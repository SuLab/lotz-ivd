#!/usr/bin/env python3
"""Module 06: Clustering with resolution optimization.

Clusters the integrated objects from Module 05 using Leiden clustering
on scANVI embeddings. Resolution is selected by multi-metric optimization
(silhouette, modularity).

Usage:
    python3 scripts/06_clustering.py                  # All objects
    python3 scripts/06_clustering.py --object NP      # Single object
    python3 scripts/06_clustering.py --validate-only  # Validation only
"""

import gc
import sys
import warnings
import traceback
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

anndata.settings.allow_write_nullable_strings = True

warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning, module='scanpy')
warnings.filterwarnings('ignore', category=DeprecationWarning)

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent.parent
INT_DIR = BASE / "data" / "integrated"
RESULTS_DIR = BASE / "results" / "integration"
CLUSTER_DIR = RESULTS_DIR / "clustering_resolution_optimization"

for d in [CLUSTER_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# ── Tier configuration ──────────────────────────────────────────────────────
# Embedding keys stored in the merged h5ad from Module 05
TIER_CONFIG = {
    "mesenchymal": {
        "embedding_key": "X_scanvi_mesenchymal",
        "cell_class_values": ["mesenchymal", "unknown"],
    },
    "non_mesenchymal": {
        "embedding_key": "X_scanvi_non_mesenchymal",
        "cell_class_values": ["non_mesenchymal"],
    },
}

OBJECTS = ["NP", "AF", "CEP", "all_cells"]


# ═════════════════════════════════════════════════════════════════════════════
# CLUSTERING RESOLUTION OPTIMIZATION
# ═════════════════════════════════════════════════════════════════════════════

def optimize_clustering_resolution(adata, embedding_key, object_name, tier_name,
                                   resolutions=None):
    """Test Leiden clustering at multiple resolutions, compute silhouette and modularity.

    Returns the selected resolution and stores results in obs.
    """
    if resolutions is None:
        # Use fewer resolutions for large datasets to reduce runtime
        if adata.shape[0] > 300000:
            resolutions = [0.4, 0.8, 1.0]  # 3 values for very large
        elif adata.shape[0] > 200000:
            resolutions = [0.2, 0.4, 0.6, 0.8, 1.0, 1.5]  # 6 values
        elif adata.shape[0] > 50000:
            resolutions = np.arange(0.2, 2.1, 0.2)  # 10 values
        else:
            resolutions = np.arange(0.1, 2.1, 0.1)  # 20 values

    print(f"  Optimizing clustering resolution ({tier_name})...")
    print(f"    Cells: {adata.shape[0]:,}, embedding: {embedding_key}")

    # Ensure neighbors are computed on this embedding
    sc.pp.neighbors(adata, use_rep=embedding_key)

    embedding = adata.obsm[embedding_key]
    results = []

    # Build igraph once for modularity computation (instead of per-resolution)
    # Skip for large datasets (>100K cells) — modularity is expensive and silhouette suffices
    ig_graph = None
    ig_weights = None
    MAX_CELLS_MODULARITY = 100000
    if adata.shape[0] <= MAX_CELLS_MODULARITY:
        try:
            import igraph as ig
            adj = adata.obsp['connectivities']
            sources, targets = adj.nonzero()
            weights = np.asarray(adj[sources, targets]).ravel()
            edgelist = np.column_stack((sources, targets))
            g = ig.Graph(n=adata.shape[0], edges=edgelist.tolist(), directed=False)
            g.es['weight'] = weights.tolist()
            ig_graph = g
            ig_weights = g.es['weight']
            print(f"    igraph built: {g.vcount()} vertices, {g.ecount()} edges")
        except Exception as e:
            print(f"    WARNING: igraph build failed: {e}, skipping modularity")
    else:
        print(f"    Skipping modularity (>{MAX_CELLS_MODULARITY:,} cells), using silhouette only")

    for res in resolutions:
        key = f'leiden_{res:.1f}'
        sc.tl.leiden(adata, resolution=res, key_added=key)
        n_clusters = adata.obs[key].nunique()

        # Silhouette score on embedding
        sil = np.nan
        if n_clusters > 1:
            try:
                from sklearn.metrics import silhouette_score
                n_sub = min(10000, adata.shape[0])
                sil = float(silhouette_score(
                    embedding, adata.obs[key].values,
                    sample_size=n_sub, random_state=42
                ))
            except Exception:
                pass

        # Modularity from the graph partition (reuse pre-built graph)
        mod_score = np.nan
        if ig_graph is not None:
            try:
                membership = adata.obs[key].astype(int).values.tolist()
                mod_score = float(ig_graph.modularity(membership, weights=ig_weights))
            except Exception:
                pass

        results.append({
            'resolution': round(res, 1),
            'n_clusters': n_clusters,
            'silhouette': sil,
            'modularity': mod_score,
        })
        print(f"    res={res:.1f}: {n_clusters} clusters, sil={sil:.3f}, mod={mod_score:.3f}")

    results_df = pd.DataFrame(results)
    results_df.to_csv(
        CLUSTER_DIR / f"{object_name}_{tier_name}_resolutions.tsv",
        sep='\t', index=False
    )

    # Select resolution: peak silhouette, or knee in modularity
    # For non-mesenchymal tier, enforce minimum resolution of 0.5 to avoid
    # too-coarse clustering that creates mixed T cell/myeloid clusters
    min_resolution = 0.5 if "non_mesenchymal" in tier_name.lower() else 0.1
    valid = results_df.dropna(subset=['silhouette'])
    if len(valid) > 0:
        valid_above_min = valid[valid['resolution'] >= min_resolution]
        if len(valid_above_min) > 0:
            best_idx = valid_above_min['silhouette'].idxmax()
            selected_res = valid_above_min.loc[best_idx, 'resolution']
        else:
            best_idx = valid['silhouette'].idxmax()
            selected_res = valid.loc[best_idx, 'resolution']
    else:
        selected_res = max(0.5, min_resolution)  # fallback

    n_sel = results_df.loc[results_df['resolution'] == selected_res, 'n_clusters'].values[0]
    print(f"    Selected resolution: {selected_res} ({n_sel} clusters)")

    # Store selected clustering
    selected_key = f'leiden_{selected_res:.1f}'
    adata.obs['leiden'] = adata.obs[selected_key].copy()
    adata.uns['selected_resolution'] = selected_res

    # Also store comparison resolutions (0.5 and 1.0)
    for ref_res in [0.5, 1.0]:
        ref_key = f'leiden_{ref_res:.1f}'
        if ref_key in adata.obs.columns:
            adata.obs[f'leiden_{ref_res}'] = adata.obs[ref_key].copy()

    # Plot
    _plot_resolution_optimization(results_df, object_name, tier_name, selected_res)

    # Clean up intermediate leiden columns
    for res in resolutions:
        key = f'leiden_{res:.1f}'
        if key in adata.obs.columns and key != selected_key:
            del adata.obs[key]

    return selected_res


def _plot_resolution_optimization(results_df, object_name, tier_name, selected_res):
    """Plot silhouette and modularity vs resolution."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    # Silhouette
    ax = axes[0]
    ax.plot(results_df['resolution'], results_df['silhouette'], 'o-', color='#2196F3')
    ax.axvline(selected_res, color='red', linestyle='--', alpha=0.7)
    ax.set_xlabel('Resolution')
    ax.set_ylabel('Silhouette Score')
    ax.set_title('Silhouette Score')

    # Modularity
    ax = axes[1]
    ax.plot(results_df['resolution'], results_df['modularity'], 'o-', color='#4CAF50')
    ax.axvline(selected_res, color='red', linestyle='--', alpha=0.7)
    ax.set_xlabel('Resolution')
    ax.set_ylabel('Modularity')
    ax.set_title('Modularity')

    # Cluster count
    ax = axes[2]
    ax.plot(results_df['resolution'], results_df['n_clusters'], 'o-', color='#FF9800')
    ax.axvline(selected_res, color='red', linestyle='--', alpha=0.7)
    ax.set_xlabel('Resolution')
    ax.set_ylabel('Number of Clusters')
    ax.set_title('Cluster Count')

    plt.suptitle(f'{object_name} — {tier_name}', fontweight='bold')
    plt.tight_layout()
    plt.savefig(CLUSTER_DIR / f"{object_name}_{tier_name}_optimization.png",
                dpi=150, bbox_inches='tight')
    plt.close()


# ═════════════════════════════════════════════════════════════════════════════
# TIER CLUSTERING
# ═════════════════════════════════════════════════════════════════════════════

def cluster_tier(adata, object_name, tier_name, tier_cfg):
    """Cluster one tier (mesenchymal or non-mesenchymal) of an object.

    Subsets cells by cell_class, runs resolution optimization on the
    tier-specific scVI embedding, and returns the subset adata with
    clustering results.
    """
    embedding_key = tier_cfg["embedding_key"]
    cell_class_vals = tier_cfg["cell_class_values"]

    # Check if the embedding exists
    if embedding_key not in adata.obsm:
        print(f"  SKIP: {object_name}/{tier_name} — embedding {embedding_key} not found")
        return None, None

    # Subset to cells of this tier
    mask = adata.obs['cell_class'].isin(cell_class_vals)
    n_cells = mask.sum()
    if n_cells == 0:
        print(f"  SKIP: {object_name}/{tier_name} — no cells with cell_class in {cell_class_vals}")
        return None, None

    print(f"\n  --- {object_name}/{tier_name}: {n_cells:,} cells ---")
    tier_adata = adata[mask].copy()

    # Drop NaN rows in the embedding (cells from the other tier have NaN)
    valid_emb = ~np.isnan(tier_adata.obsm[embedding_key]).any(axis=1)
    if valid_emb.sum() < n_cells:
        print(f"    Dropping {n_cells - valid_emb.sum()} cells with NaN embedding")
        tier_adata = tier_adata[valid_emb].copy()

    if tier_adata.shape[0] < 10:
        print(f"  SKIP: {object_name}/{tier_name} — too few cells ({tier_adata.shape[0]})")
        return None, None

    # Run resolution optimization
    selected_res = optimize_clustering_resolution(
        tier_adata, embedding_key, object_name, tier_name
    )

    return tier_adata, selected_res


# ═════════════════════════════════════════════════════════════════════════════
# PROCESS ONE OBJECT
# ═════════════════════════════════════════════════════════════════════════════

def process_object(object_name):
    """Load integrated h5ad, cluster each tier, merge clustering back, save."""
    input_path = INT_DIR / f"{object_name}.h5ad"
    if not input_path.exists():
        print(f"  ERROR: {input_path} not found, skipping {object_name}")
        return None

    print(f"\n{'='*60}")
    print(f"Clustering object: {object_name}")
    print(f"{'='*60}")

    adata = sc.read_h5ad(input_path)
    print(f"  Loaded: {adata.shape[0]:,} cells, {adata.shape[1]:,} genes")
    print(f"  obsm keys: {list(adata.obsm.keys())}")

    if 'cell_class' not in adata.obs.columns:
        print(f"  ERROR: {object_name} missing cell_class column")
        return None

    print(f"  cell_class distribution: {adata.obs['cell_class'].value_counts().to_dict()}")

    # Track tier results for merging back
    tier_results = {}

    for tier_name, tier_cfg in TIER_CONFIG.items():
        tier_adata, selected_res = cluster_tier(adata, object_name, tier_name, tier_cfg)
        if tier_adata is not None:
            tier_results[tier_name] = {
                'adata': tier_adata,
                'selected_res': selected_res,
            }

    if not tier_results:
        print(f"  WARNING: No tiers clustered for {object_name}")
        return None

    # ── Merge tier clustering results back into main adata ─────────────────
    print(f"\n  Merging tier clustering into main object...")

    # Initialize new columns as string (avoid Categorical mismatch issues)
    adata.obs['leiden'] = 'unassigned'

    for tier_name, result in tier_results.items():
        tier_adata = result['adata']
        selected_res = result['selected_res']

        # Store tier-specific leiden column
        tier_col = f'leiden_{tier_name}'
        adata.obs[tier_col] = 'NA'
        adata.obs.loc[tier_adata.obs_names, tier_col] = tier_adata.obs['leiden'].astype(str).values

        # Store comparison resolutions per tier
        for ref_res in [0.5, 1.0]:
            ref_col = f'leiden_{ref_res}'
            if ref_col in tier_adata.obs.columns:
                tier_ref_col = f'leiden_{tier_name}_{ref_res}'
                adata.obs[tier_ref_col] = 'NA'
                adata.obs.loc[tier_adata.obs_names, tier_ref_col] = \
                    tier_adata.obs[ref_col].astype(str).values

        # Build combined leiden: prefix cluster IDs with tier abbreviation
        prefix = 'M' if tier_name == 'mesenchymal' else 'NM'
        combined_labels = [f'{prefix}{c}' for c in tier_adata.obs['leiden'].astype(str).values]
        adata.obs.loc[tier_adata.obs_names, 'leiden'] = combined_labels

        # Store selected resolution in uns
        adata.uns[f'selected_resolution_{tier_name}'] = selected_res

        print(f"    {tier_name}: res={selected_res}, "
              f"{tier_adata.obs['leiden'].nunique()} clusters, "
              f"{tier_adata.shape[0]:,} cells")

    # Convert leiden to categorical
    adata.obs['leiden'] = pd.Categorical(adata.obs['leiden'])
    n_combined = adata.obs['leiden'].nunique()
    n_unassigned = (adata.obs['leiden'] == 'unassigned').sum()
    print(f"    Combined: {n_combined} clusters, {n_unassigned} unassigned cells")

    # ── Save ──────────────────────────────────────────────────────────────
    output_path = INT_DIR / f"{object_name}.h5ad"
    tmp_path = Path(str(output_path) + '.tmp')
    adata.write_h5ad(tmp_path)
    tmp_path.rename(output_path)
    print(f"  Saved: {output_path}")

    # Clean up
    for result in tier_results.values():
        del result['adata']
    del adata
    gc.collect()

    return output_path


# ═════════════════════════════════════════════════════════════════════════════
# VALIDATION
# ═════════════════════════════════════════════════════════════════════════════

def validate_clustering():
    """Validate clustering outputs. Returns (pass, messages)."""
    print("\n" + "=" * 60)
    print("Validating clustering outputs")
    print("=" * 60)

    messages = []
    all_pass = True

    for obj_name in OBJECTS:
        path = INT_DIR / f"{obj_name}.h5ad"
        if not path.exists():
            messages.append(f"FAIL: {obj_name} output missing")
            all_pass = False
            continue

        adata = sc.read_h5ad(path, backed='r')

        # Check clustering at selected resolution
        if 'leiden' in adata.obs.columns:
            n_clusters = adata.obs['leiden'].nunique()
            n_unassigned = (adata.obs['leiden'] == 'unassigned').sum()
            messages.append(f"PASS: {obj_name} has leiden column with {n_clusters} clusters")

            # Blob check
            # Count real clusters (exclude 'unassigned' and 'NA')
            real_clusters = [c for c in adata.obs['leiden'].cat.categories
                            if c not in ('unassigned', 'NA')]
            if len(real_clusters) > 1:
                messages.append(f"PASS: {obj_name} has {len(real_clusters)} real clusters (no blob)")
            else:
                messages.append(f"FAIL: {obj_name} has {len(real_clusters)} real cluster(s) (blob!)")
                all_pass = False
        else:
            messages.append(f"FAIL: {obj_name} missing leiden column")
            all_pass = False

        # Check tier-specific leiden columns
        for tier_name in TIER_CONFIG:
            tier_col = f'leiden_{tier_name}'
            if tier_col in adata.obs.columns:
                n_tier = (adata.obs[tier_col] != 'NA').sum()
                n_tier_clusters = len([c for c in adata.obs[tier_col].unique() if c != 'NA'])
                messages.append(f"PASS: {obj_name}/{tier_name} has {n_tier_clusters} clusters "
                                f"({n_tier:,} cells)")
            else:
                messages.append(f"WARNING: {obj_name} missing {tier_col} column")

        # Check comparison resolutions exist
        has_comparison = False
        for tier_name in TIER_CONFIG:
            for ref_res in [0.5, 1.0]:
                ref_col = f'leiden_{tier_name}_{ref_res}'
                if ref_col in adata.obs.columns:
                    has_comparison = True
        if has_comparison:
            messages.append(f"PASS: {obj_name} has comparison resolution columns")
        else:
            messages.append(f"WARNING: {obj_name} missing comparison resolution columns")

        # Check resolution optimization plots
        for tier_name in TIER_CONFIG:
            plot_path = CLUSTER_DIR / f"{obj_name}_{tier_name}_optimization.png"
            if plot_path.exists():
                messages.append(f"PASS: {obj_name}/{tier_name} optimization plot exists")
            else:
                messages.append(f"WARNING: {obj_name}/{tier_name} optimization plot missing")

            tsv_path = CLUSTER_DIR / f"{obj_name}_{tier_name}_resolutions.tsv"
            if tsv_path.exists():
                messages.append(f"PASS: {obj_name}/{tier_name} resolution TSV exists")
            else:
                messages.append(f"WARNING: {obj_name}/{tier_name} resolution TSV missing")

        # Check selected resolution is documented in uns
        for tier_name in TIER_CONFIG:
            res_key = f'selected_resolution_{tier_name}'
            if res_key in adata.uns:
                messages.append(f"PASS: {obj_name}/{tier_name} selected resolution = "
                                f"{adata.uns[res_key]}")
            else:
                messages.append(f"WARNING: {obj_name}/{tier_name} selected resolution not in uns")

        # ARI check: study shouldn't perfectly predict clusters
        if 'leiden' in adata.obs.columns and 'study' in adata.obs.columns:
            try:
                from sklearn.metrics import adjusted_rand_score
                # Use only real clustered cells (not unassigned)
                mask = ~adata.obs['leiden'].isin(['unassigned', 'NA'])
                obs_sub = adata.obs[mask]
                if len(obs_sub) > 0:
                    n_sub = min(10000, len(obs_sub))
                    idx = np.random.RandomState(42).choice(len(obs_sub), n_sub, replace=False)
                    sub = obs_sub.iloc[idx]
                    ari = adjusted_rand_score(sub['study'].values, sub['leiden'].values)
                    if ari < 1.0:
                        messages.append(f"PASS: {obj_name} study-cluster ARI = {ari:.3f} (< 1.0)")
                    else:
                        messages.append(f"WARNING: {obj_name} study perfectly predicts clusters "
                                        f"(ARI = {ari:.3f})")
            except Exception as e:
                messages.append(f"WARNING: {obj_name} ARI check failed: {e}")

        adata.file.close()

    for msg in messages:
        print(f"  {msg}")

    print(f"\n{'='*60}")
    print(f"OVERALL: {'ALL CHECKS PASSED' if all_pass else 'SOME CHECKS FAILED'}")
    print(f"{'='*60}")
    return all_pass, messages


# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════

def main():
    args = sys.argv[1:]

    validate_only = '--validate-only' in args

    # Parse --object
    object_filter = None
    for i, a in enumerate(args):
        if a == '--object' and i + 1 < len(args):
            object_filter = args[i + 1]

    if validate_only:
        passed, _ = validate_clustering()
        sys.exit(0 if passed else 1)

    start = datetime.now()
    print(f"Module 06: Clustering — started at {start.strftime('%Y-%m-%d %H:%M:%S')}")

    objects_to_run = [object_filter] if object_filter else OBJECTS
    results_summary = {}

    for obj_name in objects_to_run:
        try:
            output = process_object(obj_name)
            if output:
                results_summary[obj_name] = "SUCCESS"
            else:
                results_summary[obj_name] = "SKIPPED"
        except Exception as e:
            print(f"\n  ERROR processing {obj_name}: {e}")
            traceback.print_exc()
            results_summary[obj_name] = f"FAILED: {e}"

    # Run validation
    print("\n")
    passed, val_messages = validate_clustering()

    # Summary
    elapsed = datetime.now() - start
    print(f"\n{'='*60}")
    print(f"Module 06 complete — elapsed {elapsed}")
    print(f"{'='*60}")
    for obj_name, status in results_summary.items():
        print(f"  {obj_name}: {status}")

    sys.exit(0 if passed else 1)


if __name__ == '__main__':
    main()
