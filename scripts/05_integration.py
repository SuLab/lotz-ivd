#!/usr/bin/env python3
"""Module 05: Integration, clustering, and de novo annotation for IVD atlas.

Integrates cells across studies into four compartment-based objects (NP, AF,
CEP, all-cells), each with tiered scVI integration (mesenchymal and
non-mesenchymal separately). Clusters are optimized via silhouette/modularity
analysis, and cell types are annotated de novo from cluster DE markers and
canonical marker expression. CellTypist validates immune subtypes.

The all-cells object is secondary: annotations transfer from compartment-
specific objects where possible.

Usage:
    python3 scripts/05_integration.py                  # All objects
    python3 scripts/05_integration.py --object NP      # Single object
    python3 scripts/05_integration.py --validate-only  # Validation only
    python3 scripts/05_integration.py --force          # Re-run even if outputs exist
"""

import gc
import sys
import warnings
import traceback
from pathlib import Path
from datetime import datetime
from collections import Counter

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from jinja2 import Template

warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning, module='scanpy')
warnings.filterwarnings('ignore', category=UserWarning, module='scvi')
warnings.filterwarnings('ignore', category=DeprecationWarning)

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent.parent
PROC_DIR = BASE / "data" / "processed"
INT_DIR = BASE / "data" / "integrated"
MODEL_DIR = INT_DIR / "models"
RESULTS_DIR = BASE / "results" / "integration"
META_DIR = BASE / "metadata"

for d in [INT_DIR, MODEL_DIR, RESULTS_DIR,
          RESULTS_DIR / "clustering_resolution_optimization",
          RESULTS_DIR / "cluster_markers",
          RESULTS_DIR / "annotation_dotplots",
          RESULTS_DIR / "celltypist_validation"]:
    d.mkdir(parents=True, exist_ok=True)

# ── Study assignments per object ─────────────────────────────────────────────
# Each entry: (accession, compartment_filter or None for all compartments)
# compartment_filter is matched against obs['compartment']

NP_STUDIES = [
    ("GSE160756", "NP"),
    ("GSE165722", None),        # NP-only study
    ("GSE199866", "NP"),
    ("GSE205535", None),        # NP-only study
    ("GSE244889", None),        # NP-only study
    ("GSE251686", None),        # NP-only study
    ("GSE230809", "NP"),
    ("CNP0002664", None),       # NP-only study
]

AF_STUDIES = [
    ("GSE160756", "AF"),
    ("GSE199866", "AF"),        # Inner AF
    ("GSE230809", "AF"),
]

CEP_STUDIES = [
    ("GSE160756", "CEP"),
    ("GSE255768", None),        # CEP-only study
    ("GSE242443", None),        # CEP-only study (culture-expanded)
]

# All-cells: union of above plus GSE189916 (adult only, compartments not separated)
ALL_CELLS_EXTRA = [
    ("GSE189916", None),        # Adult whole IVD (neonatal excluded in Module 04)
]

# Sample-level exclusions
EXCLUDED_SAMPLES = {
    "GSE251686_NP3",  # Corrupt matrix file
}

# obs columns to keep during concatenation
OBS_COLS_KEEP = [
    "sample_id", "study", "compartment", "condition_harmonized",
    "cell_class", "pct_counts_mt",
]

# ── scVI parameters (from spec) ─────────────────────────────────────────────
SCVI_PARAMS = {
    "batch_key": "study",
    "n_latent": 20,
    "max_epochs": 200,
    "n_top_genes": 3000,
}

# ── Canonical marker panels for annotation ───────────────────────────────────
CANONICAL_MARKERS = {
    # NP subtypes
    "NP_notochordal":      ["T", "SHH", "NOG", "CD24", "KRT8", "KRT18", "KRT19"],
    "NP_mature_chondrocyte": ["ACAN", "COL2A1", "SOX9", "COMP", "PRG4"],
    "NP_stressed_degen":   ["MMP13", "ADAMTS5", "IL1B", "TNF", "VEGFA", "HIF1A"],
    "NP_fibrocartilaginous": ["COL1A1", "COL2A1", "VCAN"],
    # AF subtypes
    "AF_inner":            ["COL2A1", "ACAN", "SOX9"],
    "AF_outer":            ["COL1A1", "COL1A2", "THY1", "DCN", "LUM"],
    "AF_mechanical_stress": ["COMP", "CILP", "THBS1"],
    # Endplate
    "EP_hyaline":          ["COL2A1", "COL10A1", "SOX9"],
    "EP_ossification":     ["RUNX2", "SP7", "BGLAP"],
    # Non-mesenchymal
    "Macrophage":          ["CD68", "CD14", "CSF1R", "CD163", "CD86"],
    "T_cell":              ["CD3D", "CD3E", "CD4", "CD8A"],
    "B_cell":              ["CD79A", "MS4A1"],
    "NK_cell":             ["NKG7", "GNLY"],
    "Mast_cell":           ["KIT", "TPSAB1"],
    "Endothelial":         ["PECAM1", "VWF", "CDH5"],
    "Pericyte_SMC":        ["ACTA2", "RGS5", "PDGFRB"],
}

# Gene aliases
GENE_ALIASES = {"T": ["T", "TBXT"], "TBXT": ["TBXT", "T"], "SP7": ["SP7", "OSX"]}

# Continuous score signatures for mesenchymal continuum
CONTINUUM_SCORES = {
    "score_notochordal":   ["T", "SHH", "NOG", "CD24", "KRT8", "KRT18", "KRT19"],
    "score_degenerative":  ["MMP13", "ADAMTS5", "IL1B", "TNF", "VEGFA", "HIF1A"],
    "score_fibrotic":      ["COL1A1", "COL3A1", "FN1", "VCAN", "ACTA2", "TGFB1"],
}


# ═════════════════════════════════════════════════════════════════════════════
# GENE RESOLUTION
# ═════════════════════════════════════════════════════════════════════════════

def _resolve_gene(name, var_names):
    """Resolve gene name to form present in var_names."""
    if name in var_names:
        return name
    for alias in GENE_ALIASES.get(name, []):
        if alias in var_names:
            return alias
    return None


def _resolve_gene_list(gene_list, var_names):
    """Resolve genes, returning (resolved, missing)."""
    var_set = set(var_names)
    resolved, missing = [], []
    for g in gene_list:
        match = _resolve_gene(g, var_set)
        (resolved if match else missing).append(match or g)
    return resolved, missing


# ═════════════════════════════════════════════════════════════════════════════
# STUDY ASSIGNMENTS
# ═════════════════════════════════════════════════════════════════════════════

def get_study_assignments(object_name):
    """Return list of (accession, compartment_filter) for an object."""
    if object_name == "NP":
        return NP_STUDIES
    elif object_name == "AF":
        return AF_STUDIES
    elif object_name == "CEP":
        return CEP_STUDIES
    elif object_name == "all_cells":
        # Union of all compartment-specific studies + extra
        seen = set()
        assignments = []
        for study_list in [NP_STUDIES, AF_STUDIES, CEP_STUDIES, ALL_CELLS_EXTRA]:
            for acc, comp in study_list:
                key = (acc, comp)
                if key not in seen:
                    seen.add(key)
                    assignments.append((acc, comp))
        return assignments
    else:
        raise ValueError(f"Unknown object: {object_name}")


# ═════════════════════════════════════════════════════════════════════════════
# MEMORY-EFFICIENT LOADING & CONCATENATION
# ═════════════════════════════════════════════════════════════════════════════

def load_and_build_object(object_name):
    """Load and concatenate cells for one integrated object.

    Filters by compartment where specified, excludes excluded samples,
    and splits by cell_class.

    Returns AnnData with .X (log-norm), .layers['counts'], selected obs cols.
    """
    assignments = get_study_assignments(object_name)
    obs_cols = OBS_COLS_KEEP

    print(f"  Loading cells for {object_name} object...")
    print(f"    Studies: {[a for a, _ in assignments]}")

    # Phase 1: scan to identify matching cells and common genes
    gene_sets = []
    cell_counts = {}

    for acc, comp_filter in assignments:
        path = PROC_DIR / f"{acc}.h5ad"
        if not path.exists():
            print(f"    WARNING: {path} not found, skipping")
            continue
        adata = sc.read_h5ad(path, backed='r')
        obs_df = adata.obs

        # Build filter mask
        mask = pd.Series(True, index=obs_df.index)

        # Compartment filter
        if comp_filter is not None and 'compartment' in obs_df.columns:
            mask &= obs_df['compartment'].str.upper().str.contains(comp_filter.upper(), na=False)

        # Sample exclusions
        if 'sample_id' in obs_df.columns:
            mask &= ~obs_df['sample_id'].isin(EXCLUDED_SAMPLES)

        # Must have cell_class from Module 04
        if 'cell_class' in obs_df.columns:
            mask &= obs_df['cell_class'].isin(['mesenchymal', 'non_mesenchymal'])
        else:
            print(f"    WARNING: {acc} missing cell_class column, skipping")
            adata.file.close()
            continue

        n_match = mask.sum()
        if n_match > 0:
            cell_counts[(acc, comp_filter)] = n_match
            gene_sets.append(set(adata.var_names))
        adata.file.close()

    if not cell_counts:
        print("  WARNING: No matching cells found!")
        return None

    common_genes = sorted(set.intersection(*gene_sets))
    total_cells = sum(cell_counts.values())
    print(f"    Found {total_cells:,} cells across {len(cell_counts)} study-compartment combos")
    print(f"    Common genes: {len(common_genes):,}")

    # Phase 2: load, subset, collect
    adatas = []
    for acc, comp_filter in assignments:
        if (acc, comp_filter) not in cell_counts:
            continue
        path = PROC_DIR / f"{acc}.h5ad"
        adata = sc.read_h5ad(path)

        mask = pd.Series(True, index=adata.obs.index)
        if comp_filter is not None and 'compartment' in adata.obs.columns:
            mask &= adata.obs['compartment'].str.upper().str.contains(comp_filter.upper(), na=False)
        if 'sample_id' in adata.obs.columns:
            mask &= ~adata.obs['sample_id'].isin(EXCLUDED_SAMPLES)
        if 'cell_class' in adata.obs.columns:
            mask &= adata.obs['cell_class'].isin(['mesenchymal', 'non_mesenchymal'])

        adata_sub = adata[mask, common_genes].copy()

        # Keep only selected obs columns
        cols_available = [c for c in obs_cols if c in adata_sub.obs.columns]
        adata_sub.obs = adata_sub.obs[cols_available].copy()

        # Ensure study column
        if 'study' not in adata_sub.obs.columns:
            adata_sub.obs['study'] = acc

        # Ensure counts layer, int32
        if 'counts' in adata_sub.layers:
            if adata_sub.layers['counts'].dtype != np.int32:
                adata_sub.layers['counts'] = adata_sub.layers['counts'].astype(np.int32)
        else:
            adata_sub.layers['counts'] = adata_sub.X.copy().astype(np.int32)

        adatas.append(adata_sub)
        print(f"    {acc} ({comp_filter or 'all'}): {adata_sub.shape[0]:,} cells")
        del adata

    # Concatenate
    print("  Concatenating...")
    adata = ad.concat(adatas, merge='same')
    del adatas

    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(adata.X)
    if not sparse.issparse(adata.layers['counts']):
        adata.layers['counts'] = sparse.csr_matrix(adata.layers['counts'])

    print(f"    Final shape: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")
    return adata


# ═════════════════════════════════════════════════════════════════════════════
# PREPARE OBJECT (HVG, PCA, NEIGHBORS, UMAP BASELINE)
# ═════════════════════════════════════════════════════════════════════════════

def prepare_object(adata, n_top_genes=3000):
    """Re-derive HVGs, PCA, neighbors, UMAP on concatenated object.

    Stores unintegrated baseline in obsm['X_pca_unintegrated'] and
    obsm['X_umap_unintegrated'].
    """
    print("  Preparing object (HVG, PCA, neighbors, UMAP)...")

    n_top = min(n_top_genes, adata.shape[1] - 1)

    # Batch-aware HVG selection
    hvg_done = False
    MIN_CELLS_HVG_BATCH = 50
    if 'study' in adata.obs.columns:
        study_counts = adata.obs['study'].value_counts()
        large_studies = study_counts[study_counts >= MIN_CELLS_HVG_BATCH].index.tolist()
        if len(large_studies) > 1:
            try:
                adata_hvg = adata[adata.obs['study'].isin(large_studies)].copy()
                sc.pp.highly_variable_genes(
                    adata_hvg, n_top_genes=n_top, flavor='seurat_v3',
                    batch_key='study', layer='counts'
                )
                adata.var['highly_variable'] = False
                adata.var.loc[adata_hvg.var_names[adata_hvg.var['highly_variable']],
                              'highly_variable'] = True
                for col in ['highly_variable_rank', 'highly_variable_nbatches',
                            'means', 'variances', 'variances_norm']:
                    if col in adata_hvg.var.columns:
                        adata.var[col] = adata_hvg.var[col]
                del adata_hvg
                hvg_done = True
            except Exception as e:
                print(f"    WARNING: batch-aware HVG failed ({e}), falling back")

    if not hvg_done:
        try:
            sc.pp.highly_variable_genes(
                adata, n_top_genes=n_top, flavor='seurat_v3', layer='counts'
            )
        except ValueError:
            print("    WARNING: seurat_v3 HVG failed, falling back to cell_ranger flavor")
            sc.pp.highly_variable_genes(
                adata, n_top_genes=n_top, flavor='cell_ranger'
            )
    n_hvg = adata.var['highly_variable'].sum()
    print(f"    HVGs: {n_hvg}")

    # PCA on HVG subset
    adata_pca = adata[:, adata.var['highly_variable']].copy()
    gc.collect()
    sc.pp.scale(adata_pca, max_value=10)
    n_comps = min(50, min(adata_pca.shape) - 1)
    sc.tl.pca(adata_pca, n_comps=n_comps, svd_solver='arpack')

    adata.obsm['X_pca'] = adata_pca.obsm['X_pca']
    adata.uns['pca'] = adata_pca.uns['pca']
    adata.varm['PCs'] = np.zeros((adata.shape[1], n_comps))
    hvg_idx = np.where(adata.var['highly_variable'])[0]
    adata.varm['PCs'][hvg_idx] = adata_pca.varm['PCs']
    del adata_pca

    # Unintegrated baseline
    n_pcs = min(30, n_comps)
    sc.pp.neighbors(adata, n_pcs=n_pcs)
    sc.tl.umap(adata)
    adata.obsm['X_pca_unintegrated'] = adata.obsm['X_pca'].copy()
    adata.obsm['X_umap_unintegrated'] = adata.obsm['X_umap'].copy()
    print(f"    PCA: {n_comps} components, baseline UMAP computed")

    return adata


# ═════════════════════════════════════════════════════════════════════════════
# scVI INTEGRATION
# ═════════════════════════════════════════════════════════════════════════════

def run_scvi(adata, batch_key='study', n_latent=20, max_epochs=200,
             model_dir=None, label="scvi"):
    """Run scVI integration. Stores obsm[f'X_scvi_{label}'] and UMAP."""
    import scvi as scvi_module

    print(f"  Running scVI: batch_key={batch_key}, n_latent={n_latent}")

    scvi_module.model.SCVI.setup_anndata(
        adata, layer='counts', batch_key=batch_key,
    )

    model = scvi_module.model.SCVI(
        adata, n_latent=n_latent,
        dispersion='gene-batch', gene_likelihood='nb',
    )

    model.train(
        max_epochs=max_epochs,
        early_stopping=True,
        early_stopping_patience=10,
        early_stopping_monitor='elbo_validation',
        train_size=0.9,
        batch_size=256,
    )

    embedding_key = f'X_scvi_{label}'
    adata.obsm[embedding_key] = model.get_latent_representation()
    print(f"    scVI complete: {model.history['elbo_train'].shape[0]} epochs")

    # Neighbors and UMAP on this embedding
    sc.pp.neighbors(adata, use_rep=embedding_key)
    sc.tl.umap(adata)
    adata.obsm[f'X_umap_scvi_{label}'] = adata.obsm['X_umap'].copy()

    if model_dir is not None:
        save_path = model_dir / f"scvi_{label}"
        model.save(str(save_path), overwrite=True)

    return model, embedding_key


# ═════════════════════════════════════════════════════════════════════════════
# INTEGRATION METRICS
# ═════════════════════════════════════════════════════════════════════════════

def _quick_ilisi(embedding, batch_labels, n_sample=5000, k=30):
    """Quick iLISI estimate on a random subsample."""
    from sklearn.neighbors import NearestNeighbors

    n = embedding.shape[0]
    if n > n_sample:
        rng = np.random.RandomState(42)
        idx = rng.choice(n, n_sample, replace=False)
        emb_sub = embedding[idx]
        labels_sub = np.asarray(batch_labels)[idx]
    else:
        emb_sub = embedding
        labels_sub = np.asarray(batch_labels)

    nn = NearestNeighbors(n_neighbors=k, algorithm='auto')
    nn.fit(emb_sub)
    indices = nn.kneighbors(emb_sub, return_distance=False)

    ilisi_scores = []
    for i in range(len(emb_sub)):
        neighbor_labels = labels_sub[indices[i]]
        counts = Counter(neighbor_labels)
        total = sum(counts.values())
        p_sq = sum((c / total) ** 2 for c in counts.values())
        ilisi_scores.append(1.0 / p_sq)

    return np.mean(ilisi_scores)


def compute_integration_metrics(adata, embedding_key, batch_key='study',
                                max_cells=30_000):
    """Compute integration quality metrics (iLISI, batch-ASW, condition-ASW).

    Returns dict with metric values.
    """
    print(f"    Computing metrics for {embedding_key}...")
    metrics = {}

    if embedding_key not in adata.obsm:
        print(f"    WARNING: {embedding_key} not in obsm")
        return metrics

    embedding = adata.obsm[embedding_key]
    batch = adata.obs[batch_key].values

    # Subsample if needed
    n_cells = embedding.shape[0]
    if n_cells > max_cells:
        rng = np.random.RandomState(42)
        idx = rng.choice(n_cells, max_cells, replace=False)
        embedding = embedding[idx]
        batch = batch[idx]
        cond = adata.obs['condition_harmonized'].values[idx] if 'condition_harmonized' in adata.obs.columns else None
        print(f"    Subsampled {n_cells:,} -> {max_cells:,} for metrics")
    else:
        cond = adata.obs['condition_harmonized'].values if 'condition_harmonized' in adata.obs.columns else None

    metrics['n_cells_metrics'] = int(embedding.shape[0])

    # iLISI
    try:
        metrics['iLISI'] = float(_quick_ilisi(embedding, batch))
    except Exception as e:
        print(f"    WARNING: iLISI failed: {e}")
        metrics['iLISI'] = np.nan

    # Batch ASW
    try:
        from sklearn.metrics import silhouette_score
        if len(set(batch)) > 1:
            metrics['batch_ASW'] = float(silhouette_score(
                embedding, batch, sample_size=min(5000, len(embedding)),
                random_state=42
            ))
        else:
            metrics['batch_ASW'] = np.nan
    except Exception as e:
        print(f"    WARNING: batch_ASW failed: {e}")
        metrics['batch_ASW'] = np.nan

    # Condition ASW
    try:
        if cond is not None:
            valid = pd.notna(cond)
            if valid.sum() > 100 and len(set(cond[valid])) > 1:
                metrics['condition_ASW'] = float(silhouette_score(
                    embedding[valid], cond[valid],
                    sample_size=min(5000, valid.sum()),
                    random_state=42
                ))
            else:
                metrics['condition_ASW'] = np.nan
        else:
            metrics['condition_ASW'] = np.nan
    except Exception as e:
        print(f"    WARNING: condition_ASW failed: {e}")
        metrics['condition_ASW'] = np.nan

    return metrics


# ═════════════════════════════════════════════════════════════════════════════
# CLUSTERING WITH RESOLUTION OPTIMIZATION
# ═════════════════════════════════════════════════════════════════════════════

def optimize_clustering_resolution(adata, embedding_key, object_name, tier_name,
                                   resolutions=None):
    """Test Leiden clustering at multiple resolutions, compute silhouette and modularity.

    Returns the selected resolution and stores results in obs.
    """
    if resolutions is None:
        resolutions = np.arange(0.1, 2.1, 0.1)

    print(f"  Optimizing clustering resolution ({tier_name})...")

    # Ensure neighbors are computed on this embedding
    sc.pp.neighbors(adata, use_rep=embedding_key)

    embedding = adata.obsm[embedding_key]
    results = []

    # Build igraph once for modularity computation (instead of per-resolution)
    ig_graph = None
    ig_weights = None
    try:
        import leidenalg
        import igraph as ig
        adj = adata.obsp['connectivities']
        sources, targets = adj.nonzero()
        weights = np.asarray(adj[sources, targets]).ravel()
        g = ig.Graph(n=adata.shape[0], edges=list(zip(sources.tolist(), targets.tolist())),
                     directed=False)
        g.es['weight'] = weights.tolist()
        ig_graph = g
        ig_weights = g.es['weight']
    except Exception:
        pass

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

    results_df = pd.DataFrame(results)
    results_df.to_csv(
        RESULTS_DIR / "clustering_resolution_optimization" / f"{object_name}_{tier_name}_resolutions.tsv",
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

    print(f"    Selected resolution: {selected_res} "
          f"({results_df.loc[results_df['resolution'] == selected_res, 'n_clusters'].values[0]} clusters)")

    # Store selected clustering
    selected_key = f'leiden_{selected_res:.1f}'
    adata.obs['leiden'] = adata.obs[selected_key].copy()
    adata.uns['selected_resolution'] = selected_res

    # Also store a couple of reference resolutions
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
    plt.savefig(RESULTS_DIR / "clustering_resolution_optimization" /
                f"{object_name}_{tier_name}_optimization.png",
                dpi=150, bbox_inches='tight')
    plt.close()


# ═════════════════════════════════════════════════════════════════════════════
# DE NOVO ANNOTATION
# ═════════════════════════════════════════════════════════════════════════════

def compute_cluster_markers(adata, object_name, tier_name, n_genes=50):
    """Compute DE markers per cluster using Wilcoxon rank-sum test."""
    print(f"  Computing cluster markers ({tier_name})...")

    sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon',
                            n_genes=n_genes)

    # Save to TSV
    marker_df = sc.get.rank_genes_groups_df(adata, group=None)
    marker_df.to_csv(
        RESULTS_DIR / "cluster_markers" / f"{object_name}_{tier_name}_markers.tsv",
        sep='\t', index=False
    )
    print(f"    Saved markers: {len(marker_df)} rows")
    return marker_df


def generate_annotation_dotplots(adata, object_name, tier_name, marker_panel):
    """Generate dot plots of canonical markers per cluster."""
    var_names = set(adata.var_names)
    resolved_markers = {}
    for group, genes in marker_panel.items():
        resolved = [_resolve_gene(g, var_names) for g in genes
                     if _resolve_gene(g, var_names) is not None]
        if resolved:
            resolved_markers[group] = resolved

    if not resolved_markers and adata.obs['leiden'].nunique() <= 1:
        return

    try:
        # Flat list for dotplot, preserving group structure
        all_genes = []
        for genes in resolved_markers.values():
            for g in genes:
                if g not in all_genes:
                    all_genes.append(g)

        if all_genes:
            sc.pl.dotplot(adata, var_names=all_genes[:50], groupby='leiden', show=False)
            plt.savefig(
                RESULTS_DIR / "annotation_dotplots" /
                f"{object_name}_{tier_name}_canonical_markers.pdf",
                bbox_inches='tight'
            )
            plt.close()
    except Exception as e:
        print(f"    WARNING: Dotplot failed: {e}")


def annotate_clusters_de_novo(adata, object_name, tier_name):
    """Annotate clusters using DE markers + canonical marker expression.

    Assigns obs['cell_type'], obs['cell_type_confidence'], obs['annotation_evidence'].
    """
    print(f"  Annotating clusters de novo ({tier_name})...")

    # Select appropriate marker panel based on tier and object
    if tier_name == "non_mesenchymal":
        panel = {k: v for k, v in CANONICAL_MARKERS.items()
                 if k in ("Macrophage", "T_cell", "B_cell", "NK_cell",
                           "Mast_cell", "Endothelial", "Pericyte_SMC")}
    else:
        if object_name in ("NP", "all_cells"):
            panel = {k: v for k, v in CANONICAL_MARKERS.items()
                     if k.startswith("NP_")}
        elif object_name == "AF":
            panel = {k: v for k, v in CANONICAL_MARKERS.items()
                     if k.startswith("AF_")}
        elif object_name == "CEP":
            panel = {k: v for k, v in CANONICAL_MARKERS.items()
                     if k.startswith("EP_")}
        else:
            panel = CANONICAL_MARKERS

    # Generate dotplots
    generate_annotation_dotplots(adata, object_name, tier_name, panel)

    # Score each cluster for each canonical marker set
    var_names = set(adata.var_names)
    cluster_scores = {}

    for cluster in adata.obs['leiden'].unique():
        mask = adata.obs['leiden'] == cluster
        cluster_scores[cluster] = {}

        for ct_name, genes in panel.items():
            resolved = [_resolve_gene(g, var_names) for g in genes
                         if _resolve_gene(g, var_names) is not None]
            if not resolved:
                cluster_scores[cluster][ct_name] = 0.0
                continue

            # Expression scoring: weight marker specificity over diffuse expression.
            # For each gene, compute in-cluster vs out-of-cluster expression ratio
            # to penalize ubiquitously expressed genes (e.g., CD68 in stressed cells).
            expr_in = adata[mask, resolved].X
            expr_out = adata[~mask, resolved].X
            if sparse.issparse(expr_in):
                expr_in = expr_in.toarray()
            if sparse.issparse(expr_out):
                expr_out = expr_out.toarray()
            # Per-gene scores: fraction expressing * specificity ratio
            gene_scores = []
            for g_idx in range(expr_in.shape[1]):
                frac_in = (expr_in[:, g_idx] > 0).mean()
                mean_in = expr_in[:, g_idx].mean()
                mean_out = expr_out[:, g_idx].mean()
                # Specificity: how enriched is this gene in-cluster vs out-of-cluster
                specificity = mean_in / (mean_out + 0.01)  # avoid div by zero
                gene_scores.append(frac_in * min(specificity, 5.0))  # cap at 5x
            cluster_scores[cluster][ct_name] = float(np.mean(gene_scores)) if gene_scores else 0.0

    # Assign labels
    cell_types = pd.Series("unassigned", index=adata.obs_names, dtype=object)
    confidences = pd.Series("low", index=adata.obs_names, dtype=object)
    evidences = pd.Series("", index=adata.obs_names, dtype=object)

    for cluster in adata.obs['leiden'].unique():
        mask = adata.obs['leiden'] == cluster
        scores = cluster_scores[cluster]

        if not scores or all(v == 0 for v in scores.values()):
            continue

        sorted_scores = sorted(scores.items(), key=lambda x: x[1], reverse=True)
        best_name, best_score = sorted_scores[0]
        second_score = sorted_scores[1][1] if len(sorted_scores) > 1 else 0

        margin = best_score - second_score

        if best_score > 0.05:
            cell_types.loc[mask] = best_name
            if margin > 0.1:
                confidences.loc[mask] = "high"
            elif margin > 0.03:
                confidences.loc[mask] = "medium"
            else:
                confidences.loc[mask] = "low"
            evidences.loc[mask] = f"top_marker={best_name}(score={best_score:.3f})"

    adata.obs['cell_type'] = cell_types.values
    adata.obs['cell_type_confidence'] = confidences.values
    adata.obs['annotation_evidence'] = evidences.values

    # Summary
    ct_counts = adata.obs['cell_type'].value_counts()
    print(f"    Annotation summary ({tier_name}):")
    for ct, n in ct_counts.items():
        pct = n / adata.shape[0] * 100
        print(f"      {ct}: {n} ({pct:.1f}%)")


def compute_continuous_scores(adata):
    """Compute continuous gene signature scores for mesenchymal continuum."""
    var_names = set(adata.var_names)
    for score_name, genes in CONTINUUM_SCORES.items():
        resolved = [_resolve_gene(g, var_names) for g in genes
                     if _resolve_gene(g, var_names) is not None]
        if resolved:
            try:
                sc.tl.score_genes(adata, gene_list=resolved, score_name=score_name)
            except Exception:
                adata.obs[score_name] = 0.0
        else:
            adata.obs[score_name] = 0.0


# ═════════════════════════════════════════════════════════════════════════════
# CELLTYPIST VALIDATION (immune subtypes)
# ═════════════════════════════════════════════════════════════════════════════

def validate_immune_with_celltypist(adata, object_name):
    """Run CellTypist on non-mesenchymal cells as a validation check.

    Compares CellTypist predictions against de novo labels per cluster.
    Saves concordance table. Flags disagreements in annotation_evidence.
    """
    print(f"  Running CellTypist validation on non-mesenchymal cells...")

    try:
        import celltypist
        from celltypist import models as ct_models
    except ImportError:
        print("    WARNING: celltypist not installed, skipping validation")
        return None

    try:
        model = ct_models.Model.load(model="Immune_All_Low.pkl")
    except Exception as e:
        print(f"    WARNING: Could not load CellTypist model: {e}")
        return None

    try:
        predictions = celltypist.annotate(adata, model=model, majority_voting=True)
        result_adata = predictions.to_adata()
        adata.obs['celltypist_prediction'] = result_adata.obs['majority_voting'].values
    except Exception as e:
        print(f"    WARNING: CellTypist annotation failed: {e}")
        return None

    # Build concordance table: cluster × CellTypist majority prediction
    concordance = []
    for cluster in adata.obs['leiden'].unique():
        mask = adata.obs['leiden'] == cluster
        de_novo_label = adata.obs.loc[mask, 'cell_type'].mode()
        de_novo_label = de_novo_label.iloc[0] if len(de_novo_label) > 0 else "unknown"
        ct_labels = adata.obs.loc[mask, 'celltypist_prediction']
        ct_majority = ct_labels.mode()
        ct_majority = ct_majority.iloc[0] if len(ct_majority) > 0 else "unknown"
        ct_agreement = (ct_labels == ct_majority).mean()

        agrees = _labels_agree(de_novo_label, ct_majority)

        concordance.append({
            'cluster': cluster,
            'n_cells': mask.sum(),
            'de_novo_label': de_novo_label,
            'celltypist_majority': ct_majority,
            'celltypist_agreement_pct': f"{ct_agreement * 100:.1f}",
            'concordant': agrees,
        })

        # Flag disagreements
        if not agrees:
            mask_idx = mask[mask].index
            current_evidence = adata.obs.loc[mask_idx, 'annotation_evidence']
            adata.obs.loc[mask_idx, 'annotation_evidence'] = (
                current_evidence + f"; CELLTYPIST_DISAGREES={ct_majority}"
            )
            print(f"    WARNING: Cluster {cluster}: de novo={de_novo_label}, "
                  f"CellTypist={ct_majority} — flagged for review")

    conc_df = pd.DataFrame(concordance)
    conc_df.to_csv(
        RESULTS_DIR / "celltypist_validation" / f"{object_name}_concordance.tsv",
        sep='\t', index=False
    )
    n_agree = conc_df['concordant'].sum()
    print(f"    Concordance: {n_agree}/{len(conc_df)} clusters agree")
    return conc_df


def _labels_agree(de_novo, celltypist):
    """Check if de novo and CellTypist labels broadly agree."""
    dn = de_novo.lower()
    ct = celltypist.lower()

    # Map common label patterns
    mappings = {
        'macrophage': ['macro', 'mono', 'dc', 'dendritic'],
        't_cell': ['t cell', 'cd4', 'cd8', 'tem', 'tcm', 'treg'],
        'b_cell': ['b cell', 'plasma', 'pre-b', 'pro-b'],
        'nk_cell': ['nk'],
        'mast_cell': ['mast'],
        'endothelial': ['endothelial'],
        'pericyte_smc': ['pericyte', 'smooth muscle', 'smc'],
    }

    dn_category = None
    ct_category = None

    for key, patterns in mappings.items():
        if any(p in dn for p in patterns) or key in dn:
            dn_category = key
        if any(p in ct for p in patterns):
            ct_category = key

    # If both match a known category, they must match the same one
    if dn_category is not None and ct_category is not None:
        return dn_category == ct_category

    # If only one matches a known category, they disagree
    if dn_category is not None or ct_category is not None:
        return False

    # If neither matches known patterns, can't determine disagreement
    return True


# ═════════════════════════════════════════════════════════════════════════════
# TIER MERGING
# ═════════════════════════════════════════════════════════════════════════════

def merge_tiers(mes_adata, non_mes_adata, mes_embedding_key, non_mes_embedding_key):
    """Merge mesenchymal and non-mesenchymal tiers into a single object.

    Stores tier-specific embeddings as X_scvi_mesenchymal and X_scvi_non_mesenchymal.
    """
    print("  Merging tiers...")

    # Store tier-specific embeddings with NaN for the other tier
    n_latent_mes = mes_adata.obsm[mes_embedding_key].shape[1]
    n_latent_non = non_mes_adata.obsm[non_mes_embedding_key].shape[1]

    # Pad non-mesenchymal with NaN for mesenchymal embedding and vice versa
    mes_adata.obsm['X_scvi_mesenchymal'] = mes_adata.obsm[mes_embedding_key].copy()
    mes_adata.obsm['X_scvi_non_mesenchymal'] = np.full(
        (mes_adata.shape[0], n_latent_non), np.nan
    )

    non_mes_adata.obsm['X_scvi_non_mesenchymal'] = non_mes_adata.obsm[non_mes_embedding_key].copy()
    non_mes_adata.obsm['X_scvi_mesenchymal'] = np.full(
        (non_mes_adata.shape[0], n_latent_mes), np.nan
    )

    # Ensure both have the same obs columns
    common_obs = list(set(mes_adata.obs.columns) & set(non_mes_adata.obs.columns))
    mes_adata.obs = mes_adata.obs[common_obs].copy()
    non_mes_adata.obs = non_mes_adata.obs[common_obs].copy()

    # Ensure same var
    common_genes = sorted(set(mes_adata.var_names) & set(non_mes_adata.var_names))
    mes_sub = mes_adata[:, common_genes].copy()
    non_mes_sub = non_mes_adata[:, common_genes].copy()

    # Concatenate
    merged = ad.concat([mes_sub, non_mes_sub], merge='same')

    if not sparse.issparse(merged.X):
        merged.X = sparse.csr_matrix(merged.X)

    print(f"    Merged: {merged.shape[0]:,} cells "
          f"({mes_adata.shape[0]:,} mes + {non_mes_adata.shape[0]:,} non-mes)")

    del mes_sub, non_mes_sub
    return merged


# ═════════════════════════════════════════════════════════════════════════════
# CHECKPOINT SAVE
# ═════════════════════════════════════════════════════════════════════════════

def _save_checkpoint(adata, output_path, label):
    """Atomic save of adata checkpoint."""
    import anndata
    anndata.settings.allow_write_nullable_strings = True
    tmp_path = Path(str(output_path) + '.tmp')
    adata.write_h5ad(tmp_path)
    tmp_path.rename(output_path)
    print(f"  Checkpoint saved ({label}): {output_path}")
    sys.stdout.flush()


# ═════════════════════════════════════════════════════════════════════════════
# PROCESS ONE OBJECT (MAIN ORCHESTRATOR)
# ═════════════════════════════════════════════════════════════════════════════

def process_object(object_name, force=False):
    """Full pipeline for one integrated object.

    Steps: load → split by cell_class → integrate tiers → cluster →
    annotate → CellTypist validation → compute continuum scores → merge → save.
    """
    output_path = INT_DIR / f"{object_name}.h5ad"

    if output_path.exists() and not force:
        print(f"\n=== {object_name}: output exists, skipping (use --force) ===")
        return output_path, {}

    print(f"\n{'='*60}")
    print(f"Processing object: {object_name}")
    print(f"{'='*60}")

    # Load all cells for this object
    adata_all = load_and_build_object(object_name)
    if adata_all is None or adata_all.shape[0] < 50:
        print(f"  ERROR: Too few cells for {object_name}")
        return None, {}

    # Split by cell_class
    mes_mask = adata_all.obs['cell_class'] == 'mesenchymal'
    non_mes_mask = adata_all.obs['cell_class'] == 'non_mesenchymal'

    n_mes = mes_mask.sum()
    n_non = non_mes_mask.sum()
    print(f"  Mesenchymal: {n_mes:,} cells")
    print(f"  Non-mesenchymal: {n_non:,} cells")

    all_metrics = {}
    mes_adata = None
    non_mes_adata = None
    mes_emb_key = None
    non_mes_emb_key = None

    # ── Tier A: Mesenchymal ──────────────────────────────────────────────
    if n_mes >= 50:
        print(f"\n  --- Tier A: Mesenchymal ({object_name}) ---")
        mes_adata = adata_all[mes_mask].copy()
        mes_adata = prepare_object(mes_adata, n_top_genes=SCVI_PARAMS['n_top_genes'])

        _, mes_emb_key = run_scvi(
            mes_adata,
            batch_key=SCVI_PARAMS['batch_key'],
            n_latent=SCVI_PARAMS['n_latent'],
            max_epochs=SCVI_PARAMS['max_epochs'],
            model_dir=MODEL_DIR,
            label=f"{object_name}_mesenchymal"
        )

        # Integration metrics
        mes_metrics = compute_integration_metrics(mes_adata, mes_emb_key)
        all_metrics['mesenchymal'] = mes_metrics

        # Clustering
        optimize_clustering_resolution(mes_adata, mes_emb_key, object_name, "mesenchymal")

        # DE markers
        compute_cluster_markers(mes_adata, object_name, "mesenchymal")

        # Annotation
        annotate_clusters_de_novo(mes_adata, object_name, "mesenchymal")

        # Continuum scores
        compute_continuous_scores(mes_adata)

        gc.collect()
    else:
        print(f"  Skipping mesenchymal tier: only {n_mes} cells")

    # ── Tier B: Non-mesenchymal ──────────────────────────────────────────
    # Need enough cells for HVG + scVI to work (at least 200 cells)
    MIN_CELLS_NON_MES = 200
    if n_non >= MIN_CELLS_NON_MES:
        print(f"\n  --- Tier B: Non-mesenchymal ({object_name}) ---")
        non_mes_adata = adata_all[non_mes_mask].copy()
        non_mes_adata = prepare_object(non_mes_adata, n_top_genes=SCVI_PARAMS['n_top_genes'])

        _, non_mes_emb_key = run_scvi(
            non_mes_adata,
            batch_key=SCVI_PARAMS['batch_key'],
            n_latent=SCVI_PARAMS['n_latent'],
            max_epochs=SCVI_PARAMS['max_epochs'],
            model_dir=MODEL_DIR,
            label=f"{object_name}_non_mesenchymal"
        )

        # Integration metrics
        non_mes_metrics = compute_integration_metrics(non_mes_adata, non_mes_emb_key)
        all_metrics['non_mesenchymal'] = non_mes_metrics

        # Clustering
        optimize_clustering_resolution(non_mes_adata, non_mes_emb_key,
                                       object_name, "non_mesenchymal")

        # DE markers
        compute_cluster_markers(non_mes_adata, object_name, "non_mesenchymal")

        # Annotation
        annotate_clusters_de_novo(non_mes_adata, object_name, "non_mesenchymal")

        # CellTypist validation
        validate_immune_with_celltypist(non_mes_adata, object_name)

        gc.collect()
    else:
        print(f"  Skipping non-mesenchymal tier: only {n_non} cells")

    del adata_all
    gc.collect()

    # ── Merge tiers ──────────────────────────────────────────────────────
    if mes_adata is not None and non_mes_adata is not None:
        merged = merge_tiers(mes_adata, non_mes_adata, mes_emb_key, non_mes_emb_key)
        del mes_adata, non_mes_adata
    elif mes_adata is not None:
        merged = mes_adata
        merged.obsm['X_scvi_mesenchymal'] = merged.obsm[mes_emb_key].copy()
        del mes_adata
    elif non_mes_adata is not None:
        merged = non_mes_adata
        merged.obsm['X_scvi_non_mesenchymal'] = merged.obsm[non_mes_emb_key].copy()
        del non_mes_adata
    else:
        print(f"  ERROR: No tiers processed for {object_name}")
        return None, {}

    # ── Generate UMAPs ───────────────────────────────────────────────────
    _plot_object_umaps(merged, object_name)

    # ── Save ─────────────────────────────────────────────────────────────
    _save_checkpoint(merged, output_path, "final")
    print(f"  {object_name} complete: {merged.shape[0]:,} cells, "
          f"{merged.obs['cell_type'].nunique()} cell types")

    del merged
    gc.collect()
    return output_path, all_metrics


def process_all_cells_secondary(force=False):
    """Process all-cells as a secondary object.

    Cells from compartment-specific objects carry their annotations.
    Only GSE189916 adult cells get de novo annotation.
    """
    output_path = INT_DIR / "all_cells.h5ad"

    if output_path.exists() and not force:
        print(f"\n=== all_cells: output exists, skipping (use --force) ===")
        return output_path, {}

    print(f"\n{'='*60}")
    print(f"Processing object: all_cells (secondary)")
    print(f"{'='*60}")

    # Load all cells
    adata_all = load_and_build_object("all_cells")
    if adata_all is None:
        return None, {}

    # Transfer annotations from compartment-specific objects
    print("  Transferring annotations from compartment-specific objects...")
    adata_all.obs['cell_type'] = "unassigned"
    adata_all.obs['cell_type_confidence'] = "low"
    adata_all.obs['annotation_evidence'] = ""

    n_transferred = 0
    for obj_name in ["NP", "AF", "CEP"]:
        obj_path = INT_DIR / f"{obj_name}.h5ad"
        if not obj_path.exists():
            continue
        obj_adata = sc.read_h5ad(obj_path)
        if 'cell_type' not in obj_adata.obs.columns:
            del obj_adata
            continue

        # Match cells by obs_names
        common_cells = adata_all.obs_names.intersection(obj_adata.obs_names)
        if len(common_cells) > 0:
            adata_all.obs.loc[common_cells, 'cell_type'] = obj_adata.obs.loc[common_cells, 'cell_type'].values
            if 'cell_type_confidence' in obj_adata.obs.columns:
                adata_all.obs.loc[common_cells, 'cell_type_confidence'] = obj_adata.obs.loc[common_cells, 'cell_type_confidence'].values
            if 'annotation_evidence' in obj_adata.obs.columns:
                adata_all.obs.loc[common_cells, 'annotation_evidence'] = obj_adata.obs.loc[common_cells, 'annotation_evidence'].values
            n_transferred += len(common_cells)
            print(f"    {obj_name}: transferred {len(common_cells):,} cell annotations")
        del obj_adata

    # De novo annotate remaining unassigned cells (e.g. GSE189916)
    unassigned_mask = adata_all.obs['cell_type'] == 'unassigned'
    n_unassigned = unassigned_mask.sum()
    print(f"  Cells with transferred annotations: {n_transferred:,}")
    print(f"  Cells needing de novo annotation: {n_unassigned:,}")

    if n_unassigned > 50:
        print("  Running de novo annotation on unassigned cells...")
        unassigned_adata = adata_all[unassigned_mask].copy()

        # Quick integration + clustering for unassigned cells
        if unassigned_adata.obs['study'].nunique() > 1:
            unassigned_adata = prepare_object(unassigned_adata)
            _, emb_key = run_scvi(unassigned_adata, model_dir=MODEL_DIR,
                                  label="all_cells_unassigned")
            optimize_clustering_resolution(unassigned_adata, emb_key,
                                           "all_cells", "unassigned")
        else:
            # Single study — just cluster on PCA
            sc.pp.pca(unassigned_adata)
            sc.pp.neighbors(unassigned_adata)
            sc.tl.leiden(unassigned_adata, resolution=0.5, key_added='leiden')
            emb_key = 'X_pca'

        annotate_clusters_de_novo(unassigned_adata, "all_cells", "unassigned")

        # Transfer back
        adata_all.obs.loc[unassigned_mask, 'cell_type'] = unassigned_adata.obs['cell_type'].values
        adata_all.obs.loc[unassigned_mask, 'cell_type_confidence'] = unassigned_adata.obs['cell_type_confidence'].values
        del unassigned_adata

    # Integration of the full all-cells object for visualization
    print("  Integrating all-cells for unified visualization...")
    adata_all = prepare_object(adata_all)

    # Run scVI on mesenchymal tier for visualization
    mes_mask = adata_all.obs['cell_class'] == 'mesenchymal'
    if mes_mask.sum() > 100:
        mes = adata_all[mes_mask].copy()
        run_scvi(mes, model_dir=MODEL_DIR, label="all_cells_mesenchymal")
        # Store embedding back
        adata_all.obsm['X_scvi_mesenchymal'] = np.full(
            (adata_all.shape[0], SCVI_PARAMS['n_latent']), np.nan
        )
        adata_all.obsm['X_scvi_mesenchymal'][np.where(mes_mask)[0]] = \
            mes.obsm[f'X_scvi_all_cells_mesenchymal']
        del mes

    # Compute integration metrics on full object
    metrics = compute_integration_metrics(adata_all, 'X_pca', batch_key='study')

    # UMAPs
    _plot_object_umaps(adata_all, "all_cells")

    # Save
    _save_checkpoint(adata_all, output_path, "final")
    print(f"  all_cells complete: {adata_all.shape[0]:,} cells")

    del adata_all
    gc.collect()
    return output_path, {'all_cells': metrics}


# ═════════════════════════════════════════════════════════════════════════════
# UMAP PLOTS
# ═════════════════════════════════════════════════════════════════════════════

def _plot_object_umaps(adata, object_name):
    """Generate UMAP panels for one object."""
    color_keys = ['study', 'leiden', 'cell_type', 'condition_harmonized', 'cell_class']
    color_keys = [c for c in color_keys if c in adata.obs.columns]
    n_cols = len(color_keys)

    if n_cols == 0 or 'X_umap' not in adata.obsm:
        return

    fig, axes = plt.subplots(1, n_cols, figsize=(5 * n_cols, 4.5))
    if n_cols == 1:
        axes = [axes]

    for ax, color in zip(axes, color_keys):
        try:
            sc.pl.umap(adata, color=color, ax=ax, show=False,
                       frameon=False, s=2, alpha=0.5,
                       title=f"{object_name} — {color}")
        except Exception:
            ax.set_title(f"{object_name} — {color} (failed)")

    plt.tight_layout()
    fig.savefig(RESULTS_DIR / f"umap_{object_name}.png",
                dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: umap_{object_name}.png")


# ═════════════════════════════════════════════════════════════════════════════
# INCLUSION SUMMARY & STUDY CAVEATS
# ═════════════════════════════════════════════════════════════════════════════

def generate_inclusion_summary():
    """Generate inclusion_summary.tsv and .html for manuscript supplement."""
    print("\n  Generating inclusion summary...")

    meta_path = META_DIR / "sample_metadata.tsv"
    meta = pd.read_csv(meta_path, sep='\t') if meta_path.exists() else None

    rows = []
    for obj_name in ["NP", "AF", "CEP", "all_cells"]:
        obj_path = INT_DIR / f"{obj_name}.h5ad"
        if not obj_path.exists():
            continue
        adata = sc.read_h5ad(obj_path, backed='r')

        for study in sorted(adata.obs['study'].unique()):
            study_mask = adata.obs['study'] == study
            n_cells = int(study_mask.sum())
            n_samples = adata.obs.loc[study_mask, 'sample_id'].nunique() if 'sample_id' in adata.obs.columns else 0

            # Get metadata
            first_author = ""
            year = ""
            compartment = ""
            conditions = ""
            platform = ""
            if meta is not None:
                study_meta = meta[meta['study_accession'] == study]
                if len(study_meta) > 0:
                    first_author = study_meta.iloc[0].get('first_author', '')
                    year = str(study_meta.iloc[0].get('year', ''))
                    compartment = ', '.join(study_meta['compartment'].unique())
                    conditions = ', '.join(study_meta['condition_harmonized'].dropna().unique())
                    platform = study_meta.iloc[0].get('sequencing_platform', '')

            rows.append({
                'object': obj_name,
                'accession': study,
                'first_author': first_author,
                'year': year,
                'n_samples': n_samples,
                'n_cells': n_cells,
                'compartment': compartment,
                'conditions': conditions,
                'platform': platform,
            })
        adata.file.close()

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS_DIR / "inclusion_summary.tsv", sep='\t', index=False)

    # HTML version
    html = df.to_html(index=False, classes='summary-table')
    html_page = f"""<!DOCTYPE html>
<html><head><title>Inclusion Summary</title>
<style>
body {{ font-family: Arial; max-width: 1200px; margin: 0 auto; padding: 20px; }}
.summary-table {{ border-collapse: collapse; width: 100%; }}
.summary-table th {{ background: #3498db; color: white; padding: 8px; }}
.summary-table td {{ border: 1px solid #ddd; padding: 6px; }}
.summary-table tr:nth-child(even) {{ background: #f2f2f2; }}
</style></head><body>
<h1>Study Inclusion Summary</h1>
<p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
{html}
</body></html>"""
    (RESULTS_DIR / "inclusion_summary.html").write_text(html_page)
    print(f"  Saved: inclusion_summary.tsv and .html")


def generate_study_caveats():
    """Generate study_caveats.tsv for manuscript supplement."""
    caveats = [
        {"study": "GSE165722", "first_author": "Tu 2022",
         "caveat": "BD Rhapsody platform (not 10x)",
         "impact": "Different capture efficiency, gene detection",
         "mitigation": "Platform-aware batch correction via scVI study-level batch key"},
        {"study": "GSE205535", "first_author": "Li 2022",
         "caveat": "BD Rhapsody platform; published corrigenda",
         "impact": "See above; potential data quality issues",
         "mitigation": "Monitor for outlier behavior in integration"},
        {"study": "CNP0002664", "first_author": "Han 2022",
         "caveat": "Singleron Matrix platform (not 10x)",
         "impact": "Different capture efficiency",
         "mitigation": "Same as above"},
        {"study": "GSE242443", "first_author": "Kuchynsky 2024",
         "caveat": "Culture-expanded CEP cells",
         "impact": "Culture alters cell states; may not reflect in vivo biology",
         "mitigation": "Caveat in all CEP results; compare with non-expanded CEP from GSE160756"},
        {"study": "GSE255768", "first_author": "Shi 2024",
         "caveat": "Degenerative endplate only; no healthy control",
         "impact": "Cannot do healthy vs. degenerated comparison for this study alone",
         "mitigation": "Healthy CEP baseline from GSE160756"},
        {"study": "GSE230809", "first_author": "Swahn 2024",
         "caveat": "All-male donors; age-disease confounded",
         "impact": "Cannot separate age from degeneration effects",
         "mitigation": "Note in interpretation; sex-specific effects cannot be assessed"},
        {"study": "GSE205535 NNP", "first_author": "Li 2022",
         "caveat": "11yo spinal cord injury, classified as 'healthy'",
         "impact": "Trauma confound",
         "mitigation": "Excluded from DE comparisons"},
        {"study": "GSE189916", "first_author": "Jiang 2022",
         "caveat": "Whole IVD (compartments not separated)",
         "impact": "Cannot assign cells to NP/AF/CEP",
         "mitigation": "Included only in all-cells object"},
    ]
    df = pd.DataFrame(caveats)
    df.to_csv(RESULTS_DIR / "study_caveats.tsv", sep='\t', index=False)
    print(f"  Saved: study_caveats.tsv")


# ═════════════════════════════════════════════════════════════════════════════
# INTEGRATION REPORT
# ═════════════════════════════════════════════════════════════════════════════

INTEGRATION_REPORT_TEMPLATE = """<!DOCTYPE html>
<html>
<head>
<title>Integration Report — Module 05</title>
<style>
  body { font-family: Arial, sans-serif; max-width: 1400px; margin: 0 auto; padding: 20px; }
  h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
  h2 { color: #34495e; margin-top: 30px; }
  table { border-collapse: collapse; width: 100%; margin: 10px 0; }
  th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
  th { background-color: #3498db; color: white; }
  tr:nth-child(even) { background-color: #f2f2f2; }
  img { max-width: 100%; height: auto; margin: 10px 0; border: 1px solid #ddd; }
  .stats-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0; }
  .stat-box { background: #ecf0f1; padding: 15px; border-radius: 5px; }
  .stat-box h3 { margin: 0 0 5px 0; color: #7f8c8d; font-size: 12px; }
  .stat-box p { margin: 0; font-size: 24px; font-weight: bold; color: #2c3e50; }
  .warning { color: #e74c3c; }
  .pass { color: #27ae60; }
  .secondary { color: #95a5a6; font-style: italic; }
  .validation li { margin: 4px 0; }
</style>
</head>
<body>
<h1>Module 05: Integration, Clustering, and De Novo Annotation</h1>
<p>Generated: {{ date }}</p>

{% for obj in objects %}
<h2>{{ obj.name }}{% if obj.secondary %} <span class="secondary">(secondary)</span>{% endif %}</h2>
{% if obj.umap_exists %}
<img src="umap_{{ obj.name }}.png" alt="{{ obj.name }} UMAP">
{% endif %}
<p>Cells: {{ obj.n_cells }} | Cell types: {{ obj.n_types }}</p>
{% endfor %}

<h2>Validation</h2>
<ul class="validation">
  {% for msg in validation_messages %}
  <li class="{{ 'pass' if msg.startswith('PASS') else 'warning' if msg.startswith('FAIL') or msg.startswith('WARNING') else '' }}">{{ msg }}</li>
  {% endfor %}
</ul>

<h2>Human Checkpoint Questions</h2>
<ol>
  <li>Does the clustering resolution capture biologically meaningful groups without over-splitting?</li>
  <li>Do the cell type annotations make biological sense?</li>
  <li>Are any clusters clearly batch-driven rather than biology-driven?</li>
  <li>For the mesenchymal continuum: are discrete labels appropriate, or should some clusters be merged?</li>
  <li>Are there unexpected cell types or missing expected types?</li>
  <li>Is the all-cells object consistent with the compartment-specific objects?</li>
  <li>Are the study caveats adequately documented?</li>
  <li>Should any annotations be revised before proceeding to differential analysis?</li>
</ol>

</body>
</html>"""


def generate_integration_report(validation_messages):
    """Generate HTML integration report."""
    print("\n  Generating integration report...")

    objects = []
    for obj_name in ["NP", "AF", "CEP", "all_cells"]:
        obj_path = INT_DIR / f"{obj_name}.h5ad"
        if obj_path.exists():
            adata = sc.read_h5ad(obj_path, backed='r')
            objects.append({
                'name': obj_name,
                'n_cells': f"{adata.shape[0]:,}",
                'n_types': adata.obs['cell_type'].nunique() if 'cell_type' in adata.obs.columns else 0,
                'umap_exists': (RESULTS_DIR / f"umap_{obj_name}.png").exists(),
                'secondary': obj_name == 'all_cells',
            })
            adata.file.close()

    template = Template(INTEGRATION_REPORT_TEMPLATE)
    html = template.render(
        date=datetime.now().strftime("%Y-%m-%d %H:%M"),
        objects=objects,
        validation_messages=validation_messages,
    )
    (RESULTS_DIR / "integration_report.html").write_text(html)
    print(f"  Report: {RESULTS_DIR / 'integration_report.html'}")


# ═════════════════════════════════════════════════════════════════════════════
# VALIDATION
# ═════════════════════════════════════════════════════════════════════════════

def validate_integration():
    """Validate integration outputs. Returns (pass, messages)."""
    print("\n" + "=" * 60)
    print("Validating integration outputs")
    print("=" * 60)

    messages = []
    all_pass = True

    # Check all four output files
    for obj_name in ["NP", "AF", "CEP", "all_cells"]:
        path = INT_DIR / f"{obj_name}.h5ad"
        if path.exists():
            messages.append(f"PASS: {obj_name} output exists")
            adata = sc.read_h5ad(path, backed='r')

            # Check cell_type column
            if 'cell_type' in adata.obs.columns:
                n_unassigned = (adata.obs['cell_type'] == 'unassigned').sum()
                pct_unassigned = n_unassigned / adata.shape[0] * 100
                if pct_unassigned < 10:
                    messages.append(f"PASS: {obj_name} unassigned = {pct_unassigned:.1f}% (< 10%)")
                else:
                    messages.append(f"WARNING: {obj_name} unassigned = {pct_unassigned:.1f}% (>= 10%)")
            else:
                messages.append(f"FAIL: {obj_name} missing cell_type column")
                all_pass = False

            # Blob check
            if 'leiden' in adata.obs.columns:
                n_clusters = adata.obs['leiden'].nunique()
                if n_clusters > 1:
                    messages.append(f"PASS: {obj_name} has {n_clusters} clusters (no blob)")
                else:
                    messages.append(f"WARNING: {obj_name} has only {n_clusters} cluster (blob?)")

            # ARI check (study shouldn't perfectly predict clusters)
            if 'leiden' in adata.obs.columns and 'study' in adata.obs.columns:
                try:
                    from sklearn.metrics import adjusted_rand_score
                    ari = adjusted_rand_score(
                        adata.obs['study'].values[:10000],
                        adata.obs['leiden'].values[:10000]
                    )
                    if ari < 1.0:
                        messages.append(f"PASS: {obj_name} study-cluster ARI = {ari:.3f} (< 1.0)")
                    else:
                        messages.append(f"WARNING: {obj_name} study perfectly predicts clusters")
                except Exception:
                    pass

            adata.file.close()
        else:
            messages.append(f"FAIL: {obj_name} output missing")
            all_pass = False

    # Check supplementary outputs
    for fname in ["inclusion_summary.tsv", "study_caveats.tsv"]:
        if (RESULTS_DIR / fname).exists():
            messages.append(f"PASS: {fname} exists")
        else:
            messages.append(f"FAIL: {fname} missing")
            all_pass = False

    # Check report
    if (RESULTS_DIR / "integration_report.html").exists():
        messages.append("PASS: Integration report exists")
    else:
        messages.append("FAIL: Integration report missing")
        all_pass = False

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
    force = '--force' in args

    # Parse --object
    object_filter = None
    for i, a in enumerate(args):
        if a == '--object' and i + 1 < len(args):
            object_filter = args[i + 1]

    if validate_only:
        passed, _ = validate_integration()
        sys.exit(0 if passed else 1)

    print("=" * 60)
    print("Module 05: Integration, Clustering, and De Novo Annotation")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    all_metrics = {}

    # Process compartment-specific objects (primary)
    primary_objects = ["NP", "AF", "CEP"]
    if object_filter:
        if object_filter in primary_objects:
            primary_objects = [object_filter]
        elif object_filter == "all_cells":
            primary_objects = []
        else:
            print(f"Unknown object: {object_filter}")
            sys.exit(1)

    for obj_name in primary_objects:
        _, metrics = process_object(obj_name, force=force)
        all_metrics[obj_name] = metrics

    # Process all-cells (secondary) — after primary objects are done
    if object_filter is None or object_filter == "all_cells":
        _, metrics = process_all_cells_secondary(force=force)
        all_metrics['all_cells'] = metrics

    # Save metrics
    metrics_rows = []
    for obj_name, obj_metrics in all_metrics.items():
        for tier, tier_metrics in obj_metrics.items():
            if isinstance(tier_metrics, dict):
                row = {'object': obj_name, 'tier': tier}
                row.update(tier_metrics)
                metrics_rows.append(row)
    if metrics_rows:
        pd.DataFrame(metrics_rows).to_csv(
            RESULTS_DIR / "integration_metrics.tsv", sep='\t', index=False
        )

    # Generate supplementary tables
    generate_inclusion_summary()
    generate_study_caveats()

    # Validation
    passed, val_messages = validate_integration()

    # Report
    generate_integration_report(val_messages)

    print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Overall: {'PASSED' if passed else 'FAILED'}")


if __name__ == "__main__":
    main()
