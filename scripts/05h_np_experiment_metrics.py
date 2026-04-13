#!/usr/bin/env python3
"""NP Integration Quality Experiment — Expanded Metrics.

Computes an expanded suite of integration quality metrics across the three
experimental runs (tiered_v5, flat_v4, tiered_v4) plus the existing v5 flat
CCA baseline. Reads bridge-exported files (embedding, metadata, counts).

Metrics:
  Batch removal:    iLISI, batch_ASW
  Bio preservation: cLISI (on coarse_label), bio_ASW (on coarse_label)
  Condition signal: condition_ASW, condition_LISI
  Cluster agreement: ARI, NMI (Leiden vs coarse_label)
  Continuum signal: marker gene variance ratios (COL2A1, ACAN, SOX9, COL1A1)

Usage:
    python3 scripts/05h_np_experiment_metrics.py
    python3 scripts/05h_np_experiment_metrics.py --run tiered_v5
"""

import sys
import warnings
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread
from scipy.sparse import csc_matrix
from sklearn.metrics import silhouette_score

import scib_metrics
from scib_metrics.nearest_neighbors import pynndescent, NeighborsResults

warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)

# ── Paths ────────────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent.parent
EXP_DIR = BASE / "data" / "integrated" / "np_experiment"
BASELINE_DIR = BASE / "data" / "integrated" / "cca" / "bridge_export" / "NP"
RESULTS_DIR = BASE / "results" / "integration" / "np_experiment"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ── Parameters ───────────────────────────────────────────────────────────
MAX_CELLS_LISI = 50_000    # subsample for LISI (KNN-based)
MAX_CELLS_ASW = 30_000     # subsample for silhouette (pairwise distances)
KNN_NEIGHBORS = 90         # for LISI computation (perplexity ~30)
KNN_JOBS = 16
MARKER_GENES = ["COL2A1", "ACAN", "SOX9", "COL1A1"]
RNG_SEED = 42

# ── Run discovery ────────────────────────────────────────────────────────
# Each run has one or more scopes (tiers)
RUN_SCOPES = {
    "baseline_flat_v5": [
        {"scope": "all", "bridge_dir": BASELINE_DIR}
    ],
    "tiered_v5": [
        {"scope": "mesenchymal", "bridge_dir": EXP_DIR / "tiered_v5" / "mesenchymal"},
        {"scope": "non_mesenchymal", "bridge_dir": EXP_DIR / "tiered_v5" / "non_mesenchymal"},
    ],
    "flat_v4": [
        {"scope": "all", "bridge_dir": EXP_DIR / "flat_v4" / "all"},
    ],
    "tiered_v4": [
        {"scope": "mesenchymal", "bridge_dir": EXP_DIR / "tiered_v4" / "mesenchymal"},
        {"scope": "non_mesenchymal", "bridge_dir": EXP_DIR / "tiered_v4" / "non_mesenchymal"},
    ],
}


# ═════════════════════════════════════════════════════════════════════════
# DATA LOADING
# ═════════════════════════════════════════════════════════════════════════

def _detect_embedding_file(bridge_dir):
    """Find the best embedding file: prefer integrated.cca, then pca."""
    for name in ["embedding_integrated.cca.csv.gz", "embedding_pca.csv.gz"]:
        path = bridge_dir / name
        if path.exists():
            return path
    # Fallback: any embedding file
    emb_files = sorted(bridge_dir.glob("embedding_*.csv.gz"))
    if emb_files:
        return emb_files[0]
    return None


def load_bridge_data(bridge_dir):
    """Load embedding, metadata, and optionally counts from bridge export.

    Returns:
        embedding: np.ndarray (n_cells, n_dims)
        metadata: pd.DataFrame
        counts: scipy sparse matrix or None (cells x genes)
        gene_names: list or None
    """
    bridge_dir = Path(bridge_dir)

    # Embedding
    emb_path = _detect_embedding_file(bridge_dir)
    if emb_path is None:
        raise FileNotFoundError(f"No embedding file found in {bridge_dir}")
    print(f"    Loading embedding: {emb_path.name}")
    emb_df = pd.read_csv(emb_path, index_col=0)
    embedding = emb_df.values.astype(np.float32)

    # Metadata
    meta_path = bridge_dir / "metadata.csv.gz"
    if not meta_path.exists():
        meta_path = bridge_dir / "metadata.csv"
    print(f"    Loading metadata: {meta_path.name}")
    metadata = pd.read_csv(meta_path)

    # Align metadata with embedding rows
    if len(metadata) != embedding.shape[0]:
        print(f"    WARNING: metadata ({len(metadata)}) != embedding ({embedding.shape[0]}) rows")
        n = min(len(metadata), embedding.shape[0])
        metadata = metadata.iloc[:n]
        embedding = embedding[:n]

    # Counts (optional, for marker variance)
    counts = None
    gene_names = None
    mtx_path = bridge_dir / "counts.mtx.gz"
    if not mtx_path.exists():
        mtx_path = bridge_dir / "counts.mtx"
    genes_path = bridge_dir / "genes.csv"

    if mtx_path.exists() and genes_path.exists():
        print(f"    Loading counts matrix...")
        raw = mmread(str(mtx_path))
        # Seurat exports genes x cells; transpose to cells x genes
        if raw.shape[0] != embedding.shape[0]:
            raw = raw.T
        counts = csc_matrix(raw)
        genes_df = pd.read_csv(genes_path)
        gene_names = genes_df['gene'].tolist()
        print(f"    Counts: {counts.shape[0]} cells x {counts.shape[1]} genes")

    print(f"    Loaded: {embedding.shape[0]} cells, {embedding.shape[1]} dims")
    return embedding, metadata, counts, gene_names


# ═════════════════════════════════════════════════════════════════════════
# SUBSAMPLING
# ═════════════════════════════════════════════════════════════════════════

def _subsample(embedding, metadata, max_cells):
    """Subsample embedding and metadata with fixed seed."""
    n = embedding.shape[0]
    if n <= max_cells:
        return embedding, metadata, np.arange(n)
    rng = np.random.RandomState(RNG_SEED)
    idx = rng.choice(n, max_cells, replace=False)
    idx.sort()
    return embedding[idx], metadata.iloc[idx].reset_index(drop=True), idx


# ═════════════════════════════════════════════════════════════════════════
# METRIC COMPUTATION
# ═════════════════════════════════════════════════════════════════════════

def compute_batch_metrics(embedding, metadata):
    """Compute iLISI and batch_ASW."""
    metrics = {}
    batch = metadata['study'].values
    unique_batches = len(set(batch))

    if unique_batches < 2:
        print("    WARNING: <2 batches, skipping batch metrics")
        return {'iLISI': np.nan, 'batch_ASW': np.nan}

    # Subsample for LISI
    emb_lisi, meta_lisi, _ = _subsample(embedding, metadata, MAX_CELLS_LISI)
    batch_lisi = meta_lisi['study'].values

    # iLISI
    try:
        print("    Computing KNN for LISI...")
        nn = pynndescent(emb_lisi, n_neighbors=KNN_NEIGHBORS,
                         random_state=RNG_SEED, n_jobs=KNN_JOBS)
        metrics['iLISI'] = float(scib_metrics.ilisi_knn(nn, batch_lisi, scale=True))
        print(f"    iLISI: {metrics['iLISI']:.4f}")
    except Exception as e:
        print(f"    WARNING: iLISI failed: {e}")
        metrics['iLISI'] = np.nan

    # batch_ASW (using scib-metrics: per-label batch silhouette)
    try:
        emb_asw, meta_asw, _ = _subsample(embedding, metadata, MAX_CELLS_ASW)
        labels = meta_asw['coarse_label'].values if 'coarse_label' in meta_asw.columns else None
        batch_asw = meta_asw['study'].values
        if labels is not None and len(set(labels)) > 1:
            metrics['batch_ASW'] = float(scib_metrics.silhouette_batch(
                emb_asw, labels, batch_asw, rescale=True
            ))
        else:
            # Fallback: raw silhouette on batch
            metrics['batch_ASW'] = float(silhouette_score(
                emb_asw, batch_asw,
                sample_size=min(5000, len(emb_asw)),
                random_state=RNG_SEED
            ))
        print(f"    batch_ASW: {metrics['batch_ASW']:.4f}")
    except Exception as e:
        print(f"    WARNING: batch_ASW failed: {e}")
        metrics['batch_ASW'] = np.nan

    return metrics


def compute_bio_metrics(embedding, metadata):
    """Compute cLISI and bio_ASW on coarse_label."""
    metrics = {}

    if 'coarse_label' not in metadata.columns:
        print("    WARNING: coarse_label not in metadata, skipping bio metrics")
        return {'cLISI': np.nan, 'bio_ASW': np.nan}

    labels = metadata['coarse_label'].values
    valid = pd.notna(labels)
    if valid.sum() < 100 or len(set(labels[valid])) < 2:
        print("    WARNING: insufficient coarse_label variation, skipping bio metrics")
        return {'cLISI': np.nan, 'bio_ASW': np.nan}

    # cLISI (subsample for KNN)
    try:
        emb_lisi, meta_lisi, _ = _subsample(embedding, metadata, MAX_CELLS_LISI)
        labels_lisi = meta_lisi['coarse_label'].values
        nn = pynndescent(emb_lisi, n_neighbors=KNN_NEIGHBORS,
                         random_state=RNG_SEED, n_jobs=KNN_JOBS)
        metrics['cLISI'] = float(scib_metrics.clisi_knn(nn, labels_lisi, scale=True))
        print(f"    cLISI: {metrics['cLISI']:.4f} (lower = better cell-type purity)")
    except Exception as e:
        print(f"    WARNING: cLISI failed: {e}")
        metrics['cLISI'] = np.nan

    # bio_ASW
    try:
        emb_asw, meta_asw, _ = _subsample(embedding, metadata, MAX_CELLS_ASW)
        labels_asw = meta_asw['coarse_label'].values
        metrics['bio_ASW'] = float(scib_metrics.silhouette_label(
            emb_asw, labels_asw, rescale=True
        ))
        print(f"    bio_ASW: {metrics['bio_ASW']:.4f} (higher = better bio separation)")
    except Exception as e:
        print(f"    WARNING: bio_ASW failed: {e}")
        metrics['bio_ASW'] = np.nan

    return metrics


def compute_condition_metrics(embedding, metadata):
    """Compute condition_ASW and condition_LISI."""
    metrics = {}

    if 'condition_harmonized' not in metadata.columns:
        return {'condition_ASW': np.nan, 'condition_LISI': np.nan}

    cond = metadata['condition_harmonized'].values
    valid = pd.notna(cond)
    if valid.sum() < 100 or len(set(cond[valid])) < 2:
        return {'condition_ASW': np.nan, 'condition_LISI': np.nan}

    # condition_ASW (raw sklearn silhouette — not rescaled, for interpretability)
    try:
        emb_asw, meta_asw, _ = _subsample(embedding, metadata, MAX_CELLS_ASW)
        cond_asw = meta_asw['condition_harmonized'].values
        valid_asw = pd.notna(cond_asw)
        metrics['condition_ASW'] = float(silhouette_score(
            emb_asw[valid_asw], cond_asw[valid_asw],
            sample_size=min(10_000, int(valid_asw.sum())),
            random_state=RNG_SEED
        ))
        print(f"    condition_ASW: {metrics['condition_ASW']:.4f} (positive = conditions separable)")
    except Exception as e:
        print(f"    WARNING: condition_ASW failed: {e}")
        metrics['condition_ASW'] = np.nan

    # condition_LISI
    try:
        emb_lisi, meta_lisi, _ = _subsample(embedding, metadata, MAX_CELLS_LISI)
        cond_lisi = meta_lisi['condition_harmonized'].values
        nn = pynndescent(emb_lisi, n_neighbors=KNN_NEIGHBORS,
                         random_state=RNG_SEED, n_jobs=KNN_JOBS)
        lisi_per_cell = scib_metrics.lisi_knn(nn, cond_lisi)
        metrics['condition_LISI'] = float(np.mean(lisi_per_cell))
        print(f"    condition_LISI: {metrics['condition_LISI']:.4f} (lower = better condition preservation)")
    except Exception as e:
        print(f"    WARNING: condition_LISI failed: {e}")
        metrics['condition_LISI'] = np.nan

    return metrics


def compute_cluster_agreement(embedding, metadata):
    """Compute ARI and NMI: Leiden clusters vs coarse_label."""
    metrics = {}

    if 'coarse_label' not in metadata.columns:
        return {'ARI': np.nan, 'NMI': np.nan}

    labels = metadata['coarse_label'].values
    if len(set(labels)) < 2:
        return {'ARI': np.nan, 'NMI': np.nan}

    try:
        # Use scib-metrics Leiden-based NMI/ARI
        emb_sub, meta_sub, _ = _subsample(embedding, metadata, MAX_CELLS_LISI)
        labels_sub = meta_sub['coarse_label'].values
        nn = pynndescent(emb_sub, n_neighbors=KNN_NEIGHBORS,
                         random_state=RNG_SEED, n_jobs=KNN_JOBS)
        result = scib_metrics.nmi_ari_cluster_labels_leiden(
            nn, labels_sub, optimize_resolution=True, seed=RNG_SEED
        )
        metrics['NMI'] = float(result['nmi'])
        metrics['ARI'] = float(result['ari'])
        print(f"    NMI: {metrics['NMI']:.4f}, ARI: {metrics['ARI']:.4f}")
    except Exception as e:
        print(f"    WARNING: NMI/ARI failed: {e}")
        metrics['NMI'] = np.nan
        metrics['ARI'] = np.nan

    return metrics


def compute_marker_variance(counts, gene_names, embedding, metadata):
    """Compute marker gene within-cluster variance ratios.

    For each marker, clusters the embedding (Leiden res=0.5) and computes:
        ratio = mean(within-cluster variance) / total variance

    Higher ratio = more heterogeneity preserved within clusters (good).
    """
    metrics = {}

    if counts is None or gene_names is None:
        print("    WARNING: counts not available, skipping marker variance")
        return {f"var_ratio_{g}": np.nan for g in MARKER_GENES}

    # Build a minimal AnnData for Leiden clustering
    adata = sc.AnnData(X=counts)
    adata.obsm['X_emb'] = embedding

    # Leiden clustering
    sc.pp.neighbors(adata, use_rep='X_emb', n_neighbors=15)
    sc.tl.leiden(adata, resolution=0.5, flavor='igraph', n_iterations=2)
    clusters = adata.obs['leiden'].values

    # Normalize counts for variance computation
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    for gene in MARKER_GENES:
        key = f"var_ratio_{gene}"
        if gene not in gene_names:
            # Try with dashes (Seurat replaces _ with -)
            gene_alt = gene.replace("_", "-")
            if gene_alt in gene_names:
                gene_idx = gene_names.index(gene_alt)
            else:
                print(f"    WARNING: {gene} not found in gene list")
                metrics[key] = np.nan
                continue
        else:
            gene_idx = gene_names.index(gene)

        expr = np.asarray(adata.X[:, gene_idx].todense()).flatten()
        total_var = np.var(expr)
        if total_var < 1e-10:
            metrics[key] = np.nan
            continue

        # Weighted mean of within-cluster variances
        cluster_vars = []
        cluster_sizes = []
        for cl in np.unique(clusters):
            mask = clusters == cl
            n_cl = mask.sum()
            if n_cl < 5:
                continue
            cluster_vars.append(np.var(expr[mask]))
            cluster_sizes.append(n_cl)

        if len(cluster_vars) == 0:
            metrics[key] = np.nan
            continue

        cluster_vars = np.array(cluster_vars)
        cluster_sizes = np.array(cluster_sizes)
        weighted_within_var = np.average(cluster_vars, weights=cluster_sizes)
        ratio = weighted_within_var / total_var
        metrics[key] = float(ratio)
        print(f"    var_ratio_{gene}: {ratio:.4f}")

    del adata
    return metrics


# ═════════════════════════════════════════════════════════════════════════
# MAIN METRICS PIPELINE
# ═════════════════════════════════════════════════════════════════════════

def compute_all_metrics(bridge_dir, run_name, scope_name):
    """Compute full metrics suite for one run/scope."""
    print(f"\n{'='*60}")
    print(f"  {run_name} / {scope_name}")
    print(f"{'='*60}")

    bridge_dir = Path(bridge_dir)
    if not bridge_dir.exists():
        print(f"  SKIP: {bridge_dir} does not exist")
        return None

    # Check for required files
    emb_path = _detect_embedding_file(bridge_dir)
    if emb_path is None:
        print(f"  SKIP: no embedding file in {bridge_dir}")
        return None

    embedding, metadata, counts, gene_names = load_bridge_data(bridge_dir)

    result = {
        'run': run_name,
        'scope': scope_name,
        'n_cells': embedding.shape[0],
        'n_dims': embedding.shape[1],
    }

    # Batch metrics
    print("\n  --- Batch removal metrics ---")
    result.update(compute_batch_metrics(embedding, metadata))

    # Bio preservation metrics
    print("\n  --- Biological preservation metrics ---")
    result.update(compute_bio_metrics(embedding, metadata))

    # Condition metrics
    print("\n  --- Condition signal metrics ---")
    result.update(compute_condition_metrics(embedding, metadata))

    # Cluster agreement
    print("\n  --- Cluster agreement metrics ---")
    result.update(compute_cluster_agreement(embedding, metadata))

    # Marker variance (only if counts available)
    print("\n  --- Marker gene variance ---")
    result.update(compute_marker_variance(counts, gene_names, embedding, metadata))

    return result


# ═════════════════════════════════════════════════════════════════════════
# COMPARISON TABLE + HTML REPORT
# ═════════════════════════════════════════════════════════════════════════

HTML_TEMPLATE = """<!DOCTYPE html>
<html>
<head>
<title>NP Integration Experiment — Metrics Comparison</title>
<style>
  body {{ font-family: Arial, sans-serif; max-width: 1600px; margin: 0 auto; padding: 20px; }}
  h1 {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
  h2 {{ color: #34495e; margin-top: 30px; }}
  table {{ border-collapse: collapse; width: 100%; margin: 15px 0; }}
  th, td {{ border: 1px solid #ddd; padding: 8px 12px; text-align: right; }}
  th {{ background-color: #3498db; color: white; text-align: center; }}
  td:first-child, td:nth-child(2) {{ text-align: left; }}
  tr:nth-child(even) {{ background-color: #f9f9f9; }}
  tr.baseline {{ background-color: #fff3cd; font-weight: bold; }}
  .better {{ color: #27ae60; font-weight: bold; }}
  .worse {{ color: #e74c3c; }}
  .neutral {{ color: #7f8c8d; }}
  .summary {{ background: #eaf2f8; padding: 15px; border-radius: 5px; margin: 20px 0; }}
  .metric-desc {{ font-size: 0.85em; color: #555; }}
</style>
</head>
<body>
<h1>NP Integration Quality Experiment</h1>
<p>Generated: {date}</p>

<div class="summary">
<h3>Experiment</h3>
<p>Comparing three integration strategies for NP to address over-integration in the flat v5 CCA baseline.
The baseline (highlighted) has iLISI=3.68 and condition_ASW=-0.156, indicating aggressive batch mixing
that erases condition signal.</p>
<ul>
  <li><strong>tiered_v5:</strong> Seurat v5 IntegrateLayers + mesenchymal/non-mes split</li>
  <li><strong>flat_v4:</strong> Seurat v4 SCT + FindIntegrationAnchors (no split)</li>
  <li><strong>tiered_v4:</strong> Seurat v4 SCT + FindIntegrationAnchors + mesenchymal/non-mes split</li>
</ul>
</div>

<h2>Metric Interpretation</h2>
<table>
<tr><th>Metric</th><th>Good direction</th><th>Description</th></tr>
<tr><td>iLISI</td><td>Higher</td><td>Batch mixing in neighborhoods (1=no mixing, n_batches=perfect)</td></tr>
<tr><td>batch_ASW</td><td>Higher (rescaled)</td><td>Per-label batch silhouette (0=poor, 1=perfect mixing)</td></tr>
<tr><td>cLISI</td><td>Lower</td><td>Cell-type mixing in neighborhoods (0=pure, 1=fully mixed)</td></tr>
<tr><td>bio_ASW</td><td>Higher</td><td>Cell-type silhouette (0=poor, 1=perfect separation)</td></tr>
<tr><td>condition_ASW</td><td>More positive</td><td>Condition separability (-1 to 1; negative = over-mixed)</td></tr>
<tr><td>condition_LISI</td><td>Lower</td><td>Condition mixing in neighborhoods (1=preserved, high=erased)</td></tr>
<tr><td>NMI</td><td>Higher</td><td>Leiden vs coarse_label agreement (0 to 1)</td></tr>
<tr><td>ARI</td><td>Higher</td><td>Leiden vs coarse_label agreement (-1 to 1)</td></tr>
<tr><td>var_ratio_*</td><td>Higher</td><td>Within-cluster marker variance / total variance (heterogeneity preserved)</td></tr>
</table>

<h2>Results</h2>
{table_html}

</body>
</html>"""


def build_comparison_table(all_results):
    """Build comparison TSV and HTML report."""
    df = pd.DataFrame(all_results)

    # Order columns
    id_cols = ['run', 'scope', 'n_cells', 'n_dims']
    metric_cols = ['iLISI', 'batch_ASW', 'cLISI', 'bio_ASW',
                   'condition_ASW', 'condition_LISI', 'NMI', 'ARI']
    var_cols = [c for c in df.columns if c.startswith('var_ratio_')]
    ordered = id_cols + metric_cols + sorted(var_cols)
    ordered = [c for c in ordered if c in df.columns]
    df = df[ordered]

    # Save TSV
    tsv_path = RESULTS_DIR / "comparison_table.tsv"
    df.to_csv(tsv_path, sep='\t', index=False, float_format='%.4f')
    print(f"\n  Saved: {tsv_path}")

    # Build HTML
    # Format numeric columns
    display_df = df.copy()
    for col in display_df.columns:
        if col in id_cols:
            continue
        display_df[col] = display_df[col].apply(
            lambda x: f"{x:.4f}" if pd.notna(x) and isinstance(x, (int, float, np.floating)) else str(x)
        )

    # Mark baseline row
    table_html = display_df.to_html(index=False, classes='results', na_rep='N/A', escape=False)
    # Add baseline class to baseline rows
    table_html = table_html.replace(
        '<td>baseline_flat_v5</td>',
        '<td class="baseline">baseline_flat_v5</td>'
    )

    html = HTML_TEMPLATE.format(
        date=datetime.now().strftime("%Y-%m-%d %H:%M"),
        table_html=table_html,
    )

    html_path = RESULTS_DIR / "comparison_table.html"
    html_path.write_text(html)
    print(f"  Saved: {html_path}")

    # Print summary to console
    print(f"\n{'='*80}")
    print("COMPARISON TABLE (mesenchymal / all scopes only)")
    print(f"{'='*80}")
    focus = df[df['scope'].isin(['all', 'mesenchymal'])]
    print(focus.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    return df


# ═════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════

def main():
    import argparse
    parser = argparse.ArgumentParser(description="NP Integration Experiment Metrics")
    parser.add_argument('--run', type=str, default=None,
                        choices=['baseline_flat_v5', 'tiered_v5', 'flat_v4', 'tiered_v4'],
                        help="Compute metrics for a single run (default: all)")
    args = parser.parse_args()

    print("=" * 60)
    print("NP Integration Quality Experiment — Expanded Metrics")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    # Determine which runs to process
    if args.run:
        runs_to_process = {args.run: RUN_SCOPES[args.run]}
    else:
        runs_to_process = RUN_SCOPES

    all_results = []

    for run_name, scopes in runs_to_process.items():
        for scope_info in scopes:
            scope_name = scope_info['scope']
            bridge_dir = scope_info['bridge_dir']

            result = compute_all_metrics(bridge_dir, run_name, scope_name)
            if result is not None:
                all_results.append(result)

    if not all_results:
        print("\nNo results computed. Ensure bridge export files exist.")
        sys.exit(1)

    build_comparison_table(all_results)
    print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    main()
