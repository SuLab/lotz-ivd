#!/usr/bin/env python3
"""Tiered v4 Integration Metrics — AF, CEP, all_cells.

Generalization of 05h (NP) to the remaining compartments. Reads bridge
exports for the production v5 CCA baseline and the new tiered_v4 mes /
non-mes outputs, then writes a per-compartment comparison_table.{tsv,html}
under results/integration/{compartment}_experiment/.

Metrics, parameters, and reporting format mirror 05h exactly so the
notebook and downstream interpretation can reuse the same machinery.

Usage:
    python3 scripts/05l_compartment_metrics.py                # all 3
    python3 scripts/05l_compartment_metrics.py --compartment AF
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
from scib_metrics.nearest_neighbors import pynndescent

warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)

# ── Paths ────────────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent.parent
INTEG = BASE / "data" / "integrated"
BASELINE_BASE = INTEG / "cca" / "bridge_export"
RESULTS_BASE = BASE / "results" / "integration"

# ── Parameters (same as 05h) ─────────────────────────────────────────────
MAX_CELLS_LISI = 50_000
MAX_CELLS_ASW = 30_000
KNN_NEIGHBORS = 90
KNN_JOBS = 16
MARKER_GENES = ["COL2A1", "ACAN", "SOX9", "COL1A1"]
RNG_SEED = 42

# ── Per-compartment run definitions ──────────────────────────────────────
COMPARTMENT_RUNS = {
    "AF": {
        "exp_dir_name": "af_experiment",
        "bridge_subdir": "AF",
        "tiers_present": ["mesenchymal"],
    },
    "CEP": {
        "exp_dir_name": "cep_experiment",
        "bridge_subdir": "CEP",
        "tiers_present": ["mesenchymal", "non_mesenchymal"],
    },
    "all_cells": {
        "exp_dir_name": "all_cells_experiment",
        "bridge_subdir": "all_cells",
        "tiers_present": ["mesenchymal", "non_mesenchymal"],
    },
}


def build_run_scopes(compartment):
    """Return RUN_SCOPES dict for a given compartment."""
    cfg = COMPARTMENT_RUNS[compartment]
    exp_dir = INTEG / cfg["exp_dir_name"]
    baseline_dir = BASELINE_BASE / cfg["bridge_subdir"]

    runs = {
        "baseline_flat_v5": [
            {"scope": "all", "bridge_dir": baseline_dir}
        ],
        "tiered_v4": [
            {"scope": tier,
             "bridge_dir": exp_dir / "tiered_v4" / tier}
            for tier in cfg["tiers_present"]
        ],
    }
    return runs


# ═════════════════════════════════════════════════════════════════════════
# DATA LOADING (unchanged from 05h)
# ═════════════════════════════════════════════════════════════════════════

def _detect_embedding_file(bridge_dir):
    for name in ["embedding_integrated.cca.csv.gz", "embedding_pca.csv.gz"]:
        path = bridge_dir / name
        if path.exists():
            return path
    emb_files = sorted(bridge_dir.glob("embedding_*.csv.gz"))
    if emb_files:
        return emb_files[0]
    return None


def load_bridge_data(bridge_dir):
    bridge_dir = Path(bridge_dir)

    emb_path = _detect_embedding_file(bridge_dir)
    if emb_path is None:
        raise FileNotFoundError(f"No embedding file found in {bridge_dir}")
    print(f"    Loading embedding: {emb_path.name}")
    emb_df = pd.read_csv(emb_path, index_col=0)
    embedding = emb_df.values.astype(np.float32)

    meta_path = bridge_dir / "metadata.csv.gz"
    if not meta_path.exists():
        meta_path = bridge_dir / "metadata.csv"
    print(f"    Loading metadata: {meta_path.name}")
    metadata = pd.read_csv(meta_path)

    if len(metadata) != embedding.shape[0]:
        print(f"    WARNING: metadata ({len(metadata)}) != embedding "
              f"({embedding.shape[0]}) rows")
        n = min(len(metadata), embedding.shape[0])
        metadata = metadata.iloc[:n]
        embedding = embedding[:n]

    counts = None
    gene_names = None
    mtx_path = bridge_dir / "counts.mtx.gz"
    if not mtx_path.exists():
        mtx_path = bridge_dir / "counts.mtx"
    genes_path = bridge_dir / "genes.csv"

    if mtx_path.exists() and genes_path.exists():
        print(f"    Loading counts matrix...")
        raw = mmread(str(mtx_path))
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
    n = embedding.shape[0]
    if n <= max_cells:
        return embedding, metadata, np.arange(n)
    rng = np.random.RandomState(RNG_SEED)
    idx = rng.choice(n, max_cells, replace=False)
    idx.sort()
    return embedding[idx], metadata.iloc[idx].reset_index(drop=True), idx


def _knn_neighbors(n_cells):
    """Pick a safe k for KNN given the cell count (KNN_NEIGHBORS, capped)."""
    return max(5, min(KNN_NEIGHBORS, n_cells - 1))


# ═════════════════════════════════════════════════════════════════════════
# METRICS
# ═════════════════════════════════════════════════════════════════════════

def compute_batch_metrics(embedding, metadata):
    metrics = {}
    batch = metadata['study'].values
    unique_batches = len(set(batch))

    if unique_batches < 2:
        print("    WARNING: <2 batches, skipping batch metrics")
        return {'iLISI': np.nan, 'batch_ASW': np.nan}

    emb_lisi, meta_lisi, _ = _subsample(embedding, metadata, MAX_CELLS_LISI)
    batch_lisi = meta_lisi['study'].values
    k_lisi = _knn_neighbors(emb_lisi.shape[0])

    try:
        print("    Computing KNN for LISI...")
        nn = pynndescent(emb_lisi, n_neighbors=k_lisi,
                         random_state=RNG_SEED, n_jobs=KNN_JOBS)
        metrics['iLISI'] = float(scib_metrics.ilisi_knn(nn, batch_lisi, scale=True))
        print(f"    iLISI: {metrics['iLISI']:.4f}")
    except Exception as e:
        print(f"    WARNING: iLISI failed: {e}")
        metrics['iLISI'] = np.nan

    try:
        emb_asw, meta_asw, _ = _subsample(embedding, metadata, MAX_CELLS_ASW)
        labels = meta_asw['coarse_label'].values if 'coarse_label' in meta_asw.columns else None
        batch_asw = meta_asw['study'].values
        if labels is not None and len(set(labels)) > 1:
            metrics['batch_ASW'] = float(scib_metrics.silhouette_batch(
                emb_asw, labels, batch_asw, rescale=True
            ))
        else:
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
    metrics = {}
    if 'coarse_label' not in metadata.columns:
        return {'cLISI': np.nan, 'bio_ASW': np.nan}

    labels = metadata['coarse_label'].values
    valid = pd.notna(labels)
    if valid.sum() < 100 or len(set(labels[valid])) < 2:
        return {'cLISI': np.nan, 'bio_ASW': np.nan}

    try:
        emb_lisi, meta_lisi, _ = _subsample(embedding, metadata, MAX_CELLS_LISI)
        labels_lisi = meta_lisi['coarse_label'].values
        k_lisi = _knn_neighbors(emb_lisi.shape[0])
        nn = pynndescent(emb_lisi, n_neighbors=k_lisi,
                         random_state=RNG_SEED, n_jobs=KNN_JOBS)
        metrics['cLISI'] = float(scib_metrics.clisi_knn(nn, labels_lisi, scale=True))
        print(f"    cLISI: {metrics['cLISI']:.4f}")
    except Exception as e:
        print(f"    WARNING: cLISI failed: {e}")
        metrics['cLISI'] = np.nan

    try:
        emb_asw, meta_asw, _ = _subsample(embedding, metadata, MAX_CELLS_ASW)
        labels_asw = meta_asw['coarse_label'].values
        metrics['bio_ASW'] = float(scib_metrics.silhouette_label(
            emb_asw, labels_asw, rescale=True
        ))
        print(f"    bio_ASW: {metrics['bio_ASW']:.4f}")
    except Exception as e:
        print(f"    WARNING: bio_ASW failed: {e}")
        metrics['bio_ASW'] = np.nan

    return metrics


def compute_condition_metrics(embedding, metadata):
    metrics = {}
    if 'condition_harmonized' not in metadata.columns:
        return {'condition_ASW': np.nan, 'condition_LISI': np.nan}

    cond = metadata['condition_harmonized'].values
    valid = pd.notna(cond)
    if valid.sum() < 100 or len(set(cond[valid])) < 2:
        return {'condition_ASW': np.nan, 'condition_LISI': np.nan}

    try:
        emb_asw, meta_asw, _ = _subsample(embedding, metadata, MAX_CELLS_ASW)
        cond_asw = meta_asw['condition_harmonized'].values
        valid_asw = pd.notna(cond_asw)
        metrics['condition_ASW'] = float(silhouette_score(
            emb_asw[valid_asw], cond_asw[valid_asw],
            sample_size=min(10_000, int(valid_asw.sum())),
            random_state=RNG_SEED
        ))
        print(f"    condition_ASW: {metrics['condition_ASW']:.4f}")
    except Exception as e:
        print(f"    WARNING: condition_ASW failed: {e}")
        metrics['condition_ASW'] = np.nan

    try:
        emb_lisi, meta_lisi, _ = _subsample(embedding, metadata, MAX_CELLS_LISI)
        cond_lisi = meta_lisi['condition_harmonized'].values
        k_lisi = _knn_neighbors(emb_lisi.shape[0])
        nn = pynndescent(emb_lisi, n_neighbors=k_lisi,
                         random_state=RNG_SEED, n_jobs=KNN_JOBS)
        lisi_per_cell = scib_metrics.lisi_knn(nn, cond_lisi)
        metrics['condition_LISI'] = float(np.mean(lisi_per_cell))
        print(f"    condition_LISI: {metrics['condition_LISI']:.4f}")
    except Exception as e:
        print(f"    WARNING: condition_LISI failed: {e}")
        metrics['condition_LISI'] = np.nan

    return metrics


def compute_cluster_agreement(embedding, metadata):
    metrics = {}
    if 'coarse_label' not in metadata.columns:
        return {'ARI': np.nan, 'NMI': np.nan}

    labels = metadata['coarse_label'].values
    if len(set(labels)) < 2:
        return {'ARI': np.nan, 'NMI': np.nan}

    try:
        emb_sub, meta_sub, _ = _subsample(embedding, metadata, MAX_CELLS_LISI)
        labels_sub = meta_sub['coarse_label'].values
        k_lisi = _knn_neighbors(emb_sub.shape[0])
        nn = pynndescent(emb_sub, n_neighbors=k_lisi,
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
    metrics = {}
    if counts is None or gene_names is None:
        return {f"var_ratio_{g}": np.nan for g in MARKER_GENES}

    if counts.shape[0] < 50:
        print("    WARNING: <50 cells, skipping marker variance")
        return {f"var_ratio_{g}": np.nan for g in MARKER_GENES}

    adata = sc.AnnData(X=counts)
    adata.obsm['X_emb'] = embedding

    sc.pp.neighbors(adata, use_rep='X_emb', n_neighbors=15)
    sc.tl.leiden(adata, resolution=0.5, flavor='igraph', n_iterations=2)
    clusters = adata.obs['leiden'].values

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    for gene in MARKER_GENES:
        key = f"var_ratio_{gene}"
        if gene not in gene_names:
            gene_alt = gene.replace("_", "-")
            if gene_alt in gene_names:
                gene_idx = gene_names.index(gene_alt)
            else:
                print(f"    WARNING: {gene} not found")
                metrics[key] = np.nan
                continue
        else:
            gene_idx = gene_names.index(gene)

        expr = np.asarray(adata.X[:, gene_idx].todense()).flatten()
        total_var = np.var(expr)
        if total_var < 1e-10:
            metrics[key] = np.nan
            continue

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


def compute_all_metrics(bridge_dir, run_name, scope_name):
    print(f"\n{'='*60}")
    print(f"  {run_name} / {scope_name}")
    print(f"{'='*60}")

    bridge_dir = Path(bridge_dir)
    if not bridge_dir.exists():
        print(f"  SKIP: {bridge_dir} does not exist")
        return None

    if _detect_embedding_file(bridge_dir) is None:
        print(f"  SKIP: no embedding file in {bridge_dir}")
        return None

    embedding, metadata, counts, gene_names = load_bridge_data(bridge_dir)

    result = {
        'run': run_name,
        'scope': scope_name,
        'n_cells': embedding.shape[0],
        'n_dims': embedding.shape[1],
    }
    print("\n  --- Batch removal ---")
    result.update(compute_batch_metrics(embedding, metadata))
    print("\n  --- Bio preservation ---")
    result.update(compute_bio_metrics(embedding, metadata))
    print("\n  --- Condition ---")
    result.update(compute_condition_metrics(embedding, metadata))
    print("\n  --- Cluster agreement ---")
    result.update(compute_cluster_agreement(embedding, metadata))
    print("\n  --- Marker variance ---")
    result.update(compute_marker_variance(counts, gene_names, embedding, metadata))
    return result


# ═════════════════════════════════════════════════════════════════════════
# REPORTING
# ═════════════════════════════════════════════════════════════════════════

HTML_TEMPLATE = """<!DOCTYPE html>
<html>
<head>
<title>{compartment} Integration — Tiered v4 vs Baseline v5 CCA</title>
<style>
  body {{ font-family: Arial, sans-serif; max-width: 1600px; margin: 0 auto; padding: 20px; }}
  h1 {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
  table {{ border-collapse: collapse; width: 100%; margin: 15px 0; }}
  th, td {{ border: 1px solid #ddd; padding: 8px 12px; text-align: right; }}
  th {{ background-color: #3498db; color: white; text-align: center; }}
  td:first-child, td:nth-child(2) {{ text-align: left; }}
  tr:nth-child(even) {{ background-color: #f9f9f9; }}
  tr.baseline {{ background-color: #fff3cd; font-weight: bold; }}
  .summary {{ background: #eaf2f8; padding: 15px; border-radius: 5px; margin: 20px 0; }}
</style>
</head>
<body>
<h1>{compartment} Integration Quality — Tiered v4</h1>
<p>Generated: {date}</p>
<div class="summary">
<p>Comparing the production <strong>baseline_flat_v5</strong> CCA workflow
to the new <strong>tiered_v4</strong> (Seurat v4 SCT + CCA, mes / non-mes
split) for the {compartment} compartment.</p>
</div>
<h2>Metrics</h2>
{table_html}
</body>
</html>"""


def build_comparison_table(all_results, compartment, results_dir):
    df = pd.DataFrame(all_results)
    id_cols = ['run', 'scope', 'n_cells', 'n_dims']
    metric_cols = ['iLISI', 'batch_ASW', 'cLISI', 'bio_ASW',
                   'condition_ASW', 'condition_LISI', 'NMI', 'ARI']
    var_cols = [c for c in df.columns if c.startswith('var_ratio_')]
    ordered = id_cols + metric_cols + sorted(var_cols)
    ordered = [c for c in ordered if c in df.columns]
    df = df[ordered]

    tsv_path = results_dir / "comparison_table.tsv"
    df.to_csv(tsv_path, sep='\t', index=False, float_format='%.4f')
    print(f"\n  Saved: {tsv_path}")

    display_df = df.copy()
    for col in display_df.columns:
        if col in id_cols:
            continue
        display_df[col] = display_df[col].apply(
            lambda x: f"{x:.4f}" if pd.notna(x) and isinstance(x, (int, float, np.floating)) else str(x)
        )

    table_html = display_df.to_html(index=False, na_rep='N/A', escape=False)
    table_html = table_html.replace(
        '<td>baseline_flat_v5</td>',
        '<td class="baseline">baseline_flat_v5</td>'
    )

    html = HTML_TEMPLATE.format(
        compartment=compartment,
        date=datetime.now().strftime("%Y-%m-%d %H:%M"),
        table_html=table_html,
    )
    html_path = results_dir / "comparison_table.html"
    html_path.write_text(html)
    print(f"  Saved: {html_path}")

    print(f"\n{'='*80}")
    print(f"COMPARISON TABLE — {compartment}")
    print(f"{'='*80}")
    print(df.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    return df


def run_compartment(compartment):
    print("=" * 60)
    print(f"Tiered v4 metrics — {compartment}")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    runs = build_run_scopes(compartment)
    cfg = COMPARTMENT_RUNS[compartment]
    results_dir = RESULTS_BASE / cfg["exp_dir_name"]
    results_dir.mkdir(parents=True, exist_ok=True)

    all_results = []
    for run_name, scopes in runs.items():
        for scope_info in scopes:
            result = compute_all_metrics(
                scope_info['bridge_dir'], run_name, scope_info['scope']
            )
            if result is not None:
                all_results.append(result)

    if not all_results:
        print(f"\nNo results for {compartment}.")
        return None

    build_comparison_table(all_results, compartment, results_dir)
    print(f"\n{compartment} complete: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    return all_results


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--compartment', type=str, default=None,
                        choices=list(COMPARTMENT_RUNS.keys()),
                        help="Compartment to process (default: all)")
    args = parser.parse_args()

    targets = [args.compartment] if args.compartment else list(COMPARTMENT_RUNS.keys())
    for compartment in targets:
        run_compartment(compartment)


if __name__ == "__main__":
    main()
