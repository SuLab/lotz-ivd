#!/usr/bin/env python3
"""NP Integration Experiment — Continuum Preservation Controls.

Adds two analyses that address the cluster-count confound in the original
marker-variance metric from `05h_np_experiment_metrics.py`:

  Option 3 — Cluster-free (KNN neighborhood variance). For each cell,
      compute variance of marker expression across its k=50 nearest
      neighbors in the integrated embedding. Average over cells, divide
      by total variance. Independent of any clustering choice; fixed k
      gives the same "resolution" across arms.

  Option 2 — Resolution sweep (Leiden at 5 resolutions). For each arm,
      cluster at res in {0.1, 0.25, 0.5, 1.0, 2.0}; for each (res, marker)
      pair compute n_clusters and within-cluster var_ratio. Plotting
      var_ratio vs n_clusters lets us read off the cluster-count-matched
      comparison (Option 1) and test whether the v5/v4 ordering is stable
      across resolutions.

Outputs:
  continuum_knn_var_ratio.tsv   — one row per arm-scope, k=50 KNN variance
  continuum_sweep.tsv           — long-form: arm × scope × res × markers

Usage:
    python3 scripts/05j_continuum_control_metrics.py
    python3 scripts/05j_continuum_control_metrics.py --run tiered_v5
    python3 scripts/05j_continuum_control_metrics.py --skip-sweep  # only KNN
"""

import argparse
import sys
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread
from scipy.sparse import csc_matrix

from scib_metrics.nearest_neighbors import pynndescent

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

# ── Paths ─────────────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent.parent
EXP_DIR = BASE / "data" / "integrated" / "np_experiment"
BASELINE_DIR = BASE / "data" / "integrated" / "cca" / "bridge_export" / "NP"
RESULTS_DIR = BASE / "results" / "integration" / "np_experiment"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ── Parameters ────────────────────────────────────────────────────────────
KNN_K = 50                          # neighborhood size for cluster-free metric
KNN_JOBS = 16
SWEEP_RESOLUTIONS = [0.1, 0.25, 0.5, 1.0, 2.0]
MARKER_GENES = ["COL2A1", "ACAN", "SOX9", "COL1A1"]
RNG_SEED = 42
MIN_CELLS_PER_STUDY = 500           # skip studies smaller than this in within-study test

# Each run/scope: embedding + metadata + counts
RUN_SCOPES = {
    "baseline_flat_v5": [
        {"scope": "all", "bridge_dir": BASELINE_DIR},
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
# DATA LOADING  (mirrors 05h_np_experiment_metrics.py)
# ═════════════════════════════════════════════════════════════════════════

def _detect_embedding_file(bridge_dir):
    for name in ["embedding_integrated.cca.csv.gz", "embedding_pca.csv.gz"]:
        path = bridge_dir / name
        if path.exists():
            return path
    emb_files = sorted(bridge_dir.glob("embedding_*.csv.gz"))
    return emb_files[0] if emb_files else None


def load_bridge_data(bridge_dir):
    bridge_dir = Path(bridge_dir)

    emb_path = _detect_embedding_file(bridge_dir)
    if emb_path is None:
        raise FileNotFoundError(f"No embedding file in {bridge_dir}")
    print(f"    Embedding: {emb_path.name}")
    embedding = pd.read_csv(emb_path, index_col=0).values.astype(np.float32)

    meta_path = bridge_dir / "metadata.csv.gz"
    if not meta_path.exists():
        meta_path = bridge_dir / "metadata.csv"
    metadata = pd.read_csv(meta_path)
    if len(metadata) != embedding.shape[0]:
        n = min(len(metadata), embedding.shape[0])
        metadata = metadata.iloc[:n]
        embedding = embedding[:n]

    mtx_path = bridge_dir / "counts.mtx.gz"
    if not mtx_path.exists():
        mtx_path = bridge_dir / "counts.mtx"
    genes_path = bridge_dir / "genes.csv"

    if not (mtx_path.exists() and genes_path.exists()):
        raise FileNotFoundError(f"No counts in {bridge_dir}")

    print(f"    Loading counts...")
    raw = mmread(str(mtx_path))
    if raw.shape[0] != embedding.shape[0]:
        raw = raw.T
    counts = csc_matrix(raw)
    gene_names = pd.read_csv(genes_path)["gene"].tolist()
    print(f"    {counts.shape[0]} cells x {counts.shape[1]} genes, {embedding.shape[1]} dims")
    return embedding, metadata, counts, gene_names


def gene_idx(gene, gene_names):
    if gene in gene_names:
        return gene_names.index(gene)
    alt = gene.replace("_", "-")
    if alt in gene_names:
        return gene_names.index(alt)
    return None


def logn_expr(counts, gene_names, markers):
    """Return {gene: log1p-normalized expression vector}."""
    adata = sc.AnnData(X=counts)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    out = {}
    for g in markers:
        idx = gene_idx(g, gene_names)
        if idx is None:
            print(f"    WARNING: {g} not found")
            out[g] = None
        else:
            out[g] = np.asarray(adata.X[:, idx].todense()).flatten()
    del adata
    return out


# ═════════════════════════════════════════════════════════════════════════
# OPTION 3 — CLUSTER-FREE KNN METRICS
# ═════════════════════════════════════════════════════════════════════════

def compute_knn_spatial_metrics(embedding, expr_by_gene, k=KNN_K):
    """Two cluster-free metrics on a single k-NN graph:

      knn_var_ratio  — mean(var(expr over each cell's k-NN)) / total var.
                       Local: asks "are my neighbors similar TO EACH OTHER?"
                       Low = smooth gradient; high (→1) = no local structure.

      morans_i       — Moran's I spatial autocorrelation on the k-NN graph
                       with binary weights. Global: asks "am I similar to
                       THE MEAN OF my neighbors, relative to total variance?"
                       With binary KNN weights this reduces to
                           I = sum_i z_i · mean_neighbors(z_i) / sum_i z_i²
                       where z = x - mean(x). Range roughly [-1, 1];
                       higher = stronger spatial autocorrelation
                       = continuum preserved at graph scale.

    The two answer different questions and can diverge: neighborhood
    variance misses global pocket arrangement that Moran's I picks up.
    """
    print(f"    Building KNN (k={k}) on full embedding...")
    # pynndescent returns NeighborsResults with .indices (n_cells, k+1 incl self)
    nn = pynndescent(embedding, n_neighbors=k + 1,
                     random_state=RNG_SEED, n_jobs=KNN_JOBS)
    # Drop self (first column is always the cell itself)
    idx = nn.indices[:, 1:]                 # (n_cells, k)

    out = {}
    for g, expr in expr_by_gene.items():
        if expr is None:
            out[f"knn_var_ratio_{g}"] = np.nan
            out[f"morans_i_{g}"] = np.nan
            continue
        total_var = float(np.var(expr))
        if total_var < 1e-10:
            out[f"knn_var_ratio_{g}"] = np.nan
            out[f"morans_i_{g}"] = np.nan
            continue

        # --- Local: KNN neighborhood variance ---
        neighbor_expr = expr[idx]                     # (n_cells, k)
        neighbor_var = np.var(neighbor_expr, axis=1)  # (n_cells,)
        knn_ratio = float(np.mean(neighbor_var) / total_var)
        out[f"knn_var_ratio_{g}"] = knn_ratio

        # --- Global: Moran's I with binary KNN weights ---
        z = expr - expr.mean()
        neighbor_sum_z = z[idx].sum(axis=1)           # (n_cells,)
        denom = float(np.sum(z ** 2))
        numerator = float(np.sum(z * neighbor_sum_z))
        morans_i = numerator / (k * denom)
        out[f"morans_i_{g}"] = morans_i

        print(f"    {g}: knn_var_ratio={knn_ratio:.4f}  morans_i={morans_i:.4f}")
    return out


# ═════════════════════════════════════════════════════════════════════════
# OPTION 3c — WITHIN-STUDY MORAN'S I (batch-confound test)
# ═════════════════════════════════════════════════════════════════════════

def _morans_i_binary_knn(expr, knn_indices, k):
    """Moran's I with binary KNN weights on a pre-built neighbor index."""
    z = expr - expr.mean()
    denom = float(np.sum(z ** 2))
    if denom < 1e-12:
        return np.nan
    neighbor_sum_z = z[knn_indices].sum(axis=1)
    return float(np.sum(z * neighbor_sum_z) / (k * denom))


def compute_within_study_morans_i(embedding, metadata, expr_by_gene,
                                   k=KNN_K, min_cells=MIN_CELLS_PER_STUDY):
    """For each study with at least `min_cells` cells, build a KNN graph
    restricted to that study's cells and compute Moran's I per marker.

    Tests whether the pooled Moran's I is driven by *within-study* spatial
    structure (biology) or *between-study* structure (batch). If a method
    fails to mix studies, the pooled Moran's I captures the inter-study
    separation of cells that differ systematically in marker expression,
    inflating the value even when within-study biology is flat.

    Returns list of dicts, one per (study, gene).
    """
    rows = []
    studies = metadata["study"].value_counts()
    studies = studies[studies >= min_cells].index.tolist()
    print(f"    Studies with ≥{min_cells} cells: {len(studies)}")

    for study in studies:
        mask = (metadata["study"] == study).values
        n_study = int(mask.sum())
        emb_s = embedding[mask]
        # Safe k for this study
        k_study = min(k, n_study - 2)
        nn = pynndescent(emb_s, n_neighbors=k_study + 1,
                         random_state=RNG_SEED, n_jobs=KNN_JOBS)
        idx_s = nn.indices[:, 1:]
        for g, expr in expr_by_gene.items():
            if expr is None:
                rows.append({"study": study, "n_cells_study": n_study,
                             "k_used": k_study, "gene": g, "morans_i": np.nan})
                continue
            expr_s = expr[mask]
            mi = _morans_i_binary_knn(expr_s, idx_s, k_study)
            rows.append({"study": study, "n_cells_study": n_study,
                         "k_used": k_study, "gene": g, "morans_i": mi})
            print(f"    {study} (n={n_study}): {g} Moran's I = {mi:.4f}")
    return rows


# ═════════════════════════════════════════════════════════════════════════
# OPTION 2 — LEIDEN RESOLUTION SWEEP
# ═════════════════════════════════════════════════════════════════════════

def compute_sweep(embedding, expr_by_gene, resolutions=SWEEP_RESOLUTIONS):
    """Run Leiden at several resolutions; for each resolution × marker,
    record n_clusters and cluster-weighted within-cluster var ratio.
    """
    print(f"    Building scanpy neighbors for Leiden...")
    adata = sc.AnnData(X=np.zeros((embedding.shape[0], 1), dtype=np.float32))
    adata.obsm["X_emb"] = embedding
    sc.pp.neighbors(adata, use_rep="X_emb", n_neighbors=15)

    rows = []
    for res in resolutions:
        t0 = datetime.now()
        sc.tl.leiden(adata, resolution=res, flavor="igraph",
                     n_iterations=2, key_added=f"_lei_{res}")
        clusters = adata.obs[f"_lei_{res}"].astype(str).values
        k_clusters = len(set(clusters))
        elapsed = (datetime.now() - t0).total_seconds()
        print(f"    res={res}: {k_clusters} clusters  ({elapsed:.0f}s)")

        for g, expr in expr_by_gene.items():
            if expr is None:
                rows.append({"resolution": res, "n_clusters": k_clusters,
                             "gene": g, "var_ratio": np.nan})
                continue
            total_var = float(np.var(expr))
            if total_var < 1e-10:
                rows.append({"resolution": res, "n_clusters": k_clusters,
                             "gene": g, "var_ratio": np.nan})
                continue
            # Cluster-size-weighted mean of within-cluster variance
            cluster_vars, cluster_sizes = [], []
            for cl in np.unique(clusters):
                mask = clusters == cl
                n_cl = mask.sum()
                if n_cl < 5:
                    continue
                cluster_vars.append(np.var(expr[mask]))
                cluster_sizes.append(n_cl)
            if not cluster_vars:
                rows.append({"resolution": res, "n_clusters": k_clusters,
                             "gene": g, "var_ratio": np.nan})
                continue
            cv = np.array(cluster_vars)
            cs = np.array(cluster_sizes)
            within_var = np.average(cv, weights=cs)
            rows.append({"resolution": res, "n_clusters": k_clusters,
                         "gene": g, "var_ratio": float(within_var / total_var)})

    del adata
    return rows


# ═════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=str, default=None,
                        choices=list(RUN_SCOPES))
    parser.add_argument("--skip-sweep", action="store_true",
                        help="Skip Leiden resolution sweep")
    parser.add_argument("--skip-knn", action="store_true",
                        help="Skip pooled KNN metrics")
    parser.add_argument("--skip-within-study", action="store_true",
                        help="Skip within-study Moran's I (batch-confound test)")
    parser.add_argument("--only-within-study", action="store_true",
                        help="Only compute within-study Moran's I")
    args = parser.parse_args()
    if args.only_within_study:
        args.skip_knn = True
        args.skip_sweep = True

    print("=" * 60)
    print("NP Continuum Preservation Controls")
    print(f"Started: {datetime.now():%Y-%m-%d %H:%M:%S}")
    print(f"  k (cluster-free KNN): {KNN_K}")
    print(f"  Leiden sweep resolutions: {SWEEP_RESOLUTIONS}")
    print("=" * 60)

    runs = {args.run: RUN_SCOPES[args.run]} if args.run else RUN_SCOPES

    knn_rows = []
    sweep_rows = []
    within_study_rows = []

    for run_name, scopes in runs.items():
        for scope_info in scopes:
            scope = scope_info["scope"]
            bridge_dir = scope_info["bridge_dir"]
            print(f"\n{'='*60}\n  {run_name} / {scope}\n{'='*60}")
            if not bridge_dir.exists():
                print(f"  SKIP: missing {bridge_dir}")
                continue

            embedding, metadata, counts, gene_names = load_bridge_data(bridge_dir)

            print("  Computing log-normalized marker expression...")
            expr_by_gene = logn_expr(counts, gene_names, MARKER_GENES)
            # counts no longer needed
            del counts

            if not args.skip_knn:
                print("  --- Option 3: cluster-free KNN metrics ---")
                knn_res = compute_knn_spatial_metrics(embedding, expr_by_gene)
                knn_rows.append({
                    "run": run_name, "scope": scope,
                    "n_cells": embedding.shape[0], "k": KNN_K,
                    **knn_res,
                })

            if not args.skip_sweep:
                print("  --- Option 2: Leiden resolution sweep ---")
                for r in compute_sweep(embedding, expr_by_gene):
                    sweep_rows.append({
                        "run": run_name, "scope": scope,
                        "n_cells": embedding.shape[0], **r,
                    })

            if not args.skip_within_study:
                print("  --- Option 3c: within-study Moran's I ---")
                for r in compute_within_study_morans_i(
                        embedding, metadata, expr_by_gene):
                    within_study_rows.append({
                        "run": run_name, "scope": scope, **r,
                    })

            del embedding, expr_by_gene

    # Save outputs
    if knn_rows:
        # Filename kept as continuum_knn_var_ratio.tsv for backwards
        # compatibility; schema now also carries morans_i_<GENE> columns.
        knn_path = RESULTS_DIR / "continuum_knn_var_ratio.tsv"
        pd.DataFrame(knn_rows).to_csv(knn_path, sep="\t", index=False,
                                      float_format="%.4f")
        print(f"\nSaved: {knn_path}")
    if sweep_rows:
        sweep_path = RESULTS_DIR / "continuum_sweep.tsv"
        pd.DataFrame(sweep_rows).to_csv(sweep_path, sep="\t", index=False,
                                        float_format="%.4f")
        print(f"Saved: {sweep_path}")
    if within_study_rows:
        ws_path = RESULTS_DIR / "continuum_within_study_morans_i.tsv"
        pd.DataFrame(within_study_rows).to_csv(ws_path, sep="\t", index=False,
                                               float_format="%.4f")
        print(f"Saved: {ws_path}")

    print(f"\nCompleted: {datetime.now():%Y-%m-%d %H:%M:%S}")


if __name__ == "__main__":
    main()
