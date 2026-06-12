#!/usr/bin/env python3
"""Module 05n — Harmony integration sweep across the four tiered_v4 compartments.

Standalone re-running of the Harmony integration method on the same data as the
Seurat v4 tiered_v4 production pipeline, so the comparison table in the manuscript
Methods sub-section §Integration-method comparison includes a contemporaneous
Harmony row alongside Seurat v4/v5 (flat + tiered), scANVI, and STACAS.

Pipeline per compartment:
  1. Load data/integrated/tiered_v4/<compartment>.h5ad (raw counts in .X)
  2. CP10K + log1p; HVG selection (n_top_genes=3000, flavor=seurat_v3)
  3. Scale + PCA (50 dims) on HVGs
  4. harmonypy correction via scanpy.external.pp.harmony_integrate, group=study
  5. Compute the same metric battery as np_experiment/comparison_table.tsv:
     iLISI, batch_ASW, cLISI, bio_ASW, condition_ASW, condition_LISI, NMI, ARI,
     var_ratio_{ACAN,COL1A1,COL2A1,SOX9}.
  6. Save corrected embedding + metrics to results/integration/harmony/<compartment>/

Sequential by compartment to keep peak RAM bounded. Per-compartment results
persist as soon as they finish, so a mid-run crash leaves a partial sweep
recoverable.

Usage:
  python3 scripts/05n_harmony_integration.py
  python3 scripts/05n_harmony_integration.py --compartment NP
  python3 scripts/05n_harmony_integration.py --order CEP AF NP all_cells
"""
from __future__ import annotations
import argparse, json, os, sys, time, warnings
from pathlib import Path
import numpy as np, pandas as pd
import scanpy as sc
import scanpy.external as sce
from sklearn.metrics import silhouette_score

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

import scib_metrics
from scib_metrics.nearest_neighbors import pynndescent

BASE = Path(__file__).resolve().parent.parent
H5AD = {c: BASE / "data" / "integrated" / "tiered_v4" / f"{c}.h5ad"
        for c in ["NP", "AF", "CEP", "all_cells"]}
OUTROOT = BASE / "results" / "integration" / "harmony"
OUTROOT.mkdir(parents=True, exist_ok=True)

MAX_CELLS_LISI = 50_000
MAX_CELLS_ASW = 30_000
KNN_NEIGHBORS = 90
KNN_JOBS = 16
MARKER_GENES = ["COL2A1", "ACAN", "SOX9", "COL1A1"]
RNG_SEED = 42
HVG_N = 3000
PCA_DIMS = 50


def _ts():
    return time.strftime("%H:%M:%S")


def _subsample(emb, meta, n_max):
    n = emb.shape[0]
    if n <= n_max:
        return emb, meta, np.arange(n)
    rng = np.random.RandomState(RNG_SEED)
    idx = np.sort(rng.choice(n, n_max, replace=False))
    return emb[idx], meta.iloc[idx].reset_index(drop=True), idx


def compute_metrics(embedding, obs, counts_norm_log, gene_names, do_var=True):
    m = {}
    # ── Batch metrics
    batch = obs["study"].astype(str).values
    if len(set(batch)) >= 2:
        try:
            elisi, mlisi, _ = _subsample(embedding, obs, MAX_CELLS_LISI)
            nn = pynndescent(elisi, n_neighbors=KNN_NEIGHBORS, random_state=RNG_SEED, n_jobs=KNN_JOBS)
            m["iLISI"] = float(scib_metrics.ilisi_knn(nn, mlisi["study"].astype(str).values, scale=True))
            print(f"[{_ts()}]   iLISI={m['iLISI']:.4f}")
        except Exception as e:
            print(f"[{_ts()}]   iLISI FAILED: {e}"); m["iLISI"] = np.nan
        try:
            easw, masw, _ = _subsample(embedding, obs, MAX_CELLS_ASW)
            labels = masw["coarse_label"].astype(str).values if "coarse_label" in masw.columns else None
            ba = masw["study"].astype(str).values
            if labels is not None and len(set(labels)) > 1:
                m["batch_ASW"] = float(scib_metrics.silhouette_batch(easw, labels, ba, rescale=True))
            else:
                m["batch_ASW"] = float(silhouette_score(easw, ba,
                                                        sample_size=min(5000, len(easw)),
                                                        random_state=RNG_SEED))
            print(f"[{_ts()}]   batch_ASW={m['batch_ASW']:.4f}")
        except Exception as e:
            print(f"[{_ts()}]   batch_ASW FAILED: {e}"); m["batch_ASW"] = np.nan

    # ── Bio metrics
    if "coarse_label" in obs.columns:
        labels = obs["coarse_label"].astype(str).values
        if len(set(labels)) >= 2:
            try:
                elisi, mlisi, _ = _subsample(embedding, obs, MAX_CELLS_LISI)
                nn = pynndescent(elisi, n_neighbors=KNN_NEIGHBORS, random_state=RNG_SEED, n_jobs=KNN_JOBS)
                m["cLISI"] = float(scib_metrics.clisi_knn(nn, mlisi["coarse_label"].astype(str).values, scale=True))
                print(f"[{_ts()}]   cLISI={m['cLISI']:.4f}")
            except Exception as e:
                print(f"[{_ts()}]   cLISI FAILED: {e}"); m["cLISI"] = np.nan
            try:
                easw, masw, _ = _subsample(embedding, obs, MAX_CELLS_ASW)
                m["bio_ASW"] = float(scib_metrics.silhouette_label(
                    easw, masw["coarse_label"].astype(str).values, rescale=True))
                print(f"[{_ts()}]   bio_ASW={m['bio_ASW']:.4f}")
            except Exception as e:
                print(f"[{_ts()}]   bio_ASW FAILED: {e}"); m["bio_ASW"] = np.nan

    # ── Condition metrics
    if "condition_harmonized" in obs.columns:
        cond = obs["condition_harmonized"].astype(str).values
        if len(set(cond)) >= 2:
            try:
                easw, masw, _ = _subsample(embedding, obs, MAX_CELLS_ASW)
                ca = masw["condition_harmonized"].astype(str).values
                m["condition_ASW"] = float(silhouette_score(easw, ca,
                                                            sample_size=min(10000, len(easw)),
                                                            random_state=RNG_SEED))
                print(f"[{_ts()}]   condition_ASW={m['condition_ASW']:.4f}")
            except Exception as e:
                print(f"[{_ts()}]   condition_ASW FAILED: {e}"); m["condition_ASW"] = np.nan
            try:
                elisi, mlisi, _ = _subsample(embedding, obs, MAX_CELLS_LISI)
                nn = pynndescent(elisi, n_neighbors=KNN_NEIGHBORS, random_state=RNG_SEED, n_jobs=KNN_JOBS)
                lisi_pc = scib_metrics.lisi_knn(nn, mlisi["condition_harmonized"].astype(str).values)
                m["condition_LISI"] = float(np.mean(lisi_pc))
                print(f"[{_ts()}]   condition_LISI={m['condition_LISI']:.4f}")
            except Exception as e:
                print(f"[{_ts()}]   condition_LISI FAILED: {e}"); m["condition_LISI"] = np.nan

    # ── NMI / ARI (Leiden vs coarse_label)
    if "coarse_label" in obs.columns and len(set(obs["coarse_label"].astype(str))) >= 2:
        try:
            esub, msub, _ = _subsample(embedding, obs, MAX_CELLS_LISI)
            nn = pynndescent(esub, n_neighbors=KNN_NEIGHBORS, random_state=RNG_SEED, n_jobs=KNN_JOBS)
            r = scib_metrics.nmi_ari_cluster_labels_leiden(
                nn, msub["coarse_label"].astype(str).values, optimize_resolution=True, seed=RNG_SEED)
            m["NMI"] = float(r["nmi"]); m["ARI"] = float(r["ari"])
            print(f"[{_ts()}]   NMI={m['NMI']:.4f}  ARI={m['ARI']:.4f}")
        except Exception as e:
            print(f"[{_ts()}]   NMI/ARI FAILED: {e}"); m["NMI"] = np.nan; m["ARI"] = np.nan

    # ── Marker variance ratios (Leiden clusters on the embedding + log-norm gene matrix)
    if do_var and counts_norm_log is not None and gene_names is not None:
        try:
            tmp = sc.AnnData(X=counts_norm_log)
            tmp.obsm["X_emb"] = embedding
            sc.pp.neighbors(tmp, use_rep="X_emb", n_neighbors=15)
            sc.tl.leiden(tmp, resolution=0.5, flavor="igraph", n_iterations=2)
            clusters = tmp.obs["leiden"].astype(str).values
            for gene in MARKER_GENES:
                key = f"var_ratio_{gene}"
                if gene not in gene_names:
                    m[key] = np.nan; continue
                gidx = gene_names.index(gene)
                expr = np.asarray(tmp.X[:, gidx].todense()).flatten() if hasattr(tmp.X, "todense") else np.asarray(tmp.X[:, gidx]).flatten()
                total_var = float(np.var(expr))
                if total_var < 1e-10:
                    m[key] = np.nan; continue
                cv, cs = [], []
                for c in np.unique(clusters):
                    mask = clusters == c
                    if mask.sum() < 5: continue
                    cv.append(float(np.var(expr[mask])))
                    cs.append(int(mask.sum()))
                if cv:
                    m[key] = float(np.average(cv, weights=cs) / total_var)
                else:
                    m[key] = np.nan
            print(f"[{_ts()}]   var_ratios computed for {MARKER_GENES}")
        except Exception as e:
            print(f"[{_ts()}]   var_ratio FAILED: {e}")
            for g in MARKER_GENES: m[f"var_ratio_{g}"] = np.nan

    return m


def run_compartment(comp: str):
    out = OUTROOT / comp
    out.mkdir(parents=True, exist_ok=True)
    flag = out / "DONE"
    if flag.exists():
        print(f"[{_ts()}] [{comp}] DONE flag exists — skipping (delete {flag} to re-run)")
        return
    t0 = time.time()
    print(f"[{_ts()}] [{comp}] === START ===")
    print(f"[{_ts()}] [{comp}] loading {H5AD[comp]}")
    ad = sc.read_h5ad(H5AD[comp])
    print(f"[{_ts()}] [{comp}] shape={ad.shape}, X dtype={ad.X.dtype}")

    # 1. Normalize + log1p (preserves raw counts in a layer for var_ratios)
    print(f"[{_ts()}] [{comp}] normalize_total + log1p")
    ad.layers["counts"] = ad.X.copy()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)

    # 2. HVG (seurat_v3 flavor wants raw counts as a layer, but standard works on log1p too)
    print(f"[{_ts()}] [{comp}] HVG n_top_genes={HVG_N}")
    sc.pp.highly_variable_genes(ad, n_top_genes=HVG_N, batch_key="study", subset=False, flavor="seurat")

    # 3. Scale + PCA on HVGs
    print(f"[{_ts()}] [{comp}] PCA on {ad.var['highly_variable'].sum()} HVGs, dims={PCA_DIMS}")
    hv = ad.var["highly_variable"].values
    ad_hv = ad[:, hv].copy()
    sc.pp.scale(ad_hv, max_value=10)
    sc.tl.pca(ad_hv, n_comps=PCA_DIMS, random_state=RNG_SEED)
    ad.obsm["X_pca"] = ad_hv.obsm["X_pca"]

    # 4. Harmony grouped by study
    # NB: call harmonypy directly rather than sce.pp.harmony_integrate. The
    # PyTorch-backend harmonypy (0.2.0) returns Z_corr already oriented as
    # (n_cells, n_dims), whereas legacy harmonypy returned (n_dims, n_cells).
    # scanpy's wrapper unconditionally transposes, so it (and a naive port)
    # produces the wrong shape and anndata rejects the .obsm assignment. We
    # orient by matching n_obs, which is correct for either harmonypy version.
    print(f"[{_ts()}] [{comp}] Harmony — group=study, max_iter=20, theta=2.0")
    import harmonypy
    ho = harmonypy.run_harmony(np.asarray(ad.obsm["X_pca"]), ad.obs, ["study"],
                               max_iter_harmony=20, theta=2.0, random_state=RNG_SEED)
    Z = ho.Z_corr
    if hasattr(Z, "cpu"):          # torch tensor -> numpy
        Z = Z.cpu().numpy()
    Z = np.asarray(Z)
    n_obs = ad.n_obs
    if Z.shape[0] != n_obs and Z.shape[1] == n_obs:
        Z = Z.T                    # legacy (n_dims, n_cells) -> (n_cells, n_dims)
    if Z.shape[0] != n_obs:
        raise ValueError(f"Harmony embedding {Z.shape} does not match n_obs={n_obs}")
    emb = np.ascontiguousarray(Z)  # (n_cells, n_dims)
    ad.obsm["X_pca_harmony"] = emb
    print(f"[{_ts()}] [{comp}] Harmony done. X_pca_harmony shape={emb.shape}")

    # 5. Compute metrics
    print(f"[{_ts()}] [{comp}] computing metrics on {emb.shape[0]:,} × {emb.shape[1]} embedding")
    # For marker variance, use log1p-normalized expression matrix subset to marker genes
    counts_norm_log = ad.X  # already log1p of CP10K
    metrics = compute_metrics(emb, ad.obs, counts_norm_log, list(ad.var_names))

    # 6. Save
    np.save(out / "embedding_harmony.npy", emb)
    pd.DataFrame({"barcode": ad.obs_names.astype(str)}).to_csv(out / "cell_index.csv.gz", index=False)
    metrics_row = {"method": "harmony", "compartment": comp, "scope": "all",
                   "n_cells": int(emb.shape[0]), "n_dims": int(emb.shape[1]),
                   **metrics}
    pd.DataFrame([metrics_row]).to_csv(out / "metrics.tsv", sep="\t", index=False)
    with open(out / "params.json", "w") as f:
        json.dump({
            "hvg_n": HVG_N, "pca_dims": PCA_DIMS,
            "harmony_max_iter": 20, "harmony_theta": 2.0,
            "max_cells_lisi": MAX_CELLS_LISI, "max_cells_asw": MAX_CELLS_ASW,
            "knn_neighbors": KNN_NEIGHBORS, "rng_seed": RNG_SEED,
            "marker_genes": MARKER_GENES,
        }, f, indent=2)
    flag.touch()
    print(f"[{_ts()}] [{comp}] === DONE in {(time.time()-t0)/60:.1f} min ===")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--compartment", choices=list(H5AD), default=None)
    ap.add_argument("--order", nargs="+", default=["CEP", "AF", "NP", "all_cells"],
                    help="Run order; smallest first by default for fail-fast.")
    args = ap.parse_args()

    targets = [args.compartment] if args.compartment else args.order
    print(f"[{_ts()}] harmony sweep targets: {targets}")
    for c in targets:
        try:
            run_compartment(c)
        except Exception as e:
            print(f"[{_ts()}] [{c}] FAILED: {type(e).__name__}: {e}")
            import traceback; traceback.print_exc()
    print(f"[{_ts()}] sweep finished")


if __name__ == "__main__":
    main()
