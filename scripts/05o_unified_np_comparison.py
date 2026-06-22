#!/usr/bin/env python3
"""Module 05o — Unified NP integration-method comparison.

Builds ONE NP comparison table with directly-comparable metrics across seven
method variants, plus a uniform NP UMAP grid, by scoring every embedding through
the SAME metric battery (`compute_metrics` from 05n_harmony_integration.py, which
is byte-identical to the battery that produced the published comparison_table.tsv
CCA rows). Only scANVI and STACAS are newly computed here; the four Seurat CCA
rows and the Harmony row are reused as-is because they were produced by the
identical battery.

The seven variants (scope in parentheses):
  1. Flat CCA v5            (all)          — reuse comparison_table.tsv baseline_flat_v5
  2. Flat CCA v4            (all)          — reuse comparison_table.tsv flat_v4
  3. Tiered CCA v5          (mesenchymal)  — reuse comparison_table.tsv tiered_v5
  4. Tiered CCA v4          (mesenchymal)  — reuse comparison_table.tsv tiered_v4
  5. scANVI                 (all)          — NEW: flat whole-NP, 3000 HVG, 20-latent
  6. STACAS                 (all)          — NEW: full-NP, scored from 05c+05e export
  7. Harmony                (all)          — reuse harmony/NP/metrics.tsv

Why scANVI/STACAS were not comparable before: they were only ever scored under an
unnormalized-LISI convention (results/integration/workflow_comparison.tsv: raw
iLISI ~3.7, raw silhouette, flipped sign) and STACAS used only 16k cells. Their
embeddings no longer existed on disk, so they are regenerated here on the full NP
set (the same data/integrated/tiered_v4/NP.h5ad set Harmony used).

Stages (run independently; the scanvi stage is the multi-hour CPU step):
  --stage scanvi        Train flat whole-NP scANVI (CPU), score it, persist embedding.
  --stage stacas-score  Score the STACAS export (produced by 05c + 05e) on the battery.
  --stage assemble      Concatenate reused + new rows -> unified_comparison.tsv (+ .md).
  --stage figure        Render the NP UMAP grid (fig14).
  --stage all           stacas-score -> assemble -> figure (NOT scanvi; run that first).

Usage:
  CUDA_VISIBLE_DEVICES="" python3 scripts/05o_unified_np_comparison.py --stage scanvi
  python3 scripts/05o_unified_np_comparison.py --stage stacas-score
  python3 scripts/05o_unified_np_comparison.py --stage assemble
  python3 scripts/05o_unified_np_comparison.py --stage figure
"""
from __future__ import annotations
import argparse, importlib.util, os, sys, time, warnings
from pathlib import Path
import numpy as np, pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

BASE = Path(__file__).resolve().parent.parent
SCRIPTS = BASE / "scripts"

# ── Inputs ──────────────────────────────────────────────────────────────────
NP_H5AD = BASE / "data" / "integrated" / "tiered_v4" / "NP.h5ad"   # raw counts in .X
COMPARISON_TABLE = BASE / "results" / "integration" / "np_experiment" / "comparison_table.tsv"
HARMONY_METRICS = BASE / "results" / "integration" / "harmony" / "NP" / "metrics.tsv"
STACAS_EXPORT = BASE / "data" / "integrated" / "stacas" / "embeddings_export" / "NP_embedding.csv.gz"

# CCA / tiered embeddings for the UMAP grid
EMB_FILES = {
    "Flat CCA (v5)":      BASE / "data" / "integrated" / "cca" / "bridge_export" / "NP" / "embedding_integrated.cca.csv.gz",
    "Tiered CCA (v5)":    BASE / "data" / "integrated" / "np_experiment" / "tiered_v5" / "mesenchymal" / "embedding_integrated.cca.csv.gz",
    "Tiered CCA (v4)":    BASE / "data" / "integrated" / "np_experiment" / "tiered_v4" / "mesenchymal" / "embedding_pca.csv.gz",
    "Harmony":            BASE / "results" / "integration" / "harmony" / "NP" / "embedding_harmony.npy",
    # scANVI / STACAS filled from the regenerated outputs below
}

# ── Outputs ─────────────────────────────────────────────────────────────────
SCANVI_OUT = BASE / "results" / "integration" / "scanvi_np_flat"
STACAS_OUT = BASE / "results" / "integration" / "stacas_np_flat"
UNIFIED_TSV = BASE / "results" / "integration" / "np_experiment" / "unified_comparison.tsv"
UNIFIED_MD = BASE / "results" / "integration" / "np_experiment" / "unified_comparison.md"
FIG_DIR = BASE / "docs" / "manuscript_figures"
FIG_PATH = FIG_DIR / "fig14_np_integration_umap_grid.png"

# ── Parameters (kept in lockstep with 05n / 05h) ──────────────────────────────
HVG_N = 3000
SCANVI_N_LATENT = 20
SCVI_MAX_EPOCHS = 150
SCANVI_MAX_EPOCHS = 30
RNG_SEED = 42
MARKER_GENES = ["COL2A1", "ACAN", "SOX9", "COL1A1"]
METRIC_COLS = ["iLISI", "batch_ASW", "cLISI", "bio_ASW", "condition_ASW",
               "condition_LISI", "NMI", "ARI",
               "var_ratio_ACAN", "var_ratio_COL1A1", "var_ratio_COL2A1", "var_ratio_SOX9"]

# Display order + scope for the unified table
ROW_ORDER = [
    ("Flat CCA (v5)",   "all"),
    ("Flat CCA (v4)",   "all"),
    ("Tiered CCA (v5)", "mesenchymal"),
    ("Tiered CCA (v4)", "mesenchymal"),
    ("scANVI",          "all"),
    ("STACAS",          "all"),
    ("Harmony",         "all"),
]


def _ts():
    return time.strftime("%H:%M:%S")


def _load_module(path: Path, name: str):
    """Import a script by file path (handles leading-digit module names)."""
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _get_compute_metrics():
    return _load_module(SCRIPTS / "05n_harmony_integration.py", "harmony05n").compute_metrics


# ═════════════════════════════════════════════════════════════════════════════
# Shared: load NP expression context (log1p CP10K full-gene matrix + obs)
# ═════════════════════════════════════════════════════════════════════════════

def load_np_context():
    """Return (ad_lognorm, gene_names) for tiered_v4/NP.h5ad.

    ad_lognorm.X is log1p(CP10K) over ALL genes (for marker var_ratios), and
    ad_lognorm.obs carries study / coarse_label / condition_harmonized. Cell order
    is the canonical NP order used to align embeddings by barcode.
    """
    print(f"[{_ts()}] loading NP context {NP_H5AD}")
    ad = sc.read_h5ad(NP_H5AD)
    ad.obs_names = ad.obs_names.astype(str)
    # .X is raw counts -> log1p CP10K (same as 05n's counts_norm_log)
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    print(f"[{_ts()}]   NP context: {ad.shape[0]:,} cells x {ad.shape[1]:,} genes")
    return ad, list(ad.var_names)


def _align_obs_and_counts(ad_ctx, barcodes):
    """Reindex the NP context to a given barcode order; report coverage."""
    idx = ad_ctx.obs_names
    pos = pd.Index(idx).get_indexer(pd.Index(barcodes))
    found = pos >= 0
    cov = float(found.mean())
    if cov < 1.0:
        print(f"[{_ts()}]   WARNING: barcode join coverage {cov:.4f} "
              f"({(~found).sum():,} of {len(barcodes):,} not found)")
    pos = pos[found]
    sub = ad_ctx[pos]
    return sub.obs.copy(), sub.X, found


def _score_embedding(embedding, obs, counts_norm_log, gene_names):
    compute_metrics = _get_compute_metrics()
    return compute_metrics(embedding, obs, counts_norm_log, gene_names, do_var=True)


# ═════════════════════════════════════════════════════════════════════════════
# STAGE: scanvi  — flat whole-NP scANVI on CPU, then score
# ═════════════════════════════════════════════════════════════════════════════

def stage_scanvi(args):
    SCANVI_OUT.mkdir(parents=True, exist_ok=True)
    emb_path = SCANVI_OUT / "embedding_scanvi.npy"
    done = SCANVI_OUT / "DONE"
    if done.exists() and not args.force:
        print(f"[{_ts()}] scANVI DONE exists — skip (use --force to retrain)")
        return

    os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")  # force CPU
    from scipy import sparse
    mod05b = _load_module(SCRIPTS / "05b_integration_scanvi.py", "scanvi05b")

    print(f"[{_ts()}] loading {NP_H5AD}")
    ad = sc.read_h5ad(NP_H5AD)
    ad.obs_names = ad.obs_names.astype(str)
    print(f"[{_ts()}]   shape={ad.shape}, X dtype={ad.X.dtype}")

    # counts layer (raw integer) for scVI / seurat_v3 HVG
    counts = ad.X.copy()
    if not sparse.issparse(counts):
        counts = sparse.csr_matrix(counts)
    counts.data = np.rint(counts.data).astype(np.float32)
    ad.layers["counts"] = counts.astype(np.int32)

    # batch-aware HVG on raw counts, then SUBSET (matches the 3000-HVG convention
    # used by Flat CCA / Harmony, and keeps scVI tractable on CPU)
    print(f"[{_ts()}] HVG seurat_v3 n_top_genes={HVG_N} batch_key=study (subset)")
    sc.pp.highly_variable_genes(ad, n_top_genes=HVG_N, flavor="seurat_v3",
                                batch_key="study", layer="counts", subset=True)
    # log-norm .X (scANVI uses the counts layer; this keeps .X sane for any use)
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    print(f"[{_ts()}]   subset to {ad.shape[1]} HVGs")

    print(f"[{_ts()}] flat scANVI: scVI {SCVI_MAX_EPOCHS} + scANVI {SCANVI_MAX_EPOCHS} epochs (CPU)")
    _, emb_key = mod05b.run_scanvi(
        ad, batch_key="study", labels_key="coarse_label",
        unlabeled_category="Unknown", n_latent=SCANVI_N_LATENT,
        scvi_max_epochs=SCVI_MAX_EPOCHS, scanvi_max_epochs=SCANVI_MAX_EPOCHS,
        model_dir=SCANVI_OUT / "models", label="NP_flat",
    )
    emb = np.asarray(ad.obsm[emb_key])
    np.save(emb_path, emb)
    pd.DataFrame({"barcode": ad.obs_names.astype(str)}).to_csv(
        SCANVI_OUT / "cell_index.csv.gz", index=False)
    umap_key = f"X_umap_scanvi_NP_flat"
    if umap_key in ad.obsm:
        np.save(SCANVI_OUT / "umap_scanvi.npy", np.asarray(ad.obsm[umap_key]))
    print(f"[{_ts()}] scANVI embedding saved: {emb.shape}")

    # Score immediately (obs + full-gene log-norm are needed; reload context to get all genes)
    _score_and_write("scANVI", "all", emb, ad.obs_names.astype(str).tolist(), SCANVI_OUT)
    done.touch()
    print(f"[{_ts()}] scANVI stage DONE")


# ═════════════════════════════════════════════════════════════════════════════
# STAGE: stacas-score — score the STACAS export from 05c + 05e
# ═════════════════════════════════════════════════════════════════════════════

def _read_embedding_csv(path: Path):
    """Read an embedding csv.gz -> (matrix[n,d], barcodes[n]).

    Handles both layouts seen in this repo: numeric dims + study/condition/barcode
    columns (05e export), or an index column of barcodes + numeric dims.
    """
    df = pd.read_csv(path)
    meta_cols = [c for c in df.columns if c in
                 ("study", "condition_harmonized", "cell_barcode", "barcode")]
    bc_col = next((c for c in ("cell_barcode", "barcode") if c in df.columns), None)
    if bc_col is not None:
        barcodes = df[bc_col].astype(str).values
    else:
        # first column is the index/barcode
        barcodes = df.iloc[:, 0].astype(str).values
        meta_cols = [df.columns[0]] + meta_cols
    num = df.drop(columns=[c for c in meta_cols if c in df.columns])
    num = num.select_dtypes(include=[np.number])
    return num.values.astype(np.float32), barcodes


def stage_stacas_score(args):
    STACAS_OUT.mkdir(parents=True, exist_ok=True)
    if not STACAS_EXPORT.exists():
        print(f"[{_ts()}] STACAS export not found: {STACAS_EXPORT}")
        print("  Run first:  Rscript scripts/05c_integration_stacas.R --object NP --no-downsample --force")
        print("  then:       Rscript scripts/05e_export_rds_embeddings.R --workflow stacas")
        sys.exit(2)

    emb, barcodes = _read_embedding_csv(STACAS_EXPORT)
    print(f"[{_ts()}] STACAS export: {emb.shape[0]:,} cells x {emb.shape[1]} dims")
    _score_and_write("STACAS", "all", emb, list(barcodes), STACAS_OUT)
    print(f"[{_ts()}] STACAS scoring DONE")


def _score_and_write(method, scope, emb, barcodes, out_dir):
    """Align to NP context, score, write per-method metrics.tsv."""
    ad_ctx, gene_names = load_np_context()
    obs, counts_norm_log, found = _align_obs_and_counts(ad_ctx, barcodes)
    emb = np.asarray(emb)[found]
    print(f"[{_ts()}] scoring {method}: {emb.shape[0]:,} cells aligned")
    metrics = _score_embedding(emb, obs.reset_index(drop=True), counts_norm_log, gene_names)
    row = {"method": method, "scope": scope, "n_cells": int(emb.shape[0]),
           "n_dims": int(emb.shape[1]), **metrics}
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([row]).to_csv(out_dir / "metrics.tsv", sep="\t", index=False)
    print(f"[{_ts()}]   wrote {out_dir/'metrics.tsv'}")
    for k in METRIC_COLS:
        if k in metrics:
            print(f"      {k}={metrics[k]:.4f}")


# ═════════════════════════════════════════════════════════════════════════════
# STAGE: assemble — concatenate reused + new rows into the unified table
# ═════════════════════════════════════════════════════════════════════════════

def stage_assemble(args):
    rows = []

    # 1-4: reused Seurat CCA rows from comparison_table.tsv
    cca = pd.read_csv(COMPARISON_TABLE, sep="\t")
    cca_map = {
        ("Flat CCA (v5)", "all"):         ("baseline_flat_v5", "all"),
        ("Flat CCA (v4)", "all"):         ("flat_v4", "all"),
        ("Tiered CCA (v5)", "mesenchymal"): ("tiered_v5", "mesenchymal"),
        ("Tiered CCA (v4)", "mesenchymal"): ("tiered_v4", "mesenchymal"),
    }
    for (method, scope), (run, sc_scope) in cca_map.items():
        sub = cca[(cca["run"] == run) & (cca["scope"] == sc_scope)]
        if sub.empty:
            print(f"[{_ts()}] WARNING: {run}/{sc_scope} missing from comparison_table.tsv")
            continue
        r = sub.iloc[0]
        rows.append(_row(method, scope, r))

    # 7: reused Harmony row
    if HARMONY_METRICS.exists():
        h = pd.read_csv(HARMONY_METRICS, sep="\t").iloc[0]
        rows.append(_row("Harmony", "all", h))
    else:
        print(f"[{_ts()}] WARNING: harmony metrics missing {HARMONY_METRICS}")

    # 5-6: new scANVI / STACAS rows
    for method, out_dir in [("scANVI", SCANVI_OUT), ("STACAS", STACAS_OUT)]:
        mpath = out_dir / "metrics.tsv"
        if mpath.exists():
            rows.append(_row(method, "all", pd.read_csv(mpath, sep="\t").iloc[0]))
        else:
            print(f"[{_ts()}] WARNING: {method} metrics missing {mpath} — run its stage first")

    df = pd.DataFrame(rows)
    # order by ROW_ORDER
    order = {m: i for i, (m, _) in enumerate(ROW_ORDER)}
    df["__o"] = df["method"].map(order)
    df = df.sort_values("__o").drop(columns="__o").reset_index(drop=True)

    UNIFIED_TSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(UNIFIED_TSV, sep="\t", index=False)
    print(f"[{_ts()}] wrote {UNIFIED_TSV} ({len(df)} rows)")

    _write_markdown(df)
    print(df.to_string(index=False))


def _row(method, scope, r):
    """Extract the canonical metric columns from a source row (align by NAME)."""
    out = {"method": method, "scope": scope,
           "n_cells": int(r["n_cells"]) if "n_cells" in r and pd.notna(r["n_cells"]) else np.nan}
    for k in METRIC_COLS:
        out[k] = float(r[k]) if k in r.index and pd.notna(r[k]) else np.nan
    return out


def _write_markdown(df):
    # Manuscript-facing table: friendly marker headers, 3 decimals.
    # condition_ASW / condition_LISI surface the condition-signal axis that drove
    # the tiered-v4 selection in the §5 NP experiment (notebook 05, §5d). Direction
    # note: condition_ASW higher (closer to 0) = better; condition_LISI lower = better
    # — the only displayed column where lower is preferable.
    disp = ["iLISI", "batch_ASW", "cLISI", "bio_ASW", "condition_ASW", "condition_LISI",
            "NMI", "ARI",
            "var_ratio_ACAN", "var_ratio_COL2A1", "var_ratio_SOX9", "var_ratio_COL1A1"]
    hdr = ["Method", "Scope", "iLISI", "batch_ASW", "cLISI", "bio_ASW",
           "condition_ASW", "condition_LISI", "NMI", "ARI",
           "ACAN", "COL2A1", "SOX9", "COL1A1"]
    lines = ["| " + " | ".join(hdr) + " |",
             "|" + "|".join(["---"] * len(hdr)) + "|"]
    for _, r in df.iterrows():
        cells = [r["method"], r["scope"]] + [
            ("%.3f" % r[c]) if pd.notna(r[c]) else "—" for c in disp]
        lines.append("| " + " | ".join(str(x) for x in cells) + " |")
    UNIFIED_MD.write_text("\n".join(lines) + "\n")
    print(f"[{_ts()}] wrote {UNIFIED_MD}")


# ═════════════════════════════════════════════════════════════════════════════
# STAGE: figure — uniform NP UMAP grid (study + coarse_label), 6 methods
# ═════════════════════════════════════════════════════════════════════════════

def _umap_for(emb, barcodes, ad_ctx, precomputed=None):
    """Return (umap[n,2], obs) for an embedding aligned to NP context."""
    pos = pd.Index(ad_ctx.obs_names).get_indexer(pd.Index(barcodes))
    found = pos >= 0
    pos = pos[found]
    obs = ad_ctx.obs.iloc[pos].reset_index(drop=True)
    if precomputed is not None and precomputed.shape[0] == len(barcodes):
        return precomputed[found], obs
    a = sc.AnnData(X=np.zeros((found.sum(), 1), dtype=np.float32))
    a.obsm["X_emb"] = np.asarray(emb)[found]
    sc.pp.neighbors(a, use_rep="X_emb", n_neighbors=15, random_state=RNG_SEED)
    sc.tl.umap(a, random_state=RNG_SEED)
    return a.obsm["X_umap"], obs


def stage_figure(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # methods to plot (omit flat_v4 — embedding gone)
    methods = ["Flat CCA (v5)", "Tiered CCA (v5)", "Tiered CCA (v4)",
               "scANVI", "STACAS", "Harmony"]
    emb_files = dict(EMB_FILES)
    emb_files["scANVI"] = SCANVI_OUT / "embedding_scanvi.npy"
    emb_files["STACAS"] = STACAS_EXPORT

    ad_ctx, _ = load_np_context()

    panels = {}
    for m in methods:
        p = emb_files.get(m)
        if p is None or not Path(p).exists():
            print(f"[{_ts()}] skip {m}: embedding not found ({p})")
            continue
        if str(p).endswith(".npy"):
            emb = np.load(p)
            bc = pd.read_csv(Path(p).parent / "cell_index.csv.gz")["barcode"].astype(str).values
        else:
            emb, bc = _read_embedding_csv(Path(p))
        print(f"[{_ts()}] UMAP for {m}: {emb.shape}")
        umap, obs = _umap_for(emb, bc, ad_ctx)
        panels[m] = (umap, obs)

    if not panels:
        print(f"[{_ts()}] no panels available — nothing to render")
        return

    cols = [m for m in methods if m in panels]
    ncols = len(cols)
    fig, axes = plt.subplots(2, ncols, figsize=(3.2 * ncols, 6.8), squeeze=False)

    # consistent study palette across panels
    studies = sorted(ad_ctx.obs["study"].astype(str).unique())
    study_cmap = {s: plt.cm.tab20(i % 20) for i, s in enumerate(studies)}
    labels = sorted(ad_ctx.obs["coarse_label"].astype(str).unique())
    label_cmap = {l: plt.cm.tab10(i % 10) for i, l in enumerate(labels)}

    for j, m in enumerate(cols):
        umap, obs = panels[m]
        for row, (key, cmap) in enumerate([("study", study_cmap), ("coarse_label", label_cmap)]):
            ax = axes[row][j]
            vals = obs[key].astype(str).values
            colors = np.array([cmap[v] for v in vals])
            order = np.random.RandomState(RNG_SEED).permutation(len(umap))
            ax.scatter(umap[order, 0], umap[order, 1], c=colors[order], s=1.0,
                       linewidths=0, rasterized=True)
            ax.set_xticks([]); ax.set_yticks([])
            if row == 0:
                ax.set_title(m, fontsize=11)
            if j == 0:
                ax.set_ylabel("by study" if key == "study" else "by coarse type",
                              fontsize=10)

    # legends
    from matplotlib.lines import Line2D
    study_handles = [Line2D([0], [0], marker="o", linestyle="", markersize=5,
                            markerfacecolor=study_cmap[s], label=s) for s in studies]
    label_handles = [Line2D([0], [0], marker="o", linestyle="", markersize=5,
                            markerfacecolor=label_cmap[l], label=l) for l in labels]
    axes[0][-1].legend(handles=study_handles, fontsize=6, loc="center left",
                       bbox_to_anchor=(1.01, 0.5), frameon=False)
    axes[1][-1].legend(handles=label_handles, fontsize=6, loc="center left",
                       bbox_to_anchor=(1.01, 0.5), frameon=False)

    fig.suptitle("NP integration methods — UMAP of each embedding", fontsize=13)
    fig.tight_layout(rect=[0, 0, 0.92, 0.97])
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_PATH, dpi=150, bbox_inches="tight")
    print(f"[{_ts()}] wrote {FIG_PATH}")


# ═════════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True,
                    choices=["scanvi", "stacas-score", "assemble", "figure", "all"])
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if args.stage == "scanvi":
        stage_scanvi(args)
    elif args.stage == "stacas-score":
        stage_stacas_score(args)
    elif args.stage == "assemble":
        stage_assemble(args)
    elif args.stage == "figure":
        stage_figure(args)
    elif args.stage == "all":
        stage_stacas_score(args)
        stage_assemble(args)
        stage_figure(args)


if __name__ == "__main__":
    main()
