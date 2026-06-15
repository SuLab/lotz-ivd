#!/usr/bin/env python3
"""Module 05o — Unified NP integration-method comparison (single metric convention).

Scores ALL six integration methods for the NP compartment through ONE identical
metric battery, so the manuscript's "Integration-method comparison" table is
apples-to-apples across every row. This replaces the previous situation where:

  * flat/tiered CCA (v4, v5) were scored by 05h via the scib normalized battery
    (ilisi_knn scale=True, silhouette_batch rescale=True, 50k/30k subsample, k=90);
  * scANVI + STACAS were originally scored by 05d via a *homemade* inverse-Simpson
    LISI and a raw sklearn batch silhouette on a 5k/k=30 subsample — a DIFFERENT,
    unnormalized convention that is NOT comparable to the scib numbers;
  * Harmony was scored by 05n, which already uses the 05h scib convention.

To guarantee one convention, this script IMPORTS the exact metric functions from
05h_np_experiment_metrics.py rather than re-implementing them. Every method below
is loaded into the same (embedding, metadata, counts, gene_names) form and run
through the same compute_* functions with the same constants (MARKER_GENES,
MAX_CELLS_LISI/ASW, KNN_NEIGHBORS, RNG_SEED).

Methods (per the manuscript table; tiered methods use the MESENCHYMAL tier only —
the continuum-relevant population the marker-variance argument is about):

  flat_cca_v5    bridge  data/integrated/cca/bridge_export/NP
  tiered_cca_v5  bridge  data/integrated/np_experiment/tiered_v5/mesenchymal
  flat_cca_v4    bridge  data/integrated/np_experiment/flat_v4/all
  tiered_cca_v4  bridge  data/integrated/np_experiment/tiered_v4/mesenchymal
  scanvi         h5ad    data/integrated/scanvi/NP.h5ad           (obsm X_scanvi_*)
  harmony        npy     results/integration/harmony/NP/embedding_harmony.npy
                         (+ obs/counts from data/integrated/tiered_v4/NP.h5ad)

FAIRNESS CAVEAT (must be footnoted in the manuscript): the flat methods score all
NP cells with all coarse labels; the tiered rows score only the mesenchymal tier
(fewer label classes, fewer cells). bio_ASW / cLISI / NMI / ARI and the
marker-variance ratios are therefore computed on different cell populations for
flat vs. tiered rows. This is intrinsic to comparing flat vs. tiered integration
and is not a bug. Batch-mixing metrics (iLISI, batch_ASW) remain interpretable
within each row.

Any method whose inputs are missing on disk is SKIPPED with a loud warning and a
status='missing' row — the sweep never silently drops a method.

Usage:
    python3 scripts/05o_unified_np_comparison.py
    python3 scripts/05o_unified_np_comparison.py --method scanvi harmony
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread  # noqa: F401  (kept for parity / optional use)
from scipy.sparse import csc_matrix, issparse

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

# ── Paths ──────────────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent.parent
INT = BASE / "data" / "integrated"
RESULTS_DIR = BASE / "results" / "integration" / "np_unified_comparison"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ── Import 05h's metric battery as the single source of truth ───────────────
_05H_PATH = BASE / "scripts" / "05h_np_experiment_metrics.py"


def _load_05h():
    spec = importlib.util.spec_from_file_location("m05h", _05H_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load metric module from {_05H_PATH}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # does not run main() (guarded by __main__)
    return mod


m05h = _load_05h()

# Coarse cell-type label column: try these in order when normalizing metadata.
LABEL_CANDIDATES = ["coarse_label", "coarse_cell_type", "cell_class"]


# ═════════════════════════════════════════════════════════════════════════
# METADATA NORMALIZATION
# ═════════════════════════════════════════════════════════════════════════

def _normalize_metadata(meta: pd.DataFrame) -> pd.DataFrame:
    """Ensure a 'coarse_label' column exists, falling back to coarse_cell_type /
    cell_class. Leaves 'study' and 'condition_harmonized' untouched (the metric
    functions tolerate their absence and emit NaN)."""
    meta = meta.copy()
    if "coarse_label" not in meta.columns:
        for cand in LABEL_CANDIDATES[1:]:
            if cand in meta.columns:
                print(f"    metadata: using '{cand}' as coarse_label")
                meta["coarse_label"] = meta[cand].astype(str).values
                break
    return meta


def _as_cells_by_genes(counts, n_cells: int):
    """Return counts oriented as (cells, genes) csc, or None."""
    if counts is None:
        return None
    if not issparse(counts):
        counts = csc_matrix(counts)
    if counts.shape[0] != n_cells and counts.shape[1] == n_cells:
        counts = counts.T
    return csc_matrix(counts)


# ═════════════════════════════════════════════════════════════════════════
# LOADERS  — each returns (embedding, metadata, counts, gene_names)
# ═════════════════════════════════════════════════════════════════════════

def load_bridge(bridge_dir: Path):
    """CCA runs: reuse 05h's bridge loader verbatim (same file format)."""
    emb, meta, counts, genes = m05h.load_bridge_data(bridge_dir)
    meta = _normalize_metadata(meta)
    return emb, meta, counts, genes


def load_scanvi_h5ad(h5ad_path: Path):
    """scANVI: read the integrated h5ad, auto-detect the latent obsm key."""
    ad = sc.read_h5ad(h5ad_path)
    # Prefer a flat/whole-object latent over a tier-specific one.
    keys = list(ad.obsm.keys())
    pref = [k for k in keys if k.lower().startswith("x_scanvi")
            and "mesenchymal" not in k.lower()]
    pref += [k for k in keys if k.lower() in ("x_scanvi", "x_scvi")]
    pref += [k for k in keys if k.lower().startswith("x_scanvi")]
    if not pref:
        raise KeyError(f"No scANVI latent obsm key in {h5ad_path} (have {keys})")
    emb_key = pref[0]
    print(f"    scANVI embedding obsm key: {emb_key}")
    emb = np.asarray(ad.obsm[emb_key], dtype=np.float32)
    # Drop padding rows (tier embeddings are NaN-filled for the other tier).
    finite = np.isfinite(emb).all(axis=1)
    if not finite.all():
        print(f"    dropping {int((~finite).sum()):,} NaN-padded rows from {emb_key}")
        ad = ad[finite].copy()
        emb = emb[finite]
    meta = _normalize_metadata(ad.obs.reset_index(drop=False).rename(columns={"index": "barcode"}))
    counts = ad.layers["counts"] if "counts" in ad.layers else ad.X
    counts = _as_cells_by_genes(counts, emb.shape[0])
    genes = list(ad.var_names)
    return emb, meta, counts, genes


def load_harmony_npy(emb_npy: Path, cell_index: Path, source_h5ad: Path):
    """Harmony: embedding from 05n's .npy, obs/counts from the source h5ad,
    aligned by barcode (cell_index order == embedding row order)."""
    emb = np.load(emb_npy).astype(np.float32)
    bc = pd.read_csv(cell_index)["barcode"].astype(str).tolist()
    if len(bc) != emb.shape[0]:
        raise ValueError(f"cell_index ({len(bc)}) != embedding rows ({emb.shape[0]})")
    ad = sc.read_h5ad(source_h5ad)
    ad.obs_names = ad.obs_names.astype(str)
    missing = [b for b in bc[:5] if b not in set(ad.obs_names)]
    if missing:
        raise KeyError(f"Harmony barcodes not found in {source_h5ad} (e.g. {missing})")
    ad = ad[bc].copy()  # reorder to embedding order
    meta = _normalize_metadata(ad.obs.reset_index(drop=False).rename(columns={"index": "barcode"}))
    counts = ad.layers["counts"] if "counts" in ad.layers else ad.X
    counts = _as_cells_by_genes(counts, emb.shape[0])
    genes = list(ad.var_names)
    return emb, meta, counts, genes


# ═════════════════════════════════════════════════════════════════════════
# METHOD REGISTRY  (order = display order in the table)
# ═════════════════════════════════════════════════════════════════════════

METHODS = [
    {"key": "flat_cca_v5", "label": "Flat CCA (v5)", "scope": "all",
     "loader": lambda: load_bridge(INT / "cca" / "bridge_export" / "NP")},
    {"key": "tiered_cca_v5", "label": "Tiered CCA (v5)", "scope": "mesenchymal",
     "loader": lambda: load_bridge(INT / "np_experiment" / "tiered_v5" / "mesenchymal")},
    {"key": "flat_cca_v4", "label": "Flat CCA (v4)", "scope": "all",
     "loader": lambda: load_bridge(INT / "np_experiment" / "flat_v4" / "all")},
    {"key": "tiered_cca_v4", "label": "Tiered CCA (v4)", "scope": "mesenchymal",
     "loader": lambda: load_bridge(INT / "np_experiment" / "tiered_v4" / "mesenchymal")},
    {"key": "scanvi", "label": "scANVI", "scope": "all",
     "loader": lambda: load_scanvi_h5ad(INT / "scanvi" / "NP.h5ad")},
    {"key": "harmony", "label": "Harmony", "scope": "all",
     "loader": lambda: load_harmony_npy(
         BASE / "results" / "integration" / "harmony" / "NP" / "embedding_harmony.npy",
         BASE / "results" / "integration" / "harmony" / "NP" / "cell_index.csv.gz",
         INT / "tiered_v4" / "NP.h5ad")},
]


# ═════════════════════════════════════════════════════════════════════════
# SCORING — reuse 05h's metric functions verbatim
# ═════════════════════════════════════════════════════════════════════════

def score_method(spec) -> dict:
    key, label, scope = spec["key"], spec["label"], spec["scope"]
    print(f"\n{'='*64}\n  {label}  [{key} / {scope}]\n{'='*64}")
    row = {"method": key, "label": label, "scope": scope, "status": "complete"}
    try:
        emb, meta, counts, genes = spec["loader"]()
    except Exception as e:
        print(f"  SKIP ({type(e).__name__}): {e}")
        row["status"] = "missing"
        return row

    row["n_cells"] = int(emb.shape[0])
    row["n_dims"] = int(emb.shape[1])
    print(f"  loaded {emb.shape[0]:,} cells x {emb.shape[1]} dims")

    print("  --- batch removal ---")
    row.update(m05h.compute_batch_metrics(emb, meta))
    print("  --- biological preservation ---")
    row.update(m05h.compute_bio_metrics(emb, meta))
    print("  --- condition signal ---")
    row.update(m05h.compute_condition_metrics(emb, meta))
    print("  --- cluster agreement ---")
    row.update(m05h.compute_cluster_agreement(emb, meta))
    print("  --- marker variance ---")
    row.update(m05h.compute_marker_variance(counts, genes, emb, meta))
    return row


# ═════════════════════════════════════════════════════════════════════════
# OUTPUT
# ═════════════════════════════════════════════════════════════════════════

METRIC_COLS = ["iLISI", "batch_ASW", "cLISI", "bio_ASW",
               "condition_ASW", "condition_LISI", "NMI", "ARI",
               "var_ratio_ACAN", "var_ratio_COL2A1",
               "var_ratio_SOX9", "var_ratio_COL1A1"]


def write_outputs(rows: list[dict]):
    df = pd.DataFrame(rows)
    id_cols = ["method", "label", "scope", "status", "n_cells", "n_dims"]
    ordered = [c for c in id_cols + METRIC_COLS if c in df.columns]
    df = df[ordered]

    tsv = RESULTS_DIR / "comparison_table.tsv"
    df.to_csv(tsv, sep="\t", index=False, float_format="%.4f")
    print(f"\n  Saved: {tsv}")

    # Markdown table (manuscript-ready: the chondrogenic/fibrogenic columns
    # renamed to bare gene names to match the existing manuscript layout).
    rename = {"var_ratio_ACAN": "ACAN", "var_ratio_COL2A1": "COL2A1",
              "var_ratio_SOX9": "SOX9", "var_ratio_COL1A1": "COL1A1"}
    md_cols = ["label", "iLISI", "batch_ASW", "cLISI", "bio_ASW", "NMI", "ARI",
               "var_ratio_ACAN", "var_ratio_COL2A1", "var_ratio_SOX9", "var_ratio_COL1A1"]
    md_cols = [c for c in md_cols if c in df.columns]
    disp = df[df["status"] == "complete"][md_cols].rename(columns=rename)
    header = "| Method | " + " | ".join(c for c in disp.columns if c != "label") + " |"
    sep = "|" + "---|" * len(disp.columns)
    lines = [header, sep]
    for _, r in disp.iterrows():
        vals = [r["label"]] + [f"{r[c]:.3f}" if pd.notna(r[c]) else "—"
                               for c in disp.columns if c != "label"]
        lines.append("| " + " | ".join(vals) + " |")
    md = RESULTS_DIR / "comparison_table.md"
    md.write_text("\n".join(lines) + "\n")
    print(f"  Saved: {md}")

    print(f"\n{'='*80}\nUNIFIED NP COMPARISON (status=complete rows)\n{'='*80}")
    print(df.to_string(index=False))
    return df


def main():
    ap = argparse.ArgumentParser(description="Unified NP integration comparison")
    ap.add_argument("--method", nargs="+", default=None,
                    choices=[m["key"] for m in METHODS],
                    help="Subset of methods to score (default: all six).")
    args = ap.parse_args()

    print("=" * 64)
    print("Module 05o — Unified NP integration-method comparison")
    print(f"Metric battery imported from: {_05H_PATH.name}")
    print(f"Started: {datetime.now():%Y-%m-%d %H:%M:%S}")
    print("=" * 64)

    targets = ([m for m in METHODS if m["key"] in set(args.method)]
               if args.method else METHODS)
    rows = [score_method(m) for m in targets]

    write_outputs(rows)

    missing = [r["method"] for r in rows if r.get("status") == "missing"]
    if missing:
        print(f"\n  WARNING: {len(missing)} method(s) missing inputs and skipped: "
              f"{', '.join(missing)}")
    print(f"\nCompleted: {datetime.now():%Y-%m-%d %H:%M:%S}")
    if rows and all(r.get("status") == "missing" for r in rows):
        sys.exit(1)


if __name__ == "__main__":
    main()
