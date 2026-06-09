#!/usr/bin/env python3
"""
ML#20 — Look for notochordal / progenitor cells in the NP.

Scores NP mesenchymal cells against a notochordal-cell marker panel (and a small
progenitor panel) and asks: is there a notochordal-like sub-population, where does it
sit (which cell type, which study, where on the trajectory), and is it enriched in the
neonatal dataset (GSE189916) as expected?

Needs the tiered annotated AnnData (NP or all_cells; not in the git checkout — lives on
the analysis instance). Run:

    python3 scripts/analysis_ml20_notochordal.py \
        --h5ad data/integrated/.../all_cells_annotated.h5ad --out results/ML20

Reads raw counts from .X (or a 'counts' layer) and CP10K+log1p-normalizes for scoring.
"""
import argparse, os, sys
import numpy as np
import pandas as pd
import scanpy as sc

# Notochordal-cell markers (ML#20) + classic notochordal/NP-progenitor markers.
PANELS = {
    "notochordal":  ["KRT8", "KRT18", "KRT19", "FOXA2", "TBXT", "CD24", "CA12"],
    "progenitor":   ["PROCR", "GDF5", "CD24", "TEK", "ENG"],   # PROCR⁺ NP progenitors (Gan 2021)
}
ALIASES = {"TBXT": ["TBXT", "T", "BRACHYURY"], "T": ["T", "TBXT"]}

CELLTYPE_COLS = ["cell_type", "cell_type_v5", "celltype", "CellType"]
COMP_COLS     = ["compartment", "tissue_compartment"]
STUDY_COLS    = ["study", "study_accession", "batch", "dataset", "sample_study"]
SAMPLE_COLS   = ["sample_id", "sample", "orig.ident"]
PT_COLS       = ["dpt_pseudotime", "dpt", "pseudotime", "X_dpt"]


def first_present(cols, obs):
    return next((c for c in cols if c in obs.columns), None)


def resolve(genes, var_names):
    vn = set(var_names); out = []
    for g in genes:
        if g in vn: out.append(g); continue
        out += [a for a in ALIASES.get(g, []) if a in vn][:1]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--out", default="results/ML20")
    ap.add_argument("--score-quantile", type=float, default=0.95,
                    help="cells above this quantile of the notochordal score are 'noto-high'")
    ap.add_argument("--no-normalize", action="store_true")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    ad = sc.read_h5ad(args.h5ad)
    print(f"[load] {ad.shape[0]:,} cells × {ad.shape[1]:,} genes")

    ct = first_present(CELLTYPE_COLS, ad.obs) or sys.exit("[err] no cell_type column")
    comp = first_present(COMP_COLS, ad.obs)
    study = first_present(STUDY_COLS, ad.obs)
    pt = first_present(PT_COLS, ad.obs)
    print(f"[cols] cell_type={ct} compartment={comp} study={study} pseudotime={pt}")

    # Restrict to NP mesenchymal cells (NP_* labels), the locus of ML#20.
    np_mask = ad.obs[ct].astype(str).str.startswith("NP_")
    if comp is not None:
        np_mask |= ad.obs[comp].astype(str).str.upper().eq("NP")
    sub = ad[np_mask].copy()
    print(f"[subset] {sub.shape[0]:,} NP cells across {sub.obs[ct].nunique()} cell types")

    if not args.no_normalize:
        if "counts" in sub.layers:
            sub.X = sub.layers["counts"].copy()
        if sub.X.max() > 50:
            sc.pp.normalize_total(sub, target_sum=1e4); sc.pp.log1p(sub)
            print("[norm] CP10K + log1p applied")

    # Score panels
    score_cols = []
    for name, genes in PANELS.items():
        g = resolve(genes, sub.var_names)
        miss = sorted(set(genes) - set(g) - {"T"})
        if miss: print(f"[panel:{name}] missing {miss}")
        if g:
            sc.tl.score_genes(sub, g, score_name=f"score_{name}")
            score_cols.append(f"score_{name}")

    # Mean scores per NP cell type
    by_ct = sub.obs.groupby(ct)[score_cols].mean().round(3)
    print("\n[mean scores by NP cell type]\n", by_ct.to_string())
    by_ct.to_csv(os.path.join(args.out, "notochordal_score_by_celltype.csv"))

    # 'Noto-high' cells and where they live
    if "score_notochordal" in score_cols:
        thr = sub.obs["score_notochordal"].quantile(args.score_quantile)
        sub.obs["noto_high"] = sub.obs["score_notochordal"] > thr
        print(f"\n[noto-high] threshold (q{args.score_quantile}) = {thr:.3f}; "
              f"{int(sub.obs['noto_high'].sum()):,} cells")
        print("\n[noto-high fraction by cell type]\n",
              sub.obs.groupby(ct)["noto_high"].mean().round(4).to_string())
        if study is not None:
            tab = sub.obs.groupby(study)["noto_high"].agg(["mean", "sum", "size"]).round(4)
            print("\n[noto-high by study — expect GSE189916/neonatal enriched]\n", tab.to_string())
            tab.to_csv(os.path.join(args.out, "notohigh_by_study.csv"))
        # marker-by-celltype dotplot
        try:
            mk = resolve(PANELS["notochordal"], sub.var_names)
            sc.pl.dotplot(sub, mk, groupby=ct, standard_scale="var", show=False, save="_ml20.png")
        except Exception as e:
            print(f"[plot] dotplot skipped: {e}")

    # Relationship to trajectory, if pseudotime is present
    if pt is not None and "score_notochordal" in score_cols:
        from scipy.stats import spearmanr
        d = sub.obs[[pt, "score_notochordal"]].dropna()
        if len(d) > 100:
            rho, p = spearmanr(d[pt], d["score_notochordal"])
            print(f"\n[trajectory] Spearman(notochordal score, {pt}) rho={rho:.3f}, p={p:.2e} "
                  f"(negative rho = notochordal signal toward the trajectory root)")

    print("\n[done] outputs in", args.out)


if __name__ == "__main__":
    main()
