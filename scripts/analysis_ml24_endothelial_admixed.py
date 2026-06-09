#!/usr/bin/env python3
"""
ML#24 / ML#25 / ML#14 — Are the "endothelial-admixed" NP_fibrocartilaginous cells actually
mural cells / adventitial fibroblasts (ML#24); is that fraction enriched in
degenerated vs. healthy NP (ML#25); and how do contamination fractions (RBC +
endothelial_admixed) break down by disease state (ML#14)?

The atlas flags 1,831 NP_fibrocartilaginous cells as `endothelial_admixed`
(contamination panel CD34/EMCN/AQP1). This script tests whether those cells carry a
coherent PERIVASCULAR signature (mural/pericyte or adventitial-fibroblast) rather than a
pure endothelial-doublet signal, by scoring them against marker panels and comparing to
bona fide Endothelial, Pericyte_SMC, and "clean" NP_fibrocartilaginous cells.

This is an EMPIRICAL analysis — it needs the tiered annotated all_cells AnnData
(not in the git checkout; lives on the analysis instance). Run:

    python3 scripts/analysis_ml24_endothelial_admixed.py \
        --h5ad data/integrated/tiered_v4/all_cells_annotated.h5ad \
        --out results/ML24

Interpretation guide (printed at the end):
  * If EA cells score HIGH on Mural/Adventitial and only MODESTLY on Endothelial,
    and express COL1A1/VCAN (fibroblast identity) -> supports ML#24 (genuine
    perivascular/adventitial population, not mere contamination).
  * If EA cells score HIGH on Endothelial markers comparably to bona fide Endothelial
    (co-expressing fibroblast + endothelial at doublet-like levels) -> supports the
    contamination/doublet interpretation the manuscript currently assumes.
"""
import argparse, os, sys
import numpy as np
import pandas as pd
import scanpy as sc

# ── Marker panels ────────────────────────────────────────────────────────────
PANELS = {
    "Endothelial":            ["PECAM1", "CDH5", "VWF", "EMCN", "CLDN5", "FLT1", "CD34", "AQP1"],
    "Mural_pericyte_SMC":     ["RGS5", "PDGFRB", "ACTA2", "NOTCH3", "TAGLN", "MYH11", "MCAM", "KCNJ8", "CSPG4"],
    "Adventitial_fibroblast": ["PI16", "DPT", "MFAP5", "SEMA3C", "PCOLCE2", "ABCA8", "CD34"],
    "NP_fibrocartilaginous":  ["COL1A1", "COL2A1", "VCAN", "DCN", "FN1"],
}

CELLTYPE_COLS = ["cell_type", "cell_type_v5", "celltype", "CellType"]
SUBTYPE_COLS  = ["cell_subtype", "sub_state", "substate", "cell_substate"]
CONTAM_COLS   = ["contamination_type", "contamination", "contam_type"]
COND_COLS     = ["condition_harmonized", "condition", "degeneration_severity"]
COMP_COLS     = ["compartment", "tissue_compartment"]


def first_present(cols, obs):
    for c in cols:
        if c in obs.columns:
            return c
    return None


def resolve_genes(genes, var_names):
    aliases = {"T": ["T", "TBXT"], "TBXT": ["TBXT", "T"]}
    out = []
    vn = set(var_names)
    for g in genes:
        if g in vn:
            out.append(g); continue
        for a in aliases.get(g, []):
            if a in vn:
                out.append(a); break
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True, help="tiered annotated all_cells .h5ad")
    ap.add_argument("--out", default="results/ML24")
    ap.add_argument("--no-normalize", action="store_true",
                    help="skip CP10K+log1p (use if .X is already log-normalized)")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    ad = sc.read_h5ad(args.h5ad)
    print(f"[load] {ad.shape[0]:,} cells × {ad.shape[1]:,} genes")

    ct = first_present(CELLTYPE_COLS, ad.obs)
    if ct is None:
        sys.exit(f"[err] no cell_type column found among {CELLTYPE_COLS}; obs has {list(ad.obs.columns)}")
    sub = first_present(SUBTYPE_COLS, ad.obs)
    contam = first_present(CONTAM_COLS, ad.obs)
    cond = first_present(COND_COLS, ad.obs)
    comp = first_present(COMP_COLS, ad.obs)
    print(f"[cols] cell_type={ct} subtype={sub} contam={contam} condition={cond} compartment={comp}")

    # Identify endothelial-admixed NP_fibrocartilaginous cells (robust to where the flag lives)
    is_npf = ad.obs[ct].astype(str).eq("NP_fibrocartilaginous")
    ea_mask = pd.Series(False, index=ad.obs_names)
    if contam is not None:
        ea_mask |= ad.obs[contam].astype(str).str.contains("endothelial_admixed", case=False, na=False)
    if sub is not None:
        ea_mask |= ad.obs[sub].astype(str).str.contains("endothelial_admixed", case=False, na=False)
    if not ea_mask.any():
        sys.exit("[err] could not locate endothelial_admixed cells in contam/subtype columns")

    groups = pd.Series(index=ad.obs_names, dtype=object)
    groups[is_npf & ea_mask] = "EA_NP_fibro"
    groups[is_npf & ~ea_mask] = "clean_NP_fibro"
    groups[ad.obs[ct].astype(str).eq("Endothelial")] = "Endothelial"
    groups[ad.obs[ct].astype(str).eq("Pericyte_SMC")] = "Pericyte_SMC"
    ad.obs["ml24_group"] = groups
    print("\n[group sizes]\n", ad.obs["ml24_group"].value_counts(dropna=True).to_string())

    # Normalize for scoring
    if not args.no_normalize:
        xmax = ad.X.max()
        if xmax > 50:  # looks like raw counts
            sc.pp.normalize_total(ad, target_sum=1e4)
            sc.pp.log1p(ad)
            print("[norm] applied CP10K + log1p")
        else:
            print(f"[norm] .X max={xmax:.2f}; assuming already log-normalized")

    # Score panels
    score_cols = []
    for name, genes in PANELS.items():
        g = resolve_genes(genes, ad.var_names)
        missing = sorted(set(genes) - set(g))
        if missing:
            print(f"[panel:{name}] missing {missing}")
        if not g:
            print(f"[panel:{name}] no genes present — skipped"); continue
        sc.tl.score_genes(ad, g, score_name=f"score_{name}")
        score_cols.append(f"score_{name}")

    sub_ad = ad[ad.obs["ml24_group"].notna()]
    score_tbl = sub_ad.obs.groupby("ml24_group")[score_cols].mean().round(3)
    print("\n[mean panel scores by group]\n", score_tbl.to_string())
    score_tbl.to_csv(os.path.join(args.out, "panel_scores_by_group.csv"))

    # Per-marker mean expression (dotplot) across groups
    all_markers = resolve_genes(sorted({m for p in PANELS.values() for m in p}), ad.var_names)
    try:
        sc.pl.dotplot(sub_ad, all_markers, groupby="ml24_group", standard_scale="var",
                      show=False, save="_ml24.png")
        print("[plot] dotplot saved (scanpy figures/ dir)")
    except Exception as e:
        print(f"[plot] dotplot skipped: {e}")

    # ML#25 — is the EA fraction enriched in degenerated vs healthy NP?
    if cond is not None:
        npf = ad.obs[is_npf].copy()
        npf["ea"] = ea_mask[is_npf].values
        c = npf[cond].astype(str)
        is_healthy = c.str.contains("healthy", case=False, na=False)
        deg = np.where(is_healthy, "healthy", "degenerated")
        tab = pd.crosstab(deg, npf["ea"])
        print("\n[ML#25] NP_fibrocartilaginous EA fraction by condition:\n", tab.to_string())
        frac = (tab.get(True, 0) / tab.sum(axis=1)).round(4)
        print("\n[ML#25] EA fraction:\n", frac.to_string())
        try:
            from scipy.stats import fisher_exact
            if tab.shape == (2, 2):
                orr, p = fisher_exact(tab.values)
                print(f"[ML#25] Fisher exact OR={orr:.2f}, p={p:.2e}")
        except Exception as e:
            print(f"[ML#25] test skipped: {e}")
        tab.to_csv(os.path.join(args.out, "ea_fraction_by_condition.csv"))

    # ML#14 — full contamination_type × disease state (RBC + endothelial_admixed + clean),
    # atlas-wide and per compartment.
    if contam is not None and cond is not None:
        c = ad.obs[cond].astype(str)
        deg = np.where(c.str.contains("healthy", case=False, na=False), "healthy", "degenerated")
        ad.obs["_cond2"] = deg
        ctab = pd.crosstab(ad.obs[contam].astype(str), ad.obs["_cond2"])
        cfrac = (ctab / ctab.sum(axis=0)).round(4)  # fraction of each condition's cells
        print("\n[ML#14] contamination_type × disease state (atlas-wide counts):\n", ctab.to_string())
        print("\n[ML#14] contamination_type as fraction of each condition:\n", cfrac.to_string())
        ctab.to_csv(os.path.join(args.out, "contamination_by_condition.csv"))
        if comp is not None:
            for cp in sorted(ad.obs[comp].astype(str).unique()):
                m = ad.obs[comp].astype(str).eq(cp)
                sub = pd.crosstab(ad.obs.loc[m, contam].astype(str), ad.obs.loc[m, "_cond2"])
                print(f"\n[ML#14] {cp} contamination_type × disease state:\n", sub.to_string())

    print("\n[done] outputs in", args.out)


if __name__ == "__main__":
    main()
