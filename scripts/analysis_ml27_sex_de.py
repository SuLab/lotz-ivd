#!/usr/bin/env python3
"""
ML#27 — Address sex differences in the DE analysis.

A fully sex-stratified DE (e.g. healthy-vs-degenerated within females) is not powerable
here: sex is confounded with disease state and study, and the female/healthy cell is
near-empty (see the manuscript Caveat 6). This script instead does the rigorous,
feasible thing:

  1. Tabulate sex availability per cell-type × contrast (which contrasts can even
     support a sex term).
  2. Where a sex term is estimable, run pseudobulk DESeq2 twice — design `~group`
     (sex-naive, the manuscript result) and `~sex + group` (sex-adjusted) — and report
     whether the DEGs are robust to sex adjustment (i.e. not sex-confounded).

Needs the tiered annotated all_cells AnnData (analysis instance) + the in-repo
metadata/sample_metadata.tsv. Mirrors scripts/08_differential.py conventions. Run:

    python3 scripts/analysis_ml27_sex_de.py \
        --h5ad data/integrated/.../all_cells_annotated.h5ad --out results/ML27
"""
import argparse, os, sys
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

BASE = Path(__file__).resolve().parent.parent
DEFAULT_META = BASE / "metadata" / "sample_metadata.tsv"

EXCLUDE_SAMPLES = {"GSE205535_NNP"}      # matches 08_differential.py
MIN_SAMPLES_PER_GROUP = 3
MIN_CELLS_PER_SAMPLE = 50
LFC, PADJ = 0.5, 0.05

COMPARISONS = [
    ("healthy_vs_degenerated_all",   "healthy", ["degenerated_mild", "degenerated_severe", "degenerated_ungraded"]),
    ("healthy_vs_degenerated_mild",  "healthy", ["degenerated_mild"]),
    ("healthy_vs_degenerated_severe","healthy", ["degenerated_severe"]),
    ("mild_vs_severe", "degenerated_mild", ["degenerated_severe"]),
]
DEFAULT_CELLTYPES = ["NP_fibrocartilaginous", "NP_fibrochondrocyte_chondroid",
                     "NP_mature_chondrocyte", "AF_outer"]

CT_COLS, SAMPLE_COLS = ["cell_type", "cell_type_v5", "celltype"], ["sample_id", "sample", "orig.ident"]


def first_present(cols, obs):
    return next((c for c in cols if c in obs.columns), None)


def norm_sex(s):
    s = str(s).strip().lower()
    return "M" if s in ("m", "male") else "F" if s in ("f", "female") else "unk"


def pseudobulk(adata, sample_ids, cells_mask):
    """Sum raw counts per sample for the masked cells. Returns genes × samples DataFrame."""
    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    out = {}
    for sid in sample_ids:
        m = cells_mask & (adata.obs["_sample"].values == sid)
        if m.sum() < MIN_CELLS_PER_SAMPLE:
            continue
        sub = X[m]
        v = np.asarray(sub.sum(axis=0)).ravel()
        out[sid] = v
    if not out:
        return None
    return pd.DataFrame(out, index=adata.var_names)


def deseq(counts_df, clinical, design):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    dds = DeseqDataSet(counts=counts_df.T.astype(int), metadata=clinical,
                       design="+".join(design) if len(design) > 1 else design[0],
                       quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["group", "test", "reference"], quiet=True)
    st.summary()
    r = st.results_df.dropna(subset=["padj"])
    return r[(r["padj"] < PADJ) & (r["log2FoldChange"].abs() > LFC)]


def estimable(clinical):
    """Can we add a sex term? Need >=2 of each sex and sex not collinear with group."""
    nM = (clinical["sex"] == "M").sum(); nF = (clinical["sex"] == "F").sum()
    if nM < 2 or nF < 2:
        return False
    # sex must vary within at least one group (else perfectly confounded with group)
    varies = clinical.groupby("group")["sex"].nunique()
    return (varies >= 2).any()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--meta", default=str(DEFAULT_META))
    ap.add_argument("--out", default="results/ML27")
    ap.add_argument("--cell-types", nargs="*", default=DEFAULT_CELLTYPES)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    adata = sc.read_h5ad(args.h5ad)
    if "counts" not in adata.layers and adata.X.max() < 50:
        print("[warn] .X is not raw counts and no 'counts' layer — pseudobulk needs raw counts")
    ct = first_present(CT_COLS, adata.obs) or sys.exit("[err] no cell_type column")
    sidc = first_present(SAMPLE_COLS, adata.obs) or sys.exit("[err] no sample_id column")
    adata.obs["_sample"] = adata.obs[sidc].astype(str)

    meta = pd.read_csv(args.meta, sep="\t")
    meta["sex"] = meta["sex"].map(norm_sex)
    msex = dict(zip(meta["sample_id"], meta["sex"]))
    mcond = dict(zip(meta["sample_id"], meta["condition_harmonized"]))
    mstudy = dict(zip(meta["sample_id"], meta["study_accession"]))

    rows = []
    for ctype in args.cell_types:
        cmask = (adata.obs[ct].astype(str) == ctype).values
        if cmask.sum() == 0:
            print(f"[skip] {ctype}: not present"); continue
        for cname, ref, tests in COMPARISONS:
            # sample sets for this contrast
            samples = [s for s in adata.obs.loc[cmask, "_sample"].unique()
                       if s not in EXCLUDE_SAMPLES and mcond.get(s) in ([ref] + tests)]
            clin = pd.DataFrame({"sample": samples,
                                 "group": ["reference" if mcond[s] == ref else "test" for s in samples],
                                 "sex": [msex.get(s, "unk") for s in samples],
                                 "study": [mstudy.get(s, "NA") for s in samples]}).set_index("sample")
            nref = (clin["group"] == "reference").sum(); ntest = (clin["group"] == "test").sum()
            sexM = (clin["sex"] == "M").sum(); sexF = (clin["sex"] == "F").sum()
            row = dict(cell_type=ctype, comparison=cname, n_ref=int(nref), n_test=int(ntest),
                       M=int(sexM), F=int(sexF), unk=int((clin["sex"] == "unk").sum()))
            if nref < MIN_SAMPLES_PER_GROUP or ntest < MIN_SAMPLES_PER_GROUP:
                row["status"] = "underpowered"; rows.append(row); continue
            cdf = pseudobulk(adata, clin.index.tolist(), cmask)
            if cdf is None or cdf.shape[1] < 2 * MIN_SAMPLES_PER_GROUP:
                row["status"] = "too_few_after_cellfilter"; rows.append(row); continue
            clin = clin.loc[cdf.columns]
            try:
                naive = deseq(cdf, clin.copy(), ["group"])
                row["n_DEG_naive"] = len(naive)
            except Exception as e:
                row["status"] = f"naive_failed:{e}"; rows.append(row); continue
            if not estimable(clin):
                row["status"] = "sex_term_not_estimable"; rows.append(row); continue
            try:
                adj = deseq(cdf, clin.copy(), ["sex", "group"])
                ov = len(set(naive.index) & set(adj.index))
                row.update(n_DEG_sexadj=len(adj), n_overlap=ov,
                           frac_naive_retained=round(ov / max(len(naive), 1), 3),
                           status="sex_adjusted_OK")
            except Exception as e:
                row["status"] = f"sexadj_failed:{e}"
            rows.append(row)
            print(f"[{ctype} | {cname}] {row.get('status')}: "
                  f"naive={row.get('n_DEG_naive')} sexadj={row.get('n_DEG_sexadj')} "
                  f"overlap={row.get('n_overlap')}")

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(args.out, "sex_adjustment_summary.csv"), index=False)
    print("\n[summary]\n", df.to_string(index=False))
    print("\n[done] outputs in", args.out)


if __name__ == "__main__":
    main()
