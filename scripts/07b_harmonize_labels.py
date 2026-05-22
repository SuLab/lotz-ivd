"""Module 07b — label harmonization + contamination flagging (post-Module-07).

One-shot script that updates the tiered_v4 h5ads in place:

1. **Label renames** (compartment-prefix consistency for downstream
   pseudobulk grouping):

   - `EP_hyaline`                  → `CEP_hyaline`
   - `EP_ossification`             → `CEP_ossification`
   - `Fibroblast_like` (CEP only)  → `CEP_outer`
   - `Fibrochondrocyte_chondroid`  → `NP_fibrochondrocyte_chondroid`
   - `Fibrochondrocyte_fibroid`    → `CEP_fibrochondrocyte_fibroid`

   `Fibroblast_like` / `Chondrocyte_like` / `Fibrochondrocyte_like` /
   `Macrophage` (generic) labels coming from `compartment = "IVD_mixed"`
   (GSE189916) are intentionally left untouched — that dataset doesn't
   separate NP/AF/CEP at sampling, so the generic labels are honest.

   Renames are also applied to `coarse_cell_type` (where the same string
   may appear) and to `cell_subtype` (where labels are `<cell_type>_<state>`).

2. **Contamination flag columns** (used by Module 08+ to filter or annotate):

   - `is_contamination` (bool): True iff `cell_type == "Erythrocyte"` OR
     `cell_subtype.endswith("_endothelial_admixed")`
   - `contamination_type` (str): one of {"RBC", "endothelial_admixed", "clean"}.

3. **Regenerate** `results/integration/tiered_v4/cell_type_definitions.tsv`
   with the harmonized labels.

Idempotent — safe to re-run.
"""

from __future__ import annotations

import sys
from pathlib import Path

import anndata
import pandas as pd

anndata.settings.allow_write_nullable_strings = True

BASE = Path(__file__).resolve().parent.parent
INT_DIR = BASE / "data" / "integrated" / "tiered_v4"
RESULTS_DIR = BASE / "results" / "integration" / "tiered_v4"

# ── Cell_type renames (apply to cell_type AND coarse_cell_type columns) ────
CELL_TYPE_RENAMES = {
    "EP_hyaline":                 "CEP_hyaline",
    "EP_ossification":             "CEP_ossification",
    "Fibrochondrocyte_chondroid":  "NP_fibrochondrocyte_chondroid",
    "Fibrochondrocyte_fibroid":    "CEP_fibrochondrocyte_fibroid",
}

# `Fibroblast_like` rename is COMPARTMENT-CONDITIONAL: only cells in
# CEP-compartment cells become `CEP_outer`. Cells in `IVD_mixed`
# (GSE189916) keep the generic `Fibroblast_like` label because that dataset
# doesn't separate compartments.
CONDITIONAL_RENAMES = [
    {
        "from_label": "Fibroblast_like",
        "to_label":   "CEP_outer",
        "compartment": "CEP",
    },
]

OBJECTS = ["NP", "AF", "CEP", "all_cells"]


def _split_subtype(label: str) -> tuple[str, str]:
    """Split `cell_subtype` value into (base cell_type, sub-state).

    Sub-states are one of {proliferating, inflammatory, stressed, matrix_active,
    migratory, homeostatic, endothelial_admixed}. The cell_type prefix may
    itself contain underscores (e.g. `NP_mature_chondrocyte_matrix_active`),
    so we split on the longest known sub-state suffix.
    """
    known_suffixes = (
        "_endothelial_admixed",
        "_matrix_active",
        "_proliferating",
        "_inflammatory",
        "_homeostatic",
        "_migratory",
        "_stressed",
    )
    for suf in known_suffixes:
        if label.endswith(suf):
            return label[: -len(suf)], suf[1:]  # strip leading underscore
    return label, ""


def _apply_renames(adata: anndata.AnnData, object_name: str) -> dict:
    """Apply cell_type / coarse_cell_type / cell_subtype renames in place.

    Returns a stats dict for the run summary.
    """
    stats = {"object": object_name, "renamed_cell_type": 0,
             "renamed_subtype": 0, "renamed_coarse": 0,
             "renamed_cep_outer": 0}

    # 1. Unconditional renames on cell_type and coarse_cell_type
    for col in ("cell_type", "coarse_cell_type"):
        if col not in adata.obs.columns:
            continue
        before = adata.obs[col].astype(str)
        after = before.replace(CELL_TYPE_RENAMES)
        n_changed = int((before != after).sum())
        adata.obs[col] = after.values
        if col == "cell_type":
            stats["renamed_cell_type"] = n_changed
        else:
            stats["renamed_coarse"] = n_changed

    # 2. Conditional rename: CEP-compartment Fibroblast_like → CEP_outer
    if "cell_type" in adata.obs.columns and "compartment" in adata.obs.columns:
        for rule in CONDITIONAL_RENAMES:
            mask = (adata.obs["cell_type"].astype(str) == rule["from_label"]) & \
                   (adata.obs["compartment"].astype(str) == rule["compartment"])
            n = int(mask.sum())
            if n > 0:
                ct = adata.obs["cell_type"].astype(str).copy()
                ct.loc[mask] = rule["to_label"]
                adata.obs["cell_type"] = ct.values
                stats["renamed_cep_outer"] += n

    # 3. cell_subtype renames — split into base + sub-state, rename base, recombine
    if "cell_subtype" in adata.obs.columns:
        st_before = adata.obs["cell_subtype"].astype(str)
        new_st = []
        for label in st_before:
            base, state = _split_subtype(label)
            # Apply unconditional cell_type renames to the base
            base_new = CELL_TYPE_RENAMES.get(base, base)
            new_label = f"{base_new}_{state}" if state else base_new
            new_st.append(new_label)
        st_after = pd.Series(new_st, index=adata.obs.index)
        stats["renamed_subtype"] = int((st_before != st_after).sum())
        adata.obs["cell_subtype"] = st_after.values

        # Conditional CEP_outer rename on cell_subtype too
        if "compartment" in adata.obs.columns:
            for rule in CONDITIONAL_RENAMES:
                # Find subtypes that start with the from_label
                pattern = rule["from_label"] + "_"
                mask = (adata.obs["cell_subtype"].astype(str).str.startswith(pattern)) & \
                       (adata.obs["compartment"].astype(str) == rule["compartment"])
                n = int(mask.sum())
                if n > 0:
                    st = adata.obs["cell_subtype"].astype(str).copy()
                    st.loc[mask] = st.loc[mask].str.replace(
                        pattern, rule["to_label"] + "_", n=1, regex=False)
                    adata.obs["cell_subtype"] = st.values
                    stats["renamed_subtype"] += n

    return stats


def _add_contamination_flags(adata: anndata.AnnData) -> dict:
    """Add `is_contamination` (bool) and `contamination_type` (str) columns."""
    stats = {"n_rbc": 0, "n_endo_admixed": 0, "n_clean": 0}

    ct = adata.obs["cell_type"].astype(str) if "cell_type" in adata.obs.columns \
        else pd.Series("", index=adata.obs.index)
    st = adata.obs["cell_subtype"].astype(str) if "cell_subtype" in adata.obs.columns \
        else pd.Series("", index=adata.obs.index)

    rbc_mask = (ct == "Erythrocyte")
    endo_mask = st.str.endswith("_endothelial_admixed") & (~rbc_mask)

    contamination_type = pd.Series("clean", index=adata.obs.index, dtype=object)
    contamination_type.loc[rbc_mask] = "RBC"
    contamination_type.loc[endo_mask] = "endothelial_admixed"

    adata.obs["is_contamination"] = (rbc_mask | endo_mask).values
    adata.obs["contamination_type"] = contamination_type.values

    stats["n_rbc"] = int(rbc_mask.sum())
    stats["n_endo_admixed"] = int(endo_mask.sum())
    stats["n_clean"] = int((~(rbc_mask | endo_mask)).sum())
    return stats


def _regenerate_cell_type_definitions():
    """Rebuild cell_type_definitions.tsv from the harmonized h5ads."""
    rows = []
    for obj in OBJECTS:
        path = INT_DIR / f"{obj}.h5ad"
        if not path.exists():
            continue
        a = anndata.read_h5ad(path, backed="r")
        if "cell_type" not in a.obs.columns:
            a.file.close()
            continue
        obs = a.obs[["cell_type", "coarse_cell_type", "cell_type_confidence",
                     "cell_subtype"]].copy() if "cell_subtype" in a.obs.columns \
              else a.obs[["cell_type", "coarse_cell_type", "cell_type_confidence"]].copy()
        a.file.close()

        for ct in sorted(obs["cell_type"].astype(str).unique()):
            ct_mask = obs["cell_type"].astype(str) == ct
            n = int(ct_mask.sum())
            coarse = obs.loc[ct_mask, "coarse_cell_type"].mode().iloc[0] \
                if "coarse_cell_type" in obs.columns else ""
            conf_dist = obs.loc[ct_mask, "cell_type_confidence"].value_counts().to_dict() \
                if "cell_type_confidence" in obs.columns else {}
            conf_str = "; ".join(f"{k}={v}" for k, v in conf_dist.items())
            subtypes = sorted(obs.loc[ct_mask, "cell_subtype"].astype(str).unique()) \
                if "cell_subtype" in obs.columns else []
            rows.append({
                "object": obj,
                "cell_type": ct,
                "coarse_cell_type": coarse,
                "n_cells": n,
                "n_subtypes": len(subtypes),
                "subtypes": ", ".join(subtypes),
                "confidence_distribution": conf_str,
            })

    df = pd.DataFrame(rows)
    out_path = RESULTS_DIR / "cell_type_definitions.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    print(f"  Regenerated {out_path} ({len(df)} rows)")


def main():
    print("=" * 60)
    print("Module 07b — Label Harmonization + Contamination Flagging")
    print("=" * 60)

    all_stats = []
    for obj in OBJECTS:
        path = INT_DIR / f"{obj}.h5ad"
        if not path.exists():
            print(f"\n  {obj}: input not found at {path}, skipping")
            continue
        print(f"\n  Loading {obj}.h5ad ...")
        adata = anndata.read_h5ad(path)
        print(f"    {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

        rename_stats = _apply_renames(adata, obj)
        contam_stats = _add_contamination_flags(adata)

        # Drop the obsolete categorical categories before write (Categorical
        # columns retain the old levels otherwise).
        for col in ("cell_type", "coarse_cell_type", "cell_subtype"):
            if col in adata.obs.columns and \
                    isinstance(adata.obs[col].dtype, pd.CategoricalDtype):
                adata.obs[col] = adata.obs[col].astype(object)

        print(f"    Renames: cell_type={rename_stats['renamed_cell_type']:,}, "
              f"CEP_outer={rename_stats['renamed_cep_outer']:,}, "
              f"coarse={rename_stats['renamed_coarse']:,}, "
              f"subtype={rename_stats['renamed_subtype']:,}")
        print(f"    Contamination: RBC={contam_stats['n_rbc']:,}, "
              f"endothelial_admixed={contam_stats['n_endo_admixed']:,}, "
              f"clean={contam_stats['n_clean']:,}")
        print(f"    Writing back to {path} ...")
        adata.write_h5ad(path, compression="gzip")
        all_stats.append({**rename_stats, **contam_stats})

    print("\n  Regenerating cell_type_definitions.tsv ...")
    _regenerate_cell_type_definitions()

    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(pd.DataFrame(all_stats).to_string(index=False))
    print("\nDone.")


if __name__ == "__main__":
    main()
