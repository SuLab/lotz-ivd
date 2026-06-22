"""Render docs/manuscript_figures/fig02b_umap_substates.png — mesenchymal-tier UMAPs
per compartment, coloured by Stage 3 cell_subtype. Matches the §6g view in
notebooks/07_annotation.ipynb (Module 7 tiered_v4).

Three horizontal panels: NP (8 sub-states), AF (5), CEP (6). Full cell counts
from the post-harmonization h5ads — no downsampling. Same X_umap embedding
used in fig02 (cell_type view), so sub-states sit in the same spatial layout.
"""
from __future__ import annotations
import os, anndata as ad, numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scanpy as sc

mpl.rcParams["pdf.fonttype"] = 42
sc.settings.set_figure_params(dpi=150, facecolor="white", frameon=False)

BASE = "data/integrated/tiered_v4"
OUT = "docs/manuscript_figures"
os.makedirs(OUT, exist_ok=True)


def lite_mes_load(path: str) -> ad.AnnData:
    """Load only mesenchymal-tier cells: obs + obsm['X_umap']."""
    big = sc.read_h5ad(path, backed="r")
    obs = big.obs.to_memory() if hasattr(big.obs, "to_memory") else big.obs.copy()
    umap = np.asarray(big.obsm["X_umap"][:])
    if "tier" in obs.columns:
        keep = obs["tier"].astype(str).eq("mesenchymal").values
    elif "cell_class" in obs.columns:
        keep = obs["cell_class"].astype(str).eq("mesenchymal").values
    else:
        keep = np.ones(obs.shape[0], dtype=bool)
    obs = obs.loc[keep].copy()
    umap = umap[keep]
    return ad.AnnData(
        X=np.zeros((obs.shape[0], 1), dtype=np.float32),
        obs=obs,
        obsm={"X_umap": umap},
    )


def fig02b():
    panels = [("NP", f"{BASE}/NP.h5ad"),
              ("AF", f"{BASE}/AF.h5ad"),
              ("CEP", f"{BASE}/CEP.h5ad")]
    # Stack the three compartments vertically (3 rows x 1 col) so each UMAP
    # spans the full figure/page width — much larger per-panel than the old
    # 1x3 wide row, which rendered tiny when embedded at page width.
    fig, axes = plt.subplots(3, 1, figsize=(8, 22))
    for ax, (name, path) in zip(axes, panels):
        a = lite_mes_load(path)
        n = a.shape[0]
        # Drop "unassigned" residual from category list for display
        cs = a.obs["cell_subtype"].astype(str)
        n_unassigned = int((cs == "unassigned").sum())
        n_sub = cs[cs != "unassigned"].nunique()
        print(f"[fig02b:{name}] mesenchymal cells={n:,}  unique cell_subtypes (excl unassigned)={n_sub}  unassigned={n_unassigned:,}")
        sc.pl.umap(
            a, color="cell_subtype",
            size=4, alpha=0.55,
            legend_loc="on data",
            legend_fontsize=9, legend_fontoutline=2,
            title=f"{name} mesenchymal ({n:,} cells, {n_sub} sub-states)",
            ax=ax, show=False, frameon=False,
        )
    plt.tight_layout()
    out = f"{OUT}/fig02b_umap_substates.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig02b] wrote {out}")


if __name__ == "__main__":
    fig02b()
    print("[done]")
