"""Re-render docs/manuscript_figures/fig01_umap_all_cells.png and fig02_umap_compartments.png
from the post-harmonization tiered_v4 h5ads at full cell count (no downsampling).

Matches the existing visual layout: fig01 = single all_cells panel with side legend,
fig02 = three horizontal compartment panels (NP, AF, CEP) with in-plot labels.
"""
from __future__ import annotations
import os, anndata as ad, numpy as np, pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import scanpy as sc

mpl.rcParams["pdf.fonttype"] = 42
sc.settings.set_figure_params(dpi=150, facecolor="white", frameon=False)

BASE = "data/integrated/tiered_v4"
OUT = "docs/manuscript_figures"
os.makedirs(OUT, exist_ok=True)


def lite_load(path: str) -> ad.AnnData:
    """Load only what's needed for sc.pl.umap: obs + obsm['X_umap'] + .uns colors."""
    big = sc.read_h5ad(path, backed="r")
    obs = big.obs.to_memory() if hasattr(big.obs, "to_memory") else big.obs.copy()
    umap = np.asarray(big.obsm["X_umap"][:])
    uns = dict(big.uns) if hasattr(big, "uns") else {}
    small = ad.AnnData(
        X=np.zeros((obs.shape[0], 1), dtype=np.float32),
        obs=obs,
        obsm={"X_umap": umap},
        uns=uns,
    )
    return small


def fig01_all_cells():
    a = lite_load(f"{BASE}/all_cells.h5ad")
    n = a.shape[0]
    n_ct = a.obs["cell_type"].astype("category").cat.categories.size
    print(f"[fig01] all_cells loaded: {n:,} cells, {n_ct} cell types")
    fig, ax = plt.subplots(figsize=(12, 7.5))
    sc.pl.umap(
        a, color="cell_type",
        size=2, alpha=0.55,
        legend_loc="right margin",
        legend_fontsize=8, legend_fontoutline=0,
        title=f"Integrated atlas: {n:,} cells coloured by cell_type",
        ax=ax, show=False, frameon=False,
    )
    plt.tight_layout()
    out = f"{OUT}/fig01_umap_all_cells.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig01] wrote {out}")


def fig02_compartments():
    panels = [("NP", f"{BASE}/NP.h5ad"),
              ("AF", f"{BASE}/AF.h5ad"),
              ("CEP", f"{BASE}/CEP.h5ad")]
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    for ax, (name, path) in zip(axes, panels):
        a = lite_load(path)
        n = a.shape[0]
        n_ct = a.obs["cell_type"].astype("category").cat.categories.size
        print(f"[fig02:{name}] {n:,} cells, {n_ct} cell types")
        sc.pl.umap(
            a, color="cell_type",
            size=3, alpha=0.55,
            legend_loc="on data",
            legend_fontsize=7, legend_fontoutline=2,
            title=f"{name} ({n:,} cells, {n_ct} cell types)",
            ax=ax, show=False, frameon=False,
        )
    plt.tight_layout()
    out = f"{OUT}/fig02_umap_compartments.png"
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"[fig02] wrote {out}")


if __name__ == "__main__":
    fig02_compartments()
    fig01_all_cells()
    print("[done]")
