#!/usr/bin/env python3
"""Assemble tiered v4 bridge exports into per-compartment h5ad files.

For each compartment (NP, AF, CEP, all_cells) the tiered v4 integration
produced bridge files (counts.mtx.gz, embedding_pca.csv.gz,
embedding_umap.csv.gz, metadata.csv.gz, genes.csv) under
`data/integrated/{compartment}_experiment/tiered_v4/{mesenchymal,
non_mesenchymal}/`. This script merges the mes + non-mes tiers into a
single AnnData per compartment, padding the per-tier PCA/UMAP
embeddings with NaN across tiers, and writes
`data/integrated/tiered_v4/{compartment}.h5ad` for downstream modules
to consume.

The output schema matches what `scripts/06_clustering.py` expects:
- `obs['cell_class']` ∈ {mesenchymal, unknown, non_mesenchymal}
- `obsm['X_integrated']` (n × 50) — per-tier PCA, NaN-padded across
  tiers (clustering happens per-tier so the cross-tier coordinate
  system never matters)
- `obsm['X_umap']` (n × 2) — per-tier UMAP, NaN-padded
- `X` = raw counts (CSR sparse, gene union across tiers; missing genes
  are filled with 0 by anndata's outer concat)

The v5 h5ads at `data/integrated/{compartment}.h5ad` are NOT
overwritten; both pipelines coexist on disk.

Usage:
    python3 scripts/05m_assemble_tiered_v4.py
    python3 scripts/05m_assemble_tiered_v4.py --compartment NP
"""

import argparse
import gc
import sys
import warnings
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from scipy.io import mmread
from scipy.sparse import csr_matrix

anndata.settings.allow_write_nullable_strings = True
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)

BASE = Path(__file__).resolve().parent.parent
INTEG = BASE / "data" / "integrated"
OUT_DIR = INTEG / "tiered_v4"
OUT_DIR.mkdir(parents=True, exist_ok=True)

COMPARTMENTS = {
    "NP":        ("np_experiment",        ["mesenchymal", "non_mesenchymal"]),
    "AF":        ("af_experiment",        ["mesenchymal"]),  # non_mes not run (too few cells/study)
    "CEP":       ("cep_experiment",       ["mesenchymal", "non_mesenchymal"]),
    "all_cells": ("all_cells_experiment", ["mesenchymal", "non_mesenchymal"]),
}

EMBEDDING_DIM = 50
UMAP_DIM = 2


def _load_tier(bridge_dir):
    """Load a single tier from bridge files into an AnnData (cells x genes)."""
    bridge_dir = Path(bridge_dir)
    print(f"    Loading bridge: {bridge_dir.name}")

    metadata = pd.read_csv(bridge_dir / "metadata.csv.gz")
    barcodes = pd.read_csv(bridge_dir / "barcodes.csv")['barcode'].tolist()
    genes_df = pd.read_csv(bridge_dir / "genes.csv")
    gene_names = genes_df['gene'].tolist()

    if len(barcodes) != len(metadata):
        raise ValueError(f"barcodes ({len(barcodes)}) != metadata rows ({len(metadata)})")

    print(f"      cells={len(barcodes):,}  genes={len(gene_names):,}")
    print(f"      reading counts.mtx.gz ...")
    mtx = mmread(str(bridge_dir / "counts.mtx.gz"))
    # Bridge writes genes x cells from Seurat; transpose to cells x genes
    if mtx.shape[0] == len(gene_names) and mtx.shape[1] == len(barcodes):
        mtx = mtx.T
    elif mtx.shape[0] == len(barcodes) and mtx.shape[1] == len(gene_names):
        pass
    else:
        raise ValueError(f"counts shape {mtx.shape} matches neither "
                         f"({len(gene_names)},{len(barcodes)}) nor "
                         f"({len(barcodes)},{len(gene_names)})")
    counts = csr_matrix(mtx)

    pca = pd.read_csv(bridge_dir / "embedding_pca.csv.gz", index_col=0).values.astype(np.float32)
    umap = pd.read_csv(bridge_dir / "embedding_umap.csv.gz", index_col=0).values.astype(np.float32)

    if pca.shape[0] != len(barcodes):
        raise ValueError(f"PCA rows {pca.shape[0]} != cells {len(barcodes)}")
    if umap.shape[0] != len(barcodes):
        raise ValueError(f"UMAP rows {umap.shape[0]} != cells {len(barcodes)}")

    metadata.index = pd.Index(barcodes, name='cell_id')
    adata = anndata.AnnData(X=counts, obs=metadata,
                            var=pd.DataFrame(index=pd.Index(gene_names, name='gene')))
    adata.obsm['X_pca_tier'] = pca
    adata.obsm['X_umap_tier'] = umap
    print(f"      tier AnnData: {adata.shape}")
    return adata


def _concat_tiers(tier_adatas, tier_names):
    """Concatenate tiers cell-wise with a 'tier' obs column.

    Uses anndata.concat with join='outer' so the gene set is the union.
    Cells from one tier get 0 for genes only present in the other.
    """
    for ad, name in zip(tier_adatas, tier_names):
        ad.obs['tier'] = name
    print(f"    concat tiers: {[ad.n_obs for ad in tier_adatas]} cells per tier")
    combined = anndata.concat(tier_adatas, join='outer', axis=0,
                              merge='unique', label=None,
                              fill_value=0)
    # Make obs_names unique to avoid downstream collisions
    combined.obs_names_make_unique()
    print(f"    combined: {combined.shape}")
    return combined


def _build_padded_embeddings(combined, tier_names, tier_lengths):
    """Build NaN-padded X_integrated and X_umap from per-tier obsm.

    Each tier's PCA fills the rows that came from that tier; the
    other rows are NaN. Same for UMAP. Per-tier obsm slots
    (X_pca_tier, X_umap_tier) are concatenated by anndata as full
    matrices with the in-tier rows filled and the out-of-tier rows
    set to 0 — but we want NaN so the 06_clustering tier subset
    correctly drops them.
    """
    n = combined.n_obs
    x_int = np.full((n, EMBEDDING_DIM), np.nan, dtype=np.float32)
    x_umap = np.full((n, UMAP_DIM), np.nan, dtype=np.float32)

    cursor = 0
    for tier_name, length in zip(tier_names, tier_lengths):
        # The tiers were concatenated in order, so this tier occupies
        # rows [cursor, cursor + length).
        # The post-concat obsm slot has the tier's own values in those
        # rows and zeros for the other tiers — copy those back into
        # the NaN-padded matrices for the in-tier rows only.
        slice_idx = slice(cursor, cursor + length)
        x_int[slice_idx] = combined.obsm['X_pca_tier'][slice_idx]
        x_umap[slice_idx] = combined.obsm['X_umap_tier'][slice_idx]
        cursor += length

    combined.obsm['X_integrated'] = x_int
    combined.obsm['X_umap'] = x_umap

    # Drop per-tier scratch slots
    del combined.obsm['X_pca_tier']
    del combined.obsm['X_umap_tier']
    return combined


def assemble_compartment(compartment):
    cfg_dir, tier_names = COMPARTMENTS[compartment]
    tiered_dir = INTEG / cfg_dir / "tiered_v4"

    print(f"\n{'='*60}\n{compartment}\n{'='*60}")
    print(f"  Bridge dir: {tiered_dir}")
    print(f"  Tiers: {tier_names}")

    tier_adatas = []
    tier_lengths = []
    for tier in tier_names:
        bd = tiered_dir / tier
        if not bd.exists():
            print(f"  SKIP tier: {tier} (no bridge dir)")
            continue
        ad = _load_tier(bd)
        tier_adatas.append(ad)
        tier_lengths.append(ad.n_obs)

    if not tier_adatas:
        print(f"  ERROR: no tiers loaded for {compartment}")
        return None

    if len(tier_adatas) == 1:
        combined = tier_adatas[0]
        combined.obs['tier'] = tier_names[0]
        combined.obsm['X_integrated'] = combined.obsm.pop('X_pca_tier')
        combined.obsm['X_umap'] = combined.obsm.pop('X_umap_tier')
    else:
        combined = _concat_tiers(tier_adatas, tier_names[:len(tier_adatas)])
        combined = _build_padded_embeddings(combined, tier_names[:len(tier_adatas)],
                                            tier_lengths)

    # Sanity: cell_class column already in metadata (from bridge), but verify
    if 'cell_class' not in combined.obs.columns:
        raise ValueError(f"{compartment}: cell_class missing in metadata")

    print(f"  Final shape: {combined.shape}")
    print(f"  cell_class: {combined.obs['cell_class'].value_counts().to_dict()}")
    print(f"  tier: {combined.obs['tier'].value_counts().to_dict()}")
    print(f"  obsm keys: {list(combined.obsm.keys())}")
    print(f"  X dtype/sparse: {type(combined.X).__name__}, nnz={combined.X.nnz:,}")

    out_path = OUT_DIR / f"{compartment}.h5ad"
    tmp_path = Path(str(out_path) + ".tmp")
    print(f"  Writing {out_path} ...")
    combined.write_h5ad(tmp_path, compression='gzip')
    tmp_path.rename(out_path)
    size_gb = out_path.stat().st_size / 1e9
    print(f"  Wrote {out_path} ({size_gb:.2f} GB)")

    del combined, tier_adatas
    gc.collect()
    return out_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--compartment', type=str, default=None,
                        choices=list(COMPARTMENTS.keys()))
    args = parser.parse_args()

    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Output dir: {OUT_DIR}")

    targets = [args.compartment] if args.compartment else list(COMPARTMENTS.keys())
    written = []
    for c in targets:
        path = assemble_compartment(c)
        if path:
            written.append(path)

    print(f"\n{'='*60}\nSummary\n{'='*60}")
    for p in written:
        print(f"  {p} ({p.stat().st_size / 1e9:.2f} GB)")
    print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    main()
