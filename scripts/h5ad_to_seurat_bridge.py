#!/usr/bin/env python3
"""Convert h5ad files to a format Seurat v5 can load.

Saves sparse counts matrix (.mtx), barcodes, features, and metadata as
separate files that Seurat's Read10X + AddMetaData can consume.

Usage:
    python3 scripts/h5ad_to_seurat_bridge.py data/processed/GSE160756.h5ad /tmp/bridge/GSE160756
    python3 scripts/h5ad_to_seurat_bridge.py data/processed/GSE160756.h5ad /tmp/bridge/GSE160756 --compartment NP
"""

import sys
import os
import warnings
from pathlib import Path
from scipy.io import mmwrite
from scipy.sparse import csc_matrix
import pandas as pd
import numpy as np
import anndata

anndata.settings.allow_write_nullable_strings = True
warnings.filterwarnings('ignore')


def convert_h5ad_to_bridge(h5ad_path, out_dir, compartment_filter=None,
                           exclude_samples=None):
    """Convert an h5ad file to mtx + metadata for R loading."""
    import scanpy as sc

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(h5ad_path)

    # Filter by compartment if specified
    if compartment_filter and 'compartment' in adata.obs.columns:
        if compartment_filter == 'ALL':
            pass  # Keep all cells
        else:
            mask = adata.obs['compartment'] == compartment_filter
            if mask.sum() == 0:
                print(f"  WARNING: No cells with compartment={compartment_filter}")
                return False
            adata = adata[mask].copy()

    # Exclude samples
    if exclude_samples and 'sample_id' in adata.obs.columns:
        for s in exclude_samples:
            mask = adata.obs['sample_id'] != s
            adata = adata[mask].copy()

    if adata.shape[0] == 0:
        print(f"  WARNING: No cells remaining after filtering")
        return False

    # Get raw counts
    if 'counts' in adata.layers:
        X = adata.layers['counts']
    elif adata.raw is not None:
        X = adata.raw.X
    else:
        X = adata.X

    # Ensure sparse CSC for mmwrite
    X = csc_matrix(X)

    # Save matrix (genes x cells for Read10X compatibility)
    mmwrite(str(out_dir / "matrix.mtx"), X.T)

    # Save barcodes
    barcodes = adata.obs_names.tolist()
    pd.DataFrame(barcodes, columns=['barcode']).to_csv(
        out_dir / "barcodes.tsv", sep='\t', index=False, header=False)

    # Save features/genes
    genes = adata.var_names.tolist()
    gene_ids = adata.var['gene_ids'].tolist() if 'gene_ids' in adata.var.columns else genes
    pd.DataFrame({'gene_id': gene_ids, 'gene_name': genes}).to_csv(
        out_dir / "features.tsv", sep='\t', index=False, header=False)

    # Save metadata
    meta_cols = ['sample_id', 'study', 'compartment', 'condition_harmonized',
                 'cell_class', 'coarse_label', 'pct_counts_mt', 'n_genes_by_counts',
                 'total_counts', 'donor_id', 'age', 'sex']
    available_cols = [c for c in meta_cols if c in adata.obs.columns]
    meta = adata.obs[available_cols].copy()
    meta.index.name = 'barcode'
    meta.to_csv(out_dir / "metadata.csv")

    # Save embeddings if present
    for key in adata.obsm.keys():
        if key.startswith('X_'):
            emb = adata.obsm[key]
            if not isinstance(emb, np.ndarray):
                emb = np.array(emb)
            pd.DataFrame(emb, index=adata.obs_names).to_csv(
                out_dir / f"{key}.csv")

    print(f"  Converted: {h5ad_path} -> {out_dir} ({adata.shape[0]} cells, {adata.shape[1]} genes)")
    return True


if __name__ == "__main__":
    h5ad_path = sys.argv[1]
    out_dir = sys.argv[2]
    compartment = sys.argv[3] if len(sys.argv) > 3 else None
    exclude = sys.argv[4].split(',') if len(sys.argv) > 4 else None
    convert_h5ad_to_bridge(h5ad_path, out_dir, compartment, exclude)
