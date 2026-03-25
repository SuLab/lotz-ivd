#!/usr/bin/env python3
"""Assemble h5ad files from Seurat bridge exports.

Reads the flat files produced by seurat_to_h5ad_bridge.R and creates
proper AnnData h5ad files at data/integrated/{object}.h5ad.

Usage:
    python3 scripts/seurat_to_h5ad_assemble.py --workflow cca
    python3 scripts/seurat_to_h5ad_assemble.py --workflow cca --object NP
"""

import sys
import gc
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import anndata
import scanpy as sc
from scipy.io import mmread
from scipy.sparse import csr_matrix

anndata.settings.allow_write_nullable_strings = True
warnings.filterwarnings('ignore')

BASE = Path(__file__).resolve().parent.parent
INT_DIR = BASE / "data" / "integrated"

OBJECTS = ["NP", "AF", "CEP", "all_cells"]


def assemble_h5ad(workflow, obj_name):
    """Assemble one h5ad from bridge export files."""
    export_dir = INT_DIR / workflow / "bridge_export" / obj_name
    if not export_dir.exists():
        print(f"  SKIP: {export_dir} not found")
        return False

    print(f"\n{'='*60}")
    print(f"Assembling {workflow}/{obj_name}")
    print(f"{'='*60}")

    # 1. Load counts matrix (genes x cells in R; transpose to cells x genes)
    mtx_path = export_dir / "counts.mtx.gz"
    print(f"  Loading counts from {mtx_path.name}...")
    counts = mmread(str(mtx_path))
    counts = csr_matrix(counts.T)  # Transpose: cells x genes
    print(f"  Counts: {counts.shape[0]} cells x {counts.shape[1]} genes")

    # 2. Load gene names
    genes = pd.read_csv(export_dir / "genes.csv")['gene'].values
    print(f"  Genes: {len(genes)}")

    # 3. Load barcodes
    barcodes = pd.read_csv(export_dir / "barcodes.csv")['barcode'].values
    print(f"  Barcodes: {len(barcodes)}")

    # 4. Load metadata
    meta = pd.read_csv(export_dir / "metadata.csv.gz")
    meta.index = barcodes
    print(f"  Metadata: {meta.shape[1]} columns")

    # 5. Create AnnData
    adata = anndata.AnnData(
        X=counts,
        obs=meta,
        var=pd.DataFrame(index=genes)
    )

    # Store raw counts in layers and normalize
    adata.layers['counts'] = counts.copy()

    # Recompute log-normalization (equivalent to Seurat NormalizeData LogNormalize scale.factor=10000)
    print(f"  Computing log-normalization (matching Seurat NormalizeData)...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    print(f"  .X = log-normalized, .layers['counts'] = raw counts")

    # 7. Load embeddings
    import glob
    emb_files = sorted(export_dir.glob("embedding_*.csv.gz"))
    for emb_file in emb_files:
        red_name = emb_file.stem.replace("embedding_", "").replace(".csv", "")
        print(f"  Loading embedding: {red_name}...")
        emb_df = pd.read_csv(emb_file, index_col=0)
        emb = emb_df.values.astype(np.float32)
        print(f"    Shape: {emb.shape}")

        # Map Seurat reduction names to scanpy conventions
        if red_name in ("integrated.cca", "integrated.dr"):
            adata.obsm['X_integrated'] = emb
            print(f"    -> stored as X_integrated")
        elif red_name == "pca":
            adata.obsm['X_pca'] = emb
            print(f"    -> stored as X_pca")
        elif red_name == "umap":
            adata.obsm['X_umap'] = emb
            print(f"    -> stored as X_umap")
        else:
            key = f"X_{red_name.replace('.', '_')}"
            adata.obsm[key] = emb
            print(f"    -> stored as {key}")

    # 8. Compute neighbors and UMAP on X_integrated if no UMAP exists
    if 'X_integrated' in adata.obsm and 'X_umap' not in adata.obsm:
        print(f"  Computing neighbors + UMAP on X_integrated...")
        sc.pp.neighbors(adata, use_rep='X_integrated')
        sc.tl.umap(adata)
        print(f"  UMAP computed")

    # 9. Save
    out_path = INT_DIR / f"{obj_name}.h5ad"
    print(f"  Saving to {out_path}...")
    adata.write_h5ad(out_path)
    size_mb = out_path.stat().st_size / 1e6
    print(f"  Saved: {out_path} ({size_mb:.0f} MB)")
    print(f"  Final shape: {adata.shape}")
    print(f"  Embeddings: {list(adata.obsm.keys())}")
    print(f"  Layers: {list(adata.layers.keys())}")

    del adata
    gc.collect()
    return True


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--workflow', required=True, choices=['cca', 'stacas'])
    parser.add_argument('--object', default=None)
    args = parser.parse_args()

    objects = [args.object] if args.object else OBJECTS
    success = 0
    for obj_name in objects:
        if assemble_h5ad(args.workflow, obj_name):
            success += 1

    print(f"\n\nAssembly complete: {success}/{len(objects)} objects")


if __name__ == '__main__':
    main()
