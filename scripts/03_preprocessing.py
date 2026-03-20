#!/usr/bin/env python3
"""Module 03: Per-dataset preprocessing for IVD scRNA-seq atlas.

Applies consistent QC, normalization, and preliminary annotation to each dataset
independently. Produces clean AnnData objects, marker gene tables, and QC reports.

Usage:
    python3 scripts/03_preprocessing.py                  # Process all datasets
    python3 scripts/03_preprocessing.py GSE230809        # Process one dataset
    python3 scripts/03_preprocessing.py --validate-only  # Run validation only
"""

import sys
import os
import warnings
import tarfile
import tempfile
import gzip
import shutil
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
from scipy.io import mmread
import scrublet as scr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from jinja2 import Template

warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning, module='scanpy')

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent.parent
RAW_DIR = BASE / "data" / "raw"
PROC_DIR = BASE / "data" / "processed"
QC_DIR = BASE / "results" / "qc_reports"
ANN_DIR = BASE / "results" / "annotations"
META_FILE = BASE / "metadata" / "sample_metadata.tsv"

PROC_DIR.mkdir(parents=True, exist_ok=True)
QC_DIR.mkdir(parents=True, exist_ok=True)
ANN_DIR.mkdir(parents=True, exist_ok=True)

# ── Dataset ordering ─────────────────────────────────────────────────────────
ALL_ACCESSIONS = [
    "GSE160756", "GSE165722", "GSE189916", "GSE199866", "GSE205535",
    "CNP0002664", "GSE233666", "GSE244889", "GSE251686", "GSE255768",
    "GSE230809", "GSE242443",
]

# ── Sample-to-file mapping ───────────────────────────────────────────────────
# Maps sample_id → file identifier within the tar/directory
# Format: {accession: {sample_id: file_key}}

SAMPLE_FILE_MAP = {
    "GSE160756": {
        "GSE160756_hNP_1":  "GSM4878538_umi_hNP_1.loom.gz",
        "GSE160756_hNP_2":  "GSM4878539_umi_hNP_2.loom.gz",
        "GSE160756_hNP_3":  "GSM4878540_umi_hNP_3.loom.gz",
        "GSE160756_hCEP_1": "GSM4878541_umi_hCEP_1.loom.gz",
        "GSE160756_hCEP_2": "GSM4878542_umi_hCEP_2.loom.gz",
        "GSE160756_hAF_1":  "GSM4878543_umi_hAF_1.loom.gz",
        "GSE160756_hAF_2":  "GSM4878544_umi_hAF_2.loom.gz",
    },
    "GSE165722": {
        "GSE165722_S1": ("GSM5048708_Sample1.counts.tsv.gz", "GSM5048708_Sample1.cellname.txt.gz"),
        "GSE165722_S2": ("GSM5048709_Sample2.counts.tsv.gz", "GSM5048709_Sample2.cellname.txt.gz"),
        "GSE165722_S3": ("GSM5048710_Sample3.counts.tsv.gz", "GSM5048710_Sample3.cellname.txt.gz"),
        "GSE165722_S4": ("GSM5048711_Sample4.counts.tsv.gz", "GSM5048711_Sample4.cellname.txt.gz"),
        "GSE165722_S5": ("GSM5048712_Sample5.counts.tsv.gz", "GSM5048712_Sample5.cellname.txt.gz"),
        "GSE165722_S6": ("GSM5048713_Sample6.counts.tsv.gz", "GSM5048713_Sample6.cellname.txt.gz"),
        "GSE165722_S7": ("GSM5048714_Sample7.counts.tsv.gz", "GSM5048714_Sample7.cellname.txt.gz"),
        "GSE165722_S8": ("GSM5048715_Sample8.counts.tsv.gz", "GSM5048715_Sample8.cellname.txt.gz"),
    },
    "GSE189916": {
        "GSE189916_Neonatal_IVD_1": "GSM5709571_Sample_1",
        "GSE189916_Neonatal_IVD_2": "GSM5709572_Sample_2",
        "GSE189916_Neonatal_IVD_3": "GSM5709573_Sample_3",
        "GSE189916_Adult_IVD_1":    "GSM5709574_Sample_4",
        "GSE189916_Adult_IVD_2":    "GSM5709575_Sample_5_processed",
        "GSE189916_Adult_IVD_3":    "GSM5709576_Sample_6_processed",
    },
    "GSE199866": {
        "GSE199866_AFH": "GSM5989808_AFH_filtered_gene_bc_matrices_h5.h5",
        "GSE199866_AFD": "GSM5989809_AFD_filtered_gene_bc_matrices_h5.h5",
        "GSE199866_NPH": "GSM5989810_NPH_filtered_gene_bc_matrices_h5.h5",
        "GSE199866_NPD": "GSM5989811_NPD_filtered_gene_bc_matrices_h5.h5",
    },
    "GSE205535": {
        "GSE205535_NNP": "GSM6214392_NP-normal",
        "GSE205535_DNP": "GSM6214393_NP",
    },
    "CNP0002664": {
        "CNP0002664_Ctrl":  "ctrl_matrix.tsv.gz",
        "CNP0002664_NP2":   "NP2_matrix.tsv.gz",
        "CNP0002664_NP4":   "NP4_matrix.tsv.gz",
        "CNP0002664_NP8":   "NP8_matrix.tsv.gz",
        "CNP0002664_NP9":   "NP9_matrix.tsv.gz",
        "CNP0002664_NP10":  "NP10_matrix.tsv.gz",
    },
    "GSE233666": {
        "GSE233666_Patient_1": "GSM7432171_Patient_1",
        "GSE233666_Patient_2": "GSM7432172_Patient_2",
        "GSE233666_Patient_3": "GSM7432173_patient_3",  # lowercase in tar
        "GSE233666_Patient_4": "GSM7432174_Patient_4",
    },
    "GSE244889": {
        "GSE244889_Pa-17F":  "GSM7831813_Pa-17F",
        "GSE244889_Pb-55F":  "GSM7831814_Pb-55F",
        "GSE244889_Pab-24M": "GSM7831815_Pab-24M",
        "GSE244889_Pab-30M": "GSM7831816_Pab-30M",
        "GSE244889_Pc-41M":  "GSM7831817_Pc-41M",
        "GSE244889_Pd-62F":  "GSM7831818_Pd-62F",
        "GSE244889_Pd-59F":  "GSM7831819_Pd-59F",
    },
    "GSE251686": {
        "GSE251686_NP1": "GSM7986001_NP1_EmptyDrops_CR_matrix.tar.gz",
        "GSE251686_NP3": "GSM7986002_NP3_EmptyDrops_CR_matrix.tar.gz",
        "GSE251686_NP4": "GSM7986003_NP4_EmptyDrops_CR_matrix.tar.gz",
        "GSE251686_NP5": "GSM7986004_NP5_EmptyDrops_CR_matrix.tar.gz",
        "GSE251686_NP6": "GSM7986005_NP6_EmptyDrops_CR_matrix.tar.gz",
        "GSE251686_NP9": "GSM7986006_NP9_EmptyDrops_CR_matrix.tar.gz",
    },
    "GSE255768": {
        "GSE255768_S1": "GSM8079184_S1",
        "GSE255768_S2": "GSM8079185_S2",
    },
    "GSE230809": {
        "GSE230809_AF_SP21_015": "GSM7173748_AF_SP21_015",
        "GSE230809_AF_SP21_018": "GSM7173749_AF_SP21_018",
        "GSE230809_AF_SP22_001": "GSM7173750_AF_SP22_001",
        "GSE230809_NP_SP21_015": "GSM7173751_NP_SP21_015",
        "GSE230809_NP_SP21_018": "GSM7173752_NP_SP21_018",
        "GSE230809_NP_SP22_001": "GSM7173753_NP_SP22_001",
        "GSE230809_AF_SP20_002": "GSM7235331_AF_SP20_002",
        "GSE230809_AF_SP20_006": "GSM7235332_AF_SP20_006",
        "GSE230809_AF_SP21_007": "GSM7235333_AF_SP21_007",
        "GSE230809_AF_SP21_011": "GSM7235334_AF_SP21_011",
        "GSE230809_AF_SP21_013": "GSM7235335_AF_SP21_013",
        "GSE230809_AF_SP21_014": "GSM7235336_AF_SP21_014",
        "GSE230809_AF_SP21_016": "GSM7235337_AF_SP21_016",
        "GSE230809_AF_SP21_017": "GSM7235338_AF_SP21_017",
        "GSE230809_AF_SP22_002": "GSM7235339_AF_SP22_002",
        "GSE230809_AF_SP22_003": "GSM7235340_AF_SP22_003",
        "GSE230809_NP_SP21_007": "GSM7235341_NP_SP21_007",
        "GSE230809_NP_SP21_011": "GSM7235342_NP_SP21_011",
        "GSE230809_NP_SP21_013": "GSM7235343_NP_SP21_013",
        "GSE230809_NP_SP21_014": "GSM7235344_NP_SP21_014",
        "GSE230809_NP_SP21_016": "GSM7235345_NP_SP21_016",
        "GSE230809_NP_SP21_017": "GSM7235346_NP_SP21_017",
        "GSE230809_NP_SP22_002": "GSM7235347_NP_SP22_002",
        "GSE230809_NP_SP22_003": "GSM7235348_NP_SP22_003",
    },
    "GSE242443": {
        "GSE242443_Y8444_H1": "GSM7763584_Y8444_H1",
        "GSE242443_Y8445_D2": "GSM7763585_Y8445_D2",
    },
}

# ── QC thresholds ─────────────────────────────────────────────────────────────
QC_PARAMS = {
    "min_genes": 200,
    "max_genes": 6000,
    "max_pct_mt": 5.0,
    "min_counts": 1000,
    "max_counts": 25000,
    "doublet_score_threshold": 0.25,
    "min_cells_per_gene": 3,
}

# ── Cell type markers ─────────────────────────────────────────────────────────
CELL_TYPE_MARKERS = {
    "Chondrocyte_NP":       ["ACAN", "COL2A1", "SOX9", "KRT19"],
    "Fibroblast_AF":        ["COL1A1", "COL1A2", "THY1", "SCX"],
    "Endothelial":          ["PECAM1", "VWF", "CDH5", "FLT1"],
    "Macrophage":           ["CD68", "CD163", "CSF1R", "MRC1"],
    "T_cell":               ["CD3D", "CD3E", "CD4", "CD8A"],
    "B_cell":               ["CD79A", "MS4A1", "CD19"],
    "Mast_cell":            ["KIT", "TPSAB1", "CPA3"],
    "Pericyte_SMC":         ["ACTA2", "TAGLN", "MYH11", "RGS5"],
    "Notochordal":          ["T", "SHH", "NOG"],
    "Progenitor":           ["CD44", "PROM1", "NES"],
}

# ── Clustering resolutions ────────────────────────────────────────────────────
LEIDEN_RESOLUTIONS = [round(x * 0.1, 1) for x in range(1, 21)]  # 0.1 to 2.0 in 0.1 steps
WORKING_RESOLUTION = 0.5  # default; updated per-dataset by silhouette/modularity scoring


# ═══════════════════════════════════════════════════════════════════════════════
# DATASET LOADING
# ═══════════════════════════════════════════════════════════════════════════════

def _extract_tar(accession, extract_dir):
    """Extract the tar archive for a dataset to a temp directory."""
    tar_path = RAW_DIR / accession / f"{accession}_RAW.tar"
    if not tar_path.exists():
        raise FileNotFoundError(f"No tar file at {tar_path}")
    with tarfile.open(tar_path, 'r') as tf:
        tf.extractall(extract_dir)
    return extract_dir


def _load_10x_mtx(prefix, extract_dir):
    """Load a 10x-style MTX triplet (barcodes, features/genes, matrix).

    `prefix` is the common part before _barcodes/_features/_matrix.
    Handles inconsistent casing (e.g. Patient vs patient for features).
    """
    files = os.listdir(extract_dir)
    prefix_lower = prefix.lower()

    barcodes_file = features_file = matrix_file = None
    for f in files:
        fl = f.lower()
        if fl.startswith(prefix_lower.lower()) and "barcodes" in fl:
            barcodes_file = os.path.join(extract_dir, f)
        elif fl.startswith(prefix_lower.lower()) and ("features" in fl or "genes" in fl):
            features_file = os.path.join(extract_dir, f)
        elif fl.startswith(prefix_lower.lower()) and "matrix" in fl:
            matrix_file = os.path.join(extract_dir, f)

    if not all([barcodes_file, features_file, matrix_file]):
        raise FileNotFoundError(
            f"Missing 10x files for prefix '{prefix}' in {extract_dir}. "
            f"Found: barcodes={barcodes_file}, features={features_file}, matrix={matrix_file}"
        )

    # Read matrix
    mat = mmread(barcodes_file.replace("barcodes", "matrix").replace(
        barcodes_file, matrix_file) if False else matrix_file)
    mat = sparse.csr_matrix(mat.T)  # genes x cells → cells x genes

    # Read barcodes
    barcodes = pd.read_csv(barcodes_file, header=None, sep='\t')[0].values

    # Read features (may have 1, 2, or 3 columns)
    features = pd.read_csv(features_file, header=None, sep='\t')
    if features.shape[1] >= 2:
        gene_ids = features[0].values
        gene_names = features[1].values
    else:
        gene_ids = features[0].values
        gene_names = features[0].values

    adata = ad.AnnData(X=mat)
    adata.obs_names = pd.Index(barcodes).astype(str)
    adata.var_names = pd.Index(gene_names).astype(str)
    if not np.array_equal(gene_ids, gene_names):
        adata.var['gene_ids'] = gene_ids

    return adata


def _load_10x_h5(h5_path):
    """Load a 10x HDF5 file."""
    adata = sc.read_10x_h5(h5_path)
    adata.var_names_make_unique()
    return adata


def load_GSE160756(extract_dir):
    """Load loom.gz files for Gan 2021 using h5py (loompy validation too strict).

    Reads the dense HDF5 matrix in gene-chunks to avoid materializing the full
    ~3 GB float64 array in RAM per sample.
    """
    import h5py
    samples = {}
    for sample_id, filename in SAMPLE_FILE_MAP["GSE160756"].items():
        gz_path = os.path.join(extract_dir, filename)
        # Decompress .loom.gz to .loom
        loom_path = gz_path.replace('.gz', '')
        with gzip.open(gz_path, 'rb') as f_in, open(loom_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

        # Read loom file directly with h5py — chunked to limit memory
        with h5py.File(loom_path, 'r') as f:
            dset = f['matrix']  # genes x cells, dense float64
            n_genes, n_cells = dset.shape
            gene_names = f['row_attrs']['gene_names'][:].astype(str)
            cell_ids = f['col_attrs']['cell_names'][:].astype(str)

            # Read in gene-chunks, convert each to sparse immediately
            chunk_size = 500
            sparse_chunks = []
            for i in range(0, n_genes, chunk_size):
                chunk = dset[i:i + chunk_size, :]  # small dense slice (genes x cells)
                sparse_chunks.append(sparse.csr_matrix(chunk, dtype=np.float32))
            # Stack gives genes x cells sparse; transpose to cells x genes
            mat = sparse.vstack(sparse_chunks, format='csr').T

        adata = ad.AnnData(X=mat)
        adata.obs_names = pd.Index(cell_ids).astype(str)
        adata.var_names = pd.Index(gene_names).astype(str)
        os.remove(loom_path)

        samples[sample_id] = adata
        print(f"  {sample_id}: {adata.shape[0]} cells x {adata.shape[1]} genes")
    return samples


def load_GSE165722(extract_dir):
    """Load BD Rhapsody counts.tsv + cellname.txt for Tu 2022."""
    samples = {}
    for sample_id, (counts_file, cellname_file) in SAMPLE_FILE_MAP["GSE165722"].items():
        counts_path = os.path.join(extract_dir, counts_file)
        cellname_path = os.path.join(extract_dir, cellname_file)

        # counts.tsv: gene × cell (first col = gene name, columns = cell index C1, C2, ...)
        counts_df = pd.read_csv(counts_path, sep='\t', index_col=0)
        # cellname.txt: CellName, CellIndex mapping
        cellnames = pd.read_csv(cellname_path, sep='\t')

        # Map cell indices to real barcodes
        idx_to_name = dict(zip(cellnames['CellIndex'], cellnames['CellName']))
        real_names = [idx_to_name.get(c, c) for c in counts_df.columns]

        mat = sparse.csr_matrix(counts_df.values.T)  # transpose: cells x genes
        adata = ad.AnnData(X=mat)
        adata.obs_names = pd.Index(real_names).astype(str)
        adata.var_names = pd.Index(counts_df.index).astype(str)

        samples[sample_id] = adata
        print(f"  {sample_id}: {adata.shape[0]} cells x {adata.shape[1]} genes")
    return samples


def load_10x_mtx_dataset(accession, extract_dir):
    """Generic loader for datasets with standard 10x MTX triplet files."""
    samples = {}
    for sample_id, prefix in SAMPLE_FILE_MAP[accession].items():
        adata = _load_10x_mtx(prefix, extract_dir)
        adata.var_names_make_unique()
        samples[sample_id] = adata
        print(f"  {sample_id}: {adata.shape[0]} cells x {adata.shape[1]} genes")
    return samples


def load_GSE199866(extract_dir):
    """Load 10x H5 files for Cherif 2022."""
    samples = {}
    for sample_id, filename in SAMPLE_FILE_MAP["GSE199866"].items():
        h5_path = os.path.join(extract_dir, filename)
        adata = _load_10x_h5(h5_path)
        samples[sample_id] = adata
        print(f"  {sample_id}: {adata.shape[0]} cells x {adata.shape[1]} genes")
    return samples


def load_CNP0002664():
    """Load Singleron TSV matrices for Han 2022 (no tar, direct files)."""
    data_dir = RAW_DIR / "CNP0002664"
    samples = {}
    for sample_id, filename in SAMPLE_FILE_MAP["CNP0002664"].items():
        filepath = data_dir / filename
        # gene × cell TSV (first col unnamed = gene names, cols = barcodes)
        df = pd.read_csv(filepath, sep='\t', index_col=0)
        mat = sparse.csr_matrix(df.values.T)
        adata = ad.AnnData(X=mat)
        adata.obs_names = pd.Index(df.columns).astype(str)
        adata.var_names = pd.Index(df.index).astype(str)

        samples[sample_id] = adata
        print(f"  {sample_id}: {adata.shape[0]} cells x {adata.shape[1]} genes")
    return samples


def load_GSE251686(extract_dir):
    """Load nested tar.gz files for Jia 2024.

    Note: GSE251686_NP3 (GSM7986002) has a corrupt matrix.mtx on GEO
    (5.5 MB NULL block at byte 32193536, verified 2026-03-01). Skipped.
    """
    # Samples with known data corruption on GEO (verified identical MD5 on re-download)
    CORRUPT_SAMPLES = {"GSE251686_NP3"}

    samples = {}
    for sample_id, nested_tar_name in SAMPLE_FILE_MAP["GSE251686"].items():
        if sample_id in CORRUPT_SAMPLES:
            print(f"  {sample_id}: SKIPPED (corrupt matrix.mtx on GEO)")
            continue

        nested_path = os.path.join(extract_dir, nested_tar_name)
        # Extract nested tar.gz to a subdirectory
        sub_dir = os.path.join(extract_dir, sample_id)
        os.makedirs(sub_dir, exist_ok=True)
        with tarfile.open(nested_path, 'r:gz') as tf:
            tf.extractall(sub_dir)

        # Find the mtx/barcodes/genes files (might be in a subdirectory)
        mtx_file = None
        for root, dirs, files in os.walk(sub_dir):
            for f in files:
                if f == 'matrix.mtx' or f.endswith('matrix.mtx.gz'):
                    mtx_file = os.path.join(root, f)
                    break
            if mtx_file:
                break

        if mtx_file is None:
            raise FileNotFoundError(f"No matrix.mtx found for {sample_id} in {sub_dir}")

        mtx_dir = os.path.dirname(mtx_file)
        # Singleron uses genes.tsv / barcodes.tsv (not features.tsv)
        adata = sc.read_10x_mtx(mtx_dir, var_names='gene_symbols', make_unique=True)

        samples[sample_id] = adata
        print(f"  {sample_id}: {adata.shape[0]} cells x {adata.shape[1]} genes")
    return samples


def load_dataset(accession):
    """Load all samples for a dataset. Returns dict {sample_id: AnnData}."""
    print(f"\n{'='*60}")
    print(f"Loading {accession}")
    print(f"{'='*60}")

    if accession == "CNP0002664":
        return load_CNP0002664()

    # Extract tar
    extract_dir = tempfile.mkdtemp(prefix=f"{accession}_")
    try:
        _extract_tar(accession, extract_dir)

        if accession == "GSE160756":
            return load_GSE160756(extract_dir)
        elif accession == "GSE165722":
            return load_GSE165722(extract_dir)
        elif accession == "GSE199866":
            return load_GSE199866(extract_dir)
        elif accession == "GSE251686":
            return load_GSE251686(extract_dir)
        else:
            # Standard 10x MTX triplet format
            return load_10x_mtx_dataset(accession, extract_dir)
    finally:
        shutil.rmtree(extract_dir, ignore_errors=True)


# ═══════════════════════════════════════════════════════════════════════════════
# GENE NAME STANDARDIZATION
# ═══════════════════════════════════════════════════════════════════════════════

def standardize_gene_names(adata):
    """Ensure gene names are HGNC symbols, not Ensembl IDs."""
    # Check if var_names look like Ensembl IDs
    sample_names = adata.var_names[:20].tolist()
    ensembl_pattern = sum(1 for n in sample_names if n.startswith('ENSG'))

    if ensembl_pattern > len(sample_names) / 2:
        # Mostly Ensembl IDs — try to use gene_ids column if available
        if 'gene_ids' in adata.var.columns:
            # gene_ids has Ensembl, var_names should have symbols — but they're swapped
            print("  WARNING: var_names appear to be Ensembl IDs. Attempting to use symbol column.")
            # This shouldn't happen with our loading code, but handle it
        else:
            print("  WARNING: Gene names appear to be Ensembl IDs. No symbol mapping available.")
            # Could add mapping logic here if needed

    # Make gene names unique
    adata.var_names_make_unique()
    return adata


# ═══════════════════════════════════════════════════════════════════════════════
# PREPROCESSING PIPELINE
# ═══════════════════════════════════════════════════════════════════════════════

def preprocess_dataset(accession, samples_dict, metadata_df):
    """Full preprocessing pipeline for one dataset.

    Parameters
    ----------
    accession : str
    samples_dict : dict of {sample_id: AnnData}
    metadata_df : pd.DataFrame with sample metadata

    Returns
    -------
    adata : AnnData (processed)
    qc_stats : dict with QC statistics
    """
    qc_stats = {"accession": accession, "samples": {}}

    # ── Step 1: Format and merge samples ─────────────────────────────────────
    print(f"\n  Step 1: Loading and formatting...")
    adatas_list = []
    for sample_id, adata in samples_dict.items():
        adata = standardize_gene_names(adata)

        # Record original barcodes
        adata.obs['original_barcode'] = adata.obs_names.values.copy()
        adata.obs['sample_id'] = sample_id
        adata.obs['study'] = accession

        # Attach metadata
        meta_row = metadata_df[metadata_df['sample_id'] == sample_id]
        if len(meta_row) == 1:
            for col in metadata_df.columns:
                if col != 'sample_id':
                    adata.obs[col] = meta_row[col].values[0]
        else:
            print(f"  WARNING: No metadata match for {sample_id}")

        # Make obs_names unique across samples
        adata.obs_names = [f"{sample_id}_{bc}" for bc in adata.obs_names]

        qc_stats["samples"][sample_id] = {"cells_raw": adata.shape[0]}
        adatas_list.append(adata)

    # Find common genes across all samples
    if len(adatas_list) > 1:
        common_genes = set(adatas_list[0].var_names)
        for a in adatas_list[1:]:
            common_genes &= set(a.var_names)
        common_genes = sorted(common_genes)
        print(f"  Common genes across samples: {len(common_genes)}")
        # Subset in-place to avoid holding originals + copies simultaneously
        for i in range(len(adatas_list)):
            adatas_list[i] = adatas_list[i][:, common_genes].copy()

    adata = ad.concat(adatas_list, merge='same')
    del adatas_list  # free sample-level copies immediately
    adata.obs_names_make_unique()

    total_raw = adata.shape[0]
    qc_stats["cells_raw_total"] = total_raw
    qc_stats["genes_raw"] = adata.shape[1]
    print(f"  Combined: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Ensure counts are integers stored as sparse
    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(adata.X)
    # Round to int if float (some formats produce float)
    if adata.X.dtype != np.int32 and adata.X.dtype != np.int64:
        # Check if non-zero values are close to integers (sparse .data only)
        data = adata.X.data
        if np.allclose(data, np.round(data)):
            # Round and cast in-place on sparse data — avoids dense conversion
            adata.X.data = np.round(adata.X.data).astype(np.int32)
            adata.X = adata.X.astype(np.int32)
        else:
            print("  WARNING: Count matrix contains non-integer values")

    # Store raw counts
    adata.layers['counts'] = adata.X.copy()

    # ── Step 2: QC metrics ───────────────────────────────────────────────────
    print(f"  Step 2: Computing QC metrics...")
    adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
    adata.var['ribo'] = adata.var_names.str.upper().str.match(r'^RP[SL]')

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True
    )
    adata.obs.rename(columns={
        'pct_counts_mt': 'pct_counts_mt',
        'pct_counts_ribo': 'pct_counts_ribo',
    }, inplace=True)

    # ── Step 2b: QC filtering (per-sample) ───────────────────────────────────
    print(f"  Step 2b: Applying QC filters (per-sample)...")
    keep_mask = np.ones(adata.shape[0], dtype=bool)

    for sample_id in adata.obs['sample_id'].unique():
        sample_mask = adata.obs['sample_id'] == sample_id
        n_before = sample_mask.sum()

        # Gene count filters
        low_genes = adata.obs['n_genes_by_counts'] < QC_PARAMS['min_genes']
        high_genes = adata.obs['n_genes_by_counts'] > QC_PARAMS['max_genes']
        high_mt = adata.obs['pct_counts_mt'] > QC_PARAMS['max_pct_mt']
        low_counts = adata.obs['total_counts'] < QC_PARAMS['min_counts']
        high_counts = adata.obs['total_counts'] > QC_PARAMS['max_counts']

        sample_remove = sample_mask & (low_genes | high_genes | high_mt | low_counts | high_counts)
        keep_mask &= ~sample_remove

        n_removed = sample_remove.sum()
        n_after = n_before - n_removed
        qc_stats["samples"][sample_id]["cells_after_qc_filter"] = int(n_after)
        print(f"    {sample_id}: {n_before} → {n_after} cells "
              f"(removed {n_removed}: {(low_genes & sample_mask).sum()} low_genes, "
              f"{(high_genes & sample_mask).sum()} high_genes, "
              f"{(high_mt & sample_mask).sum()} high_mt, "
              f"{(low_counts & sample_mask).sum()} low_counts, "
              f"{(high_counts & sample_mask).sum()} high_counts)")

    adata = adata[keep_mask].copy()
    print(f"  After QC filtering: {adata.shape[0]} cells")

    # ── Step 2c: Doublet detection (per-sample) ─────────────────────────────
    print(f"  Step 2c: Running Scrublet doublet detection...")
    adata.obs['doublet_score'] = 0.0
    adata.obs['predicted_doublet'] = False

    for sample_id in adata.obs['sample_id'].unique():
        sample_mask = adata.obs['sample_id'] == sample_id
        sample_adata = adata[sample_mask]

        n_cells = sample_adata.shape[0]
        if n_cells < 50:
            print(f"    {sample_id}: Too few cells ({n_cells}) for Scrublet, skipping")
            continue

        try:
            counts = sample_adata.layers['counts']
            if sparse.issparse(counts):
                counts = counts.toarray()
            scrub = scr.Scrublet(counts, expected_doublet_rate=min(0.06, n_cells / 1000 * 0.008))
            doublet_scores, predicted = scrub.scrub_doublets(verbose=False)
            adata.obs.loc[sample_mask, 'doublet_score'] = doublet_scores
            adata.obs.loc[sample_mask, 'predicted_doublet'] = (
                doublet_scores > QC_PARAMS['doublet_score_threshold']
            )
            n_doublets = (doublet_scores > QC_PARAMS['doublet_score_threshold']).sum()
            print(f"    {sample_id}: {n_doublets}/{n_cells} predicted doublets")
        except Exception as e:
            print(f"    {sample_id}: Scrublet failed ({e}), skipping doublet detection")

    # Remove doublets
    n_before_doublet = adata.shape[0]
    adata = adata[~adata.obs['predicted_doublet']].copy()
    n_doublets_removed = n_before_doublet - adata.shape[0]
    print(f"  After doublet removal: {adata.shape[0]} cells (removed {n_doublets_removed})")

    # Update per-sample stats after doublet removal
    for sample_id in adata.obs['sample_id'].unique():
        n_after = (adata.obs['sample_id'] == sample_id).sum()
        qc_stats["samples"][sample_id]["cells_after_doublet"] = int(n_after)

    # ── Step 2d: Gene filtering ──────────────────────────────────────────────
    print(f"  Step 2d: Filtering genes...")
    genes_before = adata.shape[1]
    sc.pp.filter_genes(adata, min_cells=QC_PARAMS['min_cells_per_gene'])
    print(f"  Genes: {genes_before} → {adata.shape[1]}")

    qc_stats["cells_after_qc"] = adata.shape[0]
    qc_stats["genes_after_qc"] = adata.shape[1]

    # ── Step 3: Normalization (SCTransform equivalent) ───────────────────────
    print(f"  Step 3: Normalizing (Pearson residuals / SCTransform equivalent)...")
    # Update counts layer after gene filtering
    adata.layers['counts'] = adata.X.copy()
    # Pearson residual normalization: the scanpy equivalent of Seurat's SCTransform.
    # Variance-stabilizes counts without separate normalize → log → scale steps.
    sc.experimental.pp.normalize_pearson_residuals(adata)

    # ── Step 4: Feature selection ────────────────────────────────────────────
    print(f"  Step 4: Identifying highly variable genes...")
    n_samples = adata.obs['sample_id'].nunique()
    n_top = min(3000, adata.shape[1] - 1)

    try:
        if n_samples > 1:
            sc.pp.highly_variable_genes(
                adata, n_top_genes=n_top, flavor='seurat_v3',
                batch_key='sample_id', layer='counts'
            )
        else:
            sc.pp.highly_variable_genes(
                adata, n_top_genes=n_top, flavor='seurat_v3',
                layer='counts'
            )
    except Exception as e:
        print(f"  WARNING: seurat_v3 HVG failed ({e}), falling back to seurat method")
        if n_samples > 1:
            sc.pp.highly_variable_genes(
                adata, n_top_genes=n_top, flavor='seurat',
                batch_key='sample_id'
            )
        else:
            sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor='seurat')

    n_hvg = adata.var['highly_variable'].sum()
    print(f"  HVGs selected: {n_hvg}")
    qc_stats["n_hvgs"] = int(n_hvg)

    # ── Step 5: Dimensionality reduction ─────────────────────────────────────
    print(f"  Step 5: PCA, neighbors, UMAP...")
    # NOTE: When using Pearson residuals (SCTransform equivalent), the separate
    # scaling step is not needed since Pearson residuals are already
    # variance-stabilized. We skip sc.pp.scale() and run PCA directly.
    adata_pca = adata[:, adata.var['highly_variable']].copy()

    n_comps = min(50, min(adata_pca.shape) - 1)
    sc.tl.pca(adata_pca, n_comps=n_comps, svd_solver='arpack')

    # Transfer PCA results back
    adata.obsm['X_pca'] = adata_pca.obsm['X_pca']
    adata.uns['pca'] = adata_pca.uns['pca']
    adata.varm['PCs'] = np.zeros((adata.shape[1], n_comps))
    hvg_idx = np.where(adata.var['highly_variable'])[0]
    adata.varm['PCs'][hvg_idx] = adata_pca.varm['PCs']

    # Use all 50 PCs (or fewer if matrix is smaller)
    n_pcs = n_comps
    print(f"  PCs used: {n_pcs} (all computed components)")
    qc_stats["n_pcs"] = n_pcs

    sc.pp.neighbors(adata, n_pcs=n_pcs)
    sc.tl.umap(adata)

    del adata_pca

    # ── Step 6: Clustering ───────────────────────────────────────────────────
    print(f"  Step 6: Leiden clustering...")
    from sklearn.metrics import silhouette_score

    resolution_scores = {}
    for res in LEIDEN_RESOLUTIONS:
        key = f"leiden_res_{res}"
        sc.tl.leiden(adata, resolution=res, key_added=key, flavor='igraph',
                     n_iterations=2)
        n_clusters = adata.obs[key].nunique()
        print(f"    Resolution {res}: {n_clusters} clusters")

        # Score each resolution: silhouette on PCA embedding + modularity
        labels = adata.obs[key].astype(int).values
        if n_clusters > 1 and n_clusters < adata.shape[0]:
            try:
                # Subsample for silhouette if dataset is large (>50k cells)
                if adata.shape[0] > 50000:
                    idx = np.random.choice(adata.shape[0], 50000, replace=False)
                    sil = silhouette_score(adata.obsm['X_pca'][idx], labels[idx])
                else:
                    sil = silhouette_score(adata.obsm['X_pca'], labels)
            except Exception:
                sil = -1.0
            # Modularity from the Leiden partition (stored by scanpy)
            mod_key = f"leiden_res_{res}"
            mod = adata.uns.get(f'{mod_key}', {}).get('params', {}).get('modularity', None)
            # Combined score: silhouette is primary metric
            resolution_scores[res] = sil
            print(f"      Silhouette score: {sil:.4f}")

    # Select optimal resolution by highest silhouette score
    if resolution_scores:
        best_res = max(resolution_scores, key=resolution_scores.get)
        best_sil = resolution_scores[best_res]
        print(f"  Optimal resolution: {best_res} (silhouette={best_sil:.4f})")
        working_resolution = best_res
    else:
        print(f"  WARNING: Could not compute silhouette scores, defaulting to res=0.5")
        working_resolution = WORKING_RESOLUTION

    qc_stats["working_resolution"] = working_resolution
    qc_stats["n_clusters_working"] = int(
        adata.obs[f"leiden_res_{working_resolution}"].nunique()
    )

    # ── Step 7: Marker genes ─────────────────────────────────────────────────
    print(f"  Step 7: Computing marker genes...")
    working_key = f"leiden_res_{working_resolution}"
    sc.tl.rank_genes_groups(adata, groupby=working_key, method='wilcoxon',
                            use_raw=False)

    # Save markers to TSV
    markers_list = []
    for cluster in adata.obs[working_key].cat.categories:
        result = sc.get.rank_genes_groups_df(adata, group=str(cluster))
        result = result.head(50)
        result['cluster'] = cluster
        markers_list.append(result)

    markers_df = pd.concat(markers_list, ignore_index=True)
    markers_path = ANN_DIR / f"{accession}_markers.tsv"
    markers_df.to_csv(markers_path, sep='\t', index=False)
    print(f"  Markers saved to {markers_path}")

    # ── Step 8: Preliminary cell type labels ─────────────────────────────────
    print(f"  Step 8: Preliminary cell type scoring...")
    # Score each cell type
    for ct_name, markers in CELL_TYPE_MARKERS.items():
        available_markers = [g for g in markers if g in adata.var_names]
        if len(available_markers) >= 1:
            sc.tl.score_genes(adata, gene_list=available_markers,
                              score_name=f"score_{ct_name}")
        else:
            adata.obs[f"score_{ct_name}"] = 0.0

    # Assign preliminary labels based on highest score
    score_cols = [f"score_{ct}" for ct in CELL_TYPE_MARKERS.keys()]
    score_matrix = adata.obs[score_cols].values
    ct_names = list(CELL_TYPE_MARKERS.keys())

    max_scores = np.max(score_matrix, axis=1)
    best_idx = np.argmax(score_matrix, axis=1)
    labels = np.array([ct_names[i] for i in best_idx])

    # Confidence levels
    confidence = np.where(
        max_scores > 0.5, 'high',
        np.where(max_scores > 0.2, 'medium', 'low')
    )

    adata.obs['cell_type_preliminary'] = labels
    adata.obs['cell_type_confidence'] = confidence

    # Summary
    print(f"  Cell type distribution:")
    for ct in sorted(adata.obs['cell_type_preliminary'].unique()):
        n = (adata.obs['cell_type_preliminary'] == ct).sum()
        pct = n / adata.shape[0] * 100
        print(f"    {ct}: {n} ({pct:.1f}%)")

    qc_stats["cell_type_counts"] = adata.obs['cell_type_preliminary'].value_counts().to_dict()

    return adata, qc_stats


# ═══════════════════════════════════════════════════════════════════════════════
# QC REPORT GENERATION
# ═══════════════════════════════════════════════════════════════════════════════

QC_REPORT_TEMPLATE = """<!DOCTYPE html>
<html>
<head>
<title>QC Report: {{ accession }}</title>
<style>
  body { font-family: Arial, sans-serif; max-width: 1200px; margin: 0 auto; padding: 20px; }
  h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
  h2 { color: #34495e; margin-top: 30px; }
  table { border-collapse: collapse; width: 100%; margin: 10px 0; }
  th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
  th { background-color: #3498db; color: white; }
  tr:nth-child(even) { background-color: #f2f2f2; }
  img { max-width: 100%; height: auto; margin: 10px 0; border: 1px solid #ddd; }
  .stats-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; }
  .stat-box { background: #ecf0f1; padding: 15px; border-radius: 5px; }
  .stat-box h3 { margin: 0 0 5px 0; color: #7f8c8d; font-size: 12px; }
  .stat-box p { margin: 0; font-size: 24px; font-weight: bold; color: #2c3e50; }
  .warning { color: #e74c3c; font-weight: bold; }
</style>
</head>
<body>
<h1>QC Report: {{ accession }}</h1>
<p>Generated: {{ date }}</p>

<div class="stats-grid">
  <div class="stat-box"><h3>Cells (raw)</h3><p>{{ stats.cells_raw_total | int }}</p></div>
  <div class="stat-box"><h3>Cells (after QC)</h3><p>{{ stats.cells_after_qc | int }}</p></div>
  <div class="stat-box"><h3>Retention</h3><p>{{ "%.1f" | format(stats.cells_after_qc / stats.cells_raw_total * 100) }}%</p></div>
  <div class="stat-box"><h3>Genes (after QC)</h3><p>{{ stats.genes_after_qc | int }}</p></div>
  <div class="stat-box"><h3>HVGs</h3><p>{{ stats.n_hvgs | int }}</p></div>
  <div class="stat-box"><h3>PCs used</h3><p>{{ stats.n_pcs | int }}</p></div>
  <div class="stat-box"><h3>Clusters (res={{ working_res }})</h3><p>{{ stats.n_clusters_working | int }}</p></div>
</div>

<h2>Per-Sample Cell Counts</h2>
<table>
  <tr><th>Sample</th><th>Raw</th><th>After QC</th><th>After Doublet</th><th>Retention %</th></tr>
  {% for sid, ss in stats.samples.items() %}
  <tr>
    <td>{{ sid }}</td>
    <td>{{ ss.cells_raw }}</td>
    <td>{{ ss.get('cells_after_qc_filter', 'N/A') }}</td>
    <td>{{ ss.get('cells_after_doublet', 'N/A') }}</td>
    <td>
      {% if ss.get('cells_after_doublet') and ss.cells_raw > 0 %}
        {{ "%.1f" | format(ss.cells_after_doublet / ss.cells_raw * 100) }}%
        {% if ss.cells_after_doublet / ss.cells_raw < 0.5 %}
          <span class="warning"> ⚠ LOW</span>
        {% endif %}
      {% else %}
        N/A
      {% endif %}
    </td>
  </tr>
  {% endfor %}
</table>

<h2>QC Metric Distributions</h2>
<img src="{{ accession }}_qc_violins.png" alt="QC violin plots">

<h2>UMAP Embeddings</h2>
<img src="{{ accession }}_umap_sample.png" alt="UMAP by sample">
<img src="{{ accession }}_umap_clusters.png" alt="UMAP by clusters">
<img src="{{ accession }}_umap_celltype.png" alt="UMAP by cell type">
<img src="{{ accession }}_umap_qc.png" alt="UMAP QC metrics">

<h2>Marker Expression</h2>
<img src="{{ accession }}_dotplot_markers.png" alt="Marker dot plot">

<h2>Top 10 Markers per Cluster</h2>
<img src="{{ accession }}_top_markers.png" alt="Top markers heatmap">

<h2>Cell Type Distribution</h2>
<table>
  <tr><th>Cell Type</th><th>Count</th><th>Percentage</th></tr>
  {% for ct, count in stats.cell_type_counts.items() | sort(attribute='1', reverse=True) %}
  <tr>
    <td>{{ ct }}</td>
    <td>{{ count }}</td>
    <td>{{ "%.1f" | format(count / stats.cells_after_qc * 100) }}%</td>
  </tr>
  {% endfor %}
</table>

</body>
</html>
"""

OVERVIEW_TEMPLATE = """<!DOCTYPE html>
<html>
<head>
<title>Cross-Dataset QC Overview</title>
<style>
  body { font-family: Arial, sans-serif; max-width: 1400px; margin: 0 auto; padding: 20px; }
  h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
  table { border-collapse: collapse; width: 100%; margin: 10px 0; }
  th, td { border: 1px solid #ddd; padding: 8px; text-align: right; }
  th { background-color: #3498db; color: white; text-align: left; }
  td:first-child { text-align: left; font-weight: bold; }
  tr:nth-child(even) { background-color: #f2f2f2; }
  .warning { color: #e74c3c; font-weight: bold; }
  img { max-width: 100%; margin: 10px 0; }
</style>
</head>
<body>
<h1>Cross-Dataset QC Overview — IVD scRNA-seq Atlas</h1>
<p>Generated: {{ date }}</p>
<p>Datasets processed: {{ n_datasets }}</p>

<h2>Summary Table</h2>
<table>
  <tr>
    <th>Dataset</th><th>Samples</th><th>Cells (raw)</th><th>Cells (QC)</th>
    <th>Retention %</th><th>Genes</th><th>HVGs</th><th>PCs</th><th>Clusters</th>
  </tr>
  {% for row in summary %}
  <tr>
    <td>{{ row.accession }}</td>
    <td>{{ row.n_samples }}</td>
    <td>{{ row.cells_raw_total | int }}</td>
    <td>{{ row.cells_after_qc | int }}</td>
    <td>{{ "%.1f" | format(row.cells_after_qc / row.cells_raw_total * 100) }}%</td>
    <td>{{ row.genes_after_qc | int }}</td>
    <td>{{ row.n_hvgs | int }}</td>
    <td>{{ row.n_pcs | int }}</td>
    <td>{{ row.n_clusters_working | int }}</td>
  </tr>
  {% endfor %}
</table>

<h2>Cross-Dataset Comparison</h2>
<img src="overview_cells_comparison.png" alt="Cell count comparison">
<img src="overview_qc_comparison.png" alt="QC metrics comparison">

</body>
</html>
"""


def generate_qc_plots(adata, accession, qc_stats):
    """Generate all QC plots for one dataset."""
    sc.settings.figdir = str(QC_DIR)
    plt.rcParams['figure.dpi'] = 100

    working_key = f"leiden_res_{WORKING_RESOLUTION}"

    # 1. QC violin plots
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    for ax, metric, title in zip(
        axes,
        ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
        ['Genes per cell', 'UMI counts per cell', '% Mitochondrial']
    ):
        # Per-sample violin
        sample_ids = adata.obs['sample_id'].unique()
        data_by_sample = [adata.obs.loc[adata.obs['sample_id'] == s, metric].values
                          for s in sample_ids]
        parts = ax.violinplot(data_by_sample, positions=range(len(sample_ids)),
                              showmedians=True)
        ax.set_xticks(range(len(sample_ids)))
        short_labels = [s.replace(f"{accession}_", "") for s in sample_ids]
        ax.set_xticklabels(short_labels, rotation=45, ha='right', fontsize=8)
        ax.set_title(title)
        ax.set_ylabel(metric)
    plt.tight_layout()
    plt.savefig(QC_DIR / f"{accession}_qc_violins.png", bbox_inches='tight')
    plt.close()

    # 2. UMAP by sample
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color='sample_id', ax=ax, show=False, title=f'{accession} — Sample')
    plt.tight_layout()
    plt.savefig(QC_DIR / f"{accession}_umap_sample.png", bbox_inches='tight')
    plt.close()

    # 3. UMAP by clusters (multiple resolutions)
    n_res = len(LEIDEN_RESOLUTIONS)
    fig, axes = plt.subplots(1, n_res, figsize=(6*n_res, 5))
    if n_res == 1:
        axes = [axes]
    for ax, res in zip(axes, LEIDEN_RESOLUTIONS):
        key = f"leiden_res_{res}"
        sc.pl.umap(adata, color=key, ax=ax, show=False,
                    title=f'Leiden res={res}', legend_loc='on data',
                    legend_fontsize=6)
    plt.tight_layout()
    plt.savefig(QC_DIR / f"{accession}_umap_clusters.png", bbox_inches='tight')
    plt.close()

    # 4. UMAP by cell type
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color='cell_type_preliminary', ax=ax, show=False,
                title=f'{accession} — Preliminary Cell Type')
    plt.tight_layout()
    plt.savefig(QC_DIR / f"{accession}_umap_celltype.png", bbox_inches='tight')
    plt.close()

    # 5. UMAP by QC metrics
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    for ax, metric, title in zip(
        axes,
        ['pct_counts_mt', 'n_genes_by_counts', 'total_counts'],
        ['% MT', 'Genes', 'UMI counts']
    ):
        sc.pl.umap(adata, color=metric, ax=ax, show=False, title=title,
                    color_map='viridis')
    plt.tight_layout()
    plt.savefig(QC_DIR / f"{accession}_umap_qc.png", bbox_inches='tight')
    plt.close()

    # 6. Dot plot of canonical markers
    all_markers = []
    for ct, genes in CELL_TYPE_MARKERS.items():
        available = [g for g in genes if g in adata.var_names]
        all_markers.extend(available)
    all_markers = list(dict.fromkeys(all_markers))  # dedupe preserving order

    if all_markers:
        fig, ax = plt.subplots(figsize=(max(12, len(all_markers)*0.5), 6))
        try:
            sc.pl.dotplot(adata, var_names=all_markers, groupby=working_key,
                          ax=ax, show=False)
        except Exception:
            sc.pl.dotplot(adata, var_names=all_markers[:30], groupby=working_key,
                          show=False)
        plt.tight_layout()
        plt.savefig(QC_DIR / f"{accession}_dotplot_markers.png", bbox_inches='tight')
        plt.close()
    else:
        # Create placeholder
        fig, ax = plt.subplots(figsize=(6, 2))
        ax.text(0.5, 0.5, "No canonical markers found in gene set",
                ha='center', va='center', fontsize=14)
        ax.axis('off')
        plt.savefig(QC_DIR / f"{accession}_dotplot_markers.png", bbox_inches='tight')
        plt.close()

    # 7. Top markers heatmap
    try:
        fig, ax = plt.subplots(figsize=(12, 8))
        sc.pl.rank_genes_groups(adata, n_genes=10, ax=ax, show=False)
        plt.tight_layout()
        plt.savefig(QC_DIR / f"{accession}_top_markers.png", bbox_inches='tight')
        plt.close()
    except Exception:
        fig, ax = plt.subplots(figsize=(6, 2))
        ax.text(0.5, 0.5, "Could not generate marker heatmap",
                ha='center', va='center')
        ax.axis('off')
        plt.savefig(QC_DIR / f"{accession}_top_markers.png", bbox_inches='tight')
        plt.close()


def generate_qc_report(accession, qc_stats):
    """Generate HTML QC report for one dataset."""
    template = Template(QC_REPORT_TEMPLATE)
    html = template.render(
        accession=accession,
        date=datetime.now().strftime("%Y-%m-%d %H:%M"),
        stats=qc_stats,
        working_res=WORKING_RESOLUTION,
    )
    report_path = QC_DIR / f"{accession}_qc_report.html"
    report_path.write_text(html)
    print(f"  QC report: {report_path}")


def generate_overview(all_stats):
    """Generate cross-dataset overview report."""
    # Summary TSV
    summary_rows = []
    for stats in all_stats:
        # stats may be a dict from current run (has "samples" key) or
        # a dict from loading the TSV (has "n_samples" key directly)
        if "samples" in stats:
            n_samples = len(stats["samples"])
        else:
            n_samples = int(stats.get("n_samples", 0))
        row = {
            "accession": stats["accession"],
            "n_samples": n_samples,
            "cells_raw_total": stats["cells_raw_total"],
            "cells_after_qc": stats["cells_after_qc"],
            "retention_pct": stats["cells_after_qc"] / stats["cells_raw_total"] * 100,
            "genes_raw": stats.get("genes_raw", 0),
            "genes_after_qc": stats["genes_after_qc"],
            "n_hvgs": stats["n_hvgs"],
            "n_pcs": stats["n_pcs"],
            "n_clusters_working": stats["n_clusters_working"],
        }
        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(QC_DIR / "qc_summary.tsv", sep='\t', index=False)
    print(f"\nQC summary: {QC_DIR / 'qc_summary.tsv'}")

    # Overview plots
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    # Cell count comparison
    ax = axes[0]
    x = range(len(summary_df))
    ax.bar(x, summary_df['cells_raw_total'], alpha=0.5, label='Raw')
    ax.bar(x, summary_df['cells_after_qc'], alpha=0.8, label='After QC')
    ax.set_xticks(x)
    ax.set_xticklabels(summary_df['accession'], rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Cell count')
    ax.set_title('Cell Counts: Raw vs After QC')
    ax.legend()

    # Retention comparison
    ax = axes[1]
    colors = ['#e74c3c' if r < 50 else '#f39c12' if r < 70 else '#27ae60'
              for r in summary_df['retention_pct']]
    ax.bar(x, summary_df['retention_pct'], color=colors)
    ax.axhline(y=50, color='red', linestyle='--', alpha=0.5, label='50% threshold')
    ax.set_xticks(x)
    ax.set_xticklabels(summary_df['accession'], rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Retention %')
    ax.set_title('Cell Retention After QC')
    ax.set_ylim(0, 100)
    ax.legend()

    plt.tight_layout()
    plt.savefig(QC_DIR / "overview_cells_comparison.png", bbox_inches='tight')
    plt.close()

    # QC metrics comparison
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(x, summary_df['n_clusters_working'], color='#3498db')
    ax.set_xticks(x)
    ax.set_xticklabels(summary_df['accession'], rotation=45, ha='right', fontsize=8)
    ax.set_ylabel(f'Clusters (Leiden res={WORKING_RESOLUTION})')
    ax.set_title('Number of Clusters per Dataset')
    plt.tight_layout()
    plt.savefig(QC_DIR / "overview_qc_comparison.png", bbox_inches='tight')
    plt.close()

    # HTML overview
    template = Template(OVERVIEW_TEMPLATE)
    html = template.render(
        date=datetime.now().strftime("%Y-%m-%d %H:%M"),
        n_datasets=len(all_stats),
        summary=summary_rows,
    )
    (QC_DIR / "qc_overview.html").write_text(html)
    print(f"QC overview: {QC_DIR / 'qc_overview.html'}")


# ═══════════════════════════════════════════════════════════════════════════════
# VALIDATION
# ═══════════════════════════════════════════════════════════════════════════════

def validate_dataset(accession, metadata_df):
    """Run per-dataset validation checks. Returns (pass, messages)."""
    messages = []
    all_pass = True

    h5ad_path = PROC_DIR / f"{accession}.h5ad"

    # Check 1: File exists and loadable
    if not h5ad_path.exists():
        messages.append(f"FAIL: {h5ad_path} does not exist")
        return False, messages

    try:
        adata = sc.read_h5ad(h5ad_path)
        messages.append(f"PASS: {accession}.h5ad loaded ({adata.shape[0]} cells x {adata.shape[1]} genes)")
    except Exception as e:
        messages.append(f"FAIL: Cannot load {h5ad_path}: {e}")
        return False, messages

    # Check 2: Raw counts in layers
    if 'counts' in adata.layers:
        counts = adata.layers['counts']
        if sparse.issparse(counts):
            sample_data = counts.data[:1000]
        else:
            sample_data = counts.ravel()[:1000]
        if np.allclose(sample_data, np.round(sample_data)):
            messages.append("PASS: .layers['counts'] contains integer-like values")
        else:
            messages.append("FAIL: .layers['counts'] contains non-integer values")
            all_pass = False
    else:
        messages.append("FAIL: .layers['counts'] not found")
        all_pass = False

    # Check 3: Required metadata fields
    required_obs = ['sample_id', 'study', 'original_barcode', 'study_accession',
                    'condition_harmonized', 'compartment']
    missing = [f for f in required_obs if f not in adata.obs.columns]
    if missing:
        messages.append(f"FAIL: Missing obs fields: {missing}")
        all_pass = False
    else:
        messages.append("PASS: All required metadata fields present")

    # Check 4: Cell retention >= 50%
    sample_map = SAMPLE_FILE_MAP.get(accession, {})
    total_raw = 0
    for sample_id in sample_map:
        meta_row = metadata_df[metadata_df['sample_id'] == sample_id]
        if len(meta_row) > 0:
            n_raw = meta_row['n_cells_raw'].values[0]
            if pd.notna(n_raw) and str(n_raw) != 'NA':
                total_raw += int(float(n_raw))

    if total_raw > 0:
        retention = adata.shape[0] / total_raw * 100
        if retention >= 50:
            messages.append(f"PASS: Cell retention {retention:.1f}% (>= 50%)")
        else:
            messages.append(f"WARNING: Cell retention {retention:.1f}% (< 50%) — may indicate aggressive filtering")
    else:
        messages.append(f"INFO: Cannot check retention (raw cell counts not in metadata)")

    # Check 5: At least 200 cells per sample
    for sample_id in adata.obs['sample_id'].unique():
        n_cells = (adata.obs['sample_id'] == sample_id).sum()
        if n_cells < 200:
            messages.append(f"WARNING: {sample_id} has only {n_cells} cells (< 200)")

    # Check 6: Canonical marker expression
    chondro_markers = ['ACAN', 'COL2A1', 'COL1A1']
    working_key = f"leiden_res_{WORKING_RESOLUTION}"
    found_markers = [g for g in chondro_markers if g in adata.var_names]
    if found_markers:
        # Check if at least one cluster shows enrichment
        for gene in found_markers:
            max_expr = 0
            for cluster in adata.obs[working_key].unique():
                mask = adata.obs[working_key] == cluster
                if gene in adata.var_names:
                    expr = np.mean(adata[mask, gene].X.toarray() if sparse.issparse(adata.X)
                                   else adata[mask, gene].X)
                    max_expr = max(max_expr, expr)
            if max_expr > 0.1:
                messages.append(f"PASS: {gene} expressed (max cluster mean = {max_expr:.2f})")
                break
        else:
            messages.append(f"WARNING: Canonical markers {found_markers} not strongly expressed in any cluster")
    else:
        messages.append(f"WARNING: None of {chondro_markers} found in gene names")

    # Check 7: Immune marker localization
    immune_markers = ['CD68', 'CD3D']
    for gene in immune_markers:
        if gene in adata.var_names:
            expr_per_cluster = {}
            for cluster in adata.obs[working_key].unique():
                mask = adata.obs[working_key] == cluster
                expr = np.mean(adata[mask, gene].X.toarray() if sparse.issparse(adata.X)
                               else adata[mask, gene].X)
                expr_per_cluster[cluster] = expr
            n_expressing = sum(1 for v in expr_per_cluster.values() if v > 0.1)
            total_clusters = len(expr_per_cluster)
            if n_expressing > total_clusters * 0.5:
                messages.append(f"WARNING: {gene} diffusely expressed across {n_expressing}/{total_clusters} clusters")
            elif n_expressing > 0:
                messages.append(f"PASS: {gene} localized to {n_expressing}/{total_clusters} clusters")
            else:
                messages.append(f"INFO: {gene} not substantially expressed (absent in this tissue type)")

    # Check 8: No single cluster >80% of cells
    cluster_sizes = adata.obs[working_key].value_counts(normalize=True)
    max_cluster_pct = cluster_sizes.max() * 100
    if max_cluster_pct > 80:
        messages.append(f"FAIL: Largest cluster contains {max_cluster_pct:.1f}% of cells (>80%)")
        all_pass = False
    else:
        messages.append(f"PASS: Largest cluster = {max_cluster_pct:.1f}% (< 80%)")

    # Check 9: QC report exists
    report_path = QC_DIR / f"{accession}_qc_report.html"
    if report_path.exists():
        messages.append(f"PASS: QC report exists")
    else:
        messages.append(f"FAIL: QC report not found at {report_path}")
        all_pass = False

    # Check 10: QC summary updated
    summary_path = QC_DIR / "qc_summary.tsv"
    if summary_path.exists():
        summary_df = pd.read_csv(summary_path, sep='\t')
        if accession in summary_df['accession'].values:
            messages.append(f"PASS: {accession} in qc_summary.tsv")
        else:
            messages.append(f"FAIL: {accession} not in qc_summary.tsv")
            all_pass = False
    else:
        messages.append(f"FAIL: qc_summary.tsv not found")
        all_pass = False

    del adata
    return all_pass, messages


def validate_all(metadata_df):
    """Run all validation checks."""
    print("\n" + "="*60)
    print("VALIDATION")
    print("="*60)

    all_pass = True
    for accession in ALL_ACCESSIONS:
        print(f"\n--- {accession} ---")
        passed, messages = validate_dataset(accession, metadata_df)
        for msg in messages:
            print(f"  {msg}")
        if not passed:
            all_pass = False

    # Cross-dataset checks
    print(f"\n--- Cross-dataset ---")
    n_processed = sum(1 for a in ALL_ACCESSIONS
                      if (PROC_DIR / f"{a}.h5ad").exists())
    if n_processed == len(ALL_ACCESSIONS):
        print(f"  PASS: All {len(ALL_ACCESSIONS)} datasets processed")
    else:
        print(f"  FAIL: Only {n_processed}/{len(ALL_ACCESSIONS)} datasets processed")
        all_pass = False

    overview_path = QC_DIR / "qc_overview.html"
    if overview_path.exists():
        print(f"  PASS: Cross-dataset overview exists")
    else:
        print(f"  FAIL: qc_overview.html not found")
        all_pass = False

    print(f"\n{'='*60}")
    print(f"OVERALL: {'ALL CHECKS PASSED' if all_pass else 'SOME CHECKS FAILED'}")
    print(f"{'='*60}")

    return all_pass


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    args = sys.argv[1:]

    # Load metadata
    metadata_df = pd.read_csv(META_FILE, sep='\t')

    if args and args[0] == '--validate-only':
        validate_all(metadata_df)
        return

    # Determine which datasets to process
    if args:
        accessions = [a for a in args if a in ALL_ACCESSIONS]
        if not accessions:
            print(f"Unknown accession(s): {args}")
            print(f"Available: {ALL_ACCESSIONS}")
            sys.exit(1)
    else:
        accessions = ALL_ACCESSIONS

    all_stats = []

    for accession in accessions:
        output_path = PROC_DIR / f"{accession}.h5ad"

        try:
            # Load
            samples = load_dataset(accession)

            # Preprocess (samples dict consumed — freed inside)
            adata, qc_stats = preprocess_dataset(accession, samples, metadata_df)
            del samples

            # Generate plots
            print(f"\n  Step 9: Generating QC report...")
            generate_qc_plots(adata, accession, qc_stats)
            generate_qc_report(accession, qc_stats)

            # Save
            print(f"  Saving {output_path}...")
            adata.write_h5ad(output_path)
            print(f"  Saved: {adata.shape[0]} cells x {adata.shape[1]} genes")

            all_stats.append(qc_stats)

            del adata

        except Exception as e:
            print(f"\nERROR processing {accession}: {e}")
            import traceback
            traceback.print_exc()
            # Continue with other datasets
            continue

    # Load any previously computed stats for datasets not reprocessed
    if len(accessions) < len(ALL_ACCESSIONS):
        existing_summary = QC_DIR / "qc_summary.tsv"
        if existing_summary.exists():
            prev = pd.read_csv(existing_summary, sep='\t')
            processed_accessions = {s["accession"] for s in all_stats}
            for _, row in prev.iterrows():
                if row['accession'] not in processed_accessions:
                    all_stats.append(row.to_dict())

    # Generate overview
    if all_stats:
        generate_overview(all_stats)

    # Run validation
    validate_all(metadata_df)


if __name__ == "__main__":
    main()
