#!/usr/bin/env python3
"""Module 05: Cross-dataset integration for IVD scRNA-seq atlas.

Combines cells across 12 studies into shared representations using a tiered
strategy: Tier 1 (non-resident immune/endothelial/pericyte) gets standard
scVI integration; Tier 2 (resident NP/AF) is benchmarked across four
approaches (scVI, scANVI, Harmony, BBKNN).

Usage:
    python3 scripts/05_integration.py                  # All tiers, approaches A-D
    python3 scripts/05_integration.py --tier1-only     # Non-resident only
    python3 scripts/05_integration.py --tier2-only     # Resident only (NP + AF)
    python3 scripts/05_integration.py --approach C     # Single approach for Tier 2
    python3 scripts/05_integration.py --validate-only  # Validation only
    python3 scripts/05_integration.py --force          # Re-run even if outputs exist
"""

import gc
import sys
import os
import warnings
import traceback
from pathlib import Path
from datetime import datetime
from collections import Counter

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
from scipy.stats import pearsonr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from jinja2 import Template

warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning, module='scanpy')
warnings.filterwarnings('ignore', category=UserWarning, module='scvi')
warnings.filterwarnings('ignore', category=DeprecationWarning)

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent.parent
PROC_DIR = BASE / "data" / "processed"
INT_DIR = BASE / "data" / "integrated"
MODEL_DIR = INT_DIR / "models"
RESULTS_DIR = BASE / "results" / "integration"

INT_DIR.mkdir(parents=True, exist_ok=True)
MODEL_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ── Dataset ordering ─────────────────────────────────────────────────────────
ALL_ACCESSIONS = [
    "GSE160756", "GSE165722", "GSE189916", "GSE199866", "GSE205535",
    "CNP0002664", "GSE233666", "GSE244889", "GSE251686", "GSE255768",
    "GSE230809", "GSE242443",
]

# ── Tier classification ──────────────────────────────────────────────────────
# Non-resident cell types (CellTypist-derived or marker-based names)
NONRESIDENT_KEYWORDS = [
    "Endothelial", "Macrophage", "T cell", "T_cell", "NK cell", "NK cells",
    "B cell", "B_cell", "Mast", "Pericyte", "SMC", "DC", "monocyte",
    "Monocyte", "Plasma", "Neutrophil", "pre-B", "Pro-B",
    "Tem/", "Tcm/", "CD16+",
]

NP_KEYWORDS = ["NP_"]
AF_KEYWORDS = ["AF_"]
EP_KEYWORDS = ["EP_"]
FIBROBLAST_KEYWORDS = ["Fibroblast"]

# Signature score columns for continuum metrics
NP_SCORE_COLS = [
    "score_NP_notochordal", "score_NP_mature_chondrocyte",
    "score_NP_stressed_degenerative", "score_NP_fibrocartilaginous",
]
AF_SCORE_COLS = [
    "score_AF_inner", "score_AF_outer", "score_AF_mechanical_stress",
]

# obs columns to keep during concatenation
OBS_COLS_KEEP = [
    "sample_id", "study", "compartment", "condition_harmonized",
    "cell_type_final", "cell_type_confidence",
    "pct_counts_mt",
    "score_NP_notochordal", "score_NP_mature_chondrocyte",
    "score_NP_stressed_degenerative", "score_NP_fibrocartilaginous",
    "score_AF_inner", "score_AF_outer", "score_AF_mechanical_stress",
    "score_EP_hyaline_cartilage", "score_EP_ossification",
]


# ═══════════════════════════════════════════════════════════════════════════════
# CELL TIER CLASSIFICATION
# ═══════════════════════════════════════════════════════════════════════════════

def classify_cell_tier(label):
    """Classify a cell_type_final label into a tier.

    Returns one of: "nonresident", "NP", "AF", "EP", "fibroblast", "other".
    """
    if pd.isna(label) or label == "unassigned":
        return "other"
    s = str(label)
    for kw in NONRESIDENT_KEYWORDS:
        if kw in s:
            return "nonresident"
    for kw in NP_KEYWORDS:
        if s.startswith(kw):
            return "NP"
    for kw in AF_KEYWORDS:
        if s.startswith(kw):
            return "AF"
    for kw in EP_KEYWORDS:
        if s.startswith(kw):
            return "EP"
    for kw in FIBROBLAST_KEYWORDS:
        if kw in s:
            return "fibroblast"
    return "other"


# ═══════════════════════════════════════════════════════════════════════════════
# MEMORY-EFFICIENT LOADING & CONCATENATION
# ═══════════════════════════════════════════════════════════════════════════════

def load_subset_concat(accessions, filter_fn, obs_cols=None):
    """Load cells matching filter_fn from all datasets, concatenate efficiently.

    Args:
        accessions: list of dataset accession IDs
        filter_fn: callable(obs_df) -> boolean mask
        obs_cols: obs columns to keep (None = all in OBS_COLS_KEEP)

    Returns:
        AnnData with .X (log-norm), .layers['counts'], selected obs columns.
    """
    if obs_cols is None:
        obs_cols = OBS_COLS_KEEP

    # Phase 1: scan backed to identify matching cells and common genes
    print("  Phase 1: Scanning datasets for matching cells...")
    gene_sets = []
    cell_counts = {}
    for acc in accessions:
        path = PROC_DIR / f"{acc}.h5ad"
        if not path.exists():
            print(f"    WARNING: {path} not found, skipping")
            continue
        adata = sc.read_h5ad(path, backed='r')
        obs_df = adata.obs
        mask = filter_fn(obs_df)
        n_match = mask.sum()
        if n_match > 0:
            cell_counts[acc] = n_match
            gene_sets.append(set(adata.var_names))
        adata.file.close()

    if not cell_counts:
        print("  WARNING: No matching cells found!")
        return None

    # Common genes across datasets that have matching cells
    common_genes = sorted(set.intersection(*gene_sets))
    total_cells = sum(cell_counts.values())
    print(f"    Found {total_cells:,} matching cells across {len(cell_counts)} datasets")
    print(f"    Common genes: {len(common_genes):,}")

    # Phase 2: sequential load, subset, collect
    print("  Phase 2: Loading and subsetting...")
    adatas = []
    for acc in accessions:
        if acc not in cell_counts:
            continue
        path = PROC_DIR / f"{acc}.h5ad"
        adata = sc.read_h5ad(path)
        mask = filter_fn(adata.obs)
        adata_sub = adata[mask, common_genes].copy()

        # Keep only selected obs columns
        cols_available = [c for c in obs_cols if c in adata_sub.obs.columns]
        adata_sub.obs = adata_sub.obs[cols_available].copy()

        # Ensure counts layer exists
        if 'counts' in adata_sub.layers:
            pass  # already there
        else:
            print(f"    WARNING: {acc} missing 'counts' layer, using .X")
            adata_sub.layers['counts'] = adata_sub.X.copy()

        adatas.append(adata_sub)
        print(f"    {acc}: {adata_sub.shape[0]:,} cells")
        del adata

    # Phase 3: concatenate
    print("  Phase 3: Concatenating...")
    adata = ad.concat(adatas, merge='same')
    del adatas

    # Ensure sparse
    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(adata.X)
    if not sparse.issparse(adata.layers['counts']):
        adata.layers['counts'] = sparse.csr_matrix(adata.layers['counts'])

    print(f"    Final shape: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
    return adata


# ═══════════════════════════════════════════════════════════════════════════════
# PREPARE INTEGRATED OBJECT (HVG, PCA, NEIGHBORS, UMAP)
# ═══════════════════════════════════════════════════════════════════════════════

def prepare_integrated_object(adata, n_top_genes=3000):
    """Re-derive HVGs, PCA, neighbors, and UMAP on concatenated object.

    Stores unintegrated baseline in obsm['X_pca_unintegrated'] and
    obsm['X_umap_unintegrated'].
    """
    print("  Preparing integrated object (HVG, PCA, neighbors, UMAP)...")

    # HVG selection using counts layer (seurat_v3 requires raw counts)
    n_top = min(n_top_genes, adata.shape[1] - 1)
    n_studies = adata.obs['study'].nunique()

    # Filter out tiny batches for HVG selection (seurat_v3 loess fails on <50 cells)
    MIN_CELLS_HVG_BATCH = 50
    study_counts = adata.obs['study'].value_counts()
    large_studies = study_counts[study_counts >= MIN_CELLS_HVG_BATCH].index.tolist()

    hvg_done = False
    if n_studies > 1 and len(large_studies) > 1:
        try:
            adata_hvg = adata[adata.obs['study'].isin(large_studies)].copy()
            sc.pp.highly_variable_genes(
                adata_hvg, n_top_genes=n_top, flavor='seurat_v3',
                batch_key='study', layer='counts'
            )
            # Transfer HVG annotations to full object
            adata.var['highly_variable'] = False
            adata.var.loc[adata_hvg.var_names[adata_hvg.var['highly_variable']], 'highly_variable'] = True
            # Copy other HVG columns
            for col in ['highly_variable_rank', 'highly_variable_nbatches',
                        'highly_variable_intersection', 'means', 'variances',
                        'variances_norm']:
                if col in adata_hvg.var.columns:
                    adata.var[col] = adata_hvg.var[col]
            del adata_hvg
            hvg_done = True
            n_excluded = n_studies - len(large_studies)
            if n_excluded > 0:
                print(f"    HVG: excluded {n_excluded} studies with <{MIN_CELLS_HVG_BATCH} cells")
        except Exception as e:
            print(f"    WARNING: batch-aware HVG failed ({e}), falling back to non-batch")

    if not hvg_done:
        sc.pp.highly_variable_genes(
            adata, n_top_genes=n_top, flavor='seurat_v3', layer='counts'
        )
    n_hvg = adata.var['highly_variable'].sum()
    print(f"    HVGs: {n_hvg}")

    # PCA on HVG subset
    adata_pca = adata[:, adata.var['highly_variable']].copy()
    gc.collect()
    sc.pp.scale(adata_pca, max_value=10)
    n_comps = min(50, min(adata_pca.shape) - 1)
    sc.tl.pca(adata_pca, n_comps=n_comps, svd_solver='arpack')

    # Transfer PCA back
    adata.obsm['X_pca'] = adata_pca.obsm['X_pca']
    adata.uns['pca'] = adata_pca.uns['pca']
    adata.varm['PCs'] = np.zeros((adata.shape[1], n_comps))
    hvg_idx = np.where(adata.var['highly_variable'])[0]
    adata.varm['PCs'][hvg_idx] = adata_pca.varm['PCs']
    del adata_pca

    # Determine effective dimensionality
    var_ratio = adata.uns['pca']['variance_ratio']
    cumvar = np.cumsum(var_ratio)
    n_pcs = int(np.searchsorted(cumvar, 0.90) + 1)
    n_pcs = max(n_pcs, 10)
    n_pcs = min(n_pcs, n_comps)
    adata.uns['n_pcs_used'] = n_pcs
    print(f"    PCA: {n_comps} components, using {n_pcs} (90% variance)")

    # Unintegrated baseline
    sc.pp.neighbors(adata, n_pcs=n_pcs)
    sc.tl.umap(adata)
    adata.obsm['X_pca_unintegrated'] = adata.obsm['X_pca'].copy()
    adata.obsm['X_umap_unintegrated'] = adata.obsm['X_umap'].copy()
    print("    Unintegrated baseline computed")

    return adata


# ═══════════════════════════════════════════════════════════════════════════════
# APPROACH A: scVI
# ═══════════════════════════════════════════════════════════════════════════════

def run_scvi(adata, batch_key='study', n_latent=20, max_epochs=200,
             use_covariates=True, model_dir=None, label="scvi"):
    """Run scVI integration (Approach A).

    Returns the trained scVI model for potential scANVI initialization.
    Stores obsm['X_scvi'] and obsm['X_umap_scvi'].
    """
    import scvi as scvi_module

    print(f"  Running scVI (Approach A): batch_key={batch_key}, n_latent={n_latent}")

    # Setup anndata for scVI
    scvi_module.model.SCVI.setup_anndata(
        adata,
        layer='counts',
        batch_key=batch_key,
        categorical_covariate_keys=['compartment'] if use_covariates and 'compartment' in adata.obs.columns else None,
        continuous_covariate_keys=['pct_counts_mt'] if use_covariates and 'pct_counts_mt' in adata.obs.columns else None,
    )

    model = scvi_module.model.SCVI(
        adata,
        n_latent=n_latent,
        dispersion='gene-batch',
        gene_likelihood='nb',
    )

    model.train(
        max_epochs=max_epochs,
        early_stopping=True,
        early_stopping_patience=10,
        early_stopping_monitor='elbo_validation',
        train_size=0.9,
        batch_size=256,
    )

    # Extract latent representation
    adata.obsm['X_scvi'] = model.get_latent_representation()
    print(f"    scVI training complete: {model.history['elbo_train'].shape[0]} epochs")

    # Compute neighbors and UMAP on scVI latent space
    sc.pp.neighbors(adata, use_rep='X_scvi')
    sc.tl.umap(adata)
    adata.obsm['X_umap_scvi'] = adata.obsm['X_umap'].copy()

    # Save model
    if model_dir is not None:
        save_path = model_dir / f"scvi_{label}"
        model.save(str(save_path), overwrite=True)
        print(f"    Model saved: {save_path}")

    return model


# ═══════════════════════════════════════════════════════════════════════════════
# APPROACH B: scANVI (semi-supervised)
# ═══════════════════════════════════════════════════════════════════════════════

def run_scanvi(adata, scvi_model, max_epochs=50, model_dir=None, label="scanvi"):
    """Run scANVI integration (Approach B), initialized from scVI model.

    Uses high+medium confidence labels as seeds; low confidence → "Unknown".
    Stores obsm['X_scanvi'], obsm['X_umap_scanvi'], obs['cell_type_scanvi_predicted'].
    """
    import scvi as scvi_module

    print("  Running scANVI (Approach B): semi-supervised from scVI")

    # Create label column: use high+medium confidence, mask low as Unknown
    labels_col = '_scanvi_labels'
    labels = adata.obs['cell_type_final'].astype(str).copy()
    if 'cell_type_confidence' in adata.obs.columns:
        low_mask = adata.obs['cell_type_confidence'] == 'low'
        labels[low_mask] = 'Unknown'
        n_seed = (~low_mask).sum()
        print(f"    Seed labels: {n_seed:,} / {adata.shape[0]:,} cells "
              f"({n_seed / adata.shape[0] * 100:.1f}%)")
    else:
        print("    No confidence column; using all labels as seeds")
    adata.obs[labels_col] = labels

    scanvi_model = scvi_module.model.SCANVI.from_scvi_model(
        scvi_model,
        unlabeled_category='Unknown',
        labels_key=labels_col,
    )

    scanvi_model.train(
        max_epochs=max_epochs,
        early_stopping=True,
        early_stopping_patience=5,
        train_size=0.9,
        batch_size=256,
    )

    # Extract latent representation and predictions
    adata.obsm['X_scanvi'] = scanvi_model.get_latent_representation()
    adata.obs['cell_type_scanvi_predicted'] = scanvi_model.predict()
    print(f"    scANVI training complete: {scanvi_model.history['elbo_train'].shape[0]} epochs")

    # Compute neighbors and UMAP on scANVI latent space
    sc.pp.neighbors(adata, use_rep='X_scanvi')
    sc.tl.umap(adata)
    adata.obsm['X_umap_scanvi'] = adata.obsm['X_umap'].copy()

    # Clean up temp column
    del adata.obs[labels_col]

    # Save model
    if model_dir is not None:
        save_path = model_dir / f"scanvi_{label}"
        scanvi_model.save(str(save_path), overwrite=True)
        print(f"    Model saved: {save_path}")

    return scanvi_model


# ═══════════════════════════════════════════════════════════════════════════════
# APPROACH C: Harmony
# ═══════════════════════════════════════════════════════════════════════════════

def run_harmony(adata, batch_key='study', thetas=(0.5, 1.0, 2.0)):
    """Run Harmony integration (Approach C).

    Tests multiple theta values. Selects the best by quick iLISI check.
    Stores obsm['X_harmony'] and obsm['X_umap_harmony'].
    """
    from harmonypy import run_harmony

    print(f"  Running Harmony (Approach C): testing theta={list(thetas)}")

    n_pcs = adata.uns.get('n_pcs_used', 20)
    pca_input = adata.obsm['X_pca'][:, :n_pcs]
    batch_labels = adata.obs[batch_key].values

    best_theta = None
    best_score = -np.inf
    best_embedding = None

    for theta in thetas:
        print(f"    theta={theta}...", end=" ")
        harmony_out = run_harmony(
            pca_input, adata.obs, batch_key,
            theta=theta, max_iter_harmony=30,
        )
        embedding = harmony_out.Z_corr  # (n_cells, n_pcs) — harmonypy >= 0.0.9

        # Quick iLISI estimate on a subsample
        score = _quick_ilisi(embedding, batch_labels, n_sample=5000)
        print(f"iLISI={score:.3f}")

        if score > best_score:
            best_score = score
            best_theta = theta
            best_embedding = embedding

    print(f"    Best theta: {best_theta} (iLISI={best_score:.3f})")

    adata.obsm['X_harmony'] = best_embedding
    adata.uns['harmony_theta'] = best_theta

    # Compute neighbors and UMAP on Harmony embedding
    sc.pp.neighbors(adata, use_rep='X_harmony')
    sc.tl.umap(adata)
    adata.obsm['X_umap_harmony'] = adata.obsm['X_umap'].copy()

    return best_theta


def _quick_ilisi(embedding, batch_labels, n_sample=5000, k=30):
    """Quick iLISI estimate on a random subsample."""
    from sklearn.neighbors import NearestNeighbors

    n = embedding.shape[0]
    if n > n_sample:
        rng = np.random.RandomState(42)
        idx = rng.choice(n, n_sample, replace=False)
        emb_sub = embedding[idx]
        labels_sub = batch_labels[idx] if isinstance(batch_labels, np.ndarray) else np.array(batch_labels)[idx]
    else:
        emb_sub = embedding
        labels_sub = batch_labels if isinstance(batch_labels, np.ndarray) else np.array(batch_labels)

    nn = NearestNeighbors(n_neighbors=k, algorithm='auto')
    nn.fit(emb_sub)
    indices = nn.kneighbors(emb_sub, return_distance=False)

    # Compute inverse Simpson's index per cell
    ilisi_scores = []
    for i in range(len(emb_sub)):
        neighbor_labels = labels_sub[indices[i]]
        counts = Counter(neighbor_labels)
        total = sum(counts.values())
        p_sq = sum((c / total) ** 2 for c in counts.values())
        ilisi_scores.append(1.0 / p_sq)

    return np.mean(ilisi_scores)


# ═══════════════════════════════════════════════════════════════════════════════
# APPROACH D: BBKNN
# ═══════════════════════════════════════════════════════════════════════════════

def run_bbknn(adata, batch_key='study', n_pcs=20, neighbors_within_batch=3):
    """Run BBKNN integration (Approach D).

    Modifies the neighbor graph directly. No embedding stored — only UMAP.
    Stores obsm['X_umap_bbknn'].
    """
    import bbknn

    print(f"  Running BBKNN (Approach D): n_pcs={n_pcs}, neighbors_within_batch={neighbors_within_batch}")

    n_pcs_actual = min(n_pcs, adata.obsm['X_pca'].shape[1])

    bbknn.bbknn(
        adata,
        batch_key=batch_key,
        n_pcs=n_pcs_actual,
        neighbors_within_batch=neighbors_within_batch,
    )

    sc.tl.umap(adata)
    adata.obsm['X_umap_bbknn'] = adata.obsm['X_umap'].copy()
    print("    BBKNN complete")


# ═══════════════════════════════════════════════════════════════════════════════
# INTEGRATION METRICS
# ═══════════════════════════════════════════════════════════════════════════════

def compute_metrics(adata, embedding_key, batch_key='study', label_key='cell_type_final',
                    max_cells=30_000, random_state=42):
    """Compute integration quality metrics using scib-metrics.

    When the dataset exceeds *max_cells*, a stratified random subsample is
    drawn (preserving the proportions of *label_key*) so that metric
    computation stays within memory on 16 GB machines.  Embeddings and
    integration are still computed on the full dataset — only evaluation is
    subsampled, consistent with the scIB benchmark (Luecken et al. 2022).

    Returns dict with iLISI, cLISI, batch_ASW, celltype_ASW, isolated_label_F1,
    scib_overall, and n_cells_metrics (number of cells used for evaluation).
    """
    import gc

    print(f"    Computing metrics for {embedding_key}...")

    metrics = {}

    # Get embedding
    if embedding_key in adata.obsm:
        embedding = adata.obsm[embedding_key]
    else:
        print(f"    WARNING: {embedding_key} not in obsm, skipping metrics")
        return metrics

    batch = adata.obs[batch_key].values
    labels = adata.obs[label_key].values

    # --- Subsample if needed (stratified by cell type) ---
    n_cells = embedding.shape[0]
    if n_cells > max_cells:
        rng = np.random.RandomState(random_state)
        # Stratified sampling: keep proportions of each label
        unique_labels, label_indices = np.unique(labels, return_inverse=True)
        idx = np.arange(n_cells)
        sampled = []
        for li in range(len(unique_labels)):
            mask = label_indices == li
            group_idx = idx[mask]
            # Always keep at least all members of tiny groups
            n_take = max(1, int(round(max_cells * mask.sum() / n_cells)))
            if n_take >= len(group_idx):
                sampled.append(group_idx)
            else:
                sampled.append(rng.choice(group_idx, size=n_take, replace=False))
        sampled = np.concatenate(sampled)
        rng.shuffle(sampled)
        embedding = embedding[sampled]
        batch = batch[sampled]
        labels = labels[sampled]
        print(f"    Subsampled {n_cells:,} → {len(sampled):,} cells for metric evaluation (stratified, seed={random_state})")
    metrics['n_cells_metrics'] = int(embedding.shape[0])

    try:
        from scib_metrics import ilisi_knn, clisi_knn, silhouette_batch, silhouette_label, isolated_labels
        from scib_metrics.nearest_neighbors import pynndescent
    except ImportError:
        print("    WARNING: scib_metrics not available, using manual metrics")
        return _manual_metrics(embedding, batch, labels)

    # scib-metrics LISI functions expect a NeighborsResults object, not raw arrays
    knn_result = None
    try:
        knn_result = pynndescent(embedding, n_neighbors=90)
    except Exception as e:
        print(f"    WARNING: pynndescent kNN failed: {e}")

    try:
        # iLISI — batch mixing (higher = better mixing)
        if knn_result is not None:
            ilisi = ilisi_knn(knn_result, batch)
            metrics['iLISI'] = float(np.median(ilisi))
            del ilisi
        else:
            metrics['iLISI'] = np.nan
    except Exception as e:
        print(f"    WARNING: iLISI failed: {e}")
        metrics['iLISI'] = np.nan

    try:
        # cLISI — cell type separation (lower = better separation, closer to 1)
        if knn_result is not None:
            clisi = clisi_knn(knn_result, labels)
            metrics['cLISI'] = float(np.median(clisi))
            del clisi
        else:
            metrics['cLISI'] = np.nan
    except Exception as e:
        print(f"    WARNING: cLISI failed: {e}")
        metrics['cLISI'] = np.nan

    # Free kNN graph before silhouette computation
    del knn_result
    gc.collect()

    try:
        # Batch ASW — should be close to 0 (no batch structure)
        metrics['batch_ASW'] = float(silhouette_batch(embedding, labels, batch))
    except Exception as e:
        print(f"    WARNING: batch_ASW failed: {e}")
        metrics['batch_ASW'] = np.nan
    gc.collect()

    try:
        # Cell type ASW — should be positive (cell types separate)
        metrics['celltype_ASW'] = float(silhouette_label(embedding, labels))
    except Exception as e:
        print(f"    WARNING: celltype_ASW failed: {e}")
        metrics['celltype_ASW'] = np.nan
    gc.collect()

    try:
        # Isolated label F1
        metrics['isolated_label_F1'] = float(isolated_labels(embedding, labels, batch))
    except Exception as e:
        print(f"    WARNING: isolated_label_F1 failed: {e}")
        metrics['isolated_label_F1'] = np.nan
    gc.collect()

    # scib overall: average of normalized scores
    batch_scores = []
    bio_scores = []
    if not np.isnan(metrics.get('iLISI', np.nan)):
        # Normalize iLISI: 1 = no mixing, n_batches = perfect; scale to 0-1
        n_batches = len(set(batch))
        batch_scores.append((metrics['iLISI'] - 1) / max(n_batches - 1, 1))
    if not np.isnan(metrics.get('batch_ASW', np.nan)):
        batch_scores.append(metrics['batch_ASW'])
    if not np.isnan(metrics.get('celltype_ASW', np.nan)):
        bio_scores.append(metrics['celltype_ASW'])
    if not np.isnan(metrics.get('cLISI', np.nan)):
        # Normalize cLISI: 1 = perfect, n_types = worst; invert
        n_types = len(set(labels))
        bio_scores.append(1 - (metrics['cLISI'] - 1) / max(n_types - 1, 1))
    if not np.isnan(metrics.get('isolated_label_F1', np.nan)):
        bio_scores.append(metrics['isolated_label_F1'])

    batch_mean = np.mean(batch_scores) if batch_scores else 0.0
    bio_mean = np.mean(bio_scores) if bio_scores else 0.0
    metrics['scib_overall'] = 0.4 * batch_mean + 0.6 * bio_mean

    return metrics


def _manual_metrics(embedding, batch, labels):
    """Fallback manual metrics using sklearn if scib-metrics fails."""
    from sklearn.metrics import silhouette_score
    from sklearn.neighbors import NearestNeighbors

    metrics = {}
    try:
        # Quick iLISI
        metrics['iLISI'] = _quick_ilisi(embedding, batch, n_sample=5000)
    except Exception:
        metrics['iLISI'] = np.nan

    try:
        # Cell type ASW
        if len(set(labels)) > 1:
            metrics['celltype_ASW'] = float(silhouette_score(
                embedding[:5000] if len(embedding) > 5000 else embedding,
                labels[:5000] if len(labels) > 5000 else labels,
                sample_size=min(5000, len(embedding)),
            ))
        else:
            metrics['celltype_ASW'] = np.nan
    except Exception:
        metrics['celltype_ASW'] = np.nan

    metrics['batch_ASW'] = np.nan
    metrics['cLISI'] = np.nan
    metrics['isolated_label_F1'] = np.nan
    metrics['scib_overall'] = np.nan
    return metrics


def compute_continuum_metrics(adata, embedding_key, score_cols, label_key='cell_type_final'):
    """Compute continuum-preservation metrics for a given embedding.

    Returns dict with:
        score_variance_ratio: ratio of score variance in integrated vs unintegrated
        n_clusters_05: number of clusters at resolution 0.5 on integrated embedding
        condition_accuracy: logistic regression accuracy on integrated embedding (5-fold CV)
    """
    metrics = {}

    # Score variance ratio: how much cell state score variance is preserved
    if score_cols:
        available_scores = [c for c in score_cols if c in adata.obs.columns]
        if available_scores:
            # Variance in integrated embedding's neighbor-smoothed space
            # Simpler: compare raw score variance (which is unchanged) — use as baseline
            # Instead compare cluster-level score variance
            unint_var = 0
            int_var = 0
            for col in available_scores:
                scores = adata.obs[col].values
                unint_var += np.var(scores)
                # Smoothed variance: average within Leiden clusters
                if f'leiden_integrated' in adata.obs.columns:
                    cluster_means = adata.obs.groupby('leiden_integrated')[col].mean()
                    int_var += np.var(cluster_means)
                else:
                    int_var += np.var(scores)  # fallback
            metrics['score_variance_ratio'] = int_var / max(unint_var, 1e-10)

    # Cluster count at resolution 0.5
    try:
        if embedding_key in adata.obsm:
            # Temporarily compute neighbors on this embedding
            adata_tmp = adata.copy()
            if embedding_key.startswith('X_umap'):
                # For BBKNN: use the existing neighbor graph (already modified)
                pass
            else:
                sc.pp.neighbors(adata_tmp, use_rep=embedding_key)
            sc.tl.leiden(adata_tmp, resolution=0.5, key_added='_temp_leiden')
            n_clusters = adata_tmp.obs['_temp_leiden'].nunique()
            metrics['n_clusters_05'] = n_clusters
            adata.obs['leiden_integrated'] = adata_tmp.obs['_temp_leiden'].values
            del adata_tmp
        else:
            metrics['n_clusters_05'] = np.nan
    except Exception as e:
        print(f"    WARNING: Cluster count failed: {e}")
        metrics['n_clusters_05'] = np.nan

    # Condition classifier accuracy
    try:
        if 'condition_harmonized' in adata.obs.columns and embedding_key in adata.obsm:
            cond = adata.obs['condition_harmonized'].values
            valid = pd.notna(cond)
            unique_cond = set(cond[valid])
            if len(unique_cond) >= 2:
                from sklearn.linear_model import LogisticRegression
                from sklearn.model_selection import cross_val_score
                X = adata.obsm[embedding_key][valid]
                y = cond[valid]
                # Subsample for speed
                n_sub = min(10000, len(X))
                rng = np.random.RandomState(42)
                idx = rng.choice(len(X), n_sub, replace=False)
                clf = LogisticRegression(max_iter=500, random_state=42)
                scores = cross_val_score(clf, X[idx], y[idx], cv=5, scoring='accuracy')
                metrics['condition_accuracy'] = float(np.mean(scores))
            else:
                metrics['condition_accuracy'] = np.nan
        else:
            metrics['condition_accuracy'] = np.nan
    except Exception as e:
        print(f"    WARNING: Condition classifier failed: {e}")
        metrics['condition_accuracy'] = np.nan

    return metrics


# ═══════════════════════════════════════════════════════════════════════════════
# TIER 1: NON-RESIDENT INTEGRATION
# ═══════════════════════════════════════════════════════════════════════════════

def process_tier1(force=False):
    """Integrate non-resident cells (immune, endothelial, pericyte) with scVI."""
    output_path = INT_DIR / "tier1_nonresident.h5ad"

    if output_path.exists() and not force:
        print("\n=== Tier 1: Non-resident integration (already exists, skipping) ===")
        print(f"    Use --force to re-run")
        return output_path

    print("\n" + "=" * 60)
    print("Tier 1: Non-resident cell integration (scVI)")
    print("=" * 60)

    # Load non-resident cells
    def filter_nonresident(obs_df):
        return obs_df['cell_type_final'].map(classify_cell_tier) == 'nonresident'

    adata = load_subset_concat(ALL_ACCESSIONS, filter_nonresident)
    if adata is None or adata.shape[0] < 100:
        print("  ERROR: Too few non-resident cells for integration")
        return None

    # Prepare (HVG, PCA, baseline UMAP)
    adata = prepare_integrated_object(adata, n_top_genes=3000)

    # Run scVI
    scvi_model = run_scvi(
        adata, batch_key='study', n_latent=20, max_epochs=200,
        use_covariates=False,  # non-resident don't need compartment covariate
        model_dir=MODEL_DIR, label="tier1_nonresident"
    )
    del scvi_model

    # Compute Leiden clusters on scVI embedding
    sc.pp.neighbors(adata, use_rep='X_scvi')
    sc.tl.leiden(adata, resolution=0.5, key_added='leiden_scvi_05')
    sc.tl.leiden(adata, resolution=1.0, key_added='leiden_scvi_10')

    # Compute metrics
    metrics = compute_metrics(adata, 'X_scvi')
    print(f"  Metrics: {metrics}")

    # Generate UMAP figure
    _plot_tier1_umaps(adata)

    # Save
    adata.write_h5ad(output_path)
    print(f"  Saved: {output_path}")

    # Return metrics for summary
    metrics['tier'] = 'tier1'
    metrics['compartment'] = 'nonresident'
    metrics['approach'] = 'scVI'
    return output_path, metrics


def _plot_tier1_umaps(adata):
    """Generate UMAP plots for Tier 1 non-resident integration."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    for ax, (color_key, title) in zip(axes, [
        ('study', 'By study'),
        ('cell_type_final', 'By cell type'),
        ('condition_harmonized', 'By condition'),
    ]):
        if color_key in adata.obs.columns:
            sc.pl.umap(adata, color=color_key, ax=ax, show=False,
                       title=f"Tier 1 — {title}", frameon=False, s=5)
        else:
            ax.set_title(f"Tier 1 — {title} (N/A)")

    plt.tight_layout()
    fig.savefig(RESULTS_DIR / "umap_tier1_nonresident.png", dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {RESULTS_DIR / 'umap_tier1_nonresident.png'}")


# ═══════════════════════════════════════════════════════════════════════════════
# TIER 2: RESIDENT INTEGRATION (PER COMPARTMENT)
# ═══════════════════════════════════════════════════════════════════════════════

def _approach_key_map():
    """Map approach letter to the obsm key that marks it as complete."""
    return {
        'A': 'X_scvi',
        'B': 'X_scanvi',
        'C': 'X_harmony',
        'D': 'X_umap_bbknn',
    }


def _save_checkpoint(adata, output_path, approach_label):
    """Save adata to disk after completing an approach (resumable checkpoint).

    Writes to a temp file first, then atomically renames, so a crash
    mid-write cannot corrupt the existing checkpoint.
    """
    import anndata
    anndata.settings.allow_write_nullable_strings = True
    tmp_path = Path(str(output_path) + '.tmp')
    adata.write_h5ad(tmp_path)
    tmp_path.rename(output_path)
    print(f"  Checkpoint saved after {approach_label}: {output_path}")
    sys.stdout.flush()


def process_tier2_compartment(compartment, approaches=('A', 'B', 'C', 'D'), force=False):
    """Integrate resident cells for one compartment across 4 approaches.

    Saves a checkpoint after each approach completes so the run can be
    resumed if interrupted.  On resume, existing approaches found in the
    h5ad are skipped (unless --force is passed).

    Args:
        compartment: "NP" or "AF"
        approaches: tuple of approach letters to run
        force: re-run even if output exists

    Returns:
        (output_path, list_of_metric_dicts)
    """
    output_path = INT_DIR / f"tier2_resident_{compartment}.h5ad"
    key_map = _approach_key_map()

    # Determine which approaches still need running
    existing_keys = set()
    if output_path.exists() and not force:
        adata_check = sc.read_h5ad(output_path, backed='r')
        existing_keys = set(adata_check.obsm.keys())
        adata_check.file.close()

    needed = [a for a in approaches if force or key_map.get(a, '') not in existing_keys]
    if not needed:
        print(f"\n=== Tier 2 {compartment}: all approaches already exist, skipping ===")
        print(f"    Use --force to re-run")
        return output_path, []

    skipped = [a for a in approaches if a not in needed]
    if skipped:
        print(f"\n  Resuming Tier 2 {compartment}: skipping completed approaches {skipped}")

    print(f"\n{'=' * 60}")
    print(f"Tier 2: {compartment} resident cell integration")
    print(f"  Approaches to run: {needed}")
    print(f"{'=' * 60}")
    sys.stdout.flush()

    # Define filter function
    if compartment == "NP":
        def filter_fn(obs_df):
            tiers = obs_df['cell_type_final'].map(classify_cell_tier)
            return (tiers == 'NP') | (tiers == 'EP')  # EP included with NP
        score_cols = NP_SCORE_COLS
    elif compartment == "AF":
        def filter_fn(obs_df):
            tiers = obs_df['cell_type_final'].map(classify_cell_tier)
            return (tiers == 'AF') | (tiers == 'fibroblast')  # fibroblasts with AF
        score_cols = AF_SCORE_COLS
    else:
        raise ValueError(f"Unknown compartment: {compartment}")

    # Load from checkpoint if available, otherwise from raw data
    if output_path.exists() and not force and existing_keys:
        print("  Loading from checkpoint (sparse-backed to save memory)...")
        sys.stdout.flush()
        # Load X as sparse to keep peak memory manageable on 16 GB machines.
        # The dense X is only needed by scVI/scANVI setup; approaches C/D and
        # metric evaluation work fine with sparse or convert small slices.
        adata = sc.read_h5ad(output_path)
        # Convert X to sparse if it is dense — saves ~4 GB for 139K×11K
        import scipy.sparse as sp
        if not sp.issparse(adata.X):
            print("    Converting X to sparse CSR to reduce memory footprint...")
            adata.X = sp.csr_matrix(adata.X)
            gc.collect()
        print(f"    Loaded: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
        print(f"    Existing obsm: {list(adata.obsm.keys())}")
        sys.stdout.flush()
    else:
        adata = load_subset_concat(ALL_ACCESSIONS, filter_fn)
        if adata is None or adata.shape[0] < 100:
            print(f"  ERROR: Too few {compartment} cells for integration")
            return None, []

        # Flag EP cells
        if compartment == "NP":
            adata.obs['is_ep'] = adata.obs['cell_type_final'].map(
                lambda x: classify_cell_tier(x) == 'EP'
            )
            n_ep = adata.obs['is_ep'].sum()
            print(f"  EP cells included: {n_ep:,}")

        # Prepare (HVG, PCA, baseline UMAP)
        adata = prepare_integrated_object(adata, n_top_genes=3000)

    # Run approaches
    all_metrics = []
    scvi_model = None

    # Approach A: scVI
    if 'A' in needed:
        try:
            scvi_model = run_scvi(
                adata, batch_key='study', n_latent=20, max_epochs=200,
                use_covariates=True, model_dir=MODEL_DIR,
                label=f"tier2_{compartment}"
            )
            m = compute_metrics(adata, 'X_scvi')
            m.update(compute_continuum_metrics(adata, 'X_scvi', score_cols))
            m['tier'] = 'tier2'
            m['compartment'] = compartment
            m['approach'] = 'A_scVI'
            all_metrics.append(m)
            print(f"    Approach A metrics: {m}")
            _save_checkpoint(adata, output_path, 'A_scVI')
        except Exception as e:
            print(f"  ERROR in Approach A (scVI): {e}")
            traceback.print_exc()

    # Approach B: scANVI
    if 'B' in needed:
        if scvi_model is None:
            # Try to reload saved scVI model
            saved_model_path = MODEL_DIR / f"scvi_tier2_{compartment}"
            if saved_model_path.exists():
                print("  Loading saved scVI model for scANVI initialization...")
                import scvi as scvi_module
                scvi_module.model.SCVI.setup_anndata(
                    adata, layer='counts', batch_key='study',
                    categorical_covariate_keys=['compartment'] if 'compartment' in adata.obs.columns else None,
                    continuous_covariate_keys=['pct_counts_mt'] if 'pct_counts_mt' in adata.obs.columns else None,
                )
                scvi_model = scvi_module.model.SCVI.load(str(saved_model_path), adata=adata)
            else:
                print("  WARNING: scANVI requires scVI model; running scVI first...")
                try:
                    scvi_model = run_scvi(
                        adata, batch_key='study', n_latent=20, max_epochs=200,
                        use_covariates=True, model_dir=MODEL_DIR,
                        label=f"tier2_{compartment}"
                    )
                    _save_checkpoint(adata, output_path, 'A_scVI_for_B')
                except Exception as e:
                    print(f"  ERROR: Could not train scVI for scANVI init: {e}")
        if scvi_model is not None:
            try:
                run_scanvi(
                    adata, scvi_model, max_epochs=50,
                    model_dir=MODEL_DIR, label=f"tier2_{compartment}"
                )
                m = compute_metrics(adata, 'X_scanvi')
                m.update(compute_continuum_metrics(adata, 'X_scanvi', score_cols))
                m['tier'] = 'tier2'
                m['compartment'] = compartment
                m['approach'] = 'B_scANVI'
                all_metrics.append(m)
                print(f"    Approach B metrics: {m}")
                _save_checkpoint(adata, output_path, 'B_scANVI')
            except Exception as e:
                print(f"  ERROR in Approach B (scANVI): {e}")
                traceback.print_exc()

    del scvi_model  # free memory
    gc.collect()

    # Approach C: Harmony
    if 'C' in needed:
        try:
            best_theta = run_harmony(adata, batch_key='study', thetas=(0.5, 1.0, 2.0))
            m = compute_metrics(adata, 'X_harmony')
            m.update(compute_continuum_metrics(adata, 'X_harmony', score_cols))
            m['tier'] = 'tier2'
            m['compartment'] = compartment
            m['approach'] = 'C_Harmony'
            m['harmony_theta'] = best_theta
            all_metrics.append(m)
            print(f"    Approach C metrics: {m}")
            _save_checkpoint(adata, output_path, 'C_Harmony')
        except Exception as e:
            print(f"  ERROR in Approach C (Harmony): {e}")
            traceback.print_exc()

    # Approach D: BBKNN
    if 'D' in needed:
        try:
            run_bbknn(adata, batch_key='study', n_pcs=20, neighbors_within_batch=3)
            # BBKNN has no embedding — use PCA for metrics, UMAP for viz
            m = compute_metrics(adata, 'X_pca')
            m.update(compute_continuum_metrics(adata, 'X_pca', score_cols))
            m['tier'] = 'tier2'
            m['compartment'] = compartment
            m['approach'] = 'D_BBKNN'
            all_metrics.append(m)
            print(f"    Approach D metrics: {m}")
            _save_checkpoint(adata, output_path, 'D_BBKNN')
        except Exception as e:
            print(f"  ERROR in Approach D (BBKNN): {e}")
            traceback.print_exc()

    # Generate UMAP figures
    _plot_tier2_umaps(adata, compartment)

    # Final save (with plots info in uns)
    _save_checkpoint(adata, output_path, 'final')
    print(f"  Saved: {output_path}")

    return output_path, all_metrics


def _plot_tier2_umaps(adata, compartment):
    """Generate multi-panel UMAP comparison for Tier 2 approaches."""
    # Collect available embeddings
    umap_keys = []
    umap_labels = []
    if 'X_umap_unintegrated' in adata.obsm:
        umap_keys.append('X_umap_unintegrated')
        umap_labels.append('Unintegrated')
    if 'X_umap_scvi' in adata.obsm:
        umap_keys.append('X_umap_scvi')
        umap_labels.append('A: scVI')
    if 'X_umap_scanvi' in adata.obsm:
        umap_keys.append('X_umap_scanvi')
        umap_labels.append('B: scANVI')
    if 'X_umap_harmony' in adata.obsm:
        umap_keys.append('X_umap_harmony')
        umap_labels.append('C: Harmony')
    if 'X_umap_bbknn' in adata.obsm:
        umap_keys.append('X_umap_bbknn')
        umap_labels.append('D: BBKNN')

    n_cols = len(umap_keys)
    if n_cols == 0:
        return

    color_keys = ['study', 'cell_type_final', 'condition_harmonized', 'compartment']
    color_keys = [c for c in color_keys if c in adata.obs.columns]
    n_rows = len(color_keys)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4.5 * n_rows))
    if n_rows == 1:
        axes = axes[np.newaxis, :]
    if n_cols == 1:
        axes = axes[:, np.newaxis]

    for col_i, (umap_key, umap_label) in enumerate(zip(umap_keys, umap_labels)):
        # Temporarily set X_umap for sc.pl.umap
        adata.obsm['X_umap'] = adata.obsm[umap_key]
        for row_i, color_key in enumerate(color_keys):
            ax = axes[row_i, col_i]
            try:
                sc.pl.umap(adata, color=color_key, ax=ax, show=False,
                           frameon=False, s=2, alpha=0.5)
                if row_i == 0:
                    ax.set_title(umap_label, fontsize=11, fontweight='bold')
                if col_i == 0:
                    ax.set_ylabel(color_key.replace('_', ' '), fontsize=10)
            except Exception:
                ax.set_title(f"{umap_label}\n({color_key})")

    plt.tight_layout()
    fig.savefig(RESULTS_DIR / f"umap_tier2_{compartment}_by_approach.png",
                dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {RESULTS_DIR / f'umap_tier2_{compartment}_by_approach.png'}")


# ═══════════════════════════════════════════════════════════════════════════════
# REPORT GENERATION
# ═══════════════════════════════════════════════════════════════════════════════

INTEGRATION_REPORT_TEMPLATE = """<!DOCTYPE html>
<html>
<head>
<title>Integration Report — Module 05</title>
<style>
  body { font-family: Arial, sans-serif; max-width: 1400px; margin: 0 auto; padding: 20px; }
  h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
  h2 { color: #34495e; margin-top: 30px; }
  h3 { color: #7f8c8d; }
  table { border-collapse: collapse; width: 100%; margin: 10px 0; }
  th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
  th { background-color: #3498db; color: white; }
  tr:nth-child(even) { background-color: #f2f2f2; }
  img { max-width: 100%; height: auto; margin: 10px 0; border: 1px solid #ddd; }
  .stats-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0; }
  .stat-box { background: #ecf0f1; padding: 15px; border-radius: 5px; }
  .stat-box h3 { margin: 0 0 5px 0; color: #7f8c8d; font-size: 12px; }
  .stat-box p { margin: 0; font-size: 24px; font-weight: bold; color: #2c3e50; }
  .warning { color: #e74c3c; font-weight: bold; }
  .pass { color: #27ae60; }
  .info { color: #7f8c8d; }
  .best { background-color: #d5f5e3; }
  .validation li { margin: 4px 0; }
</style>
</head>
<body>
<h1>Module 05: Cross-Dataset Integration Report</h1>
<p>Generated: {{ date }}</p>

<div class="stats-grid">
  <div class="stat-box"><h3>Total cells integrated</h3><p>{{ total_cells }}</p></div>
  <div class="stat-box"><h3>Tier 1 (non-resident)</h3><p>{{ tier1_cells }}</p></div>
  <div class="stat-box"><h3>Tier 2 NP</h3><p>{{ tier2_np_cells }}</p></div>
  <div class="stat-box"><h3>Tier 2 AF</h3><p>{{ tier2_af_cells }}</p></div>
  <div class="stat-box"><h3>Studies</h3><p>{{ n_studies }}</p></div>
  <div class="stat-box"><h3>Approaches tested</h3><p>{{ n_approaches }}</p></div>
</div>

<h2>Tier 1: Non-Resident Cells</h2>
{% if has_tier1_umap %}
<img src="umap_tier1_nonresident.png" alt="Tier 1 UMAP">
{% endif %}

<h2>Tier 2: NP Resident Cells</h2>
{% if has_tier2_np_umap %}
<img src="umap_tier2_NP_by_approach.png" alt="Tier 2 NP UMAP">
{% endif %}

<h2>Tier 2: AF Resident Cells</h2>
{% if has_tier2_af_umap %}
<img src="umap_tier2_AF_by_approach.png" alt="Tier 2 AF UMAP">
{% endif %}

<h2>Integration Metrics</h2>
<table>
  <tr>
    <th>Tier</th><th>Compartment</th><th>Approach</th>
    <th>iLISI</th><th>cLISI</th><th>Batch ASW</th><th>Cell Type ASW</th>
    <th>Isolated Label F1</th><th>Overall Score</th>
    <th>Clusters (res 0.5)</th><th>Condition Accuracy</th>
  </tr>
  {% for row in metrics_table %}
  <tr class="{{ 'best' if row.is_best else '' }}">
    <td>{{ row.tier }}</td>
    <td>{{ row.compartment }}</td>
    <td>{{ row.approach }}</td>
    <td>{{ row.iLISI }}</td>
    <td>{{ row.cLISI }}</td>
    <td>{{ row.batch_ASW }}</td>
    <td>{{ row.celltype_ASW }}</td>
    <td>{{ row.isolated_label_F1 }}</td>
    <td>{{ row.scib_overall }}</td>
    <td>{{ row.n_clusters_05 }}</td>
    <td>{{ row.condition_accuracy }}</td>
  </tr>
  {% endfor %}
</table>

{% if has_metrics_comparison %}
<h2>Metrics Comparison</h2>
<img src="metrics_comparison_NP.png" alt="NP metrics comparison">
<img src="metrics_comparison_AF.png" alt="AF metrics comparison">
{% endif %}

{% if has_cluster_count %}
<h2>Cluster Count Comparison</h2>
<img src="cluster_count_comparison.png" alt="Cluster count comparison">
{% endif %}

<h2>Validation</h2>
<ul class="validation">
  {% for msg in validation_messages %}
  <li class="{{ 'pass' if msg.startswith('PASS') else 'warning' if msg.startswith('FAIL') or msg.startswith('WARNING') else 'info' }}">{{ msg }}</li>
  {% endfor %}
</ul>

<h2>Human Checkpoint Questions</h2>
<ol>
  <li>Which integration approach (A-D) best preserves cell state variation while adequately removing batch effects?</li>
  <li>Is any approach clearly superior, or is a combination needed (e.g., scANVI for NP, Harmony for AF)?</li>
  <li>Does the "blob" problem recur with any approach? If so, is Approach E or F the appropriate fallback?</li>
  <li>Should the analysis proceed with integrated data, per-dataset data, or both in parallel?</li>
  <li>Are there any study-specific effects that persist after integration?</li>
  <li>Does the integration reveal any new cell states not visible in per-dataset analysis?</li>
</ol>

</body>
</html>"""


def generate_integration_report(all_metrics, validation_messages):
    """Generate HTML integration report."""
    print("\n  Generating integration report...")

    # Load cell counts from output files
    tier1_cells = 0
    tier2_np_cells = 0
    tier2_af_cells = 0
    n_studies = 12

    t1_path = INT_DIR / "tier1_nonresident.h5ad"
    if t1_path.exists():
        a = sc.read_h5ad(t1_path, backed='r')
        tier1_cells = a.shape[0]
        a.file.close()
    np_path = INT_DIR / "tier2_resident_NP.h5ad"
    if np_path.exists():
        a = sc.read_h5ad(np_path, backed='r')
        tier2_np_cells = a.shape[0]
        a.file.close()
    af_path = INT_DIR / "tier2_resident_AF.h5ad"
    if af_path.exists():
        a = sc.read_h5ad(af_path, backed='r')
        tier2_af_cells = a.shape[0]
        a.file.close()

    total_cells = tier1_cells + tier2_np_cells + tier2_af_cells

    # Build metrics table rows
    metrics_table = []
    # Find best overall per compartment
    best_by_comp = {}
    for m in all_metrics:
        comp = m.get('compartment', '')
        score = m.get('scib_overall', 0)
        if comp not in best_by_comp or score > best_by_comp[comp]:
            best_by_comp[comp] = score

    for m in all_metrics:
        comp = m.get('compartment', '')
        is_best = (m.get('scib_overall', 0) == best_by_comp.get(comp, -1)
                   and not np.isnan(m.get('scib_overall', np.nan)))
        metrics_table.append({
            'tier': m.get('tier', ''),
            'compartment': comp,
            'approach': m.get('approach', ''),
            'iLISI': f"{m.get('iLISI', np.nan):.3f}" if not np.isnan(m.get('iLISI', np.nan)) else '—',
            'cLISI': f"{m.get('cLISI', np.nan):.3f}" if not np.isnan(m.get('cLISI', np.nan)) else '—',
            'batch_ASW': f"{m.get('batch_ASW', np.nan):.3f}" if not np.isnan(m.get('batch_ASW', np.nan)) else '—',
            'celltype_ASW': f"{m.get('celltype_ASW', np.nan):.3f}" if not np.isnan(m.get('celltype_ASW', np.nan)) else '—',
            'isolated_label_F1': f"{m.get('isolated_label_F1', np.nan):.3f}" if not np.isnan(m.get('isolated_label_F1', np.nan)) else '—',
            'scib_overall': f"{m.get('scib_overall', np.nan):.3f}" if not np.isnan(m.get('scib_overall', np.nan)) else '—',
            'n_clusters_05': str(m.get('n_clusters_05', '—')),
            'condition_accuracy': f"{m.get('condition_accuracy', np.nan):.3f}" if not np.isnan(m.get('condition_accuracy', np.nan)) else '—',
            'is_best': is_best,
        })

    template = Template(INTEGRATION_REPORT_TEMPLATE)
    html = template.render(
        date=datetime.now().strftime("%Y-%m-%d %H:%M"),
        total_cells=f"{total_cells:,}",
        tier1_cells=f"{tier1_cells:,}",
        tier2_np_cells=f"{tier2_np_cells:,}",
        tier2_af_cells=f"{tier2_af_cells:,}",
        n_studies=n_studies,
        n_approaches=len(all_metrics),
        has_tier1_umap=(RESULTS_DIR / "umap_tier1_nonresident.png").exists(),
        has_tier2_np_umap=(RESULTS_DIR / "umap_tier2_NP_by_approach.png").exists(),
        has_tier2_af_umap=(RESULTS_DIR / "umap_tier2_AF_by_approach.png").exists(),
        has_metrics_comparison=(RESULTS_DIR / "metrics_comparison_NP.png").exists(),
        has_cluster_count=(RESULTS_DIR / "cluster_count_comparison.png").exists(),
        metrics_table=metrics_table,
        validation_messages=validation_messages,
    )

    report_path = RESULTS_DIR / "integration_report.html"
    report_path.write_text(html)
    print(f"  Report: {report_path}")


def plot_metrics_comparison(all_metrics):
    """Generate grouped bar charts comparing metrics across approaches."""
    if not all_metrics:
        return

    df = pd.DataFrame(all_metrics)

    for compartment in df['compartment'].unique():
        comp_df = df[df['compartment'] == compartment].copy()
        if len(comp_df) < 2:
            continue

        metric_cols = ['iLISI', 'celltype_ASW', 'batch_ASW', 'isolated_label_F1', 'scib_overall']
        available = [c for c in metric_cols if c in comp_df.columns and comp_df[c].notna().any()]
        if not available:
            continue

        fig, ax = plt.subplots(figsize=(10, 5))
        x = np.arange(len(available))
        width = 0.8 / len(comp_df)

        for i, (_, row) in enumerate(comp_df.iterrows()):
            vals = [row.get(c, np.nan) for c in available]
            ax.bar(x + i * width, vals, width, label=row['approach'], alpha=0.8)

        ax.set_xticks(x + width * (len(comp_df) - 1) / 2)
        ax.set_xticklabels(available, rotation=30, ha='right')
        ax.set_ylabel('Score')
        ax.set_title(f'{compartment}: Integration Metrics Comparison')
        ax.legend()
        plt.tight_layout()
        fig.savefig(RESULTS_DIR / f"metrics_comparison_{compartment}.png",
                    dpi=150, bbox_inches='tight')
        plt.close(fig)

    # Cluster count comparison
    if 'n_clusters_05' in df.columns and df['n_clusters_05'].notna().any():
        fig, ax = plt.subplots(figsize=(8, 5))
        for compartment in df['compartment'].unique():
            comp_df = df[df['compartment'] == compartment]
            ax.bar(
                [f"{row['approach']}\n({compartment})" for _, row in comp_df.iterrows()],
                comp_df['n_clusters_05'].values,
                alpha=0.8,
            )
        ax.set_ylabel('Number of clusters (resolution 0.5)')
        ax.set_title('Cluster Count by Approach')
        plt.tight_layout()
        fig.savefig(RESULTS_DIR / "cluster_count_comparison.png",
                    dpi=150, bbox_inches='tight')
        plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════════
# VALIDATION
# ═══════════════════════════════════════════════════════════════════════════════

def validate_integration():
    """Validate integration outputs. Returns (pass, messages)."""
    print("\n" + "=" * 60)
    print("Validating integration outputs")
    print("=" * 60)

    messages = []
    all_pass = True

    # Check output files exist
    for fname, desc in [
        ("tier1_nonresident.h5ad", "Tier 1 non-resident"),
        ("tier2_resident_NP.h5ad", "Tier 2 NP"),
        ("tier2_resident_AF.h5ad", "Tier 2 AF"),
    ]:
        path = INT_DIR / fname
        if path.exists():
            messages.append(f"PASS: {desc} output exists ({path})")
        else:
            messages.append(f"FAIL: {desc} output missing ({path})")
            all_pass = False

    # Check obsm keys for Tier 2
    for fname, comp in [
        ("tier2_resident_NP.h5ad", "NP"),
        ("tier2_resident_AF.h5ad", "AF"),
    ]:
        path = INT_DIR / fname
        if not path.exists():
            continue
        adata = sc.read_h5ad(path, backed='r')
        obsm_keys = set(adata.obsm.keys())
        for key, approach in [
            ('X_scvi', 'A_scVI'),
            ('X_scanvi', 'B_scANVI'),
            ('X_harmony', 'C_Harmony'),
            ('X_umap_bbknn', 'D_BBKNN'),
        ]:
            if key in obsm_keys:
                messages.append(f"PASS: {comp} has {approach} embedding ({key})")
            else:
                messages.append(f"FAIL: {comp} missing {approach} embedding ({key})")
                all_pass = False
        adata.file.close()

    # Blob check: every approach should produce >1 cluster at res 0.5
    for fname, comp in [
        ("tier2_resident_NP.h5ad", "NP"),
        ("tier2_resident_AF.h5ad", "AF"),
    ]:
        path = INT_DIR / fname
        if not path.exists():
            continue
        adata = sc.read_h5ad(path)
        for key, approach in [
            ('X_scvi', 'A_scVI'),
            ('X_scanvi', 'B_scANVI'),
            ('X_harmony', 'C_Harmony'),
        ]:
            if key not in adata.obsm:
                continue
            try:
                sc.pp.neighbors(adata, use_rep=key)
                sc.tl.leiden(adata, resolution=0.5, key_added='_blob_check')
                n_clusters = adata.obs['_blob_check'].nunique()
                if n_clusters > 1:
                    messages.append(f"PASS: {comp} {approach} blob check: {n_clusters} clusters at res 0.5")
                else:
                    messages.append(f"WARNING: {comp} {approach} blob detected: only {n_clusters} cluster(s) at res 0.5")
            except Exception as e:
                messages.append(f"WARNING: {comp} {approach} blob check failed: {e}")

        # Study→cluster ARI check (ARI < 1.0 means integration worked)
        for key, approach in [
            ('X_scvi', 'A_scVI'),
            ('X_scanvi', 'B_scANVI'),
            ('X_harmony', 'C_Harmony'),
        ]:
            if key not in adata.obsm:
                continue
            try:
                from sklearn.metrics import adjusted_rand_score
                sc.pp.neighbors(adata, use_rep=key)
                sc.tl.leiden(adata, resolution=0.5, key_added='_ari_check')
                ari = adjusted_rand_score(
                    adata.obs['study'].values,
                    adata.obs['_ari_check'].values
                )
                if ari < 1.0:
                    messages.append(f"PASS: {comp} {approach} study-cluster ARI = {ari:.3f} (< 1.0)")
                else:
                    messages.append(f"WARNING: {comp} {approach} study perfectly predicts clusters (ARI = {ari:.3f})")
            except Exception as e:
                messages.append(f"WARNING: {comp} {approach} ARI check failed: {e}")
        del adata

    # Check metrics file
    metrics_path = RESULTS_DIR / "integration_metrics.tsv"
    if metrics_path.exists():
        mdf = pd.read_csv(metrics_path, sep='\t')
        messages.append(f"PASS: Metrics file exists with {len(mdf)} rows")
    else:
        messages.append(f"FAIL: Metrics file missing ({metrics_path})")
        all_pass = False

    # Check report
    report_path = RESULTS_DIR / "integration_report.html"
    if report_path.exists():
        messages.append(f"PASS: Integration report exists")
    else:
        messages.append(f"FAIL: Integration report missing")
        all_pass = False

    # Print results
    for msg in messages:
        status = "PASS" if msg.startswith("PASS") else "FAIL" if msg.startswith("FAIL") else "WARN"
        print(f"  [{status}] {msg}")

    return all_pass, messages


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    args = sys.argv[1:]

    # Parse flags
    validate_only = '--validate-only' in args
    tier1_only = '--tier1-only' in args
    tier2_only = '--tier2-only' in args
    force = '--force' in args

    # Parse --approach
    approach_filter = None
    for i, a in enumerate(args):
        if a == '--approach' and i + 1 < len(args):
            approach_filter = tuple(args[i + 1].upper().split(','))

    # Validation only
    if validate_only:
        passed, messages = validate_integration()
        sys.exit(0 if passed else 1)

    print("=" * 60)
    print("Module 05: Cross-Dataset Integration")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    all_metrics = []

    # Tier 1
    if not tier2_only:
        result = process_tier1(force=force)
        if result is not None and isinstance(result, tuple):
            _, tier1_metrics = result
            all_metrics.append(tier1_metrics)

    # Tier 2
    if not tier1_only:
        approaches = approach_filter if approach_filter else ('A', 'B', 'C', 'D')

        # NP
        _, np_metrics = process_tier2_compartment('NP', approaches=approaches, force=force)
        all_metrics.extend(np_metrics)

        # AF
        _, af_metrics = process_tier2_compartment('AF', approaches=approaches, force=force)
        all_metrics.extend(af_metrics)

    # Save metrics — merge with any existing from previous runs
    metrics_path = RESULTS_DIR / "integration_metrics.tsv"
    if all_metrics:
        new_df = pd.DataFrame(all_metrics)
        if metrics_path.exists():
            old_df = pd.read_csv(metrics_path, sep='\t')
            # Remove old rows that were re-computed in this run
            new_keys = set(zip(new_df['tier'], new_df['compartment'], new_df['approach']))
            old_df = old_df[~old_df.apply(
                lambda r: (r.get('tier'), r.get('compartment'), r.get('approach')) in new_keys,
                axis=1
            )]
            metrics_df = pd.concat([old_df, new_df], ignore_index=True)
        else:
            metrics_df = new_df
        metrics_df.to_csv(metrics_path, sep='\t', index=False)
        print(f"\n  Metrics saved: {metrics_path} ({len(metrics_df)} total rows)")

        # Generate comparison plots
        plot_metrics_comparison(all_metrics)

    # Validation
    passed, val_messages = validate_integration()

    # Generate report
    generate_integration_report(all_metrics, val_messages)

    print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Overall validation: {'PASSED' if passed else 'FAILED'}")


if __name__ == "__main__":
    main()
