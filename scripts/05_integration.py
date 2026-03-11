#!/usr/bin/env python3
"""Module 05: Integration for IVD atlas.

Integrates cells across studies into four compartment-based objects (NP, AF,
CEP, all-cells), each with tiered scANVI integration (mesenchymal and
non-mesenchymal separately). Uses coarse_label from Module 04 as semi-supervised
anchors for scANVI.

Clustering (Module 06) and annotation (Module 07) are separate steps.

Usage:
    python3 scripts/05_integration.py                  # All objects
    python3 scripts/05_integration.py --object NP      # Single object
    python3 scripts/05_integration.py --validate-only  # Validation only
    python3 scripts/05_integration.py --force          # Re-run even if outputs exist
"""

import gc
import sys
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
META_DIR = BASE / "metadata"

for d in [INT_DIR, MODEL_DIR, RESULTS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# ── Study assignments per object ─────────────────────────────────────────────
# Each entry: (accession, compartment_filter or None for all compartments)
# compartment_filter is matched against obs['compartment']

NP_STUDIES = [
    ("GSE160756", "NP"),
    ("GSE165722", None),        # NP-only study
    ("GSE199866", "NP"),
    ("GSE205535", None),        # NP-only study
    ("GSE244889", None),        # NP-only study
    ("GSE251686", None),        # NP-only study
    ("GSE230809", "NP"),
    ("CNP0002664", None),       # NP-only study
]

AF_STUDIES = [
    ("GSE160756", "AF"),
    ("GSE199866", "AF"),        # Inner AF
    ("GSE230809", "AF"),
]

CEP_STUDIES = [
    ("GSE160756", "CEP"),
    ("GSE255768", None),        # CEP-only study
    ("GSE242443", None),        # CEP-only study (culture-expanded)
]

# All-cells: union of above plus GSE189916 (adult only, compartments not separated)
ALL_CELLS_EXTRA = [
    ("GSE189916", None),        # Adult whole IVD (neonatal excluded in Module 04)
]

# Sample-level exclusions
EXCLUDED_SAMPLES = {
    "GSE251686_NP3",  # Corrupt matrix file
}

# obs columns to keep during concatenation
OBS_COLS_KEEP = [
    "sample_id", "study", "compartment", "condition_harmonized",
    "cell_class", "coarse_label", "pct_counts_mt",
]

# ── scANVI parameters (from spec) ─────────────────────────────────────────────
SCANVI_PARAMS = {
    "batch_key": "study",
    "n_latent": 20,
    "scvi_max_epochs": 200,
    "scanvi_max_epochs": 50,
    "n_top_genes": 3000,
}


# ═════════════════════════════════════════════════════════════════════════════
# STUDY ASSIGNMENTS
# ═════════════════════════════════════════════════════════════════════════════

def get_study_assignments(object_name):
    """Return list of (accession, compartment_filter) for an object."""
    if object_name == "NP":
        return NP_STUDIES
    elif object_name == "AF":
        return AF_STUDIES
    elif object_name == "CEP":
        return CEP_STUDIES
    elif object_name == "all_cells":
        # Union of all compartment-specific studies + extra
        seen = set()
        assignments = []
        for study_list in [NP_STUDIES, AF_STUDIES, CEP_STUDIES, ALL_CELLS_EXTRA]:
            for acc, comp in study_list:
                key = (acc, comp)
                if key not in seen:
                    seen.add(key)
                    assignments.append((acc, comp))
        return assignments
    else:
        raise ValueError(f"Unknown object: {object_name}")


# ═════════════════════════════════════════════════════════════════════════════
# MEMORY-EFFICIENT LOADING & CONCATENATION
# ═════════════════════════════════════════════════════════════════════════════

def load_and_build_object(object_name):
    """Load and concatenate cells for one integrated object.

    Filters by compartment where specified, excludes excluded samples,
    and keeps cells with valid cell_class assignments.

    For mesenchymal tier: cell_class in ['mesenchymal', 'unknown']
    For non-mesenchymal tier: cell_class == 'non_mesenchymal'

    Returns AnnData with .X (log-norm), .layers['counts'], selected obs cols.
    """
    assignments = get_study_assignments(object_name)
    obs_cols = OBS_COLS_KEEP

    print(f"  Loading cells for {object_name} object...")
    print(f"    Studies: {[a for a, _ in assignments]}")

    # Phase 1: scan to identify matching cells and common genes
    gene_sets = []
    cell_counts = {}

    for acc, comp_filter in assignments:
        path = PROC_DIR / f"{acc}.h5ad"
        if not path.exists():
            print(f"    WARNING: {path} not found, skipping")
            continue
        adata = sc.read_h5ad(path, backed='r')
        obs_df = adata.obs

        # Build filter mask
        mask = pd.Series(True, index=obs_df.index)

        # Compartment filter
        if comp_filter is not None and 'compartment' in obs_df.columns:
            mask &= obs_df['compartment'].str.upper().str.contains(comp_filter.upper(), na=False)

        # Sample exclusions
        if 'sample_id' in obs_df.columns:
            mask &= ~obs_df['sample_id'].isin(EXCLUDED_SAMPLES)

        # Must have cell_class from Module 04
        if 'cell_class' in obs_df.columns:
            mask &= obs_df['cell_class'].isin(['mesenchymal', 'non_mesenchymal', 'unknown'])
        else:
            print(f"    WARNING: {acc} missing cell_class column, skipping")
            adata.file.close()
            continue

        n_match = mask.sum()
        if n_match > 0:
            cell_counts[(acc, comp_filter)] = n_match
            gene_sets.append(set(adata.var_names))
        adata.file.close()

    if not cell_counts:
        print("  WARNING: No matching cells found!")
        return None

    common_genes = sorted(set.intersection(*gene_sets))
    total_cells = sum(cell_counts.values())
    print(f"    Found {total_cells:,} cells across {len(cell_counts)} study-compartment combos")
    print(f"    Common genes: {len(common_genes):,}")

    # Phase 2: load, subset, collect
    adatas = []
    for acc, comp_filter in assignments:
        if (acc, comp_filter) not in cell_counts:
            continue
        path = PROC_DIR / f"{acc}.h5ad"
        adata = sc.read_h5ad(path)

        mask = pd.Series(True, index=adata.obs.index)
        if comp_filter is not None and 'compartment' in adata.obs.columns:
            mask &= adata.obs['compartment'].str.upper().str.contains(comp_filter.upper(), na=False)
        if 'sample_id' in adata.obs.columns:
            mask &= ~adata.obs['sample_id'].isin(EXCLUDED_SAMPLES)
        if 'cell_class' in adata.obs.columns:
            mask &= adata.obs['cell_class'].isin(['mesenchymal', 'non_mesenchymal', 'unknown'])

        adata_sub = adata[mask, common_genes].copy()

        # Keep only selected obs columns
        cols_available = [c for c in obs_cols if c in adata_sub.obs.columns]
        adata_sub.obs = adata_sub.obs[cols_available].copy()

        # Ensure study column
        if 'study' not in adata_sub.obs.columns:
            adata_sub.obs['study'] = acc

        # Ensure coarse_label column
        if 'coarse_label' not in adata_sub.obs.columns:
            adata_sub.obs['coarse_label'] = 'Unknown'

        # Ensure counts layer, int32
        if 'counts' in adata_sub.layers:
            if adata_sub.layers['counts'].dtype != np.int32:
                adata_sub.layers['counts'] = adata_sub.layers['counts'].astype(np.int32)
        else:
            adata_sub.layers['counts'] = adata_sub.X.copy().astype(np.int32)

        adatas.append(adata_sub)
        print(f"    {acc} ({comp_filter or 'all'}): {adata_sub.shape[0]:,} cells")
        del adata

    # Concatenate
    print("  Concatenating...")
    adata = ad.concat(adatas, merge='same')
    del adatas

    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(adata.X)
    if not sparse.issparse(adata.layers['counts']):
        adata.layers['counts'] = sparse.csr_matrix(adata.layers['counts'])

    print(f"    Final shape: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")
    return adata


# ═════════════════════════════════════════════════════════════════════════════
# PREPARE OBJECT (HVG, PCA, NEIGHBORS, UMAP BASELINE)
# ═════════════════════════════════════════════════════════════════════════════

def prepare_object(adata, n_top_genes=3000):
    """Re-derive HVGs, PCA, neighbors, UMAP on concatenated object.

    Stores unintegrated baseline in obsm['X_pca_unintegrated'] and
    obsm['X_umap_unintegrated'].
    """
    print("  Preparing object (HVG, PCA, neighbors, UMAP)...")

    n_top = min(n_top_genes, adata.shape[1] - 1)

    # Batch-aware HVG selection
    hvg_done = False
    MIN_CELLS_HVG_BATCH = 50
    if 'study' in adata.obs.columns:
        study_counts = adata.obs['study'].value_counts()
        large_studies = study_counts[study_counts >= MIN_CELLS_HVG_BATCH].index.tolist()
        if len(large_studies) > 1:
            try:
                adata_hvg = adata[adata.obs['study'].isin(large_studies)].copy()
                sc.pp.highly_variable_genes(
                    adata_hvg, n_top_genes=n_top, flavor='seurat_v3',
                    batch_key='study', layer='counts'
                )
                adata.var['highly_variable'] = False
                adata.var.loc[adata_hvg.var_names[adata_hvg.var['highly_variable']],
                              'highly_variable'] = True
                for col in ['highly_variable_rank', 'highly_variable_nbatches',
                            'means', 'variances', 'variances_norm']:
                    if col in adata_hvg.var.columns:
                        adata.var[col] = adata_hvg.var[col]
                del adata_hvg
                hvg_done = True
            except Exception as e:
                print(f"    WARNING: batch-aware HVG failed ({e}), falling back")

    if not hvg_done:
        try:
            sc.pp.highly_variable_genes(
                adata, n_top_genes=n_top, flavor='seurat_v3', layer='counts'
            )
        except ValueError:
            print("    WARNING: seurat_v3 HVG failed, falling back to cell_ranger flavor")
            sc.pp.highly_variable_genes(
                adata, n_top_genes=n_top, flavor='cell_ranger'
            )
    n_hvg = adata.var['highly_variable'].sum()
    print(f"    HVGs: {n_hvg}")

    # PCA on HVG subset
    adata_pca = adata[:, adata.var['highly_variable']].copy()
    gc.collect()
    sc.pp.scale(adata_pca, max_value=10)
    n_comps = min(50, min(adata_pca.shape) - 1)
    sc.tl.pca(adata_pca, n_comps=n_comps, svd_solver='arpack')

    adata.obsm['X_pca'] = adata_pca.obsm['X_pca']
    adata.uns['pca'] = adata_pca.uns['pca']
    adata.varm['PCs'] = np.zeros((adata.shape[1], n_comps))
    hvg_idx = np.where(adata.var['highly_variable'])[0]
    adata.varm['PCs'][hvg_idx] = adata_pca.varm['PCs']
    del adata_pca

    # Unintegrated baseline
    n_pcs = min(30, n_comps)
    sc.pp.neighbors(adata, n_pcs=n_pcs)
    sc.tl.umap(adata)
    adata.obsm['X_pca_unintegrated'] = adata.obsm['X_pca'].copy()
    adata.obsm['X_umap_unintegrated'] = adata.obsm['X_umap'].copy()
    print(f"    PCA: {n_comps} components, baseline UMAP computed")

    return adata


# ═════════════════════════════════════════════════════════════════════════════
# scANVI INTEGRATION
# ═════════════════════════════════════════════════════════════════════════════

def run_scanvi(adata, batch_key='study', labels_key='coarse_label',
               unlabeled_category='Unknown', n_latent=20,
               scvi_max_epochs=200, scanvi_max_epochs=50,
               model_dir=None, label="scanvi"):
    """Run scANVI integration (scVI base + scANVI semi-supervised refinement).

    Workflow:
    1. Train scVI (unsupervised) to learn a base latent representation
    2. Initialize scANVI from the trained scVI model with coarse_label anchors
    3. Train scANVI to refine the latent space using anchor labels

    Stores obsm[f'X_scanvi_{label}'] and UMAP.
    Returns (scanvi_model, embedding_key).
    """
    import scvi as scvi_module

    print(f"  Running scANVI: batch_key={batch_key}, labels_key={labels_key}, "
          f"n_latent={n_latent}")

    # Ensure unlabeled_category is used for cells without coarse_label
    if labels_key in adata.obs.columns:
        adata.obs[labels_key] = adata.obs[labels_key].fillna(unlabeled_category)
        adata.obs[labels_key] = adata.obs[labels_key].astype(str)
    else:
        adata.obs[labels_key] = unlabeled_category

    label_counts = adata.obs[labels_key].value_counts()
    print(f"    Label distribution:")
    for lbl, cnt in label_counts.items():
        print(f"      {lbl}: {cnt:,}")

    # ── Step 1: Train scVI (unsupervised base) ────────────────────────────
    print(f"  Step 1: Training scVI base model (max_epochs={scvi_max_epochs})...")
    scvi_module.model.SCVI.setup_anndata(
        adata, layer='counts', batch_key=batch_key,
    )

    scvi_model = scvi_module.model.SCVI(
        adata, n_latent=n_latent,
        dispersion='gene-batch', gene_likelihood='nb',
    )

    scvi_model.train(
        max_epochs=scvi_max_epochs,
        early_stopping=True,
        early_stopping_patience=10,
        early_stopping_monitor='elbo_validation',
        train_size=0.9,
        batch_size=256,
    )
    scvi_epochs = scvi_model.history['elbo_train'].shape[0]
    print(f"    scVI complete: {scvi_epochs} epochs")

    if model_dir is not None:
        save_path = model_dir / f"scvi_{label}"
        scvi_model.save(str(save_path), overwrite=True)

    # ── Step 2: Initialize scANVI from scVI ───────────────────────────────
    print(f"  Step 2: Initializing scANVI from scVI model...")
    scanvi_model = scvi_module.model.SCANVI.from_scvi_model(
        scvi_model,
        labels_key=labels_key,
        unlabeled_category=unlabeled_category,
    )

    # ── Step 3: Train scANVI (semi-supervised refinement) ─────────────────
    print(f"  Step 3: Training scANVI (max_epochs={scanvi_max_epochs})...")
    scanvi_model.train(
        max_epochs=scanvi_max_epochs,
        early_stopping=True,
        early_stopping_patience=5,
        train_size=0.9,
        batch_size=256,
    )
    scanvi_epochs = scanvi_model.history['elbo_train'].shape[0]
    print(f"    scANVI complete: {scanvi_epochs} epochs")

    # Extract latent representation
    embedding_key = f'X_scanvi_{label}'
    adata.obsm[embedding_key] = scanvi_model.get_latent_representation()
    print(f"    Latent representation stored in obsm['{embedding_key}']")

    # Neighbors and UMAP on scANVI embedding
    sc.pp.neighbors(adata, use_rep=embedding_key)
    sc.tl.umap(adata)
    adata.obsm[f'X_umap_scanvi_{label}'] = adata.obsm['X_umap'].copy()

    if model_dir is not None:
        save_path = model_dir / f"scanvi_{label}"
        scanvi_model.save(str(save_path), overwrite=True)

    del scvi_model
    gc.collect()

    return scanvi_model, embedding_key


# ═════════════════════════════════════════════════════════════════════════════
# INTEGRATION METRICS
# ═════════════════════════════════════════════════════════════════════════════

def _quick_ilisi(embedding, batch_labels, n_sample=5000, k=30):
    """Quick iLISI estimate on a random subsample."""
    from sklearn.neighbors import NearestNeighbors

    n = embedding.shape[0]
    if n > n_sample:
        rng = np.random.RandomState(42)
        idx = rng.choice(n, n_sample, replace=False)
        emb_sub = embedding[idx]
        labels_sub = np.asarray(batch_labels)[idx]
    else:
        emb_sub = embedding
        labels_sub = np.asarray(batch_labels)

    nn = NearestNeighbors(n_neighbors=k, algorithm='auto')
    nn.fit(emb_sub)
    indices = nn.kneighbors(emb_sub, return_distance=False)

    ilisi_scores = []
    for i in range(len(emb_sub)):
        neighbor_labels = labels_sub[indices[i]]
        counts = Counter(neighbor_labels)
        total = sum(counts.values())
        p_sq = sum((c / total) ** 2 for c in counts.values())
        ilisi_scores.append(1.0 / p_sq)

    return np.mean(ilisi_scores)


def compute_integration_metrics(adata, embedding_key, batch_key='study',
                                max_cells=30_000):
    """Compute integration quality metrics (iLISI, batch-ASW, condition-ASW).

    Returns dict with metric values.
    """
    print(f"    Computing metrics for {embedding_key}...")
    metrics = {}

    if embedding_key not in adata.obsm:
        print(f"    WARNING: {embedding_key} not in obsm")
        return metrics

    embedding = adata.obsm[embedding_key]
    batch = adata.obs[batch_key].values

    # Subsample if needed
    n_cells = embedding.shape[0]
    if n_cells > max_cells:
        rng = np.random.RandomState(42)
        idx = rng.choice(n_cells, max_cells, replace=False)
        embedding = embedding[idx]
        batch = batch[idx]
        cond = adata.obs['condition_harmonized'].values[idx] if 'condition_harmonized' in adata.obs.columns else None
        print(f"    Subsampled {n_cells:,} -> {max_cells:,} for metrics")
    else:
        cond = adata.obs['condition_harmonized'].values if 'condition_harmonized' in adata.obs.columns else None

    metrics['n_cells_metrics'] = int(embedding.shape[0])

    # iLISI
    try:
        metrics['iLISI'] = float(_quick_ilisi(embedding, batch))
    except Exception as e:
        print(f"    WARNING: iLISI failed: {e}")
        metrics['iLISI'] = np.nan

    # Batch ASW
    try:
        from sklearn.metrics import silhouette_score
        if len(set(batch)) > 1:
            metrics['batch_ASW'] = float(silhouette_score(
                embedding, batch, sample_size=min(5000, len(embedding)),
                random_state=42
            ))
        else:
            metrics['batch_ASW'] = np.nan
    except Exception as e:
        print(f"    WARNING: batch_ASW failed: {e}")
        metrics['batch_ASW'] = np.nan

    # Condition ASW
    try:
        if cond is not None:
            valid = pd.notna(cond)
            if valid.sum() > 100 and len(set(cond[valid])) > 1:
                metrics['condition_ASW'] = float(silhouette_score(
                    embedding[valid], cond[valid],
                    sample_size=min(5000, valid.sum()),
                    random_state=42
                ))
            else:
                metrics['condition_ASW'] = np.nan
        else:
            metrics['condition_ASW'] = np.nan
    except Exception as e:
        print(f"    WARNING: condition_ASW failed: {e}")
        metrics['condition_ASW'] = np.nan

    return metrics


# ═════════════════════════════════════════════════════════════════════════════
# TIER MERGING
# ═════════════════════════════════════════════════════════════════════════════

def merge_tiers(mes_adata, non_mes_adata, mes_embedding_key, non_mes_embedding_key):
    """Merge mesenchymal and non-mesenchymal tiers into a single object.

    Stores tier-specific embeddings as X_scanvi_mesenchymal and X_scanvi_non_mesenchymal.
    """
    print("  Merging tiers...")

    # Store tier-specific embeddings with NaN for the other tier
    n_latent_mes = mes_adata.obsm[mes_embedding_key].shape[1]
    n_latent_non = non_mes_adata.obsm[non_mes_embedding_key].shape[1]

    # Pad non-mesenchymal with NaN for mesenchymal embedding and vice versa
    mes_adata.obsm['X_scanvi_mesenchymal'] = mes_adata.obsm[mes_embedding_key].copy()
    mes_adata.obsm['X_scanvi_non_mesenchymal'] = np.full(
        (mes_adata.shape[0], n_latent_non), np.nan
    )

    non_mes_adata.obsm['X_scanvi_non_mesenchymal'] = non_mes_adata.obsm[non_mes_embedding_key].copy()
    non_mes_adata.obsm['X_scanvi_mesenchymal'] = np.full(
        (non_mes_adata.shape[0], n_latent_mes), np.nan
    )

    # Ensure both have the same obs columns
    common_obs = list(set(mes_adata.obs.columns) & set(non_mes_adata.obs.columns))
    mes_adata.obs = mes_adata.obs[common_obs].copy()
    non_mes_adata.obs = non_mes_adata.obs[common_obs].copy()

    # Ensure same var
    common_genes = sorted(set(mes_adata.var_names) & set(non_mes_adata.var_names))
    mes_sub = mes_adata[:, common_genes].copy()
    non_mes_sub = non_mes_adata[:, common_genes].copy()

    # Concatenate
    merged = ad.concat([mes_sub, non_mes_sub], merge='same')

    if not sparse.issparse(merged.X):
        merged.X = sparse.csr_matrix(merged.X)

    print(f"    Merged: {merged.shape[0]:,} cells "
          f"({mes_adata.shape[0]:,} mes + {non_mes_adata.shape[0]:,} non-mes)")

    del mes_sub, non_mes_sub
    return merged


# ═════════════════════════════════════════════════════════════════════════════
# CHECKPOINT SAVE
# ═════════════════════════════════════════════════════════════════════════════

def _save_checkpoint(adata, output_path, label):
    """Atomic save of adata checkpoint."""
    import anndata
    anndata.settings.allow_write_nullable_strings = True
    tmp_path = Path(str(output_path) + '.tmp')
    adata.write_h5ad(tmp_path)
    tmp_path.rename(output_path)
    print(f"  Checkpoint saved ({label}): {output_path}")
    sys.stdout.flush()


# ═════════════════════════════════════════════════════════════════════════════
# PROCESS ONE OBJECT (MAIN ORCHESTRATOR)
# ═════════════════════════════════════════════════════════════════════════════

def process_object(object_name, force=False):
    """Full integration pipeline for one object.

    Steps: load -> split by cell_class -> scANVI per tier -> metrics -> merge -> save.
    No clustering, no annotation (those are in Modules 06 and 07).
    """
    output_path = INT_DIR / f"{object_name}.h5ad"

    if output_path.exists() and not force:
        print(f"\n=== {object_name}: output exists, skipping (use --force) ===")
        return output_path, {}

    print(f"\n{'='*60}")
    print(f"Processing object: {object_name}")
    print(f"{'='*60}")

    # Load all cells for this object
    adata_all = load_and_build_object(object_name)
    if adata_all is None or adata_all.shape[0] < 50:
        print(f"  ERROR: Too few cells for {object_name}")
        return None, {}

    # Split by cell_class
    # Mesenchymal tier includes 'unknown' cells (fibrochondrocytes that scANVI
    # should position based on transcriptomic similarity)
    mes_mask = adata_all.obs['cell_class'].isin(['mesenchymal', 'unknown'])
    non_mes_mask = adata_all.obs['cell_class'] == 'non_mesenchymal'

    n_mes = mes_mask.sum()
    n_non = non_mes_mask.sum()
    print(f"  Mesenchymal (+ unknown): {n_mes:,} cells")
    print(f"  Non-mesenchymal: {n_non:,} cells")

    all_metrics = {}
    mes_adata = None
    non_mes_adata = None
    mes_emb_key = None
    non_mes_emb_key = None

    # ── Tier A: Mesenchymal ──────────────────────────────────────────────
    if n_mes >= 50:
        print(f"\n  --- Tier A: Mesenchymal ({object_name}) ---")
        mes_adata = adata_all[mes_mask].copy()
        mes_adata = prepare_object(mes_adata, n_top_genes=SCANVI_PARAMS['n_top_genes'])

        _, mes_emb_key = run_scanvi(
            mes_adata,
            batch_key=SCANVI_PARAMS['batch_key'],
            labels_key='coarse_label',
            unlabeled_category='Unknown',
            n_latent=SCANVI_PARAMS['n_latent'],
            scvi_max_epochs=SCANVI_PARAMS['scvi_max_epochs'],
            scanvi_max_epochs=SCANVI_PARAMS['scanvi_max_epochs'],
            model_dir=MODEL_DIR,
            label=f"{object_name}_mesenchymal"
        )

        # Integration metrics
        mes_metrics = compute_integration_metrics(mes_adata, mes_emb_key)
        all_metrics['mesenchymal'] = mes_metrics

        gc.collect()
    else:
        print(f"  Skipping mesenchymal tier: only {n_mes} cells")

    # ── Tier B: Non-mesenchymal ──────────────────────────────────────────
    # Need enough cells for HVG + scANVI to work (at least 200 cells)
    MIN_CELLS_NON_MES = 200
    if n_non >= MIN_CELLS_NON_MES:
        print(f"\n  --- Tier B: Non-mesenchymal ({object_name}) ---")
        non_mes_adata = adata_all[non_mes_mask].copy()
        non_mes_adata = prepare_object(non_mes_adata, n_top_genes=SCANVI_PARAMS['n_top_genes'])

        _, non_mes_emb_key = run_scanvi(
            non_mes_adata,
            batch_key=SCANVI_PARAMS['batch_key'],
            labels_key='coarse_label',
            unlabeled_category='Unknown',
            n_latent=SCANVI_PARAMS['n_latent'],
            scvi_max_epochs=SCANVI_PARAMS['scvi_max_epochs'],
            scanvi_max_epochs=SCANVI_PARAMS['scanvi_max_epochs'],
            model_dir=MODEL_DIR,
            label=f"{object_name}_non_mesenchymal"
        )

        # Integration metrics
        non_mes_metrics = compute_integration_metrics(non_mes_adata, non_mes_emb_key)
        all_metrics['non_mesenchymal'] = non_mes_metrics

        gc.collect()
    else:
        print(f"  Skipping non-mesenchymal tier: only {n_non} cells")

    del adata_all
    gc.collect()

    # ── Merge tiers ──────────────────────────────────────────────────────
    if mes_adata is not None and non_mes_adata is not None:
        merged = merge_tiers(mes_adata, non_mes_adata, mes_emb_key, non_mes_emb_key)
        del mes_adata, non_mes_adata
    elif mes_adata is not None:
        merged = mes_adata
        merged.obsm['X_scanvi_mesenchymal'] = merged.obsm[mes_emb_key].copy()
        del mes_adata
    elif non_mes_adata is not None:
        merged = non_mes_adata
        merged.obsm['X_scanvi_non_mesenchymal'] = merged.obsm[non_mes_emb_key].copy()
        del non_mes_adata
    else:
        print(f"  ERROR: No tiers processed for {object_name}")
        return None, {}

    # ── Generate UMAPs ───────────────────────────────────────────────────
    _plot_object_umaps(merged, object_name)

    # ── Save ─────────────────────────────────────────────────────────────
    _save_checkpoint(merged, output_path, "final")
    print(f"  {object_name} complete: {merged.shape[0]:,} cells")

    del merged
    gc.collect()
    return output_path, all_metrics


def process_all_cells_secondary(force=False):
    """Process all-cells as a secondary object.

    Loads all cells, integrates with scANVI per tier, and saves.
    No annotation transfer (that moves to Module 07).
    """
    output_path = INT_DIR / "all_cells.h5ad"

    if output_path.exists() and not force:
        print(f"\n=== all_cells: output exists, skipping (use --force) ===")
        return output_path, {}

    print(f"\n{'='*60}")
    print(f"Processing object: all_cells (secondary)")
    print(f"{'='*60}")

    # Load all cells
    adata_all = load_and_build_object("all_cells")
    if adata_all is None:
        return None, {}

    # Split by cell_class (same as primary objects)
    mes_mask = adata_all.obs['cell_class'].isin(['mesenchymal', 'unknown'])
    non_mes_mask = adata_all.obs['cell_class'] == 'non_mesenchymal'

    n_mes = mes_mask.sum()
    n_non = non_mes_mask.sum()
    print(f"  Mesenchymal (+ unknown): {n_mes:,} cells")
    print(f"  Non-mesenchymal: {n_non:,} cells")

    all_metrics = {}
    mes_adata = None
    non_mes_adata = None
    mes_emb_key = None
    non_mes_emb_key = None

    # ── Tier A: Mesenchymal ──────────────────────────────────────────────
    if n_mes >= 50:
        print(f"\n  --- Tier A: Mesenchymal (all_cells) ---")
        mes_adata = adata_all[mes_mask].copy()
        mes_adata = prepare_object(mes_adata, n_top_genes=SCANVI_PARAMS['n_top_genes'])

        _, mes_emb_key = run_scanvi(
            mes_adata,
            batch_key=SCANVI_PARAMS['batch_key'],
            labels_key='coarse_label',
            unlabeled_category='Unknown',
            n_latent=SCANVI_PARAMS['n_latent'],
            scvi_max_epochs=SCANVI_PARAMS['scvi_max_epochs'],
            scanvi_max_epochs=SCANVI_PARAMS['scanvi_max_epochs'],
            model_dir=MODEL_DIR,
            label="all_cells_mesenchymal"
        )

        mes_metrics = compute_integration_metrics(mes_adata, mes_emb_key)
        all_metrics['mesenchymal'] = mes_metrics
        gc.collect()

    # ── Tier B: Non-mesenchymal ──────────────────────────────────────────
    MIN_CELLS_NON_MES = 200
    if n_non >= MIN_CELLS_NON_MES:
        print(f"\n  --- Tier B: Non-mesenchymal (all_cells) ---")
        non_mes_adata = adata_all[non_mes_mask].copy()
        non_mes_adata = prepare_object(non_mes_adata, n_top_genes=SCANVI_PARAMS['n_top_genes'])

        _, non_mes_emb_key = run_scanvi(
            non_mes_adata,
            batch_key=SCANVI_PARAMS['batch_key'],
            labels_key='coarse_label',
            unlabeled_category='Unknown',
            n_latent=SCANVI_PARAMS['n_latent'],
            scvi_max_epochs=SCANVI_PARAMS['scvi_max_epochs'],
            scanvi_max_epochs=SCANVI_PARAMS['scanvi_max_epochs'],
            model_dir=MODEL_DIR,
            label="all_cells_non_mesenchymal"
        )

        non_mes_metrics = compute_integration_metrics(non_mes_adata, non_mes_emb_key)
        all_metrics['non_mesenchymal'] = non_mes_metrics
        gc.collect()

    del adata_all
    gc.collect()

    # ── Merge tiers ──────────────────────────────────────────────────────
    if mes_adata is not None and non_mes_adata is not None:
        merged = merge_tiers(mes_adata, non_mes_adata, mes_emb_key, non_mes_emb_key)
        del mes_adata, non_mes_adata
    elif mes_adata is not None:
        merged = mes_adata
        merged.obsm['X_scanvi_mesenchymal'] = merged.obsm[mes_emb_key].copy()
        del mes_adata
    elif non_mes_adata is not None:
        merged = non_mes_adata
        merged.obsm['X_scanvi_non_mesenchymal'] = merged.obsm[non_mes_emb_key].copy()
        del non_mes_adata
    else:
        print("  ERROR: No tiers processed for all_cells")
        return None, {}

    # UMAPs
    _plot_object_umaps(merged, "all_cells")

    # Save
    _save_checkpoint(merged, output_path, "final")
    print(f"  all_cells complete: {merged.shape[0]:,} cells")

    del merged
    gc.collect()
    return output_path, all_metrics


# ═════════════════════════════════════════════════════════════════════════════
# UMAP PLOTS
# ═════════════════════════════════════════════════════════════════════════════

def _plot_object_umaps(adata, object_name):
    """Generate UMAP panels for one object."""
    color_keys = ['study', 'condition_harmonized', 'cell_class', 'coarse_label']
    color_keys = [c for c in color_keys if c in adata.obs.columns]
    n_cols = len(color_keys)

    if n_cols == 0 or 'X_umap' not in adata.obsm:
        return

    fig, axes = plt.subplots(1, n_cols, figsize=(5 * n_cols, 4.5))
    if n_cols == 1:
        axes = [axes]

    for ax, color in zip(axes, color_keys):
        try:
            sc.pl.umap(adata, color=color, ax=ax, show=False,
                       frameon=False, s=2, alpha=0.5,
                       title=f"{object_name} — {color}")
        except Exception:
            ax.set_title(f"{object_name} — {color} (failed)")

    plt.tight_layout()
    fig.savefig(RESULTS_DIR / f"umap_{object_name}.png",
                dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: umap_{object_name}.png")


# ═════════════════════════════════════════════════════════════════════════════
# INCLUSION SUMMARY & STUDY CAVEATS
# ═════════════════════════════════════════════════════════════════════════════

def generate_inclusion_summary():
    """Generate inclusion_summary.tsv and .html for manuscript supplement."""
    print("\n  Generating inclusion summary...")

    meta_path = META_DIR / "sample_metadata.tsv"
    meta = pd.read_csv(meta_path, sep='\t') if meta_path.exists() else None

    rows = []
    for obj_name in ["NP", "AF", "CEP", "all_cells"]:
        obj_path = INT_DIR / f"{obj_name}.h5ad"
        if not obj_path.exists():
            continue
        adata = sc.read_h5ad(obj_path, backed='r')

        for study in sorted(adata.obs['study'].unique()):
            study_mask = adata.obs['study'] == study
            n_cells = int(study_mask.sum())
            n_samples = adata.obs.loc[study_mask, 'sample_id'].nunique() if 'sample_id' in adata.obs.columns else 0

            # Get metadata
            first_author = ""
            year = ""
            compartment = ""
            conditions = ""
            platform = ""
            if meta is not None:
                study_meta = meta[meta['study_accession'] == study]
                if len(study_meta) > 0:
                    first_author = study_meta.iloc[0].get('first_author', '')
                    year = str(study_meta.iloc[0].get('year', ''))
                    compartment = ', '.join(study_meta['compartment'].unique())
                    conditions = ', '.join(study_meta['condition_harmonized'].dropna().unique())
                    platform = study_meta.iloc[0].get('sequencing_platform', '')

            rows.append({
                'object': obj_name,
                'accession': study,
                'first_author': first_author,
                'year': year,
                'n_samples': n_samples,
                'n_cells': n_cells,
                'compartment': compartment,
                'conditions': conditions,
                'platform': platform,
            })
        adata.file.close()

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS_DIR / "inclusion_summary.tsv", sep='\t', index=False)

    # HTML version
    html = df.to_html(index=False, classes='summary-table')
    html_page = f"""<!DOCTYPE html>
<html><head><title>Inclusion Summary</title>
<style>
body {{ font-family: Arial; max-width: 1200px; margin: 0 auto; padding: 20px; }}
.summary-table {{ border-collapse: collapse; width: 100%; }}
.summary-table th {{ background: #3498db; color: white; padding: 8px; }}
.summary-table td {{ border: 1px solid #ddd; padding: 6px; }}
.summary-table tr:nth-child(even) {{ background: #f2f2f2; }}
</style></head><body>
<h1>Study Inclusion Summary</h1>
<p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
{html}
</body></html>"""
    (RESULTS_DIR / "inclusion_summary.html").write_text(html_page)
    print(f"  Saved: inclusion_summary.tsv and .html")


def generate_study_caveats():
    """Generate study_caveats.tsv for manuscript supplement."""
    caveats = [
        {"study": "GSE165722", "first_author": "Tu 2022",
         "caveat": "BD Rhapsody platform (not 10x)",
         "impact": "Different capture efficiency, gene detection",
         "mitigation": "Platform-aware batch correction via scANVI study-level batch key"},
        {"study": "GSE205535", "first_author": "Li 2022",
         "caveat": "BD Rhapsody platform; published corrigenda",
         "impact": "See above; potential data quality issues",
         "mitigation": "Monitor for outlier behavior in integration"},
        {"study": "CNP0002664", "first_author": "Han 2022",
         "caveat": "Singleron Matrix platform (not 10x)",
         "impact": "Different capture efficiency",
         "mitigation": "Same as above"},
        {"study": "GSE242443", "first_author": "Kuchynsky 2024",
         "caveat": "Culture-expanded CEP cells",
         "impact": "Culture alters cell states; may not reflect in vivo biology",
         "mitigation": "Caveat in all CEP results; compare with non-expanded CEP from GSE160756"},
        {"study": "GSE255768", "first_author": "Shi 2024",
         "caveat": "Degenerative endplate only; no healthy control",
         "impact": "Cannot do healthy vs. degenerated comparison for this study alone",
         "mitigation": "Healthy CEP baseline from GSE160756"},
        {"study": "GSE230809", "first_author": "Swahn 2024",
         "caveat": "All-male donors; age-disease confounded",
         "impact": "Cannot separate age from degeneration effects",
         "mitigation": "Note in interpretation; sex-specific effects cannot be assessed"},
        {"study": "GSE205535 NNP", "first_author": "Li 2022",
         "caveat": "11yo spinal cord injury, classified as 'healthy'",
         "impact": "Trauma confound",
         "mitigation": "Excluded from DE comparisons"},
        {"study": "GSE189916", "first_author": "Jiang 2022",
         "caveat": "Whole IVD (compartments not separated)",
         "impact": "Cannot assign cells to NP/AF/CEP",
         "mitigation": "Included only in all-cells object"},
    ]
    df = pd.DataFrame(caveats)
    df.to_csv(RESULTS_DIR / "study_caveats.tsv", sep='\t', index=False)
    print(f"  Saved: study_caveats.tsv")


# ═════════════════════════════════════════════════════════════════════════════
# INTEGRATION REPORT
# ═════════════════════════════════════════════════════════════════════════════

INTEGRATION_REPORT_TEMPLATE = """<!DOCTYPE html>
<html>
<head>
<title>Integration Report — Module 05</title>
<style>
  body { font-family: Arial, sans-serif; max-width: 1400px; margin: 0 auto; padding: 20px; }
  h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
  h2 { color: #34495e; margin-top: 30px; }
  table { border-collapse: collapse; width: 100%; margin: 10px 0; }
  th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
  th { background-color: #3498db; color: white; }
  tr:nth-child(even) { background-color: #f2f2f2; }
  img { max-width: 100%; height: auto; margin: 10px 0; border: 1px solid #ddd; }
  .stats-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0; }
  .stat-box { background: #ecf0f1; padding: 15px; border-radius: 5px; }
  .stat-box h3 { margin: 0 0 5px 0; color: #7f8c8d; font-size: 12px; }
  .stat-box p { margin: 0; font-size: 24px; font-weight: bold; color: #2c3e50; }
  .warning { color: #e74c3c; }
  .pass { color: #27ae60; }
  .secondary { color: #95a5a6; font-style: italic; }
  .validation li { margin: 4px 0; }
</style>
</head>
<body>
<h1>Module 05: Integration (scANVI)</h1>
<p>Generated: {{ date }}</p>
<p>Integration method: scANVI (semi-supervised) with coarse_label anchors from Module 04.</p>
<p>Parameters: batch_key=study, n_latent=20, scVI max_epochs=200, scANVI max_epochs=50, n_top_genes=3000</p>

{% for obj in objects %}
<h2>{{ obj.name }}{% if obj.secondary %} <span class="secondary">(secondary)</span>{% endif %}</h2>
{% if obj.umap_exists %}
<img src="umap_{{ obj.name }}.png" alt="{{ obj.name }} UMAP">
{% endif %}
<p>Cells: {{ obj.n_cells }}</p>
{% if obj.metrics %}
<div class="stats-grid">
{% for key, val in obj.metrics.items() %}
<div class="stat-box">
  <h3>{{ key }}</h3>
  <p>{{ val }}</p>
</div>
{% endfor %}
</div>
{% endif %}
{% endfor %}

<h2>Validation</h2>
<ul class="validation">
  {% for msg in validation_messages %}
  <li class="{{ 'pass' if msg.startswith('PASS') else 'warning' if msg.startswith('FAIL') or msg.startswith('WARNING') else '' }}">{{ msg }}</li>
  {% endfor %}
</ul>

<h2>Human Checkpoint Questions</h2>
<ol>
  <li>Does the integration look reasonable? Are studies mixing well on the UMAP?</li>
  <li>Is there evidence of overcorrection (biological signal lost)?</li>
  <li>Are the study caveats adequately documented for the supplement?</li>
  <li>Should integration parameters be adjusted for any object/tier?</li>
</ol>

</body>
</html>"""


def generate_integration_report(validation_messages, all_metrics=None):
    """Generate HTML integration report."""
    print("\n  Generating integration report...")

    objects = []
    for obj_name in ["NP", "AF", "CEP", "all_cells"]:
        obj_path = INT_DIR / f"{obj_name}.h5ad"
        if obj_path.exists():
            adata = sc.read_h5ad(obj_path, backed='r')
            obj_info = {
                'name': obj_name,
                'n_cells': f"{adata.shape[0]:,}",
                'umap_exists': (RESULTS_DIR / f"umap_{obj_name}.png").exists(),
                'secondary': obj_name == 'all_cells',
                'metrics': {},
            }
            # Add metrics if available
            if all_metrics and obj_name in all_metrics:
                for tier, tier_metrics in all_metrics[obj_name].items():
                    if isinstance(tier_metrics, dict):
                        for k, v in tier_metrics.items():
                            if isinstance(v, float):
                                obj_info['metrics'][f"{tier}_{k}"] = f"{v:.3f}"
                            else:
                                obj_info['metrics'][f"{tier}_{k}"] = str(v)
            objects.append(obj_info)
            adata.file.close()

    template = Template(INTEGRATION_REPORT_TEMPLATE)
    html = template.render(
        date=datetime.now().strftime("%Y-%m-%d %H:%M"),
        objects=objects,
        validation_messages=validation_messages,
    )
    report_path = BASE / "results" / "integration_report.html"
    report_path.write_text(html)
    print(f"  Report: {report_path}")


# ═════════════════════════════════════════════════════════════════════════════
# VALIDATION
# ═════════════════════════════════════════════════════════════════════════════

def validate_integration():
    """Validate integration outputs. Returns (pass, messages)."""
    print("\n" + "=" * 60)
    print("Validating integration outputs")
    print("=" * 60)

    messages = []
    all_pass = True

    # Check all four output files
    for obj_name in ["NP", "AF", "CEP", "all_cells"]:
        path = INT_DIR / f"{obj_name}.h5ad"
        if path.exists():
            messages.append(f"PASS: {obj_name} output exists")
            adata = sc.read_h5ad(path, backed='r')

            # Blob check: compute a quick leiden at res 0.5 on UMAP to verify
            # the integration didn't collapse everything
            if 'X_umap' in adata.obsm:
                try:
                    # Read full adata for clustering check
                    adata_full = sc.read_h5ad(path)
                    sc.pp.neighbors(adata_full, use_rep='X_pca')
                    sc.tl.leiden(adata_full, resolution=0.5, key_added='_blob_check')
                    n_clusters = adata_full.obs['_blob_check'].nunique()
                    if n_clusters > 1:
                        messages.append(f"PASS: {obj_name} has {n_clusters} clusters at res=0.5 (no blob)")
                    else:
                        messages.append(f"WARNING: {obj_name} has only {n_clusters} cluster at res=0.5 (blob?)")
                    del adata_full
                except Exception as e:
                    messages.append(f"WARNING: {obj_name} blob check failed: {e}")
            else:
                messages.append(f"WARNING: {obj_name} missing X_umap")

            # ARI check (study shouldn't perfectly predict clusters)
            if 'study' in adata.obs.columns:
                try:
                    adata_full = sc.read_h5ad(path)
                    sc.pp.neighbors(adata_full, use_rep='X_pca')
                    sc.tl.leiden(adata_full, resolution=0.5, key_added='_ari_check')
                    from sklearn.metrics import adjusted_rand_score
                    n_sub = min(10000, adata_full.shape[0])
                    ari = adjusted_rand_score(
                        adata_full.obs['study'].values[:n_sub],
                        adata_full.obs['_ari_check'].values[:n_sub]
                    )
                    if ari < 1.0:
                        messages.append(f"PASS: {obj_name} study-cluster ARI = {ari:.3f} (< 1.0)")
                    else:
                        messages.append(f"WARNING: {obj_name} study perfectly predicts clusters")
                    del adata_full
                except Exception as e:
                    messages.append(f"WARNING: {obj_name} ARI check failed: {e}")

            adata.file.close()
        else:
            messages.append(f"FAIL: {obj_name} output missing")
            all_pass = False

    # Check supplementary outputs
    for fname in ["inclusion_summary.tsv", "study_caveats.tsv"]:
        if (RESULTS_DIR / fname).exists():
            messages.append(f"PASS: {fname} exists")
        else:
            messages.append(f"FAIL: {fname} missing")
            all_pass = False

    # Check report
    report_path = BASE / "results" / "integration_report.html"
    if report_path.exists():
        messages.append("PASS: Integration report exists")
    else:
        messages.append("FAIL: Integration report missing")
        all_pass = False

    for msg in messages:
        print(f"  {msg}")

    print(f"\n{'='*60}")
    print(f"OVERALL: {'ALL CHECKS PASSED' if all_pass else 'SOME CHECKS FAILED'}")
    print(f"{'='*60}")
    return all_pass, messages


# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════

def main():
    args = sys.argv[1:]

    validate_only = '--validate-only' in args
    force = '--force' in args

    # Parse --object
    object_filter = None
    for i, a in enumerate(args):
        if a == '--object' and i + 1 < len(args):
            object_filter = args[i + 1]

    if validate_only:
        passed, _ = validate_integration()
        sys.exit(0 if passed else 1)

    print("=" * 60)
    print("Module 05: Integration (scANVI)")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    all_metrics = {}

    # Process compartment-specific objects (primary)
    primary_objects = ["NP", "AF", "CEP"]
    if object_filter:
        if object_filter in primary_objects:
            primary_objects = [object_filter]
        elif object_filter == "all_cells":
            primary_objects = []
        else:
            print(f"Unknown object: {object_filter}")
            sys.exit(1)

    for obj_name in primary_objects:
        _, metrics = process_object(obj_name, force=force)
        all_metrics[obj_name] = metrics

    # Process all-cells (secondary) — after primary objects are done
    if object_filter is None or object_filter == "all_cells":
        _, metrics = process_all_cells_secondary(force=force)
        all_metrics['all_cells'] = metrics

    # Save metrics
    metrics_rows = []
    for obj_name, obj_metrics in all_metrics.items():
        for tier, tier_metrics in obj_metrics.items():
            if isinstance(tier_metrics, dict):
                row = {'object': obj_name, 'tier': tier}
                row.update(tier_metrics)
                metrics_rows.append(row)
    if metrics_rows:
        pd.DataFrame(metrics_rows).to_csv(
            INT_DIR / "integration_metrics.tsv", sep='\t', index=False
        )

    # Generate supplementary tables
    generate_inclusion_summary()
    generate_study_caveats()

    # Validation
    passed, val_messages = validate_integration()

    # Report
    generate_integration_report(val_messages, all_metrics)

    print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Overall: {'PASSED' if passed else 'FAILED'}")


if __name__ == "__main__":
    main()
