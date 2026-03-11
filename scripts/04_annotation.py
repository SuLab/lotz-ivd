#!/usr/bin/env python3
"""Module 04: Coarse Cell Classification for Integration Anchors.

Classifies cells into 5 anchor categories (Chondrocyte_like, Fibroblast_like,
Immune, Endothelial, Pericyte_SMC) plus Unknown. These coarse labels serve as
seed labels for scANVI semi-supervised integration in Module 05.
Fine-grained cell type annotation happens after integration (in Module 07).

Usage:
    python3 scripts/04_annotation.py                  # Classify all datasets
    python3 scripts/04_annotation.py GSE255768        # Classify one dataset
    python3 scripts/04_annotation.py --validate-only  # Run validation only
"""

import sys
import warnings
from pathlib import Path
from datetime import datetime

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

anndata.settings.allow_write_nullable_strings = True
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from jinja2 import Template

warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning, module='scanpy')

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent.parent
PROC_DIR = BASE / "data" / "processed"
ANN_DIR = BASE / "results" / "annotations"
META_DIR = BASE / "metadata"

ANN_DIR.mkdir(parents=True, exist_ok=True)

# ── Dataset list ─────────────────────────────────────────────────────────────
# GSE233666 excluded entirely (herniated-only, not the focus of this project)
ALL_ACCESSIONS = [
    "GSE160756", "GSE165722", "GSE189916", "GSE199866", "GSE205535",
    "CNP0002664", "GSE244889", "GSE251686", "GSE255768",
    "GSE230809", "GSE242443",
]

# ── Clustering resolution for majority voting ───────────────────────────────
CLASSIFICATION_RESOLUTION = 0.5
FALLBACK_RESOLUTION = 0.3

# ── Gene aliases ─────────────────────────────────────────────────────────────
GENE_ALIASES = {
    "T":    ["T", "TBXT"],
    "TBXT": ["TBXT", "T"],
    "SP7":  ["SP7", "OSX"],
}

# ── Marker gene panels ──────────────────────────────────────────────────────
# Chondrocyte-like markers
CHONDROCYTE_GENES = ["COL2A1", "ACAN", "SOX9"]

# Fibroblast-like markers
FIBROBLAST_GENES = ["COL1A1", "COL1A2", "DCN", "LUM"]

# Non-mesenchymal markers
IMMUNE_PRIMARY = ["PTPRC"]
IMMUNE_SUPPORTING = [
    "CD3D", "CD3E", "CD68", "CD14", "CSF1R",
    "CD79A", "MS4A1", "KIT", "TPSAB1", "NKG7", "GNLY",
]
ENDOTHELIAL_GENES = ["PECAM1", "VWF", "CDH5"]
PERICYTE_GENES = ["RGS5", "PDGFRB"]

# Rescue markers (prevent stressed disc cells from being called non-mesenchymal)
RESCUE_GENES = ["ACAN", "SOX9"]

# All classification markers (flat list for reporting)
ALL_MARKER_GENES = (
    CHONDROCYTE_GENES + FIBROBLAST_GENES +
    IMMUNE_PRIMARY + IMMUNE_SUPPORTING +
    ENDOTHELIAL_GENES + PERICYTE_GENES
)

# Coarse label palette (6 categories)
COARSE_PALETTE = {
    "Chondrocyte_like": "#2196F3",
    "Fibroblast_like": "#4CAF50",
    "Immune": "#E91E63",
    "Endothelial": "#FF9800",
    "Pericyte_SMC": "#9C27B0",
    "Unknown": "#9E9E9E",
}

VALID_COARSE_LABELS = set(COARSE_PALETTE.keys())
VALID_CELL_CLASSES = {"mesenchymal", "non_mesenchymal", "unknown"}


# ═════════════════════════════════════════════════════════════════════════════
# GENE RESOLUTION
# ═════════════════════════════════════════════════════════════════════════════

def resolve_gene(name, var_names):
    """Resolve a gene name to the form present in var_names, handling aliases."""
    if name in var_names:
        return name
    aliases = GENE_ALIASES.get(name, [])
    for alias in aliases:
        if alias in var_names:
            return alias
    return None


def resolve_gene_list(gene_list, var_names):
    """Resolve all genes in a list, returning (resolved, missing)."""
    var_set = set(var_names)
    resolved, missing = [], []
    for gene in gene_list:
        match = resolve_gene(gene, var_set)
        if match is not None:
            resolved.append(match)
        else:
            missing.append(gene)
    return resolved, missing


# ═════════════════════════════════════════════════════════════════════════════
# SAMPLE FILTERING
# ═════════════════════════════════════════════════════════════════════════════

def get_excluded_samples():
    """Return set of sample_ids to exclude (neonatal samples from GSE189916)."""
    meta_path = META_DIR / "sample_metadata.tsv"
    excluded = set()
    if meta_path.exists():
        meta = pd.read_csv(meta_path, sep='\t')
        # Exclude neonatal samples from GSE189916
        neonatal = meta[
            (meta['study_accession'] == 'GSE189916') &
            (meta['age_group'] == 'neonatal')
        ]
        excluded.update(neonatal['sample_id'].values)
    return excluded


def filter_excluded_samples(adata, accession, excluded_samples):
    """Remove excluded samples from an AnnData object. Returns filtered adata."""
    if not excluded_samples or 'sample_id' not in adata.obs.columns:
        return adata

    mask = ~adata.obs['sample_id'].isin(excluded_samples)
    n_removed = (~mask).sum()
    if n_removed > 0:
        print(f"    Excluded {n_removed} cells from {n_removed} neonatal samples")
        adata = adata[mask].copy()
    return adata


# ═════════════════════════════════════════════════════════════════════════════
# CLASSIFICATION
# ═════════════════════════════════════════════════════════════════════════════

def _get_gene_expression(adata, gene):
    """Extract expression vector for a gene, handling sparse matrices."""
    if gene not in adata.var_names:
        return None
    col = adata[:, gene].X
    if sparse.issparse(col):
        return col.toarray().ravel()
    return np.asarray(col).ravel()


def _compute_top_percentile_threshold(expr, percentile=90):
    """Compute the threshold for the top Nth percentile of nonzero values.

    Returns the value at the given percentile of nonzero expression values.
    Cells above this threshold are considered to have strong expression.
    """
    nonzero = expr[expr > 0]
    if len(nonzero) == 0:
        return np.inf  # No expression at all — threshold is unreachable
    return np.percentile(nonzero, percentile)


def score_cell_class(adata):
    """Score each cell for chondrocyte and fibroblast marker expression.

    Uses sc.tl.score_genes() for chondrocyte and fibroblast marker sets.
    Stores scores in obs['score_chondrocyte'] and obs['score_fibroblast'].
    """
    var_names = set(adata.var_names)

    # Score chondrocyte markers
    chondro_resolved, chondro_missing = resolve_gene_list(CHONDROCYTE_GENES, var_names)
    if chondro_resolved:
        sc.tl.score_genes(adata, gene_list=chondro_resolved, score_name='score_chondrocyte')
    else:
        adata.obs['score_chondrocyte'] = 0.0
    print(f"    Chondrocyte markers: {len(chondro_resolved)}/{len(CHONDROCYTE_GENES)} available"
          + (f" (missing: {', '.join(chondro_missing)})" if chondro_missing else ""))

    # Score fibroblast markers
    fibro_resolved, fibro_missing = resolve_gene_list(FIBROBLAST_GENES, var_names)
    if fibro_resolved:
        sc.tl.score_genes(adata, gene_list=fibro_resolved, score_name='score_fibroblast')
    else:
        adata.obs['score_fibroblast'] = 0.0
    print(f"    Fibroblast markers: {len(fibro_resolved)}/{len(FIBROBLAST_GENES)} available"
          + (f" (missing: {', '.join(fibro_missing)})" if fibro_missing else ""))

    # Report non-mesenchymal marker availability
    nonmes_genes = IMMUNE_PRIMARY + IMMUNE_SUPPORTING + ENDOTHELIAL_GENES + PERICYTE_GENES
    nonmes_resolved, nonmes_missing = resolve_gene_list(nonmes_genes, var_names)
    print(f"    Non-mesenchymal markers: {len(nonmes_resolved)}/{len(nonmes_genes)} available"
          + (f" (missing: {', '.join(nonmes_missing)})" if nonmes_missing else ""))

    return {
        "chondrocyte": {"available": chondro_resolved, "missing": chondro_missing},
        "fibroblast": {"available": fibro_resolved, "missing": fibro_missing},
        "non_mesenchymal": {"available": nonmes_resolved, "missing": nonmes_missing},
    }


def classify_cells(adata):
    """Assign coarse_label based on hierarchical marker rules.

    Classification hierarchy (first match wins):
    1. Immune: PTPRC in top 10th percentile, OR ≥2 immune supporting markers
       each in the top 10th percentile. Must NOT co-express ACAN or SOX9.
    2. Endothelial: PECAM1, VWF, or CDH5 in top 10th percentile.
       Must NOT co-express ACAN or SOX9.
    3. Pericyte_SMC: RGS5 AND PDGFRB co-expressed (both in top 10th percentile).
       Must NOT co-express ACAN or SOX9.
    4. Chondrocyte_like: chondrocyte score > 2x fibroblast score AND > 0.
    5. Fibroblast_like: fibroblast score > 2x chondrocyte score AND > 0.
    6. Unknown: everything else.

    Derives cell_class from coarse_label:
    - Chondrocyte_like, Fibroblast_like -> "mesenchymal"
    - Immune, Endothelial, Pericyte_SMC -> "non_mesenchymal"
    - Unknown -> "unknown"
    """
    n_cells = adata.shape[0]
    labels = np.full(n_cells, "Unknown", dtype=object)

    score_chondro = adata.obs['score_chondrocyte'].values
    score_fibro = adata.obs['score_fibroblast'].values

    # Precompute expression vectors and top-10th-percentile thresholds
    all_genes = (IMMUNE_PRIMARY + IMMUNE_SUPPORTING + ENDOTHELIAL_GENES +
                 PERICYTE_GENES + RESCUE_GENES)
    marker_data = {}
    for gene in all_genes:
        expr = _get_gene_expression(adata, gene)
        if expr is not None:
            thresh = _compute_top_percentile_threshold(expr, percentile=90)
            marker_data[gene] = (expr, thresh)

    def _is_top(gene, i):
        """Check if cell i has gene expression above the top 10th percentile threshold."""
        if gene not in marker_data:
            return False
        expr, thresh = marker_data[gene]
        return expr[i] > thresh

    def _has_rescue(i):
        """Check if cell i co-expresses rescue genes (ACAN or SOX9)."""
        return _is_top('ACAN', i) or _is_top('SOX9', i)

    # Classify each cell using the hierarchy
    for i in range(n_cells):
        # ── Rule 1: Immune ──────────────────────────────────────────────
        is_immune = False
        if _is_top('PTPRC', i):
            is_immune = True
        else:
            # Check ≥2 immune supporting markers in top 10th percentile
            n_supporting = sum(1 for g in IMMUNE_SUPPORTING if _is_top(g, i))
            if n_supporting >= 2:
                is_immune = True
        if is_immune and not _has_rescue(i):
            labels[i] = "Immune"
            continue

        # ── Rule 2: Endothelial ─────────────────────────────────────────
        is_endo = any(_is_top(g, i) for g in ENDOTHELIAL_GENES)
        if is_endo and not _has_rescue(i):
            labels[i] = "Endothelial"
            continue

        # ── Rule 3: Pericyte_SMC ────────────────────────────────────────
        is_pericyte = _is_top('RGS5', i) and _is_top('PDGFRB', i)
        if is_pericyte and not _has_rescue(i):
            labels[i] = "Pericyte_SMC"
            continue

        # ── Rule 4: Chondrocyte_like ────────────────────────────────────
        s_c = score_chondro[i]
        s_f = score_fibro[i]
        if s_c > 0 and s_c > 2 * s_f:
            labels[i] = "Chondrocyte_like"
            continue

        # ── Rule 5: Fibroblast_like ─────────────────────────────────────
        if s_f > 0 and s_f > 2 * s_c:
            labels[i] = "Fibroblast_like"
            continue

        # ── Rule 6: Unknown ─────────────────────────────────────────────
        # labels[i] already set to "Unknown"

    adata.obs['coarse_label_raw'] = labels

    # Log rescue stats
    n_rescued = 0
    for i in range(n_cells):
        # Count cells that matched immune/endo/pericyte but were rescued
        if _has_rescue(i):
            is_immune = _is_top('PTPRC', i) or sum(1 for g in IMMUNE_SUPPORTING if _is_top(g, i)) >= 2
            is_endo = any(_is_top(g, i) for g in ENDOTHELIAL_GENES)
            is_pericyte = _is_top('RGS5', i) and _is_top('PDGFRB', i)
            if is_immune or is_endo or is_pericyte:
                n_rescued += 1
    if n_rescued > 0:
        print(f"    ACAN/SOX9 rescue: {n_rescued} cells kept from non-mesenchymal assignment")


def cluster_majority_vote(adata, resolution=CLASSIFICATION_RESOLUTION):
    """Assign coarse_label by cluster-level majority voting to reduce noise.

    1. Compute Leiden clusters at the given resolution
    2. For each cluster, compute the majority coarse_label
    3. If ≥85% of cells share a label, assign that label to the entire cluster
    4. Log clusters where majority is <70%
    """
    cluster_key = f"leiden_res_{resolution}"

    # Compute Leiden if not already done at this resolution
    if cluster_key not in adata.obs.columns:
        # Check if we have a neighbor graph
        if 'neighbors' not in adata.uns:
            print(f"    Computing neighbors for clustering...")
            sc.pp.neighbors(adata, n_neighbors=15)
        print(f"    Computing Leiden clusters at resolution {resolution}...")
        sc.tl.leiden(adata, resolution=resolution, key_added=cluster_key)

    n_clusters = adata.obs[cluster_key].nunique()
    print(f"    {n_clusters} clusters at resolution {resolution}")

    labels = adata.obs['coarse_label_raw'].values.copy()
    final_labels = labels.copy()
    mixed_clusters = []

    for cluster in adata.obs[cluster_key].unique():
        mask = adata.obs[cluster_key] == cluster
        cluster_labels = labels[mask]

        counts = pd.Series(cluster_labels).value_counts()
        majority_label = counts.index[0]
        majority_pct = counts.iloc[0] / len(cluster_labels) * 100

        if majority_pct >= 85:
            final_labels[mask] = majority_label
        else:
            # Mixed cluster: keep per-cell labels
            if majority_pct < 70:
                mixed_clusters.append({
                    "cluster": cluster,
                    "majority_class": majority_label,
                    "majority_pct": f"{majority_pct:.1f}",
                    "n_cells": mask.sum(),
                })

    adata.obs['coarse_label'] = final_labels

    # Derive cell_class from coarse_label
    class_map = {
        "Chondrocyte_like": "mesenchymal",
        "Fibroblast_like": "mesenchymal",
        "Immune": "non_mesenchymal",
        "Endothelial": "non_mesenchymal",
        "Pericyte_SMC": "non_mesenchymal",
        "Unknown": "unknown",
    }
    adata.obs['cell_class'] = adata.obs['coarse_label'].map(class_map)

    if mixed_clusters:
        print(f"    WARNING: {len(mixed_clusters)} clusters with <70% majority:")
        for mc in mixed_clusters:
            print(f"      Cluster {mc['cluster']}: {mc['majority_class']} "
                  f"({mc['majority_pct']}%, {mc['n_cells']} cells)")

    return mixed_clusters


# ═════════════════════════════════════════════════════════════════════════════
# VALIDATION
# ═════════════════════════════════════════════════════════════════════════════

def validate_classification(adata, accession):
    """Run classification validation checks. Returns (pass, messages)."""
    messages = []
    all_pass = True
    n_cells = adata.shape[0]

    # Check 1: All cells have coarse_label and cell_class
    if 'coarse_label' not in adata.obs.columns:
        messages.append("FAIL: coarse_label column missing")
        return False, messages
    if 'cell_class' not in adata.obs.columns:
        messages.append("FAIL: cell_class column missing")
        return False, messages

    coarse_labels = adata.obs['coarse_label']
    cell_classes = adata.obs['cell_class']

    invalid_coarse = ~coarse_labels.isin(VALID_COARSE_LABELS)
    if invalid_coarse.sum() > 0:
        messages.append(f"FAIL: {invalid_coarse.sum()} cells have invalid coarse_label values")
        all_pass = False
    else:
        messages.append(f"PASS: All {n_cells} cells have valid coarse_label values")

    invalid_class = ~cell_classes.isin(VALID_CELL_CLASSES)
    if invalid_class.sum() > 0:
        messages.append(f"FAIL: {invalid_class.sum()} cells have invalid cell_class values")
        all_pass = False
    else:
        messages.append(f"PASS: All {n_cells} cells have valid cell_class values")

    # Check 2: Unknown proportion 5-40%
    n_unknown = (coarse_labels == 'Unknown').sum()
    pct_unknown = n_unknown / n_cells * 100
    if pct_unknown > 40:
        messages.append(f"FAIL: {pct_unknown:.1f}% Unknown (> 40%, thresholds may be too strict)")
        all_pass = False
    elif pct_unknown < 5:
        messages.append(f"WARNING: {pct_unknown:.1f}% Unknown (< 5%, thresholds may be too loose)")
    else:
        messages.append(f"PASS: {pct_unknown:.1f}% Unknown (within 5-40% range)")

    # Check 3: PTPRC enriched in Immune vs Chondrocyte_like (>10-fold)
    if 'PTPRC' in adata.var_names:
        ptprc = _get_gene_expression(adata, 'PTPRC')
        immune_mask = coarse_labels == 'Immune'
        chondro_mask = coarse_labels == 'Chondrocyte_like'

        if immune_mask.sum() > 0 and chondro_mask.sum() > 0:
            mean_immune = ptprc[immune_mask].mean() + 1e-10
            mean_chondro = ptprc[chondro_mask].mean() + 1e-10
            fold = mean_immune / mean_chondro

            if fold > 10:
                messages.append(f"PASS: PTPRC {fold:.0f}-fold enriched in Immune vs Chondrocyte_like")
            else:
                messages.append(f"WARNING: PTPRC only {fold:.1f}-fold enriched in "
                                f"Immune vs Chondrocyte_like (expected >10-fold)")
        else:
            messages.append("INFO: Cannot compute PTPRC enrichment "
                            "(missing Immune or Chondrocyte_like cells)")
    else:
        messages.append("INFO: PTPRC not in var_names")

    # Check 4: COL2A1 enriched in Chondrocyte_like vs Immune (>10-fold)
    if 'COL2A1' in adata.var_names:
        col2a1 = _get_gene_expression(adata, 'COL2A1')
        chondro_mask = coarse_labels == 'Chondrocyte_like'
        immune_mask = coarse_labels == 'Immune'

        if chondro_mask.sum() > 0 and immune_mask.sum() > 0:
            mean_chondro = col2a1[chondro_mask].mean() + 1e-10
            mean_immune = col2a1[immune_mask].mean() + 1e-10
            fold = mean_chondro / mean_immune

            if fold > 10:
                messages.append(f"PASS: COL2A1 {fold:.0f}-fold enriched in Chondrocyte_like vs Immune")
            else:
                messages.append(f"WARNING: COL2A1 only {fold:.1f}-fold enriched in "
                                f"Chondrocyte_like vs Immune (expected >10-fold)")
        else:
            messages.append("INFO: Cannot compute COL2A1 enrichment "
                            "(missing Chondrocyte_like or Immune cells)")
    else:
        messages.append("INFO: COL2A1 not in var_names")

    # Check 5: COL1A1 higher in Fibroblast_like than Chondrocyte_like
    if 'COL1A1' in adata.var_names:
        col1a1 = _get_gene_expression(adata, 'COL1A1')
        fibro_mask = coarse_labels == 'Fibroblast_like'
        chondro_mask = coarse_labels == 'Chondrocyte_like'

        if fibro_mask.sum() > 0 and chondro_mask.sum() > 0:
            mean_fibro = col1a1[fibro_mask].mean() + 1e-10
            mean_chondro = col1a1[chondro_mask].mean() + 1e-10
            fold = mean_fibro / mean_chondro

            if fold > 1:
                messages.append(f"PASS: COL1A1 {fold:.1f}-fold higher in Fibroblast_like vs Chondrocyte_like")
            else:
                messages.append(f"WARNING: COL1A1 {fold:.2f}-fold in Fibroblast_like vs "
                                f"Chondrocyte_like (expected >1)")
        else:
            messages.append("INFO: Cannot compute COL1A1 enrichment "
                            "(missing Fibroblast_like or Chondrocyte_like cells)")
    else:
        messages.append("INFO: COL1A1 not in var_names")

    # Check 6: Classification report exists
    report_path = ANN_DIR / f"{accession}_classification_report.html"
    if report_path.exists():
        messages.append("PASS: Classification report HTML generated")
    else:
        messages.append("FAIL: Classification report HTML not found")
        all_pass = False

    return all_pass, messages


# ═════════════════════════════════════════════════════════════════════════════
# PLOTS
# ═════════════════════════════════════════════════════════════════════════════

def generate_classification_plots(adata, accession):
    """Generate classification-specific plots for one dataset."""
    plt.rcParams['figure.dpi'] = 100

    # 1. UMAP colored by coarse_label
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color='coarse_label', ax=ax, show=False, palette=COARSE_PALETTE,
               title=f'{accession} — Coarse Label')
    plt.tight_layout()
    plt.savefig(ANN_DIR / f"{accession}_umap_cell_class.png",
                bbox_inches='tight', dpi=150)
    plt.close()

    # 2. Dot plot of classification markers by coarse_label
    all_markers = CHONDROCYTE_GENES + FIBROBLAST_GENES + IMMUNE_PRIMARY + IMMUNE_SUPPORTING + ENDOTHELIAL_GENES + PERICYTE_GENES
    resolved_markers = [resolve_gene(g, set(adata.var_names)) for g in all_markers
                        if resolve_gene(g, set(adata.var_names)) is not None]
    # Deduplicate preserving order
    seen = set()
    unique_markers = []
    for g in resolved_markers:
        if g not in seen:
            unique_markers.append(g)
            seen.add(g)

    if unique_markers and adata.obs['coarse_label'].nunique() > 1:
        try:
            sc.pl.dotplot(
                adata, var_names=unique_markers,
                groupby='coarse_label', show=False,
            )
            plt.savefig(ANN_DIR / f"{accession}_classification_dotplot.pdf",
                        bbox_inches='tight')
            plt.close()
        except Exception as e:
            print(f"    WARNING: dotplot failed: {e}")

    # 3. Bar plot: coarse label proportions
    counts = adata.obs['coarse_label'].value_counts()
    fig, ax = plt.subplots(figsize=(8, 5))
    colors = [COARSE_PALETTE.get(c, '#999999') for c in counts.index]
    counts.plot(kind='bar', ax=ax, color=colors)
    ax.set_title(f'{accession} — Coarse Label Proportions')
    ax.set_ylabel('Number of cells')
    ax.set_xlabel('')
    for i, (idx, val) in enumerate(counts.items()):
        ax.text(i, val + 50, f'{val/adata.shape[0]*100:.1f}%', ha='center', fontsize=9)
    plt.tight_layout()
    plt.savefig(ANN_DIR / f"{accession}_class_proportions.png",
                bbox_inches='tight', dpi=150)
    plt.close()

    # 4. Validation: marker expression distributions
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    key_markers = {
        'PTPRC': 'Immune marker',
        'COL2A1': 'Chondrocyte marker',
        'COL1A1': 'Fibroblast marker',
    }
    for ax, (gene, desc) in zip(axes, key_markers.items()):
        resolved = resolve_gene(gene, set(adata.var_names))
        if resolved and adata.obs['coarse_label'].nunique() > 1:
            sc.pl.violin(adata, keys=resolved, groupby='coarse_label',
                         ax=ax, show=False, rotation=45)
            ax.set_title(f'{resolved} ({desc})')
        else:
            ax.text(0.5, 0.5, f'{gene} not available', ha='center', va='center')
            ax.set_title(gene)
    plt.tight_layout()
    plt.savefig(ANN_DIR / f"{accession}_marker_validation.png",
                bbox_inches='tight', dpi=150)
    plt.close()

    plt.close('all')


# ═════════════════════════════════════════════════════════════════════════════
# HTML REPORT
# ═════════════════════════════════════════════════════════════════════════════

CLASSIFICATION_REPORT_TEMPLATE = """<!DOCTYPE html>
<html>
<head>
<title>Classification Report: {{ accession }}</title>
<style>
  body { font-family: Arial, sans-serif; max-width: 1200px; margin: 0 auto; padding: 20px; }
  h1 { color: #2c3e50; border-bottom: 2px solid #2ecc71; padding-bottom: 10px; }
  h2 { color: #34495e; margin-top: 30px; }
  table { border-collapse: collapse; width: 100%; margin: 10px 0; }
  th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
  th { background-color: #2ecc71; color: white; }
  tr:nth-child(even) { background-color: #f2f2f2; }
  img { max-width: 100%; height: auto; margin: 10px 0; border: 1px solid #ddd; }
  .stats-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 15px; margin: 20px 0; }
  .stat-box { background: #ecf0f1; padding: 15px; border-radius: 5px; }
  .stat-box h3 { margin: 0 0 5px 0; color: #7f8c8d; font-size: 12px; }
  .stat-box p { margin: 0; font-size: 24px; font-weight: bold; color: #2c3e50; }
  .warning { color: #e74c3c; font-weight: bold; }
  .pass { color: #27ae60; }
  .info { color: #7f8c8d; }
  .validation li { margin: 4px 0; }
</style>
</head>
<body>
<h1>Classification Report: {{ accession }}</h1>
<p>Generated: {{ date }}</p>

<div class="stats-grid">
  <div class="stat-box"><h3>Total cells</h3><p>{{ stats.n_cells }}</p></div>
  <div class="stat-box"><h3>Chondrocyte-like</h3><p>{{ stats.n_chondrocyte_like }} ({{ stats.pct_chondrocyte_like }}%)</p></div>
  <div class="stat-box"><h3>Fibroblast-like</h3><p>{{ stats.n_fibroblast_like }} ({{ stats.pct_fibroblast_like }}%)</p></div>
  <div class="stat-box"><h3>Immune</h3><p>{{ stats.n_immune }} ({{ stats.pct_immune }}%)</p></div>
  <div class="stat-box"><h3>Endothelial</h3><p>{{ stats.n_endothelial }} ({{ stats.pct_endothelial }}%)</p></div>
  <div class="stat-box"><h3>Pericyte/SMC</h3><p>{{ stats.n_pericyte_smc }} ({{ stats.pct_pericyte_smc }}%)</p></div>
  <div class="stat-box"><h3>Unknown</h3><p>{{ stats.n_unknown }} ({{ stats.pct_unknown }}%)</p></div>
  <div class="stat-box"><h3>Clustering resolution</h3><p>{{ stats.resolution }}</p></div>
  <div class="stat-box"><h3>Compartment</h3><p>{{ stats.dominant_compartment }}</p></div>
</div>

<h2>UMAP — Coarse Label</h2>
<img src="{{ accession }}_umap_cell_class.png" alt="UMAP coarse label">

<h2>Coarse Label Proportions</h2>
<img src="{{ accession }}_class_proportions.png" alt="Label proportions">

<h2>Marker Dotplot</h2>
{% if has_dotplot %}
<embed src="{{ accession }}_classification_dotplot.pdf" type="application/pdf" width="100%" height="600px">
{% else %}
<p>Dotplot not generated (insufficient data).</p>
{% endif %}

<h2>Marker Expression Validation</h2>
<img src="{{ accession }}_marker_validation.png" alt="Marker validation">

{% if mixed_clusters %}
<h2>Mixed Clusters (majority < 70%)</h2>
<table>
  <tr><th>Cluster</th><th>Majority label</th><th>Majority %</th><th>Cells</th></tr>
  {% for mc in mixed_clusters %}
  <tr>
    <td>{{ mc.cluster }}</td>
    <td>{{ mc.majority_class }}</td>
    <td>{{ mc.majority_pct }}%</td>
    <td>{{ mc.n_cells }}</td>
  </tr>
  {% endfor %}
</table>
{% endif %}

<h2>Marker Availability</h2>
<table>
  <tr><th>Marker set</th><th>Available / Total</th><th>Missing genes</th></tr>
  {% for row in marker_table %}
  <tr>
    <td>{{ row.name }}</td>
    <td>{{ row.n_avail }} / {{ row.n_total }}</td>
    <td>{{ row.missing }}</td>
  </tr>
  {% endfor %}
</table>

<h2>Validation</h2>
<ul class="validation">
  {% for msg in validation_messages %}
  <li class="{{ 'pass' if msg.startswith('PASS') else 'warning' if msg.startswith('WARNING') or msg.startswith('FAIL') else 'info' }}">{{ msg }}</li>
  {% endfor %}
</ul>

</body>
</html>
"""


def generate_classification_report(accession, stats, marker_stats,
                                   mixed_clusters, validation_messages):
    """Generate HTML classification report for one dataset."""
    marker_table = []
    for name, info in marker_stats.items():
        marker_table.append({
            "name": name,
            "n_avail": len(info["available"]),
            "n_total": len(info["available"]) + len(info["missing"]),
            "missing": ", ".join(info["missing"]) if info["missing"] else "—",
        })

    template = Template(CLASSIFICATION_REPORT_TEMPLATE)
    html = template.render(
        accession=accession,
        date=datetime.now().strftime("%Y-%m-%d %H:%M"),
        stats=stats,
        marker_table=marker_table,
        mixed_clusters=mixed_clusters,
        validation_messages=validation_messages,
        has_dotplot=(ANN_DIR / f"{accession}_classification_dotplot.pdf").exists(),
    )
    report_path = ANN_DIR / f"{accession}_classification_report.html"
    report_path.write_text(html)
    print(f"    Report: {report_path}")


# ═════════════════════════════════════════════════════════════════════════════
# AGGREGATE OUTPUTS
# ═════════════════════════════════════════════════════════════════════════════

def generate_aggregate_outputs(all_stats):
    """Generate classification_summary.tsv with coarse label proportions per dataset."""
    rows = []
    for stats in all_stats:
        rows.append({
            "accession": stats["accession"],
            "n_cells": stats["n_cells"],
            "n_chondrocyte_like": stats["n_chondrocyte_like"],
            "pct_chondrocyte_like": stats["pct_chondrocyte_like"],
            "n_fibroblast_like": stats["n_fibroblast_like"],
            "pct_fibroblast_like": stats["pct_fibroblast_like"],
            "n_immune": stats["n_immune"],
            "pct_immune": stats["pct_immune"],
            "n_endothelial": stats["n_endothelial"],
            "pct_endothelial": stats["pct_endothelial"],
            "n_pericyte_smc": stats["n_pericyte_smc"],
            "pct_pericyte_smc": stats["pct_pericyte_smc"],
            "n_unknown": stats["n_unknown"],
            "pct_unknown": stats["pct_unknown"],
            "dominant_compartment": stats["dominant_compartment"],
        })

    df = pd.DataFrame(rows)
    tsv_path = ANN_DIR / "classification_summary.tsv"
    df.to_csv(tsv_path, sep='\t', index=False)
    print(f"\nClassification summary: {tsv_path}")


# ═════════════════════════════════════════════════════════════════════════════
# PER-DATASET PROCESSING
# ═════════════════════════════════════════════════════════════════════════════

def classify_dataset(accession, excluded_samples):
    """Full classification pipeline for one dataset.

    Returns (stats_dict, validation_messages) or None on error.
    """
    h5ad_path = PROC_DIR / f"{accession}.h5ad"
    if not h5ad_path.exists():
        print(f"  SKIP {accession}: {h5ad_path} not found")
        return None

    print(f"\n{'='*60}")
    print(f"Classifying {accession}")
    print(f"{'='*60}")

    adata = sc.read_h5ad(h5ad_path)
    print(f"  Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # ── Filter excluded samples ──────────────────────────────────────────
    adata = filter_excluded_samples(adata, accession, excluded_samples)
    print(f"  After filtering: {adata.shape[0]} cells")

    # ── Step 1: Score marker sets ────────────────────────────────────────
    print(f"\n  Step 1: Scoring marker gene sets...")
    marker_stats = score_cell_class(adata)

    # ── Step 2: Cell-level classification ────────────────────────────────
    print(f"\n  Step 2: Cell-level classification...")
    classify_cells(adata)
    raw_counts = adata.obs['coarse_label_raw'].value_counts()
    for cls, n in raw_counts.items():
        print(f"    {cls}: {n} ({n/adata.shape[0]*100:.1f}%)")

    # ── Step 3: Cluster majority voting ──────────────────────────────────
    print(f"\n  Step 3: Cluster-level majority voting...")
    mixed_clusters = cluster_majority_vote(adata)
    final_counts = adata.obs['coarse_label'].value_counts()
    print(f"    Final coarse_label distribution:")
    for cls, n in final_counts.items():
        print(f"      {cls}: {n} ({n/adata.shape[0]*100:.1f}%)")
    class_counts = adata.obs['cell_class'].value_counts()
    print(f"    Derived cell_class distribution:")
    for cls, n in class_counts.items():
        print(f"      {cls}: {n} ({n/adata.shape[0]*100:.1f}%)")

    # ── Step 4: Generate plots ───────────────────────────────────────────
    print(f"\n  Step 4: Generating classification plots...")
    generate_classification_plots(adata, accession)

    # ── Step 5: Validation ───────────────────────────────────────────────
    print(f"\n  Step 5: Running validation checks...")
    passed, val_messages = validate_classification(adata, accession)
    for msg in val_messages:
        print(f"    {msg}")

    # ── Collect stats ────────────────────────────────────────────────────
    dominant_comp = "mixed"
    if 'compartment' in adata.obs.columns:
        dominant_comp = adata.obs['compartment'].value_counts().index[0]

    n_total = adata.shape[0]
    n_chondro = int((adata.obs['coarse_label'] == 'Chondrocyte_like').sum())
    n_fibro = int((adata.obs['coarse_label'] == 'Fibroblast_like').sum())
    n_immune = int((adata.obs['coarse_label'] == 'Immune').sum())
    n_endo = int((adata.obs['coarse_label'] == 'Endothelial').sum())
    n_peri = int((adata.obs['coarse_label'] == 'Pericyte_SMC').sum())
    n_unknown = int((adata.obs['coarse_label'] == 'Unknown').sum())

    stats = {
        "accession": accession,
        "n_cells": n_total,
        "n_chondrocyte_like": n_chondro,
        "pct_chondrocyte_like": f"{n_chondro / n_total * 100:.1f}",
        "n_fibroblast_like": n_fibro,
        "pct_fibroblast_like": f"{n_fibro / n_total * 100:.1f}",
        "n_immune": n_immune,
        "pct_immune": f"{n_immune / n_total * 100:.1f}",
        "n_endothelial": n_endo,
        "pct_endothelial": f"{n_endo / n_total * 100:.1f}",
        "n_pericyte_smc": n_peri,
        "pct_pericyte_smc": f"{n_peri / n_total * 100:.1f}",
        "n_unknown": n_unknown,
        "pct_unknown": f"{n_unknown / n_total * 100:.1f}",
        "resolution": CLASSIFICATION_RESOLUTION,
        "dominant_compartment": dominant_comp,
    }

    # ── Step 6: Generate report ──────────────────────────────────────────
    print(f"\n  Step 6: Generating classification report...")
    generate_classification_report(accession, stats, marker_stats,
                                   mixed_clusters, val_messages)

    # ── Step 7: Save updated h5ad ────────────────────────────────────────
    print(f"\n  Step 7: Saving updated h5ad...")
    # Clean up old annotation columns if present (from previous pipeline run)
    old_cols = ['cell_type_final', 'cell_type_marker_based', 'cell_type_celltypist',
                'cell_type_reference_based', 'cell_type_confidence',
                'cell_type_marker_confidence', 'celltypist_conf_score',
                'cell_class_raw', 'score_mesenchymal', 'score_non_mesenchymal']
    for col in old_cols:
        if col in adata.obs.columns:
            del adata.obs[col]
    # Also remove old score columns (but keep our new ones)
    keep_scores = {'score_chondrocyte', 'score_fibroblast'}
    old_score_cols = [c for c in adata.obs.columns
                      if c.startswith('score_') and c not in keep_scores]
    for col in old_score_cols:
        del adata.obs[col]

    adata.write_h5ad(h5ad_path)
    print(f"    Saved: {h5ad_path}")

    del adata
    return stats, val_messages


# ═════════════════════════════════════════════════════════════════════════════
# VALIDATION-ONLY MODE
# ═════════════════════════════════════════════════════════════════════════════

def validate_all():
    """Run validation across all datasets."""
    print("\n" + "=" * 60)
    print("CLASSIFICATION VALIDATION")
    print("=" * 60)

    all_pass = True
    for accession in ALL_ACCESSIONS:
        h5ad_path = PROC_DIR / f"{accession}.h5ad"
        if not h5ad_path.exists():
            print(f"\n--- {accession} ---")
            print(f"  SKIP: {h5ad_path} not found")
            continue

        print(f"\n--- {accession} ---")
        adata = sc.read_h5ad(h5ad_path)

        if 'coarse_label' not in adata.obs.columns:
            print(f"  SKIP: coarse_label not present (not yet classified)")
            del adata
            continue

        passed, messages = validate_classification(adata, accession)
        for msg in messages:
            print(f"  {msg}")
        if not passed:
            all_pass = False
        del adata

    # Cross-dataset checks
    print(f"\n--- Cross-dataset ---")
    tsv_path = ANN_DIR / "classification_summary.tsv"
    if tsv_path.exists():
        summary_df = pd.read_csv(tsv_path, sep='\t')
        print(f"  PASS: classification_summary.tsv has {len(summary_df)} rows")
    else:
        print(f"  FAIL: classification_summary.tsv not found")
        all_pass = False

    n_reports = sum(1 for a in ALL_ACCESSIONS
                    if (ANN_DIR / f"{a}_classification_report.html").exists())
    print(f"  Reports: {n_reports}/{len(ALL_ACCESSIONS)} datasets")

    print(f"\n{'='*60}")
    print(f"OVERALL: {'ALL CHECKS PASSED' if all_pass else 'SOME CHECKS FAILED'}")
    print(f"{'='*60}")
    return all_pass


# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════

def main():
    args = sys.argv[1:]

    if args and args[0] == '--validate-only':
        validate_all()
        return

    # Load excluded samples
    excluded_samples = get_excluded_samples()
    if excluded_samples:
        print(f"Excluded samples (neonatal): {excluded_samples}")

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
        try:
            result = classify_dataset(accession, excluded_samples)
            if result is not None:
                stats, _ = result
                all_stats.append(stats)
        except Exception as e:
            print(f"\nERROR classifying {accession}: {e}")
            import traceback
            traceback.print_exc()
            continue

    # Generate aggregate outputs
    if all_stats:
        print(f"\n{'='*60}")
        print("AGGREGATE OUTPUTS")
        print(f"{'='*60}")
        generate_aggregate_outputs(all_stats)

    # Run validation
    validate_all()


if __name__ == "__main__":
    main()
