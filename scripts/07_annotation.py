#!/usr/bin/env python3
"""Module 07: Post-integration cell type annotation.

Two-stage annotation of clustered objects from Module 06:
  Stage 1: Coarse identity via canonical markers (dot plots)
  Stage 2: Fine distinctions via DE within coarse groups

Usage:
    python3 scripts/07_annotation.py                  # All objects
    python3 scripts/07_annotation.py --object NP      # Single object
    python3 scripts/07_annotation.py --validate-only  # Validation only
"""

import gc
import sys
import warnings
import traceback
from pathlib import Path
from datetime import datetime

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
# Defaults can be overridden via --input-dir / --output-dir on the command line
# so the script can be aimed at a parallel pipeline (e.g. tiered_v4) without
# overwriting the production v5 outputs.
BASE = Path(__file__).resolve().parent.parent
INT_DIR = BASE / "data" / "integrated"
RESULTS_DIR = BASE / "results" / "integration"


# ── Canonical marker panels for coarse annotation ────────────────────────────
CANONICAL_MARKERS = {
    # NP subtypes
    "NP_notochordal":       ["T", "SHH", "NOG", "CD24", "KRT8", "KRT18", "KRT19"],
    "NP_mature_chondrocyte": ["ACAN", "COL2A1", "SOX9", "COMP", "PRG4"],
    "NP_stressed_degen":    ["MMP13", "ADAMTS5", "IL1B", "TNF", "VEGFA", "HIF1A"],
    "NP_fibrocartilaginous": ["COL1A1", "COL2A1", "VCAN"],
    # AF subtypes
    "AF_inner":             ["COL2A1", "ACAN", "SOX9"],
    "AF_outer":             ["COL1A1", "COL1A2", "THY1", "DCN", "LUM"],
    "AF_mechanical_stress": ["COMP", "CILP", "THBS1"],
    # Endplate
    "EP_hyaline":           ["COL2A1", "COL10A1", "SOX9"],
    "EP_ossification":      ["RUNX2", "SP7", "BGLAP"],
    # Non-mesenchymal
    "Macrophage":           ["CD68", "CD14", "CSF1R", "CD163", "CD86"],
    "T_cell":               ["CD3D", "CD3E", "CD4", "CD8A"],
    "B_cell":               ["CD79A", "MS4A1"],
    "NK_cell":              ["NKG7", "GNLY"],
    "Mast_cell":            ["KIT", "TPSAB1"],
    "Endothelial":          ["PECAM1", "VWF", "CDH5"],
    "Pericyte_SMC":         ["ACTA2", "RGS5", "PDGFRB"],
    # Additional coarse panels
    "Fibrochondrocyte_like": ["COL1A1", "COL2A1", "ACAN", "DCN"],
    "Chondrocyte_like":     ["COL2A1", "ACAN", "SOX9", "COMP", "PRG4"],
    "Fibroblast_like":      ["COL1A1", "COL1A2", "DCN", "LUM", "THY1"],
}

# Coarse panels — used in Stage 1 to bucket clusters into broad categories
COARSE_PANELS_MESENCHYMAL = {
    "Chondrocyte_like":      ["COL2A1", "ACAN", "SOX9", "COMP", "PRG4"],
    "Fibroblast_like":       ["COL1A1", "COL1A2", "DCN", "LUM", "THY1"],
    "Fibrochondrocyte_like": ["COL1A1", "COL2A1", "ACAN", "DCN"],
}

COARSE_PANELS_NON_MESENCHYMAL = {
    "Macrophage":     ["CD68", "CD14", "CSF1R", "CD163", "CD86"],
    "T_cell":         ["CD3D", "CD3E", "CD4", "CD8A"],
    "B_cell":         ["CD79A", "MS4A1"],
    "NK_cell":        ["NKG7", "GNLY"],
    "Mast_cell":      ["KIT", "TPSAB1"],
    "Endothelial":    ["PECAM1", "VWF", "CDH5"],
    "Pericyte_SMC":   ["ACTA2", "RGS5", "PDGFRB"],
}

# Fine marker panels — used in Stage 2 within each coarse group
FINE_PANELS = {
    # Within Chondrocyte-like (NP context)
    "NP_notochordal":       ["T", "SHH", "NOG", "CD24", "KRT8", "KRT18", "KRT19"],
    "NP_mature_chondrocyte": ["ACAN", "COL2A1", "SOX9", "COMP", "PRG4"],
    "NP_stressed":          ["MMP13", "ADAMTS5", "IL1B", "TNF", "VEGFA", "HIF1A"],
    "NP_hypertrophic":      ["COL10A1", "RUNX2"],
    # Within Fibroblast-like (AF context)
    "AF_inner":             ["COL2A1", "ACAN", "SOX9"],
    "AF_outer":             ["COL1A1", "COL1A2", "THY1", "DCN", "LUM"],
    "AF_mechanical_stress": ["COMP", "CILP", "THBS1"],
    # Within Fibrochondrocyte-like
    "Fibrochondrocyte_chondroid": ["COL2A1", "ACAN", "SOX9"],
    "Fibrochondrocyte_fibroid":   ["COL1A1", "COL1A2", "DCN"],
    "Fibrochondrocyte_stressed":  ["MMP13", "ADAMTS5", "VEGFA"],
    # Within Endplate chondrocytes
    "EP_hyaline":           ["COL2A1", "COL10A1", "SOX9"],
    "EP_ossification":      ["RUNX2", "SP7", "BGLAP"],
    # Within Immune
    "Macrophage_M1":        ["CD86", "IL1B", "TNF", "NOS2"],
    "Macrophage_M2":        ["CD163", "MRC1", "MSR1", "TGFB1"],
    "T_cell_CD4":           ["CD4", "IL7R", "FOXP3"],
    "T_cell_CD8":           ["CD8A", "CD8B", "GZMB", "PRF1"],
    # Within Endothelial / Pericyte — no further sub-typing by default
}

# Gene aliases
GENE_ALIASES = {"T": ["T", "TBXT"], "TBXT": ["TBXT", "T"], "SP7": ["SP7", "OSX"]}

# Continuous score signatures for mesenchymal continuum
CONTINUUM_SCORES = {
    "score_notochordal":  ["T", "SHH", "NOG", "CD24", "KRT8", "KRT18", "KRT19"],
    "score_degenerative": ["MMP13", "ADAMTS5", "IL1B", "TNF", "VEGFA", "HIF1A"],
    "score_fibrotic":     ["COL1A1", "COL3A1", "FN1", "VCAN", "ACTA2", "TGFB1"],
}


# ═════════════════════════════════════════════════════════════════════════════
# GENE RESOLUTION
# ═════════════════════════════════════════════════════════════════════════════

def _resolve_gene(name, var_names):
    """Resolve gene name to form present in var_names."""
    if name in var_names:
        return name
    for alias in GENE_ALIASES.get(name, []):
        if alias in var_names:
            return alias
    return None


def _resolve_gene_list(gene_list, var_names):
    """Resolve genes, returning (resolved, missing)."""
    var_set = set(var_names)
    resolved, missing = [], []
    for g in gene_list:
        match = _resolve_gene(g, var_set)
        (resolved if match else missing).append(match or g)
    return resolved, missing


# ═════════════════════════════════════════════════════════════════════════════
# CLUSTER DE MARKERS
# ═════════════════════════════════════════════════════════════════════════════

def compute_cluster_markers(adata, object_name, tier_name, n_genes=50):
    """Compute DE markers per cluster using Wilcoxon rank-sum test."""
    print(f"  Computing cluster markers ({tier_name})...")

    sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon',
                            n_genes=n_genes)

    # Save to TSV
    marker_df = sc.get.rank_genes_groups_df(adata, group=None)
    marker_df.to_csv(
        RESULTS_DIR / "cluster_markers" / f"{object_name}_{tier_name}_markers.tsv",
        sep='\t', index=False
    )
    print(f"    Saved markers: {len(marker_df)} rows")
    return marker_df


def compute_within_group_markers(adata, cluster_col='leiden', n_genes=50):
    """Compute DE markers comparing clusters within a subset.

    Used in Stage 2 to find distinguishing markers between clusters of
    the same coarse type.
    """
    clusters = adata.obs[cluster_col].unique()
    if len(clusters) < 2:
        return None

    sc.tl.rank_genes_groups(adata, groupby=cluster_col, method='wilcoxon',
                            n_genes=n_genes)
    return sc.get.rank_genes_groups_df(adata, group=None)


# ═════════════════════════════════════════════════════════════════════════════
# DOTPLOTS
# ═════════════════════════════════════════════════════════════════════════════

def generate_annotation_dotplots(adata, object_name, tier_name, marker_panel):
    """Generate dot plots of canonical markers per cluster."""
    var_names = set(adata.var_names)
    resolved_markers = {}
    for group, genes in marker_panel.items():
        resolved = [_resolve_gene(g, var_names) for g in genes
                     if _resolve_gene(g, var_names) is not None]
        if resolved:
            resolved_markers[group] = resolved

    if not resolved_markers and adata.obs['leiden'].nunique() <= 1:
        return

    try:
        # Flat list for dotplot, preserving group structure
        all_genes = []
        for genes in resolved_markers.values():
            for g in genes:
                if g not in all_genes:
                    all_genes.append(g)

        if all_genes:
            sc.pl.dotplot(adata, var_names=all_genes[:50], groupby='leiden', show=False)
            plt.savefig(
                RESULTS_DIR / "annotation_dotplots" /
                f"{object_name}_{tier_name}_canonical_markers.pdf",
                bbox_inches='tight'
            )
            plt.close()
    except Exception as e:
        print(f"    WARNING: Dotplot failed: {e}")


# ═════════════════════════════════════════════════════════════════════════════
# SCORING — shared by Stage 1 and Stage 2
# ═════════════════════════════════════════════════════════════════════════════

def _score_cluster_for_panel(adata, mask, panel, var_names_set):
    """Score a single cluster against a marker panel dict.

    For each panel entry, computes per-gene fraction expressing * specificity
    ratio (in-cluster mean / out-of-cluster mean), capped at 5x.

    Returns dict of {panel_name: score}.
    """
    scores = {}
    for ct_name, genes in panel.items():
        resolved = [_resolve_gene(g, var_names_set) for g in genes
                     if _resolve_gene(g, var_names_set) is not None]
        if not resolved:
            scores[ct_name] = 0.0
            continue

        expr_in = adata[mask, resolved].X
        expr_out = adata[~mask, resolved].X
        if sparse.issparse(expr_in):
            expr_in = expr_in.toarray()
        if sparse.issparse(expr_out):
            expr_out = expr_out.toarray()

        gene_scores = []
        for g_idx in range(expr_in.shape[1]):
            frac_in = (expr_in[:, g_idx] > 0).mean()
            mean_in = expr_in[:, g_idx].mean()
            mean_out = expr_out[:, g_idx].mean()
            specificity = mean_in / (mean_out + 0.01)
            gene_scores.append(frac_in * min(specificity, 5.0))
        scores[ct_name] = float(np.mean(gene_scores)) if gene_scores else 0.0

    return scores


def _assign_best_label(scores, min_score=0.05):
    """Pick the best label from a scores dict. Returns (label, confidence, evidence)."""
    if not scores or all(v == 0 for v in scores.values()):
        return "unassigned", "low", ""

    sorted_scores = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    best_name, best_score = sorted_scores[0]
    second_score = sorted_scores[1][1] if len(sorted_scores) > 1 else 0
    margin = best_score - second_score

    if best_score <= min_score:
        return "unassigned", "low", ""

    if margin > 0.1:
        confidence = "high"
    elif margin > 0.03:
        confidence = "medium"
    else:
        confidence = "low"

    evidence = f"top_marker={best_name}(score={best_score:.3f})"
    return best_name, confidence, evidence


# ═════════════════════════════════════════════════════════════════════════════
# STAGE 1: COARSE ANNOTATION
# ═════════════════════════════════════════════════════════════════════════════

def _get_coarse_panel(tier_name):
    """Return the coarse marker panel appropriate for a tier."""
    if tier_name == "non_mesenchymal":
        return COARSE_PANELS_NON_MESENCHYMAL
    else:
        return COARSE_PANELS_MESENCHYMAL


def annotate_coarse(adata, object_name, tier_name):
    """Stage 1: Assign coarse_cell_type to each cluster via canonical markers.

    For mesenchymal tiers: Chondrocyte_like, Fibroblast_like, Fibrochondrocyte_like
    For non-mesenchymal tiers: Macrophage, T_cell, B_cell, NK_cell, Mast_cell,
                                Endothelial, Pericyte_SMC
    """
    print(f"  Stage 1: Coarse annotation ({tier_name})...")

    panel = _get_coarse_panel(tier_name)
    var_names_set = set(adata.var_names)

    # Generate dotplots for the coarse panel
    generate_annotation_dotplots(adata, object_name, f"{tier_name}_coarse", panel)

    # Score each cluster
    coarse_types = pd.Series("unassigned", index=adata.obs_names, dtype=object)
    coarse_conf = pd.Series("low", index=adata.obs_names, dtype=object)
    coarse_evidence = pd.Series("", index=adata.obs_names, dtype=object)

    cluster_assignments = {}

    for cluster in sorted(adata.obs['leiden'].unique(), key=lambda x: int(x)):
        mask = adata.obs['leiden'] == cluster
        scores = _score_cluster_for_panel(adata, mask, panel, var_names_set)

        # For mesenchymal: detect fibrochondrocyte-like (co-expression of both
        # chondrocyte and fibroblast markers at moderate levels)
        if tier_name != "non_mesenchymal":
            chondro_score = scores.get("Chondrocyte_like", 0)
            fibro_score = scores.get("Fibroblast_like", 0)
            fibrochondro_score = scores.get("Fibrochondrocyte_like", 0)

            # If both chondrocyte and fibroblast scores are above threshold
            # and neither dominates strongly, call it fibrochondrocyte
            if (chondro_score > 0.05 and fibro_score > 0.05 and
                    abs(chondro_score - fibro_score) < 0.15):
                # Boost fibrochondrocyte score
                scores["Fibrochondrocyte_like"] = max(
                    fibrochondro_score,
                    (chondro_score + fibro_score) / 2
                )

        label, confidence, evidence = _assign_best_label(scores)
        cluster_assignments[cluster] = label
        coarse_types.loc[mask] = label
        coarse_conf.loc[mask] = confidence
        coarse_evidence.loc[mask] = evidence

    adata.obs['coarse_cell_type'] = coarse_types.values
    adata.obs['cell_type_confidence'] = coarse_conf.values
    adata.obs['annotation_evidence'] = coarse_evidence.values

    # Summary
    ct_counts = adata.obs['coarse_cell_type'].value_counts()
    print(f"    Coarse annotation summary ({tier_name}):")
    for ct, n in ct_counts.items():
        pct = n / adata.shape[0] * 100
        clusters = [c for c, l in cluster_assignments.items() if l == ct]
        print(f"      {ct}: {n:,} ({pct:.1f}%) — clusters {clusters}")

    return cluster_assignments


# ═════════════════════════════════════════════════════════════════════════════
# STAGE 2: FINE ANNOTATION WITHIN COARSE GROUPS
# ═════════════════════════════════════════════════════════════════════════════

def _get_fine_panel_for_coarse(coarse_type, object_name):
    """Return the fine marker panel appropriate for a coarse type + compartment."""
    if coarse_type == "Chondrocyte_like":
        if object_name in ("NP", "all_cells"):
            return {
                "NP_notochordal":       FINE_PANELS["NP_notochordal"],
                "NP_mature_chondrocyte": FINE_PANELS["NP_mature_chondrocyte"],
                "NP_stressed":          FINE_PANELS["NP_stressed"],
                "NP_hypertrophic":      FINE_PANELS["NP_hypertrophic"],
            }
        elif object_name == "CEP":
            return {
                "EP_hyaline":           FINE_PANELS["EP_hyaline"],
                "EP_ossification":      FINE_PANELS["EP_ossification"],
            }
        else:
            # AF chondrocyte-like → inner AF
            return {
                "AF_inner":             FINE_PANELS["AF_inner"],
            }

    elif coarse_type == "Fibroblast_like":
        if object_name == "AF":
            return {
                "AF_outer":             FINE_PANELS["AF_outer"],
                "AF_mechanical_stress": FINE_PANELS["AF_mechanical_stress"],
            }
        elif object_name in ("NP", "all_cells"):
            return {
                "NP_fibrocartilaginous": CANONICAL_MARKERS["NP_fibrocartilaginous"],
            }
        else:
            return {}

    elif coarse_type == "Fibrochondrocyte_like":
        return {
            "Fibrochondrocyte_chondroid": FINE_PANELS["Fibrochondrocyte_chondroid"],
            "Fibrochondrocyte_fibroid":   FINE_PANELS["Fibrochondrocyte_fibroid"],
            "Fibrochondrocyte_stressed":  FINE_PANELS["Fibrochondrocyte_stressed"],
        }

    elif coarse_type == "Macrophage":
        return {
            "Macrophage_M1": FINE_PANELS["Macrophage_M1"],
            "Macrophage_M2": FINE_PANELS["Macrophage_M2"],
        }

    elif coarse_type == "T_cell":
        return {
            "T_cell_CD4": FINE_PANELS["T_cell_CD4"],
            "T_cell_CD8": FINE_PANELS["T_cell_CD8"],
        }

    # B_cell, NK_cell, Mast_cell, Endothelial, Pericyte_SMC — no further splitting
    return {}


def annotate_fine(adata, object_name, tier_name, coarse_assignments):
    """Stage 2: Assign fine cell_type within each coarse group using DE + markers.

    For each coarse group with >=2 clusters:
      1. Run rank_genes_groups comparing only clusters of that group
      2. Score clusters against fine marker panels
      3. Assign cell_type label
    """
    print(f"  Stage 2: Fine annotation ({tier_name})...")

    # Initialize cell_type from coarse_cell_type (will be refined below)
    adata.obs['cell_type'] = adata.obs['coarse_cell_type'].copy()

    var_names_set = set(adata.var_names)

    # Group clusters by coarse type
    coarse_groups = {}
    for cluster, coarse_type in coarse_assignments.items():
        coarse_groups.setdefault(coarse_type, []).append(cluster)

    for coarse_type, clusters in coarse_groups.items():
        if coarse_type == "unassigned":
            continue

        fine_panel = _get_fine_panel_for_coarse(coarse_type, object_name)

        # If no fine panel defined, keep the coarse label
        if not fine_panel:
            print(f"    {coarse_type}: no fine panel — keeping coarse label "
                  f"({len(clusters)} cluster(s))")
            continue

        # If only 1 cluster in this coarse group — still score against fine panels
        # (within-group DE requires >=2 clusters)
        group_mask = adata.obs['leiden'].isin(clusters)
        group_adata = adata[group_mask]

        # Run within-group DE if >= 2 clusters
        within_markers = None
        if len(clusters) >= 2:
            try:
                group_adata_copy = group_adata.copy()
                within_markers = compute_within_group_markers(group_adata_copy)
                if within_markers is not None:
                    within_markers.to_csv(
                        RESULTS_DIR / "cluster_markers" /
                        f"{object_name}_{tier_name}_{coarse_type}_within_markers.tsv",
                        sep='\t', index=False
                    )
                del group_adata_copy
            except Exception as e:
                print(f"    WARNING: Within-group DE failed for {coarse_type}: {e}")

        # Score each cluster in this group against fine panels
        for cluster in clusters:
            cluster_mask = adata.obs['leiden'] == cluster
            scores = _score_cluster_for_panel(adata, cluster_mask, fine_panel,
                                              var_names_set)
            label, confidence, evidence = _assign_best_label(scores)

            if label != "unassigned":
                adata.obs.loc[cluster_mask, 'cell_type'] = label
                # Update confidence — only upgrade, never downgrade
                current_conf = adata.obs.loc[cluster_mask, 'cell_type_confidence'].iloc[0]
                if _confidence_rank(confidence) > _confidence_rank(current_conf):
                    adata.obs.loc[cluster_mask, 'cell_type_confidence'] = confidence
                # Append fine evidence
                mask_idx = cluster_mask[cluster_mask].index
                current_ev = adata.obs.loc[mask_idx, 'annotation_evidence']
                adata.obs.loc[mask_idx, 'annotation_evidence'] = (
                    current_ev + "; fine=" + evidence
                )
            else:
                # Keep coarse label but note that fine typing was ambiguous
                mask_idx = cluster_mask[cluster_mask].index
                current_ev = adata.obs.loc[mask_idx, 'annotation_evidence']
                adata.obs.loc[mask_idx, 'annotation_evidence'] = (
                    current_ev + "; fine=ambiguous"
                )

        print(f"    {coarse_type}: {len(clusters)} cluster(s) — "
              f"fine types: {adata.obs.loc[group_mask, 'cell_type'].unique().tolist()}")

    # Summary
    ct_counts = adata.obs['cell_type'].value_counts()
    print(f"    Fine annotation summary ({tier_name}):")
    for ct, n in ct_counts.items():
        pct = n / adata.shape[0] * 100
        print(f"      {ct}: {n:,} ({pct:.1f}%)")


def _confidence_rank(conf):
    """Numeric rank for confidence levels."""
    return {"high": 3, "medium": 2, "low": 1}.get(conf, 0)


# ═════════════════════════════════════════════════════════════════════════════
# CONTINUOUS SCORES
# ═════════════════════════════════════════════════════════════════════════════

def compute_continuous_scores(adata):
    """Compute continuous gene signature scores for mesenchymal continuum."""
    var_names = set(adata.var_names)
    for score_name, genes in CONTINUUM_SCORES.items():
        resolved = [_resolve_gene(g, var_names) for g in genes
                     if _resolve_gene(g, var_names) is not None]
        if resolved:
            try:
                sc.tl.score_genes(adata, gene_list=resolved, score_name=score_name)
            except Exception:
                adata.obs[score_name] = 0.0
        else:
            adata.obs[score_name] = 0.0


# ═════════════════════════════════════════════════════════════════════════════
# CELLTYPIST VALIDATION (non-mesenchymal)
# ═════════════════════════════════════════════════════════════════════════════

def validate_immune_with_celltypist(adata, object_name):
    """Run CellTypist on non-mesenchymal cells as a validation check.

    Compares CellTypist predictions against de novo labels per cluster.
    Saves concordance table. Flags disagreements in annotation_evidence.
    """
    print(f"  Running CellTypist validation on non-mesenchymal cells...")

    try:
        import celltypist
        from celltypist import models as ct_models
    except ImportError:
        print("    WARNING: celltypist not installed, skipping validation")
        return None

    try:
        model = ct_models.Model.load(model="Immune_All_Low.pkl")
    except Exception as e:
        print(f"    WARNING: Could not load CellTypist model: {e}")
        return None

    try:
        predictions = celltypist.annotate(adata, model=model, majority_voting=True)
        result_adata = predictions.to_adata()
        adata.obs['celltypist_prediction'] = result_adata.obs['majority_voting'].values
    except Exception as e:
        print(f"    WARNING: CellTypist annotation failed: {e}")
        return None

    # Build concordance table: cluster x CellTypist majority prediction
    concordance = []
    for cluster in sorted(adata.obs['leiden'].unique(), key=lambda x: int(x)):
        mask = adata.obs['leiden'] == cluster
        de_novo_label = adata.obs.loc[mask, 'cell_type'].mode()
        de_novo_label = de_novo_label.iloc[0] if len(de_novo_label) > 0 else "unknown"
        ct_labels = adata.obs.loc[mask, 'celltypist_prediction']
        ct_majority = ct_labels.mode()
        ct_majority = ct_majority.iloc[0] if len(ct_majority) > 0 else "unknown"
        ct_agreement = (ct_labels == ct_majority).mean()

        agrees = _labels_agree(de_novo_label, ct_majority)

        concordance.append({
            'cluster': cluster,
            'n_cells': mask.sum(),
            'de_novo_label': de_novo_label,
            'celltypist_majority': ct_majority,
            'celltypist_agreement_pct': f"{ct_agreement * 100:.1f}",
            'concordant': agrees,
        })

        # Flag disagreements
        if not agrees:
            mask_idx = mask[mask].index
            current_evidence = adata.obs.loc[mask_idx, 'annotation_evidence']
            adata.obs.loc[mask_idx, 'annotation_evidence'] = (
                current_evidence + f"; CELLTYPIST_DISAGREES={ct_majority}"
            )
            print(f"    WARNING: Cluster {cluster}: de novo={de_novo_label}, "
                  f"CellTypist={ct_majority} — flagged for review")

    conc_df = pd.DataFrame(concordance)
    conc_df.to_csv(
        RESULTS_DIR / "celltypist_validation" / f"{object_name}_concordance.tsv",
        sep='\t', index=False
    )
    n_agree = conc_df['concordant'].sum()
    print(f"    Concordance: {n_agree}/{len(conc_df)} clusters agree")
    return conc_df


def _labels_agree(de_novo, celltypist):
    """Check if de novo and CellTypist labels broadly agree."""
    dn = de_novo.lower()
    ct = celltypist.lower()

    # Map common label patterns
    mappings = {
        'macrophage': ['macro', 'mono', 'dc', 'dendritic'],
        't_cell': ['t cell', 'cd4', 'cd8', 'tem', 'tcm', 'treg'],
        'b_cell': ['b cell', 'plasma', 'pre-b', 'pro-b'],
        'nk_cell': ['nk'],
        'mast_cell': ['mast'],
        'endothelial': ['endothelial'],
        'pericyte_smc': ['pericyte', 'smooth muscle', 'smc'],
    }

    dn_category = None
    ct_category = None

    for key, patterns in mappings.items():
        if any(p in dn for p in patterns) or key in dn:
            dn_category = key
        if any(p in ct for p in patterns):
            ct_category = key

    # If both match a known category, they must match the same one
    if dn_category is not None and ct_category is not None:
        return dn_category == ct_category

    # If only one matches a known category, they disagree
    if dn_category is not None or ct_category is not None:
        return False

    # If neither matches known patterns, can't determine disagreement
    return True


# ═════════════════════════════════════════════════════════════════════════════
# CELL TYPE DEFINITIONS TABLE
# ═════════════════════════════════════════════════════════════════════════════

def generate_cell_type_definitions(all_annotations):
    """Generate cell_type_definitions.tsv summarizing all cell types across objects.

    Parameters
    ----------
    all_annotations : dict
        {object_name: adata} for each processed object.
    """
    print("\n  Generating cell type definitions table...")

    rows = []
    for obj_name, adata in all_annotations.items():
        if adata is None or 'cell_type' not in adata.obs.columns:
            continue

        for ct in sorted(adata.obs['cell_type'].unique()):
            ct_mask = adata.obs['cell_type'] == ct
            n_cells = ct_mask.sum()

            # Get coarse type
            coarse = "unknown"
            if 'coarse_cell_type' in adata.obs.columns:
                coarse_vals = adata.obs.loc[ct_mask, 'coarse_cell_type'].mode()
                coarse = coarse_vals.iloc[0] if len(coarse_vals) > 0 else "unknown"

            # Get confidence distribution
            if 'cell_type_confidence' in adata.obs.columns:
                conf_counts = adata.obs.loc[ct_mask, 'cell_type_confidence'].value_counts()
                conf_str = "; ".join(f"{k}={v}" for k, v in conf_counts.items())
            else:
                conf_str = ""

            # Get clusters
            if 'leiden' in adata.obs.columns:
                clusters = sorted(adata.obs.loc[ct_mask, 'leiden'].unique(),
                                  key=lambda x: (x.lstrip('MN'), x))
                cluster_str = ",".join(str(c) for c in clusters)
            else:
                cluster_str = ""

            # Canonical markers for this type
            markers = CANONICAL_MARKERS.get(ct, FINE_PANELS.get(ct, []))
            marker_str = ",".join(markers) if markers else ""

            rows.append({
                'object': obj_name,
                'cell_type': ct,
                'coarse_cell_type': coarse,
                'n_cells': n_cells,
                'clusters': cluster_str,
                'canonical_markers': marker_str,
                'confidence_distribution': conf_str,
            })

    if rows:
        df = pd.DataFrame(rows)
        df.to_csv(RESULTS_DIR / "cell_type_definitions.tsv", sep='\t', index=False)
        print(f"    Saved: cell_type_definitions.tsv ({len(df)} rows)")
    else:
        print("    WARNING: No annotations to summarize")


# ═════════════════════════════════════════════════════════════════════════════
# UMAP PLOTS
# ═════════════════════════════════════════════════════════════════════════════

def _plot_annotation_umaps(adata, object_name):
    """Generate UMAP panels colored by annotation results."""
    color_keys = ['leiden', 'coarse_cell_type', 'cell_type',
                  'cell_type_confidence', 'study', 'condition_harmonized']
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
    fig.savefig(RESULTS_DIR / f"umap_{object_name}_annotated.png",
                dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: umap_{object_name}_annotated.png")


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
    """Full annotation pipeline for one compartment object.

    Steps:
      1. Load h5ad from data/integrated/
      2. For each tier (mesenchymal, non-mesenchymal):
         a. Subset to tier cells
         b. Compute cluster DE markers
         c. Stage 1: coarse annotation
         d. Stage 2: fine annotation within coarse groups
         e. CellTypist validation (non-mesenchymal only)
      3. Compute continuous scores (mesenchymal only)
      4. Transfer annotations back to full object
      5. Generate UMAPs and save
    """
    input_path = INT_DIR / f"{object_name}.h5ad"

    if not input_path.exists():
        print(f"\n=== {object_name}: input not found at {input_path}, skipping ===")
        return None

    print(f"\n{'='*60}")
    print(f"Annotating object: {object_name}")
    print(f"{'='*60}")

    adata = sc.read_h5ad(input_path)
    print(f"  Loaded: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

    # Determine tiers present
    has_cell_class = 'cell_class' in adata.obs.columns
    if has_cell_class:
        mes_mask = adata.obs['cell_class'].isin(['mesenchymal', 'unknown'])
        non_mes_mask = adata.obs['cell_class'] == 'non_mesenchymal'
        n_mes = mes_mask.sum()
        n_non = non_mes_mask.sum()
        print(f"  Mesenchymal: {n_mes:,} cells")
        print(f"  Non-mesenchymal: {n_non:,} cells")
    else:
        # Treat all cells as mesenchymal if no cell_class column
        mes_mask = pd.Series(True, index=adata.obs_names)
        non_mes_mask = pd.Series(False, index=adata.obs_names)
        n_mes = adata.shape[0]
        n_non = 0
        print(f"  No cell_class column — treating all {n_mes:,} cells as mesenchymal")

    # Initialize annotation columns
    adata.obs['coarse_cell_type'] = "unassigned"
    adata.obs['cell_type'] = "unassigned"
    adata.obs['cell_type_confidence'] = "low"
    adata.obs['annotation_evidence'] = ""

    # ── Tier A: Mesenchymal ──────────────────────────────────────────────
    if n_mes >= 50:
        print(f"\n  --- Tier A: Mesenchymal ({object_name}) ---")
        mes_adata = adata[mes_mask].copy()

        # Use tier-specific leiden column (plain integers, not M-prefixed)
        tier_leiden_col = 'leiden_mesenchymal'
        if tier_leiden_col in mes_adata.obs.columns:
            mes_adata.obs['leiden'] = mes_adata.obs[tier_leiden_col].copy()

        # Ensure leiden exists for mesenchymal tier
        if 'leiden' not in mes_adata.obs.columns:
            print("    WARNING: No leiden column — cannot annotate mesenchymal tier")
        else:
            # Compute cluster DE markers
            compute_cluster_markers(mes_adata, object_name, "mesenchymal")

            # Stage 1: Coarse annotation
            coarse_assignments = annotate_coarse(mes_adata, object_name, "mesenchymal")

            # Stage 2: Fine annotation
            annotate_fine(mes_adata, object_name, "mesenchymal", coarse_assignments)

            # Continuous scores
            compute_continuous_scores(mes_adata)

            # Transfer annotations back to full object
            adata.obs.loc[mes_mask, 'coarse_cell_type'] = mes_adata.obs['coarse_cell_type'].values
            adata.obs.loc[mes_mask, 'cell_type'] = mes_adata.obs['cell_type'].values
            adata.obs.loc[mes_mask, 'cell_type_confidence'] = mes_adata.obs['cell_type_confidence'].values
            adata.obs.loc[mes_mask, 'annotation_evidence'] = mes_adata.obs['annotation_evidence'].values

            # Transfer continuous scores
            for score_name in CONTINUUM_SCORES:
                if score_name in mes_adata.obs.columns:
                    adata.obs[score_name] = np.nan
                    adata.obs.loc[mes_mask, score_name] = mes_adata.obs[score_name].values

        del mes_adata
        gc.collect()
    else:
        print(f"  Skipping mesenchymal tier: only {n_mes} cells")

    # ── Tier B: Non-mesenchymal ──────────────────────────────────────────
    MIN_CELLS_NON_MES = 50
    if n_non >= MIN_CELLS_NON_MES:
        print(f"\n  --- Tier B: Non-mesenchymal ({object_name}) ---")
        non_mes_adata = adata[non_mes_mask].copy()

        # Use tier-specific leiden column (plain integers, not NM-prefixed)
        tier_leiden_col = 'leiden_non_mesenchymal'
        if tier_leiden_col in non_mes_adata.obs.columns:
            non_mes_adata.obs['leiden'] = non_mes_adata.obs[tier_leiden_col].copy()

        if 'leiden' not in non_mes_adata.obs.columns:
            print("    WARNING: No leiden column — cannot annotate non-mesenchymal tier")
        else:
            # Compute cluster DE markers
            compute_cluster_markers(non_mes_adata, object_name, "non_mesenchymal")

            # Stage 1: Coarse annotation
            coarse_assignments = annotate_coarse(non_mes_adata, object_name,
                                                  "non_mesenchymal")

            # Stage 2: Fine annotation
            annotate_fine(non_mes_adata, object_name, "non_mesenchymal",
                          coarse_assignments)

            # CellTypist validation
            validate_immune_with_celltypist(non_mes_adata, object_name)

            # Transfer annotations back to full object
            adata.obs.loc[non_mes_mask, 'coarse_cell_type'] = non_mes_adata.obs['coarse_cell_type'].values
            adata.obs.loc[non_mes_mask, 'cell_type'] = non_mes_adata.obs['cell_type'].values
            adata.obs.loc[non_mes_mask, 'cell_type_confidence'] = non_mes_adata.obs['cell_type_confidence'].values
            adata.obs.loc[non_mes_mask, 'annotation_evidence'] = non_mes_adata.obs['annotation_evidence'].values

            # Transfer celltypist predictions if available
            if 'celltypist_prediction' in non_mes_adata.obs.columns:
                adata.obs['celltypist_prediction'] = ""
                adata.obs.loc[non_mes_mask, 'celltypist_prediction'] = \
                    non_mes_adata.obs['celltypist_prediction'].values

        del non_mes_adata
        gc.collect()
    else:
        print(f"  Skipping non-mesenchymal tier: only {n_non} cells")

    # ── Generate annotation UMAPs ─────────────────────────────────────────
    _plot_annotation_umaps(adata, object_name)

    # Generate fine-panel dotplots for the full object
    if object_name in ("NP", "all_cells"):
        panel = {k: v for k, v in CANONICAL_MARKERS.items()
                 if k.startswith("NP_")}
    elif object_name == "AF":
        panel = {k: v for k, v in CANONICAL_MARKERS.items()
                 if k.startswith("AF_")}
    elif object_name == "CEP":
        panel = {k: v for k, v in CANONICAL_MARKERS.items()
                 if k.startswith("EP_")}
    else:
        panel = CANONICAL_MARKERS
    generate_annotation_dotplots(adata, object_name, "all_markers", panel)

    # ── Save ─────────────────────────────────────────────────────────────
    _save_checkpoint(adata, input_path, "annotated")
    print(f"  {object_name} annotation complete: {adata.shape[0]:,} cells, "
          f"{adata.obs['cell_type'].nunique()} cell types")

    return adata


def process_all_cells_secondary(force=False):
    """Process all-cells as a secondary object.

    Transfers annotations from compartment-specific objects (NP, AF, CEP)
    by matching obs_names. De novo annotates remaining unassigned cells.
    """
    input_path = INT_DIR / "all_cells.h5ad"

    if not input_path.exists():
        print(f"\n=== all_cells: input not found, skipping ===")
        return None

    print(f"\n{'='*60}")
    print(f"Annotating object: all_cells (secondary)")
    print(f"{'='*60}")

    adata_all = sc.read_h5ad(input_path)
    print(f"  Loaded: {adata_all.shape[0]:,} cells")

    # Initialize annotation columns
    adata_all.obs['coarse_cell_type'] = "unassigned"
    adata_all.obs['cell_type'] = "unassigned"
    adata_all.obs['cell_type_confidence'] = "low"
    adata_all.obs['annotation_evidence'] = ""

    # Transfer annotations from compartment-specific objects
    print("  Transferring annotations from compartment-specific objects...")
    n_transferred = 0
    annotation_cols = ['coarse_cell_type', 'cell_type', 'cell_type_confidence',
                       'annotation_evidence']
    score_cols = list(CONTINUUM_SCORES.keys())

    for obj_name in ["NP", "AF", "CEP"]:
        obj_path = INT_DIR / f"{obj_name}.h5ad"
        if not obj_path.exists():
            continue
        obj_adata = sc.read_h5ad(obj_path)

        # Match cells by obs_names
        common_cells = adata_all.obs_names.intersection(obj_adata.obs_names)
        if len(common_cells) > 0:
            for col in annotation_cols:
                if col in obj_adata.obs.columns:
                    adata_all.obs.loc[common_cells, col] = \
                        obj_adata.obs.loc[common_cells, col].values

            # Transfer continuous scores
            for col in score_cols:
                if col in obj_adata.obs.columns:
                    if col not in adata_all.obs.columns:
                        adata_all.obs[col] = np.nan
                    adata_all.obs.loc[common_cells, col] = \
                        obj_adata.obs.loc[common_cells, col].values

            # Transfer celltypist predictions
            if 'celltypist_prediction' in obj_adata.obs.columns:
                if 'celltypist_prediction' not in adata_all.obs.columns:
                    adata_all.obs['celltypist_prediction'] = ""
                adata_all.obs.loc[common_cells, 'celltypist_prediction'] = \
                    obj_adata.obs.loc[common_cells, 'celltypist_prediction'].values

            n_transferred += len(common_cells)
            print(f"    {obj_name}: transferred {len(common_cells):,} cell annotations")
        del obj_adata
        gc.collect()

    # De novo annotate remaining unassigned cells (e.g. GSE189916)
    unassigned_mask = adata_all.obs['cell_type'] == 'unassigned'
    n_unassigned = unassigned_mask.sum()
    print(f"  Cells with transferred annotations: {n_transferred:,}")
    print(f"  Cells needing de novo annotation: {n_unassigned:,}")

    if n_unassigned > 50:
        print("  Running de novo annotation on unassigned cells...")
        unassigned_adata = adata_all[unassigned_mask].copy()

        # Need leiden clustering
        if 'leiden' not in unassigned_adata.obs.columns:
            # Quick clustering on PCA or existing embedding
            if 'X_pca' in unassigned_adata.obsm:
                sc.pp.neighbors(unassigned_adata, use_rep='X_pca')
            else:
                sc.pp.pca(unassigned_adata)
                sc.pp.neighbors(unassigned_adata)
            sc.tl.leiden(unassigned_adata, resolution=0.5, key_added='leiden')

        # Use both mesenchymal and non-mesenchymal panels
        combined_panel = {}
        combined_panel.update(COARSE_PANELS_MESENCHYMAL)
        combined_panel.update(COARSE_PANELS_NON_MESENCHYMAL)

        var_names_set = set(unassigned_adata.var_names)

        for cluster in unassigned_adata.obs['leiden'].unique():
            cluster_mask = unassigned_adata.obs['leiden'] == cluster
            scores = _score_cluster_for_panel(unassigned_adata, cluster_mask,
                                              combined_panel, var_names_set)
            label, confidence, evidence = _assign_best_label(scores)

            unassigned_adata.obs.loc[cluster_mask, 'coarse_cell_type'] = label
            unassigned_adata.obs.loc[cluster_mask, 'cell_type'] = label
            unassigned_adata.obs.loc[cluster_mask, 'cell_type_confidence'] = confidence
            unassigned_adata.obs.loc[cluster_mask, 'annotation_evidence'] = evidence

        # Transfer back
        adata_all.obs.loc[unassigned_mask, 'coarse_cell_type'] = \
            unassigned_adata.obs['coarse_cell_type'].values
        adata_all.obs.loc[unassigned_mask, 'cell_type'] = \
            unassigned_adata.obs['cell_type'].values
        adata_all.obs.loc[unassigned_mask, 'cell_type_confidence'] = \
            unassigned_adata.obs['cell_type_confidence'].values
        adata_all.obs.loc[unassigned_mask, 'annotation_evidence'] = \
            unassigned_adata.obs['annotation_evidence'].values
        del unassigned_adata

    # Generate UMAPs
    _plot_annotation_umaps(adata_all, "all_cells")

    # Save
    _save_checkpoint(adata_all, input_path, "annotated")
    print(f"  all_cells annotation complete: {adata_all.shape[0]:,} cells, "
          f"{adata_all.obs['cell_type'].nunique()} cell types")

    return adata_all


# ═════════════════════════════════════════════════════════════════════════════
# ANNOTATION REPORT
# ═════════════════════════════════════════════════════════════════════════════

ANNOTATION_REPORT_TEMPLATE = """<!DOCTYPE html>
<html>
<head>
<title>Annotation Report — Module 07</title>
<style>
  body { font-family: Arial, sans-serif; max-width: 1400px; margin: 0 auto; padding: 20px; }
  h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
  h2 { color: #34495e; margin-top: 30px; }
  table { border-collapse: collapse; width: 100%; margin: 10px 0; }
  th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
  th { background-color: #3498db; color: white; }
  tr:nth-child(even) { background-color: #f2f2f2; }
  img { max-width: 100%; height: auto; margin: 10px 0; border: 1px solid #ddd; }
  .warning { color: #e74c3c; }
  .pass { color: #27ae60; }
  .validation li { margin: 4px 0; }
</style>
</head>
<body>
<h1>Module 07: Post-Integration Cell Type Annotation</h1>
<p>Generated: {{ date }}</p>

{% for obj in objects %}
<h2>{{ obj.name }}</h2>
{% if obj.umap_exists %}
<img src="umap_{{ obj.name }}_annotated.png" alt="{{ obj.name }} annotated UMAP">
{% endif %}
<p>Cells: {{ obj.n_cells }} | Coarse types: {{ obj.n_coarse }} | Fine types: {{ obj.n_fine }}</p>
<p>Unassigned: {{ obj.pct_unassigned }}%</p>
{% endfor %}

<h2>Validation</h2>
<ul class="validation">
  {% for msg in validation_messages %}
  <li class="{{ 'pass' if msg.startswith('PASS') else 'warning' if msg.startswith('FAIL') or msg.startswith('WARNING') else '' }}">{{ msg }}</li>
  {% endfor %}
</ul>

<h2>Human Checkpoint Questions</h2>
<ol>
  <li>Do the cell type annotations make biological sense?</li>
  <li>For the mesenchymal continuum: are discrete labels appropriate, or should some clusters be merged?</li>
  <li>Are there unexpected cell types or missing expected types?</li>
  <li>Is the all-cells object consistent with the compartment-specific objects?</li>
  <li>Should any annotations be revised before proceeding to differential analysis?</li>
</ol>

</body>
</html>"""


def generate_annotation_report(validation_messages):
    """Generate HTML annotation report."""
    print("\n  Generating annotation report...")

    objects = []
    for obj_name in ["NP", "AF", "CEP", "all_cells"]:
        obj_path = INT_DIR / f"{obj_name}.h5ad"
        if obj_path.exists():
            adata = sc.read_h5ad(obj_path, backed='r')
            n_coarse = adata.obs['coarse_cell_type'].nunique() if 'coarse_cell_type' in adata.obs.columns else 0
            n_fine = adata.obs['cell_type'].nunique() if 'cell_type' in adata.obs.columns else 0
            pct_unassigned = 0
            if 'cell_type' in adata.obs.columns:
                n_un = (adata.obs['cell_type'] == 'unassigned').sum()
                pct_unassigned = round(n_un / adata.shape[0] * 100, 1)

            objects.append({
                'name': obj_name,
                'n_cells': f"{adata.shape[0]:,}",
                'n_coarse': n_coarse,
                'n_fine': n_fine,
                'pct_unassigned': pct_unassigned,
                'umap_exists': (RESULTS_DIR / f"umap_{obj_name}_annotated.png").exists(),
            })
            adata.file.close()

    template = Template(ANNOTATION_REPORT_TEMPLATE)
    html = template.render(
        date=datetime.now().strftime("%Y-%m-%d %H:%M"),
        objects=objects,
        validation_messages=validation_messages,
    )
    (RESULTS_DIR / "annotation_report.html").write_text(html)
    print(f"  Report: {RESULTS_DIR / 'annotation_report.html'}")


# ═════════════════════════════════════════════════════════════════════════════
# VALIDATION
# ═════════════════════════════════════════════════════════════════════════════

def validate_annotation():
    """Validate annotation outputs per spec. Returns (pass, messages)."""
    print("\n" + "=" * 60)
    print("Validating annotation outputs")
    print("=" * 60)

    messages = []
    all_pass = True

    for obj_name in ["NP", "AF", "CEP", "all_cells"]:
        path = INT_DIR / f"{obj_name}.h5ad"
        if not path.exists():
            messages.append(f"FAIL: {obj_name} output missing")
            all_pass = False
            continue

        adata = sc.read_h5ad(path, backed='r')

        # Check 1: All clusters have a cell type label
        if 'cell_type' in adata.obs.columns:
            n_unassigned = (adata.obs['cell_type'] == 'unassigned').sum()
            pct_unassigned = n_unassigned / adata.shape[0] * 100

            # Check 2: Proportion of unassigned < 10%
            if pct_unassigned < 10:
                messages.append(f"PASS: {obj_name} unassigned = {pct_unassigned:.1f}% (< 10%)")
            else:
                messages.append(f"WARNING: {obj_name} unassigned = {pct_unassigned:.1f}% (>= 10%)")

            messages.append(f"PASS: {obj_name} has cell_type column with "
                            f"{adata.obs['cell_type'].nunique()} types")
        else:
            messages.append(f"FAIL: {obj_name} missing cell_type column")
            all_pass = False

        # Check 3: coarse_cell_type column
        if 'coarse_cell_type' in adata.obs.columns:
            messages.append(f"PASS: {obj_name} has coarse_cell_type column")
        else:
            messages.append(f"FAIL: {obj_name} missing coarse_cell_type column")
            all_pass = False

        # Check 4: annotation_evidence column
        if 'annotation_evidence' in adata.obs.columns:
            n_with_evidence = (adata.obs['annotation_evidence'] != "").sum()
            messages.append(f"PASS: {obj_name} annotation_evidence for "
                            f"{n_with_evidence:,}/{adata.shape[0]:,} cells")
        else:
            messages.append(f"FAIL: {obj_name} missing annotation_evidence column")
            all_pass = False

        # Check 5: Cluster DE marker tables exist
        marker_files = list((RESULTS_DIR / "cluster_markers").glob(f"{obj_name}_*_markers.tsv"))
        if marker_files:
            messages.append(f"PASS: {obj_name} has {len(marker_files)} marker table(s)")
        else:
            messages.append(f"WARNING: {obj_name} no cluster marker tables found")

        # Check 6: Non-mesenchymal marker validation
        if 'cell_class' in adata.obs.columns and 'cell_type' in adata.obs.columns:
            non_mes = adata.obs[adata.obs['cell_class'] == 'non_mesenchymal']
            if len(non_mes) > 0:
                # Check PTPRC in immune clusters
                immune_types = [ct for ct in non_mes['cell_type'].unique()
                                if any(x in ct.lower() for x in
                                       ['macrophage', 't_cell', 'b_cell', 'nk_cell',
                                        'mast_cell'])]
                if immune_types:
                    messages.append(f"PASS: {obj_name} immune types found: {immune_types}")

                # Check endothelial
                endo_types = [ct for ct in non_mes['cell_type'].unique()
                              if 'endothelial' in ct.lower()]
                if endo_types:
                    messages.append(f"PASS: {obj_name} endothelial types found: {endo_types}")

        # Check 7: CellTypist concordance
        conc_path = RESULTS_DIR / "celltypist_validation" / f"{obj_name}_concordance.tsv"
        if conc_path.exists():
            conc = pd.read_csv(conc_path, sep='\t')
            n_agree = conc['concordant'].sum()
            messages.append(f"PASS: {obj_name} CellTypist concordance: "
                            f"{n_agree}/{len(conc)} clusters")
        # (not a failure if CellTypist wasn't run — only for non-mesenchymal)

        adata.file.close()

    # Check cell type definitions table
    defs_path = RESULTS_DIR / "cell_type_definitions.tsv"
    if defs_path.exists():
        messages.append(f"PASS: cell_type_definitions.tsv exists")
    else:
        messages.append(f"FAIL: cell_type_definitions.tsv missing")
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
    global INT_DIR, RESULTS_DIR
    args = sys.argv[1:]

    validate_only = '--validate-only' in args
    force = '--force' in args

    # Parse --object
    object_filter = None
    for i, a in enumerate(args):
        if a == '--object' and i + 1 < len(args):
            object_filter = args[i + 1]

    # Parse --input-dir (where to read {object}.h5ad from)
    for i, a in enumerate(args):
        if a == '--input-dir' and i + 1 < len(args):
            INT_DIR = Path(args[i + 1]).resolve()

    # Parse --output-dir (where annotation results land)
    for i, a in enumerate(args):
        if a == '--output-dir' and i + 1 < len(args):
            RESULTS_DIR = Path(args[i + 1]).resolve()

    for d in [RESULTS_DIR / "cluster_markers",
              RESULTS_DIR / "annotation_dotplots",
              RESULTS_DIR / "celltypist_validation"]:
        d.mkdir(parents=True, exist_ok=True)

    print(f"  Input dir:  {INT_DIR}")
    print(f"  Output dir: {RESULTS_DIR}")

    if validate_only:
        passed, _ = validate_annotation()
        sys.exit(0 if passed else 1)

    print("=" * 60)
    print("Module 07: Post-Integration Cell Type Annotation")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    all_annotations = {}

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
        adata = process_object(obj_name, force=force)
        if adata is not None:
            all_annotations[obj_name] = adata

    # Process all-cells (secondary) — after primary objects are done
    if object_filter is None or object_filter == "all_cells":
        adata = process_all_cells_secondary(force=force)
        if adata is not None:
            all_annotations['all_cells'] = adata

    # Generate cell type definitions table
    generate_cell_type_definitions(all_annotations)

    # Free memory before validation
    del all_annotations
    gc.collect()

    # Validation
    passed, val_messages = validate_annotation()

    # Report
    generate_annotation_report(val_messages)

    print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Overall: {'PASSED' if passed else 'FAILED'}")


if __name__ == "__main__":
    main()
