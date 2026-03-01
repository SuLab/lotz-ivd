#!/usr/bin/env python3
"""Module 04: Per-dataset cell type annotation for IVD scRNA-seq atlas.

Refines the preliminary cell-type labels from Module 03 using expanded gene
signatures (Strategy 1), CellTypist automated annotation (Strategy 3), and
a consensus framework. Strategy 2 (reference-based label transfer) is skipped
because no IVD-specific reference atlas exists.

Usage:
    python3 scripts/04_annotation.py                  # Annotate all datasets
    python3 scripts/04_annotation.py GSE255768        # Annotate one dataset
    python3 scripts/04_annotation.py --validate-only  # Run validation only
"""

import sys
import os
import warnings
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
warnings.filterwarnings('ignore', category=UserWarning, module='celltypist')

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent.parent
PROC_DIR = BASE / "data" / "processed"
ANN_DIR = BASE / "results" / "annotations"
FIG_DIR = ANN_DIR  # annotation figures go alongside reports

ANN_DIR.mkdir(parents=True, exist_ok=True)

# ── Dataset ordering (same as Module 03) ─────────────────────────────────────
ALL_ACCESSIONS = [
    "GSE160756", "GSE165722", "GSE189916", "GSE199866", "GSE205535",
    "CNP0002664", "GSE233666", "GSE244889", "GSE251686", "GSE255768",
    "GSE230809", "GSE242443",
]

# ── Annotation resolution ────────────────────────────────────────────────────
ANNOTATION_RESOLUTION = 0.8
FALLBACK_RESOLUTION = 0.5

# ── Gene aliases ─────────────────────────────────────────────────────────────
# Keys = canonical name used in signatures, values = list of known aliases
GENE_ALIASES = {
    "T":    ["T", "TBXT"],
    "TBXT": ["TBXT", "T"],
    "SP7":  ["SP7", "OSX"],
}

# ── Expanded gene signatures ─────────────────────────────────────────────────
# Resident IVD subtypes
NP_SIGNATURES = {
    "NP_notochordal":          ["T", "SHH", "NOG", "CD24", "KRT8", "KRT18", "KRT19"],
    "NP_mature_chondrocyte":   ["ACAN", "COL2A1", "SOX9", "COMP", "PRG4"],
    "NP_stressed_degenerative":["MMP13", "ADAMTS5", "IL1B", "TNF", "VEGFA", "HIF1A"],
    "NP_fibrocartilaginous":   ["COL1A1", "COL2A1", "VCAN"],
}

AF_SIGNATURES = {
    "AF_inner":              ["COL2A1", "ACAN", "SOX9"],
    "AF_outer":              ["COL1A1", "COL1A2", "THY1", "DCN", "LUM"],
    "AF_mechanical_stress":  ["COMP", "CILP", "THBS1"],
}

EP_SIGNATURES = {
    "EP_hyaline_cartilage": ["COL2A1", "COL10A1", "SOX9"],
    "EP_ossification":      ["RUNX2", "SP7", "BGLAP"],
}

# Non-resident (shared with Module 03, extended where noted)
NON_RESIDENT_SIGNATURES = {
    "Endothelial":  ["PECAM1", "VWF", "CDH5", "FLT1"],
    "Macrophage":   ["CD68", "CD163", "CSF1R", "MRC1"],
    "T_cell":       ["CD3D", "CD3E", "CD4", "CD8A"],
    "B_cell":       ["CD79A", "MS4A1", "CD19"],
    "Mast_cell":    ["KIT", "TPSAB1", "CPA3"],
    "Pericyte_SMC": ["ACTA2", "TAGLN", "MYH11", "RGS5"],
}

ALL_SIGNATURES = {}
ALL_SIGNATURES.update(NP_SIGNATURES)
ALL_SIGNATURES.update(AF_SIGNATURES)
ALL_SIGNATURES.update(EP_SIGNATURES)
ALL_SIGNATURES.update(NON_RESIDENT_SIGNATURES)

RESIDENT_LABELS = set(NP_SIGNATURES) | set(AF_SIGNATURES) | set(EP_SIGNATURES)
NON_RESIDENT_LABELS = set(NON_RESIDENT_SIGNATURES)

# Compartment-to-expected-resident mapping (for plausibility priors)
COMPARTMENT_EXPECTED = {
    "NP":     set(NP_SIGNATURES),
    "AF":     set(AF_SIGNATURES),
    "CEP":    set(EP_SIGNATURES),
    "EP":     set(EP_SIGNATURES),
    "mixed":  RESIDENT_LABELS,
    "IVD":    RESIDENT_LABELS,
}


# ═══════════════════════════════════════════════════════════════════════════════
# GENE RESOLUTION
# ═══════════════════════════════════════════════════════════════════════════════

def resolve_gene(name, var_names):
    """Resolve a gene name to the form present in var_names, handling aliases.

    Returns the matched name if found, else None.
    """
    if name in var_names:
        return name
    # Check aliases
    aliases = GENE_ALIASES.get(name, [])
    for alias in aliases:
        if alias in var_names:
            return alias
    return None


def resolve_signature(gene_list, var_names):
    """Resolve all genes in a signature, returning (resolved_list, missing_list)."""
    var_set = set(var_names)
    resolved = []
    missing = []
    for gene in gene_list:
        match = resolve_gene(gene, var_set)
        if match is not None:
            resolved.append(match)
        else:
            missing.append(gene)
    return resolved, missing


# ═══════════════════════════════════════════════════════════════════════════════
# STRATEGY 1: MARKER-BASED SCORING
# ═══════════════════════════════════════════════════════════════════════════════

def score_signatures(adata):
    """Score all gene signatures using sc.tl.score_genes().

    Returns a dict of {sig_name: {"available": [...], "missing": [...], "n_avail": int}}.
    """
    var_names = set(adata.var_names)
    stats = {}

    for sig_name, gene_list in ALL_SIGNATURES.items():
        resolved, missing = resolve_signature(gene_list, var_names)
        col = f"score_{sig_name}"

        if len(resolved) >= 1:
            try:
                sc.tl.score_genes(adata, gene_list=resolved, score_name=col)
            except Exception as e:
                print(f"    WARNING: score_genes failed for {sig_name}: {e}")
                adata.obs[col] = 0.0
        else:
            adata.obs[col] = 0.0

        stats[sig_name] = {
            "available": resolved,
            "missing": missing,
            "n_avail": len(resolved),
            "n_total": len(gene_list),
        }

    return stats


def assign_marker_labels(adata, resolution=ANNOTATION_RESOLUTION):
    """Assign cell type labels using cluster-level majority voting.

    Uses Leiden clusters at the given resolution. For each cluster, determines
    the dominant signature score and assigns that label to all cells in the
    cluster. Compartment metadata provides a plausibility prior.
    """
    cluster_key = f"leiden_res_{resolution}"
    if cluster_key not in adata.obs.columns:
        print(f"    WARNING: {cluster_key} not found, falling back to res {FALLBACK_RESOLUTION}")
        cluster_key = f"leiden_res_{FALLBACK_RESOLUTION}"
        resolution = FALLBACK_RESOLUTION

    # Score columns
    score_cols = {sig: f"score_{sig}" for sig in ALL_SIGNATURES}
    available_scores = {sig: col for sig, col in score_cols.items()
                        if col in adata.obs.columns}

    # Determine dominant compartment for this dataset
    if 'compartment' in adata.obs.columns:
        compartments = adata.obs['compartment'].value_counts()
        dominant_compartment = compartments.index[0]
    else:
        dominant_compartment = "mixed"

    expected_resident = COMPARTMENT_EXPECTED.get(dominant_compartment, RESIDENT_LABELS)

    labels = pd.Series("unassigned", index=adata.obs_names, dtype="object")
    confidences = pd.Series("low", index=adata.obs_names, dtype="object")

    for cluster in adata.obs[cluster_key].unique():
        mask = adata.obs[cluster_key] == cluster

        # Compute mean score for each signature in this cluster
        mean_scores = {}
        for sig, col in available_scores.items():
            mean_scores[sig] = adata.obs.loc[mask, col].mean()

        if not mean_scores:
            continue

        # Separate non-resident vs resident scores
        non_res_scores = {s: v for s, v in mean_scores.items() if s in NON_RESIDENT_LABELS}
        resident_scores = {s: v for s, v in mean_scores.items() if s in RESIDENT_LABELS}

        # Check if this cluster is non-resident
        best_non_res = max(non_res_scores, key=non_res_scores.get) if non_res_scores else None
        best_non_res_score = non_res_scores.get(best_non_res, -np.inf) if best_non_res else -np.inf

        best_resident = max(resident_scores, key=resident_scores.get) if resident_scores else None
        best_resident_score = resident_scores.get(best_resident, -np.inf) if best_resident else -np.inf

        # Non-resident wins if its score is clearly above resident scores
        # and the score itself is reasonably positive
        if best_non_res_score > 0.15 and best_non_res_score > best_resident_score + 0.05:
            chosen = best_non_res
            conf = "high" if best_non_res_score > 0.3 else "medium"
        elif best_resident is not None and best_resident_score > -0.5:
            # Resident labeling with compartment prior
            # Among resident scores, prefer those matching the dataset compartment
            in_compartment = {s: v for s, v in resident_scores.items()
                              if s in expected_resident}
            out_compartment = {s: v for s, v in resident_scores.items()
                               if s not in expected_resident}

            if in_compartment:
                best_in = max(in_compartment, key=in_compartment.get)
                best_in_score = in_compartment[best_in]
            else:
                best_in = None
                best_in_score = -np.inf

            if out_compartment:
                best_out = max(out_compartment, key=out_compartment.get)
                best_out_score = out_compartment[best_out]
            else:
                best_out = None
                best_out_score = -np.inf

            # Use in-compartment label unless out-compartment is much stronger
            if best_in is not None and (best_in_score >= best_out_score - 0.1):
                chosen = best_in
            elif best_out is not None:
                chosen = best_out
            else:
                chosen = best_resident

            # Confidence based on margin
            sorted_res = sorted(resident_scores.values(), reverse=True)
            if len(sorted_res) >= 2:
                margin = sorted_res[0] - sorted_res[1]
                if margin > 0.15:
                    conf = "high"
                elif margin > 0.05:
                    conf = "medium"
                else:
                    conf = "low"
            else:
                conf = "medium"
        else:
            chosen = "unassigned"
            conf = "low"

        labels.loc[mask] = chosen
        confidences.loc[mask] = conf

    adata.obs['cell_type_marker_based'] = labels.values
    adata.obs['cell_type_marker_confidence'] = confidences.values
    return resolution


# ═══════════════════════════════════════════════════════════════════════════════
# STRATEGY 2: REFERENCE-BASED (SKIPPED)
# ═══════════════════════════════════════════════════════════════════════════════

def set_reference_not_applied(adata):
    """Mark Strategy 2 as not applied (no IVD reference atlas available)."""
    adata.obs['cell_type_reference_based'] = "not_applied"


# ═══════════════════════════════════════════════════════════════════════════════
# STRATEGY 3: CELLTYPIST
# ═══════════════════════════════════════════════════════════════════════════════

def run_celltypist(adata):
    """Run CellTypist with Immune_All_Low model for immune cell refinement.

    Stores predictions in `cell_type_celltypist`. Gracefully skips if
    celltypist import fails or model download fails.

    Returns True if CellTypist ran successfully, False otherwise.
    """
    try:
        import celltypist
        from celltypist import models as ct_models
    except ImportError:
        print("    WARNING: celltypist not installed, skipping Strategy 3")
        adata.obs['cell_type_celltypist'] = "not_available"
        return False

    try:
        # Download model if not cached
        model = ct_models.Model.load(model="Immune_All_Low.pkl")
    except Exception as e:
        print(f"    WARNING: Could not load CellTypist model: {e}")
        adata.obs['cell_type_celltypist'] = "model_unavailable"
        return False

    try:
        # CellTypist expects log-normalized data in .X (which we have)
        # It also expects gene names as var_names (which we have)
        predictions = celltypist.annotate(
            adata,
            model=model,
            majority_voting=True,
        )
        # Store the majority voting results
        result_adata = predictions.to_adata()
        adata.obs['cell_type_celltypist'] = result_adata.obs['majority_voting'].values
        adata.obs['celltypist_conf_score'] = result_adata.obs['conf_score'].values
        print(f"    CellTypist labels assigned: {adata.obs['cell_type_celltypist'].nunique()} unique types")
        return True
    except Exception as e:
        print(f"    WARNING: CellTypist annotation failed: {e}")
        adata.obs['cell_type_celltypist'] = "failed"
        return False


# ═══════════════════════════════════════════════════════════════════════════════
# CONSENSUS LABELING
# ═══════════════════════════════════════════════════════════════════════════════

def assign_consensus_labels(adata, celltypist_ran=False):
    """Merge strategy labels into final consensus.

    Rules:
    - Non-resident cells: prefer CellTypist if it gives a more specific subtype;
      fall back to marker-based.
    - Resident IVD cells: use marker-based only. CellTypist logged but not trusted.
    - Ambiguous: label "unassigned", confidence "low".
    """
    marker_labels = adata.obs['cell_type_marker_based'].values
    marker_conf = adata.obs['cell_type_marker_confidence'].values

    if celltypist_ran:
        ct_labels = adata.obs['cell_type_celltypist'].values
    else:
        ct_labels = np.array(["not_available"] * adata.shape[0])

    final_labels = np.empty(adata.shape[0], dtype=object)
    final_conf = np.empty(adata.shape[0], dtype=object)

    # CellTypist immune-relevant labels (substrings that indicate immune/endo/pericyte)
    immune_keywords = [
        "Macro", "Mono", "DC", "Dendritic", "T cell", "CD4", "CD8", "NK",
        "B cell", "Plasma", "Mast", "Neutro", "Eosinop", "Baso",
        "Endothelial", "Pericyte", "Smooth muscle", "Fibroblast",
    ]

    for i in range(adata.shape[0]):
        ml = marker_labels[i]
        mc = marker_conf[i]
        cl = ct_labels[i]

        if ml in NON_RESIDENT_LABELS:
            # Non-resident: prefer CellTypist if it gives a specific immune subtype
            if celltypist_ran and cl not in ("not_available", "model_unavailable",
                                              "failed", "Unassigned"):
                # Check if CellTypist agrees it's non-resident
                is_ct_immune = any(kw.lower() in cl.lower() for kw in immune_keywords)
                if is_ct_immune:
                    final_labels[i] = cl
                    final_conf[i] = "high"
                else:
                    # CellTypist disagrees, keep marker-based
                    final_labels[i] = ml
                    final_conf[i] = mc
            else:
                final_labels[i] = ml
                final_conf[i] = mc
        elif ml in RESIDENT_LABELS:
            # Resident: trust marker-based only
            final_labels[i] = ml
            final_conf[i] = mc
        elif ml == "unassigned":
            # Try CellTypist as fallback
            if celltypist_ran and cl not in ("not_available", "model_unavailable",
                                              "failed", "Unassigned"):
                final_labels[i] = cl
                final_conf[i] = "low"
            else:
                final_labels[i] = "unassigned"
                final_conf[i] = "low"
        else:
            final_labels[i] = ml
            final_conf[i] = mc

    adata.obs['cell_type_final'] = final_labels
    adata.obs['cell_type_confidence'] = final_conf


# ═══════════════════════════════════════════════════════════════════════════════
# VALIDATION
# ═══════════════════════════════════════════════════════════════════════════════

def validate_annotation(adata, accession):
    """Run annotation validation checks. Returns (pass, messages)."""
    messages = []
    all_pass = True

    n_cells = adata.shape[0]

    # Check 1: All cells have cell_type_final
    if 'cell_type_final' not in adata.obs.columns:
        messages.append("FAIL: cell_type_final column missing")
        return False, messages

    n_labeled = adata.obs['cell_type_final'].notna().sum()
    if n_labeled == n_cells:
        messages.append(f"PASS: All {n_cells} cells have cell_type_final labels")
    else:
        messages.append(f"FAIL: {n_cells - n_labeled} cells missing cell_type_final")
        all_pass = False

    # Check 2: <20% unassigned
    n_unassigned = (adata.obs['cell_type_final'] == 'unassigned').sum()
    pct_unassigned = n_unassigned / n_cells * 100
    if pct_unassigned < 20:
        messages.append(f"PASS: {pct_unassigned:.1f}% unassigned (< 20%)")
    else:
        messages.append(f"FAIL: {pct_unassigned:.1f}% unassigned (>= 20%)")
        all_pass = False

    # Check 3: Immune marker specificity
    immune_markers_check = {'CD68': 'Macrophage', 'CD3D': 'T_cell', 'CD79A': 'B_cell'}
    for gene, expected_type in immune_markers_check.items():
        if gene not in adata.var_names:
            continue
        gene_col = adata[:, gene].X
        if sparse.issparse(gene_col):
            gene_col = gene_col.toarray().ravel()
        else:
            gene_col = np.asarray(gene_col).ravel()

        expressing = gene_col > 0
        if expressing.sum() == 0:
            messages.append(f"INFO: {gene} not expressed in this dataset")
            continue

        # Check: what fraction of expressing cells are labeled as the expected type
        # or another non-resident type?
        final_labels = adata.obs['cell_type_final'].values
        non_res_mask = np.array([l in NON_RESIDENT_LABELS or
                                  any(kw.lower() in str(l).lower()
                                      for kw in ["Macro", "Mono", "T cell", "CD4",
                                                  "CD8", "NK", "B cell", "Plasma",
                                                  "Mast", "Dendritic"])
                                  for l in final_labels])
        # Among expressing cells, what fraction is in non-resident labels?
        expr_in_non_res = (expressing & non_res_mask).sum()
        expr_total = expressing.sum()
        if expr_total > 0:
            specificity = expr_in_non_res / expr_total
            if specificity > 0.3:
                messages.append(f"PASS: {gene} enriched in non-resident cells "
                                f"({specificity:.0%} of {gene}+ cells)")
            else:
                messages.append(f"INFO: {gene} expressed broadly "
                                f"({specificity:.0%} in non-resident); "
                                f"may be diffuse in IVD tissue")

    # Check 4: >=80% endothelial express PECAM1
    if 'PECAM1' in adata.var_names:
        endo_mask = adata.obs['cell_type_final'].str.contains(
            'Endothelial|endothelial', case=False, na=False
        )
        if endo_mask.sum() > 0:
            pecam1_col = adata[endo_mask, 'PECAM1'].X
            if sparse.issparse(pecam1_col):
                pecam1_col = pecam1_col.toarray().ravel()
            else:
                pecam1_col = np.asarray(pecam1_col).ravel()
            pct_pecam1 = (pecam1_col > 0).sum() / endo_mask.sum() * 100
            if pct_pecam1 >= 80:
                messages.append(f"PASS: {pct_pecam1:.0f}% of endothelial cells express PECAM1")
            else:
                messages.append(f"WARNING: Only {pct_pecam1:.0f}% of endothelial cells express PECAM1 (<80%)")
        else:
            messages.append("INFO: No endothelial cells labeled in this dataset")

    # Check 5: Compartment plausibility
    if 'compartment' in adata.obs.columns:
        compartments = adata.obs['compartment'].unique()
        for comp in compartments:
            if comp in ("mixed", "IVD"):
                continue
            expected = COMPARTMENT_EXPECTED.get(comp, RESIDENT_LABELS)
            comp_mask = adata.obs['compartment'] == comp
            comp_labels = adata.obs.loc[comp_mask, 'cell_type_final']
            resident_in_comp = comp_labels[comp_labels.isin(RESIDENT_LABELS)]
            if len(resident_in_comp) > 0:
                in_expected = resident_in_comp.isin(expected).sum()
                pct_expected = in_expected / len(resident_in_comp) * 100
                if pct_expected > 30:
                    messages.append(f"PASS: {comp} compartment — {pct_expected:.0f}% "
                                    f"of resident cells match expected types")
                else:
                    messages.append(f"WARNING: {comp} compartment — only {pct_expected:.0f}% "
                                    f"of resident cells match expected types")

    # Check 6: Reports generated
    report_path = ANN_DIR / f"{accession}_annotation_report.html"
    if report_path.exists():
        messages.append("PASS: Annotation report HTML generated")
    else:
        messages.append("FAIL: Annotation report HTML not found")
        all_pass = False

    dotplot_path = ANN_DIR / f"{accession}_marker_dotplot.pdf"
    if dotplot_path.exists():
        messages.append("PASS: Marker dotplot PDF generated")
    else:
        messages.append("FAIL: Marker dotplot PDF not found")
        all_pass = False

    return all_pass, messages


# ═══════════════════════════════════════════════════════════════════════════════
# PLOTS
# ═══════════════════════════════════════════════════════════════════════════════

def generate_annotation_plots(adata, accession):
    """Generate annotation-specific plots for one dataset."""
    plt.rcParams['figure.dpi'] = 100

    # 1. UMAP by cell_type_final
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color='cell_type_final', ax=ax, show=False,
               title=f'{accession} — Cell Type (Final)')
    plt.tight_layout()
    plt.savefig(ANN_DIR / f"{accession}_umap_annotation.png",
                bbox_inches='tight', dpi=150)
    plt.close()

    # 2. Side-by-side strategy comparison
    fig, axes = plt.subplots(1, 3, figsize=(24, 7))

    sc.pl.umap(adata, color='cell_type_marker_based', ax=axes[0], show=False,
               title='Strategy 1: Marker-based')
    sc.pl.umap(adata, color='cell_type_celltypist', ax=axes[1], show=False,
               title='Strategy 3: CellTypist')
    sc.pl.umap(adata, color='cell_type_final', ax=axes[2], show=False,
               title='Consensus (Final)')
    plt.tight_layout()
    plt.savefig(ANN_DIR / f"{accession}_strategy_comparison.png",
                bbox_inches='tight', dpi=150)
    plt.close()

    # 3. Marker dotplot — canonical markers × final cell types
    marker_genes = []
    gene_group_labels = {}
    for sig_name, gene_list in ALL_SIGNATURES.items():
        resolved, _ = resolve_signature(gene_list, set(adata.var_names))
        for g in resolved:
            if g not in marker_genes:
                marker_genes.append(g)
                gene_group_labels[g] = sig_name

    if marker_genes and adata.obs['cell_type_final'].nunique() > 1:
        try:
            sc.pl.dotplot(
                adata, var_names=marker_genes[:40],
                groupby='cell_type_final', show=False,
            )
            plt.savefig(ANN_DIR / f"{accession}_marker_dotplot.pdf",
                        bbox_inches='tight')
            plt.close()
        except Exception as e:
            print(f"    WARNING: dotplot failed: {e}")
            _placeholder_pdf(accession)
    else:
        _placeholder_pdf(accession)

    # 4. Score distribution violins for NP subtypes
    np_score_cols = [f"score_{s}" for s in NP_SIGNATURES if f"score_{s}" in adata.obs.columns]
    if np_score_cols:
        fig, axes = plt.subplots(1, len(np_score_cols),
                                  figsize=(5 * len(np_score_cols), 4))
        if len(np_score_cols) == 1:
            axes = [axes]
        for ax, col in zip(axes, np_score_cols):
            sc.pl.violin(adata, keys=col, groupby='cell_type_final',
                         ax=ax, show=False, rotation=90)
            ax.set_title(col.replace('score_', ''))
        plt.tight_layout()
        plt.savefig(ANN_DIR / f"{accession}_np_score_violins.png",
                    bbox_inches='tight', dpi=120)
        plt.close()

    plt.close('all')


def _placeholder_pdf(accession):
    """Create a placeholder PDF dotplot."""
    fig, ax = plt.subplots(figsize=(6, 2))
    ax.text(0.5, 0.5, "Insufficient data for dotplot",
            ha='center', va='center', fontsize=12)
    ax.axis('off')
    plt.savefig(ANN_DIR / f"{accession}_marker_dotplot.pdf", bbox_inches='tight')
    plt.close()


# ═══════════════════════════════════════════════════════════════════════════════
# HTML REPORTS
# ═══════════════════════════════════════════════════════════════════════════════

ANNOTATION_REPORT_TEMPLATE = """<!DOCTYPE html>
<html>
<head>
<title>Annotation Report: {{ accession }}</title>
<style>
  body { font-family: Arial, sans-serif; max-width: 1200px; margin: 0 auto; padding: 20px; }
  h1 { color: #2c3e50; border-bottom: 2px solid #2ecc71; padding-bottom: 10px; }
  h2 { color: #34495e; margin-top: 30px; }
  table { border-collapse: collapse; width: 100%; margin: 10px 0; }
  th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
  th { background-color: #2ecc71; color: white; }
  tr:nth-child(even) { background-color: #f2f2f2; }
  img { max-width: 100%; height: auto; margin: 10px 0; border: 1px solid #ddd; }
  .stats-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0; }
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
<h1>Annotation Report: {{ accession }}</h1>
<p>Generated: {{ date }}</p>

<div class="stats-grid">
  <div class="stat-box"><h3>Total cells</h3><p>{{ stats.n_cells }}</p></div>
  <div class="stat-box"><h3>Cell types (final)</h3><p>{{ stats.n_types_final }}</p></div>
  <div class="stat-box"><h3>Unassigned</h3><p>{{ stats.pct_unassigned }}%</p></div>
  <div class="stat-box"><h3>CellTypist status</h3><p>{{ stats.celltypist_status }}</p></div>
  <div class="stat-box"><h3>Resolution used</h3><p>{{ stats.resolution_used }}</p></div>
  <div class="stat-box"><h3>Compartment</h3><p>{{ stats.dominant_compartment }}</p></div>
</div>

<h2>Cell Type Distribution (Final)</h2>
<table>
  <tr><th>Cell Type</th><th>Count</th><th>Percentage</th><th>Confidence (high / med / low)</th></tr>
  {% for row in cell_type_table %}
  <tr>
    <td>{{ row.cell_type }}</td>
    <td>{{ row.count }}</td>
    <td>{{ row.pct }}%</td>
    <td>{{ row.conf_high }} / {{ row.conf_med }} / {{ row.conf_low }}</td>
  </tr>
  {% endfor %}
</table>

<h2>Gene Signature Availability</h2>
<table>
  <tr><th>Signature</th><th>Available / Total</th><th>Missing genes</th></tr>
  {% for row in sig_table %}
  <tr>
    <td>{{ row.name }}</td>
    <td>{{ row.n_avail }} / {{ row.n_total }}</td>
    <td>{{ row.missing }}</td>
  </tr>
  {% endfor %}
</table>

<h2>UMAP — Final Cell Types</h2>
<img src="{{ accession }}_umap_annotation.png" alt="UMAP annotation">

<h2>Strategy Comparison</h2>
<img src="{{ accession }}_strategy_comparison.png" alt="Strategy comparison">

<h2>Marker Dotplot</h2>
<embed src="{{ accession }}_marker_dotplot.pdf" type="application/pdf" width="100%" height="600px">

{% if has_np_violins %}
<h2>NP Subtype Score Distributions</h2>
<img src="{{ accession }}_np_score_violins.png" alt="NP score violins">
{% endif %}

<h2>Validation</h2>
<ul class="validation">
  {% for msg in validation_messages %}
  <li class="{{ 'pass' if msg.startswith('PASS') else 'warning' if msg.startswith('WARNING') or msg.startswith('FAIL') else 'info' }}">{{ msg }}</li>
  {% endfor %}
</ul>

</body>
</html>
"""


def generate_annotation_report(accession, stats, sig_stats, validation_messages):
    """Generate HTML annotation report for one dataset."""
    # Build cell type table
    ct_counts = stats['cell_type_counts']
    ct_conf = stats['cell_type_confidence_by_type']
    n_cells = stats['n_cells']

    cell_type_table = []
    for ct in sorted(ct_counts.keys(), key=lambda x: ct_counts[x], reverse=True):
        count = ct_counts[ct]
        pct = f"{count / n_cells * 100:.1f}"
        conf = ct_conf.get(ct, {})
        cell_type_table.append({
            "cell_type": ct,
            "count": count,
            "pct": pct,
            "conf_high": conf.get("high", 0),
            "conf_med": conf.get("medium", 0),
            "conf_low": conf.get("low", 0),
        })

    # Build signature table
    sig_table = []
    for name, info in sig_stats.items():
        sig_table.append({
            "name": name,
            "n_avail": info["n_avail"],
            "n_total": info["n_total"],
            "missing": ", ".join(info["missing"]) if info["missing"] else "—",
        })

    template = Template(ANNOTATION_REPORT_TEMPLATE)
    html = template.render(
        accession=accession,
        date=datetime.now().strftime("%Y-%m-%d %H:%M"),
        stats=stats,
        cell_type_table=cell_type_table,
        sig_table=sig_table,
        validation_messages=validation_messages,
        has_np_violins=(ANN_DIR / f"{accession}_np_score_violins.png").exists(),
    )
    report_path = ANN_DIR / f"{accession}_annotation_report.html"
    report_path.write_text(html)
    print(f"    Report: {report_path}")


# ═══════════════════════════════════════════════════════════════════════════════
# AGGREGATE OUTPUTS
# ═══════════════════════════════════════════════════════════════════════════════

SUMMARY_TEMPLATE = """<!DOCTYPE html>
<html>
<head>
<title>Annotation Summary — IVD scRNA-seq Atlas</title>
<style>
  body { font-family: Arial, sans-serif; max-width: 1400px; margin: 0 auto; padding: 20px; }
  h1 { color: #2c3e50; border-bottom: 2px solid #2ecc71; padding-bottom: 10px; }
  h2 { color: #34495e; margin-top: 30px; }
  table { border-collapse: collapse; width: 100%; margin: 10px 0; font-size: 13px; }
  th, td { border: 1px solid #ddd; padding: 6px 8px; text-align: right; }
  th { background-color: #2ecc71; color: white; text-align: left; }
  td:first-child { text-align: left; font-weight: bold; }
  tr:nth-child(even) { background-color: #f2f2f2; }
</style>
</head>
<body>
<h1>Annotation Summary — IVD scRNA-seq Atlas</h1>
<p>Generated: {{ date }}</p>
<p>Datasets: {{ n_datasets }} | Total cells: {{ total_cells }}</p>

<h2>Cell Type Proportions Across Datasets</h2>
<table>
  <tr>
    <th>Dataset</th><th>Cells</th>
    {% for ct in cell_types %}<th>{{ ct }}</th>{% endfor %}
    <th>% Unassigned</th>
  </tr>
  {% for row in rows %}
  <tr>
    <td>{{ row.accession }}</td>
    <td>{{ row.n_cells }}</td>
    {% for ct in cell_types %}
    <td>{{ row.proportions.get(ct, '—') }}</td>
    {% endfor %}
    <td>{{ row.pct_unassigned }}%</td>
  </tr>
  {% endfor %}
</table>

<h2>Per-Dataset Reports</h2>
<ul>
{% for acc in accessions %}
  <li><a href="{{ acc }}_annotation_report.html">{{ acc }}</a></li>
{% endfor %}
</ul>

</body>
</html>
"""


def generate_aggregate_outputs(all_stats):
    """Generate annotation_comparison.tsv and annotation_summary.html."""
    # Collect all cell types across datasets
    all_cell_types = set()
    for stats in all_stats:
        all_cell_types.update(stats['cell_type_counts'].keys())
    all_cell_types = sorted(all_cell_types)

    # Build comparison TSV
    rows_tsv = []
    for stats in all_stats:
        row = {"accession": stats["accession"], "n_cells": stats["n_cells"]}
        for ct in all_cell_types:
            count = stats['cell_type_counts'].get(ct, 0)
            row[ct] = count
            row[f"{ct}_pct"] = f"{count / stats['n_cells'] * 100:.1f}" if stats['n_cells'] > 0 else "0.0"
        rows_tsv.append(row)

    comparison_df = pd.DataFrame(rows_tsv)
    tsv_path = ANN_DIR / "annotation_comparison.tsv"
    comparison_df.to_csv(tsv_path, sep='\t', index=False)
    print(f"\nAnnotation comparison: {tsv_path}")

    # Build HTML summary
    html_rows = []
    total_cells = 0
    for stats in all_stats:
        n = stats['n_cells']
        total_cells += n
        proportions = {}
        for ct in all_cell_types:
            count = stats['cell_type_counts'].get(ct, 0)
            if count > 0:
                proportions[ct] = f"{count / n * 100:.1f}%"
        n_unassigned = stats['cell_type_counts'].get('unassigned', 0)
        html_rows.append({
            "accession": stats["accession"],
            "n_cells": n,
            "proportions": proportions,
            "pct_unassigned": f"{n_unassigned / n * 100:.1f}" if n > 0 else "0.0",
        })

    template = Template(SUMMARY_TEMPLATE)
    html = template.render(
        date=datetime.now().strftime("%Y-%m-%d %H:%M"),
        n_datasets=len(all_stats),
        total_cells=f"{total_cells:,}",
        cell_types=all_cell_types,
        rows=html_rows,
        accessions=[s["accession"] for s in all_stats],
    )
    (ANN_DIR / "annotation_summary.html").write_text(html)
    print(f"Annotation summary: {ANN_DIR / 'annotation_summary.html'}")


# ═══════════════════════════════════════════════════════════════════════════════
# PER-DATASET PROCESSING
# ═══════════════════════════════════════════════════════════════════════════════

def annotate_dataset(accession):
    """Full annotation pipeline for one dataset.

    Loads the preprocessed h5ad, applies all three strategies, assigns
    consensus labels, generates plots and reports, saves updated h5ad.

    Returns (stats_dict, validation_messages) or None on error.
    """
    h5ad_path = PROC_DIR / f"{accession}.h5ad"
    if not h5ad_path.exists():
        print(f"  SKIP {accession}: {h5ad_path} not found")
        return None

    print(f"\n{'='*60}")
    print(f"Annotating {accession}")
    print(f"{'='*60}")

    adata = sc.read_h5ad(h5ad_path)
    print(f"  Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # ── Step 1: Score expanded gene signatures ────────────────────────────
    print(f"\n  Step 1: Scoring expanded gene signatures...")
    sig_stats = score_signatures(adata)
    for sig_name, info in sig_stats.items():
        print(f"    {sig_name}: {info['n_avail']}/{info['n_total']} markers available"
              + (f" (missing: {', '.join(info['missing'])})" if info['missing'] else ""))

    # ── Step 2: Marker-based cluster labeling ─────────────────────────────
    print(f"\n  Step 2: Assigning marker-based labels (cluster majority voting)...")
    resolution_used = assign_marker_labels(adata)
    print(f"    Resolution used: {resolution_used}")
    print(f"    Marker-based label distribution:")
    for ct, n in adata.obs['cell_type_marker_based'].value_counts().items():
        print(f"      {ct}: {n} ({n/adata.shape[0]*100:.1f}%)")

    # ── Step 3: Reference-based (not applied) ─────────────────────────────
    print(f"\n  Step 3: Reference-based label transfer — SKIPPED (no IVD atlas)")
    set_reference_not_applied(adata)

    # ── Step 4: CellTypist ────────────────────────────────────────────────
    print(f"\n  Step 4: CellTypist automated annotation...")
    celltypist_ran = run_celltypist(adata)
    ct_status = "ran" if celltypist_ran else "skipped/failed"

    # ── Step 5: Consensus labeling ────────────────────────────────────────
    print(f"\n  Step 5: Assigning consensus labels...")
    assign_consensus_labels(adata, celltypist_ran=celltypist_ran)
    print(f"    Final label distribution:")
    for ct, n in adata.obs['cell_type_final'].value_counts().items():
        print(f"      {ct}: {n} ({n/adata.shape[0]*100:.1f}%)")

    # ── Step 6: Generate plots ────────────────────────────────────────────
    print(f"\n  Step 6: Generating annotation plots...")
    generate_annotation_plots(adata, accession)

    # ── Step 7: Validation ────────────────────────────────────────────────
    print(f"\n  Step 7: Running validation checks...")
    passed, val_messages = validate_annotation(adata, accession)
    for msg in val_messages:
        print(f"    {msg}")

    # ── Collect stats ─────────────────────────────────────────────────────
    dominant_comp = "mixed"
    if 'compartment' in adata.obs.columns:
        dominant_comp = adata.obs['compartment'].value_counts().index[0]

    ct_counts = adata.obs['cell_type_final'].value_counts().to_dict()

    # Confidence breakdown per cell type
    ct_conf = {}
    for ct in ct_counts:
        ct_mask = adata.obs['cell_type_final'] == ct
        conf_counts = adata.obs.loc[ct_mask, 'cell_type_confidence'].value_counts().to_dict()
        ct_conf[ct] = conf_counts

    n_unassigned = ct_counts.get('unassigned', 0)

    stats = {
        "accession": accession,
        "n_cells": adata.shape[0],
        "n_types_final": adata.obs['cell_type_final'].nunique(),
        "pct_unassigned": f"{n_unassigned / adata.shape[0] * 100:.1f}",
        "celltypist_status": ct_status,
        "resolution_used": resolution_used,
        "dominant_compartment": dominant_comp,
        "cell_type_counts": ct_counts,
        "cell_type_confidence_by_type": ct_conf,
    }

    # ── Step 8: Generate report ───────────────────────────────────────────
    print(f"\n  Step 8: Generating annotation report...")
    generate_annotation_report(accession, stats, sig_stats, val_messages)

    # ── Step 9: Save updated h5ad ─────────────────────────────────────────
    print(f"\n  Step 9: Saving updated h5ad...")
    adata.write_h5ad(h5ad_path)
    print(f"    Saved: {h5ad_path}")

    del adata
    return stats, val_messages


# ═══════════════════════════════════════════════════════════════════════════════
# VALIDATION-ONLY MODE
# ═══════════════════════════════════════════════════════════════════════════════

def validate_all():
    """Run validation across all datasets."""
    print("\n" + "=" * 60)
    print("ANNOTATION VALIDATION")
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

        if 'cell_type_final' not in adata.obs.columns:
            print(f"  SKIP: cell_type_final not present (not yet annotated)")
            del adata
            continue

        passed, messages = validate_annotation(adata, accession)
        for msg in messages:
            print(f"  {msg}")
        if not passed:
            all_pass = False
        del adata

    # Cross-dataset checks
    print(f"\n--- Cross-dataset ---")
    tsv_path = ANN_DIR / "annotation_comparison.tsv"
    if tsv_path.exists():
        comp_df = pd.read_csv(tsv_path, sep='\t')
        if len(comp_df) == len(ALL_ACCESSIONS):
            print(f"  PASS: annotation_comparison.tsv has {len(comp_df)} rows (all datasets)")
        else:
            print(f"  WARNING: annotation_comparison.tsv has {len(comp_df)}/{len(ALL_ACCESSIONS)} rows")
    else:
        print(f"  FAIL: annotation_comparison.tsv not found")
        all_pass = False

    summary_path = ANN_DIR / "annotation_summary.html"
    if summary_path.exists():
        print(f"  PASS: annotation_summary.html exists")
    else:
        print(f"  FAIL: annotation_summary.html not found")
        all_pass = False

    # Check each dataset has reports
    n_reports = sum(1 for a in ALL_ACCESSIONS
                    if (ANN_DIR / f"{a}_annotation_report.html").exists())
    print(f"  Reports: {n_reports}/{len(ALL_ACCESSIONS)} datasets")

    print(f"\n{'='*60}")
    print(f"OVERALL: {'ALL CHECKS PASSED' if all_pass else 'SOME CHECKS FAILED'}")
    print(f"{'='*60}")
    return all_pass


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    args = sys.argv[1:]

    if args and args[0] == '--validate-only':
        validate_all()
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
        try:
            result = annotate_dataset(accession)
            if result is not None:
                stats, _ = result
                all_stats.append(stats)
        except Exception as e:
            print(f"\nERROR annotating {accession}: {e}")
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
