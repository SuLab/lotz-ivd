#!/usr/bin/env python3
"""Module 04: Coarse cell classification (mesenchymal vs non-mesenchymal).

Classifies cells into two broad categories using unambiguous marker genes.
This coarse classification enables tiered integration in Module 05.
Fine-grained cell type annotation happens after integration (in Module 05).

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
# Non-mesenchymal markers: any of these expressed → non-mesenchymal candidate
NON_MESENCHYMAL_MARKERS = {
    "immune": {
        "pan_immune":  ["PTPRC"],       # CD45
        "t_cell":      ["CD3D", "CD3E"],
        "macrophage":  ["CD68", "CD14", "CSF1R"],
        "b_cell":      ["CD79A", "MS4A1"],
        "mast_cell":   ["KIT", "TPSAB1"],
        "nk_cell":     ["NKG7", "GNLY"],
    },
    "endothelial": {
        "endothelial": ["PECAM1", "VWF", "CDH5"],
    },
    "pericyte": {
        "pericyte":    ["ACTA2", "RGS5", "PDGFRB"],
    },
}

# Flat list for scoring
NON_MESENCHYMAL_GENES = [
    "PTPRC", "CD3D", "CD3E", "CD68", "CD14", "CSF1R", "CD79A", "MS4A1",
    "KIT", "TPSAB1", "NKG7", "GNLY", "PECAM1", "VWF", "CDH5",
    "ACTA2", "RGS5", "PDGFRB",
]

# Mesenchymal markers: expected in IVD disc cells
MESENCHYMAL_GENES = ["COL2A1", "COL1A1", "ACAN", "SOX9", "DCN", "LUM", "VCAN"]


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


def score_cell_class(adata):
    """Score each cell for mesenchymal and non-mesenchymal marker expression.

    Computes fraction of markers detected above background for each set.
    Stores scores in obs['score_mesenchymal'] and obs['score_non_mesenchymal'].
    """
    var_names = set(adata.var_names)

    # Score mesenchymal markers
    mes_resolved, mes_missing = resolve_gene_list(MESENCHYMAL_GENES, var_names)
    if mes_resolved:
        sc.tl.score_genes(adata, gene_list=mes_resolved, score_name='score_mesenchymal')
    else:
        adata.obs['score_mesenchymal'] = 0.0
    print(f"    Mesenchymal markers: {len(mes_resolved)}/{len(MESENCHYMAL_GENES)} available"
          + (f" (missing: {', '.join(mes_missing)})" if mes_missing else ""))

    # Score non-mesenchymal markers
    non_mes_resolved, non_mes_missing = resolve_gene_list(NON_MESENCHYMAL_GENES, var_names)
    if non_mes_resolved:
        sc.tl.score_genes(adata, gene_list=non_mes_resolved, score_name='score_non_mesenchymal')
    else:
        adata.obs['score_non_mesenchymal'] = 0.0
    print(f"    Non-mesenchymal markers: {len(non_mes_resolved)}/{len(NON_MESENCHYMAL_GENES)} available"
          + (f" (missing: {', '.join(non_mes_missing)})" if non_mes_missing else ""))

    return {
        "mesenchymal": {"available": mes_resolved, "missing": mes_missing},
        "non_mesenchymal": {"available": non_mes_resolved, "missing": non_mes_missing},
    }


def _compute_expression_threshold(expr):
    """Compute a per-gene threshold for 'expressed' status.

    Uses the mean of nonzero values as a guide. Cells above 10% of that mean
    are considered to have meaningful expression. This avoids treating
    normalization artifacts or ambient RNA as true expression.
    """
    nonzero = expr[expr > 0]
    if len(nonzero) == 0:
        return 0.0
    return max(nonzero.mean() * 0.1, 0.1)


def classify_cells(adata):
    """Assign cell_class based on marker scores, with caveats for CD68 and ACTA2.

    Primary classification uses sc.tl.score_genes() scores (mesenchymal vs
    non-mesenchymal). Strong single-marker signals (PTPRC, PECAM1) override
    scores when expression is well above background.

    Caveats:
    - CD68 alone doesn't mean immune (expressed in stressed IVD cells)
    - ACTA2 alone doesn't mean pericyte (expressed in myofibroblasts)
    """
    n_cells = adata.shape[0]
    labels = np.full(n_cells, "ambiguous", dtype=object)

    score_mes = adata.obs['score_mesenchymal'].values
    score_non = adata.obs['score_non_mesenchymal'].values

    # Load key marker expression vectors and compute thresholds
    marker_data = {}
    for gene in ['PTPRC', 'PECAM1', 'CD68', 'CD14', 'ACTA2', 'RGS5', 'PDGFRB']:
        expr = _get_gene_expression(adata, gene)
        if expr is not None:
            thresh = _compute_expression_threshold(expr)
            marker_data[gene] = (expr, thresh)

    def _is_expressed(gene, i):
        if gene not in marker_data:
            return False
        expr, thresh = marker_data[gene]
        return expr[i] > thresh

    for i in range(n_cells):
        s_mes = score_mes[i]
        s_non = score_non[i]

        # Strong PTPRC expression → non-mesenchymal (immune)
        if _is_expressed('PTPRC', i):
            labels[i] = "non_mesenchymal"
            continue

        # Strong PECAM1 expression → non-mesenchymal (endothelial)
        if _is_expressed('PECAM1', i):
            labels[i] = "non_mesenchymal"
            continue

        # CD68 caveat: in IVD cells, CD68 is expressed at low levels in
        # stressed disc cells. Require PTPRC or CD14 co-expression.
        if _is_expressed('CD68', i):
            if _is_expressed('PTPRC', i) or _is_expressed('CD14', i):
                labels[i] = "non_mesenchymal"
                continue
            # CD68 alone — fall through to score-based classification

        # ACTA2 caveat: degenerated disc cells may express ACTA2 (myofibroblasts).
        # Classify as mesenchymal unless co-expressed with pericyte markers.
        if _is_expressed('ACTA2', i):
            if _is_expressed('RGS5', i) or _is_expressed('PDGFRB', i):
                labels[i] = "non_mesenchymal"
                continue
            # ACTA2 alone — likely myofibroblast, classify as mesenchymal
            labels[i] = "mesenchymal"
            continue

        # Score-based classification
        if s_mes > s_non:
            labels[i] = "mesenchymal"
        elif s_non > s_mes:
            labels[i] = "non_mesenchymal"
        # else: both near zero or tied → stays "ambiguous"

    adata.obs['cell_class_raw'] = labels


def cluster_majority_vote(adata, resolution=CLASSIFICATION_RESOLUTION):
    """Assign cell_class by cluster-level majority voting to reduce noise.

    1. Compute Leiden clusters at the given resolution
    2. For each cluster, compute the majority cell_class
    3. Assign majority class to all cells in the cluster
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

    labels = adata.obs['cell_class_raw'].values.copy()
    final_labels = np.empty(adata.shape[0], dtype=object)
    mixed_clusters = []

    for cluster in adata.obs[cluster_key].unique():
        mask = adata.obs[cluster_key] == cluster
        cluster_labels = labels[mask]

        # Count non-ambiguous votes
        non_amb = cluster_labels[cluster_labels != "ambiguous"]
        if len(non_amb) == 0:
            # All ambiguous — assign based on scores
            final_labels[mask] = "ambiguous"
            continue

        counts = pd.Series(non_amb).value_counts()
        majority_class = counts.index[0]
        majority_pct = counts.iloc[0] / len(non_amb) * 100

        final_labels[mask] = majority_class

        if majority_pct < 70:
            mixed_clusters.append({
                "cluster": cluster,
                "majority_class": majority_class,
                "majority_pct": f"{majority_pct:.1f}",
                "n_cells": mask.sum(),
            })

    adata.obs['cell_class'] = final_labels

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

    # Check 1: All cells have cell_class
    if 'cell_class' not in adata.obs.columns:
        messages.append("FAIL: cell_class column missing")
        return False, messages

    valid_classes = {"mesenchymal", "non_mesenchymal", "ambiguous"}
    cell_classes = adata.obs['cell_class']
    invalid = ~cell_classes.isin(valid_classes)
    if invalid.sum() > 0:
        messages.append(f"FAIL: {invalid.sum()} cells have invalid cell_class values")
        all_pass = False
    else:
        messages.append(f"PASS: All {n_cells} cells have valid cell_class labels")

    # Check 2: <5% ambiguous per dataset
    n_ambiguous = (cell_classes == 'ambiguous').sum()
    pct_ambiguous = n_ambiguous / n_cells * 100
    if pct_ambiguous < 5:
        messages.append(f"PASS: {pct_ambiguous:.1f}% ambiguous (< 5%)")
    else:
        messages.append(f"FAIL: {pct_ambiguous:.1f}% ambiguous (>= 5%)")
        all_pass = False

    # Check 3: PTPRC enriched in non_mesenchymal vs mesenchymal (>10-fold)
    if 'PTPRC' in adata.var_names:
        ptprc = _get_gene_expression(adata, 'PTPRC')
        mes_mask = cell_classes == 'mesenchymal'
        non_mes_mask = cell_classes == 'non_mesenchymal'

        if mes_mask.sum() > 0 and non_mes_mask.sum() > 0:
            mean_mes = ptprc[mes_mask].mean() + 1e-10
            mean_non = ptprc[non_mes_mask].mean() + 1e-10
            fold = mean_non / mean_mes

            if fold > 10:
                messages.append(f"PASS: PTPRC {fold:.0f}-fold enriched in non_mesenchymal")
            else:
                messages.append(f"WARNING: PTPRC only {fold:.1f}-fold enriched in "
                                f"non_mesenchymal (expected >10-fold)")
        else:
            messages.append("INFO: Cannot compute PTPRC enrichment "
                            "(missing mesenchymal or non_mesenchymal cells)")
    else:
        messages.append("INFO: PTPRC not in var_names")

    # Check 4: COL2A1 or COL1A1 enriched in mesenchymal
    for gene in ["COL2A1", "COL1A1"]:
        if gene not in adata.var_names:
            continue
        expr = _get_gene_expression(adata, gene)
        mes_mask = cell_classes == 'mesenchymal'
        non_mes_mask = cell_classes == 'non_mesenchymal'

        if mes_mask.sum() > 0 and non_mes_mask.sum() > 0:
            mean_mes = expr[mes_mask].mean() + 1e-10
            mean_non = expr[non_mes_mask].mean() + 1e-10
            fold = mean_mes / mean_non

            if fold > 2:
                messages.append(f"PASS: {gene} {fold:.1f}-fold enriched in mesenchymal")
            else:
                messages.append(f"INFO: {gene} only {fold:.1f}-fold enriched in mesenchymal")
        break  # Only need one of COL2A1/COL1A1

    # Check 5: Classification report exists
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

    # 1. UMAP colored by cell_class
    fig, ax = plt.subplots(figsize=(10, 8))
    palette = {"mesenchymal": "#2196F3", "non_mesenchymal": "#E91E63", "ambiguous": "#9E9E9E"}
    sc.pl.umap(adata, color='cell_class', ax=ax, show=False, palette=palette,
               title=f'{accession} — Cell Class')
    plt.tight_layout()
    plt.savefig(ANN_DIR / f"{accession}_umap_cell_class.png",
                bbox_inches='tight', dpi=150)
    plt.close()

    # 2. Dot plot of classification markers by cell_class
    all_markers = MESENCHYMAL_GENES + NON_MESENCHYMAL_GENES
    resolved_markers = [g for g in all_markers if resolve_gene(g, set(adata.var_names))]
    # Resolve to actual var_names
    resolved_markers = [resolve_gene(g, set(adata.var_names)) for g in all_markers
                        if resolve_gene(g, set(adata.var_names)) is not None]
    # Deduplicate preserving order
    seen = set()
    unique_markers = []
    for g in resolved_markers:
        if g not in seen:
            unique_markers.append(g)
            seen.add(g)

    if unique_markers and adata.obs['cell_class'].nunique() > 1:
        try:
            sc.pl.dotplot(
                adata, var_names=unique_markers,
                groupby='cell_class', show=False,
            )
            plt.savefig(ANN_DIR / f"{accession}_classification_dotplot.pdf",
                        bbox_inches='tight')
            plt.close()
        except Exception as e:
            print(f"    WARNING: dotplot failed: {e}")

    # 3. Bar plot: mesenchymal vs non-mesenchymal proportions
    counts = adata.obs['cell_class'].value_counts()
    fig, ax = plt.subplots(figsize=(6, 4))
    colors = [palette.get(c, '#999999') for c in counts.index]
    counts.plot(kind='bar', ax=ax, color=colors)
    ax.set_title(f'{accession} — Cell Class Proportions')
    ax.set_ylabel('Number of cells')
    ax.set_xlabel('')
    for i, (idx, val) in enumerate(counts.items()):
        ax.text(i, val + 50, f'{val/adata.shape[0]*100:.1f}%', ha='center', fontsize=9)
    plt.tight_layout()
    plt.savefig(ANN_DIR / f"{accession}_class_proportions.png",
                bbox_inches='tight', dpi=150)
    plt.close()

    # 4. Validation: marker expression distributions
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    key_markers = {
        'PTPRC': 'Non-mesenchymal (immune)',
        'COL2A1': 'Mesenchymal',
    }
    for ax, (gene, desc) in zip(axes, key_markers.items()):
        resolved = resolve_gene(gene, set(adata.var_names))
        if resolved and adata.obs['cell_class'].nunique() > 1:
            sc.pl.violin(adata, keys=resolved, groupby='cell_class',
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
<h1>Classification Report: {{ accession }}</h1>
<p>Generated: {{ date }}</p>

<div class="stats-grid">
  <div class="stat-box"><h3>Total cells</h3><p>{{ stats.n_cells }}</p></div>
  <div class="stat-box"><h3>Mesenchymal</h3><p>{{ stats.n_mesenchymal }} ({{ stats.pct_mesenchymal }}%)</p></div>
  <div class="stat-box"><h3>Non-mesenchymal</h3><p>{{ stats.n_non_mesenchymal }} ({{ stats.pct_non_mesenchymal }}%)</p></div>
  <div class="stat-box"><h3>Ambiguous</h3><p>{{ stats.n_ambiguous }} ({{ stats.pct_ambiguous }}%)</p></div>
  <div class="stat-box"><h3>Clustering resolution</h3><p>{{ stats.resolution }}</p></div>
  <div class="stat-box"><h3>Compartment</h3><p>{{ stats.dominant_compartment }}</p></div>
</div>

<h2>UMAP — Cell Class</h2>
<img src="{{ accession }}_umap_cell_class.png" alt="UMAP cell class">

<h2>Cell Class Proportions</h2>
<img src="{{ accession }}_class_proportions.png" alt="Class proportions">

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
  <tr><th>Cluster</th><th>Majority class</th><th>Majority %</th><th>Cells</th></tr>
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
    """Generate classification_summary.tsv with cell class proportions per dataset."""
    rows = []
    for stats in all_stats:
        rows.append({
            "accession": stats["accession"],
            "n_cells": stats["n_cells"],
            "n_mesenchymal": stats["n_mesenchymal"],
            "pct_mesenchymal": stats["pct_mesenchymal"],
            "n_non_mesenchymal": stats["n_non_mesenchymal"],
            "pct_non_mesenchymal": stats["pct_non_mesenchymal"],
            "n_ambiguous": stats["n_ambiguous"],
            "pct_ambiguous": stats["pct_ambiguous"],
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
    raw_counts = adata.obs['cell_class_raw'].value_counts()
    for cls, n in raw_counts.items():
        print(f"    {cls}: {n} ({n/adata.shape[0]*100:.1f}%)")

    # ── Step 3: Cluster majority voting ──────────────────────────────────
    print(f"\n  Step 3: Cluster-level majority voting...")
    mixed_clusters = cluster_majority_vote(adata)
    final_counts = adata.obs['cell_class'].value_counts()
    print(f"    Final class distribution:")
    for cls, n in final_counts.items():
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

    n_mes = int((adata.obs['cell_class'] == 'mesenchymal').sum())
    n_non = int((adata.obs['cell_class'] == 'non_mesenchymal').sum())
    n_amb = int((adata.obs['cell_class'] == 'ambiguous').sum())
    n_total = adata.shape[0]

    stats = {
        "accession": accession,
        "n_cells": n_total,
        "n_mesenchymal": n_mes,
        "pct_mesenchymal": f"{n_mes / n_total * 100:.1f}",
        "n_non_mesenchymal": n_non,
        "pct_non_mesenchymal": f"{n_non / n_total * 100:.1f}",
        "n_ambiguous": n_amb,
        "pct_ambiguous": f"{n_amb / n_total * 100:.1f}",
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
                'cell_type_marker_confidence', 'celltypist_conf_score']
    for col in old_cols:
        if col in adata.obs.columns:
            del adata.obs[col]
    # Also remove old score columns
    old_score_cols = [c for c in adata.obs.columns
                      if c.startswith('score_') and c not in
                      ('score_mesenchymal', 'score_non_mesenchymal')]
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

        if 'cell_class' not in adata.obs.columns:
            print(f"  SKIP: cell_class not present (not yet classified)")
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
