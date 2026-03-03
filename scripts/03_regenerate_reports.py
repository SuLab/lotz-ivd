#!/usr/bin/env python3
"""Regenerate Module 03 reporting artifacts from existing h5ad files.

The preprocessing script (03_preprocessing.py) produced complete h5ad files for
all 12 datasets, but the reporting artifacts (QC reports, plots, marker tables,
overview) were not persisted.  This script regenerates them from the h5ad files
without re-running the expensive preprocessing pipeline.

Usage:
    python3 scripts/03_regenerate_reports.py
"""

import sys
import os
import warnings
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

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent.parent
PROC_DIR = BASE / "data" / "processed"
QC_DIR = BASE / "results" / "qc_reports"
ANN_DIR = BASE / "results" / "annotations"
META_FILE = BASE / "metadata" / "sample_metadata.tsv"

QC_DIR.mkdir(parents=True, exist_ok=True)
ANN_DIR.mkdir(parents=True, exist_ok=True)

ALL_ACCESSIONS = [
    "GSE160756", "GSE165722", "GSE189916", "GSE199866", "GSE205535",
    "CNP0002664", "GSE233666", "GSE244889", "GSE251686", "GSE255768",
    "GSE230809", "GSE242443",
]

CELL_TYPE_MARKERS = {
    "Chondrocyte_NP":  ["ACAN", "COL2A1", "SOX9", "KRT19"],
    "Fibroblast_AF":   ["COL1A1", "COL1A2", "THY1", "SCX"],
    "Endothelial":     ["PECAM1", "VWF", "CDH5", "FLT1"],
    "Macrophage":      ["CD68", "CD163", "CSF1R", "MRC1"],
    "T_cell":          ["CD3D", "CD3E", "CD4", "CD8A"],
    "B_cell":          ["CD79A", "MS4A1", "CD19"],
    "Mast_cell":       ["KIT", "TPSAB1", "CPA3"],
    "Pericyte_SMC":    ["ACTA2", "TAGLN", "MYH11", "RGS5"],
    "Notochordal":     ["T", "SHH", "NOG"],
    "Progenitor":      ["CD44", "PROM1", "NES"],
}

LEIDEN_RESOLUTIONS = [0.2, 0.5, 0.8, 1.0, 1.5]
WORKING_RESOLUTION = 0.5


# ═════════════════════════════════════════════════════════════════════════════
# RECONSTRUCT QC STATS FROM H5AD + METADATA
# ═════════════════════════════════════════════════════════════════════════════

def reconstruct_qc_stats(adata, accession, metadata_df):
    """Reconstruct the qc_stats dict from a processed h5ad and metadata.

    Note: metadata n_cells_raw are paper-reported counts, which may differ from
    actual raw barcode counts in the GEO files.  When a sample's post-QC count
    exceeds its metadata raw count (i.e. the raw files had more barcodes than
    the paper reported), we use the post-QC count as the raw count floor.
    """
    stats = {"accession": accession, "samples": {}}

    working_key = f"leiden_res_{WORKING_RESOLUTION}"

    # Per-sample stats
    total_raw = 0
    for sample_id in adata.obs['sample_id'].unique():
        sample_cells = int((adata.obs['sample_id'] == sample_id).sum())
        meta_row = metadata_df[metadata_df['sample_id'] == sample_id]

        cells_raw = sample_cells  # default: use post-QC count
        if len(meta_row) > 0:
            raw_val = meta_row['n_cells_raw'].values[0]
            if pd.notna(raw_val):
                meta_raw = int(float(raw_val))
                # Use whichever is larger: metadata raw or post-QC count
                cells_raw = max(meta_raw, sample_cells)

        total_raw += cells_raw
        stats["samples"][sample_id] = {
            "cells_raw": cells_raw,
            "cells_after_qc_filter": sample_cells,
            "cells_after_doublet": sample_cells,
        }

    stats["cells_raw_total"] = total_raw
    stats["cells_after_qc"] = adata.shape[0]
    stats["genes_raw"] = adata.shape[1]  # approximate (post-gene-filtering)
    stats["genes_after_qc"] = adata.shape[1]
    stats["n_hvgs"] = int(adata.var['highly_variable'].sum()) if 'highly_variable' in adata.var else 0
    stats["n_pcs"] = adata.obsm['X_pca'].shape[1] if 'X_pca' in adata.obsm else 0
    stats["n_clusters_working"] = int(adata.obs[working_key].nunique()) if working_key in adata.obs else 0
    stats["cell_type_counts"] = adata.obs['cell_type_preliminary'].value_counts().to_dict()

    return stats


# ═════════════════════════════════════════════════════════════════════════════
# MARKER GENE EXTRACTION
# ═════════════════════════════════════════════════════════════════════════════

def extract_markers(adata, accession):
    """Extract marker genes from stored rank_genes_groups and save as TSV."""
    working_key = f"leiden_res_{WORKING_RESOLUTION}"

    if 'rank_genes_groups' not in adata.uns:
        print(f"  WARNING: rank_genes_groups not in .uns, recomputing...")
        sc.tl.rank_genes_groups(adata, groupby=working_key, method='wilcoxon',
                                use_raw=False)

    markers_list = []
    for cluster in adata.obs[working_key].cat.categories:
        result = sc.get.rank_genes_groups_df(adata, group=str(cluster))
        result = result.head(50)
        result['cluster'] = cluster
        markers_list.append(result)

    markers_df = pd.concat(markers_list, ignore_index=True)
    markers_path = ANN_DIR / f"{accession}_markers.tsv"
    markers_df.to_csv(markers_path, sep='\t', index=False)
    print(f"  Markers: {markers_path} ({len(markers_df)} rows)")


# ═════════════════════════════════════════════════════════════════════════════
# QC PLOTS (from 03_preprocessing.py)
# ═════════════════════════════════════════════════════════════════════════════

def generate_qc_plots(adata, accession):
    """Generate all QC plots for one dataset."""
    plt.rcParams['figure.dpi'] = 100
    working_key = f"leiden_res_{WORKING_RESOLUTION}"

    # 1. QC violin plots per sample
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    for ax, metric, title in zip(
        axes,
        ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
        ['Genes per cell', 'UMI counts per cell', '% Mitochondrial']
    ):
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
    fig, axes = plt.subplots(1, n_res, figsize=(6 * n_res, 5))
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

    # 5. UMAP by QC metrics (pct_mt, n_genes, total_counts)
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
    all_markers = list(dict.fromkeys(all_markers))

    if all_markers:
        try:
            sc.pl.dotplot(adata, var_names=all_markers, groupby=working_key,
                          show=False, save=False)
            # dotplot creates its own figure; save it
            plt.savefig(QC_DIR / f"{accession}_dotplot_markers.png", bbox_inches='tight')
            plt.close('all')
        except Exception as e:
            print(f"  WARNING: dotplot failed ({e}), creating placeholder")
            fig, ax = plt.subplots(figsize=(6, 2))
            ax.text(0.5, 0.5, f"Dotplot failed: {e}", ha='center', va='center')
            ax.axis('off')
            plt.savefig(QC_DIR / f"{accession}_dotplot_markers.png", bbox_inches='tight')
            plt.close()
    else:
        fig, ax = plt.subplots(figsize=(6, 2))
        ax.text(0.5, 0.5, "No canonical markers found in gene set",
                ha='center', va='center', fontsize=14)
        ax.axis('off')
        plt.savefig(QC_DIR / f"{accession}_dotplot_markers.png", bbox_inches='tight')
        plt.close()

    # 7. Top markers
    try:
        sc.pl.rank_genes_groups(adata, n_genes=10, show=False, save=False)
        plt.savefig(QC_DIR / f"{accession}_top_markers.png", bbox_inches='tight')
        plt.close('all')
    except Exception as e:
        print(f"  WARNING: rank_genes_groups plot failed ({e})")
        fig, ax = plt.subplots(figsize=(6, 2))
        ax.text(0.5, 0.5, "Could not generate marker heatmap", ha='center', va='center')
        ax.axis('off')
        plt.savefig(QC_DIR / f"{accession}_top_markers.png", bbox_inches='tight')
        plt.close()

    print(f"  Plots saved to {QC_DIR}")


# ═════════════════════════════════════════════════════════════════════════════
# HTML REPORTS (templates from 03_preprocessing.py)
# ═════════════════════════════════════════════════════════════════════════════

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
          <span class="warning"> LOW</span>
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
<h1>Cross-Dataset QC Overview &mdash; IVD scRNA-seq Atlas</h1>
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
    print(f"  HTML report: {report_path}")


def generate_overview(all_stats):
    """Generate cross-dataset overview: qc_summary.tsv, overview plots, HTML."""
    summary_rows = []
    for stats in all_stats:
        n_samples = len(stats["samples"])
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
    print(f"\nqc_summary.tsv: {QC_DIR / 'qc_summary.tsv'}")

    # Overview plots
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    x = range(len(summary_df))
    ax = axes[0]
    ax.bar(x, summary_df['cells_raw_total'], alpha=0.5, label='Raw')
    ax.bar(x, summary_df['cells_after_qc'], alpha=0.8, label='After QC')
    ax.set_xticks(list(x))
    ax.set_xticklabels(summary_df['accession'], rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Cell count')
    ax.set_title('Cell Counts: Raw vs After QC')
    ax.legend()

    ax = axes[1]
    colors = ['#e74c3c' if r < 50 else '#f39c12' if r < 70 else '#27ae60'
              for r in summary_df['retention_pct']]
    ax.bar(x, summary_df['retention_pct'], color=colors)
    ax.axhline(y=50, color='red', linestyle='--', alpha=0.5, label='50% threshold')
    ax.set_xticks(list(x))
    ax.set_xticklabels(summary_df['accession'], rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Retention %')
    ax.set_title('Cell Retention After QC')
    ax.set_ylim(0, 100)
    ax.legend()

    plt.tight_layout()
    plt.savefig(QC_DIR / "overview_cells_comparison.png", bbox_inches='tight')
    plt.close()

    # Clusters comparison
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(list(x), summary_df['n_clusters_working'], color='#3498db')
    ax.set_xticks(list(x))
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
    print(f"qc_overview.html: {QC_DIR / 'qc_overview.html'}")


# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════

def main():
    metadata_df = pd.read_csv(META_FILE, sep='\t')
    all_stats = []

    for accession in ALL_ACCESSIONS:
        h5ad_path = PROC_DIR / f"{accession}.h5ad"
        if not h5ad_path.exists():
            print(f"\nSKIP {accession}: {h5ad_path} not found")
            continue

        print(f"\n{'='*60}")
        print(f"Regenerating artifacts for {accession}")
        print(f"{'='*60}")

        adata = sc.read_h5ad(h5ad_path)
        print(f"  Loaded: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

        # 1. Reconstruct QC stats
        qc_stats = reconstruct_qc_stats(adata, accession, metadata_df)
        all_stats.append(qc_stats)

        # 2. Extract marker genes to TSV
        extract_markers(adata, accession)

        # 3. Generate QC plots
        generate_qc_plots(adata, accession)

        # 4. Generate HTML report
        generate_qc_report(accession, qc_stats)

        del adata

    # 5. Generate cross-dataset overview
    if all_stats:
        generate_overview(all_stats)

    print(f"\n{'='*60}")
    print("All artifacts regenerated.")
    print(f"  QC reports: {QC_DIR}")
    print(f"  Marker tables: {ANN_DIR}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
