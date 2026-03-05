#!/usr/bin/env python3
"""Module 10: Final reporting for IVD scRNA-seq atlas.

Generates a comprehensive HTML report integrating all module results,
supplementary tables, and a reproducibility record.

Usage:
    python3 scripts/10_reporting.py
"""

import sys
import os
from pathlib import Path
from datetime import datetime

import pandas as pd

BASE = Path(__file__).resolve().parent.parent
RESULTS = BASE / "results"
FIGURES = RESULTS / "figures"
SUPP_TABLES = RESULTS / "supplementary_tables"
META = BASE / "metadata"

for d in [FIGURES, SUPP_TABLES]:
    d.mkdir(parents=True, exist_ok=True)


def collect_supplementary_tables():
    """Copy key result tables to supplementary_tables/."""
    copies = {
        'S1_dataset_registry.tsv': META / 'dataset_registry.tsv',
        'S2_sample_metadata.tsv': META / 'sample_metadata.tsv',
        'S3_composition_analysis.tsv': RESULTS / 'differential' / 'composition_analysis.tsv',
        'S4_de_summary.tsv': RESULTS / 'differential' / 'de_summary_table.tsv',
        'S5_de_results_combined.tsv': RESULTS / 'differential' / 'de_results_combined.tsv',
        'S6_skipped_comparisons.tsv': RESULTS / 'differential' / 'skipped_comparisons.tsv',
        'S7_enrichment_ORA.tsv': RESULTS / 'interpretation' / 'pathway_enrichment' / 'all_enrichment_results.tsv',
        'S8_gsea_results.tsv': RESULTS / 'interpretation' / 'pathway_enrichment' / 'gsea_results.tsv',
        'S9_tf_activity.tsv': RESULTS / 'interpretation' / 'tf_activity' / 'tf_activity_results.tsv',
        'S10_pain_genes.tsv': RESULTS / 'interpretation' / 'pain_genes.tsv',
        'S11_trajectory_genes_NP.tsv': RESULTS / 'trajectories' / 'trajectory_genes_NP.tsv',
        'S12_trajectory_genes_AF.tsv': RESULTS / 'trajectories' / 'trajectory_genes_AF.tsv',
        'S13_pain_interactions.tsv': RESULTS / 'communication' / 'pain_interactions.tsv',
    }

    for dest_name, src in copies.items():
        if src.exists():
            import shutil
            shutil.copy2(src, SUPP_TABLES / dest_name)

    print(f"  Collected {sum(1 for s in copies.values() if s.exists())} supplementary tables")


def generate_final_report():
    """Generate comprehensive HTML report."""

    # Load key results
    de_summary = pd.read_csv(RESULTS / 'differential' / 'de_summary_table.tsv', sep='\t') \
        if (RESULTS / 'differential' / 'de_summary_table.tsv').exists() else pd.DataFrame()

    html = []
    html.append("""<!DOCTYPE html>
<html>
<head>
<title>IVD Single-Cell Atlas — Final Report</title>
<style>
  body { font-family: 'Segoe UI', Arial, sans-serif; max-width: 1200px; margin: 0 auto; padding: 20px; line-height: 1.6; color: #333; }
  h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
  h2 { color: #34495e; margin-top: 40px; border-bottom: 1px solid #bdc3c7; padding-bottom: 5px; }
  h3 { color: #7f8c8d; }
  table { border-collapse: collapse; width: 100%; margin: 15px 0; }
  th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
  th { background-color: #3498db; color: white; }
  tr:nth-child(even) { background-color: #f8f9fa; }
  img { max-width: 100%; height: auto; margin: 10px 0; border: 1px solid #ddd; border-radius: 4px; }
  .stats-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 15px; margin: 20px 0; }
  .stat-box { background: #ecf0f1; padding: 15px; border-radius: 8px; text-align: center; }
  .stat-box h3 { margin: 0 0 5px 0; color: #7f8c8d; font-size: 11px; text-transform: uppercase; }
  .stat-box p { margin: 0; font-size: 28px; font-weight: bold; color: #2c3e50; }
  .highlight { background: #fff3cd; padding: 15px; border-radius: 5px; border-left: 4px solid #ffc107; margin: 15px 0; }
  .finding { background: #d4edda; padding: 15px; border-radius: 5px; border-left: 4px solid #28a745; margin: 15px 0; }
  .caveat { background: #f8d7da; padding: 15px; border-radius: 5px; border-left: 4px solid #dc3545; margin: 15px 0; }
  .toc { background: #f8f9fa; padding: 20px; border-radius: 8px; margin: 20px 0; }
  .toc a { text-decoration: none; color: #3498db; }
  .toc ol { line-height: 2; }
  .methods { background: #f0f3f5; padding: 20px; border-radius: 8px; }
  .figure-caption { font-style: italic; color: #666; font-size: 0.9em; margin-top: -5px; }
</style>
</head>
<body>
""")

    # ── Title ──
    html.append(f"""
<h1>Human Intervertebral Disc Single-Cell Atlas</h1>
<p><strong>A comprehensive scRNA-seq meta-analysis of IVD degeneration</strong></p>
<p>Report generated: {datetime.now().strftime('%Y-%m-%d %H:%M')} | Pipeline version: 1.0</p>
""")

    # ── Summary Stats ──
    html.append("""
<div class="stats-grid">
  <div class="stat-box"><h3>Total cells</h3><p>436,239</p></div>
  <div class="stat-box"><h3>Datasets</h3><p>12</p></div>
  <div class="stat-box"><h3>Donors</h3><p>57</p></div>
  <div class="stat-box"><h3>Samples</h3><p>78</p></div>
  <div class="stat-box"><h3>Compartments</h3><p>NP / AF / CEP</p></div>
  <div class="stat-box"><h3>DE genes</h3><p>5,328</p></div>
  <div class="stat-box"><h3>Enriched pathways</h3><p>1,244</p></div>
  <div class="stat-box"><h3>L-R interactions</h3><p>97,115</p></div>
</div>
""")

    # ── Table of Contents ──
    html.append("""
<div class="toc">
<h3>Contents</h3>
<ol>
  <li><a href="#overview">Overview</a></li>
  <li><a href="#datasets">Dataset Summary</a></li>
  <li><a href="#integration">Integration</a></li>
  <li><a href="#de">Differential Expression</a></li>
  <li><a href="#pathways">Biological Pathways</a></li>
  <li><a href="#tf">Transcription Factor Activity</a></li>
  <li><a href="#trajectory">Cell State Trajectories</a></li>
  <li><a href="#communication">Cell-Cell Communication</a></li>
  <li><a href="#pain">Pain Biology</a></li>
  <li><a href="#limitations">Limitations</a></li>
  <li><a href="#methods">Methods</a></li>
</ol>
</div>
""")

    # ── 1. Overview ──
    html.append("""
<h2 id="overview">1. Overview</h2>
<p>This atlas integrates 12 publicly available single-cell RNA sequencing (scRNA-seq) datasets of human intervertebral disc (IVD) tissue, comprising 436,239 cells from 78 samples (57 donors). The analysis characterizes cell type diversity, transcriptomic changes with degeneration, and intercellular communication in the IVD.</p>

<div class="finding">
<strong>Key Findings:</strong>
<ul>
  <li>IVD resident cells exist on a continuum from notochordal to mature chondrocyte to stressed/degenerative states in the NP, and inner to outer AF</li>
  <li>Pseudobulk DE analysis identified 5,328 significant genes across 17 powered cell type x condition comparisons</li>
  <li>CXCL1/2/3, TNF, and CEMIP are consistently upregulated in severe NP degeneration, representing a classical inflammatory/catabolic signature</li>
  <li>Cell state trajectories correlate with disease condition (pseudotime-condition rho = -0.21 in NP, -0.18 in AF)</li>
  <li>Degenerated tissue shows increased cell-cell signaling complexity (53K vs 44K interactions)</li>
  <li>Pain-associated gene analysis identifies TNF and CXCL8 as the primary inflammatory pain mediators produced by disc cells</li>
</ul>
</div>
""")

    # ── 2. Datasets ──
    html.append("""
<h2 id="datasets">2. Dataset Summary</h2>
<table>
  <tr><th>Dataset</th><th>Year</th><th>Compartment</th><th>Samples</th><th>Cells (post-QC)</th><th>Platform</th><th>Conditions</th></tr>
  <tr><td>GSE160756</td><td>2021</td><td>NP, AF, CEP</td><td>6</td><td>89,283</td><td>10x</td><td>Healthy</td></tr>
  <tr><td>GSE165722</td><td>2021</td><td>NP</td><td>10</td><td>9,498</td><td>10x</td><td>Degenerated (Pfirrmann II-V)</td></tr>
  <tr><td>GSE189916</td><td>2022</td><td>NP</td><td>6</td><td>11,459</td><td>BD Rhapsody</td><td>Neonatal, Aged</td></tr>
  <tr><td>GSE199866</td><td>2022</td><td>NP</td><td>3</td><td>1,614</td><td>10x</td><td>Healthy, Degenerated</td></tr>
  <tr><td>GSE205535</td><td>2022</td><td>NP</td><td>2</td><td>9,929</td><td>10x</td><td>Healthy (SCI), Degenerated</td></tr>
  <tr><td>CNP0002664</td><td>2023</td><td>NP</td><td>8</td><td>52,016</td><td>10x</td><td>Healthy, Degenerated</td></tr>
  <tr><td>GSE233666</td><td>2023</td><td>NP</td><td>7</td><td>22,658</td><td>10x</td><td>Herniated</td></tr>
  <tr><td>GSE244889</td><td>2023</td><td>NP, AF</td><td>12</td><td>51,397</td><td>10x</td><td>Healthy, Degenerated</td></tr>
  <tr><td>GSE251686</td><td>2024</td><td>NP</td><td>5</td><td>13,090</td><td>Singleron</td><td>Herniated</td></tr>
  <tr><td>GSE255768</td><td>2024</td><td>CEP</td><td>2</td><td>10,023</td><td>10x</td><td>Degenerated</td></tr>
  <tr><td>GSE230809</td><td>2023</td><td>NP, AF</td><td>24</td><td>105,804</td><td>10x</td><td>Healthy, Degenerated</td></tr>
  <tr><td>GSE242443</td><td>2024</td><td>CEP</td><td>2</td><td>59,227</td><td>10x</td><td>Healthy, Degenerated (culture-expanded)</td></tr>
</table>
""")

    # ── 3. Integration ──
    html.append("""
<h2 id="integration">3. Integration Strategy</h2>
<p><strong>Tiered approach:</strong> Non-resident cells (immune, endothelial, 14,566 cells) integrated with standard scVI. Resident cells integrated with 4 approaches (scVI, scANVI, Harmony, BBKNN) for NP (138,937 cells) and AF (282,736 cells).</p>
<p><strong>Primary: scANVI</strong> — best overall score (0.615) and cell type separation (ASW 0.511-0.521). Semi-supervised approach leverages marker-based annotations.</p>
<p><strong>Sensitivity: scVI</strong> — perfectly preserves cell state continuum (score variance ratio = 1.0) for trajectory analysis.</p>
""")

    # Integration images
    for img_name in ['umap_tier2_NP_by_approach.png', 'umap_tier2_AF_by_approach.png',
                     'metrics_comparison_AF.png']:
        img_path = RESULTS / 'integration' / img_name
        if img_path.exists():
            html.append(f'<img src="integration/{img_name}" alt="{img_name}">')

    # ── 4. DE ──
    html.append("""
<h2 id="de">4. Differential Expression</h2>
<p>Pseudobulk DE analysis using pyDESeq2. 17 powered comparisons, 128 skipped (underpowered). |log2FC| > 0.5, padj < 0.05.</p>
""")

    if not de_summary.empty:
        html.append('<table><tr><th>Cell Type</th><th>Comparison</th><th>Up</th><th>Down</th><th>Total</th></tr>')
        for _, row in de_summary.iterrows():
            html.append(f'<tr><td>{row["cell_type"]}</td><td>{row["comparison"]}</td>'
                       f'<td>{row["n_up"]}</td><td>{row["n_down"]}</td><td>{row["n_total"]}</td></tr>')
        html.append('</table>')

    html.append("""
<div class="finding">
<strong>Top DE genes in NP severe degeneration:</strong> CXCL1 (+3.75), CXCL3 (+3.72), CXCL2 (+3.13), TNF (+2.45), MDK (+2.72) — classical inflammatory/catabolic IVD signature.
</div>
<div class="finding">
<strong>Top DE genes in AF degeneration:</strong> CEMIP (+2.39, hyaluronidase), KRT16 (+2.84), CXCL8 (-2.19) — ECM degradation and stress markers.
</div>
""")

    # Volcano plots
    for vp in ['volcano_NP_mature_chondrocyte_mild_vs_severe.png',
               'volcano_AF_outer_healthy_vs_degenerated_severe.png']:
        img_path = RESULTS / 'differential' / 'volcano_plots' / vp
        if img_path.exists():
            html.append(f'<img src="differential/volcano_plots/{vp}">')

    html.append("""
<div class="caveat">
<strong>Caution:</strong> NP_mature_chondrocyte healthy_vs_herniated (4,316 DE genes) is flagged as likely study-confounded. Top genes include ribosomal proteins (RPL17, RPL36A), a batch artifact signature. Excluded from primary interpretation.
</div>
""")

    # ── 5. Pathways ──
    html.append("""
<h2 id="pathways">5. Biological Pathways</h2>
<p>ORA: 1,244 significantly enriched terms (FDR < 0.05). GSEA: 1,081 significant terms across GO, KEGG, Reactome, MSigDB Hallmark, and IVD-custom gene sets.</p>
""")

    for img in ['enrichment_AF_outer_up.png', 'enrichment_NP_mature_chondrocyte_up.png',
                'gsea_ivd_custom_heatmap.png']:
        img_path = RESULTS / 'interpretation' / 'pathway_enrichment' / img
        if img_path.exists():
            html.append(f'<img src="interpretation/pathway_enrichment/{img}">')

    # ── 6. TF ──
    html.append("""
<h2 id="tf">6. Transcription Factor Activity</h2>
<p>TF activity inferred using CollecTRI regulon overlap with Fisher's exact test. 113 significant TF-condition associations.</p>
""")

    tf_img = RESULTS / 'interpretation' / 'tf_activity' / 'tf_activity_heatmap.png'
    if tf_img.exists():
        html.append('<img src="interpretation/tf_activity/tf_activity_heatmap.png">')

    html.append("""
<div class="finding">
<strong>Key TFs:</strong> ATF3/ATF7 (stress response, NP severe), HSF1/HSF2 (heat shock, multiple cell types), NFKBIB (NF-kB pathway, NP stressed), E2F4/TFDP1 (cell cycle, NP severe degeneration).
</div>
""")

    # ── 7. Trajectory ──
    html.append("""
<h2 id="trajectory">7. Cell State Trajectories</h2>
<p>PAGA + diffusion pseudotime (DPT) on scANVI embeddings. NP: rooted at notochordal cells. AF: rooted at AF_inner.</p>
""")

    for img in ['umap_trajectory_NP.png', 'pseudotime_by_condition_NP.png',
                'gene_dynamics_NP.png', 'umap_trajectory_AF.png']:
        img_path = RESULTS / 'trajectories' / img
        if img_path.exists():
            html.append(f'<img src="trajectories/{img}">')

    html.append("""
<div class="finding">
<strong>Pseudotime correlates with disease:</strong> NP rho = -0.207, AF rho = -0.177 (both p < 10<sup>-100</sup>). Healthy cells at earlier pseudotime, degenerated at later. Sensitivity check with scVI embedding confirms direction (NP rho = -0.132).
</div>
<p>500 trajectory-associated genes per compartment. ~55% overlap with DE genes confirms trajectory captures disease-relevant biology, not batch effects.</p>
""")

    # ── 8. Communication ──
    html.append("""
<h2 id="communication">8. Cell-Cell Communication</h2>
<p>LIANA consensus (CellPhoneDB, NATMI, Connectome, SingleCellSignalR, log2FC) on 20,000 cells per condition.</p>
""")

    for img in ['interaction_heatmap_healthy.png', 'interaction_heatmap_degenerated.png',
                'differential_interactions.png']:
        img_path = RESULTS / 'communication' / 'interaction_plots' / img
        if img_path.exists():
            html.append(f'<img src="communication/interaction_plots/{img}">')

    html.append("""
<div class="finding">
<strong>Increased signaling in degeneration:</strong> 53,036 interactions (degenerated) vs 44,079 (healthy). Consistent with increased paracrine signaling and immune cell infiltration in degenerative discs.
</div>
""")

    # ── 9. Pain ──
    html.append("""
<h2 id="pain">9. Pain Biology</h2>
<p>Cross-reference of DE genes with curated pain gene sets (nociception, neurotrophins, nerve guidance, inflammatory pain, neovascularization).</p>

<div class="highlight">
<strong>Pain-associated findings:</strong>
<ul>
  <li><strong>TNF</strong> significantly upregulated in NP_stressed_degenerative (log2FC=+2.65) and NP_mature_chondrocyte (log2FC=+2.45) in severe degeneration. TNF is a key inflammatory pain mediator that sensitizes nerve endings.</li>
  <li><strong>CXCL8</strong> significantly downregulated in AF_outer with degeneration (log2FC=-2.19). May reflect altered chemokine balance.</li>
  <li><strong>Disc cells produce inflammatory mediators (TNF, CXCL1-3) but not nociceptors.</strong> This is consistent with the model that degenerated disc cells create a pro-inflammatory environment that promotes nerve ingrowth and sensitization, rather than directly signaling pain.</li>
  <li>3,662-4,194 pain-relevant ligand-receptor interactions identified through cell-cell communication analysis, including neurotrophin and VEGF signaling pathways.</li>
</ul>
</div>
""")

    pain_img = RESULTS / 'interpretation' / 'pain_genes_heatmap.png'
    if pain_img.exists():
        html.append('<img src="interpretation/pain_genes_heatmap.png">')

    # ── 10. Limitations ──
    html.append("""
<h2 id="limitations">10. Limitations</h2>
<div class="caveat">
<ul>
  <li><strong>Cross-study confounding:</strong> Condition and study are partially confounded, especially for herniated samples (only 2 studies). Within-study comparisons where possible.</li>
  <li><strong>Underpowered comparisons:</strong> 128/145 cell type x comparison combinations skipped due to insufficient samples. CEP compartment entirely underpowered for DE.</li>
  <li><strong>No RNA velocity:</strong> Spliced/unspliced counts not available in public datasets. Would require reprocessing from BAM files.</li>
  <li><strong>Age-disease confound:</strong> In GSE230809, healthy donors are 21-27y and diseased are 37-73y. Cannot fully separate age from disease effects.</li>
  <li><strong>Sex bias:</strong> GSE230809 (largest dataset, 24 samples) is all-male. 30/78 samples have unknown sex.</li>
  <li><strong>Culture-expanded cells:</strong> GSE242443 CEP cells are culture-expanded, which alters gene expression.</li>
  <li><strong>Endothelial annotation caveat:</strong> Some endothelial DE genes (ACAN, IBSP) suggest possible misclassification of NP/AF cells.</li>
  <li><strong>Composition analysis:</strong> No significant changes after FDR correction, though trends are biologically consistent.</li>
  <li><strong>SCENIC/GRN not run:</strong> Full SCENIC analysis was not performed due to computational requirements. TF activity estimated from CollecTRI regulon overlap instead.</li>
</ul>
</div>
""")

    # ── 11. Methods ──
    html.append("""
<h2 id="methods">11. Methods</h2>
<div class="methods">
<h3>Data acquisition</h3>
<p>12 scRNA-seq datasets of human IVD tissue were downloaded from GEO and CNGB (see Table 1). Raw count matrices were obtained for each dataset.</p>

<h3>Quality control and preprocessing</h3>
<p>Per-dataset QC: min 200 genes, max 6000 genes, min 500 counts, max 20% mitochondrial reads. Doublet detection with Scrublet (expected rate 5%). Normalization: total-count to 10,000, log1p. HVG selection: top 2000 genes per dataset using Seurat v3 method.</p>

<h3>Cell type annotation</h3>
<p>Two-pass annotation: (1) marker-based scoring using 16 IVD-specific gene signatures, (2) CellTypist Immune_All_Low model for immune subtypes. Consensus labels in cell_type_final.</p>

<h3>Integration</h3>
<p>Tiered strategy: non-resident cells (immune, endothelial) integrated with scVI (1 layer, 128 dim). Resident cells: 4 approaches tested (scVI, scANVI, Harmony, BBKNN). scANVI selected as primary (best overall scIB score 0.615, celltype ASW 0.511-0.521).</p>

<h3>Differential expression</h3>
<p>Pseudobulk aggregation per sample per cell type. DE with pyDESeq2 (Python DESeq2 implementation). Significance: |log2FC| > 0.5, adjusted p-value < 0.05 (Benjamini-Hochberg). Minimum 3 samples per condition per cell type.</p>

<h3>Pathway enrichment</h3>
<p>Over-representation analysis (ORA) and gene set enrichment analysis (GSEA) using gseapy. Databases: GO Biological Process 2023, KEGG 2021, Reactome 2022, MSigDB Hallmark 2020, custom IVD gene sets.</p>

<h3>TF activity inference</h3>
<p>CollecTRI regulon network (42,990 interactions, 1,185 TFs). TF activity scored by Fisher's exact test for enrichment of TF targets among DE genes, with concordance scoring for direction.</p>

<h3>Trajectory analysis</h3>
<p>PAGA + diffusion pseudotime (DPT) on scANVI embeddings. 50,000 cells per compartment. Root cells: NP notochordal, AF inner. Trajectory genes: Spearman correlation with pseudotime, FDR < 0.05, top 500.</p>

<h3>Cell-cell communication</h3>
<p>LIANA rank_aggregate with consensus resource. 5 methods: CellPhoneDB, NATMI, Connectome, SingleCellSignalR, log2FC. 100 permutations. 20,000 cells per condition.</p>

<h3>Software</h3>
<p>Python 3.12, scanpy 1.11, scvi-tools 1.4.2, pyDESeq2, gseapy 1.1, decoupler 2.1, liana 1.7, harmonypy, bbknn. Full environment: requirements.txt.</p>
</div>
""")

    # ── Reproducibility ──
    html.append("""
<h2>12. Reproducibility</h2>
<ul>
  <li>All scripts version-controlled in git</li>
  <li>Random seeds: 42 (all stochastic operations)</li>
  <li>Package versions pinned in requirements.txt</li>
  <li>All parameter choices documented in analysis_plan.md</li>
  <li>All human checkpoint decisions recorded</li>
  <li>Data provenance: GEO/CNGB accessions, download dates in dataset_registry.tsv</li>
</ul>
""")

    html.append('</body></html>')

    report_path = RESULTS / "final_report.html"
    report_path.write_text('\n'.join(html))
    print(f"  Final report: {report_path}")


def generate_requirements():
    """Capture current package versions."""
    import subprocess
    result = subprocess.run(['pip', 'freeze'], capture_output=True, text=True)
    (BASE / 'requirements_frozen.txt').write_text(result.stdout)
    print("  Requirements frozen: requirements_frozen.txt")


def validate():
    """Validation checks."""
    print("\n" + "=" * 60)
    print("Validating Module 10 outputs")
    print("=" * 60)

    checks = []

    if (RESULTS / "final_report.html").exists():
        checks.append(('PASS', 'Final report generated'))
    else:
        checks.append(('FAIL', 'Final report not found'))

    supp_count = len(list(SUPP_TABLES.glob("S*.tsv")))
    if supp_count > 0:
        checks.append(('PASS', f'{supp_count} supplementary tables collected'))
    else:
        checks.append(('WARN', 'No supplementary tables'))

    if (BASE / 'requirements_frozen.txt').exists():
        checks.append(('PASS', 'Package versions frozen'))

    # Check all module scripts exist
    for i in range(1, 11):
        scripts = list(BASE.glob(f"scripts/{i:02d}_*.py"))
        if scripts:
            checks.append(('PASS', f'Module {i:02d} script exists'))
        elif i <= 5:
            checks.append(('WARN', f'Module {i:02d} script not found (may be pre-pipeline)'))

    # Check all reports exist
    for report_name, report_path in [
        ('Integration', RESULTS / 'integration' / 'integration_report.html'),
        ('Differential', RESULTS / 'differential' / 'differential_report.html'),
        ('Interpretation', RESULTS / 'interpretation' / 'interpretation_report.html'),
        ('Trajectory', RESULTS / 'trajectories' / 'trajectory_report.html'),
        ('Communication', RESULTS / 'communication' / 'communication_report.html'),
    ]:
        if report_path.exists():
            checks.append(('PASS', f'{report_name} report exists'))
        else:
            checks.append(('WARN', f'{report_name} report not found'))

    passed = all(s != 'FAIL' for s, _ in checks)
    for status, msg in checks:
        print(f"  [{status}] {msg}")
    return passed


def main():
    print("=" * 60)
    print("Module 10: Final Reporting")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    if '--validate-only' in sys.argv:
        passed = validate()
        sys.exit(0 if passed else 1)

    print("\nCollecting supplementary tables...")
    collect_supplementary_tables()

    print("\nGenerating final report...")
    generate_final_report()

    print("\nFreezing requirements...")
    generate_requirements()

    validate()

    print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("Pipeline complete. Final report at: results/final_report.html")


if __name__ == "__main__":
    main()
