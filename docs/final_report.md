# Human Intervertebral Disc Single-Cell Atlas — Final Report

**A comprehensive scRNA-seq meta-analysis of IVD degeneration**

| Field | Value |
|-------|-------|
| Report generated | 2026-03-25 14:23 |
| Pipeline version | v5 |
| Git commit | `eabcb1c` (branch: `main`) |
| Source of truth | `analysis_plan.md` |

## Contents

1. [Atlas Summary](#summary)
2. [Dataset Summary](#datasets)
3. [Integration](#integration)
4. [Clustering & Annotation](#clustering)
5. [Differential Expression](#de)
6. [Biological Pathways](#pathways)
7. [Transcription Factor Activity](#tf)
8. [Cell State Trajectories](#trajectory)
9. [Cell-Cell Communication](#communication)
10. [Pain Biology](#pain)
11. [Limitations & Caveats](#limitations)
12. [Methods](#methods)
13. [Reproducibility](#reproducibility)

## 1. Atlas Summary {#summary}

| Metric | Value | Source |
|--------|-------|--------|
| Samples | 78 | *[source: `metadata/sample_metadata.tsv`]* |
| Donors | 57 | *[source: `metadata/sample_metadata.tsv`]* |
| Compartments | AF, CEP, IVD_mixed, NP | *[source: `metadata/sample_metadata.tsv`]* |
| Powered DE comparisons | 10 | *[source: `docs/v5_results/de_summary_table.tsv`]* |
| Significant DE genes (total hits) | 1,198 | *[source: `docs/v5_results/de_summary_table.tsv`]* |
| L-R interactions (all conditions) | 59,745 | *[source: `docs/v5_results/communication_stats.tsv`]* |

## 2. Dataset Summary {#datasets}

**12 datasets included** (of 22 evaluated). *[source: `metadata/dataset_registry.tsv`]*

| accession | first_author | year | compartment | n_samples | technology | conditions |
| --- | --- | --- | --- | --- | --- | --- |
| GSE160756 | Gan Y | 2021 | NP, AF, CEP | 7.0 | 10x Genomics | Healthy young and adult |
| GSE165722 | Tu J | 2022 | NP | 8.0 | BD Rhapsody | Pfirrmann Grade II-V progressive |
| GSE189916 | Jiang W | 2022 | Whole IVD | 6.0 | 10x Genomics Chromium 3' v2 | Neonatal vs adult |
| GSE199866 | Cherif H | 2022 | NP, inner AF | 4.0 | 10x Genomics | Paired degen vs non-degen from same individual |
| GSE205535 | Li Z | 2022 | NP | 2.0 | BD Rhapsody | Normal vs degenerative NP |
| GSE233666 | Guo S | 2023 | NP | 4.0 | 10x Genomics | IDD diagnosis - disc herniation |
| GSE244889 | Chen F | 2024 | NP | 7.0 | 10x Genomics | Mild vs severe degeneration (Pfirrmann I-II vs III-V) |
| GSE251686 | Jia S | 2024 | NP | 6.0 | 10x Genomics | Mild vs severe degeneration |
| GSE255768 | Shi C | 2024 | CEP/Endplate | 2.0 | 10x Genomics | Degenerative endplate (Modic changes) |
| GSE230809 | Swahn H | 2024 | NP, AF (surgically separated) | 24.0 | 10x Genomics | Healthy (Thompson II) vs diseased (Thompson II-IV) |
| CNP0002664 | Han S | 2022 | NP | 6.0 | Singleron Matrix | Normal, mild, severe IVDD |
| GSE242443 | Kuchynsky K | 2024 | CEP | 2.0 | 10x Genomics 3' v3.1 | Non-degen vs degen CEP |

**Platform heterogeneity:** 2 non-10x platform(s) (BD Rhapsody, Singleron Matrix). Batch correction handles platform differences via study-level integration keys.

### Sample demographics

- Age range: 0–81 years (21 samples with unknown age). *[source: `metadata/sample_metadata.tsv`]*
- Sex distribution: M=36, unknown=30, F=12. *[source: `metadata/sample_metadata.tsv`]*

## 3. Integration {#integration}

### Workflow comparison

Three integration workflows were compared. *[source: `analysis_plan.md`]*

| Object | Workflow | Cells | Clusters | iLISI | batch_ASW | condition_ASW |
|--------|----------|-------|----------|-------|-----------|---------------|
| NP | **CCA** | 262,967 | 24 | **3.68** | -0.11 | -0.16 |
| NP | scANVI | 262,967 | 29 | 1.23 | 0.08 | 0.00 |
| NP | STACAS | 16,000* | 21 | 2.08 | -0.06 | -0.05 |
| AF | **CCA** | 84,624 | 22 | **1.49** | -0.12 | 0.05 |
| AF | scANVI | 84,568 | 18 | 1.01 | 0.16 | 0.02 |
| AF | STACAS | 84,624 | 23 | 1.06 | 0.05 | 0.01 |
| CEP | **CCA** | 50,858 | 14 | **1.63** | -0.07 | -0.09 |
| CEP | scANVI | 50,769 | 13 | 1.03 | 0.21 | 0.04 |
| CEP | STACAS | 50,858 | 15 | 1.13 | 0.05 | 0.00 |
| all | **CCA** | 410,759 | 44 | **3.18** | -0.15 | -0.14 |
| all | scANVI | 410,759 | 29 | 1.23 | 0.07 | -0.02 |
| all | STACAS | 30,000* | 17 | 2.42 | -0.06 | -0.10 |

**Rationale for CCA:**
- Label-free: no circular annotation risk (does not depend on Module 04 coarse labels)
- Full cell counts for all 4 objects (no downsampling)
- Strongest batch mixing (iLISI 1.5-3.7 vs ~1.0-1.2 for scANVI)
- Smooth embedding topology consistent with mesenchymal continuum hypothesis
- Negative batch_ASW indicates possible overcorrection, but DE uses pseudobulk on raw counts (not embeddings)


## 4. Clustering & Annotation {#clustering}

### Cell type definitions

*[source: `docs/v5_results/cell_type_definitions.tsv`]*

| object | cell_type | coarse_cell_type | n_cells | clusters | canonical_markers | confidence_distribution |
| --- | --- | --- | --- | --- | --- | --- |
| NP | Endothelial | Endothelial | 2645 | NM0,NM1,NM4 | PECAM1,VWF,CDH5 | high=1381; medium=1264 |
| NP | Macrophage_M2 | Macrophage | 181 | NM3 | CD163,MRC1,MSR1,TGFB1 | high=181; medium=0 |
| NP | NP_fibrocartilaginous | Fibroblast_like | 73764 | M2,M3,M4,M5 | COL1A1,COL2A1,VCAN | high=73764; medium=0 |
| NP | NP_mature_chondrocyte | Chondrocyte_like | 185794 | M0,M1,M6 | ACAN,COL2A1,SOX9,COMP,PRG4 | high=185794; medium=0 |
| NP | T_cell_CD8 | T_cell | 583 | NM2 | CD8A,CD8B,GZMB,PRF1 | medium=583; high=0 |
| AF | AF_inner | Chondrocyte_like | 32839 | M0,M8,M9 | COL2A1,ACAN,SOX9 | high=32839 |
| AF | AF_outer | Fibroblast_like | 51729 | M1,M2,M3,M4,M5,M6,M7 | COL1A1,COL1A2,THY1,DCN,LUM | high=51729 |
| AF | Endothelial | Endothelial | 22 | NM1 | PECAM1,VWF,CDH5 | high=22 |
| AF | Macrophage_M2 | Macrophage | 34 | NM0 | CD163,MRC1,MSR1,TGFB1 | high=34 |
| CEP | EP_hyaline | Chondrocyte_like | 12597 | M1 | COL2A1,COL10A1,SOX9 | high=12597; medium=0 |
| CEP | Endothelial | Endothelial | 38 | NM0 | PECAM1,VWF,CDH5 | high=38; medium=0 |
| CEP | Fibroblast_like | Fibroblast_like | 33582 | M0,M3,M4 | COL1A1,COL1A2,DCN,LUM,THY1 | high=33582; medium=0 |
| CEP | Fibrochondrocyte_chondroid | Fibrochondrocyte_like | 4292 | M2 | COL2A1,ACAN,SOX9 | high=4292; medium=0 |
| CEP | Fibrochondrocyte_fibroid | Fibrochondrocyte_like | 298 | M5 | COL1A1,COL1A2,DCN | medium=298; high=0 |
| CEP | NK_cell | NK_cell | 21 | NM2 | NKG7,GNLY | high=21; medium=0 |
| CEP | Pericyte_SMC | Pericyte_SMC | 30 | NM1 | ACTA2,RGS5,PDGFRB | high=30; medium=0 |
| all_cells | AF_inner | Chondrocyte_like | 32839 | M0,M1,M2,M3,M4,M5,M6,M7 | COL2A1,ACAN,SOX9 | high=32839; low=0; medium=0 |
| all_cells | AF_outer | Fibroblast_like | 51729 | M0,M1,M2,M3,M4,M5,M6,M7 | COL1A1,COL1A2,THY1,DCN,LUM | high=51729; low=0; medium=0 |
| all_cells | Chondrocyte_like | Chondrocyte_like | 8794 | M0,M1,M2,M4,M8 | COL2A1,ACAN,SOX9,COMP,PRG4 | high=6017; medium=2777; low=0 |
| all_cells | EP_hyaline | Chondrocyte_like | 12597 | M0,M1,M2,M3,M4,M5,M6,M7 | COL2A1,COL10A1,SOX9 | high=12597; low=0; medium=0 |
| all_cells | Endothelial | Endothelial | 2715 | NM0,NM1,NM2,NM3,NM4 | PECAM1,VWF,CDH5 | high=1451; medium=1264; low=0 |
| all_cells | Fibroblast_like | Fibroblast_like | 34020 | M0,M1,M2,M3,M4,M5,M6,M7 | COL1A1,COL1A2,DCN,LUM,THY1 | high=33582; medium=438; low=0 |
| all_cells | Fibrochondrocyte_chondroid | Fibrochondrocyte_like | 4292 | M0,M1,M2,M3,M4,M5,M6,M7 | COL2A1,ACAN,SOX9 | high=4292; low=0; medium=0 |
| all_cells | Fibrochondrocyte_fibroid | Fibrochondrocyte_like | 298 | M0,M1,M2,M3,M4,M5,M6 | COL1A1,COL1A2,DCN | medium=298; high=0; low=0 |
| all_cells | Fibrochondrocyte_like | Fibrochondrocyte_like | 3052 | M3,M5,M7 | COL1A1,COL2A1,ACAN,DCN | low=2304; medium=748; high=0 |
| all_cells | Macrophage | Macrophage | 16 | NM0,NM2 | CD68,CD14,CSF1R,CD163,CD86 | high=16; low=0; medium=0 |
| all_cells | Macrophage_M2 | Macrophage | 215 | NM0,NM1,NM2,NM3,NM4 | CD163,MRC1,MSR1,TGFB1 | high=215; low=0; medium=0 |
| all_cells | NK_cell | NK_cell | 21 | NM0,NM1,NM2 | NKG7,GNLY | high=21; low=0; medium=0 |
| all_cells | NP_fibrocartilaginous | Fibroblast_like | 73764 | M0,M1,M2,M3,M4,M5,M6,M7,M9 | COL1A1,COL2A1,VCAN | high=73764; low=0; medium=0 |
| all_cells | NP_mature_chondrocyte | Chondrocyte_like | 185794 | M0,M1,M2,M3,M4,M5,M6,M7,M8 | ACAN,COL2A1,SOX9,COMP,PRG4 | high=185794; low=0; medium=0 |
| all_cells | Pericyte_SMC | Pericyte_SMC | 30 | NM0,NM1,NM4 | ACTA2,RGS5,PDGFRB | high=30; low=0; medium=0 |
| all_cells | T_cell_CD8 | T_cell | 583 | NM0,NM1,NM2,NM3 | CD8A,CD8B,GZMB,PRF1 | medium=583; high=0; low=0 |

### Clustering resolution optimization

*[source: `results/integration/clustering_resolution_optimization`]*

- **TEST_mes:** best resolution 0.4, 2 clusters (silhouette=0.102)

## 5. Differential Expression {#de}

**10 powered comparisons**, 1,198 significant genes (460 up, 738 down). Thresholds: |log2FC| > 0.5, padj < 0.05. *[source: `docs/v5_results/de_summary_table.tsv`]*

| cell_type | comparison | n_up | n_down | n_total |
| --- | --- | --- | --- | --- |
| AF_inner | healthy_vs_degenerated_severe | 3 | 2 | 5 |
| AF_outer | healthy_vs_degenerated_all | 0 | 1 | 1 |
| AF_outer | healthy_vs_degenerated_mild | 4 | 114 | 118 |
| AF_outer | healthy_vs_degenerated_severe | 15 | 4 | 19 |
| AF_outer | mild_vs_severe | 3 | 1 | 4 |
| NP_fibrocartilaginous | healthy_vs_degenerated_all | 20 | 8 | 28 |
| NP_fibrocartilaginous | healthy_vs_degenerated_mild | 7 | 7 | 14 |
| NP_fibrocartilaginous | healthy_vs_degenerated_severe | 237 | 319 | 556 |
| NP_fibrocartilaginous | mild_vs_severe | 138 | 263 | 401 |
| NP_mature_chondrocyte | mild_vs_severe | 33 | 19 | 52 |

### Skipped comparisons (underpowered)

47 comparisons skipped due to insufficient samples (< 3 per condition per cell type). *[source: `docs/v5_results/skipped_comparisons.tsv`]*

## 6. Biological Pathways {#pathways}

ORA results not found. *[source: `results/interpretation/pathway_enrichment/all_enrichment_results.tsv`]*

GSEA results not found. *[source: `results/interpretation/pathway_enrichment/gsea_results.tsv`]*

## 7. Transcription Factor Activity {#tf}

**20 significant TF-condition associations** (padj < 0.05). *[source: `docs/v5_results/tf_activity_top.tsv`]*

**20 unique TFs** with significant activity changes.

| TF | Score | Cell Type | Comparison |
|----|-------|-----------|------------|
| PITX1 | -0.062 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| ID1 | 0.041 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| SMAD7 | -0.016 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| SMAD1 | 0.015 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| NFATC2 | 0.013 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| BCL6 | -0.011 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| JUND | 0.011 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| RUNX2 | -0.009 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| KLF4 | 0.008 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| ATF4 | 0.008 | NP_fibrocartilaginous | healthy_vs_degenerated_all |

## 8. Cell State Trajectories {#trajectory}

PAGA + diffusion pseudotime (DPT) analysis. *[source: `docs/v5_results/pseudotime_correlations.tsv`]*

| compartment | embedding | test | rho | pvalue | n | statistic | n_healthy | n_degen | median_healthy | median_degen |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NP | scVI_mesenchymal | pseudotime_vs_condition_ordinal | -0.0882 | 0.0000 | 5.00e+04 | nan | nan | nan | nan | nan |
| NP | scVI_mesenchymal | pseudotime_healthy_vs_degenerated | nan | 0.0000 | nan | 2.66e+08 | 1.46e+04 | 2.99e+04 | 0.0805 | 0.0503 |
| NP | scVI_mesenchymal | pseudotime_vs_condition_NP_mature_chondrocyte | -0.0018 | 0.7709 | 2.70e+04 | nan | nan | nan | nan | nan |
| NP | scVI_mesenchymal | pseudotime_vs_condition_NP_fibrocartilaginous | -0.2017 | 0.0000 | 2.30e+04 | nan | nan | nan | nan | nan |
| AF | scVI_mesenchymal | pseudotime_vs_condition_ordinal | 0.1947 | 0.0000 | 5.00e+04 | nan | nan | nan | nan | nan |
| AF | scVI_mesenchymal | pseudotime_healthy_vs_degenerated | nan | 0.0000 | nan | 2.77e+08 | 2.48e+04 | 2.52e+04 | 0.6587 | 0.6587 |
| AF | scVI_mesenchymal | pseudotime_vs_condition_AF_outer | 0.1096 | 0.0000 | 3.78e+04 | nan | nan | nan | nan | nan |
| AF | scVI_mesenchymal | pseudotime_vs_condition_AF_inner | 0.2178 | 0.0000 | 1.22e+04 | nan | nan | nan | nan | nan |
| CEP | scVI_mesenchymal | pseudotime_vs_condition_ordinal | 0.0734 | 0.0000 | 3.21e+04 | nan | nan | nan | nan | nan |
| CEP | scVI_mesenchymal | pseudotime_healthy_vs_degenerated | nan | 0.0000 | nan | 1.08e+08 | 2.05e+04 | 1.17e+04 | 0.0653 | 0.0662 |
| CEP | scVI_mesenchymal | pseudotime_vs_condition_Fibrochondrocyte_chondroid | -0.2494 | 0.0000 | 2.21e+03 | nan | nan | nan | nan | nan |
| CEP | scVI_mesenchymal | pseudotime_vs_condition_EP_hyaline | 0.1373 | 0.0000 | 7.42e+03 | nan | nan | nan | nan | nan |
| CEP | scVI_mesenchymal | pseudotime_vs_condition_Fibroblast_like | -0.3059 | 0.0000 | 2.24e+04 | nan | nan | nan | nan | nan |
| CEP | scVI_mesenchymal | pseudotime_vs_condition_Fibrochondrocyte_fibroid | -0.4296 | 0.0000 | 1.62e+02 | nan | nan | nan | nan | nan |

### Trajectory-associated genes

*[source: `docs/v5_results/trajectory_gene_counts.tsv`]*

- **AF:** 500 genes correlated with pseudotime (FDR < 0.05)
- **CEP:** 500 genes correlated with pseudotime (FDR < 0.05)
- **NP:** 500 genes correlated with pseudotime (FDR < 0.05)

## 9. Cell-Cell Communication {#communication}

LIANA rank_aggregate (CellPhoneDB, NATMI, Connectome, SingleCellSignalR, log2FC). *[source: `docs/v5_results/communication_stats.tsv`]*

| Condition | Interactions |
|-----------|-------------|
| degenerated | 34,208 |
| healthy | 25,537 |

## 10. Pain Biology {#pain}

**13 pain-associated DE genes identified.** *[source: `docs/v5_results/pain_genes.tsv`]*

| gene | padj | cell_type | comparison |
| --- | --- | --- | --- |
| IL1B | 7.83e-03 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| NTN1 | 0.031 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| PLA2G2A | 7.83e-03 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| NTN1 | 1.17e-04 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| NTN4 | 0.045 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| PDGFA | 0.030 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| PENK | 0.044 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| PLA2G2A | 6.65e-03 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| UNC5B | 6.65e-03 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| VEGFA | 0.017 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| IL6 | 5.22e-03 | NP_fibrocartilaginous | mild_vs_severe |
| CXCL8 | 7.58e-03 | AF_outer | healthy_vs_degenerated_mild |
| PLA2G2A | 8.62e-03 | AF_outer | healthy_vs_degenerated_mild |

## 11. Limitations & Caveats {#limitations}

### Data and design limitations

*[source: `analysis_plan.md`]*

- **NGDC datasets excluded:** PRJCA014236 and PRJCA007656 not downloaded. NP already well-covered.
- **GSE205535 corrigenda:** Published corrections exist — reviewed during preprocessing.
- **Platform heterogeneity:** 3 non-10x datasets (BD Rhapsody, Singleron). Handled by CCA batch correction (v5 primary). scANVI and STACAS also tested.
- **CEP underpowered:** Only 3 CEP datasets (6 samples). Compartment-specific CEP analyses are limited.
- **GSE242443 culture-expanded:** CEP cells are culture-expanded. Included with caveats.
- **GSE230809 sex bias:** All 24 samples from male donors. Limits sex-stratified analyses.
- **GSE230809 age-disease confound:** Healthy=21-27y, diseased=37-73y. Cannot separate age from disease.
- **Missing demographics:** 18/78 samples unknown age, 30/78 unknown sex.
- **GSE251686_NP3 excluded:** Corrupt matrix file (5/6 samples retained).
- **GSE165722 Pfirrmann offset:** GEO says I-IV, paper says II-V. Paper grades used.

### Underpowered comparisons

47 cell type x condition comparisons were skipped due to insufficient sample counts. *[source: `docs/v5_results/skipped_comparisons.tsv`]*

### Result sensitivity across pipeline versions

Several results are sensitive to upstream methodological choices (integration method, annotation, cell sampling). These are documented here to flag areas requiring cautious interpretation. *[source: `docs/version_history.md`]*

- **Trajectory pseudotime-condition correlations** are sensitive to integration method and root cell choice. In v5 (CCA): 
  NP rho=-0.088; AF rho=+0.195; CEP rho=+0.073. Prior versions showed sign changes (e.g., CEP: -0.163 in v2, +0.135 in v3, +0.073 in v5), indicating these correlations are not robust to upstream choices.

- **CCC interaction counts** in v5: degenerated: 34,208, healthy: 25,537. The direction of the healthy-vs-degenerated difference has varied across pipeline versions (v1: degenerated > healthy; v2: healthy > degenerated; v3: near-equal), making this result sensitive to cell type definitions and sampling.

- **CellTypist concordance** is limited for IVD-specific cell types. CellTypist lacks IVD reference data, so disagreements with de novo labels are expected for mesenchymal populations. De novo labels are retained as primary annotations; CellTypist is used for immune subtype validation only.

- **AF pseudotime sign:** AF consistently shows positive rho (degenerated cells at later pseudotime) across pipeline versions, opposite to NP. This may reflect genuine compartment-specific biology or root cell choice effects.

## 12. Methods {#methods}

Full parameter choices and rationale documented in `analysis_plan.md`. *[source: `analysis_plan.md`]*

### Data acquisition

12 scRNA-seq datasets of human IVD tissue downloaded from GEO and CNGB. Raw count matrices obtained per dataset. See `metadata/dataset_registry.tsv` for accessions and details. *[source: `scripts/01_dataset_download.py`, `metadata/dataset_registry.tsv`]*

### Quality control and preprocessing

Per-dataset QC: min 200 genes, max 6000 genes, min 500 counts, max 20% mitochondrial reads. Doublet detection with Scrublet (expected rate 5%). Normalization: total-count to 10,000, log1p. HVG selection: top 2000 genes per dataset (Seurat v3 method). *[source: `scripts/03_preprocessing.py`, `specs/03_PREPROCESSING.md`]*

### Cell classification

Coarse classification into 5 anchor categories using marker gene scoring (immune: PTPRC, CD3D, CD68, PECAM1; mesenchymal: COL2A1, COL1A1, ACAN, SOX9). Cluster-level majority voting with 85% threshold. *[source: `scripts/04_annotation.py`, `specs/04_ANNOTATION.md`]*

### Integration

Three workflows compared (CCA, scANVI, STACAS) on four compartment objects (NP, AF, CEP, all_cells). CCA (Seurat v5 `IntegrateLayers(method=CCAIntegration)`) selected as primary: label-free, full cell counts, strongest batch mixing (iLISI). *[source: `scripts/05a_integration_cca.R`, `specs/05_INTEGRATION.md`, `analysis_plan.md`]*

### Clustering

Leiden clustering with multi-resolution optimization. Resolution selected by silhouette score. *[source: `scripts/06_clustering.py`, `specs/06_CLUSTERING.md`]*

### Post-integration annotation

De novo cell type annotation from cluster DE markers and canonical marker panels. CellTypist (Immune_All_Low model) for immune subtype validation. *[source: `scripts/07_annotation.py`, `specs/07_ANNOTATION.md`]*

### Differential expression

Pseudobulk aggregation per sample per cell type. DE with pyDESeq2. Significance: |log2FC| > 0.5, padj < 0.05 (Benjamini-Hochberg). Minimum 3 samples per condition per cell type. *[source: `scripts/08_differential.py`, `specs/08_DIFFERENTIAL.md`]*

### Pathway enrichment

ORA and GSEA using gseapy. Databases: GO Biological Process 2023, KEGG 2021, Reactome 2022, MSigDB Hallmark 2020, custom IVD gene sets. *[source: `scripts/09_interpretation.py`, `specs/09_INTERPRETATION.md`]*

### TF activity inference

CollecTRI regulon network. TF activity scored by Fisher's exact test for enrichment of TF targets among DE genes. *[source: `scripts/09_interpretation.py`]*

### Trajectory analysis

PAGA + diffusion pseudotime (DPT) on mesenchymal embeddings. Root cells defined per compartment (NP: notochordal; AF: AF_inner). Trajectory genes: Spearman correlation with pseudotime, FDR < 0.05. *[source: `scripts/10_trajectory.py`, `specs/10_TRAJECTORY.md`]*

### Cell-cell communication

LIANA rank_aggregate with consensus resource. 5 methods: CellPhoneDB, NATMI, Connectome, SingleCellSignalR, log2FC. 100 permutations. *[source: `scripts/11_communication.py`, `specs/11_COMMUNICATION.md`]*

### Software

Python 3.12, scanpy, scvi-tools, pyDESeq2, gseapy, decoupler, liana. R: Seurat 5.4.0, STACAS 2.4.1. Full environment: `requirements.txt` / `requirements_frozen.txt`.

## 13. Reproducibility {#reproducibility}

- **Git commit:** `eabcb1c2d89a65322284db0e70d8a1d423c15b12` (branch: `main`)
- **Random seeds:** 42 (all stochastic operations)
- **Package versions:** pinned in `requirements.txt`, frozen in `requirements_frozen.txt`
- **Parameter choices:** documented in `analysis_plan.md`
- **Human checkpoint decisions:** recorded in `analysis_plan.md`
- **Data provenance:** GEO/CNGB accessions and download dates in `metadata/dataset_registry.tsv`
- **File checksums:** `metadata/file_checksums.json`
