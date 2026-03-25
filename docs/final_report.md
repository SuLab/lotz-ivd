# Human Intervertebral Disc Single-Cell Atlas — Final Report

**A comprehensive scRNA-seq meta-analysis of IVD degeneration**

| Field | Value |
|-------|-------|
| Report generated | 2026-03-25 23:43 |
| Pipeline version | v5 |
| Git commit | `7338828` (branch: `main`) |
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
| Powered DE comparisons | 10 | *[source: `results/differential/de_summary_table.tsv`]* |
| Significant DE genes (total hits) | 1,198 | *[source: `results/differential/de_summary_table.tsv`]* |
| Enriched pathways (ORA) | 11,407 | *[source: `results/interpretation/pathway_enrichment/all_enrichment_results.tsv`]* |
| L-R interactions (all conditions) | 59,745 | *[source: `results/communication`]* |

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


### CellTypist validation

5 concordant, 5 discordant cluster(s) across all objects. *[source: `results/integration/celltypist_validation`]*

| Object | Cluster | Cells | De Novo Label | CellTypist Label | Agreement % |
|--------|---------|-------|---------------|------------------|-------------|
| AF | 0 | 34 | Macrophage_M2 | Fibroblasts | 38.2% |
| CEP | 1 | 30 | Pericyte_SMC | NK cells | 76.7% |
| CEP | 2 | 21 | NK_cell | Classical monocytes | 52.4% |
| NP | 1 | 1,264 | Endothelial | Classical monocytes | 67.9% |
| NP | 2 | 583 | T_cell_CD8 | Classical monocytes | 46.3% |

De novo labels retained as primary (CellTypist lacks IVD-specific cell types).

## 4. Clustering & Annotation {#clustering}

### Cell type definitions

*[source: `results/integration/cell_type_definitions.tsv`]*

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

- **AF_mesenchymal:** best resolution 0.2, 10 clusters (silhouette=0.114)
- **AF_non_mesenchymal:** best resolution 0.2, 2 clusters (silhouette=0.379)
- **CEP_mesenchymal:** best resolution 0.2, 6 clusters (silhouette=0.135)
- **CEP_non_mesenchymal:** best resolution 0.3, 3 clusters (silhouette=0.414)
- **NP_mesenchymal:** best resolution 0.2, 7 clusters (silhouette=0.105)
- **NP_non_mesenchymal:** best resolution 0.2, 3 clusters (silhouette=0.223)
- **all_cells_mesenchymal:** best resolution 0.4, 10 clusters (silhouette=0.036)
- **all_cells_non_mesenchymal:** best resolution 0.2, 2 clusters (silhouette=0.634)

## 5. Differential Expression {#de}

**10 powered comparisons**, 1,198 significant genes (460 up, 738 down). Thresholds: |log2FC| > 0.5, padj < 0.05. *[source: `results/differential/de_summary_table.tsv`]*

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

47 comparisons skipped due to insufficient samples (< 3 per condition per cell type). *[source: `results/differential/skipped_comparisons.tsv`]*

## 6. Biological Pathways {#pathways}

**ORA:** 2,506 significantly enriched terms (FDR < 0.05). *[source: `results/interpretation/pathway_enrichment/all_enrichment_results.tsv`]*

**GO_Biological_Process_2023** (top 5):

- Mitotic Sister Chromatid Segregation (GO:0000070) (padj=1.14e-38, NP_fibrocartilaginous down)
- DNA Metabolic Process (GO:0006259) (padj=2.17e-27, NP_fibrocartilaginous down)
- Mitotic Sister Chromatid Segregation (GO:0000070) (padj=1.76e-26, NP_fibrocartilaginous down)
- Microtubule Cytoskeleton Organization Involved In Mitosis (GO:1902850) (padj=3.87e-26, NP_fibrocartilaginous down)
- Mitotic Spindle Organization (GO:0007052) (padj=2.20e-25, NP_fibrocartilaginous down)

**IVD_custom** (top 5):

- ECM_homeostasis (padj=2.69e-06, AF_outer down)
- Inflammatory_signaling (padj=8.11e-06, NP_fibrocartilaginous down)
- Inflammatory_signaling (padj=2.04e-05, NP_fibrocartilaginous down)
- Inflammatory_pain (padj=1.48e-03, AF_outer down)
- Nerve_guidance (padj=1.59e-03, NP_fibrocartilaginous up)

**KEGG_2021_Human** (top 5):

- Cell cycle (padj=1.26e-29, NP_fibrocartilaginous down)
- Cell cycle (padj=4.26e-27, NP_fibrocartilaginous down)
- Fanconi anemia pathway (padj=2.28e-12, NP_fibrocartilaginous down)
- DNA replication (padj=6.27e-12, NP_fibrocartilaginous down)
- Fanconi anemia pathway (padj=5.63e-10, NP_fibrocartilaginous down)

**MSigDB_Hallmark_2020** (top 5):

- E2F Targets (padj=7.74e-85, NP_fibrocartilaginous down)
- E2F Targets (padj=7.93e-84, NP_fibrocartilaginous down)
- G2-M Checkpoint (padj=6.38e-77, NP_fibrocartilaginous down)
- G2-M Checkpoint (padj=1.07e-74, NP_fibrocartilaginous down)
- Mitotic Spindle (padj=1.20e-28, NP_fibrocartilaginous down)

**Reactome_2022** (top 5):

- Cell Cycle R-HSA-1640170 (padj=1.78e-94, NP_fibrocartilaginous down)
- Cell Cycle, Mitotic R-HSA-69278 (padj=5.76e-89, NP_fibrocartilaginous down)
- Cell Cycle R-HSA-1640170 (padj=5.64e-85, NP_fibrocartilaginous down)
- Cell Cycle, Mitotic R-HSA-69278 (padj=3.99e-79, NP_fibrocartilaginous down)
- Cell Cycle Checkpoints R-HSA-69620 (padj=2.07e-51, NP_fibrocartilaginous down)

**GSEA:** 58,837 significant terms. *[source: `results/interpretation/pathway_enrichment/gsea_results.tsv`]*

## 7. Transcription Factor Activity {#tf}

**288 significant TF-condition associations** (padj < 0.05). *[source: `results/interpretation/tf_activity/tf_activity_results.tsv`]*

**185 unique TFs** with significant activity changes.

| TF | Score | Cell Type | Comparison |
|----|-------|-----------|------------|
| HCFC1 | -0.571 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| E2F7 | 0.500 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| ARID3A | -0.500 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| HCFC1 | -0.429 | NP_fibrocartilaginous | mild_vs_severe |
| STOX1 | -0.400 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| HES6 | -0.400 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| STOX1 | -0.400 | NP_fibrocartilaginous | mild_vs_severe |
| E2F5 | -0.364 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| GLIS3 | 0.333 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| ARID3A | -0.333 | NP_fibrocartilaginous | mild_vs_severe |

## 8. Cell State Trajectories {#trajectory}

PAGA + diffusion pseudotime (DPT) analysis. *[source: `results/trajectories`]*

| compartment | embedding | test | rho | pvalue | n | statistic | n_healthy | n_degen | median_healthy | median_degen |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| AF | scVI_mesenchymal | pseudotime_vs_condition_ordinal | 0.1947 | 0.0000 | 5.00e+04 | nan | nan | nan | nan | nan |
| AF | scVI_mesenchymal | pseudotime_healthy_vs_degenerated | nan | 0.0000 | nan | 2.77e+08 | 2.48e+04 | 2.52e+04 | 0.6587 | 0.6587 |
| AF | scVI_mesenchymal | pseudotime_vs_condition_AF_outer | 0.1096 | 0.0000 | 3.78e+04 | nan | nan | nan | nan | nan |
| AF | scVI_mesenchymal | pseudotime_vs_condition_AF_inner | 0.2178 | 0.0000 | 1.22e+04 | nan | nan | nan | nan | nan |
| CEP | scVI_mesenchymal | pseudotime_vs_condition_ordinal | 0.0734 | 0.0000 | 3.21e+04 | nan | nan | nan | nan | nan |
| CEP | scVI_mesenchymal | pseudotime_healthy_vs_degenerated | nan | 0.0000 | nan | 1.08e+08 | 2.05e+04 | 1.17e+04 | 0.0652 | 0.0662 |
| CEP | scVI_mesenchymal | pseudotime_vs_condition_Fibrochondrocyte_chondroid | -0.2494 | 0.0000 | 2.21e+03 | nan | nan | nan | nan | nan |
| CEP | scVI_mesenchymal | pseudotime_vs_condition_EP_hyaline | 0.1373 | 0.0000 | 7.42e+03 | nan | nan | nan | nan | nan |
| CEP | scVI_mesenchymal | pseudotime_vs_condition_Fibroblast_like | -0.3059 | 0.0000 | 2.24e+04 | nan | nan | nan | nan | nan |
| CEP | scVI_mesenchymal | pseudotime_vs_condition_Fibrochondrocyte_fibroid | -0.4296 | 0.0000 | 1.62e+02 | nan | nan | nan | nan | nan |
| NP | scVI_mesenchymal | pseudotime_vs_condition_ordinal | -0.0882 | 0.0000 | 5.00e+04 | nan | nan | nan | nan | nan |
| NP | scVI_mesenchymal | pseudotime_healthy_vs_degenerated | nan | 0.0000 | nan | 2.66e+08 | 1.46e+04 | 2.99e+04 | 0.0805 | 0.0503 |
| NP | scVI_mesenchymal | pseudotime_vs_condition_NP_mature_chondrocyte | -0.0018 | 0.7709 | 2.70e+04 | nan | nan | nan | nan | nan |
| NP | scVI_mesenchymal | pseudotime_vs_condition_NP_fibrocartilaginous | -0.2017 | 0.0000 | 2.30e+04 | nan | nan | nan | nan | nan |

### Trajectory-associated genes

*[source: `results/trajectories`]*

- **AF:** 500 genes correlated with pseudotime (FDR < 0.05)
- **CEP:** 500 genes correlated with pseudotime (FDR < 0.05)
- **NP:** 500 genes correlated with pseudotime (FDR < 0.05)

## 9. Cell-Cell Communication {#communication}

LIANA rank_aggregate (CellPhoneDB, NATMI, Connectome, SingleCellSignalR, log2FC). *[source: `results/communication`]*

| Condition | Interactions |
|-----------|-------------|
| degenerated | 34,208 |
| healthy | 25,537 |

### Differential interactions

*[source: `results/communication/differential_interactions.tsv`]*

### Pain-relevant interactions

3075 pain-relevant L-R interactions identified. *[source: `results/communication/pain_interactions.tsv`]*

## 10. Pain Biology {#pain}

**1037 pain-associated DE genes identified.** *[source: `results/interpretation/pain_genes.tsv`]*

| gene | padj | cell_type | comparison |
| --- | --- | --- | --- |
| ANGPT1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| ANGPT2 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| ASIC1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| ASIC2 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| ASIC3 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| BDKRB1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| BDKRB2 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| BDNF | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| CALCA | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| CALCB | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| CCL2 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| CXCL8 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| FGF2 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| FLT1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| GAL | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| IL1B | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| IL6 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| KDR | nan | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| KLK1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| NGF | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| NGFR | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| NPY | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| NRP1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| NRP2 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| NTF3 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| NTN1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| NTN4 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| NTRK2 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| NTRK3 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| OPRD1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| OPRK1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| P2RX4 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| P2RX7 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| PDGFA | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| PDGFB | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| PENK | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| PLA2G2A | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| PTGES | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| PTGS2 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| ROBO1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| ROBO2 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| SCN11A | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| SCN9A | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| SEMA3A | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| SEMA3E | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| SLIT2 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| SLIT3 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| TAC1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| TEK | nan | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| TNF | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| TRPV1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| TRPV4 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| UNC5B | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| VEGFA | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| VEGFB | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| VIP | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| TRPA1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| VEGFC | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| SEMA3F | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| SLIT1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| DCC | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| NTRK1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| OPRM1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_all |
| ANGPT1 | 1.000 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| ANGPT2 | 0.794 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| ASIC1 | 0.768 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| ASIC2 | nan | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| ASIC3 | 1.000 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| BDKRB1 | 0.929 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| BDKRB2 | 1.000 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| BDNF | 0.155 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| CALCA | nan | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| CALCB | nan | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| CCL2 | 1.000 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| CXCL8 | 0.226 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| FGF2 | 0.372 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| FLT1 | 1.000 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| GAL | 1.000 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| IL1B | 7.83e-03 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| IL6 | 1.000 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| KDR | 1.000 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| KLK1 | 1.000 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| NGF | 0.932 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| NGFR | 0.826 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| NPY | 0.821 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| NRP1 | 1.000 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| NRP2 | 0.714 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| NTF3 | 1.000 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| NTN1 | 0.031 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| NTN4 | 0.413 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| NTRK2 | 0.691 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| NTRK3 | 0.791 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| OPRD1 | 1.000 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| P2RX4 | 0.813 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| P2RX7 | 0.507 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| PDGFA | 0.452 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| PDGFB | 0.651 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| PENK | 0.251 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| PLA2G2A | 7.83e-03 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| PTGES | 1.000 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| PTGS2 | 1.000 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| ROBO1 | 0.578 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| ROBO2 | 1.000 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| SCN11A | 1.000 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| SCN9A | 0.850 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| SEMA3A | 1.000 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| SEMA3E | 0.879 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| SLIT2 | 0.998 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| SLIT3 | 1.000 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| TAC1 | 0.740 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| TEK | 0.682 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| TNF | 0.740 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| TRPV1 | 1.000 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| TRPV4 | 1.000 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| UNC5B | 0.507 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| VEGFA | 0.547 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| VEGFB | 0.738 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| VIP | nan | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| TRPA1 | 0.775 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| VEGFC | 0.864 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| SEMA3F | 0.907 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| SLIT1 | nan | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| DCC | nan | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| NTRK1 | 1.000 | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| OPRM1 | nan | NP_fibrocartilaginous | healthy_vs_degenerated_all |
| ANGPT1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| ANGPT2 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| ASIC1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| ASIC2 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| ASIC3 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| BDKRB1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| BDKRB2 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| BDNF | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| CALCA | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| CALCB | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| CCL2 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| CXCL8 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| FGF2 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| FLT1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| GAL | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| IL1B | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| IL6 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| KDR | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| KLK1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| NGF | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| NGFR | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| NPY | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| NRP1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| NRP2 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| NTF3 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| NTN1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| NTN4 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| NTRK2 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| NTRK3 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| OPRD1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| OPRK1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| P2RX4 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| P2RX7 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| PDGFA | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| PDGFB | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| PENK | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| PLA2G2A | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| PTGES | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| PTGS2 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| ROBO1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| ROBO2 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| SCN11A | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| SCN9A | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| SEMA3A | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| SEMA3E | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| SLIT2 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| SLIT3 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| TAC1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| TEK | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| TNF | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| TRPV1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| TRPV4 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| UNC5B | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| VEGFA | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| VEGFB | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| VIP | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| TRPA1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| VEGFC | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| SEMA3F | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| SLIT1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| DCC | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| NTRK1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| OPRM1 | 1.000 | NP_mature_chondrocyte | healthy_vs_degenerated_mild |
| ANGPT1 | 0.833 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| ANGPT2 | 0.935 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| ASIC1 | 0.813 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| ASIC2 | nan | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| ASIC3 | nan | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| BDKRB1 | 0.990 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| BDKRB2 | 0.868 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| BDNF | 0.269 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| CALCA | nan | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| CALCB | nan | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| CCL2 | 0.967 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| CXCL8 | 0.441 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| FGF2 | 0.320 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| FLT1 | 0.929 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| GAL | 0.730 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| IL1B | 0.553 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| IL6 | 0.604 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| KDR | nan | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| KLK1 | nan | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| NGF | nan | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| NGFR | 0.848 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| NPY | nan | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| NRP1 | 0.998 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| NRP2 | 0.790 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| NTF3 | 0.912 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| NTN1 | 0.191 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| NTN4 | 0.554 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| NTRK2 | 0.823 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| NTRK3 | nan | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| OPRD1 | nan | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| P2RX4 | 0.781 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| P2RX7 | 0.741 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| PDGFA | 0.582 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| PDGFB | 0.741 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| PENK | 0.580 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| PLA2G2A | 0.100 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| PTGES | 0.951 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| PTGS2 | 0.832 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| ROBO1 | 0.704 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| ROBO2 | 0.918 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| SCN11A | nan | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| SCN9A | 0.888 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| SEMA3A | 0.908 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| SEMA3E | 0.702 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| SLIT2 | 0.817 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| SLIT3 | 0.901 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| TAC1 | nan | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| TEK | nan | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| TNF | 0.491 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| TRPV1 | 0.934 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| TRPV4 | 0.962 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| UNC5B | 0.703 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| VEGFA | 0.964 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| VEGFB | 0.758 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| VIP | nan | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| TRPA1 | nan | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| VEGFC | 0.729 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| SEMA3F | 0.859 | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| SLIT1 | nan | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| DCC | nan | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| NTRK1 | nan | NP_fibrocartilaginous | healthy_vs_degenerated_mild |
| ANGPT1 | 0.779 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| ANGPT2 | 0.442 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| ASIC1 | 0.514 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| ASIC2 | nan | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| ASIC3 | nan | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| BDKRB1 | 0.265 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| BDKRB2 | 0.370 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| BDNF | 0.245 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| CALCA | nan | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| CALCB | nan | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| CCL2 | 0.622 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| CXCL8 | 0.174 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| FGF2 | 0.074 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| FLT1 | 0.508 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| GAL | 0.994 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| IL1B | 0.410 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| IL6 | 0.440 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| KDR | nan | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| KLK1 | nan | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| NGF | 0.890 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| NGFR | 0.629 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| NPY | nan | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| NRP1 | 0.872 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| NRP2 | 0.255 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| NTF3 | 0.991 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| NTN1 | 1.17e-04 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| NTN4 | 0.045 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| NTRK2 | 0.284 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| NTRK3 | 0.439 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| OPRD1 | nan | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| P2RX4 | 0.669 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| P2RX7 | 0.124 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| PDGFA | 0.030 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| PDGFB | nan | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| PENK | 0.044 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| PLA2G2A | 6.65e-03 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| PTGES | 0.484 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| PTGS2 | 0.633 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| ROBO1 | 0.078 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| ROBO2 | 0.923 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| SCN11A | nan | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| SCN9A | 0.731 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| SEMA3A | 0.800 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| SEMA3E | 0.665 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| SLIT2 | 0.963 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| SLIT3 | 0.989 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| TAC1 | 0.645 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| TEK | nan | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| TNF | 0.972 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| TRPV1 | 0.963 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| TRPV4 | 0.960 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| UNC5B | 6.65e-03 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| VEGFA | 0.017 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| VEGFB | 0.383 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| VIP | nan | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| TRPA1 | nan | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| VEGFC | 0.884 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| SEMA3F | 0.917 | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| SLIT1 | nan | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| DCC | nan | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| NTRK1 | nan | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| OPRM1 | nan | NP_fibrocartilaginous | healthy_vs_degenerated_severe |
| ANGPT1 | 0.955 | NP_mature_chondrocyte | mild_vs_severe |
| ANGPT2 | 0.939 | NP_mature_chondrocyte | mild_vs_severe |
| ASIC1 | 0.994 | NP_mature_chondrocyte | mild_vs_severe |
| ASIC2 | nan | NP_mature_chondrocyte | mild_vs_severe |
| ASIC3 | nan | NP_mature_chondrocyte | mild_vs_severe |
| BDKRB1 | 0.599 | NP_mature_chondrocyte | mild_vs_severe |
| BDKRB2 | 0.203 | NP_mature_chondrocyte | mild_vs_severe |
| BDNF | 0.985 | NP_mature_chondrocyte | mild_vs_severe |
| CALCA | nan | NP_mature_chondrocyte | mild_vs_severe |
| CALCB | nan | NP_mature_chondrocyte | mild_vs_severe |
| CCL2 | 0.291 | NP_mature_chondrocyte | mild_vs_severe |
| CXCL8 | 0.732 | NP_mature_chondrocyte | mild_vs_severe |
| FGF2 | 0.414 | NP_mature_chondrocyte | mild_vs_severe |
| FLT1 | 0.816 | NP_mature_chondrocyte | mild_vs_severe |
| GAL | 0.436 | NP_mature_chondrocyte | mild_vs_severe |
| IL1B | 0.969 | NP_mature_chondrocyte | mild_vs_severe |
| IL6 | 0.161 | NP_mature_chondrocyte | mild_vs_severe |
| KDR | 0.939 | NP_mature_chondrocyte | mild_vs_severe |
| KLK1 | nan | NP_mature_chondrocyte | mild_vs_severe |
| NGF | 0.829 | NP_mature_chondrocyte | mild_vs_severe |
| NGFR | 0.802 | NP_mature_chondrocyte | mild_vs_severe |
| NPY | nan | NP_mature_chondrocyte | mild_vs_severe |
| NRP1 | 0.502 | NP_mature_chondrocyte | mild_vs_severe |
| NRP2 | 0.721 | NP_mature_chondrocyte | mild_vs_severe |
| NTF3 | 0.386 | NP_mature_chondrocyte | mild_vs_severe |
| NTN1 | 0.239 | NP_mature_chondrocyte | mild_vs_severe |
| NTN4 | 0.606 | NP_mature_chondrocyte | mild_vs_severe |
| NTRK2 | 0.674 | NP_mature_chondrocyte | mild_vs_severe |
| NTRK3 | nan | NP_mature_chondrocyte | mild_vs_severe |
| OPRD1 | nan | NP_mature_chondrocyte | mild_vs_severe |
| OPRK1 | nan | NP_mature_chondrocyte | mild_vs_severe |
| P2RX4 | 0.886 | NP_mature_chondrocyte | mild_vs_severe |
| P2RX7 | 0.817 | NP_mature_chondrocyte | mild_vs_severe |
| PDGFA | 0.059 | NP_mature_chondrocyte | mild_vs_severe |
| PDGFB | 0.649 | NP_mature_chondrocyte | mild_vs_severe |
| PENK | 0.999 | NP_mature_chondrocyte | mild_vs_severe |
| PLA2G2A | 0.931 | NP_mature_chondrocyte | mild_vs_severe |
| PTGES | 0.776 | NP_mature_chondrocyte | mild_vs_severe |
| PTGS2 | 0.994 | NP_mature_chondrocyte | mild_vs_severe |
| ROBO1 | 0.732 | NP_mature_chondrocyte | mild_vs_severe |
| ROBO2 | 0.904 | NP_mature_chondrocyte | mild_vs_severe |
| SCN11A | nan | NP_mature_chondrocyte | mild_vs_severe |
| SCN9A | 0.992 | NP_mature_chondrocyte | mild_vs_severe |
| SEMA3A | 0.974 | NP_mature_chondrocyte | mild_vs_severe |
| SEMA3E | 0.900 | NP_mature_chondrocyte | mild_vs_severe |
| SLIT2 | 0.678 | NP_mature_chondrocyte | mild_vs_severe |
| SLIT3 | 0.816 | NP_mature_chondrocyte | mild_vs_severe |
| TAC1 | 0.962 | NP_mature_chondrocyte | mild_vs_severe |
| TEK | 0.803 | NP_mature_chondrocyte | mild_vs_severe |
| TNF | 0.566 | NP_mature_chondrocyte | mild_vs_severe |
| TRPV1 | 0.984 | NP_mature_chondrocyte | mild_vs_severe |
| TRPV4 | 0.482 | NP_mature_chondrocyte | mild_vs_severe |
| UNC5B | 0.893 | NP_mature_chondrocyte | mild_vs_severe |
| VEGFA | 0.994 | NP_mature_chondrocyte | mild_vs_severe |
| VEGFB | 0.782 | NP_mature_chondrocyte | mild_vs_severe |
| VIP | nan | NP_mature_chondrocyte | mild_vs_severe |
| TRPA1 | nan | NP_mature_chondrocyte | mild_vs_severe |
| VEGFC | 0.326 | NP_mature_chondrocyte | mild_vs_severe |
| SEMA3F | 0.686 | NP_mature_chondrocyte | mild_vs_severe |
| SLIT1 | nan | NP_mature_chondrocyte | mild_vs_severe |
| DCC | nan | NP_mature_chondrocyte | mild_vs_severe |
| NTRK1 | nan | NP_mature_chondrocyte | mild_vs_severe |
| OPRM1 | nan | NP_mature_chondrocyte | mild_vs_severe |
| ANGPT1 | 0.221 | NP_fibrocartilaginous | mild_vs_severe |
| ANGPT2 | 0.371 | NP_fibrocartilaginous | mild_vs_severe |
| ASIC1 | 0.553 | NP_fibrocartilaginous | mild_vs_severe |
| ASIC2 | nan | NP_fibrocartilaginous | mild_vs_severe |
| ASIC3 | nan | NP_fibrocartilaginous | mild_vs_severe |
| BDKRB1 | 0.122 | NP_fibrocartilaginous | mild_vs_severe |
| BDKRB2 | 0.053 | NP_fibrocartilaginous | mild_vs_severe |
| BDNF | 0.823 | NP_fibrocartilaginous | mild_vs_severe |
| CALCA | nan | NP_fibrocartilaginous | mild_vs_severe |
| CALCB | nan | NP_fibrocartilaginous | mild_vs_severe |
| CCL2 | 0.177 | NP_fibrocartilaginous | mild_vs_severe |
| CXCL8 | 0.754 | NP_fibrocartilaginous | mild_vs_severe |
| FGF2 | 0.843 | NP_fibrocartilaginous | mild_vs_severe |
| FLT1 | 0.521 | NP_fibrocartilaginous | mild_vs_severe |
| GAL | 0.856 | NP_fibrocartilaginous | mild_vs_severe |
| IL1B | 0.321 | NP_fibrocartilaginous | mild_vs_severe |
| IL6 | 5.22e-03 | NP_fibrocartilaginous | mild_vs_severe |
| KDR | nan | NP_fibrocartilaginous | mild_vs_severe |
| KLK1 | nan | NP_fibrocartilaginous | mild_vs_severe |
| NGF | nan | NP_fibrocartilaginous | mild_vs_severe |
| NGFR | 0.985 | NP_fibrocartilaginous | mild_vs_severe |
| NPY | nan | NP_fibrocartilaginous | mild_vs_severe |
| NRP1 | 0.457 | NP_fibrocartilaginous | mild_vs_severe |
| NRP2 | 0.685 | NP_fibrocartilaginous | mild_vs_severe |
| NTF3 | 0.743 | NP_fibrocartilaginous | mild_vs_severe |
| NTN1 | 0.268 | NP_fibrocartilaginous | mild_vs_severe |
| NTN4 | 0.736 | NP_fibrocartilaginous | mild_vs_severe |
| NTRK2 | 0.103 | NP_fibrocartilaginous | mild_vs_severe |
| NTRK3 | nan | NP_fibrocartilaginous | mild_vs_severe |
| OPRD1 | nan | NP_fibrocartilaginous | mild_vs_severe |
| P2RX4 | 0.788 | NP_fibrocartilaginous | mild_vs_severe |
| P2RX7 | 0.894 | NP_fibrocartilaginous | mild_vs_severe |
| PDGFA | 0.217 | NP_fibrocartilaginous | mild_vs_severe |
| PDGFB | nan | NP_fibrocartilaginous | mild_vs_severe |
| PENK | 0.709 | NP_fibrocartilaginous | mild_vs_severe |
| PLA2G2A | 0.541 | NP_fibrocartilaginous | mild_vs_severe |
| PTGES | 0.278 | NP_fibrocartilaginous | mild_vs_severe |
| PTGS2 | 0.086 | NP_fibrocartilaginous | mild_vs_severe |
| ROBO1 | 0.501 | NP_fibrocartilaginous | mild_vs_severe |
| ROBO2 | 0.993 | NP_fibrocartilaginous | mild_vs_severe |
| SCN11A | nan | NP_fibrocartilaginous | mild_vs_severe |
| SCN9A | 0.512 | NP_fibrocartilaginous | mild_vs_severe |
| SEMA3A | 0.610 | NP_fibrocartilaginous | mild_vs_severe |
| SEMA3E | 0.750 | NP_fibrocartilaginous | mild_vs_severe |
| SLIT2 | 0.399 | NP_fibrocartilaginous | mild_vs_severe |
| SLIT3 | 0.430 | NP_fibrocartilaginous | mild_vs_severe |
| TAC1 | 0.462 | NP_fibrocartilaginous | mild_vs_severe |
| TEK | nan | NP_fibrocartilaginous | mild_vs_severe |
| TNF | 0.098 | NP_fibrocartilaginous | mild_vs_severe |
| TRPV1 | 0.706 | NP_fibrocartilaginous | mild_vs_severe |
| TRPV4 | 0.659 | NP_fibrocartilaginous | mild_vs_severe |
| UNC5B | 0.401 | NP_fibrocartilaginous | mild_vs_severe |
| VEGFA | 0.055 | NP_fibrocartilaginous | mild_vs_severe |
| VEGFB | 0.335 | NP_fibrocartilaginous | mild_vs_severe |
| VIP | nan | NP_fibrocartilaginous | mild_vs_severe |
| TRPA1 | 0.313 | NP_fibrocartilaginous | mild_vs_severe |
| VEGFC | 0.749 | NP_fibrocartilaginous | mild_vs_severe |
| SEMA3F | 0.910 | NP_fibrocartilaginous | mild_vs_severe |
| SLIT1 | nan | NP_fibrocartilaginous | mild_vs_severe |
| DCC | nan | NP_fibrocartilaginous | mild_vs_severe |
| NTRK1 | nan | NP_fibrocartilaginous | mild_vs_severe |
| OPRM1 | nan | NP_fibrocartilaginous | mild_vs_severe |
| ANGPT1 | 0.996 | Endothelial | mild_vs_severe |
| ANGPT2 | 0.996 | Endothelial | mild_vs_severe |
| ASIC3 | 0.996 | Endothelial | mild_vs_severe |
| BDKRB1 | 0.996 | Endothelial | mild_vs_severe |
| BDKRB2 | 0.996 | Endothelial | mild_vs_severe |
| CCL2 | 0.996 | Endothelial | mild_vs_severe |
| CXCL8 | 0.996 | Endothelial | mild_vs_severe |
| FGF2 | 0.996 | Endothelial | mild_vs_severe |
| FLT1 | 0.996 | Endothelial | mild_vs_severe |
| GAL | 0.996 | Endothelial | mild_vs_severe |
| IL1B | 0.996 | Endothelial | mild_vs_severe |
| IL6 | 0.996 | Endothelial | mild_vs_severe |
| KDR | 0.996 | Endothelial | mild_vs_severe |
| NGF | 0.996 | Endothelial | mild_vs_severe |
| NGFR | 0.996 | Endothelial | mild_vs_severe |
| NRP1 | 0.996 | Endothelial | mild_vs_severe |
| NRP2 | 0.996 | Endothelial | mild_vs_severe |
| NTF3 | 0.996 | Endothelial | mild_vs_severe |
| NTN1 | 0.996 | Endothelial | mild_vs_severe |
| NTN4 | 0.996 | Endothelial | mild_vs_severe |
| NTRK2 | 0.996 | Endothelial | mild_vs_severe |
| P2RX4 | 0.996 | Endothelial | mild_vs_severe |
| P2RX7 | 0.996 | Endothelial | mild_vs_severe |
| PDGFA | 0.996 | Endothelial | mild_vs_severe |
| PDGFB | 0.996 | Endothelial | mild_vs_severe |
| PENK | 0.996 | Endothelial | mild_vs_severe |
| PLA2G2A | 0.996 | Endothelial | mild_vs_severe |
| PTGES | 0.996 | Endothelial | mild_vs_severe |
| PTGS2 | 0.996 | Endothelial | mild_vs_severe |
| ROBO1 | 0.996 | Endothelial | mild_vs_severe |
| SCN11A | 0.996 | Endothelial | mild_vs_severe |
| SCN9A | 0.996 | Endothelial | mild_vs_severe |
| SEMA3A | 0.996 | Endothelial | mild_vs_severe |
| SEMA3E | 0.996 | Endothelial | mild_vs_severe |
| SLIT2 | 0.996 | Endothelial | mild_vs_severe |
| SLIT3 | 0.996 | Endothelial | mild_vs_severe |
| TAC1 | 0.996 | Endothelial | mild_vs_severe |
| TEK | 0.996 | Endothelial | mild_vs_severe |
| TNF | 0.996 | Endothelial | mild_vs_severe |
| TRPV1 | 0.996 | Endothelial | mild_vs_severe |
| TRPV4 | 0.996 | Endothelial | mild_vs_severe |
| UNC5B | 0.996 | Endothelial | mild_vs_severe |
| VEGFA | 0.998 | Endothelial | mild_vs_severe |
| VEGFB | 0.996 | Endothelial | mild_vs_severe |
| VIP | 0.996 | Endothelial | mild_vs_severe |
| VEGFC | 0.996 | Endothelial | mild_vs_severe |
| SEMA3F | 0.996 | Endothelial | mild_vs_severe |
| ANGPT1 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| ANGPT2 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| ASIC1 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| ASIC2 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| ASIC3 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| BDKRB1 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| BDKRB2 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| BDNF | 1.000 | AF_inner | healthy_vs_degenerated_all |
| CALCA | 1.000 | AF_inner | healthy_vs_degenerated_all |
| CALCB | 1.000 | AF_inner | healthy_vs_degenerated_all |
| CCL2 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| CXCL8 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| FGF2 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| FLT1 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| GAL | 1.000 | AF_inner | healthy_vs_degenerated_all |
| IL1B | 1.000 | AF_inner | healthy_vs_degenerated_all |
| IL6 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| KDR | 1.000 | AF_inner | healthy_vs_degenerated_all |
| KLK1 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| NGF | 1.000 | AF_inner | healthy_vs_degenerated_all |
| NGFR | 1.000 | AF_inner | healthy_vs_degenerated_all |
| NPY | 1.000 | AF_inner | healthy_vs_degenerated_all |
| NRP1 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| NRP2 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| NTF3 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| NTN1 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| NTN4 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| NTRK2 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| NTRK3 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| OPRD1 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| OPRK1 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| P2RX4 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| P2RX7 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| PDGFA | 1.000 | AF_inner | healthy_vs_degenerated_all |
| PDGFB | 1.000 | AF_inner | healthy_vs_degenerated_all |
| PENK | 1.000 | AF_inner | healthy_vs_degenerated_all |
| PLA2G2A | 1.000 | AF_inner | healthy_vs_degenerated_all |
| PTGES | 1.000 | AF_inner | healthy_vs_degenerated_all |
| PTGS2 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| ROBO1 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| ROBO2 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| SCN11A | 1.000 | AF_inner | healthy_vs_degenerated_all |
| SCN9A | 1.000 | AF_inner | healthy_vs_degenerated_all |
| SEMA3A | 1.000 | AF_inner | healthy_vs_degenerated_all |
| SEMA3E | 1.000 | AF_inner | healthy_vs_degenerated_all |
| SLIT2 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| SLIT3 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| TAC1 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| TEK | 1.000 | AF_inner | healthy_vs_degenerated_all |
| TNF | 1.000 | AF_inner | healthy_vs_degenerated_all |
| TRPV1 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| TRPV4 | 1.000 | AF_inner | healthy_vs_degenerated_all |
| UNC5B | 1.000 | AF_inner | healthy_vs_degenerated_all |
| VEGFA | 1.000 | AF_inner | healthy_vs_degenerated_all |
| VEGFB | 1.000 | AF_inner | healthy_vs_degenerated_all |
| VIP | 1.000 | AF_inner | healthy_vs_degenerated_all |
| SEMA3F | 1.000 | AF_inner | healthy_vs_degenerated_all |
| VEGFC | 1.000 | AF_inner | healthy_vs_degenerated_all |
| ANGPT1 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| ANGPT2 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| ASIC1 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| ASIC2 | 0.997 | AF_outer | healthy_vs_degenerated_all |
| ASIC3 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| BDKRB1 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| BDKRB2 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| BDNF | 0.996 | AF_outer | healthy_vs_degenerated_all |
| CALCA | 0.996 | AF_outer | healthy_vs_degenerated_all |
| CALCB | 0.996 | AF_outer | healthy_vs_degenerated_all |
| CCL2 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| CXCL8 | nan | AF_outer | healthy_vs_degenerated_all |
| FGF2 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| FLT1 | 1.000 | AF_outer | healthy_vs_degenerated_all |
| GAL | 0.996 | AF_outer | healthy_vs_degenerated_all |
| IL1B | 0.996 | AF_outer | healthy_vs_degenerated_all |
| IL6 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| KDR | 0.996 | AF_outer | healthy_vs_degenerated_all |
| KLK1 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| NGF | 0.996 | AF_outer | healthy_vs_degenerated_all |
| NGFR | 0.996 | AF_outer | healthy_vs_degenerated_all |
| NPY | 0.997 | AF_outer | healthy_vs_degenerated_all |
| NRP1 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| NRP2 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| NTF3 | 0.997 | AF_outer | healthy_vs_degenerated_all |
| NTN1 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| NTN4 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| NTRK2 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| NTRK3 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| OPRD1 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| OPRK1 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| P2RX4 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| P2RX7 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| PDGFA | 0.996 | AF_outer | healthy_vs_degenerated_all |
| PDGFB | 0.997 | AF_outer | healthy_vs_degenerated_all |
| PENK | 0.996 | AF_outer | healthy_vs_degenerated_all |
| PLA2G2A | 0.996 | AF_outer | healthy_vs_degenerated_all |
| PTGES | 0.996 | AF_outer | healthy_vs_degenerated_all |
| PTGS2 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| ROBO1 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| ROBO2 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| SCN11A | 0.997 | AF_outer | healthy_vs_degenerated_all |
| SCN9A | 0.997 | AF_outer | healthy_vs_degenerated_all |
| SEMA3A | 0.996 | AF_outer | healthy_vs_degenerated_all |
| SEMA3E | 0.997 | AF_outer | healthy_vs_degenerated_all |
| SLIT2 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| SLIT3 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| TAC1 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| TEK | 0.996 | AF_outer | healthy_vs_degenerated_all |
| TNF | 0.996 | AF_outer | healthy_vs_degenerated_all |
| TRPV1 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| TRPV4 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| UNC5B | 0.996 | AF_outer | healthy_vs_degenerated_all |
| VEGFA | 0.996 | AF_outer | healthy_vs_degenerated_all |
| VEGFB | 0.996 | AF_outer | healthy_vs_degenerated_all |
| VIP | 0.999 | AF_outer | healthy_vs_degenerated_all |
| P2RX3 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| SEMA3F | 0.997 | AF_outer | healthy_vs_degenerated_all |
| SLIT1 | 0.997 | AF_outer | healthy_vs_degenerated_all |
| VEGFC | 0.996 | AF_outer | healthy_vs_degenerated_all |
| DCC | 0.996 | AF_outer | healthy_vs_degenerated_all |
| NTRK1 | 0.999 | AF_outer | healthy_vs_degenerated_all |
| OPRM1 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| TRPA1 | 0.996 | AF_outer | healthy_vs_degenerated_all |
| ANGPT1 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| ANGPT2 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| ASIC1 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| ASIC2 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| ASIC3 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| BDKRB1 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| BDKRB2 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| BDNF | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| CALCA | nan | AF_inner | healthy_vs_degenerated_mild |
| CALCB | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| CCL2 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| CXCL8 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| FGF2 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| FLT1 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| GAL | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| IL1B | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| IL6 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| KDR | nan | AF_inner | healthy_vs_degenerated_mild |
| KLK1 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| NGF | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| NGFR | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| NPY | nan | AF_inner | healthy_vs_degenerated_mild |
| NRP1 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| NRP2 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| NTF3 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| NTN1 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| NTN4 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| NTRK2 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| NTRK3 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| OPRD1 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| OPRK1 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| P2RX4 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| P2RX7 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| PDGFA | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| PDGFB | nan | AF_inner | healthy_vs_degenerated_mild |
| PENK | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| PLA2G2A | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| PTGES | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| PTGS2 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| ROBO1 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| ROBO2 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| SCN11A | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| SCN9A | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| SEMA3A | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| SEMA3E | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| SLIT2 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| SLIT3 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| TAC1 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| TEK | nan | AF_inner | healthy_vs_degenerated_mild |
| TNF | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| TRPV1 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| TRPV4 | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| UNC5B | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| VEGFA | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| VEGFB | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| VIP | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| SEMA3F | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| VEGFC | 1.000 | AF_inner | healthy_vs_degenerated_mild |
| ANGPT1 | 0.836 | AF_outer | healthy_vs_degenerated_mild |
| ANGPT2 | 0.160 | AF_outer | healthy_vs_degenerated_mild |
| ASIC1 | 0.956 | AF_outer | healthy_vs_degenerated_mild |
| ASIC2 | nan | AF_outer | healthy_vs_degenerated_mild |
| ASIC3 | nan | AF_outer | healthy_vs_degenerated_mild |
| BDKRB1 | 0.921 | AF_outer | healthy_vs_degenerated_mild |
| BDKRB2 | 0.748 | AF_outer | healthy_vs_degenerated_mild |
| BDNF | 0.721 | AF_outer | healthy_vs_degenerated_mild |
| CALCA | nan | AF_outer | healthy_vs_degenerated_mild |
| CALCB | nan | AF_outer | healthy_vs_degenerated_mild |
| CCL2 | 0.623 | AF_outer | healthy_vs_degenerated_mild |
| CXCL8 | 7.58e-03 | AF_outer | healthy_vs_degenerated_mild |
| FGF2 | 0.546 | AF_outer | healthy_vs_degenerated_mild |
| FLT1 | 0.156 | AF_outer | healthy_vs_degenerated_mild |
| GAL | 0.795 | AF_outer | healthy_vs_degenerated_mild |
| IL1B | nan | AF_outer | healthy_vs_degenerated_mild |
| IL6 | 0.985 | AF_outer | healthy_vs_degenerated_mild |
| KDR | nan | AF_outer | healthy_vs_degenerated_mild |
| KLK1 | nan | AF_outer | healthy_vs_degenerated_mild |
| NGF | 0.832 | AF_outer | healthy_vs_degenerated_mild |
| NGFR | 0.054 | AF_outer | healthy_vs_degenerated_mild |
| NPY | nan | AF_outer | healthy_vs_degenerated_mild |
| NRP1 | 0.766 | AF_outer | healthy_vs_degenerated_mild |
| NRP2 | 0.778 | AF_outer | healthy_vs_degenerated_mild |
| NTF3 | 0.932 | AF_outer | healthy_vs_degenerated_mild |
| NTN1 | 0.685 | AF_outer | healthy_vs_degenerated_mild |
| NTN4 | 0.757 | AF_outer | healthy_vs_degenerated_mild |
| NTRK2 | 0.993 | AF_outer | healthy_vs_degenerated_mild |
| NTRK3 | nan | AF_outer | healthy_vs_degenerated_mild |
| OPRD1 | nan | AF_outer | healthy_vs_degenerated_mild |
| OPRK1 | nan | AF_outer | healthy_vs_degenerated_mild |
| P2RX4 | 0.891 | AF_outer | healthy_vs_degenerated_mild |
| P2RX7 | 0.550 | AF_outer | healthy_vs_degenerated_mild |
| PDGFA | 0.577 | AF_outer | healthy_vs_degenerated_mild |
| PDGFB | nan | AF_outer | healthy_vs_degenerated_mild |
| PENK | 0.464 | AF_outer | healthy_vs_degenerated_mild |
| PLA2G2A | 8.62e-03 | AF_outer | healthy_vs_degenerated_mild |
| PTGES | 0.961 | AF_outer | healthy_vs_degenerated_mild |
| PTGS2 | 0.273 | AF_outer | healthy_vs_degenerated_mild |
| ROBO1 | 0.797 | AF_outer | healthy_vs_degenerated_mild |
| ROBO2 | 0.975 | AF_outer | healthy_vs_degenerated_mild |
| SCN11A | nan | AF_outer | healthy_vs_degenerated_mild |
| SCN9A | nan | AF_outer | healthy_vs_degenerated_mild |
| SEMA3A | 0.975 | AF_outer | healthy_vs_degenerated_mild |
| SEMA3E | 0.854 | AF_outer | healthy_vs_degenerated_mild |
| SLIT2 | 0.748 | AF_outer | healthy_vs_degenerated_mild |
| SLIT3 | 0.933 | AF_outer | healthy_vs_degenerated_mild |
| TAC1 | nan | AF_outer | healthy_vs_degenerated_mild |
| TEK | nan | AF_outer | healthy_vs_degenerated_mild |
| TNF | nan | AF_outer | healthy_vs_degenerated_mild |
| TRPV1 | 0.874 | AF_outer | healthy_vs_degenerated_mild |
| TRPV4 | 0.112 | AF_outer | healthy_vs_degenerated_mild |
| UNC5B | 0.845 | AF_outer | healthy_vs_degenerated_mild |
| VEGFA | 0.648 | AF_outer | healthy_vs_degenerated_mild |
| VEGFB | 0.664 | AF_outer | healthy_vs_degenerated_mild |
| VIP | nan | AF_outer | healthy_vs_degenerated_mild |
| P2RX3 | nan | AF_outer | healthy_vs_degenerated_mild |
| SEMA3F | 0.955 | AF_outer | healthy_vs_degenerated_mild |
| SLIT1 | nan | AF_outer | healthy_vs_degenerated_mild |
| VEGFC | 0.521 | AF_outer | healthy_vs_degenerated_mild |
| DCC | nan | AF_outer | healthy_vs_degenerated_mild |
| NTRK1 | nan | AF_outer | healthy_vs_degenerated_mild |
| OPRM1 | nan | AF_outer | healthy_vs_degenerated_mild |
| TRPA1 | nan | AF_outer | healthy_vs_degenerated_mild |
| ANGPT1 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| ANGPT2 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| ASIC1 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| ASIC2 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| ASIC3 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| BDKRB1 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| BDKRB2 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| BDNF | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| CALCA | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| CALCB | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| CCL2 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| CXCL8 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| FGF2 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| FLT1 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| GAL | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| IL1B | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| IL6 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| KDR | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| KLK1 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| NGF | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| NGFR | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| NPY | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| NRP1 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| NRP2 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| NTF3 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| NTN1 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| NTN4 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| NTRK2 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| NTRK3 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| OPRD1 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| OPRK1 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| P2RX4 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| P2RX7 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| PDGFA | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| PDGFB | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| PENK | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| PLA2G2A | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| PTGES | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| PTGS2 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| ROBO1 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| ROBO2 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| SCN11A | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| SCN9A | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| SEMA3A | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| SEMA3E | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| SLIT2 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| SLIT3 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| TAC1 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| TEK | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| TNF | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| TRPV1 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| TRPV4 | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| UNC5B | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| VEGFA | 0.892 | AF_inner | healthy_vs_degenerated_severe |
| VEGFB | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| VIP | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| SEMA3F | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| VEGFC | 1.000 | AF_inner | healthy_vs_degenerated_severe |
| ANGPT1 | 0.936 | AF_outer | healthy_vs_degenerated_severe |
| ANGPT2 | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| ASIC1 | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| ASIC2 | nan | AF_outer | healthy_vs_degenerated_severe |
| ASIC3 | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| BDKRB1 | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| BDKRB2 | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| BDNF | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| CALCA | 0.917 | AF_outer | healthy_vs_degenerated_severe |
| CALCB | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| CCL2 | 0.306 | AF_outer | healthy_vs_degenerated_severe |
| CXCL8 | nan | AF_outer | healthy_vs_degenerated_severe |
| FGF2 | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| FLT1 | 0.967 | AF_outer | healthy_vs_degenerated_severe |
| GAL | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| IL1B | nan | AF_outer | healthy_vs_degenerated_severe |
| IL6 | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| KDR | nan | AF_outer | healthy_vs_degenerated_severe |
| KLK1 | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| NGF | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| NGFR | 0.974 | AF_outer | healthy_vs_degenerated_severe |
| NPY | 0.963 | AF_outer | healthy_vs_degenerated_severe |
| NRP1 | 0.958 | AF_outer | healthy_vs_degenerated_severe |
| NRP2 | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| NTF3 | 0.926 | AF_outer | healthy_vs_degenerated_severe |
| NTN1 | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| NTN4 | 0.961 | AF_outer | healthy_vs_degenerated_severe |
| NTRK2 | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| NTRK3 | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| OPRD1 | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| OPRK1 | nan | AF_outer | healthy_vs_degenerated_severe |
| P2RX4 | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| P2RX7 | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| PDGFA | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| PDGFB | nan | AF_outer | healthy_vs_degenerated_severe |
| PENK | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| PLA2G2A | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| PTGES | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| PTGS2 | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| ROBO1 | 0.942 | AF_outer | healthy_vs_degenerated_severe |
| ROBO2 | 0.945 | AF_outer | healthy_vs_degenerated_severe |
| SCN11A | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| SCN9A | 0.936 | AF_outer | healthy_vs_degenerated_severe |
| SEMA3A | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| SEMA3E | 0.941 | AF_outer | healthy_vs_degenerated_severe |
| SLIT2 | 0.935 | AF_outer | healthy_vs_degenerated_severe |
| SLIT3 | 0.813 | AF_outer | healthy_vs_degenerated_severe |
| TAC1 | 0.945 | AF_outer | healthy_vs_degenerated_severe |
| TEK | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| TNF | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| TRPV1 | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| TRPV4 | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| UNC5B | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| VEGFA | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| VEGFB | 0.937 | AF_outer | healthy_vs_degenerated_severe |
| VIP | nan | AF_outer | healthy_vs_degenerated_severe |
| P2RX3 | nan | AF_outer | healthy_vs_degenerated_severe |
| SEMA3F | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| SLIT1 | nan | AF_outer | healthy_vs_degenerated_severe |
| VEGFC | 0.914 | AF_outer | healthy_vs_degenerated_severe |
| DCC | nan | AF_outer | healthy_vs_degenerated_severe |
| NTRK1 | nan | AF_outer | healthy_vs_degenerated_severe |
| OPRM1 | nan | AF_outer | healthy_vs_degenerated_severe |
| TRPA1 | nan | AF_outer | healthy_vs_degenerated_severe |
| ANGPT1 | 1.000 | AF_inner | mild_vs_severe |
| ANGPT2 | 1.000 | AF_inner | mild_vs_severe |
| ASIC1 | 1.000 | AF_inner | mild_vs_severe |
| ASIC2 | 1.000 | AF_inner | mild_vs_severe |
| ASIC3 | 1.000 | AF_inner | mild_vs_severe |
| BDKRB1 | 1.000 | AF_inner | mild_vs_severe |
| BDKRB2 | 1.000 | AF_inner | mild_vs_severe |
| BDNF | 1.000 | AF_inner | mild_vs_severe |
| CALCA | 1.000 | AF_inner | mild_vs_severe |
| CALCB | 1.000 | AF_inner | mild_vs_severe |
| CCL2 | 1.000 | AF_inner | mild_vs_severe |
| CXCL8 | nan | AF_inner | mild_vs_severe |
| FGF2 | 1.000 | AF_inner | mild_vs_severe |
| FLT1 | 1.000 | AF_inner | mild_vs_severe |
| GAL | 1.000 | AF_inner | mild_vs_severe |
| IL1B | 1.000 | AF_inner | mild_vs_severe |
| IL6 | 1.000 | AF_inner | mild_vs_severe |
| KDR | 1.000 | AF_inner | mild_vs_severe |
| KLK1 | 1.000 | AF_inner | mild_vs_severe |
| NGF | 1.000 | AF_inner | mild_vs_severe |
| NGFR | 1.000 | AF_inner | mild_vs_severe |
| NPY | 1.000 | AF_inner | mild_vs_severe |
| NRP1 | 1.000 | AF_inner | mild_vs_severe |
| NRP2 | 1.000 | AF_inner | mild_vs_severe |
| NTF3 | 1.000 | AF_inner | mild_vs_severe |
| NTN1 | 1.000 | AF_inner | mild_vs_severe |
| NTN4 | 1.000 | AF_inner | mild_vs_severe |
| NTRK2 | 1.000 | AF_inner | mild_vs_severe |
| NTRK3 | 1.000 | AF_inner | mild_vs_severe |
| OPRD1 | 1.000 | AF_inner | mild_vs_severe |
| OPRK1 | 1.000 | AF_inner | mild_vs_severe |
| P2RX4 | 1.000 | AF_inner | mild_vs_severe |
| P2RX7 | 1.000 | AF_inner | mild_vs_severe |
| PDGFA | 1.000 | AF_inner | mild_vs_severe |
| PDGFB | 1.000 | AF_inner | mild_vs_severe |
| PENK | 1.000 | AF_inner | mild_vs_severe |
| PLA2G2A | 0.957 | AF_inner | mild_vs_severe |
| PTGES | nan | AF_inner | mild_vs_severe |
| PTGS2 | 0.973 | AF_inner | mild_vs_severe |
| ROBO1 | 1.000 | AF_inner | mild_vs_severe |
| ROBO2 | 1.000 | AF_inner | mild_vs_severe |
| SCN11A | 1.000 | AF_inner | mild_vs_severe |
| SCN9A | 1.000 | AF_inner | mild_vs_severe |
| SEMA3A | 1.000 | AF_inner | mild_vs_severe |
| SEMA3E | 1.000 | AF_inner | mild_vs_severe |
| SLIT2 | 1.000 | AF_inner | mild_vs_severe |
| SLIT3 | 1.000 | AF_inner | mild_vs_severe |
| TAC1 | 1.000 | AF_inner | mild_vs_severe |
| TEK | 1.000 | AF_inner | mild_vs_severe |
| TNF | 1.000 | AF_inner | mild_vs_severe |
| TRPV1 | 1.000 | AF_inner | mild_vs_severe |
| TRPV4 | 1.000 | AF_inner | mild_vs_severe |
| UNC5B | 1.000 | AF_inner | mild_vs_severe |
| VEGFA | 1.000 | AF_inner | mild_vs_severe |
| VEGFB | 1.000 | AF_inner | mild_vs_severe |
| P2RX3 | 1.000 | AF_inner | mild_vs_severe |
| SEMA3F | 1.000 | AF_inner | mild_vs_severe |
| SLIT1 | 1.000 | AF_inner | mild_vs_severe |
| VEGFC | 1.000 | AF_inner | mild_vs_severe |
| DCC | 1.000 | AF_inner | mild_vs_severe |
| OPRM1 | 1.000 | AF_inner | mild_vs_severe |
| ANGPT1 | 1.000 | AF_outer | mild_vs_severe |
| ANGPT2 | 0.676 | AF_outer | mild_vs_severe |
| ASIC1 | 0.950 | AF_outer | mild_vs_severe |
| ASIC2 | 0.956 | AF_outer | mild_vs_severe |
| ASIC3 | 0.606 | AF_outer | mild_vs_severe |
| BDKRB1 | 0.857 | AF_outer | mild_vs_severe |
| BDKRB2 | 0.857 | AF_outer | mild_vs_severe |
| BDNF | 0.955 | AF_outer | mild_vs_severe |
| CALCA | 0.991 | AF_outer | mild_vs_severe |
| CALCB | 0.870 | AF_outer | mild_vs_severe |
| CCL2 | 0.828 | AF_outer | mild_vs_severe |
| CXCL8 | 0.684 | AF_outer | mild_vs_severe |
| FGF2 | 0.865 | AF_outer | mild_vs_severe |
| FLT1 | 0.857 | AF_outer | mild_vs_severe |
| GAL | 0.600 | AF_outer | mild_vs_severe |
| IL1B | 0.869 | AF_outer | mild_vs_severe |
| IL6 | 0.956 | AF_outer | mild_vs_severe |
| KDR | 0.965 | AF_outer | mild_vs_severe |
| KLK1 | 0.661 | AF_outer | mild_vs_severe |
| NGF | 0.828 | AF_outer | mild_vs_severe |
| NGFR | 0.778 | AF_outer | mild_vs_severe |
| NPY | 0.917 | AF_outer | mild_vs_severe |
| NRP1 | 0.917 | AF_outer | mild_vs_severe |
| NRP2 | 0.963 | AF_outer | mild_vs_severe |
| NTF3 | 0.870 | AF_outer | mild_vs_severe |
| NTN1 | 0.612 | AF_outer | mild_vs_severe |
| NTN4 | 0.927 | AF_outer | mild_vs_severe |
| NTRK2 | 0.893 | AF_outer | mild_vs_severe |
| NTRK3 | 1.000 | AF_outer | mild_vs_severe |
| OPRD1 | 0.828 | AF_outer | mild_vs_severe |
| OPRK1 | 0.917 | AF_outer | mild_vs_severe |
| P2RX4 | 0.863 | AF_outer | mild_vs_severe |
| P2RX7 | 0.925 | AF_outer | mild_vs_severe |
| PDGFA | 0.946 | AF_outer | mild_vs_severe |
| PDGFB | 0.870 | AF_outer | mild_vs_severe |
| PENK | 0.878 | AF_outer | mild_vs_severe |
| PLA2G2A | nan | AF_outer | mild_vs_severe |
| PTGES | 0.705 | AF_outer | mild_vs_severe |
| PTGS2 | 0.628 | AF_outer | mild_vs_severe |
| ROBO1 | 0.942 | AF_outer | mild_vs_severe |
| ROBO2 | 1.000 | AF_outer | mild_vs_severe |
| SCN11A | 0.840 | AF_outer | mild_vs_severe |
| SCN9A | 0.870 | AF_outer | mild_vs_severe |
| SEMA3A | 0.600 | AF_outer | mild_vs_severe |
| SEMA3E | 0.907 | AF_outer | mild_vs_severe |
| SLIT2 | 0.933 | AF_outer | mild_vs_severe |
| SLIT3 | 0.646 | AF_outer | mild_vs_severe |
| TAC1 | 0.828 | AF_outer | mild_vs_severe |
| TEK | 0.857 | AF_outer | mild_vs_severe |
| TRPV1 | 0.606 | AF_outer | mild_vs_severe |
| TRPV4 | 0.612 | AF_outer | mild_vs_severe |
| UNC5B | 1.000 | AF_outer | mild_vs_severe |
| VEGFA | 0.762 | AF_outer | mild_vs_severe |
| VEGFB | 0.951 | AF_outer | mild_vs_severe |
| VIP | 0.646 | AF_outer | mild_vs_severe |
| P2RX3 | 0.869 | AF_outer | mild_vs_severe |
| SEMA3F | 0.749 | AF_outer | mild_vs_severe |
| SLIT1 | 0.857 | AF_outer | mild_vs_severe |
| VEGFC | 0.857 | AF_outer | mild_vs_severe |
| DCC | 0.857 | AF_outer | mild_vs_severe |
| NTRK1 | 0.870 | AF_outer | mild_vs_severe |
| OPRM1 | 0.755 | AF_outer | mild_vs_severe |
| TRPA1 | 0.917 | AF_outer | mild_vs_severe |
| ANGPT1 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| ANGPT2 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| ASIC1 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| ASIC2 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| ASIC3 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| BDKRB1 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| BDKRB2 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| BDNF | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| CALCA | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| CALCB | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| CCL2 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| CXCL8 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| FGF2 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| FLT1 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| GAL | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| IL1B | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| IL6 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| KDR | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| KLK1 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| NGF | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| NGFR | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| NPY | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| NRP1 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| NRP2 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| NTF3 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| NTN1 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| NTN4 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| NTRK2 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| NTRK3 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| OPRD1 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| P2RX4 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| P2RX7 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| PDGFA | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| PDGFB | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| PENK | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| PLA2G2A | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| PTGES | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| PTGS2 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| ROBO1 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| ROBO2 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| SCN11A | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| SCN9A | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| SEMA3A | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| SEMA3E | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| SLIT2 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| SLIT3 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| TAC1 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| TEK | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| TNF | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| TRPV1 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| TRPV4 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| UNC5B | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| VEGFA | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| VEGFB | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| VIP | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| SEMA3F | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| VEGFC | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| DCC | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| NTF4 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| NTRK1 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| OPRM1 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| SCN10A | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| SLIT1 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |
| TRPA1 | 1.000 | Fibroblast_like | healthy_vs_degenerated_all |

3075 pain-relevant ligand-receptor interactions from CCC analysis. *[source: `results/communication/pain_interactions.tsv`]*

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

47 cell type x condition comparisons were skipped due to insufficient sample counts. *[source: `results/differential/skipped_comparisons.tsv`]*

### Result sensitivity across pipeline versions

Several results are sensitive to upstream methodological choices (integration method, annotation, cell sampling). These are documented here to flag areas requiring cautious interpretation. *[source: `docs/version_history.md`]*

- **Trajectory pseudotime-condition correlations** are sensitive to integration method and root cell choice. In v5 (CCA): 
  AF rho=+0.195; CEP rho=+0.073; NP rho=-0.088. Prior versions showed sign changes (e.g., CEP: -0.163 in v2, +0.135 in v3, +0.073 in v5), indicating these correlations are not robust to upstream choices.

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

- **Git commit:** `73388288c5d2696db93b3646193f7ca50d2ebf33` (branch: `main`)
- **Random seeds:** 42 (all stochastic operations)
- **Package versions:** pinned in `requirements.txt`, frozen in `requirements_frozen.txt`
- **Parameter choices:** documented in `analysis_plan.md`
- **Human checkpoint decisions:** recorded in `analysis_plan.md`
- **Data provenance:** GEO/CNGB accessions and download dates in `metadata/dataset_registry.tsv`
- **File checksums:** `metadata/file_checksums.json`
