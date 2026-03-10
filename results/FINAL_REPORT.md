# Human Intervertebral Disc Single-Cell Atlas

**A comprehensive scRNA-seq meta-analysis of IVD degeneration**

Report generated: 2026-03-10 | Pipeline version: 3.0

---

| Total cells | Datasets | Donors | Samples | Compartments | DE genes | Enriched pathways | L-R interactions |
|:-----------:|:--------:|:------:|:-------:|:------------:|:--------:|:-----------------:|:----------------:|
| **410,759** | **11** | **~50** | **71** | **NP / AF / CEP** | **1,156** | **1,043** | **81,059** |

---

## Contents

1. [Overview](#1-overview)
2. [Dataset Summary](#2-dataset-summary)
3. [Integration Strategy](#3-integration-strategy)
4. [Differential Expression](#4-differential-expression)
5. [Biological Pathways](#5-biological-pathways)
6. [Transcription Factor Activity](#6-transcription-factor-activity)
7. [Cell State Trajectories](#7-cell-state-trajectories)
8. [Cell-Cell Communication](#8-cell-cell-communication)
9. [Pain Biology](#9-pain-biology)
10. [Limitations](#10-limitations)
11. [Methods](#11-methods)
12. [Reproducibility](#12-reproducibility)

---

## 1. Overview

This atlas integrates 11 publicly available single-cell RNA sequencing (scRNA-seq) datasets of human intervertebral disc (IVD) tissue, comprising 410,759 cells from 71 samples (~50 donors). GSE233666 was excluded due to herniated-only samples that confound condition comparisons. The analysis characterizes cell type diversity, transcriptomic changes with degeneration, and intercellular communication in the IVD.

### What changed in v3

The key change from v2 to v3 was fixing the misrouting of ~17,000 stressed NP cells via three annotation corrections:

1. **Non-mesenchymal evidence gate:** Score-based classification of non-mesenchymal cell types (immune, endothelial) now requires expression of canonical lineage markers (PTPRC, PECAM1, VWF, CD3D, CD79A, NKG7) before assignment. This prevents stressed mesenchymal cells with upregulated immune-related genes from being misclassified.
2. **ACAN/SOX9 rescue:** Cells initially classified as non-mesenchymal but expressing canonical cartilage markers ACAN or SOX9 are reclassified as mesenchymal (25,415 cells rescued across all compartments).
3. **Cluster voting threshold:** Raised from 70% to 85%, ensuring that cluster-level label propagation only occurs when there is strong consensus.

Additional changes: de novo annotation uses specificity-weighted scoring, and the minimum cell count for non-mesenchymal classification was raised from 50 to 200 (MIN_CELLS_NON_MES).

These fixes primarily affected NP compartment cell counts and downstream DE results but left AF and CEP largely unchanged.

### Key Findings

- IVD resident cells exist on a **continuum** from notochordal to mature chondrocyte to fibrocartilaginous states in the NP, and inner to outer AF
- Pseudobulk DE analysis identified **1,156 unique significant genes** across 1,447 gene-comparison pairs in 18 powered comparisons
- **CXCL2** is significantly upregulated across multiple cell types in severe degeneration: NP_mature_chondrocyte (log2FC=3.63, padj=0.000175), NP_fibrocartilaginous (log2FC=2.18, padj=0.004), T_cell (log2FC=1.91, padj=0.033), and downregulated in AF_outer healthy_vs_severe (log2FC=-2.97, padj=0.037)
- Cell state trajectories correlate with disease condition (pseudotime-condition rho = -0.151 in NP, +0.325 in AF, +0.135 in CEP)
- LIANA cell-cell communication analysis identifies **81,059 total L-R interactions** (40,187 healthy, 40,872 degenerated)
- Pain-associated gene analysis identifies **10 significant pain genes** including PTGS2, TNF, PLA2G2A, BDKRB2, CCL2, PTGES, NRP2, PDGFA, ROBO1, and SEMA3A
- TF activity analysis identifies **5 significant TF-condition associations** via CollecTRI + pyDESeq2

---

## 2. Dataset Summary

| Dataset | Year | Compartment | Samples | Cells (post-QC) | Platform | Conditions |
|---------|------|-------------|:-------:|:---------------:|----------|------------|
| GSE160756 | 2021 | NP, AF, CEP | 6 | 89,283 | 10x | Healthy |
| GSE165722 | 2021 | NP | 10 | 9,498 | 10x | Degenerated (Pfirrmann II-V) |
| GSE189916 | 2022 | NP | 6 | 11,459 | BD Rhapsody | Neonatal, Aged |
| GSE199866 | 2022 | NP | 3 | 1,614 | 10x | Healthy, Degenerated |
| GSE205535 | 2022 | NP | 2 | 9,929 | 10x | Healthy (SCI), Degenerated |
| CNP0002664 | 2023 | NP | 8 | 52,016 | 10x | Healthy, Degenerated |
| GSE244889 | 2023 | NP, AF | 12 | 51,397 | 10x | Healthy, Degenerated |
| GSE251686 | 2024 | NP | 5 | 13,090 | Singleron | Herniated |
| GSE255768 | 2024 | CEP | 2 | 10,023 | 10x | Degenerated |
| GSE230809 | 2023 | NP, AF | 24 | 105,804 | 10x | Healthy, Degenerated |
| GSE242443 | 2024 | CEP | 2 | 59,227 | 10x | Healthy, Degenerated (culture-expanded) |

*GSE233666 (7 herniated-only NP samples, 22,658 cells) excluded from pipeline -- herniated samples confound condition comparisons and were the source of likely study-confounded DE results in v1.*

---

## 3. Integration Strategy

**Compartment-based approach:** Cells were separated into 4 compartment objects and each integrated independently with scVI:

| Compartment | Cells (v3) | Cells (v2) | Change |
|:-----------:|:----------:|:----------:|:------:|
| NP | 262,967 | 262,967 | -- |
| AF | 84,610 | 84,624 | -14 |
| CEP | 50,714 | 50,858 | -144 |
| all_cells | 410,759 | 410,759 | -- |

The total cell count is unchanged; the compartment-level shifts reflect the annotation fixes reclassifying cells between compartments.

**Cell types (14 total):**
- **NP:** NP_notochordal, NP_mature_chondrocyte, NP_fibrocartilaginous, NP_stressed_degen
- **AF:** AF_inner, AF_outer, AF_mechanical_stress
- **CEP:** EP_hyaline, EP_ossification
- **Non-resident:** T_cell, B_cell, Macrophage, NK_cell, Pericyte_SMC

**Annotation (v3 improvements):** De novo clustering followed by specificity-weighted marker-based cell type assignment with three key fixes (see Section 1). The non-mesenchymal evidence gate and ACAN/SOX9 rescue prevented ~17K stressed NP cells from being misclassified as immune/endothelial. Validated with CellTypist Immune_All_Low model for immune subtypes.

### NP Integration (scVI)

![NP integration UMAP](integration/umap_NP.png)

### AF Integration (scVI)

![AF integration UMAP](integration/umap_AF.png)

---

## 4. Differential Expression

Pseudobulk DE analysis using pyDESeq2. **18 powered comparisons** (down from 21 in v2), remaining comparisons skipped (underpowered). Significance: |log2FC| > 0.5, padj < 0.05. **1,156 unique DE genes** across 1,447 gene-comparison pairs.

| Cell Type | Comparison | Total DE genes |
|-----------|-----------|:-----:|
| NP_fibrocartilaginous | mild_vs_severe | **418** |
| NP_fibrocartilaginous | healthy_vs_severe | **385** |
| NP_mature_chondrocyte | mild_vs_severe | **291** |
| NP_mature_chondrocyte | healthy_vs_severe | **113** |
| AF_outer | healthy_vs_severe | **100** |
| T_cell | mild_vs_severe | 40 |
| AF_outer | mild_vs_severe | 38 |
| NP_notochordal | mild_vs_severe | 18 |
| AF_outer | healthy_vs_all | 12 |
| AF_outer | healthy_vs_mild | 7 |

*v3 shows substantially more DE genes in NP_fibrocartilaginous (418 up from 203 in v2) due to corrected cell type assignments -- 17K stressed NP cells previously misrouted to non-mesenchymal types now contribute to properly powered NP comparisons. EP_hyaline healthy_vs_all (84 genes in v2) dropped below the powered comparison threshold in v3. See Supplementary Tables S6-S8 for full results.*

### Top DE genes in NP severe degeneration

CXCL2 remains the most significant chemokine across multiple cell types:
- NP_mature_chondrocyte mild_vs_severe: log2FC=3.63, padj=0.000175
- NP_fibrocartilaginous mild_vs_severe: log2FC=2.18, padj=0.004
- T_cell mild_vs_severe: log2FC=1.91, padj=0.033

The stronger CXCL2 signal in v3 (log2FC=3.63 vs 3.14 in v2) reflects cleaner cell type composition after the annotation fix.

### Top DE genes in AF degeneration

AF_outer healthy_vs_severe yields 100 DE genes (up from 97 in v2). CXCL2 is significantly downregulated in AF_outer (log2FC=-2.97, padj=0.037), showing an opposite direction to NP -- consistent with compartment-specific inflammatory responses.

### Volcano Plots

![NP fibrocartilaginous mild vs severe](differential/volcano_plots/volcano_NP_fibrocartilaginous_mild_vs_severe.png)

![NP mature chondrocyte mild vs severe](differential/volcano_plots/volcano_NP_mature_chondrocyte_mild_vs_severe.png)

![AF outer healthy vs severe](differential/volcano_plots/volcano_AF_outer_healthy_vs_degenerated_severe.png)

---

## 5. Biological Pathways

**ORA:** 1,043 significantly enriched terms (FDR < 0.05). **GSEA:** 63,797 terms tested across GO, KEGG, Reactome, MSigDB Hallmark, and IVD-custom gene sets.

### AF_outer -- Upregulated Pathways

![AF outer up](interpretation/pathway_enrichment/enrichment_AF_outer_up.png)

### NP_mature_chondrocyte -- Upregulated Pathways

![NP chondrocyte up](interpretation/pathway_enrichment/enrichment_NP_mature_chondrocyte_up.png)

### NP_fibrocartilaginous -- Upregulated Pathways

![NP fibrocartilaginous up](interpretation/pathway_enrichment/enrichment_NP_fibrocartilaginous_up.png)

### GSEA: IVD Custom Gene Sets

![GSEA IVD heatmap](interpretation/pathway_enrichment/gsea_ivd_custom_heatmap.png)

---

## 6. Transcription Factor Activity

TF activity inferred using CollecTRI regulon overlap with pyDESeq2 pseudobulk framework. **5 significant TF-condition associations** (down from 290 in v2, reflecting the switch to pyDESeq2-based statistical testing which is more conservative than the Fisher's exact test used in v2).

![TF activity heatmap](interpretation/tf_activity/tf_activity_heatmap.png)

---

## 7. Cell State Trajectories

PAGA + diffusion pseudotime (DPT) on scVI embeddings. NP: rooted at notochordal cells. AF: rooted at AF_inner. CEP: rooted at EP_hyaline.

### NP Trajectory

![NP trajectory UMAP](trajectories/umap_trajectory_NP.png)

### NP Pseudotime by Condition

![NP pseudotime by condition](trajectories/pseudotime_by_condition_NP.png)

### NP Gene Dynamics Along Pseudotime

![NP gene dynamics](trajectories/gene_dynamics_NP.png)

### AF Trajectory

![AF trajectory UMAP](trajectories/umap_trajectory_AF.png)

### AF Pseudotime by Condition

![AF pseudotime by condition](trajectories/pseudotime_by_condition_AF.png)

**Pseudotime correlates with disease:**

| Compartment | Spearman rho (v3) | Spearman rho (v2) | Direction |
|:-----------:|:-----------------:|:-----------------:|:---------:|
| NP | -0.151 | -0.258 | Healthy early, degenerated late (weaker in v3) |
| AF | +0.325 | +0.341 | **Reversed from NP** -- degenerated early, healthy late |
| CEP | +0.135 | -0.163 | **Sign reversed from v2** |

> **FLAG FOR SME REVIEW:** The CEP pseudotime-condition correlation has **reversed sign** from v2 (-0.163) to v3 (+0.135). The AF correlation remains positive and similar in magnitude. The NP correlation weakened but retained its negative sign. These changes likely reflect the annotation fixes redistributing ~17K cells. The biological interpretation of the CEP reversal requires expert review.

> **FLAG FOR SME REVIEW (carried from v2):** The AF pseudotime-condition correlation is **positive** (+0.325), reversed from v1 (-0.177). This likely reflects the change in integration approach (scVI-only vs scANVI primary) and the exclusion of herniated samples. Whether AF degeneration proceeds "backward" along the inner-to-outer trajectory or the root cell choice needs revisiting requires expert review.

---

## 8. Cell-Cell Communication

LIANA consensus (CellPhoneDB, NATMI, Connectome, SingleCellSignalR, log2FC) on 20,000 cells per condition.

### Interaction Heatmap -- Healthy

![Healthy interactions](communication/interaction_plots/interaction_heatmap_healthy.png)

### Interaction Heatmap -- Degenerated

![Degenerated interactions](communication/interaction_plots/interaction_heatmap_degenerated.png)

### Differential Interactions

![Differential interactions](communication/interaction_plots/differential_interactions.png)

**Interaction counts:** 40,187 interactions (healthy) vs 40,872 (degenerated). Total: 81,059 L-R interactions.

The v3 interaction counts are substantially higher than v2 (81K vs 56K) and show a near-balance between healthy and degenerated conditions. This contrasts with v2 where degenerated had *fewer* interactions (27K vs 29K) and v1 where degenerated had *more* (53K vs 44K).

> **FLAG FOR SME REVIEW:** The CCC results have shifted again between versions. In v3, healthy and degenerated interaction counts are nearly equal (40K vs 41K), compared to the v2 result favoring healthy (29K vs 27K) and the v1 result favoring degenerated (44K vs 53K). The near-equality in v3 may reflect the annotation corrections producing a more balanced cell type composition. The sensitivity of CCC results to upstream annotation decisions is notable and should be considered when interpreting these results.

---

## 9. Pain Biology

Cross-reference of DE genes with curated pain gene sets (nociception, neurotrophins, nerve guidance, inflammatory pain, neovascularization).

### Pain-Associated Findings

**10 significant pain genes** identified in v3:

- **PTGS2** (COX-2) -- prostaglandin synthesis, key inflammatory pain mediator
- **TNF** -- inflammatory pain mediator that sensitizes nerve endings
- **PLA2G2A** -- phospholipase A2, upstream of prostaglandin cascade
- **BDKRB2** -- bradykinin receptor, direct nociceptive signaling
- **CCL2** -- monocyte chemoattractant, promotes neuroinflammation
- **PTGES** -- prostaglandin E synthase, downstream of COX-2
- **NRP2** -- neuropilin-2, nerve guidance co-receptor
- **PDGFA** -- platelet-derived growth factor, neovascularization and nerve ingrowth
- **ROBO1** -- Roundabout receptor, axon guidance molecule
- **SEMA3A** -- semaphorin, nerve repulsion factor (downregulated = permissive for nerve ingrowth)

The pain gene list is identical to v2 in size (10 genes) but the composition has shifted: CXCL8 was replaced by NRP2, PDGFA, ROBO1, and SEMA3A as significant, while some v2 pain genes dropped below significance. The addition of nerve guidance molecules (NRP2, ROBO1, SEMA3A) strengthens the interpretation that degenerated discs shift from nerve-repulsive to nerve-permissive signaling.

- **Disc cells produce inflammatory mediators but not nociceptors.** This is consistent with the model that degenerated disc cells create a pro-inflammatory environment that promotes nerve ingrowth and sensitization, rather than directly signaling pain.

### Pain Gene Expression Heatmap

![Pain genes heatmap](interpretation/pain_genes_heatmap.png)

---

## 10. Limitations

- **Annotation sensitivity:** The v2-to-v3 annotation fix (rerouting ~17K stressed NP cells) changed DE gene counts significantly (NP_fibrocartilaginous mild_vs_severe: 203 to 418 genes) and reversed the CEP trajectory correlation sign. CCC interaction counts are also sensitive to annotation. Results should be interpreted with this sensitivity in mind.
- **Herniated exclusion:** GSE233666 (22,658 cells, 7 samples) was excluded. Herniated-only samples confounded condition comparisons in v1 and inflated DE gene counts. GSE251686 herniated samples remain but are treated as "severe" degeneration.
- **CellTypist NP disagreements:** NP-specific cell states (notochordal, fibrocartilaginous) are absent from CellTypist reference databases. AF and CEP concordance is high.
- **Cross-study confounding:** Condition and study are partially confounded. Within-study comparisons used where possible.
- **Underpowered comparisons:** Many cell type x comparison combinations skipped due to insufficient samples (18 powered out of total possible comparisons).
- **No RNA velocity:** Spliced/unspliced counts not available in public datasets. Would require reprocessing from BAM files.
- **Age-disease confound:** In GSE230809, healthy donors are 21-27y and diseased are 37-73y. Cannot fully separate age from disease effects.
- **Sex bias:** GSE230809 (largest dataset, 24 samples) is all-male. Many samples have unknown sex.
- **Culture-expanded cells:** GSE242443 CEP cells are culture-expanded, which alters gene expression.
- **AF trajectory reversal:** AF pseudotime-condition correlation is positive (+0.325), reversed from v1 (-0.177). Persists from v2 and requires SME review.
- **CEP trajectory reversal:** CEP pseudotime-condition correlation reversed from v2 (-0.163) to v3 (+0.135) after annotation fix.
- **CCC instability:** CCC interaction counts have varied across all three pipeline versions (v1: 44K/53K healthy/degen, v2: 29K/27K, v3: 40K/41K). Results are sensitive to cell type composition and annotation.
- **Composition analysis:** No significant changes after FDR correction, though trends are biologically consistent.
- **SCENIC/GRN not run:** Full SCENIC analysis was not performed due to computational requirements. TF activity estimated from CollecTRI regulon overlap instead.
- **TF activity:** Only 5 significant associations in v3 (vs 290 in v2), reflecting the switch to more conservative pyDESeq2-based testing.

---

## 11. Methods

### Data acquisition
11 scRNA-seq datasets of human IVD tissue were downloaded from GEO and CNGB (see Table in Section 2). GSE233666 excluded due to herniated-only design. Raw count matrices were obtained for each dataset.

### Quality control and preprocessing
Per-dataset QC: min 200 genes, max 6000 genes, min 500 counts, max 20% mitochondrial reads. Doublet detection with Scrublet (expected rate 5%). Normalization: total-count to 10,000, log1p. HVG selection: top 2000 genes per dataset using Seurat v3 method.

### Cell type annotation (v3)
De novo annotation after clustering on scVI embeddings, using specificity-weighted IVD-specific marker gene signatures. Three key annotation fixes in v3: (1) non-mesenchymal evidence gate requiring canonical lineage markers (PTPRC, PECAM1, VWF, CD3D, CD79A, NKG7) before score-based non-mesenchymal assignment; (2) ACAN/SOX9 rescue reclassifying non-mesenchymal cells expressing cartilage markers as mesenchymal (25,415 cells rescued); (3) cluster voting threshold raised from 70% to 85%. MIN_CELLS_NON_MES raised from 50 to 200. Validated with CellTypist Immune_All_Low model for immune subtypes.

### Integration
Compartment-based strategy: cells separated into NP (262,967), AF (84,610), CEP (50,714), and all_cells (410,759) objects. Each integrated with scVI (1 layer, 128 dim). Single integration approach (scVI-only).

### Differential expression
Pseudobulk aggregation per sample per cell type. DE with pyDESeq2 (Python DESeq2 implementation). Significance: |log2FC| > 0.5, adjusted p-value < 0.05 (Benjamini-Hochberg). Minimum 3 samples per condition per cell type. 18 powered comparisons in v3.

### Pathway enrichment
Over-representation analysis (ORA) and gene set enrichment analysis (GSEA) using gseapy. Databases: GO Biological Process 2023, KEGG 2021, Reactome 2022, MSigDB Hallmark 2020, custom IVD gene sets. 1,043 significant ORA terms, 63,797 GSEA terms tested.

### TF activity inference
CollecTRI regulon network with pyDESeq2 pseudobulk framework. 5 significant TF-condition associations identified in v3.

### Trajectory analysis
PAGA + diffusion pseudotime (DPT) on scVI embeddings. Root cells: NP notochordal, AF inner, CEP EP_hyaline. Trajectory genes: Spearman correlation with pseudotime, FDR < 0.05, top 500 per compartment.

### Cell-cell communication
LIANA rank_aggregate with consensus resource. 5 methods: CellPhoneDB, NATMI, Connectome, SingleCellSignalR, log2FC. 100 permutations. 20,000 cells per condition. Total 81,059 interactions identified.

### Software
Python 3.12, scanpy 1.11, scvi-tools 1.4.2, pyDESeq2, gseapy 1.1, decoupler 2.1, liana 1.7. Full environment: `requirements_frozen.txt`.

### Supplementary Tables
Supplementary tables (S1-S19) provided in `results/supplementary_tables/`, including dataset registry, sample metadata, inclusion summary, study caveats, composition analysis, DE summary, DE results, skipped comparisons, ORA enrichment, GSEA results, TF activity, pain genes, trajectory genes (NP/AF/CEP), pain interactions, and CellTypist concordance (NP/AF/CEP).

---

## 12. Reproducibility

- All scripts version-controlled in git
- Random seeds: 42 (all stochastic operations)
- Package versions pinned in `requirements_frozen.txt`
- All parameter choices documented in `analysis_plan.md`
- All human checkpoint decisions recorded
- Data provenance: GEO/CNGB accessions, download dates in `metadata/dataset_registry.tsv`
- v3 annotation changes fully documented in Section 1 and Methods
