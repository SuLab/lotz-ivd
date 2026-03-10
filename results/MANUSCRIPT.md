# Inflammatory Signatures and Cell State Continua in Human Intervertebral Disc Degeneration: An 11-Dataset Single-Cell Transcriptomic Meta-Analysis

**Draft Manuscript (v2 pipeline)**
**Analysis Date: March 2026**

> **Companion Analysis:** An independent analysis of 7 of these 11 datasets (173,628 cells, 29 donors) was performed using a separate pipeline (Harmony integration, R-based DESeq2, LIANA CCC). That analysis and its full report are available at [`phylo_analysis/report_IVD_scRNAseq_analysis.md`](../phylo_analysis/report_IVD_scRNAseq_analysis.md), with a corresponding draft manuscript at [`phylo_analysis/draft_manuscript.md`](../phylo_analysis/draft_manuscript.md). The phylo analysis is referenced throughout this document where its findings converge with or diverge from the present 11-dataset analysis.

---

## Table of Contents

1. [Abstract](#1-abstract)
2. [Introduction](#2-introduction)
3. [Study Design and Datasets](#3-study-design-and-datasets)
4. [Methods](#4-methods)
5. [Results](#5-results)
   - 5.1 Integrated Cell Atlas
   - 5.2 Differential Gene Expression
   - 5.3 Pathway Enrichment
   - 5.4 Transcription Factor Activity
   - 5.5 Cell State Trajectories
   - 5.6 Cell-Cell Communication
   - 5.7 Pain Biology
6. [Biological Interpretation and Mechanistic Model](#6-biological-interpretation-and-mechanistic-model)
7. [Therapeutic Targets](#7-therapeutic-targets)
8. [Novel and Discordant Findings](#8-novel-and-discordant-findings)
9. [Limitations](#9-limitations)
10. [Conclusion](#10-conclusion)
11. [References](#11-references)

---

## 1. Abstract

Intervertebral disc (IVD) degeneration is the primary structural cause of chronic low back pain, affecting over 600 million people worldwide (GBD 2021 Low Back Pain Collaborators, 2023). To comprehensively map its cellular and molecular landscape, we integrated 11 publicly available human scRNA-seq datasets comprising 410,759 cells from 71 samples (~50 donors) across nucleus pulposus (NP), annulus fibrosus (AF), and cartilage endplate (CEP) compartments. Using scVI integration with compartment-specific objects (NP, AF, CEP, and a combined all-cells object), followed by de novo clustering and marker-based annotation, we identified cell types existing on a continuum from notochordal to mature chondrocyte to stressed/degenerative states. Pseudobulk differential expression with pyDESeq2 identified 949 unique significant genes across 1,231 gene-by-comparison pairs in 21 powered comparisons, revealing an inflammatory/catabolic signature in severe NP degeneration. Among the top upregulated genes in NP_mature_chondrocyte mild_vs_severe were HHEX (log2FC=+2.78), CXCL2 (+3.14), ICAM1 (+1.67), and UGCG (+1.06). Pathway enrichment confirmed inflammatory and chemokine-mediated signaling among upregulated programs in NP cells, alongside heat shock protein activation and mitochondrial dysfunction in AF cells. Transcription factor analysis identified 290 significant TF-condition associations. PAGA/diffusion pseudotime trajectory analysis demonstrated that pseudotime correlates with disease condition in NP (rho=-0.258, p<10^-100) and CEP (rho=-0.163), while AF showed an unexpected positive correlation (rho=+0.341) requiring further investigation. Cell-cell communication analysis (LIANA) revealed that degenerated tissue has fewer interactions than healthy tissue (27,011 vs 28,878), contrasting with some prior reports. Pain gene analysis identified 10 significant pain-relevant genes across comparisons. These findings define an inflammatory mechanism of IVD degeneration centered on NF-kB-driven chemokine and cytokine activation, with compartment-specific stress responses, and identify TNF/NF-kB inhibition, chemokine modulation, and heat shock protein targeting as candidate therapeutic strategies.

---

## 2. Introduction

### 2.1 The Clinical Problem

Low back pain (LBP) is the leading cause of years lived with disability worldwide, affecting approximately 619 million people and imposing annual costs exceeding $100 billion in the United States alone (GBD 2021 Low Back Pain Collaborators, 2023; Dieleman et al., 2020). Approximately 40% of symptomatic LBP is attributable to IVD degeneration (Wang et al., 2023a). Despite decades of research, no disease-modifying therapy exists; current treatments are limited to symptomatic relief (analgesics, physical therapy) or surgical intervention (discectomy, fusion) for refractory cases.

### 2.2 IVD Structure and Biology

The IVD is a fibrocartilaginous structure comprising three compartments: the nucleus pulposus (NP), a highly hydrated gel-like core rich in aggrecan and type II collagen that absorbs compressive loads; the annulus fibrosus (AF), a tough outer ring of concentric collagen I-rich lamellae providing tensile strength; and the cartilage endplate (CEP), thin hyaline cartilage layers that serve as the primary route for nutrient diffusion into the avascular NP (Oichi et al., 2020). The NP is the largest avascular structure in the human body, forcing its resident cells to operate under near-anoxic conditions via anaerobic glycolysis (Oichi et al., 2020).

### 2.3 Pathomechanisms of Degeneration

IVD degeneration is characterized by progressive loss of NP hydration through aggrecan degradation by ADAMTS4/5 and MMPs (Liang et al., 2022), inflammatory activation driven by TNF-alpha, IL-1beta, and NF-kB signaling (Risbud and Shapiro, 2014; Xia et al., 2024), cellular senescence and apoptosis (Song et al., 2023a), oxidative stress from mitochondrial dysfunction (Song et al., 2023b; Wang et al., 2023a), and fibrocartilaginous replacement of the NP by type I collagen-producing cells (Antoniou et al., 1996). In advanced degeneration, nerve fibers and blood vessels invade the normally avascular NP through AF fissures, contributing to discogenic pain (Freemont et al., 2002).

### 2.4 Rationale for Single-Cell Meta-Analysis

Prior single-cell studies of the IVD have been limited by small sample sizes (2-7 donors), single datasets, or focus on a single compartment (Gan et al., 2021; Fernandes et al., 2020; Li et al., 2022a). By integrating 11 datasets spanning ~50 donors, we aimed to create a comprehensive single-cell atlas of human IVD degeneration and achieve sufficient statistical power for pseudobulk differential expression analysis, the gold standard for scRNA-seq DE that avoids the inflated false positive rates of naive single-cell approaches (Squair et al., 2021; Zimmerman et al., 2021).

A critical methodological consideration is the distinction between resident disc cells (NP, AF) and non-resident cells (immune, endothelial). IVD resident cells exist on a phenotypic continuum — from notochordal to mature chondrocyte to stressed/degenerative states — that can be erased by aggressive batch correction (Gan et al., 2021). Our compartment-specific integration strategy addresses this by building separate scVI models for NP, AF, and CEP compartments, followed by de novo clustering and marker-based annotation that respects biological heterogeneity.

---

## 3. Study Design and Datasets

### 3.1 Dataset Selection

Eleven publicly available scRNA-seq datasets of human IVD tissue were identified from GEO and CNGB. Selection criteria included: (1) human IVD tissue, (2) single-cell resolution (not single-nucleus), (3) raw count matrices available. GSE233666 was excluded as it contained only herniated samples with no healthy or graded-degeneration controls.

**Table 1. Datasets included in the integrated atlas.**

| Accession | Year | Compartment | Samples | Cells (post-QC) | Platform | Conditions |
|-----------|------|-------------|:-------:|:---------------:|----------|------------|
| GSE160756 | 2021 | NP, AF, CEP | 6 | 89,283 | 10x | Healthy |
| GSE165722 | 2021 | NP | 10 | 9,498 | 10x | Degenerated (Pfirrmann II-V) |
| GSE189916 | 2022 | NP | 6 | 11,459 | BD Rhapsody | Neonatal, Aged |
| GSE199866 | 2022 | NP | 3 | 1,614 | 10x | Healthy, Degenerated |
| GSE205535 | 2022 | NP | 2 | 9,929 | 10x | Healthy, Degenerated |
| CNP0002664 | 2023 | NP | 8 | 52,016 | 10x | Healthy, Degenerated |
| GSE244889 | 2023 | NP, AF | 12 | 51,397 | 10x | Healthy, Degenerated |
| GSE251686 | 2024 | NP | 5 | 13,090 | Singleron | Herniated |
| GSE255768 | 2024 | CEP | 2 | 10,023 | 10x | Degenerated |
| GSE230809 | 2023 | NP, AF | 24 | 105,804 | 10x | Healthy, Degenerated |
| GSE242443 | 2024 | CEP | 2 | 59,227 | 10x | Healthy, Degenerated (culture-expanded) |

**Total:** 410,759 cells from 71 samples (~50 donors).

### 3.2 Condition Harmonization

Degeneration severity was harmonized across datasets using Pfirrmann grading where available: **healthy** (Pfirrmann I-II), **mild** (Pfirrmann II-III), **severe** (Pfirrmann IV-V). Herniated samples were retained but herniated comparisons were excluded from DE analysis due to single-study confounding (only GSE251686 contributes herniated samples after GSE233666 exclusion, making herniated fully confounded with study).

---

## 4. Methods

### 4.1 Quality Control and Preprocessing

Per-dataset QC applied fixed thresholds: minimum 200 genes, maximum 6,000 genes, minimum 500 counts, maximum 20% mitochondrial reads. Doublet detection used Scrublet (Wolock et al., 2019) at 5% expected rate. Normalization: total-count to 10,000, log1p transformation. HVG selection: top 2,000 genes per dataset (Seurat v3 method).

### 4.2 Cell Type Annotation

Two-stage annotation was used:

1. **Module 04 — Binary classification:** Cells were classified as mesenchymal or non-mesenchymal using marker-based scoring with IVD-specific gene signatures curated from published atlases (Gan et al., 2021; Risbud and Shapiro, 2014).

2. **Module 05 — De novo annotation:** After scVI integration and clustering, cell types were assigned by marker gene scoring. For non-mesenchymal clusters, CellTypist (Immune_All_Low model; Dominguez Conde et al., 2022) was used for validation. The resulting `cell_type` labels include: NP_mature_chondrocyte, NP_fibrocartilaginous, NP_notochordal, NP_stressed_degenerative, AF_inner, AF_outer, EP_hyaline, T_cell, B_cell, Macrophage, Endothelial_cells, Pericyte_SMC, and NK_cell.

### 4.3 Integration

All cells were integrated using scVI (Lopez et al., 2018; 1 layer, 128 dimensions) with four compartment-specific objects:
- **NP:** 262,967 cells
- **AF:** 84,624 cells
- **CEP:** 50,858 cells
- **all_cells:** 410,759 cells (combined)

This replaces the v1 two-tier structure (tier1_nonresident + tier2_resident_NP/AF). scVI was chosen for its strong batch correction while preserving biological variation, particularly for the cell state continua present in IVD resident cells.

### 4.4 Pseudobulk Differential Expression

Cells were aggregated into pseudobulk samples per donor per cell type. DE analysis used pyDESeq2 (Love et al., 2014) with Benjamini-Hochberg correction. Significance thresholds: |log2FC| > 0.5 and adjusted p-value < 0.05. Minimum 3 samples per condition per cell type were required for inclusion. LFC shrinkage was applied by default to constrain estimates to biologically plausible ranges. Herniated comparisons were excluded entirely due to single-study confounding.

### 4.5 Pathway Enrichment

Over-representation analysis (ORA) and gene set enrichment analysis (GSEA) were performed using gseapy (Fang et al., 2023) against GO Biological Process 2023, KEGG 2021, Reactome 2022, MSigDB Hallmark 2020 (Liberzon et al., 2015), and custom IVD-specific gene sets. For GSEA, genes were ranked by sign(log2FC) x -log10(p-value).

### 4.6 Transcription Factor Activity

TF activity was inferred using CollecTRI regulon networks (Garcia-Alonso et al., 2019) containing 42,990 TF-target interactions across 1,185 TFs. For each TF, enrichment of its targets among DE genes was tested using Fisher's exact test, with concordance scoring to account for activation vs. repression direction.

### 4.7 Trajectory Analysis

PAGA + diffusion pseudotime (DPT; Haghverdi et al., 2016) was computed on scVI embeddings. 50,000 cells were sampled per compartment. Root cells: NP notochordal for NP compartment, AF inner for AF compartment, EP hyaline for CEP compartment. Trajectory-associated genes were identified by Spearman correlation with pseudotime (FDR < 0.05, top 500).

### 4.8 Cell-Cell Communication

LIANA rank_aggregate (Dimitrov et al., 2022) was applied with five consensus methods (CellPhoneDB, NATMI, Connectome, SingleCellSignalR, log2FC) and 100 permutations. Analysis was run separately on healthy and degenerated subsets (20,000 cells each) using per-dataset processed files to avoid integration artifacts.

---

## 5. Results

### 5.1 Integrated Cell Atlas

The atlas comprises 410,759 cells organized into distinct populations across four compartment-specific objects. NP cells (262,967) segregate into four major states: NP_notochordal (expressing KRT8, KRT18, T/TBXT), NP_mature_chondrocyte (ACAN, COL2A1, SOX9), NP_fibrocartilaginous (COL1A1, transitional phenotype), and NP_stressed_degenerative (HSPA5, DDIT3, stress markers). AF cells (84,624) separate into AF_inner (transitional, cartilage-like) and AF_outer (COL1A1, COL1A2, fibrous). CEP cells (50,858) are annotated as EP_hyaline. Non-resident populations include T_cell, B_cell, Macrophage, Endothelial_cells, Pericyte_SMC, and NK_cell.

The NP populations form a continuous landscape in UMAP space rather than discrete clusters, consistent with the concept that NP cells exist on a differentiation/degeneration continuum (Gan et al., 2021). scVI integration preserves this continuum while correcting batch effects across the 11 datasets. The companion phylo analysis independently identified a similar set of NP cell states — including canonical, stress-response, UPR/degenerative, and metallothionein-high populations — using Harmony integration on a 7-dataset subset ([`phylo_analysis/report_IVD_scRNAseq_analysis.md`, Section 4](../phylo_analysis/report_IVD_scRNAseq_analysis.md#4-cell-type-annotation)).

**CellTypist validation for non-mesenchymal clusters:** CellTypist concordance varied by compartment. In NP, 5 of 13 non-mesenchymal clusters were concordant with CellTypist predictions (8 discordant). In AF, 1 cluster was discordant. In CEP, 3 clusters were discordant. The NP discordance is partly due to misrouted mesenchymal cells (~17K cells), scoring formula issues, and the lack of a pericyte category in the CellTypist model. This annotation quality issue is flagged as a limitation (see Section 9).

**Figure 1. NP integration — scVI UMAP by cell type.**
![NP integration UMAP](integration/umap_tier2_NP_by_approach.png)

**Figure 2. AF integration — scVI UMAP by cell type.**
![AF integration UMAP](integration/umap_tier2_AF_by_approach.png)

### 5.2 Differential Gene Expression

Pseudobulk DE identified **949 unique significant genes** across **1,231 gene-by-comparison pairs** in **21 powered comparisons** out of 74 tested (53 skipped due to insufficient samples per condition; Table 2). Notably, the CEP compartment is now powered for DE analysis, and herniated comparisons were excluded entirely due to single-study confounding.

**Table 2. Powered DE comparisons and significant genes.**

| Cell Type | Comparison | Up | Down | Total |
|-----------|-----------|:---:|:----:|:-----:|
| NP_mature_chondrocyte | mild_vs_severe | 183 | 132 | 315 |
| NP_fibrocartilaginous | mild_vs_severe | 99 | 104 | 203 |
| NP_mature_chondrocyte | healthy_vs_degenerated_severe | 94 | 78 | 172 |
| NP_fibrocartilaginous | healthy_vs_degenerated_severe | 83 | 44 | 127 |
| AF_outer | healthy_vs_degenerated_severe | 49 | 48 | 97 |
| EP_hyaline | healthy_vs_degenerated_all | 0 | 84 | 84 |
| T_cell | mild_vs_severe | 35 | 13 | 48 |
| AF_inner | healthy_vs_degenerated_severe | 33 | 14 | 47 |
| AF_outer | mild_vs_severe | 27 | 13 | 40 |
| AF_inner | healthy_vs_degenerated_all | 26 | 11 | 37 |
| NP_fibrocartilaginous | healthy_vs_degenerated_all | 18 | 5 | 23 |
| AF_outer | healthy_vs_degenerated_all | 8 | 4 | 12 |
| AF_outer | healthy_vs_degenerated_mild | 3 | 8 | 11 |
| NP_mature_chondrocyte | healthy_vs_degenerated_all | 2 | 3 | 5 |
| Macrophage | mild_vs_severe | 0 | 5 | 5 |
| B_cell | mild_vs_severe | 1 | 2 | 3 |
| NP_fibrocartilaginous | healthy_vs_degenerated_mild | 0 | 1 | 1 |
| unassigned | mild_vs_severe | 1 | 0 | 1 |

**Key finding: Inflammatory signature in NP severe degeneration.** The mild_vs_severe comparison in NP_mature_chondrocyte — which is more robust against cross-study confounding than healthy_vs_severe — reveals an inflammatory/catabolic signature, though more moderate than v1 estimates:

- **UGCG** (log2FC=+1.06, padj=6.3x10^-6): UDP-glucose ceramide glucosyltransferase, involved in glycosphingolipid metabolism
- **ICAM1** (log2FC=+1.67, padj=1.9x10^-5): intercellular adhesion molecule, a key mediator of immune cell recruitment
- **HHEX** (log2FC=+2.78, padj=1.9x10^-5): hematopoietically-expressed homeobox transcription factor
- **CXCL2** (log2FC=+3.14, padj=0.005): GRO-beta, inflammatory chemokine and neutrophil chemoattractant
- **TNC** (log2FC=+2.84, padj=0.017): tenascin-C, an ECM glycoprotein induced by inflammation and tissue damage
- **IL32** (log2FC=+2.57, padj=0.014): interleukin-32, a pro-inflammatory cytokine
- **COL12A1** (log2FC=+2.62, padj=0.026): type XII collagen, associated with fibrocartilaginous remodeling

Notably, several genes that were prominent in v1 did not reach significance in this reanalysis: CXCL1 (no significant padj), CXCL3 (padj=0.099), TNF (padj=0.22), and MDK (not significant). CEMIP was borderline (padj=0.055). The CXC chemokine signal is present but attenuated compared to v1, with only CXCL2 reaching significance among the GRO chemokines.

**CXCL2 remains a consistent signal:** Despite the attenuation of the broader CXC chemokine triad, CXCL2 (log2FC=+3.14, padj=0.005) remains robustly significant, supporting continued involvement of chemokine-mediated neutrophil recruitment in severe degeneration.

**NP_fibrocartilaginous — a new cell type with substantial DE signal.** NP_fibrocartilaginous cells, newly defined in the v2 de novo annotation, show 203 DE genes in mild_vs_severe and 127 in healthy_vs_degenerated_severe, making them the second most transcriptionally responsive NP cell type. This population likely captures transitional cells between the mature chondrocyte and fibrous phenotypes.

**CEP now powered for DE analysis.** EP_hyaline shows 84 downregulated genes in healthy_vs_degenerated_all (0 upregulated), suggesting a pattern of transcriptional silencing rather than activation in CEP degeneration. This is a notable new finding enabled by the compartment-specific integration approach.

**AF degeneration signature.** AF_outer in the healthy_vs_degenerated_severe comparison showed 97 DE genes (49 up, 48 down), a more balanced pattern than the v1 analysis.

**Figure 4. Volcano plot — NP mature chondrocyte, mild vs. severe degeneration.**
![NP mild vs severe volcano](differential/volcano_plots/volcano_NP_mature_chondrocyte_mild_vs_severe.png)

**Figure 5. Volcano plot — AF outer, healthy vs. severe degeneration.**
![AF healthy vs severe volcano](differential/volcano_plots/volcano_AF_outer_healthy_vs_degenerated_severe.png)

### 5.3 Pathway Enrichment

ORA identified **1,577 significantly enriched terms** (FDR < 0.05) and GSEA identified **1,576 significant terms** across GO, KEGG, Reactome, and MSigDB Hallmark databases.

**NP_mature_chondrocyte (mild_vs_severe, upregulated):** The dominant enriched pathways include:
- Cellular response to lipopolysaccharide
- Chemokine-mediated signaling pathway
- Neutrophil chemotaxis
- Inflammatory response
- Granulocyte chemotaxis

These pathways remain enriched even with the attenuated CXC chemokine signal, driven by CXCL2, ICAM1, and other inflammatory mediators.

**NP_mature_chondrocyte (healthy_vs_degenerated_severe, upregulated):** Cell cycle pathways dominate:
- Mitotic sister chromatid segregation
- Mitotic spindle organization

This may indicate compensatory proliferation of surviving chondrocytes, consistent with the "cluster formation" phenomenon observed histologically in degenerated discs (Johnson et al., 2001).

**Figure 6. Pathway enrichment — NP mature chondrocyte, upregulated in severe degeneration.**
![NP chondrocyte up pathways](interpretation/pathway_enrichment/enrichment_NP_mature_chondrocyte_up.png)

**Figure 7. Pathway enrichment — AF outer, upregulated in severe degeneration.**
![AF outer up pathways](interpretation/pathway_enrichment/enrichment_AF_outer_up.png)

**Figure 8. GSEA heatmap — IVD-specific custom gene sets across cell types.**
![GSEA IVD heatmap](interpretation/pathway_enrichment/gsea_ivd_custom_heatmap.png)

**AF_inner (mild_vs_severe):**
- **Upregulated:** Cellular response to heat, response to unfolded protein, TNF-mediated signaling regulation, granulocyte chemotaxis
- **Downregulated:** Oxidative phosphorylation, aerobic electron transport chain, mitochondrial ATP synthesis

The simultaneous heat shock protein upregulation and mitochondrial dysfunction in AF cells is a notable observation. It suggests that AF cells are experiencing proteotoxic stress (driving HSP induction) concurrent with metabolic failure (reduced oxidative phosphorylation), a combination that may represent an energy crisis limiting the ability of AF cells to maintain ECM homeostasis.

**Important negative finding:** Neither Wnt signaling, Notch signaling, nor cellular senescence pathways reached significance in our GSEA analysis for any cell type. This contrasts with the companion phylo analysis ([`phylo_analysis/report_IVD_scRNAseq_analysis.md`, Section 7](../phylo_analysis/report_IVD_scRNAseq_analysis.md#7-pathway-enrichment-gsea)), which reported consistent suppression of Wnt, Notch, and senescence across all cell types using its 7-dataset subset. We address this discordance in Section 8.

### 5.4 Transcription Factor Activity

TF activity inference using CollecTRI regulon overlap identified **290 significant TF-condition associations** (padj < 0.05; Fisher's exact test).

**Key TFs with strongest evidence:**

| TF | Cell Type | Comparison | padj | Targets DE | Direction |
|----|-----------|-----------|------|-----------|-----------|
| E2F4 | NP_mature_chondrocyte | healthy_vs_severe | 8.4x10^-9 | 11/149 | cell cycle |
| HSF1 | Endothelial_cells | healthy_vs_all | 4.8x10^-8 | 8/71 | heat shock |
| HSF1 | AF_inner | mild_vs_severe | 5.0x10^-6 | 5/66 | heat shock |
| E2F1 | NP_mature_chondrocyte | healthy_vs_severe | 2.0x10^-4 | 9/252 | cell cycle |
| HSF2 | Endothelial_cells | healthy_vs_all | 1.7x10^-4 | 4/20 | heat shock |
| EGR1 | NP_stressed_degenerative | mild_vs_severe | 4.7x10^-5 | 7/224 | stress |
| SP1 | NP_stressed_degenerative | mild_vs_severe | 8.7x10^-5 | 10/786 | general |
| RELA | AF_inner | mild_vs_severe | 0.002 | 5/316 | NF-kB |
| NFKB1 | AF_inner | mild_vs_severe | 8.7x10^-4 | 5/230 | NF-kB |
| STAT3 | AF_inner | mild_vs_severe | 0.001 | 5/258 | JAK-STAT |
| FOS | AF_inner | mild_vs_severe | 0.004 | 4/191 | AP-1 |
| ATF7 | NP_stressed_degenerative | mild_vs_severe | 0.003 | 2/5 | stress |
| FOXO3 | NP_stressed_degenerative | mild_vs_severe | 7.4x10^-4 | 5/153 | apoptosis |

**Interpretation:**

1. **E2F4/E2F1 in NP degeneration:** These cell cycle transcription factors are activated in NP_mature_chondrocyte severe degeneration, consistent with the cell cycle pathway enrichment in ORA (Section 5.3). E2F4 typically acts as a repressor of proliferation, and its activation alongside proliferative genes suggests dysregulated cell cycle control — a feature of chondrocyte cluster formation in degenerated discs (Johnson et al., 2001).

2. **HSF1/HSF2 across cell types:** Heat shock factors are among the most significantly activated TFs, consistent with the GSEA heat response enrichment. HSF1 is significant in Endothelial_cells, AF_inner, and NP_stressed_degenerative, suggesting tissue-wide proteotoxic stress. This is consistent with the challenging biophysical environment of the degenerated disc (increased acidity, altered osmolarity, oxidative stress; Wang et al., 2023a).

3. **RELA/NFKB1 in AF:** NF-kB pathway TFs are activated in AF_inner (RELA padj=0.002, NFKB1 padj=8.7x10^-4), directly confirming TNF/NF-kB pathway activation at the transcription factor level — not just at the gene expression level. RELA is the p65 subunit of NF-kB, and its activation drives expression of inflammatory cytokines, MMPs, and ADAMTS enzymes (Wuertz et al., 2012; Xia et al., 2024).

4. **FOXO3 in NP_stressed_degenerative:** FOXO3 (padj=7.4x10^-4) is a key mediator of apoptosis and cellular stress response. Its activation in the stressed/degenerative NP population is consistent with the extrinsic apoptotic signaling pathway enrichment in this cell type.

5. **Expanded TF landscape (v2):** The increase from 113 to 290 significant TF-condition associations in v2 likely reflects the larger number of DE genes available per comparison (particularly the NP_fibrocartilaginous and EP_hyaline cell types), providing more statistical power for TF enrichment testing.

**Figure 9. Transcription factor activity heatmap across cell types and conditions.**
![TF activity heatmap](interpretation/tf_activity/tf_activity_heatmap.png)

### 5.5 Cell State Trajectories

PAGA + diffusion pseudotime analysis revealed structured connectivity between cell states across NP, AF, and CEP compartments.

**Figure 10. NP cell state trajectory — UMAP with pseudotime overlay.**
![NP trajectory UMAP](trajectories/umap_trajectory_NP.png)

**NP trajectory:** Rooted at NP_notochordal cells, the trajectory progresses through NP_mature_chondrocyte and NP_fibrocartilaginous to NP_stressed_degenerative. Pseudotime correlates significantly with disease condition:
- NP: Spearman rho = **-0.258** (p < 10^-100)
- Healthy cells occupy earlier pseudotime; degenerated cells occupy later pseudotime
- This correlation is stronger than v1 (rho=-0.207), suggesting that the de novo annotation and scVI-only integration better resolve the disease-associated continuum.

**Figure 11. NP pseudotime distribution by disease condition.**
![NP pseudotime by condition](trajectories/pseudotime_by_condition_NP.png)

**AF trajectory:** Rooted at AF_inner, progressing toward AF_outer states. Pseudotime-condition correlation:
- AF: Spearman rho = **+0.341** (p < 10^-100)

> **SME REVIEW REQUIRED:** The AF pseudotime-condition correlation is **reversed** compared to v1 (which showed rho=-0.177). In v2, higher pseudotime is associated with healthier tissue in the AF compartment. This reversal may reflect: (1) different cell composition in the v2 AF object (84,624 cells vs. the v1 tier2_resident_AF which included ~283K cells), (2) reassignment of cells through de novo annotation, or (3) a genuine biological pattern where AF_inner-to-AF_outer maturation represents a different axis than degeneration. This finding requires careful examination of the AF object composition and root cell selection before biological interpretation.

**CEP trajectory (new in v2):** The CEP compartment now has sufficient cells for trajectory analysis:
- CEP: Spearman rho = **-0.163** (p < 10^-100)
- EP_hyaline cells show a weaker but significant pseudotime-condition correlation in the expected direction.

**Trajectory-DE overlap:** 500 trajectory-associated genes per compartment were identified. Overlap with DE genes was lower than v1: NP 96/500 genes (19%), AF 110/500 (22%), CEP 38/500 (8%). The reduced overlap compared to v1's ~55% may reflect differences in cell type composition and annotation between versions. The non-overlapping genes may represent gradual, continuous changes not captured by the binary DE framework (e.g., subtle shifts in metabolic gene programs along the continuum).

**Figure 12. Gene expression dynamics along NP pseudotime.**
![NP gene dynamics](trajectories/gene_dynamics_NP.png)

**Figure 13. AF cell state trajectory — UMAP with pseudotime overlay.**
![AF trajectory UMAP](trajectories/umap_trajectory_AF.png)

**Gene dynamics along NP pseudotime:** Notochordal markers (KRT8, KRT18) decline monotonically with pseudotime, while stress/inflammatory markers increase, consistent with the proposed continuum model. Mature chondrocyte markers (ACAN, COL2A1) peak at intermediate pseudotime and decline at the degenerative end, suggesting an initial maintenance phase followed by loss of chondrocyte identity.

### 5.6 Cell-Cell Communication

**Figure 14. Cell-cell interaction heatmap — healthy tissue.**
![Healthy interactions](communication/interaction_plots/interaction_heatmap_healthy.png)

**Figure 15. Cell-cell interaction heatmap — degenerated tissue.**
![Degenerated interactions](communication/interaction_plots/interaction_heatmap_degenerated.png)

**Figure 16. Differential interactions between healthy and degenerated tissue.**
![Differential interactions](communication/interaction_plots/differential_interactions.png)

LIANA consensus analysis identified **28,878 ligand-receptor interactions in healthy** and **27,011 in degenerated** tissue — a **6.5% decrease** in signaling interactions with degeneration.

> **SME REVIEW REQUIRED:** This finding is **reversed** compared to v1 (which showed 53,036 degenerated vs 44,079 healthy, a 20% increase). The v2 result — fewer interactions in degenerated tissue — contrasts with the v1 interpretation of "increased signaling complexity in degeneration." This reversal may reflect: (1) different cell type composition due to de novo annotation, (2) the exclusion of GSE233666, (3) different subsampling outcomes, or (4) sensitivity of LIANA results to the specific cell populations included. The companion phylo analysis found increased interaction strength in severe degeneration, which also disagrees with the v2 finding. This discrepancy underscores the fragility of CCC quantitative comparisons in cross-study meta-analyses.

**Differential interactions:** 36,014 differential interaction pairs were identified between healthy and degenerated conditions.

**Pain-relevant interactions:** 2,077 interactions were flagged as pain-relevant through cross-referencing with curated gene sets (nociception, neurotrophins, nerve guidance, inflammatory pain, neovascularization).

### 5.7 Pain Biology

Cross-referencing DE genes with curated pain gene sets identified **10 significant pain-relevant genes** across comparisons:

- **PTGS2** (COX-2): prostaglandin synthesis enzyme, a direct mediator of inflammatory pain
- **TNF** (x2 comparisons): master inflammatory cytokine and pain mediator
- **PLA2G2A**: phospholipase A2, produces arachidonic acid precursors for prostaglandin synthesis
- **BDKRB2**: bradykinin receptor B2, mediates pain signaling
- **CCL2**: monocyte chemoattractant protein-1, involved in neuroinflammation
- **PTGES**: prostaglandin E synthase, catalyzes PGE2 production
- **CXCL8**: neutrophil-recruiting chemokine with pain-modulating properties

This expanded pain gene set (compared to v1's 3 significant pain genes: TNF x2, CXCL8 x1) provides broader evidence for inflammatory pain mechanisms in disc degeneration.

**Directly supported by our DE data:**
- The prostaglandin axis (PTGS2, PLA2G2A, PTGES) constitutes a coherent pathway from arachidonic acid release to PGE2 production, a direct sensitizer of nociceptive nerve endings (Risbud and Shapiro, 2014).
- CXCL2 (significant in NP_mature_chondrocyte mild_vs_severe) recruits neutrophils and macrophages, which produce additional pain mediators.

**Not detected in our data:**
- NGF (nerve growth factor) and BDNF (brain-derived neurotrophic factor), classically associated with nerve ingrowth into degenerated discs (Freemont et al., 2002), were **not significantly upregulated** in any powered comparison. This may reflect insufficient statistical power or disease stage specificity.

**Figure 17. Pain-associated gene expression heatmap across cell types.**
![Pain genes heatmap](interpretation/pain_genes_heatmap.png)

**Model:** Degenerated disc cells create a pro-inflammatory microenvironment through chemokine and prostaglandin production that promotes nerve ingrowth and sensitization, rather than directly signaling pain. This is consistent with the two-signal model of discogenic pain: (1) structural disruption permits nerve ingrowth into the NP, and (2) the inflammatory milieu sensitizes ingrown nerves (Freemont et al., 2002; Risbud and Shapiro, 2014).

---

## 6. Biological Interpretation and Mechanistic Model

### 6.1 The Inflammatory/Catabolic Cascade

Synthesizing our DE, pathway, TF, and CCC results, we propose that inflammatory cytokine and chemokine production by NP cells contributes to the degenerative cascade:

1. **Initiation:** Mechanical stress, aging, or microinjury activates NF-kB signaling in disc cells (supported by: RELA and NFKB1 TF activation in AF_inner).

2. **Chemokine and cytokine activation:** NF-kB drives expression of inflammatory mediators including CXCL2, ICAM1, and IL32 by NP chondrocytes (supported by: CXCL2 log2FC=+3.14, padj=0.005; ICAM1 log2FC=+1.67, padj=1.9x10^-5 in NP_mature_chondrocyte mild_vs_severe; chemokine pathway enrichment).

3. **Immune cell recruitment:** Chemokines recruit neutrophils and activate macrophages (supported by: neutrophil chemotaxis pathway enrichment; ICAM1 upregulation facilitating immune cell adhesion).

4. **Catabolic cascade:** Recruited immune cells produce additional inflammatory mediators, further degrading the ECM and activating NF-kB in a feed-forward loop (supported by: inflammatory response pathway enrichment; TNC upregulation indicating ECM remodeling).

5. **Cell state deterioration:** Sustained stress drives NP cells along the notochordal → mature chondrocyte → fibrocartilaginous → stressed/degenerative trajectory (supported by: pseudotime-condition rho=-0.258; trajectory-DE overlap).

6. **Metabolic failure in AF:** AF cells experience simultaneous proteotoxic stress (HSP activation) and mitochondrial dysfunction (oxidative phosphorylation downregulation), compromising their ability to maintain the structural integrity of the outer disc (supported by: HSF1 TF activation padj=5.0x10^-6 in AF_inner).

### 6.2 The HSF1 Axis: A Novel Therapeutic Target?

Heat shock factor 1 (HSF1) emerges from our analysis as one of the most consistently activated TFs across cell types and comparisons (significant in Endothelial_cells, AF_inner, NP_stressed_degenerative). HSF1 activation drives expression of heat shock proteins (HSPA1A, HSPA1B, HSPA6, HSP90AA1) that serve as molecular chaperones to refold damaged proteins.

The dual role of HSF1 is therapeutically relevant:
- **Protective:** HSF1-driven HSP expression helps maintain protein homeostasis under stress
- **Inflammatory:** Extracellular HSPs act as damage-associated molecular patterns (DAMPs) that activate TLR2/4 on macrophages, amplifying inflammation (Asea et al., 2002)

This duality suggests that the disc's attempt to cope with proteotoxic stress (via HSF1/HSP activation) paradoxically contributes to inflammation when HSPs are released from dying cells. The timing and location of HSF1 intervention would therefore be critical.

---

## 7. Therapeutic Targets

Based on the evidence from this analysis, we propose the following therapeutic targets, ranked by strength of supporting data.

### 7.1 Tier 1: Strong Direct Evidence From This Analysis

**Target 1: TNF/NF-kB Inhibition**
- **Evidence from this analysis:** RELA TF activation in AF_inner (padj=0.002). NFKB1 TF activation (padj=8.7x10^-4). Multiple inflammatory pathway enrichments driven by NF-kB target genes. TNF is significant in pain gene analysis (2 comparisons). Note: TNF itself did not reach significance in the NP_mature_chondrocyte mild_vs_severe comparison in v2 (padj=0.22), though the NF-kB pathway activation remains robust at the TF level.
- **Mechanism:** NF-kB drives CXCL chemokines, MMPs, and ADAMTS — the entire catabolic cascade.
- **Approach:** Intradiscal anti-TNF biologics (etanercept, adalimumab) or small molecule NF-kB inhibitors. Prior literature supports this approach (Wuertz et al., 2012).
- **Status:** Intradiscal anti-TNF has been proposed; early clinical data available for epidural anti-TNF (Cohen et al., 2009).

**Target 2: HSP/Proteostasis Modulation**
- **Evidence from this analysis:** HSF1 significant in 3 cell types (padj range 5x10^-6 to 3x10^-3). HSPA1A and HSPA1B significantly upregulated in AF_inner and Endothelial_cells. Heat response is a top GSEA pathway in AF_inner.
- **Novel aspect:** HSF1 activation is among the strongest TF signals in our data, yet has received little attention as a disc degeneration therapeutic target. The simultaneous HSP activation and mitochondrial dysfunction suggests an energy crisis — cells are attempting protein rescue but lack metabolic capacity.
- **Approach:** Chemical chaperones (4-PBA, TUDCA) to reduce ER stress and alleviate the need for HSP overexpression. These have shown efficacy in cartilage models (Husa et al., 2013).

**Target 3: Chemokine Modulation (CXCL2/CXCR2)**
- **Evidence from this analysis:** CXCL2 (log2FC=+3.14, padj=0.005) in NP_mature_chondrocyte mild_vs_severe. Chemokine pathway enrichment in ORA. ICAM1 upregulation supporting immune cell recruitment axis.
- **Note:** The evidence for CXC chemokine blockade is weaker in v2 than v1. Only CXCL2 (not the full CXCL1/2/3 triad) reaches significance. This target should be considered alongside broader NF-kB inhibition rather than as a standalone strategy.
- **Mechanism:** CXCL2 signals through CXCR2 on neutrophils. CXCR2 antagonists (e.g., navarixin, AZD5069) have been tested in clinical trials for inflammatory diseases (Rennard et al., 2015).
- **Approach:** Intradiscal delivery of CXCR2 antagonist to block neutrophil/macrophage recruitment without systemic immunosuppression.

### 7.2 Tier 2: Moderate Evidence, Requires Validation

**Target 4: Prostaglandin Pathway Inhibition**
- **Evidence from this analysis:** PTGS2, PLA2G2A, and PTGES identified as significant pain-relevant DE genes. This constitutes a complete prostaglandin synthesis pathway.
- **Mechanism:** PLA2G2A → arachidonic acid → PTGS2/COX-2 → PGH2 → PTGES → PGE2, which directly sensitizes nociceptive nerve endings.
- **Approach:** Selective COX-2 inhibitors (celecoxib) or intradiscal PGE2 receptor antagonists. COX-2 inhibition is already used clinically for pain but could be targeted more precisely to the disc.

**Target 5: Mitochondrial Rescue in AF**
- **Evidence from this analysis:** GSEA shows oxidative phosphorylation, electron transport chain, and mitochondrial ATP synthesis downregulated in AF_inner.
- **Mechanism:** Restoring mitochondrial function could reduce ROS, improve energy metabolism, and support ECM maintenance (Song et al., 2023b).
- **Approach:** Mitochondria-targeted antioxidants (MitoQ, SS-31) or NAD+ precursors (NMN, NR).

**Target 6: E2F4/Cell Cycle Regulation**
- **Evidence from this analysis:** E2F4 is the most significantly activated TF in NP_mature_chondrocyte severe degeneration (padj=8.4x10^-9). Cell cycle pathways dominate ORA for this comparison.
- **Mechanism:** Dysregulated proliferation in degenerated chondrocytes produces the characteristic "cell clusters" seen histologically (Johnson et al., 2001), but these clusters are metabolically inefficient and may deplete local nutrients.
- **Approach:** CDK inhibitors to normalize cell cycle control; however, this is a high-risk target given the already low cellularity of degenerated discs.

### 7.3 Tier 3: Supported by Literature, Not Directly Demonstrated in This Data

**Target 7: ADAMTS5 Inhibition**
- **This analysis:** ADAMTS5 shows a trend toward upregulation in NP_stressed_degenerative and AF_inner in mild_vs_severe, but **does not reach significance** after FDR correction.
- **Literature:** ADAMTS5 is the primary aggrecanase in cartilaginous tissues (Stanton et al., 2005) and is consistently reported as upregulated in disc degeneration (Liang et al., 2022). Our failure to detect significance may reflect underpowering.
- **Status:** Small molecule inhibitors developed for osteoarthritis are in preclinical testing.

**Target 8: TIMP1 Restoration**
- **This analysis:** TIMP1-CD63 loss was not among the top differential interactions in our CCC analysis. However, the companion phylo analysis ([`phylo_analysis/report_IVD_scRNAseq_analysis.md`, Section 9](../phylo_analysis/report_IVD_scRNAseq_analysis.md#9-cell-cell-communication-liana)) identified TIMP1-CD63 as the dominant lost interaction in its 7-dataset CCC analysis.
- **Literature:** The MMP/TIMP balance is a well-established axis of disc degeneration (Vo et al., 2013; Cabral-Pacheco et al., 2020). AAV-TIMP1 gene therapy has shown preclinical efficacy (Han et al., 2021).

**Target 9: Senolytic Therapy**
- **This analysis:** Senescence pathways did not reach significance in our GSEA. However, senescence is well-established in IVD degeneration literature (Song et al., 2023a), and our E2F4/cell cycle TF findings may relate to senescence-associated cell cycle arrest.
- **Literature:** Dasatinib + quercetin senolytics ameliorate disc degeneration in mice (Novais et al., 2021).

### 7.4 Summary Therapeutic Target Table

| Target | Gene(s) | Evidence Level | Key Data Point | Approach |
|--------|---------|---------------|----------------|----------|
| TNF/NF-kB inhibition | RELA, NFKB1 | Strong (TF level) | RELA padj=0.002 | Biologic / small mol |
| HSP modulation | HSF1, HSPA1A/B | Strong (this study) | HSF1 padj=5.0x10^-6 | Chemical chaperone |
| CXCR2 antagonism | CXCL2 | Moderate (this study) | CXCL2 padj=0.005 | Small molecule |
| Prostaglandin inhibition | PTGS2, PLA2G2A | Moderate (pain genes) | 3 pathway genes sig | COX-2 inhibitor |
| Mitochondrial rescue | OXPHOS genes | Moderate (this study) | GSEA suppression | MitoQ / NAD+ |
| E2F4 modulation | E2F4 | Moderate (this study) | padj=8.4x10^-9 | CDK inhibitor |
| ADAMTS5 inhibition | ADAMTS5 | Literature only | Not sig in this study | Small molecule |
| TIMP1 restoration | TIMP1, CD63 | Literature only | Not primary CCC finding | Gene therapy |
| Senolytics | CDKN1A/2A | Literature only | Not sig in GSEA | D+Q |

---

## 8. Novel and Discordant Findings

### 8.1 Attenuated CXC Chemokine Signal in v2

The v1 analysis identified CXCL1/2/3 as the most significantly DE genes in NP severe degeneration — a finding that was novel in the single-cell context. In v2, only CXCL2 retains significance (padj=0.005), while CXCL1 (no significant padj), CXCL3 (padj=0.099), and TNF (padj=0.22) do not reach the significance threshold. This attenuation likely reflects: (1) the exclusion of GSE233666, (2) different cell type boundaries from de novo annotation, and (3) the use of scVI-only integration (vs. scANVI in v1). The inflammatory signal is still present — chemokine pathway enrichment persists in ORA — but the evidence for CXC chemokine dominance is weaker. The CXCL2 signal, as a GRO-beta chemokine that recruits neutrophils, still supports a model of immune cell recruitment in severe degeneration.

### 8.2 NP_fibrocartilaginous: A New Transcriptionally Active Cell Type

The NP_fibrocartilaginous cell type, newly defined through de novo annotation in v2, shows substantial transcriptomic responsiveness to degeneration (203 DE genes in mild_vs_severe, 127 in healthy_vs_severe). This population, which likely represents transitional cells producing type I collagen in the NP, warrants further characterization as a potential driver of fibrocartilaginous replacement — a hallmark of disc degeneration (Antoniou et al., 1996).

### 8.3 CEP Transcriptional Silencing

The EP_hyaline compartment shows an asymmetric DE pattern: 84 downregulated genes and 0 upregulated in healthy_vs_degenerated_all. This pattern of transcriptional silencing — rather than inflammatory activation — in CEP degeneration has not been previously described at single-cell resolution and may reflect endplate calcification and reduced metabolic activity.

### 8.4 CXCL8 Compartment Specificity: Potentially Novel

CXCL8 (IL-8) appeared among pain-relevant DE genes. Compartment-specific chemokine patterns have not been extensively characterized in the IVD literature and may reflect the distinct microenvironments of the NP (avascular, hypoxic) versus AF (partially vascularized).

### 8.5 Discordance with Companion Phylo Analysis: Wnt, Notch, and Senescence

The companion phylo analysis ([`phylo_analysis/report_IVD_scRNAseq_analysis.md`, Section 7](../phylo_analysis/report_IVD_scRNAseq_analysis.md#7-pathway-enrichment-gsea); Good, 2026) — which analyzed 7 of the same 11 datasets using Harmony integration and R-based DESeq2 — reported consistent suppression of Wnt signaling, Notch signaling, and cellular senescence pathways across all NP cell types. Our analysis did not replicate these findings. Several factors likely contribute:

1. **Histone gene artifacts:** Examination of the phylo analysis's GSEA results reveals that the top 30+ enriched pathways (including Wnt, Notch, senescence, and DNA methylation) are driven almost exclusively by the same set of replication-dependent histone genes (H4C15, H4C11, H2BC12, H2AC8, etc.). These genes appear in the core enrichment of essentially all "suppressed" pathways because many Reactome pathways include histone-related genes. Histone genes are known to be highly sensitive to cell cycle state, dissociation protocols, and ambient RNA contamination (Slyper et al., 2020), making them unreliable indicators of pathway activity in cross-study comparisons.

2. **Comparison design:** The phylo analysis used healthy_vs_severe comparisons, where study and condition are maximally confounded. Our prioritization of mild_vs_severe (within-study) comparisons reduces this confounding.

3. **LFC shrinkage:** Our pyDESeq2 analysis applies LFC shrinkage by default, constraining fold-change estimates to biologically plausible ranges. The phylo analysis reported log2FC values up to -28 (LINC01578), which likely reflect technical artifacts rather than true expression changes of 10^8-fold magnitude.

4. **Integration method:** scVI (deep learning) vs. Harmony (linear correction) may differentially preserve cell state heterogeneity, affecting which genes appear as DE.

**Our interpretation:** The Wnt, Notch, and senescence pathway suppression reported in the phylo analysis is substantially driven by a histone gene artifact that propagates through pathway databases. This does not mean these pathways are unaltered in disc degeneration — literature evidence for Wnt/Notch involvement is substantial (Li et al., 2023a; Long et al., 2019) — but our data does not independently confirm these pathway changes at the GSEA level. The inflammatory/chemokine signature we detect is more robust because it is driven by diverse, biologically coherent gene sets.

### 8.6 TIMP1-CD63 Not Replicated

The companion phylo analysis ([`phylo_analysis/report_IVD_scRNAseq_analysis.md`, Section 9](../phylo_analysis/report_IVD_scRNAseq_analysis.md#9-cell-cell-communication-liana)) identified TIMP1-CD63 loss as the dominant CCC change, alongside FN1-driven macrophage signaling gains. Our CCC analysis did not replicate this as a primary finding. This discordance likely reflects differences in: (1) cell type resolution, (2) CCC methodology (per-dataset vs. integrated), and (3) subsampling strategies. The TIMP1-CD63 finding remains biologically plausible (Vo et al., 2013) and warrants targeted investigation.

---

## 9. Limitations

1. **Cross-study confounding:** Condition and study are partially confounded. Herniated comparisons were excluded entirely due to single-study confounding (only GSE251686 after GSE233666 exclusion). Within-study comparisons (mild_vs_severe) are prioritized throughout this manuscript.

2. **Underpowered comparisons:** 53 of 74 cell type x comparison combinations were skipped due to insufficient samples (< 3 per condition). Key genes like ADAMTS5, ACAN, and COL2A1 may fail to reach significance due to donor variability rather than absence of change.

3. **Age-disease confound:** In GSE230809 (the largest dataset, 24 samples), healthy donors are 21-27 years old and diseased are 37-73 years old. Age and degeneration effects cannot be fully separated.

4. **Sex bias:** GSE230809 is all-male. Many samples have unknown sex. Sex-stratified analyses are not possible.

5. **Culture-expanded cells:** GSE242443 CEP cells are culture-expanded, which alters gene expression (particularly collagen ratios and surface markers).

6. **No RNA velocity:** Spliced/unspliced counts were not available from public deposits. RNA velocity would provide directional evidence for cell state transitions.

7. **No SCENIC/GRN:** Full SCENIC analysis (gene regulatory networks) was not performed due to computational requirements. TF activity was estimated from CollecTRI regulon overlap, which captures target enrichment but not regulatory network structure.

8. **Composition analysis underpowered:** No cell type proportion changes reached significance after FDR correction, though trends were biologically consistent (e.g., reduced NP notochordal cells, increased AF cells in degeneration).

9. **NP non-mesenchymal annotation quality:** CellTypist validation showed 8 of 13 NP non-mesenchymal clusters were discordant with automated predictions. Contributing factors include ~17K misrouted mesenchymal cells in non-mesenchymal clusters, scoring formula issues in the de novo annotation, and the CellTypist Immune_All_Low model lacking a pericyte category. Non-mesenchymal cell type labels in the NP object should be interpreted with caution pending annotation refinement.

10. **CCC methodology and fragility:** LIANA was run on per-dataset files (not integrated data), which avoids integration artifacts but fragments the analysis across datasets. The reversal of the healthy-vs-degenerated interaction count pattern between v1 and v2 (see Section 5.6) highlights the sensitivity of CCC quantification to cell type definitions and subsampling.

11. **AF trajectory reversal:** The positive pseudotime-condition correlation in AF (rho=+0.341) is reversed from v1 (rho=-0.177). This may reflect cell composition changes from de novo annotation or root cell selection artifacts. AF trajectory results should be considered preliminary pending SME review.

---

## 10. Conclusion

This 11-dataset, 410,759-cell meta-analysis of human IVD degeneration reveals an inflammatory transcriptomic signature in severe NP degeneration, with CXCL2 upregulation, ICAM1-mediated immune cell recruitment, and NF-kB pathway activation confirmed at the transcription factor level. This signature is accompanied by heat shock protein activation and mitochondrial dysfunction in AF cells, representing a tissue-wide stress response. A newly characterized NP_fibrocartilaginous population shows substantial transcriptomic responsiveness to degeneration, and the CEP compartment — now powered for DE analysis — shows a pattern of transcriptional silencing rather than activation. Cell state trajectory analysis confirms that NP cells exist on a disease-associated continuum from notochordal to degenerative states, with pseudotime correlating with clinical disease severity (rho=-0.258).

Several findings require SME review before firm biological conclusions: the AF trajectory reversal (positive pseudotime-condition correlation), the decreased cell-cell communication in degenerated tissue (contrasting with v1 and the companion phylo analysis), and the attenuated CXC chemokine signal compared to v1. These discrepancies underscore the sensitivity of single-cell meta-analysis results to methodological choices including dataset selection, integration method, cell type annotation strategy, and subsampling.

The primary therapeutic opportunities emerging from this analysis are TNF/NF-kB inhibition (supported by TF-level evidence), HSP/proteostasis modulation (supported by the strongest TF signals), and chemokine pathway modulation (supported by CXCL2 DE and pathway enrichment). The prostaglandin synthesis pathway (PTGS2, PLA2G2A, PTGES) is newly identified as a coherent pain-relevant target. Classical targets such as ADAMTS5 and TIMP1 remain valid based on extensive literature but were not independently confirmed in our powered comparisons.

Importantly, this analysis highlights the sensitivity of scRNA-seq meta-analysis results to methodological choices: integration method, comparison design, LFC shrinkage, cell type annotation, and artifact awareness all substantially impact biological conclusions. The v1-to-v2 changes documented here — including the attenuation of the CXC chemokine signal, the reversal of AF trajectory direction, and the reversal of CCC interaction counts — serve as a cautionary example of how analytical pipeline decisions shape biological narratives from the same underlying data.

---

## 11. References

Adams MA, Roughley PJ. (2006). What is intervertebral disc degeneration, and what causes it? *Spine*, 31(18):2151-2161.

Antoniou J, Steffen T, Nelson F, et al. (1996). The human lumbar intervertebral disc: evidence for changes in the biosynthesis and denaturation of the extracellular matrix. *Journal of Clinical Investigation*, 98(4):996-1003.

Asea A, Rehli M, Kabingu E, et al. (2002). Novel signal transduction pathway utilized by extracellular HSP70. *Journal of Biological Chemistry*, 277(17):15028-15034.

Cabral-Pacheco GA, Garza-Veloz I, Castruita-De la Rosa C, et al. (2020). The Roles of Matrix Metalloproteinases and Their Inhibitors in Human Diseases. *International Journal of Molecular Sciences*, 21:9739.

Cohen SP, Bogduk N, Dragovich A, et al. (2009). Randomized, double-blind, placebo-controlled, dose-response, and preclinical safety study of transforaminal epidural etanercept for the treatment of sciatica. *Anesthesiology*, 110(5):1116-1126.

Dieleman JL, Cao J, Chapin A, et al. (2020). US health care spending by payer and health condition, 1996-2016. *JAMA*, 323(9):863-884.

Dimitrov D, Turei D, Garrber M, et al. (2022). Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. *Nature Communications*, 13:3224.

Dominguez Conde C, Xu C, Jarvis LB, et al. (2022). Cross-tissue immune cell analysis reveals tissue-specific features in humans. *Science*, 376(6594):eabl5197.

Fang Z, Liu X, Peltz G. (2023). GSEApy: a comprehensive package for performing gene set enrichment analysis in Python. *Bioinformatics*, 39(1):btac757.

Fernandes LM, Khan N, Trochez CM, et al. (2020). Single-cell RNA-seq identifies unique transcriptional landscapes of human nucleus pulposus and annulus fibrosus cells. *Scientific Reports*, 10:15263.

Freemont AJ, Watkins A, Le Maitre C, et al. (2002). Nerve growth factor expression and innervation of the painful intervertebral disc. *Journal of Pathology*, 197(3):286-292.

Gan Y, He J, Zhu J, et al. (2021). Spatially defined single-cell transcriptional profiling characterizes diverse chondrocyte subtypes and nucleus pulposus progenitors in human intervertebral discs. *Bone Research*, 9:37.

Garcia-Alonso L, Holland CH, Ibrahim MM, et al. (2019). Benchmark and integration of resources for the estimation of human transcription factor activities. *Genome Research*, 29(8):1363-1375.

GBD 2021 Low Back Pain Collaborators. (2023). Global, regional, and national burden of low back pain, 1990-2020. *The Lancet Rheumatology*, 5(6):e316-e329.

Good B. (2026). Single-Cell Transcriptomic Atlas of Human Intervertebral Disc Degeneration. Draft manuscript, Phylo/Biomni analysis. Report: [`phylo_analysis/report_IVD_scRNAseq_analysis.md`](../phylo_analysis/report_IVD_scRNAseq_analysis.md); Manuscript: [`phylo_analysis/draft_manuscript.md`](../phylo_analysis/draft_manuscript.md).

Haghverdi L, Buttner M, Wolf FA, Buettner F, Theis FJ. (2016). Diffusion pseudotime robustly reconstructs lineage branching. *Nature Methods*, 13:845-848.

Han Y, Ouyang Z, Wawrose R, et al. (2021). ISSLS prize in basic science 2021: a novel inducible system to regulate transgene expression of TIMP1. *European Spine Journal*, 30:1098-1107.

Husa M, Petursson F, Loer R, et al. (2013). C/EBP homologous protein drives pro-catabolic responses in chondrocytes. *Arthritis Research & Therapy*, 15:R218.

Johnson WE, Eisenstein SM, Roberts S. (2001). Cell cluster formation in degenerate lumbar intervertebral discs is associated with increased disc cell proliferation. *Connective Tissue Research*, 42(3):197-207.

Levin D, Azar S, Engel A. (2019). A Randomized, Double-blind, Active-control, Multi-center Study of Hyaluronic Acid vs Corticosteroid for Intradiscal Injection for the Treatment of Lumbar Discogenic Pain. *Spine*, 44(16):1127-1135.

Li X, Han Y, Li G, et al. (2023a). Role of Wnt signaling pathway in joint development and cartilage degeneration. *Frontiers in Cell and Developmental Biology*, 11:1181619.

Li Z, Ye D, Dai L, et al. (2022a). Single-Cell RNA Sequencing Reveals the Difference in Human Normal and Degenerative Nucleus Pulposus Tissue Profiles and Cellular Interactions. *Frontiers in Cell and Developmental Biology*, 10:910626.

Liang H, Luo R, Li G, et al. (2022). The Proteolysis of ECM in Intervertebral Disc Degeneration. *International Journal of Molecular Sciences*, 23:1715.

Liberzon A, Birger C, Thorvaldsdottir H, et al. (2015). The Molecular Signatures Database Hallmark gene set collection. *Cell Systems*, 1(6):417-425.

Long J, Wang X, Du X, et al. (2019). JAG2/Notch2 inhibits intervertebral disc degeneration by modulating cell proliferation, apoptosis, and extracellular matrix. *Arthritis Research & Therapy*, 21:213.

Lopez R, Regier J, Cole MB, et al. (2018). Deep generative modeling for single-cell transcriptomics. *Nature Methods*, 15:1053-1058.

Love MI, Huber W, Anders S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15:550.

Luecken MD, Buttner M, Chaichoompu K, et al. (2022). Benchmarking atlas-level data integration in single-cell genomics. *Nature Methods*, 19:41-50.

Martel-Pelletier J, Barr AJ, Cicuttini FM, et al. (2020). Osteoarthritis. *Nature Reviews Disease Primers*, 2:16072.

Novais EJ, Tran VA, Johnston SN, et al. (2021). Long-term treatment with senolytic drugs dasatinib and quercetin ameliorates age-dependent intervertebral disc degeneration in mice. *Nature Communications*, 12:5213.

Oichi T, Taniguchi Y, Oshima Y, et al. (2020). Pathomechanism of intervertebral disc degeneration. *JOR Spine*, 3:e1076.

Onishi RM, Gaffen SL. (2010). Interleukin-17 and its target genes: mechanisms of interleukin-17 function in disease. *Immunology*, 129(3):311-321.

Rennard SI, Dale DC, Donohue JF, et al. (2015). CXCR2 Antagonist MK-7123. A Phase 2 Proof-of-Concept Trial for Chronic Obstructive Pulmonary Disease. *American Journal of Respiratory and Critical Care Medicine*, 191(9):1001-1011.

Risbud MV, Shapiro IM. (2014). Role of cytokines in intervertebral disc degeneration: pain and disc content. *Nature Reviews Rheumatology*, 10(1):44-56.

Slyper M, Porter CBM, Ashenberg O, et al. (2020). A single-cell and single-nucleus RNA-Seq toolbox for fresh and frozen human tumors. *Nature Medicine*, 26:792-802.

Song C, Cai W, Liu F, et al. (2022). An in-depth analysis of the immunomodulatory mechanisms of intervertebral disc degeneration. *JOR Spine*, 5:e1233.

Song C, Zhou Y, Cheng K, et al. (2023a). Cellular senescence — Molecular mechanisms of intervertebral disc degeneration from an immune perspective. *Biomedicine & Pharmacotherapy*, 162:114711.

Song C, Xu Y, Peng Q, et al. (2023b). Mitochondrial dysfunction: a new molecular mechanism of intervertebral disc degeneration. *Inflammation Research*, 72:2249-2260.

Squair JW, Gautier M, Kathe C, et al. (2021). Confronting false discoveries in single-cell differential expression. *Nature Communications*, 12:5692.

Stanton H, Rogerson FM, East CJ, et al. (2005). ADAMTS5 is the major aggrecanase in mouse cartilage in vivo and in vitro. *Nature*, 434:648-652.

Vo N, Hartman R, Yurube T, et al. (2013). Expression and regulation of metalloproteinases and their inhibitors in intervertebral disc aging and degeneration. *The Spine Journal*, 13:331-341.

Wang Y, Cheng H, Wang T, et al. (2023a). Oxidative stress in intervertebral disc degeneration: Molecular mechanisms, pathogenesis and treatment. *Cell Proliferation*, 56:e13448.

Wolock SL, Lopez R, Klein AM. (2019). Scrublet: computational identification of cell doublets in single-cell transcriptomic data. *Cell Systems*, 8(4):281-291.e9.

Wuertz K, Vo N, Kletsas D, Boos N. (2012). Inflammatory and catabolic signalling in intervertebral discs: the roles of NF-kB and MAP kinases. *European Cells and Materials*, 23:103-120.

Xia Q, Zhao Y, Dong H, et al. (2024). Progress in the study of molecular mechanisms of intervertebral disc degeneration. *Biomedicine & Pharmacotherapy*, 174:116593.

Yoshida H, Nagaoka A, Kusaka-Kikushima A, et al. (2013). KIAA1199, a deafness gene of unknown function, is a new hyaluronan binding protein involved in hyaluronan depolymerization. *Proceedings of the National Academy of Sciences*, 110(14):5612-5617.

Zimmerman KD, Espeland MA, Langefeld CD. (2021). A practical solution to pseudoreplication bias in single-cell studies. *Nature Communications*, 12:738.

---

*Analysis performed using a 10-module human-gated agentic pipeline (v2). All code version-controlled. Random seed: 42. Package versions: Python 3.12, scanpy 1.11, scvi-tools 1.4.2, pyDESeq2, gseapy 1.1, decoupler 2.1, liana 1.7. 19 supplementary tables (S1-S19, including S17-S19 for CellTypist concordance).*

*This is a computational analysis draft. All findings require experimental validation before clinical application.*
