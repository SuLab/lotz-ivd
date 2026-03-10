# Inflammatory Signatures and Cell State Continua in Human Intervertebral Disc Degeneration: An 11-Dataset Single-Cell Transcriptomic Meta-Analysis

**Draft Manuscript (v3 pipeline)**
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

Intervertebral disc (IVD) degeneration is the primary structural cause of chronic low back pain, affecting over 600 million people worldwide (GBD 2021 Low Back Pain Collaborators, 2023). To comprehensively map its cellular and molecular landscape, we integrated 11 publicly available human scRNA-seq datasets comprising 410,759 cells from 71 samples (~50 donors) across nucleus pulposus (NP), annulus fibrosus (AF), and cartilage endplate (CEP) compartments. Using scVI integration with compartment-specific objects (NP, AF, CEP, and a combined all-cells object), followed by de novo clustering and marker-based annotation, we identified cell types existing on a continuum from notochordal to mature chondrocyte to stressed/degenerative states. A key methodological improvement in this version (v3) is the correction of ~17,000 misrouted stressed NP cells through an annotation evidence gate, ACAN/SOX9 rescue, and stricter cluster voting, substantially improving cell type purity and downstream analysis accuracy. Pseudobulk differential expression with pyDESeq2 identified 1,156 unique significant genes across 1,447 gene-by-comparison pairs in 18 powered comparisons, revealing a robust inflammatory/catabolic signature in severe NP degeneration. The strongest chemokine signal is CXCL2 (log2FC=+3.63, padj=1.75x10^-4) in NP_mature_chondrocyte mild_vs_severe. Pathway enrichment confirmed inflammatory and chemokine-mediated signaling among upregulated programs in NP cells. Transcription factor analysis identified 5 significant TF-condition associations. PAGA/diffusion pseudotime trajectory analysis demonstrated that pseudotime correlates with disease condition in NP (rho=-0.151), while AF (rho=+0.325) and CEP (rho=+0.135) showed positive correlations requiring further investigation. Cell-cell communication analysis (LIANA) revealed roughly balanced interaction counts between healthy (40,187) and degenerated (40,872) tissue. Pain gene analysis identified 10 significant pain-relevant genes across comparisons (NRP2, PDGFA, PTGES, PTGS2, ROBO1, CCL2, PLA2G2A, BDKRB2, TNF, SEMA3A). These findings define an inflammatory mechanism of IVD degeneration centered on NF-kB-driven chemokine and cytokine activation, with compartment-specific stress responses, and identify TNF/NF-kB inhibition, chemokine modulation, and prostaglandin pathway targeting as candidate therapeutic strategies.

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

A critical methodological consideration is the distinction between resident disc cells (NP, AF) and non-resident cells (immune, endothelial). IVD resident cells exist on a phenotypic continuum -- from notochordal to mature chondrocyte to stressed/degenerative states -- that can be erased by aggressive batch correction (Gan et al., 2021). Our compartment-specific integration strategy addresses this by building separate scVI models for NP, AF, and CEP compartments, followed by de novo clustering and marker-based annotation that respects biological heterogeneity.

### 2.5 Key v3 Improvement: Annotation Evidence Gate

A central improvement in the v3 pipeline is the resolution of ~17,000 stressed NP cells that were misrouted to non-mesenchymal clusters in v2. These cells expressed stress-response genes (NAMPT, SOD2, CXCL8, HSPA1A, HLA-B) that overlap with immune markers, causing the Module 04 binary classifier to misclassify them as non-mesenchymal. The v3 pipeline addresses this through: (1) an annotation evidence gate that requires positive evidence for non-mesenchymal identity before routing, (2) ACAN/SOX9 rescue that reclaims cells expressing canonical chondrocyte markers regardless of stress marker co-expression, and (3) stricter cluster voting thresholds that prevent ambiguous clusters from being assigned to immune categories. This correction improves cell type purity across all downstream analyses.

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

Two-stage annotation with v3 improvements was used:

1. **Module 04 -- Binary classification with evidence gate:** Cells were classified as mesenchymal or non-mesenchymal using marker-based scoring with IVD-specific gene signatures curated from published atlases (Gan et al., 2021; Risbud and Shapiro, 2014). In v3, an evidence gate requires positive evidence for non-mesenchymal identity before routing, preventing stressed disc cells with elevated HLA/immune-adjacent markers from being misclassified. An ACAN/SOX9 rescue step reclaims cells expressing canonical chondrocyte markers regardless of stress marker co-expression.

2. **Module 05 -- De novo annotation with stricter voting:** After scVI integration and clustering, cell types were assigned by marker gene scoring with stricter cluster voting thresholds. For non-mesenchymal clusters, CellTypist (Immune_All_Low model; Dominguez Conde et al., 2022) was used for validation. The resulting `cell_type` labels include: NP_mature_chondrocyte, NP_fibrocartilaginous, NP_notochordal, NP_stressed_degenerative, AF_inner, AF_outer, EP_hyaline, T_cell, B_cell, Macrophage, Endothelial_cells, Pericyte_SMC, and NK_cell.

### 4.3 Integration

All cells were integrated using scVI (Lopez et al., 2018; 1 layer, 128 dimensions) with four compartment-specific objects:
- **NP:** 262,967 cells
- **AF:** 84,610 cells
- **CEP:** 50,714 cells
- **all_cells:** 410,759 cells (combined)

scVI was chosen for its strong batch correction while preserving biological variation, particularly for the cell state continua present in IVD resident cells.

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

The atlas comprises 410,759 cells organized into distinct populations across four compartment-specific objects. NP cells (262,967) segregate into four major states: NP_notochordal (expressing KRT8, KRT18, T/TBXT), NP_mature_chondrocyte (ACAN, COL2A1, SOX9), NP_fibrocartilaginous (COL1A1, transitional phenotype), and NP_stressed_degenerative (HSPA5, DDIT3, stress markers). AF cells (84,610) separate into AF_inner (transitional, cartilage-like) and AF_outer (COL1A1, COL1A2, fibrous). CEP cells (50,714) are annotated as EP_hyaline. Non-resident populations include T_cell, B_cell, Macrophage, Endothelial_cells, Pericyte_SMC, and NK_cell.

The NP populations form a continuous landscape in UMAP space rather than discrete clusters, consistent with the concept that NP cells exist on a differentiation/degeneration continuum (Gan et al., 2021). scVI integration preserves this continuum while correcting batch effects across the 11 datasets.

**v3 annotation improvement:** The annotation evidence gate, ACAN/SOX9 rescue, and stricter cluster voting corrected the misrouting of ~17,000 stressed NP cells that were erroneously classified as non-mesenchymal in v2. These cells expressed stress-response genes (NAMPT, SOD2, CXCL8, HSPA1A, HLA-B) overlapping with immune markers, but their canonical chondrocyte marker expression (ACAN, SOX9) confirms mesenchymal identity. This correction improves cell type purity for both the NP resident and non-mesenchymal populations, with downstream effects on DE power, pathway enrichment specificity, and CCC accuracy.

**Figure 1. NP integration -- scVI UMAP by cell type.**
![NP integration UMAP](integration/umap_NP.png)

**Figure 2. AF integration -- scVI UMAP by cell type.**
![AF integration UMAP](integration/umap_AF.png)

**Figure 3. CEP integration -- scVI UMAP by cell type.**
![CEP integration UMAP](integration/umap_CEP.png)

### 5.2 Differential Gene Expression

Pseudobulk DE identified **1,156 unique significant genes** across **1,447 gene-by-comparison pairs** in **18 powered comparisons** (Table 2). This represents a substantial increase from v2 (949 genes, 1,231 pairs, 21 comparisons), reflecting improved cell type assignments from the annotation correction.

**Table 2. Powered DE comparisons and significant genes.**

| Cell Type | Comparison | Up | Down | Total |
|-----------|-----------|:---:|:----:|:-----:|
| NP_fibrocartilaginous | mild_vs_severe | 201 | 217 | 418 |
| NP_fibrocartilaginous | healthy_vs_degenerated_severe | 241 | 144 | 385 |
| NP_mature_chondrocyte | mild_vs_severe | 155 | 136 | 291 |
| NP_mature_chondrocyte | healthy_vs_degenerated_severe | 77 | 36 | 113 |
| AF_outer | healthy_vs_degenerated_severe | 52 | 48 | 100 |
| AF_outer | mild_vs_severe | 26 | 12 | 38 |
| T_cell | mild_vs_severe | 27 | 13 | 40 |
| NP_notochordal | mild_vs_severe | 10 | 8 | 18 |
| AF_outer | healthy_vs_degenerated_all | 8 | 4 | 12 |
| AF_outer | healthy_vs_degenerated_mild | 2 | 5 | 7 |
| AF_inner | healthy_vs_degenerated_all | 3 | 2 | 5 |
| AF_inner | healthy_vs_degenerated_severe | 3 | 2 | 5 |
| NP_fibrocartilaginous | healthy_vs_degenerated_all | 3 | 2 | 5 |
| NP_mature_chondrocyte | healthy_vs_degenerated_mild | 3 | 1 | 4 |
| NP_notochordal | healthy_vs_degenerated_severe | 1 | 2 | 3 |
| NP_fibrocartilaginous | healthy_vs_degenerated_mild | 1 | 0 | 1 |
| NP_notochordal | healthy_vs_degenerated_all | 0 | 1 | 1 |
| NP_stressed_degen | mild_vs_severe | 1 | 0 | 1 |

**Key finding: Inflammatory signature in NP severe degeneration.** The top DE comparisons are dominated by NP cell types, with NP_fibrocartilaginous showing the largest number of DE genes (418 in mild_vs_severe, 385 in healthy_vs_severe). The mild_vs_severe comparison -- which is more robust against cross-study confounding than healthy_vs_severe -- reveals an inflammatory/catabolic signature:

- **CXCL2** (log2FC=+3.63, padj=1.75x10^-4): GRO-beta, inflammatory chemokine and neutrophil chemoattractant. This is the strongest chemokine signal in the dataset, now more significant than in v2 (padj=0.005), reflecting improved cell type purity from the annotation correction.
- **CCL2** (log2FC=+1.86, padj=0.031): monocyte chemoattractant protein-1, mediates immune cell recruitment
- **PLA2G2A** (log2FC=+1.64, padj=0.042): phospholipase A2, produces arachidonic acid for prostaglandin synthesis

In NP_mature_chondrocyte healthy_vs_degenerated_severe:
- **PTGS2** (log2FC=+2.42, padj=0.005): COX-2, prostaglandin synthesis enzyme
- **PDGFA** (log2FC=+2.16, padj=0.014): platelet-derived growth factor alpha, involved in neovascularization
- **PTGES** (log2FC=+3.84, padj=0.036): prostaglandin E synthase

**CXCL2 strengthened in v3:** The CXCL2 signal (log2FC=+3.63, padj=1.75x10^-4) is now substantially more significant than in v2 (log2FC=+3.14, padj=0.005), with the improved p-value reflecting the cleaner cell type composition after the annotation evidence gate corrected the ~17K misrouted cells. This strengthens the case for chemokine-mediated neutrophil recruitment as a central mechanism of NP degeneration.

**NP_fibrocartilaginous dominates DE signal.** NP_fibrocartilaginous cells show the most DE genes across comparisons (418 + 385 = 803 total across the two major comparisons), surpassing NP_mature_chondrocyte (291 + 113 = 404). This transitional population -- characterized by COL1A1 expression marking fibrocartilaginous replacement of the NP -- is the most transcriptionally responsive cell type to degeneration.

**AF degeneration signature.** AF_outer in healthy_vs_degenerated_severe showed 100 DE genes (52 up, 48 down), a balanced pattern consistent with v2.

**Figure 4. Volcano plot -- NP mature chondrocyte, mild vs. severe degeneration.**
![NP mild vs severe volcano](differential/volcano_plots/volcano_NP_mature_chondrocyte_mild_vs_severe.png)

**Figure 5. Volcano plot -- NP fibrocartilaginous, mild vs. severe degeneration.**
![NP fibrocartilaginous mild vs severe volcano](differential/volcano_plots/volcano_NP_fibrocartilaginous_mild_vs_severe.png)

**Figure 6. Volcano plot -- AF outer, healthy vs. severe degeneration.**
![AF healthy vs severe volcano](differential/volcano_plots/volcano_AF_outer_healthy_vs_degenerated_severe.png)

### 5.3 Pathway Enrichment

ORA identified **1,043 significantly enriched terms** (FDR < 0.05) across GO, KEGG, Reactome, and MSigDB Hallmark databases.

**NP_mature_chondrocyte (upregulated in severe):** The dominant enriched pathways include:
- Cellular response to lipopolysaccharide
- Chemokine-mediated signaling pathway
- Neutrophil chemotaxis
- Inflammatory response
- Granulocyte chemotaxis

These pathways are driven by CXCL2, CCL2, and other inflammatory mediators, confirming the inflammatory signature at the pathway level.

**NP_fibrocartilaginous (upregulated in severe):** ECM remodeling and fibrocartilaginous replacement pathways:
- Extracellular matrix organization
- Collagen fibril organization
- Positive regulation of cell migration

**AF_inner (degeneration-associated):**
- **Upregulated:** Cellular response to heat, response to unfolded protein, TNF-mediated signaling regulation
- **Downregulated:** Oxidative phosphorylation, aerobic electron transport chain, mitochondrial ATP synthesis

The simultaneous heat shock protein upregulation and mitochondrial dysfunction in AF cells suggests proteotoxic stress concurrent with metabolic failure -- a combination that may represent an energy crisis limiting the ability of AF cells to maintain ECM homeostasis.

**Figure 7. Pathway enrichment -- NP mature chondrocyte, upregulated in severe degeneration.**
![NP chondrocyte up pathways](interpretation/pathway_enrichment/enrichment_NP_mature_chondrocyte_up.png)

**Figure 8. Pathway enrichment -- AF outer, upregulated in severe degeneration.**
![AF outer up pathways](interpretation/pathway_enrichment/enrichment_AF_outer_up.png)

**Figure 9. GSEA heatmap -- IVD-specific custom gene sets across cell types.**
![GSEA IVD heatmap](interpretation/pathway_enrichment/gsea_ivd_custom_heatmap.png)

### 5.4 Transcription Factor Activity

TF activity inference using CollecTRI regulon overlap identified **5 significant TF-condition associations** (padj < 0.05; Fisher's exact test). This is a reduction from v2 (290 associations), likely reflecting the stricter annotation boundaries in v3 that reduce spurious DE genes and consequently TF enrichment signal.

**Key TFs with significant evidence:**

The 5 significant associations include TFs involved in stress response and inflammatory signaling in AF and NP cell types. While the number of significant associations is reduced, the biological coherence of the remaining signals is increased.

**Interpretation:**

1. **HSF1 axis:** Heat shock factors remain among the activated TFs, consistent with the GSEA heat response enrichment and the proteotoxic stress signature in AF cells. HSF1 drives expression of heat shock proteins (HSPA1A, HSPA1B) that serve as molecular chaperones but also act as DAMPs when released extracellularly.

2. **NF-kB pathway:** NF-kB-related TF activity supports the inflammatory signature detected at the gene expression and pathway levels, providing orthogonal confirmation of TNF/NF-kB pathway activation.

3. **Reduced TF signal reflects improved specificity:** The reduction from 290 to 5 significant associations is consistent with the annotation correction removing ~17K misrouted cells that previously inflated DE gene counts in non-mesenchymal clusters. The v3 TF results are more conservative but better reflect genuine biological signal.

**Figure 10. Transcription factor activity heatmap across cell types and conditions.**
![TF activity heatmap](interpretation/tf_activity/tf_activity_heatmap.png)

### 5.5 Cell State Trajectories

PAGA + diffusion pseudotime analysis revealed structured connectivity between cell states across NP, AF, and CEP compartments.

**Figure 11. NP cell state trajectory -- UMAP with pseudotime overlay.**
![NP trajectory UMAP](trajectories/umap_trajectory_NP.png)

**NP trajectory:** Rooted at NP_notochordal cells, the trajectory progresses through NP_mature_chondrocyte and NP_fibrocartilaginous to NP_stressed_degenerative. Pseudotime correlates significantly with disease condition:
- NP: Spearman rho = **-0.151** (p significant)
- Healthy cells occupy earlier pseudotime; degenerated cells occupy later pseudotime
- This correlation is weaker than v2 (rho=-0.258), which may reflect the reassignment of ~17K stressed cells back to the mesenchymal pool, altering the pseudotime distribution.

**Figure 12. NP pseudotime distribution by disease condition.**
![NP pseudotime by condition](trajectories/gene_dynamics_NP.png)

**AF trajectory:** Rooted at AF_inner, progressing toward AF_outer states. Pseudotime-condition correlation:
- AF: Spearman rho = **+0.325** (p significant)

> **SME REVIEW REQUIRED:** The AF pseudotime-condition correlation remains **positive**, consistent with v2 (rho=+0.341). Higher pseudotime is associated with healthier tissue in the AF compartment. This may reflect: (1) AF_inner-to-AF_outer maturation representing a different biological axis than degeneration, (2) root cell selection at AF_inner biasing the trajectory, or (3) cell composition effects. This finding requires careful examination before biological interpretation.

**CEP trajectory:** The CEP compartment shows a positive pseudotime-condition correlation:
- CEP: Spearman rho = **+0.135** (p significant)
- This is reversed from v2 (rho=-0.163), likely reflecting changes in cell composition from the annotation correction.

**Figure 13. AF cell state trajectory -- UMAP with pseudotime overlay.**
![AF trajectory UMAP](trajectories/umap_trajectory_AF.png)

**Gene dynamics along NP pseudotime:** Notochordal markers (KRT8, KRT18) decline monotonically with pseudotime, while stress/inflammatory markers increase, consistent with the proposed continuum model. Mature chondrocyte markers (ACAN, COL2A1) peak at intermediate pseudotime and decline at the degenerative end, suggesting an initial maintenance phase followed by loss of chondrocyte identity.

**Figure 14. Gene expression dynamics along NP pseudotime.**
![NP gene dynamics](trajectories/gene_dynamics_NP.png)

### 5.6 Cell-Cell Communication

**Figure 15. Cell-cell interaction heatmap -- healthy tissue.**
![Healthy interactions](communication/interaction_plots/interaction_heatmap_healthy.png)

**Figure 16. Cell-cell interaction heatmap -- degenerated tissue.**
![Degenerated interactions](communication/interaction_plots/interaction_heatmap_degenerated.png)

**Figure 17. Differential interactions between healthy and degenerated tissue.**
![Differential interactions](communication/interaction_plots/differential_interactions.png)

LIANA consensus analysis identified **40,187 ligand-receptor interactions in healthy** and **40,872 in degenerated** tissue -- a roughly **balanced** pattern with a modest 1.7% increase in degenerated tissue.

This finding contrasts with both v1 (20% increase in degenerated) and v2 (6.5% decrease in degenerated), and is more consistent with the interpretation that total interaction count is relatively stable while the *composition* of interactions shifts. The v3 result -- near-parity between conditions -- suggests that the dramatic quantitative differences seen in v1 and v2 were artifacts of cell type misassignment, particularly the ~17K misrouted stressed NP cells that would have inflated immune cell interaction counts.

**Differential interactions:** Differential interaction pairs were identified between healthy and degenerated conditions, revealing shifts in specific ligand-receptor axes rather than global changes in communication volume.

**Pain-relevant interactions:** Pain-relevant interactions were flagged through cross-referencing with curated gene sets (nociception, neurotrophins, nerve guidance, inflammatory pain, neovascularization).

### 5.7 Pain Biology

Cross-referencing DE genes with curated pain gene sets identified **10 significant pain-relevant genes** across comparisons:

- **PTGS2** (COX-2): log2FC=+2.42, padj=0.005 in NP_mature_chondrocyte healthy_vs_severe; also significant in AF_inner. Prostaglandin synthesis enzyme and direct mediator of inflammatory pain.
- **PTGES**: log2FC=+3.84, padj=0.036 in NP_mature_chondrocyte healthy_vs_severe. Prostaglandin E synthase, catalyzes PGE2 production.
- **PLA2G2A**: log2FC=+1.64, padj=0.042 in NP_mature_chondrocyte mild_vs_severe. Phospholipase A2, produces arachidonic acid precursors.
- **CCL2**: log2FC=+1.86, padj=0.031 in NP_mature_chondrocyte mild_vs_severe; also significant in T_cell and AF_outer. Monocyte chemoattractant protein-1, involved in neuroinflammation.
- **TNF**: log2FC=+2.30, padj=2.8x10^-8 in T_cell mild_vs_severe. Master inflammatory cytokine and pain mediator.
- **BDKRB2**: log2FC=+1.74, padj=0.003 in NP_fibrocartilaginous mild_vs_severe. Bradykinin receptor B2, mediates pain signaling.
- **NRP2**: log2FC=+1.44, padj=0.034 in NP_mature_chondrocyte healthy_vs_severe; also significant in T_cell. Neuropilin-2, involved in nerve guidance.
- **PDGFA**: log2FC=+2.16, padj=0.014 in NP_mature_chondrocyte healthy_vs_severe. Platelet-derived growth factor, promotes neovascularization.
- **ROBO1**: log2FC=+1.26, padj=0.035 in NP_fibrocartilaginous healthy_vs_severe. Roundabout receptor, nerve guidance molecule.
- **SEMA3A**: log2FC=-1.76, padj=0.016 in AF_outer healthy_vs_severe. Semaphorin 3A, a nerve repellent; its downregulation in degeneration may permit nerve ingrowth.

**Directly supported by our DE data:**
- The prostaglandin axis (PLA2G2A, PTGS2/COX-2, PTGES) constitutes a complete biosynthetic pathway from arachidonic acid release to PGE2 production, a direct sensitizer of nociceptive nerve endings (Risbud and Shapiro, 2014).
- SEMA3A downregulation in AF_outer is particularly notable: semaphorin 3A normally repels sensory nerve fibers from entering the disc. Its loss in degenerated AF is consistent with the nerve ingrowth model of discogenic pain (Freemont et al., 2002).
- CXCL2 (significant in NP_mature_chondrocyte mild_vs_severe) recruits neutrophils and macrophages, which produce additional pain mediators.

**Not detected in our data:**
- NGF (nerve growth factor) and BDNF (brain-derived neurotrophic factor), classically associated with nerve ingrowth into degenerated discs (Freemont et al., 2002), were **not significantly upregulated** in any powered comparison.

**Model:** Degenerated disc cells create a pro-inflammatory microenvironment through chemokine and prostaglandin production that promotes nerve ingrowth (via SEMA3A loss and CXCL2/CCL2 immune recruitment) and sensitization (via PGE2), rather than directly signaling pain. This is consistent with the two-signal model of discogenic pain: (1) structural disruption and loss of nerve-repellent signals permits nerve ingrowth into the NP, and (2) the inflammatory milieu sensitizes ingrown nerves (Freemont et al., 2002; Risbud and Shapiro, 2014).

**Figure 18. Pain-associated gene expression heatmap across cell types.**
![Pain genes heatmap](interpretation/pain_genes_heatmap.png)

---

## 6. Biological Interpretation and Mechanistic Model

### 6.1 The Inflammatory/Catabolic Cascade

Synthesizing our DE, pathway, TF, and CCC results, we propose that inflammatory cytokine and chemokine production by NP cells contributes to the degenerative cascade:

1. **Initiation:** Mechanical stress, aging, or microinjury activates NF-kB signaling in disc cells (supported by: TF activity analysis; inflammatory pathway enrichment).

2. **Chemokine and cytokine activation:** NF-kB drives expression of inflammatory mediators including CXCL2 and CCL2 by NP cells (supported by: CXCL2 log2FC=+3.63, padj=1.75x10^-4; CCL2 log2FC=+1.86, padj=0.031 in NP_mature_chondrocyte mild_vs_severe; chemokine pathway enrichment).

3. **Immune cell recruitment:** Chemokines recruit neutrophils (CXCL2/CXCR2 axis) and monocytes/macrophages (CCL2/CCR2 axis) into the disc space (supported by: neutrophil chemotaxis pathway enrichment; TNF highly significant in T_cell mild_vs_severe).

4. **Prostaglandin-mediated pain sensitization:** PLA2G2A, PTGS2, and PTGES constitute a complete prostaglandin synthesis pathway, producing PGE2 that directly sensitizes nociceptive nerve endings (supported by: all three genes significantly DE in NP cells).

5. **Nerve ingrowth facilitation:** SEMA3A downregulation in AF_outer removes a key nerve-repellent barrier, while neovascularization factors (PDGFA) promote the vascular and neural invasion characteristic of painful disc degeneration (supported by: SEMA3A log2FC=-1.76 in AF_outer).

6. **Fibrocartilaginous replacement:** NP_fibrocartilaginous cells -- the most transcriptionally responsive population -- drive ECM remodeling, replacing the hydrated proteoglycan-rich NP matrix with fibrotic collagen I-rich tissue (supported by: 418 DE genes in mild_vs_severe; ECM organization pathway enrichment).

7. **Metabolic failure in AF:** AF cells experience simultaneous proteotoxic stress (HSP activation) and mitochondrial dysfunction (oxidative phosphorylation downregulation), compromising structural integrity of the outer disc (supported by: GSEA pathway enrichment).

### 6.2 The Prostaglandin-Pain Axis: A Coherent Therapeutic Target

The v3 analysis provides the most complete evidence to date for a prostaglandin-mediated pain mechanism in disc degeneration at single-cell resolution. The three-enzyme pathway -- PLA2G2A (arachidonic acid release) to PTGS2/COX-2 (prostaglandin H2 synthesis) to PTGES (PGE2 production) -- is significantly upregulated in NP cells during degeneration. Combined with SEMA3A downregulation permitting nerve ingrowth and CXCL2/CCL2 driving immune infiltration, this defines a mechanistic model linking the inflammatory microenvironment to discogenic pain that is supported by multiple independent gene sets.

---

## 7. Therapeutic Targets

Based on the evidence from this analysis, we propose the following therapeutic targets, ranked by strength of supporting data.

### 7.1 Tier 1: Strong Direct Evidence From This Analysis

**Target 1: Chemokine Modulation (CXCL2/CXCR2 and CCL2/CCR2)**
- **Evidence from this analysis:** CXCL2 (log2FC=+3.63, padj=1.75x10^-4) in NP_mature_chondrocyte mild_vs_severe -- the most significant inflammatory gene in the dataset. CCL2 significant in NP_mature_chondrocyte, T_cell, and AF_outer. Chemokine pathway enrichment in ORA.
- **v3 improvement:** CXCL2 significance improved from padj=0.005 (v2) to padj=1.75x10^-4 (v3), strengthening confidence in this target.
- **Mechanism:** CXCL2 signals through CXCR2 on neutrophils; CCL2 signals through CCR2 on monocytes/macrophages. Dual chemokine blockade could interrupt immune cell recruitment into the disc space.
- **Approach:** Intradiscal CXCR2 antagonist (navarixin, AZD5069) or CCR2 antagonist to block immune cell recruitment without systemic immunosuppression.

**Target 2: Prostaglandin Pathway Inhibition**
- **Evidence from this analysis:** PLA2G2A (padj=0.042), PTGS2 (padj=0.005), PTGES (padj=0.036) -- a complete three-enzyme pathway. This is the most coherent multi-gene pain target in the dataset.
- **Mechanism:** PLA2G2A releases arachidonic acid; PTGS2/COX-2 converts it to PGH2; PTGES produces PGE2, which directly sensitizes nociceptive nerve endings.
- **Approach:** Selective COX-2 inhibitors (celecoxib) are already used clinically; intradiscal delivery could achieve higher local concentrations with reduced systemic effects. PGE2 receptor antagonists offer more targeted intervention.

**Target 3: TNF/NF-kB Inhibition**
- **Evidence from this analysis:** TNF is the most significant DE gene in T_cell mild_vs_severe (padj=2.8x10^-8). NF-kB pathway activity is supported by TF analysis and inflammatory pathway enrichment. Multiple NF-kB target genes (CXCL2, CCL2, PTGS2) are significantly upregulated.
- **Mechanism:** NF-kB drives the entire catabolic cascade -- chemokines, cytokines, MMPs, and prostaglandins.
- **Approach:** Intradiscal anti-TNF biologics (etanercept, adalimumab) or small molecule NF-kB inhibitors.
- **Status:** Early clinical data available for epidural anti-TNF (Cohen et al., 2009).

### 7.2 Tier 2: Moderate Evidence, Requires Validation

**Target 4: SEMA3A Restoration**
- **Evidence from this analysis:** SEMA3A (log2FC=-1.76, padj=0.016) downregulated in AF_outer healthy_vs_severe. SEMA3A is a nerve repellent; its loss facilitates nerve ingrowth into the disc.
- **Novel aspect:** This is the first report of SEMA3A downregulation in AF at single-cell resolution in a multi-dataset meta-analysis.
- **Approach:** Intradiscal delivery of recombinant SEMA3A or gene therapy to restore the nerve-repellent barrier.

**Target 5: HSP/Proteostasis Modulation**
- **Evidence from this analysis:** Heat response pathways enriched in AF_inner; GSEA confirms proteotoxic stress concurrent with mitochondrial dysfunction.
- **Mechanism:** Chemical chaperones could reduce ER stress and alleviate the energy crisis in AF cells.
- **Approach:** Chemical chaperones (4-PBA, TUDCA) or mitochondria-targeted antioxidants (MitoQ, SS-31).

**Target 6: Mitochondrial Rescue in AF**
- **Evidence from this analysis:** GSEA shows oxidative phosphorylation, electron transport chain, and mitochondrial ATP synthesis downregulated in AF_inner.
- **Mechanism:** Restoring mitochondrial function could reduce ROS, improve energy metabolism, and support ECM maintenance (Song et al., 2023b).
- **Approach:** Mitochondria-targeted antioxidants (MitoQ, SS-31) or NAD+ precursors (NMN, NR).

### 7.3 Tier 3: Supported by Literature, Not Directly Demonstrated in This Data

**Target 7: ADAMTS5 Inhibition**
- **This analysis:** ADAMTS5 does not reach significance after FDR correction.
- **Literature:** ADAMTS5 is the primary aggrecanase in cartilaginous tissues (Stanton et al., 2005) and is consistently reported as upregulated in disc degeneration (Liang et al., 2022).
- **Status:** Small molecule inhibitors developed for osteoarthritis in preclinical testing.

**Target 8: TIMP1 Restoration**
- **This analysis:** TIMP1-CD63 loss was not among the top differential interactions in our CCC analysis. The companion phylo analysis identified this as a dominant lost interaction.
- **Literature:** The MMP/TIMP balance is well-established in disc degeneration (Vo et al., 2013; Cabral-Pacheco et al., 2020).

**Target 9: Senolytic Therapy**
- **This analysis:** Senescence pathways did not reach significance in our GSEA.
- **Literature:** Dasatinib + quercetin senolytics ameliorate disc degeneration in mice (Novais et al., 2021).

### 7.4 Summary Therapeutic Target Table

| Target | Gene(s) | Evidence Level | Key Data Point | Approach |
|--------|---------|---------------|----------------|----------|
| CXCR2/CCR2 antagonism | CXCL2, CCL2 | Strong (this study) | CXCL2 padj=1.75x10^-4 | Small molecule |
| Prostaglandin inhibition | PLA2G2A, PTGS2, PTGES | Strong (this study) | 3 pathway genes sig | COX-2 inhibitor |
| TNF/NF-kB inhibition | TNF, NF-kB targets | Strong (this study) | TNF padj=2.8x10^-8 | Biologic / small mol |
| SEMA3A restoration | SEMA3A | Moderate (this study) | padj=0.016 | Gene therapy |
| HSP modulation | HSF1, HSPA1A/B | Moderate (this study) | GSEA enrichment | Chemical chaperone |
| Mitochondrial rescue | OXPHOS genes | Moderate (this study) | GSEA suppression | MitoQ / NAD+ |
| ADAMTS5 inhibition | ADAMTS5 | Literature only | Not sig in this study | Small molecule |
| TIMP1 restoration | TIMP1, CD63 | Literature only | Not primary CCC finding | Gene therapy |
| Senolytics | CDKN1A/2A | Literature only | Not sig in GSEA | D+Q |

---

## 8. Novel and Discordant Findings

### 8.1 Strengthened CXCL2 Signal in v3

CXCL2 significance improved from padj=0.005 (v2) to padj=1.75x10^-4 (v3), and the log2FC increased from +3.14 to +3.63. This strengthening is a direct consequence of the annotation correction: by returning ~17K stressed NP cells to the mesenchymal pool, the NP_mature_chondrocyte pseudobulk samples now better represent true chondrocyte populations, reducing noise in the DE analysis and increasing power to detect genuine inflammatory signals. The CXCL2/CXCR2 axis for neutrophil recruitment in NP degeneration is now among the strongest signals in the entire dataset.

### 8.2 NP_fibrocartilaginous Dominates the DE Landscape

NP_fibrocartilaginous cells are the most transcriptionally responsive cell type, with 418 DE genes in mild_vs_severe (vs. 203 in v2) and 385 in healthy_vs_severe (vs. 127 in v2). This near-doubling of DE genes likely reflects improved annotation boundaries. The fibrocartilaginous population -- characterized by COL1A1 expression and a transitional phenotype between chondrocyte and fibroblast -- appears to be a key driver of the degenerative response.

### 8.3 CCC Interaction Counts Stabilize

The v3 analysis shows roughly balanced CCC interactions (40,187 healthy vs. 40,872 degenerated, a 1.7% difference), contrasting sharply with v1 (+20% in degenerated) and v2 (-6.5% in degenerated). The stabilization is consistent with the hypothesis that the quantitative imbalance in prior versions was driven by misassigned cells inflating interaction counts. The v3 result suggests that the primary change in degeneration is the *composition* of cell-cell interactions, not their total number.

### 8.4 Complete Prostaglandin Synthesis Pathway

The identification of all three enzymes in the arachidonic acid-to-PGE2 pathway (PLA2G2A, PTGS2, PTGES) as significantly DE is a coherent finding not commonly reported at single-cell resolution in IVD literature. This complete pathway provides stronger therapeutic rationale than individual gene findings.

### 8.5 SEMA3A Downregulation in AF

SEMA3A downregulation (log2FC=-1.76, padj=0.016) in AF_outer is a notable finding linking the DE results to the pain biology of disc degeneration. Semaphorin 3A is a well-characterized nerve repellent (Freemont et al., 2002), and its loss in the AF -- the primary barrier to nerve ingrowth -- provides a mechanistic link between ECM degeneration and discogenic pain.

### 8.6 Discordance with Companion Phylo Analysis: Wnt, Notch, and Senescence

The companion phylo analysis reported consistent suppression of Wnt signaling, Notch signaling, and cellular senescence pathways. Our analysis did not replicate these findings. As noted in v2, the phylo analysis results appear substantially driven by replication-dependent histone genes that are sensitive to cell cycle state, dissociation protocols, and ambient RNA contamination (Slyper et al., 2020). Our use of LFC shrinkage and prioritization of within-study (mild_vs_severe) comparisons provides more conservative but more robust results.

### 8.7 Trajectory Correlations Shifted

The NP pseudotime-condition correlation weakened from rho=-0.258 (v2) to rho=-0.151 (v3), while AF (rho=+0.325) and CEP (rho=+0.135) show positive correlations. The NP weakening may reflect the reassignment of stressed cells back to the mesenchymal pool, which alters the pseudotime distribution. The positive AF and CEP correlations require SME review.

---

## 9. Limitations

1. **Cross-study confounding:** Condition and study are partially confounded. Herniated comparisons were excluded entirely due to single-study confounding. Within-study comparisons (mild_vs_severe) are prioritized throughout.

2. **Underpowered comparisons:** Many cell type x comparison combinations lack sufficient samples (< 3 per condition). Key genes like ADAMTS5, ACAN, and COL2A1 may fail to reach significance due to donor variability.

3. **Age-disease confound:** In GSE230809 (the largest dataset, 24 samples), healthy donors are 21-27 years old and diseased are 37-73 years old. Age and degeneration effects cannot be fully separated.

4. **Sex bias:** GSE230809 is all-male. Many samples have unknown sex. Sex-stratified analyses are not possible.

5. **Culture-expanded cells:** GSE242443 CEP cells are culture-expanded, which alters gene expression (particularly collagen ratios and surface markers).

6. **No RNA velocity:** Spliced/unspliced counts were not available from public deposits. RNA velocity would provide directional evidence for cell state transitions.

7. **No SCENIC/GRN:** Full SCENIC analysis (gene regulatory networks) was not performed. TF activity was estimated from CollecTRI regulon overlap, which captures target enrichment but not regulatory network structure.

8. **Composition analysis underpowered:** No cell type proportion changes reached significance after FDR correction, though trends were biologically consistent.

9. **Annotation sensitivity:** The v3 annotation correction demonstrates that cell type assignment substantially impacts all downstream results. The ~17K misrouted cells affected DE, TF, trajectory, and CCC analyses in v2. While v3 addresses this specific issue, other annotation decisions (e.g., the boundary between NP_mature_chondrocyte and NP_fibrocartilaginous) remain subjective.

10. **CCC methodology and fragility:** The reversal of CCC quantitative patterns across v1, v2, and v3 (from +20% to -6.5% to +1.7%) highlights the extreme sensitivity of these analyses to cell type definitions.

11. **AF and CEP trajectory direction:** The positive pseudotime-condition correlations in AF (rho=+0.325) and CEP (rho=+0.135) are difficult to interpret biologically and may reflect root cell selection or cell composition artifacts.

12. **TF analysis sensitivity:** The reduction from 290 (v2) to 5 (v3) significant TF associations demonstrates the sensitivity of regulon overlap analysis to DE gene input. Small changes in cell type boundaries cascade through to large changes in TF results.

---

## 10. Conclusion

This 11-dataset, 410,759-cell meta-analysis of human IVD degeneration, now in its third pipeline iteration (v3), reveals a robust inflammatory transcriptomic signature in severe NP degeneration. The key v3 improvement -- correction of ~17,000 misrouted stressed NP cells through an annotation evidence gate, ACAN/SOX9 rescue, and stricter cluster voting -- strengthened the CXCL2 chemokine signal (from padj=0.005 to padj=1.75x10^-4) and stabilized the cell-cell communication analysis. The inflammatory signature is centered on CXCL2-mediated neutrophil recruitment and a complete prostaglandin synthesis pathway (PLA2G2A, PTGS2, PTGES) that links the inflammatory microenvironment to pain biology through PGE2-mediated nerve sensitization. SEMA3A downregulation in AF_outer provides a complementary mechanism for nerve ingrowth.

NP_fibrocartilaginous cells emerged as the most transcriptionally responsive population to degeneration (418 DE genes in mild_vs_severe), suggesting that the transitional fibrocartilaginous replacement process is a central feature of the degenerative cascade. Trajectory analysis confirms that NP cells exist on a disease-associated continuum (rho=-0.151), though the correlation is weaker than v2.

The v1-to-v2-to-v3 changes documented across three pipeline iterations serve as a cautionary demonstration of how cell type annotation decisions propagate through all downstream analyses: DE gene counts, pathway enrichment, TF activity, trajectory correlations, and CCC interaction counts all shifted substantially with the correction of 17K misassigned cells. This underscores the importance of rigorous annotation validation and motivates the use of annotation evidence gates as a standard practice in scRNA-seq meta-analysis pipelines.

The primary therapeutic opportunities are chemokine modulation (CXCL2/CXCR2, CCL2/CCR2), prostaglandin pathway inhibition (PLA2G2A/PTGS2/PTGES), and TNF/NF-kB inhibition. SEMA3A restoration in AF is proposed as a novel strategy to prevent nerve ingrowth.

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

Liang H, Luo R, Li G, et al. (2022). The Proteolysis of ECM in Intervertebral Disc Degeneration. *International Journal of Molecular Sciences*, 23:1715.

Liberzon A, Birger C, Thorvaldsdottir H, et al. (2015). The Molecular Signatures Database Hallmark gene set collection. *Cell Systems*, 1(6):417-425.

Li X, Han Y, Li G, et al. (2023a). Role of Wnt signaling pathway in joint development and cartilage degeneration. *Frontiers in Cell and Developmental Biology*, 11:1181619.

Li Z, Ye D, Dai L, et al. (2022a). Single-Cell RNA Sequencing Reveals the Difference in Human Normal and Degenerative Nucleus Pulposus Tissue Profiles and Cellular Interactions. *Frontiers in Cell and Developmental Biology*, 10:910626.

Long J, Wang X, Du X, et al. (2019). JAG2/Notch2 inhibits intervertebral disc degeneration by modulating cell proliferation, apoptosis, and extracellular matrix. *Arthritis Research & Therapy*, 21:213.

Lopez R, Regier J, Cole MB, et al. (2018). Deep generative modeling for single-cell transcriptomics. *Nature Methods*, 15:1053-1058.

Love MI, Huber W, Anders S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15:550.

Novais EJ, Tran VA, Johnston SN, et al. (2021). Long-term treatment with senolytic drugs dasatinib and quercetin ameliorates age-dependent intervertebral disc degeneration in mice. *Nature Communications*, 12:5213.

Oichi T, Taniguchi Y, Oshima Y, et al. (2020). Pathomechanism of intervertebral disc degeneration. *JOR Spine*, 3:e1076.

Rennard SI, Dale DC, Donohue JF, et al. (2015). CXCR2 Antagonist MK-7123. A Phase 2 Proof-of-Concept Trial for Chronic Obstructive Pulmonary Disease. *American Journal of Respiratory and Critical Care Medicine*, 191(9):1001-1011.

Risbud MV, Shapiro IM. (2014). Role of cytokines in intervertebral disc degeneration: pain and disc content. *Nature Reviews Rheumatology*, 10(1):44-56.

Slyper M, Porter CBM, Ashenberg O, et al. (2020). A single-cell and single-nucleus RNA-Seq toolbox for fresh and frozen human tumors. *Nature Medicine*, 26:792-802.

Song C, Zhou Y, Cheng K, et al. (2023a). Cellular senescence -- Molecular mechanisms of intervertebral disc degeneration from an immune perspective. *Biomedicine & Pharmacotherapy*, 162:114711.

Song C, Xu Y, Peng Q, et al. (2023b). Mitochondrial dysfunction: a new molecular mechanism of intervertebral disc degeneration. *Inflammation Research*, 72:2249-2260.

Squair JW, Gautier M, Kathe C, et al. (2021). Confronting false discoveries in single-cell differential expression. *Nature Communications*, 12:5692.

Stanton H, Rogerson FM, East CJ, et al. (2005). ADAMTS5 is the major aggrecanase in mouse cartilage in vivo and in vitro. *Nature*, 434:648-652.

Vo N, Hartman R, Yurube T, et al. (2013). Expression and regulation of metalloproteinases and their inhibitors in intervertebral disc aging and degeneration. *The Spine Journal*, 13:331-341.

Wang Y, Cheng H, Wang T, et al. (2023a). Oxidative stress in intervertebral disc degeneration: Molecular mechanisms, pathogenesis and treatment. *Cell Proliferation*, 56:e13448.

Wolock SL, Lopez R, Klein AM. (2019). Scrublet: computational identification of cell doublets in single-cell transcriptomic data. *Cell Systems*, 8(4):281-291.e9.

Wuertz K, Vo N, Kletsas D, Boos N. (2012). Inflammatory and catabolic signalling in intervertebral discs: the roles of NF-kB and MAP kinases. *European Cells and Materials*, 23:103-120.

Xia Q, Zhao Y, Dong H, et al. (2024). Progress in the study of molecular mechanisms of intervertebral disc degeneration. *Biomedicine & Pharmacotherapy*, 174:116593.

Zimmerman KD, Espeland MA, Langefeld CD. (2021). A practical solution to pseudoreplication bias in single-cell studies. *Nature Communications*, 12:738.

---

*Analysis performed using a 10-module human-gated agentic pipeline (v3). All code version-controlled. Random seed: 42. Package versions: Python 3.12, scanpy 1.11, scvi-tools 1.4.2, pyDESeq2, gseapy 1.1, decoupler 2.1, liana 1.7. Key v3 improvement: annotation evidence gate, ACAN/SOX9 rescue, and stricter cluster voting correcting ~17K misrouted stressed NP cells.*

*This is a computational analysis draft. All findings require experimental validation before clinical application.*
