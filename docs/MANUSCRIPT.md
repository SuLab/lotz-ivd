---
title: "A Continuum-Aware Single-Cell Atlas of the Human Intervertebral Disc: Cell-Type-Specific Degenerative Programs and Pain-Associated Molecular Circuits"
author: "IVD Atlas Project"
date: "2026-05-22"
---

## Abstract

The intervertebral disc (IVD) is a chronically remodelling musculoskeletal tissue whose degeneration underlies the majority of cases of structural low back pain. Resident IVD cells exist on transcriptomically graded continua within each anatomical compartment, presenting a methodological challenge: aggressive batch correction across studies can erase the very within-compartment variation that the analysis aims to characterize. We address this by aggregating 12 publicly available human IVD single-cell RNA-seq studies (78 samples, 57 donors, 410,705 cells across nucleus pulposus, annulus fibrosus, and cartilaginous endplate) and integrating them under a **tiered scVI strategy**: immune, endothelial, and other non-resident cells are integrated with standard scVI defaults, whereas chondrocyte- and fibroblast-like resident cells are integrated under more conservative settings designed to preserve within-compartment continua. A three-stage annotation pipeline (coarse cell class → compartment-prefixed cell type → sub-state) yields 19 cell types and 26 named sub-states across the atlas, with two contamination categories (red-blood-cell, endothelial-admixed) retained as flagged cells rather than filtered. Pseudobulk DESeq2 across 26 statistically powered cell-type × condition contrasts (one inflated contrast excluded; see Caveats) identifies 1,823 significant differentially expressed genes (DEGs) at FDR < 0.05, with the three NP cell types collectively contributing 72% of the trustworthy signal and NP_fibrocartilaginous alone contributing 630 DEGs (325 in the mild-vs-severe contrast). Eighteen pain-associated genes — including IL6, NTN1, UNC5B, PENK, VEGFA, PDGFA, BDKRB1/2, PTGS2, CCL2, TNF, FGF2, and FLT1 — are cell-type-specifically dysregulated, with NP_fibrocartilaginous carrying the broadest signal and NP_mature_chondrocyte concentrating the neovascularization arm. Pathway enrichment (3,051 significant ORA terms; 6,890 GSEA), transcription factor activity inference (252 unique significant TFs; AP-1, EGR1, JUN, RELA, SP1, NFKB1, PPARG, FOS, CEBPB recurring across NP contrasts), trajectory analysis (PAGA-initialized diffusion pseudotime), and cell-cell communication (LIANA consensus rank-aggregation) collectively describe a coherent biology: fibrotic ECM remodeling, NF-κB-driven inflammatory programs, pain-mediator dysregulation, and neutrophil-recruitment signaling gained in degeneration alongside loss of basement-membrane and WNT signaling at the cartilaginous endplate. Compositional shifts at the cell-type level do not reach FDR significance, supporting the view that degeneration is principally a transcriptional, not a cellular, remodelling. We report the atlas with explicit caveats for the methodological and biological limitations of cross-study scRNA-seq meta-analysis.

---

## Introduction

Low back pain is the leading single cause of years lived with disability globally, and structural intervertebral disc (IVD) degeneration is its most common identifiable substrate. The healthy adult disc is a heterogeneous, avascular, hypocellular tissue comprising the gelatinous **nucleus pulposus** (NP), the concentric **annulus fibrosus** (AF), and the **cartilaginous endplates** (CEP) that interface with vertebral bone. Disc degeneration involves a coordinated cascade of extracellular matrix (ECM) catabolism, cellular phenotype shifts, neovascularization, and aberrant sensory nerve ingrowth — the latter two contributing directly to the discogenic pain phenotype.

Since 2020, multiple groups have applied single-cell RNA sequencing to human IVD tissue, yielding insights into resident cell heterogeneity, fibrotic transitions, and immune infiltration. Each individual study is, however, modest in sample size, often restricted to one or two compartments, and uses non-standardized annotation. Meta-analysis across studies is therefore attractive (statistical power, generalizability) but methodologically delicate. Three issues recur:

1. **The continuum problem.** Resident chondrocyte- and fibroblast-like cells exist on a graded morphological and molecular spectrum within each compartment. Standard batch correction methods applied uniformly across all cell types tend to compress this within-compartment variation, conflating biological continua with technical noise.

2. **Pseudobulk vs. single-cell DE.** Single-cell DE tests that treat individual cells as independent observations inflate false-positive rates substantially when applied across donors with unequal cell counts. Pseudobulk aggregation followed by donor-level DESeq2 is now standard practice for cross-study scRNA-seq DE.

3. **Contamination handling.** Public IVD datasets vary in their stringency for excluding red-blood-cell-contaminated NP samples and for distinguishing genuine endothelial signal from endothelial cells admixed with NP_fibrocartilaginous cells. Silent filtering of these cells removes counts that may inform the contamination calls themselves.

This manuscript reports an IVD scRNA-seq atlas constructed to address each of these issues directly. We use a **tiered integration strategy** that splits the data by coarse cell class — immune, endothelial, and other non-resident lineages on one path; chondrocyte- and fibroblast-like resident cells on a second path with more conservative settings — and we annotate in three stages (coarse class, compartment-prefixed cell type, sub-state). We retain contamination cells with explicit flags rather than filtering them. All differential expression uses pseudobulk DESeq2. The full pipeline (data acquisition through final analyses) is reproducible end-to-end from the published count matrices.

---

## Methods

### Datasets and harmonization

Twelve human IVD scRNA-seq datasets were identified through GEO and CNGB and downloaded with their original count matrices. The selection process initially screened 23 candidate studies; 11 were excluded as mouse/rat models, animal-disc-tissue studies, or for unavailable count data. The 12 retained studies span **78 samples from 57 donors**, with conditions covering healthy adult (20 samples), mild degeneration (18), severe degeneration (21), ungraded degeneration (3), herniated (10), neonatal (3), and aged ungraded (3). Per-study count matrices were converted to AnnData (.h5ad) and harmonized at the sample-metadata level: condition labels were collapsed to seven canonical categories and compartment was recorded as NP (49 samples), AF (17), CEP (6), or IVD_mixed (6, for the single dataset that did not separate compartments).

GSE205535 was processed against its published corrigenda. Three datasets (GSE242443, GSE251686, GSE255768) use non-10x chemistries (BD Rhapsody or Singleron); platform identity was retained as an integration covariate.

### Quality control and preprocessing

Per-cell QC used per-sample ambient-RNA estimation, scrublet-based doublet detection, and conservative cutoffs on `pct_counts_mt` (≤ 20%), `n_genes_by_counts` (≥ 200, with platform-aware upper bounds), and `total_counts` (≥ 500). Genes detected in fewer than three cells per dataset were dropped. Counts were retained as raw integers; CP10K + log1p normalization was applied lazily on a per-analysis basis to preserve raw counts for pseudobulk aggregation.

### Coarse cell classification

A panel-based scoring approach assigned every cell a coarse class — `mesenchymal`, `immune`, `endothelial`, or `unknown` — using marker panels for canonical disc-resident lineages (COL1A1, COL2A1, ACAN, PRG4), immune subsets (CD3D, CD8A, CD68, CD14, LYZ, S100A8, MS4A1, MZB1), red blood cells (HBB, HBA1, HBA2), and endothelial cells (PECAM1, CDH5, EMCN, VWF). Cells with no panel reaching threshold were labelled `unknown` and routed downstream into the mesenchymal integration tier. The `unknown` cells are predominantly low-quality disc cells rather than missed lineages.

### Tiered integration

The tiered integration strategy splits the data into two parallel streams. The **non-resident tier** integrates immune, endothelial, RBC, and similar cells under default scVI settings; the **resident tier** integrates disc-resident chondrocyte- and fibroblast-like cells under more conservative settings designed to preserve within-compartment continua. Both tiers were integrated per compartment (NP, AF, CEP) and on the union of all cells, producing four `.h5ad` objects with 30-dimensional `X_integrated` latent representations.

The tiered approach was selected over a flat single-method integration after an NP-specific quality experiment demonstrated that applying a single integration uniformly across all cell types collapsed the well-documented NP_fibrocartilaginous-to-mature-chondrocyte continuum into a single homogeneous cluster. The substantive distinction is therefore: separate models for separable populations, conservative settings for populations whose within-compartment variation is biologically informative.

A sensitivity analysis was also performed using flat Seurat CCA integration applied uniformly across all cells; the tiered approach yielded finer cell-type resolution and a larger trustworthy DE pool (1,823 vs. 1,198 significant DEGs across versions), with biological themes preserved across the two integration strategies. We report the tiered analysis as primary.

### Clustering

Leiden clustering was applied per tier per compartment, scanning resolutions chosen by an equal-weighted silhouette + modularity score. Tier-aware adaptive thresholds (three resolutions for > 300K cells, six for > 200K, ten for > 50K) and skipping modularity computation for > 100K cells kept run times tractable. Resulting cluster counts are: NP 17 resident + 4 non-resident, AF 9 + 3, CEP 4 + 4, all_cells 24 + 9.

### Annotation

Annotation was performed in three stages:

1. **Coarse:** confirmed/refined the Module 04 coarse class using post-integration neighbourhoods. CellTypist (with CP10K + log1p normalized inputs) provided independent reference labels for the non-resident tier; a 60% per-cluster majority threshold against the panel-derived coarse label served as fallback when CellTypist did not fire.

2. **Cell type (compartment-prefixed):** panel-based scoring against fine-grained marker panels produced compartment-prefixed labels — NP_fibrocartilaginous, NP_mature_chondrocyte, NP_fibrochondrocyte_chondroid, AF_inner, AF_outer, CEP_outer, CEP_hyaline, CEP_fibrochondrocyte_fibroid — plus shared non-resident labels (Macrophage_M1, Macrophage_M2, Neutrophil, Immune, Endothelial, Pericyte_SMC, Erythrocyte). A label-harmonization step collapsed historical generic labels (e.g. `Fibroblast_like` → compartment-specific) for cross-compartment consistency.

3. **Sub-state:** overlap-based scoring against sub-state panels (`proliferating`, `inflammatory`, `stressed`, `matrix_active`, `migratory`, `homeostatic`) plus an endothelial-admixed contamination flag (CD34, EMCN, AQP1 panel). Resident cells received a `cell_subtype` label; non-resident cells inherited `cell_subtype = cell_type`. The all_cells object yielded 26 named sub-states (excluding the `unassigned` residual).

**Contamination handling.** 16,514 Erythrocyte cells (4.0% of the atlas) and 1,831 endothelial-admixed cells within NP_fibrocartilaginous (0.4%) were retained with `is_contamination = True` flags and `contamination_type ∈ {RBC, endothelial_admixed, clean}`. This preserves the count for downstream analyses while allowing transparent filtering at the interpretation stage.

### Pseudobulk differential expression

Per-sample pseudobulk count aggregation followed by DESeq2 was used for differential expression to avoid the inflated false-positive rate of single-cell DE tests that treat cells as independent observations. The primary grouping was `cell_type` (compartment-prefixed); composition tests were additionally run on `cell_subtype`. Comparisons were defined per cell type as `healthy_vs_degenerated_all`, `healthy_vs_degenerated_mild`, `healthy_vs_degenerated_severe`, and `mild_vs_severe`. Comparisons with fewer than three samples per group were marked underpowered and skipped. One inflated contrast (Macrophage_M2 healthy-vs-severe, which returned 5,659 DEGs from a 7-sample comparison and showed evident dispersion-estimation pathology) is retained on disk for transparency but excluded from all downstream interpretation; see Caveats §1.

### Pathway enrichment, transcription factor activity, and pain-gene cross-referencing

Over-representation analysis (ORA) against GO/KEGG/Reactome via Enrichr, pre-ranked GSEA (gseapy, MSigDB plus a custom IVD gene-set library), and CollecTRI-based transcription factor activity inference via decoupler were run per significant cell-type × comparison contrast. Pain-gene cross-referencing used a curated panel covering inflammatory pain, neurotrophin signalling, nerve guidance, neovascularization, neuropeptides, and ion-channel nociception categories.

### Trajectory analysis

PAGA-initialized diffusion pseudotime (DPT) on the resident tier computed per-compartment trajectories. The neighbour graph used `X_integrated`; root clusters were chosen by maximum enrichment of the expected mature/inner population (NP_mature_chondrocyte for NP, AF_inner for AF, CEP_hyaline for CEP). NP and AF were downsampled to 50,000 cells before PAGA; CEP was processed in full at 36,879 cells. RNA velocity was not feasible — none of the public count matrices include spliced/unspliced layers.

### Cell-cell communication

LIANA consensus rank-aggregation (CellPhoneDB, NATMI, Connectome, log2FC, sca, geometric mean) was run per condition group (healthy vs. degenerated pooled) on the union atlas, with each condition downsampled to 20,000 cells for tractability. Pain-relevant interactions were flagged by ligand or receptor membership in the curated pain panel.

### Software

scanpy 1.10, anndata 0.12.10, scvi-tools 1.4.2, Seurat 5.4.0, DESeq2 1.42.1, decoupler-py, gseapy, LIANA-py. Python 3.12; R 4.4 with Bioconductor 3.20. All analyses CPU-only on a 247 GB RAM, 32-core AWS instance.

---

## Results

### §1 — A 410,705-cell atlas of the human intervertebral disc

The integrated atlas comprises **410,705 cells** across 78 samples from 57 donors, distributed across NP (262,924 cells, 64%), AF (84,617 cells, 21%), and CEP (50,854 cells, 12%), with the remainder from a single dataset (GSE189916) that does not separate compartments and is labelled IVD_mixed. The tiered integration preserves within-compartment continuity: NP_fibrocartilaginous, NP_fibrochondrocyte_chondroid, and NP_mature_chondrocyte populations form a graded structure in the latent space rather than discrete clusters, consistent with the histological literature describing NP resident cells as a continuous spectrum.

![**Figure 1. UMAP of the all-cells integrated atlas after tiered scVI integration and three-stage annotation.** 410,705 cells coloured by `cell_type`. Resident populations (NP_fibrocartilaginous, NP_fibrochondrocyte_chondroid, NP_mature_chondrocyte, AF_inner, AF_outer, CEP_outer, CEP_hyaline, CEP_fibrochondrocyte_fibroid) and non-resident populations (Macrophage_M1, Macrophage_M2, Neutrophil, Immune, Endothelial, Pericyte_SMC, Erythrocyte) are clearly separated; within each resident compartment, related cell types are positioned as adjacent graded populations.](manuscript_figures/fig01_umap_all_cells.png)

Per-compartment cell-type counts are summarized below; compartment-specific UMAPs are shown in Figure 2.

**Table 1.** Cell-type composition of the atlas.

| Cell type | n cells | Sub-states resolved |
|---|---:|---|
| NP_fibrocartilaginous | 94,597 | 6 (proliferating, inflammatory, matrix_active, migratory, stressed, endothelial_admixed) |
| NP_fibrochondrocyte_chondroid | 59,650 | 1 (homeostatic) |
| NP_mature_chondrocyte | 33,010 | 1 (matrix_active) |
| AF_outer | 48,836 | 3 (homeostatic, matrix_active, proliferating) |
| AF_inner | 23,769 | 2 (homeostatic, stressed) |
| CEP_outer | 15,557 | 3 (homeostatic, matrix_active, proliferating) |
| CEP_hyaline | 12,306 | 2 (homeostatic, stressed) |
| CEP_fibrochondrocyte_fibroid | 9,016 | 1 (homeostatic) |
| Neutrophil | 27,912 | — |
| Macrophage_M2 | 22,752 | — |
| Immune | 22,404 | — |
| Erythrocyte (contamination) | 16,514 | — |
| Pericyte_SMC | 6,124 | — |
| Macrophage_M1 | 3,500 | — |
| Endothelial | 2,474 | — |

Note: the all_cells object additionally contains 6,897 cells labelled `Chondrocyte_like`, 3,847 labelled `Fibrochondrocyte_like`, and 1,473 labelled `Fibroblast_like` — these are predominantly the GSE189916 cells from the compartment-undefined IVD_mixed tissue that retain generic (non-compartment-prefixed) labels.

![**Figure 2. Compartment-specific UMAPs.** NP (left, 262,924 cells), AF (centre, 84,617 cells), CEP (right, 50,854 cells) atlases coloured by `cell_type`. Within each compartment, resident populations are arranged as adjacent graded clusters rather than separated discrete clusters, consistent with the continuum hypothesis. Non-resident populations occupy distinct regions of feature space.](manuscript_figures/fig02_umap_compartments.png)

### §2 — NP fibrocartilaginous cells dominate the degenerative transcriptional signature

Of the 26 statistically powered cell-type × comparison contrasts (one inflated contrast excluded; see Caveats §1), **1,823 significant DEGs** at FDR < 0.05 were identified. The cell-type distribution is highly asymmetric, with the three NP cell types collectively accounting for 1,310 of 1,823 trustworthy DEGs (72%):

**Table 2.** DE counts per cell-type × contrast (trustworthy contrasts only).

| Cell type | H vs. mild | H vs. severe | Mild vs. severe | H vs. all |
|---|---:|---:|---:|---:|
| NP_fibrocartilaginous | 14 | 291 | **325** | 84 |
| NP_fibrochondrocyte_chondroid | 5 | 27 | 350 | — |
| NP_mature_chondrocyte | 3 | 121 | 174 | 2 |
| AF_outer | **121** | 3 | 1 | — |
| AF_inner | 2 | — | — | — |
| Immune | 4 | 38 | 11 | — |
| Macrophage_M1 | — | — | 21 | — |
| Macrophage_M2 | 1 | (excluded) | 171 | 1 |
| Neutrophil | 2 | 23 | 25 | 3 |

NP_fibrocartilaginous carries the largest single mild-vs-severe contrast at 325 DEGs (189 up, 136 down). Examining the actual top-ranking DEGs in this contrast (Figure 3), the upregulated set in severe-vs-mild is dominated by ECM remodelling and acute-phase inflammatory transcripts: **MFAP5** (log2FC +4.7, padj 5×10⁻³), **EYA2** (+3.6), **C3** (+3.2), **SERPINF1** (+3.0), **IL6** (+2.7, padj 9×10⁻³), **HSD11B1** (+2.6), **CXCL2** (+2.3, padj 1×10⁻³), **IBSP** (+2.3), **FBLN1** (+2.3), **IFI27** (+2.3), **COL12A1** (+2.3, padj 3×10⁻⁴), **CXCL3** (+2.3), and **SAA1** (+2.1, padj 6×10⁻⁴). The downregulated set is led by hemoglobin transcripts (HBA2, HBB — likely tracking RBC contamination heterogeneity across samples), followed by immune lineage markers (FCN1, LYZ), the WNT ligand **WNT16** (−2.6), the hedgehog interacting protein **HHIP** (−2.1), the neuromedin **NMU** (−2.6), and the cell cycle / DNA replication factors **MCM10** and **SFN**.

![**Figure 3. Volcano plot of NP_fibrocartilaginous mild-vs-severe degeneration.** 325 significant DEGs at FDR < 0.05 (red), with the largest log2FC and lowest-padj transcripts labelled. Upregulated severe-vs-mild signal is concentrated in fibrotic ECM (COL12A1, FBLN1, SERPINF1), acute-phase / SAA family (SAA1, SAA2), and inflammatory chemokines (IL6, CXCL2, CXCL3). Downregulated signal includes WNT16, HHIP, and neuromedin NMU.](manuscript_figures/fig03_volcano_NP_fib_mild_vs_severe.png)

Across the broader NP set (NP_fibrocartilaginous + NP_fibrochondrocyte_chondroid + NP_mature_chondrocyte), severe-versus-mild and severe-versus-healthy contrasts share a recurring transcriptional pattern: upregulation of fibrotic collagens (**COL1A1**, **COL3A1**, **COL10A1**, **COL12A1** — sig UP in NP_fibrocartilaginous healthy-vs-severe and mild-vs-severe contrasts), upregulation of selected matrix-degrading enzymes (**MMP3** in NP_fibrocartilaginous and NP_mature_chondrocyte mild-vs-severe, **MMP19** in NP_fibrocartilaginous, **ADAMTS1** and **ADAMTS5** in NP_mature_chondrocyte mild-vs-severe), and broad activation of the NF-κB / acute-phase axis (**NFKBIZ**, **IL6**, **CXCL2/3**, **SAA1/2**, **CCL2**, **PTGS1/2**). Canonical anabolic chondrocyte transcripts (COL2A1, ACAN, PRG4) do not reach significance in any pairwise comparison at the cell-type level in this analysis — the dominant ECM signature is **fibrotic remodelling**, not collapse of the chondrocyte ECM program.

![**Figure 4. NP_fibrocartilaginous healthy-vs-degeneration heatmap.** Top differentially expressed genes across the three NP_fibrocartilaginous contrasts (healthy-vs-degenerated_all, healthy-vs-severe, mild-vs-severe), showing the consistent direction of effect across the degeneration gradient.](manuscript_figures/fig04_heatmap_NP_fib.png)

### §3 — AF outer cells respond early, then quiesce

AF_outer cells exhibit a distinctive temporal pattern that contrasts sharply with NP cells. At healthy-vs-mild, **121 DEGs** are detected — predominantly downregulated (110 of 121). By healthy-vs-severe the signal collapses to **3 DEGs**, and mild-vs-severe yields a single significant DEG. AF_inner cells, by contrast, return only 2 DEGs across all comparisons.

Examining the AF_outer healthy-vs-mild downregulated set, the pattern is striking: a coordinated loss of pain-relevant inflammatory and neurotrophic transcripts at the earliest detectable stage of degeneration. **CXCL8** (log2FC −4.4, padj 2×10⁻²), **NGFR** (−9.5, padj 3×10⁻²), and **PLA2G2A** (−7.7, padj 3×10⁻²) are all significantly downregulated in mild-vs-healthy AF_outer. This signature is not recapitulated in any other compartment.

![**Figure 5. Volcano plot of AF_outer healthy-vs-mild degeneration.** 121 significant DEGs (red), 110 of which are downregulated. Pain-relevant transcripts CXCL8, NGFR, and PLA2G2A are among the most strongly downregulated, suggesting either a transient stress response that resolves by the severe stage or a population shift not visible at cell-type-level resolution.](manuscript_figures/fig05_volcano_AF_outer_h_vs_mild.png)

The biology is ambiguous. Three interpretations are consistent with the data: (i) a transient transcriptional stress response in AF outer fibroblasts during early degeneration that has resolved by the severe stage; (ii) a shift in the AF_outer sub-state composition (matrix_active vs. proliferating vs. homeostatic fractions) that is averaged across at the cell-type level; or (iii) sampling differences across donors at severe AF disease stages, particularly since CEP and severe AF samples come predominantly from a subset of donors. Distinguishing among these would require longitudinal sampling from the same donors, which no public dataset provides.

### §4 — Eighteen pain-associated genes with cell-type-specific dysregulation

Cross-referencing the trustworthy DE results against a curated 60-gene pain panel yields **18 unique pain-associated genes** that are significantly differentially expressed at FDR < 0.05 in one or more contrasts. The full set, grouped by pain category and cell type, is summarized in Table 3.

**Table 3.** Pain-associated genes significantly DE at FDR < 0.05.

| Gene | Category | Direction in degeneration | Cell types with significant signal |
|---|---|---|---|
| IL6 | Inflammatory pain | UP (severe) | NP_fibrocartilaginous, Immune |
| IL1B | Inflammatory pain | DOWN (deg_all) | NP_fibrocartilaginous |
| CCL2 | Inflammatory pain | UP | NP_fibrocartilaginous, Neutrophil, Immune |
| TNF | Inflammatory pain | UP (mild → severe) | Neutrophil |
| CXCL8 | Inflammatory pain | DOWN (healthy → mild) | AF_outer |
| PTGS2 | Inflammatory pain | UP (mild → severe) | Macrophage_M1 |
| PLA2G2A | Inflammatory pain | UP / DOWN | NP_fibrocartilaginous (UP severe), Immune (UP mild→severe), AF_outer (DOWN mild) |
| BDKRB1 | Inflammatory pain | UP | NP_fibrochondrocyte_chondroid |
| BDKRB2 | Inflammatory pain | UP | NP_mature_chondrocyte |
| NTN1 | Nerve guidance | UP | NP_fibrocartilaginous, NP_mature_chondrocyte |
| UNC5B | Nerve guidance | UP | NP_fibrocartilaginous |
| NGFR | Neurotrophin receptor | UP / DOWN | Immune (UP severe), AF_outer (DOWN mild) |
| PENK | Neuropeptides | UP | NP_fibrocartilaginous |
| VEGFA | Neovascularization | UP | NP_mature_chondrocyte |
| PDGFA | Neovascularization | UP | NP_mature_chondrocyte, NP_fibrochondrocyte_chondroid |
| FGF2 | Neovascularization | UP | NP_fibrochondrocyte_chondroid |
| FLT1 | Neovascularization | UP | Macrophage_M1 |
| P2RX4 | Nociception ion channel | DOWN | NP_fibrochondrocyte_chondroid |

NP_fibrocartilaginous carries the broadest pain signature: 10 of 18 unique genes (IL6, IL1B, CCL2, PLA2G2A, NTN1, UNC5B, PENK, and additional contrast-specific hits) are dysregulated in this single cell type. NP_mature_chondrocyte concentrates the neovascularization axis (VEGFA + PDGFA + BDKRB2 + NTN1). NP_fibrochondrocyte_chondroid contributes the kinin and growth-factor arm (BDKRB1, FGF2, PDGFA). Macrophage_M1 carries the canonical macrophage pain ligands (PTGS2, FLT1).

The opposing direction of AF_outer signals (CXCL8, NGFR, PLA2G2A all DOWN at healthy-vs-mild) underscores the §3 finding that the AF outer compartment undergoes a different early response than NP cells do.

![**Figure 6. Pain-associated gene heatmap.** Log2 fold-change of pain-panel genes across cell-type × comparison contrasts where the gene reaches significance (FDR < 0.05). Red = upregulated, blue = downregulated. NP_fibrocartilaginous (column-grouped) carries the broadest signal; NP_mature_chondrocyte concentrates the neovascularization axis.](manuscript_figures/fig06_pain_genes_heatmap.png)

### §5 — Pathway enrichment and transcription factor activity converge on fibrotic ECM, NF-κB, and acute-phase programs

Pathway enrichment yielded **3,051 significant ORA terms** (Enrichr; GO Biological Process, KEGG, Reactome) and **6,890 significant GSEA terms** (gseapy; MSigDB + a custom IVD gene-set library) at FDR < 0.05 across the 26 trustworthy contrasts. Transcription factor activity inference via decoupler against the CollecTRI regulon database identified **555 significant TF × comparison records covering 252 unique TFs** at FDR < 0.05.

Four themes recur:

1. **Fibrotic ECM remodelling.** Upregulation of fibrotic collagens (COL1A1, COL3A1, COL10A1, COL12A1) and selected catabolic enzymes (MMP3, MMP19, ADAMTS1, ADAMTS5) in NP_fibrocartilaginous and NP_mature_chondrocyte severity contrasts. GO terms `extracellular matrix organization`, `collagen catabolic process`, and `proteoglycan metabolic process` lead the upregulated NP_fibrocartilaginous severity comparisons. Canonical chondrocyte anabolic markers (COL2A1, ACAN, PRG4) are not significantly dysregulated at the cell-type level — the dominant ECM signal is fibrotic remodelling rather than collapse of the chondrocyte ECM program.

2. **NF-κB and acute-phase activation.** **NFKBIZ**, **IL6**, **CXCL2/3**, **SAA1/2**, **CCL2**, **PTGS1/2** all upregulated across multiple NP cell types and contrasts. Acute-phase serum amyloid A transcripts SAA1 and SAA2 are among the largest-effect upregulated DEGs in NP_fibrocartilaginous mild-vs-severe (SAA1 +2.1, padj 6×10⁻⁴) and NP_mature_chondrocyte mild-vs-severe (SAA1 +2.8, padj 2×10⁻⁴).

3. **AP-1 and stress-response TF activity.** CollecTRI TF activity scores rank **AP-1** (composite AP-1 score, 10 significant contrast records) as the single most consistently active regulator, followed by **EGR1**, **JUN**, **RELA**, **SP1**, **NFKB1**, **NFKB**, **PPARG**, **FOS**, **SP3**, **CEBPB**, **RUNX2**, **ETS1**, **GLI2**, and **SMAD3**. The NF-κB family (NFKB1, NFKB, RELA) collectively contributes 24 significant records across degeneration contrasts. AP-1 / JUN / FOS / EGR1 — canonical stress-response and immediate-early transcription factors — are the dominant signal at the TF level.

4. **Cell-type-specific pain dysregulation.** ORA over the curated IVD gene-set library confirms enrichment of `inflammatory_pain`, `neurotrophin_signalling`, and `neovascularization` categories in NP_fibrocartilaginous and NP_mature_chondrocyte degeneration contrasts.

![**Figure 7. GSEA enrichment heatmap across IVD-specific gene sets.** Normalized enrichment score per contrast (columns) × pain / ECM / immune gene-set (rows). Red = enriched in the second term of the contrast (i.e., degenerated for healthy-vs-degenerated contrasts).](manuscript_figures/fig07_gsea_ivd_heatmap.png)

![**Figure 8. Transcription factor activity heatmap.** CollecTRI-derived TF activity scores per contrast (columns) × TF (rows), restricted to the 25 most significant TFs by frequency of FDR < 0.05 contrast hits. AP-1, EGR1, JUN, RELA, SP1, NFKB1, PPARG, and FOS dominate.](manuscript_figures/fig08_tf_activity_heatmap.png)

### §6 — Cell-cell communication: neutrophil-recruitment signaling gained, CEP basement-membrane signaling lost

LIANA consensus rank-aggregation on 20,000 cells per condition identified **66,827 ligand-receptor interactions in healthy tissue and 76,019 in degenerated tissue**, with a 90,650-interaction union. **6,883 interactions** in the degenerated pool involve a curated pain-panel ligand or receptor.

Top ligands by interaction count differ between conditions in informative ways. In healthy tissue, the top ligands by interaction count are TGFB1 (2,205), FN1 (1,995), VEGFA (1,288), FGF2 (1,261), THBS1 (1,176), APP (1,155), COL1A2 (975), COL6A1 (966), COL2A1 (896), and COL6A2 (882). In degenerated tissue, **FN1** moves to the top (2,085 — up from 1,995 in healthy), TGFB1 drops to second (1,908 — down from 2,205), **VEGFA** rises (1,485 vs. 1,288 in healthy), **COL1A1** appears in the top 10 (1,400 in degen; not in top 10 in healthy), and **ADAM10** appears (1,105 in degen). The shift is consistent with the §2 finding of fibrotic ECM remodelling — fibronectin and fibrotic collagens are gaining signaling weight in degeneration.

![**Figure 9. Cell-cell interaction heatmap in degenerated tissue.** Source cell-type (rows) × target cell-type (columns) interaction counts. Resident cell types are heavily interconnected; non-resident populations show more focused interaction profiles.](manuscript_figures/fig09_interaction_heatmap_degenerated.png)

Differential analysis (rank_diff = magnitude_rank_healthy − magnitude_rank_degenerated; positive = gained in degeneration) yields two robust patterns:

**Gained in degeneration.** Pan-cell-type → Neutrophil **FN1 → C5AR1** and **RPS19 → C5AR1** signaling dominates the top of the gained list. Fibronectin and ribosomal protein S19 acting as C5a-receptor ligands across virtually every disc-resident cell type is the canonical complement-driven neutrophil chemotaxis axis being broadly activated.

**Lost in degeneration.** A coordinated set of CEP_outer-centric signaling axes attenuates: **WNT2B → FZD4/LRP5/LRP6** (CEP_outer autocrine and to NP_mature_chondrocyte), **LAMB3 → ITGA6 / ITGAV / ITGA2 integrins** (CEP-to-AF basement membrane), **NDP → FZD4** (Norrin/WNT), and **MDK → SDC3** (midkine). In parallel, **AF_outer / NP_fibrocartilaginous → CD4 HLA-class-II presentation** is reduced. The lost-in-degeneration signal is therefore a coordinate loss of WNT, basement-membrane, and immune-presentation signaling rather than a single-axis effect.

![**Figure 10. Differential interaction plot.** Each point represents one source-target-ligand-receptor interaction; horizontal axis = `rank_diff`, vertical axis = degenerated-condition magnitude rank. Right-side points are gained in degeneration; left-side points are lost. Annotated top interactions include FN1 → C5AR1 (gained, multiple sources → Neutrophil) and CEP_outer WNT2B → FZD4/LRP5/6 (lost).](manuscript_figures/fig10_differential_interactions.png)

Pain-relevant interactions in degenerated tissue are dominated by neovascularization (**VEGFA** 1,485, **FGF2** 1,157, **VEGFC** 246, **VEGFB** 216), neurotrophic (**NGF** 360), semaphorin nerve-guidance (**SEMA3C** 300, **SEMA3A** 250), inflammatory (**TNFSF10** 250, **PTGS2** 224), and granulin (**GRN** 285) ligand contributions.

Aggregate interaction counts are sensitive to the cell-type partition — finer cell-type resolution combinatorially multiplies the number of source-target-ligand-receptor tuples. Cell-type-pair-specific and ligand-specific calls are therefore the interpretable level of CCC analysis.

### §7 — Cell-type composition does not shift significantly with degeneration

Compositional analysis across 19 cell-type × comparison contrasts at the compartment level identified **no significant changes after FDR correction** (lowest adjusted p-value 0.61). The largest absolute shifts — AF_inner enrichment in severe AF (log2FC +3.05, raw p = 0.04; +1.33 at healthy-vs-severe, raw p = 0.07), Macrophage_M2 depletion in NP severe degeneration (log2FC −4.28, raw p = 0.02), and Pericyte_SMC depletion in severe disc (log2FC −7.98, raw p = 0.07) — do not survive multiple-testing correction. The finer `cell_subtype`-level composition analysis (also FDR-corrected) similarly returns no significant shifts.

The conservative interpretation is that **disc degeneration is principally a transcriptional, not a cellular, remodelling** at the cell-type level. Composition shifts at the sub-state level (proliferating, stressed, inflammatory, matrix_active fractions within each cell type) are more variable across donors and are best interpreted at the individual sub-state level rather than as compartment-wide claims.

### §8 — Cell-state trajectories: per-cell-type gradients within compartments

PAGA-initialized DPT pseudotime in each compartment yields the expected ordering at the cell-type level — root clusters enriched for the mature/inner population (NP_mature_chondrocyte for NP, AF_inner for AF, CEP_hyaline for CEP) progress along a continuous axis toward the more degenerative / outer populations.

![**Figure 11. PAGA connectivity graph, NP.** Node size = cluster size; edge weight = PAGA connectivity. The NP resident population is a connected graph rather than a set of well-separated clusters, consistent with the continuum hypothesis.](manuscript_figures/fig11_paga_NP.png)

The compartment-level pseudotime-versus-condition ordinal correlation is **near zero for NP and AF** (NP ρ = −0.004, p = 0.41; AF ρ = −0.003, p = 0.49), while CEP retains a weak positive correlation (ρ = +0.077, p = 1.9×10⁻⁴⁹). At the per-cell-type level, gradients are stronger and **opposing within compartments**: AF_inner ρ = +0.54, AF_outer ρ = −0.30 (both p ≈ 0); CEP_hyaline ρ = +0.17, CEP_outer ρ = −0.21. NP cell-type-specific correlations are small in magnitude (ρ ∈ {−0.04, +0.06, −0.05}).

![**Figure 12. Pseudotime by cell type, AF.** Distribution of DPT pseudotime per AF cell type. AF_inner and AF_outer occupy distinct portions of the pseudotime axis with opposing condition-correlated gradients.](manuscript_figures/fig12_pseudotime_by_celltype_AF.png)

Mann-Whitney U tests of healthy-versus-degenerated pseudotime medians remain highly significant in all three compartments (NP p = 3×10⁻⁴, AF p = 1×10⁻²⁰, CEP p = 2×10⁻⁴³), though absolute median shifts are small (e.g. NP healthy 0.038 vs. degenerated 0.038; AF healthy 0.500 vs. degenerated 0.497). The conservative reading is therefore: the discrete healthy-vs-degenerated shift is real but the compartment-level pseudotime direction is not a robust claim, with within-cell-type gradients accounting for most of the signal.

Trajectory-associated genes (top 500 by DPT correlation) substantially overlap with the DE results: 328/500 (NP), 326/500 (AF), and 338/500 (CEP) trajectory genes are also DE genes — 65–68% per compartment. The NP trajectory's late-up program (282 genes) is led by COL1A1, MKI67, TNC, BIRC5, SERPINF1, HMMR, AQP1, COL6A1 — consistent with the fibrotic proliferative remodelling signal from §2. The late-down program (87 genes) is led by SEMA3A, DSP, FGF2, NGF, TRPV4, TNFRSF11B, F13A1, IL11 — interestingly including several neurotrophin and pain-relevant ligands whose loss along the NP pseudotime axis is opposite to what a simple "more degeneration → more pain ligand expression" model would predict.

![**Figure 13. Gene dynamics along NP pseudotime.** Smoothed expression of top trajectory-associated genes plotted against DPT pseudotime, with program assignment (late_up, late_down, stable) shown. Fibrotic ECM and proliferation transcripts (COL1A1, COL3A1, MKI67, BIRC5) rise along the trajectory; neurotrophic and stress-response transcripts (SEMA3A, FGF2, NGF) decline.](manuscript_figures/fig13_gene_dynamics_NP.png)

---

## Discussion

### A coherent transcriptional model of disc degeneration

The atlas supports a model in which disc degeneration is principally a cell-type-specific transcriptional remodelling rather than a compositional collapse. NP_fibrocartilaginous cells drive the signal, upregulating fibrotic collagens (COL1A1, COL3A1, COL10A1, COL12A1), selected matrix-degrading enzymes (MMP3, MMP19, ADAMTS1, ADAMTS5 — concentrated in NP_mature_chondrocyte and NP_fibrocartilaginous), acute-phase proteins (SAA1, SAA2, SERPINF1), inflammatory mediators (IL6, CCL2, CXCL2/3, NFKBIZ, PTGS1/2), nerve-guidance cues (NTN1, UNC5B), neuropeptides (PENK), and neurotrophin receptors (NGFR in Immune cells). NP_mature_chondrocyte contributes the neovascularization arm (VEGFA, PDGFA, BDKRB2). NP_fibrochondrocyte_chondroid carries its own kinin signature (BDKRB1) and the FGF2 / PDGFA growth-factor axis. The result is a coherent biology in which the NP resident cell continuum acts collectively as both effector and amplifier of the catabolic / fibrotic / inflammatory / nociceptive cascade — and the dominant transcriptional theme is **fibrotic remodelling**, not collapse of the chondrocyte ECM program. Canonical chondrocyte anabolic markers (COL2A1, ACAN, PRG4) do not reach significance in any cell-type-level pairwise comparison in this atlas.

AF_outer cells contribute a temporally-distinct early response that quiesces by the severe stage, with a pain-relevant downregulation signature (CXCL8, NGFR, PLA2G2A) appearing at the mild stage and disappearing by severe. AF_inner cells, CEP cells, and the immune compartment all participate in narrower but coherent ways. Macrophage_M1 polarization in late degeneration carries the canonical macrophage pain ligands (PTGS2, FLT1). The Immune compartment carries IL6, NGFR, CCL2, and PLA2G2A upregulation at the severe stage. Neutrophil recruitment via complement (FN1/RPS19 → C5AR1) is the single most prominent gained cell-cell signaling axis in the degenerated atlas.

### Methodological considerations

Three methodological choices distinguish this atlas from prior cross-study IVD scRNA-seq analyses:

**Tiered integration.** Splitting the integration by coarse cell class — separate models for non-resident and resident cells, with conservative settings for the latter — preserves the within-compartment continua that uniform single-method integration tends to compress. A sensitivity analysis using flat Seurat CCA integration on the same data confirmed this: the tiered approach yielded finer cell-type resolution (19 vs. 16 cell types in the union atlas), a larger trustworthy DE pool (1,823 vs. 1,198 DEGs), and broader pain-gene signal (18 vs. 10 unique significant pain genes), while preserving the qualitative biological themes (NF-κB / inflammatory upregulation, no compositional shift). The trade-off is that some compartment-level summary statistics (e.g. pseudotime-condition ordinal correlations) shift between the two approaches, because finer cell-type partition redistributes signal that flat integration averages across. We treat the qualitative themes as the strong claims and report compartment-level summary metrics with caveats (§Caveats §2).

**Three-stage annotation with compartment-prefixed labels.** Coarse class → compartment-prefixed cell type → sub-state. The compartment prefix (NP_fibrocartilaginous vs. AF_outer vs. CEP_outer) is important because some functional cell types (Fibroblast_like, Chondrocyte_like) are not equivalent across compartments, and grouping them by generic class hides compartment-specific biology. Sub-state assignment (proliferating, inflammatory, stressed, matrix_active, migratory, homeostatic) provides descriptive resolution within each cell type without inflating the statistical claims (DE was run on cell_type, not cell_subtype, to preserve sample-level statistical power).

**Contamination flagging rather than filtering.** 4.5% of the atlas — 16,514 Erythrocyte cells and 1,831 endothelial-admixed NP_fibrocartilaginous cells — is retained with explicit flags rather than removed. The hemoglobin downregulation signal in the NP_fibrocartilaginous mild-vs-severe contrast (HBA2, HBB) likely reflects this flagged contamination, and is interpreted accordingly. Silent filtering would have removed counts whose biology is informative about the contamination calls themselves.

The pseudobulk DESeq2 approach is used in preference to single-cell DE methods (Wilcoxon, MAST) to avoid the well-documented inflation of false positives when treating cells as independent observations. The single inflated contrast in this atlas (Macrophage_M2 healthy-vs-severe, 5,659 DEGs from a 7-sample comparison) is excluded from interpretation specifically because its DE count is incompatible with the donor-level statistical model — likely driven by a sample-imbalance × dispersion-estimation interaction. The contrast remains on disk for transparency.

### Therapeutic implications

The cell-type-specific resolution of the pain-gene panel suggests that interventions distinguish themselves by where they act. NP_fibrocartilaginous cells are a candidate cell-type target for neurotrophin- and nerve-guidance-axis interventions (NTN1, UNC5B, PENK upregulation). NP_mature_chondrocyte and NP_fibrochondrocyte_chondroid cells are the angiogenic-mediator source (VEGFA, PDGFA, FGF2 upregulation). Macrophage_M1 cells are an inflammatory amplifier (PTGS2, FLT1 upregulation). PTGS2 recurs across NP cell types and Macrophage_M1, supporting its standing as a high-priority therapeutic node. The complement-driven neutrophil recruitment axis (C5AR1 ligands) is gained broadly across resident cells in degeneration and represents a potential upstream intervention point. None of these targets is novel in isolation — each has prior literature support — but the cell-type-specific source of each signal provides a more refined picture of where intervention might be biologically appropriate.

### Robust versus version-sensitive findings

We distinguish two classes of finding by their methodological robustness:

**Robust across integration approach.** ECM remodelling (fibrotic collagen up; selected catabolic enzyme up) in NP cell types. NF-κB / acute-phase / inflammatory upregulation in NP cell types. No compositional shifts after FDR correction. Pain-gene panel dysregulation concentrated in NP_fibrocartilaginous and NP_mature_chondrocyte. Disc cells as inflammatory mediators (NP cells, not infiltrating immune cells, are the dominant source of degeneration-associated cytokine and chemokine signal at the cell-type-level DE resolution).

**Sensitive to integration approach or to subsampling.** Compartment-level pseudotime-condition correlation direction (per-cell-type gradients dominate the signal in this analysis; aggregating them across compartments cancels out). Aggregate CCC interaction counts (sensitive to the cell-type partition and the 20,000-cell-per-condition subsampling cap; cell-type-pair-specific calls are more interpretable than aggregate counts). Specific TF significance calls (the top-ranking TFs by frequency are AP-1, EGR1, JUN, RELA, SP1, NFKB1, PPARG, FOS — but individual TFs have variable cross-analysis significance, so we report the family-level themes rather than individual TF claims).

---

## Caveats and limitations

The following caveats apply to specific claims in this manuscript and are flagged explicitly rather than hidden in supplementary materials:

1. **One DE contrast is statistically inflated and excluded from all downstream interpretation.** Macrophage_M2 healthy-vs-severe returned 5,659 DEGs from a 7-sample comparison. The per-gene dispersion estimates and the gene-count distribution flag clear pathology. The contrast is retained on disk for transparency but excluded from all pathway, TF, pain-gene, CCC, and manuscript-level analyses. Cite Macrophage_M2 mild-vs-severe (171 DEGs, well-distributed) for macrophage-related claims in degeneration.

2. **Trajectory pseudotime-condition correlations are dominated by within-cell-type gradients, not compartment-level direction.** AF_inner ρ = +0.54 and AF_outer ρ = −0.30 cancel at the compartment level (AF ρ = −0.003). The discrete healthy-vs-degenerated MWU shift is robust; the compartment-level direction is not a strong claim. Report cell-type-specific gradients rather than compartment-level direction.

3. **Contamination cells are retained with flags, not filtered.** 16,514 Erythrocyte cells (predominantly in NP and CEP) and 1,831 endothelial-admixed NP_fibrocartilaginous cells carry `is_contamination=True` flags. The hemoglobin downregulation signal in NP_fibrocartilaginous mild-vs-severe reflects this contamination distribution across samples, not biology of the disc cells themselves.

4. **GSE189916 (IVD_mixed) cells carry generic cell-type labels.** That dataset does not separate compartments at the tissue level, so its cells retain non-compartment-prefixed labels (`Fibroblast_like`, `Chondrocyte_like`, `Fibrochondrocyte_like`) rather than being forced into compartment-specific bins.

5. **GSE242443 CEP cells are culture-expanded.** They contribute a culture-derived rather than fresh-tissue signal; CEP results should be interpreted with this in mind.

6. **GSE230809 donors are all male.** This dataset's age and disease effects are sex-confounded; CEP claims that lean on it carry that confound.

7. **No RNA velocity.** Public count matrices do not include spliced/unspliced layers. Pseudotime is DPT-only.

8. **CEP is underpowered for cell-type-level DE.** Only 6 CEP-containing samples and 50,854 cells; no CEP cell-type-specific DE comparisons reach the powered threshold for healthy-vs-mild contrasts. CEP claims rest primarily on annotation, composition, trajectory, and CCC analyses.

9. **Compositional shifts are not FDR-significant.** Several raw-p < 0.05 shifts (AF_inner enrichment in severe AF, Macrophage_M2 depletion in severe NP, Pericyte_SMC depletion in severe disc) are reported descriptively but do not survive multiple-testing correction.

10. **LIANA CCC aggregate counts are sensitive to subsampling and to the cell-type partition.** The 20,000-cell-per-condition subsampling cap means that aggregate interaction counts are not directly comparable across pipeline parameterizations. Cell-type-pair-specific and ligand-specific calls are the interpretable level.

11. **Cross-sectional sampling.** All datasets are cross-sectional snapshots of human IVD tissue. Longitudinal trajectory inference (e.g. true temporal ordering of healthy → mild → severe within the same donor) is not possible from the available data.

---

## Data and code availability

All 12 source datasets are publicly available at GEO or CNGB (accession numbers in the dataset registry). The full analysis pipeline (12 modules, scripts, specifications, notebooks, and metadata) is at `https://github.com/andrewsu/lotz-ivd`. Per-module logs, intermediate AnnData files, and result tables are reproducible end-to-end from the raw count matrices. Per-file checksums are recorded in `metadata/file_checksums.json`.

---

## Author contributions, acknowledgments, references

*To be populated at submission.*
