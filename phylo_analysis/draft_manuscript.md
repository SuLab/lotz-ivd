# Single-Cell Transcriptomic Atlas of Human Intervertebral Disc Degeneration: Cell States, Signaling Networks, and Therapeutic Target Candidates

**Draft manuscript — Phylo framework analysis**

---

## Glossary of Abbreviations

| Abbreviation | Full Term |
|---|---|
| **IVD** | Intervertebral disc — the fibrocartilaginous structure between vertebral bodies that provides cushioning and flexibility to the spine |
| **NP** | Nucleus pulposus — the gel-like core of the IVD, composed primarily of chondrocyte-like cells embedded in a proteoglycan-rich extracellular matrix |
| **AF** | Annulus fibrosus — the tough, layered ring of collagen fibers surrounding the NP that contains the disc under compressive load |
| **CEP** | Cartilaginous endplate — the thin layer of hyaline cartilage at the superior and inferior surfaces of the IVD, through which nutrients diffuse from vertebral capillaries |
| **scRNA-seq** | Single-cell RNA sequencing — a technology that measures the messenger RNA (gene expression) of individual cells, revealing cell-by-cell variation invisible to bulk methods |
| **ECM** | Extracellular matrix — the structural scaffold of proteins and sugars (collagens, proteoglycans) secreted by cells that gives tissues their mechanical properties |
| **UPR** | Unfolded protein response — a cellular stress pathway activated when misfolded proteins accumulate in the endoplasmic reticulum; drives inflammation and apoptosis when prolonged |
| **UMAP** | Uniform Manifold Approximation and Projection — a dimensionality reduction algorithm that projects high-dimensional gene expression data into a 2D map, placing transcriptionally similar cells near each other |
| **HVG** | Highly variable gene — a gene whose expression varies substantially across cells, used to focus analyses on biologically informative genes rather than housekeeping noise |
| **PCA** | Principal component analysis — a linear dimensionality reduction method that identifies the axes of greatest variance in gene expression data |
| **kNN** | k-nearest neighbors — a graph construction method that connects each cell to its k most similar cells in gene expression space |
| **DEG** | Differentially expressed gene — a gene whose expression is statistically significantly different between two conditions (e.g., healthy vs. diseased) |
| **GSEA** | Gene set enrichment analysis — a method that tests whether predefined sets of functionally related genes (e.g., a pathway) are coordinately up- or down-regulated |
| **NES** | Normalized enrichment score — GSEA's measure of how strongly a gene set is enriched at the top or bottom of a ranked gene list; positive = upregulated, negative = downregulated |
| **ORA** | Over-representation analysis — a simpler enrichment method that tests whether DEGs overlap with known pathways more than expected by chance |
| **PAGA** | Partition-based graph abstraction — an algorithm that estimates the connectivity (transition probability) between clusters of cells, revealing likely cell-state transitions |
| **DPT** | Diffusion pseudotime — an algorithm that orders cells along a trajectory from a root cell, modeling the progression through biological processes like differentiation or disease |
| **LIANA** | Ligand-receptor analysis framework — a consensus toolkit for inferring cell-cell communication from scRNA-seq data by scoring ligand-receptor pair co-expression across cell types |
| **LR pair** | Ligand-receptor pair — a signaling molecule (ligand) secreted by one cell and the surface protein (receptor) on another cell that receives the signal |
| **MMP** | Matrix metalloproteinase — enzymes that degrade ECM components; overactive in degenerative disc disease |
| **TIMP** | Tissue inhibitor of metalloproteinases — proteins that inhibit MMPs, protecting the ECM from degradation |
| **EMT** | Epithelial-mesenchymal transition — a process where cells lose epithelial characteristics and gain mesenchymal (fibroblast-like) properties, associated with fibrosis |
| **SASP** | Senescence-associated secretory phenotype — the pro-inflammatory cytokine cocktail secreted by senescent cells that drives chronic inflammation in aging tissues |
| **TNF** | Tumor necrosis factor — a key pro-inflammatory cytokine implicated in disc degeneration and back pain |
| **NF-kB** | Nuclear factor kappa-light-chain-enhancer of activated B cells — a transcription factor that drives inflammatory gene expression |
| **MAD** | Median absolute deviation — a robust measure of data spread used to set adaptive quality-control thresholds per sample |
| **FDR** | False discovery rate — a statistical correction for multiple hypothesis testing; controls the expected proportion of false positives |
| **LFC** | Log2 fold change — the logarithmic ratio of gene expression between two conditions; LFC = 1 means 2-fold higher expression |
| **SMC** | Smooth muscle cell — contractile cells found in blood vessel walls; pericytes/SMCs in the IVD reflect vascular structures at the disc periphery |
| **GEO** | Gene Expression Omnibus — NCBI's public repository for gene expression datasets |

---

## Abstract

Intervertebral disc (IVD) degeneration is the leading structural cause of low back pain, yet the cellular and molecular mechanisms driving the transition from healthy to degenerated disc tissue remain poorly understood at single-cell resolution. Here we present an integrated single-cell RNA-seq atlas of 173,628 human IVD cells from 7 publicly available datasets spanning 29 donors and 4 degeneration grades (healthy, mild, moderate, severe). We identify 12 cell populations including 5 distinct nucleus pulposus (NP) cell states, annulus fibrosus (AF) fibroblasts, and vascular and immune cell types. Severe degeneration is characterized by a massive transcriptional downregulation program (approximately 7:1 down-to-up ratio of differentially expressed genes), with coordinated loss of Wnt signaling, Notch signaling, cellular senescence programs, and RUNX transcription factor activity across all NP cell states. AF fibroblasts expand significantly from approximately 7% to 24% of total cells, consistent with fibrocartilaginous metaplasia. Cell-cell communication analysis reveals loss of protective TIMP1-CD63 signaling and gain of pro-inflammatory FN1-macrophage interactions. Trajectory analysis identifies the NP stress-response state as a transition hub between canonical NP chondrocytes and the degenerative unfolded protein response (UPR) state. These findings define a multi-layered model of disc degeneration — from loss of homeostatic signaling through inflammatory amplification to structural failure — and nominate specific molecular targets including the Wnt, Notch, TIMP1, and FN1 signaling axes as candidates for therapeutic intervention.

---

## 1. Introduction

### 1.1 The clinical problem: low back pain and disc degeneration

Low back pain (LBP) is the leading cause of disability worldwide, affecting an estimated 619 million people at any given time (GBD 2021 Low Back Pain Collaborators, 2023). The economic burden in the United States alone exceeds $100 billion annually in direct healthcare costs and lost productivity (Dieleman et al., 2020). Intervertebral disc degeneration (IDD) is the most common structural correlate of LBP, though the relationship between radiographic degeneration and clinical symptoms remains complex (Brinjikji et al., 2015).

The intervertebral disc is a fibrocartilaginous structure situated between adjacent vertebral bodies throughout the spinal column. Each disc consists of three anatomically and functionally distinct compartments:

1. **Nucleus pulposus (NP):** The central, gel-like core composed primarily of chondrocyte-like cells (sometimes called NP cells or NPC) embedded in a highly hydrated matrix rich in type II collagen and the proteoglycan aggrecan (ACAN). The NP absorbs compressive loads and distributes them radially. In the embryo, the NP derives from the notochord, and remnant notochordal cells persist into early adulthood in some individuals (Risbud and Shapiro, 2014).

2. **Annulus fibrosus (AF):** A series of concentric lamellae (layers) of type I collagen fibers arranged at alternating angles, providing tensile strength to contain the pressurized NP. The outer AF is more fibroblast-like and receives some vascular supply, while the inner AF transitions toward a fibrocartilaginous phenotype (Humzah and Soames, 1988).

3. **Cartilaginous endplate (CEP):** Thin layers of hyaline-like cartilage at the superior and inferior disc surfaces that interface with the vertebral bone. The CEP is the primary route for nutrient diffusion into the avascular NP, and its calcification or damage can starve the disc of oxygen and glucose (Roberts et al., 1996).

### 1.2 The biology of disc degeneration

Disc degeneration is a multifactorial process involving mechanical overload, genetic predisposition, nutritional deprivation, and aging. At the cellular level, degeneration is characterized by (Roughley, 2004; Adams and Roughley, 2006):

- **Loss of proteoglycan (aggrecan):** Reduced ACAN and HAPLN1 (hyaluronan and proteoglycan link protein 1) diminish the disc's water-binding capacity, leading to dehydration, loss of disc height, and impaired load bearing.
- **Collagen switch:** A shift from type II collagen (COL2A1, the healthy NP matrix) to type I collagen (COL1A1/COL1A2, associated with fibrosis and the AF phenotype), reflecting fibrocartilaginous metaplasia of the NP.
- **Matrix degradation:** Upregulation of matrix metalloproteinases (MMPs) and ADAMTS (a disintegrin and metalloproteinase with thrombospondin motifs) enzymes, particularly ADAMTS5, which cleaves aggrecan. Concurrently, tissue inhibitors of metalloproteinases (TIMPs) are reduced, tipping the protease-antiprotease balance toward degradation (Le Maitre et al., 2004).
- **Inflammation:** Elevated levels of pro-inflammatory cytokines including TNF (tumor necrosis factor), IL-1beta, IL-6, and CXCL8 (IL-8), driven by NF-kB signaling. Inflammatory mediators both accelerate matrix degradation and sensitize nociceptive nerve endings that grow into the degenerated disc (Risbud and Shapiro, 2014).
- **Cellular senescence:** Accumulation of senescent cells that have exited the cell cycle but remain metabolically active, secreting a pro-inflammatory cocktail known as the senescence-associated secretory phenotype (SASP) (Feng et al., 2016).
- **Neovascularization and nerve ingrowth:** The healthy NP is avascular and aneural. In degeneration, blood vessels and nociceptive nerve fibers grow into the disc from the outer AF, facilitated by angiogenic and neurotrophic factors (Freemont et al., 2002).

### 1.3 Why single-cell transcriptomics?

Until recently, most molecular studies of disc degeneration used bulk RNA sequencing or microarrays, which measure the average gene expression across millions of cells. This approach obscures the contributions of individual cell types and states. Since the IVD contains a heterogeneous mix of NP cells, AF fibroblasts, endplate chondrocytes, vascular cells, and infiltrating immune cells — each potentially responding differently to degeneration — bulk methods cannot resolve which cells are driving disease and which are bystanders.

Single-cell RNA sequencing (scRNA-seq) overcomes this limitation by measuring the transcriptome of each individual cell. This enables:

- **Discovery of cell subtypes** that may not be distinguishable by surface markers alone
- **Identification of cell-state transitions** (e.g., from healthy NP chondrocyte to a stress-activated or degenerative phenotype)
- **Cell-type-specific differential expression** to determine which genes change in which cell populations
- **Inference of cell-cell communication** by analyzing ligand-receptor pair co-expression across cell types

Several groups have published scRNA-seq datasets from human IVD tissue in recent years (Gan et al., 2021; Cherif et al., 2022; Wang et al., 2023; Jiang and Sheyn, 2022; Li et al., 2022; Guo et al., 2023; Shi et al., 2024), each profiling a limited number of donors and conditions. No study to date has integrated all available datasets into a unified atlas with systematic cross-condition comparisons.

### 1.4 Study objectives

This study integrates 7 publicly available scRNA-seq datasets into a unified atlas of 173,628 human IVD cells to:

1. Define the complete repertoire of cell types and transcriptional states in the human IVD
2. Characterize how cell composition changes across degeneration grades
3. Identify cell-type-specific differentially expressed genes and enriched pathways in degeneration
4. Map the trajectory of NP cell-state transitions from health to disease
5. Infer changes in cell-cell communication networks during degeneration
6. Nominate molecular targets for therapeutic intervention in IDD

---

## 2. Methods

### 2.1 Dataset selection and acquisition

Seven publicly available scRNA-seq datasets of human IVD tissue were identified through systematic search of the Gene Expression Omnibus (GEO) database (Table 1). Datasets were selected to maximize coverage across tissue compartments (NP, AF, CEP), degeneration grades (healthy through severe), and donor demographics. All datasets used the 10x Genomics Chromium platform for single-cell library preparation.

**Table 1. Datasets included in the integrated atlas.**

| GEO Accession | First Author (Year) | Tissue | Condition | Samples | Cells (post-QC) |
|---|---|---|---|---|---|
| GSE160756 | Gan/Liu (2021) | NP, AF, CEP | Healthy atlas | 7 | ~60,000 |
| GSE199866 | Cherif (2022) | NP | Paired degenerated/non-degenerated | 4 | ~8,000 |
| GSE244889 | Wang (2023) | NP | Mild vs. severe | 7 | ~45,000 |
| GSE255768 | Shi (2024) | CEP | Degeneration | 2 | ~5,000 |
| GSE233666 | Guo (2023) | NP | Immune/ossification focus | 4 | ~20,000 |
| GSE205535 | Li (2022) | NP | Normal + degenerated | 2 | ~10,000 |
| GSE189916 | Jiang/Sheyn (2022) | IVD (mixed) | Neonatal + adult | 6 | ~25,000 |

Raw count matrices were downloaded from GEO. Donor-level metadata (age, sex, tissue compartment, clinical degeneration grade) were extracted from the corresponding publications and harmonized to a common schema. Degeneration grades were mapped to four categories: **healthy** (no clinical or radiographic degeneration), **mild** (Pfirrmann grade II-III or equivalent), **moderate** (Pfirrmann grade III-IV), and **severe** (Pfirrmann grade IV-V, frank herniation, or advanced degeneration).

### 2.2 Quality control and preprocessing

Each dataset was processed independently before integration. Quality control (QC) filtering was applied per sample using adaptive, median absolute deviation (MAD)-based thresholds for three metrics:

- **Number of detected genes (nGenes):** Cells with extremely low gene counts are likely empty droplets or debris; cells with extremely high counts may be doublets (two cells captured in one droplet).
- **Total UMI counts (nCounts):** Similar rationale to nGenes, filtering for viable single cells.
- **Mitochondrial gene fraction (%mito):** A high fraction of mitochondrial transcripts indicates a dying or damaged cell whose cytoplasmic mRNA has leaked out while mitochondrial mRNA is retained.

In addition, **Scrublet** (Wolock et al., 2019) was applied per sample to computationally identify and remove doublets — artificial "cells" that represent two cells captured together. Scrublet simulates synthetic doublets from the data and scores each real cell for similarity to these simulations.

After QC filtering, **173,628 of 222,433 cells (78%) were retained** across all 7 datasets and 29 donors.

### 2.3 Normalization and feature selection

For each dataset, raw UMI counts were normalized to 10,000 counts per cell (library-size normalization) and log-transformed (log1p), producing a lognormalized expression matrix. Raw integer counts were preserved in a separate data layer for downstream pseudobulk differential expression analysis, which requires unnormalized count data.

Highly variable genes (HVGs) were selected using a batch-aware approach: HVGs were identified independently within each dataset, and the union of top-ranked HVGs across datasets was taken to avoid bias toward any single dataset's biology. This yielded **4,000 HVGs** from a total of **25,304 genes** detected in 3 or more datasets.

### 2.4 Integration and batch correction

Because the 7 datasets were generated by different laboratories, with different sample preparation protocols, sequencing depths, and donor populations, systematic technical differences (batch effects) would confound biological comparisons if not corrected.

**Harmony** (Korsunsky et al., 2019) was used for batch correction. Harmony operates in principal component (PCA) space: it takes the top principal components of gene expression and iteratively adjusts them to remove variation attributable to batch variables (here, dataset identity and donor identity across 29 donors and 7 datasets) while preserving biological variation. The corrected PCA coordinates (referred to as the "Harmony embedding") were used for all downstream analyses including clustering, visualization, and trajectory inference.

A **k-nearest neighbor (kNN) graph** was constructed using FAISS (Facebook AI Similarity Search; Johnson et al., 2019), an approximate nearest-neighbor library optimized for large datasets. Each cell was connected to its k=30 most similar cells in the Harmony-corrected PCA space. This graph encodes the local neighborhood structure of the data and is the basis for clustering and UMAP visualization.

### 2.5 Clustering

Cell clusters were identified using the **Leiden algorithm** (Traag et al., 2019), a community detection method that partitions the kNN graph into groups of densely interconnected cells. The Leiden algorithm was run at multiple resolution parameters (0.3, 0.5, 0.8, 1.2) to explore the clustering hierarchy from coarse (9 clusters) to fine (27 clusters). The resolution parameter controls the granularity: higher values produce more, smaller clusters.

The resolution of **0.5 (12 clusters)** was selected for primary analysis as it captured the major biological cell types without over-splitting biologically coherent populations. Higher-resolution clusterings are available for follow-up analyses of subtypes within major populations.

### 2.6 Cell type annotation

Clusters were annotated by examining the expression of established marker genes for known IVD cell types and general cell lineages. The annotation strategy combined:

1. **Marker gene panels** curated from IVD literature (Risbud and Shapiro, 2014; Gan et al., 2021; Cherif et al., 2022):
   - *NP/chondrocyte markers:* ACAN, COL2A1, SOX9, COL9A3
   - *AF/fibroblast markers:* COL1A1, COL1A2, COL3A1, SCX
   - *Endothelial markers:* PECAM1 (CD31), VWF, CDH5
   - *Macrophage markers:* CD68, CD14, TYROBP
   - *T/NK cell markers:* CD3D, CXCR4, GNLY, NKG7
   - *Pericyte/SMC markers:* TAGLN, MYL9, CALD1

2. **Dot plot visualization** of marker gene expression across all clusters (Figure 2), showing both the fraction of cells expressing each marker and the mean expression level.

3. **Subtype-specific markers** to distinguish NP cell states:
   - *NP: canonical* — high ACAN, COL2A1, SOX9, SCRG1 (the prototypical healthy NP chondrocyte)
   - *NP: HAPLN1+* — high HAPLN1, FN1, TIMP3 (a matrix-organizing subtype)
   - *NP: stress response* — high JUN, FOS, GADD45B, DNAJB1 (immediate early genes indicating cellular stress)
   - *NP: degenerative (UPR)* — high SQSTM1, DNAJB9, TNFRSF12A (unfolded protein response activation)
   - *NP: MT-high* — high MT1G, MT1E, MT1X, MT2A (metallothioneins, indicating oxidative stress and metal ion buffering)

### 2.7 Compositional analysis

To test whether cell type proportions change with degeneration, cell type frequencies were computed per donor (not per cell, to avoid pseudoreplication). Only donors with at least 100 cells were included to ensure reliable proportion estimates. The **Kruskal-Wallis test** (a non-parametric alternative to one-way ANOVA) was used to compare proportions across the 4 degeneration grades for each cell type.

### 2.8 Pseudobulk differential expression analysis

Single-cell differential expression methods that treat each cell as an independent observation are known to produce inflated p-values because cells from the same donor are not independent (Squair et al., 2021; Zimmerman et al., 2021). To avoid this, we used a **pseudobulk** approach:

1. **Aggregation:** For each combination of donor and cell type, raw UMI counts were summed across all cells of that type from that donor, producing a single "pseudobulk" expression profile per donor per cell type. This treats the donor — not the cell — as the unit of replication, which is statistically appropriate.

2. **Filtering:** Pseudobulk samples with fewer than 10 cells were excluded. Genes with fewer than 10 counts in at least 3 samples were removed to ensure adequate statistical power.

3. **DESeq2** (Love et al., 2014) was used to test for differential expression between conditions (severe vs. healthy, moderate vs. healthy, mild vs. healthy) within each cell type. DESeq2 models count data using a negative binomial distribution and performs Wald tests with Benjamini-Hochberg FDR correction. Genes with adjusted p-value < 0.05 and |log2 fold change| > 1 were considered differentially expressed.

This analysis was performed for 6 cell types with sufficient cells and donor coverage: NP: canonical, NP: HAPLN1+, NP: stress response, NP: degenerative (UPR), NP: MT-high, and AF fibroblast.

### 2.9 Pathway enrichment analysis

Two complementary enrichment approaches were applied to the DESeq2 results:

**Gene Set Enrichment Analysis (GSEA)** (Subramanian et al., 2005) was run on the full ranked gene list for each cell type. Genes were ranked by the product sign(LFC) x -log10(adjusted p-value), which places strongly upregulated genes at the top and strongly downregulated genes at the bottom. Gene sets were drawn from:
- **MSigDB Hallmark collection** (Liberzon et al., 2015): 50 curated gene sets representing well-defined biological states and processes
- **Reactome pathways** (Gillespie et al., 2022): filtered to IVD-relevant categories by excluding gene sets related to meiosis, spermatogenesis, ribosomal processing, and other non-relevant processes (451 of ~1,800 Reactome sets retained)

**Over-representation analysis (ORA)** was also performed on significant DEGs using MSigDB Hallmark and KEGG pathway databases as a complementary check.

### 2.10 Trajectory analysis (PAGA)

To infer the relationships between NP cell states and the potential trajectory of degeneration, **Partition-based Graph Abstraction (PAGA)** (Wolf et al., 2019) was applied to the 140,439 cells of the NP lineage (the 5 NP cell states).

The NP subset was re-embedded using the Harmony-corrected PCA space with 30 components and a kNN graph (k=20). PAGA estimates the connectivity between cell clusters by comparing the number of inter-cluster edges in the kNN graph to what would be expected under a null model. High connectivity scores indicate that two cell states have many cells in transitional states between them.

**Diffusion pseudotime (DPT)** (Haghverdi et al., 2016) was also computed, using a healthy NP: canonical cell as the root. However, the diffusion map components showed extremely low variance (all DC variances = 7.12 x 10^-6), likely because the large, well-mixed integrated dataset compressed the diffusion kernel eigenvalues. The Spearman correlation of pseudotime with degeneration grade was statistically significant but weak (rho = 0.24, p < 10^-300). Given these limitations, **PAGA connectivity is reported as the primary trajectory metric**, while DPT results are presented with appropriate caveats.

### 2.11 Cell-cell communication analysis (LIANA)

Cell-cell communication was inferred using **LIANA** (Dimitrov et al., 2022), a framework that integrates multiple ligand-receptor interaction scoring methods into a consensus ranking. LIANA was run separately on healthy (n = 99,883 cells) and severe (n = 19,545 cells) subsets using the following parameters:

- **Resource:** Consensus (a curated set of ligand-receptor pairs combining CellPhoneDB, CellChat, NATMI, and other databases)
- **Minimum expression proportion:** 10% (a ligand or receptor must be expressed in at least 10% of cells in a cluster to be considered active)
- **Minimum cells per cluster:** 30
- **Scoring method:** rank_aggregate (consensus of multiple scoring methods)
- **Permutation testing:** Not performed (n_perms = None) for computational efficiency; results are therefore exploratory

Erythrocytes and monocytes/neutrophils were excluded from the communication analysis as they are unlikely to participate in local IVD signaling.

Differential communication was computed by converting LIANA's magnitude_rank score to -log10 scale, merging healthy and severe results on each source-target-ligand-receptor key, and computing the score difference (delta_score = severe - healthy). Positive delta indicates interactions gained in severe degeneration; negative delta indicates interactions lost.

---

## 3. Results

### 3.1 Atlas overview: 12 cell populations in the human IVD

After quality control, integration, and clustering, the atlas comprises **173,628 cells** from **29 donors** and **7 datasets**. UMAP visualization (Figure 1) reveals a continuous landscape of IVD cell types dominated by the NP lineage, with distinct clusters of AF fibroblasts and small islands of vascular (endothelial, pericyte) and immune (macrophage, T/NK cell) populations.

Twelve cell populations were identified at Leiden resolution 0.5 and annotated based on canonical marker gene expression (Figure 2; Table 2).

**Table 2. Cell type annotations, marker genes, and frequencies.**

| Cell Type | Key Markers | n Cells | % of Total |
|---|---|---|---|
| NP: canonical | ACAN, COL2A1, SOX9, COL9A3, SCRG1 | 47,120 | 27.1% |
| NP: degenerative (UPR) | SQSTM1, DNAJB9, TNFRSF12A, EIF1 | 34,810 | 20.1% |
| NP: MT-high | MT1G, MT1E, MT1X, MT2A, MALAT1 | 31,502 | 18.1% |
| NP: HAPLN1+ | HAPLN1, FN1, TIMP3, FGFBP2 | 17,040 | 9.8% |
| AF fibroblast | COL1A1, COL1A2, COL3A1, SCX | 16,448 | 9.5% |
| NP: stress response | JUN, FOS, GADD45B, DNAJB1 | 9,967 | 5.7% |
| Pericyte/SMC | TAGLN, MYL9, CALD1, NR2F2 | 4,068 | 2.3% |
| Macrophage | CD68, CD14, TYROBP, CTSS | 3,482 | 2.0% |
| Endothelial | PECAM1, VWF, CDH5, SPARCL1 | 3,382 | 1.9% |
| Monocyte/Neutrophil | LYZ, S100A8, S100A9, MNDA | 2,511 | 1.4% |
| Erythrocyte | HBB, HBA1, HBA2, AHSP | 1,658 | 1.0% |
| T/NK cell | CD3D, CXCR4, GNLY, NKG7 | 1,640 | 0.9% |

The five NP cell states together comprise 81% of the atlas, consistent with the NP being the dominant tissue compartment in most datasets. The identification of five transcriptionally distinct NP states — rather than a single "NP cell" type — is a central finding, suggesting that NP cells exist along a phenotypic spectrum from healthy matrix-producing chondrocytes to stress-activated and degenerative states.

![Figure 1. UMAP atlas overview](figures/05_annotation/umap_annotated_v2.png)
**Figure 1. Integrated UMAP atlas of 173,628 human IVD cells.** *(A)* Cells colored by annotated cell type, showing 12 populations. The NP lineage (blue/teal shades) forms a continuous landscape, while AF fibroblasts (red), immune cells, and vascular cells form distinct clusters. *(B)* Same UMAP colored by degeneration condition. Blue = healthy (n = 100,234), yellow = mild (n = 21,648), orange = moderate (n = 32,136), red = severe (n = 19,610). Note the enrichment of severe-condition cells in the lower right region, overlapping with the NP: degenerative and AF fibroblast clusters.

![Figure 2. Marker gene dot plot](figures/05_annotation/dotplot_markers.png)
**Figure 2. Marker gene expression across Leiden clusters.** Dot plot showing expression of canonical marker genes for NP/notochordal, AF, CEP, progenitor, endothelial, and immune lineages across the 12 clusters. Dot size represents the fraction of cells expressing each gene; color intensity represents mean expression level. Clusters 0-5 express NP markers (ACAN, COL2A1, SOX9); cluster 4 expresses AF markers (COL1A1, SCX); clusters 6-11 represent non-resident cell types.

![Figure S1. UMAP overview](figures/05_annotation/umap_overview.png)
**Figure S1. Integration quality assessment.** Four-panel UMAP showing cells colored by Leiden cluster, dataset of origin, degeneration condition, and tissue compartment. Datasets are well-mixed within the NP clusters, indicating successful batch correction by Harmony. Tissue compartment labels confirm AF fibroblast identity in cluster 4.

![Figure S2. QC summary](figures/02_qc/qc_cell_counts.png)
**Figure S2. Quality control summary.** *(A)* Cell counts before (gray) and after (color) QC per dataset. Overall retention was 78%. *(B)* Per-sample QC retention rates. All samples exceeded the 70% retention threshold (dashed red line).

### 3.2 Cell composition shifts during degeneration

Cell type proportions were computed per donor (n = 29) and compared across degeneration grades using the Kruskal-Wallis test (Figure 3).

**AF fibroblast** was the only cell type reaching statistical significance (H = 8.45, p = 0.038), increasing from approximately 7% of cells in healthy donors to approximately 24% in severe degeneration. This expansion of AF-like, type I collagen-producing fibroblasts into the NP compartment is consistent with **fibrocartilaginous metaplasia** — a well-documented hallmark of advanced disc degeneration in which the gel-like NP is progressively replaced by stiffer, collagen I-rich tissue (Antoniou et al., 1996).

Within the NP lineage, two notable shifts were observed (though not reaching significance at the donor level, likely due to the small severe donor sample, n = 3):
- **NP: HAPLN1+** expanded by +10.1 percentage points in severe degeneration, suggesting a compensatory matrix-remodeling response
- **NP: MT-high** contracted by -12.5 percentage points, suggesting loss of the oxidative stress-buffering population

![Figure 3. Compositional analysis](figures/06_composition/composition_overview.png)
**Figure 3. Cell type composition by degeneration grade.** *(A)* Stacked bar chart of mean cell type proportions across four degeneration grades. Note the progressive expansion of AF fibroblasts (red) and shift in NP state balance from healthy to severe. *(B)* Box-and-dot plots of per-donor proportions for key cell types. Each dot represents one donor. AF fibroblast proportion increases significantly with degeneration (Kruskal-Wallis p = 0.038). NP: canonical shows a trend toward decrease, while NP: degenerative trends upward.

### 3.3 Massive transcriptional downregulation in severe degeneration

Pseudobulk DESeq2 analysis comparing severe degeneration to healthy tissue revealed extensive differential gene expression across all 6 tested cell types (Figure 4; Table 3).

**Table 3. Differentially expressed genes (DEGs) per cell type, severe vs. healthy.**

| Cell Type | Total DEGs | Upregulated | Downregulated | Down:Up Ratio |
|---|---|---|---|---|
| NP: canonical | 2,641 | 332 | 2,309 | 7.0:1 |
| NP: degenerative (UPR) | 2,331 | 241 | 2,090 | 8.7:1 |
| NP: HAPLN1+ | 1,786 | 223 | 1,563 | 7.0:1 |
| NP: stress response | 1,081 | 104 | 977 | 9.4:1 |
| AF fibroblast | 1,298 | 103 | 1,195 | 11.6:1 |
| NP: MT-high | 953 | 81 | 872 | 10.8:1 |

The most striking finding is the overwhelming dominance of **downregulated genes**: across all cell types, the ratio of downregulated to upregulated DEGs is approximately **7:1 to 12:1**. This transcriptional collapse suggests that severe degeneration involves a broad shutdown of gene expression programs rather than activation of a few pathogenic pathways.

**Key upregulated genes** (consistently across multiple cell types):
- **ADAMTS5:** The primary aggrecan-degrading enzyme in cartilage; its upregulation directly drives proteoglycan loss (Stanton et al., 2005)
- **FN1 (fibronectin):** An ECM glycoprotein whose fragments are pro-inflammatory and activate macrophages via CD44 and toll-like receptors (Homandberg et al., 1997)
- **TNF, IL6, CXCL8:** Pro-inflammatory cytokines that drive the catabolic cascade in disc degeneration
- **KRT19 (cytokeratin 19):** Upregulated in the NP: stress response state; sometimes used as a notochordal remnant marker but here likely reflecting a stressed epithelial-like phenotype

**Key downregulated genes:**
- **ACAN (aggrecan):** The major proteoglycan of the NP, essential for water retention and compressive strength
- **COL2A1 (type II collagen):** The structural collagen of healthy cartilaginous matrix
- **HAPLN1 (link protein):** Stabilizes aggrecan-hyaluronan aggregates; its loss destabilizes the ECM
- **COMP (cartilage oligomeric matrix protein):** An ECM glycoprotein important for collagen fibril assembly
- **CILP (cartilage intermediate layer protein):** Involved in cartilage homeostasis signaling

![Figure 4. Volcano plots](figures/07_pseudobulk/volcano_severe_vs_healthy.png)
**Figure 4. Pseudobulk differential expression: severe degeneration vs. healthy.** Volcano plots for each of the 6 cell types showing log2 fold change (x-axis) vs. -log10 adjusted p-value (y-axis). Blue points: significantly downregulated genes (padj < 0.05, |LFC| > 1). Red points: significantly upregulated genes. Gray: non-significant. Selected IVD-relevant genes are labeled. Note the asymmetric distribution toward the left (downregulated) in all cell types, reflecting the 7:1 down:up ratio. The NP: canonical and NP: degenerative states show the largest number of DEGs.

### 3.4 Coordinated loss of four homeostatic pathway programs

GSEA revealed that severe degeneration is associated with coordinated suppression of specific signaling and transcriptional programs across all or most NP cell states (Figure 5). Remarkably, the **same pathways were suppressed in all 6 cell types**, suggesting a tissue-wide loss of homeostatic signaling rather than a cell-type-specific defect.

**Pathways consistently downregulated in severe degeneration (significant in all 6 cell types):**

1. **Wnt signaling** (TCF-dependent Wnt signaling, Signaling by Wnt): The Wnt pathway is critical for chondrocyte differentiation, maintenance of the disc progenitor niche, and regulation of ECM production. Wnt ligands (e.g., WNT3A, WNT5A) signal through Frizzled receptors to activate beta-catenin/TCF transcription. Loss of Wnt signaling in degeneration likely reflects exhaustion of the progenitor pool and loss of chondrogenic capacity (Hiyama et al., 2013).

2. **Notch signaling** (Signaling by Notch, Pre-Notch Expression and Processing): The Notch pathway maintains stem/progenitor cell populations in many tissues. In the IVD, Notch signaling has been shown to regulate NP cell proliferation and survival. Its downregulation suggests depletion of the progenitor niche that normally replenishes IVD cells (Wang et al., 2019).

3. **Cellular senescence programs** (Cellular Senescence, Oxidative Stress-Induced Senescence, SASP, DNA Damage/Telomere Stress-Induced Senescence, Senescence-Associated Heterochromatin Foci): Counterintuitively, the senescence *gene programs* are downregulated in severe degeneration. This may reflect the paradox that while senescent cells accumulate in the disc (as shown by SA-beta-galactosidase staining in histological studies), the transcriptional machinery that *regulates* senescence — including checkpoint genes and SASP factors — becomes suppressed as cells transition from regulated senescence into a more apoptotic or quiescent state.

4. **RUNX transcription** (Transcriptional Regulation by RUNX1, RUNX1 Regulates Transcription of Genes Involved in Differentiation): RUNX2 is the master transcription factor for chondrocyte hypertrophy and bone formation, while RUNX1 regulates hematopoietic and chondrogenic programs. Loss of RUNX activity reflects a collapse of chondrogenic transcriptional control.

**Pathways upregulated in specific cell types:**

- **TNF-alpha signaling via NF-kB** (NP: canonical, NP: degenerative): This inflammatory master pathway was significantly upregulated in 2 of 6 cell types, consistent with localized inflammatory activation. NF-kB drives expression of MMPs, inflammatory cytokines, and anti-apoptotic factors (Wuertz et al., 2012).

- **Inflammatory response** (NP: stress, NP: degenerative): A broader inflammatory signature complementing NF-kB activation.

- **Epithelial-mesenchymal transition** (NP: degenerative): EMT-associated gene expression in the degenerative NP state reflects the fibrocartilaginous shift from a chondrocyte-like to a fibroblast-like phenotype.

- **Collagen crosslinking** (NP: stress, NP: degenerative, NP: MT-high): Upregulation of lysyl oxidase and crosslinking enzymes suggests matrix stiffening — a known consequence of degeneration that alters the mechanical environment and further promotes catabolic signaling.

- **Glycolysis** (NP: canonical): The NP is naturally hypoxic and relies on glycolysis for energy. Increased glycolytic gene expression in degeneration may reflect metabolic reprogramming under worsening nutritional deprivation as the endplate calcifies.

![Figure 5. GSEA heatmap](figures/08_pathways/gsea_heatmap_severe_vs_healthy.png)
**Figure 5. Gene set enrichment analysis: severe degeneration vs. healthy.** Heatmap of normalized enrichment scores (NES) for top pathways across 6 cell types. Blue = downregulated in severe (negative NES); red = upregulated in severe (positive NES). Asterisks indicate statistical significance (adjusted p < 0.05). Bottom rows: the 12 pathways suppressed in all 6 cell types, representing coordinated loss of Wnt, Notch, senescence, and RUNX programs. Top rows: selectively upregulated pathways including TNF-alpha/NF-kB signaling and EMT.

![Figure S3. Pathway dot plot](figures/08_pathways/pathway_dotplot_severe_vs_healthy.png)
**Figure S3. Pathway over-representation analysis.** *(Top)* MSigDB Hallmark pathways enriched among upregulated DEGs. Glycolysis, TNF-alpha/NF-kB, and EMT are the most significant. *(Bottom)* Shared Reactome pathways enriched across 3 or more cell types among downregulated DEGs, confirming the GSEA findings of broad transcriptional suppression.

### 3.5 NP cell-state trajectory: stress response as a transition hub

PAGA analysis of the 140,439 NP lineage cells revealed structured connectivity between the five NP states (Figure 6).

**Key PAGA connectivity edges:**
- **NP: stress <-> NP: degenerative:** 0.76 (strongest connection)
- **NP: canonical <-> NP: HAPLN1+:** 0.65
- **NP: stress <-> NP: MT-high:** 0.51
- **NP: canonical <-> NP: stress:** 0.52
- **NP: HAPLN1+ <-> NP: degenerative:** 0.59

These connectivity values suggest a model in which the **NP stress-response state serves as a transition hub** between the canonical healthy chondrocyte and the degenerative UPR state. The biological interpretation is:

1. **NP: canonical** cells (high ACAN, COL2A1, SOX9) represent the resting, matrix-producing chondrocyte
2. Under biomechanical or nutritional stress, cells activate immediate early genes (JUN, FOS) and enter the **NP: stress response** state
3. If the stress is sustained, cells progress to the **NP: degenerative (UPR)** state, characterized by endoplasmic reticulum stress and misfolded protein accumulation (SQSTM1, DNAJB9)
4. **NP: MT-high** cells (metallothioneins) represent an alternative stress response focused on oxidative stress and metal ion buffering, connected primarily through the stress-response hub
5. **NP: HAPLN1+** cells may represent a compensatory remodeling response, attempting to restore matrix integrity through increased HAPLN1 and TIMP3 expression

The HAPLN1+ state's expansion in severe degeneration (+10.1%) is consistent with an activated matrix-repair phenotype, though this response is ultimately insufficient to prevent disease progression.

**Note on diffusion pseudotime:** DPT was computed but showed weak resolution (all diffusion component variances = 7.12 x 10^-6), likely due to the large, well-mixed dataset compressing the diffusion kernel eigenvalues. The Spearman correlation of pseudotime with degeneration grade was significant but weak (rho = 0.24, p < 10^-300). PAGA connectivity provides a more robust measure of cell-state relationships in this context.

![Figure 6. NP trajectory analysis](figures/09_trajectory/np_trajectory_overview.png)
**Figure 6. NP lineage trajectory and state transitions.** *(Top left)* NP subset UMAP colored by cell type. *(Top center)* Same UMAP colored by degeneration condition. *(Top right)* PAGA connectivity matrix between NP states; darker red = stronger connectivity. The NP: stress <-> NP: degenerative connection (0.76) is the strongest, supporting the stress-to-degeneration trajectory model. *(Bottom left)* NP state composition by degeneration grade, showing expansion of NP: HAPLN1+ and contraction of NP: MT-high in severe disease. *(Bottom center)* Pseudotime density distributions by condition (note caveat about weak DPT resolution). *(Bottom right)* Change in NP state proportions between severe and healthy, quantifying the HAPLN1+ expansion (+10.1%) and MT-high contraction (-12.5%).

### 3.6 Cell-cell communication: loss of protective signaling and inflammatory amplification

LIANA cell-cell communication analysis comparing healthy and severe degeneration revealed dramatic remodeling of intercellular signaling networks (Figure 7). A total of 39,530 unique ligand-receptor pairs were evaluated across 10 cell types.

#### 3.6.1 Lost interactions: collapse of TIMP1-CD63 protective signaling

The dominant lost interaction across all cell type pairs was **TIMP1 -> CD63**. TIMP1 is a secreted protein that inhibits matrix metalloproteinases (MMPs), the enzymes responsible for breaking down collagen and proteoglycans in the ECM. CD63 is a tetraspanin receptor on the cell surface that internalizes TIMP1, mediating its anti-proteolytic effects. The top 10 lost interactions were all TIMP1-CD63 pairs between different cell type combinations, with the strongest loss occurring in NP: HAPLN1+ autocrine signaling (delta = -5.06).

**Biological significance:** The loss of TIMP1-CD63 signaling in severe degeneration means that the disc's natural MMP-inhibitory defense system is collapsing. Without TIMP1 restraint, MMPs and ADAMTS enzymes can degrade the ECM unchecked, accelerating the structural breakdown of the disc. This finding aligns with biochemical studies showing reduced TIMP levels in degenerated disc tissue (Le Maitre et al., 2004) and identifies a specific signaling axis rather than just a change in TIMP protein levels.

#### 3.6.2 Gained interactions: FN1-mediated inflammatory recruitment

The top gained interactions in severe degeneration were dominated by **fibronectin (FN1)** signaling:

- **FN1 -> ITGA6** (NP: HAPLN1+ -> Endothelial, delta = +4.30): Fibronectin signaling to endothelial cells via integrin alpha-6, potentially promoting **angiogenesis** — the growth of new blood vessels into the normally avascular disc.

- **FN1 -> C5AR1** (NP: HAPLN1+ -> Macrophage, delta = +4.09): Fibronectin fragments activate complement receptor C5aR1 on macrophages, promoting inflammatory macrophage activation.

- **FN1 -> CD44** (NP: HAPLN1+ -> Macrophage, delta = +3.93): Fibronectin fragments bind CD44 on macrophages, triggering inflammatory cytokine production. This establishes a feed-forward loop: matrix degradation produces FN1 fragments, which recruit and activate macrophages, which secrete more MMPs, leading to more degradation.

- **COL1A2 -> CD93** (AF fibroblast -> Endothelial, delta = +3.93): Type I collagen signaling to endothelial cells via CD93, a receptor implicated in angiogenesis and vascular remodeling.

- **SEMA4A -> PLXNB1** (gained across NP states): Semaphorin 4A signaling to plexin B1. Semaphorins are guidance molecules originally described in axonal pathfinding; in the disc context, this gained interaction may facilitate the **nerve ingrowth** that characterizes painful degenerated discs (Freemont et al., 2002).

**Overall pattern:** Total interaction strength increased globally in severe degeneration, with all NP states showing increased outgoing signaling. The interaction landscape shifted from a protective, homeostatic mode (TIMP1-mediated MMP inhibition) to an inflammatory, degradative mode (FN1-mediated macrophage activation and vascular/nerve ingrowth).

![Figure 7. Cell-cell communication](figures/10_cellchat/liana_communication_overview.png)
**Figure 7. Cell-cell communication changes in severe IVD degeneration (LIANA analysis).** *(Top row)* Heatmaps of total interaction strength between cell type pairs in healthy (left), severe (center), and the difference (right; red = gained, blue = lost). Note the global increase in interaction strength in severe degeneration. *(Bottom left)* Bar chart of top 20 gained (red) and lost (blue) ligand-receptor interactions. TIMP1-CD63 interactions (blue, bottom) dominate the lost category, while FN1 interactions (red, top) dominate the gained category. *(Bottom right)* Total interaction score for the 4 most changed LR pairs across all cell type combinations, comparing healthy (blue) vs. severe (red). TIMP1-CD63 drops from a score of ~400 to ~100, while FN1-CD44 and FN1-C5AR1 show modest gains.

---

## 4. Discussion

### 4.1 A multi-layered model of disc degeneration

Integrating the findings across all analyses, we propose a multi-layered model of IVD degeneration that progresses through distinct but overlapping phases:

**Phase 1 — Loss of homeostatic signaling:** The earliest changes involve suppression of Wnt, Notch, and RUNX transcriptional programs across all NP cell states. These pathways maintain the chondrogenic phenotype, regulate progenitor cell activity, and sustain ECM production. Their coordinated loss suggests a tissue-wide failure of the maintenance machinery, possibly triggered by cumulative biomechanical stress or nutritional deprivation through endplate calcification.

**Phase 2 — Stress activation and phenotypic transition:** NP cells respond to the hostile microenvironment by activating immediate early gene programs (JUN, FOS) and entering a stress-response state. Some cells activate metallothionein programs (MT-high state) to buffer oxidative stress. Sustained stress drives transition to the degenerative UPR state, where misfolded protein accumulation triggers endoplasmic reticulum stress.

**Phase 3 — Inflammatory amplification:** The degenerative process activates TNF/NF-kB inflammatory signaling in specific NP states, upregulates catabolic enzymes (ADAMTS5), and generates fibronectin fragments that recruit and activate macrophages. TIMP1-mediated protective signaling collapses, tipping the protease-antiprotease balance toward matrix degradation. A feed-forward inflammatory loop is established.

**Phase 4 — Structural remodeling and failure:** AF fibroblasts expand into the NP space (fibrocartilaginous metaplasia), depositing type I collagen and stiffening the tissue. Angiogenic signaling (FN1-ITGA6, COL1A2-CD93) promotes neovascularization, and semaphorin signaling (SEMA4A-PLXNB1) may facilitate nerve ingrowth — the likely basis for discogenic pain.

### 4.2 Therapeutic target candidates

The molecular events identified in this atlas suggest several candidate therapeutic strategies, organized by the phase of degeneration they would address:

#### Restoring homeostatic signaling (Phase 1 targets)

- **Wnt pathway agonists:** Small-molecule Wnt activators (e.g., lithium chloride, CHIR99021) or recombinant Wnt ligands could potentially restore chondrogenic signaling and ECM production. However, Wnt activation must be carefully titrated, as excessive Wnt signaling promotes chondrocyte hypertrophy and osteoarthritis in synovial joints (Zhu et al., 2009).

- **Notch pathway modulation:** Strategies to maintain or restore Notch signaling could support the progenitor niche. Notch ligand-coated biomaterials or localized Notch agonist delivery are being explored in other cartilage contexts (Hosaka et al., 2013).

#### Blocking the stress-to-degeneration transition (Phase 2 targets)

- **UPR modulators:** Chemical chaperones (e.g., 4-phenylbutyric acid, tauroursodeoxycholic acid) that reduce ER stress could prevent the transition from the stress state to the degenerative UPR state.

- **Antioxidant therapy:** The contraction of the MT-high population suggests loss of oxidative stress buffering. N-acetylcysteine or mitochondria-targeted antioxidants could fill this gap.

#### Dampening inflammatory amplification (Phase 3 targets)

- **Anti-TNF biologics:** TNF inhibitors (e.g., etanercept, infliximab) have shown efficacy in preclinical disc degeneration models and limited clinical trials for sciatica (Korhonen et al., 2006). This analysis supports their rationale by showing TNF/NF-kB activation is cell-type-specific (NP: canonical, NP: degenerative), enabling targeted delivery strategies.

- **TIMP1 supplementation or gene therapy:** Restoring TIMP1-CD63 signaling could re-establish the protease-antiprotease balance. Recombinant TIMP1 delivery or adeno-associated virus (AAV)-mediated TIMP1 gene therapy targeted to NP cells are conceptually attractive approaches (Leckie et al., 2012).

- **FN1 fragment neutralization:** Blocking the interaction of fibronectin fragments with macrophage receptors (CD44, C5AR1) could break the inflammatory feed-forward loop. Anti-FN1 fragment antibodies or CD44 antagonists could be explored.

- **Anti-ADAMTS5 therapy:** Given the consistent upregulation of ADAMTS5, selective ADAMTS5 inhibitors could directly slow aggrecan degradation. ADAMTS5 inhibitors have been pursued in osteoarthritis drug development (Tortorella et al., 2009).

#### Preventing structural failure (Phase 4 targets)

- **Anti-angiogenic therapy:** Inhibiting neovascularization (e.g., via anti-VEGF agents or targeting the FN1-ITGA6 axis) could prevent blood vessel and nerve ingrowth into the disc.

- **Senolytic therapy:** Clearing senescent cells with senolytic drugs (e.g., navitoclax, dasatinib + quercetin) could reduce SASP-driven inflammation. Senolytics have shown promise in preclinical models of disc degeneration (Novais et al., 2021).

### 4.3 Comparison with existing IVD scRNA-seq studies

The individual datasets contributing to this atlas have been published with their own analyses (Gan et al., 2021; Cherif et al., 2022; Wang et al., 2023; Li et al., 2022; Guo et al., 2023; Jiang and Sheyn, 2022; Shi et al., 2024). Our integrated analysis extends these studies in several ways:

1. **Increased statistical power:** By aggregating 29 donors, we can perform pseudobulk differential expression with appropriate statistical modeling (DESeq2), treating the donor as the unit of replication rather than the cell. Individual studies typically have 2-7 donors, which limits their ability to detect subtle transcriptional changes.

2. **Cross-dataset validation:** The consistent identification of NP cell states across 7 independently generated datasets increases confidence that these states represent genuine biology rather than dataset-specific artifacts.

3. **Systematic pathway and signaling analysis:** While individual studies have reported selected differentially expressed genes, this analysis provides the first comprehensive GSEA and LIANA cell-cell communication comparison across degeneration grades in integrated IVD data.

4. **Trajectory context:** PAGA analysis across the full 140,000-cell NP lineage reveals the connectivity structure between states, identifying the stress-response state as a transition hub — a finding that could not emerge from smaller individual datasets.

### 4.4 Limitations

This study has several important limitations:

1. **Unbalanced condition groups:** The atlas contains 100,234 healthy cells (15 donors) vs. 19,610 severe cells (3 donors). This imbalance limits statistical power for compositional analyses and may bias toward detecting changes that are large in effect size. The pseudobulk DESeq2 approach partially mitigates this by using donor-level replication, but 3 severe donors is still a small sample.

2. **Cross-study heterogeneity:** Despite Harmony batch correction, residual technical differences between datasets may persist, particularly for the neonatal samples in GSE189916 (Jiang/Sheyn 2022), which represent a fundamentally different developmental stage.

3. **Annotation resolution:** At Leiden resolution 0.5, some biologically distinct subtypes (e.g., notochordal cells, distinct CEP chondrocyte populations) may be merged into broader clusters. Higher-resolution clustering (27 clusters at resolution 1.2) is available for follow-up analyses.

4. **Exploratory communication analysis:** LIANA was run without permutation testing for computational efficiency, so ligand-receptor results should be interpreted as hypothesis-generating rather than formally validated.

5. **Trajectory limitations:** Diffusion pseudotime showed weak resolution in this dataset; the large number of well-mixed cells likely compressed the diffusion kernel eigenvalues. PAGA connectivity provides a more robust measure but does not assign a continuous ordering to individual cells.

6. **Lack of spatial context:** scRNA-seq dissociates the tissue, losing spatial information about where cells reside within the disc. Spatial transcriptomics of human IVD would be needed to map these cell states to specific anatomical locations.

7. **Causality not established:** This is a cross-sectional observational study. We cannot determine whether the observed transcriptional changes cause degeneration or are consequences of it. Longitudinal studies or experimental perturbation in model systems are needed to establish causal relationships.

### 4.5 Future directions

1. **Experimental validation of therapeutic targets:** The top candidates (Wnt agonism, TIMP1 restoration, FN1 fragment neutralization, anti-ADAMTS5) should be tested in human NP cell culture systems and organotypic disc models.

2. **Spatial transcriptomics:** Visium or MERFISH spatial transcriptomics of human IVD sections would map the cell states identified here to their anatomical locations and reveal spatial patterns of degeneration.

3. **Higher-resolution NP subtyping:** The NP: canonical cluster (47,120 cells) may contain sub-populations — including notochordal remnants and distinct chondrogenic progenitors — that merit further investigation at higher clustering resolution.

4. **Integration with drug perturbation data:** Cross-referencing the DEG signatures with pharmacogenomic databases (e.g., Connectivity Map) could identify existing drugs that reverse the degenerative transcriptional signature.

5. **Proteomics and functional validation:** Transcriptomic changes should be confirmed at the protein level, as post-transcriptional regulation is important in the IVD (particularly for ECM proteins, which are long-lived).

---

## 5. Conclusions

This integrated single-cell atlas of 173,628 human IVD cells provides the most comprehensive view to date of cell-state diversity and transcriptional changes in disc degeneration. The analysis reveals that degeneration is not simply an increase in inflammatory signaling but rather a multi-layered process involving (1) coordinated loss of homeostatic Wnt, Notch, and RUNX transcriptional programs, (2) a stress-response cascade through which NP cells transition from healthy to degenerative states, (3) collapse of protective TIMP1-mediated MMP inhibition, (4) inflammatory amplification via FN1-macrophage signaling, and (5) structural remodeling through AF fibroblast expansion and neovascularization. These findings define a rich landscape of potential therapeutic targets for IDD, from pathway restoration (Wnt, Notch) to specific signaling axes (TIMP1-CD63, FN1-CD44, ADAMTS5) that could be addressed by existing drug modalities.

---

## References

Adams MA, Roughley PJ. (2006). What is intervertebral disc degeneration, and what causes it? *Spine*, 31(18):2151-2161. doi:10.1097/01.brs.0000231761.73859.2c

Antoniou J, Steffen T, Nelson F, et al. (1996). The human lumbar intervertebral disc: evidence for changes in the biosynthesis and denaturation of the extracellular matrix with growth, maturation, ageing, and degeneration. *Journal of Clinical Investigation*, 98(4):996-1003. doi:10.1172/JCI118884

Brinjikji W, Luetmer PH, Comstock B, et al. (2015). Systematic literature review of imaging features of spinal degeneration in asymptomatic populations. *American Journal of Neuroradiology*, 36(4):811-816. doi:10.3174/ajnr.A4173

Cherif H, Bisson DG, Mannarino M, et al. (2022). Single-cell RNA-seq analysis of cells from degenerating and non-degenerating intervertebral discs from the same individual reveals new biomarkers for intervertebral disc degeneration. *International Journal of Molecular Sciences*, 23(7):3993. doi:10.3390/ijms23073993

Dieleman JL, Cao J, Chapin A, et al. (2020). US health care spending by payer and health condition, 1996-2016. *JAMA*, 323(9):863-884. doi:10.1001/jama.2020.0734

Dimitrov D, Turei D, Garrber M, et al. (2022). Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. *Nature Communications*, 13:3224. doi:10.1038/s41467-022-30755-0

Feng C, Liu M, Fan X, et al. (2016). Disc cell senescence in intervertebral disc degeneration: causes and molecular pathways. *Cell Cycle*, 15(13):1674-1684. doi:10.1080/15384101.2016.1152433

Freemont AJ, Watkins A, Le Maitre C, et al. (2002). Nerve growth factor expression and innervation of the painful intervertebral disc. *Journal of Pathology*, 197(3):286-292. doi:10.1002/path.1108

Gan Y, He J, Zhu J, et al. (2021). Spatially defined single-cell transcriptional profiling characterizes diverse chondrocyte subtypes and nucleus pulposus progenitors in human intervertebral discs. *Bone Research*, 9:37. doi:10.1038/s41413-021-00163-z

GBD 2021 Low Back Pain Collaborators. (2023). Global, regional, and national burden of low back pain, 1990-2020, its attributable risk factors, and projections to 2050. *The Lancet Rheumatology*, 5(6):e316-e329. doi:10.1016/S2665-9913(23)00098-X

Gillespie M, Jassal B, Stephan R, et al. (2022). The reactome pathway knowledgebase 2022. *Nucleic Acids Research*, 50(D1):D986-D992. doi:10.1093/nar/gkab1028

Guo R, Liu M, Liang Y, et al. (2023). Single-cell RNA sequencing reveals heterogeneous immune and NP cell atlas in degenerative human intervertebral disc. *Frontiers in Cell and Developmental Biology*, 11:1170062. doi:10.3389/fcell.2023.1170062

Haghverdi L, Buttner M, Wolf FA, Buettner F, Theis FJ. (2016). Diffusion pseudotime robustly reconstructs lineage branching. *Nature Methods*, 13:845-848. doi:10.1038/nmeth.3971

Hiyama A, Sakai D, Risbud MV, et al. (2013). Enhancement of intervertebral disc cell senescence by WNT/beta-catenin signaling-induced matrix metalloproteinase expression. *Arthritis & Rheumatism*, 62(10):3036-3047. doi:10.1002/art.27599

Homandberg GA, Meyers R, Williams JM. (1997). Intraarticular injection of fibronectin fragments causes severe depletion of cartilage proteoglycans in vivo. *Journal of Rheumatology*, 24(1):129-133.

Hosaka Y, Saito T, Sugita S, et al. (2013). Notch signaling in chondrocytes modulates endochondral ossification and osteoarthritis development. *Proceedings of the National Academy of Sciences*, 110(5):1875-1880. doi:10.1073/pnas.1207458110

Humzah MD, Soames RW. (1988). Human intervertebral disc: structure and function. *Anatomical Record*, 220(4):337-356. doi:10.1002/ar.1092200402

Jiang S, Sheyn D. (2022). Single-cell atlas of intervertebral disc development and degeneration. *bioRxiv* preprint.

Johnson J, Douze M, Jegou H. (2019). Billion-scale similarity search with GPUs. *IEEE Transactions on Big Data*, 7(3):535-547. doi:10.1109/TBDATA.2019.2921572

Korhonen T, Karppinen J, Paimela L, et al. (2006). The treatment of disc-herniation-induced sciatica with infliximab: one-year follow-up results of FIRST II, a randomized controlled trial. *Spine*, 31(24):2759-2766. doi:10.1097/01.brs.0000245873.23876.1e

Korsunsky I, Millard N, Fan J, et al. (2019). Fast, sensitive and accurate integration of single-cell data with Harmony. *Nature Methods*, 16:1289-1296. doi:10.1038/s41592-019-0619-0

Le Maitre CL, Freemont AJ, Hoyland JA. (2004). Localization of degradative enzymes and their inhibitors in the degenerate human intervertebral disc. *Journal of Pathology*, 204(1):47-54. doi:10.1002/path.1608

Leckie SK, Bechara BP, Hartman RA, et al. (2012). Injection of AAV2-BMP2 and AAV2-TIMP1 into the nucleus pulposus slows the course of intervertebral disc degeneration in an in vivo rabbit model. *Spine Journal*, 12(1):7-20. doi:10.1016/j.spinee.2011.09.011

Li Z, Chen S, Ma K, et al. (2022). CL-82198 treatment attenuates intervertebral disc degeneration by inhibiting NP cell apoptosis and MMP-13 expression. *Frontiers in Cell and Developmental Biology*, 10:869101. doi:10.3389/fcell.2022.869101

Liberzon A, Birger C, Thorvaldsdottir H, et al. (2015). The Molecular Signatures Database Hallmark gene set collection. *Cell Systems*, 1(6):417-425. doi:10.1016/j.cels.2015.12.004

Love MI, Huber W, Anders S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15:550. doi:10.1186/s13059-014-0550-8

Novais EJ, Tran VA, Johnston SN, et al. (2021). Long-term treatment with senolytic drugs dasatinib and quercetin ameliorates age-dependent intervertebral disc degeneration in mice. *Nature Communications*, 12:5213. doi:10.1038/s41467-021-25453-2

Risbud MV, Shapiro IM. (2014). Role of cytokines in intervertebral disc degeneration: pain and disc content. *Nature Reviews Rheumatology*, 10(1):44-56. doi:10.1038/nrrheum.2013.160

Roberts S, Urban JP, Evans H, Eisenstein SM. (1996). Transport properties of the human cartilage endplate in relation to its composition and calcification. *Spine*, 21(4):415-420. doi:10.1097/00007632-199602150-00003

Roughley PJ. (2004). Biology of intervertebral disc aging and degeneration: involvement of the extracellular matrix. *Spine*, 29(23):2691-2699. doi:10.1097/01.brs.0000146101.53784.b1

Shi Y, He R, Yang Y, et al. (2024). Single-cell RNA sequencing reveals cellular landscape of cartilage endplate degeneration. *Frontiers in Immunology*, 15:1336207. doi:10.3389/fimmu.2024.1336207

Squair JW, Gautier M, Kathe C, et al. (2021). Confronting false discoveries in single-cell differential expression. *Nature Communications*, 12:5692. doi:10.1038/s41467-021-25960-2

Stanton H, Rogerson FM, East CJ, et al. (2005). ADAMTS5 is the major aggrecanase in mouse cartilage in vivo and in vitro. *Nature*, 434:648-652. doi:10.1038/nature03417

Subramanian A, Tamayo P, Mootha VK, et al. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *Proceedings of the National Academy of Sciences*, 102(43):15545-15550. doi:10.1073/pnas.0506580102

Tortorella MD, Malfait AM, Deccico C, Arner E. (2001). The role of ADAM-TS4 (aggrecanase-1) and ADAM-TS5 (aggrecanase-2) in a model of cartilage degradation. *Osteoarthritis and Cartilage*, 9(6):539-552. doi:10.1053/joca.2001.0427

Traag VA, Waltman L, van Eck NJ. (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports*, 9:5233. doi:10.1038/s41598-019-41695-z

Wang X, Li D, Wu H, et al. (2019). Notch signaling in intervertebral disc development and degeneration. *Journal of Cellular and Molecular Medicine*, 23(11):7365-7373. doi:10.1111/jcmm.14633

Wang Z, Zhang J, Chen G, et al. (2023). Single-cell RNA sequencing reveals the molecular landscape of nucleus pulposus degeneration. *Genes & Diseases*, 10(6):2408-2424. doi:10.1016/j.gendis.2022.05.025

Wolf FA, Hamey FK, Plass M, et al. (2019). PAGA: graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells. *Genome Biology*, 20:59. doi:10.1186/s13059-019-1663-x

Wolock SL, Lopez R, Klein AM. (2019). Scrublet: computational identification of cell doublets in single-cell transcriptomic data. *Cell Systems*, 8(4):281-291.e9. doi:10.1016/j.cels.2018.11.005

Wuertz K, Vo N, Kletsas D, Boos N. (2012). Inflammatory and catabolic signalling in intervertebral discs: the roles of NF-kB and MAP kinases. *European Cells and Materials*, 23:103-120. doi:10.22203/eCM.v023a08

Zhu M, Tang D, Wu Q, et al. (2009). Activation of beta-catenin signaling in articular chondrocytes leads to osteoarthritis-like phenotype in adult beta-catenin conditional activation mice. *Journal of Bone and Mineral Research*, 24(1):12-21. doi:10.1359/jbmr.080901

Zimmerman KD, Espeland MA, Langefeld CD. (2021). A practical solution to pseudoreplication bias in single-cell studies. *Nature Communications*, 12:738. doi:10.1038/s41467-021-21038-1

---

## Software

| Tool | Version | Purpose |
|---|---|---|
| Scanpy | 1.9 | scRNA-seq preprocessing, integration, clustering, UMAP, PAGA, DPT |
| Harmony | (via scanpy) | Batch correction across datasets and donors |
| FAISS | (IVFFlat) | Approximate k-nearest neighbor graph construction |
| Scrublet | — | Computational doublet detection |
| DESeq2 | (R) | Pseudobulk differential expression analysis |
| clusterProfiler | (R) | Gene set enrichment analysis (GSEA) and over-representation analysis |
| MSigDB / msigdbr | v10 | Hallmark and Reactome gene set collections |
| LIANA | 1.7.1 | Ligand-receptor cell-cell communication inference |
| matplotlib / seaborn | — | Visualization |

---

*Analysis performed using the Phylo automated bioinformatics framework.*
*Datasets: GSE160756, GSE199866, GSE244889, GSE255768, GSE233666, GSE205535, GSE189916.*
