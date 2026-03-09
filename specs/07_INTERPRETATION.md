# Module 07: Biological Interpretation

## Objective

Translate differential expression results into biological insight: what pathways and programs are altered in IVD degeneration, what are the effector molecules, and what regulatory networks drive the changes? Specifically address the project's interest in pain-associated changes.

## Inputs

- DE results from `results/differential/`
- Integrated and annotated AnnData objects from Module 05 (with `obs['cell_type']` from de novo annotation)
- Curated gene sets (see below)

## Outputs

- `results/interpretation/pathway_enrichment/{cell_type}_{comparison}.tsv`
- `results/interpretation/gene_programs.tsv` — coordinated gene programs identified
- `results/interpretation/pain_genes.tsv` — pain-associated DE genes with context
- `results/interpretation/regulatory_networks/` — GRN results
- `results/interpretation/interpretation_report.html`

### Notebook: `notebooks/07_interpretation.ipynb`

Produced after pathway, GRN, and pain gene analyses. Contains:

*Pathway section:*
- Dot plots: top enriched GO/KEGG/Reactome terms per cell type per comparison (up and down separately)
- REVIGO-style semantic clustering of enriched terms to reduce redundancy
- Comparison of enriched pathways across cell types (shared vs. cell-type-specific responses to degeneration)

*GRN section:*
- Heatmap: regulon activity (AUCell scores) across cell types and conditions
- Top regulons with differential activity between healthy and degenerated
- Network visualization of key TF → target gene relationships for top regulons

*Pain section:*
- Table: all pain-associated DE genes with cell type, direction, fold change, pain pathway
- Dot plot: expression of pain gene panel across cell types and conditions
- Schematic: proposed model of how degenerated IVD cells promote nerve ingrowth and pain signaling (this will need human input to draft but the notebook provides the data support)

**Manuscript mapping:** Figure 4: Pathway enrichment overview. Figure 5: Pain-associated gene expression and signaling model. Supplementary Figure S4: Full regulon analysis. Discussion section on pain mechanisms.

## Part 1: Pathway and Gene Set Enrichment

### Method

Run gene set enrichment analysis on DE results:

**Primary tool:** `decoupler` (supports multiple enrichment methods in a unified framework)

**Gene set databases:**
- GO Biological Process, Molecular Function, Cellular Component
- KEGG pathways
- Reactome pathways
- MSigDB Hallmark gene sets
- Custom IVD-relevant gene sets (see below)

**Approach:**
1. For each DE result (cell type × comparison), run over-representation analysis (ORA) on significant DE genes (up and down separately)
2. Also run GSEA (Gene Set Enrichment Analysis) on the ranked gene list (ranked by log2FC × -log10(p-value)) for a threshold-free approach
3. Report enriched terms with FDR < 0.05
4. Cluster redundant terms using semantic similarity (REVIGO or similar) to reduce noise

### Custom IVD-relevant gene sets

Compile gene sets for IVD-specific processes. These should be assembled from literature before running the analysis:

- **Extracellular matrix homeostasis:** collagens (COL1A1, COL2A1, COL9A1, COL11A1, etc.), proteoglycans (ACAN, VCAN, HAPLN1, BGN, DCN), glycoproteins (FN1, COMP, THBS1)
- **Matrix degradation:** MMPs (MMP1, MMP3, MMP7, MMP9, MMP13), ADAMTS (ADAMTS4, ADAMTS5), TIMPs (TIMP1, TIMP2, TIMP3)
- **Inflammatory signaling:** IL1B, IL6, TNF, CXCL8, CCL2, CCL5, PTGS2/COX2, NFkB pathway members
- **Pain mediators:** See Part 3 below
- **Cellular senescence:** CDKN1A/p21, CDKN2A/p16, TP53, SERPINE1, GLB1, senescence-associated secretory phenotype (SASP) components
- **Autophagy:** BECN1, MAP1LC3B, ATG5, ATG7, SQSTM1
- **Apoptosis:** BAX, BCL2, CASP3, CASP9, FASLG
- **Oxidative stress / hypoxia:** HIF1A, VEGFA, SOD2, NOS2, GPX1
- **Mechanotransduction:** YAP1, TAZ/WWTR1, PIEZO1, TRPV4, integrins
- **Notochordal/developmental:** T/TBXT, SHH, NOG, WNT signaling, Hedgehog signaling

## Part 2: Gene Regulatory Networks

### Method

**Primary tool:** pySCENIC (Single-Cell Regulatory Network Inference and Clustering)

**Steps:**
1. Run GRNBoost2 to infer transcription factor (TF) → target gene links
2. Prune using RcisTarget (motif enrichment near target genes) to retain direct regulatory links
3. Score regulon activity per cell using AUCell
4. Compare regulon activity between conditions

**Scope:** Run SCENIC on the resident IVD cell populations (NP and AF separately). Immune cells can be included if immune-IVD interactions are of interest.

**Expected outputs:**
- List of active regulons per cell type and condition
- Regulons whose activity changes significantly between healthy and degenerated
- TF candidates that may drive the degenerative phenotype

**Practical notes:**
- SCENIC is computationally expensive. Run on a subset of cells if full dataset is too large (e.g., downsample to 5000 cells per condition per cell type).
- Use the human genome databases (hg38 gene annotations, JASPAR/TRANSFAC motifs)

## Part 3: Pain-Associated Gene Analysis

### Rationale

A specific project goal is to interpret DE genes for their contributions to pain. IVD degeneration is a major cause of low back pain, but the molecular mechanisms linking disc pathology to pain are incompletely understood.

### Pain gene sets

Compile from literature and databases:

**Nociception/pain signaling:**
- Ion channels: TRPV1, TRPV4, TRPA1, SCN9A/Nav1.7, SCN10A/Nav1.8, SCN11A/Nav1.9, ASIC1-3, P2RX3
- Neuropeptides: TAC1 (substance P), CALCA/CALCB (CGRP), NPY, VIP, BDNF, NGF
- Receptors: NTRK1 (TrkA), NTRK2 (TrkB), NGFR (p75NTR), OPRD1, OPRM1

**Neurotrophin signaling:**
- NGF, BDNF, NT3/NTF3, NT4/NTF4
- NTRK1, NTRK2, NTRK3, NGFR

**Nerve growth / neurite outgrowth:**
- Semaphorins (SEMA3A, SEMA3F — typically repulsive), Netrins (NTN1)
- Guidance receptors: NRP1, NRP2, ROBO1

**Inflammatory pain mediators:**
- PTGS2/COX-2, PTGES, PLA2G2A
- IL1B, IL6, TNF, CCL2, CXCL8
- Bradykinin pathway: BDKRB1, BDKRB2, KLK1

**Neovascularization (nerves follow blood vessels into degenerated discs):**
- VEGFA, VEGFB, FGF2, PDGF, ANGPT1/2

### Analysis steps

1. Cross-reference all significant DE genes with the pain gene sets above
2. For each pain-associated DE gene, record: cell type, comparison, direction (up/down), fold change, and the pain pathway it belongs to
3. Look for coordinated upregulation of pain-related programs (e.g., simultaneous upregulation of NGF and its receptor, or co-upregulation of inflammatory mediators and neurotrophins)
4. Check whether pain-associated genes are enriched in specific cell types/states rather than globally upregulated

### Context from literature

The healthy IVD is avascular and aneural (except outer AF). In degeneration, neovascularization and nerve ingrowth occur, which is thought to contribute to discogenic pain. DE genes that promote neovascularization or nerve growth in degenerated disc cells are particularly relevant.

## Automated Validation

- [ ] Pathway enrichment results exist for all DE comparisons
- [ ] At least some expected IVD-relevant pathways are enriched (e.g., ECM organization, inflammatory response in degeneration)
- [ ] SCENIC results are generated (regulon list, AUCell scores)
- [ ] Pain gene cross-reference table is generated
- [ ] Interpretation report HTML is generated
- [ ] No enrichment result has >500 significant terms at FDR<0.05 (suggests a problem with the background or gene list)

## Human Checkpoint

### Review materials
- Top enriched pathways per cell type per comparison
- SCENIC regulon activity heatmap (conditions × regulons)
- Pain-associated DE gene table
- Interpretation report

### Questions for the reviewer
1. Are the enriched pathways consistent with known IVD degeneration biology?
2. Are there novel pathways or TFs that warrant follow-up?
3. Do the pain-associated findings suggest specific cell types as primary contributors to discogenic pain?
4. Are there actionable targets (e.g., druggable genes) among the top DE/regulatory hits?
5. Are there findings that contradict established IVD biology? If so, are they artifacts or genuinely novel?

### Potential plan revisions
- If pain-associated genes are concentrated in immune cells rather than resident cells, this shifts the narrative and may motivate additional immune-focused analysis
- If SCENIC identifies TFs not in the standard IVD literature, literature search should be expanded to assess their relevance
- Findings here may generate new hypotheses that require going back to DE analysis with different cell type groupings or comparisons
