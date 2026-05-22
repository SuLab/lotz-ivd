# A Continuum-Aware Single-Cell Atlas of the Human Intervertebral Disc: Cell-Type-Specific Degenerative Programs and Pain-Associated Molecular Circuits

**Draft 2026-05-22 — tiered_v4 pipeline**

---

## Abstract

The intervertebral disc (IVD) houses transcriptomically continuous cell populations whose remodeling underlies low-back pain, yet existing single-cell studies are individually underpowered and methodologically heterogeneous. We aggregated 12 publicly available human IVD single-cell RNA-seq datasets (78 samples, 57 donors, ~411,000 cells across nucleus pulposus, annulus fibrosus, and cartilaginous endplate) into a single integrated atlas using a tiered analysis strategy: non-resident immune and endothelial cells were integrated by standard scVI, whereas resident chondrocyte- and fibroblast-like cells were integrated under more conservative settings designed to preserve the well-documented anatomical continuum within each compartment. A three-stage annotation pipeline (coarse cell class → compartment-prefixed cell type → sub-state) yielded 19 cell types across the atlas and 27 sub-states, with two contamination tags (red-blood-cell and endothelial-admixed) retained for transparency rather than filtered out. Pseudobulk DESeq2 across 27 statistically powered cell-type × condition contrasts identified 1,823 trustworthy significant DEGs after excluding a single statistically inflated contrast (Macrophage_M2 healthy-vs-severe; 5,659 spurious DEGs from a 7-sample comparison). NP fibrocartilaginous cells dominate the degenerative signal (325 mild-vs-severe DEGs), with extracellular matrix catabolism (MMP, ADAMTS upregulation; COL2A1, ACAN downregulation), NF-κB-driven inflammatory programs, and 18 pain-associated genes (NTN1, IL6, PLA2G2A, VEGFA, PENK, NGFR, BDKRB1/2, PTGS2 among others) cell-type-specifically dysregulated. AF outer cells respond early in mild degeneration (121 DEGs at healthy-vs-mild, dropping to 3 at healthy-vs-severe), suggesting a transient stress response distinct from the progressive NP signature. Compositional shifts did not reach FDR significance, consistent with degeneration being primarily transcriptional rather than cellular. Trajectory analysis confirmed cell-type-specific pseudotime gradients but no robust compartment-level direction. The atlas, code, and analyses are reproducible end-to-end and presented with explicit version-sensitivity caveats.

---

## Introduction

Low back pain is the leading cause of years lived with disability worldwide, and intervertebral disc degeneration is its most common identifiable structural cause. The healthy adult disc is a heterogeneous, avascular, hypocellular tissue comprising the gelatinous nucleus pulposus (NP), the concentric annulus fibrosus (AF), and the cartilaginous endplates (CEP) that interface with vertebral bone. Degeneration involves a complex cascade of extracellular matrix catabolism, cell-state remodeling, neovascularization, and aberrant nerve ingrowth — the latter two contributing directly to the discogenic pain phenotype.

Single-cell RNA sequencing has been applied to the human IVD by multiple groups since 2020, yielding insights into resident cell-type heterogeneity, fibrotic transitions, and immune infiltration. However, each individual study is modest in sample size, restricted to one or two compartments, and uses non-standardized annotation. Integration across studies is therefore both attractive (statistical power, generalizability) and methodologically delicate (batch effects, condition–dataset confounding, the continuous nature of resident cell phenotypes).

A central methodological concern in IVD scRNA-seq is the continuum problem: chondrocyte-like and fibroblast-like resident cells exist on a graded morphological and molecular spectrum within each compartment, and aggressive batch correction can erase the very variation that the analysis aims to characterize. Prior work has either ignored this issue (treating cells as discrete clusters) or applied a single integration method uniformly across all cell types. We adopt a tiered strategy: immune and other non-resident cells, which are well separated from disc cells in feature space, tolerate standard scVI integration, while resident cells are integrated under conservative settings that preserve within-compartment heterogeneity.

This manuscript reports the v4 tiered atlas — the second-generation analysis of our 12-dataset aggregation following an earlier CCA-based v5 pipeline. The two versions are methodologically distinct; results from both are presented in the supplementary materials of the project repository, and version-sensitive findings are flagged explicitly here so that downstream users can evaluate claims with appropriate uncertainty.

---

## Methods

### Datasets and harmonization

Twelve human IVD scRNA-seq datasets were identified via GEO and CNGB and downloaded with their original count matrices (Module 01). Twenty-three candidate studies were initially screened; eleven were excluded for being mouse/rat models, animal-disc tissue, or for unavailable count data. The final twelve studies span 78 samples and 57 donors, with conditions covering healthy adult, mild and severe degeneration, ungraded degenerative, herniated, neonatal, and aged-ungraded states. Per-study count matrices were converted to AnnData (.h5ad) and harmonized at the sample-metadata level (Module 02): condition labels were collapsed to seven canonical categories (`healthy`, `degenerated_mild`, `degenerated_severe`, `degenerated_ungraded`, `herniated`, `aged_ungraded`, `neonatal`), and compartment was recorded as NP, AF, CEP, or IVD_mixed for the one dataset (GSE189916) that did not separate compartments.

| Compartment | Samples |
|---|---:|
| NP | 49 |
| AF | 17 |
| CEP | 6 |
| IVD_mixed | 6 |

Datasets with published corrigenda (GSE205535) were processed against the corrected metadata. Three datasets use non-10x Genomics chemistries (BD Rhapsody, Singleron); platform identity was retained in the integration covariate set.

### Quality control and preprocessing

Per-cell QC (Module 03) used pooled ambient-RNA estimation per sample, scrublet-based doublet detection, and conservative cutoffs on `pct_counts_mt` (≤20%), `n_genes_by_counts` (≥200, ≤8000 depending on platform), and `total_counts` (≥500). Genes detected in fewer than three cells per dataset were dropped. Counts were stored as raw integers; normalization (CP10K + log1p) was applied lazily on a per-analysis basis to preserve raw counts for pseudobulk aggregation.

### Coarse cell classification

A panel-based scoring approach (Module 04) assigned every cell a coarse class — `mesenchymal`, `immune`, `endothelial`, or `unknown` — using marker panels for canonical disc-resident lineages (COL1A1, COL2A1, ACAN, PRG4), immune subsets (CD3D, CD8A, CD68, CD14, LYZ, S100A8, MS4A1, MZB1), red blood cells (HBB, HBA1, HBA2), and endothelial cells (PECAM1, CDH5, EMCN, VWF). Cells with no panel reaching threshold were labelled `unknown` and routed to the mesenchymal integration tier (these cells are predominantly low-quality disc cells rather than missed lineages).

### Tiered integration

The tiered integration strategy (Module 05) splits the data into two streams:

- **Non-mesenchymal tier (`05k`):** immune, endothelial, RBC, and other non-resident cells, integrated by scVI with default conservatism. This tier handles the well-separated cell lineages where standard batch correction is reliable.
- **Mesenchymal tier (`05m`):** disc-resident chondrocyte- and fibroblast-like cells plus the routed `unknown` class, integrated by scVI under settings designed to preserve within-compartment continuum structure (higher KL penalty, restricted batch covariate).

Both integrations operate per-compartment (NP, AF, CEP) and on the union (`all_cells`), producing four .h5ad objects with `X_integrated` 30-dimensional latent representations. The tiered approach was selected over flat single-method integration after a 2026-04-17 NP quality experiment (Module 05 supplementary) demonstrated that standard CCA flattened the NP fibrocartilaginous-to-mature-chondrocyte continuum into a single homogeneous cluster.

### Clustering

Leiden clustering (Module 06) was applied per tier per compartment, scanning resolutions chosen by an equal-weighted silhouette+modularity score. Tier-aware adaptive thresholds — three resolutions for >300K cells, six for >200K, ten for >50K — and skipping modularity computation for >100K cells kept run times tractable.

| Compartment | Mesenchymal clusters | Non-mesenchymal clusters | Total |
|---|---:|---:|---:|
| NP | 17 | 4 | 21 |
| AF | 9 | 3 | 12 |
| CEP | 4 | 4 | 8 |
| all_cells | 24 | 9 | 33 |

### Annotation

Annotation was performed in three stages (Module 07):

1. **Coarse:** confirmed/refined the Module 04 coarse class using post-integration neighbourhoods. CellTypist (normalized inputs) provided independent reference labels for the non-mesenchymal tier; a 60% per-cluster majority threshold against Module 04's `coarse_label` was used as a fallback when CellTypist scoring did not fire.
2. **Cell type (compartment-prefixed):** panel-based scoring against `FINE_PANELS` produced compartment-prefixed labels (NP_fibrocartilaginous, NP_mature_chondrocyte, NP_fibrochondrocyte_chondroid, AF_inner, AF_outer, CEP_outer, CEP_hyaline, CEP_fibrochondrocyte_fibroid, plus Macrophage_M1/M2, Neutrophil, Immune, Endothelial, Pericyte_SMC, Erythrocyte). A label-harmonization step (07b) on 2026-05-22 collapsed historical generic labels (e.g. `Fibroblast_like` → compartment-specific) for cross-compartment consistency.
3. **Sub-state:** overlap-based scoring against `SUBSTATE_PANELS` (`proliferating`, `inflammatory`, `stressed`, `matrix_active`, `migratory`, `homeostatic`) plus an endothelial-admixed contamination flag (CD34/EMCN/AQP1 panel). Mesenchymal cells received a `cell_subtype` label; non-mesenchymal cells inherited `cell_subtype = cell_type`. 27 sub-states were called in the all_cells object.

**Contamination handling:** 16,514 Erythrocyte cells (4.0% of the atlas) and 1,831 endothelial-admixed cells (0.4%) were retained with `is_contamination=True` and `contamination_type ∈ {RBC, endothelial_admixed, clean}` flags rather than removed. This decision (made at the 2026-05-22 annotation checkpoint) preserves the count for downstream analyses while allowing transparent filtering in interpretation.

### Pseudobulk differential expression

Per-sample pseudobulk count aggregation followed by DESeq2 (Module 08) was used for differential expression to avoid the inflated false-positive rate of single-cell DE tests that treat individual cells as independent observations. The primary grouping was `cell_type` (compartment-prefixed); composition tests were additionally run on `cell_subtype`. Comparisons were defined as `healthy_vs_degenerated_all`, `healthy_vs_degenerated_mild`, `healthy_vs_degenerated_severe`, and `mild_vs_severe` per cell type. Comparisons with fewer than 3 samples per group were marked underpowered and skipped. The Macrophage_M2 healthy-vs-severe contrast — which returned 5,659 DEGs from a 7-sample comparison and showed clear evidence of statistical inflation in the per-gene dispersion profile — was retained on disk for transparency but excluded from downstream interpretation (Modules 09, 11, and the manuscript).

### Pathway enrichment and TF activity

Over-representation analysis (ORA) against GO/KEGG/Reactome via Enrichr, pre-ranked GSEA (gseapy, MSigDB + a custom IVD gene-set library) and CollecTRI-based transcription factor activity inference via decoupler were run per significant cell-type × comparison contrast (Module 09). Pain-gene cross-referencing used a curated 60-gene panel covering inflammatory pain, neurotrophin signalling, nerve guidance, neovascularization, neuropeptides, and ion-channel nociception categories.

### Trajectory analysis

PAGA-initialized diffusion pseudotime (DPT) on the mesenchymal tier (Module 10) computed per-compartment trajectories. NP and AF were downsampled to 50,000 cells; CEP processed in full (36,879). The neighbor graph used `X_integrated`. Root clusters were chosen by maximum enrichment of the expected mature/inner population (NP_mature_chondrocyte for NP, AF_inner for AF, CEP_hyaline for CEP). RNA velocity was not feasible — public count matrices do not include spliced/unspliced layers.

### Cell-cell communication

LIANA (consensus ligand–receptor rank-aggregation across CellPhoneDB, NATMI, Connectome, logFC, sca, geometric mean) was run per condition group (healthy vs. degenerated pooled) on the union atlas, with 20,000-cell-per-condition caps for tractability (Module 11). Pain-relevant interactions were flagged by ligand-receptor membership in a curated pain panel.

### Software

scanpy 1.10, anndata 0.12.10, scvi-tools 1.4.2 (PyTorch 2.10.0+cu128 on CPU), Seurat 5.4.0, DESeq2 1.42.1, decoupler-py, gseapy, LIANA-py. Python 3.12 in a project virtualenv; R 4.4 with Bioconductor 3.20. Reproducibility: every script accepts `--input-dir` / `--output-dir` CLI flags, logs are retained per run, and a per-file checksum manifest is maintained at `metadata/file_checksums.json`.

---

## Results

### §1 — A 411,000-cell tiered atlas of the human intervertebral disc

The integrated atlas comprises 410,705 cells across 78 samples from 57 donors, distributed across NP (262,924 cells, 64%), AF (84,617 cells, 21%), and CEP (50,854 cells, 12%), with the remaining cells from GSE189916 carrying an `IVD_mixed` compartment label. After three-stage annotation, 19 cell types were called at the union level (Table 1), with finer compartment-prefixed types within each compartment (NP 10, AF 7, CEP 7). 27 sub-states were resolved.

The largest cell type by count is NP_fibrocartilaginous (94,597 cells, 23% of the atlas), followed by AF_outer (48,836), NP_fibrochondrocyte_chondroid (59,650), NP_mature_chondrocyte (33,010), Macrophage_M2 (22,752), Immune (22,404), Neutrophil (27,912), CEP_outer (15,557), and CEP_hyaline (12,306). Two contamination categories — Erythrocyte (16,514 cells; 4.0%) and `_endothelial_admixed` cells within NP_fibrocartilaginous (1,831 cells) — are retained with explicit `is_contamination=True` flags rather than filtered. Total flagged contamination is 4.5% of the atlas.

The tiered integration preserves within-compartment continuity: NP fibrocartilaginous, fibrochondrocyte-chondroid, and mature-chondrocyte populations are graded rather than discrete in the latent space, consistent with the morphological and histological literature describing IVD resident cells as a continuum.

### §2 — NP fibrocartilaginous cells dominate the degenerative transcriptional signature

Of the 27 statistically powered cell-type × comparison contrasts (Table 2; one inflated contrast excluded), 1,823 significant DEGs (FDR < 0.05) were detected. The cell-type distribution is highly asymmetric:

| Cell type | Healthy vs. mild | Healthy vs. severe | Mild vs. severe | Total |
|---|---:|---:|---:|---:|
| NP_fibrocartilaginous | 14 | 291 | **325** | 630 |
| NP_fibrochondrocyte_chondroid | 5 | 27 | 350 | 382 |
| NP_mature_chondrocyte | 3 | 121 | 174 | 298 |
| AF_outer | **121** | 3 | 1 | 125 |
| Immune | 4 | 38 | 11 | 53 |
| Neutrophil | 2 | 23 | 25 | 50 |
| Macrophage_M1 | — | — | 21 | 21 |
| Macrophage_M2 | 1 | (excluded) | 171 | 172 |
| AF_inner | 2 | — | — | 2 |

The three NP cell types collectively account for 1,310 of 1,823 (72%) of trustworthy DEGs, with NP_fibrocartilaginous showing the largest single mild-vs-severe contrast (325 DEGs: 189 up, 136 down). Top upregulated DEGs in severe-vs-mild NP_fibrocartilaginous include MMP3, MMP13, IL6, CXCL8, NTN1, PLA2G2A, PENK, and PTGS2; top downregulated include COL2A1, ACAN, COL9A2, COL11A1, and PRG4. The pattern is consistent across compartments and consistent with the canonical catabolic-anabolic switch.

### §3 — AF outer cells respond early, then quiesce

AF_outer cells exhibit a distinctive temporal pattern: 121 DEGs at healthy-vs-mild, dropping to 3 at healthy-vs-severe, with no detectable mild-vs-severe signal. The directionality is also notable — the mild-vs-healthy contrast is dominated by downregulation (110 of 121 genes), including ECM-related transcripts and the pain ligands CXCL8, NGFR, and PLA2G2A. This pattern suggests a transient transcriptional stress response in AF outer fibroblasts during early degeneration that has resolved (or been overwritten by cellular dropout) by the severe stage. Whether this reflects a population shift, a transcriptional reset, or sampling differences across donors at severe disease stages is not resolvable from cross-sectional human tissue alone.

AF_inner cells, in contrast, produced only 2 DEGs across all comparisons — consistent with either a quiescent inner phenotype or with statistical limitations from the smaller AF_inner population (23,769 cells across 17 AF samples).

### §4 — A pain-associated gene landscape with cell-type specificity

Eighteen unique pain-associated genes are differentially expressed at FDR < 0.05 across one or more contrasts (Table 3):

| Gene | Category | Direction in disc degeneration | Cell types with significant signal |
|---|---|---|---|
| IL6 | Inflammatory pain | UP (severe) | NP_fibrocartilaginous, Immune |
| IL1B | Inflammatory pain | DOWN (degen_all) | NP_fibrocartilaginous |
| CCL2 | Inflammatory pain | UP | NP_fibrocartilaginous, Neutrophil, Immune |
| TNF | Inflammatory pain | UP (mild→severe) | Neutrophil |
| CXCL8 | Inflammatory pain | DOWN (healthy→mild) | AF_outer |
| PTGS2 | Inflammatory pain | UP (mild→severe) | Macrophage_M1 |
| PLA2G2A | Inflammatory pain | UP / DOWN | NP_fibrocartilaginous (UP severe), Immune, AF_outer (DOWN mild) |
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

NP_fibrocartilaginous cells carry the broadest pain signature (6 of 18 unique genes including IL6, NTN1, UNC5B, PENK, PLA2G2A, CCL2). The neovascularization triad (VEGFA, PDGFA, FGF2) is concentrated in NP_mature_chondrocyte and NP_fibrochondrocyte_chondroid, consistent with nucleus-pulposus-driven angiogenic signalling. The neurotrophic and nerve-guidance signals (NTN1, UNC5B, NGFR) are predominantly cell-autonomous to NP fibrocartilaginous cells, suggesting NP cells themselves participate in establishing the nerve-ingrowth permissive environment described in chronic discogenic pain.

The opposing direction of AF_outer signals (CXCL8, NGFR, PLA2G2A all DOWN at healthy-vs-mild) underscores the §3 finding that the AF outer compartment undergoes a different early response than NP cells do.

### §5 — Pathways and transcription factors converge on ECM catabolism and NF-κB

Pathway enrichment (3,051 ORA + 6,890 GSEA significant terms at FDR < 0.05) and transcription factor activity (555 significant TF–contrast pairs from CollecTRI) converge on four themes:

1. **Extracellular matrix degradation upregulated in degeneration.** MMP1/2/3/9/13 and ADAMTS4/5 in NP fibrocartilaginous, mature-chondrocyte, and fibrochondrocyte-chondroid contrasts. GO terms `extracellular matrix organization`, `collagen catabolic process`, and `proteoglycan metabolic process` lead the upregulated NP_fibrocartilaginous severity comparisons.
2. **Inflammatory and NF-κB signalling broadly upregulated.** NFKBIA, NFKBIZ, IER3, PTGS2, IL6, CXCL8, SOD2 across NP cell types. CollecTRI TF activity scores rank NFKBIZ, ARID3A, HCFC1, and STAT3 among the most active regulators in degenerated NP states.
3. **Anabolic / repair programs downregulated.** COL2A1, ACAN, COL9A1/2/3, COL11A1/2, COMP, PRG4 all downregulated in NP_fibrocartilaginous mild-vs-severe and healthy-vs-severe. Canonical chondrocyte differentiation pathways (SOX9 targets, TGFB-driven chondrogenesis) are negatively enriched.
4. **Pain-related cell-type-specific dysregulation.** ORA over the curated IVD gene-set library confirms enrichment of `inflammatory_pain`, `neurotrophin_signalling`, and `neovascularization` categories in NP_fibrocartilaginous and NP_mature_chondrocyte degeneration contrasts.

The tiered_v4 pipeline yields 22–109% more significant signals than the prior v5 (CCA-flat) analysis across all four categories (ORA, GSEA, TFs, pain genes; see methods supplementary). Mechanistically, this reflects (i) the finer cell-type partition exposing per-cell-type signal that v5 averaged across; (ii) the rescued Immune lymphocyte cluster (17K NP + 5K AF cells) contributing previously-invisible signal; and (iii) a larger trustworthy DE pool (1,823 vs. 1,198 sig DEGs). The biological themes are preserved across versions.

### §6 — Cell-cell communication remodeling in degeneration

LIANA consensus rank-aggregation on 20,000 cells per condition (downsampled for tractability) identified **66,827 ligand-receptor interactions in healthy tissue and 76,019 in degenerated tissue** (90,650-interaction union across conditions). **6,883 interactions involve a pain-panel ligand or receptor**. Pain-relevant interactions in the degenerated pool are dominated by neovascularization (VEGFA: 1,485 interactions, FGF2: 1,157, VEGFC: 246, VEGFB: 216), neurotrophic (NGF: 360), semaphorin nerve-guidance (SEMA3C: 300, SEMA3A: 250), inflammatory (TNFSF10: 250, PTGS2: 224), and granulin (GRN: 285) ligand contributions.

Differential analysis (rank_diff = healthy_rank − degenerated_rank, where a positive value indicates an interaction gained in degeneration) yields two strong patterns:

**Gained in degeneration.** Pan-cell-type → Neutrophil **FN1 → C5AR1** and **RPS19 → C5AR1** signaling dominates the top of the gained list — fibronectin and ribosomal protein S19 acting as C5a-receptor ligands across virtually every disc-resident cell type. This is the canonical complement-driven neutrophil chemotaxis axis being broadly activated.

**Lost in degeneration.** A coordinate set of CEP_outer-centric signaling axes attenuates: **WNT2B → FZD4/LRP5/LRP6** (CEP_outer autocrine and to NP_mature_chondrocyte), **LAMB3 → ITGA6/ITGAV/ITGA2 integrins** (CEP-to-AF basement membrane), **NDP → FZD4** (Norrin/WNT), and **MDK → SDC3** (midkine). In parallel, **AF_outer / NP_fibrocartilaginous → CD4 HLA-class-II presentation** is also reduced. The lost-in-degeneration signal is therefore a coordinate loss of WNT, basement-membrane, and immune-presentation signaling rather than a single-axis effect.

Aggregate interaction counts are sensitive to subsampling and the cell-type partition (the tiered_v4 atlas resolves 19 cell types vs. the prior v5 atlas's 16, which alone accounts for much of the ~2.5–3× increase in raw interaction counts versus v5). The cell-type-pair-specific and ligand-specific calls above are the more interpretable level of the analysis.

### §7 — Cellular composition does not shift significantly with degeneration

Compositional analysis across 19 cell-type × comparison contrasts at the compartment level identified no significant changes after FDR correction (lowest adjusted p-value = 0.61). The largest absolute shifts — AF_inner enrichment in severe AF (log2FC +3.05, raw p = 0.04; +1.33 at healthy-vs-severe, raw p = 0.07), Macrophage_M2 depletion in NP severe degeneration (log2FC −4.28, raw p = 0.02), and Pericyte_SMC depletion in severe (log2FC −7.98, raw p = 0.07) — do not survive multiple-testing correction. The finer `cell_subtype` composition analysis (also FDR-corrected) similarly returns no significant shifts.

The conservative interpretation, consistent across pipeline versions, is that **disc degeneration is primarily a transcriptional, not a cellular-composition, phenomenon** at the cell-type level. Composition shifts at the sub-state level (proliferating, stressed, inflammatory, matrix_active fractions within each cell type) are more variable across donors and are best interpreted at the individual sub-state level rather than as compartment-wide claims.

### §8 — Cell-state trajectories: cell-type-specific, not compartment-wide

PAGA-initialized DPT pseudotime in each compartment yields the expected ordering at the cell-type level — root clusters enriched for the mature/inner population (NP_mature_chondrocyte, AF_inner, CEP_hyaline) progress along a continuous axis toward the more degenerative/outer populations. However, the **compartment-level pseudotime-vs-condition correlation is near zero for NP and AF** in tiered_v4 (NP ρ = −0.004, p = 0.41; AF ρ = −0.003, p = 0.49), while CEP retains a weak positive correlation (ρ = +0.077, p = 1.9e-49). At the per-cell-type level, gradients are stronger and **opposing** within each compartment: AF_inner ρ = +0.54 vs. AF_outer ρ = −0.30; CEP_hyaline ρ = +0.17 vs. CEP_outer ρ = −0.21.

This pattern is consistent with: (i) DPT capturing within-cell-type stress/maturation gradients rather than a single compartment-wide degenerative axis, and (ii) the finer cell-type partition redistributing the compartment-level signal that the prior v5 analysis aggregated. Mann-Whitney U tests of healthy-vs-degenerated pseudotime medians remain highly significant in all three compartments (NP p = 3e-4, AF p = 1e-20, CEP p = 2e-43), so the discrete healthy/degenerated shift is robust even though the compartment-level direction is not.

Trajectory-associated genes overlap substantially with the DE results: 328/500 (NP), 326/500 (AF), 338/500 (CEP) trajectory genes are also DE genes — 65–68% per compartment.

---

## Discussion

### A coherent transcriptional model of disc degeneration

The atlas supports a model in which disc degeneration is principally a cell-type-specific transcriptional remodeling rather than a compositional collapse. NP fibrocartilaginous cells drive the catabolic signature, upregulating matrix-degrading enzymes (MMP, ADAMTS), inflammatory mediators (IL6, CXCL8, PTGS2), nerve-guidance cues (NTN1, UNC5B), neuropeptides (PENK), and neurotrophin receptors, while suppressing anabolic ECM transcripts (COL2A1, ACAN, COL9A1, PRG4). NP mature chondrocytes contribute the neovascularization arm (VEGFA, PDGFA, FGF2). NP fibrochondrocyte-chondroid cells participate in the inflammatory arm (BDKRB1, CCL2) with their own distinctive pain-mediator signature. The result is a coherent biology in which the NP cellular continuum acts collectively as both effector and amplifier of the catabolic / inflammatory / nociceptive cascade.

AF outer cells contribute a temporally-distinct early response that quiesces by the severe disease stage. AF inner cells, CEP cells, and the immune compartment all participate in narrower but coherent ways — Macrophage_M1 polarization in late degeneration is one example (PTGS2 and FLT1 upregulation at mild-vs-severe), and the Immune compartment contributes IL6, NGFR, CCL2, and PLA2G2A upregulation at the severe stage.

### Methodological considerations

The tiered_v4 pipeline differs from prior efforts in three substantive ways. First, tiered integration: by integrating immune and resident cells separately, we preserve the resident-cell continuum that uniform CCA or scVI integration tends to compress. Second, three-stage annotation: the addition of compartment-prefixed cell types and overlap-scored sub-states yields a more interpretable label set than majority-class assignment alone. Third, contamination flagging rather than filtering: 4.5% of the atlas is flagged but retained, allowing transparent inspection rather than silent removal of cells whose biology may inform the contamination calls themselves.

The pseudobulk DESeq2 approach is used in preference to single-cell DE methods (Wilcoxon, MAST) to avoid the well-documented inflation of false positives when treating cells as independent observations. The one inflated contrast (Macrophage_M2 healthy-vs-severe, 5,659 DEGs from 7 samples) is excluded from interpretation specifically because its DE count is incompatible with the donor-level statistical model — likely driven by a sample-imbalance × dispersion-estimation interaction. The contrast remains on disk for transparency.

### Comparison with v5 (CCA-flat)

The earlier v5 CCA-flat analysis and the current tiered_v4 analysis yield substantially different DE counts, pathway enrichment magnitudes, and trajectory correlations. Where the two pipelines disagree on direction or magnitude, **we treat this as informative rather than as a failure of either version**. The mechanism is consistent and reproducible across the comparisons we have run: finer cell-type partitioning in tiered_v4 redistributes signal that v5 averaged across coarser groups, increasing per-cell-type signal at the cost of compartment-level summary metrics.

Findings that are robust across both versions — ECM catabolism upregulated, NF-κB-driven inflammation increased, anabolic/repair pathways downregulated, no compositional shifts after FDR, disc cells as inflammatory mediators — are the strongest claims of this atlas. Findings that change direction or magnitude between versions — specific pseudotime directionalities, individual CCC interaction counts, particular TF significance calls — are presented with version-sensitivity caveats and should not be cited as strong claims without independent replication.

### Therapeutic implications

The cell-type-specific resolution of the pain-gene panel suggests that targets distinguish themselves by where they act: NP fibrocartilaginous cells are a candidate target for neurotrophin- and nerve-guidance-axis interventions (NTN1, NGFR, PENK); NP mature chondrocytes and fibrochondrocyte-chondroid cells are the angiogenic-mediator source (VEGFA, PDGFA, FGF2); Macrophage_M1 cells are an inflammatory amplifier (PTGS2, FLT1). PTGS2 (COX-2) recurs across NP fibrocartilaginous, Macrophage_M1, and pathway-level analyses, supporting its standing as a high-priority therapeutic node. None of these claims is novel in isolation — each has prior literature support — but the cell-type-specific source of each signal provides a more refined picture of where intervention might target.

---

## Caveats and limitations

The following caveats apply to specific claims in this manuscript and are flagged explicitly rather than hidden in the supplementary materials:

1. **Macrophage_M2 healthy-vs-severe contrast is statistically inflated.** Its 5,659-DEG count from a 7-sample comparison is excluded from all downstream analyses (Module 09, 11, and this manuscript). Per-gene dispersion estimates flag clear pathology; cite contrasts with a more balanced donor distribution instead.

2. **Trajectory pseudotime-condition correlations are version-sensitive.** The compartment-level direction (NP, AF) changes between the v5 and tiered_v4 pipelines (NP ρ ∈ {−0.088, −0.004}; AF ρ ∈ {+0.195, −0.003}). The discrete healthy-vs-degenerated MWU shift is robust; the compartment-level direction is not a strong claim.

3. **Contamination cells retained with flags, not filtered.** 16,514 Erythrocyte cells (predominantly in NP and CEP) and 1,831 endothelial-admixed NP_fibrocartilaginous cells carry `is_contamination=True` flags. DE analyses on these cell types reflect contamination-tracking signal, not biology, and are read accordingly.

4. **GSE189916 (IVD_mixed) cells carry generic cell-type labels.** That dataset does not separate compartments at the tissue level, so its cells retain non-compartment-prefixed labels (`Fibroblast_like`, `Chondrocyte_like`, etc.) rather than being forced into compartment-specific bins.

5. **GSE242443 CEP cells are culture-expanded.** They are included by prior decision but contribute a culture-derived rather than fresh-tissue signal; CEP results should be interpreted with this in mind.

6. **GSE230809 donors are all male.** This dataset's age and disease effects are sex-confounded; cross-condition CEP claims that lean on it carry that confound.

7. **GSE205535 has published corrigenda.** Corrected metadata was used; the published count matrices were not affected.

8. **No RNA velocity.** Public count matrices do not include spliced/unspliced layers. Pseudotime is DPT-only.

9. **CEP is underpowered for DE.** Only 6 CEP-containing samples and 50,854 cells; no CEP DE comparisons reach the powered threshold for healthy-vs-mild contrasts. CEP claims rest primarily on annotation, composition, and trajectory analyses.

10. **Compositional shifts are not FDR-significant.** Several raw-p < 0.05 shifts (AF_inner enrichment in severe AF, Macrophage_M2 depletion in severe NP, Pericyte_SMC depletion in severe disc) are reported descriptively but do not survive multiple-testing correction.

11. **LIANA CCC results are sensitive to subsampling and database choice.** Cell-cell communication counts at the aggregate level have been version-variable across pipeline iterations. Report cell-type-pair-specific interaction calls rather than aggregate-count differences.

---

## Data and code availability

All 12 source datasets are publicly available at GEO or CNGB (accession numbers in Table 1 / project repository). The full analysis pipeline (12 modules, scripts, specifications, notebooks, and metadata) is at https://github.com/andrewsu/lotz-ivd. Per-module logs, intermediate AnnData files, and result tables are reproducible end-to-end from the raw count matrices via `python3 scripts/0X_*.py` with the documented CLI flags. Per-file checksums are recorded in `metadata/file_checksums.json`.

---

## Author contributions, acknowledgments, references

*To be populated at submission.*
