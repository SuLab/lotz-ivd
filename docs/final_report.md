# Human Intervertebral Disc Single-Cell Atlas — Final Report

**A comprehensive scRNA-seq meta-analysis of IVD degeneration**

| Field | Value |
|-------|-------|
| Report generated | 2026-03-25 13:45 |
| Pipeline version | v5 |
| Git commit | `49e6992` (branch: `main`) |
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

Cell type definitions not found. *[source: `results/integration/cell_type_definitions.tsv`]*

### Clustering resolution optimization

*[source: `results/integration/clustering_resolution_optimization`]*

- **TEST_mes:** best resolution 0.4, 2 clusters (silhouette=0.102)

## 5. Differential Expression {#de}

DE summary not found. *[source: `results/differential/de_summary_table.tsv`]*

## 6. Biological Pathways {#pathways}

ORA results not found. *[source: `results/interpretation/pathway_enrichment/all_enrichment_results.tsv`]*

GSEA results not found. *[source: `results/interpretation/pathway_enrichment/gsea_results.tsv`]*

## 7. Transcription Factor Activity {#tf}

TF activity results not found. *[source: `results/interpretation/tf_activity/tf_activity_results.tsv`]*

## 8. Cell State Trajectories {#trajectory}

Trajectory correlation results not found. *[source: `results/trajectories`]*

## 9. Cell-Cell Communication {#communication}

CCC interaction files not found. *[source: `results/communication`]*

## 10. Pain Biology {#pain}

Pain gene results not found. *[source: `results/interpretation/pain_genes.tsv`]*

## 11. Limitations & Caveats {#limitations}

From `analysis_plan.md` Known Issues section: *[source: `analysis_plan.md`]*

- **NGDC datasets excluded:** PRJCA014236 and PRJCA007656 not downloaded. NP already well-covered.
- **GSE205535 corrigenda:** Published corrections exist — reviewed during preprocessing.
- **Platform heterogeneity:** 3 non-10x datasets (BD Rhapsody, Singleron). Handled by scANVI batch correction. CCA and STACAS also correct for this via study-level integration.
- **SeuratDisk incompatible with Seurat v5:** `GetAssayData(slot=...)` removed in SeuratObject 5.0. Workaround: R export to MTX/CSV + Python assembly (`scripts/seurat_to_h5ad_bridge.R` + `scripts/seurat_to_h5ad_assemble.py`).
- **CEP underpowered:** Only 3 CEP datasets (6 samples). Compartment-specific CEP analyses are limited.
- **GSE242443 culture-expanded:** CEP cells are culture-expanded. Included with caveats.
- **GSE230809 sex bias:** All 24 samples from male donors. Limits sex-stratified analyses.
- **GSE230809 age-disease confound:** Healthy=21-27y, diseased=37-73y. Cannot separate age from disease.
- **Missing demographics:** 18/78 samples unknown age, 30/78 unknown sex.
- **GSE251686_NP3 excluded:** Corrupt matrix file (5/6 samples retained).
- **GSE165722 Pfirrmann offset:** GEO says I-IV, paper says II-V. Paper grades used.

### Items requiring SME review

*[source: `analysis_plan.md`]*

1. **Trajectory instability across versions:** Pseudotime-condition correlations change sign between pipeline versions (e.g., CEP went from -0.163 in v2 to +0.135 in v3). This sensitivity to upstream annotation choices means trajectory results should be interpreted cautiously.
2. **CellTypist NP disagreements:** 8/13 de novo NP clusters were discordant with CellTypist in v3. CellTypist lacks IVD-specific cell types, so de novo labels are retained, but this should be acknowledged.
3. **CCC direction sensitivity:** v1 showed more interactions in degeneration (53K vs 44K), v2 showed fewer (27K vs 29K), v3 shows near-equal (40K vs 41K). The direction of this result is sensitive to annotation and sampling choices.
4. **AF pseudotime sign:** AF consistently shows positive rho (degenerated at later pseudotime) across v2 and v3, opposite to NP. May reflect genuine AF-specific biology or root cell choice issues.

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

- **Git commit:** `49e6992614fbe62644d1978b9fa3642e27df00d5` (branch: `main`)
- **Random seeds:** 42 (all stochastic operations)
- **Package versions:** pinned in `requirements.txt`, frozen in `requirements_frozen.txt`
- **Parameter choices:** documented in `analysis_plan.md`
- **Human checkpoint decisions:** recorded in `analysis_plan.md`
- **Data provenance:** GEO/CNGB accessions and download dates in `metadata/dataset_registry.tsv`
- **File checksums:** `metadata/file_checksums.json`
