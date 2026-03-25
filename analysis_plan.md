# IVD Analysis Plan

## Current Status

**Pipeline v5 — COMPLETE.** CCA integration, all 12 modules finished.

**Pipeline version:** v5

---

## Pipeline Summary (v5)

| Module | Status | Script | Description |
|--------|--------|--------|-------------|
| 01: Dataset Discovery | Complete (v1) | 01_dataset_download.py | 12 datasets downloaded |
| 02: Metadata Harmonization | Complete (v1) | 02_metadata_harmonization.py | Condition mappings finalized |
| 03: Preprocessing | Complete (v1) | 03_preprocessing.py | 12 per-dataset h5ad files, ~429K cells |
| 04: Coarse Classification | Complete (v4) | 04_annotation.py | 5 coarse categories + Unknown |
| 05: Integration | **COMPLETE** | 05a_integration_cca.R | CCA selected (label-free, full-cell). scANVI+STACAS for comparison. |
| 06: Clustering | **COMPLETE** | 06_clustering.py | NP 12, AF 12, CEP 9, all_cells 15 clusters |
| 07: Post-Integration Annotation | **COMPLETE** | 07_annotation.py | NP 5, AF 4, CEP 7, all 16 cell types |
| 08: Differential Analysis | **COMPLETE** | 08_differential.py | 17 powered comparisons, 1,198 sig genes |
| 09: Biological Interpretation | **COMPLETE** | 09_interpretation.py | 2,506 enrichments, 288 TFs, 10 pain genes |
| 10: Trajectory Analysis | **COMPLETE** | 10_trajectory.py | NP rho=-0.088, AF +0.195, CEP +0.073 |
| 11: Cell-Cell Communication | **COMPLETE** | 11_communication.py | 25K healthy vs 34K degenerated interactions |
| 12: Final Reporting | **COMPLETE** | 12_reporting.py | 19 supplementary tables, final report |

## Active Step

**Pipeline v5 COMPLETE (2026-03-25).** All 12 modules finished with CCA integration.

### Module 05 Workflow Selection (2026-03-25)

**Decision: CCA (Seurat v5) selected as primary integration workflow.**

Three workflows compared with full integration metrics (iLISI, batch_ASW, condition_ASW):

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

*\*STACAS downsampled for NP/all_cells (RAM-bound)*

**Rationale for CCA:**
- Label-free: no circular annotation risk (does not depend on Module 04 coarse labels)
- Full cell counts for all 4 objects (no downsampling)
- Strongest batch mixing (iLISI 1.5-3.7 vs ~1.0-1.2 for scANVI)
- Smooth embedding topology consistent with mesenchymal continuum hypothesis
- Negative batch_ASW indicates possible overcorrection, but DE uses pseudobulk on raw counts (not embeddings)

Now converting CCA RDS → h5ad and running Modules 06-12.

> CCA operational incident log: see [`docs/version_history.md`](docs/version_history.md#cca-run-incident-log-v5-2026-03-24).

---

## v5 Module 05 Integration Results

### Workflow A: CCA (Seurat v5, label-free)

Uses `IntegrateLayers(method = CCAIntegration)` with NormalizeData (log-normalization keeps layers split). Standard CCA for all objects — no downsampling (247GB RAM machine).

**Re-running full-cell CCA** (previous results used downsampled NP 15K, all_cells 11K on 62GB machine).

| Object | Cells | Status | Prior (downsampled) |
|--------|-------|--------|---------------------|
| NP | 262,967 | **Running** — at IntegrateLayers step | Was 15,000 (9 clusters) |
| AF | 84,624 | Pending | Was 84,624 (23 clusters, no change expected) |
| CEP | 50,858 | Pending | Was 50,858 (14 clusters, no change expected) |
| all_cells | 410,759 | Pending | Was 11,000 (8 clusters) |

### Workflow B: scANVI (semi-supervised, full cell counts)

Tiered scANVI with 5 coarse anchor categories. All cells integrated (no downsampling).

| Object | Cells | Clusters (res=0.5) | Time |
|--------|-------|-------------------|------|
| NP | 262,967 | 29 | ~55 min |
| AF | 84,568 | 18 | ~15 min |
| CEP | 50,769 | 13 | ~10 min |
| all_cells | 410,759 | 29 | ~90 min |

### Workflow C: STACAS (semi-supervised, R-native)

Uses `Run.STACAS()` with coarse_label anchors. Objects >100K cells downsampled.

| Object | Cells | Clusters (res=0.5) | Method |
|--------|-------|-------------------|--------|
| NP | 16,000 | 21 | Downsampled (2,000/study) |
| AF | 84,624 | 23 | Standard |
| CEP | 50,858 | 15 | Standard |
| all_cells | 30,000 | — | Downsampled (2,000/study) |

### Key observations

- **CCA re-running at full cell counts** on 247GB machine (previously downsampled on 62GB). This will enable a fair apples-to-apples comparison with scANVI across all objects.
- **STACAS still downsampled** for NP and all_cells (RAM-bound in R even on 247GB — STACAS memory footprint is higher than CCA).
- **Cluster counts are consistent across workflows** for AF (23, 18-23) and CEP (13-15), suggesting stable structure.
- **NP cluster counts diverged previously:** CCA 9 (downsampled 15K) vs scANVI 29 (full 263K). Full-cell CCA will clarify whether this was a downsampling artifact or a real methodological difference.

---

## Condition Mapping Decisions (reviewed 2026-03-05)

1. **Herniated samples:** Kept as separate category. GSE233666 excluded in v2+ (herniated-only study confounds comparisons). GSE251686 herniated samples retained, treated as "severe."
2. **GSE205535 NNP (11yo trauma):** Excluded from DE comparisons (not representative of healthy disc biology).
3. **Thompson III boundary:** II-III → degenerated_mild, III-IV → degenerated_severe. Conservative.
4. **Neonatal (GSE189916, n=3):** Separate category, not mixed into healthy.
5. **Aged ungraded (GSE189916 adult, n=3):** "aged_ungraded" — excluded from healthy vs. degenerated.
6. **Degenerated ungraded (GSE205535_DNP + GSE255768, n=3):** Included in "degenerated_all" but not mild/severe.

**DE comparison plan:**
- Primary: healthy vs. degenerated_all, healthy vs. degenerated_severe, healthy vs. degenerated_mild
- Secondary: mild vs. severe
- Per cell type per compartment where sample counts ≥ 3 per group

---

## Integration Approach (v5 — three parallel workflows)

**Three integration workflows** run in parallel on each compartment object (NP, AF, CEP, all_cells):

- **Workflow A (Seurat CCA, v5):** R-only, label-free. Seurat v5 `IntegrateLayers(method=CCAIntegration)` with `NormalizeData`. Standard CCA for all objects (no downsampling on 247GB machine).
- **Workflow B (scANVI):** Python, semi-supervised with coarse anchor labels from Module 04. Tiered (mesenchymal + non-mesenchymal). Full cell counts via GPU.
- **Workflow C (STACAS):** R-only, `Run.STACAS()` with coarse_label anchors. Large objects downsampled.

scANVI and STACAS complete. CCA re-running at full cell counts on 247GB machine (previously downsampled). Human checkpoint deferred until CCA finishes. See `specs/05_INTEGRATION.md` for details, `docs/v5_execution_dialogue.md` for execution history.

---

## Items Requiring SME Review (v5)

1. **Trajectory results (v5):** NP rho=-0.088, AF rho=+0.195, CEP rho=+0.073. Weak correlations suggest pseudotime does not strongly track degeneration severity in the CCA embedding. This result is sensitive to integration method and root cell choice — correlations changed sign across prior versions.

2. **CellTypist concordance (v5):** CellTypist lacks IVD-specific reference types, so disagreements with de novo mesenchymal labels are expected. De novo labels retained as primary; CellTypist used for immune subtype validation only.

3. **CCC direction (v5):** 25,537 healthy vs 34,208 degenerated interactions (more in degeneration). This direction has varied across pipeline versions and should be treated as uncertain.

4. **AF pseudotime sign (v5):** AF rho=+0.195 (degenerated at later pseudotime), opposite to NP (-0.088). Consistent across all pipeline versions. May reflect AF-specific biology or root cell choice effects.

> Cross-version sensitivity analysis with detailed version-by-version comparisons: see [`docs/version_history.md`](docs/version_history.md#cross-version-sensitivity-observations).

---

## Deferred Questions

- Should spatial transcriptomics data be incorporated? No human IVD spatial datasets found.
- Should the final atlas be deposited to CellxGene? No IVD data currently exists there — this would be the first.
- Cross-species validation with mouse/rat/bovine data? Datasets identified but not incorporated.

---

## Known Issues

- **NGDC datasets excluded:** PRJCA014236 and PRJCA007656 not downloaded. NP already well-covered.
- **GSE205535 corrigenda:** Published corrections exist — reviewed during preprocessing.
- **Platform heterogeneity:** 3 non-10x datasets (BD Rhapsody, Singleron). Handled by CCA batch correction (v5 primary). scANVI and STACAS also tested.
- **SeuratDisk incompatible with Seurat v5 (implementation):** Workaround in place — R export to MTX/CSV + Python assembly (`scripts/seurat_to_h5ad_bridge.R` + `scripts/seurat_to_h5ad_assemble.py`).
- **CEP underpowered:** Only 3 CEP datasets (6 samples). Compartment-specific CEP analyses are limited.
- **GSE242443 culture-expanded:** CEP cells are culture-expanded. Included with caveats.
- **GSE230809 sex bias:** All 24 samples from male donors. Limits sex-stratified analyses.
- **GSE230809 age-disease confound:** Healthy=21-27y, diseased=37-73y. Cannot separate age from disease.
- **Missing demographics:** 18/78 samples unknown age, 30/78 unknown sex.
- **GSE251686_NP3 excluded:** Corrupt matrix file (5/6 samples retained).
- **GSE165722 Pfirrmann offset:** GEO says I-IV, paper says II-V. Paper grades used.

---

## Version History

> Full changelog (v1–v5), key decisions log, CCA incident log, and cross-version
> sensitivity analysis: see [`docs/version_history.md`](docs/version_history.md).
