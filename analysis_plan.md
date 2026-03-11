# IVD Analysis Plan

## Current Status

**Pipeline v4 ready.** Scripts restructured from 10-module to 12-module pipeline. Awaiting execution of Module 04+.

**Pipeline version:** v4 (spec restructuring + scANVI integration). See `results/THREE_VERSION_SUMMARY.md` for v1-v3 history.

## Active Step

**Module 04: Coarse Cell Classification — READY TO RUN**

Scripts 04-12 have been updated to match the restructured specs. Modules 01-03 outputs are unchanged from v3.

---

## Pipeline Summary (v4 — planned)

| Module | Status | Script | Description |
|--------|--------|--------|-------------|
| 01: Dataset Discovery | Complete (v1) | 01_dataset_download.py | 11 datasets, ~410K cells, 71 samples |
| 02: Metadata Harmonization | Complete (v1) | 02_metadata_harmonization.py | Condition mappings finalized |
| 03: Preprocessing | Complete (v1) | 03_preprocessing.py | 410,759 cells post-QC across 11 datasets |
| 04: Coarse Classification | Ready | 04_annotation.py | 5 anchor categories + Unknown for scANVI |
| 05: Integration | Ready | 05_integration.py | Tiered scANVI (semi-supervised) integration |
| 06: Clustering | Ready | 06_clustering.py | Leiden with resolution optimization |
| 07: Post-Integration Annotation | Ready | 07_annotation.py | Two-stage: coarse markers → fine DE |
| 08: Differential Analysis | Ready | 08_differential.py | Composition + pseudobulk DE |
| 09: Biological Interpretation | Ready | 09_interpretation.py | Pathways, GRN, pain genes |
| 10: Trajectory Analysis | Ready | 10_trajectory.py | PAGA + DPT pseudotime |
| 11: Cell-Cell Communication | Ready | 11_communication.py | LIANA ligand-receptor |
| 12: Final Reporting | Ready | 12_reporting.py | Comprehensive report + supplements |

### v4 Key Changes from v3

1. **Module 04 restructured:** Binary mesenchymal/non-mesenchymal → 5 coarse anchor categories (Chondrocyte_like, Fibroblast_like, Immune, Endothelial, Pericyte_SMC) + Unknown. Provides richer anchor labels for scANVI.

2. **scANVI replaces scVI:** Semi-supervised integration uses coarse_label anchors from Module 04. Should produce better batch correction, especially across platforms (10x, BD Rhapsody, Singleron).

3. **Module 05 split into 3:** Old monolithic integration+clustering+annotation script split into:
   - 05: Integration only (scANVI)
   - 06: Clustering with resolution optimization (NEW)
   - 07: Two-stage post-integration annotation (NEW)

4. **Downstream modules renumbered:** 06→08, 07→09, 08→10, 09→11, 10→12

5. **Two-stage annotation (Module 07):** Stage 1 assigns coarse identity via canonical markers. Stage 2 refines within coarse groups using cluster DE markers. More principled than v3's single-pass approach.

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

## Integration Approach (v4)

**Tiered scANVI** (semi-supervised). One scANVI model per tier per compartment object (NP, AF, CEP, all_cells) with `batch_key='study'` and `labels_key='coarse_label'`.

- Mesenchymal tier anchors: Chondrocyte_like, Fibroblast_like; Unknown cells are unlabeled (scANVI positions them by similarity)
- Non-mesenchymal tier anchors: Immune, Endothelial, Pericyte_SMC

Workflow: train scVI (max_epochs=200) → initialize scANVI from scVI → train scANVI (max_epochs=50, early_stopping).

---

## Items Requiring SME Review

1. **Trajectory instability across versions:** Pseudotime-condition correlations change sign between pipeline versions (e.g., CEP went from -0.163 in v2 to +0.135 in v3). This sensitivity to upstream annotation choices means trajectory results should be interpreted cautiously.

2. **CellTypist NP disagreements:** 8/13 de novo NP clusters were discordant with CellTypist in v3. CellTypist lacks IVD-specific cell types, so de novo labels are retained, but this should be acknowledged.

3. **CCC direction sensitivity:** v1 showed more interactions in degeneration (53K vs 44K), v2 showed fewer (27K vs 29K), v3 shows near-equal (40K vs 41K). The direction of this result is sensitive to annotation and sampling choices.

4. **AF pseudotime sign:** AF consistently shows positive rho (degenerated at later pseudotime) across v2 and v3, opposite to NP. May reflect genuine AF-specific biology or root cell choice issues.

---

## Deferred Questions

- Should spatial transcriptomics data be incorporated? No human IVD spatial datasets found.
- Should the final atlas be deposited to CellxGene? No IVD data currently exists there — this would be the first.
- Cross-species validation with mouse/rat/bovine data? Datasets identified but not incorporated.

---

## Known Issues

- **NGDC datasets excluded:** PRJCA014236 and PRJCA007656 not downloaded. NP already well-covered.
- **GSE205535 corrigenda:** Published corrections exist — reviewed during preprocessing.
- **Platform heterogeneity:** 3 non-10x datasets (BD Rhapsody, Singleron). Handled by scANVI batch correction.
- **CEP underpowered:** Only 3 CEP datasets (6 samples). Compartment-specific CEP analyses are limited.
- **GSE242443 culture-expanded:** CEP cells are culture-expanded. Included with caveats.
- **GSE230809 sex bias:** All 24 samples from male donors. Limits sex-stratified analyses.
- **GSE230809 age-disease confound:** Healthy=21-27y, diseased=37-73y. Cannot separate age from disease.
- **Missing demographics:** 18/78 samples unknown age, 30/78 unknown sex.
- **GSE251686_NP3 excluded:** Corrupt matrix file (5/6 samples retained).
- **GSE165722 Pfirrmann offset:** GEO says I-IV, paper says II-V. Paper grades used.

---

## Version History

### v4 (2026-03-11): Spec restructuring + scANVI
- Pipeline restructured from 10 to 12 modules (clustering and annotation split from integration)
- Module 04: 5 coarse anchor categories replace binary classification
- Module 05: scANVI (semi-supervised) replaces scVI (unsupervised)
- Module 07: Two-stage post-integration annotation (coarse → fine)
- Modules 01-03 unchanged
- Scripts updated, ready for execution

### v3 (2026-03-10): Annotation fix
- Three fixes to Module 04 classification (evidence gate, ACAN/SOX9 rescue, 85% voting)
- Recovered 17K stressed NP cells from non-mesenchymal to mesenchymal
- Modules 04-10 rerun; Modules 01-03 unchanged
- All reports, notebooks, supplementary tables regenerated

### v2 (2026-03-09 to 2026-03-10): Pipeline restructure
- Module 04 narrowed to binary classification; Module 05 expanded to include annotation
- Four compartment objects (NP, AF, CEP, all_cells) replace two-tier structure
- scVI-only replaces 4-approach integration benchmark
- GSE233666 excluded (herniated-only confound)
- Results archived in `results_v2/` and `results/v2_archive/`

### v1 (2026-02-26 to 2026-03-05): Original pipeline
- 12 datasets (including GSE233666), two-tier integration (resident/non-resident)
- 4-approach integration benchmark (scVI, scANVI, Harmony, BBKNN)
- scANVI chosen as primary, scVI for trajectory sensitivity
- 17 powered DE comparisons, 5,328 significant genes
- Results superseded by v2 restructuring

### Key decisions log
| Date | Decision |
|------|----------|
| 2026-02-26 | Specs approved. GSE242443 included. Zhou 2023 deferred. NGDC dropped. |
| 2026-02-26 | Condition mappings tentatively approved. Revisit before Module 08. |
| 2026-03-03 | Retroactive checkpoint review of Modules 03-05. No blocking issues. |
| 2026-03-05 | scANVI primary for v1 (later superseded). Condition mappings finalized. |
| 2026-03-05 | v1 pipeline complete. All modules 01-10 done. |
| 2026-03-09 | Spec restructuring. GSE233666 excluded. scVI-only. v2 rerun initiated. |
| 2026-03-10 | v2 complete. v3 annotation fix applied. Full rerun Modules 04-10. |
| 2026-03-10 | v3 complete. Stale files cleaned. Notebooks re-executed. Reports updated. |
| 2026-03-11 | Spec restructuring: 10→12 modules. scANVI integration. Scripts updated for v4. |
