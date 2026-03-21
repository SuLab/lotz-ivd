# IVD Analysis Plan

## Current Status

**Pipeline v5 — STARTING.** Fresh run. v4 results archived in `results_v4/`. Processed per-dataset files retained in `data/processed/`; integrated data cleared.

**Pipeline version:** v5

## Active Step

**Module 01: Dataset Discovery.** Verify existing downloads and metadata.

---

## Pipeline Summary (v5)

| Module | Status | Script | Description |
|--------|--------|--------|-------------|
| 01: Dataset Discovery | Complete (v1) | 01_dataset_download.py | 12 datasets downloaded |
| 02: Metadata Harmonization | Complete (v1) | 02_metadata_harmonization.py | Condition mappings finalized |
| 03: Preprocessing | Complete (v1) | 03_preprocessing.py | 12 per-dataset h5ad files, ~429K cells |
| 04: Coarse Classification | Complete (v4) | 04_annotation.py | 5 coarse categories + Unknown |
| 05: Integration | In Progress (v5) | 05_integration.py | 3 parallel workflows: CCA, scANVI, STACAS |
| 06: Clustering | Pending | 06_clustering.py | |
| 07: Post-Integration Annotation | Pending | 07_annotation.py | |
| 08: Differential Analysis | Pending | 08_differential.py | |
| 09: Biological Interpretation | Pending | 09_interpretation.py | |
| 10: Trajectory Analysis | Pending | 10_trajectory.py | |
| 11: Cell-Cell Communication | Pending | 11_communication.py | |
| 12: Final Reporting | Pending | 12_reporting.py | |

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

- **Workflow A (Seurat CCA):** R-only, label-free reference pattern. Tests both flat and tiered. **Primary workflow.**
- **Workflow B (scANVI):** R+Python, semi-supervised with coarse anchor labels from Module 04. Tiered via batch hierarchy.
- **Workflow C (STACAS):** R-only, label-guided (optional coarse labels). Tiered integration compatible.

Human checkpoint selects which workflow to carry forward to Modules 06+. See `specs/05_INTEGRATION.md` for full details.

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

### v4 (2026-03-11): Spec restructuring + scANVI — COMPLETE
- Pipeline restructured from 10 to 12 modules (clustering and annotation split from integration)
- Module 04: 5 coarse anchor categories replace binary classification
- Module 05: scANVI (semi-supervised) replaces scVI (unsupervised); checkpoint resume added
- Module 06: Leiden clustering with resolution optimization; adaptive resolution count for large datasets
- Module 07: Two-stage post-integration annotation (coarse → fine); 10 NP types, 2 AF types, 3 CEP types
- Module 08: 23 powered DE comparisons; 925+ significant genes
- Module 09: 1,772 enriched pathways; 7 pain-related DE genes
- Module 10: Pseudotime in NP/AF/CEP; CEP rho=0.396
- Module 11: 39K healthy vs 37K degenerated CCC interactions
- Module 12: Final report + 27 supplementary tables
- Modules 01-03 unchanged from v1

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
| 2026-03-21 | v5 initiated. v4 results archived. data/integrated cleared. data/processed retained. |
