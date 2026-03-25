# IVD Pipeline Version History

> Active pipeline state is in [`analysis_plan.md`](../analysis_plan.md).
> This file captures the complete historical record across all pipeline versions.

---

## Version Changelog

### v5 (2026-03-21 to 2026-03-25): CCA integration — COMPLETE
- Spec 05 restructured for three parallel workflows (CCA, scANVI, STACAS)
- R environment installed: Seurat 5.4.0, STACAS 2.4.1, DESeq2 1.42.1, speckle 0.99.7
- Module 05: CCA selected as primary (label-free, strongest batch mixing iLISI 1.5-3.7)
- CCA full-cell on 247GB machine (migrated 2026-03-23); scANVI + STACAS retained for comparison
- Module 06: CCA produces fewer clusters (NP 12 vs v4 62) — smoother embedding
- Module 07: NP 5 types (mature_chondrocyte 72%, fibrocartilaginous 28%), AF 4, CEP 7, all 16
- Module 08: 17 powered comparisons, 1,198 sig genes (NP_fibrocartilaginous h→s: 556 genes)
- Module 09: 2,506 enrichments, 3,301 GSEA, 288 sig TFs, 10 pain genes
- Module 10: NP rho=-0.088, AF rho=+0.195, CEP rho=+0.073
- Module 11: 25,537 healthy vs 34,208 degenerated CCC interactions
- Module 12: 19 supplementary tables, final report
- Modules 01-04 reused from v1/v4

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

---

## Key Decisions Log

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
| 2026-03-21 | R installed (Seurat 5.4.0, STACAS 2.4.1, DESeq2). SeuratDisk broken with v5 — bridge workaround. |
| 2026-03-21 | scANVI workflow complete (all 4 objects). CCA v4 approach infeasible on 62GB — rewrote for Seurat v5. |
| 2026-03-22 | CCA v5 complete (all 4 objects, downsampled for NP/all_cells). STACAS v2.4 API fixed, complete. |
| 2026-03-22 | Module 05 human checkpoint ready. 3-workflow comparison report and notebooks generated. |
| 2026-03-23 | Migrated to 247GB RAM machine. CCA script updated to remove downsampling. Re-running CCA full-cell for all objects. |
| 2026-03-25 | CCA full-cell complete (all 4 objects). Metrics computed for all 3 workflows. CCA selected as primary workflow. Modules 06-12 proceeding with CCA. |

---

## CCA Run Incident Log (v5, 2026-03-24)

1. **Single-threaded NP run killed (03:00 UTC).** Original NP run started 2026-03-23 ~20:50 with reference BLAS (single-threaded) — consumed 6+ hours of CPU time on one core with 31 cores idle. Killed and restarted after compute optimization.

2. **Compute optimization applied (03:03 UTC).** Installed OpenBLAS (multi-threaded BLAS) and added `future::plan("multicore", workers=16)` during `IntegrateLayers`. Reference BLAS replaced system-wide via `update-alternatives`. `future` plan set to sequential during data loading (fork-safety with `system2` bridge calls), multicore only during integration. No functional impact on results — identical linear algebra, parallelized execution only.

3. **First optimized run crashed silently (03:03 UTC).** `future::plan("multicore")` was set globally at script startup. Forked workers conflicted with `system2()` calls in the h5ad-to-Seurat bridge conversion during data loading. Process died with no error output. Fix: moved `plan("multicore")` to activate only during `IntegrateLayers`, reverts to `plan("sequential")` after.

4. **NP, AF, CEP completed successfully (05:19–07:09 UTC).** ~2 hours total for all three objects (vs 6h+ for NP alone single-threaded).

5. **all_cells failed at IntegrateLayers (~ 07:30 UTC).** `future.globals.maxSize` was 16 GB; all_cells exported 16.25 GB of globals to workers. Error: `The total size of the 59 globals exported for future expression ('FUN()') is 16.25 GiB. This exceeds the maximum allowed size 16.00 GiB`. Fix: bumped `future.globals.maxSize` to 200 GB (machine has 247 GB RAM).

6. **all_cells restarted (15:24 UTC).** Running with 200 GB future globals limit, OpenBLAS, 16 workers.

---

## Cross-Version Sensitivity Observations

These observations document how specific results changed across pipeline versions (v1–v5). They are maintained here as evidence for the sensitivity caveats in the final report.

### Trajectory pseudotime-condition correlations

Pseudotime-condition correlations change sign across versions, indicating sensitivity to integration method, annotation choices, and root cell selection.

| Version | Integration | NP rho | AF rho | CEP rho |
|---------|-------------|--------|--------|---------|
| v2 | scVI | — | — | -0.163 |
| v3 | scVI | -0.207 | -0.177 | +0.135 |
| v4 | scANVI | — | — | +0.396 |
| v5 | CCA | -0.088 | +0.195 | +0.073 |

CEP correlation changed sign between v2 and v3. NP weakened substantially from v3 to v5. AF flipped sign from v3 (negative) to v5 (positive). These correlations are not robust to upstream methodological choices.

### CellTypist NP disagreements

| Version | Discordant NP clusters | Notes |
|---------|----------------------|-------|
| v3 | 8/13 | scVI integration, finer NP cell types |
| v4 | — | scANVI, 10 NP types |
| v5 | — | CCA, 5 NP types (broader categories) |

CellTypist lacks IVD-specific cell types, so disagreements with de novo mesenchymal labels are expected. The number of discordant clusters depends on clustering granularity. De novo labels are retained as primary in all versions.

### CCC interaction direction

The direction of the healthy-vs-degenerated interaction count difference has varied across versions.

| Version | Integration | Healthy interactions | Degenerated interactions | Direction |
|---------|-------------|---------------------|--------------------------|-----------|
| v1 | scANVI | 44K | 53K | degenerated > healthy |
| v2 | scVI | 29K | 27K | healthy > degenerated |
| v3 | scVI | 41K | 40K | near-equal |
| v4 | scANVI | 37K | 39K | near-equal |
| v5 | CCA | 25,537 | 34,208 | degenerated > healthy |

This result is sensitive to cell type definitions, integration method, and subsampling strategy. Treat the direction as uncertain.

### AF pseudotime sign

AF consistently shows positive rho (degenerated cells at later pseudotime) across v2–v5, opposite to NP in v3. In v5, NP also shows weak negative rho (-0.088), so the sign difference persists. This may reflect genuine compartment-specific biology (AF degeneration involves different processes than NP) or root cell choice effects (AF rooted at AF_inner).
