# IVD Analysis Plan

## Current Status

**Pipeline v3 complete.** All 10 modules executed. Awaiting final human checkpoint review.

**Pipeline version:** v3 (annotation fix applied 2026-03-10). See `results/THREE_VERSION_SUMMARY.md` for full version history.

## Active Step

**Module 10: Human checkpoint — final review — AWAITING**

All results, notebooks, reports, and supplementary tables are current and pushed to GitHub.

---

## Pipeline Summary (v3)

| Module | Status | Key Output |
|--------|--------|------------|
| 01: Dataset Discovery | Complete | 11 datasets, ~410K cells, 71 samples |
| 02: Metadata Harmonization | Complete | Condition mappings finalized |
| 03: Preprocessing | Complete | 410,759 cells post-QC across 11 datasets |
| 04: Cell Classification | Complete (v3) | Binary mesenchymal/non-mesenchymal with evidence gate |
| 05: Integration + Annotation | Complete (v3) | 4 compartment objects: NP (263K), AF (85K), CEP (51K), all_cells (411K) |
| 06: Differential Analysis | Complete (v3) | 18 powered comparisons, 1,156 unique DE genes |
| 07: Biological Interpretation | Complete (v3) | 1,043 ORA, 1,943 GSEA, 399 sig TFs, 10 pain genes |
| 08: Trajectory Analysis | Complete (v3) | NP rho=-0.151, AF rho=+0.325, CEP rho=+0.135 |
| 09: Cell-Cell Communication | Complete (v3) | 40.2K healthy / 40.9K degenerated interactions |
| 10: Final Reporting | Complete (v3) | FINAL_REPORT.md, MANUSCRIPT.md, 27 supplementary tables |

### v3 Annotation Fix (2026-03-10)

Three changes to Module 04 cell classification to fix 17K stressed NP cells misrouted to non-mesenchymal in v2:

1. **Non-mesenchymal evidence gate:** Cells must express at least one canonical non-mesenchymal marker (PTPRC, PECAM1, VWF, CDH5, CD68, CD163) above threshold to be classified non-mesenchymal. Prevents stressed disc cells with upregulated HLA/inflammatory genes from being misclassified.

2. **ACAN/SOX9 rescue:** Cells expressing ACAN or SOX9 (core IVD markers) are rescued to mesenchymal regardless of non-mesenchymal score. These genes are essentially never expressed in immune/endothelial cells.

3. **85% cluster voting:** After per-cell classification, if >85% of cells in a Leiden cluster share the same class, the entire cluster is assigned that class. Smooths out noisy per-cell calls using neighborhood information.

Modules 04-10 were rerun after the fix. Modules 01-03 were unchanged.

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

## Integration Approach

**scVI-only** (v2+ simplification). One scVI model per compartment (NP, AF, CEP) with `batch_key='study'`. Separate models for mesenchymal and non-mesenchymal tiers within each compartment. Combined `all_cells` object for cross-compartment analyses.

The v1 4-approach benchmark (scVI, scANVI, Harmony, BBKNN) was replaced by scVI-only in v2.

---

## Key v3 Results

### Differential Expression
- **18 powered comparisons**, 56 skipped (underpowered)
- **1,156 unique significant genes** (1,447 gene-comparison pairs)
- Top: NP_mature_chondrocyte mild_vs_severe (315 genes), NP_fibrocartilaginous mild_vs_severe (203), NP_mature_chondrocyte healthy_vs_severe (172)
- Herniated comparisons excluded (single-study confound)

### Biological Interpretation
- **1,043 significant ORA enrichments** (GO/KEGG/Reactome/MSigDB/IVD-custom)
- **1,943 significant GSEA terms** (FDR < 0.05)
- **399 significant TF-condition associations** (CollecTRI regulon overlap)
- **10 significant pain genes:** PTGS2, TNF, PLA2G2A, BDKRB2, CCL2, PTGES, CXCL8, and others

### Trajectories
- PAGA + DPT for NP, AF, CEP compartments
- NP rho=-0.151, AF rho=+0.325, CEP rho=+0.135
- Trajectory-DE overlap: NP 96/500, AF 110/500, CEP 38/500

### Cell-Cell Communication
- LIANA 5-method consensus on 20K cells per condition
- Healthy: 40,187 interactions, Degenerated: 40,872 interactions
- Near-equal counts between conditions (in contrast to v1 and v2 which showed larger differences)

---

## Items Requiring SME Review

1. **Trajectory instability across versions:** Pseudotime-condition correlations change sign between pipeline versions (e.g., CEP went from -0.163 in v2 to +0.135 in v3). This sensitivity to upstream annotation choices means trajectory results should be interpreted cautiously.

2. **CellTypist NP disagreements:** 8/13 de novo NP clusters are discordant with CellTypist. CellTypist lacks IVD-specific cell types, so de novo labels are retained, but this should be acknowledged.

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
- **Platform heterogeneity:** 3 non-10x datasets (BD Rhapsody, Singleron). Handled by scVI batch correction.
- **CEP underpowered:** Only 3 CEP datasets (6 samples). Compartment-specific CEP analyses are limited.
- **GSE242443 culture-expanded:** CEP cells are culture-expanded. Included with caveats.
- **GSE230809 sex bias:** All 24 samples from male donors. Limits sex-stratified analyses.
- **GSE230809 age-disease confound:** Healthy=21-27y, diseased=37-73y. Cannot separate age from disease.
- **Missing demographics:** 18/78 samples unknown age, 30/78 unknown sex.
- **GSE251686_NP3 excluded:** Corrupt matrix file (5/6 samples retained).
- **GSE165722 Pfirrmann offset:** GEO says I-IV, paper says II-V. Paper grades used.

---

## Version History

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
| 2026-02-26 | Condition mappings tentatively approved. Revisit before Module 06. |
| 2026-03-03 | Retroactive checkpoint review of Modules 03-05. No blocking issues. |
| 2026-03-05 | scANVI primary for v1 (later superseded). Condition mappings finalized. |
| 2026-03-05 | v1 pipeline complete. All modules 01-10 done. |
| 2026-03-09 | Spec restructuring. GSE233666 excluded. scVI-only. v2 rerun initiated. |
| 2026-03-10 | v2 complete. v3 annotation fix applied. Full rerun Modules 04-10. |
| 2026-03-10 | v3 complete. Stale files cleaned. Notebooks re-executed. Reports updated. |
