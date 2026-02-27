# IVD Analysis Plan

## Current Status

Module 02 complete. WAITING FOR HUMAN REVIEW of condition mappings and metadata harmonization.

## Active Step

**Module 02: Human checkpoint** — Review harmonized metadata, condition mappings, and summary statistics.

## Completed Steps

| Step | Date | Outcome | Notes |
|------|------|---------|-------|
| Spec writing | 2026-02-26 | Complete | 11 spec files written (00_PROJECT + 01-10 modules) |
| Spec review & approval | 2026-02-26 | Approved | Human approved all specs; proceeding to execution |
| Module 01: Dataset discovery | 2026-02-26 | Complete | 13 datasets included (11 GEO downloaded, 3 Chinese repos pending); 6 new datasets found beyond original 8; ~533K cells total |
| Module 01: Human checkpoint | 2026-02-26 | Approved | Decisions: include GSE242443 (culture-expanded CEP); defer Zhou 2023 (embryonic) to Module 08; proceed without NGDC datasets; coverage adequate |
| Module 02: Metadata harmonization | 2026-02-26 | Complete | 78 samples harmonized across 12 studies, 57 donors; cell counts from curated_metadata.xlsx; 3 low-cell-count samples flagged |

## Pending Steps

1. [x] Human review and approval of specs — DONE 2026-02-26
2. [x] Module 01: Dataset discovery & acquisition — DONE 2026-02-26
3. [x] Module 01: Human checkpoint — approve dataset list — DONE 2026-02-26
4. [x] Module 02: Metadata harmonization — DONE 2026-02-26
5. [ ] Module 02: Human checkpoint — approve condition mappings ← **ACTIVE (WAITING FOR HUMAN REVIEW)**
6. [ ] Module 03: Per-dataset preprocessing
7. [ ] Module 03: Human checkpoint — review QC reports
8. [ ] Module 04: Per-dataset annotation
9. [ ] Module 04: Human checkpoint — approve cell type labels
10. [ ] Module 05: Integration strategy (multiple approaches)
11. [ ] Module 05: Human checkpoint — choose integration approach (CRITICAL)
12. [ ] Module 06: Differential analysis
13. [ ] Module 06: Human checkpoint — review DE results
14. [ ] Module 07: Biological interpretation
15. [ ] Module 07: Human checkpoint — evaluate findings
16. [ ] Module 08: Trajectory analysis
17. [ ] Module 08: Human checkpoint — evaluate trajectory validity
18. [ ] Module 09: Cell-cell communication
19. [ ] Module 09: Human checkpoint — review interactions
20. [ ] Module 10: Reporting
21. [ ] Module 10: Human checkpoint — final review

## Revisions Log

- 2026-02-26: Module 01 execution. Searched 7 databases with 8+ query combinations. Found 6 datasets not in the original known list.
- 2026-02-26: Module 01 checkpoint. Human decisions: (1) GSE242443 included despite culture expansion, (2) Zhou 2023 embryonic data deferred to Module 08 trajectory analysis, (3) proceed without PRJCA014236 and PRJCA007656 (NGDC), (4) coverage deemed adequate.
- 2026-02-26: Module 02 execution. Harmonized metadata for 78 samples across 12 studies. Sources: GEO SOFT metadata, full-text papers (PMC), curated_metadata.xlsx from domain expert. Per-sample cell counts obtained for 53/78 samples. Key decisions: (1) GSE165722 Pfirrmann grades corrected (paper says II-V, not GEO's I-IV), (2) herniated samples classified as "herniated" not "degenerated", (3) GSE244889 Pfirrmann I reclassified as "healthy" despite authors' MDD label, (4) Thompson III alone classified as "degenerated_mild" (boundary). GSE251686 platform corrected to Singleron GEXSCOPE (was incorrectly listed as 10x). 3 low-cell-count samples flagged (<500 cells).

## Known Issues

- **NGDC datasets excluded from pipeline:** PRJCA014236 (Wang 2023) and PRJCA007656 (Ling 2022) not downloaded. Both are NP-only, which is already well-covered. Could revisit if NP coverage proves insufficient.
- **GSE205535 (Li Z 2022)** has published corrections/corrigenda — needs careful review during preprocessing.
- **Platform heterogeneity:** 3 datasets use non-10x platforms (BD Rhapsody, Singleron Matrix) which may require platform-aware batch correction during integration.
- **CEP coverage is limited:** 3 endplate datasets (GSE160756: 2 samples, GSE255768: 2 samples, GSE242443: 2 culture-expanded samples). Compartment-specific endplate analysis may be underpowered.
- **GSE242443 (Kuchynsky 2024):** Included per human decision, but CEP cells were culture-expanded — note this caveat during interpretation.
- **Low cell count samples:** CNP0002664_Ctrl (249 cells), GSE255768_S2 (423 cells), GSE230809_AF_SP20_002 (467 cells) may not survive QC. Monitor in Module 03.
- **GSE165722 GEO grade offset:** GEO lists Pfirrmann I-IV but paper Table 1 says II-V. Used paper grades (authoritative). GEO metadata has systematic off-by-one error.
- **GSE251686 platform mismatch:** Registry says "10x Genomics" but GEO metadata indicates Singleron GEXSCOPE platform. Corrected in sample_metadata.tsv.
- **GSE230809 sex bias:** ALL 11 donors are male. Combined with this being the largest dataset (24 samples), sex-stratified analyses are limited.
- **Strong age-disease confound (GSE230809):** Healthy donors 21-27y, diseased 37-73y. Cannot separate age from disease effects in this dataset alone.
- **Missing demographics:** 18/78 samples have unknown age, 30/78 have unknown sex. Limits demographic stratification.

## Deferred Questions

- Should spatial transcriptomics data (if found) be incorporated, and if so, how? **Update:** No human IVD spatial transcriptomics datasets were found. Zhou 2023 used mouse Visium only.
- Should the analysis include cross-species comparisons (e.g., mouse IVD data) for validation? Mouse/rat/bovine/goat datasets identified and logged in registry.
- Should the final atlas be deposited to CellxGene for community use? **Update:** No IVD data currently exists on CellxGene, HCA, or Single Cell Portal — this atlas would be the first.
- ~~Should Kuchynsky 2024 (GSE242443, culture-expanded CEP) be included?~~ **RESOLVED:** Yes, included.
- ~~Should Zhou 2023 (embryonic IVD) be included?~~ **RESOLVED:** Deferred to Module 08 (trajectory analysis).
- ~~How to handle Chinese repository datasets?~~ **RESOLVED:** CNP0002664 downloaded; PRJCA014236 and PRJCA007656 dropped (NP already well-covered).
