# IVD Analysis Plan

## Current Status

Module 01 complete. Awaiting human review of dataset list before proceeding to Module 02.

## Active Step

**HUMAN CHECKPOINT: Approve dataset list.** Review the dataset registry, search log, and summary report before proceeding to Module 02 (Metadata harmonization).

## Completed Steps

| Step | Date | Outcome | Notes |
|------|------|---------|-------|
| Spec writing | 2026-02-26 | Complete | 11 spec files written (00_PROJECT + 01-10 modules) |
| Spec review & approval | 2026-02-26 | Approved | Human approved all specs; proceeding to execution |
| Module 01: Dataset discovery | 2026-02-26 | Complete | 13 datasets included (11 GEO downloaded, 3 Chinese repos pending); 6 new datasets found beyond original 8; ~533K cells total |

## Pending Steps

1. [x] Human review and approval of specs — DONE 2026-02-26
2. [x] Module 01: Dataset discovery & acquisition — DONE 2026-02-26
3. [ ] Module 01: Human checkpoint — approve dataset list ← **ACTIVE**
4. [ ] Module 02: Metadata harmonization
5. [ ] Module 02: Human checkpoint — approve condition mappings
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

## Known Issues

- **Chinese repository datasets not yet downloaded:** CNP0002664 (CNGB), PRJCA014236 (GSA-Human), PRJCA007656 (NGDC) require separate access/registration. These are 3 of the 13 included datasets.
- **GSE205535 (Li Z 2022)** has published corrections/corrigenda — needs careful review during preprocessing.
- **Platform heterogeneity:** 3 datasets use non-10x platforms (BD Rhapsody, Singleron Matrix) which may require platform-aware batch correction during integration.
- **CEP coverage is limited:** Only 2 endplate datasets (GSE160756 has 2 CEP samples; GSE255768 has 2 CEP samples). Compartment-specific endplate analysis may be underpowered.
- **Borderline datasets need reviewer decision:** GSE242443 (culture-expanded CEP) and Zhou 2023 (embryonic IVD).

## Deferred Questions

- Should spatial transcriptomics data (if found) be incorporated, and if so, how? **Update:** No human IVD spatial transcriptomics datasets were found. Zhou 2023 used mouse Visium only.
- Should the analysis include cross-species comparisons (e.g., mouse IVD data) for validation? Mouse/rat/bovine/goat datasets identified and logged in registry.
- Should the final atlas be deposited to CellxGene for community use? **Update:** No IVD data currently exists on CellxGene, HCA, or Single Cell Portal — this atlas would be the first.
- Should Kuchynsky 2024 (GSE242443, culture-expanded CEP) be included despite in vitro expansion?
- Should Zhou 2023 (embryonic IVD, 177K cells) be included despite being fetal tissue?
- How to handle the 3 Chinese repository datasets (CNP0002664, PRJCA014236, PRJCA007656) — can we obtain access?
