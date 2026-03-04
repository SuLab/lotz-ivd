# IVD Analysis Plan

## Current Status

Module 05 Tier 1 (non-resident integration) approved at human checkpoint 2026-03-03. Tier 2 resident cell integration now running: approaches A–D (scVI, scANVI, Harmony, BBKNN) for NP and AF compartments.

## Active Step

**Module 05: Tier 2 Resident Cell Integration** — IN PROGRESS

Running 4 integration approaches for NP and AF compartments per spec. Will produce `data/integrated/tier2_resident_NP.h5ad` and `data/integrated/tier2_resident_AF.h5ad` with integration metrics. Next stop: Module 05 Human Checkpoint (Tier 2 review).

## Completed Steps

| Step | Date | Outcome | Notes |
|------|------|---------|-------|
| Spec writing | 2026-02-26 | Complete | 11 spec files written (00_PROJECT + 01-10 modules) |
| Spec review & approval | 2026-02-26 | Approved | Human approved all specs; proceeding to execution |
| Module 01: Dataset discovery | 2026-02-26 | Complete | 13 datasets included (11 GEO downloaded, 3 Chinese repos pending); 6 new datasets found beyond original 8; ~533K cells total |
| Module 01: Human checkpoint | 2026-02-26 | Approved | Decisions: include GSE242443 (culture-expanded CEP); defer Zhou 2023 (embryonic) to Module 08; proceed without NGDC datasets; coverage adequate |
| Module 02: Metadata harmonization | 2026-02-26 | Complete | 78 samples harmonized across 12 studies, 57 donors; cell counts from curated_metadata.xlsx; 3 low-cell-count samples flagged |
| Module 02: Human checkpoint | 2026-02-26 | Approved (tentative) | All mappings tentatively approved. **MUST revisit condition mappings before Module 06 (DE analysis)** — changes after that point require full reanalysis. |
| Module 03: Preprocessing | 2026-02-26 | Complete | 12 datasets preprocessed (436,558 cells post-QC). QC thresholds: min_genes=200, max_genes=6000, min_counts=500, max_mt=20%, Scrublet doublets. 4 datasets had 100% retention (pre-filtered input). GSE251686_NP3 excluded (corrupt matrix). |
| Module 03: Human checkpoint | 2026-03-03 | Retroactive review | Checkpoint was not properly gated during execution. QC reports reviewed 2026-03-03. Notebooks corrected and re-executed. No blocking issues found. |
| Module 04: Annotation | 2026-02-26 | Complete | Per-dataset annotation using marker-based scoring (16 signatures) + CellTypist (Immune_All_Low). Consensus labels in `cell_type_final`. No IVD reference atlas available for label transfer. |
| Module 04: Human checkpoint | 2026-03-03 | Retroactive review | Checkpoint was not properly gated during execution. Annotation notebook reviewed 2026-03-03. Notebooks corrected and re-executed. No blocking issues found. |
| Module 05: Integration (Tier 1) | 2026-03-02 | Complete | Tier 1 non-resident cells integrated with scVI: 14,566 cells from 9 studies. |
| Module 05: Human checkpoint (Tier 1) | 2026-03-03 | Approved | Human approved Tier 1 integration and retroactive review of Modules 03-04. Proceeding to Tier 2 resident cell integration. |

## Pending Steps

1. [x] Human review and approval of specs — DONE 2026-02-26
2. [x] Module 01: Dataset discovery & acquisition — DONE 2026-02-26
3. [x] Module 01: Human checkpoint — approve dataset list — DONE 2026-02-26
4. [x] Module 02: Metadata harmonization — DONE 2026-02-26
5. [x] Module 02: Human checkpoint — approve condition mappings — DONE 2026-02-26 (tentative; revisit before Module 06)
6. [x] Module 03: Per-dataset preprocessing — DONE 2026-02-26
7. [x] Module 03: Human checkpoint — retroactive review DONE 2026-03-03
8. [x] Module 04: Per-dataset annotation — DONE 2026-02-26
9. [x] Module 04: Human checkpoint — retroactive review DONE 2026-03-03
10. [x] Module 05: Integration strategy (Tier 1) — DONE 2026-03-02
11. [x] Module 05: Human checkpoint (Tier 1) — APPROVED 2026-03-03
12. [~] Module 05: Integration strategy (Tier 2) — IN PROGRESS ← **ACTIVE**
13. [ ] Module 05: Human checkpoint (Tier 2) — WAITING FOR TIER 2 COMPLETION
14. [ ] Module 06: Differential analysis
15. [ ] Module 06: Human checkpoint — review DE results
16. [ ] Module 07: Biological interpretation
17. [ ] Module 07: Human checkpoint — evaluate findings
18. [ ] Module 08: Trajectory analysis
19. [ ] Module 08: Human checkpoint — evaluate trajectory validity
20. [ ] Module 09: Cell-cell communication
21. [ ] Module 09: Human checkpoint — review interactions
22. [ ] Module 10: Reporting
23. [ ] Module 10: Human checkpoint — final review

## Revisions Log

- 2026-02-26: Module 01 execution. Searched 7 databases with 8+ query combinations. Found 6 datasets not in the original known list.
- 2026-02-26: Module 01 checkpoint. Human decisions: (1) GSE242443 included despite culture expansion, (2) Zhou 2023 embryonic data deferred to Module 08 trajectory analysis, (3) proceed without PRJCA014236 and PRJCA007656 (NGDC), (4) coverage deemed adequate.
- 2026-02-26: Module 02 checkpoint. All condition mappings tentatively approved. Human decision: revisit all mappings before Module 06 (differential expression), since changes after that point require full reanalysis. Key items to revisit: whether "herniated" should be a separate axis vs folded into degeneration severity; GSE205535 NNP (11yo spinal cord injury) classification; Thompson III boundary.
- 2026-02-26: Module 02 execution. Harmonized metadata for 78 samples across 12 studies. Sources: GEO SOFT metadata, full-text papers (PMC), curated_metadata.xlsx from domain expert. Per-sample cell counts obtained for 53/78 samples. Key decisions: (1) GSE165722 Pfirrmann grades corrected (paper says II-V, not GEO's I-IV), (2) herniated samples classified as "herniated" not "degenerated", (3) GSE244889 Pfirrmann I reclassified as "healthy" despite authors' MDD label, (4) Thompson III alone classified as "degenerated_mild" (boundary). GSE251686 platform corrected to Singleron GEXSCOPE (was incorrectly listed as 10x). 3 low-cell-count samples flagged (<500 cells).
- 2026-02-26: Modules 03-05 executed without proper checkpoint gating. The agent loop continued past Module 03 and Module 04 human checkpoints without updating analysis_plan.md or waiting for human review. Discovered 2026-03-03 during manual review.
- 2026-03-03: Retroactive review of Modules 03-05. All notebooks updated to reflect actual analysis state, re-executed with zero errors, and committed. PROMPT.md revised to enforce checkpoint gating. Shell-level gate added in run_pipeline.sh.
- 2026-03-03: Module 03 key findings: 436,558 cells post-QC across 12 datasets. 4 datasets (GSE160756, GSE165722, GSE244889, GSE242443) had 100% retention — input was pre-filtered by authors. GSE189916 had lowest retention (89.3%). Diffuse CD68 expression in 6/12 datasets (expected IVD biology). All validation checks pass.
- 2026-03-03: Module 04 key findings: consensus annotation using marker-based scoring + CellTypist. NP subtypes (notochordal, mature chondrocyte, stressed/degenerative, fibrocartilaginous), AF subtypes (inner, outer, mechanical stress), EP subtypes, and CellTypist-refined immune populations.
- 2026-03-03: Module 05 key findings: Tier 1 scVI integration of 14,566 non-resident cells from 9 studies (3 studies had no non-resident cells). Tier 2 resident cell integration not yet run — code exists but data files not generated.
- 2026-03-04: Tier 2 Approach A (scVI, NP) completed training (200 epochs, ~10.7h CPU) but OOM-killed during metric computation. Patched `compute_metrics` to use stratified subsampling (30K cells, preserving rare cell type proportions) for evaluation, with `gc.collect()` between metric steps and between approaches. Subsampling follows scIB benchmark convention (Luecken et al. 2022). Integration embeddings remain computed on all 139K cells; only metric evaluation is subsampled.

## Known Issues

- **NGDC datasets excluded from pipeline:** PRJCA014236 (Wang 2023) and PRJCA007656 (Ling 2022) not downloaded. Both are NP-only, which is already well-covered. Could revisit if NP coverage proves insufficient.
- **GSE205535 (Li Z 2022)** has published corrections/corrigenda — needs careful review during preprocessing.
- **Platform heterogeneity:** 3 datasets use non-10x platforms (BD Rhapsody, Singleron Matrix) which may require platform-aware batch correction during integration.
- **CEP coverage is limited:** 3 endplate datasets (GSE160756: 2 samples, GSE255768: 2 samples, GSE242443: 2 culture-expanded samples). Compartment-specific endplate analysis may be underpowered.
- **GSE242443 (Kuchynsky 2024):** Included per human decision, but CEP cells were culture-expanded — note this caveat during interpretation.
- **Low cell count samples:** CNP0002664_Ctrl (249 cells), GSE255768_S2 (423 cells), GSE230809_AF_SP20_002 (467 cells) — all survived QC but are small.
- **GSE165722 GEO grade offset:** GEO lists Pfirrmann I-IV but paper Table 1 says II-V. Used paper grades (authoritative). GEO metadata has systematic off-by-one error.
- **GSE251686 platform mismatch:** Registry says "10x Genomics" but GEO metadata indicates Singleron GEXSCOPE platform. Corrected in sample_metadata.tsv.
- **GSE251686_NP3 excluded:** Corrupt GEO matrix file. 5 of 6 samples processed.
- **GSE230809 sex bias:** ALL 11 donors are male. Combined with this being the largest dataset (24 samples), sex-stratified analyses are limited.
- **Strong age-disease confound (GSE230809):** Healthy donors 21-27y, diseased 37-73y. Cannot separate age from disease effects in this dataset alone.
- **Missing demographics:** 18/78 samples have unknown age, 30/78 have unknown sex. Limits demographic stratification.
- **100% QC retention in 4 datasets:** GSE160756, GSE165722, GSE244889, GSE242443 input was pre-filtered by authors. Our QC thresholds removed zero cells.
- **Diffuse CD68 expression:** CD68 expressed across many clusters in 6/12 datasets (GSE189916, GSE199866, GSE205535, GSE233666, GSE230809, GSE242443). Known IVD biology — stressed disc cells express CD68 at low levels.
- **GSE230809 metadata cell count discrepancy:** sample_metadata.tsv records 92,348 cells (from publication), but raw GEO files contain 110,556. Post-QC: 105,804. The publication numbers appear to be from a downstream analysis, not the raw data.
- **Tier 2 OOM kill (RESOLVED):** Approach A (scVI) for NP completed 200 epochs of training but the process was killed by the OOM killer during post-training metric computation (`compute_metrics`). Root cause: 139K-cell kNN graph + silhouette scores exceeded 16 GB RAM. Fix: (1) stratified subsampling to 30K cells for metric evaluation (consistent with scIB benchmark practice, Luecken et al. 2022 — embeddings are still computed on all cells), (2) explicit `gc.collect()` between metric computations and between integration approaches. scVI checkpoint and model were saved before the kill; resume skips Approach A.
- **Checkpoint gating failure (RESOLVED):** Modules 03-05 ran without proper human checkpoints. Root cause: agent loop did not enforce checkpoint stops; analysis_plan.md was not updated. Fixed 2026-03-03 with revised PROMPT.md and run_pipeline.sh gate.
- **ACTION REQUIRED BEFORE MODULE 06:** All Module 02 condition mappings were tentatively approved. Must do a final review of condition_harmonized categories, especially herniated vs degenerated classification and ambiguous cases, before running differential expression. Changes after Module 06 are expensive.

## Deferred Questions

- Should spatial transcriptomics data (if found) be incorporated, and if so, how? **Update:** No human IVD spatial transcriptomics datasets were found. Zhou 2023 used mouse Visium only.
- Should the analysis include cross-species comparisons (e.g., mouse IVD data) for validation? Mouse/rat/bovine/goat datasets identified and logged in registry.
- Should the final atlas be deposited to CellxGene for community use? **Update:** No IVD data currently exists on CellxGene, HCA, or Single Cell Portal — this atlas would be the first.
- ~~Should Kuchynsky 2024 (GSE242443, culture-expanded CEP) be included?~~ **RESOLVED:** Yes, included.
- ~~Should Zhou 2023 (embryonic IVD) be included?~~ **RESOLVED:** Deferred to Module 08 (trajectory analysis).
- ~~How to handle Chinese repository datasets?~~ **RESOLVED:** CNP0002664 downloaded; PRJCA014236 and PRJCA007656 dropped (NP already well-covered).
