# IVD Analysis Plan

## Current Status

**PIPELINE RERUN IN PROGRESS.** Specs restructured (commits 7938443, 0307f65): Module 04 narrowed to binary classification, Module 05 expanded to include clustering + de novo annotation, four compartment objects replace two-tier structure. Modules 04–10 must be rerun with updated scripts.

## Active Step

**Pipeline v2 rerun complete.** Awaiting human checkpoint review.

### Condition Mapping Review (required before Module 06, per Module 02 checkpoint)

Reviewed 2026-03-05. Decisions on flagged items:

1. **Herniated samples (10 NP, from GSE233666 + GSE251686):** Keep as separate "herniated" category, do NOT fold into degenerated. Rationale: herniated tissue is mechanically disrupted and may have distinct inflammatory/repair signatures vs. in-situ degenerated tissue. The 10 samples provide enough power for herniated vs. healthy comparisons. Can also run herniated vs. degenerated as an exploratory comparison.

2. **GSE205535 NNP (11yo spinal cord injury, classified "healthy"):** Reclassify to **exclude from DE comparisons** rather than treating as healthy. An acute spinal cord injury in an 11-year-old is not representative of "healthy" disc biology — trauma response genes may contaminate healthy baseline. Keep for annotation/integration but flag as excluded in DE. GSE205535_DNP (81yo degenerated) remains as degenerated_ungraded.

3. **Thompson III boundary (GSE230809):** Currently Thompson II-III → degenerated_mild, Thompson III → degenerated_mild, Thompson III-IV → degenerated_severe. This is reasonable — the uncertainty is at the III boundary, and having it in "mild" is conservative. No change.

4. **Neonatal samples (GSE189916, n=3):** Keep as separate category. Useful for developmental comparisons but should NOT be mixed into "healthy" (neonatal disc biology is fundamentally different from adult healthy).

5. **Aged ungraded (GSE189916 adult, n=3):** These are >65yo adults with no back pain history and unknown degeneration grade. Keep as "aged_ungraded" — they could be healthy-aged or subclinically degenerated. Useful for aging analyses but should not be in the healthy vs. degenerated comparison.

6. **Degenerated ungraded (GSE205535_DNP + GSE255768 CEP, n=3):** Include in "degenerated_all" comparisons but not in mild vs. severe since grade is unknown.

**Final comparison plan for Module 06:**
- Primary: healthy (20 samples) vs. degenerated_all (mild+severe+ungraded, 42 samples)
- Secondary: healthy vs. degenerated_mild (18), healthy vs. degenerated_severe (21), mild vs. severe
- Exploratory: healthy vs. herniated (10), herniated vs. degenerated, neonatal vs. adult-healthy
- Exclude from all DE: GSE205535_NNP (trauma confound)
- Per compartment where sample counts allow (NP has best power, AF second, CEP underpowered)

## Integration Approach Decision (Module 05 Checkpoint, 2026-03-05)

**Primary: scANVI (Approach B). Sensitivity check: scVI (Approach A) for trajectory analysis.**

Rationale: No approach clearly dominated — overall scores tightly clustered (0.60-0.62). Decision based on downstream needs:

| Metric | A: scVI | B: scANVI | C: Harmony | D: BBKNN |
|--------|---------|-----------|------------|----------|
| NP overall | 0.607 | **0.618** | 0.599 | 0.614 |
| AF overall | 0.608 | **0.615** | 0.601 | 0.611 |
| NP celltype ASW | 0.502 | **0.521** | 0.459 | 0.489 |
| AF celltype ASW | 0.496 | **0.511** | 0.484 | 0.494 |
| NP clusters (0.5) | 24 | 27 | 17 | 29 |
| AF clusters (0.5) | 34 | 34 | 19 | 27 |
| NP condition acc | 0.653 | 0.629 | 0.567 | **0.928** |
| AF condition acc | 0.648 | 0.624 | 0.581 | **0.867** |
| NP score var ratio | **1.0** | 0.655 | 0.683 | 0.527 |
| AF score var ratio | **1.0** | 0.461 | 0.523 | 0.468 |
| NP study-cluster ARI | 0.229 | 0.238 | 0.184 | — |
| AF study-cluster ARI | 0.202 | 0.218 | 0.113 | — |

Key reasoning:
1. **scANVI (B) chosen as primary** — best overall score and cell type separation for both NP and AF. Semi-supervised approach leverages Module 04 annotations, producing the most refined cell type structure. Standard atlas choice.
2. **scVI (A) retained for sensitivity** — perfectly preserves cell state continuum (score variance ratio = 1.0). Important for Module 08 trajectory analysis, where overcorrection could erase gradual cell state transitions.
3. **Harmony (C) not chosen** — most aggressive correction. Fewest clusters (17 NP, 19 AF) suggests it merges real biological groups. Lowest condition accuracy means disease signal is partially erased.
4. **BBKNN (D) not chosen as primary** — highest condition accuracy (preserves disease signal well), but no corrected embedding (only corrected neighbor graph), limiting downstream flexibility. All 4 embeddings remain available in the h5ad files.
5. **No blob problem** with any approach — no overcorrection to single-cluster outputs.
6. **For Module 06 (DE):** pseudobulk approach uses cell type labels + raw counts, not embeddings. scANVI labels are the most refined, supporting this choice.
7. **For Module 08 (Trajectory):** will run on both scANVI and scVI embeddings to assess sensitivity of continuum results to integration method.

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
| Module 05: Integration (Tier 2) | 2026-03-05 | Complete | 4 approaches (scVI, scANVI, Harmony, BBKNN) run for NP (138,937 cells) and AF (282,736 cells). All validation checks pass. |
| Module 05: Human checkpoint (Tier 2) | 2026-03-05 | Approved | Primary: scANVI (B). Sensitivity: scVI (A) for trajectory. Rationale: best overall score + cell type separation; scVI preserves continuum for Module 08. |
| Module 06: Differential analysis | 2026-03-05 | Complete | Composition analysis (Mann-Whitney U): 0/58 significant cell type proportion changes. Pseudobulk DE (pyDESeq2): 17 powered comparisons run, 128 skipped (underpowered), 5,328 significant genes (|log2FC|>0.5, padj<0.05). Top results: NP_mature_chondrocyte healthy_vs_herniated (4,316 DE genes), AF_outer healthy_vs_degenerated_severe (203), AF_outer mild_vs_severe (133), Endothelial healthy_vs_herniated (414). Volcano plots + heatmaps generated. All validation checks PASS. |
| Module 06: Human checkpoint | 2026-03-05 | Approved | Decisions: (1) Herniated comparison flagged as exploratory/likely study-confounded (RPL genes in top hits, only 2 studies); exclude from primary interpretation. (2) Endothelial annotation caveat noted (ACAN/IBSP/CYTL1 suggest misclassified NP/AF cells); no re-annotation needed. (3) No additional comparisons. (4) Composition trends biologically sensible despite FDR failure. (5) No systematic batch domination in degeneration comparisons. |
| Module 07: Biological interpretation | 2026-03-05 | Complete | Part 1 ORA: 1,244 significant enrichments (GO/KEGG/Reactome/MSigDB/IVD-custom) across 8 cell type x comparison groups. ECM, inflammatory, collagen, immune pathways confirmed. Part 1b GSEA: 1,081 significant terms. IVD custom gene set heatmap shows Inflammatory_signaling, Cellular_senescence, Matrix_degradation enriched in degeneration. Part 2 TF activity: 113 significant TFs (CollecTRI regulon overlap). Key: ATF3/ATF7 in NP severe, HSF1/HSF2 across cell types (stress), NFKBIB in NP_stressed, E2F4/TFDP1 (cell cycle) in NP severe. Part 3 Pain: only 3 significant pain gene hits (TNF x2, CXCL8 x1) — pain mediators mostly below detection in disc cells (expected: disc cells produce pro-inflammatory mediators that sensitize nerves, not nociceptors themselves). All validation PASS. |
| Module 07: Human checkpoint | 2026-03-05 | Approved | Pathways consistent with known IVD biology. Novel TF findings (ATF3/7, HSF1/2) worth highlighting. Pain analysis confirms indirect signaling model. No contradictions. Proceed to Module 08. |
| Module 08: Trajectory analysis | 2026-03-05 | Complete | PAGA + DPT pseudotime for NP (50K downsampled from 139K) and AF (50K from 281K). NP root: notochordal cluster (99% NP_notochordal). AF root: AF_inner cluster. Pseudotime-condition correlation: NP rho=-0.207, AF rho=-0.177 (both p~0, negative = healthy cells at earlier pseudotime). 500 trajectory genes per compartment (NP: 417 late_up, 83 late_down; AF: 353 late_up, 147 late_down). Overlap with DE: NP 278/500, AF 254/500 trajectory genes are also DE. RNA velocity unavailable (no spliced/unspliced layers). Sensitivity check (scVI): NP rho=-0.132, consistent direction. All validation PASS. |
| Module 08: Human checkpoint | 2026-03-05 | Approved | Trajectory biologically sensible: notochordal→mature→stressed gradient in NP, inner→outer→mechanical_stress in AF. Pseudotime aligns with disease condition. ~55% DE overlap confirms consistency. RNA velocity absence documented and acceptable. |
| Module 09: Cell-cell communication | 2026-03-05 | Complete | LIANA (CellPhoneDB+NATMI+Connectome+SingleCellSignalR+log2FC consensus) on 20K cells/condition from per-dataset files. Healthy: 44,079 interactions, 17 cell types. Degenerated: 53,036 interactions, 22 cell types. 3,662-4,194 pain-relevant interactions flagged. Differential analysis: 79,654 compared. Collagen-integrin positive controls confirmed. All validation PASS. |
| Module 09: Human checkpoint | 2026-03-05 | Approved | Interactions biologically plausible. Pain-relevant interactions include neurotrophin and VEGF pathways. More interactions in degeneration (53K vs 44K) consistent with increased paracrine signaling. Proceed to Module 10. |
| Module 10: Final reporting | 2026-03-05 | Complete | Final report (results/final_report.html) with 12 sections. 13 supplementary tables collected. Requirements frozen. All 18 validation checks PASS. All module scripts (01-10) present. All intermediate reports present. |

## Pending Steps (Rerun with Restructured Specs)

Modules 01–03 unchanged. Modules 04–10 rerun with restructured specs.

1. [x] Modules 01–03: Dataset discovery, metadata, preprocessing — DONE (unchanged)
2. [x] **Infrastructure prep:** Swap reactivated (8 GB), disk verified (68 GB free)
3. [x] Module 04 (rerun): Coarse cell classification — 11 datasets classified (GSE233666 excluded), all validation PASS, 0% ambiguous across all datasets. GSE230809 and GSE242443 100% mesenchymal.
4. [x] Module 05 (rerun): Integration + clustering + de novo annotation — 4 objects (NP 263K, AF 85K, CEP 51K, all_cells 411K). scVI-only, resolution optimization, de novo annotation, CellTypist validation. All checks PASS. Completed 2026-03-10 04:58.
5. [ ] Module 05: Human checkpoint — annotation review (most critical gate) — **DEFERRED, proceeding with downstream modules per user instruction**
   - **CellTypist disagreements to review:**
     - NP: 8 discordant clusters (out of 13). Pericyte_SMC→Fibroblasts (3,970 cells), Macrophage→Endothelial (4,612 cells), 4 "unassigned" clusters called Classical monocytes/Late erythroid by CellTypist, T_cell→Classical monocytes (3,004 cells), B_cell→Double-positive thymocytes (1,521 cells).
     - AF: 1 discordant — Pericyte_SMC→Fibroblasts (257 cells).
     - CEP: 3 discordant — B_cell→T cells (58 cells), Pericyte_SMC→Endothelial (40 cells), NK_cell→T cells (53 cells).
   - NP non-mesenchymal annotation quality is notably low — many disagreements may indicate poor marker resolution at the selected clustering resolution (0.2, 13 clusters for 40K cells).
   - **Hypotheses for NP disagreements:**
     1. **Misrouted mesenchymal cells (biggest issue):** Clusters 0, 2, 4 (~17K cells, 42% of NP non-mes) are likely stressed/inflammatory disc cells misclassified by Module 04. Top markers (NAMPT, SOD2, CXCL8, HSPA1A, HLA-B) are stress-response genes IVD cells upregulate under degeneration, not canonical immune markers. De novo correctly leaves them "unassigned"; CellTypist forces "Classical monocytes" because Immune_All_Low has no IVD cell type. These cells should probably be reclassified mesenchymal.
     2. **De novo scoring formula disadvantages endothelial markers:** Cluster 1 (4,612 cells) is clearly endothelial by markers (GNG11, SPARCL1) but scored as Macrophage by de novo. The `frac_expr * 0.5 + mean_expr * 0.5` formula may favor CD68 (diffuse low IVD expression) over sparser but specific PECAM1/VWF/CDH5. CellTypist is likely correct (99.8% agreement on "Endothelial cells").
     3. **CellTypist Immune_All_Low model has no pericyte/SMC category:** Cluster 3 (3,970 cells, Pericyte_SMC vs Fibroblasts) — CellTypist lumps pericytes as "Fibroblasts." These may also be IVD fibroblast-like cells (CALD1, TPM2, IGFBP7) misrouted from mesenchymal.
     4. **Mixed clusters from insufficient resolution:** Clusters 6, 8 (T_cell by de novo but with myeloid markers MPO, CD74, LYZ) suggest merged populations. Higher resolution would separate them.
   - **Potential fixes for rerun:**
     - Tighten Module 04 classification to reduce mesenchymal→non-mesenchymal leakage (e.g., require co-expression of ≥2 immune markers, not just score-based)
     - Fix igraph rebuild bottleneck in `optimize_clustering_resolution` (build graph once, not 20x)
     - Use higher minimum resolution for non-mesenchymal tier (e.g., 0.5 floor)
     - Fix de novo scoring to weight marker specificity, not just mean expression
     - Consider using CellTypist as primary for non-mesenchymal annotation rather than the generic marker panel scoring
6. [x] **Disk cleanup:** Old tier files removed (~14.6 GB freed), 58 GB available
7. [x] Module 06 (rerun): Pseudobulk DE — 21 powered comparisons, 53 skipped. Herniated excluded (single-study confound). Key: NP_mature_chondrocyte mild_vs_severe 315 genes, NP_fibrocartilaginous mild_vs_severe 203 genes, AF_outer healthy_vs_degenerated_severe 97 genes, EP_hyaline healthy_vs_degenerated_all 84 genes. Completed 2026-03-10 05:32.
8. [ ] Module 06: Human checkpoint — deferred
9. [x] Module 07 (rerun): 1,577 sig ORA enrichments, 55K GSEA terms tested, 10 sig pain genes (PTGS2, TNF, PLA2G2A, BDKRB2, CCL2, PTGES). Completed 2026-03-10 05:47.
10. [ ] Module 07: Human checkpoint — deferred
11. [x] Module 08 (rerun): Trajectory on NP, AF, CEP. NP rho=-0.258, AF rho=+0.341 (reversed — flag for review), CEP rho=-0.163. 500 trajectory genes each. Completed 2026-03-10 06:00.
12. [ ] Module 08: Human checkpoint — deferred
13. [x] Module 09 (rerun): LIANA on integrated objects. Healthy 28,878 interactions, degenerated 27,011. 2,077 pain-relevant. 36,014 differential. Completed 2026-03-10 06:03.
14. [ ] Module 09: Human checkpoint — deferred
15. [x] Module 10 (rerun): Final report with 19 supplementary tables (incl. S17-S19 CellTypist concordance), CellTypist disagreement table in report. Completed 2026-03-10 06:03.
16. [ ] Module 10: Human checkpoint — final review — **AWAITING**

### Previous Pipeline Run (v1, 2026-02-26 to 2026-03-05)

All modules completed; see Completed Steps table below. Results superseded by spec restructuring.

## Revisions Log

- 2026-03-09: **Spec restructuring (commits 7938443, 0307f65).** Key changes: (1) Module 04 narrowed from fine-grained annotation to binary mesenchymal/non-mesenchymal classification; (2) Module 05 expanded to include clustering resolution optimization and de novo annotation (absorbs old Module 04 annotation logic); (3) Four compartment objects (NP, AF, CEP, all_cells) replace two-tier structure (tier1_nonresident, tier2_resident_NP/AF); (4) scVI-only replaces 4-approach benchmark; (5) GSE233666 excluded (herniated-only); (6) `cell_type_final` renamed to `cell_type`; (7) "resident/non-resident" terminology changed to "mesenchymal/non-mesenchymal." All modules 04-10 require rerun. Hardware upgraded: 62 GB RAM / 16 CPUs (from 30 GB / 4 CPUs), same A10G GPU.
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
- **int64 counts bloat (RESOLVED):** `layers['counts']` stored as int64 instead of int32, wasting ~50% space. Fixed in `load_subset_concat` to cast to int32. AF file dropped from 13 GB to ~11 GB (with all 4 approach embeddings).
- **Disk full crash during AF scANVI checkpoint (RESOLVED):** 61 GB EBS was insufficient. scANVI checkpoint atomic write (10 GB temp file) filled the disk at 90% usage. EBS expanded to 123 GB. AF integration restarted from scVI checkpoint successfully.
- **Tier 2 OOM kill (RESOLVED):** Approach A (scVI) for NP completed 200 epochs of training but the process was killed by the OOM killer during post-training metric computation (`compute_metrics`). Root cause: 139K-cell kNN graph + silhouette scores exceeded 16 GB RAM. Fix: (1) stratified subsampling to 30K cells for metric evaluation (consistent with scIB benchmark practice, Luecken et al. 2022 — embeddings are still computed on all cells), (2) explicit `gc.collect()` between metric computations and between integration approaches. scVI checkpoint and model were saved before the kill; resume skips Approach A.
- **Checkpoint gating failure (RESOLVED):** Modules 03-05 ran without proper human checkpoints. Root cause: agent loop did not enforce checkpoint stops; analysis_plan.md was not updated. Fixed 2026-03-03 with revised PROMPT.md and run_pipeline.sh gate.
- 2026-03-05: Module 05 Tier 2 complete. NP: 138,937 cells (including 155 EP cells), AF: 282,736 cells (including 1,862 fibroblasts). All 4 approaches (scVI, scANVI, Harmony, BBKNN) run for both compartments. Validation: 25/26 PASS (1 cosmetic FAIL for report timing, resolved). Key metrics — NP: scVI overall=0.608, scANVI=0.615, Harmony=0.601, BBKNN=0.611. AF: scVI overall=0.608, scANVI=0.615, Harmony=0.601, BBKNN=0.611. No blob problem detected. Study-cluster ARI 0.11-0.24 (good batch mixing). Condition accuracy 0.58-0.87 (biological signal preserved). Cluster counts range 17-34 at res 0.5. EBS expanded from 61 GB to 123 GB mid-run to resolve disk space crash.
- ~~**ACTION REQUIRED BEFORE MODULE 06:**~~ **RESOLVED 2026-03-05.** Condition mappings reviewed and finalized. See "Condition Mapping Review" section above. GSE205535_NNP excluded from DE, herniated kept separate, Thompson III boundary accepted.
- 2026-03-05: Module 06 completed. Composition analysis found no significant cell type proportion changes (0/58 at padj<0.05). Pseudobulk DE: 17 comparisons powered, 128 skipped. 5,328 significant genes total. Dominant results from NP_mature_chondrocyte|healthy_vs_herniated (4,316 genes — very large, may reflect biological + technical differences in herniated tissue). Positive controls: ADAMTS5 found in NP degeneration results (log2FC=0.92) but padj=1 (underpowered). MMP13 not detected in NP results (may not pass gene filtering). Heatmaps generated for AF_outer and Endothelial cells (healthy_vs_degenerated_all comparison).

## Deferred Questions

- Should spatial transcriptomics data (if found) be incorporated, and if so, how? **Update:** No human IVD spatial transcriptomics datasets were found. Zhou 2023 used mouse Visium only.
- Should the analysis include cross-species comparisons (e.g., mouse IVD data) for validation? Mouse/rat/bovine/goat datasets identified and logged in registry.
- Should the final atlas be deposited to CellxGene for community use? **Update:** No IVD data currently exists on CellxGene, HCA, or Single Cell Portal — this atlas would be the first.
- ~~Should Kuchynsky 2024 (GSE242443, culture-expanded CEP) be included?~~ **RESOLVED:** Yes, included.
- ~~Should Zhou 2023 (embryonic IVD) be included?~~ **RESOLVED:** Deferred to Module 08 (trajectory analysis).
- ~~How to handle Chinese repository datasets?~~ **RESOLVED:** CNP0002664 downloaded; PRJCA014236 and PRJCA007656 dropped (NP already well-covered).
