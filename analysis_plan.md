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

**Tiered v4 pipeline — Module 11 complete 2026-05-22 (cell-cell communication); manuscript draft (`docs/MANUSCRIPT.md`) in progress. Module 12 is replaced by the manuscript for this version.**

Module 06 rerun on 2026-05-22 (post-rescue, 50-cell non-mes threshold) was followed the same day by Module 07 reruns with four patches to `scripts/07_annotation.py`: (a) CellTypist input is now CP10K + log1p normalized (fixes the 2026-05-15 hard failure on CEP non-mes); (b) `annotate_coarse()` for the non-mes tier falls back to Module 04's per-cell `coarse_label` majority (≥60% threshold) when stage-1 panel scoring does not fire; (c) `process_all_cells_secondary()` Categorical-assignment bug fix (cast both sides to object dtype before transfer); (d) **new Stage 3 sub-state annotation (`annotate_subtype`) using overlap-based scoring** against `SUBSTATE_PANELS` (`proliferating`, `inflammatory`, `stressed`, `matrix_active`, `migratory`, `homeostatic` fallback) plus an endothelial-admixed contamination flag (CD34/EMCN/AQP1). Writes a new `cell_subtype` column on the mesenchymal tier; non-mes cells get `cell_subtype = cell_type`.

Validator overall status: **PASSED**. **0.0 % unassigned across all four compartments**. Cell-type counts: NP 10, AF 7, CEP 7, all_cells 18. Sub-state counts: NP 15, AF 10, CEP 10, all_cells 27. **1,831 cells flagged `_endothelial_admixed`** (NP non-mes cluster 8 — held for review). See **Tiered v4 Module 07 Rerun (2026-05-22)** below.

**Module 08 grouping decision (agreed 2026-05-22):** primary pseudobulk DE runs on `cell_type` (collapsed) for statistical power and v5 cross-version comparability; composition tests run on both `cell_type` and `cell_subtype`; `cell_subtype` is the supplementary descriptive lens for within-cell-type characterization and Module 09 contextualization, not the unit of statistical claims.

**Label harmonization + contamination flagging (2026-05-22, applied by `scripts/07b_harmonize_labels.py`):** before Module 08, the cell_type labels in all four h5ads were harmonized for compartment-prefix consistency, and two new contamination-tracking columns were added.

Renames (applied to `cell_type`, `coarse_cell_type`, and `cell_subtype` base):
- `EP_hyaline` → `CEP_hyaline`
- `EP_ossification` → `CEP_ossification`
- `Fibrochondrocyte_chondroid` → `NP_fibrochondrocyte_chondroid`
- `Fibrochondrocyte_fibroid` → `CEP_fibrochondrocyte_fibroid`
- `Fibroblast_like` → `CEP_outer` *(CEP-compartment cells only; IVD_mixed cells from GSE189916 keep generic labels since that dataset doesn't separate compartments)*

New columns:
- `is_contamination` (bool) — True for `Erythrocyte` cells and `_endothelial_admixed` sub-states
- `contamination_type` (string) — one of `RBC`, `endothelial_admixed`, `clean`

Result: 19 cell_types in all_cells (8 compartment-prefixed + 4 IVD_mixed generic + 7 non-resident). 18,345 cells flagged contamination (16,514 RBC, 1,831 endothelial_admixed) = 4.5 % of total; **decision to retain with caveats rather than filter** (per 2026-05-22 checkpoint). Module 08 will filter or annotate via the contamination columns.

The `scripts/07_annotation.py` panel definitions (`CANONICAL_MARKERS`, `FINE_PANELS`, `_get_fine_panel_for_coarse`) were updated to produce the harmonized labels natively — future end-to-end reruns will not need the 07b post-step.

Historical snapshots retained: **Tiered v4 Module 06 Rerun (2026-05-22)** for the new cluster counts; **Tiered v4 Module 06 Results (2026-05-09)** for the 2026-05-14 review decisions; **Tiered v4 Module 07 Results (2026-05-15)** for the pre-rescue annotation state.

**Pipeline v5 COMPLETE (2026-03-25).** All 12 modules finished with CCA integration. v5 results retained for comparison.

---

## DONE (2026-06-15) — Unified NP integration-method comparison (`05o`)

**Resolved 2026-06-15** (branch `docs/combine-martin-ml-analyses`). The scaffold below (added earlier the same day) was executed and extended. Final outcome:

- **Scope extended to 7 methods** — STACAS added alongside flat CCA v5/v4, tiered CCA v5/v4 (mesenchymal), scANVI, Harmony — per in-session direction to re-run **both** scANVI and STACAS on full NP.
- **scANVI** regenerated flat on the full NP set (the historical `data/integrated/scanvi/NP.h5ad` was gone); **STACAS** regenerated on the full NP set (vs. the old 16k subsample) via `05c --no-downsample` + `05e` export.
- **Single convention** guaranteed: `05o` scores scANVI/STACAS through `compute_metrics` from `05n`, verified byte-identical to `05h` (same constants, `scale=True`/`rescale=True`, Leiden res=0.5 var-ratio). The 4 CCA rows + Harmony row are reused from `comparison_table.tsv` / `harmony/NP/metrics.tsv` unchanged.
- **Issue #2 (mislabel):** confirmed — the "Flat scANVI" row was the flat CCA v4 numbers. Relabelled **Flat CCA (v4)**; the stale "not placed on the scale" exclusion sentence deleted.
- **Issue #3 (iLISI vs batch_ASW):** both reported; divergence explained in prose (most pronounced for scANVI/STACAS, which also prompted a semi-supervised-circularity caveat).
- **Outputs:** `results/integration/np_experiment/unified_comparison.{tsv,md}`; new **Figure 14** (`docs/manuscript_figures/fig14_np_integration_umap_grid.png`); manuscript table + prose updated; changelog entry added. PDF left stale (no PDF engine on host — regenerate separately).

**Follow-up (noted 2026-06-22):** Figure 14 omits the **Flat CCA (v4)** UMAP panel — its metrics row is reused from `comparison_table.tsv`, but the embedding coordinate file (`data/integrated/np_experiment/flat_v4/all/`) was cleared on disk by commit `d99b45e` and, unlike scANVI/STACAS, was not regenerated for `05o` (non-production arm; re-running a full-NP Seurat v4 SCT+CCA flat integration purely to draw one panel wasn't warranted). The caption already states this. *If a reviewer flags the table-row-without-panel asymmetry:* regenerate the flat_v4 embedding via `Rscript scripts/05g_np_experiment.R --mode all` (or the flat_v4 arm specifically) on AWS, add its path to the `EMBEDDINGS` dict in `05o`, and re-run `python3 scripts/05o_unified_np_comparison.py --stage figure`. Low priority — cosmetic completeness only; the metrics comparison is unaffected.

> Original scaffold note (retained as the planning record):

**Goal:** rebuild the manuscript's "Integration-method comparison" table (NP) so all six methods are scored on ONE metric convention. Requested by Lotz lab: include **flat CCA v5, tiered CCA v5, flat CCA v4, tiered CCA v4, scANVI, Harmony** (tiered = mesenchymal tier only).

**Why this is needed (three issues found 2026-06-15, on branch `docs/combine-martin-ml-analyses`):**
1. **Mixed metric conventions.** flat/tiered CCA (v4,v5) were scored by `05h` and Harmony by `05n` using the scib normalized battery (`ilisi_knn(scale=True)`, `silhouette_batch(rescale=True)`, 50k/30k subsample, k=90). But scANVI/STACAS were originally scored by `05d` with a *homemade* inverse-Simpson LISI + raw sklearn batch silhouette on a 5k/k=30 subsample — NOT comparable. Evidence: `analysis_plan.md` workflow table records NP CCA iLISI **3.68** / scANVI **1.23** (raw 1–N), vs. the manuscript's rescaled 0.258 / 0.209.
2. **scANVI table-vs-prose contradiction in the manuscript.** The NP table includes a `Flat scANVI` row with the full normalized battery (incl. cLISI/bio_ASW/NMI/ARI/var_ratios, which only `05h` computes), yet the prose + commit `ed95768` say scANVI was "scored under a different unnormalized-LISI convention… not placed on the scale." The row appears to have been rescored normalized; the exclusion sentence is stale. Confirm against the regenerated table, then **delete that sentence**.
3. **iLISI vs batch_ASW divergence for Harmony** (iLISI 0.126 worst, batch_ASW 0.857 best). Legitimate — iLISI is global neighborhood batch-diversity; batch_ASW is within-coarse-label batch separability. The manuscript currently cites only the favorable batch_ASW. **Report both and explain the divergence honestly.**

**Action:** run `python3 scripts/05o_unified_np_comparison.py` on AWS. It imports `05h`'s exact metric functions (single source of truth — no drift) and loads each method into the same `(embedding, metadata, counts, gene_names)` form:
- `flat_cca_v5`  ← `data/integrated/cca/bridge_export/NP`
- `tiered_cca_v5` ← `data/integrated/np_experiment/tiered_v5/mesenchymal`
- `flat_cca_v4`  ← `data/integrated/np_experiment/flat_v4/all`
- `tiered_cca_v4` ← `data/integrated/np_experiment/tiered_v4/mesenchymal`
- `scanvi`       ← `data/integrated/scanvi/NP.h5ad` (obsm `X_scanvi_*`)
- `harmony`      ← `results/integration/harmony/NP/embedding_harmony.npy` (+ obs/counts from `data/integrated/tiered_v4/NP.h5ad`)

Outputs `results/integration/np_unified_comparison/comparison_table.{tsv,md}`.

**Pre-flight risk:** commit `d99b45e` ("Prepare v5 pipeline… clear v4 results") cleared `results/` in git. If the np_experiment bridge exports or the v4 scANVI embedding were also cleared on disk, the relevant rows will report `status=missing` (the script skips loudly, never silently). Regenerate as needed: bridge exports via `Rscript scripts/05g_np_experiment.R --mode all`; scANVI via `scripts/05b_integration_scanvi.py`.

**Manuscript follow-ups (after the table returns, do NOT edit prose before then):** replace the 3-row NP table in `docs/IVD_MANUSCRIPT_2026-06-12_combined.md` with the 6-row unified table; fix issues #2 and #3 in the surrounding prose; add a footnote on the flat-vs-tiered fairness caveat (tiered rows score only the mesenchymal tier → bio/cluster/var-ratio metrics are on a different cell population than the flat rows).

---

## Tiered v4 Module 06 Rerun (2026-05-22) — Post NP-Rescue + 50-Cell Non-Mes Threshold

End-to-end rerun of Module 04 → 05k → 05m → 06 after implementing the 2026-05-18 root-cause fixes. Two consolidated upstream changes:

1. **NP unassigned rescue (commit `15944f1`, executed 2026-05-18 → 2026-05-21).** Module 04 cell-class panels extended for neutrophils, plasma cells, erythrocytes; Module 07 `COARSE_PANELS_NON_MESENCHYMAL` extended in parallel. Adaptive HVG selection in 05k honors `--hvg-cap` / `--max-cfo-product` (commit `a26525c`) to keep CCA tractable on the larger tier objects.
2. **Non-mes per-study threshold raised 5 → 50 cells (uncommitted, this update; `scripts/05k_tiered_v4_compartments.R`).** Below the new floor, studies are dropped from the non-mes tier rather than forced through a simple-merge fallback that produced study-segregated UMAP lobes (no real batch correction). AF lost `GSE199866_AF` (7 cells) to the drop. Every surviving non-mes object now clears CCA's `dims=1:50` requirement, so non-mes integration runs through the same Seurat v4 SCT + CCA path as mesenchymal.

### Tiered v4 atlas re-assembly

`scripts/05m_assemble_tiered_v4.py` regenerated per-compartment AnnData under `data/integrated/tiered_v4/{NP,AF,CEP,all_cells}.h5ad` on 2026-05-21:

| Compartment | n_cells | n_genes | Mes cells | Non-mes cells | "unknown" residual | Output size |
|---|---|---|---|---|---|---|
| NP        | 262,924 | 47,663 | 134,189 | 75,667  | 53,068 | 6.45 GB after Module 06 |
| AF        |  84,617 | 37,846 |  60,593 | 12,012  | 12,012 | 2.59 GB |
| CEP       |  50,854 | 32,956 |  27,748 | 13,975  |  9,131 | 1.27 GB |
| all_cells | 410,705 | 49,623 | 236,999 | 101,747 | 71,959 | 10.82 GB |

Non-mes cell counts vs the 2026-05-09 assembly: NP 3,393 → 75,667; AF 56 → 12,012; CEP 71 → 13,975; all_cells 3,464 → 101,747. The rescue moved the bulk of previously-mis-tiered immune / endothelial / blood cells into the non-mes tier where the coarse panels can actually name them.

### Module 06 results (Leiden resolution sweep, 2026-05-22)

`scripts/06_clustering.py` re-run with the same equal-weighted silhouette+modularity selection (PR #4, commit `d44a0c4`). Full resolution sweep (0.1 → 2.0, step 0.1) per tier; argmax of `(sil_norm + mod_norm) / 2`, ties broken on silhouette.

| Compartment | Tier | Cells | Resolution | n_clusters | Combined | Notes |
|---|---|---|---|---|---|---|
| NP        | mes     | 187,257 | 0.4 |  9 | 0.668 | sil = 0.040 (genuinely positive); modularity + silhouette agree at the peak |
| NP        | non-mes |  75,667 | 0.3 |  9 | 0.628 | broad combined plateau (0.625–0.628 across res 0.3–0.4); tier is now substantive |
| AF        | mes     |  72,605 | 0.3 | 12 | 0.615 | sil ~0.013; modularity-driven peak |
| AF        | non-mes |  12,012 | 0.4 |  7 | 0.876 | clear curvature; sil = 0.085 — highest of any AF tier in the rerun |
| CEP       | mes     |  36,879 | 0.4 |  9 | 0.531 | flat combined across res 0.3–0.5; sil = 0.045 |
| CEP       | non-mes |  13,975 | 0.7 | 10 | 0.582 | sil = 0.040 at selection; sil-peak alternative is res = 0.1 (3 clusters) |
| all_cells | mes     | 308,958 | 0.7 | 19 | 0.538 | plateau across res 0.4–0.9 (0.523–0.538); modest shift from 17 (2026-05-09) |
| all_cells | non-mes | 101,747 | 0.2 |  5 | 0.834 | highest combined score of any tier; sil = 0.095 — solid selection |

Combined cluster count per compartment (after merging tiers with M / NM prefix):

| Compartment | v5 (CCA flat) | Tiered v4 (2026-05-09) | Tiered v4 (2026-05-22, this rerun) |
|---|---|---|---|
| NP        | 12 | 32 | **18** (9 mes + 9 non-mes) |
| AF        | 12 | 17 | **19** (12 mes + 7 non-mes) |
| CEP       |  9 | 20 | **19** (9 mes + 10 non-mes) |
| all_cells | 15 | 21 | **24** (19 mes + 5 non-mes) |

Mesenchymal-tier clustering became coarser in NP (27 → 9) and CEP (15 → 9): removing the mis-tiered non-mes contamination flattens modularity earlier, exposing the true mes-only signal. AF mes stayed similar (13 → 12); all_cells mes shifted modestly (17 → 19). Non-mes tiers, which previously had only 14 cells/cluster on average in AF/CEP, now have meaningful substructure for the first time.

### Validation (all PASS)

- No cluster collapses to a single blob across the four compartments.
- Study identity does not predict cluster identity: max study × leiden ARI = **0.102 (CEP)**; NP 0.074, AF 0.047, all_cells 0.056 — all well below 1.0.
- Soft warning: validator reports "missing comparison resolution columns" for all four compartments (the side-by-side res=0.5/1.0 columns are not retained). Non-blocking.

### Review materials (2026-05-22)

- `results/integration/tiered_v4/clustering_resolution_optimization/{NP,AF,CEP,all_cells}_{mesenchymal,non_mesenchymal}_optimization.png` — silhouette / modularity / cluster-count plots (overwritten)
- `results/integration/tiered_v4/clustering_resolution_optimization/*_resolutions.tsv` — full resolution sweeps (overwritten)
- `data/integrated/tiered_v4/{NP,AF,CEP,all_cells}.h5ad` — clustering written back
- `notebooks/06_clustering.ipynb` §5 — refreshed to 2026-05-22 cluster counts, tables, UMAPs, optimization plots
- Logs: `logs/05k_nonmes_rerun_2026-05-21.log`, `logs/05k_all_cells_tiered_v4_2026-05-2{0,1}*.log`, `logs/05m_assemble_tiered_v4_2026-05-21.log`, `logs/06_clustering_tiered_v4_2026-05-21.log`

**Status: Module 06 cluster counts accepted (2026-05-22) — questions 1–3 in `notebooks/06_clustering.ipynb` §5f resolved in favor of the statistical selection (keep 9 mes clusters for NP, default to the chosen resolution everywhere, don't compare per-compartment vs unified cluster sums at this stage). Phase-5 gates pending evaluation after Module 08 lands.**

---

## Tiered v4 Module 07 Rerun (2026-05-22)

Module 07 was re-run on the 2026-05-22 Module 06 outputs (NP 18 · AF 19 · CEP 19 · all_cells 24 clusters) with three patches to `scripts/07_annotation.py`:

1. **CellTypist `.X` normalization.** `validate_immune_with_celltypist()` now builds a CP10K + log1p copy of the non-mes subset (without mutating the caller's `.X`, since 05m writes raw counts and downstream DE needs them raw) before handing it to `celltypist.annotate()`. The 2026-05-15 hard failure on CEP non-mes (`Invalid expression matrix in .X`) is resolved; CellTypist runs cleanly on all three non-mes tiers.

2. **Module 04 → Module 07 label propagation (non-mes tier only).** `annotate_coarse()` falls back to Module 04's per-cell `coarse_label` majority (≥60% threshold) when stage-1 panel scoring on rank_genes_groups markers fails to fire. This catches two recurrent failure modes the pre-propagation v1 rerun exposed:
   - **Naive/resting lymphocyte clusters dominated by ribosomal genes** (NP non-mes cluster 1, 17,282 cells, 86% Module 04 `Immune`; AF non-mes cluster 0, 5,122 cells, 98% Module 04 `Immune`). `rank_genes_groups` returns RPS/RPL/EEF1A1/EEF1B2 as cluster-distinguishing markers — none of which are cell-type-specific — so the existing T_cell / B_cell / NK_cell panels never fire. Propagated to `Immune`.
   - **Tier-wide RBC contamination** (CEP non-mes clusters 0–4, 12,566 cells, 72–94% Module 04 `Erythrocyte`; AF non-mes cluster 5, 370 cells, 96% Module 04 `Erythrocyte`). When most of the tier is RBCs, hemoglobin genes don't appear as cluster-vs-rest markers because they're shared. Propagated to `Erythrocyte`.

3. **Categorical-assignment bug fix in `process_all_cells_secondary()`.** The all_cells secondary pass crashed at the `celltypist_prediction` transfer step (`TypeError: Cannot set a Categorical with another, without identical categories`) when re-running on h5ads that already contained a categorical `celltypist_prediction` column from a prior run. Fix: cast both LHS and RHS to object dtype before the assignment. The v3 resume completed in ~3.5 min.

### Results

| Object | Cells | Cell types | Unassigned % | CellTypist concordance |
|---|---|---|---|---|
| NP        | 262,924 | 10 | **0.0%** | 4/9 |
| AF        |  84,617 |  7 | **0.0%** | 2/7 |
| CEP       |  50,854 |  7 | **0.0%** | 3/10 |
| all_cells | 410,705 | 18 | 0.0% | n/a (transfer mode) |

Unassigned trajectory across runs:

| Run | NP | AF | CEP | all_cells | Notes |
|---|---|---|---|---|---|
| 2026-05-15 (pre-rescue) | 12.4% | 0.0% | 0.0% | 0.2% | NP `unassigned` was mis-tiered non-mes cells; AF/CEP non-mes too small for unassigned anything |
| 2026-05-22 v1 (post-rescue, pre-propagation) | 6.6% | 6.5% | 24.7% | 0.0% | Rescue moved cells into non-mes tier; stage-1 panel scoring couldn't name them all |
| **2026-05-22 v2 (post-propagation, current)** | **0.0%** | **0.0%** | **0.0%** | **0.0%** | Module 04 majority-vote fallback closes the gap |

CellTypist concordance unchanged between v1 and v2 — the patch only changes how unassigned-by-panel clusters get labeled, not what CellTypist itself predicts. Concordance lows (NP 4/9, AF 2/7, CEP 3/10) are driven by two systematic mismaps in `Immune_All_Low.pkl`: neutrophils called as classical monocytes (the model has no neutrophil class) and pericytes / matrix-remodeling macrophages called as fibroblasts. De novo marker calls are the trustworthy labels in these cases.

### Per-compartment cell type composition

NP (10 types): `NP_mature_chondrocyte`, `NP_fibrocartilaginous`, `Macrophage_M2`, `Macrophage_M1`, `Neutrophil`, `Immune`, `Endothelial`, `Pericyte_SMC`, `Erythrocyte`, plus one within-coarse fine type.

AF (7 types): `AF_outer`, `AF_inner`, `AF_mechanical_stress` (within Fibroblast_like), `Macrophage_M2`, `Pericyte_SMC`, `Erythrocyte`, `Immune`.

CEP (7 types): `EP_hyaline`, `EP_ossification`, `Macrophage_M2`, `Pericyte_SMC`, `Endothelial`, `Erythrocyte`, plus one additional immune subtype.

all_cells (18 types): union across compartments, transferred from per-compartment objects.

### Human checkpoint items (2026-05-22 — scientific decisions, not pipeline issues)

1. **RBC contamination filter.** CEP carries 12,566 `Erythrocyte` cells (24.7% of CEP), NP carries ~16K, AF ~0.5K. CEP samples sit next to vascularized bone, so heavy RBC contamination is expected. Pseudobulk DE on RBCs is meaningless and they can distort proportion analyses. Decision: filter before Module 08, retain as a contamination-tracking marker, or keep with caveats?

2. **`Immune` catch-all (Module 04 propagation).** Clusters labeled `Immune` rather than a specific subtype (Macrophage / Neutrophil / T_cell / etc.) are those where the fine-panel scoring did not fire — typically ribosomal-dominant naive/resting lymphocytes (NP cluster 1, 17,282 cells; AF cluster 0, 5,122 cells). Decision: leave as `Immune`, or expand fine panels (e.g., resting-T-cell rule with CD3D + low-activation markers, or a "Lymphocyte_naive" panel keyed on RPS/RPL + CD45) for one more iteration of Module 07?

3. **Coarse-type naming convention (deferred from 2026-05-15).** Compartment-specific labels (`NP_fibrocartilaginous`, `AF_inner`, `EP_hyaline`) coexist with generic labels (`Fibroblast_like`, `Chondrocyte_like`) in the `all_cells` object. Harmonize before Module 08 pseudobulk?

### Review materials

- `data/integrated/tiered_v4/{NP,AF,CEP,all_cells}.h5ad` — annotated objects (`cell_type`, `coarse_cell_type`, `cell_type_confidence`, `annotation_evidence` in `.obs`)
- `results/integration/tiered_v4/cell_type_definitions.tsv` — cell-type marker definitions
- `results/integration/tiered_v4/umap_{NP,AF,CEP,all_cells}_annotated.png`
- `results/integration/tiered_v4/cluster_markers/` — 13 per-tier and within-coarse marker tables
- `results/integration/tiered_v4/annotation_dotplots/` — canonical-marker dot-plot PDFs
- `results/integration/tiered_v4/annotation_report.html` — auto-generated annotation report
- `notebooks/07_annotation.ipynb` §6 — refreshed to 2026-05-22 v2 state, including the three checkpoint items
- Logs: `logs/07_annotation_tiered_v4_2026-05-22.log` (v1, pre-propagation), `..._v2.log` (v2 main run that crashed on all_cells), `..._v3.log` (all_cells resume after categorical fix)

**Status: pending human checkpoint review of the three scientific items above before resuming at Module 08.**

---

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

## Tiered v4 Module 08 Results (2026-05-22)

Pseudobulk DE and composition analysis run via `scripts/08_differential.py` against `data/integrated/tiered_v4/` after Module 07 + harmonization. Validator overall status: **PASSED**.

Script changes for this rerun (committed alongside results):
- Added `--input-dir` / `--output-dir` / `--exclude-contamination` CLI flags (mirroring the Module 07 pattern, so v5 and tiered_v4 outputs don't trample each other).
- Refactored `run_composition_analysis()` to iterate per-object (NP / AF / CEP / all_cells) rather than combining via a single `combined_obs` — the Module 07 transfer step duplicates per-compartment cells in all_cells, so the original pooled approach would have double-counted.
- Added a second composition pass on `cell_subtype` (in addition to `cell_type`) per the Module 08 grouping decision.
- Skipped `all_cells` in the DE pass (same double-counting concern).
- Added `--exclude-contamination` flag (default False) to honor the "retain with caveats" decision; sets are dropped via `is_contamination` column if enabled.
- **Bug fix in `run_de_analysis()` index alignment.** The original code did `obs_reset.reset_index()` then `merged.set_index('index')`, but the obs index is named `cell_id` (not unnamed), so the rename column became `cell_id` and the set_index call silently failed — every cell index then ended up as integer RangeIndex, the position lookup in `pseudobulk_aggregate` returned empty arrays, and **every comparison got mis-categorized as "no valid samples with ≥50 cells"**. Fix: capture `obs_reset.columns[0]` (whatever the original index name was) and use that for `set_index`.

### Composition results (both grouping levels)

| Grouping | Total tests | Significant (padj < 0.05) |
|---|---:|---:|
| `cell_type` | 112 | **0** |
| `cell_subtype` | 182 | **0** |

No FDR-significant composition shifts at either granularity. Consistent with v5 / v4 (logged as a "Robust Findings (all 3+ versions)" entry in memory) — degeneration acts at the gene-expression level rather than the compositional level in this dataset.

### DE results — 27 powered comparisons (significant DEGs at |log2FC| > 0.5 & padj < 0.05)

| Cell type | h_vs_all | h_vs_mild | h_vs_severe | mild_vs_severe |
|---|---:|---:|---:|---:|
| **NP_fibrocartilaginous** | 84 | 14 | 291 | **325** |
| **NP_fibrochondrocyte_chondroid** | — | 5 | 27 | **350** |
| **NP_mature_chondrocyte** | 2 | 3 | 121 | 174 |
| **AF_outer** | — | 121 | 3 | 1 |
| **AF_inner** | — | 2 | — | — |
| **Immune** (naive lymphocyte catch-all) | — | 4 | 38 | 11 |
| **Macrophage_M2** | 1 | 1 | **5,659 ⚠️** | 171 |
| **Macrophage_M1** | — | — | — | 21 |
| **Neutrophil** | 3 | 2 | 23 | 25 |
| **Erythrocyte** | 0 | 0 | — | — |

Total: 59 candidate comparisons skipped underpowered (mostly small non-mes populations and the entire CEP DE matrix). CEP had **0 powered comparisons** — matches the v5 sample-count limitation.

### Headline findings

- **NP_fibrocartilaginous + NP_fibrochondrocyte_chondroid + NP_mature_chondrocyte mild_vs_severe** = 849 sig DEGs combined; the disc degeneration response is dominated by these three NP cell types responding to severity escalation.
- **Reproducibility check vs v4:** NP_fibrocartilaginous mild_vs_severe = 325 here vs 305 in v4 — within 7%, validating the dominant NP degeneration signature across pipeline versions.
- **Immune cluster (Module 04 propagation catch-all) healthy_vs_severe = 38 sig genes** — the rescue produced a real biological signal that would have been invisible in the pre-rescue pipeline.
- **AF_outer healthy_vs_mild = 121 but healthy_vs_severe = 3** — interesting AF biology: the AF transcriptional response peaks at the mild stage rather than escalating with severity. Worth flagging at the manuscript stage.

### ⚠️ Caveat: Macrophage_M2 healthy_vs_severe is statistically inflated

`Macrophage_M2 healthy_vs_severe = 5,659 significant genes` exceeds the >5,000 inflation threshold the spec explicitly warns about. Only 7 samples in that contrast (1 healthy ref + 6 severe test, or similar split), so donor-level variance dominates condition-level variance and DESeq2 over-calls DEGs. **Do not cite the M2 healthy_vs_severe row as a M2 finding in the manuscript.** The trustworthy M2 contrast is `mild_vs_severe` (9 samples → 171 sig genes), which is in the expected range for a substantive immune response. Module 09 should down-weight or exclude the inflated row when computing pathway / TF enrichments to avoid spurious dominance.

### Contamination cell handling

Erythrocyte and `_endothelial_admixed` cells were retained in DE/composition per the 2026-05-22 checkpoint decision, with their `is_contamination` / `contamination_type` columns preserved in the obs of each h5ad. Erythrocyte healthy_vs_all = 0 and healthy_vs_mild = 0 (expected — anucleate cells have no real gene regulation; non-zero DE here would itself be a red flag). Downstream interpretation should treat any DE involving these cell types as descriptive / contamination-tracking rather than biological.

### Review materials

- `results/differential/tiered_v4/composition_analysis.tsv` — cell_type composition (112 tests, 0 sig)
- `results/differential/tiered_v4/composition_analysis_subtype.tsv` — cell_subtype composition (182 tests, 0 sig)
- `results/differential/tiered_v4/de_summary_table.tsv` — 27-row powered-comparison summary
- `results/differential/tiered_v4/de_results_combined.tsv` — 1.03 M gene-test rows (180 MB)
- `results/differential/tiered_v4/de_results/` — 27 per-comparison TSV files
- `results/differential/tiered_v4/skipped_comparisons.tsv` — 59 skipped contrasts with reasons
- `results/differential/tiered_v4/volcano_plots/` — 27 volcano PNGs
- `results/differential/tiered_v4/heatmaps/` — top-DE-gene heatmaps per powered comparison
- `notebooks/08_differential.ipynb` §2 — refreshed for 2026-05-22 tiered_v4 results
- Logs: `logs/08_differential_tiered_v4_2026-05-22.log` (v1, hit the index bug — composition only completed), `..._v2.log` (DE-only resume after the index fix)

**Status: Module 08 complete. Ready to proceed to Module 09 (pathway enrichment, TF inference, pain gene cross-referencing) on the 26 trustworthy DE comparisons (excluding M2 healthy_vs_severe).**

---

## Tiered v4 Module 09 Results (2026-05-22)

Pathway enrichment, TF activity inference, and pain gene cross-referencing via `scripts/09_interpretation.py` against the Module 08 outputs. Wall time **1 h 9 min** (20:13:14 → 21:22:47 UTC). Validator overall status: **PASSED**.

Script changes for this rerun (committed alongside results):
- Added `--input-dir` / `--output-dir` CLI flags (mirroring Modules 07 and 08).
- Added `INFLATED_CONTRASTS = {('Macrophage_M2', 'healthy_vs_degenerated_severe')}` set + filter in `load_de_results()` to drop the inflated M2 contrast before any enrichment / TF / pain analysis. Confirmed in the run log: **31,481 gene-tests removed** from the input pool before Part 1.

### Headline counts (tiered v4 vs v5)

| Part | v5 (memory) | Tiered v4 | Δ |
|---|---:|---:|---:|
| ORA (Enrichr GO/KEGG/Reactome) — significant pathways at FDR < 0.05 | 2,506 | **3,051** | **+22%** |
| GSEA (pre-ranked, gseapy + custom IVD gene sets) — significant pathways | 3,301 | **6,890** | **+109%** |
| TF activity (decoupler + CollecTRI) — significant TFs | 288 | **555** | **+93%** |
| Pain-associated DE genes — unique significant | 10 | **18** | **+80%** |

Tiered v4 systematically delivers 1.2–2× the signal of v5. The scaling reflects three contributing factors:
1. A larger and cleaner trustworthy DE pool (1,820 vs 1,198 sig DEGs).
2. Compartment-prefixed cell types resolve per-compartment biology that the v5 generic labels averaged across (e.g., AF_outer vs AF_inner now distinct, NP_fibrocartilaginous separated from generic Fibroblast_like).
3. The rescued Immune lymphocyte cluster (17K NP + 5K AF cells) contributes signal that was invisible pre-rescue.

### Validator expected-pathway sanity checks

All four expected pathway terms detected in ORA results:
- `extracellular matrix` ✓
- `inflammatory` ✓
- `collagen` ✓
- `immune` ✓

### Soft warnings worth noting

- `gseapy` logged `Duplicated values found in preranked stats: 13–20% of genes` for several comparisons. This is the gene-rank-tie warning that pseudobulk DESeq2 output triggers (many genes have padj saturating at 1.0). Not blocking — it just means some GSEA NES estimates are minor approximations.
- `gseapy [ERROR] No hits returned for all input gene sets!` appeared for small-DE comparisons (1–4 sig genes — AF_inner mild, AF_outer severe, AF_outer mild_vs_severe, Macrophage_M2 healthy_vs_all, Macrophage_M2 healthy_vs_mild, Erythrocyte healthy_vs_all, Erythrocyte healthy_vs_mild). Expected — small gene lists can't hit any library at threshold.

### Key biological themes (consistent across v5 and tiered v4)

- **ECM degradation upregulated in degeneration:** MMP1/2/3/9/13, ADAMTS4/5 across NP fibrocartilaginous + mature_chondrocyte + fibrochondrocyte clusters.
- **Inflammatory and immune signaling increased:** NF-κB axis (NFKBIA, NFKBIZ, IER3, PTGS2, IL6, CXCL8, SOD2), TNF and IL1B pathways.
- **Anabolic / repair pathways downregulated:** COL2A1, ACAN, COL9A1/2/3, COL11A1/2, COMP, PRG4.
- **Pain-related cell-type-specific dysregulation:** NTN1, IL6, PLA2G2A, VEGFA, PENK, NGF, BDNF, CALCA targets cross-referenced; 18 unique genes significant across contrasts.

### Review materials

- `results/interpretation/tiered_v4/pathway_enrichment/all_enrichment_results.tsv` — ORA full table (23,357 terms)
- `results/interpretation/tiered_v4/pathway_enrichment/enrichment_*_{up,down}.png` — per-cell-type up/down dot plots
- `results/interpretation/tiered_v4/pathway_enrichment/gsea_results.tsv` — GSEA full table (123,191 terms)
- `results/interpretation/tiered_v4/pathway_enrichment/gsea_ivd_custom_heatmap.png` — IVD custom gene-set heatmap (NES per contrast)
- `results/interpretation/tiered_v4/tf_activity/tf_activity_results.tsv` — CollecTRI inferences (11,550 TF-comparison pairs)
- `results/interpretation/tiered_v4/tf_activity/tf_activity_heatmap.png` — top TFs heatmap
- `results/interpretation/tiered_v4/pain_genes.tsv` — pain gene cross-reference table
- `results/interpretation/tiered_v4/pain_genes_heatmap.png` — log2FC heatmap per cell_type × comparison
- `results/interpretation/tiered_v4/interpretation_report.html` — auto-generated HTML report
- `notebooks/09_interpretation.ipynb` §2 — refreshed for 2026-05-22 tiered v4 results (5.2 MB, all panels rendered)
- Log: `logs/09_interpretation_tiered_v4_2026-05-22.log`

**Status: Module 09 complete. Ready to proceed to Module 10 (trajectory + pseudotime).**

---

## Tiered v4 Module 10 Results (2026-05-22)

Trajectory + DPT pseudotime rerun via `scripts/10_trajectory.py` against the tiered v4 integrated atlas + Module 08 DE results. Wall time **23 min** (21:28:10 → 21:51:31 UTC). Validator overall status: **PASSED**.

Script changes for this rerun (committed alongside results):
- Added `--input-dir` / `--output-dir` / `--de-dir` CLI flags (mirroring Modules 07–09).

### Inputs

Mesenchymal tier only — selected by `cell_class` ∈ {NP/AF/CEP mesenchymal sub-states}. NP and AF downsampled to 50,000 cells each before PAGA (originally 187,257 / 72,605); CEP processed in full (36,879). Neighbor graph built on `X_integrated` (30 dims, scVI mesenchymal embedding from Module 05k). Root cluster chosen by max enrichment of the expected mature/inner population per compartment (NP_mature_chondrocyte → cluster 3, AF_inner → cluster 5, CEP_hyaline → cluster 0). RNA velocity unavailable — no spliced/unspliced layers in public count matrices; pseudotime is DPT-only.

### Headline correlations (vs v5)

| Compartment | v5 ordinal ρ | Tiered v4 ordinal ρ | Tiered v4 H-vs-D MWU p | Tiered v4 trajectory ∩ DE |
|---|---:|---:|---:|---:|
| NP  | −0.088 | **−0.004** (p=0.41)    | 3.0e-04 | 328 / 500 (66%) |
| AF  | +0.195 | **−0.003** (p=0.49)    | 1.1e-20 | 326 / 500 (65%) |
| CEP | +0.073 | **+0.077** (p=1.9e-49) | 1.6e-43 | 338 / 500 (68%) |

Compartment-level **ordinal correlations weaken substantially** in tiered v4 (NP and AF drop to ~0). However, **healthy-vs-degenerated MWU tests remain highly significant** in all three compartments. Reconciling these: median pseudotime shifts between conditions are small but consistent across very large samples, whereas the per-cell pseudotime is no longer monotonically ordered with the categorical health/mild/severe axis at compartment level.

### Per-cell-type pseudotime correlations (the signal that moved)

The within-compartment signal redistributes into **per-cell-type gradients with opposing directions**:
- **AF**: AF_inner ρ = **+0.54** (p≈0), AF_outer ρ = **−0.30** (p≈0). Same DPT axis, opposite biology — likely an outer vs inner residual gradient rather than a degeneration axis.
- **CEP**: CEP_hyaline ρ = **+0.17**, CEP_fibrochondrocyte_fibroid ρ = **+0.10**, CEP_outer ρ = **−0.21**. CEP_outer pulls the compartment-level ρ down; the other two pull it up — net positive but small.
- **NP**: All three NP cell types have small magnitudes (ρ ∈ {−0.040, +0.063, −0.052}), consistent with no clean compartment-level direction.

### Mechanism behind the v5→tiered v4 shift

In v5 a coarser cell-type partition let one direction dominate the compartment ordinal correlation. In tiered v4 the finer cell-type partition exposes that **opposing per-cell-type gradients sum near zero at compartment level**. This is the same redistribution mechanism already observed in DE (e.g. CXCL2 signal redistributed across finer cell types between versions). **Trajectory direction is a version-sensitive finding** (memory: "AF/CEP trajectory directions flipped between versions") and should be presented with caveats in the manuscript rather than as a strong claim.

### Trajectory-associated genes

500 genes per compartment with significant DPT correlation (FDR < 0.05). Program assignments (`late_up` / `late_down` / `stable`):
- NP: 282 late_up · 87 late_down · 131 stable. NP_fibrocartilaginous degeneration markers (COL1A1, COL3A1) are `late_up`; ACAN, COL2A1 (anabolic) are `late_down`.
- AF: 281 late_up · 219 late_down · 0 stable.
- CEP: 221 late_up · 279 late_down · 0 stable.

### Validator checks

- [PASS] NP pseudotime: 50,000 cells · NP pseudotime correlates with biological variable · NP trajectory genes: 500
- [PASS] AF pseudotime: 50,000 cells · AF pseudotime correlates with biological variable · AF trajectory genes: 500
- [PASS] CEP pseudotime: 36,879 cells · CEP pseudotime correlates with biological variable · CEP trajectory genes: 500
- [PASS] Trajectory report generated

### Review materials

- `results/trajectories/tiered_v4/paga_{NP,AF,CEP}.png` — PAGA connectivity graphs
- `results/trajectories/tiered_v4/pseudotime_by_condition_{NP,AF,CEP}.png` — distributions stratified by `condition_harmonized`
- `results/trajectories/tiered_v4/pseudotime_by_celltype_{NP,AF,CEP}.png` — per-cell-type ordering
- `results/trajectories/tiered_v4/gene_dynamics_{NP,AF,CEP}.png` — top trajectory genes vs pseudotime
- `results/trajectories/tiered_v4/pseudotime_{NP,AF,CEP}.tsv` — per-cell pseudotime
- `results/trajectories/tiered_v4/pseudotime_correlations_{NP,AF,CEP}.tsv` — all correlation/MWU tests
- `results/trajectories/tiered_v4/trajectory_genes_{NP,AF,CEP}.tsv` — 500 genes per compartment with program labels
- `results/trajectories/tiered_v4/trajectory_de_overlap_{NP,AF,CEP}.tsv` — trajectory ∩ Module 08 DE
- `results/trajectories/tiered_v4/trajectory_report.html` — auto-generated HTML report
- `notebooks/10_trajectory.ipynb` §2 — refreshed for 2026-05-22 tiered v4 results
- Log: `logs/10_trajectory_tiered_v4_2026-05-22.log`

**Status: Module 10 complete. Ready to proceed to Module 11 (cell-cell communication, LIANA).**

---

## Tiered v4 Module 11 Results (2026-05-22)

LIANA cell-cell communication rerun via `scripts/11_communication.py` against the tiered v4 integrated atlas. Wall time **2 min 22 s** (22:03:46 → 22:06:08 UTC). Validator overall status: **PASSED**.

Script changes for this rerun (committed alongside results):
- Added `--input-dir` / `--output-dir` CLI flags (mirroring Modules 07–10).
- No `--de-dir` flag and no INFLATED_CONTRASTS filter — Module 11 does not cross-reference DE outputs; pain flagging uses the internal `PAIN_LIGANDS` / `PAIN_RECEPTORS` panels.

### Inputs

Per-compartment tiered_v4 h5ads (NP, AF, CEP) loaded, concatenated, and split by `condition_harmonized` into healthy vs degenerated pools. Cell types with <50 cells per pool filtered out per condition. Each pool downsampled to **20,000 cells** for LIANA tractability (default cap, same as v5). 100 permutations. LIANA consensus rank-aggregation across CellPhoneDB, NATMI, Connectome, log2FC, sca, geometric mean.

### Headline counts (tiered v4 vs v5)

| Metric | v5 (memory) | Tiered v4 | Δ |
|---|---:|---:|---:|
| Healthy interactions | 25,537 | **66,827** | **+162%** |
| Degenerated interactions | 34,208 | **76,019** | **+122%** |
| Pain-relevant interactions | 3,075 | **6,883** | **+124%** |
| Differential union | n/a | 90,650 | — |

Tiered v4 yields ~2.5–3× more interactions per condition. This is expected — tiered v4 resolves more cell types (19 vs 16), and per-cell-type interaction calls scale combinatorially. CCC aggregate counts have been version-sensitive across all prior versions (memory note); cell-type-pair-level signals are more informative than aggregate count differences.

### Top differential interactions

**Gained in degeneration** (positive `rank_diff`, dominant pattern): pan-cell-type → **Neutrophil FN1 → C5AR1** and **RPS19 → C5AR1** recruitment signals. This is the canonical complement-driven neutrophil chemotaxis axis being activated across most resident cell types in the degenerated state.

**Lost in degeneration** (negative `rank_diff`, dominant pattern):
- **CEP_outer WNT2B → FZD4/LRP5/6** (CEP_outer autocrine and to NP_mature_chondrocyte)
- **CEP_outer LAMB3 → ITGA6/ITGAV/ITGA2 integrins** (CEP-to-AF basement-membrane signaling)
- **CEP_outer NDP → FZD4** (Norrie disease protein WNT signaling)
- **AF_outer / NP_fibrocartilaginous → AF_outer / NP_mature_chondrocyte HLA-DPA1/DMA → CD4** (immune-mediated signaling losing weight in degeneration)
- **CEP_fibrochondrocyte_fibroid MDK → SDC3** (midkine signaling)

Interpretation: WNT signaling at the CEP / IVD interface and HLA-class-II → CD4 immune presentation are coordinately downregulated; the dominant gained signals are neutrophil-recruitment ligands across resident cells.

### Top pain ligands (degenerated)

VEGFA (1,485 interactions), FGF2 (1,157), NGF (360), SEMA3C (300), GRN (285), TNFSF10 (250), SEMA3A (250), VEGFC (246), PTGS2 (224), VEGFB (216). Neovascularization (VEGFA/B/C, FGF2) and neurotrophin (NGF) ligands dominate the pain-relevant interaction count.

### Validator checks

- [PASS] Healthy interactions: 66,827 · collagen positive control detected
- [PASS] Degenerated interactions: 76,019 · collagen positive control detected
- [PASS] Pain-relevant interactions flagged: 6,883
- [PASS] 5 interaction plots generated · communication report generated

### Review materials

- `results/communication/tiered_v4/interactions_healthy.tsv` (11 MB) — full LIANA output for healthy pool
- `results/communication/tiered_v4/interactions_degenerated.tsv` (13 MB) — full LIANA output for degenerated pool
- `results/communication/tiered_v4/differential_interactions.tsv` (31 MB) — 90,650 union interactions with rank_diff
- `results/communication/tiered_v4/pain_interactions.tsv` (1.4 MB) — pain-flagged subset
- `results/communication/tiered_v4/interaction_plots/{interaction_heatmap,top_interactions}_{healthy,degenerated}.png` — per-condition plots
- `results/communication/tiered_v4/interaction_plots/differential_interactions.png` — differential plot
- `results/communication/tiered_v4/communication_report.html` — auto-generated HTML report
- `notebooks/11_communication.ipynb` §2 — refreshed for 2026-05-22 tiered v4 results
- Log: `logs/11_communication_tiered_v4_2026-05-22.log`

**Status: Module 11 complete. Manuscript drafting (`docs/MANUSCRIPT.md`) is the final step for this pipeline version — Module 12 reporting is replaced by the manuscript.**

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

## Tiered v4 Pipeline — Compartment Metrics + Module 06 (2026-05-05)

The 2026-04-17 NP-only switch plan was extended to all four compartments after the tiered v4 integration was run for AF, CEP, and all_cells (commits `7d157be`, `2204069`). Per-compartment metrics are at `results/integration/{af,cep,all_cells}_experiment/comparison_table.tsv`; bar charts and tables in `notebooks/05_integration.ipynb` §6.

### Per-compartment metric direction (tiered_v4 vs baseline_flat_v5)

| Compartment | iLISI ↑ | batch_ASW ↑ | cLISI ↓ | bio_ASW ↑ | condition_ASW (closer to 0 = better) |
|---|---|---|---|---|---|
| NP (mes)    | 0.216 vs 0.258 | 0.861 vs 0.850 | **0.729 vs 0.869** | **0.510 vs 0.417** | **−0.020 vs −0.165** |
| AF (mes)    | 0.091 vs 0.056 | 0.839 vs 0.829 | 0.880 vs 0.955 | 0.482 vs 0.465 | +0.024 vs +0.049 |
| CEP (mes)   | **0.241 vs 0.089** | **0.949 vs 0.832** | **0.640 vs 0.881** | 0.496 vs 0.505 | **−0.042 vs −0.103** |
| all_cells (mes) | 0.200 vs 0.152 | 0.857 vs 0.841 | **0.690 vs 0.835** | **0.508 vs 0.425** | −0.144 vs −0.150 |

CEP shows the largest improvement; AF the smallest (only 3 studies → less batch effect to remove). all_cells improves on most axes but condition_ASW barely budges. NP confirms the original 2026-04-17 finding.

### Tiered v4 atlas assembly

`scripts/05m_assemble_tiered_v4.py` merges per-tier bridge files into per-compartment AnnData under `data/integrated/tiered_v4/{NP,AF,CEP,all_cells}.h5ad`:

- `obs['cell_class']` ∈ {mesenchymal, unknown, non_mesenchymal} (already set in bridge metadata)
- `obs['tier']` records mes vs non-mes membership
- `obsm['X_integrated']` (n × 50) — per-tier PCA, NaN-padded across tiers (clustering happens per-tier)
- `obsm['X_umap']` (n × 2) — per-tier UMAP, NaN-padded
- `X` = raw counts (CSR sparse, gene union via `anndata.concat(join='outer')`)

| Compartment | n_cells | n_genes | Output |
|---|---|---|---|
| NP        | 262,951 | 49,623 | NP.h5ad — assembled to 1.34 GB, 6.45 GB after Module 06 |
| AF        | 84,568  | 37,846 | AF.h5ad — 0.52 GB → 2.59 GB |
| CEP       | 50,840  | 32,956 | CEP.h5ad — 0.26 GB → 1.27 GB |
| all_cells | 410,643 | 49,623 | all_cells.h5ad — 2.21 GB → 10.82 GB |

(File sizes grow ~5× post-Module 06 because `06_clustering.py` writes back without compression; the original assembly used gzip on the sparse counts.)

AF has only a mesenchymal tier — non-mes cells were too few per study after splitting (the smallest objects fell below the integration anchor threshold). AF non-mesenchymal cells are present only via the all_cells non-mes tier.

### Module 06 results (Leiden resolution sweep, refreshed 2026-05-09)

`scripts/06_clustering.py` was extended with `--input-dir` / `--output-dir` flags so it could be aimed at the tiered v4 outputs without overwriting v5. The 2026-05-08/09 reruns apply the equal-weighted silhouette+modularity selection introduced in PR #4 (commit `d44a0c4`):

- Full resolution sweep (0.1 → 2.0, step 0.1) is run for every tier; the cell-count-gated coarse sweeps are gone.
- Silhouette and modularity are min-max normalized within each sweep, then averaged into a combined score. The selected resolution is `argmax combined`, with ties broken on silhouette.
- Earlier 2026-05-06 changes still in effect (commits `0886c50`, `562adcd`, `f8d5940`): the non-mesenchymal `min_resolution = 0.5` floor is gone; the 05k AF tier-export threshold was raised to ≥ 5 cells per study, which created an AF non-mesenchymal tier (56 cells) for the first time; non-mes tiers below the CCA dims threshold fall back to `integrate_simple`; 05m auto-detects on-disk tiers per compartment.

#### Cluster counts (current, 2026-05-09)

| Compartment | Tier | Cells | Resolution | n_clusters | Combined | Notes |
|---|---|---|---|---|---|---|
| NP        | mes     | 259,558 | 0.8 | **27** | 0.543 | combined flat across res 0.6–1.0 (0.515–0.543); silhouette ≲0.05 → modularity-driven |
| NP        | non-mes |   3,393 | 0.3 |  5 | 0.783 | sharp curvature past res=0.3; modularity moved selection from sil-only res=0.1 (2 clusters) |
| AF        | mes     |  84,568 | 0.4 | 13 | 0.640 | silhouette barely above zero; res=0.6 (15 clusters, 0.634) within ~1% of peak |
| AF        | non-mes |      56 | 1.0 |  4 | —     | tier created by 05k threshold raise; ~14 cells/cluster — interpret with caution |
| CEP       | mes     |  50,769 | 0.9 | 15 | —     | shifted from prior 6 clusters at res=0.2 under PR #4; integration metrics favor finer granularity |
| CEP       | non-mes |      71 | 1.0 |  5 | —     | collapses to 1 cluster at every res ≤0.9 — smallest res that splits, not a curvature minimum |
| all_cells | mes     | 407,179 | 0.5 | 17 | 0.541 | combined flat across res 0.4–0.9 (0.512–0.541); v5-comparable count |
| all_cells | non-mes |   3,464 | 0.3 |  4 | —     | sharp curvature past res=0.3; shifted from prior 5 clusters at res=0.5 |

Combined cluster count per compartment (after merging tiers with M / NM prefix):

| Compartment | v5 (CCA flat) | Tiered v4 (2026-05-09) |
|---|---|---|
| NP        | 12 | **32** |
| AF        | 12 | **17** |
| CEP       |  9 | **20** |
| all_cells | 15 | **21** |

#### Validation (all PASS)

- No cluster collapses to a single blob across the four compartments.
- Study identity does not predict cluster identity: max study × leiden ARI = 0.079 (CEP); all four compartments well below 1.0.
- Comparison resolutions (0.5, 1.0) stored per tier where computed.

### Review decisions (2026-05-14)

Cluster counts accepted as-is for all four compartments. Rationale: the goal is a statistics-based resolution choice; biological merging happens at the Module 07 annotation checkpoint where DE-evidence within coarse groups (Stage 2) can collapse near-duplicate clusters. Specific notes carried forward:

- **NP mes (27 clusters)**: combined-score plateau across res 0.6–1.0 means selection is "thin but not wrong"; 17–29 clusters all defensible within noise. Granularity question deferred to Module 07 — if many of the 27 annotate to the same `cell_type` with no distinguishing DE markers, they merge there.
- **AF mes (13 clusters)**: roughly v5-comparable count. 13-vs-15 is genuinely within noise; committing to 13 ahead of Module 07.
- **AF non-mes (56 cells, 4 clusters)** and **CEP non-mes (71 cells, 5 clusters)**: single-study tiers; will almost certainly fail Phase-5 `max_study_pct < 85%` and pseudobulk power gates. Flagged for likely drop at Module 08+; AF and CEP non-mes cells will be analyzed only through `all_cells.h5ad` in practice.
- **CEP mes (15 clusters)** and **all_cells (17 mes + 4 non-mes)**: accepted; integration metrics support CEP's finer granularity, and all_cells mes lands at v5-comparable count.

The unified-vs-compartment asymmetry (NP+AF+CEP mes sum = 55 fine clusters; all_cells mes = 17) is expected: cross-compartment cell-type differences dominate the unified atlas's degrees of freedom, while compartment-level atlases resolve within-tissue substates. Both views are kept; downstream usage is per-compartment for DE, unified for Figure 1 / proportions.

### Phase-5 gates still pending (per `docs/np_switch_to_tiered_v4_plan.md`)

The original NP plan included DE-concordance and pain-gene recovery gates that were never run because we paused at this step. Phase 5 should be widened to all four compartments before declaring tiered v4 the primary atlas. Specifically:

- ≥ 80% top-100 DE concordance with v5 per powered comparison (or biologically coherent divergence)
- ≥ 5 of 10 v5 significant pain genes recovered
- max_study_pct < 85% per cluster
- pseudobulk power ≥ v5's 17 powered comparisons

These cannot be evaluated until Modules 07–08 run.

### Review materials

- `results/integration/tiered_v4/clustering_resolution_optimization/{NP,AF,CEP,all_cells}_{mesenchymal,non_mesenchymal}_optimization.png` — silhouette / modularity / cluster-count plots
- `results/integration/tiered_v4/clustering_resolution_optimization/*_resolutions.tsv` — full resolution sweeps
- `data/integrated/tiered_v4/{NP,AF,CEP,all_cells}.h5ad` — clustering written back
- Logs: `logs/05m_assemble_*.log`, `logs/06_clustering_*_tiered_v4.log`, `logs/06_clustering_validate_tiered_v4.log`

**Human checkpoint cleared 2026-05-14. Ready to resume at Module 07.**

---

## Tiered v4 Module 07 Results (2026-05-15)

`scripts/07_annotation.py` was re-run against `data/integrated/tiered_v4/` after commit `2d183ed` added `--input-dir` / `--output-dir` flags (so the tiered v4 outputs do not overwrite the v5 production results). Run wall time 18:54 → 19:45 (~51 min). Validator overall status: **PASSED**.

### Cell-type counts (per-compartment cell_type column)

| Object | Cells | Cell types | Notes |
|--------|-------|------------|-------|
| NP        | 262,951 | 8  | 12.4% `unassigned` (32,621 cells) — soft warning, see below |
| AF        | 84,624  | 7  | 0.0% unassigned |
| CEP       | 50,840  | 6  | 0.0% unassigned |
| all_cells | 410,643 | 19 | 0.2% unassigned (908 cells); 398,359 transferred + 44,905 de novo |

Finer than v5 (NP 5 · AF 4 · CEP 7 · all_cells 16), consistent with the larger Module 06 tiered v4 cluster counts (NP 32, AF 17, CEP 20, all_cells 21).

### Validator warnings to weigh at checkpoint

1. **NP `unassigned` = 12.4%** — above the 10% spec threshold. Mesenchymal clusters that did not score above the stage-1 coarse-marker cutoff for any of Chondrocyte_like / Fibroblast_like / Fibrochondrocyte_like. Question: biologically meaningful (trajectory intermediates?) or low-quality clusters to merge into a neighbor / accept with a relaxed cutoff?
2. **CellTypist failed on CEP non-mesenchymal cells** with *"Invalid expression matrix in `.X`, expect log1p normalized expression to 10000 counts per cell."* Coarse and fine labels from the marker-based annotator are unaffected; only the automated secondary check is missing for the CEP non-mes subset. One-line fix in `07_annotation.py` (pass a normalized copy of `.X` to CellTypist).
3. **`all_cells`: no cluster marker tables found** — expected, because `all_cells` annotations are transferred from compartment-specific objects rather than recomputed. Validator does not currently know to skip this check for transfer-mode objects; consider gating the check on the presence of locally-computed clusters.

### Cell-type naming question

Compartment-specific labels (`NP_fibrocartilaginous`, `NP_mature_chondrocyte`, `AF_outer`, `AF_inner`, `EP_hyaline`) coexist with generic labels (`Fibroblast_like`, `Chondrocyte_like`, `Fibrochondrocyte_like`) in the `all_cells` object. Is this the intended cross-compartment scheme (compartment-specific identities preserved where the cell sits in its own compartment, generic where it's a transfer-mode call) or should it be harmonized before Module 08 pseudobulk?

### Review materials

- `data/integrated/tiered_v4/{NP,AF,CEP,all_cells}.h5ad` — annotated objects (`cell_type`, `coarse_cell_type`, `cell_type_confidence`, `annotation_evidence` in `.obs`)
- `results/integration/tiered_v4/cell_type_definitions.tsv` — 40 rows
- `results/integration/tiered_v4/umap_{NP,AF,CEP,all_cells}_annotated.png`
- `results/integration/tiered_v4/cluster_markers/` — 14 per-tier and within-coarse marker tables
- `results/integration/tiered_v4/annotation_dotplots/` — 9 canonical-marker dot plot PDFs
- `results/integration/tiered_v4/annotation_report.html` — auto-generated annotation report
- `notebooks/07_annotation.ipynb` §6 — embedded tables, UMAPs, dot plots, and validator notes
- Log: `logs/07_annotation_tiered_v4.log`

**Status: pending human checkpoint review of the three validator warnings + the cross-compartment naming question above before resuming at Module 08.**

### Root-cause analysis of the NP unassigned cells (2026-05-18)

Follow-up marker-scan on the 8 NP mes-tier "unassigned" Leiden clusters (`scripts/07_annotation.py` Stage-1 panel: Chondrocyte_like / Fibroblast_like / Fibrochondrocyte_like). Every cluster's top-15 markers are unambiguous **non-mesenchymal lineage markers**:

| Cluster | Cells | Top markers | Identity | Mes / Unknown (Module 04) |
|---|---|---|---|---|
| 8  | 10,779 | CXCL8, S100A8/9, FPR1, C5AR1, BCL2A1, G0S2, NAMPT, SRGN | Neutrophils / activated myeloid | 1,843 / 8,936 |
| 15 | 6,982  | LYZ, S100A8/9, LTF, BPI, CEACAM8, MNDA, CD24             | Neutrophils / granulocytes      | 793 / 6,189   |
| 18 | 3,765  | CD74, HLA-DRB1/DPA1/DPB1, FCER1G, LYZ, CTSS, PTPRC       | Macrophages / APCs              | 1,026 / 2,739 |
| 19 | 3,578  | PECAM1, EMCN, CD34, EGFL7, PODXL, AQP1, RAMP2, GNG11     | Endothelial                     | 1,256 / 2,322 |
| 21 | 2,390  | TRBC2, CD48, CCL5, RUNX3, CXCR4, REL, PTPRC              | T cells                         | 187 / 2,203   |
| 22 | 2,243  | HBB, HBA1, HBA2, HBD, AHSP, ALAS2, GYPA, HEMGN           | Erythrocytes                    | 481 / 1,762   |
| 23 | 1,014  | CD79A, IGHM, IGKC, HLA-DR/DP/DQ, SPIB, CD37              | B cells                         | 69 / 945      |
| 25 | 271    | MZB1, DERL3, FKBP11, SLAMF7, CD38, SSR4, HERPUD1         | Plasma cells                    | 13 / 258      |

74.6% of the 31,022 cells carry `coarse_label = "Unknown"` from Module 04. They are not trajectory intermediates or low-quality mesenchymal cells; they are mis-tiered non-mesenchymal cells. The mes-tier coarse panel correctly rejected them.

The same marker scan on NP non-mes cluster 0 (the 1,599 cells = 47.1% unassigned in non-mes tier) shows the **identical neutrophil signature** (CXCL8, S100A8, C5AR1, G0S2, BCL2A1, NAMPT, FPR1). AF and CEP each have **zero unassigned cells** in any tier — both compartments have much smaller immune/blood populations and the tier-routing did not leak.

### Two compounding root causes

1. **Tier-routing leak in `scripts/05k_tiered_v4_compartments.R:286`** — `mes_cells <- which(cc %in% c("mesenchymal", "unknown"))` routes every `cell_class = "unknown"` cell into the mesenchymal integration tier. This was likely intended to capture ambiguous chondrocyte/fibroblast intermediates but in practice dumps ~23K NP cells (mostly neutrophils, plus macrophages / T / B / plasma / endothelial / RBCs) into the mes tier where they cannot be labeled by chondrocyte/fibroblast panels.

2. **Panel coverage gap in `scripts/04_annotation.py` and `scripts/07_annotation.py`** — Module 04's `IMMUNE_SUPPORTING` list covers T cells (CD3D/E), monocytes/macs (CD68, CD14, CSF1R), B cells (CD79A, MS4A1), mast cells (KIT, TPSAB1), NK cells (NKG7, GNLY) — but **no neutrophil markers, no plasma-cell markers, no erythroid markers**. Neutrophils have only weak PTPRC and none of the supporting markers, so they fall through every rule into "Unknown." Module 07's `COARSE_PANELS_NON_MESENCHYMAL` has the same gap, so even if the routing fix were made in isolation, NP non-mes cluster 0 would still be "unassigned."

### Proposed fix (Module 04 + Module 07; rerun cost is low)

**A. Module 04 (`scripts/04_annotation.py`):** Extend the non-mesenchymal classification rules so neutrophils / plasma cells / erythrocytes are routed to `cell_class = non_mesenchymal` rather than `unknown`. Minimal-change option is to extend `IMMUNE_SUPPORTING`:

```python
IMMUNE_SUPPORTING = [
    "CD3D", "CD3E", "CD68", "CD14", "CSF1R",
    "CD79A", "MS4A1", "KIT", "TPSAB1", "NKG7", "GNLY",
    "S100A8", "S100A9", "FCGR3B", "CSF3R",    # neutrophil / granulocyte
    "MZB1", "JCHAIN", "SDC1",                 # plasma cell
]
```

…and add a separate Erythrocyte rule before the chondrocyte/fibroblast scoring (HBB / HBA1 / HBA2 / GYPA in top decile → `cell_class = non_mesenchymal`, `coarse_label = "Erythrocyte"`), so these clearly-contaminant cells can be filtered downstream rather than smuggled into the mes tier.

**B. Module 07 (`scripts/07_annotation.py:81-89`):** Add matching coarse panels so the non-mes tier can name them after integration:

```python
COARSE_PANELS_NON_MESENCHYMAL = {
    "Macrophage":   ["CD68", "CD14", "CSF1R", "CD163", "CD86"],
    "Neutrophil":   ["S100A8", "S100A9", "FCGR3B", "CSF3R", "FPR1"],
    "T_cell":       ["CD3D", "CD3E", "CD4", "CD8A"],
    "B_cell":       ["CD79A", "MS4A1"],
    "Plasma_cell":  ["MZB1", "DERL3", "SDC1", "CD38"],
    "NK_cell":      ["NKG7", "GNLY"],
    "Mast_cell":    ["KIT", "TPSAB1"],
    "Endothelial":  ["PECAM1", "VWF", "CDH5"],
    "Pericyte_SMC": ["ACTA2", "RGS5", "PDGFRB"],
    "Erythrocyte":  ["HBB", "HBA1", "HBA2", "GYPA"],
}
```

**C. Optional `scripts/05k_tiered_v4_compartments.R` tightening:** Once (A) re-routes the bulk of these cells out of `unknown`, the line-286 dump becomes lower-stakes. If residual `unknown` cells remain (e.g. <5% per compartment), keeping them in the mes tier as ambiguous intermediates is defensible. If a larger Unknown bucket persists, route through both tiers and pick the higher-scoring panel.

**Practical recommendation:** drop cluster 22 (2,243 erythrocytes) entirely as RBC tissue-prep contamination — they have no analytic value for the DE/CCC/trajectory questions and the hemoglobin signal can distort integration. The remaining ~28.8K cells get clean non-mesenchymal identities after the fix.

**Decision points for the human checkpoint:**
1. Approve the Module 04 + Module 07 panel extensions above? (Re-run Module 04 → re-run 05k tiered_v4 integration → re-run 06/07. Roughly half a day end-to-end; reuses checkpointed scANVI / anchorset where possible.)
2. Filter the ~2.2K erythrocyte cluster pre-integration?
3. Open question deferred from main checkpoint section: cross-compartment naming scheme (compartment-specific vs. generic) — still needs a call before Module 08.

Review materials added since the original 2026-05-15 entry:
- This memo (above)
- Marker scan output retained in conversation; can be exported to `results/integration/tiered_v4/cluster_markers/unassigned_diagnostic.tsv` if desired

---

## NP Integration Quality Experiment (2026-04-17)

Follow-up experiment addressing the over-integration concern: does flat CCA on the full NP object erase the chondrocyte ↔ fibrocartilaginous continuum? Four integration strategies compared on the 262,967-cell NP object using expanded metrics (iLISI, batch_ASW, cLISI, bio_ASW, condition metrics, Leiden-vs-coarse_label NMI/ARI, and marker-variance preservation for COL2A1/ACAN/SOX9/COL1A1).

Scripts: `scripts/05g_np_experiment.R` (runs the 3 experimental arms), `scripts/05h_np_experiment_metrics.py` (metric computation). Bridge-count export bug repaired via `scripts/05i_repair_v4_bridge_counts.R`.

| Run (mesenchymal scope) | iLISI↑ | batch_ASW | cLISI↓ | bio_ASW↑ | var_COL1A1↑ | var_COL2A1↑ |
|---|---|---|---|---|---|---|
| baseline_flat_v5 (v5 primary) | 0.258 | 0.850 | 0.869 | 0.417 | **0.839** | **0.679** |
| tiered_v5 (v5 with mes/non-mes split) | 0.209 | 0.869 | 0.799 | 0.455 | 0.799 | 0.674 |
| flat_v4 (v4 SCT + CCA, no split) | 0.209 | 0.796 | 0.888 | 0.500 | 0.558 | 0.651 |
| tiered_v4 (v4 SCT + CCA, with split) | 0.216 | 0.861 | **0.729** | **0.510** | 0.552 | 0.630 |

**Takeaways:**
- **Marker-variance preservation (the continuum signal)** separates v5 from v4: COL1A1 variance ratio collapses from ~0.80 (v5) to ~0.55 (v4). Tiering within a method has a small effect; the v5 vs v4 axis dominates. **Supports retaining v5 CCA as primary for NP.**
- **Cell-type purity (cLISI)** improves with tiering in both v4 and v5 (0.87→0.80 and 0.89→0.73), but cLISI was already adequate under baseline_flat_v5.
- **Bio_ASW and batch_ASW** favor v4 methods modestly, but at the cost of continuum preservation.
- Non-mesenchymal-scope metrics (3,393 cells): tiered_v5 and tiered_v4 give similar values; small tier so less informative for the main question.

**Initial conclusion (2026-04-17 AM):** no change to the v5 pipeline — based on cluster-based `var_ratio_*` metric showing v5 preserving marker variance better.

**Revised conclusion (2026-04-17 PM) after follow-up controls:** SWITCH NP primary from flat v5 CCA to tiered v4 CCA.

Four follow-up controls (cluster-free KNN variance, pooled Moran's I, within-study Moran's I, Leiden resolution sweep; script `scripts/05j_continuum_control_metrics.py`) falsified the original claim. v5 partially flattens within-donor marker spatial structure that v4 preserves (~20–30% lower Moran's I, confirmed not a between-study batch artifact). For the atlas's DE-between-conditions goal, tiered v4 is better on every decision-relevant metric:

- condition_ASW: −0.020 (tiered v4) vs −0.165 (flat v5) — condition signal preserved
- cLISI: 0.729 vs 0.869 — cleaner cell-type purity
- bio_ASW: 0.510 vs 0.417 — better cell-type separability
- Per-study COL1A1 Moran's I: 0.653 vs 0.491 — within-donor gradient preserved
- n_clusters at res=0.5 (mes): 18 vs 13 — more powered DE comparisons

An earlier draft claimed flat v5 also won on UMAP coherence, but close inspection of the tiered_v4 mes UMAP shows all 8 studies mixed through the main mass — the fragmentation claim was overstated for tiered_v4 (flat_v4 does fragment more visibly, which is one argument for tiered over flat within v4).

Execution plan: [`docs/np_switch_to_tiered_v4_plan.md`](docs/np_switch_to_tiered_v4_plan.md). AF and CEP remain on v5 (not in dispute). Status: **proposal pending execution approval**. Includes Phase 5 DE-concordance gate before declaring tiered v4 primary; v5 results archived under `results_v5_np_cca/` for comparison.

Raw metrics: `results/integration/np_experiment/comparison_table.tsv`, `continuum_knn_var_ratio.tsv`, `continuum_sweep.tsv`, `continuum_within_study_morans_i.tsv`.

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
