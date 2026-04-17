# Plan: Switch NP Integration from flat v5 CCA to tiered v4 CCA

**Date drafted:** 2026-04-17
**Status:** Proposed — pending execution approval
**Scope:** NP compartment only (AF / CEP / all_cells not in dispute)

## Decision

Switch the primary NP integrated object from `data/integrated/cca/bridge_export/NP/`
(flat v5 CCA produced by `scripts/05a_integration_cca.R` via Seurat v5
`IntegrateLayers`) to `data/integrated/np_experiment/tiered_v4/`
(SCTransform + Seurat v4 `FindIntegrationAnchors` / `IntegrateData`, with
mesenchymal / non-mesenchymal tier split; produced by
`scripts/05g_np_experiment.R --mode tiered_v4`).

AF and CEP are NOT affected — they weren't tested in this experiment,
and their biology (AF especially) is less dominated by a single continuum.
They remain on the v5 CCA workflow.

## Reasoning (summary; full evidence in `notebooks/05_integration.ipynb` §5)

For the project's stated goals — cross-study integration, cell clustering,
and differential expression between conditions — tiered v4 beats flat v5
on every decision-relevant axis:

| Axis | flat v5 (current) | tiered v4 (proposed) | Winner |
|---|---|---|---|
| Batch mixing (iLISI) | 0.258 | 0.216 | tied |
| Cell-type purity (cLISI ↓) | 0.869 | **0.729** | tiered v4 |
| Cell-type separability (bio_ASW ↑) | 0.417 | **0.510** | tiered v4 |
| Per-tier batch_ASW (↑) | 0.850 | **0.861** | tiered v4 |
| Condition signal (condition_ASW) | **−0.165** | −0.020 | tiered v4 |
| Condition local mixing (condition_LISI ↓) | 2.52 | **2.23** | tiered v4 |
| COL1A1 Moran's I (pooled, embedding-wide) | 0.484 | **0.683** | tiered v4 |
| COL1A1 Moran's I (within-study mean) | 0.491 | **0.653** | tiered v4 |
| COL2A1 / ACAN / SOX9 Moran's I (pooled & within-study) | 0.29–0.60 | **0.34–0.60** | tiered v4 |
| n_clusters at res=0.5 (mes tier) | 13 | 18 | tiered v4 (more DE contrasts) |
| UMAP topology | smooth single blob, studies mixed | coherent main mass, studies mixed through bulk, thin edge wisps | ~tied |

Tiered v4 wins or ties on every decision-relevant metric. (An earlier
draft of this plan claimed flat v5 won on "UMAP coherence" — that
was overstated on close inspection of the tiered_v4 mes UMAP, which
is itself coherent with studies mixed through the main mass. The
flat_v4 UMAP does show more visible fragmentation, which is one
reason to prefer tiered over flat within the v4 family.)

The original `var_ratio_*` claim that v5 "preserves the continuum better"
was falsified by four follow-up controls (cluster-free KNN variance,
pooled Moran's I, within-study Moran's I batch-confound test, and Leiden
resolution sweep). The falsification holds across all four collagen /
chondrocyte markers: v5 partially flattens within-donor marker spatial
structure by ~20–30% relative to v4 (Moran's I units), and this is not
a batch artifact (confirmed by within-study KNN reconstruction).

For pseudobulk DE specifically, tiered v4's finer clusters give more
powered comparisons (v4 archived pipeline had 23 powered comparisons
vs. v5's 17), and its preserved condition signal means the DE model
starts closer to the biology.

## Why tiered rather than flat v4

- batch_ASW 0.861 (tiered) vs 0.796 (flat) — tighter batch mixing within cell types
- cLISI 0.729 (tiered) vs 0.888 (flat) — meaningfully cleaner cell-type purity
- Non-mes tier (3.4K cells) gets its own integration, so the endothelial /
  immune / pericyte-SMC separation is clean (better than tiered_v5 non-mes,
  which has an orphan cluster)
- Module 06 clustering in the v5 pipeline already expects tiered input
  (`leiden_mesenchymal` / `leiden_non_mesenchymal` columns), so no
  downstream module needs redesign
- Conceptually correct: mes and non-mes are categorically distinct
  biology; tiering prevents CCA dimensions from being dominated by the
  obvious mes-vs-non-mes axis

## Implementation plan

### Phase 0 — Archive existing v5 NP downstream results

Preserve the current v5 NP outputs so the two pipelines can be compared
after re-running:

```bash
mkdir -p results_v5_np_cca/
for d in clustering annotation differential interpretation trajectories communication final_report.html; do
  cp -r results/$d results_v5_np_cca/ 2>/dev/null || true
done
# (copy only NP-scope files; keep AF/CEP/all_cells results as-is)
```

Note: `results_v5_np_cca/` goes into `.gitignore` (matches existing
`results_v*` convention). Record the archive manifest in
`results_v5_np_cca/README.md`.

### Phase 1 — Bridge NP tiered_v4 RDS to h5ad

The tiered_v4 RDS files already exist:
- `data/integrated/np_experiment/tiered_v4/mesenchymal.rds`
- `data/integrated/np_experiment/tiered_v4/non_mesenchymal.rds`
- Plus bridge-export directories with repaired counts
  (`data/integrated/np_experiment/tiered_v4/mesenchymal/`, `.../non_mesenchymal/`)

Build a unified `NP.h5ad` for downstream modules:
1. Load both tiers' bridge exports (metadata + counts + embedding).
2. Concatenate cells preserving `cell_class` labels.
3. Set `X_integrated` = tier's integrated PCA (padded / tiered consistently).
4. Write `data/integrated/NP.h5ad` (overwrites the v5 object after backup).
5. Checksum file.

Script to write: `scripts/05k_assemble_tiered_v4_np.py`.

### Phase 2 — Module 06 (clustering) on tiered_v4 NP

Re-run `scripts/06_clustering.py` for NP only:
- `leiden_mesenchymal` clustering on mes tier
- `leiden_non_mesenchymal` clustering on non-mes tier
- Merge as in v5 pipeline

Expected cluster counts (from the experiment's Leiden sweep):
- mes tier: ~18 at res=0.5 (vs. v5's 13)
- non-mes tier: ~6 at res=0.5

Validation:
- No cluster comprises cells from a single study only (`max_per_cluster_study_pct < 90%`)
- No cluster comprises cells from a single sample only
- Each cluster has ≥ 3 samples represented (power floor for pseudobulk)

**Checkpoint**: review UMAPs and cluster composition before Module 07.

### Phase 3 — Module 07 (post-integration annotation) on tiered_v4 NP

Re-run `scripts/07_annotation.py` with the new clustering. Expect more,
finer cell types than v5's 5 NP types. Annotations should still include
the core NP types (mature_chondrocyte, fibrocartilaginous, notochordal-like
if detected, plus immune / endothelial / pericyte-SMC in non-mes tier).

Validation:
- Annotations match known NP biology (chondrocyte-like markers: COL2A1,
  ACAN, SOX9; fibrocartilaginous: COL1A1, elevated; non-mes: standard
  PECAM1, PTPRC, RGS5)
- CellTypist concordance comparable to v5

**Checkpoint**: confirm annotations before Module 08.

### Phase 4 — Modules 08–12 (DE → reporting)

Re-run 08 (DE), 09 (interpretation), 10 (trajectory), 11 (CCC), 12
(reporting) on NP only. AF and CEP results stay on v5. The final atlas
mixes v5 (AF, CEP) and tiered_v4 (NP) — document clearly in the final
report.

### Phase 5 — Comparison to archived v5 results

Before declaring tiered_v4 primary, run a focused comparison:
- **NP_fibrocartilaginous h→s DE**: concordance in gene list top-100,
  effect-size correlation, FDR distribution.
- **Pain-gene detection**: which of the 10 v5-significant pain genes are
  recovered? which are lost? any new significant pain genes?
- **Cluster biology**: do v4's ~18 mes clusters subdivide v5's 5 broad
  NP types in interpretable ways, or do they introduce study-linked
  fragmentation?

If the DE picture is comparable or better (more powered comparisons,
similar pain genes, no study-driven clusters), declare tiered_v4
primary for NP. If study-driven clusters dominate or pain-gene
recovery collapses, fall back to v5.

## Risks and mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| Study-linked clusters in Module 06 break pseudobulk DE | medium | Phase 2 validation checks; if triggered, try higher Leiden resolution that collapses study-specific sub-clusters |
| Non-mes tier annotation quality degrades | low | tiered_v4 non-mes UMAP already looks cleaner than v5 non-mes; validate at Phase 3 |
| Cross-pipeline mix (v5 AF/CEP + tiered_v4 NP) creates comparison awkwardness | low | Document clearly; per-compartment analyses are already independent in Modules 08–12 |
| Pain-gene and core findings change materially | medium | Phase 5 comparison is the explicit gate; report both versions in manuscript supplement |
| Downstream re-runs take too long | medium | Module 06 on 263K cells ~36 min in v4; Modules 07–12 cumulatively ~2–4 hours — acceptable |

## Files to write / modify

| Path | Action |
|---|---|
| `results_v5_np_cca/` | **NEW** — archive of current NP downstream outputs |
| `scripts/05k_assemble_tiered_v4_np.py` | **NEW** — bridge RDS → h5ad |
| `data/integrated/NP.h5ad` | **MODIFY** — replace with tiered_v4 assembly |
| `metadata/file_checksums.json` | **MODIFY** — update NP.h5ad checksum |
| `analysis_plan.md` | **MODIFY** — record the switch and its rationale |
| `docs/version_history.md` | **MODIFY** — add entry for the NP workflow change |
| `notebooks/05_integration.ipynb` §5d | **MODIFY** — update conclusion to recommend tiered_v4 and link here |
| `specs/05_INTEGRATION.md` | **NO CHANGE** — the NP experiment section already allows for this revision |

## Success criteria

The switch is declared a success if, after Phase 5 comparison:
- [ ] NP_fibrocartilaginous h→s DE has ≥ 80% top-100 gene concordance with v5, OR a biologically coherent divergence explainable by finer v4 clusters.
- [ ] ≥ 5 of v5's 10 significant pain genes are recovered.
- [ ] No NP cluster is dominated by a single study (`max_study_pct < 85%`).
- [ ] Pseudobulk power (number of powered DE comparisons) is ≥ v5's 17.
- [ ] Tiered_v4 mesenchymal annotations resolve chondrocyte ↔ fibrocartilaginous axis into ≥ 3 interpretable states.

## Rollback

If Phase 5 gate fails, revert `data/integrated/NP.h5ad` and
`metadata/file_checksums.json` to the pre-switch state (git-tracked), and
restore downstream results from `results_v5_np_cca/`. Keep tiered_v4
outputs archived under `results_tiered_v4_np/` for comparison/transparency
in the manuscript supplement.

## Open questions for the SME reviewer

1. Are there domain-specific cell types (e.g., notochordal remnants,
   progenitor-like states) that v5 was detecting and v4 might split in
   an interpretable way — or split incoherently?
2. Does the tiered_v4 non-mesenchymal separation into 3 clean classes
   match the expected IVD non-resident cell composition, or are there
   known subtypes that should appear as more than 3 clusters?
3. For the final manuscript: disclose both pipelines' results, or only
   the selected pipeline with a supplement note?
