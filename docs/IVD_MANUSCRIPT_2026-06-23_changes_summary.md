# Manuscript Changes Since Original MANUSCRIPT.pdf

**Baseline:** `docs/MANUSCRIPT.pdf` (commit `45fc27f`, 2026-05-22)
**Current:** `docs/IVD_MANUSCRIPT_2026-06-23_combined.pdf`

---

## 1. Martin Lotz's 27 review edits (commit `d589034`, 2026-06-12)

Martin submitted a clean-edited version (`IVD_MANUSCRIPT_2026-05-22_clean.md`) with editorial and scientific revisions throughout. Key terminology change: "resident" → "mesenchymal" cells.

---

## 2. Four additive insertions responding to ML reviewer comments (commit `9e72c26`, 2026-06-12)

All are new content only — nothing from the clean version was removed.

| Comment | Location | Finding |
|---------|----------|---------|
| **ML#24** — Endothelial-admixed cells | Methods §Contamination | ~1,831 endothelial-admixed NP_fibrocartilaginous cells score like true endothelial cells (0.91), not pericytes/fibroblasts. Supports doublet/contamination framing. |
| **ML#20** — Notochordal cells | Results §1 | Mean notochordal panel scores ≈ 0 across all NP mesenchymal types. No coherent notochordal subpopulation detected (consistent with postnatal loss in humans). Caveat: neonatal GSE189916 cells excluded by NP-label filter. |
| **ML#14/25** — Contamination by condition | Caveat 3 | RBC contamination: 7.5% healthy vs 1.7% degenerated; endothelial-admixed: 3.5% healthy vs 1.3% degenerated (OR=2.68, p=8×10⁻⁹³). Pattern reflects sample-handling differences, not disease-driven vascular changes. |
| **ML#27** — Sex adjustment of NP DE | Caveat 6 | NP contrasts: 69–94% of DEGs retained under `~sex + group` — NP results are not sex-confounded. AF_outer contrasts: sex term is **not estimable** (all-male healthy reference). AF_outer §3 findings carry a residual sex confound — explicitly flagged. |

New supplementary tables: S20 (notochordal scores), S21 (endothelial-admixed panels), S22 (sex-adjustment summary).

---

## 3. Figure additions and fixes

| Commit | Change |
|--------|--------|
| `e731196` | Re-rendered Fig 1 and Fig 2 from post-harmonization h5ads (source data correction) |
| `4cecc25` | Added **Figure 2b** — mesenchymal sub-state UMAPs (§6g) |
| `ed95768` | Added Harmony integration comparison + fixed harmonypy orientation bug |
| `833c771` | Added **Figure 14** — NP integration UMAP grid (6 methods); corrected "Flat scANVI" mislabel → Flat CCA (v4) |
| `500e9fd` | Enlarged Fig 1 cluster points; matched Fig 2b layout to Fig 2 |
| `45f68e0` + `a3dcccd` | Added Flat CCA (v4) panel to Figure 14 (regenerated; embedding was not previously retained) |

---

## 4. Integration-method comparison section rewritten (commits `833c771`, `45df855`)

Methods table expanded from 3 rows (with a Flat scANVI/CCA mislabel) to **7 rows** — flat CCA v5/v4, tiered CCA v5/v4, scANVI, STACAS, Harmony — all scored on the same metric battery. A new "continuum-preservation controls" subsection was added showing that flat CCA retains more chondrogenic/fibrogenic marker variance than scANVI, supporting the tiered-CCA design rationale.

---

## Editorial flag

The Caveat 6 addition explicitly states that §3's AF_outer findings carry "an additional residual sex confound." This meaningfully weakens a prominent result. If Martin prefers to flag it only at the caveat level without cross-referencing §3, it is a single-sentence rollback.
