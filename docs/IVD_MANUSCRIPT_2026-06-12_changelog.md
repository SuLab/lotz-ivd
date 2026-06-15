# IVD Manuscript — clean → combined changelog

**Base:** `docs/IVD_MANUSCRIPT_2026-05-22_clean.md` (Martin Lotz's 27 review edits, committed by Hannah Swahn in `d589034`, 2026-06-12).

**Result:** `docs/IVD_MANUSCRIPT_2026-06-12_combined.md` — same content, plus four additive insertions answering ML reviewer comments #20, #24/#25/#14, and #27. No content from the clean version was removed.

The four insertions came from three standalone analysis scripts added in commit `0694162`:

- `scripts/analysis_ml20_notochordal.py`
- `scripts/analysis_ml24_endothelial_admixed.py`
- `scripts/analysis_ml27_sex_de.py`

All three ran against `data/integrated/tiered_v4/all_cells.h5ad` (410,705 cells × 49,623 genes) on 2026-06-12. Raw outputs in `results/ML20/`, `results/ML24/`, `results/ML27/`; supplementary-table copies staged as `S20*`, `S21*`, `S22*` under `results/supplementary_tables/`.

---

## Insertion 1 — Methods §Annotation > Contamination handling (ML#24)

A post-hoc marker-panel test of whether the 1,831 endothelial-admixed NP_fibrocartilaginous cells represent a genuine perivascular population versus doublets.

**Finding.** Endothelial-admixed cells score **0.91** on the Endothelial panel (≈ bona fide Endothelial 0.89), only **−0.02** on a Mural/Pericyte panel (RGS5, PDGFRB, ACTA2, MCAM, KCNJ8), **0.09** on an Adventitial-fibroblast panel (PI16, DPT, MFAP5, PCOLCE2), and the NP_fibrocartilaginous-panel score collapses from **1.31 (clean)** to **0.25**. Pattern is consistent with doublet/contamination, not a genuine perivascular cell state. Supports the manuscript's existing contamination framing.

**Where inserted.** Appended to the existing "Contamination handling" paragraph.

---

## Insertion 2 — Results §1, between IVD_mixed note and Figure 2 (ML#20)

A targeted marker-panel scan of the NP cell types against a notochordal-cell panel (KRT8, KRT18, KRT19, FOXA2, TBXT, CD24, CA12) plus a progenitor panel (PROCR, GDF5, CD24, TEK, ENG).

**Finding.** Mean notochordal scores are essentially zero across all three NP mesenchymal cell types (NP_fibrocartilaginous −0.010, NP_fibrochondrocyte_chondroid −0.005, NP_mature_chondrocyte −0.001). No coherent notochordal-like sub-population exists at the cell-type level — consistent with the well-described postnatal loss of notochordal cells from the human NP.

**Limitation noted in the inserted paragraph.** The single neonatal dataset (GSE189916) sits in IVD_mixed and its cells carry generic `Chondrocyte_like` / `Fibrochondrocyte_like` / `Fibroblast_like` labels rather than NP-prefixed labels, so the script's NP-restricted filter excluded them. The strongest version of the ML#20 question (do GSE189916 neonatal cells light up the notochordal panel?) was not tested; a dedicated neonatal-versus-adult notochordal analysis would require compartment re-annotation of GSE189916.

**Where inserted.** New paragraph between the IVD_mixed note and the Figure 2 image.

---

## Insertion 3 — Caveat 3 (ML#14, ML#25)

Both contamination categories are **concentrated in healthy rather than degenerated samples**.

**Findings.**

- Atlas-wide RBC contamination affects **7.5% of healthy cells versus 1.7% of degenerated cells**.
- Endothelial-admixed flag affects **3.5% of healthy NP_fibrocartilaginous cells versus 1.3% of degenerated NP_fibrocartilaginous cells** (Fisher OR = 2.68, p = 8 × 10⁻⁹³).

**Reading.** Consistent with cross-study sample-handling differences (some healthy donors are organ-donor cadaveric tissue with retained vascular content; degenerated samples are surgical and tend to be cleaner) rather than degeneration-driven changes in disc vascular architecture. The directional opposite of a neovascularization-with-degeneration narrative — worth pre-empting since reviewers may ask.

**Where inserted.** Appended to Caveat 3 (contamination flags).

---

## Insertion 4 — Caveat 6 (ML#27)

An empirical test of the sex-adjustment robustness of the manuscript-level NP DE results, complementing Martin's caveat that a fully sex-stratified DE is not powerable.

**Method.** For each NP cell-type × contrast, pyDESeq2 was fit twice — `~group` (the manuscript design) and `~sex + group` — and the overlap between significant DEG sets was measured.

**Findings.**

- Sex term estimable for **all 12 NP cell-type × contrast combinations**; estimable for **none of the AF_outer contrasts** (AF_outer's healthy reference is all-male, n=6).
- Across NP contrasts with ≥85 naive DEGs, **69–94% of naive DEGs are retained under sex adjustment** → the manuscript-level NP DEGs are not driven by sex confounding.
- Sex-adjusted designs return substantially more DEGs in several contrasts (e.g. NP_mature_chondrocyte healthy-vs-mild: 0 → 2,068). pyDESeq2 dispersion-trend warnings throughout the run identify this as a parametric→mean-based dispersion-trend fallback, not a clean gain in power. We report only the retention figure as the robust claim.
- The AF_outer §3 findings (the 121-DEG early downregulation including CXCL8, NGFR, PLA2G2A) **cannot be sex-adjusted with the present sample inventory** and therefore carry an additional residual sex confound. This is flagged explicitly in the merged Caveat 6.

**Where inserted.** Appended to Martin's rewritten Caveat 6.

---

## Style alignments applied to my insertions

- "resident" → "mesenchymal" throughout (matches Martin's terminology change).
- Backtick formatting for file paths (`results/ML24`, `~sex + group`, etc.).
- Sample-level (not cell-level) sex framing in Caveat 6, matching Martin's "48 of 78 samples" framing.
- Supplementary-table citations renumbered to S20 / S21 / S22 (existing tables end at S19).

## New supplementary tables (under `results/supplementary_tables/`, gitignored by project convention)

| ID | File | Source |
|----|------|--------|
| S20 | `S20_notochordal_score_by_celltype.csv` | ML#20 |
| S20b | `S20b_notohigh_by_study.csv` | ML#20 |
| S21 | `S21_endothelial_admixed_panel_scores.csv` | ML#24 |
| S21b | `S21b_endothelial_admixed_fraction_by_condition.csv` | ML#25 |
| S21c | `S21c_contamination_by_condition.csv` | ML#14 |
| S22 | `S22_sex_adjustment_summary.csv` | ML#27 |

## One judgement call worth Martin's attention

In the merged Caveat 6, I added an explicit sentence stating that §3's AF_outer findings carry **an additional residual sex confound** because all-male sample inventory made sex adjustment impossible. That's a meaningful weakening of one of the manuscript's striking results. If Martin prefers to leave §3 as-is and flag the issue only at the caveat level without explicitly attaching a confound to §3, the change is a single-sentence rollback.

---

# Update 2026-06-15 — Integration-method comparison: corrected + completed (§Methods)

Independent follow-up edit to Methods §Integration-method comparison. Produced by `scripts/05o_unified_np_comparison.py` (new) with full-NP scANVI and STACAS re-runs; outputs in `results/integration/{scanvi_np_flat,stacas_np_flat}/` and `results/integration/np_experiment/unified_comparison.tsv`.

**Correction (important).** The previous three-row table labelled its second row "Flat scANVI" with values iLISI 0.209 / batch_ASW 0.796 / … . Those numbers are in fact the **flat CCA v4** row from `comparison_table.tsv` — a mislabelled CCA run, not a scANVI result. There was no real scANVI row in the manuscript. The label is corrected to **Flat CCA (v4)**.

**Completion.** scANVI and STACAS had only ever been scored under a different (unnormalized-LISI) convention, and STACAS on a 16k-cell subsample, so they were previously excluded "from the scale." Both were regenerated on the full NP set and scored on the **identical** battery as the CCA/Harmony rows (`compute_metrics` from `05n`, byte-for-byte the battery behind `comparison_table.tsv`). The table is now seven rows — flat CCA v5/v4, tiered CCA v5/v4 (mesenchymal), scANVI, STACAS, Harmony — all directly comparable, with a `Scope` column and a cell-count footnote (flat CCA 262,967; Harmony/scANVI/STACAS 262,924; tiered 259,558).

**Reading.** All methods reach comparable global batch_ASW (0.80–0.90). scANVI/STACAS give the sharpest discrete structure (cLISI 0.99, top NMI/ARI) but are semi-supervised — guided by the same coarse labels those metrics score against — so that lead is partly self-fulfilling and is flagged as such. On the label-free discriminator (marker-variance retention), flat CCA preserves the most chondrogenic/fibrogenic variance and scANVI the least, reinforcing the tiered-CCA design rationale.

**New figure.** Figure 14 (`manuscript_figures/fig14_np_integration_umap_grid.png`) — NP UMAPs of six method embeddings (flat CCA v5, tiered CCA v5/v4, scANVI, STACAS, Harmony; flat CCA v4 omitted, embedding not retained), coloured by study and by coarse cell class.
