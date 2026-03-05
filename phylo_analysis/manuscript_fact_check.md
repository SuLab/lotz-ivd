# Manuscript Fact-Check Report

**Date:** 2026-03-05
**Manuscript:** `draft_manuscript.md`
**Checked by:** Claude (automated validation against source data)

## Summary

Systematic validation of all quantitative claims in the draft manuscript against the underlying data files (cell metadata CSV, DEG CSV, GSEA CSVs, LIANA CSVs, Kruskal-Wallis CSV, and execution trace notebook). 6 parallel validation agents checked: cell counts, DEG counts, PAGA/trajectory, LIANA/composition, figure-text concordance, and GSEA pathway claims.

**Result:** 14 discrepancies found (6 high severity, 3 moderate, 5 low). The most serious errors involve the LLM analysis agent writing biologically plausible claims about "key genes" based on its training knowledge of IVD literature rather than the actual analysis output.

---

## HIGH SEVERITY — Factually wrong

### 1. COMP claimed as key downregulated gene — actually significantly UPREGULATED

- **Manuscript (line ~405):** Lists COMP as a key downregulated gene in severe degeneration
- **Data (`all_DEGs_severe_vs_healthy.csv`):**
  - NP_chondrocyte: LFC=+2.117, padj=0.022 (significantly UP)
  - NP_degenerative_UPR: LFC=+2.139, padj=0.010 (significantly UP)
  - NP_metallothionein: LFC=+2.340, padj=0.017 (significantly UP)
  - NP_stress_response: LFC=+2.004, padj=0.016 (significantly UP)
- **Root cause:** LLM agent likely assumed COMP is downregulated in disc degeneration (as commonly reported in OA/IVD literature) without checking the DESeq2 output

### 2. ECM organization claimed as "consistently downregulated in all 6 cell types"

- **Manuscript (lines ~415-425):** Lists ECM organization and collagen formation alongside Wnt, Notch, senescence, and RUNX as pathways "consistently downregulated in severe degeneration (all 6 cell types)"
- **Data (GSEA CSVs):**
  - REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION: significantly UPREGULATED in NP_chondrocyte (NES=+1.640, padj=0.001) and NP_HAPLN1+ (NES=+1.436, padj=0.002). Not significant in other 4 cell types.
  - REACTOME_COLLAGEN_FORMATION: significantly UPREGULATED in NP_chondrocyte (NES=+1.695, padj=0.030) and NP_HAPLN1+ (NES=+1.601, padj=0.003). Not significant in other 4.
- **Root cause:** LLM conflated individual gene downregulation (ACAN, COL2A1) with pathway-level enrichment, which also captures compensatory/remodeling gene activity

### 3. Mild degeneration: 7 donors claimed, actually 4

- **Manuscript (line ~209):** "Mild degeneration: n = 21,648 cells, 7 donors"
- **Data (`cell_metadata.csv`):** mild_degen has 21,646 cells from 4 donors: Wang_Pa17F, Wang_Pab24M, Wang_Pab30M, Wang_Pb55F
- **Note:** The total donor count of 29 is correct (the per-condition counts don't add to 29 because Cherif_d1 appears in both healthy and moderate)

### 4. PAGA methods: Harmony-corrected PCA with 30 components → actually raw PCA with 20 PCs

- **Manuscript (line ~307):** "The NP subset was re-embedded using Harmony-corrected PCA with 30 components and a kNN graph (k=20)"
- **Execution trace (notebook line 3301):** `sc.pp.neighbors(adata_np, n_neighbors=20, n_pcs=20)` using fresh PCA from lognorm expression (not Harmony). The agent explicitly abandoned Harmony because it "over-compressed the embedding" (line 3279-3288).

### 5. IL6 claimed as key upregulated cytokine — not significant anywhere

- **Manuscript (line ~398):** Lists "TNF, IL6, CXCL8" as key upregulated pro-inflammatory cytokines
- **Data:** IL6 is not significant in any cell type. Missing entirely from NP_chondrocyte, NP_HAPLN1+, and NP_stress_response DESeq2 output (filtered for low expression). Where present, padj ranges 0.13-0.98.

### 6. CXCL8 claimed as key upregulated cytokine — not significant, trends downward

- **Data:** CXCL8 is not significantly DE in any cell type. Trends downward in 4/6 cell types (LFC ranges from -2.64 to +1.84, padj 0.33-0.99).

---

## MODERATE SEVERITY — Misleading

### 7. CILP claimed as key downregulated gene — actually trends upward

- **Manuscript (line ~406):** Lists CILP as a key downregulated gene
- **Data:** CILP trends upward in all 6 cell types (LFC +1.07 to +2.24), though none reach significance (padj 0.30-0.73)

### 8. ACAN, COL2A1, HAPLN1 claimed as key downregulated genes — none statistically significant

- **Data:**
  - ACAN: LFC -0.57 to -1.63, padj 0.33-0.86 (trending down but not significant)
  - COL2A1: LFC -0.33 to -1.48, padj 0.60-0.94 (trending down but not significant)
  - HAPLN1: LFC -0.28 to +0.30, padj 0.85-0.99 (no consistent direction, nowhere near significant)
- **Note:** The direction is consistent for ACAN and COL2A1 (always negative LFC), suggesting real but underpowered effects. HAPLN1 shows no consistent pattern. The issue is presenting these as definitively DE when the statistics don't support it.

### 9. Glycolysis pathway claimed as upregulated in NP: canonical

- **Manuscript (line ~437):** "Glycolysis (NP: canonical): Metabolic reprogramming under worsening hypoxic stress"
- **Data:** No HALLMARK_GLYCOLYSIS pathway appears in any GSEA output. Closest match is REACTOME_GLUCONEOGENESIS in NP_chondrocyte (NES=+1.764, padj=0.049), which is the metabolic opposite of glycolysis.

---

## LOW SEVERITY — Minor inaccuracies

### 10-11. Mild/moderate cell counts off by 2

- Mild: manuscript says 21,648, actual 21,646
- Moderate: manuscript says 32,136, actual 32,138
- Errors are complementary (sum is correct at 173,628)

### 12. Dot plot caption: "Clusters 0-5 express NP markers"

- Cluster 4 is the AF fibroblast cluster. Should say "Clusters 0-3 and 5"

### 13. Volcano plots: "all cell types"

- Only 6 of 12 cell types are shown. Should say "all six major disc cell types" or similar.

### 14. PAGA "strongest connection"

- Should specify "strongest inter-state connection" — diagonal (self-connectivity) values also reach 0.76

---

## VERIFIED CORRECT

- All 12 cell type counts and percentages (Table 2)
- All 6 DEG total/up/down counts and ratios (Table 3)
- All 10 PAGA connectivity values (rounded from 3 to 2 decimal places)
- NP state composition percentages (healthy vs severe, all 10 values exact match)
- LIANA pair counts (24,079 healthy / 35,831 severe / 39,530 merged)
- All LIANA delta scores (TIMP1-CD63, FN1-ITGA6, FN1-C5AR1, FN1-CD44, COL1A2-CD93)
- Kruskal-Wallis statistics (H=8.45, p=0.038 for AF fibroblast)
- DPT diagnostics (DC variance 7.12e-6, Spearman rho=0.24)
- Wnt, Notch, RUNX, senescence: confirmed consistently downregulated in all 6 cell types
- TNF/NF-kB: confirmed upregulated only in NP: canonical and NP: degenerative
- EMT: confirmed upregulated only in NP: degenerative
- Collagen crosslinking: confirmed upregulated in NP: stress, NP: degenerative, NP: MT-high
- ADAMTS5: confirmed upregulated (3/6 significant, all 6 trending up)
- FN1: confirmed upregulated (4/6 significant, all 6 trending up)
- TNF: confirmed upregulated in 2/6 cell types (NP canonical, NP degenerative)
- QC retention: correctly notes M7831814 below 70% threshold (fixed in prior commit)

---

## Source Files Used for Validation

| File | Contents |
|---|---|
| `figures/05_annotation/cell_metadata.csv` | 173,628 rows, cell-level annotations |
| `figures/07_pseudobulk/all_DEGs_severe_vs_healthy.csv` | 89,457 rows, all DESeq2 results |
| `figures/08_pathways/gsea_*.csv` | Per-cell-type GSEA results (6 files) |
| `figures/10_cellchat/liana_healthy_vs_severe_all.csv` | 39,530 merged LR pairs |
| `figures/10_cellchat/liana_gained_severe.csv` | Top 16 gained interactions |
| `figures/10_cellchat/liana_lost_severe.csv` | Top 16 lost interactions |
| `figures/06_composition/kruskal_wallis_results.csv` | KW test results, 12 cell types |
| `execution_trace_sess_c0d131b6d4c5.ipynb` | Full analysis notebook (76 cells) |
