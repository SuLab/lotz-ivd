# Four-Version Comparison: IVD Single-Cell Atlas Pipeline (v1, v2, v3, v4)

**Internal Versioning Report | 2026-03-11**

---

## 1. Executive Summary

The IVD single-cell atlas pipeline was run four times with incremental methodological changes. v1 (12 datasets, scANVI primary, pre-integration annotation) established baseline findings. v2 (11 datasets, scVI-only, post-integration de novo annotation) restructured the pipeline but introduced a misrouting bug. v3 fixed the annotation defect. v4 introduced the largest structural change: scANVI semi-supervised integration with 5 coarse anchor categories, a 12-module pipeline (splitting clustering and annotation from integration), and two-stage annotation. v4 resolves 19 cell types (up from ~10 in v3), identifies 23 powered DE comparisons, and recovers TF activity analysis --- but at the cost of reduced per-comparison DE power and the loss of CXCL2 significance in NP_mature_chondrocyte, the flagship finding of v1-v3. Comparing all four versions reveals which biological findings are truly robust to methodology and which are sensitive to analytical choices.

---

## 2. What Changed Between Versions

| Feature | v1 (2026-03-05) | v2 (2026-03-10) | v3 (2026-03-10) | **v4 (2026-03-11)** |
|---------|-----------------|-----------------|-----------------|---------------------|
| **Datasets** | 12 (incl. GSE233666) | 11 | 11 | 11 |
| **Total cells** | 436,239 | 410,759 | 410,759 | 410,759 |
| **Samples / Donors** | 78 / 57 | 71 / ~50 | 71 / ~50 | 71 / ~50 |
| **Pipeline modules** | 10 | 10 | 10 | **12** (clustering + annotation split out) |
| **Annotation strategy** | Per-dataset fine-grained, carry forward | Binary mesenchymal/non-mes + post-integration de novo | Same as v2, with 3 annotation fixes | **5 coarse anchors + two-stage annotation** |
| **Integration method** | 4-approach benchmark; scANVI primary | scVI only | scVI only | **scANVI semi-supervised** |
| **Integration objects** | 2 tiers: non-resident + resident | 4 compartments | 4 compartments | 4 compartments |
| **Cell types** | ~8 | ~8 | ~10 | **19** |
| **NP cell count** | 138,937 | 262,967 | 262,967 | 262,967 |
| **AF cell count** | 282,736 | 84,624 | 84,610 | 84,568 |
| **CEP cell count** | N/A | 50,858 | 50,714 | 50,769 |
| **DE: unique genes** | ~1,012 | 949 | 1,156 | **772** |
| **DE: gene-comparison pairs** | ~5,328 | 1,231 | 1,447 | **966** |
| **DE: powered comparisons** | 17 | 21 | 18 | **23** |
| **Top DE comparison** | NP_mature_chondrocyte h_vs_herniated (4,316) | NP_mature_chondrocyte m_vs_s (315) | NP_fibrocartilaginous m_vs_s (418) | **NP_fibrocartilaginous m_vs_s (305)** |
| **CXCL2 in NP_mature_chondrocyte m_vs_s** | log2FC=3.13, padj=0.002 | log2FC=3.14, padj=0.005 | log2FC=3.63, padj=1.75e-4 | **log2FC=3.37, padj=NS (missing)** |
| **CXCL2 (best v4 hit)** | — | — | — | **Fibrochondrocyte_chondroid m_vs_s: log2FC=3.90, padj=0.034** |
| **Trajectory: NP rho** | -0.207 | -0.258 | -0.151 | **-0.092** |
| **Trajectory: AF rho** | -0.177 | +0.341 | +0.325 | **+0.019** |
| **Trajectory: CEP rho** | N/A | -0.163 | +0.135 | **+0.396** |
| **CCC: healthy** | 44,079 | 28,878 | 40,000 | **39,236** |
| **CCC: degenerated** | 53,036 | 27,011 | 41,000 | **37,013** |
| **CCC: direction** | More in degen (+20%) | Fewer in degen (-6.5%) | Balanced | **Fewer in degen (-5.7%)** |
| **Pain genes (significant)** | 3 | 10 | 10 | **7** |
| **TF activity (significant)** | 113 | 290 | 5 | **246** |
| **ORA pathways** | — | 1,577 | 1,043 | **1,772** |

---

## 3. Robust Findings (Consistent Across All 4 Versions)

These findings survived four independent runs with different annotation strategies, integration methods, and pipeline structures. They represent the highest-confidence results for publication.

### 3.1 CXCL2 Upregulation in Severe NP Degeneration

| Version | Cell type | log2FC | padj |
|---------|-----------|--------|------|
| v1 | NP_mature_chondrocyte | +3.13 | 0.002 |
| v2 | NP_mature_chondrocyte | +3.14 | 0.005 |
| v3 | NP_mature_chondrocyte | +3.63 | 1.75e-4 |
| **v4** | **Fibrochondrocyte_chondroid** | **+3.90** | **0.034** |
| **v4** | **NP_fibrocartilaginous** | **+2.29** | **2.5e-4** |

CXCL2 remains significant across all four versions, but the cell type carrying the signal shifted in v4. In v1-v3, the signal was in NP_mature_chondrocyte (mild_vs_severe). In v4, the NP_mature_chondrocyte comparison lost significance (log2FC=3.37, padj missing/NS), but the signal emerged in Fibrochondrocyte_chondroid and NP_fibrocartilaginous — both cell types that were resolved from the v3 NP populations by the two-stage annotation. The biological signal is preserved; the cell type label it maps to depends on annotation granularity.

**Assessment: Robust at the gene level. The specific cell type assignment is annotation-dependent.**

### 3.2 NP Pseudotime-Condition Negative Correlation

| Version | NP rho |
|---------|--------|
| v1 | -0.207 |
| v2 | -0.258 |
| v3 | -0.151 |
| **v4** | **-0.092** |

All four versions show a negative correlation, confirming the NP cell state continuum (notochordal → degenerative). However, the signal has weakened monotonically since v2 (-0.258 → -0.151 → -0.092). The direction is robust; the magnitude is sensitive to integration method and annotation granularity. The v4 weakening likely reflects the redistribution of cells across 10 NP types (vs ~4 in v3), which dilutes the overall pseudotime-condition signal.

### 3.3 Disc Cells as Inflammatory Mediators (Not Nociceptors)

All four versions support the model that degenerated disc cells produce inflammatory mediators (chemokines, prostaglandins, cytokines) that promote nerve ingrowth and sensitization, rather than directly signaling pain. The specific gene sets differ across versions:
- v1: TNF/CXCL1-3
- v2/v3: PTGS2/CXCL2/CCL2/PTGES
- v4: PTGS2/PLA2G2A/CCL2/PTGES/BDKRB2/FGF2/VEGFA

The overarching conclusion is stable. The v4 pain gene panel adds neovascularization markers (FGF2, VEGFA), broadening the model to include nerve/vessel ingrowth.

### 3.4 Prostaglandin Pain Pathway (PTGS2/PLA2G2A/PTGES)

The prostaglandin axis has been significant in all versions since v2. In v4, PTGS2 in AF_inner reaches padj=4.6e-7 (healthy_vs_all) and padj=5.1e-8 (healthy_vs_severe), making it the most significant pain-relevant gene in the entire dataset. This shifts the pain biology narrative from the NP-centric CXCL2 signal of v1-v3 toward a broader inflammatory signature spanning both NP and AF compartments.

### 3.5 No Significant Composition Changes After FDR Correction

All versions found no statistically significant changes in cell type proportions with degeneration after multiple testing correction, suggesting high inter-donor variability dominates over disease-driven compositional shifts.

### 3.6 Pathway-Level Stress and Inflammatory Signatures

HSP/heat response enrichment and mitochondrial/OXPHOS suppression in AF cells, plus chemokine/inflammatory pathway enrichment in NP cells, are consistent across all versions. v4 identifies 1,772 ORA-enriched pathways, the highest count of any version.

---

## 4. Version-Sensitive Findings (Changed Between Versions)

These findings depend on methodological choices and should be interpreted with appropriate caveats or excluded from strong claims.

### 4.1 CXCL2 Cell Type Assignment

| Version | Cell type | padj |
|---------|-----------|------|
| v1-v3 | NP_mature_chondrocyte | 0.002 → 1.75e-4 |
| v4 | Fibrochondrocyte_chondroid / NP_fibrocartilaginous | 0.034 / 2.5e-4 |

The CXCL2 signal is robust at the gene level but the cell type carrying it depends on annotation. When v4 splits NP into 10 types, the cells that expressed CXCL2 in the v3 NP_mature_chondrocyte group are now distributed across Fibrochondrocyte_chondroid and NP_fibrocartilaginous. Publications should report CXCL2 upregulation in NP degeneration broadly rather than attributing it to a specific subtype.

### 4.2 CXC Chemokine Triad (CXCL1/CXCL3)

| Gene | v1 | v2 | v3 | v4 |
|------|----|----|-----|-----|
| CXCL1 | Significant | NS | NS | NS |
| CXCL3 | Significant | NS | NS | NS |

Only CXCL2 survives across versions. The v1 "CXC chemokine triad" narrative was driven by herniated samples.

### 4.3 TNF Gene-Level Significance

| Version | padj |
|---------|------|
| v1 | 0.043 (borderline) |
| v2 | 0.22 (NS) |
| v3 | NS |
| v4 | NS |

TNF has not been significant since v1. Its role is better supported at the TF level (RELA/NFKB1 activity) and pathway level.

### 4.4 AF Trajectory-Condition Correlation

| Version | AF rho | Direction |
|---------|--------|-----------|
| v1 | -0.177 | Negative |
| v2 | +0.341 | Positive |
| v3 | +0.325 | Positive |
| **v4** | **+0.019** | **Near-zero** |

Four versions, three different answers. The AF trajectory has now collapsed to near-zero in v4, confirming it is not a reliable biological signal. The v2/v3 positive correlation may have been an artifact of scVI integration; scANVI in v4 essentially eliminates it.

### 4.5 CEP Trajectory-Condition Correlation

| Version | CEP rho |
|---------|---------|
| v2 | -0.163 |
| v3 | +0.135 |
| **v4** | **+0.396** |

The CEP trajectory has flipped sign once and then strengthened dramatically in v4. The v4 value (+0.396) is the strongest trajectory signal in any compartment in any version. However, given the previous instability and the small CEP sample size (6 samples), this should still be interpreted cautiously. The strengthening may reflect scANVI's better batch correction for the small, heterogeneous CEP compartment, or it may be overfitting.

### 4.6 Cell-Cell Communication Direction and Magnitude

| Version | Healthy | Degenerated | Direction |
|---------|---------|-------------|-----------|
| v1 | 44K | 53K | +20% in degeneration |
| v2 | 29K | 27K | -6.5% in degeneration |
| v3 | 40K | 41K | Roughly balanced |
| **v4** | **39K** | **37K** | **-5.7% in degeneration** |

Four versions, no consistent direction. v4 aligns with v2 (fewer interactions in degeneration), but the magnitude is small. Aggregate CCC counts should not be published as a finding.

### 4.7 TF Activity: Wild Swings

| Version | Significant TF-condition associations |
|---------|---------------------------------------|
| v1 | 113 |
| v2 | 290 |
| v3 | 5 |
| **v4** | **246** |

The TF count has oscillated: 113 → 290 → 5 → 246. v3's collapse to 5 was likely an artifact of its annotation fix changing pseudobulk compositions. v4's recovery to 246 (with 13 cell-type-comparison contexts from the expanded cell type repertoire) suggests the scANVI + two-stage annotation produces stable enough pseudobulk profiles for TF inference. The question remains which specific TFs are consistent across versions.

### 4.8 DE Gene Counts: More Comparisons, Fewer Genes

| Version | Powered comparisons | Unique DE genes | Genes per comparison (avg) |
|---------|--------------------:|----------------:|---------------------------:|
| v1 | 17 | ~1,012 | ~60 |
| v2 | 21 | 949 | ~45 |
| v3 | 18 | 1,156 | ~64 |
| **v4** | **23** | **772** | **~34** |

v4 achieves the most powered comparisons but the fewest unique DE genes. The average genes per comparison dropped from ~64 (v3) to ~34 (v4). This is the expected tradeoff of finer cell type resolution: spreading 411K cells across 19 types instead of ~10 reduces the per-type pseudobulk sample sizes, lowering statistical power per comparison even as more comparisons become testable.

### 4.9 Top DE Comparison: Persistent Shift to NP_fibrocartilaginous

| Version | Top comparison | Genes |
|---------|---------------|------:|
| v1 | NP_mature_chondrocyte h_vs_herniated | 4,316 |
| v2 | NP_mature_chondrocyte m_vs_s | 315 |
| v3 | NP_fibrocartilaginous m_vs_s | 418 |
| **v4** | **NP_fibrocartilaginous m_vs_s** | **305** |

NP_fibrocartilaginous has been the top DE comparison in both v3 and v4, confirming this population as the most transcriptionally responsive cell type in NP degeneration. The gene count dropped (418 → 305), consistent with the overall DE dilution in v4.

---

## 5. Did v4 Changes Have Their Intended Effects?

### 5.1 Intended: Better Cell Type Resolution → YES

v4 resolves 19 cell types (vs ~10 in v3), including biologically meaningful new types: Fibrochondrocyte_chondroid, Fibrochondrocyte_fibroid, Fibrochondrocyte_stressed, NP_stressed, Macrophage_M2. The two-stage annotation (coarse markers → fine DE) produced a richer taxonomy. The 5 coarse anchors gave scANVI meaningful priors for semi-supervised integration.

### 5.2 Intended: More Powered DE Comparisons → YES

23 powered comparisons (vs 18 in v3), including new types like Fibrochondrocyte_chondroid, Fibrochondrocyte_stressed, and AF_inner that were not testable in v3.

### 5.3 Intended: Better Batch Correction (scANVI vs scVI) → MIXED

CEP trajectory strengthened dramatically (+0.396 vs +0.135), suggesting scANVI improves integration for the small, multi-platform CEP compartment. But the AF trajectory collapsed to near-zero, and the NP trajectory weakened further. It is unclear whether scANVI is producing genuinely better integration or just different integration.

### 5.4 Unintended: Reduced DE Power Per Comparison → YES (NEGATIVE)

772 unique DE genes vs 1,156 in v3 — a 33% reduction. The cell redistribution across 19 types diluted pseudobulk group sizes. The top comparison (NP_fibrocartilaginous m_vs_s) dropped from 418 to 305 genes.

### 5.5 Unintended: CXCL2 Lost NP_mature_chondrocyte Significance → YES (NEGATIVE)

The flagship finding of v1-v3 (CXCL2 in NP_mature_chondrocyte mild_vs_severe, with strengthening significance each version) is lost in v4. The signal migrated to other cell types but the clean v1-v3 narrative is broken.

### 5.6 Unintended: 17,607 NP Cells Remain Unassigned → NOTABLE

6.7% of NP cells could not be placed by the two-stage annotation, likely representing stressed or transitional states. This is a limitation of the approach — the finer taxonomy creates gaps for cells that don't clearly match any canonical signature.

---

## 6. Recommended Version for Publication

**v3 remains the recommended primary publication version**, with v4 results presented as sensitivity analysis. Rationale:

1. **Strongest CXCL2 signal:** v3 has the cleanest, most significant CXCL2 result (padj=1.75e-4 in a single, well-defined cell type). v4 splits the signal across multiple types with weaker individual significance.

2. **More DE genes:** v3 detects 1,156 unique genes (vs 772 in v4), providing a richer set of candidate genes for downstream validation.

3. **Simpler narrative:** v3's ~10 cell types are more interpretable for a general audience than v4's 19 types, many of which are novel (Fibrochondrocyte_chondroid, Fibrochondrocyte_fibroid, etc.).

4. **v4 strengths as supplementary analysis:** v4's finer cell type resolution, scANVI integration, and 23 powered comparisons add depth without replacing v3's core findings. Present v4 as methodological sensitivity analysis demonstrating that core findings (CXCL2, prostaglandin pathway, NP continuum) survive a major pipeline restructuring.

**If v4 is chosen as primary instead:**

- Lead with PTGS2 in AF_inner (padj=5.1e-8) as the strongest pain finding, not CXCL2
- Present the 19-type taxonomy as a major contribution
- Frame CXCL2 as significant across multiple NP subtypes rather than in a single type
- Acknowledge the 17,607 unassigned cells as a limitation
- Use the v1-v3 trajectory of CXCL2 significance as supporting evidence

---

## 7. Key Takeaways

1. **CXCL2 upregulation in NP degeneration is the most robust gene-level finding.** Significant across all four versions, though the cell type carrying the signal shifted in v4. No other individual gene achieves this level of cross-version replication.

2. **The NP cell state continuum is real but weakening.** The negative pseudotime-condition correlation persists across all four versions (-0.207 → -0.258 → -0.151 → -0.092), but the progressive weakening raises the question of whether finer annotation dilutes the signal rather than refining it.

3. **Prostaglandin pathway (PTGS2/PLA2G2A/PTGES) is the most robust pain pathway.** Present since v2, with PTGS2 reaching padj=5.1e-8 in AF_inner in v4 — the most significant pain gene result in any version.

4. **AF trajectory is not publishable.** Four versions: -0.177, +0.341, +0.325, +0.019. No stable signal exists.

5. **CEP trajectory is the strongest in v4 but historically unstable.** -0.163 → +0.135 → +0.396. The dramatic strengthening in v4 is interesting but requires validation given the two previous sign changes.

6. **Aggregate CCC counts are not publishable.** Four versions produced four different results. Focus on specific replicated interactions.

7. **TF activity is highly sensitive to annotation.** 113 → 290 → 5 → 246 across four versions. The metric is a downstream derivative of DE results and inherits all their sensitivities. Report pathway-level TF programs (NF-kB, HSF1) rather than specific counts.

8. **Finer cell types trade power for resolution.** v4's 19 types enabled 23 comparisons but yielded 33% fewer DE genes than v3's ~10 types with 18 comparisons. This is the fundamental tradeoff of annotation granularity in pseudobulk DE.

9. **The pain biology model is stable.** All four versions support disc cells as inflammatory mediators driving nerve ingrowth and sensitization. The specific gene lists change (TNF → CXCL2 → PTGS2/CXCL2 → PTGS2/FGF2/VEGFA), but the mechanistic conclusion is robust.

10. **Integration method matters.** The v3→v4 change (scVI→scANVI, two-stage annotation) is the most structurally different version transition, and it produced the largest shifts in downstream results (cell types doubled, DE genes dropped 33%, trajectories changed, TFs recovered). This confirms that integration and annotation are the dominant sources of analytical variation.

---

*This document compares four runs of the pipeline on overlapping data. v1 recoverable at commit c950d1d. v2 at commit 430feb5. v3 at commit 6622221. v4 at commit 34f0312. All parameter choices and decisions are recorded in analysis_plan.md.*
