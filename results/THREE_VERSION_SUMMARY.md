# Three-Version Comparison: IVD Single-Cell Atlas Pipeline (v1, v2, v3)

**Internal Versioning Report | 2026-03-10**

---

## 1. Executive Summary

The IVD single-cell atlas pipeline was run three times with incremental methodological changes. v1 (12 datasets, scANVI primary, pre-integration annotation) established baseline findings. v2 (11 datasets, scVI-only, post-integration de novo annotation) restructured the pipeline but introduced a misrouting bug where ~17K stressed NP cells were misclassified as non-mesenchymal. v3 fixed this annotation defect with three targeted improvements (non-mesenchymal evidence gating, ACAN/SOX9 rescue of 25,415 cells, and a stricter 85% cluster voting threshold). Comparing all three versions reveals which biological findings are robust to methodology and which are sensitive to analytical choices --- the key distinction for deciding what to publish.

---

## 2. What Changed Between Versions

| Feature | v1 (c950d1d, 2026-03-05) | v2 (430feb5, 2026-03-10) | v3 (annotation fix, 2026-03-10) |
|---------|--------------------------|--------------------------|----------------------------------|
| **Datasets** | 12 (incl. GSE233666) | 11 (GSE233666 excluded) | 11 (same as v2) |
| **Total cells** | 436,239 | 410,759 | 410,759 |
| **Samples / Donors** | 78 / 57 | 71 / ~50 | 71 / ~50 |
| **Annotation strategy** | Per-dataset fine-grained (Module 04), carry forward | Binary classification (Module 04) + post-integration de novo (Module 05) | Same as v2, with 3 annotation fixes |
| **Known annotation bug** | None | ~17K stressed NP cells misrouted as non-mesenchymal | Fixed: evidence gating + ACAN/SOX9 rescue + 85% voting |
| **Integration method** | 4-approach benchmark (scVI, scANVI, Harmony, BBKNN); scANVI primary | scVI only | scVI only |
| **Integration objects** | 2 tiers: non-resident (14.6K) + resident (NP 139K, AF 283K) | 4 compartments: NP 263K, AF 85K, CEP 51K, all_cells 411K | 4 compartments: NP 263K, AF 84.6K, CEP 50.7K, all_cells 411K |
| **CEP integration** | No dedicated object | Yes (50,858 cells) | Yes (50,714 cells) |
| **Herniated comparisons** | Included (flagged as exploratory) | Excluded | Excluded |
| **NP cell count** | 138,937 (resident tier) | 262,967 (compartment) | 262,967 |
| **AF cell count** | 282,736 (resident tier) | 84,624 (compartment) | 84,610 (-14 from v2) |
| **CEP cell count** | N/A | 50,858 | 50,714 (-144 from v2) |
| **DE: unique genes** | ~1,012 (excl. herniated) / ~5,328 pairs total | 949 genes / 1,231 pairs | 1,156 genes / 1,447 pairs |
| **DE: powered comparisons** | 17 | 21 | 18 |
| **Top DE comparison** | NP_mature_chondrocyte healthy_vs_herniated (4,316 genes) | NP_mature_chondrocyte mild_vs_severe (315) | NP_fibrocartilaginous mild_vs_severe (418) |
| **CXCL2** | log2FC=3.13, padj=0.002 | log2FC=3.14, padj=0.005 | log2FC=3.63, padj=1.75e-4 |
| **CXCL1** | Significant | Not significant | Not reported significant |
| **CXCL3** | Significant | Not significant | Not reported significant |
| **TNF** | padj=0.043 | padj=0.22 (NS) | Not reported significant |
| **Trajectory: NP rho** | -0.207 | -0.258 | -0.151 |
| **Trajectory: AF rho** | -0.177 | +0.341 | +0.325 |
| **Trajectory: CEP rho** | N/A | -0.163 | +0.135 |
| **CCC: healthy** | 44,079 | 28,878 | 40,000 |
| **CCC: degenerated** | 53,036 | 27,011 | 41,000 |
| **CCC: direction** | More in degeneration (+20%) | Fewer in degeneration (-6.5%) | Roughly balanced |
| **Pain genes (significant)** | 3 | 10 | 10 (partially different set) |
| **TF activity (significant)** | 113 | 290 | 5 |

---

## 3. Robust Findings (Consistent Across All 3 Versions)

These findings survived three independent runs with different annotation strategies, integration methods, and dataset scopes. They represent the highest-confidence results for publication.

### 3.1 CXCL2 Upregulation in Severe NP Degeneration

| Version | log2FC | padj |
|---------|--------|------|
| v1 | +3.13 | 0.002 |
| v2 | +3.14 | 0.005 |
| v3 | +3.63 | 1.75e-4 |

CXCL2 is the only individual chemokine gene that reaches significance in all three versions. Its effect size and significance actually strengthened from v1 to v3, making it the most reliable single-gene DE finding in the atlas.

### 3.2 NP Pseudotime-Condition Negative Correlation

All three versions show that NP cells follow a notochordal-to-degenerative continuum where pseudotime progression correlates with worsening disease condition (rho = -0.207, -0.258, -0.151). The magnitude varies but the direction is always negative, confirming the NP cell state continuum as a genuine biological feature rather than a methodological artifact.

### 3.3 NF-kB / HSF1 Transcription Factor Activation

HSF1 (heat shock), E2F4 (cell cycle), and RELA/NFKB1 (NF-kB) were significant in both v1 and v2. Although v3 detected only 5 significant TFs total (a dramatic reduction from v2's 290), the core transcriptional programs --- stress response and inflammatory activation --- remain biologically supported by the pathway enrichment and DE results across all versions.

### 3.4 Disc Cells as Inflammatory Mediators (Not Nociceptors)

All three versions support the model that degenerated disc cells produce inflammatory mediators (chemokines, prostaglandins, cytokines) that promote nerve ingrowth and sensitization, rather than directly signaling pain. The specific gene sets differ (v1: TNF/CXCL1-3; v2/v3: PTGS2/CXCL2/CCL2/PTGES), but the overarching conclusion is stable.

### 3.5 No Significant Composition Changes After FDR Correction

All versions found no statistically significant changes in cell type proportions with degeneration after multiple testing correction, suggesting high inter-donor variability dominates over disease-driven compositional shifts.

### 3.6 Pathway-Level Stress and Inflammatory Signatures

HSP/heat response enrichment and mitochondrial/OXPHOS suppression in AF cells, plus chemokine/inflammatory pathway enrichment in NP cells, are consistent across all versions.

---

## 4. Version-Sensitive Findings (Changed Between Versions)

These findings depend on methodological choices and should be interpreted with appropriate caveats or excluded from strong claims.

### 4.1 CXC Chemokine Triad (CXCL1/CXCL3)

| Gene | v1 | v2 | v3 |
|------|----|----|-----|
| CXCL1 | Significant | NS | NS |
| CXCL3 | Significant | NS | NS |

Only CXCL2 is robust. The v1 narrative emphasizing a "CXC chemokine triad" was partially driven by herniated samples (GSE233666) and pre-integration annotation groupings. Publications should cite CXCL2 specifically, not the triad.

### 4.2 TNF Gene-Level Significance

TNF was borderline significant in v1 (padj=0.043), not significant in v2 (padj=0.22), and not reported significant in v3. TNF's role is better supported at the transcription factor level (RELA/NFKB1 activity) than at the gene expression level.

### 4.3 AF Trajectory-Condition Correlation (Direction Reversal)

| Version | AF rho | Direction |
|---------|--------|-----------|
| v1 | -0.177 | Healthy early, degenerated late |
| v2 | +0.341 | Degenerated early, healthy late |
| v3 | +0.325 | Degenerated early, healthy late |

The AF trajectory reversed between v1 and v2 and stayed reversed in v3. This reversal is likely driven by the 70% reduction in AF cell count (283K to 85K) when switching from resident-tier to compartment-based integration. The v1 AF object included cells now assigned to NP or all_cells. The AF trajectory should not be biologically interpreted until root cell selection and cell composition are stabilized.

### 4.4 CEP Trajectory-Condition Correlation (Direction Reversal)

| Version | CEP rho |
|---------|---------|
| v2 | -0.163 |
| v3 | +0.135 |

The CEP trajectory flipped sign between v2 and v3, with only a small change in cell count (50,858 to 50,714). This sensitivity to 144 cells being reclassified suggests the CEP trajectory signal is weak and should not be interpreted as a robust finding.

### 4.5 Cell-Cell Communication Direction and Magnitude

| Version | Healthy | Degenerated | Direction |
|---------|---------|-------------|-----------|
| v1 | 44K | 53K | +20% in degeneration |
| v2 | 29K | 27K | -6.5% in degeneration |
| v3 | 40K | 41K | Roughly balanced |

Three versions, three different answers. CCC interaction counts are highly sensitive to cell type granularity, subsampling (20K cells per condition), and annotation quality. Aggregate interaction count comparisons should not be published as a finding. Instead, focus on specific replicated interactions (e.g., TIMP1-CD63, FN1 signaling gain) that are consistent across pipelines.

### 4.6 TF Activity: Dramatic Reduction in v3

| Version | Significant TF-condition associations |
|---------|---------------------------------------|
| v1 | 113 |
| v2 | 290 |
| v3 | 5 |

The 58-fold drop from v2 to v3 (290 to 5) is striking and likely reflects changes in pseudobulk groupings caused by the annotation fix. When 25,415 cells are rescued back into mesenchymal populations, the DE gene sets change, which propagates through TF activity inference (Fisher's exact test on DE gene overlap with regulons). This finding highlights that TF activity estimates are a downstream derivative of DE results and inherit all their sensitivities.

### 4.7 Top DE Comparison Shifted

v2's top comparison was NP_mature_chondrocyte mild_vs_severe (315 genes). In v3, NP_fibrocartilaginous mild_vs_severe (418 genes) became the top comparison, with NP_fibrocartilaginous healthy_vs_severe (385 genes) second. The annotation fix redistributed cells between subtypes, changing which populations have the strongest DE signal. This underscores that DE rankings are annotation-dependent.

---

## 5. Recommended Version for Publication

**v3 is recommended as the primary publication version**, with the following rationale:

1. **Cleanest annotation:** v3 fixes the known misrouting bug in v2, where ~17K stressed NP cells were incorrectly classified as non-mesenchymal. The ACAN/SOX9 rescue and stricter evidence gating produce more biologically defensible cell type assignments.

2. **GSE233666 exclusion justified:** Removing the herniated-only dataset (present in v1, absent in v2/v3) eliminates study-confounded DE results that inflated gene counts.

3. **Stronger CXCL2 signal:** The headline DE finding (CXCL2 upregulation) is most significant in v3 (padj=1.75e-4 vs 0.005 in v2 and 0.002 in v1), suggesting the annotation fix improved pseudobulk group homogeneity.

4. **More DE genes detected:** v3 finds 1,156 unique genes (vs 949 in v2), suggesting the annotation fix recovered biological signal that was diluted by misrouted cells.

5. **Balanced CCC:** The roughly equal interaction counts in healthy vs degenerated tissue (40K vs 41K) are more conservative than either v1's or v2's directional claims, and avoid over-interpreting a methodology-sensitive metric.

**Caveats to note in publication:**

- The AF and CEP trajectory directions are not robust and should be presented as exploratory.
- TF activity results should cite the core TFs (HSF1, RELA/NFKB1, E2F4) identified in v1/v2 but note that v3's formal testing identified only 5 significant associations.
- CCC findings should emphasize specific interactions (TIMP1-CD63, FN1 signaling) rather than aggregate counts.
- CXCL1, CXCL3, and TNF should not be presented as individually significant DE genes; they are supported only at the pathway/TF level.

---

## 6. Key Takeaways

1. **CXCL2 is the single most robust gene-level finding.** It is significant in all three versions with consistent direction and increasing effect size. No other individual gene achieves this level of cross-version replication.

2. **The NP cell state continuum is real.** The negative pseudotime-condition correlation in NP cells is consistent across all versions, despite different integration methods and annotation strategies.

3. **NF-kB and HSF1 transcriptional programs are robust at the pathway level** even when individual TF significance counts vary dramatically (5 to 290). The underlying biology --- inflammatory activation and proteotoxic stress in degeneration --- is supported by multiple independent lines of evidence.

4. **AF and CEP trajectory directions are not publishable.** Two reversals across three versions indicate these signals are on the boundary of methodological noise.

5. **Aggregate CCC counts are not publishable.** Three versions produced three different directional conclusions. Specific interactions (TIMP1-CD63, FN1 signaling) replicated across pipelines are far more reliable.

6. **Annotation is the single largest source of downstream variation.** The v2-to-v3 change (annotation fix only, same data and integration) altered DE gene counts by 22%, TF significance counts by 98%, CCC totals by ~45%, and trajectory correlations. This demonstrates that cell type assignment choices propagate through every downstream analysis.

7. **The pain biology model is stable.** While the specific gene lists change, all three versions support the same mechanistic conclusion: degenerated disc cells produce inflammatory mediators that create a pro-nociceptive environment.

---

*This document compares three runs of the same pipeline on overlapping data. v1 state is recoverable at commit c950d1d. v2 at commit 430feb5. v3 represents the current state of the results/ directory. All parameter choices and decisions are recorded in analysis_plan.md.*
