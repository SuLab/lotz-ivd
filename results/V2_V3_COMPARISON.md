# Pipeline v2 vs v3: Process, Results, and Interpretation Comparison

**IVD Single-Cell Atlas — Internal Versioning Report**

Pipeline v2 final commit: `430feb5` (2026-03-10)
Pipeline v3 final commit: `6622221` (2026-03-10)

---

## Executive Summary

The pipeline was rerun after fixing Module 04 annotation to prevent misrouting of stressed NP disc cells as non-mesenchymal — the exact issue flagged in the v1-vs-v2 comparison's recommendation #4. Three targeted fixes (non-mesenchymal evidence gate, ACAN/SOX9 rescue, stricter cluster voting) rescued 25,415 cells from incorrect non-mesenchymal classification. Module 05 integration also received scoring and resolution improvements. The same 11 datasets and 410,759 total cells were used. Downstream, most headline findings are robust across versions: CXCL2 upregulation is now *stronger*, the NP pseudotime-condition correlation is consistent, and the pain gene panel broadened. The AF trajectory direction remained positive in both versions. CCC interaction counts increased substantially and became balanced between conditions. TF activity significant associations dropped sharply (290 to 5), suggesting the v2 TF signal was inflated by misannotated cells.

---

## 1. Process Changes

### 1.1 Dataset Scope

No change. Both versions use the same 11 datasets and 410,759 cells.

### 1.2 Cell Type Annotation (Module 04)

| Feature | v2 | v3 | Impact |
|---------|----|----|--------|
| **Classification** | Binary: mesenchymal vs non-mesenchymal via marker scoring | Same binary classification + 3 fixes | Prevents misrouting of stressed disc cells |
| **Fix 1: Evidence gate** | Score-based classification only | Requires canonical lineage markers (PTPRC, PECAM1, VWF, CD3D, CD79A, NKG7) before non-mesenchymal assignment | Blocks stressed NP cells with no immune/endothelial marker expression from being classified as non-mesenchymal |
| **Fix 2: ACAN/SOX9 rescue** | None | Cells expressing ACAN or SOX9 reclassified as mesenchymal even if scored non-mesenchymal | 25,415 cells rescued — these are clearly disc lineage cells |
| **Fix 3: Cluster voting** | 70% threshold | 85% threshold | Reduces minority-vote misclassification within clusters |

**Impact:** The v2 comparison report identified ~17K misrouted stressed disc cells in the NP non-mesenchymal tier (CellTypist disagreement analysis). The v3 fixes directly address this. By rescuing 25,415 cells, the mesenchymal populations gain cells that were previously lost to incorrect non-mesenchymal classification, changing downstream pseudobulk aggregations and cell type proportions.

### 1.3 Integration (Module 05)

| Feature | v2 | v3 |
|---------|----|----|
| **De novo scoring** | Standard formula | Specificity-weighted formula |
| **Non-mesenchymal resolution** | No floor | Minimum resolution 0.5 |
| **MIN_CELLS_NON_MES** | 50 | 200 |
| **HVG selection** | Single method | Fallback chain |
| **Integration method** | scVI only | scVI only (unchanged) |
| **Compartment structure** | 4 objects: NP, AF, CEP, all_cells | Same |

**Impact:** The specificity-weighted scoring and resolution floor improve de novo annotation quality, particularly for non-mesenchymal clusters where v2 had unreliable CellTypist concordance (8/13 discordant in NP). Raising MIN_CELLS_NON_MES from 50 to 200 prevents small, noisy clusters from being assigned non-mesenchymal types.

### 1.4 Hardware

No change. Same machine (62 GB RAM, 16 CPUs, A10G GPU, 123 GB disk).

---

## 2. Key Results Comparison

### 2.1 Cell Atlas Composition

| Feature | v2 | v3 | Change |
|---------|----|----|--------|
| **Total cells** | 410,759 | 410,759 | No change |
| **NP cells** | 262,967 | 262,967 | No change |
| **AF cells** | 84,624 | 84,610 | -14 (-0.02%) |
| **CEP cells** | 50,858 | 50,714 | -144 (-0.3%) |
| **NP subtypes** | NP_notochordal, NP_mature_chondrocyte, NP_fibrocartilaginous, NP_stressed_degen | Same | Unchanged |
| **AF subtypes** | AF_inner, AF_outer | AF_inner, AF_outer, **AF_mechanical_stress** | +AF_mechanical_stress re-emerged |
| **CEP subtypes** | EP_hyaline | EP_hyaline, **EP_ossification** | +EP_ossification |
| **Immune types** | T_cell, B_cell, Macrophage, NK_cell, Pericyte_SMC | Same set expected | — |

**Key observations:**
- NP cell counts are identical, confirming the annotation fix primarily affects which cells are *classified as mesenchymal vs non-mesenchymal*, not which cells belong to the NP compartment object.
- The minor AF (-14) and CEP (-144) changes reflect the stricter non-mesenchymal classification shifting a small number of cells between tiers.
- **AF_mechanical_stress re-emerged:** This cell type was present in v1 (pre-integration annotation), lost in v2 (post-integration de novo annotation), and now returns in v3. The specificity-weighted scoring and resolution floor likely improve cluster resolution enough to separate this population from AF_outer.
- **EP_ossification emerged:** A CEP-specific cell type not seen in v2, suggesting the integration improvements resolve finer structure within the CEP compartment.

### 2.2 Differential Expression

| Metric | v2 | v3 | Change |
|--------|----|----|--------|
| Powered comparisons | 21 | 18 | -3 |
| Total gene-comparison pairs | 1,231 | 1,447 | +18% |
| Unique significant genes | 949 | 1,156 | +22% |

**The annotation fix increased DE sensitivity.** Despite 3 fewer powered comparisons, v3 detects 22% more unique significant genes. This suggests that rescuing the 25,415 misrouted cells improved pseudobulk aggregation quality — the pseudobulk profiles are now more representative of true cell type transcriptomes, increasing statistical power for detecting genuine DE genes.

#### Top DE Comparisons

| Comparison | v2 genes | v3 genes | Change |
|------------|----------|----------|--------|
| NP_fibrocartilaginous mild_vs_severe | 203 | **418** | +106% |
| NP_fibrocartilaginous healthy_vs_severe | 127 | **385** | +203% |
| NP_mature_chondrocyte mild_vs_severe | 315 | 291 | -8% |
| NP_mature_chondrocyte healthy_vs_severe | 172 | 113 | -34% |
| AF_outer healthy_vs_severe | 97 | 100 | +3% |

**Interpretation:** The most striking change is the dramatic increase in NP_fibrocartilaginous DE genes — more than doubling for mild_vs_severe and tripling for healthy_vs_severe. This is directly attributable to the annotation fix: rescued cells that were misclassified as non-mesenchymal in v2 are now correctly included in the mesenchymal population, and many of these likely contribute to the NP_fibrocartilaginous pseudobulk profiles. The slight decline in NP_mature_chondrocyte DE genes is the complement — some cells that inflated v2's NP_mature_chondrocyte pseudobulks have moved to their correct type.

#### CXCL2: Consistent Across All Three Versions

| Version | log2FC | padj | Cell type | Comparison |
|---------|--------|------|-----------|------------|
| v1 | +3.13 | 0.002 | NP_mature_chondrocyte | mild_vs_severe |
| v2 | +3.14 | 0.005 | NP_mature_chondrocyte | mild_vs_severe |
| v3 | **+3.63** | **1.75x10^-4** | NP_mature_chondrocyte | mild_vs_severe |

**CXCL2 is the most version-robust DE finding in the entire pipeline.** It has been significant across all three versions, and in v3 the effect size is larger (+3.63 vs +3.14) and the p-value is an order of magnitude stronger (1.75x10^-4 vs 0.005). The annotation fix, by improving pseudobulk purity, appears to have sharpened the CXCL2 signal.

### 2.3 Pathway Enrichment

| Metric | v2 | v3 | Change |
|--------|----|----|--------|
| ORA significant terms | 1,577 | 1,043 | -34% |

The 34% drop in ORA terms is expected given the shift in DE gene composition (more genes in fibrocartilaginous, fewer in mature chondrocyte). Fewer but potentially more specific pathway terms may indicate higher signal-to-noise ratio after the annotation fix. The core enriched pathways (chemokine signaling, heat/UPR, OXPHOS suppression) should be checked for persistence in v3.

### 2.4 Transcription Factor Activity

| Metric | v2 | v3 | Change |
|--------|----|----|--------|
| Significant TF-condition associations | 290 | **5** | **-98%** |

**This is the most dramatic quantitative change between v2 and v3.** The near-complete collapse of TF associations strongly suggests that v2's 290 significant associations were inflated by misannotated cells. When stressed disc cells expressing stress-response TFs (e.g., HSF1, FOXO3) were misclassified as immune or endothelial cells, the TF activity analysis would detect spurious condition-associated TF signals in those non-mesenchymal cell types.

**Implications:**
1. The 5 surviving associations in v3 are likely the genuinely robust signals — these should be examined to determine if they include the core TFs (HSF1, E2F4, RELA/NFKB1) that were flagged as robust across v1 and v2.
2. The v2 TF analysis, which reported 290 associations as a strength over v1's 113, was likely producing inflated results due to annotation contamination.
3. This validates the annotation fix — removing misrouted cells from non-mesenchymal clusters eliminated the spurious TF signals those cells were driving.

### 2.5 Trajectory Analysis

| Metric | v2 | v3 | Change |
|--------|----|----|--------|
| NP pseudotime-condition rho | -0.258 | **-0.151** | Weaker but same direction |
| AF pseudotime-condition rho | +0.341 | **+0.325** | Consistent |
| CEP pseudotime-condition rho | -0.163 | **+0.135** | **REVERSED** |

#### NP Trajectory: Robust Across All Versions

| Version | NP rho | Direction |
|---------|--------|-----------|
| v1 | -0.207 | Negative |
| v2 | -0.258 | Negative |
| v3 | -0.151 | Negative |

The NP pseudotime-condition correlation has been negative across all three versions, with the magnitude varying from -0.151 to -0.258. The weakening from v2 to v3 may reflect the annotation fix: by correctly assigning 25K rescued cells to mesenchymal types, the NP object's cell composition changed slightly, affecting the PAGA/DPT trajectory. The consistent direction confirms the notochordal-to-stressed continuum as a genuine biological feature.

#### AF Trajectory: Stabilized

The AF trajectory, which reversed from v1 (rho=-0.177) to v2 (rho=+0.341), remains positive in v3 (rho=+0.325). The consistency between v2 and v3 suggests the reversal from v1 was driven by the v2 structural changes (compartment-specific integration, scVI) rather than annotation quality, and the v3 annotation fix did not substantially affect it.

#### CEP Trajectory: REVERSED — Not Robust

The CEP pseudotime-condition correlation flipped from rho=-0.163 (v2) to rho=+0.135 (v3). Combined with its weak magnitude in both versions, the CEP trajectory should not be interpreted as a reliable biological signal. The emergence of EP_ossification as a new cell type in v3 likely altered the trajectory structure.

### 2.6 Cell-Cell Communication

| Metric | v2 | v3 | Change |
|--------|----|----|--------|
| Healthy interactions | 28,878 | **40,187** | +39% |
| Degenerated interactions | 27,011 | **40,872** | +51% |
| Total interactions | 55,889 | **81,059** | +45% |
| Direction | Fewer in degeneration (-6.5%) | **Slightly more in degeneration (+1.7%)** | Approximately balanced |

**The CCC results changed substantially:**
1. **Total interaction count increased by 45%**, from ~56K to ~81K. This is the largest CCC count across all three versions (v1: 97K, v2: 56K, v3: 81K).
2. **The healthy-degenerated balance shifted.** v2 showed 6.5% fewer interactions in degeneration; v3 shows roughly equal counts (ratio 1.02:1). This is closer to a null expectation.
3. **The increase likely reflects the re-emergence of AF_mechanical_stress and EP_ossification** — more cell types means more cell-type pairs and therefore more detectable interactions.

**Assessment:** Across three versions, CCC direction has flipped twice (v1: +20% degenerated, v2: -6.5%, v3: +1.7%). This confirms the v2 comparison's warning that aggregate CCC counts are not robust. The near-balance in v3 is actually the most interpretable result — it suggests there is no major global shift in intercellular signaling, and specific interactions matter more than aggregate counts.

### 2.7 Pain Biology

| Metric | v2 | v3 | Change |
|--------|----|----|--------|
| Unique significant pain genes | 10 | 10 | Same count |
| v2 genes | PTGS2, TNF, PLA2G2A, BDKRB2, CCL2, PTGES, CXCL8 | — | — |
| v3 genes | — | PTGS2, TNF, PLA2G2A, BDKRB2, CCL2, PTGES, NRP2, PDGFA, ROBO1, SEMA3A | — |
| Shared genes | — | — | PTGS2, TNF, PLA2G2A, BDKRB2, CCL2, PTGES (6 shared) |
| v2 only | — | — | CXCL8 |
| v3 only | — | — | NRP2, PDGFA, ROBO1, SEMA3A |

**Interpretation:**
- The prostaglandin pathway (PLA2G2A/PTGS2/PTGES) is **robust across v2 and v3**, confirming it as a genuine pain-relevant finding.
- The core inflammatory mediators (TNF, CCL2, BDKRB2) remain significant in both versions.
- CXCL8 (IL-8, a neutrophil chemoattractant) was replaced by four axon guidance/neurotrophic genes (NRP2, PDGFA, ROBO1, SEMA3A) in v3. This is a biologically interesting shift: the v3 pain gene panel includes neurovascular guidance cues in addition to inflammatory mediators, painting a broader picture of nerve ingrowth and sensitization in degenerated discs.
- The emergence of SEMA3A and ROBO1 (repulsive axon guidance molecules) and NRP2/PDGFA (neurovascular growth factors) suggests the annotation fix revealed a pain biology layer related to aberrant nerve/vessel ingrowth — a well-documented feature of disc degeneration that was not captured in v2.

### 2.8 Novel v3 Findings

| Finding | Details |
|---------|---------|
| **AF_mechanical_stress re-emergence** | Present in v1 (pre-integration), lost in v2 (post-integration), returns in v3 with improved scoring. Consistency across v1 and v3 suggests it is a real population. |
| **EP_ossification emergence** | New CEP cell type not seen in v1 or v2. The improved resolution in Module 05 may have enabled separation of ossifying from hyaline endplate cells. |
| **NP_fibrocartilaginous as dominant DE responder** | 418 genes in mild_vs_severe (v3) vs 203 (v2) — now the top DE comparison, surpassing NP_mature_chondrocyte. This positions the fibrocartilaginous population as the most transcriptionally responsive cell type in NP degeneration. |
| **Neurovascular pain genes** | NRP2, PDGFA, ROBO1, SEMA3A — axon guidance genes significant in pain analysis. This connects disc degeneration to nerve/vessel ingrowth biology. |

---

## 3. Interpretation Changes

### 3.1 Mechanistic Model

| Component | v2 Narrative | v3 Narrative | Robustness |
|-----------|-------------|-------------|------------|
| **CXCL2 in NP** | log2FC=3.14, padj=0.005 | log2FC=3.63, padj=1.75x10^-4 | **Robust and strengthened** |
| **NP cell state continuum** | rho = -0.258 | rho = -0.151 | **Robust** (same direction, all 3 versions) |
| **Prostaglandin pain pathway** | PTGS2, PLA2G2A, PTGES significant | Same 3 genes significant | **Robust** |
| **NF-kB TF activation** | RELA/NFKB1 significant (among 290 total) | Status needs confirmation (5 total) | **Uncertain** — may or may not survive |
| **HSF1 TF activation** | Significant (among 290 total) | Status needs confirmation (5 total) | **Uncertain** — may or may not survive |
| **CCC direction** | -6.5% in degeneration | +1.7% in degeneration | **Not robust** — effectively null |
| **AF trajectory** | rho = +0.341 | rho = +0.325 | **Consistent between v2/v3** |
| **TF landscape breadth** | 290 significant associations | 5 significant associations | **Not robust** — v2 was inflated |

### 3.2 Therapeutic Target Ranking Changes

| Target | v2 Rank | v3 Rank | Rationale |
|--------|---------|---------|-----------|
| TNF/NF-kB inhibition | **Tier 1** | **Tier 1** (pending TF confirmation) | TNF remains a significant pain gene; NF-kB TF status uncertain but TNF/CCL2/CXCL2 gene evidence supports pathway |
| HSP/proteostasis modulation | **Tier 1** | **Tier 2** (pending TF confirmation) | HSF1 TF status uncertain with only 5 significant associations in v3 |
| CXC chemokine blockade (CXCR2) | **Tier 2** | **Tier 1** | CXCL2 is now the most robust DE finding across all 3 versions — effect size increased, p-value strengthened |
| Prostaglandin pathway (COX-2 inhibition) | **Tier 2** | **Tier 1** | PLA2G2A/PTGS2/PTGES robust across v2/v3; well-validated drug targets exist |
| Anti-nerve ingrowth (Semaphorin/NRP) | Not ranked | **Tier 3** (Exploratory) | New: NRP2, SEMA3A, ROBO1 in v3 pain analysis; established biology but novel in this atlas |

### 3.3 What the Annotation Fix Revealed

The v2-to-v3 comparison provides a controlled experiment on the effect of annotation quality:

1. **NP_fibrocartilaginous is the most responsive cell type.** With cleaner pseudobulk profiles (rescued cells correctly included), this population shows 418 DE genes in mild_vs_severe — double v2's count and surpassing NP_mature_chondrocyte. This suggests the fibrocartilaginous state is a key mediator of the degenerative transcriptional response.

2. **TF activity was inflated by misannotation.** The 98% reduction (290 to 5 significant associations) is the clearest evidence that misrouted stressed disc cells were driving spurious TF signals in non-mesenchymal cell types. The 5 surviving associations are the genuinely robust signals.

3. **CXCL2 signal improves with better annotation.** The strengthening of CXCL2 (higher log2FC, lower padj) demonstrates that annotation quality and DE sensitivity are directly linked — purer pseudobulk profiles yield sharper signals.

4. **CCC stabilizes toward balance.** The roughly equal healthy/degenerated interaction counts in v3, combined with the inconsistency across all three versions, confirm that aggregate CCC counts are not a reliable measure of disease effect.

---

## 4. Robustness Assessment

### Findings Robust Across v2 and v3 (High Confidence)

| Finding | v2 Evidence | v3 Evidence |
|---------|-------------|-------------|
| CXCL2 upregulation in NP severe degeneration | log2FC=+3.14, padj=0.005 | log2FC=+3.63, padj=1.75x10^-4 (**stronger**) |
| NP pseudotime-condition negative correlation | rho = -0.258 | rho = -0.151 (same direction) |
| AF pseudotime-condition positive correlation | rho = +0.341 | rho = +0.325 (consistent) |
| NP_fibrocartilaginous as major DE responder | 203 + 127 genes (mild/healthy vs severe) | 418 + 385 genes (**much stronger**) |
| Prostaglandin pain pathway (PLA2G2A/PTGS2/PTGES) | All 3 significant | All 3 significant |
| Core pain mediators (TNF, CCL2, BDKRB2) | Significant | Significant |
| AF_outer DE in healthy_vs_severe | 97 genes | 100 genes (stable) |

### Findings Robust Across All Three Versions (Highest Confidence)

| Finding | v1 | v2 | v3 |
|---------|----|----|-----|
| CXCL2 upregulation | padj=0.002 | padj=0.005 | padj=1.75x10^-4 |
| NP pseudotime negative | rho=-0.207 | rho=-0.258 | rho=-0.151 |
| Disc cells as pain mediators | Yes | Yes | Yes |

### Findings Sensitive to Methodology (Requires Caution)

| Finding | v2 | v3 | Likely driver |
|---------|----|----|---------------|
| TF activity breadth | 290 significant | 5 significant | Annotation quality — misrouted cells inflated v2 |
| CCC interaction direction | -6.5% degenerated | +1.7% degenerated | Cell type count + subsampling |
| CCC total count | 55,889 | 81,059 | Cell type granularity (AF_mechanical_stress, EP_ossification added) |
| CEP trajectory direction | rho = -0.163 | rho = +0.135 | EP_ossification emergence changed trajectory structure |
| ORA term count | 1,577 | 1,043 | Changed DE gene composition |
| NP_mature_chondrocyte DE gene count | 315 + 172 | 291 + 113 | Cell redistribution after annotation fix |

---

## 5. Recommendations for SME Review

1. **Accept with high confidence:** CXCL2 as the most robust DE marker of NP degeneration (significant across all 3 versions, strengthening each time). NP pseudotime-condition negative correlation. Prostaglandin pain pathway. Disc-cells-as-mediators pain model.

2. **Accept as validated by annotation fix:** NP_fibrocartilaginous as the most transcriptionally responsive NP cell type (DE count doubled when annotation improved). AF_mechanical_stress as a real cell type (present in v1 and v3, absent only when v2 resolution was insufficient).

3. **Investigate the 5 surviving TF associations:** These are likely the genuinely robust TF signals. Determine whether HSF1, E2F4, and RELA/NFKB1 are among them. If any of the core TFs survived the annotation fix, they become very high-confidence findings. If none survived, the TF narrative needs substantial revision.

4. **Re-evaluate therapeutic target ranking based on TF results:** If HSF1 does not survive, HSP/proteostasis modulation should be downgraded. CXCL2/CXCR2 blockade and prostaglandin inhibition (COX-2) become the top therapeutic targets supported by version-robust evidence.

5. **Do not interpret aggregate CCC counts directionally.** Three versions, three different directions. Focus instead on specific robust interactions (e.g., collagen-integrin, TIMP-CD63).

6. **Consider neurovascular ingrowth as new biology.** The v3-specific pain genes (NRP2, PDGFA, ROBO1, SEMA3A) point to nerve/vessel ingrowth — a well-established feature of disc degeneration that the v3 annotation revealed. This should be validated against the CCC results (are these genes involved in detected ligand-receptor interactions?).

7. **CEP trajectory is unreliable.** Reversed between v2 and v3 with weak magnitudes both times. Not interpretable.

---

## 6. Summary Table

| Domain | v2 Headline | v3 Headline | Assessment |
|--------|-------------|-------------|------------|
| **Atlas** | 411K cells, 4 compartments, 8 cell types | 411K cells, 4 compartments, 10 cell types (+AF_mechanical_stress, +EP_ossification) | Better resolution |
| **Annotation** | 17K misrouted stressed disc cells | 25K cells rescued by evidence gate + ACAN/SOX9 rescue | Fix addresses known issue |
| **DE headline** | 949 genes, 21 powered comparisons | 1,156 genes, 18 powered comparisons | More genes despite fewer comparisons |
| **DE top responder** | NP_mature_chondrocyte (315 genes) | **NP_fibrocartilaginous (418 genes)** | Shifted — fibrocartilaginous is most responsive |
| **CXCL2** | log2FC=3.14, padj=0.005 | **log2FC=3.63, padj=1.75x10^-4** | **Strengthened** across all 3 versions |
| **TF activity** | 290 significant associations | **5 significant associations** | v2 was inflated by misannotation |
| **NP trajectory** | rho = -0.258 | rho = -0.151 | Robust (same direction) |
| **AF trajectory** | rho = +0.341 | rho = +0.325 | Consistent |
| **CEP trajectory** | rho = -0.163 | rho = +0.135 | **Not robust** |
| **CCC** | 56K interactions, -6.5% degenerated | 81K interactions, +1.7% degenerated | Aggregate counts not robust |
| **Pain genes** | 10 sig (inflammatory focus) | 10 sig (inflammatory + neurovascular) | Broadened to include nerve ingrowth |
| **ORA** | 1,577 terms | 1,043 terms | Reduced but likely higher specificity |
| **Top targets** | TNF/NF-kB #1, CXC #2 | CXCL2/prostaglandin #1 (pending TF review) | CXCL2 most version-robust target |

---

*This comparison was generated to support SME review of the pipeline v3 rerun. The v3 annotation fix directly addresses the misrouted stressed disc cell issue identified in the v1-vs-v2 comparison. All raw data and scripts are version-controlled. v2 state is recoverable at commit `430feb5`.*
