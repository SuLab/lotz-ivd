# Pipeline v3 vs v4: Process, Results, and Interpretation Comparison

**IVD Single-Cell Atlas — Internal Versioning Report**

Pipeline v3 final commit: `6622221` (2026-03-10)
Pipeline v4 final commit: `34f0312` (2026-03-11)

---

## Executive Summary

v4 is the largest structural change in the pipeline's history: scANVI replaces scVI, the pipeline expands from 10 to 12 modules, and a two-stage annotation system replaces single-pass de novo annotation. The 5 coarse anchor categories (Chondrocyte_like, Fibroblast_like, Immune, Endothelial, Pericyte_SMC) give scANVI biologically meaningful priors for semi-supervised integration. The result is 19 cell types (vs ~10 in v3), 23 powered DE comparisons (vs 18), and recovery of TF activity analysis (246 vs 5 associations). However, these gains come with tradeoffs: 33% fewer unique DE genes (772 vs 1,156), loss of the flagship CXCL2-in-NP_mature_chondrocyte finding, further weakening of the NP trajectory signal, and collapse of the AF trajectory to near-zero. The core biological conclusions — CXCL2 upregulation, prostaglandin pain pathway, disc cells as inflammatory mediators — survive the transition.

---

## 1. Process Changes

### 1.1 Dataset Scope

No change. Both versions use the same 11 datasets and 410,759 cells.

### 1.2 Pipeline Structure

| Feature | v3 | v4 | Impact |
|---------|----|----|--------|
| **Modules** | 10 (01-10) | **12** (01-12) | Clustering and annotation split from integration |
| **Module 04** | Binary mesenchymal/non-mesenchymal + 3 fixes | **5 coarse anchor categories** + Unknown | Richer semi-supervised priors for scANVI |
| **Module 05** | scVI integration + clustering + annotation | **scANVI integration only** | Semi-supervised; checkpoint resume |
| **Module 06** | (part of Module 05) | **Clustering with resolution optimization** | Adaptive resolutions by dataset size |
| **Module 07** | (part of Module 05) | **Two-stage annotation** (coarse → fine) | Stage 1: marker scoring; Stage 2: cluster DE |
| **Modules 08-12** | Modules 06-10 | **Renumbered** (DE, interpretation, trajectory, CCC, reporting) | Same analyses, new module numbers |

### 1.3 Integration Method

| Feature | v3 | v4 |
|---------|----|----|
| **Method** | scVI (unsupervised) | **scANVI (semi-supervised)** |
| **Anchor labels** | None | 5 coarse categories from Module 04 |
| **Workflow** | scVI (max_epochs=200) | scVI pre-training (200 epochs) → scANVI fine-tuning (50 epochs, early stopping) |
| **Unlabeled cells** | N/A | `Unknown` category positioned by scANVI similarity |
| **Batch key** | study | study |
| **Tiered** | Yes (mesenchymal / non-mesenchymal) | Yes (mesenchymal / non-mesenchymal) |

**Impact:** scANVI uses the coarse labels as soft constraints during integration, theoretically producing better batch correction while preserving biologically meaningful variation. The practical effect is mixed (see Section 2).

### 1.4 Annotation

| Feature | v3 | v4 |
|---------|----|----|
| **Strategy** | Single-pass de novo scoring | **Two-stage: coarse markers → fine DE** |
| **Stage 1** | Canonical marker scoring (specificity-weighted) | Canonical marker scoring (coarse identity) |
| **Stage 2** | None | **Within-group refinement using cluster DE markers** |
| **Cell types resolved** | ~10 | **19** |
| **CellTypist validation** | Yes (NP 8/13 discordant) | Not reported |

### 1.5 Clustering

| Feature | v3 | v4 |
|---------|----|----|
| **Method** | Leiden | Leiden |
| **Resolution selection** | Fixed set | **Adaptive by dataset size** (3-20 resolutions) |
| **Modularity metric** | Computed | **Skipped for >100K cells** (performance) |
| **NP clusters** | Not separately reported | 56 mesenchymal + 6 non-mesenchymal = **62** |
| **AF clusters** | Not separately reported | **14** mesenchymal |
| **CEP clusters** | Not separately reported | **9** mesenchymal |
| **all_cells clusters** | Not separately reported | 62 mesenchymal + 8 non-mesenchymal = **70** |

---

## 2. Key Results Comparison

### 2.1 Cell Atlas Composition

| Feature | v3 | v4 | Change |
|---------|----|----|--------|
| **Total cells** | 410,759 | 410,759 | No change |
| **NP cells** | 262,967 | 262,967 | No change |
| **AF cells** | 84,610 | 84,568 | -42 (-0.05%) |
| **CEP cells** | 50,714 | 50,769 | +55 (+0.1%) |
| **Total cell types** | ~10 | **19** | +9 new types |
| **NP types** | 4 (NP_notochordal, NP_mature_chondrocyte, NP_fibrocartilaginous, NP_stressed_degen) | **10** (+Fibrochondrocyte_chondroid, Fibrochondrocyte_fibroid, Fibrochondrocyte_stressed, NP_stressed, Macrophage_M2, Pericyte_SMC) | Much finer resolution |
| **AF types** | 3 (AF_inner, AF_outer, AF_mechanical_stress) | **2** (AF_inner, AF_outer) | Lost AF_mechanical_stress |
| **CEP types** | 2 (EP_hyaline, EP_ossification) | **3** (EP_hyaline, Fibroblast_like, Fibrochondrocyte_chondroid) | Different types emerged |
| **Unassigned NP** | 0 | **17,607 (6.7%)** | New unassigned category |

**Key observations:**
- NP cell types nearly tripled (4 → 10), the primary structural achievement of v4's two-stage annotation.
- The Fibrochondrocyte lineage (chondroid, fibroid, stressed) is a new concept in v4, splitting what v3 called NP_fibrocartilaginous and NP_stressed_degen into finer categories.
- **AF_mechanical_stress**, which re-emerged in v3 (present in v1, absent in v2), is lost again in v4. This cell type appears to be on the boundary of detectability.
- **17,607 unassigned NP cells (6.7%)** is a new limitation — the two-stage annotation created gaps for cells that don't clearly match any canonical signature.

### 2.2 Differential Expression

| Metric | v3 | v4 | Change |
|--------|----|----|--------|
| **Powered comparisons** | 18 | **23** | +5 (+28%) |
| **Unique significant genes** | 1,156 | **772** | -384 (-33%) |
| **Gene-comparison pairs** | 1,447 | **966** | -481 (-33%) |
| **Avg genes per comparison** | ~64 | **~34** | -47% |

**The fundamental tradeoff:** More cell types enabled more testable comparisons but reduced the number of cells per pseudobulk group, lowering per-comparison statistical power. v4 tests 28% more comparisons but detects 33% fewer genes.

#### Top DE Comparisons

| Comparison | v3 genes | v4 genes | Change |
|------------|----------|----------|--------|
| NP_fibrocartilaginous mild_vs_severe | **418** | 305 | -27% |
| NP_fibrocartilaginous healthy_vs_severe | **385** | 182 | -53% |
| NP_mature_chondrocyte mild_vs_severe | 291 | 242 | -17% |
| NP_mature_chondrocyte healthy_vs_severe | 113 | — | Not in top comparisons |
| AF_outer healthy_vs_severe | 100 | 58 | -42% |
| AF_inner healthy_vs_severe | N/A | **52** | **New** (v4 only) |
| AF_inner healthy_vs_all | N/A | **39** | **New** (v4 only) |
| Fibrochondrocyte_chondroid mild_vs_severe | N/A | **14** | **New** (v4 only) |
| Fibrochondrocyte_stressed mild_vs_severe | N/A | **14** | **New** (v4 only) |

**Interpretation:** Every comparison that existed in both versions lost genes in v4. The gains are in new comparisons that were not testable with v3's coarser cell types (AF_inner, Fibrochondrocyte_chondroid, Fibrochondrocyte_stressed). Whether the breadth of new comparisons compensates for the depth loss is a judgment call.

#### CXCL2: The Flagship Finding

| Version | Cell type | Comparison | log2FC | padj |
|---------|-----------|------------|--------|------|
| v3 | NP_mature_chondrocyte | mild_vs_severe | +3.63 | **1.75e-4** |
| v4 | NP_mature_chondrocyte | mild_vs_severe | +3.37 | **NS (missing)** |
| v4 | Fibrochondrocyte_chondroid | mild_vs_severe | +3.90 | **0.034** |
| v4 | NP_fibrocartilaginous | mild_vs_severe | +2.29 | **2.5e-4** |
| v4 | NP_fibrocartilaginous | healthy_vs_mild | -3.91 | **3.5e-4** |
| v4 | T_cell | mild_vs_severe | +1.91 | **0.033** |

**Critical finding:** CXCL2 in NP_mature_chondrocyte — the most robust individual finding across v1-v3, with strengthening significance each version — lost FDR significance in v4. The pvalue and padj columns are empty/missing in the v4 DE results, suggesting the gene was not statistically testable (likely due to reduced pseudobulk group size after cell redistribution).

The CXCL2 signal migrated to Fibrochondrocyte_chondroid and NP_fibrocartilaginous. These types were part of the v3 NP population but are now separated. The biological signal (CXCL2 upregulation in severe NP degeneration) persists, but the clean narrative of a single cell type with improving significance across versions is broken.

### 2.3 Pathway Enrichment

| Metric | v3 | v4 | Change |
|--------|----|----|--------|
| **ORA significant terms** | 1,043 | **1,772** | +70% |
| **GSEA significant terms** | — | **2,024** | — |

The 70% increase in ORA terms despite fewer DE genes is noteworthy. The expanded cell type repertoire provides more distinct pathway contexts. More cell types × more comparisons = more pathway tests, even with fewer genes per comparison.

### 2.4 Transcription Factor Activity

| Metric | v3 | v4 | Change |
|--------|----|----|--------|
| **Significant TF-condition associations** | 5 | **246** | +49x |
| **Cell-type-comparison contexts tested** | — | **13** | — |

**Recovery from v3's collapse:** v3 found only 5 significant TF associations, likely because the annotation fix changed pseudobulk compositions in ways that eliminated most TF signals. v4's 246 associations (across 13 cell-type-comparison contexts) suggests the scANVI + two-stage annotation produces sufficiently stable pseudobulk profiles for TF inference.

The 246 number is close to v2's 290 (which was deemed "inflated by misannotation"). Whether v4's 246 represents genuine signal or a different flavor of inflation remains to be determined. The expanded cell type count (19 vs ~8 in v2) means there are more contexts to test, which inflates the total even if per-context TF counts are modest.

### 2.5 Trajectory Analysis

| Compartment | v3 rho | v4 rho | Change |
|:-----------:|:------:|:------:|:------:|
| **NP** | -0.151 | **-0.092** | Weakened (same direction) |
| **AF** | +0.325 | **+0.019** | Collapsed to near-zero |
| **CEP** | +0.135 | **+0.396** | Strengthened dramatically |

#### NP Trajectory: Continuing to Weaken

The NP pseudotime-condition correlation has weakened in each version since v2 (-0.258 → -0.151 → -0.092). The direction remains negative (healthy early, degenerated late), confirming the cell state continuum, but the signal is approaching zero. The v4 weakening may reflect the distribution of cells across 10 NP types diluting the trajectory structure, or scANVI producing a latent space where the notochordal-to-degenerative axis is less dominant.

#### AF Trajectory: Collapsed

The AF trajectory, which appeared to stabilize in v2/v3 (+0.341/+0.325), has collapsed to +0.019 in v4. Across all four versions: -0.177, +0.341, +0.325, +0.019. This metric is unreliable for AF and should not be interpreted biologically.

#### CEP Trajectory: Strongest Signal Ever

The CEP trajectory strengthened from +0.135 (v3) to +0.396 (v4), the strongest trajectory signal in any compartment in any version. This is consistent with degenerated CEP cells occupying later pseudotime states. However, the CEP trajectory reversed sign between v2 (-0.163) and v3 (+0.135), so the v4 strengthening should be interpreted cautiously. If robust, it suggests that scANVI integration better resolves the CEP degeneration axis than scVI.

### 2.6 Cell-Cell Communication

| Metric | v3 | v4 | Change |
|--------|----|----|--------|
| **Healthy interactions** | 40,000 | **39,236** | -2% |
| **Degenerated interactions** | 41,000 | **37,013** | -10% |
| **Total interactions** | 81,000 | **76,249** | -6% |
| **Direction** | Balanced (ratio 1.02:1) | **Fewer in degeneration (-5.7%)** | Shifted |
| **Pain-relevant interactions** | — | **3,184** | — |

The CCC counts shifted slightly from v3's near-balance to v4 showing fewer interactions in degeneration (-5.7%). This aligns with v2's direction (-6.5%). The cross-version pattern (v1: +20%, v2: -6.5%, v3: balanced, v4: -5.7%) continues to demonstrate that aggregate interaction counts are not a stable metric.

### 2.7 Pain Biology

| Metric | v3 | v4 |
|--------|----|----|
| **Significant pain genes** | 10 | **7** |
| **Shared** | — | PTGS2, PLA2G2A, BDKRB2, CCL2, PTGES (5 shared) |
| **v3 only** | TNF, NRP2, PDGFA, ROBO1, SEMA3A | — |
| **v4 only** | — | FGF2, VEGFA |

**Core inflammatory pain axis preserved:** PTGS2/PLA2G2A/PTGES/CCL2/BDKRB2 are significant in both versions.

**Shifts:**
- v3's axon guidance genes (NRP2, ROBO1, SEMA3A, PDGFA) dropped below significance. These were novel to v3 and did not replicate in v4.
- v4 gains FGF2 and VEGFA (neovascularization), supporting the model that degenerated discs promote vascular ingrowth.
- TNF continues its decline — significant only in v1, it is now absent from v4's pain gene list entirely.

**New headline finding:** PTGS2 in AF_inner reaches padj=5.1e-8 (healthy_vs_severe), the most significant pain-relevant gene result in any version. This shifts the pain narrative from NP-centric CXCL2 (v1-v3) toward a broader NP+AF inflammatory signature.

---

## 3. Interpretation Changes

### 3.1 Mechanistic Model

| Component | v3 Narrative | v4 Narrative | Assessment |
|-----------|-------------|-------------|------------|
| **CXCL2 in NP** | log2FC=3.63, padj=1.75e-4 (NP_mature_chondrocyte) | log2FC=3.90, padj=0.034 (Fibrochondrocyte_chondroid); log2FC=2.29, padj=2.5e-4 (NP_fibrocartilaginous) | **Signal migrated** — present but in different cell types |
| **PTGS2 as top pain gene** | Significant | **padj=5.1e-8 in AF_inner** — strongest result ever | **Strengthened** |
| **NP cell state continuum** | rho = -0.151 | rho = -0.092 | **Weakened** but same direction |
| **Prostaglandin pathway** | PTGS2/PLA2G2A/PTGES significant | Same 3 genes significant | **Robust** |
| **TF activity** | 5 significant associations | 246 significant associations | **Recovered** — v3 was likely artifactually low |
| **CCC direction** | Balanced (40K/41K) | Fewer in degeneration (39K/37K) | **Not robust** |
| **AF trajectory** | rho = +0.325 | rho = +0.019 | **Collapsed** |
| **CEP trajectory** | rho = +0.135 | rho = +0.396 | **Strengthened** |
| **Pain gene narrative** | Inflammatory + axon guidance (NRP2/ROBO1/SEMA3A) | Inflammatory + neovascularization (FGF2/VEGFA) | **Shifted** — different supporting genes |

### 3.2 Therapeutic Target Ranking Changes

| Target | v3 Rank | v4 Rank | Rationale |
|--------|---------|---------|-----------|
| Prostaglandin pathway (COX-2 inhibition) | Tier 1 | **Tier 1** | PTGS2 now the strongest pain gene result ever (padj=5.1e-8). Strengthened. |
| CXC chemokine blockade (CXCR2) | Tier 1 | **Tier 1** (with caveat) | CXCL2 still significant but in different cell types; signal less clean |
| TNF/NF-kB inhibition | Tier 1 (pending TF) | **Tier 2** | TNF not significant in v4 pain genes. NF-kB supported at pathway level only. |
| Anti-nerve ingrowth (Semaphorin/NRP) | Tier 3 (exploratory) | **Dropped** | NRP2/ROBO1/SEMA3A from v3 did not replicate |
| Anti-neovascularization (FGF2/VEGFA) | Not ranked | **Tier 2** (new) | FGF2 and VEGFA significant in v4 pain analysis |

### 3.3 What scANVI + Two-Stage Annotation Revealed

1. **Fibrochondrocyte subtypes are real populations.** v4 resolves Fibrochondrocyte_chondroid (18K cells), Fibrochondrocyte_fibroid (3.6K), and Fibrochondrocyte_stressed (4.2K) within what v3 called NP_fibrocartilaginous and NP_stressed_degen. These have distinct DE signatures (14 genes each in mild_vs_severe for chondroid and stressed types). Whether this level of granularity is useful depends on validation.

2. **AF_inner is a distinct functional compartment.** v4's AF_inner (35K cells) shows 52 DE genes in healthy_vs_severe and 39 in healthy_vs_all — a previously invisible signal. The separation of AF into inner and outer with distinct DE profiles is a genuine gain from the two-stage annotation.

3. **CEP Fibroblast_like is a new population.** 17K CEP cells classified as Fibroblast_like — a population not seen in v3 (which had EP_hyaline and EP_ossification). The replacement of EP_ossification by Fibroblast_like suggests different cell type definitions rather than new biology.

4. **The NP taxonomy is overfit.** 10 NP types from 262K cells is ambitious. The 17,607 unassigned cells (6.7%) and the loss of power per comparison suggest the annotation may be too fine for the available data.

---

## 4. Robustness Assessment

### Findings Robust Across v3 and v4 (High Confidence)

| Finding | v3 Evidence | v4 Evidence |
|---------|-------------|-------------|
| CXCL2 upregulation in NP severe degeneration | padj=1.75e-4 (NP_mature_chondrocyte) | padj=0.034 (Fibrochondrocyte_chondroid), padj=2.5e-4 (NP_fibrocartilaginous) |
| Prostaglandin pain pathway | PTGS2/PLA2G2A/PTGES significant | Same 3 + PTGS2 padj=5.1e-8 in AF_inner (**stronger**) |
| NP pseudotime negative direction | rho = -0.151 | rho = -0.092 (same direction, weaker) |
| NP_fibrocartilaginous as top DE responder | 418 genes (m_vs_s) | 305 genes (m_vs_s) — still #1 |
| Core pain mediators (CCL2, BDKRB2) | Significant | Significant |
| Disc cells as inflammatory mediators | Supported | Supported |

### Findings Robust Across All Four Versions (Highest Confidence)

| Finding | v1 | v2 | v3 | v4 |
|---------|----|----|-----|-----|
| CXCL2 upregulation in NP | padj=0.002 | padj=0.005 | padj=1.75e-4 | padj=0.034 / 2.5e-4 |
| NP pseudotime negative | rho=-0.207 | rho=-0.258 | rho=-0.151 | rho=-0.092 |
| Disc cells as pain mediators | Yes | Yes | Yes | Yes |
| Prostaglandin pathway | — | Sig | Sig | Sig (strongest ever) |

### Findings Sensitive to Methodology (Requires Caution)

| Finding | v3 | v4 | Likely driver |
|---------|----|----|---------------|
| CXCL2 cell type assignment | NP_mature_chondrocyte | Fibrochondrocyte_chondroid / NP_fibrocartilaginous | Annotation granularity |
| TF activity count | 5 | 246 | Annotation + cell type count |
| AF trajectory | rho = +0.325 | rho = +0.019 | Integration method (scVI → scANVI) |
| NP trajectory magnitude | rho = -0.151 | rho = -0.092 | Cell type redistribution |
| DE gene count | 1,156 | 772 | Cell type granularity tradeoff |
| Axon guidance pain genes (NRP2/ROBO1/SEMA3A) | Significant | NS | Single-version finding (v3 only) |
| CCC direction | Balanced | -5.7% | Still unstable |

---

## 5. Recommendations for SME Review

1. **Decide on publication version.** v3 has a cleaner narrative (CXCL2 in a single cell type, more DE genes, simpler taxonomy). v4 has richer biology (19 types, AF_inner DE, PTGS2 as strongest pain gene, TF recovery). Consider v3 as primary with v4 as sensitivity analysis, or vice versa.

2. **Evaluate the 19-type taxonomy.** Are Fibrochondrocyte_chondroid, Fibrochondrocyte_fibroid, and Fibrochondrocyte_stressed biologically meaningful or overfit? The 17,607 unassigned cells suggest the annotation may be too aggressive.

3. **Lead with PTGS2 if using v4.** PTGS2 in AF_inner (padj=5.1e-8) is the most significant pain gene result across all versions. This is a stronger finding than any CXCL2 result but shifts the narrative to AF rather than NP.

4. **Do not present v3 axon guidance genes as robust.** NRP2, ROBO1, SEMA3A, and PDGFA were significant only in v3 and did not replicate in v4. These should not be published as confirmed findings.

5. **AF trajectory should be removed from the biological narrative.** Four versions, four different results (including near-zero in v4). There is no stable signal.

6. **CEP trajectory deserves cautious attention.** The +0.396 rho is the strongest trajectory result in any version, but the previous sign reversal (v2→v3) means it needs independent validation.

7. **Frame the DE power tradeoff explicitly.** If presenting v4, acknowledge that finer cell types come at the cost of per-comparison power. The 33% reduction in DE genes is the direct price of the 19-type taxonomy.

---

## 6. Summary Table

| Domain | v3 Headline | v4 Headline | Assessment |
|--------|-------------|-------------|------------|
| **Atlas** | 411K cells, ~10 types | 411K cells, **19 types** | More resolution, some unassigned |
| **Pipeline** | 10 modules, scVI | **12 modules, scANVI** | Major restructuring |
| **DE headline** | 1,156 genes, 18 comparisons | **772 genes, 23 comparisons** | More breadth, less depth |
| **Top DE** | NP_fibrocartilaginous m_vs_s (418) | NP_fibrocartilaginous m_vs_s (**305**) | Same top type, fewer genes |
| **CXCL2** | padj=1.75e-4 (NP_mature_chondrocyte) | **padj=NS** (NP_mature_chondrocyte); padj=0.034/2.5e-4 (other types) | Signal migrated |
| **PTGS2** | Significant | **padj=5.1e-8** (AF_inner) | **Strongest pain gene ever** |
| **TF activity** | 5 significant | **246 significant** | Recovered |
| **NP trajectory** | rho = -0.151 | rho = **-0.092** | Weakened |
| **AF trajectory** | rho = +0.325 | rho = **+0.019** | Collapsed |
| **CEP trajectory** | rho = +0.135 | rho = **+0.396** | Strengthened |
| **CCC** | 81K interactions, balanced | 76K interactions, **-5.7% degen** | Still unstable |
| **Pain genes** | 10 (inflammatory + axon guidance) | **7** (inflammatory + neovascularization) | Core preserved, periphery shifted |
| **ORA pathways** | 1,043 | **1,772** | More terms |
| **Unassigned** | 0 | **17,607 NP cells (6.7%)** | New limitation |

---

*This comparison was generated to support SME review of the v4 pipeline rerun. The v4 changes (scANVI, two-stage annotation, 12-module pipeline) represent the most structurally significant version transition. All raw data and scripts are version-controlled. v3 state is recoverable at commit `6622221`. v4 state at commit `34f0312`.*
