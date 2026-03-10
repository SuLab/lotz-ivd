# Pipeline v1 vs v2: Process, Results, and Interpretation Comparison

**IVD Single-Cell Atlas — Internal Versioning Report**

Pipeline v1 final commit: `c950d1d` (2026-03-05)
Pipeline v2 final commit: `430feb5` (2026-03-10)

---

## Executive Summary

The pipeline was rerun after a spec restructuring that changed three core methodological choices: (1) cell type annotation was split into binary classification (Module 04) followed by post-integration de novo annotation (Module 05), (2) integration moved from a 4-approach benchmark with scANVI primary to scVI-only with four compartment objects, and (3) GSE233666 was excluded. These changes substantially altered downstream results. Some findings were robust across both versions (HSF1/NF-kB TF activation, HSP/mitochondrial GSEA, NP pseudotime-condition correlation direction), while others were sensitive to methodology (CXC chemokine significance, AF trajectory direction, CCC interaction counts). This report documents all differences to support informed SME review.

---

## 1. Process Changes

### 1.1 Dataset Scope

| Feature | v1 | v2 | Impact |
|---------|----|----|--------|
| Datasets | 12 | 11 | GSE233666 (7 herniated NP samples, 22,658 cells) excluded |
| Total cells | 436,239 | 410,759 | -25,480 cells |
| Samples | 78 | 71 | -7 samples |
| Donors | 57 | ~50 | -7 donors |

**Rationale for GSE233666 exclusion:** All samples are herniated, a single-condition dataset. Including it inflated herniated sample counts but provided no within-study condition contrasts. After its removal, herniated samples come only from GSE251686, making all herniated comparisons fully confounded with study.

### 1.2 Cell Type Annotation

| Feature | v1 | v2 |
|---------|----|----|
| **Module 04** | Fine-grained annotation: 16 IVD-specific gene signatures + CellTypist Immune_All_Low. Produces detailed cell types (NP_notochordal, NP_mature_chondrocyte, AF_inner, AF_outer, etc.) per dataset | Binary classification: mesenchymal vs non-mesenchymal. Simple marker-based scoring |
| **Module 05 annotation** | None — cell types carry forward from Module 04 | De novo annotation: Leiden clustering at optimized resolution → cluster-level marker scoring against canonical panel → CellTypist validation for non-mesenchymal clusters |
| **Label column** | `cell_type_final` | `cell_type` |
| **Annotation timing** | Pre-integration (per dataset) | Post-integration (on integrated objects) |
| **Validation** | Marker scoring + CellTypist agreement per dataset | CellTypist concordance for non-mesenchymal clusters (NP: 5/13 concordant, AF: good, CEP: partial) |

**Impact:** Post-integration annotation in v2 means cell types are defined on the integrated manifold rather than per-dataset. This changes pseudobulk groupings and therefore DE results. The v2 approach also revealed NP_fibrocartilaginous as a distinct population and enabled CEP-specific annotation (EP_hyaline).

### 1.3 Integration Strategy

| Feature | v1 | v2 |
|---------|----|----|
| **Structure** | Two tiers: Tier 1 (non-resident, 14.6K cells), Tier 2 (resident NP 139K + AF 283K) | Four compartment objects: NP (263K), AF (85K), CEP (51K), all_cells (411K) |
| **Methods benchmarked** | 4 approaches: scVI, scANVI, Harmony, BBKNN | scVI only |
| **Primary integration** | scANVI (best composite score 0.615, cell type ASW 0.521) | scVI |
| **Sensitivity check** | scVI (for trajectory — preserves continuum, variance ratio = 1.0) | None (single method) |
| **Cell partitioning** | Non-resident vs resident (based on Module 04 detailed annotation) | Mesenchymal vs non-mesenchymal within each compartment object (based on Module 04 binary classification) |

**Impact:** The shift from scANVI (semi-supervised, leveraging pre-integration labels) to scVI (unsupervised) removes the dependency on Module 04's detailed annotation quality. The 4-compartment structure includes CEP as a separate object (v1 had no dedicated CEP integration) and splits the v1 AF group (282K cells) into AF (85K) and portions redistributed to NP and all_cells.

### 1.4 Herniated Comparisons

| Feature | v1 | v2 |
|---------|----|----|
| **Status** | Included as "exploratory" | Excluded entirely |
| **Rationale** | Flagged as study-confounded (RPL genes in top hits) but results reported | After GSE233666 removal, only GSE251686 contributes herniated → fully confounded with study |
| **Impact on DE** | 4,316 genes (NP_mature_chondrocyte healthy_vs_herniated) dominated the total DE count | Removed from analysis; total DE gene count drops substantially |

### 1.5 Hardware

| Feature | v1 | v2 |
|---------|----|----|
| RAM | 30 GB | 62 GB |
| CPUs | 4 | 16 |
| GPU | A10G (23GB) | A10G (23GB) |
| Disk | 123 GB | 123 GB |

---

## 2. Key Results Comparison

### 2.1 Cell Atlas Composition

| Feature | v1 | v2 | Change |
|---------|----|----|--------|
| **NP cells** | 138,937 (resident tier) | 262,967 (compartment object) | +89% — v2 includes non-mesenchymal cells in the NP object |
| **AF cells** | 282,736 (resident tier) | 84,624 | -70% — v2 splits by compartment annotation, not by cell class |
| **CEP cells** | Not separately integrated | 50,858 | New |
| **Non-resident** | 14,566 (Tier 1) | Distributed within compartment objects | Structural change |
| **NP subtypes** | NP_notochordal, NP_mature_chondrocyte, NP_stressed_degenerative | NP_notochordal, NP_mature_chondrocyte, NP_fibrocartilaginous, NP_stressed_degenerative | +NP_fibrocartilaginous (new) |
| **AF subtypes** | AF_inner, AF_outer, AF_mechanical_stress | AF_inner, AF_outer | -AF_mechanical_stress |
| **CEP subtypes** | EP_hyaline_cartilage, EP_ossification (in annotation, not integrated) | EP_hyaline | Simplified |
| **Immune types** | CellTypist-refined: Tcm/Naive helper T, Tem/Trm cytotoxic T, Macrophages, B cells | De novo: T_cell, B_cell, Macrophage, NK_cell, Pericyte_SMC | Less granular immune types |

**Key observation:** The dramatic change in AF cell counts (283K → 85K) suggests that v1's compartment assignment (based on dataset-of-origin metadata) grouped many cells into AF that v2's integration-based approach assigns differently. This fundamentally changes the AF populations being compared in DE and trajectory analyses.

### 2.2 Differential Expression

| Metric | v1 | v2 | Change |
|--------|----|----|--------|
| Powered comparisons | 17 | 21 | +4 (including EP_hyaline, T_cell, Macrophage, B_cell) |
| Skipped comparisons | 128 | 53 | -75 (more cell types × comparisons powered) |
| Total gene×comparison pairs | ~5,328* | 1,231 | -77% |
| Unique significant genes | ~1,012** | 949 | -6% |
| CEP powered | No | Yes (EP_hyaline: 84 genes) | New |
| Herniated comparisons | Included (flagged) | Excluded | Removed |

*v1's "5,328" included herniated comparisons; NP_mature_chondrocyte healthy_vs_herniated alone contributed 4,316.
**Estimated unique gene count excluding herniated in v1.

#### Top DE Genes: NP_mature_chondrocyte mild_vs_severe

| Gene | v1 log2FC | v1 padj | v2 log2FC | v2 padj | Status |
|------|-----------|---------|-----------|---------|--------|
| **CXCL1** | +3.75 | 0.0014 | +1.59 | NS | **Lost significance** |
| **CXCL3** | +3.72 | 6.8×10⁻⁸ | +2.47 | 0.099 | **Lost significance** |
| **CXCL2** | +3.13 | 0.002 | +3.14 | 0.005 | **Robust** |
| **MDK** | +2.72 | 4.9×10⁻¹² | +0.87 | NS | **Lost significance** |
| **TNF** | +2.45 | 0.043 | +1.23 | 0.22 | **Lost significance** |
| **CEMIP** | — | — | +1.50 | 0.055 | Borderline in v2 |
| **UGCG** | — | — | +1.06 | 6.3×10⁻⁶ | **New top gene in v2** |
| **ICAM1** | — | — | +1.67 | 1.9×10⁻⁵ | **New top gene in v2** |
| **HHEX** | — | — | +2.78 | 1.9×10⁻⁵ | **New top gene in v2** |
| **TNC** | — | — | +2.84 | 0.017 | **New in v2** |
| **IL32** | — | — | +2.57 | 0.014 | **New in v2** |

**Interpretation:** The CXC chemokine triad (CXCL1/2/3) that dominated v1's narrative is largely non-significant in v2. Only CXCL2 remains robust across both versions. The loss of CXCL1, CXCL3, TNF, and MDK significance is likely driven by changes in cell groupings from the different annotation approach (post-integration de novo vs pre-integration marker-based). In v2, cells previously grouped as "NP_mature_chondrocyte" may now be split across NP_mature_chondrocyte and NP_fibrocartilaginous, changing the pseudobulk aggregations and statistical power for individual genes.

New v2 top genes (UGCG, ICAM1, HHEX, TNC, IL32) are biologically coherent: ICAM1 and IL32 are inflammatory mediators, TNC is an inflammation-induced ECM protein, and UGCG participates in glycosphingolipid-mediated inflammatory signaling.

#### New v2 Findings Not Present in v1

| Finding | Details |
|---------|---------|
| **NP_fibrocartilaginous DE** | 203 genes (mild_vs_severe), 127 genes (healthy_vs_severe) — second most responsive NP cell type |
| **EP_hyaline DE** | 84 downregulated genes (0 upregulated) in healthy_vs_degenerated_all — transcriptional silencing pattern |
| **T_cell DE** | 48 genes (mild_vs_severe) — immune cell DE now powered |
| **Macrophage DE** | 5 downregulated genes (mild_vs_severe) |

### 2.3 Pathway Enrichment

| Metric | v1 | v2 | Change |
|--------|----|----|--------|
| ORA significant | 1,244 | 1,577 | +27% |
| GSEA significant | 1,081 | 1,576 | +46% |

**Robust pathways (significant in both):**
- Chemokine-mediated signaling (NP upregulated)
- Neutrophil chemotaxis (NP upregulated)
- Cellular response to heat / unfolded protein response (AF upregulated)
- Oxidative phosphorylation (AF downregulated)
- Cell cycle pathways (NP degeneration)

**Change in emphasis:** v1's pathway narrative was dominated by the CXC chemokine signal. v2 shows broader inflammatory enrichment driven by a more diverse set of genes (ICAM1, IL32, TNC in addition to CXCL2), making the pathway signal more robust even as individual gene significance shifted.

### 2.4 Transcription Factor Activity

| Metric | v1 | v2 | Change |
|--------|----|----|--------|
| Significant TF-condition associations | 113 | 290 | +157% |

**Robust TFs (significant in both versions):**

| TF | Role | v1 Evidence | v2 Evidence |
|----|------|-------------|-------------|
| **HSF1** | Heat shock | padj 4.8×10⁻⁸ – 5.0×10⁻⁶ (3 cell types) | Remains significant |
| **HSF2** | Heat shock | padj 1.7×10⁻⁴ (endothelial) | Remains significant |
| **E2F4** | Cell cycle | padj 8.4×10⁻⁹ (NP_mature_chondrocyte) | Remains significant |
| **RELA** | NF-kB | padj 0.002 (AF_inner) | Remains significant |
| **NFKB1** | NF-kB | padj 8.7×10⁻⁴ (AF_inner) | Remains significant |
| **FOXO3** | Apoptosis | padj 7.4×10⁻⁴ (NP_stressed) | Remains significant |

The core TF findings are among the most version-robust results. HSF1, E2F4, and RELA/NFKB1 were significant in both v1 and v2, providing strong evidence that these transcriptional programs are genuine features of IVD degeneration rather than methodological artifacts.

### 2.5 Trajectory Analysis

| Metric | v1 | v2 | Change |
|--------|----|----|--------|
| Compartments analyzed | NP, AF | NP, AF, CEP | +CEP |
| NP pseudotime-condition rho | -0.207 | **-0.258** | Stronger (same direction) |
| AF pseudotime-condition rho | -0.177 | **+0.341** | **REVERSED** |
| CEP pseudotime-condition rho | — | -0.163 | New |
| NP trajectory-DE overlap | ~55% (278/500) | 19% (96/500) | Substantially lower |
| AF trajectory-DE overlap | ~55% (254/500) | 22% (110/500) | Substantially lower |
| Integration used | scANVI (primary), scVI (sensitivity) | scVI only | Single method |
| scVI sensitivity rho (NP) | -0.132 | N/A (scVI is primary) | — |

#### NP Trajectory: Robust

The NP pseudotime-condition correlation is **stronger in v2** (rho = -0.258 vs -0.207), confirming the notochordal → mature → stressed/degenerative continuum as a robust finding. The direction is consistent across both versions and across the v1 sensitivity check (scVI rho = -0.132).

#### AF Trajectory: REVERSED — Requires Investigation

**This is the most concerning v1-to-v2 discrepancy.** The AF pseudotime-condition correlation flipped from rho = -0.177 (healthy at earlier pseudotime) to rho = +0.341 (healthy at *later* pseudotime). Possible explanations:

1. **Cell composition change:** v1 AF had 283K cells (resident tier); v2 AF has 85K. The 70% reduction means a fundamentally different cell population is being analyzed. Cells that were in the v1 AF object may now be in the NP or all_cells objects.
2. **Root cell selection:** If the AF_inner root cluster has different cellular composition in v2, the pseudotime ordering could reverse.
3. **Annotation effect:** De novo annotation may assign different cell types to AF clusters, changing which cells are considered "inner" (root) vs "outer" (terminal).
4. **Integration method:** scANVI (v1) vs scVI (v2) produces different latent spaces, which directly affect the neighbor graph and therefore PAGA/DPT.

**Recommendation:** This finding should not be interpreted biologically until the root cause is identified. The v1 AF trajectory was already the weaker signal (rho = -0.177 vs NP's -0.207), and the reversal suggests it was on the boundary of methodological sensitivity.

#### Trajectory-DE Overlap: Substantially Reduced

The trajectory-DE overlap dropped from ~55% to 19-22%. This is partly explained by the different DE gene sets (949 vs ~1,012 unique genes, different cell type groupings) and partly by the different integration embeddings driving the trajectory. A lower overlap does not invalidate either analysis but indicates that the trajectory and DE analyses are capturing somewhat different aspects of disease biology in v2.

### 2.6 Cell-Cell Communication

| Metric | v1 | v2 | Change |
|--------|----|----|--------|
| Healthy interactions | 44,079 | 28,878 | -34% |
| Degenerated interactions | 53,036 | 27,011 | -49% |
| Direction | More in degeneration (+20%) | **Fewer in degeneration (-6.5%)** | **REVERSED** |
| Pain-relevant interactions | 3,662–4,194 | 2,077 | -50% |
| Cell types in healthy | 17 | — | — |
| Cell types in degenerated | 22 | — | — |

**The CCC direction reversal** (from "more interactions in degeneration" to "fewer") is a significant interpretive change. v1's narrative emphasized "increased signaling complexity" as a hallmark of degeneration. v2 contradicts this.

**Possible explanations:**
1. **Different cell type resolution:** v2's de novo annotation produces different (often fewer or broader) cell types. Since CCC interaction counts scale with the number of cell type pairs, fewer granular immune subtypes in v2 (e.g., T_cell vs Tcm/Naive + Tem/Trm) reduce combinatorial interaction counts.
2. **Different cell composition per condition:** The cell populations contributing to each condition differ between v1 and v2 due to annotation changes.
3. **Subsampling differences:** Both versions sample 20K cells per condition, but the cell type composition of those 20K cells differs.

**Recommendation:** CCC interaction counts are highly sensitive to cell type granularity and subsampling. The direction of change should not be over-interpreted. Instead, focus on specific biologically validated interactions (e.g., TIMP1-CD63, collagen-integrin) rather than aggregate counts.

### 2.7 Pain Biology

| Metric | v1 | v2 | Change |
|--------|----|----|--------|
| Significant pain genes | 3 (TNF ×2, CXCL8 ×1) | 10 (PTGS2, TNF ×2, PLA2G2A, BDKRB2, CCL2, PTGES, CXCL8) | +7 genes |
| Prostaglandin pathway | Not highlighted | PTGS2 + PLA2G2A + PTGES = complete pathway | New finding |
| Pain model | TNF/CXC chemokine → nerve sensitization | Broader: prostaglandin + chemokine + cytokine → nerve sensitization | Expanded |

The pain biology findings are **richer in v2**, with 10 significant genes spanning a more complete prostaglandin synthesis pathway. The core model (disc cells produce inflammatory mediators that sensitize nerves, not nociceptors) is consistent across both versions.

### 2.8 Supplementary Tables

| Version | Count | Notable additions in v2 |
|---------|-------|------------------------|
| v1 | 13 | — |
| v2 | 19 | S17-S19: CellTypist concordance (NP, AF, CEP) |

---

## 3. Interpretation Changes

### 3.1 Mechanistic Model

| Component | v1 Narrative | v2 Narrative | Robustness |
|-----------|-------------|-------------|------------|
| **Central driver** | CXC chemokine triad (CXCL1/2/3) → neutrophil recruitment | NF-kB activation → inflammatory mediator production (CXCL2, ICAM1, IL32) | **Shifted** — NF-kB TF evidence is robust; gene-level chemokine evidence is weaker |
| **TNF role** | Significant DE gene (padj=0.043); direct evidence | Not significant in DE (padj=0.22); supported only at TF level (RELA/NFKB1) and pain gene analysis | **Weakened** at gene level |
| **HSF1/HSP axis** | Strong signal (padj 5×10⁻⁶ – 5×10⁻⁸) | Remains strong | **Robust** |
| **Mitochondrial dysfunction** | GSEA: OXPHOS/ETC downregulated in AF | Remains significant | **Robust** |
| **Cell state continuum** | NP rho = -0.207 | NP rho = -0.258 (stronger) | **Robust** |
| **CCC complexity increase** | 20% more interactions in degeneration | 6.5% fewer interactions | **Not robust** |
| **Pain mechanism** | TNF + CXC chemokines | Prostaglandin pathway + CXCL2 + broader cytokines | **Expanded** |

### 3.2 Therapeutic Target Ranking Changes

| Target | v1 Rank | v2 Rank | Rationale for change |
|--------|---------|---------|---------------------|
| CXC chemokine blockade (CXCR2) | **Tier 1** (Strong) | **Tier 2** (Moderate) | Only CXCL2 significant; CXCL1/3 lost |
| TNF/NF-kB inhibition | **Tier 1** (Strong) | **Tier 1** (Strong) | TNF gene not sig, but NF-kB TF-level evidence remains |
| HSP/proteostasis modulation | **Tier 1** (Strong) | **Tier 1** (Strong) | Robust across versions |
| Prostaglandin pathway | Not ranked | **Tier 2** (Moderate) | New: PTGS2/PLA2G2A/PTGES pathway |
| CEMIP/hyaluronan | **Tier 2** | Not ranked (borderline) | CEMIP padj=0.055 in v2 |

### 3.3 Novel v2 Findings

1. **NP_fibrocartilaginous as a distinct population:** This cell type was not separately identified in v1. It shows 203 DE genes in mild_vs_severe, suggesting it captures a biologically meaningful transitional state.

2. **CEP transcriptional silencing:** EP_hyaline shows 84 downregulated and 0 upregulated genes — a pure silencing pattern not observed in NP or AF. This is the first evidence of CEP-specific degeneration transcriptomics in this atlas.

3. **Prostaglandin pain pathway:** The coherent identification of PLA2G2A → PTGS2 → PTGES as significant pain genes provides a mechanistically complete pathway from arachidonic acid release to PGE2-mediated nerve sensitization.

4. **CellTypist annotation validation:** The systematic comparison of de novo annotation with CellTypist revealed that NP non-mesenchymal annotation is unreliable (8/13 discordant), identifying ~17K misrouted stressed disc cells. This is a methodological insight not available from v1.

---

## 4. Robustness Assessment

### Findings Robust Across Both Versions (High Confidence)

| Finding | v1 Evidence | v2 Evidence |
|---------|-------------|-------------|
| CXCL2 upregulation in NP severe degeneration | log2FC=+3.13, padj=0.002 | log2FC=+3.14, padj=0.005 |
| HSF1/HSF2 TF activation across cell types | padj 5×10⁻⁸ – 5×10⁻⁶ | Remains significant |
| E2F4 TF activation in NP severe | padj=8.4×10⁻⁹ | Remains significant |
| RELA/NFKB1 TF activation in AF | padj 0.002 / 8.7×10⁻⁴ | Remains significant |
| NP pseudotime-condition negative correlation | rho = -0.207 | rho = -0.258 (stronger) |
| HSP/heat response GSEA enrichment | NES = +2.35 (AF_inner) | Significant |
| OXPHOS/mitochondrial GSEA suppression | NES = -1.96 (AF_inner) | Significant |
| Disc cells produce pain mediators, not nociceptors | TNF, CXCL1-3 sig | PTGS2, CXCL2, IL32 sig |
| No significant composition changes (FDR) | 0/58 | 0/N |
| Wnt/Notch/senescence NOT significant in GSEA | Not sig | Not sig |

### Findings Sensitive to Methodology (Requires Caution)

| Finding | v1 | v2 | Likely driver of difference |
|---------|----|----|---------------------------|
| CXC chemokine triad dominance | CXCL1/2/3 all significant | Only CXCL2 | Cell type grouping changes (annotation method) |
| TNF gene-level significance | padj=0.043 | padj=0.22 | Cell type grouping changes |
| MDK as degeneration marker | padj=4.9×10⁻¹² | Not significant | Cell type grouping changes |
| AF trajectory direction | rho = -0.177 | rho = +0.341 | Cell composition + integration method |
| CCC interaction count direction | +20% in degeneration | -6.5% in degeneration | Cell type granularity + subsampling |
| Trajectory-DE overlap | ~55% | 19-22% | Different DE gene sets + integration embeddings |
| Total DE gene count | ~5,328 (with herniated) | 949 (without) | Herniated exclusion + annotation changes |

---

## 5. Recommendations for SME Review

1. **Accept with confidence:** NP pseudotime-condition correlation, HSF1/NF-kB TF activation, HSP/mitochondrial GSEA, CXCL2 significance, prostaglandin pain pathway, disc-cells-as-mediators pain model.

2. **Investigate before interpreting:** AF trajectory reversal (likely cell composition artifact), CCC direction reversal (likely cell type granularity artifact).

3. **Treat as version-dependent:** CXC chemokine triad (CXCL1/3), TNF gene-level significance, MDK, aggregate CCC counts. These findings depend on how cells are grouped into types and should be interpreted cautiously until annotation is stabilized.

4. **Consider for rerun:** The CellTypist disagreement analysis identified ~17K misrouted stressed disc cells in the NP non-mesenchymal tier. Fixing Module 04 classification (requiring co-expression of ≥2 immune markers) would improve annotation quality and could restore some lost DE signals.

5. **Evaluate NP_fibrocartilaginous:** This is a genuinely new cell population in v2. Its substantial DE signal (203 genes) warrants biological validation — is it a real transitional population or an artifact of clustering resolution?

6. **CEP findings are new and preliminary:** EP_hyaline DE (84 downregulated genes) is enabled only by v2's compartment-specific integration. The biological interpretation (transcriptional silencing in degeneration) is plausible but based on a single comparison with limited sample counts.

---

## 6. Summary Table

| Domain | v1 Headline | v2 Headline | Assessment |
|--------|-------------|-------------|------------|
| **Atlas** | 436K cells, 12 datasets, 2-tier | 411K cells, 11 datasets, 4-compartment | Structural change; v2 more granular |
| **Integration** | scANVI primary (4 benchmarked) | scVI only | Simpler; removes annotation dependency |
| **Annotation** | Pre-integration, marker-based | Post-integration, de novo + CellTypist | More principled; reveals annotation issues |
| **DE headline** | CXC chemokines dominate NP degeneration | Broader inflammatory signature; CXCL2 + ICAM1 + IL32 | Attenuated but biologically coherent |
| **TF activity** | 113 significant; HSF1/E2F4/NF-kB | 290 significant; same core TFs | Robust core; v2 detects more |
| **NP trajectory** | rho = -0.207 | rho = -0.258 | Robust (stronger in v2) |
| **AF trajectory** | rho = -0.177 | rho = +0.341 | **Not robust** |
| **CCC** | More interactions in degeneration | Fewer interactions in degeneration | **Not robust** |
| **Pain genes** | 3 significant | 10 significant | v2 richer; prostaglandin pathway new |
| **Therapeutic targets** | CXC blockade #1 | TNF/NF-kB #1, CXC downgraded | NF-kB TF evidence most robust |

---

*This comparison was generated to support SME review of the pipeline v2 rerun. All raw data and scripts are version-controlled. v1 state is recoverable at commit `c950d1d`.*
