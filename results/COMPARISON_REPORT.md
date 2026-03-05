# Comparative Analysis: Claude Code Pipeline vs. Phylo Pipeline

**IVD Single-Cell RNA-seq Meta-Analysis**

---

## Executive Summary

Two independent analytical pipelines were applied to overlapping but non-identical collections of human intervertebral disc (IVD) scRNA-seq datasets. Despite major differences in dataset scope, integration strategy, cell type resolution, and DE methodology, the two analyses converge on several core biological conclusions while diverging in important ways that reveal methodological sensitivities. This document systematically compares the two approaches and interprets their concordances and discordances.

---

## 1. Dataset Scope and Scale

| Feature | Claude Code Pipeline | Phylo Pipeline |
|---------|---------------------|----------------|
| Datasets | 12 | 7 |
| Total cells (post-QC) | 436,239 | 173,628 |
| Donors | 57 | 29 |
| Samples | 78 | ~40 |
| Compartments | NP, AF, CEP | NP, AF, CEP |
| Unique datasets | GSE189916, GSE199866, GSE205535, CNP0002664, GSE251686 | (none unique) |

The Claude Code pipeline includes all 7 Phylo datasets plus 5 additional studies, more than doubling the cell count and donor pool. The added datasets contribute more NP-specific studies (including neonatal, aged, and herniated conditions), increasing statistical power for pseudobulk DE and broadening condition coverage.

**Implication:** The larger dataset in the Claude Code pipeline provides more statistical power for rare comparisons but also increases heterogeneity and the risk of batch effects, making integration strategy more critical.

---

## 2. Integration Strategy

| Feature | Claude Code Pipeline | Phylo Pipeline |
|---------|---------------------|----------------|
| Primary method | scANVI (semi-supervised) | Harmony |
| Alternatives tested | scVI, Harmony, BBKNN | (not reported) |
| Integration tiers | Yes — non-resident vs. resident cells | No — unified integration |
| Selection metric | scIB composite score (0.615) | Visual assessment |
| Cell type ASW | 0.511-0.521 | Not reported |

**Key difference:** The Claude Code pipeline uses a tiered integration strategy that separately processes non-resident cells (immune, endothelial; ~14.5K cells) from resident disc cells (NP, AF). This preserves the biological continuum among disc-resident populations that aggressive batch correction can erase. scANVI's semi-supervised approach also leverages marker-based annotations to guide integration.

Harmony, used by the Phylo pipeline, is a linear correction method that operates in PCA space. It is effective for removing batch effects but may over-correct subtle biological gradients, particularly the notochordal-to-degenerative continuum in NP cells.

**Hypothesis:** The tiered scANVI approach may better preserve cell state gradients, which could explain why the Claude Code pipeline's trajectory analysis shows stronger pseudotime-condition correlations and more trajectory-DE gene overlap (~55%).

---

## 3. Cell Type Annotation

| Feature | Claude Code Pipeline | Phylo Pipeline |
|---------|---------------------|----------------|
| NP subtypes | NP_notochordal, NP_mature_chondrocyte, NP_stressed_degenerative | NP_chondrocyte, NP_chondrocyte_HAPLN1, NP_degenerative_UPR, NP_stress_response, NP_metallothionein |
| AF subtypes | AF_inner, AF_outer | AF_fibroblast |
| Immune | CellTypist subtypes (Tcm/Naive, Tem/Trm, etc.) | T_NK_cell, Macrophage |
| Endothelial | Endothelial cells | Endothelial |
| Total clusters | ~10 types | 12 clusters |
| Method | Marker scoring + CellTypist | Leiden clustering + manual |

The two pipelines use different cell type ontologies, making direct comparison of DE results challenging. The Phylo pipeline resolves NP cells into 5 subtypes (including HAPLN1+, UPR, metallothionein, stress response clusters), while the Claude Code pipeline uses 3 broader NP categories. Conversely, the Claude Code pipeline distinguishes AF_inner from AF_outer, while Phylo treats AF as a single fibroblast population.

**Implication:** The finer NP resolution in the Phylo pipeline may capture specialized subpopulation responses (e.g., metallothionein-high cells as a distinct stress state), while the Claude Code pipeline's AF subdivision enables inner-vs-outer AF comparisons that are anatomically meaningful.

---

## 4. Differential Expression

### 4.1 Methodology

| Feature | Claude Code Pipeline | Phylo Pipeline |
|---------|---------------------|----------------|
| Tool | pyDESeq2 (Python) | DESeq2 (R) |
| Approach | Pseudobulk per sample | Pseudobulk per sample |
| Thresholds | \|log2FC\| > 0.5, padj < 0.05 | padj < 0.05 (no LFC cutoff stated) |
| Powered comparisons | 17 | ~5-6 (severe_vs_healthy per cluster) |
| Total DE genes | 5,328 | ~1,000+ (mostly downregulated) |

### 4.2 Direction of Change

A striking divergence: the Phylo pipeline reports a ~7:1 downregulated-to-upregulated ratio in severe degeneration, while the Claude Code pipeline shows a more balanced ratio with prominent upregulated inflammatory genes.

**Phylo top DE genes (NP_chondrocyte severe vs healthy):**
- LINC01578 (-28 LFC), H3F3B (-17.5 LFC), H3F3A (-15.6 LFC), H2AFZ (-15.4 LFC), HIST1H1C (-14.3 LFC)
- Dominated by **histone genes** and lncRNAs with extreme log2FC values

**Claude Code top DE genes (NP_mature_chondrocyte mild vs severe):**
- CXCL1 (+3.75), CXCL3 (+3.72), CXCL2 (+3.13), TNF (+2.45), MDK (+2.72)
- Classical **inflammatory/catabolic** IVD signature

### 4.3 Interpretation of DE Divergence

The Phylo pipeline's extreme log2FC values (up to -28) and histone gene dominance strongly suggest **technical confounding**:

1. **Histone gene artifacts:** Replication-dependent histone genes (H3F3B, H2AFZ, HIST1H1C, H4C variants) are known to be highly sensitive to cell cycle state, dissociation protocols, and ambient RNA contamination. Their extreme downregulation likely reflects differences in tissue processing or cell viability between healthy and degenerated samples, not true transcriptomic changes.

2. **Comparison design:** The Phylo pipeline compares healthy vs. severe directly, where study and condition are maximally confounded (healthy samples come primarily from different studies than severe). The Claude Code pipeline includes mild_vs_severe comparisons that are more often within-study, reducing confounding.

3. **LFC shrinkage:** The Claude Code pipeline's pyDESeq2 applies LFC shrinkage by default, constraining estimates to biologically plausible ranges. The Phylo pipeline's extreme LFCs suggest either no shrinkage or insufficient regularization.

**Hypothesis:** The Phylo pipeline's DE results are substantially driven by batch/processing artifacts (histone genes), while the Claude Code pipeline's more conservative approach surfaces genuine disease biology (inflammatory chemokines, ECM remodeling). The Claude Code pipeline explicitly flagged its own healthy_vs_herniated comparison (4,316 DE genes with ribosomal protein enrichment) as study-confounded, demonstrating similar artifact awareness.

---

## 5. Pathway Enrichment

### 5.1 Phylo Pipeline GSEA Results

The top ~30 enriched pathways in the Phylo pipeline (NP_chondrocyte, severe vs healthy) are almost entirely driven by histone genes:
- RNA Polymerase I Transcription (NES -1.96)
- Senescence-Associated Secretory Phenotype (NES -1.91)
- Cellular Senescence (NES -1.86)
- Beta-catenin/TCF transactivating complex (NES -1.90)
- DNA methylation, PRC2, HDAC deacetylation pathways

The **core enrichment** for virtually all of these pathways consists of the same set of histone genes (H4C15, H4C11, H2BC12, H2AC8, etc.). This means the enrichment results are not independent signals — they are a single histone depletion artifact propagating through dozens of pathway gene sets that contain histones.

The Phylo report interprets these as evidence for "suppression of Wnt signaling," "loss of Notch," and "senescence." While senescence-related chromatin changes are biologically plausible in degeneration, the overwhelming dominance of replication-dependent histones suggests the signal is primarily technical.

### 5.2 Claude Code Pipeline GSEA Results

The Claude Code pipeline's top enrichments include:
- **AF_inner (mild vs severe up):** Heat acclimation, protein refolding, unfolded protein response — driven by heat shock proteins (HSPA1A/B, HSPA6, DNAJB1)
- **AF_inner (down):** Oxidative phosphorylation, electron transport chain — mitochondrial dysfunction
- **NP/AF:** Granulocyte chemotaxis (CXCL1/2/3, CCL2), TNF signaling, cytokine-mediated signaling

These pathways are driven by diverse gene sets (not a single gene family), are biologically coherent with IVD degeneration (inflammation, stress response, metabolic decline), and show pathway-specific leading edges.

### 5.3 Convergent Pathway Findings

Despite the methodological differences, both analyses identify:
- **Senescence/stress response** as a feature of degeneration (though the Phylo signal is histone-contaminated)
- **Wnt/beta-catenin pathway** changes (Phylo: suppressed; Claude Code: not among top hits but trajectory-associated)
- **Inflammatory signaling** (more prominent in Claude Code; present but secondary in Phylo)

---

## 6. Cell-Cell Communication

| Feature | Claude Code Pipeline | Phylo Pipeline |
|---------|---------------------|----------------|
| Tool | LIANA (5-method consensus) | LIANA |
| Methods | CellPhoneDB, NATMI, Connectome, SingleCellSignalR, log2FC | Similar consensus |
| Healthy interactions | 44,079 | Not quantified |
| Degenerated interactions | 53,036 | Not quantified |

### 6.1 Convergent CCC Findings

**FN1 signaling gain:** Both pipelines identify increased FN1 (fibronectin) signaling in degeneration.
- Phylo: FN1->ITGA6, FN1->C5AR1, FN1->CD44 among top gained interactions
- Claude Code: FN1 interactions among increased degeneration signaling (53K vs 44K total interactions)

This is biologically robust — FN1 is a hallmark of fibrotic ECM remodeling in degenerative discs, and its interaction with integrins and complement receptors reflects both structural change and immune activation.

### 6.2 Phylo-Specific Finding: TIMP1->CD63 Loss

The Phylo pipeline's most striking CCC result is the **loss of TIMP1->CD63 signaling** across virtually all cell pair combinations (12 of 16 top lost interactions involve TIMP1->CD63).

TIMP1 (tissue inhibitor of metalloproteinases-1) binding to CD63 (a tetraspanin) promotes cell survival and inhibits apoptosis. Its loss in degeneration is consistent with increased MMP activity and cell death. This was not highlighted in the Claude Code pipeline, possibly because:
- Different cell type granularity affects which interactions are detected
- The Claude Code pipeline's per-dataset (not integrated) approach to CCC may dilute this signal
- Threshold differences in what constitutes "significant" change

### 6.3 Phylo-Specific Finding: SEMA4A->PLXNB1 Gain

Both NP chondrocyte subtypes show gained SEMA4A->PLXNB1 signaling. Semaphorins are nerve guidance molecules, linking to the pain biology hypothesis. This interaction was not prominently featured in the Claude Code pipeline's results but the Claude Code pipeline identified pain-relevant findings through a different route (TNF, CXCL8, neurotrophin pathways).

---

## 7. Trajectory and Cell State Analysis

| Feature | Claude Code Pipeline | Phylo Pipeline |
|---------|---------------------|----------------|
| Method | PAGA + DPT | Not performed |
| NP pseudotime-condition rho | -0.207 | N/A |
| AF pseudotime-condition rho | -0.177 | N/A |
| Trajectory genes | 500 per compartment | N/A |
| Trajectory-DE overlap | ~55% | N/A |

The Claude Code pipeline's trajectory analysis provides a unique contribution: evidence that NP cells follow a notochordal -> mature chondrocyte -> stressed/degenerative continuum, and that progression along this trajectory correlates with disease severity. The ~55% overlap between trajectory-associated genes and DE genes validates that the trajectory captures disease biology, not batch artifacts.

The Phylo pipeline did not perform trajectory analysis, missing this dimension of the biology.

---

## 8. Pain Biology

| Feature | Claude Code Pipeline | Phylo Pipeline |
|---------|---------------------|----------------|
| Dedicated analysis | Yes — curated pain gene sets | Partial — SEMA4A noted |
| Key pain mediators | TNF, CXCL8, CXCL1-3 | SEMA4A->PLXNB1 (nerve guidance) |
| Pain-relevant L-R pairs | 3,662-4,194 | Not quantified |
| Conclusion | Disc cells create pro-inflammatory environment promoting nerve ingrowth | Gained semaphorin signaling may affect innervation |

Both pipelines contribute complementary views of pain biology:
- **Claude Code:** Disc cells produce inflammatory mediators (TNF, CXCLs) that sensitize nerve endings, but do not express nociceptors themselves — consistent with the "inflammatory milieu" model
- **Phylo:** Gained SEMA4A->PLXNB1 signaling suggests active nerve guidance/repulsion changes in degeneration

Together, these suggest degeneration involves both inflammation-mediated nerve sensitization AND altered nerve guidance signaling.

---

## 9. Composition Analysis

Both pipelines attempted to test whether cell type proportions change with degeneration:
- **Claude Code:** No significant changes after FDR correction, though trends are biologically consistent
- **Phylo:** Only AF fibroblast showed marginal significance (Kruskal-Wallis p=0.038)

The lack of robust composition changes likely reflects high inter-donor variability and the confounding of cell type proportions with tissue sampling (biopsy location within the disc).

---

## 10. Summary of Key Concordances

| Finding | Claude Code | Phylo | Confidence |
|---------|:-----------:|:-----:|:----------:|
| FN1 signaling increases in degeneration | Yes | Yes | **High** |
| Inflammatory gene upregulation (CXCL, TNF) | Strong | Weak | **High** (Claude Code) |
| Stress response / UPR activation | Yes (HSPs) | Yes (UPR cluster) | **High** |
| Senescence features | TF activity (E2F4) | GSEA (histone-driven) | **Moderate** |
| No robust composition changes | Yes | Yes | **High** |
| Increased CCC complexity in degeneration | 53K vs 44K | Not quantified | **Moderate** |
| ECM remodeling (CEMIP, collagens) | Yes (AF) | Yes (FN1) | **High** |

---

## 11. Summary of Key Discordances

| Feature | Claude Code | Phylo | Likely Explanation |
|---------|-------------|-------|-------------------|
| DE gene direction ratio | Balanced up/down | 7:1 down:up | Phylo dominated by histone artifacts; missing LFC shrinkage |
| Top DE genes | CXCL1/2/3, TNF, MDK | LINC01578, H3F3B, histones | Phylo's healthy vs severe comparison is cross-study confounded |
| GSEA top pathways | HSP/inflammation/mitochondria | Histone/chromatin/senescence | Same histone artifact propagating through pathways |
| TIMP1->CD63 loss | Not highlighted | Dominant finding | Different CCC approaches and cell type resolution |
| Trajectory analysis | Strong disease correlation | Not performed | Methodological difference |
| Cell type resolution | Finer AF, coarser NP | Finer NP, coarser AF | Different clustering strategies |

---

## 12. Recommendations

1. **Histone gene filtering:** Future analyses should consider excluding replication-dependent histone genes from DE and GSEA analyses, or at minimum flagging their dominance. Their extreme sensitivity to technical factors makes them unreliable disease markers in cross-study comparisons.

2. **Within-study comparisons:** When possible, DE should prioritize within-study comparisons (e.g., mild vs. severe from the same study) to reduce batch confounding. Cross-study comparisons (healthy vs. severe) require careful batch modeling.

3. **TIMP1->CD63 validation:** The Phylo pipeline's TIMP1->CD63 finding is biologically compelling and should be investigated in the Claude Code pipeline's CCC results at a more granular level.

4. **SEMA4A pain axis:** The semaphorin-plexin finding from the Phylo pipeline complements the Claude Code pipeline's inflammatory pain mediator analysis and should be integrated into a unified pain biology model.

5. **Trajectory + DE integration:** The Claude Code pipeline's demonstration that ~55% of trajectory genes overlap with DE genes provides strong evidence that the NP cell state continuum is disease-relevant. This should be explored in the Phylo pipeline's finer NP subtype resolution.

---

## 13. Conclusion

The two pipelines converge on a core biological narrative: IVD degeneration involves ECM remodeling (FN1, CEMIP), inflammatory activation (CXCLs, TNF), stress responses (HSPs, UPR), and increased intercellular signaling complexity. They diverge primarily due to technical factors — the Phylo pipeline's DE and GSEA results are substantially contaminated by histone gene artifacts from cross-study confounding, while the Claude Code pipeline's more conservative approach (LFC shrinkage, within-study comparisons, study-confound flagging) surfaces more biologically interpretable signals.

The pipelines are complementary: the Phylo pipeline's finer NP subtype resolution and TIMP1->CD63/SEMA4A findings add biological depth, while the Claude Code pipeline's larger dataset, trajectory analysis, and artifact-aware methodology provide a more robust statistical framework. An ideal analysis would combine the Claude Code pipeline's methodological rigor with the Phylo pipeline's cell type granularity, while applying histone gene filtering and within-study DE prioritization to both.
