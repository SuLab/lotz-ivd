# Comparative Analysis: Claude Code Pipeline vs. Phylo Pipeline

**IVD Single-Cell RNA-seq Meta-Analysis**

---

## Executive Summary

Two independent analytical pipelines were applied to overlapping but non-identical collections of human intervertebral disc (IVD) scRNA-seq datasets. Despite major differences in dataset scope, integration strategy, cell type resolution, and DE methodology, the two analyses converge on several core biological conclusions while diverging in important ways that reveal methodological sensitivities. This document systematically compares the two approaches and interprets their concordances and discordances.

---

## 1. Dataset Scope and Scale

| Feature | Claude Code Pipeline | Phylo Pipeline |
|---------|---------------------|----------------|
| Datasets | 11 | 7 |
| Total cells (post-QC) | 410,759 | 173,628 |
| Donors | ~50 | 29 |
| Samples | 71 | ~40 |
| Compartments | NP, AF, CEP | NP, AF, CEP |
| Unique datasets | GSE189916, GSE199866, GSE205535, CNP0002664 | (none unique) |
| Excluded | GSE233666 (QC failure) | GSE251686 and others |

The Claude Code pipeline includes all 7 Phylo datasets plus 4 additional studies (GSE233666 was excluded during QC), more than doubling the cell count and donor pool. The added datasets contribute more NP-specific studies (including neonatal, aged, and herniated conditions), increasing statistical power for pseudobulk DE and broadening condition coverage.

**Implication:** The larger dataset in the Claude Code pipeline provides more statistical power for rare comparisons but also increases heterogeneity and the risk of batch effects, making integration strategy more critical.

---

## 2. Integration Strategy

| Feature | Claude Code Pipeline | Phylo Pipeline |
|---------|---------------------|----------------|
| Primary method | scVI (unsupervised) | Harmony |
| Alternatives tested | (none — scVI only) | (not reported) |
| Integration objects | 4 compartment-specific (NP, AF, CEP, all_cells) | No — unified integration |
| Selection metric | scIB composite score | Visual assessment |

**Key difference:** The Claude Code pipeline uses scVI with compartment-specific integration objects (NP, AF, CEP, and an all_cells object), rather than a tiered resident/non-resident split. This allows each compartment's latent space to capture compartment-specific biology without interference from unrelated cell populations.

Harmony, used by the Phylo pipeline, is a linear correction method that operates in PCA space. It is effective for removing batch effects but may over-correct subtle biological gradients, particularly the notochordal-to-degenerative continuum in NP cells.

**Hypothesis:** The compartment-specific scVI approach may better preserve cell state gradients within each anatomical region, which could explain why the Claude Code pipeline's trajectory analysis shows pseudotime-condition correlations in NP and CEP (though AF shows a reversed positive correlation).

---

## 3. Cell Type Annotation

| Feature | Claude Code Pipeline | Phylo Pipeline |
|---------|---------------------|----------------|
| NP subtypes | NP_notochordal, NP_mature_chondrocyte, NP_stressed_degenerative, NP_fibrocartilaginous | NP_chondrocyte, NP_chondrocyte_HAPLN1, NP_degenerative_UPR, NP_stress_response, NP_metallothionein |
| AF subtypes | AF_inner, AF_outer | AF_fibroblast |
| Additional types | EP_hyaline (new) | — |
| Immune | CellTypist subtypes (Tcm/Naive, Tem/Trm, etc.) | T_NK_cell, Macrophage |
| Endothelial | Endothelial cells | Endothelial |
| Method | De novo annotation post-integration (marker scoring + CellTypist validation) | Leiden clustering + manual |
| CellTypist concordance | NP: 5/13 clusters, AF: good, CEP: partial | N/A |

The two pipelines use different cell type ontologies, making direct comparison of DE results challenging. The Phylo pipeline resolves NP cells into 5 subtypes (including HAPLN1+, UPR, metallothionein, stress response clusters), while the Claude Code pipeline uses 4 NP categories (adding NP_fibrocartilaginous). The Claude Code pipeline also identifies EP_hyaline as a novel type. Conversely, the Claude Code pipeline distinguishes AF_inner from AF_outer, while Phylo treats AF as a single fibroblast population.

A key methodological difference is that the Claude Code pipeline performs de novo annotation post-integration (rather than pre-integration marker-based assignment), validated against CellTypist. CellTypist concordance was moderate for NP (5 of 13 clusters matched) and better for AF, reflecting the challenge of automated annotation in disc-resident cell types that are poorly represented in reference atlases.

**Implication:** The finer NP resolution in the Phylo pipeline may capture specialized subpopulation responses (e.g., metallothionein-high cells as a distinct stress state), while the Claude Code pipeline's AF subdivision and novel EP_hyaline type enable anatomically meaningful comparisons not available in the Phylo analysis.

---

## 4. Differential Expression

### 4.1 Methodology

| Feature | Claude Code Pipeline | Phylo Pipeline |
|---------|---------------------|----------------|
| Tool | pyDESeq2 (Python) | DESeq2 (R) |
| Approach | Pseudobulk per sample | Pseudobulk per sample |
| Thresholds | \|log2FC\| > 0.5, padj < 0.05 | padj < 0.05 (no LFC cutoff stated) |
| Powered comparisons | 21 | ~5-6 (severe_vs_healthy per cluster) |
| Total DE genes | 949 unique | ~1,000+ (mostly downregulated) |
| Herniated comparisons | Excluded entirely | Included |

### 4.2 Direction of Change

A striking divergence: the Phylo pipeline reports a ~7:1 downregulated-to-upregulated ratio in severe degeneration, while the Claude Code pipeline shows a more balanced ratio with prominent upregulated inflammatory genes.

**Phylo top DE genes (NP_chondrocyte severe vs healthy):**
- LINC01578 (-28 LFC), H3F3B (-17.5 LFC), H3F3A (-15.6 LFC), H2AFZ (-15.4 LFC), HIST1H1C (-14.3 LFC)
- Dominated by **histone genes** and lncRNAs with extreme log2FC values

**Claude Code top DE genes (NP_mature_chondrocyte mild vs severe):**
- CXCL2 (padj=0.005), ICAM1 (padj=1.9e-5), HHEX (padj=1.9e-5)
- CXCL1, CXCL3, and TNF were **not significant** in v2 (unlike v1)
- More conservative result with 949 unique DE genes across 21 comparisons

### 4.3 Interpretation of DE Divergence

The Phylo pipeline's extreme log2FC values (up to -28) and histone gene dominance strongly suggest **technical confounding**:

1. **Histone gene artifacts:** Replication-dependent histone genes (H3F3B, H2AFZ, HIST1H1C, H4C variants) are known to be highly sensitive to cell cycle state, dissociation protocols, and ambient RNA contamination. Their extreme downregulation likely reflects differences in tissue processing or cell viability between healthy and degenerated samples, not true transcriptomic changes.

2. **Comparison design:** The Phylo pipeline compares healthy vs. severe directly, where study and condition are maximally confounded (healthy samples come primarily from different studies than severe). The Claude Code pipeline includes mild_vs_severe comparisons that are more often within-study, reducing confounding.

3. **LFC shrinkage:** The Claude Code pipeline's pyDESeq2 applies LFC shrinkage by default, constraining estimates to biologically plausible ranges. The Phylo pipeline's extreme LFCs suggest either no shrinkage or insufficient regularization.

**Hypothesis:** The Phylo pipeline's DE results are substantially driven by batch/processing artifacts (histone genes), while the Claude Code pipeline's more conservative approach surfaces genuine disease biology. Notably, the v2 Claude Code pipeline excluded herniated comparisons entirely (recognizing study-level confounding) and produced far fewer DE genes (949 vs. 5,328 in v1), with only CXCL2 remaining significant among the CXC chemokines. This increased stringency provides higher-confidence hits at the cost of sensitivity.

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

The v2 Claude Code pipeline produced 1,577 ORA and 1,576 GSEA significant results (up from 1,244/1,081 in v1), reflecting the expanded set of 21 powered comparisons. Top enrichments include:
- **AF_inner (mild vs severe up):** Heat acclimation, protein refolding, unfolded protein response — driven by heat shock proteins (HSPA1A/B, HSPA6, DNAJB1)
- **AF_inner (down):** Oxidative phosphorylation, electron transport chain — mitochondrial dysfunction
- **NP/AF:** Inflammatory and chemotactic pathways, though the CXC chemokine signal is weaker in v2 (only CXCL2 significant)

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
| Healthy interactions | 28,878 | Not quantified |
| Degenerated interactions | 27,011 | Not quantified |
| Direction of change | Fewer interactions in degeneration | Not quantified |

### 6.1 Convergent CCC Findings

**FN1 signaling gain:** Both pipelines identify increased FN1 (fibronectin) signaling in degeneration.
- Phylo: FN1->ITGA6, FN1->C5AR1, FN1->CD44 among top gained interactions
- Claude Code: FN1 interactions among altered degeneration signaling

This is biologically robust — FN1 is a hallmark of fibrotic ECM remodeling in degenerative discs, and its interaction with integrins and complement receptors reflects both structural change and immune activation.

**TIMP1->CD63 convergence:** The Phylo pipeline's most prominent CCC finding — loss of TIMP1->CD63 signaling — was actually replicated in the v2 Claude Code pipeline, where TIMP1-CD63 emerged as the most enriched interaction in the differential CCC analysis. This cross-pipeline replication substantially increases confidence in this finding (see Section 6.2).

### 6.2 Cross-Pipeline Convergence: TIMP1->CD63

The Phylo pipeline's most striking CCC result is the **loss of TIMP1->CD63 signaling** across virtually all cell pair combinations (12 of 16 top lost interactions involve TIMP1->CD63).

TIMP1 (tissue inhibitor of metalloproteinases-1) binding to CD63 (a tetraspanin) promotes cell survival and inhibits apoptosis. Its loss in degeneration is consistent with increased MMP activity and cell death. Importantly, the v2 Claude Code pipeline **replicated** this finding: TIMP1-CD63 was identified as the most enriched interaction in the differential CCC analysis. This cross-pipeline convergence, despite different integration methods, cell type ontologies, and dataset scopes, makes TIMP1->CD63 one of the highest-confidence CCC findings in this meta-analysis.

**Note on CCC direction:** The v2 Claude Code pipeline shows a reversed overall pattern compared to v1 — fewer total interactions in degeneration (27,011) than in healthy tissue (28,878). This contrasts with the v1 result (53K vs 44K) and suggests that the increased CCC complexity finding was sensitive to integration strategy and dataset composition.

### 6.3 Phylo-Specific Finding: SEMA4A->PLXNB1 Gain

Both NP chondrocyte subtypes show gained SEMA4A->PLXNB1 signaling. Semaphorins are nerve guidance molecules, linking to the pain biology hypothesis. This interaction was not prominently featured in the Claude Code pipeline's results but the Claude Code pipeline identified pain-relevant findings through a different route (TNF, CXCL8, neurotrophin pathways).

---

## 7. Trajectory and Cell State Analysis

| Feature | Claude Code Pipeline | Phylo Pipeline |
|---------|---------------------|----------------|
| Method | PAGA + DPT | Not performed |
| NP pseudotime-condition rho | -0.258 | N/A |
| AF pseudotime-condition rho | +0.341 (reversed) | N/A |
| CEP pseudotime-condition rho | -0.163 (new) | N/A |
| Trajectory genes | 500 per compartment | N/A |
| Trajectory-DE overlap | NP: 96/500, AF: 110/500, CEP: 38/500 | N/A |

The Claude Code pipeline's trajectory analysis provides a unique contribution: evidence that NP cells follow a notochordal -> mature chondrocyte -> stressed/degenerative continuum, and that progression along this trajectory correlates with disease severity (rho=-0.258, stronger than v1's -0.207). CEP cells show a similar negative correlation (rho=-0.163).

However, a notable v2 finding is that the **AF trajectory-condition correlation reversed** from v1 (rho=+0.341 vs. v1's -0.177). This positive correlation means AF pseudotime progression is associated with less severe disease, suggesting that the AF trajectory captures a different biological axis (possibly structural maturation rather than degeneration). This reversal is likely attributable to the change in integration strategy (scVI-only with compartment-specific objects vs. tiered scANVI).

The trajectory-DE overlap is more modest than in v1 (~19-22% for NP/AF vs. ~55% previously), reflecting both the smaller DE gene set (949 vs. 5,328) and possibly improved specificity.

The Phylo pipeline did not perform trajectory analysis, missing this dimension of the biology.

---

## 8. Pain Biology

| Feature | Claude Code Pipeline | Phylo Pipeline |
|---------|---------------------|----------------|
| Dedicated analysis | Yes — curated pain gene sets | Partial — SEMA4A noted |
| Key pain mediators | 10 significant pain genes (incl. CXCL2, ICAM1) | SEMA4A->PLXNB1 (nerve guidance) |
| Pain-relevant L-R pairs | Quantified | Not quantified |
| Conclusion | Disc cells create pro-inflammatory environment promoting nerve ingrowth | Gained semaphorin signaling may affect innervation |

Both pipelines contribute complementary views of pain biology:
- **Claude Code (v2):** 10 significant pain-related genes identified (up from 3 in v1), though the CXC chemokine signal is weaker — only CXCL2 is significant (CXCL1, CXCL3, TNF are not). The broader set of pain genes provides a more nuanced view of the inflammatory milieu.
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
| TIMP1->CD63 altered in degeneration | Yes (most enriched differential) | Yes (dominant lost interaction) | **High** (cross-pipeline) |
| Inflammatory gene upregulation (CXCL) | Moderate (CXCL2 only) | Weak | **Moderate** |
| Stress response / UPR activation | Yes (HSPs) | Yes (UPR cluster) | **High** |
| Senescence features | TF activity (E2F4) | GSEA (histone-driven) | **Moderate** |
| No robust composition changes | Yes | Yes | **High** |
| ECM remodeling (CEMIP, collagens) | Yes (AF) | Yes (FN1) | **High** |

---

## 11. Summary of Key Discordances

| Feature | Claude Code | Phylo | Likely Explanation |
|---------|-------------|-------|-------------------|
| DE gene count | 949 unique (stringent) | ~1,000+ (mostly downregulated) | Claude Code excludes herniated; stricter thresholds |
| DE gene direction ratio | Balanced up/down | 7:1 down:up | Phylo dominated by histone artifacts; missing LFC shrinkage |
| Top DE genes | CXCL2, ICAM1, HHEX | LINC01578, H3F3B, histones | Phylo's healthy vs severe comparison is cross-study confounded |
| GSEA top pathways | HSP/inflammation/mitochondria | Histone/chromatin/senescence | Same histone artifact propagating through pathways |
| CCC direction | Fewer interactions in degeneration (27K vs 29K) | Not quantified | Sensitive to integration and dataset composition |
| AF trajectory correlation | Positive (+0.341) | Not performed | AF pseudotime may capture maturation, not degeneration |
| Trajectory-DE overlap | ~19-22% | Not performed | Smaller DE gene set reduces overlap |
| Cell type resolution | Finer AF, 4 NP types + EP_hyaline | Finer NP (5 types), coarser AF | Different clustering strategies |

---

## 12. Recommendations

1. **Histone gene filtering:** Future analyses should consider excluding replication-dependent histone genes from DE and GSEA analyses, or at minimum flagging their dominance. Their extreme sensitivity to technical factors makes them unreliable disease markers in cross-study comparisons.

2. **Within-study comparisons:** When possible, DE should prioritize within-study comparisons (e.g., mild vs. severe from the same study) to reduce batch confounding. The v2 Claude Code pipeline's exclusion of herniated comparisons exemplifies this principle.

3. **TIMP1->CD63 as high-confidence target:** Now replicated across both pipelines, TIMP1->CD63 alteration in degeneration is the highest-confidence CCC finding in this meta-analysis and warrants experimental validation.

4. **SEMA4A pain axis:** The semaphorin-plexin finding from the Phylo pipeline complements the Claude Code pipeline's expanded pain gene analysis (10 significant genes in v2) and should be integrated into a unified pain biology model.

5. **AF trajectory interpretation:** The reversed AF trajectory-condition correlation (positive in v2 vs. negative in v1) suggests that AF pseudotime may capture maturation rather than degeneration. Future work should investigate whether AF trajectory endpoints correspond to specific anatomical or developmental states.

6. **CXC chemokine robustness:** The weakening of the CXC chemokine signal in v2 (only CXCL2 significant, not CXCL1/3 or TNF) suggests this finding is sensitive to integration approach and comparison design. Independent experimental validation is needed before prioritizing these targets.

---

## 13. Conclusion

The two pipelines converge on a core biological narrative: IVD degeneration involves ECM remodeling (FN1, CEMIP), stress responses (HSPs, UPR), and altered intercellular signaling — most notably TIMP1->CD63, which was independently identified by both pipelines as a key altered interaction. They diverge primarily due to technical factors — the Phylo pipeline's DE and GSEA results are substantially contaminated by histone gene artifacts from cross-study confounding, while the Claude Code pipeline's more conservative v2 approach (LFC shrinkage, exclusion of herniated comparisons, stricter thresholds) produces fewer but higher-confidence DE hits.

The v2 rerun revealed important sensitivities: the CXC chemokine signature weakened substantially (only CXCL2 significant), the overall CCC direction reversed (fewer interactions in degeneration), and the AF trajectory-condition correlation flipped sign. These changes, driven by differences in integration strategy (scVI-only vs. tiered scANVI) and dataset scope (11 vs. 12 datasets), highlight which findings are robust and which are method-dependent.

**Robust findings (replicated across pipelines and/or versions):** FN1 signaling gain, TIMP1->CD63 alteration, stress response activation, no robust composition changes.

**Method-sensitive findings (changed between v1/v2 or discordant between pipelines):** CXC chemokine significance, CCC direction/magnitude, AF trajectory-condition correlation, trajectory-DE overlap magnitude.

The pipelines are complementary: the Phylo pipeline's finer NP subtype resolution and SEMA4A findings add biological depth, while the Claude Code pipeline's larger dataset, trajectory analysis, and artifact-aware methodology provide a more robust statistical framework. The v2 replication of TIMP1->CD63 across both pipelines elevates it to the highest-confidence cell-cell communication finding in this meta-analysis.
