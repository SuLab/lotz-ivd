# Human IVD scRNA-seq Multi-Dataset Integration Analysis

**173,628 cells | 7 datasets | 29 donors | 4 degeneration grades**

---

## 1. Overview

This report summarizes an integrated single-cell RNA-seq analysis of human intervertebral disc (IVD) tissue across 7 publicly available GEO datasets, spanning healthy donors and patients with mild, moderate, and severe disc degeneration. The goal was to characterize IVD cell types, identify how cell states and gene expression change with degeneration, and map altered signaling networks.

---

## 2. Datasets

| GEO Accession | Study | Tissue | Condition | Samples | Cells (post-QC) |
|---|---|---|---|---|---|
| GSE160756 | Gan/Liu 2021 | NP, AF, CEP | Healthy atlas | 7 | ~60,000 |
| GSE199866 | Cherif 2022 | NP | Paired degen/non-degen | 4 | ~8,000 |
| GSE244889 | Wang 2023 | NP | Mild vs severe | 7 | ~45,000 |
| GSE255768 | Shi 2024 | CEP | Degeneration | 2 | ~5,000 |
| GSE233666 | Guo 2023 | NP | Immune/ossification | 4 | ~20,000 |
| GSE205535 | Li 2022 | NP | Normal + degenerated | 2 | ~10,000 |
| GSE189916 | Jiang/Sheyn 2022 | IVD | Neonatal + adult | 6 | ~25,000 |

**QC:** MAD-based filtering per sample (nGenes, nCounts, %mito); Scrublet doublet removal. 173,628 / 222,433 cells retained (78%).

---

## 3. Integration

- **Batch correction:** Harmony on dataset + donor_id (29 donors, 7 datasets)
- **Neighbor graph:** FAISS approximate kNN (k=30, IVFFlat index)
- **Clustering:** Leiden at res=0.3/0.5/0.8/1.2 (9/12/19/27 clusters)
- **Gene space:** 25,304 genes present in ≥3 datasets; 4,000 batch-aware HVGs

---

## 4. Cell Type Annotation

12 clusters identified at Leiden res=0.5, manually annotated using IVD-specific marker panels:

| Cell Type | Cluster | Key Markers | n cells | % |
|---|---|---|---|---|
| NP: canonical | 0 | ACAN, COL2A1, SOX9, COL9A3, SCRG1 | 47,120 | 27.1% |
| NP: degenerative (UPR) | 3 | SQSTM1, DNAJB9, TNFRSF12A, EIF1 | 34,810 | 20.1% |
| NP: MT-high | 5 | MT1G, MT1E, MT1X, MT2A, MALAT1 | 31,502 | 18.1% |
| NP: HAPLN1+ | 1 | HAPLN1, FN1, TIMP3, FGFBP2 | 17,040 | 9.8% |
| AF fibroblast | 4 | COL1A1, COL1A2, COL3A1, SCX | 16,448 | 9.5% |
| NP: stress response | 2 | JUN, FOS, GADD45B, DNAJB1 | 9,967 | 5.7% |
| Pericyte/SMC | 7 | TAGLN, MYL9, CALD1, NR2F2 | 4,068 | 2.3% |
| Macrophage | 10 | CD68, CD14, TYROBP, CTSS | 3,482 | 2.0% |
| Endothelial | 8 | PECAM1, VWF, CDH5, SPARCL1 | 3,382 | 1.9% |
| Monocyte/Neutrophil | 11 | LYZ, S100A8, S100A9, MNDA | 2,511 | 1.4% |
| Erythrocyte | 9 | HBB, HBA1, HBA2, AHSP | 1,658 | 1.0% |
| T/NK cell | 6 | CD3D, CXCR4, GNLY, NKG7 | 1,640 | 0.9% |

**Key figures:** `05_annotation/umap_annotated_v2.png`, `05_annotation/dotplot_markers.png`

---

## 5. Compositional Analysis

Cell type proportions computed per donor (n=29), tested with Kruskal-Wallis across 4 conditions.

- **AF fibroblast** is the only cell type with significant compositional change (p=0.038), increasing from ~7% (healthy) to ~24% (severe) — consistent with fibrocartilaginous metaplasia.
- NP state proportions show high inter-donor variability, limiting statistical power with n=29 donors.
- Within the NP lineage: **NP: HAPLN1+** expands +10.1% in severe; **NP: MT-high** contracts −12.5%.

**Key figure:** `06_composition/composition_overview.png`

---

## 6. Pseudobulk Differential Expression (DESeq2)

Pseudobulk aggregation per donor × cell type; DESeq2 on 6 cell types × 3 contrasts (severe/moderate/mild vs healthy).

### Severe vs Healthy — DEG summary

| Cell Type | Total DEGs | Up | Down |
|---|---|---|---|
| NP: canonical | 2,641 | 332 | 2,309 |
| NP: degenerative | 2,331 | 241 | 2,090 |
| NP: HAPLN1+ | 1,786 | 223 | 1,563 |
| NP: stress | 1,081 | 104 | 977 |
| AF fibroblast | 1,298 | 103 | 1,195 |
| NP: MT-high | 953 | 81 | 872 |

**Dominant pattern:** Massive transcriptional downregulation in severe degeneration across all cell types (~7:1 down:up ratio). Key upregulated genes: ADAMTS5, FN1, TNF, CXCL8, IL6. Key downregulated: ACAN, COL2A1, HAPLN1, COMP, CILP.

**Key figure:** `07_pseudobulk/volcano_severe_vs_healthy.png`  
**Data:** `07_pseudobulk/all_DEGs_severe_vs_healthy.csv`

---

## 7. Pathway Enrichment (GSEA)

GSEA on ranked gene lists (sign(LFC) × −log10(padj)) using MSigDB Hallmarks + Reactome (IVD-relevant subset). 70–85 significant pathways per cell type.

### Consistently downregulated in severe degeneration (all/most cell types)
- **Wnt signaling** (TCF-dependent Wnt, Signaling by Wnt) — loss of chondrocyte maintenance
- **Notch signaling** (Signaling by Notch, Pre-Notch processing) — loss of progenitor niche
- **Cellular senescence** (Oxidative stress-induced senescence, SASP, DNA damage senescence)
- **ECM organization** (Extracellular Matrix Organization, Collagen Formation)
- **RUNX1/2 transcription** — loss of chondrogenic transcription factor activity

### Upregulated in severe degeneration (selected cell types)
- **TNFα signaling via NF-κB** (NP: canonical, NP: degenerative) — inflammatory activation
- **Inflammatory response** (NP: stress, NP: degenerative)
- **Epithelial-mesenchymal transition** (NP: degenerative) — fibrotic shift
- **Collagen crosslinking** (NP: stress, NP: degenerative, NP: MT-high) — matrix stiffening
- **Glycolysis** (NP: canonical) — metabolic reprogramming under hypoxia

**Key figure:** `08_pathways/gsea_heatmap_severe_vs_healthy.png`

---

## 8. NP Lineage Trajectory (PAGA)

PAGA connectivity analysis on 140,439 NP lineage cells.

**PAGA connectivity (key edges):**
- NP: stress ↔ NP: degenerative: **0.76** (strongest connection)
- NP: canonical ↔ NP: HAPLN1+: **0.65**
- NP: stress ↔ NP: MT-high: **0.51**

**Interpretation:** The NP stress-response state acts as a transition hub between canonical NP chondrocytes and the degenerative UPR state. The HAPLN1+ subtype (matrix-organizing) expands in severe degeneration, possibly representing a compensatory remodeling response.

**Note:** Diffusion pseudotime (DPT) was weakly resolved in this dataset (all DC variances = 7.12e-06), likely due to the large, well-mixed dataset compressing diffusion kernel eigenvalues. Spearman correlation of pseudotime with degeneration grade was significant but weak (ρ=0.24, p<10⁻³⁰⁰). PAGA connectivity is the more reliable trajectory metric here.

**Key figure:** `09_trajectory/np_trajectory_overview.png`

---

## 9. Cell-Cell Communication (LIANA)

LIANA rank_aggregate (consensus resource) on healthy (n=99,883 cells) vs severe (n=19,545 cells) subsets. 39,530 unique ligand-receptor pairs evaluated.

### Top lost interactions in severe degeneration
- **TIMP1 → CD63**: Dominant lost interaction across all cell type pairs. TIMP1 (tissue inhibitor of metalloproteinases) signals through CD63 to suppress MMP activity. Loss of this interaction in severe degeneration is consistent with the known increase in matrix degradation.

### Top gained interactions in severe degeneration
- **FN1 → CD44/C5AR1**: Fibronectin signaling to macrophages (via CD44 and complement receptor C5AR1) — consistent with inflammatory macrophage recruitment and activation.
- **FN1 → ITGA6** (NP:HAPLN1+ → Endothelial): Fibronectin-integrin signaling promoting angiogenesis.
- **COL1A2 → CD93** (AF fibroblast → Endothelial): Collagen-endothelial interaction, potentially driving vascular ingrowth.
- **SEMA4A → PLXNB1**: Semaphorin signaling gained in NP cells — role in axonal/vascular ingrowth.

**Overall pattern:** Global increase in total interaction strength in severe degeneration (all NP states show increased outgoing signaling), driven by FN1, inflammatory cytokines, and matrix-degradation signals. Protective TIMP1 signaling is lost.

**Key figure:** `10_cellchat/liana_communication_overview.png`

---

## 10. Key Biological Conclusions

1. **Five NP cell states** coexist in the IVD, representing a spectrum from canonical chondrocytes to stress-response, UPR-driven degenerative, and metallothionein-high oxidative stress states.

2. **Degeneration is characterized by transcriptional collapse**: 7:1 down:up DEG ratio in severe vs healthy, with loss of ECM maintenance genes (ACAN, COL2A1, HAPLN1, COMP) and chondrogenic transcription factors.

3. **Four pathways are consistently suppressed across all NP states**: Wnt signaling, Notch signaling, cellular senescence programs, and RUNX transcription — suggesting a coordinated loss of tissue homeostasis machinery.

4. **Inflammatory activation is cell-state specific**: TNF/NF-κB and EMT upregulation is strongest in NP: canonical and NP: degenerative states, not uniformly across all NP cells.

5. **AF fibroblast expansion** is the most robust compositional change (p=0.038), increasing from 7% to 24% in severe degeneration — consistent with fibrocartilaginous replacement of NP tissue.

6. **TIMP1→CD63 loss** is the dominant signaling change in degeneration, implicating reduced MMP inhibition as a key driver of matrix degradation. **FN1→macrophage** signaling gain suggests inflammatory amplification via fibronectin fragments.

---

## 11. Limitations

- **Unbalanced condition groups**: 100,234 healthy vs 19,610 severe cells; 15 healthy donors vs 3 severe donors limits statistical power for compositional and DE analyses.
- **DPT trajectory**: Diffusion pseudotime was weakly resolved; PAGA connectivity is used as the primary trajectory metric.
- **Annotation resolution**: At Leiden res=0.5, some biologically distinct subtypes (e.g., notochordal cells, CEP chondrocytes) may be merged. Higher resolution clustering (res=1.2, 27 clusters) available for follow-up.
- **LIANA without permutation testing**: n_perms=None used for computational efficiency; p-values are not permutation-based. Results should be interpreted as exploratory.
- **Batch effects**: Despite Harmony correction, dataset-of-origin effects may persist, particularly for the neonatal samples in GSE189916.

---

## 12. Output Files

| Folder | Key Files | Description |
|---|---|---|
| `05_annotation/` | `umap_annotated_v2.png`, `dotplot_markers.png`, `cell_metadata.csv` | Cell type annotation |
| `06_composition/` | `composition_overview.png`, `kruskal_wallis_results.csv` | Compositional analysis |
| `07_pseudobulk/` | `volcano_severe_vs_healthy.png`, `all_DEGs_severe_vs_healthy.csv` | DESeq2 results |
| `08_pathways/` | `gsea_heatmap_severe_vs_healthy.png`, `gsea_*.csv` | GSEA results |
| `09_trajectory/` | `np_trajectory_overview.png` | PAGA trajectory |
| `10_cellchat/` | `liana_communication_overview.png`, `liana_*.csv` | LIANA CCC results |

---

*Analysis performed with: scanpy 1.9, Harmony, FAISS, DESeq2, clusterProfiler, MSigDB, LIANA 1.7.1*  
*Datasets: GSE160756, GSE199866, GSE244889, GSE255768, GSE233666, GSE205535, GSE189916*
