# Module 04: Coarse Cell Classification for Integration Anchors

## Objective

Assign each cell a coarse label from 5 anchor categories — chondrocyte-like, fibroblast-like, immune, endothelial, pericyte/SMC — or "Unknown" if the cell doesn't clearly fit any category. These coarse labels serve as seed labels for scANVI semi-supervised integration in Module 05. Fine-grained cell type annotation happens *after* integration (in Module 05), not here.

## Rationale

The goal of this project is to build a de novo IVD cell atlas from integrated data, without relying on original manuscript annotations. Fine annotation before integration risks circular logic: labels shape integration, which shapes downstream analysis.

However, purely unsupervised integration (scVI) can fail to adequately correct batch effects when datasets come from different platforms and tissue processing protocols. scANVI (semi-supervised) uses labeled cells as anchors to align similar cells across studies, producing better integration.

The key insight is that coarse labels are reliable across datasets — a PTPRC+ cell is immune regardless of study, a COL2A1+/ACAN+ cell is chondrocyte-like regardless of platform. By labeling only the confident, unambiguous cells and leaving uncertain cells as "Unknown," we give scANVI enough structure to integrate well without imposing fine-grained distinctions that should emerge from the integrated data.

Cells on the mesenchymal continuum (e.g., fibrochondrocytes that co-express both chondrocyte and fibroblast markers) are labeled "Unknown" — scANVI will position them in the latent space based on their transcriptomic similarity to the labeled anchors at either end of the continuum.

## Inputs

- Processed AnnData objects from `data/processed/{accession}.h5ad`
- **Excluded from all analysis:** GSE233666 (Guo 2023) — herniated discs only, not the focus of this project
- **Excluded samples:** Neonatal samples from GSE189916 (Jiang 2022) — neonatal disc biology is fundamentally different from adult; only adult samples are retained

## Outputs

Per dataset:
- Updated `.obs` in the h5ad file with:
  - `coarse_label` — one of: "Chondrocyte_like", "Fibroblast_like", "Immune", "Endothelial", "Pericyte_SMC", "Unknown"
  - `cell_class` — derived from coarse_label: "mesenchymal" (Chondrocyte_like, Fibroblast_like), "non_mesenchymal" (Immune, Endothelial, Pericyte_SMC), or "unknown"
- `results/annotations/{accession}_classification_report.html`

Aggregate:
- `results/annotations/classification_summary.tsv` — coarse label proportions per dataset

### Notebook: `notebooks/04_classification.ipynb`

Contains:
- Per-dataset UMAPs colored by coarse_label
- Dot plot of all classification markers by coarse_label
- Bar plot: coarse label proportions per dataset
- Validation: marker expression distributions confirming clean separation of anchor categories
- Proportion of "Unknown" cells per dataset

## Anchor Categories and Marker Genes

### Category 1: Chondrocyte-like

Cells with clearly dominant chondrocyte markers.

**Positive markers:** COL2A1, ACAN, SOX9
**Requirement:** Chondrocyte score clearly higher than fibroblast score (e.g., chondrocyte score > 2× fibroblast score)

### Category 2: Fibroblast-like

Cells with clearly dominant fibroblast markers.

**Positive markers:** COL1A1, COL1A2, DCN, LUM
**Requirement:** Fibroblast score clearly higher than chondrocyte score (e.g., fibroblast score > 2× chondrocyte score)

### Category 3: Immune

**Positive markers:** PTPRC (CD45) — primary gate
**Supporting markers:** CD3D, CD3E (T cells), CD68, CD14, CSF1R (macrophage/monocyte), CD79A, MS4A1 (B cells), KIT, TPSAB1 (mast cells), NKG7, GNLY (NK cells)
**Requirement:** PTPRC expression above background, or co-expression of ≥2 supporting immune markers

### Category 4: Endothelial

**Positive markers:** PECAM1 (CD31), VWF, CDH5
**Requirement:** Expression of ≥1 endothelial marker above background, without co-expression of mesenchymal markers (COL2A1, ACAN)

### Category 5: Pericyte/SMC

**Positive markers:** RGS5, PDGFRB
**Requirement:** Co-expression of pericyte markers without dominant mesenchymal marker expression (ACTA2 alone is insufficient — degenerated disc cells can express ACTA2)

### Unknown

Cells that don't clearly fit any anchor category. This includes:
- Fibrochondrocytes and other transitional mesenchymal cells (both chondrocyte and fibroblast scores moderate, neither clearly dominant)
- Cells with low/ambiguous marker expression
- Any cell where classification confidence is low

**The "Unknown" fraction should be generous** — better to under-label and let scANVI infer than to force a wrong anchor label. Expect 20-30% of mesenchymal cells to be Unknown. This is not a problem — scANVI is designed for semi-supervised learning with partial labels.

## Classification Logic

### Per-cell scoring

1. Compute chondrocyte score: `sc.tl.score_genes()` with [COL2A1, ACAN, SOX9]
2. Compute fibroblast score: `sc.tl.score_genes()` with [COL1A1, COL1A2, DCN, LUM]
3. Check non-mesenchymal markers individually (PTPRC, PECAM1, VWF, CDH5, RGS5, PDGFRB)

### Assignment hierarchy

Apply in order (first match wins):

1. **Immune:** PTPRC > threshold, OR ≥2 immune supporting markers expressed. Must NOT co-express ACAN or SOX9 (rescue rule — prevents stressed disc cells with upregulated HLA genes from being misclassified).
2. **Endothelial:** PECAM1, VWF, or CDH5 > threshold. Must NOT co-express ACAN or SOX9.
3. **Pericyte/SMC:** RGS5 and PDGFRB co-expressed. Must NOT co-express ACAN or SOX9.
4. **Chondrocyte-like:** Chondrocyte score > 2× fibroblast score, and chondrocyte score above a minimum threshold.
5. **Fibroblast-like:** Fibroblast score > 2× chondrocyte score, and fibroblast score above a minimum threshold.
6. **Unknown:** Everything else.

### Cluster-level smoothing

To reduce noise from individual cell-level scoring:
1. Compute Leiden clusters per dataset (resolution 0.5-1.0)
2. For each cluster, compute the majority coarse_label
3. If >85% of cells in a cluster share the same label, assign that label to the entire cluster
4. Log any clusters where the majority is <70% (potential mixed clusters or doublets)

## Known Caveats

- **CD68 in IVD cells:** CD68 is expressed at low levels in stressed disc cells in some datasets. Do NOT rely on CD68 alone — require PTPRC or co-expression with other immune markers for immune classification.
- **ACTA2 in myofibroblasts:** Some degenerated disc cells may express ACTA2. Do not use ACTA2 alone for pericyte classification — require RGS5 or PDGFRB co-expression.
- **Chondrocyte-fibroblast continuum:** The 2× score ratio threshold is intentionally conservative. Cells near the boundary will be labeled Unknown, which is the correct behavior — they are fibrochondrocytes and should not be forced into either anchor.

## Automated Validation

- [ ] All cells have a `coarse_label` (including "Unknown")
- [ ] All cells have a `cell_class` ("mesenchymal", "non_mesenchymal", or "unknown")
- [ ] Proportion of "Unknown" cells is < 40% per dataset (some Unknown is expected; >40% suggests thresholds are too strict)
- [ ] Proportion of "Unknown" cells is > 5% for datasets with NP+AF cells (some continuum cells should be Unknown; <5% suggests thresholds are too loose)
- [ ] PTPRC expression is enriched in Immune cells vs. Chondrocyte_like (>10-fold)
- [ ] COL2A1 expression is enriched in Chondrocyte_like vs. Immune (>10-fold)
- [ ] COL1A1 expression is higher in Fibroblast_like than in Chondrocyte_like
- [ ] Classification report is generated per dataset
- [ ] Summary TSV is generated

## Human Checkpoint

### Review materials
- Per-dataset UMAPs colored by coarse_label
- Marker dot plots per coarse_label
- Classification summary table with Unknown proportions
- Any clusters flagged as mixed

### Questions for the reviewer
1. Are the anchor categories cleanly separated in marker expression?
2. Is the Unknown proportion reasonable (20-30% of mesenchymal cells)?
3. Are there datasets where classification looks suspicious?
4. Is the chondrocyte vs. fibroblast 2× threshold appropriate, or should it be adjusted?

### Potential plan revisions
- If Unknown proportion is too high or too low, adjust the score ratio threshold
- If CD68 is causing immune misclassification, tighten to require PTPRC co-expression
- If any anchor category has very few cells (<100 across all datasets), it may not serve as a useful anchor for scANVI
