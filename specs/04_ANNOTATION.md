# Module 04: Coarse Cell Classification (Mesenchymal vs. Non-Mesenchymal)

## Objective

Classify cells into two broad categories — mesenchymal (chondrocytes, fibroblasts) and non-mesenchymal (immune, endothelial) — using unambiguous marker genes. This coarse classification enables tiered integration in Module 05. Fine-grained cell type annotation happens *after* integration (in Module 05), not here.

## Rationale

The goal of this project is to build a de novo IVD cell atlas from integrated data, without relying on original manuscript annotations. Fine annotation before integration risks circular logic: labels shape integration, which shapes downstream analysis. Instead, we perform only the minimal classification needed to separate mesenchymal from non-mesenchymal cells (which have fundamentally different transcriptomic profiles and require separate integration to avoid one group dominating the latent space).

This coarse split is reliable because the markers distinguishing mesenchymal from non-mesenchymal cells are unambiguous — PTPRC/CD45 marks immune cells, PECAM1 marks endothelial cells, and these are well-separated from COL2A1/ACAN-expressing mesenchymal cells.

Note: we use "mesenchymal" (chondrocytes, fibroblasts) vs. "non-mesenchymal" (immune, endothelial) rather than "resident" vs. "non-resident," because endothelial cells can be resident in vascularized IVD regions (outer AF, endplates).

## Inputs

- Processed AnnData objects from `data/processed/{accession}.h5ad`
- **Excluded from all analysis:** GSE233666 (Guo 2023) — herniated discs only, not the focus of this project
- **Excluded samples:** Neonatal samples from GSE189916 (Jiang 2022) — neonatal disc biology is fundamentally different from adult; only adult samples are retained

## Outputs

Per dataset:
- Updated `.obs` in the h5ad file with `cell_class` column ("mesenchymal" or "non_mesenchymal")
- `results/annotations/{accession}_classification_report.html`

Aggregate:
- `results/annotations/classification_summary.tsv` — cell class proportions per dataset

### Notebook: `notebooks/04_classification.ipynb`

Contains:
- Per-dataset UMAPs colored by cell_class
- Dot plot of classification markers (PTPRC, CD3D, CD68, PECAM1, COL2A1, ACAN, etc.) by cell_class
- Bar plot: mesenchymal vs. non-mesenchymal proportions per dataset
- Validation: marker expression distributions confirming clean separation

## Classification Method

### Marker genes for classification

**Non-mesenchymal markers** (any of these expressed → non-mesenchymal candidate):

*Immune:*
- Pan-immune: PTPRC (CD45)
- T cells: CD3D, CD3E
- Macrophage/monocyte: CD68, CD14, CSF1R
- B cells: CD79A, MS4A1
- Mast cells: KIT, TPSAB1
- NK cells: NKG7, GNLY

*Endothelial:*
- PECAM1 (CD31), VWF, CDH5

*Pericyte/SMC:*
- ACTA2, RGS5, PDGFRB

**Mesenchymal markers** (expected in IVD disc cells):
- COL2A1, COL1A1, ACAN, SOX9, DCN, LUM, VCAN

### Classification logic

For each cell, compute a simple score:
1. Score each cell for non-mesenchymal markers (fraction of non-mesenchymal markers detected above background)
2. Score each cell for mesenchymal markers (fraction of mesenchymal markers detected above background)
3. Assign cell_class based on which score is higher
4. Cells where both scores are near zero or tied → "ambiguous" (these should be rare; flag if >5%)

This can be implemented via `sc.tl.score_genes()` for each marker set, or simply by checking expression of key markers (PTPRC > 0 is a strong immune indicator).

### Cluster-level assignment

To reduce noise from individual cell-level scoring:
1. Compute Leiden clusters per dataset (resolution 0.5-1.0)
2. For each cluster, compute the majority cell_class
3. Assign the majority class to all cells in the cluster
4. Log any clusters where the majority is <70% (potential mixed clusters or doublets)

### Known caveats

- **CD68 in IVD cells:** CD68 is expressed at low levels in stressed disc cells in some datasets. Do NOT rely on CD68 alone — require co-expression with other immune markers (PTPRC, CD14) for immune classification.
- **ACTA2 in myofibroblasts:** Some degenerated disc cells may express ACTA2. Classify as mesenchymal unless they also express pericyte markers (RGS5, PDGFRB).

## Automated Validation

- [ ] All cells have a `cell_class` label ("mesenchymal", "non_mesenchymal", or "ambiguous")
- [ ] Proportion of "ambiguous" cells is < 5% per dataset
- [ ] PTPRC expression is enriched in non_mesenchymal cells vs. mesenchymal (>10-fold)
- [ ] COL2A1 or COL1A1 expression is enriched in mesenchymal cells vs. non_mesenchymal
- [ ] Classification report is generated per dataset
- [ ] Summary TSV is generated

## Human Checkpoint

### Review materials
- Per-dataset UMAPs colored by cell_class
- Marker dot plots confirming clean separation
- Classification summary table

### Questions for the reviewer
1. Is the mesenchymal vs. non-mesenchymal split clean? Any datasets with suspicious proportions?
2. Are there clusters flagged as mixed that need manual review?
3. Are the ambiguous cells a concern, or can they be excluded?

### Potential plan revisions
- If any dataset has >20% ambiguous cells, investigate and potentially re-threshold
- If CD68 is causing misclassification, tighten the immune classification to require PTPRC co-expression
