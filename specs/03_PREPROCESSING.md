# Module 03: Per-Dataset Preprocessing

## Objective

Apply a consistent QC and preprocessing pipeline to each dataset independently, producing clean, normalized, and annotated AnnData objects ready for downstream analysis. Per-dataset processing avoids the integration-blob problem and establishes baseline cell populations within each study before any cross-study comparison.

## Rationale

Each dataset has different quality characteristics, sequencing depth, and cell recovery. A consistent pipeline applied independently allows us to: (1) assess each dataset's quality before deciding whether to include it in integration, (2) identify cell populations within each study using only that study's variation, and (3) compare results across studies without forcing cells into a shared embedding that may erase biological signal.

## Inputs

- Raw count matrices in `data/raw/{accession}/`
- `metadata/sample_metadata.tsv` (from Module 02)

## Outputs

Per dataset:
- `data/processed/{accession}.h5ad` — QC'd, normalized, with embeddings and preliminary clusters
- `results/qc_reports/{accession}_qc_report.html` — QC summary with plots

Aggregate:
- `results/qc_reports/qc_summary.tsv` — one row per dataset with key QC metrics
- `results/qc_reports/qc_overview.html` — cross-dataset comparison of QC metrics

## Pipeline Steps

Apply the following to each dataset independently. All parameters are defaults that may be adjusted at the human checkpoint.

### Step 1: Load and format

- Load raw count matrix into AnnData
- Ensure gene names are standardized (gene symbols, not Ensembl IDs; map if necessary using a reference like HGNC)
- Attach sample-level metadata from `sample_metadata.tsv` to each cell's `.obs`
- Record original cell barcodes in `.obs['original_barcode']`
- Record study accession in `.obs['study']`

### Step 2: Quality control

Calculate per-cell metrics:
- `n_genes_by_counts` — number of genes detected
- `total_counts` — total UMI counts
- `pct_counts_mt` — percentage of mitochondrial gene counts (genes matching `^MT-`)
- `pct_counts_ribo` — percentage of ribosomal gene counts (genes matching `^RP[SL]`)

**Default QC thresholds** (apply per-sample, not per-dataset, to account for sample-level variation):
- Remove cells with `n_genes_by_counts` < 200
- Remove cells with `n_genes_by_counts` > 6000 (potential doublets)
- Remove cells with `pct_counts_mt` > 20%
- Remove cells with `total_counts` < 500

**Doublet detection:**
- Run Scrublet (or DoubletFinder equivalent) per sample
- Flag predicted doublets in `.obs['predicted_doublet']`
- Remove cells with doublet score > 0.25 (default threshold)

**Gene filtering:**
- Remove genes detected in fewer than 3 cells

Record the number of cells and genes before and after each filter step.

### Step 3: Normalization

- Normalize to 10,000 counts per cell (`sc.pp.normalize_total(target_sum=1e4)`)
- Log-transform (`sc.pp.log1p()`)
- Store raw counts in `.layers['counts']` before normalization (needed for DE analysis and scVI)

### Step 4: Feature selection

- Identify highly variable genes (HVGs) using `sc.pp.highly_variable_genes()`
  - Method: `seurat_v3` (uses raw counts, preferred for count data)
  - `n_top_genes`: 3000 (default; may be adjusted)
  - If dataset has multiple samples, use `batch_key='sample_id'` to identify HVGs that are variable across batches

### Step 5: Dimensionality reduction

- Scale to unit variance (`sc.pp.scale(max_value=10)`) — only for PCA input, do not overwrite `.X`
- PCA: compute 50 components
- Determine effective dimensionality: use elbow method or variance explained > 90% cumulative. Record the number of PCs used.
- Compute neighbor graph using the selected number of PCs (`sc.pp.neighbors()`)
- UMAP embedding (`sc.tl.umap()`)

### Step 6: Clustering

- Leiden clustering at multiple resolutions: 0.2, 0.5, 0.8, 1.0, 1.5
- Store each in `.obs['leiden_res_{resolution}']`
- The "working" resolution for initial analysis: 0.5 (default; adjusted per dataset at human checkpoint)

### Step 7: Marker genes

- For the working resolution clustering, compute marker genes using Wilcoxon rank-sum test (`sc.tl.rank_genes_groups()`)
- Record top 50 markers per cluster
- Save to `results/annotations/{accession}_markers.tsv`

### Step 8: Preliminary cell type labels

Apply broad cell type labels based on canonical markers. This is a coarse first pass, not the final annotation.

**Expected major cell types and key markers:**

| Cell type | Key markers |
|-----------|-------------|
| Chondrocyte-like (NP) | ACAN, COL2A1, SOX9, KRT19 (NP-specific) |
| Fibroblast-like (AF) | COL1A1, COL1A2, THY1, SCX |
| Endothelial | PECAM1, VWF, CDH5, FLT1 |
| Immune - Macrophage | CD68, CD163, CSF1R, MRC1 |
| Immune - T cell | CD3D, CD3E, CD4, CD8A |
| Immune - B cell | CD79A, MS4A1, CD19 |
| Immune - Mast cell | KIT, TPSAB1, CPA3 |
| Pericyte/Smooth muscle | ACTA2, TAGLN, MYH11, RGS5 |
| Notochordal | T (Brachyury), TBXT, SHH, NOG |
| Progenitor-like | CD44, PROM1, NES |

Assign labels using a scoring approach (e.g., `sc.tl.score_genes()` for each marker set), not manual inspection. Record confidence: high (clear marker enrichment), medium (markers present but not dominant), low (ambiguous).

Store in `.obs['cell_type_preliminary']`

### Step 9: QC report generation

Generate an HTML report per dataset containing:
- Sample-level cell count (before and after QC)
- Violin plots of QC metrics (n_genes, total_counts, pct_mt) per sample
- UMAP colored by: sample, leiden clusters (at multiple resolutions), preliminary cell type, pct_mt, n_genes
- Dot plot of canonical markers across clusters
- Table of top 10 markers per cluster
- Summary statistics (total cells, cells per sample, cells per preliminary cell type)

## Automated Validation

Per dataset:
- [ ] `data/processed/{accession}.h5ad` exists and is loadable
- [ ] `.layers['counts']` contains raw integer counts
- [ ] `.obs` contains all required metadata fields from sample_metadata.tsv
- [ ] Cell count after QC is at least 50% of cell count before QC (flag if not — may indicate overly aggressive filtering or bad data)
- [ ] Cell count after QC is at least 200 per sample (minimum viability threshold)
- [ ] At least one cluster expresses canonical chondrocyte/fibroblast markers (ACAN, COL2A1, COL1A1)
- [ ] Immune markers (CD68, CD3D) are either absent (acceptable for some studies) or localized to specific clusters (not diffusely spread, which would suggest contamination or poor QC)
- [ ] No single cluster contains >80% of all cells (would suggest failed clustering)
- [ ] QC report HTML is generated
- [ ] `results/qc_reports/qc_summary.tsv` is updated with metrics for this dataset

Cross-dataset:
- [ ] All included datasets have been processed
- [ ] `results/qc_reports/qc_overview.html` is generated showing cross-dataset comparison

## Human Checkpoint

### Review materials
- QC reports for each dataset (focus on UMAPs and marker expression)
- Cross-dataset QC overview
- `qc_summary.tsv`

### Questions for the reviewer
1. Are QC thresholds appropriate for each dataset? Some datasets (e.g., from tissue digestion) may need different mitochondrial thresholds.
2. Do the preliminary cell type labels make sense? Are expected cell types present where they should be (e.g., endothelial cells in tissue-derived but not in sorted NP cell datasets)?
3. Are there datasets where quality is too low to include in downstream analysis?
4. Do any datasets show strong batch effects between samples within the same study that need correction?
5. Are there unexpected cell populations (e.g., contaminating cell types that shouldn't be in IVD)?
6. For the chondrocyte/fibroblast continuum populations — do they form one large cluster or multiple subclusters? Does the clustering resolution need adjustment?

### Potential plan revisions
- If a dataset fails QC badly, remove it from the included list and update the analysis_plan.md
- If chondrocyte/fibroblast populations are already showing interesting substructure within individual datasets, this supports the per-dataset-first strategy and may reduce the need for aggressive integration
- If some datasets have very different gene detection rates or sequencing depth, this will inform the choice of integration method (scVI handles this better than Harmony)
- If immune cells are rare or absent in some datasets, cell-cell communication analysis may need to be scoped to a subset of studies
- If QC thresholds need dataset-specific adjustment, record the rationale and update this spec for reproducibility
