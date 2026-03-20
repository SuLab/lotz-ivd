# Module 08: Differential Analysis

## Objective

Identify changes in cell composition and gene expression between conditions (healthy vs. degenerated, mild vs. severe, young vs. aged) to understand how IVD cell populations change with disease.

## Rationale

Two complementary questions: (1) Do the proportions of cell types/states change with disease? (2) Within a given cell type, which genes change expression with disease? Both require careful statistical methods that account for the nested structure of the data (cells within samples within studies).

## Inputs

- Integrated AnnData objects from Module 05: `data/integrated/NP.h5ad`, `data/integrated/AF.h5ad`, `data/integrated/CEP.h5ad`, `data/integrated/all_cells.h5ad`
- `metadata/sample_metadata.tsv`
- Cell type labels from post-integration de novo annotation (`obs['cell_type']` assigned in Module 07)

**Note:** GSE233666 (herniated-only) is excluded from all analysis. Neonatal samples from GSE189916 are excluded. See Module 05 spec for full study/sample assignments per object.

## Outputs

- `results/differential/composition_analysis.tsv` — cell type proportion changes
- `results/differential/de_results/{cell_type}_{comparison}.tsv` — DE results per cell type per comparison
- `results/differential/de_summary.html` — overview report
- `results/differential/volcano_plots/` — volcano plots per comparison
- `results/differential/heatmaps/` — top DE gene heatmaps

### Notebook: `notebooks/08_differential.ipynb`

Produced after all DE and composition analyses. Contains:

*Composition section:*
- Stacked bar plots: cell type proportions per sample, grouped by condition
- Box plots: proportion of each cell type by condition, with significance annotations
- Table of significant composition changes

*DE section:*
- Volcano plots for key comparisons (healthy vs. degenerated in NP chondrocytes is the lead figure)
- Heatmap of top 50 DE genes across samples with condition annotation bars
- UpSet plot: overlap of DE genes across cell types for the same comparison
- UpSet plot: overlap of DE genes across comparisons for the same cell type
- Summary table: number of significant DE genes (up/down) per cell type per comparison
- Highlight panel: known IVD-relevant genes (MMPs, collagens, inflammatory cytokines) with their DE status across all comparisons

**Manuscript mapping:** Figure 2: Cell composition changes with degeneration. Figure 3: Differential expression (volcano, heatmap). Table 2: Top DE genes. Supplementary Tables: full DE results.

## Part 1: Cell Composition Analysis

### Comparisons

Run the following comparisons (adjust based on available samples after metadata harmonization):

1. **Healthy vs. degenerated** (all severities combined, if healthy sample count is limited)
2. **Healthy vs. mild degeneration** (if sufficient samples)
3. **Healthy vs. severe degeneration** (if sufficient samples)
4. **Mild vs. severe degeneration**
5. **Young/healthy vs. aged** (using the developmental/aging axis where degeneration grade is controlled or absent)

Each comparison should be run per object (NP, AF, CEP, all-cells) where sample counts allow. CEP is likely underpowered. The all-cells object enables cross-compartment comparisons.

### Method

Use a method that properly accounts for the compositional nature of cell type proportions and the sample-level replicate structure:

**Run both methods and compare:**
- `scCODA` (Bayesian compositional analysis) — models cell type proportions as compositional data; provides credible intervals
- `propeller` (from speckle package in R) — logit-transformed proportions with empirical Bayes variance moderation

Report results from both. Cell types called significant by both methods are high-confidence; those called by only one should be noted as method-dependent. Use propeller as the primary result for the manuscript, with scCODA as a sensitivity analysis.

**Do NOT use:** Simple proportion tests (chi-squared, Fisher's) on pooled cells — these ignore sample-level variation and dramatically inflate false positives.

**Steps:**
1. Compute cell type proportions per sample (using `cell_type` labels from Module 07 de novo annotation)
2. Build a design matrix with condition as the primary variable, adjusting for study as a covariate
3. Run both scCODA and propeller for differential abundance of each cell type
4. Report effect sizes (log fold change in proportion) and significance (FDR-corrected for propeller, credible intervals for scCODA)
5. Generate a concordance table: which cell types are significant in both, only one, or neither
6. Visualize with stacked bar plots (proportions per sample) and box plots (proportion by condition)

### Handling confounds

- **Study as covariate:** Cell type proportions vary by study due to tissue processing, compartment, etc. Include study as a covariate.
- **Study-condition confounding:** Run all comparisons, but flag confounded ones. If a comparison has <2 studies contributing to either condition, add a warning in the results table and exclude it from primary manuscript figures. Only fully cross-study comparisons (≥2 studies per condition) go in the main figures.
- **Tissue vs. cells:** If both tissue-derived and cell-isolated samples are present in a comparison, either restrict to one type or include as a covariate.
- **Age-degeneration confound:** If aged samples are also degenerated, composition changes cannot be attributed to one variable. Flag this if the data doesn't allow separation.

### R interface

DESeq2 and propeller are R packages. Use standalone R scripts (in `scripts/`) called via subprocess from Python. This is more robust and debuggable than rpy2.

## Part 2: Differential Gene Expression

### Method

**Primary method:** Pseudobulk DE analysis using DESeq2.

**Why pseudobulk:** Single-cell DE methods (Wilcoxon, MAST) treat each cell as an independent observation, inflating sample sizes and producing false positives. Pseudobulk aggregates cells per sample per cell type, producing one expression profile per sample, which is the correct replicate unit.

**Assay and layer setup:**
- DE analysis must use the RNA assay (raw counts), NOT the SCT assay
- In R/Seurat: call `JoinLayers()` before running DE to merge split count layers
- In Python: use raw counts from `.layers['counts']`

**Steps:**
1. For each cell type and each comparison:
   a. Subset cells to the relevant cell type and conditions
   b. Aggregate raw counts per sample (sum across cells within each sample for each gene)
   c. Filter genes: retain genes detected in at least 10% of cells in at least one condition
   d. Run DESeq2 with design: `~ condition + study + ivd_score` (include study as batch covariate; include IVD score as a continuous covariate)
   e. Extract results: log2 fold change, adjusted p-value (BH correction), base mean expression
   f. Filter significant DE genes: |log2FC| > 0.5 AND adjusted p-value < 0.05
2. Save results to `results/differential/de_results/{cell_type}_{comparison}.tsv`

### Comparisons to run

Same comparison list as composition analysis, crossed with cell types:
- Run DE for each comparison × each cell type that has at least 3 samples per condition with at least 50 cells per sample. Skip underpowered comparisons and log them.

### Minimum requirements

- At least 3 samples per condition per cell type (for pseudobulk to be meaningful)
- At least 50 cells per sample per cell type (for reliable aggregation)
- If these minimums aren't met, log the gap and skip rather than run underpowered analysis

### For the continuum populations

If some mesenchymal clusters exist on a continuum rather than as discrete types (as indicated by Module 07 annotation), DE can also be done by:
- Binning cells by continuous score (e.g., `score_degenerative`) into quantiles and running pseudobulk within bins
- Using the continuous score as a covariate in the DE model (interaction term: condition × cell state score)
- Document the approach used

### Visualization

- Volcano plots per comparison per cell type (highlight known IVD-relevant genes)
- Heatmap of top 50 DE genes (by adjusted p-value) across samples, with condition annotation
- Upset plot or Venn diagram showing overlap of DE genes across cell types for the same comparison

## Automated Validation

- [ ] Composition analysis results exist for all planned comparisons
- [ ] DE results exist for all cell type × comparison pairs that met the minimum sample requirements
- [ ] No DE results file contains > 5000 significant genes at FDR < 0.05 (likely indicates inflated statistics — pseudobulk should prevent this, but check)
- [ ] DE results include expected positive controls where applicable (e.g., MMP13, ADAMTS5 should be upregulated in degeneration for NP chondrocytes based on prior literature)
- [ ] Volcano plots and heatmaps are generated
- [ ] Summary report is generated
- [ ] All skipped comparisons (due to insufficient power) are logged with reasons

## Human Checkpoint

### Review materials
- Composition analysis results and plots
- DE summary report
- Volcano plots for key comparisons (healthy vs. degenerated in NP chondrocytes is the highest priority)
- Top DE gene lists

### Questions for the reviewer
1. Do the composition changes make biological sense? (e.g., increased immune infiltration with degeneration is expected)
2. Are the DE results consistent with known IVD biology? (e.g., upregulation of catabolic enzymes, inflammatory cytokines in degeneration)
3. Are there unexpected findings that warrant follow-up?
4. Is there evidence that the study covariate is dominating the results (more genes associated with study than with condition)?
5. Are there any comparisons that should be added, removed, or redefined based on the results so far?
6. Should any DE results feed back into the annotation (e.g., if a subcluster of NP cells has a very different DE profile, it may be a distinct cell state)?

### Potential plan revisions
- If DE results are dominated by batch/study effects, reconsider the integration strategy or restrict analysis to within-study comparisons
- If composition changes are minimal, the focus may shift entirely to within-cell-type DE rather than cell type proportion changes
- If the continuum populations show condition-dependent shifts in continuous scores (even without DE in individual genes), this may point to a coordinated gene program that the interpretation module should explore
- If DE results reveal pain-associated genes, prioritize these for the interpretation module
