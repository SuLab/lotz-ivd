# Module 06: Differential Analysis

## Objective

Identify changes in cell composition and gene expression between conditions (healthy vs. degenerated, mild vs. severe, young vs. aged) to understand how IVD cell populations change with disease.

## Rationale

Two complementary questions: (1) Do the proportions of cell types/states change with disease? (2) Within a given cell type, which genes change expression with disease? Both require careful statistical methods that account for the nested structure of the data (cells within samples within studies).

## Inputs

- Integrated or per-dataset annotated AnnData objects (depending on Module 05 outcome)
- `metadata/sample_metadata.tsv`
- Final cell type labels

## Outputs

- `results/differential/composition_analysis.tsv` — cell type proportion changes
- `results/differential/de_results/{cell_type}_{comparison}.tsv` — DE results per cell type per comparison
- `results/differential/de_summary.html` — overview report
- `results/differential/volcano_plots/` — volcano plots per comparison
- `results/differential/heatmaps/` — top DE gene heatmaps

## Part 1: Cell Composition Analysis

### Comparisons

Run the following comparisons (adjust based on available samples after metadata harmonization):

1. **Healthy vs. degenerated** (all severities combined, if healthy sample count is limited)
2. **Healthy vs. mild degeneration** (if sufficient samples)
3. **Healthy vs. severe degeneration** (if sufficient samples)
4. **Mild vs. severe degeneration**
5. **Young/healthy vs. aged** (using the developmental/aging axis where degeneration grade is controlled or absent)

Each comparison should be run per compartment where possible (NP, AF). CEP may lack sufficient samples.

### Method

Use a method that properly accounts for the compositional nature of cell type proportions and the sample-level replicate structure:

**Primary method:** `scCODA` (Bayesian compositional analysis) or `propeller` (from speckle package in R)

**Do NOT use:** Simple proportion tests (chi-squared, Fisher's) on pooled cells — these ignore sample-level variation and dramatically inflate false positives.

**Steps:**
1. Compute cell type proportions per sample (using `cell_type_final` labels)
2. Build a design matrix with condition as the primary variable, adjusting for study as a covariate
3. Test for differential abundance of each cell type
4. Report effect sizes (log fold change in proportion) and significance (FDR-corrected)
5. Visualize with stacked bar plots (proportions per sample) and box plots (proportion by condition)

### Handling confounds

- **Study as covariate:** Cell type proportions vary by study due to tissue processing, compartment, etc. Include study as a covariate.
- **Tissue vs. cells:** If both tissue-derived and cell-isolated samples are present in a comparison, either restrict to one type or include as a covariate.
- **Age-degeneration confound:** If aged samples are also degenerated, composition changes cannot be attributed to one variable. Flag this if the data doesn't allow separation.

## Part 2: Differential Gene Expression

### Method

**Primary method:** Pseudobulk DE analysis using DESeq2 or edgeR.

**Why pseudobulk:** Single-cell DE methods (Wilcoxon, MAST) treat each cell as an independent observation, inflating sample sizes and producing false positives. Pseudobulk aggregates cells per sample per cell type, producing one expression profile per sample, which is the correct replicate unit.

**Steps:**
1. For each cell type and each comparison:
   a. Subset cells to the relevant cell type and conditions
   b. Aggregate raw counts per sample (sum across cells within each sample for each gene)
   c. Filter genes: retain genes detected in at least 10% of cells in at least one condition
   d. Run DESeq2 (or edgeR) with design: `~ condition + study` (include study as covariate if samples come from multiple studies)
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

If the resident IVD cells were not discretely classified (i.e., continuous scores were used instead of hard labels), DE can still be done by:
- Binning cells by continuous score into quantiles and running pseudobulk within bins
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
