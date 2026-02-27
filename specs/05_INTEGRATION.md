# Module 05: Integration Strategy

## Objective

Combine cells across studies into a shared representation that preserves biological variation (cell types, cell states, disease effects) while removing unwanted technical variation (batch effects). The approach must be tailored to the known challenge that IVD resident cell populations exist on a continuum and are vulnerable to overcorrection.

## Rationale

This is the crux of the project's technical challenge. Previous attempts using scVI with batch correction by sample ID produced a "blob" for resident IVD cells. This module defines a tiered integration strategy that treats different cell populations differently, rather than applying one-size-fits-all integration.

## Inputs

- Processed and annotated AnnData objects from `data/processed/{accession}.h5ad`
- `metadata/sample_metadata.tsv`
- Results from per-dataset annotation (Module 04) — especially the cell type labels and continuous scores

## Outputs

- `data/integrated/tier1_nonresident.h5ad` — integrated non-resident cells (immune, endothelial, pericyte)
- `data/integrated/tier2_resident_{compartment}.h5ad` — integrated resident cells, per compartment
- `data/integrated/integration_metrics.tsv` — quantitative integration assessment
- `results/integration_report.html` — visualization of integration results
- Updated `analysis_plan.md` with the chosen strategy

### Notebook: `notebooks/05_integration.ipynb`

Produced after all integration approaches have been run. This is the key decision-support notebook. Contains:
- Side-by-side UMAP panels for each integration approach (A-F), colored by: study, cell type, condition, compartment
- Integration metrics table (kBET, LISI, ASW) for all approaches, with visual comparison (radar plot or grouped bar chart)
- Continuum preservation check: distribution of continuous cell state scores before vs. after integration, per approach
- Cluster count comparison: number of clusters at resolution 0.5 per approach
- Condition classifier accuracy: can we still distinguish healthy from degenerated after integration?
- Recommendation summary (to be confirmed at human checkpoint)

**Manuscript mapping:** Supplementary Figure S3: Integration benchmarking. Methods section on integration strategy and rationale for chosen approach.

## Strategy: Tiered Integration

### Tier 1: Non-resident cell populations

**Goal:** Integrate immune, endothelial, and pericyte cells across all studies.

**Approach:** Standard scVI or Harmony integration. These populations have strong, discrete transcriptomic identities that survive batch correction. Overcorrection is not a concern here.

**Steps:**
1. Subset all cells labeled as immune, endothelial, or pericyte from all datasets
2. Concatenate into a single AnnData
3. Re-identify HVGs on the concatenated object
4. Integrate using scVI with `batch_key='study'` (not sample — too many batches for a small number of cells)
5. Cluster and re-annotate at higher resolution (immune subtypes: M1/M2 macrophages, T cell subtypes, etc.)
6. Assess integration quality (see metrics below)

**This tier is relatively straightforward and should succeed.**

### Tier 2: Resident IVD cell populations

**Goal:** Create a cross-study representation of chondrocyte-like and fibroblast-like cells that preserves cell state variation.

**This is the hard part.** Multiple approaches should be tested and compared. The agent should run all of the following and generate comparison visualizations for human review.

#### Approach A: scVI with conservative batch correction

- Use `batch_key='study'` (coarser than sample ID — less aggressive correction)
- Reduce the number of latent dimensions (e.g., 20 instead of default 30)
- Consider using `categorical_covariate_keys=['compartment']` and `continuous_covariate_keys=['pct_counts_mt']` to model known covariates explicitly rather than lumping them into batch

#### Approach B: scANVI (semi-supervised)

- Use the per-dataset cell type annotations as seed labels
- scANVI leverages label information during integration, which should help preserve annotated cell states
- Use only high-confidence labels; leave ambiguous cells unlabeled for scANVI to assign

#### Approach C: Harmony with controlled parameters

- Test with `theta` parameter reduced (less aggressive correction) — try theta = 0.5, 1.0, 2.0
- Use study as the batch variable
- This is faster than scVI and may be sufficient if the batch effects are not severe for resident cells

#### Approach D: BBKNN (batch-balanced k-nearest neighbors)

- Does not produce a corrected expression matrix — operates on the neighbor graph
- Less likely to overcorrect because it only adjusts connectivity, not gene expression
- May preserve the continuum better than methods that transform the expression space

#### Approach E: No integration — metacell aggregation

- Instead of integrating at the single-cell level, compute metacells (using SEACells or MC2) within each dataset
- Compare metacell transcriptomic profiles across studies using correlation or classification
- This avoids the overcorrection problem entirely but sacrifices single-cell resolution
- Particularly useful if the continuum makes single-cell integration fundamentally unreliable

#### Approach F: Label transfer without forced integration

- Use a reference dataset (the largest or highest-quality study) to train an scANVI or CellTypist model
- Transfer labels to other datasets without embedding them in a shared space
- Each dataset retains its own UMAP; comparisons are done via label proportions and DE within shared labels
- This is the most conservative approach and is a fallback if integration consistently fails

### Compartment-specific integration

For Tier 2, integrate NP cells and AF cells separately (not all resident cells together). Rationale:
- NP and AF have genuinely different transcriptomic profiles
- Forcing them into one space adds unnecessary variation that competes with cell state variation
- CEP cells may be too few for independent integration; if so, include them with NP (closer biology) or analyze only within the two studies that have endplate data

### Integration parameter choices

For each approach, record all parameters used. Common parameters to vary:
- Batch variable: `study` vs. `sample_id` vs. `donor_id`
- Number of HVGs: 2000, 3000, 5000
- Number of latent dimensions (scVI/scANVI): 10, 20, 30
- Neighbor graph parameters: n_neighbors = 15 (default), 30

The agent should NOT exhaustively grid-search all combinations. Run each approach with reasonable defaults first, then adjust based on metrics.

## Integration Quality Metrics

For each integration result, compute:

1. **Batch mixing metrics:**
   - kBET (k-nearest neighbor batch effect test) — measures whether local neighborhoods contain proportional representation of batches
   - iLISI (integration local inverse Simpson's Index) — higher is better mixing
   - Batch-ASW (average silhouette width by batch) — should be near 0 (no batch separation)

2. **Biological conservation metrics:**
   - cLISI (cell-type local inverse Simpson's Index) — higher means cell types are less mixed (good)
   - Cell-type-ASW — should be positive (cell types should separate)
   - Isolated label F1 — checks whether rare cell types are preserved

3. **Continuum-specific metrics (custom):**
   - For resident IVD cells, compute the variance of continuous cell state scores (e.g., degenerative_score, notochordal_score) in the integrated space. If integration squashes this variance compared to per-dataset analysis, it's overcorrecting.
   - Compare the number of distinct clusters at resolution 0.5 before and after integration. A large reduction suggests overcorrection.
   - Check whether condition-associated variation (healthy vs. degenerated) is preserved by running a simple classifier on the integrated space.

4. **Scree metrics:**
   - scIB (single-cell integration benchmarking) aggregate score if implementable

## Automated Validation

- [ ] All integration approaches listed above have been executed (or explicitly skipped with documented reason)
- [ ] `integration_metrics.tsv` contains metrics for all approaches
- [ ] For each approach, UMAP visualizations are generated colored by: study, cell type, condition, compartment
- [ ] No approach produces a result where all resident cells form a single cluster at resolution 0.5 (the blob check — if this happens, flag it)
- [ ] No approach produces a result where study identity perfectly predicts cluster identity (failed integration)
- [ ] Integration report HTML is generated

## Human Checkpoint

This is the most critical decision point in the entire pipeline.

### Review materials
- Integration report with side-by-side UMAPs for all approaches
- Integration metrics table
- Comparison of cell state score distributions (integrated vs. per-dataset)
- Cluster counts at multiple resolutions for each approach

### Questions for the reviewer
1. Which integration approach (A-F) best preserves cell state variation while adequately removing batch effects?
2. Is any approach clearly superior, or is a combination needed (e.g., scANVI for NP, Harmony for AF)?
3. Does the "blob" problem recur with any approach? If so, is Approach E or F the appropriate fallback?
4. Should the analysis proceed with integrated data, per-dataset data, or both in parallel?
5. Are there any study-specific effects that persist after integration and need to be addressed as covariates in downstream DE analysis?
6. Does the integration reveal any new cell states not visible in per-dataset analysis?

### Potential plan revisions
- **If no integration approach preserves the continuum:** Switch to Approach F (label transfer) or Approach E (metacells) as the primary strategy. Downstream DE analysis would use pseudobulk per sample within shared labels rather than integrated single-cell data.
- **If integration works for some compartments but not others:** Use different strategies per compartment and document the rationale.
- **If condition-associated variation is lost after integration:** Add condition as a covariate in the integration model or abandon integration for condition-related analyses.
- **If integration reveals that tissue-vs-cells is a dominant confounder:** Consider restricting the analysis to tissue-only studies for the primary results, using cell-derived studies only for validation.
- **This checkpoint may result in a substantial rewrite of Modules 06-09.** That is expected and acceptable.
