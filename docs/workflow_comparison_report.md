# Workflow Comparison Report: Hannah's Seurat/R Pipeline vs. Current Scanpy/Python Pipeline

**Prepared for discussion between the Good and Lotz labs**
**Date: 2026-03-13**

---

## Purpose

This report provides a detailed, decision-by-decision comparison of the two single-cell RNA-seq analysis workflows under consideration for the IVD atlas project:

- **Workflow A** — Hannah's validated R/Seurat workflow (`single_nuclei_r/`), developed for single-nuclei synovium analysis
- **Workflow B** — The current IVD pipeline (`scripts/03_preprocessing.py`, `scripts/05_integration.py`), built in Python/Scanpy with scVI integration

The goal is to evaluate each methodological choice on its merits, identify where one approach is clearly superior, where it's a toss-up, and where context matters. We explicitly aim to avoid blindly adopting either workflow.

---

## Decision 1: Normalization — SCTransform vs. normalize_total + log1p

### Workflow A (Hannah): SCTransform

```r
SCTransform(object, vars.to.regress = c("percent.mt"))
```

Uses regularized negative binomial regression (Hafemeister & Satija 2019; v2 in Choudhary & Satija 2022). Models the relationship between sequencing depth and gene expression per gene, then produces Pearson residuals as the "normalized" expression values. Variance stabilization is inherent — no separate HVG selection step required (though one is still applied downstream for integration features).

**Advantages:**
- Properly handles the mean-variance relationship in count data (heteroskedasticity)
- Does not distort the relative expression of highly vs. lowly expressed genes (a known problem with log-normalization)
- Variance-stabilized residuals give better PCA/downstream results, especially for lowly-expressed genes that are biologically important (e.g., transcription factors, signaling molecules)
- The v2 update (2022) includes improved regularization and is faster
- Can regress out covariates (percent.mt) during normalization rather than as a separate step
- Well-validated in thousands of published studies

**Disadvantages:**
- R/Seurat-specific (no production-quality Python implementation). Using it commits the pipeline to R for preprocessing
- Computationally heavier than simple normalization (minutes vs. seconds per sample)
- The Pearson residuals can behave unexpectedly for genes with very high expression (a known issue documented by Lause et al. 2021, "Analytic Pearson residuals")
- Regressing out percent.mt during normalization (rather than filtering) can mask quality issues

### Workflow B (Current): normalize_total + log1p

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

Library-size normalization followed by log(x+1) transform. The standard scanpy approach.

**Advantages:**
- Simple, fast, well-understood
- Language-agnostic (works identically in any framework)
- Does not distort the data in unexpected ways — what you see is what you get
- For scVI integration (which models raw counts directly), the normalization of .X is only used for HVG selection and visualization — the integration model uses `.layers['counts']` anyway

**Disadvantages:**
- log-normalization introduces known biases: it inflates variance of lowly-detected genes and compresses variance of highly-expressed genes
- Library-size normalization assumes all cells have the same total RNA, which is biologically incorrect
- HVG selection on log-normalized data can be biased toward high-expression genes
- Does not account for the mean-variance relationship inherent in count data

### Verdict

**SCTransform is the better normalization method in general.** The statistical arguments are well-established. However, the practical impact depends on what happens downstream:

- **If integration uses raw counts (scVI):** The normalization of `.X` matters primarily for HVG selection, PCA visualization, and marker gene testing. scVI re-learns the generative model from raw counts regardless. The impact of normalization choice is attenuated.
- **If integration uses CCA on normalized data (Seurat):** The normalization directly feeds into the integration. Here, SCTransform makes a much larger difference.
- **If the downstream analysis stays in R/Seurat:** SCTransform is clearly preferred.

**Risk of adopting SCTransform:** Commits the pipeline to R for all preprocessing and integration steps. This is not inherently bad, but requires R expertise for maintenance and debugging.

---

## Decision 2: Integration Method — Seurat CCA vs. scVI

### Workflow A (Hannah): Seurat CCA with SCT

```r
Synovium.features.RNA <- SelectIntegrationFeatures(object.list = data.list.RNA, nfeatures = 3000)
data.list.RNA <- PrepSCTIntegration(object.list = data.list.RNA, anchor.features = Synovium.features.RNA)
Anchors <- FindIntegrationAnchors(object.list = data.list.RNA, dims = 1:50,
                                   normalization.method = "SCT", anchor.features = Synovium.features.RNA)
Integrated <- IntegrateData(anchorset = Anchors, dims = 1:50, normalization.method = "SCT")
```

Canonical Correlation Analysis identifies shared sources of variation across datasets, then uses mutual nearest neighbors (MNNs) to find "anchors" — cell pairs that are likely the same type across datasets. The integrated data is then corrected using these anchors.

**Advantages:**
- Extensively benchmarked in large-scale comparisons (Luecken et al. 2022 "Benchmarking atlas-level data integration"; Tran et al. 2020)
- Transparent: you can inspect anchors, see which cells are being connected across datasets
- Does not require GPU or deep learning
- Mature, stable implementation (Seurat v4/v5)
- Works directly on the corrected expression matrix — downstream clustering, DE, visualization all use the same representation
- CCA is specifically designed to find shared variation, which is exactly what integration needs
- Well-validated for combining datasets from different labs

**Disadvantages:**
- Scales poorly with very large numbers of cells (>500K). The IVD atlas has ~423K cells — this is at the edge
- CCA assumes that the primary axes of shared variation correspond to biology, not batch effects. If datasets share technical artifacts (e.g., stress signatures), CCA can lock onto these
- The anchor-based approach can overcorrect when cell types are present in very different proportions across datasets (the "rare cell type problem")
- Does not explicitly model platform-specific effects (10x vs. BD Rhapsody vs. Singleron) — these are treated the same as any batch
- Creates a "corrected" expression matrix that is a mathematical construction, not real expression values. Using this for DE analysis is inappropriate (Seurat docs themselves warn against this)

### Workflow B (Current): scVI

```python
scvi.model.SCVI.setup_anndata(adata, layer='counts', batch_key='study')
model = scvi.model.SCVI(adata, n_latent=20, dispersion='gene-batch', gene_likelihood='nb')
model.train(max_epochs=200, early_stopping=True)
embedding = model.get_latent_representation()
```

A variational autoencoder that learns a latent representation of each cell while explicitly modeling batch effects and the count distribution.

**Advantages:**
- Models raw counts directly with a proper count distribution (negative binomial) — does not require normalized data as input
- Explicitly models batch effects as a covariate, separately from biological variation
- `dispersion='gene-batch'` allows different technical noise per gene per batch — ideal for multi-platform data (10x, BD Rhapsody, Singleron)
- Scales well to large datasets (GPU acceleration available)
- Produces a probabilistic generative model that can be used for imputation, differential expression, and other downstream tasks
- scANVI extension (in the spec, not yet implemented) can use semi-supervised labels to guide integration
- Recent benchmarks (Luecken et al. 2022) rank scVI among the top methods for atlas-level integration

**Disadvantages:**
- Black-box latent space — harder to interpret what the 20 dimensions represent
- Sensitive to hyperparameters (n_latent, n_epochs, learning rate, batch_size). Poor choices can lead to underfitting or overfitting
- Requires more computational resources (GPU recommended)
- Less mature ecosystem for post-integration analysis compared to Seurat
- The latent representation is not an expression matrix — cannot be used directly for gene-level analysis without going back to the decoder
- Reproducibility can be an issue (different random seeds → somewhat different embeddings)
- Training stability: model can sometimes fail to converge or produce degenerate solutions

### Verdict

**This is the most consequential choice and genuinely a toss-up — context matters.**

**Arguments for CCA (Workflow A):**
- More transparent and interpretable
- Hannah's team has direct experience validating CCA results
- For 16 samples of the same tissue/platform (Hannah's synovium study), CCA is excellent
- Produces a corrected expression matrix usable for many downstream analyses

**Arguments for scVI (Workflow B):**
- Better theoretical fit for the IVD atlas: 12 studies, 3 platforms (10x, BD Rhapsody, Singleron), ~423K cells
- `dispersion='gene-batch'` explicitly models platform-specific technical effects — CCA does not
- Scales better to atlas-size data
- The probabilistic framework allows proper uncertainty quantification

**The critical question is whether the IVD integration challenges are due to the normalization feeding into the integration, or the integration method itself.** A fair test would be: run CCA on SCT-normalized data AND run scVI on raw counts for the same dataset, compare results.

**Recommendation:** Run both on a single IVD dataset (e.g., GSE230809, which has the most samples and NP+AF compartments) and compare integration quality metrics before committing.

---

## Decision 3: QC Thresholds

### Mitochondrial Cutoff: 5% (Workflow A) vs. 20% (Workflow B)

This is the most striking difference and has clear implications.

**Workflow A: percent.mt < 5%**
- Appropriate for **single-nuclei** data where intact nuclei should contain very little mitochondrial RNA (mitochondria are cytoplasmic)
- Strict cutoff aggressively removes damaged/lysed cells
- Risk: may remove metabolically active cell types that legitimately have higher mitochondrial content (e.g., chondrocytes in the IVD are metabolically active under hypoxia)

**Workflow B: pct_counts_mt < 20%**
- Appropriate for **single-cell** (whole-cell dissociation) data where mitochondrial content varies naturally
- Permissive cutoff retains more cells but may include damaged cells
- Risk: retains cells with compromised membranes, adding noise to downstream analysis

**Verdict:** The 20% threshold is too permissive for the IVD atlas, especially for datasets that used nuclear isolation. However, 5% may be too aggressive — some IVD cell types (NP cells under hypoxic stress, activated immune cells) genuinely express more mitochondrial genes. **A threshold of 10-15% for single-cell and 5% for single-nuclei datasets, applied per-protocol, is more defensible.** The current pipeline could benefit from protocol-aware QC thresholds.

### Count/Gene Filters

| Filter | Workflow A | Workflow B |
|--------|-----------|-----------|
| Min counts (UMI) | 1,000 | 500 |
| Max counts (UMI) | 25,000 | — |
| Min genes | — | 200 |
| Max genes | — | 6,000 |

**Workflow A** filters on total UMI counts (1K–25K). The upper bound catches doublets via library size.
**Workflow B** filters on gene counts (200–6K) and min UMI (500). The upper gene bound catches doublets; min gene bound catches empty droplets.

**Verdict:** These are complementary, not contradictory. Filtering on both total counts and gene counts is defensible. Workflow A's higher minimum (1,000 UMI) is more aggressive and may remove low-quality cells that Workflow B would retain. For IVD specifically, the quality of tissue dissociation varies — some datasets may benefit from higher minimums. **A combined approach (min_counts=1000, min_genes=200, max_genes=6000, max_counts=25000) would be more thorough.**

### Doublet Detection

**Workflow A:** Not shown in the provided code. May occur upstream or may rely solely on the max-count filter.
**Workflow B:** Scrublet per sample, with score threshold 0.25.

**Verdict:** Explicit doublet detection (Workflow B) is better practice than relying on count-based filters alone. Doublets with moderate counts will be missed by a max-count threshold. **Scrublet or DoubletFinder should be included regardless of which pipeline is adopted.**

---

## Decision 4: PCA Dimensionality

### Workflow A: All 50 PCs used for everything

```r
RunPCA(object, npcs = 50)
RunUMAP(object, dims = 1:50)
FindNeighbors(object, dims = 1:50)
FindIntegrationAnchors(object.list, dims = 1:50)
```

Uses all 50 computed PCs for neighbors, UMAP, and integration.

### Workflow B: Variance-based selection

```python
cumvar = np.cumsum(var_ratio)
n_pcs = int(np.searchsorted(cumvar, 0.90) + 1)  # 90% cumulative variance
n_pcs = max(n_pcs, 10)
```

Selects PCs based on cumulative variance explained (90% threshold, minimum 10).

**Verdict:** This is a minor difference. Using all 50 PCs includes noise dimensions, but CCA and scVI both handle this reasonably well. The variance-based approach is more principled but can be overly conservative (10 PCs may miss real biological variation in the later components). For atlas-scale integration, **30-50 PCs is standard and reasonable** — the exact number matters less than other choices. However, for per-sample preprocessing where you're just computing neighbors for clustering, variance-based selection is fine.

---

## Decision 5: Clustering Resolution Selection

### Workflow A: Systematic silhouette + modularity scoring

Hannah's second script (`Silhouette_Modularity_Scoring.R`) implements a rigorous approach:
1. Cluster at 20 resolutions (0.1–2.0)
2. Compute silhouette score (subsampled, on PCA space) — measures cluster separation
3. Compute modularity (full SNN graph via igraph) — measures graph community structure
4. Normalize and combine both metrics
5. Select the resolution maximizing the combined score

**Advantages:**
- Data-driven, reproducible, defensible
- Balances two complementary metrics (geometry vs. graph topology)
- Subsampling for silhouette handles memory constraints
- Generates interpretable diagnostic plots

**Disadvantages:**
- Computationally expensive (20 rounds of clustering + distance matrices)
- The equal weighting of silhouette and modularity is arbitrary — why not 60/40 or use a different combination?
- Neither metric directly measures biological relevance (a perfect technical clustering may not correspond to real cell types)

### Workflow B: Multiple resolutions stored, working resolution 0.5

Clusters at 5 resolutions (0.2, 0.5, 0.8, 1.0, 1.5), stores all, uses 0.5 as default. The spec mentions silhouette + modularity scoring for post-integration clustering but it's not yet implemented.

**Verdict:** Workflow A's approach is clearly more rigorous for resolution selection. **This should be adopted regardless of other choices.** The only caveat is that optimal resolution should be evaluated in the context of known biology (e.g., do the resulting clusters correspond to expected cell types?), not purely on mathematical metrics.

---

## Decision 6: Marker Gene Detection

### Workflow A: FindAllMarkers with PrepSCTFindMarkers

```r
DefaultAssay(object) <- "SCT"
object <- PrepSCTFindMarkers(object)
markers <- FindAllMarkers(object, logfc.threshold = 0.25, min.pct = 0.1,
                          only.pos = TRUE, recorrect_umi = FALSE)
```

Uses Seurat's Wilcoxon test on SCT-corrected counts. `PrepSCTFindMarkers` recorrects the SCT residuals for proper comparison.

### Workflow B (planned): Pseudobulk DESeq2

The spec calls for DESeq2 pseudobulk for differential expression (Module 08), which aggregates cells per sample before testing — treating samples (not cells) as replicates.

**Verdict:** These address different questions.

- **FindAllMarkers (Workflow A):** Answers "which genes distinguish cluster X from all others?" — a marker discovery step. Treats each cell as an independent observation, which inflates significance but is standard for marker identification.
- **Pseudobulk DESeq2 (Workflow B):** Answers "which genes differ between conditions?" — a proper statistical test for differential expression. Correctly accounts for biological replication.

**Both are needed.** Marker discovery (FindAllMarkers or scanpy's rank_genes_groups) identifies cluster-defining genes. Pseudobulk DE tests for condition effects. These are complementary steps, not alternatives. Workflow B's plan for pseudobulk DE is statistically superior for condition comparisons.

---

## Decision 7: Per-sample Processing vs. Per-dataset Processing

### Workflow A: Independent per-sample processing, then integrate

Each sample is loaded, QC'd, and SCTransform'd independently, then saved as a separate Seurat object. Integration happens on the list of independently-processed objects.

### Workflow B: Per-dataset processing (samples merged within study first), then integrate across studies

Samples within a study are merged first, then QC'd together. HVGs are identified per-dataset with batch_key='sample_id'. Integration happens across studies.

**Verdict:** Workflow A's per-sample processing is more conservative and allows per-sample SCTransform models (important because SCTransform learns per-sample regression parameters). Workflow B's per-dataset approach is more practical for 78 samples across 12 studies and allows within-study batch-aware HVG selection.

**For the IVD atlas specifically:** Per-sample SCTransform → per-study merge → across-study integration would be ideal. This is how Seurat v5 recommends handling multi-sample experiments.

---

## Decision 8: Tiered Integration (IVD-specific)

### Workflow A: Not applicable
Hannah's synovium workflow integrates all samples together — appropriate for a single tissue type.

### Workflow B: Tiered (mesenchymal vs. non-mesenchymal)
Separates cells into two tiers before integration to avoid wasting latent capacity on the major mesenchymal/non-mesenchymal axis.

**Verdict:** The tiered approach is specific to the IVD atlas and addresses a real concern — chondrocyte/fibroblast continuum biology can be masked by integration with transcriptomically very different cell types (immune, endothelial). **This should be retained regardless of which integration method is used.** It would work equally well with CCA (run FindIntegrationAnchors separately on mesenchymal and non-mesenchymal subsets).

---

## Decision 9: Framework Choice — R/Seurat vs. Python/Scanpy

This is not a technical decision per se but has practical implications.

### R/Seurat
- More mature for standard single-cell workflows (QC → normalize → integrate → cluster → annotate)
- SCTransform, Seurat v5 integration, FindMarkers are all native
- Better for interactive exploration (RStudio)
- DESeq2, propeller (compositional analysis), CellChat all native in R
- Seurat v5 has improved scaling for large datasets (BPCells, sketch-based analysis)

### Python/Scanpy
- Better for atlas-scale computation (scVI, GPU support)
- More flexible for custom analysis pipelines
- AnnData format is more portable (.h5ad)
- Better integration with deep learning methods
- Stronger ecosystem for trajectory analysis (scVelo, CellRank)

**Verdict:** The choice should follow the team's expertise and the downstream analysis needs. For this project, **many of the downstream analyses (DESeq2, propeller, CellChat) are R-native**, suggesting that an R-based preprocessing and integration step would reduce friction. The key question is whether the team has the R expertise to maintain the pipeline.

---

## Summary Table

| Decision | Recommendation | Confidence |
|----------|---------------|------------|
| Normalization | SCTransform is better; adopt from Workflow A | **High** |
| Integration method | Test both CCA and scVI head-to-head on one dataset | **Medium** — context-dependent |
| MT cutoff | Protocol-aware (5% nuclear, 10-15% cell), not a single threshold | **High** |
| Count/gene filters | Combine both approaches (more thorough) | **High** |
| Doublet detection | Keep Scrublet/DoubletFinder (Workflow B) | **High** |
| PCA dims | Minor difference; 30-50 is fine | **Low** impact |
| Resolution selection | Adopt silhouette + modularity (Workflow A) | **High** |
| Marker detection | Both needed; keep pseudobulk DE for conditions | **High** |
| Tiered integration | Keep (Workflow B, IVD-specific) | **High** |
| Framework | Follow team expertise; R has less friction for downstream | **Medium** |

---

## Proposed Validation Experiment

Before committing to either approach for the full atlas:

1. **Select one IVD dataset** from our own lab (e.g., from GSE230809 or whichever dataset Hannah's team can validate against known results)
2. **Run both workflows** on this dataset:
   - Workflow A: SCTransform → CCA (adapted for single study, multi-sample)
   - Workflow B: normalize_total + log1p → scVI
3. **Compare:**
   - Cluster number and composition
   - Known cell type recovery (do both find the expected populations?)
   - Marker gene overlap
   - UMAP structure
   - Integration metrics (if multi-sample)
4. **Use the results** to make an informed choice, then scale to the full atlas

This validation step is fast (one dataset, a few hours of compute) and eliminates guesswork.

---

## Questions for Discussion

1. Does Hannah's team have a dataset where the "ground truth" cell types are well-established, so we can validate against known biology?
2. Is the team comfortable maintaining an R-based pipeline, or should we aim for a hybrid (R for preprocessing/integration, Python for atlas-scale downstream)?
3. For the tiered integration: should we tier within CCA as well (separate FindIntegrationAnchors for mesenchymal vs. non-mesenchymal)?
4. The IVD atlas includes 3 non-10x platforms (BD Rhapsody, Singleron). Does Hannah's team have experience integrating across platforms with CCA, or has the synovium work been 10x-only?
5. Should we consider Seurat v5's newer integration methods (e.g., IntegrateLayers with CCA or Harmony) which scale better than v4's pairwise anchor approach?
