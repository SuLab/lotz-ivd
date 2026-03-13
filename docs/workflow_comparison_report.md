# Workflow Comparison Report: Hannah's Seurat/R Pipeline vs. Current Scanpy/Python Pipeline

**Prepared for discussion between the Good and Lotz labs**
**Date: 2026-03-13**

---

## Purpose

This report provides a detailed, decision-by-decision comparison of the two single-cell RNA-seq analysis workflows under consideration for the IVD atlas project:

- **Workflow A** — Hannah's validated R/Seurat workflow (`single_nuclei_r/`), developed for single-nuclei synovium analysis
- **Workflow B** — The current IVD pipeline (v4, completed 2026-03-11), built in Python/Scanpy with scANVI semi-supervised integration across a 12-module architecture

The IVD pipeline has gone through 4 iterations (v1–v4). The current v4 uses scANVI (semi-supervised) integration with 5 coarse anchor labels, adaptive resolution-optimized clustering, and two-stage annotation — resolving 19 cell types across 4 integrated objects (NP, AF, CEP, all-cells). Key findings (e.g., PTGS2 in AF inner cells, CXCL2 in fibrocartilaginous NP) have been stable across all four versions.

The goal of this report is to evaluate each methodological choice on its merits, identify where one approach is clearly superior, where it's a toss-up, and where context matters. We explicitly aim to avoid blindly adopting either workflow.

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
- For scANVI integration (which models raw counts directly), the normalization of .X is only used for HVG selection and visualization — the integration model uses `.layers['counts']` anyway

**Disadvantages:**
- log-normalization introduces known biases: it inflates variance of lowly-detected genes and compresses variance of highly-expressed genes
- Library-size normalization assumes all cells have the same total RNA, which is biologically incorrect
- HVG selection on log-normalized data can be biased toward high-expression genes
- Does not account for the mean-variance relationship inherent in count data

### Verdict

**SCTransform is the better normalization method in general.** The statistical arguments are well-established. However, the practical impact depends on what happens downstream:

- **If integration uses raw counts (scANVI/scVI):** The normalization of `.X` matters primarily for HVG selection, PCA visualization, and marker gene testing. scANVI re-learns the generative model from raw counts regardless. The impact of normalization choice is attenuated.
- **If integration uses CCA on normalized data (Seurat):** The normalization directly feeds into the integration. Here, SCTransform makes a much larger difference.
- **If the downstream analysis stays in R/Seurat:** SCTransform is clearly preferred.

**Risk of adopting SCTransform:** Commits the pipeline to R for all preprocessing and integration steps. This is not inherently bad, but requires R expertise for maintenance and debugging.

---

## Decision 2: Integration Method — Seurat CCA vs. scANVI

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

### Workflow B (Current v4): scANVI (semi-supervised)

```python
# Step 1: Train scVI base model (unsupervised)
scvi.model.SCVI.setup_anndata(adata, layer='counts', batch_key='study')
base_model = scvi.model.SCVI(adata, n_latent=20, dispersion='gene-batch', gene_likelihood='nb')
base_model.train(max_epochs=200, early_stopping=True)

# Step 2: Initialize scANVI from scVI, using coarse anchor labels
scanvi_model = scvi.model.SCANVI.from_scvi_model(base_model, labels_key='coarse_label',
                                                   unlabeled_category='Unknown')
scanvi_model.train(max_epochs=50, early_stopping=True)
embedding = scanvi_model.get_latent_representation()
```

A two-stage approach: first train an unsupervised variational autoencoder (scVI) to learn a base latent representation of each cell while explicitly modeling batch effects and count distributions, then initialize scANVI which refines the latent space using 5 coarse biological anchor labels (Chondrocyte_like, Fibroblast_like, Immune, Endothelial, Pericyte_SMC) while leaving "Unknown" cells free to be positioned by transcriptomic similarity.

**Important context:** The pipeline evolved through 4 versions — v1 benchmarked 4 integration methods (scVI, scANVI, Harmony, BBKNN), v2 simplified to scVI-only, v3 fixed annotation issues, and v4 implemented the full scANVI semi-supervised approach with coarse anchors from a dedicated classification module (Module 04). The current v4 has been run end-to-end, producing 19 resolved cell types and 23 powered DE comparisons.

**Advantages:**
- Models raw counts directly with a proper count distribution (negative binomial) — does not require normalized data as input
- Explicitly models batch effects as a covariate, separately from biological variation
- `dispersion='gene-batch'` allows different technical noise per gene per batch — ideal for multi-platform data (10x, BD Rhapsody, Singleron)
- **Semi-supervised**: scANVI uses coarse anchor labels to guide integration, conceptually similar to CCA's anchor approach but in a probabilistic framework. This constrains the latent space with biological priors without imposing fine-grained cell type identity
- Scales well to large datasets (GPU acceleration available)
- Produces a probabilistic generative model that can be used for imputation, differential expression, and other downstream tasks
- Recent benchmarks (Luecken et al. 2022) rank scVI/scANVI among the top methods for atlas-level integration
- **Already validated on IVD data**: v4 results show 19 cell types resolved, key findings stable across 4 pipeline versions, and 246 TF-activity associations recovered (vs. 5 in v3 with unsupervised scVI)

**Disadvantages:**
- Latent space is harder to interpret than CCA dimensions (though the semi-supervised anchors partially address this)
- Sensitive to hyperparameters (n_latent, n_epochs, learning rate, batch_size). Poor choices can lead to underfitting or overfitting
- Requires more computational resources (GPU recommended)
- Less mature ecosystem for post-integration analysis compared to Seurat
- The latent representation is not an expression matrix — cannot be used directly for gene-level analysis without going back to the decoder
- Reproducibility can be an issue (different random seeds → somewhat different embeddings)
- Quality of integration depends on the quality of the coarse anchor labels from the upstream classification step

### Verdict

**This is the most consequential choice, and both approaches have genuine strengths in different contexts.**

**Arguments for CCA (Workflow A):**
- More transparent and interpretable
- Hannah's team has direct experience validating CCA results
- For 16 samples of the same tissue/platform (Hannah's synovium study), CCA is excellent
- Produces a corrected expression matrix usable for many downstream analyses
- No dependency on a upstream classification step for anchor labels

**Arguments for scANVI (Workflow B):**
- Better theoretical fit for the IVD atlas: 12 studies, 3 platforms (10x, BD Rhapsody, Singleron), ~423K cells
- `dispersion='gene-batch'` explicitly models platform-specific technical effects — CCA does not
- Scales better to atlas-size data
- The probabilistic framework allows proper uncertainty quantification
- Semi-supervised approach uses biological priors (like CCA anchors) but within a generative model
- Already run end-to-end on the IVD data with stable, biologically coherent results

**The critical question is whether the IVD integration challenges are due to the normalization feeding into the integration, or the integration method itself.** A fair test would be: run CCA on SCT-normalized data AND run scANVI on raw counts for the same dataset, compare results. Note that scANVI's semi-supervised approach makes this a fairer comparison than the earlier scVI-only versions — both CCA and scANVI now use biological priors to guide integration.

**Recommendation:** Run both on a single IVD dataset (e.g., GSE230809, which has the most samples and NP+AF compartments) and compare integration quality metrics before committing. The IVD pipeline's v4 scANVI results provide a baseline to compare against.

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

### Workflow B (v4): Adaptive resolution optimization with silhouette scoring

The v4 pipeline (`scripts/06_clustering.py`) implements adaptive resolution optimization that scales the search based on dataset size:
- 300K+ cells: Test 3 resolutions (0.4, 0.8, 1.0)
- 200K–300K cells: Test 6 resolutions (0.2, 0.4, 0.6, 0.8, 1.0, 1.5)
- 50K–200K cells: Test 10 resolutions (0.2–2.0 in 0.2 steps)
- <50K cells: Test 20 resolutions (0.1–2.0 in 0.1 steps)

Uses silhouette scoring on the scANVI embedding to select the optimal resolution. Modularity scoring is skipped for objects >100K cells (too expensive).

**Advantages:**
- Adapts to dataset size (avoids wasting compute on huge objects)
- Silhouette scoring is data-driven

**Disadvantages:**
- Only uses silhouette, not the combined silhouette + modularity metric that Hannah's approach uses
- Fewer resolution candidates for large objects (3 options for 300K+ cells is quite coarse)
- No subsampling strategy documented — may be memory-limited for large objects

### Verdict

Both workflows now implement data-driven resolution selection, which is good. **Workflow A's combined silhouette + modularity approach is more rigorous** — the two metrics capture complementary aspects of cluster quality (geometric separation vs. graph community structure). Workflow B's adaptive scaling is pragmatic but loses information by dropping modularity for large objects and testing fewer resolutions. **The ideal would be to combine both: Hannah's dual-metric scoring with Workflow B's adaptive scaling for computational feasibility.** The only caveat is that optimal resolution should be evaluated in the context of known biology (e.g., do the resulting clusters correspond to expected cell types?), not purely on mathematical metrics.

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

### Workflow B (v4): Two-stage marker detection + Pseudobulk DESeq2

The v4 pipeline uses two complementary approaches:
1. **Marker discovery** (Module 07, `scripts/07_annotation.py`): Two-stage annotation — first assign coarse types via canonical markers, then discover fine-grained subtypes via DE within each coarse group. Uses scanpy's `rank_genes_groups` (Wilcoxon) for marker identification.
2. **Condition DE** (Module 08, `scripts/08_differential.py`): DESeq2 pseudobulk for differential expression between conditions (healthy vs. degenerated, mild vs. severe), aggregating cells per sample to properly account for biological replication. This is fully implemented and has produced 23 powered comparisons with 772 unique DE genes in v4.

**Verdict:** These address different questions.

- **FindAllMarkers (Workflow A):** Answers "which genes distinguish cluster X from all others?" — a marker discovery step. Treats each cell as an independent observation, which inflates significance but is standard for marker identification.
- **Pseudobulk DESeq2 (Workflow B):** Answers "which genes differ between conditions?" — a proper statistical test for differential expression. Correctly accounts for biological replication.

**Both are needed.** Marker discovery (FindAllMarkers or scanpy's rank_genes_groups) identifies cluster-defining genes. Pseudobulk DE tests for condition effects. These are complementary steps, not alternatives. Workflow B's pseudobulk DE is statistically superior for condition comparisons and has already produced results on the IVD data.

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
| Integration method | Test CCA vs. scANVI head-to-head on one dataset; scANVI has strong IVD results already | **Medium** — context-dependent |
| MT cutoff | Protocol-aware (5% nuclear, 10-15% cell), not a single threshold | **High** |
| Count/gene filters | Combine both approaches (more thorough) | **High** |
| Doublet detection | Keep Scrublet/DoubletFinder (Workflow B) | **High** |
| PCA dims | Minor difference; 30-50 is fine | **Low** impact |
| Resolution selection | Combine: Hannah's dual-metric scoring + Workflow B's adaptive scaling | **High** |
| Marker detection | Both needed; pseudobulk DE (already implemented) for conditions | **High** |
| Tiered integration | Keep (Workflow B, IVD-specific) | **High** |
| Framework | Follow team expertise; R has less friction for downstream | **Medium** |

---

## Proposed Validation Experiment

The v4 pipeline has already been run end-to-end, so we have baseline results to compare against. The validation question is: **does switching normalization and/or integration method improve the results?**

**Option 1: Targeted normalization test (fastest)**

The clearest gap between the workflows is normalization (SCTransform vs. log-normalization). Since scANVI models raw counts, the normalization primarily affects HVG selection and per-dataset QC/visualization. Test whether SCTransform-derived HVGs produce a better scANVI integration:

1. Take one multi-sample IVD dataset (e.g., GSE230809, 24 samples across NP+AF)
2. Run SCTransform per sample in R, extract the top 3000 variable features
3. Feed those HVGs (plus raw counts) into the existing scANVI pipeline
4. Compare integration metrics against the current v4 results for this dataset

**Option 2: Full head-to-head comparison (more thorough)**

1. **Select one IVD dataset** from our own lab (e.g., from GSE230809 or whichever dataset Hannah's team can validate against known results)
2. **Run both workflows** on this dataset:
   - Workflow A: SCTransform → CCA (adapted for single study, multi-sample)
   - Workflow B: Current v4 pipeline (normalize_total + log1p → scANVI)
   - Workflow C (hybrid): SCTransform → scANVI (best normalization + best integration for multi-platform)
3. **Compare:**
   - Cluster number and composition
   - Known cell type recovery (do both find the expected populations?)
   - Marker gene overlap
   - UMAP structure
   - Integration metrics (iLISI, batch-ASW, condition-ASW)
4. **Use the results** to make an informed choice, then scale to the full atlas

This comparison also tests whether the integration challenges are primarily a normalization issue or an integration method issue — or both.

---

## Questions for Discussion

1. Does Hannah's team have a dataset where the "ground truth" cell types are well-established, so we can validate against known biology?
2. The v4 pipeline has produced 19 cell types and stable findings across 4 versions. **What specifically are the integration challenges** that motivated this review? Are there specific cell types that aren't separating, studies that don't mix, or biological signals that seem lost? Understanding the specific failure modes will help us target the right fix rather than overhauling the entire pipeline.
3. Is the team comfortable maintaining an R-based pipeline, or should we aim for a hybrid (R for preprocessing, Python for integration/downstream)? A hybrid approach (SCTransform in R → export counts/HVGs → scANVI in Python) could combine the best of both.
4. For the tiered integration: should we tier within CCA as well (separate FindIntegrationAnchors for mesenchymal vs. non-mesenchymal)?
5. The IVD atlas includes 3 non-10x platforms (BD Rhapsody, Singleron). Does Hannah's team have experience integrating across platforms with CCA, or has the synovium work been 10x-only? This is a key consideration — scANVI's `dispersion='gene-batch'` explicitly models platform-specific technical effects, while CCA treats all batches identically.
6. Should we consider Seurat v5's newer integration methods (e.g., IntegrateLayers with CCA or Harmony) which scale better than v4's pairwise anchor approach?
7. The v4 pipeline already implements silhouette-based resolution optimization. Would adopting Hannah's combined silhouette + modularity scoring (from `Silhouette_Modularity_Scoring.R`) as an upgrade to the existing approach be a quick win everyone agrees on?
