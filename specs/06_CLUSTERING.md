# Module 06: Clustering with Resolution Optimization

## Objective

Cluster the integrated objects from Module 05 to identify cell populations. Each tier (mesenchymal and non-mesenchymal) within each object (NP, AF, CEP, all-cells) is clustered independently, with resolution selected by multi-metric optimization.

## Inputs

- Integrated objects from the selected Module 05 workflow: `data/integrated/{wf}/NP.rds`, `data/integrated/{wf}/AF.rds`, `data/integrated/{wf}/CEP.rds`, `data/integrated/{wf}/all_cells.rds` (where `{wf}` is `cca`, `scanvi`, or `stacas` — determined at Module 05 human checkpoint)
- If Workflow A (CCA) was selected: single integrated embedding per object (CCA reduction)
- If Workflow B (scANVI) was selected: two tier embeddings per object (`X_scanvi_mesenchymal`, `X_scanvi_non_mesenchymal`)
- If Workflow C (STACAS) was selected: single integrated embedding per object (STACAS reduction)

## Outputs

- Updated AnnData objects with clustering results in `obs['leiden']` (selected resolution) and comparison resolutions
- `results/integration/clustering_resolution_optimization/` — resolution selection diagnostics per object per tier
- Updated `analysis_plan.md` with chosen resolutions and rationale

### Notebook: `notebooks/06_clustering.ipynb`

Contains:
- Clustering resolution optimization plots (silhouette, modularity, clustree) per object per tier
- UMAP panels per object, colored by: study, cluster, condition
- Cluster size distributions per object
- Selected resolution justification

**Manuscript mapping:** Supplementary Figure S3 (continued): Resolution optimization diagnostics.

## Method

Leiden clustering on the integration neighbor graph (from whichever workflow was selected in Module 05), with resolution selected by multi-metric optimization.

### Resolution optimization

Test resolutions from 0.1 to 2.0 in steps of 0.1, and compute:

1. **Silhouette score** (on the scANVI embedding): measures how well-separated clusters are. Higher is better, but will decrease at very high resolutions as clusters become too granular.

2. **Modularity score**: measures the quality of the graph partition. Computed by the Leiden algorithm itself.

3. **Clustree stability analysis**: visualize how clusters split/merge across resolutions using `clustree` (R package) or a Python equivalent. Identify the resolution where clusters are stable (cells don't constantly reshuffle between clusters). Resolutions above which clusters fragment into noise should be avoided.

### Resolution selection logic

1. Plot silhouette score vs. resolution — look for a peak or plateau
2. Plot modularity vs. resolution — look for the knee point (modularity drops off)
3. Examine clustree — identify the resolution above which clusters start fragmenting
4. Select the resolution that balances all three: good silhouette, high modularity, stable clusters
5. If metrics disagree, prefer the resolution that gives biologically interpretable clusters (assessed in Module 07)

Store clustering results at the selected resolution in `obs['leiden']`, and also store results at 2-3 other resolutions for comparison (e.g., `obs['leiden_0.5']`, `obs['leiden_1.0']`).

### Expectations per object

- **NP:** Multiple chondrocyte-like clusters varying along the notochordal → mature → stressed/degenerative continuum
- **AF:** Clusters along the inner (chondrocyte-like) → outer (fibroblast-like) gradient, plus mechanical stress states
- **CEP:** Hyaline cartilage-like and potentially ossifying chondrocytes (limited data — expect fewer clusters)
- **All-cells:** Should recover all of the above plus cross-compartment structure; some clusters may contain cells from multiple compartments (reflecting shared biology)
- **Non-mesenchymal (all objects):** Discrete immune subtypes (macrophages, T cells, B cells, etc.) + endothelial + pericyte

## Automated Validation

- [ ] Clustering results are stored at the selected resolution and 2-3 comparison resolutions per object per tier
- [ ] Resolution optimization plots are generated per object per tier (silhouette, modularity, clustree)
- [ ] Selected resolution is documented with rationale per object per tier
- [ ] No clustering result collapses all cells into a single cluster (blob check)
- [ ] No clustering result has study identity perfectly predicting cluster identity (ARI < 1.0)

## Human Checkpoint

### Review materials
- Resolution optimization plots per object per tier
- UMAP per object colored by cluster at selected resolution
- Cluster size distributions
- Cell type proportions per cluster per dataset (to check for study-specific artifacts)

### Questions for the reviewer
1. Does the clustering resolution capture biologically meaningful groups without over-splitting?
2. Are any clusters clearly batch-driven rather than biology-driven?
3. For the mesenchymal continuum: is the resolution appropriate, or should some clusters be merged/split?
4. Are the non-mesenchymal clusters well-separated?

### Potential plan revisions
- If the optimal resolution differs substantially between metrics, test both downstream and report sensitivity
- If batch effects dominate certain clusters, consider excluding those clusters or adjusting integration parameters
- If the mesenchymal continuum resists discrete clustering, consider relying more on continuous scores for downstream analysis
