# Module 10: Trajectory & Dynamics Analysis

## Objective

Characterize the continuum of IVD cell states using trajectory inference and RNA velocity. Determine whether there are directional transitions between cell states, how these relate to degeneration, and whether aging and disease progression follow predictable paths through cell state space.

## Rationale

The chondrocyte-fibroblast continuum in IVD is a key feature of the biology. Rather than forcing discrete clusters, trajectory analysis can model this continuum explicitly and identify gene expression programs that vary along it. RNA velocity can add directionality — are cells transitioning from a healthy state toward a degenerative state?

## Inputs

- Integrated AnnData objects from Module 05: `data/integrated/NP.h5ad`, `data/integrated/AF.h5ad`
- Cell type labels from Module 07 post-integration de novo annotation (`obs['cell_type']`)
- Raw count matrices (needed for RNA velocity)

## Outputs

- `results/trajectories/pseudotime_{compartment}.tsv` — pseudotime values per cell
- `results/trajectories/velocity_{compartment}.h5ad` — RNA velocity results
- `results/trajectories/trajectory_genes_{compartment}.tsv` — genes that vary along trajectories
- `results/trajectories/trajectory_report.html`

### Notebook: `notebooks/10_trajectory.ipynb`

Produced after trajectory and velocity analyses. Contains:

*Trajectory section:*
- PAGA graph colored by cell type, overlaid on UMAP
- UMAP colored by pseudotime, with condition overlay (side-by-side or contour plot)
- Pseudotime distribution by condition (density plot or violin) — does degeneration shift cells along the trajectory?
- Pseudotime vs. age scatter (if age data is available per cell/sample)
- Heatmap: top trajectory-associated genes ordered by pseudotime, grouped by gene module
- Gene expression curves along pseudotime for key IVD genes (ACAN, COL2A1, MMP13, etc.)

*Velocity section (if available):*
- UMAP with velocity stream arrows
- Latent time distribution by condition
- Velocity confidence map
- Key velocity genes

*Gene program section:*
- Pathway enrichment for each temporal gene module (early, middle, late)
- Cross-reference with DE results: Venn diagram of trajectory genes vs. DE genes

**Manuscript mapping:** Figure 6: Cell state trajectory in NP (PAGA, pseudotime UMAP, gene dynamics). Supplementary: velocity results, AF trajectory if applicable.

## Part 1: Trajectory Inference (Pseudotime)

### Scope

Run trajectory analysis on mesenchymal IVD cells only (not immune/endothelial — those are discrete populations, not a continuum). The integrated mesenchymal objects from Module 05 (NP and AF separately) are already the correct input — no further subsetting needed beyond the Module 07 annotation.

### Methods

**Primary:** Diffusion pseudotime (DPT) via scanpy (`sc.tl.diffusion_pseudotime`)

**Alternative:** PAGA (Partition-based Graph Abstraction) to establish coarse-grained connectivity between clusters before computing pseudotime. PAGA → DPT is a well-validated workflow in scanpy.

**Steps:**
1. Subset to resident cells in the target compartment
2. Recompute HVGs, PCA, and neighbor graph on the subset
3. Run PAGA to identify cluster connectivity
4. Select root cell(s) — options:
   - If neonatal/healthy cells are present, use the cluster annotated as notochordal-like (NP) or inner AF (AF) from Module 07 de novo annotation
   - Alternatively, use the cluster with highest `score_notochordal` or highest proportion of healthy/neonatal cells
   - If no clear developmental root, test multiple roots and compare
5. Compute DPT from root
6. Identify genes that correlate with pseudotime (Spearman correlation, FDR correction)
7. Group trajectory-associated genes into early, middle, late programs

### Interpretation

- Does pseudotime correlate with condition (healthy → degenerated)?
- Does pseudotime correlate with age?
- What gene programs activate or deactivate along the trajectory?
- Is there a branch point where cells diverge into distinct fates (e.g., fibrotic vs. hypertrophic)?

## Part 2: RNA Velocity

### Prerequisites

RNA velocity requires spliced and unspliced read counts, which must be quantified from the BAM files or fastqs. Many public datasets only provide the final count matrix without spliced/unspliced information.

**Check for each dataset:**
1. Are BAM files or fastqs available? → Can run velocyto or STARsolo to quantify spliced/unspliced
2. Did the original study provide spliced/unspliced layers (some h5ad files include these)?
3. If neither is available, RNA velocity cannot be computed for that dataset. Log and skip.

### Method

**Tool:** scVelo (stochastic or dynamical mode)

**Steps:**
1. For datasets where spliced/unspliced counts are available:
   a. Load into AnnData with layers `spliced` and `unspliced`
   b. Filter and normalize as per scVelo documentation
   c. Compute moments
   d. Run velocity estimation (dynamical mode preferred, stochastic as fallback)
   e. Project velocity onto the UMAP embedding
   f. Compute latent time
2. Identify genes with significant velocity (velocity_genes)
3. Overlay velocity arrows on UMAP — do they indicate directionality in the cell state continuum?

### Caveats

- RNA velocity has known limitations in slow-turnover tissues like cartilage. IVD cell turnover is very low. Velocity results should be interpreted cautiously.
- If velocity arrows are random/incoherent, this is an expected negative result for IVD, not necessarily a pipeline failure. Report it as such.
- Dynamical mode is more robust than stochastic but slower.

## Part 3: Gene Programs Along Trajectories

### Method

For genes that significantly correlate with pseudotime or latent time:

1. Smooth gene expression curves along pseudotime (e.g., using a rolling window or spline fit)
2. Cluster genes by their smoothed expression pattern using k-means to identify gene modules: sets of genes with correlated temporal dynamics
3. Run pathway enrichment on each module
4. Cross-reference with DE results from Module 08 — are the same genes identified as trajectory-associated and DE between conditions?

This provides a richer view than DE alone: not just "gene X is up in degeneration" but "gene X activates at a specific point along the degenerative trajectory."

## Automated Validation

- [ ] Pseudotime results exist for NP (and AF if sufficient data)
- [ ] Pseudotime correlates with at least one known biological variable (condition, age, or continuous cell state score) — if it doesn't correlate with anything, it may be capturing batch variation instead
- [ ] If RNA velocity was run, velocity confidence scores are reasonable (not uniformly low)
- [ ] Trajectory gene lists are generated
- [ ] Report HTML is generated
- [ ] If RNA velocity couldn't be run (no spliced/unspliced data), this is clearly documented

## Human Checkpoint

### Review materials
- PAGA graphs showing cluster connectivity
- UMAP colored by pseudotime, with condition overlay
- RNA velocity stream plots (if available)
- Top trajectory-associated genes and their expression patterns
- Gene module pathway enrichment results

### Questions for the reviewer
1. Does the inferred trajectory make biological sense? Is the root cell choice appropriate?
2. Does pseudotime align with the expected healthy → degenerated axis, or is it capturing something else (e.g., batch, compartment)?
3. Are the velocity results coherent or noisy? Should they be included in the final analysis?
4. Do the gene programs along the trajectory reveal a staged degenerative process, or is it more like a smooth gradient?
5. Are there branch points suggesting divergent cell fates?
6. Should trajectory findings feed back into the cell type annotation (e.g., defining trajectory-based cell states)?

### Potential plan revisions
- If trajectory analysis reveals a clear staged process, this restructures the biological narrative and should be central to the final report
- If velocity is incoherent (expected for IVD), de-emphasize it in the final analysis and rely on pseudotime + static DE instead
- If trajectory-based cell states are more biologically meaningful than cluster-based annotations, consider redefining cell type labels and rerunning DE with the new groupings
