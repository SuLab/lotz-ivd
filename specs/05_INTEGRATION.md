# Module 05: Integration, Clustering, and De Novo Annotation

## Objective

Integrate cells across studies into shared representations, cluster the integrated objects, and annotate cell types de novo. This module produces the final cell atlas — cell type identities emerge from the integrated data rather than being imposed beforehand.

Four integrated objects are produced, each clustered and annotated independently:
1. **NP** — nucleus pulposus cells from studies with clearly separated NP tissue
2. **AF** — annulus fibrosus cells from studies with clearly separated AF tissue
3. **CEP** — cartilaginous endplate cells
4. **All-cells** — all IVD cells together, including studies where compartments were not separated

The all-cells object serves as a whole IVD atlas. IVD subcompartments (NP, AF, CEP) are related tissues with overlapping cell compositions — not distinct organs — so a combined analysis is biologically meaningful for understanding the full IVD cellular landscape and cross-compartment relationships.

## Rationale

Module 04 assigns coarse anchor labels (5 categories + Unknown) that are reliable across datasets but deliberately avoid fine-grained distinctions. This module uses those anchors for semi-supervised integration (scANVI), then clusters the result and discovers fine cell types from the integrated data. The coarse labels guide integration without imposing the final atlas identity.

Within each object, mesenchymal and non-mesenchymal cells are integrated separately (tiered integration) because they have fundamentally different transcriptomic profiles. Integrating them together would force the model to spend latent capacity on that major axis of variation rather than on subtler differences within each group.

## Study and Sample Assignments

### Excluded from all objects

- **GSE233666 (Guo 2023):** Excluded — contains only herniated discs, which is not the focus of this project (degeneration).
- **GSE189916 neonatal samples (Jiang 2022):** Excluded — neonatal disc biology is fundamentally different from adult. Only the adult samples from GSE189916 are retained.

### Object 1: NP-only

Studies where NP tissue was clearly defined and surgically separated:

| Study | Accession | NP Samples | Conditions | Notes |
|-------|-----------|------------|------------|-------|
| Gan 2021 | GSE160756 | NP portion | Healthy young/adult | 3-compartment atlas; NP surgically separated |
| Tu 2022 | GSE165722 | 8 | Pfirrmann II-V | BD Rhapsody platform |
| Cherif 2022 | GSE199866 | NP portion | Paired degen/non-degen | Same-patient paired design |
| Li 2022 | GSE205535 | 2 | Normal vs degenerative | BD Rhapsody; corrigenda exist |
| Chen 2024 | GSE244889 | 7 | Mild vs severe | Serglycin/fibrotic NP |
| Jia 2024 | GSE251686 | 5 (of 6) | Mild vs severe | NP3 excluded (corrupt matrix) |
| Swahn 2024 | GSE230809 | NP portion | Healthy vs diseased | Largest study; NP surgically separated |
| Han 2022 | CNP0002664 | 6 | Normal/mild/severe | Singleron platform |

**Not included in NP-only:** GSE189916 (whole IVD, compartments not separated — adult samples go to all-cells object only).

### Object 2: AF-only

Studies where AF tissue was clearly defined and surgically separated:

| Study | Accession | AF Samples | Conditions | Notes |
|-------|-----------|------------|------------|-------|
| Gan 2021 | GSE160756 | AF portion | Healthy young/adult | 3-compartment atlas |
| Cherif 2022 | GSE199866 | Inner AF portion | Paired degen/non-degen | Inner AF only |
| Swahn 2024 | GSE230809 | AF portion | Healthy vs diseased | AF surgically separated |

**Note:** AF coverage is limited to 3 studies. Clustering resolution and annotation granularity should be adjusted accordingly.

### Object 3: CEP-only

| Study | Accession | CEP Samples | Conditions | Notes |
|-------|-----------|-------------|------------|-------|
| Gan 2021 | GSE160756 | CEP portion | Healthy young/adult | 3-compartment atlas |
| Shi 2024 | GSE255768 | 2 | Degenerative endplate | No healthy control |
| Kuchynsky 2024 | GSE242443 | 2 | Non-degen vs degen | **Culture-expanded cells** |

**Note:** CEP coverage is sparse and includes culture-expanded cells (GSE242443). Results require strong caveats.

### Object 4: All-cells

All studies and samples from Objects 1-3, plus:

| Study | Accession | Samples | Conditions | Notes |
|-------|-----------|---------|------------|-------|
| Jiang 2022 | GSE189916 | Adult only (3) | Adult whole IVD | Compartments not separated; neonatal excluded |

This object contains all cells across all compartments for a unified IVD atlas view.

### Sample-level exclusions (across all objects)

- **GSE251686 NP3:** Corrupt matrix file on GEO (verified on re-download)
- **GSE205535 NNP:** Included in integration but excluded from DE analysis (Module 06) — 11yo spinal cord injury is a trauma confound

## Required Output: Inclusion Summary Table

Generate `results/integration/inclusion_summary.tsv` containing, for each object × study combination:
- Object name (NP, AF, CEP, all_cells)
- Study accession
- First author and year
- Number of samples included
- Number of cells included (post-QC)
- Compartment
- Conditions represented
- Platform

Also generate `results/integration/inclusion_summary.html` with a formatted version suitable for a manuscript supplementary table.

## Required Output: Study Caveats Table

Generate `results/integration/study_caveats.tsv` documenting per-study caveats for the manuscript supplement:

| Study | Caveat | Impact | Mitigation |
|-------|--------|--------|------------|
| GSE165722 (Tu 2022) | BD Rhapsody platform (not 10x) | Different capture efficiency, gene detection | Platform-aware batch correction via scVI study-level batch key |
| GSE205535 (Li 2022) | BD Rhapsody platform; published corrigenda | See above; potential data quality issues | Monitor for outlier behavior in integration |
| CNP0002664 (Han 2022) | Singleron Matrix platform (not 10x) | Different capture efficiency | Same as above |
| GSE242443 (Kuchynsky 2024) | Culture-expanded CEP cells | Culture alters cell states; may not reflect in vivo biology | Caveat in all CEP results; compare with non-expanded CEP from GSE160756 |
| GSE255768 (Shi 2024) | Degenerative endplate only; no healthy control | Cannot do healthy vs. degenerated comparison for this study alone | Healthy CEP baseline from GSE160756 |
| GSE230809 (Swahn 2024) | All-male donors; age-disease confounded | Cannot separate age from degeneration effects | Note in interpretation; sex-specific effects cannot be assessed |
| GSE205535 NNP sample | 11yo spinal cord injury, classified as "healthy" | Trauma confound | Excluded from DE comparisons |
| GSE189916 (Jiang 2022) | Whole IVD (compartments not separated) | Cannot assign cells to NP/AF/CEP | Included only in all-cells object |

## Inputs

- Processed AnnData objects from `data/processed/{accession}.h5ad` with `coarse_label` and `cell_class` from Module 04
- `metadata/sample_metadata.tsv`

## Outputs

- `data/integrated/NP.h5ad` — integrated NP cells
- `data/integrated/AF.h5ad` — integrated AF cells
- `data/integrated/CEP.h5ad` — integrated CEP cells
- `data/integrated/all_cells.h5ad` — integrated all IVD cells
- `data/integrated/integration_metrics.tsv` — quantitative integration assessment
- `results/integration/inclusion_summary.tsv` — study × object inclusion table
- `results/integration/inclusion_summary.html` — formatted version for manuscript supplement
- `results/integration/study_caveats.tsv` — per-study caveats for supplement
- `results/integration/clustering_resolution_optimization/` — resolution selection diagnostics per object
- `results/integration/cluster_markers/` — DE markers per cluster per object
- `results/integration/annotation_dotplots/` — canonical marker expression per cluster per object
- `results/integration_report.html` — visualization of integration, clustering, and annotation results
- Updated `analysis_plan.md` with chosen parameters and cell type definitions

### Notebook: `notebooks/05_integration.ipynb`

Contains:
- Inclusion summary table (studies × objects)
- UMAP panels per object, colored by: study, cluster, cell type (annotated), condition, compartment
- Integration quality metrics (batch mixing, biological conservation)
- Clustering resolution optimization plots (silhouette, modularity, clustree) per object
- Dot plots of canonical markers per cluster per object (for annotation support)
- Cluster DE marker heatmaps per object
- Final cell type definitions table per object
- Study caveats table

**Manuscript mapping:** Figure 1: IVD cell atlas (UMAP, dot plot, proportions). Supplementary Table S1: Study inclusion and caveats. Supplementary Figure S3: Integration benchmarking and resolution optimization. Methods section on integration, clustering, and annotation.

## Part 1: Tiered Integration with scANVI

For each of the four objects (NP, AF, CEP, all-cells), perform tiered integration using scANVI with the coarse anchor labels from Module 04.

### Integration method: scANVI (semi-supervised)

scANVI uses labeled cells as anchors to align similar cells across studies while leaving "Unknown" cells free to be positioned by transcriptomic similarity. This produces better batch correction than unsupervised scVI, especially across platforms (10x, BD Rhapsody, Singleron).

**Workflow per tier:**
1. Train scVI first (unsupervised) to learn a base latent representation
2. Initialize scANVI from the trained scVI model
3. Feed `coarse_label` as `labels_key`, with "Unknown" as `unlabeled_category`
4. scANVI refines the latent space using the anchor labels

**Parameters:**
- `batch_key='study'`
- `n_latent=20`
- scVI: `max_epochs=200`
- scANVI: `max_epochs=50`, `early_stopping=True`, `early_stopping_patience=5`
- `n_top_genes=3000`

### Tier A: Mesenchymal cells

**Anchor labels used:** Chondrocyte_like, Fibroblast_like, Unknown
**scANVI unlabeled_category:** "Unknown"

**Steps:**
1. Subset cells with `cell_class == 'mesenchymal'` or `cell_class == 'unknown'` for the relevant studies/samples (include Unknown-class cells that may be fibrochondrocytes or other transitional mesenchymal types)
2. Concatenate across datasets
3. Re-identify HVGs on the concatenated object (n_top_genes=3000)
4. Train scVI → initialize scANVI with coarse_label anchors
5. Extract scANVI latent representation
6. Compute neighbors and UMAP on scANVI embedding
7. Proceed to clustering (Part 2) and annotation (Part 3)

### Tier B: Non-mesenchymal cells

**Anchor labels used:** Immune, Endothelial, Pericyte_SMC
**scANVI unlabeled_category:** "Unknown" (should be very few in this tier)

Same workflow as Tier A. These populations have strong, discrete transcriptomic identities — scANVI should integrate them cleanly.

### Merging tiers

After independent integration, clustering, and annotation of each tier, merge them back into a single AnnData per object for downstream analysis. Store tier-specific embeddings in `obsm` (e.g., `X_scanvi_mesenchymal`, `X_scanvi_non_mesenchymal`).

### Integration Quality Metrics

For each integrated object/tier, compute:

1. **Batch mixing:**
   - iLISI (integration LISI) — higher is better mixing
   - Batch-ASW — should be near 0

2. **Biological conservation:**
   - Condition-ASW — conditions should remain distinguishable
   - Condition classifier accuracy — can we still distinguish healthy from degenerated?

3. **Continuum preservation (mesenchymal only):**
   - Compare cluster count at resolution 0.5 before and after integration (large reduction suggests overcorrection)
   - Verify that the integrated object does not collapse into a single blob

## Part 2: Clustering with Resolution Optimization

After integration, cluster each tier within each object independently.

### Method

Leiden clustering on the scANVI neighbor graph, with resolution selected by multi-metric optimization.

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
5. If metrics disagree, prefer the resolution that gives biologically interpretable clusters (assessed in Part 3)

Store clustering results at the selected resolution in `obs['leiden']`, and also store results at 2-3 other resolutions for comparison (e.g., `obs['leiden_0.5']`, `obs['leiden_1.0']`).

### Expectations per object

- **NP:** Multiple chondrocyte-like clusters varying along the notochordal → mature → stressed/degenerative continuum
- **AF:** Clusters along the inner (chondrocyte-like) → outer (fibroblast-like) gradient, plus mechanical stress states
- **CEP:** Hyaline cartilage-like and potentially ossifying chondrocytes (limited data — expect fewer clusters)
- **All-cells:** Should recover all of the above plus cross-compartment structure; some clusters may contain cells from multiple compartments (reflecting shared biology)
- **Non-mesenchymal (all objects):** Discrete immune subtypes (macrophages, T cells, B cells, etc.) + endothelial + pericyte

## Part 3: Two-Stage De Novo Cell Type Annotation

Annotate clusters in two stages: first assign a coarse identity, then refine within each coarse group. Do NOT rely on the original manuscript annotations. The Module 04 coarse labels served integration only — annotation here is driven by the integrated data.

### Stage 1: Coarse annotation via canonical markers

For each cluster, check expression of canonical marker genes to assign a coarse identity. This determines *what* the cell is.

**Canonical marker panels for coarse assignment:**

*Mesenchymal — Chondrocyte-like:*
- COL2A1, ACAN, SOX9, COMP, PRG4

*Mesenchymal — Fibroblast-like:*
- COL1A1, COL1A2, DCN, LUM, THY1

*Mesenchymal — Fibrochondrocyte-like:*
- Co-expression of both chondrocyte markers (COL2A1, ACAN) and fibroblast markers (COL1A1, DCN) at moderate levels

*Non-mesenchymal — Immune:*
- Macrophage: CD68, CD14, CSF1R, CD163 (M2), CD86 (M1)
- T cell: CD3D, CD3E, CD4, CD8A
- B cell: CD79A, MS4A1
- NK cell: NKG7, GNLY
- Mast cell: KIT, TPSAB1

*Non-mesenchymal — Endothelial:*
- PECAM1, VWF, CDH5

*Non-mesenchymal — Pericyte/SMC:*
- ACTA2, RGS5, PDGFRB

**Visualizations for Stage 1:**
- Dot plots: canonical markers × clusters (fraction expressing + mean expression)
- Feature plots (UMAPs): key markers overlaid on clusters

Each cluster receives a `coarse_cell_type` label (e.g., "Chondrocyte-like", "Fibroblast-like", "Fibrochondrocyte-like", "Macrophage", "T_cell", "Endothelial", etc.).

### Stage 2: Fine annotation via DE markers within coarse groups

Within each coarse category, compute DE genes between the clusters of that category to find what distinguishes them. This determines *what state or subtype* the cell is in.

For example, if there are 5 chondrocyte-like clusters:
1. Run `sc.tl.rank_genes_groups()` comparing only those 5 clusters against each other
2. Extract the top 20-50 DE markers per cluster
3. Use these to assign subtype/state labels

**Expected fine distinctions within coarse groups:**

*Within Chondrocyte-like (NP):*
- Notochordal-like (T/TBXT, SHH, NOG, CD24, KRT8, KRT18, KRT19)
- Mature chondrocyte (high ACAN, COL2A1, COMP)
- Stressed/degenerative (MMP13, ADAMTS5, IL1B, TNF, VEGFA, HIF1A)
- Hypertrophic (COL10A1, RUNX2)

*Within Fibroblast-like (AF):*
- Inner AF (more chondrocyte-adjacent: SOX9+)
- Outer AF (classic fibroblast: high COL1A1, COL1A2)
- Mechanical stress (COMP, CILP, THBS1)

*Within Fibrochondrocyte-like:*
- May subdivide by the ratio of chondrocyte-to-fibroblast markers, or by stress/degenerative markers

*Within Immune:*
- M1 vs. M2 macrophage, CD4 vs. CD8 T cells, etc.

*Within Endplate chondrocytes:*
- Hyaline cartilage (COL2A1, COL10A1, SOX9)
- Ossification (RUNX2, SP7/OSX, BGLAP)

**Visualizations for Stage 2:**
- Heatmap: top DE markers per cluster within each coarse group
- Dot plots: fine markers within coarse groups

### Final labels

Store annotations in `obs`:
- `coarse_cell_type` — from Stage 1 (e.g., "Chondrocyte-like", "Macrophage")
- `cell_type` — from Stage 2, combining coarse identity and fine distinction (e.g., "NP_notochordal", "NP_mature_chondrocyte", "NP_stressed", "AF_outer_fibroblast", "Macrophage_M2")

### Annotation procedure

For each cluster:
1. Check canonical marker expression (Stage 1) → assign coarse_cell_type
2. Within each coarse group, review DE markers between clusters (Stage 2) → assign cell_type
3. If a cluster's fine identity is ambiguous, check whether it splits at a higher clustering resolution into identifiable subtypes
4. If multiple clusters have the same cell_type, consider whether they represent genuine subtypes (different marker profiles) or should be merged (same markers, split by batch)

### Annotation confidence

For each cluster annotation, record:
- `cell_type` — assigned label
- `cell_type_confidence` — high (clear markers), medium (consistent but not definitive), low (ambiguous)
- `annotation_evidence` — brief note on which markers/DE genes support the label

### Continuous scores

For mesenchymal clusters that exist on a continuum (e.g., notochordal → mature → degenerative in NP), also compute continuous gene signature scores using `sc.tl.score_genes()`:
- `score_notochordal`, `score_degenerative`, `score_fibrotic`, etc.
These complement the discrete labels and preserve continuum information for downstream trajectory analysis.

## Automated Validation

### Integration
- [ ] All four integrated objects are saved (NP, AF, CEP, all_cells)
- [ ] Inclusion summary table is generated
- [ ] Study caveats table is generated
- [ ] No integration result collapses all cells into a single cluster at resolution 0.5 (blob check)
- [ ] No integration result has study identity perfectly predicting cluster identity (ARI < 1.0)
- [ ] Integration metrics are recorded per object

### Clustering
- [ ] Resolution optimization plots are generated per object (silhouette, modularity, clustree)
- [ ] Selected resolution is documented with rationale per object
- [ ] Clustering results are stored at the selected resolution and 2-3 comparison resolutions

### Annotation
- [ ] All clusters have a cell type label (including "unassigned" if ambiguous)
- [ ] Proportion of "unassigned" cells is < 10% per object
- [ ] Cluster DE marker tables are generated per object
- [ ] Canonical marker dot plots are generated per object
- [ ] For non-mesenchymal: PTPRC expressed predominantly in immune clusters, PECAM1 in endothelial
- [ ] For mesenchymal: expected compartment markers are expressed in the corresponding clusters
- [ ] Annotation evidence is recorded for each cluster

## Human Checkpoint

This is the most critical decision point in the entire pipeline — it defines the cell atlas.

### Review materials
- Inclusion summary table (which studies/samples in which objects)
- Study caveats table
- UMAP per object colored by cluster and annotated cell type
- Integration quality metrics per object
- Resolution optimization plots per object
- Cluster DE marker heatmaps per object
- Canonical marker dot plots per object
- Cell type proportions per dataset within each object (to check for study-specific artifacts)
- Annotation confidence table per object

### Questions for the reviewer
1. Does the clustering resolution capture biologically meaningful groups without over-splitting?
2. Do the cell type annotations make biological sense?
3. Are any clusters clearly batch-driven rather than biology-driven?
4. For the mesenchymal continuum: are discrete labels appropriate, or should some clusters be merged?
5. Are there unexpected cell types or missing expected types?
6. Is the all-cells object consistent with the compartment-specific objects? Do the same cell types appear?
7. Are the study caveats adequately documented for the supplement?
8. Should any annotations be revised before proceeding to differential analysis?

### Potential plan revisions
- If the optimal resolution differs substantially between metrics, test both downstream and report sensitivity
- If batch effects dominate certain clusters, consider excluding those clusters or adjusting integration parameters
- If the mesenchymal continuum resists discrete clustering, rely more heavily on continuous scores for downstream analysis
- If CEP results are unreliable due to culture expansion artifacts, note this prominently and consider excluding CEP from primary analyses
- Annotation decisions here directly affect Modules 06-09 — any changes require careful propagation
