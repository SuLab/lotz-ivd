# Module 07: Post-Integration Cell Type Annotation

## Objective

Annotate clusters from Module 06 with cell type identities using a two-stage approach: first assign a coarse identity via canonical markers, then refine within each coarse group using DE markers. This module produces the final cell atlas — cell type identities emerge from the integrated data rather than being imposed beforehand.

## Rationale

Module 04 assigned coarse anchor labels (5 categories + Unknown) that guided scANVI integration in Module 05. Those labels served integration only. Annotation here is driven by the integrated, clustered data — it discovers fine cell types de novo.

## Inputs

- Clustered AnnData objects from Module 06: `data/integrated/NP.h5ad`, `data/integrated/AF.h5ad`, `data/integrated/CEP.h5ad`, `data/integrated/all_cells.h5ad`
- Each object contains clustering results in `obs['leiden']`

## Outputs

- Updated AnnData objects with annotation in `obs['coarse_cell_type']`, `obs['cell_type']`, `obs['cell_type_confidence']`, `obs['annotation_evidence']`
- `results/integration/cluster_markers/` — DE markers per cluster per object
- `results/integration/annotation_dotplots/` — canonical marker expression per cluster per object
- `results/integration/cell_type_definitions.tsv` — final cell type definitions with markers and evidence
- Updated `analysis_plan.md` with cell type definitions

### Notebook: `notebooks/07_annotation.ipynb`

Contains:
- Dot plots of canonical markers per cluster per object (for annotation support)
- Cluster DE marker heatmaps per object
- UMAP panels per object colored by: coarse_cell_type, cell_type
- Final cell type definitions table per object
- Cell type proportions per dataset within each object

**Manuscript mapping:** Figure 1: IVD cell atlas (UMAP, dot plot, proportions). Methods section on annotation.

## Stage 1: Coarse Annotation via Canonical Markers

For each cluster, check expression of canonical marker genes to assign a coarse identity. This determines *what* the cell is.

### Canonical marker panels for coarse assignment

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

### Visualizations for Stage 1
- Dot plots: canonical markers × clusters (fraction expressing + mean expression)
- Feature plots (UMAPs): key markers overlaid on clusters

Each cluster receives a `coarse_cell_type` label (e.g., "Chondrocyte-like", "Fibroblast-like", "Fibrochondrocyte-like", "Macrophage", "T_cell", "Endothelial", etc.).

## Stage 2: Fine Annotation via DE Markers Within Coarse Groups

Within each coarse category, compute DE genes between the clusters of that category to find what distinguishes them. This determines *what state or subtype* the cell is in.

For example, if there are 5 chondrocyte-like clusters:
1. Run `sc.tl.rank_genes_groups()` comparing only those 5 clusters against each other
2. Extract the top 20-50 DE markers per cluster
3. Use these to assign subtype/state labels

### Expected fine distinctions within coarse groups

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

### Visualizations for Stage 2
- Heatmap: top DE markers per cluster within each coarse group
- Dot plots: fine markers within coarse groups

## Final Labels

Store annotations in `obs`:
- `coarse_cell_type` — from Stage 1 (e.g., "Chondrocyte-like", "Macrophage")
- `cell_type` — from Stage 2, combining coarse identity and fine distinction (e.g., "NP_notochordal", "NP_mature_chondrocyte", "NP_stressed", "AF_outer_fibroblast", "Macrophage_M2")

## Annotation Procedure

For each cluster:
1. Check canonical marker expression (Stage 1) → assign coarse_cell_type
2. Within each coarse group, review DE markers between clusters (Stage 2) → assign cell_type
3. If a cluster's fine identity is ambiguous, check whether it splits at a higher clustering resolution into identifiable subtypes
4. If multiple clusters have the same cell_type, consider whether they represent genuine subtypes (different marker profiles) or should be merged (same markers, split by batch)

## Annotation Confidence

For each cluster annotation, record:
- `cell_type` — assigned label
- `cell_type_confidence` — high (clear markers), medium (consistent but not definitive), low (ambiguous)
- `annotation_evidence` — brief note on which markers/DE genes support the label

## Continuous Scores

For mesenchymal clusters that exist on a continuum (e.g., notochordal → mature → degenerative in NP), also compute continuous gene signature scores using `sc.tl.score_genes()`:
- `score_notochordal`, `score_degenerative`, `score_fibrotic`, etc.
These complement the discrete labels and preserve continuum information for downstream trajectory analysis.

## Automated Validation

- [ ] All clusters have a cell type label (including "unassigned" if ambiguous)
- [ ] Proportion of "unassigned" cells is < 10% per object
- [ ] Cluster DE marker tables are generated per object
- [ ] Canonical marker dot plots are generated per object
- [ ] For non-mesenchymal: PTPRC expressed predominantly in immune clusters, PECAM1 in endothelial
- [ ] For mesenchymal: expected compartment markers are expressed in the corresponding clusters
- [ ] Annotation evidence is recorded for each cluster
- [ ] Cell type definitions table is generated

## Human Checkpoint

This is the most critical decision point in the entire pipeline — it defines the cell atlas.

### Review materials
- UMAP per object colored by cluster and annotated cell type
- Cluster DE marker heatmaps per object
- Canonical marker dot plots per object
- Cell type proportions per dataset within each object (to check for study-specific artifacts)
- Annotation confidence table per object
- Cell type definitions table

### Questions for the reviewer
1. Do the cell type annotations make biological sense?
2. For the mesenchymal continuum: are discrete labels appropriate, or should some clusters be merged?
3. Are there unexpected cell types or missing expected types?
4. Is the all-cells object consistent with the compartment-specific objects? Do the same cell types appear?
5. Should any annotations be revised before proceeding to differential analysis?

### Potential plan revisions
- If the mesenchymal continuum resists discrete clustering, rely more heavily on continuous scores for downstream analysis
- If CEP results are unreliable due to culture expansion artifacts, note this prominently and consider excluding CEP from primary analyses
- Annotation decisions here directly affect Modules 08-11 (DE, interpretation, trajectory, communication) — any changes require careful propagation
