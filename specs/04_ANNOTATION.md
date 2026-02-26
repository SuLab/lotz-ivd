# Module 04: Per-Dataset Cell Type Annotation

## Objective

Assign robust cell type and cell state labels to each cell in each dataset, using multiple complementary annotation strategies. This is done per-dataset before integration to avoid batch-correction artifacts influencing annotations.

## Rationale

The preliminary labels from Module 03 are coarse (based on a few markers per type). This module applies more rigorous annotation to resolve the chondrocyte/fibroblast continuum into meaningful subtypes and cell states. Multiple methods are used because no single approach is reliable for IVD tissue, where the dominant cell populations lack the sharp transcriptomic boundaries seen in, e.g., blood.

## Inputs

- Processed AnnData objects from `data/processed/{accession}.h5ad`
- Canonical marker gene lists (defined below)
- Reference datasets/atlases (if suitable ones exist)

## Outputs

Per dataset:
- Updated `.obs` in the h5ad file with annotation columns (see below)
- `results/annotations/{accession}_annotation_report.html`
- `results/annotations/{accession}_marker_dotplot.pdf`

Aggregate:
- `results/annotations/annotation_comparison.tsv` — cross-dataset comparison of cell type proportions
- `results/annotations/annotation_summary.html`

## Annotation Strategies

Apply all three strategies. Final labels are assigned by consensus or human judgment at the checkpoint.

### Strategy 1: Marker-based scoring

Use `sc.tl.score_genes()` or `decoupler` to score each cell for curated gene signatures. This extends the preliminary annotation with finer resolution.

**Gene signature sets:**

These should be compiled from IVD literature. Below is a starting framework — the agent should search PubMed for recent IVD scRNA-seq papers to refine these lists before applying them.

*NP subtypes:*
- Notochordal-like NP: T/TBXT, SHH, NOG, CD24, KRT8, KRT18, KRT19
- Mature NP chondrocyte: ACAN, COL2A1, SOX9, COMP, PRG4
- Stressed/degenerative NP: MMP13, ADAMTS5, IL1B, TNF, VEGFA, HIF1A
- Fibrocartilaginous NP: COL1A1 co-expressed with COL2A1, VCAN

*AF subtypes:*
- Inner AF: COL2A1, ACAN, SOX9 (chondrocyte-like)
- Outer AF: COL1A1, COL1A2, THY1, DCN, LUM (fibroblast-like)
- Mechanical stress AF: COMP, CILP, THBS1

*Endplate:*
- Hyaline cartilage: COL2A1, COL10A1, SOX9
- Ossification markers: RUNX2, SP7/OSX, BGLAP

*Non-resident:*
- Use the same marker sets from Module 03 for immune, endothelial, pericyte

### Strategy 2: Reference-based label transfer

If a suitable reference atlas exists (check CellxGene, published IVD atlases, or musculoskeletal references), use label transfer to annotate query cells.

Tools: `scvi-tools` (scANVI), `celltypist`, or `scanpy.tl.ingest()`

**Important caveats:**
- IVD-specific reference atlases may not exist or may be from a single study (circular if that study is in our dataset)
- Musculoskeletal or cartilage references may be more broadly available but less specific
- If no suitable reference exists, skip this strategy and note the gap

### Strategy 3: Automated annotation tools

Run `CellTypist` with the built-in human models (Immune_All_Low, Immune_All_High, and any available tissue-specific models). This is primarily useful for annotating the non-resident cell populations (immune subtypes, endothelial) with more precision.

For the resident IVD cells, automated tools trained on other tissues are unlikely to perform well. Log their predictions but do not rely on them for resident cell labels.

## Consensus Labeling

For each cell, produce a final label following this hierarchy:

1. **Non-resident cells (immune, endothelial, pericyte):** Use CellTypist or reference-based labels if available; fall back to marker scoring. These are typically well-resolved.

2. **Resident IVD cells:** Use marker-based scoring as the primary method. If substructure is evident in the clustering and supported by differential marker expression, assign subtype labels. If cells fall on a continuum without clear boundaries, assign a broad label (e.g., "NP_chondrocyte") and add continuous scores for relevant programs (e.g., "notochordal_score", "degenerative_score", "fibrotic_score") rather than forcing discrete categories.

3. **Ambiguous cells:** Label as "unassigned" rather than guessing. Record which strategies disagreed.

Store in `.obs`:
- `cell_type_final` — consensus label (discrete)
- `cell_type_confidence` — high / medium / low
- `cell_type_marker_based` — label from Strategy 1
- `cell_type_reference_based` — label from Strategy 2 (if applied)
- `cell_type_celltypist` — label from Strategy 3
- Continuous scores as applicable (e.g., `score_notochordal`, `score_degenerative`, etc.)

## Annotation Coherence Checks

After assigning labels, verify internal consistency:

- Cells labeled as immune should cluster together and separately from resident cells
- Cells labeled as NP chondrocyte should express NP markers and NOT express AF-specific markers (and vice versa), unless they are in a transitional state (in which case, both scores should be moderate)
- Within a single study, the proportion of immune cells should be consistent across samples from the same condition (large variation may indicate technical issues)
- Endothelial cells should be absent or very rare in studies that used cell sorting/selection for NP cells

## Automated Validation

Per dataset:
- [ ] All cells have a `cell_type_final` label (including "unassigned")
- [ ] Proportion of "unassigned" cells is < 20% (flag if higher)
- [ ] Immune cell markers (CD68, CD3D, CD79A) are expressed predominantly in cells labeled as immune, not in resident cell clusters
- [ ] At least 80% of cells labeled as endothelial express PECAM1
- [ ] Annotation report HTML is generated
- [ ] Marker dot plot is generated

Cross-dataset:
- [ ] `annotation_comparison.tsv` is generated
- [ ] Cell type proportions are plausible (e.g., NP-only studies shouldn't have >50% AF-labeled cells)

## Human Checkpoint

### Review materials
- Annotation reports per dataset (UMAPs colored by cell type, dot plots)
- `annotation_comparison.tsv` — do the same cell types appear across studies?
- Continuous score distributions for the resident cell continuum

### Questions for the reviewer
1. Do the cell type labels make biological sense for each dataset?
2. For the chondrocyte/fibroblast continuum: are the discrete labels meaningful, or should we rely on continuous scores for downstream analysis?
3. Are there cell populations that appear in some studies but not others? Is this biology (different compartments, conditions) or technical artifact?
4. Should any gene signatures be revised based on what the data shows?
5. Is the annotation granularity appropriate — too coarse or too fine?
6. Do the original study annotations (if provided) agree with ours? Where they disagree, which is more credible?

### Potential plan revisions
- If the continuum is genuinely continuous within most datasets, the integration module should preserve this structure rather than trying to cluster it
- If certain cell states only appear in degenerated samples, this is a key finding that should inform the DE analysis design
- If annotation quality is poor for certain datasets (high unassigned rate, inconsistent markers), consider excluding them or downweighting in integration
- If notochordal cells are found only in neonatal samples, this confirms a developmental trajectory that the trajectory module should explore
