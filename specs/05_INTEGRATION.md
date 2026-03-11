# Module 05: Integration

## Objective

Integrate cells across studies into shared representations using tiered scANVI. Clustering (Module 06) and annotation (Module 07) are separate steps.

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
- **GSE205535 NNP:** Included in integration but excluded from DE analysis (Module 08) — 11yo spinal cord injury is a trauma confound

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
- `results/integration_report.html` — visualization of integration, clustering, and annotation results
- Updated `analysis_plan.md` with integration parameters and metrics

### Notebook: `notebooks/05_integration.ipynb`

Contains:
- Inclusion summary table (studies × objects)
- UMAP panels per object, colored by: study, condition, compartment
- Integration quality metrics (batch mixing, biological conservation)
- Study caveats table

**Manuscript mapping:** Supplementary Table S1: Study inclusion and caveats. Supplementary Figure S3: Integration benchmarking. Methods section on integration.

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

## Automated Validation

- [ ] All four integrated objects are saved (NP, AF, CEP, all_cells)
- [ ] Inclusion summary table is generated
- [ ] Study caveats table is generated
- [ ] No integration result collapses all cells into a single cluster at resolution 0.5 (blob check)
- [ ] No integration result has study identity perfectly predicting cluster identity (ARI < 1.0)
- [ ] Integration metrics are recorded per object

## Human Checkpoint

### Review materials
- Inclusion summary table (which studies/samples in which objects)
- Study caveats table
- UMAP per object colored by study, condition, compartment
- Integration quality metrics per object (iLISI, batch-ASW, condition-ASW)
- Continuum preservation check (mesenchymal tier): cluster count comparison before/after integration

### Questions for the reviewer
1. Does the integration look reasonable? Are studies mixing well on the UMAP?
2. Is there evidence of overcorrection (biological signal lost)?
3. Are the study caveats adequately documented for the supplement?
4. Should integration parameters be adjusted for any object/tier?

### Potential plan revisions
- If batch effects persist, consider adjusting scANVI parameters or adding covariates
- If the mesenchymal continuum collapses into a blob, reduce integration strength or use alternative methods
- If a specific study is an outlier, consider excluding it and re-integrating
