# Module 05: Integration

## Objective

Integrate cells across studies into shared representations using three parallel integration strategies, then compare results to select the best approach. Clustering (Module 06) and annotation (Module 07) are separate steps.

Four integrated objects are produced per workflow, each clustered and annotated independently:
1. **NP** — nucleus pulposus cells from studies with clearly separated NP tissue
2. **AF** — annulus fibrosus cells from studies with clearly separated AF tissue
3. **CEP** — cartilaginous endplate cells
4. **All-cells** — all IVD cells together, including studies where compartments were not separated

The all-cells object serves as a whole IVD atlas. IVD subcompartments (NP, AF, CEP) are related tissues with overlapping cell compositions — not distinct organs — so a combined analysis is biologically meaningful for understanding the full IVD cellular landscape and cross-compartment relationships.

### Three Integration Workflows

| | A: Seurat CCA (v5) | B: scANVI | C: STACAS | Notes |
|---|---|---|---|---|
| Language | R only | R + Python | R only | Prefer A or C to avoid Python dependency |
| Cell labels? | No | Yes (coarse) | Yes (coarse; optional) | Same coarse labels shared across B and C |
| Scalability | Standard CCA on full data | Handles large datasets natively | Anchor-based (similar to CCA) | |
| Best for | Standard multi-lab; primary analysis | Large atlases; probabilistic cell type transfer | R-native label-guided; strong bio-conservation | |

**Workflow A is the primary workflow.** Workflows B and C are run for comparison. The human checkpoint selects which workflow to carry forward.

## Rationale

Running three integration strategies guards against method-specific artifacts and provides the reviewer with an informed choice:

- **Workflow A (CCA)** is label-free — it finds shared correlation structure across datasets without requiring any prior cell type information. This avoids any risk of circular annotation influencing integration.
- **Workflows B (scANVI) and C (STACAS)** leverage the coarse anchor labels from Module 04 for semi-supervised correction, which can improve batch correction across platforms (10x, BD Rhapsody, Singleron) at the cost of depending on those labels.

Module 04 assigns coarse anchor labels (5 categories + Unknown) that are reliable across datasets but deliberately avoid fine-grained distinctions. Workflows B and C use those anchors. Fine cell types emerge from post-integration clustering and annotation (Modules 06–07), not from the integration step.

Within Workflows B and C, mesenchymal and non-mesenchymal cells are integrated separately (tiered integration) because they have fundamentally different transcriptomic profiles. Integrating them together would force the model to spend latent capacity on that major axis of variation rather than on subtler differences within each group. Workflow A tests both flat and tiered approaches.

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
| GSE165722 (Tu 2022) | BD Rhapsody platform (not 10x) | Different capture efficiency, gene detection | Platform-aware batch correction via study-level batch key |
| GSE205535 (Li 2022) | BD Rhapsody platform; published corrigenda | See above; potential data quality issues | Monitor for outlier behavior in integration |
| CNP0002664 (Han 2022) | Singleron Matrix platform (not 10x) | Different capture efficiency | Same as above |
| GSE242443 (Kuchynsky 2024) | Culture-expanded CEP cells | Culture alters cell states; may not reflect in vivo biology | Caveat in all CEP results; compare with non-expanded CEP from GSE160756 |
| GSE255768 (Shi 2024) | Degenerative endplate only; no healthy control | Cannot do healthy vs. degenerated comparison for this study alone | Healthy CEP baseline from GSE160756 |
| GSE230809 (Swahn 2024) | All-male donors; age-disease confounded | Cannot separate age from degeneration effects | Note in interpretation; sex-specific effects cannot be assessed |
| GSE205535 NNP sample | 11yo spinal cord injury, classified as "healthy" | Trauma confound | Excluded from DE comparisons |
| GSE189916 (Jiang 2022) | Whole IVD (compartments not separated) | Cannot assign cells to NP/AF/CEP | Included only in all-cells object |

## Inputs

- Processed Seurat objects / AnnData objects from Module 03 (per-dataset, QC'd, SCTransform-normalized)
- `coarse_label` and `cell_class` from Module 04 — **required for Workflows B and C; not required for Workflow A**
- `metadata/sample_metadata.tsv`

## Outputs

Per workflow (`{wf}` = `cca`, `scanvi`, `stacas`):
- `data/integrated/{wf}/NP.rds` (or `.h5ad`) — integrated NP cells
- `data/integrated/{wf}/AF.rds` — integrated AF cells
- `data/integrated/{wf}/CEP.rds` — integrated CEP cells
- `data/integrated/{wf}/all_cells.rds` — integrated all IVD cells
- `data/integrated/{wf}/integration_metrics.tsv` — quantitative integration assessment

CCA parameter sweep outputs (non-default parameters):
- `data/integrated/cca_hvg{N}_dims{D}/` — integrated objects for each parameter combination
- `results/integration_cca_hvg{N}_dims{D}/` — UMAPs and metrics for each combination

Shared outputs (generated once):
- `results/integration/inclusion_summary.tsv` — study × object inclusion table
- `results/integration/inclusion_summary.html` — formatted version for manuscript supplement
- `results/integration/study_caveats.tsv` — per-study caveats for supplement
- `results/integration/workflow_comparison.tsv` — side-by-side metrics across all three workflows
- `results/integration/workflow_comparison_report.html` — visual comparison of workflows
- Updated `analysis_plan.md` with integration parameters, metrics, and workflow selection

### Notebooks

#### `notebooks/05a_integration_cca.ipynb` — Workflow A results

- UMAP panels per object (NP, AF, CEP, all-cells), colored by: study, condition, compartment
- Integration quality metrics (iLISI, batch-ASW, condition-ASW)
- Flat vs tiered comparison (side-by-side UMAPs and metrics)
- Batch mixing assessment: per-study cell density on UMAP
- Continuum preservation check: mesenchymal cell spread in embedding space

#### `notebooks/05b_integration_scanvi.ipynb` — Workflow B results

- UMAP panels per object per tier (mesenchymal, non-mesenchymal), colored by: study, condition, coarse_label
- Integration quality metrics per tier
- scANVI training diagnostics (loss curves, convergence)
- Tier merge visualization: combined UMAP with tier of origin overlay
- Continuum preservation check

#### `notebooks/05c_integration_stacas.ipynb` — Workflow C results

- UMAP panels per object, colored by: study, condition, coarse_label
- Integration quality metrics
- Batch mixing assessment
- Continuum preservation check

#### `notebooks/05d_integration_comparison.ipynb` — Workflow comparison

- Side-by-side UMAP panels: same object across all three workflows
- Integration metrics comparison table and bar plots (iLISI, batch-ASW, cLISI, condition-ASW per workflow per object)
- Continuum preservation comparison across workflows
- Inclusion summary table (studies × objects) — shared across workflows
- Study caveats table
- Summary recommendation for human checkpoint

**Manuscript mapping:** Supplementary Table S1: Study inclusion and caveats. Supplementary Figure S3: Integration benchmarking and workflow comparison. Methods section on integration.

## Shared Integration Parameters

All three workflows share these parameters (from Shared Parameters):
- **Normalization:** SCTransform per sample, regressing out percent.mt
- **HVGs:** 3,000 selected across all datasets (`SelectIntegrationFeatures(nfeatures = 3000)`) — default; configurable via `--n-hvg`
- **Dimensionality:** 50 PCs for all PCA, neighbor graph, and integration steps — default; configurable via `--n-dims`

### CCA Parameter Sweep

The CCA script (`scripts/05a_integration_cca.R`) supports `--n-hvg` and `--n-dims` flags for parameter sweeps. Non-default parameter combinations write to separate output directories (e.g., `data/integrated/cca_hvg2000_dims30/`, `results/integration_cca_hvg2000_dims30/`) to allow side-by-side comparison with the baseline.

Planned parameter combinations:

| HVGs | Dims | Directory suffix | Rationale |
|------|------|-------------------|-----------|
| 3,000 | 50 | (default) | Baseline — original v5 integration |
| 2,000 | 50 | `cca_hvg2000_dims50` | Fewer HVGs to reduce batch-driven features |
| 3,000 | 30 | `cca_hvg3000_dims30` | Fewer dims to exclude noisy components |
| 2,000 | 30 | `cca_hvg2000_dims30` | Combined reduction |

The goal is to improve batch mixing (increase iLISI, move batch ASW toward 0) without losing biological signal (condition ASW should not become more negative).

---

## Workflow A (Primary): Seurat CCA Integration

**Framework:** R / Seurat v5. Use this as the main workflow for all results.

CCA (Canonical Correlation Analysis) finds shared correlation structure across datasets without requiring cell type labels. This makes it the most assumption-free approach — integration quality depends only on the data, not on the accuracy of prior annotations. This workflow follows the integration approach from the reference R code in `single_nuclei_r/Sample_QC_Integration.R`, updated to use Seurat v5's `IntegrateLayers` API.

**Seurat v5 vs v4:** The CCA algorithm is identical. Seurat v5 changes only the execution model: data is stored in per-study layers within a single merged object, and `IntegrateLayers(method = CCAIntegration)` replaces the v4 pipeline of `FindIntegrationAnchors` → `IntegrateData`. Standard (non-sketch) CCA is used for all objects regardless of size.

### Step 1 — Load and merge into a single Seurat v5 object

Load per-study h5ad files, apply compartment/sample filters, merge into one object where each study's counts are a separate layer.

### Step 2 — Normalization

`NormalizeData()` + `FindVariableFeatures(nfeatures = N_HVG)` + `ScaleData()` — log-normalization that preserves the per-study layer structure required by `IntegrateLayers`. The CCA algorithm itself is normalization-agnostic (it finds shared correlation structure); log-normalization is the standard Seurat v5 `IntegrateLayers` input. N_HVG defaults to 3,000; configurable via `--n-hvg`.

### Step 3 — PCA and CCA integration

For all objects (NP, AF, CEP, all_cells):
1. `RunPCA(npcs = N_DIMS)` on the merged object
2. `IntegrateLayers(method = CCAIntegration, orig.reduction = "pca", dims = 1:N_DIMS)`
3. `JoinLayers()` to recombine after integration

N_DIMS defaults to 50; configurable via `--n-dims`.

### Step 4 — Dimensionality reduction

1. `RunUMAP(reduction = "integrated.cca", dims = 1:N_DIMS)`
2. `FindNeighbors(reduction = "integrated.cca", dims = 1:N_DIMS)`

---

## Workflow B (Alternative): scANVI Semi-supervised Integration

**Framework:** R (preprocessing) + Python (integration via scvi-tools). Use when label-guided integration is desired or as a comparison to Workflow A.

### Rationale

scANVI is a semi-supervised deep generative model (VAE-based) that incorporates coarse cell type labels to guide integration, helping preserve biologically meaningful structure that purely unsupervised methods may over-correct. It operates on raw counts, so SCTransform normalization does not affect the integration itself — only HVG selection upstream.

### Step 1 — Upstream processing in R (same as Workflow A)

1. Apply QC filters: percent.mt < 5%, nCount_RNA 1,000–25,000
2. Run `SCTransform()` per sample (regress percent.mt) — for HVG selection only
3. Select 3,000 HVGs using `SelectIntegrationFeatures(nfeatures = 3000)`

### Step 2 — Coarse pre-annotation (from Module 04)

Labels are assigned per-cell using marker gene scoring, not cluster-based annotation. This avoids circular dependencies between clustering and annotation.

- Score each cell using `sc.tl.score_genes()` for:
  - Chondrocyte markers: COL2A1, ACAN, SOX9
  - Fibroblast markers: COL1A1, COL1A2, DCN, LUM
  - Non-mesenchymal markers: PTPRC, PECAM1, VWF, CDH5, RGS5, PDGFRB — checked against per-gene top-10th-percentile expression thresholds
- Assign coarse labels using a hierarchical first-match-wins rule:
  1. **Immune** — PTPRC in top 10th percentile, or ≥2 immune supporting markers (unless ACAN/SOX9 co-expressed — rescue rule)
  2. **Endothelial** — PECAM1/VWF/CDH5, same rescue rule
  3. **Pericyte_SMC** — RGS5 + PDGFRB co-expression, same rescue rule
  4. **Chondrocyte_like** — chondrocyte score > 2× fibroblast score and > 0
  5. **Fibroblast_like** — fibroblast score > 2× chondrocyte score and > 0
  6. **Unknown** — everything else
- Apply cluster-level smoothing: compute Leiden clusters (resolution 0.5); if >85% of cells in a cluster share the same label, assign that label to the entire cluster
- The "Unknown" category is intentionally generous (~20–30% of mesenchymal cells), capturing fibrochondrocytes and transitional cells; scANVI positions these by transcriptomic similarity
- Derive `cell_class`: mesenchymal (Chondrocyte_like, Fibroblast_like), non_mesenchymal (Immune, Endothelial, Pericyte_SMC), unknown (Unknown)

### Step 3 — Integration in Python (scvi-tools)

1. Load AnnData; subset to the 3,000 HVGs selected in R; store raw counts in `adata.layers["counts"]`
2. Train scVI first (unsupervised VAE) using batch = sample ID; then initialize scANVI from the scVI model using the coarse cell type labels
3. Extract the scANVI latent embedding; compute neighbors and UMAP in Python (scanpy), or export the embedding back to R

**scANVI operates on raw counts — do not pass log-normalized or SCT-corrected values as model input.**

**Parameters:**
- `batch_key='study'`
- `n_latent=20`
- scVI: `max_epochs=200`
- scANVI: `max_epochs=50`, `early_stopping=True`, `early_stopping_patience=5`
- `n_top_genes=3000`

### Tiered integration

**Tier A — Mesenchymal cells:**
- Anchor labels: Chondrocyte_like, Fibroblast_like, Unknown
- scANVI `unlabeled_category`: "Unknown"
- Subset cells with `cell_class == 'mesenchymal'` or `cell_class == 'unknown'` (include Unknown-class cells that may be fibrochondrocytes)
- Concatenate across datasets → re-identify HVGs (n_top_genes=3000) → train scVI → initialize scANVI → extract embedding → compute neighbors/UMAP

**Tier B — Non-mesenchymal cells:**
- Anchor labels: Immune, Endothelial, Pericyte_SMC
- scANVI `unlabeled_category`: "Unknown" (should be very few in this tier)
- Same workflow as Tier A. These populations have strong, discrete transcriptomic identities.

**Merging tiers:** After independent integration of each tier, merge back into a single object per compartment. Store tier-specific embeddings in `obsm` (e.g., `X_scanvi_mesenchymal`, `X_scanvi_non_mesenchymal`).

### Step 4 — Export to R

Import scANVI embedding into Seurat as a custom DimReduc object; proceed with `FindNeighbors()` to build the neighbor graph.

### Key notes

- "Unknown" cells are handled natively; scANVI supports partially labeled data
- Required Python packages: scvi-tools (≥ 1.0), scanpy, anndata
- The rescue rule prevents stressed IVD cells with upregulated immune markers from being misclassified as non-mesenchymal

---

## Workflow C (Alternative): STACAS Semi-supervised Integration

**Framework:** R only (Seurat v5 + STACAS package). Use when staying fully within R is preferred and coarse cell type priors are available.

### Rationale

STACAS is an anchor-based integration method (similar in spirit to Seurat CCA) that uses prior cell type labels to weight anchors, reducing overcorrection while preserving biological variability. Fully R-native, scalable to large datasets, and robust to imprecise or incomplete labels. Recent benchmarks show it outperforms both unsupervised methods and scANVI for within-species, multi-lab integration when even rough annotations are available.

### Step 1 — Upstream processing in R (same as Workflow A)

Apply QC filters and per-sample `SCTransform()` (regress percent.mt). Select 3,000 HVGs using `SelectIntegrationFeatures(nfeatures = 3000)`.

### Step 2 — Coarse pre-annotation in R (same as Workflow B, Step 2)

Use the same per-cell marker scoring approach described in Workflow B Step 2. Store resulting labels in a metadata column (e.g., `cell_type_coarse`).

### Step 3 — STACAS Integration

Install: `remotes::install_github("carmonalab/STACAS")`

Run semi-supervised integration:
```r
SampleIntegration(seurat_list, dims = 1:50, cell.labels = "cell_type_coarse")
```

STACAS returns an integrated Seurat object. Tiered integration is compatible: run within-study integration first, then across studies, applying cell labels at both steps.

### Step 4 — Dimensionality reduction

`RunPCA(npcs = 50)` → `RunUMAP()` → `FindNeighbors(dims = 1:50)`

### Key notes

- STACAS is natively compatible with Seurat v5; no Python required
- Cells can be labeled "Unknown" — STACAS still improves integration over the unsupervised baseline

---

## Workflow Comparison

After all three workflows complete, compare them to inform the human checkpoint decision.

### Metrics (computed per workflow, per object)

1. **Batch mixing:**
   - iLISI (integration LISI) — higher is better mixing
   - Batch-ASW — should be near 0

2. **Biological conservation:**
   - cLISI (cell-type LISI) — lower is better preservation
   - Condition-ASW — conditions should remain distinguishable
   - Condition classifier accuracy — can we still distinguish healthy from degenerated?

3. **Continuum preservation (mesenchymal cells only):**
   - Compare cluster count at resolution 0.5 before and after integration (large reduction suggests overcorrection)
   - Verify that the integrated object does not collapse into a single blob
   - Spread of mesenchymal cells in embedding space

### Comparison outputs

- `results/integration/workflow_comparison.tsv` — all metrics side-by-side
- `results/integration/workflow_comparison_report.html` — side-by-side UMAPs and metric plots
- Summary recommendation (but final decision is made at human checkpoint)

### Selection criteria

Workflow A (CCA) is used unless it clearly underperforms on batch mixing or biological conservation. If CCA performs adequately, prefer it because it is label-free and avoids any circular dependency between annotation and integration. If CCA underperforms, the comparison informs which alternative to select.

## Automated Validation

Per workflow:
- [ ] All four integrated objects are saved (NP, AF, CEP, all_cells)
- [ ] No integration result collapses all cells into a single cluster at resolution 0.5 (blob check)
- [ ] No integration result has study identity perfectly predicting cluster identity (ARI < 1.0)
- [ ] Integration metrics are recorded per object

Shared:
- [ ] Inclusion summary table is generated
- [ ] Study caveats table is generated
- [ ] Workflow comparison table and report are generated
- [ ] At least one workflow passes all integration quality checks per object

## Human Checkpoint

### Review materials
- Inclusion summary table (which studies/samples in which objects)
- Study caveats table
- Per-workflow UMAP per object colored by study, condition, compartment
- Side-by-side workflow comparison (UMAPs and metrics)
- Integration quality metrics per workflow per object (iLISI, batch-ASW, condition-ASW)
- Continuum preservation check (mesenchymal cells): cluster count comparison before/after integration

### Questions for the reviewer
1. Does the integration look reasonable? Are studies mixing well on the UMAP?
2. Is there evidence of overcorrection (biological signal lost) in any workflow?
3. **Which workflow should be carried forward for downstream analysis (Modules 06–12)?**
4. For the selected workflow: should integration parameters be adjusted for any object/tier?
5. Are the study caveats adequately documented for the supplement?

### Potential plan revisions
- If Workflow A (CCA) performs adequately, proceed with it as the primary result — it is label-free and avoids circular annotation risk
- If CCA underperforms on batch mixing (especially across platforms), select the best-performing alternative (B or C)
- If the mesenchymal continuum collapses into a blob in any workflow, reduce integration strength or exclude that workflow
- If a specific study is an outlier in all workflows, consider excluding it and re-integrating
- If batch effects persist across all workflows, consider adding covariates or restricting to within-study comparisons
