# IVD Single-Cell Atlas Project Specification

## Project Goal

Identify the cell types and cell states present in the human intervertebral disc (IVD) and determine how these change with aging and degeneration.

## Biological Context

The IVD consists of three main compartments: nucleus pulposus (NP), annulus fibrosus (AF), and cartilaginous endplate (CEP). The resident cell populations are primarily chondrocyte-like and fibroblast-like cells that exist on a continuum rather than as sharply distinct types. Immune cells and endothelial cells are also present and are more clearly separable. This continuum is a central challenge: standard clustering approaches tend to either over-split noise or under-resolve real cell states, and batch correction during cross-study integration can erase the subtle variation that distinguishes these states.

## Strategy

Use a human-gated agentic pipeline. The agent executes well-defined computational steps with automated validation. At defined decision points, a human reviews results and decides whether to advance, adjust, or revise the plan. The agent is the executor, not the PI.

## Key Variables

The analysis must track and account for the following sources of variation:

- **Biological:** age (neonatal, young adult, aged), disease state (healthy, mild degeneration, severe degeneration), IVD compartment (NP, AF, CEP, whole IVD)
- **Technical:** study of origin, sequencing platform, tissue vs. dissociated cells, donor, sample processing protocol

Age and degeneration are potentially separable variables and should not be conflated without evidence.

## Pipeline Modules

Each module is defined in its own spec file. Modules execute in order but the plan is revisable at each human checkpoint.

| Module | Spec File | Automatable? | Human Checkpoint? |
|--------|-----------|-------------|-------------------|
| Dataset discovery & acquisition | 01_DATASET_DISCOVERY.md | Mostly | Yes — approve final dataset list |
| Metadata harmonization | 02_METADATA.md | Mostly | Yes — approve condition mappings |
| Per-dataset preprocessing | 03_PREPROCESSING.md | Yes | Yes — review QC reports |
| Coarse cell classification | 04_ANNOTATION.md | Yes | Yes — approve mesenchymal vs. non-mesenchymal split |
| Integration, clustering & annotation | 05_INTEGRATION.md | Partially | Yes — evaluate integration, approve cell type atlas |
| Differential analysis | 06_DIFFERENTIAL.md | Mostly | Yes — review DE results |
| Biological interpretation | 07_INTERPRETATION.md | Partially | Yes — evaluate biological claims |
| Trajectory & dynamics | 08_TRAJECTORY.md | Partially | Yes — evaluate trajectory validity |
| Cell-cell communication | 09_COMMUNICATION.md | Mostly | Yes — review interaction results |
| Reporting & reproducibility | 10_REPORTING.md | Yes | Yes — final review |

## Decision Framework

### Automated Validation (regression safeguards)

These are necessary but not sufficient conditions. They block progression if they fail. They do NOT confirm that a step is scientifically correct or complete.

Examples:
- Expected output files exist and are non-empty
- Cell counts are within plausible range after QC
- Known marker genes appear in expected clusters (e.g., immune markers in immune clusters)
- No runtime errors or warnings that indicate data corruption
- Integration metrics (kBET, LISI) are within defined thresholds

### Human Decision Points

At each checkpoint, the human evaluates:

1. **Are the results scientifically sound?** Not just "did the code run" but "does this make biological sense?"
2. **Does the downstream plan still make sense?** Given what we've learned, should we adjust the approach for subsequent modules?
3. **Should we revisit a previous step?** New information may reveal that an earlier step, even one that passed automated checks, should be rerun with different parameters or additional covariates.

### Plan Revision Triggers (non-exhaustive)

- A batch effect is discovered that wasn't accounted for in preprocessing
- Cell type annotations don't match expected biology (e.g., no immune cells found, or markers are inconsistent)
- Integration produces the "blob" problem for resident IVD cells — triggers a strategy change
- A new dataset is discovered that substantially changes the sample composition
- DE analysis reveals that a confounding variable (e.g., tissue vs. cells) dominates over the biological signal
- Trajectory analysis contradicts assumptions about cell state relationships

## Plan Document (analysis_plan.md)

A living document maintained alongside the analysis. Structure:

```
# IVD Analysis Plan

## Current Status
[Which module is active, what was the last decision point outcome]

## Completed Steps
[Step, date, outcome, any notes]

## Active Step
[What is being executed now]

## Pending Steps
[Ordered list of what comes next, with any known parameter choices]

## Revisions Log
[Record of any plan changes and the reasoning behind them]

## Known Issues
[Problems discovered that need resolution, even if not blocking current work]

## Deferred Questions
[Scientific questions that emerged during analysis but are not part of the current scope]
```

## Computational Environment

- **Language:** Python (scanpy/scverse ecosystem as primary, R/Bioconductor for specific tools like DESeq2, edgeR)
- **Key packages:** scanpy, scvi-tools, scanorama, bbknn, harmonypy, decoupler, liana, cellchat, scenic/pyscenic, scvelo, squidpy, scikit-misc
- **Data format:** AnnData (.h5ad) as the canonical format; convert from other formats as needed
- **Version control:** All scripts and specs under git. Data files tracked via manifest, not committed.
- **Reproducibility:** Random seeds set for all stochastic operations. Package versions pinned. All parameters recorded in plan document.

## File Organization

```
ivd-analysis/
├── specs/                  # Module specifications (this directory)
├── data/
│   ├── raw/                # Downloaded raw data, organized by study
│   ├── processed/          # Per-dataset processed h5ad files
│   └── integrated/         # Cross-study integrated objects
├── metadata/
│   ├── dataset_registry.tsv    # Master list of all datasets
│   └── sample_metadata.tsv     # Harmonized sample-level metadata
├── results/
│   ├── qc_reports/         # Per-dataset QC summaries and plots
│   ├── annotations/        # Cell type annotation results
│   ├── differential/       # DE and composition results
│   ├── trajectories/       # Trajectory analysis outputs
│   ├── communication/      # Cell-cell communication results
│   └── figures/            # Publication-ready figures
├── scripts/                # Analysis scripts (heavy compute, run by agent or on HPC)
├── notebooks/              # Jupyter notebooks (visualization, interpretation, manuscript figures)
│   ├── 01_datasets.ipynb
│   ├── 02_metadata.ipynb
│   ├── 03_qc.ipynb
│   ├── 04_annotation.ipynb
│   ├── 05_integration.ipynb
│   ├── 06_differential.ipynb
│   ├── 07_interpretation.ipynb
│   ├── 08_trajectory.ipynb
│   └── 09_communication.ipynb
├── analysis_plan.md        # Living plan document
└── AGENT.md                # Instructions for the agent on how to run the pipeline
```

## Scripts vs. Notebooks

The pipeline uses a deliberate split between scripts and notebooks:

**Scripts** (`scripts/`) handle compute: downloading data, running QC, training models, executing DE, running SCENIC. These are what the agent executes in the loop. They write output files to `data/` and `results/`. They can run headlessly on HPC. They should be fast to rerun and self-contained.

**Notebooks** (`notebooks/`) handle visualization and interpretation: they load outputs from scripts and produce figures, summary tables, and narrative text. These are the review artifacts at each human checkpoint. They are also the draft manuscript figures — each notebook maps to a section of the paper.

The agent should produce both: a script that does the compute, and a notebook that visualizes the results. The notebook should be executable independently of the script (it loads saved output files, not in-memory objects from the script).

## Manuscript Figure Mapping

Each notebook corresponds to one or more manuscript figures and tables. This mapping is preliminary and will be refined as the analysis progresses.

| Notebook | Manuscript Section | Likely Figures/Tables |
|----------|-------------------|----------------------|
| 01_datasets.ipynb | Methods: Data sources | Table 1: Dataset characteristics |
| 02_metadata.ipynb | Methods: Study design | Table 1 (continued): Sample metadata summary |
| 03_qc.ipynb | Supplementary | Fig S1: QC metrics per dataset |
| 04_classification.ipynb | Supplementary | Fig S2: Mesenchymal vs. non-mesenchymal classification QC |
| 05_integration.ipynb | Results: IVD cell atlas; Methods: Integration | Fig 1: UMAP atlas, marker dot plots; Fig S3: Integration benchmark, resolution optimization |
| 06_differential.ipynb | Results: Disease-associated changes | Fig 2: Composition changes; Fig 3: DE volcano/heatmaps; Table 2: Top DE genes |
| 07_interpretation.ipynb | Results: Pathways & regulation; Discussion: Pain | Fig 4: Pathway enrichment; Fig 5: Pain-associated genes; Fig S4: GRN regulons |
| 08_trajectory.ipynb | Results: Cell state transitions | Fig 6: Pseudotime trajectory, gene dynamics |
| 09_communication.ipynb | Results: Intercellular signaling | Fig 7: Cell-cell communication changes |
