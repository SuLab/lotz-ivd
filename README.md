# IVD Single-Cell Atlas

**Goal:** Identify the cell types and cell states present in the human intervertebral disc (IVD) and determine how these change with aging and degeneration.

## Approach

This project uses a human-gated agentic pipeline to analyze publicly available single-cell RNA-seq datasets of human IVD tissue. An AI agent (Claude) executes well-defined computational steps, while a human PI reviews results and makes scientific decisions at defined checkpoints.

The pipeline is driven by a loop:

```bash
while :; do cat PROMPT.md | claude; done
```

Each iteration, the agent reads the current state from `analysis_plan.md`, executes the next task as defined in the module specs, runs validation, and updates the plan. The loop halts at human checkpoints — the agent prepares review materials and stops until the human advances the plan.

## Pipeline Modules

| # | Module | Description |
|---|--------|-------------|
| 01 | Dataset Discovery | Systematic search for all human IVD scRNA-seq datasets |
| 02 | Metadata Harmonization | Standardize condition labels, demographics, covariates |
| 03 | Preprocessing | Per-dataset QC, normalization, clustering |
| 04 | Annotation | Cell type labeling using markers, references, and automated tools |
| 05 | Integration | Cross-study integration with tiered strategy for resident vs. non-resident cells |
| 06 | Differential Analysis | Cell composition changes and pseudobulk DE between conditions |
| 07 | Biological Interpretation | Pathway enrichment, GRNs, pain-associated gene analysis |
| 08 | Trajectory & Dynamics | Pseudotime, RNA velocity for the cell state continuum |
| 09 | Cell-Cell Communication | Ligand-receptor interactions between IVD cell populations |
| 10 | Reporting | Final report, figures, reproducibility documentation |

## Key Files

- `PROMPT.md` — Fed to the agent on each loop iteration
- `AGENT.md` — Execution rules and environment instructions
- `analysis_plan.md` — Living document tracking progress, decisions, and revisions
- `specs/` — Module specifications defining inputs, outputs, methods, validation, and checkpoints

## Compute Requirements

Most modules run on a standard workstation (32GB RAM). Two modules benefit from HPC:

- **Module 05 (Integration):** scVI/scANVI training is significantly faster with GPU. ~200k cells across all studies.
- **Module 07 (SCENIC):** Gene regulatory network inference is RAM-intensive (64GB+).

## Key Design Decisions

**Human-gated, not fully autonomous.** The agent executes computational steps but stops at decision points for human review. The automated validation checks are regression safeguards, not proof of correctness.

**Per-dataset first, then integrate.** Each dataset is preprocessed and annotated independently before cross-study integration. This avoids the known problem where batch correction erases the subtle cell state variation in the chondrocyte/fibroblast continuum.

**Tiered integration.** Non-resident cells (immune, endothelial) integrate easily and are handled with standard methods. Resident IVD cells require conservative integration or alternative approaches (label transfer, metacells) to preserve the biological continuum.

**The plan is revisable.** Every human checkpoint includes an evaluation of whether the downstream plan still makes sense given what's been learned. The analysis may loop back to earlier steps with revised parameters.

**Pseudobulk DE, not single-cell DE.** Differential expression uses pseudobulk aggregation (DESeq2/edgeR) to avoid inflated statistics from treating cells as independent observations.

## Directory Structure

```
ivd-analysis/
├── PROMPT.md               # Agent loop prompt
├── AGENT.md                # Agent execution rules
├── README.md               # This file
├── analysis_plan.md        # Living plan document
├── specs/                  # Module specifications
│   ├── 00_PROJECT.md
│   ├── 01_DATASET_DISCOVERY.md
│   ├── 02_METADATA.md
│   ├── 03_PREPROCESSING.md
│   ├── 04_ANNOTATION.md
│   ├── 05_INTEGRATION.md
│   ├── 06_DIFFERENTIAL.md
│   ├── 07_INTERPRETATION.md
│   ├── 08_TRAJECTORY.md
│   ├── 09_COMMUNICATION.md
│   └── 10_REPORTING.md
├── data/
│   ├── raw/                # Downloaded datasets
│   ├── processed/          # Per-dataset h5ad files
│   └── integrated/         # Cross-study integrated objects
├── metadata/               # Dataset registry, sample metadata
├── results/                # All analysis outputs
├── scripts/                # Analysis scripts
└── notebooks/              # Jupyter notebooks for checkpoint review
```

## Citation

If this analysis contributes to a publication, cite:
- The original publications for each included dataset
- The tools used (scanpy, scvi-tools, DESeq2, etc.)
- This pipeline methodology as appropriate
