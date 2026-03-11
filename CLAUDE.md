# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

IVD Single-Cell Atlas: a human-gated agentic bioinformatics pipeline that analyzes publicly available scRNA-seq datasets of human intervertebral disc (IVD) tissue. The goal is to identify cell types/states and how they change with aging and degeneration.

The pipeline is driven by an agent loop (`while :; do cat PROMPT.md | claude; done`). Each iteration, the agent reads the current state from `analysis_plan.md`, executes one task per the relevant module spec in `specs/`, runs validation, and updates the plan. The loop halts at human checkpoints.

## Environment Setup

```bash
source .venv/bin/activate
pip install -r requirements.txt
```

Requires `python3-dev` system package (for annoy/scib-metrics). Use `python3` not `python` to invoke scripts. R with DESeq2/propeller needed for Module 08+.

## Running the Pipeline

1. Read `analysis_plan.md` for the current active step
2. Read the relevant spec in `specs/` (e.g., `specs/03_PREPROCESSING.md`)
3. Read `AGENT.md` for execution rules
4. Run the script: `python3 scripts/03_preprocessing.py`
5. Generate/update the corresponding notebook in `notebooks/`
6. Update `analysis_plan.md` with outcomes
7. Commit scripts and metadata to git (not data files)

## Architecture

**12-module pipeline** (specs in `specs/01_*.md` through `specs/12_*.md`):
1. Dataset Discovery → 2. Metadata Harmonization → 3. Preprocessing → 4. Coarse Annotation → 5. Integration → 6. Clustering → 7. Post-Integration Annotation → 8. Differential Analysis → 9. Interpretation → 10. Trajectory → 11. Cell-Cell Communication → 12. Reporting

**Scripts vs Notebooks split:** Scripts (`scripts/`) do heavy compute headlessly. Notebooks (`notebooks/`) load saved results from `data/` and `results/` for visualization — they are independent of scripts and serve as checkpoint review artifacts and draft manuscript figures.

**Data flow:** `data/raw/` → `data/processed/` (per-dataset .h5ad) → `data/integrated/` → `results/`. Canonical format is AnnData (.h5ad).

## Key Execution Rules

- **One task per session.** Do not run multiple modules in one invocation.
- **Stop at human checkpoints.** Prepare review materials and stop — do not advance.
- **Do not modify spec files** without human approval. Flag issues in `analysis_plan.md`.
- **Run all automated validation checks** listed in the spec after completing a step.
- **Record everything** in `analysis_plan.md`: parameters, outputs, decisions, issues.
- `analysis_plan.md` is the source of truth for pipeline status and all decisions.

## Git Conventions

- **Commit:** scripts, specs, metadata files, `analysis_plan.md`, `AGENT.md`, notebooks
- **Do NOT commit:** data files (.h5ad, .h5, .mtx, count matrices), `results/`, `.ipynb_checkpoints/`
- Data integrity tracked via `metadata/file_checksums.json`

## Key Scientific Context

- IVD resident cells (chondrocyte-like, fibroblast-like) exist on a **continuum**, not discrete types. Standard batch correction can erase this variation — the pipeline uses per-dataset processing first, then tiered integration.
- **Tiered integration:** non-resident cells (immune, endothelial) use standard methods; resident cells use conservative approaches.
- **Pseudobulk DE** (DESeq2), not single-cell DE, to avoid treating cells as independent observations.
- 12 datasets, ~423K cells, 78 samples, 57 donors across NP/AF/CEP compartments.

## Current Dataset List

12 studies: GSE160756, GSE165722, GSE189916, GSE199866, GSE205535, CNP0002664, GSE233666, GSE244889, GSE251686, GSE255768, GSE230809, GSE242443. Full details in `AGENT.md` and `metadata/dataset_registry.tsv`.

## Key Known Issues

- GSE205535 has published corrigenda
- 3 non-10x platform datasets (BD Rhapsody, Singleron) — need platform-aware batch correction
- GSE242443 CEP cells are culture-expanded (included by decision, requires caveats)
- GSE230809: all-male donors with confounded age-disease effects
- Condition mappings must be revisited before Module 08 (DE analysis)
