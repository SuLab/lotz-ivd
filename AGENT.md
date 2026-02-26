# AGENT.md — IVD Single-Cell Analysis Pipeline

## Overview

This is a human-gated agentic bioinformatics pipeline. The agent executes one task at a time following the module specs in `specs/`. The agent does NOT autonomously decide what to do — it follows the current active step in `analysis_plan.md`.

## Execution Rules

1. **Read the active step in `analysis_plan.md` before doing anything.**
2. **Read the relevant spec file before executing any step.**
3. **One task per execution.** Do not try to run multiple pipeline modules in one session.
4. **Run automated validation checks after completing a step.** If any check fails, attempt to fix it. If the fix requires changing parameters or approach, document the issue in `analysis_plan.md` under Known Issues and STOP — do not proceed to the next module.
5. **At human checkpoints, STOP.** Generate the review materials specified in the module spec and present them. Do not advance until the human has reviewed and approved.
6. **Record everything.** All commands, parameters, outputs, and decisions go in `analysis_plan.md`.
7. **Do not modify specs without human approval.** If a spec seems wrong or incomplete, flag it in `analysis_plan.md` under Known Issues.

## Computational Environment

- Python 3.10+
- Key packages: scanpy, scvi-tools, anndata, pandas, numpy, matplotlib, seaborn
- R (for DESeq2, edgeR, propeller): accessible via rpy2 or standalone R scripts
- Install packages as needed; record versions

## Running Scripts

- All scripts go in `scripts/`
- Name scripts by module: `01_dataset_discovery.py`, `03_preprocessing.py`, etc.
- Scripts should be self-contained and re-runnable
- Use `data/raw/` for input, `data/processed/` for intermediate output, `results/` for final output

## Git Usage

- Commit after completing each task with a descriptive message
- Do NOT commit large data files (h5ad, count matrices). Add them to .gitignore.
- DO commit: scripts, specs, metadata files, analysis_plan.md, AGENT.md

## When Things Break

- If a script fails: read the error, fix the script, rerun. Do not change the spec.
- If validation fails: document in analysis_plan.md, attempt to resolve. If unresolvable, STOP.
- If the data doesn't match expectations: document in analysis_plan.md, STOP for human review.
- If a module takes longer than expected: that's fine. Don't cut corners to save time.
