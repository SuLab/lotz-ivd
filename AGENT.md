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
- R (for DESeq2, propeller): accessible via rpy2 or standalone R scripts
- Install packages as needed; record versions

## Running Scripts

- All scripts go in `scripts/`
- Name scripts by module: `01_dataset_discovery.py`, `03_preprocessing.py`, etc.
- Scripts should be self-contained and re-runnable
- Use `data/raw/` for input, `data/processed/` for intermediate output, `results/` for final output

## Notebooks

- All visualization/interpretation notebooks go in `notebooks/`
- Name notebooks by module: `01_datasets.ipynb`, `03_qc.ipynb`, etc.
- **Each module produces both a script (compute) and a notebook (visualization).**
- Notebooks load saved output files from `results/` and `data/` — they must NOT depend on in-memory objects from scripts.
- Notebooks should be executable independently: a reviewer should be able to run the notebook without re-running the compute scripts.
- Notebooks serve as the review artifact at human checkpoints.
- Notebooks are the draft manuscript figures. Each maps to specific figures/tables (see spec files for mapping).

## Git Usage

- Commit after completing each task with a descriptive message
- Do NOT commit large data files (h5ad, count matrices). Add them to .gitignore.
- DO commit: scripts, specs, metadata files, analysis_plan.md, AGENT.md

## Final Dataset List (established Module 01, approved at checkpoint 2026-02-26)

12 datasets, ~423K reported cells, ~81 samples:

| Accession | Author | Compartment | Samples | Platform | Format |
|-----------|--------|-------------|---------|----------|--------|
| GSE160756 | Gan 2021 | NP, AF, CEP | 7 | 10x | .loom.gz |
| GSE165722 | Tu 2022 | NP | 8 | BD Rhapsody | counts.tsv.gz + cellname.txt.gz |
| GSE189916 | Jiang 2022 | Whole IVD | 6 | 10x | MTX triplet |
| GSE199866 | Cherif 2022 | NP, iAF | 4 | 10x | .h5 |
| GSE205535 | Li 2022 | NP | 2 | BD Rhapsody | MTX triplet |
| CNP0002664 | Han 2022 | NP | 6 | Singleron | matrix.tsv.gz |
| GSE233666 | Guo 2023 | NP | 4 | 10x | MTX triplet |
| GSE244889 | Chen 2024 | NP | 7 | 10x | MTX triplet |
| GSE251686 | Jia 2024 | NP | 6 | 10x | nested tar.gz |
| GSE255768 | Shi 2024 | CEP | 2 | 10x | MTX triplet |
| GSE230809 | Swahn 2024 | NP, AF | 24 | 10x | MTX triplet |
| GSE242443 | Kuchynsky 2024 | CEP | 2 | 10x | MTX triplet |

Notes:
- GSE242443 CEP cells are culture-expanded (included by human decision for coverage)
- Zhou 2023 (embryonic IVD) deferred to Module 10 trajectory analysis
- PRJCA014236, PRJCA007656 excluded at checkpoint (NP well-covered, NGDC access not obtained)
- Use `python3` not `python` to invoke scripts

## When Things Break

- If a script fails: read the error, fix the script, rerun. Do not change the spec.
- If validation fails: document in analysis_plan.md, attempt to resolve. If unresolvable, STOP.
- If the data doesn't match expectations: document in analysis_plan.md, STOP for human review.
- If a module takes longer than expected: that's fine. Don't cut corners to save time.
