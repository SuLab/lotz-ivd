# IVD Single-Cell Analysis — Loop Prompt

## Instructions

You are executing one step of a bioinformatics analysis pipeline for human
intervertebral disc (IVD) single-cell RNA-seq data.

**Before doing anything:**

1. Read `analysis_plan.md` to determine the current active step.
2. If the Active Step says "WAITING FOR HUMAN REVIEW", report the current
   status and STOP. Do not do any work.
3. Read the relevant spec file in `specs/` for that step.
4. Read `AGENT.md` for execution rules.

**Then execute the active step (one task only):**

5. Run the compute script for the active module.
6. Run all automated validation checks listed in the spec. If any check
   fails, attempt to fix it. If you cannot fix it, document the failure in
   `analysis_plan.md` under Known Issues and STOP.

**After compute succeeds, complete ALL deliverables before advancing:**

7. Generate or update the corresponding notebook in `notebooks/`. The
   notebook MUST:
   - Load results from `data/` and `results/` (not from in-memory objects)
   - Execute cleanly with `jupyter nbconvert --execute` (zero errors)
   - Contain no stale or placeholder text — all markdown must describe
     what was actually computed
   - Cover every visualization listed in the spec's notebook section
   - Have all markdown cells (titles, descriptions, status summaries)
     reflect the current pipeline version and integration method from
     `analysis_plan.md`. When updating a notebook from a previous run,
     review and rewrite EVERY markdown cell — not just re-execute code.
8. Verify the notebook by running:
   `jupyter nbconvert --to notebook --execute --inplace notebooks/XX_name.ipynb`
   If it fails, fix it before proceeding.
   Then verify no stale version references remain:
   `python3 -c "import json,sys; nb=json.load(open(sys.argv[1])); [sys.exit(f'STALE: cell {i}: {line.strip()}') for i,c in enumerate(nb['cells']) if c['cell_type']=='markdown' for line in c['source'] if any(old in line for old in ['(v4)', 'scANVI-based'])]" notebooks/XX_name.ipynb`
   Update the stale-string list in this check to include previous version
   identifiers whenever the pipeline version changes.
9. Update `analysis_plan.md`:
   - Move the completed step to the Completed Steps table with today's
     date, outcome, and key parameters/decisions.
   - Set the Active Step to the module's **human checkpoint** (not the
     next module). Use this exact format:
     `**Module XX: Human checkpoint** — WAITING FOR HUMAN REVIEW`
   - Log any issues under Known Issues.
10. Commit all scripts, notebooks, and metadata to git.
11. STOP. Do not proceed to the next module. The human must review and
    approve before the next module begins.

**Rules:**

- Do not modify spec files without human approval.
- Do not skip validation checks.
- Do not advance past human checkpoints.
- Do not run multiple modules in one session.
- Do not commit notebooks that have not been re-executed successfully.
- Commit scripts and metadata to git after each completed task. Do not
  commit large data files.
- Record all parameter choices and command outputs in `analysis_plan.md`.
- If you discover something unexpected, document it in `analysis_plan.md`.
- If the active step says "WAITING FOR HUMAN REVIEW", do nothing. Report
  the current status and stop.
