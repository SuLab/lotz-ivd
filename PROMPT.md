# IVD Single-Cell Analysis — Loop Prompt

## Instructions

You are executing one step of a bioinformatics analysis pipeline for human intervertebral disc (IVD) single-cell RNA-seq data.

**Before doing anything:**

1. Read `analysis_plan.md` to determine the current active step.
2. Read the relevant spec file in `specs/` for that step.
3. Read `AGENT.md` for execution rules.

**Then:**

4. Execute the active step. One task only — do not skip ahead to future modules.
5. After completing the task, run all automated validation checks listed in the spec.
6. If all checks pass, update `analysis_plan.md`:
   - Move the completed step to the Completed Steps table with today's date and outcome.
   - Advance the Active Step to the next item in the Pending Steps list.
   - Log any issues or observations under Known Issues or Deferred Questions.
7. If a check fails, attempt to fix it. If you cannot fix it, document the failure in `analysis_plan.md` under Known Issues and STOP.
8. If the next step is a **human checkpoint**, prepare the review materials specified in the spec, summarize what you did and what needs review, then STOP. Do not proceed past a human checkpoint.

**Rules:**

- Do not modify spec files without human approval.
- Do not skip validation checks.
- Do not advance past human checkpoints.
- Do not run multiple modules in one session.
- Commit scripts and metadata to git after each completed task. Do not commit large data files.
- Record all parameter choices and command outputs in `analysis_plan.md`.
- If you discover something unexpected (a dataset issue, a biological anomaly, a tool limitation), document it in `analysis_plan.md` even if it doesn't block the current task.
- If the active step says "WAITING FOR HUMAN REVIEW", do nothing. Report the current status and stop.
