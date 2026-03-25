# Lessons Learned: Human-Gated Agentic Bioinformatics

> A retrospective on the IVD Single-Cell Atlas pipeline (v1–v5, Feb 26 – Mar 25, 2026).
> Written to inform the design of a reusable, SME-driven analysis framework.

---

## 1. Project Summary

Over 28 days, an AI agent (Claude Opus 4.6) executed a 12-module scRNA-seq
meta-analysis pipeline across 5 pipeline versions, processing ~423K cells
from 12 public datasets of human intervertebral disc tissue. Andrew Su
designed the original system architecture — the agent loop, spec format,
checkpoint protocol, and module structure. Ben Good served as the primary
human operator — running the agent loop, monitoring execution, diagnosing
stalls and crashes, driving compute environment changes, and directing the
agent through problems it could not resolve on its own. Martin Lotz
provided subject matter expertise on IVD biology, guiding scientific
decisions at checkpoints. Lisa Janssen and Hannah Swahn contributed
bioinformatics domain expertise through spec revisions that triggered the
v2 and v5 restructures respectively. The agent authored or co-authored 85
of 96 commits; the remaining 11 human commits drove every major
architectural change.

The pipeline produced cell type annotations, differential expression results,
pathway enrichments, trajectory analyses, cell-cell communication networks,
and a pain biology synthesis. But more instructively, it produced five
versions — each triggered by a domain expert identifying something the
agent could not: a flawed annotation strategy, a misclassified cell
population, an anchoring bias in result interpretation, or a missing
integration method.

---

## 2. What Worked

### 2.1 The agent-loop architecture

The `while :; do cat PROMPT.md | claude; done` loop with `analysis_plan.md`
as shared state was remarkably effective for execution. The agent could:

- Read the current pipeline state, identify the active step, execute it,
  and update the plan — all in a single invocation
- Write or rewrite 1,000+ line scripts from spec descriptions
- Debug OOM crashes, API version mismatches, and format incompatibilities
- Generate notebooks, run validation, and commit — a complete
  build-test-document cycle per iteration

### 2.2 Specs as the contract between humans and agent

The `specs/` directory served as the authoritative interface between domain
expertise and execution. When Lisa Janssen rewrote the annotation spec
(triggering v2) or Hannah Swahn restructured the integration spec
(triggering v5), the agent adapted without confusion. The spec format —
Objective, Inputs/Outputs, Method, Validation, Human Checkpoint — gave the
agent enough structure to execute while leaving scientific decisions to the
checkpoint.

### 2.3 Compartment-specific analysis

The v2 decision to create four compartment objects (NP, AF, CEP, all_cells)
instead of pooling all cells was a domain-expert insight that the agent
would never have initiated. It prevented cell-type confusion across tissue
types and enabled compartment-specific trajectories. This architectural
pattern — separate objects for biologically distinct populations, merged
only for cross-compartment comparison — should be a default in future
frameworks.

### 2.4 Multi-workflow comparison (v5)

Running three integration methods (CCA, scANVI, STACAS) in parallel and
comparing them at a human checkpoint was the most mature design pattern in
the project. It removed the "pick one method and hope" problem. The
comparison metrics (iLISI, batch_ASW, condition_ASW) gave the human
reviewer objective criteria for selection rather than requiring deep
integration expertise.

---

## 3. What Went Wrong

### 3.1 Checkpoint strategy must be an explicit design choice

During v1, the human operator intentionally pushed the agent through the
entire pipeline — all modules, including manuscript generation — without
SME review at any checkpoint. This was a deliberate decision: get a
complete end-to-end result as fast as possible so that the SMEs could see
the full pipeline output, including a draft manuscript, and have concrete
artifacts to react to. It worked — the v1 results and manuscript
stimulated discussion, demonstrated the approach's potential, and gave the
domain experts something specific to critique rather than an abstract plan.

However, this created a problem for subsequent versions: the pipeline had
been designed with checkpoints as mandatory gates, but had been operated
with checkpoints as optional. When the project moved to a mode where
checkpoints needed to be enforced (v2+), there was no mechanism to
guarantee it. The agent loop would happily continue past a checkpoint if
the human didn't intervene.

**Fix applied:** `run_pipeline.sh` was created as a shell-level gate that
greps `analysis_plan.md` for "WAITING FOR HUMAN REVIEW" and refuses to
start a new agent iteration if found. PROMPT.md was also hardened with
redundant stop directives.

**Lesson:** The framework should support multiple checkpoint modes as a
first-class concept — not just "checkpoint on" or "checkpoint off," but
explicitly selectable strategies:

- **Fast pass:** Run all modules without stopping. Useful for a first cut
  or when rerunning after minor parameter changes.
- **Gated:** Stop at every checkpoint. Required for new analyses or after
  major structural changes.
- **Selective:** Stop only at designated high-leverage checkpoints (e.g.,
  annotation, integration selection) and auto-advance past low-risk ones
  (e.g., QC, trajectory).

The v1 experience shows that all three modes are valuable at different
stages of a project. The framework should make the choice explicit and
mechanically enforced, rather than relying on the operator to manually
manage the loop.

### 3.2 The agent fabricated gene-level claims (v1)

During manuscript generation, the agent produced 14 factual errors about
specific genes and pathways. Examples:

- COMP was claimed as downregulated but was actually significantly
  *upregulated* in 4/6 cell types
- ECM pathways were claimed as "consistently downregulated in all 6 cell
  types" but were actually upregulated in two
- IL6 and CXCL8 were listed as "key upregulated genes" but neither reached
  statistical significance

These were not formatting errors or misinterpretations — the agent generated
specific, confident, wrong claims about computed results. A systematic
fact-check against the actual data tables caught all 14.

**Root cause:** The agent was asked to write a narrative synthesis from
summary statistics and pathway names. It generated plausible-sounding
biological narratives that were not grounded in the actual numbers. This is
a known failure mode of language models, but it is especially dangerous in
scientific reporting where specific claims carry weight.

**Fix applied:** `CLAUDE.md` was updated with explicit interpretation
guidelines: "default to conservative, defensible claims," "do not
editorialize," "bold guesses only when asked." A "manuscript_fact_check.md"
was added documenting each error. The pipeline spec (Module 12) was updated
to require systematic validation of all gene-level claims.

**Lesson:** The agent should never write scientific prose that makes specific
quantitative claims without those claims being mechanically verified against
the source data. In a future framework, any gene name, fold change, or
p-value in generated text should be programmatically checked against the
result tables before the text is finalized. Narrative generation and fact
verification must be separate steps.

### 3.3 Annotation was the single largest source of downstream variation

Across versions, changing the annotation strategy alone altered:
- DE gene counts by 22%
- Significant TF associations by 98% (from 290 to 5 after fixing
  misrouted cells in v3)
- CCC interaction direction (reversed between v1, v2, and v3)

The v2 binary classifier misrouted ~17K stressed NP cells to the
non-mesenchymal tier because stress markers (NAMPT, SOD2, CXCL8) confused
the scoring. These cells were then integrated with immune cells, producing
biologically meaningless DE results that propagated through every downstream
module.

**Lesson:** Cell type annotation is the highest-leverage decision in a
scRNA-seq pipeline. In v4, annotation was split into its own module (07)
with the spec explicitly noting it as "the most critical decision point in
the entire pipeline." Future frameworks should treat annotation as a
first-class human checkpoint — not something the agent resolves
autonomously.

### 3.4 Results that changed direction between versions were treated as findings

CCC interaction counts (healthy vs. degenerated) reversed direction in
v1, v2, v3, and v5. Trajectory correlations changed sign between versions.
These were eventually recognized as "method-sensitive" results and flagged
with caveats, but in early versions they were reported as findings.

**Lesson:** Cross-version instability is a useful signal, but not a
definitive one. Results that shift direction across versions *may* indicate
method sensitivity — or they may reflect genuine improvements as upstream
errors are corrected (e.g., the v2→v3 annotation fix legitimately changed
downstream results). The framework should track and flag cross-version
changes automatically, but the interpretation of those changes requires
human judgment. Running at least two integration methods in parallel and
comparing downstream results is a practical way to surface this kind of
sensitivity early, rather than discovering it across sequential pipeline
versions.

### 3.5 The agent cannot manage its own compute environment

The pipeline ran on four different compute configurations over 28 days:

| Phase | RAM | CPUs | GPU | Trigger for change |
|-------|-----|------|-----|--------------------|
| v1 start | 30 GB | 4 | A10G | Initial setup |
| v1 mid | 62 GB | 16 | A10G | OOM during integration |
| v5 start | 62 GB | 16 | A10G | CCA required downsampling |
| v5 mid | 247 GB | 32 | — | Full-cell CCA feasible |

Every one of these transitions was human-initiated. The agent never
proactively identified compute bottlenecks. On multiple occasions, the
pipeline stalled — running single-threaded for 6+ hours with 31 cores
idle, or silently consuming all available memory — and would have remained
stalled indefinitely without human intervention. The human operator had to
notice the stall, diagnose the cause, and either change the compute
environment (adding a GPU, quadrupling RAM from 62GB to 247GB) or
explicitly instruct the agent to investigate and implement parallelization
(e.g., OpenBLAS, `future::plan("multicore")`). Left to its own devices,
the agent would have waited for single-threaded jobs to finish or crashed
repeatedly without escalating.

This is a critical barrier to SME-driven operation. Without someone with
significant computing experience monitoring execution, the pipeline would
not have completed successfully — or would have taken an order of magnitude
longer. A framework intended for use by subject matter experts who are not
compute specialists must either automate environment management or provide
clear escalation paths when jobs stall or fail.

Additionally, R and Python library incompatibilities were a recurring
problem:
- SeuratDisk was broken with Seurat v5 (required a custom bridge)
- STACAS v2.4 changed its API (`SampleIntegration` → `Run.STACAS`)
- Seurat v5 changed its integration API (`FindIntegrationAnchors` →
  `IntegrateLayers`)
- `future.globals.maxSize` default of 16 GB was too small for 410K cells

**Lesson:** Computational requirements are unpredictable at the start of an
analysis. The framework should either (a) provision generously upfront
(10x the expected peak), (b) use a cloud-native architecture that scales
on demand, or (c) include an explicit "environment sizing" module that runs
a subset of the data to estimate resource requirements before committing to
the full run. The agent should be instrumented to detect stalls (e.g., job
running >N hours with <X% CPU utilization) and automatically escalate to
the human operator rather than waiting silently.

R and Python library version conflicts are near-certain. Pin exact versions
in a container image, test them before committing to a workflow, and have
a fallback plan (as the MTX/CSV bridge demonstrated).

---

## 4. Structural Observations

### 4.1 Humans drove architecture; the agent drove execution

The 11 human commits across 3 people changed *what* the pipeline does.
The 85 agent commits changed *how* it gets done. No architectural decision
— annotation timing, compartment structure, integration method selection,
module granularity — was initiated by the agent. Every such decision came
from a domain expert reviewing results and recognizing a flaw the agent
could not see.

This is not a limitation to be fixed. It is the correct division of labor.
The agent lacks the domain knowledge to evaluate whether an annotation
strategy is scientifically sound. It can execute any strategy competently,
but it cannot choose between them.

**Implication for framework design:** The framework should make it *easy*
for domain experts to express architectural decisions (via specs or
configuration) and *hard* for the agent to bypass them. The current spec
format works well for this. What needs improvement is the feedback loop —
making it faster for a domain expert to see what happened, judge whether
it's right, and redirect.

### 4.2 The analysis_plan.md grew unwieldy

As the source of truth for both operational state and historical record,
`analysis_plan.md` grew to 258 lines mixing v5 active status with v1-v4
changelogs, cross-version comparisons, and incident logs. This caused
problems:

- The agent consumed stale cross-version references and reproduced them in
  outputs (the "8/13 NP clusters discordant in v3" appearing in v5 reports)
- The file became difficult for humans to scan for the current state
- Historical context was valuable but cluttered the operational document

The eventual fix (splitting into `analysis_plan.md` for current state and
`docs/version_history.md` for history) should have been the design from
the start.

**Recommendation:** Separate operational state (what's active, what's next)
from historical record (what happened, what was learned) from the
beginning. The agent reads the operational state file; humans read the
history file. Cross-references link them.

### 4.3 Validation caught format errors, not scientific errors

The automated validation checks in each spec verify completeness: do the
output files exist, are they non-empty, do the shapes match. They never
caught:

- The 17K misrouted stressed NP cells (v2→v3)
- The fabricated gene claims in the manuscript (v1)
- The stale method references in the report (v5)
- The fact that CCC direction was unstable across versions

All scientific corrections came from human review. The specs
acknowledged this explicitly: "Passing validation does NOT confirm
scientific correctness."

**Recommendation:** Validation should be expanded to include scientific
sanity checks that *can* be automated:

- Known positive controls (e.g., "COL2A1 must be significantly expressed in
  chondrocyte clusters")
- Cross-method consistency checks (results directionally stable across
  integration methods)
- Anti-hallucination checks (every gene name/fold change in generated text
  matches a row in the result table)
- Distribution checks (cell type proportions within expected ranges for the
  tissue type)

These are not sufficient for scientific correctness, but they would have
caught at least 3 of the 4 errors above.

### 4.4 The notebook/results split created a persistent problem

Scripts write to `results/` (gitignored). Notebooks read from `results/`
for visualization. But when the pipeline runs on a remote machine, the
results aren't available locally, and the report can't be generated.

The eventual fix — extracting key tables from notebook outputs into
`docs/v5_results/` — is a workaround. The correct design would ensure that
all data needed to generate the final report is committed to the repository
(or stored in a durable, accessible location).

**Recommendation:** Summary result tables (DE summary, cell type counts,
trajectory correlations, etc.) should be committed alongside the notebooks,
even if the full result files (individual DE results per cell type) are too
large. The report should be generatable from committed files alone.

---

## 5. Recommendations for a Reusable Framework

### 5.1 Control architecture

```
specs/                  # Human-authored, agent-readable contracts
  00_project.md         # Framework rules, checkpoint protocol
  01_module.md          # Per-module: objective, method, validation, checkpoint
state/
  active.md             # Current module, status, next action (agent reads/writes)
  decisions.md          # All human checkpoint decisions (append-only log)
  issues.md             # Known issues, caveats (append-only)
history/
  changelog.md          # Version history (auto-generated from git)
  sensitivity.md        # Cross-version result stability analysis
```

**Key change from IVD project:** `active.md` should be a minimal,
machine-parseable state file — not a narrative document. Something like:

```yaml
version: 5
module: 08
status: WAITING_FOR_HUMAN_REVIEW
checkpoint_materials:
  - results/differential/de_summary_table.tsv
  - notebooks/08_differential.ipynb
```

The agent updates `status`; the orchestrator reads it. No ambiguity about
whether a checkpoint is active.

### 5.2 Checkpoint enforcement

Checkpoints must be enforced at the orchestration layer, not in the prompt.
The `run_pipeline.sh` pattern works but is fragile (depends on string
matching). A more robust approach:

1. The agent writes a structured state file after each module
2. The orchestrator reads the state file and decides whether to start the
   next iteration
3. If `status == WAITING_FOR_HUMAN_REVIEW`, the orchestrator notifies the
   human (email, Slack, dashboard) and blocks
4. The human reviews, makes decisions, updates the state file (or a
   decisions file), and unblocks
5. The orchestrator starts the next agent iteration

This removes the agent from the control loop entirely. It cannot skip
a checkpoint because it never decides whether to continue.

### 5.3 Anti-hallucination architecture

Every claim in generated text should be backed by a mechanical reference:

1. **Result tables are the source of truth.** Summary tables (DE summary,
   enrichment results, etc.) are committed to the repo.
2. **Generated text references rows in tables.** When the report says
   "COMP is upregulated (log2FC=+2.3, padj=0.001)," those numbers must
   match a specific row in `de_results_combined.tsv`.
3. **A validation step checks references.** After generating a report, a
   separate script verifies that every gene name + fold change + p-value
   in the text appears in the cited source file.
4. **Narrative and data are separated.** The agent generates a data-driven
   report (tables, numbers, source references). Narrative interpretation
   is a separate, human-reviewed step.

### 5.4 Multi-method stability as a quality gate

Instead of running one integration method and hoping it's right, the
framework should:

1. Run at least two integration methods (ideally with different assumptions,
   e.g., one label-free, one semi-supervised)
2. Run downstream analyses (DE, trajectory, CCC) on both
3. Automatically flag any result that changes direction between methods
4. Present the comparison at a human checkpoint
5. Only report results that are directionally stable

This is what v5 partially implemented with CCA vs. scANVI vs. STACAS. The
lesson is that this comparison should be built into the framework, not
added as an afterthought in the fifth iteration.

### 5.5 Computing environment

**Provision for peak, not average.** The IVD project needed 247 GB RAM for
its largest operation (CCA on 410K cells). This was discovered only after
three prior configurations proved insufficient.

Recommendations:
- **Environment sizing module:** Before the full run, execute a
  representative subset (e.g., 2 of 12 datasets) to estimate memory, disk,
  and time requirements. Scale provisioning from there.
- **Cloud-native by default.** Use spot/preemptible instances that can scale
  RAM on demand. The cost of over-provisioning for a few hours is far less
  than the cost of OOM crashes, debugging, and re-runs.
- **Container the environment.** Pin R and Python package versions in a
  Dockerfile. The SeuratDisk/v5 incompatibility, STACAS API change, and
  BLAS configuration issues would all have been avoided with a pre-built
  container.
- **Separate R and Python execution.** Rather than calling R from Python
  (or vice versa) via bridges, run R and Python modules as independent
  steps that communicate through files. This avoids fork-safety issues,
  memory sharing problems, and API incompatibilities.

### 5.6 Making the framework SME-drivable

The goal is a framework that subject matter experts can operate without
needing to write code or debug infrastructure. This requires:

**A. Spec templates, not blank pages.** Provide a library of pre-written
module specs for common scRNA-seq operations (QC, integration, clustering,
annotation, DE, enrichment, trajectory, CCC). The SME selects which modules
to include, fills in dataset-specific parameters (species, tissue,
expected cell types, conditions), and the agent executes.

**B. Visual checkpoint review.** The current checkpoint materials are
notebooks and TSV files. An SME needs a dashboard showing:
- UMAPs colored by cell type, condition, and dataset
- DE volcano plots and summary tables
- A simple approve/redirect/reject interface
- The ability to annotate ("rename cluster 5 to 'stressed_NPC'") without
  editing code

**C. Plain-language spec editing.** When an SME wants to change the
analysis (e.g., "exclude the herniated samples" or "split NP into
subtypes"), they should be able to state this in natural language. The
agent translates it into a spec change, the SME confirms, and the agent
re-executes.

**D. Provenance by default.** Every result should be traceable to the
code, parameters, and data that produced it — without the SME needing to
manage git. The framework should auto-commit after each module, tag
versions, and maintain a human-readable changelog.

**E. Graceful re-runs.** When the SME redirects (e.g., changes the
annotation), the framework should identify which downstream modules are
invalidated and offer to re-run them, rather than requiring a full pipeline
restart. The IVD project re-ran everything from scratch for each version;
a smarter framework would cache and reuse unaffected upstream results.

---

## 6. What We Would Change If Starting Over

1. **Start with two integration methods from the beginning.** The v1→v5
   journey was largely a search for the right integration approach. Running
   two methods in parallel from v1 and comparing at a checkpoint would have
   saved at least two full pipeline iterations.

2. **Enforce checkpoints mechanically from day one.** The shell gate should
   have been the first thing built, not a retroactive fix after the agent
   ran away.

3. **Commit summary result tables.** The report generation problem (results
   on a remote machine, not available locally) was entirely avoidable.
   Summary tables are small enough to commit.

4. **Use containers.** Every R/Python version conflict was avoidable with
   a pinned Docker image.

5. **Separate narrative from data in the report.** The v1 manuscript
   hallucination problem could have been prevented by never asking the
   agent to write interpretive prose. Instead: data tables with source
   references, reviewed by a human who writes the interpretation.

6. **Build the sensitivity analysis into the pipeline.** Cross-version
   result stability was tracked manually and recognized late. It should be
   an automated quality gate.

7. **Use a structured state file, not a markdown narrative.** The
   analysis_plan.md approach works but introduces ambiguity. A YAML or
   JSON state file, read by the orchestrator, eliminates the possibility of
   the agent misinterpreting its own state.

8. **Provision compute generously upfront.** The four-configuration journey
   (30GB → 62GB → 62GB+GPU → 247GB) consumed multiple days of debugging.
   Starting with 256GB and a GPU would have cost marginally more and saved
   significant time.

---

## 7. The Division of Labor That Emerged

After five iterations, a natural division crystallized:

| Role | Responsibility | How expressed |
|------|---------------|---------------|
| **Domain expert** | What to analyze, how to interpret, when to redirect | Spec edits, checkpoint decisions |
| **Agent** | How to execute, how to debug, how to validate | Scripts, notebooks, commits |
| **Orchestrator** | When to run, when to stop, when to notify | Shell gate, state file |

The key insight is that these three roles should be explicitly designed
into the framework, not discovered through failure. The domain expert
should never need to write code. The agent should never make scientific
judgment calls. The orchestrator should never depend on the agent's
self-reported state.

---

## 8. Open Questions

- **How much can annotation be automated?** The IVD project required three
  human-directed annotation revisions. Would a reference atlas (CellTypist,
  Azimuth) have reduced this, or is tissue-specific annotation inherently
  a human-judgment problem?

- **Can the framework detect when it needs human input?** The agent never
  flagged that its annotation was wrong — a human always caught it. Could
  disagreement between two methods (e.g., CellTypist vs. de novo) trigger
  an automatic checkpoint?

- **What is the minimum viable checkpoint?** Some modules (QC, trajectory)
  had checkpoints that were pro forma approvals. Could the framework
  auto-advance past low-risk modules and only block at high-leverage
  points (annotation, integration selection, DE interpretation)?

- **How do we handle the SME's time?** The framework assumes a domain
  expert is available to review checkpoints promptly. In practice, review
  delays were the longest pauses in the pipeline. How do we minimize the
  burden on the SME while maintaining scientific rigor?

---

*This document reflects the experience of one project (IVD scRNA-seq atlas,
12 datasets, 423K cells, 5 pipeline versions over 28 days). Findings may
not generalize to all analysis types, but the failure modes — checkpoint
evasion, claim fabrication, annotation fragility, compute under-provisioning
— are likely universal in human-gated agentic analysis pipelines.*
