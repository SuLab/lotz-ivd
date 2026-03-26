# v2 System Design: SME-Driven Agentic Analysis

> A generalizable framework for letting subject matter experts create and run
> complex data-intensive analyses using a coding agent — without writing code,
> managing infrastructure, or learning pipeline tooling.
>
> The framework is domain-agnostic. It was developed from an scRNA-seq atlas
> project but is designed to support any structured analysis where a subject
> matter expert defines the question and a computational pipeline produces
> the answer: proteomics, clinical trial evaluation, time series analysis,
> imaging, epidemiology, or any domain where expert judgment and
> computational execution must stay coupled.

---

## 1. Design Principles

Five principles derived from the IVD atlas project (5 pipeline versions,
12 modules, ~423K cells, 28 days). These are domain-agnostic — they apply
whether the analysis is single-cell genomics, proteomics, clinical trial
evaluation, or time series forecasting.

1. **Humans decide what; the agent decides how.** Every scientific or
   analytical judgment — what cell types exist, whether a treatment effect
   is real, which model specification to use — is made by a human at a
   checkpoint. The agent writes code, debugs errors, runs validation, and
   generates review materials. It never makes analytical calls
   autonomously.

2. **Specs are the contract.** The interface between the SME and the agent
   is a structured specification document, not code. The SME describes what
   they want in structured natural language; the agent translates that into
   executable scripts. If the spec is wrong, the SME changes the spec — not
   the code.

3. **Checkpoints are mechanical, not advisory.** The orchestrator enforces
   checkpoints. The agent cannot skip them, because the agent does not
   decide whether to continue. The orchestrator reads a structured state
   file, and if the state says "waiting for review," the loop stops.

4. **Results must be verifiable, not just plausible.** Every quantitative
   claim in any generated output must be traceable to a specific row in a
   result table. Narrative and data are separated. Fact-checking is
   automated.

5. **Instability is a signal, not a bug.** Results that change across
   methods, parameter choices, or sensitivity analyses are flagged
   automatically. Where appropriate, the framework runs multiple methods
   and reports which findings are robust and which are method-sensitive.

---

## 2. Actors

Three roles, explicitly separated. One person can fill multiple roles, but
the framework treats them as distinct.

| Role | Responsibility | Never does |
|------|---------------|------------|
| **SME** | Defines the analysis question. Selects modules. Fills in domain-specific parameters. Reviews checkpoints. Makes all scientific/analytical decisions. | Writes code. Debugs infrastructure. Manages compute. |
| **Agent** | Translates specs into scripts. Executes code. Debugs runtime errors. Generates notebooks and review materials. Runs validation. | Makes scientific judgments. Chooses between methods without being told. Writes interpretive prose without mechanical verification. |
| **Orchestrator** | Manages the agent loop. Enforces checkpoints. Monitors compute health. Notifies humans. Manages versioning and provenance. | Makes scientific decisions. Modifies code or specs. |

In the IVD project, a fourth role emerged empirically: the **operator** — a
person with computing expertise who diagnosed stalls, managed the compute
environment, and directed the agent through infrastructure problems. In v2,
the orchestrator subsumes most of this role through automated monitoring
and environment management (Section 9). But for analyses that push hardware
limits, an operator may still be needed. The framework should make it
obvious when this is the case by escalating clearly rather than stalling
silently.

---

## 3. Project Lifecycle

```
┌─────────────────────────────────────────────────────────────┐
│                     SME: PROJECT SETUP                       │
│  1. Choose analysis type (template)                          │
│  2. Describe the biological question                         │
│  3. Specify datasets, species, conditions                    │
│  4. Select modules (from template library)                   │
│  5. Set checkpoint mode (fast / gated / selective)           │
└─────────────────────┬───────────────────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────────────────┐
│                  AGENT: INITIALIZATION                        │
│  1. Generate project scaffold (dirs, state files, specs)     │
│  2. Run environment sizing on data subset                    │
│  3. Provision compute (or recommend instance size)           │
│  4. Download / validate input data                           │
│  5. Report readiness → first checkpoint                      │
└─────────────────────┬───────────────────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────────────────┐
│                 EXECUTION LOOP (per module)                   │
│                                                              │
│  ┌──────────┐    ┌───────────┐    ┌──────────────────────┐  │
│  │ Agent    │───▶│ Validate  │───▶│ Checkpoint           │  │
│  │ executes │    │ outputs   │    │ (if mode requires)   │  │
│  │ module   │    │           │    │                      │  │
│  └──────────┘    └───────────┘    └──────────┬───────────┘  │
│                                              │              │
│                        ┌─────────────────────┤              │
│                        │                     │              │
│                        ▼                     ▼              │
│                  ┌──────────┐          ┌──────────┐         │
│                  │ Advance  │          │ Redirect │         │
│                  │ to next  │          │ (SME     │         │
│                  │ module   │          │  changes │         │
│                  └──────────┘          │  spec)   │         │
│                                       └──────────┘         │
└─────────────────────────────────────────────────────────────┘
```

The loop runs until all modules are complete or the SME terminates the
analysis. At any checkpoint, the SME can:

- **Advance:** Accept results, proceed to the next module.
- **Redirect:** Change parameters, swap methods, add/remove modules. The
  agent regenerates affected specs and re-executes.
- **Restart module:** Re-run the current module with different settings.
- **Fork:** Branch the analysis to compare two approaches in parallel
  (e.g., two integration methods), then reconverge at a later checkpoint.

---

## 4. Project Structure

```
project/
├── project.yaml              # Project-level config (analysis type, species,
│                              # conditions, checkpoint mode, compute profile)
├── specs/
│   ├── 00_project.md          # Framework rules, generated from template
│   ├── 01_module.md           # Per-module specs, generated from template
│   ├── ...                    #   + SME customizations
│   └── custom_module.md       # SME-added modules (optional)
├── state/
│   ├── pipeline.yaml          # Machine-readable pipeline state (Section 5)
│   ├── decisions.yaml         # Append-only log of checkpoint decisions
│   └── issues.yaml            # Known issues, caveats, flags
├── scripts/                   # Agent-generated compute scripts
├── notebooks/                 # Agent-generated review notebooks
├── data/
│   ├── raw/                   # Input data (immutable after import)
│   ├── processed/             # Intermediate outputs
│   └── derived/               # Cross-dataset / derived results
├── results/
│   ├── tables/                # Summary result tables (committed to repo)
│   └── figures/               # Generated figures
├── history/
│   ├── versions.md            # Human-readable version history
│   └── sensitivity.yaml       # Cross-method / cross-version stability
├── AGENT.md                   # Agent execution rules (generated)
├── PROMPT.md                  # Loop prompt (generated)
└── orchestrator.yaml          # Orchestrator config (notifications, compute)
```

### Key changes from v1

- **`state/pipeline.yaml` replaces `analysis_plan.md`** as the
  machine-readable state. The agent reads and writes structured YAML, not
  free-form markdown. This eliminates the ambiguity that caused the agent
  to misread its own state in the IVD project.

- **`state/decisions.yaml` is append-only.** Every human checkpoint
  decision is logged with a timestamp, the decision, and the rationale.
  This is the audit trail. The agent never modifies past entries.

- **`results/tables/` is committed.** Summary tables (DE summaries, cell
  counts, enrichment results) are small enough to commit. The final report
  can be generated from committed files alone — solving the "results on a
  remote machine" problem.

- **`history/sensitivity.yaml`** tracks cross-method and cross-version
  result stability automatically. Any result that changes direction between
  methods is flagged before it reaches the SME.

---

## 5. State Management

The pipeline state is a structured YAML file, not a markdown narrative.
The orchestrator is the only component that decides whether to start a new
agent iteration — based solely on this file.

```yaml
# state/pipeline.yaml
project: "IVD Single-Cell Atlas"
version: 5
checkpoint_mode: selective    # fast | gated | selective

modules:
  - id: "01_dataset_discovery"
    status: complete           # pending | running | complete | failed |
                               # waiting_for_review | skipped
    completed_at: "2026-02-26"
    outputs:
      - data/raw/manifest.tsv
    checkpoint: false          # not a designated checkpoint in selective mode

  - id: "05_integration"
    status: waiting_for_review
    started_at: "2026-03-23"
    checkpoint: true
    checkpoint_materials:
      - notebooks/05d_integration_comparison.ipynb
      - results/tables/integration_metrics.tsv
      - results/figures/umap_by_method.png
    review_prompt: >
      Three integration methods ran on all compartments. Review the
      comparison notebook and metrics table. Select which workflow to
      carry forward, or request additional methods.

active_module: "05_integration"
next_module: "06_clustering"
```

### Checkpoint modes

| Mode | Behavior | When to use |
|------|----------|-------------|
| **Fast** | Run all modules without stopping. All checkpoints auto-advance. | First pass on a new dataset. Rerunning after minor parameter tweaks. |
| **Gated** | Stop at every module's checkpoint. | New analyses. After major spec changes. When the SME wants maximum control. |
| **Selective** | Stop only at designated high-leverage checkpoints (annotation, integration selection, DE review). Auto-advance past low-risk modules (QC, trajectory, CCC). | Established pipelines where the SME knows which decisions matter most. |

The SME sets the mode at project creation and can override it per-module.
The orchestrator enforces it — the agent has no knowledge of checkpoint
modes.

---

## 6. Spec Format

Specs are the contract between the SME and the agent. They define *what*
the agent should do, not *how* to code it. The format is structured enough
for the agent to execute but readable enough for the SME to author and
review.

### Template

```markdown
# Module NN: [Module Name]

## Objective
[One paragraph: what this module accomplishes and why it matters.]

## Inputs
- [List of input files / objects with expected format]

## Method
[Structured description of what the agent should do. Not pseudocode —
describe the analytical steps, parameters, and decision points.]

### Parameters
| Parameter | Value | Rationale |
|-----------|-------|-----------|
| min_genes | 200 | Standard QC threshold |
| ...       | ...   | ...       |

### Alternatives (optional)
[If the SME wants multiple methods compared, list them here with the
comparison criteria.]

## Outputs
- [List of expected output files with format and location]

## Validation
[Automated checks the agent must pass before declaring the module
complete. These are necessary but not sufficient for scientific
correctness.]

- [ ] Output files exist and are non-empty
- [ ] Cell counts within expected range
- [ ] [Domain-specific sanity checks]

## Scientific Sanity Checks (optional)
[Domain-specific checks that go beyond format validation — positive
controls, expected distributions, known relationships. These do not
block progression but are reported at the checkpoint. Examples:]

- [ ] [Expected entity] present in [expected context]
- [ ] [Control condition] yields [expected result]
- [ ] [Distribution/proportion] within [expected range]

## Human Checkpoint
**Checkpoint level:** required | recommended | auto-advance
**Review materials:**
- [List of notebooks, tables, and figures for the SME to review]
**Decision required:**
- [What the SME needs to decide at this checkpoint]
```

### Spec authoring workflow

1. The SME selects a module from the template library (or creates a
   custom module).
2. The framework presents a questionnaire for the module parameters:
   species, tissue type, expected cell types, conditions, etc.
3. The agent generates a draft spec from the template + answers.
4. The SME reviews and edits the spec (in plain language — the agent
   handles any technical translation).
5. The spec is committed and becomes the contract.

At any checkpoint, the SME can edit a spec. The agent detects the change,
identifies affected downstream modules, and offers to re-execute them.

---

## 7. Analysis Types and Module Libraries

The framework is organized around **analysis types** — domain-specific
templates that define the module set, dependency graph, default
parameters, validation checks, and checkpoint recommendations for a class
of analysis. The SME selects an analysis type, customizes it, and the
framework generates the project scaffold.

An analysis type is a reusable package containing:

- A default module sequence with dependency graph
- Spec templates for each module
- Domain-specific validation checks (Tier 2)
- Domain-specific claim types for the verifier
- A project.yaml template with the right domain fields
- A container image (or environment spec) with pinned dependencies

### Defining a new analysis type

Any completed project can be generalized into a new analysis type. The
framework extracts the module structure, dependency graph, and spec
templates from the project's specs — replacing project-specific values
with parameterized placeholders. This is how the scRNA-seq template was
created from the IVD project, and how future templates would be created
from the first successful proteomics or clinical trial project.

### Example: scRNA-seq atlas

| Module | Purpose | Checkpoint level |
|--------|---------|-----------------|
| Data acquisition | Download and validate input datasets | Recommended |
| Metadata harmonization | Standardize sample/condition annotations | Required |
| QC & preprocessing | Filter cells/genes, normalize, find HVGs | Recommended |
| Coarse annotation | Broad cell type classification | Required |
| Integration | Merge datasets, correct batch effects | Required |
| Clustering | Identify cell communities | Recommended |
| Fine annotation | Assign cell type labels to clusters | Required |
| Differential expression | Compare conditions (pseudobulk) | Required |
| Pathway / functional enrichment | ORA, GSEA, TF activity | Recommended |
| Trajectory analysis | Pseudotime, RNA velocity | Recommended |
| Cell-cell communication | Ligand-receptor interactions | Recommended |
| Reporting | Summary tables, figures, narrative | Required |

### Example: quantitative proteomics

| Module | Purpose | Checkpoint level |
|--------|---------|-----------------|
| Data import | Load instrument output (MaxQuant, DIA-NN, Spectronaut) | Recommended |
| Quality assessment | Missingness, CV, intensity distributions, batch effects | Required |
| Normalization & imputation | Median normalization, KNN/MinProb imputation | Required |
| Differential abundance | Limma, MSstats, or mixed models | Required |
| Functional enrichment | ORA/GSEA on protein-level results | Recommended |
| PTM analysis | Phospho/acetyl site quantification (if applicable) | Recommended |
| Network / interaction analysis | STRING, protein complex enrichment | Recommended |
| Reporting | Summary tables, volcano plots, heatmaps | Required |

### Example: clinical trial evaluation

| Module | Purpose | Checkpoint level |
|--------|---------|-----------------|
| Data import & CDISC mapping | Load CRF/EDC data, map to SDTM/ADaM | Required |
| Population definition | ITT, mITT, per-protocol populations | Required |
| Baseline characteristics | Demographics table (Table 1), balance checks | Required |
| Primary endpoint analysis | Pre-specified primary analysis per SAP | Required |
| Secondary endpoints | Multiplicity-adjusted secondary analyses | Required |
| Sensitivity analyses | Per-protocol, subgroup, missing data handling | Required |
| Safety analysis | AE tables, exposure, lab shifts | Required |
| Subgroup analyses | Forest plots, interaction tests | Recommended |
| Reporting | CSR tables, figures, listings (TFLs) | Required |

### Example: time series / longitudinal analysis

| Module | Purpose | Checkpoint level |
|--------|---------|-----------------|
| Data import & cleaning | Load time series, handle missing/irregular intervals | Recommended |
| Exploratory analysis | Trends, seasonality, autocorrelation, stationarity tests | Required |
| Feature engineering | Lags, rolling statistics, external covariates | Recommended |
| Model specification | ARIMA, state space, mixed effects, or ML models | Required |
| Model fitting & diagnostics | Fit model, check residuals, cross-validation | Required |
| Forecasting / inference | Point forecasts, intervals, or causal estimates | Required |
| Sensitivity analysis | Alternative specifications, holdout periods | Recommended |
| Reporting | Forecast plots, model comparisons, diagnostics | Required |

### Module structure: DAGs vs. exploratory workflows

The scRNA-seq and proteomics examples are naturally **linear DAGs** —
each module depends on prior modules in a mostly sequential chain.
Clinical trial evaluation and time series analysis are more often
**exploratory** — the SME may want to run a subgroup analysis that wasn't
planned, or re-specify a model after seeing diagnostics.

The framework handles this through two mechanisms:

1. **Optional modules.** The SME can add or skip modules at any
   checkpoint. The dependency graph determines what's invalidated, but
   new modules can be inserted without restarting.

2. **Branching.** At any checkpoint, the SME can fork the analysis into
   parallel branches (e.g., "run both a Cox model and a parametric AFT
   model"), compare results, and select one to carry forward. This is
   the same "multi-method comparison" concept from Section 11, but
   applied to the workflow structure itself.

The framework does not force linearity. The module dependency graph is a
DAG, but it can have multiple roots, parallel branches, and optional
leaves. The orchestrator handles any valid DAG.

---

## 8. Confirmatory vs. Exploratory Analysis Modes

Different domains have fundamentally different relationships between the
analysis plan and the analysis execution. The framework must support
both.

### Exploratory mode

In exploratory analyses (scRNA-seq atlases, proteomics discovery, time
series modeling), the SME defines a starting plan but expects to revise
it based on intermediate results. The plan is a hypothesis — the results
may contradict it. This is the mode the IVD project operated in: five
pipeline versions, each triggered by a human recognizing something that
needed to change.

In exploratory mode:
- The SME can add, remove, or reorder modules at any checkpoint
- Method comparisons are encouraged — run multiple approaches and select
- The dependency graph handles invalidation and re-execution
- The final report describes what was done and why, including revisions

### Confirmatory mode

In confirmatory analyses (clinical trials, pre-registered studies,
regulated submissions), the analysis plan is fixed before execution
begins. Deviations from the plan are permitted but must be documented,
justified, and clearly distinguished from pre-specified analyses. The
Statistical Analysis Plan (SAP) is a legal document, not a hypothesis.

In confirmatory mode:
- The module sequence and methods are locked at project creation
- The SME can still redirect at checkpoints, but deviations are flagged
  as **post-hoc** in the audit trail and the final report
- Sensitivity analyses are pre-specified in the plan, not added ad hoc
- The framework tracks which analyses were pre-specified vs. exploratory
- The audit trail (decisions.yaml) records the timestamp and rationale
  for every deviation from the original plan

### Hybrid mode

Many real projects are hybrid: a confirmatory primary analysis with
exploratory secondary objectives. The framework supports this by
allowing per-module mode flags:

```yaml
modules:
  - id: primary_endpoint
    mode: confirmatory        # locked to SAP specification
  - id: subgroup_analysis
    mode: confirmatory        # pre-specified subgroups only
  - id: biomarker_exploration
    mode: exploratory         # SME can redirect freely
```

---

## 9. Orchestrator

The orchestrator is the control plane. It manages the agent loop,
enforces checkpoints, monitors compute health, and handles notifications.
It is a separate process from the agent — the agent cannot influence the
orchestrator's decisions.

### Core responsibilities

```
┌─────────────────────────────────────────────────────┐
│                   ORCHESTRATOR                       │
│                                                      │
│  ┌─────────────┐  ┌──────────────┐  ┌────────────┐ │
│  │ Loop        │  │ Checkpoint   │  │ Compute    │ │
│  │ controller  │  │ enforcer     │  │ monitor    │ │
│  │             │  │              │  │            │ │
│  │ Reads state │  │ Blocks on    │  │ Detects    │ │
│  │ Starts/     │  │ review.      │  │ stalls,    │ │
│  │ stops agent │  │ Notifies SME │  │ OOM, idle  │ │
│  │ iterations  │  │ via configured│  │ CPUs.      │ │
│  │             │  │ channels     │  │ Escalates. │ │
│  └─────────────┘  └──────────────┘  └────────────┘ │
│                                                      │
│  ┌─────────────┐  ┌──────────────┐  ┌────────────┐ │
│  │ Provenance  │  │ Environment  │  │ Sensitivity│ │
│  │ manager     │  │ manager      │  │ tracker    │ │
│  │             │  │              │  │            │ │
│  │ Auto-commit │  │ Sizes compute│  │ Compares   │ │
│  │ after each  │  │ from pilot   │  │ results    │ │
│  │ module. Tags│  │ run. Scales  │  │ across     │ │
│  │ versions.   │  │ on demand.   │  │ methods.   │ │
│  └─────────────┘  └──────────────┘  └────────────┘ │
└─────────────────────────────────────────────────────┘
```

### Loop control

```python
# Pseudocode — the actual implementation may be a shell script,
# a Python daemon, or a cloud-native workflow engine.

while not pipeline_complete(state):
    state = read_state("state/pipeline.yaml")

    if state.active_module.status == "waiting_for_review":
        notify_sme(state.active_module)     # email, Slack, dashboard
        block_until_review_complete()
        state = read_state("state/pipeline.yaml")  # SME updated it

    if state.active_module.status == "failed":
        if retry_count < max_retries:
            retry_module()
        else:
            notify_operator(state.active_module)
            block_until_resolved()

    # Compute health check
    if detect_stall(cpu_util < 5%, duration > 30min):
        notify_operator("Possible stall detected")

    # Start next agent iteration
    run_agent(prompt="PROMPT.md", context=state)

    # Auto-commit after module completion
    auto_commit(state.active_module)
```

### Compute monitoring

The orchestrator monitors resource utilization during agent execution and
acts on anomalies:

| Signal | Threshold | Action |
|--------|-----------|--------|
| CPU utilization < 5% for > 30 min | Possible stall | Alert operator |
| Memory > 90% for > 5 min | Approaching OOM | Alert operator, suggest scale-up |
| GPU idle while training | Misconfigured job | Alert with fix suggestion |
| Job running > 2x expected duration | Possible inefficiency | Alert with profiling suggestion |
| Disk > 90% | Running out of space | Alert, suggest cleanup |

These thresholds are configurable. The point is that silent stalls — the
single biggest time sink in the IVD project — are detected automatically.

### Environment sizing

Before the full pipeline runs, the orchestrator executes an **environment
sizing module** on a subset of the data (e.g., 2 of 12 datasets, or a
10% subsample). This estimates:

- Peak memory per module
- GPU memory requirements (if applicable)
- Estimated wall time per module
- Disk space requirements

The orchestrator then provisions or recommends an instance accordingly.
The IVD project went through four compute configurations (30GB → 62GB →
62GB+GPU → 247GB) — all human-initiated. The sizing module prevents this.

---

## 10. Anti-Hallucination Architecture

The IVD project's v1 manuscript contained 14 fabricated gene-level claims.
This architecture prevents that class of error.

### Principle: separate data from narrative

The agent produces two distinct types of output:

1. **Data outputs:** Tables, figures, metrics — mechanically generated
   from computation. These are the source of truth.
2. **Narrative outputs:** Text that describes, summarizes, or interprets
   the data outputs. Every claim in narrative text must reference a
   specific data output.

### Claim verification pipeline

```
Agent generates report text
         │
         ▼
┌─────────────────────────────────┐
│  CLAIM EXTRACTOR                │
│                                 │
│  Parses text for domain-        │
│  specific verifiable entities:  │
│  - Named entities (genes,       │
│    proteins, drugs, endpoints)  │
│  - Effect sizes (log2FC, HR,    │
│    beta coefficients, AUC)      │
│  - Statistical measures         │
│    (p-values, CIs, NNT)         │
│  - Direction claims (up/down,   │
│    improved/worsened, sig/NS)   │
│  - Category names (cell types,  │
│    treatment arms, subgroups)   │
│  - Counts and proportions       │
└──────────┬──────────────────────┘
           │
           ▼
┌─────────────────────────────────┐
│  CLAIM VERIFIER                 │
│                                 │
│  For each extracted claim:      │
│  1. Look up the cited source    │
│     table                       │
│  2. Find the matching row       │
│  3. Compare extracted values    │
│     to table values             │
│  4. Flag mismatches             │
└──────────┬──────────────────────┘
           │
           ▼
┌─────────────────────────────────┐
│  VERIFICATION REPORT            │
│                                 │
│  - N claims checked             │
│  - N verified                   │
│  - N mismatches (with details)  │
│  - N unverifiable (no source)   │
└─────────────────────────────────┘
```

The claim extractor is configured per analysis type. Each analysis type
defines its **verifiable entity types** — the categories of quantitative
claims that can be mechanically checked against result tables:

| Analysis type | Entity types | Effect size types | Direction vocabulary |
|---------------|-------------|-------------------|---------------------|
| scRNA-seq | Gene names, cell types, pathways | log2FC, enrichment score | upregulated/downregulated, enriched/depleted |
| Proteomics | Protein names, PTM sites, complexes | log2FC, intensity ratio | abundant/depleted, modified/unmodified |
| Clinical trial | Treatment arms, endpoints, subgroups | Hazard ratio, odds ratio, mean difference, NNT | improved/worsened, superior/non-inferior |
| Time series | Variables, time periods, model terms | Coefficients, forecast values, RMSE | increasing/decreasing, significant/NS |

### Rules for agent-generated text

1. Every named entity + statistic pair in generated text must
   correspond to a row in a committed result table.
2. Direction claims ("upregulated," "superior," "increasing") must
   match the sign/direction of the corresponding statistic in the table.
3. Aggregation claims ("consistently upregulated across all cell types,"
   "significant in all subgroups") must be verified against every
   applicable row, not sampled.
4. The claim verifier runs as a validation step after report generation
   — before the checkpoint. If any claim fails verification, the report
   is flagged and the mismatches are shown to the SME.

### What the agent should and should not write

| Agent writes | Human writes |
|-------------|-------------|
| Data-driven summaries: "N entities significant at threshold X" | Domain interpretation: what the results mean |
| Table references: "Top results by effect size (Table SN)" | Mechanistic or causal hypotheses |
| Method descriptions: parameters, software versions, assumptions | Discussion, limitations, and conclusions |
| Comparisons between methods: metrics, overlap, concordance | Claims about real-world significance or actionability |

---

## 11. Robustness and Sensitivity Analysis

Results that depend on a single analytical choice are fragile. The
framework builds sensitivity analysis into the pipeline as a first-class
concept — but the form of that analysis varies by domain.

### Two modes of robustness checking

**Mode 1: Multi-method comparison (exploratory analyses).** When the
analysis involves method choice (which integration algorithm, which
clustering resolution, which normalization), the spec lists 2+ methods.
The agent runs all of them and compares. This is the default for
exploratory analyses like scRNA-seq, proteomics, and time series, where
the "right" method is not known a priori.

**Mode 2: Pre-specified sensitivity analyses (confirmatory analyses).**
When the analysis has a pre-specified primary method — as in a clinical
trial governed by a Statistical Analysis Plan (SAP) — the framework does
not run multiple methods and pick the best one. Instead, it runs the
primary analysis exactly as specified, then runs pre-defined sensitivity
analyses (per-protocol population, alternative missing-data handling,
tipping-point analysis) and reports whether the primary result is robust
to these perturbations.

The SME specifies which mode applies at the analysis-type level or
per-module. Both modes produce the same output: a classification of each
result as robust, concordant, or discordant.

### How it works

1. **At designated comparison/sensitivity points,** the spec lists the
   methods or sensitivity analyses to run.
2. The agent runs all listed analyses and produces standardized outputs.
3. A comparison module runs automatically:
   - Overlap of significant results between methods/analyses
   - Directional concordance (do results go in the same direction?)
   - Effect size correlation
   - Higher-level concordance (pathways, endpoints, subgroups)
4. Results are classified:

| Classification | Criterion | Reporting |
|---------------|-----------|-----------|
| **Robust** | Same direction and significance across all methods/analyses | Report as a finding |
| **Concordant** | Same direction, significance varies | Report with caveat |
| **Discordant** | Direction or conclusion changes | Flag as method-sensitive or assumption-sensitive |

5. The stability report is presented at the checkpoint.

### Domain-specific sensitivity examples

| Domain | Primary analysis | Sensitivity analyses |
|--------|-----------------|---------------------|
| scRNA-seq | CCA integration | scANVI, STACAS — compare cell types, DE results |
| Proteomics | Median normalization + limma | Quantile normalization, MSstats — compare protein hits |
| Clinical trial | ITT with MMRM (per SAP) | Per-protocol, LOCF, tipping point, subgroup interactions |
| Time series | ARIMA(p,d,q) | Alternative order selection, structural breaks, holdout validation |

### Cross-version tracking

If the pipeline is re-run (new version), the orchestrator automatically
compares the new results to the prior version's results using the same
concordance framework. Changes are flagged and explained in context. The
agent generates mechanistic explanations for the changes; the SME decides
whether they matter.

---

## 12. SME Interface

The framework must be operable by someone who has domain expertise but
does not write code, manage git, or configure cloud infrastructure. This
section describes what the SME sees and does.

### 12.1 Project creation

The SME interacts with a setup wizard (CLI, web, or chat-based). The
wizard adapts its questions based on the selected analysis type. The
first question is always the same; everything after is domain-specific.

**Common questions (all analysis types):**
```
1. What type of analysis?
   [scRNA-seq atlas / proteomics / clinical trial / time series / ...]

2. What is the core question?
   [Free text — the SME's research question in their own words]

3. Checkpoint mode?
   [Fast / Gated / Selective]
```

**Domain-specific questions (examples):**

*scRNA-seq:*
```
- Species and tissue?
- Dataset accessions?
- Conditions to compare?
- Expected cell types and markers? (for sanity checks)
```

*Proteomics:*
```
- Instrument and acquisition mode? (DDA, DIA, TMT, label-free)
- Search engine output format? (MaxQuant, DIA-NN, Spectronaut, FragPipe)
- Experimental design? (conditions, replicates, batches)
- Proteins or PTM sites of interest? (for sanity checks)
```

*Clinical trial:*
```
- Study design? (parallel, crossover, factorial, adaptive)
- Phase and regulatory context?
- Primary endpoint and analysis method? (from SAP)
- Secondary endpoints?
- Key populations? (ITT, mITT, per-protocol)
- Subgroups of interest?
```

*Time series:*
```
- Data frequency? (daily, monthly, irregular)
- Forecast horizon or comparison period?
- Known external factors / covariates?
- Stationarity expectations?
```

From these answers, the framework generates:
- `project.yaml` with all parameters
- A full set of module specs from templates, customized with the SME's
  parameters
- `AGENT.md` and `PROMPT.md`
- The directory scaffold

### 12.2 Checkpoint review

When a checkpoint fires, the SME receives a notification (email, Slack,
dashboard — configurable) with:

- A plain-language summary of what the module did
- Links to the review notebook (rendered as HTML — no Jupyter needed)
- Key metrics in a summary table
- The specific decision the SME needs to make

The SME responds through a review interface that supports:

```
┌─────────────────────────────────────────────────────┐
│  MODULE 07: Post-Integration Annotation              │
│                                                      │
│  Status: WAITING FOR REVIEW                          │
│                                                      │
│  Summary: Identified 5 cell types in NP, 4 in AF,   │
│  7 in CEP. Mature chondrocytes dominate NP (72%).    │
│                                                      │
│  [View Notebook]  [View Marker Table]  [View UMAPs]  │
│                                                      │
│  ┌───────────────────────────────────────────────┐   │
│  │ Decision:                                      │   │
│  │ ○ Approve and advance                         │   │
│  │ ○ Approve with modifications:                 │   │
│  │   [Rename cluster 3 to "stressed_NPC"]        │   │
│  │   [Merge clusters 7 and 8]                    │   │
│  │   [Add a custom marker: ACAN]                 │   │
│  │ ○ Reject — rerun with different parameters    │   │
│  │ ○ Reject — change approach (edit spec)        │   │
│  └───────────────────────────────────────────────┘   │
│                                                      │
│  Notes for the record:                               │
│  [Free text — logged in decisions.yaml]              │
│                                                      │
│  [Submit Decision]                                   │
└─────────────────────────────────────────────────────┘
```

### 12.3 Spec editing via natural language

When the SME wants to change the analysis, they state it in natural
language:

> "Exclude the herniated samples from the degeneration comparison."
> "Split NP mesenchymal cells into subtypes before DE."
> "Add STACAS as a third integration method."
> "Use Milo instead of propeller for composition analysis."

The agent translates this into a spec change, shows the diff to the SME
for confirmation, and re-executes affected modules. The SME never edits
YAML or code directly.

### 12.4 Provenance (invisible to the SME)

The framework auto-manages provenance:

- Auto-commits after each module with a descriptive message
- Tags each pipeline version (`v1`, `v2`, etc.)
- Maintains `history/versions.md` with a human-readable changelog
- Tracks file checksums for data integrity
- Logs every checkpoint decision with timestamp and rationale

The SME gets reproducibility for free. Any result can be traced to the
exact code, parameters, data, and decisions that produced it.

---

## 13. Dependency-Aware Re-Execution

When the SME redirects at a checkpoint, the framework must know which
downstream modules are invalidated. The IVD project re-ran everything
from scratch for each version — a smarter framework caches and reuses
unaffected upstream results.

### Module dependency graph

```yaml
# Defined per analysis type. Example for scRNA-seq:
dependencies:
  01_data_acquisition: []
  02_metadata: [01_data_acquisition]
  03_preprocessing: [01_data_acquisition, 02_metadata]
  04_coarse_annotation: [03_preprocessing]
  05_integration: [03_preprocessing, 04_coarse_annotation]
  06_clustering: [05_integration]
  07_annotation: [06_clustering]
  08_differential: [07_annotation, 02_metadata]
  09_interpretation: [08_differential]
  10_trajectory: [05_integration, 07_annotation]
  11_communication: [07_annotation, 02_metadata]
  12_reporting: [08_differential, 09_interpretation, 10_trajectory, 11_communication]
```

### Re-execution logic

When a module is re-executed (due to spec change or SME redirect):

1. Mark the module as `pending`.
2. Walk the dependency graph forward: mark all downstream modules as
   `invalidated`.
3. Present the invalidation list to the SME: "Re-running Module 05 will
   invalidate Modules 06–12. Proceed?"
4. On confirmation, archive current results (tagged with version) and
   re-execute from the changed module forward.
5. Modules upstream of the change are not re-executed — their cached
   outputs are reused.

### Partial invalidation

Some changes only affect a subset of downstream modules. For example,
changing the annotation spec (Module 07) invalidates DE (08), pathways
(09), trajectory (10), CCC (11), and reporting (12) — but not integration
(05) or clustering (06). The dependency graph captures this precisely.

---

## 14. Compute Environment Strategy

### Containerized by default

The framework ships pre-built container images for each analysis type.
The container pins exact versions of all dependencies — eliminating the
version conflicts that consumed days in the IVD project (Seurat
v5/SeuratDisk incompatibility, STACAS API changes, BLAS configuration).

Each analysis type has a corresponding container that is tested before
release — all cross-language bridges and tool interactions verified.

### Sizing strategy

Compute requirements vary enormously across analysis types. An scRNA-seq
atlas with 400K cells needs 247 GB RAM; a clinical trial with 500
patients needs 4 GB. The sizing strategy adapts:

| Phase | Approach |
|-------|----------|
| **Template defaults** | Each analysis type has baseline compute recommendations (small/medium/large profile) based on data size heuristics |
| **Pilot** | Run the first 2–3 modules on the default profile to measure actual resource usage |
| **Sizing report** | Orchestrator projects peak requirements from pilot metrics |
| **Provisioning** | Provision at 1.5x projected peak (headroom prevents OOM) |
| **Monitoring** | Orchestrator tracks utilization and recommends right-sizing |

Analysis types with predictable compute (clinical trials, most time
series) can skip the pilot and use template defaults directly. Analysis
types with data-dependent compute (scRNA-seq, proteomics with large
cohorts) should always run the pilot.

### Cloud-native scaling

For cloud deployments, the orchestrator can:
- Start with a smaller instance for lightweight modules (data import,
  metadata)
- Scale up for compute-heavy modules (integration, model fitting)
- Scale back down for reporting
- Use spot/preemptible instances (with checkpointing) for cost efficiency

### Multi-language execution

When an analysis requires multiple languages (R + Python, SAS + R,
Python + Julia), modules in different languages execute as independent
steps that communicate through files — not through in-process bridges
like rpy2 or reticulate. This avoids fork-safety issues, memory sharing
problems, and API incompatibilities.

```
Language A module → writes .tsv / standard format → Language B module reads
```

Cross-language data exchange uses TSV/CSV for tabular data and
domain-standard formats (h5ad, .rds, parquet) where both languages have
reliable readers. The analysis type template specifies which interchange
formats are tested.

---

## 15. Validation Architecture

Three tiers of validation, each catching a different class of error.

### Tier 1: Format validation (automated, blocks progression)

- Output files exist and are non-empty
- Expected columns present in output tables
- Cell/gene counts within plausible range
- No NaN/Inf in critical fields
- File checksums match expected formats

These are the checks that the IVD project already implemented. They catch
coding errors and data corruption. They do not catch scientific errors.

### Tier 2: Domain-specific sanity checks (automated, reported at checkpoint)

Tier 2 checks are defined per analysis type. They encode domain
knowledge that can be tested computationally — things the SME would
immediately notice if they looked, but that should be flagged
automatically so the checkpoint review is focused.

**scRNA-seq examples:**
- Known marker genes expressed in expected clusters (COL2A1 in
  chondrocytes, PTPRC in immune cells)
- Markers absent from unexpected clusters
- Cell type proportions within literature-reported ranges

**Proteomics examples:**
- Known abundant proteins (albumin, actin) present in expected
  fractions
- Coefficient of variation within replicates below threshold
- No systematic intensity bias between batches after normalization

**Clinical trial examples:**
- Baseline characteristics balanced between arms (standardized
  differences < 0.1)
- Randomization ratio matches protocol
- Primary endpoint analysis matches pre-specified method in SAP
- AE rates within plausible ranges for the indication

**Time series examples:**
- Residuals approximately white noise (Ljung-Box test)
- No data leakage (future values not used in features)
- Forecast intervals have correct coverage on holdout data

**Cross-domain checks (always active):**
- Cross-method concordance (Section 11)
- Claim verification (Section 10)

These checks would have caught 3 of the 4 major errors in the IVD
project (the misrouted stressed NP cells, the fabricated gene claims,
and the stale method references).

### Tier 3: Human review (checkpoint)

- Is the biology plausible?
- Do the results match domain knowledge?
- Are the caveats appropriate?
- Should the approach change?

This is irreplaceable. The framework optimizes the human's time by
presenting pre-validated, well-organized review materials — not raw
outputs.

---

## 16. Reporting

### Result tables (committed)

Every module produces summary tables that are committed to the repository.
The specific tables are defined by the analysis type template. Examples:

**scRNA-seq:**

| Module | Table | Contents |
|--------|-------|----------|
| QC | `qc_summary.tsv` | Per-dataset cell counts, median genes, % mito |
| Annotation | `cell_type_counts.tsv` | Cell type × dataset × condition counts |
| DE | `de_summary.tsv` | Per-comparison: n_sig, top genes, direction |
| Enrichment | `enrichment_summary.tsv` | Significant pathways per cell type |

**Clinical trial:**

| Module | Table | Contents |
|--------|-------|----------|
| Baseline | `demographics.tsv` | Table 1: baseline characteristics by arm |
| Primary | `primary_endpoint.tsv` | Effect estimate, CI, p-value, sensitivity results |
| Safety | `ae_summary.tsv` | AE rates by SOC, SAEs, discontinuations |
| Subgroups | `subgroup_summary.tsv` | Forest plot data: subgroup × effect estimate |

**Proteomics:**

| Module | Table | Contents |
|--------|-------|----------|
| QC | `qc_summary.tsv` | Identified proteins, peptides, missingness |
| DA | `da_summary.tsv` | Per-comparison: n_sig, top proteins, direction |
| Enrichment | `enrichment_summary.tsv` | Significant GO/KEGG terms |

### Final report

The report is generated from committed tables — it can be built on any
machine without access to the large result files. It contains:

1. **Methods:** Automatically generated from specs and parameters
2. **Results:** Data-driven summaries with table references (verified by
   the claim verification pipeline)
3. **Stability analysis:** Which results are robust vs. method-sensitive
4. **Supplementary tables:** All committed result tables
5. **Figures:** Generated from notebooks

The report is the SME's primary deliverable. It is complete enough to
hand to a collaborator or submit as supplementary material.

---

## 17. What This Framework Does Not Do

Being explicit about scope prevents overdesign:

- **It does not replace the SME.** The framework automates execution,
  not judgment. Every analytical decision goes through a human.
- **It does not auto-select methods.** The agent can compare methods the
  SME specifies, but it does not decide which method is "best." That is
  an analytical judgment.
- **It does not generate interpretations.** The agent produces data
  summaries. Scientific, clinical, or business interpretation is the
  SME's job.
- **It does not handle novel analysis types without a template.** New
  analysis types require a new set of module specs. The framework makes
  it easy to create templates from prior analyses (Section 7), but
  someone with domain expertise must define the module structure, default
  parameters, and validation checks for the first instance.
- **It does not guarantee correctness.** Tier 1 and Tier 2 validation
  catch many errors, but scientific correctness ultimately requires human
  judgment. The framework is designed to make that judgment efficient, not
  to replace it.
- **It does not replace regulatory review.** For regulated analyses
  (clinical trials, GxP), the framework provides provenance, audit
  trails, and reproducibility — but it does not substitute for
  regulatory-compliant validation (IQ/OQ/PQ), 21 CFR Part 11 compliance,
  or GxP-qualified environments. These would need to be layered on top.

---

## 18. Implementation Roadmap

### Phase 1: Core infrastructure (domain-agnostic)

- Structured state management (YAML state file, orchestrator loop)
- Checkpoint enforcement with notification system
- Claim verification pipeline (configurable entity types)
- Auto-commit and provenance tracking
- Analysis type plugin architecture (template registry)

### Phase 2: First analysis type + SME interface

- Complete module template set for one domain (scRNA-seq as the
  reference implementation, since it has been validated through 5
  pipeline versions)
- Container image for that domain
- Project setup wizard (CLI first, web later)
- Checkpoint review interface (rendered notebooks + decision form)
- Natural-language spec editing

### Phase 3: Second and third analysis types

- Proteomics module template set + container
- Clinical trial module template set + container (with regulatory
  provenance extensions)
- Template-from-project extraction tool (generalize a completed project
  into a reusable analysis type)
- Cross-domain validation of the core infrastructure

### Phase 4: Advanced features

- Multi-method/sensitivity analysis framework
- Dependency-aware re-execution with caching
- Compute auto-scaling and environment sizing
- Cross-version sensitivity tracking
- Multi-user support (multiple SMEs reviewing different checkpoints)
- Integration with lab notebook / ELN / CTMS systems

---

## Appendix A: Lessons from the IVD Project That Shaped This Design

| IVD failure mode | v2 design response | Section |
|-----------------|-------------------|---------|
| Agent skipped checkpoints | Orchestrator enforces checkpoints mechanically | 5, 9 |
| 14 fabricated gene claims | Claim verification pipeline; narrative/data separation | 10 |
| 4 compute configurations, all human-initiated | Environment sizing module; compute monitoring | 9, 14 |
| `analysis_plan.md` grew unwieldy | Structured YAML state; separate history | 4, 5 |
| CCC direction reversed across versions | Multi-method stability analysis | 11 |
| Annotation errors propagated to all downstream modules | Annotation as required checkpoint; dependency-aware re-execution | 6, 13 |
| R/Python version conflicts | Containerized environment; multi-language separation | 14 |
| Results on remote machine, report not generatable locally | Summary tables committed to repo | 4, 16 |
| Silent compute stalls | Orchestrator monitors CPU/memory utilization | 9 |
| SME needed operator to run pipeline | Orchestrator automates loop; SME interface for review | 9, 12 |

## Appendix B: Example `project.yaml` files

### B.1 scRNA-seq atlas

```yaml
project:
  name: "IVD Single-Cell Atlas"
  type: scrna_atlas
  description: >
    Identify cell types in human intervertebral disc and characterize
    changes with degeneration.

domain:                              # domain-specific section — schema
  species: human                     # defined by the analysis type template
  tissue: intervertebral_disc
  compartments: [NP, AF, CEP]
  conditions:
    - name: healthy
      aliases: [normal, control, Pfirrmann_1, Pfirrmann_2]
    - name: mild_degeneration
      aliases: [mild, Pfirrmann_3]
    - name: severe_degeneration
      aliases: [severe, Pfirrmann_4, Pfirrmann_5, herniated]
  sanity_checks:
    expected_cell_types:
      - chondrocyte_like: {markers: [COL2A1, ACAN, SOX9]}
      - immune: {markers: [PTPRC, CD68, CD3D]}
  known_issues:
    - "GSE242443 CEP cells are culture-expanded"

data:
  sources:
    - accession: GSE160756
      compartments: [NP, AF, CEP]
      platform: 10x
    - accession: GSE165722
      compartments: [NP]
      platform: BD_Rhapsody

pipeline:
  checkpoint_mode: selective
  checkpoints_required: [04_annotation, 05_integration, 07_annotation]
  method_comparisons:
    integration: [CCA, scANVI]

compute:
  profile: cloud_auto
  gpu_required: true
  min_ram_gb: 128
  container: "bioagent/scrna:1.0"

notifications:
  checkpoint: [email, slack]
  stall: [slack]
```

### B.2 Clinical trial evaluation

```yaml
project:
  name: "AURORA Phase III Analysis"
  type: clinical_trial
  description: >
    Primary and secondary endpoint analysis for a Phase III randomized
    controlled trial of Drug X vs. placebo in moderate-to-severe RA.

domain:
  study_design: parallel_group
  phase: 3
  indication: rheumatoid_arthritis
  arms:
    - name: treatment
      drug: "Drug X 200mg"
    - name: placebo
      drug: "Placebo"
  populations:
    primary: ITT
    sensitivity: [mITT, per_protocol]
  endpoints:
    primary:
      name: "ACR20 response at Week 24"
      type: binary
      method: logistic_regression
      covariates: [site, baseline_DAS28, prior_biologic]
    secondary:
      - name: "DAS28-CRP change from baseline"
        type: continuous
        method: MMRM
      - name: "HAQ-DI change from baseline"
        type: continuous
        method: MMRM
  multiplicity: hierarchical_testing
  missing_data: [MAR_primary, tipping_point_sensitivity]
  sanity_checks:
    - "Randomization ratio approximately 1:1"
    - "Baseline DAS28 > 5.1 (moderate-to-severe inclusion criterion)"
    - "Primary endpoint matches SAP specification"
  known_issues:
    - "Site 042 had high dropout — monitor per-protocol exclusion rate"

data:
  sources:
    - path: data/raw/adsl.sas7bdat
      type: ADaM
    - path: data/raw/adae.sas7bdat
      type: ADaM

pipeline:
  checkpoint_mode: gated            # all checkpoints required for Phase III
  method_comparisons:
    primary_endpoint: [logistic_regression, CMH_stratified]

compute:
  profile: local                    # clinical trial data is small
  min_ram_gb: 8
  container: "bioagent/clinical:1.0"

notifications:
  checkpoint: [email]
```

### B.3 Quantitative proteomics

```yaml
project:
  name: "Synovial Fluid Proteome in OA"
  type: proteomics
  description: >
    Differential protein abundance in synovial fluid between early and
    late-stage osteoarthritis using DIA-MS.

domain:
  species: human
  sample_type: synovial_fluid
  acquisition: DIA
  search_engine: DIA-NN
  conditions:
    - name: early_OA
      aliases: [KL_1, KL_2]
    - name: late_OA
      aliases: [KL_3, KL_4]
  design:
    replicates: biological           # biological vs. technical
    batches: 3
  sanity_checks:
    expected_proteins: [ALB, A2M, FGA, FGB, FGG]   # abundant in SF
    max_cv_within_replicates: 0.25

data:
  sources:
    - path: data/raw/diann_report.tsv
      format: DIA-NN_main_output
    - path: data/raw/sample_annotation.tsv
      format: metadata

pipeline:
  checkpoint_mode: selective
  checkpoints_required: [normalization, differential_abundance]
  method_comparisons:
    normalization: [median, quantile]
    differential: [limma, MSstats]

compute:
  profile: cloud_fixed
  min_ram_gb: 32
  container: "bioagent/proteomics:1.0"
```

Note the common structure across all three examples: `project` (universal
metadata), `domain` (analysis-type-specific parameters including sanity
checks), `data` (input sources), `pipeline` (checkpoint and comparison
config), `compute` (resources), `notifications`. The `domain` section's
schema is defined by the analysis type template — it contains whatever
that domain needs. The rest is framework-standard.
