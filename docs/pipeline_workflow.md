# IVD Single-Cell Atlas: Pipeline Workflow

## Overview

Human-gated agentic bioinformatics pipeline analyzing 12 scRNA-seq datasets (436K cells, 78 samples, 57 donors) of human intervertebral disc tissue. Each module produces results that are reviewed at a human checkpoint before the pipeline advances.

## Workflow Diagram

```mermaid
flowchart TD
    %% ── Styles ──
    classDef module fill:#4A90D9,stroke:#2C5F8A,color:#fff,stroke-width:2px
    classDef checkpoint fill:#E8A838,stroke:#B07D20,color:#fff,stroke-width:3px
    classDef decision fill:#FFFACD,stroke:#B07D20,color:#333,stroke-width:1px,font-size:11px
    classDef data fill:#7BC47F,stroke:#4A8A4E,color:#fff,stroke-width:2px

    %% ── Spec Writing ──
    SPECS["Spec Writing<br/>10 module specifications"]:::module
    SPECS --> CP_SPECS

    CP_SPECS{{"HUMAN CHECKPOINT<br/>Spec Review & Approval"}}:::checkpoint
    CP_SPECS --> M01

    %% ── Module 01 ──
    M01["Module 01: Dataset Discovery<br/>13 candidates across 7 databases<br/>~533K cells total"]:::module
    M01 --> CP01

    CP01{{"HUMAN CHECKPOINT<br/>Approve Dataset List"}}:::checkpoint
    N01["Decisions:<br/>- Include GSE242443 (culture-expanded CEP)<br/>- Defer Zhou 2023 embryonic data to Module 08<br/>- Drop 2 NGDC datasets (NP well-covered)<br/>- 12 datasets, coverage adequate"]:::decision
    CP01 -.- N01
    CP01 --> M02

    %% ── Module 02 ──
    M02["Module 02: Metadata Harmonization<br/>78 samples, 57 donors, 12 studies<br/>Condition mapping, demographics, compartment labels"]:::module
    M02 --> CP02

    CP02{{"HUMAN CHECKPOINT<br/>Approve Condition Mappings"}}:::checkpoint
    N02["Decisions (tentative):<br/>- Herniated kept separate from degenerated<br/>- GSE165722 Pfirrmann grades corrected (GEO off-by-one)<br/>- GSE244889 Pfirrmann I reclassified as healthy<br/>- Thompson III boundary: degenerated_mild<br/>- MUST revisit mappings before Module 06"]:::decision
    CP02 -.- N02
    CP02 --> M03

    %% ── Module 03 ──
    M03["Module 03: Preprocessing<br/>436,558 cells post-QC across 12 datasets<br/>QC: min_genes=200, max_genes=6000, max_mt=20%<br/>Scrublet doublet detection"]:::module
    M03 --> CP03

    CP03{{"HUMAN CHECKPOINT<br/>QC Review (retroactive)"}}:::checkpoint
    N03["Notes:<br/>- 4 datasets had 100% retention (pre-filtered)<br/>- GSE251686_NP3 excluded (corrupt matrix)<br/>- Diffuse CD68 expected IVD biology<br/>- No blocking issues"]:::decision
    CP03 -.- N03
    CP03 --> M04

    %% ── Module 04 ──
    M04["Module 04: Cell Type Annotation<br/>Marker-based scoring (16 signatures)<br/>+ CellTypist (Immune_All_Low)<br/>Consensus labels in cell_type_final"]:::module
    M04 --> CP04

    CP04{{"HUMAN CHECKPOINT<br/>Annotation Review (retroactive)"}}:::checkpoint
    N04["Notes:<br/>- NP subtypes: notochordal, mature chondrocyte, stressed, fibrocartilaginous<br/>- AF subtypes: inner, outer, mechanical stress<br/>- CellTypist refined immune populations<br/>- No IVD reference atlas available"]:::decision
    CP04 -.- N04
    CP04 --> M05

    %% ── Module 05 ──
    M05["Module 05: Integration<br/>Tier 1: Non-resident cells (14.6K cells, scVI)<br/>Tier 2: Resident cells (NP 139K, AF 283K)<br/>4 approaches: scVI, scANVI, Harmony, BBKNN"]:::module
    M05 --> CP05a

    CP05a{{"HUMAN CHECKPOINT<br/>Tier 1 Integration Review"}}:::checkpoint
    N05a["Decisions:<br/>- Tier 1 scVI integration approved<br/>- Retroactive approval of Modules 03-04<br/>- Proceed to Tier 2"]:::decision
    CP05a -.- N05a
    CP05a --> CP05b

    CP05b{{"HUMAN CHECKPOINT<br/>Tier 2 Integration Selection"}}:::checkpoint
    N05b["Decisions:<br/>- Primary: scANVI (best overall + cell type ASW)<br/>- Sensitivity: scVI (preserves cell state continuum)<br/>- Harmony rejected (overcorrects, merges clusters)<br/>- BBKNN not primary (no corrected embedding)<br/>- Pseudobulk DE uses scANVI labels + raw counts"]:::decision
    CP05b -.- N05b

    %% ── Condition Mapping Revisit ──
    CP05b --> COND_REVIEW

    COND_REVIEW{{"HUMAN CHECKPOINT<br/>Condition Mapping Revisit<br/>(required before DE analysis)"}}:::checkpoint
    N_COND["Decisions:<br/>- Herniated: separate category (10 samples), exploratory<br/>- GSE205535_NNP (11yo trauma): exclude from DE<br/>- Thompson III boundary: accepted as mild<br/>- Neonatal (n=3): separate, not mixed into healthy<br/>- Aged ungraded (n=3): separate category<br/>- Primary: healthy (20) vs degenerated_all (42)"]:::decision
    COND_REVIEW -.- N_COND
    COND_REVIEW --> M06

    %% ── Module 06 ──
    M06["Module 06: Differential Analysis<br/>Composition: Mann-Whitney U (0/58 significant)<br/>Pseudobulk DE: pyDESeq2, 17 powered comparisons<br/>5,328 significant genes"]:::module
    M06 --> CP06

    CP06{{"HUMAN CHECKPOINT<br/>DE Results Review"}}:::checkpoint
    N06["Decisions:<br/>- Herniated comparison exploratory only (RPL genes in top hits, likely study-confounded)<br/>- Endothelial annotation caveat (some may be misclassified NP/AF)<br/>- Composition trends sensible despite no FDR hits<br/>- No systematic batch domination"]:::decision
    CP06 -.- N06
    CP06 --> M07

    %% ── Module 07 ──
    M07["Module 07: Biological Interpretation<br/>ORA: 1,244 enrichments (GO/KEGG/Reactome/MSigDB)<br/>GSEA: 1,081 significant terms<br/>TF activity: 113 significant TFs<br/>Pain gene analysis: 3 significant hits"]:::module
    M07 --> CP07

    CP07{{"HUMAN CHECKPOINT<br/>Interpretation Review"}}:::checkpoint
    N07["Decisions:<br/>- Pathways consistent with known IVD biology<br/>- Novel TF findings (ATF3/7, HSF1/2) to highlight<br/>- Pain analysis confirms indirect signaling model<br/>- No contradictions found"]:::decision
    CP07 -.- N07
    CP07 --> M08

    %% ── Module 08 ──
    M08["Module 08: Trajectory Analysis<br/>PAGA + diffusion pseudotime<br/>NP root: notochordal cluster<br/>AF root: AF_inner cluster<br/>500 trajectory genes per compartment"]:::module
    M08 --> CP08

    CP08{{"HUMAN CHECKPOINT<br/>Trajectory Review"}}:::checkpoint
    N08["Decisions:<br/>- NP: notochordal -> mature -> stressed<br/>- AF: inner -> outer -> mechanical_stress<br/>- Pseudotime correlates with disease condition<br/>- ~55% DE gene overlap confirms consistency<br/>- RNA velocity absent (acceptable)<br/>- Sensitivity check (scVI) confirms direction"]:::decision
    CP08 -.- N08
    CP08 --> M09

    %% ── Module 09 ──
    M09["Module 09: Cell-Cell Communication<br/>LIANA (4 methods consensus)<br/>Healthy: 44K interactions | Degenerated: 53K<br/>3,662-4,194 pain-relevant interactions"]:::module
    M09 --> CP09

    CP09{{"HUMAN CHECKPOINT<br/>Communication Review"}}:::checkpoint
    N09["Decisions:<br/>- Interactions biologically plausible<br/>- Pain-relevant: neurotrophin + VEGF pathways<br/>- More interactions in degeneration (53K vs 44K)<br/>- Collagen-integrin positive controls confirmed"]:::decision
    CP09 -.- N09
    CP09 --> M10

    %% ── Module 10 ──
    M10["Module 10: Final Reporting<br/>12-section report, 13 supplementary tables<br/>All 18 validation checks PASS"]:::module
    M10 --> CP10

    CP10{{"HUMAN CHECKPOINT<br/>Final Review"}}:::checkpoint
    N10["AWAITING REVIEW"]:::decision
    CP10 -.- N10

    %% ── Data flow annotations ──
    D1[("data/raw/<br/>12 datasets")]:::data
    D2[("data/processed/<br/>per-dataset .h5ad")]:::data
    D3[("data/integrated/<br/>NP + AF atlases")]:::data
    D4[("results/<br/>DE, enrichments, trajectory, CCC")]:::data
    D5[("results/final_report.html")]:::data

    D1 -.-> M03
    M03 -.-> D2
    D2 -.-> M05
    M05 -.-> D3
    D3 -.-> M06
    M06 -.-> D4
    D4 -.-> M10
    M10 -.-> D5
```

## Legend

| Element | Meaning |
|---------|---------|
| Blue boxes | Computational modules (agent-executed) |
| Orange hexagons | Human checkpoints (pipeline pauses for expert review) |
| Yellow notes | Key decisions made at each checkpoint |
| Green cylinders | Data artifacts |
| Dashed arrows | Data flow |
| Solid arrows | Pipeline sequence |

## Key Pipeline Characteristics

- **Human-in-the-loop**: 12 human checkpoints across 10 modules ensure scientific rigor
- **Tiered integration**: Non-resident cells (immune, endothelial) integrated separately from resident cells (NP, AF) to preserve the IVD cell state continuum
- **Pseudobulk DE**: Avoids treating cells as independent observations (a common single-cell pitfall)
- **Multi-method validation**: Integration compared 4 approaches; CCC used 4 method consensus
- **Condition mapping revisited**: Explicit checkpoint before DE analysis to finalize disease categories
