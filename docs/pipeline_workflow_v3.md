# IVD Single-Cell Atlas: Pipeline Workflow (v3)

## Overview

Human-gated agentic bioinformatics pipeline analyzing 11 scRNA-seq datasets (410,759 cells, 71 samples, ~50 donors) of human intervertebral disc tissue. Each module produces results reviewed at a human checkpoint before the pipeline advances.

**Pipeline version:** v3 (annotation fix, 2026-03-10). Commit: `6622221`.

## What Changed from v2

- **Three annotation fixes** in Module 04 to address the ~17K misrouted stressed NP cells:
  1. Non-mesenchymal evidence gate (requires canonical immune/endothelial markers)
  2. ACAN/SOX9 rescue (25,415 cells reclassified to mesenchymal)
  3. Stricter 85% cluster voting threshold (up from 70%)
- **Specificity-weighted scoring** in Module 05 de novo annotation
- **Minimum non-mesenchymal resolution** floor of 0.5 in Module 05
- **Modules 04-10 rerun**; Modules 01-03 unchanged

## Workflow Diagram

```mermaid
flowchart TD
    classDef module fill:#4A90D9,stroke:#2C5F8A,color:#fff,stroke-width:2px
    classDef checkpoint fill:#E8A838,stroke:#B07D20,color:#fff,stroke-width:3px
    classDef decision fill:#FFFACD,stroke:#B07D20,color:#333,stroke-width:1px,font-size:11px
    classDef data fill:#7BC47F,stroke:#4A8A4E,color:#fff,stroke-width:2px

    SPECS["Spec Writing<br/>10 module specifications"]:::module
    SPECS --> CP_SPECS
    CP_SPECS{{"HUMAN CHECKPOINT<br/>Spec Review & Approval"}}:::checkpoint
    CP_SPECS --> M01

    M01["Module 01: Dataset Discovery<br/>11 datasets selected<br/>410,759 cells, 71 samples"]:::module
    M01 --> CP01
    CP01{{"HUMAN CHECKPOINT<br/>Approve Dataset List"}}:::checkpoint
    N01["11 datasets (GSE233666 excluded)<br/>Same as v2"]:::decision
    CP01 -.- N01
    CP01 --> M02

    M02["Module 02: Metadata Harmonization<br/>71 samples, ~50 donors, 11 studies"]:::module
    M02 --> CP02
    CP02{{"HUMAN CHECKPOINT<br/>Approve Condition Mappings"}}:::checkpoint
    N02["Same mappings as v2<br/>Herniated excluded from DE"]:::decision
    CP02 -.- N02
    CP02 --> M03

    M03["Module 03: Preprocessing<br/>410,759 cells post-QC<br/>Unchanged from v2"]:::module
    M03 --> CP03
    CP03{{"HUMAN CHECKPOINT<br/>QC Review"}}:::checkpoint
    CP03 --> M04

    M04["Module 04: Binary Classification (v3 fix)<br/>Mesenchymal vs non-mesenchymal<br/>+ Evidence gate + ACAN/SOX9 rescue + 85% voting<br/>25,415 cells rescued to mesenchymal"]:::module
    M04 --> CP04
    CP04{{"HUMAN CHECKPOINT<br/>Classification Review"}}:::checkpoint
    N04["v3 annotation fix applied:<br/>1. Non-mes evidence gate (PTPRC/PECAM1/CD68)<br/>2. ACAN/SOX9 rescue (25,415 cells)<br/>3. 85% cluster voting (up from 70%)<br/>Fixed 17K misrouted stressed NP cells"]:::decision
    CP04 -.- N04
    CP04 --> M05

    M05["Module 05: Integration + Clustering + Annotation<br/>scVI-only, 4 compartments:<br/>NP (263K), AF (84.6K), CEP (50.7K), all_cells (411K)<br/>Specificity-weighted scoring, min resolution 0.5"]:::module
    M05 --> CP05
    CP05{{"HUMAN CHECKPOINT<br/>Integration + Annotation Review"}}:::checkpoint
    N05["MOST CRITICAL GATE<br/>~10 cell types resolved<br/>AF_mechanical_stress re-emerged (was in v1, lost in v2)<br/>EP_ossification emerged (new CEP type)<br/>NP_fibrocartilaginous now top DE responder"]:::decision
    CP05 -.- N05
    CP05 --> M06

    M06["Module 06: Differential Analysis<br/>Pseudobulk DE: pyDESeq2<br/>18 powered comparisons<br/>1,156 unique sig genes, 1,447 pairs"]:::module
    M06 --> CP06
    CP06{{"HUMAN CHECKPOINT<br/>DE Results Review"}}:::checkpoint
    N06["18 powered comparisons<br/>1,156 unique genes (+22% from v2)<br/>CXCL2 log2FC=3.63, padj=1.75e-4 (strongest)<br/>NP_fibrocartilaginous m_vs_s: 418 genes (top)"]:::decision
    CP06 -.- N06
    CP06 --> M07

    M07["Module 07: Biological Interpretation<br/>ORA: 1,043 enrichments<br/>5 significant TF associations (down from 290)<br/>10 significant pain genes"]:::module
    M07 --> CP07
    CP07{{"HUMAN CHECKPOINT<br/>Interpretation Review"}}:::checkpoint
    N07["TF associations collapsed (290 → 5)<br/>Prostaglandin pathway preserved<br/>New pain genes: NRP2, ROBO1, SEMA3A (axon guidance)<br/>TF collapse likely artifact of annotation change"]:::decision
    CP07 -.- N07
    CP07 --> M08

    M08["Module 08: Trajectory Analysis<br/>PAGA + diffusion pseudotime<br/>NP, AF, and CEP compartments"]:::module
    M08 --> CP08
    CP08{{"HUMAN CHECKPOINT<br/>Trajectory Review"}}:::checkpoint
    N08["NP rho=-0.151 (consistent, weaker)<br/>AF rho=+0.325 (consistent with v2)<br/>CEP rho=+0.135 (REVERSED from v2)<br/>CEP instability flagged"]:::decision
    CP08 -.- N08
    CP08 --> M09

    M09["Module 09: Cell-Cell Communication<br/>LIANA (5 methods consensus)<br/>Healthy: 40K | Degenerated: 41K<br/>Roughly balanced"]:::module
    M09 --> CP09
    CP09{{"HUMAN CHECKPOINT<br/>Communication Review"}}:::checkpoint
    N09["Near-balanced (40K vs 41K)<br/>Third directional result in 3 versions<br/>CCC aggregate counts deemed unreliable"]:::decision
    CP09 -.- N09
    CP09 --> M10

    M10["Module 10: Final Reporting<br/>Full report + supplementary tables"]:::module
    M10 --> CP10
    CP10{{"HUMAN CHECKPOINT<br/>Final Review"}}:::checkpoint

    D1[("data/raw/<br/>11 datasets")]:::data
    D2[("data/processed/<br/>per-dataset .h5ad")]:::data
    D3[("data/integrated/<br/>NP, AF, CEP, all_cells")]:::data
    D4[("results/<br/>DE, enrichments, trajectory, CCC")]:::data

    D1 -.-> M03
    M03 -.-> D2
    D2 -.-> M05
    M05 -.-> D3
    D3 -.-> M06
    M06 -.-> D4
    D4 -.-> M10
```

## Key v3 Characteristics

- **Same 11 datasets and 410,759 cells** as v2
- **Annotation fix** rescued 25,415 misrouted cells back to mesenchymal
- **~10 cell types** including AF_mechanical_stress (re-emerged) and EP_ossification (new)
- **scVI-only integration** (unchanged from v2)
- **Strongest CXCL2 result** across v1-v3 (padj=1.75e-4)
- **TF activity collapsed** from 290 to 5 significant associations

## v3 Key Results

| Metric | Value |
|--------|-------|
| Datasets | 11 |
| Total cells | 410,759 |
| Cell types | ~10 |
| Powered DE comparisons | 18 |
| Unique DE genes | 1,156 |
| Top DE comparison | NP_fibrocartilaginous m_vs_s (418 genes) |
| CXCL2 | log2FC=3.63, padj=1.75e-4 |
| NP trajectory rho | -0.151 |
| AF trajectory rho | +0.325 |
| CEP trajectory rho | +0.135 (reversed from v2) |
| CCC healthy/degen | 40K / 41K (balanced) |
| TF associations | 5 |
| Pain genes | 10 |

## What the Annotation Fix Revealed

1. **NP_fibrocartilaginous became the top DE responder** — 418 genes in mild_vs_severe (doubled from v2's 203)
2. **CXCL2 signal strengthened** with purer pseudobulk profiles
3. **TF associations collapsed** (290 → 5) — v2's 290 were likely inflated by misannotated cells
4. **Axon guidance pain genes emerged** (NRP2, ROBO1, SEMA3A) — later did not replicate in v4
5. **Annotation is the single largest source of downstream variation** — changing annotation alone altered DE counts by 22%, TFs by 98%, CCC by ~45%

---

*Superseded by v4 (scANVI semi-supervised, 12-module pipeline, two-stage annotation). See `results/VERSION_DIFFERENCES_SUMMARY.md` for cross-version comparison.*
