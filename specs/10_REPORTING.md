# Module 10: Reporting & Reproducibility

## Objective

Produce a comprehensive, reproducible record of the entire analysis: methods, parameters, results, and interpretations. Generate publication-quality figures and a final integrated report.

## Outputs

- `results/final_report.html` — comprehensive analysis report
- `results/figures/` — publication-quality figures
- `results/supplementary_tables/` — all result tables in a single organized directory
- `analysis_plan.md` — finalized with complete revision history
- Environment specification (`environment.yml` or `requirements.txt`)
- Script archive with all analysis code

## Report Structure

1. **Overview:** Project goal, datasets included, total cells analyzed
2. **Dataset summary:** Table of all datasets with key characteristics, QC outcomes
3. **Cell type atlas:** Final cell type definitions, marker genes, UMAP visualizations
4. **Integration:** Strategy chosen and rationale, quality metrics
5. **Composition changes:** Cell type proportion shifts with disease/aging
6. **Differential expression:** Key DE findings per cell type, volcano plots, heatmaps
7. **Biological pathways:** Enriched pathways and gene programs
8. **Regulatory networks:** Key TFs and regulons
9. **Trajectories:** Cell state transitions and gene dynamics
10. **Cell communication:** Intercellular signaling changes
11. **Pain biology:** Synthesis of pain-relevant findings
12. **Limitations:** Known caveats, confounders, underpowered comparisons
13. **Methods:** Complete methods section suitable for a manuscript

## Reproducibility Requirements

- [ ] All scripts are version-controlled and can recreate the analysis from raw data
- [ ] All random seeds are recorded
- [ ] Package versions are pinned
- [ ] All parameter choices are documented with rationale
- [ ] All human checkpoint decisions are recorded in analysis_plan.md
- [ ] Data provenance: download dates, checksums, accessions all recorded

## Human Checkpoint (Final)

### Questions for the reviewer
1. Does the report accurately represent the findings?
2. Are the conclusions supported by the evidence?
3. Are limitations adequately described?
4. What follow-up analyses or experiments are suggested by the results?
5. Is this ready for presentation to collaborators or for manuscript preparation?
