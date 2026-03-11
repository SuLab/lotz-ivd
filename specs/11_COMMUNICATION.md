# Module 11: Cell-Cell Communication Analysis

## Objective

Infer ligand-receptor interactions between IVD cell populations and determine how these change with degeneration. Of particular interest: immune cell signaling to resident IVD cells, and whether degenerated cells produce signals that promote nerve ingrowth and pain.

## Inputs

- Integrated all-cells AnnData object: `data/integrated/all_cells.h5ad` — with de novo cell type labels (`obs['cell_type']`) from Module 07
- DE results from Module 08

**Scope:** Run cell-cell communication on the **all-cells object only**. Compartment-specific objects (NP, AF, CEP) may have very few non-mesenchymal cells, making interaction inference unreliable.

## Outputs

- `results/communication/interactions_{comparison}.tsv`
- `results/communication/interaction_plots/`
- `results/communication/communication_report.html`

### Notebook: `notebooks/11_communication.ipynb`

Produced after cell-cell communication analysis. Contains:
- Chord diagram: number/strength of interactions between cell type pairs, per condition
- Dot plot: top significant ligand-receptor pairs, source → target, per condition
- Differential interaction heatmap: gained/lost interactions in degeneration
- Focused panels:
  - Immune → NP signaling (inflammatory axis)
  - NP autocrine/paracrine signaling changes
  - Pain-relevant interactions: neurotrophin, semaphorin, VEGF pathways highlighted
- Network visualization: cell types as nodes, interaction strength as edges, comparing healthy vs. degenerated

**Manuscript mapping:** Figure 7: Cell-cell communication in healthy vs. degenerated IVD. Supplementary: full interaction tables.

## Method

### Primary tool: LIANA+

LIANA aggregates multiple ligand-receptor inference methods (CellPhoneDB, NATMI, Connectome, etc.) and provides a consensus score. This is more robust than any single method.

### CellChat

CellChat is not run. LIANA already aggregates multiple methods internally, making a second tool redundant.

### Steps

1. For each condition (or comparison):
   a. Subset the annotated AnnData to the relevant samples
   b. Run LIANA with default methods (CellPhoneDB, NATMI, SingleCellSignalR at minimum)
   c. Filter for significant interactions using LIANA's aggregate rank p-value < 0.05
   d. Record: source cell type, target cell type, ligand, receptor, significance score
2. Compare interactions between conditions using LIANA+'s built-in differential analysis (`liana.multi.test`):
   a. Which interactions are gained or lost in degeneration?
   b. Which interactions change in magnitude?
3. Focus analyses:
   a. **Immune → resident:** What are macrophages/T cells signaling to NP/AF cells? (e.g., IL1B-IL1R1, TNF-TNFR)
   b. **Resident → resident:** Do healthy NP cells communicate differently than degenerated NP cells?
   c. **Pain-relevant interactions:** Are there interactions that promote neurotrophin signaling (NGF-NTRK1), nerve growth (SEMA, NTN pathways), or neovascularization (VEGF signaling)?

### Minimum requirements

- Cell-cell communication requires multiple cell types to be present in the same dataset/sample
- Datasets that only contain sorted NP cells (no immune/endothelial) cannot be used for this analysis
- Minimum 50 cells per cell type per sample for reliable interaction inference

## Automated Validation

- [ ] Interaction results exist for at least one condition
- [ ] Results include known positive interactions (e.g., collagen-integrin interactions between matrix-producing cells and their neighbors)
- [ ] Pain-relevant interactions are specifically flagged
- [ ] Visualization plots are generated (chord diagrams, dot plots)
- [ ] Communication report is generated

## Human Checkpoint

### Questions for the reviewer
1. Are the top interactions biologically plausible?
2. Do condition-specific interactions suggest mechanisms by which immune cells drive degeneration?
3. Are there pain-relevant interactions that could be therapeutic targets?
4. Are any interactions likely artifacts of ambient RNA or doublets?

### Potential plan revisions
- If immune-resident interactions are a major finding, expand immune subtype analysis
- If pain-relevant interactions are identified, prioritize these in the final report
