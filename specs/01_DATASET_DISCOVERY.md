# Module 01: Dataset Discovery & Acquisition

## Objective

Systematically identify, evaluate, and acquire all publicly available human IVD single-cell RNA-seq datasets suitable for the project goal. This includes the 8 studies already identified and any additional datasets that may have been missed or recently published.

## Rationale

The known dataset list (8 studies, ~53 samples) was compiled manually and may be incomplete. New studies are published regularly. A systematic search ensures comprehensive coverage and makes the inclusion criteria explicit and reproducible.

## Inputs

- Known dataset list (see Known Datasets section below)
- Inclusion/exclusion criteria (defined in this spec)

## Outputs

- `metadata/dataset_registry.tsv` — master list of all candidate datasets with accession, status (included/excluded), and reason
- `metadata/search_log.md` — record of all searches performed, with dates, queries, and result counts
- Raw data files downloaded to `data/raw/{study_accession}/`
- A summary report for human review

## Inclusion Criteria

A dataset is included if ALL of the following are true:

1. **Species:** Human (Homo sapiens)
2. **Tissue:** Intervertebral disc tissue or cells derived from IVD (NP, AF, CEP, or whole IVD). Adjacent tissues (e.g., facet joint, ligamentum flavum, paraspinal muscle) are excluded unless they include IVD compartments.
3. **Technology:** Single-cell or single-nucleus RNA-seq (10x Genomics, Drop-seq, Smart-seq2, or similar). Spatial transcriptomics datasets should be flagged separately but not excluded — they may be useful for validation.
4. **Data availability:** Raw count matrix (or fastq files that can be processed to count matrix) is publicly accessible. Datasets described in publications but without deposited data are logged but excluded.
5. **Minimum cells:** At least 200 cells per sample after standard QC. Datasets with very low cell counts may introduce noise without adding power.

## Exclusion Criteria

A dataset is excluded if ANY of the following are true:

1. Non-human species (mouse, rat, bovine IVD studies are logged for reference but excluded from the main analysis)
2. Bulk RNA-seq, microarray, or proteomics only
3. Data is not publicly available or requires special access that cannot be obtained
4. The study is a reanalysis of an already-included dataset (avoid double-counting cells). If a reanalysis adds new samples, include only the new samples.
5. Solely in vitro cultured cells with no primary tissue comparison (culture conditions alter cell states substantially)

## Search Strategy

### Databases to search

1. **GEO (Gene Expression Omnibus)** — primary repository for most scRNA-seq studies
   - Search URL: https://www.ncbi.nlm.nih.gov/geo/
   - Also search via GEO DataSets and GEO Profiles
2. **ArrayExpress / BioStudies** — European equivalent
   - Search URL: https://www.ebi.ac.uk/biostudies/arrayexpress
3. **CellxGene** (Chan Zuckerberg Initiative) — curated single-cell data portal
   - Search URL: https://cellxgene.cziscience.com/
4. **Human Cell Atlas Data Portal**
   - Search URL: https://data.humancellatlas.org/
5. **Single Cell Portal** (Broad Institute)
   - Search URL: https://singlecell.broadinstitute.org/single_cell
6. **PubMed / Google Scholar** — to find publications that may reference deposited data not yet indexed in the above
7. **bioRxiv / medRxiv** — preprints may have associated data

### Search queries

Use combinations of the following terms. Adapt syntax to each database.

Primary terms:
- "intervertebral disc" OR "intervertebral disk"
- "nucleus pulposus"
- "annulus fibrosus" OR "anulus fibrosus"
- "endplate" AND ("spine" OR "spinal" OR "vertebral")

Combined with:
- "single-cell" OR "single cell" OR "scRNA-seq" OR "scRNA" OR "single-nucleus" OR "snRNA-seq"

Additional queries to catch edge cases:
- "disc degeneration" AND "single-cell"
- "IVD" AND "single-cell RNA"
- "notochordal" AND "single-cell" (notochordal cells are NP precursors)

### Search procedure

1. Execute all queries across all databases
2. Deduplicate results by accession number and publication DOI
3. For each candidate, record: accession, first author, year, title, journal, DOI, compartment(s), conditions, species, technology, number of samples, number of cells (if reported), data availability status
4. Apply inclusion/exclusion criteria
5. For any borderline cases, flag for human review with a note explaining the ambiguity

### Cross-referencing

After the primary search, cross-reference:
- Check the references and citations of all included publications for additional datasets
- Check review articles on IVD single-cell studies for any missed datasets
- Check CellxGene and HCA for curated IVD collections that may aggregate multiple studies

## Known Datasets

These have been previously identified and should be confirmed during the search:

| # | First Author | Year | Accession | Compartment | Conditions |
|---|-------------|------|-----------|-------------|------------|
| 1 | Jia S | 2024 | GSE251686 | NP | mild/severe degeneration |
| 2 | Cherif H | 2022 | GSE199866 | NP & inner AF | non-degenerate/degenerate |
| 3 | Jiang W | 2022 | GSE189916 | IVD | neonatal/adult |
| 4 | Li Z | 2022 | GSE205535 | NP | normal/degenerative |
| 5 | Han S | 2022 | PMID35265617 | NP | normal/mild/severe IVDD |
| 6 | Guo S | 2023 | GSE233666 | NP | IDD diagnosis |
| 7 | Shi C | 2024 | GSE255768 | Endplate | degenerative modic changes |
| 8 | Swahn H | 2024 | PMID38403470 | NP, AF | young healthy/aged diseased |

## Data Acquisition

For each included dataset:

1. Download the count matrix in the most raw available form (prefer raw counts over normalized data)
2. Download any provided metadata (sample annotations, cell barcodes, cluster labels from original study)
3. Record the exact download source, date, and file checksums (md5)
4. Organize as `data/raw/{accession}/` with a README noting the source and any preprocessing already applied by the original authors
5. If only fastq files are available, log this — alignment and quantification will be handled in the preprocessing module

### File format expectations

- 10x Genomics: expect matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz (or genes.tsv.gz), or .h5 files
- Some studies provide pre-made Seurat objects (.rds) or AnnData (.h5ad) — these are acceptable but we should verify they contain raw counts
- If only normalized/scaled data is available, log this as a limitation

## Automated Validation

These checks must pass before advancing:

- [ ] All search queries have been executed and results logged in `metadata/search_log.md`
- [ ] `metadata/dataset_registry.tsv` exists and contains all candidates with inclusion/exclusion status and reason
- [ ] All included datasets have raw data successfully downloaded to `data/raw/`
- [ ] File checksums are recorded
- [ ] No included dataset is a duplicate of another (same cells counted twice)
- [ ] Each known dataset from the table above is either included or has an explicit exclusion reason
- [ ] A summary report is generated listing: total candidates found, included, excluded (with breakdown by reason), total expected cells, coverage by compartment, coverage by condition

## Human Checkpoint

### Review materials
- The summary report
- `metadata/dataset_registry.tsv` — full candidate list
- `metadata/search_log.md` — to verify search was comprehensive

### Questions for the reviewer
1. Are the inclusion/exclusion criteria appropriate? Should any be adjusted?
2. Are there any borderline datasets that should be reconsidered?
3. Is the coverage across compartments (NP, AF, CEP) and conditions (healthy, degenerated, aged, neonatal) adequate for the project goals? If not, are there known datasets that the search missed?
4. Are there any datasets where the condition labels are ambiguous and need clarification before proceeding?
5. Are there spatial transcriptomics datasets that should be earmarked for later validation?

### Potential plan revisions triggered by this checkpoint
- If very few healthy/normal samples are found, the comparative analysis strategy may need adjustment (e.g., using "mild" as the reference instead of "normal")
- If endplate or AF data is sparse, compartment-specific analysis may not be feasible for all compartments
- If new datasets are large enough to substantially change the analysis design, downstream module parameters may need updating
- If any datasets lack raw counts, a decision is needed on whether to include normalized data with appropriate caveats
