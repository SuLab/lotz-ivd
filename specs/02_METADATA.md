# Module 02: Metadata Harmonization

## Objective

Create a unified, standardized metadata schema across all included datasets so that downstream analysis can stratify, compare, and correct for biological and technical variables consistently.

## Rationale

The existing studies use inconsistent terminology for conditions (e.g., "non-degenerate" vs. "normal" vs. "healthy" vs. "Pfirrmann I"), compartments, and grading systems. Without harmonization, any cross-study comparison is unreliable. This step is deceptively important — errors here propagate silently through every downstream analysis.

## Inputs

- `metadata/dataset_registry.tsv` (from Module 01)
- Original metadata files from each study in `data/raw/{accession}/`
- Published papers for each study (for condition definitions, grading criteria, donor demographics)

## Outputs

- `metadata/sample_metadata.tsv` — one row per sample, harmonized columns
- `metadata/harmonization_notes.md` — documenting all mapping decisions and ambiguities
- `metadata/ontology_mappings.tsv` — mapping original terms to standardized terms

## Standardized Schema

Each sample in `sample_metadata.tsv` must have the following fields:

| Field | Description | Allowed values / format |
|-------|-------------|------------------------|
| sample_id | Unique identifier (study_accession + original sample ID) | string |
| study_accession | GEO/ArrayExpress accession | string |
| first_author | First author of the publication | string |
| year | Publication year | integer |
| donor_id | Donor identifier (unique within study; flag if same donor contributes multiple samples) | string |
| n_cells_raw | Number of cells before QC | integer |
| compartment | IVD compartment | NP, AF, CEP, IVD_mixed |
| condition_original | Condition label as reported in the study | string (verbatim) |
| condition_harmonized | Harmonized condition category | see Condition Harmonization below |
| degeneration_grade_original | Grading as reported | string (verbatim) |
| degeneration_grade_system | Which grading system was used | Pfirrmann, Thompson, other, none |
| degeneration_severity | Standardized severity | none, mild, moderate, severe, ungraded |
| age_group | Age category | neonatal, young_adult, middle_aged, aged, unknown |
| age_years | Donor age in years if reported | numeric or NA |
| sex | Donor sex if reported | M, F, unknown |
| tissue_or_cells | Whether the input was fresh tissue or isolated/cultured cells | tissue, cells, unknown |
| sequencing_platform | Sequencing platform | 10x_3prime_v2, 10x_3prime_v3, drop_seq, smart_seq2, other |
| species | Should always be human for included datasets | human |
| notes | Any additional relevant information | free text |

## Condition Harmonization

This is the most critical and ambiguous mapping. The following standardized categories should be used:

- **healthy**: Explicitly described as healthy, normal, no degeneration, no pain history. Must have supporting grading (Pfirrmann I, Thompson I-II) or explicit clinical confirmation.
- **degenerated_mild**: Mild degeneration. Pfirrmann II-III, Thompson II-III, or described as "mild" or "early" degeneration by the authors.
- **degenerated_severe**: Severe degeneration. Pfirrmann IV-V, Thompson III-IV or higher, or described as "severe" or "advanced" degeneration.
- **degenerated_ungraded**: Described as degenerated but without grading or clear severity classification.
- **herniated**: Specifically described as herniated, which may or may not overlap with degeneration categories. If both herniation and degeneration grade are reported, record both.
- **neonatal**: Neonatal tissue (not a disease state, but a distinct developmental stage).
- **aged_ungraded**: Aged donors without explicit degeneration grading. Important: "aged" is not synonymous with "degenerated" — some aged samples may be relatively healthy.

### Ambiguous cases to flag

- GSE251686: Samples are described as "mildly degenerative" (Pfirrmann II-III) AND have herniation. Should they be categorized as degenerated_mild, herniated, or both? Record both attributes.
- GSE189916 Sample_4-6: "Adult samples >65 years old, no pain history." These are aged but not explicitly degenerated. Categorize as aged_ungraded unless grading is available.
- GSE205535: "Normal" is from an 11-year-old with acute spinal cord injury. Is this truly "normal" IVD tissue? Flag for review.
- GSE233666: "IDD diagnosis" with no grading. Categorize as degenerated_ungraded.
- PMID35265617 Ctrl: Only 249 cells. Normal tissue (Pfirrmann I) but very low cell count. Flag potential QC concern.

## Age Harmonization

Where donor age in years is available, assign age_group:
- neonatal: < 1 year
- young_adult: 1-35 years
- middle_aged: 36-55 years
- aged: > 55 years
- unknown: age not reported

These cutoffs are approximate and may be revised at the human checkpoint.

## Automated Validation

- [ ] `metadata/sample_metadata.tsv` exists with one row per sample across all included studies
- [ ] All required fields are populated (no empty cells except where explicitly allowed as NA/unknown)
- [ ] Every sample in the registry marked as "included" appears in sample_metadata.tsv
- [ ] No duplicate sample_ids
- [ ] condition_harmonized values are all from the allowed set
- [ ] compartment values are all from the allowed set
- [ ] Donor counts and sample counts match what's reported in the original publications (cross-check against dataset_registry.tsv)
- [ ] `metadata/harmonization_notes.md` documents every non-trivial mapping decision

## Human Checkpoint

### Review materials
- `metadata/sample_metadata.tsv`
- `metadata/harmonization_notes.md`
- Summary statistics: number of samples per condition_harmonized, per compartment, per age_group

### Questions for the reviewer
1. Are the condition mappings accurate? Especially the ambiguous cases flagged above.
2. Is the condition hierarchy appropriate? Should "herniated" be a separate axis or folded into degeneration severity?
3. Is the age group binning appropriate for the scientific questions?
4. Are there any donor-level confounds (e.g., same donor contributing samples to multiple conditions) that need special handling?
5. Should GSE205535 "normal" (11 y/o with spinal cord injury) be reclassified or excluded?
6. Given the sample distribution across conditions and compartments, is the analysis plan still viable? Are any comparisons underpowered?

### Potential plan revisions
- If the healthy/normal sample count is very low, consider whether "mild degeneration" can serve as a comparison group for severe cases
- If compartment coverage is highly uneven, decide whether to focus the main analysis on NP (best covered) and treat AF/CEP as secondary
- If age and degeneration are strongly confounded in the available data, the plan for separating their effects may need adjustment
