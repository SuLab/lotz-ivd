# Metadata Harmonization Notes

Module 02 — IVD Single-Cell Atlas

Generated: 2026-02-26

Total samples: 78
Total studies: 12
Total unique donors: 57

## Summary Statistics

### Samples per condition_harmonized

| Condition | N samples |
|-----------|----------|
| aged_ungraded | 3 |
| degenerated_mild | 18 |
| degenerated_severe | 21 |
| degenerated_ungraded | 3 |
| healthy | 20 |
| herniated | 10 |
| neonatal | 3 |

### Samples per compartment

| Compartment | N samples |
|-------------|----------|
| AF | 17 |
| CEP | 6 |
| IVD_mixed | 6 |
| NP | 49 |

### Samples per age_group

| Age group | N samples |
|-----------|----------|
| aged | 26 |
| middle_aged | 11 |
| neonatal | 3 |
| unknown | 18 |
| young_adult | 20 |

### Donors per study

| Study | Donors | Samples |
|-------|--------|--------|
| CNP0002664 | 6 | 6 |
| GSE160756 | 2 | 7 |
| GSE165722 | 8 | 8 |
| GSE189916 | 4 | 6 |
| GSE199866 | 1 | 4 |
| GSE205535 | 2 | 2 |
| GSE230809 | 13 | 24 |
| GSE233666 | 4 | 4 |
| GSE242443 | 2 | 2 |
| GSE244889 | 7 | 7 |
| GSE251686 | 6 | 6 |
| GSE255768 | 2 | 2 |

---

## Condition Harmonization Decisions

### General Rules

- **Pfirrmann I** → `healthy` (normal disc)
- **Pfirrmann II-III** → `degenerated_mild`
- **Pfirrmann IV-V** → `degenerated_severe`
- **Thompson I-II** → `healthy`
- **Thompson II-III** → `degenerated_mild`
- **Thompson III alone** → `degenerated_mild` (boundary; conservative choice)
- **Thompson III-IV or higher** → `degenerated_severe`
- **Herniated discs with grading** → `herniated` (herniation takes priority as condition)
- **Neonatal tissue** → `neonatal`
- **Aged without grading** → `aged_ungraded`

### Ambiguous Cases (Flagged for Human Review)

#### 1. GSE251686 (Jia 2024): Herniation + degeneration grading

All 6 samples described as having both degeneration (mild: Pfirrmann II-III, severe: Pfirrmann IV) AND herniation (lumbar disc herniation). **Decision:** classified as `herniated` with degeneration_severity recorded separately (mild or severe). The spec explicitly flags this case: 'Should they be categorized as degenerated_mild, herniated, or both?' We chose `herniated` as the primary condition.

#### 2. GSE189916 (Jiang 2022): Adult samples — aged but ungraded

Adult samples are from donors >65 years old with no back pain history. No degeneration grading available. **Decision:** classified as `aged_ungraded` per spec guidance ('Aged donors without explicit degeneration grading'). These are NOT classified as 'healthy' because no grading confirmation exists, nor as 'degenerated' because no degeneration was diagnosed.

#### 3. GSE205535 (Li 2022): 'Normal' from spinal cord injury patient

The 'normal' (NNP) sample is from an 11-year-old with acute spinal cord injury. The disc itself was not degenerated, but the clinical context (spinal trauma) is unusual. **Decision:** classified as `healthy` since the disc was not degenerated, but **FLAGGED** for human review. The 11-year-old's disc biology may differ from adult healthy discs due to developmental stage, and the spinal cord injury may have affected the local environment.

#### 4. GSE233666 (Guo 2023): IDD with herniation, has Pfirrmann grading

All 4 samples are from IDD patients undergoing lumbar discectomy for herniation. Pfirrmann grades II-III are available. **Decision:** classified as `herniated` because the primary clinical presentation and reason for surgery was disc herniation. Degeneration grades recorded in degeneration_grade_original and degeneration_severity.

#### 5. GSE244889 (Chen 2024): Author vs. spec classification conflict

Authors classify grade 1-2 as 'MDD' (mild disc degeneration) and grade 3-4 as 'SDD' (severe disc degeneration). However, our spec classifies Pfirrmann I as 'healthy.' **Decision:** followed the spec (Pfirrmann I = healthy, II-III = mild, IV-V = severe) rather than the authors' grouping. Two samples (Pa-17F, Pb-55F) with Pfirrmann 1 are classified as `healthy` despite authors calling them MDD.

#### 6. CNP0002664 (Han 2022): Ctrl sample with very low cell count

The control sample (Pfirrmann I, normal) has only ~249 cells reported. This is extremely low and may not survive QC filtering. **Decision:** included in metadata but flagged for QC review in Module 03.

## Grading System Discrepancies

### GSE165722 GEO vs. Paper Grade Offset

GEO metadata lists Pfirrmann grades I-IV for the 8 samples, but the paper's Table 1 (PMC8787427) lists grades II-V. The paper explicitly states 'Pfirrmann grade I disc tissues were difficult to obtain.' **Decision:** used the paper's grading (II-V), which is authoritative. The GEO metadata has a systematic off-by-one error in Pfirrmann grades.

### Thompson III boundary classification

Thompson Grade III appears in both the 'mild' (II-III) and 'severe' (III-IV) ranges in the spec. For samples graded as Thompson III alone (not III-IV), we classified as `degenerated_mild` (conservative, lower boundary). For Thompson III-IV, we classified as `degenerated_severe`. This affects GSE230809 (Swahn 2024) samples SP21.011 and SP21.013 (Thompson III → mild) vs. SP20.002 etc. (Thompson III-IV → severe).

## Platform Notes

- **10x_3prime_v3**: 8 datasets (GSE160756, GSE189916, GSE233666, GSE244889, GSE255768, GSE230809, GSE242443; GSE199866 inferred as v2)
- **10x_3prime_v2**: GSE189916 (Chromium 3' v2), GSE199866 (inferred from HiSeq 4000 sequencer and CellRanger h5 format)
- **other (BD Rhapsody)**: GSE165722, GSE205535
- **other (Singleron Matrix)**: CNP0002664
- **other (Singleron GEXSCOPE)**: GSE251686 — NOTE: registry incorrectly lists as '10x Genomics'

## Missing Metadata Summary

| Dataset | Age Missing | Sex Missing | Notes |
|---------|-------------|-------------|-------|
| GSE160756 | No | Yes | Sex in Supplementary Table 1 (not accessible) |
| GSE165722 | No | No | Full demographics from Table 1 |
| GSE189916 | Partial | Partial | Neonatal age=0, adult >65 (exact unknown); neonatal sex=M, adult unknown |
| GSE199866 | Yes | Yes | Single donor; demographics in Supplementary Table S26 (not accessible) |
| GSE205535 | No | Yes | Ages from paper (11y, 81y) |
| CNP0002664 | Yes | Yes | Demographics in Supplementary Table S1 (not accessible) |
| GSE233666 | No | No | Full demographics from Table 1 |
| GSE244889 | No | No | Encoded in sample names |
| GSE251686 | Yes | Yes | Not in paper text or GEO; may be in supplementary |
| GSE255768 | No | No | From paper text |
| GSE230809 | No | No | Excellent GEO metadata |
| GSE242443 | Yes | Yes | Not available in paper or GEO |

## Donor-Level Confounds

### Multi-sample donors

- **GSE160756 (Gan 2021):** Donor Gan_D02 (31y) contributes 6 of 7 samples (NP_2, NP_3, CEP_1, CEP_2, AF_1, AF_2). Donor Gan_D01 (18y) contributes only NP_1.
- **GSE189916 (Jiang 2022):** Donor Jiang_Neo01 contributes all 3 neonatal samples (different spinal levels from same donor).
- **GSE199866 (Cherif 2022):** Single donor (Cherif_D01) contributes all 4 samples (paired healthy/degenerated × NP/AF).
- **GSE230809 (Swahn 2024):** 8 of 11 donors contribute paired NP+AF samples. 2 donors (SP20.002, SP20.006) contribute AF only.

### Age-disease confounding

- **GSE230809 (Swahn 2024):** Healthy donors are 21-27 years old; diseased donors are 37-73 years old. This is a strong age-disease confound that must be addressed during integration.
- **GSE205535 (Li 2022):** Normal donor is 11y, degenerative is 81y. Extreme age difference.

### Sex bias

- **GSE230809 (Swahn 2024):** ALL 11 donors are male. This is the largest dataset (24 samples) and will dominate sex-unaware analyses.
- Most other datasets have no sex information, limiting sex-stratified analyses.

## Special Considerations

### Culture-expanded cells (GSE242443)

Both samples from Kuchynsky 2024 are culture-expanded CEP cells (passage 1). Gene expression may be altered by in vitro conditions. Included per human decision at Module 01 checkpoint to improve sparse CEP coverage. tissue_or_cells field set to 'cells' to distinguish from fresh tissue.

### Normalized vs. raw counts (GSE165722)

GEO supplementary files for Tu 2022 may contain normalized counts rather than raw UMI counts. This must be verified during Module 03 preprocessing. If only normalized counts are available, raw data may need to be obtained from SRA.

### Published corrections (GSE205535)

Li Z 2022 has published corrections/corrigenda. These should be reviewed during preprocessing to determine if any data corrections are needed.

## Data Sources

### curated_metadata.xlsx (domain expert)

Per-sample cell counts for 8 of 12 datasets were obtained from
`data/curated_metadata.xlsx`, a spreadsheet prepared by a domain expert on the
project. This provided cell counts for datasets not covered by GEO metadata:
CNP0002664, GSE251686, GSE199866, GSE189916, GSE205535, GSE233666, GSE255768,
GSE230809. Cell counts for GSE160756, GSE165722, GSE244889, and GSE242443 were
not in the curated file (these were added in Module 01 and may not have been
reviewed by the expert yet).

Notable discrepancy: The curated file lists only 5 of 6 GSE251686 samples
(MI1/NP1, MI3/NP3, SE1/NP5, SE2/NP6, SE3/NP9). The third mild sample
(NP4/GSM7986003) is absent. Its cell count is marked NA.

### Low cell count samples flagged

Three samples have <500 cells and may not survive QC:
- **CNP0002664_Ctrl**: 249 cells (Pfirrmann I control)
- **GSE255768_S2**: 423 cells (degenerative endplate)
- **GSE230809_AF_SP20_002**: 467 cells (diseased AF, Thompson III-IV)

Additional samples with <1000 cells:
- **GSE230809_AF_SP21_016**: 587 cells
- **GSE230809_NP_SP21_016**: 1,205 cells

These should be monitored during QC in Module 03.

