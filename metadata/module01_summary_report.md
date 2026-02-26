# Module 01: Dataset Discovery & Acquisition — Summary Report

**Date:** 2026-02-26
**Status:** Complete (GEO datasets downloaded; non-GEO datasets pending)

## Search Summary

| Metric | Count |
|--------|-------|
| Databases searched | 7 (GEO, ArrayExpress, CellxGene, HCA, Single Cell Portal, PubMed, Google Scholar) |
| Total unique candidate datasets identified | 21 |
| Included | 13 |
| Borderline (flagged for review) | 2 |
| Excluded | 6 |

## Included Datasets (13 total)

| # | Accession | First Author | Year | Compartment | Conditions | Cells (reported) | Samples | Technology | Repository |
|---|-----------|-------------|------|-------------|------------|-------------------|---------|-----------|------------|
| 1 | GSE160756 | Gan Y | 2021 | NP, AF, CEP | Healthy young/adult | 108,108 | 7 | 10x Genomics | GEO |
| 2 | GSE165722 | Tu J | 2022 | NP | Pfirrmann II-V | 39,732 | 8 | BD Rhapsody | GEO |
| 3 | GSE189916 | Jiang W | 2022 | Whole IVD | Neonatal vs adult | ~18,000 | 6 | 10x Genomics | GEO |
| 4 | GSE199866 | Cherif H | 2022 | NP, inner AF | Paired degen/non-degen | 13,636 | 4 | 10x Genomics | GEO |
| 5 | GSE205535 | Li Z | 2022 | NP | Normal vs degenerative | 10,152 | 2 | BD Rhapsody | GEO |
| 6 | CNP0002664 | Han S | 2022 | NP | Normal/mild/severe | 30,300 | 6 | Singleron Matrix | CNGB |
| 7 | PRJCA007656 | Ling Z | 2022 | NP | Pfirrmann II-IV | 36,196 | 3 | 10x Genomics | NGDC |
| 8 | GSE233666 | Guo S | 2023 | NP | IDD (disc herniation) | N/A | 4 | 10x Genomics | GEO |
| 9 | GSE244889 | Chen F | 2024 | NP | Mild vs severe degen | 55,264 | 7 | 10x Genomics | GEO |
| 10 | PRJCA014236 | Wang D | 2023 | NP, AF | Pfirrmann I-V | 73,562 | 14 | 10x Genomics | GSA |
| 11 | GSE251686 | Jia S | 2024 | NP | Mild vs severe degen | 47,610 | 6 | 10x Genomics | GEO |
| 12 | GSE255768 | Shi C | 2024 | CEP/Endplate | Degenerative (Modic) | 8,534 | 2 | 10x Genomics | GEO |
| 13 | GSE230809 | Swahn H | 2024 | NP, AF | Healthy vs diseased | ~92,334 | 24 | 10x Genomics | GEO |

**Total estimated cells across all included datasets:** ~533,000+
**Total samples:** ~93

## Borderline Datasets (flagged for human review)

| Accession | First Author | Year | Issue |
|-----------|-------------|------|-------|
| GSE242443 | Kuchynsky K | 2024 | CEP cells were culture-expanded before scRNA-seq; not primary tissue profiling. Exclusion criterion #5 says "solely in vitro cultured cells with no primary tissue comparison" are excluded. These CEP cells were isolated from tissue but expanded in culture. Flagged for reviewer decision. |
| N/A | Zhou T | 2023 | Embryonic/fetal IVD (gestational weeks 7-11); 177,725 cells. Highly relevant for understanding NC-to-NP transitions but may not be appropriate for a degeneration-focused atlas. Accession not confirmed from abstracts. |

## Excluded Datasets (6)

| Candidate | Reason |
|-----------|--------|
| Fernandes LM (2020) | Data not publicly deposited (available from corresponding author on request only) |
| Zhang Y (2021) | Shares data with Han 2022 (both use CNP0002664) — excluded to avoid double-counting |
| Liu (2024) | Re-analysis of existing datasets (GSE233666, GSE205535, GSE189916, CNP0002664); no new primary data |
| Gao B (2022) GSE192789 | Primary scRNA-seq is mouse; human used only for validation |
| GSE211407 (2022) | Non-human species (rat) |
| GSE235198 (2024) | Non-human species (mouse) |

## Coverage Analysis

### By IVD Compartment

| Compartment | # Datasets | # Samples (approx) | Notes |
|-------------|-----------|---------------------|-------|
| Nucleus pulposus (NP) | 12 | ~75 | Well covered across conditions |
| Annulus fibrosus (AF) | 4 | ~20 | GSE160756, GSE199866, GSE230809, PRJCA014236 |
| Cartilage endplate (CEP) | 2-3 | ~4 | GSE160756 (2 samples), GSE255768 (2 samples), GSE242443 (borderline, 2 samples) |
| Whole IVD | 1 | 6 | GSE189916 (neonatal + adult) |

### By Condition

| Condition | # Datasets | Notes |
|-----------|-----------|-------|
| Healthy / non-degenerated | 4 | GSE160756, GSE189916 (adult), GSE199866 (non-degen), GSE230809 (healthy) |
| Mild degeneration (Pfirrmann I-II) | 5 | GSE165722, GSE244889, GSE251686, CNP0002664, PRJCA014236 |
| Severe degeneration (Pfirrmann III-V) | 7 | GSE165722, GSE244889, GSE251686, GSE233666, CNP0002664, PRJCA014236, PRJCA007656 |
| Neonatal / developmental | 1 | GSE189916 |
| Endplate-specific pathology | 1 | GSE255768 (Modic changes) |

### By Technology Platform

| Platform | # Datasets |
|----------|-----------|
| 10x Genomics | 10 |
| BD Rhapsody | 2 (GSE165722, GSE205535) |
| Singleron Matrix | 1 (CNP0002664) |

## Data Acquisition Status

### Successfully Downloaded (11 GEO datasets)

| Accession | File | Size | MD5 |
|-----------|------|------|-----|
| GSE160756 | GSE160756_RAW.tar | 433 MB | 54aa105250ba6770c147932ce44be3a8 |
| GSE165722 | GSE165722_RAW.tar | 67 MB | c5c6b7b5f00e47e3e1b8165fd129756b |
| GSE189916 | GSE189916_RAW.tar | 244 MB | 945a0c96de4f70650da276b8e991edb8 |
| GSE199866 | GSE199866_RAW.tar | 51 MB | 8b64d6152d7409867f30f223bed4cc38 |
| GSE205535 | GSE205535_RAW.tar | 61 MB | ad92b07aec19c0b454bcc523c3b3f62f |
| GSE233666 | GSE233666_RAW.tar | 230 MB | 3f598eb1f360e913959754c4a4dce9b8 |
| GSE244889 | GSE244889_RAW.tar | 521 MB | b8e532986390a8f046744bd5a83ff359 |
| GSE251686 | GSE251686_RAW.tar | 275 MB | 8f000dd3852fdd137b1c970540c5e0b2 |
| GSE255768 | GSE255768_RAW.tar | 46 MB | d0928d887a9b6c9e3a47dff49cc3e10d |
| GSE230809 | GSE230809_RAW.tar | 1,146 MB | 4e27be9a9e6e459e59dd3fa88e65e07c |
| GSE242443 | GSE242443_RAW.tar | 255 MB | 3e31f4340fac05a8e9323f34e6b4f62b |

**Total GEO download size:** ~3.3 GB

### Downloaded from CNGB (1 dataset)

| Accession | Files | Total Size | Notes |
|-----------|-------|-----------|-------|
| CNP0002664 (Han S 2022) | 6 count matrices (ctrl, NP2, NP4, NP8, NP9, NP10) | 55 MB | All MD5 checksums verified |

### Pending Download (2 Chinese repository datasets)

| Accession | Repository | Status |
|-----------|-----------|--------|
| PRJCA014236 (Wang D 2023) | Genome Sequence Archive (GSA-Human) | Requires NGDC account/access |
| PRJCA007656 (Ling Z 2022) | NGDC BioProject | Requires NGDC account/access |

## File Format Summary

| Format | Datasets |
|--------|----------|
| 10x MTX (barcodes.tsv.gz + features.tsv.gz + matrix.mtx.gz) | GSE189916, GSE205535, GSE230809, GSE233666, GSE242443, GSE244889, GSE255768 |
| 10x MTX in nested tar.gz | GSE251686 |
| Loom (.loom.gz) | GSE160756 |
| HDF5 (.h5) | GSE199866 |
| Count matrices (counts.tsv.gz + cellname.txt.gz) | GSE165722 |
| Count matrices (genes x cells TSV) | CNP0002664 |

## Automated Validation Checklist

- [x] All search queries executed and results logged in `metadata/search_log.md`
- [x] `metadata/dataset_registry.tsv` exists with all candidates, inclusion/exclusion status and reason
- [x] All included GEO datasets have raw data downloaded to `data/raw/`
- [ ] **PARTIAL:** 2 datasets in Chinese repositories not yet downloaded (PRJCA014236, PRJCA007656). CNP0002664 downloaded and verified.
- [x] File checksums recorded in `metadata/file_checksums.json`
- [x] No included dataset is a duplicate of another (Zhang Y 2021 excluded as duplicate of Han 2022)
- [x] Each known dataset from the original table is either included or has explicit exclusion reason
- [x] Summary report generated (this document)

## Key Observations

1. **6 new datasets discovered** beyond the 8 originally known: GSE160756 (Gan 2021), GSE165722 (Tu 2022), GSE244889 (Chen 2024), GSE242443 (Kuchynsky 2024), PRJCA014236 (Wang 2023), PRJCA007656 (Ling 2022)

2. **No IVD datasets on curated portals** — CellxGene, HCA, and Single Cell Portal have no IVD data. This atlas would be the first curated IVD resource on such a platform.

3. **NP is heavily overrepresented** — nearly all datasets profile NP. AF coverage is moderate (4 datasets), and CEP coverage is limited (2-3 datasets). This asymmetry should inform compartment-specific analyses.

4. **Platform heterogeneity** — Most use 10x Genomics, but 3 datasets use different platforms (BD Rhapsody, Singleron Matrix). Integration may require platform-aware batch correction.

5. **Chinese repository datasets** — 3 datasets are only available through Chinese national databases. These require separate access procedures and may need manual download.

6. **GSE205535 (Li Z 2022)** has published corrections/corrigenda — should be examined carefully during preprocessing.

7. **No snRNA-seq or spatial transcriptomics** datasets for human IVD were identified.
