# Dataset Discovery Search Log

## Search Date: 2026-02-26

## Databases Searched

1. **GEO (Gene Expression Omnibus)** — via NCBI web search and direct accession lookup
2. **ArrayExpress / BioStudies** — via EBI web search
3. **CellxGene** — via API query and web search
4. **Human Cell Atlas Data Portal** — via web search
5. **Single Cell Portal (Broad Institute)** — via web search
6. **PubMed** — systematic literature search
7. **Google Scholar** — supplementary literature search

## Search Queries and Results

### GEO / NCBI Searches

| # | Query | Results Found |
|---|-------|---------------|
| 1 | "intervertebral disc" single-cell RNA-seq | GSE251686, GSE160756 |
| 2 | "nucleus pulposus" single-cell RNA-seq | GSE165722, GSE160756 |
| 3 | "annulus fibrosus" single-cell RNA-seq | GSE160756 |
| 4 | "intervertebral disc" scRNA-seq | GSE160756, GSE244889 |
| 5 | "disc degeneration" single-cell | GSE165722, GSE244889 |
| 6 | "IVD" single-cell RNA | GSE199866, GSE205535 |
| 7 | "notochordal" single-cell | GSE189916 (human), GSE235198 (mouse) |
| 8 | "endplate" spine single-cell | GSE255768, GSE242443, GSE160756 |

Additional GEO datasets found through cross-referencing: GSE229711/GSE230808/GSE230809, GSE233666

### CellxGene

Searched via CellxGene API filtering by tissue labels including "intervertebral disc," "nucleus pulposus," "annulus fibrosus." **No intervertebral disc datasets currently hosted on CellxGene.** Closest tissue labels ("bone spine", "vertebral column") correspond to cancer/lemur datasets, not IVD studies.

### Human Cell Atlas Data Portal

**Zero results** for intervertebral disc. HCA does not currently include IVD projects.

### Single Cell Portal (Broad Institute)

**No direct hits** for intervertebral disc datasets.

### ArrayExpress / BioStudies

**No IVD-specific scRNA-seq datasets found** that were not already indexed in GEO.

### PubMed Searches

| # | Query | Unique Publications Found |
|---|-------|--------------------------|
| 1 | "intervertebral disc" AND ("single-cell" OR "scRNA-seq" OR "single-nucleus" OR "snRNA-seq") | ~25 results (includes reviews, non-human) |
| 2 | "nucleus pulposus" AND "single-cell" | ~15 results |
| 3 | "annulus fibrosus" AND "single-cell" | ~5 results |
| 4 | "disc degeneration" AND "scRNA-seq" | ~10 results |
| 5 | "IVD" AND "single-cell transcriptom*" | ~8 results |
| 6 | "notochordal cells" AND "single-cell" | ~5 results |

### Google Scholar

Supplementary search confirmed PubMed findings. No additional datasets beyond those found via PubMed/GEO.

## Cross-Referencing

- Checked citation lists of all included publications
- Checked review articles on IVD single-cell studies (including the fibrocyte enrichment meta-analysis, Bone Research 2024)
- Confirmed no curated IVD collections exist on CellxGene or HCA
- Identified Liu 2024 (Frontiers) as a re-analysis of existing datasets (GSE233666, GSE205535, GSE189916, CNP0002664) — no new primary data

## Key Findings

1. **No IVD datasets on major curated portals** (CellxGene, HCA, Single Cell Portal) — all data resides in primary repositories (GEO, CNGB, NGDC/GSA)
2. **No snRNA-seq studies** of human IVD identified — all use standard scRNA-seq
3. **No spatial transcriptomics** datasets for human IVD found (Zhou 2023 used mouse Visium only; human was scRNA-seq)
4. **Three datasets in Chinese repositories** (CNGB, NGDC/GSA) that are not mirrored in GEO
5. **Search identified 6 new datasets** beyond the 8 originally known: GSE160756 (Gan 2021), GSE165722 (Tu 2022), GSE244889 (Chen 2024), GSE242443 (Kuchynsky 2024), PRJCA014236 (Wang 2023), PRJCA007656 (Ling 2022)
