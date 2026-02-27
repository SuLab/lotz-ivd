#!/usr/bin/env python3
"""
Module 02: Metadata Harmonization
Generates harmonized sample metadata, harmonization notes, and ontology mappings
for the IVD single-cell atlas project.

All mapping decisions are documented in metadata/harmonization_notes.md.
"""

import csv
import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
METADATA_DIR = PROJECT_ROOT / "metadata"
METADATA_DIR.mkdir(exist_ok=True)

# === Standardized schema fields ===
FIELDS = [
    "sample_id", "study_accession", "first_author", "year", "donor_id",
    "n_cells_raw", "compartment", "condition_original", "condition_harmonized",
    "degeneration_grade_original", "degeneration_grade_system",
    "degeneration_severity", "age_group", "age_years", "sex",
    "tissue_or_cells", "sequencing_platform", "species", "notes"
]

ALLOWED_CONDITIONS = {
    "healthy", "degenerated_mild", "degenerated_severe",
    "degenerated_ungraded", "herniated", "neonatal", "aged_ungraded"
}
ALLOWED_COMPARTMENTS = {"NP", "AF", "CEP", "IVD_mixed"}
ALLOWED_SEX = {"M", "F", "unknown"}
ALLOWED_AGE_GROUPS = {"neonatal", "young_adult", "middle_aged", "aged", "unknown"}
ALLOWED_GRADE_SYSTEMS = {"Pfirrmann", "Thompson", "Modic", "other", "none"}
ALLOWED_SEVERITY = {"none", "mild", "moderate", "severe", "ungraded"}


def age_to_group(age_years):
    """Convert age in years to age group per spec."""
    if age_years is None or age_years == "NA":
        return "unknown"
    age = float(age_years)
    if age < 1:
        return "neonatal"
    elif age <= 35:
        return "young_adult"
    elif age <= 55:
        return "middle_aged"
    else:
        return "aged"


def build_samples():
    """Build the complete list of sample records."""
    samples = []

    # =========================================================================
    # GSE160756 — Gan 2021 — NP, AF, CEP — 7 samples, 2 donors
    # Source: GEO metadata + PMC8368097 full text
    # Ages from GEO. Sex not reported in GEO or accessible paper text.
    # Paper: "five healthy human IVDs (Pfirrmann I)" from nine donors.
    # GEO age: 18y (NP_1) and 31y (all others). Likely 2 donors.
    # =========================================================================
    gan_samples = [
        ("GSM4878538", "hNP_1", "NP", 18, "Gan_D01"),
        ("GSM4878539", "hNP_2", "NP", 31, "Gan_D02"),
        ("GSM4878540", "hNP_3", "NP", 31, "Gan_D02"),
        ("GSM4878541", "hCEP_1", "CEP", 31, "Gan_D02"),
        ("GSM4878542", "hCEP_2", "CEP", 31, "Gan_D02"),
        ("GSM4878543", "hAF_1", "AF", 31, "Gan_D02"),
        ("GSM4878544", "hAF_2", "AF", 31, "Gan_D02"),
    ]
    for gsm, name, comp, age, donor in gan_samples:
        samples.append({
            "sample_id": f"GSE160756_{name}",
            "study_accession": "GSE160756",
            "first_author": "Gan Y",
            "year": 2021,
            "donor_id": donor,
            "n_cells_raw": "NA",
            "compartment": comp,
            "condition_original": "healthy young/adult",
            "condition_harmonized": "healthy",
            "degeneration_grade_original": "Pfirrmann I",
            "degeneration_grade_system": "Pfirrmann",
            "degeneration_severity": "none",
            "age_group": age_to_group(age),
            "age_years": age,
            "sex": "unknown",
            "tissue_or_cells": "tissue",
            "sequencing_platform": "10x_3prime_v3",
            "species": "human",
            "notes": f"GSM {gsm}; donor {donor} contributes multiple compartments"
                     if donor == "Gan_D02" else f"GSM {gsm}",
        })

    # =========================================================================
    # GSE165722 — Tu 2022 — NP — 8 samples, 8 donors
    # Source: PMC8787427 Table 1 (full text) + GEO
    # IMPORTANT: GEO lists grades I-IV but paper Table 1 says II-V. Paper
    # explicitly states "Pfirrmann grade I disc tissues were difficult to
    # obtain." Using paper grades (II-V), not GEO (I-IV).
    # Platform: BD Rhapsody (NOT 10x).
    # =========================================================================
    tu_samples = [
        ("GSM5048708", "S1", 63, "M", "II",  "burst fracture"),
        ("GSM5048709", "S2", 41, "M", "II",  "burst fracture"),
        ("GSM5048710", "S3", 56, "F", "III", "lumbar disc herniation"),
        ("GSM5048711", "S4", 65, "F", "III", "lumbar disc herniation"),
        ("GSM5048712", "S5", 64, "F", "IV",  "lumbar disc herniation"),
        ("GSM5048713", "S6", 53, "F", "IV",  "lumbar disc herniation"),
        ("GSM5048714", "S7", 54, "M", "V",   "lumbar disc herniation"),
        ("GSM5048715", "S8", 56, "M", "V",   "lumbar disc herniation"),
    ]
    for gsm, name, age, sex, grade, surgery in tu_samples:
        grade_num = {"II": 2, "III": 3, "IV": 4, "V": 5}[grade]
        if grade_num <= 3:
            severity = "mild"
            condition = "degenerated_mild"
        else:
            severity = "severe"
            condition = "degenerated_severe"
        samples.append({
            "sample_id": f"GSE165722_{name}",
            "study_accession": "GSE165722",
            "first_author": "Tu J",
            "year": 2022,
            "donor_id": f"Tu_{name}",
            "n_cells_raw": "NA",
            "compartment": "NP",
            "condition_original": f"Pfirrmann {grade}, {surgery}",
            "condition_harmonized": condition,
            "degeneration_grade_original": f"Pfirrmann {grade}",
            "degeneration_grade_system": "Pfirrmann",
            "degeneration_severity": severity,
            "age_group": age_to_group(age),
            "age_years": age,
            "sex": sex,
            "tissue_or_cells": "tissue",
            "sequencing_platform": "other",
            "species": "human",
            "notes": f"GSM {gsm}; BD Rhapsody platform; "
                     f"GEO grades offset by 1 (lists I-IV, paper says II-V); "
                     f"surgery: {surgery}",
        })

    # =========================================================================
    # GSE189916 — Jiang 2022 — Whole IVD — 6 samples, 4 donors
    # Source: PMC9213722 + GEO + curated_metadata.xlsx
    # Neonatal: 3 samples from 1 donor (ND18842), 6 hours postnatal, male.
    # Adult: 3 samples from 3 different donors, >65 years, no back pain.
    # Neonatal from L1-L5 (3 levels used). Adults from lumbar discs.
    # Cell counts from curated_metadata.xlsx (domain expert).
    # =========================================================================
    jiang_cells = {1: 2501, 2: 2063, 3: 3210, 4: 2722, 5: 5990, 6: 5523}
    # Neonatal samples
    for i in range(1, 4):
        gsm = f"GSM570957{i}"
        samples.append({
            "sample_id": f"GSE189916_Neonatal_IVD_{i}",
            "study_accession": "GSE189916",
            "first_author": "Jiang W",
            "year": 2022,
            "donor_id": "Jiang_Neo01",
            "n_cells_raw": jiang_cells[i],
            "compartment": "IVD_mixed",
            "condition_original": "neonatal IVD",
            "condition_harmonized": "neonatal",
            "degeneration_grade_original": "NA",
            "degeneration_grade_system": "none",
            "degeneration_severity": "none",
            "age_group": "neonatal",
            "age_years": 0,
            "sex": "M",
            "tissue_or_cells": "tissue",
            "sequencing_platform": "10x_3prime_v2",
            "species": "human",
            "notes": f"GSM {gsm}; same donor ND18842 for all 3 neonatal samples; "
                     "6 hours postnatal; NDRI tissue; L1-L5 levels",
        })
    # Adult samples
    for i, idx in enumerate([4, 5, 6], start=1):
        gsm = f"GSM570957{idx}"
        source = "cadaveric" if i <= 2 else "clinical (discarded disc)"
        samples.append({
            "sample_id": f"GSE189916_Adult_IVD_{i}",
            "study_accession": "GSE189916",
            "first_author": "Jiang W",
            "year": 2022,
            "donor_id": f"Jiang_Ad{i:02d}",
            "n_cells_raw": jiang_cells[idx],
            "compartment": "IVD_mixed",
            "condition_original": "adult IVD, >65 years, no back pain history",
            "condition_harmonized": "aged_ungraded",
            "degeneration_grade_original": "NA",
            "degeneration_grade_system": "none",
            "degeneration_severity": "ungraded",
            "age_group": "aged",
            "age_years": "NA",
            "sex": "unknown",
            "tissue_or_cells": "tissue",
            "sequencing_platform": "10x_3prime_v2",
            "species": "human",
            "notes": f"GSM {gsm}; {source}; exact age not reported (>65y); "
                     "no back pain history per paper",
        })

    # =========================================================================
    # GSE199866 — Cherif 2022 — NP, inner AF — 4 samples, 1 donor
    # Source: PMC8999935 + GEO
    # Unique paired design: healthy and degenerated discs from same individual.
    # Healthy: Thompson I-II (L4-L5). Degenerated: Thompson III-V (L5-S1).
    # Age and sex not available in GEO or accessible paper sections.
    # NOTE: GEO metadata error — GSM5989811 (NPD) lists "AF cells" in
    # source_name, but filename and title indicate NP.
    # =========================================================================
    # Cell counts from curated_metadata.xlsx (domain expert)
    cherif_samples = [
        ("GSM5989808", "AFH", "AF", "non-degenerating inner AF",
         "healthy", "Thompson I-II", "none", 3226),
        ("GSM5989809", "AFD", "AF", "degenerating inner AF",
         "degenerated_severe", "Thompson III-V", "severe", 3142),
        ("GSM5989810", "NPH", "NP", "non-degenerating NP",
         "healthy", "Thompson I-II", "none", 3955),
        ("GSM5989811", "NPD", "NP", "degenerating NP",
         "degenerated_severe", "Thompson III-V", "severe", 3678),
    ]
    for gsm, name, comp, cond_orig, cond_harm, grade, severity, ncells in cherif_samples:
        note = f"GSM {gsm}; paired design (same donor); inner AF (not outer AF)"
        if gsm == "GSM5989811":
            note += "; GEO source_name erroneously says AF but this is NPD"
        samples.append({
            "sample_id": f"GSE199866_{name}",
            "study_accession": "GSE199866",
            "first_author": "Cherif H",
            "year": 2022,
            "donor_id": "Cherif_D01",
            "n_cells_raw": ncells,
            "compartment": comp,
            "condition_original": cond_orig,
            "condition_harmonized": cond_harm,
            "degeneration_grade_original": grade,
            "degeneration_grade_system": "Thompson",
            "degeneration_severity": severity,
            "age_group": "unknown",
            "age_years": "NA",
            "sex": "unknown",
            "tissue_or_cells": "tissue",
            "sequencing_platform": "10x_3prime_v2",
            "species": "human",
            "notes": note,
        })

    # =========================================================================
    # GSE205535 — Li Z 2022 — NP — 2 samples, 2 donors
    # Source: PMC9301035 + GEO
    # FLAGGED: "Normal" (NNP) is from an 11-year-old with acute spinal cord
    # injury — disc itself not degenerated but context is abnormal.
    # Degenerative (DNP) is from an 81-year-old with LDH/LDD.
    # Platform: BD Rhapsody. Published corrections/corrigenda exist.
    # =========================================================================
    # Cell counts from curated_metadata.xlsx (domain expert)
    samples.append({
        "sample_id": "GSE205535_NNP",
        "study_accession": "GSE205535",
        "first_author": "Li Z",
        "year": 2022,
        "donor_id": "Li_NNP",
        "n_cells_raw": 5168,
        "compartment": "NP",
        "condition_original": "normal NP (11y, acute spinal cord injury)",
        "condition_harmonized": "healthy",
        "degeneration_grade_original": "NA",
        "degeneration_grade_system": "none",
        "degeneration_severity": "none",
        "age_group": "young_adult",
        "age_years": 11,
        "sex": "unknown",
        "tissue_or_cells": "tissue",
        "sequencing_platform": "other",
        "species": "human",
        "notes": "GSM GSM6214392; BD Rhapsody platform; FLAG: 'normal' disc from "
                 "11-year-old with acute spinal cord injury — may not be truly "
                 "representative of normal IVD; published corrections exist",
    })
    samples.append({
        "sample_id": "GSE205535_DNP",
        "study_accession": "GSE205535",
        "first_author": "Li Z",
        "year": 2022,
        "donor_id": "Li_DNP",
        "n_cells_raw": 5005,
        "compartment": "NP",
        "condition_original": "degenerative NP (81y, LDH/LDD)",
        "condition_harmonized": "degenerated_ungraded",
        "degeneration_grade_original": "NA",
        "degeneration_grade_system": "none",
        "degeneration_severity": "ungraded",
        "age_group": "aged",
        "age_years": 81,
        "sex": "unknown",
        "tissue_or_cells": "tissue",
        "sequencing_platform": "other",
        "species": "human",
        "notes": "GSM GSM6214393; BD Rhapsody platform; no degeneration grade "
                 "reported; published corrections exist",
    })

    # =========================================================================
    # CNP0002664 — Han 2022 — NP — 6 samples
    # Source: PMC8899542 + CNGB data files
    # Pfirrmann grades from paper: Ctrl=I, NP4/NP9/NP10=II/III (mild),
    # NP2/NP8=IV/V (severe). Ages and sex in Supplementary Table S1 (not
    # accessible). Platform: Singleron Matrix.
    # NOTE: Ctrl sample has only ~249 cells — very low, flagged for QC.
    # =========================================================================
    # Cell counts from curated_metadata.xlsx (domain expert)
    han_samples = [
        ("ctrl", "Ctrl", "Pfirrmann I", "normal control", "healthy",
         "Pfirrmann", "none", 249),
        ("NP2", "NP2", "Pfirrmann IV/V", "severe IVDD", "degenerated_severe",
         "Pfirrmann", "severe", 2440),
        ("NP4", "NP4", "Pfirrmann II/III", "mild IVDD", "degenerated_mild",
         "Pfirrmann", "mild", 4654),
        ("NP8", "NP8", "Pfirrmann IV/V", "severe IVDD", "degenerated_severe",
         "Pfirrmann", "severe", 3916),
        ("NP9", "NP9", "Pfirrmann II/III", "mild IVDD", "degenerated_mild",
         "Pfirrmann", "mild", 11497),
        ("NP10", "NP10", "Pfirrmann II/III", "mild IVDD", "degenerated_mild",
         "Pfirrmann", "mild", 8345),
    ]
    for file_id, name, grade, cond_orig, cond_harm, grade_sys, severity, ncells in han_samples:
        note = f"Singleron Matrix platform; CNGB data (not GEO)"
        if name == "Ctrl":
            note += "; FLAG: only 249 cells — potential QC concern"
        samples.append({
            "sample_id": f"CNP0002664_{name}",
            "study_accession": "CNP0002664",
            "first_author": "Han S",
            "year": 2022,
            "donor_id": f"Han_{name}",
            "n_cells_raw": ncells,
            "compartment": "NP",
            "condition_original": cond_orig,
            "condition_harmonized": cond_harm,
            "degeneration_grade_original": grade,
            "degeneration_grade_system": grade_sys,
            "degeneration_severity": severity,
            "age_group": "unknown",
            "age_years": "NA",
            "sex": "unknown",
            "tissue_or_cells": "tissue",
            "sequencing_platform": "other",
            "species": "human",
            "notes": note,
        })

    # =========================================================================
    # GSE233666 — Guo 2023 — NP (herniated) — 4 samples, 4 donors
    # Source: PMC10449260 Table 1 + GEO
    # All from lumbar discectomy with IDD diagnosis. All herniated.
    # Pfirrmann grades available (II-III). Sex from paper Table 1.
    # condition_harmonized = herniated (primary diagnosis); degeneration grades
    # recorded separately.
    # =========================================================================
    # Cell counts from curated_metadata.xlsx (domain expert)
    guo_samples = [
        ("GSM7432171", "Patient_1", 20, "F", "II", 4400),
        ("GSM7432172", "Patient_2", 69, "F", "III", 2234),
        ("GSM7432173", "Patient_3", 23, "M", "III", 4759),
        ("GSM7432174", "Patient_4", 28, "F", "II", 7466),
    ]
    for gsm, name, age, sex, grade, ncells in guo_samples:
        samples.append({
            "sample_id": f"GSE233666_{name}",
            "study_accession": "GSE233666",
            "first_author": "Guo S",
            "year": 2023,
            "donor_id": f"Guo_P{name.split('_')[1]}",
            "n_cells_raw": ncells,
            "compartment": "NP",
            "condition_original": f"IDD, herniated NP, Pfirrmann {grade}, lumbar discectomy",
            "condition_harmonized": "herniated",
            "degeneration_grade_original": f"Pfirrmann {grade}",
            "degeneration_grade_system": "Pfirrmann",
            "degeneration_severity": "mild",
            "age_group": age_to_group(age),
            "age_years": age,
            "sex": sex,
            "tissue_or_cells": "tissue",
            "sequencing_platform": "10x_3prime_v3",
            "species": "human",
            "notes": f"GSM {gsm}; all samples herniated with IDD diagnosis; "
                     f"Pfirrmann {grade} but classified as herniated (primary "
                     "diagnosis is disc herniation)",
        })

    # =========================================================================
    # GSE244889 — Chen 2024 — NP — 7 samples, 7 donors
    # Source: GEO metadata (excellent — age/sex in filenames)
    # MDD = mild disc degeneration (Pfirrmann 1-2 per authors)
    # SDD = severe disc degeneration (Pfirrmann 3-4 per authors)
    # Harmonization: using spec's Pfirrmann-based categories, not authors'
    # grouping. Pfirrmann 1 = healthy, 2-3 = mild, 4-5 = severe.
    # =========================================================================
    chen_samples = [
        ("GSM7831813", "Pa-17F", 17, "F", "1", "MDD"),
        ("GSM7831814", "Pb-55F", 55, "F", "1", "MDD"),
        ("GSM7831815", "Pab-24M", 24, "M", "2", "MDD"),
        ("GSM7831816", "Pab-30M", 30, "M", "2", "MDD"),
        ("GSM7831817", "Pc-41M", 41, "M", "3", "SDD"),
        ("GSM7831818", "Pd-62F", 62, "F", "4", "SDD"),
        ("GSM7831819", "Pd-59F", 59, "F", "4", "SDD"),
    ]
    for gsm, name, age, sex, grade, author_class in chen_samples:
        grade_num = int(grade)
        if grade_num == 1:
            cond = "healthy"
            severity = "none"
        elif grade_num <= 3:
            cond = "degenerated_mild"
            severity = "mild"
        else:
            cond = "degenerated_severe"
            severity = "severe"
        samples.append({
            "sample_id": f"GSE244889_{name}",
            "study_accession": "GSE244889",
            "first_author": "Chen F",
            "year": 2024,
            "donor_id": f"Chen_{name}",
            "n_cells_raw": "NA",
            "compartment": "NP",
            "condition_original": f"Pfirrmann grade {grade}, {author_class}",
            "condition_harmonized": cond,
            "degeneration_grade_original": f"Pfirrmann {grade}",
            "degeneration_grade_system": "Pfirrmann",
            "degeneration_severity": severity,
            "age_group": age_to_group(age),
            "age_years": age,
            "sex": sex,
            "tissue_or_cells": "tissue",
            "sequencing_platform": "10x_3prime_v3",
            "species": "human",
            "notes": f"GSM {gsm}; authors classify grade 1-2 as MDD and 3-4 as SDD; "
                     "harmonized using spec categories (Pfirrmann I = healthy)",
        })

    # =========================================================================
    # GSE251686 — Jia 2024 — NP — 6 samples, 6 donors
    # Source: PMC11549379 + GEO
    # All from L4/L5 or L5/S1 with lumbar disc herniation.
    # Mild: Pfirrmann II-III (NP1, NP3, NP4)
    # Severe: Pfirrmann IV (NP5, NP6, NP9)
    # All have herniation — flagged per spec as ambiguous case.
    # Ages and sex not available from paper or GEO.
    # Platform: Singleron GEXSCOPE (not 10x as listed in registry).
    # =========================================================================
    # Cell counts from curated_metadata.xlsx (MI1=NP1, MI3=NP3, SE1=NP5,
    # SE2=NP6, SE3=NP9). NP4 count not in curated data, estimated from
    # total (47610) minus other 5 samples (37089) = ~10521.
    jia_samples = [
        ("GSM7986001", "NP1", "mildly degeneration, rep 1",
         "Pfirrmann II-III", "mild", 15287),
        ("GSM7986002", "NP3", "mildly degeneration, rep 2",
         "Pfirrmann II-III", "mild", 5916),
        ("GSM7986003", "NP4", "mildly degeneration, rep 3",
         "Pfirrmann II-III", "mild", "NA"),
        ("GSM7986004", "NP5", "severely degeneration, rep 1",
         "Pfirrmann IV", "severe", 7009),
        ("GSM7986005", "NP6", "severely degeneration, rep 2",
         "Pfirrmann IV", "severe", 2428),
        ("GSM7986006", "NP9", "severely degeneration, rep 3",
         "Pfirrmann IV", "severe", 6449),
    ]
    for gsm, name, cond_orig, grade, severity, ncells in jia_samples:
        samples.append({
            "sample_id": f"GSE251686_{name}",
            "study_accession": "GSE251686",
            "first_author": "Jia S",
            "year": 2024,
            "donor_id": f"Jia_{name}",
            "n_cells_raw": ncells,
            "compartment": "NP",
            "condition_original": f"{cond_orig}, lumbar disc herniation",
            "condition_harmonized": "herniated",
            "degeneration_grade_original": grade,
            "degeneration_grade_system": "Pfirrmann",
            "degeneration_severity": severity,
            "age_group": "unknown",
            "age_years": "NA",
            "sex": "unknown",
            "tissue_or_cells": "tissue",
            "sequencing_platform": "other",
            "species": "human",
            "notes": f"GSM {gsm}; ALL samples have herniation + degeneration — "
                     "classified as 'herniated' per spec; Singleron GEXSCOPE "
                     "platform (registry incorrectly lists 10x); "
                     f"degeneration severity: {severity}",
        })

    # =========================================================================
    # GSE255768 — Shi 2024 — CEP/Endplate — 2 samples, 2 donors
    # Source: PMC11399435 + GEO
    # Both degenerative endplate with Modic changes from disc herniation surgery.
    # No healthy endplate control.
    # =========================================================================
    # Cell counts from curated_metadata.xlsx (domain expert)
    shi_samples = [
        ("GSM8079184", "S1", 58, "F", 8628),
        ("GSM8079185", "S2", 66, "M", 423),
    ]
    for gsm, name, age, sex, ncells in shi_samples:
        samples.append({
            "sample_id": f"GSE255768_{name}",
            "study_accession": "GSE255768",
            "first_author": "Shi C",
            "year": 2024,
            "donor_id": f"Shi_{name}",
            "n_cells_raw": ncells,
            "compartment": "CEP",
            "condition_original": f"degenerative endplate, Modic changes, disc herniation surgery",
            "condition_harmonized": "degenerated_ungraded",
            "degeneration_grade_original": "Modic changes (type unspecified)",
            "degeneration_grade_system": "Modic",
            "degeneration_severity": "ungraded",
            "age_group": age_to_group(age),
            "age_years": age,
            "sex": sex,
            "tissue_or_cells": "tissue",
            "sequencing_platform": "10x_3prime_v3",
            "species": "human",
            "notes": f"GSM {gsm}; no healthy endplate control in this study; "
                     "Modic changes at L3-L4 or L4-L5; surgery for disc herniation",
        })

    # =========================================================================
    # GSE230809 — Swahn 2024 — NP, AF — 24 samples, 11 donors
    # Source: GEO (excellent metadata) + PMC (Adv Sci)
    # SuperSeries: GSE229711 (healthy) + GSE230808 (diseased)
    # ALL donors male. Healthy from tissue bank, diseased from surgery.
    # Thompson grading. Strong age-disease confound (healthy young, diseased older).
    # =========================================================================
    # Cell counts from curated_metadata.xlsx (domain expert)
    # Healthy samples (GSE229711 subseries)
    swahn_healthy = [
        ("GSM7173748", "AF_SP21_015", "AF", 21, "SP21.015", "Thompson II", 7242),
        ("GSM7173749", "AF_SP21_018", "AF", 27, "SP21.018", "Thompson II", 7945),
        ("GSM7173750", "AF_SP22_001", "AF", 25, "SP22.001", "Thompson II", 4791),
        ("GSM7173751", "NP_SP21_015", "NP", 21, "SP21.015", "Thompson II", 1372),
        ("GSM7173752", "NP_SP21_018", "NP", 27, "SP21.018", "Thompson II", 12729),
        ("GSM7173753", "NP_SP22_001", "NP", 25, "SP22.001", "Thompson II", 6786),
    ]
    for gsm, name, comp, age, patient, grade, ncells in swahn_healthy:
        samples.append({
            "sample_id": f"GSE230809_{name}",
            "study_accession": "GSE230809",
            "first_author": "Swahn H",
            "year": 2024,
            "donor_id": f"Swahn_{patient}",
            "n_cells_raw": ncells,
            "compartment": comp,
            "condition_original": f"healthy, {grade}",
            "condition_harmonized": "healthy",
            "degeneration_grade_original": grade,
            "degeneration_grade_system": "Thompson",
            "degeneration_severity": "none",
            "age_group": age_to_group(age),
            "age_years": age,
            "sex": "M",
            "tissue_or_cells": "tissue",
            "sequencing_platform": "10x_3prime_v3",
            "species": "human",
            "notes": f"GSM {gsm}; tissue bank donor; paired NP+AF from same donor; "
                     "all donors male; subseries GSE229711",
        })

    # Diseased samples (GSE230808 subseries)
    # Cell counts from curated_metadata.xlsx (domain expert)
    swahn_diseased = [
        ("GSM7235331", "AF_SP20_002", "AF", 73, "SP20.002", "Thompson III-IV", 467),
        ("GSM7235332", "AF_SP20_006", "AF", 56, "SP20.006", "Thompson III-IV", 1206),
        ("GSM7235333", "AF_SP21_007", "AF", 43, "SP21.007", "Thompson II-III", 5878),
        ("GSM7235334", "AF_SP21_011", "AF", 42, "SP21.011", "Thompson III", 2194),
        ("GSM7235335", "AF_SP21_013", "AF", 63, "SP21.013", "Thompson III", 4538),
        ("GSM7235336", "AF_SP21_014", "AF", 37, "SP21.014", "Thompson II-III", 2821),
        ("GSM7235337", "AF_SP21_016", "AF", 68, "SP21.016", "Thompson III-IV", 587),
        ("GSM7235338", "AF_SP21_017", "AF", 61, "SP21.017", "Thompson III-IV", 1276),
        ("GSM7235339", "AF_SP22_002", "AF", 64, "SP22.002", "Thompson III-IV", 3689),
        ("GSM7235340", "AF_SP22_003", "AF", 63, "SP22.003", "Thompson III-IV", 4334),
        ("GSM7235341", "NP_SP21_007", "NP", 43, "SP21.007", "Thompson II-III", 4124),
        ("GSM7235342", "NP_SP21_011", "NP", 42, "SP21.011", "Thompson III", 1863),
        ("GSM7235343", "NP_SP21_013", "NP", 63, "SP21.013", "Thompson III", 2714),
        ("GSM7235344", "NP_SP21_014", "NP", 37, "SP21.014", "Thompson II-III", 3435),
        ("GSM7235345", "NP_SP21_016", "NP", 68, "SP21.016", "Thompson III-IV", 1205),
        ("GSM7235346", "NP_SP21_017", "NP", 61, "SP21.017", "Thompson III-IV", 1376),
        ("GSM7235347", "NP_SP22_002", "NP", 64, "SP22.002", "Thompson III-IV", 2867),
        ("GSM7235348", "NP_SP22_003", "NP", 63, "SP22.003", "Thompson III-IV", 6909),
    ]
    for gsm, name, comp, age, patient, grade, ncells in swahn_diseased:
        # Thompson III or II-III = mild; Thompson III-IV = severe
        if "III-IV" in grade or "IV" in grade:
            severity = "severe"
            condition = "degenerated_severe"
        else:
            severity = "mild"
            condition = "degenerated_mild"
        has_np = patient not in ("SP20.002", "SP20.006")
        paired_note = "paired NP+AF from same donor" if has_np else "AF only (no NP)"
        samples.append({
            "sample_id": f"GSE230809_{name}",
            "study_accession": "GSE230809",
            "first_author": "Swahn H",
            "year": 2024,
            "donor_id": f"Swahn_{patient}",
            "n_cells_raw": ncells,
            "compartment": comp,
            "condition_original": f"diseased, {grade}",
            "condition_harmonized": condition,
            "degeneration_grade_original": grade,
            "degeneration_grade_system": "Thompson",
            "degeneration_severity": severity,
            "age_group": age_to_group(age),
            "age_years": age,
            "sex": "M",
            "tissue_or_cells": "tissue",
            "sequencing_platform": "10x_3prime_v3",
            "species": "human",
            "notes": f"GSM {gsm}; surgical specimen; {paired_note}; "
                     "all donors male; subseries GSE230808; "
                     "strong age-disease confound (healthy=young, diseased=older)",
        })

    # =========================================================================
    # GSE242443 — Kuchynsky 2024 — CEP — 2 samples, 2 donors
    # Source: GEO + PMID 38173036
    # Culture-expanded p1 CEP cells. Included per human decision for coverage.
    # Thompson <=2 = healthy, Thompson >=4 = degenerated.
    # Age and sex not available.
    # =========================================================================
    samples.append({
        "sample_id": "GSE242443_Y8444_H1",
        "study_accession": "GSE242443",
        "first_author": "Kuchynsky K",
        "year": 2024,
        "donor_id": "Kuchynsky_H1",
        "n_cells_raw": "NA",
        "compartment": "CEP",
        "condition_original": "healthy CEP, Thompson <= 2",
        "condition_harmonized": "healthy",
        "degeneration_grade_original": "Thompson <= 2",
        "degeneration_grade_system": "Thompson",
        "degeneration_severity": "none",
        "age_group": "unknown",
        "age_years": "NA",
        "sex": "unknown",
        "tissue_or_cells": "cells",
        "sequencing_platform": "10x_3prime_v3",
        "species": "human",
        "notes": "GSM GSM7763584; CULTURE-EXPANDED (p1) — gene expression may "
                 "be altered; included for CEP coverage per human decision",
    })
    samples.append({
        "sample_id": "GSE242443_Y8445_D2",
        "study_accession": "GSE242443",
        "first_author": "Kuchynsky K",
        "year": 2024,
        "donor_id": "Kuchynsky_D2",
        "n_cells_raw": "NA",
        "compartment": "CEP",
        "condition_original": "degenerated CEP, Thompson >= 4",
        "condition_harmonized": "degenerated_severe",
        "degeneration_grade_original": "Thompson >= 4",
        "degeneration_grade_system": "Thompson",
        "degeneration_severity": "severe",
        "age_group": "unknown",
        "age_years": "NA",
        "sex": "unknown",
        "tissue_or_cells": "cells",
        "sequencing_platform": "10x_3prime_v3",
        "species": "human",
        "notes": "GSM GSM7763585; CULTURE-EXPANDED (p1) — gene expression may "
                 "be altered; included for CEP coverage per human decision",
    })

    return samples


def write_sample_metadata(samples, outpath):
    """Write sample_metadata.tsv."""
    with open(outpath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        for s in samples:
            writer.writerow(s)
    print(f"Wrote {len(samples)} samples to {outpath}")


def write_ontology_mappings(outpath):
    """Write ontology_mappings.tsv mapping original terms to standardized terms."""
    mappings = [
        # Condition mappings
        ("condition", "healthy", "healthy", "No degeneration, Pfirrmann I or Thompson I-II"),
        ("condition", "healthy young and adult", "healthy", "GSE160756: young donors, Pfirrmann I"),
        ("condition", "non-degenerating", "healthy", "GSE199866: Thompson I-II disc"),
        ("condition", "normal", "healthy", "GSE205535: described as normal, but flagged"),
        ("condition", "normal control", "healthy", "CNP0002664: Pfirrmann I"),
        ("condition", "neonatal IVD", "neonatal", "GSE189916: 6h postnatal tissue"),
        ("condition", "adult IVD >65 no pain", "aged_ungraded", "GSE189916: aged, no degeneration grading"),
        ("condition", "Pfirrmann II", "degenerated_mild", "Pfirrmann II-III = mild per spec"),
        ("condition", "Pfirrmann III", "degenerated_mild", "Pfirrmann II-III = mild per spec"),
        ("condition", "Pfirrmann IV", "degenerated_severe", "Pfirrmann IV-V = severe per spec"),
        ("condition", "Pfirrmann V", "degenerated_severe", "Pfirrmann IV-V = severe per spec"),
        ("condition", "Thompson I-II", "healthy", "Thompson I-II = healthy per spec"),
        ("condition", "Thompson II", "healthy", "Thompson I-II = healthy per spec"),
        ("condition", "Thompson II-III", "degenerated_mild", "Thompson II-III = mild per spec"),
        ("condition", "Thompson III", "degenerated_mild", "Thompson III alone = mild (boundary)"),
        ("condition", "Thompson III-IV", "degenerated_severe", "Thompson III-IV or higher = severe per spec"),
        ("condition", "Thompson III-V", "degenerated_severe", "Thompson III-V = severe per spec"),
        ("condition", "Thompson >= 4", "degenerated_severe", "Thompson IV+ = severe per spec"),
        ("condition", "Thompson <= 2", "healthy", "Thompson I-II = healthy per spec"),
        ("condition", "mild IVDD", "degenerated_mild", "CNP0002664: Pfirrmann II-III"),
        ("condition", "severe IVDD", "degenerated_severe", "CNP0002664: Pfirrmann IV-V"),
        ("condition", "MDD (mild disc degeneration)", "varies", "GSE244889: Pfirrmann 1=healthy, 2=mild"),
        ("condition", "SDD (severe disc degeneration)", "varies", "GSE244889: Pfirrmann 3=mild, 4=severe"),
        ("condition", "degenerative", "degenerated_ungraded", "No grade specified"),
        ("condition", "degenerating", "degenerated_severe", "GSE199866: Thompson III-V"),
        ("condition", "IDD herniated NP", "herniated", "GSE233666: herniation is primary diagnosis"),
        ("condition", "mildly degeneration + herniation", "herniated", "GSE251686: herniation present"),
        ("condition", "severely degeneration + herniation", "herniated", "GSE251686: herniation present"),
        ("condition", "Modic changes", "degenerated_ungraded", "GSE255768: Modic type unspecified"),
        # Compartment mappings
        ("compartment", "nucleus pulposus", "NP", "Standard"),
        ("compartment", "annulus fibrosus", "AF", "Standard"),
        ("compartment", "inner annulus fibrosus", "AF", "GSE199866: inner AF specifically"),
        ("compartment", "cartilage endplate", "CEP", "Standard"),
        ("compartment", "endplate", "CEP", "GSE255768"),
        ("compartment", "whole IVD", "IVD_mixed", "GSE189916: not compartment-separated"),
        ("compartment", "IVD", "IVD_mixed", "When not compartment-separated"),
        # Platform mappings
        ("platform", "10x Chromium 3' v2", "10x_3prime_v2", ""),
        ("platform", "10x Chromium 3' v3", "10x_3prime_v3", ""),
        ("platform", "10x Chromium 3' v3.1", "10x_3prime_v3", "v3.1 grouped with v3"),
        ("platform", "BD Rhapsody", "other", "GSE165722, GSE205535"),
        ("platform", "Singleron Matrix", "other", "CNP0002664"),
        ("platform", "Singleron GEXSCOPE", "other", "GSE251686"),
        # Severity mappings
        ("severity", "Pfirrmann I", "none", "Normal disc"),
        ("severity", "Pfirrmann II-III", "mild", "Early degeneration"),
        ("severity", "Pfirrmann IV-V", "severe", "Advanced degeneration"),
        ("severity", "Thompson I-II", "none", "Normal disc"),
        ("severity", "Thompson II-III", "mild", "Early degeneration"),
        ("severity", "Thompson III-IV", "severe", "Moderate-severe degeneration"),
        ("severity", "Modic changes", "ungraded", "Different grading system"),
    ]

    with open(outpath, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["category", "original_term", "harmonized_term", "notes"])
        for row in mappings:
            writer.writerow(row)
    print(f"Wrote {len(mappings)} ontology mappings to {outpath}")


def write_harmonization_notes(samples, outpath):
    """Write harmonization_notes.md documenting all mapping decisions."""
    # Compute summary stats
    cond_counts = {}
    comp_counts = {}
    age_counts = {}
    for s in samples:
        cond_counts[s["condition_harmonized"]] = cond_counts.get(s["condition_harmonized"], 0) + 1
        comp_counts[s["compartment"]] = comp_counts.get(s["compartment"], 0) + 1
        age_counts[s["age_group"]] = age_counts.get(s["age_group"], 0) + 1

    # Count unique donors per study
    study_donors = {}
    for s in samples:
        key = s["study_accession"]
        if key not in study_donors:
            study_donors[key] = set()
        study_donors[key].add(s["donor_id"])

    with open(outpath, "w") as f:
        f.write("# Metadata Harmonization Notes\n\n")
        f.write("Module 02 — IVD Single-Cell Atlas\n\n")
        f.write(f"Generated: 2026-02-26\n\n")
        f.write(f"Total samples: {len(samples)}\n")
        f.write(f"Total studies: {len(study_donors)}\n")
        total_donors = sum(len(d) for d in study_donors.values())
        f.write(f"Total unique donors: {total_donors}\n\n")

        # Summary tables
        f.write("## Summary Statistics\n\n")
        f.write("### Samples per condition_harmonized\n\n")
        f.write("| Condition | N samples |\n|-----------|----------|\n")
        for k in sorted(cond_counts.keys()):
            f.write(f"| {k} | {cond_counts[k]} |\n")

        f.write("\n### Samples per compartment\n\n")
        f.write("| Compartment | N samples |\n|-------------|----------|\n")
        for k in sorted(comp_counts.keys()):
            f.write(f"| {k} | {comp_counts[k]} |\n")

        f.write("\n### Samples per age_group\n\n")
        f.write("| Age group | N samples |\n|-----------|----------|\n")
        for k in sorted(age_counts.keys()):
            f.write(f"| {k} | {age_counts[k]} |\n")

        f.write("\n### Donors per study\n\n")
        f.write("| Study | Donors | Samples |\n|-------|--------|--------|\n")
        for study in sorted(study_donors.keys()):
            n_donors = len(study_donors[study])
            n_samples = sum(1 for s in samples if s["study_accession"] == study)
            f.write(f"| {study} | {n_donors} | {n_samples} |\n")

        # Mapping decisions
        f.write("\n---\n\n")
        f.write("## Condition Harmonization Decisions\n\n")

        f.write("### General Rules\n\n")
        f.write("- **Pfirrmann I** → `healthy` (normal disc)\n")
        f.write("- **Pfirrmann II-III** → `degenerated_mild`\n")
        f.write("- **Pfirrmann IV-V** → `degenerated_severe`\n")
        f.write("- **Thompson I-II** → `healthy`\n")
        f.write("- **Thompson II-III** → `degenerated_mild`\n")
        f.write("- **Thompson III alone** → `degenerated_mild` (boundary; conservative choice)\n")
        f.write("- **Thompson III-IV or higher** → `degenerated_severe`\n")
        f.write("- **Herniated discs with grading** → `herniated` (herniation takes priority as condition)\n")
        f.write("- **Neonatal tissue** → `neonatal`\n")
        f.write("- **Aged without grading** → `aged_ungraded`\n\n")

        f.write("### Ambiguous Cases (Flagged for Human Review)\n\n")

        f.write("#### 1. GSE251686 (Jia 2024): Herniation + degeneration grading\n\n")
        f.write("All 6 samples described as having both degeneration (mild: Pfirrmann II-III, "
                "severe: Pfirrmann IV) AND herniation (lumbar disc herniation). "
                "**Decision:** classified as `herniated` with degeneration_severity "
                "recorded separately (mild or severe). The spec explicitly flags this "
                "case: 'Should they be categorized as degenerated_mild, herniated, "
                "or both?' We chose `herniated` as the primary condition.\n\n")

        f.write("#### 2. GSE189916 (Jiang 2022): Adult samples — aged but ungraded\n\n")
        f.write("Adult samples are from donors >65 years old with no back pain history. "
                "No degeneration grading available. **Decision:** classified as "
                "`aged_ungraded` per spec guidance ('Aged donors without explicit "
                "degeneration grading'). These are NOT classified as 'healthy' because "
                "no grading confirmation exists, nor as 'degenerated' because no "
                "degeneration was diagnosed.\n\n")

        f.write("#### 3. GSE205535 (Li 2022): 'Normal' from spinal cord injury patient\n\n")
        f.write("The 'normal' (NNP) sample is from an 11-year-old with acute spinal "
                "cord injury. The disc itself was not degenerated, but the clinical "
                "context (spinal trauma) is unusual. **Decision:** classified as "
                "`healthy` since the disc was not degenerated, but **FLAGGED** for "
                "human review. The 11-year-old's disc biology may differ from adult "
                "healthy discs due to developmental stage, and the spinal cord injury "
                "may have affected the local environment.\n\n")

        f.write("#### 4. GSE233666 (Guo 2023): IDD with herniation, has Pfirrmann grading\n\n")
        f.write("All 4 samples are from IDD patients undergoing lumbar discectomy for "
                "herniation. Pfirrmann grades II-III are available. **Decision:** "
                "classified as `herniated` because the primary clinical presentation "
                "and reason for surgery was disc herniation. Degeneration grades "
                "recorded in degeneration_grade_original and degeneration_severity.\n\n")

        f.write("#### 5. GSE244889 (Chen 2024): Author vs. spec classification conflict\n\n")
        f.write("Authors classify grade 1-2 as 'MDD' (mild disc degeneration) and "
                "grade 3-4 as 'SDD' (severe disc degeneration). However, our spec "
                "classifies Pfirrmann I as 'healthy.' **Decision:** followed the spec "
                "(Pfirrmann I = healthy, II-III = mild, IV-V = severe) rather than "
                "the authors' grouping. Two samples (Pa-17F, Pb-55F) with Pfirrmann 1 "
                "are classified as `healthy` despite authors calling them MDD.\n\n")

        f.write("#### 6. CNP0002664 (Han 2022): Ctrl sample with very low cell count\n\n")
        f.write("The control sample (Pfirrmann I, normal) has only ~249 cells reported. "
                "This is extremely low and may not survive QC filtering. "
                "**Decision:** included in metadata but flagged for QC review in "
                "Module 03.\n\n")

        # Grading discrepancies
        f.write("## Grading System Discrepancies\n\n")

        f.write("### GSE165722 GEO vs. Paper Grade Offset\n\n")
        f.write("GEO metadata lists Pfirrmann grades I-IV for the 8 samples, but the "
                "paper's Table 1 (PMC8787427) lists grades II-V. The paper explicitly "
                "states 'Pfirrmann grade I disc tissues were difficult to obtain.' "
                "**Decision:** used the paper's grading (II-V), which is authoritative. "
                "The GEO metadata has a systematic off-by-one error in Pfirrmann grades.\n\n")

        f.write("### Thompson III boundary classification\n\n")
        f.write("Thompson Grade III appears in both the 'mild' (II-III) and 'severe' "
                "(III-IV) ranges in the spec. For samples graded as Thompson III alone "
                "(not III-IV), we classified as `degenerated_mild` (conservative, "
                "lower boundary). For Thompson III-IV, we classified as "
                "`degenerated_severe`. This affects GSE230809 (Swahn 2024) samples "
                "SP21.011 and SP21.013 (Thompson III → mild) vs. SP20.002 etc. "
                "(Thompson III-IV → severe).\n\n")

        # Platform notes
        f.write("## Platform Notes\n\n")
        f.write("- **10x_3prime_v3**: 8 datasets (GSE160756, GSE189916, GSE233666, "
                "GSE244889, GSE255768, GSE230809, GSE242443; GSE199866 inferred as v2)\n")
        f.write("- **10x_3prime_v2**: GSE189916 (Chromium 3' v2), GSE199866 (inferred "
                "from HiSeq 4000 sequencer and CellRanger h5 format)\n")
        f.write("- **other (BD Rhapsody)**: GSE165722, GSE205535\n")
        f.write("- **other (Singleron Matrix)**: CNP0002664\n")
        f.write("- **other (Singleron GEXSCOPE)**: GSE251686 — NOTE: registry "
                "incorrectly lists as '10x Genomics'\n\n")

        # Missing data
        f.write("## Missing Metadata Summary\n\n")
        f.write("| Dataset | Age Missing | Sex Missing | Notes |\n")
        f.write("|---------|-------------|-------------|-------|\n")
        f.write("| GSE160756 | No | Yes | Sex in Supplementary Table 1 (not accessible) |\n")
        f.write("| GSE165722 | No | No | Full demographics from Table 1 |\n")
        f.write("| GSE189916 | Partial | Partial | Neonatal age=0, adult >65 (exact unknown); neonatal sex=M, adult unknown |\n")
        f.write("| GSE199866 | Yes | Yes | Single donor; demographics in Supplementary Table S26 (not accessible) |\n")
        f.write("| GSE205535 | No | Yes | Ages from paper (11y, 81y) |\n")
        f.write("| CNP0002664 | Yes | Yes | Demographics in Supplementary Table S1 (not accessible) |\n")
        f.write("| GSE233666 | No | No | Full demographics from Table 1 |\n")
        f.write("| GSE244889 | No | No | Encoded in sample names |\n")
        f.write("| GSE251686 | Yes | Yes | Not in paper text or GEO; may be in supplementary |\n")
        f.write("| GSE255768 | No | No | From paper text |\n")
        f.write("| GSE230809 | No | No | Excellent GEO metadata |\n")
        f.write("| GSE242443 | Yes | Yes | Not available in paper or GEO |\n")

        # Donor-level notes
        f.write("\n## Donor-Level Confounds\n\n")
        f.write("### Multi-sample donors\n\n")
        f.write("- **GSE160756 (Gan 2021):** Donor Gan_D02 (31y) contributes 6 of 7 "
                "samples (NP_2, NP_3, CEP_1, CEP_2, AF_1, AF_2). Donor Gan_D01 (18y) "
                "contributes only NP_1.\n")
        f.write("- **GSE189916 (Jiang 2022):** Donor Jiang_Neo01 contributes all 3 "
                "neonatal samples (different spinal levels from same donor).\n")
        f.write("- **GSE199866 (Cherif 2022):** Single donor (Cherif_D01) contributes "
                "all 4 samples (paired healthy/degenerated × NP/AF).\n")
        f.write("- **GSE230809 (Swahn 2024):** 8 of 11 donors contribute paired NP+AF "
                "samples. 2 donors (SP20.002, SP20.006) contribute AF only.\n\n")

        f.write("### Age-disease confounding\n\n")
        f.write("- **GSE230809 (Swahn 2024):** Healthy donors are 21-27 years old; "
                "diseased donors are 37-73 years old. This is a strong age-disease "
                "confound that must be addressed during integration.\n")
        f.write("- **GSE205535 (Li 2022):** Normal donor is 11y, degenerative is 81y. "
                "Extreme age difference.\n\n")

        f.write("### Sex bias\n\n")
        f.write("- **GSE230809 (Swahn 2024):** ALL 11 donors are male. This is the "
                "largest dataset (24 samples) and will dominate sex-unaware analyses.\n")
        f.write("- Most other datasets have no sex information, limiting "
                "sex-stratified analyses.\n\n")

        # Culture expansion
        f.write("## Special Considerations\n\n")
        f.write("### Culture-expanded cells (GSE242443)\n\n")
        f.write("Both samples from Kuchynsky 2024 are culture-expanded CEP cells "
                "(passage 1). Gene expression may be altered by in vitro conditions. "
                "Included per human decision at Module 01 checkpoint to improve "
                "sparse CEP coverage. tissue_or_cells field set to 'cells' to "
                "distinguish from fresh tissue.\n\n")

        f.write("### Normalized vs. raw counts (GSE165722)\n\n")
        f.write("GEO supplementary files for Tu 2022 may contain normalized counts "
                "rather than raw UMI counts. This must be verified during Module 03 "
                "preprocessing. If only normalized counts are available, raw data may "
                "need to be obtained from SRA.\n\n")

        f.write("### Published corrections (GSE205535)\n\n")
        f.write("Li Z 2022 has published corrections/corrigenda. These should be "
                "reviewed during preprocessing to determine if any data corrections "
                "are needed.\n\n")

    print(f"Wrote harmonization notes to {outpath}")


def validate(samples):
    """Run automated validation checks."""
    errors = []

    # Check all required fields populated
    for s in samples:
        for field in FIELDS:
            if field not in s or s[field] == "" or s[field] is None:
                errors.append(f"Empty field '{field}' in sample {s.get('sample_id', '?')}")

    # Check no duplicate sample_ids
    ids = [s["sample_id"] for s in samples]
    if len(ids) != len(set(ids)):
        dupes = [x for x in ids if ids.count(x) > 1]
        errors.append(f"Duplicate sample_ids: {set(dupes)}")

    # Check allowed values
    for s in samples:
        if s["condition_harmonized"] not in ALLOWED_CONDITIONS:
            errors.append(f"Invalid condition_harmonized '{s['condition_harmonized']}' "
                          f"in {s['sample_id']}")
        if s["compartment"] not in ALLOWED_COMPARTMENTS:
            errors.append(f"Invalid compartment '{s['compartment']}' in {s['sample_id']}")
        if s["sex"] not in ALLOWED_SEX:
            errors.append(f"Invalid sex '{s['sex']}' in {s['sample_id']}")
        if s["age_group"] not in ALLOWED_AGE_GROUPS:
            errors.append(f"Invalid age_group '{s['age_group']}' in {s['sample_id']}")

    # Check species
    for s in samples:
        if s["species"] != "human":
            errors.append(f"Non-human species '{s['species']}' in {s['sample_id']}")

    # Check expected sample counts per study
    expected_counts = {
        "GSE160756": 7, "GSE165722": 8, "GSE189916": 6, "GSE199866": 4,
        "GSE205535": 2, "CNP0002664": 6, "GSE233666": 4, "GSE244889": 7,
        "GSE251686": 6, "GSE255768": 2, "GSE230809": 24, "GSE242443": 2,
    }
    actual_counts = {}
    for s in samples:
        acc = s["study_accession"]
        actual_counts[acc] = actual_counts.get(acc, 0) + 1
    for acc, expected in expected_counts.items():
        actual = actual_counts.get(acc, 0)
        if actual != expected:
            errors.append(f"Sample count mismatch for {acc}: expected {expected}, got {actual}")

    # Check all included studies present
    for acc in expected_counts:
        if acc not in actual_counts:
            errors.append(f"Study {acc} missing from metadata")

    return errors


if __name__ == "__main__":
    samples = build_samples()

    # Validate before writing
    errors = validate(samples)
    if errors:
        print("VALIDATION ERRORS:")
        for e in errors:
            print(f"  - {e}")
        print(f"\n{len(errors)} error(s) found. Writing files anyway for inspection.")
    else:
        print("All validation checks passed!")

    # Write outputs
    write_sample_metadata(samples, METADATA_DIR / "sample_metadata.tsv")
    write_ontology_mappings(METADATA_DIR / "ontology_mappings.tsv")
    write_harmonization_notes(samples, METADATA_DIR / "harmonization_notes.md")

    print(f"\nSummary: {len(samples)} samples across "
          f"{len(set(s['study_accession'] for s in samples))} studies")
    print("\nCondition distribution:")
    cond_counts = {}
    for s in samples:
        c = s["condition_harmonized"]
        cond_counts[c] = cond_counts.get(c, 0) + 1
    for k, v in sorted(cond_counts.items()):
        print(f"  {k}: {v}")
