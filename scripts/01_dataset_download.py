#!/usr/bin/env python3
"""
Module 01: Dataset Discovery & Acquisition - Download Script

Downloads raw count matrices and metadata for all included human IVD
single-cell RNA-seq datasets from GEO.

Datasets in Chinese repositories (CNGB, NGDC/GSA) are handled separately
as they require different access procedures.

Usage:
    python scripts/01_dataset_download.py
"""

import os
import subprocess
import hashlib
import json
from pathlib import Path
from datetime import datetime

BASE_DIR = Path(__file__).resolve().parent.parent
RAW_DIR = BASE_DIR / "data" / "raw"
METADATA_DIR = BASE_DIR / "metadata"

# GEO datasets to download (accession -> metadata)
GEO_DATASETS = {
    "GSE160756": {
        "first_author": "Gan Y",
        "year": 2021,
        "compartment": "NP, AF, CEP",
        "description": "Healthy IVD atlas - all 3 compartments",
    },
    "GSE165722": {
        "first_author": "Tu J",
        "year": 2022,
        "compartment": "NP",
        "description": "NP progressive degeneration (Pfirrmann II-V)",
    },
    "GSE189916": {
        "first_author": "Jiang W",
        "year": 2022,
        "compartment": "Whole IVD",
        "description": "Neonatal vs adult IVD",
    },
    "GSE199866": {
        "first_author": "Cherif H",
        "year": 2022,
        "compartment": "NP, inner AF",
        "description": "Paired degen vs non-degen from same individual",
    },
    "GSE205535": {
        "first_author": "Li Z",
        "year": 2022,
        "compartment": "NP",
        "description": "Normal vs degenerative NP",
    },
    "GSE233666": {
        "first_author": "Guo S",
        "year": 2023,
        "compartment": "NP",
        "description": "NP from IDD patients - immune involvement",
    },
    "GSE244889": {
        "first_author": "Chen F",
        "year": 2024,
        "compartment": "NP",
        "description": "NP mild vs severe degeneration - serglycin biomarker",
    },
    "GSE251686": {
        "first_author": "Jia S",
        "year": 2024,
        "compartment": "NP",
        "description": "NP mild vs severe degeneration",
    },
    "GSE255768": {
        "first_author": "Shi C",
        "year": 2024,
        "compartment": "CEP/Endplate",
        "description": "Degenerative endplate (Modic changes)",
    },
    "GSE230809": {
        "first_author": "Swahn H",
        "year": 2024,
        "compartment": "NP, AF",
        "description": "SuperSeries: healthy + diseased NP and AF",
    },
    "GSE242443": {
        "first_author": "Kuchynsky K",
        "year": 2024,
        "compartment": "CEP",
        "description": "CEP non-degen vs degen (BORDERLINE: culture-expanded)",
    },
}


def md5sum(filepath):
    """Compute MD5 checksum of a file."""
    h = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def get_geo_file_list(accession):
    """Parse the GEO FTP directory listing to get supplementary file URLs."""
    import re
    import urllib.request

    prefix = accession[:-3] + "nnn"
    base_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{accession}/suppl/"

    try:
        req = urllib.request.Request(base_url)
        with urllib.request.urlopen(req, timeout=30) as resp:
            html = resp.read().decode("utf-8")
    except Exception as e:
        print(f"  ERROR fetching directory listing: {e}")
        return []

    # Parse href links from the HTML directory listing
    # Exclude parent directory and filelist.txt
    files = re.findall(r'href="([^"]+)"', html)
    data_files = [
        f for f in files
        if not f.startswith("/")
        and not f.startswith("http")
        and f != "filelist.txt"
    ]
    return [(f, base_url + f) for f in data_files]


def download_geo_supplementary(accession, output_dir):
    """Download supplementary files from GEO."""
    output_dir.mkdir(parents=True, exist_ok=True)

    prefix = accession[:-3] + "nnn"
    base_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{accession}/suppl/"

    print(f"\n{'='*60}")
    print(f"Downloading {accession}")
    print(f"  URL: {base_url}")
    print(f"  Output: {output_dir}")
    print(f"{'='*60}")

    file_list = get_geo_file_list(accession)
    if not file_list:
        print(f"  No files found in directory listing for {accession}")
        return False

    print(f"  Found {len(file_list)} file(s): {[f[0] for f in file_list]}")

    success = True
    for fname, url in file_list:
        outpath = output_dir / fname
        print(f"  Downloading {fname}...")
        cmd = [
            "wget", "-q", "--show-progress", "--no-check-certificate",
            "-O", str(outpath),
            url,
        ]
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
            if result.returncode != 0:
                print(f"    FAILED: {result.stderr[:200]}")
                success = False
            else:
                size = outpath.stat().st_size
                if size < 1000:
                    print(f"    WARNING: file suspiciously small ({size} bytes)")
                else:
                    size_mb = size / (1024 * 1024)
                    print(f"    OK ({size_mb:.1f} MB)")
        except subprocess.TimeoutExpired:
            print(f"    TIMEOUT downloading {fname}")
            success = False
        except Exception as e:
            print(f"    ERROR: {e}")
            success = False

    return success


def compute_checksums(directory):
    """Compute MD5 checksums for all files in directory."""
    checksums = {}
    for f in sorted(directory.rglob("*")):
        if f.is_file():
            checksums[str(f.relative_to(directory))] = md5sum(f)
    return checksums


def write_readme(accession, info, output_dir, checksums):
    """Write a README for each downloaded dataset."""
    readme_path = output_dir / "README.md"
    with open(readme_path, "w") as f:
        f.write(f"# {accession} - {info['first_author']} ({info['year']})\n\n")
        f.write(f"**Compartment:** {info['compartment']}\n\n")
        f.write(f"**Description:** {info['description']}\n\n")
        f.write(f"**Download date:** {datetime.now().strftime('%Y-%m-%d')}\n\n")
        f.write(f"**Source:** GEO ({accession})\n\n")
        f.write("## File Checksums (MD5)\n\n")
        for fname, checksum in checksums.items():
            if fname != "README.md":
                f.write(f"- `{fname}`: `{checksum}`\n")
        f.write("\n")


def main():
    print("=" * 60)
    print("IVD Single-Cell Atlas - Dataset Download")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    RAW_DIR.mkdir(parents=True, exist_ok=True)

    results = {}
    all_checksums = {}

    for accession, info in GEO_DATASETS.items():
        output_dir = RAW_DIR / accession
        if output_dir.exists() and any(output_dir.iterdir()):
            print(f"\n{accession}: Already downloaded, skipping.")
            # Still compute checksums for existing files
            checksums = compute_checksums(output_dir)
            all_checksums[accession] = checksums
            results[accession] = "already_present"
            continue

        success = download_geo_supplementary(accession, output_dir)
        if success:
            checksums = compute_checksums(output_dir)
            all_checksums[accession] = checksums
            write_readme(accession, info, output_dir, checksums)
            results[accession] = "downloaded"
            print(f"  SUCCESS: {len(checksums)} files downloaded")
        else:
            results[accession] = "failed"

    # Save checksums to metadata
    checksums_path = METADATA_DIR / "file_checksums.json"
    with open(checksums_path, "w") as f:
        json.dump(all_checksums, f, indent=2)
    print(f"\nChecksums saved to {checksums_path}")

    # Print summary
    print("\n" + "=" * 60)
    print("DOWNLOAD SUMMARY")
    print("=" * 60)
    for accession, status in results.items():
        info = GEO_DATASETS[accession]
        print(f"  {accession} ({info['first_author']} {info['year']}): {status}")

    failed = [a for a, s in results.items() if s == "failed"]
    if failed:
        print(f"\nWARNING: {len(failed)} dataset(s) failed to download: {', '.join(failed)}")
        print("These may need manual download or alternative access methods.")

    # Note about non-GEO datasets
    print("\n" + "=" * 60)
    print("NON-GEO DATASETS (require separate download)")
    print("=" * 60)
    print("  CNP0002664 (Han S 2022) - China National GeneBank (CNGB)")
    print("    URL: https://db.cngb.org/search/project/CNP0002664/")
    print("  PRJCA014236 (Wang D 2023) - Genome Sequence Archive")
    print("    URL: https://ngdc.cncb.ac.cn/gsa-human/browse/HRA004080")
    print("  PRJCA007656 (Ling Z 2022) - NGDC BioProject")
    print("    URL: https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA007656")
    print("\nThese datasets are in Chinese national repositories and may require")
    print("registration or special access procedures.")


if __name__ == "__main__":
    main()
