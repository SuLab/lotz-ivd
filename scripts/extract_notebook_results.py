#!/usr/bin/env python3
"""Extract key result tables from executed notebook outputs.

Parses .ipynb JSON to recover structured data from cell outputs, writing
them as TSVs to docs/v5_results/. This allows 12_reporting.py to generate
a complete report from committed files alone, without needing results/.

Usage:
    python3 scripts/extract_notebook_results.py
"""

import json
import re
import io
from pathlib import Path
from datetime import datetime

import pandas as pd

BASE = Path(__file__).resolve().parent.parent
NOTEBOOKS = BASE / "notebooks"
OUTDIR = BASE / "docs" / "v5_results"
OUTDIR.mkdir(parents=True, exist_ok=True)


def load_notebook(name):
    """Load a notebook and return its cells."""
    path = NOTEBOOKS / name
    if not path.exists():
        print(f"  SKIP: {name} not found")
        return []
    with open(path) as f:
        nb = json.load(f)
    return nb.get("cells", [])


def get_stream_text(cell):
    """Extract all stream (stdout) text from a cell's outputs."""
    parts = []
    for out in cell.get("outputs", []):
        if out.get("output_type") == "stream":
            parts.extend(out.get("text", []))
    return "".join(parts)


def get_html_tables(cell):
    """Extract HTML table content from display_data outputs."""
    tables = []
    for out in cell.get("outputs", []):
        if out.get("output_type") in ("display_data", "execute_result"):
            html = "".join(out.get("data", {}).get("text/html", []))
            if "<table" in html.lower():
                try:
                    dfs = pd.read_html(io.StringIO(html))
                    tables.extend(dfs)
                except Exception:
                    pass
    return tables


def get_text_plain(cell):
    """Extract text/plain from display_data outputs."""
    parts = []
    for out in cell.get("outputs", []):
        if out.get("output_type") in ("display_data", "execute_result"):
            txt = "".join(out.get("data", {}).get("text/plain", []))
            if txt:
                parts.append(txt)
    return parts


def _parse_interlaced_compartment_counts(cell, count_pattern):
    """Parse outputs where compartment labels (in markdown display_data) and
    counts (in stream text) are interlaced across separate outputs."""
    results = {}
    current_comp = None
    for out in cell.get("outputs", []):
        if out.get("output_type") == "display_data":
            md = "".join(out.get("data", {}).get("text/markdown", []))
            for comp in ["NP", "AF", "CEP"]:
                if comp in md:
                    current_comp = comp
                    break
        elif out.get("output_type") == "stream":
            text = "".join(out.get("text", []))
            m = re.search(count_pattern, text)
            if m and current_comp:
                results[current_comp] = int(m.group(1))
                current_comp = None
    return results


# ---------------------------------------------------------------------------
# Extractors — one per notebook
# ---------------------------------------------------------------------------

def extract_08_differential(cells):
    """Extract DE summary, top genes, skipped comparisons, composition."""
    extracted = {}

    # Cell 3: DE summary table (HTML Styler)
    if len(cells) > 3:
        tables = get_html_tables(cells[3])
        if tables:
            extracted["de_summary_table"] = tables[0]

    # Cell 4: stream has aggregate counts; after fix, will also have top genes
    # For now, parse stream for gene counts
    if len(cells) > 4:
        stream = get_stream_text(cells[4])
        # Try to get top gene tables from HTML outputs (added by our fix)
        tables = get_html_tables(cells[4])
        if len(tables) >= 2:
            extracted["de_top_upregulated"] = tables[0]
            extracted["de_top_downregulated"] = tables[1]
        elif len(tables) == 1:
            extracted["de_top_upregulated"] = tables[0]

    # Cell 6: skipped comparisons
    if len(cells) > 6:
        tables = get_html_tables(cells[6])
        # Also try text/plain (skipped is displayed as plain DataFrame)
        if tables:
            extracted["skipped_comparisons"] = tables[0]
        else:
            text_parts = get_text_plain(cells[6])
            for txt in text_parts:
                if "comparison" in txt and "reason" in txt:
                    try:
                        df = pd.read_csv(io.StringIO(txt), sep=r"\s{2,}", engine="python")
                        extracted["skipped_comparisons"] = df
                    except Exception:
                        pass

    # Cell 16: composition analysis
    if len(cells) > 16:
        tables = get_html_tables(cells[16])
        if tables:
            extracted["composition_analysis"] = tables[0]

    return extracted


def extract_09_interpretation(cells):
    """Extract ORA, GSEA, TF, pain gene results."""
    extracted = {}

    # Cell 3: ORA summary + top terms
    if len(cells) > 3:
        stream = get_stream_text(cells[3])
        m = re.search(r"Total enrichment results:\s*(\d+)", stream)
        if m:
            extracted["_ora_total"] = int(m.group(1))

        tables = get_html_tables(cells[3])
        if len(tables) >= 1:
            extracted["ora_summary"] = tables[0]
        if len(tables) >= 2:
            extracted["ora_top_terms"] = tables[1]

    # Cell 9: GSEA
    if len(cells) > 9:
        stream = get_stream_text(cells[9])
        m = re.search(r"Total GSEA terms tested:\s*(\d+)", stream)
        if m:
            extracted["_gsea_total"] = int(m.group(1))
        m = re.search(r"Significant GSEA terms.*?:\s*(\d+)", stream)
        if m:
            extracted["_gsea_significant"] = int(m.group(1))

        tables = get_html_tables(cells[9])
        if tables:
            extracted["gsea_top"] = tables[0]

    # Cell 12: TF activity
    if len(cells) > 12:
        stream = get_stream_text(cells[12])
        m = re.search(r"Significant TFs.*?:\s*(\d+)", stream)
        if m:
            extracted["_tf_significant"] = int(m.group(1))

        tables = get_html_tables(cells[12])
        if tables:
            extracted["tf_activity_top"] = tables[0]

    # Cell 15: Pain genes
    if len(cells) > 15:
        stream = get_stream_text(cells[15])
        m = re.search(r"Significant pain genes:\s*(\d+)", stream)
        if m:
            extracted["_pain_significant"] = int(m.group(1))

        tables = get_html_tables(cells[15])
        if tables:
            extracted["pain_genes"] = tables[0]

    return extracted


def extract_10_trajectory(cells):
    """Extract pseudotime correlations and trajectory gene counts."""
    extracted = {}

    # Cell 3: pseudotime correlations table
    if len(cells) > 3:
        tables = get_html_tables(cells[3])
        if tables:
            extracted["pseudotime_correlations"] = tables[0]

    # Cell 13: trajectory gene counts — outputs are interlaced:
    #   display_data (markdown: "### NP — ...") → stream ("Total trajectory genes: 500")
    if len(cells) > 13:
        gene_counts = _parse_interlaced_compartment_counts(
            cells[13], r"Total trajectory genes[=:]\s*(\d+)")
        if gene_counts:
            extracted["trajectory_gene_counts"] = gene_counts

    # Cell 15: DE overlap counts — same interlaced pattern
    if len(cells) > 15:
        overlaps = _parse_interlaced_compartment_counts(
            cells[15], r"Overlapping genes[=:]\s*(\d+)")
        if overlaps:
            extracted["trajectory_de_overlap"] = overlaps

    return extracted


def extract_11_communication(cells):
    """Extract CCC interaction counts and pain interactions."""
    extracted = {}

    # Cell 3: interaction counts
    if len(cells) > 3:
        stream = get_stream_text(cells[3])
        m = re.search(r"Healthy interactions:\s*([\d,]+)", stream)
        if m:
            extracted["_ccc_healthy"] = int(m.group(1).replace(",", ""))
        m = re.search(r"Degenerated interactions:\s*([\d,]+)", stream)
        if m:
            extracted["_ccc_degenerated"] = int(m.group(1).replace(",", ""))

    # Cell 10: differential interactions
    if len(cells) > 10:
        stream = get_stream_text(cells[10])
        m = re.search(r"Total differential interactions\s*=\s*([\d,]+)", stream)
        if m:
            extracted["_ccc_differential"] = int(m.group(1).replace(",", ""))

    # Cell 13: pain interactions
    if len(cells) > 13:
        stream = get_stream_text(cells[13])
        m = re.search(r"Total pain-relevant interactions:\s*([\d,]+)", stream)
        if m:
            extracted["_pain_interactions"] = int(m.group(1).replace(",", ""))

        # Parse category breakdown
        categories = {}
        for m in re.finditer(r"(pain_\w+:\w+)\s*:\s*(\d+)", stream):
            categories[m.group(1)] = int(m.group(2))
        if categories:
            extracted["pain_interaction_categories"] = categories

    return extracted


def extract_07_annotation(cells):
    """Extract cell type definitions and CellTypist concordance."""
    extracted = {}

    # Cell 3: cell type definitions table
    if len(cells) > 3:
        tables = get_html_tables(cells[3])
        if tables:
            extracted["cell_type_definitions"] = tables[0]

    # Cell 9: cell type counts (from stream)
    if len(cells) > 9:
        stream = get_stream_text(cells[9])
        extracted["_annotation_stream"] = stream

    # Cell 11: CellTypist concordance tables (text/plain DataFrames)
    if len(cells) > 11:
        text_parts = get_text_plain(cells[11])
        stream = get_stream_text(cells[11])
        # Combine for parsing
        all_text = stream + "\n".join(text_parts)
        extracted["_celltypist_text"] = all_text

        tables = get_html_tables(cells[11])
        if tables:
            for i, t in enumerate(tables):
                extracted[f"celltypist_concordance_{i}"] = t

    return extracted


def extract_06_clustering(cells):
    """Extract clustering resolution optimization tables."""
    extracted = {}

    # Cell 5: resolution scan tables (multiple HTML Styler tables)
    if len(cells) > 5:
        tables = get_html_tables(cells[5])
        if tables:
            extracted["resolution_tables"] = tables

    return extracted


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 60)
    print("Extracting results from notebook outputs")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Output: {OUTDIR}")
    print("=" * 60)

    written = 0

    # ── 06 Clustering ──
    print("\n06_clustering.ipynb...")
    cells = load_notebook("06_clustering.ipynb")
    if cells:
        data = extract_06_clustering(cells)
        if "resolution_tables" in data:
            combined = pd.concat(data["resolution_tables"], ignore_index=True)
            combined.to_csv(OUTDIR / "clustering_resolutions.tsv", sep="\t", index=False)
            print(f"  clustering_resolutions.tsv ({len(combined)} rows)")
            written += 1

    # ── 07 Annotation ──
    print("\n07_annotation.ipynb...")
    cells = load_notebook("07_annotation.ipynb")
    if cells:
        data = extract_07_annotation(cells)
        if "cell_type_definitions" in data:
            data["cell_type_definitions"].to_csv(
                OUTDIR / "cell_type_definitions.tsv", sep="\t", index=False)
            print(f"  cell_type_definitions.tsv ({len(data['cell_type_definitions'])} rows)")
            written += 1
        for key in sorted(data):
            if key.startswith("celltypist_concordance_"):
                df = data[key]
                df.to_csv(OUTDIR / f"{key}.tsv", sep="\t", index=False)
                print(f"  {key}.tsv ({len(df)} rows)")
                written += 1

    # ── 08 Differential ──
    print("\n08_differential.ipynb...")
    cells = load_notebook("08_differential.ipynb")
    if cells:
        data = extract_08_differential(cells)
        for key in ["de_summary_table", "de_top_upregulated", "de_top_downregulated",
                     "skipped_comparisons", "composition_analysis"]:
            if key in data:
                df = data[key]
                df.to_csv(OUTDIR / f"{key}.tsv", sep="\t", index=False)
                print(f"  {key}.tsv ({len(df)} rows)")
                written += 1

    # ── 09 Interpretation ──
    print("\n09_interpretation.ipynb...")
    cells = load_notebook("09_interpretation.ipynb")
    if cells:
        data = extract_09_interpretation(cells)
        for key in ["ora_summary", "ora_top_terms", "gsea_top",
                     "tf_activity_top", "pain_genes"]:
            if key in data:
                df = data[key]
                df.to_csv(OUTDIR / f"{key}.tsv", sep="\t", index=False)
                print(f"  {key}.tsv ({len(df)} rows)")
                written += 1

        # Write scalar stats as a simple key-value TSV
        scalars = {k: v for k, v in data.items() if k.startswith("_")}
        if scalars:
            pd.DataFrame([scalars]).to_csv(
                OUTDIR / "interpretation_stats.tsv", sep="\t", index=False)
            print(f"  interpretation_stats.tsv ({len(scalars)} scalars)")
            written += 1

    # ── 10 Trajectory ──
    print("\n10_trajectory.ipynb...")
    cells = load_notebook("10_trajectory.ipynb")
    if cells:
        data = extract_10_trajectory(cells)
        if "pseudotime_correlations" in data:
            data["pseudotime_correlations"].to_csv(
                OUTDIR / "pseudotime_correlations.tsv", sep="\t", index=False)
            print(f"  pseudotime_correlations.tsv ({len(data['pseudotime_correlations'])} rows)")
            written += 1
        if "trajectory_gene_counts" in data:
            df = pd.DataFrame(
                [{"compartment": k, "n_genes": v}
                 for k, v in data["trajectory_gene_counts"].items()])
            df.to_csv(OUTDIR / "trajectory_gene_counts.tsv", sep="\t", index=False)
            print(f"  trajectory_gene_counts.tsv ({len(df)} rows)")
            written += 1
        if "trajectory_de_overlap" in data:
            df = pd.DataFrame(
                [{"compartment": k, "n_overlap": v}
                 for k, v in data["trajectory_de_overlap"].items()])
            df.to_csv(OUTDIR / "trajectory_de_overlap.tsv", sep="\t", index=False)
            print(f"  trajectory_de_overlap.tsv ({len(df)} rows)")
            written += 1

    # ── 11 Communication ──
    print("\n11_communication.ipynb...")
    cells = load_notebook("11_communication.ipynb")
    if cells:
        data = extract_11_communication(cells)
        scalars = {k: v for k, v in data.items()
                   if isinstance(v, (int, float, str))}
        if scalars:
            pd.DataFrame([scalars]).to_csv(
                OUTDIR / "communication_stats.tsv", sep="\t", index=False)
            print(f"  communication_stats.tsv ({len(scalars)} scalars)")
            written += 1
        if "pain_interaction_categories" in data:
            df = pd.DataFrame(
                [{"category": k, "count": v}
                 for k, v in data["pain_interaction_categories"].items()])
            df.to_csv(OUTDIR / "pain_interaction_categories.tsv", sep="\t", index=False)
            print(f"  pain_interaction_categories.tsv ({len(df)} rows)")
            written += 1

    print(f"\n{'=' * 60}")
    print(f"Done. {written} files written to {OUTDIR}")


if __name__ == "__main__":
    main()
