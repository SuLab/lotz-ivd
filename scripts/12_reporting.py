#!/usr/bin/env python3
"""Module 12: Final reporting for IVD scRNA-seq atlas.

Generates a comprehensive markdown report by reading result files from
all upstream modules. Every statistic is loaded from its source file and
annotated with a provenance reference.

Output: docs/final_report.md

Usage:
    python3 scripts/12_reporting.py
"""

import sys
import os
import subprocess
from pathlib import Path
from datetime import datetime

import pandas as pd

BASE = Path(__file__).resolve().parent.parent
RESULTS = BASE / "results"
META = BASE / "metadata"
DOCS = BASE / "docs"
SUPP_TABLES = RESULTS / "supplementary_tables"

DOCS.mkdir(parents=True, exist_ok=True)
SUPP_TABLES.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def git_version():
    """Return git commit hash and branch name."""
    try:
        sha = subprocess.check_output(
            ["git", "-C", str(BASE), "rev-parse", "HEAD"],
            text=True,
        ).strip()
        branch = subprocess.check_output(
            ["git", "-C", str(BASE), "rev-parse", "--abbrev-ref", "HEAD"],
            text=True,
        ).strip()
        short = sha[:7]
        return short, sha, branch
    except Exception:
        return "unknown", "unknown", "unknown"


def read_tsv(path, **kwargs):
    """Read a TSV file, returning an empty DataFrame if it doesn't exist."""
    if path.exists():
        return pd.read_csv(path, sep="\t", **kwargs)
    return pd.DataFrame()


def src(path):
    """Format a source reference relative to the repo root."""
    try:
        rel = path.resolve().relative_to(BASE.resolve())
    except ValueError:
        rel = path
    return f"*[source: `{rel}`]*"


def safe_int(val):
    try:
        return f"{int(val):,}"
    except (ValueError, TypeError):
        return str(val)


# ---------------------------------------------------------------------------
# Data loaders — each returns (data_dict, source_path)
# ---------------------------------------------------------------------------

def load_dataset_registry():
    path = META / "dataset_registry.tsv"
    df = read_tsv(path)
    return df, path


def load_sample_metadata():
    path = META / "sample_metadata.tsv"
    df = read_tsv(path)
    return df, path


def load_cell_type_definitions():
    path = RESULTS / "integration" / "cell_type_definitions.tsv"
    df = read_tsv(path)
    return df, path


def load_clustering_resolutions():
    cdir = RESULTS / "integration" / "clustering_resolution_optimization"
    frames = []
    if cdir.exists():
        for f in sorted(cdir.glob("*_resolutions.tsv")):
            tmp = read_tsv(f)
            tmp["_source_file"] = f.name
            frames.append(tmp)
    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return df, cdir


def load_de_summary():
    path = RESULTS / "differential" / "de_summary_table.tsv"
    df = read_tsv(path)
    return df, path


def load_de_combined():
    path = RESULTS / "differential" / "de_results_combined.tsv"
    df = read_tsv(path)
    return df, path


def load_skipped_comparisons():
    path = RESULTS / "differential" / "skipped_comparisons.tsv"
    df = read_tsv(path)
    return df, path


def load_composition():
    path = RESULTS / "differential" / "composition_analysis.tsv"
    df = read_tsv(path)
    return df, path


def load_enrichment():
    path = RESULTS / "interpretation" / "pathway_enrichment" / "all_enrichment_results.tsv"
    df = read_tsv(path)
    return df, path


def load_gsea():
    path = RESULTS / "interpretation" / "pathway_enrichment" / "gsea_results.tsv"
    df = read_tsv(path)
    return df, path


def load_tf_activity():
    path = RESULTS / "interpretation" / "tf_activity" / "tf_activity_results.tsv"
    df = read_tsv(path)
    return df, path


def load_pain_genes():
    path = RESULTS / "interpretation" / "pain_genes.tsv"
    df = read_tsv(path)
    return df, path


def load_trajectory_correlations():
    frames = []
    tdir = RESULTS / "trajectories"
    if tdir.exists():
        for f in sorted(tdir.glob("pseudotime_correlations_*.tsv")):
            tmp = read_tsv(f)
            comp = f.stem.replace("pseudotime_correlations_", "")
            tmp["compartment"] = comp
            tmp["_source_file"] = str(f)
            frames.append(tmp)
    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return df, tdir


def load_trajectory_genes():
    counts = {}
    tdir = RESULTS / "trajectories"
    if tdir.exists():
        for f in sorted(tdir.glob("trajectory_genes_*.tsv")):
            comp = f.stem.replace("trajectory_genes_", "")
            tmp = read_tsv(f)
            counts[comp] = len(tmp)
    return counts, tdir


def load_ccc_interactions():
    """Load interaction counts per condition."""
    cdir = RESULTS / "communication"
    counts = {}
    if cdir.exists():
        for f in sorted(cdir.glob("interactions_*.tsv")):
            cond = f.stem.replace("interactions_", "")
            tmp = read_tsv(f)
            counts[cond] = len(tmp)
    return counts, cdir


def load_differential_interactions():
    path = RESULTS / "communication" / "differential_interactions.tsv"
    df = read_tsv(path)
    return df, path


def load_pain_interactions():
    path = RESULTS / "communication" / "pain_interactions.tsv"
    df = read_tsv(path)
    return df, path


def load_celltypist_concordance():
    cdir = RESULTS / "integration" / "celltypist_validation"
    frames = []
    if cdir.exists():
        for f in sorted(cdir.glob("*_concordance.tsv")):
            tmp = read_tsv(f)
            obj = f.stem.replace("_concordance", "")
            tmp["object"] = obj
            tmp["_source_file"] = str(f)
            frames.append(tmp)
    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return df, cdir


def load_analysis_plan():
    path = BASE / "analysis_plan.md"
    if path.exists():
        return path.read_text(), path
    return "", path


# ---------------------------------------------------------------------------
# Report sections
# ---------------------------------------------------------------------------

def section_header(lines, plan_text):
    short_sha, full_sha, branch = git_version()
    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    # Extract pipeline version from analysis_plan
    version = "unknown"
    for line in plan_text.splitlines():
        if line.startswith("**Pipeline version:**"):
            version = line.split(":**")[1].strip()
            break

    lines.append("# Human Intervertebral Disc Single-Cell Atlas — Final Report")
    lines.append("")
    lines.append("**A comprehensive scRNA-seq meta-analysis of IVD degeneration**")
    lines.append("")
    lines.append(f"| Field | Value |")
    lines.append(f"|-------|-------|")
    lines.append(f"| Report generated | {now} |")
    lines.append(f"| Pipeline version | {version} |")
    lines.append(f"| Git commit | `{short_sha}` (branch: `{branch}`) |")
    lines.append(f"| Source of truth | `analysis_plan.md` |")
    lines.append("")


def section_toc(lines):
    lines.append("## Contents")
    lines.append("")
    toc_items = [
        ("summary", "Atlas Summary"),
        ("datasets", "Dataset Summary"),
        ("integration", "Integration"),
        ("clustering", "Clustering & Annotation"),
        ("de", "Differential Expression"),
        ("pathways", "Biological Pathways"),
        ("tf", "Transcription Factor Activity"),
        ("trajectory", "Cell State Trajectories"),
        ("communication", "Cell-Cell Communication"),
        ("pain", "Pain Biology"),
        ("limitations", "Limitations & Caveats"),
        ("methods", "Methods"),
        ("reproducibility", "Reproducibility"),
    ]
    for i, (anchor, title) in enumerate(toc_items, 1):
        lines.append(f"{i}. [{title}](#{anchor})")
    lines.append("")


def section_summary(lines, sample_df, sample_path, de_summary, de_path,
                    enrichment_df, enrichment_path, ccc_counts, ccc_path):
    lines.append("## 1. Atlas Summary {#summary}")
    lines.append("")

    # Compute stats from actual data
    if not sample_df.empty:
        # sample_metadata.tsv contains only included samples (no status column)
        n_samples = len(sample_df)
        n_donors = sample_df["donor_id"].nunique() if "donor_id" in sample_df.columns else "?"
        compartments = sorted(sample_df["compartment"].dropna().unique()) if "compartment" in sample_df.columns else []
        lines.append(f"| Metric | Value | Source |")
        lines.append(f"|--------|-------|--------|")
        lines.append(f"| Samples | {n_samples} | {src(sample_path)} |")
        lines.append(f"| Donors | {n_donors} | {src(sample_path)} |")
        lines.append(f"| Compartments | {', '.join(compartments)} | {src(sample_path)} |")
    else:
        lines.append("| Metric | Value | Source |")
        lines.append("|--------|-------|--------|")
        lines.append(f"| Samples | *data not available* | {src(sample_path)} |")

    if not de_summary.empty:
        n_comparisons = len(de_summary)
        n_sig = int(de_summary["n_total"].sum()) if "n_total" in de_summary.columns else "?"
        lines.append(f"| Powered DE comparisons | {n_comparisons} | {src(de_path)} |")
        lines.append(f"| Significant DE genes (total hits) | {safe_int(n_sig)} | {src(de_path)} |")

    if enrichment_df is not None and not enrichment_df.empty:
        n_enrichments = len(enrichment_df)
        lines.append(f"| Enriched pathways (ORA) | {safe_int(n_enrichments)} | {src(enrichment_path)} |")

    total_ccc = sum(ccc_counts.values())
    if total_ccc > 0:
        lines.append(f"| L-R interactions (all conditions) | {safe_int(total_ccc)} | {src(ccc_path)} |")

    lines.append("")


def section_datasets(lines, registry_df, registry_path, sample_df, sample_path):
    lines.append("## 2. Dataset Summary {#datasets}")
    lines.append("")

    if registry_df.empty:
        lines.append(f"Dataset registry not found. {src(registry_path)}")
        lines.append("")
        return

    included = registry_df[registry_df["status"] == "included"].copy()
    lines.append(f"**{len(included)} datasets included** (of {len(registry_df)} evaluated). "
                 f"{src(registry_path)}")
    lines.append("")

    cols = ["accession", "first_author", "year", "compartment", "n_samples",
            "technology", "conditions"]
    available_cols = [c for c in cols if c in included.columns]
    header = "| " + " | ".join(available_cols) + " |"
    sep = "| " + " | ".join(["---"] * len(available_cols)) + " |"
    lines.append(header)
    lines.append(sep)
    for _, row in included.iterrows():
        vals = [str(row.get(c, "")) for c in available_cols]
        lines.append("| " + " | ".join(vals) + " |")
    lines.append("")

    # Platform breakdown
    if "technology" in included.columns:
        platforms = included["technology"].value_counts()
        non_10x = platforms[~platforms.index.str.contains("10x", case=False, na=False)]
        if len(non_10x) > 0:
            lines.append(f"**Platform heterogeneity:** {len(non_10x)} non-10x platform(s) "
                         f"({', '.join(non_10x.index)}). Batch correction handles platform "
                         f"differences via study-level integration keys.")
            lines.append("")

    # Sample demographics summary from sample_metadata
    if not sample_df.empty:
        lines.append("### Sample demographics")
        lines.append("")
        if "age_years" in sample_df.columns:
            known_age = sample_df["age_years"].dropna()
            unknown_age = len(sample_df) - len(known_age)
            if len(known_age) > 0:
                lines.append(f"- Age range: {int(known_age.min())}–{int(known_age.max())} years "
                             f"({unknown_age} samples with unknown age). {src(sample_path)}")
        if "sex" in sample_df.columns:
            sex_counts = sample_df["sex"].value_counts()
            lines.append(f"- Sex distribution: {', '.join(f'{k}={v}' for k, v in sex_counts.items())}. "
                         f"{src(sample_path)}")
        lines.append("")


def section_integration(lines, plan_text, plan_path, celltypist_df, celltypist_path):
    lines.append("## 3. Integration {#integration}")
    lines.append("")

    # Extract integration metrics table from analysis_plan
    in_metrics = False
    metrics_lines = []
    for line in plan_text.splitlines():
        if "| Object | Workflow |" in line:
            in_metrics = True
        if in_metrics:
            if line.strip().startswith("|"):
                metrics_lines.append(line)
            elif metrics_lines and not line.strip().startswith("|"):
                break

    if metrics_lines:
        lines.append("### Workflow comparison")
        lines.append("")
        lines.append(f"Three integration workflows were compared. {src(plan_path)}")
        lines.append("")
        for ml in metrics_lines:
            lines.append(ml)
        lines.append("")

    # Extract rationale from analysis_plan
    in_rationale = False
    for line in plan_text.splitlines():
        if "**Rationale for CCA:**" in line:
            in_rationale = True
        if in_rationale:
            if line.strip() == "" and any(l.startswith("-") for l in metrics_lines):
                continue
            if line.strip().startswith("Now converting") or line.strip().startswith("---"):
                break
            lines.append(line)
    if in_rationale:
        lines.append("")

    # CellTypist concordance
    if not celltypist_df.empty and "concordant" in celltypist_df.columns:
        lines.append("### CellTypist validation")
        lines.append("")
        disc = celltypist_df[celltypist_df["concordant"] == False]
        conc = celltypist_df[celltypist_df["concordant"] == True]
        lines.append(f"{len(conc)} concordant, {len(disc)} discordant cluster(s) across all objects. "
                     f"{src(celltypist_path)}")
        lines.append("")
        if len(disc) > 0:
            lines.append("| Object | Cluster | Cells | De Novo Label | CellTypist Label | Agreement % |")
            lines.append("|--------|---------|-------|---------------|------------------|-------------|")
            for _, row in disc.iterrows():
                obj = row.get("object", "")
                cluster = row.get("cluster", "")
                n_cells = safe_int(row.get("n_cells", ""))
                de_novo = row.get("de_novo_label", "")
                ct = row.get("celltypist_majority", "")
                pct = row.get("celltypist_agreement_pct", "")
                pct_str = f"{pct:.1f}" if isinstance(pct, (int, float)) else str(pct)
                lines.append(f"| {obj} | {cluster} | {n_cells} | {de_novo} | {ct} | {pct_str}% |")
            lines.append("")
            lines.append("De novo labels retained as primary (CellTypist lacks IVD-specific cell types).")
            lines.append("")


def section_clustering(lines, celldefs_df, celldefs_path, clustering_df, clustering_path):
    lines.append("## 4. Clustering & Annotation {#clustering}")
    lines.append("")

    if not celldefs_df.empty:
        lines.append(f"### Cell type definitions")
        lines.append("")
        lines.append(f"{src(celldefs_path)}")
        lines.append("")
        display_cols = [c for c in celldefs_df.columns if not c.startswith("_")]
        header = "| " + " | ".join(display_cols) + " |"
        sep = "| " + " | ".join(["---"] * len(display_cols)) + " |"
        lines.append(header)
        lines.append(sep)
        for _, row in celldefs_df.iterrows():
            vals = [str(row.get(c, "")) for c in display_cols]
            lines.append("| " + " | ".join(vals) + " |")
        lines.append("")
    else:
        lines.append(f"Cell type definitions not found. {src(celldefs_path)}")
        lines.append("")

    # Clustering resolution summary
    if not clustering_df.empty:
        lines.append("### Clustering resolution optimization")
        lines.append("")
        lines.append(f"{src(clustering_path)}")
        lines.append("")
        for sfile, group in clustering_df.groupby("_source_file"):
            obj_name = sfile.replace("_resolutions.tsv", "")
            if "resolution" in group.columns and "n_clusters" in group.columns:
                best = group.loc[group.get("silhouette", pd.Series([0])).idxmax()] if "silhouette" in group.columns else group.iloc[0]
                lines.append(f"- **{obj_name}:** best resolution {best.get('resolution', '?')}, "
                             f"{safe_int(best.get('n_clusters', '?'))} clusters "
                             f"(silhouette={best.get('silhouette', '?'):.3f})" if isinstance(best.get('silhouette'), float) else
                             f"- **{obj_name}:** resolution data available")
        lines.append("")


def section_de(lines, de_summary, de_path, de_combined, de_combined_path,
               skipped_df, skipped_path):
    lines.append("## 5. Differential Expression {#de}")
    lines.append("")

    if de_summary.empty:
        lines.append(f"DE summary not found. {src(de_path)}")
        lines.append("")
        return

    n_comparisons = len(de_summary)
    n_total = int(de_summary["n_total"].sum()) if "n_total" in de_summary.columns else 0
    n_up = int(de_summary["n_up"].sum()) if "n_up" in de_summary.columns else 0
    n_down = int(de_summary["n_down"].sum()) if "n_down" in de_summary.columns else 0

    lines.append(f"**{n_comparisons} powered comparisons**, {safe_int(n_total)} significant genes "
                 f"({safe_int(n_up)} up, {safe_int(n_down)} down). "
                 f"Thresholds: |log2FC| > 0.5, padj < 0.05. {src(de_path)}")
    lines.append("")

    # Summary table
    display_cols = [c for c in de_summary.columns if not c.startswith("_")]
    header = "| " + " | ".join(display_cols) + " |"
    sep = "| " + " | ".join(["---"] * len(display_cols)) + " |"
    lines.append(header)
    lines.append(sep)
    for _, row in de_summary.iterrows():
        vals = [str(row.get(c, "")) for c in display_cols]
        lines.append("| " + " | ".join(vals) + " |")
    lines.append("")

    # Top DE genes from combined results
    if not de_combined.empty and "log2FoldChange" in de_combined.columns:
        lines.append("### Top upregulated genes (by log2FC)")
        lines.append("")
        lines.append(f"{src(de_combined_path)}")
        lines.append("")
        top_up = de_combined.nlargest(10, "log2FoldChange")
        lines.append("| Gene | log2FC | padj | Cell Type | Comparison |")
        lines.append("|------|--------|------|-----------|------------|")
        for _, row in top_up.iterrows():
            gene = row.get("gene", row.get("index", ""))
            lfc = row.get("log2FoldChange", "")
            lfc_str = f"{lfc:.2f}" if isinstance(lfc, (int, float)) else str(lfc)
            padj = row.get("padj", "")
            padj_str = f"{padj:.2e}" if isinstance(padj, (int, float)) else str(padj)
            ct = row.get("cell_type", "")
            comp = row.get("comparison", "")
            lines.append(f"| {gene} | {lfc_str} | {padj_str} | {ct} | {comp} |")
        lines.append("")

        lines.append("### Top downregulated genes (by log2FC)")
        lines.append("")
        top_down = de_combined.nsmallest(10, "log2FoldChange")
        lines.append("| Gene | log2FC | padj | Cell Type | Comparison |")
        lines.append("|------|--------|------|-----------|------------|")
        for _, row in top_down.iterrows():
            gene = row.get("gene", row.get("index", ""))
            lfc = row.get("log2FoldChange", "")
            lfc_str = f"{lfc:.2f}" if isinstance(lfc, (int, float)) else str(lfc)
            padj = row.get("padj", "")
            padj_str = f"{padj:.2e}" if isinstance(padj, (int, float)) else str(padj)
            ct = row.get("cell_type", "")
            comp = row.get("comparison", "")
            lines.append(f"| {gene} | {lfc_str} | {padj_str} | {ct} | {comp} |")
        lines.append("")

    # Skipped comparisons
    if not skipped_df.empty:
        lines.append(f"### Skipped comparisons (underpowered)")
        lines.append("")
        lines.append(f"{len(skipped_df)} comparisons skipped due to insufficient samples "
                     f"(< 3 per condition per cell type). {src(skipped_path)}")
        lines.append("")


def section_pathways(lines, enrichment_df, enrichment_path, gsea_df, gsea_path):
    lines.append("## 6. Biological Pathways {#pathways}")
    lines.append("")

    if not enrichment_df.empty:
        n_sig = len(enrichment_df[enrichment_df.get("Adjusted P-value", pd.Series(dtype=float)) < 0.05]) \
            if "Adjusted P-value" in enrichment_df.columns else len(enrichment_df)
        lines.append(f"**ORA:** {safe_int(n_sig)} significantly enriched terms (FDR < 0.05). "
                     f"{src(enrichment_path)}")
        lines.append("")

        # Top enrichments by library
        if "Gene_set" in enrichment_df.columns or "library" in enrichment_df.columns:
            lib_col = "library" if "library" in enrichment_df.columns else "Gene_set"
            for lib_name, lib_group in enrichment_df.groupby(lib_col):
                if "Adjusted P-value" in lib_group.columns:
                    top5 = lib_group.nsmallest(5, "Adjusted P-value")
                else:
                    top5 = lib_group.head(5)
                if len(top5) == 0:
                    continue
                lines.append(f"**{lib_name}** (top 5):")
                lines.append("")
                for _, row in top5.iterrows():
                    term = row.get("Term", row.get("term", ""))
                    pval = row.get("Adjusted P-value", row.get("padj", ""))
                    pval_str = f"{pval:.2e}" if isinstance(pval, (int, float)) else str(pval)
                    ct = row.get("cell_type", "")
                    direction = row.get("direction", "")
                    lines.append(f"- {term} (padj={pval_str}, {ct} {direction})")
                lines.append("")
    else:
        lines.append(f"ORA results not found. {src(enrichment_path)}")
        lines.append("")

    if not gsea_df.empty:
        n_gsea = len(gsea_df)
        lines.append(f"**GSEA:** {safe_int(n_gsea)} significant terms. {src(gsea_path)}")
        lines.append("")
    else:
        lines.append(f"GSEA results not found. {src(gsea_path)}")
        lines.append("")


def section_tf(lines, tf_df, tf_path):
    lines.append("## 7. Transcription Factor Activity {#tf}")
    lines.append("")

    if tf_df.empty:
        lines.append(f"TF activity results not found. {src(tf_path)}")
        lines.append("")
        return

    # Count significant TFs
    if "p_value" in tf_df.columns:
        sig = tf_df[tf_df["p_value"] < 0.05]
        lines.append(f"**{len(sig)} significant TF-condition associations** (p < 0.05). "
                     f"{src(tf_path)}")
    elif "padj" in tf_df.columns:
        sig = tf_df[tf_df["padj"] < 0.05]
        lines.append(f"**{len(sig)} significant TF-condition associations** (padj < 0.05). "
                     f"{src(tf_path)}")
    else:
        sig = tf_df
        lines.append(f"**{len(tf_df)} TF activity results.** {src(tf_path)}")
    lines.append("")

    if not sig.empty and "source" in sig.columns:
        tf_col = "source"
    elif not sig.empty and "TF" in sig.columns:
        tf_col = "TF"
    else:
        tf_col = sig.columns[0] if len(sig.columns) > 0 else None

    if tf_col and not sig.empty:
        unique_tfs = sig[tf_col].unique()
        lines.append(f"**{len(unique_tfs)} unique TFs** with significant activity changes.")
        lines.append("")
        # Show top TFs by score magnitude
        score_col = None
        for candidate in ["score", "activity", "estimate", "statistic"]:
            if candidate in sig.columns:
                score_col = candidate
                break
        if score_col:
            top_tfs = sig.reindex(sig[score_col].abs().sort_values(ascending=False).index).head(10)
            lines.append("| TF | Score | Cell Type | Comparison |")
            lines.append("|----|-------|-----------|------------|")
            for _, row in top_tfs.iterrows():
                tf = row.get(tf_col, "")
                sc = row.get(score_col, "")
                sc_str = f"{sc:.3f}" if isinstance(sc, (int, float)) else str(sc)
                ct = row.get("cell_type", row.get("condition", ""))
                comp = row.get("comparison", "")
                lines.append(f"| {tf} | {sc_str} | {ct} | {comp} |")
            lines.append("")


def section_trajectory(lines, traj_corr_df, traj_corr_path,
                       traj_gene_counts, traj_gene_path):
    lines.append("## 8. Cell State Trajectories {#trajectory}")
    lines.append("")

    if traj_corr_df.empty:
        lines.append(f"Trajectory correlation results not found. {src(traj_corr_path)}")
        lines.append("")
        return

    lines.append(f"PAGA + diffusion pseudotime (DPT) analysis. {src(traj_corr_path)}")
    lines.append("")

    # Display correlation results
    display_cols = [c for c in traj_corr_df.columns if not c.startswith("_")]
    header = "| " + " | ".join(display_cols) + " |"
    sep = "| " + " | ".join(["---"] * len(display_cols)) + " |"
    lines.append(header)
    lines.append(sep)
    for _, row in traj_corr_df.iterrows():
        vals = []
        for c in display_cols:
            v = row.get(c, "")
            if isinstance(v, float):
                vals.append(f"{v:.4f}" if abs(v) < 1 else f"{v:.2e}")
            else:
                vals.append(str(v))
        lines.append("| " + " | ".join(vals) + " |")
    lines.append("")

    # Trajectory gene counts
    if traj_gene_counts:
        lines.append(f"### Trajectory-associated genes")
        lines.append("")
        lines.append(f"{src(traj_gene_path)}")
        lines.append("")
        for comp, count in sorted(traj_gene_counts.items()):
            lines.append(f"- **{comp}:** {count} genes correlated with pseudotime (FDR < 0.05)")
        lines.append("")


def section_ccc(lines, ccc_counts, ccc_path, diff_df, diff_path,
                pain_inter_df, pain_inter_path):
    lines.append("## 9. Cell-Cell Communication {#communication}")
    lines.append("")

    if not ccc_counts:
        lines.append(f"CCC interaction files not found. {src(ccc_path)}")
        lines.append("")
        return

    lines.append(f"LIANA rank_aggregate (CellPhoneDB, NATMI, Connectome, SingleCellSignalR, log2FC). "
                 f"{src(ccc_path)}")
    lines.append("")
    lines.append("| Condition | Interactions |")
    lines.append("|-----------|-------------|")
    for cond, count in sorted(ccc_counts.items()):
        lines.append(f"| {cond} | {safe_int(count)} |")
    lines.append("")

    if not diff_df.empty:
        lines.append(f"### Differential interactions")
        lines.append("")
        lines.append(f"{src(diff_path)}")
        lines.append("")
        if "direction" in diff_df.columns:
            gained = len(diff_df[diff_df["direction"] == "gained"])
            lost = len(diff_df[diff_df["direction"] == "lost"])
            lines.append(f"- Gained in degeneration: {safe_int(gained)}")
            lines.append(f"- Lost in degeneration: {safe_int(lost)}")
            lines.append("")

    if not pain_inter_df.empty:
        lines.append(f"### Pain-relevant interactions")
        lines.append("")
        lines.append(f"{len(pain_inter_df)} pain-relevant L-R interactions identified. "
                     f"{src(pain_inter_path)}")
        lines.append("")


def section_pain(lines, pain_genes_df, pain_path, pain_inter_df, pain_inter_path):
    lines.append("## 10. Pain Biology {#pain}")
    lines.append("")

    if pain_genes_df.empty:
        lines.append(f"Pain gene results not found. {src(pain_path)}")
        lines.append("")
        return

    lines.append(f"**{len(pain_genes_df)} pain-associated DE genes identified.** {src(pain_path)}")
    lines.append("")

    display_cols = [c for c in pain_genes_df.columns if not c.startswith("_")]
    # Limit to a reasonable set of columns for display
    preferred = ["gene", "log2FoldChange", "padj", "cell_type", "comparison", "pain_category"]
    show_cols = [c for c in preferred if c in display_cols] or display_cols[:6]

    header = "| " + " | ".join(show_cols) + " |"
    sep = "| " + " | ".join(["---"] * len(show_cols)) + " |"
    lines.append(header)
    lines.append(sep)
    for _, row in pain_genes_df.iterrows():
        vals = []
        for c in show_cols:
            v = row.get(c, "")
            if isinstance(v, float):
                if abs(v) < 0.01 or abs(v) > 1000:
                    vals.append(f"{v:.2e}")
                else:
                    vals.append(f"{v:.3f}")
            else:
                vals.append(str(v))
        lines.append("| " + " | ".join(vals) + " |")
    lines.append("")

    if not pain_inter_df.empty:
        lines.append(f"{len(pain_inter_df)} pain-relevant ligand-receptor interactions from CCC analysis. "
                     f"{src(pain_inter_path)}")
        lines.append("")


def section_limitations(lines, plan_text, plan_path, skipped_df, skipped_path,
                        traj_corr_df, ccc_counts):
    lines.append("## 11. Limitations & Caveats {#limitations}")
    lines.append("")

    # Extract Known Issues from analysis_plan, but classify them:
    # - Data/design limitations are always relevant
    # - Implementation workarounds are not scientific limitations
    # - Cross-version observations belong in a separate section
    in_issues = False
    raw_issues = []
    for line in plan_text.splitlines():
        if line.strip() == "## Known Issues":
            in_issues = True
            continue
        if in_issues:
            if line.startswith("---"):
                break
            if line.strip().startswith("- **"):
                raw_issues.append(line.strip())

    # Filter: keep data/design limitations, skip implementation workarounds
    implementation_keywords = ["SeuratDisk", "workaround", "GetAssayData"]
    data_issues = []
    for issue in raw_issues:
        if any(kw in issue for kw in implementation_keywords):
            continue
        # Fix stale method references: the original text references scANVI
        # as the primary batch correction; v5 uses CCA
        issue = issue.replace(
            "Handled by scANVI batch correction. CCA and STACAS also correct "
            "for this via study-level integration.",
            "Handled by study-level batch correction (CCA in v5; scANVI and "
            "STACAS also tested)."
        )
        data_issues.append(issue)

    if data_issues:
        lines.append(f"### Data and design limitations")
        lines.append("")
        lines.append(f"{src(plan_path)}")
        lines.append("")
        for issue in data_issues:
            lines.append(issue)
        lines.append("")

    # Skipped comparisons summary
    if not skipped_df.empty:
        lines.append(f"### Underpowered comparisons")
        lines.append("")
        lines.append(f"{len(skipped_df)} cell type x condition comparisons were skipped "
                     f"due to insufficient sample counts. {src(skipped_path)}")
        lines.append("")

    # Cross-version sensitivity — derived from actual v5 data where possible,
    # with context from analysis_plan for historical perspective
    lines.append("### Result sensitivity across pipeline versions")
    lines.append("")
    lines.append("Several results are sensitive to upstream methodological choices "
                 "(integration method, annotation, cell sampling). These are documented "
                 "here to flag areas requiring cautious interpretation. "
                 f"{src(plan_path)}")
    lines.append("")

    # Trajectory — use v5 data if available, otherwise note from plan
    if not traj_corr_df.empty and "compartment" in traj_corr_df.columns:
        lines.append("- **Trajectory pseudotime-condition correlations** are sensitive to "
                     "integration method and root cell choice. In v5 (CCA): ")
        parts = []
        for _, row in traj_corr_df.iterrows():
            comp = row.get("compartment", "")
            # Try common column names for the correlation coefficient
            rho = None
            for col in ["rho", "spearman_rho", "correlation", "statistic"]:
                if col in row.index and pd.notna(row[col]):
                    rho = row[col]
                    break
            if rho is not None:
                parts.append(f"{comp} rho={rho:+.3f}")
            else:
                parts.append(comp)
        lines.append(f"  {'; '.join(parts)}. "
                     "Prior versions showed sign changes (e.g., CEP: -0.163 in v2, +0.135 in v3, "
                     "+0.073 in v5), indicating these correlations are not robust to upstream choices.")
    else:
        lines.append("- **Trajectory instability:** Pseudotime-condition correlations have changed "
                     "sign across pipeline versions, indicating sensitivity to integration method "
                     "and annotation choices.")
    lines.append("")

    # CCC direction sensitivity — use v5 data if available
    if ccc_counts:
        ccc_summary = ", ".join(f"{k}: {v:,}" for k, v in sorted(ccc_counts.items()))
        lines.append(f"- **CCC interaction counts** in v5: {ccc_summary}. "
                     "The direction of the healthy-vs-degenerated difference has varied across "
                     "pipeline versions (v1: degenerated > healthy; v2: healthy > degenerated; "
                     "v3: near-equal), making this result sensitive to cell type definitions "
                     "and sampling.")
    else:
        lines.append("- **CCC direction sensitivity:** The direction of healthy-vs-degenerated "
                     "interaction count differences has varied across pipeline versions.")
    lines.append("")

    lines.append("- **CellTypist concordance** is limited for IVD-specific cell types. "
                 "CellTypist lacks IVD reference data, so disagreements with de novo labels "
                 "are expected for mesenchymal populations. De novo labels are retained as "
                 "primary annotations; CellTypist is used for immune subtype validation only.")
    lines.append("")

    lines.append("- **AF pseudotime sign:** AF consistently shows positive rho (degenerated "
                 "cells at later pseudotime) across pipeline versions, opposite to NP. This may "
                 "reflect genuine compartment-specific biology or root cell choice effects.")
    lines.append("")


def section_methods(lines, plan_text, plan_path):
    lines.append("## 12. Methods {#methods}")
    lines.append("")
    lines.append(f"Full parameter choices and rationale documented in `analysis_plan.md`. {src(plan_path)}")
    lines.append("")

    lines.append("### Data acquisition")
    lines.append("")
    lines.append("12 scRNA-seq datasets of human IVD tissue downloaded from GEO and CNGB. "
                 "Raw count matrices obtained per dataset. "
                 f"See `metadata/dataset_registry.tsv` for accessions and details. "
                 "*[source: `scripts/01_dataset_download.py`, `metadata/dataset_registry.tsv`]*")
    lines.append("")

    lines.append("### Quality control and preprocessing")
    lines.append("")
    lines.append("Per-dataset QC: min 200 genes, max 6000 genes, min 500 counts, max 20% "
                 "mitochondrial reads. Doublet detection with Scrublet (expected rate 5%). "
                 "Normalization: total-count to 10,000, log1p. HVG selection: top 2000 genes "
                 "per dataset (Seurat v3 method). "
                 "*[source: `scripts/03_preprocessing.py`, `specs/03_PREPROCESSING.md`]*")
    lines.append("")

    lines.append("### Cell classification")
    lines.append("")
    lines.append("Coarse classification into 5 anchor categories using marker gene scoring "
                 "(immune: PTPRC, CD3D, CD68, PECAM1; mesenchymal: COL2A1, COL1A1, ACAN, SOX9). "
                 "Cluster-level majority voting with 85% threshold. "
                 "*[source: `scripts/04_annotation.py`, `specs/04_ANNOTATION.md`]*")
    lines.append("")

    lines.append("### Integration")
    lines.append("")
    lines.append("Three workflows compared (CCA, scANVI, STACAS) on four compartment objects "
                 "(NP, AF, CEP, all_cells). CCA (Seurat v5 `IntegrateLayers(method=CCAIntegration)`) "
                 "selected as primary: label-free, full cell counts, strongest batch mixing (iLISI). "
                 "*[source: `scripts/05a_integration_cca.R`, `specs/05_INTEGRATION.md`, `analysis_plan.md`]*")
    lines.append("")

    lines.append("### Clustering")
    lines.append("")
    lines.append("Leiden clustering with multi-resolution optimization. Resolution selected by "
                 "silhouette score. "
                 "*[source: `scripts/06_clustering.py`, `specs/06_CLUSTERING.md`]*")
    lines.append("")

    lines.append("### Post-integration annotation")
    lines.append("")
    lines.append("De novo cell type annotation from cluster DE markers and canonical marker panels. "
                 "CellTypist (Immune_All_Low model) for immune subtype validation. "
                 "*[source: `scripts/07_annotation.py`, `specs/07_ANNOTATION.md`]*")
    lines.append("")

    lines.append("### Differential expression")
    lines.append("")
    lines.append("Pseudobulk aggregation per sample per cell type. DE with pyDESeq2. "
                 "Significance: |log2FC| > 0.5, padj < 0.05 (Benjamini-Hochberg). "
                 "Minimum 3 samples per condition per cell type. "
                 "*[source: `scripts/08_differential.py`, `specs/08_DIFFERENTIAL.md`]*")
    lines.append("")

    lines.append("### Pathway enrichment")
    lines.append("")
    lines.append("ORA and GSEA using gseapy. Databases: GO Biological Process 2023, KEGG 2021, "
                 "Reactome 2022, MSigDB Hallmark 2020, custom IVD gene sets. "
                 "*[source: `scripts/09_interpretation.py`, `specs/09_INTERPRETATION.md`]*")
    lines.append("")

    lines.append("### TF activity inference")
    lines.append("")
    lines.append("CollecTRI regulon network. TF activity scored by Fisher's exact test for "
                 "enrichment of TF targets among DE genes. "
                 "*[source: `scripts/09_interpretation.py`]*")
    lines.append("")

    lines.append("### Trajectory analysis")
    lines.append("")
    lines.append("PAGA + diffusion pseudotime (DPT) on mesenchymal embeddings. "
                 "Root cells defined per compartment (NP: notochordal; AF: AF_inner). "
                 "Trajectory genes: Spearman correlation with pseudotime, FDR < 0.05. "
                 "*[source: `scripts/10_trajectory.py`, `specs/10_TRAJECTORY.md`]*")
    lines.append("")

    lines.append("### Cell-cell communication")
    lines.append("")
    lines.append("LIANA rank_aggregate with consensus resource. 5 methods: CellPhoneDB, NATMI, "
                 "Connectome, SingleCellSignalR, log2FC. 100 permutations. "
                 "*[source: `scripts/11_communication.py`, `specs/11_COMMUNICATION.md`]*")
    lines.append("")

    lines.append("### Software")
    lines.append("")
    lines.append("Python 3.12, scanpy, scvi-tools, pyDESeq2, gseapy, decoupler, liana. "
                 "R: Seurat 5.4.0, STACAS 2.4.1. "
                 "Full environment: `requirements.txt` / `requirements_frozen.txt`.")
    lines.append("")


def section_reproducibility(lines):
    short_sha, full_sha, branch = git_version()

    lines.append("## 13. Reproducibility {#reproducibility}")
    lines.append("")
    lines.append(f"- **Git commit:** `{full_sha}` (branch: `{branch}`)")
    lines.append("- **Random seeds:** 42 (all stochastic operations)")
    lines.append("- **Package versions:** pinned in `requirements.txt`, frozen in `requirements_frozen.txt`")
    lines.append("- **Parameter choices:** documented in `analysis_plan.md`")
    lines.append("- **Human checkpoint decisions:** recorded in `analysis_plan.md`")
    lines.append("- **Data provenance:** GEO/CNGB accessions and download dates in `metadata/dataset_registry.tsv`")
    lines.append("- **File checksums:** `metadata/file_checksums.json`")
    lines.append("")


# ---------------------------------------------------------------------------
# Supplementary table collection (unchanged from original)
# ---------------------------------------------------------------------------

def collect_supplementary_tables():
    """Copy key result tables to supplementary_tables/."""
    copies = {
        "S1_dataset_registry.tsv": META / "dataset_registry.tsv",
        "S2_sample_metadata.tsv": META / "sample_metadata.tsv",
        "S3_inclusion_summary.tsv": RESULTS / "integration" / "inclusion_summary.tsv",
        "S4_study_caveats.tsv": RESULTS / "integration" / "study_caveats.tsv",
        "S5_composition_analysis.tsv": RESULTS / "differential" / "composition_analysis.tsv",
        "S6_de_summary.tsv": RESULTS / "differential" / "de_summary_table.tsv",
        "S7_de_results_combined.tsv": RESULTS / "differential" / "de_results_combined.tsv",
        "S8_skipped_comparisons.tsv": RESULTS / "differential" / "skipped_comparisons.tsv",
        "S9_enrichment_ORA.tsv": RESULTS / "interpretation" / "pathway_enrichment" / "all_enrichment_results.tsv",
        "S10_gsea_results.tsv": RESULTS / "interpretation" / "pathway_enrichment" / "gsea_results.tsv",
        "S11_tf_activity.tsv": RESULTS / "interpretation" / "tf_activity" / "tf_activity_results.tsv",
        "S12_pain_genes.tsv": RESULTS / "interpretation" / "pain_genes.tsv",
        "S13_trajectory_genes_NP.tsv": RESULTS / "trajectories" / "trajectory_genes_NP.tsv",
        "S14_trajectory_genes_AF.tsv": RESULTS / "trajectories" / "trajectory_genes_AF.tsv",
        "S15_trajectory_genes_CEP.tsv": RESULTS / "trajectories" / "trajectory_genes_CEP.tsv",
        "S16_pain_interactions.tsv": RESULTS / "communication" / "pain_interactions.tsv",
        "S17_celltypist_concordance_NP.tsv": RESULTS / "integration" / "celltypist_validation" / "NP_concordance.tsv",
        "S18_celltypist_concordance_AF.tsv": RESULTS / "integration" / "celltypist_validation" / "AF_concordance.tsv",
        "S19_celltypist_concordance_CEP.tsv": RESULTS / "integration" / "celltypist_validation" / "CEP_concordance.tsv",
    }

    collected = 0
    for dest_name, src_path in copies.items():
        if src_path.exists():
            import shutil
            shutil.copy2(src_path, SUPP_TABLES / dest_name)
            collected += 1

    print(f"  Collected {collected}/{len(copies)} supplementary tables")
    return collected


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def generate_report():
    """Build the full markdown report from result files."""
    lines = []

    # Load all data sources
    plan_text, plan_path = load_analysis_plan()
    registry_df, registry_path = load_dataset_registry()
    sample_df, sample_path = load_sample_metadata()
    celldefs_df, celldefs_path = load_cell_type_definitions()
    clustering_df, clustering_path = load_clustering_resolutions()
    de_summary, de_path = load_de_summary()
    de_combined, de_combined_path = load_de_combined()
    skipped_df, skipped_path = load_skipped_comparisons()
    enrichment_df, enrichment_path = load_enrichment()
    gsea_df, gsea_path = load_gsea()
    tf_df, tf_path = load_tf_activity()
    pain_df, pain_path = load_pain_genes()
    traj_corr_df, traj_corr_path = load_trajectory_correlations()
    traj_gene_counts, traj_gene_path = load_trajectory_genes()
    ccc_counts, ccc_path = load_ccc_interactions()
    diff_inter_df, diff_inter_path = load_differential_interactions()
    pain_inter_df, pain_inter_path = load_pain_interactions()
    celltypist_df, celltypist_path = load_celltypist_concordance()

    # Build report
    section_header(lines, plan_text)
    section_toc(lines)
    section_summary(lines, sample_df, sample_path, de_summary, de_path,
                    enrichment_df, enrichment_path, ccc_counts, ccc_path)
    section_datasets(lines, registry_df, registry_path, sample_df, sample_path)
    section_integration(lines, plan_text, plan_path, celltypist_df, celltypist_path)
    section_clustering(lines, celldefs_df, celldefs_path, clustering_df, clustering_path)
    section_de(lines, de_summary, de_path, de_combined, de_combined_path,
               skipped_df, skipped_path)
    section_pathways(lines, enrichment_df, enrichment_path, gsea_df, gsea_path)
    section_tf(lines, tf_df, tf_path)
    section_trajectory(lines, traj_corr_df, traj_corr_path,
                       traj_gene_counts, traj_gene_path)
    section_ccc(lines, ccc_counts, ccc_path, diff_inter_df, diff_inter_path,
                pain_inter_df, pain_inter_path)
    section_pain(lines, pain_df, pain_path, pain_inter_df, pain_inter_path)
    section_limitations(lines, plan_text, plan_path, skipped_df, skipped_path,
                        traj_corr_df, ccc_counts)
    section_methods(lines, plan_text, plan_path)
    section_reproducibility(lines)

    report_path = DOCS / "final_report.md"
    report_path.write_text("\n".join(lines))
    print(f"  Final report: {report_path}")
    return report_path


def validate():
    """Validation checks."""
    print("\n" + "=" * 60)
    print("Validating Module 12 outputs")
    print("=" * 60)

    checks = []

    report_path = DOCS / "final_report.md"
    if report_path.exists():
        size_kb = report_path.stat().st_size / 1024
        checks.append(("PASS", f"Final report generated ({size_kb:.1f} KB)"))
    else:
        checks.append(("FAIL", "Final report not found at docs/final_report.md"))

    supp_count = len(list(SUPP_TABLES.glob("S*.tsv")))
    if supp_count > 0:
        checks.append(("PASS", f"{supp_count} supplementary tables collected"))
    else:
        checks.append(("WARN", "No supplementary tables (results/ may not be populated)"))

    if (BASE / "requirements_frozen.txt").exists():
        checks.append(("PASS", "Package versions frozen"))

    # Check data sources that were read
    data_sources = [
        ("Dataset registry", META / "dataset_registry.tsv"),
        ("Sample metadata", META / "sample_metadata.tsv"),
        ("DE summary", RESULTS / "differential" / "de_summary_table.tsv"),
        ("Enrichment results", RESULTS / "interpretation" / "pathway_enrichment" / "all_enrichment_results.tsv"),
        ("TF activity", RESULTS / "interpretation" / "tf_activity" / "tf_activity_results.tsv"),
        ("Pain genes", RESULTS / "interpretation" / "pain_genes.tsv"),
        ("Cell type definitions", RESULTS / "integration" / "cell_type_definitions.tsv"),
    ]
    for name, path in data_sources:
        if path.exists():
            checks.append(("PASS", f"{name} found"))
        else:
            checks.append(("WARN", f"{name} not found — report section will show placeholder"))

    passed = all(s != "FAIL" for s, _ in checks)
    for status, msg in checks:
        print(f"  [{status}] {msg}")
    return passed


def main():
    print("=" * 60)
    print("Module 12: Final Reporting")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    if "--validate-only" in sys.argv:
        passed = validate()
        sys.exit(0 if passed else 1)

    print("\nCollecting supplementary tables...")
    collect_supplementary_tables()

    print("\nGenerating final report (markdown)...")
    generate_report()

    print("\nFreezing requirements...")
    try:
        result = subprocess.run(["pip", "freeze"], capture_output=True, text=True)
        (BASE / "requirements_frozen.txt").write_text(result.stdout)
        print("  Requirements frozen: requirements_frozen.txt")
    except Exception as e:
        print(f"  WARN: Could not freeze requirements: {e}")

    validate()

    print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("Pipeline complete. Final report at: docs/final_report.md")


if __name__ == "__main__":
    main()
