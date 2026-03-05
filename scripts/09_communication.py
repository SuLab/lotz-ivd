#!/usr/bin/env python3
"""Module 09: Cell-cell communication analysis for IVD scRNA-seq atlas.

Uses LIANA to infer ligand-receptor interactions between IVD cell populations,
comparing healthy vs degenerated conditions.

Usage:
    python3 scripts/09_communication.py
    python3 scripts/09_communication.py --validate-only
"""

import gc
import sys
import os
import warnings
import traceback
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent.parent
INT_DIR = BASE / "data" / "integrated"
PROC_DIR = BASE / "data" / "processed"
META_PATH = BASE / "metadata" / "sample_metadata.tsv"
RESULTS_DIR = BASE / "results" / "communication"
PLOT_DIR = RESULTS_DIR / "interaction_plots"

for d in [RESULTS_DIR, PLOT_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# ── Configuration ────────────────────────────────────────────────────────────
# Minimum cells per cell type to include in analysis
MIN_CELLS_PER_TYPE = 50

# Conditions to compare
CONDITIONS = {
    'healthy': ['healthy'],
    'degenerated': ['degenerated_mild', 'degenerated_severe', 'degenerated_ungraded'],
}

# Max cells per condition for LIANA (memory constraint)
MAX_CELLS_PER_CONDITION = 20000

# Pain-relevant ligand-receptor pairs to highlight
PAIN_LIGANDS = ['NGF', 'BDNF', 'NTF3', 'NTF4', 'VEGFA', 'VEGFB', 'VEGFC',
                'SEMA3A', 'SEMA3F', 'NTN1', 'IL1B', 'IL6', 'TNF', 'CCL2',
                'CXCL8', 'CXCL1', 'CXCL2', 'CXCL3', 'FGF2', 'PDGFA', 'PDGFB',
                'PTGS2', 'TAC1', 'CALCA', 'ANGPT1', 'ANGPT2']

PAIN_RECEPTORS = ['NTRK1', 'NTRK2', 'NTRK3', 'NGFR', 'KDR', 'FLT1', 'NRP1',
                  'NRP2', 'IL1R1', 'IL6R', 'TNFRSF1A', 'TNFRSF1B', 'CCR2',
                  'CXCR1', 'CXCR2', 'ROBO1', 'ROBO2', 'TEK']


# ═══════════════════════════════════════════════════════════════════════════════
# Data Loading
# ═══════════════════════════════════════════════════════════════════════════════

def load_data_for_communication():
    """Load and combine resident + non-resident cells per study for CCC analysis.

    For CCC we need multiple cell types in the same dataset. We'll use per-dataset
    processed files which contain all cell types.
    """
    meta = pd.read_csv(META_PATH, sep='\t')

    # Find per-dataset processed files
    processed_files = sorted(PROC_DIR.glob("*.h5ad"))
    print(f"  Found {len(processed_files)} processed datasets")

    # Load and label by condition
    healthy_adatas = []
    degen_adatas = []

    for f in processed_files:
        acc = f.stem
        # Get condition for this dataset's samples
        dataset_meta = meta[meta['study_accession'] == acc]
        if dataset_meta.empty:
            continue

        conditions = dataset_meta['condition_harmonized'].unique()
        has_healthy = any(c in CONDITIONS['healthy'] for c in conditions)
        has_degen = any(c in CONDITIONS['degenerated'] for c in conditions)

        if not (has_healthy or has_degen):
            continue

        print(f"  Loading {acc}...")
        try:
            adata = sc.read_h5ad(f)
        except Exception as e:
            print(f"    ERROR loading {acc}: {e}")
            continue

        # Must have cell_type_final annotation
        if 'cell_type_final' not in adata.obs.columns:
            print(f"    Skipping {acc}: no cell_type_final")
            continue

        # Must have sample_id to link to conditions
        if 'sample_id' not in adata.obs.columns:
            print(f"    Skipping {acc}: no sample_id")
            continue

        # Tag with condition from metadata
        sample_cond = dataset_meta.set_index('sample_id')['condition_harmonized'].to_dict()
        adata.obs['condition_harmonized'] = adata.obs['sample_id'].map(sample_cond)
        adata.obs['study_id'] = acc

        # Split into healthy and degenerated
        for group_name, group_conds in CONDITIONS.items():
            mask = adata.obs['condition_harmonized'].isin(group_conds)
            if mask.sum() < MIN_CELLS_PER_TYPE:
                continue

            sub = adata[mask].copy()
            ct_counts = sub.obs['cell_type_final'].value_counts()
            valid_cts = ct_counts[ct_counts >= MIN_CELLS_PER_TYPE].index
            if len(valid_cts) < 2:
                print(f"    {acc} {group_name}: only {len(valid_cts)} cell type(s) with >={MIN_CELLS_PER_TYPE} cells, skipping")
                continue

            sub = sub[sub.obs['cell_type_final'].isin(valid_cts)].copy()
            print(f"    {acc} {group_name}: {sub.n_obs:,} cells, {len(valid_cts)} cell types")

            if group_name == 'healthy':
                healthy_adatas.append(sub)
            else:
                degen_adatas.append(sub)

        del adata
        gc.collect()

    # Concatenate
    result = {}
    for name, adatas in [('healthy', healthy_adatas), ('degenerated', degen_adatas)]:
        if not adatas:
            print(f"  WARNING: No {name} datasets suitable for CCC")
            continue

        combined = ad.concat(adatas, join='inner')
        print(f"  {name}: {combined.n_obs:,} cells, {combined.n_vars} genes, "
              f"{combined.obs['cell_type_final'].nunique()} cell types")

        # Downsample if too large
        if combined.n_obs > MAX_CELLS_PER_CONDITION:
            print(f"    Downsampling to {MAX_CELLS_PER_CONDITION:,}")
            sc.pp.subsample(combined, n_obs=MAX_CELLS_PER_CONDITION, random_state=42)

        # Ensure normalized for LIANA
        sc.pp.normalize_total(combined, target_sum=1e4)
        sc.pp.log1p(combined)

        result[name] = combined

    return result


# ═══════════════════════════════════════════════════════════════════════════════
# LIANA Analysis
# ═══════════════════════════════════════════════════════════════════════════════

def run_liana(adata, condition_name):
    """Run LIANA cell-cell communication analysis."""
    import liana as li

    print(f"\n  Running LIANA for {condition_name} ({adata.n_obs:,} cells)...")

    # Set cell type column
    adata.obs['cell_type'] = adata.obs['cell_type_final'].astype(str)

    try:
        li.mt.rank_aggregate(
            adata,
            groupby='cell_type',
            use_raw=False,
            verbose=True,
            resource_name='consensus',
            n_perms=100,
        )
    except Exception as e:
        print(f"    LIANA failed: {e}")
        traceback.print_exc()
        return pd.DataFrame()

    # Extract results
    if hasattr(adata, 'uns') and 'liana_res' in adata.uns:
        results = adata.uns['liana_res'].copy()
        results['condition'] = condition_name
        print(f"    LIANA results: {len(results):,} interactions")

        # Save
        results.to_csv(RESULTS_DIR / f"interactions_{condition_name}.tsv",
                       sep='\t', index=False)
        return results
    else:
        print("    WARNING: No LIANA results found")
        return pd.DataFrame()


# ═══════════════════════════════════════════════════════════════════════════════
# Differential Interactions
# ═══════════════════════════════════════════════════════════════════════════════

def compare_interactions(healthy_res, degen_res):
    """Compare interactions between healthy and degenerated."""
    if healthy_res.empty or degen_res.empty:
        return pd.DataFrame()

    print("\n  Comparing healthy vs degenerated interactions...")

    # Create interaction keys
    for df in [healthy_res, degen_res]:
        df['interaction_key'] = (df['source'].astype(str) + '|' +
                                  df['target'].astype(str) + '|' +
                                  df['ligand_complex'].astype(str) + '|' +
                                  df['receptor_complex'].astype(str))

    # Merge
    merged = healthy_res.merge(degen_res, on='interaction_key', how='outer',
                                suffixes=('_healthy', '_degen'))

    # Score columns vary by LIANA version; use magnitude_rank or similar
    score_col_h = None
    score_col_d = None
    for col in merged.columns:
        if 'magnitude_rank' in col and '_healthy' in col:
            score_col_h = col
        if 'magnitude_rank' in col and '_degen' in col:
            score_col_d = col

    if score_col_h and score_col_d:
        merged['rank_diff'] = merged[score_col_h].fillna(1) - merged[score_col_d].fillna(1)
        # Positive rank_diff = gained in degeneration (lower rank = more significant)
        merged = merged.sort_values('rank_diff')

    merged.to_csv(RESULTS_DIR / "differential_interactions.tsv", sep='\t', index=False)
    print(f"  Differential interactions: {len(merged):,}")

    return merged


# ═══════════════════════════════════════════════════════════════════════════════
# Pain-Relevant Interaction Flagging
# ═══════════════════════════════════════════════════════════════════════════════

def flag_pain_interactions(results):
    """Flag pain-relevant ligand-receptor interactions."""
    if results.empty:
        return results

    pain_flags = []
    for _, row in results.iterrows():
        ligand = str(row.get('ligand_complex', ''))
        receptor = str(row.get('receptor_complex', ''))

        is_pain = False
        categories = []

        for pl in PAIN_LIGANDS:
            if pl in ligand:
                is_pain = True
                categories.append(f'pain_ligand:{pl}')
        for pr in PAIN_RECEPTORS:
            if pr in receptor:
                is_pain = True
                categories.append(f'pain_receptor:{pr}')

        pain_flags.append({
            'is_pain_relevant': is_pain,
            'pain_categories': '; '.join(categories) if categories else '',
        })

    pain_df = pd.DataFrame(pain_flags)
    results = pd.concat([results.reset_index(drop=True), pain_df], axis=1)

    pain_results = results[results['is_pain_relevant']]
    if not pain_results.empty:
        pain_results.to_csv(RESULTS_DIR / "pain_interactions.tsv", sep='\t', index=False)
        print(f"  Pain-relevant interactions: {len(pain_results):,}")

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# Visualization
# ═══════════════════════════════════════════════════════════════════════════════

def plot_interaction_summary(results, condition_name):
    """Plot interaction summary for a condition."""
    if results.empty:
        return

    # Number of interactions per cell type pair
    pair_counts = results.groupby(['source', 'target']).size().reset_index(name='n_interactions')
    pivot = pair_counts.pivot_table(index='source', columns='target',
                                     values='n_interactions', fill_value=0)

    fig, ax = plt.subplots(figsize=(max(8, pivot.shape[1] * 0.8),
                                     max(6, pivot.shape[0] * 0.6)))
    pivot = pivot.astype(int)
    sns.heatmap(pivot, annot=True, fmt='d', cmap='YlOrRd', ax=ax,
                linewidths=0.5, cbar_kws={'label': 'Number of interactions'})
    ax.set_title(f'Cell-Cell Interactions: {condition_name}')
    plt.tight_layout()
    fig.savefig(PLOT_DIR / f"interaction_heatmap_{condition_name}.png",
                dpi=150, bbox_inches='tight')
    plt.close(fig)

    # Top interactions dot plot
    # Use specificity_rank if available
    rank_col = None
    for col in results.columns:
        if 'specificity_rank' in col:
            rank_col = col
            break
        elif 'magnitude_rank' in col:
            rank_col = col

    if rank_col:
        top = results.nsmallest(30, rank_col)
        top['interaction'] = (top['ligand_complex'].astype(str) + ' → ' +
                              top['receptor_complex'].astype(str))
        top['pair'] = top['source'].astype(str) + ' → ' + top['target'].astype(str)

        fig, ax = plt.subplots(figsize=(12, 8))
        for i, (_, row) in enumerate(top.iterrows()):
            color = 'red' if row.get('is_pain_relevant', False) else 'steelblue'
            ax.scatter(row[rank_col], i, s=100, c=color, edgecolors='black', linewidth=0.5)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels([f"{row['pair']}: {row['interaction']}"
                            for _, row in top.iterrows()], fontsize=7)
        ax.set_xlabel('Rank (lower = more significant)')
        ax.set_title(f'Top 30 Interactions: {condition_name}\n(red = pain-relevant)')
        ax.invert_yaxis()
        plt.tight_layout()
        fig.savefig(PLOT_DIR / f"top_interactions_{condition_name}.png",
                    dpi=150, bbox_inches='tight')
        plt.close(fig)

    print(f"  Plots saved for {condition_name}")


def plot_differential(diff_df):
    """Plot gained/lost interactions in degeneration."""
    if diff_df.empty:
        return

    # Top gained and lost
    if 'rank_diff' in diff_df.columns:
        diff_df_clean = diff_df.dropna(subset=['rank_diff'])
        if diff_df_clean.empty:
            return

        top_gained = diff_df_clean.nsmallest(15, 'rank_diff')  # most negative = gained in degen
        top_lost = diff_df_clean.nlargest(15, 'rank_diff')  # most positive = lost in degen

        combined = pd.concat([top_gained, top_lost])

        fig, ax = plt.subplots(figsize=(10, 8))

        # Use source/target from whichever condition has data
        for col_prefix in ['source', 'target', 'ligand_complex', 'receptor_complex']:
            h_col = f'{col_prefix}_healthy'
            d_col = f'{col_prefix}_degen'
            if h_col in combined.columns and d_col in combined.columns:
                combined[col_prefix] = combined[h_col].fillna(combined[d_col])

        if 'source' in combined.columns:
            labels = (combined.get('source', '').astype(str) + '→' +
                      combined.get('target', '').astype(str) + ': ' +
                      combined.get('ligand_complex', '').astype(str) + '→' +
                      combined.get('receptor_complex', '').astype(str))

            colors = ['#e74c3c' if v < 0 else '#3498db' for v in combined['rank_diff']]
            ax.barh(range(len(combined)), combined['rank_diff'], color=colors)
            ax.set_yticks(range(len(combined)))
            ax.set_yticklabels([l[:60] for l in labels], fontsize=7)
            ax.set_xlabel('Rank difference (negative = gained in degeneration)')
            ax.set_title('Differential Interactions: Healthy vs Degenerated')
            ax.axvline(0, color='black', linewidth=0.5)
            plt.tight_layout()
            fig.savefig(PLOT_DIR / "differential_interactions.png",
                        dpi=150, bbox_inches='tight')
        plt.close(fig)
        print("  Differential interaction plot saved")


# ═══════════════════════════════════════════════════════════════════════════════
# Report Generation
# ═══════════════════════════════════════════════════════════════════════════════

def generate_report(condition_results, diff_df, pain_results):
    """Generate HTML communication report."""
    html = ['<!DOCTYPE html><html><head>',
            '<title>Communication Report — Module 09</title>',
            '<style>',
            '  body { font-family: Arial, sans-serif; max-width: 1400px; margin: 0 auto; padding: 20px; }',
            '  h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }',
            '  h2 { color: #34495e; margin-top: 30px; }',
            '  h3 { color: #7f8c8d; }',
            '  table { border-collapse: collapse; width: 100%; margin: 10px 0; }',
            '  th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }',
            '  th { background-color: #3498db; color: white; }',
            '  tr:nth-child(even) { background-color: #f2f2f2; }',
            '  img { max-width: 100%; height: auto; margin: 10px 0; border: 1px solid #ddd; }',
            '  .stats-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0; }',
            '  .stat-box { background: #ecf0f1; padding: 15px; border-radius: 5px; }',
            '  .stat-box h3 { margin: 0 0 5px 0; color: #7f8c8d; font-size: 12px; }',
            '  .stat-box p { margin: 0; font-size: 24px; font-weight: bold; color: #2c3e50; }',
            '  .pain { background-color: #fde8e8; }',
            '</style></head><body>',
            '<h1>Module 09: Cell-Cell Communication Report</h1>',
            f'<p>Generated: {datetime.now().strftime("%Y-%m-%d %H:%M")}</p>']

    # Stats
    total_interactions = sum(len(r) for r in condition_results.values())
    n_pain = len(pain_results) if not pain_results.empty else 0
    html.append('<div class="stats-grid">')
    html.append(f'<div class="stat-box"><h3>Total interactions</h3><p>{total_interactions:,}</p></div>')
    html.append(f'<div class="stat-box"><h3>Pain-relevant</h3><p>{n_pain}</p></div>')
    html.append(f'<div class="stat-box"><h3>Conditions analyzed</h3><p>{len(condition_results)}</p></div>')
    html.append('</div>')

    # Per condition
    for cond, results in condition_results.items():
        html.append(f'<h2>{cond.title()} Condition</h2>')

        heatmap = PLOT_DIR / f"interaction_heatmap_{cond}.png"
        if heatmap.exists():
            html.append(f'<img src="interaction_plots/interaction_heatmap_{cond}.png">')

        top_plot = PLOT_DIR / f"top_interactions_{cond}.png"
        if top_plot.exists():
            html.append(f'<img src="interaction_plots/top_interactions_{cond}.png">')

    # Differential
    if not diff_df.empty:
        html.append('<h2>Differential Interactions</h2>')
        diff_plot = PLOT_DIR / "differential_interactions.png"
        if diff_plot.exists():
            html.append('<img src="interaction_plots/differential_interactions.png">')

    # Pain interactions
    if not pain_results.empty:
        html.append('<h2>Pain-Relevant Interactions</h2>')
        html.append('<table><tr><th>Source</th><th>Target</th><th>Ligand</th>'
                   '<th>Receptor</th><th>Condition</th><th>Pain Category</th></tr>')
        for _, row in pain_results.head(50).iterrows():
            html.append(f'<tr class="pain"><td>{row.get("source", "")}</td>'
                       f'<td>{row.get("target", "")}</td>'
                       f'<td>{row.get("ligand_complex", "")}</td>'
                       f'<td>{row.get("receptor_complex", "")}</td>'
                       f'<td>{row.get("condition", "")}</td>'
                       f'<td>{row.get("pain_categories", "")}</td></tr>')
        html.append('</table>')

    # Checkpoint
    html.append('<h2>Human Checkpoint Questions</h2><ol>')
    html.append('<li>Are the top interactions biologically plausible?</li>')
    html.append('<li>Do condition-specific interactions suggest mechanisms of degeneration?</li>')
    html.append('<li>Are there pain-relevant interactions that could be therapeutic targets?</li>')
    html.append('<li>Are any interactions likely artifacts of ambient RNA or doublets?</li>')
    html.append('</ol>')

    html.append('</body></html>')
    (RESULTS_DIR / "communication_report.html").write_text('\n'.join(html))
    print(f"  Report saved: {RESULTS_DIR / 'communication_report.html'}")


# ═══════════════════════════════════════════════════════════════════════════════
# Validation
# ═══════════════════════════════════════════════════════════════════════════════

def validate():
    """Automated validation checks."""
    print("\n" + "=" * 60)
    print("Validating Module 09 outputs")
    print("=" * 60)

    checks = []

    # Check interaction results
    for cond in ['healthy', 'degenerated']:
        path = RESULTS_DIR / f"interactions_{cond}.tsv"
        if path.exists():
            df = pd.read_csv(path, sep='\t')
            checks.append(('PASS', f'{cond} interactions: {len(df):,}'))

            # Check for collagen-integrin (positive control)
            has_collagen = df['ligand_complex'].astype(str).str.contains('COL', na=False).any()
            if has_collagen:
                checks.append(('PASS', f'{cond}: collagen interactions found (positive control)'))
        else:
            checks.append(('WARN', f'{cond} interaction results not found'))

    # Pain interactions
    pain_path = RESULTS_DIR / "pain_interactions.tsv"
    if pain_path.exists():
        pain = pd.read_csv(pain_path, sep='\t')
        checks.append(('PASS', f'Pain-relevant interactions flagged: {len(pain):,}'))
    else:
        checks.append(('WARN', 'Pain interactions not found'))

    # Plots
    plot_count = len(list(PLOT_DIR.glob("*.png")))
    if plot_count > 0:
        checks.append(('PASS', f'{plot_count} interaction plots generated'))
    else:
        checks.append(('WARN', 'No interaction plots found'))

    # Report
    if (RESULTS_DIR / "communication_report.html").exists():
        checks.append(('PASS', 'Communication report generated'))
    else:
        checks.append(('FAIL', 'Communication report not generated'))

    passed = all(s != 'FAIL' for s, _ in checks)
    for status, msg in checks:
        print(f"  [{status}] {msg}")
    return passed, checks


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 60)
    print("Module 09: Cell-Cell Communication")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    if '--validate-only' in sys.argv:
        passed, _ = validate()
        sys.exit(0 if passed else 1)

    # Load data
    print("\nLoading data for communication analysis...")
    condition_data = load_data_for_communication()

    if not condition_data:
        print("ERROR: No suitable data for communication analysis")
        sys.exit(1)

    # Run LIANA per condition
    condition_results = {}
    all_pain = []

    for cond_name, adata in condition_data.items():
        results = run_liana(adata, cond_name)
        if not results.empty:
            results = flag_pain_interactions(results)
            condition_results[cond_name] = results
            plot_interaction_summary(results, cond_name)

            pain = results[results.get('is_pain_relevant', False)]
            if not pain.empty:
                all_pain.append(pain)

        del adata
        gc.collect()

    # Compare conditions
    diff_df = pd.DataFrame()
    if 'healthy' in condition_results and 'degenerated' in condition_results:
        diff_df = compare_interactions(condition_results['healthy'],
                                        condition_results['degenerated'])
        plot_differential(diff_df)

    pain_combined = pd.concat(all_pain) if all_pain else pd.DataFrame()

    # Generate report
    print(f"\n{'=' * 60}")
    print("Generating Report")
    print('=' * 60)
    generate_report(condition_results, diff_df, pain_combined)

    # Validate
    passed, _ = validate()
    print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Overall validation: {'PASSED' if passed else 'NEEDS REVIEW'}")


if __name__ == "__main__":
    main()
