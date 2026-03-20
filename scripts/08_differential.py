#!/usr/bin/env python3
"""Module 08: Differential analysis for IVD scRNA-seq atlas.

Part 1: Cell composition analysis (propeller-style proportion tests)
Part 2: Pseudobulk differential expression using pyDESeq2

Usage:
    python3 scripts/08_differential.py
    python3 scripts/08_differential.py --composition-only
    python3 scripts/08_differential.py --de-only
    python3 scripts/08_differential.py --validate-only
"""

import gc
import sys
import os
import warnings
import traceback
from pathlib import Path
from datetime import datetime
from collections import defaultdict

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
from scipy.stats import mannwhitneyu, fisher_exact
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
META_PATH = BASE / "metadata" / "sample_metadata.tsv"
RESULTS_DIR = BASE / "results" / "differential"
VOLCANO_DIR = RESULTS_DIR / "volcano_plots"
HEATMAP_DIR = RESULTS_DIR / "heatmaps"

for d in [RESULTS_DIR, VOLCANO_DIR, HEATMAP_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# ── Configuration ────────────────────────────────────────────────────────────

# Samples to EXCLUDE from DE (per condition mapping review)
EXCLUDE_SAMPLES = {'GSE205535_NNP'}  # 11yo spinal cord injury — trauma confound

# Minimum requirements for pseudobulk DE
MIN_SAMPLES_PER_CONDITION = 3
MIN_CELLS_PER_SAMPLE = 50

# DE significance thresholds
LFC_THRESHOLD = 0.5
PADJ_THRESHOLD = 0.05

# Comparisons to run: (name, condition_A, condition_B)
# condition_A is the reference (denominator in fold change)
COMPARISONS = [
    ("healthy_vs_degenerated_all", "healthy", ["degenerated_mild", "degenerated_severe", "degenerated_ungraded"]),
    ("healthy_vs_degenerated_mild", "healthy", ["degenerated_mild"]),
    ("healthy_vs_degenerated_severe", "healthy", ["degenerated_severe"]),
    ("mild_vs_severe", "degenerated_mild", ["degenerated_severe"]),
    # ("healthy_vs_herniated", "healthy", ["herniated"]),  # Excluded: only GSE251686 contributes herniated samples after GSE233666 removal — fully confounded with study
]

# IVD-relevant genes for highlighting in volcano plots
IVD_GENES = [
    'MMP1', 'MMP3', 'MMP9', 'MMP13', 'ADAMTS4', 'ADAMTS5',
    'COL1A1', 'COL1A2', 'COL2A1', 'COL3A1', 'COL6A1', 'COL9A1', 'COL11A1',
    'ACAN', 'BGN', 'DCN', 'VCAN', 'LUM',
    'IL1B', 'IL6', 'TNF', 'TNFRSF11B', 'CCL2', 'CXCL8',
    'SOX9', 'FOXF1', 'KRT19', 'CD24', 'COMP',
    'NGF', 'BDNF', 'CALCA', 'TAC1', 'TRPV1',  # Pain-related
    'HIF1A', 'VEGFA', 'CA12',  # Hypoxia
    'TP53', 'CDKN1A', 'CDKN2A',  # Senescence
]


# ═══════════════════════════════════════════════════════════════════════════════
# DATA LOADING
# ═══════════════════════════════════════════════════════════════════════════════

def load_integrated_data():
    """Load integrated h5ad files and merge with sample metadata."""
    meta = pd.read_csv(META_PATH, sep='\t')

    datasets = {}
    for fname, label in [
        ("NP.h5ad", "NP"),
        ("AF.h5ad", "AF"),
        ("CEP.h5ad", "CEP"),
    ]:
        path = INT_DIR / fname
        if not path.exists():
            print(f"  WARNING: {path} not found, skipping")
            continue
        a = sc.read_h5ad(path)
        # Ensure counts layer exists
        if 'counts' not in a.layers:
            if sparse.issparse(a.X):
                a.layers['counts'] = a.X.copy()
            else:
                a.layers['counts'] = sparse.csr_matrix(a.X)
        datasets[label] = a
        print(f"  Loaded {label}: {a.shape[0]:,} cells × {a.shape[1]:,} genes")

    return datasets, meta


def filter_to_comparison(obs, meta, condition_ref, condition_test, exclude=EXCLUDE_SAMPLES):
    """Filter obs to cells in the comparison, return merged DataFrame."""
    # Normalize condition_test to list
    if isinstance(condition_test, str):
        condition_test = [condition_test]

    # Merge with metadata
    merged = obs.merge(meta[['sample_id', 'condition_harmonized', 'study_accession']],
                       on='sample_id', how='left', suffixes=('', '_meta'))

    # Use metadata condition if available, fall back to obs
    if 'condition_harmonized_meta' in merged.columns:
        merged['condition'] = merged['condition_harmonized_meta'].fillna(merged.get('condition_harmonized', ''))
    else:
        merged['condition'] = merged['condition_harmonized']

    # Exclude flagged samples
    merged = merged[~merged['sample_id'].isin(exclude)]

    # Filter to comparison conditions
    ref_mask = merged['condition'] == condition_ref
    test_mask = merged['condition'].isin(condition_test)
    merged = merged[ref_mask | test_mask].copy()

    # Assign binary group label
    merged['group'] = 'reference'
    merged.loc[test_mask[merged.index], 'group'] = 'test'

    return merged


# ═══════════════════════════════════════════════════════════════════════════════
# PART 1: COMPOSITION ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════════

def run_composition_analysis(datasets, meta):
    """Test for differential cell type proportions between conditions."""
    print("\n" + "=" * 60)
    print("Part 1: Cell Composition Analysis")
    print("=" * 60)

    all_results = []

    # Combine obs from all datasets
    all_obs = []
    for label, adata in datasets.items():
        obs = adata.obs[['sample_id', 'cell_type']].copy()
        obs['tier'] = label
        all_obs.append(obs)
    combined_obs = pd.concat(all_obs, axis=0)

    for comp_name, cond_ref, cond_test in COMPARISONS:
        print(f"\n  --- {comp_name} ---")
        merged = filter_to_comparison(combined_obs.reset_index(), meta, cond_ref, cond_test)

        if merged.empty:
            print(f"    No cells for this comparison, skipping")
            continue

        # Compute proportions per sample
        ct_counts = merged.groupby(['sample_id', 'cell_type', 'group']).size().reset_index(name='n_cells')
        sample_totals = merged.groupby('sample_id').size().reset_index(name='total')
        ct_counts = ct_counts.merge(sample_totals, on='sample_id')
        ct_counts['proportion'] = ct_counts['n_cells'] / ct_counts['total']

        # Test each cell type
        for ct in ct_counts['cell_type'].unique():
            ct_data = ct_counts[ct_counts['cell_type'] == ct]
            ref_props = ct_data[ct_data['group'] == 'reference']['proportion'].values
            test_props = ct_data[ct_data['group'] == 'test']['proportion'].values

            if len(ref_props) < 2 or len(test_props) < 2:
                continue

            # Mann-Whitney U test on proportions
            try:
                stat, pval = mannwhitneyu(ref_props, test_props, alternative='two-sided')
            except ValueError:
                continue

            lfc = np.log2((np.mean(test_props) + 1e-6) / (np.mean(ref_props) + 1e-6))

            all_results.append({
                'comparison': comp_name,
                'cell_type': ct,
                'mean_prop_ref': np.mean(ref_props),
                'mean_prop_test': np.mean(test_props),
                'log2FC_proportion': lfc,
                'n_samples_ref': len(ref_props),
                'n_samples_test': len(test_props),
                'pvalue': pval,
            })

    if not all_results:
        print("  No composition results generated")
        return pd.DataFrame()

    results_df = pd.DataFrame(all_results)

    # FDR correction (Benjamini-Hochberg)
    from statsmodels.stats.multitest import multipletests
    _, results_df['padj'], _, _ = multipletests(results_df['pvalue'], method='fdr_bh')

    # Sort and save
    results_df = results_df.sort_values('padj')
    results_df.to_csv(RESULTS_DIR / "composition_analysis.tsv", sep='\t', index=False)
    print(f"\n  Saved: {RESULTS_DIR / 'composition_analysis.tsv'}")
    print(f"  Significant (padj < 0.05): {(results_df['padj'] < 0.05).sum()} / {len(results_df)}")

    # Plot
    _plot_composition(results_df, combined_obs, meta)

    return results_df


def _plot_composition(results_df, combined_obs, meta):
    """Generate composition analysis plots."""
    sig = results_df[results_df['padj'] < 0.05]
    if sig.empty:
        print("  No significant composition changes to plot")
        return

    # Box plots of significant cell type proportions by condition
    for _, row in sig.iterrows():
        comp_name = row['comparison']
        ct = row['cell_type']
        comp_info = next((c for c in COMPARISONS if c[0] == comp_name), None)
        if comp_info is None:
            continue

        _, cond_ref, cond_test = comp_info
        merged = filter_to_comparison(combined_obs.reset_index(), meta, cond_ref, cond_test)
        ct_data = merged[merged['cell_type'] == ct]
        sample_props = ct_data.groupby(['sample_id', 'group']).size().reset_index(name='n')
        totals = merged.groupby('sample_id').size().reset_index(name='total')
        sample_props = sample_props.merge(totals, on='sample_id')
        sample_props['proportion'] = sample_props['n'] / sample_props['total']

        fig, ax = plt.subplots(figsize=(4, 5))
        sns.boxplot(data=sample_props, x='group', y='proportion', ax=ax)
        sns.stripplot(data=sample_props, x='group', y='proportion', ax=ax,
                      color='black', size=4, alpha=0.6)
        ax.set_title(f"{ct}\n{comp_name}\npadj={row['padj']:.3g}")
        ax.set_ylabel("Proportion")
        plt.tight_layout()
        fig.savefig(RESULTS_DIR / f"composition_{ct}_{comp_name}.png", dpi=150, bbox_inches='tight')
        plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════════
# PART 2: PSEUDOBULK DIFFERENTIAL EXPRESSION
# ═══════════════════════════════════════════════════════════════════════════════

def pseudobulk_aggregate(adata, obs_merged, cell_type):
    """Aggregate raw counts per sample for a given cell type."""
    # Filter to cell type
    ct_mask = obs_merged['cell_type'] == cell_type
    ct_obs = obs_merged[ct_mask]

    # Count cells per sample
    cells_per_sample = ct_obs.groupby('sample_id').size()
    valid_samples = cells_per_sample[cells_per_sample >= MIN_CELLS_PER_SAMPLE].index.tolist()

    if len(valid_samples) == 0:
        return None, None

    # Aggregate counts
    counts_list = []
    sample_info = []
    for sample_id in valid_samples:
        sample_mask = ct_obs['sample_id'] == sample_id
        cell_indices = ct_obs.index[sample_mask]

        # Get positions in adata
        pos = np.where(adata.obs.index.isin(cell_indices))[0]
        if len(pos) == 0:
            continue

        # Sum counts across cells
        if 'counts' in adata.layers:
            counts_matrix = adata.layers['counts'][pos]
        else:
            counts_matrix = adata.X[pos]

        if sparse.issparse(counts_matrix):
            sample_counts = np.asarray(counts_matrix.sum(axis=0)).flatten()
        else:
            sample_counts = counts_matrix.sum(axis=0).flatten()

        counts_list.append(sample_counts)
        group = ct_obs.loc[ct_obs['sample_id'] == sample_id, 'group'].iloc[0]
        study = ct_obs.loc[ct_obs['sample_id'] == sample_id, 'study_accession'].iloc[0]
        sample_info.append({
            'sample_id': sample_id,
            'group': group,
            'study': study,
            'n_cells': len(pos),
        })

    if not counts_list:
        return None, None

    # Build count matrix (samples × genes)
    count_matrix = np.vstack(counts_list).astype(int)
    sample_df = pd.DataFrame(sample_info)

    # Create DataFrame with gene names
    counts_df = pd.DataFrame(count_matrix, columns=adata.var_names, index=sample_df['sample_id'])

    return counts_df, sample_df


def run_deseq2(counts_df, sample_df, cell_type, comparison_name):
    """Run pyDESeq2 on pseudobulk counts."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    # Check power
    ref_n = (sample_df['group'] == 'reference').sum()
    test_n = (sample_df['group'] == 'test').sum()
    if ref_n < MIN_SAMPLES_PER_CONDITION or test_n < MIN_SAMPLES_PER_CONDITION:
        return None, f"Underpowered: {ref_n} ref, {test_n} test (need >= {MIN_SAMPLES_PER_CONDITION})"

    # Filter lowly expressed genes: keep genes detected in >= 10% of samples
    gene_detection = (counts_df > 0).mean(axis=0)
    keep_genes = gene_detection[gene_detection >= 0.1].index
    counts_filtered = counts_df[keep_genes]

    if counts_filtered.shape[1] < 100:
        return None, f"Too few genes after filtering: {counts_filtered.shape[1]}"

    # Prepare metadata
    # DE: RNA assay; correct for study (batch) and IVD score (continuous covariate)
    cols_to_keep = ['group', 'study']
    if 'ivd_score' in sample_df.columns:
        cols_to_keep.append('ivd_score')
    clinical = sample_df.set_index('sample_id')[cols_to_keep].copy()

    # If all samples from same study, don't include study in design
    n_studies = clinical['study'].nunique()
    if n_studies > 1 and all(clinical.groupby('group')['study'].nunique() > 1):
        design_factors = ['group', 'study']
    else:
        design_factors = ['group']

    # Include ivd_score as continuous covariate when available and variable
    if 'ivd_score' in clinical.columns:
        ivd_vals = pd.to_numeric(clinical['ivd_score'], errors='coerce')
        if ivd_vals.notna().sum() >= 2 and ivd_vals.nunique() > 1:
            clinical['ivd_score'] = ivd_vals
            design_factors.append('ivd_score')
        else:
            clinical = clinical.drop(columns=['ivd_score'])

    try:
        dds = DeseqDataSet(
            counts=counts_filtered,
            metadata=clinical,
            design="+".join(design_factors) if len(design_factors) > 1 else "group",
        )
        dds.deseq2()

        stat_res = DeseqStats(dds, contrast=["group", "test", "reference"])
        stat_res.summary()

        results = stat_res.results_df.copy()
        results['gene'] = results.index
        results = results.rename(columns={
            'log2FoldChange': 'log2FC',
            'pvalue': 'pvalue',
            'padj': 'padj',
            'baseMean': 'baseMean',
        })

        # Add metadata
        results['cell_type'] = cell_type
        results['comparison'] = comparison_name
        results['n_ref'] = ref_n
        results['n_test'] = test_n
        results['n_genes_tested'] = len(keep_genes)
        results['design'] = "+".join(design_factors)

        return results, None

    except Exception as e:
        return None, f"DESeq2 error: {str(e)}"


def run_de_analysis(datasets, meta):
    """Run pseudobulk DE for all cell types × comparisons."""
    print("\n" + "=" * 60)
    print("Part 2: Pseudobulk Differential Expression")
    print("=" * 60)

    all_de_results = []
    skipped = []

    for tier_label, adata in datasets.items():
        print(f"\n  === Tier: {tier_label} ===")

        # Get cell types in this tier
        cell_types = adata.obs['cell_type'].unique()

        for comp_name, cond_ref, cond_test in COMPARISONS:
            # Build merged obs for this comparison
            obs_reset = adata.obs[['sample_id', 'cell_type']].copy()
            obs_reset.index = adata.obs.index  # preserve index for alignment
            merged = filter_to_comparison(obs_reset.reset_index(), meta, cond_ref, cond_test)

            if merged.empty:
                continue

            # Set index back for alignment with adata
            if 'index' in merged.columns:
                merged = merged.set_index('index')

            for ct in cell_types:
                label = f"{ct} | {comp_name}"

                # Pseudobulk aggregate
                counts_df, sample_df = pseudobulk_aggregate(adata, merged, ct)
                if counts_df is None:
                    skipped.append((label, "No valid samples with >= 50 cells"))
                    continue

                # Run DESeq2
                print(f"    {label}: {sample_df.shape[0]} samples...", end=" ", flush=True)
                results, skip_reason = run_deseq2(counts_df, sample_df, ct, comp_name)

                if results is not None:
                    n_sig = ((results['padj'] < PADJ_THRESHOLD) & (results['log2FC'].abs() > LFC_THRESHOLD)).sum()
                    print(f"{n_sig} DE genes")
                    all_de_results.append(results)

                    # Save individual results
                    safe_ct = ct.replace(' ', '_').replace('/', '_')
                    outpath = RESULTS_DIR / "de_results" / f"{safe_ct}_{comp_name}.tsv"
                    outpath.parent.mkdir(parents=True, exist_ok=True)
                    results.to_csv(outpath, sep='\t', index=False)
                else:
                    print(f"SKIPPED: {skip_reason}")
                    skipped.append((label, skip_reason))

                gc.collect()

    # Log skipped comparisons
    if skipped:
        skip_df = pd.DataFrame(skipped, columns=['comparison', 'reason'])
        skip_df.to_csv(RESULTS_DIR / "skipped_comparisons.tsv", sep='\t', index=False)
        print(f"\n  Skipped {len(skipped)} comparisons (see skipped_comparisons.tsv)")

    if not all_de_results:
        print("  No DE results generated")
        return pd.DataFrame()

    # Combine all results
    combined = pd.concat(all_de_results, ignore_index=True)
    combined.to_csv(RESULTS_DIR / "de_results_combined.tsv", sep='\t', index=False)
    print(f"\n  Combined DE results: {combined.shape[0]:,} gene-tests")

    sig = combined[(combined['padj'] < PADJ_THRESHOLD) & (combined['log2FC'].abs() > LFC_THRESHOLD)]
    print(f"  Significant (|log2FC| > {LFC_THRESHOLD}, padj < {PADJ_THRESHOLD}): {sig.shape[0]:,}")

    # Summary table
    summary = sig.groupby(['cell_type', 'comparison']).agg(
        n_up=('log2FC', lambda x: (x > 0).sum()),
        n_down=('log2FC', lambda x: (x < 0).sum()),
        n_total=('log2FC', 'count'),
    ).reset_index()
    summary.to_csv(RESULTS_DIR / "de_summary_table.tsv", sep='\t', index=False)
    print(f"\n  DE Summary:")
    print(summary.to_string(index=False))

    # Generate plots
    _plot_volcanos(combined)
    _plot_heatmap(combined, datasets)

    return combined


def _plot_volcanos(combined):
    """Generate volcano plots for each cell type × comparison with DE genes."""
    sig_combos = combined.groupby(['cell_type', 'comparison']).apply(
        lambda x: ((x['padj'] < PADJ_THRESHOLD) & (x['log2FC'].abs() > LFC_THRESHOLD)).any()
    )
    sig_combos = sig_combos[sig_combos].index.tolist()

    for ct, comp in sig_combos:
        data = combined[(combined['cell_type'] == ct) & (combined['comparison'] == comp)].copy()
        data = data.dropna(subset=['log2FC', 'padj'])
        if data.empty:
            continue

        data['neg_log10_padj'] = -np.log10(data['padj'].clip(lower=1e-300))

        fig, ax = plt.subplots(figsize=(8, 6))

        # Color by significance
        sig_up = (data['padj'] < PADJ_THRESHOLD) & (data['log2FC'] > LFC_THRESHOLD)
        sig_down = (data['padj'] < PADJ_THRESHOLD) & (data['log2FC'] < -LFC_THRESHOLD)
        ns = ~(sig_up | sig_down)

        ax.scatter(data.loc[ns, 'log2FC'], data.loc[ns, 'neg_log10_padj'],
                   c='gray', s=3, alpha=0.3, label='NS')
        ax.scatter(data.loc[sig_up, 'log2FC'], data.loc[sig_up, 'neg_log10_padj'],
                   c='red', s=5, alpha=0.5, label=f'Up ({sig_up.sum()})')
        ax.scatter(data.loc[sig_down, 'log2FC'], data.loc[sig_down, 'neg_log10_padj'],
                   c='blue', s=5, alpha=0.5, label=f'Down ({sig_down.sum()})')

        # Highlight IVD genes
        for gene in IVD_GENES:
            gene_data = data[data['gene'] == gene]
            if not gene_data.empty:
                row = gene_data.iloc[0]
                if row['padj'] < PADJ_THRESHOLD:
                    ax.annotate(gene, (row['log2FC'], row['neg_log10_padj']),
                                fontsize=7, ha='center', va='bottom',
                                fontweight='bold', color='darkred' if row['log2FC'] > 0 else 'darkblue')

        ax.axhline(-np.log10(PADJ_THRESHOLD), ls='--', c='gray', lw=0.5)
        ax.axvline(LFC_THRESHOLD, ls='--', c='gray', lw=0.5)
        ax.axvline(-LFC_THRESHOLD, ls='--', c='gray', lw=0.5)
        ax.set_xlabel('log2 Fold Change')
        ax.set_ylabel('-log10(adjusted p-value)')
        ax.set_title(f'{ct}\n{comp}')
        ax.legend(loc='upper right', fontsize=8)
        plt.tight_layout()

        safe_ct = ct.replace(' ', '_').replace('/', '_')
        fig.savefig(VOLCANO_DIR / f"volcano_{safe_ct}_{comp}.png", dpi=150, bbox_inches='tight')
        plt.close(fig)

    print(f"  Volcano plots saved to {VOLCANO_DIR}")


def _plot_heatmap(combined, datasets):
    """Generate heatmap of top DE genes for the primary comparison."""
    # Focus on the primary comparison: healthy_vs_degenerated_all
    primary = combined[combined['comparison'] == 'healthy_vs_degenerated_all']
    if primary.empty:
        return

    for ct in primary['cell_type'].unique():
        ct_data = primary[primary['cell_type'] == ct]
        sig = ct_data[(ct_data['padj'] < PADJ_THRESHOLD) & (ct_data['log2FC'].abs() > LFC_THRESHOLD)]
        if sig.empty or sig.shape[0] < 5:
            continue

        # Top 50 genes by adjusted p-value
        top_genes = sig.nsmallest(50, 'padj')['gene'].tolist()

        # Find the adata containing this cell type
        target_adata = None
        for label, adata in datasets.items():
            if ct in adata.obs['cell_type'].values:
                target_adata = adata
                break

        if target_adata is None:
            continue

        # Get pseudobulk expression for top genes (log-normalized X, averaged per sample)
        ct_mask = target_adata.obs['cell_type'] == ct
        ct_adata = target_adata[ct_mask]
        available_genes = [g for g in top_genes if g in ct_adata.var_names]
        if len(available_genes) < 5:
            continue

        # Average expression per sample
        expr_data = {}
        meta_full = pd.read_csv(META_PATH, sep='\t')
        for sample_id in ct_adata.obs['sample_id'].unique():
            if sample_id in EXCLUDE_SAMPLES:
                continue
            sample_mask = (ct_adata.obs['sample_id'] == sample_id).values  # .values for numpy bool array
            n_cells = sample_mask.sum()
            if n_cells < MIN_CELLS_PER_SAMPLE:
                continue

            gene_idx = [ct_adata.var_names.get_loc(g) for g in available_genes]
            X_sub = ct_adata.X[sample_mask][:, gene_idx]
            if sparse.issparse(X_sub):
                X_sub = X_sub.toarray()
            expr_data[sample_id] = X_sub.mean(axis=0)

        if len(expr_data) < 4:
            continue

        expr_df = pd.DataFrame(expr_data, index=available_genes).T

        # Add condition annotation
        sample_conditions = meta_full.set_index('sample_id')['condition_harmonized']
        conditions = expr_df.index.map(lambda x: sample_conditions.get(x, 'unknown'))

        # Color map for conditions
        cond_colors = {
            'healthy': '#2ecc71', 'degenerated_mild': '#f39c12',
            'degenerated_severe': '#e74c3c', 'degenerated_ungraded': '#e67e22',
            'herniated': '#9b59b6', 'neonatal': '#3498db',
            'aged_ungraded': '#95a5a6',
        }
        row_colors = conditions.map(lambda x: cond_colors.get(x, 'gray'))

        try:
            g = sns.clustermap(expr_df, row_colors=row_colors.values,
                               cmap='RdBu_r', center=0, figsize=(12, max(6, len(expr_df) * 0.3)),
                               xticklabels=True, yticklabels=True,
                               dendrogram_ratio=0.1, cbar_pos=(0.02, 0.8, 0.03, 0.15))
            g.ax_heatmap.set_title(f'Top DE genes: {ct}\nhealthy_vs_degenerated_all', fontsize=10)
            g.savefig(HEATMAP_DIR / f"heatmap_{ct.replace(' ', '_')}_healthy_vs_deg.png",
                      dpi=150, bbox_inches='tight')
            plt.close()
        except Exception as e:
            print(f"    Heatmap error for {ct}: {e}")

    print(f"  Heatmaps saved to {HEATMAP_DIR}")


# ═══════════════════════════════════════════════════════════════════════════════
# VALIDATION
# ═══════════════════════════════════════════════════════════════════════════════

def validate():
    """Run validation checks per spec."""
    print("\n" + "=" * 60)
    print("Validating Module 08 outputs")
    print("=" * 60)

    messages = []

    # Composition results exist
    comp_path = RESULTS_DIR / "composition_analysis.tsv"
    if comp_path.exists():
        messages.append("PASS: Composition analysis results exist")
    else:
        messages.append("FAIL: Composition analysis results missing")

    # DE results exist
    de_dir = RESULTS_DIR / "de_results"
    if de_dir.exists():
        de_files = list(de_dir.glob("*.tsv"))
        messages.append(f"PASS: {len(de_files)} DE result files generated")
    else:
        messages.append("FAIL: No DE results directory")

    # Check for inflated statistics (>5000 sig genes)
    combined_path = RESULTS_DIR / "de_results_combined.tsv"
    if combined_path.exists():
        combined = pd.read_csv(combined_path, sep='\t')
        for (ct, comp), group in combined.groupby(['cell_type', 'comparison']):
            sig = group[(group['padj'] < 0.05) & (group['log2FC'].abs() > 0.5)]
            if sig.shape[0] > 5000:
                messages.append(f"WARNING: {ct} | {comp}: {sig.shape[0]} significant genes (>5000, possible inflation)")
            else:
                messages.append(f"PASS: {ct} | {comp}: {sig.shape[0]} significant genes")

        # Check for positive controls (MMP13, ADAMTS5 in NP degeneration)
        np_deg = combined[(combined['cell_type'].str.contains('NP_mature')) &
                          (combined['comparison'] == 'healthy_vs_degenerated_all')]
        if not np_deg.empty:
            for gene in ['MMP13', 'ADAMTS5']:
                gene_row = np_deg[np_deg['gene'] == gene]
                if not gene_row.empty:
                    row = gene_row.iloc[0]
                    if row['padj'] < 0.05 and row['log2FC'] > 0:
                        messages.append(f"PASS: Positive control {gene} upregulated in NP degeneration (log2FC={row['log2FC']:.2f}, padj={row['padj']:.2g})")
                    else:
                        messages.append(f"INFO: Positive control {gene} in NP degeneration: log2FC={row['log2FC']:.2f}, padj={row['padj']:.2g}")
                else:
                    messages.append(f"INFO: Positive control {gene} not found in NP results")

    # Volcano plots exist
    volcanos = list(VOLCANO_DIR.glob("*.png"))
    if volcanos:
        messages.append(f"PASS: {len(volcanos)} volcano plots generated")
    else:
        messages.append("WARNING: No volcano plots generated")

    # Skipped comparisons logged
    skip_path = RESULTS_DIR / "skipped_comparisons.tsv"
    if skip_path.exists():
        skipped = pd.read_csv(skip_path, sep='\t')
        messages.append(f"PASS: {len(skipped)} skipped comparisons logged")

    # Print results
    all_pass = True
    for msg in messages:
        prefix = "PASS" if msg.startswith("PASS") else ("FAIL" if msg.startswith("FAIL") else "INFO")
        if msg.startswith("FAIL"):
            all_pass = False
        print(f"  [{prefix}] {msg}")

    return all_pass, messages


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 60)
    print("Module 08: Differential Analysis")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    args = sys.argv[1:]
    composition_only = '--composition-only' in args
    de_only = '--de-only' in args
    validate_only = '--validate-only' in args

    if validate_only:
        passed, messages = validate()
        sys.exit(0 if passed else 1)

    # Load data
    print("\nLoading integrated data...")
    datasets, meta = load_integrated_data()

    if not datasets:
        print("ERROR: No integrated datasets found")
        sys.exit(1)

    # Part 1: Composition
    if not de_only:
        comp_results = run_composition_analysis(datasets, meta)

    # Part 2: DE
    if not composition_only:
        de_results = run_de_analysis(datasets, meta)

    # Validate
    passed, messages = validate()

    print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Overall validation: {'PASSED' if passed else 'FAILED'}")


if __name__ == "__main__":
    main()
