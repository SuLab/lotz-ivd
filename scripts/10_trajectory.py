#!/usr/bin/env python3
"""Module 10: Trajectory & dynamics analysis for IVD scRNA-seq atlas.

Part 1: Trajectory inference (PAGA + DPT pseudotime)
Part 2: RNA velocity check (document availability)
Part 3: Gene programs along trajectories

Usage:
    python3 scripts/10_trajectory.py
    python3 scripts/10_trajectory.py --validate-only
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
from scipy import sparse, stats
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
RESULTS_DIR = BASE / "results" / "trajectories"
DE_DIR = BASE / "results" / "differential"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ── Configuration ────────────────────────────────────────────────────────────
# Cell types to include per compartment (mesenchymal/resident only)
# These will be updated after Module 07 de novo annotation assigns final labels.
# For now, filter to mesenchymal cells and let the de novo labels drive the analysis.
NP_RESIDENT = None  # Will use all mesenchymal cell types in NP object
AF_RESIDENT = None  # Will use all mesenchymal cell types in AF object

# Embeddings to use
PRIMARY_EMBEDDING = 'X_integrated'    # CCA integrated embedding
SENSITIVITY_EMBEDDING = None  # No second embedding in new pipeline

PSEUDOTIME_CORR_FDR = 0.05
N_TOP_TRAJECTORY_GENES = 500
DOWNSAMPLE_FOR_PAGA = 50000  # Max cells for PAGA computation


# ═══════════════════════════════════════════════════════════════════════════════
# Data Loading
# ═══════════════════════════════════════════════════════════════════════════════

def load_compartment(compartment):
    """Load integrated data for a compartment.

    Subsets to mesenchymal-tier cells — for tiered_v4 this is the cells with
    `obs['tier'] == 'mesenchymal'` (covers both cell_class=='mesenchymal' AND
    the 'unknown'-class residual that 05m routed into the mes integration
    tier). For v5 the `tier` column may not exist, so we fall back to
    cell_class=='mesenchymal'.
    """
    path = INT_DIR / f"{compartment}.h5ad"

    if not path.exists():
        print(f"  WARNING: {path} not found")
        return None

    print(f"  Loading {path.name}...")
    adata = sc.read_h5ad(path)
    print(f"    Full: {adata.n_obs:,} cells × {adata.n_vars} genes")

    if 'tier' in adata.obs.columns:
        mask = adata.obs['tier'].astype(str) == 'mesenchymal'
        adata = adata[mask].copy()
        print(f"    Mesenchymal tier only: {adata.n_obs:,} cells")
    elif 'cell_class' in adata.obs.columns:
        mask = adata.obs['cell_class'] == 'mesenchymal'
        adata = adata[mask].copy()
        print(f"    cell_class=='mesenchymal' only: {adata.n_obs:,} cells")

    print(f"    Cell types: {dict(adata.obs['cell_type'].value_counts())}")

    return adata


def load_metadata():
    """Load sample metadata."""
    meta = pd.read_csv(META_PATH, sep='\t')
    return meta


# ═══════════════════════════════════════════════════════════════════════════════
# Part 1: Trajectory Inference
# ═══════════════════════════════════════════════════════════════════════════════

def run_trajectory(adata, compartment, embedding_key, embedding_name):
    """Run PAGA + DPT trajectory analysis."""
    print(f"\n  === {compartment} trajectory ({embedding_name}) ===")

    if embedding_key not in adata.obsm:
        print(f"    WARNING: {embedding_key} not found, skipping")
        return None

    # Use the integration embedding for neighbors
    emb = adata.obsm[embedding_key]
    n_dims = min(emb.shape[1], 30)
    print(f"    Using {embedding_key} ({n_dims} dims) for neighbor graph")

    # Downsample if needed for memory
    if adata.n_obs > DOWNSAMPLE_FOR_PAGA:
        print(f"    Downsampling from {adata.n_obs:,} to {DOWNSAMPLE_FOR_PAGA:,} for PAGA")
        sc.pp.subsample(adata, n_obs=DOWNSAMPLE_FOR_PAGA, random_state=42, copy=False)
        emb = adata.obsm[embedding_key]

    # Compute neighbors on the embedding
    sc.pp.neighbors(adata, use_rep=embedding_key, n_neighbors=15, random_state=42)

    # Leiden clustering at moderate resolution for PAGA
    sc.tl.leiden(adata, resolution=0.5, random_state=42, key_added='leiden_traj')
    print(f"    Leiden clusters: {adata.obs['leiden_traj'].nunique()}")

    # PAGA
    sc.tl.paga(adata, groups='leiden_traj')
    # Must plot PAGA to compute layout positions needed for UMAP init
    sc.pl.paga(adata, show=False)
    plt.close('all')

    # Compute UMAP using PAGA initialization
    sc.tl.umap(adata, init_pos='paga', random_state=42)

    # Diffusion map for pseudotime
    sc.tl.diffmap(adata, n_comps=15)

    # Select root cell: highest proportion of notochordal/healthy cells
    root_cluster = _select_root(adata, compartment)
    print(f"    Root cluster: {root_cluster}")

    # Set root cell
    adata.uns['iroot'] = _get_root_cell(adata, root_cluster)

    # DPT
    sc.tl.dpt(adata, n_dcs=10)
    print(f"    Pseudotime range: [{adata.obs['dpt_pseudotime'].min():.3f}, "
          f"{adata.obs['dpt_pseudotime'].max():.3f}]")

    # Save pseudotime
    pt_df = adata.obs[['cell_type', 'dpt_pseudotime']].copy()
    if 'condition_harmonized' in adata.obs.columns:
        pt_df['condition'] = adata.obs['condition_harmonized']
    if 'study_id' in adata.obs.columns:
        pt_df['study'] = adata.obs['study_id']

    suffix = f"_{embedding_name}" if embedding_name != 'scVI_mesenchymal' else ""
    pt_df.to_csv(RESULTS_DIR / f"pseudotime_{compartment}{suffix}.tsv", sep='\t')

    return adata


def _select_root(adata, compartment):
    """Select root cluster for DPT based on biology."""
    ct_col = 'cell_type'

    if compartment == 'NP':
        # Root = notochordal cells (most primitive/healthy)
        preferred = ['NP_notochordal', 'NP_mature_chondrocyte']
    elif compartment == 'AF':
        # Root = AF_inner (less exposed to mechanical stress)
        preferred = ['AF_inner', 'AF_outer']
    else:
        # CEP or other — use healthy-enriched cluster as fallback
        preferred = []

    # Find leiden cluster with highest proportion of preferred cell types
    cross = pd.crosstab(adata.obs['leiden_traj'], adata.obs[ct_col], normalize='index')
    for pref in preferred:
        if pref in cross.columns:
            best_cluster = cross[pref].idxmax()
            pref_prop = cross.loc[best_cluster, pref]
            if pref_prop > 0.3:
                print(f"    Root selection: cluster {best_cluster} "
                      f"({pref_prop:.0%} {pref})")
                return best_cluster

    # Fallback: cluster with most healthy cells
    if 'condition_harmonized' in adata.obs.columns:
        cond_cross = pd.crosstab(adata.obs['leiden_traj'],
                                  adata.obs['condition_harmonized'], normalize='index')
        for cond in ['healthy', 'neonatal']:
            if cond in cond_cross.columns:
                best = cond_cross[cond].idxmax()
                if cond_cross.loc[best, cond] > 0.2:
                    print(f"    Root selection (fallback): cluster {best} "
                          f"({cond_cross.loc[best, cond]:.0%} {cond})")
                    return best

    # Last resort: cluster 0
    print("    Root selection: defaulting to cluster 0")
    return '0'


def _get_root_cell(adata, cluster):
    """Get the index of the root cell (most central in cluster)."""
    mask = adata.obs['leiden_traj'] == cluster
    if not mask.any():
        return 0

    # Use the diffmap center of the cluster
    cluster_cells = np.where(mask)[0]
    if 'X_diffmap' in adata.obsm:
        diffmap = adata.obsm['X_diffmap'][cluster_cells]
        center = diffmap.mean(axis=0)
        dists = np.sum((diffmap - center) ** 2, axis=1)
        return cluster_cells[np.argmin(dists)]
    return cluster_cells[0]


def plot_trajectory(adata, compartment, embedding_name):
    """Generate trajectory visualization plots."""
    suffix = f"_{embedding_name}" if embedding_name != 'scVI_mesenchymal' else ""

    fig, axes = plt.subplots(2, 2, figsize=(16, 14))

    # 1. UMAP by cell type
    sc.pl.umap(adata, color='cell_type', ax=axes[0, 0], show=False,
               title=f'{compartment} Cell Types')

    # 2. UMAP by pseudotime
    sc.pl.umap(adata, color='dpt_pseudotime', ax=axes[0, 1], show=False,
               title=f'{compartment} Pseudotime', cmap='viridis')

    # 3. UMAP by condition
    if 'condition_harmonized' in adata.obs.columns:
        sc.pl.umap(adata, color='condition_harmonized', ax=axes[1, 0], show=False,
                   title=f'{compartment} Condition')
    else:
        axes[1, 0].set_visible(False)

    # 4. UMAP by study
    if 'study_id' in adata.obs.columns:
        sc.pl.umap(adata, color='study_id', ax=axes[1, 1], show=False,
                   title=f'{compartment} Study')
    else:
        axes[1, 1].set_visible(False)

    plt.tight_layout()
    fig.savefig(RESULTS_DIR / f"umap_trajectory_{compartment}{suffix}.png",
                dpi=150, bbox_inches='tight')
    plt.close(fig)

    # PAGA graph
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    sc.pl.paga(adata, ax=axes[0], show=False, title=f'{compartment} PAGA (Leiden)')
    sc.tl.paga(adata, groups='cell_type')
    sc.pl.paga(adata, ax=axes[1], show=False, title=f'{compartment} PAGA (Cell Type)')
    plt.tight_layout()
    fig.savefig(RESULTS_DIR / f"paga_{compartment}{suffix}.png",
                dpi=150, bbox_inches='tight')
    plt.close(fig)

    # Pseudotime distributions by condition
    if 'condition_harmonized' in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(10, 6))
        pt_data = adata.obs[['dpt_pseudotime', 'condition_harmonized']].dropna()
        conditions = pt_data['condition_harmonized'].unique()
        for cond in sorted(conditions):
            vals = pt_data[pt_data['condition_harmonized'] == cond]['dpt_pseudotime']
            ax.hist(vals, bins=50, alpha=0.5, label=f'{cond} (n={len(vals)})', density=True)
        ax.set_xlabel('Pseudotime')
        ax.set_ylabel('Density')
        ax.set_title(f'{compartment} Pseudotime by Condition')
        ax.legend()
        plt.tight_layout()
        fig.savefig(RESULTS_DIR / f"pseudotime_by_condition_{compartment}{suffix}.png",
                    dpi=150, bbox_inches='tight')
        plt.close(fig)

    # Pseudotime by cell type (violin)
    if adata.obs['cell_type'].nunique() > 1:
        fig, ax = plt.subplots(figsize=(10, 6))
        sc.pl.violin(adata, keys='dpt_pseudotime', groupby='cell_type',
                     ax=ax, show=False, rotation=45)
        ax.set_title(f'{compartment} Pseudotime by Cell Type')
        plt.tight_layout()
        fig.savefig(RESULTS_DIR / f"pseudotime_by_celltype_{compartment}{suffix}.png",
                    dpi=150, bbox_inches='tight')
        plt.close(fig)

    print(f"    Plots saved for {compartment} ({embedding_name})")


def pseudotime_condition_correlation(adata, compartment, embedding_name):
    """Test if pseudotime correlates with condition."""
    suffix = f"_{embedding_name}" if embedding_name != 'scVI_mesenchymal' else ""
    results = []

    if 'condition_harmonized' not in adata.obs.columns:
        return pd.DataFrame()

    pt = adata.obs['dpt_pseudotime'].values
    cond = adata.obs['condition_harmonized'].values

    # Encode condition as ordinal: healthy=0, mild=1, severe=2
    cond_order = {'neonatal': -1, 'healthy': 0, 'degenerated_mild': 1,
                  'degenerated_ungraded': 1.5, 'degenerated_severe': 2,
                  'aged_ungraded': 1, 'herniated': 2}

    ordinal = np.array([cond_order.get(c, np.nan) for c in cond])
    valid = ~np.isnan(ordinal) & ~np.isnan(pt)

    if valid.sum() > 10:
        rho, pval = stats.spearmanr(pt[valid], ordinal[valid])
        results.append({
            'compartment': compartment, 'embedding': embedding_name,
            'test': 'pseudotime_vs_condition_ordinal',
            'rho': rho, 'pvalue': pval, 'n': valid.sum(),
        })
        print(f"    Pseudotime vs condition: rho={rho:.3f}, p={pval:.1e}")

    # Mann-Whitney: healthy vs degenerated
    healthy_mask = (cond == 'healthy') & ~np.isnan(pt)
    degen_mask = np.isin(cond, ['degenerated_mild', 'degenerated_severe', 'degenerated_ungraded']) & ~np.isnan(pt)

    if healthy_mask.sum() > 5 and degen_mask.sum() > 5:
        stat, pval = stats.mannwhitneyu(pt[healthy_mask], pt[degen_mask], alternative='two-sided')
        results.append({
            'compartment': compartment, 'embedding': embedding_name,
            'test': 'pseudotime_healthy_vs_degenerated',
            'statistic': stat, 'pvalue': pval,
            'n_healthy': healthy_mask.sum(), 'n_degen': degen_mask.sum(),
            'median_healthy': np.median(pt[healthy_mask]),
            'median_degen': np.median(pt[degen_mask]),
        })
        print(f"    Healthy vs degenerated pseudotime: "
              f"median {np.median(pt[healthy_mask]):.3f} vs {np.median(pt[degen_mask]):.3f}, "
              f"p={pval:.1e}")

    # Per cell type
    for ct in adata.obs['cell_type'].unique():
        ct_mask = adata.obs['cell_type'] == ct
        ct_pt = pt[ct_mask]
        ct_cond = cond[ct_mask]

        ct_ordinal = np.array([cond_order.get(c, np.nan) for c in ct_cond])
        ct_valid = ~np.isnan(ct_ordinal) & ~np.isnan(ct_pt)
        if ct_valid.sum() > 10:
            rho, pval = stats.spearmanr(ct_pt[ct_valid], ct_ordinal[ct_valid])
            results.append({
                'compartment': compartment, 'embedding': embedding_name,
                'test': f'pseudotime_vs_condition_{ct}',
                'rho': rho, 'pvalue': pval, 'n': ct_valid.sum(),
            })

    if results:
        res_df = pd.DataFrame(results)
        res_df.to_csv(RESULTS_DIR / f"pseudotime_correlations_{compartment}{suffix}.tsv",
                      sep='\t', index=False)
        return res_df
    return pd.DataFrame()


# ═══════════════════════════════════════════════════════════════════════════════
# Part 2: RNA Velocity Check
# ═══════════════════════════════════════════════════════════════════════════════

def check_velocity_data(adata, compartment):
    """Check if spliced/unspliced counts are available."""
    print(f"\n  === RNA Velocity Data Check ({compartment}) ===")

    has_spliced = 'spliced' in adata.layers
    has_unspliced = 'unspliced' in adata.layers

    print(f"    Spliced layer: {'YES' if has_spliced else 'NO'}")
    print(f"    Unspliced layer: {'YES' if has_unspliced else 'NO'}")
    print(f"    Available layers: {list(adata.layers.keys())}")

    if has_spliced and has_unspliced:
        print("    RNA velocity CAN be computed")
        return True
    else:
        print("    RNA velocity CANNOT be computed (no spliced/unspliced counts)")
        print("    This is expected for most public scRNA-seq datasets that only provide")
        print("    final count matrices. Would require reprocessing from BAM/fastq files.")
        return False


# ═══════════════════════════════════════════════════════════════════════════════
# Part 3: Trajectory Gene Programs
# ═══════════════════════════════════════════════════════════════════════════════

def find_trajectory_genes(adata, compartment, embedding_name):
    """Find genes correlated with pseudotime."""
    suffix = f"_{embedding_name}" if embedding_name != 'scVI_mesenchymal' else ""
    print(f"\n  Finding trajectory-associated genes ({compartment}, {embedding_name})...")

    pt = adata.obs['dpt_pseudotime'].values
    valid = ~np.isnan(pt)
    if valid.sum() < 100:
        print("    Too few cells with valid pseudotime")
        return pd.DataFrame()

    pt_valid = pt[valid]

    # Get expression matrix
    X = adata.X[valid]
    if sparse.issparse(X):
        X = X.toarray()

    results = []
    for i, gene in enumerate(adata.var_names):
        expr = X[:, i]
        if expr.std() == 0:
            continue
        rho, pval = stats.spearmanr(pt_valid, expr)
        results.append({'gene': gene, 'rho': rho, 'pvalue': pval})

    if not results:
        return pd.DataFrame()

    gene_df = pd.DataFrame(results)

    # FDR correction
    from scipy.stats import false_discovery_control
    gene_df['padj'] = false_discovery_control(gene_df['pvalue'].values, method='bh')

    sig = gene_df[gene_df['padj'] < PSEUDOTIME_CORR_FDR].copy()
    sig = sig.reindex(sig['rho'].abs().nlargest(N_TOP_TRAJECTORY_GENES).index)

    # Assign temporal programs
    sig['program'] = pd.cut(sig['rho'],
                            bins=[-1, -0.1, 0.1, 1],
                            labels=['late_down', 'stable', 'late_up'])

    sig.to_csv(RESULTS_DIR / f"trajectory_genes_{compartment}{suffix}.tsv",
               sep='\t', index=False)
    print(f"    Found {len(sig)} trajectory-associated genes (top {N_TOP_TRAJECTORY_GENES})")
    print(f"    Programs: {dict(sig['program'].value_counts())}")

    return sig


def plot_trajectory_genes(adata, traj_genes, compartment, embedding_name):
    """Plot key gene expression along pseudotime."""
    suffix = f"_{embedding_name}" if embedding_name != 'scVI_mesenchymal' else ""

    if traj_genes.empty:
        return

    # Key IVD genes to highlight
    ivd_genes = ['ACAN', 'COL2A1', 'COL1A1', 'MMP13', 'MMP3', 'ADAMTS5',
                 'SOX9', 'TNF', 'IL6', 'CDKN1A', 'CDKN2A', 'SERPINE1',
                 'TBXT', 'KRT19', 'FN1', 'COMP', 'VEGFA', 'HIF1A',
                 'CXCL1', 'CXCL2', 'CXCL3', 'CEMIP']

    # Top trajectory genes + IVD genes that are in the data
    available_ivd = [g for g in ivd_genes if g in adata.var_names]
    top_traj = traj_genes.head(10)['gene'].tolist()
    genes_to_plot = list(dict.fromkeys(available_ivd + top_traj))[:20]

    pt = adata.obs['dpt_pseudotime'].values
    valid = ~np.isnan(pt)

    n_genes = len(genes_to_plot)
    n_cols = 4
    n_rows = (n_genes + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 3 * n_rows))
    axes = axes.flatten()

    for i, gene in enumerate(genes_to_plot):
        if gene not in adata.var_names:
            axes[i].set_visible(False)
            continue

        ax = axes[i]
        idx = adata.var_names.get_loc(gene)
        expr = adata.X[valid, idx]
        if sparse.issparse(expr):
            expr = expr.toarray().flatten()
        else:
            expr = np.asarray(expr).flatten()

        pt_v = pt[valid]

        # Scatter with lowess trend
        ax.scatter(pt_v, expr, s=1, alpha=0.1, c='gray')

        # Bin and compute mean
        bins = np.linspace(pt_v.min(), pt_v.max(), 50)
        bin_idx = np.digitize(pt_v, bins)
        bin_means = [expr[bin_idx == b].mean() if (bin_idx == b).any() else np.nan
                     for b in range(1, len(bins))]
        ax.plot(bins[:-1] + np.diff(bins)/2, bin_means, 'r-', linewidth=2)

        # Correlation info
        in_traj = traj_genes[traj_genes['gene'] == gene]
        if not in_traj.empty:
            rho = in_traj.iloc[0]['rho']
            ax.set_title(f'{gene} (rho={rho:.2f})', fontsize=9)
        else:
            ax.set_title(gene, fontsize=9)
        ax.set_xlabel('Pseudotime', fontsize=8)
        ax.set_ylabel('Expression', fontsize=8)
        ax.tick_params(labelsize=7)

    for i in range(len(genes_to_plot), len(axes)):
        axes[i].set_visible(False)

    plt.suptitle(f'{compartment} Gene Expression Along Pseudotime ({embedding_name})',
                 fontsize=14)
    plt.tight_layout()
    fig.savefig(RESULTS_DIR / f"gene_dynamics_{compartment}{suffix}.png",
                dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Gene dynamics plot saved")


def cross_reference_de(traj_genes, compartment, embedding_name):
    """Cross-reference trajectory genes with DE results."""
    suffix = f"_{embedding_name}" if embedding_name != 'scVI_mesenchymal' else ""

    de_path = DE_DIR / "de_results_combined.tsv"
    if not de_path.exists() or traj_genes.empty:
        return

    de = pd.read_csv(de_path, sep='\t')
    sig_de = de[(de['padj'] < 0.05) & (de['log2FC'].abs() > 0.5)]

    traj_set = set(traj_genes['gene'])
    de_set = set(sig_de['gene'])
    overlap = traj_set & de_set

    print(f"\n    Trajectory genes: {len(traj_set)}")
    print(f"    DE genes: {len(de_set)}")
    print(f"    Overlap: {len(overlap)}")

    if overlap:
        overlap_df = traj_genes[traj_genes['gene'].isin(overlap)].merge(
            sig_de[sig_de['gene'].isin(overlap)][['gene', 'cell_type', 'comparison', 'log2FC', 'padj']],
            on='gene', how='inner'
        )
        overlap_df.to_csv(RESULTS_DIR / f"trajectory_de_overlap_{compartment}{suffix}.tsv",
                          sep='\t', index=False)

    # Venn-like summary
    summary = {
        'trajectory_only': len(traj_set - de_set),
        'de_only': len(de_set - traj_set),
        'both': len(overlap),
    }
    print(f"    Trajectory only: {summary['trajectory_only']}, "
          f"DE only: {summary['de_only']}, Both: {summary['both']}")

    return summary


# ═══════════════════════════════════════════════════════════════════════════════
# Report Generation
# ═══════════════════════════════════════════════════════════════════════════════

def generate_report(compartments_run, velocity_available):
    """Generate HTML trajectory report."""
    html = ['<!DOCTYPE html><html><head>',
            '<title>Trajectory Report — Module 10</title>',
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
            '  .section { margin: 20px 0; padding: 15px; background: #f8f9fa; border-radius: 5px; border-left: 4px solid #3498db; }',
            '  .warning { border-left-color: #e74c3c; }',
            '</style></head><body>',
            '<h1>Module 10: Trajectory & Dynamics Report</h1>',
            f'<p>Generated: {datetime.now().strftime("%Y-%m-%d %H:%M")}</p>']

    for comp in compartments_run:
        html.append(f'<h2>{comp} Compartment</h2>')

        # UMAP + pseudotime
        umap_path = f"umap_trajectory_{comp}.png"
        if (RESULTS_DIR / umap_path).exists():
            html.append(f'<img src="{umap_path}" alt="{comp} trajectory UMAP">')

        # PAGA
        paga_path = f"paga_{comp}.png"
        if (RESULTS_DIR / paga_path).exists():
            html.append(f'<img src="{paga_path}" alt="{comp} PAGA">')

        # Pseudotime by condition
        pt_cond = f"pseudotime_by_condition_{comp}.png"
        if (RESULTS_DIR / pt_cond).exists():
            html.append(f'<img src="{pt_cond}" alt="{comp} pseudotime by condition">')

        # Pseudotime by cell type
        pt_ct = f"pseudotime_by_celltype_{comp}.png"
        if (RESULTS_DIR / pt_ct).exists():
            html.append(f'<img src="{pt_ct}" alt="{comp} pseudotime by cell type">')

        # Gene dynamics
        gene_dyn = f"gene_dynamics_{comp}.png"
        if (RESULTS_DIR / gene_dyn).exists():
            html.append(f'<img src="{gene_dyn}" alt="{comp} gene dynamics">')

        # Correlation results
        corr_path = RESULTS_DIR / f"pseudotime_correlations_{comp}.tsv"
        if corr_path.exists():
            corr = pd.read_csv(corr_path, sep='\t')
            html.append('<h3>Pseudotime Correlations</h3><table>')
            html.append('<tr><th>Test</th><th>Statistic</th><th>P-value</th><th>N</th></tr>')
            for _, row in corr.iterrows():
                stat_val = f"rho={row['rho']:.3f}" if pd.notna(row.get('rho')) else f"U={row.get('statistic', 0):.0f}"
                n_val = row.get('n', 0)
                if pd.isna(n_val):
                    n_val = int(row.get('n_healthy', 0) or 0) + int(row.get('n_degen', 0) or 0)
                html.append(f'<tr><td>{row["test"]}</td><td>{stat_val}</td>'
                           f'<td>{row["pvalue"]:.1e}</td><td>{int(n_val)}</td></tr>')
            html.append('</table>')

        # Trajectory genes
        traj_path = RESULTS_DIR / f"trajectory_genes_{comp}.tsv"
        if traj_path.exists():
            traj = pd.read_csv(traj_path, sep='\t')
            html.append(f'<h3>Trajectory Genes: {len(traj)} total</h3>')
            html.append('<table><tr><th>Gene</th><th>rho</th><th>padj</th><th>Program</th></tr>')
            for _, row in traj.head(20).iterrows():
                html.append(f'<tr><td>{row["gene"]}</td><td>{row["rho"]:.3f}</td>'
                           f'<td>{row["padj"]:.1e}</td><td>{row["program"]}</td></tr>')
            html.append('</table>')

        # (Sensitivity check removed — single scVI embedding per new pipeline)

    # RNA velocity
    html.append('<h2>RNA Velocity</h2>')
    if any(velocity_available.values()):
        html.append('<p>Velocity data available for some datasets. Results below.</p>')
    else:
        html.append('<div class="section warning">')
        html.append('<strong>RNA velocity could not be computed.</strong> ')
        html.append('No spliced/unspliced count layers found in the integrated data. ')
        html.append('This is expected for public scRNA-seq datasets that only provide final count matrices. ')
        html.append('Computing velocity would require reprocessing from BAM/fastq files using velocyto or STARsolo. ')
        html.append('Given that IVD is a slow-turnover tissue, velocity results may be uninformative even if computed.')
        html.append('</div>')

    # Checkpoint
    html.append('<h2>Human Checkpoint Questions</h2><ol>')
    html.append('<li>Does the inferred trajectory make biological sense?</li>')
    html.append('<li>Does pseudotime align with healthy → degenerated?</li>')
    html.append('<li>Are there branch points suggesting divergent cell fates?</li>')
    html.append('<li>Do gene programs along the trajectory reveal staged degeneration?</li>')
    html.append('<li>Should trajectory findings redefine cell states?</li>')
    html.append('</ol>')

    html.append('</body></html>')
    (RESULTS_DIR / "trajectory_report.html").write_text('\n'.join(html))
    print(f"  Report saved: {RESULTS_DIR / 'trajectory_report.html'}")


# ═══════════════════════════════════════════════════════════════════════════════
# Validation
# ═══════════════════════════════════════════════════════════════════════════════

def validate():
    """Run automated validation checks."""
    print("\n" + "=" * 60)
    print("Validating Module 10 outputs")
    print("=" * 60)

    checks = []

    for comp in ['NP', 'AF', 'CEP']:
        pt_path = RESULTS_DIR / f"pseudotime_{comp}.tsv"
        if pt_path.exists():
            pt = pd.read_csv(pt_path, sep='\t')
            checks.append(('PASS', f'{comp} pseudotime: {len(pt):,} cells'))

            # Check correlation with condition
            corr_path = RESULTS_DIR / f"pseudotime_correlations_{comp}.tsv"
            if corr_path.exists():
                corr = pd.read_csv(corr_path, sep='\t')
                sig_corr = corr[corr['pvalue'] < 0.05]
                if not sig_corr.empty:
                    checks.append(('PASS', f'{comp} pseudotime correlates with biological variable'))
                else:
                    checks.append(('WARN', f'{comp} pseudotime does not significantly correlate with any variable'))
        else:
            checks.append(('WARN', f'{comp} pseudotime not found'))

        traj_path = RESULTS_DIR / f"trajectory_genes_{comp}.tsv"
        if traj_path.exists():
            traj = pd.read_csv(traj_path, sep='\t')
            checks.append(('PASS', f'{comp} trajectory genes: {len(traj)}'))
        else:
            checks.append(('WARN', f'{comp} trajectory genes not found'))

    report = RESULTS_DIR / "trajectory_report.html"
    if report.exists():
        checks.append(('PASS', 'Trajectory report generated'))
    else:
        checks.append(('FAIL', 'Trajectory report not generated'))

    passed = all(s != 'FAIL' for s, _ in checks)
    for status, msg in checks:
        print(f"  [{status}] {msg}")

    return passed, checks


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    global INT_DIR, RESULTS_DIR, DE_DIR

    args = sys.argv[1:]

    # Parse --input-dir (where to read h5ads from)
    for i, a in enumerate(args):
        if a == '--input-dir' and i + 1 < len(args):
            INT_DIR = Path(args[i + 1]).resolve()

    # Parse --output-dir (where Module 10 results land)
    for i, a in enumerate(args):
        if a == '--output-dir' and i + 1 < len(args):
            RESULTS_DIR = Path(args[i + 1]).resolve()
            RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Parse --de-dir (where Module 08 DE results live, for cross-reference)
    for i, a in enumerate(args):
        if a == '--de-dir' and i + 1 < len(args):
            DE_DIR = Path(args[i + 1]).resolve()

    print("=" * 60)
    print("Module 10: Trajectory & Dynamics Analysis")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  Input dir:  {INT_DIR}")
    print(f"  Output dir: {RESULTS_DIR}")
    print(f"  DE dir:     {DE_DIR}")
    print("=" * 60)

    if '--validate-only' in args:
        passed, _ = validate()
        sys.exit(0 if passed else 1)

    meta = load_metadata()
    compartments_run = []
    velocity_available = {}

    for compartment in ['NP', 'AF', 'CEP']:
        print(f"\n{'=' * 60}")
        print(f"Processing {compartment}")
        print('=' * 60)

        adata = load_compartment(compartment)
        if adata is None:
            continue

        # Check velocity data
        velocity_available[compartment] = check_velocity_data(adata, compartment)

        # Primary trajectory (scVI mesenchymal)
        adata_traj = run_trajectory(adata, compartment, PRIMARY_EMBEDDING, 'scVI_mesenchymal')
        if adata_traj is not None:
            plot_trajectory(adata_traj, compartment, 'scVI_mesenchymal')
            pseudotime_condition_correlation(adata_traj, compartment, 'scVI_mesenchymal')
            traj_genes = find_trajectory_genes(adata_traj, compartment, 'scVI_mesenchymal')
            plot_trajectory_genes(adata_traj, traj_genes, compartment, 'scVI_mesenchymal')
            cross_reference_de(traj_genes, compartment, 'scVI_mesenchymal')
            compartments_run.append(compartment)

        del adata_traj
        gc.collect()

        del adata
        gc.collect()

    # Generate report
    print(f"\n{'=' * 60}")
    print("Generating Report")
    print('=' * 60)
    generate_report(compartments_run, velocity_available)

    # Validate
    passed, _ = validate()
    print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Overall validation: {'PASSED' if passed else 'NEEDS REVIEW'}")


if __name__ == "__main__":
    main()
