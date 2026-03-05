#!/usr/bin/env python3
"""Module 07: Biological interpretation for IVD scRNA-seq atlas.

Part 1: Pathway and gene set enrichment (ORA + GSEA via gseapy/decoupler)
Part 2: TF activity inference (decoupler, CollecTRI regulons)
Part 3: Pain-associated gene analysis

Usage:
    python3 scripts/07_interpretation.py
    python3 scripts/07_interpretation.py --pathways-only
    python3 scripts/07_interpretation.py --pain-only
    python3 scripts/07_interpretation.py --tf-only
    python3 scripts/07_interpretation.py --validate-only
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
from scipy.stats import fisher_exact, false_discovery_control
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent.parent
DE_DIR = BASE / "results" / "differential"
INT_DIR = BASE / "data" / "integrated"
RESULTS_DIR = BASE / "results" / "interpretation"
PATHWAY_DIR = RESULTS_DIR / "pathway_enrichment"
TF_DIR = RESULTS_DIR / "tf_activity"
PAIN_DIR = RESULTS_DIR

for d in [RESULTS_DIR, PATHWAY_DIR, TF_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# ── Configuration ────────────────────────────────────────────────────────────
PADJ_THRESHOLD = 0.05
LFC_THRESHOLD = 0.5
ENRICHMENT_FDR = 0.05
MIN_GENE_SET_SIZE = 10
MAX_GENE_SET_SIZE = 500

# Comparisons to EXCLUDE from primary interpretation (flagged as confounded)
CONFOUNDED_COMPARISONS = {'healthy_vs_herniated'}

# ── IVD-Specific Gene Sets ──────────────────────────────────────────────────
IVD_GENE_SETS = {
    'ECM_homeostasis': [
        'COL1A1', 'COL1A2', 'COL2A1', 'COL9A1', 'COL9A2', 'COL9A3',
        'COL11A1', 'COL11A2', 'ACAN', 'VCAN', 'HAPLN1', 'BGN', 'DCN',
        'LUM', 'FMOD', 'FN1', 'COMP', 'THBS1', 'THBS2', 'SPARC',
    ],
    'Matrix_degradation': [
        'MMP1', 'MMP2', 'MMP3', 'MMP7', 'MMP9', 'MMP13', 'MMP14',
        'ADAMTS4', 'ADAMTS5', 'ADAMTSL1', 'TIMP1', 'TIMP2', 'TIMP3',
    ],
    'Inflammatory_signaling': [
        'IL1B', 'IL6', 'TNF', 'CXCL1', 'CXCL2', 'CXCL3', 'CXCL8',
        'CCL2', 'CCL5', 'PTGS2', 'NFKB1', 'NFKB2', 'RELA', 'RELB',
        'IL1A', 'IL18', 'TGFB1', 'TGFB2', 'TGFB3',
    ],
    'Cellular_senescence': [
        'CDKN1A', 'CDKN2A', 'CDKN2B', 'TP53', 'SERPINE1', 'GLB1',
        'IL1A', 'IL6', 'CXCL8', 'MMP1', 'MMP3', 'IGFBP3', 'IGFBP5',
        'IGFBP7', 'CCL2', 'GDF15',
    ],
    'Oxidative_stress_hypoxia': [
        'HIF1A', 'EPAS1', 'VEGFA', 'SOD1', 'SOD2', 'SOD3', 'CAT',
        'GPX1', 'GPX4', 'NOS2', 'NOS3', 'HMOX1', 'NFE2L2', 'NQO1',
    ],
    'Apoptosis': [
        'BAX', 'BAK1', 'BCL2', 'BCL2L1', 'MCL1', 'CASP3', 'CASP8',
        'CASP9', 'FASLG', 'FAS', 'PMAIP1', 'BID', 'BIRC5', 'XIAP',
    ],
    'Autophagy': [
        'BECN1', 'MAP1LC3B', 'ATG5', 'ATG7', 'ATG12', 'SQSTM1',
        'ULK1', 'LAMP1', 'LAMP2', 'BNIP3', 'BNIP3L',
    ],
    'Mechanotransduction': [
        'YAP1', 'WWTR1', 'PIEZO1', 'PIEZO2', 'TRPV4', 'ITGA5',
        'ITGB1', 'ITGAV', 'CTGF', 'CYR61', 'ANKRD1',
    ],
    'Notochordal_developmental': [
        'TBXT', 'SHH', 'NOG', 'WNT3A', 'WNT5A', 'WNT11', 'BMP2',
        'BMP4', 'BMP7', 'SOX9', 'SOX5', 'SOX6', 'FOXA1', 'FOXA2',
        'KRT8', 'KRT18', 'KRT19', 'CD24',
    ],
}

# Pain gene sets (from spec)
PAIN_GENE_SETS = {
    'Nociception_ion_channels': [
        'TRPV1', 'TRPV4', 'TRPA1', 'SCN9A', 'SCN10A', 'SCN11A',
        'ASIC1', 'ASIC2', 'ASIC3', 'P2RX3', 'P2RX4', 'P2RX7',
    ],
    'Neuropeptides': [
        'TAC1', 'CALCA', 'CALCB', 'NPY', 'VIP', 'BDNF', 'NGF',
        'PENK', 'PDYN', 'GAL',
    ],
    'Neurotrophin_receptors': [
        'NTRK1', 'NTRK2', 'NTRK3', 'NGFR', 'OPRD1', 'OPRM1', 'OPRK1',
    ],
    'Neurotrophin_signaling': [
        'NGF', 'BDNF', 'NTF3', 'NTF4', 'NTRK1', 'NTRK2', 'NTRK3', 'NGFR',
    ],
    'Nerve_guidance': [
        'SEMA3A', 'SEMA3F', 'SEMA3E', 'NTN1', 'NTN4',
        'NRP1', 'NRP2', 'ROBO1', 'ROBO2', 'SLIT1', 'SLIT2', 'SLIT3',
        'DCC', 'UNC5B',
    ],
    'Inflammatory_pain': [
        'PTGS2', 'PTGES', 'PLA2G2A', 'IL1B', 'IL6', 'TNF', 'CCL2',
        'CXCL8', 'BDKRB1', 'BDKRB2', 'KLK1',
    ],
    'Neovascularization': [
        'VEGFA', 'VEGFB', 'VEGFC', 'FGF2', 'PDGFA', 'PDGFB',
        'ANGPT1', 'ANGPT2', 'TEK', 'KDR', 'FLT1',
    ],
}


# ═══════════════════════════════════════════════════════════════════════════════
# Part 1: Pathway Enrichment
# ═══════════════════════════════════════════════════════════════════════════════

def load_de_results():
    """Load combined DE results from Module 06."""
    path = DE_DIR / "de_results_combined.tsv"
    if not path.exists():
        print(f"ERROR: DE results not found at {path}")
        sys.exit(1)
    df = pd.read_csv(path, sep='\t')
    print(f"  Loaded {len(df):,} DE results")

    # Filter out confounded comparisons for primary analysis
    primary = df[~df['comparison'].isin(CONFOUNDED_COMPARISONS)]
    confounded = df[df['comparison'].isin(CONFOUNDED_COMPARISONS)]
    print(f"  Primary comparisons: {len(primary):,} gene-tests")
    print(f"  Confounded (excluded from primary): {len(confounded):,} gene-tests")
    return primary, confounded


def run_pathway_enrichment(de_results):
    """Run ORA and GSEA on DE results."""
    import gseapy as gp

    sig = de_results[(de_results['padj'] < PADJ_THRESHOLD) &
                     (de_results['log2FC'].abs() > LFC_THRESHOLD)]

    all_enrichment = []
    comparisons = sig.groupby(['cell_type', 'comparison'])

    for (ct, comp), group in comparisons:
        ct_clean = ct.replace('/', '_').replace(' ', '_')
        print(f"  {ct} | {comp}: {len(group)} DE genes")

        up_genes = group[group['log2FC'] > 0]['gene'].tolist()
        down_genes = group[group['log2FC'] < 0]['gene'].tolist()

        # Get background genes for this cell type
        ct_all = de_results[(de_results['cell_type'] == ct) &
                            (de_results['comparison'] == comp)]
        background = ct_all['gene'].tolist()

        for direction, genes in [('up', up_genes), ('down', down_genes)]:
            if len(genes) < 5:
                continue

            # ORA against GO, KEGG, Reactome via Enrichr
            for lib in ['GO_Biological_Process_2023', 'KEGG_2021_Human',
                        'Reactome_2022', 'MSigDB_Hallmark_2020']:
                try:
                    enr = gp.enrich(
                        gene_list=genes,
                        gene_sets=lib,
                        background=background if len(background) > 100 else None,
                        outdir=None,
                        no_plot=True,
                        verbose=False,
                    )
                    if enr.results is not None and len(enr.results) > 0:
                        res = enr.results.copy()
                        res['cell_type'] = ct
                        res['comparison'] = comp
                        res['direction'] = direction
                        res['library'] = lib
                        all_enrichment.append(res)
                except Exception as e:
                    pass  # Some libraries may not be available offline

            # Also test IVD custom gene sets
            try:
                combined_ivd = {**IVD_GENE_SETS, **PAIN_GENE_SETS}
                enr = gp.enrich(
                    gene_list=genes,
                    gene_sets=combined_ivd,
                    background=background if len(background) > 100 else None,
                    outdir=None,
                    no_plot=True,
                    verbose=False,
                )
                if enr.results is not None and len(enr.results) > 0:
                    res = enr.results.copy()
                    res['cell_type'] = ct
                    res['comparison'] = comp
                    res['direction'] = direction
                    res['library'] = 'IVD_custom'
                    all_enrichment.append(res)
            except Exception:
                pass

    if all_enrichment:
        enrichment_df = pd.concat(all_enrichment, ignore_index=True)
        enrichment_df.to_csv(PATHWAY_DIR / "all_enrichment_results.tsv",
                             sep='\t', index=False)
        print(f"\n  Total enrichment results: {len(enrichment_df):,}")

        sig_enr = enrichment_df[enrichment_df['Adjusted P-value'] < ENRICHMENT_FDR]
        print(f"  Significant (FDR < {ENRICHMENT_FDR}): {len(sig_enr):,}")
        return enrichment_df
    else:
        print("  WARNING: No enrichment results obtained")
        return pd.DataFrame()


def run_gsea_ranked(de_results):
    """Run pre-ranked GSEA on full ranked gene lists."""
    import gseapy as gp

    all_gsea = []
    comparisons = de_results.groupby(['cell_type', 'comparison'])

    for (ct, comp), group in comparisons:
        # Only run for comparisons with enough genes
        if len(group) < 100:
            continue

        # Rank by signed significance: sign(LFC) * -log10(pvalue)
        ranked = group[['gene', 'log2FC', 'pvalue']].dropna().copy()
        ranked['rank_metric'] = np.sign(ranked['log2FC']) * -np.log10(ranked['pvalue'].clip(1e-300))
        ranked = ranked.sort_values('rank_metric', ascending=False)
        rank_series = ranked.set_index('gene')['rank_metric']

        ct_clean = ct.replace('/', '_').replace(' ', '_')
        print(f"  GSEA: {ct} | {comp} ({len(rank_series)} genes)")

        for lib in ['GO_Biological_Process_2023', 'MSigDB_Hallmark_2020']:
            try:
                gsea_res = gp.prerank(
                    rnk=rank_series,
                    gene_sets=lib,
                    outdir=None,
                    no_plot=True,
                    verbose=False,
                    min_size=MIN_GENE_SET_SIZE,
                    max_size=MAX_GENE_SET_SIZE,
                    permutation_num=1000,
                    seed=42,
                )
                if gsea_res.res2d is not None and len(gsea_res.res2d) > 0:
                    res = gsea_res.res2d.copy()
                    res['cell_type'] = ct
                    res['comparison'] = comp
                    res['library'] = lib
                    all_gsea.append(res)
            except Exception as e:
                print(f"    GSEA failed for {lib}: {e}")

        # Custom IVD sets
        try:
            combined_ivd = {**IVD_GENE_SETS, **PAIN_GENE_SETS}
            # Filter sets to genes present in data
            filtered_sets = {k: [g for g in v if g in rank_series.index]
                            for k, v in combined_ivd.items()}
            filtered_sets = {k: v for k, v in filtered_sets.items() if len(v) >= 5}

            if filtered_sets:
                gsea_res = gp.prerank(
                    rnk=rank_series,
                    gene_sets=filtered_sets,
                    outdir=None,
                    no_plot=True,
                    verbose=False,
                    min_size=5,
                    max_size=500,
                    permutation_num=1000,
                    seed=42,
                )
                if gsea_res.res2d is not None and len(gsea_res.res2d) > 0:
                    res = gsea_res.res2d.copy()
                    res['cell_type'] = ct
                    res['comparison'] = comp
                    res['library'] = 'IVD_custom'
                    all_gsea.append(res)
        except Exception:
            pass

    if all_gsea:
        gsea_df = pd.concat(all_gsea, ignore_index=True)
        gsea_df.to_csv(PATHWAY_DIR / "gsea_results.tsv", sep='\t', index=False)
        print(f"\n  Total GSEA results: {len(gsea_df):,}")

        sig_gsea = gsea_df[gsea_df['FDR q-val'].astype(float) < ENRICHMENT_FDR]
        print(f"  Significant (FDR < {ENRICHMENT_FDR}): {len(sig_gsea):,}")
        return gsea_df
    else:
        print("  WARNING: No GSEA results obtained")
        return pd.DataFrame()


def plot_enrichment_summary(enrichment_df):
    """Plot top enriched pathways per cell type."""
    if enrichment_df.empty:
        return

    sig = enrichment_df[enrichment_df['Adjusted P-value'] < ENRICHMENT_FDR].copy()
    if sig.empty:
        print("  No significant enrichment results to plot")
        return

    # Top pathways per cell type (by combined score or p-value)
    for ct in sig['cell_type'].unique():
        ct_data = sig[sig['cell_type'] == ct].copy()
        ct_clean = ct.replace('/', '_').replace(' ', '_')

        for direction in ['up', 'down']:
            dir_data = ct_data[ct_data['direction'] == direction]
            if dir_data.empty:
                continue

            # Top 20 by adjusted p-value
            top = dir_data.nsmallest(20, 'Adjusted P-value')

            fig, ax = plt.subplots(figsize=(10, max(4, len(top) * 0.4)))
            top_sorted = top.sort_values('Adjusted P-value', ascending=False)

            colors = {'GO_Biological_Process_2023': '#3498db',
                      'KEGG_2021_Human': '#e74c3c',
                      'Reactome_2022': '#2ecc71',
                      'MSigDB_Hallmark_2020': '#9b59b6',
                      'IVD_custom': '#f39c12'}

            barcolors = [colors.get(lib, '#95a5a6') for lib in top_sorted['library']]
            ax.barh(range(len(top_sorted)), -np.log10(top_sorted['Adjusted P-value']),
                    color=barcolors)
            ax.set_yticks(range(len(top_sorted)))
            ax.set_yticklabels([t[:60] for t in top_sorted['Term']], fontsize=8)
            ax.set_xlabel('-log10(adjusted p-value)')
            ax.set_title(f'{ct} — {direction}regulated pathways')
            ax.axvline(-np.log10(ENRICHMENT_FDR), color='gray', linestyle='--', alpha=0.5)

            # Legend
            from matplotlib.patches import Patch
            handles = [Patch(color=c, label=l) for l, c in colors.items()
                       if l in top_sorted['library'].values]
            if handles:
                ax.legend(handles=handles, loc='lower right', fontsize=7)

            plt.tight_layout()
            fig.savefig(PATHWAY_DIR / f"enrichment_{ct_clean}_{direction}.png",
                        dpi=150, bbox_inches='tight')
            plt.close(fig)

    print(f"  Enrichment plots saved to {PATHWAY_DIR}")


def plot_gsea_summary(gsea_df):
    """Plot GSEA IVD custom gene set results as a heatmap."""
    if gsea_df.empty:
        return

    ivd = gsea_df[gsea_df['library'] == 'IVD_custom'].copy()
    if ivd.empty:
        return

    # Pivot to NES matrix
    ivd['label'] = ivd['cell_type'] + ' | ' + ivd['comparison']
    ivd['FDR q-val'] = ivd['FDR q-val'].astype(float)
    ivd['NES'] = ivd['NES'].astype(float)

    pivot = ivd.pivot_table(index='Term', columns='label', values='NES', aggfunc='first')
    fdr_pivot = ivd.pivot_table(index='Term', columns='label', values='FDR q-val', aggfunc='first')

    if pivot.empty or pivot.shape[0] < 2:
        return

    fig, ax = plt.subplots(figsize=(max(8, pivot.shape[1] * 1.2),
                                     max(4, pivot.shape[0] * 0.5)))
    sns.heatmap(pivot, cmap='RdBu_r', center=0, ax=ax,
                annot=True, fmt='.2f', annot_kws={'fontsize': 7},
                linewidths=0.5, cbar_kws={'label': 'NES'})

    # Mark significant cells
    for i, row in enumerate(pivot.index):
        for j, col in enumerate(pivot.columns):
            fdr = fdr_pivot.loc[row, col] if pd.notna(fdr_pivot.loc[row, col]) else 1.0
            if fdr < 0.05:
                ax.text(j + 0.5, i + 0.85, '*', ha='center', va='center',
                        fontsize=12, fontweight='bold', color='black')

    ax.set_title('IVD Gene Set Enrichment (GSEA NES)\n* = FDR < 0.05')
    plt.tight_layout()
    fig.savefig(PATHWAY_DIR / "gsea_ivd_custom_heatmap.png", dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  GSEA IVD heatmap saved")


# ═══════════════════════════════════════════════════════════════════════════════
# Part 2: TF Activity Inference
# ═══════════════════════════════════════════════════════════════════════════════

def run_tf_activity(de_results):
    """Infer TF activity from DE results using CollecTRI regulon overlap + Fisher's exact test."""
    import decoupler as dc

    print("\n  Loading CollecTRI regulon network...")
    try:
        net = dc.op.collectri(organism='human')
        print(f"  CollecTRI: {len(net)} interactions, {net['source'].nunique()} TFs")
    except Exception as e:
        print(f"  WARNING: Could not load CollecTRI: {e}")
        print("  Falling back to manual TF scoring...")
        return _run_tf_activity_manual(de_results)

    # Build TF -> target dict with direction
    tf_targets = {}
    for _, row in net.iterrows():
        tf = row['source']
        target = row['target']
        weight = row.get('weight', 1)
        if tf not in tf_targets:
            tf_targets[tf] = {}
        tf_targets[tf][target] = weight

    all_tf_results = []
    comparisons = de_results.groupby(['cell_type', 'comparison'])

    for (ct, comp), group in comparisons:
        sig = group[(group['padj'] < PADJ_THRESHOLD) & (group['log2FC'].abs() > LFC_THRESHOLD)]
        if len(sig) < 5:
            continue

        all_genes = set(group['gene'])
        sig_up = set(sig[sig['log2FC'] > 0]['gene'])
        sig_down = set(sig[sig['log2FC'] < 0]['gene'])
        n_total = len(all_genes)

        for tf, targets in tf_targets.items():
            tf_target_genes = set(targets.keys()) & all_genes
            if len(tf_target_genes) < 5:
                continue

            # Compute concordance: how many TF targets are DE in expected direction?
            activated_targets = {g for g, w in targets.items() if g in all_genes and w > 0}
            repressed_targets = {g for g, w in targets.items() if g in all_genes and w < 0}

            # TF "up" score: activated targets UP + repressed targets DOWN
            concordant_up = len(activated_targets & sig_up) + len(repressed_targets & sig_down)
            # TF "down" score: activated targets DOWN + repressed targets UP
            concordant_down = len(activated_targets & sig_down) + len(repressed_targets & sig_up)

            n_targets = len(tf_target_genes)
            n_sig = len(sig)

            # Fisher's exact: are TF targets enriched in DE genes?
            de_and_target = len(tf_target_genes & set(sig['gene']))
            table = [[de_and_target, n_sig - de_and_target],
                     [n_targets - de_and_target, n_total - n_sig - n_targets + de_and_target]]
            if any(v < 0 for row in table for v in row):
                continue
            _, pval = fisher_exact(table, alternative='greater')

            # Activity = net concordance score
            activity = (concordant_up - concordant_down) / max(n_targets, 1)

            all_tf_results.append({
                'cell_type': ct, 'comparison': comp, 'tf': tf,
                'activity': activity, 'pvalue': pval,
                'n_targets': n_targets, 'n_de_targets': de_and_target,
                'concordant_up': concordant_up, 'concordant_down': concordant_down,
            })

    if all_tf_results:
        tf_df = pd.DataFrame(all_tf_results)
        # FDR correction per comparison
        for (ct, comp), g in tf_df.groupby(['cell_type', 'comparison']):
            mask = (tf_df['cell_type'] == ct) & (tf_df['comparison'] == comp)
            pvals = tf_df.loc[mask, 'pvalue'].values
            if len(pvals) > 0:
                fdr = false_discovery_control(pvals, method='bh')
                tf_df.loc[mask, 'padj'] = fdr

        tf_df.to_csv(TF_DIR / "tf_activity_results.tsv", sep='\t', index=False)

        sig_tf = tf_df[tf_df['padj'] < 0.05]
        print(f"\n  TF activity results: {len(tf_df):,} total, {len(sig_tf):,} significant")

        # Top TFs per comparison
        for (ct, comp), g in sig_tf.groupby(['cell_type', 'comparison']):
            top = g.reindex(g['activity'].abs().nlargest(5).index)
            print(f"    {ct} | {comp}: top TFs = {', '.join(top['tf'].values)}")

        return tf_df
    else:
        print("  WARNING: No TF activity results")
        return pd.DataFrame()


def _run_tf_activity_manual(de_results):
    """Manual TF activity scoring without CollecTRI (fallback)."""
    # Just identify known IVD-relevant TFs in DE results
    ivd_tfs = ['NFKB1', 'RELA', 'HIF1A', 'EPAS1', 'SOX9', 'SOX5', 'SOX6',
               'FOXA1', 'FOXA2', 'TP53', 'TBXT', 'NFE2L2', 'STAT3', 'JUN',
               'FOS', 'MYC', 'RUNX2', 'SP7', 'PPARG', 'ATF3', 'ATF4',
               'XBP1', 'DDIT3', 'EGR1']

    sig = de_results[(de_results['padj'] < PADJ_THRESHOLD) &
                     (de_results['log2FC'].abs() > LFC_THRESHOLD)]
    tf_in_de = sig[sig['gene'].isin(ivd_tfs)]

    if not tf_in_de.empty:
        tf_in_de.to_csv(TF_DIR / "tf_activity_results.tsv", sep='\t', index=False)
        print(f"  Found {len(tf_in_de)} IVD-relevant TFs in DE results")
    return tf_in_de


def plot_tf_activity(tf_df):
    """Plot TF activity heatmap."""
    if tf_df.empty or 'activity' not in tf_df.columns:
        return

    sig = tf_df[tf_df['padj'] < 0.05].copy()
    if sig.empty:
        return

    # Top TFs by mean absolute activity
    top_tfs = sig.groupby('tf')['activity'].apply(lambda x: x.abs().mean()).nlargest(30).index

    sig_top = sig[sig['tf'].isin(top_tfs)]
    sig_top['label'] = sig_top['cell_type'] + ' | ' + sig_top['comparison']

    pivot = sig_top.pivot_table(index='tf', columns='label', values='activity', aggfunc='first')
    if pivot.shape[0] < 2 or pivot.shape[1] < 1:
        return

    fig, ax = plt.subplots(figsize=(max(8, pivot.shape[1] * 1.5),
                                     max(6, pivot.shape[0] * 0.4)))
    sns.heatmap(pivot.fillna(0), cmap='RdBu_r', center=0, ax=ax,
                annot=True, fmt='.1f', annot_kws={'fontsize': 7},
                linewidths=0.5, cbar_kws={'label': 'TF Activity (ULM)'})
    ax.set_title('Top TF Activity Changes (padj < 0.05)')
    plt.tight_layout()
    fig.savefig(TF_DIR / "tf_activity_heatmap.png", dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  TF activity heatmap saved")


# ═══════════════════════════════════════════════════════════════════════════════
# Part 3: Pain Gene Analysis
# ═══════════════════════════════════════════════════════════════════════════════

def run_pain_analysis(de_results):
    """Cross-reference DE genes with pain gene sets."""
    all_pain_genes = {}
    for category, genes in PAIN_GENE_SETS.items():
        for gene in genes:
            if gene not in all_pain_genes:
                all_pain_genes[gene] = []
            all_pain_genes[gene].append(category)

    # Find pain genes in DE results
    pain_in_de = []
    for _, row in de_results.iterrows():
        gene = row['gene']
        if gene in all_pain_genes:
            pain_in_de.append({
                'gene': gene,
                'cell_type': row['cell_type'],
                'comparison': row['comparison'],
                'log2FC': row['log2FC'],
                'padj': row['padj'],
                'baseMean': row['baseMean'],
                'direction': 'UP' if row['log2FC'] > 0 else 'DOWN',
                'pain_categories': '; '.join(all_pain_genes[gene]),
                'significant': (abs(row['log2FC']) > LFC_THRESHOLD) and (row['padj'] < PADJ_THRESHOLD),
            })

    pain_df = pd.DataFrame(pain_in_de)
    if pain_df.empty:
        print("  No pain genes found in DE results")
        return pain_df

    pain_df.to_csv(PAIN_DIR / "pain_genes.tsv", sep='\t', index=False)

    sig_pain = pain_df[pain_df['significant']]
    print(f"\n  Pain genes in DE results: {len(pain_df):,} total, {len(sig_pain):,} significant")
    print(f"  Unique significant pain genes: {sig_pain['gene'].nunique()}")

    if not sig_pain.empty:
        print("\n  Significant pain-associated DE genes:")
        for _, row in sig_pain.sort_values('padj').iterrows():
            print(f"    {row['gene']:12s} {row['direction']:4s} log2FC={row['log2FC']:+.2f} "
                  f"padj={row['padj']:.1e} | {row['cell_type']} | {row['comparison']} "
                  f"[{row['pain_categories']}]")

    return pain_df


def plot_pain_genes(pain_df, de_results):
    """Plot pain gene expression dot plot."""
    if pain_df.empty:
        return

    sig_pain = pain_df[pain_df['significant']].copy()
    if sig_pain.empty:
        print("  No significant pain genes to plot")
        return

    # Create a matrix of pain gene LFCs across cell types
    unique_genes = sig_pain['gene'].unique()
    if len(unique_genes) > 30:
        # Take top 30 by significance
        top_genes = sig_pain.nsmallest(30, 'padj')['gene'].unique()
    else:
        top_genes = unique_genes

    # Also show all pain genes (significant or not) from key comparisons
    key_comparisons = ['healthy_vs_degenerated_all', 'healthy_vs_degenerated_severe', 'mild_vs_severe']
    key_data = de_results[
        (de_results['comparison'].isin(key_comparisons)) &
        (de_results['gene'].isin(list(set().union(*PAIN_GENE_SETS.values()))))
    ].copy()

    if key_data.empty:
        return

    # Heatmap of pain gene LFCs
    key_data['label'] = key_data['cell_type'] + '\n' + key_data['comparison']
    pivot = key_data.pivot_table(index='gene', columns='label', values='log2FC', aggfunc='first')
    padj_pivot = key_data.pivot_table(index='gene', columns='label', values='padj', aggfunc='first')

    # Filter to genes present in at least one comparison
    pivot = pivot.dropna(how='all')
    if pivot.shape[0] < 2:
        return

    # Cluster genes by category
    gene_cats = {}
    for cat, genes in PAIN_GENE_SETS.items():
        for g in genes:
            if g in pivot.index:
                gene_cats[g] = cat
    ordered_genes = sorted(pivot.index, key=lambda g: gene_cats.get(g, 'zzz'))
    pivot = pivot.loc[ordered_genes]

    fig, ax = plt.subplots(figsize=(max(8, pivot.shape[1] * 1.5),
                                     max(6, pivot.shape[0] * 0.35)))
    sns.heatmap(pivot.fillna(0), cmap='RdBu_r', center=0, ax=ax,
                annot=True, fmt='.1f', annot_kws={'fontsize': 6},
                linewidths=0.5, cbar_kws={'label': 'log2FC'})

    # Mark significant cells with asterisks
    for i, gene in enumerate(pivot.index):
        for j, label in enumerate(pivot.columns):
            padj = padj_pivot.loc[gene, label] if gene in padj_pivot.index and label in padj_pivot.columns else 1.0
            lfc = pivot.loc[gene, label]
            if pd.notna(padj) and padj < 0.05 and pd.notna(lfc) and abs(lfc) > 0.5:
                ax.text(j + 0.5, i + 0.85, '*', ha='center', va='center',
                        fontsize=10, fontweight='bold', color='black')

    ax.set_title('Pain-Associated Gene Expression Changes\n* = |log2FC| > 0.5 & padj < 0.05')
    plt.tight_layout()
    fig.savefig(PAIN_DIR / "pain_genes_heatmap.png", dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Pain gene heatmap saved")


# ═══════════════════════════════════════════════════════════════════════════════
# Report Generation
# ═══════════════════════════════════════════════════════════════════════════════

def generate_report(enrichment_df, gsea_df, tf_df, pain_df):
    """Generate HTML interpretation report."""
    html = ['<!DOCTYPE html><html><head>',
            '<title>Interpretation Report — Module 07</title>',
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
            '  .warning { color: #e74c3c; }',
            '  .pass { color: #27ae60; }',
            '  .section { margin: 20px 0; padding: 15px; background: #f8f9fa; border-radius: 5px; border-left: 4px solid #3498db; }',
            '</style></head><body>',
            '<h1>Module 07: Biological Interpretation Report</h1>',
            f'<p>Generated: {datetime.now().strftime("%Y-%m-%d %H:%M")}</p>']

    # Stats
    n_enrich = len(enrichment_df) if not enrichment_df.empty else 0
    n_sig_enrich = len(enrichment_df[enrichment_df['Adjusted P-value'] < 0.05]) if n_enrich > 0 else 0
    n_gsea = len(gsea_df) if not gsea_df.empty else 0
    n_sig_gsea = len(gsea_df[gsea_df['FDR q-val'].astype(float) < 0.05]) if n_gsea > 0 else 0
    n_tf = len(tf_df[tf_df['padj'] < 0.05]) if not tf_df.empty and 'padj' in tf_df.columns else 0
    n_pain = len(pain_df[pain_df['significant']]) if not pain_df.empty else 0

    html.append('<div class="stats-grid">')
    html.append(f'<div class="stat-box"><h3>ORA enriched terms (FDR&lt;0.05)</h3><p>{n_sig_enrich:,}</p></div>')
    html.append(f'<div class="stat-box"><h3>GSEA significant terms</h3><p>{n_sig_gsea:,}</p></div>')
    html.append(f'<div class="stat-box"><h3>Significant TFs</h3><p>{n_tf:,}</p></div>')
    html.append(f'<div class="stat-box"><h3>Pain DE genes</h3><p>{n_pain}</p></div>')
    html.append('</div>')

    # Pathway enrichment section
    html.append('<h2>1. Pathway Enrichment (ORA)</h2>')
    if not enrichment_df.empty:
        sig = enrichment_df[enrichment_df['Adjusted P-value'] < 0.05]
        for ct in sig['cell_type'].unique():
            ct_data = sig[sig['cell_type'] == ct]
            ct_clean = ct.replace('/', '_').replace(' ', '_')
            html.append(f'<h3>{ct}</h3>')

            # Show enrichment plots if they exist
            for direction in ['up', 'down']:
                plot_path = PATHWAY_DIR / f"enrichment_{ct_clean}_{direction}.png"
                if plot_path.exists():
                    html.append(f'<img src="pathway_enrichment/enrichment_{ct_clean}_{direction}.png" '
                               f'alt="{ct} {direction}">')

            # Top terms table
            top = ct_data.nsmallest(15, 'Adjusted P-value')
            html.append('<table><tr><th>Term</th><th>Direction</th><th>Library</th>'
                       '<th>Overlap</th><th>Adj P</th><th>Comparison</th></tr>')
            for _, row in top.iterrows():
                html.append(f'<tr><td>{row["Term"][:80]}</td><td>{row["direction"]}</td>'
                           f'<td>{row["library"]}</td><td>{row.get("Overlap", "")}</td>'
                           f'<td>{row["Adjusted P-value"]:.1e}</td><td>{row["comparison"]}</td></tr>')
            html.append('</table>')

    # GSEA section
    html.append('<h2>2. Gene Set Enrichment Analysis (GSEA)</h2>')
    gsea_heatmap = PATHWAY_DIR / "gsea_ivd_custom_heatmap.png"
    if gsea_heatmap.exists():
        html.append('<h3>IVD Custom Gene Sets</h3>')
        html.append('<img src="pathway_enrichment/gsea_ivd_custom_heatmap.png" alt="GSEA IVD heatmap">')

    if not gsea_df.empty:
        sig_gsea = gsea_df[gsea_df['FDR q-val'].astype(float) < 0.05]
        if not sig_gsea.empty:
            html.append('<h3>Top GSEA Terms</h3>')
            html.append('<table><tr><th>Term</th><th>NES</th><th>FDR</th>'
                       '<th>Cell Type</th><th>Comparison</th><th>Library</th></tr>')
            sig_gsea_sorted = sig_gsea.copy()
            sig_gsea_sorted['FDR q-val'] = sig_gsea_sorted['FDR q-val'].astype(float)
            for _, row in sig_gsea_sorted.nsmallest(30, 'FDR q-val').iterrows():
                html.append(f'<tr><td>{str(row["Term"])[:80]}</td><td>{float(row["NES"]):.2f}</td>'
                           f'<td>{float(row["FDR q-val"]):.1e}</td><td>{row["cell_type"]}</td>'
                           f'<td>{row["comparison"]}</td><td>{row["library"]}</td></tr>')
            html.append('</table>')

    # TF activity section
    html.append('<h2>3. Transcription Factor Activity</h2>')
    tf_heatmap = TF_DIR / "tf_activity_heatmap.png"
    if tf_heatmap.exists():
        html.append('<img src="tf_activity/tf_activity_heatmap.png" alt="TF activity heatmap">')

    if not tf_df.empty and 'padj' in tf_df.columns:
        sig_tf = tf_df[tf_df['padj'] < 0.05]
        if not sig_tf.empty:
            html.append('<h3>Top Differential TF Activities</h3>')
            html.append('<table><tr><th>TF</th><th>Activity</th><th>padj</th>'
                       '<th>Cell Type</th><th>Comparison</th></tr>')
            for _, row in sig_tf.reindex(sig_tf['activity'].abs().nlargest(30).index).iterrows():
                direction = 'UP' if row['activity'] > 0 else 'DOWN'
                html.append(f'<tr><td><strong>{row["tf"]}</strong></td>'
                           f'<td>{row["activity"]:.2f} ({direction})</td>'
                           f'<td>{row["padj"]:.1e}</td><td>{row["cell_type"]}</td>'
                           f'<td>{row["comparison"]}</td></tr>')
            html.append('</table>')

    # Pain gene section
    html.append('<h2>4. Pain-Associated Gene Analysis</h2>')
    pain_heatmap = PAIN_DIR / "pain_genes_heatmap.png"
    if pain_heatmap.exists():
        html.append('<img src="pain_genes_heatmap.png" alt="Pain genes heatmap">')

    if not pain_df.empty:
        sig_pain = pain_df[pain_df['significant']]
        if not sig_pain.empty:
            html.append('<h3>Significant Pain-Associated DE Genes</h3>')
            html.append('<table><tr><th>Gene</th><th>Direction</th><th>log2FC</th><th>padj</th>'
                       '<th>Cell Type</th><th>Comparison</th><th>Pain Category</th></tr>')
            for _, row in sig_pain.sort_values('padj').iterrows():
                color = '#e74c3c' if row['direction'] == 'UP' else '#2980b9'
                html.append(f'<tr><td><strong>{row["gene"]}</strong></td>'
                           f'<td style="color:{color}">{row["direction"]}</td>'
                           f'<td>{row["log2FC"]:+.2f}</td><td>{row["padj"]:.1e}</td>'
                           f'<td>{row["cell_type"]}</td><td>{row["comparison"]}</td>'
                           f'<td>{row["pain_categories"]}</td></tr>')
            html.append('</table>')

    # Checkpoint questions
    html.append('<h2>5. Human Checkpoint Questions</h2>')
    html.append('<ol>')
    html.append('<li>Are the enriched pathways consistent with known IVD degeneration biology?</li>')
    html.append('<li>Are there novel pathways or TFs that warrant follow-up?</li>')
    html.append('<li>Do the pain-associated findings suggest specific cell types as primary contributors to discogenic pain?</li>')
    html.append('<li>Are there actionable targets (druggable genes) among the top hits?</li>')
    html.append('<li>Are there findings that contradict established IVD biology?</li>')
    html.append('</ol>')

    html.append('</body></html>')

    report_path = RESULTS_DIR / "interpretation_report.html"
    report_path.write_text('\n'.join(html))
    print(f"  Report saved: {report_path}")


# ═══════════════════════════════════════════════════════════════════════════════
# Validation
# ═══════════════════════════════════════════════════════════════════════════════

def validate():
    """Run automated validation checks."""
    print("\n" + "=" * 60)
    print("Validating Module 07 outputs")
    print("=" * 60)

    checks = []

    # Check pathway enrichment
    enrich_path = PATHWAY_DIR / "all_enrichment_results.tsv"
    if enrich_path.exists():
        enr = pd.read_csv(enrich_path, sep='\t')
        sig = enr[enr['Adjusted P-value'] < 0.05]
        checks.append(('PASS', f'ORA enrichment results: {len(enr):,} total, {len(sig):,} significant'))

        # Check for expected pathways
        if not sig.empty:
            terms_lower = sig['Term'].str.lower()
            for expected in ['extracellular matrix', 'inflammatory', 'collagen', 'immune']:
                found = terms_lower.str.contains(expected, na=False).any()
                if found:
                    checks.append(('PASS', f'Expected pathway found: "{expected}"'))

            # Sanity: no comparison with >500 significant terms
            per_comp = sig.groupby(['cell_type', 'comparison']).size()
            max_terms = per_comp.max()
            if max_terms > 500:
                checks.append(('WARN', f'Comparison with {max_terms} significant terms (may indicate issue)'))
            else:
                checks.append(('PASS', f'Max enriched terms per comparison: {max_terms} (< 500)'))
    else:
        checks.append(('WARN', 'ORA enrichment results not found'))

    # GSEA
    gsea_path = PATHWAY_DIR / "gsea_results.tsv"
    if gsea_path.exists():
        gsea = pd.read_csv(gsea_path, sep='\t')
        checks.append(('PASS', f'GSEA results exist: {len(gsea):,} terms tested'))
    else:
        checks.append(('WARN', 'GSEA results not found'))

    # TF activity
    tf_path = TF_DIR / "tf_activity_results.tsv"
    if tf_path.exists():
        checks.append(('PASS', 'TF activity results exist'))
    else:
        checks.append(('WARN', 'TF activity results not found'))

    # Pain genes
    pain_path = PAIN_DIR / "pain_genes.tsv"
    if pain_path.exists():
        pain = pd.read_csv(pain_path, sep='\t')
        sig_pain = pain[pain['significant'] == True]
        checks.append(('PASS', f'Pain gene analysis: {len(pain):,} total, {sig_pain["gene"].nunique()} unique significant'))
    else:
        checks.append(('WARN', 'Pain gene results not found'))

    # Report
    report_path = RESULTS_DIR / "interpretation_report.html"
    if report_path.exists():
        checks.append(('PASS', 'Interpretation report generated'))
    else:
        checks.append(('FAIL', 'Interpretation report not generated'))

    passed = all(s != 'FAIL' for s, _ in checks)
    for status, msg in checks:
        color = {'PASS': 'pass', 'FAIL': 'warning', 'WARN': 'info'}.get(status, '')
        print(f"  [{status}] {msg}")

    return passed, checks


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 60)
    print("Module 07: Biological Interpretation")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    args = sys.argv[1:]
    pathways_only = '--pathways-only' in args
    pain_only = '--pain-only' in args
    tf_only = '--tf-only' in args
    validate_only = '--validate-only' in args

    if validate_only:
        passed, checks = validate()
        sys.exit(0 if passed else 1)

    # Load DE results
    print("\nLoading DE results...")
    de_primary, de_confounded = load_de_results()

    enrichment_df = pd.DataFrame()
    gsea_df = pd.DataFrame()
    tf_df = pd.DataFrame()
    pain_df = pd.DataFrame()

    # Part 1: Pathway enrichment
    if not pain_only and not tf_only:
        print("\n" + "=" * 60)
        print("Part 1: Pathway Enrichment (ORA)")
        print("=" * 60)
        enrichment_df = run_pathway_enrichment(de_primary)
        plot_enrichment_summary(enrichment_df)

        print("\n" + "=" * 60)
        print("Part 1b: Gene Set Enrichment Analysis (GSEA)")
        print("=" * 60)
        gsea_df = run_gsea_ranked(de_primary)
        plot_gsea_summary(gsea_df)

    # Part 2: TF activity
    if not pathways_only and not pain_only:
        print("\n" + "=" * 60)
        print("Part 2: TF Activity Inference")
        print("=" * 60)
        tf_df = run_tf_activity(de_primary)
        plot_tf_activity(tf_df)

    # Part 3: Pain genes
    if not pathways_only and not tf_only:
        print("\n" + "=" * 60)
        print("Part 3: Pain-Associated Gene Analysis")
        print("=" * 60)
        pain_df = run_pain_analysis(de_primary)
        plot_pain_genes(pain_df, de_primary)

    # Generate report
    print("\n" + "=" * 60)
    print("Generating Report")
    print("=" * 60)
    # Reload results if needed
    if enrichment_df.empty:
        p = PATHWAY_DIR / "all_enrichment_results.tsv"
        if p.exists():
            enrichment_df = pd.read_csv(p, sep='\t')
    if gsea_df.empty:
        p = PATHWAY_DIR / "gsea_results.tsv"
        if p.exists():
            gsea_df = pd.read_csv(p, sep='\t')
    if tf_df.empty:
        p = TF_DIR / "tf_activity_results.tsv"
        if p.exists():
            tf_df = pd.read_csv(p, sep='\t')
    if pain_df.empty:
        p = PAIN_DIR / "pain_genes.tsv"
        if p.exists():
            pain_df = pd.read_csv(p, sep='\t')

    generate_report(enrichment_df, gsea_df, tf_df, pain_df)

    # Validate
    passed, checks = validate()
    print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Overall validation: {'PASSED' if passed else 'NEEDS REVIEW'}")


if __name__ == "__main__":
    main()
