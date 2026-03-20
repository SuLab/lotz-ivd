#!/usr/bin/env python3
"""Module 05 — Workflow comparison: compare integration results across CCA, scANVI, and STACAS.

Loads integration results from each workflow subdirectory, computes integration
metrics (iLISI, batch-ASW, condition-ASW) for each, and generates a side-by-side
comparison table and HTML report.

This script can be run independently after all workflows complete.

Usage:
    python3 scripts/05d_integration_comparison.py
    python3 scripts/05d_integration_comparison.py --object NP
"""

import sys
import warnings
from pathlib import Path
from datetime import datetime
from collections import Counter

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from jinja2 import Template

warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning, module='scanpy')
warnings.filterwarnings('ignore', category=DeprecationWarning)

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent.parent
INT_BASE = BASE / "data" / "integrated"
RESULTS_DIR = BASE / "results" / "integration"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

WORKFLOWS = {
    "cca": {
        "dir": INT_BASE / "cca",
        "ext": ".rds",
        "label": "CCA (Seurat v5)",
        "language": "R",
    },
    "scanvi": {
        "dir": INT_BASE / "scanvi",
        "ext": ".h5ad",
        "label": "scANVI",
        "language": "Python",
    },
    "stacas": {
        "dir": INT_BASE / "stacas",
        "ext": ".rds",
        "label": "STACAS",
        "language": "R",
    },
}

OBJECTS = ["NP", "AF", "CEP", "all_cells"]


# ═════════════════════════════════════════════════════════════════════════════
# METRICS COMPUTATION
# ═════════════════════════════════════════════════════════════════════════════

def _quick_ilisi(embedding, batch_labels, n_sample=5000, k=30):
    """Quick iLISI estimate on a random subsample."""
    from sklearn.neighbors import NearestNeighbors

    n = embedding.shape[0]
    if n > n_sample:
        rng = np.random.RandomState(42)
        idx = rng.choice(n, n_sample, replace=False)
        emb_sub = embedding[idx]
        labels_sub = np.asarray(batch_labels)[idx]
    else:
        emb_sub = embedding
        labels_sub = np.asarray(batch_labels)

    nn = NearestNeighbors(n_neighbors=k, algorithm='auto')
    nn.fit(emb_sub)
    indices = nn.kneighbors(emb_sub, return_distance=False)

    ilisi_scores = []
    for i in range(len(emb_sub)):
        neighbor_labels = labels_sub[indices[i]]
        counts = Counter(neighbor_labels)
        total = sum(counts.values())
        p_sq = sum((c / total) ** 2 for c in counts.values())
        ilisi_scores.append(1.0 / p_sq)

    return np.mean(ilisi_scores)


def compute_metrics_for_adata(adata, batch_key='study', max_cells=30_000):
    """Compute integration quality metrics for an AnnData object.

    Looks for embeddings in obsm: X_pca, X_scanvi_*, or similar.
    Returns dict with metric values.
    """
    metrics = {'n_cells': adata.shape[0]}

    # Find the best embedding to use
    embedding_key = None
    for key in ['X_pca', 'X_scANVI', 'X_scanvi', 'X_integrated']:
        if key in adata.obsm:
            embedding_key = key
            break
    # Fallback: use any scanvi embedding
    if embedding_key is None:
        for key in adata.obsm.keys():
            if 'scanvi' in key.lower() or 'pca' in key.lower():
                embedding_key = key
                break

    if embedding_key is None:
        print(f"    WARNING: No suitable embedding found in obsm keys: {list(adata.obsm.keys())}")
        return metrics

    embedding = adata.obsm[embedding_key]
    metrics['embedding_key'] = embedding_key

    if batch_key not in adata.obs.columns:
        print(f"    WARNING: batch_key '{batch_key}' not in obs columns")
        return metrics

    batch = adata.obs[batch_key].values

    # Subsample if needed
    n_cells = embedding.shape[0]
    if n_cells > max_cells:
        rng = np.random.RandomState(42)
        idx = rng.choice(n_cells, max_cells, replace=False)
        embedding = embedding[idx]
        batch = batch[idx]
        cond = adata.obs['condition_harmonized'].values[idx] if 'condition_harmonized' in adata.obs.columns else None
    else:
        cond = adata.obs['condition_harmonized'].values if 'condition_harmonized' in adata.obs.columns else None

    # iLISI
    try:
        metrics['iLISI'] = float(_quick_ilisi(embedding, batch))
    except Exception as e:
        print(f"    WARNING: iLISI failed: {e}")
        metrics['iLISI'] = np.nan

    # Batch ASW
    try:
        from sklearn.metrics import silhouette_score
        if len(set(batch)) > 1:
            metrics['batch_ASW'] = float(silhouette_score(
                embedding, batch, sample_size=min(5000, len(embedding)),
                random_state=42
            ))
        else:
            metrics['batch_ASW'] = np.nan
    except Exception as e:
        print(f"    WARNING: batch_ASW failed: {e}")
        metrics['batch_ASW'] = np.nan

    # Condition ASW
    try:
        if cond is not None:
            valid = pd.notna(cond)
            if valid.sum() > 100 and len(set(cond[valid])) > 1:
                metrics['condition_ASW'] = float(silhouette_score(
                    embedding[valid], cond[valid],
                    sample_size=min(5000, int(valid.sum())),
                    random_state=42
                ))
            else:
                metrics['condition_ASW'] = np.nan
        else:
            metrics['condition_ASW'] = np.nan
    except Exception as e:
        print(f"    WARNING: condition_ASW failed: {e}")
        metrics['condition_ASW'] = np.nan

    return metrics


def load_rds_metrics(rds_path):
    """Load metrics from an RDS file's associated metrics TSV, or return basic info.

    Since we can't load RDS directly in Python, we check for the per-workflow
    integration_metrics.tsv file.
    """
    # Check for pre-computed metrics TSV in the same directory
    metrics_tsv = rds_path.parent / "integration_metrics.tsv"
    obj_name = rds_path.stem

    if metrics_tsv.exists():
        df = pd.read_csv(metrics_tsv, sep='\t')
        obj_rows = df[df['object'] == obj_name]
        if len(obj_rows) > 0:
            row = obj_rows.iloc[0].to_dict()
            return row

    # If no metrics file, return basic info
    return {
        'n_cells': 'N/A (RDS file — run R metrics script)',
        'iLISI': np.nan,
        'batch_ASW': np.nan,
        'condition_ASW': np.nan,
    }


# ═════════════════════════════════════════════════════════════════════════════
# COMPARISON TABLE
# ═════════════════════════════════════════════════════════════════════════════

def build_comparison_table(object_filter=None):
    """Build a comparison table across all workflows and objects.

    Returns a DataFrame with columns: object, workflow, n_cells, iLISI,
    batch_ASW, condition_ASW, embedding_key, status.
    """
    objects = [object_filter] if object_filter else OBJECTS
    rows = []

    for obj_name in objects:
        for wf_name, wf_info in WORKFLOWS.items():
            wf_dir = wf_info['dir']
            ext = wf_info['ext']
            path = wf_dir / f"{obj_name}{ext}"

            row = {
                'object': obj_name,
                'workflow': wf_name,
                'workflow_label': wf_info['label'],
                'language': wf_info['language'],
                'status': 'missing',
            }

            if not path.exists():
                # Check alternative extensions
                alt_ext = '.h5ad' if ext == '.rds' else '.rds'
                alt_path = wf_dir / f"{obj_name}{alt_ext}"
                if alt_path.exists():
                    path = alt_path
                    ext = alt_ext
                else:
                    rows.append(row)
                    continue

            row['status'] = 'complete'

            if ext == '.h5ad':
                # Load and compute metrics directly
                try:
                    print(f"  Loading {wf_name}/{obj_name}...")
                    adata = sc.read_h5ad(path)
                    metrics = compute_metrics_for_adata(adata)
                    row.update(metrics)
                    del adata
                except Exception as e:
                    print(f"    ERROR loading {path}: {e}")
                    row['status'] = 'error'
            elif ext == '.rds':
                # Load from pre-computed metrics
                metrics = load_rds_metrics(path)
                row.update(metrics)

            rows.append(row)

    df = pd.DataFrame(rows)
    return df


# ═════════════════════════════════════════════════════════════════════════════
# COMPARISON PLOTS
# ═════════════════════════════════════════════════════════════════════════════

def plot_metric_comparison(df):
    """Generate bar plots comparing metrics across workflows."""
    metric_cols = ['iLISI', 'batch_ASW', 'condition_ASW']
    available_metrics = [c for c in metric_cols if c in df.columns and df[c].notna().any()]

    if not available_metrics:
        print("  No metrics available for plotting")
        return

    n_metrics = len(available_metrics)
    fig, axes = plt.subplots(1, n_metrics, figsize=(6 * n_metrics, 5))
    if n_metrics == 1:
        axes = [axes]

    for ax, metric in zip(axes, available_metrics):
        plot_df = df[df[metric].notna()].copy()
        if plot_df.empty:
            ax.set_title(f"{metric} (no data)")
            continue

        sns.barplot(data=plot_df, x='object', y=metric, hue='workflow', ax=ax)
        ax.set_title(metric)
        ax.set_xlabel('Object')
        ax.set_ylabel(metric)
        ax.legend(title='Workflow', fontsize=8)

    plt.tight_layout()
    fig.savefig(RESULTS_DIR / "workflow_comparison_metrics.png",
                dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: workflow_comparison_metrics.png")


# ═════════════════════════════════════════════════════════════════════════════
# COMPARISON REPORT
# ═════════════════════════════════════════════════════════════════════════════

COMPARISON_REPORT_TEMPLATE = """<!DOCTYPE html>
<html>
<head>
<title>Integration Workflow Comparison — Module 05</title>
<style>
  body { font-family: Arial, sans-serif; max-width: 1400px; margin: 0 auto; padding: 20px; }
  h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
  h2 { color: #34495e; margin-top: 30px; }
  table { border-collapse: collapse; width: 100%; margin: 10px 0; }
  th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
  th { background-color: #3498db; color: white; }
  tr:nth-child(even) { background-color: #f2f2f2; }
  .missing { color: #e74c3c; font-style: italic; }
  .complete { color: #27ae60; }
  .best { font-weight: bold; background-color: #d5f5e3; }
  img { max-width: 100%; height: auto; margin: 10px 0; border: 1px solid #ddd; }
  .summary { background: #eaf2f8; padding: 15px; border-radius: 5px; margin: 20px 0; }
</style>
</head>
<body>
<h1>Module 05: Integration Workflow Comparison</h1>
<p>Generated: {{ date }}</p>
<p>Comparing three integration workflows: CCA (label-free), scANVI (semi-supervised), STACAS (semi-supervised).</p>

<div class="summary">
<h3>Summary</h3>
<p>Workflow A (CCA) is the primary workflow. Workflows B (scANVI) and C (STACAS) are run for comparison.
The human checkpoint selects which workflow to carry forward for Modules 06-12.</p>
<ul>
  <li><strong>Workflows available:</strong> {{ n_complete }} / {{ n_total }} (object x workflow combinations)</li>
  <li><strong>Workflows missing:</strong> {{ n_missing }}</li>
</ul>
</div>

<h2>Metrics Comparison</h2>
{{ comparison_table }}

{% if metrics_plot_exists %}
<h2>Metric Plots</h2>
<img src="workflow_comparison_metrics.png" alt="Workflow comparison metrics">
{% endif %}

{% for obj_name in objects %}
<h2>{{ obj_name }} — Side-by-side UMAPs</h2>
<div style="display: flex; gap: 10px; flex-wrap: wrap;">
{% for wf in workflows %}
{% if umap_exists[obj_name][wf] %}
<div>
  <h4>{{ wf }}</h4>
  <img src="umap_{{ wf }}_{{ obj_name }}.png" alt="{{ wf }} {{ obj_name }} UMAP" style="max-width: 400px;">
</div>
{% else %}
<div>
  <h4>{{ wf }}</h4>
  <p class="missing">UMAP not available</p>
</div>
{% endif %}
{% endfor %}
</div>
{% endfor %}

<h2>Selection Criteria</h2>
<ol>
  <li>Workflow A (CCA) is used unless it clearly underperforms on batch mixing or biological conservation.</li>
  <li>If CCA performs adequately, prefer it because it is label-free and avoids circular annotation dependencies.</li>
  <li>If CCA underperforms, the comparison informs which alternative to select.</li>
  <li>If the mesenchymal continuum collapses into a blob in any workflow, reduce integration strength or exclude that workflow.</li>
</ol>

<h2>Human Checkpoint Questions</h2>
<ol>
  <li>Does the integration look reasonable? Are studies mixing well on the UMAP?</li>
  <li>Is there evidence of overcorrection (biological signal lost) in any workflow?</li>
  <li><strong>Which workflow should be carried forward for downstream analysis (Modules 06-12)?</strong></li>
  <li>For the selected workflow: should integration parameters be adjusted for any object/tier?</li>
</ol>

</body>
</html>"""


def generate_comparison_report(df):
    """Generate the HTML comparison report."""
    print("\n  Generating comparison report...")

    # Build comparison table HTML
    display_cols = ['object', 'workflow_label', 'language', 'status', 'n_cells']
    metric_cols = ['iLISI', 'batch_ASW', 'condition_ASW']
    for c in metric_cols:
        if c in df.columns:
            display_cols.append(c)

    display_df = df[[c for c in display_cols if c in df.columns]].copy()

    # Format numeric columns
    for c in metric_cols:
        if c in display_df.columns:
            display_df[c] = display_df[c].apply(
                lambda x: f"{x:.3f}" if pd.notna(x) and isinstance(x, (int, float)) else str(x)
            )

    table_html = display_df.to_html(index=False, classes='comparison-table',
                                     na_rep='N/A', escape=False)

    # Check which UMAPs exist
    umap_exists = {}
    for obj_name in OBJECTS:
        umap_exists[obj_name] = {}
        for wf_name in WORKFLOWS:
            umap_path = RESULTS_DIR / f"umap_{wf_name}_{obj_name}.png"
            umap_exists[obj_name][wf_name] = umap_path.exists()

    n_complete = len(df[df['status'] == 'complete'])
    n_total = len(df)
    n_missing = len(df[df['status'] == 'missing'])

    template = Template(COMPARISON_REPORT_TEMPLATE)
    html = template.render(
        date=datetime.now().strftime("%Y-%m-%d %H:%M"),
        comparison_table=table_html,
        metrics_plot_exists=(RESULTS_DIR / "workflow_comparison_metrics.png").exists(),
        objects=OBJECTS,
        workflows=list(WORKFLOWS.keys()),
        umap_exists=umap_exists,
        n_complete=n_complete,
        n_total=n_total,
        n_missing=n_missing,
    )

    report_path = RESULTS_DIR / "workflow_comparison_report.html"
    report_path.write_text(html)
    print(f"  Report: {report_path}")


# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════

def main():
    args = sys.argv[1:]

    # Parse --object
    object_filter = None
    for i, a in enumerate(args):
        if a == '--object' and i + 1 < len(args):
            object_filter = args[i + 1]

    print("=" * 60)
    print("Module 05 — Integration Workflow Comparison")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    # Build comparison table
    print("\nCollecting metrics from all workflows...")
    df = build_comparison_table(object_filter)

    # Save comparison table
    tsv_path = RESULTS_DIR / "workflow_comparison.tsv"
    df.to_csv(tsv_path, sep='\t', index=False)
    print(f"  Saved: {tsv_path}")

    # Print summary
    print("\n  Comparison summary:")
    for _, row in df.iterrows():
        status_icon = "[OK]" if row['status'] == 'complete' else "[--]"
        print(f"    {status_icon} {row['workflow']:8s} / {row['object']:10s} — "
              f"n_cells={row.get('n_cells', 'N/A')}, "
              f"iLISI={row.get('iLISI', 'N/A')}, "
              f"batch_ASW={row.get('batch_ASW', 'N/A')}")

    # Generate plots
    plot_metric_comparison(df)

    # Generate HTML report
    generate_comparison_report(df)

    print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    main()
