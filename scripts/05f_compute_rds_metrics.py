#!/usr/bin/env python3
"""Compute iLISI, batch_ASW, condition_ASW from exported RDS embeddings.

Reads CSV.gz files produced by 05e_export_rds_embeddings.R and computes
the same metrics as the scANVI pipeline, then updates the workflow's
integration_metrics.tsv.

Usage:
    python3 scripts/05f_compute_rds_metrics.py --workflow cca
    python3 scripts/05f_compute_rds_metrics.py --workflow stacas
"""

import sys
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
INT_BASE = BASE / "data" / "integrated"


def quick_ilisi(embedding, batch_labels, n_sample=5000, k=30):
    """Quick iLISI on a random subsample."""
    from sklearn.neighbors import NearestNeighbors

    n = embedding.shape[0]
    if n > n_sample:
        rng = np.random.RandomState(42)
        idx = rng.choice(n, n_sample, replace=False)
        emb = embedding[idx]
        labels = np.asarray(batch_labels)[idx]
    else:
        emb = embedding
        labels = np.asarray(batch_labels)

    nn = NearestNeighbors(n_neighbors=k, algorithm='auto')
    nn.fit(emb)
    indices = nn.kneighbors(emb, return_distance=False)

    scores = []
    for i in range(len(emb)):
        counts = Counter(labels[indices[i]])
        total = sum(counts.values())
        p_sq = sum((c / total) ** 2 for c in counts.values())
        scores.append(1.0 / p_sq)

    return float(np.mean(scores))


def compute_metrics(emb_path):
    """Compute metrics from an exported embedding CSV."""
    print(f"  Loading {emb_path.name}...")
    df = pd.read_csv(emb_path)

    # Separate embedding columns from metadata
    meta_cols = ['study', 'condition_harmonized', 'cell_barcode']
    emb_cols = [c for c in df.columns if c not in meta_cols]
    embedding = df[emb_cols].values
    n_cells = embedding.shape[0]
    print(f"    {n_cells} cells, {len(emb_cols)} dims")

    metrics = {'n_cells': n_cells}

    batch = df['study'].values
    cond = df['condition_harmonized'].values

    # Subsample for metric computation
    max_cells = 30_000
    if n_cells > max_cells:
        rng = np.random.RandomState(42)
        idx = rng.choice(n_cells, max_cells, replace=False)
        embedding_sub = embedding[idx]
        batch_sub = batch[idx]
        cond_sub = cond[idx]
    else:
        embedding_sub = embedding
        batch_sub = batch
        cond_sub = cond

    # iLISI
    try:
        metrics['iLISI'] = quick_ilisi(embedding_sub, batch_sub)
        print(f"    iLISI: {metrics['iLISI']:.4f}")
    except Exception as e:
        print(f"    WARNING: iLISI failed: {e}")
        metrics['iLISI'] = np.nan

    # batch ASW
    try:
        from sklearn.metrics import silhouette_score
        unique_batches = len(set(batch_sub))
        if unique_batches > 1:
            metrics['batch_ASW'] = float(silhouette_score(
                embedding_sub, batch_sub,
                sample_size=min(5000, len(embedding_sub)),
                random_state=42
            ))
            print(f"    batch_ASW: {metrics['batch_ASW']:.4f}")
        else:
            metrics['batch_ASW'] = np.nan
    except Exception as e:
        print(f"    WARNING: batch_ASW failed: {e}")
        metrics['batch_ASW'] = np.nan

    # condition ASW
    try:
        valid = pd.notna(cond_sub)
        if valid.sum() > 100 and len(set(cond_sub[valid])) > 1:
            metrics['condition_ASW'] = float(silhouette_score(
                embedding_sub[valid], cond_sub[valid],
                sample_size=min(5000, int(valid.sum())),
                random_state=42
            ))
            print(f"    condition_ASW: {metrics['condition_ASW']:.4f}")
        else:
            metrics['condition_ASW'] = np.nan
    except Exception as e:
        print(f"    WARNING: condition_ASW failed: {e}")
        metrics['condition_ASW'] = np.nan

    return metrics


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--workflow', required=True, choices=['cca', 'stacas'])
    args = parser.parse_args()

    wf = args.workflow
    export_dir = INT_BASE / wf / "embeddings_export"

    if not export_dir.exists():
        print(f"ERROR: {export_dir} not found. Run 05e_export_rds_embeddings.R first.")
        sys.exit(1)

    # Load existing metrics
    metrics_path = INT_BASE / wf / "integration_metrics.tsv"
    if metrics_path.exists():
        existing = pd.read_csv(metrics_path, sep='\t')
    else:
        existing = pd.DataFrame()

    objects = ["NP", "AF", "CEP", "all_cells"]
    all_metrics = []

    for obj_name in objects:
        emb_path = export_dir / f"{obj_name}_embedding.csv.gz"
        if not emb_path.exists():
            print(f"  SKIP: {emb_path} not found")
            continue

        print(f"\n{'='*60}")
        print(f"Computing metrics for {wf}/{obj_name}")
        print(f"{'='*60}")
        m = compute_metrics(emb_path)
        m['object'] = obj_name
        m['workflow'] = wf
        all_metrics.append(m)

    if not all_metrics:
        print("No metrics computed.")
        return

    new_df = pd.DataFrame(all_metrics)

    # Merge with existing metrics (keep n_clusters_res05 from existing)
    if not existing.empty:
        merged = existing.copy()
        for _, row in new_df.iterrows():
            obj = row['object']
            mask = merged['object'] == obj
            if mask.any():
                for col in ['iLISI', 'batch_ASW', 'condition_ASW', 'n_cells']:
                    merged.loc[mask, col] = row[col]
            else:
                merged = pd.concat([merged, pd.DataFrame([row])], ignore_index=True)
        merged.to_csv(metrics_path, sep='\t', index=False)
    else:
        new_df.to_csv(metrics_path, sep='\t', index=False)

    print(f"\nUpdated: {metrics_path}")
    print(pd.read_csv(metrics_path, sep='\t').to_string())


if __name__ == '__main__':
    main()
