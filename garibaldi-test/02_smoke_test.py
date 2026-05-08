#!/usr/bin/env python3
"""Garibaldi smoke test for the lotz-ivd stack.

Checks: package imports, compute-node outbound HTTPS, and a tiny end-to-end
scVI training run on the staged pbmc3k dataset. Run via the CPU/GPU sbatch
scripts in this directory; not meant to be invoked directly on a login node.
"""
import argparse
import os
import socket
import sys
import time
import urllib.error
import urllib.request


def banner(msg: str) -> None:
    print(f"\n=== {msg} ===", flush=True)


def check_egress() -> None:
    banner("Outbound HTTPS from compute node")
    targets = [
        "https://api.anthropic.com",
        "https://github.com",
        "https://ftp.ncbi.nlm.nih.gov",
    ]
    for url in targets:
        try:
            req = urllib.request.Request(url, method="HEAD")
            with urllib.request.urlopen(req, timeout=5) as r:
                print(f"  OK    {url} -> {r.status}")
        except (urllib.error.URLError, TimeoutError, OSError) as e:
            print(f"  FAIL  {url} -> {type(e).__name__}: {e}")


def check_imports():
    banner("Package versions")
    import anndata
    import numpy
    import pandas
    import scanpy as sc
    import scipy
    import scvi
    import torch

    print(f"  python      {sys.version.split()[0]}")
    print(f"  numpy       {numpy.__version__}")
    print(f"  scipy       {scipy.__version__}")
    print(f"  pandas      {pandas.__version__}")
    print(f"  anndata     {anndata.__version__}")
    print(f"  scanpy      {sc.__version__}")
    print(f"  scvi-tools  {scvi.__version__}")
    print(f"  torch       {torch.__version__}")
    print(f"  CUDA avail  {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        for i in range(torch.cuda.device_count()):
            print(f"    GPU {i}: {torch.cuda.get_device_name(i)}")
    return torch


def run_scvi(adata_path: str, use_gpu: bool, max_epochs: int) -> None:
    banner(f"scVI training ({'GPU' if use_gpu else 'CPU'}, {max_epochs} epochs)")
    import scanpy as sc
    import scvi
    import torch

    adata = sc.read_h5ad(adata_path)
    print(f"  loaded {adata.shape[0]} cells x {adata.shape[1]} genes from {adata_path}")

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.highly_variable_genes(
        adata, n_top_genes=1000, flavor="seurat_v3", layer="counts"
    )
    adata = adata[:, adata.var.highly_variable].copy()
    print(f"  after HVG: {adata.shape[0]} cells x {adata.shape[1]} genes")

    scvi.model.SCVI.setup_anndata(adata, layer="counts")
    model = scvi.model.SCVI(adata)

    accelerator = "gpu" if use_gpu and torch.cuda.is_available() else "cpu"
    t0 = time.time()
    model.train(max_epochs=max_epochs, accelerator=accelerator)
    print(f"  trained on {accelerator} in {time.time() - t0:.1f}s")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data", required=True, help="path to staged pbmc3k.h5ad")
    p.add_argument("--gpu", action="store_true")
    p.add_argument("--epochs", type=int, default=5)
    p.add_argument("--skip-egress", action="store_true")
    args = p.parse_args()

    print(f"host:      {socket.gethostname()}")
    print(f"slurm job: {os.environ.get('SLURM_JOB_ID', '(none)')}")
    print(f"partition: {os.environ.get('SLURM_JOB_PARTITION', '(none)')}")

    if not args.skip_egress:
        check_egress()
    check_imports()
    run_scvi(args.data, use_gpu=args.gpu, max_epochs=args.epochs)
    banner("DONE")


if __name__ == "__main__":
    main()
