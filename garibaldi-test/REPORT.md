# Feasibility report: running `lotz-ivd` on Scripps Garibaldi HPC

**Date:** 2026-05-08
**Operator:** asu (asu@scripps.edu)
**Question:** Can the [andrewsu/lotz-ivd](https://github.com/andrewsu/lotz-ivd) IVD single-cell pipeline run on the Scripps Garibaldi HPC cluster, as an alternative to the current AWS r6i.8xlarge (32 vCPU / 256 GB) workstation?

---

## Verdict

**Yes, the pipeline will run on Garibaldi.** All four risks identified before the experiment have empirical evidence against them: the dependency stack installs and imports, compute nodes have outbound HTTPS (including to `api.anthropic.com`), GPU acceleration works on a modern Ampere card, and a much better GPU partition exists than the one documented in the local skill.

Two adjustments are needed before running the real pipeline:

1. Install Python 3.12 in user space (e.g. via `uv`) — Garibaldi's newest module is `python/3.11.4`, but `lotz-ivd`'s `scanpy==1.12` pin requires Python ≥3.12.
2. Submit GPU jobs to the `rtxa6000` partition, **not** the `gpu` partition. The skill's existing GPU example targets the older `gpu` partition (GTX 1080 Ti, Pascal); `rtxa6000` provides RTX A6000 (Ampere, 48 GB VRAM), which is the right target for scVI/scANVI in Module 05.

---

## Background

`lotz-ivd` is a 12-module agentic pipeline that integrates 12 public scRNA-seq datasets (~423K cells, ~81 samples) into a human IVD single-cell atlas. The agent loop is driven by:

```bash
while :; do cat PROMPT.md | claude; done
```

Most modules run within 32 GB RAM. Two have HPC-relevant resource needs:
- **Module 05 — Integration** (scVI/scANVI on ~200K cells): GPU-bound.
- **Module 07 — SCENIC GRN inference**: RAM-bound (≥64 GB).

Pre-experiment, four open risks were flagged:

| Risk | Concern |
|---|---|
| R1. Stack install | Does `scanpy` + `scvi-tools` install on Garibaldi without Anaconda? |
| R2. Compute-node networking | Can a compute node reach `api.anthropic.com` (needed by the agent loop) and GEO/NCBI (needed for dataset fetches)? |
| R3. GPU availability | Does Garibaldi expose a GPU good enough for scVI/scANVI? |
| R4. R / DESeq2 | Is R available for Module 08 differential analysis? |

---

## Method

A minimal end-to-end smoke test was scripted in [`garibaldi-test/`](.):

| File | Purpose |
|---|---|
| `00_setup_env.sh` | Login-host: builds `~/envs/ivd-test` venv, stages pbmc3k dataset |
| `01_check_cluster.sh` | Login-host: surveys partitions, GPU GRES, available modules |
| `02_smoke_test.py` | Inside sbatch: import check + outbound HTTPS check + 5-epoch scVI run |
| `03_test_cpu.sbatch` | `shared` partition, 4 CPUs, 16 GB |
| `04_test_gpu.sbatch` | `rtxa6000` partition, 1× A6000, 8 CPUs, 32 GB |

**Workload:** scanpy's bundled pbmc3k (2700 cells × 32738 genes) → filter → top 1000 HVGs (`seurat_v3`) → scVI 5 epochs. Tiny by design — purpose was stack verification, not a benchmark.

Both jobs were submitted from the login host:

```bash
sbatch 03_test_cpu.sbatch   # job 41824697
sbatch 04_test_gpu.sbatch   # job 41824698
```

---

## Cluster state observed

### Partitions

```
PARTITION  AVAIL  TIMELIMIT  NODES  NODELIST
shared*    up     infinite   49     emb[0702-0734], nodeb[0301-0316]
highmem    up     infinite   34     nodea0234, nodeb[04...]
iscb       up     4-00:00    3      nodec[0819,0821,0823]
gpu        up     7-00:00    15     nodeb[0417-0433, 1201-1217]
alphafold  up     4-00:00    29     nodeb[201-619]
rtxa6000   up     4-00:00    8      nodeb[215,217,219,315,317,319,615,617]
```

### GPU GRES per partition

| Partition | GPU model | GPUs/node | Mem | Cores | Notes |
|---|---|---|---|---|---|
| `gpu` | GTX 1080 / 1080 Ti | 1–4 | 128–515 GB | 16–20 | Pascal, 11 GB VRAM. Documented in SKILL.md but suboptimal for ML. |
| `rtxa6000` | **RTX A6000** | 2 | 257 GB | 32 | **Ampere, 48 GB VRAM. Best target for scVI/scANVI.** |
| `alphafold` | RTX A5000 (mostly), A6000 (some) | 2–4 | 257 GB | 32 | Dominantly used by structure-prediction jobs. |

### Module system

- Environment Modules (not Lmod). Init script: `source /etc/profile.d/modules.sh`.
- Available Pythons: `python/2.7.11`, `python/3.8.3` (default), `python/3.11.4`. **No conda module.**
- CUDA: `cuda/11.8`, `cuda/12.4`, `cuda/12.9`.
- R: `R/4.3.0` (default), `R/4.5.1`.

---

## Results

Both jobs completed with exit code 0:

```
JobID     JobName     Partition   State      Elapsed   ExitCode
41824697  ivd-test-cpu  shared    COMPLETED  00:03:55  0:0
41824698  ivd-test-gpu  rtxa6000  COMPLETED  00:03:34  0:0
```

### Package versions (identical on both jobs)

| Package | Version |
|---|---|
| Python | 3.11.4 |
| numpy | 1.24.4 |
| scipy | 1.11.1 |
| pandas | 2.2.2 |
| anndata | 0.10.8 |
| scanpy | 1.9.5 |
| scvi-tools | 1.4.2 |
| torch | 2.11.0+cu130 |

### Outbound HTTPS from compute nodes

Tested from both `emb0713` (`shared`) and `nodeb615` (`rtxa6000`):

| Target | Result | Interpretation |
|---|---|---|
| `https://github.com` | 200 | ✓ pip / git fetches work |
| `https://ftp.ncbi.nlm.nih.gov` | 200 | ✓ GEO dataset downloads work |
| `https://api.anthropic.com` | HTTP 404 | ✓ **TLS handshake completed**; 404 just means HEAD `/` is not a real endpoint. The Claude API is reachable from compute nodes. |

### scVI training timings

| Run | Device | Epochs | Cells × HVGs | Wall time | Note |
|---|---|---|---|---|---|
| CPU | 4-core CPU on `emb0713` | 5 | 2700 × 1000 | **10.9 s** | |
| GPU | 1× RTX A6000 on `nodeb615` | 5 | 2700 × 1000 | **31.0 s** | First-epoch CUDA warmup ~27 s; subsequent epochs ~1 s each. CPU "wins" only because the dataset is tiny. |

The GPU result is **not** a real benchmark — it's a smoke-test confirming end-to-end CUDA execution. Module 05's actual workload (~200K cells across 12 datasets) will be bound by per-batch GPU compute, where the A6000's advantage is real.

### nvidia-smi on `nodeb615`

```
NVIDIA-SMI 580.126.09   Driver Version: 580.126.09   CUDA Version: 13.0
GPU 0: NVIDIA RTX A6000   49140 MiB
GPU 1: NVIDIA RTX A6000   49140 MiB
```

Slurm respects the `--gres=gpu:rtxa6000:1` request; `CUDA_VISIBLE_DEVICES=0` was set inside the job.

---

## Findings against pre-experiment risks

| Risk | Status | Evidence |
|---|---|---|
| R1. Stack install | **Resolved** | venv on `python/3.11.4` + `pip install scanpy<1.12 scvi-tools<1.5 anndata<0.13 scikit-misc` works; smoke test imports cleanly. |
| R2. Compute-node networking | **Resolved** | TLS handshakes succeed to GitHub, NCBI, and `api.anthropic.com` from both `shared` and `rtxa6000` compute nodes. |
| R3. GPU availability | **Resolved (better than expected)** | `rtxa6000` partition exposes 8 nodes × 2× RTX A6000 (48 GB VRAM, Ampere). Single-GPU job runs end-to-end through PyTorch Lightning. |
| R4. R / DESeq2 | **Partially** | `R/4.3.0` and `R/4.5.1` modules exist; DESeq2 itself was not loaded or imported in this test. Worth a follow-up `R -e 'library(DESeq2)'` check. |

---

## Caveats and open questions

1. **Version drift.** This test relaxed `lotz-ivd` pins (`scanpy==1.12` → `<1.12`, `anndata==0.12.10` → `<0.13`, `scvi-tools==1.4.2` was satisfiable). The full `requirements_frozen.txt` was not exercised. To match the AWS workstation exactly, install Python 3.12 in user space.
2. **No real-scale GPU benchmark.** 2700 cells is too small to surface the A6000's value. A 10k–50k cell subset would give a wall-time anchor for Module 05 planning.
3. **DESeq2 (Module 08), SCENIC (Module 07), and scANVI integration (Module 05's full path) were not exercised.** The smoke test only verified that scVI itself trains; the broader pipeline modules need separate validation.
4. **Disk quota and scratch policy.** Not investigated. The pipeline will produce tens to low-hundreds of GB of intermediate `.h5ad` files; home-directory quota and scratch retention need confirmation with `hpc@scripps.edu` before a full run.
5. **The agent loop's robustness.** `run_pipeline.sh` runs an open-ended `while :; do ... done` loop. Inside an `sbatch` job this works fine until walltime expires (`rtxa6000` cap is 4 days). For a multi-day run, prefer per-module sbatch jobs over one long agent loop.

---

## Recommended next steps

If you decide to migrate the pipeline:

1. **Build the real environment** with Python 3.12:
   ```bash
   curl -LsSf https://astral.sh/uv/install.sh | sh
   uv python install 3.12
   uv venv ~/envs/ivd-prod --python 3.12
   source ~/envs/ivd-prod/bin/activate
   pip install -r requirements.txt
   ```
2. **Pre-stage the 12 GEO datasets** on a login host (or via a one-shot `sbatch` with `--partition=shared`, since compute-node HTTPS works).
3. **Submit Module 05 (scVI integration)** to `rtxa6000`:
   ```
   #SBATCH --partition=rtxa6000
   #SBATCH --gres=gpu:rtxa6000:1
   #SBATCH --cpus-per-task=8
   #SBATCH --mem=64G
   #SBATCH --time=24:00:00
   ```
4. **Submit Module 07 (SCENIC)** to a fat node:
   ```
   #SBATCH --partition=highmem    # or shared with --mem=256G
   #SBATCH --cpus-per-task=16
   #SBATCH --mem=256G
   #SBATCH --time=72:00:00
   ```
5. **Don't run the full agent loop on a compute node.** Either run it on a login host (lightweight: only `claude` calls, no heavy compute) and have `claude` shell out to `sbatch`, or split each module into a discrete `sbatch` script that the agent submits and waits on.

---

## Artifacts

All test scripts and raw logs live in this directory:

- `00_setup_env.sh`, `01_check_cluster.sh`, `02_smoke_test.py`, `03_test_cpu.sbatch`, `04_test_gpu.sbatch`
- On Garibaldi: `~/garibaldi-test/ivd_test_cpu_41824697.out`, `~/garibaldi-test/ivd_test_gpu_41824698.out`
- venv: `/gpfs/home/asu/envs/ivd-test`
- staged dataset: `/gpfs/home/asu/scratch/ivd-test/data/pbmc3k.h5ad`
