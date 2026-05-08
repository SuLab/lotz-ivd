# Garibaldi smoke test for lotz-ivd

Five small scripts that verify the lotz-ivd stack runs on Garibaldi without
committing to the full pipeline. Total runtime: a few minutes once queued.

## What it checks

1. Conda env builds with `scanpy==1.12` + `scvi-tools==1.4.2` (the load-bearing deps for Module 05)
2. Compute nodes can reach the outside internet (`api.anthropic.com`, GitHub, NCBI)
3. The stack imports cleanly inside an `sbatch` job
4. scVI trains end-to-end on a tiny dataset on CPU
5. scVI trains on a Garibaldi GPU, and reports the device name + wall time

## Usage

Copy this folder to your Garibaldi home, then on a login host:

```bash
bash 00_setup_env.sh        # builds env, stages pbmc3k (~30 MB)
bash 01_check_cluster.sh    # surveys partitions, GPUs, modules
sbatch 03_test_cpu.sbatch
sbatch 04_test_gpu.sbatch
squeue -u $USER
```

Read `ivd_test_cpu_<jobid>.out` and `ivd_test_gpu_<jobid>.out` when the jobs
finish. The "Outbound HTTPS" section is the single most important line to read
— if it shows `FAIL` for `api.anthropic.com`, the agent loop in
`run_pipeline.sh` cannot run on a compute node and you'll need to split the
pipeline into per-module sbatch jobs.

## Things you may need to edit

- `00_setup_env.sh`, `03_test_cpu.sbatch`, `04_test_gpu.sbatch`: the conda
  module name. Run `01_check_cluster.sh` first to see what's actually called
  on Garibaldi (`miniconda3`, `anaconda3`, `python/3.11`, etc.).
- `04_test_gpu.sbatch`: `--gres=gpu:1` requests any GPU. If you want to pin a
  specific card type, change to e.g. `--gres=gpu:gtx1080ti:1`.
- `$SCRATCH`: defaults to `$HOME/scratch/ivd-test`. Override by exporting
  `SCRATCH` to wherever Scripps points your scratch space.
