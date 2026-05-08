#!/bin/bash
# Run on a Garibaldi login host. Builds a Python venv that exercises the same
# stack as lotz-ivd Module 05 (scanpy + scvi-tools), and stages a tiny
# bundled dataset for the smoke test. Garibaldi has no conda module, so we
# use python/3.11.4 + venv.
set -euo pipefail

source /etc/profile.d/modules.sh

VENV_DIR="$HOME/envs/ivd-test"
DATA_DIR="$HOME/scratch/ivd-test/data"

mkdir -p "$DATA_DIR" "$(dirname "$VENV_DIR")"

echo "=== Loading python module ==="
module purge
module load python/3.11.4

echo "=== Creating venv: $VENV_DIR ==="
if [[ -d "$VENV_DIR" ]]; then
    echo "venv already exists, skipping create"
else
    python -m venv "$VENV_DIR"
fi
source "$VENV_DIR/bin/activate"
python -m pip install --quiet --upgrade pip

echo "=== Installing test deps ==="
# scanpy 1.12 requires Python >=3.12; Garibaldi has 3.11.4, so use the latest
# 3.11-compatible release. scvi-tools and anndata follow likewise. The full
# lotz-ivd requirements_frozen.txt pins exact versions for reproducibility,
# but for a smoke test we just need the stack to import + train.
pip install --quiet \
    "scanpy<1.12" \
    "scvi-tools<1.5" \
    "anndata<0.13" \
    "scikit-misc" \
    "numpy" "scipy" "pandas"

echo "=== Staging pbmc3k dataset ==="
python - <<EOF
import os, scanpy as sc
out = os.path.join("$DATA_DIR", "pbmc3k.h5ad")
if os.path.exists(out):
    print(f"already staged: {out}")
else:
    adata = sc.datasets.pbmc3k()
    adata.write_h5ad(out)
    print(f"wrote {out}: {adata.shape[0]} cells x {adata.shape[1]} genes")
EOF

cat <<MSG

setup complete
  venv: $VENV_DIR
  data: $DATA_DIR/pbmc3k.h5ad

next:
  sbatch 03_test_cpu.sbatch
  sbatch 04_test_gpu.sbatch
MSG
