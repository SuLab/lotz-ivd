#!/bin/bash
#SBATCH --job-name=ivd-test-cpu
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=ivd_test_cpu_%j.out

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

source /etc/profile.d/modules.sh
module purge
module load python/3.11.4

source "$HOME/envs/ivd-test/bin/activate"

python 02_smoke_test.py --data "$HOME/scratch/ivd-test/data/pbmc3k.h5ad" --epochs 5
