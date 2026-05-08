#!/bin/bash
#SBATCH --job-name=ivd-test-gpu
#SBATCH --partition=rtxa6000
#SBATCH --gres=gpu:rtxa6000:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=ivd_test_gpu_%j.out

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

source /etc/profile.d/modules.sh
module purge
module load python/3.11.4
module load cuda/12.4

source "$HOME/envs/ivd-test/bin/activate"

nvidia-smi || echo "(nvidia-smi unavailable)"

python 02_smoke_test.py --data "$HOME/scratch/ivd-test/data/pbmc3k.h5ad" --gpu --epochs 5
