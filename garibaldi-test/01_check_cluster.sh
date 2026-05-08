#!/bin/bash
# Survey what's on Garibaldi. Run on a login host. Read-only, no allocations.
source /etc/profile.d/modules.sh 2>/dev/null || true
set +e

echo "=== Partitions ==="
sinfo -o "%P %a %l %D %N"

echo
echo "=== GPU GRES per partition ==="
for p in gpu rtxa6000 alphafold; do
    echo "--- $p ---"
    sinfo -p "$p" -o "%N %G %m %c" 2>/dev/null
done

echo
echo "=== Python / CUDA modules ==="
module avail python cuda 2>&1 | head -20

echo
echo "=== R modules ==="
module avail R 2>&1 | head -10

echo
echo "=== Your queue ==="
squeue -u "$USER"
