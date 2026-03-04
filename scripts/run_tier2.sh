#!/bin/bash
# Tier 2 resident cell integration — session-resilient wrapper
# Runs one approach at a time (A→D) for both NP and AF compartments.
# The script has checkpoint logic — completed approaches are auto-skipped on re-run.
# Log output goes to results/integration/tier2_run.log
set -euo pipefail

cd /home/ubuntu/lotz-ivd
source .venv/bin/activate

LOG="results/integration/tier2_run.log"
mkdir -p results/integration

echo "===== Tier 2 integration started: $(date) =====" | tee -a "$LOG"

for APPROACH in A B C D; do
    echo "" | tee -a "$LOG"
    echo "====== Approach $APPROACH: $(date) ======" | tee -a "$LOG"
    python3 scripts/05_integration.py --tier2-only --approach "$APPROACH" 2>&1 | tee -a "$LOG"
    RC=$?
    echo "====== Approach $APPROACH finished (exit=$RC): $(date) ======" | tee -a "$LOG"
done

# Final validation and report generation
echo "" | tee -a "$LOG"
echo "====== Validation: $(date) ======" | tee -a "$LOG"
python3 scripts/05_integration.py --validate-only 2>&1 | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "===== Tier 2 integration completed: $(date) =====" | tee -a "$LOG"
