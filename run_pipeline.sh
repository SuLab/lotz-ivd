#!/bin/bash
# IVD Analysis Pipeline Runner
#
# Wraps the agent loop with a gate that refuses to start a new iteration
# if analysis_plan.md indicates we are waiting for human review.
#
# Usage: ./run_pipeline.sh

set -euo pipefail
cd "$(dirname "$0")"

while :; do
  # Gate: refuse to run if waiting for human review
  if grep -q "WAITING FOR HUMAN REVIEW" analysis_plan.md; then
    echo "=== PIPELINE HALTED: waiting for human review ==="
    echo ""
    echo "Current active step:"
    grep "Active Step" analysis_plan.md -A 2 | head -3
    echo ""
    echo "Review the notebooks and analysis_plan.md, then update the"
    echo "Active Step to the next module to continue."
    exit 0
  fi

  cat PROMPT.md | claude
done
