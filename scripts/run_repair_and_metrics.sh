#!/usr/bin/env bash
# Chain: repair v4 bridge counts → relaunch metrics comparison.
# Run under nohup so it survives session close. Writes to two dedicated logs
# and drops a sentinel file at the end so completion is easy to spot.

set -uo pipefail

BASE=/home/ubuntu/lotz-ivd
LOG_DIR=$BASE/results/integration/np_experiment
REPAIR_LOG=$LOG_DIR/repair.log
METRICS_LOG=$LOG_DIR/metrics.log
DONE_MARKER=$LOG_DIR/pipeline_done.txt
FAIL_MARKER=$LOG_DIR/pipeline_failed.txt

rm -f "$DONE_MARKER" "$FAIL_MARKER"

cd "$BASE"

echo "===== STAGE 1: repair ($(date -u '+%Y-%m-%d %H:%M:%S UTC')) =====" > "$REPAIR_LOG"
Rscript scripts/05i_repair_v4_bridge_counts.R >> "$REPAIR_LOG" 2>&1
RC=$?
if [ $RC -ne 0 ]; then
  echo "REPAIR FAILED rc=$RC at $(date -u)" | tee -a "$REPAIR_LOG" > "$FAIL_MARKER"
  exit $RC
fi
echo "===== STAGE 1 complete ($(date -u '+%Y-%m-%d %H:%M:%S UTC')) =====" >> "$REPAIR_LOG"

source "$BASE/.venv/bin/activate"

echo "===== STAGE 2: metrics ($(date -u '+%Y-%m-%d %H:%M:%S UTC')) =====" > "$METRICS_LOG"
python3 scripts/05h_np_experiment_metrics.py >> "$METRICS_LOG" 2>&1
RC=$?
if [ $RC -ne 0 ]; then
  echo "METRICS FAILED rc=$RC at $(date -u)" | tee -a "$METRICS_LOG" > "$FAIL_MARKER"
  exit $RC
fi
echo "===== STAGE 2 complete ($(date -u '+%Y-%m-%d %H:%M:%S UTC')) =====" >> "$METRICS_LOG"

echo "Pipeline complete at $(date -u '+%Y-%m-%d %H:%M:%S UTC')" > "$DONE_MARKER"
