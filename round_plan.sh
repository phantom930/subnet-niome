#!/usr/bin/env bash
# round_plan.sh — once an hour: resolve live predictions, then write the next round's window plan.
#
# Order matters. seed_window_model.py resolves any pending prediction against the task that has
# since appeared and emits a fresh one, so the shadow log keeps accruing evidence; seed_model is
# then refreshed if a round has stamped since its last run, and window_plan.py writes
# data/window_plan.json last so it picks up whichever prediction was just promoted.
#
# Nothing here restarts a miner. The plan is read fresh on every build, so a hotkey picks up a new
# window on its next task with no restart, and a missing or stale plan leaves it on its
# NIOME_HDR_WINDOW pin (TTL 6h, so one missed run changes nothing).
set -u
cd /root/workspace/subnet-niome || exit 1
PY=/root/workspace/subnet-niome/.venv/bin/python
PY_ML=/root/workspace/subnet-niome/.venv-ml/bin/python
echo "===== $(date -Is) ====="
"$PY" seed_window_model.py 2>&1 | tail -12 | sed 's/^/[live ] /'

# seed_model refresh — walk-forward backtest, per-cell finetune, then a fresh next_prediction.json
# that window_plan.py reads through JOINED_SOURCE. Guarded on new data: the cron is hourly and
# rounds stamp every ~2h24m, so refreshing unconditionally would re-finetune on identical data and
# add a same-sample-count entry to state.json's history for nothing.
#
# Niced and thread-capped deliberately. The miners' all-HDR and all-cut builds are CPU-bound across
# 15 cores (see CLAUDE.md, "Memory, not the GPU, is the binding constraint"), and an unrestrained
# 15-thread torch run landing inside a prefetch window would slow a build that has a deadline.
# A decline here is harmless: window_plan.py falls back to the previous prediction file, and
# _seed_model_windows falls back to _rank_freq_windows if that file is unusable too.
if [ ! -x "$PY_ML" ]; then
  echo "[model] .venv-ml missing — skipping refresh (window_plan uses the existing prediction)"
else
  # Capture rather than pipe: a pipeline's status is the LAST command's, so piping the guard
  # straight into sed would always look like success and the refresh would run every hour.
  guard_out=$("$PY" seed_refresh_guard.py 2>&1); guard_rc=$?
  printf '%s\n' "$guard_out" | sed 's/^/[model] /'
  if [ "$guard_rc" -eq 0 ]; then
    OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 \
      nice -n 10 "$PY_ML" -m seed_model.update --end now 2>&1 | tail -16 | sed 's/^/[model] /'
  fi
fi

"$PY" window_plan.py 2>&1 | tail -14 | sed 's/^/[plan ] /'
