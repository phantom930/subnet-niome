#!/usr/bin/env bash
# round_plan.sh — once an hour: resolve live predictions, then write the next round's window plan.
#
# Order matters. seed_window_model.py resolves any pending prediction against the task that has
# since appeared and emits a fresh one, so the shadow log keeps accruing evidence; window_plan.py
# then writes data/window_plan.json, which the miners read per build.
#
# Nothing here restarts a miner. The plan is read fresh on every build, so a hotkey picks up a new
# window on its next task with no restart, and a missing or stale plan leaves it on its
# NIOME_HDR_WINDOW pin (TTL 6h, so one missed run changes nothing).
set -u
cd /root/workspace/subnet-niome || exit 1
PY=/root/workspace/subnet-niome/.venv/bin/python
echo "===== $(date -Is) ====="
"$PY" seed_window_model.py 2>&1 | tail -12 | sed 's/^/[live ] /'
"$PY" window_plan.py 2>&1 | tail -14 | sed 's/^/[plan ] /'
