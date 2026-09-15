#!/usr/bin/env bash
# plan_refresh.sh — rewrite data/window_plan.json as soon as a round STAMPS, not on the hour.
#
# Why this exists, measured rather than assumed. A round's seeds are stamped ~2h21m after the task
# is created, and the next task is created ~2h24m after it -- so a stamp lands about THREE MINUTES
# before the task that should use it is published and built. The hourly cron cannot meet that
# deadline: fe564247 was built at 00:13:03 off the 23:17 plan while the stamp it needed had landed
# at 00:10:40, so it played the round-BEFORE-last's classes. Every history-based strategy in
# strategy_rank.py (repeat_last, hot_hand, cold_hand, marginal) reads the last stamped round, so all
# of them were one round behind roughly always.
#
# What this is worth: FIDELITY, not payout. Band position is free under a uniform generator, so
# playing the round-before-last costs nothing measurable -- it just means the selector was not
# running the strategy it claimed to. Do not expect a score change from this.
#
# Cost is why it can run every minute: the guard is one HTTP GET, and the whole refresh is ~16s
# (walk_update --refresh 11.4s + window_plan.py 3.5s), which fires once per stamp (~every 2h24m).
#
# Ordering: walk_update BEFORE window_plan. A newly stamped round is also a new scoreable sample for
# strategy_rank's `model` arm, and if walk_record.json does not cover it the model is dropped from
# that cell's ranking pool with a note. Dropping it exactly on the cells that just stamped would be
# a systematic exclusion rather than a random one, and 11s of budget is cheap enough to avoid it.
set -u
cd /root/workspace/subnet-niome || exit 1
PY=/root/workspace/subnet-niome/.venv/bin/python
PY_ML=/root/workspace/subnet-niome/.venv-ml/bin/python
LOCK=/root/workspace/subnet-niome/data/.plan.lock

# Non-blocking, and shared with round_plan.sh: both write data/window_plan.json, and the hourly job
# can hold the lock for minutes while it retrains. Skipping is correct there -- that job writes a
# plan from fresher data at the end of its own run, so waiting would only duplicate it.
exec 9>"$LOCK" || exit 0
if ! flock -n 9; then
  echo "$(date -Is) plan-refresh: another plan writer holds the lock; skipping"
  exit 0
fi

# Capture rather than print: this runs 1440 times a day and fires on ~10 of them, so echoing the
# routine "plan is current" line would bury the ~10 real events under a megabyte of noise a month.
# rc 2 (feed unreachable) IS logged -- a guard that cannot see the feed is a silent failure of the
# whole trigger, and looks identical to "nothing stamped" if it is not surfaced.
guard_out=$("$PY" plan_guard.py 2>&1); rc=$?
if [ "$rc" -eq 1 ]; then
  exit 0
fi
printf '%s\n' "$guard_out"
[ "$rc" -ne 0 ] && exit 0

echo "===== $(date -Is) stamp-triggered plan refresh ====="
if [ -x "$PY_ML" ]; then
  OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 \
    nice -n 10 "$PY_ML" -m seed_model.walk_update --refresh 2>&1 | tail -6 | sed 's/^/[walk ] /'
else
  echo "[walk ] .venv-ml missing — skipping; strategy_rank will drop the model arm from the pool"
fi
"$PY" window_plan.py 2>&1 | tail -14 | sed 's/^/[plan ] /'
