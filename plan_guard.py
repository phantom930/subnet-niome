#!/usr/bin/env python3
"""plan_guard.py — is data/window_plan.json built on an older stamp than the feed now carries?

Exit 0 = the plan is stale and should be rewritten, 1 = current, 2 = cannot tell (skip).

**This is a different question from `seed_refresh_guard.py`'s and the two are not interchangeable.**
That guard asks "has a round stamped since seed_model/seeds.json was written", i.e. is there new
TRAINING data; it gates a ~65s finetune that is fine to defer an hour. This one asks "is the plan
itself built on a stale stamp", and the answer has a hard three-minute deadline:

    task N stamps at  T + 2h21m   (measured: 2:20:21 / 2:21:05 / 2:21:28 / 2:22:03 over four rounds)
    task N+1 is born  T + 2h24m   (measured median cadence 144 min)

so a stamp lands about **three minutes** before the next task is published and built. The hourly
cron cannot meet that: `fe564247` was built at 00:13:03 off the 23:17 plan, while the stamp it
needed had landed at 00:10:40 — 53 minutes after that plan was written. It therefore played the
round-BEFORE-last's classes, which is not the strategy `repeat_last` was measured as.

The comparison is against the plan's own `history_through`, which `window_plan.py` sets to
`rows[-1]["at"]` where `rows` comes from `load_tasks()` — and that helper keeps only tasks whose
contract carries three comma-joined seeds. So `history_through` is exactly "the newest STAMPED task
this plan saw", which is the input every history-based strategy reads. Comparing against
`generated_at` instead would fire every time the plan aged, and against seeds.json would answer the
other guard's question.

Worth being clear about the value: under a uniform generator band position is free, so playing the
round-before-last costs nothing measurable. This buys FIDELITY to the strategy the selector claims
to be running, not payout.
"""
import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PLAN = ROOT / "data" / "window_plan.json"
TASKS_URL = "https://niome-api.genomes.io/api/v3/tasks?limit=500"


def newest_stamped():
    with urllib.request.urlopen(TASKS_URL, timeout=60) as handle:
        doc = json.load(handle)
    items = doc if isinstance(doc, list) else doc.get("data") or doc.get("items") or []
    best = ""
    for task in items:
        contract = (task.get("content") or {}).get("contract") or {}
        seeds = [s for s in str(contract.get("seed") or "").split(",") if s.strip().isdigit()]
        at = task.get("created_at") or ""
        if len(seeds) == 3 and at > best:
            best = at
    return best


def main():
    try:
        live = newest_stamped()
    except Exception as exc:
        print(f"plan-guard: cannot reach the task feed ({exc})")
        return 2
    if not live:
        print("plan-guard: feed carries no stamped task; skipping")
        return 2
    # A missing plan is stale by definition -- that is the cold-start case, and writing one is
    # exactly what should happen.
    if not PLAN.exists():
        print(f"plan-guard: no plan on disk; newest stamp {live[:19]} -> rewrite")
        return 0
    try:
        have = json.loads(PLAN.read_text()).get("history_through") or ""
    except Exception as exc:
        print(f"plan-guard: plan unreadable ({exc}) -> rewrite")
        return 0
    if live > have:
        print(f"plan-guard: {live[:19]} has stamped, plan was built through "
              f"{have[:19] or 'nothing'} -> rewrite")
        return 0
    print(f"plan-guard: plan is current through {have[:19]}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
