#!/usr/bin/env python3
"""seed_refresh_guard.py — has a round stamped since seed_model/seeds.json was written?

Exit 0 = refresh needed, exit 1 = nothing new, exit 2 = cannot tell (treated as "skip" by
round_plan.sh, because re-finetuning on unchanged data is worse than waiting an hour).

round_plan.sh runs hourly but rounds stamp every ~2h24m, so without this guard ~60% of cycles
would re-finetune on identical data: 65s of CPU for nothing, a fresh entry in state.json's
history with the same sample count, and a pointless promotion decision. The model's own
walk-forward number is the thing to watch over time, so keeping that series one-entry-per-new-
sample is what makes it readable.

A task counts as stamped only when its contract carries three comma-joined seeds. The listing
reports `seed: 0` for rounds that were in fact stamped later, so an unstamped row here means
"not yet usable as training data", not "never will be".
"""
import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SEEDS_JSON = ROOT / "seed_model" / "seeds.json"
TASKS_URL = "https://niome-api.genomes.io/api/v3/tasks?limit=500"


def stamped_by_cell():
    with urllib.request.urlopen(TASKS_URL, timeout=60) as handle:
        doc = json.load(handle)
    items = doc if isinstance(doc, list) else doc.get("data") or doc.get("items") or []
    newest = {}
    for task in items:
        contract = (task.get("content") or {}).get("contract") or {}
        seeds = [s for s in str(contract.get("seed") or "").split(",") if s.strip().isdigit()]
        if len(seeds) != 3:
            continue
        cell, at = contract.get("cell_type"), task.get("created_at") or ""
        if cell and at > newest.get(cell, ""):
            newest[cell] = at
    return newest


def have_by_cell():
    doc = json.loads(SEEDS_JSON.read_text())
    rounds = (doc.get("rounds") or {}).get("by_cell_type") or {}
    return {cell: max((r.get("created_at") or "") for r in rows) if rows else ""
            for cell, rows in rounds.items()}


def main():
    try:
        live, have = stamped_by_cell(), have_by_cell()
    except Exception as exc:
        print(f"guard: cannot tell ({exc})")
        return 2
    fresh = [(c, have.get(c, ""), a) for c, a in sorted(live.items()) if a > have.get(c, "")]
    for cell, old, new in fresh:
        print(f"guard: {cell} has {new[:19]} stamped, seeds.json holds {old[:19] or 'nothing'}")
    if not fresh:
        newest = max(live.values(), default="")
        print(f"guard: no new stamped round (newest {newest[:19]}); skipping refresh")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
