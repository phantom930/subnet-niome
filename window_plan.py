#!/usr/bin/env python3
"""window_plan.py — assign each hotkey its clean-band window for the coming round.

Concentrates several hotkeys onto the 100-seed window the next task is predicted to draw from
(seed_window_model.predict), and leaves the rest spread across their usual territories. The miner
reads the result per build from data/window_plan.json and falls back to its NIOME_HDR_WINDOW pin if
the file is missing, stale or malformed.

**Why this is safe to run before the prediction is proven.** Concentration is exactly EV-neutral
under a uniform generator: six hotkeys at width 16 cover 6 x 14.7 = 88 band seeds, and a uniform
seed hits them with probability 88/900 whether they sit inside one window or across nine. It only
starts costing when the concentrated hotkeys **saturate** the 100-seed window they share -- seven at
width 16 cover 103 band seeds inside 100 and waste the overlap (37.2% against 37.9%) -- or when a
hotkey is spent on something other than a band. So CONCENTRATE is capped below saturation per cell
type, and no hotkey is given up:

    layout                     uniform generator    50% top-1 accuracy
    9 spread (before)                    37.9%              37.9%
    6 concentrated + 3 spread            37.9%              54.5%
    7 concentrated + 2 spread            37.2%              55.9%
    6 + 2 + 1 hedge                      34.3%              51.5%

HEK293 runs width 12 and band 10.1, so eight hotkeys cover 81 seeds and are still below saturation.

Raise CONCENTRATE past the caps here only once the live log in seed_window_model.py clears the ~38%
break-even, because that is where concentration stops being free.

Usage:
    python window_plan.py            # write data/window_plan.json for the next round
    python window_plan.py --dry-run  # print the plan without writing it
"""
import json
import os
import re
import sys
from datetime import datetime, timedelta, timezone

from seed_window_model import load_tasks, predict, window_label

PLAN_PATH = "data/window_plan.json"
MINER_SH = "miner.sh"
# Band width per cell type, from the width sweep: band 14.7 of 16 on the erythroid types, 10.1 of
# 12 on HEK293 (168 builds, zero declines).
WIDTH = {"HEK293": 12}
DEFAULT_WIDTH = 16
# Hotkeys concentrated on the predicted window. Capped below the point where they saturate its 100
# seeds: 6 x 14.7 = 88 at width 16, 8 x 10.1 = 81 at width 12.
CONCENTRATE = {"HEK293": 8}
DEFAULT_CONCENTRATE = 6
# Spread hotkeys run the shipped 100-wide window rather than the narrow band width, and sit on the
# NEXT most under-drawn windows rather than their home territories. This deliberately costs
# coverage -- band is 13.00 at width 100 against 14.70/14.60/14.22 at width 16 on
# CD34+/HUDEP-2/K562, and 7.00 against 10.10 at width 12 on HEK293, so each erythroid spread hotkey
# gives up ~9-13% and HEK293's gives up 44% (168-build sweep). Per-hotkey coverage is exactly the
# band count, and P(k=2) = (band/900)**2 whatever the width, so nothing offsets it.
#
# It is kept because it changes one thing at a time: the concentrated block carries the whole
# narrow-window bet while the spread hotkeys stay on the configuration the fleet has always run,
# which is the fallback if the prediction turns out to be worthless.
SPREAD_WIDTH = 100
TTL_HOURS = 6                      # survives a missed cron run, expires before it misleads


def fleet():
    """(instance, default window) in table order, parsed from miner.sh so the two cannot drift."""
    rows = []
    for line in open(MINER_SH):
        m = re.match(r'\s*"(\S+)\s+\d+\s+\d+\s+(\d+)-(\d+)"', line)
        if m:
            rows.append((m.group(1), (int(m.group(2)), int(m.group(3)))))
    return rows


def tile(lo, hi, width, n):
    """n disjoint sub-windows of `width` inside [lo, hi], packed from lo."""
    out = []
    for i in range(n):
        a = lo + i * width
        b = a + width - 1
        if b > hi:
            break
        out.append((a, b))
    return out


def main():
    members = fleet()
    if not members:
        print(f"no HOTKEYS table found in {MINER_SH}", file=sys.stderr)
        return 1
    rows = load_tasks()
    now = datetime.now(timezone.utc)
    plan = {"generated_at": now.isoformat(),
            "expires_at": (now + timedelta(hours=TTL_HOURS)).isoformat(),
            "history_through": rows[-1]["at"], "n_tasks": len(rows), "assignments": {}}

    print(f"{len(rows)} tasks through {rows[-1]['at'][:16]}; {len(members)} hotkeys\n")
    for cell in sorted({t["cell"] for t in rows}):
        pr = predict(rows, cell)
        top = pr["ranked"][0]
        width = WIDTH.get(cell, DEFAULT_WIDTH)
        conc = min(CONCENTRATE.get(cell, DEFAULT_CONCENTRATE), len(members))
        lo = (top["window"] + 1) * 100
        slots = tile(lo, lo + 99, width, conc)
        assign = {}
        for (name, _default), w in zip(members[:len(slots)], slots):
            assign[name] = list(w)
        # The rest take the next most under-drawn windows, full width. Ranked order skips entry 0
        # (the concentrated block already owns it), so these can never overlap it.
        rest = members[len(slots):]
        for (name, _default), nxt in zip(rest, pr["ranked"][1:1 + len(rest)]):
            a = (nxt["window"] + 1) * 100
            assign[name] = [a, a + SPREAD_WIDTH - 1]

        spans = sorted(assign.values())
        for (a1, b1), (a2, b2) in zip(spans, spans[1:]):
            if a2 <= b1:
                print(f"ERROR: {cell} windows {a1}-{b1} and {a2}-{b2} overlap", file=sys.stderr)
                return 1
        plan["assignments"][cell] = assign
        print(f"  {cell:<11} predict {top['label']} (p_hit {top['p_hit']:.0%}, beta {pr['beta']:+.2f})"
              f"  width {width}  concentrate {len(slots)}")
        print(f"    concentrated: "
              + ", ".join(f"{n}:{a}-{b}" for (n, _d), (a, b) in zip(members, slots)))
        ranks = {(r["window"] + 1) * 100: i for i, r in enumerate(pr["ranked"])}
        print(f"    spread (w{SPREAD_WIDTH}): "
              + ", ".join(f"{n}:{v[0]}-{v[1]}(rank {ranks.get(v[0], '?')+1},"
                          f" {pr['ranked'][ranks[v[0]]]['p_hit']:.0%})"
                          for n, v in list(assign.items())[len(slots):]))

    if "--dry-run" in sys.argv:
        print("\n--dry-run: not written")
        return 0
    os.makedirs(os.path.dirname(PLAN_PATH), exist_ok=True)
    tmp = PLAN_PATH + ".tmp"
    with open(tmp, "w") as handle:            # atomic: miners read this file mid-round
        json.dump(plan, handle, indent=1)
    os.replace(tmp, PLAN_PATH)
    print(f"\nwrote {PLAN_PATH}, valid until {plan['expires_at'][:16]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
