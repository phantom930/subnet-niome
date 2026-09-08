#!/usr/bin/env python3
"""sd_variants.py — do sibling seed-depend variants submit *similar* rows and *similar* scores?

Four hotkeys on seed-depend are only worth 0.85 of the payout curve (ranks 1-4) if all four still
out-score the field's best. The build is deterministic in (contract, seed, variant), so siblings
must carry distinct variant indices or they submit byte-identical rows from four uids — and the
open question is what that index costs. Two things are measured here, on a real round:

  1. score spread across variants, against that round's actual field
  2. row overlap between variants — the reason for the index in the first place

Both matter in opposite directions: too much overlap and the siblings are the same submission, too
little and the later variants are worse builds.

    python sd_variants.py [task_id] [n_variants]
"""
import os
import sys

_ARGV = list(sys.argv)                 # before bittensor's argparse eats it
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")

# Per-task DATA_DIR, set before settings is imported (it resolves the paths at import time).
# Two of these running at once under one instance name write contract.json, submission.json and
# every stage output to the SAME files: the runs interleave and score one task's rows against
# another's contract. That produced a fidelity of 0.0010 and looked exactly like a build failure.
_TASK_ARG = _ARGV[1] if len(_ARGV) > 1 else "138a47f7-f842-4490-91b6-7d63466e89d0"
os.environ["NIOME_INSTANCE"] = f"sdvar_{_TASK_ARG[:8]}"

import dataclasses  # noqa: E402
import json  # noqa: E402
import logging  # noqa: E402
import time  # noqa: E402

logging.basicConfig(level=logging.ERROR)

import genExp as G  # noqa: E402
from niome_subnet.genomics import seed_depend as SD  # noqa: E402
from sd_task import DIST, field, score, task_content  # noqa: E402

TASK = _TASK_ARG
NVAR = int(_ARGV[2]) if len(_ARGV) > 2 else 4
BUDGET = float(os.getenv("SD_BUDGET", "900"))
OUT = os.getenv("SD_OUT", f"sd_variants_{TASK[:8]}.json")


def rows_for(variant, contract, reference, cell_types):
    """Build (or reuse) one variant's rows. Cached per task+variant — builds are ~100s each."""
    cache = f"sd_task_rows_{TASK[:8]}.json" if variant == 1 else \
            f"sd_var_rows_{TASK[:8]}_v{variant}.json"
    if os.path.exists(cache):
        got = json.load(open(cache))
        return got["rows"], got["meta"], 0.0
    t0 = time.monotonic()
    cfg = dataclasses.replace(SD.SeedDependConfig(), variant=variant)
    rows, meta = SD.build(contract, reference, cell_types, seed=0, cfg=cfg, budget_s=BUDGET)
    if not rows:
        raise SystemExit(f"variant {variant} declined: {meta.get('reason', meta)}")
    json.dump({"rows": rows, "meta": meta}, open(cache, "w"))
    return rows, meta, time.monotonic() - t0


def main():
    task, contract, reference = task_content(TASK)
    seeds = str(contract.get("seed", "")).split(",")
    unstamped = len(seeds) == 1 and seeds[0].strip() in ("0", "")
    print(f"task {TASK}  {contract.get('cell_type')}  opened {task.get('created_at','?')[:19]}")
    print(f"  seed field {contract.get('seed')!r} -> "
          f"{'never stamped' if unstamped else 'stamped'}\n")
    if not unstamped:
        print("  NOTE: this round was stamped, so every seed-0 build scores the ~0.10 floor and the")
        print("  comparison below is not the case seed-depend is for.\n")

    cell_types = G.fetch_cell_types()
    G.load_sequence()

    got = {}
    for v in range(1, NVAR + 1):
        rows, meta, secs = rows_for(v, contract, reference, cell_types)
        s = score(rows, contract, reference, cell_types, seed=0)
        s["guides"] = {r["guideRNA"] for r in rows}
        s["heavy"] = meta.get("heavy")
        got[v] = s
        print(f"  variant {v}: {len(rows)} rows in {secs:>5.0f}s | weighted {s['weighted']:>7.1f} "
              f"x cons {s['consistency']:.4f} x fid {s['fidelity']:.4f} = {s['final']:>7.2f}")

    finals = [got[v]["final"] for v in got]
    print(f"\nscore spread across {len(got)} variants: {min(finals):.2f} to {max(finals):.2f} "
          f"({max(finals) - min(finals):.2f} points)")

    print("\nrow overlap (shared guides, % of 250):")
    print("        " + "".join(f"  v{v:<5}" for v in got))
    for a in got:
        cells = "".join(f"  {len(got[a]['guides'] & got[b]['guides']) / 250:>5.0%}" for b in got)
        print(f"    v{a}  " + cells)

    rows_field = field(TASK)
    tied = [r for r in rows_field if r[1] >= 0.999]
    print(f"\nthe field: {len(rows_field)} miners, {len(tied)} at consistency 1.000, "
          f"best {rows_field[0][0]:.2f}")

    # Place all our siblings into the field at once: they compete with each other for ranks.
    combined = sorted([(f, "OURS") for f in finals]
                      + [(r[0], r[4][:12]) for r in rows_field], reverse=True)
    share = 0.0
    print(f"\n  {'rank':>4} {'final':>8}  who")
    for i, (f, who) in enumerate(combined[:12], 1):
        mark = "  <<" if who == "OURS" else ""
        if who == "OURS" and i <= len(DIST):
            share += DIST[i - 1]
        print(f"  {i:>4} {f:>8.2f}  {who}{mark}")
    solo = DIST[sum(1 for r in rows_field if r[0] > max(finals))] \
        if sum(1 for r in rows_field if r[0] > max(finals)) < len(DIST) else 0.0
    print(f"\n  {len(got)} siblings take {share:.1%} of this round's curve "
          f"(one hotkey alone would take {solo:.1%})")

    json.dump({"task": TASK, "cell": contract.get("cell_type"), "unstamped": unstamped,
               "variants": {v: {k: got[v][k] for k in
                                ("weighted", "consistency", "fidelity", "final", "heavy")}
                            for v in got},
               "overlap": {f"v{a}-v{b}": len(got[a]["guides"] & got[b]["guides"]) / 250
                           for a in got for b in got if a < b},
               "field_best": rows_field[0][0], "n_tied": len(tied),
               "share_fleet": share, "share_solo": solo}, open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
