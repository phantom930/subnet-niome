#!/usr/bin/env python3
"""hybrid2.py — is the field's 0.92-0.99 consistency a designed mixture, or a leaky pure build?

hybrid_test.py established that a cas-aligned no-cut/HDR mixture scores cons 0.9483 against a pure
pin's 1.0000, at identical weighted — so the mechanism is real but costs 5% for nothing at a single
seed. That leaves two questions it could not answer.

**1. The simpler explanation.** A pure build with a few non-compliant rows produces the same
signature: the forest predicts 245 of 250 rows perfectly and lands just under 1.0. If n=5 strays
give ~0.98, then the miners at 0.9890 are not doing anything clever — they are slightly worse at the
search than a clean pin, and the whole "deliberate mixture" reading is unnecessary. This measures
the curve cons(n) directly.

**2. A mixture the forest CAN learn.** hybrid_test's dist_score arm degenerated: it ranked rows by
weighted, and dist_score is a component of weighted, so every selected row fell on the HDR side of
the threshold. Here the two halves are filled separately (stratified), so the mixture is forced to
exist, and the split rides a feature stage 4 is handed explicitly.

    python hybrid2.py [task_id] [seed]
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")

TASK = _ARGV[1] if len(_ARGV) > 1 else None
SEED = int(_ARGV[2]) if len(_ARGV) > 2 else 150
os.environ["NIOME_INSTANCE"] = "hybrid"

import json  # noqa: E402
import logging  # noqa: E402
import time  # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np  # noqa: E402

import genExp as G  # noqa: E402
from hybrid_test import ROWS, allocate, candidates, rows_of  # noqa: E402
from nocut100 import pick_task  # noqa: E402
from sd_task import score, task_content  # noqa: E402

CACHE = os.getenv("HY_POOLS", f"hybrid_pools_{SEED}.json")
OUT = os.getenv("HY_OUT", "hybrid2.json")
STRAYS = (0, 1, 2, 3, 5, 10, 20, 40)


def pools_for(ctx, sites, contract, cell_types):
    """Compliant guides per (rule, cas) at this seed, cached — the screen costs ~370s."""
    if os.path.exists(CACHE):
        got = json.load(open(CACHE))
        return {tuple(k.split("|")): v for k, v in got.items()}
    out, t0 = {}, time.monotonic()
    for rule in ("hdr", "nocut"):
        for cas in ("Cas9", "Cas12a"):
            out[(rule, cas)] = candidates(ctx, sites, contract, cell_types, SEED, rule, cas)
            print(f"  {rule:<6} {cas:<7} {len(out[(rule, cas)]):>8} "
                  f"({time.monotonic() - t0:.0f}s)")
    json.dump({"|".join(k): v for k, v in out.items()}, open(CACHE, "w"))
    return out


def stratified(hdr_pool, nocut_pool, n_hdr):
    """n_hdr rows from the HDR side and the rest from the no-cut side, each with its cell floors.

    Filling the halves separately is the point: ranking one merged pool by weighted put every row
    on the HDR side of the threshold, because dist_score is part of weighted.
    """
    a, _ = allocate(hdr_pool, n_hdr)
    b, _ = allocate(nocut_pool, ROWS - n_hdr)
    return a + b


def main():
    task, contract, reference = task_content(TASK) if TASK else pick_task()
    cell = contract.get("cell_type")
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    ctx = G.build_context({**contract, "seed": 0}, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    print(f"task {task['id'][:8]}  {cell}  seed {SEED}\n")
    pools = pools_for(ctx, sites, contract, cell_types)
    hdr_all = pools[("hdr", "Cas9")] + pools[("hdr", "Cas12a")]
    nocut_all = pools[("nocut", "Cas9")] + pools[("nocut", "Cas12a")]
    print(f"\npools: hdr {len(hdr_all)}, nocut {len(nocut_all)}\n")

    out = {"strays": {}, "stratified": {}}

    # 1. cons(n) for a pure HDR build with n rows swapped for no-cut-compliant ones. Only the
    #    outcome mix changes; the swapped rows are chosen by weighted like every other row.
    base, _ = allocate(hdr_all, ROWS)
    print("a pure HDR build with n non-compliant (no-cut) rows swapped in:")
    print(f"  {'n strays':>9} {'% of rows':>10} {'weighted':>9} {'cons':>7} "
          f"{'fid':>7} {'final':>8}")
    swap_pool, _ = allocate(nocut_all, max(STRAYS) + 20)
    for n in STRAYS:
        rows = rows_of(base[:ROWS - n] + swap_pool[:n], cell)
        if len(rows) < ROWS:
            continue
        s = score(rows, contract, reference, cell_types, seed=SEED)
        out["strays"][n] = {k: s[k] for k in ("weighted", "consistency", "fidelity", "final")}
        print(f"  {n:>9} {n / ROWS:>9.1%} {s['weighted']:>9.1f} {s['consistency']:>7.4f} "
              f"{s['fidelity']:>7.4f} {s['final']:>8.2f}")

    # 2. A forced mixture: n_hdr rows on HDR, the rest on no-cut, both halves filled properly.
    print("\na forced mixture, both halves filled with their own cell floors:")
    print(f"  {'hdr rows':>9} {'nocut':>6} {'cas9':>5} {'weighted':>9} {'cons':>7} "
          f"{'fid':>7} {'final':>8}")
    for n_hdr in (250, 200, 170, 125, 80, 40, 0):
        chosen = stratified(hdr_all, nocut_all, n_hdr)
        if len(chosen) < ROWS:
            print(f"  {n_hdr:>9} {ROWS - n_hdr:>6}  only {len(chosen)} rows")
            continue
        rows = rows_of(chosen, cell)
        n9 = sum(1 for r in chosen if r["cas_system"] == "Cas9")
        s = score(rows, contract, reference, cell_types, seed=SEED)
        out["stratified"][n_hdr] = {**{k: s[k] for k in
                                       ("weighted", "consistency", "fidelity", "final")},
                                    "cas9_rows": n9}
        print(f"  {n_hdr:>9} {ROWS - n_hdr:>6} {n9:>5} {s['weighted']:>9.1f} "
              f"{s['consistency']:>7.4f} {s['fidelity']:>7.4f} {s['final']:>8.2f}")

    json.dump({"task": task["id"], "cell": cell, "seed": SEED, **out}, open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
