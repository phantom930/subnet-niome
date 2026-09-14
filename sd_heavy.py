#!/usr/bin/env python3
"""sd_heavy.py — the best seed-0 submission at each mutation-weight spread.

`allocate` re-solves the heavy/light split per contract, and it moves the right way: over the six
never-stamped rounds the chosen heavy share tracks the contract's weight spread at **r = +0.975**
(197-199 rows at spread 1.94, 205-207 at 2.54-2.62). What that correlation does NOT establish is
whether the two-stage search lands on the PEAK or merely near it -- the coarse pass steps by 4 with
the k-mer price off, and the priced refinement only looks +/- 4 rows around the unpriced argmax.

So this pins the share with `force_heavy` and maps `weighted x fidelity` directly against n_heavy,
one point per row count, and compares the argmax to what the shipped search picks.

**One screen per contract.** `enumerate_candidates` is the ~300s half of a build and does not
depend on the split, so it runs once and every n_heavy reuses it. That turns a 20-point curve from
100 minutes into about 300s + 20 cheap allocations.

The objective is `weighted x fidelity` and nothing else, because on a seed-0 round every miner that
pins the seed reaches consistency exactly 1.000 and the ranking collapses to that product.

    SH_TASKS=a66f01fa,...  SH_LO=185 SH_HI=225 SH_STEP=2  python -u sd_heavy.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "sdheavy")

import dataclasses as _dc   # noqa: E402
import json                 # noqa: E402
import logging              # noqa: E402
import time                 # noqa: E402

logging.basicConfig(level=logging.ERROR)

import genExp as G          # noqa: E402
from niome_subnet.genomics import seed_depend as SD   # noqa: E402
from sd_task import DIST, field, score, task_content  # noqa: E402

LO = int(os.getenv("SH_LO", "185"))
HI = int(os.getenv("SH_HI", "229"))
STEP = int(os.getenv("SH_STEP", "2"))
OUT = os.getenv("SH_OUT", "sd_heavy.json")
TASKS = [t for t in os.getenv("SH_TASKS", ",".join([
    "a66f01fa-9795-4044-a243-15e33cdf9f81",
    "138a47f7-f842-4490-91b6-7d63466e89d0",
    "67bdd18a-9aa1-43b3-9bc7-13c23fadfed1",
    "e32ec7b3-33b9-468e-8894-a47313dabbe4",
    "8f02f1a4-0658-4590-8011-615b9f1f43f3",
    "19018a0a-2cfa-4cc2-ac3e-7a35b17dc5b9"])).split(",") if t]


def rows_of(chosen, contract):
    rows, seen = [], set()
    for i, rec in enumerate(chosen):
        key = (rec["cas_system"], rec["start"], rec["strand"], rec["guide"])
        if key in seen:
            continue
        seen.add(key)
        rows.append({"experiment_id": f"exp-{i:05d}", "guideRNA": rec["guide"],
                     "target_alignment_start": rec["start"],
                     "target_alignment_end": rec["start"] + rec["length"],
                     "strand": rec["strand"], "mutation": rec["mutation"],
                     "cas_system": rec["cas_system"],
                     "cell_type": contract.get("cell_type")})
    return rows


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    out = {"lo": LO, "hi": HI, "step": STEP, "contracts": []}
    print(f"seed-0 heavy-share curve   n_heavy {LO}..{HI} step {STEP}\n")

    for tid in TASKS:
        task, contract, reference = task_content(tid)
        w = contract.get("mutation_weights") or {}
        vals = sorted(w.values(), reverse=True)
        spread = vals[0] / vals[-1] if vals and vals[-1] else float("nan")
        pinned = dict(contract)
        pinned["seed"] = 0
        ctx = G.build_context(pinned, reference, cell_types)
        sites = G.enumerate_sites(ctx, SD.SeedDependConfig().flank, SD.SeedDependConfig().lengths)
        n_rows = pinned["rules"].get("max_experiments") or ctx.max_experiments

        t0 = time.monotonic()
        by_cell = SD.enumerate_candidates(ctx, sites, SD.SeedDependConfig(), None)
        screen_s = time.monotonic() - t0
        if len(by_cell) < 8:
            print(f"  {tid[:8]} only {len(by_cell)}/8 cells — skipped")
            continue

        # What the shipped search picks, for comparison.
        _ch, ship = SD.allocate(by_cell, pinned, ctx, n_rows, SD.SeedDependConfig())
        fld = field(tid)

        print(f"{tid[:8]}  {contract.get('cell_type'):<11} spread {spread:.3f}  "
              f"screen {screen_s:.0f}s  shipped picks heavy {ship.get('heavy')} "
              f"-> {ship.get('product', 0):.2f}")
        curve = []
        for n in range(LO, HI + 1, STEP):
            cfg = _dc.replace(SD.SeedDependConfig(), force_heavy=n)
            ch, m = SD.allocate(by_cell, pinned, ctx, n_rows, cfg)
            if not ch:
                continue
            curve.append({"heavy": n, "weighted": m["weighted"],
                          "fidelity": m["fidelity"], "product": m["product"]})
        if not curve:
            continue
        peak = max(curve, key=lambda x: x["product"])
        gain = peak["product"] - ship.get("product", 0)
        # Score the true peak against the round's real field.
        cfg = _dc.replace(SD.SeedDependConfig(), force_heavy=peak["heavy"])
        ch, _m = SD.allocate(by_cell, pinned, ctx, n_rows, cfg)
        rows = rows_of(ch, pinned)
        s = score(rows, pinned, reference, cell_types, seed=0)
        rank = sum(1 for r in fld if r[0] > s["final"]) + 1
        share = DIST[rank - 1] if rank <= len(DIST) else 0.0
        ship_rank = sum(1 for r in fld if r[0] > ship.get("product", 0)) + 1
        ship_share = DIST[ship_rank - 1] if ship_rank <= len(DIST) else 0.0
        print(f"    true peak heavy {peak['heavy']:>3} -> {peak['product']:.2f} "
              f"({gain:+.2f} vs shipped)   rank {rank} ({share:.1%}) "
              f"vs shipped rank {ship_rank} ({ship_share:.1%})")
        top = sorted(curve, key=lambda x: -x["product"])[:5]
        print("    top 5: " + "  ".join(f"{c['heavy']}:{c['product']:.1f}" for c in top))
        out["contracts"].append({
            "task": tid, "cell": contract.get("cell_type"), "spread": spread,
            "shipped_heavy": ship.get("heavy"), "shipped_product": ship.get("product"),
            "peak_heavy": peak["heavy"], "peak_product": peak["product"], "gain": gain,
            "peak_rank": rank, "peak_share": share,
            "shipped_rank": ship_rank, "shipped_share": ship_share,
            "screen_s": round(screen_s, 1), "curve": curve})
        json.dump(out, open(OUT, "w"), indent=1)
        print()
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
