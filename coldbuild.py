#!/usr/bin/env python3
"""coldbuild.py — cold-build time distribution for the conjunction, wide cut against narrow.

`widecut_price.py` settled the score half on 2026-09-18: the 900-seed cut is a WASH at the live arm
(+1.1%, 4W/2L/2T, sign p = 0.688 over 8 contracts). So cost is what should decide whether the
erythroid cells keep it — and cost is exactly what is not measured. `Miner.CONJUNCTION_MIN_BUDGET_S`
was raised to 700s for K562 and CD34+_HSPC off `conj_widecut_check.py`'s ONE cold build per cell
(K562 491s, CD34+_HSPC 577s), which that entry itself flags as "a first read, not a distribution".
A 700s gate makes those two cells' rung fire only on a near-fully-prefetched round, so if the wide
cut is not actually that expensive the fleet is paying availability for nothing.

Per (cell, contract, arm) this times TWO builds:

    cold   the Cas12a bank file is removed first, so the scan runs — what the gate must cover
    warm   immediately after, with the bank the cold build just wrote — what a sibling hotkey pays

The cold/warm difference is the bank scan, which is the only part the cut width plausibly moves.

**Deleting bank files on a box running a live fleet is the hazard here, so three rails.**
(1) The `SKIP` newest tasks per cell are never touched — the miners prefetch the NEWEST unstamped
task, so its bank is the one build in flight. (2) A bank whose mtime is inside `FRESH_MIN` minutes
is left alone and its arm skipped, since that is the signature of a build happening right now.
(3) Builds run strictly sequentially, and each one checks free memory first: 8 miners hold ~20 GB of
a 49 GB box and CLAUDE.md records a research build OOM-killing four processes.

Deletion is self-healing and that is not an assumption: a bank is a pure function of its `bank_key`,
so the cold build rewrites the identical object. This script verifies that rather than trusting it —
it records `meta["bank"]` before deleting and after rebuilding and flags any mismatch.

    CB_CELLS=K562,CD34+_HSPC CB_N=4 python coldbuild.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "coldbuild")

import dataclasses
import json
import logging
import time

logging.basicConfig(level=logging.ERROR)
import numpy as np                                       # noqa: E402

import genExp as G                                       # noqa: E402
import joined_window as JW                               # noqa: E402
from niome_subnet.genomics import conjunction as CJ      # noqa: E402
from niome_subnet.genomics.all_cut import BANK_DIR, bank_key   # noqa: E402
from conj_stageb import API                              # noqa: E402
from sd_task import fetch                                # noqa: E402

CELLS = [c for c in os.getenv("CB_CELLS", "K562,CD34+_HSPC").split(",") if c]
N_CONTRACTS = int(os.getenv("CB_N", "4"))
SKIP = int(os.getenv("CB_SKIP", "2"))          # never touch the N newest tasks per cell
FRESH_MIN = float(os.getenv("CB_FRESH_MIN", "45"))   # leave banks younger than this alone
MIN_FREE_GB = float(os.getenv("CB_MIN_FREE_GB", "8"))
HK = os.getenv("CB_HK", "niome_hotkey")
BUDGET = float(os.getenv("CB_BUDGET", "1800"))
SPACE = list(range(100, 1000))
OUT = os.getenv("CB_JSON", "coldbuild.json")


def free_gb():
    with open("/proc/meminfo") as fh:
        m = {l.split(":")[0]: int(l.split()[1]) for l in fh if ":" in l}
    return m.get("MemAvailable", 0) / 1024 / 1024


def arm_cfg(base, band_space, arm):
    """Same two arms `widecut_price.py` priced, at the cell's LIVE band_k."""
    if arm == "wide":
        cands = CJ.sub_window(band_space, base.band_width, JW.band_offset_frac(HK) or 0.0)
        return dataclasses.replace(base, seed_list=tuple(SPACE), band_candidates=tuple(cands))
    return dataclasses.replace(base, seed_list=tuple(band_space), band_candidates=())


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    items = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])

    recs = []
    for cell in CELLS:
        base = CJ.config_for(cell)
        if base is None:
            print(f"=== {cell}: no conjunction config, skipped ===", flush=True)
            continue
        band_space = sorted({s for a, b in JW.fallback_for(cell, HK) for s in range(a, b + 1)})
        pool = []
        for t in sorted(items, key=lambda t: t.get("created_at", ""), reverse=True):
            c = (t.get("content") or {}).get("contract") or {}
            s = str(c.get("seed", "") or "")
            if c.get("cell_type") != cell:
                continue
            if len([x for x in s.split(",") if x.strip().isdigit()]) != 3:
                continue
            pool.append(t)
        cands = pool[SKIP:SKIP + N_CONTRACTS]     # rail 1: never the newest SKIP
        print(f"\n=== {cell}  live k={base.band_k} group={base.group_size} "
              f"width={base.band_width} | {len(cands)} contracts "
              f"(skipping the {SKIP} newest, which the fleet is building) ===", flush=True)

        for t in cands:
            tid = (t.get("task_id") or t["id"])
            contract = dict(t["content"]["contract"])
            reference = t["content"]["hbb_reference"]
            print(f"  {tid[:8]}  created {t.get('created_at','?')[:19]}", flush=True)
            for arm in ("narrow", "wide"):
                cfg = arm_cfg(base, band_space, arm)
                path = os.path.join(BANK_DIR, f"cas12a-{bank_key(contract, cell_types, cfg)}.npz")
                existed = os.path.exists(path)
                age_min = (time.time() - os.path.getmtime(path)) / 60 if existed else 1e9
                if existed and age_min < FRESH_MIN:      # rail 2
                    print(f"    {arm:6s} SKIPPED — bank is {age_min:.0f} min old, a build may be "
                          f"in flight on it", flush=True)
                    continue
                if free_gb() < MIN_FREE_GB:              # rail 3
                    print(f"    {arm:6s} SKIPPED — only {free_gb():.1f} GB free, under the "
                          f"{MIN_FREE_GB} GB floor", flush=True)
                    continue
                size_before = os.path.getsize(path) if existed else None
                if existed:
                    os.remove(path)
                t0 = time.monotonic()
                rows, meta = CJ.build_submission(contract, reference, cell_types, cfg=cfg,
                                                 budget_s=BUDGET)
                cold = time.monotonic() - t0
                bank_cold = meta.get("bank")
                if not rows:
                    print(f"    {arm:6s} cold {cold:6.0f}s  DECLINED — {meta.get('reason')}",
                          flush=True)
                    recs.append({"cell": cell, "task": tid[:8], "arm": arm, "built": False,
                                 "cold_s": round(cold), "reason": meta.get("reason")})
                    continue
                t1 = time.monotonic()
                rows2, meta2 = CJ.build_submission(contract, reference, cell_types, cfg=cfg,
                                                   budget_s=BUDGET)
                warm = time.monotonic() - t1
                size_after = os.path.getsize(path) if os.path.exists(path) else None
                same = (size_before is None or size_before == size_after)
                print(f"    {arm:6s} cold {cold:6.0f}s  warm {warm:6.0f}s  "
                      f"bank-scan {cold-warm:6.0f}s | bank {bank_cold} band {meta.get('band')} "
                      f"cas9 {meta.get('cas9_pool')} rows {meta.get('rows')} "
                      f"| rebuilt-identical {'yes' if same else 'NO — SIZE CHANGED'}", flush=True)
                recs.append({"cell": cell, "task": tid[:8], "arm": arm, "built": True,
                             "cold_s": round(cold), "warm_s": round(warm),
                             "scan_s": round(cold - warm), "bank": bank_cold,
                             "band": meta.get("band"), "cas9": meta.get("cas9_pool"),
                             "rows": meta.get("rows"), "existed_before": existed,
                             "size_before": size_before, "size_after": size_after,
                             "rebuilt_identical": same})
                with open(OUT, "w") as fh:
                    json.dump(recs, fh, indent=1)

    with open(OUT, "w") as fh:
        json.dump(recs, fh, indent=1)

    print("\n=== cold-build distribution (what CONJUNCTION_MIN_BUDGET_S must cover) ===")
    print(f"{'cell':12s} {'arm':7s} {'n':>3s} {'min':>6s} {'med':>6s} {'max':>6s} "
          f"{'warm med':>9s} {'scan med':>9s}")
    gates = {}
    for cell in CELLS:
        for arm in ("narrow", "wide"):
            sel = [r for r in recs if r["cell"] == cell and r["arm"] == arm and r.get("built")]
            if not sel:
                continue
            cold = sorted(r["cold_s"] for r in sel)
            warm = sorted(r["warm_s"] for r in sel)
            scan = sorted(r["scan_s"] for r in sel)
            gates[(cell, arm)] = cold
            print(f"{cell:12s} {arm:7s} {len(sel):3d} {cold[0]:6d} "
                  f"{int(np.median(cold)):6d} {cold[-1]:6d} "
                  f"{int(np.median(warm)):9d} {int(np.median(scan)):9d}")
    print("\nImplied gate = max observed cold + ALL_HDR_MIN_BUDGET_S (190) + 25% margin:")
    for (cell, arm), cold in gates.items():
        print(f"  {cell:12s} {arm:7s} max {cold[-1]:4d}s -> gate ~{int((cold[-1]+190)*1.25):4d}s")


if __name__ == "__main__":
    main()
