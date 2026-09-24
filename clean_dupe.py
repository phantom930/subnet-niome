#!/usr/bin/env python3
"""clean_dupe.py — how much of the fleet's CUT-CLEAN set is duplicated across the 7 hotkeys?

`hdr_saturate7.py` measured duplication in the BAND (11.5% of band seats land on a seed a sibling
already took). The off-band floor is the other half of the conjunction and it is much larger — 128-146
cut-clean seeds per hotkey against a band of 11 — so whether those overlap decides how much of the
fleet's floor is actually distinct.

**This needs full builds, which is why it was not folded into the saturation sweep.** The clean set
comes from `FastGreedy`'s min-union over the assembled group, so `choose_band` alone cannot produce
it; only `build_submission` does, via `meta["clean_seeds"]`.

All seven hotkeys share one bank here (erythroid is back on the 900-seed wide cut, so one `bank_key`
per contract): the FIRST build pays the cold scan and the rest load it warm. What differs between
them is the band window, hence the band, hence the surviving pool, hence the group — so any overlap
in the clean sets is a property of the construction rather than of shared inputs.

Reported per contract:
  * per-hotkey band and clean-set sizes
  * seats = sum of the clean-set sizes, distinct = |union|, dup rate = 1 - distinct/seats
  * mean pairwise overlap, and the per-seed coverage histogram

    CD_N=4 CD_CELL=K562 python clean_dupe.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "cleandupe")

import dataclasses                                       # noqa: E402
import itertools                                         # noqa: E402
import json                                              # noqa: E402
import logging                                           # noqa: E402
import time                                              # noqa: E402
from collections import Counter                          # noqa: E402

logging.basicConfig(level=logging.ERROR)
import numpy as np                                       # noqa: E402

import genExp as G                                       # noqa: E402
import joined_window as JW                               # noqa: E402
from niome_subnet.genomics import conjunction as CJ      # noqa: E402
from niome_subnet.genomics.validation import _parse_seeds  # noqa: E402
from conj_stageb import API                              # noqa: E402
from sd_task import fetch                                # noqa: E402

CELL = os.getenv("CD_CELL", "K562")
N = int(os.getenv("CD_N", "4"))
HKS = JW.BAND_HK
OUT = os.getenv("CD_JSON", "clean_dupe.json")


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    items = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])
    base = CJ.config_for(CELL)
    tasks = []
    for t in sorted(items, key=lambda x: x.get("created_at", ""), reverse=True):
        c = (t.get("content") or {}).get("contract") or {}
        if c.get("cell_type") != CELL or len(_parse_seeds(c.get("seed"))) != 3:
            continue
        tasks.append(t)
        if len(tasks) >= N:
            break

    out = []
    if os.path.exists(OUT):
        try:
            out = json.load(open(OUT))
        except Exception:
            out = []
    done = {r["task"] for r in out}

    print(f"=== {CELL}  k={base.band_k} group={base.group_size} width={base.band_width} "
          f"| {len(HKS)} hotkeys | {len(tasks)} contracts ===", flush=True)
    for t in tasks:
        tid = (t.get("task_id") or t["id"])
        if tid[:8] in done:
            continue
        contract, reference = dict(t["content"]["contract"]), t["content"]["hbb_reference"]
        cleans, bands, rows_ok = {}, {}, 0
        print(f"\n  {tid[:8]}", flush=True)
        for hk in HKS:
            cands = CJ.sub_window(JW.band_space(hk, None), base.band_width,
                                  JW.band_offset_frac(hk) or 0.0)
            cut = JW.conjunction_cut_seeds(hk, CELL, predicted=None, band_candidates=cands)
            cfg = dataclasses.replace(base, seed_list=tuple(sorted(cut or JW.FULL_SPACE)),
                                      band_candidates=tuple(sorted(cands)))
            t0 = time.monotonic()
            rows, meta = CJ.build_submission(dict(contract, seed=0), reference, cell_types,
                                             cfg=cfg, budget_s=1800)
            if not rows:
                print(f"    {hk[-2:]:>3} DECLINED {meta.get('reason')}", flush=True)
                continue
            rows_ok += 1
            cleans[hk] = set(int(x) for x in meta["clean_seeds"])
            bands[hk] = set(int(x) for x in meta["band_seeds"])
            print(f"    {hk[-2:]:>3} band {len(bands[hk]):>3}  clean {len(cleans[hk]):>4}  "
                  f"pool {meta.get('pool')}  {time.monotonic()-t0:.0f}s", flush=True)
        if rows_ok < 2:
            continue

        def stats(sets):
            seats = sum(len(s) for s in sets.values())
            uni = set().union(*sets.values())
            pair = [len(a & b) / max(1, min(len(a), len(b)))
                    for a, b in itertools.combinations(sets.values(), 2)]
            hist = Counter(Counter(x for s in sets.values() for x in s).values())
            return seats, len(uni), (1 - len(uni) / seats if seats else 0), \
                float(np.mean(pair)) if pair else 0.0, dict(sorted(hist.items()))

        cs, cu, cd, cp, ch = stats(cleans)
        bs, bu, bd, bp, bh = stats(bands)
        print(f"    CLEAN  seats {cs:>4}  distinct {cu:>4}  dup {cd:>6.1%}  "
              f"mean pairwise overlap {cp:.1%}", flush=True)
        print(f"           seeds held by n hotkeys: {ch}", flush=True)
        print(f"    BAND   seats {bs:>4}  distinct {bu:>4}  dup {bd:>6.1%}  "
              f"mean pairwise overlap {bp:.1%}", flush=True)
        out.append({"task": tid[:8], "cell": CELL, "built": rows_ok,
                    "clean": {"seats": cs, "distinct": cu, "dup": cd, "pair": cp, "hist": ch},
                    "band": {"seats": bs, "distinct": bu, "dup": bd, "pair": bp, "hist": bh},
                    "per_hk_clean": {h: len(v) for h, v in cleans.items()}})
        json.dump(out, open(OUT, "w"), indent=1)

    if out:
        print(f"\n=== summary, {len(out)} contracts ===", flush=True)
        print(f"  CLEAN  dup {np.mean([r['clean']['dup'] for r in out]):.1%}   "
              f"distinct {np.mean([r['clean']['distinct'] for r in out]):.0f} of 900   "
              f"seats {np.mean([r['clean']['seats'] for r in out]):.0f}", flush=True)
        print(f"  BAND   dup {np.mean([r['band']['dup'] for r in out]):.1%}   "
              f"distinct {np.mean([r['band']['distinct'] for r in out]):.0f}", flush=True)
    print(f"\nwrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
