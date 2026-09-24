#!/usr/bin/env python3
"""clean_dupe_loops.py — does ONE hotkey's cut-clean set change when it re-bands?

`clean_dupe.py` measured duplication ACROSS the 7 hotkeys in a single round: 36.5% of clean-seed
seats are shared, 608 distinct of 900. This is the same question along the other axis — ONE hotkey,
re-banded loop after loop with its OWN previous band seeds excluded, exactly as `hdr_saturate7.py`
drives it. Each loop gets a different band, hence a different surviving pool, hence a different
group, hence a different clean set. Whether those clean sets differ decides what a LOOP buys:

  * if the clean sets barely overlap, extra loops (= extra hotkeys) widen the floor as well as the
    band, and headcount compounds on both halves of the construction;
  * if they largely coincide, loops buy band coverage ONLY and the floor saturates immediately.

The band side is 0% duplicated by construction here -- each loop excludes its predecessors' band
seeds -- so the clean set is the whole question, and any overlap found is a real property of the
min-union rather than an artefact of shared candidates.

Full `build_submission` per loop: the clean set comes from `FastGreedy`'s min-union over the
assembled group and `choose_band` alone cannot produce it. All loops share one bank (one cut space
per contract), so only the first pays the cold scan.

    CDL_N=2 CDL_HK=niome_hotkey CDL_LOOPS=10 python clean_dupe_loops.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "cdloops")

import dataclasses, itertools, json, logging, time                      # noqa: E402
from collections import Counter                                        # noqa: E402

logging.basicConfig(level=logging.ERROR)
import numpy as np                                                     # noqa: E402

import genExp as G                                                     # noqa: E402
import joined_window as JW                                             # noqa: E402
from niome_subnet.genomics import conjunction as CJ                    # noqa: E402
from niome_subnet.genomics.validation import _parse_seeds              # noqa: E402
from conj_stageb import API                                            # noqa: E402
from sd_task import fetch                                              # noqa: E402

CELL = os.getenv("CDL_CELL", "K562")
N = int(os.getenv("CDL_N", "2"))
HK = os.getenv("CDL_HK", "niome_hotkey")
MAX_LOOPS = int(os.getenv("CDL_LOOPS", "10"))
OUT = os.getenv("CDL_JSON", "clean_dupe_loops.json")


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    items = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])
    base = CJ.config_for(CELL)
    window = CJ.sub_window(JW.band_space(HK, None), base.band_width,
                           JW.band_offset_frac(HK) or 0.0)
    cut = JW.conjunction_cut_seeds(HK, CELL, predicted=None, band_candidates=window)

    tasks = []
    for t in sorted(items, key=lambda x: x.get("created_at", ""), reverse=True):
        c = (t.get("content") or {}).get("contract") or {}
        if c.get("cell_type") != CELL or len(_parse_seeds(c.get("seed"))) != 3:
            continue
        tasks.append(t)
        if len(tasks) >= N:
            break

    out = json.load(open(OUT)) if os.path.exists(OUT) else []
    done = {r["task"] for r in out}
    print(f"=== {CELL} {HK} window {min(window)}-{max(window)} ({len(window)} seeds) "
          f"cut {len(cut)} | k={base.band_k} group={base.group_size} | "
          f"{len(tasks)} contracts x <={MAX_LOOPS} loops ===", flush=True)

    for t in tasks:
        tid = (t.get("task_id") or t["id"])
        if tid[:8] in done:
            continue
        contract, reference = dict(t["content"]["contract"]), t["content"]["hbb_reference"]
        print(f"\n  {tid[:8]}", flush=True)
        banned, cleans, bands = set(), [], []
        for loop in range(1, MAX_LOOPS + 1):
            cands = [s for s in window if s not in banned]
            if len(cands) < base.band_k:
                print(f"    loop {loop:>2} candidates exhausted ({len(cands)} left)", flush=True)
                break
            cfg = dataclasses.replace(base, seed_list=tuple(sorted(cut or JW.FULL_SPACE)),
                                      band_candidates=tuple(sorted(cands)))
            t0 = time.monotonic()
            rows, meta = CJ.build_submission(dict(contract, seed=0), reference, cell_types,
                                             cfg=cfg, budget_s=1800)
            if not rows:
                print(f"    loop {loop:>2} DECLINED — {meta.get('reason')}", flush=True)
                break
            b = set(int(x) for x in meta["band_seeds"])
            c = set(int(x) for x in meta["clean_seeds"])
            banned |= b
            bands.append(b)
            cleans.append(c)
            newc = len(c - set().union(*cleans[:-1])) if len(cleans) > 1 else len(c)
            print(f"    loop {loop:>2}  band {len(b):>3}  clean {len(c):>4}  "
                  f"new clean {newc:>4}  cand left {len(cands)-len(b):>3}  "
                  f"{time.monotonic()-t0:.0f}s", flush=True)
        if len(cleans) < 2:
            continue
        seats = sum(len(c) for c in cleans)
        uni = set().union(*cleans)
        pair = [len(a & b) / max(1, min(len(a), len(b)))
                for a, b in itertools.combinations(cleans, 2)]
        adj = [len(cleans[i] & cleans[i+1]) / max(1, min(len(cleans[i]), len(cleans[i+1])))
               for i in range(len(cleans)-1)]
        hist = dict(sorted(Counter(Counter(x for c in cleans for x in c).values()).items()))
        bseats = sum(len(b) for b in bands); buni = set().union(*bands)
        print(f"    CLEAN over {len(cleans)} loops: seats {seats}  distinct {len(uni)}  "
              f"dup {1-len(uni)/seats:.1%}  mean pairwise {np.mean(pair):.1%}  "
              f"adjacent-loop {np.mean(adj):.1%}", flush=True)
        print(f"           seeds held by n loops: {hist}", flush=True)
        print(f"    BAND  seats {bseats} distinct {len(buni)} dup {1-len(buni)/bseats:.1%} "
              f"(0% expected — loops exclude their predecessors)", flush=True)
        out.append({"task": tid[:8], "hk": HK, "loops": len(cleans),
                    "clean_sizes": [len(c) for c in cleans],
                    "seats": seats, "distinct": len(uni), "dup": 1-len(uni)/seats,
                    "pairwise": float(np.mean(pair)), "adjacent": float(np.mean(adj)),
                    "hist": hist, "band_dup": 1-len(buni)/bseats})
        json.dump(out, open(OUT, "w"), indent=1)

    if out:
        print(f"\n=== summary, {len(out)} contracts ===", flush=True)
        print(f"  CLEAN dup across loops {np.mean([r['dup'] for r in out]):.1%}   "
              f"distinct {np.mean([r['distinct'] for r in out]):.0f} of 900   "
              f"pairwise {np.mean([r['pairwise'] for r in out]):.1%}", flush=True)
    print(f"\nwrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
