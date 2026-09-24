#!/usr/bin/env python3
"""loop_tile8.py — EIGHT hotkeys tiling two specific seed classes at 4x overlap, one band each.

The arm under test:

    joined window / cut   900 (wide) -- HEK293's pre-2026-09-22 span, forced here for every task
    band space            200 seeds: 200-299 + 400-499 (a JOINED, non-contiguous space)
    sub window            100 seeds
    stride                25  ->  8 x 25 = 200, a CIRCULAR tiling, every seed in exactly 4 windows
    loops                 1 per hotkey -- this is the OFFSET layout, not the loop axis
    everything else       HEK293's live CELL_CONFIG (k=8, group 80, light 12, cell-aware)

**Tasks are selected, not sampled**: only HEK293 rounds whose three seeds include at least one in
200-299 AND at least one in 400-499. Nine of 63 qualify. That selection is the whole point -- it
asks what eight concentrated hotkeys catch on precisely the rounds their two classes were built
for, which is the best case for a concentrated layout and an upper bound on it.

**Read the result as conditional.** Chance puts a round in this set 6.6% of the time
(1 - 2*(800/900)^3 + (700/900)^3), so whatever these nine show has to be discounted by that before
it is compared against a layout that covers the whole 900. A hit rate here is P(hit | both classes
drawn), never P(hit).

Cheap by construction: the cut is the same 900 for all eight, so ONE bank per task, and
`choose_band` restricts to its `candidates` by column lookup, so one `hdr_compliance` pass over the
whole 200-seed space serves every sub-window of it.

    LT_JSON=loop_tile8.json python loop_tile8.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "looptile8")

import dataclasses, gc, json, logging, time                            # noqa: E402

logging.basicConfig(level=logging.ERROR)
import numpy as np                                                    # noqa: E402

import genExp as G                                                    # noqa: E402
import joined_window as JW                                            # noqa: E402
from niome_subnet.genomics import conjunction as CJ                   # noqa: E402
from niome_subnet.genomics import mt19937 as MT                       # noqa: E402
from niome_subnet.genomics.validation import _parse_seeds             # noqa: E402
from conj_stageb import API                                           # noqa: E402
from sd_task import fetch                                             # noqa: E402

CELL = os.getenv("LT_CELL", "HEK293")
NHK = int(os.getenv("LT_HK", "8"))
WIDTH = int(os.getenv("LT_WIDTH", "100"))
STRIDE = int(os.getenv("LT_STRIDE", "25"))
OUT = os.getenv("LT_JSON", "loop_tile8.json")
BLOCK_A = (200, 300)
BLOCK_B = (400, 500)
SPACE = list(range(*BLOCK_A)) + list(range(*BLOCK_B))
assert NHK * STRIDE == len(SPACE), f"{NHK}x{STRIDE} != {len(SPACE)}: not an exact circular tiling"
V_CLEAN = {"HEK293": 0.212}
V_DIRTY = {"HEK293": 0.082}


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    items = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])
    A, B = set(range(*BLOCK_A)), set(range(*BLOCK_B))
    tasks = []
    for t in sorted(items, key=lambda x: x.get("created_at", ""), reverse=True):
        c = (t.get("content") or {}).get("contract") or {}
        if c.get("cell_type") != CELL:
            continue
        sd = _parse_seeds(c.get("seed"))
        if len(sd) != 3 or not (any(s in A for s in sd) and any(s in B for s in sd)):
            continue
        tasks.append(t)

    base = CJ.config_for(CELL)
    windows = [CJ.sub_window(SPACE, WIDTH, (i * STRIDE) / float(len(SPACE))) for i in range(NHK)]
    out = json.load(open(OUT)) if os.path.exists(OUT) else []
    done = {r["task"] for r in out}
    print(f"=== {CELL} | {NHK} hotkeys x width {WIDTH} at stride {STRIDE} over "
          f"{BLOCK_A[0]}-{BLOCK_A[1]-1} + {BLOCK_B[0]}-{BLOCK_B[1]-1} | cut 900 wide | "
          f"k={base.band_k} group={base.group_size} | {len(tasks)} qualifying tasks ===", flush=True)
    for i, w in enumerate(windows):
        print(f"    h{i}  {len(w)} seeds  {min(w)}..{max(w)}", flush=True)
    print(flush=True)

    cut = list(JW.FULL_SPACE)
    for ti, t in enumerate(tasks, 1):
        tid = (t.get("task_id") or t["id"])[:8]
        if tid in done:
            continue
        contract, reference = dict(t["content"]["contract"]), t["content"]["hbb_reference"]
        seeds = _parse_seeds(contract["seed"])
        t0 = time.monotonic()
        mf = max(1, round(base.cas12a_max_fail * len(cut) / 900))
        cfg0 = dataclasses.replace(base, seed_list=tuple(cut), start_seed=cut[0],
                                   end_seed=cut[-1], cas12a_max_fail=mf)
        ctx = G.build_context(dict(contract, seed=0), reference, cell_types)
        sites = G.enumerate_sites(ctx, 3000, (20, 23))
        path = os.path.join(CJ.BANK_DIR,
                            f"cas12a-{CJ.bank_key(dict(contract, seed=0), cell_types, cfg0)}.npz")
        if not os.path.exists(path):
            os.makedirs(CJ.BANK_DIR, exist_ok=True)
            bank = CJ.build_bank(dict(contract, seed=0), reference, cell_types, ctx, sites,
                                 cfg0, None)
            MT.free_gpu_memory()
            if not bank:
                print(f"  {tid} bank empty, skipped", flush=True); continue
            CJ.save_bank(path, bank); del bank; gc.collect()
        records = CJ.load_bank(path, limit=cfg0.bank_keep)
        n_bank = len(records)
        ok = CJ.hdr_compliance(records, dict(contract, seed=0), cell_types, ctx, SPACE,
                               cfg0.band_rule)
        MT.free_gpu_memory()
        ok = {j: np.fromiter(v, dtype=np.int32) for j, v in ok.items()}
        del records; gc.collect()
        cell_ok = None
        if cfg0.band_cell_aware:
            cell_ok = CJ.cas9_cell_probe(dict(contract, seed=0), cell_types, ctx, sites, cfg0,
                                         SPACE, cfg0.band_rule)
            MT.free_gpu_memory()

        hks = []
        for i, window in enumerate(windows):
            band, alive = CJ.choose_band(ok, cfg0.band_k, window, cfg0.group_size, cell_ok, {},
                                         cfg0.band_quota)
            b = sorted(int(x) for x in band)
            hit = [s for s in seeds if s in set(b)]
            hks.append({"hk": i, "window": [min(window), max(window)], "n_window": len(window),
                        "in_window": [s for s in seeds if s in set(window)],
                        "band": b, "formed": len(b) == cfg0.band_k, "alive": len(alive),
                        "hits": hit, "builds": []})
        del ok, cell_ok; gc.collect()

        for h in hks:
            if not h["hits"]:
                continue
            bcfg = dataclasses.replace(cfg0, band_candidates=tuple(h["band"]))
            rows, meta = CJ.build_submission(dict(contract, seed=0), reference, cell_types,
                                             cfg=bcfg, budget_s=1800)
            MT.free_gpu_memory()
            if not rows:
                h["builds"].append({"declined": meta.get("reason")}); continue
            got = sorted(int(x) for x in meta["band_seeds"])
            clean = set(int(x) for x in meta["clean_seeds"])
            bh = [s for s in seeds if s in set(got)]
            ch = [s for s in seeds if s in clean and s not in set(got)]
            off = [s for s in seeds if s not in set(got)]
            vc = V_CLEAN.get(CELL, 0.212); vd = V_DIRTY.get(CELL, 0.104)
            h["builds"].append({"band": got, "band_match": got == h["band"],
                                "clean_n": len(clean), "band_hits": bh, "clean_hits": ch,
                                "n_offband": len(off), "n_offband_in_cut": len(off),
                                "cons": (len(bh) + len(ch) * vc + (len(off) - len(ch)) * vd) / 3,
                                "pool": meta.get("pool"), "cas9": meta.get("cas9_pool")})

        allb = sorted({s for h in hks for s in h["band"]})
        caught = sorted({s for h in hks for s in h["hits"]})
        print(f"  [{ti}/{len(tasks)}] {tid} seeds {seeds}  bank {n_bank}  "
              f"band union {len(allb)}/200  caught {len(caught)}/3 {caught}  "
              f"({time.monotonic()-t0:.0f}s)", flush=True)
        out.append({"cell": CELL, "task": tid, "seeds": seeds, "bank": n_bank,
                    "space": [BLOCK_A, BLOCK_B], "width": WIDTH, "stride": STRIDE,
                    "k": cfg0.band_k, "cut": len(cut), "hotkeys": hks,
                    "band_union": len(allb), "caught": caught})
        json.dump(out, open(OUT, "w"), indent=1)
        gc.collect()
    print(f"\nwrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
