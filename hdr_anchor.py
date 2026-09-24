#!/usr/bin/env python3
"""hdr_anchor.py — all-HDR on a JOINED 20-seed band space anchored on the round's own seeds.

`hdr_both.py` found the one recipe that let a single all-HDR hotkey band TWO drawn seeds on
f18ee409: a joined narrow band space plus a reduced group. Density is what does it -- 22% at a
contiguous width 36 caught nothing, 83% at a joined 12 caught one, 92% caught both -- and joined
beats contiguous because two seeds 35 apart force any contiguous window to 36 seeds while two
narrow blocks reach them in 12.

This generalises the recipe and tests it across rounds rather than on one:

    band space   two width-10 BLOCKS, one containing the round's seed in 200-299 and one
                 containing its seed in 400-499. Seed 481 -> block 480-489; seed 208 -> 200-209;
                 joined space = those 20 seeds.
    group        80 (live) and 42, paired within task -- the second lever from `hdr_both.py`
    cut / clean  200-299 + 400-499, enumerated from the shipped rows
    tasks        every erythroid round (K562, HUDEP-2, CD34+_HSPC) with a seed in BOTH ranges

**This is block-level ORACLE knowledge and is not buildable blind.** It needs to know which
width-10 block each seed fell in, which is far finer than the 100-seed windows the seed model
predicts. It is measured to establish what the construction can do when aimed -- the ceiling the
prediction side would have to reach for -- not to propose a live config.

Where a round has several seeds in one range, the LOWEST is the anchor and the others are recorded
as incidental: a second seed in the same block is caught free, one in a different block is out of
reach by construction.
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "hdranchor")

import dataclasses, gc, json, logging, time                           # noqa: E402

logging.basicConfig(level=logging.ERROR)
import numpy as np                                                    # noqa: E402

import genExp as G                                                    # noqa: E402
from niome_subnet.genomics import all_hdr as AH                       # noqa: E402
from niome_subnet.genomics import conjunction as CJ                   # noqa: E402
from niome_subnet.genomics import mt19937 as MT                       # noqa: E402
from niome_subnet.genomics.validation import _parse_seeds             # noqa: E402
from conj_stageb import API                                           # noqa: E402
from sd_task import fetch, score                                      # noqa: E402
from widecut_price import records_of                                  # noqa: E402

CELLS = [c for c in os.getenv("HA_CELLS", "K562,HUDEP-2,CD34+_HSPC").split(",") if c]
GROUPS = [int(x) for x in os.getenv("HA_GROUPS", "80,42").split(",")]
BLOCK = int(os.getenv("HA_BLOCK", "10"))
OUT = os.getenv("HA_JSON", "hdr_anchor.json")
MIN_FREE_GB = float(os.getenv("HA_MIN_FREE_GB", "7"))
RA, RB = (200, 300), (400, 500)
SPACE = sorted(list(range(*RA)) + list(range(*RB)))


def free_gb() -> float:
    m = {}
    for line in open("/proc/meminfo"):
        p = line.split()
        m[p[0].rstrip(":")] = int(p[1])
    return m.get("MemAvailable", 0) / 1048576.0


def wait_for_memory(tag=""):
    for _ in range(120):
        if free_gb() >= MIN_FREE_GB:
            return
        print(f"      waiting for memory ({free_gb():.1f} GB) {tag}", flush=True)
        time.sleep(30)


def compliant(rows, contract, cell_types, ctx, rule):
    ok = CJ.hdr_compliance(records_of(rows), contract, cell_types, ctx, SPACE, rule)
    MT.free_gpu_memory()
    if not ok:
        return []
    keep = set(SPACE)
    for i in sorted(ok):
        keep &= set(int(x) for x in ok[i])
        if not keep:
            break
    return sorted(keep)


def main():
    cell_types = G.fetch_cell_types()
    bad = [c for c in CELLS if (cell_types or {}).get(c, {}).get("accessibility", 1.0) >= 0.99]
    if not cell_types or bad:
        raise SystemExit(f"cell-types unusable (accessibility 1.0 fallback on {bad or 'all'})")
    G.load_sequence()
    items = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])
    tasks = []
    for t in sorted(items, key=lambda x: x.get("created_at", ""), reverse=True):
        c = (t.get("content") or {}).get("contract") or {}
        if c.get("cell_type") not in CELLS:
            continue
        sd = _parse_seeds(c.get("seed"))
        if len(sd) != 3:
            continue
        inA = sorted(s for s in sd if RA[0] <= s < RA[1])
        inB = sorted(s for s in sd if RB[0] <= s < RB[1])
        if inA and inB:
            tasks.append((t, inA, inB))

    out = json.load(open(OUT)) if os.path.exists(OUT) else []
    done = {(r["task"], r["group"]) for r in out}
    print(f"=== all-HDR | joined 2 x width-{BLOCK} blocks anchored on the round's own seeds | "
          f"groups {GROUPS} | cut+clean {len(SPACE)} (200-299 + 400-499) | "
          f"cells {','.join(CELLS)} | {len(tasks)} qualifying tasks ===\n", flush=True)

    for ti, (t, inA, inB) in enumerate(tasks, 1):
        tid = (t.get("task_id") or t["id"])[:8]
        contract, reference = dict(t["content"]["contract"]), t["content"]["hbb_reference"]
        seeds = _parse_seeds(contract["seed"])
        cell = contract["cell_type"]
        tA, tB = inA[0], inB[0]
        blocks = [(tA // BLOCK) * BLOCK, (tB // BLOCK) * BLOCK]
        sl = sorted({s for b in blocks for s in range(b, b + BLOCK)})
        extra = [s for s in seeds if s in set(sl) and s not in (tA, tB)]
        c0 = dict(contract, seed=0)
        ctx = None
        print(f"[{ti}/{len(tasks)}] {tid} {cell:11s} {t.get('created_at','')[:16]} seeds {seeds}"
              f"  anchors {tA},{tB} -> blocks {blocks[0]}-{blocks[0]+BLOCK-1},"
              f"{blocks[1]}-{blocks[1]+BLOCK-1}" + (f"  (+{extra} in space)" if extra else ""),
              flush=True)
        base = AH.config_for(cell)
        for grp in GROUPS:
            if (tid, grp) in done:
                continue
            if ctx is None:
                ctx = G.build_context(c0, reference, cell_types)
            wait_for_memory(f"{tid} g{grp}")
            cfg = dataclasses.replace(
                base, hdr_range=(sl[0], sl[-1]), seed_list=tuple(sl), group_size=grp,
                main_max_fail=AH._scaled_max_fail(cell, base.main_max_fail, len(sl),
                                                  AH._native_span(cell)))
            t0 = time.monotonic()
            rows, meta = AH.build_submission(c0, reference, cell_types, cfg=cfg, budget_s=1800)
            MT.free_gpu_memory()
            if not rows:
                print(f"    g{grp:<3d} DECLINED  {meta.get('reason')}", flush=True)
                out.append({"task": tid, "cell": cell, "group": grp, "seeds": seeds,
                            "anchors": [tA, tB], "space": sl,
                            "declined": meta.get("reason"), "cells": meta.get("cells")})
                continue
            band = compliant(rows, contract, cell_types, ctx, "hdr")
            cut = compliant(rows, contract, cell_types, ctx, "cut")
            bs = set(band); cs = set(cut) - bs
            per = [score(rows, contract, reference, cell_types, seed=s) for s in seeds]
            cons = float(np.mean([p["consistency"] for p in per]))
            wxf = per[0]["weighted"] * per[0]["fidelity"]
            hits = [s for s in seeds if s in bs]
            both = tA in bs and tB in bs
            print(f"    g{grp:<3d} band {len(band):3d} ({len(band)/len(sl):4.0%} of {len(sl)})  "
                  f"clean {len(cs):3d}  {tA}:{'Y' if tA in bs else 'n'} "
                  f"{tB}:{'Y' if tB in bs else 'n'}  {'BOTH' if both else ('one' if hits else '-'):4s}"
                  f"  cons {cons:.4f} final {wxf*cons:6.1f}  cells {meta.get('cells')}/8  "
                  f"({time.monotonic()-t0:.0f}s)", flush=True)
            out.append({"task": tid, "cell": cell, "group": grp, "seeds": seeds,
                        "anchors": [tA, tB], "space": sl, "band_n": len(band), "band": band,
                        "clean_n": len(cs), "hits": hits, "both": both,
                        "has_A": tA in bs, "has_B": tB in bs, "round_cons": cons,
                        "weighted": per[0]["weighted"], "fidelity": per[0]["fidelity"],
                        "round_final": wxf * cons, "cells": meta.get("cells"),
                        "cas9": meta.get("cas9_pool")})
            del rows; gc.collect()
            json.dump(out, open(OUT, "w"), indent=1)
        print(flush=True)
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
