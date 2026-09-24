#!/usr/bin/env python3
"""hdr_vs_conj.py — all-HDR spike against conjunction+floor, on the SAME ten disjoint windows.

    cell          K562
    seed space    200 seeds: 200-299 + 400-499 (JOINED, non-contiguous). Nothing outside it is
                  considered -- not the cut, not the clean set, not the dirty regime.
    windows       10 x width 20 at stride 20  ->  10*20 = 200, a CIRCULAR tiling with
                  width == stride, so the ten windows PARTITION the space. No overlap, no gap.
    arm A         conjunction + floor: Cas12a min-unions CUT over the whole 200-seed space and the
                  band is pinned inside this hotkey's 20 seeds (`band_candidates`). Narrow cut, so
                  `band` is a subset of `cut` by construction.
    arm B         all-HDR spike: the Cas12a group min-unions HDR over the 20-seed window itself,
                  and the clean band is what survives. `build_for_cell(seed_list=...)` scales
                  `main_max_fail` to span 20 exactly as the live miner would.
    tasks         every K562 round whose three seeds include one in 200-299 AND one in 400-499.

**Both arms are priced identically and from the SHIPPED ROWS**, which is the only honest way to
compare two constructions that build over different spaces (CLAUDE.md: `meta["clean"]` counts
within each arm's own space and is biased across arms). For each build:

  * `band`      = every seed of the 200 where EVERY row satisfies `hdr`  -> stage 4 pins all three
                  targets, consistency exactly 1.0
  * `cut-clean` = the same on `cut`, minus the band -> only `is_cut` pinned, the floor regime
  * the round's three REAL seeds are then scored through the validator's own stages, so the
    reported consistency is measured per seed rather than assumed from the value ladder.

That last point matters here: CLAUDE.md measures a cut-clean seed at 0.237 in an all-cut row
composition and 0.162 in an all-HDR one, so assuming one v_clean for both arms would hand the
spike arm a floor it does not have.

    HC_JSON=hdr_vs_conj.json python hdr_vs_conj.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "hdrvsconj")

import dataclasses, gc, json, logging, time                            # noqa: E402

logging.basicConfig(level=logging.ERROR)
import numpy as np                                                    # noqa: E402

import genExp as G                                                    # noqa: E402
from niome_subnet.genomics import conjunction as CJ                   # noqa: E402
from niome_subnet.genomics import all_hdr as AH                       # noqa: E402
from niome_subnet.genomics import mt19937 as MT                       # noqa: E402
from niome_subnet.genomics.validation import _parse_seeds             # noqa: E402
from conj_stageb import API                                           # noqa: E402
from sd_task import fetch, score                                      # noqa: E402
from widecut_price import records_of                                  # noqa: E402

CELL = os.getenv("HC_CELL", "K562")
WIN_W = int(os.getenv("HC_WIDTH", "20"))
WIN_S = int(os.getenv("HC_STRIDE", "20"))
OUT = os.getenv("HC_JSON", "hdr_vs_conj.json")
MIN_FREE_GB = float(os.getenv("HC_MIN_FREE_GB", "7"))
BLOCK_A, BLOCK_B = (200, 300), (400, 500)
SPACE = sorted(list(range(*BLOCK_A)) + list(range(*BLOCK_B)))
NWIN = len(SPACE) // WIN_S
assert NWIN * WIN_S == len(SPACE), f"{NWIN}x{WIN_S} != {len(SPACE)}: not an exact tiling"
assert WIN_W == WIN_S, "width must equal stride for the windows to partition the space"


def free_gb() -> float:
    m = {}
    for line in open("/proc/meminfo"):
        p = line.split()
        m[p[0].rstrip(":")] = int(p[1])
    return m.get("MemAvailable", 0) / 1048576.0


def wait_for_memory(tag: str = "") -> None:
    for _ in range(120):
        if free_gb() >= MIN_FREE_GB:
            return
        print(f"      waiting for memory ({free_gb():.1f} GB free) {tag}", flush=True)
        time.sleep(30)


def compliant(rows, contract, cell_types, ctx, rule):
    """Every seed of SPACE on which EVERY row satisfies `rule` — arm-agnostic, from the rows."""
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


def price(rows, contract, reference, cell_types, ctx, seeds):
    """Band, cut-clean and the round's REAL per-seed consistency, all from the shipped rows."""
    band = compliant(rows, contract, cell_types, ctx, "hdr")
    cut = compliant(rows, contract, cell_types, ctx, "cut")
    bs, cs = set(band), set(cut) - set(band)
    per = []
    for s in seeds:
        r = score(rows, contract, reference, cell_types, seed=s)
        per.append({"seed": s, "cons": r["consistency"], "weighted": r["weighted"],
                    "fidelity": r["fidelity"],
                    "regime": "band" if s in bs else ("clean" if s in cs else "dirty")})
    cons = float(np.mean([p["cons"] for p in per]))
    wxf = per[0]["weighted"] * per[0]["fidelity"]
    return {"band_n": len(band), "band": band, "cut_n": len(cut), "clean_n": len(cs),
            "band_hits": [s for s in seeds if s in bs], "clean_hits": [s for s in seeds if s in cs],
            "per_seed": per, "round_cons": cons, "weighted": per[0]["weighted"],
            "fidelity": per[0]["fidelity"], "round_final": wxf * cons}


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
    windows = [SPACE[i * WIN_S:i * WIN_S + WIN_W] for i in range(NWIN)]
    out = json.load(open(OUT)) if os.path.exists(OUT) else []
    done = {r["task"] for r in out}
    print(f"=== {CELL} | {NWIN} windows x width {WIN_W} at stride {WIN_S} partitioning "
          f"{BLOCK_A[0]}-{BLOCK_A[1]-1} + {BLOCK_B[0]}-{BLOCK_B[1]-1} ({len(SPACE)} seeds) | "
          f"conjunction cut {len(SPACE)} narrow, k={base.band_k} group={base.group_size} | "
          f"all-HDR band space = the window | {len(tasks)} qualifying tasks ===", flush=True)
    for i, w in enumerate(windows):
        print(f"    w{i}  {w[0]}-{w[-1]}", flush=True)
    print(flush=True)

    mf = max(1, round(base.cas12a_max_fail * len(SPACE) / 900))
    for ti, t in enumerate(tasks, 1):
        tid = (t.get("task_id") or t["id"])[:8]
        if tid in done:
            continue
        contract, reference = dict(t["content"]["contract"]), t["content"]["hbb_reference"]
        seeds = _parse_seeds(contract["seed"])
        created = t.get("created_at", "")
        t0 = time.monotonic()
        c0 = dict(contract, seed=0)
        ctx = G.build_context(c0, reference, cell_types)
        print(f"[{ti}/{len(tasks)}] {tid} {created[:16]} seeds {seeds}", flush=True)
        wins = []
        for i, w in enumerate(windows):
            rec = {"win": i, "range": [w[0], w[-1]], "in_window": [s for s in seeds if s in set(w)]}
            # --- arm A: conjunction + floor -------------------------------------------------
            wait_for_memory(f"{tid} w{i} conj")
            cfgC = dataclasses.replace(base, seed_list=tuple(SPACE), start_seed=SPACE[0],
                                       end_seed=SPACE[-1], cas12a_max_fail=mf,
                                       band_candidates=tuple(w))
            rows, meta = CJ.build_submission(c0, reference, cell_types, cfg=cfgC, budget_s=1800)
            MT.free_gpu_memory()
            if not rows:
                rec["conj"] = {"declined": meta.get("reason")}
            else:
                rec["conj"] = dict(price(rows, contract, reference, cell_types, ctx, seeds),
                                   pool=meta.get("pool"), cas9=meta.get("cas9_pool"),
                                   band_cfg=sorted(int(x) for x in meta.get("band_seeds", [])))
            del rows; gc.collect()
            # --- arm B: all-HDR spike -------------------------------------------------------
            wait_for_memory(f"{tid} w{i} hdr")
            rows, meta = AH.build_for_cell(c0, reference, cell_types, budget_s=1800,
                                           seed_list=tuple(w))
            MT.free_gpu_memory()
            if not rows:
                rec["hdr"] = {"declined": meta.get("reason")}
            else:
                rec["hdr"] = dict(price(rows, contract, reference, cell_types, ctx, seeds),
                                  cas9=meta.get("cas9_pool"), bank=meta.get("bank"),
                                  clean_cfg=meta.get("clean"))
            del rows; gc.collect()
            c, h = rec["conj"], rec["hdr"]
            print(f"    w{i} {w[0]}-{w[-1]}  "
                  f"conj " + (f"DECLINED" if "declined" in c else
                              f"band {c['band_n']:3d} clean {c['clean_n']:3d} "
                              f"hits {len(c['band_hits'])}b/{len(c['clean_hits'])}c "
                              f"cons {c['round_cons']:.4f} final {c['round_final']:6.1f}")
                  + "  |  hdr " + (f"DECLINED" if "declined" in h else
                                   f"band {h['band_n']:3d} clean {h['clean_n']:3d} "
                                   f"hits {len(h['band_hits'])}b/{len(h['clean_hits'])}c "
                                   f"cons {h['round_cons']:.4f} final {h['round_final']:6.1f}"),
                  flush=True)
            wins.append(rec)
        out.append({"cell": CELL, "task": tid, "created_at": created, "seeds": seeds,
                    "space": [BLOCK_A, BLOCK_B], "width": WIN_W, "stride": WIN_S,
                    "k": base.band_k, "cut": len(SPACE), "windows": wins})
        json.dump(out, open(OUT, "w"), indent=1)
        print(f"    ({time.monotonic()-t0:.0f}s)\n", flush=True)
        gc.collect()
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
