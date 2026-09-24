#!/usr/bin/env python3
"""conj_tile.py — conjunction+floor on N narrow windows of ONE 100-seed class.

    cell          K562
    seed space    CT_CLASSES, default the single class 200-299 (100 seeds). Nothing outside it is
                  considered -- not the cut, not the clean set, not the dirty regime.
    windows       CT_NWIN slices of width CT_WIDTH at stride CT_STRIDE over that space, CIRCULAR
                  (they wrap), so the union covers the space even when count*stride overshoots it.
    cut           the whole space (NARROW), so `band` is a subset of `cut` by construction and the
                  Cas12a min-union optimises cut-cleanness INSIDE the window rather than across 900.
    arm           conjunction only, at the cell's live CELL_CONFIG (k, group, light, cell-aware).
    tasks         every K562 round whose three seeds include one in 200-299 AND one in 400-499 --
                  the SAME six rounds every other arm in this series ran, so the comparison is
                  paired even though this space is half the size.

The question this asks that no earlier arm did: `band_k` 11 has to be found inside a **15-seed**
candidate window. `choose_band` is a greedy prefix, so a narrow window costs it the freedom to
pick easy seeds -- CLAUDE.md records the band GROWING toward the window as the window shrinks
(band 10 of 10 at width 10 on HEK293), which makes the conditional Cas9 fill harder, not easier.
Whether k=11 of 15 forms at all is the first result; the clean set is the second.

Priced from the SHIPPED ROWS, identically to `hdr_vs_conj.py` and `hdr_tile.py`: `band` is every
seed of the space on which EVERY row satisfies `hdr`, `cut-clean` the same on `cut`, and the
round's three real seeds are scored through the validator's own stages.

    CT_WIDTH=15 CT_STRIDE=15 CT_NWIN=7 CT_JSON=conj_tile_w15.json python conj_tile.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "conjtile")

import dataclasses, gc, json, logging, time                            # noqa: E402

logging.basicConfig(level=logging.ERROR)
import itertools                                                      # noqa: E402
import numpy as np                                                    # noqa: E402

import genExp as G                                                    # noqa: E402
from niome_subnet.genomics import conjunction as CJ                   # noqa: E402
from niome_subnet.genomics import mt19937 as MT                       # noqa: E402
from niome_subnet.genomics.validation import _parse_seeds             # noqa: E402
from conj_stageb import API                                           # noqa: E402
from sd_task import fetch, score                                      # noqa: E402
from widecut_price import records_of                                  # noqa: E402

CELL = os.getenv("CT_CELL", "K562")
WIN_W = int(os.getenv("CT_WIDTH", "15"))
WIN_S = int(os.getenv("CT_STRIDE", "15"))
NWIN = int(os.getenv("CT_NWIN", "7"))
OUT = os.getenv("CT_JSON", "conj_tile_w15.json")
MIN_FREE_GB = float(os.getenv("CT_MIN_FREE_GB", "7"))
# Band depth and group, overriding the cell's live CELL_CONFIG. Both default to the shipped value,
# so an unset env reproduces the live arm exactly. `band_k` is the quantity CLAUDE.md puts AT the
# erythroid wall of 12 -- one step past 11 the surviving pool falls ~0.57x, which is the whole
# question this sweep asks at a 20-seed candidate window.
K_OVERRIDE = os.getenv("CT_K")
GROUP_OVERRIDE = os.getenv("CT_GROUP")
# LOOP AXIS mode. Every hotkey shares the WHOLE band space and hotkey i plays loop i+1, excluding
# every seed loops 1..i banded -- so the bands are disjoint by construction and the union is
# exactly `hotkeys * band_k` until the chain runs out of candidates. Off by default, in which case
# the hotkeys take offset window slices as usual.
LOOP_AXIS = os.getenv("CT_LOOP_AXIS", "") not in ("", "0", "false")
# An explicit task list, as in hdr_tile.py: comma-separated id prefixes.
_ONLY = tuple(x.strip() for x in os.getenv("CT_TASK", "").split(",") if x.strip())
# TWO spaces, deliberately different.
#
#   BAND_SPACE  the seeds the windows tile -- 200-299 by default, one class.
#   SPACE       the seeds the CUT min-unions over and the clean set is counted in -- 200-299 PLUS
#               400-499 by default, both classes.
#
# Separating them is the point of this arm. The band lives in one class, but the floor is a
# property of the CUT, and a cut restricted to 200-299 can never make a 400-499 seed clean --
# which would leave a third of each round's seeds unreachable by either regime. Counting clean
# over both classes is also what makes this comparable to every earlier 200-seed arm.
BAND_SPACE = sorted(int(x) for c in os.getenv("CT_BAND_CLASSES", "1").split(",")
                    for x in range(int(c) * 100 + 100, int(c) * 100 + 200))
SPACE = sorted(int(x) for c in os.getenv("CT_CLASSES", "1,3").split(",")
               for x in range(int(c) * 100 + 100, int(c) * 100 + 200))
BLOCK_A, BLOCK_B = (200, 300), (400, 500)
assert WIN_W <= len(BAND_SPACE), "the window is wider than the band space"
assert set(BAND_SPACE) <= set(SPACE), "the band space must sit inside the cut space"


def _spans(seeds) -> str:
    """"200-299 + 400-499 + 600-699" for a possibly non-contiguous seed list."""
    xs = sorted(seeds)
    out, lo, prev = [], xs[0], xs[0]
    for x in xs[1:]:
        if x != prev + 1:
            out.append(f"{lo}-{prev}"); lo = x
        prev = x
    out.append(f"{lo}-{prev}")
    return " + ".join(out)


# Task selection. "pairs" keeps every K562 round whose seeds contain one of the five class-pairs
# 200/200, 200/400, 200/600, 400/600, 400/400 -- any two seeds drawn from classes 2xx/4xx/6xx
# except the 600/600 pair, which is deliberately not in the set. A round can match several pairs.
# "both" is the older filter (one seed in 200-299 AND one in 400-499), kept so the earlier runs in
# this series stay reproducible against their own task list.
TASK_FILTER = os.getenv("TASK_FILTER", "pairs")
_PAIRS = {(2, 2), (2, 4), (2, 6), (4, 6), (4, 4)}


def _qualifies(sd) -> bool:
    if len(sd) != 3:
        return False
    if TASK_FILTER == "both":
        return any(200 <= s < 300 for s in sd) and any(400 <= s < 500 for s in sd)
    cls = sorted(s // 100 for s in sd if s // 100 in (2, 4, 6))
    return any(tuple(sorted(p)) in _PAIRS for p in itertools.combinations(cls, 2))


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
            "n_offband_in_cut": len([s for s in seeds if s in set(SPACE) and s not in bs]),
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
        if _ONLY:
            tid_full = t.get("task_id") or t["id"]
            if not any(tid_full.startswith(x) for x in _ONLY):
                continue
        elif not _qualifies(sd):
            continue
        tasks.append(t)

    base = CJ.config_for(CELL)
    if K_OVERRIDE:
        base = dataclasses.replace(base, band_k=int(K_OVERRIDE))
    if GROUP_OVERRIDE:
        base = dataclasses.replace(base, group_size=int(GROUP_OVERRIDE))
    windows = ([sorted(BAND_SPACE)] * NWIN if LOOP_AXIS else
               [sorted({BAND_SPACE[(i * WIN_S + k) % len(BAND_SPACE)] for k in range(WIN_W)})
                for i in range(NWIN)])
    seats = NWIN * WIN_W
    union = len({s for w in windows for s in w})
    nb = len(BAND_SPACE)
    out = json.load(open(OUT)) if os.path.exists(OUT) else []
    done = {r["task"] for r in out}
    print(f"=== {CELL} conjunction+floor | {NWIN} band windows x width {WIN_W} at stride "
          f"{WIN_S} over {BAND_SPACE[0]}-{BAND_SPACE[-1]} ({nb} seeds) | CUT and CLEAN over "
          f"{len(SPACE)} seeds ({_spans(SPACE)}) | "
          f"k={base.band_k} group={base.group_size} | "
          + (f"LOOP AXIS: {NWIN} loops x k={base.band_k} = {NWIN*base.band_k} seats from {nb} "
             f"candidates" if LOOP_AXIS else f"seats {seats}, union {union}/{nb}, "
             f"{seats - union} duplicated")
          + f" | {len(tasks)} qualifying tasks ===", flush=True)
    for i, w in enumerate(windows):
        print(f"    w{i}  {len(w):2d} seeds  {w[0]}..{w[-1]}", flush=True)
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
            rec = {"win": i, "n_window": len(w), "lo": w[0], "hi": w[-1],
                   "in_window": [s for s in seeds if s in set(w)]}
            wait_for_memory(f"{tid} w{i}")
            cfg = dataclasses.replace(base, seed_list=tuple(SPACE), start_seed=SPACE[0],
                                      end_seed=SPACE[-1], cas12a_max_fail=mf,
                                      band_candidates=tuple(w),
                                      band_loop=(i + 1) if LOOP_AXIS else 1)
            rows, meta = CJ.build_submission(c0, reference, cell_types, cfg=cfg, budget_s=1800)
            MT.free_gpu_memory()
            if not rows:
                rec["conj"] = {"declined": meta.get("reason"), "band_reached": meta.get("band")}
                print(f"    w{i} {w[0]}-{w[-1]} DECLINED  {meta.get('reason')}", flush=True)
            else:
                rec["conj"] = dict(price(rows, contract, reference, cell_types, ctx, seeds),
                                   pool=meta.get("pool"), cas9=meta.get("cas9_pool"),
                                   band_cfg=sorted(int(x) for x in meta.get("band_seeds", [])))
                c = rec["conj"]
                print(f"    w{i} {w[0]}-{w[-1]} ({len(w)})  band {c['band_n']:3d} "
                      f"clean {c['clean_n']:3d}  hits {len(c['band_hits'])}b/"
                      f"{len(c['clean_hits'])}c  cons {c['round_cons']:.4f} "
                      f"final {c['round_final']:6.1f}  pool {c['pool']}", flush=True)
            del rows; gc.collect()
            wins.append(rec)
        out.append({"cell": CELL, "task": tid, "created_at": created, "seeds": seeds,
                    "space": [SPACE[0], SPACE[-1]], "n_space": len(SPACE),
                    "band_space": [BAND_SPACE[0], BAND_SPACE[-1]], "n_band_space": nb,
                    "width": WIN_W,
                    "stride": WIN_S, "nwin": NWIN, "method": "conjunction",
                    "k": base.band_k, "group": base.group_size, "cut": len(SPACE),
                    "windows": wins})
        json.dump(out, open(OUT, "w"), indent=1)
        print(f"    ({time.monotonic()-t0:.0f}s)\n", flush=True)
        gc.collect()
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
