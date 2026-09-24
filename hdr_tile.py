#!/usr/bin/env python3
"""hdr_tile.py — the all-HDR SPIKE alone, on N circular windows of one 200-seed joined space.

    cell          K562
    seed space    200 seeds: 200-299 + 400-499 (JOINED). Nothing outside it is considered --
                  not the band space, not the clean set, not the dirty regime.
    windows       HT_NWIN slices of width HT_WIDTH at stride HT_STRIDE, CIRCULAR
                  (`count * stride == span`), so they wrap the space exactly once.
    arm           all-HDR only. `build_for_cell(seed_list=...)` scales `main_max_fail` to the
                  window's span exactly as the live miner does.
    tasks         every K562 round whose three seeds include one in 200-299 AND one in 400-499.

Companion to `hdr_vs_conj.py`, which ran width 20 at stride 20 (disjoint) against the
conjunction. Here only the width changes, so the question is what the SPIKE alone is worth as its
window grows: a wider window gives the min-union more seeds to fail on, which SHRINKS the clean
band (CLAUDE.md measures 13 at width 100 against this test's 14-15 at width 20), while the
overlap it buys means each seed is covered by several hotkeys at once.

Priced from the SHIPPED ROWS, identically to `hdr_vs_conj.py`: `band` is every seed of the 200 on
which EVERY row satisfies `hdr`, `cut-clean` the same on `cut`, and the round's three real seeds
are scored through the validator's own stages rather than assumed from the value ladder.

    HT_WIDTH=100 HT_STRIDE=20 HT_JSON=hdr_tile_w100.json python hdr_tile.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "hdrtile")

import gc, json, logging, time                                        # noqa: E402

logging.basicConfig(level=logging.ERROR)
import itertools                                                      # noqa: E402
import numpy as np                                                    # noqa: E402

import genExp as G                                                    # noqa: E402
from niome_subnet.genomics import all_hdr as AH                       # noqa: E402
from niome_subnet.genomics import conjunction as CJ                   # noqa: E402
from niome_subnet.genomics import mt19937 as MT                       # noqa: E402
from niome_subnet.genomics.validation import _parse_seeds             # noqa: E402
from conj_stageb import API                                           # noqa: E402
from sd_task import fetch, score                                      # noqa: E402
from widecut_price import records_of                                  # noqa: E402

CELL = os.getenv("HT_CELL", "K562")
WIN_W = int(os.getenv("HT_WIDTH", "100"))
WIN_S = int(os.getenv("HT_STRIDE", "20"))
OUT = os.getenv("HT_JSON", "hdr_tile.json")
MIN_FREE_GB = float(os.getenv("HT_MIN_FREE_GB", "7"))
BLOCK_A, BLOCK_B = (200, 300), (400, 500)
# TWO spaces, as in `conj_tile.py`. BAND_SPACE is what the windows tile; SPACE is what `band` and
# `cut-clean` are enumerated over. Separating them lets an all-HDR arm sit on one class while
# still being priced across the cut the rest of the fleet uses.
BAND_SPACE = sorted(int(x) for c in os.getenv("HT_BAND_CLASSES", "1,3").split(",")
                    for x in range(int(c) * 100 + 100, int(c) * 100 + 200))
SPACE = sorted(int(x) for c in os.getenv("HT_CLASSES", "1,3").split(",")
               for x in range(int(c) * 100 + 100, int(c) * 100 + 200))
NWIN = int(os.getenv("HT_NWIN", str(len(BAND_SPACE) // WIN_S)))
assert NWIN * WIN_S >= len(BAND_SPACE), \
    (f"{NWIN} windows x stride {WIN_S} does not reach {len(BAND_SPACE)} seeds: some seed sits in "
     f"no window")
assert (NWIN - 1) * WIN_S < len(BAND_SPACE), \
    "the last window wraps a full revolution onto another window's slice"
assert WIN_W <= len(BAND_SPACE), "the window is wider than the band space"
assert set(BAND_SPACE) <= set(SPACE), "the band space must sit inside the cut space"


# Task selection. "pairs" keeps every K562 round whose seeds contain one of the five class-pairs
# 200/200, 200/400, 200/600, 400/600, 400/400 -- any two seeds drawn from classes 2xx/4xx/6xx
# except the 600/600 pair, which is deliberately not in the set. A round can match several pairs.
# "both" is the older filter (one seed in 200-299 AND one in 400-499), kept so the earlier runs in
# this series stay reproducible against their own task list.
TASK_FILTER = os.getenv("TASK_FILTER", "pairs")
_PAIRS = {(2, 2), (2, 4), (2, 6), (4, 6), (4, 4)}
# An explicit task list overrides the class-pair filter entirely: comma-separated id prefixes,
# matched against the full task id, for probing one round rather than a sampled set.
_ONLY = tuple(x.strip() for x in os.getenv("HT_TASK", "").split(",") if x.strip())


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

    windows = [sorted({BAND_SPACE[(i * WIN_S + k) % len(BAND_SPACE)] for k in range(WIN_W)})
               for i in range(NWIN)]
    out = json.load(open(OUT)) if os.path.exists(OUT) else []
    done = {r["task"] for r in out}
    print(f"=== {CELL} all-HDR | {NWIN} band windows x width {WIN_W} at stride {WIN_S} over "
          f"{BAND_SPACE[0]}-{BAND_SPACE[-1]} ({len(BAND_SPACE)} seeds) | band+clean priced over "
          f"{len(SPACE)} seeds | seats {NWIN*WIN_W} | {len(tasks)} qualifying tasks ===",
          flush=True)
    for i, w in enumerate(windows):
        print(f"    w{i}  {len(w):3d} seeds  {w[0]}..{w[-1]}", flush=True)
    print(flush=True)

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
            rows, meta = AH.build_for_cell(c0, reference, cell_types, budget_s=1800,
                                           seed_list=tuple(w))
            MT.free_gpu_memory()
            if not rows:
                rec["hdr"] = {"declined": meta.get("reason")}
                print(f"    w{i} DECLINED {meta.get('reason')}", flush=True)
            else:
                rec["hdr"] = dict(price(rows, contract, reference, cell_types, ctx, seeds),
                                  cas9=meta.get("cas9_pool"), bank=meta.get("bank"),
                                  clean_cfg=meta.get("clean"))
                h = rec["hdr"]
                print(f"    w{i} {w[0]}-{w[-1]} ({len(w)})  band {h['band_n']:3d} "
                      f"clean {h['clean_n']:3d}  hits {len(h['band_hits'])}b/"
                      f"{len(h['clean_hits'])}c  cons {h['round_cons']:.4f} "
                      f"final {h['round_final']:6.1f}", flush=True)
            del rows; gc.collect()
            wins.append(rec)
        out.append({"cell": CELL, "task": tid, "created_at": created, "seeds": seeds,
                    "space": [BLOCK_A, BLOCK_B], "width": WIN_W, "stride": WIN_S,
                    "nwin": NWIN, "method": "all-hdr", "windows": wins})
        json.dump(out, open(OUT, "w"), indent=1)
        print(f"    ({time.monotonic()-t0:.0f}s)\n", flush=True)
        gc.collect()
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
