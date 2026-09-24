#!/usr/bin/env python3
"""loop_axis.py — one fixed 300-seed window, re-banded N times: when are the round's seeds reached?

The configuration under test, as specified:
    joined window / cut   900 (wide) -- FORCED here for every cell, including HEK293, whose live
                          `CUT_MODE` is "union"; the point is to hold the cut constant across cells
    sub window            300 seeds
    stride                0  -> every hotkey shares ONE window, so the fleet collapses to a single
                          distinct band per loop and this is the pure LOOP axis
    everything else       each cell's live CELL_CONFIG (k, group_size, light_cell_rows, cell-aware)

**A structural limit to read before the results.** The band is confined to the 300-seed window, so a
seed outside it can NEVER be found however many loops run. Seeds are uniform over 100-999, so on
average only 1 of 3 is reachable and `P(>=2 reachable) = 7/27 ~ 26%`. Expect most tasks to report the
second seed as "never" -- that is the cost of concentrating on 300 seeds, not a failure of the loop.

Cheap by design: only the BAND is needed to answer "was this seed reached", so `choose_band` is
enough and no full `build_submission` runs. Compliance is computed once over the 300-seed window
(not 900), which also keeps memory ~10x below the sweeps that were OOM-killed earlier today.

    LA_CELLS=HUDEP-2,HEK293,K562,CD34+_HSPC LA_N=10 LA_LOOPS=10 python loop_axis.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "loopaxis")

import dataclasses, json, logging, random, time                        # noqa: E402

logging.basicConfig(level=logging.ERROR)
import numpy as np                                                    # noqa: E402

import genExp as G                                                    # noqa: E402
import joined_window as JW                                            # noqa: E402
from niome_subnet.genomics import conjunction as CJ                   # noqa: E402
from niome_subnet.genomics import mt19937 as MT                       # noqa: E402
from niome_subnet.genomics.validation import _parse_seeds             # noqa: E402
from conj_stageb import API                                           # noqa: E402
from sd_task import fetch                                             # noqa: E402

CELLS = [c for c in os.getenv("LA_CELLS", "HUDEP-2,HEK293,K562,CD34+_HSPC").split(",") if c]
N = int(os.getenv("LA_N", "10"))
LOOPS = int(os.getenv("LA_LOOPS", "10"))
WIDTH = int(os.getenv("LA_WIDTH", "300"))
# "oracle" = the REAL width-100 classes the round's three seeds drew from, unioned. That is
# not available at build time -- it is the perfect-prediction case -- and it exists to remove
# the reachability limit: under "uniform" the window is a fixed 100-399 and a seed outside it
# can never be reached at any loop count, so ~74% of tasks could not report a 2-seed number at
# all. With the oracle window all three seeds are in-window by construction and every task
# yields a loop count, which is what makes the cells comparable.
# Note the window is NOT always 300 seeds: two seeds sharing a class give 200, three give 100.
# `classes` is recorded per task so a short window is never mistaken for a fast one.
MODE = os.getenv("LA_WINDOW", "oracle")
OUT = os.getenv("LA_JSON", "loop_axis.json")


def prepare(contract, reference, cell_types, cell, window):
    base = CJ.config_for(cell)
    cut = list(JW.FULL_SPACE)                 # 900, forced for every cell
    cfg = dataclasses.replace(base, seed_list=tuple(cut), band_candidates=tuple(window))
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    path = os.path.join(CJ.BANK_DIR, f"cas12a-{CJ.bank_key(contract, cell_types, cfg)}.npz")
    if not os.path.exists(path):
        os.makedirs(CJ.BANK_DIR, exist_ok=True)
        bank = CJ.build_bank(contract, reference, cell_types, ctx, sites, cfg, None)
        MT.free_gpu_memory()
        if not bank:
            return None
        CJ.save_bank(path, bank)
    records = CJ.load_bank(path, limit=cfg.bank_keep)
    ok = CJ.hdr_compliance(records, contract, cell_types, ctx, window, cfg.band_rule)
    MT.free_gpu_memory()
    # numpy-backed: a set of ints costs ~16x an int32 array and OOM-killed three runs today.
    ok = {i: np.fromiter(v, dtype=np.int32) for i, v in ok.items()}
    cell_ok = None
    if cfg.band_cell_aware:
        cell_ok = CJ.cas9_cell_probe(contract, cell_types, ctx, sites, cfg, window, cfg.band_rule)
        MT.free_gpu_memory()
    return {"cfg": cfg, "ok": ok, "cell_ok": cell_ok, "bank": len(records)}


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    items = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])
    fixed = CJ.sub_window(list(JW.FULL_SPACE), WIDTH, 0.0)
    out = json.load(open(OUT)) if os.path.exists(OUT) else []
    done = {(r["cell"], r["task"]) for r in out}
    print(f"=== loop axis | window mode {MODE!r} | cut 900 forced | {LOOPS} loops | "
          f"cells {CELLS} ===", flush=True)
    if MODE == "oracle":
        print("    window = the real width-100 classes the 3 seeds drew from (perfect prediction);"
              "\n    all 3 seeds in-window by construction, so every task yields a loop count.\n",
              flush=True)
    else:
        print(f"    window fixed at {min(fixed)}-{max(fixed)}; P(>=2 of 3 reachable) = 26%\n",
              flush=True)

    for cell in CELLS:
        base = CJ.config_for(cell)
        tasks = []
        for t in sorted(items, key=lambda x: x.get("created_at", ""), reverse=True):
            c = (t.get("content") or {}).get("contract") or {}
            if c.get("cell_type") != cell or len(_parse_seeds(c.get("seed"))) != 3:
                continue
            tasks.append(t)
            if len(tasks) >= N:
                break
        print(f"--- {cell}  k={base.band_k} group={base.group_size} "
              f"light={base.light_cell_rows} | {len(tasks)} tasks ---", flush=True)
        for t in tasks:
            tid = (t.get("task_id") or t["id"])
            if (cell, tid[:8]) in done:
                continue
            contract = dict(t["content"]["contract"])
            seeds = _parse_seeds(contract["seed"])
            if MODE == "oracle":
                real = sorted({(s - 100) // 100 for s in seeds})
                # PAD to three classes when the seeds share one. Two seeds in the same width-100
                # class happens on ~24% of rounds and three on ~1.2%, and without padding those
                # tasks would run on a 200- or 100-seed window -- covering far faster per loop and
                # reporting a low loop count for a reason that has nothing to do with the band. The
                # padding classes are drawn from the 6-8 the seeds did NOT use, so they add window
                # size without adding reachable seeds, which is exactly the control wanted.
                # Deterministic per task so a re-run reproduces the same window.
                rng = random.Random(int(tid[:8], 16))
                pad = [w for w in range(9) if w not in real]
                rng.shuffle(pad)
                classes = sorted(real + pad[:max(0, 3 - len(real))])
                padded = sorted(set(classes) - set(real))
                window = sorted(s for w in classes for s in range(w*100+100, w*100+200))
            else:
                real, classes, padded, window = None, None, None, fixed
            inwin = [s for s in seeds if s in set(window)]
            t0 = time.monotonic()
            prep = prepare(dict(contract, seed=0), t["content"]["hbb_reference"],
                           cell_types, cell, window)
            if prep is None:
                print(f"  {tid[:8]} bank empty, skipped", flush=True)
                continue
            cfg = prep["cfg"]
            covered, found = set(), {}
            for loop in range(1, LOOPS + 1):
                cands = [s for s in window if s not in covered]
                if len(cands) < cfg.band_k:
                    break
                band, alive = CJ.choose_band(prep["ok"], cfg.band_k, cands, cfg.group_size,
                                             prep["cell_ok"], {}, cfg.band_quota)
                if len(band) < cfg.band_k or len(alive) < cfg.group_size:
                    break
                covered |= {int(x) for x in band}
                for s in seeds:
                    if s in covered and s not in found:
                        found[s] = loop
            v = sorted(found.values())
            two = v[1] if len(v) >= 2 else None
            print(f"  {tid[:8]} seeds {seeds}  win {len(window)}sd"
                  f"{' real' + str(len(real)) + '+pad' + str(len(padded)) if real is not None else ''}  "
                  f"found {({s: found[s] for s in seeds if s in found}) or '{}'}  "
                  f"2-seed L{two if two else '-'}  covered {len(covered)}/{len(window)}  "
                  f"({time.monotonic()-t0:.0f}s)", flush=True)
            out.append({"cell": cell, "task": tid[:8], "seeds": seeds, "in_window": inwin,
                        "classes": classes, "real_classes": real, "padded_classes": padded,
                        "window_size": len(window), "mode": MODE,
                        "found": {str(s): found.get(s) for s in seeds}, "two_seed_loop": two,
                        "covered": len(covered), "window": [min(window), max(window)]})
            json.dump(out, open(OUT, "w"), indent=1)

        cr = [r for r in out if r["cell"] == cell]
        if cr:
            twos = [r["two_seed_loop"] for r in cr if r["two_seed_loop"]]
            ones = [min(x for x in r["found"].values() if x) for r in cr
                    if any(r["found"].values())]
            print(f"    {cell}: 2 seeds reached on {len(twos)}/{len(cr)} tasks"
                  + (f", loops {sorted(twos)}, median {int(np.median(twos))}" if twos else "")
                  + f" | >=1 seed on {len(ones)}/{len(cr)}"
                  + (f", median L{int(np.median(ones))}" if ones else ""), flush=True)
    print(f"\nwrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
