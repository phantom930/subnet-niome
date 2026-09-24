#!/usr/bin/env python3
"""loop_clean.py — when a loop-axis band HITS a round seed, what does that same loop's CLEAN set hold?

`loop_axis.py` answered "which loop reaches a seed" using `choose_band` alone, which is all the BAND
question needs. It cannot answer the floor question: a hotkey that catches one seed in its band
still scores the other two somewhere, and whether those land in its cut-clean set (~0.21) or dirty
(~0.10) is the difference between a round that places and one that does not.

    cons = (band_hits*1.0 + clean_hits*v_clean + dirty*v_dirty) / 3

So this replays each task's loop chain -- deterministically, the same `choose_band` calls -- and at
every loop where a seed was FOUND runs ONE full `build_submission` to enumerate that loop's clean
set, then counts all three round seeds against both regimes.

**How the band is pinned, since `build_submission` has no force parameter.** `band_candidates` is
narrowed to exactly the k seeds that loop's band came out as. `choose_band` is a greedy prefix that
takes a candidate while the surviving pool holds `group_size`; the original chain already proved the
intersection over all k clears that bar, and intersections only shrink, so every prefix clears it
too and the greedy takes all k. The band set is therefore reproduced exactly -- asserted per build,
not assumed.

One bank per (contract, cell): the cut space is 900 for every loop and `bank_key` excludes the band
params, so only the first build of a task can pay a cold scan, and the loop_axis run already warmed
most of them.

    LC_JSON=loop_clean.json python loop_clean.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "loopclean")

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
from loop_axis import prepare                                         # noqa: E402

SRC = os.getenv("LC_SRC", "loop_axis_oracle.json")
OUT = os.getenv("LC_JSON", "loop_clean.json")
LOOPS = int(os.getenv("LC_LOOPS", "10"))
CELLS = [c for c in os.getenv("LC_CELLS", "HUDEP-2,HEK293,K562,CD34+_HSPC").split(",") if c]

# The measured per-seed value ladder (CLAUDE.md, "What the competition actually does"): a band seed
# pins all three stage-4 targets and scores exactly 1.0; a cut-clean seed pins `is_cut` only.
V_CLEAN = {"HEK293": 0.212, "_": 0.212}
V_DIRTY = {"HEK293": 0.082, "_": 0.104}


def main():
    src = {(r["cell"], r["task"]): r for r in json.load(open(SRC))}
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    items = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])
    by_id = {(t.get("task_id") or t["id"])[:8]: t for t in items}

    out = json.load(open(OUT)) if os.path.exists(OUT) else []
    done = {(r["cell"], r["task"], r["loop"]) for r in out}

    targets = []
    for (cell, tid), r in src.items():
        if cell not in CELLS:
            continue
        loops = sorted({l for l in r["found"].values() if l})
        if loops:
            targets.append((cell, tid, loops, r))
    targets.sort(key=lambda x: (CELLS.index(x[0]), x[1]))
    n_pts = sum(len(t[2]) for t in targets)
    print(f"=== loop clean sets | {len(targets)} tasks, {n_pts} hit-loops | cut 900 forced ===",
          flush=True)

    for cell, tid, loops, r in targets:
        if all((cell, tid, l) in done for l in loops):
            continue
        t = by_id.get(tid)
        if t is None:
            print(f"  {cell} {tid} not in listing, skipped", flush=True)
            continue
        contract, reference = dict(t["content"]["contract"]), t["content"]["hbb_reference"]
        seeds = _parse_seeds(contract["seed"])
        window = sorted(s for w in r["classes"] for s in range(w * 100 + 100, w * 100 + 200))
        print(f"\n  {cell} {tid}  seeds {seeds}  win {len(window)}sd  hit-loops {loops}", flush=True)

        prep = prepare(dict(contract, seed=0), reference, cell_types, cell, window)
        if prep is None:
            print("    bank empty, skipped", flush=True)
            continue
        cfg = prep["cfg"]

        # Replay the chain with `choose_band` only -- cheap, and bit-identical to the loop_axis run.
        bands, covered = {}, set()
        for loop in range(1, LOOPS + 1):
            cands = [s for s in window if s not in covered]
            if len(cands) < cfg.band_k:
                break
            band, alive = CJ.choose_band(prep["ok"], cfg.band_k, cands, cfg.group_size,
                                         prep["cell_ok"], {}, cfg.band_quota)
            if len(band) < cfg.band_k or len(alive) < cfg.group_size:
                break
            bands[loop] = sorted(int(x) for x in band)
            covered |= set(bands[loop])

        for loop in loops:
            if (cell, tid, loop) in done or loop not in bands:
                continue
            band = bands[loop]
            bcfg = dataclasses.replace(cfg, band_candidates=tuple(band))
            t0 = time.monotonic()
            rows, meta = CJ.build_submission(dict(contract, seed=0), reference, cell_types,
                                             cfg=bcfg, budget_s=1800)
            MT.free_gpu_memory()
            if not rows:
                print(f"    L{loop:<2} DECLINED — {meta.get('reason')}", flush=True)
                continue
            got = sorted(int(x) for x in meta["band_seeds"])
            clean = set(int(x) for x in meta["clean_seeds"])
            match = got == band
            bh = [s for s in seeds if s in set(got)]
            ch = [s for s in seeds if s in clean and s not in set(got)]
            dirty = [s for s in seeds if s not in clean and s not in set(got)]
            vc = V_CLEAN.get(cell, V_CLEAN["_"]); vd = V_DIRTY.get(cell, V_DIRTY["_"])
            cons = (len(bh) * 1.0 + len(ch) * vc + len(dirty) * vd) / 3
            print(f"    L{loop:<2} band {len(got)}{'' if match else ' MISMATCH'}  "
                  f"clean {len(clean):>4}/900  |  band hits {len(bh)} {bh}  "
                  f"clean hits {len(ch)} {ch}  dirty {len(dirty)}  ->  cons {cons:.4f}  "
                  f"({time.monotonic()-t0:.0f}s)", flush=True)
            out.append({"cell": cell, "task": tid, "loop": loop, "seeds": seeds,
                        "band": got, "band_match": match, "clean_n": len(clean),
                        "band_hits": bh, "clean_hits": ch, "n_dirty": len(dirty),
                        "cons": cons, "wxfid": meta.get("weighted_x_fidelity"),
                        "window_size": len(window)})
            json.dump(out, open(OUT, "w"), indent=1)

    if out:
        print(f"\n=== summary, {len(out)} hit-loops ===", flush=True)
        for cell in CELLS:
            cr = [x for x in out if x["cell"] == cell]
            if not cr:
                continue
            print(f"  {cell:<12} clean {np.mean([x['clean_n'] for x in cr]):>6.1f}/900  "
                  f"clean hits {sum(len(x['clean_hits']) for x in cr)}/{2*len(cr)} off-band seeds  "
                  f"cons {np.mean([x['cons'] for x in cr]):.4f}", flush=True)
    print(f"\nwrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
