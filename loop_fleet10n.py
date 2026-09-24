#!/usr/bin/env python3
"""loop_fleet10n.py — the NARROW-cut arm of `loop_fleet10.py`: cut == each hotkey's OWN sub-window.

Identical to `loop_fleet10.py` in every other respect -- ten virtual hotkeys, sub window 270 at
stride 90 (circular, 3x overlap, wrapping), ten loops EACH, the same latest 10 tasks, each cell's
live CELL_CONFIG -- so the pair isolates the CUT SPAN and nothing else:

    loop_fleet10.py   cut 900, shared by all ten          -> ONE bank per (contract, cell)
    loop_fleet10n.py  cut 270, each hotkey's own window   -> TEN banks per (contract, cell)

**That bank multiplication is the real cost and it is structural.** `bank_key` folds in `seed_list`,
so ten different cut windows are ten different banks, every one of them a cold Cas12a scan that
fills to the 300,000 `bank_keep` cap (a narrow cut is an easier constraint than a wide one, so more
guides qualify and the scan runs to the cap rather than collapsing below it). Nothing can be shared:
`hdr_compliance` is computed against a hotkey's OWN bank, so the "one pass over 900 serves every
sub-window" trick that made the wide arm cheap does not apply here at all. This is the same
arithmetic that made the live fleet's 2026-09-20 per-hotkey cut cost seven concurrent 566-580s
scans; here they are sequential rather than concurrent, so it is slow rather than contended.

**One asymmetry to carry into the comparison, because it is easy to read wrong.** The clean set is
a subset of the CUT space, so in this arm a round seed outside a hotkey's 270-seed window cannot be
a clean hit for that hotkey at any density -- where in the wide arm all three seeds were eligible.
`meta["clean"]` is therefore out of 270 here and out of 900 there, which CLAUDE.md flags as biased
across arms. The quantity that IS comparable is the hit rate among off-band seeds THAT LIE IN THE
CUT SPACE, so `n_offband_in_cut` is recorded per build alongside the raw count.

    LF_N=10 LF_LOOPS=10 LF_WIDTH=270 LF_STRIDE=90 python loop_fleet10n.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "loopfleet10n")

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

N = int(os.getenv("LF_N", "10"))
LOOPS = int(os.getenv("LF_LOOPS", "10"))
WIDTH = int(os.getenv("LF_WIDTH", "270"))
STRIDE = int(os.getenv("LF_STRIDE", "90"))
NHK = int(os.getenv("LF_HK", "10"))
OUT = os.getenv("LF_JSON", "loop_fleet10n.json")
MIN_FREE_GB = float(os.getenv("LF_MIN_FREE_GB", "7"))
SPAN = len(JW.FULL_SPACE)
assert NHK * STRIDE == SPAN, f"{NHK} x {STRIDE} != {SPAN}: the circular tiling is not exact"

V_CLEAN = {"HEK293": 0.212}
V_DIRTY = {"HEK293": 0.082}


def wait_for_memory():
    for _ in range(60):
        avail = int([l for l in open("/proc/meminfo") if l.startswith("MemAvailable")][0]
                    .split()[1]) / 1048576.0
        if avail >= MIN_FREE_GB:
            return
        print(f"    waiting for memory: {avail:.1f} GB < {MIN_FREE_GB} GB", flush=True)
        time.sleep(60)


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    items = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])
    # **Pinned to the WIDE run's task list, not re-derived as "the latest 10".** The task feed moves
    # -- two new rounds landed between the two arms -- so re-deriving would silently swap tasks in
    # and out and leave the arms unpaired on the one thing the comparison is about. `LF_PAIR` names
    # the wide run's output; set it empty to take the latest 10 instead.
    pair = os.getenv("LF_PAIR", "loop_fleet10.json")
    want = None
    if pair and os.path.exists(pair):
        want = [r["task"] for r in json.load(open(pair))]
        print(f"    task list PINNED to {pair}: {want}", flush=True)
    tasks = []
    if want:
        by = {(t.get("task_id") or t["id"])[:8]: t for t in items}
        tasks = [by[x] for x in want if x in by]
        missing = [x for x in want if x not in by]
        if missing:
            print(f"    WARNING {len(missing)} pinned task(s) not in the feed: {missing}",
                  flush=True)
    else:
        for t in sorted(items, key=lambda x: x.get("created_at", ""), reverse=True):
            c = (t.get("content") or {}).get("contract") or {}
            if not c.get("cell_type") or len(_parse_seeds(c.get("seed"))) != 3:
                continue
            if CJ.config_for(c["cell_type"]) is None:
                continue
            tasks.append(t)
            if len(tasks) >= N:
                break

    out = json.load(open(OUT)) if os.path.exists(OUT) else []
    done = {r["task"] for r in out}
    full = list(JW.FULL_SPACE)
    windows = [CJ.sub_window(full, WIDTH, (i * STRIDE) / float(SPAN)) for i in range(NHK)]
    print(f"=== {NHK} test hotkeys | sub window {WIDTH} stride {STRIDE} | "
          f"CUT = own window ({WIDTH}) | {LOOPS} loops EACH | {len(tasks)} tasks ===", flush=True)
    print(f"    ten distinct cut spans -> TEN banks per (contract, cell); nothing is shared\n",
          flush=True)

    for t in tasks:
        tid = (t.get("task_id") or t["id"])[:8]
        if tid in done:
            continue
        wait_for_memory()
        contract, reference = dict(t["content"]["contract"]), t["content"]["hbb_reference"]
        cell = contract["cell_type"]
        seeds = _parse_seeds(contract["seed"])
        base = CJ.config_for(cell)
        ctx = G.build_context(dict(contract, seed=0), reference, cell_types)
        sites = G.enumerate_sites(ctx, 3000, (20, 23))
        print(f"  {tid} {cell:<11} seeds {seeds}  k={base.band_k} group={base.group_size}",
              flush=True)

        hks = []
        for i, window in enumerate(windows):
            w0 = time.monotonic()
            wset = set(window)
            # `cas12a_max_fail` is calibrated at the 900-seed span, so it scales with the cut --
            # exactly what `build_for_cell` does. Without it a 270-seed cut silently demands a far
            # stricter guide than the arm was measured at.
            cfg = dataclasses.replace(
                base, seed_list=tuple(sorted(window)), band_candidates=tuple(sorted(window)),
                start_seed=min(window), end_seed=max(window),
                cas12a_max_fail=max(1, round(base.cas12a_max_fail * len(window) / 900)))
            path = os.path.join(
                CJ.BANK_DIR, f"cas12a-{CJ.bank_key(dict(contract, seed=0), cell_types, cfg)}.npz")
            if not os.path.exists(path):
                os.makedirs(CJ.BANK_DIR, exist_ok=True)
                bank = CJ.build_bank(dict(contract, seed=0), reference, cell_types, ctx, sites,
                                     cfg, None)
                MT.free_gpu_memory()
                if not bank:
                    print(f"    h{i} bank empty, skipped", flush=True)
                    continue
                CJ.save_bank(path, bank)
                del bank
                gc.collect()
            records = CJ.load_bank(path, limit=cfg.bank_keep)
            n_bank = len(records)
            ok = CJ.hdr_compliance(records, dict(contract, seed=0), cell_types, ctx, window,
                                   cfg.band_rule)
            MT.free_gpu_memory()
            ok = {j: np.fromiter(v, dtype=np.int32) for j, v in ok.items()}
            del records
            gc.collect()
            cell_ok = None
            if cfg.band_cell_aware:
                cell_ok = CJ.cas9_cell_probe(dict(contract, seed=0), cell_types, ctx, sites, cfg,
                                             window, cfg.band_rule)
                MT.free_gpu_memory()

            bands, covered, found = {}, set(), {}
            for loop in range(1, LOOPS + 1):
                cands = [s for s in window if s not in covered]
                if len(cands) < cfg.band_k:
                    break
                band, alive = CJ.choose_band(ok, cfg.band_k, cands, cfg.group_size, cell_ok, {},
                                             cfg.band_quota)
                if len(band) < cfg.band_k or len(alive) < cfg.group_size:
                    break
                bands[loop] = sorted(int(x) for x in band)
                covered |= set(bands[loop])
                for s in seeds:
                    if s in covered and s not in found:
                        found[s] = loop
            del ok, cell_ok
            gc.collect()

            inwin = [s for s in seeds if s in wset]
            print(f"    h{i} win {min(window)}-{max(window)} bank {n_bank} inwin {inwin} "
                  f"cov {len(covered)}/{len(window)} loops {len(bands)} "
                  f"found {({s: found[s] for s in seeds if s in found}) or '{}'} "
                  f"({time.monotonic()-w0:.0f}s)", flush=True)

            builds = []
            for loop in sorted({l for l in found.values() if l}):
                bcfg = dataclasses.replace(cfg, band_candidates=tuple(bands[loop]))
                b0 = time.monotonic()
                rows, meta = CJ.build_submission(dict(contract, seed=0), reference, cell_types,
                                                 cfg=bcfg, budget_s=1800)
                MT.free_gpu_memory()
                if not rows:
                    print(f"      h{i} L{loop} DECLINED — {meta.get('reason')}", flush=True)
                    continue
                got = sorted(int(x) for x in meta["band_seeds"])
                clean = set(int(x) for x in meta["clean_seeds"])
                bh = [s for s in seeds if s in set(got)]
                ch = [s for s in seeds if s in clean and s not in set(got)]
                off = [s for s in seeds if s not in set(got)]
                # Only an off-band seed INSIDE the cut can ever be clean here; the wide arm had no
                # such restriction. This is the denominator the two arms must be compared on.
                off_in_cut = [s for s in off if s in wset]
                vc = V_CLEAN.get(cell, 0.212); vd = V_DIRTY.get(cell, 0.104)
                cons = (len(bh) * 1.0 + len(ch) * vc + (len(off) - len(ch)) * vd) / 3
                print(f"      h{i} L{loop:<2} band {len(got)}"
                      f"{'' if got == bands[loop] else ' MISMATCH'}  clean {len(clean):>4}/"
                      f"{len(window)} ({len(clean)/len(window):>5.1%})  band hits {len(bh)} {bh}  "
                      f"clean hits {len(ch)} {ch}  off-in-cut {len(off_in_cut)}  "
                      f"-> cons {cons:.4f} ({time.monotonic()-b0:.0f}s)", flush=True)
                builds.append({"loop": loop, "band": got, "band_match": got == bands[loop],
                               "clean_n": len(clean), "clean_frac": len(clean) / len(window),
                               "clean_seeds": sorted(clean), "band_hits": bh, "clean_hits": ch,
                               "n_offband": len(off), "n_offband_in_cut": len(off_in_cut),
                               "cons": cons})
            hks.append({"hk": i, "window": [min(window), max(window)], "n_window": len(window),
                        "bank": n_bank, "in_window": inwin, "bands": bands,
                        "covered": len(covered),
                        "found": {str(s): found.get(s) for s in seeds}, "builds": builds})
            gc.collect()

        out.append({"task": tid, "cell": cell, "seeds": seeds, "k": base.band_k,
                    "group": base.group_size, "width": WIDTH, "stride": STRIDE,
                    "cut": WIDTH, "arm": "narrow", "loops": LOOPS, "hotkeys": hks})
        json.dump(out, open(OUT, "w"), indent=1)
        gc.collect()
    print(f"\nwrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
