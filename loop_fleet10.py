#!/usr/bin/env python3
"""loop_fleet10.py — TEN virtual hotkeys on rotated width-270 windows, each running its OWN 10 loops.

The arm under test, as specified:

    joined window / cut   900 (wide) for every cell
    sub window            270 seeds
    stride                90  ->  10 x 90 = 900, a CIRCULAR tiling at 3x overlap (270/90), wrapping
    loops                 10 per hotkey, each hotkey chained SEPARATELY
    everything else       each cell's live CELL_CONFIG (k, group_size, light_cell_rows, cell-aware)

**This is the cross product, not the diagonal.** The live fleet gives hotkey i loop i+1 (one band
each, ten disjoint bands). Here every one of the ten hotkeys runs its own full ten-loop chain, so
there are 100 (hotkey, loop) band draws per task -- which is what a "hotkey rows x loop columns"
table needs. Bands are disjoint WITHIN a hotkey's chain by construction and free to collide ACROSS
hotkeys, because the windows overlap 3x.

**Why this is affordable.** One cut space (900) means ONE bank per (contract, cell) shared by all
ten. `choose_band` restricts both `ok` and `cell_ok` to its `candidates` by column lookup, so one
`hdr_compliance` + `cas9_cell_probe` pass over the whole 900 is valid for every sub-window of it --
the same trick `width30_test.py` used. The per-hotkey slice is then taken from that one pass, one
hotkey at a time so only a single slice is resident (pre-slicing all ten as Python sets of ints is
what OOM-killed three runs on 2026-09-20; these are numpy int32).

The clean-set half needs real rows -- `meta["clean_seeds"]` does not exist until `build_submission`
runs -- so one full build is done at each (hotkey, loop) that CAUGHT a seed, with the band pinned by
narrowing `band_candidates` to exactly that loop's k seeds. The greedy provably takes all k (the
chain already proved the intersection clears `group_size`, and intersections only shrink), and the
resulting set is asserted equal to the chain's band rather than assumed.

    LF_N=10 LF_LOOPS=10 LF_WIDTH=270 LF_STRIDE=90 python loop_fleet10.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "loopfleet10")

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
OUT = os.getenv("LF_JSON", "loop_fleet10.json")
MIN_FREE_GB = float(os.getenv("LF_MIN_FREE_GB", "7"))
SPAN = len(JW.FULL_SPACE)
assert NHK * STRIDE == SPAN, f"{NHK} x {STRIDE} != {SPAN}: the circular tiling is not exact"

V_CLEAN = {"HEK293": 0.212}
V_DIRTY = {"HEK293": 0.082}


def slice_ok(ok, window):
    """`ok` restricted to one hotkey's window, numpy-backed.

    Not cosmetic: `choose_band` walks every compliant seed of every record through a dict lookup,
    so handing it compliance over all 900 when only 270 are candidates triples that walk on every
    one of the 100 calls. Sliced once per hotkey and discarded before the next.
    """
    w = np.fromiter(sorted(window), dtype=np.int32)
    out = {}
    for i, v in ok.items():
        keep = v[np.isin(v, w)]
        if keep.size:
            out[i] = keep
    return out


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    items = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])
    tasks = []
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
    cut = list(JW.FULL_SPACE)
    windows = [CJ.sub_window(cut, WIDTH, (i * STRIDE) / float(SPAN)) for i in range(NHK)]
    print(f"=== {NHK} test hotkeys | sub window {WIDTH} stride {STRIDE} (circular, "
          f"{WIDTH//STRIDE}x overlap) | cut {len(cut)} wide | {LOOPS} loops EACH | "
          f"{len(tasks)} tasks ===", flush=True)
    for i, w in enumerate(windows):
        print(f"    h{i}  {min(w)}-{max(w)}  {len(w)} seeds"
              + ("  (wraps)" if max(w) - min(w) + 1 != len(w) else ""), flush=True)
    print(flush=True)

    for t in tasks:
        tid = (t.get("task_id") or t["id"])[:8]
        if tid in done:
            continue
        # The ten live miners hold ~3.3 GB each at idle, so this run has ~9 GB of headroom and the
        # first attempt died at 10.97 GB. Wait rather than race them into another OOM kill.
        for _ in range(60):
            avail = int([l for l in open("/proc/meminfo") if l.startswith("MemAvailable")][0]
                        .split()[1]) / 1048576.0
            if avail >= MIN_FREE_GB:
                break
            print(f"    waiting for memory: {avail:.1f} GB available < {MIN_FREE_GB} GB",
                  flush=True)
            time.sleep(60)
        contract, reference = dict(t["content"]["contract"]), t["content"]["hbb_reference"]
        cell = contract["cell_type"]
        seeds = _parse_seeds(contract["seed"])
        base = CJ.config_for(cell)
        cfg = dataclasses.replace(base, seed_list=tuple(cut),
                                  cas12a_max_fail=max(1, round(base.cas12a_max_fail
                                                               * len(cut) / 900)),
                                  start_seed=cut[0], end_seed=cut[-1])
        t0 = time.monotonic()
        ctx = G.build_context(dict(contract, seed=0), reference, cell_types)
        sites = G.enumerate_sites(ctx, 3000, (20, 23))
        path = os.path.join(CJ.BANK_DIR,
                            f"cas12a-{CJ.bank_key(dict(contract, seed=0), cell_types, cfg)}.npz")
        if not os.path.exists(path):
            os.makedirs(CJ.BANK_DIR, exist_ok=True)
            bank = CJ.build_bank(dict(contract, seed=0), reference, cell_types, ctx, sites,
                                 cfg, None)
            MT.free_gpu_memory()
            if not bank:
                print(f"  {tid} {cell} bank empty, skipped", flush=True)
                continue
            CJ.save_bank(path, bank)
        records = CJ.load_bank(path, limit=cfg.bank_keep)
        n_bank = len(records)
        ok_all = CJ.hdr_compliance(records, dict(contract, seed=0), cell_types, ctx, cut,
                                   cfg.band_rule)
        MT.free_gpu_memory()
        ok_all = {i: np.fromiter(v, dtype=np.int32) for i, v in ok_all.items()}
        # `records` is the single biggest object here -- 300,000 banked guides with their `fails`
        # arrays -- and NOTHING below needs it: the chains run off `ok_all`, `cas9_cell_probe` takes
        # `sites`, and `build_submission` re-loads the bank itself. Holding it through the build
        # phase meant two copies resident at once, which is what OOM-killed the first run.
        del records
        gc.collect()
        cell_all = None
        if cfg.band_cell_aware:
            cell_all = CJ.cas9_cell_probe(dict(contract, seed=0), cell_types, ctx, sites, cfg,
                                          cut, cfg.band_rule)
            MT.free_gpu_memory()
        print(f"  {tid} {cell:<11} seeds {seeds}  k={cfg.band_k} group={cfg.group_size}  "
              f"bank {n_bank}  prep {time.monotonic()-t0:.0f}s", flush=True)

        hks = []
        for i, window in enumerate(windows):
            ok = slice_ok(ok_all, window)
            cok = None
            if cell_all is not None:
                ws = set(window)
                cok = {k: [{s for s in ss if s in ws} for ss in v] for k, v in cell_all.items()}
            bands, covered, found = {}, set(), {}
            for loop in range(1, LOOPS + 1):
                cands = [s for s in window if s not in covered]
                if len(cands) < cfg.band_k:
                    break
                band, alive = CJ.choose_band(ok, cfg.band_k, cands, cfg.group_size, cok, {},
                                             cfg.band_quota)
                if len(band) < cfg.band_k or len(alive) < cfg.group_size:
                    break
                bands[loop] = sorted(int(x) for x in band)
                covered |= set(bands[loop])
                for s in seeds:
                    if s in covered and s not in found:
                        found[s] = loop
            inwin = [s for s in seeds if s in set(window)]
            hks.append({"hk": i, "window": [min(window), max(window)], "n_window": len(window),
                        "in_window": inwin, "bands": bands, "covered": len(covered),
                        "found": {str(s): found.get(s) for s in seeds}, "builds": []})
            print(f"    h{i} win {min(window)}-{max(window)} inwin {inwin} "
                  f"cov {len(covered)}/{len(window)} loops {len(bands)}  "
                  f"found {({s: found[s] for s in seeds if s in found}) or '{}'}", flush=True)
            del ok, cok

        # **Free the compliance structures BEFORE the builds.** `build_submission` re-loads the
        # bank itself, so `records` / `ok_all` / `cell_all` are dead weight from here on -- and at a
        # 300,000-record erythroid bank over 900 seeds they are the bulk of the process. Holding
        # them through the build phase is what OOM-killed the first run of this script at 10.97 GB
        # with the ten live miners already holding ~27 GB of the box's 49.
        del ok_all, cell_all
        ok_all = cell_all = None
        gc.collect()

        # Clean sets, only where a band caught a seed.
        for h in hks:
            for loop in sorted({l for l in h["found"].values() if l}):
                band = h["bands"][loop]
                bcfg = dataclasses.replace(cfg, band_candidates=tuple(band))
                b0 = time.monotonic()
                rows, meta = CJ.build_submission(dict(contract, seed=0), reference, cell_types,
                                                 cfg=bcfg, budget_s=1800)
                MT.free_gpu_memory()
                if not rows:
                    print(f"      h{h['hk']} L{loop} DECLINED — {meta.get('reason')}", flush=True)
                    continue
                got = sorted(int(x) for x in meta["band_seeds"])
                clean = set(int(x) for x in meta["clean_seeds"])
                bh = [s for s in seeds if s in set(got)]
                ch = [s for s in seeds if s in clean and s not in set(got)]
                dirty = [s for s in seeds if s not in clean and s not in set(got)]
                vc = V_CLEAN.get(cell, 0.212); vd = V_DIRTY.get(cell, 0.104)
                cons = (len(bh) * 1.0 + len(ch) * vc + len(dirty) * vd) / 3
                print(f"      h{h['hk']} L{loop:<2} band {len(got)}"
                      f"{'' if got == band else ' MISMATCH'}  clean {len(clean):>4}/900 "
                      f"({len(clean)/900:>5.1%})  band hits {len(bh)} {bh}  "
                      f"clean hits {len(ch)} {ch}  -> cons {cons:.4f} "
                      f"({time.monotonic()-b0:.0f}s)", flush=True)
                h["builds"].append({"loop": loop, "band": got, "band_match": got == band,
                                    "clean_n": len(clean), "clean_seeds": sorted(clean),
                                    "band_hits": bh, "clean_hits": ch, "n_dirty": len(dirty),
                                    "cons": cons})
        out.append({"task": tid, "cell": cell, "seeds": seeds, "k": cfg.band_k,
                    "group": cfg.group_size, "width": WIDTH, "stride": STRIDE,
                    "cut": len(cut), "loops": LOOPS, "hotkeys": hks})
        json.dump(out, open(OUT, "w"), indent=1)
        gc.collect()
    print(f"\nwrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
