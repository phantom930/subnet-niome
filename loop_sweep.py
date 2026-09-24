#!/usr/bin/env python3
"""loop_sweep.py — the loop axis across CUT SPAN x BAND DEPTH, to pick the best arm per cell.

Six arms per cell, the full cross product of the two knobs under test:

    cut    wide (900)  |  narrow (the task's own 300-seed oracle window)
    band_k HEK293 6/8/9      erythroid 10/11/12

11 loops per arm per task, chains run to completion so the loop at which EACH of the three seeds
is first reached is recorded separately (`found` maps seed -> loop, or None).

**Window: the REAL three seed classes, not `baseline_uniform`.** The width-100 classes the round's
seeds actually drew from, unioned. When two or three seeds share a class the window would be 200 or
100 seeds wide and would cover far faster per loop for a reason that has nothing to do with the
arm, so it is PADDED to three classes with classes the seeds did NOT use -- adding window size
without adding reachable seeds. Deterministic per task, so a re-run reproduces the same window.

**What makes six arms affordable: `bank_key` excludes every band parameter.** So the three k values
of one cut span share ONE bank and ONE `hdr_compliance` pass -- two banks per task total, not six.
`choose_band` restricts to its `candidates` by column lookup, so compliance computed over the
300-seed window serves both spans (the band candidates are the window in BOTH arms; only the CUT
differs).

**The clean-set half is deliberately run on FEWER tasks than the band half.** `meta["clean_seeds"]`
needs a real `build_submission`, ~60-130s each, and at 6 arms x 50 tasks x 4 cells that is days.
The band half (chains only) is cheap and runs on all `LS_N` tasks; builds run on the newest
`LS_BUILD_N`. Clean density is a property of (arm, task) that barely moves across loops within a
task -- measured 97-103 of 270 across ten loops on 2026-09-21 -- so a smaller n estimates it well,
while the loop-count question genuinely wants the larger n.

**Cross-arm clean counts are NOT comparable and the script records what makes them so.** A narrow
arm's clean set is counted within 300 seeds and a wide arm's within 900, and clean is a subset of
the cut -- so in the narrow arm a seed outside the window can never be a clean hit. `n_offband_in_cut`
is stored per build; compute hit rates on that, never on all off-band seeds.

    LS_CELLS=HEK293 LS_N=50 LS_LOOPS=11 LS_BUILD_N=10 python loop_sweep.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "loopsweep")

import dataclasses, gc, json, logging, random, time                    # noqa: E402

logging.basicConfig(level=logging.ERROR)
import numpy as np                                                    # noqa: E402

import genExp as G                                                    # noqa: E402
import joined_window as JW                                            # noqa: E402
from niome_subnet.genomics import conjunction as CJ                   # noqa: E402
from niome_subnet.genomics import mt19937 as MT                       # noqa: E402
from niome_subnet.genomics.validation import _parse_seeds             # noqa: E402
from conj_stageb import API                                           # noqa: E402
from sd_task import fetch                                             # noqa: E402

CELLS = [c for c in os.getenv("LS_CELLS", "HEK293,K562,HUDEP-2,CD34+_HSPC").split(",") if c]
N = int(os.getenv("LS_N", "50"))
LOOPS = int(os.getenv("LS_LOOPS", "11"))
BUILD_N = int(os.getenv("LS_BUILD_N", "10"))
OUT = os.getenv("LS_JSON", "loop_sweep.json")
MIN_FREE_GB = float(os.getenv("LS_MIN_FREE_GB", "7"))
KS = {"HEK293": (6, 8, 9), "K562": (10, 11, 12),
      "HUDEP-2": (10, 11, 12), "CD34+_HSPC": (10, 11, 12)}
V_CLEAN = {"HEK293": 0.212}
V_DIRTY = {"HEK293": 0.082}


def wait_mem():
    for _ in range(90):
        a = int([l for l in open("/proc/meminfo") if l.startswith("MemAvailable")][0]
                .split()[1]) / 1048576.0
        if a >= MIN_FREE_GB:
            return
        print(f"    waiting for memory: {a:.1f} GB < {MIN_FREE_GB}", flush=True)
        time.sleep(60)


def oracle_window(tid, seeds):
    """The real width-100 classes the seeds drew from, padded to three."""
    real = sorted({(s - 100) // 100 for s in seeds})
    rng = random.Random(int(tid[:8], 16))
    pad = [w for w in range(9) if w not in real]
    rng.shuffle(pad)
    classes = sorted(real + pad[:max(0, 3 - len(real))])
    win = sorted(s for w in classes for s in range(w * 100 + 100, w * 100 + 200))
    return classes, real, sorted(set(classes) - set(real)), win


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    items = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])
    out = json.load(open(OUT)) if os.path.exists(OUT) else []
    done = {(r["cell"], r["task"]) for r in out}

    for cell in CELLS:
        base = CJ.config_for(cell)
        ks = KS[cell]
        tasks = []
        for t in sorted(items, key=lambda x: x.get("created_at", ""), reverse=True):
            c = (t.get("content") or {}).get("contract") or {}
            if c.get("cell_type") != cell or len(_parse_seeds(c.get("seed"))) != 3:
                continue
            tasks.append(t)
            if len(tasks) >= N:
                break
        build_ids = {(t.get("task_id") or t["id"])[:8] for t in tasks[:BUILD_N]}
        print(f"\n=== {cell} | k {ks} x cut {{wide 900, narrow 300}} | {LOOPS} loops | "
              f"{len(tasks)} tasks, builds on newest {len(build_ids)} ===", flush=True)

        for ti, t in enumerate(tasks, 1):
            tid = (t.get("task_id") or t["id"])[:8]
            if (cell, tid) in done:
                continue
            wait_mem()
            contract, reference = dict(t["content"]["contract"]), t["content"]["hbb_reference"]
            seeds = _parse_seeds(contract["seed"])
            classes, real, padded, window = oracle_window(tid, seeds)
            t0 = time.monotonic()
            ctx = G.build_context(dict(contract, seed=0), reference, cell_types)
            sites = G.enumerate_sites(ctx, 3000, (20, 23))
            rec = {"cell": cell, "task": tid, "seeds": seeds, "classes": classes,
                   "real_classes": real, "padded_classes": padded, "loops": LOOPS, "arms": {}}

            for span, cut in (("wide", list(JW.FULL_SPACE)), ("narrow", sorted(window))):
                mf = max(1, round(base.cas12a_max_fail * len(cut) / 900))
                cfg0 = dataclasses.replace(base, seed_list=tuple(sorted(cut)),
                                           start_seed=min(cut), end_seed=max(cut),
                                           cas12a_max_fail=mf,
                                           band_candidates=tuple(sorted(window)))
                path = os.path.join(CJ.BANK_DIR,
                                    f"cas12a-{CJ.bank_key(dict(contract, seed=0), cell_types, cfg0)}.npz")
                cold = not os.path.exists(path)
                if cold:
                    os.makedirs(CJ.BANK_DIR, exist_ok=True)
                    bank = CJ.build_bank(dict(contract, seed=0), reference, cell_types, ctx,
                                         sites, cfg0, None)
                    MT.free_gpu_memory()
                    if not bank:
                        print(f"  {tid} {span} bank empty, skipped", flush=True); continue
                    CJ.save_bank(path, bank); del bank; gc.collect()
                records = CJ.load_bank(path, limit=cfg0.bank_keep)
                n_bank = len(records)
                ok = CJ.hdr_compliance(records, dict(contract, seed=0), cell_types, ctx,
                                       window, cfg0.band_rule)
                MT.free_gpu_memory()
                ok = {i: np.fromiter(v, dtype=np.int32) for i, v in ok.items()}
                del records; gc.collect()
                cell_ok = None
                if cfg0.band_cell_aware:
                    cell_ok = CJ.cas9_cell_probe(dict(contract, seed=0), cell_types, ctx, sites,
                                                 cfg0, window, cfg0.band_rule)
                    MT.free_gpu_memory()

                for k in ks:
                    cfg = dataclasses.replace(cfg0, band_k=k)
                    bands, covered, found = {}, set(), {}
                    for loop in range(1, LOOPS + 1):
                        cands = [s for s in window if s not in covered]
                        if len(cands) < k:
                            break
                        band, alive = CJ.choose_band(ok, k, cands, cfg.group_size, cell_ok, {},
                                                     cfg.band_quota)
                        if len(band) < k or len(alive) < cfg.group_size:
                            break
                        bands[loop] = sorted(int(x) for x in band)
                        covered |= set(bands[loop])
                        for s in seeds:
                            if s in covered and s not in found:
                                found[s] = loop
                    rec["arms"][f"{span}/k{k}"] = {
                        "span": span, "cut": len(cut), "k": k, "bank": n_bank, "cold": cold,
                        "loops_built": len(bands), "covered": len(covered),
                        "found": {str(s): found.get(s) for s in seeds},
                        "bands": bands, "builds": []}
                del ok, cell_ok; gc.collect()

                # clean sets: newest tasks only, one build per arm at its FIRST found loop
                if tid in build_ids:
                    for k in ks:
                        a = rec["arms"].get(f"{span}/k{k}")
                        if not a: continue
                        hits = sorted({l for l in a["found"].values() if l})
                        if not hits: continue
                        loop = hits[0]
                        bcfg = dataclasses.replace(cfg0, band_k=k,
                                                   band_candidates=tuple(a["bands"][loop]))
                        rows, meta = CJ.build_submission(dict(contract, seed=0), reference,
                                                         cell_types, cfg=bcfg, budget_s=1800)
                        MT.free_gpu_memory()
                        if not rows:
                            a["builds"].append({"loop": loop, "declined": meta.get("reason")})
                            continue
                        got = sorted(int(x) for x in meta["band_seeds"])
                        clean = set(int(x) for x in meta["clean_seeds"])
                        bh = [s for s in seeds if s in set(got)]
                        ch = [s for s in seeds if s in clean and s not in set(got)]
                        off = [s for s in seeds if s not in set(got)]
                        oic = [s for s in off if s in set(cut)]
                        vc = V_CLEAN.get(cell, 0.212); vd = V_DIRTY.get(cell, 0.104)
                        cons = (len(bh) + len(ch) * vc + (len(off) - len(ch)) * vd) / 3
                        a["builds"].append({
                            "loop": loop, "band": got, "band_match": got == a["bands"][loop],
                            "clean_n": len(clean), "clean_of": len(cut),
                            "band_hits": bh, "clean_hits": ch, "n_offband": len(off),
                            "n_offband_in_cut": len(oic), "cons": cons,
                            "pool": meta.get("pool"), "cas9": meta.get("cas9_pool")})
                gc.collect()

            line = "  ".join(
                f"{a}:{sum(1 for v in rec['arms'][a]['found'].values() if v)}/3"
                f"@{sorted(v for v in rec['arms'][a]['found'].values() if v) or '-'}"
                for a in sorted(rec["arms"]))
            print(f"  [{ti}/{len(tasks)}] {tid} seeds {seeds} cls {classes}"
                  f"{'+pad' + str(len(padded)) if padded else ''}  ({time.monotonic()-t0:.0f}s)\n"
                  f"      {line}", flush=True)
            out.append(rec)
            json.dump(out, open(OUT, "w"), indent=1)
            gc.collect()
    print(f"\nwrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
