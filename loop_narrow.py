#!/usr/bin/env python3
"""loop_narrow.py — the loop axis at a NARROW cut: joined window == sub window == the same 300 seeds.

`loop_axis.py` + `loop_clean.py` ran the loop axis at cut 900 (wide) and measured the floor as
worthless: 4 clean hits across 84 in-window off-band seeds (4.8%), against a clean set that is
15.2% of 900 globally. That gap is the point of this run. Those 84 seeds are direct membership
tests INSIDE the band window, so the wide arm's clean set looks depleted exactly where the round's
seeds sit -- plausible, because a min-union over 900 seeds has no reason to spend its freedom on
the 300 that matter.

Here the cut is the SAME 300 seeds as the band window, so:

  * `band` is a subset of `cut` by construction (the `narrow`/`union` property, CLAUDE.md);
  * the min-union optimises cut-cleanness INSIDE the window rather than across all 900;
  * `meta["clean"]` is the within-window count directly, so the density needs no inference.

**The bands change too, and that is why this is a full re-run.** `bank_key` folds in `seed_list`,
so a different cut is a different bank, different surviving records, different `choose_band` output.
The loop at which each seed is reached is therefore NOT carried over from the wide run -- it is
re-derived here, and the two arms are compared on rates, not on matched loops.

One pass per task: prepare (bank + compliance + probe) once, run the whole `choose_band` chain, then
one full `build_submission` at each loop that caught a seed, with the band pinned by narrowing
`band_candidates` to exactly that loop's k seeds (the greedy provably takes all k; asserted).

**Comparing clean counts across the two arms needs care and is deliberately not done here.**
CLAUDE.md records that `meta["clean"]` is biased across arms with different cut spans -- 900 seeds
of opportunity against 300. The quantity that IS comparable, and the one the floor question turns
on, is the HIT RATE: what fraction of the round's in-window off-band seeds land in the clean set.
Both arms have all three seeds in-window by construction, so that rate is measured identically.

    LN_CELLS=HUDEP-2,HEK293,K562,CD34+_HSPC LN_LOOPS=10 python loop_narrow.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "loopnarrow")

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

CELLS = [c for c in os.getenv("LN_CELLS", "HUDEP-2,HEK293,K562,CD34+_HSPC").split(",") if c]
N = int(os.getenv("LN_N", "10"))
LOOPS = int(os.getenv("LN_LOOPS", "10"))
OUT = os.getenv("LN_JSON", "loop_narrow.json")

V_CLEAN = {"HEK293": 0.212}
V_DIRTY = {"HEK293": 0.082}


def prepare(contract, reference, cell_types, cell, window):
    """Same as `loop_axis.prepare` except the cut is the WINDOW, not the full 900."""
    base = CJ.config_for(cell)
    cut = sorted(window)                       # <-- the whole point of this script
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
    out = json.load(open(OUT)) if os.path.exists(OUT) else []
    done = {(r["cell"], r["task"]) for r in out}
    print(f"=== loop axis, NARROW cut (cut == window == 300) | {LOOPS} loops | cells {CELLS} ===",
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
        print(f"\n--- {cell}  k={base.band_k} group={base.group_size} "
              f"light={base.light_cell_rows} | {len(tasks)} tasks ---", flush=True)

        for t in tasks:
            tid = (t.get("task_id") or t["id"])
            if (cell, tid[:8]) in done:
                continue
            contract, reference = dict(t["content"]["contract"]), t["content"]["hbb_reference"]
            seeds = _parse_seeds(contract["seed"])
            # Identical oracle window to the wide run: same per-task RNG, same padding rule, so the
            # two arms differ in the CUT alone.
            real = sorted({(s - 100) // 100 for s in seeds})
            rng = random.Random(int(tid[:8], 16))
            pad = [w for w in range(9) if w not in real]
            rng.shuffle(pad)
            classes = sorted(real + pad[:max(0, 3 - len(real))])
            padded = sorted(set(classes) - set(real))
            window = sorted(s for w in classes for s in range(w * 100 + 100, w * 100 + 200))

            t0 = time.monotonic()
            prep = prepare(dict(contract, seed=0), reference, cell_types, cell, window)
            if prep is None:
                print(f"  {tid[:8]} bank empty, skipped", flush=True)
                continue
            cfg = prep["cfg"]
            print(f"  {tid[:8]} seeds {seeds}  win {len(window)}sd  bank {prep['bank']}  "
                  f"prep {time.monotonic()-t0:.0f}s", flush=True)

            bands, covered, found = {}, set(), {}
            for loop in range(1, LOOPS + 1):
                cands = [s for s in window if s not in covered]
                if len(cands) < cfg.band_k:
                    break
                band, alive = CJ.choose_band(prep["ok"], cfg.band_k, cands, cfg.group_size,
                                             prep["cell_ok"], {}, cfg.band_quota)
                if len(band) < cfg.band_k or len(alive) < cfg.group_size:
                    print(f"    loop {loop} STOPPED band {len(band)}/{cfg.band_k} "
                          f"alive {len(alive)}/{cfg.group_size}", flush=True)
                    break
                bands[loop] = sorted(int(x) for x in band)
                covered |= set(bands[loop])
                for s in seeds:
                    if s in covered and s not in found:
                        found[s] = loop

            v = sorted(found.values())
            two = v[1] if len(v) >= 2 else None
            print(f"    found {({s: found[s] for s in seeds if s in found}) or '{}'}  "
                  f"2-seed L{two if two else '-'}  covered {len(covered)}/{len(window)}",
                  flush=True)

            builds = []
            for loop in sorted(set(found.values())):
                band = bands[loop]
                bcfg = dataclasses.replace(cfg, band_candidates=tuple(band))
                b0 = time.monotonic()
                rows, meta = CJ.build_submission(dict(contract, seed=0), reference, cell_types,
                                                 cfg=bcfg, budget_s=1800)
                MT.free_gpu_memory()
                if not rows:
                    print(f"    L{loop:<2} DECLINED — {meta.get('reason')}", flush=True)
                    continue
                got = sorted(int(x) for x in meta["band_seeds"])
                clean = set(int(x) for x in meta["clean_seeds"])
                bh = [s for s in seeds if s in set(got)]
                ch = [s for s in seeds if s in clean and s not in set(got)]
                dirty = [s for s in seeds if s not in clean and s not in set(got)]
                vc = V_CLEAN.get(cell, 0.212); vd = V_DIRTY.get(cell, 0.104)
                cons = (len(bh) * 1.0 + len(ch) * vc + len(dirty) * vd) / 3
                print(f"    L{loop:<2} band {len(got)}{'' if got == band else ' MISMATCH'}  "
                      f"clean {len(clean):>4}/{len(window)} ({len(clean)/len(window):>5.1%})  |  "
                      f"band hits {len(bh)} {bh}  clean hits {len(ch)} {ch}  dirty {len(dirty)}"
                      f"  ->  cons {cons:.4f}  ({time.monotonic()-b0:.0f}s)", flush=True)
                builds.append({"loop": loop, "band": got, "band_match": got == band,
                               "clean_n": len(clean), "clean_frac": len(clean) / len(window),
                               "clean_seeds": sorted(clean), "band_hits": bh, "clean_hits": ch,
                               "n_dirty": len(dirty), "cons": cons,
                               "cas9": meta.get("cas9_pool"), "pool": meta.get("pool")})

            out.append({"cell": cell, "task": tid[:8], "seeds": seeds, "classes": classes,
                        "real_classes": real, "padded_classes": padded, "arm": "narrow",
                        "window_size": len(window), "window": [min(window), max(window)],
                        "bank": prep["bank"], "found": {str(s): found.get(s) for s in seeds},
                        "two_seed_loop": two, "covered": len(covered),
                        "bands": bands, "builds": builds})
            json.dump(out, open(OUT, "w"), indent=1)

    print(f"\n=== summary ===", flush=True)
    for cell in CELLS:
        cr = [r for r in out if r["cell"] == cell]
        if not cr:
            continue
        bl = [b for r in cr for b in r["builds"]]
        twos = [r["two_seed_loop"] for r in cr if r["two_seed_loop"]]
        off = sum(3 - len(b["band_hits"]) for b in bl)
        ch = sum(len(b["clean_hits"]) for b in bl)
        print(f"  {cell:<12} 2-seed {len(twos)}/{len(cr)}  seeds "
              f"{sum(1 for r in cr for v in r['found'].values() if v)}/{3*len(cr)}  "
              f"| clean {np.mean([b['clean_n'] for b in bl]) if bl else 0:>5.1f}"
              f"/{cr[0]['window_size']} "
              f"({np.mean([b['clean_frac'] for b in bl]) if bl else 0:>5.1%})  "
              f"clean hits {ch}/{off}  cons {np.mean([b['cons'] for b in bl]) if bl else 0:.4f}",
              flush=True)
    print(f"\nwrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
