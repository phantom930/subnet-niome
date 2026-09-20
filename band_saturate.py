#!/usr/bin/env python3
"""band_saturate.py — how many times can h0-h5 re-band a contract before the window is covered?

The procedure, as specified: take the union of h0-h5's bands; re-run the band step with those seeds
EXCLUDED from every hotkey's candidate window; union the new bands with the old; repeat. Two things
are read off it — the loop at which the round's REAL in-window seeds first get covered, and how far
the cumulative union ever gets before band formation collapses.

**What the loop count actually means, because it is not a per-round strategy.** A hotkey submits ONE
band per round, so N loops is not something one fleet does over N rounds -- the seeds are redrawn and
the bank is rebuilt every round, so nothing accumulates across them. N loops is what **6N hotkeys**
would do on ONE contract, each handed a candidate window with its predecessors' seeds removed. So
this measures the marginal value of fleet headcount on the predicted window, which CLAUDE.md calls
the one lever with headroom left, and it prices it exactly: the loop that first covers a drawn seed
is the headcount at which that seed would have been caught.

**Why it is cheap.** `choose_band` and `cas9_cell_probe` both restrict to their `candidates`
argument by column lookup (`if seed in col`), so compliance computed ONCE over the whole 300-seed
window is valid for every subset of it. One bank scan and one `hdr_compliance` pass per task, then
each loop is just the greedy -- which is why this can sweep several tasks instead of one.

**The exclusion makes the wall bite, and that is the point rather than a flaw.** Each loop hands the
greedy a thinner candidate list (150 minus whatever is already covered), and `choose_band` stops
when the survivors fall under `group_size`. A hotkey that cannot reach `band_k` DECLINES and
contributes nothing, exactly as `build_submission` would -- so the union saturates rather than
growing forever, and where it saturates is the answer.

    BS_N=4 python band_saturate.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "saturate")

import dataclasses                                       # noqa: E402
import json                                              # noqa: E402
import logging                                           # noqa: E402
import time                                              # noqa: E402

logging.basicConfig(level=logging.ERROR)

import genExp as G                                       # noqa: E402
import joined_window as JW                               # noqa: E402
from niome_subnet.genomics import conjunction as CJ      # noqa: E402
from niome_subnet.genomics import mt19937 as MT          # noqa: E402
from niome_subnet.genomics.validation import _parse_seeds  # noqa: E402
from conj_stageb import API                              # noqa: E402
from sd_task import fetch                                # noqa: E402

N_TASKS = int(os.getenv("BS_N", "4"))
MAX_LOOPS = int(os.getenv("BS_MAX_LOOPS", "40"))
HKS = JW.BAND_HK                                          # h0-h5, the predicted-window group
OUT = os.getenv("BS_JSON", "band_saturate.json")


def prepare(contract, reference, cell_types, cell, predicted):
    """Bank + per-guide HDR compliance + the cell probe, once, over the WHOLE 300-seed window."""
    base = CJ.config_for(cell)
    # h0-h5 all share one cut space within a cell, so one bank covers every hotkey and every loop.
    cut = JW.conjunction_cut_seeds(HKS[0], cell, predicted=predicted,
                                   band_candidates=predicted) or predicted
    cfg = dataclasses.replace(base, seed_list=tuple(sorted(cut)),
                              band_candidates=tuple(predicted))
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    path = os.path.join(CJ.BANK_DIR,
                        f"cas12a-{CJ.bank_key(contract, cell_types, cfg)}.npz")
    if not os.path.exists(path):
        os.makedirs(CJ.BANK_DIR, exist_ok=True)
        bank = CJ.build_bank(contract, reference, cell_types, ctx, sites, cfg, None)
        MT.free_gpu_memory()
        if not bank:
            return None
        CJ.save_bank(path, bank)
    records = CJ.load_bank(path, limit=cfg.bank_keep)
    ok = CJ.hdr_compliance(records, contract, cell_types, ctx, predicted, cfg.band_rule)
    MT.free_gpu_memory()
    cell_ok = None
    if cfg.band_cell_aware:
        cell_ok = CJ.cas9_cell_probe(contract, cell_types, ctx, sites, cfg, predicted,
                                     cfg.band_rule)
        MT.free_gpu_memory()
    return {"cfg": cfg, "ok": ok, "cell_ok": cell_ok, "bank": len(records), "cut": len(cut)}


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    items = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])
    plan = json.load(open("data/window_plan.json"))

    tasks = []
    for t in sorted(items, key=lambda x: x.get("created_at", ""), reverse=True):
        c = (t.get("content") or {}).get("contract") or {}
        if len(_parse_seeds(c.get("seed"))) != 3:
            continue
        tasks.append(t)
        if len(tasks) >= N_TASKS:
            break

    out = []
    for t in tasks:
        tid = (t.get("task_id") or t["id"])
        contract = dict(t["content"]["contract"])
        reference = t["content"]["hbb_reference"]
        cell = contract["cell_type"]
        seeds = _parse_seeds(contract["seed"])
        raw = plan["assignments"][cell][HKS[0]]
        pairs = raw if isinstance(raw[0], (list, tuple)) else [raw]
        predicted = sorted({s for a, b in pairs for s in range(a, b + 1)})
        target = [s for s in seeds if s in set(predicted)]

        t0 = time.monotonic()
        # The build must not see the stamped seeds: the miner prefetches before they exist.
        prep = prepare(dict(contract, seed=0), reference, cell_types, cell, predicted)
        if prep is None:
            print(f"\n=== {tid[:8]} {cell}: bank empty, skipped ===", flush=True)
            continue
        cfg, ok, cell_ok = prep["cfg"], prep["ok"], prep["cell_ok"]
        print(f"\n=== {tid[:8]}  {cell}  seeds {seeds}  k={cfg.band_k} group={cfg.group_size} "
              f"| window {pairs} ({len(predicted)}) | bank {prep['bank']} cut {prep['cut']} "
              f"| prep {time.monotonic()-t0:.0f}s ===", flush=True)
        print(f"    real seeds inside the window: {target or 'NONE'}", flush=True)

        covered, loops, found = set(), [], {}
        for loop in range(1, MAX_LOOPS + 1):
            new, built, declined = set(), 0, 0
            for hk in HKS:
                win = CJ.sub_window(predicted, cfg.band_width,
                                    JW.band_offset_frac(hk) or 0.0)
                cands = [s for s in win if s not in covered]
                if len(cands) < cfg.band_k:
                    declined += 1
                    continue
                band, _alive = CJ.choose_band(ok, cfg.band_k, cands, cfg.group_size, cell_ok,
                                              {}, cfg.band_quota)
                if len(band) < cfg.band_k:
                    declined += 1          # build_submission would decline here
                    continue
                built += 1
                new |= set(int(x) for x in band)
            covered |= new
            for s in target:
                if s in covered and s not in found:
                    found[s] = loop
            loops.append({"loop": loop, "built": built, "declined": declined,
                          "new": len(new), "covered": len(covered)})
            print(f"    loop {loop:>2}  built {built}/6  new {len(new):>3}  "
                  f"covered {len(covered):>3}/{len(predicted)} "
                  f"({len(covered)/len(predicted):5.1%})"
                  + (f"   <- caught {[s for s in target if found.get(s) == loop]}"
                     if any(found.get(s) == loop for s in target) else ""), flush=True)
            if built == 0:
                print(f"    SATURATED: no hotkey can form a {cfg.band_k}-seed band from what is "
                      f"left ({len(predicted)-len(covered)} seeds uncovered)", flush=True)
                break
            if len(covered) >= len(predicted):
                print(f"    FULL COVERAGE at loop {loop}", flush=True)
                break
        miss = [s for s in target if s not in found]
        print(f"    -> caught {found or '{}'}"
              + (f", NEVER caught {miss}" if miss else "")
              + f" | saturation {len(covered)}/{len(predicted)} after {len(loops)} loops",
              flush=True)
        out.append({"task": tid[:8], "cell": cell, "seeds": seeds, "window": pairs,
                    "k": cfg.band_k, "target": target, "found": found, "missed": miss,
                    "loops": loops, "covered": len(covered), "window_size": len(predicted)})

    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\n=== summary ===", flush=True)
    for r in out:
        cov = r["covered"] / r["window_size"]
        fnd = "/".join(f"{s}@L{r['found'][str(s)] if str(s) in r['found'] else r['found'].get(s)}"
                       for s in r["target"]) or "no in-window seed"
        print(f"  {r['task']} {r['cell']:12} {len(r['loops']):>3} loops  "
              f"saturate {cov:5.1%}  seeds {fnd}", flush=True)
    print(f"\nwrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
