#!/usr/bin/env python3
"""conj_test.py — all-cut AND an HDR band on the SAME rows: does it keep all-cut's floor?

The target is the floor, not the spike. With a 0.10 floor a pure pin can only pay 0.40 / 0.70 /
1.00 on a 3-seed round, yet the leaders sit at 0.45-0.65 — which a floor near 0.25-0.35 generates
directly:

    all-HDR today   (1.00 + 0.13 + 0.13) / 3 = 0.42     measured off-band mean 0.1256
    all-cut today   (0.25 + 0.25 + 0.25) / 3 = 0.25     measured 0.2544 on 62% of seeds
    conjunction     (1.00 + 0.25 + 0.25) / 3 = 0.50     if both hold at once

HDR is a SUBSET of cut, so these are nested requirements on one row rather than competing ones:
a row that repairs by HDR at a seed necessarily cut at that seed. The construction is therefore a
sequential filter, not a second objective — which is what distinguishes it from the "combined
construction" already falsified in CLAUDE.md (that measured all-HDR's existing rows for
cut-cleanliness, and pitted two min-union objectives against each other).

The cost is pool shrinkage: requiring HDR on k band seeds keeps ~0.57**k of the all-cut bank, and a
smaller pool min-unions to a LARGER union of failed seeds, i.e. fewer cut-clean seeds and a lower
floor. That trade is the whole experiment, so k is swept with k=0 as the pure all-cut reference.

    python conj_test.py [task_id] [k1,k2,...]
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")

TASK = _ARGV[1] if len(_ARGV) > 1 else None
KS = [int(x) for x in (_ARGV[2].split(",") if len(_ARGV) > 2 else ["0", "1", "2", "4"])]
os.environ.setdefault("NIOME_INSTANCE", "conj")

import json  # noqa: E402
import logging  # noqa: E402
import random  # noqa: E402
import statistics as st  # noqa: E402
import time  # noqa: E402
from collections import Counter, defaultdict  # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np  # noqa: E402

import genExp as G  # noqa: E402
from niome_subnet.genomics import all_cut as AC  # noqa: E402
from niome_subnet.genomics import fastgreedy as FG  # noqa: E402
from niome_subnet.genomics import mt19937 as MT  # noqa: E402
from niome_subnet.genomics import seed_agnostic as SA  # noqa: E402
from niome_subnet.genomics.validation import stage3  # noqa: E402
from nocut100 import pick_task  # noqa: E402
from sd_task import score, task_content  # noqa: E402

OUT = os.getenv("CONJ_OUT", "conj_test.json")
# The cut window. Narrowing it makes the cut requirement cheaper per guide (fewer seeds to satisfy)
# and so leaves a bigger pool for the HDR filter — but it also shrinks the cut-clean SET, which is
# what actually pays, and everything outside the window falls back to the ~0.11 floor. cas12a_max_fail
# is scaled with the span: 100 is calibrated for 900 seeds, and left unscaled it would accept a guide
# failing every seed in a 100-wide window.
WINDOW = os.getenv("CONJ_WINDOW", "100-999")
# Group size and the range the band is chosen from, both overridable: the all-cut default of 42 is
# tuned for a whole-window build (83/17 cas mix), and a nested band needs its candidate seeds
# confined to the sub-window it is meant to live in rather than scattered across the cut window.
GROUP = int(os.getenv("CONJ_GROUP", "0")) or None
BAND_RANGE = os.getenv("CONJ_BAND_RANGE", "")
N_SAMPLE = int(os.getenv("CONJ_SAMPLE", "18"))     # off-band seeds scored per arm
BAND_POOL = int(os.getenv("CONJ_BAND_POOL", "40"))  # candidate seeds the band is chosen from
# Pin the band to named seeds instead of letting the greedy choose. This is ORACLE knowledge — the
# round's stamped seeds are not known at build time — so it prices what a correct seed-window
# prediction would be worth, and must never be read as an achievable score.
BAND_FORCE = [int(x) for x in os.getenv("CONJ_BAND_FORCE", "").replace(",", " ").split()]


def hdr_compliance(records, contract, cell_types, ctx, band_candidates):
    """For each bank guide, which of the candidate band seeds it repairs by HDR on.

    Grouped by (site, mutation) because the screen is per target: every guide of one target shares
    gc/energy/cut_p columns, which is what makes the batched kernel worth using.
    """
    acc = cell_types.get(contract.get("cell_type"), {}).get("accessibility", 1.0)
    regions = contract.get("mutation_regions") or {}
    seeds = np.asarray(band_candidates, dtype=np.int64)
    groups = defaultdict(list)
    for i, rec in enumerate(records):
        groups[(rec["start"], rec["strand"], rec["cas_system"], rec["length"],
                rec["mutation"])].append(i)
    ok = {}
    for (start, strand, cas, length, mutation), idxs in groups.items():
        site = type("S", (), {"start": start, "strand": strand, "cas": cas, "length": length})()
        distance = abs(start - ctx.mutation_map[mutation])
        offset = stage3.REGION_ENERGY_OFFSETS.get(regions.get(mutation), 0.0)
        params_of = AC._params_fn(site, distance, acc, offset)
        guides = [records[i]["guide"] for i in idxs]
        got = MT.screen_guides_rule_gpu(guides, seeds, mutation, cas, start, strand,
                                        params_of, "hdr", len(seeds))
        for i, guide in zip(idxs, guides):
            fails = got.get(guide)
            bad = set(int(x) for x in fails) if fails is not None else set(int(s) for s in seeds)
            ok[i] = set(int(s) for s in seeds) - bad
    return ok


def choose_band(ok, k, band_candidates, need):
    """Greedily pick k seeds keeping the largest set of guides HDR-compliant on ALL of them."""
    if k == 0:
        return [], set(ok)
    alive = set(ok)
    band = []
    for _ in range(k):
        best, best_keep = None, None
        for seed in band_candidates:
            if seed in band:
                continue
            keep = {i for i in alive if seed in ok[i]}
            if best_keep is None or len(keep) > len(best_keep):
                best, best_keep = seed, keep
        if best is None or len(best_keep) < need:
            break
        band.append(best)
        alive = best_keep
    return band, alive


def own_field(task_id):
    """The real scored finals of every OTHER miner on this exact task, best row per hotkey.

    Matched to the contract rather than pooled across tasks: the rank-10 cutoff swings 45.9-135.1
    across K562 rounds, so placing a contract-specific score against a pool of other contracts'
    fields measures field softness instead of the construction.
    """
    from sd_task import OURS
    try:
        rows = json.load(open("sd_task_scores.json"))
    except (OSError, json.JSONDecodeError):
        return []
    rows = rows if isinstance(rows, list) else (rows.get("data") or rows.get("items") or [])
    best = {}
    for r in rows:
        if r.get("task_id") == task_id and r.get("miner_hotkey") not in OURS:
            h = r["miner_hotkey"]
            best[h] = max(float(r.get("final_score") or 0), best.get(h, 0.0))
    return sorted(best.values(), reverse=True)


def round_score(rows, contract, reference, cell_types, field):
    """The ACTUAL round score on this contract's stamped seeds, and where it lands in the field.

    The validator averages every breakdown field over the round's seeds, so each seed's own
    weighted/consistency/fidelity are multiplied and only then averaged.
    """
    seeds = [int(x) for x in str(contract.get("seed", "")).split(",") if x.strip().isdigit()]
    if not seeds:
        return None
    per = [score(rows, contract, reference, cell_types, seed=s) for s in seeds]
    final = st.mean(p["final"] for p in per)
    rank = sum(1 for f in field if f > final) + 1 if field else None
    dist = [0.30, 0.20, 0.20, 0.15, 0.05, 0.03, 0.025, 0.02, 0.015, 0.01]
    return {"seeds": seeds, "per_seed_cons": [p["consistency"] for p in per],
            "consistency": st.mean(p["consistency"] for p in per),
            "final": final, "rank": rank, "field_n": len(field),
            "share": (dist[rank - 1] if rank and rank <= len(dist) else 0.0)}


def main():
    task, contract, reference = task_content(TASK) if TASK else pick_task()
    cell = contract.get("cell_type")
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    cfg = AC.config_for(cell) or AC.AllCutConfig()
    w_lo, w_hi = (int(x) for x in WINDOW.split("-"))
    span = w_hi - w_lo + 1
    import dataclasses as _dc
    # ``seeds`` is a property derived from start_seed/end_seed, not a field — setting it here
    # would raise, and setting the bounds is what actually moves the window.
    cfg = _dc.replace(cfg, start_seed=w_lo, end_seed=w_hi,
                      cas12a_max_fail=max(0, round(cfg.cas12a_max_fail * span / 900)),
                      **({"group_size": GROUP} if GROUP else {}))
    print(f"cut window {w_lo}-{w_hi} (span {span}), "
          f"cas12a_max_fail scaled to {cfg.cas12a_max_fail}")
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
    print(f"task {task['id'][:8]}  {cell}  group {cfg.group_size}  rows {n_rows}  "
          f"window {cfg.start_seed}-{cfg.end_seed}\n")

    path = os.path.join(AC.BANK_DIR, f"cas12a-{AC.bank_key(contract, cell_types, cfg)}.npz")
    if not os.path.exists(path):
        print("building the all-cut Cas12a bank (the slow half, ~120-300s)")
        bank = AC.build_bank(contract, reference, cell_types, ctx, sites, cfg)
        if not bank:
            raise SystemExit("bank scan produced nothing")
        AC.save_bank(path, bank)
    records = AC.load_bank(path)
    print(f"bank: {len(records)} Cas12a guides\n")

    rng = random.Random(11)
    field = own_field(task["id"])
    if field:
        cut10 = field[9] if len(field) >= 10 else 0.0
        print(f"own field: {len(field)} other miners, rank-1 {field[0]:.1f}, "
              f"rank-10 cutoff {cut10:.1f}")
    if BAND_RANGE:
        b_lo, b_hi = (int(x) for x in BAND_RANGE.split("-"))
    else:
        b_lo, b_hi = cfg.start_seed, cfg.end_seed
    band_candidates = sorted(rng.sample(range(b_lo, b_hi + 1),
                                        min(BAND_POOL, b_hi - b_lo + 1)))
    print(f"band candidates drawn from {b_lo}-{b_hi}")
    t0 = time.monotonic()
    ok = hdr_compliance(records, contract, cell_types, ctx, band_candidates)
    per_seed = [len([i for i in ok if s in ok[i]]) for s in band_candidates]
    print(f"HDR screen over {BAND_POOL} candidate band seeds in {time.monotonic()-t0:.0f}s")
    print(f"  bank guides HDR-compliant per single seed: median {int(st.median(per_seed))} "
          f"of {len(records)} ({st.median(per_seed)/len(records):.0%})\n")

    print(f"  {'k':>3} {'pool':>7} {'union':>6} {'clean':>6} {'cas9':>6} {'rows':>5} "
          f"{'cells':>5} {'weighted':>9} {'band':>7} {'clean cons':>10} {'dirty':>7} "
          f"{'offwin':>7} {'E[cons]':>8} {'E[final]':>9}")
    out = {}
    for k in (KS if not BAND_FORCE else [len(BAND_FORCE)]):
        if BAND_FORCE:
            band = list(BAND_FORCE)
            alive = {i for i in ok if all(s in ok[i] for s in band)}
            print(f"  band pinned to {band} (oracle: these are the round's own seeds)")
        else:
            band, alive = choose_band(ok, k, band_candidates, cfg.group_size)
        if k and len(band) < k:
            print(f"  {k:>3} {len(alive):>7}  band of {k} unreachable (pool fell below the group)")
            out[k] = {"reason": "band unreachable"}
            continue
        pool = [records[i] for i in sorted(alive)]
        if len(pool) < cfg.group_size:
            print(f"  {k:>3} {len(pool):>7}  pool below group size {cfg.group_size}")
            out[k] = {"reason": "pool below group size"}
            continue
        selector = FG.FastGreedy(pool, window_lo=cfg.start_seed, window_hi=cfg.end_seed)
        index, _union = selector.best(cfg.group_size, restarts=12)
        group = [pool[i] for i in index]
        bad = set()
        for rec in group:
            bad.update(int(x) for x in rec["fails"])
        clean = np.array(sorted(set(range(cfg.start_seed, cfg.end_seed + 1)) - bad),
                         dtype=np.int64)
        cas9 = AC.scan_cas9(clean, contract, cell_types, ctx, sites, cfg,
                            n_rows - cfg.group_size)
        # The Cas9 half must ALSO repair by HDR on the band, or the band seed is not a spike.
        if band and cas9:
            ok9 = hdr_compliance(cas9, contract, cell_types, ctx, band)
            cas9 = [cas9[i] for i in sorted(ok9) if len(ok9[i]) == len(band)]
        if len(cas9) < n_rows - cfg.group_size:
            print(f"  {k:>3} {len(pool):>7} {len(bad):>6} {len(clean):>6} {len(cas9):>6}"
                  f"  Cas9 half short of {n_rows - cfg.group_size}")
            out[k] = {"reason": f"cas9 pool {len(cas9)}", "pool": len(pool),
                      "clean": len(clean)}
            continue
        rows = AC.assemble(group, cas9, contract, ctx, cfg, n_rows)
        cells = len(Counter((r["mutation"], r["cas_system"], r["strand"]) for r in rows))
        cleanset = set(int(x) for x in clean)
        band_seed = band[0] if band else None
        s_band = score(rows, contract, reference, cell_types,
                       seed=band_seed) if band_seed else None
        pool_clean = [s for s in cleanset if s not in set(band)]
        pool_dirty = sorted(bad)
        cs = [score(rows, contract, reference, cell_types, seed=s)["consistency"]
              for s in rng.sample(pool_clean, min(N_SAMPLE, len(pool_clean)))]
        ds = [score(rows, contract, reference, cell_types, seed=s)["consistency"]
              for s in rng.sample(pool_dirty, min(N_SAMPLE, len(pool_dirty)))] or [0.10]
        # Off-window seeds are part of the round's draw and were never screened, so they must be
        # measured rather than assumed. Dividing by the window span instead of the full 100-999
        # space would silently credit a narrow window with the coverage of a wide one.
        offw = [s for s in range(100, 1000)
                if not (cfg.start_seed <= s <= cfg.end_seed)]
        os_ = [score(rows, contract, reference, cell_types, seed=s)["consistency"]
               for s in rng.sample(offw, min(N_SAMPLE, len(offw)))] if offw else []
        e_cons = ((len(band) * (s_band["consistency"] if s_band else 0.0)
                   + (len(cleanset) - len(band)) * st.mean(cs)
                   + len(bad) * st.mean(ds)
                   + len(offw) * (st.mean(os_) if os_ else 0.0)) / 900.0)
        base = s_band or score(rows, contract, reference, cell_types, seed=int(clean[0]))
        e_final = base["weighted"] * e_cons * base["fidelity"]
        rnd = round_score(rows, contract, reference, cell_types, field)
        out[k] = {"round": rnd,
                  "pool": len(pool), "union": len(bad), "clean": len(cleanset),
                  "cas9": len(cas9), "rows": len(rows), "cells": cells,
                  "weighted": base["weighted"], "fidelity": base["fidelity"],
                  "band_cons": s_band["consistency"] if s_band else None,
                  "clean_cons": st.mean(cs), "dirty_cons": st.mean(ds),
                  "offwindow_cons": (st.mean(os_) if os_ else None),
                  "offwindow_seeds": len(offw),
                  "e_cons": e_cons, "e_final": e_final, "band": band}
        print(f"  {k:>3} {len(pool):>7} {len(bad):>6} {len(cleanset):>6} {len(cas9):>6} "
              f"{len(rows):>5} {cells:>5} {base['weighted']:>9.1f} "
              f"{(s_band['consistency'] if s_band else float('nan')):>7.4f} "
              f"{st.mean(cs):>10.4f} {st.mean(ds):>7.4f} "
              f"{(st.mean(os_) if os_ else float('nan')):>7.4f} "
              f"{e_cons:>8.4f} {e_final:>9.2f}")
        if rnd:
            hits = [x for x in band if x in rnd["seeds"]]
            print(f"      round seeds {rnd['seeds']} -> cons "
                  f"{'/'.join(f'{c:.3f}' for c in rnd['per_seed_cons'])} = {rnd['consistency']:.4f}"
                  f"   final {rnd['final']:.1f}   rank {rnd['rank']} of {rnd['field_n'] + 1}"
                  f"   share {rnd['share']:.3f}"
                  + (f"   BAND HIT on {hits}" if hits else "   (band hit none)"))
        MT.free_gpu_memory()

    print(f"\n  E[cons] weights every regime by its share of the FULL 100-999 draw space (900"
          f" seeds), including the {900 - span} outside the cut window.")
    print(f"  all-HDR for comparison: band 12-16 at 1.000, off-band 0.1256 -> E[cons] ~0.15.")
    json.dump({"task": task["id"], "cell": cell, "group": cfg.group_size,
               "cut_window": [cfg.start_seed, cfg.end_seed],
               "band_range": [b_lo, b_hi], "bank": len(records), "arms": out},
              open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
