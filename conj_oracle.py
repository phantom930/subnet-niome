#!/usr/bin/env python3
"""conj_oracle.py — the conjunction under ORACLE window knowledge, across a fleet of 10 hotkeys.

Assumes what no live miner knows: the three 100-seed windows the round's seeds fall in. The joined
cut space is exactly those windows (200 seeds when two seeds share one, 300 otherwise), every
hotkey min-unions on cut over the WHOLE of it, and the hotkeys are decorrelated only by where their
HDR band sub-window sits — rotated at `STRIDE` so 10 slices of `HDR_WIDTH` tile the space.

    h0: indices   0..149      h1: indices  30..179      h2: indices  60..209   ...

What this measures is coverage, not payout: per task it unions the 10 hotkeys' bands and the 10
hotkeys' clean sets, then asks how many of the round's three real seeds each union contains. A band
seed scores `consistency_factor` 1.0 and a clean seed ~0.21-0.24, so band>=2 is the round that wins
outright and clean>=1 is the round that is merely elevated.

**`light_cell_rows` is deliberately NOT swept.** It governs only how `all_cut.assemble` apportions
rows across (mutation, Cas9, strand) cells; the band comes from the greedy HDR prefix and the clean
set from the Cas12a min-union, both settled before `assemble` is called. It cannot move either
number, so sweeping it would have quadrupled the run for four identical answers.

**The Cas9 half is not built here either**, so these are coverage upper bounds: a (group, k) whose
Cas9 conditional fill cannot reach `n_rows - group` would decline in the miner and ship nothing.
Feasibility is a separate, much more expensive question -- see `--feasible`.

    python conj_oracle.py --tasks 1              # pilot, times each phase
    CO_GROUPS=80,100,125 CO_KMAX=9 python conj_oracle.py --tasks 40
"""
import argparse
import dataclasses as _dc
import json
import logging
import os
import sys
import time

sys.argv_full = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "conjoracle")
logging.basicConfig(level=logging.ERROR)

import numpy as np                                          # noqa: E402

import genExp as G                                          # noqa: E402
from niome_subnet.genomics import all_cut as AC             # noqa: E402
from niome_subnet.genomics import conjunction as CJ         # noqa: E402
from niome_subnet.genomics import fastgreedy as FG          # noqa: E402
from niome_subnet.genomics import mt19937 as MT             # noqa: E402
from joined300 import fetch_tasks                           # noqa: E402

CELL = os.getenv("CO_CELL", "HEK293")
GROUPS = [int(x) for x in os.getenv("CO_GROUPS", "80,100,125").split(",")]
KMAX = int(os.getenv("CO_KMAX", "9"))
KMIN = int(os.getenv("CO_KMIN", "6"))
HDR_WIDTH = int(os.getenv("CO_HDR_WIDTH", "150"))
STRIDE = int(os.getenv("CO_STRIDE", "30"))
NHK = int(os.getenv("CO_NHK", "10"))
MF_SPAN900 = int(os.getenv("CO_MF", "100"))     # scaled to the joined span, as build_for_cell does
RESTARTS = int(os.getenv("CO_RESTARTS", "12"))
OUT = os.getenv("CO_OUT", "conj_oracle_{cell}.json")


def slice_at(seeds, offset, width):
    """The circular width-`width` slice of `seeds` starting at index `offset`."""
    n = len(seeds)
    return [seeds[(offset + j) % n] for j in range(min(width, n))]


def greedy_prefix(ok, kmax, candidates, need):
    """Greedy HDR band to depth kmax, returning (band, alive) at EVERY prefix length.

    Run once per hotkey at the SMALLEST group in the sweep: the picks a greedy makes do not depend
    on `need`, only where it stops, so the prefixes of a need=80 run contain the need=100 and
    need=125 runs exactly. The caller truncates by `len(alive) >= group` instead of re-running.
    """
    cand = sorted(set(int(x) for x in candidates))
    col = {s: j for j, s in enumerate(cand)}
    idx = sorted(ok)
    M = np.zeros((len(idx), len(cand)), dtype=bool)
    for r, i in enumerate(idx):
        for s in ok[i]:
            if s in col:
                M[r, col[s]] = True
    alive = np.ones(len(idx), dtype=bool)
    band, taken = [], set()
    steps = [([], [idx[r] for r in np.flatnonzero(alive)])]
    for _ in range(kmax):
        counts = M[alive].sum(axis=0)
        for j in taken:
            counts[j] = -1
        j = int(np.argmax(counts))
        if counts[j] < need:
            break
        taken.add(j)
        band.append(cand[j])
        alive &= M[:, j]
        steps.append((list(band), [idx[r] for r in np.flatnonzero(alive)]))
    return steps


def windows_of(seeds):
    """The distinct 100-seed classes the round's seeds fall in -> the sorted joined seed space."""
    cls = sorted({(int(s) - 100) // 100 for s in seeds})
    return sorted(s for c in cls for s in range(c * 100 + 100, c * 100 + 200)), len(cls)


def run_task(task, cell_types, timing=None):
    """-> {(group, k): {'band': set, 'clean': set, 'per_hk': [...]}} plus meta, for one task."""
    content = task["content"]
    contract, reference = content["contract"], content["hbb_reference"]
    seeds = sorted(int(x) for x in str(contract.get("seed", "")).split(",") if x.strip().isdigit())
    joined, nwin = windows_of(seeds)
    t0 = time.monotonic()
    cfg = CJ.config_for(cell) if (cell := contract.get("cell_type")) else None
    mf = max(1, round(MF_SPAN900 * len(joined) / 900))
    cfg = _dc.replace(cfg, seed_list=tuple(joined), start_seed=joined[0], end_seed=joined[-1],
                      cas12a_max_fail=mf)
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    os.makedirs(AC.BANK_DIR, exist_ok=True)
    path = os.path.join(AC.BANK_DIR, f"cas12a-{AC.bank_key(contract, cell_types, cfg)}.npz")
    if not os.path.exists(path):
        bank = AC.build_bank(contract, reference, cell_types, ctx, sites, cfg)
        MT.free_gpu_memory()
        if not bank:
            return None, {"reason": f"cut bank empty at mf {mf}"}
        AC.save_bank(path, bank)
    records = AC.load_bank(path, limit=cfg.bank_keep)
    t_bank = time.monotonic() - t0

    t0 = time.monotonic()
    ok = CJ.hdr_compliance(records, contract, cell_types, ctx, joined)
    MT.free_gpu_memory()
    t_screen = time.monotonic() - t0

    t0 = time.monotonic()
    out = {}
    per_hk = {}
    jset = set(joined)
    for h in range(NHK):
        cand = set(slice_at(joined, (h * STRIDE) % len(joined), HDR_WIDTH))
        okw = {i: (v & cand) for i, v in ok.items()}
        steps = greedy_prefix(okw, KMAX, sorted(cand), min(GROUPS))
        for g in GROUPS:
            for k in range(KMIN, KMAX + 1):
                key = (g, k)
                rec = out.setdefault(key, {"band": set(), "clean": set(), "hk": 0})
                if k >= len(steps):
                    continue
                band, alive = steps[k]
                if len(alive) < g:
                    continue
                pool = [records[i] for i in alive]
                sel = FG.FastGreedy(pool, window_lo=joined[0], window_hi=joined[-1],
                                    seeds=np.asarray(joined, dtype=np.int64))
                idx, _u = sel.best(g, restarts=RESTARTS)
                bad = set()
                for r in (pool[i] for i in idx):
                    bad.update(int(x) for x in r["fails"])
                rec["band"].update(band)
                rec["clean"].update(jset - bad)
                rec["hk"] += 1
                per_hk.setdefault(key, []).append({"hk": h, "band": len(band),
                                                   "clean": len(jset - bad)})
    t_build = time.monotonic() - t0
    if timing is not None:
        timing.update(bank=t_bank, screen=t_screen, build=t_build, bank_n=len(records))
    meta = {"task": (task.get("task_id") or task["id"])[:8], "at": task.get("created_at", "")[:16],
            "seeds": seeds, "n_windows": nwin, "span": len(joined), "mf": mf,
            "bank": len(records)}
    return (out, per_hk), meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tasks", type=int, default=40)
    ap.add_argument("--out", default=None)
    args = ap.parse_args(sys.argv_full[1:])

    items = fetch_tasks(limit=500)
    cands = []
    for t in sorted(items, key=lambda x: x.get("created_at", ""), reverse=True):
        c = (t.get("content") or {}).get("contract") or {}
        if c.get("cell_type") != CELL:
            continue
        sd = [x for x in str(c.get("seed", "") or "").split(",") if x.strip().isdigit()]
        if len(sd) != 3:
            continue
        cands.append(t)
        if len(cands) >= args.tasks:
            break
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    print(f"{CELL}: {len(cands)} tasks | {NHK} hotkeys, HDR width {HDR_WIDTH} stride {STRIDE} "
          f"| groups {GROUPS} k {KMIN}-{KMAX} | mf {MF_SPAN900}@900\n")
    rows = []
    for n, t in enumerate(cands, 1):
        timing = {}
        got, meta = run_task(t, cell_types, timing)
        if got is None:
            print(f"  {meta.get('reason')}")
            continue
        out, per_hk = got
        rec = {**meta, "timing": {k: round(v, 1) for k, v in timing.items() if k != "bank_n"},
               "arms": []}
        sset = set(meta["seeds"])
        for (g, k), v in sorted(out.items()):
            b, c = v["band"], v["clean"]
            rec["arms"].append({
                "group": g, "k": k, "hotkeys_built": v["hk"],
                "band_total": len(b), "clean_total": len(c),
                "band_hits": len(b & sset), "clean_hits": len(c & sset),
                "per_hk_band": [x["band"] for x in per_hk.get((g, k), [])],
                "per_hk_clean": [x["clean"] for x in per_hk.get((g, k), [])],
            })
        rows.append(rec)
        best = max(rec["arms"], key=lambda a: (a["band_hits"], a["clean_hits"]))
        print(f"  [{n}/{len(cands)}] {meta['task']} {meta['at']} seeds {meta['seeds']} "
              f"win {meta['n_windows']} span {meta['span']} | bank {meta['bank']} "
              f"| bank {timing['bank']:.0f}s screen {timing['screen']:.0f}s "
              f"build {timing['build']:.0f}s | best band hits {best['band_hits']}/3 "
              f"(g{best['group']} k{best['k']}) clean {best['clean_hits']}/3")
        json.dump(rows, open(args.out or OUT.format(cell=CELL.replace("+", "")), "w"), indent=1)
    print(f"\nwrote {args.out or OUT.format(cell=CELL.replace('+', ''))}")


if __name__ == "__main__":
    main()
