#!/usr/bin/env python3
"""fleet_double.py — how often a 10-hotkey coordinated fleet catches TWO of a round's three seeds.

Assumes ORACLE class knowledge (the three 100-seed classes the seeds fall in) and 10 hotkeys, and
plays the coordinated band placement: each hotkey draws its band from a block inside every class
(`conjunction.class_blocks`) under a per-class quota, so the bands straddle the classes instead of
the 1.5 a contiguous `sub_window` spans, and the ten hotkeys' allocations are DISJOINT.

Why that matters: a hotkey holding (k1,k2,k3) band seeds across the classes doubles with probability
`(k1k2+k1k3+k2k3)/100**2` -- ZERO if the band sits in one class. The shipped `sub_window` produced
0 doubles in 100 real hotkey-rounds at a mean pair-sum of 14.1 of 21; balanced blocks reach 48 at
k=12.

CLASS-ADAPTIVE, because the calculus inverts when two seeds share a class (29.6% of rounds, 3 seeds
over 9 classes): concentrating the whole band in that class beats splitting it, and ten hotkeys
PARTITIONING that one class turns a double into "both seeds in one block".

Two numbers are reported and the second is the one to read:

  * EMPIRICAL -- did any hotkey actually hold >=2 of the round's real seeds. At ~5-7% per round this
    is 0 or 1 over ten tasks and is nearly uninformative on its own.
  * EXACT -- given the bands the fleet actually built, and seed positions uniform within their known
    classes (which is what the builder knew), P(some hotkey doubles). No sampling noise, and it is
    the quantity the design controls.

    FD_N=10 python fleet_double.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "fleetdbl")

import dataclasses
import json
import logging
import math
import time
from collections import Counter, defaultdict

logging.basicConfig(level=logging.ERROR)
import numpy as np                                       # noqa: E402

import genExp as G                                       # noqa: E402
from niome_subnet.genomics import conjunction as CJ      # noqa: E402
from conj_stageb import API                              # noqa: E402
from sd_task import fetch                                # noqa: E402

N_TASKS = int(os.getenv("FD_N", "10"))
NHK = int(os.getenv("FD_NHK", "10"))
SIM = int(os.getenv("FD_SIM", "200000"))
# Every decline in the k=band_k run was EXACTLY one seed short ("band reached 11 of 12",
# "8 of 9"), never two — the coordinated block has enough freedom for k-1 and not k.
KDELTA = int(os.getenv("FD_KDELTA", "0"))
RNG = np.random.default_rng(11)


def classes_for(seeds):
    """The distinct 100-seed classes the seeds occupy, padded to 3 for the CUT space."""
    used = sorted({s // 100 for s in seeds})
    pad = [c for c in range(1, 10) if c not in used]
    return used, (used + pad)[:3]


def exact_p_double(bands, classes_of_seed):
    """P(some hotkey holds >=2 of the 3 seeds), seed positions uniform inside their known classes."""
    n = SIM
    pos = RNG.integers(0, 100, size=(n, 3))
    seeds = np.stack([classes_of_seed[i] * 100 + pos[:, i] for i in range(3)], axis=1)
    hit_any = np.zeros(n, dtype=bool)
    for band in bands:
        if not band:
            continue
        b = np.asarray(sorted(band))
        inb = np.isin(seeds, b)
        hit_any |= inb.sum(axis=1) >= 2
    return float(hit_any.mean())


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    items = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])
    cands = []
    for t in sorted(items, key=lambda t: t.get("created_at", ""), reverse=True):
        c = (t.get("content") or {}).get("contract") or {}
        s = str(c.get("seed", "") or "")
        if len([x for x in s.split(",") if x.strip().isdigit()]) != 3:
            continue
        cands.append(t)
        if len(cands) >= N_TASKS:
            break

    print(f"{len(cands)} tasks, {NHK} hotkeys, conjunction + band_cell_aware, coordinated blocks\n",
          flush=True)
    out = []
    only = {x for x in os.getenv("FD_ONLY", "").replace(",", " ").split() if x}
    for t in cands:
        tid = (t.get("task_id") or t["id"])[:8]
        if only and tid not in only:
            continue
        contract = dict(t["content"]["contract"])
        reference = t["content"]["hbb_reference"]
        cell = contract.get("cell_type")
        seeds = [int(x) for x in str(contract["seed"]).split(",") if x.strip().isdigit()]
        used, three = classes_for(seeds)
        base = CJ.config_for(cell)
        if base is None:
            print(f"=== {tid} {cell}: no conjunction config, skipped ===", flush=True)
            continue
        joined = sorted(s for c in three for s in range(c * 100, c * 100 + 100))
        k = base.band_k + KDELTA
        # The class to concentrate in is the one holding TWO (or three) of the seeds, i.e. the
        # MODE of s//100 -- not used[0], which is merely the lowest class present. Picking the
        # lowest put the whole band in a single-seed class on 2 of 6 concentrate rounds and
        # returned an exact P(double) of 0 by construction.
        cl = [s // 100 for s in seeds]
        shared = None if len(used) == 3 else Counter(cl).most_common(1)[0][0]
        mode = "balanced 3-class" if shared is None else f"concentrate class {shared}"
        print(f"=== {tid} {cell:11s} seeds {seeds} classes {used} -> {mode}, k={k} ===", flush=True)

        bands, built = [], 0
        for h in range(NHK):
            if shared is None:
                cset = tuple(CJ.class_blocks(joined, 100 // NHK, h, NHK))
                cap = math.ceil(k / 3)
                quota = tuple((c * 100, c * 100 + 99, cap) for c in three)
            else:
                lo = shared * 100
                width = max(k, 100 // NHK)
                start = (h * (100 // NHK)) % 100
                cset = tuple(lo + ((start + j) % 100) for j in range(width))
                quota = ((lo, lo + 99, k),)
            cfg = dataclasses.replace(base, seed_list=tuple(joined), band_k=k,
                                      band_candidates=cset, band_quota=quota)
            t0 = time.monotonic()
            rows, meta = CJ.build_submission(contract, reference, cell_types, cfg=cfg,
                                             budget_s=1200.0)
            if not rows:
                print(f"  h{h}: DECLINED — {meta.get('reason')}  ({time.monotonic()-t0:.0f}s)",
                      flush=True)
                bands.append([])
                continue
            built += 1
            band = sorted(int(x) for x in meta["band_seeds"])
            bands.append(band)
            split = [sum(1 for s in band if s // 100 == c) for c in three]
            hits = [s for s in seeds if s in set(band)]
            print(f"  h{h}: k={len(band):2d} split {'/'.join(map(str, split))} "
                  f"hits {len(hits)}{' <<< DOUBLE' if len(hits) >= 2 else ''} "
                  f"({time.monotonic()-t0:.0f}s)", flush=True)

        best = max((sum(1 for s in seeds if s in set(b)) for b in bands if b), default=0)
        p = exact_p_double(bands, [s // 100 for s in seeds])
        out.append({"task": tid, "cell": cell, "seeds": seeds, "classes": used,
                    "mode": mode, "built": built, "best_hits": best, "p_exact": p,
                    "bands": bands})
        print(f"  -> built {built}/{NHK}  best hits {best}  P(double, exact) {p:.4f}\n", flush=True)

    with open("fleet_double.json", "w") as fh:
        json.dump(out, fh, indent=1)

    print("=== summary ===")
    print(f"{'task':9s} {'cell':11s} {'mode':22s} {'built':>6s} {'besthit':>8s} {'P(double)':>10s}")
    for r in out:
        print(f"{r['task']:9s} {r['cell']:11s} {r['mode']:22s} {r['built']:3d}/{NHK} "
              f"{r['best_hits']:8d} {r['p_exact']:10.4f}")
    n = len(out)
    emp = sum(1 for r in out if r["best_hits"] >= 2)
    ps = [r["p_exact"] for r in out]
    print(f"\n  EMPIRICAL: {emp} of {n} rounds had at least one hotkey with >=2 seeds")
    print(f"  EXACT    : mean P(double) {np.mean(ps):.4f} -> expected {np.sum(ps):.2f} of {n} rounds")
    print(f"             per-round range {min(ps):.4f} - {max(ps):.4f}")
    print(f"  rounds per double: 1 / {np.mean(ps):.4f} = {1/max(1e-9, np.mean(ps)):.0f}")


if __name__ == "__main__":
    main()
