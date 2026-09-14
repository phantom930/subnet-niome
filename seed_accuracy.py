#!/usr/bin/env python3
"""seed_accuracy.py — do the REAL seeds land in our bands, given a correct window prediction?

Everything about the fleet's coverage so far is modelled: band sizes are measured, the union is
measured, and P(hit) is then computed as `1 - C(S-U,3)/C(S,3)` on the assumption that the three
seeds are uniform on the predicted space. This replays that against the seeds the backend actually
drew, on every three-seed task of a cell type.

**The prediction is GRANTED, not tested.** For each task the joined space is built from the classes
that really hold the three seeds, so the window model is correct by construction and what is left
is exactly the question asked: given the right three windows, how often is a drawn seed one of the
clean seeds our eleven hotkeys select out of them?

Two details that decide whether the number means anything:

* **Two-class rounds.** 9 of 31 K562 rounds put two seeds in one width-100 class, so only two
  classes hold seeds and the third predicted window is empty. Padding it with an arbitrary unused
  class spreads a third of the fleet's band onto seeds that cannot be drawn -- but it dilutes union
  and space by the same factor, so `union/space` and therefore P(hit) are unchanged in expectation.
  The padding is deterministic (lowest unused class) and is recorded per task.
* **Band, not submission.** A band seed is one where the Cas12a min-union group is clean; the Cas9
  fill lands on that band and cannot alter it (band_rebuild.py). Same definition as band_dupe.py.

    SA_CELL=K562 SA_LAYOUT=plan|tile27 python -u seed_accuracy.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "seedacc")

import dataclasses as _dc   # noqa: E402
import json                 # noqa: E402
import logging              # noqa: E402
import math                 # noqa: E402
import time                 # noqa: E402
from collections import Counter  # noqa: E402

logging.basicConfig(level=logging.ERROR)

import genExp as G          # noqa: E402
import joined_window as JW  # noqa: E402
from niome_subnet.genomics import all_hdr as AH     # noqa: E402
from niome_subnet.genomics import fastgreedy as FG  # noqa: E402
from niome_subnet.genomics import mt19937 as MT     # noqa: E402
from sd_task import task_content                    # noqa: E402

CELL = os.getenv("SA_CELL", "K562")
LAYOUT = os.getenv("SA_LAYOUT", "plan")
TILE_W = int(os.getenv("SA_TILE_W", "27"))
TILE_G = int(os.getenv("SA_TILE_G", "42"))
OUT = os.getenv("SA_OUT", f"seed_accuracy_{CELL.replace('+','plus').replace('-','_')}.json")
H, NS = 11, 3


def tasks():
    d = json.load(open("sd_task_listing.json"))
    items = d if isinstance(d, list) else (d.get("items") or [])
    got = []
    for t in items:
        c = (t.get("content") or {}).get("contract") or {}
        if c.get("cell_type") != CELL:
            continue
        sd = [int(x) for x in str(c.get("seed", "")).split(",") if x.strip().isdigit()]
        if len(sd) == 3:
            got.append((t.get("created_at", ""), t["id"], sd))
    return sorted(got)


def classes_of(seeds):
    """The three width-100 classes a correct prediction would name, padded if the seeds use fewer."""
    used = sorted({(s - 100) // 100 for s in seeds})
    pad = [c for c in range(9) if c not in used]
    return sorted(used + pad[:3 - len(used)]), len(used)


def windows_for(space):
    """{hotkey: seed list} under the layout being tested."""
    if LAYOUT == "tile27":
        return {f"h{i}": space[i * TILE_W:(i + 1) * TILE_W] for i in range(H)}
    width = JW.sub_width(CELL)
    got = {f"h{i+1}": JW.slice_at(space, (i * JW.STRIDE) % len(space), width) for i in range(H - 1)}
    got["h0"] = JW.half_stride(space, width)
    return {k: sorted({s for a, b in v for s in range(a, b + 1)}) for k, v in got.items()}


def band_for(seeds, contract, reference, cell_types, ctx, sites):
    base = AH.config_for(CELL)
    span = len(seeds)
    group = TILE_G if LAYOUT == "tile27" else base.group_size
    cfg = _dc.replace(base, hdr_range=(seeds[0], seeds[-1]), seed_list=tuple(seeds),
                      group_size=group,
                      main_max_fail=AH._scaled_max_fail(CELL, base.main_max_fail, span,
                                                        AH._native_span(CELL)))
    if span > 100:
        cfg = _dc.replace(cfg, variants=min(cfg.variants, AH.WIDE_WINDOW_VARIANTS))
    cfg = _dc.replace(cfg, light_cell_rows=AH.resolve_light(cfg.light_cell_rows, contract))
    p = os.path.join(AH.HDR_BANK_DIR, f"cas12a-{AH.bank_key(contract, cell_types, cfg)}.npz")
    if not os.path.exists(p):
        bk = AH.build_bank(contract, reference, cell_types, ctx, sites, cfg, None)
        if not bk:
            return None
        AH.save_bank(p, bk)
    recs = AH.load_bank(p)
    if len(recs) < cfg.group_size:
        return None
    sp = cfg.band_seeds()
    sel = FG.FastGreedy(recs, window_lo=cfg.start_seed, window_hi=cfg.end_seed,
                        per_cell_min=cfg.per_cell_min, caps=AH._group_caps(contract, ctx, cfg),
                        seeds=sp)
    idx, _u = sel.best(cfg.group_size, restarts=cfg.restarts)
    bad = set()
    for r in (recs[i] for i in idx):
        bad.update(int(x) for x in r["fails"])
    MT.free_gpu_memory()
    return sorted(set(int(x) for x in sp) - bad)


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    todo = tasks()
    print(f"{CELL}: {len(todo)} three-seed tasks, layout={LAYOUT}, "
          f"{H} hotkeys — window prediction GRANTED\n")
    out = {"cell": CELL, "layout": LAYOUT, "rounds": []}
    for i, (created, tid, seeds) in enumerate(todo, 1):
        cls, n_used = classes_of(seeds)
        space = sorted(s for c in cls for s in range(c * 100 + 100, c * 100 + 200))
        wins = windows_for(space)
        _t, contract, reference = task_content(tid)
        ctx = G.build_context(contract, reference, cell_types)
        sites = G.enumerate_sites(ctx, 3000, (20, 23))
        t0 = time.monotonic()
        bands, per_hk = {}, {}
        for hk, sd in wins.items():
            b = band_for(sd, contract, reference, cell_types, ctx, sites)
            if b is None:
                continue
            bands[hk] = b
            per_hk[hk] = sum(1 for s in seeds if s in set(b))
        if not bands:
            print(f"  {i:>2}/{len(todo)} {tid[:8]} all hotkeys declined")
            continue
        union = {s for b in bands.values() for s in b}
        covered = [s for s in seeds if s in union]
        kmax = max(per_hk.values())
        rec = {"task": tid, "created": created, "seeds": seeds, "classes": cls,
               "distinct_classes": n_used, "union": len(union),
               "band_sizes": {k: len(v) for k, v in bands.items()},
               "seeds_covered": covered, "n_covered": len(covered),
               "k_max": kmax, "per_hotkey": per_hk,
               "elapsed_s": round(time.monotonic() - t0, 1)}
        out["rounds"].append(rec)
        print(f"  {i:>2}/{len(todo)} {tid[:8]} seeds {seeds} cls {cls}{'*' if n_used < 3 else ''}"
              f"  union {len(union):>3}  covered {len(covered)}/3  kmax {kmax}"
              f"  {rec['elapsed_s']:.0f}s", flush=True)
        json.dump(out, open(OUT, "w"), indent=1)

    r = out["rounds"]
    if r:
        n = len(r)
        seeds_tot = 3 * n
        hit_seeds = sum(x["n_covered"] for x in r)
        p_round = sum(1 for x in r if x["n_covered"] >= 1) / n
        p_dbl = sum(1 for x in r if x["k_max"] >= 2) / n
        mean_u = sum(x["union"] for x in r) / n
        c = math.comb
        model1 = 1 - c(300 - round(mean_u), NS) / c(300, NS)
        print(f"\n{CELL}  {n} rounds, layout {LAYOUT}")
        print(f"  mean union                 {mean_u:.1f} of 300")
        print(f"  per-SEED accuracy          {hit_seeds}/{seeds_tot} = {hit_seeds/seeds_tot:.1%}"
              f"   (model {mean_u/300:.1%})")
        print(f"  round hit  (>=1 seed)      {p_round:.1%}   (model {model1:.1%})")
        print(f"  double hit (one hk >=2)    {p_dbl:.1%}")
        print(f"  k_max distribution         {dict(sorted(Counter(x['k_max'] for x in r).items()))}")
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
