#!/usr/bin/env python3
"""both_seeds.py — can ONE hotkey band two specific seeds, 214 and 249, on f18ee409?

Seven windows across M11-M13 and M16 have contained two drawn seeds and none banded both. This
asks whether that is bad luck or structure, on one concrete HEK293 round (seeds 214, 249, 784).

The obstacle is arithmetic. 214 and 249 are 35 apart, so the narrowest window holding both is 36
wide, and HEK293's `band_k` is 8 against a wall of 8-9. A band of 8 in a window of 36 is 22%
density, so a blind greedy picks both with probability ~4-6%. Reaching 50% would need k >= 26 at
width 36 -- three times past the wall.

Arms, all on the live cut (200-299 + 400-499), everything else HEK293's live CELL_CONFIG:

    oracle8   band_candidates = {214, 249} + 6 fillers   k=8   -- the greedy has no choice
    oracle2   band_candidates = {214, 249}               k=2   -- the pure two-seed band
    w36/k8    band_candidates = 214..249                 k=8   -- narrowest blind window
    w36/k9    band_candidates = 214..249                 k=9   -- and at HEK293's other depth
    w40/k8    band_candidates = 210..249                 k=8
    w100/k8   band_candidates = 200..299                 k=8   -- the whole class

The oracle arms are NOT buildable blind: they need the round's seeds, which are stamped after
broadcast. They are here to price the ceiling -- what a k=2 round is worth on this contract -- so
the blind arms can be judged against something real rather than against each other.
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "bothseeds")

import dataclasses, gc, json, logging, random, time                   # noqa: E402

logging.basicConfig(level=logging.ERROR)
import numpy as np                                                    # noqa: E402

import genExp as G                                                    # noqa: E402
from niome_subnet.genomics import conjunction as CJ                   # noqa: E402
from niome_subnet.genomics import mt19937 as MT                       # noqa: E402
from niome_subnet.genomics.validation import _parse_seeds             # noqa: E402
from conj_stageb import API                                           # noqa: E402
from sd_task import fetch, score                                      # noqa: E402
from widecut_price import records_of                                  # noqa: E402

TASK = os.getenv("BS_TASK", "f18ee409")
OUT = os.getenv("BS_JSON", "both_seeds.json")
SPACE = sorted(list(range(200, 300)) + list(range(400, 500)))   # cut + clean, the live span
A, B = 214, 249


def compliant(rows, contract, cell_types, ctx, rule):
    ok = CJ.hdr_compliance(records_of(rows), contract, cell_types, ctx, SPACE, rule)
    MT.free_gpu_memory()
    if not ok:
        return []
    keep = set(SPACE)
    for i in sorted(ok):
        keep &= set(int(x) for x in ok[i])
        if not keep:
            break
    return sorted(keep)


def price(rows, contract, reference, cell_types, ctx, seeds):
    band = compliant(rows, contract, cell_types, ctx, "hdr")
    cut = compliant(rows, contract, cell_types, ctx, "cut")
    bs, cs = set(band), set(cut) - set(band)
    per = [dict(seed=s, **{k: v for k, v in
                           score(rows, contract, reference, cell_types, seed=s).items()})
           for s in seeds]
    cons = float(np.mean([p["consistency"] for p in per]))
    wxf = per[0]["weighted"] * per[0]["fidelity"]
    return {"band_n": len(band), "band": band, "clean_n": len(cs),
            "band_hits": [s for s in seeds if s in bs],
            "clean_hits": [s for s in seeds if s in cs],
            "has_A": A in bs, "has_B": B in bs, "both": (A in bs and B in bs),
            "round_cons": cons, "weighted": per[0]["weighted"],
            "fidelity": per[0]["fidelity"], "round_final": wxf * cons,
            "per_seed": [{"seed": p["seed"], "cons": p["consistency"]} for p in per]}


def main():
    cell_types = G.fetch_cell_types()
    # The cell-types endpoint 504s intermittently and `fetch_cell_types` then falls back to
    # accessibility 1.0, which is ~3x HEK293's real 0.35 and inflates every stage-3 draw -- the
    # first run of this script did exactly that and produced clean sets of 124-198 of 200 against
    # the 9-15 the same cell gives with real data. Refuse to build on the fallback.
    acc = (cell_types or {}).get("HEK293", {}).get("accessibility")
    if not cell_types or acc is None or acc >= 0.99:
        raise SystemExit(f"cell-types unusable (HEK293 accessibility {acc!r}); not building on "
                         f"the accessibility-1.0 fallback")
    G.load_sequence()
    items = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])
    t = next(x for x in items if (x.get("task_id") or x["id"]).startswith(TASK))
    contract, reference = dict(t["content"]["contract"]), t["content"]["hbb_reference"]
    seeds = _parse_seeds(contract["seed"])
    cell = contract["cell_type"]
    base = CJ.config_for(cell)
    c0 = dict(contract, seed=0)
    ctx = G.build_context(c0, reference, cell_types)
    mf = max(1, round(base.cas12a_max_fail * len(SPACE) / 900))

    rng = random.Random(0)
    fill = rng.sample([s for s in range(200, 300) if s not in (A, B)], 6)
    ARMS = [
        ("oracle8", sorted([A, B] + fill), 8),
        ("oracle2", [A, B], 2),
        ("w36/k8", list(range(A, B + 1)), 8),
        ("w36/k9", list(range(A, B + 1)), 9),
        ("w40/k8", list(range(210, 250)), 8),
        ("w100/k8", list(range(200, 300)), 8),
    ]
    print(f"=== {cell} {TASK} seeds {seeds} | cut {len(SPACE)} (200-299 + 400-499) | "
          f"live k={base.band_k} group={base.group_size} | targets {A} and {B} ===\n", flush=True)
    out = []
    for name, cands, k in ARMS:
        cfg = dataclasses.replace(base, seed_list=tuple(SPACE), start_seed=SPACE[0],
                                  end_seed=SPACE[-1], cas12a_max_fail=mf,
                                  band_candidates=tuple(sorted(cands)), band_k=k)
        rows, meta = CJ.build_submission(c0, reference, cell_types, cfg=cfg, budget_s=1800)
        MT.free_gpu_memory()
        if not rows:
            print(f"  {name:9s} cands {len(cands):3d} k={k}  DECLINED  {meta.get('reason')}",
                  flush=True)
            out.append({"arm": name, "k": k, "n_cands": len(cands),
                        "declined": meta.get("reason")})
            continue
        r = price(rows, contract, reference, cell_types, ctx, seeds)
        r.update(arm=name, k=k, n_cands=len(cands), pool=meta.get("pool"),
                 cas9=meta.get("cas9_pool"), band_cfg=sorted(int(x) for x in meta["band_seeds"]))
        mark = "BOTH" if r["both"] else ("one" if (r["has_A"] or r["has_B"]) else "-")
        print(f"  {name:9s} cands {len(cands):3d} k={k}  band {r['band_n']:2d} clean "
              f"{r['clean_n']:3d}  {A}:{'Y' if r['has_A'] else 'n'} {B}:"
              f"{'Y' if r['has_B'] else 'n'}  {mark:4s}  cons {r['round_cons']:.4f} "
              f"final {r['round_final']:6.1f}  pool {r['pool']}", flush=True)
        print(f"            band seeds {r['band_cfg']}", flush=True)
        out.append(r)
        del rows; gc.collect()
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
