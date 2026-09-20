#!/usr/bin/env python3
"""replay_round.py — what would the fleet have scored on one round under the CURRENT config?

Rebuilds all ten hotkeys for a given task with whatever `joined_window` / `conjunction.CELL_CONFIG`
say today, scores each through the validator's own five stages at the round's REAL stamped seeds,
and ranks the result in that round's REAL field. The point is to price a config change against a
round that actually happened rather than against an own-field Monte Carlo.

**Two fidelity details that decide whether this replicates the fleet or merely resembles it.**

* **The build uses `seed: 0`, the scoring uses the stamped seeds.** The miner prefetches hours
  before the seeds are stamped, so its build never sees them -- and `Context.seed` feeds the variant
  RNG, so building against a stamped contract would produce different rows than the fleet could
  ever have shipped. The contract is therefore copied with `seed: 0` for the build and restored for
  the scoring, which is the live sequence exactly.
* **The window comes from the LIVE plan**, not from `fallback_for`, because that is what the miner
  reads. Under `JOINED_SOURCE = "uniform"` that is a constant 100-399 on every cell.

Round final is the mean over the three seeds of `weighted x consistency x fidelity`, which is what
`validation/__init__.py` computes -- weighted and fidelity are seed-independent, so it is equally
`w x fid x mean(cons)`; both are printed so a divergence is visible rather than assumed.

    RR_TASK=cc0458cc python replay_round.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "replay")

import dataclasses                                       # noqa: E402
import json                                              # noqa: E402
import logging                                           # noqa: E402
import time                                              # noqa: E402

logging.basicConfig(level=logging.ERROR)
import numpy as np                                       # noqa: E402

import genExp as G                                       # noqa: E402
import joined_window as JW                               # noqa: E402
from niome_subnet.genomics import conjunction as CJ      # noqa: E402
from niome_subnet.genomics.validation import _parse_seeds  # noqa: E402
from conj_stageb import API, DIST, OURS                  # noqa: E402
from sd_task import fetch, score                         # noqa: E402
from widecut_price import own_fields                     # noqa: E402

TASK = os.getenv("RR_TASK", "cc0458cc")
HKS = JW.BAND_HK + JW.REST_HK
OUT = os.getenv("RR_JSON", "replay_round.json")
# Override `band_k` to replay the same round at another depth. `choose_band` is a greedy PREFIX, so
# a shallower band is the same seeds minus the last picks -- which makes k an unusually clean thing
# to vary: the arms are NESTED rather than merely comparable, and every other input is identical.
K = int(os.getenv("RR_K", "0")) or None


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    items = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])
    task = next(t for t in items if (t.get("task_id") or t["id"]).startswith(TASK))
    tid = task.get("task_id") or task["id"]
    contract = dict(task["content"]["contract"])
    reference = task["content"]["hbb_reference"]
    cell = contract["cell_type"]
    seeds = _parse_seeds(contract["seed"])
    # What the miner actually held at prefetch time.
    build_contract = dict(contract, seed=0)

    plan = json.load(open("data/window_plan.json"))
    raw = plan["assignments"][cell][HKS[0]]
    pairs = raw if isinstance(raw[0], (list, tuple)) else [raw]
    predicted = sorted({s for a, b in pairs for s in range(a, b + 1)})

    field = own_fields(cell).get(tid)
    base = CJ.config_for(cell)
    print(f"=== {tid[:8]}  {cell}  seeds {seeds} ===", flush=True)
    print(f"    plan window {pairs} ({len(predicted)} seeds) | k={K or base.band_k}"
          f"{' (OVERRIDE, config says %d)' % base.band_k if K and K != base.band_k else ''} "
          f"group={base.group_size} width={base.band_width} | cut mode "
          f"{JW.CUT_MODE.get(cell)}", flush=True)
    if field:
        print(f"    real field: {len(field)} miners, top {field[0]:.1f}, "
              f"rank-10 cutoff {field[9]:.1f}", flush=True)
    in_pred = [s for s in seeds if s in set(predicted)]
    print(f"    of the 3 seeds, {len(in_pred)} sit in the predicted window: {in_pred}", flush=True)

    recs = []
    print(f"\n{'hk':4} {'grp':10} {'band window':>13} {'cut':>5} {'band':>5} {'hits':>4} "
          f"{'weighted':>9} {'fid':>6} {'cons/seed':>22} {'final':>7} {'rank':>5} {'share':>6}",
          flush=True)
    for hk in HKS:
        bspace = JW.band_space(hk, predicted)
        cands = CJ.sub_window(bspace, base.band_width, JW.band_offset_frac(hk) or 0.0)
        cut = JW.conjunction_cut_seeds(hk, cell, predicted=predicted, band_candidates=cands)
        cfg = dataclasses.replace(base, seed_list=tuple(sorted(cut)),
                                  band_candidates=tuple(sorted(cands)),
                                  **({"band_k": K} if K else {}))
        t0 = time.monotonic()
        rows, meta = CJ.build_submission(build_contract, reference, cell_types, cfg=cfg,
                                         budget_s=1800)
        grp = "predicted" if hk in JW.BAND_HK else "complement"
        if not rows:
            print(f"{hk[-2:]:4} {grp:10} {'-':>13} {len(cut):>5} DECLINED {meta.get('reason')}",
                  flush=True)
            recs.append({"hk": hk, "group": grp, "built": False,
                         "reason": meta.get("reason"), "final": 0.0})
            continue
        band = sorted(int(x) for x in meta["band_seeds"])
        hits = [s for s in seeds if s in set(band)]
        per = [score(rows, contract, reference, cell_types, seed=s) for s in seeds]
        cons = [p["consistency"] for p in per]
        final = float(np.mean([p["final"] for p in per]))
        alt = per[0]["weighted"] * per[0]["fidelity"] * float(np.mean(cons))
        rank = 1 + sum(1 for x in field if x > final) if field else None
        share = DIST[rank - 1] if (rank and rank <= 10) else 0.0
        print(f"{hk[-2:]:4} {grp:10} {min(cands)}-{max(cands):<8} {len(cut):>5} "
              f"{len(band):>5} {len(hits):>4} {per[0]['weighted']:>9.1f} "
              f"{per[0]['fidelity']:>6.4f} "
              f"{'/'.join(f'{c:.3f}' for c in cons):>22} {final:>7.1f} "
              f"{str(rank):>5} {share:>6.1%}"
              + (f"   BAND HIT {hits}" if hits else "")
              + (f"   [w*fid*mean(cons)={alt:.1f}]" if abs(alt - final) > 0.05 else ""),
              flush=True)
        recs.append({"hk": hk, "group": grp, "built": True, "band": len(band),
                     "band_seeds": band, "hits": hits, "cut": len(cut),
                     "weighted": per[0]["weighted"], "fidelity": per[0]["fidelity"],
                     "cons": cons, "final": final, "rank": rank, "share": share,
                     "clean_in_cut": meta.get("clean"), "pool": meta.get("pool"),
                     "elapsed": round(time.monotonic() - t0, 1)})

    with open(OUT, "w") as fh:
        json.dump({"task": tid, "cell": cell, "seeds": seeds, "window": pairs,
                   "field_top": field[0] if field else None,
                   "field_cut10": field[9] if field else None, "hotkeys": recs}, fh, indent=1)

    built = [r for r in recs if r.get("built")]
    print(f"\n=== fleet result, {len(built)}/{len(recs)} built ===", flush=True)
    if built:
        tot = sum(r["share"] for r in built)
        hitters = [r for r in built if r["hits"]]
        print(f"  band hits: {len(hitters)} hotkey(s)"
              + (f" -> {[(r['hk'][-2:], r['hits']) for r in hitters]}" if hitters else ""),
              flush=True)
        print(f"  best final {max(r['final'] for r in built):.1f} "
              f"(rank {min(r['rank'] for r in built if r['rank'])})", flush=True)
        print(f"  TOTAL curve share {tot:.2%}", flush=True)

    # What the fleet ACTUALLY scored on this round, for the same task, from the score feed.
    sc = fetch(f"{API}/miners/scores?limit=40000", cache="sd_task_scores.json")
    sc = sc if isinstance(sc, list) else (sc.get("data") or sc.get("items") or [])
    ours = {}
    for r in sc:
        if r.get("task_id") != tid or r.get("miner_hotkey") not in OURS:
            continue
        hk, f = r["miner_hotkey"], float(r.get("final_score") or 0)
        ours[hk] = max(ours.get(hk, 0.0), f)
    if ours:
        print(f"\n=== what actually shipped on this round (old config) ===", flush=True)
        for hk, f in sorted(ours.items(), key=lambda kv: -kv[1]):
            rk = 1 + sum(1 for x in field if x > f) if field else None
            print(f"  {hk[:10]}... final {f:7.2f}  rank {rk}", flush=True)
        print(f"  actual total curve share "
              f"{sum(DIST[1 + sum(1 for x in field if x > f) - 1] for f in ours.values() if field and 1 + sum(1 for x in field if x > f) <= 10):.2%}",
              flush=True)
    print(f"\nwrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
