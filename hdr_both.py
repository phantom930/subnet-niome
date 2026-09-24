#!/usr/bin/env python3
"""hdr_both.py — can ONE all-HDR hotkey band both 214 and 249 on f18ee409?

`both_seeds.py` answered this for the conjunction: no. Its band is exactly `band_k` seeds chosen
by a greedy whose picks depend on which guides survive the cut min-union, and across widths 36/40
/100 and k=8/9 it never took both -- it took 214 once, 249 never, and nothing twice.

all-HDR is a different mechanism and has two levers the conjunction lacks:

  * **the band is EMERGENT, not a fixed k.** It is `band_space` minus the union of the group's HDR
    failures, so it GROWS toward the window as the window shrinks -- CLAUDE.md records band 10 of
    10 at width 10 on HEK293, against 8 at width 40 (M16) and 7 at width 100.
  * **`group_size` sets the failure union.** Fewer guides fail on fewer seeds, so a smaller group
    leaves a wider clean band -- measured 14.6 at group 42 against 12.1 at group 100 on the
    erythroid cells. It is paid for in stage-5 cas-coverage entropy.

and one more that matters here: `seed_list` takes an explicit NON-CONTIGUOUS band space, so the
two targets can be bracketed by two narrow blocks instead of one 36-seed span. 214 and 249 are 35
apart, which forces any contiguous window to 36; a joined space reaches both in 12-20 seeds.

Arms (HEK293 live config otherwise: group 80, main_max_fail scaled to the span):

    oracle2    seed_list = {214, 249}                    -- the band space IS the targets
    w36/g80    214..249            (36 seeds)  group 80  -- narrowest contiguous, live group
    w36/g42    214..249            (36)        group 42  -- smaller group, wider band
    w36/g20    214..249            (36)        group 20
    j20/g80    210-219 + 245-254   (20)        group 80  -- joined, brackets both
    j12/g80    212-217 + 247-252   (12)        group 80  -- joined, tighter
    j12/g42    212-217 + 247-252   (12)        group 42

`oracle2` needs the stamped seeds and is not buildable blind; it prices the ceiling. The joined
arms need only to know roughly WHERE the seeds are -- window-level knowledge, which is the
resolution the seed model actually operates at -- so they are the interesting ones.
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "hdrboth")

import dataclasses, gc, json, logging                                 # noqa: E402

logging.basicConfig(level=logging.ERROR)
import numpy as np                                                    # noqa: E402

import genExp as G                                                    # noqa: E402
from niome_subnet.genomics import all_hdr as AH                       # noqa: E402
from niome_subnet.genomics import conjunction as CJ                   # noqa: E402
from niome_subnet.genomics import mt19937 as MT                       # noqa: E402
from niome_subnet.genomics.validation import _parse_seeds             # noqa: E402
from conj_stageb import API                                           # noqa: E402
from sd_task import fetch, score                                      # noqa: E402
from widecut_price import records_of                                  # noqa: E402

TASK = os.getenv("HB_TASK", "f18ee409")
OUT = os.getenv("HB_JSON", "hdr_both.json")
SPACE = sorted(list(range(200, 300)) + list(range(400, 500)))
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


def main():
    cell_types = G.fetch_cell_types()
    acc = (cell_types or {}).get("HEK293", {}).get("accessibility")
    if not cell_types or acc is None or acc >= 0.99:
        raise SystemExit(f"cell-types unusable (HEK293 accessibility {acc!r})")
    G.load_sequence()
    items = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])
    t = next(x for x in items if (x.get("task_id") or x["id"]).startswith(TASK))
    contract, reference = dict(t["content"]["contract"]), t["content"]["hbb_reference"]
    seeds = _parse_seeds(contract["seed"])
    cell = contract["cell_type"]
    base = AH.config_for(cell)
    c0 = dict(contract, seed=0)
    ctx = G.build_context(c0, reference, cell_types)

    j20 = list(range(210, 220)) + list(range(245, 255))
    j12 = list(range(212, 218)) + list(range(247, 253))
    ARMS = [("oracle2", [A, B], 80),
            ("w36/g80", list(range(A, B + 1)), 80),
            ("w36/g42", list(range(A, B + 1)), 42),
            ("w36/g20", list(range(A, B + 1)), 20),
            ("j20/g80", j20, 80),
            ("j12/g80", j12, 80),
            ("j12/g42", j12, 42)]
    print(f"=== all-HDR | {cell} {TASK} seeds {seeds} | accessibility {acc} | "
          f"live group {base.group_size} | targets {A} and {B} ===\n", flush=True)
    out = []
    for name, sl, grp in ARMS:
        sl = sorted(set(int(x) for x in sl))
        span = len(sl)
        cfg = dataclasses.replace(
            base, hdr_range=(sl[0], sl[-1]), seed_list=tuple(sl), group_size=grp,
            main_max_fail=AH._scaled_max_fail(cell, base.main_max_fail, span,
                                              AH._native_span(cell)))
        rows, meta = AH.build_submission(c0, reference, cell_types, cfg=cfg, budget_s=1800)
        MT.free_gpu_memory()
        if not rows:
            print(f"  {name:8s} space {span:3d} group {grp:3d}  DECLINED  {meta.get('reason')}",
                  flush=True)
            out.append({"arm": name, "span": span, "group": grp,
                        "declined": meta.get("reason")})
            continue
        band = compliant(rows, contract, cell_types, ctx, "hdr")
        cut = compliant(rows, contract, cell_types, ctx, "cut")
        bs = set(band); cs = set(cut) - bs
        per = [score(rows, contract, reference, cell_types, seed=s) for s in seeds]
        cons = float(np.mean([p["consistency"] for p in per]))
        wxf = per[0]["weighted"] * per[0]["fidelity"]
        both = (A in bs and B in bs)
        mark = "BOTH" if both else ("one" if (A in bs or B in bs) else "-")
        print(f"  {name:8s} space {span:3d} group {grp:3d}  band {len(band):3d} "
              f"({len(band)/span:4.0%} of space)  clean {len(cs):3d}  {A}:{'Y' if A in bs else 'n'} "
              f"{B}:{'Y' if B in bs else 'n'}  {mark:4s}  cons {cons:.4f} final {wxf*cons:6.1f}  "
              f"cfg-clean {meta.get('clean')}  cells {meta.get('cells')}/8", flush=True)
        print(f"           band {band}", flush=True)
        out.append({"arm": name, "span": span, "group": grp, "band_n": len(band), "band": band,
                    "clean_n": len(cs), "has_A": A in bs, "has_B": B in bs, "both": both,
                    "round_cons": cons, "weighted": per[0]["weighted"],
                    "fidelity": per[0]["fidelity"], "round_final": wxf * cons,
                    "cells": meta.get("cells"), "cas9": meta.get("cas9_pool")})
        del rows; gc.collect()
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
