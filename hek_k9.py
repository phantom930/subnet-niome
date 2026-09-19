#!/usr/bin/env python3
"""hek_k9.py — does band_k 9 pay at cut 450 the way it pays at cut 300?

`hek_cutspan.py` established that HEK293's band-formation wall is 9 on the 300-seed and 450-seed cut
spans and 8 on the 900-seed one, so k=9 is reachable on the two narrower arms and not on the shipped
wide one. The entire argument against the wide cut then rests on ONE number that has never been
measured at these spans: CLAUDE.md's band-depth table prices k=9 at **0.000163 against k=8's
0.000136, +20%** -- and that was measured by `band_hit.py` at the NARROW cut, with a different
pricing method, on 12 contracts that are not these four.

So this sweeps k INSIDE each span rather than across spans, which makes the k-step the only variable:

    cut 300 (h0, band ⊆ cut)   k=8 vs k=9     the control -- does +20% reproduce under this method?
    cut 450 (h7, band ⊆ cut)   k=8 vs k=9     the question

Both k=8 arms are rebuilt here rather than spliced in from `hek_cutspan.json`, so every pair is
produced in one process against identical banks and identical sampling. `bank_key` excludes the band
parameters, so k=9 loads the same bank k=8 built and the pairs cost one scan per span, not two.

k=9 sits exactly AT the wall on both spans. That is the same position `band_k` 8 occupies on the
wide cut, and CLAUDE.md's depth table records k=9 +floor building 12/12 there with `band_cell_aware`
on (it is, live), so at-the-wall is expected to build -- but whether it does on 4/4 here is itself
part of the answer, and an arm that declines is scored 0.0 rather than dropped.

    HK9_N=4 python hek_k9.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "hekk9")

import dataclasses                                       # noqa: E402
import json                                              # noqa: E402
import logging                                           # noqa: E402
import time                                              # noqa: E402
from collections import defaultdict                      # noqa: E402

logging.basicConfig(level=logging.ERROR)
import numpy as np                                       # noqa: E402

import genExp as G                                       # noqa: E402
import joined_window as JW                               # noqa: E402
from niome_subnet.genomics import conjunction as CJ      # noqa: E402
from conj_stageb import API                              # noqa: E402
from sd_task import fetch                                # noqa: E402
from widecut_price import SPACE, evaluate, own_fields    # noqa: E402

CELL = "HEK293"
N_CONTRACTS = int(os.getenv("HK9_N", "4"))
KS = [int(x) for x in os.getenv("HK9_KS", "8,9").split(",")]
BUDGET = float(os.getenv("HK9_BUDGET", "1800"))
OUT = os.getenv("HK9_JSON", "hek_k9.json")


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    items = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])
    base = CJ.config_for(CELL)
    ftasks = own_fields(CELL)
    predicted = sorted({s for a, b in JW.fallback_for(CELL, "niome_hotkey")
                        for s in range(a, b + 1)})

    tasks = []
    for t in sorted(items, key=lambda x: x.get("created_at", ""), reverse=True):
        c = (t.get("content") or {}).get("contract") or {}
        s = str(c.get("seed", "") or "")
        if c.get("cell_type") != CELL or len(
                [x for x in s.split(",") if x.strip().isdigit()]) != 3:
            continue
        if (t.get("task_id") or t["id"]) not in ftasks:
            continue
        tasks.append(t)
        if len(tasks) >= N_CONTRACTS:
            break

    # (label, hotkey, cut) -- both spans carry band ⊆ cut, which is what makes them comparable.
    spans = []
    for hk, label in (("niome_hotkey", "cut300"), ("niome_hotkey7", "cut450")):
        bspace = JW.band_space(hk, predicted)
        cands = CJ.sub_window(bspace, base.band_width, JW.band_offset_frac(hk) or 0.0)
        cut = predicted if set(cands) <= set(predicted) else sorted(set(predicted) | set(cands))
        spans.append((label, cands, cut))

    print(f"=== {CELL}  group {base.group_size} width {base.band_width} "
          f"light {base.light_cell_rows} ca {base.band_cell_aware} | k in {KS} "
          f"| {len(tasks)} contracts ===", flush=True)
    for label, cands, cut in spans:
        print(f"    {label}: cut {len(cut)} seeds, band cands {min(cands)}-{max(cands)} "
              f"({len(cands)}), band ⊆ cut = {set(cands) <= set(cut)}", flush=True)

    recs, agg = [], defaultdict(list)
    for t in tasks:
        tid = (t.get("task_id") or t["id"])
        contract, reference = dict(t["content"]["contract"]), t["content"]["hbb_reference"]
        field = ftasks[tid]
        print(f"\n  {tid[:8]}  own cut10 {field[9]:.1f}  top {field[0]:.1f}", flush=True)
        ctx = G.build_context(contract, reference, cell_types)
        for label, cands, cut in spans:
            for k in KS:
                cfg = dataclasses.replace(base, seed_list=tuple(sorted(cut)),
                                          band_candidates=tuple(sorted(cands)), band_k=k)
                key = f"{label}/k{k}"
                t0 = time.monotonic()
                rows, meta = CJ.build_submission(contract, reference, cell_types, cfg=cfg,
                                                 budget_s=BUDGET)
                if not rows:
                    print(f"    {key:12} DECLINED — {meta.get('reason')} "
                          f"({time.monotonic()-t0:.0f}s)", flush=True)
                    agg[key].append(0.0)
                    recs.append({"task": tid[:8], "arm": key, "cut": len(cut), "k": k,
                                 "built": False, "reason": meta.get("reason"), "share": 0.0})
                    continue
                pool, grp = meta.get("pool"), meta.get("group_size")
                r = evaluate(rows, contract, reference, cell_types, ctx, field, key, t0)
                agg[key].append(r["share"] if r else 0.0)
                recs.append(dict(r or {}, task=tid[:8], arm=key, cut=len(cut), k=k, built=True,
                                 pool=pool, margin=pool / grp, cas9=meta.get("cas9_pool"),
                                 bank=meta.get("bank"),
                                 elapsed=round(time.monotonic() - t0, 1), cut10=field[9]))
                print(f"                 pool {pool}/{grp} = {pool/grp:.2f}x  "
                      f"cas9 {meta.get('cas9_pool')}  bank {meta.get('bank')}", flush=True)

    with open(OUT, "w") as fh:
        json.dump(recs, fh, indent=1)

    print(f"\n=== mean E[own-field share], {len(tasks)} contracts, 0.0 for a decline ===",
          flush=True)
    for key in sorted(agg):
        v = agg[key]
        built = sum(1 for x in v if x > 0)
        print(f"  {key:12} {np.mean(v):.6f}  built {built}/{len(v)}   "
              f"{' '.join(f'{x:.6f}' for x in v)}", flush=True)

    print("\n=== the k-step, paired WITHIN span and contract ===", flush=True)
    for label, _c, _cut in spans:
        a, b = f"{label}/k8", f"{label}/k9"
        d = [(x, y) for x, y in zip(agg[a], agg[b])]
        w = sum(1 for x, y in d if y > x)
        lo = sum(1 for x, y in d if y < x)
        ma, mb = float(np.mean([x for x, _ in d])), float(np.mean([y for _, y in d]))
        print(f"  {label}:  k8 {ma:.6f} -> k9 {mb:.6f}  = "
              f"{(mb/ma if ma else float('nan')):.3f}x   {w}W/{lo}L/{len(d)-w-lo}T", flush=True)
    print(f"\nwrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
