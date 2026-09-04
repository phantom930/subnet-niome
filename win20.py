#!/usr/bin/env python3
"""win20.py — does a narrower screening window buy a wider band *fraction*?

The fleet tiles 100-999 as nine 100-seed windows and each hotkey's clean band is ~13 of its 100.
Coverage is therefore ``n_hotkeys x band``, and the only way to raise it without more hotkeys is a
band that is a larger *fraction* of the window it was searched in.

Theory says that is possible. A band of size B needs ``group_size`` bank guides all HDR-clean on
all B seeds, so with a pool of N enumerated guides B is bounded by ``N x P(HDR)**B >= group_size``
— a bound with no window width in it. The measured 13/12/11/9 at widths 100/150/200/300 is
consistent with that: roughly constant, not proportional. Extrapolated down to width 20 the same
bound permits a band far past 13 of 20, which would be a coverage win.

``main_max_fail`` cannot be carried over, and neither rate-scaling nor the shipped 45 is right:

* 45 of 20 seeds is not a screen at all — every guide passes and the bank is the whole pool.
* Rate-scaling to 9 keeps ~68% of the pool, so the bank_keep truncation decides the bank anyway.

What matters is that the screen not exclude wide-band participants. A guide in a band of size B
fails at most ``width - B`` seeds, so ``max_fail = width - B_target`` is the loosest screen that
excludes nothing a band of >= B_target could use; anything looser only adds guides that cannot
join one. ``max_fail 6`` over 20 seeds admits every guide compatible with a band of 14 or more,
and still passes enough of the pool that ``bank_keep`` (fewest fails first) sets the bank.

Runs the full ``build_submission`` rather than the bank/min-union alone, because a wider band is
not free: ``scan_cas9`` must reach HDR on *every* band seed, which costs ``P(HDR)`` per seed added
(a 29-31 seed band yielded zero Cas9 on all three erythroid types). A band that widens into a Cas9
decline is worth nothing — the build falls through to all-cut and never spikes.
"""
import dataclasses
import json
import sys
import time

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
import logging
logging.basicConfig(level=logging.ERROR)

import genExp as G
from niome_subnet.genomics import all_hdr as AH
from niome_subnet.genomics import mt19937 as MT

TASK = "task-k562.json"
BASE = (300, 399)
NARROW = [(lo, lo + 19) for lo in range(300, 400, 20)]
NARROW_MAX_FAIL = 6
OUT = "win20_results.json"


def run(contract, reference, cell_types, window, max_fail):
    """build_submission's steps, inlined so the band *seeds* are visible, not just the count."""
    from collections import Counter
    import os

    import numpy as np

    from niome_subnet.genomics import fastgreedy as FG
    from niome_subnet.genomics.all_cut import bank_key, load_bank, save_bank

    cfg = dataclasses.replace(AH.config_for(contract["cell_type"]),
                              hdr_range=window, main_max_fail=max_fail)
    started = time.monotonic()
    rec = {"window": list(window), "max_fail": max_fail, "width": window[1] - window[0] + 1}
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    os.makedirs(AH.HDR_BANK_DIR, exist_ok=True)
    path = os.path.join(AH.HDR_BANK_DIR, f"cas12a-{bank_key(contract, cell_types, cfg)}.npz")
    if not os.path.exists(path):
        bank = AH.build_bank(contract, reference, cell_types, ctx, sites, cfg, None)
        MT.free_gpu_memory()
        if not bank:
            return {**rec, "built": False, "reason": "bank scan produced nothing",
                    "wall_s": round(time.monotonic() - started, 1)}
        save_bank(path, bank)
    records = load_bank(path)
    rec["bank"] = len(records)
    rec["bank_min_fails"] = int(min(len(r["fails"]) for r in records))
    if len(records) < cfg.group_size:
        return {**rec, "built": False, "reason": f"bank short of group {cfg.group_size}",
                "wall_s": round(time.monotonic() - started, 1)}

    sel = FG.FastGreedy(records, window_lo=cfg.start_seed, window_hi=cfg.end_seed,
                        per_cell_min=cfg.per_cell_min,
                        caps=AH._group_caps(contract, ctx, cfg))
    index, _union = sel.best(cfg.group_size, restarts=cfg.restarts)
    group = [records[i] for i in index]
    bad = set()
    for r in group:
        bad.update(int(x) for x in r["fails"])
    clean = np.array(sorted(set(range(cfg.start_seed, cfg.end_seed + 1)) - bad), dtype=np.int64)
    rec.update(clean=int(clean.size), band=[int(x) for x in clean],
               clean_fraction=round(clean.size / rec["width"], 4))
    if clean.size == 0:
        return {**rec, "built": False, "reason": "failures cover the whole window",
                "wall_s": round(time.monotonic() - started, 1)}

    n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
    want = n_rows - cfg.group_size
    cas9 = AH.scan_cas9(clean, contract, cell_types, ctx, sites, cfg, want, None)
    MT.free_gpu_memory()
    rec["cas9_pool"] = len(cas9)
    cells4 = len({(r["mutation"], r["strand"]) for r in cas9})
    if len(cas9) < want or cells4 < 4:
        return {**rec, "built": False,
                "reason": f"Cas9 pool {len(cas9)} over {cells4} cells short of {want}",
                "wall_s": round(time.monotonic() - started, 1)}
    rows = AH.assemble(group, cas9, contract, ctx, cfg, n_rows)
    rec.update(rows=len(rows),
               cells=len(Counter((r["mutation"], r["cas_system"], r["strand"]) for r in rows)),
               cas_mix=dict(Counter(r["cas_system"] for r in rows)))
    rec["built"] = rec["rows"] >= n_rows and rec["cells"] >= 8
    if not rec["built"]:
        rec["reason"] = f"assembled {rec['rows']} rows over {rec['cells']} cells"
    rec["wall_s"] = round(time.monotonic() - started, 1)
    return rec


def main():
    task = json.load(open(TASK))
    contract = task["content"]["contract"]
    reference = task["content"]["hbb_reference"]
    seeds = [int(x) for x in str(contract["seed"]).split(",")]
    cell_types = G.fetch_cell_types()
    G.load_sequence()

    jobs = [(BASE, AH.config_for(contract["cell_type"]).main_max_fail)]
    jobs += [(w, NARROW_MAX_FAIL) for w in NARROW]

    out = {"task_id": task["id"], "cell_type": contract["cell_type"], "seeds": seeds,
           "results": []}
    for window, max_fail in jobs:
        rec = run(contract, reference, cell_types, window, max_fail)
        out["results"].append(rec)
        hits = [x for x in seeds if x in rec.get("band", [])]
        print(f"{window[0]}-{window[1]:>3} mf{max_fail:<3} bank {rec.get('bank', 0):>7} "
              f"clean {rec.get('clean', 0):>3}/{rec['width']:<3} "
              f"({rec.get('clean_fraction', 0) * 100:>5.1f}%) "
              f"cas9 {rec.get('cas9_pool', 0):>5} "
              f"rows {rec.get('rows', 0):>3} cells {rec.get('cells', 0)} "
              f"{rec['wall_s']:>6.1f}s  {'built' if rec['built'] else rec.get('reason')}",
              flush=True)
        print(f"      band {rec.get('band')}", flush=True)
        if hits:
            print(f"      HITS drawn seed(s) {hits}", flush=True)
        json.dump(out, open(OUT, "w"), indent=1)

    base = out["results"][0]
    narrow = [r for r in out["results"] if r["width"] == 20]
    cov_base = len(base.get("band", []))
    cov_narrow = sum(len(r.get("band", [])) for r in narrow if r["built"])
    print(f"\ncoverage of 300-399, one hotkey per window:")
    print(f"  1 x width 100 : {cov_base} of 100 seeds   ({base['built'] and 'built' or 'DECLINED'})")
    print(f"  5 x width  20 : {cov_narrow} of 100 seeds "
          f"({sum(1 for r in narrow if r['built'])}/5 built)")
    print(f"  per-hotkey    : {cov_base} vs {cov_narrow / 5:.1f} band seeds")
    out["summary"] = {"base_band": cov_base, "narrow_band_total": cov_narrow,
                      "narrow_built": sum(1 for r in narrow if r["built"])}
    json.dump(out, open(OUT, "w"), indent=1)


if __name__ == "__main__":
    main()
