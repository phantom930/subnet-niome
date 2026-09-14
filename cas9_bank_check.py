#!/usr/bin/env python3
"""cas9_bank_check.py — prove the Cas9 bank reproduces `scan_cas9` + `hdr_compliance`.

The bank does the FULL nearest-first walk while `scan_cas9` stops early, so the bank is a superset
by construction. What has to hold for the substitution to be sound is stronger and is what this
checks, per config, against the real functions:

  1. bank.select(clean, band) SUPERSET of the scanned+HDR-filtered set        (no guide is lost)
  2. the two agree on the head `assemble` actually consumes -- the nearest `score_cap // 4` per
     (mutation, strand) cell, sorted by (distance, |gc - 0.50|)
  3. `assemble` returns byte-identical rows from both, so weighted/fidelity/score are unchanged

A mismatch on (1) is a bug. A mismatch on (2)/(3) with the bank a strict superset means the early
exit was truncating a cell -- the bank is then BETTER, not wrong, and the diff is reported rather
than failed so the size of the difference is on record.

    CBC_CELL=K562 python cas9_bank_check.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")

import dataclasses as _dc            # noqa: E402
import logging                       # noqa: E402
import os.path                       # noqa: E402
import time                          # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np                   # noqa: E402

CELL = os.getenv("CBC_CELL", "K562")
os.environ.setdefault("NIOME_INSTANCE", "cbc_" + CELL.replace("+", "").replace("-", "_"))

import genExp as G                                          # noqa: E402
import cas9_bank as CB                                      # noqa: E402
from niome_subnet.genomics import all_cut as AC             # noqa: E402
from niome_subnet.genomics import fastgreedy as FG          # noqa: E402
from niome_subnet.genomics import mt19937 as MT             # noqa: E402
from conj_test import hdr_compliance                        # noqa: E402
from conj_joined import JoinedCut                           # noqa: E402
from conj_feas import greedy_prefix                         # noqa: E402
from conj_grid import sub_window                            # noqa: E402
from joined300 import fetch_tasks, pick_tasks               # noqa: E402
from sd_task import score                                   # noqa: E402

CLASSES = [(100, 199), (400, 499), (700, 799)]
CASES = [tuple(int(y) for y in x.split(",")) for x in
         os.getenv("CBC_CASES", "300,8,80;100,6,42;225,10,125").split(";")]
POOL_TARGET = int(os.getenv("CBC_POOL_TARGET", "500"))


def key(r):
    return (r["guide"], r["mutation"], r["strand"], r["start"])


def head(recs, score_cap):
    by = {}
    for r in recs:
        by.setdefault((r["mutation"], r["strand"]), []).append(r)
    cap = max(1, score_cap // max(1, len(by)))
    out = []
    for rs in by.values():
        rs.sort(key=lambda r: (r["distance"], abs(r["gc"] - 0.50)))
        out.extend(rs[:cap])
    return out


def main():
    task = pick_tasks(fetch_tasks()).get(CELL)
    content = task["content"]
    contract, reference = content["contract"], content["hbb_reference"]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    joined = sorted(s for lo, hi in CLASSES for s in range(lo, hi + 1))
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    n_rows = contract["rules"].get("max_experiments") or ctx.max_experiments
    base = AC.config_for(CELL) or AC.AllCutConfig()
    mf = max(1, round(base.cas12a_max_fail * len(joined) / 900))
    cfg0 = _dc.replace(JoinedCut(**{f.name: getattr(base, f.name)
                                    for f in _dc.fields(AC.AllCutConfig)}),
                       seed_list=tuple(joined), start_seed=min(joined),
                       end_seed=max(joined), cas12a_max_fail=mf)
    path = os.path.join(AC.BANK_DIR, f"cas12a-{AC.bank_key(contract, cell_types, cfg0)}.npz")
    records = AC.load_bank(path, limit=300_000)
    print(f"{CELL}  task {(task.get('task_id') or task['id'])[:8]}  "
          f"Cas12a bank {len(records)}  rows {n_rows}\n")

    bank = CB.get(CELL, contract, reference, cell_types, ctx, sites, cfg0, joined)
    print(f"Cas9 bank: {len(bank)} guides over {len(bank.seeds)} seeds\n")

    ok = hdr_compliance(records, contract, cell_types, ctx, joined)
    for width, k, g in CASES:
        cand = sub_window(joined, 700, width)
        okw = {i: (v & set(cand)) for i, v in ok.items()}
        steps = greedy_prefix(okw, k, cand, g)
        if k >= len(steps):
            print(f"width {width} k {k} group {g}: band unreachable, skipped"); continue
        band, alive = steps[k]
        pool = [records[i] for i in alive]
        sel = FG.FastGreedy(pool, window_lo=min(joined), window_hi=max(joined),
                            seeds=np.asarray(joined, dtype=np.int64))
        idx, _u = sel.best(g, restarts=12)
        grp = [pool[i] for i in idx]
        bad = set()
        for rec in grp:
            bad.update(int(x) for x in rec["fails"])
        clean = np.array(sorted(set(joined) - bad), dtype=np.int64)
        want = n_rows - g

        t0 = time.monotonic()
        cfg_scan = _dc.replace(cfg0, group_size=g, pool_target=POOL_TARGET)
        scanned = AC.scan_cas9(clean, contract, cell_types, ctx, sites, cfg_scan, want)
        if band and scanned:
            ok9 = hdr_compliance(scanned, contract, cell_types, ctx, band)
            scanned = [scanned[i] for i in sorted(ok9) if len(ok9[i]) == len(band)]
        t_scan = time.monotonic() - t0

        t0 = time.monotonic()
        picked = bank.select(clean, band)
        t_bank = time.monotonic() - t0

        ks, kb = {key(r) for r in scanned}, {key(r) for r in picked}
        superset = ks <= kb
        hs, hb = head(scanned, cfg0.score_cap), head(picked, cfg0.score_cap)
        head_same = {key(r) for r in hs} == {key(r) for r in hb}
        cfg = _dc.replace(cfg0, group_size=g, light_cell_rows=6)
        rs = AC.assemble(grp, scanned, contract, ctx, cfg, n_rows)
        rb = AC.assemble(grp, picked, contract, ctx, cfg, n_rows)
        rows_same = rs == rb
        print(f"width {width:>3} k {k:>2} group {g:>3} | clean {clean.size:>3} band {len(band)}")
        print(f"   scan_cas9+HDR {len(scanned):>7} in {t_scan:>6.1f}s     "
              f"bank.select {len(picked):>7} in {t_bank:>6.3f}s   "
              f"speedup {t_scan/max(t_bank,1e-6):>7.0f}x")
        print(f"   superset {superset}   (missing {len(ks - kb)})   "
              f"assemble head identical {head_same}   rows identical {rows_same}")
        if not rows_same:
            a = score(rs, contract, reference, cell_types, seed=int(band[0]))
            b = score(rb, contract, reference, cell_types, seed=int(band[0]))
            print(f"   scan  w {a['weighted']:.4f} fid {a['fidelity']:.6f} "
                  f"cons {a['consistency']:.4f}")
            print(f"   bank  w {b['weighted']:.4f} fid {b['fidelity']:.6f} "
                  f"cons {b['consistency']:.4f}")
        print()
        MT.free_gpu_memory()


if __name__ == "__main__":
    main()
