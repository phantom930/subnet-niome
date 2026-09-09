#!/usr/bin/env python3
"""band_overlap.py — do stride-10 sibling windows produce independent bands, or correlated ones?

The proposed layout puts 11 hotkeys on width-200 windows at stride 10. Adjacent windows then share
**190 of 200 seeds**, so their Cas12a candidate pools are nearly identical and `FastGreedy` is
deterministic given its bank — the bands could come out highly correlated, collapsing the union
toward a single band and making 11 hotkeys worth little more than one. CLAUDE.md is explicit that a
re-correlated sibling "costs the entire decorrelation", so this has to be measured, not assumed.

Coverage is `|union of bands| / 900`, and what the layout is worth is
`1 - (1 - |union|/900)**3` — the chance a round's three seeds touch some sibling's band.

Two layouts at the same hotkey count and width, so concentration is the only difference:

    concentrated  offsets 0,10,...,100 over a 300-seed span   (the proposal)
    spread        offsets 0,70,...,700 over the full 900      (the control)

Only the Cas12a group is needed to fix the band — the Cas9 half must comply on whatever clean set
the group leaves — so this builds banks and runs the min-union, skipping `scan_cas9`/`assemble`.

    python band_overlap.py [task_id]
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "bandov")

import dataclasses as _dc     # noqa: E402
import itertools              # noqa: E402
import json                   # noqa: E402
import logging                # noqa: E402
import statistics as st       # noqa: E402
import time                   # noqa: E402

logging.basicConfig(level=logging.ERROR)

import genExp as G            # noqa: E402
from niome_subnet.genomics import all_cut as AC      # noqa: E402
from niome_subnet.genomics import all_hdr as AH      # noqa: E402
from niome_subnet.genomics import fastgreedy as FG   # noqa: E402
from niome_subnet.genomics import mt19937 as MT      # noqa: E402
from sd_task import task_content                      # noqa: E402

TASK = _ARGV[1] if len(_ARGV) > 1 else "1eca4bf1-c91f-4291-93ff-4391d058d3dc"
WIDTH = int(os.getenv("BO_WIDTH", "200"))
N_HK = int(os.getenv("BO_N", "11"))
LAYOUTS = [("concentrated stride 10", 100, 10), ("spread stride 70", 100, 70)]


def band_for(window, contract, reference, cell_types, ctx, sites, base, cell):
    """The clean band of the min-union Cas12a group over `window`, as a set of seeds."""
    lo, hi = window
    span = hi - lo + 1
    cfg = _dc.replace(base, hdr_range=(lo, hi),
                      main_max_fail=AH._scaled_max_fail(cell, base.main_max_fail, span),
                      variants=(min(base.variants, AH.WIDE_WINDOW_VARIANTS)
                                if span > 100 else base.variants))
    path = os.path.join(AH.HDR_BANK_DIR, f"cas12a-{AC.bank_key(contract, cell_types, cfg)}.npz")
    if not os.path.exists(path):
        bank = AH.build_bank(contract, reference, cell_types, ctx, sites, cfg)
        if not bank:
            return None, cfg
        AH.save_bank(path, bank)
    recs = AC.load_bank(path)
    if len(recs) < cfg.group_size:
        return None, cfg
    sel = FG.FastGreedy(recs, window_lo=lo, window_hi=hi, per_cell_min=cfg.per_cell_min,
                        caps=AH._group_caps(contract, ctx, cfg))
    idx, _u = sel.best(cfg.group_size, restarts=cfg.restarts)
    bad = set()
    for i in idx:
        bad.update(int(x) for x in recs[i]["fails"])
    MT.free_gpu_memory()
    return set(range(lo, hi + 1)) - bad, cfg


def main():
    task, contract, reference = task_content(TASK)
    cell = contract["cell_type"]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    base = AH.config_for(cell)
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    print(f"task {TASK[:8]}  {cell}  {N_HK} hotkeys x width {WIDTH}  group {base.group_size}\n")
    out = {}
    for label, start, stride in LAYOUTS:
        wins = [(start + h * stride, start + h * stride + WIDTH - 1) for h in range(N_HK)]
        if wins[-1][1] > 999:
            wins = [(lo, min(hi, 999)) for lo, hi in wins]
        print(f"  {label}: {wins[0][0]}-{wins[0][1]} … {wins[-1][0]}-{wins[-1][1]}")
        bands, t0 = [], time.monotonic()
        for i, w in enumerate(wins):
            b, cfg = band_for(w, contract, reference, cell_types, ctx, sites, base, cell)
            bands.append(b)
            print(f"    h{i:<2} {w[0]:>3}-{w[1]:<3} mf {cfg.main_max_fail:>3}  "
                  + (f"band {len(b):>3}  {sorted(b)[:6]}{'…' if len(b) > 6 else ''}"
                     if b else "declined"))
        ok = [b for b in bands if b]
        if len(ok) < 2:
            print("    too few bands to compare\n")
            continue
        union = set().union(*ok)
        pair = [len(a & b) for a, b in itertools.combinations(ok, 2)]
        mean_band = st.mean(len(b) for b in ok)
        # if bands were independent draws inside their own windows, two windows overlapping in
        # `shared` seeds would share about band_a*band_b*shared/(w_a*w_b)
        exp_pair = []
        for (wa, ba), (wb, bb) in itertools.combinations(list(zip(wins, ok)), 2):
            shared = len(set(range(wa[0], wa[1] + 1)) & set(range(wb[0], wb[1] + 1)))
            exp_pair.append(len(ba) * len(bb) * shared /
                            ((wa[1] - wa[0] + 1) * (wb[1] - wb[0] + 1)))
        cov = 1 - (1 - len(union) / 900.0) ** 3
        out[label] = {"windows": wins, "bands": [sorted(b) for b in ok],
                      "band_sizes": [len(b) for b in ok], "union": len(union),
                      "sum_bands": sum(len(b) for b in ok),
                      "mean_pair_overlap": st.mean(pair), "exp_pair_overlap": st.mean(exp_pair),
                      "p_hit": cov, "build_s": round(time.monotonic() - t0, 1)}
        print(f"    -> bands {[len(b) for b in ok]}  sum {sum(len(b) for b in ok)}  "
              f"UNION {len(union)}")
        print(f"       mean pairwise overlap {st.mean(pair):.2f} seeds "
              f"(independent-draw expectation {st.mean(exp_pair):.2f})")
        print(f"       coverage {len(union)}/900 -> P(>=1 of 3 seeds) {cov:.1%}"
              f"   [{time.monotonic()-t0:.0f}s]\n")
    if len(out) == 2:
        a, b = [out[k] for k, _s, _t in [(l, s, t) for l, s, t in LAYOUTS] if k in out][:2]
        print(f"  concentrated union {a['union']} vs spread union {b['union']}  "
              f"-> P(hit) {a['p_hit']:.1%} vs {b['p_hit']:.1%}")
        print(f"  current fleet for reference: 7 hotkeys, ~77 band seeds, 24.3%")
        eff_a = a["union"] / a["sum_bands"]
        eff_b = b["union"] / b["sum_bands"]
        print(f"  band independence (union / sum): concentrated {eff_a:.2f}, spread {eff_b:.2f}"
              f"   (1.00 = perfectly disjoint)")
    json.dump({"task": TASK, "cell": cell, "width": WIDTH, "n": N_HK, "layouts": out},
              open(os.getenv("BO_OUT", "band_overlap.json"), "w"), indent=1)
    print(f"\nwrote {os.getenv('BO_OUT', 'band_overlap.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
