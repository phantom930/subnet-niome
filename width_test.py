#!/usr/bin/env python3
"""width_test.py — band window width vs weighted/fidelity, on one contract.

On HEK293 round a8b9f1bb our band hotkeys held weighted 163-167 and fidelity 0.90 while the field
ran 221 and 0.962 -- and ranks 7-10 did it at cons 0.386, the same single-band-hit signature as our
h3. So the field is getting a band hit AND good structural terms; we are not.

Hypothesis: the width-300 layout is the cause. A wider window makes a guide's fail count a tighter
estimate of its true P(HDR), so selecting the extreme tail selects harder on P(HDR) -- which on
HEK293 IS GC, because accessibility 0.35 leaves `energy` under the clamp (measured: P(HDR) rises
0.346 -> 0.416 across GC 0.30 -> 0.75, where K562 is flat above GC 0.40). At width 100 the draw is
noisier, so low-GC guides survive the cut and gc_score stays high.

CLAUDE.md records HEK293 all-HDR at width 100 with weighted 318 / fidelity 0.924 against our 166 /
0.90 at width 300, but on a different contract -- this measures both on the same one.

Coverage is reported for a 7-hotkey fleet at each width, since a narrower window means a smaller
band and the comparison is frequency x value, not weighted alone.

    python width_test.py [task_id]
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "wtest")

import dataclasses as _dc     # noqa: E402
import json                   # noqa: E402
import logging                # noqa: E402
import statistics as st       # noqa: E402
import time                   # noqa: E402

logging.basicConfig(level=logging.ERROR)

import genExp as G            # noqa: E402
from niome_subnet.genomics import all_cut as AC      # noqa: E402
from niome_subnet.genomics import all_hdr as AH      # noqa: E402
from niome_subnet.genomics import mt19937 as MT      # noqa: E402
from niome_subnet.utils import settings              # noqa: E402
from sd_task import score, task_content              # noqa: E402

TASK = _ARGV[1] if len(_ARGV) > 1 else "a8b9f1bb-f25d-42af-92a1-e0aff16bee3a"
FLOOR = 0.10
N_HOTKEYS = 7
# Two width-100 windows and two width-300, so a difference cannot be blamed on window position.
WINDOWS = [(300, 399), (600, 699), (300, 599), (100, 399)]


def main():
    task, contract, reference = task_content(TASK)
    cell = contract["cell_type"]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    base = AH.config_for(cell)
    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    cons_k1 = (1 + 2 * FLOOR) / 3
    print(f"task {TASK[:8]}  {cell}  accessibility "
          f"{cell_types.get(cell, {}).get('accessibility')}  group {base.group_size}")
    print(f"k=1 round consistency {cons_k1:.4f}; fleet = {N_HOTKEYS} band hotkeys\n")
    print(f"  {'window':<10} {'span':>5} {'mf':>4} {'var':>6} {'band':>5} {'cas9':>6} "
          f"{'weighted':>9} {'gc_s':>6} {'dist_s':>7} {'fid':>7} {'k=1 final':>10} "
          f"{'fleet cov':>10} {'freq x val':>11} {'s':>5}")
    out = []
    for lo, hi in WINDOWS:
        span = hi - lo + 1
        cfg = _dc.replace(base, hdr_range=(lo, hi),
                          main_max_fail=AH._scaled_max_fail(cell, base.main_max_fail, span),
                          variants=(min(base.variants, AH.WIDE_WINDOW_VARIANTS)
                                    if span > 100 else base.variants))
        p = os.path.join(AH.HDR_BANK_DIR, f"cas12a-{AC.bank_key(contract, cell_types, cfg)}.npz")
        if os.path.exists(p):
            os.remove(p)
        t0 = time.monotonic()
        rows, meta = AH.build_submission(contract, reference, cell_types, cfg=cfg, budget_s=800)
        dt = time.monotonic() - t0
        MT.free_gpu_memory()
        if not rows:
            print(f"  {f'{lo}-{hi}':<10} {span:>5} {cfg.main_max_fail:>4} {cfg.variants:>6}"
                  f"   declined: {meta.get('reason')}")
            out.append({"window": [lo, hi], "reason": meta.get("reason")})
            continue
        s = score(rows, contract, reference, cell_types, seed=lo)
        feats = [d["features"] for d in json.load(open(settings.VALID_EXPERIMENTS_PATH))]
        band = meta["clean"]
        cov = 1 - (1 - N_HOTKEYS * band / 900.0) ** 3
        k1 = s["weighted"] * cons_k1 * s["fidelity"]
        rec = {"window": [lo, hi], "span": span, "mf": cfg.main_max_fail,
               "variants": cfg.variants, "band": band, "cas9": meta.get("cas9_pool"),
               "weighted": s["weighted"], "fidelity": s["fidelity"],
               "gc_score": st.mean(f["gc_score"] for f in feats),
               "dist_score": st.mean(f["dist_score"] for f in feats),
               "gc_mean": st.mean(f["gc"] for f in feats),
               "k1_final": k1, "fleet_cov": cov, "freq_x_val": cov * k1,
               "build_s": round(dt, 1)}
        out.append(rec)
        print(f"  {f'{lo}-{hi}':<10} {span:>5} {cfg.main_max_fail:>4} {cfg.variants:>6} "
              f"{band:>5} {meta.get('cas9_pool'):>6} {s['weighted']:>9.2f} "
              f"{rec['gc_score']:>6.3f} {rec['dist_score']:>7.3f} {s['fidelity']:>7.4f} "
              f"{k1:>10.2f} {cov:>9.1%} {rec['freq_x_val']:>11.2f} {dt:>5.0f}")

    good = [r for r in out if r.get("band")]
    if good:
        w100 = [r for r in good if r["span"] == 100]
        w300 = [r for r in good if r["span"] == 300]
        print()
        if w100 and w300:
            def m(rs, k): return st.mean(r[k] for r in rs)
            print(f"  width 100 (n={len(w100)}): band {m(w100,'band'):.1f}  "
                  f"weighted {m(w100,'weighted'):.1f}  gc_s {m(w100,'gc_score'):.3f}  "
                  f"fid {m(w100,'fidelity'):.4f}  k=1 {m(w100,'k1_final'):.1f}  "
                  f"cov {m(w100,'fleet_cov'):.1%}  freq x val {m(w100,'freq_x_val'):.2f}")
            print(f"  width 300 (n={len(w300)}): band {m(w300,'band'):.1f}  "
                  f"weighted {m(w300,'weighted'):.1f}  gc_s {m(w300,'gc_score'):.3f}  "
                  f"fid {m(w300,'fidelity'):.4f}  k=1 {m(w300,'k1_final'):.1f}  "
                  f"cov {m(w300,'fleet_cov'):.1%}  freq x val {m(w300,'freq_x_val'):.2f}")
            r = m(w100, 'freq_x_val') / m(w300, 'freq_x_val')
            print(f"\n  width 100 is {r:.2f}x width 300 on frequency x value  "
                  f"-> {'WIDTH 100' if r > 1 else 'WIDTH 300'} for {cell}")
        print(f"\n  this round's rank-10 cutoff was 81.92; a k=1 final below that never places.")
    json.dump({"task": TASK, "cell": cell, "arms": out},
              open(os.getenv("WT_OUT", "width_test.json"), "w"), indent=1)
    print(f"\nwrote {os.getenv('WT_OUT', 'width_test.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
