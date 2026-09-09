#!/usr/bin/env python3
"""gcdist.py — the distance/GC route to `total_weighted_score`, measured against the band it costs.

CLAUDE.md records our weighted gap as the `base structural` half and notes the distance/GC route as
"untested and does not trade against fidelity". This tests it. On task 12be08f9 our eight hotkeys
sat at the 10-14th percentile of the field on weighted while holding the 76-79th on fidelity, and
h7 lost a curve slot by 8.7 points split evenly between weighted and fidelity.

    gc_score   = max(0, 1 - |gc - 0.5|*2)      peaks at GC 0.50, linear falloff
    dist_score = exp(-dist / base_padding)     base_padding 400
    base       = 0.625*gc_score + 0.375*dist_score

`all_hdr` ships cas12a_gc/cas9_gc = (0.40, 0.95) and max_distance 400. A guide at GC 0.95 scores
gc_score 0.10, so the pool admits rows that score badly: measured gc_score 0.938 and dist_score
0.840 give base 0.901 against a ceiling of 1.0.

**The trade this has to clear.** Tightening either bound shrinks the Cas12a candidate pool, and the
min-union turns pool size into BAND size -- which is the payoff. So weighted alone decides nothing;
each arm is scored on the round value that actually places (one band seed of three) and on the
fleet coverage the band buys.

    python gcdist.py [task_id] [lo-hi]
"""
import os
import sys

_ARGV = list(sys.argv)
sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "gcdist")

import dataclasses as _dc      # noqa: E402
import json                    # noqa: E402
import logging                 # noqa: E402
import statistics as st        # noqa: E402
import time                    # noqa: E402

logging.basicConfig(level=logging.ERROR)

import genExp as G             # noqa: E402
from niome_subnet.genomics import all_cut as AC       # noqa: E402
from niome_subnet.genomics import all_hdr as AH       # noqa: E402
from niome_subnet.genomics import mt19937 as MT       # noqa: E402
from niome_subnet.utils import settings                # noqa: E402
from sd_task import score, task_content               # noqa: E402

TASK = _ARGV[1] if len(_ARGV) > 1 else "12be08f9-8e88-455d-aa18-3a995626df20"
WINDOW = _ARGV[2] if len(_ARGV) > 2 else "100-399"
FLOOR = 0.10          # measured off-band consistency, 0.097-0.106

# (label, gc bounds, max_distance). Baseline first so every later arm is read against it.
ARMS_ALL = [
    ("baseline  (0.40-0.95, d400)", (0.40, 0.95), 400),
    ("tight GC  (0.45-0.55, d400)", (0.45, 0.55), 400),
    ("mid GC    (0.42-0.58, d400)", (0.42, 0.58), 400),
    ("tight d   (0.40-0.95, d150)", (0.40, 0.95), 150),
    ("both      (0.45-0.55, d150)", (0.45, 0.55), 150),
    ("both mid  (0.42-0.58, d250)", (0.42, 0.58), 250),
]
# GCD_ARMS selects a subset by substring, so a verification run can carry just the two arms that
# matter without rebuilding the four that were only there to locate the optimum.
_want = [x.strip() for x in (os.getenv("GCD_ARMS") or "").split(",") if x.strip()]
ARMS = ([a for a in ARMS_ALL if any(w in a[0] for w in _want)] if _want else ARMS_ALL)
OUT = os.getenv("GCD_OUT", "gcdist.json")


def main():
    task, contract, reference = task_content(TASK)
    cell = contract["cell_type"]
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    lo, hi = (int(x) for x in WINDOW.split("-"))
    span = hi - lo + 1
    base = AH.config_for(cell)
    mf = AH._scaled_max_fail(cell, base.main_max_fail, span)
    print(f"task {TASK[:8]}  {cell}  window {lo}-{hi} (span {span})  mf {mf}  "
          f"group {base.group_size}")
    print(f"base_padding {contract['rules'].get('base_padding')}  "
          f"k=1 round consistency = (1 + 2*{FLOOR})/3 = {(1+2*FLOOR)/3:.4f}\n")
    cons_k1 = (1 + 2 * FLOOR) / 3

    print(f"  {'arm':<28} {'band':>5} {'cas9':>6} {'weighted':>9} {'gc_s':>6} {'dist_s':>7} "
          f"{'fid':>7} {'k=1 final':>10} {'P(hit)':>7} {'fleet':>7} {'s':>5}")
    out, ref = [], None
    for label, gc, dist in ARMS:
        cfg = _dc.replace(base, hdr_range=(lo, hi), main_max_fail=mf,
                          variants=min(base.variants, AH.WIDE_WINDOW_VARIANTS),
                          cas12a_gc=gc, cas9_gc=gc, max_distance=dist)
        p = os.path.join(AH.HDR_BANK_DIR, f"cas12a-{AC.bank_key(contract, cell_types, cfg)}.npz")
        if os.path.exists(p):
            os.remove(p)
        t0 = time.monotonic()
        try:
            rows, meta = AH.build_submission(contract, reference, cell_types, cfg=cfg,
                                             budget_s=700)
        except Exception as exc:
            print(f"  {label:<28}  raised {type(exc).__name__}: {exc}")
            out.append({"arm": label, "reason": f"{type(exc).__name__}: {exc}"})
            MT.free_gpu_memory()
            continue
        dt = time.monotonic() - t0
        MT.free_gpu_memory()
        if not rows:
            print(f"  {label:<28} {'-':>5} {str(meta.get('cas9_pool') or '-'):>6}"
                  f"   declined: {meta.get('reason')}")
            out.append({"arm": label, "reason": meta.get("reason"),
                        "band": meta.get("clean"), "cas9": meta.get("cas9_pool")})
            continue
        # Score once on a band seed for weighted/fidelity, and read the structural terms back
        # from stage 12's own per-row output rather than recomputing them.
        s = score(rows, contract, reference, cell_types, seed=lo)
        detail = json.load(open(settings.VALID_EXPERIMENTS_PATH))
        # gc_score/dist_score are per-row under "features"; structural_score is stage 2's own
        # 0.625*gc + 0.375*dist, so reading it back verifies the arm rather than assuming it.
        feats = [d["features"] for d in detail]
        gcs = st.mean(f["gc_score"] for f in feats)
        ds = st.mean(f["dist_score"] for f in feats)
        struct = st.mean(d["stage2"]["structural_score"] for d in detail)
        offt = st.mean(f["offtarget_factor"] for f in feats)
        gc_mean = st.mean(f["gc"] for f in feats)
        dist_mean = st.mean(f["distance_to_mutation"] for f in feats)
        band = meta["clean"]
        p_hit = 1 - (1 - band / 900.0) ** 3
        fleet = 1 - (1 - 7 * band / 900.0) ** 3
        k1 = s["weighted"] * cons_k1 * s["fidelity"]
        rec = {"arm": label, "gc": gc, "max_distance": dist, "band": band,
               "cas9": meta.get("cas9_pool"), "weighted": s["weighted"],
               "gc_score": gcs, "dist_score": ds, "structural": struct,
               "offtarget": offt, "gc_mean": gc_mean, "dist_mean": dist_mean,
               "fidelity": s["fidelity"],
               "k1_final": k1, "p_hit": p_hit, "fleet_cov": fleet, "build_s": round(dt, 1)}
        out.append(rec)
        if ref is None:
            ref = rec
            print(f"      (baseline row detail: mean GC {gc_mean:.3f}, mean distance "
                  f"{dist_mean:.0f}bp, structural {struct:.4f}, offtarget {offt:.4f})")
        print(f"  {label:<28} {band:>5} {meta.get('cas9_pool'):>6} {s['weighted']:>9.2f} "
              f"{gcs:>6.3f} {ds:>7.3f} {s['fidelity']:>7.4f} {k1:>10.2f} {p_hit:>6.2%} "
              f"{fleet:>6.1%} {dt:>5.0f}")

    good = [r for r in out if r.get("band")]
    if ref and len(good) > 1:
        print(f"\n  against baseline (weighted {ref['weighted']:.1f}, band {ref['band']}, "
              f"k=1 {ref['k1_final']:.1f}):")
        for r in good[1:]:
            dw = (r["weighted"] / ref["weighted"] - 1) * 100
            dk = (r["k1_final"] / ref["k1_final"] - 1) * 100
            db = r["band"] - ref["band"]
            print(f"    {r['arm']:<28} weighted {dw:+6.1f}%  band {db:+3d}  "
                  f"k=1 final {dk:+6.1f}%  fleet cov {r['fleet_cov']-ref['fleet_cov']:+.1%}")
        best = max(good, key=lambda r: r["k1_final"])
        print(f"\n  best k=1 round score: {best['arm']}  ({best['k1_final']:.2f})")
        print(f"  NOTE k=1 final is what places -- rank-10 cutoffs run 45.9-135.1 on K562, "
              f"and this round's was 99.08.")
    json.dump({"task": TASK, "cell": cell, "window": [lo, hi], "arms": out},
              open(OUT, "w"), indent=1)
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
