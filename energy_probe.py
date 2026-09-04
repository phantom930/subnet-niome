#!/usr/bin/env python3
"""energy_probe.py — can feature-outcome correlation lift consistency? Per-target R2, measured.

The proposal: make guide features correlate with outcomes so stage 4's forest can predict them,
rather than pinning outcomes constant on a band. stage 3 says how far that can go:

    energy = clamp(accessibility * (1.8*gc + 0.6*exp(-dist/1500) + region_offset), 0, 1)
    cut_p  = min(0.99, max(0.4, base + 0.18*energy))      base 0.86 Cas9 / 0.78 Cas12a
    P(HDR) = (hdr_base + 0.35*energy) / (that + mh_nhej + blunt)

Energy is the only lever a designer controls, through gc and distance. But cut_p spans just
[0.78, 0.99] over the whole energy range and P(HDR) just [0.40, 0.59], so the *predictable* share
of outcome variance is small by construction: for is_hdr, Var(p) ~ 0.003 against E[p(1-p)] ~ 0.245,
an R2 ceiling near 0.01.

This measures it rather than trusting that algebra. Four row sets on the same contract, chosen only
by where they sit in energy, each scored on three seeds with the per-target R2 and normalised MAE
pulled out of stage 4:

  high     maximum energy   -- cut_p and P(HDR) both at ceiling, near-zero spread
  low      minimum energy   -- both at floor, near-zero spread
  bimodal  half at each end -- maximises Var(energy), the best case for learnability
  closed   minimum distance + k-mer diverse -- the "close the cut site, keep diversity" proposal

consistency = 0.7 * max(avg_r2, 0) + 0.3 * (1 - avg_nmae), so a design can only win here by
driving avg_r2 positive. If bimodal cannot, no arrangement of guides can.
"""
import os

os.environ["NIOME_INSTANCE"] = "win20"

import json
import math
import statistics as st
import sys
from collections import defaultdict

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
import logging
logging.basicConfig(level=logging.ERROR)

import genExp as G
from niome_subnet.genomics import seed_agnostic as SA
from niome_subnet.genomics.validation import run_stage12, run_stage3, run_stage4, run_stage5
from niome_subnet.genomics.validation import stage5
from niome_subnet.utils import settings

TASKS = ["task-k562.json"]
SEEDS = [271]
GC_BAND = (0.15, 0.95)          # deliberately wider than any shipped builder
MAX_DISTANCE = 4000             # deliberately far, to reach the low-energy end
VARIANTS = 900
PER_BUCKET = 220
OUT = "energy_probe_r2.json"


def energy_of(f):
    return max(0.0, min(1.0, f.get("cell_type_accessibility", 1.0) * (
        1.8 * f["gc"] + 0.6 * math.exp(-f["distance_to_mutation"] / 1500)
        + f.get("region_energy_offset", 0.0))))


def candidates(contract, reference, cell_types):
    c = dict(contract)
    c["seed"] = 500
    ctx = G.build_context(c, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    jobs = [(s, m) for s in sites for m in ctx.mutations
            if abs(s.start - ctx.mutation_map[m]) <= MAX_DISTANCE]
    buckets = defaultdict(list)
    for site, mutation in jobs:
        for guide in SA.enumerate_variants(site, ctx, GC_BAND[0], GC_BAND[1],
                                           ctx.max_mismatches, True, VARIANTS):
            entry = G.build_valid_entry(G.make_experiment(site, guide, mutation, ctx, "c"), ctx)
            if entry is None:
                continue
            f = entry["features"]
            e = energy_of(f)
            key = ((mutation, site.cas, site.strand), min(9, int(e * 10)))
            if len(buckets[key]) >= PER_BUCKET:
                continue
            buckets[key].append({"guide": guide, "mutation": mutation, "cas_system": site.cas,
                                 "strand": site.strand, "start": site.start,
                                 "length": site.length, "energy": e, "gc": f["gc"],
                                 "distance": f["distance_to_mutation"],
                                 "weighted_score": entry["stage2"]["weighted_score"],
                                 "kmers": frozenset(stage5.extract_kmers(guide, 12))})
    by_cell = defaultdict(list)
    for (cell, _e), recs in buckets.items():
        by_cell[cell] += recs
    return by_cell


def pick(by_cell, n_rows, mode):
    cells = sorted(by_cell)
    if len(cells) < 8:
        return None
    out = []
    base, extra = divmod(n_rows, len(cells))
    for i, key in enumerate(cells):
        take = base + (1 if i < extra else 0)
        recs = by_cell[key]
        if mode == "high":
            sel = sorted(recs, key=lambda r: -r["energy"])[:take]
        elif mode == "low":
            sel = sorted(recs, key=lambda r: r["energy"])[:take]
        elif mode == "bimodal":
            hi = sorted(recs, key=lambda r: -r["energy"])
            lo = sorted(recs, key=lambda r: r["energy"])
            sel = hi[:take // 2] + lo[:take - take // 2]
        else:                                   # closed distance + k-mer diversity
            near = sorted(recs, key=lambda r: r["distance"])[:max(take * 6, 200)]
            used, sel = set(), []
            for _ in range(take):
                best = max((r for r in near if r not in sel),
                           key=lambda r: len(r["kmers"] - used), default=None)
                if best is None:
                    break
                sel.append(best)
                used |= best["kmers"]
        out += sel
    return out[:n_rows] if len(out) >= n_rows else None


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    out = []
    for tf in TASKS:
        task = json.load(open(tf))
        contract = dict(task["content"]["contract"])
        reference = task["content"]["hbb_reference"]
        n_rows = contract["rules"].get("max_experiments") or 250
        print(f"=== {tf}  {contract['cell_type']} ===", flush=True)
        by = candidates(contract, reference, cell_types)
        print(f"  {sum(len(v) for v in by.values())} candidates over {len(by)} cells; "
              f"energy range "
              f"{min(r['energy'] for v in by.values() for r in v):.3f}-"
              f"{max(r['energy'] for v in by.values() for r in v):.3f}", flush=True)
        for mode in ("high", "low", "bimodal", "closed"):
            rows = pick(by, n_rows, mode)
            if rows is None:
                print(f"  {mode:<8} short pool; skipped", flush=True)
                continue
            es = [r["energy"] for r in rows]
            doc = dict(contract)
            doc["seed"] = ",".join(str(s) for s in SEEDS)
            json.dump(doc, open(settings.CONTRACT_PATH, "w"))
            json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))
            json.dump([{"experiment_id": f"exp-{i:05d}", "guideRNA": r["guide"],
                        "target_alignment_start": r["start"],
                        "target_alignment_end": r["start"] + r["length"],
                        "strand": r["strand"], "mutation": r["mutation"],
                        "cas_system": r["cas_system"],
                        "cell_type": contract.get("cell_type")}
                       for i, r in enumerate(rows)],
                      open(settings.MINER_SUBMISSION_PATH, "w"))
            run_stage12(cell_types)
            per = []
            for s in SEEDS:
                run_stage3(seed=s)
                run_stage4(seed=s)
                fr = json.load(open(settings.FINAL_REWARD_PATH))
                res = run_stage5()
                tgt = fr.get("model_results") or {}
                per.append({"cons": res["consistency_factor"],
                            "r2": {k: round(v.get("r2_mean", 0), 4) for k, v in tgt.items()},
                            "mae": {k: round(v.get("mae_mean", 0), 4) for k, v in tgt.items()}})
            avg = {k: st.mean(p["r2"].get(k, 0) for p in per)
                   for k in (per[0]["r2"] if per else {})}
            rec = {"task": task["id"], "cell_type": contract["cell_type"], "mode": mode,
                   "energy_mean": st.mean(es), "energy_sd": st.pstdev(es),
                   "cons": st.mean(p["cons"] for p in per), "r2_per_target": avg,
                   "detail": per}
            out.append(rec)
            json.dump(out, open(OUT, "w"), indent=1)
            print(f"  {mode:<8} energy {st.mean(es):.3f}+-{st.pstdev(es):.3f}  "
                  f"cons {st.mean(p['cons'] for p in per):.4f}  "
                  f"R2 " + " ".join(f"{k}={v:+.3f}" for k, v in avg.items()), flush=True)
        print(flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
