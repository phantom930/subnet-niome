#!/usr/bin/env python3
"""partial_probe.py — is consistency a smooth function of the *fraction* of complying rows?

Everything else about the 0.45-0.65 regime has been falsified: the leaders have no better floor
(their medians equal ours, their p10s are worse), fidelity does not cause it (0.999 fidelity scored
0.088-0.102 on 64 random seeds), and feature spread does not create learnable outcomes. What
survives is arithmetic: the observed values read as k=2 with per-seed compliance ~0.9 rather than
1.0, and reaching k>=2 on ~12% of rounds needs a band near 180 seeds -- an order of magnitude past
the ~15 our pool bound allows for *full* compliance.

Partial compliance is the only remaining way to get both. Pinning 200 of 250 rows over a wide
window is a far weaker constraint than pinning all 250 over a narrow one, because the binomial
requirement scales with the number of rows that must agree.

This measures the payoff curve before attempting the search. For one contract and one seed, rows
are drawn so that exactly ``f`` of them satisfy the ``hdr`` rule at that seed and the rest do not,
and the submission is scored at that seed. That isolates the question -- no window, no min-union,
no search -- and asks what consistency a given compliant fraction is worth.

The shape decides the next step:

* smooth, passing through 0.45-0.65 near f ~ 0.85  -> the regime is partial compliance, and the
  target is exact; the search problem becomes "f of the rows over a wide window".
* a step from ~0.10 to 1.0 with nothing between    -> partial compliance is worthless and the
  regime is something not yet imagined.

stage 4 fits a forest per target under KFold, so an intermediate curve requires the complying rows
to be *identifiable from the features* -- the non-compliant rows are ordinary candidates here, not
synthetic noise, so the answer reflects real guide populations.
"""
import os

os.environ["NIOME_INSTANCE"] = "win20"

import json
import statistics as st
import sys
import time
from collections import defaultdict

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
import logging
logging.basicConfig(level=logging.ERROR)

import genExp as G
from niome_subnet.genomics import seed_agnostic as SA
from niome_subnet.genomics.validation import run_stage12, run_stage3, run_stage4, run_stage5
from niome_subnet.utils import settings

TASKS = ["task-k562.json", "task-cd34.json"]
PROBE_SEEDS = [271, 604, 833]
# Row *counts* that break the rule, not fractions: the first sweep showed consistency climbing
# only to 0.21-0.40 by f=0.95 and then jumping to 1.000 at f=1.00, so everything interesting lives
# in the last 12 rows. The observed 0.45-0.65 regime sits inside that jump.
BAD_ROWS = [0, 1, 2, 3, 5, 8, 12, 25]
MAX_DISTANCE = 600
GC_BAND = (0.30, 0.70)
VARIANTS = 1200
PER_BUCKET = 400
OUT = "partial_probe_fine.json"


def pools(contract, reference, cell_types, seed):
    """Candidates for one seed, split per stage-5 cell into rule-satisfying and rule-breaking."""
    pinned = dict(contract)
    pinned["seed"] = seed
    ctx = G.build_context(pinned, reference, cell_types)
    sites = G.enumerate_sites(ctx, 3000, (20, 23))
    jobs = [(s, m) for s in sites for m in ctx.mutations
            if abs(s.start - ctx.mutation_map[m]) <= MAX_DISTANCE]
    jobs.sort(key=lambda j: abs(j[0].start - ctx.mutation_map[j[1]]))
    yes, no = defaultdict(list), defaultdict(list)
    for site, mutation in jobs:
        for guide in SA.enumerate_variants(site, ctx, GC_BAND[0], GC_BAND[1],
                                           ctx.max_mismatches, True, VARIANTS):
            entry = G.build_valid_entry(G.make_experiment(site, guide, mutation, ctx, "c"), ctx)
            if entry is None:
                continue
            rec = G.simulate(entry, ctx)
            key = (mutation, site.cas, site.strand)
            bucket = yes if rec["outcome"] == "HDR" else no
            if len(bucket[key]) >= PER_BUCKET:
                continue
            bucket[key].append({"guide": guide, "mutation": mutation, "cas_system": site.cas,
                                "strand": site.strand, "start": site.start, "length": site.length,
                                "weighted_score": entry["stage2"]["weighted_score"]})
    for d in (yes, no):
        for v in d.values():
            v.sort(key=lambda r: -r["weighted_score"])
    return ctx, yes, no


def compose(yes, no, n_rows, frac):
    """n_rows split so ``frac`` satisfy the rule, both halves spread evenly over the 8 cells."""
    cells = sorted(set(yes) | set(no))
    if len(cells) < 8:
        return None
    n_yes = int(round(n_rows * frac))
    out, take = [], {}
    for i, key in enumerate(cells):
        base, extra = divmod(n_yes, len(cells))
        take[key] = base + (1 if i < extra else 0)
    for key in cells:
        out += yes.get(key, [])[:take[key]]
    need = n_rows - len(out)
    per = defaultdict(int)
    for i, key in enumerate(cells):
        base, extra = divmod(need, len(cells))
        per[key] = base + (1 if i < extra else 0)
    for key in cells:
        out += no.get(key, [])[:per[key]]
    if len(out) < n_rows:
        for key in cells:                       # top up from whichever pool still has rows
            for src in (no, yes):
                while len(out) < n_rows and len(src.get(key, [])) > per[key] + take[key]:
                    out.append(src[key][per[key] + take[key]])
                    per[key] += 1
    return out[:n_rows]


def score_at(rows, contract, reference, cell_types, seed):
    doc = dict(contract)
    doc["seed"] = str(seed)
    json.dump(doc, open(settings.CONTRACT_PATH, "w"))
    json.dump(reference, open(settings.HBB_REFERENCE_PATH, "w"))
    json.dump([{"experiment_id": f"exp-{i:05d}", "guideRNA": r["guide"],
                "target_alignment_start": r["start"],
                "target_alignment_end": r["start"] + r["length"], "strand": r["strand"],
                "mutation": r["mutation"], "cas_system": r["cas_system"],
                "cell_type": contract.get("cell_type")} for i, r in enumerate(rows)],
              open(settings.MINER_SUBMISSION_PATH, "w"))
    run_stage12(cell_types)
    run_stage3(seed=seed)
    run_stage4(seed=seed)
    return run_stage5()


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    out = []
    for tf in TASKS:
        task = json.load(open(tf))
        contract = dict(task["content"]["contract"])
        reference = task["content"]["hbb_reference"]
        n_rows = contract["rules"].get("max_experiments") or 250
        print(f"=== {tf}  {contract['cell_type']}  {task['id'][:8]} ===", flush=True)
        for seed in PROBE_SEEDS:
            t0 = time.monotonic()
            _ctx, yes, no = pools(contract, reference, cell_types, seed)
            ny, nn = sum(len(v) for v in yes.values()), sum(len(v) for v in no.values())
            print(f"  seed {seed}: {ny} rule-satisfying / {nn} rule-breaking candidates "
                  f"over {len(set(yes) | set(no))} cells ({time.monotonic() - t0:.0f}s)", flush=True)
            print(f"    {'bad':>5} {'actual':>7} {'weighted':>9} {'fid':>6} {'cons':>7} "
                  f"{'final':>8}", flush=True)
            for bad in BAD_ROWS:
                frac = (n_rows - bad) / n_rows
                rows = compose(yes, no, n_rows, frac)
                if rows is None or len(rows) < n_rows:
                    print(f"    bad={bad:<3} short pool; skipped", flush=True)
                    continue
                res = score_at(rows, contract, reference, cell_types, seed)
                # how many of the chosen rows actually satisfy the rule
                ykeys = {(r["guide"], r["start"], r["strand"], r["mutation"])
                         for v in yes.values() for r in v}
                actual = sum(1 for r in rows
                             if (r["guide"], r["start"], r["strand"],
                                 r["mutation"]) in ykeys) / len(rows)
                rec = {"task": task["id"], "cell_type": contract["cell_type"], "seed": seed,
                       "bad_rows": bad, "target_frac": frac, "actual_frac": actual,
                       "rows": len(rows),
                       "weighted": res["total_weighted_score"],
                       "fidelity": res["distribution_fidelity_factor"],
                       "consistency": res["consistency_factor"], "final": res["final_score"]}
                out.append(rec)
                json.dump(out, open(OUT, "w"), indent=1)
                print(f"    {bad:>5} {actual:>7.3f} {res['total_weighted_score']:>9.1f} "
                      f"{res['distribution_fidelity_factor']:>6.3f} "
                      f"{res['consistency_factor']:>7.4f} {res['final_score']:>8.2f}", flush=True)
        print(flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
