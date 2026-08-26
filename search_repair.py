#!/usr/bin/env python3
"""search_repair.py — Cas9 candidates at the cut_p clamp whose repair mode tracks the mh coin.

Two conditions, and they are of different kinds:

* **cut_p >= 0.99** is a *design* property, no seed involved. ``cut_probability`` for Cas9 is
  ``min(0.99, max(0.4, 0.86 + 0.18*energy))``, so the floor is reachable only at the clamp, and
  only for ``energy >= (0.99 - 0.86) / 0.18 = 0.7222``. Since
  ``energy = accessibility * (1.8*gc + 0.6*exp(-d/1500) + region_offset)``, whether a cell type can
  reach it at sane GC is decided by accessibility alone — see ``--report-gate``.
* **mh -> HDR, else BLUNT_NHEJ** is a *seed* property: the repair mode is the third draw of
  ``random.Random(experiment_seed(seed, design))``, after the microhomology coin and the cut coin.
  So it holds for a given guide under some seeds and not others, and this tool tests one seed --
  the contract's, by default.

This is ``genExp``'s "mh" construction minus its ``indel_length == 1`` pin (registered there as
``mh_any``), which is why the yield is high: no 0.70 pin to hit, so roughly

    P(cut) * [ P(mh)*P(HDR|mh) + P(!mh)*P(BLUNT|!mh) ] ~= 0.99 * 0.41 ~= 40%

of enumerated variants qualify, against ~29% for "mh". The cost is stage 4's ``indel_length``
target, which stops being learnable.

    python search_repair.py                          # task.json, the contract's own seed
    python search_repair.py --seed 0                 # what a broadcast contract looks like
    python search_repair.py --gc-min 0.30 --gc-max 0.70 --variants 40000
    python search_repair.py --seed-window 100 999    # per-candidate hit rate across a seed range

``--seed-window`` exists to answer the obvious follow-up: this construction cannot be made
seed-agnostic. At a ~40% per-seed hit rate, surviving n seeds costs 0.40**n, so the best guide over
a 900-seed window holds for a few hundred of them and none hold for all. The tool measures it
rather than asserting it.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import logging
import math
import os
import random
import sys
import time
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import niome_subnet.utils.settings as settings  # noqa: E402  (must precede bittensor imports)
import genExp as G  # noqa: E402

from niome_subnet.genomics import seed_agnostic as SA  # noqa: E402
from niome_subnet.genomics.validation import stage3  # noqa: E402

logger = logging.getLogger("search_repair")

STATE: dict = {}


def energy_of(gc: float, distance: int, accessibility: float, region_offset: float) -> float:
    """stage3.sequence_energy for an explicit gc, which is the one input that varies per variant."""
    return max(0.0, min(1.0, accessibility * (
        1.8 * gc + 0.6 * math.exp(-distance / 1500) + region_offset)))


def gc_floor_for_clamp(accessibility: float, distance: int, region_offset: float,
                       cut_p_floor: float) -> float:
    """The GC a Cas9 guide needs before cut_p reaches ``cut_p_floor``. >1 means unreachable.

    Inverts cut_probability then sequence_energy. This is why HEK293 is hopeless here: at
    accessibility 0.35 the clamp needs GC ~0.81, far outside any band that also scores well on
    stage 2's gc_score (which peaks at 0.50).
    """
    energy_needed = (cut_p_floor - 0.86) / 0.18
    if energy_needed <= 0:
        return 0.0
    inner = energy_needed / accessibility - 0.6 * math.exp(-distance / 1500) - region_offset
    return inner / 1.8


def draw(guide: str, gc: float, energy: float, cut_p: float, mutation: str, cas: str,
         start: int, strand: str, seed: int) -> tuple[bool, str]:
    """Replay stage 3's draw sequence for one (design, seed). Returns (mh, outcome).

    The seed derivation is replicated (it is a hash, and MT.verify holds it bit-exact against the
    real pipeline), but every *probability* comes from the validator's own stage-3 functions, so a
    formula change there re-prices this search automatically rather than silently diverging.
    """
    digest = hashlib.sha256(
        (str(seed) + f"|{mutation}|{cas}|{guide}|{start}|{strand}").encode()).digest()
    rng = random.Random(int.from_bytes(digest[-4:], "big"))
    mh = stage3.microhomology_trigger({"gc": gc}, rng)      # draw 1
    if rng.random() > cut_p:                                # draw 2
        return mh, "no_cut"
    return mh, stage3.repair_mode(cas, energy, mh, rng)     # draw 3


# Repair-mode constructions, keyed the way genExp's CONSTRUCTIONS are. Both make is_hdr an exact
# function of design-derived features; they differ in how hard they are per seed. For Cas12a at the
# energy clamp: "hdr" holds 51.2% of the time (HDR is the single most likely mode, weight 0.59 vs
# blunt 0.35), "mh_any" only 41.0%, because it has to track the mh coin into whichever branch.
RULES = {
    "mh_any": lambda mh, outcome: outcome == ("HDR" if mh else "BLUNT_NHEJ"),
    "hdr": lambda mh, outcome: outcome == "HDR",
}


def rule_holds(mh: bool, outcome: str, rule: str = "mh_any") -> bool:
    return RULES[rule](mh, outcome)


def window_fails(guide: str, gc: float, energy: float, cut_p: float, mutation: str, cas: str,
                 start: int, strand: str, lo: int, hi: int,
                 max_cut_fail: int | None = None,
                 max_rule_fail: int | None = None,
                 rule: str = "mh_any") -> tuple[int, int] | None:
    """Count (cut failures, rule failures) for one guide across seeds [lo, hi].

    A ``no_cut`` seed counts against *both*: the row did not cut, and a row that did not cut
    satisfies neither branch of ``mh_any``. So rule_fails >= cut_fails always, and the rule is the
    binding condition. ``None`` is an early-out once either allowance is blown -- which is what
    makes a bank scan over millions of guides cheap, since a guide that fails ~60% of seeds dies
    within the first few dozen.
    """
    cut_fails = rule_fails = 0
    for sd in range(lo, hi + 1):
        mh, outcome = draw(guide, gc, energy, cut_p, mutation, cas, start, strand, sd)
        if outcome == "no_cut":
            cut_fails += 1
            rule_fails += 1
        elif not rule_holds(mh, outcome, rule):
            rule_fails += 1
        if max_cut_fail is not None and cut_fails > max_cut_fail:
            return None
        if max_rule_fail is not None and rule_fails > max_rule_fail:
            return None
    return cut_fails, rule_fails


def _init(payload: dict) -> None:
    logging.getLogger().setLevel(logging.ERROR)
    ctx = G.build_context(payload["contract"], payload["reference"], payload["cell_types"])
    contract = payload["contract"]
    regions = contract.get("mutation_regions") or {}
    STATE.update(
        ctx=ctx,
        args=payload["args"],
        accessibility=payload["cell_types"].get(contract.get("cell_type"), {})
                                           .get("accessibility", 1.0),
        region_offset={m: stage3.REGION_ENERGY_OFFSETS.get(regions.get(m), 0.0)
                       for m in ctx.mutations},
        sites=G.enumerate_sites(ctx, payload["args"].flank, (20, 23)),
        profile_per_target=payload.get("profile_per_target", 0),
    )


def scan_target(job: tuple[int, str]) -> dict:
    """Enumerate one (site, mutation) target and keep every variant meeting both conditions."""
    site_index, mutation = job
    ctx, args = STATE["ctx"], STATE["args"]
    site = STATE["sites"][site_index]
    distance = abs(site.start - ctx.mutation_map[mutation])
    if args.max_distance and distance > args.max_distance:
        return {"kept": [], "examined": 0, "gated": 0}

    accessibility = STATE["accessibility"]
    region_offset = STATE["region_offset"][mutation]
    guides = SA.enumerate_variants(site, ctx, args.gc_min, args.gc_max, ctx.max_mismatches,
                                   True, args.variants)
    if not guides:
        return {"kept": [], "examined": 0, "gated": 0}

    # cut_p depends on GC only through energy, so one lookup per distinct GC count serves the
    # whole target rather than recomputing per guide.
    by_gc: dict[int, tuple[float, float, float]] = {}

    def params(guide: str):
        n = sum(b in "GC" for b in guide)
        got = by_gc.get(n)
        if got is None:
            gc = n / site.length
            energy = energy_of(gc, distance, accessibility, region_offset)
            got = by_gc[n] = (gc, energy, stage3.cut_probability(site.cas, energy))
        return got

    kept, examined, gated = [], 0, 0
    for guide in guides:
        examined += 1
        gc, energy, cut_p = params(guide)
        if cut_p < args.cut_p_floor - 1e-12:
            gated += 1
            continue
        mh, outcome = draw(guide, gc, energy, cut_p, mutation, site.cas,
                           site.start, site.strand, args.seed)
        if not rule_holds(mh, outcome, args.rule):
            continue
        kept.append({
            "guide": guide, "mutation": mutation, "cas_system": site.cas, "strand": site.strand,
            "target_alignment_start": site.start, "length": site.length,
            "gc": round(gc, 4), "distance": distance, "energy": round(energy, 4),
            "cut_p": round(cut_p, 4), "mh": mh, "outcome": outcome,
        })
    return {"kept": kept, "examined": examined, "gated": gated}


def scan_target_window(job: tuple[int, str]) -> dict:
    """Window mode: keep guides whose cut and rule failures both stay inside their allowances.

    ``--profile`` switches off the early-out for a sample so the *achievable* frontier can be
    measured -- the minimum rule-failure count anywhere in the pool -- rather than only counting
    how many clear a threshold that may well be unreachable.
    """
    site_index, mutation = job
    ctx, args = STATE["ctx"], STATE["args"]
    site = STATE["sites"][site_index]
    distance = abs(site.start - ctx.mutation_map[mutation])
    if args.max_distance and distance > args.max_distance:
        return {"kept": [], "examined": 0, "gated": 0, "profile": []}

    accessibility = STATE["accessibility"]
    region_offset = STATE["region_offset"][mutation]
    guides = SA.enumerate_variants(site, ctx, args.gc_min, args.gc_max, ctx.max_mismatches,
                                   True, args.variants)
    if not guides:
        return {"kept": [], "examined": 0, "gated": 0, "profile": []}

    lo, hi = args.window
    by_gc: dict[int, tuple[float, float, float]] = {}

    def params(guide: str):
        n = sum(b in "GC" for b in guide)
        got = by_gc.get(n)
        if got is None:
            gc = n / site.length
            energy = energy_of(gc, distance, accessibility, region_offset)
            got = by_gc[n] = (gc, energy, stage3.cut_probability(site.cas, energy))
        return got

    kept, examined, gated, profile = [], 0, 0, []
    profile_quota = STATE["profile_per_target"]
    for guide in guides:
        examined += 1
        gc, energy, cut_p = params(guide)
        if cut_p < args.cut_p_floor - 1e-12:
            gated += 1
            continue
        if profile_quota and len(profile) < profile_quota:
            full = window_fails(guide, gc, energy, cut_p, mutation, site.cas,
                                site.start, site.strand, lo, hi, rule=args.rule)
            profile.append(full)                       # no early-out: the true frontier
            counts = full
        else:
            counts = window_fails(guide, gc, energy, cut_p, mutation, site.cas,
                                  site.start, site.strand, lo, hi,
                                  args.max_fail, args.rule_max_fail, rule=args.rule)
        if counts is None:
            continue
        cut_fails, rule_fails = counts
        if cut_fails > args.max_fail or rule_fails > args.rule_max_fail:
            continue
        kept.append({
            "guide": guide, "mutation": mutation, "cas_system": site.cas, "strand": site.strand,
            "target_alignment_start": site.start, "length": site.length,
            "gc": round(gc, 4), "distance": distance, "energy": round(energy, 4),
            "cut_p": round(cut_p, 4), "cut_fails": cut_fails, "rule_fails": rule_fails,
        })
    return {"kept": kept, "examined": examined, "gated": gated, "profile": profile}


def seed_window_profile(records: list[dict], lo: int, hi: int, sample: int) -> dict:
    """How many seeds in [lo, hi] each of a sample of candidates satisfies the rule under.

    The point is the maximum: if no candidate clears the whole window, no all-seed construction on
    this rule exists, and the count says how far off it is.
    """
    seeds = range(lo, hi + 1)
    span = len(list(seeds))
    hits = []
    for rec in records[:sample]:
        n = sum(1 for sd in range(lo, hi + 1)
                if rule_holds(*draw(rec["guide"], rec["gc"], rec["energy"], rec["cut_p"],
                                    rec["mutation"], rec["cas_system"],
                                    rec["target_alignment_start"], rec["strand"], sd)))
        hits.append(n)
    return {"window": [lo, hi], "span": span, "sampled": len(hits),
            "max": max(hits) if hits else 0, "mean": (sum(hits) / len(hits)) if hits else 0.0,
            "all_clear": sum(1 for n in hits if n == span)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--task", default="test/task.json")
    ap.add_argument("--cell-types", default="test/cell_types.json")
    ap.add_argument("--out", default="test/repair_candidates.json")
    ap.add_argument("--seed", type=int, default=None,
                    help="round seed to test the rule at (default: the contract's own seed)")
    ap.add_argument("--cas", default="Cas9", choices=("Cas9", "Cas12a"))
    ap.add_argument("--cut-p-floor", type=float, default=0.99)
    ap.add_argument("--gc-min", type=float, default=0.40)
    ap.add_argument("--gc-max", type=float, default=0.60)
    ap.add_argument("--variants", type=int, default=20000, help="cap on variants per target")
    ap.add_argument("--flank", type=int, default=3000)
    ap.add_argument("--max-distance", type=int, default=0,
                    help="0 = every enumerated target, no distance cap")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    ap.add_argument("--seed-window", type=int, nargs=2, metavar=("LO", "HI"), default=None,
                    help="also profile how many seeds in [LO, HI] each candidate holds under")
    ap.add_argument("--window-sample", type=int, default=2000)
    ap.add_argument("--window", type=int, nargs=2, metavar=("LO", "HI"), default=None,
                    help="window mode: require BOTH conditions across seeds LO..HI, subject to "
                         "--max-fail and --rule-max-fail, instead of testing one seed")
    ap.add_argument("--max-fail", type=int, default=22,
                    help="window mode: seeds the guide may fail to cut under")
    ap.add_argument("--rule-max-fail", type=int, default=None,
                    help="window mode: seeds the guide may break mh_any under "
                         "(default: same as --max-fail). A no_cut seed counts against both.")
    ap.add_argument("--rule", default="mh_any", choices=tuple(RULES),
                    help="repair-mode construction the candidate must satisfy")
    ap.add_argument("--cell-type", default=None,
                    help="override the contract's cell type (changes accessibility, hence energy)")
    ap.add_argument("--profile", type=int, default=0,
                    help="window mode: fully evaluate this many guides with no early-out, to "
                         "measure the achievable frontier rather than only counting threshold hits")
    ap.add_argument("--max-records", type=int, default=250000,
                    help="cap the records written; kept slice is ranked by what stage 2 rewards "
                         "(distance asc, then |gc - 0.50| asc). 0 = write everything")
    ap.add_argument("--score-top", type=int, default=0,
                    help="compute the exact stage-2 weighted_score for this many of the kept "
                         "records (needs the real entry per guide, so it is the slow part)")
    ap.add_argument("--verify", type=int, default=200,
                    help="cross-check this many candidates against stage3.simulate")
    ap.add_argument("--report-gate", action="store_true",
                    help="print the GC each cell type needs to reach the cut_p floor, then exit")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(message)s")

    task = json.load(open(args.task))
    contract = task["content"]["contract"]
    reference = task["content"]["hbb_reference"]
    cell_types = json.load(open(args.cell_types))
    if args.seed is None:
        args.seed = contract.get("seed", 0)
    if args.cell_type:
        contract = dict(contract, cell_type=args.cell_type)
    if args.rule_max_fail is None:
        args.rule_max_fail = args.max_fail
    if args.window and args.cut_p_floor == 0.99 and args.cas == "Cas12a":
        # Cas12a's cut_p caps at 0.96, so the Cas9 default floor would gate the whole pool out.
        # In window mode the per-seed cut test replaces the static floor anyway.
        args.cut_p_floor = 0.0

    if args.report_gate:
        print(f"\ncut_p floor {args.cut_p_floor} needs "
              f"energy >= {(args.cut_p_floor - 0.86) / 0.18:.4f} (Cas9)\n")
        print(f"  {'cell type':<14}{'access':>8}{'GC needed':>12}   reachable in a 40-60% band?")
        for name, info in sorted(cell_types.items()):
            need = gc_floor_for_clamp(info["accessibility"], 0, 0.0, args.cut_p_floor)
            verdict = "yes" if need <= 0.60 else "NO"
            print(f"  {name:<14}{info['accessibility']:>8.2f}{max(need, 0.0):>12.3f}   {verdict}")
        return

    ctx = G.build_context(contract, reference, cell_types)
    sites = G.enumerate_sites(ctx, args.flank, (20, 23))
    jobs = [(i, m) for i, s in enumerate(sites) if s.cas == args.cas for m in ctx.mutations
            if not args.max_distance
            or abs(s.start - ctx.mutation_map[m]) <= args.max_distance]

    logger.info("task %s | cell_type=%s seed=%s | %d %s targets from %d sites",
                task.get("id", args.task), contract.get("cell_type"), args.seed,
                len(jobs), args.cas, len(sites))
    rule_text = {"hdr": "outcome == HDR",
                 "mh_any": "mh -> HDR, else BLUNT_NHEJ"}[args.rule]
    logger.info("conditions: cut_p >= %.3f and (%s) | GC %.2f-%.2f, <= %d variants/target",
                args.cut_p_floor, rule_text, args.gc_min, args.gc_max, args.variants)

    payload = {"contract": contract, "reference": reference, "cell_types": cell_types,
               "args": args,
               "profile_per_target": (max(1, -(-args.profile // max(1, len(jobs))))
                                     if (args.window and args.profile) else 0)}
    worker = scan_target_window if args.window else scan_target
    t0 = time.time()
    records, examined, gated, profile = [], 0, 0, []
    with Pool(args.jobs, initializer=_init, initargs=(payload,)) as pool:
        for i, res in enumerate(pool.imap_unordered(worker, jobs, chunksize=4), 1):
            records.extend(res["kept"])
            examined += res["examined"]
            gated += res["gated"]
            profile.extend(res.get("profile", ()))
            if i % 200 == 0:
                logger.info("  %d/%d targets | %d candidates from %d variants",
                            i, len(jobs), len(records), examined)
    elapsed = time.time() - t0

    by_cell = Counter((r["mutation"], r["strand"]) for r in records)
    by_branch = Counter("mh -> HDR" if r["mh"] else "!mh -> BLUNT_NHEJ"
                        for r in records if "mh" in r)
    per_target = defaultdict(int)
    for rec in records:
        per_target[(rec["target_alignment_start"], rec["mutation"])] += 1

    scope = (f"seeds {args.window[0]}-{args.window[1]}" if args.window else f"seed {args.seed}")
    print(f"\n{'=' * 78}")
    print(f"{args.cas} candidates satisfying both conditions | "
          f"{contract.get('cell_type')} (acc {cell_types.get(contract.get('cell_type'), {}).get('accessibility')}) | {scope}")
    print(f"{'=' * 78}")
    print(f"  variants examined        {examined:,}")
    print(f"  rejected on cut_p        {gated:,}"
          f"{'  (none — the whole GC band clears the clamp)' if not gated else ''}")
    survivors = examined - gated
    print(f"  passed the cut_p gate    {survivors:,}")
    print(f"  CANDIDATES               {len(records):,}"
          f"   ({100 * len(records) / survivors:.1f}% of those that passed)" if survivors else "")
    print(f"  distinct guides          {len({r['guide'] for r in records}):,}")
    print(f"  targets with >=1         {len(per_target):,} of {len(jobs):,}")
    print(f"  scan time                {elapsed:.1f}s on {args.jobs} cores")
    if args.window:
        lo, hi = args.window
        span = hi - lo + 1
        print(f"\n  window [{lo}, {hi}] ({span} seeds), allowances: "
              f"cut <= {args.max_fail} fails, {args.rule} <= {args.rule_max_fail} fails")
        if profile:
            cut_f = sorted(pf[0] for pf in profile)
            rule_f = sorted(pf[1] for pf in profile)
            def pct(xs, q):
                return xs[min(len(xs) - 1, int(q * len(xs)))]
            print(f"\n  measured frontier over {len(profile):,} fully-evaluated guides "
                  f"(no early-out)")
            print(f"    {'':<14}{'best':>8}{'p1':>8}{'median':>8}{'worst':>8}")
            print(f"    {'cut fails':<14}{cut_f[0]:>8}{pct(cut_f, 0.01):>8}"
                  f"{pct(cut_f, 0.5):>8}{cut_f[-1]:>8}")
            print(f"    {args.rule + ' fails':<14}{rule_f[0]:>8}{pct(rule_f, 0.01):>8}"
                  f"{pct(rule_f, 0.5):>8}{rule_f[-1]:>8}")
            print(f"\n    best {args.rule}: {rule_f[0]} fails => holds under "
                  f"{span - rule_f[0]}/{span} seeds")
            print(f"\n  candidates by allowance (from the profiled sample)")
            print(f"    {'rule max-fail':>14}{'candidates':>14}")
            for thresh in (0, 22, 50, 100, 200, 300, 400, 450, 500, 550):
                n = sum(1 for pf in profile if pf[1] <= thresh and pf[0] <= max(thresh, args.max_fail))
                print(f"    {thresh:>14}{n:>14,}")
            need = rule_f[0]
            print(f"\n    smallest allowance with any candidate: {need}"
                  f"  (you asked for {args.rule_max_fail})")
        print(f"\n  CANDIDATES at the requested allowances: {len(records):,}")
        if records:
            best = min(records, key=lambda r: r["rule_fails"])
            print(f"    best: cut_fails={best['cut_fails']} rule_fails={best['rule_fails']} "
                  f"gc={best['gc']} d={best['distance']}")

    if not args.window:
        print("\n  by branch of the rule")
        for branch, n in by_branch.most_common():
            print(f"    {branch:<22}{n:>10,}")
    # Distance matters twice over: it feeds energy (so the cut_p gate tightens with it) and it is
    # stage 2's dist_score, which is what a distant candidate actually costs. Raw count alone would
    # read as more usable supply than there is.
    print("\n  by distance to the mutation")
    bands = [(0, 200), (201, 500), (501, 1000), (1001, 2000), (2001, 10 ** 9)]
    for lo, hi in bands:
        n = sum(1 for r in records if lo <= r["distance"] <= hi)
        label = f"{lo}-{hi}" if hi < 10 ** 9 else f"{lo}+"
        bar = "#" * int(40 * n / max(1, len(records)))
        print(f"    {label:>10} bp {n:>10,}  {bar}")

    print("\n  by cell (mutation x strand)")
    for (mut, strand), n in sorted(by_cell.items()):
        print(f"    {mut} {strand:<4}{n:>12,}")

    if args.verify and records and not args.window:
        step = max(1, len(records) // args.verify)
        checked = mismatched = unresolved = 0
        for rec in records[::step][:args.verify]:
            # Length is part of the identity: a start can carry both a 20- and a 23-mer site, and
            # rebuilding a 20-mer guide on the 23-mer site fails stage 1 for reasons that have
            # nothing to do with the draw. Kept separate from a real mismatch on purpose.
            site = next((s for s in sites
                         if s.start == rec["target_alignment_start"] and s.strand == rec["strand"]
                         and s.cas == rec["cas_system"] and s.length == rec["length"]), None)
            entry = None if site is None else G.build_valid_entry(
                G.make_experiment(site, rec["guide"], rec["mutation"], ctx, "verify"), ctx)
            if entry is None:
                unresolved += 1
                continue
            result = stage3.simulate(entry, args.seed)
            checked += 1
            if not (result["mh"] == rec["mh"] and result["outcome"] == rec["outcome"]
                    and G.CONSTRUCTIONS["mh_any"](result, entry)):
                mismatched += 1
        print(f"\n  verified against stage3.simulate: {checked} simulated, "
              f"{mismatched} mismatch(es), {unresolved} unresolved")

    if args.seed_window:
        lo, hi = args.seed_window
        prof = seed_window_profile(records, lo, hi, args.window_sample)
        print(f"\n  seeds in [{lo}, {hi}] each candidate holds under "
              f"({prof['sampled']:,} sampled of {prof['span']})")
        print(f"    best      {prof['max']}/{prof['span']}")
        print(f"    mean      {prof['mean']:.1f}/{prof['span']}")
        print(f"    all-clear {prof['all_clear']}  "
              f"-> {'a seed-agnostic build exists' if prof['all_clear'] else 'no seed-agnostic build on this rule'}")

    # Rank before capping. Both keys are what stage 2 actually pays for -- dist_score falls with
    # distance and gc_score peaks at 0.50 -- so the kept slice is the usable end of the supply
    # rather than an arbitrary 250k of it. The full count is reported either way.
    records.sort(key=lambda r: (r["distance"], abs(r["gc"] - 0.50)))
    kept = records if not args.max_records else records[:args.max_records]

    if args.score_top:
        by_key = {(s.start, s.strand, s.cas, s.length): s for s in sites}
        # Round-robin the cells. A straight top-N slice of a (distance, gc)-sorted list is all one
        # target -- the nearest one -- so its scores would describe that target, not the bank.
        per_cell: dict[tuple, list] = {}
        for rec in kept:
            per_cell.setdefault((rec["mutation"], rec["strand"]), []).append(rec)
        sample = [rec for group in itertools.zip_longest(*per_cell.values())
                  for rec in group if rec is not None][:args.score_top]
        scored = 0
        for rec in sample:
            site = by_key.get((rec["target_alignment_start"], rec["strand"],
                               rec["cas_system"], rec["length"]))
            entry = None if site is None else G.build_valid_entry(
                G.make_experiment(site, rec["guide"], rec["mutation"], ctx, "score"), ctx)
            if entry is not None:
                rec["weighted_score"] = entry["stage2"]["weighted_score"]
                scored += 1
        top = sorted((r for r in kept if "weighted_score" in r),
                     key=lambda r: -r["weighted_score"])
        print(f"\n  exact stage-2 score on the top {scored:,}: "
              f"best {top[0]['weighted_score']:.4f}, median "
              f"{top[len(top) // 2]['weighted_score']:.4f}" if top else "")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump({"seed": args.seed, "cas": args.cas, "cell_type": contract.get("cell_type"),
                   "cut_p_floor": args.cut_p_floor, "gc_band": [args.gc_min, args.gc_max],
                   "variants_per_target": args.variants, "examined": examined,
                   "gated_on_cut_p": gated, "candidates": len(records),
                   "records_written": len(kept),
                   "by_cell": {f"{m}|{s}": n for (m, s), n in by_cell.items()},
                   "by_distance_band": {f"{lo}-{hi}": sum(1 for r in records
                                                          if lo <= r["distance"] <= hi)
                                        for lo, hi in bands},
                   "records": kept}, fh)
    print(f"\n  {len(records):,} candidates found; wrote the best {len(kept):,} to {args.out}")


if __name__ == "__main__":
    main()
