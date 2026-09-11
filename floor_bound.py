#!/usr/bin/env python3
"""floor_bound.py -- the analytic ceiling on `consistency_factor` from row DISTRIBUTION alone.

The question this settles, before any GPU build. The field lifted its off-band per-seed value to
0.20-0.25 while still reaching k=2 (CLAUDE.md, "Nobody has a better floor"). Our 0.101 is attributed
there to `avg_nmae ~ 0.63` and a forest that overfits to negative R2. Two mechanisms could produce
their floor and they have opposite consequences:

  PINNING     every row cuts on a seed -> `is_cut` constant -> that seed is worth 0.237. Needs a
              wide cut-clean set, which competes with the HDR band for the same guides: measured
              dead twice (the conjunction caps the clean set at 52-80 of 900; the hybrid split
              scores below the pure floor).
  DISTRIBUTION shape the rows so the three targets are genuinely PREDICTABLE from
              [gc, distance, gc_score, dist_score, consistency, energy, mh] -- positive R2 and low
              nmae with nothing pinned. This does NOT compete with the band, so if it reaches
              0.20-0.25 it is the whole answer.

This computes the DISTRIBUTION ceiling exactly, with no simulation and no build, by giving the
regressor the best it could possibly do: predicting each row's true conditional mean. For a target
y with per-row success probability p_i, an oracle regressor gives

    r2   = Var(p) / (E[p(1-p)] + Var(p))            (explained / total variance)
    nmae = mean(2 p_i (1-p_i)) / sqrt(E[p(1-p)] + Var(p))

so heterogeneity in p across rows is the only thing that buys either term -- a homogeneous
submission has Var(p) = 0, r2 = 0 and nmae = 2 sqrt(p(1-p)), which is where our 0.101 comes from.
`indel_length` is handled by its exact conditional mean and variance over the three repair modes.

Every probability is read straight off stage3: `sequence_energy`, `cut_probability`,
`microhomology_trigger`, `repair_mode`, `sample_indel_length`. The real forest OVERFITS (measured
avg_r2 -0.249 off band), so this is an upper bound and a generous one: if it cannot reach 0.24, the
distributional route is dead and their floor must be pinning, which puts the problem back on the
conjunction.
"""
import itertools
import math

import numpy as np

ACCESS = {"HEK293": 0.35, "K562": 0.77, "HUDEP-2": 0.82, "CD34+_HSPC": 0.87}
# stage3.REGION_ENERGY_OFFSETS -- the reachable span, not a per-site value
OFFSETS = (-0.03, 0.05)


def energy(gc, dist, access, offset=0.0):
    return max(0.0, min(1.0, access * (1.8 * gc + 0.6 * math.exp(-dist / 1500.0) + offset)))


def cut_p(cas, e):
    base = 0.86 if cas == "Cas9" else 0.78
    return min(0.99, max(0.4, base + 0.18 * e))


def p_mh(gc):
    return min(0.6, gc * (1 - gc) * 2.2)


def p_hdr_given_cut(cas, e, mh):
    hdr = (0.32 if cas == "Cas9" else 0.24) + 0.35 * e
    return hdr / (hdr + (0.30 if mh else 0.12) + 0.35)


_IM_CACHE = {}


def indel_moments(cas, e, mh_p, cp, n=400000):
    """E[L], Var(L|row) and the conditional MAD of indel_length, by drawing the exact process in
    stage3.sample_indel_length. Closed form is available but the floor()-then-max(1,.) on a gamma
    and an exponential makes it fiddly; sampling is exact enough at n=4e5 and is cached per row
    type, and unlike a closed form it also gives the MAD the oracle MAE needs."""
    key = (cas, round(e, 4), round(mh_p, 4), round(cp, 4))
    if key in _IM_CACHE:
        return _IM_CACHE[key]
    rng = np.random.default_rng(7)
    hdr_base = 0.32 if cas == "Cas9" else 0.24
    hdr = hdr_base + 0.35 * e
    mh = rng.random(n) < mh_p                       # the microhomology coin, drawn first
    mh_nhej = np.where(mh, 0.30, 0.12)
    tot = hdr + mh_nhej + 0.35
    cut = rng.random(n) <= cp                       # stage3: no_cut iff draw > cut_p
    r = rng.random(n) * tot
    mode = np.where(r < hdr, 0, np.where(r < hdr + mh_nhej, 1, 2))   # 0 HDR, 1 MH_NHEJ, 2 BLUNT
    L = np.zeros(n)
    g = np.maximum(1, rng.gamma(2.2, 2.8, n).astype(int))
    x = np.maximum(1, rng.exponential(1 / 0.6, n).astype(int))
    L = np.where(~cut, 0.0, np.where(mode == 0, 0.0, np.where(mode == 1, g, x)))
    m = float(L.mean())
    out = (m, float(L.var()), float(np.abs(L - m).mean()))
    _IM_CACHE[key] = out
    return out


def cons_bound(rows):
    """rows: list of (cas, gc, dist, access, offset, weight). Oracle-regressor cons."""
    r2s, nmaes = [], []
    W = sum(r[5] for r in rows)
    # --- is_cut and is_hdr: Bernoulli with per-row p
    for target in ("is_cut", "is_hdr"):
        ps, ws = [], []
        for cas, gc, dist, acc, off, w in rows:
            e = energy(gc, dist, acc, off)
            cp = cut_p(cas, e)
            if target == "is_cut":
                p = cp
            else:
                mp = p_mh(gc)
                p = cp * (mp * p_hdr_given_cut(cas, e, 1) + (1 - mp) * p_hdr_given_cut(cas, e, 0))
            ps.append(p); ws.append(w)
        mu = sum(p * w for p, w in zip(ps, ws)) / W
        e_var = sum(p * (1 - p) * w for p, w in zip(ps, ws)) / W
        var_p = sum((p - mu) ** 2 * w for p, w in zip(ps, ws)) / W
        tot = e_var + var_p
        r2s.append(var_p / tot if tot > 0 else 0.0)
        mae = sum(2 * p * (1 - p) * w for p, w in zip(ps, ws)) / W
        nmaes.append(mae / math.sqrt(tot) if tot > 0 else 0.0)
    # --- indel_length: continuous, exact conditional moments
    m1s, cvars, mads, ws = [], [], [], []
    for cas, gc, dist, acc, off, w in rows:
        e = energy(gc, dist, acc, off)
        cp = cut_p(cas, e)
        a, v, mad = indel_moments(cas, e, p_mh(gc), cp)
        m1s.append(a); cvars.append(v); mads.append(mad); ws.append(w)
    mu = sum(m * w for m, w in zip(m1s, ws)) / W
    e_var = sum(v * w for v, w in zip(cvars, ws)) / W
    var_p = sum((m - mu) ** 2 * w for m, w in zip(m1s, ws)) / W
    tot = e_var + var_p
    r2s.append(var_p / tot if tot > 0 else 0.0)
    # oracle MAE = the mean conditional MAD about the conditional mean, drawn exactly above
    mae = sum(d * w for d, w in zip(mads, ws)) / W
    nmaes.append(mae / math.sqrt(tot) if tot > 0 else 0.0)
    avg_r2 = sum(r2s) / 3
    avg_nmae = sum(nmaes) / 3
    return 0.7 * max(avg_r2, 0) + 0.3 * (1 - avg_nmae), r2s, nmaes


def homogeneous(cell, cas_frac=0.32, gc=0.50, dist=300):
    acc = ACCESS[cell]
    n9 = round(250 * (1 - cas_frac)); n12 = 250 - n9
    return ([("Cas9", gc, dist, acc, 0.0, 1.0)] * n9 +
            [("Cas12a", gc, dist, acc, 0.0, 1.0)] * n12)


def main():
    print("=== 1. where our 0.101 comes from: the shipped, homogeneous submission ===")
    print("  cell         energy  cons_bound   r2(cut/hdr/indel)        nmae(cut/hdr/indel)")
    for cell in ("CD34+_HSPC", "HUDEP-2", "K562", "HEK293"):
        rows = homogeneous(cell)
        c, r2, nm = cons_bound(rows)
        e = energy(0.50, 300, ACCESS[cell])
        print("  %-11s %6.3f  %8.4f    %s   %s" % (cell, e, c,
              " ".join("%+.3f" % x for x in r2), " ".join("%.3f" % x for x in nm)))
    print("  (measured live: 0.095-0.113 erythroid, 0.076-0.085 HEK293 -- the bound tracks it,")
    print("   which is the check that these formulas are the right ones)")

    print("\n=== 2. energy is CLAMPED on the erythroid types, which kills gc as a lever ===")
    print("  cell        gc=0.30 gc=0.40 gc=0.50 gc=0.60   (distance 300, offset 0)")
    for cell in ("CD34+_HSPC", "HUDEP-2", "K562", "HEK293"):
        print("  %-11s %s" % (cell, "   ".join("%5.3f" % energy(g, 300, ACCESS[cell])
                                               for g in (0.30, 0.40, 0.50, 0.60))))
    print("  -> at accessibility 0.77-0.87, energy hits its 1.0 clamp above gc ~0.33, so every")
    print("     erythroid row sits at cut_p 0.99 / 0.96 regardless of gc. That is exactly why the")
    print("     falsified `gc_spread` sweep moved the floor 0.095 -> 0.099: it varied a clamped")
    print("     quantity. DISTANCE is the only lever that unclamps it.")

    print("\n=== 3. the distribution ceiling: sweep heterogeneity in distance and cas mix ===")
    print("  the most heterogeneous 250 rows the design space allows, per cell")
    print("  cell         arm                                 cons_bound  avg_r2  avg_nmae")
    for cell in ("CD34+_HSPC", "K562", "HEK293"):
        acc = ACCESS[cell]
        arms = {}
        arms["shipped (gc .50, d 300)"] = homogeneous(cell)
        # split the rows across the widest reachable energy span
        for lo, hi in ((300, 3000), (300, 6000), (300, 12000)):
            rows = []
            for i in range(250):
                cas = "Cas9" if i % 100 >= 32 else "Cas12a"
                d = lo if i % 2 else hi
                rows.append((cas, 0.50, d, acc, 0.0, 1.0))
            arms["distance split %d/%d" % (lo, hi)] = rows
        # gc spread as well, for contrast
        rows = []
        for i in range(250):
            cas = "Cas9" if i % 100 >= 32 else "Cas12a"
            rows.append((cas, 0.20 if i % 2 else 0.60, 300, acc, 0.0, 1.0))
        arms["gc split 0.20/0.60"] = rows
        # everything at once
        rows = []
        for i, (cas, gc, d) in enumerate(itertools.islice(itertools.cycle([
                ("Cas9", 0.20, 12000), ("Cas9", 0.60, 300), ("Cas12a", 0.20, 12000),
                ("Cas12a", 0.60, 300)]), 250)):
            rows.append((cas, gc, d, acc, 0.0, 1.0))
        arms["gc AND distance, both extremes"] = rows
        for lab, rows in arms.items():
            c, r2, nm = cons_bound(rows)
            print("  %-11s %-35s %8.4f  %+.3f  %.4f"
                  % (cell, lab, c, sum(r2) / 3, sum(nm) / 3))
    print("\n  target to beat: the field's off-band per-seed value of 0.20-0.25.")
    print("  all-cut's pinned-is_cut value, for reference: 0.237.")


if __name__ == "__main__":
    main()
