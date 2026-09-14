#!/usr/bin/env python3
"""conj_mc.py — reprice the replication with PER-SEED sampling instead of regime means.

`price()` in conj_stageb.py values each regime at its mean, so a round takes one of ten fixed point
values. That is fine when a round sits clear of the cutoff and wrong when it straddles one: the
curve is a step function, so a round modelled two points short scores zero where in reality it
places about half the time. This resamples every seed from its regime's measured distribution and
Monte-Carlos the round -- the same expectation, with the straddle handled correctly.

Two regimes carry real spread and they differ in size (stage B, 20 samples x 18 arms per cell):
the dirty floor has sd ~0.012, the cut-clean regime ~0.045 -- 3.6x wider. So this matters far more
for the conjunction, which HAS a clean regime, than for all-HDR, which does not. The band regime is
exactly 1.0 and carries no spread (`band_cons` verified 1.0000 on every arm).

Note what this does NOT rescue: all-HDR's k=1 round falls 4.5-15.8 points short of these contracts'
own cutoffs, against a per-seed sd worth 1.3-1.6 points -- a 3.3-5.6 sigma gap. Sampling cannot
close it, so the point model was not biased there. The open question this answers is the opposite
one: whether the CONJUNCTION's wider clean regime was being undervalued.

The replication stored regime means only, so spreads are imputed per cell from stage B's arms of
the same construction on the same joined window. Sampling is Gaussian truncated to [0, 1]; the
underlying values are bounded and mildly skewed, so read sub-5% differences as model choice.

    python conj_mc.py
    CMC_N=500000 python conj_mc.py
"""
import glob
import json
import os
import statistics as st
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "conjmc")

import numpy as np

from conj_replicate import fields_by_task          # noqa: E402

N = int(os.getenv("CMC_N", "200000"))
DIST = np.array([0.30, 0.20, 0.20, 0.15, 0.05, 0.03, 0.025, 0.02, 0.015, 0.01])
RNG = np.random.default_rng(int(os.getenv("CMC_SEED", "20260914")))
CELLS = [("HEK293", (6, 80, 300, 12)), ("K562", (8, 80, 150, 6)),
         ("CD34+_HSPC", (8, 80, 100, 12)), ("HUDEP-2", (8, 80, 225, 12))]


def spreads(cell):
    tag = cell.replace("+", "").replace("-", "_")
    vc, vr, fl = [], [], []
    for p in glob.glob(f"conj_stageb_{tag}*.json"):
        d = json.load(open(p))
        arms = d["arms"] if isinstance(d, dict) else d
        vc += [a["v_clean_sd"] for a in arms]
        vr += [a["v_rest_sd"] for a in arms]
        ah = d.get("all_hdr_matched") if isinstance(d, dict) else None
        if ah and ah.get("floor_sd") is not None:
            fl.append(ah["floor_sd"])
    return (st.mean(vc) if vc else 0.045, st.mean(vr) if vr else 0.012,
            st.mean(fl) if fl else (st.mean(vr) if vr else 0.012))


def mc_finals(band, clean_excl, v_clean, sd_clean, v_rest, sd_rest, wxf, n=N):
    p_b, p_c = band / 900.0, clean_excl / 900.0
    u = RNG.random((n, 3))
    is_b, is_c = u < p_b, (u >= p_b) & (u < p_b + p_c)
    v = np.clip(RNG.normal(v_rest, sd_rest, (n, 3)), 0.0, 1.0)
    if clean_excl:
        v = np.where(is_c, np.clip(RNG.normal(v_clean, sd_clean, (n, 3)), 0.0, 1.0), v)
    v = np.where(is_b, 1.0, v)
    return wxf * v.mean(axis=1)


def e_share(finals, fields):
    tot = 0.0
    for f in fields:
        asc = np.sort(np.asarray(f, dtype=float))
        rank = len(asc) - np.searchsorted(asc, finals, side="right") + 1
        tot += np.where(rank <= 10, DIST[np.clip(rank - 1, 0, 9)], 0.0).mean()
    return tot / len(fields)


def main():
    print(f"Monte Carlo repricing: {N:,} rounds per arm per contract, per-seed regime sampling\n")
    rows = []
    for cell, arm in CELLS:
        sc, sr, sf = spreads(cell)
        tag = cell.replace("+", "").replace("-", "_")
        d = json.load(open(f"conj_replicate_{tag}.json"))
        ftasks = fields_by_task(cell)
        pool = list(ftasks.values())
        by_short = {t[:8]: v for t, v in ftasks.items()}
        pc = pb = pcp = pbp = 0.0
        mc_c = mc_b = mcp_c = mcp_b = 0.0
        wins_pt = wins_mc = 0
        for r in d:
            a = next(x for x in r["arms"]
                     if (x["k"], x["group"], x["width"], x["light"]) == arm)
            ah = r["all_hdr"]
            own = [by_short[r["task"]]]
            fc = mc_finals(a["band"], a["clean"] - a["band"], a["v_clean"], sc, a["v_rest"], sr,
                           a["wxfid"])
            fb = mc_finals(ah["band"], 0, 0.0, 0.0, ah["floor"], sf, ah["wxfid"])
            c_own, b_own = e_share(fc, own), e_share(fb, own)
            c_pool, b_pool = e_share(fc, pool), e_share(fb, pool)
            pc += a["e_own"]; pb += ah["e_own"]; pcp += a["e_pool"]; pbp += ah["e_pool"]
            mc_c += c_own; mc_b += b_own; mcp_c += c_pool; mcp_b += b_pool
            wins_pt += int(a["e_own"] > ah["e_own"])
            wins_mc += int(c_own > b_own)
        rows.append((cell, arm, sc, sr, sf, pc / pb if pb else 0, mc_c / mc_b if mc_b else 0,
                     pcp / pbp if pbp else 0, mcp_c / mcp_b if mcp_b else 0,
                     wins_pt, wins_mc, len(d)))
    print(f"{'cell':<12}{'sdCln':>7}{'sdFlr':>7} | {'own pt':>8}{'own MC':>8}{'win pt':>8}"
          f"{'win MC':>8} | {'pool pt':>9}{'pool MC':>9}")
    for (cell, arm, sc, sr, sf, o_pt, o_mc, p_pt, p_mc, wp, wm, n) in rows:
        print(f"{cell:<12}{sc:>7.3f}{sf:>7.3f} | {o_pt:>7.2f}x{o_mc:>7.2f}x{wp:>6}/{n}"
              f"{wm:>6}/{n} | {p_pt:>8.2f}x{p_mc:>8.2f}x")
    json.dump([{"cell": c, "arm": list(a), "sd_clean": sc, "sd_rest": sr, "sd_floor": sf,
                "own_point": float(op), "own_mc": float(om),
                "pool_point": float(pp), "pool_mc": float(pm),
                "wins_point": int(wp), "wins_mc": int(wm), "n": int(n)}
               for (c, a, sc, sr, sf, op, om, pp, pm, wp, wm, n) in rows],
              open("conj_mc.json", "w"), indent=1)
    print("\nwrote conj_mc.json")


if __name__ == "__main__":
    main()
