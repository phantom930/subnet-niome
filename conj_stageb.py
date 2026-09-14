#!/usr/bin/env python3
"""conj_stageb.py — stage B: measure the conjunction's per-seed REGIMES and price them.

Stage A (`conj_grid.py`) ranks configs on quantities that need no seed sampling. What it cannot
give is the thing that decides payout: what an off-band seed is actually worth. The conjunction has
three regimes rather than two --

    band   (k seeds)          all three stage-4 targets pinned  -> exactly 1.000
    clean  (`clean` - k)      `is_cut` pinned                   -> ~0.24-0.30, measured here
    rest   (the remainder)    nothing pinned                    -> ~0.10, measured here

-- and a round averages three independent draws, so the expectation is an exact enumeration over
the 10 ways three seeds split across the regimes. No Monte-Carlo error enters; the only measurement
error is in `v_clean` / `v_rest`, which is why this samples 20 seeds per regime rather than the 3
earlier sweeps used. At 3 samples `conj_bandwidth.json` reported `e_final` 33.8-45.8 across arms
whose `w x fid` was flat within 1%, i.e. the ranking was sampling noise.

Each arm is priced against the REAL fields of its own cell type -- the rank-10 cutoff and the real
`SCORE_DISTRIBUTION` position, one field at a time -- because a construction's score and its field
move together and pooling across contracts prices field softness instead of the construction.

    CSB_CELL=K562 python conj_stageb.py
    CSB_CELL=HEK293 CSB_TOP=10 CSB_NS=12 python conj_stageb.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")

import dataclasses as _dc            # noqa: E402
import json                          # noqa: E402
import logging                       # noqa: E402
import os.path                       # noqa: E402
import random                        # noqa: E402
import statistics as st              # noqa: E402
import time                          # noqa: E402
import urllib.request                # noqa: E402
from collections import defaultdict  # noqa: E402
from math import factorial           # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np                   # noqa: E402

CELL = os.getenv("CSB_CELL", "K562")
os.environ.setdefault("NIOME_INSTANCE", "csb_" + CELL.replace("+", "").replace("-", "_"))

import genExp as G                                          # noqa: E402
from niome_subnet.genomics import all_cut as AC             # noqa: E402
from niome_subnet.genomics import fastgreedy as FG          # noqa: E402
from niome_subnet.genomics import mt19937 as MT             # noqa: E402
from conj_test import hdr_compliance                        # noqa: E402
from conj_joined import JoinedCut                           # noqa: E402
from conj_feas import greedy_prefix                         # noqa: E402
from conj_grid import sub_window                            # noqa: E402
from joined300 import fetch_tasks, pick_tasks               # noqa: E402
from sd_task import score, fetch                            # noqa: E402

CLASSES = [(100, 199), (400, 499), (700, 799)]
TOP = int(os.getenv("CSB_TOP", "12"))
NS = int(os.getenv("CSB_NS", "20"))
POOL_TARGET = int(os.getenv("CSB_POOL_TARGET", "500"))
NSEED = 900
DIST = [0.30, 0.20, 0.20, 0.15, 0.05, 0.03, 0.025, 0.02, 0.015, 0.01]
API = "https://niome-api.genomes.io/api/v3"
# CLAUDE.md's live per-cell all-HDR `w x fid`, kept only to show how far a CROSS-CONTRACT
# comparison misleads. The real baseline is built on THIS contract by `all_hdr_arm` below --
# stage A already shows why: the conjunction reaches w x fid 308 on K562's contract and 202 on
# HUDEP-2's, a 52% gap that is the contract, not the construction. Comparing either against a
# number measured on other contracts prices field softness, which is this file's own standing rule.
AH_LIVE = {"K562": 220.1, "HUDEP-2": 298.3, "CD34+_HSPC": 233.3, "HEK293": 179.0}
# Every hotkey of ours, by ss58: uids are recycled, so identifying our own rows by uid
# contaminates the field with another operator's submissions.
OURS = {
    "5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW", "5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb",
    "5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU", "5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf",
    "5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2", "5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn",
    "5H5v45M2i6cFtS3Di4abh2zjWrPmuJZzNPoAYtVAXVFa6sp3", "5CS8FdHr8Ddv14QHr75E5K7e9xv3wre9m66ywRHZMR8zhntN",
    "5Cr9gJ3ukDDdxnhpyRM58gdMz7Sjuw8u58iiuYphVz1ZUGD2", "5Ge6UdGJE7BabD3bv8L9AR79rW47K3XQetcTHgaLD9NF2Vf5",
    "5D7H8V67nQiFm8K8d7U6dAtZnD5XeKkbv9VJ8GkqyeUXRZzd",
}


def fields(cell):
    """Every current-regime 3-seed field this cell played, ours removed, as sorted finals."""
    sc = fetch(f"{API}/miners/scores?limit=40000", cache="sd_task_scores.json")
    sc = sc if isinstance(sc, list) else (sc.get("data") or sc.get("items") or [])
    tk = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    tk = tk if isinstance(tk, list) else (tk.get("items") or tk.get("data") or [])
    meta = {}
    for t in tk:
        c = (t.get("content") or {}).get("contract", {})
        s = str(c.get("seed", "") or "")
        meta[t.get("task_id") or t.get("id")] = (
            c.get("cell_type"), len([x for x in s.split(",") if x.strip().isdigit()]))
    best = defaultdict(dict)
    for r in sc:
        c, n = meta.get(r.get("task_id"), (None, 0))
        if c != cell or n != 3 or r.get("miner_hotkey") in OURS:
            continue
        hk, f = r["miner_hotkey"], float(r.get("final_score") or 0)
        if hk not in best[r["task_id"]] or f > best[r["task_id"]][hk]:
            best[r["task_id"]][hk] = f
    return [sorted(v.values(), reverse=True) for v in best.values()
            if len(v) >= 10 and max(v.values()) > 0]


def share(final, field):
    rank = sum(1 for f in field if f > final) + 1
    return DIST[rank - 1] if rank <= 10 else 0.0


def price(band, clean_excl_band, v_clean, v_rest, wxf, fs):
    """Exact E[curve share]: enumerate the 10 ways three seeds split across the three regimes."""
    pb, pc = band / NSEED, clean_excl_band / NSEED
    pr = max(0.0, 1.0 - pb - pc)
    tot = p_place = 0.0
    for nb in range(4):
        for nc in range(4 - nb):
            nr = 3 - nb - nc
            w = (factorial(3) / (factorial(nb) * factorial(nc) * factorial(nr))
                 * pb ** nb * pc ** nc * pr ** nr)
            if w <= 0:
                continue
            final = wxf * (nb + nc * v_clean + nr * v_rest) / 3.0
            tot += w * sum(share(final, f) for f in fs) / len(fs)
            p_place += w * sum(1 for f in fs if share(final, f) > 0) / len(fs)
    return tot, p_place


def hdr_clean_seeds(lo=100, hi=999):
    """Which seeds the whole submission repairs by HDR on -- all-HDR's band, hunted not assumed.

    Reads `valid_experiments.json`, which the preceding `score()` call has already written for this
    row set: `assemble` returns SUBMISSION rows and `stage3.simulate` needs stage 12's `features`
    block, so the band cannot be scanned off the rows directly. Pure CPU, ~900 x 250 draws.
    """
    from niome_subnet.genomics.validation import stage3 as S3
    from niome_subnet.utils import settings as ST
    valid = json.load(open(ST.VALID_EXPERIMENTS_PATH))
    band = []
    for seed in range(lo, hi + 1):
        if all(S3.simulate(e, seed)["outcome"] == "HDR" for e in valid):
            band.append(seed)
    return band


def all_hdr_arm(contract, reference, cell_types, ctx, joined, fs, rng, ns):
    """The incumbent, built on THIS contract and the SAME joined window, then priced identically.

    Without this the comparison is cross-contract: CLAUDE.md's live `w x fid` per cell was measured
    on the rounds the fleet actually played, and `weighted` moves up to 54% with the contract.
    """
    from niome_subnet.genomics import all_hdr as AH
    t0 = time.monotonic()
    rows, meta = AH.build_for_cell(contract, reference, cell_types, seed_list=joined)
    if rows is None:
        print(f"  matched all-HDR declined: {meta.get('reason')}")
        return None
    # One score() first: it writes valid_experiments.json, which the band hunt reads.
    bs = score(rows, contract, reference, cell_types, seed=500)
    band = hdr_clean_seeds()
    if not band:
        print("  matched all-HDR built but has no HDR-clean seed in 100-999")
        return None
    off = [s for s in range(100, 1000) if s not in set(band)]
    fl = [score(rows, contract, reference, cell_types, seed=s)["consistency"]
          for s in rng.sample(off, min(ns, len(off)))]
    bs = score(rows, contract, reference, cell_types, seed=int(band[0]))
    wxf = bs["weighted"] * bs["fidelity"]
    e, p = price(len(band), 0, 0.0, st.mean(fl), wxf, fs)
    print(f"  matched all-HDR: band {len(band)}  w x fid {wxf:.1f}  floor {st.mean(fl):.4f}  "
          f"band_cons {bs['consistency']:.4f}  E[share] {e:.5f}  P(place) {p:.2%}  "
          f"({time.monotonic()-t0:.0f}s)")
    return {"band": len(band), "wxfid": wxf, "floor": st.mean(fl), "floor_sd": st.pstdev(fl),
            "band_cons": bs["consistency"], "e_share": e, "p_place": p, "meta_reason": None}


def shortlist(rows):
    """Top `TOP` by w x fid, but guaranteeing every k and every group is represented once.

    `CSB_FORCE` ("k:group:width:light,...", light "n" for None) prices exactly those arms instead.
    Needed because selecting by `w x fid` is ANTI-correlated with k -- the k axis falls in w x fid
    while E[share] rises with band size -- so the natural shortlist misses the high-k arms at the
    best group, which is where the conjunction is strongest.
    """
    ok = [r for r in rows if "wxfid" in r]
    force = os.getenv("CSB_FORCE", "").strip()
    if force:
        picked = []
        for spec in force.split(","):
            k, g, w, li = spec.split(":")
            li = None if li.lower() in ("n", "none", "") else int(li)
            m = [r for r in ok if r["k"] == int(k) and r["group"] == int(g)
                 and r["width"] == int(w) and r["light"] == li]
            if m:
                picked.append(m[0])
            else:
                print(f"  forced arm {spec} not in the grid, skipped")
        return picked
    ok.sort(key=lambda r: -r["wxfid"])
    picked, seen = [], set()
    for keyf in (lambda r: ("k", r["k"]), lambda r: ("g", r["group"])):
        for r in ok:
            if keyf(r) not in seen:
                seen.add(keyf(r))
                if r not in picked:
                    picked.append(r)
    for r in ok:
        if len(picked) >= TOP:
            break
        if r not in picked:
            picked.append(r)
    return picked


def main():
    tag = CELL.replace("+", "").replace("-", "_")
    stem = f"conj_stageb_{tag}" + ("_forced" if os.getenv("CSB_FORCE", "").strip() else "")
    grid = f"conj_grid_{tag}.json"
    rows = shortlist(json.load(open(grid)))
    fs = fields(CELL)
    cuts = sorted(f[9] for f in fs)
    print(f"{CELL}: {len(fs)} current-regime 3-seed fields (ours excluded)")
    print(f"rank-10 cutoff  min {cuts[0]:.1f}  p25 {cuts[len(cuts)//4]:.1f}  "
          f"median {cuts[len(cuts)//2]:.1f}  p75 {cuts[3*len(cuts)//4]:.1f}  max {cuts[-1]:.1f}")
    print(f"pricing {len(rows)} shortlisted configs at {NS} samples per regime\n")

    # PIN the task to the one stage A built on. `pick_tasks` returns the newest STAMPED task, and
    # rounds stamp every ~2h24m -- long enough that a grid and its pricing run can straddle one.
    # Measured: the HUDEP-2 grid built on bbdece32 at 05:52, 06b78fe1 stamped at ~07:18, and stage
    # B at 08:12 then asked for a Cas12a bank that was never built. Silent drift would be worse
    # than the crash: both arms must be priced on the SAME contract as each other AND as stage A.
    want = next((r["task"] for r in json.load(open(grid)) if r.get("task")), None)
    items = fetch_tasks()
    task = next((t for t in items
                 if (t.get("task_id") or t.get("id", "")).startswith(want)), None)
    if task is None:
        raise SystemExit(f"stage A built {CELL} on task {want}, which is no longer in the listing")
    newest = pick_tasks(items).get(CELL)
    newest_id = (newest.get("task_id") or newest["id"])[:8] if newest else "-"
    print(f"task {want} (pinned from {grid}"
          + (f"; newest stamped is now {newest_id})" if newest_id != want else ")"))
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
    ok = hdr_compliance(records, contract, cell_types, ctx, joined)
    rng = random.Random(20260913)
    jset = set(joined)
    offw = [s for s in range(100, 1000) if s not in jset]

    out = []
    cache = {}
    print(f"{'k':>3}{'grp':>5}{'wid':>5}{'light':>7}{'clean':>7}{'w x fid':>9}"
          f"{'v_clean':>9}{'v_rest':>8}{'E[share]':>10}{'P(place)':>10}{'s':>6}")
    for r in rows:
        t0 = time.monotonic()
        width, k, g = r["width"], r["k"], r["group"]
        ck = (width, k, g)
        if ck not in cache:
            cand = sub_window(joined, 700, width)
            okw = {i: (v & set(cand)) for i, v in ok.items()}
            steps = greedy_prefix(okw, k, cand, g)
            band, alive = steps[k]
            pool = [records[i] for i in alive]
            sel = FG.FastGreedy(pool, window_lo=min(joined), window_hi=max(joined),
                                seeds=np.asarray(joined, dtype=np.int64))
            idx, _u = sel.best(g, restarts=12)
            grp = [pool[i] for i in idx]
            bad = set()
            for rec in grp:
                bad.update(int(x) for x in rec["fails"])
            clean = sorted(jset - bad)
            cfg_scan = _dc.replace(cfg0, group_size=g, pool_target=POOL_TARGET)
            cas9 = AC.scan_cas9(np.array(clean, dtype=np.int64), contract, cell_types, ctx,
                                sites, cfg_scan, n_rows - g)
            if band and cas9:
                ok9 = hdr_compliance(cas9, contract, cell_types, ctx, band)
                cas9 = [cas9[i] for i in sorted(ok9) if len(ok9[i]) == len(band)]
            cache[ck] = (band, grp, cas9, clean, bad)
            MT.free_gpu_memory()
        band, grp, cas9, clean, bad = cache[ck]
        cfg = _dc.replace(cfg0, group_size=g, light_cell_rows=r["light"])
        rw = AC.assemble(grp, cas9, contract, ctx, cfg, n_rows)
        pc = [s for s in clean if s not in set(band)]
        rest = sorted(bad) + offw
        vc = [score(rw, contract, reference, cell_types, seed=s)["consistency"]
              for s in rng.sample(pc, min(NS, len(pc)))] or [0.0]
        vr = [score(rw, contract, reference, cell_types, seed=s)["consistency"]
              for s in rng.sample(rest, min(NS, len(rest)))]
        bs = score(rw, contract, reference, cell_types, seed=int(band[0]))
        wxf = bs["weighted"] * bs["fidelity"]
        e, p = price(len(band), len(pc), st.mean(vc), st.mean(vr), wxf, fs)
        dt = time.monotonic() - t0
        print(f"{k:>3}{g:>5}{width:>5}{str(r['light']):>7}{len(clean):>7}{wxf:>9.2f}"
              f"{st.mean(vc):>9.4f}{st.mean(vr):>8.4f}{e:>10.5f}{p:>9.2%}{dt:>6.0f}")
        out.append({**{x: r[x] for x in ("cell", "task", "k", "group", "width", "light")},
                    "clean": len(clean), "band": len(band), "wxfid": wxf,
                    "v_clean": st.mean(vc), "v_clean_sd": st.pstdev(vc), "n_clean": len(vc),
                    "v_rest": st.mean(vr), "v_rest_sd": st.pstdev(vr), "n_rest": len(vr),
                    "band_cons": bs["consistency"], "e_share": e, "p_place": p})
        json.dump(out, open(stem + ".json", "w"), indent=1)
    print("\nmatched baseline, same contract and same joined window:")
    ah = all_hdr_arm(contract, reference, cell_types, ctx, joined, fs, rng, NS)
    ah_e = ah["e_share"] if ah else 0.0
    print(f"  (for reference, this cell's LIVE all-HDR w x fid is {AH_LIVE[CELL]} -- measured on "
          f"other contracts, so not comparable to the arms above)")
    out.sort(key=lambda r: -r["e_share"])
    b = out[0]
    print(f"BEST {CELL}: k={b['k']} group={b['group']} width={b['width']} light={b['light']}  "
          f"E[share] {b['e_share']:.5f}  P(place) {b['p_place']:.2%}  "
          f"w x fid {b['wxfid']:.1f}  clean {b['clean']}  v_clean {b['v_clean']:.3f}")
    print(f"  vs matched all-HDR: {b['e_share']/ah_e:.2f}x" if ah_e > 0 else
          "  vs matched all-HDR: the incumbent earns 0 on these fields")
    json.dump({"cell": CELL, "arms": out, "all_hdr_matched": ah, "all_hdr_live_wxfid": AH_LIVE[CELL]},
              open(stem + ".json", "w"), indent=1)


if __name__ == "__main__":
    main()
