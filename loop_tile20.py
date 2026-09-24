#!/usr/bin/env python3
"""loop_tile20.py — TWENTY hotkeys on ONE shared 200-seed window, separated by the LOOP AXIS.

The arm under test:

    joined window / cut   900 (wide)
    band space            200 seeds: 200-299 + 400-499 (a JOINED, non-contiguous space)
    sub window            200 seeds  -- the WHOLE space
    stride                0          -- every hotkey searches the IDENTICAL candidate set
    loops                 1..20, one per hotkey
    everything else       HEK293's live CELL_CONFIG (k, group, light, cell-aware)

This is the direct counterpart to `loop_tile8.py`, which gave the same two classes to eight
hotkeys at stride 25 and let WINDOW OFFSETS do the separating. There the 4x window overlap
duplicated ~25% of the band seats. Here the windows are identical and the LOOP AXIS does the
separating instead: loop L bans every seed banded by loops 1..L-1, so the bands are disjoint by
construction and the union is exactly headcount x k until the chain runs out of candidates.

Twenty loops at k=8 asks for 160 of the 200 seeds, so this also measures where the chain BREAKS --
`choose_band` is a greedy prefix and the surviving pool decays ~P(hdr) per band seed, so a late
loop can fail to form even with candidates left. Once it does, every later loop collapses onto the
same remainder (that is what `band_loop_short` records in the live builder), and this script
reproduces that behaviour exactly rather than declining.

**Tasks are selected, not sampled**: only HEK293 rounds whose three seeds include at least one in
200-299 AND at least one in 400-499. Chance puts a round in that set 6.6% of the time, so every
rate here is P(... | both classes drawn), never P(...).

    LT_JSON=loop_tile20.json python loop_tile20.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "looptile20")

import dataclasses, gc, json, logging, time                            # noqa: E402

logging.basicConfig(level=logging.ERROR)
import numpy as np                                                    # noqa: E402

import genExp as G                                                    # noqa: E402
import joined_window as JW                                            # noqa: E402
from niome_subnet.genomics import conjunction as CJ                   # noqa: E402
from niome_subnet.genomics import mt19937 as MT                       # noqa: E402
from niome_subnet.genomics.validation import _parse_seeds             # noqa: E402
from conj_stageb import API                                           # noqa: E402
from sd_task import fetch                                             # noqa: E402

CELL = os.getenv("LT_CELL", "HEK293")
NHK = int(os.getenv("LT_HK", "20"))
OUT = os.getenv("LT_JSON", "loop_tile20.json")
# "wide": the Cas12a cut min-unions over the full 100-999. "narrow": it min-unions over the SAME
# 200 seeds the band is drawn from, so `band` is a subset of `cut` by construction. A different cut
# is a different `bank_key`, hence a different bank, different survivors and a different band --
# the two arms are not the same build with one number changed.
CUTMODE = os.getenv("LT_CUT", "wide")
assert CUTMODE in ("wide", "narrow"), CUTMODE
MIN_FREE_GB = float(os.getenv("LT_MIN_FREE_GB", "7"))
BLOCK_A = (200, 300)
BLOCK_B = (400, 500)
SPACE = list(range(*BLOCK_A)) + list(range(*BLOCK_B))
V_CLEAN = {"HEK293": 0.212, "K562": 0.212, "HUDEP-2": 0.212, "CD34+_HSPC": 0.212}
V_DIRTY = {"HEK293": 0.082, "K562": 0.104, "HUDEP-2": 0.104, "CD34+_HSPC": 0.104}


def free_gb() -> float:
    m = {}
    for line in open("/proc/meminfo"):
        p = line.split()
        m[p[0].rstrip(":")] = int(p[1])
    return m.get("MemAvailable", 0) / 1048576.0


def wait_for_memory(tag: str = "") -> None:
    for _ in range(120):
        if free_gb() >= MIN_FREE_GB:
            return
        print(f"    waiting for memory ({free_gb():.1f} GB free) {tag}", flush=True)
        time.sleep(30)


def main():
    cell_types = G.fetch_cell_types()
    G.load_sequence()
    items = fetch(f"{API}/tasks?limit=500", cache="sd_task_listing.json")
    items = items if isinstance(items, list) else (items.get("items") or items.get("data") or [])
    A, B = set(range(*BLOCK_A)), set(range(*BLOCK_B))
    tasks = []
    for t in sorted(items, key=lambda x: x.get("created_at", ""), reverse=True):
        c = (t.get("content") or {}).get("contract") or {}
        if c.get("cell_type") != CELL:
            continue
        sd = _parse_seeds(c.get("seed"))
        if len(sd) != 3 or not (any(s in A for s in sd) and any(s in B for s in sd)):
            continue
        tasks.append(t)

    base = CJ.config_for(CELL)
    out = json.load(open(OUT)) if os.path.exists(OUT) else []
    done = {r["task"] for r in out}
    print(f"=== {CELL} | {NHK} hotkeys, loops 1..{NHK}, ONE shared window "
          f"{BLOCK_A[0]}-{BLOCK_A[1]-1} + {BLOCK_B[0]}-{BLOCK_B[1]-1} ({len(SPACE)} seeds), "
          f"stride 0 | cut {900 if CUTMODE == 'wide' else len(SPACE)} {CUTMODE} | "
          f"k={base.band_k} group={base.group_size} | "
          f"{len(tasks)} qualifying tasks | seats {NHK * base.band_k} of {len(SPACE)} ===",
          flush=True)

    cut = list(JW.FULL_SPACE) if CUTMODE == "wide" else sorted(SPACE)
    cutset = set(cut)
    for ti, t in enumerate(tasks, 1):
        tid = (t.get("task_id") or t["id"])[:8]
        if tid in done:
            continue
        contract, reference = dict(t["content"]["contract"]), t["content"]["hbb_reference"]
        seeds = _parse_seeds(contract["seed"])
        created = t.get("created_at", "")
        t0 = time.monotonic()
        mf = max(1, round(base.cas12a_max_fail * len(cut) / 900))
        cfg0 = dataclasses.replace(base, seed_list=tuple(cut), start_seed=cut[0],
                                   end_seed=cut[-1], cas12a_max_fail=mf)
        wait_for_memory(tid)
        ctx = G.build_context(dict(contract, seed=0), reference, cell_types)
        sites = G.enumerate_sites(ctx, 3000, (20, 23))
        path = os.path.join(CJ.BANK_DIR,
                            f"cas12a-{CJ.bank_key(dict(contract, seed=0), cell_types, cfg0)}.npz")
        if not os.path.exists(path):
            os.makedirs(CJ.BANK_DIR, exist_ok=True)
            bank = CJ.build_bank(dict(contract, seed=0), reference, cell_types, ctx, sites,
                                 cfg0, None)
            MT.free_gpu_memory()
            if not bank:
                print(f"  {tid} bank empty, skipped", flush=True); continue
            CJ.save_bank(path, bank); del bank; gc.collect()
        records = CJ.load_bank(path, limit=cfg0.bank_keep)
        n_bank = len(records)
        ok = CJ.hdr_compliance(records, dict(contract, seed=0), cell_types, ctx, SPACE,
                               cfg0.band_rule)
        MT.free_gpu_memory()
        ok = {j: np.fromiter(v, dtype=np.int32) for j, v in ok.items()}
        del records; gc.collect()
        cell_ok = None
        if cfg0.band_cell_aware:
            cell_ok = CJ.cas9_cell_probe(dict(contract, seed=0), cell_types, ctx, sites, cfg0,
                                         SPACE, cfg0.band_rule)
            MT.free_gpu_memory()

        # ---- the loop axis, reproduced exactly as `build_submission` runs it -----------------
        # Loop L bans loops 1..L-1. While the chain is intact that is one fresh choose_band per
        # loop; once a loop fails to form (short band, or survivors under the group) every later
        # loop collapses onto the same remainder, which is what `band_loop_short` records live.
        hks, banned, chain_ok, short_at, short_res = [], set(), True, None, None
        for i in range(NHK):
            loop = i + 1
            if chain_ok:
                left = [s for s in SPACE if s not in banned]
                if len(left) < cfg0.band_k:
                    chain_ok = False
                    short_at = loop
                    short_res = {"band": [], "alive": 0, "declined": True, "n_left": len(left)}
                else:
                    band, alive = CJ.choose_band(ok, cfg0.band_k, left, cfg0.group_size, cell_ok,
                                                 {}, cfg0.band_quota)
                    b = sorted(int(x) for x in band)
                    res = {"band": b, "alive": len(alive), "declined": False, "n_left": len(left)}
                    if len(b) < cfg0.band_k or len(alive) < cfg0.group_size:
                        chain_ok = False
                        short_at = loop
                        short_res = res
                    else:
                        banned.update(b)
            r = dict(short_res) if not chain_ok else res
            hit = [s for s in seeds if s in set(r["band"])]
            hks.append({"hk": i, "loop": loop, "band": r["band"], "alive": r["alive"],
                        "formed": len(r["band"]) == cfg0.band_k and not r["declined"],
                        "n_left": r["n_left"], "shared_from": short_at if not chain_ok else None,
                        "declined": r["declined"], "hits": hit, "builds": []})
        del ok; gc.collect()

        # ---- full builds, one per DISTINCT hitting band --------------------------------------
        seen: dict[tuple, dict] = {}
        for h in hks:
            if not h["hits"] or not h["band"]:
                continue
            key = tuple(h["band"])
            if key in seen:
                h["builds"].append(dict(seen[key], reused=True)); continue
            wait_for_memory(f"{tid} h{h['hk']}")
            bcfg = dataclasses.replace(cfg0, band_candidates=key)
            rows, meta = CJ.build_submission(dict(contract, seed=0), reference, cell_types,
                                             cfg=bcfg, budget_s=1800)
            MT.free_gpu_memory()
            if not rows:
                rec = {"declined": meta.get("reason")}
            else:
                got = sorted(int(x) for x in meta["band_seeds"])
                clean = set(int(x) for x in meta["clean_seeds"])
                bh = [s for s in seeds if s in set(got)]
                ch = [s for s in seeds if s in clean and s not in set(got)]
                off = [s for s in seeds if s not in set(got)]
                vc = V_CLEAN.get(CELL, 0.212); vd = V_DIRTY.get(CELL, 0.104)
                rec = {"band": got, "band_match": got == h["band"], "clean_n": len(clean),
                       "band_hits": bh, "clean_hits": ch, "n_offband": len(off),
                       "n_offband_in_cut": len([s for s in off if s in cutset]),
                       "cut_n": len(cut),
                       "cons": (len(bh) + len(ch) * vc + (len(off) - len(ch)) * vd) / 3,
                       "pool": meta.get("pool"), "cas9": meta.get("cas9_pool")}
            seen[key] = rec
            h["builds"].append(rec)
            gc.collect()

        formed = [h for h in hks if h["formed"]]
        allb = sorted({s for h in hks for s in h["band"]})
        caught = sorted({s for h in hks for s in h["hits"]})
        in_space = [s for s in seeds if s in set(SPACE)]
        print(f"  [{ti}/{len(tasks)}] {tid} {created[:16]} seeds {seeds} (in space {in_space})  "
              f"bank {n_bank}  loops formed {len(formed)}/{NHK}"
              f"{'' if short_at is None else f' short at {short_at}'}  "
              f"union {len(allb)}/{len(SPACE)}  caught {len(caught)}/{len(in_space)} {caught}  "
              f"({time.monotonic()-t0:.0f}s)", flush=True)
        out.append({"cell": CELL, "task": tid, "created_at": created, "seeds": seeds,
                    "bank": n_bank, "space": [BLOCK_A, BLOCK_B], "n_space": len(SPACE),
                    "stride": 0, "loops": NHK, "k": cfg0.band_k, "group": cfg0.group_size,
                    "cut": len(cut), "cut_mode": CUTMODE, "hotkeys": hks, "band_union": len(allb),
                    "loops_formed": len(formed), "short_at": short_at, "caught": caught})
        json.dump(out, open(OUT, "w"), indent=1)
        gc.collect()
    print(f"\nwrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
