#!/usr/bin/env python3
"""conj_k14.py — is the conjunction's band worth pushing to 14, by shrinking the group to 37?

CLAUDE.md computes the wall arithmetic and then declines to measure it: "Reaching k=13 everywhere
needs group ~57, and k=14 needs ~33, both far under the sizes the cas-mix/fidelity work settled
on." This session's own data reproduces that from the other end -- 90 builds at k=11 / group 100 /
narrow 200-seed cut held a surviving Cas12a pool of 199-249, and at ~0.57x decay per band seed that
implies pool(14) = 37-46, i.e. group must be <= 37 for k=14 to form at all.

So the two arms, paired within (task, loop), everything else the cell's live config:

    A  k=11, group 100   the live arm, control
    B  k=14, group  37   the untested arm

Band candidates are the WHOLE 200-seed joined space (200-299 + 400-499) and the cut is the same
200 seeds (narrow), so `band` is a subset of `cut` on both arms. Loops 1..CK_LOOPS give several
paired builds per task via the loop axis; note the arms do NOT cover equal seed space
(loops x k = 55 vs 70 at 5 loops), so COVERAGE is not the comparison -- the per-build quantities
are.

**What decides this is not the band count.** Group size sets the cas mix directly, and CLAUDE.md
measures an empty (mutation x cas x strand) stage-5 cell as roughly a 0.03x multiplier on the whole
score; `build_submission` declines outright at `cells < 8`. So the first result is whether group 37
assembles 8/8 cells at all, the second is what it costs on `weighted x fidelity`, and only the
third is the extra band seats.

Priced from the SHIPPED ROWS exactly as every other arm in this series.

    CK_LOOPS=5 CK_JSON=conj_k14.json python conj_k14.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "conjk14")

import dataclasses, gc, json, logging, time                            # noqa: E402

logging.basicConfig(level=logging.ERROR)
import numpy as np                                                    # noqa: E402

import genExp as G                                                    # noqa: E402
from niome_subnet.genomics import conjunction as CJ                   # noqa: E402
from niome_subnet.genomics import mt19937 as MT                       # noqa: E402
from niome_subnet.genomics.validation import _parse_seeds             # noqa: E402
from conj_stageb import API                                           # noqa: E402
from sd_task import fetch, score                                      # noqa: E402
from widecut_price import records_of                                  # noqa: E402

CELL = os.getenv("CK_CELL", "K562")
LOOPS = int(os.getenv("CK_LOOPS", "5"))
OUT = os.getenv("CK_JSON", "conj_k14.json")
MIN_FREE_GB = float(os.getenv("CK_MIN_FREE_GB", "7"))
BLOCK_A, BLOCK_B = (200, 300), (400, 500)
SPACE = sorted(list(range(*BLOCK_A)) + list(range(*BLOCK_B)))
ARMS = [("A k11/g100", 11, 100), ("B k14/g37", 14, 37)]


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
        print(f"      waiting for memory ({free_gb():.1f} GB free) {tag}", flush=True)
        time.sleep(30)


def compliant(rows, contract, cell_types, ctx, rule):
    ok = CJ.hdr_compliance(records_of(rows), contract, cell_types, ctx, SPACE, rule)
    MT.free_gpu_memory()
    if not ok:
        return []
    keep = set(SPACE)
    for i in sorted(ok):
        keep &= set(int(x) for x in ok[i])
        if not keep:
            break
    return sorted(keep)


def price(rows, contract, reference, cell_types, ctx, seeds):
    band = compliant(rows, contract, cell_types, ctx, "hdr")
    cut = compliant(rows, contract, cell_types, ctx, "cut")
    bs, cs = set(band), set(cut) - set(band)
    per = []
    for s in seeds:
        r = score(rows, contract, reference, cell_types, seed=s)
        per.append({"seed": s, "cons": r["consistency"], "weighted": r["weighted"],
                    "fidelity": r["fidelity"],
                    "regime": "band" if s in bs else ("clean" if s in cs else "dirty")})
    cons = float(np.mean([p["cons"] for p in per]))
    wxf = per[0]["weighted"] * per[0]["fidelity"]
    return {"band_n": len(band), "band": band, "cut_n": len(cut), "clean_n": len(cs),
            "band_hits": [s for s in seeds if s in bs], "clean_hits": [s for s in seeds if s in cs],
            "n_offband_in_cut": len([s for s in seeds if s in set(SPACE) and s not in bs]),
            "per_seed": per, "round_cons": cons, "weighted": per[0]["weighted"],
            "fidelity": per[0]["fidelity"], "wxf": wxf, "round_final": wxf * cons}


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
    print(f"=== {CELL} conjunction | narrow cut {len(SPACE)} ({BLOCK_A[0]}-{BLOCK_A[1]-1} + "
          f"{BLOCK_B[0]}-{BLOCK_B[1]-1}) | band candidates = the whole space | loops 1..{LOOPS} | "
          f"arms " + " vs ".join(f"{n} (k={k} group={g})" for n, k, g in ARMS) +
          f" | {len(tasks)} qualifying tasks ===", flush=True)
    print(f"    live arm is k={base.band_k} group={base.group_size}\n", flush=True)

    mf = max(1, round(base.cas12a_max_fail * len(SPACE) / 900))
    for ti, t in enumerate(tasks, 1):
        tid = (t.get("task_id") or t["id"])[:8]
        if tid in done:
            continue
        contract, reference = dict(t["content"]["contract"]), t["content"]["hbb_reference"]
        seeds = _parse_seeds(contract["seed"])
        created = t.get("created_at", "")
        t0 = time.monotonic()
        c0 = dict(contract, seed=0)
        ctx = G.build_context(c0, reference, cell_types)
        print(f"[{ti}/{len(tasks)}] {tid} {created[:16]} seeds {seeds}", flush=True)
        recs = []
        for name, k, grp in ARMS:
            for loop in range(1, LOOPS + 1):
                wait_for_memory(f"{tid} {name} loop{loop}")
                cfg = dataclasses.replace(base, seed_list=tuple(SPACE), start_seed=SPACE[0],
                                          end_seed=SPACE[-1], cas12a_max_fail=mf,
                                          band_candidates=tuple(SPACE), band_k=k,
                                          group_size=grp, band_loop=loop)
                rows, meta = CJ.build_submission(c0, reference, cell_types, cfg=cfg,
                                                 budget_s=1800)
                MT.free_gpu_memory()
                r = {"arm": name, "k": k, "group": grp, "loop": loop}
                if not rows:
                    r.update(declined=meta.get("reason"), band_reached=meta.get("band"),
                             pool=meta.get("pool"), cells=meta.get("cells"))
                    print(f"    {name} loop{loop}  DECLINED  band reached {meta.get('band')} "
                          f"| {meta.get('reason')}", flush=True)
                else:
                    r.update(price(rows, contract, reference, cell_types, ctx, seeds),
                             pool=meta.get("pool"), cas9=meta.get("cas9_pool"),
                             cells=meta.get("cells"), cas_mix=meta.get("cas_mix"),
                             band_cfg=sorted(int(x) for x in meta.get("band_seeds", [])),
                             elapsed=meta.get("elapsed_s"))
                    print(f"    {name} loop{loop}  band {r['band_n']:3d} clean {r['clean_n']:3d} "
                          f"cells {r['cells']}/8  wxf {r['wxf']:6.1f} "
                          f"(w {r['weighted']:.1f} fid {r['fidelity']:.4f})  "
                          f"hits {len(r['band_hits'])}b/{len(r['clean_hits'])}c  "
                          f"cons {r['round_cons']:.4f} final {r['round_final']:6.1f}  "
                          f"pool {r['pool']} cas9 {r['cas9']}", flush=True)
                del rows; gc.collect()
                recs.append(r)
        out.append({"cell": CELL, "task": tid, "created_at": created, "seeds": seeds,
                    "space": [BLOCK_A, BLOCK_B], "n_space": len(SPACE), "cut": len(SPACE),
                    "loops": LOOPS, "arms": [{"name": n, "k": k, "group": g} for n, k, g in ARMS],
                    "builds": recs})
        json.dump(out, open(OUT, "w"), indent=1)
        print(f"    ({time.monotonic()-t0:.0f}s)\n", flush=True)
        gc.collect()
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
