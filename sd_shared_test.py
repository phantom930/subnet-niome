#!/usr/bin/env python3
"""sd_shared_test.py — does the fleet-shared seed-depend build coordinate correctly?

Four hotkeys reach the same task within ~15s of each other and now share one build. If the lock
logic is wrong the failure is not subtle: either every process builds (12.7 GB each, which is what
the sharing exists to prevent) or every process waits forever and the whole fleet misses the round.

``SD.build`` is stubbed to a 6s sleep so this exercises the coordination, not the builder. Run as
separate PROCESSES, because that is what the miners are — a thread test would share the module
globals and prove less.

    python sd_shared_test.py            # parent: spawns the workers, checks the outcome
    python sd_shared_test.py worker N   # one worker
"""
import json
import os
import subprocess
import sys
import time
from pathlib import Path

CACHE = "data/seed_depend_test"
CONTRACT = {"cell_type": "K562", "seed": 0, "rules": {"max_experiments": 250},
            "active_mutations": ["m1"], "mutation_weights": {"m1": 1.0}}


def worker(idx: int) -> None:
    os.environ["NIOME_INSTANCE"] = f"shtest{idx}"
    sys.argv = ["x"]
    sys.path.insert(0, "/root/workspace/subnet-niome")
    import logging
    logging.basicConfig(level=logging.CRITICAL)
    from niome_subnet.genomics import seed_depend as SD

    built = {"n": 0}

    def fake_build(contract, reference, cell_types, seed=0, cfg=None, budget_s=None):
        built["n"] += 1
        time.sleep(6.0)                       # a build long enough for the others to pile up
        return ([{"experiment_id": f"exp-{i:05d}", "guideRNA": "ACGT"} for i in range(250)],
                {"rule": "mh", "weighted": 1.0, "fidelity": 1.0, "product": 1.0,
                 "elapsed_s": 6.0, "rows": 250, "heavy": 1, "cells": 8, "candidates": 1})

    SD.build = fake_build
    import neurons.miner as M

    class Stub:
        SEED_DEPEND_SHARED = True
        SEED_DEPEND_VARIANT = idx
        SEED_DEPEND_SHARED_VARIANT = 1
        SEED_DEPEND_SEED = 0
        SEED_DEPEND_CACHE_DIR = CACHE
        SEED_DEPEND_LOCK_STALE_S = 1200.0
        SEED_DEPEND_POLL_S = 0.5
        SEED_DEPEND_CACHE_KEEP = 40
        _seed_depend_key = M.Miner._seed_depend_key
        _seed_depend_cached = M.Miner._seed_depend_cached
        _prune_seed_depend_cache = M.Miner._prune_seed_depend_cache
        _build_seed_depend = M.Miner._build_seed_depend

    t0 = time.monotonic()
    rows, meta = Stub()._build_seed_depend(CONTRACT, {}, {}, 120.0, "task-shared-test")
    print(json.dumps({"idx": idx, "rows": len(rows or []), "shared": meta.get("shared"),
                      "built_here": built["n"], "secs": round(time.monotonic() - t0, 1)}))


def main() -> None:
    Path(CACHE).mkdir(parents=True, exist_ok=True)
    for stale in Path(CACHE).glob("*"):
        stale.unlink()
    procs = [subprocess.Popen([sys.executable, __file__, "worker", str(i)],
                              stdout=subprocess.PIPE, text=True,
                              cwd="/root/workspace/subnet-niome") for i in range(4)]
    got = []
    for p in procs:
        out, _ = p.communicate(timeout=180)
        line = [x for x in out.splitlines() if x.startswith("{")]
        got.append(json.loads(line[-1]) if line else {"error": out[-400:]})

    for g in sorted(got, key=lambda x: x.get("idx", 9)):
        print(f"  worker {g.get('idx')}: rows={g.get('rows')} shared={g.get('shared')} "
              f"built_here={g.get('built_here')} {g.get('secs')}s"
              + (f"  ERROR {g['error']}" if "error" in g else ""))

    builds = sum(g.get("built_here", 0) for g in got)
    reads = sum(1 for g in got if g.get("shared") == "read")
    ok = (builds == 1 and reads == 3 and all(g.get("rows") == 250 for g in got))
    print(f"\n  builds={builds} (want 1)   reads={reads} (want 3)   "
          f"all 250 rows={all(g.get('rows') == 250 for g in got)}")
    print("  PASS" if ok else "  FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "worker":
        worker(int(sys.argv[2]))
    else:
        main()
