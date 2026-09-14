#!/usr/bin/env python3
"""cas9_bank.py — a reusable Cas9 bank, so the conjunction grid needs no GPU per config.

`all_cut.scan_cas9` re-screens Cas9 guides for every (clean set, band) pair, which is 15-146s a
call and the dominant cost of a config sweep. But what it computes is config-independent: for each
candidate guide, WHICH seeds it fails to cut on. Screen that once over the whole joined window and
every later question is set arithmetic:

    eligible(clean, band) = {g : cutfail(g) INTERSECT clean = {} AND band SUBSET hdr_ok(g)}

which is two numpy reductions over a bool matrix -- microseconds instead of a GPU scan. The same
bank also answers the HDR half, so `hdr_compliance` drops out of the inner loop too.

Stored per guide: the design fields `assemble` needs (`mutation`, `strand`, `start`, `length`,
`gc`, `distance`) plus two bool matrices over the window's seeds -- `cutfail` and `hdr_ok`.

Two things that make this equivalent to `scan_cas9` rather than merely similar:

  * `scan_cas9` walks targets NEAREST-FIRST and stops early once it has `want * pool_target`
    candidates with every (mutation, strand) cell past `cas9_cell_target`. This bank does the full
    walk, so it is a SUPERSET of whatever that scan returned.
  * `assemble` re-sorts each cell by `(distance, |gc - 0.50|)` and keeps only the nearest
    `score_cap // 4` = 625, so a superset that shares the near head produces the same head. That
    is the equivalence `--validate` checks, config by config, against the real `scan_cas9`.

`CB9_MAX_FAIL` bounds what is stored: a guide is usable only if EVERY one of its cut failures lies
outside the clean set, so one with more failures than the window has dirty seeds can never qualify.
60 is far above that bound on all four cells (Cas9 `cut_p` caps at 0.99, so the mean over 300 seeds
is ~3).

!! FALSIFIED 2026-09-13 -- DO NOT USE THIS TO RANK CONFIGS. !!

`cas9_bank_check.py` measured it against the real `scan_cas9` on K562 and the superset property
FAILS on 2 of 3 configs (missing 55 and 1711 guides; a third was a clean superset). The cause is
that the fail-count distribution is Binomial(300, 0.01), **mean 3.0** -- 99.99% of guides have
F < 12 -- so `CB9_MAX_FAIL` at any usable value filters nothing, and the per-cell cap then spends
its whole budget on the NEAREST guides regardless of fail count. Only 4.87% of the 600,000 stored
(29,227) are strict over all 300 seeds. `scan_cas9` instead walks DEEPER into distant targets and
keeps only strict-on-clean guides, so it harvests strict guides from far more targets than a
near-head bank covers.

There is no cheap fix. For a clean set of 232 (dirty 68) a guide is eligible only if every failure
lands in the dirty seeds: F=0 covers 50% of eligible guides, F<=1 85%, F<=2 97%. A true superset
therefore means storing F<=2 from the FULL walk -- ~6.7M guides per cell, ~4 GB resident, ~2h of
rebuild -- and is still not exact. Worse, the error is CONFIG-DEPENDENT (0 to 1711 guides missing
across three configs), so it can reorder arms that differ by less than its own bias, which is
exactly what a config sweep measures. Measured cost of the bias where rows differed: weighted
-0.3% on one arm and +1.7% on another.

What it IS still good for: HEK293, where the clean set is only 77-108 of 300, so ~34% of near-head
guides qualify and a near-head bank is a superset in practice. And any question where the clean set
is FIXED across arms. Neither is the conjunction grid. Use `all_cut.scan_cas9` there; it is 15-146s
a call and correct.

    python cas9_bank.py                     # build for every cell type
    CB9_CELLS=K562 python cas9_bank.py --validate
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")
os.environ.setdefault("NIOME_INSTANCE", "cas9bank")

import dataclasses as _dc            # noqa: E402
import logging                       # noqa: E402
import os.path                       # noqa: E402
import time                          # noqa: E402

logging.basicConfig(level=logging.ERROR)

import numpy as np                   # noqa: E402

import genExp as G                                          # noqa: E402
from niome_subnet.genomics import all_cut as AC             # noqa: E402
from niome_subnet.genomics import mt19937 as MT             # noqa: E402
from niome_subnet.genomics import seed_agnostic as SA       # noqa: E402
from niome_subnet.genomics.validation import stage3         # noqa: E402

BANK_DIR = os.getenv("CB9_DIR", "data/cas9_bank")
MAX_FAIL = int(os.getenv("CB9_MAX_FAIL", "60"))
PER_CELL_CAP = int(os.getenv("CB9_PER_CELL", "150000"))


def bank_path(contract, cfg, seeds):
    key = AC.bank_key(contract, {}, cfg)
    tag = f"{len(seeds)}x{int(seeds[0])}-{int(seeds[-1])}-mf{MAX_FAIL}"
    return os.path.join(BANK_DIR, f"cas9-{key}-{tag}.npz")


def build(contract, reference, cell_types, ctx, sites, cfg, seeds, verbose=True):
    """Full nearest-first Cas9 walk, screened for cut AND HDR over `seeds`."""
    cell = contract.get("cell_type")
    acc = cell_types.get(cell, {}).get("accessibility", 1.0)
    regions = contract.get("mutation_regions") or {}
    sd = np.asarray(seeds, dtype=np.int64)
    jobs = [(i, m, abs(sites[i].start - ctx.mutation_map[m]))
            for i, s in enumerate(sites) if s.cas == "Cas9" for m in ctx.mutations
            if abs(s.start - ctx.mutation_map[m]) <= cfg.max_distance]
    jobs.sort(key=lambda j: j[2])
    recs, t0 = [], time.monotonic()
    for n, (si, mutation, distance) in enumerate(jobs, 1):
        site = sites[si]
        offset = stage3.REGION_ENERGY_OFFSETS.get(regions.get(mutation), 0.0)
        guides = SA.enumerate_variants(site, ctx, cfg.cas9_gc[0], cfg.cas9_gc[1],
                                       ctx.max_mismatches, True, cfg.variants)
        if not guides:
            continue
        params_of = AC._params_fn(site, distance, acc, offset)
        cut = MT.screen_guides_rule_gpu(guides, sd, mutation, "Cas9", site.start,
                                        site.strand, params_of, "cut", MAX_FAIL)
        if not cut:
            continue
        hdr = MT.screen_guides_rule_gpu(list(cut), sd, mutation, "Cas9", site.start,
                                        site.strand, params_of, "hdr", len(sd))
        for guide, fails in cut.items():
            gc, _e, _c = params_of(guide)
            recs.append({"guide": guide, "mutation": mutation, "strand": site.strand,
                         "start": site.start, "length": site.length, "gc": gc,
                         "distance": distance, "cutfail": np.asarray(fails, dtype=np.int64),
                         "hdrfail": np.asarray(hdr.get(guide, sd), dtype=np.int64)})
        if n % 200 == 0:
            MT.free_gpu_memory()
            if verbose:
                print(f"    {n}/{len(jobs)} targets | {len(recs)} guides | "
                      f"{time.monotonic()-t0:.0f}s", flush=True)
    MT.free_gpu_memory()
    # Keep the nearest-first head PER (mutation, strand) cell: `assemble` apportions per cell and
    # a global cap silently empties one, which costs ~0.03x on the whole score.
    by_cell = {}
    for r in recs:
        by_cell.setdefault((r["mutation"], r["strand"]), []).append(r)
    kept = []
    for key, rs in by_cell.items():
        rs.sort(key=lambda r: (r["distance"], abs(r["gc"] - 0.50)))
        kept.extend(rs[:PER_CELL_CAP])
        if verbose:
            print(f"    cell {key}: {len(rs)} -> {min(len(rs), PER_CELL_CAP)}")
    return kept, sd


def save(path, recs, seeds):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    pos = {int(s): j for j, s in enumerate(seeds)}
    n = len(recs)
    cutfail = np.zeros((n, len(seeds)), dtype=bool)
    hdrok = np.zeros((n, len(seeds)), dtype=bool)
    for i, r in enumerate(recs):
        for s in r["cutfail"]:
            cutfail[i, pos[int(s)]] = True
        hdrok[i] = True
        for s in r["hdrfail"]:
            hdrok[i, pos[int(s)]] = False
    np.savez_compressed(
        path, seeds=np.asarray(seeds, dtype=np.int64),
        cutfail=np.packbits(cutfail, axis=1), hdrok=np.packbits(hdrok, axis=1),
        nseeds=np.asarray([len(seeds)]),
        guide=np.array([r["guide"] for r in recs]),
        mutation=np.array([r["mutation"] for r in recs]),
        strand=np.array([r["strand"] for r in recs]),
        start=np.array([r["start"] for r in recs], dtype=np.int64),
        length=np.array([r["length"] for r in recs], dtype=np.int16),
        gc=np.array([r["gc"] for r in recs], dtype=np.float32),
        distance=np.array([r["distance"] for r in recs], dtype=np.int32))


class Cas9Bank:
    """Loaded bank; `select` answers one config with two bool reductions."""

    def __init__(self, path):
        d = np.load(path, allow_pickle=False)
        self.seeds = np.asarray(d["seeds"])
        ns = int(np.asarray(d["nseeds"])[0])
        self.cutfail = np.unpackbits(np.asarray(d["cutfail"]), axis=1)[:, :ns].astype(bool)
        self.hdrok = np.unpackbits(np.asarray(d["hdrok"]), axis=1)[:, :ns].astype(bool)
        self.guide = np.asarray(d["guide"])
        self.mutation = np.asarray(d["mutation"])
        self.strand = np.asarray(d["strand"])
        self.start = np.asarray(d["start"])
        self.length = np.asarray(d["length"])
        self.gc = np.asarray(d["gc"])
        self.distance = np.asarray(d["distance"])
        self.pos = {int(s): j for j, s in enumerate(self.seeds)}

    def __len__(self):
        return self.guide.shape[0]

    def select(self, clean, band=()):
        """Guides strict on cut over `clean` and HDR-compliant on every seed of `band`."""
        ci = np.array([self.pos[int(s)] for s in clean], dtype=np.int64)
        alive = ~self.cutfail[:, ci].any(axis=1)
        if len(band):
            bi = np.array([self.pos[int(s)] for s in band], dtype=np.int64)
            alive &= self.hdrok[:, bi].all(axis=1)
        idx = np.flatnonzero(alive)
        return [{"guide": str(self.guide[i]), "mutation": str(self.mutation[i]),
                 "cas_system": "Cas9", "strand": str(self.strand[i]),
                 "start": int(self.start[i]), "length": int(self.length[i]),
                 "gc": float(self.gc[i]), "distance": int(self.distance[i])} for i in idx]


def get(cell, contract, reference, cell_types, ctx, sites, cfg, seeds, verbose=True):
    path = bank_path(contract, cfg, seeds)
    if not os.path.exists(path):
        t0 = time.monotonic()
        if verbose:
            print(f"  building the Cas9 bank for {cell} (full nearest-first walk)...", flush=True)
        recs, sd = build(contract, reference, cell_types, ctx, sites, cfg, seeds, verbose)
        if not recs:
            raise SystemExit(f"Cas9 bank scan produced nothing for {cell}")
        save(path, recs, sd)
        if verbose:
            print(f"  wrote {path}: {len(recs)} guides in {time.monotonic()-t0:.0f}s "
                  f"({os.path.getsize(path)/1e6:.0f} MB)", flush=True)
    return Cas9Bank(path)
