"""Vectorised min-union greedy, provably the same algorithm as ``SA.min_union_group``.

SA's inner loop recomputes every candidate's marginal cost at every pick, in Python: that is
O(group * n * |fails|) and hits ~814s at n=400k, so a 1.6M-candidate pool is an hour and 3.9M is
two. The cost vector is a gather-and-sum, so it vectorises exactly:

    cost[j] = |{ sd in fails[j] : sd not yet covered }|
            = uncovered[ idx[j] ].sum()      with padding pointed at an always-covered sentinel

Equivalence, not approximation:
  * SA scans j in index order and updates only on ``cost < best_cost``, so it keeps the *first*
    minimum. ``argmin`` also returns the first occurrence -- identical pick, given identical
    eligibility.
  * The eligibility gate (per-cell floors, and the ``free = slots > unmet`` slack reservation) is
    reproduced verbatim.
  * ``assert_matches_sa`` checks the un-jittered build produces the byte-identical chosen set on a
    real bank. That is the only build SA makes deterministically, so it is the only one that *can*
    be compared exactly -- restarts here re-break ties uniformly among argmin ties rather than
    replaying SA's sequential 50% coin, so restart results are comparable in kind, not bit-exact.
"""
from collections import Counter

import numpy as np

try:
    import cupy as _cp
    _cp.arange(1)
    xp = _cp
    GPU = True
except Exception:
    xp = np
    GPU = False


class FastGreedy:
    def __init__(self, candidates, window_lo=100, window_hi=999, per_cell_min=8):
        self.n = len(candidates)
        self.span = window_hi - window_lo + 1
        self.per_cell_min = per_cell_min
        maxf = max(1, max(len(c["fails"]) for c in candidates))
        # Padding points at the sentinel slot `span`, which is held covered forever, so a short
        # fail list contributes exactly its real length and nothing more.
        idx = np.full((self.n, maxf), self.span, dtype=np.int32)
        for i, c in enumerate(candidates):
            f = c["fails"]
            idx[i, :len(f)] = np.asarray(f, dtype=np.int32) - window_lo
        self.idx = xp.asarray(idx)
        cells = [(c["mutation"], c["cas_system"], c["strand"]) for c in candidates]
        order = {c: k for k, c in enumerate(sorted(set(cells)))}
        self.cell_id = xp.asarray(np.asarray([order[c] for c in cells], dtype=np.int32))
        self.n_cells = len(order)
        self.available = Counter(order[c] for c in cells)
        self.fail_lists = [np.asarray(c["fails"], dtype=np.int32) - window_lo for c in candidates]

    def _floors(self, group_size):
        per = min(self.per_cell_min, max(1, group_size // max(1, self.n_cells)))
        return {k: min(per, self.available[k]) for k in self.available}

    def build(self, group_size, rng=None):
        floor = self._floors(group_size)
        uncovered = xp.ones(self.span + 1, dtype=xp.bool_)
        uncovered[self.span] = False                      # the sentinel: padding costs nothing
        chosen = []
        taken = xp.zeros(self.n, dtype=xp.bool_)
        cell_count = Counter()
        BIG = np.int32(1 << 20)

        for _ in range(min(group_size, self.n)):
            slots = group_size - len(chosen)
            unmet = sum(max(0, floor.get(c, 0) - cell_count[c]) for c in floor)
            cost = uncovered[self.idx].sum(axis=1).astype(xp.int32)
            cost = xp.where(taken, BIG, cost)
            if slots <= unmet:
                # No slack left: every remaining pick must land in a cell below its floor.
                need = [c for c in floor if cell_count[c] < floor[c]]
                if not need:
                    break
                allowed = xp.zeros(self.n, dtype=xp.bool_)
                for c in need:
                    allowed |= (self.cell_id == c)
                cost = xp.where(allowed, cost, BIG)
            if bool((cost >= BIG).all()):
                break
            if rng is None:
                pick = int(xp.argmin(cost))               # first minimum, exactly as SA does
            else:
                m = cost.min()
                ties = xp.flatnonzero(cost == m)
                pick = int(ties[rng.integers(len(ties))])
            chosen.append(pick)
            taken[pick] = True
            cell_count[int(self.cell_id[pick])] += 1
            uncovered[xp.asarray(self.fail_lists[pick])] = False
        return chosen

    def union(self, chosen):
        u = set()
        for i in chosen:
            u |= set(self.fail_lists[i].tolist())
        return len(u)

    def best(self, group_size, restarts=12, seed=0):
        best_idx = self.build(group_size, rng=None)        # restart 0 is the deterministic build
        best_u = self.union(best_idx)
        rng = np.random.default_rng(seed)
        for _ in range(max(0, restarts - 1)):
            c = self.build(group_size, rng=rng)
            u = self.union(c)
            if u < best_u:
                best_idx, best_u = c, u
        return best_idx, best_u


def assert_matches_sa(candidates, group_size, cfg, window=900):
    """The un-jittered build must reproduce SA's un-jittered chosen set exactly."""
    import random as _random
    from niome_subnet.genomics import seed_agnostic as SA

    fails = [frozenset(c["fails"]) for c in candidates]
    cell = [(c["mutation"], c["cas_system"], c["strand"]) for c in candidates]
    available = Counter(cell)
    per_cell = min(cfg.per_cell_min, max(1, group_size // max(1, len(available))))
    floor = {k: min(per_cell, available[k]) for k in available}

    # SA.min_union_group's build(jitter=False), lifted verbatim so the reference is the real thing.
    chosen, cover, cell_count = set(), Counter(), Counter()
    while len(chosen) < min(group_size, len(candidates)):
        slots = group_size - len(chosen)
        unmet = sum(max(0, floor.get(c, 0) - cell_count[c]) for c in floor)
        free = slots > unmet
        best, best_cost = None, 10 ** 9
        for j in range(len(candidates)):
            if j in chosen:
                continue
            if not (cell_count[cell[j]] < floor.get(cell[j], 0) or free):
                continue
            cost = sum(1 for sd in fails[j] if cover[sd] == 0)
            if cost < best_cost:
                best_cost, best = cost, j
        if best is None:
            break
        chosen.add(best)
        cell_count[cell[best]] += 1
        for sd in fails[best]:
            cover[sd] += 1

    fg = FastGreedy(candidates, per_cell_min=cfg.per_cell_min)
    mine = fg.build(group_size, rng=None)
    same = sorted(mine) == sorted(chosen)
    return same, fg.union(mine), len({sd for j in chosen for sd in fails[j]})
