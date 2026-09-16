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
    def __init__(self, candidates, window_lo=100, window_hi=999, per_cell_min=8, caps=None,
                 seeds=None, prefer=None):
        """``caps`` optionally bounds how many picks a (mutation, cas_system, strand) cell may take.

        The min-union objective is blind to ``mutation_weight``, so on a contract with a heavy and a
        light mutation it lands near a 50/50 split — and every light row costs
        ``weight_heavy - weight_light`` of ``total_weighted_score`` for nothing stage 5 needs beyond
        a non-empty cell. Capping the light cells pushes the remainder onto the heavy mutation.
        Floors still win: a cap below a cell's floor is raised to it, so a cap can never empty a
        stage-5 cell. Passing None keeps the unconstrained selection.
        """
        self.n = len(candidates)
        # ``seeds`` makes the band space an explicit, possibly NON-CONTIGUOUS set (the joined
        # window of several disjoint ranges). Fails are densified through a seed -> position map
        # instead of a `- window_lo` shift; a fail outside the set lands on the sentinel slot and
        # therefore costs nothing, which is correct because such a seed is not part of the band.
        self.seeds = None if seeds is None else np.asarray(sorted(set(int(x) for x in seeds)),
                                                           dtype=np.int64)
        self.span = (window_hi - window_lo + 1) if self.seeds is None else int(self.seeds.size)
        if self.seeds is not None:
            hi = int(self.seeds.max()) + 1
            self._pos = np.full(max(hi, 1000), self.span, dtype=np.int32)
            self._pos[self.seeds] = np.arange(self.span, dtype=np.int32)
        self.per_cell_min = per_cell_min
        self.caps_by_key = dict(caps) if caps else None
        maxf = max(1, max(len(c["fails"]) for c in candidates))
        # Padding points at the sentinel slot `span`, which is held covered forever, so a short
        # fail list contributes exactly its real length and nothing more.
        idx = np.full((self.n, maxf), self.span, dtype=np.int32)
        for i, c in enumerate(candidates):
            f = np.asarray(c["fails"], dtype=np.int64)
            if len(f):
                idx[i, :len(f)] = (self._pos[f] if self.seeds is not None
                                   else f.astype(np.int32) - window_lo)
        self.idx = xp.asarray(idx)
        cells = [(c["mutation"], c["cas_system"], c["strand"]) for c in candidates]
        order = {c: k for k, c in enumerate(sorted(set(cells)))}
        self.cell_id = xp.asarray(np.asarray([order[c] for c in cells], dtype=np.int32))
        self.n_cells = len(order)
        self.available = Counter(order[c] for c in cells)
        if self.seeds is not None:
            self.fail_lists = [
                (self._pos[np.asarray(c["fails"], dtype=np.int64)]
                 if len(c["fails"]) else np.empty(0, dtype=np.int32)) for c in candidates]
        else:
            self.fail_lists = [np.asarray(c["fails"], dtype=np.int32) - window_lo
                               for c in candidates]
        self.caps = ({order[k]: v for k, v in self.caps_by_key.items() if k in order}
                     if self.caps_by_key else None)
        # ``prefer``: one score per candidate, HIGHER is better, used ONLY to break ties among
        # equal-cost picks. The min-union objective is untouched -- a candidate is eligible for the
        # tie-break only when its marginal cost already equals the minimum -- so this cannot trade
        # cut coverage for row quality; it spends slack the greedy was otherwise resolving by index
        # order. CLAUDE.md measures that slack as large: the argmin is tied on 71 of 80 picks,
        # median 12 candidates, max 577.
        #
        # None (the default) keeps ``argmin``'s first-minimum pick, which is what
        # ``assert_matches_sa`` compares against SA. Passing an array changes restart 0 only;
        # restarts 1..n still re-break ties uniformly, so the random arms remain what they were.
        #
        # **MEASURED AND DEAD -- nothing passes this, and it should stay that way.** Built to chase
        # the one structural term with real headroom: on a live HEK293 build ``dist_score`` averaged
        # 0.8497 (17.7% below ceiling) against ``gc_score``'s 0.9505 (5.2%) and ``offtarget_factor``
        # at a perfect 1.0000, and the Cas12a half carried it -- dist 0.7061 there against Cas9's
        # 0.9173, because the min-union selects for cut-failure coincidence and is blind to
        # distance. Two arms over three HEK293 contracts, everything else identical:
        #
        #     arm    dist_score   gc_score   weighted   w x fid   wins
        #     dist     +0.87%      +0.05%     +1.07%    +0.39%    1/3
        #     base     +0.09%      +0.09%     +0.20%    +0.43%    1/3
        #
        # The targeted term moves and the PRODUCT does not, which is the same outcome CLAUDE.md
        # records for GC tie-break, mutation-weight tie-break, band width, group_size,
        # max_distance and light_cell_rows. ``ba815f07`` shows the mechanism cleanly: ``dist``
        # raised weighted 348.11 -> 359.50 (+3.3%) and fidelity fell 0.8787 -> 0.8503 (-3.2%), so
        # ``w x fid`` went 305.9 -> 305.7. Preferring nearer guides shifts the mutation mix, and
        # the mutation coverage entropy term pays back exactly what weighted gains. On the third
        # contract both arms reproduced the baseline byte-for-byte -- no tie ever preferred a
        # different guide.
        #
        # What DID hold is the design guarantee: ``clean`` and ``union`` were identical on 3/3
        # contracts, so the tie-break never traded cut coverage for row quality. The mechanism is
        # sound; the gain is not there.
        self._prefer = None if prefer is None else xp.asarray(np.asarray(prefer, dtype=np.float32))

    def _floors(self, group_size):
        per = min(self.per_cell_min, max(1, group_size // max(1, self.n_cells)))
        return {k: min(per, self.available[k]) for k in self.available}

    def _caps(self, group_size, floor):
        """Effective per-cell ceilings: never below the cell's floor, and never so tight in total
        that ``group_size`` becomes unreachable — an under-filled group would change the
        construction silently rather than decline."""
        if not self.caps:
            return None
        caps = {c: min(self.available[c], max(self.caps.get(c, self.available[c]), floor.get(c, 0)))
                for c in self.available}
        if sum(caps.values()) < group_size:
            return None
        return caps

    def build(self, group_size, rng=None):
        floor = self._floors(group_size)
        caps = self._caps(group_size, floor)
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
            if caps is not None:
                full = [c for c, k in caps.items() if cell_count[c] >= k]
                for c in full:
                    cost = xp.where(self.cell_id == c, BIG, cost)
            if bool((cost >= BIG).all()):
                break
            if rng is None and self._prefer is None:
                pick = int(xp.argmin(cost))               # first minimum, exactly as SA does
            elif rng is None:
                # Deterministic quality tie-break: among the equal-cost minima, take the best
                # ``prefer``. Ties within ``prefer`` itself fall back to the lowest index, so the
                # build stays deterministic.
                m = cost.min()
                ties = xp.flatnonzero(cost == m)
                pick = int(ties[int(xp.argmax(self._prefer[ties]))])
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

    def _prefer_total(self, chosen):
        """Summed ``prefer`` over a chosen set, or 0.0 when no preference is configured."""
        if self._prefer is None or not chosen:
            return 0.0
        return float(sum(float(self._prefer[i]) for i in chosen))

    def best(self, group_size, restarts=12, seed=0):
        """Lowest union wins; ``prefer`` only settles builds that tie on union.

        Union stays the sole objective -- a higher-quality build is never taken over a build that
        covers more seeds. That matters because the clean set is what the whole construction is
        for, and a tie-break that could trade it away would be a different algorithm rather than a
        better-resolved one.
        """
        best_idx = self.build(group_size, rng=None)        # restart 0 is the deterministic build
        best_u = self.union(best_idx)
        best_p = self._prefer_total(best_idx)
        rng = np.random.default_rng(seed)
        for _ in range(max(0, restarts - 1)):
            c = self.build(group_size, rng=rng)
            u = self.union(c)
            if u < best_u or (u == best_u and self._prefer_total(c) > best_p):
                best_idx, best_u, best_p = c, u, self._prefer_total(c)
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
