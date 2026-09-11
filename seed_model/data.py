"""Windows over the per-cell-type round sequence, and the walk-forward splits.

A sample is one round boundary: the features describe rounds 1..t (the three
classes each of them drew, and the class counts standing before each), and the
target is the class multiset that round t+1 actually drew. That is exactly the
question asked - predict the next task's three seed windows from the previous
ones plus the cumulative counts of the nine windows.

The rounds come out of seeds.json, which scripts/seed_bins.py writes. That file
is the single definition of the nine classes, the cumulative counts and the
strict class ranks, so there is no second implementation here to drift from
what seeds.html shows. Call refresh_seeds() (or --refresh) to regenerate it
before training; that step, not this module, is what talks to the API.
"""

import json
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SEEDS_JSON = ROOT / "seed_model" / "seeds.json"
SUBNET_PYTHON = ROOT / ".venv" / "bin" / "python"     # has the subnet deps

# Class width is a parameter, set once per process by SM_WIDTH. The default 100 reproduces the
# original nine classes (100-199 ... 900-999) exactly, including reading seeds.json's stored
# cumulative_counts/class_ranks/class_last_updated rather than recomputing them -- so the live
# refresh path in round_plan.sh is byte-identical to before this became configurable.
#
# Any other width is REBINNED from each round's raw `seeds`, because those stored fields are
# 100-wide and cannot be reused. 900 must divide by the width: 90/75/60/50/30/20/10 give
# 10/12/15/18/30/45/90 classes.
SEED_LO, SEED_HI = 100, 999
N_SEEDS = SEED_HI - SEED_LO + 1
WIDTH = int(os.environ.get("SM_WIDTH", "100"))
if WIDTH <= 0 or N_SEEDS % WIDTH:
    raise SystemExit(f"SM_WIDTH={WIDTH} does not divide {N_SEEDS} seeds evenly")
N_CLASSES = N_SEEDS // WIDTH
SEEDS_PER_TASK = 3
CONTEXT = 12                   # rounds of history fed to the model
# per-round feature block, in order; this list is the single source of truth
# for FEATURE_DIM, so the model and its checkpoints cannot disagree
FEATURE_BLOCKS = (
    ("drawn", N_CLASSES),      # multiset the round drew, /3
    ("share", N_CLASSES),      # counts before the round, as a distribution
    ("logcount", N_CLASSES),   # log1p(counts), scaled - keeps absolute volume
    ("rank", N_CLASSES),       # strict class rank, mapped to [0, 1]
    ("recency", N_CLASSES),    # rounds since the class was last drawn
    ("meta", 2),               # round index (scaled), gap to the next round
)
FEATURE_DIM = sum(n for _, n in FEATURE_BLOCKS)
CELL_TYPES = ("CD34+_HSPC", "HEK293", "HUDEP-2", "K562")
DRAWN = slice(0, N_CLASSES)
SHARE = slice(N_CLASSES, 2 * N_CLASSES)


def refresh_seeds(start, end, rounds_start=None, rounds_end=None,
                  path=SEEDS_JSON):
    """Re-run seed_bins.py so seeds.json covers the window, then return it."""
    cmd = [str(SUBNET_PYTHON) if SUBNET_PYTHON.exists() else sys.executable,
           str(ROOT / "seed_model" / "seed_bins.py"),
           "--start", start, "--end", end,
           "--rounds-start", rounds_start or start,
           "--rounds-end", rounds_end or end,
           "--json-out", str(path)]
    subprocess.run(cmd, check=True, cwd=ROOT,
                   stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    return load_rounds(path)


def class_name(i):
    """Canonical name of class i at the active width. '100-199' at WIDTH 100."""
    lo = SEED_LO + i * WIDTH
    return f"{lo}-{lo + WIDTH - 1}"


def class_of(seed):
    """Class index of a raw seed, or None if it falls outside 100-999."""
    i = (int(seed) - SEED_LO) // WIDTH
    return i if 0 <= i < N_CLASSES else None


def rebin(by_cell):
    """Recompute the per-round class fields at the active WIDTH from raw `seeds`.

    seeds.json stores seed_classes/cumulative_counts/class_ranks/class_last_updated at width 100
    only, so every other width has to derive them. Semantics are copied from what seed_bins.py
    writes, verified against a stored row: counts and last-drawn are the state BEFORE the round
    (K562 round 26's counts sum to 75 = 25 prior rounds x 3, and its own classes do not appear in
    its class_last_updated), and rank 1 is the highest count.

    The tie-break among equal counts is ours, not seed_bins.py's -- that file's ordering could not
    be reproduced from its output and does not need to be, because a model trained at this width
    only ever sees this function's ranks. It is deterministic (count desc, then most recently drawn,
    then index) which is what matters. Round 1 emits all-zero ranks, matching the stored files, so
    _round_features maps it to the flat middle.
    """
    out = {}
    for cell, rows in by_cell.items():
        rows = sorted(rows, key=lambda r: (r["created_at"], r.get("round", 0)))
        counts = np.zeros(N_CLASSES, dtype=np.int64)
        last = np.zeros(N_CLASSES, dtype=np.int64)
        new = []
        for row in rows:
            drawn = [c for c in (class_of(s) for s in row.get("seeds") or []) if c is not None]
            if counts.sum() == 0:
                ranks = np.zeros(N_CLASSES, dtype=np.int64)
            else:
                order = sorted(range(N_CLASSES),
                               key=lambda i: (-int(counts[i]), -int(last[i]), i))
                ranks = np.zeros(N_CLASSES, dtype=np.int64)
                for pos, i in enumerate(order, start=1):
                    ranks[i] = pos
            new.append({**row,
                        "seed_classes": [class_name(c) for c in drawn],
                        "cumulative_counts": counts.tolist(),
                        "class_ranks": ranks.tolist(),
                        "class_last_updated": last.tolist()})
            for c in drawn:
                counts[c] += 1
                last[c] = int(row.get("round", 0))
        out[cell] = new
    return out


def load_rounds(path=SEEDS_JSON):
    """-> ({cell_type: [round, ...]}, meta) from seeds.json's rounds block."""
    payload = json.loads(Path(path).read_text())
    block = payload.get("rounds")
    if not block:
        raise SystemExit(f"{path} carries no 'rounds' block - regenerate it with "
                         "scripts/seed_bins.py (needs --rounds-start/--rounds-end)")
    by_cell = block["by_cell_type"]
    if WIDTH != 100:
        by_cell = rebin(by_cell)
    return by_cell, {"window": block["window"],
                     "generated_at": payload.get("generated_at"),
                     "width": WIDTH, "n_classes": N_CLASSES}


def _class_index(name):
    """'400-499' -> 3 at WIDTH 100; None for a seed that fell outside the classes."""
    if not name:
        return None
    i = (int(name.split("-")[0]) - SEED_LO) // WIDTH
    return i if 0 <= i < N_CLASSES else None


def _drawn_vector(row):
    v = np.zeros(N_CLASSES, dtype=np.float32)
    for name in row["seed_classes"]:
        b = _class_index(name)
        if b is not None:
            v[b] += 1.0
    return v


def _round_features(row, next_created, n_rounds):
    counts = np.asarray(row["cumulative_counts"], dtype=np.float32)
    total = counts.sum()
    share = (counts / total if total > 0
             else np.full(N_CLASSES, 1.0 / N_CLASSES, dtype=np.float32))
    logcount = np.log1p(counts) / np.log1p(max(total, 1.0))
    ranks = np.asarray(row["class_ranks"], dtype=np.float32)
    # round 1 carries all-zero ranks; treat that as the flat middle
    rank_feat = (np.where(ranks > 0, ranks, (N_CLASSES + 1) / 2.0) - 1.0) / (N_CLASSES - 1.0)
    last = np.asarray(row["class_last_updated"], dtype=np.float32)
    stale = np.where(last > 0, row["round"] - last, float(row["round"]))
    recency = np.minimum(stale.astype(np.float32) / 10.0, 3.0)
    gap = 0.0
    if next_created:
        a = datetime.fromisoformat(row["created_at"])
        b = datetime.fromisoformat(next_created)
        gap = min((b - a).total_seconds() / 86400.0, 3.0)
    meta = np.array([row["round"] / max(n_rounds, 1), gap], dtype=np.float32)
    return np.concatenate([_drawn_vector(row) / SEEDS_PER_TASK, share, logcount,
                           rank_feat, recency, meta])


def _target(row):
    """Normalised class multiset of a round: three draws summing to 1."""
    y = _drawn_vector(row)
    s = y.sum()
    return y / s if s > 0 else np.full(N_CLASSES, 1.0 / N_CLASSES, dtype=np.float32)


def _history(feats, upto, context):
    """Left-padded window of feats[..upto], plus its validity mask."""
    hist = feats[max(0, upto - context + 1): upto + 1]
    block = np.zeros((context, FEATURE_DIM), dtype=np.float32)
    mask = np.zeros(context, dtype=np.float32)
    if hist:
        block[context - len(hist):] = np.stack(hist)
        mask[context - len(hist):] = 1.0
    return block, mask


def _feature_sequence(rows):
    n = len(rows)
    return [_round_features(r, rows[i + 1]["created_at"] if i + 1 < n else None, n)
            for i, r in enumerate(rows)]


def build_samples(rounds, context=CONTEXT):
    """-> dict of arrays, one sample per predictable round boundary.

    x       (n, context, FEATURE_DIM)  history, oldest last, left-padded
    mask    (n, context)               1 where x carries a real round
    y       (n, N_CLASSES)             the next round's class multiset, /3
    cell    (n,)                       cell-type index
    round   (n,)                       the round being predicted
    """
    xs, masks, ys, cells, target_rounds, created = [], [], [], [], [], []
    for cell, rows in sorted(rounds.items()):
        if cell not in CELL_TYPES or len(rows) < 2:
            continue
        feats = _feature_sequence(rows)
        for t in range(len(rows) - 1):
            block, mask = _history(feats, t, context)
            xs.append(block)
            masks.append(mask)
            ys.append(_target(rows[t + 1]))
            cells.append(CELL_TYPES.index(cell))
            target_rounds.append(rows[t + 1]["round"])
            created.append(rows[t + 1]["created_at"])
    if not xs:
        raise SystemExit("no samples - widen the window")
    return {"x": np.stack(xs), "mask": np.stack(masks), "y": np.stack(ys),
            "cell": np.asarray(cells), "round": np.asarray(target_rounds),
            "created": created}


def predict_input(rows, context=CONTEXT):
    """The one unlabelled sample that asks 'what does the next round draw?'."""
    feats = _feature_sequence(rows)
    block, mask = _history(feats, len(rows) - 1, context)
    return block[None], mask[None]


def walk_forward(sample, n_folds=5, min_train=8):
    """Expanding-window splits in time order: never train on a later round.

    Split each cell type's own sequence, then merge, so one fold tests the same
    relative position in every cell type rather than letting the longest series
    dominate the early folds.
    """
    folds = [([], []) for _ in range(n_folds)]
    for c in range(len(CELL_TYPES)):
        idx = np.where(sample["cell"] == c)[0]
        if len(idx) <= min_train:
            continue
        idx = idx[np.argsort(sample["round"][idx], kind="stable")]
        edges = np.linspace(min_train, len(idx), num=n_folds + 1).astype(int)
        for f in range(n_folds):
            lo, hi = edges[f], edges[f + 1]
            if hi <= lo:
                continue
            folds[f][0].extend(idx[:lo].tolist())
            folds[f][1].extend(idx[lo:hi].tolist())
    for tr, te in folds:
        if tr and te:
            yield np.asarray(sorted(tr)), np.asarray(sorted(te))
