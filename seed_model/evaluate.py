"""Metrics and baselines. The baselines are the point of this file.

A model that predicts three classes out of nine gets about one of the three
right by luck alone, so "we hit 1 of 3" is not a result. Everything here is
reported against explicit reference strategies, and the headline comparison
carries a permutation test so an improvement of a tenth of a hit is not read
as skill.
"""

import numpy as np

from .data import DRAWN, N_CLASSES, SEEDS_PER_TASK, SHARE


def top3_from_probs(p):
    """The three classes a distribution would bet on, highest mass first."""
    return np.argsort(-p, kind="stable")[:SEEDS_PER_TASK]


def hits(pred_classes, y):
    """How many of the round's three draws the three predicted classes cover.

    y is the normalised multiset, so y*3 is the count per class. A class drawn
    twice can contribute two hits only if it was predicted - and it can be
    predicted once, so it contributes at most one. That matches how a person
    would score the guess: three named windows against three actual seeds.
    """
    counts = np.rint(y * SEEDS_PER_TASK).astype(int)
    return int(sum(min(counts[c], 1) for c in pred_classes))


def log_loss(p, y, eps=1e-9):
    return float(-(y * np.log(np.clip(p, eps, 1))).sum())


def brier(p, y):
    return float(((p - y) ** 2).sum())


def evaluate_probs(probs, ys):
    """-> dict of aggregate metrics over a set of predictions."""
    h = [hits(top3_from_probs(p), y) for p, y in zip(probs, ys)]
    return {
        "n": len(ys),
        "mean_hits": float(np.mean(h)) if h else 0.0,
        "hit_histogram": {k: int(sum(1 for v in h if v == k)) for k in range(4)},
        "top1_accuracy": float(np.mean([
            np.rint(y * SEEDS_PER_TASK)[np.argmax(p)] > 0 for p, y in zip(probs, ys)
        ])) if len(ys) else 0.0,
        "log_loss": float(np.mean([log_loss(p, y) for p, y in zip(probs, ys)])),
        "brier": float(np.mean([brier(p, y) for p, y in zip(probs, ys)])),
        "per_sample_hits": h,
    }


# ---- reference strategies ------------------------------------------------
# Each takes the training targets seen so far plus the test sample's history
# and returns a distribution over the nine classes.

def baseline_uniform(train_y, x, mask):
    return np.full(N_CLASSES, 1.0 / N_CLASSES)


def baseline_marginal(train_y, x, mask):
    """The empirical class frequency over everything drawn so far."""
    if len(train_y) == 0:
        return baseline_uniform(train_y, x, mask)
    p = train_y.mean(0)
    return p / p.sum()


def baseline_hot(train_y, x, mask):
    """Bet on the classes with the highest cumulative count - 'follow the
    leader'. Reads the count share straight out of the last history row."""
    share = _last_share(x, mask)
    return share if share is not None else baseline_uniform(train_y, x, mask)


def baseline_cold(train_y, x, mask):
    """The gambler's-fallacy strategy: bet on the least-drawn classes."""
    share = _last_share(x, mask)
    if share is None:
        return baseline_uniform(train_y, x, mask)
    inv = share.max() - share + 1e-3
    return inv / inv.sum()


def baseline_repeat(train_y, x, mask):
    """Bet that the next round repeats the previous round's three classes."""
    idx = np.where(mask > 0)[0]
    if len(idx) == 0:
        return baseline_uniform(train_y, x, mask)
    drawn = x[idx[-1], DRAWN]                            # "drawn" block, /3
    if drawn.sum() <= 0:
        return baseline_uniform(train_y, x, mask)
    p = drawn + 1e-3
    return p / p.sum()


def _last_share(x, mask):
    idx = np.where(mask > 0)[0]
    if len(idx) == 0:
        return None
    share = x[idx[-1], SHARE].astype(np.float64)         # "share" block
    s = share.sum()
    return share / s if s > 0 else None


BASELINES = {
    "uniform": baseline_uniform,
    "marginal": baseline_marginal,
    "hot_hand": baseline_hot,
    "cold_hand": baseline_cold,
    "repeat_last": baseline_repeat,
}


def run_baselines(sample, train_idx, test_idx):
    out = {}
    train_y = sample["y"][train_idx]
    for name, fn in BASELINES.items():
        probs = [fn(train_y, sample["x"][i], sample["mask"][i]) for i in test_idx]
        out[name] = evaluate_probs(np.asarray(probs), sample["y"][test_idx])
    return out


def permutation_test(model_hits, ref_hits, n_iter=20000, seed=0):
    """Paired sign-flip test on the per-sample hit difference.

    H0: the model and the reference are exchangeable on every sample. Returns
    the one-sided p for 'model better'. With a few dozen samples this is the
    honest way to talk about a small edge.
    """
    d = np.asarray(model_hits, dtype=np.float64) - np.asarray(ref_hits, dtype=np.float64)
    if len(d) == 0:
        return {"observed": 0.0, "p_value": 1.0, "n": 0}
    obs = d.mean()
    rng = np.random.default_rng(seed)
    signs = rng.choice([-1.0, 1.0], size=(n_iter, len(d)))
    null = (signs * d).mean(axis=1)
    p = float(((null >= obs).sum() + 1) / (n_iter + 1))
    return {"observed": float(obs), "p_value": p, "n": int(len(d))}
