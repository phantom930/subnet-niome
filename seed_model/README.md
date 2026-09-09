# seed_model — next-3-seed-window prediction, per cell type

Predicts the three seed classes (100–199 … 900–999) of a cell type's **next**
3-seed task from that cell type's previous rounds: the classes each earlier
round drew, and the cumulative counts of all nine classes standing before it.

One model per cell type, trained and finetuned separately, refreshed
continuously from the live task feed.

## The result, first

**The model does not beat chance, and neither does anything else.** Over 79
held-out rounds in walk-forward backtest:

| strategy | hits / 3 |
|---|---|
| cold_hand (bet the least-drawn classes) | 1.000 |
| repeat_last (bet the previous round's classes) | 0.975 |
| uniform | 0.949 |
| hot_hand (bet the most-drawn classes) | 0.886 |
| **SeedFormer** | **0.886** |
| marginal (empirical class frequency) | 0.848 |
| chance (any 3 of 9) | 1.000 |

`model − cold_hand = −0.114 hits/round, p = 0.869` (paired sign-flip,
one-sided). Per cell type the p-values are 0.88–1.00. Every strategy sits on
top of chance at 1.000, which is what three guesses out of nine classes gets
you for free.

The model's *log-loss* is 2.228 against uniform's 2.197 (= ln 9) — it learned
to output almost exactly the uniform distribution, deviating by at most 0.027
in any class. That is the correct answer to an unpredictable target, and it is
consistent with the rest of this repo's findings: `seed_pred/` found no key
recipe over ~395M candidates, and seeds.html's strict class ranks average
4.82–5.26 against an even-draw value of exactly 5.0.

Treat `next_prediction.json` as three arbitrary class names with a ~1-in-3 hit
rate each. The file says so in its own `calibration.reading` field.

**Why it is still worth having:** it is the negative control. If the generator
ever changes — a seeded RNG, a drifting distribution, a bug that correlates
consecutive tasks — the walk-forward number moves off chance and the promotion
guard in `update.py` records it. Right now it is measuring that nothing is
there.

## Data

22–32 usable rounds per cell type (115 total, ~111 predictable boundaries).
**This is the binding constraint** — no architecture recovers signal from 30
samples that the data does not contain, and the numbers above are what an
honest protocol on this much data looks like.

Rounds are read from `seeds.json`'s `rounds.by_cell_type` block, written by
`scripts/seed_bins.py`. That script is the single definition of the classes,
the cumulative counts and the strict ranks, so the model and `seeds.html`
cannot drift apart. `--refresh` re-runs it before training.

Per-round features (65 dims, `data.FEATURE_BLOCKS`): the round's own class
multiset, the prior counts as a distribution, log1p counts (absolute volume),
the strict class rank scaled to [0,1], rounds since each class was last drawn,
and two meta values (round index, gap to the next task).

Target: the next round's three draws normalised to sum to 1, so a class drawn
twice carries twice the mass. Loss is soft cross-entropy with label smoothing.

## Model

`SeedFormer` — 3-block pre-norm Transformer, dim 64, 4 heads, 152k params:

- **RMSNorm** pre-norm blocks, no bias terms — stable at batch sizes this small
- **SwiGLU** feed-forward
- **Rotary position embeddings** — relative positions, no learned table to overfit
- **Masked attention** over the left-padded history (context 12 rounds)
- **Attention-pooled readout** rather than last-token, so a short history is
  not dominated by whichever round sits at the end
- **Cell-type embedding**, so one model pretrains across all four sequences
- **Stochastic depth + dropout 0.2 + weight decay 0.05**, and an **EMA shadow**
  used for every evaluation
- AdamW, cosine schedule with warmup, gradient clipping, early stopping on the
  training split's own time-ordered tail, 3 independent inits averaged

**Pretrain pooled → finetune per cell type.** Fitting a transformer on 22
samples memorises them. Pooling the four sequences gives ~111 samples to learn
the shape of the problem from; the per-cell-type finetune is a short
low-learning-rate pass on that cell type alone, with its own checkpoint. The
deliverable is still one separately-finetuned model per cell type.

## Protocol

Walk-forward only: a fold trains on rounds 1…k of each cell type and tests on
the rounds after k, so the model never sees a task earlier than one it
predicts. Splits are computed per cell type then merged, so a fold tests the
same relative position in every sequence. Five reference strategies are scored
on the identical folds, and the headline comparison carries a paired
sign-flip permutation test — with a few dozen samples that is the only honest
way to talk about a tenth of a hit.

## Usage

```bash
# one full cycle: refresh data, backtest, finetune per cell type, predict
.venv-ml/bin/python -m seed_model.update --end now

# keep cycling hourly
.venv-ml/bin/python -m seed_model.update --end now --loop 3600

# pieces
.venv-ml/bin/python -m seed_model.train --refresh --end now   # -> report.json
.venv-ml/bin/python -m seed_model.predict                     # -> next_prediction.json
```

`update.py` warm-starts from the promoted checkpoints, and **rejects** a run
whose walk-forward score regresses by more than `--tolerance` (default 0.10
hits/round), restoring `checkpoints_prev/`. Every cycle appends to
`state.json`, so the score's history over time is the thing to watch.

Torch lives in `.venv-ml` (CPU build), kept separate from the subnet's `.venv`
so the miner runtime is untouched. Full run ≈ 7 min on 16 CPU threads.

## Files

| path | role |
|---|---|
| `data.py` | rounds → samples, features, walk-forward splits |
| `model.py` | SeedFormer, EMA, soft cross-entropy |
| `train.py` | pooled pretrain + per-cell finetune + backtest → `report.json` |
| `evaluate.py` | metrics, the five baselines, permutation test |
| `predict.py` | next-round distribution → `next_prediction.json` |
| `update.py` | continuous cycle, promotion guard → `state.json` |
| `checkpoints/` | `pooled.pt` + one `.pt` per cell type |
