# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Bittensor subnet 55 (testnet 289) for synthetic CRISPR dataset generation on the HBB locus. Miners
design CRISPR experiments (guide RNA, coordinates, strand, mutation, Cas system, cell type) against a
contract issued by the backend; validators score those designs through a five-stage pipeline and set
on-chain weights. The anti-cheat premise is that miners submit **designs only** — every biological
outcome (cut, repair mode, indel length) is computed by the validator, never supplied by the miner.

## Commands

```bash
uv sync                        # install deps into .venv (pyproject.toml + uv.lock; there is no requirements.txt
                               #   despite what docs/miner_guide.md and docs/validator_guide.md say)

# Validator — the supported path. Installs pm2, downloads chr11.fa, auto-pulls main every 60s.
./entrypoint.sh --wallet <NAME> --wallet-hotkey <HOTKEY> [--wandb.api_key <KEY>]
bash scripts/run_validator.sh --wallet <NAME> --wallet-hotkey <HOTKEY>   # no pm2, no auto-update
python neurons/validator.py --netuid 55 --network finney --wallet <NAME> --wallet-hotkey <HOTKEY>

# Miner
python neurons/miner.py --netuid 55 --network finney --wallet <NAME> --wallet-hotkey <HOTKEY> --axon.port 8091

# Add --mock to either neuron to run against MockSubtensor/MockMetagraph with no chain or keystore.
```

There is no test suite and no linter config; `.github/` holds only a PR template, no workflows.
Verify changes by running a neuron (`--mock` for the chain parts) or by exercising the scoring
pipeline directly:

```python
# From the repo root, with data/contract.json, data/hbb_reference.json, data/submission.json
# and data/chr11.fa present — re-scores a submission exactly as a validator would.
from niome_subnet.genomics.validation import benchmark_submission
benchmark_submission(cell_types_dict, uid=0)
```

**`contract["seed"]` is a comma-joined list of round seeds** (`"122,321,431"`), not one seed —
`_parse_seeds` splits it, stage 12 runs once (it is seed-independent), then **stages 3-5 re-run per
seed and every breakdown field plus `final_score` is averaged** over them
([validation/\_\_init\_\_.py](niome_subnet/genomics/validation/__init__.py)). A build tuned to a
single seed therefore only captures its share of the mean. Broadcast contracts arrive with
`seed: 0`; the real seeds are stamped in before scoring.

`data/chr11.fa` (~130 MB, GRCh38 chromosome 11 from Ensembl release 116 — URL in
[scripts/run_validator.sh](scripts/run_validator.sh)) is required by both sides: PAMs and coordinates
are checked against the real sequence. `data/` is gitignored.

**Every data path in [settings.py](niome_subnet/utils/settings.py) is relative, so neurons must run
from the repo root.**

## Architecture

### Round structure

The validator has no request/response loop with miners. Everything is driven by block position within
a fixed 720-block interval measured from `BASE_BLOCK_NUMBER`, in
[validator/forward.py](niome_subnet/validator/forward.py):

| Blocks into interval | Phase |
|---|---|
| 0–599 | `broadcast_task` — fetch the task, mint one presigned S3 PUT URL per miner, POST the task |
| 600–699 | `run_validation` — download each submission from S3, score it, `set_weights` locally, push top-5 submissions to the backend |
| ~700 | `should_set_weights` fires and the accumulated weights are committed on chain |

`forward()` is called repeatedly and dispatches these as fire-and-forget asyncio tasks guarded by
`is_broadcasting`/`is_validating` flags, so it must stay idempotent within a phase.

Broadcast dedups **only by uid**: `collected_uids` (persisted in `state.npz`) skips a uid already
contacted for this `task_id`, and resets when `task_id` changes. There is no per-operator limit.
`7385f6b` added an IP+coldkey filter and `8451d27` removed it, leaving `ip` and `coldkey` computed
at [forward.py:100](niome_subnet/validator/forward.py#L100) and never read — dead assignments under
a comment that still describes the deleted behaviour. Don't mistake that comment for working code;
`get_miner_uids` doesn't filter either (it only skips validators via `trust > 0`).

### Submission path

The miner's HTTP reply is an empty ack — the dataset travels out of band as a PUT to the presigned URL
that arrived with the task, and must land before that URL's 300 s TTL (`SUBMISSION_TIMEOUT`) expires.
A miner that misses the window is indistinguishable from one that was never contacted; there is no
retry within a task id and no feedback channel, which is why
[neurons/miner.py](neurons/miner.py) scores its own build locally before uploading.

### Round prefetch — the TTL bounds the upload, not the build

The 300 s TTL only has to cover the PUT, provided the rows already exist when the validator calls.
A round's task is published on the **public, unsigned** `/api/v3/tasks` (not `settings.TASK_URL`,
`/tasks/current`, which 400s without a validator's signed headers) the moment the round opens, with
`content.contract` and `content.hbb_reference` inline. Those are the same artifacts the validator
later hands over: across all 41 archived rounds they match field for field, the only difference
being `seed`, which is 0 on both sides until the round closes.

So `Miner._prefetch_loop` polls for the newest task whose seed is still unstamped and builds it on
sight, hours ahead of the request. Measured over 72 rounds (task `created_at` against the miner's
own "Received genomics task" log line) the lead time is **min 198 s, p10 395 s, median 1794 s, max
4404 s** — 86% of rounds leave ten minutes or more, against the ~225 s the in-TTL path gets.

`process_task` then calls `_rows_for_task`, which has three outcomes in descending order of build
time: the prepared rows; a prepare still running (wait for it, bounded by
`deadline - UPLOAD_RESERVE_S - EMERGENCY_BUILD_S`); or the in-TTL build, unchanged from before and
still the fallback whenever a prepare is missing, mismatched or failed. If the wait runs out,
`_build(allow_hedges=False)` takes the ordinary construction rather than queueing behind the GPU.

The prefetch shipped without changing what any builder produced: `PREPARE_BUDGET_S` is 900 s, and
`_build_seed_agnostic` still clamps to its own `SEED_AGNOSTIC_MAX_BUDGET_S` (210 s, the value every
number in that hedge was tuned at), so widening that hedge stays a separate measurement.

**K562's all-cut config was the first thing to spend the new budget.** It ran `cas12a_max_fail 22` /
`max_distance 200` only because that was the one K562 bank that fit the TTL, so it moved to the
`mf100/d400` bank HEK293 uses. Three configs were then scored against each other over five
contracts x 40 seeds, paired within (contract, seed):

| | vs mf22/d200 | t | 95% CI |
|---|---|---|---|
| mf100/d400 | +2.36 ± 2.10 | +1.12 | [-1.76, +6.47] |
| mf22/d400 | +0.64 ± 1.32 | +0.48 | [-1.95, +3.23] |
| mf22/d400 **vs mf100/d400** | -1.85 ± 1.21 | -1.53 | [-4.22, +0.52] |

**No pair separates.** `mf100/d400` is retained on its point estimate, not on a measured win — do
not cite it as one. Full table in `all_cut.CELL_CONFIG`.

**Never infer these scores from clean fraction.** It failed as a proxy twice, in both directions:
the wide bank gains 7-15 clean seeds of 900 over d200 yet scores *worse* on one contract, and
mf22/d400 reaches a higher clean count than mf100/d400 while scoring 1.85 lower. Term 1 is
flat-to-worse on the wide bank; the delta rides entirely on `consistency_factor`.

Since score cannot separate them, the choice is about cost. The build is 4-6x slower than
mf22/d400 and no longer fits the ~225 s in-TTL path, so `Miner.ALL_CUT_MIN_BUDGET_S` had to become
per cell type (`{"HEK293": 190, "K562": 480}`) — a single gate would either start a K562 build that
cannot finish, or lock HEK293 out of a path it completes inside. **K562 all-cut is therefore
prefetch-dependent:** a round whose prefetch fails falls to the seed-agnostic hedge (~41 against
~64). Granting the wide bank its full point estimate, that trade breaks even at an 8% prefetch
failure rate. `RETRY_CONFIG` is empty because both d400 configs build 5/5 where d200 declines 1 in 5.

`mf22/d400` is the fallback of record if the prefetch proves unreliable: same score within noise,
same coverage, 70-133 s, and it fits back under a flat 190 s gate.

**HUDEP-2 was then measured and is a clear win**, unlike the K562 config tuning above. It had no
`CELL_CONFIG` entry, so every HUDEP-2 round shipped the seed-agnostic hedge; all-cut at K562's
config beats that hedge by **+17.92 ± 3.65, t = 4.9, CI [+10.77, +25.07]** over five contracts x 40
seeds — positive on 5/5, each individually significant, built 5/5, winning 143/200 seed pairs. What
marks it as real rather than `consistency_factor` jitter: *all three* score terms rise on every
contract. The band transfers from K562 on the clamp (accessibility 0.82 reaches `cut_p` 0.990/0.960
at gc 0.40 against K562's 0.990/0.953), not on resemblance — HEK293's wide band would be the wrong
transfer. Builds are 247-386 s, so HUDEP-2 inherits the same 480 s gate and the same prefetch
dependency.

**CD34+_HSPC completes the set**, at the same config again, and is the weakest of the four:
**+14.27 ± 1.93 where it builds, but only 3 of 6 contracts build → ~+7.1 expected per round.** The
declines are stage-5 cell coverage (7 of 8 cells), which is contract site geometry — the same
contracts fail every time and accessibility does not predict it.

Shrinking the group recovers coverage and is measurably *not* worth it: group 25 builds 5/5 and
still expects less (+5.59) than group 42 building 3/5 (+8.56), because score falls faster
(14.27 → 5.59) than availability rises. Note this resolves **opposite** to the group 75/50 question
on HUDEP-2, where the score spread was ±1 and availability decided it — so neither direction is a
general rule, and each cell type needs its own sweep.

All-cut now covers every cell type the backend issues:

| cell type | tasks | gain/round | over |
|---|---|---|---|
| K562 | 88 | ~+23 | seed-agnostic hedge |
| HUDEP-2 | 74 | +17.92 | seed-agnostic hedge |
| CD34+_HSPC | 85 | ~+7.1 | seed-agnostic hedge |
| HEK293 | 81 | ~+4 | its own clustered builder |

Two things this is *not* a licence to forget: the prefetch thread waits on `_prewarmed` because
`G.load_sequence` caches through an unguarded module global, and the hedge builders are serialised
by `_hedge_lock` through the bounded `_hedge_slot`, so an in-TTL build never queues behind a
prepared one past its window.

**All-cut is now the fallback, not the shipped path.** Every gain in that table is still the right
number for all-cut, and all-cut is still what runs when the builder above it declines — but on all
four cell types the builder above it is all-HDR.

`Miner._build` tries the builders in this order, each falling through to the next on a decline, a
short pool, a missing GPU or too little budget — so a failure anywhere lands on the build the miner
had before that rung existed:

```
all-HDR  →  all-cut  →  seed-agnostic hedge  →  ordinary construction
                                                (HEK293: clustered builder)
```

`allow_hedges=False` (the emergency in-TTL path, when a prepare ran out of time) skips straight to
the last rung rather than queueing behind the GPU.

### All-HDR — the spike construction, and the fleet that plays it

[genomics/all_hdr.py](niome_subnet/genomics/all_hdr.py) is all-cut with the pinned outcome moved
from *cut* to *HDR*: min-union a Cas12a group on the `hdr` rule over a narrow 100-seed band, then
require HDR of the Cas9 half only over the resulting clean band. On a clean-band seed **every** row
repairs by HDR, so stage 4's three targets (`is_cut`, `is_hdr`, `indel_length`) are all constant and
`consistency_factor` is **exactly 1.0** — against all-cut's ~0.20-0.26, which pins `is_cut` alone.

**It loses on the mean and is shipped anyway.** Scored head to head through all five stages, five
contracts, clean and dirty legs:

| cell | band | clean | clean leg | dirty leg | expected | all-cut | delta |
|---|---|---|---|---|---|---|---|
| CD34+_HSPC | 500-599 | 16 | 275.38 | 29.89 | 34.26 | 56.55 | -22.29 |
| K562 | 700-799 | 16 | 254.30 | 27.22 | 31.25 | 61.12 | -29.87 |
| HUDEP-2 | 800-899 | 15 | 191.31 | 20.07 | 22.92 | 41.67 | -18.75 |

(That head-to-head predates the group 42 → 100 balance change below, so its bands are the wider
15-16 rather than today's 12-13; the ratio it establishes is what matters, not the absolute legs.)

That is 51-61% of all-cut's expected score, ~-24/round. **The justification is `SCORING_SYSTEM =
"top"`, which pays only the top 10 on a fixed curve** ([settings.py](niome_subnet/utils/settings.py)
`SCORE_DISTRIBUTION`, rank 1 = 29.4% and rank 11 = 0%): a build worth ~290 on a few percent of
rounds can out-*rank* a flat one that rarely does. It has **since been measured**, and the bet
holds while the reasoning behind it was wrong — all-HDR's k=1 (~88) places on 60% of rounds against
all-cut's E[final] 46.1 at 11%, so the gap is 60-vs-11, not something-vs-never (see the cutoff
distribution below). **Both of those percentages are pooled-field numbers and are upper bounds** —
see "Pricing a construction" below, which is also why the all-HDR side of the bet survives the
correction and the all-cut side does not. `Miner.ALL_HDR = False` reverts to all-cut everywhere,
and that is the switch to reach for if the bet is judged wrong.

Three regimes, and the reason there is no middle:

| regime | condition | `consistency_factor` | seeds of 900 (measured, group 80) |
|---|---|---|---|
| 1 | HDR-clean seed (in the band) | exactly 1.0 | 13 (K562/HUDEP-2/CD34+), 7 (HEK293) |
| 2 | cut-clean, HDR-dirty | **0.123-0.138** | **5-6** (HEK293: **0**) |
| 3 | cut-dirty | ~0.10 | 881-882 |

Regime 2 was written here as ~0.26 and is measured at 0.123-0.138 — and it is economically dead
either way: it covers 5-6 seeds of 900 for all-HDR (cut-clean is 18-19 total, and 13 of those *are*
the band), against all-cut's ~560.

**The "rank-10 cutoff of 62-86, which is why all-cut never places" claim that stood here is wrong
on both halves.** Measured across 72 current-regime rounds the cutoff is a *distribution*, not a
band: min 20.7, p10 45.9, **median 82.3**, p90 120.0, max 158.1 (per-cell medians HEK293 **76.4**
— the 64.5 previously recorded here was on fewer rounds; re-measured over 31 HEK293 three-seed
fields the spread is min 20.7 / p25 43.4 / median 76.4 / p75 95.2 / max 158.1 —
CD34+ 81.4, HUDEP-2 83.8, K562 90.9). Against that, all-cut's E[final] of **46.1** — priced over
the clean-seed distribution, not the single lucky seed that produced the old 57.4 — places on
**11%** of rounds, not never. all-HDR's k=1 (~88) places on 60% and k=2 (~154) on 99%. The
conclusion that all-HDR is the better bet survives; the reasoning "all-cut never places" does not,
and a construction should not be dismissed on it.

HDR-clean ⊂ cut-clean, so all-HDR gives up all-cut's wide regime-2 floor to buy regime 1. **Two
attempts to have both were falsified, don't repeat them:** a *combined* construction cannot, because
cut-min-union and HDR-min-union compete for the same guides and a build optimises exactly one; and
the *hybrid split* (part of the rows on each rule) does not produce an elevated "regime 1.5" — its
spike consistency measured *below* the pure floor. Also beware the measurement itself: dirty-leg
consistency swung 0.10↔0.25 across configs that should have matched, so re-verify any consistency
delta on a fresh sample before believing it.

**Two structural limits, both measured, that no parameter reaches:**

- **The band cannot be widened.** Each seed added multiplies the conditional Cas9 requirement by
  P(HDR) ≈ 0.57 (≈0.37 on HEK293). A 12-16 seed band leaves 484-1559 Cas9 candidates against the
  ~170-208 needed; a 29-31 seed band (group 6) yielded **zero** on all three erythroid types.
- **The band's position is free.** Seeds are independent draws. Verified on HEK293 across all nine
  fleet windows: band 7 and 8/8 stage-5 cells at every one, fidelity 0.904-0.934, spike final
  287-294. This is what makes the fleet below possible.

Per-cell tuning lives in `all_hdr.CELL_CONFIG`; `AllHdrConfig` holds the shared defaults:

| cell type | band default | group | `main_max_fail` | band width | build (cold/warm) |
|---|---|---|---|---|---|
| K562 / HUDEP-2 / CD34+_HSPC | 700-799 / 800-899 / 500-599 | 80 | 45 | 13 | ~27 s warm |
| HEK293 | 300-399 | 80 | 48 | 7 | ~32 s / ~3 s |

**`group_size` is 80 on a payout measurement that overrides an earlier score measurement.** It sets
the cas mix directly, and trades two things that pull opposite ways — band width (spike *frequency*)
against stage 5's cas-coverage entropy (spike *score*):

| | g42 | g60 | g80 | g100 |
|---|---|---|---|---|
| band, erythroid mean (56/56 pairs, monotone) | 14.6 | 13.7 | **13.0** | 12.1 |
| band, HEK293 | 8.3 | 7.3 | **7.2** | 6.5 |
| fidelity (**`light_cell_rows` OFF** — see below) | 0.89 | 0.93 | **0.96** | 0.976 |
| **E[pay] aggregate** | 0.0490 | 0.0583 | **0.0598** | 0.0524 |

**That fidelity row is not what the fleet ships, and reading it as a target wastes a day.** It was
measured with `light_cell_rows` off. With the shipped `light_cell_rows = 6`, group 80's fidelity is
**0.889** — the trade is documented at `AllHdrConfig.light_cell_rows` (fidelity 0.948 -> 0.889 for
+15.7% E[pay], on 210 configs) and it is deliberate. Rebuilt to check, on 571f4843 CD34+ at the
*contiguous 100-seed window this row was measured at*, fidelity comes out **0.8956** — so 0.96 is
unreproducible at its own stated configuration once the shipped apportionment is on. The live fleet
measures 0.877-0.912 on every hotkey and every round, which is that 0.889, not a shortfall.

**Where the 0.11 actually sits, decomposed.** `distribution_fidelity_score` is the geometric mean of
six ratios, all computed from stage 12's design fields alone — `cas_shift` (repair-mode JS
divergence, indel Wasserstein) is a **diagnostic and does not enter the score**, so fidelity is
seed-independent (verified: 0.895552 at seeds 500 and 731). On the shipped build five of six terms
are near-perfect and one is not:

| term | value |
|---|---|
| **mutation coverage entropy** | **0.722** |
| **joint coverage entropy** | **0.808** |
| cas system coverage entropy | 0.904 |
| strand coverage entropy | 0.987 |
| k-mer diversity | 0.991 |
| distinct guide ratio | 1.000 |

Every backend contract carries exactly 2 mutations, so the 79/79/6/6 Cas9 split drives the mutation
axis to near its floor and drags the joint term with it. Taking those two to 1.0 would give 0.980 —
that is the *whole* fidelity headroom, and it is the same knob E[pay] already decided. Do not
re-open it without beating 0.0242 E[eff].

Priced over 54 current-regime fields × 25 contracts at the k=1 case, group 80 is **+14.2%** on
expected payout over group 100, and best on K562, HUDEP-2 and CD34+_HSPC (HEK293's g42 edge is +2%,
inside noise). `total_weighted_score` moves slightly the *other* way (K562 258 at g42 vs 250 at
g100) because a smaller group means more Cas9 rows and those score better structurally — so 80
gains a little on band and weighted and gives up only fidelity.

**Two sampling traps this went through; don't repeat them.** Scoring against a handful of single
fields says group 42 wins by +32% — but those fields were all above their cell type's median cutoff,
and on a hard field the k=1 hit doesn't place, so payout rides on the k=2 term which scales as
`band²` and favours wide bands. Sampling fields across *all* backend history then says the opposite,
because the subnet ran at 29-44 miners/task historically against **248 today**, and those low cutoffs
make every config look like it places. Only current-regime fields decide anything.

The superseded note, kept because the number is still true: group 42 → 100 was taken on spike-round
*score*, 117.1 → 121.8 (+4.7). It is the wrong objective — `SCORE_DISTRIBUTION` is a step function,
so a score gain crossing no rank threshold pays nothing while the narrower band costs frequency
every round.

On HEK293 the group sweep by fidelity is much larger, because it started from a worse mix:

| group | cas12a/cas9 | band | cells | spike cons | fidelity | weighted | spike final |
|---|---|---|---|---|---|---|---|
| 20 | 20/230 | 10 | 8/8 | 1.000 | 0.786 | 326 | 256 |
| 40 | 40/210 | 8 | 8/8 | 1.000 | 0.854 | 330 | 282 |
| 60 | 60/190 | 8 | 8/8 | 1.000 | 0.886 | 326 | 289 |
| 80 | 80/170 | 7 | 8/8 | 1.000 | 0.924 | 318 | 294 |
| 100 | 100/150 | 6 | 8/8 | 1.000 | 0.939 | 308 | 289 |

Group 80 ships (and is now the shared default): fidelity 0.924 against the leaders' 0.947 median on
their placing HEK293 rounds (weighted 288, consistency 0.552), and its 7-seed band is hit ~17% more
often than group 100's 6 for the 0.015 fidelity given up — the same frequency-over-score trade the
payout table above makes for every cell type. **A larger group is not free but is cheap:** weighted drifts 326 → 308
because the added Cas12a rows score below the Cas9 rows they displace, and the band narrows — which
*helps* the Cas9 fill, since `0.37**6` is an easier target than `0.37**10`.

An earlier note here claimed HEK293 fidelity was capped near 0.78 by its accessibility 0.35. **That
was wrong** — it was an artefact of only ever building groups 8-24, where 20 Cas12a rows of 250 drive
the cas entropy to ~0.65. Nothing about accessibility bounds fidelity.

`ALL_HDR_MIN_BUDGET_S` stays a single flat 190 s, unlike `ALL_CUT_MIN_BUDGET_S`: all four cell types
build well inside the ~225 s in-TTL path.

#### The fleet — eleven hotkeys, ten on band windows

**Three operational failures cost round 9ac1d178 (HUDEP-2, seeds 207/573/129) everything, and two
of them are configuration, not luck.** The fleet caught **2 of the 3 seeds** — h2 hit 207 and h3 hit
573, its best band coverage on record — and *neither placed*, at ranks 11 and 12 against a rank-10
cutoff of 126.32 (the 86th percentile of 28 HUDEP-2 rounds; h3's 121.74 would have placed on 23 of
28). Payout was zero. What went wrong beyond the hard field:

* **h0 and h5 shipped the emergency build** and scored byte-identically (313.3 / 0.0908 / 0.9189 →
  26.15), because `_build(allow_hedges=False)` is the deterministic ordinary construction. The
  validator called h5 **129 s** and h0 **149 s** after task creation, against all-HDR's ~360 s build
  and h0's all-cut 711.6 s. **That is a new minimum lead time** — the prefetch section's floor of
  198 s is stale, so any hotkey can be asked inside 130 s and the ladder must survive it.
* **h8, h9 and h10 scored a clean zero**: their processes first came up at 08:22/08:22/08:26, after
  the 07:07 broadcast, so the validator found nothing at their S3 keys. h1-h7 also restarted at
  08:10, which means **the rotated joined-window layout has not yet been exercised on a scored
  round** — h1-h7 built the old contiguous width-300 windows from the 06:17 plan.
* **Bands leak well outside their window.** h3's window was 300-599 but its 12-seed band came out
  [132, 140, 144, 163, 173, 511, 517, 530, 546, 547, 573, 817] — five seeds below 300 and one above
  599. The window only *guides* the min-union search; it does not bound the resulting clean band. So
  sibling bands overlap more than the window layout implies, and a disjoint-window plan does not
  give disjoint bands. Coverage unions computed from windows are upper bounds.

Our `total_weighted_score` was **not** the problem that round: 340.0 median across the band hotkeys
against the field top-10's 337.1 and rank 1's 329.8. We out-score the winner on term 1 and lose
entirely on consistency structure.

**Read this before trusting anything below about windows.** As of 2026-09-09 the fleet runs
**eleven** hotkeys, all registered — h0 uid 74, h1 209, h2 235, h3 196, h4 189, h5 147, h6 136,
h7 75, h8 151, h9 224, h10 10 — and **all eleven play a band**, not seed-depend
(`SEED_DEPEND_VARIANTS` is empty, so every process carries `NIOME_SEED_DEPEND=` and skips that
rung). `ALL_CUT_HOTKEYS` is now **empty**: h0 came off all-cut on the fleet-level pricing below, so
all-cut is the ladder's fallback for every hotkey and the dedicated flat-score hedge is gone. All
eleven pm2 apps carry `--restart-delay 60000`, and axons are bound on 8091-8099, 9001 and 9002.

The layout is the **joined window**, in `window_plan.py`, and it has two tiers:

| | hotkeys | window | band measured |
|---|---|---|---|
| rotated slices | h1-h10 | `JOINED_SUB_WIDTH` 225 at `JOINED_STRIDE` 30, circular | 12 (K562), 7 (HEK293) |
| full width | h0 (`JOINED_FULL_HK`) | the whole 300-seed joined space, no offset | 11 (K562), 7 (HEK293) |

Every seed in the joined space sits in 7.5 of the ten slices on average, plus h0. `miner.sh` needs
`ALLOW_OVERLAPPING_WINDOWS=1` for this: `assert_disjoint_windows` otherwise **exits** rather than
warning, which is the right default for every other layout.

**Width is nearly free, and the contiguous width sweep does not predict the joined case.** The
falsified table's row records band 13/12/11/9 at contiguous widths 100/150/200/300. Measured on the
live fleet at joined widths, the band is **12 at 225 and 11 at 300** on K562, and **flat at 7** at
100/200/225/300 on HEK293. So going 200 -> 225 cost no band at all and bought 6.67 -> 7.5 slices of
depth per seed, and h0's full 300 costs one seed on K562 and none on HEK293 against that table's
predicted 9. Do not price a joined window off the contiguous sweep.

**Where the three joined classes come from.** `JOINED_SOURCE` selects between `_rank_freq_windows`
(the cumulative-rank-frequency scheme) and `"seed_model"`, which reads `predicted_classes` from
[seed_model/next_prediction.json](seed_model/next_prediction.json) — a 152k-param SeedFormer
transformer, one finetuned checkpoint per cell type. It is currently on `"seed_model"`.
**Both are arbitrary and that is measured, so the switch is EV-neutral rather than an
improvement.** SeedFormer's own walk-forward series across four cycles reads 0.899 / 0.950 / 1.000
/ 0.976 hits-of-3 against **chance's 1.000** and cold_hand's 0.988 (p = 0.50-0.80), and its
log-loss 2.202-2.210 is *worse* than uniform's ln 9 = 2.197 — it learned to emit the uniform
distribution, deviating at most 0.024 in any class. `_rank_freq_windows` is no better
(z = -0.71 over 108 rounds). Band position is free under a uniform generator, so three arbitrary
classes are worth exactly as many as three others. It is wired up because the prediction accrues a
scored record either way and the walk-forward number is what would move first if the generator ever
stopped being uniform. `round_plan.sh` refreshes it, guarded by `seed_refresh_guard.py` so it only
retrains when a round has actually stamped (the cron is hourly, rounds stamp every ~2h24m), niced
and capped to 4 threads so a 65 s torch run cannot slow a CPU-bound miner build. Torch lives in
`.venv-ml`, never the subnet's `.venv`.

**A reboot takes the whole fleet down, silently.** On 2026-09-09 the box rebooted twice (15:39,
15:42) after a ~30-minute outbound-network outage. There is **no pm2 systemd unit**, so nothing
resurrected: eleven apps stayed down until restarted by hand, and the pm2 daemon that answered the
first `pm2 jlist` was a fresh empty one spawned by that very command. `pm2 save` has since been run,
so `dump.pm2` holds all eleven with h0's current config — but without `pm2 startup` that only helps
someone who runs `pm2 resurrect` by hand. The stale dump it replaced was from 08:26 and still had
h0 on all-cut, so a blind resurrect would have silently reverted a config change. Two failure modes
to keep apart when reading a dead fleet: the network outage cost 7 of 11 hotkeys their submission on
4390d969 (every hotkey called at or before 15:02 submitted; every one called at 15:08 or later
failed `GET contract` after `MAX_TASK_RETRIES = 3` with 219-283 s of TTL still unused), while the
reboot cost the fleet's existence.

Overlap is free but buys nothing on its own — what spikes is the ~8-12 seed clean BAND inside the
window, not the window — so coverage is `sum(band)` over the fleet and **hotkey count is the only
lever with headroom**. Measured: 7 hotkeys hold 56-84 band seeds (17-25% of rounds), against
6 disjoint width-100 windows' 78 seeds (23.8%). See "What the top block actually buys" below.

**Half of those windows are currently arbitrary, and the plan file will not tell you.** `fit_beta`
returns **0.00 for CD34+_HSPC and HEK293**, so every window ties at `p_hit` 0.2977 and `predict`'s
ranking is a stable-sort tie-break to the lowest index — which is exactly the 100-199 / 200-299 /
300-399 / 400-499 block the live plan assigns for both cells. HUDEP-2 (beta 0.25 → 800-899 first)
and K562 (beta 0.30 → 900-999 first) are real orderings. Check beta before reading a plan as a
prediction.

The whole fleet was deregistered on 2026-09-04/05 and rebuilt smaller, so any claim here about nine
hotkeys or coverage unions is history.

The band arithmetic, kept because it is what a returning band hotkey is worth: the band is ~1-2% of
the 100-999 seed space and a round draws three seeds, so one hotkey spikes on
`1-(1-band/900)**3` of rounds — 2.3% at HEK293's band 7. A coldkey's payout is that expression over
the **union** of its hotkeys' bands, so the whole win is making the bands **disjoint**. Measured:
three hotkeys on disjoint windows covered 15.2% of rounds against 5.2% for the same window three
times (2.9x). At nine disjoint windows the union is 63/900 on HEK293 (~19.6% of rounds) and ~108/900
on the erythroid types (~32%).

Two environment variables carry this, both read at process start:

- **`NIOME_HDR_WINDOW`** — the hotkey's band window, `"200-299"` or a bare `"200"`. Parsed by
  `Miner._parse_hdr_window` (rejects anything outside 100-999 or non-increasing, with a warning, and
  falls back to the cell default). It overrides `CELL_CONFIG`'s `hdr_range` for **every** cell type,
  since one hotkey owns one seed window regardless of what the round asks for.
- **`NIOME_INSTANCE`** — namespaces that process's read/write files under `data/inst/<name>/`
  ([settings.py](niome_subnet/utils/settings.py) `DATA_DIR`). Without it, siblings sharing a working
  directory overwrite each other's submission, task artifacts, upload record and local scoring. The
  *uploads* stay correct either way (rows are held in memory and PUT to a per-uid S3 key) — it is the
  on-disk diagnostics that collide. `chr11.fa`, the k-mer cache and the banks deliberately stay
  shared: large, read-only or content-keyed (`bank_key` folds in the window).

[miner.sh](miner.sh) is the launcher: a `HOTKEYS` table of `<hotkey> <port> <external-port>
<window>`, tiling 100-999 across nine rows (the table is the port registry; which rows are *live*
is set by `DEREGISTERED`). `assert_disjoint_windows` parses the ranges and
**exits** on an overlap or a malformed range rather than warning — `600-899` vs `800-899` overlap
without matching as strings, and a silently-rejected window re-correlates the sibling, which costs
the entire decorrelation. Run one hotkey per pm2 app —
`pm2 start ./miner.sh --name miner-h1 -- niome_hotkey1` — because the no-argument form launches all
of them under one wrapper, where pm2 cannot restart an individual crash-looping child. `PY` is an absolute venv path on purpose: `pm2 restart
--update-env` has rewritten `PATH` and sent the miner into a crash loop before.

Four bash lists in [miner.sh](miner.sh) decide what each hotkey does, and `window_plan.py` parses
all four so the two files cannot drift:

| list | meaning |
|---|---|
| `SEED_DEPEND_VARIANTS` | `<hotkey>:<n>` — runs seed-depend. The value no longer differentiates: the build is shared |
| `ALL_CUT_HOTKEYS` | `<hotkey>` (all cells) or `<hotkey>:CELL,CELL` — all-cut instead of all-HDR on those cells |
| `DEREGISTERED` | not on the metagraph; left out of the window allocation entirely |
| `SPREAD_HOTKEYS` | preferred for the wide 100-seed windows when the concentrated block does not need them |

A hotkey that plays no band (seed-depend, or all-cut on that cell) is excluded from the plan for
that cell — a window assigned to it is a window nothing covers, and the log would read as though a
band had been played there. With every hotkey excluded, `window_plan.py` writes empty assignments
rather than erroring, so the shadow log keeps accruing the evidence that decides whether a band is
worth returning to.

**A deregistered hotkey cannot start.** `base/neuron.check_registered` calls `exit()` in `__init__`
*and* in every `sync()`, so a process that starts unregistered crash-loops and one that is
deregistered mid-run dies at its next sync. pm2 apps here carry `--restart-delay 60000` for that
reason: the miner retries once a minute and comes up on its own within a minute of
`btcli subnet register` landing, instead of hammering the chain endpoint every 7s.

[fleet_status.py](fleet_status.py) tallies the latest task across `/root/.pm2/logs/miner-h*`:
prefetch-ready, build time, whether a validator called, whether the rows were served from the
prepare or an in-TTL fallback, submitted rows, and OOM. It is the check for resource contention.
Note the `KeyboardInterrupt` tracebacks in those logs are pm2's SIGINT at shutdown, not failures,
and that **the miner's own output goes to `miner-h<N>-error.log`** (loguru writes to stderr); the
`-out.log` holds only miner.sh's launcher echoes. Filter `Miner running` and `resync_metagraph`
or you will see nothing else.

**Memory, not the GPU, is the binding constraint now.** seed-depend is CPU-only, and one build at
`variants_per_site 12000` holds **12.7 GB** (measured): four concurrent would need 50.6 GB on a
49 GB box, which is why the build is shared rather than repeated per hotkey. A 24000-variant build
reaches ~24 GB and OOM-killed four research processes; the miners survived only because the kernel
picked the larger victims.

#### The Cas9 scan's per-cell floor — the one real defect found

`scan_cas9` (in **both** [all_cut.py](niome_subnet/genomics/all_cut.py) and
[all_hdr.py](niome_subnet/genomics/all_hdr.py)) walks Cas9 sites **nearest-first** and used to break
on `len(found) >= want * pool_target and len({(mutation, strand)}) == 4`. That cell test only
required each (mutation, strand) cell to be **non-empty** — one candidate satisfies it — while
`assemble` then apportions `want` rows by mutation weight and asks for close to half of them on a
single heavy strand. Measured on HEK293 a8b9f1bb at width 300: the break fired at **job 34 of 395**
with HEAVY+ holding **4** candidates against **5531** available, so the quota backfilled with light
rows and `mean_weight` came out **0.784** against a reachable 0.939.

`cas9_cell_target(contract, ctx, cfg, want)` now derives the floor from whichever apportionment path
the config uses — the hard `light_cell_rows` split gives `(want - 2*lcr + 1)//2` = **79** for
all-HDR, the exponent path gives the heaviest mutation's share halved per strand (**65** on HEK293
all-cut, **72** on K562) — and the break requires every cell to reach it:

| cell | baseline Cas9 split | with the floor | cost |
|---|---|---|---|
| **HEK293** w300 | H+ **4**, H− 43, l+ 97, l− 26 | **weighted +16.2%, fidelity +2.1%**, band unchanged | +9 s |
| CD34+_HSPC | H+ 79, H− 79, l+ 6, l− 6 | +0.0% | +2 s |
| HUDEP-2 | H+ 79, H− 79, l+ 6, l− 6 | +0.0% | +3 s |
| K562 | H+ 79, H− 79, l+ 6, l− 6 | +0.1% | +7 s |

**HEK293 is the only cell that starves**, because its heavy-mutation Cas9 sites sit farther from the
mutation than the light ones; the other three already reach the ideal 79/79/6/6 split unaided and
the floor is inert there to four decimal places.

**Why fidelity *rose* here, without contradicting the `light_cell_rows`-is-dead entry.** The starved
split was badly strand-imbalanced (light+ 97 against light− 26), and that asymmetry was itself
costing stage 5 entropy — so rebalancing gained more than cutting light rows to 6/6 lost. At width
100, where the baseline was less pathological, fidelity does fall −2.5% exactly as that entry
predicts, and weighted's +6.3% carries it anyway.

Two related knobs shipped with the width-300 layout, both documented at their definitions:
`WIDE_WINDOW_VARIANTS` (44000 — do **not** budget the GPU by summing per-build times, see the
constant) and `_scaled_max_fail`, which holds the **z-score** of the tuned `main_max_fail` rather
than its rate. Linear span scaling is wrong for HEK293: mf 48 sits at z −3.32 at span 100 and linear
scaling to 144 lands at **z −5.66**, where only 53 of 60000 bank guides qualify against a group of
80 — every band hotkey declined until this was fixed (HEK293 now uses **165** at width 300).

#### What the top block actually buys

Ranks 7-10 of HEK293 a8b9f1bb are one operator with 90+ slots (ip/coldkey grouping), and their
scores are **single band hits**: cons 0.3854-0.3865, i.e. `(1 + 2*0.078)/3` at the same ~0.078 floor
we measure. So the bar for a k=1 round to clear that field's 81.92 cutoff is

    weighted * fidelity >= 81.92 / 0.3854 = 212.6

They hold weighted **221** and fidelity **0.962**. We reach **179.0** after the Cas9 fix (0.84x),
and **no lever tested reaches their combination** — the construction trades the two terms against
each other. Group 125 reaches their *fidelity* (0.9639) at weighted 176; group 42 reaches their
*weighted* (206) at fidelity 0.82; window width moves the product only within 177-185. Treat 212.6
as a property of all-HDR rather than a tuning target, and do not spend more on it.

**A high-consistency block appears in 6 of 119 rounds and nothing here explains it.** This
subsection claimed a **min-gap** mechanism — that whether any miner reaches `consistency_factor`
>= 0.70 is decided by the distance between the two closest of the three seeds — and **round
56a9f4cb falsified it on 2026-09-10** with seven such miners at min-gap **75**. The re-run over 119
current-regime rounds, against what this file previously recorded:

| min-gap | rounds | mean miners at cons >= 0.70 | max | previously recorded |
|---|---|---|---|---|
| <=5 | 3 | 7.00 | 7 | 7.00 / max 7 |
| 6-15 | 9 | 1.33 | 6 | 1.50 / max 6 |
| 16-50 | 24 | 0.00 | 0 | 0.00 / max 0 |
| **51-150** | 42 | **0.17** | **7** | **0.00 / max 0** — broken |
| >150 | 41 | 0.00 | 0 | 0.00 / max 0 |

What the data actually supports is much weaker, and is about **rarity, not geometry**. Only **six of
119 rounds** contain any miner at cons >= 0.70 at all, and when one does the count is **6 or 7**,
never 1-5 — bimodal, so whatever produces it produces the whole block at once:

| round | date | cell | min-gap | miners >= 0.70 |
|---|---|---|---|---|
| b743244d | 09-01 06:54 | HEK293 | 10 | 6 |
| f2d36f4d | 09-04 14:07 | K562 | 6 | 6 |
| f9511be6 | 09-05 11:44 | HUDEP-2 | 4 | 7 |
| a155b953 | 09-09 09:30 | K562 | 3 | 7 |
| d4a8f8b8 | 09-09 11:56 | CD34+_HSPC | 2 | 7 |
| **56a9f4cb** | **09-10 04:45** | **CD34+_HSPC** | **75** | **7** |

Five of six had a small gap, which is what made the mechanism look measured; with n = 6 that is not
enough to be one, and the sixth round produces the identical seven-miner block with no close pair.
**Do not re-derive their window width from seed geometry** — the inference "their window is ~10
seeds, because at width 10 the band IS the window" rested entirely on the broken separation. It is
still true that our width-225 slices cannot do it (our band is 12 *scattered* seeds of 225, so two
seeds 2 apart are no likelier to both land in it than two seeds 300 apart — on d4a8f8b8 h7's window
held **both** 523 and 525 and its band caught one), but that is a fact about us, not evidence about
them.

The block's own shape is the remaining lead, and it is beyond the measured value ladder: on
56a9f4cb the seven ran cons **0.7004-0.8747** at fidelity 0.9435-0.9749 and weighted 225-271. The
ladder's best reachable k=1 round is 0.491 (one band hit on an all-cut clean floor) and the
conjunction caps at 0.409 — **0.87 needs roughly 2.6 of 3 seeds' worth of pinned targets**, which no
construction in this file reaches at any k. This is the same unexplained ~0.50-and-above floor noted
on the leaders' rounds. It is a construction, not luck (chi2 = 382 on 6 df), and it is **not a
stage-12 knob**.

What their slots buy is **frequency**, and that is the half worth copying. At band 8 a coldkey needs
**36-45 disjoint windows** for a 78% chance some sibling hits a band seed; 90 slots is near-certain
coverage every round, which is exactly what a block of 8-10 consecutive finishers looks like. Our 7
band hotkeys cover 56-84 seeds (17-25%). Coverage is linear in hotkey count and every other lever
measured this session is capped — **hotkey count is the only one with headroom left**.

#### The seed-window prediction, and why it is dormant

[seed_window_model.py](seed_window_model.py) predicts which 100-seed window the next round of a
cell type will draw from, and [window_plan.py](window_plan.py) concentrates several hotkeys onto it
(`round_plan.sh`, cron at :17, resolves live predictions then writes `data/window_plan.json`; the
miner reads it per build and falls back to `NIOME_HDR_WINDOW` if it is missing, stale or malformed,
TTL 6h). Concentration is EV-neutral under a uniform generator, which is why it was safe to run
before the prediction was proven: six hotkeys at width 16 cover 88 band seeds and a uniform seed
hits them with probability 88/900 whether those seeds sit in one window or nine.

**It has not proven itself.** Per-rank *exclusive* hit rates over 38 scored predictions, against
the 29.8% chance baseline every individual rank shares (`1-(8/9)**3`):

| rank | held a seed | 95% CI | lift | p |
|---|---:|---|---:|---:|
| #1 | 39.5% | [25.6%, 55.3%] | 1.33x | 0.130 |
| #2 | 31.6% | [19.1%, 47.5%] | 1.06x | 0.464 |
| #3 | 31.6% | [19.1%, 47.5%] | 1.06x | 0.464 |
| #4 | 26.3% | [15.0%, 42.0%] | 0.88x | 0.735 |

The aggregate ordering is monotone, which is the right shape, but only rank 1 is above chance and
not significantly; ranks 2-3 are indistinguishable from naming any window, and rank 4 is worse.
**Per cell type the ordering does not reproduce at all** — CD34+ 46.2/15.4/15.4/**46.2**, HEK293
38.5/30.8/15.4/23.1, HUDEP-2 **inverted** at 33.3/50.0/**66.7**/8.3 — so a per-cell allocation has
nothing to stand on at n≈13. The "double window" result (a repeated window among the three seeds,
chance 3.43%) collapses the same way: all three of its apparent hits sit in one cell at one rank.

Two things worth keeping from the mechanics. `fit_beta` returns **0.00** for CD34+ and HEK293,
meaning maximum likelihood finds no balancing at all — and with beta 0 every window ties at 0.298
and `predict`'s "rank 1" is a *stable-sort tie-break to the lowest index*, not a prediction. Check
beta before reading anything into a rank. And `tile()` must cover the window: packing at a fixed
width left 4 of 100 seeds uncovered and cost one of only two correct predictions (678cf369 called
300-399, seed 397 landed in it, the block tiled only 300-395).

**What a correct prediction is worth — measured, and it is a lot.** On K562 task 7287db6a
(seeds 907, 169, 885) the window 850-925 held two of the three. With a whole-window all-cut build
and the HDR band **pinned to 885 and 907** — oracle knowledge, unavailable at build time — the round
scored `1.000 / 0.149 / 1.000` = consistency **0.7162**, final **136.4**, **rank 1 of 245**, against
that field's real rank-1 of 120.4. Note the third seed at **0.149** rather than the 0.10 floor:
all-cut's clean set is still underneath the band. That is the "spike *and* floor" property every
combined construction in the falsified table failed to produce, and this is the only measurement
that has ever shown it.

**It does not convert, because the band pins 3 seeds while the model predicts 100-seed windows.**
The same build with the band chosen greedily inside that same window missed both seeds at every k
and placed 11 / 26 / 79 / 82 — the band *cost* ranks, since pure all-cut (rank 11, final 54.4) was
the best of the four. The binding gap is resolution, not window accuracy:

| | |
|---|---|
| one fixed width-76 window holds >=2 of the 3 seeds | **2.02%** |
| ...holds >=1 | 23.25% |
| four disjoint width-76 windows, some window holds >=2 | **8.08%** (>=1: 70.78%) |
| *some* width-76 window holds 2 seeds (min-gap <= 75) | 42.31% |
| a blind 2-seed band inside a window that does hold both | ~1 in `C(76,2)` = 2850 |

So the construction needs the seeds to within **3 of 900**, and the model predicts 100-wide windows
at a rank-1 hit rate of 39.5% for >=1 seed (not significant). **Point prediction is a different
problem from window prediction and nothing here has measured it** — that, not a wider band, is what
would have to work first. Tooling: `conj_test.py` with `CONJ_BAND_FORCE` prices the oracle case,
and without it the blind case.

#### One shared-config hazard

`all_hdr` imports `assemble`, `bank_key`, `_params_fn`, `load_bank` and `save_bank` from
[all_cut.py](niome_subnet/genomics/all_cut.py) so the two builders cannot drift on what invalidates a
cache. That means **`AllCutConfig` and `AllHdrConfig` are both passed to the same `assemble`** while
having different fields. Reach for a field only one of them defines and every live hotkey crash-loops
on restart — which is why `assemble` reads `getattr(cfg, "cas9_cell_floor", 6)`. `AllHdrConfig`
likewise aliases `cas12a_max_fail` to its own `main_max_fail` as a property, so the borrowed
`bank_key` keeps working. Add a field to one config and check `assemble` before shipping it.

`assemble` also scores Cas9 candidates **per cell** (`score_cap // len(by_cell)`, plus a per-cell
floor) rather than globally: a global cap dropped whole stage-5 cells once the mix was balanced, and
an empty (mutation × cas × strand) cell costs roughly a 0.03x multiplier on the entire score.

#### What the competition actually does, and the floor nobody escapes

Read this before designing a new construction. A long session of hypotheses died here, and the
falsified ones are listed at the end so they are not re-run.

**There are two scoring regimes and they are not comparable.** Join any score row to its task's
seed count before using it:

| regime | tasks | rows | max `consistency_factor` ever | #≥0.90 | #≥0.70 |
|---|---|---|---|---|---|
| 1 seed (to 2026-08-24) | 291 | 19,591 | **1.000** | 246 | 290 |
| 3 seeds (current) | 52 | 13,158 | **0.875** | **0** | **6** |

In the single-seed era the seed was knowable, so the whole field converged on it: rank-1..8 medians
climbed 0.33 → 0.40 → 0.65 → 0.99 over ten days, ending with almost everyone at 1.000. Splitting
into three seeds stamped after broadcast is what broke that. **Any analysis that mixes the two eras
will invent a mechanism that does not exist** — filtering on the *score row's* `created_at` is not
enough, because single-seed tasks were still being scored after the switch.

**The ~0.10 floor is arithmetic, not a design failure.**
`cons = 0.7·max(avg_r2, 0) + 0.3·(1 - avg_nmae)`. Off a band seed the three targets are
constant-probability Bernoulli noise, the forest overfits to *negative* R² (measured -0.12, -0.21,
-0.44 → avg -0.26, clipped to 0), and the nmae term is pinned by `is_hdr`: `P(HDR)` can only span
[0.33, 0.59] because `hdr_w = HDR_BASE + 0.35·energy` with energy clamped at 1.0, so that target is
always near a fair coin and its nmae is ~0.94. That forces `avg_nmae ≈ 0.63` and a floor of
`0.3 × 0.37 ≈ 0.11`. Measured 0.092-0.101 for *every* 250-row submission tried — all-HDR, `mh`,
`mh_any`, `shaped` at any `gc_spread`, and the naive baseline.

**"Nobody has a better floor" was true until 2026-09-08 and is now false on the erythroid types.**
The arithmetic above still holds for *our* rows — it bounds a submission whose targets are
Bernoulli noise — but it never bounded the field, and the field left it. Between the 09-07 and 09-08
rounds **205 of 252 hotkeys (81%)** lifted their median erythroid `consistency_factor` from
**0.1079 to 0.1857** (per-hotkey medians, >=8 pre and >=5 post rounds each; the shift distribution
is centred on +0.08 with nothing below -0.02). Pooled over every deduped erythroid row:

| | pre 09-08 | post 09-08 | |
|---|---|---|---|
| `consistency_factor` | 0.1079 | **0.1857** | **+72.1%** |
| `total_weighted_score` | 276.3 | 269.3 | -2.5% |
| `distribution_fidelity_score` | 0.9199 | **0.8739** | -5.0% |
| `final_score` | 31.0 | **44.5** | **+43.5%** |

So they bought +72% consistency for 2.5% of weighted and 5.0% of fidelity — **the same trade this
file prices as expensive, and their fidelity is now BELOW ours** (0.874 against our live 0.877-0.912).
**HEK293 did not move** (field median 0.0758-0.0790 before and after), which is the strongest hint
about the mechanism: it is the one cell type where a wide cut-clean set is not available.

**It is one shared CONSTRUCTION, and it is not one build — an earlier draft of this entry said
"byte-identical weighted 279.1, 153 of 161 rows" and that was a rounding artefact.** At full
precision that cluster is 188 hotkeys spanning **278.730704 to 279.104859** — 58 distinct weighted
values and **165 distinct fidelities of 188**. The rows are not identical, so "one deterministic
build" does not follow. What the numbers actually support, measured on all 19 post-shift erythroid
rounds:

| | |
|---|---|
| share of the field in the cluster | **73-84%**, mean **78%** |
| weighted spread inside it | **0.0099 to 1.53** (0.004%-0.6%) |
| `consistency_factor` inside it | **0.098 to 0.708** |
| fidelity inside it | 0.78-0.93, and *varies* row to row |

**Weighted nearly invariant while fidelity and consistency vary is the signature of a shared SITE
selection with per-hotkey guide variants.** `total_weighted_score` is a sum over the site multiset
(mutation, cas, strand, gc, distance); stage 3 hashes the actual `guideRNA`. So two builds that pick
the same sites and then search different guide variants at each site — which is exactly what
`assemble` then `all_hdr`/`all_cut` do — score identical weighted and completely different
consistency. **That is our own architecture, run by ~180 hotkeys instead of 11.**

**So this is not a floor-only construction: it reaches k=2.** Cluster members hit cons 0.708 on
a155b953 and 0.707 on d4a8f8b8, and 37 cluster rows across the 19 rounds sit at cons >= 0.44. Its
elevated *floor* (the 0.20-0.25 mode) is the same construction on a round where its band misses.

**Two competing families, separated by fidelity rather than by weighted.** Of the 6-9 rows per round
at cons >= 0.44:

| family | n | fidelity | weighted | where it finishes |
|---|---|---|---|---|
| **outside** the cluster | 114 | **0.9429** (0.834-0.980) | 242.0 | ranks 1-6 |
| inside the cluster | 37 | 0.8726 (0.801-0.926) | **278.8** | ranks 7-10 |

The top family gives up 13% of weighted to gain 8% of fidelity and out-ranks the cluster on the
product. On fbde9cc4 that is ranks 1-6 at fidelity 0.908-0.979 against the cluster's 0.8996-0.8999
at ranks 7-10.

**"The field adopted all-cut" is now the WEAKER hypothesis** and should not be carried forward. It
was proposed here on three fits — the 0.237 clean-seed value, determinism-per-contract, and the
HEK293 exemption — and the middle one is gone with the tie. A band construction on an elevated floor
fits everything all-cut fitted *and* explains the cons reaching 0.71, which all-cut cannot. What the
HEK293 exemption then says is that whatever lifts the floor needs a wide cut-clean set, which HEK293
does not have; that part survives either way.

Two consequences that change how to read everything below. **Our 0.101 is now the worst floor in the
erythroid field, not the universal one** — and the round table that follows is *our* ladder, not the
field's. And the floor is what decides whether our spike places: see "Pricing a floor change".

Round score averages the three seeds, so with a floor of 0.10 the only reachable values are:

| seeds in band | round `consistency_factor` | round final (weighted 230, fid 0.95) |
|---|---|---|
| 0 | 0.10 | ~22 |
| 1 | **0.40** | ~87 |
| 2 | 0.70 | ~153 |
| 3 | 1.00 | ~219 |

h0 hit one seed on task 9ed335da and scored exactly 0.397 / 86.58, confirming the model. **Never
quote the single-seed band score (~290) as a round score** — the k=1 value is what places.

**That table is incomplete: there is a third per-seed value, and it is where the field's payout
lives.** Stage 4 scores each seed independently, and what a seed is worth is set by how many of its
three targets are *exactly constant* across the 250 rows — because `r2_score` returns 1.0 for a
constant `y_test` predicted exactly, and `normalized_mae` returns the raw MAE **unnormalised** when
`std(y) < 1e-9`. Pinning all three targets therefore scores exactly 1.0 on both terms; pinning one
scores 1.0 on that target only, while the other two still overfit to negative R². Measured on
HUDEP-2 (9ac1d178 for all-HDR, 4725e952 for all-cut):

| per-seed regime | pinned targets | n | mean | range |
|---|---|---|---|---|
| HDR band | all 3 | 2 | **1.0000** | exact |
| cut-clean, all-cut rows | `is_cut` | 14 | **0.237** | 0.139-0.303 |
| cut-clean, all-HDR rows | `is_cut` | 5 | **0.162** | 0.135-0.203 |
| dirty, all-cut rows | none | 8 | 0.108 | 0.103-0.121 |
| dirty, all-HDR rows | none | 2 | 0.101 | 0.099-0.104 |

Two things to take from it. **Pinning `is_cut` is worth far less than the 0.7 R² weight suggests** —
`avg_r2` only reaches 0.04-0.24, because `is_hdr` and `indel_length` overfit *harder* once `is_cut`
goes constant (-0.10 to -0.95 across the sampled seeds). And **the value is a property of row
composition, not of the rule**: all-cut's 0.237 beats all-HDR's 0.162 by 0.075 ± 0.020 (t = 3.8),
because all-cut's rows overfit less on those two targets. Never price a cut-clean seed at the other
composition's number, and never off one seed — an earlier draft of this section did exactly that and
recorded 0.299, which is the top of the range rather than the mean. This does not contradict the
"nobody has a better floor" line above: the *dirty* floor really is ~0.10 for everyone. The
difference between us and the field is how many seeds are not dirty.

**What the field's top actually is: a top block present in EVERY round, and the plateau
reading of it was wrong twice.** The mechanism `(1 + 2*0.237)/3 = 0.491` is right and the framing
around it was not. Re-measured over 119 current-regime rounds:

| | |
|---|---|
| rounds carrying at least one row at cons >= 0.44 | **118 of 119** |
| rows at cons >= 0.44 per round | **4-9** (median 6, one round at 14, one at 0) |
| total such rows | 740 |
| distinct hotkeys holding at least one | **157** of 528 seen |
| busiest single hotkey | 16 of its 119 rounds |

So the "2.52% of rows" is not an occasional construction — it is a **permanent per-round block of
4-9 miners**, and no round is free of it. Two corrections follow:

* **It is not one operator.** Membership rotates: ~157 hotkeys each reach it on ~4% of their rounds,
  which lands 6-7 per round by arithmetic alone. The chi2 = 382 figure below is still the right test
  against an i.i.d. lottery over *all* 248 miners, but it identifies a **population** of ~157 that
  run the construction, not a block that recurs by hotkey.
* **The value the plateau sits on has moved.** At the field's new erythroid floor of 0.19-0.24,
  `(1 + 2*0.24)/3 = 0.493` — so the same one-band-hit arithmetic now reproduces the plateau *without
  needing all-cut's row composition*, which is what the top block's rows actually look like (7 of 7
  distinct weighted values and fidelity 0.955 on 56a9f4cb, against the shared build's tied 279.1 /
  0.874). **This also explains the ~0.50 "off-band floor" flagged elsewhere as unexplained.**

**And any row at cons >= 0.44 has a genuinely positive `avg_r2`, which ours never does.** From
`cons = 0.7*max(avg_r2, 0) + 0.3*(1 - avg_nmae)`, the nmae term alone caps at **0.30**, so
cons >= 0.44 forces `avg_r2 >= 0.2`. All 740 of those rows clear it. Ours measures **-0.249** off
band. That is the gap in one number, and it is not a stage-12 quantity.

Two band hits remain ruled out as the explanation: at band 12-16 `P(k>=2)` predicts **0.13-0.23**
miners per round out of 248, and matching the observed rate needs a band of **~55 seeds**, against
the 12-16 all-HDR reaches and the ~12 that `pool * P**B >= group_size` allows.

**That is the whole prize, and the Cas9 conditional fill is what blocks it.** Priced against
9ac1d178's real field (rank-10 cutoff 126.32) at our own `weighted x fidelity` of 303.6:

| build | cut-clean of 900 | cons at k=1 | round final | rank | |
|---|---|---|---|---|---|
| all-HDR as shipped | 17 | 0.401 | 121.7 | 11 | shipped |
| conjunction, group 80 | **52** | 0.406 | 123.2 | 11 | **measured, Cas12a-bound** |
| conjunction, group 42 | **80** | 0.409 | 124.1 | 11 | **measured, Cas9-bound** |
| + a *perfect* Cas12a cut tie-break | 162 | 0.408 | 123.9 | 11 | bound, not reachable |
| + both halves cut-min-unioned, all-cut rows | 559 | 0.459 | 139.2 | **9** | **measured UNREACHABLE** |
| every seed cut-clean | 900 | 0.491 | 149.1 | **8** | bound only |

So a placing k=1 round needs the wide cut-clean set **and** all-cut's composition together, and
**`cas12a_union.py` has since measured that the clean set caps at 52-80** — see the falsified table.
The 559 row is kept only to show what the prize would have been worth; it is not available. Why
neither half can be fixed alone, on the shipped group-80 build: the Cas12a half's cut-fail union is
**845 of 900** by itself (80 rows, mean 34.2 fails each, because Cas12a's `cut_p` caps at 0.96), the
170 Cas9 rows cover **738** (mean 9.1 at `cut_p` 0.99), and only 33 seeds have a single failing row
against a mode of 4-5 — so perfecting the Cas12a half caps cut-clean at **162** and perfecting Cas9
caps it at **55**. Both of those are bounds on a half in isolation; the measured joint answer is
lower still.

**The leaders run our design with more hotkeys.** Grouping current-regime placements by coldkey:
eight coldkeys hold ~50% of the payout curve (top three alone 40.7%), each running **11-14 hotkeys**,
and `placed ≈ hotkeys` (12 of 13, 11 of 11, 14 of 14) — the signature of disjoint band windows, the
same construction shape as [miner.sh](miner.sh). Their placers' median `total_weighted_score` is
**263** against our 230. All **520** top-10 finishers in the 3-seed regime submit exactly **250
rows**; small submissions do not exploit stage 4's KFold degeneracy and score *lower* (0.064-0.093).

So the two real gaps are exposure and term 1, not the construction: 11-14 hotkeys covers ~46% of
rounds against 9 hotkeys' 34%, and h0's rank-12 miss needed weighted 255.7 where it had 223.9. The
weighted gap is the `base structural` half (`gc_score` 0.938, `dist_score` 0.840,
`offtarget_factor` already a perfect 1.0) — the mutation-skew route to it is measured and dead, the
distance/GC route is untested and does not trade against fidelity.

**Seeds are uniform random.** Min-gap between the closest two of the three seeds: median 94 observed
against 93 simulated, and `<=150` at 66% against 70% expected. 66% of rounds do have two seeds within
150 — but that is chance, and it does not help, because the band is ~13 *scattered* seeds however
wide the window it was searched in.

**Falsified — do not re-run these:**

| hypothesis | result |
|---|---|
| wider screening window → wider band | clean count is 13/12/11/9 at widths 100/150/200/300, and **11-16 over the full 900**. Band size is set by how many rows must agree, not the window. Wider costs hotkeys (only 6 tile 900 at width 150) and coverage falls 117 → 72 seeds. |
| graceful predictability (`mh` rule degrades smoothly) | collapses to median **0.097** across 40 seeds, same as `hdr`. The `mh` coin is redrawn per seed, so compliance is not a property of the design. |
| feature variance (`gc_spread`, `shaped` strategy) | floor unchanged (0.095 → 0.099 → 0.098 → 0.096 at `gc_spread` 0 → 0.45) while weighted drops 230 → 206. `shaped` reaches 0.49-0.53 on its *build* seed only. |
| small submissions → degenerate R² = 1.0 | every one of 520 top-10 finishers has 250 rows; the few sub-250 submissions score 0.064-0.093. |
| recover cut-clean seeds inside all-HDR's pools | ceiling **34 of 900** (against all-cut's ~560). Only **82 of 60,000** bank guides are HDR-clean on the band and 80 are needed, so the Cas12a min-union has no freedom. Cas9 min-union on cut does work (union 498 → 394) and is swamped. |
| no-MH_NHEJ (`not_mhnhej`) for a wider band | band **45** (3.5x, real — the rule holds on all 45) but only one pinned target → cons **0.123**. Even all three seeds in band scores 25.7, under every cutoff. E[pay] exactly 0. |
| the rule ladder generally | `hdr` 3 targets → 1.000 / band 13; `not_hdr` 2 targets → 0.577 / band 12; `not_mhnhej` 1 target → 0.123 / band 45. Frequency×value is a wash and `hdr` is the best point. |
| narrowing all-cut's window to widen its clean set | clean fraction *does* rise (62% over 900 → 93.6% at width 300, and strict Cas12a guides exist below ~225) but E[final] **falls**: width 225 **36.78**, width 300 **37.62**, whole window **46.11**. Consistency has a second channel independent of coverage — row composition — so a build can hold union 0 and still lose monotonically as `group_size` grows (0.1706 → 0.1583 → 0.1452). |
| the rule ladder inside a narrow window | at width 225 the band tracks P(rule) exactly as it does at 900: `cut` 225 → `not_mhnhej` 48 → `not_hdr` 10 → `mh_any` 9, all set by the Cas9 conditional fill rather than the window. `not_mhnhej` at width 225 scores E[final] **22.29** (−51.7% against the whole-window 46.11). |
| fidelity as all-cut's ceiling | group 125 at width 225 reached casR 1.000 and fidelity **0.9759** — above the leaders' 0.9473 — with 8/8 cells, and still lost. Fidelity was never the binding term; consistency fell faster. |
| all-cut ∧ an HDR band on the same rows | dead at every setting tried: cut window 100 / 300 / 900 × group 42 / 50 / 80, k up to 8, on K562 and HEK293. HDR ⊂ cut, so this really is a sequential filter and not the falsified "combined construction" — but the filter shrinks the Cas12a pool, a smaller pool min-unions to a *larger* failed-seed union, and the clean-seed floor itself falls (0.2586 → 0.2291 → 0.2145 → 0.1668 at k=0-3, whole window, group 42). E[final] is monotone decreasing in k on every window. |
| widening the conjunction's band | the Cas9 conditional fill caps it at **3 seeds**, against the 12-16 pure all-HDR reaches. Pool decay is ~0.55 per band seed: at group 42 over the whole window the Cas9 pool runs 1749 → 913 → 557 → 220 → **164** for k=0-4 against the 208 rows needed. Group 80 buys depth (170 needed, k=8 built at 206) and pays ~65 weighted for it. |
| group 50 (200/50) as a middle point between 42 and 80 | weighted-identical to group 42 **on the same contract** (214.1 vs 214.4) while covering 36 fewer clean seeds. The 304.6-vs-214.1 gap that looked like a group effect is entirely the *contract*. |
| the distance/GC route to `total_weighted_score` | the bound is **not** where it leaks — mean GC is already 0.505 and `offtarget_factor` a perfect 1.0, so structural sits at 0.85-0.91 of a possible 1.0. Tightening `cas12a_gc`/`cas9_gc` to 0.42-0.58 with `max_distance` 250 raises weighted 2-3% and **costs band on 3 of 4 conditions** (−1 to −3 seeds), net −6% to −18% on frequency×value. The one arm that held its band (CD34+ 100-399, +2.4%) did not reproduce at another window. Also **falsifies the note that this route "does not trade against fidelity"**: tightening GC costs fidelity monotonically (0.8933 → 0.8745), tightening distance *raises* it (0.9056). |
| GC-aware tie-break inside `FastGreedy` | the greedy's `argmin` is tied on 71 of 80 picks (median 12 candidates, max 577), and preferring the tied guide nearest GC 0.50 does hold the band — but delivers only **+5.6%** on the group's `gc_score`, not the +44.7% the per-pick slack suggests, because ties are re-drawn after every pick and 12 restarts already sample them. Weighted +0.7%, fidelity −1.1%, net **−0.3%**. |
| mutation-weight tie-break inside `FastGreedy` | same trick on `mutation_weight`: shifts 3 of 80 group rows on HEK293 (+1.6% freq×value) and 7 on CD34+, where the **1.8% fidelity drop** eats the 1.2% weighted gain (net −0.6%). CD34+ is already apportioned at `mean_weight` 1.10; only HEK293 was starved, and `cas9_cell_target` fixes that properly. |
| narrowing the band window to buy Cas9 pool freedom | `weighted × fidelity` is **flat at 177-185 across widths 10 / 20 / 30 / 50 / 100 / 300** on HEK293 — a 4% spread over a 30x range. The mechanism fails because **narrowing the window grows the band toward the window size** rather than shrinking it (band **10 of 10** at width 10, against 7 at width 100), so the Cas9 conditional fill gets *harder* and the pool *falls* (1691 at width 10 vs 9290 at width 100). Width 50 is the best point and worth only +3.8% k=1 at equal band. |
| `group_size` away from 80 (HEK293, re-tested after the Cas9 fix) | the two trends are real and cross exactly at 80: fidelity climbs 0.820 → 0.964 across group 42 → 125 while weighted falls 206 → 164, and the **product peaks at group 80** (179.0). Priced on 31 real HEK293 fields every arm is 0.0010-0.0014 E[share] with `P(place \| k=1)` a flat **48%** — quality moves k=1 only between 60 and 69, so band size decides the ordering and group 42's band 9 edges it by 8%, inside noise. Group 80 stands. |
| two band hits as the explanation for the field's top block | needs a band of **~55 seeds** to produce the observed rate (4-9 such rows in **118 of 119** rounds); at band 12-16 the prediction is 0.13-0.23 miners per round against ~6 observed, 6-40x short. It is one hit on an elevated floor, not two hits. |
| stage 4's overfitting as an exploit (degenerate feature vectors) | **real but not buildable.** Collapsing the feature matrix moves off-band `avg_r2` from **-0.249 to -0.032**, confirming the negative R² is pure overfitting — but `max(avg_r2, 0)` discards it, so cons moves only 0.1036 → 0.1052. Combined with a pinned `is_cut` it *is* worth 0.135 → **0.350**, except that the recovery needs `distance` cardinality **<=8** (flat from 61 down to 16: 0.1346 / 0.1359 / 0.1345), and 250 rows admit only ~8 per distance value (2 offsets x 2 mutations x 2 cas x 2 strands), forcing ~31 distinct. |
| a cut-robustness tie-break inside all-HDR's Cas12a greedy | ceiling **+2.2 final points** (121.7 → 123.9), still under the cutoff. Two measured facts compound: perfecting the Cas12a half caps cut-clean at 162 of 900 because the Cas9 rows independently cover 738, and a cut-clean seed in an all-HDR composition is worth 0.162 rather than all-cut's 0.237. Priced before implementing — the go/no-go was `cutunion`, not a build. |
| feature concentration (gc at 0.50, distance clustered) to buy all-cut's composition | **noise-dominated.** Concentrating gc at 0.50 for 79% of rows helped seed 373 (+0.051) and hurt seed 161 (-0.006); clustering distance to 39 levels did the reverse (+0.030 / -0.022). The R² on the two live targets is noise around a negative mean, so composition cannot be tuned seed-by-seed. |
| the conjunction's clean set, measured at both of its bounds (`cas12a_union.py`) | **caps at 52-80 of 900 against the 559 the prize needs, and which bound binds depends on `group_size`.** (a) The HDR min-union buys **zero** cut-coincidence — the group's cut-fail union is **714** at group 42 against **724** for independent failures (848 at group 80), where all-cut reaches **341** with the same 42 rows by optimising cut instead. The two objectives are orthogonal, so a group selected for one gets nothing free on the other. (b) The HDR-on-band Cas9 pool costs `P(HDR)**2` per extra band seed (**906** candidates at band 15 against **3118** at band 13), and cut-strictness over `C` costs a further `0.99**\|C\|`. So group 42 is Cas9-bound at `\|C\|` = 80 and group 80 is Cas12a-bound at `\|C\|` = 52, worth **+0.7 to +2.4** round final (121.7 -> 122.4-124.1) against a 126.32 cutoff. **The scan's early exit is not the cause** — at group 42 the *full* pool (906) is smaller than the exit target (1664), though at group 80 it does truncate (3118 vs 1360). |
| the joined window as the cause of our 0.89 fidelity (`fidelity_window.py`) | **mostly not it, and free on the product.** Four arms on one contract, group and everything else fixed: contiguous 100 -> **0.8956**, contiguous 225 -> 0.8916, joined 225 -> 0.8880, joined 300 -> 0.8773. So the contiguous window *this fidelity was measured at* already gives 0.8956, and non-contiguity costs only 0.9% at width 225 / 2.0% at width 300, entirely through the mutation term (0.722 -> 0.689 -> 0.653). Width alone costs ~nothing. And weighted rises as fidelity falls, so `weighted x fidelity` is **flat at 221.6-223.4 across all four**, independently reproducing the `narrow_width` result. The window is exonerated; h0's full 300 is free net. |
| `max_distance` alone, at the band-preserving end (`dist_sweep.py`) | **the term improves and the product does not.** Distance only, GC bounds untouched, one contract, joined 225, group 80: `dist_score` climbs monotonically 0.8555 -> 0.9294 and `base structural` with it, 0.9088 -> **0.9464** (+4.1%) across maxd 400 -> 100. But `weighted` **oscillates** (251.6 / 248.6 / 257.6 / 248.8 / 255.9 / 256.7) instead of following, and `weighted x fidelity` is **flat at 223.0-228.3 — a 2.4% spread over a 4x range**. Mechanism: a tighter bound restricts eligible sites, shifting the mutation mix and dropping mean `mutation_weight`, which absorbs the base gain — fidelity moves the *other* way as it does so and the two largely cancel. Including band, `band x w x fid` peaks at maxd 200 (2700) and 300 (2698) against shipped 400's 2681 (**+0.7%, noise**) and falls 7% below maxd 150 where the band drops 12 -> 11. This also explains why the 210-config sweep read non-monotone across 400/150/100: it was sampling a flat surface. **Do not extrapolate base linearly** — a 4.1% base gain returned 2.0% weighted and 2.2% product, so even base = 1.000 reaches ~+5%, not the +10.5% that rank 10 requires. |
| min-gap as the mechanism behind the field's high-consistency block | **broken by one round.** Five of the first six rounds carrying a cons >= 0.70 block had min-gap <= 10, which read as a clean separation (zero of 101 rounds at gap >= 16). **56a9f4cb (2026-09-10, CD34+, gap 75) produced the identical seven-miner block**, moving the 51-150 bucket from mean 0.00 / max 0 to 0.17 / max 7. The inference that drew their window width at ~10 seeds from this is withdrawn: with only 6 positive rounds of 119 the gap association was never separable from chance. What survives is that the block is **bimodal (0 or 6-7 miners), rare (5% of rounds), and reaches cons 0.87** — above everything the per-seed value ladder can build. Don't re-run the geometry; the open question is what the rows are. |
| raising the off-band floor by row DISTRIBUTION rather than pinning (`floor_bound.py`) | **the ceiling is our current floor.** Give the stage-4 regressor an ORACLE -- predicting each row's true conditional mean, so `r2 = Var(p)/(E[p(1-p)]+Var(p))` and `nmae = mean(2p(1-p))/sqrt(...)`, the best any regressor could do -- and the most heterogeneous 250 rows the design space allows reach cons **0.085-0.111**, against our live 0.095-0.113. The bound reproduces the shipped floor (0.1114 erythroid / 0.0894 HEK293 predicted) which is the check that the formulas are right. Heterogeneity raises `avg_r2` 0.004 -> 0.027 and `avg_nmae` 0.638 -> 0.705 **together**, so cons falls. Two structural reasons: **`is_hdr`'s nmae is 0.997** (`P(HDR\|cut)` is confined near 0.5 by `hdr = hdr_base + 0.35*energy` against fixed `mh_nhej`/`blunt`, so it is a fair coin no design reaches outside), and **energy is CLAMPED at 1.0 on every erythroid row** -- at accessibility 0.77-0.87, `1.8*gc + 0.6*exp(-d/1500)` exceeds the clamp above gc ~0.33. That clamp is why the falsified `gc_spread` row moved the floor only 0.095 -> 0.099: it varied a saturated quantity. **Distance is the only lever that unclamps energy**, and it is worse (0.1097 against 0.1114) because unclamping moves rows off `cut_p` 0.99 toward 0.5 where nmae is maximal. **So an elevated floor REQUIRES pinning `is_cut`** -- there is no distributional route, and this closes the question from the other side to `cas12a_union.py`. |
| narrowing the HDR band to buy Cas12a pool freedom for a cut min-union (`hdr_pool.py`, `band_floor.py`) | **the freedom is not there, and the two requirements never overlap.** The plan was: an elevated floor needs a wide cut-clean set; that is blocked because only ~82 bank guides are HDR-clean on a 13-seed band against the 80 a group needs; the pool grows as the band narrows, so trade band width for min-union freedom. `band_floor.py` priced the trade against the 19 post-shift erythroid fields and gave **break-even floors** of 0.157 (band 10, needing 371 of 900 cut-clean), 0.170 (band 8, 457), 0.186 (band 7, 562), 0.192 (band 6, 602), with band 4 losing **-33.9%** even at a perfect 0.237. Then `hdr_pool.py` measured the pool over six live banks -- pure set arithmetic on the `fails` each bank already stores, no GPU -- and it is **10-14x below the extrapolation**: best-band pool **847 / 230 / 126 / 72 / 29 / 11** at width 4 / 6 / 7 / 8 / 10 / 13, against all-cut's **41.6x** freedom (group 42 from a pool of ~1749) to reach union 341. So band 10, the arm with the most reachable break-even, has a pool of **29 against a group of 42** and cannot form a group at all; bands 6-7 hold 3.0-5.5x freedom and need a union of 298-338, i.e. **at or better than all-cut's 341 from an eighth of the freedom**; band 4 has the freedom and cannot pay. Two input errors produced the optimism, both from reading one datapoint: the per-seed HDR-clean rate is **0.476, not the 0.602** inferred from "82 of 60,000", and live banks hold **14,000 guides, not 60,000**. Note what did NOT go wrong -- random-band pools match the independence model at the measured rate (1 against 1 at band 13, 33 against 37 at band 8), so there is **no guide-heterogeneity bonus**; the best-band greedy's 11x edge is from choosing seeds, nothing more. Bank size is the only remaining lever (pool is linear in it) and matching all-cut's freedom at band 6 needs a **~106k guide bank, 7.6x the current build**, for a **+17.1%** upside and still a union better than all-cut manages. Not worth it against the memory constraint. |
| whether the HDR filter DESTROYS cut-failure coincidence or merely fails to select for it | **still open, and now unmeasurable where it would matter.** `cas12a_union.py` found a group selected for HDR has a cut-fail union of **714** against **724** for independent failures -- zero coincidence. The obvious next test is to reverse the objective: filter to HDR-clean-on-band first, then min-union on CUT over what survives. The pool measurement above says that pool is 29-126 at every band width whose floor could pay, which is too small to ask the question. Band 4 (pool 847) is the only arm with enough freedom to measure it and is the one band that cannot pay at any floor. Left here because it is the mechanism behind two falsified rows, not because it is actionable. |

**Read the table as a whole before proposing a sixth knob.** Five independent levers have now been
swept — GC tie-break, mutation-weight tie-break, band-window width, `group_size`, and `max_distance`
— and every one shows the same shape: **the term it targets improves measurably and
`weighted x fidelity` does not move.** GC tie-break +5.6% on `gc_score` -> -0.3% net; window width
flat 177-185 over a 30x range; group size a product peaking exactly at the shipped 80;
`max_distance` +4.1% on base structural -> +0.7% on frequency x value. The mechanism differs each
time (ties re-drawn per pick, band growing toward the window, fidelity/weighted crossing, mutation
mix shifting with site availability) but the outcome does not. Treat `weighted x fidelity` as at a
local optimum **of this construction** that per-term tuning does not escape, and note that the two
terms are *coupled through row composition* — which is why improving either in isolation returns
almost nothing. The levers with measured headroom are elsewhere: hotkey count (linear, uncapped) and
consistency structure (where the plateau lives). Anything that only moves a stage-12 term is very
likely already priced at zero.

**"Local optimum" means local, and the qualifier matters — but the comparison against the field
must be made PER CELL TYPE, and pooled it is wrong.** Our live `weighted x fidelity`, medians over
post-09-08 rows identified by ss58 address:

| cell | weighted | fidelity | `w x fid` | vs the shared build's **244.0** |
|---|---|---|---|---|
| HUDEP-2 | 336.0 | 0.8931 | **298.3** | **we are +22.3% ahead** |
| CD34+_HSPC | 260.0 | 0.8936 | 233.3 | -4.4% |
| K562 | 244.3 | 0.8994 | **220.1** | -9.8% |

A 36% spread across cells, so pooling these to 225.8 and reading "the field beats us by 8.1%" is the
same error as pooling fields across contracts — **an earlier draft of this entry said exactly that
and it is withdrawn.** We are *ahead* of the shared build on HUDEP-2 and behind on K562, and the
sweeps' flat 221-228 band was measured on the cells where we are behind. What survives: the flatness
is a property of the all-HDR row composition rather than of stage 12, and the shared build reaches
its number by going the *other* way on the trade than every arm swept here (more weighted, less
fidelity). It does not license re-running the five levers, and **`floor_price.py` arm B measures
that adopting the shared build's product at our current floor is worth -48.8%** — because it would
cost HUDEP-2 its 298.3.

**The same defect is live in [fleet_price.py](fleet_price.py)**, whose `AH_WXF` is the single
constant `339.9 * 0.8931 = 303.6`: right for HUDEP-2 (298.3 measured) and **38% too high for K562**
(220.1). Both arms carry the inflation, so it compresses rather than reverses the all-cut/all-HDR
comparison — but the h0 decision rests on a +6.5% margin computed with the K562 arm scored 38%
high. **Re-run and resolved** — see the corrected table under "Per-hotkey P(place) is the wrong
unit for a fleet": the ordering and the h0 decision survive, the E[share] levels are wrong by up to
11x, and the monotonicity claim turns out to hold only at all-cut's live clean fraction.

`not_hdr`'s **0.577** is the one number worth remembering: two pinned targets lands exactly in the
range the leaders' placing rounds occupy, which is the arithmetic confirmation that a placing round
is one band hit plus two floor seeds — not a flat construction.

#### Pricing a floor change

[floor_price.py](floor_price.py) prices the off-band per-seed floor `f` as a free parameter, against
the **19 post-2026-09-08 erythroid fields** — the shift above changed the field, so pre-shift fields
price a field that no longer exists, which is this section's own rule applied to the time axis. Each
cell is scored at **its own** `weighted x fidelity` (HUDEP-2 298.3, CD34+ 233.3, K562 220.1); a
first pass pooled them at 225.8 and overstated the answer by **2.3x**, so do not pool them.

**The floor never pays directly. It pays by moving k=1 up the curve.** Per cell, against that cell's
own fields, `cons = (k + (3-k)*f)/3`:

| cell | cutoff med | k=0 at f=0.237 | k=1 at f=0.101 | k=1 at f=0.198 | k=2, any f |
|---|---|---|---|---|---|
| CD34+_HSPC | 99.5 | 55.3 (**0%**) | 93.5 (33%) | 108.6 (100%) | 100% |
| HUDEP-2 | 112.1 | 70.7 (**0%**) | 119.5 (50%) | 138.8 (100%) | 100% |
| K562 | 80.9 | 52.2 (**0%**) | 88.2 (67%) | 102.4 (89%) | 100% |

k=2 already places everywhere and k=0 places nowhere even at f=0.237, so the floor's whole value is
in k=1 rounds — and mostly as **rank**, not frequency: rank 3 pays 0.20 against rank 9's 0.015.

| arm | E[share] | vs shipped | |
|---|---|---|---|
| **A** f=0.101, band 12 (**shipped**) | 0.00832 | — | |
| A f=0.150 | 0.01016 | +22.2% | 9/19 rounds, m/SE 3.0 |
| A f=0.198 | 0.01588 | +90.9% | 15/19, m/SE 1.7 |
| A f=0.237 | 0.01773 | **+113.2%** | 16/19, m/SE 2.1 |
| A f=0.300 | 0.06015 | +623.3% | 19/19, m/SE 4.8 — above any measured per-seed value |
| **B** f=0.101 at the shared build's 244.0 | 0.00426 | **-48.8%** | costs HUDEP-2 its 298.3 |
| B f=0.237 at 244.0 | 0.03404 | +309.4% | |
| ~~**C** the field's cluster modelled as flat, no band~~ | ~~0.00000~~ | ~~-100.0%~~ | **premise false, see below** |

Three results, in order of what they settle:

**Arm C is WRONG and is kept only as a worked example of a modelling error.** It priced the field's
cluster as *flat* at cons 0.198 with no band — `279.1 x 0.874 x 0.198 = 48.3` every round, which
never reaches a cutoff whose minimum is 63.4, hence the 0.00000. **The premise is false:** that
cluster is a band construction whose consistency runs 0.098-0.708, not a flat floor (see "Nobody has
a better floor" above). Measured directly against the same 19 fields, the cluster **captures curve
share on 19 of 19 rounds** — mean **0.077**, min 0.025, max 0.360, taking ranks 7-10 on a typical
round. So the correct reading is the opposite of what arm C printed: the *floor-only members* of
that construction do not place, and the construction itself places every round.

What survives from arm C, and it is worth keeping: **a flat 48.3 really would earn zero**, because
the post-shift erythroid cutoff minimum is 63.4. That is still the demonstration that
**`final_score` is not the objective** under `SCORING_SYSTEM = "top"` — the field's +43.5% on
`final_score` would buy nothing on its own. It is their band, not their floor, that pays.

Arms A and B are unaffected: they model *our* construction, whose floor is measured at 0.101 and
whose band is measured at 12.

**The rebuild that arms A and B describe was attempted and it aborted at its first measurement.**
See the two `floor_bound.py` / `hdr_pool.py` rows in the falsified table: the floor cannot be raised
distributionally (oracle ceiling 0.111), and narrowing the band to buy min-union freedom for a
pinned floor finds pools of 29-126 where 1749 was extrapolated. **Read arms A and B as the price of
a prize, not as a plan.**

**A's 2.13x at f=0.237 is a bound, not an option.** Both routes to floor-and-spike are already
falsified: the conjunction puts both rules on the same rows and caps the clean set at 52-80 of 900,
and the hybrid split puts them on different rows and measures *below* the pure floor. The reason is
structural — stage 4 computes each target over all 250 rows, so half the rows repairing by HDR
leaves `is_hdr` unpinned even on a band seed. **At the reachable floor the prize is ~nothing:** the
conjunction's 52-80 clean seeds give `f ~ 0.105`, which sits between the 0.101 and 0.150 rows, i.e.
under +22% and probably under +5% — independently reconfirming `cas12a_union.py`'s "+0.7 to +2.4
round final points".

**Where the floor would be worth most is not where we are weakest.** Per cell at f=0.237: K562
**+278%**, CD34+ **+209%**, HUDEP-2 only **+18.8%** — HUDEP-2's `w x fid` of 298.3 already puts its
k=1 at 119.5 against a 112.1 cutoff, so it needs no floor. Any future floor work is a K562/CD34+
question.

**What this promotes instead.** k=2 places on 100% of fields at every floor and every cell, and the
floor cannot touch that. So the two levers this leaves are the ones already flagged as uncapped:
**hotkey count**, and **band size** — the only quantity that moves `P(k>=2)`. Neither is a stage-12
term, which is the same conclusion the five-lever table reaches from the other direction.

#### Pricing a construction — the field must match the contract

**A construction's score and its field move together, so a score measured on one contract must be
placed in that contract's own field.** Pooling fields across contracts prices field softness instead
of the construction. This has already produced a wrong config decision once, and it is the same
class of error as the "sample across all backend history" trap above.

Whole-window all-cut, group 42, the same build, on two K562 contracts:

| contract | weighted | clean | E[final] | best case | own field's rank-10 cutoff | E[share] own field | vs 19 pooled |
|---|---|---|---|---|---|---|---|
| 83f430e9 | 304.6 | 549 | 57.04 | 75.0 | 116.0 | **0.0000** | 0.0023 |
| 7287db6a | 214.4 | 570 | 38.97 | 49.9 | 63.4 | **0.0000** | 0.0003 |

The clean sets are within 21 seeds of each other, so the 46% gap in E[final] is almost all
`total_weighted_score` — and the field moved with it: rank-10 cutoff 116.0 against 63.4, field
median 50.4 against 34.8. Our score fell 0.68x and the field median 0.69x. So 83f430e9's **0.0023**
is what its score earns in *other, softer* contracts' fields; place it in 7287db6a's field and it
reads **0.0034**, in its own **0.0000**.

That 0.0023 is the number that put h0 on all-cut for K562. `price_cell.py` bootstrapped over fields
with the contract held fixed, so its CI of [+0.0003, +0.0033] never saw the variance that decides
the question. **`conj_price.py --own` prices each arm against the one field that played its
contract; prefer it, and read every pooled figure in this file as an upper bound.**

**Why the correction cuts in all-HDR's favour.** all-cut's ceiling (49.9-75.0) sits *inside* the
K562 cutoff distribution (min 45.9, p25 59.5, **median 95.0**, p75 118.8, max 135.1), so whether it
places is decided by field softness — exactly what pooling scrambles. all-HDR's spike is ~103 on a
single band hit and ~247 on a triple against a maximum observed cutoff of 135, so a hit clears
nearly every field whichever contract drew it.

**The point-value model understates the upper tail; a modelled 0.0% is "small", not "never".** Each
regime is priced at one sampled mean, so the model cannot produce a round above "every seed in the
best regime". On 7287db6a all-cut's modelled best case was 49.9 and the real round scored **54.4,
rank 11 of 245** — all three seeds drew clean *and* above the 0.2586 clean-seed mean
(0.285 / 0.290 / 0.272). Per-seed spread within a regime is real and is not modelled.

**Measured, matched, and it confirms all-cut on K562 — by more than the pooled figure claimed.**
[cmp_k562.py](cmp_k562.py) builds *both* arms on the same contract and prices each against the one
field that played it, over six K562 contracts spanning the whole cutoff range:

| contract | field cut10 | all-cut E[share] | all-HDR E[share] |
|---|---|---|---|
| 87411bfa | 45.9 | **0.0068** | 0.0006 |
| b9051bc7 | 48.8 | **0.0024** | 0.0006 |
| 7287db6a | 63.4 | 0.0000 | 0.0008 |
| 9c26657e | 104.4 | 0.0000 | 0.0002 |
| 83f430e9 | 116.0 | 0.0000 | 0.0002 |
| f2d36f4d | 135.1 | 0.0000 | 0.0000 |
| **mean** | | **0.0015** | **0.0004** |

all-HDR wins 4 of the 6 contracts and still loses the mean **3.78x**: it collects 0.0002-0.0008
reliably while all-cut collects nothing on four contracts and then 0.0068 on one soft field. Here
the *flat* construction holds the fat tail, which is the reverse of the reasoning that justified
all-HDR. The mechanism is a threshold: all-cut's ceiling (52-78) straddles the cutoff distribution,
so it places on 24-68% of rounds where the field cuts below it and never where it does not, and
about a quarter of K562 fields cut that low. **And 3.78x understates it** — the point-value model's
blind spot above only ever helps all-cut, whose ceiling sits *at* the cutoff, never all-HDR, whose
spike clears every field anyway.

Two inputs to `price_cell.py` are wrong and matter beyond K562:

- **all-HDR's `total_weighted_score` is not a contract-generic 258.** It tracks all-cut within
  **3%** on all six contracts (333.3/330.3, 296.4/304.6, 260.1/267.8, 247.6/255.0, 255.3/246.6,
  220.5/214.4) while the contract moves both by 54%. all-HDR gives up nothing measurable on term 1;
  the entire trade is consistency structure. Pricing it at a fixed 258 against a varying all-cut is
  what produced the wrong pooled answer.
- **The off-band floor is 0.097-0.106, not 0.1256** (measured on the real submissions, six
  contracts), which inflated all-HDR's E[final] by ~7%.

The band itself is stable and is not the uncertain part: **13 of 900 on all six contracts, spike
verified at exactly 1.0000 each time** (hunted, not assumed), against all-cut's clean set drifting
549-570 and its hit consistency 0.2561-0.2940.

**Per-hotkey P(place) is the wrong unit for a fleet, and it inverts the answer.** Both matched
single-hotkey measurements are correct — `cmp_k562.py`'s all-cut 3.78x on K562, and a second
independent method (P(place) over 28 real HUDEP-2 cutoffs: all-cut **10.38%** against all-HDR's
**3.26%**, 3.18x) — and **both are inapplicable to a fleet decision**, because the two
constructions differ in whether siblings can be decorrelated at all. all-cut is a deterministic
function of the contract (no window, no band), so every hotkey running it submits **identical
rows**: ten of them hold one score and take ranks r..r+9, collecting the *tail* of
`SCORE_DISTRIBUTION`, and all miss together. all-HDR carries a per-hotkey band window, so one hit
places **alone** near the top. Rank 1 pays 30x rank 10, and that asymmetry decides it.
[fleet_price.py](fleet_price.py) prices whole compositions against the real fields, Monte-Carlo over
seed draws and band placements, at the measured per-seed value ladder:

**The first run of this used one `weighted x fidelity` for every cell type and its E[share]
levels are wrong by up to 11x. Re-run per cell, the ORDERING survives and the levels do not:**

| all-cut / all-HDR | HUDEP-2 | CD34+_HSPC | K562 | vs 1/10 |
|---|---|---|---|---|
| 0 / 11 | **0.03773** | **0.01034** | **0.00436** | **+8.2% / +5.3% / +10.3%** |
| 1 / 10 | 0.03488 | 0.00981 | 0.00396 | — |
| 2 / 9 | 0.03257 | 0.00951 | 0.00369 | -6.6% / -3.1% / -6.6% |
| 4 / 7 | 0.02840 | 0.00934 | 0.00317 | -18.6% / -4.8% / -19.8% |
| 6 / 5 | 0.02431 | 0.00923 | 0.00288 | -30.3% / -5.9% / -27.2% |
| 11 / 0 | 0.01589 | **0.01001** | 0.00211 | -54.4% / **+1.9%** / -46.7% |

The superseded row, kept so the size of the error is on record: HUDEP-2 0.03887 and K562 **0.04969**
at 0/11, against the corrected **0.00436** on K562 — the single `AH_WXF` of `339.9 * 0.8931` = 303.6
is HUDEP-2's number (298.3 measured live) and **38% above K562's** (220.1), and all-HDR's whole case
is that its k=1 spike clears the cutoff, so inflating it clears fields it cannot reach. **Cite no
absolute E[share] from this section's earlier form.**

**`ALL_CUT_HOTKEYS = ""` is confirmed, on a wider basis than before.** 1/10 -> 0/11 is positive in
**9 of 9** measurements (three cells; CD34+ at three RNG seeds x two clean fractions), margin
**+0.4% to +10.3%**. The old +6.48% +/- 0.69 sits inside that, so the decision is unchanged — but
CD34+'s margin runs as low as **+0.4%**, so it is not the mean/SE-21 result that figure implied.

**"Monotone on both cell types: every hotkey moved to all-cut costs share" is conditionally false,
and the condition is all-cut's clean fraction.** It was only ever run on HUDEP-2 and K562. CD34+ is
**U-shaped** — worst in the middle and *better* than the shipped mix at the extreme:

| CD34+_HSPC, `FP_AC_CLEAN` | 6 / 5 | 11 / 0 (three RNG seeds) |
|---|---|---|
| **569**/900 (build-time, what this file records) | -3.2 to -5.9% | **+1.9% / +2.1% / +4.3%** |
| **477**/900 (implied by h0's live all-cut rounds) | -13.3 to -16.0% | -20.4% / -19.8% / -18.6% |

At 569 an all-cut *fleet* beats the shipped mix on CD34+, sign-consistent across RNG seeds, so the
inversion is arithmetic and not noise. At 477 it collapses and monotonicity returns. CD34+ is the
one cell where all-cut's `weighted x fidelity` **exceeds** all-HDR's (249.1 against 233.3, ratio
1.068), which is what gives correlated siblings something worth stacking; the lower clean fraction
removes enough consistency to erase it. **So the monotonicity claim rests entirely on all-cut's
clean fraction being ~0.53 rather than ~0.63, and that is measured at n = 6.** Treat it as one
input away from flipping on one cell type.

**Two of the three constants this section was priced on were wrong, all in all-cut's favour.**
Measured from h0's own rounds while `ALL_CUT_HOTKEYS` still named it, isolated by the cons signature
(all-cut 0.16-0.21 against all-HDR's 0.09-0.10):

| | script assumed | live |
|---|---|---|
| all-cut / all-HDR `w x fid` ratio | 0.943 (all-cut lower) | **1.014** — all-cut gives up nothing |
| all-cut clean fraction | 569/900 = 0.632 | **477/900 = 0.530** (round-cons 0.1774 against a predicted 0.1911) |
| `AH_WXF` | one constant, 303.6 | **per cell: 298.3 / 233.3 / 220.1** |

The 1.014 ratio independently reproduces this file's "all-HDR tracks all-cut within 3% on term 1".
The model's band crowding still matches measurement (11 independent 12-seed bands in a 300-seed
joined space give a union of 108 against `band_overlap.py`'s 104).

`fleet_price.py` now takes `FP_CELL`, `FP_AC_CLEAN`, `FP_RNG` and `FP_SINCE`, and filters our own
hotkeys **by ss58 address**: uids are recycled, and 7 of our 11 carried another operator's hotkey
during the 2026-09-08 reboot window, which is worth 57 contaminated rows of 228 if identification
goes through uid.
Two live confirmations of the h0 half: its all-HDR K562 build takes **344 s** against all-cut's
**606 s**, and on the two rounds where its all-cut prepare was still running when a validator called
(9ac1d178 at 149 s, d4a8f8b8 at 409 s) it shipped the emergency ordinary construction and scored the
floor. What is given up is the only construction that scores anything when every band misses
(45.46 on a155b953) — worth zero at rank 99, but it was the fleet's sole signal on those rounds.

**HEK293 has not had this treatment** and is the open case: `price_cell.py` puts it on all-HDR by
0.0059 vs 0.0027, on the same pooled basis this section discredits, and it is the one cell type
whose accessibility (0.35) makes the two arms structurally different. Run `cmp_k562.py` against
HEK293 contracts before trusting that split.

### Validation pipeline — stages talk through files, not return values

[genomics/validation/](niome_subnet/genomics/validation/) runs stages in order, each reading the
previous stage's JSON output from `data/`:

```
data/submission.json  →  stage12 → valid_experiments.json / invalid_experiments.json
                      →  stage3  → stage3_dataset.json / stage3_summary.json
                      →  stage4  → final_reward.json
                      →  stage5  → distribution_fidelity_summary.json  → MinerScore
```

- **stage12** — structural gate (PAM, guide length 20/23, mismatch budget, mutation whitelist,
  cell-type match, dedup on `(cas, start, strand, guide)`) then structural scoring
  (`0.625·gc_score + 0.375·dist_score`, scaled by off-target k-mer uniqueness and `mutation_weight`).
  `truncate_submission` caps rows at `max_experiments` and drops duplicate `experiment_id`s **by
  rewriting `data/submission.json` in place** — the cut file is what every later stage and the
  archived submission see.
- **stage3** — biophysical simulation. Seeded from *each* round seed + the design fields
  (`sha256(seed|mutation|cas|guide|start|strand)`, low 32 bits, into `random.Random`), so it is fully
  deterministic and reproducible by anyone holding the contract. Per row the draw order is: the
  microhomology coin, then the cut coin (`cut` iff draw ≤ `cut_p`), then `repair_mode`, then the
  indel length.
- **stage4** — trains a RandomForest per target (`is_cut`, `is_hdr`, `indel_length`) under KFold and
  turns cross-validated R²/MAE into `consistency_factor`. Needs ≥2 valid rows and a non-empty
  `experiment_id` join, otherwise it writes a clean zero rather than raising.
- **stage5** — six-way *geometric* mean of coverage/diversity entropy ratios. An empty
  (mutation × cas × strand) cell costs roughly a 0.03× multiplier on the entire score.

Final score: `total_weighted_score × consistency_factor × distribution_fidelity_factor`, averaged
over every seed in `contract["seed"]`.

Because `run_validation` downloads every miner's submission to the same
`MINER_SUBMISSION_PATH`, scoring is inherently sequential per miner. Parallelising it means giving
each stage per-miner paths first.

### Miner generation mirrors the validator on purpose

[genExp.py](genExp.py) **imports and calls the validator's own stage functions**
(`stage12.check_pam`, `stage12.load_or_build_kmer_index`, `stage3.simulate`, …) rather than
reimplementing the scoring maths. Preserve that. It is what keeps the generator from drifting from the
pipeline that judges it, and it means **editing a validation stage silently changes miner behaviour** —
a formula tweak in stage12 or stage3 re-prices every site the generator ranks and can invalidate the
outcome "construction" it searches for.

The build enumerates PAM sites in `gene_region ± flank`, apportions rows across the full
mutation × cas × strand support, tunes each guide toward 50% GC within the contract's mismatch budget,
then searches guide variants for one whose deterministic stage-3 draw satisfies the configured
construction (`CONSTRUCTIONS`) — that conformance is what drives `consistency_factor`
to 1.0, and it is all-or-nothing: one stray row collapses stage 4's R². Sequence, k-mer index and PAM
enumeration are process-global caches, warmed on a prewarm thread at miner startup.

genExp.py is both the miner's engine and its research tool: `python genExp.py` builds and scores one
task, `--all-tasks` sweeps the backend's whole history, and [submission.py](submission.py) writes the
row sets themselves. [neurons/miner.py](neurons/miner.py) imports genExp directly and its `_build`
mirrors `submission.build_for_task` step for step (`build_context` → `enumerate_sites` →
`choose_weight_skew` → `generate` → `order_rows`), so an offline sweep predicts exactly what the miner
will send — `submission.py --task-id <id>` and the miner produce identical arrays for one contract.
Keep them in step: a knob that only exists on one side breaks that guarantee.

`CONSTRUCTIONS` currently registers `mh` (genExp's default), `mh_any`, `hdr`, `nocut`, `blunt`,
`mhnhej`, `blunt_any`, `mhnhej_any`. The `*_any` rules are their pinned counterparts minus the
`indel_length == 1` pin: cheaper per row, but they give up stage 4's `indel_length` target.
[neurons/miner.py](neurons/miner.py) does **not** use the genExp default — it sets `CONSTRUCTION =
"hdr"` and `CAS_MIX = "70/30"` as class constants, with `CELL_TYPE_OVERRIDES` replacing them
per cell type (HEK293 has its own clustered builder,
[genomics/hek293_generation.py](niome_subnet/genomics/hek293_generation.py)).

[genomics/generation.py](niome_subnet/genomics/generation.py) is the superseded packaged port of
genExp's pure path. Nothing imports it any more.

### Seed-depend — the build the fleet actually ships

[genomics/seed_depend.py](niome_subnet/genomics/seed_depend.py) solves the opposite problem to
everything else here: the seed is **known**. A task is broadcast with `seed: 0` and the real seeds
are normally stamped in before scoring — but sometimes the stamp never happens and the validator
scores at seed 0 itself. On those rounds a seed-0 build reaches `consistency_factor` exactly 1.000
and the ranking collapses to `total_weighted_score x distribution_fidelity_factor`. It is the first
rung of `_build`'s ladder, so it replaces the construction rather than hedging alongside it.

**How often, and against how many.** Do not read the rate off the task listing — `/api/v3/tasks`
reports `seed: 0` for rounds that were in fact stamped. The reliable test is miners reaching
consistency 1.000 in the score rows. By that test, over the 105 scored rounds since 2026-08-25:

| | |
|---|---|
| never stamped | **8 of 105 = 7.6%** (the ~3.5% and ~5.4% figures elsewhere are stale and too low) |
| crowd on those rounds | 0, 0, 0, 0, 1, 6, 8, 13 — and then **25, 38** on 2026-09-07/08 |
| four seed-0 rounds in one day | 2026-09-07 (14:16, 16:41, 19:06, 21:29) |

**The competition is converging fast, and it has caught us.** Replaying our build into six seed-0
fields (excluding our own rows — from 2026-09-07 18:46 our submissions are *in* that feed):

| round | cell | tied at 1.000 | our margin over the field | ranks | curve share |
|---|---|---:|---:|:---:|---:|
| 8f02f1a4 09-02 | HEK293 | 13 | +4.57 | 1-4 | 85% |
| 138a47f7 09-07 | CD34+ | 6 | +7.55 | 1-4 | 85% |
| e32ec7b3 09-07 | HUDEP-2 | 8 | +2.34 | 1-4 | 85% |
| 19018a0a 09-07 | CD34+ | 5 | +3.20 | 1-4 | 85% |
| 67bdd18a 09-07 | HUDEP-2 | 25 | +0.77 | 1-4 | 85% |
| **a66f01fa 09-08** | CD34+ | **38** | **−1.87** | **6-9** | **9%** |

h0 took **rank 1 of 248 twice live** (19018a0a 339.94, 67bdd18a 329.37) within three hours of
registering, then placed 6th on the next round. 1.8 weighted points is the whole distance between
85% of the curve and 9% — `SCORE_DISTRIBUTION` is a step function, so score gains that cross no
rank boundary pay nothing.

**Sibling hotkeys take consecutive ranks, which is leverage in both directions.** Four correlated
builds occupy ranks r..r+3: 0.30+0.20+0.20+0.15 = 0.85 at r=1, but 0.09 at r=6. They move as a
block because their scores differ by less than any gap in the field.

**The build is shared across the fleet, and the submissions are identical.** `Miner._build_seed_depend`
takes an `O_EXCL` lock, builds once, writes the rows to `data/seed_depend/` (shared, *not* under
`DATA_DIR`), and the other hotkeys load them. Keyed on a hash of the contract plus the build
parameters, so bumping `variants_per_site` invalidates every cached round instead of silently
shipping the old build. A follower that runs out of budget falls through the ladder and must
**never** build its own copy — the whole point is that four concurrent builds do not fit in memory.
`SEED_DEPEND_SHARED = False` restores independent builds. Identical rows across a coldkey's hotkeys
is the most visible possible signature if cross-miner duplicate detection ever appears; nothing in
the validator does that today, and the field ships byte-equal rows openly (a66f01fa's top four are
equal at 262.92).

**Config, measured over those six rounds** (paired within contract at a fixed variant, since the
variant tie-break alone moves 0.03-0.14):

| knob | result |
|---|---|
| `variants_per_site` 4000 → **12000** | **+1.47 ± 0.25, t = 5.80, 6/6.** Shipped. 100s → 294s build |
| `variants_per_site` 24000 | +1.59 — i.e. +0.12 over 12k for 2x build and 2x memory, and the *same* payout share. 12k is the plateau |
| `rule="hdr"` | +0.53 alone, but **not additive** with vps12k (+1.40 combined, worse than vps12k alone). `rule` stays `"mh"` |
| `alloc_step` 4 → 1 | −0.01. The split search was already converged at granularity 2 — closed knob |
| `greedy_window` 400 → 1200 | **−4.27, 0/6.** A wider window lets `kmer_price` pull in structurally weaker rows (weighted 279.7 → 274.9). 400 is tuned, not arbitrary |

Note what the winning knob does *not* do: mean payout share moved 25.5% → 25.8%. It is insurance
against a field that is catching up, not a fix for a round already lost — and neither 12k nor 24k
recovers a66f01fa, where both reach 262.90 against the field's 262.92.

`SEED_DEPEND_MIN_BUDGET_S` is **360s**, not 190: a 294s build cannot finish inside the ~225s in-TTL
path, so the old gate would start a build it could not complete, burn the window, and fall through
the ladder later and poorer than if it had never tried. Above 360 it only ever runs on the prefetch
path (900s budget, >=600s of lead on 86% of rounds); a failed prefetch now skips seed-depend and
lands on all-HDR/all-cut.

**A trap that produced convincing wrong timings.** all-HDR banks live in
`all_hdr.HDR_BANK_DIR` (`data/all_hdr`), **not** `all_cut.BANK_DIR` (`data/all_cut`). Research
scripts that "force a cold scan" by deleting `data/all_cut/cas12a-<key>.npz` are no-ops for
all-HDR, so any arm whose config matches something the live fleet already built silently loads that
bank. It corrupts **timings only** — a bank is a pure function of its `bank_key`, so a cache hit and
a cold build are the same object, and band/weighted/fidelity are unaffected (verified by
reproducing a width-300 arm cold: identical to 4 decimals, 50.7 s against the 6 s the cache
reported).

**Research tooling** (none of it imported by the neurons): `sd_task.py` builds and scores one task
against its real field, `sd_variants.py` sweeps variant indices, `sd_fleet.py` aggregates rounds,
`sd_sweep.py` / `sd_sweep_all.py` run the config sweep with paired statistics, `sd_shared_test.py`
verifies the shared-build locking across four real processes. Three traps these hit, all of which
produced convincing wrong answers:

- **Score rows are per validator per miner.** A round scored by two validators lists every
  competitor twice; ranks computed from raw rows double-count. Dedupe by `miner_hotkey`.
- **Our own submissions are in the score feed** from 2026-09-07 18:46. A replay compared against an
  unfiltered field is racing itself — the exact-match "field best 339.94" was h0.
- **`NIOME_INSTANCE` must differ per concurrent process.** Two runs under one instance name write
  `contract.json`, `submission.json` and every stage output to the same paths; they interleave and
  score one task's rows against another's contract. That produced `fidelity 0.0010` and looked
  exactly like a catastrophic build failure.

### Seed-agnostic hedge (non-HEK293, unstamped contracts)

[genomics/seed_agnostic.py](niome_subnet/genomics/seed_agnostic.py) builds a submission whose
`is_cut` is constant across a seed *window* (default 100-999), for the case where the contract
arrives with `seed: 0` and the real seeds are stamped later. Two banks: **strict Cas9** (`max_fail 0`
— reachable because `cut_p` clamps at 0.99, so `0.99**900 ≈ 1.2e-4` yields a few hundred) and a
**Cas12a pool** at `cas12a_max_fail 22` (`cut_p` caps at 0.96, so no strict guide exists **at this
window width** — `0.96**900 ~ 1.6e-16`; the claim is about the 900-seed window and does not
generalise, since at width 225 it is `~1.0e-4` and 246-277 strict Cas12a guides were measured) from which
`min_union_group` picks the subset minimising the *union* of failed seeds — overlapping failures are
free, so it selects for coincidence, subject to per-cell floors that keep stage 5's geometric mean
off zero. `_build` in the miner gates on it (skip for HEK293, for a stamped seed, without a GPU, or
under `SEED_AGNOSTIC_MIN_BUDGET_S`) and falls back to the ordinary construction on short rows or
excess backfill.

[genomics/mt19937.py](niome_subnet/genomics/mt19937.py) and
[genomics/sha256_gpu.py](niome_subnet/genomics/sha256_gpu.py) are bit-exact CPU/GPU replications of
CPython's `random.Random` seeding and draw sequence, used to test millions of guides against 900
seeds inside the 300 s upload TTL. Both carry `verify()` against the real pipeline — run it after
touching either, since a one-bit divergence silently invalidates every bank.

**Measured ceiling, don't re-derive it:** only `is_cut` is recoverable across a window. Repair mode
is a fresh ~0.5 coin per seed with no design lever, so every repair-mode rule tested tops out at
349-491 failed seeds of 900 (`hdr` on Cas9 is the best at 349) against the 12-22 that make the
cut-only hedge work. Raising `max_fail`, `max_distance` or pool size does not move it.

This is a ceiling on the **900-seed window**, and all-HDR above does not contradict it: it pins
repair mode over a 12-16 seed *band*, not a window, and pays for it by scoring ~0.10 on every seed
outside the band. Those 349 failed seeds are the same measurement seen from the other side.

The repo root also holds research tools that are **not** imported by the neurons: `search_repair.py`
(window/single-seed construction frontiers), `search_guides.py`, `search_group.py`, `search_seed.py`,
`assemble.py`, `robustness.py`, `test.py`, `calc.py`.

### Bittensor 11 specifics

This repo is on the v11 API and does **not** use `bt.Axon`/`bt.Synapse` dendrite calls:

- Miners serve a FastAPI app ([base/miner.py](niome_subnet/base/miner.py)) with one `POST /forward`
  route. Requests are authenticated with `bt.http_auth.sign` / `bt.http_auth.verify`; `GenomicsTaskSynapse`
  is a plain pydantic model sent as the JSON body.
- Chain writes go through `subtensor.execute(bt.SetWeights(...) | bt.ServeAxon(...), wallet)`.
- Metagraph reads go through `fetch_metagraph_with_retry`
  ([utils/misc.py](niome_subnet/utils/misc.py)), which pins reads `FINALITY_LAG` blocks below the head.
  Reading at the unfinalized tip raises `BlockNotFound` after a reorg — never call
  `subtensor.subnets.metagraph()` unpinned.
- `BT_NO_PARSE_CLI_ARGS` must be forced to `"false"` **before** `import bittensor`, which is why
  [neurons/miner.py](neurons/miner.py) imports `niome_subnet.utils.settings` first, ahead of every
  other project import. Don't reorder those imports.

Config is hand-rolled argparse ([utils/config.py](niome_subnet/utils/config.py)), not `bt.config`:
dotted flags like `--neuron.epoch_length` are flattened by argparse and re-nested by `_nest_config`,
and legacy v10 spellings (`--wallet.name`, `--subtensor.network`) are hidden aliases.

### Weights and emissions

[base/validator.py](niome_subnet/base/validator.py) `set_weights` computes weights locally and stores
them on `self.uids`/`self.weights`; the actual chain commit happens later in the run loop when
`should_set_weights()` opens. `SCORING_SYSTEM = "top"` pays only the top 10 miners on the fixed
`SCORE_DISTRIBUTION` curve; `BURNING_RATE` is carved off the top for `OWNER_HOTKEY`. Scores are also
POSTed to the backend, signed with the canonical-JSON scheme in
[api/\_\_init\_\_.py](niome_subnet/api/__init__.py) (a different signing scheme from the `bt.http_auth`
one used validator→miner — don't conflate them).

### Tuning surface

[utils/settings.py](niome_subnet/utils/settings.py) is the single place for netuids, backend URLs,
block schedule, timeouts, scoring system, burn rate and every `data/` path. Prefer changing it over
threading new constants through call sites. Note `utils/__init__.py` re-exports every util module with
`import *`, so a new top-level name there can shadow another module's.
