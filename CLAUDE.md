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
| fidelity | 0.89 | 0.93 | **0.96** | 0.976 |
| **E[pay] aggregate** | 0.0490 | 0.0583 | **0.0598** | 0.0524 |

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

#### The fleet — eight hotkeys, seven on overlapping width-300 windows

**Read this before trusting anything below about windows.** As of 2026-09-09 the fleet runs
**eight** hotkeys — h0 uid 74, h1 209, h2 235, h3 196, h4 189, h5 147, h6 136, h7 75; only h8 is
unregistered — and they are on **bands**, not seed-depend (`SEED_DEPEND_VARIANTS` is empty, so every
process carries `NIOME_SEED_DEPEND=` and skips that rung).

The layout is `window_plan.FIXED_WINDOWS`, which overrides both `RANK_BY_HOTKEY` and the
concentrate/spread code and consults no prediction: h1-h7 take **overlapping width-300 windows
stepping by 100** (100-399, 200-499, … 700-999), and `ALL_CUT_HOTKEYS="niome_hotkey"` puts h0 on
all-cut for **every** cell type, so it plays no window at all. `miner.sh` needs
`ALLOW_OVERLAPPING_WINDOWS=1` for this: `assert_disjoint_windows` otherwise **exits** rather than
warning, which is the right default for every other layout.

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
`mh_any`, `shaped` at any `gc_spread`, and the naive baseline. The population median for 250-row
submissions in the 3-seed regime is **0.101**, i.e. identical to ours. **Nobody has a better floor.**

Round score averages the three seeds, so with a floor of 0.10 the only reachable values are:

| seeds in band | round `consistency_factor` | round final (weighted 230, fid 0.95) |
|---|---|---|
| 0 | 0.10 | ~22 |
| 1 | **0.40** | ~87 |
| 2 | 0.70 | ~153 |
| 3 | 1.00 | ~219 |

h0 hit one seed on task 9ed335da and scored exactly 0.397 / 86.58, confirming the model. **Never
quote the single-seed band score (~290) as a round score** — the k=1 value is what places.

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
| `group_size` away from 80 (HEK293, re-tested after the Cas9 fix) | the two trends are real and cross exactly at 80: fidelity climbs 0.820 → 0.964 across group 42 → 125 while weighted falls 206 → 164, and the **product peaks at group 80** (179.0). Priced on 31 real HEK293 fields every arm is 0.0010-0.0014 E[share] with `P(place | k=1)` a flat **48%** — quality moves k=1 only between 60 and 69, so band size decides the ordering and group 42's band 9 edges it by 8%, inside noise. Group 80 stands. |

`not_hdr`'s **0.577** is the one number worth remembering: two pinned targets lands exactly in the
range the leaders' placing rounds occupy, which is the arithmetic confirmation that a placing round
is one band hit plus two floor seeds — not a flat construction.

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
