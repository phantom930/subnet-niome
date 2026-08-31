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

### Seed-agnostic hedge (non-HEK293, unstamped contracts)

[genomics/seed_agnostic.py](niome_subnet/genomics/seed_agnostic.py) builds a submission whose
`is_cut` is constant across a seed *window* (default 100-999), for the case where the contract
arrives with `seed: 0` and the real seeds are stamped later. Two banks: **strict Cas9** (`max_fail 0`
— reachable because `cut_p` clamps at 0.99, so `0.99**900 ≈ 1.2e-4` yields a few hundred) and a
**Cas12a pool** at `cas12a_max_fail 22` (`cut_p` caps at 0.96, so no strict guide exists) from which
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
