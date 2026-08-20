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

There is no test suite, linter config, or CI. Verify changes by running a neuron (`--mock` for the
chain parts) or by exercising the scoring pipeline directly:

```python
# From the repo root, with data/contract.json, data/hbb_reference.json, data/submission.json
# and data/chr11.fa present — re-scores a submission exactly as a validator would.
from niome_subnet.genomics.validation import benchmark_submission
benchmark_submission(cell_types_dict, uid=0)
```

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

Broadcast dedups by IP **and** coldkey (`seen_ips`, `seen_coldkeys`, persisted in `state.npz`) so one
operator running many hotkeys gets one shot per task. That state resets when `task_id` changes.

### Submission path

The miner's HTTP reply is an empty ack — the dataset travels out of band as a PUT to the presigned URL
that arrived with the task, and must land before that URL's 300 s TTL (`SUBMISSION_TIMEOUT`) expires.
A miner that misses the window is indistinguishable from one that was never contacted; there is no
retry within a task id and no feedback channel, which is why
[neurons/miner.py](neurons/miner.py) scores its own build locally before uploading.

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
- **stage3** — biophysical simulation. Seeded from `contract.seed` + the design fields, so it is fully
  deterministic and reproducible by anyone holding the contract — but **the miner is not holding it**:
  see "The seed is not in the broadcast contract" below.
- **stage4** — trains a RandomForest per target (`is_cut`, `is_hdr`, `indel_length`) under KFold and
  turns cross-validated R²/MAE into `consistency_factor`. Needs ≥2 valid rows and a non-empty
  `experiment_id` join, otherwise it writes a clean zero rather than raising.
  `consistency_factor = 0.7·max(avg_r2, 0) + 0.3·(1 − avg_nmae)`, and the term that dominates it is a
  degeneracy, not a fit: when no row draws `no_cut`, `is_cut` is a constant column, every fold's
  `r2_score` hits its 0-numerator/0-denominator case and returns 1.0, and `normalized_mae`
  short-circuits on `std < 1e-9`. That single event is the whole gap between the field's ~0.10 cluster
  and its ~0.33 one. Nothing else in stage 4 is predictable from `X`, so a design's remaining job is
  to keep the forest from overfitting the two targets that are not.
- **stage5** — six-way *geometric* mean of coverage/diversity entropy ratios. An empty
  (mutation × cas × strand) cell costs roughly a 0.03× multiplier on the entire score.

Final score: `total_weighted_score × consistency_factor × distribution_fidelity_factor`.

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

Sequence, k-mer index and PAM enumeration are process-global caches, warmed on a prewarm thread at
miner startup.

### The seed is not in the broadcast contract

The `Task` a validator sends is `{id, contract_url, hbb_ref_url}`; the miner downloads the contract
from that presigned URL and it carries **`seed: 0`**. `run_validation` re-fetches the task at
validation time, and *that* copy carries the seed the backend stamped in between. The public
`/api/v3/tasks` list shows the same thing — the newest task reads `seed: 0` until its round closes.
Checked for structure and there is none to exploit: no hash of the task id or timestamp reproduces it,
no LCG fits the sequence, and the lag-92 repeat visible in the history is a single replayed block that
stopped around index 199. What *is* knowable is the support — 243 of the 261 stamped seeds land in
`[100, 999]` (`SEED_CANDIDATES`).

So `genExp.CONSTRUCTIONS` and the `pure`/`shaped` strategies, which pick guides by their stage-3 draw
under `ctx.seed`, are research tools only: on task e824bae7 the `mh` construction predicted
`consistency_factor` 1.0 and was paid **0.1048**. They are kept because they still measure the
pipeline honestly when you hand them a stamped task.

`strategy="robust"` (the default, genExp section 7b) is the production build and never reads the seed:

- **One coordinate per (mutation, cas, strand) cell.** Stage 1 dedups on `(cas, start, strand, guide)`
  — the guide is in the key — so one coordinate carries many rows, and every guide within the mismatch
  budget at an unchanged GC count is the *same row* to stages 1, 2 and 5. That makes `X` eight distinct
  feature vectors instead of ~230, and the forest can only return group means: measured r2 on
  is_hdr/indel_length goes from −0.29/−0.39 to −0.06/−0.06.
- **Guides chosen for seed survival.** A row's cut draw is a pure function of
  `(seed, mutation, cas, guide, start, strand)`, so `cut_failure_seeds` computes, for each candidate
  guide, the exact set of `SEED_CANDIDATES` it would draw `no_cut` under. Rows are then assigned to
  minimise the *union* of those sets (Cas12a cells first — `cut_p` caps at 0.96 against Cas9's 0.99 —
  then block-coordinate refinement passes). Every seed outside the union is a seed where all 250 rows
  cut and `is_cut` is constant. Measured ~460–540 of 900 survive on the high-accessibility cell types.

The scan is the build's entire cost (~50 s on 12 cores) and `variant_pool` is trimmed by core count so
the wall clock stays inside the 300 s presigned-URL TTL. Two things silently break it: a fixed
failure-count `cap` (must follow `cut_p` — HEK293's accessibility of 0.35 keeps energy under the clamp
and quadruples the failure rate), and a coordinate whose reference GC is off target, which spends
mismatch budget on GC and collapses its variant pool from thousands to ~190.

genExp.py is both the miner's engine and its research tool: `python genExp.py` builds and scores one
task, `--all-tasks` sweeps the backend's whole history, and [submission.py](submission.py) writes the
row sets themselves. [neurons/miner.py](neurons/miner.py) imports genExp directly and its `_build`
mirrors `submission.build_for_task` step for step (`build_context` → `enumerate_sites` →
`choose_weight_skew` → `generate` → `order_rows`), so an offline sweep predicts exactly what the miner
will send — `submission.py --task-id <id>` and the miner produce identical arrays for one contract.
Keep them in step: a knob that only exists on one side breaks that guarantee.

[genomics/generation.py](niome_subnet/genomics/generation.py) is the superseded packaged port of
genExp's pure path. Nothing imports it any more.

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
