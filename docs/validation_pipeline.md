# NIOME Validation Pipeline — Reference

How a validator turns one miner's `submission.json` into a single `final_score`, stage by stage, with every parameter that enters the arithmetic.

Code: [niome_subnet/genomics/validation/](../niome_subnet/genomics/validation/)
Constants: [niome_subnet/utils/settings.py](../niome_subnet/utils/settings.py)
Orchestration: [niome_subnet/validator/forward.py](../niome_subnet/validator/forward.py)

---

## 0. At a glance

```
submission.json (miner)
      │
      ├─ truncate_submission ......... cap rows, drop duplicate/blank experiment_id
      │
   ┌──▼─────────────────────────────────────────────┐
   │ Stage 1  structural gate      (pass / fail)    │  run ONCE — seed independent
   │ Stage 2  structural scoring   (weighted_score) │
   └──┬─────────────────────────────────────────────┘
      │ valid_experiments.json / invalid_experiments.json
      │
   ┌──▼─────────────────────────────────────────────┐
   │ Stage 3  biophysical simulation                │  re-run PER SEED
   │ Stage 4  cross-consistency ML (consistency_f)  │
   │ Stage 5  distributional fidelity (fidelity_f)  │
   └──┬─────────────────────────────────────────────┘
      │
      ▼
final_score = total_weighted_score × consistency_factor × distribution_fidelity_factor
             (averaged over all round seeds)
```

The three factors are **multiplicative**, so each is a veto:

| Factor | Range | Answers |
|---|---|---|
| `total_weighted_score` | 0 … `max_experiments × max(mutation_weight)` | *How many good designs did you submit?* |
| `consistency_factor` | 0 … 1 | *Are the simulated outcomes learnable from your features?* |
| `distribution_fidelity_factor` | 0 … 1 | *Is the design space covered evenly and non-redundantly?* |

---

## 1. Files the pipeline reads and writes

All paths come from [settings.py](../niome_subnet/utils/settings.py) and are relative to the validator's working directory.

| Constant | Path | Direction | Written by |
|---|---|---|---|
| `CONTRACT_PATH` | `data/contract.json` | in | backend (presigned URL, `fetch_task`) |
| `HBB_REFERENCE_PATH` | `data/hbb_reference.json` | in | backend (presigned URL) |
| `CHR11_PATH` | `data/chr11.fa` | in | static reference genome |
| `MINER_SUBMISSION_PATH` | `data/submission.json` | in / **rewritten** | miner via S3; rewritten by `truncate_submission` |
| `KMER_CACHE_DIR` | `data/kmer_cache/` | cache | stage 1–2 |
| `VALID_EXPERIMENTS_PATH` | `data/valid_experiments.json` | out | stage 1–2 |
| `INVALID_EXPERIMENTS_PATH` | `data/invalid_experiments.json` | out | stage 1–2 |
| `STAGE3_DATASET` | `data/stage3_dataset.json` | out | stage 3 |
| `STAGE3_SUMMARY_PATH` | `data/stage3_summary.json` | out (diagnostic) | stage 3 |
| `FINAL_REWARD_PATH` | `data/final_reward.json` | out | stage 4 |
| `DISTRIBUTION_FIDELITY_PATH` | `data/distribution_fidelity_summary.json` | out | stage 5 |

Because every stage communicates through files, the pipeline is **not** safe to run concurrently for two miners in the same working directory.

### The contract

```json
{
  "active_mutations": ["NC_000011.10:g.5226784G>C", "NC_000011.10:g.5225906G>T"],
  "cell_type": "HEK293",
  "mutation_regions": { "NC_000011.10:g.5226784G>C": null, ... },
  "mutation_weights": { "NC_000011.10:g.5226784G>C": 1.5, "NC_000011.10:g.5225906G>T": 0.7 },
  "rules": {
    "base_padding": 400,
    "cas_systems": ["Cas9", "Cas12a"],
    "max_experiments": 250,
    "max_mismatches": 3,
    "proximity_gate": false
  },
  "seed": 628,
  "version": "v1"
}
```

| Contract field | Consumed by | Effect |
|---|---|---|
| `active_mutations` | stage 1, stage 5 | Whitelist for `mutation`; also the support set for mutation-coverage entropy. |
| `cell_type` | stage 1, stage 2 | Every row's `cell_type` must match exactly (skipped only if the field is `null`). Also selects the accessibility multiplier. |
| `mutation_regions` | stage 2 | Maps a mutation to a gene region name → `region_energy_offset`. `null` ⇒ offset 0. |
| `mutation_weights` | stage 2, 3, 4 | Multiplies the structural score, and becomes the RF `sample_weight`. Default `1.0` for a mutation not listed. |
| `rules.base_padding` | stage 1, 2 | Distance scale (bp): decay constant of `dist_score`, the `consistency` cutoff, and the `proximity_gate` radius. |
| `rules.cas_systems` | stage 5 | Support set for Cas coverage entropy and for the Cas-shift diagnostic. Defaults to `["Cas9","Cas12a"]`. |
| `rules.max_experiments` | truncation | Hard row cap. `null` ⇒ uncapped. |
| `rules.max_mismatches` | stage 1 | Max Hamming distance between guide and genomic target. |
| `rules.proximity_gate` | stage 1 | When `true`, a guide farther than `base_padding` from the mutation is rejected outright rather than merely scored down. |
| `seed` | stages 3, 4 | Round seed(s). Parsed by `_parse_seeds` as a **comma-joined string** (`"122,321,431"`), so a plain int is a single-seed round. |

### The HBB reference

| Field | Used for |
|---|---|
| `mutation_map` | mutation ID → absolute chr11 coordinate. Drives every distance calculation. |
| `gene_region.start` / `.end` | Centre of the off-target k-mer index window (HBB: 5,225,464–5,227,064, i.e. 1,600 bp). |
| `chromosome`, `window_id`, `challenge` | Metadata; `challenge` mirrors the contract and is not read by the pipeline. |

### A submission row

```json
{
  "experiment_id": "exp-00000",
  "guideRNA": "GTAGACCATCAGCTGCCTAA",
  "target_alignment_start": 5226783,
  "target_alignment_end": 5226803,
  "strand": "+",
  "mutation": "NC_000011.10:g.5226784G>C",
  "cas_system": "Cas9",
  "cell_type": "HEK293"
}
```

Coordinates are **absolute 0-based offsets into the chr11 FASTA sequence**, not gene-relative.

---

## 2. Orchestration — `benchmark_submission`

[`validation/__init__.py`](../niome_subnet/genomics/validation/__init__.py)

```python
def benchmark_submission(cell_types: dict, uid: int) -> MinerScore
```

| Parameter | Type | Meaning |
|---|---|---|
| `cell_types` | `dict` | Fetched from `CELL_TYPES_URL` (`/api/v3/data/cell-types`). Shape `{"HEK293": {"accessibility": 0.35, ...}, ...}`. Only `accessibility` is read. |
| `uid` | `int` | Miner UID; a label on the returned `MinerScore`, it does not affect scoring. |

Flow:

1. Read `contract["seed"]` and split it into a list of round seeds.
2. `run_stage12(cell_types)` — **once**. Stages 1–2 are deterministic and seed independent.
3. For each seed: `run_stage3(seed)` → `run_stage4(seed)` → `run_stage5()`, collecting the stage-5 dict.
4. Arithmetic-mean every metric across seeds; `n_valid_experiments` is rounded to an int.

Multi-seed averaging exists so a submission cannot be tuned to one lucky RNG draw. Note that stage 4's `fold_seed` is also the round seed, so **each seed reshuffles the cross-validation folds as well as the simulation**.

### Feeds

| Source | What it supplies |
|---|---|
| `data/contract.json` → `seed` | The list of round seeds to iterate |
| `fetch_cell_types()` | `cell_types` dict, passed straight through to stage 2 |
| stage 5 return value, once per seed | The six metrics that get averaged |

### Emits — `MinerScore`

Defined in [`genomics/model.py`](../niome_subnet/genomics/model.py). Returned in memory only; `forward.run_validation` collects one per miner.

| Field | Type | Meaning |
|---|---|---|
| `uid` | `int` | Miner UID, echoed from the argument. |
| `final_score` | `float` | Mean of stage 5's `final_score` across seeds. **This is the only number `set_weights` reads.** |
| `log` | `str` | Always `""` in the current code — reserved for a per-miner diagnostic string. |
| `breakdown` | `dict` | Six averaged metrics, below. Forwarded verbatim to the backend and shown in the miner dashboard. |

| `breakdown` key | Type | Source | Meaning |
|---|---|---|---|
| `n_valid_experiments` | `int` | stage 4 (`len(df)`) | Rows that passed the gate **and** survived the stage 3 ⨝ stage 1–2 join. Averaged then rounded — a fractional mean means some seeds dropped rows. |
| `total_weighted_score` | `float` | stage 4 | Sum of `weighted_score` over merged rows. Seed-independent in practice, so the average equals any single seed's value. |
| `consistency_score` | `float` | stage 4 | The 0–100 form: `(0.7·max(R²,0) + 0.3·(1−nMAE)) × 100`. |
| `consistency_factor` | `float` | stage 4 | `consistency_score / 100`, clipped to `[0, 1]`. The multiplier actually used. |
| `distribution_fidelity_score` | `float` | stage 5 | Geometric mean of the six coverage ratios, before clipping. |
| `distribution_fidelity_factor` | `float` | stage 5 | The same value clipped to `[0, 1]`. |

`score`/`factor` pairs differ only by scale and clipping — the `_score` is the diagnostic, the `_factor` is what multiplies.

A `MinerScore` is later widened into a `MinerScoreDto` for the backend, adding `task_id`, `hotkey` (looked up from the metagraph) and `weight` (the miner's final normalised on-chain weight).

---

## 3. Pre-stage — `truncate_submission`

[`stage12.py:214`](../niome_subnet/genomics/validation/stage12.py#L214)

Applied before stage 1 sees anything, and it **rewrites `submission.json` in place** so that the archived copy uploaded back to S3 matches what was scored.

| Rule | Parameter | Behaviour on violation |
|---|---|---|
| Row cap | `rules.max_experiments` (250) | Iteration stops; the tail is silently dropped. |
| `experiment_id` present | — | Non-string or whitespace-only id ⇒ row dropped (not even counted as invalid). |
| `experiment_id` unique | — | Second and later occurrences dropped. |

The uniqueness rule matters more than it looks: stage 4 inner-joins stage 3 to stage 1–2 **on `experiment_id`**, so a repeated id would become a cross product — 250 rows sharing 50 ids would report 1,250 valid experiments and 5× the weighted score.

### Feeds and emits

| | File | Note |
|---|---|---|
| feed | `data/submission.json` | Raw miner upload, downloaded from S3 by `run_validation`. |
| feed | `data/contract.json` → `rules.max_experiments` | The row cap. |
| emit | `data/submission.json` | **Rewritten in place, only if something was cut.** Same list, same field names, fewer elements. |
| emit | return value | The kept list, handed straight to the stage 1 loop. |

### Field reference — a submission row

Every field is miner-supplied. Nothing else is permitted to influence scoring; a row carrying predicted efficiencies or repair outcomes is not rejected for it, but those extra keys are ignored everywhere.

| Field | Type | Required by | Meaning and constraints |
|---|---|---|---|
| `experiment_id` | `str` | truncation, stage 4 join | Non-empty after `strip()`, unique within the submission. The **only** key linking stage 1–2 output to stage 3 output. |
| `guideRNA` | `str` | stages 1, 2, 3, 5 | The spacer, written 5'→3' **on its own strand**, so a `-` strand guide is compared against `revcomp` of the genome slice. Length must be 20 or 23. |
| `target_alignment_start` | `int` | stages 1, 2, 3 | Absolute 0-based offset into `chr11.fa` of the first base of the protospacer on the `+` strand. Drives PAM slicing, the mismatch check and `distance`. |
| `target_alignment_end` | `int` | stage 1 | Must equal `target_alignment_start + len(guideRNA)` exactly. Redundant by construction — it exists as a self-consistency check. |
| `strand` | `str` | stages 1, 3, 5 | `"+"` or `"-"`. Selects which side of the protospacer the PAM sits on, and is a coverage dimension in stage 5. |
| `mutation` | `str` | stages 1, 2, 3, 5 | HGVS id; must be in `active_mutations`. Selects `mutation_weight`, `mutation_regions` and the `mutation_map` coordinate. |
| `cas_system` | `str` | stages 1, 2, 3, 5 | `"Cas9"` or `"Cas12a"`. Selects the PAM motif, which end of the guide is the off-target seed, and both stage 3 base rates. |
| `cell_type` | `str` | stage 1 | Must equal the contract's `cell_type`. Note the *accessibility* used downstream is looked up from the **contract's** cell type, not this field — so this is purely a gate. |

---

## 4. Stages 1 & 2 — `run_stage12`

[`stage12.py:252`](../niome_subnet/genomics/validation/stage12.py#L252)

```python
def run_stage12(cell_types: dict, offtarget_flank: int = 50000) -> tuple[list, list]
```

| Parameter | Default | Meaning |
|---|---|---|
| `cell_types` | — | Accessibility lookup (see above). |
| `offtarget_flank` | `50000` | bp added either side of `gene_region` to build the off-target k-mer index. Window = `[start − 50000, end + 50000]` clipped to the chromosome — 101,600 bp for HBB. Widening it makes off-target hits more likely and scores strictly lower. |

The k-mer index is built with `k = 12` (hard-coded at the call site), cached to `data/kmer_cache/kmer_<md5(seq+k)>.pkl`, and keyed on the sliced sequence — so it is reused across miners and across rounds as long as the gene region and flank are unchanged.

### Feeds

| File | Fields read | Used for |
|---|---|---|
| `data/submission.json` | the eight row fields above | Everything. |
| `data/contract.json` | `active_mutations`, `cell_type`, `mutation_weights`, `mutation_regions`, `rules.{max_mismatches, base_padding, proximity_gate, max_experiments}` | Gate thresholds and score weights. |
| `data/hbb_reference.json` | `mutation_map`, `gene_region.start`, `gene_region.end` | Distance calculations; the k-mer index window. |
| `data/chr11.fa` | first FASTA record only | The genomic truth: PAM slices, mismatch comparison, bounds. |
| `cell_types` argument | `[contract.cell_type]["accessibility"]` | Recorded for stage 3; missing cell type falls back to `1.0`. |
| `data/kmer_cache/*.pkl` | — | Reuses a prior 12-mer index if the window hash matches. |

### Stage 1 — structural gate

[`stage12.py:111`](../niome_subnet/genomics/validation/stage12.py#L111) — returns `1.0` (pass) or `0.0` plus a reason string. Checks run in order and **short-circuit on the first failure**:

| # | Check | Parameter | Reason string |
|---|---|---|---|
| 1 | `mutation ∈ active_mutations` | `contract.active_mutations` | `mutation_not_allowed` |
| 2 | `exp.cell_type == contract.cell_type` (skipped if contract value is `null`) | `contract.cell_type` | `cell_type_mismatch` |
| 3 | `len(guideRNA) ∈ {20, 23}` | fixed | `invalid_length` |
| 4 | `target_alignment_end == target_alignment_start + len(guideRNA)` | — | `invalid_alignment_end` |
| 5 | `0 ≤ start` and `start + len(guide) < len(chr11)` | — | `out_of_bounds` |
| 6 | PAM present at the right offset | `cas_system`, `strand` | `pam_ok` / `pam_out_of_bounds` / `pam_invalid_strand` / `pam_invalid_cas` |
| 7 | `hamming(guide, target) ≤ max_mismatches` | `rules.max_mismatches` (3) | `too_many_mismatches` |
| 8 | `\|start − mutation_pos\| ≤ base_padding`, only if gate enabled | `rules.proximity_gate`, `rules.base_padding` | `mutation_too_far` |
| 9 | `(cas_system, start, strand, guideRNA)` not seen before | — | `duplicate_experiment` |

Note on #7: for `strand == "-"` the comparison is `hamming(reverse_complement(guide), target)` — the guide is written 5'→3' on its own strand, the genome slice is read on the `+` strand.

Note on #9: the dedup key is the *biological* design, not the `experiment_id`. Two rows with different ids but the same guide at the same locus and strand — the second one is invalid.

#### PAM rules — `check_pam` ([`stage12.py:78`](../niome_subnet/genomics/validation/stage12.py#L78))

With `L = len(guideRNA)` and `s = target_alignment_start`:

| Cas | Strand | PAM slice | Motif required |
|---|---|---|---|
| Cas9 | `+` | `seq[s+L : s+L+3]` (3'-adjacent) | `NGG` — implemented as `pam[1:] == "GG"` |
| Cas9 | `−` | `revcomp(seq[s−3 : s])` | `NGG` |
| Cas12a | `+` | `seq[s−4 : s]` (5'-adjacent) | `TTTV` — implemented as `pam[:3] == "TTT"` |
| Cas12a | `−` | `revcomp(seq[s+L : s+L+4])` | `TTT…` |

Any other `cas_system` value fails with `invalid_cas`; any strand other than `+`/`−` fails with `invalid_strand`.

Every failing row is written to `invalid_experiments.json` as `{"experiment": …, "stage1_pass": false, "reason": …}` — that file is the first thing to read when a submission scores zero.

### Stage 2 — structural scoring

[`stage12.py:170`](../niome_subnet/genomics/validation/stage12.py#L170) — runs only on rows that passed stage 1.

| Quantity | Formula | Parameters | Notes |
|---|---|---|---|
| `gc` | fraction of `G`/`C` in the guide | — | Raw feature, reused in stage 3 energy. |
| `distance` | `\|target_alignment_start − mutation_map[mutation]\|` | `mutation_map` | Absolute bp on chr11. |
| `gc_score` | `max(0, 1 − \|gc − 0.5\| × 2)` | — | Triangular peak: 1.0 at 50 % GC, 0.0 at 0 % or 100 %. |
| `dist_score` | `exp(−distance / base_padding)` | `rules.base_padding` (400) | 1.0 on the mutation, 0.37 at 400 bp, 0.08 at 1 kb. |
| `consistency` | `1.0 if distance < base_padding else 0.3` | `rules.base_padding` | **Recorded as a stage-4 model feature only**; it does not enter the stage-2 score. |
| `base_structural_score` | `0.625 × gc_score + 0.375 × dist_score` | fixed 0.625 / 0.375 | GC is weighted ~1.67× more than proximity. |
| `offtarget_factor` | step function of seed-region hits | see below | Multiplicative penalty. |
| `structural_score` | `base_structural_score × offtarget_factor` | — | The value stored as `stage2.structural_score`. |
| `mutation_weight` | `contract.mutation_weights[mutation]`, default `1.0` | contract | 1.5 vs 0.7 in the current contract. |
| **`weighted_score`** | `structural_score × mutation_weight` | — | **The per-row currency of the whole pipeline** — stage 4 sums exactly this column. |
| `cell_type_accessibility` | `cell_types[contract.cell_type]["accessibility"]`, default `1.0` | cell-types API | HEK293 = **0.35**. Used only by stage 3. |
| `region_energy_offset` | lookup on `mutation_regions[mutation]` | table below | Used only by stage 3. |

#### `offtarget_uniqueness` ([`stage12.py:157`](../niome_subnet/genomics/validation/stage12.py#L157))

The seed region is the part of the guide nearest the PAM: **last 12 nt for Cas9**, **first 12 nt for Cas12a**. Its occurrence count in the ±50 kb index gives:

| Hits of the 12-nt seed | `offtarget_factor` |
|---|---|
| 0 | 1.0 |
| 1 – 5 | 0.7 |
| 6 – 20 | 0.4 |
| > 20 | 0.1 |

The on-target site itself sits inside the index window, so a perfectly-matching guide normally registers ≥ 1 hit — 1.0 is reachable mainly by guides carrying mismatches in the seed.

#### `REGION_ENERGY_OFFSETS`

| Region | Offset |
|---|---|
| `5UTR_or_upstream` | +0.05 |
| `exon1`, `exon2`, `exon3` | +0.03 |
| `intron1`, `intron2` | −0.03 |
| `3UTR` | +0.02 |
| unknown / `null` | 0.00 |

The table is duplicated verbatim in [stage3.py:15](../niome_subnet/genomics/validation/stage3.py#L15); stage 3 in fact reads the offset that stage 2 already stored, so the copy there is dead weight — keep the two in sync if either is edited.

### Emits — `valid_experiments.json`

A JSON list, one object per surviving row, in submission order. This is the pipeline's central artifact: stage 3 and stage 5 read it directly, and stage 4 reads it a second time to recover the weighted score.

```json
{
  "experiment": { ...the eight submitted fields, verbatim... },
  "features": {
    "gc": 0.5,
    "distance_to_mutation": 1,
    "gc_score": 1.0,
    "dist_score": 0.9975031223974601,
    "consistency": 1.0,
    "offtarget_factor": 1.0,
    "mutation_weight": 1.5,
    "cell_type": "HEK293",
    "cell_type_accessibility": 0.35,
    "mutation_region": null,
    "region_energy_offset": 0.0
  },
  "stage1": { "valid": true },
  "stage2": { "structural_score": 0.9990636708990475, "weighted_score": 1.4985955063485712 }
}
```

| Field | Type | Meaning | Read downstream by |
|---|---|---|---|
| `experiment` | `object` | The submitted row, unmodified — the audit trail back to what the miner actually sent. | stage 3 (RNG key, `cas`, `mutation`), stage 4 (`experiment_id`, `guideRNA`, `start`), stage 5 (`mutation`, `cas_system`, `strand`, `guideRNA`) |
| `features.gc` | `float` | GC fraction of the guide, 0–1. | stages 3, 4 |
| `features.distance_to_mutation` | `int` | Absolute bp between `target_alignment_start` and the mutation coordinate. **Renamed to `distance`** once stage 3 copies it. | stages 3, 4 |
| `features.gc_score` | `float` | Triangular GC fitness, 0–1. Deterministic function of `gc`. | stages 3, 4 |
| `features.dist_score` | `float` | `exp(−distance / base_padding)`, 0–1. | stages 3, 4 |
| `features.consistency` | `float` | `1.0` or `0.3` — a proximity flag, not a score. Carried purely to become a model feature. | stages 3, 4 |
| `features.offtarget_factor` | `float` | 1.0 / 0.7 / 0.4 / 0.1. **Already folded into `structural_score`; dropped at the stage 3 boundary and never seen by the model.** | nothing — diagnostic only |
| `features.mutation_weight` | `float` | The contract weight for this mutation. | stages 3, 4 (`sample_weight`) |
| `features.cell_type` | `str` \| `null` | Echo of the contract's cell type, for auditing. | nothing — diagnostic only |
| `features.cell_type_accessibility` | `float` | Chromatin accessibility multiplier, 0–1. The single largest lever on stage 3 energy. | stage 3 |
| `features.mutation_region` | `str` \| `null` | Region name from `mutation_regions`. | nothing — diagnostic only |
| `features.region_energy_offset` | `float` | −0.03 … +0.05, resolved from the region name. | stage 3 |
| `stage1.valid` | `bool` | Always `true` in this file — present so a row is self-describing when read in isolation. | nothing |
| `stage2.structural_score` | `float` | `(0.625·gc_score + 0.375·dist_score) × offtarget_factor`, 0–1. | stage 4 (as `stage2_score`, carried but **not** a model feature) |
| `stage2.weighted_score` | `float` | `structural_score × mutation_weight`. | stage 4 — summed into `total_weighted_score` |

### Emits — `invalid_experiments.json`

Same list shape, one object per rejected row. Empty when every row passes.

| Field | Type | Meaning |
|---|---|---|
| `experiment` | `object` | The rejected row, verbatim. |
| `stage1_pass` | `bool` | Always `false`. |
| `reason` | `str` | The **first** failing check — one of `mutation_not_allowed`, `cell_type_mismatch`, `invalid_length`, `invalid_alignment_end`, `out_of_bounds`, `pam_ok`, `pam_out_of_bounds`, `pam_invalid_strand`, `pam_invalid_cas`, `too_many_mismatches`, `mutation_too_far`, `duplicate_experiment`. |

Two quirks worth knowing. First, `reason` is `"pam_ok"` when the PAM lookup succeeded but the motif did not match — the `ok` refers to the bounds check inside `check_pam`, not to the row. Second, rows dropped by `truncate_submission` (over the cap, blank or duplicate `experiment_id`) never reach stage 1, so they appear in **neither** output file; the only evidence is that `submission.json` got shorter.

---

## 5. Stage 3 — biophysical simulation

[`stage3.py:172`](../niome_subnet/genomics/validation/stage3.py#L172)

```python
def run_stage3(seed=None) -> tuple[list, dict]
```

| Parameter | Default | Meaning |
|---|---|---|
| `seed` | `None` → `contract["seed"]` | The round seed. Mixed into every per-experiment RNG. `benchmark_submission` always passes it explicitly. |

This is the anti-gaming core: the miner supplies only the *design*, and all biology is rolled here, by the validator, under a seed the miner did not know when building the submission.

### Feeds

| Source | Fields read | Used for |
|---|---|---|
| `data/valid_experiments.json` | `experiment.{experiment_id, cas_system, mutation, guideRNA, target_alignment_start, strand}` | RNG key and record identity |
| | `features.{gc, distance_to_mutation, gc_score, dist_score, consistency, cell_type_accessibility, mutation_weight, region_energy_offset}` | `extract_features` — the eight fields that survive into the simulation |
| `data/contract.json` | `seed` | Fallback only; `benchmark_submission` passes the seed explicitly |

`extract_features` is a deliberate narrowing: `offtarget_factor`, `cell_type` and `mutation_region` are **not** copied, and `distance_to_mutation` is renamed to `distance`. Anything the model is meant to see has to survive this function.

### Per-experiment RNG — `experiment_seed`

```
sha256("round_seed|mutation|cas_system|guideRNA|target_alignment_start|strand") mod 2^32
```

Each row gets its own `random.Random`, so outcomes are independent of row order and fully reproducible. Two identical designs would draw identically — which is why stage 1 rejects duplicates.

### The simulation chain

| Step | Formula | Constants |
|---|---|---|
| `sequence_energy` | `clip(accessibility × (1.8·gc + 0.6·exp(−distance/1500) + region_offset), 0, 1)` | 1.8 GC coefficient, 0.6 proximity coefficient, **1500 bp** decay (note: *not* `base_padding`) |
| `microhomology_trigger` | `rng.random() < min(0.6, 2.2 · gc·(1−gc))` | cap 0.6, slope 2.2 → peaks at 0.55 for gc = 0.5 |
| `cut_probability` | `clip(base + 0.18 × energy, 0.4, 0.99)` | `base` = **0.86** (Cas9) / **0.78** (Cas12a) |
| cut draw | `rng.random() > cut_p` ⇒ `outcome = "no_cut"`, `indel_length = 0` | — |
| `repair_mode` | sample proportional to unnormalised weights: `hdr = hdr_base + 0.35·energy`; `mh_nhej = 0.30 if mh else 0.12`; `blunt = 0.35` | `hdr_base` = **0.32** (Cas9) / **0.24** (Cas12a) |
| `sample_indel_length` | `HDR → 0`; `MH_NHEJ → max(1, int(gamma(2.2, 2.8)))`; `BLUNT_NHEJ → max(1, int(expo(0.6)))` | mean ≈ 6.2 and ≈ 1.7 bp |

The accessibility multiplier is the dominant term. With HEK293 at 0.35, energy is capped at 0.84 even for a 100 %-GC guide, and a typical 50 %-GC on-target guide lands at **0.525** — so `cut_p` ≈ 0.954 for Cas9 and 0.874 for Cas12a. **Cut rates are high and nearly flat**, which is exactly why the `is_cut` target is hard for stage 4 to learn (see below).

RNG call order is fixed and load-bearing: `microhomology_trigger` consumes one draw *before* the cut test, so changing the order of those two lines changes every outcome.

### Emits — `stage3_dataset.json`

A JSON list, one record per valid experiment, in `valid_experiments.json` order.

```json
{
  "experiment_id": "exp-00000",
  "mutation": "NC_000011.10:g.5226784G>C",
  "mutation_weight": 1.5,
  "cas": "Cas9",
  "outcome": "BLUNT_NHEJ",
  "indel_length": 2,
  "features": {
    "gc": 0.5, "distance": 1, "gc_score": 1.0, "dist_score": 0.9975031223974601,
    "consistency": 1.0, "cell_type_accessibility": 0.35,
    "mutation_weight": 1.5, "region_energy_offset": 0.0
  },
  "energy": 0.5248600466562979,
  "mh": true
}
```

| Field | Type | Meaning | Read downstream by |
|---|---|---|---|
| `experiment_id` | `str` | Join key back to stage 1–2. | stage 4 (merge key) |
| `mutation` | `str` | Echoed for grouping. | stage 3 summary; carried into stage 4's frame but unused by the model |
| `mutation_weight` | `float` | Echoed to the top level so grouping does not have to reach into `features`. | stage 3 summary |
| `cas` | `str` | **Renamed** from `cas_system`. Watch this when writing code that spans the two files. | stage 4 (`cas_system` column), stage 5 diagnostic |
| `outcome` | `str` | `no_cut` \| `HDR` \| `MH_NHEJ` \| `BLUNT_NHEJ`. The primary simulated label. | stage 4 (`is_cut`, `is_hdr`), stage 5 diagnostic |
| `indel_length` | `int` | 0 for `no_cut` and `HDR`; ≥ 1 otherwise. | stage 4 (regression target), stage 5 diagnostic |
| `features` | `object` | The eight surviving design features, copied unchanged from stage 2. | stage 4 — seven of them become `X` |
| `energy` | `float` | `sequence_energy`, 0–1. Computed here, and **fed back in as a model feature** — the one column in `X` the miner does not control directly. | stage 4 |
| `mh` | `bool` | Whether the microhomology draw fired. Cast to `int` for the model. | stage 4 |

Note that `no_cut` records are kept, not filtered. They carry `indel_length = 0` and are what makes `is_cut` a non-constant target at all — currently 15 of 250.

### Emits — `stage3_summary.json`

**Diagnostics only. Nothing downstream reads this file** — it exists for operators debugging why a submission scored the way it did.

| Field | Type | Meaning |
|---|---|---|
| `n` | `int` | Number of simulated records = number of valid experiments. |
| `cut_rate` | `float` | `1 − no_cut/n`. Currently 0.94 — see the stage 4 note on why a high, flat cut rate suppresses `consistency_factor`. |
| `mean_indel_length` | `float` | Mean over **all** records including the zeros from `HDR` and `no_cut`, so it understates the size of actual indels. |
| `mean_energy` | `float` | Mean `sequence_energy`. Tracks accessibility and GC. |
| `outcomes` | `object` | Raw counts per outcome label, e.g. `{"BLUNT_NHEJ": 84, "HDR": 110, "MH_NHEJ": 41, "no_cut": 15}`. |
| `mutation_weight_breakdown` | `object` | Per mutation, sorted by `mutation_weight` descending: `mutation_weight`, `n`, `cut_rate`, `mean_energy`, `mean_indel_length`. The quickest read on whether the submission is skewed toward the high-weight mutation. |
| `weight_correlations.weight_vs_cut` | `float` | Pearson r between `mutation_weight` and the 0/1 cut indicator. |
| `weight_correlations.weight_vs_energy` | `float` | Pearson r between `mutation_weight` and `energy`. |
| `weight_correlations.weight_vs_indel_length` | `float` | Pearson r between `mutation_weight` and `indel_length`. |

The three correlations are a bias check on the simulator, not on the miner: they should sit near zero, because mutation weight is an economic parameter and must not leak into simulated biology. The current run reads −0.129 / −0.050 / +0.004, which is the expected noise floor for two weight levels at n = 250.

---

## 6. Stage 4 — cross-consistency

[`stage4.py:118`](../niome_subnet/genomics/validation/stage4.py#L118)

```python
def run_stage4(n_folds: int = 5, seed=None) -> dict
```

| Parameter | Default | Meaning |
|---|---|---|
| `n_folds` | `5` | Requested KFold splits. Effective splits are `max(min(n_folds, n), 2)`, so a tiny dataset degrades to 2 folds rather than crashing. |
| `seed` | `None` → `contract["seed"]` | Becomes `fold_seed`, the `random_state` of a `shuffle=True` KFold. **The round seed therefore decides the fold partition**, so the same rows in a different order can score differently. |

The question this stage asks: *given only the design features, can a model predict the simulated outcomes?* A submission of internally coherent designs yields a learnable feature→outcome mapping; noise does not.

### Feeds

| Source | Fields read |
|---|---|
| `data/stage3_dataset.json` | `experiment_id`, `mutation`, `cas`, `outcome`, `indel_length`, `energy`, `mh`, and `features.{gc, distance, gc_score, dist_score, consistency}` |
| `data/valid_experiments.json` | `experiment.{experiment_id, mutation, cas_system, guideRNA, target_alignment_start}`, `features.{gc, distance_to_mutation, gc_score, dist_score, consistency, mutation_weight}`, `stage2.{structural_score, weighted_score}` |
| `data/contract.json` | `seed` — fallback only |

Note both files are flattened with overlapping column names, then **`stage12_slim` deliberately subsets to six columns** — `experiment_id`, `guideRNA`, `start`, `stage2_score`, `mutation_weight`, `weighted_score` — so the merge produces no `_x`/`_y` suffixes. The consequence: `gc`, `distance`, `gc_score`, `dist_score` and `consistency` in the merged frame all come from the **stage 3 copy**, and `flatten_stage12` computing them a second time is wasted work.

### Assembly

1. Flatten `stage3_dataset.json` and `valid_experiments.json` into DataFrames.
2. **Guard A** — `len(stage12) < 2 or len(stage3) < 2` ⇒ write an all-zero result and return. (Empty frames have no columns to subset; one row cannot be split into folds.)
3. Inner-join on `experiment_id`, keeping `stage2_score`, `mutation_weight`, `weighted_score`, `guideRNA`, `start` from the stage 1–2 side.
4. **Guard B** — empty merge ⇒ same all-zero result.

Both guards score a clean zero rather than raising, so a malformed submission is penalised, not crashed — and the validator keeps going to the next miner.

#### Merged frame — column roles

| Column | Origin | Role |
|---|---|---|
| `experiment_id` | both | Join key. Not a feature. |
| `gc`, `distance`, `gc_score`, `dist_score`, `consistency` | stage 3 copy | **Features 1–5** |
| `energy`, `mh` | stage 3 | **Features 6–7** — `mh` cast to `int` |
| `outcome` | stage 3 | Source of targets `is_cut` and `is_hdr`; itself dropped |
| `indel_length` | stage 3 | **Target 3**, used raw |
| `mutation_weight` | stage 1–2 | `sample_weight` for fitting *and* scoring |
| `weighted_score` | stage 1–2 | Summed into `total_weighted_score`. Never a feature. |
| `mutation`, `cas_system`, `guideRNA`, `start`, `stage2_score` | carried | Present in the frame, read by nothing. Available if you extend `build_X`. |

### Model

| Setting | Value | Note |
|---|---|---|
| Estimator | `RandomForestRegressor` | Regression even for the two binary targets. |
| `n_estimators` | 200 | |
| `max_depth` | 12 | |
| `random_state` | **42, hard-coded** | Independent of the round seed — only the fold split varies per seed. |
| `sample_weight` | `df["mutation_weight"]` | Applied in `fit` **and** in `r2_score` / `mean_absolute_error`, so high-weight mutations dominate the metric too. |

**Feature matrix `X`** (7 columns, order fixed by `build_X`): `gc`, `distance`, `gc_score`, `dist_score`, `consistency`, `energy`, `mh`. A missing column raises `ValueError`.

Two of these are near-degenerate in practice: `gc_score` is a deterministic function of `gc`, and `consistency` is constant 1.0 whenever every guide sits within `base_padding`.

**Targets `y`** (each fitted separately):

| Target | Definition |
|---|---|
| `is_cut` | `1` if `outcome != "no_cut"` |
| `is_hdr` | `1` if `outcome == "HDR"` |
| `indel_length` | integer bp |

### Metrics

Per target, `evaluate` returns `r2_mean`, `r2_std`, `mae_mean`, `mae_std`, `residual_std_mean`, `n_folds` averaged over folds. Then:

```
avg_r2   = mean(r2_mean          for the 3 targets)
avg_nmae = mean(mae_mean / std(y_full)  for the 3 targets)   # normalized_mae, raw MAE if std < 1e-9

consistency_score  = (0.7 × max(avg_r2, 0) + 0.3 × (1 − avg_nmae)) × 100
consistency_factor = clip(consistency_score / 100, 0, 1)      # NaN → 0

total_weighted_score = df["weighted_score"].sum()
final_reward         = total_weighted_score × consistency_factor
```

`final_reward` in `final_reward.json` is stage 4's own view of the score; the authoritative number is stage 5's `final_score`, which multiplies in fidelity as well.

### Emits — `final_reward.json`

A single object. **Stage 5 reads this file rather than a return value**, so it must be fresh for the current seed.

| Field | Type | Meaning | Read by |
|---|---|---|---|
| `n_valid_experiments` | `int` | `len(df)` after the merge — rows that passed the gate *and* joined. Differs from the stage 1 pass count only when ids fail to match. | stage 5 (echoed), `MinerScore` |
| `total_weighted_score` | `float` | `df["weighted_score"].sum()`. Unbounded above; the scale term of the final score. | stage 5 (multiplicand) |
| `consistency_score` | `float` | The 0–100 blend. Can exceed 100 or go negative in principle; `NaN` is coerced to `0.0` on write. | `MinerScore` |
| `consistency_factor` | `float` | `consistency_score / 100` clipped to `[0, 1]`. | stage 5 (multiplicand) |
| `final_reward` | `float` | `total_weighted_score × consistency_factor`. **Superseded** — stage 5 recomputes with fidelity and never reads this. | nothing |
| `model_results` | `object` | Per-target cross-validation detail, below. Empty `{}` when a guard fired. | nothing — diagnostic |

`model_results` has one entry per target (`is_cut`, `is_hdr`, `indel_length`):

| Field | Type | Meaning |
|---|---|---|
| `r2_mean` | `float` | Mean weighted R² across folds. **Routinely negative** — the model predicts worse than the weighted mean. Clamped at 0 before it reaches the score. |
| `r2_std` | `float` | Spread of R² across folds. Large values mean the fold partition is doing the work, not the features. |
| `mae_mean` | `float` | Mean weighted absolute error. In target units — probability for the two binary targets, bp for `indel_length`. |
| `mae_std` | `float` | Fold-to-fold spread of MAE. |
| `residual_std_mean` | `float` | Mean unweighted standard deviation of `y_test − pred`. Reported only; not part of the score. |
| `n_folds` | `int` | Folds actually run — 5 normally, 2 for a very small submission. |

Only `r2_mean` and `mae_mean` feed the score. The rest are there to tell you *why* a `consistency_factor` came out where it did: compare `mae_mean` against `std(y)` for that target to see the normalised error the 0.3 term actually consumes.

An all-zero output from either guard omits nothing — it writes the same six keys with zeros and `model_results: {}`, so downstream readers never have to special-case it.

**Practical behaviour.** `max(avg_r2, 0)` clamps the R² term, and with cut rates pinned near 0.95 the outcome targets are close to constant, so `avg_r2` is routinely *negative* and the 70 % term contributes **nothing**. What survives is `0.3 × (1 − avg_nmae)` — a normalised-error term, and `is_cut`'s tiny standard deviation inflates its normalised MAE. Empirically this pins `consistency_factor` low (≈0.08 in the current `final_reward.json`) with a practical ceiling around 0.39. Treat it as a nearly-fixed haircut, not a lever: it is far cheaper to move `total_weighted_score` and `distribution_fidelity_factor`.

---

## 7. Stage 5 — distributional fidelity

[`stage5.py:188`](../niome_subnet/genomics/validation/stage5.py#L188)

```python
def run_stage5(k: int = 12) -> dict
```

| Parameter | Default | Meaning |
|---|---|---|
| `k` | `12` | k-mer length for guide-diversity entropy. Coincidentally equal to the off-target seed length but a separate knob. |

The question here: *is the submission a genuine sweep of the design space, or 250 near-copies of one good guide?* This is the anti-monoculture stage.

### Feeds

| Source | Fields read | Used for |
|---|---|---|
| `data/valid_experiments.json` | `experiment.{mutation, cas_system, strand, guideRNA}` | All six coverage ratios. Nothing from `features` or `stage2` is touched. |
| `data/stage3_dataset.json` | `cas`, `outcome`, `indel_length` | The Cas-shift diagnostic only. |
| `data/contract.json` | `active_mutations`, `rules.cas_systems` | Declared support sets for the entropy denominators. |
| `data/final_reward.json` | `total_weighted_score`, `consistency_factor`, `n_valid_experiments` | Multiplied into `final_score` and echoed. |

The support sets come from the **contract, not the data**. A mutation nobody targeted still counts as an empty bucket and drags the entropy ratio down — that is the mechanism, not a bug.

### The six ratios

All coverage ratios use `coverage_entropy_ratio` = Shannon entropy (base 2) of observed counts over the **full declared support**, divided by `log2(|support|)`. Perfectly uniform ⇒ 1.0; everything in one bucket ⇒ 0.0. A support of size ≤ 1 returns 1.0 by convention.

| Ratio | Support | Size (current contract) |
|---|---|---|
| `mutation_coverage_entropy_ratio` | `contract.active_mutations` | 2 |
| `cas_system_coverage_entropy_ratio` | `rules.cas_systems` | 2 |
| `strand_coverage_entropy_ratio` | `["+", "−"]` | 2 |
| `joint_coverage_entropy_ratio` | mutation × cas × strand | 8 |
| `kmer_diversity_entropy_ratio` | all `k`-mers pooled across every guide | `H(pool) / log2(total_kmers)` |
| `distinct_guide_ratio` | — | `len(set(guides)) / len(guides)` |

`kmer_diversity_entropy_ratio` normalises by `log2(total)` where `total` is the count of k-mer *instances*, not distinct k-mers — so the theoretical 1.0 requires every k-mer across the whole submission to be unique. It returns 0.0 when `total ≤ 1`.

```
distribution_fidelity_score  = geometric_mean(the six ratios)     # values clipped up to 1e-9
distribution_fidelity_factor = clip(score, 0, 1)
```

The **geometric** mean is the point: one collapsed dimension drags the whole factor down in a way an arithmetic mean would hide. Submitting only Cas9, or only the high-weight mutation, costs far more than the missing rows would suggest.

Note the tension with `mutation_weights`: piling rows onto the 1.5-weight mutation raises `total_weighted_score` but lowers `mutation_coverage_entropy_ratio`. The current run's 198/52 split gives a mutation ratio of 0.738 — that imbalance is a deliberate trade, and the optimum is not the uniform split.

### Diagnostics (computed, written, **not scored**)

`cas_specific_shift_diagnostic` compares the first two Cas systems that have at least `MIN_PER_CAS = 5` records:

- `repair_mode_jensen_shannon_divergence` — JSD (base 2, range 0…1) between the two repair-outcome distributions.
- `indel_length_wasserstein_distance` — mean absolute gap between the two indel-length distributions across `n_quantiles = 200` evenly spaced quantiles.

Fewer than two qualifying systems ⇒ `{"insufficient_data": true}`.

### Emits — `distribution_fidelity_summary.json`

| Field | Type | Meaning |
|---|---|---|
| `n_valid_experiments` | `int` | `len(stage12_valid)` — counted here from the file, *not* from stage 4's merge, so the two can disagree if ids failed to join. |
| `mutation_coverage_entropy_ratio` | `float` | 0–1. Uniformity across `active_mutations`. |
| `cas_system_coverage_entropy_ratio` | `float` | 0–1. Uniformity across `rules.cas_systems`. |
| `strand_coverage_entropy_ratio` | `float` | 0–1. Uniformity across `+` / `−`. |
| `joint_coverage_entropy_ratio` | `float` | 0–1 over the full mutation × cas × strand grid — catches the case where each margin looks balanced but the combinations do not. |
| `kmer_diversity_entropy_ratio` | `float` | 0–1. Sequence-level redundancy across all guides pooled. |
| `distinct_guide_ratio` | `float` | 0–1. Exact-duplicate guide strings only. Stage 1 has already rejected rows matching on `(cas, start, strand, guide)`, so this drops below 1.0 only when the same sequence appears at a different locus, strand, or Cas system. |
| `distribution_fidelity_score` | `float` | Geometric mean of the six above. |
| `cas_specific_shift_diagnostic` | `object` | See below. Not scored. |
| `coverage_detail` | `object` | Raw `mutation_counts`, `cas_system_counts`, `strand_counts` — the fastest way to see which dimension collapsed. |

When there are no valid experiments the file shrinks to three keys — `n_valid_experiments: 0`, `distribution_fidelity_score: 0.0`, and `note: "no valid experiments"` — so a reader must treat the ratio fields as optional.

`cas_specific_shift_diagnostic` fields:

| Field | Type | Meaning |
|---|---|---|
| `insufficient_data` | `bool` | `true` when fewer than two Cas systems have ≥ 5 records; the remaining fields are then absent. |
| `compared` | `array` | The two Cas systems actually compared, in `rules.cas_systems` order. |
| `repair_mode_jensen_shannon_divergence` | `float` | 0–1. Near 0 means the two nucleases produce indistinguishable repair-mode mixes. |
| `indel_length_wasserstein_distance` | `float` | In bp. Mean quantile gap between the two indel-length distributions. |

Both are sanity checks on the simulator: the two Cas systems have different base rates, so a value pinned at exactly 0 would suggest the `cas` branch is not firing.

### Emits — the return value

Stage 5 returns a dict rather than writing it; `benchmark_submission` averages these across seeds.

| Field | Source | Meaning |
|---|---|---|
| `n_valid_experiments` | stage 4 file | Merge-survivor count, passed through. |
| `total_weighted_score` | stage 4 file | Passed through. |
| `consistency_score` | stage 4 file | Passed through. |
| `consistency_factor` | stage 4 file | Passed through, and used as a multiplicand. |
| `distribution_fidelity_score` | computed here | Pre-clipping geometric mean. |
| `distribution_fidelity_factor` | computed here | Clipped to `[0, 1]`. |
| `final_score` | computed here | The product of the three factors. |

### The final number

```python
final_score = stage4["total_weighted_score"] × stage4["consistency_factor"] × distribution_fidelity_factor
```

Stage 5 reads `total_weighted_score` and `consistency_factor` back out of `final_reward.json`, so **stage 4 must have run for the same seed immediately before**.

---

## 8. From `final_score` to on-chain weight

[`base/validator.py:214`](../niome_subnet/base/validator.py#L214), [`utils/weight_utils.py`](../niome_subnet/utils/weight_utils.py)

| Constant | Value | Effect |
|---|---|---|
| `SCORING_SYSTEM` | `"top"` | `"top"` = winner-takes-most ladder; `"linear"` = weight proportional to score. |
| `TOP_MINER_COUNT` | 10 | Only the top 10 positive scorers receive any weight under `"top"`. |
| `SCORE_DISTRIBUTION` | `[0.3, 0.2, 0.2, 0.15, 0.05, 0.03, 0.025, 0.02, 0.015, 0.01]` | Rank ladder, renormalised if fewer than 10 miners score positive. Rank matters, margin does not. |
| `BURNING_RATE` | 0.02 | Miners split 98 % of weight; `OWNER_HOTKEY` receives 2 %. At `1.0` the validator burns everything and skips validation. |

After the ladder, `process_weights_for_netuid` applies the chain's `max_weight_limit` and `min_allowed_weights`, then the vector is renormalised and emitted.

Round timing ([forward.py:188](../niome_subnet/validator/forward.py#L188)), on `blocks = (current_block − BASE_BLOCK_NUMBER) mod INTERVAL_BLOCKS` with `INTERVAL_BLOCKS = 720`:

| Block offset | Phase |
|---|---|
| 0 – 599 | `broadcast_task` — fetch contract, hand each miner a presigned S3 PUT URL (`SUBMISSION_TIMEOUT` 300 s) |
| 600 – 699 | `run_validation` — download each submission, `benchmark_submission`, `set_weights`, upload the top `FINAL_SUBMISSION_COUNT = 5` submissions to the backend |
| 700 – 719 | idle |

One submission per IP **and** per coldkey — a coldkey fanning out many hotkeys behind one machine does not get extra attempts.

---

## 9. Worked example

Row `exp-00000` from the current `data/` snapshot: 20-nt Cas9 guide, `+` strand, 50 % GC, 1 bp from the 1.5-weight mutation, in HEK293 (accessibility 0.35).

**Stage 2**

```
gc_score   = 1 − |0.50 − 0.5|·2                = 1.0
dist_score = exp(−1 / 400)                     = 0.99750
base       = 0.625·1.0 + 0.375·0.99750         = 0.99906
offtarget_factor                               = 1.0
structural_score                               = 0.99906
weighted_score = 0.99906 × 1.5                 = 1.49860   ✓ matches valid_experiments.json
```

**Stage 3**

```
energy = 0.35 · (1.8·0.5 + 0.6·exp(−1/1500))   = 0.52486
cut_p  = 0.86 + 0.18 · 0.52486                 = 0.95447
p_mh   = min(0.6, 2.2 · 0.5 · 0.5)             = 0.55
repair weights (mh drawn true): hdr 0.5037 | mh_nhej 0.30 | blunt 0.35
  ⇒ P(HDR | cut, mh) = 0.5037 / 1.1537         = 0.437
```

**Whole submission (250 rows, single seed)**

```
total_weighted_score          = 280.841
consistency_factor            =   0.0806
distribution_fidelity_factor  =   0.9315
final_score = 280.841 × 0.0806 × 0.9315 = 21.10
```

The fidelity factor is nearly free here (0.93); the consistency factor is throwing away 92 % of the weighted score. That ratio is typical.

---

## 10. Field propagation map

Where each per-row field is born, and where it dies. Reading down a column tells you what that stage can actually see.

| Field | submission | valid_experiments | stage3_dataset | stage 4 `X`/`y` | stage 5 |
|---|:--:|:--:|:--:|:--:|:--:|
| `experiment_id` | ● born | ● | ● | join key | — |
| `guideRNA` | ● born | ● | — | carried, unused | ● k-mer + distinct ratios |
| `target_alignment_start` | ● born | ● | — | carried as `start`, unused | — |
| `target_alignment_end` | ● born | ● | — | — | — |
| `strand` | ● born | ● | — | — | ● coverage |
| `mutation` | ● born | ● | ● | carried, unused | ● coverage |
| `cas_system` | ● born | ● | ● as `cas` | carried, unused | ● coverage + diagnostic |
| `cell_type` | ● born | ● | — | — | — |
| `gc` | | ● born | ● | **feature** | — |
| `distance_to_mutation` | | ● born | ● as `distance` | **feature** | — |
| `gc_score` | | ● born | ● | **feature** | — |
| `dist_score` | | ● born | ● | **feature** | — |
| `consistency` | | ● born | ● | **feature** | — |
| `offtarget_factor` | | ● born | ✗ dropped | — | — |
| `cell_type_accessibility` | | ● born | ● | — | — |
| `mutation_region` | | ● born | ✗ dropped | — | — |
| `region_energy_offset` | | ● born | ● | — | — |
| `mutation_weight` | | ● born | ● | `sample_weight` | — |
| `structural_score` | | ● born | — | carried as `stage2_score`, unused | — |
| `weighted_score` | | ● born | — | **summed** | — |
| `energy` | | | ● born | **feature** | — |
| `mh` | | | ● born | **feature** | — |
| `outcome` | | | ● born | **targets** `is_cut`, `is_hdr` | ● diagnostic |
| `indel_length` | | | ● born | **target** | ● diagnostic |

Three things fall out of this table:

- **Only seven columns ever reach the model.** Five of them (`gc`, `distance`, `gc_score`, `dist_score`, `consistency`) are pure functions of the miner's design, and two (`energy`, `mh`) are validator-computed. `gc_score` is redundant with `gc`, and `consistency` is constant whenever every guide sits inside `base_padding` — so the effective feature space is smaller than it looks.
- **`offtarget_factor` and `mutation_region` are terminal.** They shape `structural_score` and `region_energy_offset` at birth and are then discarded, so the model cannot learn around them.
- **`weighted_score` never meets the model.** Scale (stage 4's sum) and learnability (stage 4's R²/MAE) are computed from disjoint columns, which is what makes the three final factors genuinely independent axes.

---

## 11. Parameter quick reference

**Contract-driven** (changes round to round)

| Parameter | Current | Stage |
|---|---|---|
| `active_mutations` | 2 HBB variants | 1, 5 |
| `cell_type` | `HEK293` (accessibility 0.35) | 1, 2, 3 |
| `mutation_weights` | 1.5 / 0.7 | 2, 3, 4 |
| `mutation_regions` | both `null` | 2, 3 |
| `rules.base_padding` | 400 | 1, 2 |
| `rules.max_mismatches` | 3 | 1 |
| `rules.max_experiments` | 250 | truncation |
| `rules.proximity_gate` | `false` | 1 |
| `rules.cas_systems` | `Cas9`, `Cas12a` | 5 |
| `seed` | 628 | 3, 4 |

**Hard-coded** (constant unless the code changes)

| Parameter | Value | Location |
|---|---|---|
| Guide lengths accepted | 20 or 23 nt | `stage1` |
| Off-target flank | 50,000 bp | `run_stage12` |
| k-mer size (off-target + fidelity) | 12 | `run_stage12`, `run_stage5` |
| Off-target factor steps | 1.0 / 0.7 / 0.4 / 0.1 at 0 / ≤5 / ≤20 / >20 hits | `offtarget_uniqueness` |
| Structural blend | 0.625 GC + 0.375 distance | `stage2` |
| Energy coefficients | 1.8 GC, 0.6 proximity, 1500 bp decay | `sequence_energy` |
| Microhomology | slope 2.2, cap 0.6 | `microhomology_trigger` |
| Cut base rate | 0.86 Cas9 / 0.78 Cas12a, +0.18·energy, clipped [0.4, 0.99] | `cut_probability` |
| HDR base | 0.32 Cas9 / 0.24 Cas12a, +0.35·energy | `repair_mode` |
| MH-NHEJ / blunt weights | 0.30 (or 0.12) / 0.35 | `repair_mode` |
| Indel distributions | Gamma(2.2, 2.8) / Exp(0.6) | `sample_indel_length` |
| RF hyperparameters | 200 trees, depth 12, `random_state=42` | `evaluate` |
| KFold | 5 splits, shuffled, `random_state = round seed` | `evaluate` |
| Consistency blend | 0.7 R² + 0.3 (1 − nMAE), ×100 | `run_stage4` |
| Fidelity aggregation | geometric mean of 6 ratios | `compute_distribution_fidelity` |
| `MIN_PER_CAS` | 5 | `compute_distribution_fidelity` |
| Wasserstein quantiles | 200 | `wasserstein_1d` |

---

## 12. Implementation notes and edge cases

Behaviour worth knowing before changing anything here:

- **`build_kmer_index` skips the final k-mer.** The loop is `range(len(seq) - k)`, so the k-mer starting at `len(seq) - k` is never indexed. One position out of ~101,600; harmless in practice, but the cache key does not encode this, so fixing the loop requires clearing `data/kmer_cache/`.
- **Bounds check is strict.** `start + len(guide) >= len(seq)` rejects a guide that ends exactly on the last base.
- **Stage 5 depends on stage 4's file, not its return value.** Running stage 5 against a stale `final_reward.json` silently mixes seeds.
- **The row order you upload is the order that gets scored.** Stage 4 shuffles from the round seed, but the shuffle is applied to the merged frame in file order, so two permutations of the same 250 rows can produce different fold partitions and different `consistency_score`s. Benchmark the exact array you submit.
- **`truncate_submission` mutates the input file.** Re-running the pipeline on an already-truncated `submission.json` is a no-op, but the original oversized submission is gone.
- **`consistency` and `region_energy_offset` are computed in stage 2 and consumed in stage 3/4** — stage 2's own score ignores both.
- **`REGION_ENERGY_OFFSETS` is defined twice**, in `stage12.py` and `stage3.py`; only the stage-2 copy is ever read.
- **`_parse_seeds` expects a string.** A contract with `"seed": 628` works (`str(628).split(",")`), but a JSON list would not.
- **No concurrency.** Every stage round-trips through fixed paths in `data/`; validating two miners in parallel in one directory corrupts both.
