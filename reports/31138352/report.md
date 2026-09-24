# Submission audit: 31138352-9ccf-4fa0-907e-dd010395e27a

**HEK293**, created **2026-09-20T13:13:22.503212 UTC**. Scoring seeds: **431, 189, 566**.

**Best of h0–h9: h4, 118.836826, rank 13/248. Leader: UID 238, 171.489132.** h4 achieved 69.30% of the leader, a gap of 52.652306. Top-10 cutoff: 126.506891; h4 was 7.670065 below it.

**Seven uploads were found: h0 and h4–h9**, each with **250 valid rows**. Six uploaded prepared conjunctions with **80 Cas12a + 170 Cas9** and **eight-seed HDR bands**. **h6 uploaded an ordinary fallback with 75 Cas12a + 175 Cas9, zero clean seeds and no HDR band across 100–999.** All seven have zero reward weight in this snapshot. h1–h3 have no task build/window/upload records or official score entries; their results are unavailable, not zero.

**h4 hit HDR seed 431: all 250 rows returned HDR**, scoring **307.934233** at that seed. Its other seed scores were **24.383393** at 189 and **24.192852** at 566, giving the mean **118.836826**. This was the only uploaded clean or HDR hit. Seeds 189 and 566 were outside every uploaded clean set.

[Official task score feed](https://niome-api.genomes.io/api/v3/miners/scores?task_id=31138352-9ccf-4fa0-907e-dd010395e27a&limit=40000), fetched 2026-09-20T16:06:45.759634+00:00. 248 records, 248 miners, 1 validator. Score timestamp: 2026-09-20T15:36:36.646596 UTC. Audit generated 2026-09-20T16:09:06.830487+00:00.

## Official scores

| Miner | UID | Rank / 248 | Final score | Weighted score | Consistency | Fidelity | Valid rows | Reward weight |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Top miner | 238 | 1 | 171.489132 | 328.830612 | 0.553975 | 0.941400 | 250 | 0.294000 |
| h0 | 248 | 162 | 23.508451 | 336.321824 | 0.076240 | 0.916830 | 250 | 0.000000 |
| h1 | N/A | N/A | No score entry | N/A | N/A | N/A | N/A | N/A |
| h2 | N/A | N/A | No score entry | N/A | N/A | N/A | N/A | N/A |
| h3 | N/A | N/A | No score entry | N/A | N/A | N/A | N/A | N/A |
| h4 | 190 | 13 | 118.836826 | 336.532329 | 0.385916 | 0.915021 | 250 | 0.000000 |
| h5 | 234 | 19 | 29.421729 | 338.309375 | 0.094829 | 0.917092 | 250 | 0.000000 |
| h6 | 175 | 17 | 30.410496 | 343.014938 | 0.098137 | 0.903392 | 250 | 0.000000 |
| h7 | 65 | 166 | 22.801175 | 332.979989 | 0.074139 | 0.923616 | 250 | 0.000000 |
| h8 | 211 | 167 | 22.636215 | 338.572713 | 0.073246 | 0.912786 | 250 | 0.000000 |
| h9 | 245 | 120 | 25.484327 | 338.054636 | 0.082284 | 0.916161 | 250 | 0.000000 |

Rank is 1 plus the number of strictly higher scores; ties share rank. Final score and breakdown factors are means over three seeds. Multiplying mean factors is not generally equivalent to averaging per-seed products. Missing hotkeys have no score or rank.

Leader hotkey: `5G3bTVRofa3vib2U6idTDWXFFyMc9GPo5jxpjECxZMZtHNG4`.

## Joined windows, actual band space and cut windows

The **12:24 UTC plan** recorded width-100 windows **100 / 200 / 300**, joined as **100–399**. Actual windows were **400 / 100 / 500** for scoring seeds **431 / 189 / 566**: one matched that prediction.

For these conjunctions, the plan's predicted space did **not** determine band placement. The active layout used the full **100–999** space, with **300-seed contiguous windows at stride 100** for h0, h4, h5, h6, h7, h8 and h9. The source's `band_space()` ignores the predicted argument for these hotkeys. The logs retain the label `predicted space of 900`; that denotes this full space, not the plan's 300 predicted seeds.

For **HEK293**, `cut mode=union` resolves to each hotkey's **own 300-seed band candidate window**. The original build logs confirm `cut 300 seeds (union) | k=8`; this task did not use a 900-seed cut search. h6 prepared that configuration but uploaded an ordinary fallback instead.

| Hotkey | Actual construction | Assigned band candidates / conjunction cut window | Offset in 100–999 | Actual clean /900 | Clean inside assigned window |
|---|---|---|---:|---:|---:|
| h0 | conjunction | 100–399 | 0 | 26 | 26 |
| h1 | No task build or upload | N/A | N/A | N/A | N/A |
| h2 | No task build or upload | N/A | N/A | N/A | N/A |
| h3 | No task build or upload | N/A | N/A | N/A | N/A |
| h4 | conjunction | 200–499 | 100 | 28 | 28 |
| h5 | conjunction | 300–599 | 200 | 24 | 24 |
| h6 | ordinary fallback | 400–699 (prepared only) | 300 | 0 | 0 |
| h7 | conjunction | 500–799 | 400 | 26 | 26 |
| h8 | conjunction | 600–899 | 500 | 28 | 28 |
| h9 | conjunction | 700–999 | 600 | 26 | 26 |

The six uploaded conjunction candidate windows together still cover all **100–999**, despite h6's fallback. They overlap. All recovered full-row clean seeds lie within their corresponding conjunction cut windows; replay found no additional full-row clean seeds outside them.

## Actual uploaded HDR bands

**Clean:** all 250 uploaded rows cut at that seed. **HDR band:** all 250 uploaded rows return HDR. Both were recovered by scanning all 900 seeds (100–999), independently of the 300-seed build search.

| Hotkey | Clean /900 | Exact uploaded HDR band | Clean scoring-seed hits | HDR hits |
|---|---:|---|---|---|
| h0 | 26 | 119, 128, 139, 159, 243, 268, 367, 393 | none | none |
| h1 | N/A | No uploaded rows available | N/A | N/A |
| h2 | N/A | No uploaded rows available | N/A | N/A |
| h3 | N/A | No uploaded rows available | N/A | N/A |
| h4 | 28 | 251, 350, 351, 388, 395, 421, 431, 460 | 431 | 431 |
| h5 | 24 | 310, 373, 432, 434, 458, 460, 502, 542 | none | none |
| h6 | 0 | none | none | none |
| h7 | 26 | 502, 522, 541, 570, 641, 702, 720, 724 | none | none |
| h8 | 28 | 605, 646, 665, 709, 745, 766, 796, 802 | none | none |
| h9 | 26 | 702, 720, 763, 906, 914, 957, 961, 992 | none | none |

The uploaded fleet union contains **148/900 clean seeds** and **44 distinct HDR-band seeds** (48 band memberships before overlaps). h6's prepared band is excluded from these totals.

## Per-seed results

Each entry is **per-seed final score / HDR rows / no-cut rows**, from 250 rows per seed.

| Hotkey | Seed 431 | Seed 189 | Seed 566 | Mean final score |
|---|---:|---:|---:|---:|
| h0 | 22.478165 / 98 / 26 | 26.053131 / 106 / 17 | 21.994057 / 100 / 15 | 23.508451 |
| h4 | 307.934233 / 250 / 0 | 24.383393 / 105 / 15 | 24.192852 / 115 / 12 | 118.836826 |
| h5 | 27.506201 / 112 / 17 | 29.981201 / 102 / 15 | 30.777785 / 109 / 16 | 29.421729 |
| h6 | 30.619901 / 112 / 20 | 30.125827 / 103 / 15 | 30.485760 / 111 / 12 | 30.410496 |
| h7 | 22.801544 / 95 / 19 | 21.759910 / 110 / 21 | 23.842070 / 118 / 15 | 22.801175 |
| h8 | 24.400054 / 112 / 20 | 24.082081 / 120 / 18 | 19.426511 / 108 / 24 | 22.636215 |
| h9 | 27.543337 / 111 / 15 | 20.327154 / 96 / 23 | 28.582491 / 116 / 15 | 25.484327 |

## Comparison with the leader

h4's weighted score was **336.532329**, versus the leader's **328.830612**. Mean consistency was **0.385916 versus 0.553975**, and fidelity was **0.915021 versus 0.941400**. The largest factor gap is consistency. h4's strong score at its HDR hit was averaged with two low seed scores. The leader's submitted rows, clean sets, HDR bands and per-seed outcomes were not available in the inspected data; its aggregate score does not identify those sets.

## h6: uploaded fallback versus prepared band

h6 began preparation at **13:13:36 UTC**. The validator requested its submission at **13:15:06**. At **13:15:09**, it waited up to **177 seconds** for the prepared build. At **13:18:06**, that build was still running (270 seconds elapsed), so the miner used ordinary construction while preparation still held the GPU. It uploaded at **13:18:22.734**, with **88 seconds** of URL TTL left. The prepared conjunction finished at **13:19:43.059**, about **80 seconds after the upload**.

That unsubmitted preparation logged **26 clean seeds within 400–699** and recorded band **427, 462, 508, 511, 542, 570, 609, 632**. The actual fallback has **zero clean seeds and no HDR band across 100–999**. Its window record combines `source=all_hdr_not_attempted` / `window=null` with a band from the completed background build; the band does not describe uploaded rows. All six other uploaded bands match their task window records.

## h1–h3: no task activity

No task preparation, received-task/upload event, window record or upload archive was found for h1–h3. None appears in the complete official score feed. Their earlier runtime logs reported them unregistered on netuid 55 on 2026-09-19. No submitted clean set, band, score or rank is available for these hotkeys on this task.

## Exact uploaded clean sets

Each list covers all 100–999 seeds; band seeds are included. h6's set is empty.

### h0: 26 clean seeds

115, 119, 128, 139, 152, 159, 168, 171, 200, 205, 221, 234, 235, 243, 250, 251, 253, 258, 268, 289, 332, 341, 367, 377, 386, 393.

### h4: 28 clean seeds

218, 226, 244, 251, 263, 273, 280, 282, 288, 299, 350, 351, 356, 364, 372, 384, 388, 391, 395, 401, 421, 431, 445, 460, 468, 477, 478, 482.

### h5: 24 clean seeds

308, 310, 321, 326, 344, 365, 366, 373, 376, 382, 383, 390, 432, 434, 458, 460, 480, 492, 496, 502, 529, 542, 571, 586.

### h6: 0 clean seeds

none.

### h7: 26 clean seeds

500, 502, 522, 523, 525, 532, 541, 543, 570, 577, 604, 617, 639, 641, 670, 675, 702, 712, 713, 720, 724, 734, 745, 747, 750, 767.

### h8: 28 clean seeds

601, 604, 605, 624, 626, 646, 665, 673, 695, 699, 709, 714, 721, 745, 753, 766, 774, 779, 796, 802, 829, 853, 858, 866, 874, 888, 890, 893.

### h9: 26 clean seeds

701, 702, 713, 720, 733, 739, 763, 802, 810, 814, 819, 826, 845, 885, 888, 906, 907, 909, 914, 939, 944, 957, 961, 980, 992, 995.

## Verification and artifacts

- All seven archives identify this task; contracts and references match the published task after the seed stamp.
- All 1,750 submitted rows passed stage 1/2, with accessibility 0.35. Local three-seed scores exactly match all seven official scores.
- Reused h6/h8/h9 archived validations and generated h0/h4/h5/h7 validations in report-local directories.
- Replayed outcomes and indel lengths match local stage-3 detail for every uploaded row at all three scoring seeds.
- Full-row clean counts inside each conjunction's cut window match the original logged Cas12a clean counts.
- h7 has a Cas12a-only clean seed outside its cut window, **385**, where the full 250-row submission is not clean. It is excluded from the reported full-row clean set.
- Recovered sets over **900 × 250 × 7 = 1,575,000 row/seed combinations**.
- [Full analysis](analysis.json), [official scores](scores.json), [task](task.json), [sanitized build/upload and registration evidence](build_evidence.json), [layout source evidence](layout_evidence.json), [manifest](manifest.json), [CSV](summary.csv).
- Reproduce with `.venv/bin/python reports/31138352/validate.py`, then `.venv/bin/python reports/31138352/audit.py`, then `.venv/bin/python reports/31138352/render_report.py`. The analysis includes validator-code and submission SHA-256 hashes.

| Hotkey | Public hotkey | Upload time UTC | Submission archive |
|---|---|---|---|
| h0 | `5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW` | 2026-09-20 13:58:38,470 | [data/inst/niome_hotkey/result/2026-09-20T13:59:07](../../data/inst/niome_hotkey/result/2026-09-20T13:59:07/submission.json) |
| h1 | `5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb` | N/A | No upload archive |
| h2 | `5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU` | N/A | No upload archive |
| h3 | `5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf` | N/A | No upload archive |
| h4 | `5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2` | 2026-09-20 13:22:08,150 | [data/inst/niome_hotkey4/result/2026-09-20T13:22:52](../../data/inst/niome_hotkey4/result/2026-09-20T13:22:52/submission.json) |
| h5 | `5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn` | 2026-09-20 14:03:02,167 | [data/inst/niome_hotkey5/result/2026-09-20T14:03:30](../../data/inst/niome_hotkey5/result/2026-09-20T14:03:30/submission.json) |
| h6 | `5H5v45M2i6cFtS3Di4abh2zjWrPmuJZzNPoAYtVAXVFa6sp3` | 2026-09-20 13:18:22,734 | [data/inst/niome_hotkey6/result/2026-09-20T13:18:53](../../data/inst/niome_hotkey6/result/2026-09-20T13:18:53/submission.json) |
| h7 | `5CS8FdHr8Ddv14QHr75E5K7e9xv3wre9m66ywRHZMR8zhntN` | 2026-09-20 13:26:48,174 | [data/inst/niome_hotkey7/result/2026-09-20T13:28:29](../../data/inst/niome_hotkey7/result/2026-09-20T13:28:29/submission.json) |
| h8 | `5Cr9gJ3ukDDdxnhpyRM58gdMz7Sjuw8u58iiuYphVz1ZUGD2` | 2026-09-20 13:39:38,028 | [data/inst/niome_hotkey8/result/2026-09-20T13:40:08](../../data/inst/niome_hotkey8/result/2026-09-20T13:40:08/submission.json) |
| h9 | `5Ge6UdGJE7BabD3bv8L9AR79rW47K3XQetcTHgaLD9NF2Vf5` | 2026-09-20 14:11:46,286 | [data/inst/niome_hotkey9/result/2026-09-20T14:12:15](../../data/inst/niome_hotkey9/result/2026-09-20T14:12:15/submission.json) |
