# Submission audit: 14dcd43f-6f42-4d7b-8349-475662dd7a85

**HEK293**, created **2026-09-18T19:53:29.058938 UTC**. Scoring seeds: **124, 479, 237**.

**Best of h0–h9: h6, 120.796799, rank 14/248. Leader: UID 9, 174.450894.** h6 achieved 69.24% of the leader, a gap of 53.654096. Top-10 cutoff: 133.869791; h6 was 13.072993 below it.

**All ten hotkeys uploaded 250 valid rows: 80 Cas12a + 170 Cas9.** Each used a prepared conjunction with an eight-seed HDR band and a cut search over all 900 seeds (100–999). All ten have zero reward weight in this score snapshot.

**h6 hit HDR seed 124: all 250 uploaded rows returned HDR, scoring 318.719399 at that seed.** Its other two seed scores were 25.200971 at 479 and 18.470026 at 237, giving the mean 120.796799. No other hotkey hit a clean seed or an HDR band at any of the three scoring seeds.

[Official task score feed](https://niome-api.genomes.io/api/v3/miners/scores?task_id=14dcd43f-6f42-4d7b-8349-475662dd7a85&limit=40000), fetched 2026-09-18T22:20:32.080077+00:00. 248 records, 248 miners, 1 validator. Leader score timestamp: 2026-09-18T22:18:23.856534 UTC. Audit generated 2026-09-18T22:25:46.800592+00:00.

## Official scores

| Miner | UID | Rank / 248 | Final score | Weighted score | Consistency | Fidelity | Valid rows | Reward weight |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Top miner | 9 | 1 | 174.450894 | 313.518890 | 0.571458 | 0.973700 | 250 | 0.294000 |
| h0 | 248 | 174 | 23.024611 | 349.251851 | 0.072509 | 0.909202 | 250 | 0.000000 |
| h1 | 41 | 114 | 27.341189 | 349.309116 | 0.086184 | 0.908199 | 250 | 0.000000 |
| h2 | 163 | 136 | 26.838529 | 343.097950 | 0.085145 | 0.918714 | 250 | 0.000000 |
| h3 | 118 | 160 | 24.920463 | 345.032146 | 0.079138 | 0.912662 | 250 | 0.000000 |
| h4 | 190 | 20 | 30.859702 | 345.742833 | 0.097437 | 0.916043 | 250 | 0.000000 |
| h5 | 234 | 177 | 22.435874 | 341.461931 | 0.070988 | 0.925590 | 250 | 0.000000 |
| h6 | 175 | 14 | 120.796799 | 350.923018 | 0.379007 | 0.908232 | 250 | 0.000000 |
| h7 | 65 | 171 | 23.660454 | 350.977577 | 0.074271 | 0.907665 | 250 | 0.000000 |
| h8 | 211 | 163 | 24.798591 | 342.511995 | 0.078432 | 0.923120 | 250 | 0.000000 |
| h9 | 245 | 154 | 25.518931 | 343.614494 | 0.080776 | 0.919412 | 250 | 0.000000 |

Rank is 1 plus the number of strictly higher scores; ties share rank. Final score is the mean of the three per-seed scores. Breakdown factors are also means, so multiplying the displayed mean factors is not generally equivalent to averaging the per-seed products.

Leader hotkey: `5CAbCnSQgkxbBpBFTVAaJpcszhJTF4SmrQVYCDkg2yuSAFQW`.

## Joined seed windows and per-hotkey candidates

The **19:17 UTC plan** selected width-100 windows **200 / 300 / 900**, joined as **200–399 ∪ 900–999** (300 seeds). The actual seed windows were **100 / 400 / 200**, so one of the three actual windows matched the prediction: seed **237** was inside the predicted space; **124 and 479** were outside it.

**h0–h5** took 150-seed circular slices of the predicted 300 seeds, at stride 50. **h6–h9** took 150-seed slices of its 600-seed complement, **100–199 ∪ 400–899**, at stride 150. Their four candidate slices partition that complement. All ten used the same full cut search, **100–999**.

| Hotkey | Band group | Exact 150-seed candidate slice | Offset within group | Clean /900 | Clean in candidate slice |
|---|---|---|---:|---:|---:|
| h0 | predicted | 200–349 | 0 | 18 | 9 |
| h1 | predicted | 250–399 | 50 | 21 | 12 |
| h2 | predicted | 300–399 ∪ 900–949 | 100 | 20 | 8 |
| h3 | predicted | 350–399 ∪ 900–999 | 150 | 19 | 10 |
| h4 | predicted | 200–249 ∪ 900–999 | 200 | 19 | 9 |
| h5 | predicted | 200–299 ∪ 950–999 | 250 | 19 | 10 |
| h6 | complement | 100–199 ∪ 400–449 | 0 | 20 | 9 |
| h7 | complement | 450–599 | 150 | 21 | 10 |
| h8 | complement | 600–749 | 300 | 19 | 8 |
| h9 | complement | 750–899 | 450 | 18 | 10 |

The build logs confirm the predicted/complement group sizes and offsets. Their candidate-window messages sometimes show only the minimum and maximum seed; the table expands the actual joined slices and preserves gaps. For example, h6's logged `100-449 (150 seeds)` is **100–199 ∪ 400–449**, not the entire contiguous 100–449 range.

## Actual submitted HDR bands

**Clean set:** all 250 uploaded rows cut at that seed. **HDR band:** all 250 uploaded rows return HDR. Both were recovered by replay over seeds 100–999. The full-row clean sets equal the Cas12a group clean sets, and their counts match the original build logs. Every HDR band is a subset of its clean set and its candidate slice.

| Hotkey | Clean /900 | Exact uploaded HDR band (8 seeds) | Clean scoring-seed hits | HDR scoring-seed hits |
|---|---:|---|---|---|
| h0 | 18 | 205, 206, 260, 318, 325, 332, 333, 346 | none | none |
| h1 | 21 | 267, 270, 330, 332, 378, 380, 384, 393 | none | none |
| h2 | 20 | 300, 302, 334, 350, 378, 381, 903, 935 | none | none |
| h3 | 19 | 359, 362, 378, 382, 906, 920, 951, 955 | none | none |
| h4 | 19 | 208, 221, 929, 944, 945, 949, 951, 961 | none | none |
| h5 | 19 | 208, 241, 243, 252, 256, 289, 297, 951 | none | none |
| h6 | 20 | 111, 124, 136, 145, 154, 174, 426, 433 | 124 | 124 |
| h7 | 21 | 450, 458, 532, 542, 550, 552, 574, 583 | none | none |
| h8 | 19 | 676, 687, 694, 720, 727, 728, 732, 749 | none | none |
| h9 | 18 | 756, 762, 768, 829, 852, 866, 882, 892 | none | none |

The fleet covers **180 distinct clean seeds** and **74 distinct HDR-band seeds** across 100–999 (80 band memberships before overlaps). Its only clean and HDR scoring-seed hit is **124 on h6**, in the complementary group. The predicted group had no clean or HDR hit even though seed 237 lay in its candidate space. This describes this task's outcome; it does not establish a long-run advantage for either group.

## Per-seed results

Each entry is **per-seed final score / HDR rows / no-cut rows**. The total is 250 rows per seed.

| Hotkey | Seed 124 | Seed 479 | Seed 237 | Mean final score |
|---|---:|---:|---:|---:|
| h0 | 25.496365 / 109 / 18 | 22.843216 / 115 / 16 | 20.734250 / 104 / 19 | 23.024611 |
| h1 | 25.759306 / 109 / 21 | 27.458823 / 102 / 21 | 28.805439 / 107 / 10 | 27.341189 |
| h2 | 24.749757 / 110 / 16 | 30.559559 / 90 / 17 | 25.206270 / 105 / 13 | 26.838529 |
| h3 | 23.348095 / 112 / 18 | 23.095810 / 107 / 19 | 28.317484 / 114 / 9 | 24.920463 |
| h4 | 31.582131 / 105 / 13 | 24.426686 / 109 / 17 | 36.570289 / 118 / 9 | 30.859702 |
| h5 | 20.972214 / 107 / 18 | 21.217788 / 112 / 22 | 25.117620 / 116 / 12 | 22.435874 |
| h6 | 318.719399 / 250 / 0 | 25.200971 / 119 / 21 | 18.470026 / 106 / 20 | 120.796799 |
| h7 | 20.803110 / 97 / 25 | 27.218157 / 110 / 9 | 22.960095 / 100 / 19 | 23.660454 |
| h8 | 20.438291 / 102 / 18 | 27.072258 / 112 / 19 | 26.885225 / 105 / 13 | 24.798591 |
| h9 | 23.805971 / 103 / 17 | 23.734511 / 105 / 16 | 29.016310 / 101 / 14 | 25.518931 |

## Comparison with the leader

h6's weighted score is **350.923018**, versus the leader's **313.518890**. Mean consistency is **0.379007 versus 0.571458**, and fidelity is **0.908232 versus 0.973700**. h6's consistency factor was 1.0 at its HDR hit, but low scores at the other two seeds reduced its task average. The mean consistency difference is the largest factor gap in the published breakdown. The leader's submitted rows, clean set, bands and individual seed outcomes were not available in the inspected data; its aggregate score does not identify those sets.

## Exact clean sets per hotkey

These are the full submitted clean sets over 100–999, including the HDR band seeds.

### h0: 18 clean seeds

164, 205, 206, 241, 260, 318, 325, 332, 333, 346, 432, 462, 468, 530, 537, 549, 723, 826.

### h1: 21 clean seeds

207, 256, 267, 270, 330, 332, 333, 377, 378, 380, 384, 393, 398, 481, 501, 505, 609, 649, 652, 844, 888.

### h2: 20 clean seeds

107, 287, 300, 302, 334, 350, 378, 381, 485, 489, 504, 543, 564, 618, 636, 731, 812, 879, 903, 935.

### h3: 19 clean seeds

203, 206, 311, 323, 359, 362, 371, 378, 382, 402, 460, 483, 617, 745, 906, 908, 920, 951, 955.

### h4: 19 clean seeds

194, 208, 221, 541, 599, 602, 753, 759, 800, 861, 863, 885, 929, 944, 945, 949, 951, 961, 986.

### h5: 19 clean seeds

166, 186, 203, 208, 241, 243, 252, 256, 289, 297, 357, 502, 508, 534, 644, 674, 827, 951, 967.

### h6: 20 clean seeds

111, 124, 132, 136, 145, 154, 174, 240, 279, 322, 326, 361, 426, 433, 475, 486, 696, 749, 887, 995.

### h7: 21 clean seeds

222, 394, 437, 447, 450, 458, 499, 532, 542, 550, 552, 560, 574, 583, 618, 679, 708, 750, 842, 912, 990.

### h8: 19 clean seeds

118, 303, 304, 403, 421, 428, 447, 448, 577, 676, 687, 694, 720, 727, 728, 732, 749, 809, 924.

### h9: 18 clean seeds

122, 144, 417, 488, 513, 677, 702, 741, 756, 762, 768, 792, 829, 852, 865, 866, 882, 892.

## Build and upload provenance

All ten original builds began at **19:53:30–36 UTC** and finished at **19:56:21–28 UTC**. Each successful upload used its prepared 250-row conjunction. The ten recorded bands match the bands recovered from the actual uploads. No fallback or later replacement build was needed for these submissions.

| Hotkey | Upload time UTC | Submission archive |
|---|---|---|
| h0 | 2026-09-18 19:56:27,996 | [data/inst/niome_hotkey/result/2026-09-18T19:57:00](../../data/inst/niome_hotkey/result/2026-09-18T19:57:00/submission.json) |
| h1 | 2026-09-18 20:13:04,181 | [data/inst/niome_hotkey1/result/2026-09-18T20:13:33](../../data/inst/niome_hotkey1/result/2026-09-18T20:13:33/submission.json) |
| h2 | 2026-09-18 20:34:55,782 | [data/inst/niome_hotkey2/result/2026-09-18T20:35:24](../../data/inst/niome_hotkey2/result/2026-09-18T20:35:24/submission.json) |
| h3 | 2026-09-18 20:45:55,345 | [data/inst/niome_hotkey3/result/2026-09-18T20:46:24](../../data/inst/niome_hotkey3/result/2026-09-18T20:46:24/submission.json) |
| h4 | 2026-09-18 20:39:53,553 | [data/inst/niome_hotkey4/result/2026-09-18T20:40:21](../../data/inst/niome_hotkey4/result/2026-09-18T20:40:21/submission.json) |
| h5 | 2026-09-18 20:40:08,878 | [data/inst/niome_hotkey5/result/2026-09-18T20:40:36](../../data/inst/niome_hotkey5/result/2026-09-18T20:40:36/submission.json) |
| h6 | 2026-09-18 20:01:37,204 | [data/inst/niome_hotkey6/result/2026-09-18T20:02:05](../../data/inst/niome_hotkey6/result/2026-09-18T20:02:05/submission.json) |
| h7 | 2026-09-18 20:47:55,303 | [data/inst/niome_hotkey7/result/2026-09-18T20:48:23](../../data/inst/niome_hotkey7/result/2026-09-18T20:48:23/submission.json) |
| h8 | 2026-09-18 20:17:08,032 | [data/inst/niome_hotkey8/result/2026-09-18T20:17:58](../../data/inst/niome_hotkey8/result/2026-09-18T20:17:58/submission.json) |
| h9 | 2026-09-18 20:37:50,788 | [data/inst/niome_hotkey9/result/2026-09-18T20:38:19](../../data/inst/niome_hotkey9/result/2026-09-18T20:38:19/submission.json) |

## Verification and artifacts

- All ten upload archives identify this task; contracts and references match the published task after the seed stamp.
- All ten have 250 stage-1/2 valid rows, accessibility 0.35, and 80 Cas12a + 170 Cas9 rows.
- Independent local three-seed final scores match all ten official API scores exactly.
- Replayed outcomes and indel lengths match local stage-3 details for every uploaded row at all three scoring seeds.
- Recovered clean sets and HDR bands over **900 × 250 × 10 = 2,250,000 row/seed combinations**.
- [Full analysis](analysis.json), [official score snapshot](scores.json), [published task](task.json), [sanitized build/upload evidence](build_evidence.json), [archive manifest](manifest.json), [summary CSV](summary.csv).
- Reproduce from the repository root: `.venv/bin/python reports/14dcd43f/validate.py`, then `.venv/bin/python reports/14dcd43f/audit.py`, then `.venv/bin/python reports/14dcd43f/render_report.py`. Validations write only to report-local directories. The analysis includes validator-code and submission SHA-256 hashes.

Public hotkey mapping:

- **h0 / UID 248**: `5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW`.
- **h1 / UID 41**: `5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb`.
- **h2 / UID 163**: `5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU`.
- **h3 / UID 118**: `5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf`.
- **h4 / UID 190**: `5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2`.
- **h5 / UID 234**: `5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn`.
- **h6 / UID 175**: `5H5v45M2i6cFtS3Di4abh2zjWrPmuJZzNPoAYtVAXVFa6sp3`.
- **h7 / UID 65**: `5CS8FdHr8Ddv14QHr75E5K7e9xv3wre9m66ywRHZMR8zhntN`.
- **h8 / UID 211**: `5Cr9gJ3ukDDdxnhpyRM58gdMz7Sjuw8u58iiuYphVz1ZUGD2`.
- **h9 / UID 245**: `5Ge6UdGJE7BabD3bv8L9AR79rW47K3XQetcTHgaLD9NF2Vf5`.
