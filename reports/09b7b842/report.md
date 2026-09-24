# Submission audit: 09b7b842-ccff-4ea5-ab63-7486437e1b59

**HEK293**, created **2026-09-21T06:05:14.343405 UTC**. Scoring seeds: **423, 999, 189**.

**Best of h0–h9: h4, 18.979722, rank 99/248. Leader: UID 148, 100.732020.** h4 reached 18.84% of the leader, a gap of 81.752298. Top-10 cutoff: 86.640083.

All ten hotkeys uploaded prepared conjunctions with **250 valid rows: 80 Cas12a + 170 Cas9**. Each uploaded an **8-seed HDR band**. No fallback was used. **None of the three scoring seeds was clean or all-HDR on any hotkey.**

**This task spans a configuration change. h0–h3 and h9 uploaded the original narrow-cut builds; h4–h8 uploaded the later shared-window, wide-cut builds.** The latest `window_used.json` entries describe the rebuild for all ten and therefore contain the wrong uploaded bands for h0–h3 and h9. This audit derives bands and clean sets from the actual archived submission rows.

[Official score feed](https://niome-api.genomes.io/api/v3/miners/scores?task_id=09b7b842-ccff-4ea5-ab63-7486437e1b59&limit=40000), fetched 2026-09-21T08:36:24.476070+00:00. 248 records, 248 miners, 1 validator. Score timestamp: 2026-09-21T08:29:34.387183 UTC. Audit generated 2026-09-21T08:42:09.508604+00:00.

## Official scores

| Miner | UID | Rank / 248 | Final score | Weighted score | Consistency | Fidelity | Valid rows | Reward weight |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Top miner | 148 | 1 | 100.732020 | 208.103488 | 0.504374 | 0.959700 | 250 | 0.294000 |
| h0 | 248 | 220 | 15.512824 | 224.629010 | 0.075823 | 0.910806 | 250 | 0.000000 |
| h1 | 189 | 205 | 16.857116 | 223.814243 | 0.082482 | 0.913139 | 250 | 0.000000 |
| h2 | 45 | 224 | 14.657636 | 224.406155 | 0.072004 | 0.907135 | 250 | 0.000000 |
| h3 | 136 | 218 | 16.246405 | 223.892013 | 0.079106 | 0.917295 | 250 | 0.000000 |
| h4 | 190 | 99 | 18.979722 | 221.712582 | 0.092307 | 0.927390 | 250 | 0.000000 |
| h5 | 234 | 164 | 17.822489 | 224.696983 | 0.086199 | 0.920173 | 250 | 0.000000 |
| h6 | 175 | 207 | 16.729445 | 221.774975 | 0.081394 | 0.926786 | 250 | 0.000000 |
| h7 | 65 | 195 | 17.064984 | 225.286484 | 0.082609 | 0.916949 | 250 | 0.000000 |
| h8 | 211 | 190 | 17.279071 | 222.757032 | 0.084032 | 0.923094 | 250 | 0.000000 |
| h9 | 245 | 221 | 15.323244 | 221.362761 | 0.074778 | 0.925702 | 250 | 0.000000 |

Ranks count strictly higher scores, so ties share rank. All identities were matched by full public hotkeys. h1–h3 appear at UIDs **189, 45 and 136** for this task; their older UIDs are not used. Final scores and breakdown factors are means across three seeds. The product of mean factors need not equal the mean final score.

## Joined seed windows and cut search

The original round plan (logged **06:04 UTC**) predicted width-100 windows **100 / 200 / 300**, joined as **100–399**. The rebuild's plan (generated **06:17:51 UTC**) records the same prediction. Actual scoring windows were **400 / 900 / 100** for seeds **423 / 999 / 189**: **one matched window**.

The later shared layout instead fixes its band candidates to **100–299 ∪ 500–599** (classes **100 / 200 / 500**, 300 seeds), independently of the plan. It likewise covers one of the three actual windows. The latest record's `window: [100,399]` describes the plan, not the rebuild's actual shared band candidates. h4–h8 submitted shared loops **5–9 of 10**.

Original h1–h3 each use a 100-seed band slice, while their cut search spans the full predicted 300 seeds. Original h0/h9 each use 300 band candidates. Rebuilt h4–h8 use all 900 seeds for cut selection. Band candidates and cut-search seeds are separate quantities.

| Hotkey | Uploaded build | Joined band candidates | Candidate count | Shared loop | Cut search | Clean in cut search | Clean over 100–999 |
|---|---|---|---:|---:|---|---:|---:|
| h0 | original / narrow | 100–399 | 300 | — | 100–399 | 23/300 | 23 |
| h1 | original / narrow | 100–199 | 100 | — | 100–399 | 23/300 | 23 |
| h2 | original / narrow | 200–299 | 100 | — | 100–399 | 24/300 | 24 |
| h3 | original / narrow | 300–399 | 100 | — | 100–399 | 25/300 | 25 |
| h4 | rebuilt / wide | 100–299 ∪ 500–599 | 300 | 5 | 100–999 | 18/900 | 18 |
| h5 | rebuilt / wide | 100–299 ∪ 500–599 | 300 | 6 | 100–999 | 18/900 | 18 |
| h6 | rebuilt / wide | 100–299 ∪ 500–599 | 300 | 7 | 100–999 | 19/900 | 19 |
| h7 | rebuilt / wide | 100–299 ∪ 500–599 | 300 | 8 | 100–999 | 20/900 | 20 |
| h8 | rebuilt / wide | 100–299 ∪ 500–599 | 300 | 9 | 100–999 | 17/900 | 17 |
| h9 | original / narrow | 700–999 | 300 | — | 700–999 | 22/300 | 22 |

## Uploaded HDR bands

**Clean** means all 250 rows cut. **HDR band** means all 250 rows return HDR. Sets were recovered by replay over every seed from 100 to 999. Every band is contained in its uploaded candidate window and clean set. The build's Cas12a clean count is checked within its cut-search space, separately from the full-row clean definition.

| Hotkey | Clean /900 | Exact uploaded HDR band | Clean scoring-seed hits | HDR scoring-seed hits | Latest band record matches upload |
|---|---:|---|---|---|---|
| h0 | 23 | 111, 140, 152, 169, 205, 260, 262, 336 | none | none | NO — later rebuild |
| h1 | 23 | 124, 127, 129, 137, 140, 141, 175, 185 | none | none | NO — later rebuild |
| h2 | 24 | 203, 205, 215, 228, 235, 262, 272, 298 | none | none | NO — later rebuild |
| h3 | 25 | 303, 316, 325, 329, 339, 345, 360, 394 | none | none | NO — later rebuild |
| h4 | 18 | 101, 123, 127, 236, 244, 262, 284, 591 | none | none | yes |
| h5 | 18 | 179, 203, 209, 217, 289, 513, 538, 557 | none | none | yes |
| h6 | 19 | 105, 160, 208, 249, 254, 271, 282, 593 | none | none | yes |
| h7 | 20 | 100, 190, 222, 280, 529, 549, 559, 587 | none | none | yes |
| h8 | 17 | 143, 167, 215, 233, 238, 285, 548, 558 | none | none | yes |
| h9 | 22 | 719, 723, 737, 796, 812, 931, 946, 952 | none | none | NO — later rebuild |

Fleet union: **173 distinct clean seeds** and **73 distinct HDR-band seeds**, from 80 band memberships. The five uploaded shared-loop bands are mutually disjoint, but bands across the mixed fleet overlap. Overlapping seeds: 127 (h1/h4); 140 (h0/h1); 203 (h2/h5); 205 (h0/h2); 215 (h2/h8); 262 (h0/h2/h4).

## Per-seed results

Each entry is **per-seed final score / HDR rows / no-cut rows**, from 250 rows.

| Hotkey | Seed 423 | Seed 999 | Seed 189 | Mean final score |
|---|---:|---:|---:|---:|
| h0 | 13.288583 / 104 / 23 | 18.029162 / 120 / 18 | 15.220727 / 105 / 16 | 15.512824 |
| h1 | 17.639009 / 101 / 16 | 17.345952 / 102 / 20 | 15.586387 / 107 / 18 | 16.857116 |
| h2 | 16.421931 / 118 / 18 | 14.719974 / 112 / 19 | 12.831002 / 102 / 19 | 14.657636 |
| h3 | 16.094782 / 96 / 19 | 16.481080 / 89 / 13 | 16.163353 / 101 / 18 | 16.246405 |
| h4 | 20.662804 / 109 / 11 | 18.360999 / 105 / 16 | 17.915363 / 103 / 11 | 18.979722 |
| h5 | 21.637369 / 114 / 12 | 16.571093 / 117 / 16 | 15.259004 / 103 / 16 | 17.822489 |
| h6 | 15.844859 / 98 / 20 | 18.090712 / 98 / 14 | 16.252766 / 112 / 14 | 16.729445 |
| h7 | 18.724063 / 106 / 20 | 17.016649 / 103 / 21 | 15.454241 / 101 / 22 | 17.064984 |
| h8 | 20.403976 / 116 / 10 | 15.947268 / 104 / 14 | 15.485970 / 109 / 19 | 17.279071 |
| h9 | 11.758505 / 102 / 17 | 20.015210 / 107 / 9 | 14.196016 / 93 / 20 | 15.323244 |

## Comparison with the leader

h4's weighted score was **221.712582**, versus the leader's **208.103488**. Mean consistency was **0.092307 versus 0.504374**, and fidelity was **0.927390 versus 0.959700**. The largest factor gap is consistency. Our best weighted score exceeds the leader's, but its much lower consistency accompanies a much lower final score. Every uploaded hotkey missed all three seeds with its full clean set and HDR band. This single mixed-configuration task does not establish which layout performs better across tasks.

Leader hotkey: `5GxfgRsx17PhVECHZjhjYzPA2qSKXXXsVHW49gJA41SDYSho`. Its submitted rows, clean sets, HDR bands and per-seed outcomes were unavailable in the inspected data. Its aggregate score does not identify those sets.

## Exact uploaded clean sets

These lists cover the full 100–999 range and include the corresponding HDR bands.

### h0: 23 clean seeds

111, 115, 117, 119, 136, 140, 152, 156, 166, 169, 171, 188, 205, 223, 227, 236, 260, 262, 266, 297, 336, 358, 394.

### h1: 23 clean seeds

106, 112, 119, 124, 125, 127, 129, 130, 137, 140, 141, 169, 175, 185, 220, 239, 265, 285, 314, 340, 366, 367, 389.

### h2: 24 clean seeds

131, 152, 188, 192, 203, 205, 211, 212, 215, 218, 219, 221, 228, 235, 241, 258, 259, 262, 272, 284, 298, 328, 375, 378.

### h3: 25 clean seeds

148, 160, 176, 180, 185, 192, 215, 223, 227, 239, 259, 265, 300, 303, 311, 316, 325, 327, 329, 339, 342, 345, 360, 381, 394.

### h4: 18 clean seeds

101, 123, 125, 127, 170, 180, 236, 244, 262, 284, 313, 368, 530, 582, 591, 693, 719, 820.

### h5: 18 clean seeds

179, 203, 209, 217, 289, 308, 398, 409, 450, 513, 538, 557, 726, 749, 890, 901, 956, 991.

### h6: 19 clean seeds

102, 105, 110, 160, 208, 249, 254, 259, 271, 282, 356, 358, 486, 566, 593, 834, 874, 878, 942.

### h7: 20 clean seeds

100, 127, 190, 222, 223, 236, 257, 274, 280, 398, 448, 529, 549, 559, 587, 687, 719, 764, 796, 913.

### h8: 17 clean seeds

143, 147, 167, 215, 233, 238, 285, 421, 432, 548, 558, 689, 713, 772, 837, 924, 991.

### h9: 22 clean seeds

714, 719, 723, 737, 739, 768, 796, 812, 892, 909, 926, 931, 946, 949, 952, 953, 957, 958, 978, 983, 992, 993.

## Build and upload provenance

Original builds completed around 06:12–06:13 UTC. Rebuilds began at 06:39 and completed around 06:42. Each log explicitly reports using the prepared 250-row submission. The timestamps below select the build actually available at upload time, and replay confirms the uploaded band and clean count.

| Hotkey | Uploaded build completed UTC | Uploaded UTC | Latest window record UTC |
|---|---|---|---|
| h0 | 2026-09-21 06:13:21,145 | 2026-09-21 06:23:37,329 | 2026-09-21T06:39:06.002694+00:00 |
| h1 | 2026-09-21 06:12:47,163 | 2026-09-21 06:19:17,521 | 2026-09-21T06:39:07.958830+00:00 |
| h2 | 2026-09-21 06:12:47,418 | 2026-09-21 06:30:54,436 | 2026-09-21T06:39:10.182799+00:00 |
| h3 | 2026-09-21 06:12:47,416 | 2026-09-21 06:33:01,297 | 2026-09-21T06:39:11.353513+00:00 |
| h4 | 2026-09-21 06:42:34,939 | 2026-09-21 06:43:18,007 | 2026-09-21T06:39:12.747146+00:00 |
| h5 | 2026-09-21 06:42:38,827 | 2026-09-21 06:59:48,410 | 2026-09-21T06:39:17.275310+00:00 |
| h6 | 2026-09-21 06:42:40,599 | 2026-09-21 06:53:25,319 | 2026-09-21T06:39:17.649850+00:00 |
| h7 | 2026-09-21 06:42:43,177 | 2026-09-21 06:44:37,121 | 2026-09-21T06:39:18.584192+00:00 |
| h8 | 2026-09-21 06:42:46,733 | 2026-09-21 06:57:13,904 | 2026-09-21T06:39:20.665598+00:00 |
| h9 | 2026-09-21 06:13:27,267 | 2026-09-21 06:21:08,063 | 2026-09-21T06:39:24.812127+00:00 |

## Verification and artifacts

- All 10 archive upload records identify this task and a successful upload. Contracts and references match the published task after substituting its seed stamp.
- All 2,500 rows passed validation with HEK293 accessibility 0.35; all 10 local three-seed final scores match official scores within 1e-10.
- Reused h7's archived validation; generated the other nine validations under this report directory.
- Replayed every row at all 900 seeds: **2,250,000 row/seed combinations**. At the three scoring seeds, every outcome and indel length matches local stage-3 detail.
- Clean counts within cut-search windows match all selected build logs. Recovered bands match the latest records for h4–h8 and differ for h0–h3/h9, as their upload timing predicts.
- [Analysis JSON](analysis.json), [CSV](summary.csv), [official scores](scores.json), [task](task.json), [manifest](manifest.json), [sanitized build evidence](build_evidence.json), [layout evidence](layout_evidence.json).
- Reproduce with `.venv/bin/python reports/09b7b842/validate.py`, then `.venv/bin/python reports/09b7b842/audit.py`, then `.venv/bin/python reports/09b7b842/render_report.py`. Validator and submission SHA-256 hashes are included in the analysis.

| Hotkey | Public hotkey | Submission archive |
|---|---|---|
| h0 | `5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW` | [data/inst/niome_hotkey/result/2026-09-21T06:24:07](../../data/inst/niome_hotkey/result/2026-09-21T06:24:07/submission.json) |
| h1 | `5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb` | [data/inst/niome_hotkey1/result/2026-09-21T06:19:47](../../data/inst/niome_hotkey1/result/2026-09-21T06:19:47/submission.json) |
| h2 | `5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU` | [data/inst/niome_hotkey2/result/2026-09-21T06:31:24](../../data/inst/niome_hotkey2/result/2026-09-21T06:31:24/submission.json) |
| h3 | `5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf` | [data/inst/niome_hotkey3/result/2026-09-21T06:33:30](../../data/inst/niome_hotkey3/result/2026-09-21T06:33:30/submission.json) |
| h4 | `5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2` | [data/inst/niome_hotkey4/result/2026-09-21T06:43:46](../../data/inst/niome_hotkey4/result/2026-09-21T06:43:46/submission.json) |
| h5 | `5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn` | [data/inst/niome_hotkey5/result/2026-09-21T07:00:17](../../data/inst/niome_hotkey5/result/2026-09-21T07:00:17/submission.json) |
| h6 | `5H5v45M2i6cFtS3Di4abh2zjWrPmuJZzNPoAYtVAXVFa6sp3` | [data/inst/niome_hotkey6/result/2026-09-21T06:53:54](../../data/inst/niome_hotkey6/result/2026-09-21T06:53:54/submission.json) |
| h7 | `5CS8FdHr8Ddv14QHr75E5K7e9xv3wre9m66ywRHZMR8zhntN` | [data/inst/niome_hotkey7/result/2026-09-21T06:45:06](../../data/inst/niome_hotkey7/result/2026-09-21T06:45:06/submission.json) |
| h8 | `5Cr9gJ3ukDDdxnhpyRM58gdMz7Sjuw8u58iiuYphVz1ZUGD2` | [data/inst/niome_hotkey8/result/2026-09-21T06:57:42](../../data/inst/niome_hotkey8/result/2026-09-21T06:57:42/submission.json) |
| h9 | `5Ge6UdGJE7BabD3bv8L9AR79rW47K3XQetcTHgaLD9NF2Vf5` | [data/inst/niome_hotkey9/result/2026-09-21T06:21:36](../../data/inst/niome_hotkey9/result/2026-09-21T06:21:36/submission.json) |
