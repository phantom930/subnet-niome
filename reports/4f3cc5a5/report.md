# Submission audit: 4f3cc5a5-cbb5-4a1e-86c4-2c88d2e51697

HEK293 · task opened 2026-09-18T00:35:26.785386 UTC · scoring seeds **758, 865, 350**.

**Best of h0–h5: h3, 26.608822, rank 96/248. Top miner: UID 90, 184.277221.** h3 reached 14.44% of the top score, a gap of 157.668399 points. It missed the top-10 cutoff (129.592319) by **102.983498** points. All six submissions were uploaded successfully, contained 250 valid rows, and received zero reward weight.

Scores are from the [task-specific public API](https://niome-api.genomes.io/api/v3/miners/scores?task_id=4f3cc5a5-cbb5-4a1e-86c4-2c88d2e51697&limit=40000); 248 records, 248 distinct miners, 1 validator. Score timestamp: 2026-09-18T03:02:18.241953 UTC. Audit generated: 2026-09-18T04:18:05.807656+00:00.

## Official score comparison

| Miner | UID this round | Rank / 248 | Final score | Weighted score | Consistency | Fidelity | Reward weight |
|---|---:|---:|---:|---:|---:|---:|---:|
| Top miner | 90 | 1 | 184.277221 | 351.164010 | 0.554774 | 0.945900 | 0.294000 |
| h0 | 248 | 177 | 23.531965 | 338.134383 | 0.077172 | 0.901793 | 0.000000 |
| h1 | 41 | 190 | 21.574193 | 333.382188 | 0.071041 | 0.910926 | 0.000000 |
| h2 | 163 | 154 | 24.913252 | 344.606031 | 0.081205 | 0.890280 | 0.000000 |
| h3 | 118 | 96 | 26.608822 | 342.152833 | 0.087072 | 0.893154 | 0.000000 |
| h4 | 190 | 189 | 21.652004 | 336.437605 | 0.070762 | 0.909476 | 0.000000 |
| h5 | 234 | 170 | 24.147755 | 334.505421 | 0.079128 | 0.912316 | 0.000000 |

Top miner hotkey: `5DkDdLMvYdts4ADidWhTvAWqjZeMXnS3kqcPByobU4XbvFz6`.

Final score is the mean of the per-seed scores. For these submissions weighted score and fidelity are constant across seeds, so final = weighted score × mean consistency × fidelity.

## Joined seed windows and recovered bands

All six original builds used joined band space **200–399 ∪ 900–999** (300 seeds), from the plan generated at 2026-09-18 00:17 UTC. **The cut search also used this same 300-seed joined window.** HEK293 was on the narrower cut configuration for this task; the logs report clean counts out of 300, not 900. Band candidate slices were 150 seeds wide, at offsets 0/50/100/150/200/250 for h0–h5. Every submission has 80 Cas12a and 170 Cas9 rows and an eight-seed HDR band.

**Clean set:** every submitted row cuts at that seed. **HDR band:** every submitted row returns HDR. The band is a subset of the clean set. Both sets were recovered by simulating the exact uploaded rows at every seed 100–999; the recovered bands match the original recorded bands exactly. Inside the joined window the full-row clean sets equal the Cas12a group clean sets, and their sizes match the build logs. The full 900-seed replay found **no additional clean seeds outside the joined window**, either for the full submission or for the Cas12a group alone.

| Hotkey | Candidate slice within joined space | Clean / 300 | Full scan clean / 900 | Actual HDR band seeds |
|---|---|---:|---:|---|
| h0 | 200–349 | 23 | 23 | 213, 311, 313, 314, 317, 328, 333, 347 |
| h1 | 250–399 | 25 | 25 | 255, 261, 263, 279, 313, 327, 351, 398 |
| h2 | 300–399, 900–949 | 20 | 20 | 309, 313, 314, 370, 374, 386, 398, 900 |
| h3 | 350–399, 900–999 | 24 | 24 | 362, 363, 392, 911, 958, 963, 980, 983 |
| h4 | 200–249, 900–999 | 18 | 18 | 212, 219, 237, 934, 939, 942, 958, 982 |
| h5 | 200–299, 950–999 | 21 | 21 | 213, 219, 237, 257, 265, 958, 959, 993 |

Candidate slices are reconstructed from the original logged joined space and the six-hotkey offset layout; every recovered band lies inside its corresponding slice.

Fleet unions: **98 distinct clean seeds** and **39 distinct HDR band seeds** (48 band memberships before overlap). The scoring seeds covered by the fleet's clean union are none; the band union hits none.

## Hits and per-seed scores

| Hotkey | Clean seed hits (includes HDR hits) | HDR band hits |
|---|---|---|
| h0 | none | none |
| h1 | none | none |
| h2 | none | none |
| h3 | none | none |
| h4 | none | none |
| h5 | none | none |

Seeds **758 and 865** are outside the joined window. Seed **350** is inside the joined window and the candidate slices of h1/h2/h3, but it is absent from every clean set and HDR band. **None of the three scoring seeds hits a clean set or HDR band for any of h0–h5.**

Each cell below is **per-seed final score / no-cut rows**.

| Hotkey | Seed 758 | Seed 865 | Seed 350 | Mean final score |
|---|---:|---:|---:|---:|
| h0 | 22.419914 / 17 | 26.024110 / 11 | 22.151872 / 19 | 23.531965 |
| h1 | 17.900947 / 28 | 25.999028 / 18 | 20.822603 / 13 | 21.574193 |
| h2 | 25.591676 / 14 | 22.029980 / 22 | 27.118101 / 14 | 24.913252 |
| h3 | 24.600902 / 19 | 26.895464 / 16 | 28.330099 / 17 | 26.608822 |
| h4 | 20.046574 / 19 | 26.618661 / 22 | 18.290776 / 19 | 21.652004 |
| h5 | 26.466664 / 18 | 25.131812 / 20 | 20.844788 / 23 | 24.147755 |

All 18 hotkey/seed combinations contain no-cut rows. h3's per-seed scores are **24.600902** at 758, **26.895464** at 865, and **28.330099** at 350, giving the fleet's best three-seed mean of **26.608822**.

## Comparison with the leader

h3's weighted score is only modestly lower than the leader's (**342.152833 versus 351.164010**). Its fidelity is also lower (**0.893154 versus 0.945900**), but the dominant gap is mean consistency: **0.087072 versus 0.554774**. Across all six hotkeys, consistency ranges from 0.070762 to 0.087072. Missing both the clean sets and the HDR bands explains the fleet's per-seed outcomes on this task; one task does not establish long-run strategy performance. The leader's clean sets, bands, and per-seed outcomes are not available from its aggregate score and were not reconstructed.

## Exact clean sets

These lists include the HDR band seeds and cover the full 100–999 seed space. Machine-readable lists are also in [analysis.json](analysis.json).

### h0: 23 clean seeds

210, 213, 237, 253, 260, 284, 311, 313, 314, 317, 320, 321, 328, 332, 333, 347, 382, 397, 913, 914, 926, 945, 953.

### h1: 25 clean seeds

202, 244, 251, 255, 256, 259, 261, 263, 279, 313, 314, 325, 327, 337, 338, 351, 357, 370, 398, 900, 957, 969, 990, 991, 998.

### h2: 20 clean seeds

223, 277, 309, 310, 313, 314, 370, 374, 386, 396, 398, 900, 903, 927, 940, 942, 945, 954, 980, 994.

### h3: 24 clean seeds

219, 258, 259, 297, 328, 333, 351, 362, 363, 392, 397, 905, 911, 913, 921, 943, 957, 958, 963, 972, 974, 980, 983, 993.

### h4: 18 clean seeds

212, 219, 237, 261, 309, 316, 323, 356, 369, 381, 904, 905, 934, 939, 942, 958, 959, 982.

### h5: 21 clean seeds

213, 219, 231, 237, 257, 259, 265, 284, 341, 381, 382, 916, 924, 948, 958, 959, 962, 968, 976, 987, 993.

## Submission provenance and verification

The six builds started between 00:35:34 and 00:35:35 UTC and finished between 00:38:12 and 00:38:18 UTC (158–163 seconds). All six hotkeys served the validator's requests from prepared submissions. **h3 was called before its build finished and waited about 70 seconds**, then uploaded the completed conjunction submission with 213 seconds of TTL remaining. It did not use an emergency fallback. Uploads succeeded between 00:38:13 and 01:34:28 UTC, with 213–289 seconds of TTL remaining. Each hotkey has one recorded build and one recorded successful upload for this task. The recorded bands agree with the archived rows' full-seed replay.

- All archived contracts and reference files match the published task after accounting for the seed changing from 0 at submission time to the three stamped scoring seeds.
- All 250 archived stage-1/2 valid-row entries match each submission exactly; accessibility is 0.35.
- All six archived three-seed final scores match the official API exactly.
- Replayed stage-3 outcomes and indel lengths match the archived detail for every row at all three real seeds.
- Scanned 900 seeds × 250 rows × 6 hotkeys = 1,350,000 row simulations. No miner configuration, running process, or uploaded submission was changed.
- [Full analysis](analysis.json), [official score snapshot](scores.json), [published task](task.json), [sanitized build/upload evidence](build_evidence.json).
- Reproduce with `.venv/bin/python reports/4f3cc5a5/audit.py` and `.venv/bin/python reports/4f3cc5a5/render_report.py` from the repository root.

Uploaded archives:

- **h0** `5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW` — [data/inst/niome_hotkey/result/2026-09-18T01:34:00](../../data/inst/niome_hotkey/result/2026-09-18T01:34:00/submission.json)
- **h1** `5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb` — [data/inst/niome_hotkey1/result/2026-09-18T01:34:58](../../data/inst/niome_hotkey1/result/2026-09-18T01:34:58/submission.json)
- **h2** `5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU` — [data/inst/niome_hotkey2/result/2026-09-18T01:09:41](../../data/inst/niome_hotkey2/result/2026-09-18T01:09:41/submission.json)
- **h3** `5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf` — [data/inst/niome_hotkey3/result/2026-09-18T00:38:42](../../data/inst/niome_hotkey3/result/2026-09-18T00:38:42/submission.json)
- **h4** `5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2` — [data/inst/niome_hotkey4/result/2026-09-18T00:58:07](../../data/inst/niome_hotkey4/result/2026-09-18T00:58:07/submission.json)
- **h5** `5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn` — [data/inst/niome_hotkey5/result/2026-09-18T01:12:21](../../data/inst/niome_hotkey5/result/2026-09-18T01:12:21/submission.json)
