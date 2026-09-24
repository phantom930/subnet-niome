# Submission audit: 5a453e21-5516-415c-a2f9-2ad2d04a75e2

**HUDEP-2**, created **2026-09-23T08:29:24.352781 UTC**.

**Best of h0–h6: h3, 118.755026, rank 10/248. Leader: UID 4, 181.923257.** h3 reached 65.28% of the leader, a gap of 63.168231. Top-10 cutoff: 118.755026.

All seven hotkeys uploaded prepared conjunctions with **250 valid rows: 100 Cas12a + 150 Cas9**. Each upload has an **11-seed HDR band**. No fallback occurred.

**Scoring seeds are unavailable in the inspected sources.** The public task feed still carries `seed: 0`, although the official score feed has results. Fresh requests with pagination and different query parameters also returned 0. The public S3 contract returned HTTP 403. The 0 is an unstamped placeholder, not a reported scoring seed. **Clean/HDR hits on the actual scoring seeds and an exact reproduction of final scores cannot be confirmed from this snapshot.** They are reported as unavailable, not zero.

[Official score feed](https://niome-api.genomes.io/api/v3/miners/scores?task_id=5a453e21-5516-415c-a2f9-2ad2d04a75e2&limit=40000), fetched 2026-09-23T10:50:52.580992+00:00. 248 records, 248 miners, 1 validator. Score timestamp: 2026-09-23T10:43:22.581380 UTC. Audit generated 2026-09-23T10:54:51.842563+00:00.

## Official scores

| Miner | UID | Rank / 248 | Final score | Weighted score | Consistency | Fidelity | Valid rows | Reward weight |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Top miner | 4 | 1 | 181.923257 | 320.974992 | 0.579533 | 0.978000 | 250 | 0.29400003 |
| h0 | 226 | 168 | 39.614837 | 314.719984 | 0.133469 | 0.943090 | 250 | 0.00000000 |
| h1 | 189 | 240 | 29.108803 | 317.680017 | 0.097484 | 0.939944 | 250 | 0.00000000 |
| h2 | 45 | 199 | 35.875951 | 320.726165 | 0.119628 | 0.935051 | 250 | 0.00000000 |
| h3 | 136 | 10 | 118.755026 | 320.178962 | 0.395557 | 0.937670 | 250 | 0.00980000 |
| h4 | 228 | 205 | 35.644777 | 314.586365 | 0.120059 | 0.943755 | 250 | 0.00000000 |
| h5 | 49 | 166 | 45.455242 | 311.311069 | 0.154222 | 0.946765 | 250 | 0.00000000 |
| h6 | 227 | 232 | 32.739887 | 306.680657 | 0.112331 | 0.950363 | 250 | 0.00000000 |

Ranks count strictly higher scores; ties share rank. All identities were matched by full public hotkeys, using the UIDs in this task's official score snapshot. h3 is the only requested hotkey with nonzero reward weight: **0.0098 (0.98%)**. The leader's weight is approximately **0.294 (29.4%)**. These are reported task weights, not a calculation of on-chain payout.

## Joined seed windows and actual layout

The task used a **fixed layout**, independent of the prediction plan. The fleet's band space is **200–299**, divided into seven circular windows of **15 seeds at stride 15**. Every hotkey uses **loop 1**. All seven share the **200-seed cut window 200–299 ∪ 400–499**. The log text `shared window, loop 1 of 7` does not mean every hotkey uses the same band candidates: each has its own 15-seed slice.

**h6 wraps:** its candidates are **200–204 ∪ 290–299**, not all 200–299. Its latest record's `[200,299]` field is only the bounding range. The build log and layout source preserve the exact joined window. h0 and h6 share five candidate seeds (200–204).

| Hotkey | Joined band candidates | Offset | Loop | Clean in shared cut /200 | Full clean /900 | Exact uploaded HDR band |
|---|---|---:|---:|---:|---:|---|
| h0 | 200–214 | 0 | 1 | 86 | 88 | 200, 201, 202, 205, 206, 207, 208, 209, 212, 213, 214 |
| h1 | 215–229 | 15 | 1 | 85 | 91 | 215, 216, 217, 218, 220, 221, 222, 223, 224, 226, 229 |
| h2 | 230–244 | 30 | 1 | 85 | 87 | 230, 231, 233, 235, 236, 237, 238, 239, 240, 241, 244 |
| h3 | 245–259 | 45 | 1 | 84 | 86 | 246, 247, 248, 249, 250, 251, 252, 253, 254, 257, 259 |
| h4 | 260–274 | 60 | 1 | 81 | 84 | 260, 261, 262, 263, 265, 267, 268, 269, 271, 272, 274 |
| h5 | 275–289 | 75 | 1 | 88 | 92 | 275, 277, 279, 280, 281, 282, 284, 286, 287, 288, 289 |
| h6 | 200–204 ∪ 290–299 | 90 | 1 | 86 | 86 | 200, 201, 203, 290, 291, 293, 294, 295, 296, 297, 298 |

**Clean** means all 250 uploaded rows cut at that seed. **HDR band** means all 250 rows return HDR. Both sets were recovered by replay over every seed from 100 through 999. The build logs count clean seeds inside the 200-seed cut window; full replay can find additional clean seeds outside it. Cas12a-only cleanliness is checked separately and is not counted as full-row cleanliness when a Cas9 row fails to cut.

Fleet union: **214 distinct clean seeds** and **75 distinct HDR-band seeds** from 77 band memberships. Clean union inside the shared cut window: **195/200**. Repeated HDR-band seeds: 200 (h0/h6); 201 (h0/h6).

## Comparison with the top miner

h3's weighted score was **320.178962**, versus the leader's **320.974992**. Mean consistency was **0.395557 versus 0.579533**, and fidelity was **0.937670 versus 0.978000**. Weighted scores are close; the largest factor gap is consistency. An aggregate final score alone does not establish which seeds hit a band. For multi-seed scoring, multiplying mean factors need not reproduce the mean final score.

Leader hotkey: `5GTfCDV5fxTpEWMChqYdXT4cCQv2VwJ5wWJ8idA2g2dGTnbA`. The leader's submission rows and exact clean/HDR sets were not available in the inspected data.

## Exact uploaded clean sets

Each list covers the entire 100–999 range and includes its uploaded HDR band.

### h0: 88 clean seeds

200, 201, 202, 203, 205, 206, 207, 208, 209, 210, 212, 213, 214, 215, 224, 225, 226, 227, 228, 229, 235, 236, 237, 238, 240, 244, 246, 247, 248, 250, 251, 252, 254, 256, 257, 261, 263, 264, 266, 267, 269, 271, 274, 281, 282, 283, 284, 289, 291, 292, 294, 297, 298, 400, 401, 406, 408, 415, 419, 425, 430, 432, 436, 437, 441, 443, 447, 448, 451, 454, 460, 464, 468, 469, 471, 472, 477, 481, 485, 487, 488, 489, 496, 497, 498, 499, 827, 861.

Clean seeds outside the shared cut window: 827, 861.

Cas12a-only clean seeds excluded from the full-row set: 154, 161, 322, 352, 559, 599, 640, 667, 701, 806, 819, 842, 851, 994.

### h1: 91 clean seeds

204, 205, 206, 210, 211, 213, 214, 215, 216, 217, 218, 220, 221, 222, 223, 224, 226, 229, 230, 231, 232, 236, 237, 248, 250, 251, 255, 259, 263, 265, 266, 267, 269, 271, 273, 274, 275, 278, 279, 283, 285, 287, 292, 294, 296, 326, 402, 403, 405, 409, 410, 411, 413, 418, 419, 422, 423, 424, 425, 426, 429, 430, 438, 444, 446, 447, 448, 453, 459, 461, 464, 466, 471, 472, 477, 478, 479, 483, 485, 488, 490, 491, 494, 495, 496, 498, 517, 661, 681, 806, 915.

Clean seeds outside the shared cut window: 326, 517, 661, 681, 806, 915.

Cas12a-only clean seeds excluded from the full-row set: 107, 364, 390, 501, 518, 522, 585, 776, 788, 820, 975.

### h2: 87 clean seeds

201, 203, 205, 206, 207, 208, 213, 216, 217, 218, 220, 221, 224, 227, 229, 230, 231, 232, 233, 235, 236, 237, 238, 239, 240, 241, 243, 244, 246, 248, 253, 254, 255, 263, 265, 266, 273, 274, 275, 277, 280, 282, 285, 286, 291, 292, 296, 297, 298, 299, 406, 407, 409, 410, 411, 415, 416, 419, 420, 421, 422, 428, 430, 431, 432, 436, 437, 442, 448, 449, 451, 454, 456, 461, 462, 466, 468, 469, 470, 478, 480, 484, 485, 495, 499, 853, 866.

Clean seeds outside the shared cut window: 853, 866.

Cas12a-only clean seeds excluded from the full-row set: 541, 571, 604, 753.

### h3: 86 clean seeds

200, 206, 207, 210, 212, 213, 214, 215, 219, 224, 229, 231, 234, 235, 237, 239, 240, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 257, 259, 261, 263, 264, 265, 266, 268, 269, 274, 275, 278, 279, 281, 283, 284, 286, 289, 293, 297, 298, 299, 403, 404, 408, 411, 416, 419, 421, 422, 423, 424, 428, 431, 434, 441, 443, 444, 452, 458, 459, 462, 465, 467, 471, 473, 474, 481, 485, 488, 489, 490, 492, 494, 495, 496, 499, 648, 739.

Clean seeds outside the shared cut window: 648, 739.

Cas12a-only clean seeds excluded from the full-row set: 377, 600, 835, 926, 970.

### h4: 84 clean seeds

172, 200, 203, 205, 211, 219, 220, 221, 222, 224, 225, 229, 231, 239, 244, 248, 250, 252, 253, 256, 257, 258, 260, 261, 262, 263, 265, 266, 267, 268, 269, 271, 272, 273, 274, 280, 281, 289, 290, 294, 297, 404, 405, 409, 410, 416, 419, 420, 421, 424, 425, 429, 431, 433, 435, 439, 443, 449, 450, 451, 454, 458, 460, 461, 462, 465, 466, 467, 468, 469, 478, 479, 481, 482, 483, 484, 485, 487, 490, 495, 496, 499, 843, 918.

Clean seeds outside the shared cut window: 172, 843, 918.

Cas12a-only clean seeds excluded from the full-row set: 303, 306, 355, 514, 539, 613, 786, 870, 972, 986.

### h5: 92 clean seeds

207, 208, 210, 211, 214, 217, 219, 221, 223, 224, 229, 237, 241, 243, 247, 248, 249, 252, 253, 257, 258, 259, 261, 262, 263, 265, 266, 272, 273, 274, 275, 277, 278, 279, 280, 281, 282, 284, 286, 287, 288, 289, 290, 294, 295, 333, 403, 404, 405, 406, 414, 417, 418, 420, 423, 424, 426, 433, 438, 439, 444, 445, 447, 448, 451, 454, 455, 457, 458, 462, 463, 465, 466, 468, 470, 471, 472, 474, 475, 476, 479, 480, 483, 485, 490, 491, 492, 493, 496, 500, 565, 836.

Clean seeds outside the shared cut window: 333, 500, 565, 836.

Cas12a-only clean seeds excluded from the full-row set: 302, 732, 735, 799, 849, 904, 969.

### h6: 86 clean seeds

200, 201, 203, 205, 206, 207, 209, 212, 218, 221, 225, 227, 228, 233, 234, 235, 240, 245, 246, 248, 249, 252, 261, 263, 266, 267, 269, 272, 273, 274, 275, 276, 277, 279, 280, 281, 284, 285, 286, 287, 289, 290, 291, 293, 294, 295, 296, 297, 298, 401, 404, 406, 407, 409, 411, 412, 413, 415, 420, 421, 428, 429, 432, 434, 438, 446, 447, 448, 450, 451, 453, 457, 459, 461, 464, 470, 473, 474, 476, 477, 481, 483, 485, 488, 494, 498.

Clean seeds outside the shared cut window: none.

Cas12a-only clean seeds excluded from the full-row set: 164, 180, 343, 363, 520, 651, 656, 701, 783, 822, 935, 951.

## Verification and provenance

- All seven successful upload archives identify this task. Contracts and references match the task snapshot field for field.
- All 1,750 rows pass stage 1/2 with HUDEP-2 accessibility 0.82, from the current public cell-types table.
- Local weighted scores and distribution-fidelity factors match all seven official records within 1e-9. These quantities do not depend on the scoring seed.
- Because the task seed stamp is unavailable, `validate.py` uses **diagnostic seed 100** to prepare validation features. Its local final scores are probes and are not presented as the official task scores.
- Replayed **900 × 250 × 7 = 1,575,000 row/seed combinations**. Diagnostic-seed outcomes and indel lengths match the local pipeline. All seven recovered bands match the recorded bands and lie inside their exact candidate windows.
- Full-row and Cas12a clean counts inside the shared cut window match all seven build log counts. All uploads occurred after their prepared builds completed.
- [Analysis JSON](analysis.json), [CSV](summary.csv), [official scores](scores.json), [task snapshot](task.json), [seed availability checks](seed_availability.json), [manifest](manifest.json), [sanitized build evidence](build_evidence.json), [layout evidence](layout_evidence.json).
- Reproduce coverage with `.venv/bin/python reports/5a453e21/validate.py`, then `.venv/bin/python reports/5a453e21/audit.py`, then `.venv/bin/python reports/5a453e21/render_report.py`. Validator and submission SHA-256 hashes are recorded in the analysis.

| Hotkey | Upload UTC | Prepared build completed UTC | Public hotkey | Submission archive |
|---|---|---|---|---|
| h0 | 2026-09-23 09:13:51,914 | 2026-09-23 08:35:29,601 | `5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW` | [data/inst/niome_hotkey/result/2026-09-23T09:14:06](../../data/inst/niome_hotkey/result/2026-09-23T09:14:06/submission.json) |
| h1 | 2026-09-23 09:18:47,873 | 2026-09-23 08:35:30,646 | `5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb` | [data/inst/niome_hotkey1/result/2026-09-23T09:19:02](../../data/inst/niome_hotkey1/result/2026-09-23T09:19:02/submission.json) |
| h2 | 2026-09-23 09:34:52,638 | 2026-09-23 08:35:30,673 | `5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU` | [data/inst/niome_hotkey2/result/2026-09-23T09:35:06](../../data/inst/niome_hotkey2/result/2026-09-23T09:35:06/submission.json) |
| h3 | 2026-09-23 09:11:50,609 | 2026-09-23 08:35:29,155 | `5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf` | [data/inst/niome_hotkey3/result/2026-09-23T09:12:05](../../data/inst/niome_hotkey3/result/2026-09-23T09:12:05/submission.json) |
| h4 | 2026-09-23 09:23:22,397 | 2026-09-23 08:35:28,017 | `5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2` | [data/inst/niome_hotkey4/result/2026-09-23T09:23:36](../../data/inst/niome_hotkey4/result/2026-09-23T09:23:36/submission.json) |
| h5 | 2026-09-23 09:33:51,661 | 2026-09-23 08:35:29,090 | `5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn` | [data/inst/niome_hotkey5/result/2026-09-23T09:34:06](../../data/inst/niome_hotkey5/result/2026-09-23T09:34:06/submission.json) |
| h6 | 2026-09-23 09:02:31,103 | 2026-09-23 08:35:29,482 | `5H5v45M2i6cFtS3Di4abh2zjWrPmuJZzNPoAYtVAXVFa6sp3` | [data/inst/niome_hotkey6/result/2026-09-23T09:02:51](../../data/inst/niome_hotkey6/result/2026-09-23T09:02:51/submission.json) |
