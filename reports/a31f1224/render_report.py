"""Render the a31f1224 audit, including distinct builds and missing h9 upload evidence."""
import csv
import json
from pathlib import Path

OUT = Path(__file__).resolve().parent
d = json.loads((OUT / "analysis.json").read_text())
ours, top, seeds = d["hotkeys"], d["top_miner"], d["seeds"]
uploaded = {h: r for h, r in ours.items() if r["submitted"]}
best_h = max(uploaded, key=lambda h: uploaded[h]["score"]["final_score"])
best = uploaded[best_h]


def nums(values):
    if values is None:
        return "not recoverable"
    return ", ".join(map(str, values)) or "none"


def score_line(label, rank, score):
    b = score["breakdown"]
    return (f"| {label} | {score['miner_uid']} | {rank} | {score['final_score']:.6f} | "
            f"{b['total_weighted_score']:.6f} | {b['consistency_factor']:.6f} | "
            f"{b['distribution_fidelity_factor']:.6f} | {b['n_valid_experiments']} | {score['weight']:.6f} |")


labels = {"conjunction": "conjunction", "ordinary_hdr_fallback": "ordinary fallback",
          "all_hdr": "all-HDR", "no_uploaded_archive": "no upload archive"}
lines = [
    f"# Submission audit: {d['task_id']}", "",
    f"**{d['cell_type']}**, created **{d['task_created_at']} UTC**. Scoring seeds: **{nums(seeds)}**.", "",
    f"**Best of h0–h9: {best_h}, {best['score']['final_score']:.6f}, rank {best['rank']}/{d['unique_miners']}. "
    f"Leader: UID {top['miner_uid']}, {top['final_score']:.6f}.** "
    f"{best_h} achieved {best['percent_of_top']:.2f}% of the top score, a gap of {best['gap_to_top']:.6f}. "
    f"The top-10 cutoff was {d['top10_cutoff']:.6f}, {d['top10_cutoff']-best['score']['final_score']:.6f} above {best_h}.", "",
    "**Nine successful uploads were found, h0–h8, each with 250 valid rows. h9 has no task upload archive or "
    "received-task/upload event in its log; its official score is 0 with 0 valid rows.** "
    "All ten hotkeys have zero reward weight. A prepared h9 build exists in the log, but its clean set and band "
    "must not be represented as submitted.", "",
    "h5 uploaded an **ordinary fallback**; h7 uploaded an **on-demand all-HDR** build; "
    "h0–h4, h6 and h8 uploaded **conjunctions**. Later rebuilds overwrote several per-task band records. "
    "All actual clean sets and bands below come from the uploaded archives.", "",
    f"[Official task score feed]({d['source_api']}), fetched {d['fetched_at']}. "
    f"{d['score_records']} score records, {d['unique_miners']} miners, {len(d['validators'])} validator. "
    f"Leader score timestamp: {top['created_at']} UTC. Audit generated {d['generated_at']}.", "",
    "## Official score comparison", "",
    f"| Miner | UID | Rank / {d['unique_miners']} | Final score | Weighted score | Consistency | Fidelity | Valid rows | Reward weight |",
    "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    score_line("Top miner", 1, top),
]
lines += [score_line(h, r["rank"], r["score"]) for h, r in ours.items()]
lines += ["", "Ranks use competition ranking: 1 plus the count of strictly higher scores. "
          "h9 shares the zero-score rank with other zero-score miners; its displayed rank does not imply a valid submission.", "",
          f"Leader hotkey: `{top['miner_hotkey']}`.", "",
          "## Original joined windows and candidate slices", "",
          "**h0–h7:** the original 14:17 plan supplied **200–299 ∪ 400–499 ∪ 800–899** "
          "(windows **200/400/800**). Seeds **287 and 486** are inside this joined space; 998 is outside. "
          "h5 attempted this plan but ultimately uploaded ordinary fallback rows.", "",
          "**h8/h9:** their logs show the built-in fallback joined layout **100–199 ∪ 300–399 ∪ 800–899**, "
          "source `env_pin_no_entry_joined`. All three scoring seeds are outside this space. "
          "h8 uploaded its conjunction; h9 only has preparation evidence.", "",
          "Submitted conjunctions searched for clean seeds over the entire **100–999** range, while their HDR "
          "bands came from 150-seed slices of the joined band space. Their rows are **100 Cas12a + 150 Cas9**. "
          "h7's all-HDR build used the full joined 300-seed candidate space and the same Cas mix. "
          "h5's ordinary fallback is **76 Cas12a + 174 Cas9**.", "",
          "| Hotkey | Uploaded construction | Joined band space used/planned | Band candidate slice | Offset | Clean /900 |",
          "|---|---|---|---|---:|---:|"]
for h, r in ours.items():
    candidate = r["candidate_window_ranges"]
    if h == "h5":
        candidate += " (unused)"
    elif h == "h8":
        candidate += " (reconstructed)"
    elif h == "h9":
        candidate += " (prepared only; reconstructed)"
    count = str(len(r["clean_seeds"])) if r["clean_seeds"] is not None else "N/A"
    offset = r["candidate_offset"] if r["candidate_offset"] is not None else "full 300"
    lines.append(f"| {h} | {labels[r['construction']]} | {r['joined_band_ranges']} | {candidate} | {offset} | {count} |")
lines += ["", "h0–h5's original offsets were **0/50/100/150/200/250**, with **25 for h6**, "
          "as supported by the historical stride-50 plan logs. h7's uploaded all-HDR build used all 300 candidates. "
          "h8/h9's 240/270 offsets are reconstructed under the new ten-hotkey stride-30 layout; "
          "their exact candidate slices were not directly logged. h8's recovered band is consistent with that slice. "
          "The exact uploaded bands below are independently established by replay, regardless of that reconstruction.", "",
          "## Actual clean sets and HDR bands", "",
          "**Clean set:** all 250 uploaded rows cut. **HDR band:** all 250 uploaded rows return HDR. "
          "Both sets were scanned over seeds 100–999. An HDR-band seed is also clean. "
          "All seven uploaded conjunctions have exactly 11 HDR seeds; their full-row clean sets equal their "
          "Cas12a clean sets and match the original build-log counts.", "",
          "| Hotkey | Clean /900 | Clean inside joined /300 | Clean outside joined | Actual uploaded HDR band seeds | Clean hits | HDR hits |",
          "|---|---:|---:|---:|---|---|---|"]
for h, r in ours.items():
    if not r["submitted"]:
        lines.append(f"| {h} | N/A | N/A | N/A | No uploaded rows available | N/A | N/A |")
        continue
    lines.append(f"| {h} | {len(r['clean_seeds'])} | {len(r['clean_seeds_in_joined_space'])} | "
                 f"{len(r['clean_seeds_outside_joined_space'])} | {nums(r['band_seeds'])} | "
                 f"{nums(r['clean_hits'])} | {nums(r['band_hits'])} |")
lines += ["", "**h4 hits clean seed 998; h8 hits clean seed 486. No other clean hits and no HDR-band hits occur.** "
          "Both clean hits are outside those hotkeys' joined band spaces, demonstrating that the wider cut search "
          "provided coverage outside the band prediction. Seed 287 is absent from every uploaded clean set.", "",
          f"Across the **nine actual uploads**, the clean union contains **{len(d['fleet_clean_union'])}/900** seeds "
          f"and the HDR-band union contains **{len(d['fleet_band_union'])} distinct seeds** "
          f"({sum(len(r['band_seeds']) for r in uploaded.values())} band memberships before overlap). "
          "h9's prepared build is excluded from these unions.", "",
          "## Per-seed results", "",
          "Each cell is **per-seed final score / no-cut rows**. Zero no-cut rows denotes a clean hit.", "",
          "| Hotkey | " + " | ".join(f"Seed {s}" for s in seeds) + " | Mean score |",
          "|---|---:|---:|---:|---:|"]
for h, r in uploaded.items():
    values = [f"{r['per_seed'][str(s)]['final_score']:.6f} / {r['per_seed'][str(s)]['outcomes'].get('no_cut', 0)}" for s in seeds]
    lines.append(f"| {h} | " + " | ".join(values) + f" | {r['score']['final_score']:.6f} |")
b, t = best["score"]["breakdown"], top["breakdown"]
lines += ["", "## Comparison with the leader", "",
          f"{best_h}'s weighted score (**{b['total_weighted_score']:.6f}**) and fidelity "
          f"(**{b['distribution_fidelity_factor']:.6f}**) are above the leader's "
          f"**{t['total_weighted_score']:.6f}** and **{t['distribution_fidelity_factor']:.6f}**. "
          f"The decisive gap is mean consistency: **{b['consistency_factor']:.6f} versus {t['consistency_factor']:.6f}**. "
          "The leader's submitted rows, clean set, bands and individual seed outcomes were not available in the "
          "inspected data and cannot be reconstructed from its aggregate score.", "",
          "## Build and upload provenance", "",
          "- **h0–h4 and h6:** prepared conjunctions started around 15:09 and finished around 15:21 UTC, "
          "taking 711–719 seconds. Later validator requests used those prepared rows.",
          "- **h5:** the validator handler waited from 15:16:29 to 15:19:28, then started ordinary construction "
          "because the prepared build still held the GPU. It uploaded the fallback at **15:19:32**, with 101 seconds "
          "of TTL remaining. Its prepared conjunction finished at 15:21:05 and was not submitted.",
          "- **h7:** at 15:06:45 its prepared cache was for another task, so it built this task on demand. "
          "It completed an **all-HDR** build and uploaded at **15:09:45**, with 103 seconds of TTL remaining. "
          "Its later prepared conjunction was not the upload.",
          "- **h8:** it prepared against the fallback joined layout at 15:11:30, finished at 15:21:17, "
          "and uploaded that prepared conjunction at 15:44:40.",
          "- **h9:** its log records a prepared conjunction from 15:11:30 to 15:21:14, but no received-task or "
          "successful upload event was found for this task. There is no matching upload archive. "
          "The API records **0 valid rows and score 0**. The reason the validator received no valid rows is not established by these records.", "",
          "### Later rebuilds and overwritten band records", "",
          "At approximately **15:50–15:56 UTC**, h0–h7 rebuilt the same task under the newer stride-30 layout. "
          "This happened after every archived upload. The current `window_used.json` entries for these hotkeys "
          "therefore refer to later builds; h0 retained the same band because its offset remained zero. "
          "h1–h7's recorded bands differ from their actual uploaded bands. h8's record still matches its upload.", "",
          "| Hotkey | Latest task-record time (UTC) | Latest recorded band matches actual upload? |",
          "|---|---|---|"]
for h, r in uploaded.items():
    lines.append(f"| {h} | {r['latest_record_at']} | {'yes' if r['latest_record_band_matches_upload'] else 'NO'} |")
lines += ["", "### h9 preparation only — not submitted evidence", "",
          f"h9's prepared build logged **{ours['h9']['logged_prepared_clean_count']}/900 clean seeds** and the recorded band "
          f"**{nums(ours['h9']['latest_record_band'])}**. These describe preparation only. "
          "No uploaded h9 rows are available to recover its actual submitted clean set or band, so the main tables use N/A.", "",
          "## Exact uploaded clean sets", "",
          "The lists include HDR-band seeds and cover the full 100–999 seed space. "
          "There is no reconstructable submitted clean set for h9.", ""]
for h, r in uploaded.items():
    lines += [f"### {h}: {len(r['clean_seeds'])} clean seeds", "", nums(r["clean_seeds"]) + ".", ""]
lines += ["## Verification and artifacts", "",
          "- All nine archived uploads identify this task. Their contracts and references match the public task after accounting for its later seed stamp.",
          "- All nine submissions contain 250 stage-1/2 valid rows, with cell accessibility 0.87, and their local three-seed scores match the API exactly.",
          "- Saved validations were reused for h0–h3, h7 and h8; h4–h6 were validated with calc.py into this report directory.",
          "- Stage-3 outcomes and indel lengths match saved/regenerated detail at all three scoring seeds for every uploaded row.",
          "- Recovered sets by evaluating **900 × 250 × 9 = 2,025,000 row/seed combinations**; h9 was not simulated without uploaded rows.",
          "- [Full analysis](analysis.json), [official scores](scores.json), [task](task.json), [sanitized miner logs](build_evidence.json), "
          "[historical layout logs](plan_history.json), [archive manifest](manifest.json), [summary CSV](summary.csv).",
          "- Reproduce from the repository root with `.venv/bin/python reports/a31f1224/validate.py`, "
          "`.venv/bin/python reports/a31f1224/audit.py`, and `.venv/bin/python reports/a31f1224/render_report.py`.", "",
          "Uploaded archives:", ""]
for h, r in uploaded.items():
    lines.append(f"- **{h}** `{r['hotkey']}` — [{r['archive']}](../../{r['archive']}/submission.json); "
                 f"SHA-256 `{r['submission_sha256']}`.")
lines += [f"- **h9** `{ours['h9']['hotkey']}` — no upload archive for this task.", ""]
(OUT / "report.md").write_text("\n".join(lines))
with (OUT / "summary.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["Hotkey", "UID", "Rank (ties shared)", "Final score", "Valid rows", "Construction",
                     "Joined space", "Candidate slice", "Clean count /900", "Clean hits", "Actual HDR band",
                     "HDR hits", "Percent of top", "Gap to top"])
    for h, r in ours.items():
        writer.writerow([h, r["uid"], r["rank"], r["score"]["final_score"], r["score"]["breakdown"]["n_valid_experiments"],
                         r["construction"], r["joined_band_ranges"], r["candidate_window_ranges"],
                         len(r["clean_seeds"]) if r["clean_seeds"] is not None else "N/A", nums(r["clean_hits"]),
                         nums(r["band_seeds"]), nums(r["band_hits"]), r["percent_of_top"], r["gap_to_top"]])
print(OUT / "report.md")
