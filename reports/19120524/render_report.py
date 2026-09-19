"""Render the six uploaded submissions and four preparation-only records for 19120524."""
import csv
import json
from pathlib import Path

OUT = Path(__file__).resolve().parent
d = json.loads((OUT / "analysis.json").read_text())
ours, top, seeds = d["hotkeys"], d["top_miner"], d["seeds"]
uploaded = {h: r for h, r in ours.items() if r["submitted"]}
missing = {h: r for h, r in ours.items() if not r["submitted"]}
best_h = max(uploaded, key=lambda h: uploaded[h]["score"]["final_score"])
best = uploaded[best_h]
assert set(uploaded) == {"h1", "h2", "h3", "h7", "h8", "h9"}


def nums(values):
    if values is None:
        return "N/A"
    return ", ".join(map(str, values)) or "none"


def score_line(label, rank, score):
    b = score["breakdown"]
    return (f"| {label} | {score['miner_uid']} | {rank} | {score['final_score']:.6f} | "
            f"{b['total_weighted_score']:.6f} | {b['consistency_factor']:.6f} | "
            f"{b['distribution_fidelity_factor']:.6f} | {b['n_valid_experiments']} | {score['weight']:.6f} |")


lines = [
    f"# Submission audit: {d['task_id']}", "",
    f"**{d['cell_type']}**, created **{d['task_created_at']} UTC**. Scoring seeds: **{nums(seeds)}**.", "",
    f"**Best of h0–h9: {best_h}, {best['score']['final_score']:.6f}, rank {best['rank']}/{d['unique_miners']}. "
    f"Leader: UID {top['miner_uid']}, {top['final_score']:.6f}.** "
    f"{best_h} achieved {best['percent_of_top']:.2f}% of the leader, a gap of {best['gap_to_top']:.6f}. "
    f"Top-10 cutoff: {d['top10_cutoff']:.6f}; {best_h} was {d['top10_cutoff']-best['score']['final_score']:.6f} below it.", "",
    "**Six uploaded conjunctions were found: h1, h2, h3, h7, h8 and h9. Each has 250 valid rows. "
    "h0, h4, h5 and h6 have no matching upload archive or received-task/upload event; each has official score 0 and 0 valid rows.** "
    "All ten hotkeys have zero reward weight. Prepared builds for the four missing uploads are documented separately below.", "",
    f"[Official task score feed]({d['source_api']}), fetched {d['fetched_at']}. "
    f"{d['score_records']} records, {d['unique_miners']} miners, {len(d['validators'])} validator. "
    f"Leader score timestamp: {top['created_at']} UTC. Audit generated {d['generated_at']}.", "",
    "## Official scores", "",
    f"| Miner | UID | Rank / {d['unique_miners']} | Final score | Weighted score | Consistency | Fidelity | Valid rows | Reward weight |",
    "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    score_line("Top miner", 1, top),
]
lines += [score_line(h, r["rank"], r["score"]) for h, r in ours.items()]
lines += ["", "Rank is 1 plus the number of strictly higher scores. All four zero-score hotkeys share rank "
          f"**{ours['h0']['rank']}** with other zero-score miners; that rank does not imply a valid upload.", "",
          f"Leader hotkey: `{top['miner_hotkey']}`.", "",
          "## Original joined windows and candidate slices", "",
          "The original **16:40 plan** selected windows **100/200/300**, joined as **100–399**. "
          "All ten hotkeys prepared against that band space. Seeds **236 and 156** are inside it; **596** is outside. "
          "The submitted builds used the original ten-hotkey **stride-30** layout, with 150-seed circular candidate slices. "
          "Every submitted conjunction used a **900-seed cut search (100–999)**, **100 Cas12a + 150 Cas9 rows**, "
          "and an **11-seed HDR band**.", "",
          "| Hotkey | Upload evidence | Band candidate slice | Offset | Actual clean /900 |",
          "|---|---|---|---:|---:|"]
for h, r in ours.items():
    status = "uploaded conjunction" if r["submitted"] else "prepared only; no upload found"
    count = str(len(r["clean_seeds"])) if r["submitted"] else "N/A"
    lines.append(f"| {h} | {status} | {r['candidate_window_ranges']} | {r['candidate_offset']} | {count} |")
lines += ["", "Candidate slices are reconstructed from the historical plan's offsets and the logged 100–399 band space. "
          "Every actual uploaded band lies inside its corresponding original slice. "
          "The later split between predicted and complementary windows was not the layout used for these uploads.", "",
          "## Actual uploaded clean sets and HDR bands", "",
          "**Clean set:** all 250 uploaded rows cut at that seed. **HDR band:** all 250 uploaded rows return HDR. "
          "Both were recovered over seeds 100–999. The full-row clean sets equal the Cas12a group clean sets, "
          "and their sizes match the original conjunction build logs.", "",
          "| Hotkey | Clean /900 | Clean inside joined /300 | Clean outside joined | Actual uploaded HDR band seeds | Clean hits | HDR hits |",
          "|---|---:|---:|---:|---|---|---|"]
for h, r in ours.items():
    if not r["submitted"]:
        lines.append(f"| {h} | N/A | N/A | N/A | No uploaded rows available | N/A | N/A |")
        continue
    lines.append(f"| {h} | {len(r['clean_seeds'])} | {len(r['clean_seeds_in_joined_space'])} | "
                 f"{len(r['clean_seeds_outside_joined_space'])} | {nums(r['band_seeds'])} | "
                 f"{nums(r['clean_hits'])} | {nums(r['band_hits'])} |")
lines += ["", "**h1 hit clean seed 236; h2 and h8 hit clean seed 156. No other uploaded clean hits or HDR-band hits occurred.** "
          "Seed 596 was outside all six uploaded clean sets. "
          f"The six-upload fleet union contains **{len(d['fleet_clean_union'])}/900 clean seeds** and "
          f"**{len(d['fleet_band_union'])} distinct HDR-band seeds** "
          f"({sum(len(r['band_seeds']) for r in uploaded.values())} band memberships before overlap). "
          "Preparation-only records are excluded from these totals.", "",
          "## Per-seed results", "",
          "Each entry is **per-seed final score / no-cut rows**. Zero no-cut rows identifies a clean hit.", "",
          "| Hotkey | " + " | ".join(f"Seed {s}" for s in seeds) + " | Mean final score |",
          "|---|---:|---:|---:|---:|"]
for h, r in uploaded.items():
    values = [f"{r['per_seed'][str(s)]['final_score']:.6f} / {r['per_seed'][str(s)]['outcomes'].get('no_cut', 0)}" for s in seeds]
    lines.append(f"| {h} | " + " | ".join(values) + f" | {r['score']['final_score']:.6f} |")
b, t = best["score"]["breakdown"], top["breakdown"]
lines += ["", "## Comparison with the leader", "",
          f"{best_h}'s weighted score is **{b['total_weighted_score']:.6f}**, versus the leader's "
          f"**{t['total_weighted_score']:.6f}**. Fidelity is **{b['distribution_fidelity_factor']:.6f} versus "
          f"{t['distribution_fidelity_factor']:.6f}**. The largest gap is mean consistency: "
          f"**{b['consistency_factor']:.6f} versus {t['consistency_factor']:.6f}**. "
          "h8's clean hit at 156 scored **60.475469**, but its other seed scores of **27.862710 at 236** "
          "and **34.670339 at 596** reduced the three-seed mean to **41.002839**. "
          "The leader's submitted rows, clean set, bands and individual seed outcomes were not available in the inspected data; "
          "aggregate scores do not identify those sets.", "",
          "## Build and upload provenance", "",
          "The logs show preparation attempts around **17:30**, **17:34**, and **17:47 UTC**. "
          "The completed 17:47 builds were used for all six successful uploads, which occurred between "
          "**17:50:23 and 18:23:34 UTC**, with 264–289 seconds of URL TTL remaining. "
          "No ordinary or all-HDR fallback was used by these six submissions.", "",
          "At **18:47–18:49 UTC**, those same six hotkeys rebuilt the task after their uploads. "
          "The newer configuration changed h1–h3's offsets and moved h7–h9 into complementary band space. "
          "Their latest `window_used.json` bands therefore differ from all six actual uploaded bands. "
          "The report uses uploaded-row replay to avoid crediting these later bands.", "",
          "| Hotkey | Actual uploaded band | Later recorded band — not submitted |",
          "|---|---|---|"]
for h, r in uploaded.items():
    assert not r["latest_record_band_matches_upload"]
    lines.append(f"| {h} | {nums(r['band_seeds'])} | {nums(r['latest_record_band'])} |")
lines += ["", "## h0/h4/h5/h6: preparation only", "",
          "These four hotkeys have completed prepared builds in their logs, but no matching upload archives or "
          "received-task/successful-upload events were found. The API gives each 0 valid experiments and score 0. "
          "The reason the validator obtained no valid rows is not established by the available records.", "",
          "The following are **logged prepared-build values**, not recovered submitted sets:", "",
          "| Hotkey | Logged prepared clean count /900 | Recorded prepared band |",
          "|---|---:|---|"]
for h, r in missing.items():
    lines.append(f"| {h} | {r['logged_prepared_clean_count']} | {nums(r['latest_record_band'])} |")
lines += ["", "h4's prepared-band record includes scoring seed **236**, but there is no uploaded submission "
          "available to validate or credit that potential band hit. Its official score remains zero.", "",
          "## Exact uploaded clean sets", "",
          "These lists include each uploaded HDR band and cover the full 100–999 seed space. "
          "No submitted clean set can be recovered for h0/h4/h5/h6.", ""]
for h, r in uploaded.items():
    lines += [f"### {h}: {len(r['clean_seeds'])} clean seeds", "", nums(r["clean_seeds"]) + ".", ""]
lines += ["## Verification and artifacts", "",
          "- All six upload archives identify this task, with contracts and references matching the published task after the seed stamp.",
          "- All six submissions have 250 stage-1/2 valid rows and accessibility 0.82. Archived local three-seed scores match the API exactly.",
          "- Replayed outcomes and indel lengths match archived stage-3 detail for every uploaded row at all three scoring seeds.",
          "- Recovered clean sets and bands over **900 × 250 × 6 = 1,350,000 row/seed combinations**.",
          "- [Full analysis](analysis.json), [official scores](scores.json), [task](task.json), [sanitized build/upload logs](build_evidence.json), "
          "[historical layout](plan_history.json), [archive manifest](manifest.json), [summary CSV](summary.csv).",
          "- Reproduce with `.venv/bin/python reports/19120524/audit.py` and "
          "`.venv/bin/python reports/19120524/render_report.py` from the repository root, using the archived validations.", "",
          "Uploaded archives:", ""]
for h, r in uploaded.items():
    lines.append(f"- **{h}** `{r['hotkey']}` — [{r['archive']}](../../{r['archive']}/submission.json); "
                 f"SHA-256 `{r['submission_sha256']}`.")
for h, r in missing.items():
    lines.append(f"- **{h}** `{r['hotkey']}` — no upload archive for this task.")
(OUT / "report.md").write_text("\n".join(lines) + "\n")
with (OUT / "summary.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["Hotkey", "UID", "Rank (ties shared)", "Final score", "Valid rows", "Uploaded",
                     "Original joined space", "Original candidate slice", "Clean count /900", "Clean hits",
                     "Actual HDR band", "HDR hits", "Percent of leader", "Gap to leader"])
    for h, r in ours.items():
        writer.writerow([h, r["uid"], r["rank"], r["score"]["final_score"], r["score"]["breakdown"]["n_valid_experiments"],
                         bool(r["submitted"]), r["joined_band_ranges"], r["candidate_window_ranges"],
                         len(r["clean_seeds"]) if r["submitted"] else "N/A", nums(r["clean_hits"]), nums(r["band_seeds"]),
                         nums(r["band_hits"]), r["percent_of_top"], r["gap_to_top"]])
print(OUT / "report.md")
