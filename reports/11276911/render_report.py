"""Render task 11276911's audited uploads and recovered seed sets."""
import csv
import json
from pathlib import Path

OUT = Path(__file__).resolve().parent
d = json.loads((OUT / "analysis.json").read_text())
evidence = json.loads((OUT / "build_evidence.json").read_text())
ours, top, seeds = d["hotkeys"], d["top_miner"], d["seeds"]
best_h = max(ours, key=lambda h: ours[h]["score"]["final_score"])
best = ours[best_h]
fallbacks = [h for h, r in ours.items() if r["construction"] == "ordinary_hdr_fallback"]
conjunctions = [h for h in ours if h not in fallbacks]
assert fallbacks == ["h0", "h6"]


def nums(values):
    return ", ".join(map(str, values)) or "none"


def score_line(label, rank, score):
    b = score["breakdown"]
    return (f"| {label} | {score['miner_uid']} | {rank} | {score['final_score']:.6f} | "
            f"{b['total_weighted_score']:.6f} | {b['consistency_factor']:.6f} | "
            f"{b['distribution_fidelity_factor']:.6f} | {score['weight']:.6f} |")


def event_time(h, text):
    matches = [line for line in evidence[h]["events"] if text in line]
    assert len(matches) == 1, (h, text, matches)
    return matches[0][:23].replace(",", ".")


lines = [
    f"# Submission audit: {d['task_id']}", "",
    f"**{d['cell_type']}**, created **{d['task_created_at']} UTC**. Scoring seeds: **{nums(seeds)}**.", "",
    f"**Best of h0–h6: {best_h}, {best['score']['final_score']:.6f}, rank {best['rank']}/{d['unique_miners']}. "
    f"Leader: UID {top['miner_uid']}, {top['final_score']:.6f}.** "
    f"{best_h} reached {best['percent_of_top']:.2f}% of the leader, with a gap of {best['gap_to_top']:.6f}. "
    f"The top-10 cutoff was {d['top10_cutoff']:.6f}, {d['top10_cutoff'] - best['score']['final_score']:.6f} above {best_h}. "
    "All seven uploaded successfully, have 250 valid rows, and received zero reward weight.", "",
    "**h0 and h6 submitted identical ordinary HDR-construction fallbacks.** Their prepared conjunctions "
    "finished after their wait deadlines and uploads. h1–h5 submitted their prepared conjunction rows. "
    "The prepared-band records for h0/h6 do not describe their actual uploads.", "",
    f"Scores: [task-specific public API]({d['source_api']}), fetched {d['fetched_at']}; "
    f"{d['score_records']} records, {d['unique_miners']} distinct miners, {len(d['validators'])} validator. "
    f"Leader's score timestamp: {top['created_at']} UTC. Audit generated {d['generated_at']}.", "",
    "## Official score comparison", "",
    f"| Miner | UID | Rank / {d['unique_miners']} | Final score | Weighted score | Consistency | Fidelity | Reward weight |",
    "|---|---:|---:|---:|---:|---:|---:|---:|",
    score_line("Top miner", 1, top),
]
lines += [score_line(h, r["rank"], r["score"]) for h, r in ours.items()]
lines += ["", "Ranks are calculated from official scores using competition ranking: exact ties share rank. "
          f"h0 and h6 both rank **{ours['h0']['rank']}**; they occupy sorted positions "
          f"{ours['h0']['sorted_api_position']} and {ours['h6']['sorted_api_position']} when API order is retained within the tie. "
          "The submitted files have identical SHA-256 hashes.", "",
          f"Leader hotkey: `{top['miner_hotkey']}`.", "",
          "Final score is the mean of the three per-seed scores. For these seven submissions, weighted score "
          "and fidelity are constant across seeds, so final = weighted score × mean consistency × fidelity.", "",
          "## Joined windows, clean sets and actual bands", "",
          "The plan generated **2026-09-18 05:17:11 UTC** selected joined band space **300–399 ∪ 600–799**, "
          "or 100-wide windows **300/600/700**. Conjunction cut search covered **100–999 (900 seeds)**. "
          "The HDR candidate slices were 150 seeds wide, at offsets **0/50/100/150/200/250/25** for h0–h6. "
          "Each of the five submitted conjunctions has **100 Cas12a + 150 Cas9 rows** and an **11-seed HDR band**. "
          "The two ordinary fallbacks each have **74 Cas12a + 176 Cas9 rows**.", "",
          "**Clean set:** every uploaded row cuts at that seed. **HDR band:** every uploaded row returns HDR. "
          "Both are measured from the exact uploaded rows over all seeds 100–999. For submitted conjunctions, "
          "the full-row clean set equals the Cas12a group clean set and its size agrees with the build log. "
          "Their recovered HDR bands exactly match the original band records.", "",
          "| Hotkey | Uploaded construction | Planned candidate slice | Offset | Clean /900 | Clean inside joined /300 | Clean outside joined |",
          "|---|---|---|---:|---:|---:|---:|"]
for h, r in ours.items():
    fallback = h in fallbacks
    label = "ordinary fallback" if fallback else "conjunction"
    candidate = r["candidate_window_ranges"] + (" (unused)" if fallback else "")
    lines.append(f"| {h} | {label} | {candidate} | {r['candidate_offset']} | {len(r['clean_seeds'])} | "
                 f"{len(r['clean_seeds_in_joined_space'])} | {len(r['clean_seeds_outside_joined_space'])} |")
lines += ["", "| Hotkey | Actual uploaded HDR band seeds | Clean scoring-seed hits | HDR scoring-seed hits |",
          "|---|---|---|---|"]
for h, r in ours.items():
    lines.append(f"| {h} | {nums(r['band_seeds'])} | {nums(r['clean_hits'])} | {nums(r['band_hits'])} |")
lines += ["", "Candidate slices are reconstructed from the logged joined space and configured offsets. "
          "All submitted conjunction bands lie inside the corresponding candidate slice.", "",
          f"Fleet union: **{len(d['fleet_clean_union'])}/900 clean seeds**, **{len(d['fleet_band_union'])} distinct HDR-band seeds** "
          f"({sum(len(r['band_seeds']) for r in ours.values())} band memberships before overlap). "
          f"Clean scoring-seed hits: **{nums(d['fleet_clean_hits'])}**; HDR-band hits: **{nums(d['fleet_band_hits'])}**.", "",
          "All three scoring seeds **448/931/493** are outside the predicted joined band space. "
          "The table above identifies any all-cut hits recovered by the wider cut search.", "",
          "## Per-seed scores and outcomes", "",
          "Each cell is **per-seed final score / no-cut rows**. Zero no-cut rows means an all-cut clean hit.", "",
          "| Hotkey | " + " | ".join(f"Seed {s}" for s in seeds) + " | Mean final score |",
          "|---|---:|---:|---:|---:|"]
for h, r in ours.items():
    values = [f"{r['per_seed'][str(s)]['final_score']:.6f} / {r['per_seed'][str(s)]['outcomes'].get('no_cut', 0)}"
              for s in seeds]
    lines.append(f"| {h} | " + " | ".join(values) + f" | {r['score']['final_score']:.6f} |")
b, t = best["score"]["breakdown"], top["breakdown"]
lines += ["", "## Comparison with the leader", "",
          f"{best_h}'s weighted score is **{b['total_weighted_score']:.6f}**, compared with the leader's "
          f"**{t['total_weighted_score']:.6f}**. Fidelity is **{b['distribution_fidelity_factor']:.6f} versus "
          f"{t['distribution_fidelity_factor']:.6f}**. The main difference is mean consistency: "
          f"**{b['consistency_factor']:.6f} versus {t['consistency_factor']:.6f}**. "
          "The leader's clean set, HDR band and per-seed outcomes cannot be recovered from its public aggregate score; "
          "its submitted rows were not available in the inspected data.", "",
          "## h0 and h6 fallback provenance", "",
          "| Hotkey | Prepared wait expired (UTC) | Ordinary upload succeeded (UTC) | Prepared conjunction finished (UTC) | Unused prepared clean count | Uploaded clean count |",
          "|---|---|---|---|---:|---:|"]
for h in fallbacks:
    r = ours[h]
    lines.append(f"| {h} | {event_time(h, 'prepared round unusable')} | {event_time(h, 'Submitted 250 rows')} | "
                 f"{event_time(h, 'Build: conjunction ')} | {r['logged_prepared_clean_count']} | {len(r['clean_seeds'])} |")
lines += ["", "Both handlers waited about 179–180 seconds for their prepared builds. When those waits expired, "
          "the GPU was still occupied, so they used ordinary construction. Both uploads succeeded with about 101 seconds "
          "of URL TTL remaining. Their prepared conjunctions completed later, so the completed bands were never submitted.", "",
          "h5 also received a request before its build finished, but its wait lasted only about 36 seconds; "
          "it successfully uploaded the prepared conjunction with 247 seconds of TTL remaining.", "",
          "The final `window_used.json` records for h0/h6 combine `all_hdr_not_attempted` with `all_hdr_built: true` "
          "and bands written by the later prepared builds. Uploaded-row replay and exact official score agreement establish "
          "which rows were actually submitted.", "",
          "| Hotkey | Recorded prepared band — not submitted |",
          "|---|---|"]
for h in fallbacks:
    lines.append(f"| {h} | {nums(ours[h]['recorded_prepared_band'])} |")
lines += ["", "## Exact uploaded clean sets", "",
          "These lists include HDR-band seeds and cover all seeds 100–999. "
          "Full arrays are also saved in [analysis.json](analysis.json).", ""]
for h, r in ours.items():
    lines += [f"### {h}: {len(r['clean_seeds'])} clean seeds", "", nums(r["clean_seeds"]) + ".", ""]
lines += ["## Verification and artifacts", "",
          "- All seven upload records confirm this task; contracts and reference files match the published task after the scoring-seed stamp.",
          "- All seven submissions were scored through calc.py into this report directory. Each has 250/250 valid rows, "
          "cell accessibility 0.87, and a three-seed score matching the official API exactly.",
          "- Stage-3 outcomes and indel lengths agree with the full local validation for all 250 rows at all three real scoring seeds.",
          "- Recovered actual clean sets and HDR bands with **900 × 250 × 7 = 1,575,000 row simulations**.",
          "- [Official scores](scores.json), [published task](task.json), [sanitized logs](build_evidence.json), "
          "[archive manifest](manifest.json), [summary CSV](summary.csv), [full analysis](analysis.json).",
          "- Reproduce from the repository root with `.venv/bin/python reports/11276911/validate.py`, "
          "then `.venv/bin/python reports/11276911/audit.py`, then `.venv/bin/python reports/11276911/render_report.py`.", "",
          "Uploaded archives:", ""]
for h, r in ours.items():
    lines.append(f"- **{h}** `{r['hotkey']}` — [{r['archive']}](../../{r['archive']}/submission.json); "
                 f"SHA-256 `{r['submission_sha256']}`.")
(OUT / "report.md").write_text("\n".join(lines) + "\n")
with (OUT / "summary.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["Hotkey", "UID", "Rank (ties shared)", "Score", "Construction", "Planned candidate slice",
                     "Clean count /900", "Clean hits", "HDR band", "HDR hits", "Percent of leader", "Gap to leader"])
    for h, r in ours.items():
        writer.writerow([h, r["uid"], r["rank"], r["score"]["final_score"], r["construction"],
                         r["candidate_window_ranges"], len(r["clean_seeds"]), nums(r["clean_hits"]),
                         nums(r["band_seeds"]), nums(r["band_hits"]), r["percent_of_top"], r["gap_to_top"]])
print(OUT / "report.md")
