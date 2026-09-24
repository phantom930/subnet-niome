"""Render seven uploaded sets and three prepared-only records for af6b11ee."""
import csv
import json
from pathlib import Path

OUT = Path(__file__).resolve().parent
d = json.loads((OUT / "analysis.json").read_text())
evidence = json.loads((OUT / "build_evidence.json").read_text())
ours, top, seeds = d["hotkeys"], d["top_miner"], d["seeds"]
uploaded = {h: r for h, r in ours.items() if r["submitted"]}
missing = {h: r for h, r in ours.items() if not r["submitted"]}
best_h = max(uploaded, key=lambda h: uploaded[h]["score"]["final_score"])
best = uploaded[best_h]


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
    "**Seven uploads were found: h0 and h4–h9.** Each has 250 valid rows, **100 Cas12a + 150 Cas9**, "
    "and an **11-seed HDR band**. All seven have zero reward weight in the official snapshot. "
    "**h1–h3 have no task upload archive, received-task/upload event or official score entry.** Their results "
    "are unavailable, not zero. Prepared builds for these hotkeys are documented separately.", "",
    "**h8 hit HDR seed 353: all 250 uploaded rows returned HDR**, scoring **235.918078** at that seed. "
    "Its scores at 218 and 143 were **21.320120** and **22.058984**, giving the task mean **93.099061**. "
    "h7 hit clean seed **218**, scoring **33.948359** there. No other uploaded clean or HDR hits occurred; "
    "143 was outside every uploaded clean set.", "",
    f"[Official task score feed]({d['source_api']}), fetched {d['fetched_at']}. "
    f"{d['score_records']} records, {d['unique_miners']} miners, {len(d['validators'])} validator. "
    f"Score timestamp: {top['created_at']} UTC. Audit generated {d['generated_at']}.", "",
    "## Official scores", "",
    f"| Miner | UID | Rank / {d['unique_miners']} | Final score | Weighted score | Consistency | Fidelity | Valid rows | Reward weight |",
    "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    score_line("Top miner", 1, top),
]
for h, r in ours.items():
    lines.append(score_line(h, r["rank"], r["score"]) if r["submitted"] else f"| {h} | N/A | N/A | No score entry | N/A | N/A | N/A | N/A | N/A |")
lines += ["", "Rank is 1 plus the number of strictly higher scores; ties share rank. Final score and breakdown "
          "factors are means over three seeds. Multiplying mean factors is not generally equivalent to averaging "
          "per-seed products. Missing hotkeys have no rank or score; they are not assigned a zero-score rank.", "",
          f"Leader hotkey: `{top['miner_hotkey']}`.", "",
          "## Historical joined windows and candidate slices", "",
          "The **03:17 UTC plan** selected width-100 windows **100 / 200 / 300**, joined as **100–399**. "
          "Actual windows were **300 / 200 / 100** for seeds **353 / 218 / 143**: **all three matched**.", "",
          "The task's original build logs show **ten disjoint 30-seed candidate slices**, at offsets "
          "0, 30, …, 270 within the predicted 300 seeds. Every prepared conjunction used **k=11, group=100** "
          "and a full **100–999 cut search**. All seven uploaded bands lie inside their corresponding slices. "
          "This report uses the task's logged layout; the repository's later seven-hotkey width-300 layout "
          "was not used for these submissions.", "",
          "| Hotkey | Upload status | Assigned 30-seed candidate slice | Offset | Actual clean /900 |",
          "|---|---|---|---:|---:|"]
for h, r in ours.items():
    status = "uploaded conjunction" if r["submitted"] else "prepared only; no upload found"
    count = len(r["clean_seeds"]) if r["submitted"] else "N/A"
    lines.append(f"| {h} | {status} | {r['candidate_window_ranges']} | {r['candidate_offset']} | {count} |")
lines += ["", "The seven uploaded candidate slices cover **100–129 ∪ 220–399** (210 seeds). "
          "The assigned slices **130–219** belong to h1–h3, for which no upload exists. Scoring seeds "
          "**143 and 218** fell into those unsubmitted slices; **353** fell into h8's uploaded slice. "
          "Clean sets can extend beyond a candidate slice because the cut search spans all 900 seeds.", "",
          "## Actual uploaded clean sets and HDR bands", "",
          "**Clean:** all 250 uploaded rows cut at that seed. **HDR band:** all 250 uploaded rows return HDR. "
          "Both sets were recovered over 100–999. For all seven uploads, the full-row clean set equals "
          "the Cas12a group clean set, its count matches the build log, and the HDR band matches the window record.", "",
          "| Hotkey | Clean /900 | Exact uploaded HDR band | Clean scoring-seed hits | HDR hits |",
          "|---|---:|---|---|---|"]
for h, r in ours.items():
    if not r["submitted"]:
        lines.append(f"| {h} | N/A | No uploaded rows available | N/A | N/A |")
    else:
        lines.append(f"| {h} | {len(r['clean_seeds'])} | {nums(r['band_seeds'])} | {nums(r['clean_hits'])} | {nums(r['band_hits'])} |")
lines += ["", f"The uploaded fleet union contains **{len(d['fleet_clean_union'])}/900 clean seeds** and "
          f"**{len(d['fleet_band_union'])} distinct HDR-band seeds**. Its seven bands are disjoint "
          "(7 × 11 = 77 seeds). All totals exclude h1–h3's prepared builds.", "",
          "## Per-seed results", "",
          "Each entry is **per-seed final score / HDR rows / no-cut rows**, from 250 rows per seed.", "",
          "| Hotkey | " + " | ".join(f"Seed {s}" for s in seeds) + " | Mean final score |",
          "|---|---:|---:|---:|---:|"]
for h, r in uploaded.items():
    values = []
    for s in seeds:
        sr = r["per_seed"][str(s)]
        values.append(f"{sr['final_score']:.6f} / {sr['outcomes'].get('HDR', 0)} / {sr['outcomes'].get('no_cut', 0)}")
    lines.append(f"| {h} | " + " | ".join(values) + f" | {r['score']['final_score']:.6f} |")
b, t = best["score"]["breakdown"], top["breakdown"]
lines += ["", "## Comparison with the leader", "",
          f"{best_h}'s weighted score was **{b['total_weighted_score']:.6f}**, versus the leader's "
          f"**{t['total_weighted_score']:.6f}**. Mean consistency was **{b['consistency_factor']:.6f} versus "
          f"{t['consistency_factor']:.6f}**, and fidelity was **{b['distribution_fidelity_factor']:.6f} versus "
          f"{t['distribution_fidelity_factor']:.6f}**. The largest factor gap is consistency. "
          "h8's strong score at its HDR hit was averaged with two low seed scores. "
          "The leader's submitted rows, clean sets, HDR bands and per-seed outcomes were not available in the "
          "inspected data; its aggregate score does not identify those sets.", "",
          "## h1–h3: prepared builds without uploads", "",
          "Runtime logs reported these three hotkeys unregistered on netuid 55 on **2026-09-19**, before this task. "
          "Their processes still prepared datasets for this task, but no received-task/upload event or matching "
          "upload archive was found, and none of their public hotkeys appears in this task's complete 248-record score feed.", "",
          "The following counts and bands are **prepared-build records, not reconstructed submitted sets**:", "",
          "| Hotkey | Unregistered log time UTC | Logged prepared clean /900 | Recorded prepared band | Scoring seeds in prepared band |",
          "|---|---|---:|---|---|"]
for h, r in missing.items():
    stamp = evidence[h]["registration_events"][-1][:23]
    lines.append(f"| {h} | {stamp} | {r['logged_prepared_clean_count']} | {nums(r['prepared_band_record'])} | {nums(r['prepared_band_scoring_seed_overlap'])} |")
lines += ["", "**h3's prepared band includes 218**, but it was not uploaded. It cannot be counted as an actual "
          "HDR hit or assigned a hypothetical official score. Exact submitted clean sets cannot be recovered "
          "for h1–h3 because there are no submitted rows to replay.", "",
          "## Exact uploaded clean sets", "",
          "Each list includes its HDR band and covers the entire 100–999 seed space.", ""]
for h, r in uploaded.items():
    lines += [f"### {h}: {len(r['clean_seeds'])} clean seeds", "", nums(r["clean_seeds"]) + ".", ""]
lines += ["## Build provenance and verification", "",
          "- All seven successful uploads used prepared conjunctions. No ordinary fallback appears in their task logs. "
          "h0 and h5 waited for their ongoing builds, then uploaded with 199 and 117 seconds of URL TTL remaining.",
          "- All seven archives identify this task; contracts and references match the published task after the seed stamp.",
          "- All 1,750 submitted rows passed stage 1/2, with accessibility 0.82. Local three-seed scores exactly match all seven official scores.",
          "- Reused h4's archived validation and generated the other six in report-local directories.",
          "- Replayed outcomes and indel lengths match local stage-3 detail for every uploaded row at all three scoring seeds.",
          "- Recovered sets over **900 × 250 × 7 = 1,575,000 row/seed combinations**.",
          "- [Full analysis](analysis.json), [official scores](scores.json), [task](task.json), "
          "[sanitized build/upload and registration evidence](build_evidence.json), [manifest](manifest.json), [CSV](summary.csv).",
          "- Reproduce with `.venv/bin/python reports/af6b11ee/validate.py`, then "
          "`.venv/bin/python reports/af6b11ee/audit.py`, then `.venv/bin/python reports/af6b11ee/render_report.py`. "
          "The analysis includes validator-code and submission SHA-256 hashes.", "",
          "| Hotkey | Public hotkey | Upload time UTC | Submission archive |",
          "|---|---|---|---|"]
for h, r in ours.items():
    archive = f"[{r['archive']}](../../{r['archive']}/submission.json)" if r["submitted"] else "No upload archive"
    lines.append(f"| {h} | `{r['hotkey']}` | {r.get('uploaded_at', 'N/A')} | {archive} |")
(OUT / "report.md").write_text("\n".join(lines) + "\n")
with (OUT / "summary.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["Hotkey", "UID", "Rank", "Official score", "Uploaded", "Valid rows", "Joined band space",
                     "Candidate slice", "Clean count /900", "Exact submitted clean set", "Clean hits", "Uploaded HDR band",
                     "HDR hits", "Percent of leader", "Gap to leader", "Prepared-only clean count", "Prepared-only band"])
    for h, r in ours.items():
        yes = r["submitted"]
        writer.writerow([h, r["uid"], r["rank"], r["score"]["final_score"] if yes else "N/A", yes,
                         r["rows"], r["band_space_ranges"], r["candidate_window_ranges"], len(r["clean_seeds"]) if yes else "N/A",
                         nums(r["clean_seeds"]), nums(r["clean_hits"]), nums(r["band_seeds"]), nums(r["band_hits"]),
                         r["percent_of_top"], r["gap_to_top"], r.get("logged_prepared_clean_count", ""),
                         nums(r["prepared_band_record"]) if not yes else ""])
print(OUT / "report.md")
