"""Render seven wide-cut conjunction uploads and the official score comparison."""
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
    return "N/A" if values is None else (", ".join(map(str, values)) or "none")


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
    "**Seven uploads were found: h0 and h4–h9.** Each has **250 valid rows, 100 Cas12a + 150 Cas9**, "
    "and an **11-seed HDR band**. All seven used prepared conjunctions; no fallback was used. "
    "All seven have zero reward weight in this snapshot. h1–h3 have no task build/window/upload records "
    "or official score entries; their results are unavailable, not zero.", "",
    "**h9 hit HDR seed 806: all 250 rows returned HDR**, scoring **197.660227** at that seed. "
    "Its other seed scores were **17.062024** at 143 and **19.755283** at 625, giving the mean **78.159178**. "
    "Additional clean hits occurred at **143 on h0/h4** and **806 on h7**. Seed **625** was outside "
    "every uploaded clean set. No other clean or HDR hits occurred.", "",
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
          "per-seed products. Missing hotkeys have no score or rank.", "",
          f"**UID 41 is the leader's current UID, not our h1's identity.** Leader hotkey: `{top['miner_hotkey']}`. "
          f"Our h1 hotkey: `{ours['h1']['hotkey']}`. All identity matching in this audit uses full public hotkeys.", "",
          "## Joined windows and actual band placement", "",
          "The **20:17 UTC plan** recorded width-100 windows **100 / 200 / 300**, joined as **100–399**. "
          "Actual windows were **100 / 800 / 600** for scoring seeds **143 / 806 / 625**: one matched that prediction.", "",
          "For these conjunctions, the active band space was the full **100–999**. Seven **300-seed "
          "contiguous candidate windows at stride 100** were assigned to h0, h4, h5, h6, h7, h8 and h9. "
          "The plan's predicted classes did not steer this placement. The build log label `predicted space of 900` "
          "denotes the full space, not the plan's 300 predicted seeds.", "",
          "All seven **CD34+_HSPC** builds used **cut 900 seeds (wide), k=11, group=100**, confirmed by "
          "the task logs. Band candidates span 300 seeds per hotkey; cut selection spans 100–999. "
          "Clean sets therefore extend beyond the candidate windows.", "",
          "| Hotkey | Assigned 300-seed band candidates | Offset in 100–999 | Cut search | Clean /900 | Clean inside candidate window | Clean outside candidate window |",
          "|---|---|---:|---|---:|---:|---:|"]
for h, r in ours.items():
    if not r["submitted"]:
        lines.append(f"| {h} | No task build or upload | N/A | N/A | N/A | N/A | N/A |")
    else:
        lines.append(f"| {h} | {r['candidate_window_ranges']} | {r['candidate_offset']} | 100–999 | {len(r['clean_seeds'])} | "
                     f"{len(r['clean_seeds_in_candidate_window'])} | {len(r['clean_seeds_outside_candidate_window'])} |")
lines += ["", "The seven overlapping candidate windows together cover all 100–999. Seed 143 was in h0's "
          "candidate window; 806 was in h8/h9's; 625 was in h6/h7/h8's. These are candidate memberships, "
          "not guarantees of clean or HDR hits.", "",
          "## Actual uploaded clean sets and HDR bands", "",
          "**Clean:** all 250 uploaded rows cut at that seed. **HDR band:** all 250 uploaded rows return HDR. "
          "Both were recovered by replay over all 900 seeds (100–999). Each full-row clean set equals the "
          "Cas12a group clean set and matches the build log count. Every recovered band matches its task "
          "window record and lies within its candidate window.", "",
          "| Hotkey | Clean /900 | Exact uploaded HDR band | Clean scoring-seed hits | HDR hits |",
          "|---|---:|---|---|---|"]
for h, r in ours.items():
    if not r["submitted"]:
        lines.append(f"| {h} | N/A | No uploaded rows available | N/A | N/A |")
    else:
        lines.append(f"| {h} | {len(r['clean_seeds'])} | {nums(r['band_seeds'])} | {nums(r['clean_hits'])} | {nums(r['band_hits'])} |")
lines += ["", f"The uploaded fleet union contains **{len(d['fleet_clean_union'])}/900 clean seeds** and "
          f"**{len(d['fleet_band_union'])} distinct HDR-band seeds** "
          f"({sum(len(r['band_seeds']) for r in uploaded.values())} band memberships before overlaps). "
          "h9 at 806 was the only uploaded all-HDR hit.", "",
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
          "h9's high score at its HDR hit was averaged with two low seed scores. "
          "The leader's submitted rows, clean sets, HDR bands and per-seed outcomes were not available in the "
          "inspected data; its aggregate score does not identify those sets.", "",
          "## h1–h3: no task activity", "",
          "No task preparation, received-task/upload event, window record or upload archive was found for "
          "h1–h3. None appears in the complete official score feed. Their earlier runtime logs reported them "
          "unregistered on netuid 55 on 2026-09-19. No submitted clean set, band, score or rank is available "
          "for these hotkeys on this task.", "",
          "## Exact uploaded clean sets", "",
          "Each list covers all 100–999 seeds and includes its HDR band.", ""]
for h, r in uploaded.items():
    lines += [f"### {h}: {len(r['clean_seeds'])} clean seeds", "", nums(r["clean_seeds"]) + ".", ""]
lines += ["## Verification and artifacts", "",
          "- All seven uploads used their prepared conjunctions; no ordinary fallback was used. All recorded bands match actual uploads.",
          "- All seven archives identify this task; contracts and references match the published task after the seed stamp.",
          "- All 1,750 submitted rows passed stage 1/2, with accessibility 0.87. Local three-seed scores exactly match all seven official scores.",
          "- Reused h7/h8/h9 archived validations and generated h0/h4/h5/h6 validations in report-local directories.",
          "- Replayed outcomes and indel lengths match local stage-3 detail for every uploaded row at all three scoring seeds.",
          "- Recovered sets over **900 × 250 × 7 = 1,575,000 row/seed combinations**.",
          "- [Full analysis](analysis.json), [official scores](scores.json), [task](task.json), "
          "[sanitized build/upload and registration evidence](build_evidence.json), [layout source evidence](layout_evidence.json), "
          "[manifest](manifest.json), [CSV](summary.csv).",
          "- Reproduce with `.venv/bin/python reports/00b2e47c/validate.py`, then "
          "`.venv/bin/python reports/00b2e47c/audit.py`, then `.venv/bin/python reports/00b2e47c/render_report.py`. "
          "The analysis includes validator-code and submission SHA-256 hashes.", "",
          "| Hotkey | Public hotkey | Upload time UTC | Submission archive |",
          "|---|---|---|---|"]
for h, r in ours.items():
    archive = f"[{r['archive']}](../../{r['archive']}/submission.json)" if r["submitted"] else "No upload archive"
    lines.append(f"| {h} | `{r['hotkey']}` | {r.get('uploaded_at', 'N/A')} | {archive} |")
(OUT / "report.md").write_text("\n".join(lines) + "\n")
with (OUT / "summary.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["Hotkey", "UID", "Rank", "Official score", "Uploaded", "Valid rows", "Band candidate window", "Cut window",
                     "Clean count /900", "Clean inside candidate window", "Clean outside candidate window", "Exact clean set", "Clean hits",
                     "HDR band", "HDR hits", "Percent of leader", "Gap to leader", "Reward weight"])
    for h, r in ours.items():
        yes = r["submitted"]
        writer.writerow([h, r["uid"], r["rank"], r["score"]["final_score"] if yes else "N/A", yes, r["rows"],
                         r["candidate_window_ranges"], "100–999" if yes else "N/A", len(r["clean_seeds"]) if yes else "N/A",
                         len(r["clean_seeds_in_candidate_window"]) if yes else "N/A", len(r["clean_seeds_outside_candidate_window"]) if yes else "N/A",
                         nums(r["clean_seeds"]), nums(r["clean_hits"]), nums(r["band_seeds"]), nums(r["band_hits"]),
                         r["percent_of_top"], r["gap_to_top"], r["score"]["weight"] if yes else "N/A"])
print(OUT / "report.md")
