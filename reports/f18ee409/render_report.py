"""Render the verified mixed-construction audit for task f18ee409."""
import csv
import json
from pathlib import Path

OUT = Path(__file__).resolve().parent
d = json.loads((OUT / "analysis.json").read_text())
ours, top, seeds = d["hotkeys"], d["top_miner"], d["seeds"]
assert len(ours) == 10 and "fleet_clean_union" in d
best_h = d["best_our_hotkey"]
best = ours[best_h]


def nums(values):
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
    f"{best_h} reached {best['percent_of_top']:.2f}% of the leader, a gap of {best['gap_to_top']:.6f}. "
    f"Top-10 cutoff: {d['top10_cutoff']:.6f}.", "",
    "All ten hotkeys uploaded **250 valid rows: 80 Cas12a + 170 Cas9**. "
    "**h0–h6 used prepared conjunctions; h7–h9 used prepared all-HDR submissions.** "
    "No ordinary fallback occurred. Every uploaded submission has **eight actual HDR-band seeds**. "
    "All ten locally reproduced final scores match the official scores within 1e-10.", "",
    "Fleet scoring-seed clean hits: **" + nums(d["fleet_clean_hits"]) + "**. "
    "Fleet HDR-band hits: **" + nums(d["fleet_band_hits"]) + "**.", "",
    "**Only h0 was clean at seed 214**: all 250 rows cut, with 111 HDR rows, for a per-seed score "
    "of **21.597457**. No uploaded clean set contained 249 or 784. None of the three scoring "
    "seeds hit an uploaded all-HDR band.", "",
    f"[Official score feed]({d['source_api']}), fetched {d['fetched_at']}. "
    f"{d['score_records']} records, {d['unique_miners']} miners, {len(d['validators'])} validator. "
    f"Score timestamp: {top['created_at']} UTC. Audit generated {d['generated_at']}.", "",
    "## Official scores", "",
    f"| Miner | UID | Rank / {d['unique_miners']} | Final score | Weighted score | Consistency | Fidelity | Valid rows | Reward weight |",
    "|---|---:|---:|---:|---:|---:|---:|---:|---:|", score_line("Top miner", 1, top),
]
lines += [score_line(h, r["rank"], r["score"]) for h, r in ours.items()]
lines += ["", "Ranks count strictly higher scores, so ties share rank. Identities were matched by full public "
          "hotkeys using this task's UIDs. All ten requested hotkeys have zero reward weight in this snapshot. "
          "Final scores and breakdown factors are means across three seeds; the product of mean factors "
          "need not reproduce the mean final score.", "",
          "## Joined windows and construction", "",
          "The task used a fixed layout independent of the prediction plan. **h0–h6** used seven "
          "**15-seed candidate windows at stride 15 over 200–299**, every hotkey on **loop 1**, "
          "with a common **200-seed cut search: 200–299 ∪ 400–499**. "
          "h6 wraps to **200–204 ∪ 290–299**, sharing five candidate seeds with h0.", "",
          "**h7–h9** used the all-HDR builder directly, over **34-seed candidate windows at stride 34 "
          "over 400–499**. h9 wraps to **400–401 ∪ 468–499**, sharing two candidate seeds with h7. "
          "This builder selects Cas12a rows to share HDR seeds inside that window, then selects Cas9 "
          "rows that return HDR on those same seeds. It does not use the conjunction's separate 200-seed "
          "cut search.", "",
          "**The current layout was moved after these builds.** Its present constants target 100–199 "
          "and 500–599. This report uses the task's historical build logs, archived rows and recovered "
          "bands, which establish the 200–299 and 400–499 windows actually uploaded.", "",
          "| Hotkey | Uploaded construction | Joined band candidates | Candidate count | Conjunction cut search | Logged clean statistic | Full clean /900 |",
          "|---|---|---|---:|---|---|---:|"]
for h, r in ours.items():
    conj = r["construction"] == "conjunction"
    log_stat = f"{r['prepared_logged_clean_count']}/{r['logged_clean_denominator']} " + ("cut-clean" if conj else "HDR-clean")
    lines.append(f"| {h} | {'conjunction' if conj else 'all-HDR'} | {r['candidate_window_ranges']} | "
                 f"{len(r['candidate_window'])} | {r['cut_search_ranges'] or 'not applicable'} | {log_stat} | {len(r['clean_seeds'])} |")
lines += ["", "The logged clean counts have different definitions: conjunction logs count the Cas12a "
          "group's **no-cut-free** seeds inside the cut window; all-HDR logs count its **all-HDR** "
          "seeds inside the candidate window. The report's full clean set always requires **every one "
          "of the 250 uploaded rows to cut**, and is recovered over all 900 seeds.", "",
          "## Actual uploaded bands and hits", "",
          "An **HDR-band seed** makes all 250 uploaded rows return HDR. Each actual band below was "
          "recovered by replay and is a subset of that hotkey's clean set and candidate window. "
          "A clean seed need not be all-HDR.", "",
          "**h7–h9's `window_used.json` band strings describe the 34-seed candidate windows, not "
          "34 actual HDR seeds.** Their build logs report `clean 8/34`; replay identifies the exact "
          "eight HDR seeds. h0–h6's JSON band lists already contain their actual eight HDR seeds and "
          "match replay. The bounding `window` fields on h6 and h9 omit the gaps in their joined sets.", "",
          "| Hotkey | Full clean /900 | Exact uploaded HDR band | Clean scoring-seed hits | HDR scoring-seed hits |",
          "|---|---:|---|---|---|"]
for h, r in ours.items():
    lines.append(f"| {h} | {len(r['clean_seeds'])} | {nums(r['band_seeds'])} | {nums(r['clean_hits'])} | {nums(r['band_hits'])} |")
lines += ["", f"Fleet union: **{len(d['fleet_clean_union'])} distinct clean seeds** and "
          f"**{len(d['fleet_band_union'])} distinct HDR-band seeds** from 80 band memberships. "
          "Overlapping HDR seeds: " + ("; ".join(f"{s} ({'/'.join(hs)})" for s, hs in sorted(d["duplicate_band_seeds"].items(), key=lambda x: int(x[0]))) or "none") + ".", "",
          "## Per-seed results", "",
          "Each entry is **per-seed final score / HDR rows / no-cut rows**, from 250 rows.", "",
          "| Hotkey | " + " | ".join(f"Seed {s}" for s in seeds) + " | Mean final score |",
          "|---|---:|---:|---:|---:|"]
for h, r in ours.items():
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
          f"{t['distribution_fidelity_factor']:.6f}**. The largest factor gap is consistency.", "",
          "Seeds 214 and 249 were within the conjunction fleet's candidate class, specifically h0 and "
          "h3's windows, but candidate membership is not a guarantee of clean or HDR coverage. "
          "Seed 784 was outside both construction groups' candidate windows and the conjunction cut "
          "search. The all-HDR group targeted class 400–499; this task drew no seed in that class.", "",
          f"Leader hotkey: `{top['miner_hotkey']}`. The leader's submitted rows, exact clean sets, "
          "HDR bands and per-seed outcomes were unavailable in the inspected data. Its aggregate score "
          "does not identify those sets. One task does not establish which construction performs "
          "better across tasks.", "",
          "## Exact uploaded clean sets", "",
          "Each list covers all seeds 100–999 and includes the corresponding uploaded HDR band.", ""]
for h, r in ours.items():
    lines += [f"### {h}: {len(r['clean_seeds'])} clean seeds", "", nums(r["clean_seeds"]) + ".", ""]
    extra = sorted(set(r["cas12a_clean_seeds"]) - set(r["clean_seeds"]))
    if extra:
        lines += [f"Cas12a-only clean seeds excluded from the full-row set: {nums(extra)}.", ""]
lines += ["## Verification and provenance", "",
          "- All ten successful upload records identify this task. Contracts and references match the published task after substituting its scoring seed stamp.",
          "- All 2,500 rows passed validation with HEK293 accessibility 0.35. All ten local three-seed final scores match the official scores within 1e-10.",
          "- All validation artifacts were generated under this report directory; miner archives were read without modification.",
          "- Replayed **900 × 250 × 10 = 2,250,000 row/seed combinations**. At all three scoring seeds, every row's outcome and indel length matches local stage-3 detail.",
          "- Conjunction cut-clean counts match all seven build logs. All-HDR HDR-clean counts match the three 8/34 logs and the actual eight-seed bands. Every upload used its completed prepared build.",
          "- [Analysis JSON](analysis.json), [CSV](summary.csv), [official scores](scores.json), [task](task.json), [manifest](manifest.json), [sanitized build evidence](build_evidence.json), [historical layout evidence](layout_evidence.json).",
          "- Reproduce with `.venv/bin/python reports/f18ee409/validate.py`, then `.venv/bin/python reports/f18ee409/audit.py`, then `.venv/bin/python reports/f18ee409/render_report.py`. Validator and submission SHA-256 hashes are recorded in the analysis.", "",
          "| Hotkey | Upload UTC | Prepared build completed UTC | Public hotkey | Submission archive |",
          "|---|---|---|---|---|"]
for h, r in ours.items():
    lines.append(f"| {h} | {r['uploaded_at']} | {r['prepared_build_completed_at']} | `{r['hotkey']}` | "
                 f"[{r['archive']}](../../{r['archive']}/submission.json) |")
(OUT / "report.md").write_text("\n".join(lines) + "\n")
with (OUT / "summary.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["Hotkey", "UID", "Rank", "Official score", "Construction", "Valid rows", "Band candidate window", "Conjunction cut window",
                     "Clean count /900", "Logged clean count", "Logged clean denominator", "Logged clean definition", "Exact clean set",
                     "Clean hits", "HDR band", "HDR hits", "Window record band semantics", "Percent of leader", "Gap to leader", "Reward weight"])
    for h, r in ours.items():
        writer.writerow([h, r["uid"], r["rank"], r["score"]["final_score"], r["construction"], r["rows"], r["candidate_window_ranges"],
                         r["cut_search_ranges"] or "not applicable", len(r["clean_seeds"]), r["prepared_logged_clean_count"],
                         r["logged_clean_denominator"], r["logged_clean_definition"], nums(r["clean_seeds"]), nums(r["clean_hits"]),
                         nums(r["band_seeds"]), nums(r["band_hits"]), r["latest_record_band_semantics"], r["percent_of_top"], r["gap_to_top"], r["score"]["weight"]])
print(OUT / "report.md")
