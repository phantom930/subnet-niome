"""Render verified uploaded coverage and official scores for task ba385efa."""
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
    "All ten hotkeys uploaded **250 valid rows** each. Nine used prepared conjunctions with "
    "**100 Cas12a + 150 Cas9**, each with an **11-seed HDR band**. "
    "**h7 uploaded an ordinary fallback with 76 Cas12a + 174 Cas9.** Its prepared conjunction finished "
    "after upload, so the late 11-seed band in its window record is not submitted coverage.", "",
    f"[Official score feed]({d['source_api']}), fetched {d['fetched_at']}. "
    f"{d['score_records']} records, {d['unique_miners']} miners, {len(d['validators'])} validator. "
    f"Score timestamp: {top['created_at']} UTC. Audit generated {d['generated_at']}.", "",
    "## Official scores", "",
    f"| Miner | UID | Rank / {d['unique_miners']} | Final score | Weighted score | Consistency | Fidelity | Valid rows | Reward weight |",
    "|---|---:|---:|---:|---:|---:|---:|---:|---:|", score_line("Top miner", 1, top),
]
lines += [score_line(h, r["rank"], r["score"]) for h, r in ours.items()]
lines += ["", "Ranks count strictly higher scores; ties share rank. Identities were matched using full public "
          "hotkeys. All ten have zero reward weight in this snapshot. Final scores and breakdown factors "
          "are means across three seeds; multiplying mean factors need not reproduce the mean final score.", "",
          "## Joined seed windows and cut search", "",
          "The round plan generated **08:17:19 UTC** predicted width-100 windows **100 / 200 / 300**, "
          "joined as **100–399**. The actual conjunction layout instead uses fixed windows "
          "**100 / 200 / 500**, joined as **100–299 ∪ 500–599**: 300 seeds shared by all loops. "
          "Both band candidates and narrow cut search use this fixed shared set. The plan does not steer "
          "these band candidates. The plan and shared set both contain scoring seeds **263 and 269**, "
          "while **486** is outside. The actual width-100 window labels are **200 / 400 / 200**.", "",
          "h0–h9 were assigned loops 1–10, each excluding earlier loops' band seeds. "
          "The uploaded conjunctions are loops **1–7 and 9–10**. h7's assigned loop 8 completed too late "
          "and its uploaded fallback has no conjunction candidate or cut-search window.", "",
          "| Hotkey | Uploaded construction | Uploaded loop | Band candidates and cut search | Clean in cut search | Full clean set /900 |",
          "|---|---|---:|---|---:|---:|"]
for h, r in ours.items():
    conj = r["construction"] == "conjunction"
    lines.append(f"| {h} | {'conjunction' if conj else 'ordinary fallback'} | {r['band_loop'] or '—'} | "
                 f"{r['candidate_window_ranges'] if conj else 'not applicable'} | "
                 f"{str(len(r['clean_seeds_in_cut_space'])) + '/300' if conj else '—'} | {len(r['clean_seeds'])} |")
lines += ["", "The full clean counts exceed the logged narrow-window counts because replay also finds "
          "some clean seeds outside the 300 construction seeds. Both counts require every uploaded row "
          "to cut; Cas12a-only clean seeds outside that window are excluded where any Cas9 row fails to cut.", "",
          "## Uploaded HDR bands and scoring-seed hits", "",
          "**Clean:** all 250 rows cut. **HDR band:** all 250 rows return HDR. Both sets were recovered "
          "from actual archived submissions by replay over all seeds 100–999. A clean hit can have "
          "mixed HDR/NHEJ outcomes. Candidate membership alone guarantees neither clean nor HDR.", "",
          "| Hotkey | Clean /900 | Exact uploaded HDR band | Clean scoring-seed hits | HDR scoring-seed hits |",
          "|---|---:|---|---|---|"]
for h, r in ours.items():
    lines.append(f"| {h} | {len(r['clean_seeds'])} | {nums(r['band_seeds'])} | {nums(r['clean_hits'])} | {nums(r['band_hits'])} |")
lines += ["", f"The uploaded fleet union contains **{len(d['fleet_clean_union'])} distinct clean seeds** "
          f"and **{len(d['fleet_band_union'])} distinct HDR-band seeds**. "
          f"Of the clean union, **{len(set(d['fleet_clean_union']) & set(d['shared_joined']))} are inside "
          f"the shared 300-seed window** and **{len(set(d['fleet_clean_union']) - set(d['shared_joined']))} are outside**. "
          "All nine uploaded bands are mutually disjoint. Fleet clean hits: " + nums(d["fleet_clean_hits"]) + "; "
          "HDR hits: " + nums(d["fleet_band_hits"]) + ".", "",
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
          f"{t['distribution_fidelity_factor']:.6f}**. The largest factor gap is consistency.", ""]
for s in best["band_hits"]:
    sr = best["per_seed"][str(s)]
    lines += [f"{best_h} returned HDR on all 250 rows at **{s}**, scoring **{sr['final_score']:.6f}** at that seed. "
              "Its final score averages this result with the other two seed scores shown above.", ""]
lines += [f"Leader hotkey: `{top['miner_hotkey']}`. The leader's submitted rows, clean sets, HDR bands and "
          "per-seed outcomes were unavailable in the inspected data; its aggregate score does not identify those sets.", "",
          "## h7 fallback and late prepared band", "",
          "h7 received the task at **08:31:31 UTC** while prefetch was still building. After waiting, "
          "at **08:34:35** the runtime reported the prepared round unusable (still building at 187 seconds). "
          "It skipped hedges because the prepared build still held the GPU and generated an ordinary "
          "submission, which uploaded at **08:34:40.849**. The prepared conjunction finished only at "
          "**08:37:06.650**, with **99/300 clean seeds and an 11-seed band**. Those are preparation "
          "statistics, not statistics of h7's upload.", "",
          "Late, unsubmitted h7 band: **" + nums(ours["h7"]["latest_record_band"]) + "**. "
          f"The actual fallback has **{len(ours['h7']['clean_seeds'])} clean seeds over 100–999** and "
          f"**{len(ours['h7']['band_seeds'])} all-HDR band seeds**. "
          "Its window record combines `source: all_hdr_not_attempted` with `all_hdr_built: true` and the late "
          "band. This audit uses the upload archive and timing to distinguish them.", "",
          "## Exact uploaded clean sets", "",
          "Each list covers all 100–999 and includes the corresponding uploaded HDR band. "
          "The full-row definition excludes seeds where only the Cas12a group cuts cleanly.", ""]
for h, r in ours.items():
    lines += [f"### {h}: {len(r['clean_seeds'])} clean seeds", "", nums(r["clean_seeds"]) + ".", ""]
    extra = sorted(set(r["cas12a_clean_seeds"]) - set(r["clean_seeds"]))
    if extra:
        lines += [f"Cas12a-only clean seeds excluded from the full-row set: {nums(extra)}.", ""]
lines += ["## Verification and provenance", "",
          "- All ten successful upload archives identify this task. Contracts and references match the published task after substituting its seed stamp.",
          "- All 2,500 rows passed validation with CD34+_HSPC accessibility 0.87. All ten local three-seed final scores match official scores within 1e-10.",
          "- Generated all ten validations under this report directory without modifying miner archives.",
          "- Replayed **900 × 250 × 10 = 2,250,000 row/seed combinations**. Every outcome and indel length at the three scoring seeds matches local stage-3 detail.",
          "- For nine conjunctions, Cas12a and full-row clean counts inside the cut-search window agree with build logs; recovered uploaded bands match the window records. h7's late prepared band is excluded.",
          "- [Analysis JSON](analysis.json), [CSV](summary.csv), [official scores](scores.json), [task](task.json), [manifest](manifest.json), [sanitized build evidence](build_evidence.json), [layout evidence](layout_evidence.json).",
          "- Reproduce with `.venv/bin/python reports/ba385efa/validate.py`, then `.venv/bin/python reports/ba385efa/audit.py`, then `.venv/bin/python reports/ba385efa/render_report.py`. Validator and submission SHA-256 hashes are recorded in the analysis.", "",
          "| Hotkey | Upload UTC | Prepared build completed UTC | Prepared build uploaded | Public hotkey | Submission archive |",
          "|---|---|---|---|---|---|"]
for h, r in ours.items():
    lines.append(f"| {h} | {r['uploaded_at']} | {r['prepared_build_completed_at']} | {r['prepared_build_uploaded']} | "
                 f"`{r['hotkey']}` | [{r['archive']}](../../{r['archive']}/submission.json) |")
(OUT / "report.md").write_text("\n".join(lines) + "\n")
with (OUT / "summary.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["Hotkey", "UID", "Rank", "Official score", "Uploaded construction", "Valid rows", "Band candidate window", "Uploaded loop",
                     "Cut window", "Clean count /900", "Clean count in cut space", "Exact clean set", "Clean hits", "HDR band", "HDR hits",
                     "Latest record matches uploaded band", "Percent of leader", "Gap to leader", "Reward weight"])
    for h, r in ours.items():
        writer.writerow([h, r["uid"], r["rank"], r["score"]["final_score"], r["construction"], r["rows"],
                         r["candidate_window_ranges"] or "not applicable", r["band_loop"], r["cut_search_ranges"] or "not applicable",
                         len(r["clean_seeds"]), len(r["clean_seeds_in_cut_space"]) if r["cut_search_space"] else "not applicable",
                         nums(r["clean_seeds"]), nums(r["clean_hits"]), nums(r["band_seeds"]), nums(r["band_hits"]),
                         r["latest_record_band_matches_upload"], r["percent_of_top"], r["gap_to_top"], r["score"]["weight"]])
print(OUT / "report.md")
