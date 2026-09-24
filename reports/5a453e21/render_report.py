"""Render official score comparison and recovered coverage, preserving the seed-stamp limitation."""
import csv
import json
from pathlib import Path

OUT = Path(__file__).resolve().parent
d = json.loads((OUT / "analysis.json").read_text())
ours, top = d["hotkeys"], d["top_miner"]
assert len(ours) == 7 and "fleet_clean_union" in d
best_h = d["best_our_hotkey"]
best = ours[best_h]


def nums(values):
    return "Unavailable" if values is None else (", ".join(map(str, values)) or "none")


def score_line(label, rank, score):
    b = score["breakdown"]
    return (f"| {label} | {score['miner_uid']} | {rank} | {score['final_score']:.6f} | "
            f"{b['total_weighted_score']:.6f} | {b['consistency_factor']:.6f} | "
            f"{b['distribution_fidelity_factor']:.6f} | {b['n_valid_experiments']} | {score['weight']:.8f} |")


lines = [
    f"# Submission audit: {d['task_id']}", "",
    f"**{d['cell_type']}**, created **{d['task_created_at']} UTC**.", "",
    f"**Best of h0–h6: {best_h}, {best['score']['final_score']:.6f}, rank {best['rank']}/{d['unique_miners']}. "
    f"Leader: UID {top['miner_uid']}, {top['final_score']:.6f}.** "
    f"{best_h} reached {best['percent_of_top']:.2f}% of the leader, a gap of {best['gap_to_top']:.6f}. "
    f"Top-10 cutoff: {d['top10_cutoff']:.6f}.", "",
    "All seven hotkeys uploaded prepared conjunctions with **250 valid rows: 100 Cas12a + 150 Cas9**. "
    "Each upload has an **11-seed HDR band**. No fallback occurred.", "",
    "**Scoring seeds are unavailable in the inspected sources.** The public task feed still carries "
    "`seed: 0`, although the official score feed has results. Fresh requests with pagination and different "
    "query parameters also returned 0. The public S3 contract returned HTTP 403. "
    "The 0 is an unstamped placeholder, not a reported scoring seed. **Clean/HDR hits on the actual "
    "scoring seeds and an exact reproduction of final scores cannot be confirmed from this snapshot.** "
    "They are reported as unavailable, not zero.", "",
    f"[Official score feed]({d['source_api']}), fetched {d['fetched_at']}. "
    f"{d['score_records']} records, {d['unique_miners']} miners, {len(d['validators'])} validator. "
    f"Score timestamp: {top['created_at']} UTC. Audit generated {d['generated_at']}.", "",
    "## Official scores", "",
    f"| Miner | UID | Rank / {d['unique_miners']} | Final score | Weighted score | Consistency | Fidelity | Valid rows | Reward weight |",
    "|---|---:|---:|---:|---:|---:|---:|---:|---:|", score_line("Top miner", 1, top),
]
lines += [score_line(h, r["rank"], r["score"]) for h, r in ours.items()]
lines += ["", "Ranks count strictly higher scores; ties share rank. All identities were matched by full public "
          "hotkeys, using the UIDs in this task's official score snapshot. h3 is the only requested hotkey "
          "with nonzero reward weight: **0.0098 (0.98%)**. The leader's weight is approximately **0.294 (29.4%)**. "
          "These are reported task weights, not a calculation of on-chain payout.", "",
          "## Joined seed windows and actual layout", "",
          "The task used a **fixed layout**, independent of the prediction plan. "
          "The fleet's band space is **200–299**, divided into seven circular windows of **15 seeds at "
          "stride 15**. Every hotkey uses **loop 1**. All seven share the **200-seed cut window "
          "200–299 ∪ 400–499**. The log text `shared window, loop 1 of 7` does not mean every hotkey "
          "uses the same band candidates: each has its own 15-seed slice.", "",
          "**h6 wraps:** its candidates are **200–204 ∪ 290–299**, not all 200–299. "
          "Its latest record's `[200,299]` field is only the bounding range. The build log and layout "
          "source preserve the exact joined window. h0 and h6 share five candidate seeds (200–204).", "",
          "| Hotkey | Joined band candidates | Offset | Loop | Clean in shared cut /200 | Full clean /900 | Exact uploaded HDR band |",
          "|---|---|---:|---:|---:|---:|---|"]
for h, r in ours.items():
    lines.append(f"| {h} | {r['candidate_window_ranges']} | {r['candidate_offset']} | {r['band_loop']} | "
                 f"{len(r['clean_seeds_in_cut_space'])} | {len(r['clean_seeds'])} | {nums(r['band_seeds'])} |")
lines += ["", "**Clean** means all 250 uploaded rows cut at that seed. **HDR band** means all 250 rows "
          "return HDR. Both sets were recovered by replay over every seed from 100 through 999. "
          "The build logs count clean seeds inside the 200-seed cut window; full replay can find additional "
          "clean seeds outside it. Cas12a-only cleanliness is checked separately and is not counted as "
          "full-row cleanliness when a Cas9 row fails to cut.", "",
          f"Fleet union: **{len(d['fleet_clean_union'])} distinct clean seeds** and "
          f"**{len(d['fleet_band_union'])} distinct HDR-band seeds** from 77 band memberships. "
          f"Clean union inside the shared cut window: **{len(set(d['fleet_clean_union']) & set(d['shared_cut_joined']))}/200**. "
          "Repeated HDR-band seeds: " + "; ".join(f"{s} ({'/'.join(hs)})" for s, hs in sorted(d["duplicate_band_seeds"].items(), key=lambda x: int(x[0]))) + ".", "",
          "## Comparison with the top miner", ""]
b, t = best["score"]["breakdown"], top["breakdown"]
lines += [f"{best_h}'s weighted score was **{b['total_weighted_score']:.6f}**, versus the leader's "
          f"**{t['total_weighted_score']:.6f}**. Mean consistency was **{b['consistency_factor']:.6f} versus "
          f"{t['consistency_factor']:.6f}**, and fidelity was **{b['distribution_fidelity_factor']:.6f} versus "
          f"{t['distribution_fidelity_factor']:.6f}**. Weighted scores are close; the largest factor gap "
          "is consistency. An aggregate final score alone does not establish which seeds hit a band. "
          "For multi-seed scoring, multiplying mean factors need not reproduce the mean final score.", "",
          f"Leader hotkey: `{top['miner_hotkey']}`. The leader's submission rows and exact clean/HDR sets "
          "were not available in the inspected data.", "",
          "## Exact uploaded clean sets", "",
          "Each list covers the entire 100–999 range and includes its uploaded HDR band.", ""]
for h, r in ours.items():
    lines += [f"### {h}: {len(r['clean_seeds'])} clean seeds", "", nums(r["clean_seeds"]) + ".", "",
              "Clean seeds outside the shared cut window: " + nums(sorted(set(r["clean_seeds"]) - set(r["cut_search_space"]))) + ".", ""]
    extra = sorted(set(r["cas12a_clean_seeds"]) - set(r["clean_seeds"]))
    if extra:
        lines += [f"Cas12a-only clean seeds excluded from the full-row set: {nums(extra)}.", ""]
lines += ["## Verification and provenance", "",
          "- All seven successful upload archives identify this task. Contracts and references match the task snapshot field for field.",
          "- All 1,750 rows pass stage 1/2 with HUDEP-2 accessibility 0.82, from the current public cell-types table.",
          "- Local weighted scores and distribution-fidelity factors match all seven official records within 1e-9. These quantities do not depend on the scoring seed.",
          "- Because the task seed stamp is unavailable, `validate.py` uses **diagnostic seed 100** to prepare validation features. Its local final scores are probes and are not presented as the official task scores.",
          "- Replayed **900 × 250 × 7 = 1,575,000 row/seed combinations**. Diagnostic-seed outcomes and indel lengths match the local pipeline. All seven recovered bands match the recorded bands and lie inside their exact candidate windows.",
          "- Full-row and Cas12a clean counts inside the shared cut window match all seven build log counts. All uploads occurred after their prepared builds completed.",
          "- [Analysis JSON](analysis.json), [CSV](summary.csv), [official scores](scores.json), [task snapshot](task.json), [seed availability checks](seed_availability.json), [manifest](manifest.json), [sanitized build evidence](build_evidence.json), [layout evidence](layout_evidence.json).",
          "- Reproduce coverage with `.venv/bin/python reports/5a453e21/validate.py`, then `.venv/bin/python reports/5a453e21/audit.py`, then `.venv/bin/python reports/5a453e21/render_report.py`. Validator and submission SHA-256 hashes are recorded in the analysis.", "",
          "| Hotkey | Upload UTC | Prepared build completed UTC | Public hotkey | Submission archive |",
          "|---|---|---|---|---|"]
for h, r in ours.items():
    lines.append(f"| {h} | {r['uploaded_at']} | {r['prepared_build_completed_at']} | `{r['hotkey']}` | "
                 f"[{r['archive']}](../../{r['archive']}/submission.json) |")
(OUT / "report.md").write_text("\n".join(lines) + "\n")
with (OUT / "summary.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["Hotkey", "UID", "Rank", "Official score", "Valid rows", "Band candidate window", "Loop", "Cut window",
                     "Clean count /900", "Clean count in cut space", "Exact clean set", "Clean scoring-seed hits", "HDR band",
                     "HDR scoring-seed hits", "Percent of leader", "Gap to leader", "Reward weight", "Scoring seed status"])
    for h, r in ours.items():
        writer.writerow([h, r["uid"], r["rank"], r["score"]["final_score"], r["rows"], r["candidate_window_ranges"], r["band_loop"],
                         r["cut_search_ranges"], len(r["clean_seeds"]), len(r["clean_seeds_in_cut_space"]), nums(r["clean_seeds"]),
                         nums(r["clean_hits"]), nums(r["band_seeds"]), nums(r["band_hits"]), r["percent_of_top"], r["gap_to_top"],
                         r["score"]["weight"], "Unavailable: public task seed is 0"])
print(OUT / "report.md")
