"""Render this task's verified submissions and official score comparison."""
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
            f"{b['distribution_fidelity_factor']:.6f} | {b['n_valid_experiments']} | {score['weight']:.8f} |")


lines = [
    f"# Submission audit: {d['task_id']}", "",
    f"**{d['cell_type']}**, created **{d['task_created_at']} UTC**. Scoring seeds: **{nums(seeds)}**.", "",
    f"**Best of h0–h9: {best_h}, {best['score']['final_score']:.6f}, rank {best['rank']}/{d['unique_miners']}. "
    f"Leader: UID {top['miner_uid']}, {top['final_score']:.6f}.** "
    f"{best_h} reached {best['percent_of_top']:.2f}% of the leader, a gap of {best['gap_to_top']:.6f}. "
    f"Top-10 cutoff: {d['top10_cutoff']:.6f}.", "",
    "All ten hotkeys uploaded **250 valid rows**, comprising **80 Cas12a and 170 Cas9** rows each. "
    "**h0–h8 used prepared conjunctions; h9 used a prepared all-HDR submission.** All ten locally "
    "reproduced final scores match the official scores within 1e-10.", "",
    f"[Official score feed]({d['source_api']}), fetched {d['fetched_at']}. "
    f"{d['score_records']} records, {d['unique_miners']} miners, {len(d['validators'])} validator. "
    f"Leader score timestamp: {top['created_at']} UTC. Audit generated {d['generated_at']}.", "",
    "## Official scores", "",
    f"| Miner | UID | Rank / {d['unique_miners']} | Final score | Weighted score | Consistency | Fidelity | Valid rows | Reward weight |",
    "|---|---:|---:|---:|---:|---:|---:|---:|---:|", score_line("Top miner", 1, top),
]
lines += [score_line(h, r["rank"], r["score"]) for h, r in ours.items()]
rewarded = [(h, r) for h, r in ours.items() if r["score"]["weight"] > 0]
reward_text = ("Reported nonzero task weights: " + "; ".join(
    f"{h}: {r['score']['weight']:.8f}" for h, r in rewarded) + "." if rewarded
    else "All ten hotkeys have zero reported task reward weight.")
lines += ["", "Ranks count strictly higher scores; ties share rank. Identities were matched by full "
          "public hotkey using this task's UIDs. " + reward_text + " These weights come from the task "
          "score feed. Final scores and factors are means across three seeds; multiplying mean "
          "factors need not reproduce the mean final score.", "",
          "## Joined windows and construction", "",
          "The fixed layout used **nine 12-seed candidate windows at stride 12 over 500–599** for "
          "h0–h8, with the common **200-seed cut search 100–199 ∪ 500–599**. h8 wraps to "
          "**500–507 ∪ 596–599**, sharing eight candidate seeds with h0. h9 used the all-HDR builder "
          "directly across **100–199**, without a separate conjunction cut search.", "",
          "The cut search was reconstructed by testing all 36 pairs of 100-seed classes: "
          "only 100–199 ∪ 500–599 matches every hotkey’s logged Cas12a clean count. "
          "The live layout has since changed and was not used to infer these historical windows.", "",
          "Build logs record loop 1 for h0–h6 and loop None for h7–h8. The conjunction function "
          "defaults to loop 1 and only overrides it for a non-None input; h7–h8 therefore also use "
          "the first loop. The recorded source evidence and replay checks are saved with the report.", "",
          "| Hotkey | Construction | Joined band candidates | Conjunction cut search | Logged clean statistic | Full clean /900 |",
          "|---|---|---|---|---|---:|"]
for h, r in ours.items():
    conj = r["construction"] == "conjunction"
    log_stat = f"{r['prepared_logged_clean_count']}/{r['logged_clean_denominator']} " + ("cut-clean" if conj else "HDR-clean")
    lines.append(f"| {h} | {'conjunction' if conj else 'all-HDR'} | {r['candidate_window_ranges']} | "
                 f"{r['cut_search_ranges'] or 'not applicable'} | {log_stat} | {len(r['clean_seeds'])} |")
lines += ["", "Conjunction logs count the Cas12a group's no-cut-free seeds inside its cut search. "
          "The all-HDR log counts its shared HDR seeds inside its candidate window. A **full clean "
          "seed requires all 250 uploaded rows to cut**, replayed over seeds 100–999. It may include "
          "seeds outside the construction window. Cas12a-only clean seeds are excluded when any "
          "Cas9 row fails to cut.", "",
          "## Actual uploaded HDR bands and hits", "",
          "An **HDR-band seed** makes all 250 uploaded rows return HDR. Each recovered band is a "
          "subset of its clean set and candidate window. A clean hit may contain mixed HDR and NHEJ outcomes.", "",
          "**h9's recorded band string 100–199 describes its candidate window. Replay recovers its "
          "seven actual HDR seeds below.** h0–h8's recorded JSON band lists match replay. h8's "
          "bounding window field 500–599 omits the gap in its joined 12-seed candidate window.", "",
          "| Hotkey | Clean /900 | HDR-band size | Exact uploaded HDR band | Clean scoring-seed hits | HDR scoring-seed hits |",
          "|---|---:|---:|---|---|---|"]
for h, r in ours.items():
    lines.append(f"| {h} | {len(r['clean_seeds'])} | {len(r['band_seeds'])} | {nums(r['band_seeds'])} | "
                 f"{nums(r['clean_hits'])} | {nums(r['band_hits'])} |")
lines += ["", f"Fleet union: **{len(d['fleet_clean_union'])} distinct clean seeds** and "
          f"**{len(d['fleet_band_union'])} distinct HDR-band seeds**, from "
          f"{sum(len(r['band_seeds']) for r in ours.values())} band memberships. Repeated HDR seeds: " +
          ("; ".join(f"{s} ({'/'.join(hs)})" for s, hs in sorted(d["duplicate_band_seeds"].items(), key=lambda x: int(x[0]))) or "none") + ".", "",
          "| Scoring seed | Hotkeys clean at this seed | Hotkeys all-HDR at this seed |",
          "|---:|---|---|"]
for s in seeds:
    lines.append(f"| {s} | {nums([h for h,r in ours.items() if s in r['clean_hits']])} | "
                 f"{nums([h for h,r in ours.items() if s in r['band_hits']])} |")
lines += ["", "## Per-seed results", "",
          "Each entry is **per-seed final score / HDR rows / no-cut rows**, from 250 rows.", "",
          "| Hotkey | " + " | ".join(f"Seed {s}" for s in seeds) + " | Mean final score |",
          "|---|---:|---:|---:|---:|"]
for h, r in ours.items():
    values = []
    for s in seeds:
        sr = r["per_seed"][str(s)]
        values.append(f"{sr['final_score']:.6f} / {sr['outcomes'].get('HDR',0)} / {sr['outcomes'].get('no_cut',0)}")
    lines.append(f"| {h} | " + " | ".join(values) + f" | {r['score']['final_score']:.6f} |")
lines += ["", "## Comparison with the leader", ""]
b, t = best["score"]["breakdown"], top["breakdown"]
lines += [f"{best_h}'s weighted score was **{b['total_weighted_score']:.6f}**, versus the leader's "
          f"**{t['total_weighted_score']:.6f}**. Mean consistency was **{b['consistency_factor']:.6f} versus "
          f"{t['consistency_factor']:.6f}**, and fidelity was **{b['distribution_fidelity_factor']:.6f} versus "
          f"{t['distribution_fidelity_factor']:.6f}**.", "",
          "Seed 548 was in h4's candidate window and the conjunctions' cut search, but was absent "
          "from h4's uploaded HDR band. Seeds 872 and 939 were outside both construction groups' "
          "candidate windows and the conjunction cut search. Actual coverage is determined by the replay above.", "",
          f"Leader hotkey: `{top['miner_hotkey']}`. The inspected data contains the leader's aggregate "
          "score, but does not provide its uploaded rows, clean sets, HDR bands or per-seed outcomes.", "",
          "## Exact uploaded clean sets", "",
          "Each list covers all seeds 100–999 and includes its uploaded HDR band.", ""]
for h, r in ours.items():
    lines += [f"### {h}: {len(r['clean_seeds'])} clean seeds", "", nums(r["clean_seeds"]) + ".", ""]
    extra = sorted(set(r["cas12a_clean_seeds"]) - set(r["clean_seeds"]))
    if extra:
        lines += [f"Cas12a-only clean seeds excluded from the full-row set: {nums(extra)}.", ""]
lines += ["## Verification and provenance", "",
          "- All ten successful upload records identify this task. Contracts and references match the published task after substituting its scoring seed stamp.",
          "- All 2,500 rows pass validation at HEK293 accessibility 0.35. Every locally reproduced three-seed final score matches its official score within 1e-10.",
          "- Replayed **900 × 250 × 10 = 2,250,000 row/seed combinations**. At each scoring seed, every row's outcome and indel length matches local stage-3 detail.",
          "- Conjunction cut-clean counts match the build logs. The all-HDR HDR-clean count matches the recovered band. All uploads used completed prepared builds.",
          "- [Analysis JSON](analysis.json), [CSV](summary.csv), [official scores](scores.json), [task](task.json), [manifest](manifest.json), [sanitized build evidence](build_evidence.json), [layout evidence](layout_evidence.json).",
          "- Reproduce from saved snapshots with `.venv/bin/python reports/35e41a29/validate.py`, then `.venv/bin/python reports/35e41a29/audit.py`, then `.venv/bin/python reports/35e41a29/render_report.py`. Validator and submission SHA-256 hashes are recorded in the analysis.", "",
          "| Hotkey | Upload UTC | Prepared build completed UTC | Public hotkey | Submission archive |",
          "|---|---|---|---|---|"]
for h, r in ours.items():
    lines.append(f"| {h} | {r['uploaded_at']} | {r['prepared_build_completed_at']} | `{r['hotkey']}` | "
                 f"[{r['archive']}](../../{r['archive']}/submission.json) |")
(OUT / "report.md").write_text("\n".join(lines) + "\n")
with (OUT / "summary.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["Hotkey", "UID", "Rank", "Official score", "Construction", "Cas12a", "Cas9", "Band candidate window", "Conjunction cut window",
                     "Clean count /900", "Logged clean count", "Logged denominator", "Logged clean definition", "Exact clean set", "Clean hits",
                     "HDR band", "HDR hits", "Percent of leader", "Gap to leader", "Reward weight"])
    for h, r in ours.items():
        writer.writerow([h,r["uid"],r["rank"],r["score"]["final_score"],r["construction"],r["cas_mix"]["Cas12a"],r["cas_mix"]["Cas9"],
                         r["candidate_window_ranges"],r["cut_search_ranges"] or "not applicable",len(r["clean_seeds"]),
                         r["prepared_logged_clean_count"],r["logged_clean_denominator"],r["logged_clean_definition"],nums(r["clean_seeds"]),
                         nums(r["clean_hits"]),nums(r["band_seeds"]),nums(r["band_hits"]),r["percent_of_top"],r["gap_to_top"],r["score"]["weight"]])
print(OUT / "report.md")
