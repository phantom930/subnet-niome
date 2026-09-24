"""Replay seven wide-cut conjunctions for 00b2e47c; h1–h3 have no build/upload."""
import hashlib
import importlib.util
import json
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
TASK_ID = "00b2e47c-f674-4f13-939d-cb4ad624eacd"
sys.path.insert(0, str(ROOT))


def read(path):
    return json.loads(Path(path).read_text())


def ranges(values):
    groups = []
    for value in sorted(set(values)):
        if groups and value == groups[-1][-1] + 1:
            groups[-1].append(value)
        else:
            groups.append([value])
    return " ∪ ".join(str(g[0]) if len(g) == 1 else f"{g[0]}–{g[-1]}" for g in groups)


def main():
    spec = importlib.util.spec_from_file_location("audit_stage3", ROOT / "niome_subnet/genomics/validation/stage3.py")
    stage3 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(stage3)
    task, manifest, evidence = (read(OUT / f"{f}.json") for f in ("task", "manifest", "build_evidence"))
    seeds = [int(x) for x in task["content"]["contract"]["seed"].split(",")]
    snapshot = read(OUT / "scores.json")
    scores = [r for r in snapshot["items"] if r["task_id"] == TASK_ID]
    assert len(scores) == snapshot["pagination"]["total"]
    best = {}
    for row in scores:
        hk = row["miner_hotkey"]
        if hk not in best or row["final_score"] > best[hk]["final_score"]:
            best[hk] = row
    ranked = sorted(best.values(), key=lambda r: -r["final_score"])
    ranks = {hk: 1 + sum(x["final_score"] > r["final_score"] for x in ranked) for hk, r in best.items()}
    predicted = list(range(100, 400))
    complement = sorted(set(range(100, 1000)) - set(predicted))
    full = list(range(100, 1000))
    active = ["h0", "h4", "h5", "h6", "h7", "h8", "h9"]
    result = {
        "task_id": TASK_ID, "task_created_at": task["created_at"],
        "cell_type": task["content"]["contract"]["cell_type"], "seeds": seeds,
        "generated_at": datetime.now(timezone.utc).isoformat(), "fetched_at": manifest["fetched_at"],
        "source_api": f"https://niome-api.genomes.io/api/v3/miners/scores?task_id={TASK_ID}&limit=40000",
        "score_records": len(scores), "unique_miners": len(ranked),
        "validators": sorted({r["validator_hotkey"] for r in scores}), "top_miner": ranked[0],
        "top10_cutoff": ranked[9]["final_score"], "top10": ranked[:10],
        "rank_definition": "1 plus the number of strictly higher scores; ties share rank.",
        "predicted_windows": [100, 200, 300], "predicted_joined": predicted,
        "actual_band_space": full, "active_layout": active,
        "complement_joined": complement, "plan_generated_at": "2026-09-20T20:17:15.356931+00:00",
        "code_sha256": {f: hashlib.sha256((ROOT / f).read_bytes()).hexdigest() for f in
                        ("calc.py", "niome_subnet/genomics/validation/stage12.py", "niome_subnet/genomics/validation/stage3.py",
                         "niome_subnet/genomics/validation/stage4.py", "niome_subnet/genomics/validation/stage5.py")},
        "definitions": {
            "clean_seeds": "All 250 uploaded rows cut at this seed, replayed over 100–999.",
            "band_seeds": "All 250 uploaded rows return HDR at this seed, replayed over 100–999.",
            "candidate_window": "Historical 300-seed slice of full 100–999 at offsets 0/100/200/300/400/500/600 for h0/h4/h5/h6/h7/h8/h9. The cut-search window is the full 100–999 for all seven.",
            "missing_upload": "h1–h3 have no task build/window/upload records or official score entries. Their submitted sets and scores are unavailable, not zero.",
        },
        "hotkeys": {},
    }
    for i in range(10):
        h = f"h{i}"
        entry, ev = manifest["hotkeys"][h], evidence[h]
        score = best.get(entry["hotkey"])
        if not entry["archive"]:
            assert score is None and ev["current_window_record"] is None and not ev["events"]
            result["hotkeys"][h] = {
                "hotkey": entry["hotkey"], "uid": None, "construction": "no_build_or_upload",
                "archive": None, "validation_dir": None, "submitted": False, "rows": None,
                "score": None, "rank": None, "gap_to_top": None, "percent_of_top": None,
                "band_space": None, "candidate_window": None, "candidate_window_ranges": None,
                "cut_search_space": None, "clean_seeds": None, "band_seeds": None,
                "clean_hits": None, "band_hits": None,
                "status": "No task build/window/upload records or official score entry.",
            }
            print(h, "no task build, upload or score entry", flush=True)
            continue
        assert score is not None
        joined = full
        offset = active.index(h) * 100
        group = "full_900"
        candidates = joined[offset:offset + 300]
        events = ev["events"]
        window_lines = [e for e in events if "Build: window " in e]
        assert window_lines and all("100-399" in e and "2026-09-20T20:17" in e for e in window_lines)
        candidate_lines = [e for e in events if "Build: band candidates " in e]
        assert candidate_lines and all(f"(300 seeds, predicted space of 900, offset {offset}) | cut 900 seeds (wide) | k=11" in e for e in candidate_lines)
        assert any("using the prepared submission (250 rows" in e for e in events)
        assert not any("using the ordinary construction" in e for e in events)
        uploaded_events = [e for e in events if f"Submitted 250 rows for task {TASK_ID}" in e]
        assert len(uploaded_events) == 1
        folder, validation_dir = ROOT / entry["archive"], ROOT / entry["validation_dir"]
        rows = read(folder / "submission.json")
        valid = read(validation_dir / "validation/valid_experiments.json")
        validation = read(validation_dir / "validation.json")
        upload = read(folder / "last_upload.json")
        assert upload["submitted"] and upload["task_id"] == TASK_ID and len(rows) == 250
        assert rows == [e["experiment"] for e in valid]
        assert {e["features"]["cell_type_accessibility"] for e in valid} == {0.87}
        assert validation["scored_seeds"] == seeds and validation["contract_pairing"]["verified"]
        contract = read(folder / "contract.json")
        assert {k: v for k, v in contract.items() if k != "seed"} == {k: v for k, v in task["content"]["contract"].items() if k != "seed"}
        reference = read(folder / "hbb_reference.json")
        reference["challenge"]["seed"] = task["content"]["hbb_reference"]["challenge"]["seed"]
        assert reference == task["content"]["hbb_reference"]
        assert abs(validation["score"]["final_score"] - score["final_score"]) < 1e-10
        assert abs(sum(validation["per_seed_final_score"].values()) / len(seeds) - score["final_score"]) < 1e-10
        clean, band, ca_clean, seed_results = [], [], [], {}
        for seed in range(100, 1000):
            simulated = [stage3.simulate(e, seed) for e in valid]
            counts = Counter(r["outcome"] for r in simulated)
            if counts["no_cut"] == 0:
                clean.append(seed)
            if counts["HDR"] == len(rows):
                band.append(seed)
            if all(r["outcome"] != "no_cut" for r in simulated if r["cas"] == "Cas12a"):
                ca_clean.append(seed)
            if seed in seeds:
                detail = read(validation_dir / f"stage3_detail.seed{seed}.json")
                old = {r["experiment_id"]: r for r in detail["experiments"]}
                assert all(r["outcome"] == old[r["experiment_id"]]["outcome"] and r["indel_length"] == old[r["experiment_id"]]["indel_length"] for r in simulated)
                seed_results[str(seed)] = {
                    "outcomes": dict(counts),
                    "outcomes_by_cas": {cas: dict(Counter(r["outcome"] for r in simulated if r["cas"] == cas)) for cas in ("Cas9", "Cas12a")},
                    "in_predicted_space": seed in predicted, "in_band_space": seed in joined,
                    "clean": seed in clean, "band": seed in band,
                    "final_score": validation["per_seed_final_score"][str(seed)],
                }
        assert set(band) <= set(clean) <= set(ca_clean)
        cas_mix = dict(Counter(r["cas_system"] for r in rows))
        builds = [e for e in events if "Build: conjunction " in e]
        assert len(builds) == 1
        logged = tuple(map(int, re.search(r"band (\d+) k=(\d+) group (\d+) clean (\d+)/(\d+)", builds[0]).groups()))
        assert logged == (11, 11, 100, len(ca_clean), 900)
        assert clean == ca_clean
        assert len(band) == 11 and set(band) <= set(candidates)
        assert cas_mix == {"Cas12a": 100, "Cas9": 150}
        assert uploaded_events[0][:23] > builds[0][:23]
        record_band = sorted(json.loads(ev["current_window_record"]["band"]))
        assert len(record_band) == 11 and record_band == band
        rec = {
            "hotkey": entry["hotkey"], "uid": score["miner_uid"], "construction": entry["construction"],
            "archive": entry["archive"], "validation_dir": entry["validation_dir"],
            "submitted": True, "uploaded_at": uploaded_events[0][:23],
            "rows": len(rows), "cas_mix": cas_mix, "score": score, "rank": ranks[entry["hotkey"]],
            "gap_to_top": ranked[0]["final_score"] - score["final_score"],
            "percent_of_top": 100 * score["final_score"] / ranked[0]["final_score"],
            "band_group": group, "band_space": joined, "band_space_ranges": ranges(joined),
            "candidate_offset": offset, "candidate_window": candidates, "candidate_window_ranges": ranges(candidates),
            "candidate_window_applies_to_upload": True,
            "cut_search_space": full,
            "clean_seeds": clean, "cas12a_clean_seeds": ca_clean, "band_seeds": band,
            "clean_seeds_in_predicted_space": sorted(set(clean) & set(predicted)),
            "clean_seeds_in_band_space": sorted(set(clean) & set(joined)),
            "clean_seeds_in_candidate_window": sorted(set(clean) & set(candidates)),
            "clean_seeds_outside_candidate_window": sorted(set(clean) - set(candidates)),
            "clean_hits": sorted(set(clean) & set(seeds)), "band_hits": sorted(set(band) & set(seeds)),
            "latest_record_at": ev["current_window_record"]["at"], "latest_record_band": record_band,
            "latest_record_band_matches_upload": record_band == band,
            "logged_prepared_clean_count": logged[3], "prepared_completed_at": builds[0][:23],
            "submission_sha256": hashlib.sha256((folder / "submission.json").read_bytes()).hexdigest(),
            "per_seed": seed_results, "local_score_matches_api": True, "stage3_replay_matches_local_detail": True,
        }
        result["hotkeys"][h] = rec
        (OUT / "analysis.json").write_text(json.dumps(result, indent=2) + "\n")
        print(f"{h}: rank={rec['rank']} final={score['final_score']:.8f} clean={len(clean)} band={band} clean_hits={rec['clean_hits']} HDR_hits={rec['band_hits']}", flush=True)
    records = [r for r in result["hotkeys"].values() if r["submitted"]]
    result["uploaded_count"] = len(records)
    result["fleet_clean_union"] = sorted(set().union(*(set(r["clean_seeds"]) for r in records)))
    result["fleet_band_union"] = sorted(set().union(*(set(r["band_seeds"]) for r in records)))
    result["fleet_clean_hits"] = sorted(set(seeds) & set(result["fleet_clean_union"]))
    result["fleet_band_hits"] = sorted(set(seeds) & set(result["fleet_band_union"]))
    (OUT / "analysis.json").write_text(json.dumps(result, indent=2) + "\n")
    print("Fleet unions:", len(result["fleet_clean_union"]), "clean;", len(result["fleet_band_union"]), "HDR band", flush=True)


if __name__ == "__main__":
    main()
