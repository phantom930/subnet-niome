"""Replay the actual uploads for a31f1224; distinguish later builds and h9's absent archive."""
import hashlib
import importlib.util
import json
import re
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
TASK_ID = "a31f1224-969b-4195-bb49-8220aa8800e5"


def read(path):
    return json.loads(Path(path).read_text())


def ranges(values):
    groups = []
    for value in sorted(set(values)):
        if groups and value == groups[-1][-1] + 1:
            groups[-1].append(value)
        else:
            groups.append([value])
    return ", ".join(str(g[0]) if len(g) == 1 else f"{g[0]}–{g[-1]}" for g in groups)


def main():
    spec = importlib.util.spec_from_file_location("audit_stage3", ROOT / "niome_subnet/genomics/validation/stage3.py")
    stage3 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(stage3)
    task, manifest, evidence = (read(OUT / f"{f}.json") for f in ("task", "manifest", "build_evidence"))
    seeds = [int(x) for x in task["content"]["contract"]["seed"].split(",")]
    scores = [r for r in read(OUT / "scores.json")["items"] if r["task_id"] == TASK_ID]
    best = {}
    for r in scores:
        if r["miner_hotkey"] not in best or r["final_score"] > best[r["miner_hotkey"]]["final_score"]:
            best[r["miner_hotkey"]] = r
    ranked = sorted(best.values(), key=lambda r: -r["final_score"])
    ranks = {hk: 1 + sum(x["final_score"] > r["final_score"] for x in ranked) for hk, r in best.items()}
    result = {
        "task_id": TASK_ID, "task_created_at": task["created_at"], "cell_type": task["content"]["contract"]["cell_type"],
        "seeds": seeds, "generated_at": datetime.now(timezone.utc).isoformat(), "fetched_at": manifest["fetched_at"],
        "source_api": f"https://niome-api.genomes.io/api/v3/miners/scores?task_id={TASK_ID}&limit=40000",
        "score_records": len(scores), "unique_miners": len(ranked),
        "validators": sorted({r["validator_hotkey"] for r in scores}), "top_miner": ranked[0],
        "top10_cutoff": ranked[9]["final_score"], "rank_definition": "1 plus number of strictly higher scores; exact ties share rank.",
        "code_sha256": {f: hashlib.sha256((ROOT / f).read_bytes()).hexdigest() for f in
                        ("calc.py", "niome_subnet/genomics/validation/stage12.py", "niome_subnet/genomics/validation/stage3.py",
                         "niome_subnet/genomics/validation/stage4.py", "niome_subnet/genomics/validation/stage5.py")},
        "definitions": {"clean_seeds": "All 250 uploaded rows cut at this seed, replayed over 100-999.",
                        "band_seeds": "All 250 uploaded rows return HDR, replayed over 100-999.",
                        "candidate_window": "Historical 150-seed slice: h0-h5 offsets 0/50/100/150/200/250; h6 override 25; h8 reconstructed at 240 under the new stride-30 layout. h7 all-HDR uses the full joined 300.",
                        "missing_archive": "h9 has no task upload archive and zero valid experiments in the API. Actual submitted clean set and band cannot be reconstructed."},
        "hotkeys": {},
    }
    for i in range(10):
        h = f"h{i}"
        entry, ev = manifest["hotkeys"][h], evidence[h]
        hk, construction = entry["hotkey"], entry["construction"]
        score = best[hk]
        joined = ([*range(100, 200), *range(300, 400), *range(800, 900)] if i >= 8
                  else [*range(200, 300), *range(400, 500), *range(800, 900)])
        offset = i * 50 if i < 6 else 25 if i == 6 else i * 30
        candidates = joined if i == 7 else [joined[(offset + j) % 300] for j in range(150)]
        original_events = [e for e in ev["events"] if "2026-09-18 15:04" <= e[:19] < "2026-09-18 15:23"]
        window_lines = [e for e in original_events if "Build: window " in e]
        assert window_lines
        expected_window = "100-199,300-399,800-899" if i >= 8 else "200-299,400-499,800-899"
        assert all(expected_window in e for e in window_lines)
        rec = {
            "hotkey": hk, "uid": score["miner_uid"], "construction": construction, "archive": entry["archive"],
            "validation_dir": entry["validation_dir"], "score": score, "rank": ranks[hk],
            "gap_to_top": ranked[0]["final_score"] - score["final_score"],
            "percent_of_top": 100 * score["final_score"] / ranked[0]["final_score"],
            "joined_band_space": joined, "joined_band_ranges": ranges(joined),
            "joined_source": "env_pin_no_entry_joined" if i >= 8 else "plan_joined (14:17 plan)",
            "candidate_offset": None if i == 7 else offset, "candidate_window": candidates,
            "candidate_window_basis": "New ten-hotkey stride-30 layout, consistent with recovered band; reconstructed, not logged directly." if i >= 8 else "Original stride-50 layout with h6 offset 25; h7 all-HDR used all 300 seeds.",
            "candidate_window_ranges": ranges(candidates), "candidate_window_applies_to_upload": i not in (5, 9),
            "cut_search_space": list(range(100, 1000)) if construction == "conjunction" else None,
            "latest_record_band": sorted(json.loads(ev["current_window_record"]["band"])),
            "latest_record_at": ev["current_window_record"]["at"],
        }
        if not entry["archive"]:
            assert score["final_score"] == 0 and score["breakdown"]["n_valid_experiments"] == 0
            assert not any("Submitted 250 rows" in e or "Received genomics task" in e for e in ev["events"])
            prepared = [e for e in original_events if "Build: conjunction " in e]
            assert len(prepared) == 1
            rec.update({"rows": None, "clean_seeds": None, "band_seeds": None, "clean_hits": None,
                        "band_hits": None, "submitted": None, "logged_prepared_clean_count": int(re.search(r"clean (\d+)/900", prepared[0])[1]),
                        "status": "No received-task/upload event or task archive found; API reports zero score and zero valid rows. Prepared band is unverified against submitted rows."})
            result["hotkeys"][h] = rec
            print(h, "no uploaded archive; official score 0, valid rows 0", flush=True)
            continue
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
        assert {k:v for k,v in contract.items() if k != "seed"} == {k:v for k,v in task["content"]["contract"].items() if k != "seed"}
        reference = read(folder / "hbb_reference.json")
        reference["challenge"]["seed"] = task["content"]["hbb_reference"]["challenge"]["seed"]
        assert reference == task["content"]["hbb_reference"]
        assert abs(validation["score"]["final_score"] - score["final_score"]) < 1e-10
        clean, band, ca_clean, seed_results = [], [], [], {}
        started = time.monotonic()
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
                seed_results[str(seed)] = {"outcomes": dict(counts),
                    "outcomes_by_cas": {cas: dict(Counter(r["outcome"] for r in simulated if r["cas"] == cas)) for cas in ("Cas9", "Cas12a")},
                    "in_joined_space": seed in joined, "clean": seed in clean, "band": seed in band,
                    "final_score": validation["per_seed_final_score"][str(seed)]}
        assert set(band) <= set(clean) <= set(ca_clean)
        cas_mix = dict(Counter(r["cas_system"] for r in rows))
        if construction == "conjunction":
            builds = [e for e in original_events if "Build: conjunction " in e]
            assert len(builds) == 1
            logged = tuple(map(int, re.search(r"band (\d+) k=(\d+) group (\d+) clean (\d+)/(\d+)", builds[0]).groups()))
            assert logged == (11, 11, 100, len(ca_clean), 900)
            assert clean == ca_clean
            assert len(band) == 11 and set(band) <= set(candidates), (h, band, ranges(candidates))
            assert cas_mix == {"Cas12a": 100, "Cas9": 150}
        elif construction == "all_hdr":
            assert any("Build: all-HDR " in e and "clean 11/300" in e for e in original_events)
            assert len(band) == 11 and set(band) <= set(joined)
        else:
            assert cas_mix == {"Cas9": 174, "Cas12a": 76}
            assert any("prepared round unusable" in e for e in original_events)
        rec.update({"rows": len(rows), "cas_mix": cas_mix, "submitted": True,
            "submission_sha256": hashlib.sha256((folder / "submission.json").read_bytes()).hexdigest(),
            "clean_seeds": clean, "band_seeds": band, "cas12a_clean_seeds": ca_clean,
            "clean_seeds_in_joined_space": sorted(set(clean) & set(joined)),
            "clean_seeds_outside_joined_space": sorted(set(clean) - set(joined)),
            "clean_hits": sorted(set(clean) & set(seeds)), "band_hits": sorted(set(band) & set(seeds)),
            "latest_record_band_matches_upload": rec["latest_record_band"] == band,
            "per_seed": seed_results, "local_score_matches_api": True, "stage3_replay_matches_local_detail": True,
            "scan_seconds": round(time.monotonic() - started, 3)})
        result["hotkeys"][h] = rec
        (OUT / "analysis.json").write_text(json.dumps(result, indent=2) + "\n")
        print(f"{h}: rank={rec['rank']} final={score['final_score']:.8f} clean={len(clean)} band={band} clean_hits={rec['clean_hits']} band_hits={rec['band_hits']} latest_band_match={rec['latest_record_band_matches_upload']}", flush=True)
    uploaded = [r for r in result["hotkeys"].values() if r["submitted"]]
    result["uploaded_count"] = len(uploaded)
    result["fleet_clean_union"] = sorted(set().union(*(set(r["clean_seeds"]) for r in uploaded)))
    result["fleet_band_union"] = sorted(set().union(*(set(r["band_seeds"]) for r in uploaded)))
    result["fleet_clean_hits"] = sorted(set(seeds) & set(result["fleet_clean_union"]))
    result["fleet_band_hits"] = sorted(set(seeds) & set(result["fleet_band_union"]))
    (OUT / "analysis.json").write_text(json.dumps(result, indent=2) + "\n")
    print("Fleet unions:", len(result["fleet_clean_union"]), "clean;", len(result["fleet_band_union"]), "HDR band", flush=True)


if __name__ == "__main__":
    main()
