"""Generate and verify this task's report without changing miner data or production code."""
import ast
from collections import Counter
import csv
from datetime import datetime, timedelta, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
TASK_ID = "35e41a29-5658-494b-9baa-7c83a18c25bf"
SEEDS = [872, 548, 939]


def read(path):
    return json.loads(Path(path).read_text())


def write(name, data):
    (OUT / name).write_text(json.dumps(data, indent=2) + "\n")


def prepare():
    listing = read("/tmp/niome-35e41a29-tasks.json")["items"]
    task = next(t for t in listing if t["id"] == TASK_ID)
    assert task["content"]["contract"]["seed"] == "872,548,939"
    write("task.json", task)
    score_path = Path("/tmp/niome-35e41a29-scores.json")
    shutil.copyfile(score_path, OUT / "scores.json")
    shutil.copyfile("/tmp/niome-35e41a29-cell-types.json", OUT / "cell_types.json")
    previous = read(ROOT / "reports/e6caf784/manifest.json")["hotkeys"]
    start = task["created_at"].replace("T", " ")[:19]
    later = [t["created_at"] for t in listing if t["created_at"] > task["created_at"]]
    end = (min(later).replace("T", " ")[:19] if later else
           (datetime.fromisoformat(task["created_at"]) + timedelta(hours=3)).strftime("%Y-%m-%d %H:%M:%S"))
    manifest = {"task_id": TASK_ID,
                "fetched_at": datetime.fromtimestamp(score_path.stat().st_mtime, timezone.utc).isoformat(),
                "scoring_seeds_available": True, "hotkeys": {}}
    evidence = {}
    for i in range(10):
        h = f"h{i}"
        base = ROOT / "data/inst" / ("niome_hotkey" + (str(i) if i else ""))
        archives = [p.parent for p in (base / "result").glob("*/last_upload.json")
                    if read(p).get("task_id") == TASK_ID]
        assert len(archives) == 1, (h, archives)
        folder = archives[0]
        assert read(folder / "last_upload.json")["submitted"]
        validation = folder / "validation.json"
        existing = (validation.exists() and read(validation).get("task_id") == TASK_ID
                    and read(validation).get("scored_seeds") == SEEDS)
        manifest["hotkeys"][h] = {
            "hotkey": previous[h]["hotkey"], "archive": str(folder.relative_to(ROOT)),
            "validation_dir": str((folder if existing else OUT / h).relative_to(ROOT)),
            "existing_validation": existing, "submitted": True,
            "construction": "conjunction" if i < 9 else "all_hdr",
        }
        log = Path(f"/root/.pm2/logs/miner-h{i}-error.log")
        log_lines = log.read_text(errors="replace").splitlines()
        ready = [line for line in log_lines if f"Prefetch: task {TASK_ID} ready" in line]
        assert len(ready) == 1, (h, ready)
        events = [re.sub(r"https?://\S+", "[URL redacted]", line)
                  for line in log_lines
                  if TASK_ID in line or (start <= line[:19] and line[:23] <= ready[0][:23]
                                        and "Build:" in line)]
        evidence[h] = {"log_path": str(log), "events": events,
                       "current_window_record": read(base / "window_used.json").get(TASK_ID)}
        mix = dict(Counter(r["cas_system"] for r in read(folder / "submission.json")))
        print(h, "archive", folder.name, "cas_mix", mix, "existing_validation", existing, flush=True)
        for line in events:
            if "Build: conjunction " in line or "Build: all-HDR " in line:
                print(line, flush=True)
    write("manifest.json", manifest)
    write("build_evidence.json", evidence)
    layout = ROOT / "joined_window.py"
    source = layout.read_text()
    constants = {}
    for node in ast.parse(source).body:
        targets = node.targets if isinstance(node, ast.Assign) else [node.target] if isinstance(node, ast.AnnAssign) else []
        for target in targets:
            if isinstance(target, ast.Name) and target.id in ("BAND_CLASSES", "CUT_CLASSES", "HDR_ONLY_CLASSES", "BAND_SUB_WIDTH", "BAND_STRIDE", "CELL_LAYOUT"):
                constants[target.id] = ast.literal_eval(node.value)
    print("Current layout constants", constants, flush=True)
    assert constants["CELL_LAYOUT"]["HEK293"] == {"conj_hk": 9, "width": 12, "stride": 12, "hdr_width": 100, "hdr_stride": 100}, constants
    cj = ROOT / "niome_subnet/genomics/conjunction.py"
    cj_source = cj.read_text()
    assert re.search(r"band_loop:\s*int\s*=\s*1\b", cj_source)
    assert "if band_loop is not None:" in cj_source
    write("layout_evidence.json", {
        "note": "Current constants differ from this historical task and are retained only as audit-time evidence. Historical candidate windows come from task-specific logs. The 200-seed cut is reconstructed by uniquely matching all nine logged Cas12a clean counts across all possible pairs of 100-seed classes in the replay.",
        "constants": constants, "layout_sha256": hashlib.sha256(layout.read_bytes()).hexdigest(),
        "conjunction_sha256": hashlib.sha256(cj.read_bytes()).hexdigest(),
        "layout_band_loop": [ast.get_source_segment(source, node) for node in ast.parse(source).body
                             if isinstance(node, ast.FunctionDef) and node.name == "band_loop"],
        "conjunction_loop_lines": [f"{i+1}: {line}" for i, line in enumerate(cj_source.splitlines())
                                   if "band_loop" in line],
    })
    shutil.copyfile(ROOT / "reports/f18ee409/validate.py", OUT / "validate.py")
    audit = (ROOT / "reports/f18ee409/audit.py").read_text()
    replacements = {
        "f18ee409-bdb0-4e80-9818-e25f9f456f46": TASK_ID,
        "seven conjunctions and three all-HDR uploads for task f18ee409": "nine conjunctions and one all-HDR upload for task 35e41a29",
        "[214, 249, 784]": "[872, 548, 939]",
        "conjunction_space = list(range(200, 300))": "conjunction_space = list(range(500, 600))",
        "hdr_space = list(range(400, 500))": "hdr_space = list(range(100, 200))",
        "conjunction_cut = conjunction_space + hdr_space": "conjunction_cut = sorted(conjunction_space + hdr_space)",
        "is_conjunction = i < 7": 'is_conjunction = i < 9\n        loop_label = "1" if i < 7 else "None"',
        "(15 * i, 15) if is_conjunction else (34 * (i - 7), 34)": "(12 * i, 12) if is_conjunction else (0, 100)",
        'assert "15 seeds (shared window, loop 1 of 7) | cut 200 seeds (narrow) | k=8"': 'assert f"12 seeds (shared window, loop {loop_label} of 9) | cut 200 seeds (narrow) | k=8"',
        "logged == (80, len(set(ca_hdr) & set(candidates)), 34)": "logged == (80, len(set(ca_hdr) & set(candidates)), 100)",
        'expected_record = label if i < 9 else "joined 34 seeds " + label': "expected_record = label",
        "assert len(band) == 8 and set(band) <= set(candidates)": "assert len(band) == (8 if is_conjunction else 7) and set(band) <= set(candidates)",
        '"conjunction_15_seed_window" if is_conjunction else "all_hdr_34_seed_window"': '"conjunction_12_seed_window" if is_conjunction else "all_hdr_100_seed_window"',
        '"candidate_offset": offset, "band_loop": 1 if is_conjunction else None,': '"candidate_offset": offset, "band_loop": 1 if is_conjunction else None, "logged_band_loop": 1 if i < 7 else None,',
        "h0–h6: width-15 windows at stride 15 over 200–299, loop 1, cut 200–299 ∪ 400–499. h7–h9: all-HDR search over width-34 windows at stride 34 over 400–499. h6 and h9 wrap.": "h0–h8: width-12 windows at stride 12 over 500–599; h8 wraps to 500–507 ∪ 596–599. Shared cut 100–199 ∪ 500–599. h9: all-HDR search across 100–199.",
        "For h7–h9, log clean 8/34": "For h9, log clean 7/100",
        "in the 34-seed candidate window": "in the 100-seed candidate window",
        "those eight seeds": "those seven seeds",
    }
    for before, after in replacements.items():
        assert before in audit, before
        audit = audit.replace(before, after)
    insertion = '    matching_cut_classes = []\n    for a in range(9):\n        for b in range(a + 1, 9):\n            candidate_cut = set(range(100 + 100*a, 200 + 100*a)) | set(range(100 + 100*b, 200 + 100*b))\n            if all(len(candidate_cut & set(r["cas12a_clean_seeds"])) == r["prepared_logged_clean_count"]\n                   for r in result["hotkeys"].values() if r["construction"] == "conjunction"):\n                matching_cut_classes.append([a, b])\n    assert matching_cut_classes == [[0, 4]], matching_cut_classes\n    result["cut_space_inference"] = {"matching_class_pairs": matching_cut_classes,\n        "tested_class_pairs": 36, "matched_hotkeys": 9,\n        "basis": "Unique pair of 100-seed classes matching all nine logged Cas12a clean counts; current layout differs."}\n'
    marker = '    records = result["hotkeys"].values()'
    assert marker in audit
    audit = audit.replace(marker, insertion + marker)
    ast.parse(audit)
    (OUT / "audit.py").write_text(audit)


def check():
    d = read(OUT / "analysis.json")
    ours = d["hotkeys"]
    assert len(ours) == 10 and d["seeds"] == SEEDS
    assert d["best_our_hotkey"] == "h3" and ours["h3"]["rank"] == 37
    assert ours["h8"]["candidate_window"] == list(range(500, 508)) + list(range(596, 600))
    assert ours["h9"]["candidate_window"] == list(range(100, 200))
    csv_rows = list(csv.DictReader((OUT / "summary.csv").open()))
    assert len(csv_rows) == 10 and {r["Hotkey"] for r in csv_rows} == set(ours)
    for row in csv_rows:
        h = row["Hotkey"]
        r = ours[h]
        assert float(row["Official score"]) == r["score"]["final_score"]
        assert int(row["Clean count /900"]) == len(r["clean_seeds"])
        assert r["local_score_matches_api"] and r["stage3_replay_matches_local_detail"]
        assert r["seed_independent_factors_match_api"] and r["rows"] == 250
        assert r["cas_mix"] == {"Cas12a": 80, "Cas9": 170}
        assert set(r["band_seeds"]) <= set(r["clean_seeds"])
        assert len(r["band_seeds"]) == (8 if h != "h9" else 7)
        if h != "h9":
            assert r["latest_record_band_matches_upload"]
        else:
            assert r["latest_record_band_semantics"] == "candidate window; not the actual HDR seeds"
    report = (OUT / "report.md").read_text()
    for target in re.findall(r"\]\(([^)]+)\)", report):
        if not target.startswith("https://"):
            assert (OUT / target).exists(), target
    for filename in ("manifest.json", "build_evidence.json", "layout_evidence.json", "audit.py", "render_report.py", "report.md"):
        text = (OUT / filename).read_text()
        assert not re.search(r"X-Amz-|AWSAccessKeyId|Signature=", text), filename
        assert "f18ee409" not in text, filename
    print("All score, replay, CSV, identity and report-link checks passed.", flush=True)
    print("\n".join(report.splitlines()[:115]), flush=True)
    subprocess.run(["git", "status", "--short"], cwd=ROOT, check=True)


if __name__ == "__main__":
    if sys.argv[1:] == ["--check-only"]:
        check()
    else:
        assert sys.argv[1:] in ([], ["--render-only"])
        if not sys.argv[1:]:
            prepare()
            for name in ("validate.py", "audit.py"):
                subprocess.run([str(ROOT / ".venv/bin/python"), str(OUT / name)], cwd=ROOT, check=True)
        subprocess.run([str(ROOT / ".venv/bin/python"), str(OUT / "render_report.py")], cwd=ROOT, check=True)
        check()
