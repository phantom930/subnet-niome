"""Read this task's upload archives and write its reproducible audit under reports/e6caf784."""
import ast
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
TASK_ID = "e6caf784-f2af-48cc-938a-1377c441a1dd"


def read(path):
    return json.loads(Path(path).read_text())


def write(name, data):
    (OUT / name).write_text(json.dumps(data, indent=2) + "\n")


def prepare():
    task = next(t for t in read("/tmp/niome-e6caf784-tasks.json")["items"] if t["id"] == TASK_ID)
    assert task["content"]["contract"]["seed"] == "481,208,809"
    write("task.json", task)
    shutil.copyfile("/tmp/niome-e6caf784-scores.json", OUT / "scores.json")
    shutil.copyfile("/tmp/niome-e6caf784-cell-types.json", OUT / "cell_types.json")
    previous = read(ROOT / "reports/f18ee409/manifest.json")["hotkeys"]
    manifest = {"task_id": TASK_ID, "fetched_at": datetime.now(timezone.utc).isoformat(),
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
                    and read(validation).get("scored_seeds") == [481, 208, 809])
        manifest["hotkeys"][h] = {
            "hotkey": previous[h]["hotkey"], "archive": str(folder.relative_to(ROOT)),
            "validation_dir": str((folder if existing else OUT / h).relative_to(ROOT)),
            "existing_validation": existing, "submitted": True,
            "construction": "conjunction" if i < 7 else "all_hdr",
        }
        log = Path(f"/root/.pm2/logs/miner-h{i}-error.log")
        events = [re.sub(r"https?://\S+", "[URL redacted]", line)
                  for line in log.read_text(errors="replace").splitlines()
                  if TASK_ID in line or ("2026-09-24 03:35" <= line[:16] <= "2026-09-24 04:50"
                                        and any(w in line for w in ("Build:", "Prefetch:")))]
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
    write("layout_evidence.json", {
        "note": "Historical windows are pinned from this task's build logs and checked by replay. Current layout constants are recorded for context, not substituted for the historical windows.",
        "historical_conjunction_band": [[200, 299]], "historical_conjunction_cut": [[200, 299], [400, 499]],
        "historical_all_hdr_band": [[400, 499]],
        "current_layout_sha256": hashlib.sha256(layout.read_bytes()).hexdigest(),
        "current_layout_constants": [f"{i+1}: {line}" for i, line in enumerate(source.splitlines())
                                     if re.match(r"^(BAND_CLASSES|CUT_CLASSES|BAND_STRIDE|BAND_SUB_WIDTH|HDR_ONLY_CLASSES|HDR_ONLY_SUB_WIDTH|HDR_ONLY_STRIDE)\b", line)],
    })
    shutil.copyfile(ROOT / "reports/f18ee409/validate.py", OUT / "validate.py")
    audit = (ROOT / "reports/f18ee409/audit.py").read_text()
    replacements = {
        "f18ee409-bdb0-4e80-9818-e25f9f456f46": TASK_ID,
        "task f18ee409": "task e6caf784",
        "[214, 249, 784]": "[481, 208, 809]",
        "== {0.35}": "== {0.82}",
        "| k=8": "| k=11",
        "logged == (8, 8, 80, len(ca_clean_in_cut), 200)": "logged == (11, 11, 100, len(ca_clean_in_cut), 200)",
        "logged == (80, len(set(ca_hdr) & set(candidates)), 34)": "logged == (cas_mix['Cas12a'], len(set(ca_hdr) & set(candidates)), 34)",
        "assert len(band) == 8 and set(band) <= set(candidates)": "assert len(band) == (11 if is_conjunction else logged_count) and set(band) <= set(candidates)",
        'assert cas_mix == {"Cas12a": 80, "Cas9": 170}': "assert sum(cas_mix.values()) == 250 and cas_mix['Cas12a'] == (100 if is_conjunction else logged[0])",
        "log clean 8/34": "log clean count/34",
        "those eight seeds": "those HDR seeds",
    }
    for before, after in replacements.items():
        assert before in audit, before
        audit = audit.replace(before, after)
    ast.parse(audit)
    (OUT / "audit.py").write_text(audit)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        assert sys.argv[1:] == ["--render-only"]
        scripts = ["render_report.py"]
    else:
        prepare()
        scripts = ["validate.py", "audit.py", "render_report.py"]
    for name in scripts:
        subprocess.run([str(ROOT / ".venv/bin/python"), str(OUT / name)], cwd=ROOT, check=True)
