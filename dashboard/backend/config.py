"""Settings for the dashboard backend.

The backend does not import niome_subnet. Refreshing runs
``scripts/bench_task.py --fetch`` as a subprocess, so the harness stays the one
implementation of fetching and merging, and its heavy imports stay out of this
process.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

BACKEND_ROOT = Path(__file__).resolve().parent
REPO_ROOT = BACKEND_ROOT.parents[1]

# The snapshot bench_task.py writes and reads. Single source of truth, shared
# with the CLI harness, so a refresh from the dashboard is visible to
# `scripts/bench_task.py` and the other way round.
SNAPSHOT_PATH = Path(os.getenv("DASHBOARD_SNAPSHOT_PATH", REPO_ROOT / "testing" / "task.json"))

# Written by the same --fetch run, alongside the task snapshot.
CELL_TYPES_PATH = Path(
    os.getenv("DASHBOARD_CELL_TYPES_PATH", REPO_ROOT / "testing" / "cell_types.json")
)

BENCH_SCRIPT = Path(os.getenv("DASHBOARD_BENCH_SCRIPT", REPO_ROOT / "scripts" / "bench_task.py"))

# The harness needs its own dependencies, so the repo venv rather than
# whatever interpreter happens to be serving this app.
_default_python = REPO_ROOT / ".venv" / "bin" / "python"
PYTHON = Path(os.getenv("DASHBOARD_PYTHON", _default_python if _default_python.exists() else sys.executable))

# A fetch is network I/O over a paginated endpoint, normally a few seconds.
FETCH_TIMEOUT_SECONDS = float(os.getenv("DASHBOARD_FETCH_TIMEOUT", "120"))

CORS_ORIGINS = [
    origin.strip()
    for origin in os.getenv(
        "DASHBOARD_CORS_ORIGINS", "http://localhost:4200,http://127.0.0.1:4200"
    ).split(",")
    if origin.strip()
]
