"""Settings for the dashboard backend.

Standalone on purpose: this service does not import niome_subnet or the
scripts/ harness, so it can be run, restarted and deployed without touching
the subnet code or its data/ directory.
"""

from __future__ import annotations

import os
from pathlib import Path

BACKEND_ROOT = Path(__file__).resolve().parent
REPO_ROOT = BACKEND_ROOT.parents[1]

# The backend's own snapshot. Refresh writes here and nowhere else.
SNAPSHOT_PATH = Path(
    os.getenv("DASHBOARD_SNAPSHOT_PATH", BACKEND_ROOT / "data" / "task.json")
)

# Read-only fallback so the dashboard has something to show before the first
# refresh. Never written to; this is the subnet harness's own snapshot.
SEED_SNAPSHOT_PATH = Path(
    os.getenv("DASHBOARD_SEED_SNAPSHOT_PATH", REPO_ROOT / "testing" / "task.json")
)

# Closed-round history. bench_task.py documents this as the one task endpoint
# that needs no hotkey signature: a plain public read that sends nothing
# identifying. /current, by contrast, rejects an unsigned GET.
UPSTREAM_BASE_URL = os.getenv("DASHBOARD_UPSTREAM_BASE_URL", "https://niome-api.genomes.io")
TASK_HISTORY_URL = f"{UPSTREAM_BASE_URL}/api/v3/tasks"
LEADERBOARD_URL = os.getenv(
    "DASHBOARD_LEADERBOARD_URL", "https://niome-leaderboard.genomes.io/tasks"
)

UPSTREAM_TIMEOUT_SECONDS = float(os.getenv("DASHBOARD_UPSTREAM_TIMEOUT", "30"))
UPSTREAM_MAX_PAGES = int(os.getenv("DASHBOARD_UPSTREAM_MAX_PAGES", "100"))

CORS_ORIGINS = [
    origin.strip()
    for origin in os.getenv(
        "DASHBOARD_CORS_ORIGINS", "http://localhost:4200,http://127.0.0.1:4200"
    ).split(",")
    if origin.strip()
]
