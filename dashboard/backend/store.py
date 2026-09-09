"""Reading, merging and writing the task snapshot.

The snapshot shape matches what the frontend already consumes:

    {source, leaderboard, fetched_at, order, note, count, unstamped, tasks}
"""

from __future__ import annotations

import json
import os
import tempfile
from datetime import datetime, timezone
from typing import Any

import config

ORDER = "newest first"
NOTE = (
    "Closed-round snapshot: each task carries the contract a validator broadcast plus the "
    "seed stamped after that round closed. Merged across refreshes, so tasks that have aged "
    "out of the backend's window are kept."
)


class SnapshotMissing(RuntimeError):
    """Neither the backend's snapshot nor the read-only fallback exists."""


def normalize_seed(seed: Any) -> int | None:
    """The stamped seed as a number, or None when the round never closed.

    The upstream sends this field as both a number and a comma-grouped string,
    and the grouping is not always correct, so commas come out before parsing.
    Zero is the 'not stamped yet' placeholder rather than a seed of zero.
    """
    if seed is None or seed == "":
        return None
    try:
        value = seed if isinstance(seed, (int, float)) else float(str(seed).replace(",", ""))
    except (TypeError, ValueError):
        return None
    if value == 0:
        return None
    return int(value)


def count_unstamped(tasks: list[dict]) -> int:
    return sum(
        1
        for task in tasks
        if normalize_seed(task.get("content", {}).get("contract", {}).get("seed")) is None
    )


def _sort_key(task: dict) -> str:
    return str(task.get("created_at", ""))


def build_snapshot(tasks: list[dict], fetched_at: str | None = None) -> dict[str, Any]:
    ordered = sorted(tasks, key=_sort_key, reverse=True)
    return {
        "source": config.TASK_HISTORY_URL,
        "leaderboard": config.LEADERBOARD_URL,
        "fetched_at": fetched_at or datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S%z"),
        "order": ORDER,
        "note": NOTE,
        "count": len(ordered),
        "unstamped": count_unstamped(ordered),
        "tasks": ordered,
    }


def read() -> dict[str, Any]:
    """The backend's snapshot, falling back to the harness's read-only copy."""
    for path in (config.SNAPSHOT_PATH, config.SEED_SNAPSHOT_PATH):
        if path.exists():
            with open(path) as snapshot_file:
                return json.load(snapshot_file)
    raise SnapshotMissing(
        f"no snapshot at {config.SNAPSHOT_PATH} and no fallback at "
        f"{config.SEED_SNAPSHOT_PATH}. Refresh to fetch one."
    )


def read_tasks() -> list[dict]:
    try:
        return read().get("tasks", [])
    except (SnapshotMissing, ValueError):
        return []


def summary() -> dict[str, Any]:
    """Counts without shipping the whole task array."""
    try:
        snapshot = read()
    except (SnapshotMissing, ValueError):
        return {"count": 0, "unstamped": 0, "fetched_at": None}
    tasks = snapshot.get("tasks", [])
    return {
        "count": snapshot.get("count", len(tasks)),
        "unstamped": snapshot.get("unstamped", count_unstamped(tasks)),
        "fetched_at": snapshot.get("fetched_at"),
    }


def merge(fetched: list[dict], replace: bool = False) -> dict[str, Any]:
    """Fold freshly fetched tasks into what is already stored.

    Merging rather than overwriting is deliberate: the upstream history is a
    window, so a task that has aged out of it stays available instead of
    disappearing. Fetched entries win on a shared id, because a round that was
    unstamped when it was last recorded carries its real seed now.
    """
    existing = [] if replace else read_tasks()
    by_id = {task["id"]: task for task in existing if "id" in task}

    added = 0
    restamped = 0
    for task in fetched:
        task_id = task.get("id")
        if task_id is None:
            continue
        previous = by_id.get(task_id)
        if previous is None:
            added += 1
        else:
            before = previous.get("content", {}).get("contract", {}).get("seed")
            after = task.get("content", {}).get("contract", {}).get("seed")
            if normalize_seed(before) != normalize_seed(after):
                restamped += 1
        by_id[task_id] = task

    snapshot = build_snapshot(list(by_id.values()))
    write(snapshot)
    return {
        "added": added,
        "restamped": restamped,
        "fetched": len(fetched),
        "count": snapshot["count"],
        "unstamped": snapshot["unstamped"],
        "fetched_at": snapshot["fetched_at"],
        "replaced": replace,
    }


def write(snapshot: dict[str, Any]) -> None:
    """Write the snapshot atomically, so a crash cannot truncate it."""
    config.SNAPSHOT_PATH.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary_path = tempfile.mkstemp(
        dir=config.SNAPSHOT_PATH.parent, prefix=".task-", suffix=".json"
    )
    try:
        with os.fdopen(handle, "w") as snapshot_file:
            json.dump(snapshot, snapshot_file)
        os.replace(temporary_path, config.SNAPSHOT_PATH)
    except BaseException:
        os.unlink(temporary_path)
        raise
