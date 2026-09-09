"""Reading the task snapshot that scripts/bench_task.py writes.

Read-only. The harness owns the file: it fetches, merges and writes, and this
module only reads what is there and diffs two reads to describe what a refresh
changed.
"""

from __future__ import annotations

import json
from typing import Any

import config


class SnapshotMissing(RuntimeError):
    """No snapshot yet. A refresh creates one."""


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


def read() -> dict[str, Any]:
    """The snapshot as the harness wrote it."""
    if not config.SNAPSHOT_PATH.exists():
        raise SnapshotMissing(
            f"no snapshot at {config.SNAPSHOT_PATH}. Refresh to fetch one, or run "
            "scripts/bench_task.py --fetch"
        )
    with open(config.SNAPSHOT_PATH) as snapshot_file:
        return json.load(snapshot_file)


def read_cell_types() -> dict[str, Any]:
    """The accessibility table, written by the same --fetch run."""
    if not config.CELL_TYPES_PATH.exists():
        return {}
    try:
        with open(config.CELL_TYPES_PATH) as cell_types_file:
            return json.load(cell_types_file)
    except (OSError, ValueError):
        return {}


def seeds_by_id() -> dict[str, int | None]:
    """Each task's normalized seed, for diffing across a refresh.

    Normalized rather than raw, because the same seed arrives as 100000 and as
    '100,000' on different tasks. Comparing raw values would report restamps
    that never happened.
    """
    try:
        snapshot = read()
    except (SnapshotMissing, ValueError):
        return {}
    return {
        task["id"]: normalize_seed(task.get("content", {}).get("contract", {}).get("seed"))
        for task in snapshot.get("tasks", [])
        if "id" in task
    }


def summary() -> dict[str, Any]:
    """Counts without shipping the whole task array."""
    # The cell-type table is written by the same --fetch run but is a separate
    # file, so it is reported whether or not the task snapshot is there.
    cell_types = len(read_cell_types())
    try:
        snapshot = read()
    except (SnapshotMissing, ValueError):
        return {"count": 0, "unstamped": 0, "fetched_at": None, "cell_types": cell_types}
    tasks = snapshot.get("tasks", [])
    return {
        "count": snapshot.get("count", len(tasks)),
        "unstamped": snapshot.get(
            "unstamped",
            sum(
                1
                for task in tasks
                if normalize_seed(task.get("content", {}).get("contract", {}).get("seed")) is None
            ),
        ),
        "fetched_at": snapshot.get("fetched_at"),
        "cell_types": cell_types,
    }


def describe_change(before: dict[str, int | None]) -> dict[str, Any]:
    """What changed between a snapshot read earlier and the one on disk now.

    Computed from the files rather than parsed out of the harness's stdout, so
    a change to its output format cannot silently break these numbers.
    """
    after = seeds_by_id()
    added = [task_id for task_id in after if task_id not in before]
    restamped = [
        task_id
        for task_id, seed in after.items()
        if task_id in before and before[task_id] != seed
    ]
    removed = [task_id for task_id in before if task_id not in after]

    return {
        "added": len(added),
        "restamped": len(restamped),
        "removed": len(removed),
        **summary(),
    }
