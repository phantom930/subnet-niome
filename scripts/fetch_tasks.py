#!/usr/bin/env python3
"""Snapshot the backend's whole task history into a task file ``scripts/run_task.py`` can read.

``/api/v3/tasks`` returns every task the backend has ever issued in one page, newest first, each one
carrying the contract and HBB reference a validator broadcast **plus the seed the backend stamped
afterwards**. That last part is the reason this script exists: a live task reaches a miner with
``seed: 0``, so the only way to score a build against a real seed is to replay a task that has
already closed.

The newest few rows are exactly the ones that have not closed yet, which is what ``--drop-newest``
is for. Its default of 3 is deliberately generous — a round is 720 blocks and the seed is stamped
somewhere inside it, so the top of the list is a mix of stamped and unstamped and dropping a couple
of good tasks costs nothing against 260 others. Tasks that still read ``seed: 0`` after the drop are
counted in the output (``unstamped``) and left in place rather than filtered: ``run_task.py`` warns
about them on its own, and a zero seed is still a legitimate thing to build against.

    python scripts/fetch_tasks.py                     # -> testing/task.json, newest 3 dropped
    python scripts/fetch_tasks.py --drop-newest 0     # keep everything, unstamped tasks included
    python scripts/fetch_tasks.py --out /tmp/t.json   # somewhere else

The file it writes is ``{"tasks": [...]}`` with the provenance around it, which is one of the shapes
``run_task.py --task-index/--task-id`` accepts; ``genExp.py --all-tasks`` reads the same history
straight off the backend and does not need this file.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from niome_subnet.utils import settings                                       # noqa: E402

TASKS_URL = f"{settings.BASE_URL}/api/v3/tasks"
LEADERBOARD = "https://niome-leaderboard.genomes.io/tasks"
DEFAULT_OUT = REPO_ROOT / "testing" / "task.json"


def fetch(url: str, timeout: float) -> list[dict]:
    """The task list, newest first. Accepts the paged envelope or a bare array."""
    with urllib.request.urlopen(url, timeout=timeout) as response:
        document = json.load(response)
    if isinstance(document, list):
        return document
    for key in ("items", "tasks"):
        if isinstance(document.get(key), list):
            return document[key]
    raise SystemExit(f"{url}: no task array under 'items' or 'tasks' (keys: {list(document)})")


def check_order(tasks: list[dict]) -> bool:
    """Whether the backend really handed them back newest first.

    ``--drop-newest`` cuts from the front, so the wrong order would silently drop the *oldest*
    tasks — the stamped ones — and keep the unstamped ones it exists to remove. Worth asserting.
    """
    stamps = [task.get("created_at") or "" for task in tasks]
    return all(a >= b for a, b in zip(stamps, stamps[1:]))


def seed_of(task: dict) -> object:
    return task.get("content", {}).get("contract", {}).get("seed")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--url", default=TASKS_URL)
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--drop-newest", type=int, default=3,
                        help="skip this many from the front of the list — the rounds still open, "
                             "whose seeds the backend has not stamped yet (default: 3)")
    parser.add_argument("--limit", type=int, default=0,
                        help="keep at most this many tasks after the drop (0 = all)")
    parser.add_argument("--timeout", type=float, default=120.0)
    args = parser.parse_args()

    print(f"[1/3] GET {args.url}")
    tasks = fetch(args.url, args.timeout)
    ordered = check_order(tasks)
    newest = (tasks[0].get("created_at") or "?")[:19] if tasks else "?"
    oldest = (tasks[-1].get("created_at") or "?")[:19] if tasks else "?"
    print(f"  {len(tasks)} tasks, {newest} .. {oldest}"
          f"{'' if ordered else '  ! NOT newest-first'}")
    if not ordered:
        raise SystemExit("  refusing to write: --drop-newest cuts from the front and assumes "
                         "newest-first, so a different order would drop the wrong end")

    kept = tasks[max(0, args.drop_newest):]
    if args.limit:
        kept = kept[:args.limit]
    dropped = [f"{t.get('id', '?')[:8]} ({(t.get('created_at') or '?')[:19]}, seed {seed_of(t)})"
               for t in tasks[:max(0, args.drop_newest)]]
    print(f"[2/3] dropped {len(dropped)} newest, kept {len(kept)}")
    for line in dropped:
        print(f"  - {line}")
    if not kept:
        raise SystemExit("  nothing left to write")

    unstamped = [t.get("id", "?")[:8] for t in kept if not seed_of(t)]
    if unstamped:
        print(f"  ! {len(unstamped)} kept task(s) still read seed 0: {', '.join(unstamped[:8])}"
              f"{' ...' if len(unstamped) > 8 else ''}")

    out = Path(args.out).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    document = {
        "source": args.url,
        "leaderboard": LEADERBOARD,
        "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "order": "newest first",
        "note": (f"snapshot of the backend's task history, newest {max(0, args.drop_newest)} "
                 "dropped as still-open rounds whose seed is not stamped yet. Each task carries "
                 "the contract a validator broadcast plus the seed stamped after that round "
                 "closed, which a live task does not."),
        "count": len(kept),
        "unstamped": len(unstamped),
        "tasks": kept,
    }
    out.write_text(json.dumps(document, indent=2))
    size = out.stat().st_size
    print(f"[3/3] wrote {out} ({len(kept)} tasks, {size / 1000:.0f} kB)")
    print(f"  newest kept {kept[0].get('id', '?')[:8]} "
          f"({(kept[0].get('created_at') or '?')[:19]}, seed {seed_of(kept[0])})")
    print(f"  scripts/run_task.py --task {out} --task-index 0   # or --task-id <prefix>")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
