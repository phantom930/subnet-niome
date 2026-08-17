#!/usr/bin/env python
"""Re-upload a task's submission to the validator's presigned URL, by hand.

The miner gets one shot per task id and no feedback channel, so when an upload fails the only
recovery is to PUT the rows again before the presigned URL expires — 300 s (``SUBMISSION_TIMEOUT``)
from the moment the validator minted it. That window is the whole constraint: this script is worth
running the moment a failure appears in the log, and worth nothing several minutes later.

What it sends is the record ``submission.py`` archived for that exact task id in ``submission.json``
— never a positional "latest" entry, because uploading another task's rows would be scored as this
task's and there is no way to take it back. If the archive has no record for the id, that is an
error with a build command, not a fallback.

The URL comes from ``data/last_upload.json``, which the miner writes when it starts a task: the
signature is deliberately kept out of the log, so that record is the only thing a failed upload
leaves behind to retry with.

Usage, from the repo root:

    python scripts/resubmit.py                          # retry the last task the miner handled
    python scripts/resubmit.py --task-id <uuid>         # pick the record explicitly
    python scripts/resubmit.py --url '<PRESIGNED>' --task-id <uuid>
"""

import argparse
import json
import logging
import os
import sys
import time

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

import niome_subnet.utils.settings as settings  # noqa: E402

from neurons.miner import Miner  # noqa: E402

logger = logging.getLogger("resubmit")

# submission.py's --out default, relative to the repo root like every other path here.
DEFAULT_ARCHIVE = "submission.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--archive",
        default=DEFAULT_ARCHIVE,
        help=f"submission.py's output, holding one record per task (default: {DEFAULT_ARCHIVE}).",
    )
    parser.add_argument(
        "--task-id",
        default=None,
        help="Task whose record to upload. Defaults to the task the miner last handled, per "
        f"{settings.LAST_UPLOAD_PATH}.",
    )
    parser.add_argument(
        "--url",
        default=None,
        help=f"Presigned PUT URL. Defaults to the one recorded in {settings.LAST_UPLOAD_PATH}. "
        "Quote it — the query string contains '&'.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Upload even when the recorded attempt already succeeded, or the URL looks expired. "
        "An expired URL is refused by S3, not by us, so this only costs a failed request.",
    )
    return parser.parse_args()


def load_recorded_upload() -> dict:
    """What the miner recorded about its most recent task, or an empty record."""
    if not os.path.exists(settings.LAST_UPLOAD_PATH):
        return {}
    try:
        with open(settings.LAST_UPLOAD_PATH) as handle:
            return json.load(handle)
    except Exception as e:
        logger.warning(f"Could not read {settings.LAST_UPLOAD_PATH} ({e})")
        return {}


def load_task_record(archive: str, task_id: str) -> dict:
    """The archived record for exactly this task id.

    Matching is on the id and nothing else. The archive is ordered newest-first and spans the
    backend's whole history, so "the last entry" is a different task on most runs — and rows built
    against another contract would still upload cleanly and then score as this task's.
    """
    if not os.path.exists(archive):
        raise SystemExit(
            f"No {archive} to read. Build this task's rows first:\n"
            f"    python submission.py --task-id {task_id}"
        )
    with open(archive) as handle:
        document = json.load(handle)

    records = document.get("tasks") if isinstance(document, dict) else document
    if not isinstance(records, list):
        raise SystemExit(
            f"{archive} is not submission.py output — expected a document with a 'tasks' list."
        )

    matches = [r for r in records if isinstance(r, dict) and r.get("task_id") == task_id]
    if not matches:
        raise SystemExit(
            f"{archive} holds {len(records)} record(s), none for task {task_id}. Build it:\n"
            f"    python submission.py --task-id {task_id}\n"
            "(that writes the same rows the miner would have sent, so the retry is not a guess)"
        )
    if len(matches) > 1:
        # Two records for one id means two different builds; there is no basis for picking one.
        raise SystemExit(
            f"{archive} holds {len(matches)} records for task {task_id} — remove the stale one, "
            "or rebuild the archive."
        )
    return matches[0]


def describe(record: dict) -> None:
    """Log what this record is, so the operator sees what is about to be sent."""
    expected = record.get("expected") or {}
    logger.info(
        f"Record for task {record.get('task_id')} | cell_type={record.get('cell_type')} "
        f"rows={record.get('rows')} construction={record.get('construction')} "
        f"weight_skew={record.get('weight_skew')}"
        + (f" expected_score={expected['final_score']:.3f}" if "final_score" in expected else "")
    )
    if record.get("seed_provisional"):
        logger.warning(
            "This record was built against an unstamped contract (seed 0). The rows stay valid, "
            "but stage 3 is keyed on the seed, so the construction — and the consistency factor "
            "it buys — will not hold if the backend has since restamped the task."
        )
    for problem in record.get("problems") or []:
        logger.warning(f"Record violates a submission invariant: {problem}")


def main() -> int:
    args = parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s | %(levelname)s | %(name)s | %(message)s"
    )

    now = time.time()
    recorded = load_recorded_upload()

    task_id = args.task_id or recorded.get("task_id")
    if not task_id:
        raise SystemExit(
            f"No task id: {settings.LAST_UPLOAD_PATH} is missing or has none, so pass --task-id."
        )

    url = args.url or recorded.get("presigned_url")
    if not url:
        raise SystemExit(
            f"No presigned URL: {settings.LAST_UPLOAD_PATH} is missing or has none. Pass --url "
            "with the URL the validator sent — it cannot be reconstructed, and the log prints it "
            "without its signature on purpose."
        )

    # The recorded outcome only describes the recorded URL; with --url this is someone else's
    # attempt and says nothing about whether that upload landed.
    if not args.url and recorded:
        logger.info(
            f"Miner's last upload for task {task_id}: "
            f"{'succeeded' if recorded.get('submitted') else 'FAILED or never completed'}"
            + (f", {recorded['rows']} rows" if recorded.get("rows") else "")
        )
        if recorded.get("submitted") and not args.force:
            logger.info(
                "That upload already succeeded, so there is nothing to recover. Re-run with "
                "--force to overwrite the object anyway."
            )
            return 0

    deadline = (
        recorded.get("expires_at")
        if not args.url and recorded.get("expires_at")
        else Miner._upload_deadline(url, now)
    )
    remaining = deadline - now
    logger.info(f"Target {Miner._url_summary(url)} | {remaining:.0f}s of TTL left")
    if remaining <= 0 and not args.force:
        logger.error(
            "The presigned URL has already expired — S3 will reject this and the round is lost; "
            "the next task is the next opportunity. Re-run with --force to try regardless."
        )
        return 1

    record = load_task_record(args.archive, task_id)
    describe(record)
    rows = record.get("submission")
    if not rows:
        raise SystemExit(
            f"The record for task {task_id} in {args.archive} carries no rows — rebuild it with "
            f"`python submission.py --task-id {task_id}`."
        )

    logger.info(f"Uploading {len(rows)} rows from {args.archive}")

    # Miner._upload is a plain method that touches no instance state, so an uninitialised Miner is
    # enough to borrow it — and borrowing it is the point: the retry then carries the same headers,
    # retry count and deadline handling as the miner's own upload, with no second copy to drift.
    # The headers in particular are load-bearing; see Miner._upload_headers.
    miner = Miner.__new__(Miner)
    try:
        # --force past the deadline still needs a window to make the attempt in.
        Miner._upload(miner, url, rows, max(deadline, now + 30) if args.force else deadline)
    except Exception as e:
        logger.error(f"Manual resubmit failed: {e}")
        return 1

    logger.info(f"Resubmitted {len(rows)} rows for task {task_id}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
