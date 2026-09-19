#!/usr/bin/env python3
"""seed_watch.py — when does a task's contract actually carry its REAL seeds, and on which source?

The whole question for a seed-aware build is timing: our presigned upload URL lives 300s from the
validator's call, and calls land at a median of +30 min after task creation (p90 +52 min, n=946).
So a seed-pinned build is only possible if the real seeds are readable BEFORE that URL dies. Two
sources disagree about when that happens and neither is documented:

  API  `/api/v3/tasks` -> content.contract.seed      -- public, unsigned, pollable at any time
  S3   the per-task contract.json the validator hands us at call time (presigned; 403 otherwise)

CLAUDE.md claims the seed is "0 on both sides until the round closes". Two live observations on
2026-09-18 do NOT fit that cleanly and are the reason this exists:

  * h7 was called at +145 min and its S3 contract carried `448,931,493`, while every hotkey called
    between +2 and +56 min that same round got `seed: 0`. So S3 flips somewhere in +56..+145 min.
  * The API reported `448,931,493` for that task at +232 min and `341,410,564` at +246 min — the
    SAME task_id and created_at, a CHANGED seed. Stable across 8 polls 4s apart, so it is not
    per-request noise; something rewrote it. Only the newest task moved; older ones were frozen.

That second point matters more than the first: a seed that is rewritten after we read it cannot be
pinned against, however early it appears. This logs both sources on a fixed cadence so the flip
time, the stability, and any API-vs-S3 lead can be read off rather than inferred from two samples.

Writes one JSON object per poll to `seed_watch.jsonl` (append-only, safe to tail while running):

    ts, task, created_at, age_min, cell, api_seed, api_changed, s3_status, s3_seed

Read-only: it fetches public endpoints and writes one local file. It never touches miner state.

    SW_EVERY=60 SW_HOURS=12 python seed_watch.py
"""
import os
import sys

sys.argv = ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")

import datetime
import json
import time

import requests                                          # noqa: E402

from conj_stageb import API                               # noqa: E402

EVERY = float(os.getenv("SW_EVERY", "60"))
HOURS = float(os.getenv("SW_HOURS", "12"))
NTASK = int(os.getenv("SW_NTASK", "3"))       # newest N tasks tracked each poll
OUT = os.getenv("SW_OUT", "seed_watch.jsonl")
S3 = "https://niome-data-bucket.s3.amazonaws.com/CRISPR/{tid}/contract.json"


def utcnow():
    return datetime.datetime.now(datetime.UTC).replace(tzinfo=None)


def parse_ts(s):
    try:
        return datetime.datetime.fromisoformat(str(s).replace("Z", "")[:26])
    except Exception:
        return None


def is_stamped(seed):
    """A real 3-seed stamp, as opposed to the 0 / empty placeholder a broadcast carries."""
    parts = [p for p in str(seed or "").split(",") if p.strip().isdigit()]
    return len(parts) >= 2 or (len(parts) == 1 and int(parts[0]) != 0)


def main():
    deadline = time.monotonic() + HOURS * 3600
    last = {}          # task -> last api seed seen, for change detection
    first_stamp = {}   # task -> (age_min, seed) when it first read as stamped
    print(f"polling every {EVERY:.0f}s for {HOURS:.1f}h -> {OUT}", flush=True)
    while time.monotonic() < deadline:
        now = utcnow()
        try:
            r = requests.get(f"{API}/tasks?limit=8", timeout=30)
            items = r.json()
            items = items if isinstance(items, list) else (items.get("items")
                                                           or items.get("data") or [])
        except Exception as exc:
            print(f"{now.isoformat()[:19]} API error {exc}", flush=True)
            time.sleep(EVERY)
            continue

        items = sorted(items, key=lambda t: t.get("created_at", ""), reverse=True)[:NTASK]
        for t in items:
            tid = t.get("task_id") or t["id"]
            c = (t.get("content") or {}).get("contract") or {}
            api_seed = c.get("seed")
            ca = parse_ts(t.get("created_at"))
            age = (now - ca).total_seconds() / 60 if ca else None

            s3_status, s3_seed = None, None
            try:
                rr = requests.get(S3.format(tid=tid), timeout=15)
                s3_status = rr.status_code
                if rr.ok:
                    s3_seed = rr.json().get("seed")
            except Exception:
                s3_status = -1

            changed = tid in last and str(last[tid]) != str(api_seed)
            if changed:
                print(f"  !! {tid[:8]} API seed CHANGED {last[tid]!r} -> {api_seed!r} "
                      f"at age {age:.1f} min", flush=True)
            if is_stamped(api_seed) and tid not in first_stamp:
                first_stamp[tid] = (age, api_seed)
                print(f"  ** {tid[:8]} API first reads STAMPED at age {age:.1f} min "
                      f"-> {api_seed!r}", flush=True)
            last[tid] = api_seed

            rec = {"ts": now.isoformat()[:19], "task": tid[:8], "created_at": str(ca)[:19],
                   "age_min": None if age is None else round(age, 1),
                   "cell": c.get("cell_type"), "api_seed": api_seed,
                   "api_stamped": is_stamped(api_seed), "api_changed": changed,
                   "s3_status": s3_status, "s3_seed": s3_seed,
                   "s3_stamped": is_stamped(s3_seed) if s3_seed is not None else None}
            with open(OUT, "a") as fh:
                fh.write(json.dumps(rec) + "\n")

        newest = items[0] if items else None
        if newest:
            c = (newest.get("content") or {}).get("contract") or {}
            ca = parse_ts(newest.get("created_at"))
            age = (now - ca).total_seconds() / 60 if ca else -1
            print(f"{now.isoformat()[:19]}  newest {(newest.get('task_id') or newest['id'])[:8]} "
                  f"age {age:6.1f}m  api_seed {str(c.get('seed'))!r:22} "
                  f"stamped={is_stamped(c.get('seed'))}", flush=True)
        time.sleep(EVERY)


if __name__ == "__main__":
    main()
