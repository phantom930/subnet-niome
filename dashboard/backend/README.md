# Dashboard Backend

FastAPI service that serves the closed-round task snapshot to the frontend and
refreshes it from the public NIOME task history.

It is deliberately standalone. It does not import `niome_subnet` or the
`scripts/` harness, and it never writes to the subnet's `data/` or `testing/`
directories. You can start, restart and deploy it without touching the subnet
code.

## Run

```bash
cd dashboard/backend
uvicorn main:app --reload --port 8000
```

The repo's own `.venv` already has FastAPI, uvicorn and httpx, so you can reuse
it rather than making a second environment:

```bash
../../.venv/bin/uvicorn main:app --reload --port 8000
```

Otherwise `pip install -r requirements.txt`.

The frontend's dev server proxies `/api` here, so start this first. If it is
not running, the Tasks page says so rather than showing an empty table.

## Endpoints

| Method | Path                 | Does                                          |
| ------ | -------------------- | --------------------------------------------- |
| GET    | `/api/health`        | Status, snapshot counts, upstream URL         |
| GET    | `/api/tasks`         | The whole snapshot, newest first              |
| POST   | `/api/tasks/refresh` | Fetch upstream history, merge, report changes |

`POST /api/tasks/refresh` takes `{"replace": false}`. With `replace` true it
discards the stored snapshot instead of merging into it. It answers 409 if a
refresh is already running, and 502 if the upstream is unreachable.

## Where the data lives

Two paths, and only the first is ever written:

1. `dashboard/backend/data/task.json` is the backend's own snapshot. Refresh
   writes here, atomically via a temp file and rename, so a crash cannot leave
   a truncated file. The directory is gitignored.
2. `../../testing/task.json` is a read-only fallback, so the dashboard has
   something to show before the first refresh. The subnet harness owns that
   file and this service never modifies it.

Both are overridable with `DASHBOARD_SNAPSHOT_PATH` and
`DASHBOARD_SEED_SNAPSHOT_PATH`. See `config.py` for the rest, including
`DASHBOARD_UPSTREAM_BASE_URL` and `DASHBOARD_CORS_ORIGINS`.

## Refresh merges rather than overwrites

The upstream history is a window, so a task that has aged out of it would
disappear on a plain overwrite. Refresh merges instead, and fetched entries win
on a shared id, because a round that was unstamped when it was last recorded
carries its real seed now. The response reports both effects:

```json
{ "added": 6, "restamped": 1, "fetched": 441, "count": 441, "unstamped": 25 }
```

`added` is tasks not previously stored, `restamped` is stored tasks whose seed
changed. A refresh that reports `added: 0, restamped: 1` did do something: one
task just got its seed.

## The seed field is mixed-type

The upstream sends the stamped seed as a raw number for most tasks and as a
comma-grouped string for others, and the grouping is not always correct. One
observed value reads `328,371,1000` for what is really 3,283,711,000. A seed of
`0` is the "round has not closed yet" placeholder, not a seed of zero.

`store.normalize_seed` strips commas, parses, and maps zero to `None`, which is
what the `unstamped` count is built from. The frontend normalizes the same way
for display and sorting. Comparing raw seed values instead would both miscount
the unstamped tasks and miss restamps.

## What is not here

Running a benchmark. `scripts/bench_task.py` scores a submission through the
validator's own stages, which needs `niome_subnet`, the 130 MB chromosome, and
the miner's design code. Wiring that into this service would couple the
dashboard to the subnet internals, which is what keeping it standalone avoids.

If you do want a benchmark button later, the cleanest shape is a separate
worker that runs `scripts/bench_task.py` as a subprocess and reports results,
kept apart from this read-only API. Note that the harness mutates process-global
state: `redirect_stage_paths` rewrites the stages' module-level path constants,
so two concurrent runs in one process would return wrong scores rather than
failing. Any such worker has to serialize runs or give each its own process.
