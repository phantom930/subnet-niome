# Dashboard Backend

FastAPI service that serves the closed-round task snapshot to the frontend and
refreshes it by running `scripts/bench_task.py --fetch`.

The harness stays the one implementation of fetching and merging. This service
does not reimplement it and does not import `niome_subnet`: it shells out, so
the snapshot in `testing/` is a single source of truth shared with the CLI.

## Run

```bash
cd dashboard/backend
../../.venv/bin/uvicorn main:app --reload --port 8000
```

The repo's `.venv` already has FastAPI and uvicorn. The harness needs the
project's own dependencies, so `config.PYTHON` defaults to `../../.venv/bin/python`
regardless of which interpreter is serving this app. Override with
`DASHBOARD_PYTHON` if your environment lives elsewhere.

The frontend's dev server proxies `/api` here, so start this first. If it is
not running, the Tasks page says so rather than showing an empty table.

## Endpoints

| Method | Path                 | Does                                        |
| ------ | -------------------- | ------------------------------------------- |
| GET    | `/api/health`        | Status, snapshot counts, resolved paths     |
| GET    | `/api/tasks`         | The whole snapshot, newest first            |
| GET    | `/api/cell-types`    | The accessibility table                     |
| POST   | `/api/tasks/refresh` | Runs `--fetch`, reports what changed        |

`POST /api/tasks/refresh` takes `{"replace": false}`; `replace` passes
`--replace`, discarding the local snapshot instead of merging into it. It
answers 409 if a refresh is already running, 502 if the harness fails, and 504
if it outruns `DASHBOARD_FETCH_TIMEOUT` (120s default).

## Why a subprocess and not an import

`bench_task` pulls in `niome_subnet`, the validation stages and a 130 MB
chromosome, and it mutates module-level state: `redirect_stage_paths` rewrites
the stages' path constants with `setattr`, and `memoise_chr11` monkey-patches
`stage12.load_chr11` with a closure holding its own cache. Importing that into
a long-lived web process would drag all of it in and make concurrent operations
share the global state. A subprocess gets its own copy and exits when done.

A refresh is also serialized behind an asyncio lock, because two `--fetch` runs
would both write the snapshot and the second merge would be built from a file
the first was still replacing.

## What a refresh touches

`--fetch` writes two files, both under `testing/`:

- `testing/task.json`, the merged task snapshot
- `testing/cell_types.json`, the accessibility table

Both are the harness's own paths, so a refresh from the dashboard is visible to
`scripts/bench_task.py` and the other way round. Override with
`DASHBOARD_SNAPSHOT_PATH` and `DASHBOARD_CELL_TYPES_PATH` if you want the
dashboard to keep a separate copy instead.

## The counts come from the files, not from stdout

The endpoint reads the snapshot's seeds before the run and diffs them against
the snapshot afterwards, so a change to the harness's output format cannot
silently break the numbers. The script's own summary is passed through in
`output` for the details it prints, including its per-page fetch progress and
any cell types in the history that are missing from the table.

```json
{ "added": 0, "restamped": 1, "removed": 0, "count": 441, "unstamped": 25 }
```

A refresh reporting `added: 0, restamped: 1` did do something: one task just
got its seed. That is why `restamped` is reported separately rather than folded
into a single "changed" count.

## The seed field is mixed-type

The upstream sends the stamped seed as a raw number for most tasks and as a
comma-grouped string for others, and the grouping is not always correct. One
observed value reads `328,371,1000` for what is really 3,283,711,000. A seed of
`0` is the "round has not closed yet" placeholder, not a seed of zero.

`store.normalize_seed` strips commas, parses, and maps zero to `None`. The diff
compares normalized seeds, because the same seed arriving as `100000` on one
fetch and `100,000` on the next would otherwise be reported as a restamp that
never happened. The frontend normalizes the same way for display and sorting.

## What is not here

A benchmark button. `--fetch` is network I/O that finishes in seconds, so it
runs inline. A default `bench_task.py` run scores a submission through all five
validation stages, which takes roughly 7 seconds for one seed and longer for
more, and a synchronous request would time out.

The shape for that is a job API rather than a blocking call: POST starts a run
and returns an id, GET polls for the result. The subprocess approach here
carries over unchanged, and it already handles the global-state problem that
makes concurrent runs unsafe.
