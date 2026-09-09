# Dashboard Frontend

Angular 22 single-page app with [PrimeNG](https://primeng.org) 22. Browses the
closed-round task snapshot, shows any task's raw JSON, and runs benchmarks
against individual tasks.

## Requirements

- Node.js 22 (built and tested on v22.23.2)
- npm 12 for a fresh dependency resolve, or npm 10 with `npm ci`

## Install and run

Start the backend first, since the dev server proxies `/api` to it:

```bash
cd dashboard/backend
../../.venv/bin/uvicorn main:app --reload --port 8000
```

Then, in a second terminal:

```bash
cd dashboard/frontend
npm ci
npm start
```

The dev server listens on <http://localhost:4200>. The port is set on the serve
target in `angular.json`; override it per run with `npm start -- --port 9001`.
If the backend is not running, the Tasks page says so instead of showing an
empty table.

Use `npm ci` if your global npm is v10. A fresh `npm install` under npm 10.9.8
crashes with `Cannot read properties of null (reading 'edgesOut')` while
resolving the Vitest peer graph. That is an npm bug, not a problem with these
dependencies, and it affects the stock Angular 22 scaffold too. `npm ci`
installs straight from `package-lock.json` and is unaffected. To re-resolve
dependencies when adding or upgrading a package, use npm 12:

```bash
npx npm@12 install
```

## Commands

| Command           | What it does                            |
| ----------------- | --------------------------------------- |
| `npm start`       | Dev server on port 4200 with hot reload |
| `npm run build`   | Production build into `dist/frontend`   |
| `npm run watch`   | Development build, rebuilding on change |
| `npm test`        | Unit tests (Vitest); watches in a TTY   |
| `npm run test:ci` | Unit tests once, no watch               |

## Layout

```text
proxy.conf.json                    Forwards /api to the backend on port 8000
src/
  app/
    app.ts / app.html / app.scss   Root shell: toolbar, nav, theme toggle
    app.config.ts                  Providers, including providePrimeNG()
    app.routes.ts                  Lazy feature routes
    core/task.models.ts            Shapes of the API payloads and table row
    core/task.service.ts           Calls the backend, flattens the snapshot
    core/benchmark.service.ts      Starts a benchmark job and polls it
    core/theme.service.ts          Light/dark mode, persisted per browser
    pages/tasks/tasks.*            Paginated task table
    pages/tasks/benchmark-dialog.* Per-row benchmark result popup
    pages/tasks/task-json-dialog.* Per-row raw JSON viewer
    theme/app-preset.ts            PrimeNG theme preset (Aura + emerald)
  styles.scss                      Global styles, CSS layer order, PrimeIcons
```

## Task data

`TaskService` calls the backend in `dashboard/backend` and flattens each task
into one table row. The dev server proxies `/api` there via `proxy.conf.json`,
so both run same-origin in development and no CORS is involved.

The Refresh button posts to `/api/tasks/refresh`, which makes the backend run
`scripts/bench_task.py --fetch` as a subprocess. That writes `testing/task.json`
and `testing/cell_types.json`, so a refresh from the dashboard is visible to the
CLI harness and the other way round.

The toast reports how many tasks were added and how many were newly stamped,
then the table reloads. Both are reported because a refresh that adds nothing
can still have done something: a task that was unstamped last time carries its
real seed now.

Two things about the data are worth knowing, both handled in `TaskService` and
mirrored in the backend:

- **The seed is mixed-type.** Most tasks send it as a number, some as a
  comma-grouped string, and the grouping is not always right. One observed
  value reads `328,371,1000`. Seeds are normalized to numbers so the column
  sorts correctly, with the original kept for the cell's tooltip.
- **A seed of `0` means unstamped**, not a seed of zero. Those rows show a tag
  instead of a value, which matches the `unstamped` count the backend reports.

The five rule fields are identical across every task, so they appear once above
the table rather than as five identical columns. If they ever stop being
uniform, the panel hides itself and the rules belong back in the rows.

## How PrimeNG is wired

Three pieces have to agree with each other:

1. `app.config.ts` calls `providePrimeNG()` with the theme preset, sets
   `darkModeSelector: '.app-dark'`, and declares the CSS layer order.
2. `styles.scss` repeats that layer order in an `@layer` statement and imports
   the PrimeIcons stylesheet. The statement runs before PrimeNG injects its
   theme at runtime, which is what lets app styles override component styles
   without specificity hacks.
3. `theme.service.ts` toggles the `.app-dark` class on `<html>`, the selector
   configured in step 1.

To re-brand, edit the `primary` palette in `src/app/theme/app-preset.ts`, or
pass `Aura` directly to `providePrimeNG()` for the stock theme.

Components are imported per feature, not globally. Add the module you need to a
component's `imports` array, for example `import { TableModule } from
'primeng/table'`.

## Notes for PrimeNG 22

Two API changes bite when copying older PrimeNG examples:

- `styleClass` was removed from most components, including Card, Toolbar,
  Select, Tag, Skeleton, and ProgressBar. Use a plain `class` attribute. Button
  still accepts `styleClass`, and Table still accepts `tableStyleClass`.
- Content templates are matched by template ref, not by the `pTemplate`
  directive. Write `<ng-template #header>` and `<ng-template #body let-row>`
  rather than `<ng-template pTemplate="header">`.

`provideHttpClient(withFetch())` is already configured for when you start
calling an API.

## Benchmarking a task

Each row has a play button that runs `scripts/bench_task.py` against that task
and shows the result in a dialog. The backend runs it as a subprocess and the
dialog polls until it finishes, roughly 7 seconds for one seed and 15 for the
default three.

The dialog leads with the final score and the product it comes from, then the
validator's fields, then a row per seed when more than one was used. The
harness's full printed report is behind "Show full report".

The numbers are parsed out of that report, because `bench_task.py` has no JSON
mode. If its format ever changes, parsed fields go missing rather than wrong,
and the raw report is still there.

`BenchmarkDialog` starts its run from its `row` input rather than from a method
the parent calls, so there is no ordering to get right between binding the row
and starting the job.

## Viewing a task as JSON

The code button on each row opens the task exactly as the backend sent it,
pretty-printed, with a toggle between the whole task and the contract alone,
and a copy button.

Nothing is refetched for this. `TaskRow` keeps the raw task, which costs
nothing because it is already loaded: a task is about 1.2 kB, so all of them
together are the half-megabyte the snapshot response already carried.

The JSON is rendered as text, not markup. Syntax highlighting would mean
building HTML out of backend data, and that injection risk is not worth it for
a viewer.

Worth knowing when comparing the two views: the JSON shows the seed exactly as
sent, so a task whose seed arrived as the string `546,343,346` reads that way
here while the table column shows the normalized number. A test pins that
difference so neither view drifts.
