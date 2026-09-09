# Dashboard Frontend

Angular 22 single-page app with [PrimeNG](https://primeng.org) 22. Browses the
closed-round task snapshot, shows any task's raw JSON, runs benchmarks against
individual tasks, and documents the mining process.

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
    pages/mining/mining.*          The mining process, step by step
    pages/mining/mining.data.ts    Its content, with a source file per step
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

The page itself is one step off the card background rather than equal to it, so
cards, the toolbar and the tables read as surfaces laid on the page instead of
one flat sheet. `styles.scss` holds that as `--app-background`, which the body
and `.app-shell` both use: `--p-surface-50` in light, `--p-surface-950` in dark.
Light stops at 50 rather than 100 because Aura gives a secondary `p-message` a
surface-100 background and those callouts sit directly on the page. Each theme
therefore keeps three distinct levels: page, card, callout.

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

## The Mining page

A step-by-step account of what happens between a validator's broadcast and a
score: the request path, the design pipeline inside the build, what each
scoring factor responds to, the contract's obligations, and the rules that cost
rows without reporting anything.

The prose lives in `mining.data.ts` rather than in the template, so there is one
place to edit, and every step names the file it was drawn from. The sources are
`docs/miner_flow.md`, `docs/validation_pipeline.md`, `neurons/miner.py`,
`niome_subnet/genomics/design.py`, `niome_subnet/validator/forward.py` and the
validation stages themselves. Numeric results quoted from those docs are
offline measurements, and the page says so rather than presenting them as
guarantees.

## The three figures

Three things on that page are pictures because the prose version of them needs
the reader to hold several parts in mind at once:

| Figure                            | What it shows                                                               |
| --------------------------------- | --------------------------------------------------------------------------- |
| The round trip, in the first card | That the score comes from an S3 object, not from the reply to `/forward`    |
| The pipeline, above the factors   | Which stage emits which factor, and which half of it the round seed reaches |
| Occupancy, in the coverage card   | Which kind of missing coverage is a haircut and which is a cliff            |

They are hand-authored inline SVG in `mining.html`: no chart library, no
runtime, no images. Five conventions hold across all three, and are worth
keeping if you add a fourth figure:

- **Each `<svg>` carries `width` and `height` attributes matching its
  `viewBox`.** That is what makes one user unit render as one CSS pixel, so the
  font sizes in `mining.scss` are the sizes that reach the screen. Left to
  stretch, an svg fills the card and scales its type up with it, which on this
  1200px layout renders a 12px label at 20px. The single `.dg` rule handles
  the rest: `max-width: 100%` to shrink on a narrow screen, `min-width` to
  scroll instead of shrinking past legibility.
- **One hue per kind of component, shared across the figures.** Indigo is the
  validator and the stages it runs, amber the bucket, sky the miner and the
  array it builds, emerald the score and its three factors. A box wears its
  role by class: `class="dg-box dg-validator"`. The role sets `--dg-hue` and
  the fill follows from it as `color-mix(… var(--dg-hue) 20%,
var(--p-content-background))`, so neither theme is hand-coded: the same
  declaration is a pale wash on a light card and a deep one on a dark card. Add
  a role by adding one line. The colour is always a shortcut, never the only
  clue, since every box is labelled and the first figure's caption states the
  key.
- Everything unroled stays `currentColor` so it falls out of the card's own
  foreground, and `--p-red-500` marks the one number that is a cliff.
  `mining.scss` groups these rules by what a mark means rather than restating
  them per class, so the ordering comment there matters: the overrides have to
  land after what they override.
- Each figure is a `<figure>` with a `<figcaption>` carrying its claim, and the
  `<svg>` repeats that claim in `aria-label` for anyone who cannot see it. A
  test asserts all three still have both.
- Geometry that needs arithmetic is computed in `mining.ts`, because Angular
  template expressions have no floor or bitwise operator. Labels that cross a
  line knock a gap in it via `paint-order`, rather than being nudged by hand.

The diagram rules take `mining.scss` past the 4 kB `anyComponentStyle` budget
that `angular.json` shipped with, so the warning threshold there is now 5 kB.
The 8 kB error is untouched and is still the real guard.

The occupancy figure is derived, not drawn: `mining.ts` ports stage 5's
`coverage_entropy_ratio` and `geometric_mean`, and the factor under each grid
comes out of them. That is how the figure came to correct the page. Stage 5
floors each of the six ratios at 1e-9 and then takes a sixth root, so a zeroed
ratio caps the factor at 0.032 rather than zeroing it, and one empty cell out
of eight does not zero any ratio at all: it costs about two percent. The 1e-9
cliff belongs to a collapsed dimension, every row on one Cas system, one strand
or one mutation. `design.py`'s own docstring and `docs/miner_flow.md` attribute
it to a single empty cell, which the arithmetic does not support.

One section is computed rather than written. The accessibility table is fetched
from `/api/cell-types` and run through stage 3's own formulas to show, per cell
type, the energy and the Cas9 and Cas12a cut probabilities:

```text
energy          = clamp(0, 1, accessibility x (1.8*gc + 0.6*exp(-dist/1500) + offset))
cut_probability = clamp(0.4, 0.99, base + 0.18 * energy)    base 0.86 Cas9, 0.78 Cas12a
```

Substituting the miner's own design, GC at 50 percent and the nearest PAM, puts
the inner term near 1.5, so energy saturates once accessibility passes about
0.67. That is worth seeing live rather than describing, because three of the
four cell types the backend currently issues are already past it, which changes
what `consistency_factor` can reach. The tests pin the derived numbers against
the formulas so this section cannot drift from the code it describes.
