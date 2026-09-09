# Dashboard Frontend

Angular 22 single-page app with [PrimeNG](https://primeng.org) 22, initialized
and ready for features.

## Requirements

- Node.js 22 (built and tested on v22.23.2)
- npm 12 for a fresh dependency resolve, or npm 10 with `npm ci`

## Install and run

```bash
cd dashboard/frontend
npm ci
npm start
```

The dev server listens on <http://localhost:4200>. The port is set on the serve
target in `angular.json`; override it per run with `npm start -- --port 9001`.

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

| Command           | What it does                                 |
| ----------------- | -------------------------------------------- |
| `npm start`       | Dev server on port 4200 with hot reload      |
| `npm run build`   | Production build into `dist/frontend`        |
| `npm run watch`   | Development build, rebuilding on change      |
| `npm test`        | Unit tests (Vitest); watches in a TTY        |
| `npm run test:ci` | Unit tests once, no watch                    |

## Layout

```text
src/
  app/
    app.ts / app.html / app.scss   Root shell: toolbar and theme toggle
    app.config.ts                  Providers, including providePrimeNG()
    app.routes.ts                  Empty; add feature routes here
    core/theme.service.ts          Light/dark mode, persisted per browser
    theme/app-preset.ts            PrimeNG theme preset (Aura + emerald)
  styles.scss                      Global styles, CSS layer order, PrimeIcons
```

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
