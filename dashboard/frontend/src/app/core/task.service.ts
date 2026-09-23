import { HttpClient, HttpErrorResponse } from '@angular/common/http';
import { Injectable, inject } from '@angular/core';
import { Observable, catchError, map, of, throwError } from 'rxjs';

import {
  CellType,
  MutationEntry,
  RawTask,
  RefreshResult,
  SeedOccurrence,
  SnapshotMeta,
  TaskRow,
  TaskRules,
  TaskSnapshot,
} from './task.models';

export interface TaskData {
  meta: SnapshotMeta;
  rows: TaskRow[];
}

/**
 * The round's stamped seeds, or [] when it never closed.
 *
 * The comma separates seeds; it is not digit grouping. A round has carried
 * three seeds since 2026-08-27 and one before that, and the backend joins them
 * ('263,486,269'), so stripping the commas and parsing one number turned three
 * seeds into 263486269 — which the number pipe then re-rendered with the same
 * commas, hiding the mistake on every round whose seeds were all three digits.
 * A seed of 1000 is what exposed it: '328,371,1000' surfaced as 3,283,711,000.
 *
 * Zero is the 'not stamped yet' placeholder rather than a seed, so it drops
 * out, and an unparseable field yields [] rather than a partial round.
 */
function parseSeeds(seed: string | number | null | undefined): number[] {
  if (seed === null || seed === undefined || seed === '') return [];
  const parts = String(seed).split(',');
  const seeds: number[] = [];
  for (const part of parts) {
    const trimmed = part.trim();
    if (!trimmed) continue;
    const value = Number(trimmed);
    if (!Number.isFinite(value)) return [];
    if (value !== 0) seeds.push(value);
  }
  return seeds;
}

/**
 * Talks to the dashboard backend in dashboard/backend.
 *
 * The dev server proxies /api to it, see proxy.conf.json. Start the backend
 * first or these calls fail with a connection error, which the Tasks page
 * shows rather than swallowing.
 */
@Injectable({ providedIn: 'root' })
export class TaskService {
  private readonly http = inject(HttpClient);

  /** Not cached: the refresh button needs the next load to see new data. */
  load(): Observable<TaskData> {
    return this.http.get<TaskSnapshot>('/api/tasks').pipe(
      map((snapshot) => this.toTaskData(snapshot)),
      catchError((error) => throwError(() => new Error(describeError(error)))),
    );
  }

  /** The accessibility table, written by the same --fetch run as the snapshot. */
  loadCellTypes(): Observable<Record<string, CellType>> {
    return this.http
      .get<Record<string, CellType>>('/api/cell-types')
      .pipe(catchError((error) => throwError(() => new Error(describeError(error)))));
  }

  /**
   * How often each seed has been stamped, for the Tasks page's colouring.
   *
   * Failure is reported as an unavailable ledger rather than thrown: the
   * colouring is a layer over the table, and losing it must not take the
   * table's own load down with it.
   */
  loadSeedOccurrence(): Observable<SeedOccurrence> {
    return this.http
      .get<SeedOccurrence>('/api/seed-occurrence')
      .pipe(
        catchError((error) =>
          of({ available: false, counts: {}, reason: describeError(error) } as SeedOccurrence),
        ),
      );
  }

  /** Pull the upstream task history and merge it into the stored snapshot. */
  refresh(replace = false): Observable<RefreshResult> {
    return this.http
      .post<RefreshResult>('/api/tasks/refresh', { replace })
      .pipe(catchError((error) => throwError(() => new Error(describeError(error)))));
  }

  private toTaskData(snapshot: TaskSnapshot): TaskData {
    const rows = snapshot.tasks.map((task) => this.toRow(task));
    return {
      meta: {
        source: snapshot.source,
        fetchedAt: new Date(snapshot.fetched_at),
        count: snapshot.count ?? rows.length,
        unstamped: snapshot.unstamped ?? rows.filter((r) => r.seeds.length === 0).length,
        sharedRules: this.findSharedRules(snapshot.tasks),
      },
      rows,
    };
  }

  private toRow(task: RawTask): TaskRow {
    const contract = task.content.contract;
    const mutations: MutationEntry[] = contract.active_mutations.map((name) => ({
      name,
      region: contract.mutation_regions?.[name] ?? null,
      weight: contract.mutation_weights?.[name] ?? null,
    }));

    const seeds = parseSeeds(contract.seed);

    return {
      id: task.id,
      shortId: task.id.split('-')[0],
      createdAt: new Date(task.created_at),
      cellType: contract.cell_type,
      mutations,
      seeds,
      seedSort: seeds.length ? seeds[0] : null,
      seedRaw: seeds.length ? String(contract.seed) : null,
      version: contract.version,
      raw: task,
    };
  }

  /**
   * Every task in the snapshot carries the same rules, so they are worth
   * showing once above the table rather than as five identical columns.
   * Returns null if that ever stops being true, and the UI hides the panel.
   */
  private findSharedRules(tasks: RawTask[]): TaskRules | null {
    if (!tasks.length) return null;
    const first = tasks[0].content.contract.rules;
    const fingerprint = JSON.stringify(first);
    const uniform = tasks.every((t) => JSON.stringify(t.content.contract.rules) === fingerprint);
    return uniform ? first : null;
  }
}

/** A message worth showing the user, rather than 'Http failure response'. */
function describeError(error: unknown): string {
  if (!(error instanceof HttpErrorResponse)) {
    return error instanceof Error ? error.message : 'Unknown error';
  }
  if (error.status === 0) {
    return 'Cannot reach the dashboard backend. Start it with "uvicorn main:app --port 8000" in dashboard/backend.';
  }
  const detail = (error.error as { detail?: string } | null)?.detail;
  if (detail) return detail;
  return `Backend answered ${error.status} ${error.statusText}`;
}
