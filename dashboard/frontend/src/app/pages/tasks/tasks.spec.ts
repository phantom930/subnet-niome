import { TestBed } from '@angular/core/testing';
import { provideHttpClient } from '@angular/common/http';
import { HttpTestingController, provideHttpClientTesting } from '@angular/common/http/testing';
import { providePrimeNG } from 'primeng/config';

import { Tasks } from './tasks';
import { AppPreset } from '../../theme/app-preset';
import { RawTask, RefreshResult, SeedOccurrence, TaskSnapshot } from '../../core/task.models';

const CELL_TYPES = ['HEK293', 'K562', 'CD34+_HSPC', 'HUDEP-2'];

/** The seeds task `i` closed under: none, one, or three. */
function seedsFor(i: number): number[] {
  if (i % 3 === 0) return [];
  if (i % 3 === 1) return [100 + i];
  return [100 + i, 200 + i, 300 + i];
}

/** Those seeds in the shape the backend sends them. */
function seedFieldFor(i: number): string | number {
  const seeds = seedsFor(i);
  if (!seeds.length) return 0;
  return seeds.length === 1 ? seeds[0] : seeds.join(',');
}

/**
 * A ledger that gives seed 101 four stampings and 302 two, so the ramp's ends
 * are both exercised. Everything else counts once, which is the baseline.
 */
function makeOccurrence(): SeedOccurrence {
  const counts: Record<string, number> = {};
  for (let i = 0; i < TOTAL; i++) {
    for (const seed of seedsFor(i)) counts[String(seed)] = 1;
  }
  counts['101'] = 4;
  counts['302'] = 2;
  return {
    available: true,
    counts,
    // Before every task in the fixture, so none are excluded as too old.
    since: '2026-01-01T00:00:00',
    tasks_recorded: TOTAL,
    undrawn: 413,
  };
}

/** Mirrors the shape the backend serves, including its mixed-type seed. */
function makeTask(i: number): RawTask {
  const mutation = `HBB:c.${100 + i}A>G`;
  const other = `NC_000011.10:g.${5225000 + i}C>T`;
  return {
    id: `${i.toString(16).padStart(8, '0')}-8e88-455d-aa18-3a995626df20`,
    created_at: new Date(Date.UTC(2026, 8, 8, 12, 0, i)).toISOString(),
    content: {
      contract: {
        active_mutations: [mutation, other],
        cell_type: CELL_TYPES[i % CELL_TYPES.length],
        mutation_regions: { [mutation]: 'exon2', [other]: null },
        mutation_weights: { [mutation]: 1.2, [other]: 0.65 },
        rules: {
          base_padding: 400,
          cas_systems: ['Cas9', 'Cas12a'],
          max_experiments: 250,
          max_mismatches: 3,
          proximity_gate: false,
        },
        // As the upstream really sends it: 0 for a round that never closed, a
        // bare number for a single-seed round, and comma-joined seeds for the
        // three-seed rounds the backend has issued since 2026-08-27.
        seed: seedFieldFor(i),
        version: 'v1',
      },
      hbb_reference: {},
    },
  };
}

const TOTAL = 60;

function makeSnapshot(): TaskSnapshot {
  const tasks = Array.from({ length: TOTAL }, (_, i) => makeTask(i));
  return {
    source: 'https://example.invalid/api/v3/tasks',
    leaderboard: 'https://example.invalid/tasks',
    fetched_at: '2026-09-08T21:04:49+0000',
    order: 'newest first',
    note: 'test snapshot',
    count: tasks.length,
    unstamped: tasks.filter((t) => t.content.contract.seed === 0).length,
    tasks,
  };
}

describe('Tasks', () => {
  let http: HttpTestingController;

  /** `ledger` false serves an unavailable one, the degraded path. */
  const create = async (ledger: SeedOccurrence | false = makeOccurrence()) => {
    const fixture = TestBed.createComponent(Tasks);
    http
      .expectOne('/api/seed-occurrence')
      .flush(ledger === false ? { available: false, counts: {}, reason: 'no ledger' } : ledger);
    http.expectOne('/api/tasks').flush(makeSnapshot());
    await fixture.whenStable();
    return fixture;
  };

  const bodyRows = (el: HTMLElement) =>
    Array.from(el.querySelectorAll('tbody tr')).filter((r) => !r.querySelector('.empty-cell'));

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [Tasks],
      providers: [
        provideHttpClient(),
        provideHttpClientTesting(),
        providePrimeNG({ theme: { preset: AppPreset } }),
      ],
    });
    http = TestBed.inject(HttpTestingController);
  });

  afterEach(() => http.verify());

  it('should show only the first page of rows', async () => {
    const fixture = await create();
    const compiled = fixture.nativeElement as HTMLElement;

    // 60 tasks loaded, default page size 25.
    expect(bodyRows(compiled).length).toBe(25);
    expect(compiled.querySelector('.p-paginator')).toBeTruthy();
  });

  it('should report the page range over the full result set', async () => {
    const fixture = await create();
    const report = (fixture.nativeElement as HTMLElement).querySelector('.p-paginator-current');
    expect(report?.textContent?.trim()).toBe(`1 to 25 of ${TOTAL}`);
  });

  it('should advance to the next page', async () => {
    const fixture = await create();
    const compiled = fixture.nativeElement as HTMLElement;

    (compiled.querySelector('.p-paginator-next') as HTMLButtonElement).click();
    await fixture.whenStable();

    expect(compiled.querySelector('.p-paginator-current')?.textContent?.trim()).toBe(
      `26 to 50 of ${TOTAL}`,
    );
    expect(bodyRows(compiled).length).toBe(25);
  });

  it('should treat a numeric zero seed as unstamped', async () => {
    const fixture = await create();
    const component = fixture.componentInstance as unknown as {
      setSeedFilter: (v: 'unstamped') => void;
      rows: () => Array<{ seeds: number[] }>;
    };

    component.setSeedFilter('unstamped');
    await fixture.whenStable();

    const rows = component.rows();
    expect(rows.length).toBe(20);
    expect(rows.every((r) => r.seeds.length === 0)).toBe(true);
  });

  it('should split a comma-joined field into separate seeds', async () => {
    const fixture = await create();
    const component = fixture.componentInstance as unknown as {
      rows: () => Array<{ id: string; seeds: number[]; seedRaw: string | null }>;
    };

    const byRaw = new Map(component.rows().map((r) => [r.seedRaw, r.seeds]));
    // A three-seed round is three seeds, not the 9-digit number that stripping
    // the commas produced — task 2 closed under 102, 202 and 302.
    expect(byRaw.get('102,202,302')).toEqual([102, 202, 302]);
    // And a single-seed round is still one seed.
    expect(byRaw.get('101')).toEqual([101]);
  });

  it('should shade each seed by how often it has been stamped', async () => {
    const fixture = await create();
    const component = fixture.componentInstance as unknown as {
      rows: () => Array<{
        seeds: number[];
        chips: Array<{ seed: number; count: number | null; level: number }>;
      }>;
    };

    const chips = component.rows().flatMap((r) => r.chips);
    const bySeed = new Map(chips.map((c) => [c.seed, c]));

    expect([bySeed.get(101)?.count, bySeed.get(101)?.level]).toEqual([4, 4]);
    expect([bySeed.get(302)?.count, bySeed.get(302)?.level]).toEqual([2, 2]);
    expect([bySeed.get(202)?.count, bySeed.get(202)?.level]).toEqual([1, 1]);
  });

  it('should render one shaded chip per seed, with the count beside it', async () => {
    const fixture = await create();
    const compiled = fixture.nativeElement as HTMLElement;
    const component = fixture.componentInstance as unknown as { setSearch: (v: string) => void };

    // Task 2 closed under 102/202/302 but sorts onto a later page, so search
    // it up — which also exercises matching a single seed of a round.
    component.setSearch('302');
    await fixture.whenStable();

    const cells = Array.from(compiled.querySelectorAll('td.seed-cell'));
    expect(cells.length).toBe(1);
    const chips = Array.from(cells[0].querySelectorAll('.seed-chip'));

    expect(chips.length).toBe(3);
    expect(chips.map((c) => c.className.match(/occ-\d/)?.[0])).toEqual(['occ-1', 'occ-1', 'occ-2']);
    // The colour is never the only cue: a repeat also carries its count.
    expect(chips[2].querySelector('.occ-count')?.textContent?.trim()).toBe('×2');
    expect(chips[0].querySelector('.occ-count')).toBeNull();
    expect(chips[2].getAttribute('title')).toContain('stamped twice');

    // And the ramp has a key.
    expect(compiled.querySelectorAll('.seed-legend .legend-swatch').length).toBe(4);
  });

  it('should filter to rounds holding a repeated seed', async () => {
    const fixture = await create();
    const component = fixture.componentInstance as unknown as {
      setSeedFilter: (v: 'repeated') => void;
      rows: () => Array<{ maxOccurrence: number }>;
    };

    component.setSeedFilter('repeated');
    await fixture.whenStable();

    const rows = component.rows();
    // Task 1 stamped 101 (4x) and task 2 stamped 302 (2x); nothing else repeats.
    expect(rows.length).toBe(2);
    expect(rows.every((r) => r.maxOccurrence >= 2)).toBe(true);
  });

  it('should render seeds uncoloured when the ledger is unavailable', async () => {
    const fixture = await create(false);
    const component = fixture.componentInstance as unknown as {
      rows: () => Array<{ chips: Array<{ count: number | null; level: number }> }>;
    };

    const chips = component.rows().flatMap((r) => r.chips);
    expect(chips.length).toBeGreaterThan(0);
    expect(chips.every((c) => c.level === 0 && c.count === null)).toBe(true);
  });

  it('should not count seeds from rounds older than the ledger', async () => {
    const ledger = { ...makeOccurrence(), since: '2099-01-01T00:00:00' };
    const fixture = await create(ledger);
    const component = fixture.componentInstance as unknown as {
      rows: () => Array<{ chips: Array<{ count: number | null; level: number }> }>;
    };

    // Every fixture task predates that cutover, so none of their seeds carry
    // an occurrence — which must not read as "never stamped".
    const chips = component.rows().flatMap((r) => r.chips);
    expect(chips.every((c) => c.count === null && c.level === 0)).toBe(true);
  });

  it('should filter by cell type and reset to the first page', async () => {
    const fixture = await create();
    const compiled = fixture.nativeElement as HTMLElement;
    const component = fixture.componentInstance as unknown as {
      setCellType: (v: string) => void;
      first: () => number;
      rows: () => unknown[];
    };

    (compiled.querySelector('.p-paginator-next') as HTMLButtonElement).click();
    await fixture.whenStable();
    expect(component.first()).toBe(25);

    component.setCellType('K562');
    await fixture.whenStable();

    expect(component.first()).toBe(0);
    expect(component.rows().length).toBe(TOTAL / CELL_TYPES.length);
  });

  it('should surface the rules shared by every task', async () => {
    const fixture = await create();
    const note = (fixture.nativeElement as HTMLElement).querySelector('.rules-note');
    expect(note?.textContent).toContain('max 250 experiments');
    expect(note?.textContent).toContain('Cas9 and Cas12a');
  });

  it('should post a refresh and reload the snapshot afterwards', async () => {
    const fixture = await create();
    const component = fixture.componentInstance as unknown as { refresh: () => void };

    component.refresh();
    const refreshRequest = http.expectOne('/api/tasks/refresh');
    expect(refreshRequest.request.method).toBe('POST');
    expect(refreshRequest.request.body).toEqual({ replace: false });

    const result: RefreshResult = {
      added: 6,
      restamped: 1,
      removed: 0,
      count: 66,
      unstamped: 22,
      fetched_at: '2026-09-09T11:06:23+0000',
      cell_types: 4,
      replaced: false,
      output: 'testing/task.json: 66 tasks (+6 new, 1 newly stamped, 22 still unstamped)',
    };
    refreshRequest.flush(result);
    await fixture.whenStable();

    // A refresh is followed by a reload, so the table shows the new snapshot —
    // and the ledger is reloaded with it, because --fetch rewrites both.
    http.expectOne('/api/seed-occurrence').flush(makeOccurrence());
    http.expectOne('/api/tasks').flush(makeSnapshot());
    await fixture.whenStable();

    expect(bodyRows(fixture.nativeElement as HTMLElement).length).toBe(25);
  });

  it('should explain a backend that is not running', async () => {
    const fixture = TestBed.createComponent(Tasks);
    // Both calls fail when the backend is down. The ledger's failure is
    // swallowed into an unavailable ledger; the snapshot's is what shows.
    http
      .expectOne('/api/seed-occurrence')
      .error(new ProgressEvent('error'), { status: 0, statusText: 'Unknown Error' });
    http
      .expectOne('/api/tasks')
      .error(new ProgressEvent('error'), { status: 0, statusText: 'Unknown Error' });
    await fixture.whenStable();

    const compiled = fixture.nativeElement as HTMLElement;
    expect(compiled.querySelector('.load-error')?.textContent).toContain(
      'Cannot reach the dashboard backend',
    );
    // The table is replaced by the error, not shown empty.
    expect(compiled.querySelector('.p-datatable')).toBeNull();
  });
});
