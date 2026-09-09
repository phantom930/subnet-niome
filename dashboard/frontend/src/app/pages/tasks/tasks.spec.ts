import { TestBed } from '@angular/core/testing';
import { provideHttpClient } from '@angular/common/http';
import { HttpTestingController, provideHttpClientTesting } from '@angular/common/http/testing';
import { providePrimeNG } from 'primeng/config';

import { Tasks } from './tasks';
import { AppPreset } from '../../theme/app-preset';
import { RawTask, RefreshResult, TaskSnapshot } from '../../core/task.models';

const CELL_TYPES = ['HEK293', 'K562', 'CD34+_HSPC', 'HUDEP-2'];

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
        // Numeric 0 for unstamped, and stamped seeds arriving as both raw
        // numbers and comma-grouped strings, as the upstream really sends them.
        seed:
          i % 3 === 0 ? 0 : i % 3 === 1 ? 100000 + i : `${(100000 + i).toLocaleString('en-US')}`,
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

  const create = async () => {
    const fixture = TestBed.createComponent(Tasks);
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
      rows: () => Array<{ seed: number | null }>;
    };

    component.setSeedFilter('unstamped');
    await fixture.whenStable();

    const rows = component.rows();
    expect(rows.length).toBe(20);
    expect(rows.every((r) => r.seed === null)).toBe(true);
  });

  it('should normalize numeric and comma-grouped seeds to the same number', async () => {
    const fixture = await create();
    const component = fixture.componentInstance as unknown as {
      rows: () => Array<{ seed: number | null; seedRaw: string | null }>;
    };

    const byRaw = new Map(component.rows().map((r) => [r.seedRaw, r.seed]));
    expect(byRaw.get('100001')).toBe(100001);
    expect(byRaw.get('100,002')).toBe(100002);
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

    // A refresh is followed by a reload, so the table shows the new snapshot.
    http.expectOne('/api/tasks').flush(makeSnapshot());
    await fixture.whenStable();

    expect(bodyRows(fixture.nativeElement as HTMLElement).length).toBe(25);
  });

  it('should explain a backend that is not running', async () => {
    const fixture = TestBed.createComponent(Tasks);
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
