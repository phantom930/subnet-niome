import { TestBed } from '@angular/core/testing';
import { provideHttpClient } from '@angular/common/http';
import { HttpTestingController, provideHttpClientTesting } from '@angular/common/http/testing';
import { providePrimeNG } from 'primeng/config';

import { BenchmarkDialog } from './benchmark-dialog';
import { AppPreset } from '../../theme/app-preset';
import { BenchmarkJob, RawTask, TaskRow } from '../../core/task.models';

const TASK_ID = 'a155b953-bb0b-423d-8a3a-2189c571a3a2';

const RAW: RawTask = {
  id: TASK_ID,
  created_at: '2026-09-09T09:30:55.284775',
  content: {
    contract: {
      active_mutations: ['HBB:c.236T>C'],
      cell_type: 'K562',
      mutation_regions: { 'HBB:c.236T>C': 'exon2' },
      mutation_weights: { 'HBB:c.236T>C': 1.2 },
      rules: {
        base_padding: 400,
        cas_systems: ['Cas9', 'Cas12a'],
        max_experiments: 250,
        max_mismatches: 3,
        proximity_gate: false,
      },
      seed: '546,343,346',
      version: 'v1',
    },
    hbb_reference: { window_id: 'w1' },
  },
};

const ROW: TaskRow = {
  id: TASK_ID,
  shortId: 'a155b953',
  createdAt: new Date('2026-09-09T09:30:55Z'),
  cellType: 'K562',
  mutations: [{ name: 'HBB:c.236T>C', region: 'exon2', weight: 1.2 }],
  seed: 546343346,
  seedRaw: '546,343,346',
  version: 'v1',
  raw: RAW,
};

function job(status: BenchmarkJob['status'], extra: Partial<BenchmarkJob> = {}): BenchmarkJob {
  return {
    id: 'dabe0fd9c7c2',
    status,
    request: { task: ROW.id, seeds: 2, rng: null, task_seed: false, per_seed: true, uid: 0 },
    created_at: 0,
    started_at: 0,
    finished_at: null,
    duration_seconds: null,
    output: null,
    result: null,
    error: null,
    ...extra,
  };
}

/** Shaped like what the backend parses out of the harness's report. */
const DONE = job('done', {
  finished_at: 15,
  duration_seconds: 15.06,
  output: 'task a155b953...\n  final_score  26.3847',
  result: {
    task: { id: ROW.id, cell_type: 'K562', accessibility: 0.77 },
    miner: { rows: '250/250' },
    validator: {
      n_valid_experiments: 250,
      total_weighted_score: 267.724,
      consistency_score: 10.582,
      consistency_factor: 0.1058,
      distribution_fidelity_score: 0.9313,
      distribution_fidelity_factor: 0.9313,
      final_score: 26.3847,
    },
    seeds: [214, 754],
    per_seed: [
      {
        seed: 214,
        total_weighted_score: 267.724,
        consistency_factor: 0.1049,
        distribution_fidelity_factor: 0.9313,
        final_score: 26.1564,
      },
      {
        seed: 754,
        total_weighted_score: 267.724,
        consistency_factor: 0.1067,
        distribution_fidelity_factor: 0.9313,
        final_score: 26.613,
      },
    ],
    spread: 0.4566,
  },
});

describe('BenchmarkDialog', () => {
  let http: HttpTestingController;

  const open = async (row: TaskRow | null = ROW) => {
    const fixture = TestBed.createComponent(BenchmarkDialog);
    fixture.componentRef.setInput('row', row);
    await fixture.whenStable();
    return fixture;
  };

  /** The dialog renders into an overlay, so query the document. */
  const text = () => document.body.textContent ?? '';

  beforeEach(() => {
    TestBed.configureTestingModule({
      providers: [
        provideHttpClient(),
        provideHttpClientTesting(),
        providePrimeNG({ theme: { preset: AppPreset } }),
      ],
    });
    http = TestBed.inject(HttpTestingController);
  });

  afterEach(() => http.verify());

  it('should start a run as soon as a row is set', async () => {
    await open();

    const request = http.expectOne('/api/benchmarks');
    expect(request.request.method).toBe('POST');
    expect(request.request.body.task).toBe(ROW.id);
    request.flush(job('queued'));
  });

  it('should not start a run when there is no row', async () => {
    await open(null);
    // afterEach's verify() proves no request was made.
    expect(true).toBe(true);
  });

  it('should show the final score once the job is done', async () => {
    const fixture = await open();
    http.expectOne('/api/benchmarks').flush(DONE);
    await fixture.whenStable();

    expect(text()).toContain('26.3847');
    expect(text()).toContain('Final score');
  });

  it('should list a row per seed', async () => {
    const fixture = await open();
    http.expectOne('/api/benchmarks').flush(DONE);
    await fixture.whenStable();

    const seedRows = document.querySelectorAll('.p-datatable tbody tr');
    expect(seedRows.length).toBe(2);
    expect(text()).toContain('214');
    expect(text()).toContain('754');
  });

  it('should surface a failed run', async () => {
    const fixture = await open();
    http
      .expectOne('/api/benchmarks')
      .flush(job('failed', { error: 'the harness exited 1: task is unstamped' }));
    await fixture.whenStable();

    expect(text()).toContain('task is unstamped');
  });

  it('should explain a backend that is not running', async () => {
    const fixture = await open();
    http
      .expectOne('/api/benchmarks')
      .error(new ProgressEvent('error'), { status: 0, statusText: 'Unknown Error' });
    await fixture.whenStable();

    expect(text()).toContain('Cannot reach the dashboard backend');
  });

  it('should poll until the job settles', async () => {
    const fixture = await open();
    http.expectOne('/api/benchmarks').flush(job('running'));
    await fixture.whenStable();

    // Still running, so the progress state shows rather than a score.
    expect(text()).toContain('Running the miner');

    // The poll is on a timer; advance past it.
    await new Promise((resolve) => setTimeout(resolve, 1600));
    http.expectOne(`/api/benchmarks/${DONE.id}`).flush(DONE);
    await fixture.whenStable();

    expect(text()).toContain('26.3847');
  });
});
