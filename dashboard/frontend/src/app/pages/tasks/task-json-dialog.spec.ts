import { TestBed } from '@angular/core/testing';
import { providePrimeNG } from 'primeng/config';

import { TaskJsonDialog } from './task-json-dialog';
import { AppPreset } from '../../theme/app-preset';
import { RawTask, TaskRow } from '../../core/task.models';

const TASK_ID = 'a155b953-bb0b-423d-8a3a-2189c571a3a2';

const RAW: RawTask = {
  id: TASK_ID,
  created_at: '2026-09-09T09:30:55.284775',
  content: {
    contract: {
      active_mutations: ['HBB:c.236T>C', 'NC_000011.10:g.5225715G>A'],
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
    hbb_reference: { window_id: 'w1', chromosome: 'chr11' },
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

describe('TaskJsonDialog', () => {
  const open = async (row: TaskRow | null = ROW) => {
    const fixture = TestBed.createComponent(TaskJsonDialog);
    fixture.componentRef.setInput('row', row);
    await fixture.whenStable();
    return fixture;
  };

  /** The dialog renders into an overlay, so query the document. */
  const jsonBlock = () => document.querySelector('.json');

  beforeEach(() => {
    TestBed.configureTestingModule({
      providers: [providePrimeNG({ theme: { preset: AppPreset } })],
    });
  });

  it('should render the whole task as valid, pretty-printed JSON', async () => {
    await open();
    const text = jsonBlock()?.textContent ?? '';

    // Round-trips, so what is shown is the task and not a lossy rendering.
    expect(JSON.parse(text)).toEqual(RAW);
    expect(text).toContain('\n  "id"');
  });

  it('should show nothing when there is no row', async () => {
    await open(null);
    expect(jsonBlock()).toBeNull();
  });

  it('should narrow to the contract when asked', async () => {
    const fixture = await open();
    const component = fixture.componentInstance as unknown as {
      setScope: (s: 'task' | 'contract') => void;
    };

    component.setScope('contract');
    await fixture.whenStable();

    const text = jsonBlock()?.textContent ?? '';
    expect(JSON.parse(text)).toEqual(RAW.content.contract);
    // hbb_reference belongs to the task, not the contract.
    expect(text).not.toContain('hbb_reference');
  });

  it('should preserve the seed exactly as sent, commas included', async () => {
    await open();
    const parsed = JSON.parse(jsonBlock()?.textContent ?? '{}');
    // The table normalizes the seed to a number; the JSON view must not.
    expect(parsed.content.contract.seed).toBe('546,343,346');
  });

  it('should report the line and byte count of what is shown', async () => {
    const fixture = await open();
    const component = fixture.componentInstance as unknown as {
      lineCount: () => number;
      byteCount: () => number;
    };

    const text = jsonBlock()?.textContent ?? '';
    expect(component.lineCount()).toBe(text.split('\n').length);
    expect(component.byteCount()).toBeGreaterThan(0);
  });

  it('should copy the shown JSON to the clipboard', async () => {
    const written: string[] = [];
    Object.defineProperty(navigator, 'clipboard', {
      value: { writeText: (text: string) => (written.push(text), Promise.resolve()) },
      configurable: true,
    });

    const fixture = await open();
    const component = fixture.componentInstance as unknown as {
      copy: () => Promise<void>;
      copied: () => boolean;
    };

    await component.copy();
    await fixture.whenStable();

    expect(written.length).toBe(1);
    expect(JSON.parse(written[0])).toEqual(RAW);
    expect(component.copied()).toBe(true);
  });

  it('should report a blocked clipboard instead of failing silently', async () => {
    Object.defineProperty(navigator, 'clipboard', {
      value: {
        writeText: () => Promise.reject(new Error('Document is not focused')),
      },
      configurable: true,
    });

    const fixture = await open();
    const component = fixture.componentInstance as unknown as {
      copy: () => Promise<void>;
      copied: () => boolean;
      copyFailed: () => boolean;
    };

    await component.copy();
    await fixture.whenStable();

    expect(component.copyFailed()).toBe(true);
    expect(component.copied()).toBe(false);
  });
});
