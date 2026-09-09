import { TestBed } from '@angular/core/testing';
import { provideHttpClient } from '@angular/common/http';
import { HttpTestingController, provideHttpClientTesting } from '@angular/common/http/testing';
import { providePrimeNG } from 'primeng/config';

import { Mining } from './mining';
import { AppPreset } from '../../theme/app-preset';
import { BUILD_STEPS, REQUEST_STEPS } from './mining.data';

/** The live table as the backend serves it. */
const CELL_TYPES = {
  'CD34+_HSPC': { accessibility: 0.87, basis: 'sourced' },
  HEK293: { accessibility: 0.35, basis: 'estimated' },
  'HUDEP-2': { accessibility: 0.82, basis: 'sourced' },
  K562: { accessibility: 0.77, basis: 'sourced' },
};

interface Internals {
  accessibilityRows: () => Array<{
    cellType: string;
    accessibility: number;
    energy: number;
    saturated: boolean;
    cas9CutProbability: number;
    cas12aCutProbability: number;
  }>;
  saturationPoint: () => number;
  saturatedCount: () => number;
  cellTypesError: () => string | null;
  occupancyPanels: Array<{
    label: string;
    fidelity: number;
    floored: boolean;
    squares: Array<{ occupied: boolean; x: number; y: number }>;
    offset: number;
  }>;
  flooredCeiling: number;
}

describe('Mining', () => {
  let http: HttpTestingController;

  const create = async (table: Record<string, { accessibility: number }> = CELL_TYPES) => {
    const fixture = TestBed.createComponent(Mining);
    http.expectOne('/api/cell-types').flush(table);
    await fixture.whenStable();
    return fixture;
  };

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

  it('should render every documented step', async () => {
    const fixture = await create();
    const text = (fixture.nativeElement as HTMLElement).textContent ?? '';

    for (const step of [...REQUEST_STEPS, ...BUILD_STEPS]) {
      expect(text).toContain(step.title);
    }
  });

  it('should name the source file for each step', async () => {
    const fixture = await create();
    const sources = Array.from(
      (fixture.nativeElement as HTMLElement).querySelectorAll('.step-source'),
    ).map((el) => el.textContent?.trim());

    expect(sources.length).toBe(REQUEST_STEPS.length + BUILD_STEPS.length);
    expect(sources).toContain('Miner._upload, _upload_headers');
  });

  it('should show the scoring formula with all three factors', async () => {
    const fixture = await create();
    const formula = (fixture.nativeElement as HTMLElement).querySelector('.formula');

    expect(formula?.textContent).toContain('total_weighted_score');
    expect(formula?.textContent).toContain('consistency_factor');
    expect(formula?.textContent).toContain('distribution_fidelity_factor');

    // Prettier puts the content on its own line inside <pre>. The parser drops
    // a newline immediately after the tag, so this must not render blank-first.
    expect(formula?.textContent?.startsWith('final_score')).toBe(true);
  });

  /**
   * The numbers below come from stage3.py: energy is accessibility times
   * (1.8*0.5 + 0.6) = 1.5, clamped to 1, and cut_probability is
   * clamp(0.4, 0.99, base + 0.18*energy).
   */
  it('should derive HEK293 energy and cut probabilities from the real formulas', async () => {
    const fixture = await create();
    const component = fixture.componentInstance as unknown as Internals;

    const hek = component.accessibilityRows().find((r) => r.cellType === 'HEK293');
    expect(hek).toBeDefined();
    // 0.35 * 1.5 = 0.525
    expect(hek!.energy).toBeCloseTo(0.525, 6);
    expect(hek!.saturated).toBe(false);
    // 0.86 + 0.18*0.525 = 0.9545, and 0.78 + 0.18*0.525 = 0.8745
    expect(hek!.cas9CutProbability).toBeCloseTo(0.9545, 6);
    expect(hek!.cas12aCutProbability).toBeCloseTo(0.8745, 6);
  });

  it('should saturate energy and cap Cas9 at 0.99 for a high-accessibility cell type', async () => {
    const fixture = await create();
    const component = fixture.componentInstance as unknown as Internals;

    const hspc = component.accessibilityRows().find((r) => r.cellType === 'CD34+_HSPC');
    // 0.87 * 1.5 = 1.305, clamped to 1.0
    expect(hspc!.energy).toBe(1);
    expect(hspc!.saturated).toBe(true);
    // 0.86 + 0.18 = 1.04, clamped to the 0.99 ceiling
    expect(hspc!.cas9CutProbability).toBe(0.99);
    expect(hspc!.cas12aCutProbability).toBeCloseTo(0.96, 6);
  });

  it('should report the saturation point and how many cell types reach it', async () => {
    const fixture = await create();
    const component = fixture.componentInstance as unknown as Internals;

    // 1 / 1.5
    expect(component.saturationPoint()).toBeCloseTo(0.6667, 4);
    // 0.87, 0.82 and 0.77 saturate; 0.35 does not.
    expect(component.saturatedCount()).toBe(3);
  });

  it('should order the accessibility table lowest first', async () => {
    const fixture = await create();
    const component = fixture.componentInstance as unknown as Internals;

    const values = component.accessibilityRows().map((r) => r.accessibility);
    expect(values).toEqual([...values].sort((a, b) => a - b));
  });

  /**
   * The three expected values are computed from stage5.py itself — entropy over
   * the declared support, then a geometric mean with each ratio floored at 1e-9
   * — with the k-mer ratio at 0.985 and the distinct-guide ratio at 1.0. They
   * are the point of the occupancy figure, so they are pinned here: an empty
   * cell is a haircut and a collapsed dimension is the cliff.
   */
  it('should derive the fidelity factor for each occupancy pattern', async () => {
    const fixture = await create();
    const panels = (fixture.componentInstance as unknown as Internals).occupancyPanels;

    expect(panels.map((p) => p.label)).toEqual([
      'all eight occupied',
      'one cell empty',
      'one Cas system only',
    ]);
    expect(panels[0].fidelity).toBeCloseTo(0.9975, 4);
    expect(panels[1].fidelity).toBeCloseTo(0.9792, 4);
    expect(panels[2].fidelity).toBeCloseTo(0.0295, 4);
  });

  it('should mark only the collapsed dimension as floored', async () => {
    const fixture = await create();
    const panels = (fixture.componentInstance as unknown as Internals).occupancyPanels;

    expect(panels.map((p) => p.floored)).toEqual([false, false, true]);
    // 1e-9 through a sixth root: the most a zeroed ratio can leave behind.
    expect((fixture.componentInstance as unknown as Internals).flooredCeiling).toBeCloseTo(
      0.0316,
      4,
    );
  });

  it('should lay every panel out as eight cells on one grid', async () => {
    const fixture = await create();
    const panels = (fixture.componentInstance as unknown as Internals).occupancyPanels;

    for (const panel of panels) {
      expect(panel.squares.length).toBe(8);
      // Two rows of four, so exactly two distinct y values and four x values.
      expect(new Set(panel.squares.map((s) => s.y)).size).toBe(2);
      expect(new Set(panel.squares.map((s) => s.x)).size).toBe(4);
    }

    expect(panels.map((p) => p.squares.filter((s) => !s.occupied).length)).toEqual([0, 1, 4]);
    expect(panels.map((p) => p.offset)).toEqual([0, 190, 380]);
  });

  /**
   * The figures tint their boxes from --p-indigo-400 and friends, which exist
   * only because the preset carries those primitive palettes. Swapping the
   * preset for one without them would leave the fills unresolved, so this fails
   * next to the diagrams rather than in the browser.
   */
  it('should expose every palette the diagrams tint with', () => {
    const primitive = (AppPreset as { primitive?: Record<string, Record<string, string>> })
      .primitive;

    for (const hue of ['indigo', 'amber', 'sky', 'emerald']) {
      expect(primitive?.[hue]?.['400']).toMatch(/^#[0-9a-f]{6}$/i);
    }
  });

  it('should give the roled boxes their hue and the rest the neutral default', async () => {
    const fixture = await create();
    const root = fixture.nativeElement as HTMLElement;

    // Every box in the two flow figures says which component it belongs to.
    const roled = root.querySelectorAll(
      'rect.dg-validator, rect.dg-store, rect.dg-miner, rect.dg-score',
    );
    // Five in the sequence figure, eleven in the pipeline.
    expect(roled.length).toBe(16);
    expect(root.querySelectorAll('line.dg-life[class*="dg-"]').length).toBe(3);
    // The two result boxes keep the accent outline on top of their role.
    expect(root.querySelectorAll('rect.dg-box-key.dg-score').length).toBe(2);
  });

  it('should give every diagram a text alternative', async () => {
    const fixture = await create();
    const figures = Array.from(
      (fixture.nativeElement as HTMLElement).querySelectorAll('figure.diagram'),
    );

    expect(figures.length).toBe(3);
    for (const figure of figures) {
      const svg = figure.querySelector('svg');
      expect(svg?.getAttribute('role')).toBe('img');
      // A label, not a placeholder: these carry the figure's whole claim.
      expect((svg?.getAttribute('aria-label') ?? '').length).toBeGreaterThan(80);
      expect(figure.querySelector('figcaption')?.textContent?.trim()).toBeTruthy();
    }
  });

  it('should render the occupancy grid as SVG rather than an image', async () => {
    const fixture = await create();
    const root = fixture.nativeElement as HTMLElement;

    // Three panels of eight cells, drawn as real shapes.
    expect(root.querySelectorAll('rect.dg-cell').length).toBe(24);
    expect(root.querySelectorAll('rect.dg-cell-empty').length).toBe(5);
    expect(root.querySelectorAll('img').length).toBe(0);
  });

  it('should warn rather than mislead when the cell-type table is empty', async () => {
    const fixture = await create({});
    const text = (fixture.nativeElement as HTMLElement).textContent ?? '';

    expect(text).toContain('cell-type table is empty');
    expect(text).toContain('defaults to 1.0');
  });

  it('should surface a failure to load the table', async () => {
    const fixture = TestBed.createComponent(Mining);
    http
      .expectOne('/api/cell-types')
      .error(new ProgressEvent('error'), { status: 0, statusText: 'Unknown Error' });
    await fixture.whenStable();

    const component = fixture.componentInstance as unknown as Internals;
    expect(component.cellTypesError()).toContain('Cannot reach the dashboard backend');
  });
});
