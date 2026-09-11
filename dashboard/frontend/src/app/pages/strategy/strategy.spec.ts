import { TestBed } from '@angular/core/testing';
import { providePrimeNG } from 'primeng/config';

import { Strategy } from './strategy';
import { AppPreset } from '../../theme/app-preset';
import { BUILD_PHASES, GLOSSARY, LADDER, SCORE_DISTRIBUTION } from './strategy.data';

/**
 * The page derives its arithmetic rather than restating it, so the tests check
 * the derivations against the published formulas instead of against literals
 * copied out of the component.
 */
interface Internals {
  payoutBars: Array<{ rank: number; share: number; cumulative: number; height: number }>;
  topBlockShare: number;
  roundRows: Array<{
    hits: number;
    floor: number;
    consistency: number;
    final: number;
    places: boolean;
  }>;
  coverageCurves: Array<{ band: number; points: Array<{ hotkeys: number; probability: number }> }>;
  fleetCoverage: Array<{ band: number; atFleet: number; atLeaders: number }>;
  singleHotkey: Array<{ band: number; probability: number }>;
  inTtlBudget: number;
  glossaryTotal: number;
  glossaryMatches: () => number;
  glossaryGroups: () => Array<{ title: string; terms: Array<{ term: string }> }>;
  setGlossaryQuery: (value: string) => void;
}

/** 1 - (1 - union/900)^3, the page's own coverage formula. */
const spikeRate = (union: number) => 1 - Math.pow(1 - union / 900, 3);

describe('Strategy', () => {
  const create = async () => {
    const fixture = TestBed.createComponent(Strategy);
    await fixture.whenStable();
    return fixture;
  };

  const text = () => document.body.textContent ?? '';

  beforeEach(() => {
    TestBed.configureTestingModule({
      providers: [providePrimeNG({ theme: { preset: AppPreset } })],
    });
  });

  it('should render without a backend', async () => {
    const fixture = await create();
    expect(fixture.componentInstance).toBeTruthy();
    // No HTTP is involved: every number on the page is derived or recorded.
    expect(text()).toContain('The mining strategy');
  });

  it('should list every rung of the ladder in order', async () => {
    await create();
    for (const rung of LADDER) {
      expect(text()).toContain(rung.title);
    }
    expect(LADDER.map((rung) => rung.tag)).toEqual(['1', '2', '3', '4', '5']);
  });

  it('should list every build phase', async () => {
    await create();
    for (const phase of BUILD_PHASES) {
      expect(text()).toContain(phase.title);
    }
  });

  it('should carry the payout curve cumulatively', async () => {
    const fixture = await create();
    const bars = (fixture.componentInstance as unknown as Internals).payoutBars;

    expect(bars.length).toBe(SCORE_DISTRIBUTION.length);
    expect(bars[0].share).toBe(SCORE_DISTRIBUTION[0]);
    // The last bar's cumulative is the whole curve, which is what the top ten share.
    expect(bars[bars.length - 1].cumulative).toBeCloseTo(
      SCORE_DISTRIBUTION.reduce((a, b) => a + b, 0),
      10,
    );
    // Ranks 1 to 4, the block a lucky run of siblings collects.
    expect((fixture.componentInstance as unknown as Internals).topBlockShare).toBeCloseTo(0.85, 10);
  });

  it('should derive the round ladder from the per-seed floor', async () => {
    const fixture = await create();
    const rows = (fixture.componentInstance as unknown as Internals).roundRows;

    // Four rows for 0 to 3 hits, plus the all-cut-floor variant.
    expect(rows.length).toBe(5);
    for (const row of rows.slice(0, 4)) {
      expect(row.consistency).toBeCloseTo((row.hits * 1 + (3 - row.hits) * row.floor) / 3, 10);
    }
    // A round with no seed in the band never places; one seed does.
    expect(rows[0].places).toBe(false);
    expect(rows[1].places).toBe(true);
    // Consistency is monotone in hits, and a full sweep reaches exactly 1.0.
    expect(rows[3].consistency).toBeCloseTo(1, 10);
  });

  it('should compute coverage as one minus the miss probability cubed', async () => {
    const fixture = await create();
    const internals = fixture.componentInstance as unknown as Internals;

    for (const curve of internals.coverageCurves) {
      for (const point of curve.points) {
        expect(point.probability).toBeCloseTo(spikeRate(curve.band * point.hotkeys), 10);
      }
      // Coverage rises with every hotkey added — the claim the figure is making.
      const probabilities = curve.points.map((point) => point.probability);
      for (let i = 1; i < probabilities.length; i += 1) {
        expect(probabilities[i]).toBeGreaterThan(probabilities[i - 1]);
      }
    }

    for (const row of internals.fleetCoverage) {
      expect(row.atLeaders).toBeGreaterThan(row.atFleet);
    }
    for (const row of internals.singleHotkey) {
      expect(row.probability).toBeCloseTo(spikeRate(row.band), 10);
    }
  });

  it('should state the in-TTL budget as the deadline less the reserve', async () => {
    const fixture = await create();
    expect((fixture.componentInstance as unknown as Internals).inTtlBudget).toBe(225);
  });

  it('should render the glossary with every term', async () => {
    const fixture = await create();
    const internals = fixture.componentInstance as unknown as Internals;
    const terms = GLOSSARY.flatMap((group) => group.terms);

    expect(internals.glossaryTotal).toBe(terms.length);
    expect(internals.glossaryMatches()).toBe(terms.length);
    expect(document.querySelectorAll('.glossary-entry').length).toBe(terms.length);
    // A glossary that defines a word twice sends the reader to the wrong one.
    expect(new Set(terms.map((entry) => entry.term)).size).toBe(terms.length);
  });

  it('should filter the glossary by term and by definition', async () => {
    const fixture = await create();
    const internals = fixture.componentInstance as unknown as Internals;

    internals.setGlossaryQuery('min-union');
    const byTerm = internals.glossaryGroups().flatMap((group) => group.terms);
    expect(byTerm.some((entry) => entry.term === 'min-union')).toBe(true);
    expect(byTerm.length).toBeLessThan(internals.glossaryTotal);

    // The definitions are searched too: a reader usually has the wording from
    // the page rather than the headword.
    internals.setGlossaryQuery('geometric mean');
    expect(internals.glossaryGroups().flatMap((group) => group.terms)).toEqual([
      expect.objectContaining({ term: 'stage-5 cells' }),
    ]);

    internals.setGlossaryQuery('   ');
    expect(internals.glossaryMatches()).toBe(internals.glossaryTotal);
  });
});
