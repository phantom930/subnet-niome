import { Component, computed, inject, signal } from '@angular/core';
import { DecimalPipe } from '@angular/common';
import { AccordionModule } from 'primeng/accordion';
import { ButtonModule } from 'primeng/button';
import { CardModule } from 'primeng/card';
import { DividerModule } from 'primeng/divider';
import { MessageModule } from 'primeng/message';
import { TableModule } from 'primeng/table';
import { TagModule } from 'primeng/tag';
import { TimelineModule } from 'primeng/timeline';

import { TaskService } from '../../core/task.service';
import {
  BUILD_STEPS,
  CONTRACT_RULES,
  FACTORS,
  FIXED_CONSTRAINTS,
  HARNESS_COMMANDS,
  OCCUPANCY_SCENARIOS,
  OccupancyScenario,
  PITFALLS,
  REQUEST_STEPS,
  ROW_FIELDS,
} from './mining.data';

/** One drawn square of a panel's 2 x 4 grid. */
interface OccupancySquare {
  occupied: boolean;
  x: number;
  y: number;
}

/** One panel of the occupancy figure, with its factor derived from the cells. */
interface OccupancyPanel extends OccupancyScenario {
  /** distribution_fidelity_factor under this pattern. */
  fidelity: number;
  /** True when a ratio hit zero and stage 5's floor did the damage. */
  floored: boolean;
  /**
   * Grid geometry, precomputed because Angular template expressions have no
   * bitwise or floor operator and the row index needs one.
   */
  squares: OccupancySquare[];
  /** Where this panel sits along the figure's x axis. */
  offset: number;
}

/** One row of the derived accessibility table. */
interface AccessibilityRow {
  cellType: string;
  accessibility: number;
  basis?: string;
  /** stage 3's energy under the assumptions stated on the page. */
  energy: number;
  /** Whether energy saturates at 1.0 under those assumptions. */
  saturated: boolean;
  cas9CutProbability: number;
  cas12aCutProbability: number;
}

/**
 * The miner pins GC to 50 percent and picks the nearest usable PAM, so the
 * inner term of stage 3's energy is about 1.8 * 0.5 + 0.6 = 1.5 before
 * accessibility scales it. Stated on the page, because the real per-row value
 * moves with the actual GC and distance.
 */
const GC_AT_OPTIMUM = 0.5;
const DISTANCE_AT_OPTIMUM = 0;

/** stage5.geometric_mean's eps: every ratio is floored at this before the log. */
const FIDELITY_FLOOR = 1e-9;

/** How many ratios stage 5 takes the geometric mean of. */
const FIDELITY_RATIO_COUNT = 6;

/**
 * The two ratios in stage 5 that the occupancy pattern does not determine, held
 * at the values design.py measures for the current design so that the panels
 * differ only in which cells are occupied. Stated on the page.
 */
const KMER_DIVERSITY_RATIO = 0.985;
const DISTINCT_GUIDE_RATIO = 1;

/**
 * Occupancy figure geometry, in the SVG's own user units, which the figure's
 * max-width holds at 1:1 with CSS pixels.
 */
const PANEL_PITCH = 190;
const GRID_ORIGIN_X = 8;
const GRID_ORIGIN_Y = 22;
const SQUARE_PITCH = 37;
const SQUARE_COLUMNS = 4;

@Component({
  selector: 'app-mining',
  imports: [
    DecimalPipe,
    AccordionModule,
    ButtonModule,
    CardModule,
    DividerModule,
    MessageModule,
    TableModule,
    TagModule,
    TimelineModule,
  ],
  templateUrl: './mining.html',
  styleUrl: './mining.scss',
})
export class Mining {
  private readonly tasks = inject(TaskService);

  protected readonly requestSteps = REQUEST_STEPS;
  protected readonly buildSteps = BUILD_STEPS;
  protected readonly factors = FACTORS;
  protected readonly pitfalls = PITFALLS;
  protected readonly rowFields = ROW_FIELDS;
  protected readonly contractRules = CONTRACT_RULES;
  protected readonly fixedConstraints = FIXED_CONSTRAINTS;
  protected readonly harnessCommands = HARNESS_COMMANDS;

  protected readonly gcAtOptimum = GC_AT_OPTIMUM;

  protected readonly cellTypesError = signal<string | null>(null);
  private readonly cellTypes = signal<Record<
    string,
    { accessibility: number; basis?: string }
  > | null>(null);

  constructor() {
    this.tasks.loadCellTypes().subscribe({
      next: (table) => this.cellTypes.set(table),
      error: (error: Error) => this.cellTypesError.set(error.message),
    });
  }

  /**
   * stage3.sequence_energy, with the miner's own design substituted in.
   *
   *   energy = clamp(0, 1, accessibility * (1.8*gc + 0.6*exp(-dist/1500) + offset))
   *
   * The region offset is 0 for every mutation region observed so far, so it is
   * left out rather than guessed at.
   */
  private energyFor(accessibility: number): number {
    const inner = 1.8 * GC_AT_OPTIMUM + 0.6 * Math.exp(-DISTANCE_AT_OPTIMUM / 1500);
    return Math.max(0, Math.min(1, accessibility * inner));
  }

  /** stage3.cut_probability: clamp(0.4, 0.99, base + 0.18*energy). */
  private cutProbability(base: number, energy: number): number {
    return Math.min(0.99, Math.max(0.4, base + 0.18 * energy));
  }

  protected readonly accessibilityRows = computed<AccessibilityRow[]>(() => {
    const table = this.cellTypes();
    if (!table) return [];
    return Object.entries(table)
      .map(([cellType, entry]) => {
        const accessibility = entry.accessibility;
        const energy = this.energyFor(accessibility);
        return {
          cellType,
          accessibility,
          basis: entry.basis,
          energy,
          saturated: energy >= 1,
          cas9CutProbability: this.cutProbability(0.86, energy),
          cas12aCutProbability: this.cutProbability(0.78, energy),
        };
      })
      .sort((a, b) => a.accessibility - b.accessibility);
  });

  protected readonly loadingCellTypes = computed(
    () => this.cellTypes() === null && this.cellTypesError() === null,
  );

  /** The accessibility at which energy saturates under the stated assumptions. */
  protected readonly saturationPoint = computed(() => {
    const inner = 1.8 * GC_AT_OPTIMUM + 0.6 * Math.exp(-DISTANCE_AT_OPTIMUM / 1500);
    return 1 / inner;
  });

  protected readonly saturatedCount = computed(
    () => this.accessibilityRows().filter((row) => row.saturated).length,
  );

  /**
   * stage5.coverage_entropy_ratio: Shannon entropy of the observed counts over
   * the full declared support, divided by log2 of the support size. Uniform is
   * 1.0, everything in one bucket is 0.0, and a support of one is 1.0 by
   * convention.
   */
  private coverageEntropyRatio(counts: number[]): number {
    if (counts.length <= 1) return 1;
    const total = counts.reduce((sum, count) => sum + count, 0);
    if (total <= 0) return 0;

    let entropy = 0;
    for (const count of counts) {
      if (count <= 0) continue;
      const p = count / total;
      entropy -= p * Math.log2(p);
    }
    return entropy / Math.log2(counts.length);
  }

  /** stage5.geometric_mean, including the 1e-9 floor it applies first. */
  private geometricMean(values: number[]): number {
    const clipped = values.map((value) => Math.max(value, FIDELITY_FLOOR));
    const logSum = clipped.reduce((sum, value) => sum + Math.log(value), 0);
    return Math.exp(logSum / clipped.length);
  }

  /**
   * The six ratios for one occupancy pattern, with the rows spread evenly over
   * the occupied cells so that only occupancy varies between the panels.
   */
  private fidelityRatiosFor(cells: boolean[]): number[] {
    const rowsIn = (predicate: (index: number) => boolean) =>
      cells.filter((occupied, index) => occupied && predicate(index)).length;
    const bothOf = (axis: (index: number) => number) => [
      rowsIn((index) => axis(index) === 0),
      rowsIn((index) => axis(index) === 1),
    ];

    return [
      this.coverageEntropyRatio(bothOf((index) => index >> 2)),
      this.coverageEntropyRatio(bothOf((index) => (index >> 1) & 1)),
      this.coverageEntropyRatio(bothOf((index) => index & 1)),
      this.coverageEntropyRatio(cells.map((occupied) => (occupied ? 1 : 0))),
      KMER_DIVERSITY_RATIO,
      DISTINCT_GUIDE_RATIO,
    ];
  }

  protected readonly occupancyPanels: OccupancyPanel[] = OCCUPANCY_SCENARIOS.map(
    (scenario, panel) => {
      const ratios = this.fidelityRatiosFor(scenario.cells);
      return {
        ...scenario,
        fidelity: this.geometricMean(ratios),
        floored: ratios.some((ratio) => ratio <= 0),
        squares: scenario.cells.map((occupied, index) => ({
          occupied,
          x: GRID_ORIGIN_X + (index % SQUARE_COLUMNS) * SQUARE_PITCH,
          y: GRID_ORIGIN_Y + Math.floor(index / SQUARE_COLUMNS) * SQUARE_PITCH,
        })),
        offset: panel * PANEL_PITCH,
      };
    },
  );

  /**
   * The most a floored ratio can leave behind: 1e-9 through a sixth root is
   * 0.032, whatever the other five ratios do. This is the actual cliff, and it
   * is a collapsed dimension rather than an empty cell.
   */
  protected readonly flooredCeiling = Math.pow(FIDELITY_FLOOR, 1 / FIDELITY_RATIO_COUNT);
}
