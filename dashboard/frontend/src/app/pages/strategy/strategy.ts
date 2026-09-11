import { Component } from '@angular/core';
import { DecimalPipe, PercentPipe } from '@angular/common';
import { AccordionModule } from 'primeng/accordion';
import { CardModule } from 'primeng/card';
import { MessageModule } from 'primeng/message';
import { TableModule } from 'primeng/table';
import { TagModule } from 'primeng/tag';
import { TimelineModule } from 'primeng/timeline';

import {
  BAND_FACTS,
  BUDGET_GATES,
  BUILD_PHASES,
  CONSTRUCTION_COMPARISON,
  CUTOFF,
  DRIFT,
  ENVELOPE,
  FALSIFIED,
  FILES,
  LADDER,
  OPEN,
  SCORE_DISTRIBUTION,
  SEED_MODEL,
  SEED_REGIMES,
  SETTLED,
  SHARP_EDGES,
  SOURCE_REVISION,
  WINDOW_LAYOUTS,
} from './strategy.data';

/** One bar of the payout figure. */
interface PayoutBar {
  rank: number;
  share: number;
  /** Cumulative share through this rank, which is what a block of siblings collects. */
  cumulative: number;
  x: number;
  y: number;
  height: number;
}

/** One row of the derived round-score ladder. */
interface RoundRow {
  label: string;
  /** Seeds of the round's three that hit the construction's clean set. */
  hits: number;
  /** The per-seed factor on a hit. */
  spike: number;
  /** The per-seed factor on a miss. */
  floor: number;
  consistency: number;
  final: number;
  places: boolean;
}

/** One point of the fleet coverage curve. */
interface CoveragePoint {
  hotkeys: number;
  probability: number;
  x: number;
  y: number;
}

/** One plotted coverage curve, for a band size. */
interface CoverageCurve {
  label: string;
  band: number;
  points: CoveragePoint[];
  path: string;
}

/**
 * A round draws three seeds, independently and uniformly, from 100-999. Every
 * frequency number on this page comes out of that.
 */
const SEEDS_PER_ROUND = 3;
const SEED_SPACE = 900;

/**
 * The stated design point for the round-score ladder, so the four rows differ
 * only in how many seeds hit. Both are this branch's own measured values for a
 * shipped build, and they are stated on the page rather than derived, because
 * they move with the contract.
 */
const STATED_WEIGHTED = 230;
const STATED_FIDELITY = 0.95;

/** Hotkeys plotted on the coverage curve, and the two counts worth marking. */
const MAX_HOTKEYS = 14;
const FLEET_HOTKEYS = 10;
const LEADER_HOTKEYS = 14;

/** Payout figure geometry, in the SVG's own user units. */
const PAYOUT_WIDTH = 560;
const PAYOUT_HEIGHT = 150;
const PAYOUT_BASELINE = 120;
const PAYOUT_BAR_PITCH = 44;
const PAYOUT_BAR_WIDTH = 30;
const PAYOUT_ORIGIN_X = 10;
/** The tallest bar, rank 1's 30%, drawn at this height. */
const PAYOUT_SCALE = 300;

/** Coverage figure geometry, in the SVG's own user units. */
const COVERAGE_WIDTH = 560;
const COVERAGE_HEIGHT = 220;
const COVERAGE_LEFT = 44;
const COVERAGE_RIGHT = 540;
const COVERAGE_TOP = 16;
const COVERAGE_BOTTOM = 178;
/** The y axis runs 0 to this, so the curves fill the panel. */
const COVERAGE_Y_MAX = 0.6;

@Component({
  selector: 'app-strategy',
  imports: [
    DecimalPipe,
    PercentPipe,
    AccordionModule,
    CardModule,
    MessageModule,
    TableModule,
    TagModule,
    TimelineModule,
  ],
  templateUrl: './strategy.html',
  styleUrl: './strategy.scss',
})
export class Strategy {
  protected readonly revision = SOURCE_REVISION;
  protected readonly ladder = LADDER;
  protected readonly buildPhases = BUILD_PHASES;
  protected readonly comparison = CONSTRUCTION_COMPARISON;
  protected readonly seedRegimes = SEED_REGIMES;
  protected readonly windowLayouts = WINDOW_LAYOUTS;
  protected readonly bandFacts = BAND_FACTS;
  protected readonly falsified = FALSIFIED;
  protected readonly drift = DRIFT;
  protected readonly cutoff = CUTOFF;
  protected readonly envelope = ENVELOPE;
  protected readonly budgetGates = BUDGET_GATES;
  protected readonly sharpEdges = SHARP_EDGES;
  protected readonly seedModel = SEED_MODEL;
  protected readonly settled = SETTLED;
  protected readonly open = OPEN;
  protected readonly files = FILES;

  /**
   * The prediction against the two things it has to beat, ordered so the model
   * sits between them and the ordering makes the verdict readable without the
   * p-value.
   */
  protected readonly seedModelRows = [
    {
      label: `best reference (${SEED_MODEL.bestReferenceName})`,
      hits: SEED_MODEL.bestReferenceHits,
    },
    { label: 'the model', hits: SEED_MODEL.modelHits },
    { label: 'chance, against each round’s own draw', hits: SEED_MODEL.chanceHits },
  ];

  /**
   * The budget a build racing the upload TTL actually gets: the deadline less
   * the reserve held back for the scan's tail and the PUT. Stated as a round
   * number on the page because the artifact fetches come out of it too.
   */
  protected readonly inTtlBudget = ENVELOPE.ttlSeconds - ENVELOPE.hedgeReserve;

  protected readonly seedsPerRound = SEEDS_PER_ROUND;
  protected readonly seedSpace = SEED_SPACE;
  protected readonly statedWeighted = STATED_WEIGHTED;
  protected readonly statedFidelity = STATED_FIDELITY;
  protected readonly fleetHotkeys = FLEET_HOTKEYS;
  protected readonly leaderHotkeys = LEADER_HOTKEYS;

  protected readonly payoutWidth = PAYOUT_WIDTH;
  protected readonly payoutHeight = PAYOUT_HEIGHT;
  protected readonly payoutBaseline = PAYOUT_BASELINE;
  protected readonly payoutBarWidth = PAYOUT_BAR_WIDTH;
  protected readonly coverageWidth = COVERAGE_WIDTH;
  protected readonly coverageHeight = COVERAGE_HEIGHT;
  protected readonly coverageLeft = COVERAGE_LEFT;
  protected readonly coverageRight = COVERAGE_RIGHT;
  protected readonly coverageTop = COVERAGE_TOP;
  protected readonly coverageBottom = COVERAGE_BOTTOM;

  /** The floor every submission scores on a seed it does not pin. */
  private readonly floor = SEED_REGIMES.find((regime) => regime.pinned === 'none')?.factor ?? 0.101;

  /** all-cut's per-seed value where its own rows hold is_cut constant. */
  private readonly allCutFloor =
    SEED_REGIMES.find((regime) => regime.label.includes('all-cut'))?.factor ?? 0.237;

  /**
   * The payout curve as bars, with the cumulative share carried alongside.
   *
   * The cumulative column is the one that matters for a fleet: correlated
   * siblings score within noise of each other and so take consecutive ranks,
   * collecting a run of this curve rather than n times its top.
   */
  protected readonly payoutBars: PayoutBar[] = SCORE_DISTRIBUTION.map((share, index) => {
    const cumulative = SCORE_DISTRIBUTION.slice(0, index + 1).reduce((sum, s) => sum + s, 0);
    const height = share * PAYOUT_SCALE;
    return {
      rank: index + 1,
      share,
      cumulative,
      x: PAYOUT_ORIGIN_X + index * PAYOUT_BAR_PITCH,
      y: PAYOUT_BASELINE - height,
      height,
    };
  });

  /** What ranks 1 to 4 collect together, which is what a lucky block takes. */
  protected readonly topBlockShare = SCORE_DISTRIBUTION.slice(0, 4).reduce((a, b) => a + b, 0);

  /** What ranks 6 to 9 collect together, which is what the same block takes one rank later. */
  protected readonly nextBlockShare = SCORE_DISTRIBUTION.slice(5, 9).reduce((a, b) => a + b, 0);

  /**
   * The round score for k of three seeds hitting, at the stated design point.
   *
   * `consistency_factor` is averaged over the round's seeds, so with a floor
   * near 0.10 and a spike at 1.0 the reachable values are almost discrete —
   * which is why the question is how often a seed is hit, not how good the
   * build is. The last row is the same arithmetic with all-cut's own clean-seed
   * value as the floor, and it is where the field's top block actually lives.
   */
  protected readonly roundRows: RoundRow[] = [
    ...[0, 1, 2, 3].map((hits) => ({
      label:
        hits === 0
          ? 'no seed in the band'
          : hits === 1
            ? 'one seed in the band'
            : `${hits} seeds in the band`,
      hits,
      spike: 1,
      floor: this.floor,
      ...this.roundScore(hits, 1, this.floor),
    })),
    {
      label: 'one band seed, on an all-cut floor',
      hits: 1,
      spike: 1,
      floor: this.allCutFloor,
      ...this.roundScore(1, 1, this.allCutFloor),
    },
  ];

  /** The single-hotkey spike rate at each cell type's band size. */
  protected readonly singleHotkey = [
    { label: 'erythroid cell types, band 13', band: 13, probability: this.spikeRate(13) },
    { label: 'HEK293, band 7', band: 7, probability: this.spikeRate(7) },
  ];

  /**
   * Coverage against hotkey count.
   *
   * The union of the fleet's bands is taken as `band × hotkeys`, which assumes
   * the bands are disjoint. They are not quite: a band leaks outside the window
   * it was searched in, so sibling bands collide more than a disjoint plan
   * implies. Read these curves as upper bounds — the shape is the point, and
   * the shape is why hotkey count is the lever with headroom left.
   */
  protected readonly coverageCurves: CoverageCurve[] = [
    { label: 'band 13 (erythroid)', band: 13 },
    { label: 'band 7 (HEK293)', band: 7 },
  ].map(({ label, band }) => {
    const points: CoveragePoint[] = [];
    for (let hotkeys = 1; hotkeys <= MAX_HOTKEYS; hotkeys += 1) {
      const probability = this.spikeRate(band * hotkeys);
      points.push({
        hotkeys,
        probability,
        x: COVERAGE_LEFT + ((hotkeys - 1) / (MAX_HOTKEYS - 1)) * (COVERAGE_RIGHT - COVERAGE_LEFT),
        y:
          COVERAGE_BOTTOM -
          Math.min(1, probability / COVERAGE_Y_MAX) * (COVERAGE_BOTTOM - COVERAGE_TOP),
      });
    }
    return {
      label,
      band,
      points,
      path: points
        .map((p, i) => `${i === 0 ? 'M' : 'L'}${p.x.toFixed(1)} ${p.y.toFixed(1)}`)
        .join(' '),
    };
  });

  /** Gridlines for the coverage figure, at every 20 percent. */
  protected readonly coverageGrid = [0, 0.2, 0.4, 0.6].map((value) => ({
    value,
    y: COVERAGE_BOTTOM - (value / COVERAGE_Y_MAX) * (COVERAGE_BOTTOM - COVERAGE_TOP),
  }));

  /** x positions for the hotkey axis labels. */
  protected readonly coverageTicks = [1, 4, 7, 10, 14].map((hotkeys) => ({
    hotkeys,
    x: COVERAGE_LEFT + ((hotkeys - 1) / (MAX_HOTKEYS - 1)) * (COVERAGE_RIGHT - COVERAGE_LEFT),
  }));

  /** What the fleet's own count buys, at each band size. */
  protected readonly fleetCoverage = this.coverageCurves.map((curve) => ({
    label: curve.label,
    band: curve.band,
    atFleet: this.spikeRate(curve.band * FLEET_HOTKEYS),
    atLeaders: this.spikeRate(curve.band * LEADER_HOTKEYS),
  }));

  /**
   * `1 - (1 - union/900)³` — the chance that at least one of a round's three
   * independent seeds lands in a set of `union` seeds.
   */
  private spikeRate(union: number): number {
    const share = Math.min(1, Math.max(0, union / SEED_SPACE));
    return 1 - Math.pow(1 - share, SEEDS_PER_ROUND);
  }

  /**
   * The round's `consistency_factor` and final score for `hits` of its three
   * seeds pinned, and whether it clears the median rank-10 cutoff.
   */
  private roundScore(
    hits: number,
    spike: number,
    floor: number,
  ): { consistency: number; final: number; places: boolean } {
    const consistency = (hits * spike + (SEEDS_PER_ROUND - hits) * floor) / SEEDS_PER_ROUND;
    const final = STATED_WEIGHTED * consistency * STATED_FIDELITY;
    return { consistency, final, places: final >= CUTOFF.median };
  }
}
