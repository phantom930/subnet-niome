/**
 * Content for the Strategy page — the mining approach on the `develop` branch.
 *
 * The Mining page documents one miner answering one validator and maximising
 * one submission. This page documents what `develop` does with the same
 * pipeline: a portfolio of constructions played across a fleet of hotkeys,
 * betting on a rare consistency spike rather than on a good average.
 *
 * Sourced from the code on `origin/develop` at 527588a (2026-09-09), read
 * directly rather than from its notes:
 *   neurons/miner.py                          the ladder, the budgets, the prefetch
 *   niome_subnet/genomics/all_cut.py           the whole-window cut pin
 *   niome_subnet/genomics/all_hdr.py           the narrow-band HDR pin
 *   niome_subnet/genomics/seed_depend.py       the seed-0 build
 *   niome_subnet/genomics/seed_agnostic.py     the cut-only hedge
 *   niome_subnet/genomics/fastgreedy.py        the min-union selector
 *   niome_subnet/genomics/mt19937.py           the bit-exact RNG replica
 *   window_plan.py, miner.sh, round_plan.sh    the fleet and its window plan
 *   niome_subnet/utils/settings.py             the payout curve
 *   CLAUDE.md                                  the measurement record
 *
 * Kept in one file so the numbers have a single place to be updated, and every
 * item names the file it came from. Where the branch's own CLAUDE.md disagrees
 * with the code, the code wins and the disagreement is listed in DRIFT below.
 */

/** The branch state this page describes. */
export const SOURCE_REVISION = {
  branch: 'develop',
  commit: '527588a',
  subject: 'seed model updated',
  date: '2026-09-09',
};

/**
 * settings.SCORE_DISTRIBUTION with SCORING_SYSTEM = "top": the fraction of the
 * round's emission each rank is paid. Rank 11 and below are paid nothing, which
 * is the fact the whole approach is organised around.
 */
export const SCORE_DISTRIBUTION = [0.3, 0.2, 0.2, 0.15, 0.05, 0.03, 0.025, 0.02, 0.015, 0.01];

/**
 * Per-seed `consistency_factor` by how many of stage 4's three targets the rows
 * hold exactly constant.
 *
 * This is the pivot of the whole approach, and it is a property of the scoring
 * code rather than of biology: `r2_score` returns exactly 1.0 for a constant
 * `y_test` predicted exactly, and `normalized_mae` returns the raw MAE
 * unnormalised when `std(y) < 1e-9`. So what a seed pays is set by how many
 * targets are pinned, not by how learnable they are.
 *
 * Measured on HUDEP-2 (9ac1d178 for all-HDR rows, 4725e952 for all-cut rows),
 * recorded in CLAUDE.md. The two `is_cut` rows differ only in row composition,
 * which is why the same pin is worth more under all-cut's rows: all-cut's rows
 * overfit less on the two targets that are still free.
 */
export interface SeedRegime {
  label: string;
  /** Which of is_cut / is_hdr / indel_length are exactly constant. */
  pinned: string;
  factor: number;
  /** Range across the sampled seeds, where more than one was sampled. */
  range?: string;
  samples?: number;
}

export const SEED_REGIMES: SeedRegime[] = [
  {
    label: 'every row repairs by HDR',
    pinned: 'all three',
    factor: 1.0,
    range: 'exact',
    samples: 2,
  },
  {
    label: 'every row cuts, all-cut rows',
    pinned: 'is_cut',
    factor: 0.237,
    range: '0.139 – 0.303',
    samples: 14,
  },
  {
    label: 'every row cuts, all-HDR rows',
    pinned: 'is_cut',
    factor: 0.162,
    range: '0.135 – 0.203',
    samples: 5,
  },
  {
    label: 'neither holds',
    pinned: 'none',
    factor: 0.101,
    range: '0.099 – 0.121',
    samples: 10,
  },
];

/**
 * The measured rank-10 cutoff, per cell type, over current-regime rounds. A
 * round final at or above the cutoff is paid; below it is worth nothing.
 *
 * Recorded in CLAUDE.md over 72 rounds. Quoted as a distribution on purpose:
 * an earlier draft there recorded a single "62-86 band" and drew the wrong
 * conclusion from it, because whether a construction places is decided by how
 * soft that round's field happened to be.
 */
export const CUTOFF = {
  min: 20.7,
  p10: 45.9,
  median: 82.3,
  p90: 120.0,
  max: 158.1,
  rounds: 72,
  byCell: [
    { cellType: 'HEK293', median: 76.4 },
    { cellType: 'CD34+_HSPC', median: 81.4 },
    { cellType: 'HUDEP-2', median: 83.8 },
    { cellType: 'K562', median: 90.9 },
  ],
};

/**
 * A rung of `Miner._build`'s ladder, in the order the code tries them. Each
 * falls through to the next on a decline, a short pool, a missing GPU or too
 * little budget, so a failure anywhere lands on the build the miner had before
 * that rung existed.
 *
 * Every builder signals a decline the same way — `(None, meta)` with a
 * human-readable `meta["reason"]`, never an exception and never an empty list —
 * which is what makes the ladder safe to fall down.
 */
export interface Rung {
  tag: string;
  title: string;
  source: string;
  /** What the rows hold constant, in the validator's terms. */
  pins: string;
  /** The gate in `_build`, as the code tests it. */
  gate: string;
  summary: string;
  details: string[];
  icon: string;
}

export const LADDER: Rung[] = [
  {
    tag: '1',
    title: 'Seed-depend',
    source: 'genomics/seed_depend.py',
    pins: 'all three targets, at seed 0 only',
    gate: 'NIOME_SEED_DEPEND set and budget ≥ 360 s',
    summary:
      'For the rounds the backend never stamps. The validator then scores at seed 0 — the value the miner was handed — so the seed is known and consistency reaches exactly 1.000.',
    details: [
      'It replaces the construction rather than hedging alongside it: the rows are pinned to seed 0 and score the ~0.10 floor on any round that is stamped. That is why it is first, and why listing a hotkey here turns its band off.',
      'The objective has no frequency term at all. It maximises sum(weighted_score) × fidelity subject to every row satisfying the rule at the target seed, and the two pull against each other — so allocate searches the mutation split and scores each candidate with stage 5 itself rather than a proxy.',
      'Never-stamped rounds are 8 of 105 since 2026-08-25, 7.6 percent. Do not read that rate off the task listing, which reports seed 0 for rounds that were in fact stamped; the reliable test is miners reaching consistency 1.000 in the score rows.',
      'variants_per_site 12000 is the one knob that paid: +1.47 ± 0.25, t = 5.80, on 6 of 6 rounds. It moved mean payout share 25.5 to 25.8 percent, which is to say nothing — the points crossed no rank boundary. 24000 adds +0.12 for twice the build and twice the memory.',
      'The build is deterministic, so several hotkeys would submit byte-identical files. variant and variant_epsilon break ties among near-equal candidates by hashing (guide, variant), changing which guide is taken without changing what is optimised.',
    ],
    icon: 'pi pi-key',
  },
  {
    tag: '2',
    title: 'All-HDR',
    source: 'genomics/all_hdr.py',
    pins: 'all three targets, on a ~13-seed band',
    gate: 'ALL_HDR and cell type measured and budget ≥ 190 s',
    summary:
      'The spike. Min-union a Cas12a group on the HDR rule over a 100-seed band, then require strict HDR of the Cas9 half over the resulting clean band. On a clean-band seed every row repairs by HDR, so all three of stage 4’s targets are constant and consistency_factor is exactly 1.0.',
    details: [
      'It loses on the mean and is shipped anyway: 51 to 61 percent of all-cut’s expected score, about −24 points a round. The justification is the payout curve — a build worth ~290 on a few percent of rounds can out-rank a flat one that rarely places.',
      'The clean band is 15 to 16 seeds of the 100 screened on the erythroid cell types and 6 to 10 on HEK293, so roughly 1.8 percent of the seed space. Verified against stage3.simulate itself, not only against the GPU replica.',
      'group_size 80 is a payout optimum that deliberately overrides an earlier score optimum. It trades band width — spike frequency — against stage 5’s cas-coverage entropy, which is spike score, and priced over 54 current-regime fields it is +14.2 percent on expected payout over group 100.',
      'main_max_fail 45 of 100 band seeds is deliberately loose: the min-union step is what produces the clean band, and a tighter screen shrinks the pool it selects from without improving the group. HEK293 uses 48, where the bank holds only 0.041 percent of guides.',
      'A band override rescales main_max_fail by z-score rather than linearly, because the median guide fails 51 of 100 seeds but 153 of 300. Linear scaling put HEK293 at z −5.66, where 53 of 60000 bank guides qualified against a group of 80, and every band hotkey declined.',
    ],
    icon: 'pi pi-bolt',
  },
  {
    tag: '3',
    title: 'All-cut',
    source: 'genomics/all_cut.py',
    pins: 'is_cut, across the whole 900-seed window',
    gate: 'ALL_CUT and cell type measured and budget ≥ 190 s (HEK293) or 480 s',
    summary:
      'The same machine with one string changed: the screening rule is "cut", which admits every repair mode. Only is_cut goes constant, but it does so on 12.7 to 29.2 percent of all 900 seeds rather than on 1.8 percent.',
    details: [
      'It was moved ahead of the seed-agnostic hedge because it is that hedge with the Cas9 constraint relaxed from strict-over-900-seeds to strict-over-the-Cas12a-clean-set, and it measured +23.45 (t = 8.4) over it on K562. Placed after, it was dead code.',
      'group_size 42 gives an 83/17 cas mix, the peak of the clean-payout sweep.',
      'The budget gate is per cell type — 190 s for HEK293, 480 s for the three erythroid types — because a single flat gate would either start a build that cannot finish or lock HEK293 out of a path it completes. That makes erythroid all-cut prefetch-dependent.',
      'Its measured per-seed value is 0.20 to 0.23, not 1.0: pinning is_cut alone is worth far less than stage 4’s 0.7 R² weight suggests, because is_hdr and indel_length overfit harder once is_cut goes constant.',
      'On K562 the matched pricing — both arms built on the same contract, each placed against the one field that played it — puts all-cut ahead 3.78x, reversing the pooled estimate. HEK293 has never had that treatment and is the open case.',
    ],
    icon: 'pi pi-scissors',
  },
  {
    tag: '4',
    title: 'Seed-agnostic hedge',
    source: 'genomics/seed_agnostic.py',
    pins: 'is_cut, across the whole 900-seed window',
    gate: 'SEED_AGNOSTIC and not HEK293 and budget ≥ 150 s',
    summary:
      'The older cut-only hedge: a strict Cas9 bank that cuts under every seed in the window, plus a min-union Cas12a group, assembled at an 85/15 mix.',
    details: [
      'It works because cut_p clamps at 0.99 for Cas9 on a high-accessibility cell type, so 0.99^900 ≈ 1.2e-4 still leaves a few hundred strict guides. Cas12a caps at 0.96, where 0.96^900 ≈ 1.6e-16 and no strict guide exists at that width — hence a group whose failures coincide.',
      'Skipped for HEK293, where accessibility 0.35 holds Cas9 to about 0.95 and no strict guide exists at any width.',
      'Its ceiling bounds every hedge of this shape: only is_cut is recoverable across a window. Repair mode is a fresh ~0.5 coin per seed with no design lever, so every repair-mode rule tested tops out at 349 to 491 failed seeds of 900 against the 12 to 22 that make the cut-only hedge work.',
      'The 210 s time budget stops the scan even with its pool cap unmet, so a build always returns inside the URL’s TTL rather than overrunning it and losing the round.',
    ],
    icon: 'pi pi-shield',
  },
  {
    tag: '5',
    title: 'Ordinary construction',
    source: 'genExp.py, genomics/hek293_generation.py',
    pins: 'nothing — it scores the floor',
    gate: 'always available',
    summary:
      'The deterministic build the miner had before any of the rungs above existed. HEK293 uses its own clustered builder instead.',
    details: [
      'This is what allow_hedges=False jumps straight to. That path exists for the case where a prepared build is still holding the GPU and only the ordinary construction will finish in time.',
      'It is not hypothetical: two hotkeys were called 129 s and 149 s after task creation on one round and shipped byte-identical emergency builds, scoring 26.15. Any hotkey can be asked cold, so the last rung has to always finish.',
    ],
    icon: 'pi pi-file',
  },
];

/**
 * The two spike constructions side by side. Both are the same machine with one
 * rule string changed, which is the fact worth carrying away: the difference in
 * outcome is entirely in how many targets the rule pins and over how many
 * seeds it can hold them.
 */
export interface ConstructionRow {
  property: string;
  allCut: string;
  allHdr: string;
}

export const CONSTRUCTION_COMPARISON: ConstructionRow[] = [
  { property: 'screening rule', allCut: '"cut" — any repair mode', allHdr: '"hdr" — HDR only' },
  { property: 'seed axis', allCut: '900 seeds, 100–999', allHdr: '100-seed band, per cell type' },
  {
    property: 'stage-4 targets pinned',
    allCut: 'is_cut',
    allHdr: 'is_cut, is_hdr, indel_length',
  },
  {
    property: 'consistency on a clean seed',
    allCut: '0.20 – 0.23 measured',
    allHdr: 'exactly 1.0',
  },
  {
    property: 'clean seeds',
    allCut: '12.7 – 29.2% of 900',
    allHdr: '15–16 of 100 (6–10 on HEK293)',
  },
  { property: 'group size and cas mix', allCut: '42, an 83/17 mix', allHdr: '80, an 80/170 mix' },
  { property: 'Cas12a GC band', allCut: '0.40 – 0.60', allHdr: '0.40 – 0.95' },
  { property: 'bank scan cost', allCut: '~120 s at d400', allHdr: '~35 – 40 s' },
  { property: 'decline gates', allCut: 'four', allHdr: 'five' },
];

/** One step of the shared construction machine. */
export interface BuildPhase {
  tag: string;
  title: string;
  source: string;
  summary: string;
  details: string[];
  icon: string;
}

/**
 * The construction both spike builders run, in order. The two-stage shape —
 * min-union one half, then require strictness of the other half over only what
 * survived — is what makes the construction reachable at all.
 */
export const BUILD_PHASES: BuildPhase[] = [
  {
    tag: '0',
    title: 'Reproduce the validator’s draw, exactly',
    source: 'genomics/mt19937.py, genomics/sha256_gpu.py',
    summary:
      'Everything above rests on one question asked millions of times: does this guide, at this seed, produce this outcome? Answering it means reproducing the validator’s own RNG bit for bit, which is why two replica modules exist.',
    details: [
      'Stage 3 seeds one generator per row from sha256 of the round seed joined to the design fields, keeps the low 32 bits, and then draws in a fixed order: the microhomology coin, the cut coin, the repair mode, the indel length. One guide against one seed therefore costs a full Mersenne Twister initialisation.',
      'The replica batches those initialisations through numpy — the same work as about 1900 array operations rather than 1900 times the batch in scalar ones. Exactness is the contract: CPython seeds an integer with init_by_array, not init_genrand, and those are not interchangeable with numpy’s legacy generator.',
      'The hashing then becomes the binding cost, so it moves to the GPU. The trick is not the hash but what crosses the bus: the key splits at its first separator, so only a per-seed digit table and a per-guide suffix table are uploaded, and each thread assembles its own message in registers.',
      'Both carry a verify() against the real pipeline. Run it after touching either — a one-bit divergence would not crash anything, it would quietly return the wrong guides and silently invalidate every bank built since.',
    ],
    icon: 'pi pi-hashtag',
  },
  {
    tag: '1',
    title: 'Screen the Cas12a candidates into a bank',
    source: 'all_cut.build_bank, mt19937.screen_guides_rule_gpu',
    summary:
      'Every (Cas12a site, mutation) pair within max_distance is enumerated into guide variants, and each variant is tested against every seed in the window. A guide that breaks the rule on at most max_fail seeds is banked together with its failed-seed array.',
    details: [
      'The bank is the expensive half and it is seed-independent — it is defined over the whole window — so one bank serves every task of the same contract shape. bank_key hashes the cell type, accessibility, mutation set and regions, GC bands, distance, variant cap, max_fail and the seed range.',
      'In practice the disk cache rarely hits across rounds, because every archived task has a distinct mutation set and the experiment key hashes the mutation string. It is kept because it costs nothing and covers a re-broadcast.',
      'A completed bank is still discarded if the deadline passes on the final target, because the deadline test runs after the append for every index including the last.',
    ],
    icon: 'pi pi-database',
  },
  {
    tag: '2',
    title: 'Min-union a group out of the bank',
    source: 'fastgreedy.FastGreedy.best',
    summary:
      'Choose group_size guides minimising the size of the union of their failed seeds — not the sum. The complement of that union is the clean set: the seeds on which the whole submission’s pin holds.',
    details: [
      'Overlapping failures are free, and that is the mechanism. Once a seed enters the union it is marked covered forever, so a later candidate failing on that same seed pays nothing for it. A guide failing 60 already-lost seeds is picked ahead of a guide failing one fresh seed.',
      'That is why the construction scales at all: the group’s cost is not the sum of about 60 failures each across 42 guides, roughly 2500, but the size of their union, about 340 of 900.',
      'Restart 0 is the deterministic first-minimum build and the remaining restarts break argmin ties randomly, keeping the smallest union. Per-cell floors stop any (mutation, cas, strand) cell being starved.',
      'The selector is blind to mutation_weight, so on a contract with a heavy and a light mutation it lands near a 50/50 split. Both attempts to correct that inside the greedy — a GC tie-break and a weight tie-break — measured net negative.',
    ],
    icon: 'pi pi-filter',
  },
  {
    tag: '3',
    title: 'Fill the Cas9 half, conditional on the clean set',
    source: 'all_cut.scan_cas9, cas9_cell_target',
    summary:
      'Scan Cas9 sites nearest-first and keep only guides that satisfy the rule on every seed of the clean set — max_fail 0, strictness rather than tolerance.',
    details: [
      'Strict rows add nothing to the union, so the clean set is preserved exactly rather than eroded. That is the entire reason for the two-stage shape.',
      'Conditioning on the clean band rather than the whole band is also what makes it affordable: HDR over 15 to 16 seeds is about 1.1e-4 per guide, against about 1e-8 over the full band.',
      'The scan early-exits once it has enough candidates and every cell has reached its quota. Requiring cells to be merely non-empty was a real defect: on one HEK293 round the exit fired at job 34 of 395 with the heavy-plus cell holding 4 candidates against 5531 available, and mean mutation weight came out 0.784 against a reachable 0.939.',
    ],
    icon: 'pi pi-lock',
  },
  {
    tag: '4',
    title: 'Assemble the rows',
    source: 'all_cut.assemble',
    summary:
      'Bucket the Cas9 pool by (mutation, strand), score a per-cell head with the validator’s own stage-2 code, apportion the rows across mutations by weight, and emit deduplicated on (cas, start, strand, guide).',
    details: [
      'Candidates are scored per cell rather than globally. A global cap dropped whole stage-5 cells once the mix was balanced, and an empty (mutation × cas × strand) cell costs roughly a 0.03x multiplier on the entire score.',
      'The quota can run dry, so the row set is topped up from the unused head by score. Without that top-up a short submission loses term 1 linearly, which silently produced a 137-row build during development.',
      'assemble is shared by both builders while their config dataclasses have different fields, so it reads optional fields through getattr with a default. A hard attribute access here once broke every all-HDR build, and would have taken the live hotkeys down at their next restart.',
    ],
    icon: 'pi pi-th-large',
  },
];

/**
 * The four window layouts in window_plan.py, checked in this order, first
 * non-empty one wins. Only the first ships at 527588a.
 */
export interface WindowLayout {
  mode: string;
  condition: string;
  layout: string;
  active: boolean;
}

export const WINDOW_LAYOUTS: WindowLayout[] = [
  {
    mode: 'joined',
    condition: 'JOINED_MODE = true',
    layout:
      'Three predicted 100-seed windows are concatenated into one 300-seed space. Ten hotkeys take circular width-225 slices at stride 30, so every seed sits in 7.5 slices on average, and one hotkey takes the whole 300.',
    active: true,
  },
  {
    mode: 'fixed',
    condition: 'FIXED_WINDOWS non-empty',
    layout:
      'Literal per-hotkey windows, overlap permitted, the prediction not consulted. Eleven width-200 windows at stride 10 over the span 100–399.',
    active: false,
  },
  {
    mode: 'explicit',
    condition: 'RANK_BY_HOTKEY non-empty',
    layout: 'One hotkey per ranked window at the full width of 100.',
    active: false,
  },
  {
    mode: 'concentrate / spread',
    condition: 'otherwise',
    layout:
      'Several hotkeys tile the top-ranked window at band width; the rest spread over the other ranks at width 100.',
    active: false,
  },
];

/** Facts about bands that a window layout hides. */
export interface BandFact {
  claim: string;
  evidence: string;
}

export const BAND_FACTS: BandFact[] = [
  {
    claim: 'A band leaks outside its window',
    evidence:
      'The window only guides the min-union search; it does not bound the resulting clean band. A hotkey on window 300–599 produced the band 132, 140, 144, 163, 173, 511, 517, 530, 546, 547, 573, 817 — five seeds below the window and one above it. So sibling bands overlap more than a disjoint plan implies, and any coverage union computed from windows is an upper bound.',
  },
  {
    claim: 'Band size falls as the window widens',
    evidence:
      'Measured 13, 12, 11 and 9 seeds at widths 100, 150, 200 and 300. Band size is set by how many rows must agree, not by how wide a range was searched, so widening a window to cover more seeds buys less than it appears to.',
  },
  {
    claim: 'Overlapping windows do not re-correlate siblings',
    evidence:
      'At stride 10 adjacent width-200 windows share 190 of their 200 seeds, yet measured pairwise band overlap is 0.56 seeds against an independent-draw expectation of 0.58. A 10-seed shift changes which guides survive the max_fail tail cut, so the greedy lands somewhere else.',
  },
];

/** Hypotheses tested and dropped, so they are not re-run. */
export interface Falsified {
  hypothesis: string;
  result: string;
}

export const FALSIFIED: Falsified[] = [
  {
    hypothesis: 'A wider screening window gives a wider band',
    result:
      'Band is 13 / 12 / 11 / 9 at widths 100 / 150 / 200 / 300, and 11 to 16 over the full 900. Band size is set by how many rows must agree, and a wider window costs hotkeys.',
  },
  {
    hypothesis: 'One construction that pins cut and HDR together',
    result:
      'Dead at every setting tried. HDR is a subset of cut, so it is a sequential filter rather than a combination, and the filter shrinks the Cas12a pool — a smaller pool min-unions to a larger failed-seed union.',
  },
  {
    hypothesis: 'A hybrid split, part of the rows on each rule',
    result: 'No elevated middle regime. Its spike consistency measured below the pure floor.',
  },
  {
    hypothesis: 'A weaker rule for a wider band',
    result:
      'not_mhnhej reaches band 45, three and a half times wider and real, but pins only one target, so consistency is 0.123. All three seeds in band still scores 25.7, under every cutoff. Expected payout exactly zero.',
  },
  {
    hypothesis: 'Narrowing all-cut’s window to widen its clean set',
    result:
      'The clean fraction does rise, to 93.6 percent at width 300, but expected final falls — 37.62 against 46.11 for the whole window. Consistency has a second channel independent of coverage, which is row composition.',
  },
  {
    hypothesis: 'Fidelity as the binding term',
    result:
      'Group 125 reached fidelity 0.976, above the leaders’ 0.947, with all eight cells occupied, and still lost. Consistency fell faster than fidelity rose.',
  },
  {
    hypothesis: 'Stage 4’s overfitting as an exploit',
    result:
      'Real but not buildable. Collapsing the feature matrix moves off-band avg_r2 from −0.249 to −0.032, confirming the negative R² is pure overfitting, but max(avg_r2, 0) discards it. The recovery needs distance cardinality of 8 or less, and 250 rows force about 31 distinct values.',
  },
  {
    hypothesis: 'Two band hits explain the field’s top block',
    result:
      'That needs a band of about 55 seeds to produce the observed 2.52 percent plateau rate; at band 12 to 16 the prediction is 0.13 to 0.23 miners per round against about 6 observed. The plateau is one band hit sitting on an all-cut clean floor.',
  },
  {
    hypothesis: 'Small submissions for a degenerate R² of 1.0',
    result:
      'Every one of 520 top-10 finishers in the three-seed regime submits exactly 250 rows. The few smaller submissions score 0.064 to 0.093.',
  },
];

/**
 * The seed-window prediction, as its own committed artifacts measure it.
 *
 * Two attempts exist. The newest is a 151,561-parameter transformer, one
 * finetuned checkpoint per cell type, predicting the three 100-seed classes the
 * next round of that cell type will draw. Its output is live and load-bearing:
 * the window plan reads it and the three predicted classes become the joined
 * space the fleet tiles.
 *
 * Numbers here come from the package's own `report.json`, which supersedes both
 * its README and the caller's comment — those quote an earlier run (0.886 hits,
 * p 0.869) that appears in no committed file. The conclusion is unchanged: no
 * detectable edge. What changes is the sign. The "chance = 1.000" those places
 * compare against is mis-specified, because a round can draw the same class
 * twice while the scorer caps a class at one hit; on the actual test set chance
 * is 0.9024, which puts the model nominally *above* it rather than below.
 */
export const SEED_MODEL = {
  parameters: 151561,
  testRounds: 82,
  /** Mean classes hit, out of the three a round draws. */
  modelHits: 0.9756,
  bestReferenceName: 'cold hand',
  bestReferenceHits: 0.9878,
  /** Chance recomputed against each round's own multiset, rather than assumed. */
  chanceHits: 0.9024,
  /** Model against the best reference, one-sided. */
  pValue: 0.579,
  logLoss: 2.2158,
  uniformLogLoss: 2.1972,
  /** Rounds of the test set that drew a duplicated class, capping the score at 2 of 3. */
  duplicateRounds: 23,
};

/**
 * The timing envelope. Every construction on the branch is far too slow to
 * build inside the presigned URL's 300-second TTL, so the build was moved off
 * that path entirely and the TTL now only has to cover the upload.
 */
export const ENVELOPE = {
  /** settings.SUBMISSION_TIMEOUT — the TTL the deadline is clamped to. */
  ttlSeconds: 300,
  /** Miner.PREPARE_BUDGET_S — what a prefetched build is allowed to take. */
  prepareBudget: 900,
  /** Reserved out of the deadline for the PUT itself. */
  uploadReserve: 45,
  /** Reserved on top of that for one ordinary construction, if a prepare overruns. */
  emergencyBuild: 60,
  /** Subtracted from the deadline to form the in-TTL hedge budget. */
  hedgeReserve: 75,
  /**
   * Measured lead time from task creation to the validator's call, in seconds,
   * over 72 rounds. `observedFloor` is a later single observation that broke
   * the measured minimum: two hotkeys were called 129 s and 149 s after task
   * creation, which is inside every construction's build time. That is why the
   * last rung of the ladder has to always finish.
   */
  lead: { min: 198, p10: 395, median: 1794, max: 4404, rounds: 72, observedFloor: 129 },
};

/** The per-rung budget gates, which is what decides who can build in the TTL. */
export interface BudgetGate {
  rung: string;
  seconds: string;
  reachableInTtl: boolean;
  note: string;
}

/**
 * The in-TTL budget is the deadline less the 75-second reserve, so roughly
 * 215-220 seconds once the artifact fetches are paid for. Three gates sit above
 * that and are therefore prefetch-only — by design, not by accident.
 */
export const BUDGET_GATES: BudgetGate[] = [
  {
    rung: 'seed-depend',
    seconds: '360',
    reachableInTtl: false,
    note: 'The build measures 285 to 333 s at the shipped variant count, so a flat 190 s gate would start a build that could not finish, burn the window and fall through the ladder poorer than if it had never tried.',
  },
  {
    rung: 'all-HDR',
    seconds: '190',
    reachableInTtl: true,
    note: 'One flat number, because all four cell types build well inside the in-TTL path.',
  },
  {
    rung: 'all-cut, HEK293',
    seconds: '190',
    reachableInTtl: true,
    note: 'Kept low deliberately so HEK293 is not locked out of a path it completes inside.',
  },
  {
    rung: 'all-cut, erythroid',
    seconds: '480',
    reachableInTtl: false,
    note: 'The wider bank these cell types moved to is 4 to 6 times slower. That makes their all-cut prefetch-dependent: a round whose prefetch fails falls to the next rung instead.',
  },
  {
    rung: 'seed-agnostic',
    seconds: '150',
    reachableInTtl: true,
    note: 'Its own scan is capped at 210 s regardless of the budget offered, so building earlier cannot silently widen it.',
  },
  {
    rung: 'ordinary construction',
    seconds: 'none',
    reachableInTtl: true,
    note: 'One to five seconds with the caches warm. This is what the emergency reserve pays for.',
  },
];

/**
 * Wiring hazards found by reading the code, each of which changes what an
 * operator should expect. None of them is in the branch's notes.
 */
export interface SharpEdge {
  title: string;
  detail: string;
}

export const SHARP_EDGES: SharpEdge[] = [
  {
    title: 'The prefetch budget is not a wall clock',
    detail:
      'After a rung declines, the ladder refreshes its budget — but on the prefetch path that refresh returns the configured 900 s verbatim rather than what is left of it. A prepared round that spends 700 s failing all-HDR still tells all-cut it has 900 s. The budget is a per-builder allowance, not a total, and nothing enforces the total.',
  },
  {
    title: 'Local scoring covers only the fallback path',
    detail:
      'The predicted-score line comes from replaying the validator’s five stages over the built rows, and it sits inside the ordinary construction. All four hedge rungs return before reaching it, so a round that ships all-HDR, all-cut, seed-depend or seed-agnostic rows gets no predicted score at all — only its builder’s own metadata. It is log-only either way and can never change what is uploaded.',
  },
  {
    title: 'One reserve protects a fallback that will refuse',
    detail:
      'The seed-agnostic minimum budget is reused as the reserve all-cut must leave behind for its caller, which suppresses all-cut’s internal retry. But seed-agnostic is skipped for HEK293 and on any host without a usable GPU, so in those cases all-cut gives up a recoverable retry to protect a rung that would decline outright.',
  },
  {
    title: 'A transient GPU error is permanent',
    detail:
      'The GPU probe memoises its result in a module global for the life of the process. A transient CUDA failure at the first probe therefore disables the seed-agnostic rung until restart — and a transient CUDA error is one of the causes the prepare-retry logic was written for. The other rungs have no GPU gate.',
  },
  {
    title: 'The plan file’s TTL is a string comparison',
    detail:
      'Staleness is decided by comparing ISO-8601 timestamps as text, which is correct only while both sides are the same shape of UTC string. A plan written with an offset suffix or a space separator would compare wrongly. An absent or empty expiry disables the check entirely, and every failure path falls back to the hotkey’s own window pin — a prediction can never stop a build.',
  },
  {
    title: 'A stale prediction is used, not rejected',
    detail:
      'The window plan checks whether the seed model’s prediction has been overtaken by a round that has since drawn, and prints STALE when it has — then uses the prediction anyway. Only an unreadable or malformed file causes a fallback. Harmless while the layout is expected-value-neutral, but the plan output must not be read as though a fallback occurred.',
  },
  {
    title: 'The model’s promotion guard cannot detect "no edge"',
    detail:
      'After each refresh the new checkpoint is kept unless it scores more than 0.10 hits below the previous cycle’s own score. It never consults the baselines, chance, or the p-value — so a model that has never beaten chance is promoted every cycle. The guard detects degradation against yesterday’s noise, which is not what the model is there to measure.',
  },
  {
    title: 'The submission file is only authoritative after the upload',
    detail:
      'A prefetched build and an in-TTL fallback can both write the same submission path, since the build lock cannot cover the prefetch thread. The upload path repairs it by re-persisting what was actually sent, so a crash between build and upload can leave rows on disk that no validator ever saw.',
  },
];

/** What the branch treats as decided, and what it knows it has not answered. */
export interface StatusItem {
  claim: string;
  evidence: string;
}

export const SETTLED: StatusItem[] = [
  {
    claim: 'The off-band floor is universal',
    evidence:
      '0.092 to 0.101 for every 250-row submission tried, and the population median in the three-seed regime is 0.101. The difference between this fleet and the field is how many seeds are not off-band, not how high the floor is.',
  },
  {
    claim: '250 rows, always',
    evidence:
      'All 520 top-10 finishers in the three-seed regime submit exactly 250 rows; the few smaller submissions score 0.064 to 0.093.',
  },
  {
    claim: 'The band cannot be widened',
    evidence:
      'Each seed added multiplies the conditional Cas9 requirement by the probability of HDR, about 0.57 and about 0.37 on HEK293. A 12 to 16 seed band leaves enough candidates; a 29 to 31 seed band yielded zero on all three erythroid cell types.',
  },
  {
    claim: 'The band’s position is free',
    evidence:
      'Verified across all nine 100-seed windows on HEK293: band 7 and all eight stage-5 cells at every one, fidelity 0.904 to 0.934. This is what makes an unproven window prediction safe to run, and what makes the fleet layout arbitrary but sound.',
  },
  {
    claim: 'The two scoring eras are not comparable',
    evidence:
      'Single-seed rounds ran until 2026-08-24, when the seed was knowable and almost the whole field reached 1.000. No three-seed round has exceeded 0.875. Pooling the two invents a mechanism that does not exist, and filtering on a score row’s own timestamp is not enough, because single-seed tasks were still being scored after the switch.',
  },
];

export const OPEN: StatusItem[] = [
  {
    claim: 'HEK293 has never had the matched pricing treatment',
    evidence:
      'The honest method builds both arms on the same contract and places each against the one field that played that contract. Done for K562, it reversed the pooled answer by 3.78x. HEK293 is still on a pooled estimate, and it is the cell type whose accessibility makes the two arms structurally different.',
  },
  {
    claim: 'Point prediction is a different problem, and nothing has measured it',
    evidence:
      'The construction needs the seeds to within about 3 of 900; both models predict 100-seed windows. An oracle build pinned to two of a round’s three seeds scored rank 1 of 245, so the prize is real — but a blind band inside a correctly predicted window is worth about one chance in 2850.',
  },
  {
    claim: 'Hotkey count is the only lever with headroom',
    evidence:
      'Every other lever measured is capped. The leading coldkeys run 11 to 14 hotkeys each with placements roughly equal to their hotkey count, which is the signature of the same construction played wider.',
  },
  {
    claim: 'The weighted-score gap to the leaders is unexplained',
    evidence:
      'Their placers’ median total weighted score is 263 against this branch’s 230, and the gap sits in the structural half where the off-target factor is already a perfect 1.0 and mean GC is already 0.505. The mutation-skew route is measured and dead; tightening GC and distance raises weighted 2 to 3 percent and costs band on three of four conditions.',
  },
];

/** Where to look in the branch, for a reader who wants the code. */
export interface FileRef {
  path: string;
  role: string;
}

export const FILES: FileRef[] = [
  { path: 'neurons/miner.py', role: 'request path, the ladder, prefetch, upload' },
  { path: 'niome_subnet/genomics/all_hdr.py', role: 'the narrow-band HDR pin' },
  {
    path: 'niome_subnet/genomics/all_cut.py',
    role: 'the whole-window cut pin, and the shared assembly',
  },
  { path: 'niome_subnet/genomics/seed_depend.py', role: 'the seed-0 build' },
  {
    path: 'niome_subnet/genomics/seed_agnostic.py',
    role: 'the cut-only hedge over a whole window',
  },
  { path: 'niome_subnet/genomics/fastgreedy.py', role: 'the min-union group selector' },
  { path: 'niome_subnet/genomics/hek293_generation.py', role: 'HEK293’s own clustered builder' },
  {
    path: 'niome_subnet/genomics/mt19937.py, sha256_gpu.py',
    role: 'bit-exact CPU and GPU replicas of the validator’s RNG',
  },
  {
    path: 'genExp.py',
    role: 'site, guide and scoring primitives, calling the validator’s own stages',
  },
  { path: 'miner.sh', role: 'the fleet launcher and its four hotkey lists' },
  { path: 'window_plan.py, round_plan.sh', role: 'window allocation, refreshed hourly' },
  { path: 'seed_model/', role: 'the seed-window transformer and its evaluation' },
  {
    path: 'niome_subnet/utils/settings.py',
    role: 'the payout curve, the schedule, every data path',
  },
];

/** Where the branch's own notes no longer match its code. */
export interface Drift {
  notes: string;
  code: string;
}

export const DRIFT: Drift[] = [
  {
    notes:
      'The fleet layout is FIXED_WINDOWS: seven hotkeys on overlapping width-300 windows stepping by 100.',
    code: 'JOINED_MODE is checked first, so FIXED_WINDOWS is dead configuration. Ten hotkeys take circular width-225 slices at stride 30 of a joined 300-seed space.',
  },
  {
    notes:
      'ALL_CUT_HOTKEYS="niome_hotkey" puts h0 on all-cut for every cell type, so ten hotkeys play a band.',
    code: 'ALL_CUT_HOTKEYS is empty in miner.sh. No hotkey is forced onto all-cut, and JOINED_FULL_HK puts h0 on the whole joined window instead — a change window_plan.py explains and the fleet notes do not.',
  },
  {
    notes: 'h8, h9 and h10 complete a stride-10 tiling of the window space.',
    code: 'The HOTKEYS table pins h8, h9 and h10 all to 900–999, and h0 and h1 both to 100–299. Those literals only matter as the fallback when the plan file is missing, but as a fallback they leave three siblings fully correlated.',
  },
  {
    notes: 'Nothing in CLAUDE.md mentions the seed_model package.',
    code: 'The branch HEAD adds seed_model/, a per-cell-type transformer whose prediction window_plan.py consumes through JOINED_SOURCE. Its own measurements put it at chance.',
  },
  {
    notes:
      'The seed model scores 0.886 hits of three against chance’s 1.000, p = 0.869, and its log loss is worse than uniform’s.',
    code: 'Those figures come from the package README and appear in no committed artifact. Its own report.json says 0.976 hits against the best reference’s 0.988, p = 0.579, over 82 rounds. The "chance = 1.000" both places compare against is also mis-specified, because a round can draw one class twice while the scorer caps a class at one hit — computed against each round’s own draw it is 0.902, which puts the model nominally above chance rather than below. The verdict of no edge survives; the sign in the prose does not.',
  },
];

/** One entry of the glossary. */
export interface GlossaryTerm {
  term: string;
  /** Other spellings the page uses for the same thing, searched alongside the term. */
  aka?: string;
  definition: string;
}

export interface GlossaryGroup {
  title: string;
  blurb: string;
  terms: GlossaryTerm[];
}

/**
 * Every term this page uses, defined for a reader who knows neither the subnet
 * nor the biology.
 *
 * Definitions carry as few measured numbers as possible on purpose. The figures
 * live in the constants above and in the prose that cites them; repeating one
 * here would be a second place to update and a second place to be wrong. A
 * number appears only where it is part of the definition — three seeds, 900
 * values, 250 rows, ten paid ranks, a 300-second link.
 */
export const GLOSSARY: GlossaryGroup[] = [
  {
    title: 'The round, and who gets paid',
    blurb: 'The shape of the game, before any of the strategy.',
    terms: [
      {
        term: 'subnet',
        definition:
          'One market inside Bittensor. Miners do a job, validators score the results, and the network pays out on those scores. This subnet asks for CRISPR gene-editing experiment designs.',
      },
      {
        term: 'miner',
        definition:
          'One registered identity that uploads a file each round. Here one miner is one hotkey running one process.',
      },
      {
        term: 'validator',
        definition:
          'The scoring side. It collects each miner’s file, replays a fixed five-stage pipeline over it, and reports one number per miner.',
      },
      {
        term: 'round',
        definition:
          'One cycle: a task is published, miners upload their designs, the validator scores them, the top ten are paid.',
      },
      {
        term: 'task',
        aka: 'contract',
        definition:
          'The job sheet for a round — which cell type, which mutations to target, what each mutation is worth, and the rules everyone builds under. Every miner in the round is handed the same one.',
      },
      {
        term: 'round seed',
        definition:
          'The three numbers, each between 100 and 999, that the validator feeds to its simulation. They are stamped onto the task only after the round has closed, so no miner can design against them. Everything on this page follows from that one fact.',
      },
      {
        term: 'seed space',
        definition:
          'The 900 values a seed can take, 100 through 999. A round draws three of them independently and uniformly.',
      },
      {
        term: 'stamped',
        aka: 'unstamped',
        definition:
          'A task carries seed 0 while its round is open — unstamped — and its real seeds afterwards. A round the backend never stamps is scored at seed 0, which is the one case where the miner knew the seed all along.',
      },
      {
        term: 'submission',
        definition: 'The file a miner uploads: up to 250 rows, one proposed experiment per row.',
      },
      {
        term: 'rank',
        aka: 'placement, places',
        definition:
          'Where a miner’s final score falls among everyone who submitted that round. Only the first ten are paid anything at all.',
      },
      {
        term: 'payout curve',
        aka: 'SCORE_DISTRIBUTION, SCORING_SYSTEM top',
        definition:
          'The fixed split of a round’s pay by rank: a large share for first, falling away to a token share for tenth, and nothing from eleventh down. Because it is a step rather than a slope, a score gain that crosses no rank boundary pays exactly nothing.',
      },
      {
        term: 'emission',
        definition: 'The tokens a round pays out, divided among the top ten by the payout curve.',
      },
      {
        term: 'rank-10 cutoff',
        aka: 'the cutoff',
        definition:
          'The tenth-placed score in a round, and so the bar to clear to be paid. It is not a constant — it moves with how strong that round’s field happened to be, which is why the page quotes it as a distribution.',
      },
      {
        term: 'three-seed regime',
        aka: 'the two scoring eras',
        definition:
          'Rounds since 2026-08-24, which draw three seeds. Before that a round drew one, it was knowable, and almost the whole field scored a perfect consistency — so scores from the two periods cannot be compared or pooled.',
      },
    ],
  },
  {
    title: 'How a submission is scored',
    blurb: 'The five stages, and the three factors that multiply into the ranked number.',
    terms: [
      {
        term: 'the five stages',
        definition:
          'The validator’s pipeline. In short: stage 2 scores each row’s design quality, stage 3 simulates the experiment at each seed, stage 4 checks the simulated outcomes are learnable, stage 5 checks the rows cover the design space.',
      },
      {
        term: 'total_weighted_score',
        aka: 'weighted',
        definition:
          'Stage 2’s sum over the rows: design quality, weighted by how much the contract values each mutation. The first of the three factors multiplied into the final score.',
      },
      {
        term: 'stage 4',
        definition:
          'Trains a small random forest to predict three outcome columns — is_cut, is_hdr, indel_length — from each row’s design fields under a five-fold split, then scores those predictions. This is where consistency_factor comes from.',
      },
      {
        term: 'R²',
        aka: 'r2_score, avg_r2',
        definition:
          'How much of a column’s variation a prediction explains: 1.0 is perfect, 0 is no better than always guessing the average, and below 0 is worse than that. The quirk the whole strategy turns on: a column that never changes, predicted exactly, returns exactly 1.0.',
      },
      {
        term: 'normalized MAE',
        aka: 'avg_nmae, normalized_mae',
        definition:
          'The average prediction error divided by how spread out the column is. When that spread is essentially zero the code returns the raw error instead of dividing — the second half of the same quirk.',
      },
      {
        term: 'consistency_factor',
        aka: 'consistency',
        definition:
          'Stage 4’s two scores combined and clamped into 0 to 1, then averaged over the round’s three seeds. In practice it measures how many of the three outcome columns every row holds identical, not how learnable they are — which is why it can be engineered rather than earned.',
      },
      {
        term: 'distribution_fidelity_factor',
        aka: 'fidelity',
        definition:
          'How closely the spread of the submission’s simulated outcomes matches what the validator expects to see. Measured as not the binding constraint: builds with better fidelity than the leaders still lost on consistency.',
      },
      {
        term: 'final_score',
        definition:
          'total_weighted_score × consistency_factor × distribution_fidelity_factor. The one number miners are ranked on.',
      },
      {
        term: 'the floor',
        aka: 'off-band floor, floor seed',
        definition:
          'The consistency a full 250-row submission scores on a seed where it has pinned nothing — near 0.10 for everything anyone has tried, this fleet and the wider field alike. So the difference between miners is how many seeds are not floor seeds, not how high their floor is.',
      },
      {
        term: 'overfitting',
        definition:
          'A model memorising noise instead of finding a pattern, which shows up as a negative R². Stage 4 floors R² at 0, so overfitting quietly costs the score rather than showing up as an error.',
      },
      {
        term: 'stage-5 cells',
        aka: 'cells, empty cells',
        definition:
          'The buckets stage 5 wants covered — each combination of mutation, enzyme and strand. Leaving one empty is close to fatal: stage 5 takes a geometric mean, so one empty bucket crushes the whole score.',
      },
      {
        term: 'off-target factor',
        definition:
          'The penalty a guide carries for also matching somewhere else in the genome. Already perfect on this branch, so it is not where the remaining gap to the leaders is.',
      },
    ],
  },
  {
    title: 'What a row is made of',
    blurb: 'The biology, only as far as the scoring actually uses it.',
    terms: [
      {
        term: 'row',
        aka: 'experiment',
        definition:
          'One proposed edit: which enzyme, which guide sequence, where it lands, which mutation it targets. 250 rows make a submission.',
      },
      {
        term: 'guide',
        definition:
          'The short RNA sequence that steers the cutting enzyme to one specific place in the DNA. Choosing guides is the entire design job.',
      },
      {
        term: 'site',
        definition:
          'A place in the DNA where an enzyme can land. Each site yields many candidate guides.',
      },
      {
        term: 'variant',
        definition:
          'One candidate guide generated at a site. The builders enumerate thousands per site and screen every one of them against every seed.',
      },
      {
        term: 'Cas9 / Cas12a',
        aka: 'cas systems, cas mix',
        definition:
          'The two cutting enzymes the rules allow. They behave differently in the simulation: Cas9 can be pushed to cut on almost every seed, Cas12a cannot. The constructions exploit that by treating the two halves of the submission differently.',
      },
      {
        term: 'strand',
        definition:
          'Which of DNA’s two complementary strands a guide binds to. Part of the coverage buckets stage 5 checks.',
      },
      {
        term: 'distance',
        aka: 'max_distance',
        definition:
          'How far a guide’s cut sits from the mutation it is meant to repair. The rules cap it, and the builders scan nearest-first.',
      },
      {
        term: 'GC content',
        aka: 'GC band',
        definition:
          'The share of G and C letters in a guide, which affects how tightly it binds. The builders accept only guides inside a chosen band.',
      },
      {
        term: 'cut',
        aka: 'is_cut',
        definition:
          'Whether the simulation says the enzyme actually cut, for one row at one seed. The first of stage 4’s three target columns, and the only one a design can hold constant across many seeds.',
      },
      {
        term: 'HDR',
        aka: 'is_hdr, homology-directed repair',
        definition:
          'One of the ways a cell repairs a cut — copying from a template, rather than gluing the ends back together. The simulator picks a repair mode per row per seed at roughly a coin flip and no design choice moves it, which is why pinning HDR needs the right seed rather than the right guide.',
      },
      {
        term: 'MH-NHEJ',
        aka: 'not_mhnhej, microhomology end joining',
        definition:
          'Another repair mode. It appears here only as a weaker screening rule that was tried, gave a much wider band, and still paid nothing.',
      },
      {
        term: 'indel_length',
        definition:
          'How many DNA letters the repair inserted or deleted. The third of stage 4’s target columns.',
      },
      {
        term: 'cell type',
        definition:
          'Which cell line the round is for: HEK293, K562, HUDEP-2 or CD34+_HSPC. It is fixed by the task, and it changes which construction can work.',
      },
      {
        term: 'erythroid',
        definition:
          'Shorthand on this page for the three blood-lineage cell types, as against HEK293 — which has much tighter chromatin and so behaves differently at almost every step.',
      },
      {
        term: 'accessibility',
        definition:
          'How open the DNA is in that cell type, on a 0 to 1 scale. It scales the chance of a cut, so it sets how much of the score is reachable at all before any design choice is made.',
      },
      {
        term: 'mutation weight',
        definition:
          'How much the contract values each mutation in the round. Rows are apportioned across mutations in proportion to it.',
      },
    ],
  },
  {
    title: 'Pinning, bands and constructions',
    blurb: 'The vocabulary of the strategy itself. Most of the page is written in these words.',
    terms: [
      {
        term: 'pin',
        definition:
          'To make one of stage 4’s target columns identical in every row, at a given seed. A pinned column scores exactly 1.0 — that is the whole mechanism this page is built on.',
      },
      {
        term: 'spike',
        definition:
          'A round where a pin holds on at least one of the three seeds, so consistency jumps from the floor to something that can place.',
      },
      {
        term: 'clean seed',
        definition: 'A seed on which a given build’s pin holds.',
      },
      {
        term: 'band',
        aka: 'clean band, clean set',
        definition:
          'All the clean seeds for one build. Its width is the strategy’s central quantity: the wider the band, the more often the build spikes. Widening it is settled as impossible — each extra seed multiplies the requirement by another coin flip.',
      },
      {
        term: 'screening rule',
        definition:
          'The one string that separates the two spike constructions: "cut" accepts any repair mode, "hdr" demands HDR. Everything else about them is the same code.',
      },
      {
        term: 'bank',
        definition:
          'Every candidate guide that survived the screen, stored with the list of seeds it failed on. It is the expensive half of a build and does not depend on the round’s seeds, so it can be cached and reused.',
      },
      {
        term: 'max_fail',
        definition:
          'How many seeds a guide is allowed to break the rule on and still be banked. Deliberately loose, because the next step is what actually produces the clean band.',
      },
      {
        term: 'strict',
        definition:
          'Breaking the rule on no seed at all of the set in question. The Cas9 half is built strict, so it adds nothing to the group’s failures and leaves the band exactly as it found it.',
      },
      {
        term: 'min-union',
        definition:
          'The trick at the centre of the construction: pick the group of guides whose failed-seed lists overlap as much as possible. Overlapping failures are free, so a group loses far fewer seeds between them than the sum of their individual failures. What is left over is the clean band.',
      },
      {
        term: 'greedy selector',
        aka: 'fastgreedy, restarts',
        definition:
          'The algorithm that does the min-union: take the best-looking candidate, mark its failures as already lost, repeat. The first pass is deterministic and later passes break ties at random, keeping whichever result was smallest.',
      },
      {
        term: 'group_size',
        definition:
          'How many Cas12a guides go into the min-union group. Bigger group, wider coverage of the design space, narrower band — the two constructions sit at different points of that trade.',
      },
      {
        term: 'the ladder',
        definition:
          'The ordered list of builders the miner tries, best-paying first. Each one either produces rows or declines, so a failure anywhere lands on the next rung down rather than on nothing.',
      },
      {
        term: 'rung',
        definition: 'One builder on that ladder.',
      },
      {
        term: 'decline',
        definition:
          'A builder’s way of saying "not this round" — no GPU, too little time left, too small a candidate pool. It returns no rows plus a readable reason instead of raising an error, which is what makes the ladder safe to fall down.',
      },
      {
        term: 'hedge',
        definition:
          'A rung that gives up the exact-1.0 spike for something worth less per hit but hit far more often.',
      },
      {
        term: 'all-HDR',
        definition:
          'The spike build. It pins all three target columns at once, so a clean seed scores a perfect consistency — but its band covers only a sliver of the seed space. It loses on the average and is shipped anyway, because the payout curve pays placement rather than averages.',
      },
      {
        term: 'all-cut',
        definition:
          'The same machine run on the "cut" rule. It pins one column instead of three, worth much less per clean seed, but holds it across a far larger share of the seed space.',
      },
      {
        term: 'seed-depend',
        definition:
          'The build for a round that never gets stamped. The validator then scores at seed 0 — the value the miner was handed — so the seed is known and consistency is exactly 1.000. It is worthless on any round that does get stamped.',
      },
      {
        term: 'seed-agnostic',
        definition:
          'The older cut-only hedge: Cas9 guides that cut under every seed in the window, with no band involved. Skipped on HEK293, where tight chromatin means no such guide exists.',
      },
      {
        term: 'ordinary construction',
        definition:
          'The plain build from before any of this existed. It pins nothing and scores the floor, but it always finishes in seconds, which is why it is the last rung and the emergency path.',
      },
      {
        term: 'oracle build',
        definition:
          'An experiment built with the round’s real seeds known. Not playable — it exists only to price what perfect prediction would be worth, and it took rank 1.',
      },
    ],
  },
  {
    title: 'The fleet, and the seed prediction',
    blurb: 'Why the strategy is run across many identities at once, and what the model does.',
    terms: [
      {
        term: 'hotkey',
        definition:
          'One registered miner identity: one process, one seed window, one submission per round.',
      },
      {
        term: 'coldkey',
        definition:
          'The wallet that owns a group of hotkeys. Payout is what matters across the whole group, so the useful quantity is the union of what its hotkeys cover — not how good any single one is.',
      },
      {
        term: 'fleet',
        definition:
          'The hotkeys run together as one strategy, each searching a different part of the seed space so their bands land in different places.',
      },
      {
        term: 'window',
        definition:
          'The slice of the seed space one hotkey searches in. It steers the search but does not fence in the result: measured bands leak outside their own window.',
      },
      {
        term: 'window plan',
        definition:
          'The file, rewritten hourly, that assigns each hotkey its window. The miner re-reads it for every build and falls back to the hotkey’s own pinned window when it is missing, stale or malformed.',
      },
      {
        term: 'joined mode',
        definition:
          'The window layout that actually ships: the seed model’s three predicted classes are concatenated into one space, and the hotkeys take overlapping slices of that.',
      },
      {
        term: 'correlated siblings',
        definition:
          'Hotkeys that build near-identical rows. They score within noise of each other and so take consecutive ranks, which means the fleet moves as a block — leverage in both directions, given how steep the payout curve is.',
      },
      {
        term: 'seed model',
        definition:
          'A small transformer that tries to predict which parts of the seed space the next round will draw from. Its own evaluation puts it at chance. Running it anyway costs nothing, because where a band sits does not affect what it is worth.',
      },
      {
        term: 'band position is free',
        definition:
          'The measured fact that makes an unproven prediction safe to act on: seeds are drawn uniformly, so a band in one part of the seed space is worth exactly as much as a band anywhere else.',
      },
      {
        term: 'walk-forward',
        definition:
          'Testing a predictor round by round in time order, training only on rounds that came before the one being predicted. The honest way to score a forecast.',
      },
      {
        term: 'chance baseline',
        definition:
          'What guessing would score on the same test. Easy to get wrong here: a round can draw the same class twice while the scorer counts it once, so chance is not a perfect score.',
      },
      {
        term: 'log loss',
        definition:
          'A score for a probability forecast, where lower is better and confidence in a wrong answer is punished. A model scoring worse than a flat "no idea" has learned that there is nothing to learn.',
      },
      {
        term: 'p-value',
        definition:
          'Roughly, how often a gap this large would appear by luck alone. A large one means the two things being compared are indistinguishable on the evidence available.',
      },
      {
        term: 'cold hand',
        definition:
          'The best of the simple reference strategies the seed model is scored against in its own evaluation. It beats the model.',
      },
    ],
  },
  {
    title: 'Timing and delivery',
    blurb: 'Why a build that cannot finish in time is worth nothing, however good it is.',
    terms: [
      {
        term: 'presigned URL',
        aka: 'TTL, the deadline',
        definition:
          'The upload link the validator hands over when it asks for a submission, alive for 300 seconds. Every construction here takes longer than that to build, so the build had to move off this path entirely.',
      },
      {
        term: 'lead time',
        definition:
          'The seconds between a task appearing publicly and the validator asking that hotkey for its file. Usually generous, occasionally shorter than any real build — which is why the bottom rung must always finish.',
      },
      {
        term: 'prefetch',
        aka: 'prepare',
        definition:
          'Building against the public, unstamped task the moment a round opens, before any validator asks. The task is the same contract the validator will hand over later; only the seed differs, and it is 0 on both sides until the round closes.',
      },
      {
        term: 'in-TTL build',
        definition:
          'The fallback: building inside the 300-second upload link because no prepared build was ready. Only the cheaper rungs fit here, once time is reserved for the upload itself.',
      },
      {
        term: 'budget gate',
        definition:
          'The number of seconds a rung demands before it will start at all. Set per rung so a build that cannot finish is never begun — starting one and running out is worse than never trying.',
      },
      {
        term: 'emergency build',
        definition:
          'The ordinary construction, held in reserve for the case where a prepared build is not ready and the clock is nearly gone.',
      },
      {
        term: 'build lock',
        definition:
          'The exclusive lock that stops sibling hotkeys building the memory-heavy shared rows at the same time. One process builds, the others load what it wrote.',
      },
      {
        term: 'drift',
        definition:
          'Where the branch’s own engineering notes no longer describe what its code does. The last section of this page lists the cases found.',
      },
    ],
  },
];
