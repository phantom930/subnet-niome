# The MIT License (MIT)
# Copyright © 2023 Yuma Rao
# Copyright © 2025 Genomes.io
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
# documentation files (the "Software"), to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of
# the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO
# THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.


import asyncio
import contextlib
import hashlib
import json
import logging
import os
import requests
import shutil
import subprocess
import sys
import threading
import time
import traceback

from collections import Counter
from dataclasses import dataclass, field, replace
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import parse_qs, urlparse

# Add project root to Python path. This has to happen before the niome_subnet imports below —
# running this file as a script puts neurons/ on sys.path, not the repo root — and settings has to
# be the first of them, because it fixes BT_NO_PARSE_CLI_ARGS before bittensor is imported.
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

import niome_subnet.utils.settings as settings  # noqa: E402

# genExp is the common generator and validator-replica layer. Most cell types mirror
# submission.py's flow directly; HEK293 reuses its context, site, validation and scoring primitives
# through a dedicated clustered builder below. Note the module chdir()s to the repo root on import —
# harmless here (every settings.py path is relative to it and neurons must run from there anyway)
# but it is an import-time side effect, so it stays below settings and above nothing that cares
# about cwd.
import genExp as G  # noqa: E402

from niome_subnet.base.miner import BaseMinerNeuron  # noqa: E402
from niome_subnet.genomics.hek293_generation import (  # noqa: E402
    generate_seed_agnostic_clustered,
)
from niome_subnet.genomics import all_cut as AC  # noqa: E402
from niome_subnet.genomics import all_hdr as AH  # noqa: E402
from niome_subnet.genomics import seed_agnostic as SA  # noqa: E402
from niome_subnet.genomics.model import Task  # noqa: E402
from niome_subnet.protocol import GenomicsTaskSynapse  # noqa: E402

logger = logging.getLogger(__name__)


@dataclass
class PreparedRound:
    """One round's submission, built when its task appeared rather than when a validator asked.

    Published to ``Miner._prepared`` *before* the build starts, so a validator arriving mid-build
    finds the round and can wait on it. ``done`` is set exactly once, on success or failure, and
    which it was is then readable from ``rows`` / ``error``.
    """

    task_id: str
    key: str
    contract: dict
    reference: dict
    created_at: str
    done: threading.Event = field(default_factory=threading.Event)
    rows: list[dict] | None = None
    error: str | None = None
    started: float = field(default_factory=time.time)
    elapsed: float | None = None
    attempt: int = 1

    def summary(self) -> str:
        if not self.done.is_set():
            return f"still building ({time.time() - self.started:.0f}s so far)"
        if self.rows:
            return f"{len(self.rows)} rows in {self.elapsed:.0f}s"
        return f"failed after {self.elapsed:.0f}s: {self.error}"


class Miner(BaseMinerNeuron):
    """
    Miner neuron. Receives genomics tasks from validators via HTTP and processes them.

    The reply to a validator is an empty ack; the dataset travels out of band as a PUT to the
    presigned S3 URL that came with the task, and must land before that URL expires
    (``SUBMISSION_TIMEOUT``, 300 s from the moment the validator minted it). Generation itself is a
    second or so once the caches are warm, so the budget is not the constraint — see :mod:`genExp`
    for what is built and why it scores.
    """

    MAX_RETRIES = 3

    # Generation knobs, the same ones submission.py builds every task's dataset with. They are
    # class attributes rather than CLI flags because they are design decisions measured over the
    # backend's whole task history, not per-run operator settings.
    #
    # "pure" forces every row's outcome to satisfy CONSTRUCTION and drops rows that will not; the
    # "shaped" strategy instead fits outcomes to a distance ramp and takes whatever the draw gives.
    # Pure is what reaches consistency_factor 1.0.
    STRATEGY = "pure"
    # "packed" ranks each cell's sites by the stage-2 score they can actually reach (proximity is
    # pure total_weighted_score under the pure strategy); "stratified" spreads them over distance
    # bands to give stage 4 a feature axis, which pure does not need.
    SELECTION = "nearest"
    # Treat any guide this close to 50% GC as equally good so distance breaks the tie. Read only by
    # the "nearest" selection. Without the tier the ranking chases the last thousandth of gc_score
    # and reaches for sites hundreds of bp further out; gc_score is flat to first order around 50%,
    # dist_score is not.
    GC_TOLERANCE = 0.03
    # Which construction every row's simulated outcome is forced to satisfy. "hdr" makes every row
    # repair by HDR: all three stage-4 targets go constant, r2_score hits its zero/zero branch and
    # normalized_mae short-circuits on std < 1e-9, so consistency_factor reaches 1.0 without the
    # forest learning anything. It is also easier to hit than "mh" (P(HDR|variant) ~ 0.5 against
    # ~0.37 for the compound rule), so fewer near-mutation sites are abandoned to the reserve —
    # measured +2.3 on term 1 across every cas mix. The cost is identifiability: three degenerate
    # target columns and a cas-shift diagnostic of exactly zero. "mh" gives up ~0.8% to avoid that.
    # See G.CONSTRUCTIONS for the alternatives.
    CONSTRUCTION = "hdr"
    # Hold every Cas system to its own cut_p ceiling — 0.99 for Cas9, 0.96 for Cas12a. These differ
    # because stage 3's base rate does (0.86 vs 0.78) and energy is clamped at 1.0, so one shared
    # floor above 0.96 is satisfiable by Cas9 alone: it empties every Cas12a cell, takes
    # cas_system_coverage_entropy_ratio to 0, and the 1e-9 clip turns that into a 0.0316x
    # multiplier. Measured 8.81 against 289.45 for the per-system form.
    CUT_P_CEILING = True
    # Row share per Cas system, in rules.cas_systems order. Under a cut_p floor this acts as a cap
    # on Cas12a rather than a target: only ~168 Cas12a sites can reach 0.96, so whatever the quota
    # cannot fill is topped up by BACKFILL_CAS.
    CAS_MIX = "70/30"
    # Top the submission up to max_experiments from this Cas system, highest cut_p first, when a
    # filter empties a cell or the construction drops a row. Only Cas9 has spare sites at its own
    # ceiling. Leaving the cap short costs whole rows off term 1, which is a sum — measured +29 at
    # a 50/50 mix. None disables it.
    BACKFILL_CAS = "Cas9"
    # Guide variants searched per site. Each is a distinct stage-3 draw at an identical feature
    # vector, so this is the budget for hitting the construction without paying anything in stage 2.
    VARIANTS = 24
    # Site-enumeration window either side of gene_region. Wider only helps if cells run short of
    # coordinates, and costs a longer enumeration.
    FLANK = 3000
    # Guide lengths enumerated. Stage 1 accepts 20 and 23 only.
    LENGTHS = (20, 23)
    # Predict the validator's score after each build (~5 s on top of a ~1 s build). Worth it:
    # validators never report back, so this is the only signal a miner gets before the next task.
    SCORE_LOCALLY = True

    # ---- Post-upload archive and scoring ---------------------------------------------------
    # Where each round's artifacts are kept, one folder per upload. A submission can only be
    # scored *after* its own round: the contract is broadcast with seed 0 and the backend stamps
    # the real seed once the task closes, so stages 3 and 4 have nothing to run against until
    # then. Each folder is therefore scored one round later, when calc.py can fetch that seed
    # from /api/v3/tasks — which is what makes the archive worth keeping at all.
    # Under the per-instance data dir so several hotkeys sharing this working directory each archive
    # (and locally re-score) their own rounds instead of interleaving into one data/result.
    ARCHIVE_ROOT = f"{settings.DATA_DIR}/result"
    # Copied out of data/ into the round's folder, in this order. The settings names rather than
    # literal paths, so a relocated data/ directory takes the archive with it. last_upload.json is
    # the important one: it records the task id, which is how the folder is paired back to its
    # task — the folder name alone cannot be, since it is the archive time, not the task's.
    ARCHIVED_PATHS = ("MINER_SUBMISSION_PATH", "LAST_UPLOAD_PATH", "CONTRACT_PATH",
                      "HBB_REFERENCE_PATH")
    # calc.py's summary file. Its presence is what marks a folder as already scored.
    VALIDATION_MARKER = "validation.json"
    # calc.py runs out of process on purpose: it repoints every settings path at its own output
    # directory, which in-process would leave this miner writing its next submission into the
    # archive.
    VALIDATION_SCRIPT = "calc.py"
    # It loads the 135 MB FASTA, builds a k-mer index and fits 15 random forests under KFold —
    # a couple of minutes warm. The cap is loose because it runs long after the upload, and only
    # exists so a hung subprocess cannot hold the worker thread for ever.
    VALIDATION_TIMEOUT_S = 1800

    # Per-cell-type replacements for the constants above. Everything the tuned configuration is
    # built on assumes stage 3's energy clamp: at accessibility 0.77 and above,
    # accessibility * (1.8*gc + 0.6*exp(-d/1500) + offset) exceeds 1.0 near the mutation, so every
    # Cas9 site sits at cut_p 0.99 and every Cas12a site at 0.96. HEK293's 0.35 caps energy near
    # 0.52 — the best Cas9 site reaches 0.9545 and the best Cas12a 0.8745, so the ceilings are
    # unreachable and the cut_p filter has nothing to select on. HEK293 therefore has a dedicated
    # seed-agnostic two-pool builder. These overrides keep the generic config compatible for local
    # reporting and offline tools; the dedicated builder owns its final Cas/mutation distribution.
    #
    # 59 of the backend's 258 tasks are HEK293, so this path is a quarter of all rounds.
    CELL_TYPE_OVERRIDES = {
        "HEK293": {
            "SELECTION": "packed",       # dedicated builder applies its own Cas-specific ranking
            "GC_TOLERANCE": 0.0,         # unread by packed selection
            "CONSTRUCTION": "mh",        # config compatibility only; outcomes are not forced
            "CUT_P_CEILING": False,      # unreachable at HEK293 accessibility
            "CAS_MIX": None,             # dedicated builder balances exact global Cas totals
        },
    }

    # --- seed-agnostic bank builder (non-HEK293) ------------------------------------------------
    # Where the contract carries no seed, the outcome the ordinary construction engineers against
    # seed 0 does not hold when the validator scores under the seed it is later assigned. The bank
    # builder instead picks guides whose *cut* survives every seed in a window: strict Cas9 (cut_p
    # 0.99, so all-window guides exist) plus a min-union Cas12a group (cut_p 0.96, where none do).
    # Measured on a K562 task: mean final ~41-46 across in-window seeds against ~29 for the
    # ordinary seed-0 build, with is_cut's r2 recovering to 1.0 on ~62% of seeds.
    SEED_AGNOSTIC = True
    # HEK293 keeps its own clustered builder: at accessibility 0.35 no strict Cas9 guide exists
    # (best reachable cut_p 0.9545, so 0.9545**900 ~ 6e-19), so this method has nothing to select.
    SEED_AGNOSTIC_SKIP_CELL_TYPES = ("HEK293",)
    # Seconds held back from the upload deadline. This has to cover more than the upload: the
    # builder's time budget bounds only the two scans, while the min-union selection, assembly and
    # coverage check add ~20s after them (measured: a 254s budget produced a 274s build). So the
    # reserve is ~20s post-scan + ~15s for the PUT + margin. Getting this wrong costs the round,
    # not just the hedge.
    SEED_AGNOSTIC_RESERVE_S = 75.0
    # Below this the scan cannot finish, and a truncated scan is *worse* than not trying: the rows
    # it fails to hedge get backfilled from ordinary guides, and a handful of unhedged rows destroy
    # the clean-seed property the whole method rests on.
    SEED_AGNOSTIC_MIN_BUDGET_S = 150.0
    # Unhedged rows tolerated before the ordinary construction is the safer submission. Each fails
    # ~1% (Cas9) to ~4% (Cas12a) of window seeds, so a few cost a little clean fraction while dozens
    # cost all of it.
    SEED_AGNOSTIC_MAX_BACKFILL = 8
    # Ceiling on the scan even when a prepared round could afford more. Every number this hedge was
    # tuned against was measured at ~210 s, and its scans stop on "time" rather than exhaustion, so
    # a prefetch budget would silently change what the pool contains on the ~half of rounds
    # (CD34+_HSPC, HUDEP-2) that reach it. Raising this should be a measurement, not a side effect
    # of building earlier — a longer scan should only ever help, but that has not been shown.
    SEED_AGNOSTIC_MAX_BUDGET_S = 210.0

    # all-cut: the conditional hedge in genomics/all_cut.py. Tried ahead of the seed-agnostic
    # builder for HEK293, and only when its Cas12a bank is already cached — that scan is ~10 min
    # against a 300 s TTL, so a cache miss must fall through rather than miss the upload window.
    # Prebuild with: python -m niome_subnet.genomics.all_cut
    ALL_CUT = True
    # Cell types with a measured all-cut config (genomics/all_cut.CELL_CONFIG). The gain is
    # largest where the method replaces the seed-agnostic hedge outright rather than replacing an
    # earlier all-cut config: HUDEP-2 +17.92 (t=4.9 clustered, 5/5 contracts) and K562's original
    # +23.45 (t=8.4) are both against the hedge, while HEK293 gains ~+4 over its clustered builder.
    # Relaxing Cas9 from "strict over all 900 seeds" to "strict over the Cas12a-clean set" is what
    # frees the rows that were over-constrained.
    #
    # All four cell types now have a measured config, so all-cut covers every task the backend
    # issues. CD34+_HSPC is the weakest at +8.56/round — not because the method works less well
    # there (+14.27 where it builds) but because 2 of 5 contracts cannot fill stage 5's eight
    # cells at any group size that scores well. See all_cut.CELL_CONFIG for that sweep.
    ALL_CUT_CELL_TYPES = ("HEK293", "K562", "HUDEP-2", "CD34+_HSPC")
    # Least budget worth starting an all-cut build with, per cell type. The bank is built per task
    # (contract shapes never repeat, so a cache cannot be relied on) and its cost differs by an
    # order of magnitude between the two: HEK293 measured ~127 s cold, K562 271-407 s (median ~350)
    # over six cold builds at the mf100/d400 bank it moved to once the prefetch made that
    # affordable — 1.9M candidates over ~230 targets, against a few thousand at mf22/d200.
    #
    # Starting a build that cannot finish is worse than not starting one: build_bank returns []
    # on its deadline rather than a partial bank, so the whole scan is spent and all-cut declines
    # anyway, having eaten the seed-agnostic hedge's window. Per cell type rather than one number
    # because a single 400 s gate would lock HEK293 out of the ~225 s in-TTL fallback path it
    # currently completes inside.
    ALL_CUT_MIN_BUDGET_S = {"HEK293": 190.0, "K562": 480.0, "HUDEP-2": 480.0,
                            "CD34+_HSPC": 480.0}

    # --- all-HDR: the tail bet ------------------------------------------------------------------
    # Tried ahead of all-cut for the three clamped cell types. It pins the *repair mode* over a
    # narrow seed band rather than the cut over the whole window, so on a band seed every row is
    # HDR, stage 4's three targets go constant and consistency_factor is exactly 1.0.
    #
    # **This loses to all-cut on the mean and is shipped deliberately.** Measured head to head,
    # clean/dirty legs through all five stages: CD34+ 34.26 vs 56.55, K562 31.25 vs 61.12,
    # HUDEP-2 22.92 vs 41.67 — 51-61% of all-cut's expected score, about -24/round across three
    # quarters of task volume. The band is ~1.8% of the seed space, so the 95% of rounds that miss
    # it pay for the 5% that hit.
    #
    # The bet is on rank, not mean: SCORING_SYSTEM = "top" pays only the top 10 on a fixed curve,
    # so a build worth 111-194 on ~5% of rounds may out-rank a flat 56 that never spikes. No rank
    # measurement supports that here — it is the operator's call, recorded as such. Set ALL_HDR
    # False to revert to all-cut everywhere.
    ALL_HDR = True
    # HEK293 is excluded on measurement, not oversight: at accessibility 0.35 energy never clamps,
    # P(HDR) falls to ~0.37 where the guides actually sit, and the bank collapses to 1,206
    # candidates from 2.9M guides (0.041% against ~14.5%) with a 7-seed band at group 20.
    ALL_HDR_CELL_TYPES = ("CD34+_HSPC", "K562", "HUDEP-2")
    # Cheaper than all-cut despite the same target count — the band is 100 seeds, not 900.
    # Measured 31-41s for the Cas12a bank plus 59-80s for the conditional Cas9 scan.
    ALL_HDR_MIN_BUDGET_S = 190.0
    # Per-hotkey clean-band window, the decorrelation lever. all-HDR's clean band is Cas9-capped at
    # ~15 seeds and lands wherever this window is placed; a coldkey's payout is
    # 1-(1-union/900)^3, so the win comes from making sibling hotkeys' bands DISJOINT. Measured: 3
    # hotkeys on 200-299/500-599/800-899 cover 15.2% of rounds against 5.2% for the same window
    # thrice (2.9x). Set NIOME_HDR_WINDOW="200-299" (or a bare "200" for 200-299) per hotkey
    # process; unset keeps the cell-type default in all_hdr.CELL_CONFIG. Applies across all cell
    # types, since one hotkey owns one seed window regardless of the round's cell type.
    HDR_WINDOW = os.getenv("NIOME_HDR_WINDOW")

    # --- round prefetch -------------------------------------------------------------------------
    # The presigned URL's 300 s TTL bounds the *upload*, not the build — provided the rows are
    # already in hand when the validator calls. A round's task is published on /api/v3/tasks the
    # moment the round opens, hours before this miner is contacted, and the contract there is the
    # same one the validator later hands over: checked across all 41 archived rounds, contract and
    # hbb_reference match field for field, the only difference being `seed`, which is 0 on both
    # sides until the round closes. So the miner watches for its own task and builds ahead of the
    # request; ``process_task`` then only has to PUT.
    #
    # Measured over 72 rounds — task created_at against this miner's "Received genomics task" log
    # line — the lead time is: min 198 s, p10 395 s, median 1794 s, max 4404 s. 86% of rounds leave
    # 10 minutes or more, against the ~225 s the in-TTL path gets.
    PREFETCH = True
    # Public and unsigned, unlike settings.TASK_URL (/tasks/current, which 400s without the signed
    # headers a validator sends). Returns the whole task history in one page, contract and
    # hbb_reference inline, so a prepared round needs no further fetch.
    TASKS_URL = f"{settings.BASE_URL}/api/v3/tasks"
    # One round, from settings. Only used to predict when the next task is due so the poll can idle
    # in between; a schedule change costs a slower pickup, never a missed task, because
    # PREFETCH_IDLE_POLL_S keeps checking regardless.
    ROUND_SECONDS = settings.INTERVAL_BLOCKS * 12
    # Poll cadence around the expected task time, and the floor elsewhere. 15 s costs ~4 requests
    # per round on a 380 KB endpoint.
    PREFETCH_POLL_S = 15.0
    PREFETCH_IDLE_POLL_S = 300.0
    # Start polling fast this long before the next task is due, to absorb backend jitter.
    PREFETCH_LEAD_S = 180.0
    # Build budget for a prepared round. Deliberately *not* spent by anything today: all-cut
    # finishes in 30-130 s and _build_seed_agnostic keeps its own 210 s cap, so shipping the
    # prefetch changes only *when* the submission is built, never what it contains. Raising what
    # the builders do with it — a d1500 all-cut bank is ~10 min — is the next step, and this is the
    # ceiling for it. The cost of that ceiling: if a builder ever did spend the full 900 s, the 19%
    # of rounds contacted sooner than that would fall to the wait, then the emergency build below.
    PREPARE_BUDGET_S = 900.0
    # Held back from the upload deadline for the PUT itself. The upload is a ~200 KB body.
    UPLOAD_RESERVE_S = 45.0
    # And held back on top of that for one ordinary-construction build, in case a validator arrives
    # while the prepare is still running and the wait runs out. That build is ~1-5 s with the caches
    # warm; the margin is for doing it while the prepare still has the GPU and its worker pool.
    EMERGENCY_BUILD_S = 60.0
    # A failed prepare is retried rather than written off: the round still has hours left, and the
    # usual causes (a backend blip, the GPU busy, a transient CUDA error) clear on their own.
    # Spaced and capped so a contract that genuinely cannot be built does not spin all round.
    PREPARE_RETRY_S = 120.0
    PREPARE_MAX_ATTEMPTS = 3

    # The constants a contract may override. Anything absent from this tuple is global.
    TUNABLE = ("STRATEGY", "SELECTION", "GC_TOLERANCE", "CONSTRUCTION", "CUT_P_CEILING",
               "CAS_MIX", "BACKFILL_CAS", "VARIANTS", "FLANK", "LENGTHS")

    @staticmethod
    def _parse_hdr_window(raw: str | None) -> tuple[int, int] | None:
        """Parse NIOME_HDR_WINDOW into a (lo, hi) seed band, or None to keep the cell default.

        Accepts "200-299" or a bare "200" (taken as 200-299, a 100-wide window). Any malformed or
        out-of-range value returns None with a warning: a bad env var must not stop the miner, only
        cost this hotkey its decorrelation.
        """
        if not raw:
            return None
        try:
            if "-" in raw:
                lo, hi = (int(x) for x in raw.split("-", 1))
            else:
                lo = int(raw)
                hi = lo + 99
            if not (100 <= lo < hi <= 999):
                raise ValueError(f"window {lo}-{hi} outside 100-999 or non-increasing")
            return lo, hi
        except (ValueError, TypeError) as exc:
            logger.warning(f"NIOME_HDR_WINDOW={raw!r} is not a valid seed window ({exc}); "
                           "using the cell-type default band")
            return None

    def __init__(self, config=None):
        super(Miner, self).__init__(config=config)

        self.gen_config = self.base_gen_config()

        # One build at a time, and at most one per task. Every validator broadcasts the same task id
        # with its own presigned URL, and the rows are a deterministic function of the contract, so
        # later broadcasts should reuse the first build and repeat only the upload.
        self._build_lock = asyncio.Lock()
        self._built: tuple[str, list[dict]] | None = None
        # asyncio holds only a weak reference to a running task, so a fire-and-forget
        # ``create_task`` result can be collected mid-flight and the round vanishes with no log
        # line at all. Keeping the handle here until it completes is what makes the submission
        # survive — and what lets the done-callback below report a crash that nothing awaits.
        self._inflight: set[asyncio.Task] = set()
        # Every validator broadcasts the same task, so process_task can run several times per
        # round. One archive/scoring pass at a time, and the losers skip rather than queue.
        self._archive_lock = threading.Lock()
        # The round the prefetch thread is working on, or has finished. Written by that thread and
        # read by process_task on the event loop; a single attribute rebind is atomic, and every
        # field of a PreparedRound is either set before publication or guarded by ``done``.
        self._prepared: PreparedRound | None = None
        # Serialises the hedge builders across the prefetch thread and process_task's worker
        # thread. They own the GPU and a worker pool, so two at once is slower than either alone —
        # and _build_lock cannot cover it, being an asyncio lock the prefetch thread cannot take.
        self._hedge_lock = threading.Lock()
        # chr11 is loaded through an unguarded module global, so the prefetch waits for the prewarm
        # rather than racing it into a second 130 MB read.
        self._prewarmed = threading.Event()

        # This hotkey's all-HDR clean-band window, parsed once. A bad value logs and falls back to
        # the cell-type default rather than taking the process down — a misconfigured window should
        # cost decorrelation, not the miner.
        self.hdr_window = self._parse_hdr_window(self.HDR_WINDOW)
        if self.hdr_window:
            logger.info(f"all-HDR clean-band window pinned to {self.hdr_window[0]}-"
                        f"{self.hdr_window[1]} for this hotkey (NIOME_HDR_WINDOW)")

        logger.info(
            f"Generation config: strategy={self.STRATEGY} selection={self.SELECTION} "
            f"construction={self.CONSTRUCTION} variants={self.VARIANTS} flank={self.FLANK} "
            f"lengths={self.LENGTHS} score_locally={self.SCORE_LOCALLY}"
        )
        threading.Thread(target=self._prewarm, name="niome-prewarm", daemon=True).start()
        if self.PREFETCH:
            threading.Thread(target=self._prefetch_loop, name="niome-prefetch",
                             daemon=True).start()
        else:
            logger.info("Prefetch disabled; submissions will be built inside the upload TTL")

    @classmethod
    def settings_for(cls, contract: dict) -> dict:
        """The generation constants that apply to this contract, after any cell-type override."""
        overrides = cls.CELL_TYPE_OVERRIDES.get(contract.get("cell_type"), {})
        return {name: overrides.get(name, getattr(cls, name)) for name in cls.TUNABLE}

    @classmethod
    def gen_config_for(cls, contract: dict) -> G.GenConfig:
        """The fully resolved config for one contract: constants, overrides, then the knobs that
        can only be keyed once ``rules.cas_systems`` is known.

        The single entry point for anything that needs to reproduce what the miner will build —
        ``_build`` here, and the offline harnesses — so a cell-type override cannot apply on one
        side and not the other.
        """
        chosen = cls.settings_for(contract)
        cfg = G.GenConfig(
            strategy=chosen["STRATEGY"],
            selection=chosen["SELECTION"],
            construction=chosen["CONSTRUCTION"],
            gc_tolerance=chosen["GC_TOLERANCE"],
            backfill_cas=chosen["BACKFILL_CAS"],
            variants=chosen["VARIANTS"],
            flank=chosen["FLANK"],
            lengths=tuple(chosen["LENGTHS"]),
        )
        return G.config_for_contract(cfg, contract, cut_p_ceiling=chosen["CUT_P_CEILING"],
                                     cas_mix=chosen["CAS_MIX"])

    @classmethod
    def base_gen_config(cls) -> G.GenConfig:
        """The static half of the generation config, built from the class constants.

        The contract-dependent knobs — the per-Cas cut_p floors and the cas mix — are keyed by Cas
        system name, and that roster comes from ``rules.cas_systems``, so neither can be resolved
        until a task arrives; ``_build`` fills them in through ``G.config_for_contract``, the same
        call submission.py makes. Everything the prewarm needs (flank, lengths) is already here.

        Exposed as a classmethod so an offline harness can reproduce the miner's exact
        configuration without standing up a wallet and a chain connection.
        """
        return G.GenConfig(
            strategy=cls.STRATEGY,
            selection=cls.SELECTION,
            construction=cls.CONSTRUCTION,
            gc_tolerance=cls.GC_TOLERANCE,
            backfill_cas=cls.BACKFILL_CAS,
            variants=cls.VARIANTS,
            flank=cls.FLANK,
            lengths=tuple(cls.LENGTHS),
        )

    async def forward(self, body: bytes, caller_hotkey: str) -> dict:
        """
        Processes an incoming genomics task request.

        Args:
            body: Raw JSON body bytes from the validator.
            caller_hotkey: Verified hotkey ss58 of the calling validator.

        Returns:
            dict: Response payload (empty acknowledgement).
        """
        try:
            synapse = GenomicsTaskSynapse.model_validate_json(body)
            if synapse.task is None:
                logger.error(f"Task missing from {caller_hotkey}'s request; ignoring")
                return {}

            task = synapse.task
            logger.info(f"Received genomics task {task.id} from {caller_hotkey}")
            logger.info(
                f"[{task.id}] contract_url={task.contract_url} hbb_ref_url={task.hbb_ref_url}"
            )
            if not synapse.presigned_url:
                logger.error(
                    f"[{task.id}] no presigned URL in the request — there is nowhere to upload "
                    "the submission, so this round cannot be scored"
                )
                return {}
            logger.info(f"[{task.id}] upload target {self._url_summary(synapse.presigned_url)}")

            # Fire and forget - run process_task asynchronously without waiting
            handle = asyncio.create_task(self.process_task(task, synapse.presigned_url))
            self._inflight.add(handle)
            handle.add_done_callback(self._on_task_done)
            logger.info(f"[{task.id}] build+upload scheduled in the background; acking validator")

            return {}
        except Exception as e:
            logger.error(f"Forward error: {e}")
            logger.debug(traceback.format_exc())
            return {"error": str(e)}

    def _on_task_done(self, handle: asyncio.Task) -> None:
        """Drop the strong reference and surface anything ``process_task`` failed to catch."""
        self._inflight.discard(handle)
        if handle.cancelled():
            logger.error("Background submission task was cancelled before it finished")
            return
        error = handle.exception()
        if error is not None:
            logger.error(f"Background submission task crashed: {error!r}")
            logger.debug(
                "".join(traceback.format_exception(type(error), error, error.__traceback__))
            )

    async def process_task(self, task: Task, presigned_url: str) -> None:
        """Build this task's dataset and upload it to the validator's presigned URL.

        Called fire-and-forget from ``forward``, so it must neither raise into the request path nor
        block the event loop: the build and the HTTP calls are synchronous and CPU-bound, and go to
        worker threads so ``/forward`` stays answerable while a build is in flight.
        """
        if task is None or not presigned_url:
            logger.error("No task or no presigned URL — nothing to submit")
            return

        started = time.time()
        deadline = self._upload_deadline(presigned_url, started)
        tag = f"[{task.id}]"
        logger.info(
            f"{tag} step 1/5 starting submission; upload deadline in "
            f"{deadline - started:.0f}s"
        )
        # The signed URL is the only way back to this key, and the log deliberately prints it
        # without its signature. Recording it here is what makes scripts/resubmit.py possible at
        # all: a failed upload can be retried by hand while the TTL lasts, which is the only
        # recovery there is — the validator never asks twice within a task id.
        self._record_upload(task.id, presigned_url, deadline, submitted=False)

        try:
            logger.info(f"{tag} step 2/5 fetching contract and HBB reference")
            fetch_started = time.time()
            contract, reference = await asyncio.to_thread(self._fetch_artifacts, task)
            logger.info(
                f"{tag} step 2/5 done in {time.time() - fetch_started:.1f}s "
                f"| seed={contract.get('seed')} cell_type={contract.get('cell_type')} "
                f"mutations={len(contract.get('active_mutations') or [])} "
                f"rules={contract.get('rules')}"
            )

            logger.info(f"{tag} step 3/5 fetching cell-type accessibility table")
            cell_types = await asyncio.to_thread(self._fetch_cell_types)
            logger.info(f"{tag} step 3/5 done | {len(cell_types)} cell types")

            logger.info(f"{tag} step 4/5 obtaining the dataset")
            build_started = time.time()
            rows = await self._rows_for_task(task, contract, reference, cell_types, deadline, tag)
            logger.info(
                f"{tag} step 4/5 done in {time.time() - build_started:.1f}s | {len(rows)} rows"
            )

            if not rows:
                logger.error(f"Built no rows for task {task.id} — nothing to upload")
                return

            logger.info(
                f"{tag} step 5/5 uploading {len(rows)} rows, "
                f"{deadline - time.time():.0f}s of TTL left"
            )
            await asyncio.to_thread(self._upload, presigned_url, rows, deadline)
            self._record_upload(task.id, presigned_url, deadline, submitted=True, rows=len(rows))
            logger.info(
                f"Submitted {len(rows)} rows for task {task.id} in {time.time() - started:.1f}s "
                f"({deadline - time.time():.0f}s of the URL's TTL to spare)"
            )

            # Housekeeping, not submission: score the previous round now that its seed exists,
            # then archive this one. The PUT has landed and the URL is spent, so nothing below
            # here can cost the round — but it is minutes of CPU, so it goes to a worker thread.
            #
            # Re-persist first: _build writes this file, but a prepared round and a fallback build
            # can both have written it since, and the archive must hold what was actually sent.
            self._persist(settings.MINER_SUBMISSION_PATH, rows)
            await asyncio.to_thread(self._post_upload, task.id)
        except Exception as e:
            # Nothing downstream reports a miner failure — a missed upload is indistinguishable from
            # never having been contacted, and there is no retry within a task id. So this log line
            # is the only evidence the round was lost.
            logger.error(f"Failed to submit task {getattr(task, 'id', '?')}: {e}")
            logger.debug(traceback.format_exc())
            logger.error(
                f"{deadline - time.time():.0f}s of the URL's TTL remain — retry by hand with "
                "`python scripts/resubmit.py` while that is positive"
            )

    def _record_upload(
        self,
        task_id: str,
        presigned_url: str,
        deadline: float,
        submitted: bool,
        rows: int | None = None,
    ) -> None:
        """Note the current task's upload target and its outcome, for scripts/resubmit.py.

        Best-effort: a failure to write this must not be what loses the round, so it is logged and
        swallowed. The file holds a signed URL, which is a write capability on one bucket key until
        it expires — ``data/`` is gitignored, and it is worthless a few minutes later.
        """
        try:
            self._persist(
                settings.LAST_UPLOAD_PATH,
                {
                    "task_id": task_id,
                    "presigned_url": presigned_url,
                    "expires_at": deadline,
                    "submitted": submitted,
                    "rows": rows,
                },
            )
        except Exception as e:
            logger.warning(f"Could not record the upload target ({e}); manual retry will need "
                           "the URL from the validator")

    def _post_upload(self, task_id: str) -> None:
        """Score the previous round, then archive this one. Never raises.

        The order is the point. When this runs, ``data/result/<latest>`` is the *previous* round:
        its task has closed, so the backend has stamped its seed and calc.py can score it against
        the stream the validator actually used. This round's folder is written afterwards, which is
        what makes it the next run's "previous round". Scoring a folder before its task closes would
        only measure it at seed 0.
        """
        if not self._archive_lock.acquire(blocking=False):
            logger.info("post-upload: an archive pass is already running; skipping this one")
            return
        try:
            self._score_previous_round(task_id)
            self._archive_round(task_id)
        except Exception as e:
            # A failed archive costs nothing that was not already uploaded, so it is logged and
            # dropped rather than raised into process_task's handler, which reports lost rounds.
            logger.warning(f"post-upload housekeeping failed: {e}")
            logger.debug(traceback.format_exc())
        finally:
            self._archive_lock.release()

    def _latest_archive(self) -> Path | None:
        """The newest round folder. Names are timestamps, so lexicographic order is chronological."""
        root = Path(self.ARCHIVE_ROOT)
        if not root.is_dir():
            return None
        folders = sorted((path for path in root.iterdir() if path.is_dir()), key=lambda p: p.name)
        return folders[-1] if folders else None

    def _score_previous_round(self, task_id: str) -> None:
        """Run calc.py over the newest archive folder, unless it has already been scored.

        calc.py is invoked as a subprocess, pinned to the task id the folder recorded when it was
        archived: it fetches that task from ``/api/v3/tasks`` for the stamped seed, re-runs all five
        validator stages and writes the per-stage detail back into the folder.
        """
        folder = self._latest_archive()
        if folder is None:
            logger.info(f"post-upload: nothing archived under {self.ARCHIVE_ROOT} yet, so there "
                        "is nothing to score")
            return
        if (folder / self.VALIDATION_MARKER).exists():
            logger.info(f"post-upload: {folder} already holds {self.VALIDATION_MARKER}; "
                        "not scoring it again")
            return

        submission = folder / Path(settings.MINER_SUBMISSION_PATH).name
        if not submission.exists():
            logger.warning(f"post-upload: {folder} has no {submission.name} to score")
            return

        recorded = self._read_json(folder / Path(settings.LAST_UPLOAD_PATH).name) or {}
        archived_task = recorded.get("task_id")
        if archived_task and archived_task == task_id:
            # This round's own folder, from an earlier broadcast of the same task or a restart.
            # Its seed is stamped only after the task closes, so scoring it now would measure the
            # submission against a stream no validator will use. It becomes scorable next round.
            logger.info(f"post-upload: {folder} is task {task_id}'s own archive — its seed is not "
                        "stamped yet, so scoring waits for the next round")
            return

        command = [sys.executable, os.path.join(PROJECT_ROOT, self.VALIDATION_SCRIPT),
                   "--folder", str(folder), "--quiet"]
        if archived_task:
            command += ["--task-id", archived_task]
        else:
            # No recorded id: calc.py falls back to the folder name against the task's created_at,
            # then to the contract fingerprint, and refuses outright if what it finds disagrees
            # with the archived contract.
            logger.info(f"post-upload: {folder} records no task id; letting calc.py match it")

        log_path = folder / "validation.log"
        logger.info(f"post-upload: scoring {folder} with {self.VALIDATION_SCRIPT}"
                    + (f" against task {archived_task}" if archived_task else "")
                    + f"; output -> {log_path}")
        started = time.time()
        try:
            with open(log_path, "w") as handle:
                completed = subprocess.run(
                    command,
                    cwd=PROJECT_ROOT,          # every settings path is relative to the repo root
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                    timeout=self.VALIDATION_TIMEOUT_S,
                    check=False,
                )
        except subprocess.TimeoutExpired:
            logger.error(f"post-upload: {self.VALIDATION_SCRIPT} exceeded "
                         f"{self.VALIDATION_TIMEOUT_S}s on {folder}; see {log_path}")
            return

        elapsed = time.time() - started
        if completed.returncode != 0:
            logger.error(f"post-upload: {self.VALIDATION_SCRIPT} exited {completed.returncode} "
                         f"after {elapsed:.0f}s — see {log_path}")
            return

        summary = self._read_json(folder / self.VALIDATION_MARKER) or {}
        score = summary.get("score") or {}
        rows = summary.get("rows") or {}
        logger.info(
            f"post-upload: scored {folder} in {elapsed:.0f}s | task {summary.get('task_id')} "
            f"seed {summary.get('scored_seed')} | {rows.get('valid')}/{rows.get('scored')} rows "
            f"valid | final_score {self._number(score.get('final_score'))} = weighted "
            f"{self._number(score.get('total_weighted_score'))} x consistency "
            f"{self._number(score.get('consistency_factor'))} x fidelity "
            f"{self._number(score.get('distribution_fidelity_factor'))}"
        )

    def _archive_round(self, task_id: str) -> Path | None:
        """Copy this round's artifacts into a new ``data/result/<archived at>`` folder.

        Everything needed to re-score the upload later, and nothing else: the rows exactly as sent,
        the broadcast contract and reference they were designed against, and the upload record that
        names the task. Skipped when the newest folder already records this task id, so repeated
        broadcasts of one task archive once.
        """
        latest = self._latest_archive()
        if latest is not None:
            recorded = self._read_json(latest / Path(settings.LAST_UPLOAD_PATH).name) or {}
            if recorded.get("task_id") == task_id:
                logger.info(f"post-upload: task {task_id} is already archived in {latest}")
                return None

        stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S")
        folder = Path(self.ARCHIVE_ROOT) / stamp
        folder.mkdir(parents=True, exist_ok=True)

        copied, missing = [], []
        for name in self.ARCHIVED_PATHS:
            source = Path(getattr(settings, name))
            if not source.exists():
                missing.append(str(source))
                continue
            shutil.copyfile(source, folder / source.name)
            copied.append(source.name)

        if copied:
            logger.info(f"post-upload: archived {', '.join(copied)} for task {task_id} -> {folder}")
        if missing:
            # A missing last_upload.json costs the folder its task pairing; a missing submission
            # makes it unscorable altogether. Either way the round itself is already uploaded.
            logger.warning(f"post-upload: not in data/, so not archived: {', '.join(missing)}")
        return folder

    @staticmethod
    def _read_json(path: Path | str):
        """Read a JSON file, or None if it is absent or unreadable. Used on archive metadata,
        where a half-written or hand-edited file must not take the miner down."""
        try:
            with open(path) as handle:
                return json.load(handle)
        except (OSError, json.JSONDecodeError):
            return None

    @staticmethod
    def _number(value, digits: int = 4) -> str:
        """Format a score that calc.py may not have produced."""
        return f"{value:.{digits}f}" if isinstance(value, (int, float)) else "?"

    def _prewarm(self) -> None:
        """Load the reference — and, if a previous task left its artifacts behind, the k-mer index
        and PAM sites too — before the first task arrives.

        All three caches are task-independent (every task issued so far shares one gene_region and
        one rules block), so this is pure critical-path removal: a warm process only pays for the
        build itself.
        """
        started = time.time()
        try:
            logger.info("Prewarm: loading chr11 reference sequence")
            self._load_sequence()
            logger.info(f"Prewarm: chr11 loaded in {time.time() - started:.1f}s")
            if os.path.exists(settings.CONTRACT_PATH) and os.path.exists(
                settings.HBB_REFERENCE_PATH
            ):
                logger.info("Prewarm: replaying the last task's artifacts to warm the site cache")
                with open(settings.CONTRACT_PATH) as handle:
                    contract = json.load(handle)
                with open(settings.HBB_REFERENCE_PATH) as handle:
                    reference = json.load(handle)
                # cell_types is left empty on purpose: it is carried on the Context but nothing
                # cached here reads it, so warming does not need the backend to answer.
                context = G.build_context(contract, reference, {})
                sites = G.enumerate_sites(context, self.gen_config.flank, self._lengths())
                logger.info(
                    f"Prewarmed chr11 and {len(sites)} PAM sites in {time.time() - started:.1f}s"
                )
            else:
                logger.info(
                    f"Prewarmed chr11 in {time.time() - started:.1f}s; PAM sites will be "
                    "enumerated on the first task"
                )
        except Exception as e:
            logger.error(f"Prewarm failed ({e}); the first task will pay the cold start")
            logger.debug(traceback.format_exc())
        finally:
            # Released on failure too: the prefetch has its own error handling, and blocking it
            # forever on a bad prewarm would silently disable the whole path.
            self._prewarmed.set()

    # ---- Round prefetch ------------------------------------------------------------------------

    def _prefetch_loop(self) -> None:
        """Watch /api/v3/tasks for each round's task and build its submission on sight.

        Runs for the life of the process in its own thread. Every failure is caught and retried on
        the next poll: this path is an optimisation, and ``process_task`` still builds in-TTL when
        it finds nothing prepared, so a backend outage must cost latency and not the round.
        """
        logger.info(
            f"Prefetch: watching {self.TASKS_URL} for seed-0 tasks; build budget "
            f"{self.PREPARE_BUDGET_S:.0f}s (against ~{settings.SUBMISSION_TIMEOUT}s in-TTL)"
        )
        self._prewarmed.wait()
        while not self.should_exit:
            try:
                item = self._newest_unstamped_task()
                prepared = self._prepared
                if item is None:
                    logger.debug("Prefetch: newest task already carries a seed; nothing to prepare")
                elif prepared is None or prepared.task_id != item["id"]:
                    self._prepare(item)
                elif (prepared.done.is_set() and not prepared.rows
                        and prepared.attempt < self.PREPARE_MAX_ATTEMPTS
                        and time.time() - prepared.started >= self.PREPARE_RETRY_S):
                    logger.info(
                        f"Prefetch: retrying {item['id']} (attempt {prepared.attempt + 1} of "
                        f"{self.PREPARE_MAX_ATTEMPTS}) after: {prepared.error}"
                    )
                    self._prepare(item, attempt=prepared.attempt + 1)
            except Exception as e:
                logger.warning(f"Prefetch poll failed ({e}); retrying")
                logger.debug(traceback.format_exc())
            self._sleep(self._next_poll_delay())

    def _sleep(self, seconds: float) -> None:
        """Sleep in slices so shutdown stays responsive across a five-minute idle poll."""
        until = time.time() + seconds
        while not self.should_exit and time.time() < until:
            time.sleep(min(5.0, until - time.time()))

    def _next_poll_delay(self) -> float:
        """How long until the next poll: fast around the time a task is due, slow in between.

        The due time is the current task's ``created_at`` plus one round, which is where the
        backend has put every task so far. It is only a hint — an idle poll runs regardless, so a
        schedule that moves costs a slower pickup rather than a missed round.
        """
        prepared = self._prepared
        if prepared is None or not prepared.created_at:
            return self.PREFETCH_POLL_S
        if (prepared.done.is_set() and not prepared.rows
                and prepared.attempt < self.PREPARE_MAX_ATTEMPTS):
            # A retry is outstanding. Come back for it on the retry interval instead of idling
            # until the next round is due, which would delay it by up to PREFETCH_IDLE_POLL_S.
            return max(self.PREFETCH_POLL_S,
                       self.PREPARE_RETRY_S - (time.time() - prepared.started))
        try:
            created = datetime.fromisoformat(prepared.created_at).replace(tzinfo=timezone.utc)
        except ValueError:
            return self.PREFETCH_IDLE_POLL_S
        due = created.timestamp() + self.ROUND_SECONDS
        until_window = due - self.PREFETCH_LEAD_S - time.time()
        if until_window <= 0:
            return self.PREFETCH_POLL_S             # inside the window the next task is due in
        return min(self.PREFETCH_IDLE_POLL_S, until_window)

    def _newest_unstamped_task(self) -> dict | None:
        """The current round's task, or None when the newest one has already been stamped.

        Only the newest item is considered. The backend stamps a task's seed when its round closes,
        so exactly one unstamped task exists at a time and it is always the newest — an older
        unstamped one would be a task whose round is over and whose upload window is long gone.
        """
        response = requests.get(self.TASKS_URL, timeout=settings.TASK_REQUEST_TIMEOUT)
        response.raise_for_status()
        items = response.json().get("items") or []
        if not items:
            raise RuntimeError("backend returned no tasks")
        newest = max(items, key=lambda item: item.get("created_at") or "")
        contract = (newest.get("content") or {}).get("contract") or {}
        return newest if self._is_unstamped(contract.get("seed")) else None

    @staticmethod
    def _is_unstamped(seed) -> bool:
        """Whether a contract seed is the broadcast placeholder rather than a scoring seed.

        ``seed`` is a comma-joined list of round seeds once stamped ("130,507,441"); before that it
        is 0. Parsed rather than compared literally because the placeholder has been seen as int 0,
        "0" and absent, and a build keyed to the wrong one of those is a lost round either way.
        """
        if seed is None or seed == "":
            return True
        try:
            return all(int(part) == 0 for part in str(seed).split(","))
        except ValueError:
            return False

    def _prepare(self, item: dict, attempt: int = 1) -> None:
        """Build the submission for a task that no validator has asked for yet.

        Publishes the round before building, so a validator arriving mid-build waits on this rather
        than starting a second one.
        """
        task_id = item["id"]
        content = item.get("content") or {}
        contract, reference = content.get("contract"), content.get("hbb_reference")
        if not contract or not reference:
            logger.warning(f"Prefetch: task {task_id} carries no contract/reference; skipping")
            return

        prepared = PreparedRound(
            task_id=task_id,
            key=self._build_key(task_id, contract),
            contract=contract,
            reference=reference,
            created_at=item.get("created_at") or "",
            attempt=attempt,
        )
        self._prepared = prepared
        logger.info(
            f"Prefetch: preparing task {task_id} (created {prepared.created_at}) | "
            f"cell_type={contract.get('cell_type')} "
            f"mutations={len(contract.get('active_mutations') or [])} "
            f"budget={self.PREPARE_BUDGET_S:.0f}s"
            + (f" | attempt {attempt}" if attempt > 1 else "")
        )
        try:
            # Same files _fetch_artifacts writes, so a restart's prewarm and any offline re-score
            # see this round's artifacts whether or not a validator has been in touch yet.
            self._persist(settings.CONTRACT_PATH, contract)
            self._persist(settings.HBB_REFERENCE_PATH, reference)
            cell_types = self._fetch_cell_types()
            prepared.rows = self._build(contract, reference, cell_types,
                                        budget_s=self.PREPARE_BUDGET_S)
        except Exception as e:
            prepared.error = str(e)
            logger.error(f"Prefetch: build for {task_id} failed ({e}); the validator's request "
                         "will fall back to an in-TTL build")
            logger.debug(traceback.format_exc())
        finally:
            prepared.elapsed = time.time() - prepared.started
            prepared.done.set()
        if prepared.rows:
            logger.info(
                f"Prefetch: task {task_id} ready — {len(prepared.rows)} rows in "
                f"{prepared.elapsed:.0f}s, waiting for a validator to ask for it"
            )

    async def _rows_for_task(self, task: Task, contract: dict, reference: dict, cell_types: dict,
                             deadline: float, tag: str) -> list[dict]:
        """The rows to upload for this task: the prepared ones where possible, built ones if not.

        Three paths, in descending order of how much time the build got:
        prepared (hours), prepared-but-still-running (wait for it), and in-TTL (what the miner did
        before the prefetch existed, and still the fallback whenever the other two miss).
        """
        key = self._build_key(task.id, contract)
        prepared = self._prepared
        if prepared is not None and prepared.key == key:
            if not prepared.done.is_set():
                # Leave room for the PUT and for one ordinary build, so a prepare that overruns
                # costs quality rather than the round.
                patience = deadline - self.UPLOAD_RESERVE_S - self.EMERGENCY_BUILD_S - time.time()
                logger.info(
                    f"{tag} prepared build is still running ({time.time() - prepared.started:.0f}s "
                    f"in); waiting up to {patience:.0f}s for it"
                )
                if patience > 0:
                    await asyncio.to_thread(prepared.done.wait, patience)
            # ``done`` is the only thing that makes ``rows`` trustworthy: it is set in the
            # builder's finally, after the rows are complete. Reading rows without it would be
            # right only by accident of assignment order.
            if prepared.done.is_set() and prepared.rows:
                logger.info(
                    f"{tag} using the prepared submission ({prepared.summary()}) — no build needed"
                )
                return prepared.rows
            logger.warning(f"{tag} prepared round unusable ({prepared.summary()}); building now")
        elif prepared is not None:
            # Same round, different contract, or a round we never saw. Either way the prepared rows
            # were designed against different rules and must not be sent.
            logger.warning(
                f"{tag} prepared round is for task {prepared.task_id} under a different contract "
                f"key; building this one instead"
            )
        else:
            logger.info(f"{tag} nothing prepared for this task; building inside the upload TTL")

        # Hedges are skipped when a prepare is still holding the GPU and its worker pool: the
        # ordinary construction is the build that reliably finishes in what is left of the TTL.
        still_preparing = prepared is not None and not prepared.done.is_set()
        async with self._build_lock:
            if self._built is not None and self._built[0] == key:
                logger.info(f"{tag} reusing the {len(self._built[1])}-row build for this task")
                return self._built[1]
            rows = await asyncio.to_thread(
                self._build, contract, reference, cell_types, deadline,
                None, not still_preparing,
            )
            # Only the current task's rows are worth keeping: a new task means a new contract, and
            # the old rows can never be submitted again.
            self._built = (key, rows)
        return rows

    @contextlib.contextmanager
    def _hedge_slot(self, timeout: float):
        """Serialise the hedge builders, yielding False rather than queueing past ``timeout``.

        A prepared round can hold this for the whole of ``PREPARE_BUDGET_S``. An in-TTL build that
        cannot get in must say so and use the ordinary construction: waiting its turn would spend
        the upload window, and a late submission scores nothing at all.
        """
        acquired = self._hedge_lock.acquire(timeout=timeout)
        try:
            yield acquired
        finally:
            if acquired:
                self._hedge_lock.release()

    def _seed_agnostic_applies(self, contract: dict, budget_s: float) -> tuple[bool, str]:
        """Whether the bank builder should run for this contract, and why not when it should not."""
        if not self.SEED_AGNOSTIC:
            return False, "disabled by SEED_AGNOSTIC"
        if contract.get("cell_type") in self.SEED_AGNOSTIC_SKIP_CELL_TYPES:
            return False, f"{contract.get('cell_type')} has its own clustered builder"
        if contract.get("seed"):
            # A real seed makes the ordinary construction exact — consistency_factor 1.0 and a
            # score around 290, against ~41 for a hedge against an unknown seed. Never trade that
            # away; the hedge exists only for the case the seed is not known at build time.
            return False, f"contract carries seed {contract['seed']}; the matched build scores far higher"
        if SA.MT._gpu() is None:
            # CPU-only the scan takes ~540s against a 300s TTL, so it would always truncate.
            return False, "no usable GPU; the scan cannot finish inside the upload window"
        if budget_s < self.SEED_AGNOSTIC_MIN_BUDGET_S:
            return False, f"only {budget_s:.0f}s of budget; need {self.SEED_AGNOSTIC_MIN_BUDGET_S:.0f}s"
        return True, ""

    def _build_seed_agnostic(self, contract: dict, reference: dict, cell_types: dict,
                             budget_s: float) -> list[dict] | None:
        """Build from the strict-Cas9 / min-union-Cas12a banks. None means fall back."""
        task = {"content": {"contract": contract, "hbb_reference": reference}}
        budget_s = min(budget_s, self.SEED_AGNOSTIC_MAX_BUDGET_S)
        cfg = replace(
            SA.SeedAgnosticConfig(),
            time_budget_s=budget_s,
            max_experiments=contract["rules"].get("max_experiments"),
        )
        logger.info(
            f"Build: seed-agnostic bank scan | window {cfg.start_seed}-{cfg.end_seed} "
            f"cas_mix={cfg.cas_mix} budget={budget_s:.0f}s gpu=on"
        )
        started = time.time()
        rows, meta = SA.build_submission(contract, reference, cell_types, task, cfg)
        logger.info(
            f"Build: seed-agnostic done in {time.time() - started:.0f}s | rows={meta['rows']} "
            f"hedged={meta['hedged_rows']} backfilled={meta['backfilled']} "
            f"strict_cas9={meta['strict_cas9_available']} cas12a_pool={meta['cas12a_pool']} "
            f"group={meta['cas12a_group']} union={meta['group_failed_seed_union']} "
            f"clean~{meta['clean_fraction_estimate']:.1%} cas={meta['cas_mix']}"
        )

        expected = cfg.max_experiments or 0
        if expected and meta["rows"] < expected:
            logger.warning(
                f"Build: seed-agnostic produced {meta['rows']}/{expected} rows; "
                "falling back to the ordinary construction"
            )
            return None
        if meta["backfilled"] > self.SEED_AGNOSTIC_MAX_BACKFILL:
            # The hedge is all-or-nothing: measured, 66 unhedged rows took the clean fraction to 0%
            # and the score below the ordinary build. Prefer the build we can reason about.
            logger.warning(
                f"Build: {meta['backfilled']} rows outside the hedge exceeds "
                f"{self.SEED_AGNOSTIC_MAX_BACKFILL}; falling back to the ordinary construction"
            )
            return None
        return rows

    def _build(self, contract: dict, reference: dict, cell_types: dict,
               deadline: float | None = None, budget_s: float | None = None,
               allow_hedges: bool = True) -> list[dict]:
        """Generate this task's submission and log what it should be worth.

        Most cell types follow the same sequence as ``submission.build_for_task``. HEK293 instead
        uses the dedicated clustered builder because its provisional seed can be replaced before
        scoring; that builder does not inspect simulated outcomes while choosing designs.

        ``budget_s`` states the time the hedge builders may take outright, for a prepared round
        that is not racing an upload TTL; without it the budget is derived from ``deadline`` as
        before. ``allow_hedges=False`` skips them entirely — the emergency path, for when a
        prepare already holds the GPU and only the ordinary construction will finish in time.
        """
        logger.info("Build: loading chr11 (cached after the first call)")
        self._load_sequence()

        def remaining() -> float:
            if budget_s is not None:
                return budget_s
            if deadline:
                return deadline - time.time() - self.SEED_AGNOSTIC_RESERVE_S
            return float(self.SEED_AGNOSTIC_MIN_BUDGET_S)

        budget = remaining()
        # How long to queue for the hedge slot. A prepared round can afford to wait out a
        # straggler; a build racing the upload TTL cannot, and takes the ordinary construction
        # instead of spending its window in a queue.
        hedge_wait = 60.0 if budget_s is not None else 10.0
        if not allow_hedges:
            logger.info("Build: hedges skipped (a prepared build is still holding the GPU); "
                        "using the ordinary construction")

        # All-cut goes first, ahead of the seed-agnostic hedge: it *is* that hedge with the Cas9
        # constraint relaxed from "strict over all 900 seeds" to "strict over the Cas12a-clean
        # set", and it measured +23.45 (t=8.4) over it on K562. Placed after, it was dead code —
        # _build_seed_agnostic returns, so K562 never reached it.
        # All-HDR goes ahead of all-cut for the cell types it is configured for. It declines to
        # None on an unmeasured cell type or a short pool, and all-cut below is then the fallback —
        # so HEK293, and any failure on the other three, still gets the build it had before.
        cell_type = contract.get("cell_type")
        if (allow_hedges and self.ALL_HDR and cell_type in self.ALL_HDR_CELL_TYPES
                and budget >= self.ALL_HDR_MIN_BUDGET_S):
            try:
                with self._hedge_slot(hedge_wait) as slot:
                    if slot:
                        hdr_rows, hdr_meta = AH.build_for_cell(
                            contract, reference, cell_types, budget_s=budget,
                            hdr_range=self.hdr_window)
                    else:
                        hdr_rows, hdr_meta = None, {"reason": "another build holds the hedge slot"}
            except Exception as exc:
                logger.warning(f"Build: all-HDR failed ({exc}); falling through to all-cut")
                logger.debug(traceback.format_exc())
                hdr_rows, hdr_meta = None, {"reason": str(exc)}
            if hdr_rows:
                logger.info(
                    f"Build: all-HDR ({cell_type}) | band {hdr_meta['band']} "
                    f"group {hdr_meta['group_size']} "
                    f"clean {hdr_meta['clean']}/{hdr_meta['clean'] + hdr_meta['union']} "
                    f"| cas9 pool {hdr_meta['cas9_pool']} | rows {hdr_meta['rows']} "
                    f"cells {hdr_meta['cells']}/8 | {hdr_meta['elapsed_s']}s"
                )
                self._persist(settings.MINER_SUBMISSION_PATH, hdr_rows)
                return hdr_rows
            logger.info("Build: all-HDR declined (%s); falling through to all-cut",
                        hdr_meta.get("reason", "unknown"))
            budget = remaining()

        min_budget = self.ALL_CUT_MIN_BUDGET_S.get(contract.get("cell_type"), 0.0)
        if (allow_hedges and self.ALL_CUT
                and contract.get("cell_type") in self.ALL_CUT_CELL_TYPES
                and budget >= min_budget):
            try:
                with self._hedge_slot(hedge_wait) as slot:
                    if slot:
                        allcut_rows, allcut_meta = AC.build_for_cell(
                            contract, reference, cell_types, budget_s=budget,
                            reserve_s=float(self.SEED_AGNOSTIC_MIN_BUDGET_S))
                    else:
                        allcut_rows = None
                        allcut_meta = {"reason": "another build holds the hedge slot"}
            except Exception as exc:
                logger.warning(f"Build: all-cut failed ({exc}); falling through")
                logger.debug(traceback.format_exc())
                allcut_rows, allcut_meta = None, {"reason": str(exc)}
            if allcut_rows:
                logger.info(
                    f"Build: all-cut ({contract.get('cell_type')}) | "
                    f"group {allcut_meta['group_size']} "
                    f"clean {allcut_meta['clean']}/900 ({100*allcut_meta['clean_fraction']:.1f}%) "
                    f"| cas9 pool {allcut_meta['cas9_pool']} | rows {allcut_meta['rows']} "
                    f"cells {allcut_meta['cells']}/8 | {allcut_meta['elapsed_s']}s"
                    + ("  (retried)" if allcut_meta.get("retried") else "")
                )
                self._persist(settings.MINER_SUBMISSION_PATH, allcut_rows)
                return allcut_rows
            logger.info("Build: all-cut declined (%s); falling through",
                        allcut_meta.get("reason", "unknown"))
            budget = remaining()

        applies, why_not = self._seed_agnostic_applies(contract, budget) if allow_hedges else (
            False, "hedges skipped")
        if applies:
            try:
                with self._hedge_slot(hedge_wait) as slot:
                    if slot:
                        rows = self._build_seed_agnostic(contract, reference, cell_types, budget)
                    else:
                        logger.warning("Build: another build holds the hedge slot; using the "
                                       "ordinary construction")
                        rows = None
            except Exception as exc:
                logger.warning(f"Build: seed-agnostic builder failed ({exc}); using the "
                               "ordinary construction")
                logger.debug(traceback.format_exc())
                rows = None
            if rows:
                self._persist(settings.MINER_SUBMISSION_PATH, rows)
                logger.info(f"Build: wrote the submission to {settings.MINER_SUBMISSION_PATH}")
                return rows
        else:
            logger.info(f"Build: seed-agnostic path not used ({why_not})")

        if not contract.get("seed"):
            if contract.get("cell_type") == "HEK293":
                logger.info(
                    "Contract carries no seed — HEK293 row selection is seed-agnostic; only the "
                    "local diagnostic outcome score is provisional"
                )
            else:
                # Stage 3 is seeded from contract.seed, so an unstamped task's construction stops
                # holding the moment the backend assigns a real one. The rows stay valid;
                # consistency_factor does not survive.
                logger.warning(
                    "Contract carries no seed — the outcome construction is provisional and the "
                    "predicted consistency will not hold if the backend restamps the task"
                )

        logger.info("Build: building context (k-mer index over gene_region +/- 50 kb)")
        context = G.build_context(contract, reference, cell_types)
        logger.info(
            f"Build: context ready | mutations={len(context.mutations)} "
            f"cas={context.cas_systems} max_experiments={context.max_experiments}"
        )

        # Resolve the per-contract settings before enumeration so a future cell-type override of
        # FLANK or LENGTHS affects the site pool as well as the downstream generator.
        cfg = self.gen_config_for(contract)
        lengths = tuple(sorted(set(cfg.lengths)))
        logger.info(f"Build: enumerating PAM sites in gene_region +/- {cfg.flank}")
        sites = G.enumerate_sites(context, cfg.flank, lengths)
        logger.info(f"Build: {len(sites)} PAM sites enumerated")

        override = self.CELL_TYPE_OVERRIDES.get(contract.get("cell_type"))
        if override:
            logger.info(
                f"Build: cell type {contract.get('cell_type')} overrides {sorted(override)} — "
                f"its accessibility cannot reach the cut_p ceilings the tuned config assumes"
            )
        logger.info(
            f"Build: config resolved | construction={cfg.construction} selection={cfg.selection} "
            f"gc_tolerance={cfg.gc_tolerance} cut_p_floors={cfg.cut_p_floors} "
            f"cas_share={cfg.cas_share} backfill_cas={cfg.backfill_cas}"
        )
        is_hek293 = contract.get("cell_type") == "HEK293"
        if is_hek293:
            logger.info(
                "Build: generating HEK293 rows from independent seed-agnostic Cas pools "
                "(adaptive fallback plus typed x170 candidate; exact 70/30 Cas and 50/50 strand "
                "margins; Cas-specific length/GC; strict TTTV; global guide uniqueness/entropy; "
                "anonymous seed-free candidate gate)"
            )
            clustered = generate_seed_agnostic_clustered(context, sites, cfg)
            rows, valid, results = clustered.rows, clustered.valid, clustered.results
            cfg = replace(cfg, weight_skew=clustered.weight_skew)
            logger.info(
                f"Build: HEK293 clusters ready | site_clusters={clustered.site_clusters} "
                f"feature_clusters={clustered.feature_clusters} cas={clustered.cas_counts} "
                f"strands={clustered.strand_counts} mutations={clustered.mutation_totals} "
                f"quota_objective={clustered.deterministic_objective:.3f} "
                f"(tws={clustered.deterministic_tws:.3f}, "
                f"fidelity={clustered.deterministic_fidelity:.6f}, "
                f"candidates={clustered.quota_candidates_evaluated}) "
                f"selected={clustered.selected_candidate}"
            )
            gate = clustered.gate_diagnostics
            if gate is not None:
                logger.info(
                    f"Build: HEK293 anonymous gate | estimator={gate.estimator} "
                    f"replicates={gate.replicates} folds={gate.folds} "
                    f"mean_C(fallback={gate.fallback_mean_consistency}, "
                    f"typed={gate.typed_mean_consistency}) "
                    f"mean_final(fallback={gate.fallback_mean_final}, "
                    f"typed={gate.typed_mean_final}) "
                    f"paired_delta_mean={gate.paired_final_delta_mean} "
                    f"se={gate.paired_final_delta_standard_error} "
                    f"lcb={gate.paired_final_delta_lcb} "
                    f"threshold={gate.minimum_lcb_gain} choice={gate.selected_candidate}"
                )
                if gate.typed_candidate_failure:
                    logger.warning(
                        "Build: HEK293 typed candidate unavailable; retained adaptive fallback: "
                        f"{gate.typed_candidate_failure}"
                    )
            for failure in clustered.quota_candidate_failures:
                logger.warning(f"Build: skipped infeasible HEK293 quota candidate: {failure}")
        else:
            if cfg.strategy == "pure":
                # The optimal skew depends on this contract's mutation-weight ratio, which moves
                # task to task; fitting it costs a few ms of surrogate scoring against selection.
                cfg = replace(cfg, weight_skew=G.choose_weight_skew(context, sites, cfg))
                logger.info(f"Build: fitted weight_skew={cfg.weight_skew}")

            logger.info(
                f"Build: generating rows (construction={cfg.construction}, "
                f"{cfg.variants} variants/site)"
            )
            rows, valid, results = G.generate(context, sites, cfg)
            rows = G.order_rows(rows, valid)
        logger.info(f"Build: generated {len(rows)} rows, {len(valid)} of them stage-1 valid")

        strategy_label = "hek293_clustered" if is_hek293 else cfg.strategy
        construction_label = "seed-agnostic" if is_hek293 else cfg.construction
        logger.info(
            f"Built {len(rows)}/{context.max_experiments} rows from {len(sites)} sites "
            f"| strategy={strategy_label} construction={construction_label} "
            f"skew={cfg.weight_skew} cas={dict(Counter(r['cas'] for r in results))} "
            f"outcomes={dict(Counter(r['outcome'] for r in results))}"
        )
        if len(rows) < context.max_experiments:
            # Term 1 is a sum, so a short submission is a straight loss. The backfill should have
            # prevented this; reaching here means the eligible pool is genuinely exhausted.
            logger.warning(
                f"Submission is {context.max_experiments - len(rows)} row(s) under the cap even "
                f"after the backfill — widen FLANK, or relax CUT_P_CEILING / CAS_MIX"
            )
        # The normal generator prints these too, but only to stdout; a miner's evidence is its log.
        # HEK293 deliberately does not force a simulated construction, so only apply the common
        # ID/design checks there. Stage-1 validity is guaranteed by the aligned `valid` row count.
        invariant_results = [] if is_hek293 else results
        for problem in G.check_invariants(rows, invariant_results, cfg, valid):
            logger.warning(f"Submission invariant violated: {problem}")

        if self.SCORE_LOCALLY and len(valid) >= 2:
            logger.info("Build: scoring locally through the validator's own pipeline")
            report = G.score_rows(valid, results, context)
            logger.info(
                f"Predicted score {report['final_score']:.3f} = "
                f"weighted {report['total_weighted_score']:.3f} "
                f"x consistency {report['consistency_factor']:.4f} "
                f"x fidelity {report['distribution_fidelity_factor']:.4f}"
            )

        # Keep a copy of exactly what was sent. This is also the path benchmark_submission reads, so
        # an operator can re-score the upload against the persisted contract and reference.
        self._persist(settings.MINER_SUBMISSION_PATH, rows)
        logger.info(f"Build: wrote the submission to {settings.MINER_SUBMISSION_PATH}")
        return rows

    def _lengths(self) -> tuple[int, ...]:
        """Guide lengths in the form the site cache is keyed on — sorted and deduped, so a prewarm
        and a build hit the same entry."""
        return tuple(sorted(set(self.gen_config.lengths)))

    @staticmethod
    def _load_sequence() -> None:
        """Load chr11, turning a missing FASTA into an ordinary error.

        ``G.load_sequence`` raises SystemExit when ``data/chr11.fa`` is absent — the right reflex for
        a CLI, wrong for a long-lived miner, where it would escape both the prewarm thread's and
        ``process_task``'s ``except Exception`` and take the round down without a log line.
        """
        try:
            G.load_sequence()
        except SystemExit as missing:
            raise RuntimeError(str(missing)) from missing

    def _fetch_artifacts(self, task: Task) -> tuple[dict, dict]:
        """GET the contract and the HBB reference.

        Plain unsigned GETs — the presigning is already in the URL. Both URLs are short-lived, so
        they are fetched per task and never cached; the *contents* are persisted, which is what lets
        a restart prewarm the site cache and an operator re-score a submission locally.
        """
        contract = self._get_json(task.contract_url, "contract")
        reference = self._get_json(task.hbb_ref_url, "hbb_reference")
        self._persist(settings.CONTRACT_PATH, contract)
        self._persist(settings.HBB_REFERENCE_PATH, reference)
        logger.info(
            f"Persisted contract to {settings.CONTRACT_PATH} and reference to "
            f"{settings.HBB_REFERENCE_PATH}"
        )
        return contract, reference

    def _fetch_cell_types(self) -> dict:
        """The accessibility table, read unsigned from the backend.

        Accessibility scales stage 3's energy, which drives the cut and repair draws, so falling back
        to the default 1.0 risks simulating different outcomes than the validator and breaking the
        construction. How much it costs depends on the clamp: energy saturates at 1.0 for
        near-mutation rows at 50% GC, so a real 0.77 was measured to change nothing, while 0.5 drops
        below the clamp and takes consistency_factor from 1.0 to 0.18. Rows stay valid either way —
        it is term 2 of the score that is at risk.
        """
        try:
            logger.info(f"Fetching cell types from {settings.CELL_TYPES_URL}")
            response = requests.get(settings.CELL_TYPES_URL, timeout=settings.TASK_REQUEST_TIMEOUT)
            response.raise_for_status()
            return response.json()
        except Exception as e:
            logger.warning(
                f"Cell-types fetch failed ({e}); falling back to accessibility 1.0, which will "
                "mis-predict stage 3 and likely cost the whole consistency factor"
            )
            return {}

    def _get_json(self, url: str, what: str) -> dict:
        last_error = None
        for attempt in range(self.MAX_RETRIES):
            try:
                logger.info(f"GET {what} ({self._url_summary(url)})")
                response = requests.get(url, timeout=settings.TASK_REQUEST_TIMEOUT)
                response.raise_for_status()
                logger.info(
                    f"GET {what} -> {response.status_code}, {len(response.content) / 1024:.1f} KB"
                )
                return response.json()
            except Exception as e:
                last_error = e
                logger.warning(
                    f"Fetching {what} failed (attempt {attempt + 1}/{self.MAX_RETRIES}): {e}"
                )
                if attempt + 1 < self.MAX_RETRIES:
                    time.sleep(settings.BASE_DELAY_SECONDS * (attempt + 1))
        raise RuntimeError(f"could not fetch {what}: {last_error}")

    def _upload(self, presigned_url: str, rows: list[dict], deadline: float) -> None:
        """PUT the bare JSON array to the validator's bucket.

        The URL must be sent exactly as received — a re-encoded query string is a
        ``SignatureDoesNotMatch`` — and so must the header set, which is *not* free to vary: see
        ``_upload_headers`` for why sending a Content-Type can break the signature outright.
        """
        payload = json.dumps(rows).encode()
        headers = self._upload_headers(presigned_url)
        last_error = None
        logger.info(
            f"PUT {len(payload) / 1024:.1f} KB to {self._url_summary(presigned_url)} "
            f"with headers {headers or '{}'}"
        )

        for attempt in range(self.MAX_RETRIES):
            remaining = deadline - time.time()
            if remaining <= 0:
                raise RuntimeError(
                    "presigned URL expired before the upload went through "
                    f"(last error: {last_error})"
                )
            try:
                logger.info(
                    f"Upload attempt {attempt + 1}/{self.MAX_RETRIES}, {remaining:.0f}s of TTL left"
                )
                response = requests.put(
                    presigned_url,
                    data=payload,
                    headers=headers,
                    timeout=min(remaining, 60),
                )
                response.raise_for_status()
                logger.info(
                    f"Uploaded {len(rows)} rows ({len(payload) / 1024:.1f} KB) -> "
                    f"{response.status_code}, etag {response.headers.get('ETag', '?')}"
                )
                return
            except Exception as e:
                last_error = e
                # S3 answers with XML, and the body is what separates an expired URL from an altered
                # one — the status is 403 either way.
                body = getattr(getattr(e, "response", None), "text", "") or ""
                logger.warning(
                    f"Upload attempt {attempt + 1}/{self.MAX_RETRIES} failed: {e} {body[:300]}"
                )
                if attempt + 1 < self.MAX_RETRIES:
                    time.sleep(settings.BASE_DELAY_SECONDS * (attempt + 1))

        raise RuntimeError(f"upload failed after {self.MAX_RETRIES} attempts: {last_error}")

    @staticmethod
    def _upload_headers(presigned_url: str) -> dict:
        """The headers this PUT may carry, decided by which signing scheme minted the URL.

        The validator presigns ``put_object`` with no ``ContentType``, and botocore answers that
        with a **SigV2** URL (``AWSAccessKeyId`` / ``Signature`` / ``Expires``). SigV2's
        string-to-sign is ``VERB\\n Content-MD5\\n Content-Type\\n Expires\\n resource`` — the
        Content-Type is *in the signature*, signed as the empty string. Sending
        ``Content-Type: application/json`` therefore makes S3 hash a different string than the
        validator did and reject the upload with ``SignatureDoesNotMatch``, which reads like a
        credentials problem and is in fact this header. So V2 URLs get no headers at all.

        SigV4 URLs (``X-Amz-SignedHeaders``) cover only the headers they name, so there a
        Content-Type is required exactly when it was signed — and forbidden otherwise, for the
        same reason. ``requests`` adds no Content-Type of its own for a bytes body, so an empty
        dict really does send none.
        """
        try:
            query = parse_qs(urlparse(presigned_url).query)
        except Exception:
            return {}

        signed = query.get("X-Amz-SignedHeaders")
        if signed:
            names = [name.strip().lower() for name in signed[0].split(";")]
            return {"Content-Type": "application/json"} if "content-type" in names else {}
        return {}

    @staticmethod
    def _upload_deadline(presigned_url: str, now: float) -> float:
        """When S3 stops accepting the PUT.

        The URL carries its own answer, which beats assuming it was minted the instant it arrived:
        ``Expires`` (an absolute epoch) on the SigV2 URLs the validator currently mints, and
        ``X-Amz-Date`` + ``X-Amz-Expires`` on SigV4 ones. Either comes off the *validator's* clock,
        so the result is clamped to our own SUBMISSION_TIMEOUT budget: a slow local clock must not
        talk us into a deadline that has in fact already passed.
        """
        fallback = now + settings.SUBMISSION_TIMEOUT
        try:
            query = parse_qs(urlparse(presigned_url).query)
            if "Expires" in query:
                expires_at = float(query["Expires"][0])
                scheme = "SigV2"
            else:
                signed_at = datetime.strptime(query["X-Amz-Date"][0], "%Y%m%dT%H%M%SZ").replace(
                    tzinfo=timezone.utc
                )
                expires_at = signed_at.timestamp() + int(query["X-Amz-Expires"][0])
                scheme = "SigV4"
            logger.info(
                f"{scheme} presigned URL expires at "
                f"{datetime.fromtimestamp(expires_at, timezone.utc).isoformat()}; "
                f"{expires_at - now:.0f}s of it left on our clock"
            )
            return min(expires_at, fallback)
        except Exception as e:
            logger.warning(
                f"Could not read the URL's expiry ({e}); assuming the local "
                f"SUBMISSION_TIMEOUT of {settings.SUBMISSION_TIMEOUT}s"
            )
            return fallback

    @staticmethod
    def _url_summary(url: str) -> str:
        """host + path of a presigned URL. The query carries the signature, so it stays out of the
        log — the bucket and key are what identify a submission, and the expiry is logged
        separately by ``_upload_deadline``."""
        try:
            parsed = urlparse(url)
            return f"{parsed.scheme}://{parsed.netloc}{parsed.path}"
        except Exception:
            return "<unparseable url>"

    @staticmethod
    def _build_key(task_id: str, contract: dict) -> str:
        """Memo key. The contract is hashed in as well, so a contract that changes under a task id
        rebuilds instead of re-uploading rows designed against the old rules."""
        canonical = json.dumps(contract, sort_keys=True, separators=(",", ":")).encode()
        return f"{task_id}:{hashlib.sha256(canonical).hexdigest()[:16]}"

    @staticmethod
    def _persist(path: str, document) -> None:
        """Write via a temp file and rename: the prewarm on the next restart reads these back, and a
        task arriving mid-write must not leave half a JSON document behind."""
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        temporary = f"{path}.tmp"
        with open(temporary, "w") as handle:
            json.dump(document, handle)
        os.replace(temporary, path)

    async def blacklist(self, caller_hotkey: str) -> bool:
        """
        Determines whether an incoming request should be blacklisted.

        Args:
            caller_hotkey: ss58 hotkey of the caller (already verified by http_auth).

        Returns:
            bool: True if the request should be rejected.
        """
        if caller_hotkey not in self.metagraph.hotkeys:
            if not self.config.blacklist.allow_non_registered:
                logger.debug(f"Blacklisting un-registered hotkey {caller_hotkey}")
                return True

        uid = self.metagraph.hotkeys.index(caller_hotkey)

        if self.config.blacklist.force_validator_permit:
            if not self.metagraph.neurons[uid].validator_permit:
                logger.warning(f"Blacklisting non-validator hotkey {caller_hotkey}")
                return True

        logger.info(f"Allowing recognized hotkey {caller_hotkey} (uid {uid})")
        return False


# This is the main function, which runs the miner.
if __name__ == "__main__":
    with Miner() as miner:
        while True:
            logger.info(f"Miner running... {time.time()}")
            time.sleep(5)
