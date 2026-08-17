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
import hashlib
import json
import logging
import os
import requests
import sys
import threading
import time
import traceback

from collections import Counter
from dataclasses import replace
from datetime import datetime, timezone
from urllib.parse import parse_qs, urlparse

# Add project root to Python path. This has to happen before the niome_subnet imports below —
# running this file as a script puts neurons/ on sys.path, not the repo root — and settings has to
# be the first of them, because it fixes BT_NO_PARSE_CLI_ARGS before bittensor is imported.
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

import niome_subnet.utils.settings as settings  # noqa: E402

# genExp is the generator, and submission.py's per-task flow (build_context -> enumerate_sites ->
# choose_weight_skew -> generate -> order_rows) is what _build mirrors. Importing it rather than
# keeping a second copy in the package is deliberate: genExp is where the design is tuned and swept
# against the whole task history, so anything measured there is what the miner submits, with no port
# in between. Note the module chdir()s to the repo root on import — harmless here (every settings.py
# path is relative to it and neurons must run from there anyway) but it is an import-time side
# effect, so it stays below settings and above nothing that cares about cwd.
import genExp as G  # noqa: E402

from niome_subnet.base.miner import BaseMinerNeuron  # noqa: E402
from niome_subnet.genomics.model import Task  # noqa: E402
from niome_subnet.protocol import GenomicsTaskSynapse  # noqa: E402

logger = logging.getLogger(__name__)


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
    SELECTION = "packed"
    # Which construction every row's simulated outcome is forced to satisfy. "mh" is the HDR/NHEJ
    # mix: microhomology rows repair by HDR, the rest are BLUNT_NHEJ pinned to a 1 bp indel. Stage 4
    # recovers that mapping exactly (mh is a column in its feature matrix), so consistency_factor
    # reaches 1.0 without the degenerate single-outcome dataset "hdr" would submit. See
    # G.CONSTRUCTIONS for the alternatives.
    CONSTRUCTION = "mh"
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

    def __init__(self, config=None):
        super(Miner, self).__init__(config=config)

        self.gen_config = G.GenConfig(
            strategy=self.STRATEGY,
            selection=self.SELECTION,
            construction=self.CONSTRUCTION,
            variants=self.VARIANTS,
            flank=self.FLANK,
            lengths=tuple(self.LENGTHS),
        )

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

        logger.info(
            f"Generation config: strategy={self.STRATEGY} selection={self.SELECTION} "
            f"construction={self.CONSTRUCTION} variants={self.VARIANTS} flank={self.FLANK} "
            f"lengths={self.LENGTHS} score_locally={self.SCORE_LOCALLY}"
        )
        threading.Thread(target=self._prewarm, name="niome-prewarm", daemon=True).start()

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

            logger.info(f"{tag} step 4/5 building the dataset")
            build_started = time.time()
            async with self._build_lock:
                key = self._build_key(task.id, contract)
                if self._built is not None and self._built[0] == key:
                    rows = self._built[1]
                    logger.info(f"{tag} reusing the {len(rows)}-row build for this task")
                else:
                    rows = await asyncio.to_thread(
                        self._build, contract, reference, cell_types
                    )
                    # Only the current task's rows are worth keeping: a new task means a new
                    # contract, and the old rows can never be submitted again.
                    self._built = (key, rows)
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
            logger.info(
                f"Submitted {len(rows)} rows for task {task.id} in {time.time() - started:.1f}s "
                f"({deadline - time.time():.0f}s of the URL's TTL to spare)"
            )
        except Exception as e:
            # Nothing downstream reports a miner failure — a missed upload is indistinguishable from
            # never having been contacted, and there is no retry within a task id. So this log line
            # is the only evidence the round was lost.
            logger.error(f"Failed to submit task {getattr(task, 'id', '?')}: {e}")
            logger.debug(traceback.format_exc())

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

    def _build(self, contract: dict, reference: dict, cell_types: dict) -> list[dict]:
        """Generate this task's submission and log what it should be worth.

        The same sequence ``submission.build_for_task`` runs, so a row set built here is byte-for-byte
        the one ``python submission.py --task-id <id>`` produces for the same contract — which is what
        makes an offline sweep a prediction of what the miner will actually send.
        """
        logger.info("Build: loading chr11 (cached after the first call)")
        self._load_sequence()
        if not contract.get("seed"):
            # Stage 3 is seeded from contract.seed, so an unstamped task's construction stops holding
            # the moment the backend assigns a real one. The rows stay valid; consistency_factor does
            # not survive.
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

        logger.info(f"Build: enumerating PAM sites in gene_region +/- {self.gen_config.flank}")
        sites = G.enumerate_sites(context, self.gen_config.flank, self._lengths())
        logger.info(f"Build: {len(sites)} PAM sites enumerated")

        cfg = self.gen_config
        if cfg.strategy == "pure":
            # The optimal skew depends on this contract's mutation-weight ratio, which moves task to
            # task; fitting it costs a few ms of surrogate scoring against the selection alone.
            cfg = replace(cfg, weight_skew=G.choose_weight_skew(context, sites, cfg))
            logger.info(f"Build: fitted weight_skew={cfg.weight_skew}")

        logger.info(f"Build: generating rows (construction={cfg.construction}, {cfg.variants} variants/site)")
        rows, valid, results = G.generate(context, sites, cfg)
        rows = G.order_rows(rows, valid)
        logger.info(f"Build: generated {len(rows)} rows, {len(valid)} of them stage-1 valid")

        logger.info(
            f"Built {len(rows)}/{context.max_experiments} rows from {len(sites)} sites "
            f"| strategy={cfg.strategy} construction={cfg.construction} "
            f"skew={cfg.weight_skew} outcomes={dict(Counter(r['outcome'] for r in results))}"
        )
        # G.generate prints these too, but only to stdout; a miner's evidence is its log.
        for problem in G.check_invariants(rows, results, cfg, valid):
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

        The URL must be sent exactly as received: only ``host`` is covered by the signature, so
        extra headers are fine, but a re-encoded query string is a ``SignatureDoesNotMatch``.
        """
        payload = json.dumps(rows).encode()
        last_error = None
        logger.info(
            f"PUT {len(payload) / 1024:.1f} KB to {self._url_summary(presigned_url)}"
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
                    headers={"Content-Type": "application/json"},
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
    def _upload_deadline(presigned_url: str, now: float) -> float:
        """When S3 stops accepting the PUT.

        ``X-Amz-Date`` plus ``X-Amz-Expires`` in the URL are S3's own answer, which beats assuming
        the URL was minted the instant it arrived. They come off the *validator's* clock though, so
        the result is clamped to our own SUBMISSION_TIMEOUT budget: a slow local clock must not talk
        us into a deadline that has in fact already passed.
        """
        fallback = now + settings.SUBMISSION_TIMEOUT
        try:
            query = parse_qs(urlparse(presigned_url).query)
            signed_at = datetime.strptime(query["X-Amz-Date"][0], "%Y%m%dT%H%M%SZ").replace(
                tzinfo=timezone.utc
            )
            expires_at = signed_at.timestamp() + int(query["X-Amz-Expires"][0])
            logger.info(
                f"Presigned URL signed at {signed_at.isoformat()} for "
                f"{query['X-Amz-Expires'][0]}s; {expires_at - now:.0f}s of it left on our clock"
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
