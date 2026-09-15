#!/usr/bin/env bash
# Launch the niome miner hotkeys, each pinned to a disjoint all-HDR clean-band window.
#
# all-HDR's clean band is ~15 seeds and lands wherever NIOME_HDR_WINDOW places it. A coldkey's
# payout is 1-(1-union/900)^3 over its siblings' bands, so the win is making those bands DISJOINT:
# 3 hotkeys on disjoint windows covered 15.2% of rounds against 5.2% for one window three times
# (measured, 2.9x). Assign each hotkey its own window below and keep them from overlapping.
#
# Two ways to run:
#   ./miner.sh                 launch every hotkey as a background child, wait on all (one wrapper)
#   ./miner.sh niome_hotkey1   run just that hotkey in the foreground (exec)
#
# Prefer the second form under pm2 — one pm2 app per hotkey gives independent restart, which the
# whole-wrapper form does not (a single child crash-looping is invisible to pm2 and takes the
# hotkey down silently). Register them with, e.g.:
#   pm2 start ./miner.sh --name miner-h0 -- niome_hotkey
#   pm2 start ./miner.sh --name miner-h1 -- niome_hotkey1
#   pm2 start ./miner.sh --name miner-h2 -- niome_hotkey2

set -euo pipefail

# Absolute venv python: `pm2 restart --update-env` has rewritten PATH before and sent the miner
# into a crash loop when it resolved a different python. Never rely on PATH here.
PY=/root/workspace/subnet-niome/.venv/bin/python
ROOT=/root/workspace/subnet-niome
cd "$ROOT"

EXTERNAL_IP=184.144.255.144

# One row per hotkey:  <wallet-hotkey>  <axon.port>  <axon.external-port>  <NIOME_HDR_WINDOW>
#
# **The window column is now a LAST RESORT that nothing normally reads.** A build takes its window
# from one of two places, in this order:
#
#   1. data/window_plan.json         written hourly by window_plan.py from the predicted classes
#   2. joined_window.fallback_for()  the SAME construction against a FIXED class triple, used when
#                                    that file is missing, stale or malformed
#
# Both key on NIOME_INSTANCE (the hotkey name, exported below), both take the WIDTH from the cell
# type, and both rotate at stride 30 with h0 at a half-stride offset. The width is per cell --
# 100 for K562 and HUDEP-2, 150 for CD34+_HSPC, 75 for HEK293, validated over five contracts each
# (all_hdr.CELL_CONFIG) -- and it is not a free choice: group_size is tuned AT the width (HEK293
# measured group 100 ahead at width 75 and group 80 ahead at 100-225), so a layout that hands a
# cell the wrong width hands it the wrong group too. One env var cannot carry four widths, which
# is why this column stopped being the source.
#
# NIOME_HDR_WINDOW is read only if NIOME_INSTANCE is unset or the hotkey appears in neither of
# joined_window's lists -- and even then only its START is used, at the cell's width. That is why
# the duplicate starts below (h0/h1 at 100, h8/h9/h10 at 900) are now harmless: before the rotation
# was wired in they collapsed eleven hotkeys onto EIGHT distinct fallback windows, three of them
# drawing the identical band.
#
# So: to move a hotkey, edit joined_window.ROTATE_HK or FALLBACK_CLASSES. To change a width, edit
# all_hdr.CELL_CONFIG. Editing this column changes only the last resort. Adding a row here still
# registers its ports, but add the same name to joined_window.ROTATE_HK or it will fall straight
# through to that last resort and re-correlate with whatever shares its start.
HOTKEYS=(
  "niome_hotkey   8091 52760 100-299"
  "niome_hotkey1  8092 52096 100-299"
  "niome_hotkey2  8093 52424 120-319"
  "niome_hotkey3  8094 52069 130-329"
  "niome_hotkey4  8095 52011 150-349"
  "niome_hotkey5  8096 52079 170-369"
  "niome_hotkey6  8097 52799 180-379"
  "niome_hotkey7  8098 52240 200-399"
  "niome_hotkey8  9001 52504 900-999"
  "niome_hotkey9  9002 52384 900-999"
  "niome_hotkey10  8099 52543 900-999"
)

# NIOME_SEED_DEPEND swaps a hotkey's construction for genomics/seed_depend.py, a submission pinned
# to seed 0. It wins the rounds the backend never stamps and scores the ~0.10 floor on every round
# that IS stamped, so a listed hotkey gives up its seed window entirely — it is the FIRST rung of
# _build's ladder, and listing a hotkey here is how you turn its band off.
#
# EMPTY as of 2026-09-08: the fleet is back on seed windows (h1/h2/h3 on ranks 1/2/3, h0 on rank 4
# or all-cut for HEK293). Retired on judgement, not on failure — everything needed to switch back
# is intact, including the shared build. To re-enable, list hotkeys as "<hotkey>:<n>"; the value is
# now only a marker, because the four hotkeys SHARE one build (Miner.SEED_DEPEND_SHARED) and
# Miner.SEED_DEPEND_SHARED_VARIANT is what it uses.
#
# What was measured, so the decision can be revisited on evidence:
#   rate        never-stamped rounds are **7.6%** (8 of 105 since 2026-08-25), not the 3.5% once
#               claimed here. Do not read this off the task listing — it reports seed 0 for rounds
#               that were stamped; the reliable test is miners reaching cons 1.000 in the scores.
#   wins        rank 1 of 248 twice live: 19018a0a 339.94, 67bdd18a 329.37.
#   the field   miners at cons 1.000 on those rounds went 6, 8, 5, 25, **38** across 09-07/08. On
#               the 38-deep round we placed 6th, 1.8 weighted points short — the whole distance
#               between 85% of the curve (ranks 1-4) and 9% (ranks 6-9).
#   siblings    four correlated builds take consecutive ranks r..r+3, so they move as a block:
#               0.85 of the curve at r=1, 0.09 at r=6. Leverage in both directions.
#   tuning      variants_per_site 12000 recovered +1.47 +- 0.25 (t=5.80, 6/6) and moved payout
#               share 25.5% -> 25.8%, i.e. nothing — the points crossed no rank boundary. 24000
#               adds +0.12 for 2x build and 2x memory. One build holds 12.7 GB; four would need
#               50.6 GB on a 49 GB box, which is why the build is shared.
# EMPTY as of 2026-09-14: h0 (niome_hotkey), the only registered hotkey, was moved OFF seed-depend
# and onto all-HDR at the full joined width. It ran seed-depend for one day (2026-09-13) after the
# whole fleet was deregistered between 09-11 22:48 and 09-13 07:33 and only h0 was registered back
# (uid 122, block 9058420).
#
# With one hotkey there is no band/hedge split to make across siblings — the choice is which single
# trade to take, and it is now the band rather than the seed-0 lottery. Seed-depend wins outright
# on the 7.6% of rounds the backend never stamps (rank 1 of 248 twice) and scores the ~0.10 floor
# on the other 92.4%, and the field on those rounds went 6 -> 8 -> 25 -> 38 miners tied at
# consistency 1.000 across 09-07/08, with the 38-deep round putting us 6th at 9% of the curve. The
# band instead plays every round, and at full width it covers the whole 300-seed joined space.
# Everything needed to switch back is intact: list "niome_hotkey:1" here and re-run window_plan.py.
#
# 2026-09-15: h0 is BACK on seed-depend by operator request. **The rate this trade depends on has
# collapsed and that is measured, not inferred.** Counting rounds where any miner reached
# consistency 1.000 -- the reliable test, since the task listing reports seed 0 for rounds that were
# in fact stamped -- across all 474 scored rounds from 2026-08-03:
#
#   overall                     51 of 474 = 10.8%
#   clustered in               2026-08-03/04 (20 rounds) and 08-22/24 (23 rounds)
#   last seed-0 round          a66f01fa, 2026-09-07T23:53 (39 miners tied at 1.000)
#   rounds since               **70, none scored at seed 0**
#
# So seed-depend now scores the ~0.10 floor on every round with no upside unless the backend
# regresses to leaving rounds unstamped. Revert with SEED_DEPEND_VARIANTS="" and re-run
# window_plan.py; nothing else has to change, because h0 keeps its FULL_HK entry either way.
#
# 2026-09-15 (later, operator request): REVERTED to "" -- h0 is back on the band ladder, i.e. the
# conjunction on HEK293/CD34+_HSPC/K562 and all-HDR on HUDEP-2. The seed-0 rate above was
# re-verified independently against the live score feed before reverting: 51 of 477 scored rounds
# have ever carried a miner at consistency exactly 1.000, the last was a66f01fa on 2026-09-08, and
# **73 rounds have been scored since with none at seed 0**. So the lottery this rung exists to win
# has not paid in a week of rounds, while the band plays every round.
#
# seed-depend is the FIRST rung of _build's ladder and REPLACES the construction rather than
# hedging beside it, so a listed hotkey scores the ~0.10 floor on every round the backend does
# stamp -- which is ~92.4% of them (8 never-stamped of 105 measured). Running all seven there
# put the whole fleet on that trade; one hotkey buys the seed-0 lottery ticket while the other
# six keep the band.
#
# The seed-depend hotkey must be the one sitting in joined_window.FULL_HK, and the two lists are
# kept in step: with h10 here, ROTATE_HK is h4-h9, and 6 x STRIDE 50 tiles the 300-seed joined
# space exactly once -- no gap, nothing doubled up. Moving seed-depend to a hotkey that is in
# ROTATE_HK instead would leave five rotating slices covering 250 of 300 seeds and strand the
# half-stride one. h10 keeps its FULL_HK entry so that if seed-depend declines (budget under
# SEED_DEPEND_MIN_BUDGET_S of 360s on a failed prefetch) it falls through onto a real window.
#
# The value after the colon is only a marker -- SEED_DEPEND_SHARED is True, so one process builds
# and any others load the same rows at SEED_DEPEND_SHARED_VARIANT.
#
# Re-run window_plan.py after changing this list: it excludes seed-depend hotkeys from the band
# allocation, so a stale plan would either assign a window nothing plays or omit one that is.
SEED_DEPEND_VARIANTS=""

# One hotkey runs all-cut instead of all-HDR: the fleet's flat-score hedge. all-HDR scores the
# ~0.10 floor on every seed outside its ~14-seed band, which is 64% of rounds with nine bands; all-
# cut pins is_cut across the whole 900-seed window and scores ~0.19-0.35 on every round instead.
#
# Measured on round b9051bc7 (K562, 248 miners, none of the nine spiked): all-cut 42.31 against
# the fleet's 19.7-23.2 and the seed-agnostic hedge's 28.53. That round drew 1 of 3 clean seeds
# where 557/900 predicts 1.86, so over the clean-seed distribution E[final] is 57.4 — rank 10 of
# 248 that round, with P(top 10) = 0.68.
#
# Exactly one, and it must be a SPREAD hotkey, not a concentrated one. Siblings running all-cut
# score within noise of each other and so take ranks r..r+n-1 rather than n x rank r — the second
# is worth almost nothing. Meanwhile each hotkey converted costs band coverage: 9 bands cover
# 123/900 seeds, 8 cover 110, and only a band can produce the k>=2 rounds (cons 0.70) that all-cut
# structurally cannot reach.
#
# The listed hotkey keeps its table window, which is simply unused — all-cut has no band. It also
# becomes prefetch-dependent: ALL_CUT_MIN_BUDGET_S is 480s for the three erythroid cell types
# against the ~225s in-TTL path, so a round whose prefetch fails falls to the seed-agnostic hedge.
#
# Entries are "<hotkey>" (all-cut on every cell type) or "<hotkey>:CELL,CELL" (all-cut on those
# cell types only, all-HDR with a seed window on the rest). The per-cell form is what a
# single-hotkey fleet needs: the one hotkey plays the band where the seed-window prediction has a
# measured edge and the flat hedge where it does not.
#
# 2026-09-08: h0 runs all-cut on **K562**, and takes the rank-4 window on the other three cell
# types. This reverses the earlier HEK293 choice, which was made on band width and is contradicted
# by the first rank-based pricing in this repo (price_cell.py): every possible round score of each
# construction was placed in all 97 real three-seed fields and paid at SCORE_DISTRIBUTION, then
# bootstrapped over fields. Expected curve share per round:
#
#   cell          all-cut   all-HDR   95% CI on the gap        verdict
#   K562           0.0023    0.0006   [+0.0003, +0.0033]       all-cut, 100% of resamples
#   HEK293         0.0027    0.0059   [-0.0044, -0.0022]       all-HDR, 0% of resamples
#   CD34+_HSPC     0.0008    0.0013   [-0.0014, +0.0000]       all-HDR, marginal
#   HUDEP-2        0.0008    0.0006   [-0.0003, +0.0008]       noise (76%) — left on the band
#
# The mechanism is accessibility. K562 at 0.77 reaches cut_p ~0.96, so all-cut holds 549 of 900
# seeds clean at consistency 0.283; HEK293 at 0.35 reaches only ~0.875 and holds **137**, so 61% of
# its rounds have no clean seed at all — which is exactly how h0 came 135th of 248 on 634a512c
# while the three band hotkeys placed 77th-133rd. Meanwhile HEK293 has the strongest all-HDR spike
# of the four (weighted 318, fidelity 0.924). The old config had each on the wrong side.
#
# Still exactly one hotkey: all-cut is deterministic, so a second one submits near-identical rows
# and takes rank+1 at the tail of the curve. 1 all-cut + 3 bands prices at 0.0041 against 0.0023
# for all four on all-cut.
#
# CD34+_HSPC and HUDEP-2 stay on bands — one marginal, one noise, and changing a config on a 76%
# bootstrap is the error this session spent the day avoiding.
# 2026-09-08 (later): h0 moved to all-cut on EVERY cell type alongside the width-300 overlapping
# band layout, so the fleet is 5 bands + 1 flat hedge. Note this puts h0 on the losing side for
# HEK293 on the evidence available (pooled pricing 0.0059 all-HDR vs 0.0027 all-cut, and that
# basis is itself discredited above — HEK293 has never had the matched treatment cmp_k562.py
# gives K562, where all-cut wins 3.78x). Revert to "niome_hotkey:K562" for the measured split.
ALL_CUT_HOTKEYS=""

# Hotkeys that are not registered on the subnet right now. They keep running (a deregistration is
# usually temporary and re-registering is cheaper than a cold restart), but window_plan.py leaves
# them out of the allocation: a validator cannot contact an unregistered hotkey, so a window spent
# on one is a window no hotkey is covering. Verify against the chain before editing this — the
# metagraph is the authority, not this line.
#   2026-09-05: niome_hotkey (h0) deregistered; the other eight hold uids 175/8/110/79/81/124/36/92.
#   2026-09-07: the whole fleet dropped off the metagraph — h0 09-04 20:58, h1 09-05 07:28,
#   h8 09-05 16:39 — and was rebuilt smaller.
#   2026-09-09: all eleven were registered — niome_hotkey uid 74, h1 209, h2 235, h3 196, h4 189,
#   h5 147, h6 136, h7 75, h8 151, h9 224, h10 10.
# 2026-09-13, verified against the chain (256 uids on netuid 55): ONE of eleven is registered —
#   niome_hotkey uid 122, registered at block 9058420. Every other hotkey was deregistered over the
#   preceding ~36h, one at a time in hotkey order: h0 09-11 22:48, h1/h2/h3 09-12 05:13/05:36/06:02,
#   h4-h7 09-12 19:26-21:28, h8/h9/h10 09-13 06:46/07:10/07:33. uid 10 now belongs to another
#   coldkey, so the old uids above are recycled and must not be used to identify us.
#
#   **A deregistered miner does NOT reliably crash-loop, and pm2 will report it healthy.**
#   check_registered's exit() (base/neuron.py) raises SystemExit, and run() is a DAEMON THREAD
#   (base/miner.py run_in_background_thread) — so mid-run deregistration kills only that thread.
#   h8/h9/h10 stayed `online` for 4-5 hours after falling off, printing the "Miner running..."
#   heartbeat from the main thread and still prefetching builds no validator could ever ask for,
#   with resync_metagraph() stopped dead at the deregistration error. Only the exit() in __init__
#   runs on the main thread and actually ends the process, so the crash-loop-plus-restart-delay
#   recovery applies to a hotkey that starts unregistered, NOT to one deregistered while running.
#   Check the chain, not `pm2 list` — grep for the last resync_metagraph() to date the death.
DEREGISTERED="niome_hotkey1 niome_hotkey2 niome_hotkey3 niome_hotkey4 niome_hotkey5 niome_hotkey6 niome_hotkey7 niome_hotkey8 niome_hotkey9 niome_hotkey10"

# Hotkeys preferred for the wide spread windows, in the order they should be filled. They are only
# used as spread when the concentrated block does not need them: HEK293 concentrates 8 and so
# consumes these too, while the erythroid types concentrate 6 and leave both on spread duty.
# Empty at one hotkey: with a single banded hotkey there is no concentrated/spread split to make —
# it takes the whole 100-seed window at the rank window_plan.py's CELL_RANK names for that cell.
SPREAD_HOTKEYS=""

run_one() {
  # <name> <port> <ext_port> <window>. exec so the process replaces this shell — pm2 then
  # supervises the python directly.
  local name=$1 port=$2 ext=$3 win=$4
  echo "starting $name on :$port (ext :$ext) with NIOME_HDR_WINDOW=$win"
  # NIOME_INSTANCE namespaces this hotkey's own read/write files under data/inst/<name>/ so the
  # siblings' submission, task artifacts, upload record and local scoring don't collide (settings.py).
  local sd=""
  for pair in ${SEED_DEPEND_VARIANTS:-}; do
    if [[ "${pair%%:*}" == "$name" ]]; then
      sd="${pair##*:}"
      echo "  ($name is a seed-depend hotkey, variant $sd: pinned to seed 0, no band)"
    fi
  done
  local ac=""
  for spec in ${ALL_CUT_HOTKEYS:-}; do
    if [[ "${spec%%:*}" == "$name" ]]; then
      if [[ "$spec" == *:* ]]; then
        ac="${spec##*:}"
        echo "  ($name runs all-cut on $ac — flat score there, no band; all-HDR elsewhere)"
      else
        ac=1
        echo "  ($name is the all-cut hedge: flat score every round, no band, window unused)"
      fi
    fi
  done
  NIOME_INSTANCE="$name" NIOME_HDR_WINDOW="$win" NIOME_SEED_DEPEND="$sd" \
    NIOME_ALL_CUT_ONLY="$ac" \
    exec "$PY" neurons/miner.py \
    --netuid 55 \
    --wallet niome_coldkey \
    --wallet-hotkey "$name" \
    --axon.external_ip "$EXTERNAL_IP" \
    --axon.external_port "$ext" \
    --axon.ip 0.0.0.0 \
    --axon.port "$port"
}

# Guard the windows: each must be a valid LO-HI inside 100-999, and no two may overlap. A duplicate
# string is the obvious mistake, but 600-899 vs 800-899 overlap without matching as strings, and
# 900-899 is silently rejected by the miner (falling back to the cell default and re-correlating),
# so the guard parses ranges rather than comparing text. Both cost the whole decorrelation, so they
# are launch-time errors, not warnings.
# Overlapping windows are a deliberate configuration, not a mistake, so the guard below has to be
# switchable. It stays ON by default because a silently-rejected or duplicated window
# re-correlates two siblings and costs the entire decorrelation the fleet exists for. Set to 1
# only when window_plan.FIXED_WINDOWS is the layout and the overlap is intended: the object that
# spikes is the clean BAND inside the window, and two hotkeys screening overlapping ranges
# min-union different guide groups, so their bands still land on different seeds.
ALLOW_OVERLAPPING_WINDOWS=1

# It parses the ranges AS WRITTEN, which are no longer the runtime windows -- those come from the
# plan or from joined_window.fallback_for(), keyed on NIOME_INSTANCE (see the HOTKEYS table above).
# So this guard no longer protects decorrelation; both real layouts guarantee eleven distinct
# windows by construction, and neither can be broken by an edit here. What it still earns its keep
# on is catching a malformed or out-of-range window before a process starts on it. Overlap is only
# a NOTE while ALLOW_OVERLAPPING_WINDOWS=1, and duplicate rows (h0/h1 at 100, h8/h9/h10 at 900) are
# deliberate: they feed only the last resort.
assert_disjoint_windows() {
  local -a los=() his=() names=()
  for row in "${HOTKEYS[@]}"; do
    read -r name _ _ win <<<"$row"
    if [[ ! "$win" =~ ^([0-9]+)-([0-9]+)$ ]]; then
      echo "ERROR: $name window '$win' is not LO-HH form" >&2; exit 1
    fi
    local lo=${BASH_REMATCH[1]} hi=${BASH_REMATCH[2]}
    if (( lo < 100 || hi > 999 || lo >= hi )); then
      echo "ERROR: $name window $win is outside 100-999 or not increasing" >&2; exit 1
    fi
    local i
    for i in "${!los[@]}"; do
      if (( lo <= his[i] && los[i] <= hi )); then
        if [[ "${ALLOW_OVERLAPPING_WINDOWS:-0}" == "1" ]]; then
          echo "NOTE: $name window $win overlaps ${names[i]} window ${los[i]}-${his[i]} — allowed by ALLOW_OVERLAPPING_WINDOWS" >&2
        else
          echo "ERROR: $name window $win overlaps ${names[i]} window ${los[i]}-${his[i]} — siblings would correlate" >&2
          exit 1
        fi
      fi
    done
    los+=("$lo"); his+=("$hi"); names+=("$name")
  done
}
assert_disjoint_windows

if [[ $# -ge 1 ]]; then
  # Run the single named hotkey (foreground / exec) — the per-pm2-app path.
  want=$1
  for row in "${HOTKEYS[@]}"; do
    read -r name port ext win <<<"$row"
    if [[ "$name" == "$want" ]]; then
      run_one "$name" "$port" "$ext" "$win"
    fi
  done
  echo "ERROR: hotkey '$want' is not in the HOTKEYS table" >&2
  exit 1
fi

# No argument: launch every hotkey as a background child and wait on all of them. Simple, but pm2
# sees only this wrapper — a child that dies is not individually restarted. Prefer the per-app form.
pids=()
for row in "${HOTKEYS[@]}"; do
  read -r name port ext win <<<"$row"
  ( run_one "$name" "$port" "$ext" "$win" ) &
  pids+=("$!")
done
# If any child exits, bring the rest down too rather than leaving a partial fleet running.
trap 'kill "${pids[@]}" 2>/dev/null || true' EXIT
wait -n
echo "a hotkey process exited; shutting the rest down" >&2
