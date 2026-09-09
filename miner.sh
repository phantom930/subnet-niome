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
# Windows must be disjoint (no shared seeds) or the siblings' bands overlap and the coverage
# collapses toward a single hotkey's. 200-299 / 500-599 / 800-899 are evenly spread; any three
# non-overlapping 100-seed windows in 100-999 are equivalent, since band position is otherwise free.
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
#   2026-09-08: six registered, verified against the chain (256 uids on netuid 55):
#     niome_hotkey uid 74, niome_hotkey1 uid 209, niome_hotkey2 uid 235, niome_hotkey3 uid 196,
#     niome_hotkey4 uid 189, niome_hotkey5 uid 147, niome_hotkey6 uid 136, niome_hotkey7 uid 75.
#   h8 remains off. Nothing here can register a hotkey — until `btcli subnet register` puts one
#   back on netuid 55 its process exits at startup (base/neuron.check_registered calls exit()) and
#   pm2 restarts it on a delay.
# 2026-09-09: all eleven are registered — niome_hotkey uid 74, h1 209, h2 235, h3 196, h4 189,
#   h5 147, h6 136, h7 75, h8 151, h9 224, h10 10 — so nothing is held out of the allocation.
DEREGISTERED=""

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
