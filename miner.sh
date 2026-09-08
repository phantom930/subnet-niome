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
  "niome_hotkey   8091 52760 100-199"
  "niome_hotkey1  8092 52096 200-299"
  "niome_hotkey2  8093 52424 300-399"
  "niome_hotkey3  8094 52069 400-499"
  "niome_hotkey4  8095 52011 500-599"
  "niome_hotkey5  8096 52079 600-699"
  "niome_hotkey6  8097 52799 700-799"
  "niome_hotkey7  8098 52240 800-899"
  "niome_hotkey8  9001 52504 900-999"
)

# EXPERIMENTAL: NIOME_SEED_DEPEND swaps a hotkey's construction for genomics/seed_depend.py, a
# submission pinned to seed 0. It wins the rounds the backend never stamps and scores the ~0.10
# floor on every round that IS stamped, so each listed hotkey gives up its seed window entirely.
#
# The value is NOT read off the task listing: /api/v3/tasks reports `seed: 0` for rounds that were
# in fact stamped (9 list that way, only 4 were scored on it). The reliable test is miners reaching
# cons 1.000 in the score rows. By that test the rate is 2 of 57 rounds since 2026-08-25 = 3.5%.
#
# Siblings compete for the same slots: n hotkeys take ranks 1..n, not n x rank 1. So 1 hotkey is
# worth 0.035 x 0.30 = 0.0105/round against the ~0.0027 a band hotkey contributes, and 3 hotkeys
# 0.035 x 0.70 = 0.0245 against ~0.0081 — each one added is worth less than the last, and the
# fleet's band coverage falls from 9 windows to 6 (117 -> 78 seeds of 900).
#
# The VALUE is the variant index, not a flag. The build is deterministic, so without distinct
# variants siblings would submit byte-identical rows. Measured at variants 0-3: finals spanned
# 338.29-338.39 (0.10 points) while row sets shared only 28-36% of their guides.
# The fleet numbers below were priced for nine hotkeys sharing the curve; with ONE hotkey the trade
# is much simpler, because there are no siblings to crowd ranks 2..n:
#   seed-depend   ~3.5% of rounds x rank 1 (0.300 of the curve) = 0.0105/round
#   one band      band 13 of 900 over 3 seeds = 4.3% of rounds, at rank ~10 = ~0.0027/round
# So at one hotkey seed-depend is worth ~3.9x the band it replaces, and unlike the band it does not
# depend on the seed-window prediction being worth anything. What it does depend on is the field:
# see the competition warning above — the payout is rank 1 of however many miners also pin seed 0.
#
# It is the FIRST rung of _build's ladder, so it replaces the construction rather than hedging
# alongside it. ALL_CUT_HOTKEYS and the window plan below stay configured and simply become the
# fallback for the rounds where seed-depend declines (hedge slot taken, build error, short budget).
#
# FOUR hotkeys on it, and here the "each one is worth less than the last" arithmetic argues FOR
# them rather than against. The earlier note assumed our build lands mid-field, where siblings
# split one rank's worth of curve. Measured on both of 2026-09-07's never-stamped rounds, this
# build beat the entire field (138a47f7 +7.53 over the best of 6 tied miners; e32ec7b3 +2.29 over
# 8), so four near-identical siblings take ranks 1-4 outright:
#   1 hotkey   0.300 of the curve      3 hotkeys  0.700
#   2 hotkeys  0.500                   4 hotkeys  0.850   (a 5th would add only 0.05)
# That holds only while we out-score the field's top. If a competitor passes us the whole block
# slides down together — the siblings are correlated, so this is leverage in both directions.
#
# The variant index no longer differentiates anything: the four hotkeys now SHARE one build
# (Miner.SEED_DEPEND_SHARED), so the first process to reach a round builds it and the others load
# the same rows and submit them verbatim. The values below are kept only to mark which hotkeys run
# seed-depend at all; Miner.SEED_DEPEND_SHARED_VARIANT is what the build actually uses.
#
# Why sharing rather than four distinct builds: at variants_per_site 12000 one build holds
# **12.7 GB** (measured), so four concurrent would need 50.6 GB on a 49 GB box — the sharing is
# what makes that config deployable, not an optimisation. The variant index cost nothing to give
# up: it moved finals by 0.03-0.14 points and never moved a rank, and identical submissions take
# the same consecutive ranks that distinct-but-equal ones did. The field does this openly — on
# a66f01fa its top four rows are byte-equal at 262.92.
SEED_DEPEND_VARIANTS="niome_hotkey:1 niome_hotkey1:1 niome_hotkey2:1 niome_hotkey3:1"

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
# 2026-09-07, one hotkey: all-cut on HEK293 and K562, band elsewhere. perrank.py's exclusive
# per-rank hit rates over 38 scored predictions put CD34+_HSPC's edge at rank 1 (46.2% against the
# 29.8% chance baseline) and HUDEP-2's at rank 3 (66.7%); K562 had no scored prediction at all in
# that window and HEK293's ranks decline roughly monotonically from 38.5%, i.e. no rank worth
# concentrating a lone hotkey on. n is 12-13 rounds per cell, so this is the best available read,
# not a significant result — all-cut's flat ~46 E[final] is the defensible play where the ranking
# has nothing.
ALL_CUT_HOTKEYS="niome_hotkey:HEK293,K562"

# Hotkeys that are not registered on the subnet right now. They keep running (a deregistration is
# usually temporary and re-registering is cheaper than a cold restart), but window_plan.py leaves
# them out of the allocation: a validator cannot contact an unregistered hotkey, so a window spent
# on one is a window no hotkey is covering. Verify against the chain before editing this — the
# metagraph is the authority, not this line.
#   2026-09-05: niome_hotkey (h0) deregistered; the other eight hold uids 175/8/110/79/81/124/36/92.
#   2026-09-07: the whole fleet dropped off the metagraph — h0 09-04 20:58, h1 09-05 07:28,
#   h8 09-05 16:39 — and was rebuilt smaller.
#   2026-09-08: four registered, verified against the chain at block 9019325:
#     niome_hotkey uid 74, niome_hotkey1 uid 209, niome_hotkey2 uid 235, niome_hotkey3 uid 196.
#   h4-h8 remain off. Nothing here can register a hotkey — until `btcli subnet register` puts one
#   back on netuid 55 its process exits at startup (base/neuron.check_registered calls exit()) and
#   pm2 restarts it on a delay.
DEREGISTERED="niome_hotkey4 niome_hotkey5 niome_hotkey6 niome_hotkey7 niome_hotkey8"

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
        echo "ERROR: $name window $win overlaps ${names[i]} window ${los[i]}-${his[i]} — siblings would correlate" >&2
        exit 1
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
