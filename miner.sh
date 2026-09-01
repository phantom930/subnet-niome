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
  "niome_hotkey1  8092 52096 300-399"
  "niome_hotkey2  8093 52424 500-599"
  "niome_hotkey3  8094 52069 600-699"
  "niome_hotkey4  8095 52011 700-799"
  "niome_hotkey5  8096 52079 800-899"
  "niome_hotkey6  8097 52799 900-999"
)

run_one() {
  # <name> <port> <ext_port> <window>. exec so the process replaces this shell — pm2 then
  # supervises the python directly.
  local name=$1 port=$2 ext=$3 win=$4
  echo "starting $name on :$port (ext :$ext) with NIOME_HDR_WINDOW=$win"
  NIOME_HDR_WINDOW="$win" exec "$PY" neurons/miner.py \
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
