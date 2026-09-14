#!/usr/bin/env bash
# Stage A of the conjunction sweep, all four cell types, sequentially.
# Sequential on purpose: two GPU screens on an 8 GB card risk an OOM that also takes the live miners.
# NOTE: never edit this file while it is running -- bash reads a script incrementally by byte
# offset, so an in-place edit can corrupt a later loop iteration. Write a new file instead.
set -u
cd /root/workspace/subnet-niome
LOG=/tmp/claude-0/-root-workspace-subnet-niome/28230e9a-8a8e-4365-bd4a-ba4d01caff4a/scratchpad
for spec in "K562:8,9,10" "HUDEP-2:8,9,10" "CD34+_HSPC:8,9,10" "HEK293:6,7,8"; do
  cell="${spec%%:*}"; ks="${spec##*:}"
  tag="$(echo "$cell" | tr -d '+' | tr '-' '_')"
  echo "=== $(date -Is)  starting $cell  k=$ks ==="
  CG_CELL="$cell" CG_KS="$ks" .venv/bin/python conj_grid.py > "$LOG/grid_$tag.log" 2>&1
  rc=$?; echo "=== $(date -Is)  $cell exit $rc  ==="
done
echo "=== $(date -Is)  ALL CELLS DONE ==="
