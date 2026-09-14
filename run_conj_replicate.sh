#!/usr/bin/env bash
# Six-contract replication, all four cells, sequential (one GPU user at a time).
# Do not edit while running: bash reads scripts incrementally by byte offset.
set -u
cd /root/workspace/subnet-niome
LOG=/tmp/claude-0/-root-workspace-subnet-niome/28230e9a-8a8e-4365-bd4a-ba4d01caff4a/scratchpad
for cell in HEK293 K562 CD34+_HSPC HUDEP-2; do
  tag="$(echo "$cell" | tr -d '+' | tr '-' '_')"
  echo "=== $(date -Is) replicate $cell ==="
  CR_CELL="$cell" CR_N=6 CR_NS=12 .venv/bin/python conj_replicate.py > "$LOG/rep_$tag.log" 2>&1
  rc=$?; echo "=== $(date -Is) replicate $cell exit $rc ==="
done
echo "=== $(date -Is) REPLICATION DONE ==="
