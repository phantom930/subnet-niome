#!/usr/bin/env bash
# Extend the remaining three cells to 12 contracts, matching what HUDEP-2 just had.
# Sequential: one GPU user at a time.
set -u
cd /root/workspace/subnet-niome
LOG=/tmp/claude-0/-root-workspace-subnet-niome/28230e9a-8a8e-4365-bd4a-ba4d01caff4a/scratchpad
for cell in HEK293 CD34+_HSPC K562; do
  tag="$(echo "$cell" | tr -d '+' | tr '-' '_')"
  echo "=== $(date -Is) $cell n=12 ==="
  CR_CELL="$cell" CR_N=12 CR_NS=12 .venv/bin/python conj_replicate.py > "$LOG/rep12_$tag.log" 2>&1
  rc=$?; echo "=== $(date -Is) $cell n=12 exit $rc ==="
done
echo "=== $(date -Is) ALL n=12 DONE ==="
