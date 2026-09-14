#!/usr/bin/env bash
set -u
cd /root/workspace/subnet-niome
LOG=/tmp/claude-0/-root-workspace-subnet-niome/28230e9a-8a8e-4365-bd4a-ba4d01caff4a/scratchpad
for cell in K562 HUDEP-2 CD34+_HSPC; do
  tag="$(echo "$cell" | tr -d '+' | tr '-' '_')"
  echo "=== $(date -Is)  stage B $cell ==="
  CSB_CELL="$cell" CSB_TOP=12 CSB_NS=20 .venv/bin/python conj_stageb.py > "$LOG/stageb_$tag.log" 2>&1
  rc=$?; echo "=== $(date -Is)  stage B $cell exit $rc ==="
done
echo "=== $(date -Is)  STAGE B (3 cells) DONE ==="
