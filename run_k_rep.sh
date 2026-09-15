#!/usr/bin/env bash
# HEK293 k sweep: 6 vs 8 vs 9 at group 80 / width 300 / light 12, over 12 contracts, own-field.
# Writes conj_replicate_HEK293_k.json -- NOT the canonical stem. Do not edit while running.
set -u
cd /root/workspace/subnet-niome
L=/tmp/claude-0/-root-workspace-subnet-niome/28230e9a-8a8e-4365-bd4a-ba4d01caff4a/scratchpad
echo "=== $(date -Is) HEK293 k sweep 6/8/9 ==="
CR_CELL=HEK293 CR_N=12 CR_NS=12 CR_SUFFIX="_k" \
  CR_ARMS="6/80/300/12,8/80/300/12,9/80/300/12" \
  .venv/bin/python conj_replicate.py > "$L/krep_HEK293.log" 2>&1
rc=$?; echo "=== $(date -Is) exit $rc ==="
