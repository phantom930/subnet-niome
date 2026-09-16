#!/usr/bin/env bash
# group 100 / k=10 and group 80 / k=10 on the SAME 10 tasks as conj_fleet.json.
# CF_ONLY pins the task set -- the default "latest N" is a moving window and would drift again.
set -u
cd /root/workspace/subnet-niome
L=/tmp/claude-0/-root-workspace-subnet-niome/28230e9a-8a8e-4365-bd4a-ba4d01caff4a/scratchpad
T="1ab688d3,2bb728d0,33c266d5,56506921,800d04c4,8457e626,920e1837,9f14cd8d,aa071e5e,fe564247"
for spec in "100 g100k10" "80 g80k10"; do
  set -- $spec
  echo "=== $(date -Is) group $1 k=10 ==="
  CF_ONLY="$T" CF_GROUP=$1 CF_K=10 CF_OUT=conj_fleet_$2.json \
    .venv/bin/python conj_fleet.py > "$L/fleet_$2.log" 2>&1
  rc=$?; echo "=== $(date -Is) group $1 k=10 exit $rc ==="
done
echo "=== $(date -Is) BOTH DONE ==="
