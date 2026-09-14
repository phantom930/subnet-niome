#!/usr/bin/env bash
# HEK293 stage B, plus its high-k arms. Separate file rather than an edit to a running script.
set -u
cd /root/workspace/subnet-niome
LOG=/tmp/claude-0/-root-workspace-subnet-niome/28230e9a-8a8e-4365-bd4a-ba4d01caff4a/scratchpad
while pgrep -f "conj_stageb.py" >/dev/null || pgrep -f "conj_grid.py" >/dev/null \
   || pgrep -f "run_conj_stageb2.sh" >/dev/null; do sleep 30; done
echo "=== $(date -Is) HEK293 stage B ==="
CSB_CELL=HEK293 CSB_TOP=12 CSB_NS=20 .venv/bin/python conj_stageb.py > "$LOG/stageb_HEK293.log" 2>&1
rc=$?; echo "=== $(date -Is) HEK293 exit $rc ==="
# k=8 is HEK293's ceiling and beats all-HDR's band of 7, but the w x fid shortlist selects k=6.
echo "=== $(date -Is) HEK293 high-k fill ==="
CSB_CELL=HEK293 CSB_NS=20 \
  CSB_FORCE="8:80:30:12,8:80:100:12,8:100:30:12,8:100:100:12,8:100:225:12,7:100:100:12" \
  .venv/bin/python conj_stageb.py > "$LOG/stagebk_HEK293.log" 2>&1
rc=$?; echo "=== $(date -Is) HEK293 high-k fill exit $rc ==="
echo "=== $(date -Is) HEK293 DONE ==="
