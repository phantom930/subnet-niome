#!/usr/bin/env python
"""Tally the 9-hotkey fleet's handling of the latest round: did prefetch win, or did any fall to
the in-TTL fallback (the sign of GPU contention)?

Reads each hotkey's pm2 error log (INFO lines land on stderr) and, for the newest task each saw,
reports whether the prefetch had rows ready before the validator called, how long the build took,
what path served the upload, and any all-HDR failure / OOM / retry. Run it any time after a round.

  python fleet_status.py            # newest task per hotkey
  python fleet_status.py <task8hex> # a specific task id prefix
"""
import re, sys, glob
from datetime import datetime
from collections import defaultdict

LOGDIR = "/root/.pm2/logs"
TS = re.compile(r"^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d),\d+")
def ts(line):
    m = TS.match(line)
    return datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S") if m else None

PATS = {
    "prepare":   re.compile(r"Prefetch: preparing task (\S+)"),
    "ready":     re.compile(r"Prefetch: task (\S+) ready . (\d+) rows in (\d+)s"),
    "pf_fail":   re.compile(r"Prefetch: build for (\S+) failed"),
    "received":  re.compile(r"Received genomics task (\S+)"),
    "used_prep": re.compile(r"using the prepared submission"),
    "waiting":   re.compile(r"prepared build is still running"),
    "in_ttl":    re.compile(r"nothing prepared for this task; building inside the upload TTL"),
    "unusable":  re.compile(r"prepared round unusable"),
    "emergency": re.compile(r"hedges skipped .a prepared build"),
    "submitted": re.compile(r"Submitted (\d+) rows for task (\S+) in ([\d.]+)s"),
    "hdr_ok":    re.compile(r"Build: all-HDR \(([^)]+)\).*group (\d+).*clean (\d+)"),
    "hdr_fail":  re.compile(r"Build: all-HDR failed \(([^)]+)\)"),
    "hdr_decl":  re.compile(r"Build: all-HDR declined"),
    "oom":       re.compile(r"[Oo]ut ?[Oo]f ?[Mm]emory|OutOfMemory|CUDA_ERROR_OUT_OF_MEMORY|cudaErrorMemoryAllocation"),
}

want = sys.argv[1] if len(sys.argv) > 1 else None
rows = []
for logf in sorted(glob.glob(f"{LOGDIR}/miner-h*-error.log"),
                   key=lambda p: int(re.search(r"miner-h(\d+)-", p).group(1))):
    h = re.search(r"(miner-h\d+)-", logf).group(1)
    ev = defaultdict(dict); tasks_seen = []
    for line in open(logf, errors="ignore"):
        t = ts(line)
        for kind, pat in PATS.items():
            m = pat.search(line)
            if not m: continue
            tid = None
            if kind in ("prepare", "ready", "pf_fail", "received"): tid = m.group(1)
            elif kind == "submitted": tid = m.group(2)
            if kind == "oom": ev["_oom"] = t
            if tid:
                if want and not tid.startswith(want): continue
                if tid not in tasks_seen: tasks_seen.append(tid)
                ev[tid][kind] = (t, m.groups())
            elif kind in ("used_prep","waiting","in_ttl","unusable","emergency","hdr_ok","hdr_fail","hdr_decl"):
                # attach to the most recent task seen
                if tasks_seen: ev[tasks_seen[-1]].setdefault(kind, (t, m.groups()))
    if not tasks_seen: rows.append((h, None, {})); continue
    tid = tasks_seen[-1] if not want else next((x for x in tasks_seen if x.startswith(want)), tasks_seen[-1])
    rows.append((h, tid, ev[tid] | ({"_oom": ev["_oom"]} if "_oom" in ev else {})))

tid_common = next((t for _,t,_ in rows if t), None)
print(f"Fleet status — task {tid_common[:8] if tid_common else '?'}\n")
print(f"  {'hk':<4}{'prefetch':>10}{'build':>7}{'validator':>11}{'served':>13}{'submit':>8}  notes")
tally = defaultdict(int)
for h, tid, ev in rows:
    if tid is None:
        print(f"  {h[-2:]:<4}{'—':>10}   no task in log"); continue
    ready = ev.get("ready"); prep = ev.get("prepare"); pf_fail = ev.get("pf_fail")
    recv = ev.get("received"); sub = ev.get("submitted")
    build = f"{ready[1][2]}s" if ready else ("FAIL" if pf_fail else "—")
    prefetch = "ready" if ready else ("FAILED" if pf_fail else ("preparing" if prep else "—"))
    # served path
    if ev.get("used_prep"): served, lane = "prepared", "prefetch"
    elif ev.get("in_ttl"): served, lane = "in-TTL", "fallback"
    elif ev.get("waiting"): served, lane = "waited", "fallback"
    elif ev.get("unusable"): served, lane = "rebuilt", "fallback"
    else: served, lane = ("prepared" if ready else "—"), ("prefetch" if ready else "?")
    # ready before validator?
    intime = ""
    if ready and recv:
        margin = (recv[0] - ready[0]).total_seconds()
        intime = f"+{margin/60:.0f}m" if margin >= 0 else f"LATE {margin:.0f}s"
    elif recv and not ready:
        intime = "no-prep"
    notes = []
    if ev.get("hdr_fail"): notes.append(f"all-HDR FAIL: {ev['hdr_fail'][1][0][:24]}")
    if ev.get("hdr_decl"): notes.append("all-HDR declined")
    if "_oom" in ev: notes.append("OOM seen")
    if ev.get("emergency"): notes.append("emergency build")
    subtxt = f"{sub[1][0]}rows" if sub else ("—" if recv else "no-call")
    print(f"  {h[-2:]:<4}{prefetch:>10}{build:>7}{('called '+intime) if recv else 'no-call':>11}"
          f"{served:>13}{subtxt:>8}  {'; '.join(notes)}")
    tally["ready" if ready else "not_ready"] += 1
    if lane == "prefetch": tally["via_prefetch"] += 1
    elif lane == "fallback": tally["via_fallback"] += 1
    if "_oom" in ev: tally["oom"] += 1
    if sub: tally["submitted"] += 1
n = sum(1 for _,t,_ in rows if t)
print(f"\n  TALLY over {n} hotkeys: {tally['ready']} prefetch-ready | "
      f"{tally['via_prefetch']} served from prefetch | {tally['via_fallback']} fell to fallback | "
      f"{tally['oom']} OOM | {tally['submitted']} submitted")
