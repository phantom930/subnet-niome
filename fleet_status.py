#!/usr/bin/env python
"""Tally the 9-hotkey fleet's handling of the latest round: did prefetch win, or did any fall to
the in-TTL fallback (the sign of GPU contention)?

Reads each hotkey's pm2 error log (INFO lines land on stderr) and, for the newest task each saw,
reports whether the prefetch had rows ready before the validator called, how long the build took,
what path served the upload, and any all-HDR failure / OOM / retry. Run it any time after a round.

  python fleet_status.py            # the task the fleet most recently worked on
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
seen_all = []
for logf in sorted(glob.glob(f"{LOGDIR}/miner-h*-error.log"),
                   key=lambda p: int(re.search(r"miner-h(\d+)-", p).group(1))):
    h = re.search(r"(miner-h\d+)-", logf).group(1)
    ev = defaultdict(dict); tasks_seen = []; last_ts = None
    for line in open(logf, errors="ignore"):
        t = ts(line)
        if t: last_ts = t
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
    if not tasks_seen: seen_all.append((h, None, {}, None)); continue
    own = tasks_seen[-1]
    seen_all.append((h, own, ev, last_ts))

# --- scope every row to ONE task -----------------------------------------------------------
# **This is the defect that made a dead hotkey read as healthy.** Each hotkey's row used to be
# scoped to `tasks_seen[-1]` -- its OWN newest task -- while the header printed the FIRST row's
# task as though it covered the table. A hotkey that is stopped or deregistered keeps its last
# successful round in the log forever, so it kept reporting "250rows" under the current round's
# header, and the tally counted it as a submission. Observed 2026-09-21: h4 showed 250 rows for
# task 512919ff, which never appears in its log at all -- its real last task was 3917ba07, the
# round before it was stopped, and its final line is `Wallet ... is not registered`.
#
# That is the same trap CLAUDE.md records for the window plan ("a stopped hotkey needs excluding
# from the allocation for exactly the same reason a deregistered one does; the plan file will not
# tell you"), reappearing in the status tool. Coverage read off either overstates the fleet by
# however many are down.
#
# The fleet-wide target is now the task the MOST hotkeys saw most recently, and a hotkey that
# never saw it is reported as such rather than silently falling back to its own stale round.
if want:
    target = next((t for _, t, _, _ in seen_all if t and t.startswith(want)), None)
else:
    freshest = {}
    for _, t, ev_, _ in seen_all:
        if not t: continue
        when = max((v[0] for v in ev_[t].values() if v and v[0]), default=None)
        if when and (t not in freshest or when > freshest[t]): freshest[t] = when
    target = max(freshest, key=freshest.get) if freshest else None

rows = []
for h, own, ev_, last_ts in seen_all:
    if target and ev_.get(target):
        rows.append((h, target, ev_[target] | ({"_oom": ev_["_oom"]} if "_oom" in ev_ else {})))
    else:
        rows.append((h, None, {"_stale": (own, last_ts)}))

print(f"Fleet status — task {target[:8] if target else '?'}\n")
print(f"  {'hk':<4}{'prefetch':>10}{'build':>7}{'validator':>11}{'served':>13}{'submit':>8}  notes")
tally = defaultdict(int)
for h, tid, ev in rows:
    if tid is None:
        st = ev.get("_stale") or (None, None)
        if st[0]:
            age = f", last log {st[1]:%H:%M}" if st[1] else ""
            print(f"  {h[-2:]:<4}{'—':>10}{'':>7}{'':>11}{'':>13}{'':>8}  "
                  f"DID NOT SEE THIS TASK (newest it saw: {st[0][:8]}{age})")
            tally["absent"] += 1
        else:
            print(f"  {h[-2:]:<4}{'—':>10}   no task in log")
        continue
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
print(f"\n  TALLY over {n} hotkeys ON THIS TASK: {tally['ready']} prefetch-ready | "
      f"{tally['via_prefetch']} served from prefetch | {tally['via_fallback']} fell to fallback | "
      f"{tally['oom']} OOM | {tally['submitted']} submitted")
if tally["absent"]:
    print(f"  {tally['absent']} hotkey(s) never saw this task and are EXCLUDED from the tally — "
          f"stopped, deregistered or dead. Do not read the fleet size off the row count.")
