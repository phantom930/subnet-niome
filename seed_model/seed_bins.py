"""Seed distribution of tasks: 9 hundred-wide classes (100-199 ... 900-999).

Counts every seed a task's contract carries into its class, split by cell type
with a TOTAL row, over a task window. Two seed formats live in the API: a bare
int, and a "a,b,c" triple (newer tasks) - a triple contributes three seeds, so
the unit here is a seed, not a task. Seeds outside 100-999 are reported
separately rather than folded into an edge class.

Also measures how tightly a triple is packed: for seeds (a, b, c) it takes the
three pairwise |differences| and keeps the smallest - the closest any two seeds
in that task come to each other - and compares the spread of those against what
three uniform draws would give.

Writes seeds.json and a seeds.html report: class matrix, class totals against
the uniform expectation, a strip plot of the raw seed values, and the
closest-pair distribution.
"""

import argparse
import html
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(ROOT))

from fetch_rank import TASKS_URL, fetch_items, parse_ts

LO, HI, STEP = 100, 999, 100
GAP_STEP = 50                                  # min-gap histogram width
GAP_BINS = [(g, g + GAP_STEP - 1) for g in range(0, 450, GAP_STEP)]


def pair_gaps(seeds):
    """The three |differences| of a triple, smallest first."""
    return sorted(abs(a - b) for i, a in enumerate(seeds) for b in seeds[i + 1:])


def strict_ranks(counts, last_round):
    """Ranks 1..9 over the nine classes, every rank used exactly once.

    Most-selected class first. A tie on the count is broken by recency - the
    class updated in the later round takes the lower rank number - and a tie on
    that by class value, the lower class taking the lower rank number. Classes
    never drawn share last_round 0 and so fall back to class order.

    -> list of nine ranks, indexed by class.
    """
    order = sorted(range(len(counts)),
                   key=lambda b: (-counts[b], -last_round[b], b))
    ranks = [0] * len(counts)
    for position, b in enumerate(order, start=1):
        ranks[b] = position
    return ranks


def build_rounds(tasks, start, end):
    """-> {cell_type: [round, ...]} over the 3-seed tasks, oldest first.

    Round N carries the class counts accumulated over rounds 1..N-1 and, for
    each of round N's own three seeds, the rank its class held in exactly that
    prior state. Keeping both halves on the same prior state is what makes a
    joined row self-checking: the three ranks can be read straight off the nine
    counts beside them.

    The nine classes are ranked 1..9 with no rank repeated - see strict_ranks
    for the tie-break. Round 1 has nothing behind it to rank, so its three rank
    cells are 0.

    Round N also carries how often each rank value had been handed out over
    those same previous rounds: ranks 1..9, plus a tenth bucket for the seeds
    that took no rank at all - round 1's three, and any seed outside 100-999.
    That row sums to 3x(N-1), one entry per seed of every previous round.

    A seed outside 100-999 belongs to no class, so it neither counts nor ranks.
    """
    per_cell = defaultdict(list)
    for t in tasks:
        created = parse_ts(t["created_at"])
        if not (start <= created <= end):
            continue
        contract = t.get("content", {}).get("contract", {}) or {}
        vals = parse_seeds(contract.get("seed"))
        if len(vals) < 2:                      # 3-seed style only
            continue
        per_cell[contract.get("cell_type") or "UNKNOWN"].append(
            (created, t["id"], [v for v, _ in sorted(vals, key=lambda x: x[1])]))

    out = {}
    for cell, entries in per_cell.items():
        entries.sort(key=lambda e: e[0])
        counts, rows = [0] * len(BINS), []
        last_round = [0] * len(BINS)           # round each class last took a seed
        rank_tally = [0] * (len(BINS) + 1)     # ranks 1..9, then "no rank"
        for i, (created, task_id, seeds) in enumerate(entries, start=1):
            bins_of = [bin_of(v) for v in seeds]
            class_ranks = strict_ranks(counts, last_round)
            if i == 1:
                # nothing drawn yet, so there is no distribution to rank against
                ranks = [0, 0, 0]
            else:
                ranks = [class_ranks[b] if b is not None else None
                         for b in bins_of]
            rows.append({
                "round": i,
                "created_at": created.isoformat(),
                "task_id": task_id,
                "seeds": seeds,
                "seed_classes": [f"{BINS[b][0]}-{BINS[b][1]}" if b is not None
                                 else None for b in bins_of],
                "cumulative_counts": list(counts),      # rounds 1..i-1
                "class_ranks": class_ranks if i > 1 else [0] * len(BINS),
                "class_last_updated": list(last_round),
                "seed_ranks": ranks,
                "rank_cumulative_counts": list(rank_tally),
            })
            for b, k in zip(bins_of, ranks):
                if b is not None:
                    counts[b] += 1
                    last_round[b] = i
                # rank 0 (round 1) and None (no class) share the last bucket
                rank_tally[(k - 1) if k else len(BINS)] += 1
        out[cell] = rows
    return out


def rounds_stats(rows):
    """Where the drawn seeds sat in the prior distribution, round 2 onwards.

    Each round's nine ranks are a permutation of 1..9, so a generator ignoring
    the running tally would average exactly 5.0 here.
    """
    ranks = [r for row in rows[1:] for r in row["seed_ranks"] if r]
    if not ranks:
        return {"ranked_seeds": 0, "mean_rank": None, "top3": 0, "bottom3": 0}
    final = rows[-1]["rank_cumulative_counts"] if rows else []
    return {
        "ranked_seeds": len(ranks),
        "mean_rank": round(sum(ranks) / len(ranks), 3),
        "top3": sum(1 for r in ranks if r <= 3),
        "bottom3": sum(1 for r in ranks if r >= 7),
        # the last row's tally covers every round but the last one
        "rank_tally_through_penultimate_round": final,
    }


def p_min_gap_above(g, span=HI - LO + 1):
    """P(all adjacent gaps > g) for 3 uniform draws: (1 - 2g/span)^3."""
    frac = 1 - 2 * g / span
    return frac ** 3 if frac > 0 else 0.0
BINS = [(lo, lo + STEP - 1) for lo in range(LO, HI, STEP)]   # 9 classes
TOTAL = "TOTAL"
# categorical slots 1-4, validated; light aqua/yellow are sub-3:1, so lanes
# carry visible labels and the task table below is the relief view
CELL_HUE_LIGHT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
CELL_HUE_DARK = ["#3987e5", "#d95926", "#199e70", "#c98500"]


def parse_seeds(value):
    """-> [(seed, position)]. Accepts 739 and "150,124,400" alike."""
    if value is None:
        return []
    if isinstance(value, int):
        return [(value, 0)]
    parts = [p.strip() for p in str(value).split(",") if p.strip()]
    out = []
    for i, part in enumerate(parts):
        try:
            out.append((int(part), i))
        except ValueError:
            continue
    return out


def bin_of(seed):
    for i, (lo, hi) in enumerate(BINS):
        if lo <= seed <= hi:
            return i
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", default="2026-08-27T12:00:00")
    ap.add_argument("--end", default="2026-08-31T23:59:59")
    ap.add_argument("--json-out", default=str(ROOT / "seed_model" / "seeds.json"))
    ap.add_argument("--html-out", default=str(ROOT / "seed_model" / "seeds.html"))
    ap.add_argument("--rounds-start", default=None,
                    help="window for the per-cell-type round tables; defaults "
                         "to --start")
    ap.add_argument("--rounds-end", default=None,
                    help="window for the per-cell-type round tables; defaults "
                         "to --end")
    args = ap.parse_args()

    start, end = parse_ts(args.start), parse_ts(args.end)
    if end.hour == end.minute == end.second == 0:
        end = end.replace(hour=23, minute=59, second=59, microsecond=999999)
    r_start = parse_ts(args.rounds_start) if args.rounds_start else start
    r_end = parse_ts(args.rounds_end) if args.rounds_end else end
    if r_end.hour == r_end.minute == r_end.second == 0:
        r_end = r_end.replace(hour=23, minute=59, second=59, microsecond=999999)

    tasks = fetch_items(TASKS_URL)
    print(f"fetched {len(tasks)} tasks")

    records, out_of_range, missing = [], [], 0
    task_seeds = {}
    tasks_in_window, formats = 0, Counter()
    for t in tasks:
        created = parse_ts(t["created_at"])
        if not (start <= created <= end):
            continue
        tasks_in_window += 1
        contract = t.get("content", {}).get("contract", {}) or {}
        raw = contract.get("seed")
        cell = contract.get("cell_type") or "UNKNOWN"
        parsed = parse_seeds(raw)
        if not parsed:
            missing += 1
            continue
        formats["triple" if len(parsed) > 1 else "single"] += 1
        if len(parsed) > 1:
            # the gap view asks how close a task's own seeds come to each other,
            # so it keeps every seed the contract carries - a seed outside
            # 100-999 is excluded from the class matrix but still a real seed
            task_seeds[t["id"]] = {
                "created_at": created.isoformat(),
                "cell_type": cell,
                "seeds": [v for v, _ in sorted(parsed, key=lambda x: x[1])],
            }
        for seed, pos in parsed:
            b = bin_of(seed)
            entry = {
                "task_id": t["id"],
                "created_at": created.isoformat(),
                "cell_type": cell,
                "seed": seed,
                "position": pos,
                "seeds_in_task": len(parsed),
            }
            if b is None:
                out_of_range.append(entry)
                continue
            records.append(dict(entry, bin=f"{BINS[b][0]}-{BINS[b][1]}", bin_index=b))
    records.sort(key=lambda r: (r["created_at"], r["position"]))

    # task-level view: how close the three seeds of one task come to each other
    gaps = []
    for task_id, meta in task_seeds.items():
        seeds_of_task = meta["seeds"]
        d = pair_gaps(seeds_of_task)
        gaps.append({
            "task_id": task_id,
            "created_at": meta["created_at"],
            "cell_type": meta["cell_type"],
            "seeds": seeds_of_task,
            "distances": d,
            "min_gap": d[0],
            "max_gap": d[-1],
        })
    gaps.sort(key=lambda g: g["min_gap"])

    cell_types = sorted({r["cell_type"] for r in records}
                        | {g["cell_type"] for g in gaps})
    rounds = build_rounds(tasks, r_start, r_end)
    counts = defaultdict(lambda: [0] * len(BINS))
    for r in records:
        counts[r["cell_type"]][r["bin_index"]] += 1
        counts[TOTAL][r["bin_index"]] += 1

    seeds = [r["seed"] for r in records]
    repeats = {s: n for s, n in Counter(seeds).items() if n > 1}
    summary = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": TASKS_URL,
        "window": {"start": start.isoformat(), "end": end.isoformat()},
        "classes": [f"{lo}-{hi}" for lo, hi in BINS],
        "totals": {
            "tasks_in_window": tasks_in_window,
            "task_seed_formats": dict(formats),
            "seeds_binned": len(records),
            "seeds_out_of_range": len(out_of_range),
            "tasks_without_seed": missing,
            "distinct_seeds": len(set(seeds)),
            "repeated_seeds": repeats,
            "seed_min": min(seeds) if seeds else None,
            "seed_max": max(seeds) if seeds else None,
        },
        "min_gap": {
            "definition": "per task: min(|a-b|, |a-c|, |b-c|) over its 3 seeds",
            "tasks_measured": len(gaps),
            "min": min((g["min_gap"] for g in gaps), default=None),
            "max": max((g["min_gap"] for g in gaps), default=None),
            "mean": round(sum(g["min_gap"] for g in gaps) / len(gaps), 2) if gaps else None,
            "median": (sorted(g["min_gap"] for g in gaps)[len(gaps) // 2]
                       if gaps else None),
            "uniform_expectation": {"mean": round((HI - LO + 1) / 8, 1),
                                    "median": round((HI - LO + 1) * (1 - 0.5 ** (1 / 3)) / 2, 1)},
            "histogram": {f"{lo}-{hi}": sum(1 for g in gaps if lo <= g["min_gap"] <= hi)
                          for lo, hi in GAP_BINS},
        },
        "task_gaps": gaps,
        "matrix": {ct: {f"{BINS[i][0]}-{BINS[i][1]}": counts[ct][i]
                        for i in range(len(BINS))}
                   for ct in cell_types + [TOTAL]},
        "out_of_range": out_of_range,
        "records": records,
        "rounds": {
            "definition": ("per cell type, one row per 3-seed task in order: the "
                           "nine class counts accumulated over all previous "
                           "rounds, and the rank each of this round's seeds' "
                           "classes held in that same prior state (rank 1 = "
                           "most-selected class so far, ties shared)"),
            "window": {"start": r_start.isoformat(), "end": r_end.isoformat()},
            "stats": {ct: rounds_stats(rows) for ct, rows in rounds.items()},
            "by_cell_type": rounds,
        },
    }
    Path(args.json_out).write_text(json.dumps(summary, indent=2) + "\n")

    write_html(summary, counts, cell_types, records, gaps, rounds, args.html_out)

    print(f"window {start.isoformat()} .. {end.isoformat()}: {len(records)} tasks binned, "
          f"{len(set(seeds))} distinct seeds, range {min(seeds, default='-')}..{max(seeds, default='-')}, "
          f"formats {dict(formats)}")
    hdr = "cell type".ljust(12) + "".join(f"{lo}-{hi}".rjust(10) for lo, hi in BINS) + "   total"
    print(hdr)
    for ct in cell_types + [TOTAL]:
        row = counts[ct]
        print(ct.ljust(12) + "".join(str(v or "-").rjust(10) for v in row)
              + str(sum(row)).rjust(8))
    if gaps:
        mg = summary["min_gap"]
        print(f"closest pair per task: min {mg['min']}, median {mg['median']}, "
              f"mean {mg['mean']} (uniform expects mean {mg['uniform_expectation']['mean']}, "
              f"median {mg['uniform_expectation']['median']})")
        print("  tightest:", ", ".join(
            f"{g['min_gap']} ({'/'.join(str(v) for v in g['seeds'])})" for g in gaps[:5]))
    if repeats:
        print("repeated seeds:", repeats)
    for ct in sorted(rounds):
        st = s_rounds = summary["rounds"]["stats"][ct]
        print(f"rounds {ct:<11} {len(rounds[ct]):>3} rounds, "
              f"mean rank of drawn class {st['mean_rank']}, "
              f"top3 {st['top3']}/{st['ranked_seeds']}, "
              f"bottom3 {st['bottom3']}/{st['ranked_seeds']}")
    if out_of_range:
        print("outside 100-999:", out_of_range)
    print(f"wrote {args.json_out} and {args.html_out}")


def write_html(s, counts, cell_types, records, gaps, rounds, path):
    tot = s["totals"]
    n = tot["seeds_binned"]
    labels = s["classes"]
    rows_all = cell_types + [TOTAL]
    cmax = max((counts[ct][i] for ct in cell_types for i in range(len(BINS))), default=1)
    tmax = max(counts[TOTAL]) if n else 1
    expected = n / len(BINS) if n else 0

    def cell(v, is_total):
        if not v:
            return '<td class="c zero">-</td>'
        mix = 18 + 58 * (v / (tmax if is_total else cmax))
        cls = "c tot" if is_total else "c"
        return (f'<td class="{cls}" style="background:color-mix(in srgb, var(--seq) '
                f'{mix:.0f}%, var(--surface))">{v}</td>')

    matrix_rows = "".join(
        f'<tr><th class="rowh">{"" if ct != TOTAL else ""}'
        + (f'<span class="dot" style="background:var(--c{cell_types.index(ct)})"></span>'
           if ct in cell_types else "")
        + f'{html.escape(ct)}</th>'
        + "".join(cell(counts[ct][i], ct == TOTAL) for i in range(len(BINS)))
        + f'<td class="num rowtot">{sum(counts[ct])}</td></tr>'
        for ct in rows_all)

    # class totals, with the uniform expectation as a reference line
    bars = "".join(
        f'<div class="bar-row"><span class="bar-label">{html.escape(labels[i])}</span>'
        f'<span class="bar-track"><span class="bar-fill" style="width:'
        f'{(counts[TOTAL][i] / tmax * 100) if tmax else 0:.1f}%"></span>'
        f'<span class="expect" style="left:{(expected / tmax * 100) if tmax else 0:.1f}%"></span>'
        f'</span><span class="bar-value">{counts[TOTAL][i] or "-"}</span>'
        f'<span class="bar-sub">{counts[TOTAL][i] - expected:+.1f} vs expected</span></div>'
        for i in range(len(BINS)))

    # strip plot: one lane per cell type, dot per task at its seed value
    span = HI + 1 - LO
    lanes = []
    for idx, ct in enumerate(cell_types):
        pts = [r for r in records if r["cell_type"] == ct]
        dots = "".join(
            f'<span class="dot-pt" style="left:{(r["seed"] - LO) / span * 100:.2f}%;'
            f'background:var(--c{idx})" title="seed {r["seed"]} · {html.escape(r["cell_type"])} '
            f'· {html.escape(r["created_at"][:16])} · task {html.escape(r["task_id"][:8])}"></span>'
            for r in pts)
        lanes.append(
            f'<div class="lane"><span class="lane-label">'
            f'<span class="dot" style="background:var(--c{idx})"></span>{html.escape(ct)}'
            f'<i>{len(pts)}</i></span><span class="lane-track">{dots}</span></div>')
    ticks = "".join(
        f'<span class="tick" style="left:{(v - LO) / span * 100:.2f}%">{v}</span>'
        for v in range(LO, HI + 2, STEP))

    table_rows = "".join(
        f'<tr><td class="mono">{html.escape(r["created_at"][:16])}</td>'
        f'<td><span class="dot" style="background:var(--c{cell_types.index(r["cell_type"])})">'
        f'</span>{html.escape(r["cell_type"])}</td>'
        f'<td class="num"><b>{r["seed"]}</b></td>'
        f'<td class="mono">{html.escape(r["bin"])}</td>'
        f'<td class="num muted">{r["position"] + 1}/{r["seeds_in_task"]}</td>'
        f'<td class="mono muted">{html.escape(r["task_id"])}</td></tr>'
        for r in records)

    rep = tot["repeated_seeds"]
    notes = []
    if rep:
        notes.append("repeated seeds: " + ", ".join(f"{k}×{v}" for k, v in sorted(rep.items())))
    if tot["seeds_out_of_range"]:
        vals = ", ".join(str(o["seed"]) for o in s["out_of_range"][:6])
        notes.append(f"{tot['seeds_out_of_range']} seed(s) outside {LO}–{HI} "
                     f"(value{'s' if tot['seeds_out_of_range'] > 1 else ''} {vals}) — "
                     f"excluded from the classes")
    if tot["tasks_without_seed"]:
        notes.append(f"{tot['tasks_without_seed']} task(s) carried no seed")
    fmt = tot["task_seed_formats"]
    notes.insert(0, f"{tot['tasks_in_window']} tasks in window — "
                    f"{fmt.get('triple', 0)} carry a 3-seed triple, "
                    f"{fmt.get('single', 0)} a single seed")
    note_line = " · ".join(notes)


    # ---- closest pair per task -------------------------------------------
    gap_html = ""
    if gaps:
        mg = s["min_gap"]
        hist = [(lo, hi, sum(1 for g in gaps if lo <= g["min_gap"] <= hi))
                for lo, hi in GAP_BINS]
        gmax = max(c for _, _, c in hist) or 1
        n_g = len(gaps)
        rows_g = []
        for lo, hi, c in hist:
            exp = n_g * (p_min_gap_above(lo) - p_min_gap_above(hi + 1))
            rows_g.append(
                f'<div class="bar-row"><span class="bar-label">{lo}–{hi}</span>'
                f'<span class="bar-track">'
                f'<span class="bar-fill" style="width:{c / gmax * 100:.1f}%"></span>'
                f'<span class="expect" style="left:{min(exp / gmax * 100, 100):.1f}%"></span>'
                f'</span><span class="bar-value">{c or "-"}</span>'
                f'<span class="bar-sub">{c - exp:+.1f} vs uniform ({exp:.1f})</span></div>')

        per_cell = []
        for ct in cell_types:
            vals = sorted(g["min_gap"] for g in gaps if g["cell_type"] == ct)
            if not vals:
                continue
            per_cell.append(
                f'<tr><th class="rowh"><span class="dot" '
                f'style="background:var(--c{cell_types.index(ct)})"></span>'
                f'{html.escape(ct)}</th>'
                f'<td class="num">{len(vals)}</td><td class="num">{vals[0]}</td>'
                f'<td class="num">{vals[len(vals) // 2]}</td>'
                f'<td class="num">{sum(vals) / len(vals):.1f}</td>'
                f'<td class="num">{vals[-1]}</td></tr>')

        gap_rows = [{
            "created": g["created_at"][:16],
            "cell_type": g["cell_type"],
            "cell_i": cell_types.index(g["cell_type"]),
            "seeds": g["seeds"],
            "min": g["distances"][0],
            "mid": g["distances"][len(g["distances"]) // 2],
            "max": g["distances"][-1],
            "task": g["task_id"],
        } for g in gaps]

        gap_html = f"""
<section class="card">
  <h2>Closest pair inside each task</h2>
  <p class="note lead">For every task the three seeds give three distances —
    |a−b|, |a−c|, |b−c| — and the smallest is kept. Over {n_g} tasks that minimum runs
    <b>{mg['min']}</b> to <b>{mg['max']}</b>, median <b>{mg['median']}</b>,
    mean <b>{mg['mean']}</b>. Three uniform draws from {LO}–{HI} would give
    mean {mg['uniform_expectation']['mean']} and median
    {mg['uniform_expectation']['median']}.</p>
  <div>{''.join(rows_g)}</div>
  <p class="note">Dashed line marks the uniform expectation for that band,
    from P(min gap &gt; g) = (1 − 2g/{HI - LO + 1})³.</p>
  <table style="margin-top:18px">
    <thead><tr><th>cell type</th><th class="num">tasks</th><th class="num">min</th>
      <th class="num">median</th><th class="num">mean</th><th class="num">max</th></tr></thead>
    <tbody>{''.join(per_cell)}</tbody>
  </table>
  <h2 style="margin:22px 0 8px">Every task ({n_g})</h2>
  <p class="note" style="margin:0 0 12px">Ordered by closest pair; click a column to
    reorder, click again to reverse.</p>
  <div class="scroll"><table id="gaps">
    <thead><tr><th data-k="created">created</th><th data-k="cell_type">cell type</th>
      <th data-k="seed0">seeds</th>
      <th data-k="min" class="num">min</th><th data-k="mid" class="num">mid</th>
      <th data-k="max" class="num">max</th><th data-k="task">task</th></tr></thead>
    <tbody></tbody>
  </table></div>
</section>

<script>
const GAPS = {json.dumps(gap_rows)};
GAPS.forEach(g => g.seed0 = g.seeds[0]);
let gk = "min", gd = "asc";
const gbody = document.querySelector("#gaps tbody");
const gesc = t => String(t).replace(/[&<>"]/g, c => ({{"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}}[c]));
function gapRender() {{
  const rows = GAPS.map((r, i) => [r, i]).sort(([a, ia], [b, ib]) => {{
    const x = a[gk], y = b[gk];
    let c = typeof x === "string" ? x.localeCompare(y) : x - y;
    if (c === 0) return ia - ib;
    return gd === "asc" ? c : -c;
  }}).map(([r]) => r);
  gbody.innerHTML = rows.map(r => `<tr>
    <td class="mono">${{gesc(r.created)}}</td>
    <td><span class="dot" style="background:var(--c${{r.cell_i}})"></span>${{gesc(r.cell_type)}}</td>
    <td class="mono seed">${{r.seeds.map(v => `<span class="s">${{v}}</span>`).join("")}}</td>
    <td class="num hit">${{r.min}}</td>
    <td class="num">${{r.mid}}</td>
    <td class="num">${{r.max}}</td>
    <td class="mono muted">${{gesc(r.task.slice(0, 8))}}</td></tr>`).join("");
  document.querySelectorAll("#gaps thead th").forEach(th =>
    th.dataset.dir = th.dataset.k === gk ? gd : "");
}}
document.querySelectorAll("#gaps thead th").forEach(th => th.onclick = () => {{
  const k = th.dataset.k;
  const numeric = ["min", "mid", "max", "seed0"].includes(k);
  gd = k === gk ? (gd === "asc" ? "desc" : "asc") : (numeric ? "asc" : "asc");
  gk = k;
  gapRender();
}});
gapRender();
</script>"""

    lv = "".join(f"  --c{i}: {CELL_HUE_LIGHT[i % 4]};\n" for i in range(len(cell_types)))
    dv = "".join(f"  --c{i}: {CELL_HUE_DARK[i % 4]};\n" for i in range(len(cell_types)))

    # ---- per cell type: cumulative class counts joined to this round's ranks
    rnd = s.get("rounds", {}) or {}
    round_cards = []
    for ct in sorted(rounds):
        rows = rounds[ct]
        if not rows:
            continue
        st = rnd.get("stats", {}).get(ct, {})
        peak = max((max(r["cumulative_counts"]) for r in rows), default=1) or 1
        rpeak = max((max(r["rank_cumulative_counts"]) for r in rows), default=1) or 1
        i_ct = cell_types.index(ct) if ct in cell_types else 0
        body = []
        for r in rows:
            cnt = "".join(
                (f'<td class="c" style="background:color-mix(in srgb, var(--seq) '
                 f'{14 + 58 * c / peak:.0f}%, var(--surface))">{c}</td>')
                if c else '<td class="c zero">-</td>'
                for c in r["cumulative_counts"])
            rk = "".join(
                (f'<td class="c rk"><b>{k}</b><i>{v}</i></td>' if k
                 else f'<td class="c rk off"><b>{"0" if k == 0 else "—"}</b>'
                      f'<i>{v}</i></td>')
                for v, k in zip(r["seeds"], r["seed_ranks"]))
            rcnt = "".join(
                (f'<td class="c" style="background:color-mix(in srgb, var(--seq) '
                 f'{14 + 58 * c / rpeak:.0f}%, var(--surface))">{c}</td>')
                if c else '<td class="c zero">-</td>'
                for c in r["rank_cumulative_counts"])
            body.append(f'<tr><td class="num rnd">{r["round"]}</td>'
                        f'<td class="mono muted">{html.escape(r["created_at"][:16])}</td>'
                        f'{cnt}<td class="gapcol"></td>{rk}'
                        f'<td class="gapcol"></td>{rcnt}</tr>')
        note = ""     # round-1 handling is spelled out in the table footnote
        stat_line = ""
        if st.get("ranked_seeds"):
            stat_line = (
                f'Over rounds 2–{len(rows)}, {st["ranked_seeds"]} drawn seeds sat at '
                f'mean rank <b>{st["mean_rank"]}</b> against an even-draw value of '
                f'exactly 5.0: <b>{st["top3"]}</b> fell in a top-3 class '
                f'({st["ranked_seeds"] / 3:.0f} expected), <b>{st["bottom3"]}</b> in a '
                f'bottom-3 one.')
        round_cards.append(f"""
<section class="card">
  <h2><span class="dot" style="background:var(--c{i_ct})"></span>{html.escape(ct)}
    — {len(rows)} rounds</h2>
  <p class="note lead">{stat_line}</p>
  <div class="scroll"><table class="rounds">
    <thead>
      <tr><th class="num" rowspan="2">round</th><th rowspan="2">created</th>
        <th class="grp" colspan="{len(BINS)}">cumulative class counts — rounds 1…N−1</th>
        <th class="gapcol" rowspan="2"></th>
        <th class="grp" colspan="3">rank of this round's own seeds</th>
        <th class="gapcol" rowspan="2"></th>
        <th class="grp" colspan="{len(BINS) + 1}">cumulative rank counts —
          rounds 1…N−1</th></tr>
      <tr>{''.join(f'<th class="c">{lo}–{hi}</th>' for lo, hi in BINS)}
        <th class="c">seed 1</th><th class="c">seed 2</th><th class="c">seed 3</th>
        {''.join(f'<th class="c">{k}</th>' for k in range(1, len(BINS) + 1))}
        <th class="c" title="round 1 and seeds outside the classes">0/—</th></tr>
    </thead>
    <tbody>{''.join(body)}</tbody>
  </table></div>
  <p class="note">Each rank cell shows the class rank with its seed value beneath.
    Round 1 has no prior distribution to rank against, so its three cells are
    <b>0</b>; a seed outside {LO}–{HI} belongs to no class and shows “—”. Both
    land in the last of the ten rank columns, which counts seeds that took no
    rank. The right block tallies every rank the middle block handed out over all
    previous rounds, so its row sums to 3×(N−1) — {3 * (len(rows) - 1)} by the
    last round.</p>
</section>""")

    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Task seed classes</title>
<style>
:root {{
  color-scheme: light;
  --plane: #f9f9f7; --surface: #fcfcfb; --ink: #0b0b0b; --ink-2: #52514e;
  --muted: #898781; --grid: #e1e0d9; --line: #c3c2b7; --track: #eeede8;
  --border: rgba(11,11,11,.10); --seq: #2a78d6;
{lv}}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    color-scheme: dark;
    --plane: #0d0d0d; --surface: #1a1a19; --ink: #fff; --ink-2: #c3c2b7;
    --muted: #898781; --grid: #2c2c2a; --line: #383835; --track: #262625;
    --border: rgba(255,255,255,.10); --seq: #3987e5;
{dv}  }}
}}
:root[data-theme="dark"] {{
  color-scheme: dark;
  --plane: #0d0d0d; --surface: #1a1a19; --ink: #fff; --ink-2: #c3c2b7;
  --muted: #898781; --grid: #2c2c2a; --line: #383835; --track: #262625;
  --border: rgba(255,255,255,.10); --seq: #3987e5;
{dv}}}
* {{ box-sizing: border-box; }}
body {{ margin: 0; padding: 28px clamp(16px, 4vw, 44px) 64px; background: var(--plane);
  color: var(--ink); font: 14px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }}
h1 {{ font-size: 21px; margin: 0 0 6px; letter-spacing: -.01em; }}
h2 {{ font-size: 14px; margin: 0 0 16px; }}
.sub {{ color: var(--ink-2); margin: 0 0 22px; }}
.sub code {{ background: var(--track); padding: 1px 6px; border-radius: 4px; font-size: 12px; }}
.card {{ background: var(--surface); border: 1px solid var(--border); border-radius: 12px;
  padding: 16px 18px; margin-bottom: 16px; }}
.tiles {{ display: grid; gap: 12px 20px; grid-template-columns: repeat(auto-fit, minmax(104px, 1fr)); }}
.tile b {{ display: block; font-size: 22px; line-height: 1.15; }}
.tile span {{ color: var(--muted); font-size: 11.5px; }}
.dot {{ width: 9px; height: 9px; border-radius: 50%; display: inline-block; margin-right: 6px;
  vertical-align: baseline; }}
table {{ width: 100%; border-collapse: collapse; font-variant-numeric: tabular-nums; }}
th, td {{ padding: 6px 8px; border-bottom: 1px solid var(--grid); font-size: 12.5px; text-align: left; }}
thead th {{ color: var(--ink-2); font-size: 11.5px; letter-spacing: .03em; text-transform: uppercase;
  white-space: nowrap; border-bottom: 1px solid var(--line); }}
td.c, th.c {{ text-align: center; width: 62px; }}
td.zero {{ color: var(--muted); }}
td.tot {{ font-weight: 600; }}
th.rowh {{ white-space: nowrap; font-weight: 600; }}
#gaps thead th {{ cursor: pointer; user-select: none; }}
#gaps thead th:hover {{ color: var(--ink); }}
#gaps thead th[data-dir]:after {{ content: " ↕"; color: var(--muted); }}
#gaps thead th[data-dir="asc"]:after {{ content: " ↑"; color: var(--ink); }}
#gaps thead th[data-dir="desc"]:after {{ content: " ↓"; color: var(--ink); }}
td.num, th.num {{ text-align: right; }}
td.rowtot {{ font-weight: 600; border-left: 1px solid var(--grid); }}
tbody tr:last-child th, tbody tr:last-child td {{ border-top: 1px solid var(--line); }}
.mono {{ font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 11.5px; }}
.muted {{ color: var(--muted); }}
.bar-row {{ display: grid; grid-template-columns: 86px minmax(120px, 1fr) 28px 118px;
  align-items: center; gap: 10px; padding: 3px 0; }}
.bar-label {{ font-size: 12px; color: var(--ink-2); font-variant-numeric: tabular-nums; }}
.bar-track {{ position: relative; background: var(--track); border-radius: 4px; height: 12px; }}
.bar-fill {{ display: block; height: 100%; border-radius: 0 4px 4px 0; background: var(--seq);
  min-width: 2px; }}
.expect {{ position: absolute; top: -3px; bottom: -3px; width: 0;
  border-left: 2px dashed var(--line); }}
.bar-value {{ text-align: right; font-weight: 600; font-variant-numeric: tabular-nums; font-size: 12px; }}
.bar-sub {{ color: var(--muted); font-size: 11px; font-variant-numeric: tabular-nums; }}
.lane {{ display: grid; grid-template-columns: 130px 1fr; align-items: center; gap: 10px;
  padding: 5px 0; }}
.lane-label {{ font-size: 12px; color: var(--ink-2); white-space: nowrap; }}
.lane-label i {{ font-style: normal; color: var(--muted); margin-left: 6px; }}
.lane-track {{ position: relative; height: 20px; border-bottom: 1px solid var(--grid); }}
.dot-pt {{ position: absolute; top: 5px; width: 10px; height: 10px; margin-left: -5px;
  border-radius: 50%; box-shadow: 0 0 0 2px var(--surface); }}
.axis {{ position: relative; height: 18px; margin: 4px 0 0 140px; }}
.tick {{ position: absolute; transform: translateX(-50%); color: var(--muted); font-size: 10.5px;
  font-variant-numeric: tabular-nums; }}
.grid-lines {{ position: absolute; inset: 0; pointer-events: none; }}
.note {{ color: var(--muted); font-size: 11.5px; margin: 10px 0 0; }}
.note.lead {{ color: var(--ink-2); font-size: 13px; margin: -4px 0 16px; max-width: 82ch; }}
td.seed {{ white-space: nowrap; }}
td.seed .s {{ display: inline-block; min-width: 30px; text-align: right; color: var(--ink-2);
  background: var(--track); border-radius: 4px; padding: 0 4px; margin-right: 2px; }}
td.num.hit {{ font-weight: 600; color: var(--ink);
  background: color-mix(in srgb, var(--seq) 22%, var(--surface)); }}
th.rowh {{ white-space: nowrap; font-weight: 600; }}
#gaps thead th {{ cursor: pointer; user-select: none; }}
#gaps thead th:hover {{ color: var(--ink); }}
#gaps thead th[data-dir]:after {{ content: " ↕"; color: var(--muted); }}
#gaps thead th[data-dir="asc"]:after {{ content: " ↑"; color: var(--ink); }}
#gaps thead th[data-dir="desc"]:after {{ content: " ↓"; color: var(--ink); }}
.legend {{ display: flex; gap: 16px; flex-wrap: wrap; color: var(--ink-2); font-size: 12px;
  margin: 0 0 16px; }}
.scroll {{ max-height: 60vh; overflow: auto; border: 1px solid var(--border); border-radius: 10px; }}
.scroll thead th {{ position: sticky; top: 0; background: var(--surface); z-index: 1; }}
table.rounds th.grp {{ text-align: center; font-size: 11px; color: var(--muted);
  border-bottom: 1px solid var(--grid); }}
table.rounds td.rnd {{ font-weight: 600; }}
table.rounds td.c {{ width: 40px; }}
table.rounds td.rk {{ line-height: 1.15; }}
table.rounds td.rk b {{ display: block; font-size: 13px; }}
table.rounds td.rk i {{ display: block; font-style: normal; font-size: 10.5px;
  color: var(--muted); font-family: ui-monospace, monospace; }}
table.rounds td.rk.off b, table.rounds td.rk.off {{ color: var(--muted); }}
table.rounds .gapcol {{ width: 10px; padding: 0; border-bottom: 0;
  border-left: 1px solid var(--line); }}
table.rounds thead th {{ position: sticky; background: var(--surface); z-index: 1; }}
table.rounds thead tr:first-child th {{ top: 0; }}
table.rounds thead tr:last-child th {{ top: 26px; }}
</style></head>
<body>
<h1>Task seed classes — {LO}–{HI} in {STEP}s</h1>
<p class="sub">{html.escape(s['window']['start'][:19])} → {html.escape(s['window']['end'][:19])} ·
  source <code>/api/v3/tasks</code> · fetched {html.escape(s['generated_at'][:19])}</p>

<div class="legend">{''.join(
  f'<span><i class="dot" style="background:var(--c{i})"></i>{html.escape(ct)}</span>'
  for i, ct in enumerate(cell_types))}</div>

<section class="card">
  <div class="tiles">
    <div class="tile"><b>{n}</b><span>seeds binned</span></div>
    <div class="tile"><b>{tot['tasks_in_window']}</b><span>tasks in window</span></div>
    <div class="tile"><b>{tot['distinct_seeds']}</b><span>distinct seeds</span></div>
    <div class="tile"><b>{tot['seed_min']}</b><span>lowest seed</span></div>
    <div class="tile"><b>{tot['seed_max']}</b><span>highest seed</span></div>
    <div class="tile"><b>{len(BINS)}</b><span>classes</span></div>
    <div class="tile"><b>{expected:.1f}</b><span>expected per class</span></div>
    <div class="tile"><b>{sum(1 for i in range(len(BINS)) if counts[TOTAL][i] == 0)}</b>
      <span>empty classes</span></div>
    <div class="tile"><b>{max(counts[TOTAL]) if n else 0}</b><span>busiest class</span></div>
    <div class="tile"><b>{s['min_gap']['min'] if gaps else '—'}</b>
      <span>tightest pair</span></div>
    <div class="tile"><b>{s['min_gap']['median'] if gaps else '—'}</b>
      <span>median closest pair</span></div>
  </div>
  <p class="note">{html.escape(note_line)}</p>
</section>

<section class="card">
  <h2>Seeds per class, by cell type</h2>
  <table>
    <thead><tr><th>cell type</th>{''.join(f'<th class="c">{html.escape(l)}</th>' for l in labels)}
      <th class="num">total</th></tr></thead>
    <tbody>{matrix_rows}</tbody>
  </table>
  <p class="note">Shade depth tracks the count within its own row group (cell types share one
    scale, TOTAL has its own); “-” means no task fell in that class.</p>
</section>

<section class="card">
  <h2>Class totals against a uniform draw</h2>
  <div>{bars}</div>
  <p class="note">Dashed line marks {expected:.1f} — what each class would hold if the
    {n} seeds were spread evenly over {LO}–{HI}.</p>
</section>

<section class="card">
  <h2>Every seed drawn</h2>
  <div>{''.join(lanes)}</div>
  <div class="axis">{ticks}</div>
  <p class="note">One dot per task at its seed value; hover for the task. Class boundaries are the
    tick marks, so a cluster straddling a boundary is visible here but split in the matrix above.</p>
</section>

{gap_html}

<section class="card">
  <h2>Round by round — class counts so far against the seeds just drawn</h2>
  <p class="note lead">One table per cell type over the 3-seed tasks in
    {html.escape(rnd.get('window', {}).get('start', '')[:19])} →
    {html.escape(rnd.get('window', {}).get('end', '')[:19])}, oldest round first.
    Three blocks share every round row. The first is the running tally of how
    many seeds each of the nine classes had taken <i>before</i> this round. The
    second gives the rank that class held at that moment for each of this
    round's three seeds. All nine classes are ranked <b>1–9 with no rank
    repeated</b>: most-selected first, and where counts tie the class updated in
    the later round takes the lower rank number, falling back to the lower class
    when even that ties. Round 1 has nothing behind it to rank, so its three
    cells are <b>0</b>. The third block tallies how often each rank had been
    handed out over those same previous rounds: ten columns for ranks 1–9 plus
    one for the seeds that took no rank — round 1's three, and any seed outside
    {LO}–{HI}. All three blocks read off the same prior state, so a row checks
    against itself: the seed marked rank 1 belongs to the class holding the
    largest count in its own row, and the third block's row sums to three seeds
    times the rounds behind it. Because each round's ranks are a permutation of
    1–9, a generator that ignored the running tally would average exactly
    <b>5.0</b>.</p>
</section>

{''.join(round_cards)}

<section class="card">
  <h2>Seeds ({n})</h2>
  <div class="scroll"><table>
    <thead><tr><th>created</th><th>cell type</th><th class="num">seed</th><th>class</th>
      <th class="num">pos</th><th>task id</th></tr></thead>
    <tbody>{table_rows}</tbody>
  </table></div>
</section>
</body></html>
"""
    Path(path).write_text(doc)


if __name__ == "__main__":
    main()
