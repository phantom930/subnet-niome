#!/usr/bin/env python3
"""band_report.py — fold validator scores into band_efficiency.json and render band_efficiency.html.

Run band_scan.py first. This adds each hotkey's actual score/rank/payout for every scanned round,
then draws one 100-999 strip per round: nine window segments, every band seed ticked, the three
drawn seeds overlaid, hits marked. The point of the picture is that a hit or miss is visible
directly rather than inferred from a table of numbers.
"""
import json, html, urllib.request
from collections import defaultdict

DIST = [0.3, 0.2, 0.2, 0.15, 0.05, 0.03, 0.025, 0.02, 0.015, 0.01]
HOTKEY_SS58 = {
    "5HT66iVw1UPgQa73toQ3PhNKQ6FvL2z1NWk2dC1EdnX5wYHW": "h0",
    "5GFE8UJcTjEW7QsdHvQDUxbPsVKLfNzTcwGUPbf6Nc5o1hvb": "h1",
    "5FP4o2SSosZbCB71TzMQC2WPxTsKQUEzghWvkcc4B4PEbUeU": "h2",
    "5Fjzzbaf6q1fQfiprNKZtv8Twxm4J8C94aev4egCFYyrTXdf": "h3",
    "5GBWGSM6ZTk1hzrf3cgrAc9vmkA7p1x8oEVz6QLk6wACsuC2": "h4",
    "5GNz6g47q45YN471GpfdSrMSeJrzynbCgBBic759nsjLaPNn": "h5",
    "5H5v45M2i6cFtS3Di4abh2zjWrPmuJZzNPoAYtVAXVFa6sp3": "h6",
    "5CS8FdHr8Ddv14QHr75E5K7e9xv3wre9m66ywRHZMR8zhntN": "h7",
    "5Cr9gJ3ukDDdxnhpyRM58gdMz7Sjuw8u58iiuYphVz1ZUGD2": "h8",
}


def attach_scores(data):
    sc = json.load(urllib.request.urlopen(
        "https://niome-api.genomes.io/api/v3/miners/scores?limit=40000", timeout=120))
    sc = sc if isinstance(sc, list) else (sc.get("data") or sc.get("items") or [])
    by = defaultdict(list)
    for x in sc:
        by[x["task_id"]].append(x)
    for task in data["tasks"]:
        rows = by.get(task["task_id"], [])
        if not rows:
            task["scored"] = False
            continue
        best = {}
        for x in rows:
            k = x["miner_hotkey"]
            if k not in best or x["final_score"] > best[k]["final_score"]:
                best[k] = x
        ranked = sorted(best.values(), key=lambda y: -y["final_score"])
        task["scored"] = True
        task["field"] = {"miners": len(ranked), "rank1": ranked[0]["final_score"],
                         "rank10": ranked[9]["final_score"] if len(ranked) > 9 else 0.0}
        payout = 0.0
        for hk, tag in HOTKEY_SS58.items():
            x = best.get(hk)
            if not x or tag not in task["hotkeys"]:
                continue
            pos = next(i for i, y in enumerate(ranked, 1) if y["miner_hotkey"] == hk)
            pay = DIST[pos - 1] if pos <= 10 and x["final_score"] > 0 else 0.0
            payout += pay
            task["hotkeys"][tag].update(
                final=round(x["final_score"], 2), rank=pos, payout=pay,
                consistency=round(x["breakdown"]["consistency_factor"], 3),
                weighted=round(x["breakdown"]["total_weighted_score"], 1),
                fidelity=round(x["breakdown"]["distribution_fidelity_factor"], 3))
        task["fleet_payout"] = round(payout, 4)
    data["summary"]["fleet_payout_total"] = round(
        sum(t.get("fleet_payout", 0.0) for t in data["tasks"]), 4)
    data["summary"]["rounds_scored"] = sum(1 for t in data["tasks"] if t.get("scored"))
    return data


def strip(task):
    """One 100-999 SVG strip: window segments, band ticks, drawn seeds, hits."""
    W, H, PAD = 1000, 132, 8
    def x(seed):
        return PAD + (seed - 100) / 899.0 * (W - 2 * PAD)
    parts = [f'<svg viewBox="0 0 {W} {H}" class="strip" role="img" '
             f'aria-label="band coverage for {html.escape(task["task_id"][:8])}">']
    for i in range(9):
        lo = 100 + i * 100
        x0, x1 = x(lo), x(min(lo + 99, 999))
        parts.append(f'<rect x="{x0:.1f}" y="30" width="{x1-x0:.1f}" height="46" '
                     f'class="win {"win-alt" if i % 2 else ""}"/>')
        parts.append(f'<text x="{(x0+x1)/2:.1f}" y="94" class="wlab">h{i}</text>')
        parts.append(f'<text x="{(x0+x1)/2:.1f}" y="107" class="wsub">{lo}-{lo+99}</text>')
    for tag, hk in task["hotkeys"].items():
        for b in hk["band"]:
            parts.append(f'<line x1="{x(b):.2f}" y1="34" x2="{x(b):.2f}" y2="72" class="band"/>')
    for s in task["seeds"]:
        owner = next((t for t, hk in task["hotkeys"].items()
                      if hk["window"][0] <= s <= hk["window"][1]), None)
        hit = bool(owner and s in task["hotkeys"][owner]["hits"])
        cls = "hit" if hit else "miss"
        parts.append(f'<line x1="{x(s):.2f}" y1="18" x2="{x(s):.2f}" y2="80" class="seed {cls}"/>')
        parts.append(f'<circle cx="{x(s):.2f}" cy="14" r="5" class="dot {cls}"/>')
        parts.append(f'<text x="{x(s):.2f}" y="127" class="slab {cls}">{s}</text>')
    parts.append("</svg>")
    return "".join(parts)


def render(data):
    s = data["summary"]
    rows = []
    for t in data["tasks"]:
        hit_hk = [f'{tag}&nbsp;<span class="mono">{hk["hits"][0]}</span>'
                  for tag, hk in sorted(t["hotkeys"].items()) if hk["hits"]]
        owners = []
        for seed in t["seeds"]:
            tag = next((k for k, hk in t["hotkeys"].items()
                        if hk["window"][0] <= seed <= hk["window"][1]), "—")
            hk = t["hotkeys"].get(tag, {})
            got = seed in hk.get("hits", [])
            # band_scan stores nearest as {seed_str: value}; band_rebuild stores a bare int.
            nz = hk.get("nearest")
            near = nz.get(str(seed)) if isinstance(nz, dict) else nz
            d = abs(near - seed) if near is not None else None
            owners.append(
                f'<tr><td class="mono">{seed}</td><td>{tag}</td>'
                f'<td class="{"y" if got else "n"}">{"HIT" if got else "miss"}</td>'
                f'<td class="mono dim">{"—" if got else f"nearest {near} (&Delta;{d})"}</td>'
                f'<td class="mono">{hk.get("consistency","—")}</td>'
                f'<td class="mono">{hk.get("final","—")}</td>'
                f'<td class="mono">{hk.get("rank","—")}</td></tr>')
        field = t.get("field") or {}
        rows.append(f'''<section class="round">
  <header>
    <div><h3>{html.escape(t["cell_type"] or "?")} <span class="dim mono">{t["task_id"][:8]}</span>
      <span class="tag {t.get("source","rebuilt")}">{t.get("source","rebuilt")}</span></h3>
      <div class="dim">{html.escape(t["created_at"][:16].replace("T"," "))} &middot; seeds
        <span class="mono">{", ".join(str(x) for x in t["seeds"])}</span></div></div>
    <div class="tally">
      <div class="big {"y" if t["hits"] else "n"}">{t["hits"]}<span class="dim">/{t["opportunities"]}</span></div>
      <div class="dim">hits</div>
    </div>
  </header>
  {strip(t)}
  <div class="cols">
    <table class="seeds"><thead><tr><th>seed</th><th>owner</th><th>result</th><th>detail</th>
      <th>cons</th><th>final</th><th>rank</th></tr></thead><tbody>{"".join(owners)}</tbody></table>
    <div class="side">
      <div><span class="dim">field rank&nbsp;1</span> <b class="mono">{field.get("rank1","—") if not field else f'{field["rank1"]:.1f}'}</b></div>
      <div><span class="dim">rank&nbsp;10 cutoff</span> <b class="mono">{f'{field["rank10"]:.1f}' if field else "—"}</b></div>
      <div><span class="dim">fleet payout</span> <b class="mono {"y" if t.get("fleet_payout") else ""}">{t.get("fleet_payout", 0):.3f}</b></div>
      <div><span class="dim">hit hotkeys</span> <b>{", ".join(hit_hk) if hit_hk else "—"}</b></div>
    </div>
  </div>
</section>''')

    band_rows = []
    agg = defaultdict(lambda: [0, 0, 0])
    for t in data["tasks"]:
        for tag, hk in t["hotkeys"].items():
            agg[tag][0] += hk["band_size"]
            agg[tag][1] += len(hk["owned_seeds"])
            agg[tag][2] += len(hk["hits"])
    n = len(data["tasks"])
    for tag in sorted(agg):
        tot, opp, hits = agg[tag]
        band_rows.append(f'<tr><td>{tag}</td><td class="mono">{data["windows"][tag][0]}-'
                         f'{data["windows"][tag][1]}</td><td class="mono">{tot/n:.1f}</td>'
                         f'<td class="mono">{opp}</td><td class="mono {"y" if hits else ""}">{hits}</td></tr>')

    return f'''<title>Band Efficiency</title>
<style>
:root {{
  --bg:#fbfaf8; --panel:#fff; --ink:#1c1a17; --dim:#6b6560; --line:#e7e2db;
  --band:#c9a227; --hit:#1f9d55; --miss:#c2410c; --win:#f3efe9; --win2:#eee9e1; --accent:#2563eb;
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --bg:#14130f; --panel:#1c1a16; --ink:#ece7de; --dim:#9a938a; --line:#2e2a24;
    --band:#e0b93c; --hit:#3fbc76; --miss:#f0803c; --win:#232019; --win2:#1e1b15; --accent:#6ea8fe;
  }}
}}
:root[data-theme="dark"] {{
  --bg:#14130f; --panel:#1c1a16; --ink:#ece7de; --dim:#9a938a; --line:#2e2a24;
  --band:#e0b93c; --hit:#3fbc76; --miss:#f0803c; --win:#232019; --win2:#1e1b15; --accent:#6ea8fe;
}}
body {{ background:var(--bg); color:var(--ink); font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif; margin:0; padding:32px 20px 64px; }}
.wrap {{ max-width:1060px; margin:0 auto; }}
h1 {{ font-size:26px; margin:0 0 4px; letter-spacing:-.02em; }}
h3 {{ font-size:16px; margin:0 0 2px; }}
.dim {{ color:var(--dim); }}
.mono {{ font-family:ui-monospace,SFMono-Regular,Menlo,monospace; font-variant-numeric:tabular-nums; }}
.y {{ color:var(--hit); font-weight:600; }} .n {{ color:var(--miss); }}
.lede {{ color:var(--dim); max-width:70ch; margin:0 0 24px; }}
.cards {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr)); gap:12px; margin:0 0 28px; }}
.card {{ background:var(--panel); border:1px solid var(--line); border-radius:10px; padding:14px 16px; }}
.card .v {{ font-size:24px; font-weight:650; letter-spacing:-.02em; }}
.card .k {{ color:var(--dim); font-size:12px; text-transform:uppercase; letter-spacing:.06em; }}
.round {{ background:var(--panel); border:1px solid var(--line); border-radius:12px; padding:18px 20px; margin:0 0 18px; }}
.round header {{ display:flex; justify-content:space-between; align-items:flex-start; gap:16px; margin-bottom:6px; }}
.tally {{ text-align:right; }} .tally .big {{ font-size:26px; font-weight:650; }}
.strip {{ width:100%; height:auto; display:block; margin:6px 0 10px; overflow:visible; }}
.win {{ fill:var(--win); }} .win-alt {{ fill:var(--win2); }}
.band {{ stroke:var(--band); stroke-width:2.2; opacity:.85; }}
.seed {{ stroke-width:2; }} .seed.hit {{ stroke:var(--hit); }} .seed.miss {{ stroke:var(--miss); }}
.dot.hit {{ fill:var(--hit); }} .dot.miss {{ fill:var(--miss); }}
.wlab {{ fill:var(--dim); font-size:11px; text-anchor:middle; font-weight:600; }}
.wsub {{ fill:var(--dim); font-size:9px; text-anchor:middle; opacity:.75; }}
.slab {{ font-size:10px; text-anchor:middle; font-family:ui-monospace,monospace; }}
.slab.hit {{ fill:var(--hit); font-weight:700; }} .slab.miss {{ fill:var(--miss); }}
.cols {{ display:grid; grid-template-columns:1fr 210px; gap:20px; align-items:start; }}
@media (max-width:760px) {{ .cols {{ grid-template-columns:1fr; }} }}
table {{ border-collapse:collapse; width:100%; font-size:13px; }}
th {{ text-align:left; color:var(--dim); font-weight:600; font-size:11px; text-transform:uppercase; letter-spacing:.05em; border-bottom:1px solid var(--line); padding:4px 8px 4px 0; }}
td {{ padding:5px 8px 5px 0; border-bottom:1px solid var(--line); }}
.side div {{ display:flex; justify-content:space-between; gap:10px; padding:5px 0; border-bottom:1px solid var(--line); }}
.note {{ background:var(--panel); border:1px solid var(--line); border-left:3px solid var(--accent); border-radius:8px; padding:14px 18px; margin:26px 0 0; color:var(--dim); }}
.note b {{ color:var(--ink); }}
.legend {{ display:flex; gap:18px; flex-wrap:wrap; margin:0 0 22px; color:var(--dim); font-size:13px; }}
.tag {{ font-size:10px; text-transform:uppercase; letter-spacing:.06em; padding:2px 7px; border-radius:99px; vertical-align:2px; margin-left:6px; }}
.tag.archived {{ background:var(--hit); color:var(--bg); }}
.tag.rebuilt {{ background:var(--line); color:var(--dim); }}
.legend i {{ display:inline-block; width:12px; height:12px; border-radius:2px; vertical-align:-1px; margin-right:6px; }}
</style>
<div class="wrap">
<h1>Band efficiency</h1>
<p class="lede">Every hotkey pins an all-HDR clean band inside its own 100-seed window. A round draws
three seeds; each falls in exactly one hotkey's window, and that hotkey scores
<span class="mono">consistency 1.000</span> on it only if the seed is one of its band seeds. This is
what actually happened, per round.</p>

<div class="cards">
  <div class="card"><div class="v">{s["hits"]}<span class="dim">/{s["opportunities"]}</span></div><div class="k">seeds hit</div></div>
  <div class="card"><div class="v">{s["hit_rate"]*100:.1f}%</div><div class="k">observed rate</div></div>
  <div class="card"><div class="v">{s["expected_hit_rate"]*100:.1f}%</div><div class="k">expected (band/100)</div></div>
  <div class="card"><div class="v">{s["mean_band_size"]}</div><div class="k">mean band size</div></div>
  <div class="card"><div class="v">{s["rounds_with_a_hit"]}<span class="dim">/{data["tasks_scanned"]}</span></div><div class="k">rounds with a hit</div></div>
  <div class="card"><div class="v">{s.get("fleet_payout_total",0):.3f}</div><div class="k">payout earned</div></div>
  <div class="card"><div class="v">{s["rounds"]}<span class="dim">/{data["tasks_in_range"]}</span></div><div class="k">rounds covered</div></div>
</div>

<div class="legend">
  <span><i style="background:var(--band)"></i>band seed (consistency 1.000 if drawn)</span>
  <span><i style="background:var(--hit)"></i>drawn seed &mdash; hit</span>
  <span><i style="background:var(--miss)"></i>drawn seed &mdash; miss</span>
</div>

{"".join(rows)}

<h3 style="margin:28px 0 8px">Per cell type</h3>
<table><thead><tr><th>cell type</th><th>mean band</th><th>seeds owned</th><th>hits</th><th>rate</th></tr></thead>
<tbody>{"".join(f'<tr><td>{html.escape(c)}</td><td class="mono">{v["mean_band"]}</td>'
                f'<td class="mono">{v["opportunities"]}</td><td class="mono y">{v["hits"]}</td>'
                f'<td class="mono">{v["hit_rate"]*100:.1f}%</td></tr>'
                for c, v in s.get("by_cell_type", {}).items())}</tbody></table>

<h3 style="margin:28px 0 8px">Per hotkey</h3>
<table><thead><tr><th>hotkey</th><th>window</th><th>mean band</th><th>seeds owned</th><th>hits</th></tr></thead>
<tbody>{"".join(band_rows)}</tbody></table>

<div class="note">
<b>Read the sample size before the rate.</b> This covers {data["tasks_scanned"]} of the
{data["tasks_in_range"]} three-seed tasks in range &mdash; the fleet only began archiving on
2026-09-01, and rebuilding the rest would use a config that has since changed
(<span class="mono">group_size</span> 100&nbsp;&rarr;&nbsp;80,
<span class="mono">light_cell_rows</span> off&nbsp;&rarr;&nbsp;6), so their bands would not be the
ones we submitted. {s["opportunities"]} drawn seeds cannot separate a {s["hit_rate"]*100:.0f}% hit
rate from {s["expected_hit_rate"]*100:.0f}% &mdash; the interval is wide. What the data does show
cleanly is the mechanism: a hit yields round consistency &asymp;0.40
<span class="mono">(1.0&nbsp;+&nbsp;0.1&nbsp;+&nbsp;0.1)/3</span>, and every miss is a miss by
chance, not by proximity &mdash; band membership comes from
<span class="mono">sha256(seed|design)</span>, so seed 159 and 160 are unrelated draws.
</div>
</div>'''


if __name__ == "__main__":
    data = attach_scores(json.load(open("band_efficiency.json")))
    json.dump(data, open("band_efficiency.json", "w"), indent=2)
    open("band_efficiency.html", "w").write(render(data))
    print(f"band_efficiency.json  {len(json.dumps(data))/1024:.0f} KB")
    print(f"band_efficiency.html  {len(open('band_efficiency.html').read())/1024:.0f} KB")
    print(f"  {data['summary']['hits']}/{data['summary']['opportunities']} hits, "
          f"payout {data['summary'].get('fleet_payout_total')}")
