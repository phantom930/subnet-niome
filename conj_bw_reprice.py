#!/usr/bin/env python3
"""conj_bw_reprice.py — re-price a CR_ARMS sweep with the sampling noise removed.

`conj_replicate.py` estimates `v_clean`/`v_rest` by drawing NS seeds per regime from a SHARED
`random.Random` that advances between arms. So two arms are never scored on the same sampled
seeds, and arms whose constructions are byte-identical still price differently — measured on the
CD34+ band_width sweep, two arms at band 8 / clean 183 / wxfid 317.9 came out at own-field 0.00021
and 0.00009, a 2.3x gap from nothing at all.

That is fine for the conjunction-vs-all-HDR question the canonical runs ask (the effect there is
~1.8x and the two arms differ structurally), and it is fatal for ARM-vs-ARM sweeps, where the
constructions differ by ~1% and the noise is ~18%.

This re-prices every arm of a sweep against its own field using a COMMON `v_clean`/`v_rest` per
contract — the mean over that contract's arms — so the only quantities that still separate the arms
are the deterministic ones the build actually produced: `band`, `clean` and `wxfid`.

    python conj_bw_reprice.py conj_replicate_CD34_HSPC_bw.json
"""
import sys

sys.argv = sys.argv[:2] or ["x"]
sys.path.insert(0, "/root/workspace/subnet-niome")

import collections   # noqa: E402
import json          # noqa: E402
import math          # noqa: E402
import os            # noqa: E402
import statistics as st  # noqa: E402

os.environ.setdefault("NIOME_INSTANCE", "bwreprice")

from conj_stageb import price, OURS, API   # noqa: E402
from sd_task import fetch                  # noqa: E402

PATH = sys.argv[1] if len(sys.argv) > 1 else "conj_replicate_CD34_HSPC_bw.json"


def own_fields():
    """{task_id[:8]: sorted finals} — our own hotkeys removed, best row per miner."""
    sc = fetch(f"{API}/miners/scores?limit=40000", cache="sd_task_scores.json")
    sc = sc if isinstance(sc, list) else (sc.get("data") or sc.get("items") or [])
    best = collections.defaultdict(dict)
    for r in sc:
        if r.get("miner_hotkey") in OURS:
            continue
        hk, f = r["miner_hotkey"], float(r.get("final_score") or 0)
        if hk not in best[r["task_id"]] or f > best[r["task_id"]][hk]:
            best[r["task_id"]][hk] = f
    return {t[:8]: sorted(v.values(), reverse=True) for t, v in best.items()}


def main():
    d = json.load(open(PATH))
    fields = own_fields()
    widths = sorted({a["width"] for r in d for a in r["arms"]})

    print(f"{PATH}: {len(d)} contracts, widths {widths}\n")
    print("spread ACROSS ARMS within a contract (construction vs sampled):")
    for key in ("clean", "wxfid", "v_clean", "v_rest"):
        sp = [max(a[key] for a in r["arms"]) - min(a[key] for a in r["arms"]) for r in d]
        lvl = [st.mean([a[key] for a in r["arms"]]) for r in d]
        print(f"  {key:<8} mean spread {st.mean(sp):8.4f} = {st.mean(sp)/st.mean(lvl)*100:5.2f}% "
              f"of level")

    lo, hi = widths[0], widths[-1]
    diffs = [{a["width"]: a for a in r["arms"]}[lo]["v_clean"]
             - {a["width"]: a for a in r["arms"]}[hi]["v_clean"] for r in d]
    m, sd = st.mean(diffs), st.stdev(diffs)
    se = sd / math.sqrt(len(diffs))
    print(f"\npaired v_clean, w={lo} minus w={hi}: mean {m:+.4f} sd {sd:.4f} t = {m/se:.2f} "
          f"(df {len(diffs)-1}), positive {sum(1 for x in diffs if x > 0)}/{len(diffs)}")

    print("\nRE-PRICED with a common v_clean/v_rest per contract:")
    print(f"{'contract':>9} " + " ".join(f"{w:>10}" for w in widths) + f" {'all-HDR':>10}")
    tot = collections.defaultdict(float)
    ahtot = 0.0
    for r in d:
        a = {x["width"]: x for x in r["arms"]}
        vc = st.mean(x["v_clean"] for x in r["arms"])
        vr = st.mean(x["v_rest"] for x in r["arms"])
        fs = [fields[r["task"]]]
        row = []
        for w in widths:
            e, _ = price(a[w]["band"], a[w]["clean"] - a[w]["band"], vc, vr, a[w]["wxfid"], fs)
            tot[w] += e
            row.append(e)
        ah = r["all_hdr"]
        eah, _ = price(ah["band"], 0, 0.0, ah["floor"], ah["wxfid"], fs)
        ahtot += eah
        print(f"{r['task']:>9} " + " ".join(f"{e:>10.5f}" for e in row) + f" {eah:>10.5f}")
    print(f"\n{'aggregate':>9} " + " ".join(f"{tot[w]:>10.5f}" for w in widths)
          + f" {ahtot:>10.5f}")
    print(f"{'vs w=%d' % widths[0]:>9} "
          + " ".join(f"{tot[w]/tot[widths[0]]:>9.3f}x" for w in widths))
    print(f"{'vs allHDR':>9} " + " ".join(f"{tot[w]/ahtot:>9.2f}x" for w in widths))


if __name__ == "__main__":
    main()
