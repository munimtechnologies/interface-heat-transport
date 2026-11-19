# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Made with AI agents working in MV2's name.
"""Phonon-side conductance implied by the Li et al. (2022) measurements under three G_ep(T) laws
(PLAN-amendment-0.3.2.md). Writes <out>/gep-li2022.json and <out>/run.txt; last line 'TBG:'."""
import argparse
import json
import os

import tbc_gep as tg

# Li et al.'s DMM curve in Fig. 1, read by MV2's hands on 2026-11-18 (MW m^-2 K^-1).
LI_DMM = {"Al/Si": {"300": 460e6, "400": 475e6, "700": 490e6},
          "Al/GaN": {"300": 440e6, "400": 455e6, "700": 470e6}}
LAWS = [("constant 1.21 GW (amendment 0.3.1)", 1.21e9, 0.0),
        ("T^-1/2 from 1.26 GW at 400 K (Li et al.'s implied; bulk-input expectation)", 1.26e9, 0.5),
        ("T^-1/4 from 1.26 GW at 400 K (half the bulk T-dependence)", 1.26e9, 0.25)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meas", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    meas = json.load(open(a.meas))
    res = {"tool": tg.VERSION, "made_by": "MV2 of Munim, Inc. (made with AI agents)",
           "source": meas.get("source", ""), "laws": []}
    lines = []
    for name, g400, n in LAWS:
        entry = {"law": name, "gep400": g400, "n": n, "pairs": {}}
        for pair in ("Al/Si", "Al/GaN"):
            r = tg.infer(meas[pair], g400, n)
            dmm = {}
            for t, gd in LI_DMM[pair].items():
                gep = tg.gep_power_law(g400, float(t), 400.0, n)
                dmm[t] = {"dmm": gd, "dmm_plus_ep": tg.series(gd, gep),
                          "meas": float(meas[pair][t])}
            r["li_dmm_series"] = dmm
            entry["pairs"][pair] = r
            lines.append("%s | %s | G_pp MW: %s | convex %s | E(400-700) %.1f MW" % (
                name, pair, " ".join("%.1f" % (x / 1e6) for x in r["G_pp"]),
                r["convex_everywhere"], r["linear_400_700"]["E"] / 1e6))
        res["laws"].append(entry)
    json.dump(res, open(os.path.join(a.out, "gep-li2022.json"), "w"), indent=1)
    half = res["laws"][1]["pairs"]
    last = "TBG: T^-1/2 G_pp Al/Si %s | Al/GaN %s | convex %s/%s" % (
        "/".join("%.1f" % (x / 1e6) for x in half["Al/Si"]["G_pp"]),
        "/".join("%.1f" % (x / 1e6) for x in half["Al/GaN"]["G_pp"]),
        half["Al/Si"]["convex_everywhere"], half["Al/GaN"]["convex_everywhere"])
    lines.append(last)
    with open(os.path.join(a.out, "run.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
