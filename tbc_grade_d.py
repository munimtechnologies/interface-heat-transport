# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Made with AI agents working in MV2's name.
"""tbc_grade_d: grade the blind test D1-D4 of PLAN-amendment-0.3.3.md (section 5).

Sergeev's direct electron-to-substrate channel depends only on the metal and grows as s T. Given a
measured G(T1) on an interface whose conductance stays flat, the least G(T2) the circuit allows
(phonon side not falling between T1 and T2) is tbc_direct.min_rise_prediction. A channel is
EXCLUDED when the measured G(T2) < 0.85 G_min(T2) (threshold fixed in amendment 0.3.3).

Also (exploratory, not pre-registered): the largest slope s each flat interface allows, with the
0.85 margin and at central values, and the phonon side of Li et al.'s Sample 1 left over after
removing the largest interface-independent channel the flat interfaces allow.

python tbc_grade_d.py --read data/read-d1d3.json --meas data/meas-li2022.json --out results/dtest
Writes grade-d1d4.json; the last printed line starts with 'TBDG:'.
"""
import argparse
import json
import os

import tbc_direct as td
from run_tbc_direct import AL, GEP400, N, TEMPS

VERSION = "mv2 tbc-grade-d 0.1.0"
MW = 1e6
CUT = 0.85          # amendment 0.3.3: excluded if measured G(T2) < 0.85 G_min(T2)
MIN_SPAN = 100.0    # D3 graded only if T2 - T1 >= 100 K


def slopes():
    """Direct-channel slopes s (W m^-2 K^-2) named in amendment 0.3.3."""
    a = AL
    return {
        "sergeev_g_DFT": td.sergeev_from_g(a["gamma_vol"], a["g_DFT"], a["u_L"], a["u_T"], a["n_e"], 300.0) / 300.0,
        "sergeev_table1": 264e6 / 300.0,
        "sergeev_lambda": td.sergeev_from_lambda(a["gamma_vol"], a["lam_w2_meV2"], a["u_L"], a["u_T"], a["n_e"],
                                                 300.0) / 300.0,
        "sergeev_eq4_printed": td.sergeev_from_lambda(a["gamma_vol"], a["lam_w2_meV2"], a["u_L"], a["u_T"],
                                                      a["n_e"], 300.0, zeta=False) / 300.0,
        "DSF_107_linear": 107e6 / 300.0,
    }


def g_min(g1, t1, t2, s):
    """Least G(T2) (W m^-2 K^-1) allowed with slope s; None if s T1 >= G(T1)."""
    return td.min_rise_prediction(g1, t1, t2, s, GEP400, N)["G_min"]


def grade(g1, t1, g2, t2, s, cut=CUT):
    """One pre-registered comparison."""
    gm = g_min(g1, t1, t2, s)
    if gm is None:
        return {"G_min": None, "excluded": True, "reason": "channel alone exceeds G(T1)"}
    return {"G_min": gm, "threshold": cut * gm, "ratio": g2 / gm, "excluded": g2 < cut * gm}


def s_allowed(g1, t1, g2, t2, cut=1.0, iters=200):
    """Largest s with cut * G_min(T2; s) <= G(T2). G_min grows monotonically with s (amendment 0.4,
    section 2), so bisection is exact. Returns 0.0 if even s = 0 fails, None if no bound."""
    lo, hi = 0.0, g1 / t1 * (1.0 - 1e-12)
    if cut * g_min(g1, t1, t2, lo) > g2:
        return 0.0
    if cut * g_min(g1, t1, t2, hi) <= g2:
        return None
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if cut * g_min(g1, t1, t2, mid) <= g2:
            lo = mid
        else:
            hi = mid
    return lo


def _mw(x):
    return None if x is None else x / MW


def _sample(read, key):
    d = read[key]
    gs = {float(t): float(v) * MW for t, v in d["G"].items()}
    return gs, float(d["reading_uncertainty"]) * MW


def run(read, meas):
    sl = slopes()
    li, u_li = _sample(read, "li2022_sample2")
    ch, u_ch = _sample(read, "cheng2020_al_sapphire")
    t2_li = float(read["li2022_sample2"]["highest_T"])
    t2_ch = float(read["cheng2020_al_sapphire"]["highest_T"])
    cases = {"D1_D2_li_sample2": (li[400.0], 400.0, li[t2_li], t2_li, u_li),
             "D3_cheng_al_sapphire": (ch[300.0], 300.0, ch[t2_ch], t2_ch, u_ch)}
    out = {"version": VERSION, "cut": CUT, "slopes_MW_per_K": {k: v / MW for k, v in sl.items()},
           "G_ep_MW": {"G400": GEP400 / MW, "n": N}, "cases": {}}
    for name, (g1, t1, g2, t2, u) in cases.items():
        c = {"T1": t1, "T2": t2, "G1_MW": g1 / MW, "G2_MW": g2 / MW, "span_K": t2 - t1,
             "graded": t2 - t1 >= MIN_SPAN, "by_slope": {}, "worst_case_reading": {}}
        for k, s in sl.items():
            r = grade(g1, t1, g2, t2, s)
            c["by_slope"][k] = {"G_min_MW": _mw(r["G_min"]), "threshold_MW": _mw(r.get("threshold")),
                                "ratio": r.get("ratio"), "excluded": r["excluded"],
                                "sigma_e_T1_MW": s * t1 / MW, "note": r.get("reason", "")}
            # reading corners least favourable to the verdict just found
            if r["excluded"]:
                w1, w2 = g1 - u, g2 + u
            else:
                w1, w2 = g1 + u, g2 - u
            w = grade(w1, t1, w2, t2, s)
            c["worst_case_reading"][k] = {"G1_MW": w1 / MW, "G2_MW": w2 / MW, "G_min_MW": _mw(w["G_min"]),
                                          "excluded": w["excluded"],
                                          "verdict_holds": w["excluded"] == r["excluded"]}
        c["s_allowed_cut085_MW_per_K"] = s_allowed(g1, t1, g2, t2, CUT) / MW
        c["s_allowed_central_MW_per_K"] = s_allowed(g1, t1, g2, t2, 1.0) / MW
        c["sigma_e_300_allowed_central_MW"] = c["s_allowed_central_MW_per_K"] * 300.0
        out["cases"][name] = c

    d12 = out["cases"]["D1_D2_li_sample2"]
    d3 = out["cases"]["D3_cheng_al_sapphire"]
    grades = {
        "D1": "EXCLUDED" if d12["by_slope"]["sergeev_g_DFT"]["excluded"] else "NOT EXCLUDED",
        "D2": "EXCLUDED" if d12["by_slope"]["DSF_107_linear"]["excluded"] else "NOT EXCLUDED",
        "D3_sergeev_g_DFT": ("UNDETERMINED" if not d3["graded"] else
                             "EXCLUDED" if d3["by_slope"]["sergeev_g_DFT"]["excluded"] else "NOT EXCLUDED"),
        "D3_DSF_107_linear": ("UNDETERMINED" if not d3["graded"] else
                              "EXCLUDED" if d3["by_slope"]["DSF_107_linear"]["excluded"] else "NOT EXCLUDED"),
    }
    grades["D4"] = ("TRIGGERED: no bulk-only direct channel of Sergeev's size at Al interfaces"
                    if grades["D1"] == "EXCLUDED" and grades["D3_sergeev_g_DFT"] == "EXCLUDED"
                    else "NOT TRIGGERED")
    grades["all_sergeev_versions_excluded_both"] = all(
        out["cases"][c]["by_slope"][k]["excluded"] for c in out["cases"] for k in sl if k.startswith("sergeev"))
    out["grades"] = grades

    # exploratory: phonon side of Sample 1 after the largest interface-independent channel allowed
    s_cap = min(d12["s_allowed_central_MW_per_K"], d3["s_allowed_central_MW_per_K"]) * MW
    s1 = {k: v for k, v in meas["Al/Si"].items()}
    inf = td.infer_with_direct(s1, GEP400, N, s_cap, TEMPS)
    gpp = [None if x is None else x / MW for x in inf["G_pp"]]
    out["sample1_phonon_side_after_cap"] = {"s_cap_MW_per_K": s_cap / MW, "T": inf["T"], "G_pp_MW": gpp,
                                            "ratio_700_400": gpp[inf["T"].index(700)] / gpp[inf["T"].index(400)]}
    out["cheng_ep_share_300K"] = ch[300.0] / (GEP400 * (400.0 / 300.0) ** N)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--read", required=True)
    ap.add_argument("--meas", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    with open(a.read) as f:
        read = json.load(f)
    with open(a.meas) as f:
        meas = json.load(f)
    res = run(read, meas)
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "grade-d1d4.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))
    g = res["grades"]
    c1, c3 = res["cases"]["D1_D2_li_sample2"], res["cases"]["D3_cheng_al_sapphire"]
    print("TBDG: D1 %s (G_min %.1f, meas %.0f) | D2 %s (G_min %.1f) | D3 sergeev %s (G_min %.1f, meas %.0f) DSF %s "
          "(G_min %.1f) | D4 %s | s_central S2 %.4f sapphire %.4f | S1 G_pp 700/400 %.3f"
          % (g["D1"], c1["by_slope"]["sergeev_g_DFT"]["G_min_MW"], c1["G2_MW"], g["D2"],
             c1["by_slope"]["DSF_107_linear"]["G_min_MW"], g["D3_sergeev_g_DFT"],
             c3["by_slope"]["sergeev_g_DFT"]["G_min_MW"], c3["G2_MW"], g["D3_DSF_107_linear"],
             c3["by_slope"]["DSF_107_linear"]["G_min_MW"], g["D4"][:9], c1["s_allowed_central_MW_per_K"],
             c3["s_allowed_central_MW_per_K"], res["sample1_phonon_side_after_cap"]["ratio_700_400"]))


if __name__ == "__main__":
    main()
