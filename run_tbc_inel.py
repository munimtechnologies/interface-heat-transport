# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""run_tbc_inel: stage 2b of the interface-heat-transport study (PLAN-amendment-0.3.md, with the
solver correction and recorded diagnostics of PLAN-amendment-0.3.1.md).

Inputs (all published in the repo before this runs): stage 2a's per-material lifetimes
results/anh/anh-<model>-<mat>.npz and results/anh/tba-summary.json (elastic two-stream G_BTE, relaxed
lattice constants), and the measured conductances read from Li et al. (2022) Fig. 1 in
data/meas-li2022.json. Computes the interfacial three-phonon channel G_inel(T) of tbc_inel for
Al/Si, Al/GaN and AlN/GaN, adds it in parallel to G_BTE (Al electron-phonon resistance in series),
and grades I0-I4. Last line of run.txt starts with 'TBI:'.
"""
import argparse
import json
import os

import numpy as np

import tbc_inel as ti

VERSION = "mv2 run-tbc-inel 0.1.1"

PAIRS = [("Al", "Si"), ("Al", "GaN"), ("AlN", "GaN")]
GEOM = {"Al": "fcc111", "Si": "diamond111", "C": "diamond111", "GaN": "wurtzite0001",
        "AlN": "wurtzite0001"}
NCONV = {"Al": 4, "Si": 8, "C": 8, "GaN": 4, "AlN": 4}   # atoms in the cell whose lengths are stored
G_EP_AL = float(np.sqrt(2.45e17 * 6.0))                    # stage 1/2a Al electron-phonon term
MEAS_300 = {"Al/Si": 379e6, "Al/GaN": 423e6}
ALNGAN_RANGE = (0.8 * 320e6, 1.2 * 490e6)
S_HI = 1000.0


def load_side(anh_dir, model, name, cell_lengths):
    z = np.load(os.path.join(anh_dir, f"anh-{model}-{name}.npz"))
    L = [x * 1e-10 for x in cell_lengths]
    if GEOM[name] == "wurtzite0001":
        a, c = L[0], L[2]
        vconv = np.sqrt(3.0) / 2.0 * a * a * c
    else:
        a, c = L[0], None
        vconv = a ** 3
    nb = z["freq"].shape[1]
    vcell = vconv * (nb / 3.0) / NCONV[name]
    return {"freq": z["freq"], "weights": z["weights"], "gamma3": z["gamma3"],
            "temps": [float(t) for t in z["temps"]], "vcell": vcell,
            "t": ti.t_eff(GEOM[name], a, c), "a": a, "c": c}


def side_at(s, T):
    it = s["temps"].index(T)
    return {"freq": s["freq"], "weights": s["weights"], "gamma": s["gamma3"][it],
            "vcell": s["vcell"], "t": s["t"]}


def total(pair, g_bte, g_in, s=1.0):
    g = g_bte + s * g_in
    return ti.series(G_EP_AL, g) if pair.startswith("Al/") else g


def phonon_side(g_meas, g_ep=G_EP_AL):
    """Phonon-side conductance implied by a measured total with g_ep in series."""
    return 1.0 / (1.0 / g_meas - 1.0 / g_ep)


def classify(vals):
    """Reading of PLAN-amendment-0.3 sec. 4 on a dict of scale factors."""
    lo, hi = min(vals.values()), max(vals.values())
    return ("too-weak" if lo > 8.0 else "too-strong" if hi < 0.25 else
            "right-order" if 0.5 <= lo and hi <= 2.0 else "intermediate")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--anh-dir", required=True)
    ap.add_argument("--meas", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    log(f"{VERSION} | {ti.VERSION}")
    tba = json.load(open(os.path.join(args.anh_dir, "tba-summary.json")))
    meas = json.load(open(args.meas))
    models = tba["models"]
    prim = models[0]
    res = {}
    for m in models:
        res[m] = {}
        for A, B in PAIRS:
            pair = f"{A}/{B}"
            key = f"oriented {pair}"
            try:
                sa = load_side(args.anh_dir, m, A, tba["info"][m][A]["info"]["cell_lengths"])
                sb = load_side(args.anh_dir, m, B, tba["info"][m][B]["info"]["cell_lengths"])
                rows_bte = tba["tables"][m][key]
            except (KeyError, FileNotFoundError, OSError) as e:
                log(f"{m} {pair}: SKIPPED ({type(e).__name__}: {e})")
                continue
            caches = {"classical": {}, "jdos": {}}
            rows = {}
            for T in sa["temps"]:
                k = str(int(T))
                if k not in rows_bte:
                    continue
                a_T, b_T = side_at(sa, T), side_at(sb, T)
                prim_r = ti.g_inel(a_T, b_T, T, kernel="classical", q_cache=caches["classical"])
                jd = ti.g_inel(a_T, b_T, T, kernel="jdos", q_cache=caches["jdos"])
                cl = ti.g_inel(a_T, b_T, T, kernel="classical", classical_c=True,
                               q_cache=caches["classical"])
                gb = float(rows_bte[k]["iso"]["G"])
                rows[k] = {"G_bte": gb, "G_inel": prim_r["G"], "G_inel_jdos": jd["G"],
                           "G_inel_classical_c": cl["G"], "parts": [prim_r["G_side1"], prim_r["G_side2"]],
                           "Qbar": prim_r["Qbar"], "rate_mean": prim_r["rate_mean"], "C": prim_r["C"],
                           "zero_nodes": prim_r["zero_nodes"], "capped_nodes": prim_r["capped_nodes"],
                           "G_total": total(pair, gb, prim_r["G"]),
                           "G_total_elastic": total(pair, gb, 0.0)}
            res[m][pair] = {"t_eff": [sa["t"], sb["t"]], "rows": rows}
            for k in ("300", "400", "700"):
                if k in rows:
                    r = rows[k]
                    log(f"{m} {pair} {k} K: G_bte {r['G_bte']/1e6:.0f} G_inel {r['G_inel']/1e6:.1f} "
                        f"({r['parts'][0]/1e6:.1f}+{r['parts'][1]/1e6:.1f}; jdos {r['G_inel_jdos']/1e6:.1f}; "
                        f"classical-c {r['G_inel_classical_c']/1e6:.1f}) Qbar {r['Qbar'][0]:.3f}/{r['Qbar'][1]:.3f} "
                        f"total {r['G_total']/1e6:.0f} (elastic only {r['G_total_elastic']/1e6:.0f}) MW/m2K")

    g = {}
    R = res.get(prim, {})

    def rr(pair, k, f):
        return R[pair]["rows"][k][f]
    try:
        rat = {p: rr(p, "700", "G_inel_classical_c") / rr(p, "400", "G_inel_classical_c")
               for p in ("Al/Si", "Al/GaN")}
        g["I0"] = ("PASS" if all(1.55 <= v <= 1.80 for v in rat.values()) else "FAIL") + \
            "(" + ",".join(f"{p}:{v:.3f}" for p, v in rat.items()) + ")"
    except KeyError as e:
        g["I0"] = f"UNDETERMINED(missing {e})"
    sstar = {}
    detail = {}
    for p in ("Al/Si", "Al/GaN"):
        try:
            el = (rr(p, "400", "G_bte"), rr(p, "700", "G_bte"))
            inn = (rr(p, "400", "G_inel"), rr(p, "700", "G_inel"))
            m400, m700 = float(meas[p]["400"]), float(meas[p]["700"])
        except KeyError:
            continue
        rise = m700 - m400
        sstar[p] = ti.solve_scale(el, inn, rise, g_ep=G_EP_AL, hi=S_HI)
        # recorded diagnostics (PLAN-amendment-0.3.1 sec. 3)
        s_max, r_max = ti.max_rise(el, inn, g_ep=G_EP_AL, hi=S_HI)
        gpp = (phonon_side(m400), phonon_side(m700))
        d_in = inn[1] - inn[0]
        s_pp = ((gpp[1] - gpp[0]) - (el[1] - el[0])) / d_in if d_in > 0 else float("inf")
        s_abs400 = (gpp[0] - el[0]) / inn[0] if inn[0] > 0 else float("inf")
        ep_t = (G_EP_AL, G_EP_AL * np.sqrt(400.0 / 700.0))
        s_ep_t = ti.solve_scale(el, inn, rise, g_ep=ep_t, hi=S_HI)
        detail[p] = {"rise_meas": rise, "max_model_rise": r_max, "s_at_max_rise": s_max,
                     "G_pp_meas": list(gpp), "s_pp": s_pp, "s_abs400": s_abs400,
                     "s_star_ep_falling": s_ep_t,
                     "no_fit": bool((not np.isfinite(sstar[p])) and s_max < 0.999 * S_HI)}
        log(f"{p}: s* {sstar[p]:.3g} | max model rise {r_max/1e6:.0f} at s {s_max:.3g} vs measured "
            f"{rise/1e6:.0f} | G_pp meas {gpp[0]/1e6:.0f}->{gpp[1]/1e6:.0f} | s_pp {s_pp:.3g} | "
            f"s_abs400 {s_abs400:.3g} | s* with falling G_ep {s_ep_t:.3g}")
    if len(sstar) == 2:
        g["I1"] = ("PASS" if all(0.5 <= v <= 2.0 for v in sstar.values()) else "FAIL") + \
            "(" + ",".join(f"s*{p}:{v:.2f}" for p, v in sstar.items()) + ")"
        s = sstar["Al/Si"]
        if np.isfinite(s):
            pred = total("Al/GaN", rr("Al/GaN", "700", "G_bte"), rr("Al/GaN", "700", "G_inel"), s) - \
                total("Al/GaN", rr("Al/GaN", "400", "G_bte"), rr("Al/GaN", "400", "G_inel"), s)
            mrise = float(meas["Al/GaN"]["700"]) - float(meas["Al/GaN"]["400"])
            g["I2"] = ("PASS" if abs(pred / mrise - 1.0) <= 0.35 else "FAIL") + \
                f"(pred rise {pred/1e6:.0f} vs meas {mrise/1e6:.0f})"
        else:
            g["I2"] = "FAIL(s* infinite)"
        nofit = [p for p in sstar if detail[p]["no_fit"]]
        g["reading"] = ("no-fit(" + ",".join(nofit) + ")") if nofit else classify(sstar)
        g["reading_pp"] = classify({p: detail[p]["s_pp"] for p in detail})
    else:
        g["I1"] = "UNDETERMINED(no measured rise)"
        g["I2"] = "UNDETERMINED(no measured rise)"
    try:
        v = {p: rr(p, "300", "G_total") for p in ("Al/Si", "Al/GaN")}
        g["I3"] = ("PASS" if all(abs(v[p] / MEAS_300[p] - 1.0) <= 0.20 for p in v) else "FAIL") + \
            "(" + ",".join(f"{p}:{v[p]/1e6:.0f}" for p in v) + ")"
    except KeyError as e:
        g["I3"] = f"UNDETERMINED(missing {e})"
    try:
        v = rr("AlN/GaN", "300", "G_total")
        g["I4"] = ("PASS" if ALNGAN_RANGE[0] <= v <= ALNGAN_RANGE[1] else "FAIL") + f"({v/1e6:.0f})"
    except KeyError as e:
        g["I4"] = f"UNDETERMINED(missing {e})"
    g["primary"] = prim
    summary = {"version": VERSION, "inel": ti.VERSION, "grades": g, "s_star": sstar,
               "s_star_detail": detail, "G_ep_Al": G_EP_AL, "meas": meas, "results": res}
    with open(os.path.join(args.out, "tbi-summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1, default=float)
    for k, v in g.items():
        log(f"{k}: {v}")
    log("TBI: " + " | ".join(f"{k} {str(v).split('(')[0]}" for k, v in g.items()))
    with open(os.path.join(args.out, "run.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
