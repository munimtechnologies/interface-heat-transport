# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""run_tbc_tin: stage 3b of the interface-heat-transport study (PLAN-amendment-0.5.2.md).

TiN on diamond, Si and Ge, 100-1000 K, against Khan et al. (APL 125, 2024). For each open
universal potential (orb-v3 first, then MACE-MPA-0) it computes bulk three-phonon lifetimes of
TiN, C, Si and Ge (tbc_anh, phono3py RTA), the elastic two-stream Boltzmann (dressed) conductance
G_BTE (tbc_bte) and the no-fit interfacial three-phonon channel G_inel (tbc_inel) for (001)
interfaces. G_total = G_BTE + G_inel is graded (T0-T6, tbc_tin.grade_t, written before this job
ran). A series electron-phonon term inside TiN, sqrt(g Lambda_p) with g = 1e18 W m^-3 K^-1 and the
model's own TiN phonon conductivity, is reported as a sensitivity only. Last line of run.txt
starts with 'TBT:'.

Nothing is copied: the stage 2a and 2b code is imported and TiN and Ge are registered in its tables.
"""
import argparse
import gc
import json
import math
import os
import time
import traceback

import tbc_tin as tz

VERSION = "mv2 run-tbc-tin 0.1.0"


def inel_rows(ra, ri, ti, tt, out, model, data, tables, log):
    res = {}
    for A, B in tz.PAIR_TUPLES:
        pair = f"{A}/{B}"
        key = f"oriented {pair}"
        if A not in data or B not in data or key not in tables:
            log(f"{model} {pair}: SKIPPED (material missing)")
            continue
        sa = ri.load_side(out, model, A, data[A]["info"]["cell_lengths"])
        sb = ri.load_side(out, model, B, data[B]["info"]["cell_lengths"])
        cache = {}
        rows = {}
        for T in sa["temps"]:
            k = str(int(T))
            if k not in tables[key]:
                continue
            a_T, b_T = ri.side_at(sa, T), ri.side_at(sb, T)
            pr = ti.g_inel(a_T, b_T, T, kernel="classical", q_cache=cache)
            el = tables[key][k]["iso"]
            gb = float(el["G"])
            g_tot = (gb + pr["G"]) / 1e6
            lam = data[A]["kappa"].get(k, {}).get("iso", {}).get("phono3py_iso")
            g_ep = tt.g_ep(tz.G_EP_TIN, float(lam)) / 1e6 if lam and lam > 0 else float("nan")
            rows[k] = {"G_bte": gb / 1e6, "G_inel": pr["G"] / 1e6, "G_total": g_tot,
                       "G_dmm": float(el["G_dmm"]) / 1e6, "G_simons": float(el["G_simons"]) / 1e6,
                       "G_bte_3ph": float(tables[key][k]["3ph"]["G"]) / 1e6,
                       "G_ep_TiN": g_ep,
                       "G_total_series_ep": (1.0 / (1.0 / g_tot + 1.0 / g_ep)
                                             if g_ep == g_ep else float("nan")),
                       "Qbar": pr["Qbar"], "parts": [pr["G_side1"] / 1e6, pr["G_side2"] / 1e6],
                       "zero_nodes": pr["zero_nodes"], "capped_nodes": pr["capped_nodes"]}
        res[pair] = {"t_eff": [sa["t"], sb["t"]], "rows": rows}
        log(f"G {model} {pair} T:bte/inel/total (dmm) MW/m2K " + " ".join(
            f"{k}:{r['G_bte']:.0f}/{r['G_inel']:.1f}/{r['G_total']:.0f}({r['G_dmm']:.0f})"
            for k, r in rows.items()))
    return res


def material_info(data):
    out = {}
    for name, d in data.items():
        k = d["kappa"]["300"]
        out[name] = {"relax_converged": d["info"]["relax_converged"],
                     "cell_lengths": d["info"]["cell_lengths"],
                     "imag_fraction": d["info"]["imag_fraction_dense"],
                     "max_thz": d["info"]["max_thz_dense"],
                     "kappa300": k["iso"]["phono3py_iso"]}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--models", default="orb-v3,mace-mpa-0")
    ap.add_argument("--budget-min", type=float, default=170.0)
    ap.add_argument("--nodes", type=int, default=400)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    tz.register_tin()
    import tbc_anh as ta
    import tbc_inel as ti
    import tbc_phonons as tp
    import tbc_transport as tt
    import run_tbc_anh as ra
    import run_tbc_inel as ri
    log(f"{VERSION} | {tz.VERSION} | {ra.VERSION} | {ri.VERSION} | {ti.VERSION} | {ta.VERSION} | "
        f"{tt.VERSION} | {tp.VERSION}")
    log("registered: MATS " + ",".join(ra.MATS) + " | TEMPS " + ",".join(str(int(t)) for t in ra.TEMPS)
        + " | GEOM " + ",".join(f"{m}:{ri.GEOM[m]}" for m in ra.MATS))
    edges_f = tt.make_edges(48.0, 0.02)
    edges_c = tt.make_edges(48.0, 0.2)
    per_model = {}
    models = []
    for m in [x.strip() for x in args.models.split(",") if x.strip()]:
        if time.time() - t0 > args.budget_min * 60.0:
            log(f"MODEL {m}: SKIPPED (time budget)")
            continue
        try:
            calc, _ = tp.load_calculator(m)
            data = ra.per_material(m, calc, edges_f, edges_c, args.out, log, t0, args.budget_min * 60.0)
            del calc
            gc.collect()
            if not data:
                continue
            tables = ra.pair_tables(data, edges_f, edges_c, args.nodes)
            inel = inel_rows(ra, ri, ti, tt, args.out, m, data, tables, log)
            per_model[m] = {"materials": material_info(data),
                            "pairs": {p: v["rows"] for p, v in inel.items()},
                            "t_eff_m": {p: v["t_eff"] for p, v in inel.items()}}
            models.append(m)
        except Exception as e:
            log(f"MODEL {m} FAILED: {type(e).__name__}: {e}")
            log(traceback.format_exc())
    prim = models[0] if models else None
    grades = tz.grade_t(per_model[prim]) if prim else {
        k: "UNDETERMINED(no model)" for k in ("T0", "T1", "T2", "T3", "T4", "T5", "T6")}
    grades["primary"] = str(prim)
    cross = {}
    if len(models) >= 2:
        cross = {m: {k: v.split("(")[0] for k, v in tz.grade_t(per_model[m]).items()} for m in models}
    summary = {"version": VERSION, "grades": grades, "per_model_grades": cross, "models": models,
               "results": per_model, "seconds": round(time.time() - t0, 1),
               "measured": {p: {str(k): v for k, v in tz.KHAN[p].items()} for p in tz.PAIRS},
               "thresholds": {"T_HIGH": tz.T_HIGH, "T2_TOL": tz.T2_TOL, "T3_TOL": tz.T3_TOL,
                              "T6_RANGE": tz.T6_RANGE, "T0_MAX_THZ": tz.T0_MAX_THZ}}
    with open(os.path.join(args.out, "tbt-summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1, default=float)
    for m in models:
        for name, v in per_model[m]["materials"].items():
            log(f"{m} {name}: kappa300 iso {v['kappa300']:.1f} W/mK, max {v['max_thz']:.2f} THz, "
                f"imag {v['imag_fraction']:.4f}, relaxed {v['relax_converged']}")
        for p, rows in per_model[m]["pairs"].items():
            meas = tz.KHAN.get(p, {})
            log(f"{m} {p} model vs measured: " + " ".join(
                f"{k}:{r['G_total']:.0f}/{meas[int(k)]:.0f}" for k, r in rows.items()
                if int(k) in meas) + " | series-ep 300: " + (
                f"{rows['300']['G_total_series_ep']:.0f}" if "300" in rows and
                not math.isnan(rows["300"]["G_total_series_ep"]) else "n/a"))
    for k, v in grades.items():
        log(f"{k}: {v}")
    if cross:
        log("by model: " + "; ".join(f"{m} " + ",".join(f"{k} {v}" for k, v in g.items()) for m, g in cross.items()))
    log("TBT: " + " | ".join(f"{k} {str(v).split('(')[0]}" for k, v in grades.items()))
    with open(os.path.join(args.out, "run.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
