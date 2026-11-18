# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""run_tbc_anh: stage 2a of the interface-heat-transport study (PLAN-amendment-0.2.md).

For each universal potential (orb-v3 first, then MACE-MPA-0): bulk three-phonon and isotope
lifetimes of Al, Si, GaN, AlN and diamond (phono3py, RTA), then the multiband two-stream Boltzmann
interface model (tbc_bte) for Al/Si, Al/GaN and AlN/GaN at 100-700 K, next to the conventional and
gray-corrected DMM. Grades A0, A1, S2, S3, S5, S6, written before this job ran. Last line of
run.txt starts with 'TBA:'.
"""
import argparse
import gc
import json
import os
import time
import traceback

import numpy as np

import tbc_bte as tb
import tbc_transport as tt

VERSION = "mv2 run-tbc-anh 0.1.0"

TEMPS = [100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0]
MATS = ["Al", "Si", "GaN", "AlN", "C"]
NORMALS = {"Al": [1, 1, 1], "Si": [1, 1, 1], "GaN": [0, 0, 1], "AlN": [0, 0, 1], "C": [1, 1, 1]}
KDIR = {"Al": "iso", "Si": "iso", "C": "iso", "GaN": "z", "AlN": "z"}  # conductivity along normal
PAIRS = [("Al", "Si"), ("Al", "GaN"), ("AlN", "GaN")]
EPH_AL = (2.45e17, 6.0)   # same Al electron-phonon inputs as stage 1 (run_tbc.EPH)
MEAS = {"Al/Si": 379e6, "Al/GaN": 423e6, "AlN/GaN": (320e6, 490e6)}


def per_material(model, calc, edges_f, edges_c, out, log, t_start, budget_s):
    import tbc_anh as ta
    data = {}
    for name in MATS:
        if time.time() - t_start > budget_s:
            log(f"{model} {name}: SKIPPED (time budget)")
            continue
        try:
            r = ta.anharmonic(name, calc, TEMPS, log=log)
            fm = ta.flux_mesh(r)
        except Exception as e:
            log(f"{model} {name}: FAILED {type(e).__name__}: {e}")
            log(traceback.format_exc())
            continue
        M_or = tt.flux_spectrum(fm["freq"], fm["vel"], fm["weights"], fm["vcell"], edges_f,
                                normal=NORMALS[name])
        M_iso = tt.flux_spectrum(fm["freq"], fm["vel"], fm["weights"], fm["vcell"], edges_f)
        d = {"info": {k: r[k] for k in ("relax_converged", "cell_lengths", "n_fc3", "seconds",
                                        "seconds_forces")}}
        d["info"]["min_thz_dense"] = float(fm["freq"].min())
        d["info"]["max_thz_dense"] = float(fm["freq"].max())
        d["info"]["imag_fraction_dense"] = float(np.mean(fm["freq"] < -0.05))
        d["kappa"] = {}
        bins = {}
        for it, T in enumerate(TEMPS):
            key = str(int(T))
            kd = {}
            for dirn in sorted({"iso", KDIR[name]}):
                k_iso = ta.kappa_bins(r, it, T, edges_c, dirn, with_isotope=True)
                k_3 = ta.kappa_bins(r, it, T, edges_c, dirn, with_isotope=False)
                k_cl = ta.kappa_bins(r, it, T, edges_c, dirn, with_isotope=False, classical=True)
                bins[(key, dirn)] = (k_iso, k_3, k_cl)
                kd[dirn] = {"sum_bins_iso": float(k_iso.sum()), "sum_bins_3ph": float(k_3.sum()),
                            "phono3py_iso": ta.kappa_ph3_dir(r, it, dirn)}
            d["kappa"][key] = kd
        d["M_or"], d["M_iso"], d["bins"] = M_or, M_iso, bins
        k300 = d["kappa"]["300"]
        log(f"{model} {name}: kappa300 " + ", ".join(
            f"{dn} bins {v['sum_bins_iso']:.1f} (3ph only {v['sum_bins_3ph']:.1f}) phono3py "
            f"{v['phono3py_iso']:.1f}" for dn, v in k300.items()) + f" W/mK; {r['seconds']} s")
        np.savez_compressed(os.path.join(out, f"anh-{model}-{name}.npz"), edges_c=edges_c,
                            freq=r["freq"], weights=r["weights"], gamma3=r["gamma3"],
                            gamma_iso=r["gamma_iso"], v2_iso=r["v2_iso"], v2_z=r["v2_z"],
                            kappa_ph3=r["kappa_ph3"], temps=np.array(TEMPS), M_or=M_or,
                            M_iso=M_iso, edges_f=edges_f)
        data[name] = d
        del r, fm
        gc.collect()
    return data


def pair_tables(data, edges_f, edges_c, nodes):
    res = {}
    for a, b in PAIRS:
        if a not in data or b not in data:
            continue
        for geo in ("oriented", "iso"):
            rows = {}
            for T in TEMPS:
                key = str(int(T))
                da, db = data[a], data[b]
                Ma = da["M_or"] if geo == "oriented" else da["M_iso"]
                Mb = db["M_or"] if geo == "oriented" else db["M_iso"]
                dira = KDIR[a] if geo == "oriented" else "iso"
                dirb = KDIR[b] if geo == "oriented" else "iso"
                row = {}
                for var, idx, classical in (("iso", 0, False), ("3ph", 1, False), ("classical", 2, True)):
                    Ka = K_b(Ma, edges_f, edges_c, T, classical)
                    Kb = K_b(Mb, edges_f, edges_c, T, classical)
                    la, fa = tb.lengths_from_kappa(da["bins"][(key, dira)][idx], Ka)
                    lb, fb = tb.lengths_from_kappa(db["bins"][(key, dirb)][idx], Kb)
                    s = tb.solve_two_stream(Ka, Kb, la, lb, n_nodes=nodes)
                    s["filled_bins"] = [fa, fb]
                    row[var] = s
                rows[key] = row
            res[f"{geo} {a}/{b}"] = rows
    return res


def K_b(M, edges_f, edges_c, T, classical):
    import tbc_anh as ta
    return ta.K_bins(M, edges_f, edges_c, T, classical=classical)


def grade(per_model, tables, models):
    g = {}
    prim = models[0] if models else None
    for m in models:
        bad = []
        for name, d in per_model[m].items():
            for dn, v in d["kappa"]["300"].items():
                if v["phono3py_iso"] > 0 and abs(v["sum_bins_iso"] / v["phono3py_iso"] - 1.0) >= 0.03:
                    bad.append(f"{name}-{dn}:{v['sum_bins_iso'] / v['phono3py_iso']:.3f}")
        g[f"A0_{m}"] = ("PASS" if not bad and per_model[m] else "FAIL") + "(" + ";".join(bad) + ")"
        parts, ok = [], True
        for name, (lo, hi) in tb_ranges().items():
            if name not in per_model[m]:
                ok = False
                parts.append(f"{name}:missing")
                continue
            k = per_model[m][name]["kappa"]["300"]["iso"]["phono3py_iso"]
            ok &= lo <= k <= hi
            parts.append(f"{name}:{k:.0f}")
        g[f"A1_{m}"] = ("PASS" if ok else "FAIL") + "(" + ",".join(parts) + ")"
    for m in models:
        t = tables.get(m, {})
        try:
            rat = {p: t[f"oriented {p}"]["700"]["iso"]["G"] / t[f"oriented {p}"]["400"]["iso"]["G"]
                   for p in ("Al/Si", "Al/GaN")}
            g[f"S2_{m}"] = ("PASS" if all(v < 1.10 for v in rat.values()) else "FAIL") + "(" + \
                ",".join(f"{p}:{v:.3f}" for p, v in rat.items()) + ")"
        except KeyError as e:
            g[f"S2_{m}"] = f"UNDETERMINED(missing {e})"
        try:
            s = t["oriented Al/Si"]["300"]["iso"]
            g[f"S3_{m}"] = ("PASS" if s["G_dmm"] < s["G"] < s["G_simons"] else "FAIL") + \
                f"(dmm {s['G_dmm']/1e6:.0f} < bte {s['G']/1e6:.0f} < simons {s['G_simons']/1e6:.0f})"
        except KeyError as e:
            g[f"S3_{m}"] = f"UNDETERMINED(missing {e})"
        try:
            rr = {p: t[f"oriented {p}"]["300"]["3ph"]["G"] / t[f"oriented {p}"]["300"]["iso"]["G"]
                  for p in ("Al/Si", "Al/GaN")}
            g[f"S5_{m}"] = ("PASS" if all(abs(v - 1) < 0.05 for v in rr.values()) else "FAIL") + "(" + \
                ",".join(f"{p}:{v:.3f}" for p, v in rr.items()) + ")"
        except KeyError as e:
            g[f"S5_{m}"] = f"UNDETERMINED(missing {e})"
    if len(models) >= 2:
        try:
            v = [tables[m]["oriented Al/Si"]["300"]["iso"]["G"] for m in models[:2]]
            g["S6"] = ("PASS" if abs(v[0] / v[1] - 1) < 0.15 else "FAIL") + f"({v[0]/1e6:.0f} vs {v[1]/1e6:.0f})"
        except KeyError as e:
            g["S6"] = f"UNDETERMINED(missing {e})"
    else:
        g["S6"] = "UNDETERMINED(one model)"
    g["primary"] = str(prim)
    return g


def tb_ranges():
    import tbc_anh as ta
    return ta.KAPPA_RANGE_300


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--models", default="orb-v3,mace-mpa-0")
    ap.add_argument("--budget-min", type=float, default=150.0)
    ap.add_argument("--nodes", type=int, default=400)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    logf = open(os.path.join(args.out, "progress.txt"), "w")

    def log(s):
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    import tbc_anh as ta
    import tbc_phonons as tp
    log(f"{VERSION} | {tb.VERSION} | {ta.VERSION} | {tt.VERSION} | {tp.VERSION}")
    edges_f = tt.make_edges(48.0, 0.02)
    edges_c = tt.make_edges(48.0, 0.2)
    per_model, tables, models = {}, {}, []
    for m in [x.strip() for x in args.models.split(",") if x.strip()]:
        if time.time() - t0 > args.budget_min * 60.0:
            log(f"MODEL {m}: SKIPPED (time budget)")
            continue
        try:
            calc, _ = tp.load_calculator(m)
            data = per_material(m, calc, edges_f, edges_c, args.out, log, t0, args.budget_min * 60.0)
            del calc
            gc.collect()
        except Exception as e:
            log(f"MODEL {m} FAILED: {type(e).__name__}: {e}")
            log(traceback.format_exc())
            continue
        if not data:
            continue
        models.append(m)
        per_model[m] = data
        tables[m] = pair_tables(data, edges_f, edges_c, args.nodes)
        # grid check on one case
        try:
            d = tables[m]["oriented Al/Si"]["300"]["iso"]
            da, db = data["Al"], data["Si"]
            Ka = K_b(da["M_or"], edges_f, edges_c, 300.0, False)
            Kb = K_b(db["M_or"], edges_f, edges_c, 300.0, False)
            la, _ = tb.lengths_from_kappa(da["bins"][("300", "iso")][0], Ka)
            lb, _ = tb.lengths_from_kappa(db["bins"][("300", "iso")][0], Kb)
            gh = tb.solve_two_stream(Ka, Kb, la, lb, n_nodes=args.nodes // 2)["G"]
            log(f"GRID {m} Al/Si 300 K: G(n={args.nodes}) {d['G']/1e6:.2f}, G(n={args.nodes // 2}) {gh/1e6:.2f} MW/m2K")
        except Exception as e:
            log(f"GRID check skipped: {type(e).__name__}: {e}")
    g = grade(per_model, tables, models)
    gep = tt.g_ep(*EPH_AL)
    summary = {"version": VERSION, "bte": tb.VERSION, "grades": g, "models": models,
               "info": {m: {n: {"info": d["info"], "kappa": d["kappa"]} for n, d in per_model[m].items()}
                        for m in models},
               "tables": tables, "g_ep_Al": gep, "seconds": round(time.time() - t0, 1)}
    with open(os.path.join(args.out, "tba-summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1, default=float)
    for m in models:
        for p, rows in tables[m].items():
            line = " ".join(
                f"{T}:{rows[T]['iso']['G_dmm']/1e6:.0f}/{rows[T]['iso']['G']/1e6:.0f}/{rows[T]['3ph']['G']/1e6:.0f}/"
                f"{rows[T]['classical']['G']/1e6:.0f}/{rows[T]['iso']['G_simons']/1e6:.0f}"
                for T in ("100", "300", "400", "500", "700") if T in rows)
            log(f"G {m} {p:18s} T:dmm/bte/bte3ph/bteclassical/simons MW/m2K {line}")
            if p.startswith("oriented Al/") and "300" in rows:
                s = tt.series(rows["300"]["iso"]["G"], gep)
                log(f"  {m} {p} 300 K with Al e-ph in series: {s/1e6:.0f} MW/m2K (measured {MEAS[p.split()[1]]/1e6:.0f})")
    for k, v in g.items():
        log(f"{k}: {v}")
    log("TBA: " + " | ".join(f"{k} {v.split('(')[0]}" for k, v in g.items()))


if __name__ == "__main__":
    main()
