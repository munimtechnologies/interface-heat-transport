# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""run_tbc_zno: stage 3a of the interface-heat-transport study (PLAN-amendment-0.5.md).

Two non-metal interfaces with no electrons in the way: ZnO/GaN (Gaskins et al., 77-500 K) and
AlN/GaN (UCLA bonded wafers to 300 K; Walwil et al. to 400 K). For each open universal potential
(orb-v3 first, then MACE-MPA-0, as in stage 2a) it computes bulk three-phonon lifetimes of ZnO, GaN
and AlN (tbc_anh, phono3py RTA), the elastic two-stream Boltzmann conductance G_BTE (tbc_bte), and
the no-fit interfacial three-phonon channel G_inel (tbc_inel, stage 2b) at 100-700 K, for the
[0001] interface. G_total = G_BTE + G_inel (no metal, so no electron-phonon term). Grades Z0-Z5
(tbc_nonmetal.grade_z), written before this job ran. Last line of run.txt starts with 'TBZ:'.

Nothing is copied: the stage 2a and 2b code is imported and ZnO is registered in its tables.
"""
import argparse
import gc
import json
import os
import time
import traceback

import tbc_anh as ta
import tbc_inel as ti
import tbc_nonmetal as tn
import tbc_phonons as tp
import tbc_transport as tt
import run_tbc_anh as ra
import run_tbc_inel as ri

VERSION = "mv2 run-tbc-zno 0.1.0"

# Room-temperature wurtzite ZnO (a, c in Angstrom, internal parameter u); relaxed by each potential.
ZNO = dict(kind="wurtzite", el="ZnO", a=3.2496, c=5.2042, u=0.3819, sc=(4, 4, 3), mesh=24)
MATS = ["ZnO", "GaN", "AlN"]
PAIRS = [("ZnO", "GaN"), ("AlN", "GaN")]


def register_zno():
    tp.MATERIALS["ZnO"] = dict(ZNO)
    tp.REF_MAX_THZ["ZnO"] = tn.ZNO_REF_MAX_THZ
    ta.ANH["ZnO"] = dict(sc3=(3, 3, 2), cut=4.0, mesh=(12, 12, 8), lattice="hex")
    ta.KAPPA_RANGE_300["ZnO"] = tn.ZNO_KAPPA_RANGE
    ra.NORMALS["ZnO"] = [0, 0, 1]
    ra.KDIR["ZnO"] = "z"
    ra.MATS = list(MATS)
    ra.PAIRS = list(PAIRS)
    ri.GEOM["ZnO"] = "wurtzite0001"
    ri.NCONV["ZnO"] = 4


def inel_rows(out, model, data, tables, log):
    res = {}
    for A, B in PAIRS:
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
            cl = ti.g_inel(a_T, b_T, T, kernel="classical", classical_c=True, q_cache=cache)
            el = tables[key][k]["iso"]
            gb = float(el["G"])
            rows[k] = {"G_bte": gb / 1e6, "G_inel": pr["G"] / 1e6,
                       "G_inel_classical_c": cl["G"] / 1e6, "G_total": (gb + pr["G"]) / 1e6,
                       "G_dmm": float(el["G_dmm"]) / 1e6, "G_simons": float(el["G_simons"]) / 1e6,
                       "G_bte_3ph": float(tables[key][k]["3ph"]["G"]) / 1e6,
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
                     "kappa300": k["iso"]["phono3py_iso"],
                     "kappa300_z": k.get("z", {}).get("phono3py_iso")}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--models", default="orb-v3,mace-mpa-0")
    ap.add_argument("--budget-min", type=float, default=150.0)
    ap.add_argument("--nodes", type=int, default=400)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    register_zno()
    log(f"{VERSION} | {tn.VERSION} | {ra.VERSION} | {ri.VERSION} | {ti.VERSION} | {ta.VERSION} | "
        f"{tt.VERSION} | {tp.VERSION}")
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
            inel = inel_rows(args.out, m, data, tables, log)
            per_model[m] = {"materials": material_info(data),
                            "pairs": {p: v["rows"] for p, v in inel.items()},
                            "t_eff_m": {p: v["t_eff"] for p, v in inel.items()}}
            models.append(m)
        except Exception as e:
            log(f"MODEL {m} FAILED: {type(e).__name__}: {e}")
            log(traceback.format_exc())
    prim = models[0] if models else None
    grades = tn.grade_z(per_model[prim]) if prim else {
        k: "UNDETERMINED(no model)" for k in ("Z0", "Z1", "Z2", "Z3", "Z4", "Z5")}
    grades["primary"] = str(prim)
    cross = {}
    if len(models) >= 2:
        cross = {m: {k: v.split("(")[0] for k, v in tn.grade_z(per_model[m]).items()} for m in models}
    summary = {"version": VERSION, "grades": grades, "per_model_grades": cross, "models": models,
               "results": per_model, "seconds": round(time.time() - t0, 1),
               "thresholds": {"Z1_MAX": tn.Z1_MAX, "Z2_MAX": tn.Z2_MAX, "Z3_RANGE": tn.Z3_RANGE,
                              "Z4_RANGE": tn.Z4_RANGE}}
    with open(os.path.join(args.out, "tbz-summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1, default=float)
    for m in models:
        for name, v in per_model[m]["materials"].items():
            log(f"{m} {name}: kappa300 iso {v['kappa300']:.1f} z {v['kappa300_z']} W/mK, max "
                f"{v['max_thz']:.2f} THz, imag {v['imag_fraction']:.4f}, relaxed {v['relax_converged']}")
    for k, v in grades.items():
        log(f"{k}: {v}")
    if cross:
        log("by model: " + "; ".join(f"{m} " + ",".join(f"{k} {v}" for k, v in g.items()) for m, g in cross.items()))
    log("TBZ: " + " | ".join(f"{k} {str(v).split('(')[0]}" for k, v in grades.items()))
    with open(os.path.join(args.out, "run.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
