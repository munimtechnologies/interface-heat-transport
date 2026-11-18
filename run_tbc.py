# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""run_tbc: stage-1 baselines for Grand Challenge 05 (thermal boundary conductance).

For each universal potential: relax and phonon the benchmark materials, build spectral fluxes,
then compute elastic G(T) (full-band DMM, gray-Simons-corrected DMM, radiation limits, Debye AMM
and DMM, electron-phonon series) for the benchmark interfaces, and grade the expectations B1-B3,
E1-E4 written in PLAN.md before this job ran. Last line of run.txt starts with 'TBC:'.
"""
import argparse
import gc
import json
import os
import time
import traceback

import numpy as np

import tbc_transport as tt

VERSION = "mv2 run-tbc 0.1.0"

# Isotropic (polycrystal) density and sound speeds for the Debye models (kg/m^3, m/s); rounded
# handbook values. TiC from K = 242 GPa, G = 190 GPa. GaN, AlN: approximate isotropic averages.
DEBYE = {
    "Al": (2700.0, 6420.0, 3040.0), "Si": (2329.0, 8433.0, 5843.0), "GaN": (6150.0, 7960.0, 4130.0),
    "AlN": (3255.0, 11000.0, 6300.0), "Cu": (8960.0, 4760.0, 2325.0), "C": (3515.0, 17500.0, 12800.0),
    "TiC": (4930.0, 10000.0, 6200.0), "Au": (19300.0, 3240.0, 1200.0), "Ti": (4506.0, 6070.0, 3125.0),
    "Ni": (8908.0, 5810.0, 3080.0),
}
# Electron-phonon inputs of metals (g in W m^-3 K^-1, phonon conductivity Lambda_p in W m^-1 K^-1).
# Approximate literature values (Al: g from Waldecker et al. 2016 as used by Li et al. 2022, with
# 4.9e17 as the high variant; Lambda_p from first-principles phonon calculations). They are bulk
# metal properties, not fitted to any interface.
EPH = {"Al": (2.45e17, 6.0), "Cu": (1.0e17, 2.6), "Au": (2.2e16, 2.0)}
EPH_HIGH = {"Al": (4.9e17, 6.0), "Cu": (2.0e17, 2.6), "Au": (4.4e16, 2.0)}

NORMALS = {"Al": [1, 1, 1], "Si": [1, 1, 1], "GaN": [0, 0, 1], "AlN": [0, 0, 1]}
TEMPS = [80.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 500.0, 600.0, 700.0]
ORIENTED_PAIRS = [("Al", "Si"), ("Al", "GaN"), ("Al", "AlN"), ("AlN", "GaN")]
ISO_PAIRS = [("Al", "Si"), ("Al", "GaN"), ("Cu", "C"), ("Cu", "TiC"), ("TiC", "C"), ("Au", "Si"),
             ("Ti", "Si"), ("Ni", "Si"), ("Al", "AlN")]
MATS = ["Al", "Si", "GaN", "AlN", "Cu", "C", "TiC", "Au", "Ti", "Ni"]

MEAS = {"AlSi300": 379e6, "AlGaN300": 423e6, "TiC_ratio": 1.48}
KAPPA_TIC_FILM = 20.0   # W/m/K, 10 nm TiC film (order of magnitude; enters only as L/kappa)
L_TIC = 10e-9
G_METAL_METAL = 4e9     # W/m^2/K, electronic Cu/TiC contact (order of Gundrum et al. 2005 Al/Cu)


def debye_dict(name, n_atoms):
    rho, vl, vt = DEBYE[name]
    return {"rho": rho, "vL": vl, "vT": vt, "n": n_atoms}


def pair_G(Ma, Mb, edges, T):
    g = tt.g_dmm(Ma, Mb, edges, T)
    ra = tt.g_radiation(Ma, edges, T)
    rb = tt.g_radiation(Mb, edges, T)
    return {"dmm": g, "dmm_simons": tt.simons_correct(g, ra, rb), "rad1": ra, "rad2": rb}


def run_model(model, edges, out, budget_s, t_start, log):
    import tbc_phonons as tp
    calc, disp = tp.load_calculator(model)
    spectra, info = {}, {}
    for name in MATS:
        if time.time() - t_start > budget_s:
            log(f"{model} {name}: SKIPPED (time budget)")
            continue
        try:
            r = tp.material_phonons(name, calc, disp)
        except Exception as e:  # keep going; record the failure
            log(f"{model} {name}: FAILED {type(e).__name__}: {e}")
            info[name] = {"failed": f"{type(e).__name__}: {e}"}
            continue
        sp = {"iso": tt.flux_spectrum(r["freq"], r["vel"], r["weights"], r["vcell"], edges)}
        if name in NORMALS:
            sp["or"] = tt.flux_spectrum(r["freq"], r["vel"], r["weights"], r["vcell"], edges,
                                        normal=NORMALS[name])
        spectra[name] = sp
        ref = tp.REF_MAX_THZ.get(name)
        info[name] = {k: r[k] for k in ("relax_converged", "cell_lengths", "n_atoms_per_m3",
                                        "max_thz", "min_thz", "imag_fraction", "seconds",
                                        "n_disp", "n_prim")}
        info[name]["ref_max_thz"] = ref
        info[name]["max_ratio"] = (r["max_thz"] / ref) if ref else None
        log(f"{model} {name}: a,b,c={['%.4f' % x for x in r['cell_lengths']]} max={r['max_thz']:.2f} "
            f"THz (ref {ref}) min={r['min_thz']:.3f} imag={r['imag_fraction']:.4f} "
            f"conv={r['relax_converged']} {r['seconds']} s")
        del r
        gc.collect()
    np.savez_compressed(os.path.join(out, f"spectra-{model}.npz"), edges=edges,
                        **{f"{k}_{v}": a for k, d in spectra.items() for v, a in d.items()})
    del calc
    gc.collect()
    return spectra, info


def tables(spectra, info, edges):
    res = {"oriented": {}, "iso": {}, "debye": {}}
    for kind, pairs in (("oriented", ORIENTED_PAIRS), ("iso", ISO_PAIRS)):
        for a, b in pairs:
            key = "or" if kind == "oriented" else "iso"
            if a not in spectra or b not in spectra or key not in spectra[a] or key not in spectra[b]:
                continue
            res[kind][f"{a}/{b}"] = {str(int(T)): pair_G(spectra[a][key], spectra[b][key], edges, T)
                                    for T in TEMPS}
    for a, b in ISO_PAIRS:
        if a in info and b in info and "n_atoms_per_m3" in info[a] and "n_atoms_per_m3" in info[b]:
            da = debye_dict(a, info[a]["n_atoms_per_m3"])
            db = debye_dict(b, info[b]["n_atoms_per_m3"])
            res["debye"][f"{a}/{b}"] = {str(int(T)): {"amm": tt.g_amm_debye(da, db, edges, T),
                                                     "dmm": tt.g_dmm_debye(da, db, edges, T)}
                                        for T in (100.0, 300.0, 700.0)}
    return res


def grade(res_by_model, info_by_model):
    models = list(res_by_model)
    g = {}
    lines = []
    # B2 per model
    for m in models:
        ok, bad = True, []
        for name, d in info_by_model[m].items():
            if "failed" in d:
                ok = False
                bad.append(f"{name}:failed")
                continue
            if d["imag_fraction"] >= 0.01:
                ok = False
                bad.append(f"{name}:imag{d['imag_fraction']:.3f}")
            if d["max_ratio"] is not None and abs(d["max_ratio"] - 1.0) >= 0.15:
                ok = False
                bad.append(f"{name}:max{d['max_ratio']:.2f}")
        missing = [n for n in MATS if n not in info_by_model[m]]
        if missing:
            ok = False
            bad.append("missing:" + ",".join(missing))
        g[f"B2_{m}"] = "PASS" if ok else "FAIL(" + ";".join(bad) + ")"
    # B3
    try:
        v = [res_by_model[m]["iso"]["Al/Si"]["300"]["dmm"] for m in models]
        g["B3"] = ("PASS" if abs(v[0] / v[1] - 1.0) < 0.15 else "FAIL") + f"({v[0]/1e6:.0f} vs {v[1]/1e6:.0f})"
    except (KeyError, IndexError):
        g["B3"] = "UNDETERMINED"
    gal, gal_hi = tt.g_ep(*EPH["Al"]), tt.g_ep(*EPH_HIGH["Al"])
    gcu = tt.g_ep(*EPH["Cu"])
    for m in models:
        r = res_by_model[m]
        try:
            o = r["oriented"]
            e1 = all(o[p]["700"]["dmm"] / o[p]["400"]["dmm"] < 1.10 for p in ("Al/Si", "Al/GaN"))
            g[f"E1_{m}"] = ("PASS" if e1 else "FAIL") + "(" + ",".join(
                f"{p}:{o[p]['700']['dmm'] / o[p]['400']['dmm']:.3f}" for p in ("Al/Si", "Al/GaN")) + ")"
            s = o["Al/Si"]["300"]
            lo = tt.series(s["dmm"], gal)
            hi = tt.series(s["dmm_simons"], gal)
            g[f"E2_{m}"] = ("PASS" if lo < MEAS["AlSi300"] < hi else "FAIL") + \
                f"(conv+ep {lo/1e6:.0f} < 379 < simons+ep {hi/1e6:.0f}; ep {gal/1e6:.0f}, hi-g ep {gal_hi/1e6:.0f})"
            rat = o["Al/GaN"]["300"]["dmm"] / o["Al/Si"]["300"]["dmm"]
            rat_s = o["Al/GaN"]["300"]["dmm_simons"] / o["Al/Si"]["300"]["dmm_simons"]
            g[f"E3_{m}"] = ("PASS" if rat < 1.12 else "FAIL") + f"(GaN/Si conv {rat:.3f}, simons {rat_s:.3f}; meas 1.116)"
        except KeyError as e:
            for k in ("E1", "E2", "E3"):
                g.setdefault(f"{k}_{m}", f"UNDETERMINED(missing {e})")
        try:
            i = r["iso"]
            outs = []
            for var in ("dmm", "dmm_simons"):
                with_ = tt.series(G_METAL_METAL, i["TiC/C"]["300"][var], KAPPA_TIC_FILM / L_TIC)
                without = tt.series(gcu, i["Cu/C"]["300"][var])
                outs.append(with_ / without)
            g[f"E4_{m}"] = ("PASS" if all(x > MEAS["TiC_ratio"] for x in outs) else "FAIL") + \
                f"(ratio conv {outs[0]:.2f}, simons {outs[1]:.2f}; meas 1.48)"
        except KeyError as e:
            g[f"E4_{m}"] = f"UNDETERMINED(missing {e})"
    return g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--models", default="mace-mpa-0,orb-v3")
    ap.add_argument("--budget-min", type=float, default=150.0)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    logf = open(os.path.join(args.out, "progress.txt"), "w")

    def log(s):
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    log(f"{VERSION} | {tt.VERSION}")
    edges = tt.make_edges(48.0, 0.02)
    res_by_model, info_by_model = {}, {}
    for m in [x.strip() for x in args.models.split(",") if x.strip()]:
        try:
            spectra, info = run_model(m, edges, args.out, args.budget_min * 60.0, t0, log)
        except Exception as e:
            log(f"MODEL {m} FAILED: {type(e).__name__}: {e}")
            log(traceback.format_exc())
            continue
        info_by_model[m] = info
        res_by_model[m] = tables(spectra, info, edges)
    g = grade(res_by_model, info_by_model)
    summary = {"version": VERSION, "transport": tt.VERSION, "grades": g, "info": info_by_model,
               "results": res_by_model, "eph": {"nominal": EPH, "high": EPH_HIGH},
               "seconds": round(time.time() - t0, 1)}
    with open(os.path.join(args.out, "tbc-summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    for m, r in res_by_model.items():
        for kind in ("oriented", "iso"):
            for p, d in r[kind].items():
                row = " ".join(f"{T}:{d[T]['dmm']/1e6:.0f}/{d[T]['dmm_simons']/1e6:.0f}" for T in
                               ("100", "300", "400", "700"))
                log(f"G {m} {kind:8s} {p:8s} conv/simons MW/m2K {row} rad {d['300']['rad1']/1e6:.0f},{d['300']['rad2']/1e6:.0f}")
        for p, d in r["debye"].items():
            log(f"G {m} debye    {p:8s} 300K amm {d['300']['amm']/1e6:.0f} dmm {d['300']['dmm']/1e6:.0f}")
    for k, v in g.items():
        log(f"{k}: {v}")
    short = " | ".join(f"{k} {v.split('(')[0]}" for k, v in g.items())
    log("TBC: " + short)


if __name__ == "__main__":
    main()
