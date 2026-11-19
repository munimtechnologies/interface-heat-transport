# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Made with AI agents working in MV2's name.
"""Reproduce PLAN-amendment-0.3.3 (direct electron-phonon channel vs Li et al. 2022).

python run_tbc_direct.py --meas data/meas-li2022.json --out results/direct
Writes direct-li2022.json; the last printed line starts with 'TBD:'.
"""
import argparse
import json
import os

import tbc_direct as td
from tbc_gep import gep_power_law, series

MW = 1e6
GEP400, N = 1.26e9, 0.5          # amendment 0.3.2: G_ep(T) = 1.26 GW (400/T)^0.5
FLOOR_AL_SI = 210e6              # ab initio harmonic phonon side, Al/Si, 300 K (arXiv:2510.27499)
AL = {"gamma_vol": 130.0, "lam_w2_meV2": 185.9, "u_L": 6240.0, "u_T": 3040.0, "n_e": 1.81e29,
      "g_DFT": 2.46e17, "k_e": 228.0, "k_p": 8.95,
      "sources": "gamma, lambda<w^2>, u_L, u_T, n_e: Lombard et al. arXiv:1501.03176 Table 1; "
                 "g, k_p, k: De San Feliciano et al. arXiv:2510.27499 (k_e = 237 - 8.95)"}
TEMPS = (300, 400, 500, 600, 700)


def r6(x):
    return None if x is None else round(x / MW, 3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meas", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    with open(a.meas) as f:
        meas = json.load(f)
    os.makedirs(a.out, exist_ok=True)
    pairs = {"Al/Si": meas["Al/Si"], "Al/GaN": meas["Al/GaN"]}

    serg = {
        "from_g_DFT": td.sergeev_from_g(AL["gamma_vol"], AL["g_DFT"], AL["u_L"], AL["u_T"], AL["n_e"], 300.0),
        "from_lambda": td.sergeev_from_lambda(AL["gamma_vol"], AL["lam_w2_meV2"], AL["u_L"], AL["u_T"], AL["n_e"], 300.0),
        "from_lambda_printed_eq4": td.sergeev_from_lambda(AL["gamma_vol"], AL["lam_w2_meV2"], AL["u_L"], AL["u_T"],
                                                          AL["n_e"], 300.0, zeta=False),
        "Lombard_table1": 264e6,
    }
    slopes = {k: v / 300.0 for k, v in serg.items()}
    slopes["DSF_107_at_300K_linear"] = 107e6 / 300.0

    res = {"version": td.VERSION, "inputs": AL, "G_ep": {"G400": GEP400, "n": N},
           "sergeev_sigma_e_300K_MW": {k: r6(v) for k, v in serg.items()},
           "g_from_lambda": td.g_from_lambda(AL["gamma_vol"], AL["lam_w2_meV2"]),
           "circuit_check_MW": {
               "exact_sigma_e_107_sigma_p_210": r6(td.circuit_exact(AL["g_DFT"], AL["k_e"], AL["k_p"], 107e6, 210e6)),
               "limit_sigma_e_107_sigma_p_210": r6(td.circuit_limit((AL["g_DFT"] * AL["k_p"]) ** 0.5, 107e6, 210e6))},
           "bounds": {}, "scan": {}, "sergeev_with_floor_Al_Si": {}, "sample2_min_G600_MW": {}}

    for name, m in pairs.items():
        s, t = td.linear_bound(m, TEMPS)
        res["bounds"][name] = {"s_max_no_floor_MW_per_K": round(s / MW, 5), "binds_at": t,
                               "sigma_e_300_max_MW": round(s * 300 / MW, 2)}
    s, t = td.floor_bound(pairs["Al/Si"], GEP400, N, FLOOR_AL_SI, TEMPS)
    res["bounds"]["Al/Si"].update({"s_max_floor210_MW_per_K": round(s / MW, 5), "floor_binds_at": t,
                                   "sigma_e_300_max_floor210_MW": round(s * 300 / MW, 2)})

    for label, s in [("0", 0.0), ("0.25", 0.25e6), ("0.30", 0.30e6), ("DSF", slopes["DSF_107_at_300K_linear"]),
                     ("0.50", 0.5e6)]:
        res["scan"][label] = {name: [r6(x) for x in td.infer_with_direct(m, GEP400, N, s, TEMPS)["G_pp"]]
                              for name, m in pairs.items()}

    for k, s in slopes.items():
        res["sergeev_with_floor_Al_Si"][k] = {
            str(T): r6(s * T + series(gep_power_law(GEP400, T, 400.0, N), FLOOR_AL_SI)) for T in TEMPS}
        mr = td.min_rise_prediction(330e6, 400.0, 600.0, s, GEP400, N)
        res["sample2_min_G600_MW"][k] = r6(mr["G_min"])

    with open(os.path.join(a.out, "direct-li2022.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))
    b = res["bounds"]
    print("TBD: sergeev Al 300K MW g/lam/eq4/table %.1f/%.1f/%.1f/264 | s_max AlSi %.4f AlGaN %.4f floor210 %.4f "
          "(sigma_e300 <= %.1f) | S2 G600 min DSF %.1f sergeev-g %.1f"
          % (serg["from_g_DFT"] / MW, serg["from_lambda"] / MW, serg["from_lambda_printed_eq4"] / MW,
             b["Al/Si"]["s_max_no_floor_MW_per_K"], b["Al/GaN"]["s_max_no_floor_MW_per_K"],
             b["Al/Si"]["s_max_floor210_MW_per_K"], b["Al/Si"]["sigma_e_300_max_floor210_MW"],
             res["sample2_min_G600_MW"]["DSF_107_at_300K_linear"], res["sample2_min_G600_MW"]["from_g_DFT"]))


if __name__ == "__main__":
    main()
