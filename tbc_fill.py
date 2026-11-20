# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""tbc_fill: are the TiN conductances 5-6x below DMM an artefact of unsampled low-frequency bins?
(amendment 0.5.8, section 4)

Job tbc-diag-v1 printed, for TiN at 300 K, a two-stream bulk conductivity sum_b 2 K_b l_b of
132.3 W/m/K, while the binned phono3py conductivity it was built from is 44.6 W/m/K. The model is
meant to reproduce the bulk value exactly (tbc_bte docstring). The difference comes from the 11
TiN bins that have a radiation conductance K_b > 0 on the dense harmonic mesh but no phono3py mode
on the 16^3 lifetime mesh: tbc_bte.lengths_from_kappa gives them log-interpolated lengths (constant
extrapolation below the first sampled bin), and those invented lengths then carry at least
(132.3 - 44.6) / 132.3 = 66 % of the model's TiN conductivity.

This module reloads the tbc-tin-v1 npz files and re-solves every TiN pair with:

  base    the job's lengths (closure against tbc-diag-v1)                              -> F0
  merge   all bins from 0 up to the first bin sampled on BOTH sides' lifetime meshes merged into
          one band per side (K and kappa summed), so every low-frequency length is a real
          kappa / (2 K) and the bulk conductivity is conserved                          -> F2-F4
  nofill  filled bins kept but given a vanishing length (1e-12 m): they carry no bulk heat and
          stay at the local bath temperature (information)
  uni_m   merged bins with one common length per side (information, as tbc-diag's 'uni')

It also prints a per-bin table of the lowest 30 bins of each material at 300 K. tbc_bte, tbc_anh
and tbc_transport are imported unchanged.
"""
import math

VERSION = "mv2 tbc-fill 0.1.0"

PAIRS = (("TiN", "C"), ("TiN", "Si"), ("TiN", "Ge"))
T_HIGH = {"TiN/C": 900, "TiN/Si": 700, "TiN/Ge": 700}
CELL_A = {"TiN": 4.2538, "C": 3.5721, "Si": 5.4608, "Ge": 5.7581}
NODES = 400
L_TINY = 1e-12  # m

# tbc-diag-v1 (run.txt): base G and own DMM, MW m^-2 K^-1, one decimal
DIAG = {"TiN/C": {300: (43.3, 264.3), 900: (52.8, 361.8)},
        "TiN/Si": {300: (52.3, 315.3), 700: (61.5, 340.2)},
        "TiN/Ge": {300: (51.5, 272.4), 700: (59.3, 288.1)}}
# tbc-tin-v1 orb-v3 printed three-phonon inelastic channel G_inel (information only)
G_INEL = {"TiN/C": {300: 5.7, 900: 44.1}, "TiN/Si": {300: 56.4, 700: 172.4},
          "TiN/Ge": {300: 61.8, 700: 188.2}}
# Khan et al. 2024 measured values used in amendment 0.5.7 (information only)
MEASURED = {"TiN/C": {300: 449.0, 900: 910.0}, "TiN/Si": {300: 507.0, 700: 916.0},
            "TiN/Ge": {300: 190.0, 700: 281.0}}

# pre-registered thresholds (amendment 0.5.8, section 4)
F0_ABS, F0_REL = 0.5, 0.01
F1_SHARE = 0.60
F1_EDGE_THZ = 3.0
F2_FACTOR = 2.0
F3_FACTOR = 1.6
F4_REL = 0.05


def fcc_prim_volume_m3(a_ang):
    return (a_ang * 1e-10) ** 3 / 4.0


def first_have(K, kap):
    """First bin index with K > 0 and kappa > 0 (None if there is none)."""
    for i, (k, c) in enumerate(zip(K, kap)):
        if k > 0.0 and c > 0.0:
            return i
    return None


def merge_low(Ka, kapa, Kb, kapb):
    """Merge bins 0..j on both sides into bin j, j = max of the two sides' first sampled bins.
    K and kappa are summed, so sum K and sum kappa are conserved on each side.
    Returns (Ka', kapa', Kb', kapb', j) as lists of floats."""
    ja, jb = first_have(Ka, kapa), first_have(Kb, kapb)
    if ja is None or jb is None:
        raise ValueError("a side has no sampled bin")
    j = max(ja, jb)
    out = []
    for K, kap in ((Ka, kapa), (Kb, kapb)):
        K = [float(x) for x in K]
        kap = [float(x) for x in kap]
        K2 = [0.0] * j + [sum(K[:j + 1])] + K[j + 1:]
        kap2 = [0.0] * j + [sum(kap[:j + 1])] + kap[j + 1:]
        out.append((K2, kap2))
    return out[0][0], out[0][1], out[1][0], out[1][1], j


def filled_flags(K, kap):
    return [bool(k > 0.0 and not c > 0.0) for k, c in zip(K, kap)]


def two_stream_kappa(K, l):
    return sum(2.0 * k * x for k, x in zip(K, l) if k > 0.0)


def filled_share(K, l, kap):
    """Share of the model conductivity sum 2 K l carried by filled bins."""
    den = two_stream_kappa(K, l)
    num = sum(2.0 * k * x for k, x, f in zip(K, l, filled_flags(K, kap)) if f)
    return num / den if den > 0.0 else float("nan")


def uniform_lengths(K, kappa_total):
    tot = sum(k for k in K if k > 0.0)
    return [kappa_total / (2.0 * tot) if (k > 0.0 and tot > 0.0) else 0.0 for k in K]


def within(v, ref, abs_tol, rel_tol):
    return abs(v - ref) <= max(abs_tol, rel_tol * abs(ref))


def grade(rows, tin300):
    """rows[pair][T]: dict with base, dmm, merge, dmm_merge, kap1_merge, kap1_bins (MW, W/mK).
    tin300: dict with fill_share (TiN, 300 K, job lengths) and fill_lo_max_thz."""
    g = {}
    bad = []
    for p, by in DIAG.items():
        for T, (gb, _gd) in by.items():
            r = rows.get(p, {}).get(T)
            if r is None:
                bad.append(f"{p} {T} missing")
            elif not within(r["base"], gb, F0_ABS, F0_REL):
                bad.append(f"{p} {T} G {r['base']:.1f} vs {gb:.1f}")
    g["F0"] = "PASS" if not bad else "FAIL(" + "; ".join(bad) + ")"
    allr = [(p, T, r) for p, by in rows.items() for T, r in by.items()]
    if not allr or tin300 is None:
        for k in ("F1", "F2", "F3", "F4"):
            g[k] = "UNDETERMINED(no rows)"
        return g
    s, e = tin300["fill_share"], tin300["fill_lo_max_thz"]
    ok1 = s >= F1_SHARE and e < F1_EDGE_THZ
    g["F1"] = ("PASS" if ok1 else "FAIL") + f"(TiN 300 K filled share {s:.3f}, highest filled bin " \
        f"starts at {e:.2f} THz)"
    b2 = [f"{p} {T} {r['merge'] / r['base']:.2f}" for p, T, r in allr
          if not r["merge"] >= F2_FACTOR * r["base"]]
    g["F2"] = ("PASS" if not b2 else "FAIL(" + ", ".join(b2) + ")") + "(" + ", ".join(
        f"{p} {T} merge/base {r['merge'] / r['base']:.2f}" for p, T, r in allr) + ")"
    b3 = [f"{p} {T}" for p, T, r in allr if not r["merge"] * F3_FACTOR >= r["dmm_merge"]]
    g["F3"] = ("PASS" if not b3 else "FAIL(" + ", ".join(b3) + ")") + "(" + ", ".join(
        f"{p} {T} merge/dmm {r['merge'] / r['dmm_merge']:.2f}" for p, T, r in allr) + ")"
    b4 = [f"{p} {T} {r['kap1_merge']:.1f} vs {r['kap1_bins']:.1f}" for p, T, r in allr
          if not abs(r["kap1_merge"] / r["kap1_bins"] - 1.0) <= F4_REL]
    g["F4"] = "PASS" if not b4 else "FAIL(" + ", ".join(b4) + ")"
    return g


def main():
    import argparse
    import json
    import os
    import time

    import numpy as np

    import tbc_anh as ta
    import tbc_bte as tb

    ap = argparse.ArgumentParser()
    ap.add_argument("--npz-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--model", default="orb-v3")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    log(f"{VERSION} | {ta.VERSION} | {tb.VERSION}")
    mats = {}
    edges = None
    for m in ("TiN", "C", "Si", "Ge"):
        z = np.load(os.path.join(args.npz_dir, f"anh-{args.model}-{m}.npz"))
        r = {k: z[k] for k in z.files}
        temps = [float(x) for x in r["temps"]]
        it = temps.index(300.0)
        r1 = dict(r)
        r1["vcell"] = 1.0
        s = float(ta.kappa_bins(r1, it, 300.0, r["edges_c"], "iso", with_isotope=True).sum())
        k_ph3 = float(np.mean(r["kappa_ph3"][it][:3]))
        r["vcell"] = s / k_ph3
        log(f"MAT {m}: vcell/(a^3/4) {r['vcell'] / fcc_prim_volume_m3(CELL_A[m]):.4f}; kappa300 "
            f"ph3 {k_ph3:.1f} W/mK; lifetime-mesh lowest positive freq "
            f"{float(np.min(r['freq'][r['freq'] > 0.001])):.3f} THz")
        if edges is None:
            edges = np.asarray(r["edges_c"], dtype=float)
        elif not np.allclose(edges, r["edges_c"]):
            raise SystemExit(f"edges_c differ for {m}")
        mats[m] = (r, temps)

    def bins(m, T):
        r, temps = mats[m]
        i = temps.index(float(T))
        K = ta.K_bins(r["M_or"], r["edges_f"], r["edges_c"], float(T), classical=False)
        kap = ta.kappa_bins(r, i, float(T), r["edges_c"], "iso", with_isotope=True)
        return np.asarray(K, float), np.asarray(kap, float)

    # per-bin table, lowest 30 bins at 300 K
    tin300 = None
    for m in ("TiN", "C", "Si", "Ge"):
        K, kap = bins(m, 300)
        ell, nf = tb.lengths_from_kappa(kap, K)
        fl = filled_flags(list(K), list(kap))
        sh = filled_share(list(K), list(ell), list(kap))
        lo = [edges[i] / 1e12 for i, f in enumerate(fl) if f]
        log(f"BINS {m} 300 K: filled {nf} at lower edges (THz) {','.join(f'{x:.1f}' for x in lo)} | "
            f"model kappa {two_stream_kappa(list(K), list(ell)):.1f} vs bins {float(kap.sum()):.1f} "
            f"W/mK | filled share {sh:.3f} | sum K {float(K.sum()) / 1e6:.0f} MW/m2K")
        for i in range(30):
            log(f"BIN {m} {edges[i] / 1e12:5.1f} K {K[i] / 1e6:9.4f} kappa {kap[i]:9.4f} "
                f"l_nm {ell[i] * 1e9:10.2f} {'FILLED' if fl[i] else ''}")
        if m == "TiN":
            tin300 = {"fill_share": sh, "fill_lo_max_thz": max(lo) if lo else 0.0,
                      "n_filled": nf, "filled_lo_thz": lo}

    rows = {}
    for A, B in PAIRS:
        pair = f"{A}/{B}"
        rows[pair] = {}
        for T in (300, T_HIGH[pair]):
            Ka, kapa = bins(A, T)
            Kb, kapb = bins(B, T)
            la, fa = tb.lengths_from_kappa(kapa, Ka)
            lb, fb = tb.lengths_from_kappa(kapb, Kb)
            base = tb.solve_two_stream(Ka, Kb, la, lb, n_nodes=NODES)
            Ka2, kapa2, Kb2, kapb2, j = merge_low(list(Ka), list(kapa), list(Kb), list(kapb))
            Ka2, kapa2, Kb2, kapb2 = (np.array(x) for x in (Ka2, kapa2, Kb2, kapb2))
            la2, fa2 = tb.lengths_from_kappa(kapa2, Ka2)
            lb2, fb2 = tb.lengths_from_kappa(kapb2, Kb2)
            mg = tb.solve_two_stream(Ka2, Kb2, la2, lb2, n_nodes=NODES)
            la3 = np.where(np.array(filled_flags(list(Ka), list(kapa))), L_TINY, la)
            lb3 = np.where(np.array(filled_flags(list(Kb), list(kapb))), L_TINY, lb)
            nof = tb.solve_two_stream(Ka, Kb, la3, lb3, n_nodes=NODES)
            ua = np.array(uniform_lengths(list(Ka2), mg["kappa1"]))
            ub = np.array(uniform_lengths(list(Kb2), mg["kappa2"]))
            um = tb.solve_two_stream(Ka2, Kb2, ua, ub, n_nodes=NODES)
            mg800 = tb.solve_two_stream(Ka2, Kb2, la2, lb2, n_nodes=2 * NODES)
            r = {"base": base["G"] / 1e6, "dmm": base["G_dmm"] / 1e6, "merge": mg["G"] / 1e6,
                 "dmm_merge": mg["G_dmm"] / 1e6, "simons_merge": mg["G_simons"] / 1e6,
                 "linbath_merge": mg["G_linear_bath"] / 1e6, "nofill": nof["G"] / 1e6,
                 "uni_merge": um["G"] / 1e6, "merge_n800": mg800["G"] / 1e6,
                 "kap1_base": base["kappa1"], "kap2_base": base["kappa2"],
                 "kap1_merge": mg["kappa1"], "kap2_merge": mg["kappa2"],
                 "kap1_bins": float(kapa.sum()), "kap2_bins": float(kapb.sum()),
                 "merge_upto_thz": float(edges[j + 1] / 1e12), "filled_base": [fa, fb],
                 "filled_merge": [fa2, fb2],
                 "fill_share_merge": [filled_share(list(Ka2), list(la2), list(kapa2)),
                                      filled_share(list(Kb2), list(lb2), list(kapb2))],
                 "residual": max(x["residual"] for x in (base, mg, nof, um, mg800)),
                 "inel": G_INEL[pair][T], "measured": MEASURED[pair][T]}
            r["total_merge"] = r["merge"] + r["inel"]
            rows[pair][T] = r
            log(f"FILL {pair} {T} K | merged 0-{r['merge_upto_thz']:.1f} THz | G base {r['base']:.1f} "
                f"dmm {r['dmm']:.1f} | merge {r['merge']:.1f} (n800 {r['merge_n800']:.1f}) dmm_m "
                f"{r['dmm_merge']:.1f} simons_m {r['simons_merge']:.1f} linbath_m {r['linbath_merge']:.1f}"
                f" uni_m {r['uni_merge']:.1f} | nofill {r['nofill']:.1f} | kappa1 base/merge/bins "
                f"{r['kap1_base']:.1f}/{r['kap1_merge']:.1f}/{r['kap1_bins']:.1f} kappa2 "
                f"{r['kap2_base']:.1f}/{r['kap2_merge']:.1f}/{r['kap2_bins']:.1f} | filled base "
                f"{fa}/{fb} merge {fa2}/{fb2} share_m {r['fill_share_merge'][0]:.3f}/"
                f"{r['fill_share_merge'][1]:.3f} | info merge+inel {r['total_merge']:.0f} vs measured "
                f"{r['measured']:.0f} | resid {r['residual']:.1e}")
    g = grade(rows, tin300)
    for k, v in g.items():
        log(f"{k}: {v}")
    with open(os.path.join(args.out, "fill-summary.json"), "w") as fh:
        json.dump({"version": VERSION, "grades": g, "tin300": tin300,
                   "rows": {p: {str(T): v for T, v in by.items()} for p, by in rows.items()},
                   "seconds": round(time.time() - t0, 1)}, fh, indent=1, default=float)
    log("FIL: " + " | ".join(f"{k} {str(v).split('(')[0]}" for k, v in g.items()) + " | TiN300 share "
        f"{tin300['fill_share']:.2f} | " + " ; ".join(
            f"{p} {T} base/merge/dmm_m/nofill {r['base']:.0f}/{r['merge']:.0f}/{r['dmm_merge']:.0f}/"
            f"{r['nofill']:.0f}" for p, by in rows.items() for T, r in by.items()))
    with open(os.path.join(args.out, "run.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
