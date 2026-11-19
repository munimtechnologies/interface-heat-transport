# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""tbc_diag: why does the dressed multiband conductance fall below its own DMM? (amendment 0.5.7)

Job tbc-tin-v1 (orb-v3) gave G_BTE = 43, 52 and 51 MW m^-2 K^-1 at 300 K for TiN on diamond, Si
and Ge (001), while the band-summed DMM of the same band conductances is 264, 315 and 272. In the
one-band (gray) limit tbc_bte reproduces the Simons / Zeng-Chen dressed conductance, which is
always ABOVE the DMM value. So the factor 5-6 deficit must come either from the spread of band
lengths (non-equilibrium between bands near the interface, which the multiband model allows) or
from a defect in the solver or its inputs. This module reloads the per-material outputs that
tbc-tin-v1 saved (anh-orb-v3-<mat>.npz, published in results/tin/orb-v3/), rebuilds the exact band
conductances and two-stream lengths the job used, and re-solves with:

  base      the job's lengths (closure check D0 against the printed numbers)
  uni       every band on BOTH sides given the same length, same bulk conductivity (gray lengths,
            multiband DMM transmission kept)                                       -> D3
  uni1      uniform lengths on side 1 (TiN) only                                   -> D4
  uni2      uniform lengths on side 2 (C, Si, Ge) only                             -> D4
  x0.5, x2  every length scaled by 0.5 and 2 (the model is scale invariant)        -> D1
  n800      800 nodes instead of 400                                               -> D2
  s1x0.01   TiN lengths x 0.01 only (crude stand-in for electron-pinned TiN phonons; information)

It also prints the share of each side's conductivity carried by bands that the other side can
receive (shared bands), and the kappa-weighted mean length of all bands and of shared bands.

The tbc_bte, tbc_anh and tbc_transport modules are imported unchanged from the repository; vcell,
which the npz files do not store, is recovered from the stored phono3py conductivity (the job
printed bins = phono3py to 0.1 W/mK) and checked against a^3/4 of the relaxed cell.
"""
import math

VERSION = "mv2 tbc-diag 0.1.0"

PAIRS = (("TiN", "C"), ("TiN", "Si"), ("TiN", "Ge"))
T_HIGH = {"TiN/C": 900, "TiN/Si": 700, "TiN/Ge": 700}
# tbc-tin-v1 stdout (orb-v3): (G_bte, G_dmm) in MW m^-2 K^-1, printed rounded to integers
PRINTED = {"TiN/C": {300: (43.0, 264.0), 900: (53.0, 362.0)},
           "TiN/Si": {300: (52.0, 315.0), 700: (61.0, 340.0)},
           "TiN/Ge": {300: (51.0, 272.0), 700: (59.0, 288.0)}}
# relaxed cubic lattice constants (Angstrom) printed by tbc-tin-v1 for orb-v3
CELL_A = {"TiN": 4.2538, "C": 3.5721, "Si": 5.4608, "Ge": 5.7581}
NODES = 400

# pre-registered tolerances (amendment 0.5.7, section 4)
D0_ABS, D0_REL = 1.0, 0.02
D1_REL = 0.02
D2_REL = 0.01


def fcc_prim_volume_m3(a_ang):
    """Volume of the fcc primitive cell (rocksalt and diamond both have it), m^3."""
    return (a_ang * 1e-10) ** 3 / 4.0


def uniform_lengths(K, kappa_total):
    """Same two-stream length on every band with K > 0, chosen so that sum 2 K l = kappa_total."""
    out = []
    tot = sum(k for k in K if k > 0.0)
    for k in K:
        out.append(kappa_total / (2.0 * tot) if (k > 0.0 and tot > 0.0) else 0.0)
    return out


def shared_share(K_self, l_self, K_other):
    """Share of this side's conductivity (sum 2 K l) carried by bands the other side also has."""
    num = den = 0.0
    for k, l, ko in zip(K_self, l_self, K_other):
        if k > 0.0:
            den += 2.0 * k * l
            if ko > 0.0:
                num += 2.0 * k * l
    return num / den if den > 0.0 else float("nan")


def mean_length(K, l, K_other=None):
    """Conductivity-weighted mean length (all bands, or only bands shared with K_other)."""
    num = den = 0.0
    for i, (k, x) in enumerate(zip(K, l)):
        if k > 0.0 and (K_other is None or K_other[i] > 0.0):
            num += 2.0 * k * x * x
            den += 2.0 * k * x
    return num / den if den > 0.0 else float("nan")


def within(v, ref, abs_tol, rel_tol):
    return abs(v - ref) <= max(abs_tol, rel_tol * abs(ref))


def grade(rows):
    """rows[pair][T] = dict of G values in MW m^-2 K^-1. Returns {D0..D4: text}."""
    g = {}
    bad = []
    for p, by_t in PRINTED.items():
        for T, (gb, gd) in by_t.items():
            r = rows.get(p, {}).get(T)
            if r is None:
                bad.append(f"{p} {T} missing")
                continue
            if not within(r["base"], gb, D0_ABS, D0_REL):
                bad.append(f"{p} {T} G {r['base']:.1f} vs {gb:.0f}")
            if not within(r["dmm"], gd, D0_ABS, D0_REL):
                bad.append(f"{p} {T} DMM {r['dmm']:.1f} vs {gd:.0f}")
    g["D0"] = "PASS" if not bad else "FAIL(" + "; ".join(bad) + ")"
    allr = [(p, T, r) for p, by_t in rows.items() for T, r in by_t.items()]
    if not allr:
        for k in ("D1", "D2", "D3", "D4"):
            g[k] = "UNDETERMINED(no rows)"
        return g
    b1 = [f"{p} {T}" for p, T, r in allr
          if not (abs(r["x0.5"] / r["base"] - 1.0) <= D1_REL and abs(r["x2"] / r["base"] - 1.0) <= D1_REL)]
    g["D1"] = "PASS" if not b1 else "FAIL(" + ", ".join(b1) + ")"
    b2 = [f"{p} {T}" for p, T, r in allr if not abs(r["n800"] / r["base"] - 1.0) <= D2_REL]
    g["D2"] = "PASS" if not b2 else "FAIL(" + ", ".join(b2) + ")"
    b3 = [f"{p} {T} {r['uni']:.0f}<{r['dmm']:.0f}" for p, T, r in allr if not r["uni"] >= r["dmm"]]
    g["D3"] = ("PASS" if not b3 else "FAIL(" + ", ".join(b3) + ")") + "(" + ", ".join(
        f"{p} {T} uni/dmm {r['uni'] / r['dmm']:.2f}" for p, T, r in allr) + ")"
    parts = []
    for p, T, r in allr:
        gap = math.log(r["uni"] / r["base"]) if r["uni"] > 0 and r["base"] > 0 else float("nan")
        s1 = math.log(r["uni1"] / r["base"]) / gap if gap and gap == gap and gap != 0 else float("nan")
        s2 = math.log(r["uni2"] / r["base"]) / gap if gap and gap == gap and gap != 0 else float("nan")
        side = "TiN" if s1 > s2 else "partner"
        parts.append(f"{p} {T} side1 {s1:.2f} side2 {s2:.2f} -> {side}")
    g["D4"] = "INFO(" + "; ".join(parts) + ")"
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
        vgeo = fcc_prim_volume_m3(CELL_A[m])
        f = r["freq"]
        log(f"MAT {m}: vcell from kappa {r['vcell']:.4e} m^3, a^3/4 {vgeo:.4e} (ratio "
            f"{r['vcell'] / vgeo:.4f}); kappa300 ph3 {k_ph3:.1f} W/mK; anh-mesh freq max "
            f"{float(f.max()):.2f} THz, min {float(f.min()):.3f}, imag share "
            f"{float(np.mean(f < -0.05)):.4f}; M_or shape {tuple(r['M_or'].shape)}; temps "
            f"{','.join(str(int(x)) for x in temps)}")
        mats[m] = (r, temps)

    rows = {}
    for A, B in PAIRS:
        pair = f"{A}/{B}"
        ra, ta_temps = mats[A]
        rb, tb_temps = mats[B]
        rows[pair] = {}
        for T in (300, T_HIGH[pair]):
            ia, ib = ta_temps.index(float(T)), tb_temps.index(float(T))
            Ka = ta.K_bins(ra["M_or"], ra["edges_f"], ra["edges_c"], float(T), classical=False)
            Kb = ta.K_bins(rb["M_or"], rb["edges_f"], rb["edges_c"], float(T), classical=False)
            kap_a = ta.kappa_bins(ra, ia, float(T), ra["edges_c"], "iso", with_isotope=True)
            kap_b = ta.kappa_bins(rb, ib, float(T), rb["edges_c"], "iso", with_isotope=True)
            la, fa = tb.lengths_from_kappa(kap_a, Ka)
            lb, fb = tb.lengths_from_kappa(kap_b, Kb)

            def solve(l1, l2, n=NODES):
                return tb.solve_two_stream(Ka, Kb, l1, l2, n_nodes=n)

            base = solve(la, lb)
            ua = np.array(uniform_lengths(list(Ka), base["kappa1"]))
            ub = np.array(uniform_lengths(list(Kb), base["kappa2"]))
            res = {"base": base, "uni": solve(ua, ub), "uni1": solve(ua, lb), "uni2": solve(la, ub),
                   "x0.5": solve(0.5 * la, 0.5 * lb), "x2": solve(2.0 * la, 2.0 * lb),
                   "n800": solve(la, lb, 2 * NODES),
                   # information only: TiN phonon lengths x 0.01, a crude stand-in for TiN phonons
                   # held near equilibrium by electrons (not modelled), section 4 of 0.5.7
                   "s1x0.01": solve(0.01 * la, lb)}
            r = {k: v["G"] / 1e6 for k, v in res.items()}
            r["dmm"] = base["G_dmm"] / 1e6
            r["simons"] = base["G_simons"] / 1e6
            r["linear_bath"] = base["G_linear_bath"] / 1e6
            r["uni_simons"] = res["uni"]["G_simons"] / 1e6
            r["kappa1"], r["kappa2"] = base["kappa1"], base["kappa2"]
            r["residual"] = max(v["residual"] for v in res.values())
            r["capped"] = [base["capped_share1"], base["capped_share2"]]
            r["filled"] = [fa, fb]
            Kal, Kbl, lal, lbl = list(Ka), list(Kb), list(la), list(lb)
            r["shared1"] = shared_share(Kal, lal, Kbl)
            r["shared2"] = shared_share(Kbl, lbl, Kal)
            r["len1_nm"] = [mean_length(Kal, lal) * 1e9, mean_length(Kal, lal, Kbl) * 1e9]
            r["len2_nm"] = [mean_length(Kbl, lbl) * 1e9, mean_length(Kbl, lbl, Kal) * 1e9]
            r["sumK"] = [float(np.sum(Ka)) / 1e6, float(np.sum(Kb)) / 1e6]
            rows[pair][T] = r
            log(f"DIAG {pair} {T} K | G base {r['base']:.1f} dmm {r['dmm']:.1f} simons {r['simons']:.1f} "
                f"linbath {r['linear_bath']:.1f} | uni {r['uni']:.1f} uni1 {r['uni1']:.1f} uni2 "
                f"{r['uni2']:.1f} s1x0.01 {r['s1x0.01']:.1f} | x0.5 {r['x0.5']:.1f} x2 {r['x2']:.1f} n800 {r['n800']:.1f} | "
                f"kappa {r['kappa1']:.1f}/{r['kappa2']:.1f} W/mK sumK {r['sumK'][0]:.0f}/{r['sumK'][1]:.0f} "
                f"| shared {r['shared1']:.3f}/{r['shared2']:.3f} | len nm all/shared "
                f"{r['len1_nm'][0]:.3g}/{r['len1_nm'][1]:.3g} and {r['len2_nm'][0]:.3g}/{r['len2_nm'][1]:.3g} "
                f"| capped {r['capped'][0]:.3f}/{r['capped'][1]:.3f} filled {fa}/{fb} resid {r['residual']:.2e}")
    g = grade(rows)
    for k, v in g.items():
        log(f"{k}: {v}")
    with open(os.path.join(args.out, "tdg-summary.json"), "w") as fh:
        json.dump({"version": VERSION, "grades": g,
                   "rows": {p: {str(T): v for T, v in by.items()} for p, by in rows.items()},
                   "seconds": round(time.time() - t0, 1)}, fh, indent=1, default=float)
    log("TDG: " + " | ".join(f"{k} {str(v).split('(')[0]}" for k, v in g.items()) + " | " + " ; ".join(
        f"{p} {T} base/dmm/uni {r['base']:.0f}/{r['dmm']:.0f}/{r['uni']:.0f}"
        for p, by in rows.items() for T, r in by.items()))
    with open(os.path.join(args.out, "run.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
