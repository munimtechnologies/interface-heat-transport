# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""tbc_nonmetal: what measured G(T) at non-metal interfaces allows for an inelastic channel, and the
pre-registered grades Z0-Z5 of stage 3a (PLAN-amendment-0.5.md). Pure Python (math only).

The bound
---------
Write the measured conductance as G(T) = E(T) + I(T), with
- E the elastic part. With fixed harmonic force constants it can only grow with T (every mode's
  heat capacity grows), so E(T2) >= E(T1) for T2 > T1;
- I an inelastic (three-phonon) channel. Its T-dependence is that of C(nu, T) * Gamma(nu, T) for
  every mode, summed with T-independent weights. Gamma is a sum of decay terms (1 + n' + n'') and
  combination terms (n' - n''). For T2 > T1 each such term, times the mode's heat capacity, grows
  by at least the classical factor rho = T2 / T1 (checked on a dense grid in the unit tests; in the
  classical limit it is exactly T2 / T1). Hence I(T2) >= rho I(T1).
Subtracting, G(T2) - G(T1) >= (rho - 1) I(T1), so

    I(T1) <= max(0, G(T2) - G(T1)) / (rho - 1).

`growth` < 1 replaces rho by 1 + growth (rho - 1), a cautious variant for channels that grow more
slowly than the classical three-phonon rates.
"""
import math

VERSION = "mv2 tbc-nonmetal 0.1.0"

H_OVER_KB_THZ = 47.99243073366221   # h * 1 THz / k_B in K

# Readings (data/read-nonmetal.json); the test checks they match the file.
ZNO_GAN = {300: 491.0, 400: 520.0, 500: 502.0, 100: 201.0, 200: 366.0}
ZNO_GAN_REL = 0.04
ZNO_GAN_AGF = {300: 300.0, 500: 305.0}
ZNO_GAN_DMM = {300: 240.0, 500: 278.0}
ZNO_GAN_RT_RANGE = (490.0 - 110.0, 490.0 + 150.0)          # abstract: 490 [+150, -110]
UCLA_GAN_ALN_300 = 324.0
WALWIL = {"GaN/AlN": {300: 483.0, 400: 382.0}, "AlN/SiC": {300: 475.0, 350: 551.0},
          "GaN/SiC": {300: 307.0, 350: 373.0, 400: 408.0, 450: 447.0}}

# ZnO sanity references for Z0: highest non-polar-corrected mode E2(high) 438 cm^-1 = 13.1 THz;
# bulk thermal conductivity at 300 K roughly 40-100 W/m/K in the literature, range kept wide.
ZNO_REF_MAX_THZ = 13.1
ZNO_KAPPA_RANGE = (20.0, 150.0)


def occ(x):
    """Bose occupation 1 / (e^x - 1)."""
    return 1.0 / math.expm1(x)


def cap(x):
    """Mode heat capacity in units of k_B: x^2 e^x / (e^x - 1)^2."""
    if x < 1e-6:
        return 1.0
    e = math.expm1(x)
    return x * x * (e + 1.0) / (e * e)


def x_of(nu_thz, T):
    return H_OVER_KB_THZ * nu_thz / T


def growth_decay(nu_thz, frac, T1, T2):
    """[C Gamma](T2) / [C Gamma](T1) for a mode nu decaying into frac*nu + (1-frac)*nu."""
    def term(T):
        return cap(x_of(nu_thz, T)) * (1.0 + occ(x_of(frac * nu_thz, T)) +
                                       occ(x_of((1.0 - frac) * nu_thz, T)))
    return term(T2) / term(T1)


def growth_combination(nu_thz, nu1_thz, T1, T2):
    """Same for a mode nu combining with nu1 into nu + nu1 (rate ~ n(nu1) - n(nu + nu1))."""
    def term(T):
        return cap(x_of(nu_thz, T)) * (occ(x_of(nu1_thz, T)) - occ(x_of(nu_thz + nu1_thz, T)))
    return term(T2) / term(T1)


def rho(T1, T2, growth=1.0):
    return 1.0 + growth * (T2 / T1 - 1.0)


def channel_bound(g1, g2, T1, T2, growth=1.0):
    """Largest I(T1) allowed by G(T1) = g1, G(T2) = g2 (module docstring)."""
    return max(0.0, g2 - g1) / (rho(T1, T2, growth) - 1.0)


def worst_case(g1, g2, rel):
    """Readings moved by the reading error in the direction that loosens the bound."""
    return g1 * (1.0 - rel), g2 * (1.0 + rel)


def exploratory():
    """Data-only numbers of PLAN-amendment-0.5 section 2 (MW m^-2 K^-1)."""
    z = ZNO_GAN
    out = {}
    lo300, hi500 = worst_case(z[300], z[500], ZNO_GAN_REL)
    lo400, hi500b = worst_case(z[400], z[500], ZNO_GAN_REL)
    b35 = channel_bound(z[300], z[500], 300, 500)
    b35w = channel_bound(lo300, hi500, 300, 500)
    b45 = channel_bound(z[400], z[500], 400, 500)
    b45w = channel_bound(lo400, hi500b, 400, 500)
    out["zno_gan"] = {
        "I300_max_central_300_500": b35,
        "share300_central": b35 / z[300],
        "I300_max_worst_300_500": b35w,
        "share300_worst": b35w / lo300,
        "I400_max_central_400_500": b45,
        "I400_max_worst_400_500": b45w,
        "I300_max_worst_from_400_500": b45w / rho(300, 400),
        "I300_max_worst_combined": min(b35w, b45w / rho(300, 400)),
        "I300_max_central_half_growth": channel_bound(z[300], z[500], 300, 500, growth=0.5),
        "I300_max_worst_half_growth": channel_bound(lo300, hi500, 300, 500, growth=0.5),
        "ratio_500_300_central": z[500] / z[300],
        "ratio_500_300_worst": hi500 / lo300,
        "agf_gap300": z[300] - ZNO_GAN_AGF[300],
        "G500_min_if_agf_gap_is_channel": ZNO_GAN_AGF[500] + (z[300] - ZNO_GAN_AGF[300]) * rho(300, 500),
        "dmm_gap300": z[300] - ZNO_GAN_DMM[300],
        "G500_min_if_dmm_gap_is_channel": ZNO_GAN_DMM[500] + (z[300] - ZNO_GAN_DMM[300]) * rho(300, 500),
    }
    out["walwil_ratios"] = {
        "GaN/AlN_400_300": WALWIL["GaN/AlN"][400] / WALWIL["GaN/AlN"][300],
        "AlN/SiC_350_300": WALWIL["AlN/SiC"][350] / WALWIL["AlN/SiC"][300],
        "GaN/SiC_450_300": WALWIL["GaN/SiC"][450] / WALWIL["GaN/SiC"][300],
    }
    out["ucla_over_agf_300"] = UCLA_GAN_ALN_300 / 346.0
    return out


# ----------------------------------------------------------------------------------------------
# Pre-registered grades of stage 3a (PLAN-amendment-0.5.md section 3), on the job's primary model.
_lo300, _hi500 = worst_case(ZNO_GAN[300], ZNO_GAN[500], ZNO_GAN_REL)
Z1_MAX = channel_bound(_lo300, _hi500, 300, 500)          # 76.08
Z2_MAX = _hi500 / _lo300                                   # 1.10760
Z3_RANGE = ZNO_GAN_RT_RANGE                                # (380, 640)
Z4_RANGE = (0.75 * UCLA_GAN_ALN_300, 1.25 * UCLA_GAN_ALN_300)   # (243, 405)


def _row(prim, pair, T, key):
    return float(prim["pairs"][pair][str(T)][key])


def grade_z(prim):
    """prim: {"materials": {"ZnO": {"relax_converged", "imag_fraction", "max_thz", "kappa300"}},
    "pairs": {"ZnO/GaN" | "AlN/GaN": {"300": {"G_bte", "G_inel", "G_total"}, ...}}}, MW m^-2 K^-1."""
    g = {}
    z = prim.get("materials", {}).get("ZnO")
    if not z:
        g["Z0"] = "UNDETERMINED(ZnO missing)"
    else:
        ok = (bool(z["relax_converged"]) and z["imag_fraction"] < 0.01 and
              abs(z["max_thz"] / ZNO_REF_MAX_THZ - 1.0) <= 0.15 and
              ZNO_KAPPA_RANGE[0] <= z["kappa300"] <= ZNO_KAPPA_RANGE[1])
        g["Z0"] = ("PASS" if ok else "FAIL") + \
            f"(conv {bool(z['relax_converged'])}, imag {z['imag_fraction']:.4f}, max {z['max_thz']:.2f} THz, " \
            f"kappa300 {z['kappa300']:.1f})"
    try:
        v = _row(prim, "ZnO/GaN", 300, "G_inel")
        g["Z1"] = ("PASS" if v <= Z1_MAX else "FAIL") + f"(G_inel300 {v:.1f} vs max {Z1_MAX:.2f})"
    except (KeyError, TypeError) as e:
        g["Z1"] = f"UNDETERMINED(missing {e})"
    try:
        r = _row(prim, "ZnO/GaN", 500, "G_total") / _row(prim, "ZnO/GaN", 300, "G_total")
        g["Z2"] = ("PASS" if r <= Z2_MAX else "FAIL") + f"(ratio {r:.4f} vs max {Z2_MAX:.4f})"
    except (KeyError, TypeError, ZeroDivisionError) as e:
        g["Z2"] = f"UNDETERMINED(missing {e})"
    try:
        v = _row(prim, "ZnO/GaN", 300, "G_total")
        side = "in" if Z3_RANGE[0] <= v <= Z3_RANGE[1] else ("low" if v < Z3_RANGE[0] else "high")
        g["Z3"] = ("PASS" if side == "in" else "FAIL") + f"({side}: {v:.0f} vs [{Z3_RANGE[0]:.0f}, {Z3_RANGE[1]:.0f}])"
    except (KeyError, TypeError) as e:
        g["Z3"] = f"UNDETERMINED(missing {e})"
    try:
        v = _row(prim, "AlN/GaN", 300, "G_total")
        side = "in" if Z4_RANGE[0] <= v <= Z4_RANGE[1] else ("low" if v < Z4_RANGE[0] else "high")
        g["Z4"] = ("PASS" if side == "in" else "FAIL") + f"({side}: {v:.0f} vs [{Z4_RANGE[0]:.0f}, {Z4_RANGE[1]:.0f}])"
    except (KeyError, TypeError) as e:
        g["Z4"] = f"UNDETERMINED(missing {e})"
    try:
        r = _row(prim, "AlN/GaN", 400, "G_total") / _row(prim, "AlN/GaN", 300, "G_total")
        g["Z5"] = ("PASS" if r >= 1.0 else "FAIL") + f"(ratio {r:.4f}; Walwil measured 0.791)"
    except (KeyError, TypeError, ZeroDivisionError) as e:
        g["Z5"] = f"UNDETERMINED(missing {e})"
    return g
