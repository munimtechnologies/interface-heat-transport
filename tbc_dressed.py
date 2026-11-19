# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""tbc_dressed: amendment 0.5.1 of interface-heat-transport. Pure Python (math only).

Why this module exists
----------------------
Amendment 0.5 bounded any three-phonon channel I at ZnO/GaN from the flat measured G(T), assuming
the elastic part E cannot fall with temperature. Two things weaken that bound:

1. The published error bars. Gaskins et al. (Nano Lett. 18, 7469, 2018) give 95 % confidence
   contours of a separate TDTR fit at each temperature. Treated as independent, they allow a far
   larger channel than the reading error does (`profile_upper_limit`, `two_point_upper_limit`).

2. The conductance a TDTR fit reports is referred to the local (diffusive) temperatures on each
   side, i.e. the "dressed" conductance (Simons 1974; Zeng and Chen 2001; Shi, Yang, Fisher and
   Ruan, arXiv:1812.07910). Shi et al. state that their dressed ZnO/GaN conductance falls at high
   temperature. For a DMM interface the dressed conductance is exactly
       1/G_d = 1/G_0 - 1/(2 R_1) - 1/(2 R_2),
   with G_0 the conventional Landauer-DMM conductance and R_i the phonon radiation limits of the
   two sides (`dressed_dmm` computes it both ways: Shi et al.'s weighted equivalent-temperature
   form and the closed form). The radiation limit of the stiffer side keeps growing after the
   softer side has saturated, so G_d can fall while every mode's heat capacity grows
   (`test_tbc_dressed.TestDressed.test_two_bin_fall`). The premise of amendment 0.5 then fails.

What replaces it
----------------
If the elastic part changes by a factor kappa_E = E(T2)/E(T1) (from a model) and I grows by at
least rho = T2/T1:  G2 - G1 >= (kappa_E - 1) E1 + (rho - 1) I1 with E1 = G1 - I1, so
    I(T1) <= max(0, G2 - kappa_E G1) / (rho - kappa_E)          (kappa_E < rho).
kappa_E = 1 gives back amendment 0.5. `grade_supp` applies it with the stage 3a model's own
G_BTE(500)/G_BTE(300) (supplementary grades Z1g, Z6, Z7, filed before any output of tbc-zno-v1).
"""
import math

import tbc_nonmetal as tn

VERSION = "mv2 tbc-dressed 0.1.0"

# Blind reading of the published Nano Lett. 2018 figure (data/read-gaskins2018.json), MW m^-2 K^-1:
# T: (G, upper bar, lower bar). The bars are 95 % confidence contours of each temperature's fit.
PUB = {300: (490.0, 107.0, 111.0), 385: (524.0, 121.0, 124.0), 500: (503.0, 127.0, 131.0)}
PUB_400_INTERP = (521.0, 122.0, 125.0)
PUB_READ_ERR = 10.0          # reading uncertainty on points
PUB_DMM = {300: 250.0, 500: 280.0}
PUB_AGF = {300: 300.0, 500: 305.0}
PUB_TEXT_300 = (490.0, 150.0, 110.0)   # the text's room-temperature value

Z95_TWO_SIDED = 1.959963984540054
Z95_ONE_SIDED = 1.6448536269514722
DCHI2_95_ONE_SIDED = Z95_ONE_SIDED ** 2   # 2.7055

GAN_TOP_THZ = 22.0           # GaN A1(LO), about 735 cm^-1: the highest phonon frequency at ZnO/GaN
DRESS_RATIO = 2.0            # G_d / G_0 for a DMM interface well below both radiation limits


def sigma_from_95(up, down):
    """1-sigma from a two-sided 95 % bar (mean of the two arms)."""
    return 0.5 * (up + down) / Z95_TWO_SIDED


def isotonic(y, w):
    """Weighted least-squares non-decreasing fit (pool-adjacent-violators)."""
    blocks = []
    for yi, wi in zip(y, w):
        blocks.append([float(yi), float(wi), 1])
        while len(blocks) > 1 and blocks[-2][0] > blocks[-1][0]:
            m2, w2, n2 = blocks.pop()
            m1, w1, n1 = blocks.pop()
            wt = w1 + w2
            blocks.append([(m1 * w1 + m2 * w2) / wt, wt, n1 + n2])
    fit = []
    for m, _, n in blocks:
        fit.extend([m] * n)
    return fit


def chi2_profile(a, temps, g, sig, T0=300.0):
    """chi^2 of G(T) = E(T) + a T/T0 with E non-decreasing, profiled over E."""
    r = [gi - a * T / T0 for gi, T in zip(g, temps)]
    w = [1.0 / s ** 2 for s in sig]
    fit = isotonic(r, w)
    return sum(wi * (ri - fi) ** 2 for wi, ri, fi in zip(w, r, fit))


def profile_upper_limit(temps, g, sig, dchi2=DCHI2_95_ONE_SIDED, amax=5000.0):
    """Approximate 95 % upper limit on I(T0) = a for a classical channel I = a T/T0 (a >= 0), with
    independent Gaussian errors. chi2_profile is convex in a, so a ternary search finds the minimum
    and a bisection the crossing. Returns (a_upper, chi2_min, a_at_min)."""
    lo, hi = 0.0, amax
    for _ in range(200):
        m1, m2 = lo + (hi - lo) / 3.0, hi - (hi - lo) / 3.0
        if chi2_profile(m1, temps, g, sig) <= chi2_profile(m2, temps, g, sig):
            hi = m2
        else:
            lo = m1
    amin = 0.5 * (lo + hi)
    if chi2_profile(0.0, temps, g, sig) <= chi2_profile(amin, temps, g, sig):
        amin = 0.0
    cmin = chi2_profile(amin, temps, g, sig)
    target = cmin + dchi2
    lo, hi = amin, amax
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if chi2_profile(mid, temps, g, sig) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), cmin, amin


def two_point_upper_limit(g1, s1, g2, s2, T1, T2, z=Z95_ONE_SIDED):
    """One-sided 95 % version of amendment 0.5's bound from two independent points."""
    return max(0.0, g2 - g1 + z * math.hypot(s1, s2)) / (T2 / T1 - 1.0)


def bound_general(g1, g2, T1, T2, kappa_e, growth=1.0):
    """I(T1) <= max(0, G2 - kappa_E G1) / (rho - kappa_E); inf if kappa_E >= rho."""
    r = tn.rho(T1, T2, growth)
    if kappa_e >= r:
        return float("inf")
    return max(0.0, g2 - kappa_e * g1) / (r - kappa_e)


def cap_ratio(nu_thz, T1, T2):
    """C(nu, T2) / C(nu, T1); grows with nu (unit test), so the top frequency gives the maximum."""
    return tn.cap(tn.x_of(nu_thz, T2)) / tn.cap(tn.x_of(nu_thz, T1))


def min_kappa_e_dmm(nu_top_thz, T1, T2, dress_ratio=DRESS_RATIO):
    """Crude floor on G_d(T2)/G_d(T1) for a DMM interface: G_0 held fixed, both radiation limits
    grown by the largest heat-capacity ratio kappa of any mode. From the closed form,
    G_d(T1)/G_d(T2) <= 1 + (1 - 1/kappa)(G_d/G_0 - 1)."""
    k = cap_ratio(nu_top_thz, T1, T2)
    return 1.0 / (1.0 + (1.0 - 1.0 / k) * (dress_ratio - 1.0))


def dressed_dmm(M1, M2, C):
    """Per frequency bin: numbers of modes M1, M2 (side 1, side 2) and mode heat capacity C.
    DMM transmission t12 = M2/(M1+M2), t21 = M1/(M1+M2) where both sides have modes, else 0.
    Returns G_0 (conventional), R_1, R_2 (radiation limits), G_d from Shi et al.'s weighted
    equivalent-equilibrium temperatures (weights C M_i), and G_d from the closed form."""
    g0 = r1 = r2 = 0.0
    n12 = n21 = 0.0
    for m1, m2, c in zip(M1, M2, C):
        r1 += c * m1
        r2 += c * m2
        if m1 > 0 and m2 > 0:
            t12 = m2 / (m1 + m2)
            t21 = m1 / (m1 + m2)
            g0 += c * m1 * t12
            n12 += c * m1 * t12
            n21 += c * m2 * t21
    mean_t12 = n12 / r1
    mean_t21 = n21 / r2
    g_shi = g0 / (1.0 - 0.5 * mean_t12 - 0.5 * mean_t21)
    g_closed = 1.0 / (1.0 / g0 - 0.5 / r1 - 0.5 / r2)
    return {"G0": g0, "R1": r1, "R2": r2, "Gd_shi": g_shi, "Gd_closed": g_closed}


def numbers():
    """Every number of PLAN-amendment-0.5.1 section 2 (results/nonmetal/bounds-0.5.1.json)."""
    temps = sorted(PUB)
    g = [PUB[T][0] for T in temps]
    sig = [sigma_from_95(PUB[T][1], PUB[T][2]) for T in temps]
    a_up, cmin, amin = profile_upper_limit(temps, g, sig)
    s300, s500 = sig[0], sig[-1]
    tp = two_point_upper_limit(PUB[300][0], s300, PUB[500][0], s500, 300.0, 500.0)
    k_top = cap_ratio(GAN_TOP_THZ, 300.0, 500.0)
    ke_min = min_kappa_e_dmm(GAN_TOP_THZ, 300.0, 500.0)
    return {
        "sigma_1s": {str(T): s for T, s in zip(temps, sig)},
        "profile_I300_upper95": a_up,
        "profile_chi2_min": cmin,
        "profile_a_at_min": amin,
        "profile_share300": a_up / PUB[300][0],
        "two_point_I300_upper95": tp,
        "correlated_I300_central_pub": tn.channel_bound(PUB[300][0], PUB[500][0], 300, 500),
        "correlated_I300_reading_worst_pub": tn.channel_bound(PUB[300][0] - PUB_READ_ERR,
                                                              PUB[500][0] + PUB_READ_ERR, 300, 500),
        "pub_minus_arxiv": {"300": PUB[300][0] - tn.ZNO_GAN[300], "400": PUB_400_INTERP[0] - tn.ZNO_GAN[400],
                            "500": PUB[500][0] - tn.ZNO_GAN[500]},
        "dmm_pub_minus_arxiv_300": PUB_DMM[300] - tn.ZNO_GAN_DMM[300],
        "gan_top_cap_ratio_300_500": k_top,
        "kappa_e_floor_dmm_crude": ke_min,
        "I300_bound_at_crude_floor_central_pub": bound_general(PUB[300][0], PUB[500][0], 300, 500, ke_min),
        "twice_dmm_pub": {"300": 2 * PUB_DMM[300], "500": 2 * PUB_DMM[500]},
    }


# ----------------------------------------------------------------------------------------------
# Supplementary grades of tbc-zno-v1 (PLAN-amendment-0.5.1 section 3), on the job's primary model,
# read from tbz-summary.json["results"][primary]. They add to Z0-Z5 and change none of them.
_W300, _W500 = tn.worst_case(tn.ZNO_GAN[300], tn.ZNO_GAN[500], tn.ZNO_GAN_REL)   # 471.36, 522.08
Z7_RANGE = (1.5, 2.5)


def grade_supp(prim):
    g = {}
    try:
        rows = prim["pairs"]["ZnO/GaN"]
        b300, b500 = float(rows["300"]["G_bte"]), float(rows["500"]["G_bte"])
        i300 = float(rows["300"]["G_inel"])
        ke = b500 / b300
        bg = bound_general(_W300, _W500, 300, 500, ke)
        g["Z1g"] = ("PASS" if i300 <= bg else "FAIL") + \
            f"(G_inel300 {i300:.1f} vs max {bg:.2f}; kappa_E {ke:.4f})"
    except (KeyError, TypeError, ZeroDivisionError) as e:
        g["Z1g"] = f"UNDETERMINED(missing {e})"
    try:
        rows = prim["pairs"]["ZnO/GaN"]
        b400, b500 = float(rows["400"]["G_bte"]), float(rows["500"]["G_bte"])
        g["Z6"] = ("PASS" if b500 < b400 else "FAIL") + f"(G_bte 400 {b400:.0f}, 500 {b500:.0f})"
    except (KeyError, TypeError) as e:
        g["Z6"] = f"UNDETERMINED(missing {e})"
    try:
        rows = prim["pairs"]["ZnO/GaN"]
        r = float(rows["300"]["G_bte"]) / float(rows["300"]["G_dmm"])
        g["Z7"] = ("PASS" if Z7_RANGE[0] <= r <= Z7_RANGE[1] else "FAIL") + \
            f"(G_bte/G_dmm at 300 K {r:.3f} vs [{Z7_RANGE[0]}, {Z7_RANGE[1]}])"
    except (KeyError, TypeError, ZeroDivisionError) as e:
        g["Z7"] = f"UNDETERMINED(missing {e})"
    return g
