# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Made with AI agents working in MV2's name.
"""tbc_direct: a direct electron-to-substrate-phonon channel at a metal/non-metal interface.

Circuit (steady state, semi-infinite metal, two-temperature model). Electrons (k_e) and phonons
(k_p) of the metal exchange heat through the bulk coupling g. At the interface the electrons pass
heat straight to the substrate through sigma_e, and the metal phonons through sigma_p. The
effective conductance is defined with the far-field (equilibrated) metal temperature extrapolated
to the interface, as a TDTR fit on a film much thicker than the coupling length does. Solving the
two coupled equations gives (P = k_e/l, l = sqrt(k_e k_p / (g (k_e + k_p))), r = k_e/k_p,
f_e = k_e/(k_e+k_p), f_p = k_p/(k_e+k_p)):

    G = [P (sigma_p + sigma_e) + sigma_p sigma_e (1 + r)] / [P + f_p sigma_e + f_e r sigma_p]

For sigma_e = 0 this is the Majumdar-Reddy series form, 1/G = 1/sigma_p + f_e^(3/2)/sqrt(g k_p).
For k_p << k_e it becomes sigma_e + series(sqrt(g k_p), sigma_p): the direct channel in parallel
with the usual series path. This is the mixed series-parallel circuit of Lu, Wang and Ruan
(Phys. Rev. B 93, 064302, 2016), which we re-derived; credit is theirs.

Sergeev's model of the direct channel (inelastic electron-boundary scattering; A. Sergeev,
Phys. Rev. B 58, R10199, 1998), as written by Lombard, Detcheverry and Merabia (arXiv:1501.03176,
Eqs. 1-3):
    sigma_e = [3 pi hbar / (35 zeta(3) k_B)] (gamma u_L / tau_ep) [1 + 2 (u_L/u_T)^3]
    G_bulk  = pi^2 m n_e c_s^2 / (6 tau_ep T)          (we take c_s = u_L)
    G_bulk  = 3 gamma (hbar^2 lambda<omega^2>) / (pi hbar k_B)
It depends only on bulk properties of the metal (not on the substrate or the interface) and grows
in proportion to T. Their printed Eq. 4 drops zeta(3) and has c_s where we get u_L, and their
Table 1 values do not all follow from it, so we evaluate the model three ways (see run file).

numpy is not needed. SI units throughout (W m^-2 K^-1 for conductances; s in W m^-2 K^-2).
"""
import math

from tbc_gep import gep_power_law, series, phonon_side

VERSION = "mv2 tbc-direct 0.1.0"

HBAR = 1.054571817e-34
KB = 1.380649e-23
ME = 9.1093837015e-31
MEV = 1.602176634e-22
ZETA3 = 1.2020569031595942


def circuit_exact(g, k_e, k_p, sigma_e, sigma_p):
    """Effective conductance of the two-temperature circuit with a direct channel."""
    if g <= 0 or k_e <= 0 or k_p <= 0 or sigma_e < 0 or sigma_p < 0:
        raise ValueError("g, k_e, k_p must be positive; sigma_e, sigma_p non-negative")
    ell = math.sqrt(k_e * k_p / (g * (k_e + k_p)))
    p = k_e / ell
    r = k_e / k_p
    fe = k_e / (k_e + k_p)
    fp = k_p / (k_e + k_p)
    num = p * (sigma_p + sigma_e) + sigma_p * sigma_e * (1.0 + r)
    den = p + fp * sigma_e + fe * r * sigma_p
    return num / den


def circuit_limit(gep, sigma_e, sigma_p):
    """k_p << k_e limit: sigma_e in parallel with series(G_ep, sigma_p)."""
    if sigma_p == 0:
        return sigma_e
    return sigma_e + series(gep, sigma_p)


def g_from_lambda(gamma_vol, lam_w2_mev2):
    """Bulk electron-phonon coupling g = 3 gamma lambda<(hbar omega)^2> / (pi hbar k_B)."""
    return 3.0 * gamma_vol * lam_w2_mev2 * MEV ** 2 / (math.pi * HBAR * KB)


def _bracket(u_l, u_t):
    return 1.0 + 2.0 * (u_l / u_t) ** 3


def sergeev_from_g(gamma_vol, g, u_l, u_t, n_e, T, m=ME):
    """Sergeev's sigma_e with tau_ep from g (c_s = u_L):
    sigma_e = 18 hbar gamma g T [1 + 2(u_L/u_T)^3] / (35 zeta(3) pi k_B m n_e u_L)."""
    return (18.0 * HBAR * gamma_vol * g * T * _bracket(u_l, u_t)
            / (35.0 * ZETA3 * math.pi * KB * m * n_e * u_l))


def sergeev_from_lambda(gamma_vol, lam_w2_mev2, u_l, u_t, n_e, T, m=ME, zeta=True):
    """Sergeev's sigma_e with g from lambda<omega^2>. zeta=False reproduces the printed Eq. 4 of
    Lombard et al. (no zeta(3)), which is larger by zeta(3) = 1.202."""
    val = sergeev_from_g(gamma_vol, g_from_lambda(gamma_vol, lam_w2_mev2), u_l, u_t, n_e, T, m)
    return val if zeta else val * ZETA3


def _pairs(meas, temps):
    ts = [t for t in temps if str(t) in meas or t in meas]
    return ts, [float(meas.get(str(t), meas.get(t))) for t in ts]


def linear_bound(meas, temps=(300, 400, 500, 600, 700)):
    """Largest slope s of a direct channel sigma_e = s T with sigma_e(T) < G(T) at every T
    (any non-negative series path). Returns (s_max, T where it binds)."""
    ts, gs = _pairs(meas, temps)
    vals = [(g / t, t) for t, g in zip(ts, gs)]
    return min(vals)


def floor_bound(meas, gep400, n, floor, temps=(300, 400, 500, 600, 700)):
    """Largest s if the phonon side is at least `floor` at every T (a phonon side that does not
    fall below an ab initio harmonic value): s T <= G(T) - series(G_ep(T), floor)."""
    ts, gs = _pairs(meas, temps)
    vals = [((g - series(gep_power_law(gep400, t, 400.0, n), floor)) / t, t)
            for t, g in zip(ts, gs)]
    return min(vals)


def infer_with_direct(meas, gep400, n, s, temps=(300, 400, 500, 600, 700)):
    """Phonon side implied by G(T) when a direct channel s T runs in parallel:
    G_pp = phonon_side(G - s T, G_ep). None where G - s T <= 0."""
    ts, gs = _pairs(meas, temps)
    out = []
    for t, g in zip(ts, gs):
        rest = g - s * t
        gep = gep_power_law(gep400, t, 400.0, n)
        out.append(phonon_side(rest, gep) if 0 < rest < gep else None)
    return {"T": ts, "G_meas": gs, "sigma_e": [s * t for t in ts], "G_pp": out}


def min_rise_prediction(g1, t1, t2, s, gep400, n):
    """Lowest G(T2) the circuit allows given G(T1), a direct channel s T, and a phonon side that
    does not fall between T1 and T2: s T2 + series(G_ep(T2), phonon_side(G1 - s T1, G_ep(T1)))."""
    gep1 = gep_power_law(gep400, t1, 400.0, n)
    gep2 = gep_power_law(gep400, t2, 400.0, n)
    rest = g1 - s * t1
    if rest <= 0:
        return {"G_min": None, "reason": "G(T1) below the direct channel alone"}
    spp = phonon_side(rest, gep1)
    return {"G_min": s * t2 + series(gep2, spp), "sigma_p_T1": spp, "sigma_e_T2": s * t2}
