# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Made with AI agents working in MV2's name.
"""tbc_gep: the electron-phonon side of a metal/non-metal interface, and what it does to the
phonon-side conductance inferred from a measured G(T).

In a metal on a non-metal, heat carried by the metal's electrons must first reach the metal's
phonons. Majumdar and Reddy (Appl. Phys. Lett. 84, 4768, 2004) give the resulting series
conductance as G_ep = sqrt(g * Lambda_p), with g the bulk electron-phonon coupling factor
(W m^-3 K^-1) and Lambda_p the metal's lattice thermal conductivity (W m^-1 K^-1). The measured
conductance is then G = [1/G_ep + 1/G_pp]^-1, and the phonon side follows as
G_pp = 1/(1/G - 1/G_ep).

Temperature dependence used here:
- g of Al is close to constant near equilibrium (Waldecker et al., PRX 6, 021003, 2016: weak
  dependence in DFT).
- Lambda_p of Al is three-phonon limited and falls as ~1/T above the Debye temperature;
  phonon-electron scattering leaves it almost unchanged in Al (Wang, Lu and Ruan, J. Appl. Phys.
  119, 225109, 2016). An optional T-independent phonon-electron term enters by Matthiessen's rule.
So G_ep ~ T^-1/2 (n = 0.5 in gep_power_law) is the bulk-input expectation; n = 0 is a constant G_ep.

The inferred G_pp amplifies relative errors of G by 1/(1 - G/G_ep) (see amplification()).

numpy is not needed. Units are SI throughout (W m^-2 K^-1, K).
"""
import math

VERSION = "mv2 tbc-gep 0.1.0"


def gep_mr(g, lam_p):
    """Majumdar-Reddy series conductance sqrt(g * Lambda_p), W m^-2 K^-1."""
    if g <= 0 or lam_p <= 0:
        raise ValueError("g and lam_p must be positive")
    return math.sqrt(g * lam_p)


def lam_p_T(lam_pp_ref, T, T_ref=300.0, lam_pe=None):
    """Lattice conductivity at T. The phonon-phonon part scales as T_ref/T; an optional
    T-independent phonon-electron limit lam_pe is added by Matthiessen's rule."""
    if T <= 0 or T_ref <= 0:
        raise ValueError("temperatures must be positive")
    lpp = lam_pp_ref * T_ref / T
    if lam_pe is None:
        return lpp
    return 1.0 / (1.0 / lpp + 1.0 / lam_pe)


def gep_power_law(gep_ref, T, T_ref=400.0, n=0.5):
    """G_ep(T) = G_ep(T_ref) * (T_ref/T)^n. n = 0.5 is Lambda_p ~ 1/T with constant g."""
    return gep_ref * (T_ref / T) ** n


def series(*gs):
    """Conductances in series."""
    return 1.0 / sum(1.0 / g for g in gs)


def phonon_side(G, gep):
    """G_pp = 1/(1/G - 1/G_ep). Needs G < G_ep."""
    if G >= gep:
        raise ValueError("measured G must be below G_ep")
    return 1.0 / (1.0 / G - 1.0 / gep)


def amplification(G, gep):
    """d ln G_pp / d ln G = 1/(1 - G/G_ep): how much a relative error in G grows in G_pp."""
    if G >= gep:
        raise ValueError("measured G must be below G_ep")
    return 1.0 / (1.0 - G / gep)


def diffs(vals):
    """First differences of a list."""
    return [b - a for a, b in zip(vals[:-1], vals[1:])]


def linear_floor(g400, g700):
    """For G_pp = E + a T fixed by two points at 400 and 700 K: (E, a, share of aT at 400 K)."""
    a = (g700 - g400) / 300.0
    e = (7.0 * g400 - 4.0 * g700) / 3.0
    return e, a, a * 400.0 / g400


def infer(meas_pair, gep400, n, temps=(300, 400, 500, 600, 700)):
    """Phonon side, its increments and second differences, and the error amplification, for
    one measured pair {T: G} and G_ep(T) = gep400 * (400/T)^n."""
    ts = [t for t in temps if str(t) in meas_pair or t in meas_pair]
    g = [float(meas_pair.get(str(t), meas_pair.get(t))) for t in ts]
    gep = [gep_power_law(gep400, t, 400.0, n) for t in ts]
    gpp = [phonon_side(a, b) for a, b in zip(g, gep)]
    inc = diffs(gpp)
    sec = diffs(inc)
    out = {
        "T": ts,
        "G_meas": g,
        "G_ep": gep,
        "G_pp": gpp,
        "increments": inc,
        "second_differences": sec,
        "convex_everywhere": bool(sec) and all(s > 0 for s in sec),
        "amplification": [amplification(a, b) for a, b in zip(g, gep)],
    }
    if 400 in ts and 700 in ts:
        e, a, share = linear_floor(gpp[ts.index(400)], gpp[ts.index(700)])
        out["linear_400_700"] = {"E": e, "a": a, "share_at_400": share}
    return out
