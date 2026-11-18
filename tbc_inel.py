# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""tbc_inel: an inelastic (three-phonon) channel across a solid-solid interface whose strength
comes only from bulk quantities: phono3py three-phonon rates, phonon densities of states and the
crystallography of the terminating atomic layers. Stage 2b of the interface-heat-transport study
(PLAN-amendment-0.3.md). numpy only; SI units unless a name says THz.

Model (interfacial anharmonic bridge, IAB)
------------------------------------------
The atoms of the last atomic layer on side i have a fraction f_i of their nearest-neighbour bonds
across the interface. On those bonds a three-phonon process can take energy from a phonon of side
i and leave it in a phonon of side j. We assume the cross bonds are as anharmonic as the bulk
bonds of side i, so a side-i phonon of frequency nu in that layer scatters across at the rate

    Gamma_i(nu, T) * Q_ij(nu),

where Gamma_i is its bulk three-phonon relaxation rate and Q_ij the share of its (weighted)
two-phonon phase space in which its energy ends up on side j, when each partner phonon on a cross
bond sits on either side with the per-atom density of states of that side:

    rho_x = (rho_i + rho_j) / 2                                    (per-atom DOS, each integrates to 1)
    D_i(nu)  = dec[rho_i, rho_i](nu) + 2 com[rho_i, rho_i](nu)     (bulk phase space of side i)
    X_ij(nu) = dec_E[rho_j, rho_x](nu) + com[rho_x, rho_j](nu)      (energy lands on side j)
    Q_ij = X_ij / D_i  (0 where D_i vanishes, capped at Q_CAP)

dec[a, b](nu) = int a(nu') b(nu - nu') w_d dnu' is the decay nu -> nu' + nu''; dec_E weights it by
the energy share nu'/nu of the product taken from a; com[a, b](nu) = int a(nu') b(nu + nu') w_c dnu'
is the combination nu + nu' -> nu''. With kernel 'classical' the weights are the classical-limit
occupation factors over the frequency product, w_d = (1/nu' + 1/nu'')/(nu' nu''),
w_c = (1/nu' - 1/nu'')/(nu' nu''); with kernel 'jdos' they are 1. For two identical sides Q = 1/2
exactly (unit test).

In linear response with a temperature jump dT across the interface, the energy that the cross
processes move per unit area and time is

    G_inel = 1/2 * sum_i t_i * sum_nu C_i(nu, T) Gamma_i(nu, T) Q_ij(nu),     t_i = f_i * sigma_i / n_i,

with C_i the volumetric spectral heat capacity, sigma_i the areal atom density of the terminating
layer and n_i the bulk atom density (so sigma_i / n_i is a length, a fraction of a layer spacing).
The prefactor 1/2 removes the double count: every cross event appears in the relaxation rate of a
side-1 phonon and of a side-2 phonon.

In the classical limit C is k_B per mode and Gamma grows in proportion to T, while Q does not
depend on T, so G_inel grows in proportion to T. This channel breaks the scale invariance of
tbc_bte (stage 2a): its length t_i is fixed by the lattice, not by the mean free paths. It adds in
parallel to the elastic two-stream conductance of tbc_bte.

Every number here is a bulk or crystallographic input. No interface conductance is fitted.
"""
import numpy as np

VERSION = "mv2 tbc-inel 0.1.0"

KB = 1.380649e-23
H = 6.62607015e-34
Q_CAP = 1.0
H_THZ = 0.1          # default frequency grid step (THz)

# terminating layer: fraction of nearest-neighbour bonds that cross, and sigma / n in units of the
# lattice constants (a for cubic, c for hexagonal); bulk-terminated, unreconstructed, ideal epitaxy.
#   fcc(111):       12 neighbours, 3 across; sigma = 4/(sqrt3 a^2), n = 4/a^3   -> sigma/n = a/sqrt3
#   diamond(111):   4 bonds, 1 across (single-bond termination); sigma = 4/(sqrt3 a^2), n = 8/a^3
#                   -> sigma/n = a/(2 sqrt3)
#   wurtzite(0001): 4 bonds, 1 across (along c); sigma = 2/(sqrt3 a^2), n = 8/(sqrt3 a^2 c)
#                   -> sigma/n = c/4
GEOMETRY = {
    "fcc111": dict(f=3.0 / 12.0, sigma_over_n=lambda a, c: a / np.sqrt(3.0)),
    "diamond111": dict(f=1.0 / 4.0, sigma_over_n=lambda a, c: a / (2.0 * np.sqrt(3.0))),
    "wurtzite0001": dict(f=1.0 / 4.0, sigma_over_n=lambda a, c: c / 4.0),
}


def t_eff(geometry, a, c=None):
    """Effective thickness t = f * sigma / n (m) of the terminating layer; a, c in m."""
    g = GEOMETRY[geometry]
    return float(g["f"] * g["sigma_over_n"](a, c))


def mode_heat_capacity(nu_hz, T):
    """Heat capacity of one mode (J/K); k_B at nu -> 0."""
    nu = np.asarray(nu_hz, dtype=float)
    x = H * np.clip(nu, 0.0, None) / (KB * T)
    out = np.full(x.shape, KB)
    big = x > 1e-6
    xb = np.minimum(x[big], 700.0)
    em1 = np.expm1(xb)
    out[big] = KB * xb * xb * np.exp(xb) / (em1 * em1)
    return out


def grid(nmax_thz, h=H_THZ):
    """Grid nodes nu_k = k h (THz), k = 0 .. n-1, covering [0, nmax_thz]."""
    n = int(np.ceil(nmax_thz / h)) + 2
    return np.arange(n) * h


def node_index(freq_thz, h, n):
    k = np.rint(np.asarray(freq_thz, dtype=float) / h).astype(int)
    return np.clip(k, 0, n - 1)


def dos(freq_thz, weights, h, n):
    """Per-atom DOS shape on the nodes (THz^-1), normalized so h * sum = 1. freq (gp, band),
    weights (gp,). Modes at or below 0 THz are dropped."""
    f = np.asarray(freq_thz, dtype=float)
    w = np.broadcast_to(np.asarray(weights, dtype=float)[:, None], f.shape)
    ok = f > 1e-6
    rho = np.bincount(node_index(f[ok], h, n), weights=w[ok], minlength=n)[:n].astype(float)
    s = rho.sum() * h
    return rho / s if s > 0 else rho


def _kernels(n, kernel):
    """Weight matrices on node indices (m, m'), frequencies in units of h (only ratios matter)."""
    k = np.arange(n, dtype=float)
    inv = np.zeros(n)
    inv[1:] = 1.0 / k[1:]
    if kernel == "jdos":
        return lambda a, b: np.ones(np.broadcast(a, b).shape), lambda a, c: np.ones(np.broadcast(a, c).shape)
    if kernel != "classical":
        raise ValueError(kernel)

    def wd(a, b):
        return (inv[a] + inv[b]) * inv[a] * inv[b]

    def wc(a, c):
        return (inv[a] - inv[c]) * inv[a] * inv[c]
    return wd, wc


def phase_space(rho_a, rho_b, h, kernel="classical"):
    """dec[a,b], dec_E[a,b] (energy share of the product drawn from a) and com[a,b] on the nodes."""
    ra = np.asarray(rho_a, dtype=float)
    rb = np.asarray(rho_b, dtype=float)
    n = ra.size
    wd, wc = _kernels(n, kernel)
    dec = np.zeros(n)
    dec_e = np.zeros(n)
    com = np.zeros(n)
    for k in range(1, n):
        m = np.arange(1, k)                     # products m and k - m, both > 0
        if m.size:
            t = ra[m] * rb[k - m] * wd(m, k - m)
            dec[k] = h * t.sum()
            dec_e[k] = h * (t * (m / k)).sum()
        mm = np.arange(1, n - k)                # partner m, product k + m
        if mm.size:
            com[k] = h * (ra[mm] * rb[k + mm] * wc(mm, k + mm)).sum()
    return dec, dec_e, com


def cross_fraction(rho_i, rho_j, h, kernel="classical"):
    """Q_ij on the nodes, the bulk phase space D_i, and the number of nodes where side i has
    states but D_i vanishes (Q set to 0 there) or Q hit the cap."""
    ri = np.asarray(rho_i, dtype=float)
    rj = np.asarray(rho_j, dtype=float)
    rx = 0.5 * (ri + rj)
    dii, _, cii = phase_space(ri, ri, h, kernel)
    D = dii + 2.0 * cii
    _, dje, _ = phase_space(rj, rx, h, kernel)
    _, _, cxj = phase_space(rx, rj, h, kernel)
    X = dje + cxj
    tiny = 1e-12 * max(float(np.max(D)), 1e-300)
    Q = np.where(D > tiny, X / np.where(D > tiny, D, 1.0), 0.0)
    n_zero = int(np.sum((ri > 0) & ~(D > tiny)))
    n_cap = int(np.sum(Q > Q_CAP))
    return np.minimum(Q, Q_CAP), D, n_zero, n_cap


def spectral_cg(freq_thz, weights, gamma_thz, T, vcell, h, n, classical=False):
    """Per node: volumetric heat capacity C (J m^-3 K^-1) and C * Gamma (W m^-3 K^-1), with
    Gamma = 4 pi gamma (phono3py's gamma in THz is the half linewidth: tau = 1 / (4 pi gamma))."""
    f = np.asarray(freq_thz, dtype=float)
    g = np.asarray(gamma_thz, dtype=float)
    w = np.broadcast_to(np.asarray(weights, dtype=float)[:, None], f.shape)
    ok = (f > 1e-6) & np.isfinite(g) & (g >= 0.0)
    c = np.full(f.shape, KB) if classical else mode_heat_capacity(f * 1e12, T)
    norm = float(np.sum(weights)) * vcell
    rate = 4.0 * np.pi * g * 1e12
    idx = node_index(f[ok], h, n)
    C = np.bincount(idx, weights=(w * c)[ok], minlength=n)[:n] / norm
    CG = np.bincount(idx, weights=(w * c * rate)[ok], minlength=n)[:n] / norm
    return C, CG


def g_inel_side(t, CG_i, Q_ij):
    return 0.5 * t * float(np.sum(np.asarray(CG_i) * np.asarray(Q_ij)))


def g_inel(side1, side2, T, h=H_THZ, kernel="classical", classical_c=False, q_cache=None):
    """G_inel (W m^-2 K^-1) between two sides. Each side is a dict with freq (gp, band) THz,
    weights (gp,), gamma (gp, band) THz at temperature T, vcell (m^3), t (m). Returns a dict with
    the total, the two one-sided parts, and the Q diagnostics."""
    nmax = max(float(np.max(side1["freq"])), float(np.max(side2["freq"])))
    n = grid(nmax, h).size
    if q_cache is not None and "Q12" in q_cache:
        Q12, Q21, diag = q_cache["Q12"], q_cache["Q21"], q_cache["diag"]
    else:
        r1 = dos(side1["freq"], side1["weights"], h, n)
        r2 = dos(side2["freq"], side2["weights"], h, n)
        Q12, _, z12, c12 = cross_fraction(r1, r2, h, kernel)
        Q21, _, z21, c21 = cross_fraction(r2, r1, h, kernel)
        diag = {"zero_nodes": [z12, z21], "capped_nodes": [c12, c21]}
        if q_cache is not None:
            q_cache.update(Q12=Q12, Q21=Q21, diag=diag)
    _, CG1 = spectral_cg(side1["freq"], side1["weights"], side1["gamma"], T, side1["vcell"], h, n,
                         classical_c)
    _, CG2 = spectral_cg(side2["freq"], side2["weights"], side2["gamma"], T, side2["vcell"], h, n,
                         classical_c)
    g1 = g_inel_side(side1["t"], CG1, Q12)
    g2 = g_inel_side(side2["t"], CG2, Q21)
    C1, _ = spectral_cg(side1["freq"], side1["weights"], side1["gamma"], T, side1["vcell"], h, n,
                        classical_c)
    C2, _ = spectral_cg(side2["freq"], side2["weights"], side2["gamma"], T, side2["vcell"], h, n,
                        classical_c)
    qbar1 = float(np.sum(CG1 * Q12) / max(np.sum(CG1), 1e-300))
    qbar2 = float(np.sum(CG2 * Q21) / max(np.sum(CG2), 1e-300))
    return {"G": g1 + g2, "G_side1": g1, "G_side2": g2, "Qbar": [qbar1, qbar2],
            "C": [float(C1.sum()), float(C2.sum())],
            "rate_mean": [float(CG1.sum() / max(C1.sum(), 1e-300)),
                          float(CG2.sum() / max(C2.sum(), 1e-300))], **diag}


def series(*gs):
    return 1.0 / sum(1.0 / g for g in gs)


def solve_scale(g_el, g_in, rise_meas, g_ep=None, lo=0.0, hi=1000.0):
    """Factor s such that the model's G(T_hi) - G(T_lo) equals the measured rise, where
    G(T) = series(g_ep, g_el(T) + s g_in(T)) (g_ep None: no series term). g_el, g_in: pairs
    (value at T_lo, value at T_hi). Returns s (0 if the elastic rise alone already reaches the
    measured rise, inf if even s = hi does not)."""
    def total(s, k):
        g = g_el[k] + s * g_in[k]
        return series(g_ep, g) if g_ep else g

    def rise(s):
        return total(s, 1) - total(s, 0)
    if rise(lo) >= rise_meas:
        return 0.0
    if rise(hi) < rise_meas:
        return float("inf")
    a, b = lo, hi
    for _ in range(200):
        m = 0.5 * (a + b)
        if rise(m) < rise_meas:
            a = m
        else:
            b = m
    return 0.5 * (a + b)
