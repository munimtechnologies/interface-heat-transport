# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""tbc_bte: thermal boundary conductance from a multiband two-stream Boltzmann transport model
with mode-resolved bulk lifetimes on both sides of an elastic (DMM) interface.

Only numpy. SI units. Stage 2a of the interface-heat-transport study (PLAN-amendment-0.2.md).

Model
-----
Frequencies are cut into bins (bands) b. On side i (1: x < 0, 2: x > 0) band b has
- K_ib: its phonon radiation conductance toward the interface, int_bin c(nu, T) M_i(nu) dnu
  (W m^-2 K^-1), with M_i the spectral flux of tbc_transport.flux_spectrum;
- l_ib: its two-stream length, l_ib = kappa_ib / (2 K_ib), where kappa_ib is the band's share of
  the bulk thermal conductivity along the interface normal (W m^-1 K^-1). With this choice the
  two-stream bulk conductivity sum_b 2 K_ib l_ib equals the input kappa_i exactly. For an isotropic
  band l = 2 Lambda / 3 (Lambda = v tau).
Each band has a stream moving toward the interface ("in") and one moving away ("out"), each
described by a temperature. A stream relaxes toward the local bath temperature theta(x) over its
length l: dT/ds = -(T - theta)/l along its direction of travel. theta(x) follows from energy
conservation: the total heat flux sum_b K_b (T_in - T_out) is the same at every x.
At x = 0 the interface is elastic and diffuse (DMM): of band b's incident flux on side 1 a
fraction alpha_1b = K_2b / (K_1b + K_2b) is transmitted into band b of side 2, the rest reflected
into band b of side 1 (and symmetrically). Detailed balance K_1b alpha_1b = K_2b alpha_2b holds.
Far from the interface both sides carry the bulk linear profile. G is the heat flux divided by the
jump between the two bulk profiles extrapolated to x = 0 (what a thermoreflectance fit of a
two-layer model with bulk conductivities would call the interface conductance).

The problem is linear. We write T = theta_lin(x) + (bulk stream offset) + u and solve for the
deviations u and phi = theta - theta_lin, so the far-field profiles never enter numerically (no
cancellation). theta is piecewise linear between grid nodes, and each stream is integrated
exactly across each cell for that theta. Grid: x = 0 and a geometric grid from x_min to x_max on
each side.

Checks (test_tbc_bte.py): one band (gray) gives exactly the Simons/Zeng-Chen corrected DMM,
1/G = 1/G_dmm - 1/(2 K_1) - 1/(2 K_2); swapping the two sides gives the same G; multiplying every
length l by the same factor leaves G unchanged.

The last property is the point of stage 2a. In the classical limit (all heat capacities k_B)
three-phonon scattering rates all grow in proportion to T, so every length scales by the same
factor and the model's G cannot change with temperature. More generally, for semi-infinite
solids and a linearized Boltzmann equation whose collision operator is multiplied by a constant,
x -> x / lambda maps the solution onto itself, so G is unchanged; this holds beyond the two-stream
and relaxation-time approximations. Bulk anharmonic thermalization near the interface can
therefore only make G depend on T through (a) quantum heat capacities, (b) scattering that does
not grow with T (isotopes, defects, boundaries, finite films), (c) four-phonon and higher
processes, or (d) temperature-dependent harmonic properties.
"""
import numpy as np

VERSION = "mv2 tbc-bte 0.1.0"

L_MIN = 1e-12   # m, floor on two-stream lengths
L_MAX = 1e-2    # m, cap on two-stream lengths (longer ones are treated as 1 cm)


def dmm_alpha(K1, K2):
    """Per-band transmitted fractions (alpha_1, alpha_2) and band conductances Gamma."""
    K1 = np.asarray(K1, dtype=float)
    K2 = np.asarray(K2, dtype=float)
    tot = K1 + K2
    ok = (K1 > 0.0) & (K2 > 0.0)
    safe = np.where(tot > 0.0, tot, 1.0)
    a1 = np.where(ok, K2 / safe, 0.0)
    a2 = np.where(ok, K1 / safe, 0.0)
    gam = np.where(ok, K1 * K2 / safe, 0.0)
    return a1, a2, gam


def g_dmm_bands(K1, K2):
    return float(np.sum(dmm_alpha(K1, K2)[2]))


def g_simons_gray(K1, K2):
    """Gray Simons/Zeng-Chen correction of the band-summed DMM; inf if not positive."""
    g = g_dmm_bands(K1, K2)
    inv = 1.0 / g - 0.5 / float(np.sum(K1)) - 0.5 / float(np.sum(K2))
    return float(1.0 / inv) if inv > 0.0 else float("inf")


def g_linear_bath(K1, K2, l1, l2):
    """Closed form if the bath temperature stayed linear up to the interface:
    G = G_dmm / (1 - sum_b Gamma_b (l_1b / kappa_1 + l_2b / kappa_2)); inf if not positive."""
    K1, K2 = np.asarray(K1, float), np.asarray(K2, float)
    l1, l2 = np.asarray(l1, float), np.asarray(l2, float)
    _, _, gam = dmm_alpha(K1, K2)
    k1 = float(np.sum(2.0 * K1 * l1))
    k2 = float(np.sum(2.0 * K2 * l2))
    s = float(np.sum(gam * (l1 / k1 + l2 / k2)))
    g = float(np.sum(gam))
    return g / (1.0 - s) if s < 1.0 else float("inf")


def _step(e, p, coef, iA, iB):
    """Propagate stream coefficient rows (B, n) from node A to node B (direction of travel):
    u_B = (1 - p) phi_B + (p - e) phi_A + e u_A."""
    out = coef * e[:, None]
    out[:, iB] += 1.0 - p
    out[:, iA] += p - e
    return out


def _cell(delta, ell):
    a = delta / ell
    e = np.exp(-a)
    p = np.where(a > 1e-12, -np.expm1(-a) / np.where(a > 1e-12, a, 1.0), 1.0 - 0.5 * a)
    return e, p


def solve_two_stream(K1, K2, l1, l2, n_nodes=160, x_min=1e-11, x_max=0.1, return_profile=False):
    """Interface conductance G (W m^-2 K^-1) of the multiband two-stream model.

    K1, K2: (B,) band radiation conductances of sides 1 and 2 on common frequency bins.
    l1, l2: (B,) two-stream lengths (m); only bands with K > 0 on a side are used there.
    Returns a dict with G, the conventional DMM G_dmm, the gray Simons value, the linear-bath
    value, the bulk conductivities, the residual of the linear solve and the share of each side's
    conductivity carried by lengths capped at L_MAX."""
    K1 = np.asarray(K1, dtype=float).copy()
    K2 = np.asarray(K2, dtype=float).copy()
    l1 = np.clip(np.asarray(l1, dtype=float), L_MIN, None)
    l2 = np.clip(np.asarray(l2, dtype=float), L_MIN, None)
    cap1 = float(np.sum(2 * K1 * l1 * (l1 > L_MAX)) / max(np.sum(2 * K1 * l1), 1e-300))
    cap2 = float(np.sum(2 * K2 * l2 * (l2 > L_MAX)) / max(np.sum(2 * K2 * l2), 1e-300))
    l1 = np.minimum(l1, L_MAX)
    l2 = np.minimum(l2, L_MAX)
    m1 = K1 > 0.0
    m2 = K2 > 0.0
    a1, a2, gam = dmm_alpha(K1, K2)
    kap1 = float(np.sum(2.0 * K1 * l1))
    kap2 = float(np.sum(2.0 * K2 * l2))
    Q = 1.0
    off1 = np.where(m1, l1 * Q / kap1, 0.0)   # bulk stream offsets a_1b, a_2b
    off2 = np.where(m2, l2 * Q / kap2, 0.0)

    X = np.concatenate([[0.0], np.geomspace(x_min, x_max, n_nodes - 1)])
    N = n_nodes
    nun = 2 * N + 1          # phi on side 1 (0..N-1), side 2 (N..2N-1), constant column 2N
    C = 2 * N
    i1 = np.arange(N)
    i2 = N + np.arange(N)

    B1 = int(np.sum(m1))
    B2 = int(np.sum(m2))
    K1a, l1a, a1a, off1a = K1[m1], l1[m1], a1[m1], off1[m1]
    K2a, l2a, a2a, off2a = K2[m2], l2[m2], a2[m2], off2[m2]

    flux1 = np.zeros((N, nun))   # sum_b K_1b (u_in - u_out) at side-1 nodes
    flux2 = np.zeros((N, nun))   # sum_b K_2b (u_out - u_in) at side-2 nodes

    def inward(nb, ell, idx, Kw, flux, sign):
        coef = np.zeros((nb, nun))
        coef[:, idx[N - 1]] = 1.0          # far boundary: incoming stream at bath temperature
        flux[N - 1] += sign * (Kw @ coef)
        for k in range(N - 2, -1, -1):
            e, p = _cell(X[k + 1] - X[k], ell)
            coef = _step(e, p, coef, idx[k + 1], idx[k])
            flux[k] += sign * (Kw @ coef)
        return coef

    def outward(coef, ell, idx, Kw, flux, sign):
        flux[0] += sign * (Kw @ coef)
        for k in range(0, N - 1):
            e, p = _cell(X[k + 1] - X[k], ell)
            coef = _step(e, p, coef, idx[k], idx[k + 1])
            flux[k + 1] += sign * (Kw @ coef)

    in1 = inward(B1, l1a, i1, K1a, flux1, +1.0)
    in2 = inward(B2, l2a, i2, K2a, flux2, -1.0)
    # full-length in-stream coefficient rows (zeros for inactive bands) to couple the two sides
    full1 = np.zeros((K1.size, nun))
    full1[m1] = in1
    full2 = np.zeros((K2.size, nun))
    full2[m2] = in2
    out1 = (1.0 - a1a)[:, None] * in1 + a1a[:, None] * full2[m1]
    out1[:, C] += (2.0 - a1a) * off1a - a1a * off2[m1]
    out2 = (1.0 - a2a)[:, None] * in2 + a2a[:, None] * full1[m2]
    out2[:, C] += -(2.0 - a2a) * off2a + a2a * off1[m2]
    outward(out1, l1a, i1, K1a, flux1, -1.0)
    outward(out2, l2a, i2, K2a, flux2, +1.0)

    pin = np.zeros((1, nun))
    pin[0, i1[N - 1]] = 1.0
    A = np.vstack([flux1, flux2, pin])
    # scale rows so that flux rows and the pin row are comparable
    scale = np.max(np.abs(A[:, :C]), axis=1)
    scale[scale == 0.0] = 1.0
    A = A / scale[:, None]
    M = A[:, :C]
    rhs = -A[:, C]
    phi, *_ = np.linalg.lstsq(M, rhs, rcond=None)
    resid = float(np.max(np.abs(M @ phi - rhs)))
    dT = float(phi[i1[N - 1]] - phi[i2[N - 1]])
    out = {
        "G": Q / dT if dT > 0.0 else float("inf"),
        "G_dmm": float(np.sum(gam)),
        "G_simons": g_simons_gray(K1, K2),
        "G_linear_bath": g_linear_bath(K1, K2, l1, l2),
        "kappa1": kap1, "kappa2": kap2, "residual": resid,
        "capped_share1": cap1, "capped_share2": cap2,
    }
    if return_profile:
        out["x"] = X
        out["phi1"] = phi[i1]
        out["phi2"] = phi[i2]
    return out


def lengths_from_kappa(kappa_b, K_b):
    """Two-stream lengths l_b = kappa_b / (2 K_b) on bins with K_b > 0. Bins with K_b > 0 but no
    kappa (a coarser lifetime mesh left them empty) get log-interpolated lengths from their
    neighbours. Returns (l_b, number of filled bins)."""
    kb = np.asarray(kappa_b, dtype=float)
    Kb = np.asarray(K_b, dtype=float)
    ell = np.zeros_like(Kb)
    have = (Kb > 0.0) & (kb > 0.0)
    ell[have] = kb[have] / (2.0 * Kb[have])
    need = (Kb > 0.0) & ~have
    if np.any(need) and np.any(have):
        idx = np.arange(Kb.size)
        ell[need] = np.exp(np.interp(idx[need], idx[have], np.log(ell[have])))
    return ell, int(np.sum(need))


def rebin(values_fine, edges_fine, edges_coarse):
    """Sum a per-bin quantity from fine bins into coarse bins (coarse edges on the fine grid)."""
    v = np.asarray(values_fine, dtype=float)
    cf = 0.5 * (np.asarray(edges_fine)[1:] + np.asarray(edges_fine)[:-1])
    j = np.searchsorted(edges_coarse, cf, side="right") - 1
    ok = (j >= 0) & (j < len(edges_coarse) - 1)
    return np.bincount(j[ok], weights=v[ok], minlength=len(edges_coarse) - 1)
