# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""tbc_transport: elastic phonon models of thermal boundary conductance G(T) across a solid-solid
interface, written as baselines for the Discovered Materials Grand Challenge 05 ("Predicting Heat
Transport across Solid Interfaces").

Only numpy. Everything is in SI units: frequencies nu in Hz (ordinary, not angular), spectral fluxes
M(nu) in 1/(m^2 s Hz), conductances in W m^-2 K^-1.

Definitions
- Spectral phonon flux of side i toward the interface (per unit area, per unit frequency):
      M_i(nu) = (1 / (N_q V_cell)) * sum_{q,s} max(v_qs . n, 0) * delta(nu - nu_qs)
  For random orientation (polycrystal) max(v.n, 0) averages to |v| / 4 over directions.
- Phonon radiation limit of side i:  G_rad,i(T) = int c(nu, T) M_i(nu) dnu,
  with c(nu, T) = h nu df/dT the heat capacity of one mode (Bose-Einstein f).
- Diffuse mismatch model (DMM):  tau_12(nu) = M_2 / (M_1 + M_2),  G = int c M_1 tau_12 dnu.
  (detailed balance M_1 tau_12 = M_2 tau_21 holds by construction.)
- Gray non-equilibrium (Simons 1974; Zeng and Chen 2001) correction, which refers G to the local
  equilibrium temperatures on both sides instead of the emitted-phonon temperatures:
      1/G_corr = 1/G - (1/2) (1/G_rad,1 + 1/G_rad,2).
- Debye models: branch j of side i with sound speed v_ij carries pi nu^2 / v_ij^2 of isotropic flux up
  to its own cutoff nu_c,ij = v_ij (3 n_i / (4 pi))^(1/3) (n_i = atoms per m^3).
- Acoustic mismatch model (AMM, Little 1959) per branch, no mode conversion:
      tau_j(theta) = 4 Z1 Z2 cos(theta1) cos(theta2) / (Z1 cos(theta1) + Z2 cos(theta2))^2,
      Z = rho v, Snell sin(theta2) = (v2/v1) sin(theta1), tau = 0 beyond the critical angle,
      Gamma_j = 2 int_0^{pi/2} tau_j cos sin dtheta,  G = sum_j Gamma_j int c M_1j dnu,
  where M_1j is cut at min(nu_c,1j, nu_c,2j) (elastic: the mode must exist on side 2).
- Electron-phonon resistance inside a metal film in series with the phonon conductance
  (Majumdar and Reddy 2004): G_ep = sqrt(g Lambda_p), g the electron-phonon coupling factor and
  Lambda_p the phonon thermal conductivity of the metal.
"""
import numpy as np

VERSION = "mv2 tbc-transport 0.1.0"

H = 6.62607015e-34  # J s
KB = 1.380649e-23   # J / K


def mode_heat_capacity(nu_hz, T):
    """h nu df/dT of one phonon mode (J/K) for Bose-Einstein occupation; -> k_B as nu -> 0."""
    nu = np.asarray(nu_hz, dtype=float)
    x = H * nu / (KB * float(T))
    out = np.full(x.shape, KB, dtype=float)
    m = x > 1e-8
    xm = np.minimum(x[m], 700.0)
    em = np.exp(-xm)
    out[m] = KB * xm * xm * em / (1.0 - em) ** 2
    return out


def make_edges(nu_max_thz=48.0, dnu_thz=0.02):
    """Common frequency bin edges in Hz."""
    n = int(round(nu_max_thz / dnu_thz))
    return np.linspace(0.0, nu_max_thz * 1e12, n + 1)


def flux_spectrum(freq_thz, vel_ms, weights, vcell_m3, edges_hz, normal=None):
    """Spectral flux M(nu) (1/(m^2 s Hz)) on the given bins from mesh frequencies (THz, shape
    (nq, nb)), group velocities (m/s, shape (nq, nb, 3)) and q-point weights (nq,).
    normal=None: isotropic average |v|/4; else the half-space flux along the unit normal."""
    f = np.asarray(freq_thz, dtype=float) * 1e12
    v = np.asarray(vel_ms, dtype=float)
    w = np.asarray(weights, dtype=float)
    if normal is None:
        vn = np.linalg.norm(v, axis=-1) / 4.0
    else:
        nv = np.asarray(normal, dtype=float)
        nv = nv / np.linalg.norm(nv)
        vn = np.clip(v @ nv, 0.0, None)
    wt = (w / w.sum())[:, None] / float(vcell_m3) * vn
    mask = f > 0.0
    hist = np.histogram(f[mask], bins=edges_hz, weights=wt[mask])[0]
    return hist / np.diff(edges_hz)


def centers_widths(edges_hz):
    e = np.asarray(edges_hz, dtype=float)
    return 0.5 * (e[1:] + e[:-1]), np.diff(e)


def g_radiation(M, edges_hz, T):
    c, w = centers_widths(edges_hz)
    return float(np.sum(mode_heat_capacity(c, T) * M * w))


def dmm_transmission(M1, M2):
    tot = M1 + M2
    return np.where(tot > 0.0, M2 / np.where(tot > 0.0, tot, 1.0), 0.0)


def g_dmm(M1, M2, edges_hz, T):
    """Conventional (emitted-temperature) DMM conductance, side 1 -> side 2."""
    c, w = centers_widths(edges_hz)
    tau = dmm_transmission(M1, M2)
    return float(np.sum(mode_heat_capacity(c, T) * M1 * tau * w))


def simons_correct(G, grad1, grad2):
    """Gray non-equilibrium correction; returns inf if the correction exceeds 1/G."""
    inv = 1.0 / G - 0.5 * (1.0 / grad1 + 1.0 / grad2)
    return float(1.0 / inv) if inv > 0.0 else float("inf")


def series(*gs):
    """Conductances in series (W m^-2 K^-1)."""
    return float(1.0 / sum(1.0 / g for g in gs))


def g_ep(g_coupling, lambda_p):
    """Electron-phonon conductance inside a metal, sqrt(g * Lambda_p) (W m^-2 K^-1)."""
    return float(np.sqrt(g_coupling * lambda_p))


def debye_cutoff(v, n_atoms):
    return v * (3.0 * n_atoms / (4.0 * np.pi)) ** (1.0 / 3.0)


def debye_branches(m):
    """(speed, cutoff) of the L, T, T branches of a material dict with rho, vL, vT, n."""
    out = []
    for v in (m["vL"], m["vT"], m["vT"]):
        out.append((v, debye_cutoff(v, m["n"])))
    return out


def debye_flux(m, edges_hz, cut_with=None):
    """Isotropic Debye spectral flux; if cut_with is another material, each branch is cut at the
    lower of the two materials' cutoffs for that branch."""
    c, _ = centers_widths(edges_hz)
    M = np.zeros_like(c)
    b1 = debye_branches(m)
    b2 = debye_branches(cut_with) if cut_with is not None else b1
    for (v, nc), (_, nc2) in zip(b1, b2):
        cut = min(nc, nc2)
        M += np.where(c < cut, np.pi * c * c / (v * v), 0.0)
    return M


def amm_gamma(Z1, Z2, v1, v2, nth=4000):
    """Angle-averaged AMM transmission 2 int tau cos sin dtheta for one branch."""
    th = (np.arange(nth) + 0.5) * (0.5 * np.pi / nth)
    s2 = (v2 / v1) * np.sin(th)
    ok = s2 < 1.0
    c1 = np.cos(th)
    c2 = np.sqrt(np.clip(1.0 - s2 * s2, 0.0, None))
    tau = np.where(ok, 4.0 * Z1 * Z2 * c1 * c2 / (Z1 * c1 + Z2 * c2) ** 2, 0.0)
    return float(2.0 * np.sum(tau * c1 * np.sin(th)) * (0.5 * np.pi / nth))


def g_amm_debye(m1, m2, edges_hz, T):
    c, w = centers_widths(edges_hz)
    cap = mode_heat_capacity(c, T)
    G = 0.0
    for (v1, nc1), (v2, nc2) in zip(debye_branches(m1), debye_branches(m2)):
        gam = amm_gamma(m1["rho"] * v1, m2["rho"] * v2, v1, v2)
        Mj = np.where(c < min(nc1, nc2), np.pi * c * c / (v1 * v1), 0.0)
        G += gam * float(np.sum(cap * Mj * w))
    return G


def g_dmm_debye(m1, m2, edges_hz, T):
    return g_dmm(debye_flux(m1, edges_hz), debye_flux(m2, edges_hz), edges_hz, T)


def g_rad_debye_classical(m):
    """High-temperature limit of the Debye radiation limit: k_B sum_j pi nu_c,j^3 / (3 v_j^2)."""
    return float(sum(KB * np.pi * nc ** 3 / (3.0 * v * v) for v, nc in debye_branches(m)))
