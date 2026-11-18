# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""Unit tests of tbc_inel (numpy only)."""
import unittest

import numpy as np

import tbc_inel as ti

G_EP = float(np.sqrt(2.45e17 * 6.0))   # the Al electron-phonon term used by run_tbc_inel


def debye_side(numax, ngp=400, nb=3, vcell=2e-29, t=5e-11, gamma_per_thz2=0.002, seed=0):
    """A synthetic side: Debye-like spectrum up to numax (THz), gamma ~ nu^2 (THz)."""
    rng = np.random.default_rng(seed)
    u = rng.uniform(0.0, 1.0, (ngp, nb))
    freq = numax * np.cbrt(u)                       # density ~ nu^2
    return {"freq": freq, "weights": np.ones(ngp), "vcell": vcell, "t": t,
            "gamma": gamma_per_thz2 * freq ** 2}


class TestInel(unittest.TestCase):
    def test_identical_sides_half(self):
        s = debye_side(10.0)
        n = ti.grid(10.0).size
        rho = ti.dos(s["freq"], s["weights"], ti.H_THZ, n)
        for kernel in ("classical", "jdos"):
            Q, D, nz, nc = ti.cross_fraction(rho, rho, ti.H_THZ, kernel)
            ok = D > 1e-9 * D.max()
            np.testing.assert_allclose(Q[ok], 0.5, rtol=1e-12)
            self.assertEqual(nc, 0)

    def test_dos_normalized(self):
        s = debye_side(8.0)
        n = ti.grid(8.0).size
        rho = ti.dos(s["freq"], s["weights"], ti.H_THZ, n)
        self.assertAlmostEqual(rho.sum() * ti.H_THZ, 1.0, places=12)

    def test_classical_linear_in_T(self):
        # classical heat capacity and rates proportional to T: G_inel proportional to T
        a = debye_side(9.0, seed=1)
        b = debye_side(15.0, seed=2, t=3e-11)
        out = []
        for T in (300.0, 600.0):
            a["gamma"] = 0.002 * a["freq"] ** 2 * T / 300.0
            b["gamma"] = 0.001 * b["freq"] ** 2 * T / 300.0
            out.append(ti.g_inel(a, b, T, classical_c=True)["G"])
        self.assertAlmostEqual(out[1] / out[0], 2.0, places=10)

    def test_zero_rates_zero_G(self):
        a = debye_side(9.0, seed=1)
        b = debye_side(15.0, seed=2)
        a["gamma"] = np.zeros_like(a["freq"])
        b["gamma"] = np.zeros_like(b["freq"])
        self.assertEqual(ti.g_inel(a, b, 300.0)["G"], 0.0)

    def test_swap_symmetry(self):
        a = debye_side(9.0, seed=1)
        b = debye_side(15.0, seed=2, t=3e-11)
        g1 = ti.g_inel(a, b, 300.0)["G"]
        g2 = ti.g_inel(b, a, 300.0)["G"]
        self.assertAlmostEqual(g1 / g2, 1.0, places=12)

    def test_mismatch_suppresses_cross(self):
        # a soft side next to a very stiff one: the soft side's energy rarely lands in the stiff
        # side (few stiff-side states at low frequency), so its Q is far below 1/2
        a = debye_side(2.0, seed=3)
        b = debye_side(40.0, seed=4)
        n = ti.grid(40.0).size
        ra = ti.dos(a["freq"], a["weights"], ti.H_THZ, n)
        rb = ti.dos(b["freq"], b["weights"], ti.H_THZ, n)
        Q, D, _, _ = ti.cross_fraction(ra, rb, ti.H_THZ)
        cg = ti.spectral_cg(a["freq"], a["weights"], a["gamma"], 300.0, a["vcell"], ti.H_THZ, n)[1]
        qbar = float(np.sum(cg * Q) / np.sum(cg))
        self.assertLess(qbar, 0.1)

    def test_quantum_below_classical(self):
        a = debye_side(9.0, seed=1)
        b = debye_side(15.0, seed=2)
        gq = ti.g_inel(a, b, 100.0)["G"]
        gc = ti.g_inel(a, b, 100.0, classical_c=True)["G"]
        self.assertLess(gq, gc)

    def test_t_eff(self):
        a = 4.04e-10
        self.assertAlmostEqual(ti.t_eff("fcc111", a), 0.25 * a / np.sqrt(3.0))
        self.assertAlmostEqual(ti.t_eff("diamond111", 5.43e-10), 0.25 * 5.43e-10 / (2 * np.sqrt(3.0)))
        self.assertAlmostEqual(ti.t_eff("wurtzite0001", 3.19e-10, 5.19e-10), 5.19e-10 / 16.0)

    def test_heat_capacity_limits(self):
        self.assertAlmostEqual(ti.mode_heat_capacity(np.array([1e9]), 300.0)[0] / ti.KB, 1.0, places=4)
        self.assertLess(ti.mode_heat_capacity(np.array([40e12]), 100.0)[0] / ti.KB, 1e-3)

    def test_single_frequency_no_phase_space(self):
        # every mode at one frequency: no nu' + nu'' = nu exists, so no channel and G = 0
        f = np.full((50, 3), 5.0)
        side = {"freq": f, "weights": np.ones(50), "vcell": 2e-29, "t": 1e-10, "gamma": np.full(f.shape, 0.1)}
        r = ti.g_inel(side, side, 300.0)
        self.assertEqual(r["G"], 0.0)
        self.assertEqual(r["zero_nodes"], [1, 1])

    def test_magnitude_identical_sides(self):
        # identical sides, Q = 1/2 wherever there is phase space: G = 2 * 1/2 * t * sum(C Gamma Q)
        s = debye_side(10.0, seed=5)
        n = ti.grid(10.0).size
        r = ti.g_inel(s, s, 300.0)
        rho = ti.dos(s["freq"], s["weights"], ti.H_THZ, n)
        Q, _, _, _ = ti.cross_fraction(rho, rho, ti.H_THZ)
        cg = ti.spectral_cg(s["freq"], s["weights"], s["gamma"], 300.0, s["vcell"], ti.H_THZ, n)[1]
        self.assertAlmostEqual(r["G"] / (s["t"] * float(np.sum(cg * Q))), 1.0, places=12)
        self.assertLessEqual(r["Qbar"][0], 0.5 + 1e-12)

    def test_solve_scale(self):
        g_el = (300e6, 320e6)
        g_in = (20e6, 40e6)
        s = ti.solve_scale(g_el, g_in, rise_meas=60e6)
        self.assertAlmostEqual(s, 2.0, places=6)
        self.assertEqual(ti.solve_scale(g_el, g_in, rise_meas=10e6), 0.0)
        s2 = ti.solve_scale(g_el, g_in, rise_meas=30e6, g_ep=1.2e9)
        tot = lambda k: ti.series(1.2e9, g_el[k] + s2 * g_in[k])
        self.assertAlmostEqual((tot(1) - tot(0)) / 30e6, 1.0, places=6)

    def test_solve_scale_series_peak(self):
        # With a series electron-phonon term the rise peaks and falls back towards 0 at large s.
        # Here the rise at s = 1000 (~36 MW) is below the target (60 MW) although s ~ 10-20
        # matches it: 0.1.0 returned inf for this case.
        g_el = (400e6, 400e6)
        g_in = (15e6, 26.25e6)
        self.assertLess(ti.model_rise(g_el, g_in, 1000.0, G_EP), 60e6)
        s = ti.solve_scale(g_el, g_in, rise_meas=60e6, g_ep=G_EP)
        self.assertTrue(10.0 < s < 20.0, s)
        self.assertAlmostEqual(ti.model_rise(g_el, g_in, s, G_EP) / 60e6, 1.0, places=6)
        self.assertLess(ti.model_rise(g_el, g_in, 0.99 * s, G_EP), 60e6)   # smallest root

    def test_solve_scale_unreachable(self):
        g_el = (400e6, 400e6)
        g_in = (15e6, 26.25e6)
        s_max, r_max = ti.max_rise(g_el, g_in, g_ep=G_EP)
        self.assertTrue(100e6 < r_max < 140e6, r_max)
        self.assertTrue(1.0 < s_max < 1000.0, s_max)
        self.assertEqual(ti.solve_scale(g_el, g_in, rise_meas=200e6, g_ep=G_EP), float("inf"))
        s = ti.solve_scale(g_el, g_in, rise_meas=0.99 * r_max, g_ep=G_EP)
        self.assertTrue(np.isfinite(s) and s <= s_max * 1.01, s)

    def test_ep_pair(self):
        g_el = (300e6, 320e6)
        g_in = (20e6, 40e6)
        a = ti.solve_scale(g_el, g_in, rise_meas=30e6, g_ep=G_EP)
        b = ti.solve_scale(g_el, g_in, rise_meas=30e6, g_ep=(G_EP, G_EP))
        self.assertAlmostEqual(a, b, places=10)
        # a falling series conductance at T_hi lowers the rise for the same s
        self.assertLess(ti.model_rise(g_el, g_in, 2.0, (G_EP, 0.8 * G_EP)),
                        ti.model_rise(g_el, g_in, 2.0, G_EP))


if __name__ == "__main__":
    unittest.main()
