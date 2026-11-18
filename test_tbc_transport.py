# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""Unit tests of tbc_transport (numpy only)."""
import math
import unittest

import numpy as np

import tbc_transport as tt

AL = {"rho": 2700.0, "vL": 6420.0, "vT": 3040.0, "n": 6.026e28}
SI = {"rho": 2329.0, "vL": 8433.0, "vT": 5843.0, "n": 4.994e28}


class TestTransport(unittest.TestCase):
    def test_heat_capacity_limits(self):
        c = tt.mode_heat_capacity(np.array([0.0, 1e6, 1e15]), 300.0)
        self.assertAlmostEqual(c[0] / tt.KB, 1.0, places=9)
        self.assertAlmostEqual(c[1] / tt.KB, 1.0, places=6)
        self.assertLess(c[2], 1e-30)
        # x = 1: x^2 e^x / (e^x - 1)^2 = 0.920673...
        nu = tt.KB * 300.0 / tt.H
        self.assertAlmostEqual(tt.mode_heat_capacity(np.array([nu]), 300.0)[0] / tt.KB,
                               math.e / (math.e - 1.0) ** 2, places=9)

    def test_debye_radiation_classical_limit(self):
        edges = tt.make_edges(30.0, 0.0015)
        M = tt.debye_flux(AL, edges)
        g = tt.g_radiation(M, edges, 1e6)  # classical limit
        self.assertAlmostEqual(g / tt.g_rad_debye_classical(AL), 1.0, delta=2e-3)

    def test_flux_spectrum_isotropic_vs_oriented(self):
        # a single branch with velocity along +z for half the q points and -z for the other half
        nq = 1000
        rng = np.random.default_rng(1)
        f = rng.uniform(1.0, 2.0, size=(nq, 1))
        v = np.zeros((nq, 1, 3))
        v[: nq // 2, 0, 2] = 1000.0
        v[nq // 2:, 0, 2] = -1000.0
        edges = tt.make_edges(3.0, 0.01)
        w = np.ones(nq)
        Mz = tt.flux_spectrum(f, v, w, 1e-29, edges, normal=[0, 0, 1])
        Miso = tt.flux_spectrum(f, v, w, 1e-29, edges, normal=None)
        _, wid = tt.centers_widths(edges)
        # oriented: half the modes at 1000 m/s -> 500 m/s / V ; isotropic: 1000/4 = 250 m/s / V
        self.assertAlmostEqual(np.sum(Mz * wid) * 1e-29, 500.0, delta=1e-6)
        self.assertAlmostEqual(np.sum(Miso * wid) * 1e-29, 250.0, delta=1e-6)
        # negative frequencies (imaginary modes) are dropped
        Mneg = tt.flux_spectrum(-f, v, w, 1e-29, edges, normal=None)
        self.assertEqual(float(np.sum(Mneg)), 0.0)

    def test_dmm_identical_and_detailed_balance(self):
        edges = tt.make_edges(20.0, 0.01)
        M1 = tt.debye_flux(AL, edges)
        M2 = tt.debye_flux(SI, edges)
        g11 = tt.g_dmm(M1, M1, edges, 300.0)
        self.assertAlmostEqual(g11 / tt.g_radiation(M1, edges, 300.0), 0.5, places=9)
        for T in (100.0, 300.0, 700.0):
            self.assertAlmostEqual(tt.g_dmm(M1, M2, edges, T) / tt.g_dmm(M2, M1, edges, T), 1.0,
                                   places=9)

    def test_simons(self):
        # identical, diffuse: tau = 1/2 -> corrected G equals the radiation limit
        self.assertAlmostEqual(tt.simons_correct(0.5, 1.0, 1.0), 1.0, places=12)
        # identical, perfect transmission -> infinite
        self.assertEqual(tt.simons_correct(1.0, 1.0, 1.0), float("inf"))
        # gray two-sided check: tau12 = 0.2 with grad1 = 1, grad2 = 0.5 -> tau21 = 0.4
        g = 0.2
        exp = 0.2 / (1.0 - (0.2 + 0.4) / 2.0)
        self.assertAlmostEqual(tt.simons_correct(g, 1.0, 0.5), exp, places=12)

    def test_amm(self):
        # identical media: tau = 1 for all angles -> Gamma = 1
        self.assertAlmostEqual(tt.amm_gamma(1.0, 1.0, 5000.0, 5000.0), 1.0, delta=1e-6)
        # normal-incidence-dominated check: Gamma <= 4 Z1 Z2 / (Z1 + Z2)^2 when v1 = v2
        Z1, Z2 = 1.0, 3.0
        self.assertAlmostEqual(tt.amm_gamma(Z1, Z2, 5000.0, 5000.0), 4 * Z1 * Z2 / (Z1 + Z2) ** 2,
                               delta=1e-6)
        edges = tt.make_edges(20.0, 0.01)
        g_same = tt.g_amm_debye(AL, AL, edges, 300.0)
        self.assertAlmostEqual(g_same / tt.g_radiation(tt.debye_flux(AL, edges), edges, 300.0), 1.0,
                               delta=1e-5)
        # AMM conductance between different media is below side 1's radiation limit
        self.assertLess(tt.g_amm_debye(AL, SI, edges, 300.0),
                        tt.g_radiation(tt.debye_flux(AL, edges), edges, 300.0))

    def test_series_and_ep(self):
        self.assertAlmostEqual(tt.series(2.0, 2.0), 1.0, places=12)
        self.assertAlmostEqual(tt.g_ep(1e17, 4.0), math.sqrt(4e17), places=3)

    def test_al_si_debye_magnitude(self):
        # sanity: Debye-model Al/Si DMM at 300 K lies within 100-2000 MW m^-2 K^-1
        edges = tt.make_edges(20.0, 0.01)
        g = tt.g_dmm_debye(AL, SI, edges, 300.0)
        self.assertGreater(g, 1e8)
        self.assertLess(g, 2e9)


if __name__ == "__main__":
    unittest.main()
