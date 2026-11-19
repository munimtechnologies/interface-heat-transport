# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Made with AI agents working in MV2's name. Expected numbers below were worked out by hand
# (PLAN-amendment-0.3.3.md shows the arithmetic).
import math
import unittest

import tbc_direct as td
from tbc_gep import series

AL_SI = {"300": 379e6, "400": 400e6, "500": 446e6, "600": 478e6, "700": 530e6}
AL_GAN = {"300": 422e6, "400": 432e6, "500": 448e6, "600": 479e6, "700": 546e6}
MW = 1e6
# Al inputs of De San Feliciano et al. (arXiv:2510.27499): g, k_p; k_e = 237 - 8.95 ~ 228
G_AL, KE_AL, KP_AL = 2.46e17, 228.0, 8.95
# Al inputs of Lombard et al. Table 1: gamma = 1.3 mJ/mol/K^2 x 1e5 mol/m^3, lambda<w^2>, u_L, u_T, n_e
GAM_AL, LW2_AL, UL_AL, UT_AL, NE_AL = 130.0, 185.9, 6240.0, 3040.0, 1.81e29
S_DSF = 107e6 / 300.0  # a T-linear channel equal to De San Feliciano's 107 MW at 300 K


class TestCircuit(unittest.TestCase):
    def test_no_direct_channel_is_majumdar_reddy(self):
        fe = KE_AL / (KE_AL + KP_AL)
        gmr = math.sqrt(G_AL * KP_AL) / fe ** 1.5
        got = td.circuit_exact(G_AL, KE_AL, KP_AL, 0.0, 210e6)
        self.assertAlmostEqual(got / series(gmr, 210e6), 1.0, places=10)
        self.assertAlmostEqual(got / MW, 185.25, delta=0.1)

    def test_with_direct_channel(self):
        self.assertAlmostEqual(td.circuit_exact(G_AL, KE_AL, KP_AL, 107e6, 210e6) / MW, 293.24, delta=0.3)
        self.assertAlmostEqual(td.circuit_limit(math.sqrt(G_AL * KP_AL), 107e6, 210e6) / MW, 290.96, delta=0.1)

    def test_only_direct_channel(self):
        self.assertAlmostEqual(td.circuit_exact(G_AL, KE_AL, KP_AL, 107e6, 0.0) / MW, 106.989, delta=0.02)
        self.assertEqual(td.circuit_limit(1e9, 107e6, 0.0), 107e6)

    def test_limit_approaches_exact_when_kp_small(self):
        a = td.circuit_exact(G_AL, 1e5 * KP_AL, KP_AL, 107e6, 210e6)
        b = td.circuit_limit(math.sqrt(G_AL * KP_AL), 107e6, 210e6)
        self.assertAlmostEqual(a / b, 1.0, places=3)

    def test_bad_input(self):
        with self.assertRaises(ValueError):
            td.circuit_exact(0.0, KE_AL, KP_AL, 0.0, 1e8)


class TestSergeev(unittest.TestCase):
    def test_g_from_lambda(self):
        self.assertAlmostEqual(td.g_from_lambda(GAM_AL, LW2_AL) / 1e17, 4.0687, delta=0.0005)

    def test_sigma_e_values(self):
        self.assertAlmostEqual(td.sergeev_from_g(GAM_AL, G_AL, UL_AL, UT_AL, NE_AL, 300.0) / MW, 177.48, delta=0.2)
        lam = td.sergeev_from_lambda(GAM_AL, LW2_AL, UL_AL, UT_AL, NE_AL, 300.0)
        self.assertAlmostEqual(lam / MW, 293.5, delta=0.5)
        noz = td.sergeev_from_lambda(GAM_AL, LW2_AL, UL_AL, UT_AL, NE_AL, 300.0, zeta=False)
        self.assertAlmostEqual(noz / lam, td.ZETA3, places=12)

    def test_linear_in_T(self):
        a = td.sergeev_from_g(GAM_AL, G_AL, UL_AL, UT_AL, NE_AL, 300.0)
        b = td.sergeev_from_g(GAM_AL, G_AL, UL_AL, UT_AL, NE_AL, 700.0)
        self.assertAlmostEqual(b / a, 7.0 / 3.0, places=12)


class TestBounds(unittest.TestCase):
    def test_linear_bound(self):
        s, t = td.linear_bound(AL_SI)
        self.assertAlmostEqual(s / MW, 530.0 / 700.0, places=9)
        self.assertEqual(t, 700)
        s, t = td.linear_bound(AL_GAN)
        self.assertAlmostEqual(s / MW, 0.78, places=9)

    def test_floor_bound(self):
        s, t = td.floor_bound(AL_SI, 1.26e9, 0.5, 210e6)
        self.assertAlmostEqual(s / MW, 0.50600, delta=0.0005)
        self.assertEqual(t, 600)


class TestInfer(unittest.TestCase):
    def test_dsf_size_channel_flattens_al_si(self):
        r = td.infer_with_direct(AL_SI, 1.26e9, 0.5, S_DSF)
        want = [334.544, 323.378, 351.043, 355.130, 397.254]
        for got, w in zip(r["G_pp"], want):
            self.assertAlmostEqual(got / MW, w, delta=0.1)

    def test_dsf_size_channel_al_gan(self):
        r = td.infer_with_direct(AL_GAN, 1.26e9, 0.5, S_DSF)
        want = [402.045, 375.577, 354.490, 356.941, 430.167]
        for got, w in zip(r["G_pp"], want):
            self.assertAlmostEqual(got / MW, w, delta=0.15)

    def test_large_channel_makes_phonon_side_fall(self):
        r = td.infer_with_direct(AL_SI, 1.26e9, 0.5, 0.5e6)
        self.assertAlmostEqual(r["G_pp"][0] / MW, 271.78, delta=0.1)
        self.assertAlmostEqual(r["G_pp"][3] / MW, 215.24, delta=0.1)

    def test_channel_larger_than_G_gives_none(self):
        r = td.infer_with_direct(AL_SI, 1.26e9, 0.5, 0.8e6)
        self.assertIsNone(r["G_pp"][4])


class TestMinRise(unittest.TestCase):
    def test_sample2_predictions(self):
        a = td.min_rise_prediction(330e6, 400.0, 600.0, S_DSF, 1.26e9, 0.5)
        self.assertAlmostEqual(a["G_min"] / MW, 395.28, delta=0.1)
        self.assertAlmostEqual(a["sigma_p_T1"] / MW, 220.05, delta=0.05)
        b = td.min_rise_prediction(330e6, 400.0, 600.0, 0.59e6, 1.26e9, 0.5)
        self.assertAlmostEqual(b["G_min"] / MW, 446.45, delta=0.1)

    def test_no_room(self):
        self.assertIsNone(td.min_rise_prediction(100e6, 400.0, 600.0, 0.59e6, 1.26e9, 0.5)["G_min"])


if __name__ == "__main__":
    unittest.main()
