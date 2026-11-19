# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Made with AI agents working in MV2's name. Expected numbers below were worked out by hand.
import math
import unittest

import tbc_gep as tg

AL_SI = {"300": 379e6, "400": 400e6, "500": 446e6, "600": 478e6, "700": 530e6}
AL_GAN = {"300": 422e6, "400": 432e6, "500": 448e6, "600": 479e6, "700": 546e6}
MW = 1e6


class TestGep(unittest.TestCase):
    def test_gep_mr(self):
        # sqrt(2.5e17 * 6) = sqrt(1.5e18) = 1.22474e9
        self.assertAlmostEqual(tg.gep_mr(2.5e17, 6.0) / 1e9, 1.224745, places=5)
        with self.assertRaises(ValueError):
            tg.gep_mr(0.0, 6.0)

    def test_lam_p_T(self):
        self.assertAlmostEqual(tg.lam_p_T(6.0, 600.0), 3.0, places=12)
        # Matthiessen with a 60 W/m/K phonon-electron limit: 1/(1/3 + 1/60) = 1/0.35
        self.assertAlmostEqual(tg.lam_p_T(6.0, 600.0, lam_pe=60.0), 1.0 / 0.35, places=12)

    def test_power_law(self):
        self.assertAlmostEqual(tg.gep_power_law(1.26e9, 400.0), 1.26e9, places=3)
        self.assertAlmostEqual(tg.gep_power_law(1.26e9, 700.0) / 1e9, 1.26 * math.sqrt(4 / 7), places=12)
        self.assertAlmostEqual(tg.gep_power_law(1.21e9, 700.0, n=0.0), 1.21e9, places=3)

    def test_phonon_side_and_series(self):
        gpp = tg.phonon_side(400e6, 1.26e9)
        self.assertAlmostEqual(gpp / MW, 586.047, places=2)
        self.assertAlmostEqual(tg.series(gpp, 1.26e9) / MW, 400.0, places=6)
        with self.assertRaises(ValueError):
            tg.phonon_side(2e9, 1.26e9)

    def test_amplification(self):
        gep = tg.gep_power_law(1.26e9, 700.0)  # 952.47 MW
        self.assertAlmostEqual(tg.amplification(530e6, gep), 2.2545, places=3)
        self.assertAlmostEqual(tg.amplification(400e6, 1.26e9), 1.0 / (1.0 - 400 / 1260), places=12)

    def test_linear_floor_matches_amendment_031(self):
        # No series term: a = 130/300 MW/K, E = (2800 - 2120)/3 = 226.67 MW, share 0.4333
        e, a, share = tg.linear_floor(400.0, 530.0)
        self.assertAlmostEqual(e, 226.6667, places=3)
        self.assertAlmostEqual(a, 0.43333, places=4)
        self.assertAlmostEqual(share, 0.43333, places=4)

    def test_infer_al_si_t_half(self):
        r = tg.infer(AL_SI, 1.26e9, 0.5)
        want = [512.505, 586.047, 738.10, 892.83, 1194.90]
        for got, w in zip(r["G_pp"], want):
            self.assertAlmostEqual(got / MW, w, delta=0.05)
        self.assertTrue(r["convex_everywhere"])
        sec = [s / MW for s in r["second_differences"]]
        self.assertAlmostEqual(sec[0], 78.52, delta=0.1)
        self.assertAlmostEqual(sec[1], 2.67, delta=0.1)
        self.assertAlmostEqual(sec[2], 147.33, delta=0.1)
        self.assertLess(r["linear_400_700"]["E"], 0.0)

    def test_infer_al_gan_t_half(self):
        r = tg.infer(AL_GAN, 1.26e9, 0.5)
        want = [594.41, 657.39, 743.60, 896.33, 1279.43]
        for got, w in zip(r["G_pp"], want):
            self.assertAlmostEqual(got / MW, w, delta=0.05)
        self.assertTrue(r["convex_everywhere"])
        inc = [d / MW for d in r["increments"]]
        self.assertAlmostEqual(inc[3], 383.10, delta=0.1)

    def test_infer_constant_gep_not_convex(self):
        # Amendment 0.3.1 values with a constant 1.21 GW: 551.85 at 300 K; increments 46, 109, 83, 153
        r = tg.infer(AL_SI, 1.21e9, 0.0)
        self.assertAlmostEqual(r["G_pp"][0] / MW, 551.85, delta=0.05)
        self.assertFalse(r["convex_everywhere"])


if __name__ == "__main__":
    unittest.main()
