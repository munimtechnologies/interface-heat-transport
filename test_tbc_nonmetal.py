# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""Unit tests of tbc_nonmetal. Every number below was first worked out by hand
(PLAN-amendment-0.5.md); the tests recompute it."""
import json
import os
import unittest

import tbc_nonmetal as tn

HERE = os.path.dirname(os.path.abspath(__file__))


class TestReadings(unittest.TestCase):
    def test_constants_match_data_file(self):
        d = json.load(open(os.path.join(HERE, "data", "read-nonmetal.json")))
        z = d["zno_gan"]
        for T in (100, 200, 300, 400, 500):
            self.assertEqual(float(z["G"][str(T)]), tn.ZNO_GAN[T])
        self.assertEqual(z["reading_uncertainty_rel"], tn.ZNO_GAN_REL)
        self.assertEqual(float(z["model_curves"]["AGF"]["300"]), tn.ZNO_GAN_AGF[300])
        self.assertEqual(float(z["model_curves"]["DMM"]["500"]), tn.ZNO_GAN_DMM[500])
        self.assertEqual(float(d["gan_aln_ucla"]["G"]["300"]), tn.UCLA_GAN_ALN_300)
        for pair, pts in tn.WALWIL.items():
            for T, v in pts.items():
                self.assertEqual(float(d["walwil"]["G"][pair][str(T)]), v)


class TestGrowthFloor(unittest.TestCase):
    """C * Gamma of every three-phonon term grows at least like T (the bound's premise)."""

    def grid(self):
        return [0.05 * k for k in range(1, 601)]          # 0.05 .. 30 THz

    def test_decay(self):
        for T1, T2 in ((300.0, 500.0), (300.0, 400.0), (400.0, 500.0)):
            floor = T2 / T1 * (1.0 - 1e-9)
            worst = 1e9
            for nu in self.grid()[::3]:
                for f in (0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5):
                    worst = min(worst, tn.growth_decay(nu, f, T1, T2))
            self.assertGreaterEqual(worst, floor)
            # the floor is approached in the classical limit
            self.assertAlmostEqual(tn.growth_decay(0.01, 0.5, T1, T2), T2 / T1, places=3)

    def test_combination(self):
        for T1, T2 in ((300.0, 500.0), (300.0, 400.0), (400.0, 500.0)):
            floor = T2 / T1 * (1.0 - 1e-9)
            worst = 1e9
            for nu in self.grid()[::5]:
                for nu1 in (0.05, 0.5, 2.0, 5.0, 10.0, 20.0):
                    worst = min(worst, tn.growth_combination(nu, nu1, T1, T2))
            self.assertGreaterEqual(worst, floor)

    def test_hand_values(self):
        # x = 2 at 300 K (nu = 12.5 THz): decay into halves, hand value 1.946 (cap 1.2266 x coth 1.586)
        nu = 2.0 * 300.0 / tn.H_OVER_KB_THZ
        self.assertAlmostEqual(tn.growth_decay(nu, 0.5, 300.0, 500.0), 1.946, delta=0.003)
        self.assertAlmostEqual(tn.cap(2.0), 0.7241, delta=1e-4)
        self.assertAlmostEqual(tn.cap(1.2), 0.8882, delta=1e-4)


class TestBound(unittest.TestCase):
    def test_simple(self):
        self.assertAlmostEqual(tn.channel_bound(491, 502, 300, 500), 16.5, places=9)
        self.assertEqual(tn.channel_bound(520, 502, 400, 500), 0.0)
        self.assertAlmostEqual(tn.channel_bound(491, 502, 300, 500, growth=0.5), 33.0, places=9)

    def test_bound_is_tight(self):
        # E constant, I = c T: the bound returns I(T1) exactly.
        E, c = 300.0, 0.2
        self.assertAlmostEqual(tn.channel_bound(E + c * 300, E + c * 500, 300, 500), c * 300, places=9)

    def test_exploratory_hand_numbers(self):
        x = tn.exploratory()
        z = x["zno_gan"]
        hand = {"I300_max_central_300_500": 16.5, "share300_central": 0.033605,
                "I300_max_worst_300_500": 76.08, "share300_worst": 0.161405,
                "I400_max_central_400_500": 0.0, "I400_max_worst_400_500": 91.52,
                "I300_max_worst_from_400_500": 68.64, "I300_max_worst_combined": 68.64,
                "I300_max_central_half_growth": 33.0, "I300_max_worst_half_growth": 152.16,
                "ratio_500_300_central": 1.022403, "ratio_500_300_worst": 1.107604,
                "agf_gap300": 191.0, "G500_min_if_agf_gap_is_channel": 623.3333,
                "dmm_gap300": 251.0, "G500_min_if_dmm_gap_is_channel": 696.3333}
        for k, v in hand.items():
            self.assertAlmostEqual(z[k], v, delta=1e-4 if v < 2 else 0.01, msg=k)
        w = x["walwil_ratios"]
        self.assertAlmostEqual(w["GaN/AlN_400_300"], 0.790890, delta=1e-5)
        self.assertAlmostEqual(w["AlN/SiC_350_300"], 1.16, delta=1e-9)
        self.assertAlmostEqual(w["GaN/SiC_450_300"], 1.456026, delta=1e-5)
        self.assertAlmostEqual(x["ucla_over_agf_300"], 0.936416, delta=1e-5)

    def test_results_file_matches(self):
        r = json.load(open(os.path.join(HERE, "results", "nonmetal", "bounds-exploratory.json")))
        x = tn.exploratory()
        for k, v in r["zno_gan"].items():
            self.assertAlmostEqual(x["zno_gan"][k], v, delta=1e-4 if abs(v) < 2 else 0.01, msg=k)
        for k, v in r["walwil_ratios"].items():
            self.assertAlmostEqual(x["walwil_ratios"][k], v, delta=1e-5, msg=k)
        self.assertAlmostEqual(x["ucla_over_agf_300"], r["ucla_over_agf_300"], delta=1e-5)
        for k in ("Z1_MAX", "Z2_MAX"):
            self.assertAlmostEqual(getattr(tn, k), r["thresholds"][k], delta=1e-4)


class TestGrades(unittest.TestCase):
    def test_thresholds(self):
        self.assertAlmostEqual(tn.Z1_MAX, 76.08, places=6)
        self.assertAlmostEqual(tn.Z2_MAX, 1.107604, delta=1e-6)
        self.assertEqual(tn.Z3_RANGE, (380.0, 640.0))
        self.assertEqual(tn.Z4_RANGE, (243.0, 405.0))

    def made_up(self, inel300=20.0, tot=(350.0, 370.0, 380.0), aln=(300.0, 310.0)):
        return {"materials": {"ZnO": {"relax_converged": True, "imag_fraction": 0.0, "max_thz": 12.8,
                                      "kappa300": 45.0}},
                "pairs": {"ZnO/GaN": {"300": {"G_inel": inel300, "G_total": tot[0]},
                                      "400": {"G_inel": inel300 * 1.4, "G_total": tot[1]},
                                      "500": {"G_inel": inel300 * 1.8, "G_total": tot[2]}},
                          "AlN/GaN": {"300": {"G_total": aln[0]}, "400": {"G_total": aln[1]}}}}

    def test_made_up_pass_fail(self):
        g = tn.grade_z(self.made_up())
        self.assertTrue(g["Z0"].startswith("PASS"))
        self.assertTrue(g["Z1"].startswith("PASS"))
        self.assertTrue(g["Z2"].startswith("PASS"))          # 380/350 = 1.0857
        self.assertTrue(g["Z3"].startswith("FAIL(low"))
        self.assertTrue(g["Z4"].startswith("PASS"))
        self.assertTrue(g["Z5"].startswith("PASS"))
        g = tn.grade_z(self.made_up(inel300=80.0, tot=(500.0, 540.0, 560.0), aln=(420.0, 400.0)))
        self.assertTrue(g["Z1"].startswith("FAIL"))
        self.assertTrue(g["Z2"].startswith("FAIL"))          # 1.12
        self.assertTrue(g["Z3"].startswith("PASS"))
        self.assertTrue(g["Z4"].startswith("FAIL(high"))
        self.assertTrue(g["Z5"].startswith("FAIL"))

    def test_missing(self):
        g = tn.grade_z({"materials": {}, "pairs": {}})
        for k in ("Z0", "Z1", "Z2", "Z3", "Z4", "Z5"):
            self.assertTrue(g[k].startswith("UNDETERMINED"), k)

    def test_z0_fails_on_bad_zno(self):
        p = self.made_up()
        p["materials"]["ZnO"]["max_thz"] = 16.0
        self.assertTrue(tn.grade_z(p)["Z0"].startswith("FAIL"))


if __name__ == "__main__":
    unittest.main()
