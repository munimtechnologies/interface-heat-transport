# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""Unit tests of tbc_dressed. Every number below was first worked out by hand
(PLAN-amendment-0.5.1.md); the tests recompute it."""
import json
import os
import random
import unittest

import tbc_dressed as td
import tbc_nonmetal as tn

HERE = os.path.dirname(os.path.abspath(__file__))


class TestReading(unittest.TestCase):
    def test_constants_match_data_file(self):
        d = json.load(open(os.path.join(HERE, "data", "read-gaskins2018.json")))
        for T, (g, up, dn) in td.PUB.items():
            p = d["points_used"][str(T)]
            self.assertEqual((float(p["G"]), float(p["up"]), float(p["down"])), (g, up, dn))
        self.assertEqual(float(d["model_curves"]["DMM"]["300"]), td.PUB_DMM[300])
        self.assertEqual(float(d["model_curves"]["DMM"]["500"]), td.PUB_DMM[500])
        self.assertEqual(float(d["model_curves"]["AGF"]["500"]), td.PUB_AGF[500])
        self.assertEqual(float(d["reading_uncertainty_points"]), td.PUB_READ_ERR)


class TestPieces(unittest.TestCase):
    def test_isotonic(self):
        self.assertEqual(td.isotonic([1, 3, 2], [1, 1, 1]), [1.0, 2.5, 2.5])
        self.assertEqual(td.isotonic([3, 2, 1], [1, 1, 2]), [1.75] * 3)
        self.assertEqual(td.isotonic([1, 2, 3], [5, 1, 1]), [1.0, 2.0, 3.0])

    def test_sigma(self):
        self.assertAlmostEqual(td.sigma_from_95(107, 111), 55.6133, delta=1e-3)
        self.assertAlmostEqual(td.sigma_from_95(121, 124), 62.5011, delta=1e-3)
        self.assertAlmostEqual(td.sigma_from_95(127, 131), 65.8175, delta=1e-3)

    def test_bound_general_reduces(self):
        for g1, g2 in ((491, 502), (471.36, 522.08), (520, 502)):
            self.assertAlmostEqual(td.bound_general(g1, g2, 300, 500, 1.0),
                                   tn.channel_bound(g1, g2, 300, 500), places=9)
        # tight: E falls by kappa, I = c T
        E, c, k = 400.0, 0.3, 0.9
        g1, g2 = E + c * 300, k * E + c * 500
        self.assertAlmostEqual(td.bound_general(g1, g2, 300, 500, k), c * 300, places=9)
        self.assertEqual(td.bound_general(400, 500, 300, 500, 2.0), float("inf"))

    def test_cap_ratio_grows_with_frequency(self):
        prev = 0.0
        for k in range(1, 400):
            r = td.cap_ratio(0.1 * k, 300.0, 500.0)
            self.assertGreater(r, prev)
            prev = r
        self.assertAlmostEqual(td.cap_ratio(0.01, 300.0, 500.0), 1.0, places=6)
        self.assertAlmostEqual(td.cap_ratio(22.0, 300.0, 500.0), 1.793239, delta=1e-4)


class TestDressed(unittest.TestCase):
    def test_identity_random(self):
        rnd = random.Random(20261119)
        for _ in range(50):
            n = 30
            cut = rnd.randint(5, 25)
            M1 = [rnd.uniform(0.1, 3.0) if i < cut else 0.0 for i in range(n)]
            M2 = [rnd.uniform(0.1, 3.0) for _ in range(n)]
            C = [rnd.uniform(0.05, 1.0) for _ in range(n)]
            d = td.dressed_dmm(M1, M2, C)
            self.assertAlmostEqual(d["Gd_shi"] / d["Gd_closed"], 1.0, places=12)
            self.assertGreater(d["Gd_closed"], d["G0"])

    def test_gray_doubles(self):
        # one shared band, equal modes: G0 = 0.5, R1 = R2 = 1, G_d = 1 / (2 - 0.5 - 0.5) = 1
        d = td.dressed_dmm([1.0], [1.0], [1.0])
        self.assertAlmostEqual(d["G0"], 0.5, places=12)
        self.assertAlmostEqual(d["Gd_closed"], 1.0, places=12)

    def test_two_bin_fall(self):
        # bin A (1 THz) on both sides, bin B (22 THz) only on side 2. Every C grows, G0 grows,
        # yet G_d falls: hand values 0.8750 (300 K) and 0.8287 (500 K).
        out = {}
        for T in (300.0, 500.0):
            cA = tn.cap(tn.x_of(1.0, T))
            cB = tn.cap(tn.x_of(22.0, T))
            out[T] = td.dressed_dmm([1.0, 0.0], [1.0, 1.0], [cA, cB])
        self.assertGreaterEqual(out[500.0]["G0"], out[300.0]["G0"])
        self.assertAlmostEqual(out[300.0]["Gd_closed"], 0.8750, delta=2e-4)
        self.assertAlmostEqual(out[500.0]["Gd_closed"], 0.8287, delta=2e-4)
        self.assertLess(out[500.0]["Gd_closed"], out[300.0]["Gd_closed"])

    def test_crude_floor(self):
        self.assertAlmostEqual(td.min_kappa_e_dmm(22.0, 300.0, 500.0), 0.693314, delta=1e-4)
        self.assertAlmostEqual(td.min_kappa_e_dmm(22.0, 300.0, 500.0, dress_ratio=1.0), 1.0, places=12)


class TestNumbers(unittest.TestCase):
    hand = {"profile_I300_upper95": (229.05, 0.3), "profile_chi2_min": (0.05353, 2e-4),
            "profile_a_at_min": (0.0, 1e-3), "profile_share300": (0.46745, 1e-3),
            "two_point_I300_upper95": (232.10, 0.1), "correlated_I300_central_pub": (19.5, 1e-9),
            "correlated_I300_reading_worst_pub": (49.5, 1e-9), "dmm_pub_minus_arxiv_300": (10.0, 1e-9),
            "gan_top_cap_ratio_300_500": (1.793239, 1e-4), "kappa_e_floor_dmm_crude": (0.693314, 1e-4),
            "I300_bound_at_crude_floor_central_pub": (167.746, 0.05)}

    def test_hand_numbers(self):
        x = td.numbers()
        for k, (v, tol) in self.hand.items():
            self.assertAlmostEqual(x[k], v, delta=tol, msg=k)
        self.assertEqual(x["pub_minus_arxiv"], {"300": -1.0, "400": 1.0, "500": 1.0})
        self.assertEqual(x["twice_dmm_pub"], {"300": 500.0, "500": 560.0})

    def test_profile_pooled_closed_form(self):
        # in the all-pooled regime chi^2(a) = 0.165844 - 0.00239097 a + 5.98672e-5 a^2 (hand)
        temps = sorted(td.PUB)
        g = [td.PUB[T][0] for T in temps]
        sig = [109 / 1.96, 62.5, 129 / 1.96]
        for a in (100.0, 227.0, 300.0):
            hand = 0.16584359 - 0.00239097 * a + 5.986721e-5 * a * a
            self.assertAlmostEqual(td.chi2_profile(a, temps, g, sig), hand, delta=2e-4)

    def test_profile_recovers_planted_channel(self):
        # E flat at 300, I = 60 T/300, sigma 0.5. Every a <= 60 fits exactly (a rising E absorbs
        # the rest), so chi2_min = 0; above 60, chi^2 = 8 (a - 60)^2 / 9, crossing 2.7055 at 61.745.
        temps = [300.0, 400.0, 500.0]
        g = [300 + 60 * T / 300 for T in temps]
        a, cmin, amin = td.profile_upper_limit(temps, g, [0.5, 0.5, 0.5])
        self.assertAlmostEqual(cmin, 0.0, delta=1e-9)
        self.assertAlmostEqual(a, 61.745, delta=0.01)

    def test_results_file_matches(self):
        r = json.load(open(os.path.join(HERE, "results", "nonmetal", "bounds-0.5.1.json")))
        x = td.numbers()
        for k, (v, tol) in self.hand.items():
            self.assertAlmostEqual(r[k], x[k], delta=tol, msg=k)


class TestSupplementaryGrades(unittest.TestCase):
    def prim(self, b=(450.0, 460.0, 455.0), inel=20.0, dmm=250.0):
        return {"pairs": {"ZnO/GaN": {"300": {"G_bte": b[0], "G_inel": inel, "G_dmm": dmm},
                                      "400": {"G_bte": b[1]}, "500": {"G_bte": b[2]}}}}

    def test_pass(self):
        g = td.grade_supp(self.prim())
        self.assertTrue(g["Z1g"].startswith("PASS"), g["Z1g"])     # max 69.38 (hand)
        self.assertIn("69.38", g["Z1g"])
        self.assertTrue(g["Z6"].startswith("PASS"))
        self.assertTrue(g["Z7"].startswith("PASS"))                # 1.8

    def test_fail(self):
        g = td.grade_supp(self.prim(b=(300.0, 320.0, 340.0), inel=60.0, dmm=300.0))
        self.assertTrue(g["Z1g"].startswith("FAIL"), g["Z1g"])     # max 0
        self.assertTrue(g["Z6"].startswith("FAIL"))
        self.assertTrue(g["Z7"].startswith("FAIL"))                # 1.0

    def test_missing(self):
        g = td.grade_supp({"pairs": {}})
        for k in ("Z1g", "Z6", "Z7"):
            self.assertTrue(g[k].startswith("UNDETERMINED"), k)


if __name__ == "__main__":
    unittest.main()
