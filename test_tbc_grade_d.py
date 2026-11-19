# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Made with AI agents working in MV2's name. The numbers in results/dtest/grade-d1d4.json were worked
# out by hand (PLAN-amendment-0.4.md shows the arithmetic); this test recomputes them in code.
import json
import os
import unittest

import tbc_grade_d as gd

HERE = os.path.dirname(os.path.abspath(__file__))
MW = 1e6


def _load(*p):
    with open(os.path.join(HERE, *p)) as f:
        return json.load(f)


READ = _load("data", "read-d1d3.json")
MEAS = _load("data", "meas-li2022.json")
HAND = _load("results", "dtest", "grade-d1d4.json")
RES = gd.run(READ, MEAS)
CASES = ("D1_D2_li_sample2", "D3_cheng_al_sapphire")
TOL = {"sergeev_g_DFT": 0.3, "DSF_107_linear": 0.3, "sergeev_table1": 0.5, "sergeev_lambda": 2.0,
       "sergeev_eq4_printed": 2.0}


class TestSlopes(unittest.TestCase):
    def test_slopes(self):
        for k, v in HAND["slopes_MW_per_K"].items():
            self.assertAlmostEqual(RES["slopes_MW_per_K"][k], v, delta=0.002, msg=k)


class TestGrades(unittest.TestCase):
    def test_g_min_and_verdicts(self):
        for c in CASES:
            for k in TOL:
                got, want = RES["cases"][c]["by_slope"][k], HAND[c][k]
                self.assertEqual(got["excluded"], want["excluded"], msg=(c, k))
                if want["G_min"] is None:
                    self.assertIsNone(got["G_min_MW"], msg=(c, k))
                else:
                    self.assertAlmostEqual(got["G_min_MW"], want["G_min"], delta=TOL[k], msg=(c, k))
                if "ratio" in want:
                    self.assertAlmostEqual(got["ratio"], want["ratio"], delta=0.001, msg=(c, k))
                    self.assertAlmostEqual(got["threshold_MW"], want["threshold"], delta=0.3, msg=(c, k))

    def test_worst_case_reading(self):
        for c in CASES:
            for k, want in HAND[c]["worst_case_reading"].items():
                got = RES["cases"][c]["worst_case_reading"][k]
                self.assertEqual(got["G1_MW"], want["G1"])
                self.assertEqual(got["G2_MW"], want["G2"])
                self.assertAlmostEqual(got["G_min_MW"], want["G_min"], delta=0.3, msg=(c, k))
                self.assertEqual(got["excluded"], want["excluded"])
                self.assertTrue(got["verdict_holds"])

    def test_allowed_slopes(self):
        for c in CASES:
            got = RES["cases"][c]
            self.assertAlmostEqual(got["s_allowed_central_MW_per_K"], HAND[c]["s_allowed_central"], delta=0.002)
            self.assertAlmostEqual(got["s_allowed_cut085_MW_per_K"], HAND[c]["s_allowed_cut085"], delta=0.003)
            self.assertAlmostEqual(got["sigma_e_300_allowed_central_MW"], HAND[c]["sigma_e_300_allowed_central"],
                                   delta=0.6)

    def test_grade_words(self):
        g = RES["grades"]
        self.assertEqual(g["D1"], "EXCLUDED")
        self.assertEqual(g["D2"], "NOT EXCLUDED")
        self.assertEqual(g["D3_sergeev_g_DFT"], "EXCLUDED")
        self.assertEqual(g["D3_DSF_107_linear"], "NOT EXCLUDED")
        self.assertTrue(g["D4"].startswith("TRIGGERED"))
        self.assertTrue(g["all_sergeev_versions_excluded_both"])
        self.assertTrue(RES["cases"]["D3_cheng_al_sapphire"]["graded"])
        for k in HAND["grades"]:
            self.assertTrue(HAND["grades"][k] if isinstance(HAND["grades"][k], bool)
                            else HAND["grades"][k].startswith(g[k].split(":")[0]), msg=k)


class TestExploratory(unittest.TestCase):
    def test_sample1_phonon_side(self):
        got = RES["sample1_phonon_side_after_cap"]
        want = HAND["exploratory_not_preregistered"]
        self.assertAlmostEqual(got["s_cap_MW_per_K"], want["s_cap"], delta=0.002)
        for a, b in zip(got["G_pp_MW"], want["sample1_phonon_side"]["G_pp"]):
            self.assertAlmostEqual(a, b, delta=1.0)
        self.assertAlmostEqual(got["ratio_700_400"], want["sample1_ratio_700_400"], delta=0.01)
        self.assertAlmostEqual(RES["cheng_ep_share_300K"], want["cheng_ep_share_300K"], delta=0.001)


class TestSolver(unittest.TestCase):
    def test_g_min_grows_with_s(self):
        # the bisection in s_allowed needs G_min(T2) to rise monotonically with s
        for g1, t1, t2 in [(329e6, 400.0, 600.0), (346e6, 300.0, 480.0), (379e6, 300.0, 700.0)]:
            prev = None
            for i in range(0, 200):
                s = i / 200.0 * g1 / t1 * 0.999
                v = gd.g_min(g1, t1, t2, s)
                if prev is not None:
                    self.assertGreater(v, prev)
                prev = v

    def test_s_allowed_inverts_g_min(self):
        s = 0.3e6
        g2 = gd.g_min(329e6, 400.0, 600.0, s)
        self.assertAlmostEqual(gd.s_allowed(329e6, 400.0, g2, 600.0, 1.0) / MW, 0.3, delta=1e-6)
        self.assertEqual(gd.s_allowed(329e6, 400.0, 200e6, 600.0, 1.0), 0.0)

    def test_none_when_channel_exceeds_g1(self):
        r = gd.grade(329e6, 400.0, 347e6, 600.0, 0.9e6)
        self.assertIsNone(r["G_min"])
        self.assertTrue(r["excluded"])


if __name__ == "__main__":
    unittest.main()
