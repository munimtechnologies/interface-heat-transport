# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""Tests for tbc_diag (pure-Python parts; the tbc_bte check is skipped if numpy or tbc_bte is absent).
Every expected number is worked by hand in the comments."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import tbc_diag as d  # noqa: E402


class TestHelpers(unittest.TestCase):
    def test_fcc_volume(self):
        # Si a = 5.4608 A: 5.4608^3 = 162.842... ; /4 = 40.7106 A^3 = 4.07106e-29 m^3
        self.assertAlmostEqual(d.fcc_prim_volume_m3(5.4608) / 4.07106e-29, 1.0, delta=1e-4)

    def test_uniform_lengths(self):
        # K = (1, 3, 0); kappa 8 -> l = 8 / (2 x 4) = 1 on bands with K > 0
        self.assertEqual(d.uniform_lengths([1.0, 3.0, 0.0], 8.0), [1.0, 1.0, 0.0])
        # check sum 2 K l = kappa: 2 (1 + 3) x 1 = 8

    def test_shared_share(self):
        # side: K (1, 1), l (1, 3) -> contributions 2, 6; other side has only band 0 -> 2 / 8
        self.assertAlmostEqual(d.shared_share([1.0, 1.0], [1.0, 3.0], [5.0, 0.0]), 0.25)

    def test_mean_length(self):
        # weights 2Kl = 2, 6; lengths 1, 3 -> (2x1 + 6x3) / 8 = 20 / 8 = 2.5
        self.assertAlmostEqual(d.mean_length([1.0, 1.0], [1.0, 3.0]), 2.5)
        # shared only band 0 -> 1
        self.assertAlmostEqual(d.mean_length([1.0, 1.0], [1.0, 3.0], [5.0, 0.0]), 1.0)

    def test_within(self):
        self.assertTrue(d.within(43.9, 43.0, 1.0, 0.02))      # |0.9| <= 1
        self.assertFalse(d.within(45.0, 43.0, 1.0, 0.02))     # 2 > max(1, 0.86)
        self.assertTrue(d.within(368.0, 362.0, 1.0, 0.02))    # 6 <= 7.24


def _row(base, dmm, uni, uni1, uni2, x05=None, x2=None, n800=None):
    return {"base": base, "dmm": dmm, "uni": uni, "uni1": uni1, "uni2": uni2,
            "x0.5": base if x05 is None else x05, "x2": base if x2 is None else x2,
            "n800": base if n800 is None else n800}


class TestGrade(unittest.TestCase):
    def _rows(self, uni_factor=1.5):
        rows = {}
        for p, by in d.PRINTED.items():
            rows[p] = {}
            for T, (gb, gd) in by.items():
                rows[p][T] = _row(gb, gd, gd * uni_factor, gb * 1.2, gd * 1.1)
        return rows

    def test_all_pass(self):
        g = d.grade(self._rows())
        self.assertEqual(g["D0"], "PASS")
        self.assertEqual(g["D1"], "PASS")
        self.assertEqual(g["D2"], "PASS")
        self.assertTrue(g["D3"].startswith("PASS"))
        # D4: uni2 recovers more than uni1 -> partner side named
        self.assertIn("-> partner", g["D4"])

    def test_d3_fail(self):
        g = d.grade(self._rows(uni_factor=0.5))
        self.assertTrue(g["D3"].startswith("FAIL"))

    def test_d0_fail_and_d1_fail(self):
        rows = self._rows()
        rows["TiN/C"][300]["base"] = 50.0          # 7 away from 43 -> D0 FAIL
        rows["TiN/Si"][700]["x2"] = 61.0 * 1.05     # 5 % -> D1 FAIL
        g = d.grade(rows)
        self.assertTrue(g["D0"].startswith("FAIL"))
        self.assertIn("TiN/C 300 G 50.0", g["D0"])
        self.assertTrue(g["D1"].startswith("FAIL"))

    def test_missing(self):
        g = d.grade({})
        self.assertTrue(g["D0"].startswith("FAIL"))
        self.assertTrue(g["D3"].startswith("UNDETERMINED"))


class TestGrayLimit(unittest.TestCase):
    """One band: tbc_bte must give the Simons value, above DMM (needs numpy and tbc_bte)."""

    def test_gray(self):
        try:
            import numpy as np
            import tbc_bte as tb
        except ImportError:
            self.skipTest("numpy or tbc_bte not present")
        # K1 = 1e9, K2 = 1e9 W/m2K, lengths 1e-7 m: DMM = 5e8; kappa = 2 x 1e9 x 1e-7 = 200 each;
        # Simons: 1/G = 1/5e8 - 1/(2e9) - 1/(2e9) = 2e-9 - 1e-9 = 1e-9 -> G = 1e9
        s = tb.solve_two_stream(np.array([1e9]), np.array([1e9]), np.array([1e-7]), np.array([1e-7]),
                                n_nodes=400)
        self.assertAlmostEqual(s["G_dmm"] / 5e8, 1.0, delta=1e-9)
        self.assertGreater(s["G"], s["G_dmm"])
        self.assertAlmostEqual(s["G"] / s["G_simons"], 1.0, delta=0.03)


if __name__ == "__main__":
    unittest.main()
