# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""Tests for tbc_fill. Pure-Python parts always run; the tbc_bte check is skipped if numpy or
tbc_bte is absent. Every expected number is worked by hand in the comments."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import tbc_fill as f  # noqa: E402


class TestHelpers(unittest.TestCase):
    def test_first_have(self):
        self.assertEqual(f.first_have([0.0, 1.0, 2.0], [0.0, 0.0, 3.0]), 2)
        self.assertIsNone(f.first_have([1.0, 0.0], [0.0, 0.0]))

    def test_merge_low(self):
        # side a first sampled bin 2, side b first sampled bin 1 -> j = 2; bins 0..2 summed
        Ka, ka, Kb, kb, j = f.merge_low([1, 2, 3, 4], [0, 0, 6, 8], [0, 5, 6, 7], [0, 1, 2, 3])
        self.assertEqual(j, 2)
        self.assertEqual(Ka, [0.0, 0.0, 6.0, 4.0])     # 1 + 2 + 3 = 6
        self.assertEqual(ka, [0.0, 0.0, 6.0, 8.0])     # 0 + 0 + 6 = 6
        self.assertEqual(Kb, [0.0, 0.0, 11.0, 7.0])    # 0 + 5 + 6 = 11
        self.assertEqual(kb, [0.0, 0.0, 3.0, 3.0])     # 0 + 1 + 2 = 3
        # totals conserved: 10, 14, 18, 6
        self.assertEqual((sum(Ka), sum(ka), sum(Kb), sum(kb)), (10.0, 14.0, 18.0, 6.0))

    def test_merge_low_no_sample(self):
        with self.assertRaises(ValueError):
            f.merge_low([1, 1], [0, 0], [1, 1], [1, 1])

    def test_filled_flags(self):
        self.assertEqual(f.filled_flags([1.0, 1.0, 0.0], [0.0, 2.0, 0.0]), [True, False, False])

    def test_two_stream_kappa(self):
        # 2 x 1 x 3 + 2 x 2 x 4 = 6 + 16 = 22 (K = 0 band ignored)
        self.assertAlmostEqual(f.two_stream_kappa([1.0, 2.0, 0.0], [3.0, 4.0, 5.0]), 22.0)

    def test_filled_share(self):
        # bin 0 filled: 2 x 1 x 3 = 6; bin 1: 2 x 1 x 1 = 2 -> 6 / 8
        self.assertAlmostEqual(f.filled_share([1.0, 1.0], [3.0, 1.0], [0.0, 2.0]), 0.75)

    def test_uniform_lengths(self):
        # K (1, 3, 0), kappa 8 -> 8 / (2 x 4) = 1
        self.assertEqual(f.uniform_lengths([1.0, 3.0, 0.0], 8.0), [1.0, 1.0, 0.0])

    def test_within(self):
        self.assertTrue(f.within(43.7, 43.3, 0.5, 0.01))     # 0.4 <= 0.5
        self.assertFalse(f.within(44.0, 43.3, 0.5, 0.01))    # 0.7 > max(0.5, 0.433)


def _rows(merge_of, kap_merge=44.6):
    rows = {}
    for p, by in f.DIAG.items():
        rows[p] = {}
        for T, (gb, gd) in by.items():
            rows[p][T] = {"base": gb, "dmm": gd, "merge": merge_of(gb, gd), "dmm_merge": gd,
                          "kap1_merge": kap_merge, "kap1_bins": 44.6}
    return rows


class TestGrade(unittest.TestCase):
    def test_all_pass(self):
        g = f.grade(_rows(lambda gb, gd: gd), {"fill_share": 0.66, "fill_lo_max_thz": 2.0})
        self.assertEqual(g["F0"], "PASS")
        self.assertTrue(g["F1"].startswith("PASS"))
        self.assertTrue(g["F2"].startswith("PASS"))   # gd / gb >= 5 everywhere
        self.assertTrue(g["F3"].startswith("PASS"))   # merge = dmm_merge
        self.assertEqual(g["F4"], "PASS")

    def test_fails(self):
        # merge 1.5 x base: F2 fails; 1.5 x 43.3 x 1.6 = 103.9 < 264.3: F3 fails;
        # kappa 50 vs 44.6 is +12 %: F4 fails; a filled bin at 12 THz: F1 fails
        g = f.grade(_rows(lambda gb, gd: 1.5 * gb, kap_merge=50.0),
                    {"fill_share": 0.66, "fill_lo_max_thz": 12.0})
        self.assertTrue(g["F1"].startswith("FAIL"))
        self.assertTrue(g["F2"].startswith("FAIL"))
        self.assertTrue(g["F3"].startswith("FAIL"))
        self.assertTrue(g["F4"].startswith("FAIL"))

    def test_f0_fail_and_missing(self):
        rows = _rows(lambda gb, gd: gd)
        rows["TiN/C"][300]["base"] = 44.0             # 0.7 away from 43.3
        g = f.grade(rows, {"fill_share": 0.66, "fill_lo_max_thz": 2.0})
        self.assertTrue(g["F0"].startswith("FAIL"))
        g = f.grade({}, None)
        self.assertTrue(g["F0"].startswith("FAIL"))
        self.assertTrue(g["F1"].startswith("UNDETERMINED"))


class TestMechanism(unittest.TestCase):
    """The filled-bin inflation and its removal by merging (needs numpy and tbc_bte)."""

    def test_inflation_and_merge(self):
        try:
            import numpy as np
            import tbc_bte as tb
        except ImportError:
            self.skipTest("numpy or tbc_bte not present")
        Ka = np.array([1e6, 2e6, 3e6])
        ka = np.array([0.0, 0.0, 10.0])
        # job rule: bins 0 and 1 get bin 2's length 10 / 6e6 (constant extrapolation)
        la, nf = tb.lengths_from_kappa(ka, Ka)
        self.assertEqual(nf, 2)
        # model kappa 2 x (1 + 2 + 3)e6 x 10 / 6e6 = 20, twice the input 10; filled share 0.5
        self.assertAlmostEqual(f.two_stream_kappa(list(Ka), list(la)), 20.0, places=9)
        self.assertAlmostEqual(f.filled_share(list(Ka), list(la), list(ka)), 0.5, places=12)
        Kb = [0.0, 1e6, 1e6]
        kb = [0.0, 5.0, 5.0]
        Ka2, ka2, Kb2, kb2, j = f.merge_low(list(Ka), list(ka), Kb, kb)
        self.assertEqual(j, 2)
        la2, nf2 = tb.lengths_from_kappa(np.array(ka2), np.array(Ka2))
        self.assertEqual(nf2, 0)
        # merged: one band K 6e6, kappa 10 -> model kappa exactly 10
        self.assertAlmostEqual(f.two_stream_kappa(Ka2, list(la2)), 10.0, places=9)


if __name__ == "__main__":
    unittest.main()
