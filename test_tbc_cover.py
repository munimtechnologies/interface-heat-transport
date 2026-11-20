# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""Tests for tbc_cover. Pure-Python parts always run; the tbc_bte check is skipped if numpy or
tbc_bte is absent. Every expected number is worked by hand in the comments."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import tbc_cover as c  # noqa: E402


class TestBands(unittest.TestCase):
    def test_interleaved(self):
        # a sampled at 1, 3; b sampled at 0, 2 -> (0,1) closes at 1, (2,3) closes at 3
        Ka, ka = [1, 1, 1, 1], [0, 1, 0, 1]
        Kb, kb = [1, 1, 1, 1], [1, 0, 1, 0]
        bands = c.cover_bands(Ka, ka, Kb, kb)
        self.assertEqual(bands, [(0, 1), (2, 3)])
        # b's last sampled bin is 2; its bin 3 (K 1, kappa 0) joins band 1 -> K 1 + 1, kappa 1
        self.assertEqual(c.banded(Kb, kb, bands), ([2.0, 2.0], [1.0, 1.0]))
        self.assertEqual(c.banded(Ka, ka, bands), ([2.0, 2.0], [1.0, 1.0]))

    def test_first_band_is_merge_range(self):
        # the tbc_fill merge_low example: a first sampled at 2, b at 1 -> first band 0..2
        bands = c.cover_bands([1, 2, 3, 4], [0, 0, 6, 8], [0, 5, 6, 7], [0, 1, 2, 3])
        self.assertEqual(bands[0], (0, 2))
        self.assertEqual(bands, [(0, 2), (3, 3)])

    def test_cutoff_spec_vs_literal(self):
        # side b has states only in bins 0-1 (bin 1 unsampled); side a sampled everywhere
        Ka, ka = [1, 1, 1, 1, 1], [1, 1, 1, 1, 1]
        Kb, kb = [1, 1, 0, 0, 0], [1, 0, 0, 0, 0]
        spec = c.cover_bands(Ka, ka, Kb, kb)
        self.assertEqual(spec, [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4)])
        # b's bin 1 (above its last sampled bin 0) joins band 0: K 2, kappa 1; nothing else
        self.assertEqual(c.banded(Kb, kb, spec), ([2.0, 0.0, 0.0, 0.0, 0.0],
                                                  [1.0, 0.0, 0.0, 0.0, 0.0]))
        lit = c.cover_bands(Ka, ka, Kb, kb, literal=True)
        self.assertEqual(lit, [(0, 4)])       # whole tail folded into the last band
        self.assertEqual(c.banded(Ka, ka, lit, per_side_tail=False), ([5.0], [5.0]))

    def test_min_bins(self):
        K, k = [1] * 5, [1] * 5
        self.assertEqual(c.cover_bands(K, k, K, k, min_bins=2), [(0, 1), (2, 4)])
        self.assertEqual(c.cover_bands(K[:4], k[:4], K[:4], k[:4], min_bins=2), [(0, 1), (2, 3)])

    def test_conserves_totals(self):
        Ka, ka = [1, 2, 3, 4, 5, 6], [0, 0, 2, 0, 1, 0]
        Kb, kb = [3, 0, 1, 1, 0, 2], [1, 0, 0, 3, 0, 0]
        bands = c.cover_bands(Ka, ka, Kb, kb)
        for K, k in ((Ka, ka), (Kb, kb)):
            K2, k2 = c.banded(K, k, bands)
            self.assertAlmostEqual(sum(K2), sum(K))
            self.assertAlmostEqual(sum(k2), sum(k))
            for x, y in zip(K2, k2):           # no band with K > 0 is left without kappa
                self.assertFalse(x > 0 and not y > 0)

    def test_no_sample(self):
        with self.assertRaises(ValueError):
            c.cover_bands([1, 1], [0, 0], [1, 1], [1, 1])

    def test_widths(self):
        self.assertEqual(c.band_widths([(0, 1), (2, 2)], [0.0, 0.2, 0.4, 0.6]), [0.4, 0.6 - 0.4])


def _rows(cover_of, dmm_factor=1.0, kap_rel=1.0, coarse_rel=1.0, hi_shift=0.0):
    rows = {}
    for p, by in c.FILL.items():
        rows[p] = {}
        for T, (gb, gm, hi) in by.items():
            cv = cover_of(p, gm)
            rows[p][T] = {"base": gb, "merge": gm, "first_band_hi_thz": hi + hi_shift,
                          "cover": cv, "dmm_cover": cv * dmm_factor,
                          "kap1_cover": 44.6 * kap_rel, "kap1_bins": 44.6,
                          "kap2_cover": 20.0, "kap2_bins": 20.0, "cover_coarse": cv * coarse_rel}
    return rows


class TestGrade(unittest.TestCase):
    def test_all_pass(self):
        # cover = 1.1 x merge everywhere: TiN/C within 25 %, Si/Ge >= merge; G/DMM 1/1.25 = 0.8
        g = c.grade(_rows(lambda p, gm: 1.1 * gm, dmm_factor=1.25),
                    earlier={"ZnO/GaN": {"base": 190.0, "cover": 240.0}})   # 1.26 >= 1.2
        for k in ("V0", "V1", "V2", "V3", "V4", "V5"):
            self.assertTrue(g[k].startswith("PASS"), (k, g[k]))

    def test_fails(self):
        # Si/Ge cover 0.9 x merge: V2 fails; G/DMM 1/3 = 0.33 < 0.4: V3 fails; kappa +2 %: V1;
        # coarse +15 %: V4; ZnO/GaN 200/190 = 1.05: V5; first band shifted: V0
        g = c.grade(_rows(lambda p, gm: 0.9 * gm, dmm_factor=3.0, kap_rel=1.02, coarse_rel=1.15,
                          hi_shift=0.2),
                    earlier={"ZnO/GaN": {"base": 190.0, "cover": 200.0}})
        for k in ("V0", "V1", "V2", "V3", "V4", "V5"):
            self.assertTrue(g[k].startswith("FAIL"), (k, g[k]))

    def test_tin_c_band(self):
        # TiN/C 1.3 x merge fails V2 (outside 25 %) even though Si/Ge rows pass
        g = c.grade(_rows(lambda p, gm: (1.3 if p == "TiN/C" else 1.0) * gm))
        self.assertTrue(g["V2"].startswith("FAIL"))
        self.assertTrue(g["V5"].startswith("UNDETERMINED"))

    def test_missing(self):
        g = c.grade({})
        self.assertTrue(g["V0"].startswith("FAIL"))
        self.assertTrue(g["V1"].startswith("UNDETERMINED"))


class TestMechanism(unittest.TestCase):
    """Interior unsampled bins inflate the model conductivity; cover bands conserve it."""

    def test_cover_conserves(self):
        try:
            import numpy as np
            import tbc_bte as tb
        except ImportError:
            self.skipTest("numpy or tbc_bte not present")
        Ka = np.array([1e6, 1e6, 1e6, 1e6])
        ka = np.array([0.0, 5.0, 0.0, 5.0])
        # interpolated lengths: every bin 5 / 2e6 = 2.5e-6 m -> model kappa 2 x 4e6 x 2.5e-6 = 20
        la, nf = tb.lengths_from_kappa(ka, Ka)
        self.assertEqual(nf, 2)
        self.assertAlmostEqual(float(np.sum(2 * Ka * la)), 20.0, places=9)
        Kb = [1e6] * 4
        kb = [1.0] * 4
        bands = c.cover_bands(list(Ka), list(ka), Kb, kb)
        self.assertEqual(bands, [(0, 1), (2, 3)])
        Ka2, ka2 = c.banded(Ka, ka, bands)
        Kb2, kb2 = c.banded(Kb, kb, bands)
        la2, nf2 = tb.lengths_from_kappa(np.array(ka2), np.array(Ka2))
        lb2, _ = tb.lengths_from_kappa(np.array(kb2), np.array(Kb2))
        self.assertEqual(nf2, 0)
        s = tb.solve_two_stream(np.array(Ka2), np.array(Kb2), la2, lb2, n_nodes=80)
        # bands K 2e6, kappa 5 -> l 1.25e-6 -> 2 x 2e6 x 1.25e-6 x 2 = 10, the input
        self.assertAlmostEqual(s["kappa1"], 10.0, places=9)
        self.assertAlmostEqual(s["kappa2"], 4.0, places=9)
        # DMM on the bands: 2 x (2e6 x 2e6 / 4e6) = 2e6
        self.assertAlmostEqual(s["G_dmm"], 2e6, places=3)


if __name__ == "__main__":
    unittest.main()
