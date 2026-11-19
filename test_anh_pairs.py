# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""Tests for anh_pairs.py and run_anh_safe.py. Skipped where numpy or the repository modules
(run_tbc_anh, tbc_anh, tbc_bte, tbc_transport) are not importable."""
import os
import tempfile
import unittest

try:
    import numpy as np
    import anh_pairs as ap
    import run_anh_safe as ras
    import tbc_anh as ta
    import tbc_transport as tt
except ImportError:  # numpy or a repository module is absent
    np = None


@unittest.skipIf(np is None, "numpy or repository modules not available")
class TestAnhPairs(unittest.TestCase):
    def test_prim_volume_cubic(self):
        # diamond / fcc: a^3 / 4; Si a = 5.431 A -> 160.19148 / 4 = 40.04787 A^3
        self.assertAlmostEqual(ap.prim_volume("Si", (5.431,) * 3) * 1e30, 40.0479, places=3)

    def test_prim_volume_hex(self):
        # wurtzite: (sqrt 3 / 2) a^2 c; GaN 3.189, 5.185 -> 45.66552 A^3
        self.assertAlmostEqual(ap.prim_volume("GaN", (3.189, 3.189, 5.185)) * 1e30, 45.6655,
                               places=3)

    def test_kappa_bins_one_classical_mode(self):
        # one mode at 5 THz, gamma 0.1 THz -> tau = 1/(4 pi 1e11) s; kB * 1e6 m^2/s^2 * tau / 1e-29 m^3
        r = self._one_mode()
        edges_c = tt.make_edges(48.0, 0.2)
        k = ta.kappa_bins(r, 0, 100.0, edges_c, "iso", with_isotope=False, classical=True)
        self.assertAlmostEqual(float(k.sum()), 1.0986856, delta=1e-4)

    def test_material_entry_matches_direct_formula(self):
        r = self._one_mode()
        edges_c = tt.make_edges(48.0, 0.2)
        d = ap.material_entry(r, "Si", edges_c)
        self.assertEqual(sorted(d["kappa"].keys(), key=int), ["100", "200", "300", "400", "500",
                                                                 "600", "700"])
        tau = 1.0 / (4.0 * np.pi * 0.1e12)
        for T in (100.0, 300.0, 700.0):
            c = float(tt.mode_heat_capacity(np.array([5e12]), T)[0])
            want = c * 1e6 * tau / 1e-29
            got = d["kappa"][str(int(T))]["iso"]
            self.assertAlmostEqual(got["sum_bins_3ph"] / want, 1.0, places=9)
            self.assertAlmostEqual(got["sum_bins_iso"] / want, 1.0, places=9)  # no isotope rate

    def test_load_material_volume_from_v1_cells(self):
        r = self._one_mode()
        del r["vcell"]
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "anh-orb-v3-GaN.npz")
            np.savez(path, **r)
            got, src = ap.load_material(path, "orb-v3", "GaN")
        self.assertEqual(src, "printed cell lengths")
        want = 0.5 * np.sqrt(3.0) * 3.2180 * 3.2180 * 5.2454 * 1e-30
        self.assertAlmostEqual(got["vcell"] / want, 1.0, places=12)

    def test_parse_sources(self):
        self.assertEqual(ap.parse_sources(["orb-v3=/a", "mace-mpa-0=/b/c"]),
                         [("orb-v3", "/a"), ("mace-mpa-0", "/b/c")])
        with self.assertRaises(ValueError):
            ap.parse_sources(["orb-v3"])

    def test_safe_driver_retry_rules(self):
        self.assertEqual(ras.attempts_for("mace-mpa-0"), ["float64", "float32"])
        self.assertEqual(ras.attempts_for("orb-v3"), ["default"])
        self.assertTrue(ras.should_retry(-9))
        self.assertFalse(ras.should_retry(1))
        self.assertFalse(ras.should_retry("timeout"))
        self.assertEqual(ras.describe_exit(-9), "killed (SIGKILL, likely out of memory)")

    @staticmethod
    def _one_mode():
        temps = np.array([100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0])
        edges_f = tt.make_edges(48.0, 0.02)
        return {"freq": np.array([[5.0]]), "weights": np.array([1.0]),
                "gamma3": np.full((7, 1, 1), 0.1), "gamma_iso": np.zeros((1, 1)),
                "v2_iso": np.array([[1e6]]), "v2_z": np.array([[1e6]]),
                "kappa_ph3": np.zeros((7, 6)), "temps": temps, "vcell": 1e-29,
                "edges_c": tt.make_edges(48.0, 0.2), "edges_f": edges_f,
                "M_or": np.zeros(edges_f.size - 1), "M_iso": np.zeros(edges_f.size - 1)}


if __name__ == "__main__":
    unittest.main()
