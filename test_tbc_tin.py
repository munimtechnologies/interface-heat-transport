# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""Unit tests of tbc_tin. Every number below was worked out by hand first
(PLAN-amendment-0.5.2.md); the tests recompute it. Pure Python: no stage 2 modules needed."""
import json
import os
import unittest

import tbc_tin as tt

HERE = os.path.dirname(os.path.abspath(__file__))


class TestReading(unittest.TestCase):
    def test_constants_match_data_file(self):
        d = json.load(open(os.path.join(HERE, "data", "read-hight.json")))
        for pair, row in tt.KHAN.items():
            got = {int(k): float(v) for k, v in d["khan2024"]["G_at"][pair].items()
                   if v is not None}
            self.assertEqual(got, row, pair)
            T, g = d["khan2024"]["last_point"][pair]
            self.assertEqual((float(T), float(g)), tt.KHAN_LAST[pair])
        for pair, c in tt.KHAN_CURVES.items():
            for name, v in c.items():
                self.assertEqual(tuple(float(x) for x in d["khan2024"]["curves_300_1000"][pair][name]),
                                 v, (pair, name))
        for pair, v in tt.KOH.items():
            self.assertEqual(tuple(float(x) for x in d["koh2020"]["G_300_400_500"][pair]), v)
            self.assertEqual(tuple(float(x) for x in d["koh2020"]["nel_dmm_300_500"][pair]),
                             tt.KOH_NEL[pair])


class TestPieces(unittest.TestCase):
    def test_cap(self):
        self.assertAlmostEqual(tt.cap(1e-9, 300.0), 1.0, places=9)
        self.assertAlmostEqual(tt.cap(18.0, 300.0), 0.522727, delta=2e-5)
        self.assertAlmostEqual(tt.cap(18.0, 996.0), 0.939601, delta=2e-5)
        self.assertAlmostEqual(tt.cap(15.6, 300.0), 0.609908, delta=2e-5)
        self.assertAlmostEqual(tt.cap(9.1, 771.0), 0.973679, delta=2e-5)

    def test_dressed_gray(self):
        # equal sides, perfect transmission: the interface disappears
        self.assertEqual(tt.dressed_gray(1.0, 1.0, 1.0), float("inf"))
        # G0 = 0.5, R1 = R2 = 1 (gray DMM): G_d = 1
        self.assertAlmostEqual(tt.dressed_gray(0.5, 1.0, 1.0), 1.0, places=12)
        self.assertAlmostEqual(tt.dressed_gray(364.0, 1290.0, 3200.0), 453.843, delta=0.01)

    def test_floors(self):
        # E = 400 -> 1.2 x 400, I = 0.3 T: G1 = 490, G2(3.32x) = 480 + 298.8
        g1, g2 = 400 + 90.0, 480 + 0.3 * 996
        self.assertAlmostEqual(tt.channel_floor_linear(g1, g2, 3.32, 1.2), 90.0, places=9)
        self.assertAlmostEqual(tt.floor_nonelastic(g1, g2, 1.2), 0.3 * 996 - 1.2 * 90, places=9)

    def test_t_eff_001(self):
        self.assertAlmostEqual(tt.t_eff_001("rocksalt001", 4.2417e-10), 4.2417e-10 / 12, places=20)
        self.assertAlmostEqual(tt.t_eff_001("diamond001", 5.431e-10), 5.431e-10 / 8, places=20)

    def test_interp(self):
        self.assertEqual(tt.interp(tt.KHAN["TiN/C"], 300), 449.0)
        self.assertAlmostEqual(tt.interp(tt.KHAN["TiN/C"], 350), 510.5, places=9)


class TestNumbers(unittest.TestCase):
    hand = {
        "TiN/C dressed elastic_limit 300": (453.843, 0.01),
        "TiN/C dressed elastic_limit 1000": (629.685, 0.01),
        "TiN/C dressed DMM 300": (301.463, 0.01),
        "TiN/C dressed DMM 1000": (408.362, 0.01),
        "TiN/Si dressed elastic_limit 300": (635.736, 0.01),
        "TiN/Si dressed elastic_limit 1000": (647.976, 0.01),
        "TiN/Si dressed DMM 300": (388.152, 0.01),
        "TiN/Si dressed DMM 1000": (403.943, 0.01),
        "TiN/C kappa dressed EL 300-1000": (1.38745, 1e-4),
        "TiN/Si kappa dressed EL 300-1000": (1.019253, 1e-5),
        "TiN/C excess over dressed EL at 996": (318.315, 0.01),
        "TiN/C floor nonelastic 996 (dressed EL trend)": (325.035, 0.05),
        "TiN/Si floor nonelastic 771 (dressed EL trend, 300-1000 ratio)": (425.239, 0.01),
        "TiN/C linear channel floor I300 (dressed EL trend)": (168.19, 0.05),
        "TiN/C kappa conv max (TiN 18.0 THz) 300-996": (1.797499, 2e-4),
        "TiN/C floor nonelastic 996 (model-free conventional)": (140.92, 0.1),
        "TiN/Si kappa conv max (Si 15.6 THz) 300-771": (1.516610, 2e-4),
        "TiN/Si floor nonelastic 771 (model-free conventional)": (173.08, 0.1),
        "TiN/Ge kappa conv max (Ge 9.1 THz) 300-771": (1.158255, 2e-4),
        "TiN/Ge floor nonelastic 771 (model-free conventional)": (78.93, 0.05),
        "TiN/C measured ratio 900/300": (2.026726, 1e-5),
        "TiN/Si measured ratio 700/300": (1.806706, 1e-5),
        "TiN/Ge measured ratio 700/300": (1.478947, 1e-5),
    }

    def test_hand_numbers(self):
        x = tt.numbers()
        for k, (v, tol) in self.hand.items():
            self.assertIn(k, x)
            self.assertAlmostEqual(x[k], v, delta=tol, msg=k)

    def test_results_file_matches(self):
        r = json.load(open(os.path.join(HERE, "results", "tin", "bounds-0.5.2.json")))
        x = tt.numbers()
        for k, (v, tol) in self.hand.items():
            self.assertAlmostEqual(r["numbers"][k], x[k], delta=tol, msg=k)

    def test_worst_reading_floors_stay_positive(self):
        # model-free floors with every reading moved by its error against the claim
        for pair, side, T2 in (("TiN/C", "TiN", 996.0), ("TiN/Si", "Si", 771.0), ("TiN/Ge", "Ge", 771.0)):
            e = tt.READ_ERR[pair]
            k = tt.cap_ratio(tt.NU_MAX_REF[side], 300.0, T2)
            f = tt.floor_nonelastic(tt.KHAN[pair][300] + e, tt.KHAN_LAST[pair][1] - e, k)
            self.assertGreater(f, 40.0, pair)


class TestGrades(unittest.TestCase):
    def prim(self, scale_hi=1.0, bte_r=1.2, tot300=None):
        tot300 = tot300 or {"TiN/C": 449.0, "TiN/Si": 507.0, "TiN/Ge": 190.0}
        mats = {m: {"max_thz": tt.NU_MAX_REF[m], "imag_fraction": 0.0} for m in tt.NU_MAX_REF}
        pairs = {}
        for p in tt.PAIRS:
            t300 = tot300[p]
            thi = tt.KHAN[p][tt.T_HIGH[p]] * scale_hi
            pairs[p] = {"300": {"G_bte": 300.0, "G_total": t300, "G_dmm": 250.0},
                        str(tt.T_HIGH[p]): {"G_bte": 300.0 * bte_r, "G_total": thi, "G_dmm": 280.0}}
        return {"materials": mats, "pairs": pairs}

    def test_all_pass(self):
        g = tt.grade_t(self.prim())
        for k in ("T0", "T1", "T2", "T3", "T4", "T5", "T6"):
            self.assertTrue(g[k].startswith("PASS"), (k, g[k]))
        self.assertIn("0.375", g["T6"])          # 190 / 507

    def test_fails(self):
        g = tt.grade_t(self.prim(scale_hi=0.6, bte_r=2.5,
                                 tot300={"TiN/C": 600.0, "TiN/Si": 500.0, "TiN/Ge": 400.0}))
        self.assertTrue(g["T1"].startswith("FAIL"))   # elastic rises 2.5x, more than the data
        self.assertTrue(g["T2"].startswith("FAIL"))
        self.assertTrue(g["T3"].startswith("FAIL"))
        self.assertTrue(g["T5"].startswith("FAIL"))   # C > Si
        self.assertTrue(g["T6"].startswith("FAIL"))   # 0.8

    def test_t0_and_missing(self):
        p = self.prim()
        p["materials"]["TiN"]["max_thz"] = 30.0
        self.assertTrue(tt.grade_t(p)["T0"].startswith("FAIL"))
        p = self.prim()
        del p["pairs"]["TiN/Ge"]
        g = tt.grade_t(p)
        self.assertTrue(g["T0"].startswith("FAIL"))
        self.assertTrue(g["T3"].startswith("UNDETERMINED"))


if __name__ == "__main__":
    unittest.main()
