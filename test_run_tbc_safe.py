# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""Tests for run_tbc_safe (no torch needed: the phonon child is replaced by a fake)."""
import multiprocessing as mp
import os
import shutil
import signal
import tempfile
import time
import unittest

import numpy as np

import run_tbc_safe as rs

INFO = dict(relax_converged=True, cell_lengths=[4.0, 4.0, 4.0], n_atoms_per_m3=6.0e28, max_thz=21.10,
            min_thz=0.1, imag_fraction=0.0, seconds=0.5, n_disp=1, n_prim=4)


def fake_child(model, name, dtype, edges, normal, stem):
    """Stands in for rs.child: 'C' is killed in float64 (as by the OOM killer), 'X' always exits 3."""
    if name == "C" and dtype == "float64":
        os.kill(os.getpid(), signal.SIGKILL)
    if name == "X":
        raise SystemExit(3)
    n = len(edges) - 1
    arrays = {"iso": np.full(n, 1.0)}
    if normal is not None:
        arrays["or"] = np.full(n, 2.0)
    info = dict(INFO)
    info["dtype"] = dtype
    info["peak_rss_gb"] = 0.05
    rs.write_part(stem, arrays, info)


class TestSmall(unittest.TestCase):
    def test_describe_exit(self):
        self.assertEqual(rs.describe_exit(0), "ok")
        self.assertIn("out of memory", rs.describe_exit(-9))
        self.assertEqual(rs.describe_exit(-15), "signal 15")
        self.assertEqual(rs.describe_exit(3), "exit 3")
        self.assertEqual(rs.describe_exit("timeout"), "timeout")

    def test_should_retry(self):
        self.assertTrue(rs.should_retry(-9))
        self.assertFalse(rs.should_retry(3))
        self.assertFalse(rs.should_retry(0))
        self.assertFalse(rs.should_retry("timeout"))
        self.assertFalse(rs.should_retry(None))

    def test_attempts(self):
        self.assertEqual(rs.attempts_for("mace-mpa-0"), ["float64", "float32"])
        self.assertEqual(rs.attempts_for("orb-v3"), ["default"])

    def test_finish_info_gan(self):
        # 21.10 / 16.8 = 1.25595 (pre-registered E1(TO) reference); 21.10 / 20.67 = 1.02080 (B1 high)
        d = rs.finish_info(dict(INFO), "GaN", 16.8)
        self.assertAlmostEqual(d["max_ratio"], 1.25595, places=4)
        self.assertAlmostEqual(d["max_ratio_b1high_info"], 1.02080, places=4)
        d2 = rs.finish_info(dict(INFO), "Ti", None)
        self.assertIsNone(d2["max_ratio"])
        self.assertNotIn("max_ratio_b1high_info", d2)

    def test_part_roundtrip(self):
        tmp = tempfile.mkdtemp()
        try:
            stem = os.path.join(tmp, "m-Al-float64")
            rs.write_part(stem, {"iso": np.arange(3.0)}, {"a": 1})
            arrays, info = rs.read_part(stem)
            self.assertEqual(info, {"a": 1})
            self.assertEqual(list(arrays["iso"]), [0.0, 1.0, 2.0])
            self.assertFalse(os.path.exists(stem + ".json.tmp"))
        finally:
            shutil.rmtree(tmp)

    def test_safe_line(self):
        ibm = {"mace-mpa-0": {
            "Al": {"dtype": "float64", "peak_rss_gb": 1.2},
            "GaN": {"dtype": "float32", "peak_rss_gb": 2.5, "max_ratio_b1high_info": 1.0208},
            "C": {"failed": "float64:killed; float32:killed"}}}
        self.assertEqual(rs.safe_line(ibm), "SAFE: mace-mpa-0 done 2/3 f32 [GaN] failed [C] peak 2.50 GB "
                                            "GaN/B1high 1.021")


@unittest.skipUnless(hasattr(signal, "SIGKILL") and os.name == "posix", "needs POSIX signals")
class TestOrchestration(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.ctx = mp.get_context("spawn")
        self.edges = np.linspace(0.0, 1e12, 5)
        self.lines = []

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_kill_then_float32_and_no_retry_after_exit(self):
        spectra, info = rs.run_model_safe("mace-mpa-0", self.edges, self.tmp, time.time() + 300,
                                          self.lines.append, self.ctx, ["Al", "C", "X"],
                                          {"Al": [1, 1, 1]}, {"Al": 9.7}, target=fake_child)
        self.assertEqual(sorted(spectra), ["Al", "C"])
        self.assertIn("or", spectra["Al"])
        self.assertNotIn("or", spectra["C"])
        self.assertEqual(len(spectra["Al"]["iso"]), 4)
        self.assertEqual(info["Al"]["dtype"], "float64")
        self.assertEqual(info["C"]["dtype"], "float32")
        self.assertTrue(info["C"]["attempts"][0].startswith("float64:killed"))
        self.assertEqual(info["C"]["attempts"][1], "float32:ok")
        self.assertEqual(info["X"], {"failed": "float64:exit 3"})
        self.assertAlmostEqual(info["Al"]["max_ratio"], 21.10 / 9.7, places=9)
        self.assertTrue(os.path.exists(os.path.join(self.tmp, "spectra-mace-mpa-0.npz")))
        self.assertTrue(any(s.startswith("mace-mpa-0 X: FAILED") for s in self.lines))

    def test_orb_one_attempt(self):
        spectra, info = rs.run_model_safe("orb-v3", self.edges, self.tmp, time.time() + 300,
                                          self.lines.append, self.ctx, ["X", "Al"], {}, {},
                                          target=fake_child)
        self.assertEqual(info["X"], {"failed": "default:exit 3"})
        self.assertEqual(info["Al"]["dtype"], "default")
        self.assertEqual(list(spectra), ["Al"])

    def test_no_time_left(self):
        spectra, info = rs.run_model_safe("orb-v3", self.edges, self.tmp, time.time() + 5,
                                          self.lines.append, self.ctx, ["Al"], {}, {},
                                          target=fake_child)
        self.assertEqual(info["Al"], {"failed": "default:no time left"})
        self.assertEqual(spectra, {})


if __name__ == "__main__":
    unittest.main()
