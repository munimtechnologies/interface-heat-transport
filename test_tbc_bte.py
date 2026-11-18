# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""Unit tests of tbc_bte (numpy only)."""
import unittest

import numpy as np

import tbc_bte as tb


def simons(G, K1, K2):
    return 1.0 / (1.0 / G - 0.5 / K1 - 0.5 / K2)


class TestBTE(unittest.TestCase):
    def test_dmm_alpha_detailed_balance(self):
        K1 = np.array([1.0, 2.0, 0.0, 3.0])
        K2 = np.array([4.0, 0.0, 1.0, 3.0])
        a1, a2, g = tb.dmm_alpha(K1, K2)
        np.testing.assert_allclose(K1 * a1, K2 * a2)
        np.testing.assert_allclose(g, [0.8, 0.0, 0.0, 1.5])

    def test_gray_equals_simons(self):
        for K1, K2, l1, l2 in ((3e8, 5e8, 2e-9, 7e-8), (1e9, 1e9, 1e-7, 1e-7), (2e8, 9e8, 5e-6, 3e-10)):
            r = tb.solve_two_stream([K1], [K2], [l1], [l2])
            G = K1 * K2 / (K1 + K2)
            self.assertAlmostEqual(r["G"] / simons(G, K1, K2), 1.0, places=7)
            self.assertAlmostEqual(r["G_simons"] / simons(G, K1, K2), 1.0, places=9)
            self.assertAlmostEqual(r["G_linear_bath"] / simons(G, K1, K2), 1.0, places=9)
            self.assertLess(r["residual"], 1e-8)

    def test_equal_lengths_reduce_to_gray(self):
        K = np.array([1e8, 3e8, 2e8, 5e7])
        ell = np.full(4, 4e-8)
        r = tb.solve_two_stream(K, K, ell, ell)
        self.assertAlmostEqual(r["G"] / float(np.sum(K)), 1.0, places=6)

    def test_kappa_matches_lengths(self):
        K1 = np.array([1e8, 2e8])
        l1 = np.array([1e-8, 3e-7])
        r = tb.solve_two_stream(K1, K1, l1, l1)
        self.assertAlmostEqual(r["kappa1"], 2 * (1e8 * 1e-8 + 2e8 * 3e-7), places=6)

    def _multiband(self):
        rng = np.random.default_rng(7)
        B = 30
        K1 = rng.uniform(1e6, 3e7, B)
        K2 = rng.uniform(1e6, 3e7, B)
        K1[-6:] = 0.0          # bands that exist only on side 2
        K2[:3] = 0.0           # and only on side 1
        l1 = np.exp(rng.uniform(np.log(1e-9), np.log(1e-6), B))
        l2 = np.exp(rng.uniform(np.log(1e-9), np.log(3e-6), B))
        return K1, K2, l1, l2

    def test_swap_symmetry(self):
        K1, K2, l1, l2 = self._multiband()
        g12 = tb.solve_two_stream(K1, K2, l1, l2)["G"]
        g21 = tb.solve_two_stream(K2, K1, l2, l1)["G"]
        self.assertAlmostEqual(g12 / g21, 1.0, places=6)
        self.assertTrue(np.isfinite(g12) and g12 > 0.0)

    def test_scale_invariance(self):
        # multiplying every length by the same factor leaves G unchanged (stage 2a, S1)
        K1, K2, l1, l2 = self._multiband()
        g = tb.solve_two_stream(K1, K2, l1, l2, n_nodes=240)["G"]
        for lam in (0.1, 3.0, 20.0):
            gl = tb.solve_two_stream(K1, K2, lam * l1, lam * l2, n_nodes=240)["G"]
            self.assertAlmostEqual(gl / g, 1.0, delta=1e-2)

    def test_grid_convergence(self):
        K1, K2, l1, l2 = self._multiband()
        g1 = tb.solve_two_stream(K1, K2, l1, l2, n_nodes=160)["G"]
        g2 = tb.solve_two_stream(K1, K2, l1, l2, n_nodes=400)["G"]
        self.assertAlmostEqual(g1 / g2, 1.0, delta=1e-2)

    def test_non_gray_differs_from_gray(self):
        # disparate lengths: the spectral answer is not the gray Simons value, and is finite
        K1, K2, l1, l2 = self._multiband()
        r = tb.solve_two_stream(K1, K2, l1, l2)
        self.assertGreater(abs(r["G"] / r["G_simons"] - 1.0), 1e-3)
        self.assertGreater(r["G"], 0.0)

    def test_lengths_from_kappa_fill(self):
        kb = np.array([0.0, 2.0, 0.0, 8.0, 0.0])
        Kb = np.array([0.0, 1.0, 1.0, 1.0, 0.0])
        ell, nfill = tb.lengths_from_kappa(kb, Kb)
        self.assertEqual(nfill, 1)
        self.assertAlmostEqual(ell[1], 1.0)
        self.assertAlmostEqual(ell[3], 4.0)
        self.assertAlmostEqual(ell[2], 2.0)   # geometric mean of neighbours
        self.assertEqual(ell[0], 0.0)

    def test_rebin(self):
        ef = np.linspace(0, 10, 11)
        ec = np.array([0.0, 5.0, 10.0])
        np.testing.assert_allclose(tb.rebin(np.ones(10), ef, ec), [5.0, 5.0])


if __name__ == "__main__":
    unittest.main()
