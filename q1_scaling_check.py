# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""Pre-registered check Q1 (PLAN-amendment-0.5.5.md), written before it was ever run.

Claim: in the multiband two-stream model (tbc_bte.solve_two_stream) the interface conductance G
does not change when the lengths of ONE side only are multiplied by a constant. Reason: each
half-space is coupled to the other only at x = 0; x -> x / lambda on that side maps its solution
onto itself, and its bulk conductivity, its stream offsets l Q / kappa and its extrapolated bulk
profile at x = 0 are all unchanged. If true, G depends only on the spectral SHAPE of each side's
lengths, so a potential that gets every lifetime of a material wrong by the same factor still
gives the right G. Expected PASS (90 %). Run in job tbc-anhpairs-v1 as information (it does not
stop the job); graded in amendment 0.5.6. Not named test_*.py on purpose: it is a pre-registered
check whose outcome is unknown, not a unit test that guards a commit.
Run: python -m unittest -v q1_scaling_check
"""
import unittest

try:
    import numpy as np
    import tbc_bte as tb
except ImportError:
    np = None


@unittest.skipIf(np is None, "numpy or tbc_bte not available")
class TestOneSideScaling(unittest.TestCase):
    K1 = [1e8, 3e8, 2e8, 5e7]
    K2 = [2e8, 1e8, 4e8, 3e8]
    L1 = [1e-8, 3e-7, 2e-6, 5e-9]
    L2 = [5e-8, 1e-7, 1e-6, 2e-8]

    def g(self, f1, f2, n=300):
        l1 = np.array(self.L1) * f1
        l2 = np.array(self.L2) * f2
        return tb.solve_two_stream(self.K1, self.K2, l1, l2, n_nodes=n)["G"]

    def test_scale_side_one(self):
        g0 = self.g(1.0, 1.0)
        for f in (0.2, 7.0):
            self.assertLess(abs(self.g(f, 1.0) / g0 - 1.0), 0.01, f"side 1 x{f}")

    def test_scale_side_two(self):
        g0 = self.g(1.0, 1.0)
        for f in (0.1, 5.0):
            self.assertLess(abs(self.g(1.0, f) / g0 - 1.0), 0.01, f"side 2 x{f}")

    def test_shape_matters(self):
        # control: changing the relative lengths within one side does change G (else the model
        # would be gray and the check above empty); expected change above 1 %
        g0 = self.g(1.0, 1.0)
        l1 = np.array(self.L1) * np.array([30.0, 1.0, 1.0, 30.0])
        g1 = tb.solve_two_stream(self.K1, self.K2, l1, np.array(self.L2), n_nodes=300)["G"]
        self.assertGreater(abs(g1 / g0 - 1.0), 0.01)


if __name__ == "__main__":
    unittest.main()
