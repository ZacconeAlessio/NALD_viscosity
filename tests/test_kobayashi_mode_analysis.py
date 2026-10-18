import unittest

import numpy as np

from benchmarks.kobayashi_2026.mode_analysis import (
    branch_counts,
    null_gamma_fraction,
    signed_frequency_thz,
    structural_prefactor,
)


class TestKobayashiModeAnalysis(unittest.TestCase):
    def test_signed_frequency(self):
        w = np.array([-4.0, 0.0, 9.0]) * (2.0 * np.pi * 1e12) ** 2
        f = signed_frequency_thz(w)
        self.assertTrue(np.allclose(f, [-2.0, 0.0, 3.0]))

    def test_structural_prefactor_keeps_negative_modes(self):
        lam = np.array([-4.0e24, 1.0e-15, 9.0e24])
        gamma = np.array([2.0, 100.0, 3.0])
        got = structural_prefactor(
            lam,
            gamma,
            volume_m3=2.0,
            cutoff_thz=0.0,
            null_lambda_s2=1.0e12,
        )
        expected = 0.5 * (2.0 / (4.0e24) ** 2 + 3.0 / (9.0e24) ** 2)
        self.assertAlmostEqual(got / expected, 1.0, places=13)

    def test_cutoff_is_symmetric_in_frequency(self):
        unit = (2.0 * np.pi * 1e12) ** 2
        lam = np.array([-0.01 * unit, -0.04 * unit, 0.01 * unit, 0.04 * unit])
        gamma = np.ones(4)
        got = structural_prefactor(
            lam,
            gamma,
            volume_m3=1.0,
            cutoff_thz=0.15,
            null_lambda_s2=0.0,
        )
        expected = 2.0 / (0.04 * unit) ** 2
        self.assertAlmostEqual(got / expected, 1.0, places=13)

    def test_branch_and_null_coupling(self):
        lam = np.array([-2.0e13, -2.0, 0.0, 3.0, 5.0e13])
        gamma = np.array([1.0, 0.0, 0.0, 0.0, 3.0])
        self.assertEqual(
            branch_counts(lam),
            {"negative": 1, "positive": 1, "null": 3},
        )
        self.assertEqual(null_gamma_fraction(lam, gamma), 0.0)


if __name__ == "__main__":
    unittest.main()
