import unittest

import numpy as np

from benchmarks.kobayashi_2026.acoustic_finite_size import (
    classify_projected_dynamical_matrix,
    plane_wave_rayleigh,
)
from benchmarks.kobayashi_2026.pressure_relaxation import cumulative_trapezoid


class TestKobayashiAcousticAndRelaxation(unittest.TestCase):
    def test_plane_wave_q_definition(self):
        coords = np.array([
            [0.20e-9, 0.20e-9, 0.20e-9],
            [0.54e-9, 0.20e-9, 0.20e-9],
            [0.37e-9, 0.49e-9, 0.31e-9],
        ])
        length = 2.0e-9
        out = plane_wave_rayleigh(coords, length, 1.0, harmonic_index=1)
        self.assertAlmostEqual(out["q_nm-1"], np.pi, places=12)
        self.assertTrue(np.isfinite(out["f_longitudinal_THz"]))
        self.assertTrue(np.isfinite(out["f_transverse_THz"]))

    def test_projected_matrix_accounts_for_polarization_mixing(self):
        rotation = np.array([
            [np.sqrt(0.8), np.sqrt(0.2), 0.0],
            [-np.sqrt(0.2), np.sqrt(0.8), 0.0],
            [0.0, 0.0, 1.0],
        ])
        eigenvalues = np.array([1.0, 4.0, 9.0]) * (2.0 * np.pi * 1.0e12) ** 2
        dyn = rotation @ np.diag(eigenvalues) @ rotation.T

        f_long, f_trans, weight = classify_projected_dynamical_matrix(dyn, 0)
        self.assertAlmostEqual(f_long, 1.0, places=12)
        np.testing.assert_allclose(np.sort(f_trans), [2.0, 3.0], rtol=1e-12)
        self.assertAlmostEqual(weight, 0.8, places=12)

    def test_cumulative_trapezoid(self):
        t = np.array([0.0, 1.0, 2.0])
        y = np.array([2.0, 2.0, 2.0])
        got = cumulative_trapezoid(t, y)
        self.assertTrue(np.allclose(got, [0.0, 2.0, 4.0]))


if __name__ == "__main__":
    unittest.main()
