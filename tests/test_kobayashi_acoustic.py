import unittest

import numpy as np

from benchmarks.kobayashi_2026.acoustic_finite_size import plane_wave_rayleigh
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

    def test_cumulative_trapezoid(self):
        t = np.array([0.0, 1.0, 2.0])
        y = np.array([2.0, 2.0, 2.0])
        got = cumulative_trapezoid(t, y)
        self.assertTrue(np.allclose(got, [0.0, 2.0, 4.0]))


if __name__ == "__main__":
    unittest.main()
