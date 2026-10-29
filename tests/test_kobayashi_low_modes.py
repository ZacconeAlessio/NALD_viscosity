import unittest

import numpy as np

from benchmarks.kobayashi_2026.low_mode_character import (
    participation_ratio,
    qmin_plane_wave_weights,
    summarize_subcutoff_modes,
)


class TestKobayashiLowModeCharacter(unittest.TestCase):
    def test_participation_ratio_limits(self):
        # Two atoms, one Cartesian degree per atom used.
        localized = np.zeros((6, 1))
        localized[0, 0] = 1.0
        self.assertAlmostEqual(
            participation_ratio(localized, 2)[0], 0.5, places=14
        )

        extended = np.zeros((6, 1))
        extended[0, 0] = 1.0 / np.sqrt(2.0)
        extended[3, 0] = 1.0 / np.sqrt(2.0)
        self.assertAlmostEqual(
            participation_ratio(extended, 2)[0], 1.0, places=14
        )

    def test_qmin_longitudinal_plane_wave_overlap(self):
        length = 2.0
        coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.5, 0.7]])
        mode = np.zeros((6, 1))
        mode[0, 0] = 1.0 / np.sqrt(2.0)
        mode[3, 0] = -1.0 / np.sqrt(2.0)

        _, longitudinal, _ = qmin_plane_wave_weights(coords, mode, length)
        self.assertAlmostEqual(longitudinal[0], 1.0, places=14)

    def test_subcutoff_summary_is_contribution_weighted(self):
        unit = (2.0 * np.pi * 1.0e12) ** 2
        lam = np.array([0.01 * unit, 0.04 * unit, 1.0e-15])
        gamma = np.array([4.0, 1.0, 100.0])
        pr = np.array([0.2, 0.8, 1.0])
        qw = np.array([0.1, 0.9, 1.0])

        out = summarize_subcutoff_modes(
            lam,
            gamma,
            pr,
            qw,
            cutoff_thz=0.15,
            volume_m3=1.0,
        )
        self.assertEqual(out["subcutoff_modes"], 1)
        self.assertAlmostEqual(out["A_weighted_participation_ratio"], 0.2)
        self.assertAlmostEqual(out["A_weighted_qmin_plane_wave_weight"], 0.1)


if __name__ == "__main__":
    unittest.main()
