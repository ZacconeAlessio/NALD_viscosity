import unittest

import numpy as np

from benchmarks.kobayashi_2026.static_nald import (
    affine_force,
    dense_mass_weighted_hessian,
    forces_and_pair_metrics,
    source_box_length_m,
)


class TestKobayashiStaticNALD(unittest.TestCase):
    def setUp(self):
        self.coords = np.array(
            [
                [0.20e-9, 0.20e-9, 0.20e-9],
                [0.54e-9, 0.20e-9, 0.20e-9],
                [0.37e-9, 0.49e-9, 0.31e-9],
            ],
            dtype=float,
        )
        self.cell = 2.0e-9

    def test_source_box_length(self):
        self.assertAlmostEqual(source_box_length_m() * 1e9, 4.3250816376, places=9)

    def test_newton_third_law(self):
        forces, pair_count, _ = forces_and_pair_metrics(self.coords, self.cell, 1.0)
        self.assertEqual(pair_count, 3)
        scale = max(float(np.linalg.norm(forces)), 1.0e-30)
        self.assertLess(float(np.linalg.norm(forces.sum(axis=0))), 1.0e-12 * scale)

    def test_hessian_symmetry_and_translations(self):
        h = dense_mass_weighted_hessian(self.coords, self.cell, 1.0)
        self.assertTrue(np.allclose(h, h.T, rtol=1.0e-13, atol=0.0))

        hnorm = max(float(np.linalg.norm(h)), 1.0)
        for axis in range(3):
            translation = np.zeros(3 * len(self.coords))
            translation[axis::3] = 1.0
            self.assertLess(float(np.linalg.norm(h @ translation)), 1.0e-12 * hnorm)

    def test_affine_force_has_zero_total_force(self):
        for channel in ("bulk", "shear"):
            xi = affine_force(self.coords, self.cell, 1.0, channel, delta=1.0e-6)
            scale = max(float(np.linalg.norm(xi)), 1.0e-30)
            self.assertLess(float(np.linalg.norm(xi.sum(axis=0))), 1.0e-9 * scale)


if __name__ == "__main__":
    unittest.main()
