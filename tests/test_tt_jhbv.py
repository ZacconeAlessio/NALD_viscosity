import unittest

import numpy as np

from benchmarks.kobayashi_2026.tt_jhbv import (
    ARGON_MASS_KG,
    ATM_COEFFICIENT_J_M9,
    C_K_NM,
    FULL_PAIR_CUTOFF_M,
    MD_SIGMA_M,
    RMIN_M,
    RMIN_NM,
    base_potential_K,
    source_force_curvature_SI,
    source_pair_cutoff_m,
)


class TestTTJHBV(unittest.TestCase):
    def test_source_constants(self):
        self.assertAlmostEqual(ARGON_MASS_KG, 6.6335e-26, places=38)
        self.assertAlmostEqual(MD_SIGMA_M, 3.3952e-10, places=22)
        self.assertAlmostEqual(FULL_PAIR_CUTOFF_M, 1.6976e-9, places=21)
        self.assertAlmostEqual(RMIN_M, 3.7618e-10, places=22)
        self.assertAlmostEqual(ATM_COEFFICIENT_J_M9, 7.32e-108, places=120)

    def test_production_cutoff_switch(self):
        self.assertEqual(source_pair_cutoff_m(0.0), RMIN_M)
        self.assertEqual(source_pair_cutoff_m(1.0), FULL_PAIR_CUTOFF_M)
        self.assertEqual(source_pair_cutoff_m(0.5), FULL_PAIR_CUTOFF_M)

    def test_potential_minimum_value(self):
        u, du, _ = base_potential_K(RMIN_NM)
        self.assertAlmostEqual(u, -143.1232013, places=5)
        self.assertLess(abs(du), 2.0)

    def test_source_dispersion_recurrence(self):
        c12 = C_K_NM[6] * (C_K_NM[10] / C_K_NM[8]) ** 3
        c14 = C_K_NM[8] * (c12 / C_K_NM[10]) ** 3
        c16 = C_K_NM[10] * (c14 / c12) ** 3
        self.assertAlmostEqual(C_K_NM[12], c12, places=18)
        self.assertAlmostEqual(C_K_NM[14], c14, places=18)
        self.assertAlmostEqual(C_K_NM[16], c16, places=18)
        self.assertAlmostEqual(C_K_NM[16], 1.170063435858e-6, places=16)

    def test_s0_force_is_zero_beyond_rmin(self):
        radii = np.array([0.95 * RMIN_M, 1.05 * RMIN_M])
        du, d2u = source_force_curvature_SI(radii, 0.0)
        self.assertNotEqual(float(du[0]), 0.0)
        self.assertNotEqual(float(d2u[0]), 0.0)
        self.assertEqual(float(du[1]), 0.0)
        self.assertEqual(float(d2u[1]), 0.0)


if __name__ == "__main__":
    unittest.main()
