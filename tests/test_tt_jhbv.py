import unittest

from benchmarks.kobayashi_2026.tt_jhbv import (
    C_K_NM,
    RMIN_NM,
    SIGMA_NM,
    base_potential_K,
)


class TestTTJHBV(unittest.TestCase):
    def test_potential_minimum_value(self):
        u, du, _ = base_potential_K(RMIN_NM)
        self.assertAlmostEqual(u, -143.1231887, places=5)
        self.assertLess(abs(du), 0.2)

    def test_sigma_is_near_zero_crossing(self):
        u, _, _ = base_potential_K(SIGMA_NM)
        self.assertLess(abs(u), 0.5)

    def test_corrigendum_c16(self):
        self.assertAlmostEqual(C_K_NM[16], 1.17006343e-6, places=15)


if __name__ == "__main__":
    unittest.main()
