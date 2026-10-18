import unittest

import numpy as np

from benchmarks.kobayashi_2026.pressure_memory import (
    estimate_omega0_squared,
    second_order_memory_laplace,
)


class TestPressureMemory(unittest.TestCase):
    def test_even_curvature_fit(self):
        omega0 = 4.0e12
        t = np.linspace(0.0, 0.08e-12, 81)
        c = 1.0 - 0.5 * omega0**2 * t**2 + 0.03 * omega0**4 * t**4
        got = estimate_omega0_squared(t, c, window_s=0.04e-12)
        self.assertAlmostEqual(got / omega0**2, 1.0, places=12)

    def test_markovian_oscillator_recovers_constant_memory(self):
        omega0 = 4.0
        nu = 1.5
        wd = np.sqrt(omega0**2 - (nu / 2.0) ** 2)
        t = np.linspace(0.0, 30.0, 60001)
        c = np.exp(-nu * t / 2.0) * (
            np.cos(wd * t) + (nu / (2.0 * wd)) * np.sin(wd * t)
        )
        for s in (0.0, 0.2, 1.0, 3.0):
            got = second_order_memory_laplace(
                t, c, s, omega0_sq_s2=omega0**2
            )
            self.assertAlmostEqual(got, nu, places=4)


if __name__ == "__main__":
    unittest.main()
