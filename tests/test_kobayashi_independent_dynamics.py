import unittest

import numpy as np

from benchmarks.kobayashi_2026.green_kubo_spectrum import (
    complex_viscosity_from_integrand,
)
from benchmarks.kobayashi_2026.vacf_memory import (
    vacf_memory_laplace,
)


class TestKobayashiIndependentDynamics(unittest.TestCase):
    def test_exponential_vacf_recovers_markovian_memory(self):
        nu = 2.5
        t = np.linspace(0.0, 20.0, 200001)
        c = np.exp(-nu * t)
        for s in (0.0, 0.3, 1.2, 2.0 + 0.7j):
            got = vacf_memory_laplace(t, c, s)
            self.assertAlmostEqual(got.real, nu, places=5)
            self.assertAlmostEqual(got.imag, 0.0, places=5)

    def test_green_kubo_zero_frequency_is_time_integral(self):
        t = np.linspace(0.0, 10.0, 100001)
        c = 3.0 * np.exp(-2.0 * t)
        got = complex_viscosity_from_integrand(t, c, np.array([0.0]))
        self.assertAlmostEqual(got[0].real, 1.5, places=7)
        self.assertAlmostEqual(got[0].imag, 0.0, places=14)


if __name__ == "__main__":
    unittest.main()
