import unittest
import numpy as np

from nald_viscosity.finite_size import shear_cutoff


class TestFiniteSize(unittest.TestCase):
    def test_jcp_shear_cutoff(self):
        self.assertAlmostEqual(
            shear_cutoff(2.0 * np.pi, 4.0, 9.0),
            1.5,
            places=14,
        )


if __name__ == "__main__":
    unittest.main()
