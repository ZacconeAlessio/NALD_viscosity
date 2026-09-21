import unittest
import numpy as np

from nald_viscosity.core import (
    GAMMA_REAL_MW_TO_SI,
    loss_modulus_si,
    select_modes,
    zero_frequency_viscosity_si,
)


class TestCore(unittest.TestCase):
    def test_single_mode_zero_frequency(self):
        lam = np.array([4.0])
        gamma = np.array([2.0])
        eta = zero_frequency_viscosity_si(lam, gamma, 5.0, 3.0)
        self.assertAlmostEqual(eta, 0.075, places=14)

    def test_negative_inm_contributes(self):
        lam = np.array([-4.0, 4.0])
        gamma = np.array([2.0, 2.0])
        eta = zero_frequency_viscosity_si(lam, gamma, 5.0, 3.0)
        self.assertAlmostEqual(eta, 0.15, places=14)

    def test_low_frequency_ratio(self):
        lam = np.array([-7.0, 4.0, 11.0])
        gamma = np.array([0.5, 2.0, 1.5])
        volume = 5.0
        nu = 3.0
        eta0 = zero_frequency_viscosity_si(lam, gamma, volume, nu)
        omega = np.array([1.0e-7])
        loss = loss_modulus_si(omega, lam, gamma, volume, nu)
        self.assertAlmostEqual(float(loss[0] / omega[0]), eta0, places=12)

    def test_cutoff_symmetric_in_lambda_sign(self):
        lam = np.array([-9.0, -1.0, 0.0, 1.0, 9.0])
        gamma = np.ones(5)
        kept, _, info = select_modes(
            lam, gamma, zero_tol=1.0e-12, cutoff_frequency=2.0
        )
        np.testing.assert_allclose(kept, [-9.0, 9.0])
        self.assertEqual(info["removed_zero"], 1)
        self.assertEqual(info["removed_cutoff"], 2)

    def test_real_unit_gamma_conversion(self):
        self.assertAlmostEqual(
            GAMMA_REAL_MW_TO_SI,
            2906915.7802365357,
            places=8,
        )


if __name__ == "__main__":
    unittest.main()
