import tempfile
import unittest
from pathlib import Path
import numpy as np

from nald_viscosity.modes import (
    affine_vector,
    gamma_from_eigenvectors,
    read_affine,
    read_hessian,
)


class TestModes(unittest.TestCase):
    def test_affine_reader_sorts_ids(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "af.data"
            path.write_text(
                "2 16.0 4.0 5.0 6.0\n"
                "1 12.0 1.0 2.0 3.0\n",
                encoding="utf-8",
            )
            mass, xi = read_affine(path, skip_header=0)
            np.testing.assert_allclose(mass, [12.0, 16.0])
            np.testing.assert_allclose(
                xi, [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
            )

    def test_hessian_reader_allows_trailing_blank(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "hessian.data"
            path.write_text(
                "1 0 0\n"
                "0 2 0\n"
                "0 0 3\n\n",
                encoding="utf-8",
            )
            h = read_hessian(path, natoms=1)
            np.testing.assert_allclose(h, np.diag([1.0, 2.0, 3.0]))

    def test_mass_weighted_projection(self):
        masses = np.array([4.0])
        xi = np.array([[2.0, 4.0, 6.0]])
        vec = affine_vector(masses, xi, mass_weighted=True)
        np.testing.assert_allclose(vec, [1.0, 2.0, 3.0])
        gamma = gamma_from_eigenvectors(np.eye(3), vec)
        np.testing.assert_allclose(gamma, [1.0, 4.0, 9.0])


if __name__ == "__main__":
    unittest.main()
