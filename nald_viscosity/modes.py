"""Hessian I/O and affine-force projection."""

from __future__ import annotations

from pathlib import Path
import numpy as np


def read_affine(path: str | Path, skip_header: int = 9) -> tuple[np.ndarray, np.ndarray]:
    """Read an affine-force dump.

    Supported numeric columns:
      id mass Xi_x Xi_y Xi_z  (recommended)
      mass Xi_x Xi_y Xi_z     (historical)
    """
    path = Path(path)
    data = np.genfromtxt(path, skip_header=skip_header)
    if data.ndim == 1:
        data = data[None, :]
    if data.shape[1] >= 5:
        data = data[np.argsort(data[:, 0])]
        mass = data[:, 1]
        xi = data[:, 2:5]
    elif data.shape[1] == 4:
        mass = data[:, 0]
        xi = data[:, 1:4]
    else:
        raise ValueError(f"{path}: expected 4 or 5 numeric columns")
    if np.any(mass <= 0):
        raise ValueError(f"{path}: all atomic masses must be positive")
    return mass.astype(float), xi.astype(float)


def read_hessian(path: str | Path, natoms: int) -> np.ndarray:
    """Read the repository/LAMMPS block format and symmetrize it."""
    path = Path(path)
    dim = 3 * natoms
    hessian = np.empty((dim, dim), dtype=float)
    with path.open("r", encoding="utf-8") as handle:
        for row in range(dim):
            values = np.empty(dim, dtype=float)
            k = 0
            for _ in range(natoms):
                line = handle.readline()
                if not line:
                    raise ValueError(f"{path}: ended early while reading row {row}")
                block = np.fromstring(line, sep=" ")
                if block.size != 3:
                    raise ValueError(f"{path}: expected 3 Hessian entries per line")
                values[k : k + 3] = block
                k += 3
            hessian[row] = values
        if any(line.strip() for line in handle):
            raise ValueError(f"{path}: extra nonblank Hessian lines found")
    return 0.5 * (hessian + hessian.T)


def affine_vector(
    masses_amu: np.ndarray,
    xi: np.ndarray,
    mass_weighted: bool = True,
) -> np.ndarray:
    """Flatten Xi, optionally using Xi_i/sqrt(m_i)."""
    if xi.shape != (len(masses_amu), 3):
        raise ValueError("Xi must have shape (N,3)")
    if mass_weighted:
        return (xi / np.sqrt(masses_amu)[:, None]).reshape(-1)
    return xi.reshape(-1)


def gamma_from_eigenvectors(eigenvectors: np.ndarray, xi_vector: np.ndarray) -> np.ndarray:
    """Return Gamma_p=(e_p.Xi)^2 for columns e_p of eigenvectors."""
    if eigenvectors.shape[0] != xi_vector.size:
        raise ValueError("eigenvector dimension and affine-force vector differ")
    return (eigenvectors.T @ xi_vector) ** 2
