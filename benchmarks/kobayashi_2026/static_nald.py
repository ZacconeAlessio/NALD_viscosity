"""Static NALD inputs for the Kobayashi 124 K argon benchmark.

The routines here use the collaborator-supplied MD conventions documented in
`source_protocol.md`: cubic periodic boundaries, the production TT/JHBV
two-body force law, and the model-dependent cutoff.

Only the pair-potential s=1.00 and s=0.00 benchmark models are handled here.
The ATM three-body term is intentionally excluded from this first benchmark.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from benchmarks.kobayashi_2026.tt_jhbv import (
    ARGON_MASS_KG,
    RMIN_M,
    source_force_curvature_SI,
    source_pair_cutoff_m,
)

N_ATOMS = 1372
DENSITY_KG_M3 = 1124.9


def source_box_length_m(
    natoms: int = N_ATOMS,
    density_kg_m3: float = DENSITY_KG_M3,
    mass_kg: float = ARGON_MASS_KG,
) -> float:
    """Cubic box length implied by the collaborator's MD source convention."""
    return (natoms * mass_kg / density_kg_m3) ** (1.0 / 3.0)


def cubic_cell_matrix(length_m: float) -> np.ndarray:
    return np.eye(3, dtype=float) * float(length_m)


def _cell_matrix(cell) -> np.ndarray:
    arr = np.asarray(cell, dtype=float)
    if arr.ndim == 0:
        return cubic_cell_matrix(float(arr))
    if arr.shape != (3, 3):
        raise ValueError("cell must be a scalar cubic length or a 3x3 matrix")
    return arr


def load_rxyz_frames(path: str | Path, natoms: int = N_ATOMS):
    """Load collaborator rxyz.dat blocks as [(step, xyz_m), ...]."""
    path = Path(path)
    lines = path.read_text(encoding="utf-8").splitlines()
    stride = natoms + 1
    if len(lines) % stride:
        raise ValueError(f"{path}: line count is not a multiple of {stride}")

    frames = []
    for start in range(0, len(lines), stride):
        step = int(lines[start].strip())
        xyz = np.loadtxt(lines[start + 1 : start + 1 + natoms], dtype=float)
        if xyz.shape != (natoms, 3):
            raise ValueError(f"{path}: malformed coordinate block at step {step}")
        frames.append((step, xyz))
    return frames


def _iter_pair_blocks(coords_m, cell, attraction_scale: float):
    """Yield vectorized pair blocks (i, js, dr, r, dUdr, d2Udr2)."""
    coords = np.asarray(coords_m, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError("coords_m must have shape (N,3)")

    h = _cell_matrix(cell)
    inv_h = np.linalg.inv(h)
    frac = coords @ inv_h.T
    cutoff = source_pair_cutoff_m(attraction_scale)

    for i in range(len(coords) - 1):
        ds = frac[i] - frac[i + 1 :]
        ds -= np.rint(ds)
        dr = ds @ h.T
        r = np.linalg.norm(dr, axis=1)
        keep = (r > 0.0) & (r <= cutoff)
        if not np.any(keep):
            continue

        js = np.nonzero(keep)[0] + i + 1
        dr = dr[keep]
        r = r[keep]
        u1, u2 = source_force_curvature_SI(r, attraction_scale)
        yield i, js, dr, r, np.asarray(u1), np.asarray(u2)


def forces_and_pair_metrics(coords_m, cell, attraction_scale: float):
    """Return forces, included-pair count, and mass-weighted Hessian trace."""
    coords = np.asarray(coords_m, dtype=float)
    forces = np.zeros_like(coords)
    pair_count = 0
    hessian_trace = 0.0

    for i, js, dr, r, u1, u2 in _iter_pair_blocks(coords, cell, attraction_scale):
        nhat = dr / r[:, None]
        fi = -u1[:, None] * nhat
        forces[i] += fi.sum(axis=0)
        np.add.at(forces, js, -fi)

        pair_count += len(js)
        # Each pair contributes two 3x3 diagonal blocks.  For a central
        # potential tr(A) = U'' + 2 U'/r.
        hessian_trace += float(np.sum(2.0 * (u2 + 2.0 * u1 / r)))

    return forces, pair_count, hessian_trace / ARGON_MASS_KG


def pair_forces(coords_m, cell, attraction_scale: float) -> np.ndarray:
    return forces_and_pair_metrics(coords_m, cell, attraction_scale)[0]


def dense_mass_weighted_hessian(coords_m, cell, attraction_scale: float) -> np.ndarray:
    """Construct the dense equal-mass dynamical matrix in s^-2.

    For one central pair, with n = r_ij/r,

        A = (U'' - U'/r) n n^T + (U'/r) I,

    and the Cartesian Hessian receives +A on ii/jj and -A on ij/ji.
    Division by the argon mass gives the mass-normalized Hessian used by NALD.
    """
    coords = np.asarray(coords_m, dtype=float)
    n_atoms = len(coords)
    hess = np.zeros((3 * n_atoms, 3 * n_atoms), dtype=float)
    eye = np.eye(3)

    for i, js, dr, r, u1, u2 in _iter_pair_blocks(coords, cell, attraction_scale):
        nhat = dr / r[:, None]
        outer = nhat[:, :, None] * nhat[:, None, :]
        blocks = (
            (u2 - u1 / r)[:, None, None] * outer
            + (u1 / r)[:, None, None] * eye
        )

        ii = slice(3 * i, 3 * i + 3)
        hess[ii, ii] += blocks.sum(axis=0)

        for a in range(3):
            jj_rows = 3 * js + a
            ii_rows = np.full_like(jj_rows, 3 * i + a)
            for b in range(3):
                jj_cols = 3 * js + b
                ii_cols = np.full_like(jj_cols, 3 * i + b)
                vals = blocks[:, a, b]

                np.add.at(hess, (jj_rows, jj_cols), vals)
                np.add.at(hess, (ii_rows, jj_cols), -vals)
                np.add.at(hess, (jj_rows, ii_cols), -vals)

    return hess / ARGON_MASS_KG


def _affine_deform(coords_m, cell, channel: str, strain: float):
    coords = np.asarray(coords_m, dtype=float)
    h = _cell_matrix(cell)

    if channel == "bulk":
        # epsilon_v = Delta V / V.  Isotropic linear scale is (1+epsilon_v)^(1/3).
        scale = (1.0 + strain) ** (1.0 / 3.0)
        a = np.eye(3) * scale
    elif channel == "shear":
        # xz engineering shear: x -> x + gamma z.
        a = np.eye(3)
        a[0, 2] = strain
    else:
        raise ValueError("channel must be 'bulk' or 'shear'")

    return coords @ a.T, a @ h


def affine_force(
    coords_m,
    cell,
    attraction_scale: float,
    channel: str,
    delta: float = 1.0e-6,
) -> np.ndarray:
    """Centered finite-difference affine force Xi = -df/d(strain)."""
    if delta <= 0.0:
        raise ValueError("delta must be positive")

    minus_coords, minus_cell = _affine_deform(coords_m, cell, channel, -0.5 * delta)
    plus_coords, plus_cell = _affine_deform(coords_m, cell, channel, +0.5 * delta)

    f_minus = pair_forces(minus_coords, minus_cell, attraction_scale)
    f_plus = pair_forces(plus_coords, plus_cell, attraction_scale)
    return (f_minus - f_plus) / delta


def mass_weight_affine_force(xi: np.ndarray) -> np.ndarray:
    """Return Xi_i/sqrt(m_Ar), flattened for projection onto mass-normalized modes."""
    return (np.asarray(xi, dtype=float) / np.sqrt(ARGON_MASS_KG)).reshape(-1)
