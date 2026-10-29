"""Independent acoustic finite-size diagnostics for the Kobayashi argon benchmark.

At each box wavevector q=n*2*pi/L we project the instantaneous Hessian onto
the three Cartesian plane-wave polarizations.  For a disordered snapshot the
projected 3x3 dynamical matrix need not be diagonal in the Cartesian basis, so
we diagonalize that matrix and identify the mode with the largest projection
onto q as longitudinal.  The other two modes are classified as transverse.

This avoids treating pure Cartesian Rayleigh quotients as exact L/T modes when
longitudinal-transverse mixing is present in an individual disordered snapshot.
"""

from __future__ import annotations

import numpy as np

from benchmarks.kobayashi_2026.static_nald import _cell_matrix, _iter_pair_blocks
from benchmarks.kobayashi_2026.tt_jhbv import ARGON_MASS_KG


def _signed_frequency_thz(eigenvalues_s2) -> np.ndarray:
    lam = np.asarray(eigenvalues_s2, dtype=float)
    return np.sign(lam) * np.sqrt(np.abs(lam)) / (2.0 * np.pi * 1.0e12)


def classify_projected_dynamical_matrix(dynamical_s2, q_axis: int):
    """Diagonalize a 3x3 projected dynamical matrix and classify L/T modes.

    Returns (f_longitudinal, f_transverse[2], longitudinal_projection).
    Frequencies are signed and expressed in THz.
    """
    dyn = np.asarray(dynamical_s2, dtype=float)
    if dyn.shape != (3, 3):
        raise ValueError("dynamical_s2 must be a 3x3 matrix")
    if q_axis not in (0, 1, 2):
        raise ValueError("q_axis must be 0, 1, or 2")

    dyn = 0.5 * (dyn + dyn.T)
    eigenvalues, eigenvectors = np.linalg.eigh(dyn)
    frequencies = _signed_frequency_thz(eigenvalues)

    longitudinal_weights = eigenvectors[q_axis, :] ** 2
    i_long = int(np.argmax(longitudinal_weights))
    i_trans = [i for i in range(3) if i != i_long]

    return (
        float(frequencies[i_long]),
        np.asarray(frequencies[i_trans], dtype=float),
        float(longitudinal_weights[i_long]),
    )


def plane_wave_rayleigh(
    coords_m,
    cell,
    attraction_scale: float,
    harmonic_index: int = 1,
):
    """Return projected longitudinal/transverse frequencies at n*2*pi/L."""
    if harmonic_index < 1:
        raise ValueError("harmonic_index must be >= 1")

    coords = np.asarray(coords_m, dtype=float)
    h = _cell_matrix(cell)
    lengths = np.diag(h)
    if not np.allclose(h, np.diag(lengths)):
        raise ValueError("plane_wave_rayleigh currently expects an orthogonal cell")
    if not np.allclose(lengths, lengths[0]):
        raise ValueError("plane_wave_rayleigh currently expects a cubic cell")

    length = float(lengths[0])
    projected = np.zeros((3, 3, 3), dtype=float)
    eye = np.eye(3)

    for _, _, dr, r, u1, u2 in _iter_pair_blocks(coords, h, attraction_scale):
        nhat = dr / r[:, None]
        outer = nhat[:, :, None] * nhat[:, None, :]
        blocks = (
            (u2 - u1 / r)[:, None, None] * outer
            + (u1 / r)[:, None, None] * eye
        )

        for q_axis in range(3):
            phase_weight = 2.0 - 2.0 * np.cos(
                2.0 * np.pi * harmonic_index * dr[:, q_axis] / length
            )
            projected[q_axis] += np.sum(
                blocks * phase_weight[:, None, None],
                axis=0,
            )

    projected /= len(coords) * ARGON_MASS_KG

    longitudinal = []
    transverse = []
    longitudinal_projection = []
    for q_axis in range(3):
        f_long, f_trans, weight = classify_projected_dynamical_matrix(
            projected[q_axis],
            q_axis,
        )
        longitudinal.append(f_long)
        transverse.extend(f_trans)
        longitudinal_projection.append(weight)

    longitudinal = np.asarray(longitudinal)
    transverse = np.asarray(transverse)
    longitudinal_projection = np.asarray(longitudinal_projection)

    f_long = float(np.mean(longitudinal))
    f_trans = float(np.mean(transverse))

    return {
        "harmonic_index": harmonic_index,
        "q_nm-1": (2.0 * np.pi * harmonic_index / length) * 1.0e-9,
        "f_longitudinal_THz": f_long,
        "f_longitudinal_std_THz": float(np.std(longitudinal, ddof=1)),
        "f_transverse_THz": f_trans,
        "f_transverse_std_THz": float(np.std(transverse, ddof=1)),
        "longitudinal_projection_mean": float(np.mean(longitudinal_projection)),
        "longitudinal_projection_min": float(np.min(longitudinal_projection)),
        "c_longitudinal_m_s": f_long * 1.0e12 * length / harmonic_index,
        "c_transverse_m_s": f_trans * 1.0e12 * length / harmonic_index,
    }
