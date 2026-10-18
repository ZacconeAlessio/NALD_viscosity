"""Independent acoustic finite-size diagnostics for the Kobayashi argon benchmark."""

from __future__ import annotations

import numpy as np

from benchmarks.kobayashi_2026.static_nald import _cell_matrix, _iter_pair_blocks
from benchmarks.kobayashi_2026.tt_jhbv import ARGON_MASS_KG


def plane_wave_rayleigh(coords_m, cell, attraction_scale: float, harmonic_index: int = 1):
    """Return longitudinal/transverse plane-wave Rayleigh frequencies at n*2*pi/L."""
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
    numerator = np.zeros((3, 3), dtype=float)

    for _, _, dr, r, u1, u2 in _iter_pair_blocks(coords, h, attraction_scale):
        nhat = dr / r[:, None]
        block_diag = (u2 - u1 / r)[:, None] * nhat**2 + (u1 / r)[:, None]
        for q_axis in range(3):
            phase_weight = 2.0 - 2.0 * np.cos(
                2.0 * np.pi * harmonic_index * dr[:, q_axis] / length
            )
            numerator[q_axis] += np.sum(block_diag * phase_weight[:, None], axis=0)

    lam = numerator / (len(coords) * ARGON_MASS_KG)
    freq_thz = np.sign(lam) * np.sqrt(np.abs(lam)) / (2.0 * np.pi * 1.0e12)

    longitudinal = np.array([freq_thz[a, a] for a in range(3)])
    transverse = np.array([freq_thz[a, p] for a in range(3) for p in range(3) if p != a])

    f_long = float(np.mean(longitudinal))
    f_trans = float(np.mean(transverse))

    return {
        "harmonic_index": harmonic_index,
        "q_nm-1": (2.0 * np.pi * harmonic_index / length) * 1.0e-9,
        "f_longitudinal_THz": f_long,
        "f_longitudinal_std_THz": float(np.std(longitudinal, ddof=1)),
        "f_transverse_THz": f_trans,
        "f_transverse_std_THz": float(np.std(transverse, ddof=1)),
        "c_longitudinal_m_s": f_long * 1.0e12 * length / harmonic_index,
        "c_transverse_m_s": f_trans * 1.0e12 * length / harmonic_index,
    }
