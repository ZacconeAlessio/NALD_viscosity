"""Low-frequency mode diagnostics for the Kobayashi benchmark.

These helpers characterize modes that fall below an acoustic finite-size
frequency without assuming that every such mode is an acoustic finite-size
artifact.  Two independent descriptors are reported:

- participation ratio (PR), measuring spatial delocalization;
- overlap with the q_min plane-wave subspace of the cubic simulation box.

The descriptors are diagnostics only.  No universal threshold is imposed.
"""

from __future__ import annotations

import numpy as np

from benchmarks.kobayashi_2026.mode_analysis import numerical_null_mask


def participation_ratio(eigenvectors: np.ndarray, natoms: int) -> np.ndarray:
    """Return atom-wise participation ratio for each normalized eigenmode."""
    vec = np.asarray(eigenvectors, dtype=float)
    if vec.shape[0] != 3 * natoms:
        raise ValueError("eigenvectors must have shape (3*natoms, nmodes)")
    atom_amp2 = np.sum(vec.reshape(natoms, 3, -1) ** 2, axis=1)
    denom = natoms * np.sum(atom_amp2**2, axis=0)
    return 1.0 / denom


def qmin_plane_wave_weights(
    coords_m: np.ndarray,
    eigenvectors: np.ndarray,
    box_length_m: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return total, longitudinal, and transverse q_min plane-wave weights.

    The three smallest reciprocal vectors q=2*pi/L along x,y,z are used.
    For each direction, the complex plane-wave overlap is resolved into
    Cartesian polarizations.  The reported weights are summed over the three
    q directions.  Because particle positions are disordered, this is a
    diagnostic projection rather than an exact orthogonal mode decomposition.
    """
    coords = np.asarray(coords_m, dtype=float)
    vec = np.asarray(eigenvectors, dtype=float)
    natoms = len(coords)
    if coords.shape != (natoms, 3):
        raise ValueError("coords_m must have shape (N,3)")
    if vec.shape[0] != 3 * natoms:
        raise ValueError("eigenvectors must have shape (3*N, nmodes)")
    if box_length_m <= 0.0:
        raise ValueError("box_length_m must be positive")

    modes = vec.reshape(natoms, 3, -1)
    longitudinal = np.zeros(vec.shape[1], dtype=float)
    transverse = np.zeros(vec.shape[1], dtype=float)

    for q_axis in range(3):
        phase = np.exp(
            -1j * 2.0 * np.pi * coords[:, q_axis] / box_length_m
        ) / np.sqrt(natoms)
        projection = np.einsum("i,iam->am", phase, modes, optimize=True)
        weight = np.abs(projection) ** 2
        longitudinal += weight[q_axis]
        transverse += np.sum(weight, axis=0) - weight[q_axis]

    return longitudinal + transverse, longitudinal, transverse


def structural_mode_contribution(
    eigenvalues_s2,
    gamma,
    volume_m3: float,
    null_lambda_s2: float = 1.0e12,
) -> np.ndarray:
    """Return per-mode contribution Gamma/(V*lambda^2) in kg/m."""
    lam = np.asarray(eigenvalues_s2, dtype=float)
    gam = np.asarray(gamma, dtype=float)
    if lam.shape != gam.shape:
        raise ValueError("eigenvalues and gamma must have the same shape")
    if volume_m3 <= 0.0:
        raise ValueError("volume_m3 must be positive")

    out = np.zeros_like(lam, dtype=float)
    keep = ~numerical_null_mask(lam, null_lambda_s2)
    out[keep] = gam[keep] / (volume_m3 * lam[keep] ** 2)
    return out


def summarize_subcutoff_modes(
    eigenvalues_s2,
    gamma,
    participation,
    qmin_weight,
    cutoff_thz: float,
    volume_m3: float,
    null_lambda_s2: float = 1.0e12,
) -> dict[str, float]:
    """Summarize static structural weight carried below a frequency cutoff."""
    lam = np.asarray(eigenvalues_s2, dtype=float)
    pr = np.asarray(participation, dtype=float)
    qw = np.asarray(qmin_weight, dtype=float)
    if not (lam.shape == pr.shape == qw.shape):
        raise ValueError("mode arrays must have matching shape")
    if cutoff_thz < 0.0:
        raise ValueError("cutoff_thz must be non-negative")

    contribution = structural_mode_contribution(
        lam, gamma, volume_m3, null_lambda_s2
    )
    frequency = np.sqrt(np.abs(lam)) / (2.0 * np.pi * 1.0e12)
    valid = ~numerical_null_mask(lam, null_lambda_s2)
    low = valid & (frequency < cutoff_thz)

    total_low = float(np.sum(contribution[low]))
    if total_low > 0.0:
        weights = contribution[low] / total_low
        weighted_pr = float(np.sum(weights * pr[low]))
        weighted_q = float(np.sum(weights * qw[low]))
        low_indices = np.flatnonzero(low)
        top = int(low_indices[np.argmax(contribution[low])])
        top_fraction = float(contribution[top] / total_low)
        top_pr = float(pr[top])
        top_q = float(qw[top])
        top_frequency = float(frequency[top])
        top_negative = float(lam[top] < 0.0)
    else:
        weighted_pr = weighted_q = top_fraction = 0.0
        top_pr = top_q = top_frequency = 0.0
        top_negative = 0.0

    high = valid & ~low
    high_sum = float(np.sum(contribution[high]))
    all_sum = total_low + high_sum

    return {
        "subcutoff_modes": int(np.count_nonzero(low)),
        "A_below_cut_kg_m": total_low,
        "A_above_cut_kg_m": high_sum,
        "A_all_over_above": float(all_sum / high_sum) if high_sum else np.inf,
        "A_weighted_participation_ratio": weighted_pr,
        "A_weighted_qmin_plane_wave_weight": weighted_q,
        "top1_low_A_fraction": top_fraction,
        "top1_participation_ratio": top_pr,
        "top1_qmin_plane_wave_weight": top_q,
        "top1_frequency_THz": top_frequency,
        "top1_negative": top_negative,
    }
