"""Mode-resolved analysis helpers for the Kobayashi argon benchmark.

For the Markovian zero-frequency NALD expression,

    eta_M = nu * A_M,

with

    A_M = (1/V) sum_p Gamma_{M,p} / lambda_p^2.

A_M isolates the structural/mode contribution from the (still unresolved for
this NVE liquid benchmark) damping or memory-kernel contribution.

Negative instantaneous normal modes are retained.  A symmetric cutoff is
applied to |f_p| = sqrt(|lambda_p|)/(2*pi), matching the convention used in
the existing finite-size treatment.  Numerical null modes are removed
separately and must not be confused with a physical finite-size cutoff.
"""

from __future__ import annotations

import numpy as np


DEFAULT_NULL_LAMBDA_S2 = 1.0e12


def signed_frequency_thz(eigenvalues_s2) -> np.ndarray:
    """Signed cyclic frequency in THz from Hessian eigenvalues in s^-2."""
    lam = np.asarray(eigenvalues_s2, dtype=float)
    return np.sign(lam) * np.sqrt(np.abs(lam)) / (2.0 * np.pi * 1.0e12)


def numerical_null_mask(
    eigenvalues_s2,
    null_lambda_s2: float = DEFAULT_NULL_LAMBDA_S2,
) -> np.ndarray:
    """Mask numerical/translational/floppy null modes."""
    if null_lambda_s2 < 0.0:
        raise ValueError("null_lambda_s2 must be non-negative")
    return np.abs(np.asarray(eigenvalues_s2, dtype=float)) < null_lambda_s2


def structural_prefactor(
    eigenvalues_s2,
    gamma,
    volume_m3: float,
    cutoff_thz: float = 0.0,
    null_lambda_s2: float = DEFAULT_NULL_LAMBDA_S2,
) -> float:
    """Return A_M=(1/V) sum Gamma/lambda^2 in kg/m.

    Modes in the numerical nullspace are always removed.  If cutoff_thz > 0,
    both stable and unstable modes with |f_p| below the requested cyclic
    frequency are removed symmetrically.
    """
    lam = np.asarray(eigenvalues_s2, dtype=float)
    gam = np.asarray(gamma, dtype=float)
    if lam.shape != gam.shape:
        raise ValueError("eigenvalues and gamma must have the same shape")
    if volume_m3 <= 0.0:
        raise ValueError("volume_m3 must be positive")
    if cutoff_thz < 0.0:
        raise ValueError("cutoff_thz must be non-negative")

    null = numerical_null_mask(lam, null_lambda_s2)
    keep = ~null
    if cutoff_thz > 0.0:
        freq = np.abs(signed_frequency_thz(lam))
        keep &= freq >= cutoff_thz

    return float(np.sum(gam[keep] / lam[keep] ** 2) / volume_m3)


def branch_counts(
    eigenvalues_s2,
    null_lambda_s2: float = DEFAULT_NULL_LAMBDA_S2,
) -> dict[str, int]:
    """Count unstable, stable, and numerical-null modes."""
    lam = np.asarray(eigenvalues_s2, dtype=float)
    return {
        "negative": int(np.sum(lam < -null_lambda_s2)),
        "positive": int(np.sum(lam > null_lambda_s2)),
        "null": int(np.sum(np.abs(lam) < null_lambda_s2)),
    }


def null_gamma_fraction(
    eigenvalues_s2,
    gamma,
    null_lambda_s2: float = DEFAULT_NULL_LAMBDA_S2,
) -> float:
    """Fraction of total affine coupling carried by the numerical nullspace."""
    lam = np.asarray(eigenvalues_s2, dtype=float)
    gam = np.asarray(gamma, dtype=float)
    total = float(np.sum(gam))
    if total == 0.0:
        return 0.0
    return float(np.sum(gam[numerical_null_mask(lam, null_lambda_s2)]) / total)
