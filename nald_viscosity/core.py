"""Modal NALD response and unit conversion.

For a deformation channel M (G for shear, K for bulk),

    M''(Omega) = (1/V) sum_p Gamma_p nu Omega /
                 ((lambda_p - Omega^2)^2 + (nu Omega)^2)

and, for a constant Markovian damping rate,

    lim_(Omega->0) M''(Omega)/Omega
        = (1/V) sum_p Gamma_p nu / lambda_p^2.

The bulk channel uses exactly the same modal kernel; only the affine-force
field (and therefore Gamma_p) changes.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np

AVOGADRO = 6.02214076e23
AMU_KG = 1.66053906660e-27
KCAL_J = 4184.0
ANGSTROM_M = 1.0e-10
THZ_SI = 1.0e12
ANGSTROM3_TO_M3 = 1.0e-30

FORCE_REAL_TO_N = (KCAL_J / AVOGADRO) / ANGSTROM_M
GAMMA_REAL_MW_TO_SI = FORCE_REAL_TO_N**2 / AMU_KG


def _load_vector(path: str | Path) -> np.ndarray:
    values = np.loadtxt(path, dtype=float)
    return np.asarray(values, dtype=float).reshape(-1)


def load_modes(
    eigenvalues_path: str | Path,
    gamma_path: str | Path,
) -> tuple[np.ndarray, np.ndarray]:
    """Load matching one-column eigenvalue and Gamma arrays."""
    lam = _load_vector(eigenvalues_path)
    gamma = _load_vector(gamma_path)
    if lam.size != gamma.size:
        raise ValueError(
            f"eigenvalues ({lam.size}) and gamma ({gamma.size}) differ in length"
        )
    if not (np.all(np.isfinite(lam)) and np.all(np.isfinite(gamma))):
        raise ValueError("mode data contain non-finite values")
    if np.any(gamma < -1.0e-12):
        raise ValueError("Gamma must be non-negative because Gamma=(e.Xi)^2")
    return lam, np.maximum(gamma, 0.0)


def convert_to_si(
    lam: np.ndarray,
    gamma: np.ndarray,
    volume: float,
    input_units: str,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Convert supported mode conventions to SI."""
    if input_units == "si":
        return lam, gamma, volume
    if input_units == "lammps-real-massweighted":
        return (
            lam * THZ_SI**2,
            gamma * GAMMA_REAL_MW_TO_SI,
            volume * ANGSTROM3_TO_M3,
        )
    raise ValueError(f"unknown input unit convention: {input_units}")


def select_modes(
    lam_native: np.ndarray,
    gamma: np.ndarray,
    zero_tol: float = 1.0e-10,
    cutoff_frequency: float = 0.0,
) -> tuple[np.ndarray, np.ndarray, dict[str, int]]:
    """Remove numerical zero modes and, optionally, a finite-size band.

    Cutoffs are applied in the native eigenvalue convention before SI
    conversion. Negative INMs outside the excluded band are retained.
    """
    if zero_tol < 0 or cutoff_frequency < 0:
        raise ValueError("zero_tol and cutoff_frequency must be non-negative")
    abs_lam = np.abs(lam_native)
    keep_zero = abs_lam > zero_tol
    if cutoff_frequency > 0.0:
        keep_cutoff = np.sqrt(abs_lam) >= cutoff_frequency
    else:
        keep_cutoff = np.ones_like(keep_zero, dtype=bool)
    keep = keep_zero & keep_cutoff
    info = {
        "total": int(lam_native.size),
        "removed_zero": int(np.count_nonzero(~keep_zero)),
        "removed_cutoff": int(np.count_nonzero(keep_zero & ~keep_cutoff)),
        "kept": int(np.count_nonzero(keep)),
    }
    return lam_native[keep], gamma[keep], info


def loss_modulus_si(
    omega_s1: np.ndarray,
    lam_s2: np.ndarray,
    gamma_si: np.ndarray,
    volume_m3: float,
    nu_s1: float,
) -> np.ndarray:
    """Return the loss modulus in Pa."""
    omega_s1 = np.asarray(omega_s1, dtype=float)
    if volume_m3 <= 0 or nu_s1 <= 0:
        raise ValueError("volume and damping rate must be positive")
    out = np.empty_like(omega_s1)
    for j, omega in enumerate(omega_s1):
        den = (lam_s2 - omega**2) ** 2 + (nu_s1 * omega) ** 2
        out[j] = np.sum(gamma_si * nu_s1 * omega / den) / volume_m3
    return out


def zero_frequency_viscosity_si(
    lam_s2: np.ndarray,
    gamma_si: np.ndarray,
    volume_m3: float,
    nu_s1: float,
) -> float:
    """Return lim M''(Omega)/Omega in Pa s."""
    if volume_m3 <= 0 or nu_s1 <= 0:
        raise ValueError("volume and damping rate must be positive")
    if np.any(lam_s2 == 0.0):
        raise ValueError("zero eigenvalues must be removed before evaluating viscosity")
    return float(np.sum(gamma_si * nu_s1 / (lam_s2**2)) / volume_m3)
