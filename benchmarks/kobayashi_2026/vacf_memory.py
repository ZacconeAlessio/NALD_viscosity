"""Independent GLE memory reconstruction from a velocity autocorrelation.

For a normalized VACF C_v(t), the first-order Mori/GLE equation is

    dC_v/dt = - integral_0^t K(t-t') C_v(t') dt'.

Its one-sided Laplace transform obeys

    K_tilde(s) = 1/C_tilde(s) - s.

K_tilde has units s^-1 and can be supplied directly to the generalized NALD
friction kernel when the same Fourier/Laplace convention is used.
"""

from __future__ import annotations

import numpy as np


def _normalized(time_s, corr):
    t = np.asarray(time_s, dtype=float)
    c = np.asarray(corr, dtype=float)
    if t.ndim != 1 or c.ndim != 1 or t.shape != c.shape:
        raise ValueError("time and correlation must be matching 1D arrays")
    if len(t) < 2 or t[0] != 0.0 or np.any(np.diff(t) <= 0.0):
        raise ValueError("time must start at zero and increase strictly")
    if c[0] == 0.0:
        raise ValueError("correlation C(0) must be non-zero")
    return t, c / c[0]


def one_sided_laplace(time_s, corr, s_s1: complex) -> complex:
    """Numerically evaluate integral C(t) exp(-s t) dt."""
    if np.real(s_s1) < 0.0:
        raise ValueError("Laplace transform requires Re(s) >= 0")
    t, c = _normalized(time_s, corr)
    y = c * np.exp(-complex(s_s1) * t)
    return complex(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(t)))


def vacf_memory_laplace(time_s, corr, s_s1: complex) -> complex:
    """Return K_tilde(s)=1/C_tilde(s)-s in s^-1."""
    chat = one_sided_laplace(time_s, corr, s_s1)
    if abs(chat) == 0.0:
        raise ValueError("VACF Laplace transform is zero")
    return 1.0 / chat - complex(s_s1)


def vacf_memory_frequency(
    time_s,
    corr,
    omega_s1,
    regularization_s1: float = 0.0,
) -> np.ndarray:
    """Return causal complex memory K_tilde(i*Omega+epsilon).

    The returned convention is compatible with the NALD denominator

        lambda - Omega^2 + i Omega K_tilde(Omega).
    """
    if regularization_s1 < 0.0:
        raise ValueError("regularization_s1 must be non-negative")
    omega = np.asarray(omega_s1, dtype=float)
    return np.asarray(
        [
            vacf_memory_laplace(
                time_s,
                corr,
                complex(regularization_s1, w),
            )
            for w in omega
        ],
        dtype=complex,
    )
