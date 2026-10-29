"""Second-order Mori memory diagnostics from pressure correlations.

For a normalized correlation C(t), consider the generalized oscillator equation

    C''(t) + omega0^2 C(t)
      + integral_0^t nu(t-t') C'(t') dt' = 0,

with C(0)=1 and C'(0)=0.  Its Laplace-domain memory is

    nu_tilde(s) =
      [s - (s^2 + omega0^2) C_tilde(s)] / [s C_tilde(s) - 1].

At s=0,

    nu_tilde(0) = omega0^2 * integral_0^infinity C(t) dt.

This has the same dimensions as the generalized NALD friction transform.
The mapping from the pressure collective variable to the microscopic NALD
memory is only a diagnostic at this stage; it is not assumed to be exact.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


def _normalized(time_s, corr):
    t = np.asarray(time_s, dtype=float)
    c = np.asarray(corr, dtype=float)
    if t.ndim != 1 or c.ndim != 1 or t.shape != c.shape:
        raise ValueError("time and correlation must be matching 1D arrays")
    if len(t) < 3 or t[0] != 0.0 or np.any(np.diff(t) <= 0.0):
        raise ValueError("time must start at zero and increase strictly")
    if c[0] == 0.0:
        raise ValueError("correlation C(0) must be non-zero")
    return t, c / c[0]


def estimate_omega0_squared(
    time_s,
    corr,
    window_s: float = 4.0e-14,
) -> float:
    """Estimate omega0^2=-C''(0) with an even quartic short-time fit."""
    if window_s <= 0.0:
        raise ValueError("window_s must be positive")
    t, c = _normalized(time_s, corr)
    keep = t <= window_s
    if np.count_nonzero(keep) < 3:
        raise ValueError("curvature window contains too few samples")

    x = t[keep] / window_s
    design = np.column_stack((x**2, x**4))
    b2, _ = np.linalg.lstsq(design, c[keep] - 1.0, rcond=None)[0]
    omega0_sq = -2.0 * b2 / window_s**2
    if omega0_sq <= 0.0:
        raise ValueError("estimated -C''(0) is not positive")
    return float(omega0_sq)


def laplace_correlation(time_s, corr, s_s1: float) -> float:
    """Return the one-sided Laplace transform of normalized C(t)."""
    if s_s1 < 0.0:
        raise ValueError("Laplace rate must be non-negative")
    t, c = _normalized(time_s, corr)
    y = c * np.exp(-s_s1 * t)
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(t)))


def second_order_memory_laplace(
    time_s,
    corr,
    s_s1: float,
    omega0_sq_s2: float | None = None,
    curvature_window_s: float = 4.0e-14,
) -> float:
    """Return nu_tilde(s) [s^-1] for the second-order Mori diagnostic."""
    if s_s1 < 0.0:
        raise ValueError("Laplace rate must be non-negative")
    if omega0_sq_s2 is None:
        omega0_sq_s2 = estimate_omega0_squared(
            time_s, corr, window_s=curvature_window_s
        )
    if omega0_sq_s2 <= 0.0:
        raise ValueError("omega0_sq_s2 must be positive")

    chat = laplace_correlation(time_s, corr, s_s1)
    denominator = s_s1 * chat - 1.0
    if abs(denominator) < 1.0e-14:
        raise ValueError("memory transform denominator is numerically singular")
    return float(
        (s_s1 - (s_s1**2 + omega0_sq_s2) * chat) / denominator
    )


def pressure_memory_summary(
    path: str | Path,
    curvature_window_s: float = 4.0e-14,
) -> dict[str, float]:
    """Return robust zero-frequency pressure-memory diagnostics."""
    data = np.loadtxt(path, dtype=float)
    if data.ndim != 2 or data.shape[1] != 2:
        raise ValueError("pressure correlation file must have two columns")

    time_s = data[:, 0]
    corr = data[:, 1]
    omega0_sq = estimate_omega0_squared(
        time_s, corr, window_s=curvature_window_s
    )
    tau_s = laplace_correlation(time_s, corr, 0.0)
    nu0 = omega0_sq * tau_s

    return {
        "curvature_window_ps": curvature_window_s * 1.0e12,
        "omega0_sq_s-2": omega0_sq,
        "bare_frequency_THz": np.sqrt(omega0_sq) / (2.0 * np.pi * 1.0e12),
        "tau_pressure_ps": tau_s * 1.0e12,
        "nu_pressure_zero_s-1": nu0,
    }
