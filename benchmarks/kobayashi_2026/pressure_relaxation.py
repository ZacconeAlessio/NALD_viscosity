"""Pressure-relaxation diagnostics for collaborator-supplied bulk_corr.dat."""

from __future__ import annotations

from pathlib import Path

import numpy as np


def load_bulk_corr(path: str | Path) -> np.ndarray:
    corr = np.loadtxt(path, dtype=float)
    if corr.ndim != 2 or corr.shape[1] != 2:
        raise ValueError("bulk_corr.dat must contain two columns")
    if not np.all(np.diff(corr[:, 0]) > 0.0):
        raise ValueError("time column must be strictly increasing")
    return corr


def cumulative_trapezoid(time_s, integrand_pa) -> np.ndarray:
    t = np.asarray(time_s, dtype=float)
    y = np.asarray(integrand_pa, dtype=float)
    increments = 0.5 * (y[1:] + y[:-1]) * np.diff(t)
    return np.concatenate(([0.0], np.cumsum(increments)))


def relaxation_summary(path: str | Path) -> dict[str, float]:
    corr = load_bulk_corr(path)
    time_s = corr[:, 0]
    integrand = corr[:, 1]
    cumulative = cumulative_trapezoid(time_s, integrand)
    zeta = float(cumulative[-1])
    c0 = float(integrand[0])

    result = {
        "corr_rows": len(corr),
        "dt_ps": float((time_s[1] - time_s[0]) * 1.0e12),
        "tmax_ps": float(time_s[-1] * 1.0e12),
        "C0_GPa": c0 / 1.0e9,
        "zeta_integral_mPa_s": zeta * 1.0e3,
        "tau_integral_ps": zeta / c0 * 1.0e12,
    }

    for fraction in (0.50, 0.80, 0.90, 0.95):
        target = fraction * zeta
        index = int(np.flatnonzero(cumulative >= target)[0])
        result[f"t{int(100*fraction)}_cumulative_ps"] = float(time_s[index] * 1.0e12)

    return result
