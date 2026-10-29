"""Frequency-resolved Green-Kubo transforms for stress correlations."""

from __future__ import annotations

import numpy as np


def complex_viscosity_from_integrand(time_s, integrand_pa, omega_s1):
    """Return integral C(t) exp(-i*Omega*t) dt in Pa s.

    The input integrand should already include the Green-Kubo prefactor, e.g.
    V <delta P(t) delta P(0)>/(k_B T) for bulk viscosity.
    """
    t = np.asarray(time_s, dtype=float)
    c = np.asarray(integrand_pa, dtype=float)
    omega = np.asarray(omega_s1, dtype=float)
    if t.ndim != 1 or c.ndim != 1 or t.shape != c.shape:
        raise ValueError("time and integrand must be matching 1D arrays")
    if len(t) < 2 or t[0] != 0.0 or np.any(np.diff(t) <= 0.0):
        raise ValueError("time must start at zero and increase strictly")

    out = np.empty(omega.shape, dtype=complex)
    dt = np.diff(t)
    for index, w in np.ndenumerate(omega):
        y = c * np.exp(-1j * float(w) * t)
        out[index] = np.sum(0.5 * (y[1:] + y[:-1]) * dt)
    return out


def dissipative_viscosity_spectrum(time_s, integrand_pa, omega_s1):
    """Return Re[eta*(Omega)] = integral C(t) cos(Omega t) dt."""
    return np.real(
        complex_viscosity_from_integrand(time_s, integrand_pa, omega_s1)
    )
