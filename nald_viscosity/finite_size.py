"""Finite-size cutoff helpers."""

from __future__ import annotations
import math


def shear_cutoff(box_length: float, density: float, shear_modulus: float) -> float:
    """Return omega_min=(2*pi/L)*sqrt(G_s/rho).

    In SI, use L [m], rho [kg m^-3], G_s [Pa] and the result is [s^-1].
    The same expression can be used in a self-consistent reduced unit system.

    This prescription is documented for the shear channel in Singh et al.,
    J. Chem. Phys. 162, 244504 (2025). This function intentionally does not
    invent a corresponding bulk/longitudinal cutoff.
    """
    if box_length <= 0 or density <= 0 or shear_modulus < 0:
        raise ValueError("L and density must be positive and G_s non-negative")
    return (2.0 * math.pi / box_length) * math.sqrt(shear_modulus / density)
