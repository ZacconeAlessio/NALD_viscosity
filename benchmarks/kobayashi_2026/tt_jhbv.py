"""Tang-Toennies/JHBV pair potential for the Kobayashi argon benchmark.

The published potential is written in kelvin and nanometres:

U(r)/k_B = A exp(a1 r + a2 r^2 + a_-1/r + a_-2/r^2)
           - sum_{p=6,8,...,16} C_p r^{-p} f_p(b r),

where
    f_p(x) = 1 - exp(-x) sum_{k=0}^p x^k/k!.

Parameters follow Alboul & Lishchuk, Phys. Rev. E 105, 054135 (2022),
Table I, together with the Jager et al. higher dispersion coefficients.
C16 uses the published corrigendum value 1.17006343e-6 K nm^16.

For the attraction-scaled model U_{2,s}, forces are unchanged below the
potential minimum and multiplied by s above it, as in the WCA-like
decomposition used by Kobayashi et al.
"""

from __future__ import annotations

import math

KB = 1.380649e-23  # J/K
NM = 1.0e-9        # m

A_K = 4.61330146e7
A1_NM1 = -2.98337630e1
A2_NM2 = -9.71208881
AM1_NM = 2.75206827e-2
AM2_NM2 = -1.01489050e-2
B_NM1 = 4.02517211e1

C_K_NM = {
    6: 4.42812017e-1,
    8: 3.26707684e-2,
    10: 2.45656537e-3,
    12: 1.88246247e-4,
    14: 1.47012192e-5,
    16: 1.17006343e-6,
}

RMIN_NM = 0.3762
EPSILON_OVER_KB_K = 143.12
SIGMA_NM = 0.3357


def _damping(p: int, x: float) -> tuple[float, float, float]:
    """Return f_p(x), df_p/dr divided by b, and d2f_p/dr2 divided by b^2.

    The second and third outputs are derivatives with respect to x.
    """
    series = sum(x**k / math.factorial(k) for k in range(p + 1))
    ex = math.exp(-x)
    f = 1.0 - ex * series
    fp_x = ex * x**p / math.factorial(p)
    if x == 0.0:
        fpp_x = 0.0
    else:
        fpp_x = ex * x ** (p - 1) * (p - x) / math.factorial(p)
    return f, fp_x, fpp_x


def base_potential_K(r_nm: float) -> tuple[float, float, float]:
    """Return U/kB [K], d(U/kB)/dr [K/nm], d2(U/kB)/dr2 [K/nm^2]."""
    if r_nm <= 0.0:
        raise ValueError("r must be positive")

    r = r_nm
    phi = A1_NM1 * r + A2_NM2 * r*r + AM1_NM / r + AM2_NM2 / (r*r)
    rep = A_K * math.exp(phi)

    phi1 = A1_NM1 + 2.0*A2_NM2*r - AM1_NM/(r*r) - 2.0*AM2_NM2/(r**3)
    phi2 = 2.0*A2_NM2 + 2.0*AM1_NM/(r**3) + 6.0*AM2_NM2/(r**4)
    rep1 = rep * phi1
    rep2 = rep * (phi1*phi1 + phi2)

    disp = 0.0
    disp1 = 0.0
    disp2 = 0.0
    x = B_NM1 * r

    for p, c in C_K_NM.items():
        f, fx, fxx = _damping(p, x)
        fr = B_NM1 * fx
        frr = B_NM1**2 * fxx

        rp = r ** (-p)
        term = c * rp * f
        term1 = c * (-p * r ** (-p - 1) * f + rp * fr)
        term2 = c * (
            p * (p + 1) * r ** (-p - 2) * f
            - 2.0 * p * r ** (-p - 1) * fr
            + rp * frr
        )
        disp += term
        disp1 += term1
        disp2 += term2

    return rep - disp, rep1 - disp1, rep2 - disp2


def scaled_potential_K(r_nm: float, attraction_scale: float) -> tuple[float, float, float]:
    """Return WCA-like U_s/kB and its first two radial derivatives.

    The additive shift below r_min affects energy only. Force and curvature
    are the original TT/JHBV values below r_min and are multiplied by s above
    r_min.
    """
    s = float(attraction_scale)
    if s < 0.0:
        raise ValueError("attraction_scale must be non-negative")
    u, u1, u2 = base_potential_K(r_nm)
    if r_nm < RMIN_NM:
        return u + (1.0 - s) * EPSILON_OVER_KB_K, u1, u2
    return s*u, s*u1, s*u2


def scaled_potential_SI(r_m: float, attraction_scale: float) -> tuple[float, float, float]:
    """Return U [J], dU/dr [N], and d2U/dr2 [N/m]."""
    u, u1, u2 = scaled_potential_K(r_m / NM, attraction_scale)
    return KB*u, KB*u1/NM, KB*u2/(NM*NM)
