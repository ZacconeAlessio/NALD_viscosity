#!/usr/bin/env python3
"""Reproduce the separately supplied direct epoxy zero-frequency formula."""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np

TU = 1.0e12
NU = 5.0e13
LEGACY_CONVERSION = 2.906e6


def legacy_eta(
    eigenvalues: np.ndarray,
    gamma: np.ndarray,
    volume_a3: float,
    cutoff: float = 1.0,
    zero_tol: float = 1.0e-5,
    omit_last_mode: bool = True,
) -> tuple[float, int]:
    lam = np.asarray(eigenvalues, dtype=float).reshape(-1)
    gam = np.asarray(gamma, dtype=float).reshape(-1).copy()
    if lam.size != gam.size:
        raise ValueError("eigenvalue and gamma lengths differ")

    gam[np.abs(lam) < zero_tol] = 0.0
    cut = (lam > -cutoff) & (lam < cutoff)
    gam[cut] = 0.0

    stop = lam.size - 1 if omit_last_mode else lam.size
    den = (TU**2 * lam[:stop]) ** 2
    eta_gpa_s = (
        np.sum(LEGACY_CONVERSION * gam[:stop] * NU / den)
        * (1.0e30 / volume_a3)
        * 1.0e-9
    )
    return float(eta_gpa_s * 1.0e9), int(np.count_nonzero(cut))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("reference_dir", type=Path)
    p.add_argument("--cutoff-frequency", type=float, default=1.0)
    p.add_argument("--zero-tol", type=float, default=1.0e-5)
    p.add_argument("--include-last-mode", action="store_true")
    args = p.parse_args()

    ref = args.reference_dir
    lam = np.loadtxt(ref / "eigenvalues.data", dtype=float)
    gamma = np.loadtxt(ref / "gamma.data", dtype=float)
    volume = float(np.loadtxt(ref / "volume.dat", dtype=float))

    eta_pa_s, ncut = legacy_eta(
        lam,
        gamma,
        volume,
        cutoff=args.cutoff_frequency,
        zero_tol=args.zero_tol,
        omit_last_mode=not args.include_last_mode,
    )
    print(f"cutoff |lambda| < {args.cutoff_frequency:g}: removed {ncut} modes")
    print(f"eta = {eta_pa_s:.12g} Pa s")


if __name__ == "__main__":
    main()
