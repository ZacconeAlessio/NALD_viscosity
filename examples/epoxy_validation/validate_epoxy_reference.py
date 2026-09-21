#!/usr/bin/env python3
"""Reproduce the historical T=300 K epoxy G_dp.f90 calculation."""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np

TU = 1.0e12
NU = 5.0e13
LEGACY_CONVERSION = 2.906e6
LAMBDA_CUTOFF = 2.57


def reproduce(reference_dir: Path) -> tuple[np.ndarray, float, int]:
    lam = np.loadtxt(reference_dir / "eigenvalues.data", dtype=float).reshape(-1)
    gamma = np.loadtxt(reference_dir / "gamma.data", dtype=float).reshape(-1)
    volume = float(np.loadtxt(reference_dir / "volume.dat", dtype=float))
    if lam.size != gamma.size:
        raise ValueError("eigenvalue and gamma lengths differ")

    gamma = gamma.copy()
    cut = (lam > -LAMBDA_CUTOFF) & (lam < LAMBDA_CUTOFF)
    gamma[cut] = 0.0

    rows = []
    for p in range(1, 236):
        omega = TU * np.exp(-15.0 + 0.1 * p)
        den = ((TU**2) * lam[:-1] - omega**2) ** 2 + (NU * omega) ** 2
        f = LEGACY_CONVERSION * omega * gamma[:-1] * NU / den
        gpp_gpa = np.sum(f) * (1.0e30 / volume) * 1.0e-9
        rows.append((omega / TU, gpp_gpa))

    active = ~cut[:-1]
    eta0 = (
        LEGACY_CONVERSION
        * (1.0e30 / volume)
        * np.sum(gamma[:-1][active] * NU / (((TU**2) * lam[:-1][active]) ** 2))
    )
    return np.asarray(rows), float(eta0), int(np.count_nonzero(cut))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("reference_dir", type=Path)
    p.add_argument("--rtol", type=float, default=2.0e-8)
    args = p.parse_args()

    calc, eta0, ncut = reproduce(args.reference_dir)
    ref = np.loadtxt(args.reference_dir / "G_dprime_T300.data", dtype=float)
    if calc.shape != ref.shape:
        raise ValueError(f"shape mismatch: calculated {calc.shape}, reference {ref.shape}")

    freq_err = np.max(np.abs(calc[:, 0] - ref[:, 0]))
    rel = np.abs(calc[:, 1] - ref[:, 1]) / np.maximum(np.abs(ref[:, 1]), 1.0e-300)
    max_rel = float(np.max(rel))

    print(f"historical cutoff removed {ncut} modes")
    print(f"max frequency absolute error: {freq_err:.6e}")
    print(f"max G'' relative error: {max_rel:.6e}")
    print(f"historical eta_s: {eta0:.12g} Pa s")
    if max_rel > args.rtol:
        raise SystemExit(f"FAIL: {max_rel:.3e} > {args.rtol:.3e}")
    print("PASS")


if __name__ == "__main__":
    main()
