#!/usr/bin/env python3
"""Diagonalize one Hessian and project multiple affine-force channels."""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np

from nald_viscosity.modes import (
    affine_vector,
    gamma_from_eigenvectors,
    read_affine,
    read_hessian,
)


def parse_channel(spec: str) -> tuple[str, Path]:
    if "=" not in spec:
        raise argparse.ArgumentTypeError("use NAME=PATH, e.g. shear=AF_shear.data")
    name, path = spec.split("=", 1)
    if not name.strip():
        raise argparse.ArgumentTypeError("channel name cannot be empty")
    return name.strip(), Path(path)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hessian", type=Path, required=True)
    p.add_argument("--affine", action="append", type=parse_channel, required=True)
    p.add_argument("--skip-header", type=int, default=9)
    p.add_argument("--outdir", type=Path, default=Path("."))
    p.add_argument(
        "--legacy-unweighted-af",
        action="store_true",
        help="project raw affine forces; only for reproducing historical data",
    )
    args = p.parse_args()

    channels: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    natoms = None
    masses_ref = None

    for name, path in args.affine:
        masses, xi = read_affine(path, args.skip_header)
        if natoms is None:
            natoms = len(masses)
            masses_ref = masses
        elif len(masses) != natoms:
            raise ValueError("all affine-force files must contain the same atoms")
        elif not np.allclose(masses, masses_ref, rtol=1e-12, atol=1e-12):
            raise ValueError("atomic masses differ among affine-force channels")
        channels[name] = (masses, xi)

    assert natoms is not None
    print(f"Reading {3*natoms} x {3*natoms} Hessian for N={natoms}...")
    hessian = read_hessian(args.hessian, natoms)
    print("Diagonalizing...")
    eigenvalues, eigenvectors = np.linalg.eigh(hessian)

    args.outdir.mkdir(parents=True, exist_ok=True)
    np.savetxt(args.outdir / "eigenvalues.data", eigenvalues)

    for name, (masses, xi) in channels.items():
        mass_weighted = not args.legacy_unweighted_af
        xi_vec = affine_vector(masses, xi, mass_weighted=mass_weighted)
        gamma = gamma_from_eigenvectors(eigenvectors, xi_vec)
        convention = "massweighted" if mass_weighted else "legacy-unweighted"
        np.savetxt(args.outdir / f"gamma_{name}.data", gamma)
        np.savetxt(
            args.outdir / f"eigen_gamma_{name}.data",
            np.column_stack((eigenvalues, gamma)),
            header=f"eigenvalue Gamma_{name}; affine_projection={convention}",
        )
        print(f"Wrote {name} channel ({convention}).")


if __name__ == "__main__":
    main()
