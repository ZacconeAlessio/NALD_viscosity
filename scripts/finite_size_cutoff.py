#!/usr/bin/env python3
"""Evaluate the published finite-size shear cutoff."""

from __future__ import annotations

import argparse
from nald_viscosity.finite_size import shear_cutoff


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--box-length", type=float, required=True)
    p.add_argument("--density", type=float, required=True)
    p.add_argument("--shear-modulus", type=float, required=True)
    p.add_argument("--units", choices=("si", "reduced"), default="si")
    args = p.parse_args()

    omega = shear_cutoff(args.box_length, args.density, args.shear_modulus)
    if args.units == "si":
        print(f"omega_min = {omega:.12g} s^-1")
        print(f"omega_min = {omega/1.0e12:.12g} THz")
        print(f"lambda_cut = {(omega/1.0e12)**2:.12g} THz^2")
    else:
        print(f"omega_min = {omega:.12g} (reduced frequency units)")
        print(f"lambda_cut = {omega**2:.12g} (reduced eigenvalue units)")


if __name__ == "__main__":
    main()
