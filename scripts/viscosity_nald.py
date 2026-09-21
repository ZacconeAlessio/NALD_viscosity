#!/usr/bin/env python3
"""Compute shear or bulk NALD loss modulus and zero-frequency viscosity."""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np

from nald_viscosity.core import (
    THZ_SI,
    convert_to_si,
    load_modes,
    loss_modulus_si,
    select_modes,
    zero_frequency_viscosity_si,
)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--eigenvalues", required=True)
    p.add_argument("--gamma", required=True)
    p.add_argument("--channel", choices=("shear", "bulk"), required=True)
    p.add_argument("--volume", type=float, required=True)
    p.add_argument("--nu", type=float, required=True, help="Markovian damping rate [s^-1]")
    p.add_argument(
        "--input-units",
        choices=("si", "lammps-real-massweighted"),
        default="lammps-real-massweighted",
    )
    p.add_argument("--zero-tol", type=float, default=1.0e-10)
    p.add_argument(
        "--cutoff-frequency",
        type=float,
        default=0.0,
        help="optional cutoff on sqrt(|lambda|) in native frequency units",
    )
    p.add_argument("--wmin", type=float, default=1.0e-6)
    p.add_argument("--wmax", type=float, default=1.0e2)
    p.add_argument("--nfrequency", type=int, default=300)
    p.add_argument("--output")
    return p


def main() -> None:
    args = build_parser().parse_args()
    if args.volume <= 0 or args.nu <= 0 or args.wmin <= 0 or args.wmax <= args.wmin:
        raise ValueError("positive volume/nu/frequencies required, with wmax>wmin")
    if args.nfrequency < 2:
        raise ValueError("nfrequency must be at least 2")

    lam_native, gamma = load_modes(args.eigenvalues, args.gamma)
    lam_native, gamma, info = select_modes(
        lam_native, gamma, args.zero_tol, args.cutoff_frequency
    )
    lam_si, gamma_si, volume_m3 = convert_to_si(
        lam_native, gamma, args.volume, args.input_units
    )

    omega_native = np.geomspace(args.wmin, args.wmax, args.nfrequency)
    if args.input_units == "si":
        omega_si = omega_native
        frequency_label = "omega_s^-1"
        modulus_scale = 1.0
        modulus_label = "loss_modulus_Pa"
    else:
        omega_si = omega_native * THZ_SI
        frequency_label = "omega_THz_repo_convention"
        modulus_scale = 1.0e-9
        modulus_label = "loss_modulus_GPa"

    loss_pa = loss_modulus_si(omega_si, lam_si, gamma_si, volume_m3, args.nu)
    viscosity_omega = loss_pa / omega_si
    viscosity_zero = zero_frequency_viscosity_si(
        lam_si, gamma_si, volume_m3, args.nu
    )

    viscosity_label = "eta_s_Pa_s" if args.channel == "shear" else "zeta_Pa_s"
    output = Path(args.output or f"{args.channel}_viscosity.data")
    header = (
        f"{frequency_label} {modulus_label} {viscosity_label}\n"
        f"channel={args.channel}; input_units={args.input_units}; "
        f"nu_s^-1={args.nu:.16g}; zero_tol={args.zero_tol:.16g}; "
        f"cutoff_frequency={args.cutoff_frequency:.16g}; "
        f"modes_total={info['total']}; modes_kept={info['kept']}; "
        f"removed_zero={info['removed_zero']}; "
        f"removed_cutoff={info['removed_cutoff']}; "
        f"zero_frequency_viscosity_Pa_s={viscosity_zero:.16g}"
    )
    np.savetxt(
        output,
        np.column_stack((omega_native, loss_pa * modulus_scale, viscosity_omega)),
        header=header,
    )

    print(f"Wrote {output}")
    print(
        f"Modes: {info['kept']}/{info['total']} kept "
        f"({info['removed_zero']} zero; {info['removed_cutoff']} cutoff)"
    )
    print(f"Zero-frequency {viscosity_label}: {viscosity_zero:.12g}")
    if np.any(lam_native < 0):
        print("Negative instantaneous-normal-mode eigenvalues retained.")


if __name__ == "__main__":
    main()
