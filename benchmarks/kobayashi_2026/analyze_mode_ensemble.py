#!/usr/bin/env python3
"""Run the 10-snapshot Kobayashi mode-resolved benchmark locally.

The collaborator coordinate files are not required to live in this public
repository. Point this script at an authorized local copy containing

    s_1.00/rxyz.dat
    s_0.00/rxyz.dat

For every snapshot the script builds the source-matched dense mass-normalized
Hessian, computes bulk/shear affine-force vectors, diagonalizes the 4116x4116
matrix, and writes a cutoff scan of the structural prefactors

    A_K = (1/V) sum Gamma_K / lambda^2
    A_G = (1/V) sum Gamma_G / lambda^2.

The zero-frequency Markovian viscosities would be zeta=nu*A_K and eta=nu*A_G.
This script does NOT fit or assume nu.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from benchmarks.kobayashi_2026.mode_analysis import (
    branch_counts,
    null_gamma_fraction,
    signed_frequency_thz,
    structural_prefactor,
)
from benchmarks.kobayashi_2026.static_nald import (
    affine_force,
    dense_mass_weighted_hessian,
    load_rxyz_frames,
    mass_weight_affine_force,
    source_box_length_m,
)


DEFAULT_CUTOFFS = (
    0.0,
    0.005,
    0.010,
    0.015,
    0.020,
    0.025,
    0.030,
    0.040,
    0.050,
    0.075,
    0.100,
    0.125,
    0.150,
    0.175,
    0.200,
    0.250,
    0.300,
    0.400,
)


def _model_scale(model: str) -> float:
    return 1.0 if model == "s_1.00" else 0.0


def _key(cutoff: float) -> str:
    return f"{cutoff:.3f}".rstrip("0").rstrip(".").replace(".", "p")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("data_dir", type=Path)
    p.add_argument("--outdir", type=Path, default=Path("kobayashi_mode_results"))
    p.add_argument("--delta", type=float, default=1.0e-6)
    p.add_argument("--max-frames", type=int, default=None)
    args = p.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []

    for model in ("s_1.00", "s_0.00"):
        scale = _model_scale(model)
        frames = load_rxyz_frames(args.data_dir / model / "rxyz.dat")
        if args.max_frames is not None:
            frames = frames[: args.max_frames]

        length = source_box_length_m(natoms=len(frames[0][1]))
        volume = length**3

        for frame_index, (step, xyz) in enumerate(frames):
            xi_bulk = mass_weight_affine_force(
                affine_force(xyz, length, scale, "bulk", args.delta)
            )
            xi_shear = mass_weight_affine_force(
                affine_force(xyz, length, scale, "shear", args.delta)
            )
            hess = dense_mass_weighted_hessian(xyz, length, scale)
            eigenvalues, eigenvectors = np.linalg.eigh(hess)

            gamma_bulk = (eigenvectors.T @ xi_bulk) ** 2
            gamma_shear = (eigenvectors.T @ xi_shear) ** 2
            counts = branch_counts(eigenvalues)
            freq = np.abs(signed_frequency_thz(eigenvalues))
            nonnull = freq[eigenvalues != 0.0]

            row: dict[str, object] = {
                "model": model,
                "frame_index": frame_index,
                "step": step,
                "negative_modes": counts["negative"],
                "positive_modes": counts["positive"],
                "null_modes": counts["null"],
                "min_abs_frequency_THz": float(np.min(nonnull)),
                "gamma_bulk_null_fraction": null_gamma_fraction(
                    eigenvalues, gamma_bulk
                ),
                "gamma_shear_null_fraction": null_gamma_fraction(
                    eigenvalues, gamma_shear
                ),
            }

            for cutoff in DEFAULT_CUTOFFS:
                key = _key(cutoff)
                a_bulk = structural_prefactor(
                    eigenvalues, gamma_bulk, volume, cutoff_thz=cutoff
                )
                a_shear = structural_prefactor(
                    eigenvalues, gamma_shear, volume, cutoff_thz=cutoff
                )
                row[f"A_bulk_{key}_kg_m"] = a_bulk
                row[f"A_shear_{key}_kg_m"] = a_shear
                row[f"A_bulk_over_A_shear_{key}"] = a_bulk / a_shear

            rows.append(row)
            print(
                model,
                frame_index,
                step,
                counts,
                flush=True,
            )

    with (args.outdir / "frame_metrics.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
