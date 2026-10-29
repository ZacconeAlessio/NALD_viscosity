#!/usr/bin/env python3
"""Audit modes below the independent acoustic cutoff.

This is intentionally separate from the viscosity calculation.  It asks
whether modes removed by the hard q_min cutoff look like extended q_min
plane waves or like other instantaneous soft modes.

Example:
    python benchmarks/kobayashi_2026/analyze_low_modes.py \
        /path/to/data_from_kobayashi --frames 0,5,9 --out low_mode_audit.csv
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from benchmarks.kobayashi_2026.acoustic_finite_size import plane_wave_rayleigh
from benchmarks.kobayashi_2026.low_mode_character import (
    participation_ratio,
    qmin_plane_wave_weights,
    summarize_subcutoff_modes,
)
from benchmarks.kobayashi_2026.static_nald import (
    affine_force,
    dense_mass_weighted_hessian,
    load_rxyz_frames,
    mass_weight_affine_force,
    source_box_length_m,
)


def _parse_frames(spec: str, nframes: int) -> list[int]:
    if spec.strip().lower() == "all":
        return list(range(nframes))
    indices = sorted({int(x.strip()) for x in spec.split(",") if x.strip()})
    if any(i < 0 or i >= nframes for i in indices):
        raise ValueError(f"frame indices must be in [0,{nframes-1}]")
    return indices


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("data_dir", type=Path)
    p.add_argument("--frames", default="0,5,9")
    p.add_argument("--delta", type=float, default=1.0e-6)
    p.add_argument("--out", type=Path, default=Path("low_mode_audit.csv"))
    args = p.parse_args()

    rows = []
    for model, scale in (("s_1.00", 1.0), ("s_0.00", 0.0)):
        frames = load_rxyz_frames(args.data_dir / model / "rxyz.dat")
        length = source_box_length_m(natoms=len(frames[0][1]))
        volume = length**3

        for frame_index in _parse_frames(args.frames, len(frames)):
            step, xyz = frames[frame_index]
            acoustic = plane_wave_rayleigh(xyz, length, scale, 1)

            xi_bulk = mass_weight_affine_force(
                affine_force(xyz, length, scale, "bulk", args.delta)
            )
            xi_shear = mass_weight_affine_force(
                affine_force(xyz, length, scale, "shear", args.delta)
            )
            hessian = dense_mass_weighted_hessian(xyz, length, scale)
            eigenvalues, eigenvectors = np.linalg.eigh(hessian)

            gamma_bulk = (eigenvectors.T @ xi_bulk) ** 2
            gamma_shear = (eigenvectors.T @ xi_shear) ** 2
            pr = participation_ratio(eigenvectors, len(xyz))
            qtotal, _, _ = qmin_plane_wave_weights(xyz, eigenvectors, length)

            bulk = summarize_subcutoff_modes(
                eigenvalues,
                gamma_bulk,
                pr,
                qtotal,
                acoustic["f_longitudinal_THz"],
                volume,
            )
            shear = summarize_subcutoff_modes(
                eigenvalues,
                gamma_shear,
                pr,
                qtotal,
                acoustic["f_transverse_THz"],
                volume,
            )

            row = {
                "model": model,
                "frame_index": frame_index,
                "step": step,
                "f_L_qmin_THz": acoustic["f_longitudinal_THz"],
                "f_T_qmin_THz": acoustic["f_transverse_THz"],
            }
            row.update({f"bulk_{k}": v for k, v in bulk.items()})
            row.update({f"shear_{k}": v for k, v in shear.items()})
            rows.append(row)
            print(model, frame_index, step, flush=True)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
