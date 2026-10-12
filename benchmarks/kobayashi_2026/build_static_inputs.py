#!/usr/bin/env python3
"""Build Hessian and shear/bulk affine-force inputs from a Kobayashi snapshot.

Example:
    python benchmarks/kobayashi_2026/build_static_inputs.py \
        /path/to/data_from_kobayashi --model s_1.00 --frame-index 0 \
        --outdir kobayashi_s1_frame0

By default the script writes the mass-normalized Hessian and the two
mass-weighted affine-force vectors as .npy files.  Full diagonalization is
optional because a 4116 x 4116 dense eigensolve is comparatively expensive.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from benchmarks.kobayashi_2026.static_nald import (
    affine_force,
    dense_mass_weighted_hessian,
    forces_and_pair_metrics,
    load_rxyz_frames,
    mass_weight_affine_force,
    source_box_length_m,
)
from benchmarks.kobayashi_2026.tt_jhbv import (
    ARGON_MASS_KG,
    source_pair_cutoff_m,
)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("data_dir", type=Path)
    p.add_argument("--model", choices=("s_1.00", "s_0.00"), required=True)
    p.add_argument("--frame-index", type=int, default=0)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--delta", type=float, default=1.0e-6)
    p.add_argument("--diagonalize", action="store_true")
    args = p.parse_args()

    s = 1.0 if args.model == "s_1.00" else 0.0
    frames = load_rxyz_frames(args.data_dir / args.model / "rxyz.dat")
    if not 0 <= args.frame_index < len(frames):
        raise SystemExit(f"frame-index must be in [0,{len(frames)-1}]")

    step, xyz = frames[args.frame_index]
    length = source_box_length_m(natoms=len(xyz))

    forces, pair_count, trace_mw = forces_and_pair_metrics(xyz, length, s)
    xi_bulk = affine_force(xyz, length, s, "bulk", args.delta)
    xi_shear = affine_force(xyz, length, s, "shear", args.delta)
    xi_bulk_mw = mass_weight_affine_force(xi_bulk)
    xi_shear_mw = mass_weight_affine_force(xi_shear)

    args.outdir.mkdir(parents=True, exist_ok=True)
    np.save(args.outdir / "affine_bulk_mass_weighted.npy", xi_bulk_mw)
    np.save(args.outdir / "affine_shear_mass_weighted.npy", xi_shear_mw)

    hess = dense_mass_weighted_hessian(xyz, length, s)
    np.save(args.outdir / "hessian_mass_weighted.npy", hess)

    metadata = {
        "model": args.model,
        "attraction_scale": s,
        "frame_index": args.frame_index,
        "step": step,
        "natoms": len(xyz),
        "box_length_m": length,
        "argon_mass_kg": ARGON_MASS_KG,
        "pair_cutoff_m": source_pair_cutoff_m(s),
        "pair_count": pair_count,
        "affine_delta": args.delta,
        "force_rms_component_N": float(np.sqrt(np.mean(forces**2))),
        "force_max_atom_N": float(np.max(np.linalg.norm(forces, axis=1))),
        "mass_weighted_hessian_trace_s-2": trace_mw,
        "bulk_affine_norm_N": float(np.linalg.norm(xi_bulk)),
        "shear_affine_norm_N": float(np.linalg.norm(xi_shear)),
        "bulk_affine_mass_weighted_norm": float(np.linalg.norm(xi_bulk_mw)),
        "shear_affine_mass_weighted_norm": float(np.linalg.norm(xi_shear_mw)),
    }

    if args.diagonalize:
        eigenvalues, eigenvectors = np.linalg.eigh(hess)
        np.save(args.outdir / "eigenvalues_s-2.npy", eigenvalues)
        np.save(args.outdir / "eigenvectors.npy", eigenvectors)

        gamma_bulk = (eigenvectors.T @ xi_bulk_mw) ** 2
        gamma_shear = (eigenvectors.T @ xi_shear_mw) ** 2
        np.save(args.outdir / "gamma_bulk.npy", gamma_bulk)
        np.save(args.outdir / "gamma_shear.npy", gamma_shear)

        metadata["negative_modes"] = int(np.sum(eigenvalues < 0.0))
        metadata["near_zero_modes_abs_lt_1e12_s-2"] = int(
            np.sum(np.abs(eigenvalues) < 1.0e12)
        )

    (args.outdir / "metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(json.dumps(metadata, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
