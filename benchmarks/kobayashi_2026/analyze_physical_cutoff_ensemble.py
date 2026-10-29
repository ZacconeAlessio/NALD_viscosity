#!/usr/bin/env python3
"""Recompute the physical-cutoff Kobayashi ensemble under the current code.

This is the authoritative generator for frame-level and ensemble summaries
that use the independently estimated longitudinal q_min cutoff for bulk and
the transverse q_min cutoff for shear.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from benchmarks.kobayashi_2026.acoustic_finite_size import plane_wave_rayleigh
from benchmarks.kobayashi_2026.mode_analysis import (
    branch_counts,
    numerical_null_mask,
    signed_frequency_thz,
)
from benchmarks.kobayashi_2026.static_nald import (
    affine_force,
    dense_mass_weighted_hessian,
    load_rxyz_frames,
    mass_weight_affine_force,
    source_box_length_m,
)

PUBLISHED = {
    "s_1.00": {"zeta": 3.79e-4, "eta": 0.958e-4, "kdiff": 0.673e9},
    "s_0.00": {"zeta": 0.344e-4, "eta": 0.790e-4, "kdiff": 0.3908e9},
}


def _structural(eigenvalues, gamma, volume, cutoff_thz):
    freq = np.abs(signed_frequency_thz(eigenvalues))
    keep = (~numerical_null_mask(eigenvalues)) & (freq >= cutoff_thz)
    a = float(np.sum(gamma[keep] / eigenvalues[keep] ** 2) / volume)
    d = float(np.sum(gamma[keep] / eigenvalues[keep]) / volume)
    return a, d, int(np.count_nonzero(keep))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("data_dir", type=Path)
    p.add_argument("--outdir", type=Path, default=Path("physical_cutoff_results"))
    p.add_argument("--delta", type=float, default=1.0e-6)
    args = p.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    for model, scale in (("s_1.00", 1.0), ("s_0.00", 0.0)):
        frames = load_rxyz_frames(args.data_dir / model / "rxyz.dat")
        length = source_box_length_m(natoms=len(frames[0][1]))
        volume = length**3

        for frame_index, (step, xyz) in enumerate(frames):
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

            a_bulk, d_bulk, n_bulk = _structural(
                eigenvalues, gamma_bulk, volume, acoustic["f_longitudinal_THz"]
            )
            a_shear, d_shear, n_shear = _structural(
                eigenvalues, gamma_shear, volume, acoustic["f_transverse_THz"]
            )
            counts = branch_counts(eigenvalues)

            rows.append(
                {
                    "model": model,
                    "frame_index": frame_index,
                    "step": step,
                    "L_nm": length * 1e9,
                    "f_L_qmin_THz": acoustic["f_longitudinal_THz"],
                    "f_T_qmin_THz": acoustic["f_transverse_THz"],
                    "f_L_axis_std_THz": acoustic["f_longitudinal_std_THz"],
                    "f_T_pol_std_THz": acoustic["f_transverse_std_THz"],
                    "longitudinal_projection_mean": acoustic[
                        "longitudinal_projection_mean"
                    ],
                    "longitudinal_projection_min": acoustic[
                        "longitudinal_projection_min"
                    ],
                    "negative_modes": counts["negative"],
                    "positive_modes": counts["positive"],
                    "null_modes": counts["null"],
                    "A_bulk_phys_kg_m": a_bulk,
                    "D_bulk_phys_Pa": d_bulk,
                    "A_shear_phys_kg_m": a_shear,
                    "D_shear_phys_Pa": d_shear,
                    "modes_bulk_phys": n_bulk,
                    "modes_shear_phys": n_shear,
                }
            )
            print(model, frame_index, step, flush=True)

    frame_path = args.outdir / "physical_cutoff_frame_metrics.csv"
    with frame_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    summary = []
    for model in ("s_1.00", "s_0.00"):
        mr = [r for r in rows if r["model"] == model]
        pub = PUBLISHED[model]

        def stats(key):
            x = np.asarray([float(r[key]) for r in mr])
            return float(x.mean()), float(x.std(ddof=1))

        f_l, f_l_sd = stats("f_L_qmin_THz")
        f_t, f_t_sd = stats("f_T_qmin_THz")
        a_b, a_b_sd = stats("A_bulk_phys_kg_m")
        d_b, d_b_sd = stats("D_bulk_phys_Pa")
        a_s, a_s_sd = stats("A_shear_phys_kg_m")
        d_s, d_s_sd = stats("D_shear_phys_Pa")
        neg, neg_sd = stats("negative_modes")
        pos, pos_sd = stats("positive_modes")
        nul, nul_sd = stats("null_modes")

        summary.append(
            {
                "model": model,
                "nframes": len(mr),
                "f_L_qmin_THz_mean": f_l,
                "f_L_qmin_THz_std": f_l_sd,
                "f_T_qmin_THz_mean": f_t,
                "f_T_qmin_THz_std": f_t_sd,
                "A_bulk_phys_kg_m_mean": a_b,
                "A_bulk_phys_kg_m_std": a_b_sd,
                "D_bulk_phys_Pa_mean": d_b,
                "D_bulk_phys_Pa_std": d_b_sd,
                "A_shear_phys_kg_m_mean": a_s,
                "A_shear_phys_kg_m_std": a_s_sd,
                "D_shear_phys_Pa_mean": d_s,
                "D_shear_phys_Pa_std": d_s_sd,
                "negative_modes_mean": neg,
                "negative_modes_std": neg_sd,
                "positive_modes_mean": pos,
                "positive_modes_std": pos_sd,
                "null_modes_mean": nul,
                "null_modes_std": nul_sd,
                "published_zeta_Pa_s": pub["zeta"],
                "published_eta_Pa_s": pub["eta"],
                "published_Kdiff_Pa": pub["kdiff"],
                "nu_bulk_required_s-1": pub["zeta"] / a_b,
                "nu_shear_required_s-1": pub["eta"] / a_s,
                "D_bulk_over_published_Kdiff": d_b / pub["kdiff"],
            }
        )

    summary_path = args.outdir / "physical_cutoff_summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0].keys()))
        writer.writeheader()
        writer.writerows(summary)


if __name__ == "__main__":
    main()
