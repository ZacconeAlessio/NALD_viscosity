#!/usr/bin/env python3
"""Validate a locally supplied Kobayashi 124 K benchmark package.

The collaborator-supplied data are intentionally NOT stored in this public
repository. Point this script to an unpacked directory containing

    s_1.00/bulk_corr.dat
    s_1.00/rxyz.dat
    s_0.00/bulk_corr.dat
    s_0.00/rxyz.dat

The bulk-correlation second column has units Pa and numerically corresponds to
the Green-Kubo integrand V <dP(0)dP(t)>/(k_B T). Its value at t=0 is therefore
K_inf-K_0, and integration over time gives zeta in Pa s.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np

N_ATOMS = 1372
DENSITY_KG_M3 = 1124.9
ARGON_MASS_U = 39.948
AMU_KG = 1.66053906660e-27

PUBLISHED = {
    "s_1.00": {
        "zeta_pa_s": 3.79e-4,
        "zeta_stderr_pa_s": 0.08e-4,
        "kdiff_pa": 0.673e9,
    },
    "s_0.00": {
        "zeta_pa_s": 0.344e-4,
        "zeta_stderr_pa_s": 0.009e-4,
        "kdiff_pa": 0.3908e9,
    },
}


def load_frames(path: Path, natoms: int = N_ATOMS) -> list[tuple[int, np.ndarray]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    stride = natoms + 1
    if len(lines) % stride:
        raise ValueError(f"{path}: line count is not a multiple of {stride}")
    frames = []
    for start in range(0, len(lines), stride):
        step = int(lines[start].strip())
        xyz = np.loadtxt(lines[start + 1:start + 1 + natoms], dtype=float)
        if xyz.shape != (natoms, 3):
            raise ValueError(f"{path}: bad coordinate block at step {step}")
        frames.append((step, xyz))
    return frames


def inferred_box_length_m(
    natoms: int = N_ATOMS,
    density_kg_m3: float = DENSITY_KG_M3,
    argon_mass_u: float = ARGON_MASS_U,
) -> float:
    volume = natoms * argon_mass_u * AMU_KG / density_kg_m3
    return volume ** (1.0/3.0)


def validate_state(root: Path, name: str) -> None:
    corr = np.loadtxt(root / name / "bulk_corr.dat", dtype=float)
    if corr.ndim != 2 or corr.shape[1] != 2:
        raise ValueError(f"{name}: bulk_corr.dat must have two columns")
    if not np.all(np.diff(corr[:, 0]) > 0):
        raise ValueError(f"{name}: time column must be strictly increasing")

    zeta = float(np.sum(0.5 * (corr[1:, 1] + corr[:-1, 1]) * np.diff(corr[:, 0])))
    k0 = float(corr[0, 1])
    pub = PUBLISHED[name]

    frames = load_frames(root / name / "rxyz.dat")
    steps = np.asarray([step for step, _ in frames], dtype=int)
    if len(frames) > 1:
        step_interval = int(np.median(np.diff(steps)))
    else:
        step_interval = 0

    box = inferred_box_length_m()
    max_coord = max(float(xyz.max()) for _, xyz in frames)
    min_coord = min(float(xyz.min()) for _, xyz in frames)

    print(f"[{name}]")
    print(f"correlation rows: {len(corr)}; final time: {corr[-1,0]*1e12:.6g} ps")
    print(f"C(0): {k0/1e9:.9g} GPa; published K_inf-K_0: {pub['kdiff_pa']/1e9:.9g} GPa")
    print(f"integral zeta: {zeta:.12g} Pa s ({zeta*1e3:.9g} mPa s)")
    print(f"published zeta: {pub['zeta_pa_s']:.12g} +/- {pub['zeta_stderr_pa_s']:.2g} Pa s")
    print(f"difference/published: {(zeta/pub['zeta_pa_s']-1.0)*100:.3f}%")
    print(f"coordinate frames: {len(frames)}; steps {steps[0]}..{steps[-1]}; interval {step_interval}")
    print(f"inferred cubic L (m_Ar={ARGON_MASS_U} u): {box*1e9:.9g} nm")
    print(f"coordinate range: {min_coord*1e9:.6g}..{max_coord*1e9:.6g} nm")
    if max_coord >= box * 1.001 or min_coord < -1e-15:
        raise ValueError(f"{name}: coordinates appear inconsistent with inferred box")
    print()


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("data_dir", type=Path)
    args = p.parse_args()

    for state in ("s_1.00", "s_0.00"):
        validate_state(args.data_dir, state)


if __name__ == "__main__":
    main()
