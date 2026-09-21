# LAMMPS affine-force templates

The files under `atomistic_real/` are templates for systems using LAMMPS
`units real`. They retain the force-field section of the atomistic epoxy
example so that the historical workflow can be reproduced, but **the
deformation logic is the reusable part**.

For a different material, replace the force-field, atom-style, data-file, and
electrostatics sections while keeping the same centered finite-difference
construction.

- `in.AF_shear`: xz shear, `Xi=-df/dgamma_xz`.
- `in.AF_bulk`: isotropic volume strain,
  `epsilon_v=Delta V/V`, `Xi=-df/d epsilon_v`.

Both files write atom ID, mass, and the three components of `Xi`. The mass is
required by the default mass-weighted projection in
`scripts/diagonalize_channels.py`.

For quantitative work, repeat each calculation at several strain amplitudes and
verify a strain-independent affine-force field within numerical noise.
