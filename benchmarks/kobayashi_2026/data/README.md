# Collaborator-supplied Kobayashi benchmark data

These data originate from H. Kobayashi's calculations for the 124 K,
1124.9 kg/m^3 liquid-argon state point in:

H. Kobayashi, Y. Ishii, and N. Ohtori,
"Structural origin of the strong effect of attraction on bulk viscosity in
simple liquids," *J. Chem. Phys.* **165**, 104506 (2026).
DOI: 10.1063/5.0351887.

H. Kobayashi gave permission in October 2026 for benchmark data supplied to
this project to be included in the public NALD_viscosity repository, with the
paper cited.

## Files

- `benchmark_summary.csv` records properties of the full collaborator-supplied
  bulk-correlation files and the viscosities obtained by direct trapezoidal
  integration of those full files.
- `bulk_corr_sampled.csv` is a compact exact-row sample of the supplied
  correlation functions for `s=1.00` and `s=0.00`. The
  `cumulative_zeta_Pa_s` column was evaluated from the **full-resolution**
  supplied data up to that original row; do not re-integrate only the sparse
  sampled rows and expect the same value.

The supplied correlation second column is treated here as the Green-Kubo bulk
integrand

    V <delta P(0) delta P(t)> / (k_B T)

in Pa. This interpretation is numerically confirmed because its t=0 value
matches the reported `K_inf-K_0`, while integration over time gives the
reported bulk viscosity.

For the full collaborator-supplied files:

- `s=1.00`: zeta = 0.380497267485 mPa s; published 0.379 +/- 0.008 mPa s.
- `s=0.00`: zeta = 0.035109834144 mPa s; published 0.0344 +/- 0.0009 mPa s.

The full package also contains ten 1372-particle coordinate snapshots for each
interaction model. Those configurations are being used for the NALD
Hessian/affine-force benchmark. A complete reproduction workflow will be added
after the collaborator's full MD program and accompanying manual have been
checked, so that the box, cutoff and force conventions can be mirrored exactly.

The journal PDFs and collaborator's Fortran source are **not redistributed**
here.

See `DATA_NOTICE.md` for reuse/provenance information.
