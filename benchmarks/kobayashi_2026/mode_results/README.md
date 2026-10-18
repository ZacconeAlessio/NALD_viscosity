# First mode-resolved Kobayashi benchmark results

These are **exploratory structural results, not yet a validated viscosity
prediction**.

The calculation uses all ten supplied 124 K snapshots for each pair model
(`s=1.00` and `s=0.00`), the source-matched TT/JHBV force law, the
mass-normalized 4116 x 4116 instantaneous Hessian, and the bulk/shear affine
forces evaluated on each same snapshot.

For the Markovian zero-frequency form it is useful to separate the structural
factor

```text
A_M = (1/V) sum_p Gamma_M,p / lambda_p^2
eta_M = nu A_M
```

from the damping parameter `nu`.  Both stable and unstable instantaneous
normal modes are retained. Numerical null modes are removed separately from
any physical low-frequency cutoff.

## Numerical checks

On the first supplied snapshot of each model, the analytic central-force
Hessian was checked against a centered finite difference of the force along a
random non-translational displacement.  The relative discrepancy is about
1e-8.  The bulk and shear affine-force vectors at strain step 1e-6 also agree
with a 3e-7 calculation to about 1e-8 relative accuracy.  See
`numerical_checks.csv`.

## Mode counts


Averaged over ten snapshots:

| model | negative modes | positive modes | numerical null modes |
| --- | ---: | ---: | ---: |
| s=1.00 | 1448.8 +/- 23.7 | 2664.2 +/- 23.7 | 3.0 |
| s=0.00 | 2224.7 +/- 21.9 | 1609.0 +/- 14.0 | 282.3 +/- 14.3 |

The three `s=1.00` null modes are the expected translations.  The
purely-repulsive instantaneous network contains a much larger floppy/null
subspace because only contacts inside `r_min` carry force and curvature.

Importantly, the null subspace is essentially orthogonal to the affine-force
fields: the mean fraction of total `Gamma_bulk` in the numerical nullspace is
about 9.7e-24 for `s=1.00` and 7.9e-25 for `s=0.00`.  Thus these modes can
be removed as numerical/floppy zero modes without creating the large
zero-frequency response seen from genuinely soft, coupled modes.

## Bulk attraction sensitivity

Without a physical low-frequency regularization, `A_K` is extremely
sensitive to rare very-soft modes and has large snapshot-to-snapshot variance.
We therefore scanned a **symmetric exploratory cutoff** in

```text
|f_p| = sqrt(|lambda_p|)/(2*pi).
```

The published 124 K bulk-viscosity ratio is

```text
zeta(s=1.00) / zeta(s=0.00) = 11.017.
```

At a 0.100 THz symmetric cutoff, the ten-snapshot structural ratio is

```text
A_K(s=1.00) / A_K(s=0.00) = 10.424.
```

A 5000-resample snapshot bootstrap gives a 95% interval of approximately
9.58--11.24 for that ratio.  Interpolating the scan, the structural ratio
crosses the published value near **0.095 THz**; the bootstrap interval for the
crossing is approximately **0.0887--0.1024 THz**.

This is encouraging: the Hessian plus bulk affine-force coupling contains a
strong attraction dependence on the same order as the Green-Kubo bulk
viscosity.  It is **not** yet a prediction, because the cutoff was identified
from the comparison rather than fixed independently.

## Why this is not yet a complete validation

The same 0.100 THz cutoff gives

```text
A_G(s=1.00) / A_G(s=0.00) = 5.368,
```

whereas the published shear-viscosity ratio is only 1.213.  Likewise, a single
channel-independent constant `nu` does not reproduce both `zeta/eta`
ratios.

The next scientific step is therefore to determine the low-frequency
regularization independently (for example from the appropriate finite-size
longitudinal/transverse scale) and to test a frequency-dependent memory kernel
rather than fitting a common Markovian damping constant.  This is especially
relevant for bulk viscosity because the Kobayashi paper identifies slow
pressure relaxation as the dominant attraction-sensitive mechanism.

## Files

- `mode_count_summary.csv`: ensemble mode/null-space statistics.
- `numerical_checks.csv`: finite-difference Hessian and affine-force checks.
- `cross_model_predictions.csv`: bulk and shear attraction-ratio scan.
- `bootstrap_cutoff_summary.csv`: bootstrap sensitivity of the bulk crossing.
- `../analyze_mode_ensemble.py`: reproducible local analysis from the
  collaborator-authorized coordinate package.

All numerical values here should be treated as research-stage benchmark output.
