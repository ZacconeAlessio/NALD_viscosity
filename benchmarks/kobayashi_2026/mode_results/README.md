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

## Exploratory cutoff scan

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

This fitted crossing is useful only as a sensitivity diagnostic. It is **not**
a prediction, because the cutoff was identified from the Green-Kubo comparison.

## Independent acoustic finite-size scale

We next estimated the fundamental propagating scale without using viscosity
data. For each snapshot, longitudinal and transverse plane-wave Rayleigh
quotients of the source-matched instantaneous Hessian were evaluated at

```text
q_min = 2*pi/L = 1.453 nm^-1.
```

This is essentially the same smallest wavevector (~1.5 nm^-1) emphasized in
the Kobayashi paper for the low-q pressure-density coupling.

Averaged over ten snapshots:

| model | f_L(q_min) (THz) | f_T(q_min) (THz) |
| --- | ---: | ---: |
| s=1.00 | 0.2746 +/- 0.0031 | 0.1619 +/- 0.0022 |
| s=0.00 | 0.2817 +/- 0.0015 | 0.1495 +/- 0.0020 |

Therefore the earlier fitted crossing near 0.095 THz lies **below** both the
independently estimated transverse and longitudinal fundamental propagating
frequencies. It cannot be justified as the simple q_min acoustic cutoff of this
box.

Using the model-specific longitudinal q_min scale for the bulk channel and the
transverse q_min scale for shear gives

```text
A_K(s=1.00) / A_K(s=0.00) = 2.191
A_G(s=1.00) / A_G(s=0.00) = 3.174
```

compared with the published Green-Kubo ratios

```text
zeta(s=1.00) / zeta(s=0.00) = 11.017
eta(s=1.00)  / eta(s=0.00)  = 1.213.
```

Thus the static Hessian plus affine-force coupling, with an independently
estimated finite-size scale, does **not** by itself reproduce the attraction
dependence of either transport coefficient.

If each model is forced into the constant-Markovian form eta_M=nu*A_M, the
required damping ratios are

```text
nu_bulk(s=1.00) / nu_bulk(s=0.00) = 5.03
nu_shear(s=1.00) / nu_shear(s=0.00) = 0.382.
```

The opposite trends in the two channels rule out a single transferable scalar
damping constant as a satisfactory explanation.

## What lies below the acoustic cutoff?

The hard q_min cutoff was audited on frames 0, 5, and 9 of both models using
two descriptors that do not enter the viscosity fit: atom-wise participation
ratio and overlap with the q_min plane-wave subspace.

For the bulk subcutoff structural weight, the three-frame averages are:

| model | A-weighted participation ratio | A-weighted q_min plane-wave weight |
| --- | ---: | ---: |
| s=1.00 | 0.274 +/- 0.049 | 0.0352 +/- 0.0051 |
| s=0.00 | 0.0573 +/- 0.0073 | 0.0348 +/- 0.0033 |

The q_min overlap is essentially unchanged while the spatial character is very
different.  In s=1.00, very soft extended instantaneous modes can dominate the
unregularized structural prefactor; in s=0.00 the corresponding subcutoff
weight is distributed over much more localized modes.

This does not license retaining all subcutoff modes: instantaneous soft modes
can be transient, so their lifetime belongs to the memory problem.  It does,
however, show that a frequency-only cutoff is not a complete mode
classification. See `low_mode_audit_sample.csv`.

## Pressure-relaxation target

The supplied Green-Kubo correlations independently show where the missing bulk
physics resides. Direct integration gives pressure-relaxation times of
0.565 ps (s=1.00) and 0.0898 ps (s=0.00). More strikingly, 90% of the
s=1.00 viscosity integral is accumulated only after 7.80 ps, versus 0.120 ps
for s=0.00.

The published Gaussian-plus-double-exponential fit is consistent with this:
for s=1.00 it contains an intermediate exponential (tau2=0.43 ps) and a
slow component (tau3=3.69 ps), whereas s=0.00 is described by the Gaussian
component alone. This agrees with the paper's conclusion that attraction
creates slow low-q pressure/density relaxation that is absent from shear.

The benchmark therefore points to the next theory step: a non-Markovian and
possibly channel-dependent memory treatment, tested against the supplied
pressure-relaxation data rather than tuned only to the final viscosity.


## Second-order pressure-memory diagnostic

To move beyond a scalar Markovian damping constant, the supplied pressure
correlations were also analyzed with a second-order Mori form,

```text
C''(t) + omega0^2 C(t)
  + integral_0^t nu(t-t') C'(t') dt' = 0.
```

For the normalized correlation, the Laplace-domain memory is

```text
nu_tilde(s) =
  [s - (s^2 + omega0^2) C_tilde(s)] / [s C_tilde(s) - 1].
```

The short-time curvature gives `omega0^2=-C''(0)`.  Using a 0.04 ps even
quartic fit to the raw supplied correlations gives:

| model | bare frequency (THz) | tau_pressure (ps) | nu_tilde(0) (ps^-1) |
| --- | ---: | ---: | ---: |
| s=1.00 | 1.649 | 0.5655 | 60.73 |
| s=0.00 | 2.238 | 0.08984 | 17.77 |

Thus the zero-frequency pressure-memory scale changes by a factor of 3.418
between the two interaction models. Combining that dynamical ratio with the
independently regularized bulk structural ratio, 2.191, gives

```text
2.191 * 3.418 = 7.49,
```

which moves substantially toward the Green-Kubo bulk-viscosity ratio 11.017
but does not fully reproduce it.

This is a **diagnostic, not a prediction**: the pressure-memory information is
derived from the same pressure-correlation data used by Green-Kubo. Its value is
that it identifies the sign, scale, and frequency dependence of the missing
dynamical contribution. The remaining model-dependent factor differs by about
1.47 between s=1.00 and s=0.00.

The Laplace scan also shows that the attraction effect is strongest at long
times: `nu_tilde_s1/nu_tilde_s0` is about 3.42 at s=0 and approaches unity
as the Laplace rate increases. This is qualitatively consistent with the slow
pressure-relaxation tail being the attraction-sensitive part of the dynamics.

## Why this is not yet a complete validation

The independent acoustic calculation resolves the finite-size question but
makes the remaining limitation clearer: a constant Markovian damping model is
not transferable between the two interaction models or between bulk and shear.

The next scientific step is to formulate and test a frequency-dependent memory
kernel and to examine whether the low-q density-relaxation times supplied by
Kobayashi can constrain that kernel without fitting the final viscosity.

## Files

- `mode_count_summary.csv`: ensemble mode/null-space statistics.
- `numerical_checks.csv`: finite-difference Hessian and affine-force checks.
- `cross_model_predictions.csv`: bulk and shear attraction-ratio scan.
- `bootstrap_cutoff_summary.csv`: bootstrap sensitivity of the exploratory crossing.
- `acoustic_qmin_summary.csv`: independent q_min longitudinal/transverse scales.
- `physical_cutoff_summary.csv`: structural prefactors using those q_min scales.
- `low_mode_audit_sample.csv`: verified frames 0/5/9 participation-ratio and q_min-overlap audit.
- `pressure_relaxation_summary.csv`: direct integral/tail diagnostics from supplied correlations.
- `pressure_memory_summary.csv`: second-order Mori zero-frequency memory diagnostic.
- `pressure_memory_laplace.csv`: real-Laplace pressure-memory scan.
- `memory_nald_bridge.csv`: structural x pressure-memory comparison with Green-Kubo.
- `published_relaxation_fit_endpoints.csv`: published s=1 and s=0 relaxation-fit endpoints.
- `../analyze_mode_ensemble.py`: mode-resolved local analysis.
- `../analyze_acoustic_ensemble.py`: independent plane-wave q_min analysis.
- `../analyze_physical_cutoff_ensemble.py`: authoritative physical-cutoff regeneration.
- `../analyze_low_modes.py`: participation-ratio/q_min subcutoff mode audit.
- `../analyze_pressure_relaxation.py`: correlation-tail analysis.

All numerical values here should be treated as research-stage benchmark output.
