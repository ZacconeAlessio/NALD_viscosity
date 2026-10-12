# Kobayashi et al. JCP 2026 benchmark

Reference: H. Kobayashi, Y. Ishii, and N. Ohtori,
"Structural origin of the strong effect of attraction on bulk viscosity in
simple liquids," *J. Chem. Phys.* **165**, 104506 (2026).
DOI: 10.1063/5.0351887.

This directory provides the liquid-argon target for quantitative validation of
the NALD bulk-viscosity channel.

`reference.csv` transcribes the published bulk viscosity `zeta`, shear
viscosity `eta`, `K_inf-K_0`, and pressure-relaxation time `tau_zeta`
from Table I and Table S1. Published viscosity values in those tables use units
of `10^-1 mPa s`.

## Why this benchmark is useful

At 124 K and 1124.9 kg/m^3, reducing the pair-attraction scale from `s=1.00`
to `s=0.00` changes the published quantities approximately as follows:

- bulk viscosity: 3.79 -> 0.344 in units of 10^-1 mPa s (about -91%);
- shear viscosity: 0.958 -> 0.790 (about -18%);
- `K_inf-K_0`: 0.673 -> 0.3908 GPa (about -42%);
- `tau_zeta`: 0.56 -> 0.088 ps (about -84%).

The strong contrast between bulk and shear channels makes this a stringent test
of whether NALD modal couplings distinguish attraction-sensitive pressure
relaxation from shear relaxation.

## Collaborator-supplied 124 K benchmark

H. Kobayashi subsequently supplied benchmark material for the 124 K,
1124.9 kg/m^3 state point:

- ten equilibrated 1372-particle coordinate snapshots for `s=1.00`;
- ten equilibrated coordinate snapshots for `s=0.00`;
- full bulk Green-Kubo correlation data for both models;
- the TT/JHBV force routine used in the in-house MD code;
- references containing the potential parameters.

H. Kobayashi explicitly gave permission in October 2026 for benchmark data to
be included in this public repository, provided the JCP paper is cited.

A compact public subset is in [`data/`](data/). It contains a summary of the
full supplied correlation files plus exact sampled rows and the cumulative
full-resolution Green-Kubo integral. The full coordinate package is being used
for the NALD calculation.

Kobayashi has now also supplied the complete in-house MD package and a manual.
We checked the force, box, cutoff, periodic-boundary and run conventions needed
for the static benchmark. They are summarized in
[`source_protocol.md`](source_protocol.md); the collaborator source itself is
not redistributed.

Direct trapezoidal integration of the full supplied correlation files gives:

| model | supplied integral (mPa s) | published (mPa s) |
| --- | ---: | ---: |
| s=1.00 | 0.3804973 | 0.379 +/- 0.008 |
| s=0.00 | 0.0351098 | 0.0344 +/- 0.0009 |

The t=0 correlation values are 0.6729 and 0.3908 GPa, respectively, matching
the reported `K_inf-K_0` values.

Run

```bash
python benchmarks/kobayashi_2026/validate_public_data.py
```

to verify the public benchmark subset.

## TT/JHBV implementation

`tt_jhbv.py` now mirrors the supplied production source at force/Hessian
level. Parameters follow the cited Jager and Alboul-Lishchuk tables, while
`C12`, `C14`, and `C16` are generated with the same recurrence used in
the Fortran code.

The attraction-scaled `U_{2,s}` model keeps the full TT force below
`r_min = 3.7618 Angstrom` and scales the branch above the minimum by `s`.
The exact production cutoff is also reproduced: `s=0` uses `r_min`, while
every nonzero `s` uses `5 sigma_MD = 1.6976 nm`.

## Reported simulation protocol

The production system contains 1372 Ar atoms. The three state points reported
in the paper are (140 K, 968 kg/m^3), (124 K, 1124.9 kg/m^3), and
(90 K, 1390 kg/m^3).

The paper reports 2.5e7 NVE production steps, normally with a 6.45 fs timestep;
3 fs is used for the full-pair+three-body and purely repulsive models. The pair
cutoff is 1.7 nm. An Axilrod-Teller-Muto three-body term is also studied in the
paper, but the present collaborator benchmark focuses on the pair-potential
`s=1.00` and `s=0.00` models.

## Static NALD reproduction

`static_nald.py` constructs pair forces, the dense equal-mass dynamical
matrix, and centered finite-difference bulk/shear affine-force fields directly
from the supplied coordinates. `build_static_inputs.py` writes the
mass-normalized Hessian and mass-weighted affine vectors as NumPy arrays and can
optionally perform the full eigendecomposition.

The first supplied snapshot of each pair model has also been run through a
force/Hessian/affine preflight; see `preflight_first_snapshot.csv`. This is a
code-path check, not yet a viscosity validation.

For each supplied configuration we will:

1. construct the instantaneous mass-normalized Hessian;
2. compute shear and volumetric affine-force fields on the same snapshot;
3. obtain `Gamma_shear,p` and `Gamma_bulk,p`;
4. establish a physically justified low-frequency/finite-size treatment;
5. compare NALD `eta` and `zeta` with the Green-Kubo targets;
6. test how attraction redistributes `Gamma_bulk,p` relative to
   `Gamma_shear,p`, especially among low-frequency collective modes.

Because the production simulations are NVE, the effective NALD damping/memory
treatment is a central part of the benchmark rather than an externally imposed
Langevin parameter.

The paper also reports a finite-size warning for the full-pair+three-body
model: bulk viscosity changes by about 8% when the system size is increased to
32,000 particles. System-size effects therefore remain an important part of a
rigorous comparison.

## Data provenance

The collaborator-supplied benchmark data are separate from the MIT-licensed
software. See [`data/DATA_NOTICE.md`](data/DATA_NOTICE.md). Journal PDFs and
collaborator source code are not redistributed here.
