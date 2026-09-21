# NALD viscosity

Non-affine lattice-dynamics (NALD) tools for computing **shear viscosity** and
developing a corresponding **bulk-viscosity** formulation from atomistic
Hessian eigenmodes and affine-force fields.

For a deformation channel (M) (shear: (M=G); bulk: (M=K)), the Markovian
modal response used here is

```text
M''(Omega) = (1/V) sum_p Gamma_p nu Omega /
             [ (lambda_p - Omega^2)^2 + (nu Omega)^2 ]
```

and the zero-frequency viscosity is

```text
eta_s or zeta = lim_(Omega->0) M''(Omega)/Omega
              = (1/V) sum_p Gamma_p nu / lambda_p^2 .
```

For shear, `Gamma_p=(e_p.Xi_shear)^2`. For bulk, the code uses the volumetric
strain `epsilon_v=Delta V/V` and `Gamma_p=(e_p.Xi_bulk)^2`, with
`Xi_bulk=-df/d epsilon_v`.

## Status

**Shear viscosity:** regression-tested against the historical 300 K atomistic
epoxy calculation and consistent with the NALD viscosity formulation used in
Singh et al., *J. Chem. Phys.* **162**, 244504 (2025).

**Bulk viscosity:** implemented as a deformation-channel extension and prepared
for quantitative validation. Published argon bulk-viscosity values from
Kobayashi, Ishii, and Ohtori, *J. Chem. Phys.* **165**, 104506 (2026), are
included as a benchmark target. The bulk channel should be regarded as
**research/development code until an end-to-end liquid benchmark is completed**.

## Repository layout

```text
nald_viscosity/                 core modal response and unit conversions
scripts/                        command-line tools
lammps/atomistic_real/          shear/bulk affine-force templates
examples/epoxy_validation/      historical shear-viscosity regression
benchmarks/kobayashi_2026/      published simple-liquid reference data
tests/                          numerical unit tests
```

## Quick start

Install locally:

```bash
python -m pip install -e .
```

Given eigenvalues and an affine-force correlator:

```bash
python scripts/viscosity_nald.py \
  --eigenvalues eigenvalues.data \
  --gamma gamma_shear.data \
  --channel shear \
  --volume 96003.3 \
  --nu 5.0e13 \
  --input-units lammps-real-massweighted
```

For bulk, use `gamma_bulk.data` and `--channel bulk`.

To diagonalize one Hessian and project both channels:

```bash
python scripts/diagonalize_channels.py \
  --hessian Hessian.dat \
  --affine shear=AF_shear.data \
  --affine bulk=AF_bulk.data \
  --outdir modes
```

The new workflow mass-weights the affine force as `Xi_i/sqrt(m_i)` before
projection, consistent with a mass-normalized dynamical matrix.

## Low-frequency modes

Negative instantaneous-normal modes are retained. Numerical zero modes are
handled separately from any physical finite-size cutoff.

For the shear channel, Singh et al. (2025) use

```text
omega_min = (2*pi/L) * sqrt(G_s/rho)
```

as the minimum propagating shear-mode frequency supported by the finite box.
Use `scripts/finite_size_cutoff.py` to evaluate it. No universal
bulk/longitudinal analogue is assumed in this repository; that prescription
must be validated independently.

## Validation

The epoxy regression reproduces the historical `G_dp.f90` loss-modulus curve
to about (10^{-8}) relative accuracy when the historical damping, cutoff,
unit conversion, and final-mode omission are reproduced exactly.

The legacy epoxy calculations are kept only as regression references. Their
hard-coded low-frequency cutoffs are **not** used automatically in the general
workflow.

The Kobayashi benchmark contains published values of bulk viscosity, shear
viscosity, (K_\infty-K_0), and pressure-relaxation time at three liquid-argon
state points. The complete Tang-Toennies/JHBV potential parameter set and
instantaneous configurations are not contained in the supplied publication
files, so the benchmark is currently a target dataset rather than a complete
reproduction deck.

## References

1. A. Zaccone, “General theory of the viscosity of liquids and solids from
   nonaffine particle motions,” *Phys. Rev. E* **108**, 044101 (2023).
   DOI: 10.1103/PhysRevE.108.044101
2. A. Singh, V. Vaibhav, T. W. Sirk, and A. Zaccone, “Viscosity of polymer
   melts using non-affine theory based on vibrational modes,”
   *J. Chem. Phys.* **162**, 244504 (2025). DOI: 10.1063/5.0272171
3. V. Vaibhav, T. W. Sirk, and A. Zaccone, “Time-Scale Bridging in Atomistic
   Simulations of Epoxy Polymer Mechanics Using Nonaffine Deformation Theory,”
   *Macromolecules* **57**, 10885–10893 (2024).
   DOI: 10.1021/acs.macromol.4c01360
4. H. Kobayashi, Y. Ishii, and N. Ohtori, “Structural origin of the strong
   effect of attraction on bulk viscosity in simple liquids,”
   *J. Chem. Phys.* **165**, 104506 (2026). DOI: 10.1063/5.0351887

## License

MIT. See [LICENSE](LICENSE).
