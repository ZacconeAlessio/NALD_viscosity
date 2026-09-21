# Kobayashi et al. JCP 2026 benchmark

Reference: H. Kobayashi, Y. Ishii, and N. Ohtori, “Structural origin of the
strong effect of attraction on bulk viscosity in simple liquids,”
*J. Chem. Phys.* **165**, 104506 (2026), DOI 10.1063/5.0351887.

`reference.csv` transcribes the reported bulk viscosity `zeta`, shear
viscosity `eta`, (K_\infty-K_0), and pressure-relaxation time
`tau_zeta` from Table I and Table S1. Viscosity values are tabulated in units
of (10^{-1}) mPa s; multiply by (10^{-4}) for Pa s.

The paper uses

```text
zeta = (K_inf - K_0) * tau_zeta .
```

At 124 K and 1124.9 kg/m^3, reducing the pair-attraction scale from (s=1) to
(s=0) changes the reported quantities approximately as follows:

- `zeta`: 3.79 -> 0.344 (about -91%)
- `eta`: 0.958 -> 0.790 (about -18%)
- `K_inf-K_0`: 0.673 -> 0.3908 GPa (about -42%)
- `tau_zeta`: 0.56 -> 0.088 ps (about -84%)

The especially strong change in the relaxation time makes this a useful target
for the NALD bulk channel.

## Reported simulation protocol

The main production system contains 1372 Ar atoms at three state points:
(140 K, 968 kg/m^3), (124 K, 1124.9 kg/m^3), and
(90 K, 1390 kg/m^3). The pair model is the Tang-Toennies/JHBV argon potential,
with a WCA-like decomposition used to vary attraction strength. An optional
Axilrod-Teller-Muto three-body term is also studied.

The paper reports (2.5e7) NVE production steps, normally with a 6.45 fs
timestep; 3 fs is used for the full-pair+three-body and purely repulsive models.
The pair cutoff is 1.7 nm and the three-body cutoff is one quarter of the box
length.

## Reproducibility limit

The complete Tang-Toennies/JHBV numerical parameter set, instantaneous
configurations, and in-house MD source are not contained in the supplied
article or Supplementary Material. Therefore this directory is a target
benchmark, not a complete reproduction deck.

The paper also notes a finite-size effect for the full-pair+three-body model:
the bulk viscosity changes by about 8% when the particle number is increased to
32,000. Any NALD comparison should therefore examine system-size dependence
rather than treating the 1372-particle value as a thermodynamic-limit result.

## NALD comparison target

For each interaction model/configuration:

1. compute one instantaneous mass-normalized Hessian;
2. compute shear and volumetric affine-force fields on the same snapshot;
3. obtain `Gamma_shear,p` and `Gamma_bulk,p`;
4. establish a physically justified low-frequency finite-size treatment;
5. compare NALD `eta` and `zeta` against `reference.csv`;
6. examine whether attraction preferentially enhances low-frequency
   `Gamma_bulk,p` relative to `Gamma_shear,p`.

The final point provides a modal-space test of the paper's finding that bulk
viscosity is coupled strongly to slow, low-q density fluctuations whereas the
shear channel is much less sensitive.
