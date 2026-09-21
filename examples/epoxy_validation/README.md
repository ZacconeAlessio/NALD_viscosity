# Epoxy shear-viscosity regression

This directory preserves two historical calculations as regression references
for the general NALD viscosity implementation.

The reference T=300 K atomistic epoxy dataset contains 30,222 modes
(`N=10,074`) and is intentionally not committed here. With the historical
`G_dp.f90` choices reproduced exactly, the loss-modulus curve is recovered to
about (10^{-8}) relative error and the zero-frequency shear viscosity is

```text
eta_s = 0.033905666086 Pa s
```

for the historical `|lambda|<2.57` cutoff.

A separately supplied direct-viscosity script used `|lambda|<1.0`, giving

```text
eta_s = 0.147762800487 Pa s
```

on the same run. With no finite cutoff beyond the numerical zero-mode
tolerance, the result rises to roughly (7.44e3) Pa s. This illustrates the
extreme sensitivity of the zero-frequency sum to the lowest modes.

These numerical windows are **legacy regression choices, not universal physical
cutoffs**. For the shear channel, use the finite-size criterion documented in
Singh et al., J. Chem. Phys. 162, 244504 (2025) when the required (L),
(ho), and (G_s) are available.

To validate against a local reference-data directory:

```bash
python examples/epoxy_validation/validate_epoxy_reference.py /path/to/run1
python examples/epoxy_validation/legacy_epoxy_zero_frequency.py /path/to/run1
```
