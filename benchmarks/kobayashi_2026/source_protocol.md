# Collaborator MD source/manual: reproduction notes

In October 2026 H. Kobayashi supplied the complete in-house MD package and
`MD_Simulation_Manual.docx` used for the benchmark system.  The source and
manual are not redistributed here.  This note records only the conventions
needed to reproduce the pair-force/Hessian calculation in NALD.

## 124 K benchmark state

- particles: 1372 Ar atoms
- temperature: 124 K
- density: 1.1249 g/cm^3 = 1124.9 kg/m^3
- argon mass used by the source: 6.6335e-26 kg
- resulting cubic cell length: 4.3250816376 nm
- periodic minimum-image boundaries

The supplied input template uses a 3 fs timestep.  After run01 the shell script
switches off velocity scaling and subsequent runs are NVE.  The published
production protocol should be followed when reproducing time-correlation
functions; the timestep does not enter the static NALD Hessian itself.

## Pair potential and attraction scaling

The production source calls the TT/JHBV two-body force.  The benchmark varies
only the attractive scaling factor; the repulsive scaling factor remains 1.

The TT minimum used by the code is

    r_min = 3.7618 Angstrom

with

    epsilon_TT = 1.976041929857538e-21 J.

For r < r_min the force is the full TT force.  For r >= r_min the TT force is
multiplied by the attraction scale s.

The source generates C12, C14 and C16 recursively from C6, C8 and C10.  This
gives C16 = 1.170063435858...e-6 K nm^16, consistent with the corrected
published coefficient.

## Production two-body cutoff

This is an important detail that was not encoded in the first public
implementation:

- for s != 0: r_c = 5 * sigma_MD, with sigma_MD = 3.3952 Angstrom,
  hence r_c = 1.6976 nm;
- for s = 0: r_c = r_min = 0.37618 nm.

The force is simply omitted outside the selected cutoff.  The NALD Hessian and
affine-force implementation now mirrors this convention.

## ATM term

The source contains an Axilrod-Teller-Muto term with coefficient

    nu_ATM = 7.32e-108 J m^9

and a three-body cutoff of L/4.  The ATM call is disabled in the supplied pair
template unless explicitly uncommented.  The present NALD benchmark therefore
starts with the pair-only s=1.00 and s=0.00 models.

## Stress/pressure outputs

The manual identifies:

- `prss.dat`: pressure needed for bulk Green-Kubo viscosity;
- `tensol_ab.dat`: off-diagonal stress components for shear viscosity;
- `rxyz.dat`: coordinates in metres.

These conventions are sufficient to build the static NALD Hessian and
volumetric/shear affine-force vectors on exactly the supplied coordinates.
