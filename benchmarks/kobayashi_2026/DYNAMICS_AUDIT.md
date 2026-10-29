# Kobayashi benchmark: remaining static-vs-dynamic audit

The public Green-Kubo integration check and the first-snapshot
force/Hessian/affine preflight are already validated and are not repeated here.

## Evidence table

| Question | Current repository evidence | Smallest remaining test |
| --- | --- | --- |
| Is the source-matched static calculation complete? | `analyze_mode_ensemble.py`, `physical_cutoff_summary.csv` | Treat the existing 10-snapshot eigenvalues and affine-force couplings as fixed inputs. |
| Is the acoustic finite-size scale independent of viscosity? | `acoustic_finite_size.py`, `acoustic_qmin_summary.csv` | Yes. Recompute summaries with `analyze_physical_cutoff_ensemble.py` after any acoustic-projection change. |
| Does a hard q_min cutoff remove only acoustic finite-size artifacts? | `low_mode_character.py`, `low_mode_audit_sample.csv` | No clear evidence. The 0/5/9-frame audit shows similar q_min plane-wave weight for s=1 and s=0 subcutoff bulk weight, but very different participation ratios. Run all ten frames before adopting a mode-character filter. |
| Is the missing factor demonstrably a memory effect? | `pressure_memory.py`, `memory_nald_bridge.csv` | Not yet. The current pressure-memory diagnostic uses the same pressure correlation that defines the Green-Kubo target. |
| What independent dynamical input is needed? | `vacf_memory.py` | Supply a VACF or force autocorrelation from the same NVE state points and reconstruct the memory without pressure/stress data. |
| How should the full comparison be made? | `green_kubo_spectrum.py`, generalized memory support in `nald_viscosity/core.py` | Compare frequency-resolved NALD K''/Omega and G''/Omega against bulk- and shear-stress Green-Kubo spectra, not only zero-frequency ratios. |

## Verified subcutoff sample

Frames 0, 5, and 9 were audited for both interaction models.

For the **bulk** modes below the independently estimated longitudinal
q_min cutoff:

- s=1.00: A-weighted participation ratio = 0.274 +/- 0.049;
  A-weighted q_min plane-wave weight = 0.0352 +/- 0.0051.
- s=0.00: A-weighted participation ratio = 0.0573 +/- 0.0073;
  A-weighted q_min plane-wave weight = 0.0348 +/- 0.0033.

Thus the strong difference in subcutoff structural weight is not accompanied
by a corresponding increase in q_min plane-wave character.  The s=1 modes
carrying that weight are substantially more extended than the s=0 modes.

This does **not** prove that the subcutoff s=1 modes should be retained in a
transport calculation: instantaneous soft modes can be transient and their
lifetime belongs to the dynamical/memory problem.  It does show that a blanket
hard frequency cutoff cannot by itself distinguish finite-size acoustic modes
from other soft instantaneous modes.

## Minimal additional collaborator input

Preferred, in order:

1. **VACF or force autocorrelation** for pair-only s=1.00 and s=0.00 at
   124 K and 1124.9 kg/m^3, from the same NVE production protocol.  A
   precomputed correlation is sufficient; raw velocities/forces are not needed
   if the normalization, units, timestep, and sampling interval are supplied.
2. **Off-diagonal stress autocorrelation / shear Green-Kubo integrand** for the
   same two models, ideally averaged over xy, xz, yz.  This provides the
   frequency-resolved shear control channel.
3. If already available, **F(q,t) or S(q,omega)** at q approximately
   2*pi/L = 1.453 nm^-1 (and optionally the next q shell).  This would test the
   proposed low-q density/pressure mechanism independently of the pressure
   autocorrelation.

If precomputed correlations are unavailable, a short contiguous NVE trajectory
containing velocities, off-diagonal stress, and positions at sufficiently fine
sampling can generate all three diagnostics locally.
