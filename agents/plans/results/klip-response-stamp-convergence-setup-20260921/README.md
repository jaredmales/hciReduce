# KLIP response-stamp convergence setup

## Question

Stage A found that the archived 11-by-11 exact response leaves 1.35--3.63%
median mode-200 squared energy on its border at every primary radius. This
experiment determines the smallest larger footprint that contains those lobes
without spending one full exact-response campaign on each candidate size.

## Fixed design

The experiment uses the completed canonical Stage-A run as its immutable
parent. At each of radii 7.5, 10, 12, 16, 20, and 24 pixels, it selects 12
sites from the complete five-pixel geometry. Selection uses only coordinates:
candidates are sorted by angle and the midpoint of each of 12 equal-count
blocks is retained. No baseline, response-energy, or source-recovery value is
used in site selection.

At each of the 72 sites, one reduction receives a positive perturbation and
one receives an equal negative perturbation. Both reductions first subtract
the fitted known planet. The perturbation half-contrast is
0.0045743624964148452, matching the archived paired response calculation.
Their full-image central difference supplies candidate response stamps of 11,
15, 19, 23, and 31 pixels. Diagnostic extractions at 39, 47, 55, and 63 pixels
measure longer-range response energy and test whether a candidate-size border
minimum is only a local zero crossing. All nine sizes share the same 144 KLIP
reductions.

All eight KL modes are analyzed. Mode 200 controls the footprint decision; the
other modes are mandatory diagnostics. For every site, mode, and size, report:

- finite support fraction;
- total squared response energy;
- squared energy on the outermost border;
- negative-lobe energy fraction; and
- energy relative to the 63-by-63 diagnostic extraction.

The external 11-by-11 finite difference must reproduce the archived native
exact response with cosine at least 0.999 and projection within 0.01 of unity
at every site and mode. A larger size passes the original footprint rule only
if every selected mode-200 stamp has complete finite support, every radial bin
has median border energy no larger than 1%, and no individual site has border
energy above 5%. The smallest passing tested size is the target for one new
full-field response campaign. If none passes, the result identifies whether
support at the inner boundary or still-extended response energy is responsible.
Sizes above 31 pixels diagnose nonlocal tails; they are not automatically
promoted because their covariance dimension and radial support depart from the
planned local-filter experiment.

At the measured 5.06 seconds per canonical Stage-A baseline reduction, the
144 reductions should take about 12--15 minutes on the fixed 16-CPU ROC
allocation.

## ROC commands

```bash
root=working/roc/klip_response_stamp_convergence_20260921

taskset -c 12-27 python3 \
  agents/plans/scripts/run_klip_response_stamp_convergence.py check

taskset -c 12-27 python3 \
  agents/plans/scripts/run_klip_response_stamp_convergence.py prepare "$root"

taskset -c 12-27 python3 \
  "$root/software/run_klip_response_stamp_convergence.py" run "$root"
```

The run is resumable. It verifies completed FITS products before reuse and
refuses to overwrite a partial task directory. When all 144 reductions are
complete, it automatically writes `results.json`, `results.md`, and a
fingerprinted completion receipt.
