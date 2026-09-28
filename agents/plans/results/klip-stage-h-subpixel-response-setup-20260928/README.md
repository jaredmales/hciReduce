# KLIP Stage-H subpixel-response setup

## Purpose

Stage G showed that the integer-grid exact response is accurate for an
integer-centered companion but loses fidelity at the planet's subpixel phase.
Stage H tests whether this is primarily a registration error that can be
repaired by shifting the existing response, or a phase-dependent change in the
KLIP response shape that requires newly measured templates.

The experiment reuses the twelve completed nominal optimized-phase Stage-G
positive injections. It adds one matching negative reduction at each site and
forms the paired central response

\[
t_{\mathrm{pair}} =
\frac{I(+c)-I(-c)}{2c}.
\]

This is the same finite-amplitude response definition used by the promoted
integer-grid `refitDifference` library, now evaluated at the planet-like
fractional phase.

## Frozen comparison

All twelve sites use the Stage-G optimized phase
`(-0.27707,+0.48596)` pixels and contrast
`0.004574362496414845`. Mode 200 is primary; all eight retained KL modes are
reported. Three templates are compared on the same complete 11-by-11 support:

| Template | Construction | Role |
| :--- | :--- | :--- |
| Current integer | Promoted 47-pixel exact response at the nearest native pixel, cropped to 11 pixels | Existing method |
| Cubic shifted | The same 47-pixel response shifted by the known row/column phase with third-order interpolation, then cropped | Cheap subpixel model |
| Paired exact | Full positive-minus-negative KLIP response at the actual fractional source location | Regenerated subpixel reference |

The shift is applied before cropping so that the 11-pixel result does not lose
edge information. Every task must have a complete 15-by-15 validity guard
around the retained stamp.

The paired response is also compared with the one-sided positive and negative
responses,

\[
t_+ = \frac{I(+c)-I(0)}{c},
\qquad
t_- = \frac{I(0)-I(-c)}{c},
\]

to separate subpixel-template error from finite-amplitude asymmetry.

## Endpoints

For identity, raw rectangular covariance, and radial-Hann covariance, the
runner reports:

- template cosine;
- matched-filter projection scale; and
- best-scaled relative residual.

The covariance models and radial normalization come from the unchanged
signal-free Stage-F baseline and use the same source-safe training geometry.
Cosine is the fraction of optimal matched-filter SNR retained by the template
under the stated covariance metric. Projection scale is the noiseless contrast
gain. These response endpoints isolate template fidelity; this experiment
does not redefine or re-estimate annular detection SNR.

The predeclared mode-200 comparisons are shifted-minus-integer cosine,
integer-minus-shifted residual, and the reduction in absolute projection-scale
error. If the shifted model approaches the paired-exact ceiling, a cubic phase
shift is sufficient for the next planet analysis. If a substantial residual
remains, the next implementation should measure a grid of fractional response
phases through the full KLIP reduction.

## Workload and provenance

Only twelve new KLIP reductions are required. The corresponding positives,
sites, phase, contrast, calibration products, exact-response field, and
software boundary come from completed Stage G. The runner fingerprints every
reused positive and all response products, freezes the twelve negative
commands, and is resumable at reduction and analysis granularity.

Each completed site retains the central, positive, negative, integer, and
shifted mode-200 response stamps for review. The aggregate JSON covers all
modes; CSV and Markdown reports present the primary mode.

## ROC commands

After pulling the setup commit, run from the repository root:

```bash
root=working/roc/klip_stage_c_development_20260925
python3 agents/plans/scripts/run_klip_stage_h_subpixel_response.py check
taskset -c 12-27 python3 \
  agents/plans/scripts/run_klip_stage_h_subpixel_response.py prepare "$root"
```

Then launch the frozen command printed by `prepare`:

```bash
root=working/roc/klip_stage_c_development_20260925
taskset -c 12-27 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 python3 \
  "$root/stage_h_subpixel_response/software/run_klip_stage_h_subpixel_response.py" \
  run "$root" \
  > "$root/stage_h_subpixel_response/driver.log" 2>&1
```

The command can be left to completion and safely rerun. Progress is available
with:

```bash
tail -f working/roc/klip_stage_c_development_20260925/stage_h_subpixel_response/driver.log
```

The primary report will be
`working/roc/klip_stage_c_development_20260925/stage_h_subpixel_response/results.md`.
