# P4 Negative-Companion Optimization Plan

## Goal

Fit and remove one astrophysical point source by repeatedly evaluating the exact finite-amplitude P4 pixel-local
path against one pristine, already-preprocessed observation cube. The optimizer must first solve the bounded
one-dimensional contrast problem at fixed separation and position angle. Joint position and contrast fitting follows
only after the contrast-only result is deterministic and agrees with a dense trial scan.

The optimizer is a scientific layer above `P4Reduction` rather than a new subtraction algorithm. Every merit
evaluation must use the ordinary P4 eigensolve, the complete temporal sample set, the configured final combination,
and the existing local validity rules. Intermediate evaluations remain in memory and do not write FITS products.

## Established prerequisites

The pixel-local implementation provides the required numerical oracle:

- `p4.localStampSize>0` refits every detector regression needed by the local sky stamp after applying the finite
  trial to the target and predictor samples;
- repeated `reduce()` calls retain the loaded target cube because `m_filesRead` remains true;
- the trial is sampled on demand and never added to `m_tgtIms`, so trial flux cannot accumulate;
- each evaluation exposes the local residual cube, local validity cube, continuous source coordinate, full-image
  origin, mode ordering, timing, fit count, and sparse-sample count; and
- output writing can be disabled while the in-memory products remain available.

An A-to-B-to-A contrast regression test now verifies that the first and third local results agree within `1e-6`, the
intermediate contrast changes the result, and every loaded input pixel remains unchanged. A separate file-backed
test verifies that two reductions on one instance invoke the input post-read hook only once.

## Initial performance measurement

The authoritative real-data configuration is the ignored working file
`working/roc/p4Reduce_afLepNaco_pl.conf`. It reads 621 already-preprocessed 256-by-256 images from the coadd5 `pp`
subdirectory, uses 19 contiguous annuli spanning 0 to 40 pixels, 15 mode fractions, skipped preprocessing, an
11-by-11 local stamp, and a trial at separation 12.7 pixels, PA 260.1 degrees, and contrast -0.1.

The turnaround measurement retained the first 96 images from that exact input and all 19 annuli, reduced the output
to three mode fractions, and used four OpenMP workers plus one BLAS worker. `/usr/bin/time -v` measured independent
complete and local processes. The complete comparison did not materialize a trial source, so it slightly favors the
complete run.

| Measurement | Complete field | Finite-amplitude local | Local / complete |
|---|---:|---:|---:|
| Total elapsed | 6.030045 s | 2.278965 s | 0.378 |
| P4 algorithm elapsed | 4.914753 s | 1.639154 s | 0.334 |
| Geometry elapsed | 0.598039 s | 0.614883 s | 1.028 |
| Aggregate measured worker time | 19.084028 s | 5.527049 s | 0.290 |
| Unique detector fits | 5,024 | 256 | 0.051 |
| Retained/requested residual samples | 482,304 | 20,803 | 0.043 |
| Peak RSS | 987,556 KiB | 940,372 KiB | 0.952 |

The local run is about 2.65 times faster end to end and 3.00 times faster in the P4 algorithm even at 96 frames.
Sparse trial-source sampling still accounts for 4.993576 worker seconds, or 90.35% of measured local worker time,
and geometry takes the same approximately 0.6 seconds in both paths because annular geometry is currently rebuilt in
full. The loaded input cube dominates peak RSS, so sparse output processing lowers process peak memory by only about
4.8% in this configuration.

The full 621-frame `finim_0038.fits` trial from the same configuration used 487 unique detector fits and 131,180
sparse residual samples. Its local geometry retained 31,806,752 bytes, all 121 output pixels were valid in all modes,
and the negative response aligned with the continuous source coordinate. The user observed that run to be very fast;
its console wall-time and peak-RSS values were not captured and remain an integration measurement.

The first real-data contrast-only optimizer run completed 26 evaluations in 109.125764 seconds of accumulated local
evaluation time, or 4.197 seconds per evaluation. Its 21-point dense scan placed the minimum at -0.005 with merit
0.104920010; Brent refinement converged at contrast -0.004117980077 with merit 0.09685237217. The unperturbed
contrast-zero merit was 0.2726269414, so the fitted negative source reduced the selected-mode aperture merit by
64.5%. The summary, merit CSV, 11-by-11-by-15 best residual/validity cubes, FITS provenance, and fitted configuration
fragment agreed with the console result.

The first real-data joint run used a 192-evaluation budget and finished 181 distinct evaluations: 26 seed-contrast,
129 simplex, and 26 final-contrast calls. It reported separation 12.4508435 pixels, PA 261.3958089 degrees, contrast
-0.004552707737, and merit 0.07414968197, but did not satisfy position convergence. The merit history diagnosed a
hard-aperture discontinuity rather than a simple lack of budget: nearby trials alternated between 79-pixel merit
near 0.07415 and 80-pixel merit near 0.07967. The selected point placed full-image pixel `(135,127)` only
0.000004 pixel outside the radius-5 aperture. The moving aperture had therefore allowed the optimizer to improve its
score by excluding a boundary pixel. Joint fitting now fixes the aperture on the initial continuous sky coordinate,
records that center in its products, and uses a separate practical Cartesian tolerance.

## Optimizer boundary

### Milestone 1: contrast only

Hold separation and PA fixed at the configured single `fake.sep` and `fake.PA`. Optimize the signed injected
`fake.contrast` over an explicit bounded interval whose upper bound is nonpositive. Evaluate one selected P4 output
mode. This isolates the photometric sign, merit function, state reuse, and convergence behavior from astrometric
degeneracy.

Use a deterministic bounded scalar method. Brent minimization is preferred after evaluating the two bounds and one
interior point; fall back to golden-section steps if a parabolic proposal is unsafe. Do not assume the P4 merit is
differentiable because rank thresholds and sigma-clipped combination can introduce discrete changes.

Before accepting the optimizer result, evaluate a dense contrast grid spanning the same bounds and require the
reported minimum to agree with the sampled basin. The grid is a validation product, not the long-term optimizer.

### Milestone 2: position and contrast

The implemented full optimizer fits two local sky offsets and signed negative contrast. Cartesian sky offsets are the
internal parameters, so PA wrap and the singularity at zero separation do not affect the search. The exact conversion
uses `row=separation*sin(-PA)` and `column=separation*cos(-PA)` and converts the fitted coordinate back to separation
and PA for configuration and FITS provenance, preserving the `0.5*(N-1)` center convention used by local processing.

A bounded three-dimensional Nelder-Mead simplex is seeded by the contrast-only result. Every trial reflects back into
the inclusive square Cartesian bounds and signed contrast interval. Position changes rebuild the sparse output
dependencies while repeated reductions reuse the loaded pristine target cube. The optimizer reserves a second
dense+Brent contrast fit at the selected position, so the published photometry and dense-basin check apply at the
fitted astrometry rather than only at the initial position. A combined absolute/relative merit tolerance uses
`meritTolerance*max(1,abs(merit))`, avoiding impossible relative convergence when a synthetic or signal-free merit is
near zero. `positionTolerance` independently controls the Cartesian simplex diameter; `parameterTolerance` remains
the contrast-width tolerance.

## Merit function

The initial merit for selected output plane `m` is the mean squared valid residual inside a circular aperture:

```text
J(theta) = sum_i valid_i * w_i * R_m(theta, i)^2 / sum_i valid_i * w_i.
```

Here `theta` is the current trial and `R_m` is the combined local residual. The aperture is centered on the initial
configured continuous sky coordinate rather than the nearest integer stamp pixel or the moving trial coordinate. A
pixel is included by its full-image center distance from that fixed coordinate, so every trial uses identical sky
pixels. Reject a trial if the fixed aperture extends outside its local stamp or if any required validity is absent;
do not reinterpret invalid samples as zero.

Start with uniform weights and an aperture radius of `2*p4.psfRadius`. For the AF Lep configuration this is five
pixels. An 11-pixel stamp supports contrast-only optimization; a one-pixel joint position bound requires at least a
13-pixel stamp so the same aperture remains present throughout the search. Noise weights or an L1/Huber merit can be
added only as explicit alternatives with injection-recovery tests; they must not silently change the baseline
definition.

Use exactly one configured P4 mode fraction initially. Highly correlated mode planes are not independent
measurements, so summing them without a covariance model would give a misleading merit and uncertainty. Report the
merit curve for every available mode as a diagnostic, but optimize the user-selected mode only. The selected
fraction must match one configured `p4.modeFractions` entry; never silently choose the nearest plane.

## Proposed configuration

Keep the existing `[fake]` values as the initial trial and add a separate optimizer group:

```ini
[p4Optimize]
enabled=true
modeFraction=0.15
apertureRadius=5
fitPosition=false
contrastLower=-0.05
contrastUpper=0
maxEvaluations=64
parameterTolerance=1e-5
positionTolerance=0.001
meritTolerance=1e-6
validationSamples=21
positionBound=1
uncertaintyBlocks=8
outputPrefix=p4Negative_
```

For the joint milestone set `fitPosition=true`, retain `positionBound=1`, set `p4.localStampSize=13`, and raise
`maxEvaluations` to 192. With 21 validation samples the checked mathematical minimum is 78, but that leaves only four
simplex vertices and is not a useful convergence budget.

Proposed meanings:

| Target | Contract |
|---|---|
| `enabled` | False preserves the ordinary one-shot `p4Reduce` path; true runs the selected optimizer mode. |
| `modeFraction` | Exact member of `p4.modeFractions` used for optimization. |
| `apertureRadius` | Positive circular merit radius in pixels; default `2*p4.psfRadius`. |
| `fitPosition` | False for contrast-only milestone; true enables two Cartesian sky offsets plus contrast. |
| `contrastLower`, `contrastUpper` | Finite ordered signed bounds with `lower < upper <= 0`. |
| `maxEvaluations` | Hard upper bound on local P4 calls, excluding optional uncertainty resamples. Joint mode requires at least `2*(validationSamples+16)+4`; 192 is recommended for 21 validation samples. |
| `parameterTolerance` | Absolute contrast-width tolerance for both contrast stages and the simplex contrast coordinate. |
| `positionTolerance` | Absolute Cartesian simplex-diameter tolerance in pixels for milestone 2; default `0.001`. |
| `meritTolerance` | Combined absolute/relative numerical allowance used for simplex convergence and dense-basin agreement. |
| `validationSamples` | Uniform bounded contrast samples, including both endpoints, used to validate the fitted basin. |
| `positionBound` | Positive symmetric per-coordinate offset bound in pixels around the initial configured source when `fitPosition=true`. |
| `uncertaintyBlocks` | Number of contiguous time blocks used by the optional delete-one-block jackknife. |
| `outputPrefix` | Prefix for the summary, sampled merit curve, and best local products. |

The configuration loader must reject an initial contrast outside the bounds, a missing selected mode, an aperture
that cannot remain inside the stamp over the position bounds, fewer than the required valid pixels, and every local
configuration already rejected by `P4Reduction`.

## Reusable evaluation API

Add a small typed trial/result seam instead of having optimizer code mutate several vectors directly:

```text
P4LocalTrial { separation, positionAngle, contrast }
P4LocalEvaluation { residual, validity, origin, source, timing, counts }
P4Reduction::evaluateLocal(trial)
```

The first implementation may internally use the existing reduction lifecycle, but it must:

1. require a loaded or loadable pristine observation;
2. disable intermediate FITS writing;
3. replace, rather than append to, the single trial tuple;
4. clear only per-evaluation results and timing;
5. preserve the loaded cube, headers, angles, input scaling, and file selection; and
6. return an owning or lifetime-safe result that cannot be overwritten unexpectedly by the next call.

Measure the fraction of each evaluation spent rebuilding annular geometry, sparse dependencies, the prepared source,
and PCA fits. Cache only immutable state with explicit invalidation when geometry, position, template, input frames,
or P4 settings change. Contrast-only optimization should ultimately reuse both annular and sparse geometry because
only the template amplitude changes.

## Output and provenance

Write no FITS product per optimizer evaluation. At completion write:

- a machine-readable summary containing the initial and fitted trial, selected mode, aperture, bounds, convergence
  status, evaluation count, final merit, timing totals, and software provenance;
- a sampled merit table sufficient to reproduce the contrast fit or joint seed/simplex/final-validation sequence;
- the best local residual and validity FITS cubes using the existing local headers plus optimizer provenance; and
- a configuration fragment containing the fitted negative source for the final complete reduction.

An optimizer failure or exhausted evaluation budget must retain the best finite trial but mark it unconverged. It
must not publish that trial as a completed fit.

## Uncertainty method

Do not interpret the optimizer curvature as a statistically calibrated uncertainty: local pixels and P4 mode planes
are correlated, and sigma-clipped combination may not be smooth. Use a contiguous delete-one-time-block jackknife as
the initial uncertainty method. Refit after omitting each of at least eight time blocks, then calculate the ordinary
jackknife covariance of contrast and, when enabled, Cartesian position.

Treat the curvature and dense-grid width as convergence diagnostics only. Validate the jackknife uncertainty later
with negative-source injection/recovery trials at nearby empty locations and comparable separation.

## Final signal-free validation

After convergence, run one complete P4 reduction with the fitted negative source applied at the same
post-preprocessing P4-input stage. This requires an explicit complete-field post-preprocessing injection path; the
current inherited fake hook is intentionally inactive under `preProcess.skip=true`. Implement that path without
changing the local optimizer contract and without accumulating flux in the optimizer's pristine cube.

From that complete signal-free reduction:

1. regenerate the spatially variable P4 PSF field;
2. regenerate the normalized filtered product;
3. compare the response-energy denominator with the original trough;
4. compare the filtered radial structure/ring with the original product; and
5. verify that ordinary unfiltered SNR away from the fitted aperture is unchanged within the expected numerical and
   statistical variation.

## Work sequence

1. [x] Add the A-to-B-to-A pristine-state and one-read regression tests.
2. [x] Record local/complete benchmark measurements and the large-data integration geometry.
3. [x] Add `P4LocalTrial`, `P4LocalEvaluation`, and `evaluateLocal()` without changing the CLI result.
4. [x] Implement and test the aperture merit independently of the optimizer.
5. [x] Implement bounded contrast minimization and dense-grid validation.
6. [x] Add summary/table/best-local outputs with complete provenance.
7. Benchmark repeated in-process evaluations on the remote 1-second-coadd dataset.
8. [x] Add joint Cartesian position/contrast fitting.
9. Add delete-one-block jackknife uncertainty.
10. Add the one-shot complete-field post-preprocessing injection and perform the signal-free PSF/filter comparison.

## Acceptance criteria

- An A-to-B-to-A trial sequence reproduces A and leaves the loaded cube bitwise unchanged.
- File discovery, loading, and preprocessing happen once per optimizer session.
- Contrast-only optimization agrees with a dense bounded scan on synthetic and AF Lep/NACO cases.
- Invalid or partially supported aperture samples cause a checked failure rather than a biased merit.
- Intermediate evaluations do not write files and cannot accumulate trial flux.
- The selected mode, merit aperture, bounds, convergence state, evaluation count, and best trial are explicit in
  persisted provenance.
- Joint astrometry uses continuous coordinates and preserves the `0.5*(N-1)` center convention.
- Reported uncertainty comes from the documented resampling method, not uncalibrated optimizer curvature.
- One complete fitted-negative reduction provides the final science and PSF/filter validation products.

## Confirmed implementation decisions

1. Confirm one selected mode rather than a joint-mode merit for the initial optimizer.
2. Confirm uniform L2 inside `2*p4.psfRadius` as the baseline merit.
3. Confirm signed negative contrast bounds, with the current `[fake]` contrast as the initial point.
4. Confirm contrast-only first, followed by bounded Cartesian position plus contrast.
5. Confirm delete-one-time-block jackknife as the initial reported uncertainty.
