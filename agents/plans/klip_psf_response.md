# KLIP Sparse PSF Response and Matched-Filter Plan

## Objective

Produce a useful KLIP matched-filter response directly from the original science reduction. The normal path must not
require a negative-companion fit followed by a second zero-signal fake-injection reduction merely to obtain the
filter. Full fake injection and refitting remain validation oracles for throughput and nonlinear basis changes.

The preferred KLIP response is the exact derivative for a frozen, target-specific KL basis. Unlike P4's
`detectorLocal` approximation, KLIP does not have a separately fitted operator at every search pixel: a sparse sky
probe can be propagated directly through each already-computed regional basis. A detector-local approximation should
therefore be considered only if measurement cost remains unacceptable after bounded accumulation and profiling.

## Existing implementation

The first sparse radial estimator landed in hciReduce commit `3e0db39`.

- `src/common/KLIPPSFModel.hpp` and `.cpp` construct a centered detector-frame probe, apply supported linear regional
  centering, calculate `p - Z Z^T p` for the live target basis, and accumulate each regional contribution directly
  into a compact derotated response stamp.
- `src/common/KLIPreduction.hpp` maps each requested radius/angle to the nearest eligible exact final-image detector
  pixel, records that pixel's actual radius/angle, invokes the probe while each target-specific basis is resident,
  combines target-frame stamps, fits one `RadialPSFModel` per requested mode count, and optionally writes canonical
  radial response and validity cubes. No post-measurement registration is applied.
- `src/common/RadialPSFModel.hpp` and `.cpp` resolve fixed-count or maximum-arc angular sampling, rotate measurements
  to the positive-column radial orientation, average samples at the same radius, and evaluate arbitrary positions by
  linearly interpolating bracketing radial responses.
- `tests/common/KLIPPSFModel_test.cpp`, `tests/common/KLIPreduction_test.cpp`, and
  `tests/common/RadialPSFModel_test.cpp` cover probe preprocessing, frozen-basis projection, compact derotation,
  radial alignment, linear radial evaluation, orchestration, and persisted schema-1 products.
- `doc/klip.dox`, `doc/klip_algorithm.dox`, and `doc/klip_config.dox` describe the public controls and initial
  restrictions.

P4 and KLIP now share the same configuration surface. The former response keys in `[p4]` and `[klip]` are not
registered or loaded; response controls live only in `[psfResponse]`:

```ini
[psfResponse]
file=/path/to/centered_psf.fits
stampSize=11
sampleRadii=1,3,5,7
samplesPerRadius=16
method=refitDifference
sampleAvoidRadius=5
refitContrast=0.0048
outputModels=true
filter=true
filterMinGoodFract=1
outputPrefix=klipPSF_
```

`src/common/PSFResponseConfig.hpp` owns this vocabulary, the method enumeration, and explicit-versus-region radial
grid resolution. `P4Reduction` and `KLIPreduction` inherit it and apply algorithm-specific validation. Both support
`skyExact` and `refitDifference`; P4 additionally supports `detectorLocal`.

Maintained ROC experiment assets are under `agents/plans/scripts`:

- `klipReduce_afLepNaco_psf_response.conf` is the response-compatible AF Lep/NACO KLIP base configuration. It follows
  `working/kr.conf`, uses the current ROC angle constant, and deliberately selects supported `imageMean` centering.
- `run_klip_psf_response_experiment.sh` runs a science-only control, a 1-pixel/16-angle radial reference, and several
  coarser radial/angular grids while recording exact commands, hashes, wall time, peak RSS, and completion state.
- `compare_klip_psf_response.py` linearly evaluates candidate response/validity cubes at the fine-reference radii,
  verifies the response
  run leaves the final science cube unchanged, and reports template-level matched-filter amplitude/cosine proxies.
  Those proxies do not replace the later negative-fit/zero-signal-injection inference oracle.
- `fit_klip_matched_response.py` fits a bounded local likelihood peak directly from the production filtered,
  normalization, support, and validity cubes for an exact configured KL mode count.
- `run_klip_matched_response_validation.sh` runs the accepted fixed-16 lambda/D-scale response grid with normalized
  filtering enabled, verifies the science cube against a science-only control, and fits every configured AF Lep KL
  mode by default. The first response-backed run is recorded below; the negative-companion oracle remains
  outstanding.
- `run_klip_finite_response_validation.sh` reuses that completed response run and performs only one inexpensive
  end-to-end KLIP reduction after subtracting the accepted P4 exact-negative companion. It is restartable under an
  exact settings check and records executable, input, template, configuration, and reference-product provenance.
- `compare_klip_finite_response.py` forms the original-minus-signal-cancelled response per unit contrast, reconstructs
  the production cubic/radial frozen response at the planet pixel, verifies that reconstruction against the persisted
  filtered cube, and measures response cosine, projection scale, and best-scaled residual for every KL mode. It also
  writes the empirical, frozen, and difference stamp cubes for inspection.
- `run_klip_central_response_validation.sh` measures paired positive/negative end-to-end KLIP differences at 15
  bounded positions and three perturbation amplitudes. The default clear samples cover the inner and outer region
  boundaries, the planet's radial neighborhood, and a mid-radius control; all remain more than 13 pixels from the
  known candidate. A separately labelled candidate sample connects the local derivative to the completed one-sided
  finite-amplitude test. The 90 default reductions should take approximately 4.5 minutes at the measured ROC rate.
- `compare_klip_central_response.py` compares every central difference with the production frozen response and, at
  the candidate, with the one-sided finite secant. It reports amplitude dependence, positive/negative secant
  asymmetry, response scale and shape, and the scatter remaining after clear samples at a common radius are rotated
  and averaged. Diagnostic FITS stacks preserve each central response, frozen comparison, and best-scaled residual.
- `run_klip_adapted_grid_validation.sh` measures paired central differences at four uniformly spaced angles on all 15
  accepted 3.6-pixel radial nodes, omitting any sample within five pixels of the known candidate. With the defaults,
  three of 60 locations are omitted, leaving 57 clear samples and 114 complete KLIP reductions.
- `build_klip_adapted_grid.py` rotates those clear measurements to a common angle, averages them independently at each
  radius with per-pixel validity, linearly interpolates the radial response, and applies the same signed normalized
  filter used by the production KLIP path. It writes a response/filter case consumed unchanged by
  `fit_klip_matched_response.py`, plus an independent comparison with the previously measured candidate response.

From the repository root on ROC, run the initial bounded set with:

```bash
OMP_NUM_THREADS=32 agents/plans/scripts/run_klip_psf_response_experiment.sh \
    science_only reference_dr1_a16 radial_ld_fixed16 radial_ld_arc
```

The current implementation supports FP32 calculation storage, complete non-overlapping annuli, enabled derotation,
no post-median subtraction, no pixel time-series normalization, and `none` or `imageMean` regional centering.

From the repository root on ROC, run the first response-backed inference test with:

```bash
OMP_NUM_THREADS=48 agents/plans/scripts/run_klip_matched_response_validation.sh
```

Run the finite-amplitude basis-adaptation test without repeating the sparse response calculation with:

```bash
OMP_NUM_THREADS=48 nohup agents/plans/scripts/run_klip_finite_response_validation.sh \
  > klip_finite_response_driver.log 2>&1 &
```

Run the bounded paired-central-difference oracle with:

```bash
OMP_NUM_THREADS=48 nohup agents/plans/scripts/run_klip_central_response_validation.sh \
  > klip_central_response_driver.log 2>&1 &
```

Run the complete candidate-avoiding numerical-response grid and matched-response fits with:

```bash
OMP_NUM_THREADS=48 nohup agents/plans/scripts/run_klip_adapted_grid_validation.sh \
  > klip_adapted_grid_driver.log 2>&1 &
```

The paired reductions should take approximately 5.6 minutes at the measured 2.934 seconds per ROC reduction; grid
construction, full-image filtering, and the response-backed fits follow automatically. The primary results are
`adapted_grid/adapted_grid_summary.md` and `klip_adapted_grid_fit_summary.md`. The run is restartable without
overwriting completed reductions or fits. It uses the full fiducial perturbation by default because the central test
found stable response shape and better small-radius numerical behavior there.

The native C++ paired-grid validation uses the same 15 radii, four requested angles, five-pixel candidate avoidance,
and finite-amplitude response as the accepted scripted grid. Each requested polar location selects the nearest clear
exact detector pixel; candidate avoidance therefore substitutes a nearby clear pixel instead of deleting the angular
sample. From the repository root on ROC, run:

```bash
OMP_NUM_THREADS=48 RESPONSE_CASE=radial_ld_refit4_filter \
  EXPERIMENT_DIR=working/roc/klip_cpp_response_$(date -u +%Y%m%dT%H%M%SZ) \
  nohup agents/plans/scripts/run_klip_matched_response_validation.sh \
  > klip_cpp_response_driver.log 2>&1 &
```

The driver requires the native output headers to report 60 retained measurements and 120 signed trial reductions,
checks that response estimation leaves the ordinary science cube unchanged, and repeats the matched-response fit for
every configured KL mode.

The default perturbation magnitudes are 0.25, 0.5, and 1.0 times the P4 exact-negative contrast. Set
`AMPLITUDE_FRACTIONS`, `SAMPLE_SPECS`, or `MODE_COUNTS` only for a deliberately reduced or expanded diagnostic; exact
settings are persisted and a partially completed directory is safely restartable.

## Current gaps

### Response combination semantics

The science residual must continue to use its configured combination. Under the accepted faint-source assumption,
active science `sigmaMean` must map to an unclipped mean for the analytic response while preserving configured
weights, ordinary geometric validity, and `combine.minGoodFract`. Independently sigma clipping a known response is
not the derivative of the science estimator. Product provenance must distinguish the configured science estimator
and threshold from the effective response estimator and threshold.

Median combination remains a nonlinear filter-template approximation. It may continue to retain frame samples until
a separately accepted bounded median approximation exists.

### Retained response memory

Before the bounded accumulator, `KLIPreduction` retained one float response cube and one float validity cube for
every sample-by-mode pair, each with one plane per target frame. That storage scaled as

```text
2 * sampleCount * modeCount * frameCount * stampSize^2 * sizeof(float).
```

For 30 radii, 16 angles, 15 mode outputs, 621 frames, and an 11-by-11 stamp, this is 4,328,121,600 bytes, or
approximately 4.03 GiB, before container overhead and combination scratch.

Mean and weighted-mean responses now retain worker-local weighted response sums plus sample-local weight sums and
valid-frame counts, then reduce workers in fixed index order after all regions. Their dominant response storage scales
with worker count rather than frame count, and the exact retained byte count is persisted. Median intentionally
preserves the frame-stack path until its bounded contract is decided.

### Products and matched filtering

The production writer now optionally evaluates the fitted radial model over each final science mode and applies
P4's signed normalized filter mathematics. It publishes amplitude beside the final image and response-energy
normalization, support, validity, and a completion manifest in the auxiliary-product directory. The maintained
response-backed fitting script consumes those products without rerunning KLIP. Comparison with a full
negative-companion oracle remains the scientific validation step.

### Scientific and performance validation

The maintained ROC driver now compares multiple sparse grids with a fine radial reference and a science-only control,
including response similarity, template-level matched-filter proxies, wall time, peak RSS, retained memory, and
product size. The first representative run is recorded below. A later experiment must compare the accepted
sparse grid with the current negative-companion plus zero-signal fake-injection construction, including matched-filter
position and contrast ranking and uncertainty. Full reinjection/refitting need not agree as a calibrated throughput
measurement, but it must establish that the sparse filter does not materially degrade the accepted inference metrics.

Known-planet avoidance is not initially required for the exact KLIP probe because it does not estimate a response
from the measured residual at that source pixel. The basis is frozen from the original science data, and the probe is
analytic. If validation shows that basis contamination by the known companion matters, compare a neighboring-angle
surrogate or masked-basis construction explicitly rather than silently moving the exact requested probe.

The 2026-09-06 run in `working/roc/klip_psf_response_20260906T163025Z` used commit `55eae39`, 32 workers, a
science-only control, a 1-pixel/16-angle reference with 864 exact probes, and a 2-pixel/16-angle candidate with 432
probes. The response overhead fell from 708.07 to 347.97 seconds (2.03x), while the ordinary science product remained
bitwise unchanged. Across modes, the candidate had median/worst relative L2 errors of 0.0641/0.0644, mean cosine
similarity 0.99794, and worst unit-signal amplitude-proxy error 0.224%. This establishes useful sparsity, but the next
ROC run must evaluate the requested 3.6-pixel radial grid with production linear radial interpolation and compare a
fixed 16-angle grid with the 3.6-pixel maximum-arc grid.

The follow-up run in `working/roc/klip_psf_response_20260906T202414Z` used commit `2c7babc` and completed that
comparison. The 3.6-pixel/fixed-16 grid used 240 measurements, reduced response overhead from 700.18 to 191.73
seconds (3.65x), and retained 35.4 MiB instead of the fine grid's approximately 127.8 MiB. Across KLIP modes its
median/worst relative response errors were 0.10336/0.10391, mean cosine similarity was 0.99467, and the worst
unit-signal matched-filter amplitude proxy error was 0.889%. The 3.6-pixel arc-spaced case expanded to 872
measurements at these radii, provided no accuracy benefit, and cost essentially the same as the fine reference.
The fixed-16 lambda/D-scale grid is therefore the accepted first matched-filter candidate. Its 192-second response
overhead is suitable for the first inference validation and does not justify a detector-local KLIP approximation
before measuring the actual response-backed position and contrast.

The first response-backed inference run in `working/roc/klip_matched_response_20260913T012040Z` used commit
`e78e948`, 48 workers, and the accepted 240-measurement fixed-16 grid. The response-enabled reduction took 177.53
seconds versus 3.00 seconds for the science-only control, retained 53.17 MiB of response accumulators, and increased
peak RSS by 80.05 MiB. Enabling the response left all eight planes of the ordinary science cube elementwise
identical. All eight local fits converged away from their position bounds with full filter-stamp support. From 125
through 350 modes, the fitted separation spanned 11.919--12.141 pixels, PA spanned 261.273--261.623 degrees, and S/N
spanned 5.264--5.476; the across-mode means and sample standard deviations were respectively 12.0359 +/- 0.0716
pixels, 261.441 +/- 0.115 degrees, and 5.360 +/- 0.076. Relative to the P4 exact-negative position, the KLIP
positions differ by 0.218--0.429 pixels and only 0.46--0.75 of each fit's local curvature sigma. This accepts the
production products and the sparse response as a stable detection/localization filter.

The frozen-basis photometry is not accepted. The fitted KLIP contrasts span `0.00114124`--`0.00119344`, with mean
`0.00117027`, or only 0.2396--0.2505 of the P4 exact-negative contrast `0.00476393`. That factor agrees strikingly
with the earlier P4 frozen-response result: the 125-mode KLIP contrast `0.00116746` differs from the P4 frozen sparse
contrast `0.00116759` by only 0.012%. The fixed-16 sampling error relative to dense KLIP was previously below 0.9%, so
radial sparsity cannot explain the factor-of-four scale discrepancy. The leading hypothesis is the same omitted
operator-adaptation term found in P4: the current KLIP response propagates the probe through a basis frozen from the
planet-bearing data, while an end-to-end negative injection changes the reference covariance and KL modes.

The finite-amplitude ROC test in `klip_finite_response_20260913T150323Z` used commit `76f7ed3` and subtracted the P4
exact-negative companion before a complete KLIP refit. The refit took 2.93 seconds with 48 workers and 2.69 GiB peak
RSS. Across all eight mode counts, reconstruction of the persisted frozen-response filter agreed with the production
amplitude to at worst `1.11e-10` in absolute contrast, and both the empirical and frozen stamps had complete support.
The end-to-end removed signal projected onto the frozen response by 0.2271--0.2304, with mean 0.2284. The prior
response-fit contrast fractions were 0.2396--0.2505, with mean 0.2457, so the independently measured response scale
accounts for 90.7--95.5% of the factor-of-four photometric discrepancy. Subtracting the full P4 companion left only
1.38--6.17% of the original frozen-filter amplitude, providing a second direct check that the adopted companion
nearly cancels the signal seen by that statistic.

This confirms that the dominant missing contribution is KL-basis adaptation at this finite companion amplitude and
rules out sparse radial interpolation as its cause. It is not merely a scale correction: the empirical/frozen stamp
cosine is only 0.8023--0.8081 and the best-scaled empirical residual is 0.5891--0.5970. A corrected matched filter must
therefore include the changed response shape rather than multiplying the present frozen response by approximately
four. Because this one-sided difference is a finite-amplitude secant, it does not by itself establish the
infinitesimal KL-basis derivative. The next validation oracle is a bounded paired central-difference test at several
representative sparse samples and perturbation amplitudes. Once its local regime is established, implement and test
the analytic KL-mode perturbation term so production response estimation does not require two complete KLIP refits
per sample.

The paired-central-difference run in `klip_central_response_20260913T154659Z` used commit `a345ed3`, 48 workers, 15
positions, and perturbations of 0.25, 0.5, and 1.0 times the P4 exact-negative contrast. Its 90 complete KLIP
reductions took 264.05 seconds in total, or 2.934 seconds each. At the candidate, the smallest-perturbation derivative
projected onto the frozen response by 0.2204--0.2252 and had cosine 0.7966--0.8017, reproducing the independently
measured finite-response scale and shape. The derivative instead projected onto the previous one-sided finite secant
by 0.9780--0.9838 with cosine 0.9986--0.9999. From the smallest to largest perturbation its mean relative change was
only 2.59%, so the finite secant is already an excellent shape oracle at this source strength and the central test
establishes a usable local response regime.

The clear-sample derivative is strongly radius dependent, with mean projection onto the frozen response near 0.21 at
7.8 and 11.4 pixels, 0.40 at 33 pixels, and 0.54 at 58.2 pixels. This reinforces the need to sample and interpolate in
radius. After rotation to a common angle, response scatter about the radial mean decreases from 56% at 7.8 pixels and
43% at 11.4 pixels to 24% at 33 pixels and 4.4% at 58.2 pixels. The corresponding radial-mean matched-filter behavior
is substantially better than those image-domain residuals suggest: per-angle unit-signal amplitudes span
0.804--1.184 and cosine 0.739--0.952 at 7.8 pixels, 0.823--1.168 and cosine 0.821--0.966 at 11.4 pixels,
0.971--1.028 and cosine 0.969--0.974 at 33 pixels, and 0.992--1.009 and cosine 0.999 at 58.2 pixels. Small-radius
azimuthal averaging is therefore an approximation that must be assessed through filter loss rather than L2 response
error alone.

Most importantly for candidate avoidance, the three clear 11.4-pixel samples predict the independently measured
candidate derivative at 11.74 pixels with projection 1.029--1.043 and cosine 0.952--0.958 for the smallest
perturbation. The result is stable across all three perturbation amplitudes: the mean projection is 1.038--1.042 and
mean cosine is 0.955--0.956. Filtering the original candidate with that clear radial mean yields contrasts
`0.004340`--`0.004704`, with across-mode means `0.004502`--`0.004525`, or 94.5--95.0% of the P4 exact-negative
contrast. This accepts sparse, candidate-avoiding radial averaging as a useful approximate matched filter at the
demonstrated location. A complete 57-sample paired grid is estimated to take about 5.6 minutes, only about 1.9 times
the existing 177.5-second frozen-response calculation. The numerical paired grid is therefore the next practical
candidate estimator; an analytic KL-mode derivative remains a potential production optimization after its accuracy
and runtime are benchmarked against this grid.

The complete adapted-grid run in `working/roc/klip_adapted_grid_20260913T163652Z` retained 57 of 60 requested
locations. Across 125--350 modes, the independently measured candidate response projected onto the clear radial grid
by 1.026--1.039 with cosine 0.9541--0.9557. All eight response-backed fits converged: contrast ranged from
`0.0047386` to `0.0049611`, or 0.9947--1.0414 of the P4 exact-negative contrast, with separation 11.910--12.039
pixels and PA 259.51--260.00 degrees. This accepts the paired candidate-avoiding radial grid as the numerical KLIP
response target. The C++ `refitDifference` implementation now executes that construction inside one `klipReduce`
run; the next ROC test is a direct implementation and performance check, not another algorithm-selection experiment.

The first native runs in `working/roc/klip_cpp_response_20260913T203850Z` and
`working/roc/klip_cpp_response_20260913T210611Z` reproduced the accepted scripted-grid fit closely, but a direct
response-stamp comparison exposed a coordinate-phase discrepancy: native paired trials were injected at continuous
requested polar coordinates and their final stamps were rounded to the nearest pixel. At 200 modes, translating the
native response by about 0.28 pixel improved its cosine with the independently measured candidate response from
0.955 to 0.983, while the registered central 9-by-9 comparison gave cosine 0.991. This translation is diagnostic
only; registering a response after measurement would violate the matched filter's definition at pixel `(x,y)`. The
implementation now instead selects the nearest eligible exact detector pixel before injection, uses its actual radius
and angle for the paired trial, extracts at that same pixel, and rotates by the actual angle before radial averaging.
The next ROC run must validate that coordinate-preserving correction.

## Implementation sequence

### 1. Align response semantics and provenance

- [x] Map active science `sigmaMean` to mean for KLIP analytic responses while retaining weights, geometric validity,
  and minimum-good-frame support.
- [x] Leave the science combination unchanged.
- [x] Record configured science combination/threshold and effective response combination/threshold separately in
  final and response-product FITS headers.
- [x] Add unweighted and weighted active-sigma regressions against direct mean-combined response oracles.

### 2. Bound linear-response storage

- [x] Accumulate mean and weighted-mean regional response contributions in worker-local stamp sums.
- [x] Count valid target frames and accumulate valid weights once per sample/frame, independent of the number of
  regions and mode outputs.
- [x] Reduce worker accumulators after every region is complete and preserve the existing validity and
  `combine.minGoodFract` rules.
- [x] Keep median on an explicitly documented frame-stack path until a bounded nonlinear approximation is accepted.
- [x] Add retained-byte accounting and provenance identifying `WORKER_SUM` versus `FRAME_STACK` accumulation.
- [x] Test multi-region equivalence, weighted support, invalid boundaries, and independence from target-frame count.

### 3. Validate on ROC

- [x] Add maintained runner and analyzer scripts under `agents/plans/scripts` using a response-compatible AF Lep/NACO
  KLIP ROC configuration.
- [x] Compare a fine one-pixel radial reference with multiple radial/angular grids.
- [x] Record response error, matched-filter proxy error, worker/wall time, peak RSS, and stored-product size.
- [x] Decide whether exact sparse sky propagation is fast enough before considering a detector-local approximation.

### 3a. Adopt lambda/D-scale sampling

- [x] Resolve per-radius angular counts from either a fixed count or a maximum circular arc step.
- [x] Average all common-angle samples independently at each configured radius.
- [x] Replace nearest-radius lookup with linear interpolation, intersecting endpoint validity and clamping only
  outside the sampled interval.
- [x] Repeat the P4 and KLIP ROC experiments at 3.6-pixel radial spacing with both fixed and arc-spaced angular grids.

### 4. Add matched-filter and fitting integration

- [x] Evaluate the fitted radial response at every requested science location and apply the shared signed normalized
  matched-filter convention.
- [x] Publish filtered value, normalization, support, and validity products with deterministic coordinates and mode
  identity.
- [x] Expose a bounded response-backed position/contrast merit calculation so a fit does not rerun KLIP merely to
  refresh its filter.
- [x] Run the first response-backed fit across the configured AF Lep mode counts and record detection, localization,
  photometry, product-validity, and performance behavior.
- [ ] Validate recovered position, contrast ranking, and curvature- or resampling-derived uncertainty against the
  negative-fit/zero-signal-injection oracle.

### 5. Include KL-basis adaptation

- [x] Measure the finite-amplitude end-to-end KLIP response at the accepted companion and compare its scale and shape
  with the frozen-basis sparse response.
- [x] Measure paired central differences at representative radius/angle samples and multiple perturbation amplitudes
  to identify the locally linear regime and provide an implementation oracle.
- [x] Run the complete candidate-avoiding paired grid, construct its radial response, and repeat matched-filter
  photometry for every configured KL mode.
- [x] Implement the paired candidate-avoiding grid, region-aware radial nodes, common-angle averaging, radial
  interpolation, and filtering in the KLIP C++ reduction path using the shared `[psfResponse]` interface.
- [x] Center every KLIP response measurement on the nearest eligible exact detector pixel and retain its actual
  radius/angle through injection, extraction, and common-angle rotation.
- [ ] Validate the native C++ paired-grid products and runtime on ROC against the accepted scripted-grid result.
- [ ] Add the analytic covariance/eigenmode perturbation contribution to the sparse KLIP response calculation.
- [ ] Validate the adapted sparse response against the paired-refit oracle before repeating matched-filter photometry
  and uncertainty tests.

## Acceptance criteria

- Enabling KLIP PSF estimation does not change the ordinary science residual or final science combination.
- Mean, weighted mean, and active science `sigmaMean` produce the documented frozen-basis analytic response with
  truthful estimator provenance.
- Linear-response retained memory is independent of target-frame count and is predicted before allocation.
- Sparse response products remain deterministic in sample/radius/mode ordering and retain explicit validity.
- A ROC comparison shows that the sparse response provides an acceptable matched filter without the second full
  reduction in the fitting loop.

## Known non-blocking follow-ups

- Median response combination still needs a bounded accepted approximation or external-storage strategy.
- Median regional centering, pixel time-series normalization, wedges/overlapping regions, and post-median subtraction
  require separate perturbation contracts.
- The 2026-09-06 radial-sampling follow-up rechecked `/home/jrmales/Source/mxlib/_build/coverage_filtered.info` for
  every mxlib API directly called by the edited KLIP configuration, preparation, region, and header functions. The
  exact app-configuration, exception, finite-check, FITS file/header, `eigenCube<float>`, `radAngImage`,
  `annulusIndices`, `cutImageRegion`, degree/radian, and time-utility paths have 100% executable-line coverage. No new
  mxlib ownership follow-up is required.
- The 2026-09-12 matched-filter integration rechecked the same current LCOV report for every mxlib API called by the
  edited configuration, validation, header, response-product, and P4 product-writing functions. The exercised
  `appConfigurator::add` and scalar/vector extraction paths, exception construction, `math::isFinite<float>`,
  `invalidNumber<float>`, `ioutils::createDirectories` and `parentPath`, `eigenCube<float>` construction/access,
  FITS header append, and FP32 image/cube FITS write paths all have 100% executable-line coverage. The containing
  `fitsFile.hpp`, `fitsHeader.hpp`, `eigenCube.hpp`, and `floatUtils.hpp` reports are respectively 491/491, 259/259,
  185/185, and 20/20 executable lines. No mxlib ownership follow-up is required.
- The 2026-09-13 exact-pixel sampling change rechecked the current
  `/home/jrmales/Source/mxlib/_build/coverage_filtered.info` report for the mxlib APIs called by the edited preparation
  and test paths. `exception.hpp`, `fitsFile.hpp`, `eigenCube.hpp`, `eigenImage.hpp`, `floatUtils.hpp`, and `geo.hpp`
  respectively report 42/42, 491/491, 185/185, 4/4, 20/20, and 16/16 executable lines. No mxlib ownership follow-up
  is required.
