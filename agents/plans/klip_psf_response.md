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
- `src/common/KLIPreduction.hpp` configures exact sky samples at `klip.psfSampleRadii` and uniformly spaced angles,
  invokes the probe while each target-specific basis is resident, combines target-frame stamps, fits one
  `RadialPSFModel` per requested mode count, and optionally writes canonical radial response and validity cubes.
- `src/common/RadialPSFModel.hpp` and `.cpp` rotate measurements to the positive-column radial orientation, average
  samples at the same radius, and evaluate arbitrary positions with nearest-radius interpolation.
- `tests/common/KLIPPSFModel_test.cpp`, `tests/common/KLIPreduction_test.cpp`, and
  `tests/common/RadialPSFModel_test.cpp` cover probe preprocessing, frozen-basis projection, compact derotation,
  radial alignment, nearest-radius evaluation, orchestration, and persisted schema-1 products.
- `doc/klip.dox`, `doc/klip_algorithm.dox`, and `doc/klip_config.dox` describe the public controls and initial
  restrictions.

The existing configuration surface is:

```ini
[klip]
psfFile=/path/to/centered_psf.fits
psfStampSize=11
psfSampleRadii=1,3,5,7
psfSamplesPerRadius=16
outputPSFModels=true
psfOutputPrefix=klipPSF_
```

Maintained ROC experiment assets are under `agents/plans/scripts`:

- `klipReduce_afLepNaco_psf_response.conf` is the response-compatible AF Lep/NACO KLIP base configuration. It follows
  `working/kr.conf`, uses the current ROC angle constant, and deliberately selects supported `imageMean` centering.
- `run_klip_psf_response_experiment.sh` runs a science-only control, a 1-pixel/16-angle radial reference, and several
  coarser radial/angular grids while recording exact commands, hashes, wall time, peak RSS, and completion state.
- `compare_klip_psf_response.py` compares canonical response/validity cubes by nearest radius, verifies the response
  run leaves the final science cube unchanged, and reports template-level matched-filter amplitude/cosine proxies.
  Those proxies do not replace the later negative-fit/zero-signal-injection inference oracle.

From the repository root on ROC, run the initial bounded set with:

```bash
OMP_NUM_THREADS=32 agents/plans/scripts/run_klip_psf_response_experiment.sh \
    science_only reference_dr1_a16 radial_dr2_a16
```

The current implementation supports FP32 calculation storage, complete non-overlapping annuli, enabled derotation,
no post-median subtraction, no pixel time-series normalization, and `none` or `imageMean` regional centering.

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

The current writer publishes only canonical radial response and validity cubes. It does not apply a KLIP matched
filter to the final science modes, publish filter normalization/support products, or expose a response-backed merit
function to a position/contrast optimizer. Prefer sharing P4's signed normalized filter mathematics and output
contract where the KLIP radial product supplies the local template.

### Scientific and performance validation

The maintained ROC driver now compares multiple sparse grids with a fine radial reference and a science-only control,
including response similarity, template-level matched-filter proxies, wall time, peak RSS, retained memory, and
product size. The representative run still needs to be performed. A later experiment must compare the accepted
sparse grid with the current negative-companion plus zero-signal fake-injection construction, including matched-filter
position and contrast ranking and uncertainty. Full reinjection/refitting need not agree as a calibrated throughput
measurement, but it must establish that the sparse filter does not materially degrade the accepted inference metrics.

Known-planet avoidance is not initially required for the exact KLIP probe because it does not estimate a response
from the measured residual at that source pixel. The basis is frozen from the original science data, and the probe is
analytic. If validation shows that basis contamination by the known companion matters, compare a neighboring-angle
surrogate or masked-basis construction explicitly rather than silently moving the exact requested probe.

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
- [ ] Compare dense or full-injection references with multiple radial/angular grids.
- [ ] Record response error, matched-filter error, worker/wall time, peak RSS, and stored-product size.
- [ ] Decide whether exact sparse sky propagation is fast enough before considering a detector-local approximation.

### 4. Add matched-filter and fitting integration

- [ ] Evaluate the fitted radial response at every requested science location and apply the shared signed normalized
  matched-filter convention.
- [ ] Publish filtered value, normalization, support, and validity products with deterministic coordinates and mode
  identity.
- [ ] Expose a bounded response-backed position/contrast merit calculation so a fit does not rerun KLIP merely to
  refresh its filter.
- [ ] Validate recovered position, contrast ranking, and curvature- or resampling-derived uncertainty against the
  negative-fit/zero-signal-injection oracle.

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
- Apply the repository mxlib coverage gate to every edited function that calls mxlib and record any ownership gap in
  `agents/plans/mxlib_cleanup.md`.
