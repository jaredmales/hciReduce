# P4 PSF Calculation and Post-Processing Plan

## Goal

Calculate a spatially variable forward-model PSF for every valid P4 output pixel and mode, then optionally filter the
final image with the PSF assigned to that pixel. The calculation must reuse the P4 model already learned from the
science data, remain bounded for thousands of frames, preserve the existing reduction when disabled, and distinguish
an analytic frozen-model response from a complete fake-source reinjection and refit.

## Scientific contract to establish first

The recommended initial definition is a **frozen P4 forward model**:

- learn the ordinary P4 coefficients from the science sequence exactly once;
- hold those coefficients fixed while propagating a supplied reference PSF template at its stored normalization;
- apply the same predictor interpolation, temporal-image selection, detector-to-sky transformation, validity, and
  linear combination weights used for the science residuals; and
- preserve the response's positive core and any negative subtraction lobes without taking an absolute value, clipping
  negative values, or renormalizing its integral or peak.

QUESTION: what does "signed response" mean?

ANSWER: P4 subtraction can give its effective PSF negative sidelobes. "Signed" means retaining those negative values
along with the positive core; it does not refer to a separate sign convention or a positive/negative source option.

This measures algorithmic attenuation and distortion for the model that produced the science image. The forward PSF
retains the template's scale; a later normalized matched filter estimates a multiplicative amplitude relative to that
template. The calculation does not
include the change in the learned PCA basis or coefficients that would occur if a finite-brightness source were added
before fitting. Existing fake injection followed by a complete rerun remains the validation path for that nonlinear
effect and for inherited preprocessing.

CONFIRMED: the only goal is to calculate an estimate of the PSF for optimal filtering.

RESOLUTION: the calculated PSF and any persisted PSF field are filter templates, not photometric-throughput
calibrations. Full fake injection is useful for validating detection/filter behavior, but matching its calibrated flux
response is not a requirement of this task.

The supplied PSF template should initially be defined at the input to P4, after any inherited preprocessing. P4 must
not silently apply the science cube's radial-profile, median-filter, mean-image, or temporal-RMS operations to a single
template: several are data-dependent or nonlinear, and a single centered image is insufficient to reproduce their
time-cube response. The template's units, normalization, center, and stage must be recorded in output provenance.

### Combination boundary

The final estimator determines whether one linear forward-model PSF exists:

| Final operation | Frozen-model interpretation |
|---|---|
| Unweighted or weighted mean | Linear; use the same per-frame validity and normalized weights as the science image. |
| `sigmaMean` with nonpositive threshold | The current implementation falls back to mean and is linear. |
| Active `sigmaMean` | Apply the configured sigma-clipped mean directly to the propagated PSF samples as the filtering approximation. |
| Median combination | Median-combine the propagated PSF samples, matching the configured science estimator as the filtering approximation. |
| `adi.postMedSub=true` | Nonlinear temporal median step before derotation. Reject initially for analytic PSF calculation. |

COMMENT: we're focused on the regime where noise >> signal in individual images, so I think we approximate the effect of the PSF perturbation on the combination as linear regardless of which combo we use.  That is, we'll just median combine it if we're using median combination.  This is justified b/c we're only using this as a filter, not a calibration.

RESOLUTION: accepted. The PSF response uses the same configured final estimator, including median or active
sigma-clipped mean, without attempting to reuse the science data's median selection or sigma-clip membership. Output
provenance describes this as a filtering approximation rather than a calibrated linear response.

The ordinary science reduction and analytic PSF calculation both support every existing final-combination method.
The initial restriction on `adi.postMedSub=true` remains separate because that operation modifies each detector pixel's
temporal sequence before derotation rather than only combining already aligned samples.

## Initial findings

### The required forward model is not retained

For one search pixel, P4 builds a target time series `y` and predictor matrix `X`, then writes residuals for requested
mode counts. `P4PCAResult` currently returns those residuals and rank support but discards the predictor-space
coefficient vectors.

For a supported truncation `m`, the compact learned model is

```text
beta_m = (truncated pseudoinverse of X) * y
r_m    = y - X * beta_m.
```

The existing eigensolve already has everything needed to accumulate `beta_m` without another decomposition:

```text
predictor-Gram mode q:  delta_beta = q * (q^T X^T y) / lambda
temporal-Gram mode u:   delta_beta = X^T u * (u^T y) / lambda
```

For rotated-frame centered fitting, `X` and `y` in these formulas are centered, while the resulting `beta_m` is
applied to the original uncentered predictor samples. Coefficient output must be optional so ordinary P4 runs do not
allocate an additional `K * M` doubles per worker.

### Detector-frame response with same-image predictors

For detector search pixel `p`, mode `m`, a post-preprocessing PSF template `F`, and same-image predictor samples `q_k`,
the frozen local response is the target template minus the learned predictor combination. In schematic continuous
coordinates,

```text
G_p,m(delta) = F(delta) - sum_k beta_p,m,k * F(delta + q_k - p).
```

Production evaluation must use the exact stored P4 interpolation records and cubic kernels rather than treating
fractional `q_k` as integer shifts. A small explicit injected-template calculation will define sign, row/column,
center, and interpolation conventions.

With `p4.numberImages=0` and one stable template, this local stamp is independent of target frame. It can be calculated
as soon as one search-pixel regression completes and retained as `stamp pixels * requested modes` floats rather than
retaining its coefficient vectors.

### Additional temporal images make the response frame-dependent

With `p4.numberImages>0`, a predictor row contains the same-image OR plus PSF-disk samples from selected earlier and
later images. The selected image identities vary with the central target row, including one-sided endpoint fallback.
A sky-stationary source occupies a different detector coordinate in each selected image. Consequently, one
frame-independent stamp per detector search pixel is not sufficient.

The bounded extension should retain the compact `beta` vectors for that annulus, synthesize one target frame's local
responses at a time using `m_temporalSelections`, derotate/accumulate them, and release the frame scratch. It must not
allocate `search pixels * target frames * modes * stamp pixels` by default.

### Output-pixel reconstruction needs an exact oracle

The final PSF at a sky pixel is not obtained by merely rotating one nearest detector stamp. The science residual
derotation uses cubic convolution, so a sky output pixel receives contributions from a detector-frame footprint and
from the validity associated with those pixels. The local PSF coordinate axes must also rotate into the sky frame.

The authoritative small-case oracle should therefore:

1. freeze explicitly known P4 coefficient vectors;
2. construct a detector cube containing a unit PSF at one fixed sky location in every frame;
3. apply the frozen local P4 operator directly at every required detector search pixel;
4. run the existing validity-aware derotation and combination; and
5. crop the resulting final response about the source position.

The compact stamp-field reconstruction must match this oracle. This test determines the required transpose, stamp
rotation, fractional shift, and detector-field interpolation order without relying on an informal convolution
argument.

The first nonzero-rotation oracle exposed an additional sampling-phase constraint: an odd local response stamp cannot
represent an even-centered template on the template's native integer sample phase (and vice versa). The retained
local response dimensions therefore preserve the corresponding template row/column parity and include enough radial
padding for both cubic footprints at any rotation. With mixed template parity (15-by-16), mixed local parity
(17-by-18), an asymmetric PSF, spatially varying coefficients, and a 0.23-radian rotation, the compact reconstruction
matched direct full-detector injection to about `6e-7` absolute float precision. The production reconstructor now
composes the two cubic footprints explicitly and propagates invalid local models instead of interpreting them as zero.

### Memory and output scale

Let `S` be the number of owned P4 search pixels, `K` the predictor count, `M` the requested output count, and `A` the
number of pixels in one PSF stamp.

Potential retained representations have leading sizes:

```text
all coefficients:                 8 * S * K * M bytes
frame-independent local stamps:   4 * S * M * A bytes
all frame-dependent stamps:       4 * S * T * M * A bytes  (prohibited by default)
one-frame stamp scratch:          4 * S * M * A bytes
one final compact PSF mode:       4 * S * A bytes
```

For the measured 2,828-pixel, 11-mode geometry and a 10-by-10 stamp, frame-independent stamps occupy about 11.9 MiB,
while 8,182-frame stamps would occupy about 95 GiB. If `K=800`, all-double coefficients occupy about 190 MiB. These
terms must be added to the existing compact-residual and worker estimates before enabling the feature.

### AF Lep/NACO test reference

The available integration reference is
`/home/jrmales/Source/mxWork/NACO/AFLep/2011-10-21/out/psf_reg_median.fits`. It is a header-minimal, 256-by-256,
single-precision FITS image with all 65,536 pixels finite. Its maximum is at zero-based `(127,127)`, and all four
pixels surrounding the 127.5-by-127.5 geometric center are bright. The file has no normalization or center metadata;
its stored sum is approximately 31,264.21, so tests must preserve the stored scale and explicitly use the geometric
center rather than infer unit flux or a peak-pixel center. Its current SHA-256 is
`bc43daa82a433e32a756af0ba82cd0061c8d3f01adb73bc2fed97ed832880194`.

Focused unit tests must remain self-contained and use synthetic asymmetric templates. The AF Lep file is for local and
remote real-data integration/oracle comparisons and should be supplied through configuration rather than embedded as
an absolute production default.

### Initial implementation measurement

The first coefficient-and-local-stamp implementation was measured on 2026-08-21 with the 96-frame, one-annulus
`working/p4Reduce_afLepNaco_profile_local.conf` configuration, 20 OpenMP workers, one OpenBLAS thread, the AF Lep PSF
above, an 11-by-11 stamp, and three output modes. The PSF-enabled run retained 0.130 MiB of local stamps, 0.258 MiB
of prepared template state, estimated 2.034 MiB per worker, and reached 311,968 KiB maximum RSS. The disabled run
estimated 1.995 MiB per worker and reached 311,732 KiB maximum RSS. On this small case the PSF work added about
0.40 seconds of wall time; its 6.11 aggregate worker seconds were 91.6% of the regression worker total. Both the
enabled and disabled runs completed successfully after rebuilding the application against the changed library ABI.

This result does not justify temporary products for the frame-independent local state: the retained model is tiny
relative to the science cube and materialization. Final-field reconstruction and the temporal-image extension still
need their own peak-memory measurements before the no-temporary-file decision can be generalized.

After enforcing phase parity and rotation-safe support, the same local run used 26-by-26 internal response stamps.
It retained 0.743 MiB of local responses, estimated 2.036 MiB per worker, and reached 312,284 KiB maximum RSS versus
311,636 KiB for a contemporaneous disabled run. Local-PSF calculation took 26.43 aggregate worker seconds and added
about 1.90 seconds of wall time. The larger internal stamp fixes the interpolation error at a modest memory cost, but
also confirms that local-response calculation is the first performance target once the full-field algorithm is
correct.

The initial `+4` radial-support rule was subsequently shown to be conservative. Exact enumeration of the composed
four-by-four PSF-shift and derotation footprints over arbitrary subpixel phase gives a minimum 22-by-22 local response
for an 11-by-11 output stamp and an even template, or 23-by-23 for an odd template. The parity-aware production bound
now uses those exact dimensions; the historical timing and memory measurements above remain measurements of the
earlier 26-by-26 implementation.

A broader four-annulus `[6,14)` output run used 96 frames, an 11-by-11 final stamp, three modes, and 20 workers. It
reconstructed `S=504` final positions in 0.262 seconds after a 13.568-second regression. Local response storage was
0.003809 GiB, maximum RSS was 372,580 KiB, and each 248-KiB mode cube had 25,936 finite samples. Exact source-center
validity was 292 of 504 positions in every mode. The omitted inner/outer calculation halo accounts for the remaining
invalid positions and confirms that a narrow isolated annulus cannot provide complete derotation-footprint support.

## Proposed configuration and products

Use opt-in P4-specific configuration so all existing controls and outputs remain unchanged when no PSF template is
configured:

| Proposed target | Contract |
|---|---|
| `p4.psfFile` | Centered FITS template at the post-preprocessing P4 input stage; unset disables analytic PSF calculation. |
| `p4.psfStampSize` | Positive square stamp width. Both odd and even sizes are allowed; the center is the geometric `(size-1)/2` coordinate. |
| `p4.outputPSFModels` | Write the compact spatially variable PSF products in addition to any in-process filtering. |
| `p4.psfFilter` | Apply the spatially variable PSF filter to the final science cube; default `false`. |
| `p4.psfFilterMinGoodFract` | Minimum usable fraction of the full odd-sized filter stamp; default `1`. |
| `p4.psfOutputPrefix` | Prefix for compact PSF-model products and their manifest inside the common output directory. |

The names and grouping should be finalized during implementation review rather than overloading `fake.fileName` or
`p4.psfRadius`. The physical exclusion radius and the forward-model stamp size are different quantities.

The recommended compact persisted schema is:

- one `stamp rows * stamp columns * owned output pixels` FITS cube per P4 output mode;
- one coordinate image/table mapping the deterministic plane index to final `(row,column)`;
- one mode-by-output-position validity product; and
- a schema version plus PSF source, template center, stamp size, regression frame, realized mode, combination policy,
  normalization, and ownership metadata in every product.

Writing one mode at a time avoids retaining the complete `M * S * A` final field. The ordinary `finim` remains the
unfiltered P4 product. Filtering should first write a distinct, clearly named cube so science and filtered products
can be compared; replacing the in-memory/output `finim` can be considered only as a separate explicit option.
Implemented filter products derive their names from the resolved ordinary final-image path, insert the product role
before its shared four-digit sequence (or before the extension for an exact name), and mirror the complete final-image
FITS header before appending product-specific P4 provenance. Compact model products continue to use
`p4.psfOutputPrefix`.

## Work sequence

### 1. Lock the frozen-model numerical seam

- [x] Add an optional coefficient result to `P4PCA` without changing the allocation or arithmetic path when it is not
  requested. Preserve mode order and use the existing rank status for unsupported coefficient columns.
- [x] Accumulate predictor coefficients in the existing retained-mode loop for both temporal- and predictor-Gram
  branches. Do not form a dense `K * K` pseudoinverse or projector.
- [x] Cover uncentered, centered, temporal-Gram, predictor-Gram, rank-insufficient, zero-rank, and eigensolver-failure
  cases. Verify `X * beta_m == y - residual_m` at accepted all-double tolerances, including centered fitting applied to
  the original uncentered predictors.
- [x] Measure the opt-in `K * M` worker allocation and update `estimatedWorkerBytes()` and the memory-policy report
  before integrating stamps.

### 2. Define and test the local PSF operator

- [x] Implement a pure numerical/geometry component that accepts a centered PSF template, one search pixel's exact
  predictor geometry, one supported coefficient vector, and a requested stamp grid.
- [x] Evaluate fractional predictor locations with the same cubic-convolution convention used by P4 sampling. Define
  and test row/column sign, geometric centers for odd and even templates/stamps, boundary support, and zero padding.
- [x] Build a brute-force detector-image oracle by shifting/injecting the template, sampling it through
  `P4PixelGrid`, and applying the same frozen coefficients. Compare every stamp element to the compact operator.
- [x] Add an AF Lep/NACO integration case using `psf_reg_median.fits` at its stored normalization and geometric center;
  record a content checksum with the test result so a changed external reference is not mistaken for a regression.
- [x] Require finite template values and checked dimensions. Reject a requested stamp whose required template support
  cannot be evaluated under the selected padding contract.
- [x] Initially implement and validate the frame-independent `numberImages=0` path; keep positive `numberImages` an
  explicit rejection until the streamed temporal path below is complete.

### 3. Capture local models during detector-frame regression

- [x] Enable analytic PSF calculation only for detector-frame P4 in the first implementation. Rotated-frame P4 remains
  a numerical reference and should reject this option until its forward-model and final-frame meaning are separately
  accepted.
- [x] For every valid, rank-supported search pixel, calculate compact local stamps while `beta` is resident and write
  them into the annulus's deterministic owned-pixel storage. Invalid/rank-insufficient modes receive explicit
  validity rather than zero-valued PSFs.
- [x] Preserve annulus ownership and non-overlap rules. Do not retain a stamp for pixels that cannot contribute to a
  valid science output.
- [x] Add local-stamp storage and scratch to checked memory calculations. Allow a configured run to fail before
  regression if compact residuals, one residual materialization, local PSF state, and one worker cannot fit.
- [ ] If local stamps plus residuals become material, write completed annulus/model blocks to a versioned temporary or
  final FITS product and release them. Use explicit completion metadata so interrupted products are not treated as
  complete.

### 4. Reconstruct one final sky-pixel PSF exactly

- [x] Implement the small direct frozen-injection oracle first for a few frames, pixels, and one mode.
- [x] Construct the compact reconstruction using the exact inverse derotation coordinates, cubic footprint weights,
  local-stamp coordinate rotation, science validity, and combination weights.
- [x] Match the oracle for zero and nonzero angles, integer and fractional inverse coordinates, mask/annulus boundaries,
  weighted and unweighted means, odd/even image sizes, and odd/even stamp sizes.
- [x] Define final PSF validity from the same `combine.minGoodFract` and effective frame support as the corresponding
  science output pixel. Do not treat an invalid local model as a zero response.
- [x] Validate source positions at multiple radii and position angles so a transpose or sign error cannot pass a
  symmetric on-axis test.

### 5. Generate the complete compact spatially variable PSF field

- [x] Process one P4 mode and a bounded tile of final output positions at a time. Retain at most the local input-model
  store, one output PSF tile/mode, and bounded frame scratch.
- [x] Parallelize over independent final output positions or tiles only after a serial oracle match. Respect the P4
  memory budget and prevent nested BLAS/OpenMP oversubscription.
- [x] Write the compact mode cube, coordinate mapping, validity, and complete provenance transactionally. Verify that
  plane ordering is deterministic across worker counts.
- [x] Add PSF timing and storage estimates to the high-level timing report only after measuring meaningful stages;
  version the timing diagnostic if its persisted columns change.

### 6. Support `p4.numberImages>0` without frame-count storage

- [ ] Retain or stream the supported `beta` vectors for one annulus instead of a frame-independent local stamp.
- [ ] For each retained central target row, use its exact `m_temporalSelections` entry and predictor-column ordering to
  evaluate the source in the central and selected neighboring detector images. Include one-sided endpoint selection
  and the annulus's realized temporal exclusion radius.
- [ ] Accumulate one target frame directly into bounded final-PSF tiles, then release its response scratch. Never
  allocate `S * T * M * A` stamps.
- [ ] Compare against a direct frozen-model cube containing a rotating unit source for `numberImages=1` and a larger
  `N`, including beginning/end targets and non-monotonic angle sequences.
- [ ] Recalculate worker and retained-memory estimates for coefficient storage; allow writing/reloading coefficient
  blocks by annulus if that is smaller than retaining all regions.

### 7. Apply the configured combination approximation

- [x] Route propagated per-frame PSF samples through the configured mean, weighted mean, median, or sigma-clipped mean
  implementation with the same validity, weights, sigma threshold, and `combine.minGoodFract` settings.
- [x] For active `sigmaMean`, calculate clipping from the propagated PSF samples themselves. Do not retain or reuse the
  science residual's clip membership.
- [x] Document median and active-sigma results as filter-template approximations justified by the individual-frame
  noise-dominated regime, not as finite-amplitude-independent photometric responses.
- [x] Keep `adi.postMedSub=true` rejected until its separate detector-frame temporal operation is explicitly accepted
  under a comparable filtering approximation.

### 8. Add spatially variable PSF filtering

- [x] Define the filter mathematically before naming the output. The implemented flux-estimator form at output pixel
  `q` is the validity-weighted local matched filter

  ```text
  a_hat(q) = sum_delta H_q(delta) * I(q+delta)
             / sum_delta H_q(delta)^2,
  ```

  with an optional raw-correlation diagnostic. The denominator and edge/invalid support must be reported rather than
  silently absorbed.
- [x] Apply one compact local PSF to the corresponding final-image neighborhood without constructing a dense
  space-variant convolution matrix. Process modes and output tiles independently.
- [x] Define a minimum usable PSF fraction and output validity at image/annulus boundaries. Do not replace invalid
  filtered pixels with zero.
- [x] Write the filtered cube separately with the source final-image identity, PSF schema/version, normalization,
  support fraction, and mode mapping. Preserve the original final image by default.
- [ ] Validate filter response, localization, and detection ranking with isolated synthetic sources, multiple
  contrasts, neighboring sources, and noise-only data. Compare against full fake injection and refitting to determine
  whether the estimated PSF materially improves filtering; calibrated throughput agreement is not required.

### 9. Documentation, provenance, and coverage

- [x] Add the accepted configuration to `doc/p4_config.dox` and the frozen-model equations, limitations, product
  schema, and filtering definition to `doc/p4.dox`.
- [x] Add expressive FITS provenance cards for the template, template stage, stamp geometry, PSF schema, frozen-model
  convention, combination policy, temporal-image support, filter normalization, and validity thresholds.
- [ ] Add application-level real-FITS tests proving that the feature is opt-in, output directories are created, errors
  are actionable, and disabled runs retain their existing data and headers.
- [x] Run the complete P4, shared ADI/HCI, application, and documentation suites. Compare final science products from
  disabled and enabled-but-nonreplacing runs.
- [x] For every edited function containing mxlib calls, verify all directly called APIs in the current mxlib LCOV
  report. Record any gap under `Known non-blocking ownership follow-ups` in `agents/plans/mxlib_cleanup.md`.

The 2026-08-21 filtered-header/naming follow-up rechecked the current mxlib LCOV report. The exact executable lines for
`getSequentialFilename`, `parentPath`, and `createDirectories` are covered, while the called FITS file/header,
`ompLoopWatcher`, finite-check, and time-utility implementation files have 100% executable-line coverage. No new mxlib
ownership follow-up is required.

## Acceptance criteria

- Enabling PSF calculation does not change the science residual or final-image data, validity, headers, or filenames
  unless a separately explicit replacement option is selected.
- P4PCA coefficient reconstruction matches the existing residual calculation in both Gram branches and centered mode.
- Compact local and final PSFs match the direct frozen-injection oracle across asymmetric geometry, nonzero rotations,
  fractional interpolation, masks, rank failure, weights, and temporal-image selections.
- Peak memory is predicted and bounded. No default path retains frame-dependent stamps, and worker selection includes
  all opt-in coefficient/stamp scratch.
- Every compact PSF plane has a deterministic output coordinate, mode, validity state, template identity, and schema
  version sufficient for later filtering without rerunning regression.
- Every configured final-combination method is applied to the propagated PSF with the science configuration's
  validity, weights, sigma threshold, and minimum-good-fraction settings. Median and active-sigma products are labeled
  as filtering approximations.
- The spatially variable filter has a documented normalization and invalid-support rule and recovers accepted
  synthetic-source flux/position tolerances.
- A representative remote run records PSF calculation time, peak RSS, persisted product size, and filtering behavior
  on full-reinjection sources before the feature is considered scientifically accepted.

## Deferred decisions

- Post-median subtraction still requires an explicit filtering approximation; final median combination is supported
  by directly median-combining the propagated PSF samples.
- Modeling inherited preprocessing belongs to a future cube-level forward operator. The initial PSF template is an
  input-stage response, and full fake reinjection remains the end-to-end comparison.
- Per-image/list PSF templates, wavelength-dependent templates, and per-frame flux scales should follow the stable
  single-template path.
- GPU acceleration is premature until the serial oracle, storage layout, and CPU tile profile are established.
- Replacing the ordinary `finim` with the filtered cube remains opt-in and separate from producing a comparison
  product.
