# P4 Pixel-Local Processing Plan

## Goal

Add an exact, bounded P4 evaluation path for a small sky region around a trial companion. The path will refit every
P4 regression needed by that region after adding a trial source to the already-preprocessed input samples, but will
not process or retain the remainder of the output field. Its first scientific use is repeated negative-companion
evaluation: load the prepared observation once, vary separation, position angle, and contrast, and drive the
local residual toward the corresponding signal-free result.

This is a finite-amplitude refit, not the frozen-coefficient response implemented by the P4 PSF calculation. It must
therefore capture trial-source cross-talk into both the target time series and the learned predictor basis.

The pixel-local evaluator should stop at a reusable in-memory result. Selecting an optimizer and merit function is a
separate layer. After an optimum is found, one ordinary full-field negative-injection reduction remains the final
validation and product-generation step.

## Motivation from the initial P4 filter result

The first spatially variable P4-filter test showed a radial structure in the filtered amplitude that is absent from
the ordinary P4 image after comparable Gaussian smoothing. At the companion radius, the raw filter correlation did
not have a corresponding excess, while the response-energy denominator had a trough. Dividing by that denominator
amplified the existing residual structure. A `sum(H*I)/sqrt(sum(H*H))` detection-statistic proxy reduced the visible
effect only slightly and did not outperform the ordinary smoothed image.

A pure frozen-template autocorrelation did not reproduce the structure. The next decisive scientific comparison is
therefore a P4 model learned from a signal-free sequence. An efficient negative-companion fit is the practical route
to that state, and P4's independent detector-pixel regressions make a local finite-amplitude refit possible.

## Scientific contract

### Exactness boundary

For a fixed post-preprocessing science cube, trial template, P4 configuration, derotation sequence, and final
combination policy, the pixel-local result must equal the matching pixels cropped from a complete P4 rerun in which
the same trial template was added at the P4 input stage. “Local” changes only which independent regressions and
residual samples are evaluated; it does not freeze coefficients, linearize the PCA fit, reduce the temporal sample
set, or approximate the predictor interpolation.

For detector search pixel `p`, input frame `t`, and trial perturbation `D_t`, the local path must build

```text
y'_p(t)   = I_t(p) + D_t(p)
X'_p(t,k) = S_p,k[I_t] + S_p,k[D_t]
```

using the same P4 sampling records as the full reduction, and then run the ordinary `P4PCA` calculation on `y'_p`
and `X'_p`. With positive `p4.numberImages`, every selected neighboring-image sample receives the perturbation from
that neighboring physical frame. This is what allows the source to change the local eigensystem and coefficients.

The initial production path should support detector-frame P4. Rotated-frame P4 is not the current science path and
has different locality: a fixed sky search pixel is already sampled through every detector frame. It should be
rejected until a separate contract and need are established.

### Preprocessing boundary

Ordinary `[fake]` injection currently occurs in `ADIobservation::postReadFiles()`, before coaddition and shared
preprocessing. P4 regression consumes the resulting post-preprocessing cube. Several supported preprocessing steps
are nonlinear or data-dependent, including median filtering, radial-profile subtraction, temporal median
subtraction, and temporal RMS normalization. Adding a cropped template after those operations is not generally
equivalent to injecting before them.

The exact initial local contract should therefore require `preProcess.skip=true`: inputs must already be coadded and
preprocessed, and the configured trial template must represent that same P4-input stage. The current inherited fake
hook does not inject when preprocessing is skipped, so the local evaluator can reuse `[fake]` position, contrast,
template, and optional scale settings without double injection. The output must record that the injection stage is
`POST_PREPROCESS_P4_INPUT`.

Support for an internally preprocessed local perturbation can be considered only after each enabled preprocessing
operation has an explicit response rule. It must not be inferred silently from the science configuration.

### Trial and result scope

The numerical seam should evaluate one trial source at a time and return:

- the final combined local residual stamp for every configured P4 mode fraction;
- an exact validity stamp for every mode;
- the integer full-image origin and continuous trial center of the stamp;
- the trial separation, position angle, contrast, template identity, and per-frame scaling provenance; and
- timing and geometry counts sufficient to measure the speedup.

The first command-line integration may iterate over the existing equally sized `fake.sep`, `fake.PA`, and
`fake.contrast` vectors, but each tuple must be evaluated independently against the same pristine post-preprocessing
cube. Trial perturbations must never accumulate.

The local result should remain available in memory so a future optimizer can call it repeatedly without FITS reload,
preprocessing, or process startup. Optional FITS persistence is useful for validation, but writing a product for every
optimizer step must not be required.

## Current implementation findings

### The local regression seam already exists conceptually

`P4Reduction::regions()` currently performs four responsibilities in one method:

1. build annular `P4PixelGrid` geometry and temporal-image selections;
2. assemble one full target vector and predictor matrix for each detector search pixel;
3. call `P4PCA` and retain either compact or full residual state; and
4. materialize, derotate, combine, and write complete output planes.

Each detector search-pixel fit is numerically independent. The worker already owns all of the scratch needed for one
local refit. Pixel-local processing therefore does not require a new PCA algorithm; it requires extracting the
geometry/session setup and the single-search-pixel sampling/fit operation so both complete and sparse runners use one
authoritative implementation.

### Final-stamp locality is sparse in both pixel and frame

For output sky pixel `q` and frame `t`, ordinary cubic derotation reads a four-by-four footprint of integer detector
residual pixels around the inverse-rotated coordinate. For a `W`-by-`W` sky stamp, one frame therefore needs at most
approximately `(W+3)^2` detector residual samples before masks and duplicates are considered.

Those detector footprints rotate through the sequence. The complete work description is not one fixed detector
patch, but a sparse map

```text
detector search pixel p -> target frame indices whose local derotation footprint uses p.
```

The number of requested `(p,t)` samples is bounded by approximately `(W+3)^2*T`. The number of expensive regressions
is the number of unique detector pixels in that map, generally a swept annular band and often much smaller than the
complete P4 search field. Each unique `p` is fit once using all `T` rows, but only its requested residual rows are
retained.

This distinction is required for exactness: calculating an independent fit for each `(p,t)` or shortening its time
series would change P4.

### Existing geometry can support the sparse map

The configured annular grids already provide every search coordinate, exact predictor interpolation record,
common-mask validity, and annulus ownership. `m_ownership` maps an image coordinate to an owning annulus, but a local
runner also needs a checked reverse lookup from `(row,column)` to the annulus-local search index. Building that lookup
alongside ownership avoids repeated linear searches.

The sparse output map must preserve annulus boundaries. A derotation footprint may use detector pixels belonging to
different configured annuli, whose realized integer mode counts can differ even though their output-plane fractions
are shared.

### A full shifted trial cube is unnecessary

The current fake injector materializes a full shifted image before adding it to each cube plane. A local evaluator
needs the perturbation only at target pixels, P4 predictor footprints, and any positive-`numberImages` direct temporal
predictors used by a selected regression. A bounded trial-source sampler can evaluate the centered PSF at those
detector coordinates with the same cubic convention and zero exterior padding.

The trial-source sampler must first be validated against `ADIobservation::injectFake()` on complete small images,
including asymmetric templates, even and odd template dimensions, fractional separation/PA, positive and negative
contrast, angle wrap, and per-frame scale factors. Only after that oracle agrees should the sparse sampler become the
production path.

### Final processing needs a reusable local form

The current compact finalizer rematerializes a full detector cube for one mode, calls the shared full-image
derotation/combination lifecycle, and then discards it. Pixel-local processing should instead:

1. use the sparse map to reconstruct only the requested `W`-by-`W` sky residual frame for each target image;
2. propagate validity through the exact cubic footprint;
3. apply the configured final estimator and weights to those local frame samples; and
4. retain only the combined local mode stamp and validity.

The arithmetic and validity rules should be factored from the shared finalization code rather than copied. Active
post-median subtraction is a detector-pixel temporal operation; it can be exact because one local regression produces
the complete residual time series transiently, but it should remain disabled until the sparse finalizer has a direct
equivalence test.

## Proposed configuration

Use one opt-in P4 flag and the existing fake-source settings for the first command-line surface:

| Proposed target | Contract |
|---|---|
| `p4.localStampSize` | Positive square final sky-stamp width; `0` disables local processing and preserves the complete reduction. |
| `fake.fileName` | Trial PSF template at the post-preprocessing P4 input stage. Initially require `fake.method=single`. |
| `fake.sep`, `fake.PA`, `fake.contrast` | Equal-length independent trial tuples. Negative contrast is the normal negative-companion case. |
| `fake.scaleFileName` | Optional existing per-input-frame scaling, matched after input selection exactly as in ordinary injection. |

Initial validation should require:

- `preProcess.skip=true`;
- detector-frame P4;
- a positive odd `p4.localStampSize` so the returned stamp has one unambiguous integer center pixel;
- a configured final combination method;
- one single-image PSF template whose crop can cover the requested local support; and
- finite trial values, nonnegative separation, and nonzero finite contrast.

The trial source itself may have a fractional sky center. Define the local integer output lattice explicitly from that
continuous center and test the rounding rule at half-pixel boundaries. Do not make the template crop size implicit
from a large full-frame PSF: `p4.localStampSize` is the requested science/result region, while the source sampler must
retain enough template support to reproduce all nonzero template values selected by a separately established crop
rule. If one size cannot express both needs cleanly, add `p4.localTemplateSize` before implementation rather than
overloading the result-stamp width.

QUESTION: Should the first CLI accept multiple trial tuples in one run, or require exactly one tuple while the
reusable C++ session API supports repeated calls? Requiring one initially gives simpler product identity; vector
support is useful for a coarse negative-contrast scan without a new optimizer.

ANSWER: let's start with the simplest form that let's us evaluate performance.

RESOLUTION: the first CLI accepts exactly one trial tuple. The reusable in-memory evaluation seam must still support
repeated calls so performance and state isolation can be measured before an optimizer is added.

QUESTION: Is `preProcess.skip=true` an acceptable initial requirement for negative optimization? This is the only
simple exact contract with the current inherited fake-injection point and is also the efficient workflow for repeated
trials.

ANSWER: yes for now. However I think that pre-processing is generally low-cost compared to the rest of the algorithm so I think we can consider full pre-processing followed by pixel-local-only analysis.

RESOLUTION: initially require `preProcess.skip=true` and a template at the P4-input stage. A subsequent exact path may
start each trial from pristine input, perform ordinary full-image fake injection and preprocessing, and then run only
the pixel-local P4 calculation. Nonlinear preprocessing remains full-cube work for each trial, but its measured cost
is small relative to P4 and it preserves the existing injection semantics.

QUESTION: Should the local result and injected-template crop use one configured width, or should the result use
`p4.localStampSize` and the template have an independent `p4.localTemplateSize`?

ANSWER: I think they can be the same.  I take this to mean that the number of reprocessed pixels is set by the size of the input?

RESOLUTION: use one `p4.localStampSize` initially for both the centered PSF-template crop and the returned sky stamp.
The sparse dependency builder still expands the requested output through cubic derotation support and unions the
resulting detector footprints over field rotation, so this width does not directly equal the number of detector-pixel
regressions. Split template and result sizes later only if PSF wings or the fitting aperture require different support.

## Implementation update (2026-08-21)

The initial detector-frame command-line path is implemented. `P4TrialSource` prepares one phase-preserving centered
template and samples it only where a local fit needs target or predictor flux. `P4LocalGeometry` constructs the exact
sparse inverse-derotation dependencies, and the complete and local runners now share one authoritative
`fitDetectorSearch()` implementation. The local runner fits each unique detector pixel once with the complete temporal
sequence, retains only requested residual rows, materializes only the local sky stack, and uses the existing final
combination implementation.

The center convention is explicit throughout. The image and source-template center is always `0.5*(N-1)`. The
continuous trial coordinate is retained independently of the nearest-integer output-lattice anchor. An odd requested
result therefore does not force an odd source crop: with the supplied 256-by-256 template and 128-by-128 detector, an
11-by-11 result uses a centered 12-by-12 internal template crop and preserves the original half-pixel phase. The sparse
shift and rotation arithmetic also uses the same float specialization and operation ordering as production
`imageShift()` and `imageRotate()` so interpolation anchors cannot differ at a precision boundary.

Focused tests compare sparse source sampling against a fully materialized shifted template; compare sparse positive,
negative, and zero-angle derotation footprints against full rotation; exercise exact half-pixel anchors, masks, odd
and even template centers, both contrast signs, per-frame scales, and positive `p4.numberImages`; and compare every
valid output of a finite-amplitude local refit against the matching full-rerun crop. A synthetic positive source plus
the exact local negative trial also reproduces a separately reduced signal-free control.

A real-data smoke test used `/home/jrmales/Source/mxWork/NACO/AFLep/2011-10-21/out/psf_reg_median.fits` with an 11-pixel
stamp. It completed 377 unique detector fits and retained 53,236 sparse detector/frame samples with a 12-by-12 internal
source crop. This confirms the intended path and resource reporting; it is not yet a local-versus-complete benchmark.

Verification passes the optimized build, all 25 CTest targets, the documentation build with no warnings from the new
or modified P4 files, and focused ASan/UBSan runs of the sparse helpers and complete finite-amplitude equivalence test.
The aggregate hciReduce coverage target was also attempted: 24/25 instrumented tests passed, but the unrelated
`P4PSFReconstructor` test exceeded both the target's 300-second timeout and an isolated 900-second retry (the same test
passes in 8.45 seconds in the optimized suite), so a fresh aggregate LCOV artifact was not produced. The separately
maintained current mxlib LCOV trace records 100% executable-line coverage for every mxlib API called by the edited
functions, including the exact float shift/rotation/kernel ranges.

The first CLI deliberately remains one trial per run and requires `preProcess.skip=true`. Repeated-trial optimizer
control, broader inverse predictor-footprint closure, post-median subtraction, rotated regression, and optional
full-image patch products remain deferred as described below.

## Proposed reusable design

### 1. `P4TrialSource`

An immutable prepared source containing the centered/cropped template, continuous sky position, contrast, and
per-frame scales. It exposes checked sampling of the trial perturbation at detector coordinates for a requested input
frame. Preparation performs all template I/O and coordinate validation outside worker loops.

### 2. `P4LocalRequest`

The requested sky output lattice plus its continuous source center. It builds the deterministic sparse
detector-pixel/frame dependency map using the production derotation geometry and cubic footprint rules.

### 3. `P4ReductionSession`

The pristine post-preprocessing cube, P4 annular grids, reverse coordinate lookup, temporal selections, mode mapping,
mask state, final weights, and memory policy. Construction happens once. A complete runner and a local runner both
call one extracted search-pixel fit routine so their PCA arithmetic cannot drift.

### 4. `P4LocalResult`

Mode-major combined residual and validity stamps, full-image origin, source center, trial provenance, unique-fit and
sparse-sample counts, and timing. FITS writing is a separate optional operation.

These names describe responsibilities rather than mandatory public API names. Keep the first implementation inside
hciReduce unless a second consumer establishes a stable public interface.

## Work sequence

### 1. Lock source-sampling and local-output oracles

- [ ] Add a self-contained asymmetric synthetic PSF and compare sparse trial-source samples with a full
  `ADIobservation::injectFake()` image over integer and fractional source positions, both contrast signs, multiple
  derotation angles, even/odd template dimensions, edge truncation, and per-frame scales.
- [x] Define the continuous source center to integer output-stamp mapping, including exact half-pixel cases.
- [x] On a small synthetic residual cube, enumerate the production cubic derotation footprints for a requested sky
  stamp and prove the sparse `(detector pixel, frame)` reconstruction matches a cropped full derotation.
- [x] Cover mask and image-edge validity before introducing PCA.

### 2. Extract one authoritative detector-pixel fit

- [x] Refactor the existing P4 worker body into a helper that assembles `y` and `X`, applies an optional trial-source
  perturbation, calls `P4PCA`, and returns residuals plus rank status.
- [x] Preserve the disabled/complete reduction's allocation and arithmetic behavior when no perturbation is supplied.
- [x] Support both same-image predictors and positive `p4.numberImages`, including endpoint fallback selections.
- [x] Unit-test the helper against the existing complete runner for unperturbed and finite-amplitude perturbed fits.
- [x] Apply the mxlib LCOV gate to every edited function containing mxlib API calls and record any gap in
  `agents/plans/mxlib_cleanup.md` as required by `AGENTS.md`.

### 3. Build sparse dependency bookkeeping

- [x] Add a checked `(row,column) -> (region,search index)` reverse lookup alongside `m_ownership`.
- [x] Generate and deduplicate all detector residual footprints needed by every local sky pixel and target frame.
- [x] Store the required frame indices for each unique detector search pixel in deterministic order.
- [x] Reject requests outside configured annular ownership or report partial validity according to one explicit rule;
  never reinterpret an uncalculated detector pixel as zero.
- [ ] Test rotations in both directions, angle wrap, duplicate footprints, annulus crossings, masks, and calculation
  halos.

### 4. Implement sparse refit and local finalization

- [x] Refit each unique detector search pixel using the complete temporal sequence and trial-perturbed target and
  predictor samples.
- [x] Retain only residual rows requested by the sparse dependency map; keep full per-fit residuals only in worker
  scratch.
- [x] Reconstruct only the local sky-frame stack with exact cubic validity, then apply the shared combination policy
  and weights.
- [x] Return in-memory mode and validity stamps with geometry, trial, rank, and timing metadata.
- [x] Add the local persistent and per-worker storage to the existing P4 memory-budget report.

### 5. Prove equivalence to a complete finite-amplitude rerun

- [x] With `preProcess.skip=true`, compare every local mode pixel and validity bit against a crop from an ordinary
  full-field reduction after full-image trial injection.
- [ ] Cover zero perturbation, positive injection, negative injection, fractional position, multiple annuli and mode
  mappings, positive `p4.numberImages`, weights, and each supported final combination estimator.
- [x] Construct a synthetic cube containing a known source, inject the exact negative source locally, and verify the
  local result matches a separately generated signal-free control.
- [ ] Run ASan/UBSan, focused Catch2 tests, complete CTest, and coverage.

### 6. Add the opt-in command-line path

- [x] Register and validate `p4.localStampSize`; preserve all current behavior and products at zero.
- [x] Prevent inherited pre-preprocessing fake injection in local mode by requiring the established
  `preProcess.skip=true` contract, then prepare the existing fake configuration as local trials.
- [x] Optionally write local residual and validity products with deterministic trial indices and complete FITS
  provenance, while keeping in-memory evaluation independent of output.
- [x] Document configuration, output geometry, resource scaling, preprocessing-stage semantics, and a reproducible
  one-trial comparison command.

### 7. Measure and connect a negative-companion optimizer

- [ ] On the small local configuration, record unique detector fits, requested sparse samples, wall time, worker time,
  and peak RSS relative to a complete rerun.
- [ ] Repeat on the remote 1-second-coadd dataset and confirm repeated trials reuse loaded/preprocessed state.
- [ ] Define the optimizer's merit function, fitting aperture, selected P4 mode or joint-mode rule, parameter bounds,
  convergence rule, and uncertainty method in a separate plan.
- [ ] Use the fitted negative source in one complete reduction, regenerate the spatially variable P4 PSF/filter, and
  determine whether the response-normalization trough and filtered ring remain in the signal-free state.

## Deferred broader influence products

A trial source can change the regression at an output detector pixel outside the requested sky stamp when the source
overlaps one of that pixel's OR or additional-image predictor samples. Finding every such pixel is an inverse
predictor-footprint problem. It matters if a local trial is used to patch a complete final image, but not for an exact
merit function evaluated only inside the requested local stamp: every regression consumed by that stamp is already
refit with the source present.

Defer the broader closure initially. If added, expose it explicitly and default it off as proposed in
`developer_todo.md`. A full-image replacement product must include all affected regressions or state clearly that it
replaces only the requested output region; it must never imply global equivalence from a local patch.

## Acceptance criteria

- With local processing disabled, complete P4 arithmetic, products, memory selection, and timing remain unchanged.
- A local finite-amplitude result agrees with the same crop from a complete post-preprocessing injection/rerun at the
  established float tolerance, including validity.
- Every local regression uses the complete configured temporal sample set and the ordinary P4 eigensolve.
- Trial flux affects target, same-image OR, and additional-image predictor samples exactly where their production
  interpolation footprints overlap it.
- Repeated trials start from the same pristine in-memory cube and do not accumulate injected flux or reduction state.
- Memory scales with the prepared input cube, unique local fits, sparse residual requests, and worker scratch—not with
  a new full shifted trial cube or full residual output cube per trial.
- Output and in-memory results identify the full-image stamp origin, continuous trial position, contrast, template,
  P4 modes, injection stage, support, and validity.
- An exact synthetic negative source reproduces a separately constructed signal-free control.
