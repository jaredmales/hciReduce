# Consolidated hciReduce Issues

## Purpose

This file is the central index of unresolved work described across `agents/plans`. The linked planning documents
remain the source of detailed mathematics, contracts, measurements, and implementation sequences; this index records
which work is still live and prevents the same issue from being rediscovered in several plans.

The inventory was reconciled on 2026-08-31. Newer implementation/status prose takes precedence over older unchecked
checklists. An unchecked source checkbox is therefore not automatically an open issue, and completed historical plans
are not copied here merely because their original work sequence was not edited in place.

Use these states:

- **Open**: concrete work can proceed under an existing contract.
- **Decision required**: a scientific, API, or product contract must be selected first.
- **Deferred**: intentionally postponed until its stated prerequisite or priority changes.
- **Conditional**: do the work only if measurement demonstrates the stated need.
- **External**: owned by mxlib or another dependency, but relevant to hciReduce.
- **Untriaged**: the source note does not yet contain enough context or acceptance criteria.

When resolving an issue, update this index and the linked source plan together. Keep issue IDs stable; mark an entry
resolved rather than reusing its ID.

## Current issue index

| ID | State | Issue |
|---|---|---|
| HCI-001 | Open | Implement bounded KLIP memory management |
| HCI-002 | Open | Complete and archive the P4 golden scientific comparison |
| HCI-003 | Decision required | Select the P4 final-derotation boundary and interpolation policy |
| HCI-004 | Decision required | Define a lossless P4 postprocess checkpoint and input-header provenance |
| HCI-005 | Open | Run the fitted-negative complete-field signal-free validation |
| HCI-006 | Open | Complete spatially variable P4 PSF/filter validation |
| HCI-007 | Open | Calibrate negative-companion jackknife uncertainty |
| HCI-008 | Open | Capture the full remote pixel-local/optimizer reuse benchmark |
| HCI-009 | Deferred | Compute broader pixel-local influence for full-image patch products |
| HCI-010 | Open | Close the exact factor-deletion promotion gates |
| HCI-011 | Deferred | Establish a bounded-rank projected exclusion backend |
| HCI-012 | Deferred | Extend structured SVD deletion beyond one removed row |
| HCI-013 | Deferred | Stream P4 per-frame residual products one mode at a time |
| HCI-014 | Open / Conditional | Calibrate P4 memory policy; add out-of-core input only if needed |
| HCI-015 | Conditional | Evaluate a broader float-precision P4 path |
| HCI-016 | Conditional | Add bounded P4 interpolation-kernel caching |
| HCI-017 | Decision required | Decide whether rotated P4 prunes mask-affected predictor columns |
| HCI-018 | Open | Complete rotated-P4 scientific validation and release gates |
| HCI-019 | Deferred | Investigate the original globally coupled rotated coefficient field |
| HCI-020 | Deferred | Research improved temporal predictors |
| HCI-021 | Deferred | Adapt factor deletion to KLIP image exclusion |
| HCI-022 | Deferred | Evaluate alternate P4 solvers and Tikhonov regularization |
| HCI-023 | Deferred | Evaluate grouping neighboring P4 target pixels |
| HCI-024 | Decision required | Inherit preprocessing provenance when preprocessing is skipped |
| HCI-025 | Deferred | Evaluate covariance-weighted and locally estimated matched filters |
| HCI-026 | Deferred | Build the external multi-algorithm grid runner |
| HCI-027 | Deferred | Define a P4 RDI or alternate-predictor extension |
| HCI-028 | Deferred | Extend the P4 forward model beyond the current single-template contract |
| HCI-029 | External / Deferred | Repair and reassess the optional CUDA eigensolver path |
| HCI-030 | Open | Reject malformed saved-reduction header vectors deterministically |
| HCI-031 | External / Decision required | Correct `appConfigurator::get` conversion-error reporting |
| HCI-032 | Deferred | Close residual platform-fault and Doxygen-warning debt |

## Core reduction and product issues

### HCI-001 — Implement bounded KLIP memory management

**State:** Open.

Implement the planned compact residual store and one-mode finalization, sparse reference-selection bookkeeping,
resource-aware worker selection, and selected-covariance/copy reductions. Preserve the existing residual, validity,
diagnostic, FITS, and ordering contracts. Streaming requested residual products and changing input retention follow
only after the bounded core is measured.

Sources: [developer TODO](developer_todo.md#klip-memory-management),
[KLIP memory work sequence](klip_memory_management_optimization.md#work-sequence), and
[KLIP acceptance criteria](klip_memory_management_optimization.md#acceptance-criteria).

### HCI-002 — Complete and archive the P4 golden scientific comparison

**State:** Open.

Archive the expected AF Lep/NACO geometry, residual, combined-image, and provenance products outside the unit-test
repository. Compare coordinates, valid pixels, predictor counts, realized modes, residuals, and combined products
against the prototype; explain intentional differences and set platform-tolerant rules for nearly degenerate mode
cutoffs.

Sources: [P4 status](P4.md#status), [Phase 0 golden-reference gate](P4.md#phase-0-dependency-and-contract-gate), and
[Phase 6 scientific parity](P4.md#phase-6-scientific-parity-and-production-hardening).

### HCI-003 — Select the P4 final-derotation boundary and interpolation policy

**State:** Decision required.

Separate calculation/support annuli from retained science-output bounds, then choose and validate the boundary
sampler. Complete-support cubic is the conservative baseline; bilinear, adaptive, masked positive-weight, or
nearest-valid behavior requires explicit provenance and photometry, astrometry, noise, angle, and deterministic
validity tests. Do not let a calculation halo silently change the configured science annulus or predictor set.

Source: [P4 Phase 6 boundary-policy work](P4.md#phase-6-scientific-parity-and-production-hardening).

### HCI-004 — Define a lossless P4 postprocess checkpoint and input-header provenance

**State:** Decision required.

Specify a versioned raw-residual checkpoint with one authoritative validity cube before enabling P4 `postprocess`.
The checkpoint must prevent final stages from being applied twice and invalid samples from being consumed. In the
same product-contract review, decide how arbitrary AF Lep input-header provenance is carried into the final combined
product.

Sources: [P4 shared final lifecycle](P4.md#phase-3-shared-adi-final-lifecycle) and
[P4 application phase](P4.md#phase-5-p4reduce-application).

### HCI-005 — Run the fitted-negative complete-field signal-free validation

**State:** Open.

Add the one-shot complete-field post-preprocessing injection path, apply the fitted negative source, and run one
ordinary full reduction. Regenerate the spatially variable PSF and filtered product, then determine whether the
response-normalization trough and filtered ring disappear while unrelated unfiltered science remains stable.

Sources: [developer TODO](developer_todo.md#p4-psf-calculation-and-post-processing),
[optimizer final validation](p4_negative_companion_optimization.md#final-signal-free-validation), and
[pixel-local measurement step](p4_pixel_local_processing.md#7-measure-and-connect-a-negative-companion-optimizer).

### HCI-006 — Complete spatially variable P4 PSF/filter validation

**State:** Open.

Validate response localization, detection ranking, support at annulus/image boundaries, isolated sources at several
contrasts, neighboring sources, and noise-only data. Add application-level real-FITS tests for opt-in behavior,
directory creation, actionable errors, and unchanged disabled outputs. Record a representative remote run including
positive `p4.numberImages`, peak RSS, product size, and comparison with complete fake injection/refitting.

Sources: [developer TODO](developer_todo.md#p4-psf-calculation-and-post-processing),
[PSF/filter work sequence](p4_psf_calculation_post_processing.md#8-add-spatially-variable-psf-filtering), and
[PSF/filter acceptance criteria](p4_psf_calculation_post_processing.md#acceptance-criteria).

### HCI-007 — Calibrate negative-companion jackknife uncertainty

**State:** Open.

Add an interleaved-block alternative to the implemented contiguous delete-one-block jackknife. Use negative-source
injection/recovery at nearby empty locations and comparable separations to compare block counts and contiguous versus
interleaved covariance before treating either result as calibrated science uncertainty.

Sources: [developer TODO](developer_todo.md#p4-pixel-local-processing) and
[optimizer uncertainty contract](p4_negative_companion_optimization.md#uncertainty-method).

### HCI-008 — Capture the full remote pixel-local/optimizer reuse benchmark

**State:** Open.

Benchmark repeated evaluations in one process on the remote 621-frame, one-second-coadd dataset. Capture console wall
time, phase timing, peak RSS, and evidence that loading/preprocessing state is reused. The existing full-size local
trial recorded geometry but not these integration measurements.

Sources: [optimizer initial measurement](p4_negative_companion_optimization.md#initial-performance-measurement),
[optimizer work sequence](p4_negative_companion_optimization.md#work-sequence), and
[pixel-local measurement step](p4_pixel_local_processing.md#7-measure-and-connect-a-negative-companion-optimizer).

### HCI-009 — Compute broader pixel-local influence for full-image patch products

**State:** Deferred.

If local results are ever used to patch a complete output image, compute the inverse predictor-footprint closure so
every regression changed by the trial is included. Keep this off by default. The present local merit contract is
already exact inside its requested stamp and must not be described as globally equivalent.

Source: [deferred broader influence products](p4_pixel_local_processing.md#deferred-broader-influence-products).

## Exclusion solvers and performance issues

### HCI-010 — Close the exact factor-deletion promotion gates

**State:** Open.

Complete the remaining one-versus-many-worker comparison with inner BLAS threading controlled, forced
application-level rank-boundary fallback, forced memory-limit behavior, and real-data science comparison. Capture any
still-missing production-scale evidence, and then decide whether `factorDowndateExact` or `rankOneSecular` should
remain opt-in or become eligible for automatic/default selection. A oneMKL compile/link smoke test remains
conditional on oneMKL availability.

Sources: [current P4 downdate checkpoint](p4_svd_downdate.md#current-p4-integration-and-review-checkpoint),
[ROC structured-backend result](p4_svd_downdate.md#roc-structured-backend-result),
[promotion step](p4_svd_downdate.md#9-promote-only-after-science-and-operations-acceptance).

### HCI-011 — Establish a bounded-rank projected exclusion backend

**State:** Deferred.

Sweep guard rank, quantify the approximation caused by selecting the global truncated subspace with the target row
present, and define residual, leakage, injection-recovery, and science-metric tolerances. Keep projected provenance
distinct from exact exclusion and do not expose a public guard-rank control until convergence data requires it.

Source: [bounded working-rank step](p4_svd_downdate.md#6-add-bounded-working-rank-and-convergence-checks).

### HCI-012 — Extend structured SVD deletion beyond one removed row

**State:** Deferred.

Implement and compare a stable small-block deletion method with within-target rank-one steps, always starting from
the immutable base factor. Retain a stable dense/complement oracle and explicit fallback. This is needed for wider
temporal exclusion sets; the current `rankOneSecular` backend intentionally supports exactly one removed row.

Sources: [structured backend step](p4_svd_downdate.md#8-add-a-structured-deletion-backend-if-dense-solves-remain-limiting)
and [target-exclusion adaptation](p4_target_image_exclusion.md#5-adapt-the-long-males-svd-downdate).

### HCI-013 — Stream P4 per-frame residual products one mode at a time

**State:** Deferred.

Stream the existing `outputPSFSub` filenames and `REDUCTION` indices without retaining the legacy full residual cube
path. Preserve filenames, headers, sidecars, ordering, and recoverable partial-write behavior.

Source: [P4 compact finalization step](p4_memory_management_optimization.md#4-store-residuals-compactly-and-finalize-one-mode-at-a-time).

### HCI-014 — Calibrate P4 memory policy; add out-of-core input only if needed

**State:** Open / Conditional.

Measure remote peak RSS at one and several workers before finalizing the automatic-memory fraction and safety margin.
Only if those measurements show the loaded preprocessed cube remains material after compact output, prototype bounded
annulus/image-plane loading. Account for repeated I/O and overlapping annuli, and preserve all metadata, fake-source,
and preprocessing contracts. Keep general out-of-core preprocessing separate.

Sources: [P4 input-retention step](p4_memory_management_optimization.md#6-reconsider-input-cube-retention-only-if-still-necessary)
and [P4 deferred memory decisions](p4_memory_management_optimization.md#deferred-decisions).

### HCI-015 — Evaluate a broader float-precision P4 path

**State:** Conditional.

Measure whether moving more non-decomposition P4 arithmetic from double to float provides a useful performance or
memory benefit. Any change requires numerical, validity, and science-product equivalence bounds; the current
all-double PCA contract remains authoritative unless separately changed.

Source: [developer TODO](developer_todo.md#evaluate-extend-of-double-vs-float-in-p4).

### HCI-016 — Add bounded P4 interpolation-kernel caching

**State:** Conditional.

Consider on-demand generation or bounded caching only if real configurations show the precomputed
`Nsearch * K * width^2` footprint is material. Preserve exact kernel values and deterministic geometry.

Source: [P4 deferred performance work](P4.md#phase-7-deferred-performance-and-algorithm-extensions).

## Rotated-P4 issues

### HCI-017 — Decide whether rotated P4 prunes mask-affected predictor columns

**State:** Decision required.

The reopened proposal would retain a target only when it has complete all-frame support while pruning an entire
predictor column if any mapped sample is invalid. Before implementation, settle the resulting per-target predictor
count, structural-rank, annulus-consistency, provenance, and empty-predictor contracts.

Source: [rotated P4 Phase 0](P4_rotated.md#phase-0-freeze-the-scientific-definition).

### HCI-018 — Complete rotated-P4 scientific validation and release gates

**State:** Open.

Finish the direct/materialized/shift-then-rotate reference checks, finite-difference objective check, baseline-treatment
matrix, remaining geometry cases, and stage-specific progress/timing. Run injected companions across separation,
position angle, contrast, rotation, and modes; measure throughput, astrometric bias, noise, pedestal, and support under
each baseline treatment. Complete coverage, ASan/UBSan, OpenMP determinism, AF Lep comparison, and controlled
performance benchmarks before calling the variant scientifically ready.

Sources: [rotated P4 implementation phases](P4_rotated.md#implementation-phases),
[verification/performance gate](P4_rotated.md#phase-4-verification-and-performance-gate), and
[release criteria](P4_rotated.md#phase-5-documentation-and-release-criteria).

### HCI-019 — Investigate the original globally coupled rotated coefficient field

**State:** Deferred.

This is a separate research program, not completion of the current per-sky-pixel rotated mode. It first needs the
local-versus-global coefficient, baseline/output, regularization, mask/support, and corrected AF Lep angle contracts;
then tiny explicit oracles, local temporal-dual and global matrix-free prototypes, solver comparisons, scientific
sweeps, and production-selection gates.

Source: [full rotated-frame status and scope](P4_rotated_full.md#status-and-scope) and
[staged investigation plan](P4_rotated_full.md#staged-investigation-plan).

## Algorithm and workflow extensions

### HCI-020 — Research improved temporal predictors

**State:** Deferred.

Compare smaller additional-image patches, mean/median patch summaries, endpoint-specific coefficient sets, PCAT-like
temporal bases, and irregular-cadence linear prediction using an autocorrelation or PSD model. Each alternative needs
an explicit missing-time, interpolation/extrapolation, validity, and provenance contract.

Source: [developer TODO](developer_todo.md#better-temporal-prediction).

### HCI-021 — Adapt factor deletion to KLIP image exclusion

**State:** Deferred.

Define KLIP's image-exclusion semantics and adapt the reusable factor-deletion machinery without assuming P4's
predictor/response layout. Preserve KLIP mode, reference-selection, and output contracts.

Source: [developer TODO](developer_todo.md#exclusion-of-target-images).

### HCI-022 — Evaluate alternate P4 solvers and Tikhonov regularization

**State:** Deferred.

Compare direct SVD and SVD plus Tikhonov regularization with the current Gram eigensolver under explicit rank,
regularization-selection, residual, and provenance contracts.

Source: [developer TODO](developer_todo.md#other-solvers).

### HCI-023 — Evaluate grouping neighboring P4 target pixels

**State:** Deferred.

Test whether fitting a shared solution for a small block such as 2-by-2 pixels improves turnaround or fit quality
without introducing unacceptable spatial bias. Define how the shared solution maps back to individual pixels.

Source: [developer TODO](developer_todo.md#pixel-grouping).

### HCI-024 — Inherit preprocessing provenance when preprocessing is skipped

**State:** Decision required.

Define an opt-in inheritance contract for already-preprocessed inputs, validate that all input files agree, and write
the inherited preprocessing settings to every downstream FITS product. Decide mismatch and missing-card behavior
before implementation.

Source: [developer TODO](developer_todo.md#propagation-of-pre-processing).

### HCI-025 — Evaluate covariance-weighted and locally estimated matched filters

**State:** Deferred.

Evaluate covariance weighting from the eigenvector projection and a signal-free response estimated from nearby PSFs
outside an exclusion radius at the same separation. The latter likely needs a rotation-aware rather than shift-only
operator.

Source: [developer TODO](developer_todo.md#matched-filtering-updates).

### HCI-026 — Build the external multi-algorithm grid runner

**State:** Deferred.

Settle the nine science/coordinate/grid/output/restart/parallelism decisions, then implement a standalone runner with
a pure grid generator, common trial schema, isolated KLIP/P4 process adapters, deterministic products, a flushed
manifest, restart/error policy, scientific equivalence tests, documentation, and resource measurements. Add
process-unique temporary diagnostic names before allowing concurrent writers to share a directory.

Sources: [future runner status](grid_mode.md#status), [implementation phases](grid_mode.md#future-implementation-phases),
and [P4 distributed-writer prerequisite](P4.md#phase-7-deferred-performance-and-algorithm-extensions).

### HCI-027 — Define a P4 RDI or alternate-predictor extension

**State:** Deferred.

P4 currently remains target-only. Any RDI or alternate predictor must enter behind the pure PCA/predictor seam and
define training rows, scaling, geometry, rank, invalidity, output meaning, and provenance as a separately named
algorithm contract.

Sources: [P4 RDI decision](P4.md#7-rdi-behavior) and
[P4 deferred extensions](P4.md#phase-7-deferred-performance-and-algorithm-extensions).

### HCI-028 — Extend the P4 forward model beyond the current single-template contract

**State:** Deferred.

Open extensions include a post-median-subtraction filtering approximation, inherited preprocessing as a cube-level
forward operator, per-image/list or wavelength-dependent templates, per-frame flux scales, and optional disk-backed
model blocks only if measured retained state becomes material. Keep replacement of ordinary `finim` explicitly
opt-in.

Source: [PSF deferred decisions](p4_psf_calculation_post_processing.md#deferred-decisions).

## Dependency and maintenance issues

### HCI-029 — Repair and reassess the optional CUDA eigensolver path

**State:** External / Deferred.

Before any hciReduce CUDA caller, mxlib must provide safe noncopyable RAII handle/parameter ownership, a reusable
all-double `syevdx` wrapper, behavioral device tests, a CUDA-enabled coverage lane, and machine-readable installed
CUDA capability. Current measurements show that direct repeated GPU eigensolves are not a useful production path;
reconsider CUDA only for a materially different workload such as batched structured downdates, and require measured
end-to-end benefit before integration.

Sources: [mxlib CUDA prerequisite](mxlib_cleanup.md#additional-audit-items),
[P4 CUDA proposal](P4.md#phase-7-deferred-performance-and-algorithm-extensions), and
[target-exclusion recommendation](p4_target_image_exclusion.md#current-status-and-recommendation).

### HCI-030 — Reject malformed saved-reduction header vectors deterministically

**State:** Open.

`parseStringVector` returns no status and appends a converted final token, so current size-only checks cannot reliably
detect empty or malformed saved KLIP/P4 header vectors. Validate raw header text/tokens in hciReduce or design a
separate status-returning mxlib API, then add malformed-header regression tests.

Source: [mxlib cleanup HCI follow-up](mxlib_cleanup.md#additional-audit-items).

### HCI-031 — Correct `appConfigurator::get` conversion-error reporting

**State:** External / Decision required.

The documented API returns `-1` on conversion errors, but the implementation uses an exception-free conversion
without requesting its error code and can return success with a fallback value. Decide whether malformed input
preserves the destination or stores the fallback, then repair and cover the mxlib API before relying on the contract.

Source: [mxlib cleanup API follow-up](mxlib_cleanup.md#additional-audit-items).

### HCI-032 — Close residual platform-fault and Doxygen-warning debt

**State:** Deferred.

The remaining `HCIobservation::readPSFSub` coverage miss requires deterministic filesystem fault injection for an
error while incrementing a directory iterator after successful construction. The documentation build also retains
pre-existing warnings; classify and eliminate them without weakening warning policy. These do not currently block
the tested reduction paths.

Sources: [initial test Phase 11](initial_tests.md#phase-11-hciobservation-branch-closure-complete-except-platform-fault-injection)
and [mxlib documentation debt record](mxlib_cleanup.md#phase-5-repair-documentation-infrastructure-and-debt).

## Reconciled as complete or intentionally excluded

The following source material was reviewed but does not create a current issue:

- `documentation_update.md`, `hciAnalyze.md`, `fake_planet_and_config_cleanup.md`, `snr_small_sample.md`, and
  `p4_profiling.md` are superseded by their completion records and current implementation evidence.
- `Adding_More_Images.md` is implemented by the positive-`p4.numberImages` detector-frame path and its current tests.
- The unchecked test/verification bullets in `p4_pixel_local_processing.md` are superseded by its 2026-08-21
  implementation update; only the remote benchmark, final signal-free run, and optional broader influence closure
  remain live above.
- The old `vectorVariance` RMS follow-up in `mxlib_cleanup.md` is superseded by the local true-RMS implementation
  recorded in `initial_tests.md`; mxlib's sample-variance API need not be redefined.
- The items under `mxlib_cleanup.md`'s **Not related to hciReduce** section belong in mxlib's own issue tracker and are
  not duplicated here.
