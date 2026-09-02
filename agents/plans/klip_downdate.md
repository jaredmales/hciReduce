# KLIP Factor-Deletion Integration Plan

## Goal

Replace KLIP's repeated target-specific reference-library eigensolves for ordinary ADI image exclusion with one
immutable base factorization per search region and a target-specific column-deletion update. The first production
path must reproduce the current direct selected-library KLIP result within explicit numerical tolerances; it is an
optimization, not a change in subtraction semantics.

The reusable mxlib operation already exists in `mx::math::svdRemoveColumns()`. This plan adapts that operation to
KLIP's image-library orientation and preserves the direct solve as the correctness oracle and fallback. The generic
linear-algebra contract is documented in [p4_svd_downdate.md](p4_svd_downdate.md); this plan does not duplicate or
change mxlib's ownership.

## Current KLIP behavior

For one region, let `X` be the mean-subtracted, optionally normalized reference matrix with `P` region pixels in rows
and `R` reference images in columns. For each target `t`, `collapseCovar()`:

1. applies `adi.excludeMethod` / `adi.excludeMethodMax`;
2. optionally applies `klip.includeMethod` / `klip.includeRefNum` ranking;
3. gathers the selected columns `X_S(t)` and their covariance; and
4. calls `calcKLModesAdaptive()` on that gathered library.

The result is exact for the current selected columns but repeats a covariance eigensolve for every target whenever
selection or exclusion is active. Existing centering and pixel time-series normalization occur before selection, so
they are properties of the complete reference library and must remain so in the factor-deletion path.

## Mathematical contract

Use a complete numerical-rank thin factorization of the already processed base library,

$$
X = U\Sigma V^T.
$$

For the target-specific excluded reference-image indices `E`, `mx::math::svdRemoveColumns()` receives the base
singular values and the corresponding rows `V_E`. It returns updated singular values and a preserved-left rotation
`W_E`. The retained KL directions are

$$
U_{\mathrm{ret}} = U W_E.
$$

The KLIP worker will store the requested leading rows of `(U_{\mathrm{ret}})^T` in the current ascending-eigenvalue row ordering,
so existing coefficient calculation, mode counting, right-reason projection, and output cubes remain unchanged.

This is exact for the represented matrix only when the base factorization contains its complete numerical rank. A
bounded working-rank or literal reduced Long--Males core is a separate approximate method and is out of scope for the
first implementation.

## Decisions and initial scope

- Support only same-cube ADI (`rims == tims`) with a real image-exclusion criterion. Independent RDI already disables
  ADI exclusion in region orchestration and retains its current master/direct paths.
- Apply factor deletion only when `includeMethod=all` or `includeRefNum=0`; limited correlation/time/angle/image-number
  selection stays on the direct gather-and-solve path. A small retained library is often cheaper to solve directly
  than to delete most of a large base library.
- Keep all current exclusion semantics, including inclusive `adi.minDPx` boundaries, endpoint truncation, and the
  current empty-selected-library error.
- Build an immutable FP64 base factorization after `meanSubtract()`. Each target is derived directly from that base;
  deletion updates are never chained.
- Use `leadingCovariance` for arbitrary multi-image deletion initially, because current KLIP itself is a covariance
  eigensolve. Use `rankOneSecular` only when exactly one image is deleted and only after it matches the dense backend
  and direct oracle. A deletion backend failure, factor-validation failure, or rank-boundary ambiguity falls back to
  the existing direct selected-library calculation.
- Do not add a user-facing solver switch in the first patch. Emit method/backend/fallback provenance and timing so the
  performance and numerical promotion decision is evidence-based.

## Factor construction and target update

1. Preserve the existing `meanSubtract()` call, including mask behavior, reference norms, and target preprocessing.
2. Construct a complete FP64 thin SVD of the processed `P x R` base library. Retain both `U` and `V`, singular values
   in descending order, and the effective numerical rank. Validate the deletion-side `V` factor once with
   `validateSvdDeletionFactor()`.
3. Refactor reference selection so one helper can produce both the current inclusion flags and sorted unique removed
   indices. The direct path gathers the retained complement; the factor path passes removed indices to the deletion
   API. `m_imsIncluded` remains the diagnostic/public record of the final selection.
4. Allocate one `svdDeletionResult<double>` and `svdDeletionWorkspace<double>` per OpenMP worker. Prepare capacities
   outside the target loop using the base rank, requested maximum mode count, and the maximum observed deletion count.
5. For each target, obtain its removed index set. If the retained library is empty, retain the current exception. If
   the target is eligible, invoke column deletion, apply the returned rotation to the base `U`, and assemble rows in
   KLIP's established ordering. Otherwise or on a recoverable numerical outcome, call the existing direct path.
6. Continue using the current coefficient/projection and right-reason code with the resulting mode matrix. The
   factor-deletion branch must not special-case science output values.

Rank-deficient input requires particular care. A zero or non-finite singular direction cannot be reconstructed by
division when forming factors. The base-factor routine must either produce a valid complete positive-rank thin system
and prove factor orthogonality, or select the direct fallback. It must not silently discard finite positive directions
or publish a truncated system as exact.

## Direct oracle and fallback policy

The existing `collapseCovar()` plus `calcKLModesAdaptive()` implementation remains available throughout development.
It is the oracle for tests and the production fallback for:

- RDI or limited-reference inclusion;
- empty selected libraries;
- base-factor validation failure;
- a non-success mxlib deletion status, including an unsupported structured deletion count;
- rank/mode ambiguity where the update and direct solve choose different valid bases near a truncation boundary; and
- any non-finite updated singular system or mode matrix.

Fallback must be per target, must preserve `m_imsIncluded`, and must be counted in timing/provenance. It must not
alter selection or cause other targets to use a chained or rebuilt base factorization.

## Verification matrix

Add focused `KLIPreduction_test.cpp` coverage before enabling the production branch.

- Compare factor-deletion and direct residuals/projectors for wide (`P < R`) and tall/equal (`R <= P`) libraries.
- Cover no deletion, delete-one, contiguous delete-block, noncontiguous angle-like deletion, and beginning/end targets.
- Cover exact-rank-deficient, repeated, zero, and near-rank-boundary matrices; verify direct fallback where a complete
  valid factor cannot be established.
- Verify every requested mode count: supported, greater than retained rank, duplicate requests, and zero/negative
  maximum-mode configuration.
- Verify all exclusion methods and their inclusive boundaries through the common selection helper.
- Verify `includeRefNum` and every non-`all` inclusion method remain on the direct path with byte-equivalent selection
  flags and residuals.
- Verify RDI remains unchanged.
- Compare serial and multi-worker results, including mixed deletion counts and one target forced through fallback.
- Verify right-reason projection, diagnostics, FITS provenance, and all timing fields remain finite and correctly
  attributed.

Use projector/residual comparisons rather than signed individual mode-vector comparisons, because eigensolver signs
and bases inside a degenerate eigenspace are not unique.

## Performance and promotion gates

The optimization is not assumed to be faster in every regime. Measure wall time, aggregate worker time, base-factor
time, deletion time, direct-fallback time, peak memory, effective outer OpenMP thread count, and inner BLAS thread
count for representative regions and real data.

Expected behavior differs by exclusion size:

- Delete-one can use the quadratic `rankOneSecular` update and is the principal expected speedup case.
- Multi-image angular/pixel exclusion initially uses a dense small deletion core, so it may not beat a direct solve;
  benchmark it rather than assuming a gain.
- When only a few references survive due to `includeRefNum`, direct solving is expected to remain preferable.

Before making factor deletion automatic beyond the guarded initial scope, require:

1. agreement with the direct oracle on the full verification matrix and representative injected-companion recovery;
2. no unexplained fallback rate or rank-boundary discrepancy;
3. serial/multi-worker reproducibility with controlled BLAS threading; and
4. a demonstrated memory and wall-time benefit for at least the delete-one production regime.

Document the selected method, backend, fallback count, and measured timing in science FITS provenance once the path is
enabled. A future bounded-rank projected backend requires separate algorithm naming, approximation provenance, and
science acceptance tests.

## Benchmark record — 2026-09-02

The NACO AFLep configuration at
`/home/jrmales/Source/mxWork/NACO/AFLep/2011-10-21/kr/kr.conf` is an appropriate real-data delete-one reference:
it uses same-cube ADI image-number exclusion (`adi.excludeMethod=imno`, `adi.minDPx=0`) on 621 preprocessed frames.
For a bounded, repeatable measurement, the run used the first 128 frames (`--input.deleteBack 493`) and 64 KL modes.
OpenBLAS was held at one thread; outer OpenMP was measured at one and 20 threads.

| Solver | Outer workers | Wall time | Peak RSS | Factor deletion worker time |
| --- | ---: | ---: | ---: | ---: |
| Direct selected-library baseline | 1 | 1.48 s | 272 MiB | n/a |
| `leadingCovariance` factor deletion | 1 | 3.34 s | 300 MiB | 2.08 CPU s |
| Direct selected-library baseline | 20 | 0.27 s | 421 MiB | n/a |
| `leadingCovariance` factor deletion | 20 | 0.97 s | 673 MiB | 4.51 aggregate CPU s |
| `rankOneSecular` factor deletion | 1 | 2.84 s | not recorded | 1.79 CPU s |
| `rankOneSecular` factor deletion | 20 | 0.88 s | not recorded | 4.23 aggregate CPU s |

The dense factor-deletion results agree scientifically with the direct baseline: among finite final-image pixels, the
maximum absolute residual difference was `3.48e-5`, RMS difference was `4.57e-6`, and the maximum difference was
`1.58e-4` of the direct-image RMS. Serial and 20-worker outputs were reproducible to `rtol=1e-5`, `atol=1e-6`.

The generic `leadingCovariance` backend therefore fails this plan's performance and memory promotion gate: it was
2.3x slower with one worker and 3.6x slower with 20 workers. The delete-one `rankOneSecular` backend also fails the
wall-time gate: it was 2.0x slower with one worker and 3.2x slower with 20 workers. Its final image matched the
direct baseline within the same `3.48e-5` maximum absolute difference observed for the dense backend.

An early temporary dual-backend integration run ended with `*** stack smashing detected ***` after writing its output.
That failure is not reproducible in the isolated `rankOneSecular` path: the exact 128-frame/64-mode NACO run now
completes in Release, and a 384-by-128 synthetic 64-mode rank-one test is clean under ASan/UBSan. Keep the backend
unpromoted for performance, but do not attribute the earlier canary failure to mxlib without a reproducer.

### Delete-one crossover map

The same NACO frames were then measured with one outer worker and one BLAS thread. The table reports rank-one
algorithm time divided by direct algorithm time; values below one favor rank-one. Radii are a convenient proxy for
the geometric region size: 10, 20, 40, and 60 pixels contain approximately 208, 1,148, 4,916, and 11,180 annular
pixels before the input mask. Requested modes were half the image count through 128 frames and 64 thereafter.

| Images | Modes | ~208 px | ~1,148 px | ~4,916 px | ~11,180 px |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 16 | 8 | 0.79x | 1.76x | 3.44x | 4.28x |
| 32 | 16 | 0.56x | 1.45x | 3.09x | 4.73x |
| 64 | 32 | 0.41x | 1.00x | 2.34x | 3.40x |
| 128 | 64 | 0.32x | 0.78x | 2.25x | 3.04x |
| 256 | 64 | 0.30x | 0.83x | 2.02x | not measured |
| 512 | 64 | 0.29x | not measured | not measured | not measured |

With 20 outer workers and one BLAS thread, overhead shifts the boundary toward even smaller regions. At 128 images
and 64 modes, rank-one/direct ratios were 1.11x, 1.95x, 4.48x, and 5.56x for the four regions above. At ~208 pixels,
rank-one becomes beneficial again at larger libraries: 0.81x for 256 images and 0.65x for 512 images. At ~1,148
pixels and 256 images it was still 2.57x slower.

The practical conclusion is not a universal image-count threshold. A future automatic delete-one selector must use
at least base numerical rank, selected output-mode count, region pixel count, and outer worker count. The present
implementation should remain direct by default; rank-one is a credible targeted optimization for compact regions
with a large image library.

## Implementation sequence

1. [x] Extract a selection result from `collapseCovar()` that records selected and deleted indices without changing the
       direct algorithm.
2. [x] Add an FP64 base-SVD/factor-deletion helper and establish direct-oracle equivalence.
3. [x] Add worker-private deletion storage, per-target fallback, timing, and diagnostics behind the narrow ADI exclusion
       eligibility condition.
4. [x] Run the focused sanitizer suite and controlled serial/multi-worker benchmarks. The sanitizer suite was clean;
   both dense and rank-one factor deletion are numerically reproducible but fail the performance gate. An earlier
   stack-smashing report from temporary dual-backend integration is not reproducible in the isolated rank-one path.
5. [x] Do not promote factor deletion to the default eligible delete-one solver. Keep direct solving as the production
   choice until it provides a demonstrated performance benefit.
   Broader multi-delete and bounded-rank work remain deferred.

## Non-goals

- Changing KLIP centering, pixel normalization, exclusion boundaries, reference-selection semantics, or mode ordering.
- Applying ADI image exclusion to independent RDI libraries.
- Replacing direct KLIP solving for `includeRefNum`-limited reference sets.
- Introducing a projected/truncated Long--Males implementation without explicit scientific approval.
- Extending the mxlib SVD-deletion API; this work consumes the existing row/column-deletion interface.
