# KLIP Memory Management and Optimization Plan

## Goal

Bound KLIP memory use for large target and reference sequences without changing region geometry, reference selection,
PCA precision, requested mode counts, residual values, or final science products. Prefer removing retained full-frame
storage before reducing concurrency, because the existing KLIP timing baseline is eigensolver-dominated and therefore
has a real performance cost when fewer target images can be processed concurrently.

## Initial findings

KLIP currently has three materially different kinds of memory:

1. Full target, optional RDI reference, and residual cubes retained for the complete reduction.
2. Compact pixel-by-image arrays and shared PCA state retained for one region at a time.
3. Private selection and eigensolver scratch retained by each OpenMP worker.

These components scale differently. A useful memory policy must account for the configured output count, region pixel
count, reference-selection path, and worker count rather than using only the number of input frames.

### Full-sequence storage

Let:

- `P` be the number of pixels in a full image;
- `T` be the number of target images;
- `R` be the number of RDI references, or `T` for ADI bookkeeping;
- `M` be the number of entries in `klip.Nmodes`;
- `S_union` be the number of distinct full-image pixels covered by all configured regions.

`KLIPreduction::regions()` currently retains one full float target cube and allocates one zeroed full float residual cube
for every requested mode count before processing any region. The leading ADI full-cube term is therefore

```text
4 * P * T * (1 + M) bytes.
```

RDI adds a full `4 * P * R` reference cube. Unlike P4, current KLIP residuals do not have matching retained validity
cubes.

KLIP also always allocates `m_imsIncluded` as an `int` matrix with `T * R` entries, even when all references are used
and diagnostics are disabled. Its leading size is

```text
4 * T * R bytes.
```

For 8,182 ADI images this matrix is about 255.4 MiB. It is not needed by the science calculation after reference
selection; when diagnostics are enabled, the final file describes the selection made for the last processed region.

### Per-region shared storage

For a region containing `S` pixels, KLIP copies the target data into a compact `S * T` float array and mean-subtracts
it in place. RDI also copies the reference data into an `S * R` float array; ADI aliases the target array. These arrays
are released before the next region.

The covariance and shared-mode behavior depends on the library configuration:

| Path | Shared basis behavior | Important retained region allocations |
|---|---|---|
| All references, `R <= S` | Reference-space basis is calculated once | `R * R` float covariance and roughly `L * S` float KL modes |
| All references, `S < R` | Pixel-space basis is calculated once | `S * S` solve storage during construction and roughly `L * S` float KL modes; the `R * R` covariance is already skipped unless diagnostics require it |
| Exclusion or ranked inclusion | A basis is calculated for every target | Full `R * R` float covariance plus private selected-library and eigensolver storage per active worker |
| Right-reason projection | Adds dense spatial masks/projections | Shared `S * S` float mask and, on the shared-basis path, a shared `S * S` float projection matrix |

Here `L = min(klip maximum mode count, S, realized reference count)`. Enabling diagnostics can force the historical
full reference covariance even when the adaptive pixel-space solve would not otherwise need it.

### Per-worker storage

The OpenMP region is currently unconstrained by memory. On the target-specific basis path, each worker retains
reusable arrays at the largest dimensions encountered by that thread. For a target with `K` admissible references and
`D = min(S, K)`, the leading private allocations include:

- `S * K` floats for copied selected references;
- `K * K` floats for the selected reference covariance;
- approximately `L * S` floats for the KL modes;
- a double-precision `syevrMem` containing a `D * D` covariance, up to `D * L` eigenvectors, eigenvalues, and LAPACK
  workspace;
- additional float eigenvector and Gram arrays used by the higher-level adaptive solve; and
- another `S * S` float projection matrix when right-reason projection is enabled.

The two adaptive branches need separate estimates. Ignoring lower-order vectors and LAPACK's linear workspace, a
conservative description of the arrays simultaneously live around the solve is:

```text
K <= S: 4*S*K + 4*K*K + 4*L*S + 8*K*K + 12*K*L bytes
S < K:  4*S*K + 4*K*K + 12*S*S + 16*S*L bytes
```

Right-reason processing adds approximately `4 * S * S` private bytes on the target-specific path. These formulas are
a starting point for checked estimates, not a substitute for measured peak RSS: allocator state, Eigen temporaries,
BLAS/LAPACK workspaces, thread stacks, and the common reduction state must also fit.

On the all-reference path, KL modes are calculated once before the OpenMP target loop. Worker-local matrices remain
empty except for approximately `L` coefficients and `S` PSF values, so reducing target-loop concurrency there would
save little memory. The automatic policy must distinguish this case.

## Existing performance baseline

The existing `showTiming` report from a representative KLIP reduction measured:

| Stage | Measured time | Share of measured worker work |
|---|---:|---:|
| KLIP algorithm elapsed | 117.95 s | -- |
| Eigen decomposition | 2,213.45 worker s | 96.44% |
| KL image calculation | 25.12 worker s | 1.09% |
| PSF calculation/subtraction | 56.55 worker s | 2.46% |

The data geometry and peak RSS were not recorded with that run, so it is a timing baseline rather than a complete
memory baseline. It nevertheless shows why automatic worker reduction should follow fixed-storage reductions and
should activate only when the configured memory budget requires it.

## Representative large-sequence scaling

A 256-by-256 float cube with 8,000 planes occupies 1.953 GiB. Current leading ADI full-cube storage is therefore:

| Requested output mode counts | Target plus residual cubes |
|---:|---:|
| 1 | 3.91 GiB |
| 5 | 11.72 GiB |
| 10 | 21.48 GiB |
| 15 | 31.25 GiB |

At 8,182 planes, one such cube is 1.998 GiB. The current target-plus-residual term is about 23.97 GiB for 11 outputs
and 31.96 GiB for 15 outputs, before region arrays, reference selection, eigensolvers, or an RDI cube are counted.

The recent P4 large-data geometry provides an illustration of the potential compact-output reduction, not an assumed
KLIP configuration: it had 2,828 distinct processed pixels out of 65,536. At 8,182 planes, compact KLIP residuals for
that union would use about 0.086 GiB per output mode. Retaining the input cube, compact residuals, and one full residual
cube for finalization would give leading peaks of approximately 4.94 GiB for 11 outputs and 5.29 GiB for 15 outputs.
Those are conservative estimates with the target cube retained through finalization. If lifecycle review permits the
input cube to be released after the last region, the corresponding phase-by-phase fixed peaks fall to approximately
2.95 GiB and 3.29 GiB. Actual KLIP savings must be calculated from its configured region union.

## Constraints

- Preserve the current float input/residual representation and double-precision eigensolve.
- Preserve the adaptive choice between reference-space and pixel-space Gram matrices and the ordering of returned
  modes.
- Preserve configured region order. Regions are not currently required to be disjoint; when they overlap, a later
  region overwrites an earlier region's residual pixels. Compact storage must retain that last-region-wins behavior.
- Preserve zeros outside configured regions, post-median subtraction, derotation, combination, final plane order,
  headers, timing semantics, and per-frame PSF-subtracted filenames.
- Preserve the exact `imsIncluded.fits` diagnostic type and contents when diagnostics are enabled. Its in-memory
  representation when diagnostics are disabled is an implementation detail and need not remain a dense matrix.
- Treat `klip.writeDiagnostics=true`, RDI, and right-reason projection as explicit high-memory paths in estimates and
  tests; do not silently disable them.
- Respect a lower user-specified OpenMP limit. An automatic policy may reduce concurrency but must not exceed the
  process's existing maximum.
- Reject configurations before regression when checked size arithmetic overflows or the minimum one-worker path
  cannot fit the selected budget.
- Do not change reference inclusion, exclusion, or requested mode counts in response to memory pressure.

## Work sequence

### 1. Establish memory baselines and checked geometry

- [ ] Add a planning/reporting pass after input and mask geometry are known but before residual allocation. Build each
  region's index list once, record `S` for every region, and construct a full-index-to-compact-index map for the union.
- [ ] Record `P`, `T`, `R`, `M`, `S_union`, the maximum per-region `S`, ADI/RDI mode, selection path, right-reason and
  diagnostic flags, OpenMP and BLAS thread limits, and peak RSS for representative runs.
- [ ] Establish rapid local cases for the all-reference pixel-Gram path and a target-specific exclusion/inclusion path.
  Use fewer images and one or two annuli for turnaround, then validate accepted changes on the complete remote data.
- [ ] Retain the 117.95-second timing report as the pre-change performance reference and capture its missing geometry
  and `/usr/bin/time -v` peak RSS on the next comparable run.

### 2. Store residuals compactly and finalize one mode at a time

- [ ] When a final combination is configured and `output.outputPSFSub=false`, replace `M` full residual cubes with
  `M * S_union * T` compact floats. Workers continue writing disjoint target-image columns; sequential regions write
  through the union map so overlap behavior remains last-region-wins.
- [ ] Keep pixels outside the region union implicitly zero. Add explicit overlapping-region and masked-pixel tests
  comparing compact and retained residual materializations.
- [ ] After all regions finish, allocate one zeroed full residual cube, materialize one requested mode, and run the
  unchanged post-median, derotation, and combination lifecycle. Retain only its combined two-dimensional plane and
  release the temporary cube before the next mode.
- [ ] Preserve final output-plane order, headers, accumulated stage timings, and one aggregate progress message per
  lifecycle stage.
- [ ] Audit the post-regression lifetime of the full target and RDI cubes. Release them before one-mode finalization
  when no configured output or later lifecycle step consumes them; retain them when compatibility requires it.
- [ ] Keep the retained-cube path initially for `combine.method=none` and `output.outputPSFSub=true`, because those
  configurations explicitly request the residual sequence. Streaming the existing indexed products is a later step.
- [ ] Add checked byte calculations for compact residuals and one-mode materialization. Report and fail before KLIP
  work if the fixed compact path alone exceeds the configured budget.

### 3. Avoid unused dense inclusion bookkeeping

- [ ] Do not allocate `m_imsIncluded` when diagnostics are disabled and reference selection has no consumer for a
  retained matrix. Let target selection use a worker-local inclusion/index vector instead.
- [ ] When diagnostics are enabled, preserve the existing integer FITS product and its last-region selection semantics.
  Allocate or materialize the dense matrix only for that explicit diagnostic path.
- [ ] Measure whether a byte or packed internal representation followed by an exact integer diagnostic materialization
  lowers peak RSS. Do not change the diagnostic file's data type merely to save memory.
- [ ] Cover all-references, exclusion, each ranked-inclusion ordering, ADI, and RDI selection results in equivalence
  tests.

### 4. Add resource-aware worker selection

- [ ] Add `klip.memoryFraction` with the same high-level contract as `p4.memoryFraction`: default to a conservative
  fraction of Linux `MemAvailable`, let zero disable automatic limiting, and clearly report an unavailable discovery
  fallback.
- [ ] Before each region, account for full/compact retained state, target and RDI region arrays, shared covariance and
  modes, diagnostics, right-reason matrices, and a safety margin. Estimate regression and one-mode finalization as
  separate lifetime phases and budget their maximum rather than summing allocations that never coexist.
- [ ] Estimate target-specific worker scratch separately for `K <= S` and `S < K`. Use a conservative admissible
  reference bound derived from `R` and `includeRefNum`; exclusions may reduce actual `K` but must not be assumed to do
  so before selection.
- [ ] Apply `num_threads(effectiveWorkers)` only where target-specific bases create large private scratch. Keep the
  normal OpenMP maximum on the shared all-reference path unless measured memory identifies another private cost.
- [ ] Report available memory, budget, fixed/shared estimate, worst one-worker estimate, requested maximum, and
  effective worker count for every region to the console and FITS provenance.
- [ ] Fail with an actionable estimate if one worker cannot fit. Never change the science configuration to make the
  estimate fit.

### 5. Remove avoidable covariance and eigensolver copies

- [ ] Split reference selection from covariance extraction so it first returns the ordered admissible indices and
  realized `K`.
- [ ] On target-specific `S < K` solves, avoid forming and copying a selected `K * K` covariance that the adaptive
  pixel-space branch does not use. Form the `S * S` pixel Gram from the selected reference columns instead.
- [ ] Decide by measurement whether a shared full `R * R` covariance is worthwhile when realized libraries use the
  `K <= S` reference-space branch. Preserve the historical full covariance whenever diagnostics request it.
- [ ] Audit the nested `calcKLModes`/`calcEigenVecs` conversion path. KLIP intentionally solves in double precision
  from float inputs, but the reusable workspace and local float outputs can coexist with avoidable duplicate matrices.
  Remove only copies shown unnecessary by lifetime analysis and numerical equivalence tests.
- [ ] Consider indexed or bounded-tile pixel-Gram construction only after the dense selected-reference version is
  validated. Any summation-order change requires accepted residual and final-image tolerances.
- [ ] Re-measure peak RSS and elapsed/eigensolver worker time at one and several workers after each change; do not trade
  a large eigensolver regression for a small allocation reduction without evidence.

### 6. Stream explicitly requested residual products

- [ ] For `output.outputPSFSub=true`, process or materialize one mode at a time and write the existing
  `<prefix>_<reduction>_<image><extension>` products with unchanged `REDUCTION` indices, headers, and weight sidecar.
- [ ] Determine whether `combine.method=none` can use the same streaming path without changing the documented meaning
  of retaining PSF-subtracted products for postprocessing.
- [ ] Keep failure behavior recoverable and explicit when only part of a streamed sequence has been written.

### 7. Reconsider full input-cube retention only if still necessary

- [ ] Measure the compact-output and bounded-worker implementation before changing input I/O.
- [ ] If the target cube remains material, prototype a preprocessed-input mode that extracts the region union while
  reading each frame and retains only `S_union * T` target values. Apply the same design independently to RDI inputs.
- [ ] Preserve header, angle, timestamp, coadd, fake-injection, and preprocessing contracts. Do not silently substitute
  an annulus-only loader when full-frame preprocessing is still requested.
- [ ] Compare a union-once loader with per-region rereads; account for overlapping regions and repeated FITS I/O before
  selecting either design.

### 8. Document and validate the bounded path

- [ ] Add `klip.memoryFraction` to the user guide and KLIP configuration dictionary, including its FITS provenance
  mapping and explicit behavior on non-Linux hosts or when disabled.
- [ ] Document which output configurations use compact storage and which deliberately retain or stream residuals.
- [ ] Update the persisted timing schema only if new memory-planning or covariance stages are added to timing output;
  keep the existing console hierarchy otherwise.
- [ ] Run the complete KLIP and shared ADI/HCI test suites, the application tests, documentation build, and focused
  retained-versus-bounded FITS comparisons.
- [ ] For every edited function that calls mxlib, verify all directly called mxlib APIs against the current mxlib LCOV
  report. Record any gap under `Known non-blocking ownership follow-ups` in `agents/plans/mxlib_cleanup.md` as required
  by the repository coverage gate.

## Acceptance criteria

- A large ADI reduction reports a checked fixed/shared/per-worker memory estimate, selected worker count, and measured
  peak RSS with enough headroom to avoid swapping or allocation failure under the configured budget.
- The default final-combination path retains compact region residuals and at most one full residual cube during
  finalization, independent of `M`.
- Automatic worker selection distinguishes shared and target-specific basis paths, respects a lower user limit, and
  never changes reference selection or mode counts.
- Tests demonstrate compact-versus-retained equivalence for residual materialization and final products, including
  overlapping regions, masks, ADI/RDI, exclusions, ranked inclusion, repeated/clamped mode counts, post-median
  subtraction, derotation, and each combination method.
- Diagnostic tests preserve `imsIncluded.fits`, covariance, projection, and right-reason outputs when enabled.
- Existing `outputPSFSub` names, reduction indices, headers, and postprocess compatibility remain unchanged.
- Remote validation records complete configuration geometry, OpenMP and BLAS limits, timing, predicted memory, and
  measured peak RSS before the bounded path is considered complete.

## Deferred decisions

- The exact automatic-memory fraction and safety margin require measured KLIP RSS on both shared-basis and
  target-specific-basis runs; they must not be inferred only from P4.
- A byte or packed inclusion matrix is useful only if exact diagnostic compatibility can be retained without replacing
  it with an equally large long-lived integer materialization.
- Out-of-core or union-only input loading is deferred until compact residual storage and worker bounds reveal the
  remaining fixed-memory requirement.
- GPU eigensolvers, changes to PCA precision, fewer references, fewer modes, or altered region geometry are separate
  performance/science decisions rather than memory-management substitutions.
