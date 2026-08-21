# P4 Memory Management and Optimization Plan

## Goal

Reduce P4 wall-clock time and bound its memory use for large reductions without changing the regression definition,
interpolation contract, PCA precision, chosen modes, or science products. Work begins from measured detector-frame P4
behavior rather than assuming that eigensolves or additional temporal images are the limiting cost.

## Initial findings

The `showTiming` report and P4 timing schema 3 were measured on one representative remote-machine reduction. Times
below are aggregate OpenMP-worker seconds unless labeled elapsed; worker totals can exceed elapsed time because they
sum work across concurrent threads.

| Stage | Measured time | Share |
|---|---:|---:|
| P4 algorithm elapsed | 953.51 s | — |
| Total sampling | 40,730.28 worker s | 98.63% of measured worker work |
| Same-image target/OR sampling | 40,728.83 worker s | 99.996% of sampling |
| Additional-image PSF-disk sampling | 1.45 worker s | 0.004% of sampling |
| Gram construction | 476.42 worker s | 1.15% |
| Eigensolve/rank selection | 85.15 worker s | 0.21% |
| Projection/residual application | 2.56 worker s | 0.01% |

The prior schema-2 run had nearly the same elapsed time (953.47 s P4 algorithm) and sampling total (40,628.58 worker
s), so the sub-timing instrumentation is not a material source of overhead at this scale.

### Interpretation

- The configured additional-image `p4.numberImages` predictor disks are not the performance problem for this data.
- Detector-frame same-image OR sampling dominates. For every search pixel and retained target image, P4 constructs a
  local predictor row. Each same-image OR predictor currently uses a scalar cubic interpolation over a four-by-four
  footprint, so this stage repeatedly performs image loads and weighted arithmetic before the all-double PCA core.
- The approximately 42.7 aggregate sampling-worker seconds per P4 elapsed second initially appeared to show productive
  parallel interpolation. The local drill-down below instead found that those workers were repeatedly copying whole
  image planes at the sampling API boundary; the original locality/SIMD interpretation is superseded.

## Local reduced baseline (2026-08-20)

A short-turnaround baseline now uses the same AF Lep detector-frame cubic sampler and OR geometry as the full
configuration while retaining 96 of 621 frames and only the `[6,8)` annulus. This gives `T=96`, 96 search pixels,
and `K=800`; mode fractions `0.05,0.1,0.2` keep all requested counts positive without affecting sampling work.
`OPENBLAS_NUM_THREADS=1` was fixed for every run. The ignored local configuration is
`working/p4Reduce_afLepNaco_profile_local.conf`; outputs and `/usr/bin/time -v` records were written under `/tmp`.

The local host is an Intel Core i9-12900HK hybrid CPU with 20 logical CPUs, 14 physical cores, 24 MiB shared L3, and
62 GiB RAM. These single-run results establish turnaround and expose scaling behavior; they are not yet a variance-
controlled thread-count recommendation.

| OpenMP threads | P4 elapsed | Sampling worker time | Process CPU | Peak RSS |
|---:|---:|---:|---:|---:|
| 5 | 8.50 s | 40.82 s | 461% | 388 MiB |
| 8 | 7.81 s | 47.43 s | 585% | 393 MiB |
| 10 | 7.28 s | 54.99 s | 718% | 396 MiB |
| 10, repeat | 7.75 s | 55.85 s | 690% | 396 MiB |
| 20 | 7.68 s | 119.01 s | 1381% | 413 MiB |

All five final FITS products were byte-identical. Twenty threads did not improve wall time over ten and more than
doubled aggregate sampling time relative to five, consistent with contention, hybrid-core effects, or both. Repeated
runs and affinity-controlled P-core/E-core cases are still needed before selecting a local optimum.

Linux hardware-counter profiling is unavailable on this host because `perf_event_paranoid=4`. A GCC 15 Release
vectorization report shows that `P4PixelGrid::sample()` vectorized zero loops, and the fixed four-by-four accumulation
loop was rejected (`no vectype for stmt`) at the image/kernel multiply. The focused benchmark below showed that this
scalar arithmetic was secondary to an unintended full-image copy at the call boundary.

### Root cause and first fix (2026-08-20)

`eigenCube::image()` returns an `Eigen::Map` view, while `P4PixelGrid::sample()` accepted a concrete
`const Eigen::Array &`. Binding the map to that parameter implicitly evaluated it into a complete 256-by-256 temporary
array for every predictor sample. The reduced case therefore copied about 256 KiB before each four-by-four
interpolation. Changing the parameter to `const Eigen::Ref<const imageT> &` accepts both owned arrays and mapped cube
planes without materializing the image.

An ignored local harness, `working/P4PixelGrid_sampling_benchmark.cpp`, reproduced the production geometry with 96
synthetic images, 12 search pixels, 800 predictors, and 921,600 samples. Best-of-three single-thread Release results
were:

| Sampling path | Throughput | Equivalence to checked owned-array path |
|---|---:|---:|
| Checked API, owned `Eigen::Array` | 125.52 Msample/s | reference |
| Checked API, `Eigen::Map` through old concrete-array boundary | 0.227 Msample/s | bit-for-bit |
| Direct records, search-major | 200.95 Msample/s | bit-for-bit |
| Direct records, source-major blocks 1-12 | 198.06-200.12 Msample/s | bit-for-bit |
| Direct records, column-major footprint order | 195.47 Msample/s | not bit-for-bit |

The mapped-input path is approximately 553 times slower than the same checked API with an already-owned array. The
small source-major blocks provide no material benefit after the copy is absent, while changing footprint accumulation
order loses exact equivalence and is not faster.

The first reduced five-thread end-to-end run after the `Eigen::Ref` change produced:

| Measurement | Before | After | Change |
|---|---:|---:|---:|
| Total reduction | 9.126 s | 0.650 s | 14.0x faster |
| P4 algorithm elapsed | 8.504 s | 0.0568 s | 149.7x faster |
| Sampling worker time | 40.817 s | 0.0863 s | 473.1x lower |
| Peak RSS | 388 MiB | 388 MiB | no material change |

The output FITS file is byte-identical to all five pre-fix thread-sweep products. Sampling is now 35.6% of measured
worker work in this small case, versus 36.9% for Gram construction and 27.1% for eigendecomposition. The complete
remote configuration must be reprofiled before choosing another optimization target.

Five immediate follow-up runs at the same five-thread setting gave a mean P4 elapsed time of 0.0530 s with a 0.0040 s
sample standard deviation. Mean sampling, Gram, and eigendecomposition worker times were 0.0834 s, 0.0864 s, and
0.0651 s respectively. Their narrow ranges confirm that all three are now comparable at this deliberately small scale;
the final product hash remained unchanged.

### Post-fix cache/batching decision

Each interpolated predictor value is currently consumed once when constructing one search pixel's local PCA matrix;
retaining that value in a large cache does not itself remove interpolation work. A useful bounded tile must instead
provide at least one of these benefits:

- keep one source image and overlapping cubic footprints hot while sampling several nearby search pixels;
- arrange independent interpolations so the compiler can use SIMD without changing each interpolation's accumulation
  order beyond an accepted numerical contract; or
- reduce repeated checked lookup/address work by validating a block once and operating on compact records.

The focused benchmark compared source-image reuse directly and did not find a throughput improvement for blocks up to
12 search pixels. Do not add a sampled-value tile from the pre-fix profile. Reconsider bounded batching only if the
corrected complete-configuration profile again identifies sampling as the limiting stage.

## Corrected remote profile (2026-08-20)

The complete remote configuration was rerun after the zero-copy sampling change with the same data and science
configuration as the schema-3 `finim_0024.fits` baseline. The output product remained unchanged.

| Measurement | `finim_0024.fits` before | `finim_0025.fits` after | Change |
|---|---:|---:|---:|
| Total reduction | 966.54 s | 31.80 s | 30.4x faster |
| P4 algorithm elapsed | 953.51 s | 18.66 s | 51.1x faster |
| Sampling worker time | 40,730.28 s | 109.97 s | 370.4x lower |
| Gram worker time | 476.42 s | 524.05 s | now 72.28% of worker work |
| Eigensolve worker time | 85.15 s | 88.48 s | now 12.20% of worker work |

About 8.4 seconds remain after subtracting the individually reported loading, preprocessing, geometry, regression,
derotation, and combination stages. The current full-cube allocation/initialization and output lifecycle are among the
unreported work; they are not represented by the PCA worker subtimers. The zero-copy sampler is accepted. Sampling
batching is no longer the next optimization target.

## Memory model and finer-coadd requirement

The motivating data set improves significantly when its temporal coadd is reduced from five seconds to one second,
but the one-second reduction consumes about 27% of the remote workstation's memory. A further factor of five, to
approximately 0.2-second coadds, must not rely on manually selecting a sufficiently small OpenMP thread count.

Memory has two materially different components:

1. Retained data and output storage scales with frame count and is independent of worker count.
2. P4PCA scratch storage scales with both the annulus dimensions and the number of concurrent OpenMP workers.

For `P` pixels per frame, `F` retained frames, and `M` requested mode fractions, the current float storage includes one
input cube, `M` PSF-subtracted cubes, and `M` full float validity cubes. Its leading fixed term is therefore

```text
4 * P * F * (1 + 2*M) bytes.
```

The AF Lep configuration uses 256-by-256 frames and 11 mode fractions, so this term is exactly 5.75 MiB per retained
frame. The expected scaling is:

| Coadd | Approximate frames | Current fixed cube term |
|---:|---:|---:|
| 5 s | 621 | 3.49 GiB |
| 1 s | 3,105 | 17.43 GiB |
| 0.2 s | 15,525 | 87.17 GiB |

These values exclude PCA workspaces, masks, headers, final images, library allocations, and other process state. The
fixed term alone shows that thread limiting cannot make the next factor-of-five coadd reduction safe.

The search annuli cover radii 6 through 20 pixels, or only about 1,144 geometric pixels before discrete-grid and mask
effects. Retaining residual values only for owned search pixels changes their storage from `M*P*F` floats to roughly
`M*S*F` floats, where `S` is the number of owned search pixels. Rank validity is constant over the temporal samples of
one local fit and can be represented once per search-pixel/output pair rather than as another full float cube.

After regression, materialize one output mode at a time into one full science cube and one full validity cube, run the
unchanged post-median/derotation/combination/output lifecycle, retain its combined two-dimensional image, and release
the temporary cubes before materializing the next mode. At 0.2-second coadds, the leading peak for this design is
approximately:

| Storage | Approximate size |
|---|---:|
| Input float cube | 3.79 GiB |
| Compact residuals for 11 modes and about 1,144 search pixels | 0.73 GiB |
| One full science/validity cube pair | 7.58 GiB |
| Leading fixed peak | 12.10 GiB |

This estimate is deliberately incomplete, but it reduces the dominant fixed term by more than sevenfold and makes
worker concurrency a useful remaining control.

For one detector-frame worker, let `T` be the target-image count, `K` the predictor count, and `D=min(T,K)`. After
removing avoidable eigensolver copies, the leading double-precision worker storage is approximately
`8*(T*K + 2*D*D + T*M)` bytes plus vectors and LAPACK work arrays. The predictor matrix grows linearly with finer
coadds even after `D` reaches `K`, so worker count still must be bounded from an explicit memory estimate.

### Initial implementation and local validation

The first compact path is active when a final combination is requested and `output.outputPSFSub=false`. It retains
annulus-only residuals and one byte of rank support per search-pixel/output pair, materializes one full output pair at
a time, and releases the pair after combination. `combine.method=none` and `output.outputPSFSub=true` deliberately
retain the legacy full cubes because those configurations request the residual sequence itself. Streaming indexed
per-frame output remains a follow-up rather than changing its filename contract in this first slice.

The 96-frame, one-annulus, three-mode local AF Lep case produced the same final science data and header as the retained-
cube zero-copy baseline before new memory-provenance cards were added. Peak RSS fell from about 388 MiB to 297,212 KiB
(about 290 MiB). The reduction now records the intentionally variable memory snapshot and worker choice in FITS, so
subsequent complete-file hashes differ by those new cards while `fitsdiff` reports no data differences.

The compact-versus-retained unit regression also exercises `numberImages=1` with one-sided endpoint selections and
compares every combined output value, including invalid pixels. The complete local suite passes all 21 test programs,
and the generated documentation builds successfully with only the repository's pre-existing warnings.

`p4.memoryFraction` defaults to `0.8`. On Linux, P4 budgets that fraction of `MemAvailable`, verifies that compact
residuals plus one materialized output pair fit, estimates private PCA storage for each annulus, and limits the P4
OpenMP region with `num_threads(effectiveWorkers)`. A lower `OMP_NUM_THREADS` value is respected. Zero disables the
policy; unavailable Linux discovery warns and leaves the OpenMP maximum unchanged. The local validation selected all
five requested workers with a 1.99 MiB estimate per worker.

### Full remote validation (2026-08-21)

The compact implementation completed the AF Lep reduction with more than 8,000 images at 0.2-second coadds, 256 by
256 pixels, 14 annuli, and 2,828 owned search pixels. At regression start it reported 99.3809 GiB available, a
79.5047 GiB automatic budget, 1.24867 GiB of compact residuals, and a 3.99512 GiB one-mode materialization. The
largest per-worker estimate was 657.726 MiB, so all 48 permitted workers fit and no automatic reduction was needed.
P4 regression took 14,700.7 seconds and the complete reduction took 14,934.4 seconds. This validates that the bounded
path can process the finer-coadd data set that motivated the work. The 0.2-second science product had lower SNR than
the one-second product, so one-second coadds remain preferable for that data independently of memory feasibility.

Because bounded finalization applies the shared ADI lifecycle once per output mode, its original progress reporting
printed `derotating` and `combining` repeatedly. The lifecycle now supports suppressing per-call progress, and compact
P4 reports each aggregate stage once while retaining per-mode processing and accumulated timing.

## Constraints

- Preserve the current detector-frame interpolation values and predictor-column ordering bit-for-bit where practical;
  any floating-point reordering must be verified against accepted numerical tolerances and science metrics.
- Keep P4PCA's all-double normal equations and eigensolve unchanged unless separately justified.
- Preserve final FITS planes, per-frame PSF-subtracted products, validity semantics, post-median subtraction,
  derotation, and combination behavior when compact residual storage changes intermediate object state.
- Do not make a full sampled `search pixels × target images × predictors` cache by default: it can exceed memory for
  the very data sets this work targets.
- Do not retain full-frame residual and validity cubes for every mode. Temporary materialization must be bounded to one
  mode unless an explicitly selected memory policy permits more.
- Do not infer a beneficial OpenMP or BLAS thread count from one machine; measure it.

## Work sequence

### 1. Establish reproducible baselines

- Record complete configuration, image geometry, annulus search-pixel counts, `K` predictor counts, retained target
  count `T`, OpenMP thread count, BLAS thread count, elapsed timing report, and peak RSS for representative datasets.
- Repeat each benchmark enough times to report variance. Run a small thread-count sweep with BLAS constrained to one
  thread while sampling dominates.
- Use the existing P4 timing schema 3 to retain same-image and additional-image sampling contributions per run.
- Use the 96-frame, one-annulus local configuration for rapid iteration, but revalidate any selected optimization on
  the complete remote configuration before acceptance.

### 2. Characterize same-image sampling

- [x] Reproduce the dominant cost in a focused sampler harness and separate mapped-image conversion from cubic
  arithmetic.
- [x] Inspect compiler vectorization; hardware counters remain unavailable locally because of host policy.
- [x] Remove the implicit full-image temporary by accepting an `Eigen::Ref` view in `P4PixelGrid::sample()`.
- [x] Re-run the complete remote configuration with the zero-copy interface before adding finer timers.
- Retain the current timing granularity unless finer instrumentation changes measured throughput materially.
- [x] Compare interpolated values per second for current search-major and source-image/block-major paths without
  invoking P4PCA or final processing.

### 3. Remove avoidable per-worker PCA copies

- [x] Call the underlying double `eigenSYEVR` directly on the already-double Gram matrix. The higher-level mxlib helper
  copies the Gram matrix and retains duplicate eigenvector storage in the reusable workspace even though P4 requests
  neither type conversion nor normalization.
- [x] Reuse the rotated-frame caller predictor matrix for temporal centering because it is refilled before every local fit.
  This removes one `T*K` allocation from that path; detector-frame P4 is already uncentered and receives no benefit.
- [x] Add direct equivalence tests and verify the production detector- and rotated-frame outputs before accepting the
  changes. Record peak RSS, not only wall time.

### 4. Store residuals compactly and finalize one mode at a time

- [x] Allocate annulus-local compact residual arrays indexed by output mode, owned search pixel, and retained target image.
  Parallel workers write disjoint search-pixel columns.
- [x] Store rank support once per local fit and mode. Reconstruct temporal validity from each annulus's existing central
  target-image selections; do not allocate a float validity value for every full-frame pixel and image.
- [x] After every annulus is complete, allocate one zeroed full-frame residual cube and validity cube. Materialize one
  output mode, write its optional detector-frame validity diagnostic, and apply the existing ADI final lifecycle.
- [x] Preserve output-plane order while processing one mode at a time. Accumulate derotation and combination timing
  across modes and write the multi-plane final FITS cube only after all combined images have been collected.
- [ ] Stream the existing `outputPSFSub` filenames and `REDUCTION` indices one mode at a time. The current implementation
  deliberately retains the legacy full-cube path whenever these per-frame products are requested.
- [x] Release each materialized pair before allocating the next. Retain only the small combined output cube when
  combination is enabled.
- [x] Add overflow checks for every compact and full-cube size calculation and fail before regression when even the fixed
  one-mode path cannot fit the configured budget.

### 5. Add resource-aware worker selection

- [x] Add a high-level P4 memory policy with a conservative automatic default and an explicit disable/override. Derive an
  effective worker count separately for every annulus from available memory, its `T`, `K`, `D`, mode count, compact
  residual allocation, and the measured/estimated one-worker scratch requirement.
- [x] Respect a lower user-specified OpenMP limit. The automatic policy may reduce concurrency but must never exceed the
  process's existing OpenMP maximum.
- [x] Apply the selected value with `num_threads(effectiveWorkers)` on the P4 regression parallel region. A single worker
  remains the minimum; fail early with an actionable estimate if one worker cannot fit.
- [x] Report the memory budget, fixed-storage estimate, per-worker estimate, requested maximum, and effective worker count
  for every annulus to the console, FITS provenance, and checked diagnostics.
- [x] Prefer Linux `MemAvailable` for automatic headroom with a clearly reported fallback when automatic discovery is
  unavailable. A separate explicit byte/GiB override is deferred until measurements show a need beyond fraction zero.

### 6. Reconsider input-cube retention only if still necessary

- Measure the compact-output implementation at 1-second and 0.2-second coadds before changing input loading.
- If the remaining input cube is material, prototype annulus-local or bounded image-plane loading from already
  preprocessed files. Account for repeated I/O across annuli and do not combine this with the first compact-output
  change.
- Keep fake-planet pixel-local processing and a general out-of-core preprocessing pipeline as separate designs.

## Acceptance criteria

- The 1-second and proposed 0.2-second configurations report predicted and measured peak RSS, effective worker count,
  and enough headroom to avoid host swapping or allocation failure under the selected budget.
- The default path retains compact annulus residuals and at most one full-frame residual/validity pair at a time; it
  never allocates a full sampled predictor tensor.
- Tests demonstrate equivalence of residual cubes when materialized, validity cubes, final products, and optional
  per-frame products, including masked/edge, one-sided temporal fallback, and `numberImages > 0` cases.
- Timing schema/version and user documentation are updated if persisted timing fields change.
- Automatic thread/resource selection is explicitly reported and never changes requested science settings.

## Deferred decisions

- GPU acceleration is not justified by the current measurements until a CPU locality/vectorization baseline is known.
- The exact automatic-memory fraction and safety margin require remote peak-RSS measurements at one and several
  workers; they must not be selected solely from local-machine throughput.
- Bounded sampling tiles are deferred because the corrected profile made Gram construction dominant and compact
  residual storage addresses the actual frame-count memory scaling.
- Reducing OR geometry or `p4.numberImages` can reduce work but changes the regression model; it is a science/config
  choice, not an implementation optimization.
