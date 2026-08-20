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

## Constraints

- Preserve the current detector-frame interpolation values and predictor-column ordering bit-for-bit where practical;
  any floating-point reordering must be verified against accepted numerical tolerances and science metrics.
- Keep P4PCA's all-double normal equations and eigensolve unchanged unless separately justified.
- Do not make a full sampled `search pixels × target images × predictors` cache by default: it can exceed memory for
  the very data sets this work targets.
- Continue to support a bounded-memory path. Cache/tile capacity must be explicit and must degrade safely to the
  existing streaming sampler.
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
- Re-run the complete remote configuration with the zero-copy interface before adding finer timers.
- Retain the current timing granularity unless finer instrumentation changes measured throughput materially.
- [x] Compare interpolated values per second for current search-major and source-image/block-major paths without
  invoking P4PCA or final processing.

### 3. Evaluate bounded annulus-local batching

- Defer implementation until the corrected remote baseline shows that sampling still warrants batching.
- Prototype processing a bounded block of search pixels at a time, organized so an input image and predictor-offset
  footprint are reused across multiple local fits before eviction.
- Consider a compact tile of sampled same-image predictor values, with its capacity selected from an explicit memory
  budget. The tile must be released before the next block and must not retain a full-annulus predictor tensor by
  default.
- Compare two layouts: search-pixel-major storage for direct P4PCA matrix assembly and image/predictor-major staging
  for contiguous sampling. Select from measurements, not intuition.
- Reject a tile that only stores one-use sampled results without improving source-footprint reuse, SIMD utilization,
  or checked-address overhead.
- Verify exact predictor ordering, mask behavior, edge-footprint rejection, and final residual/validity products
  against the streaming implementation.

### 4. Optimize the selected sampler

- If batching provides locality but remains arithmetic-limited, vectorize the fixed four-by-four cubic interpolation
  across a block of search pixels or images while preserving its kernel and footprint semantics.
- If it is bandwidth-limited, prefer cache blocking, compact addresses/weights, and prefetch-friendly traversal over
  more threads.
- Keep additional-image PSF-disk sampling separate; it is currently negligible and should not drive the design.

### 5. Resource-aware operation

- Add an optional memory-budget policy that chooses a safe tile size and reduces OpenMP workers when per-worker
  workspace plus the selected tile would exceed available resources.
- Report the chosen tile size, memory estimate, and effective worker count in diagnostics without silently changing
  scientific configuration.
- Defer out-of-core input loading and fake-planet pixel-local processing to separate designs after bounded in-memory
  tiling is measured.

## Acceptance criteria

- A benchmark record identifies whether the optimized sampler improves wall time, aggregate worker time, or both for
  at least one representative reduction, including peak RSS and variability.
- The default path remains bounded and does not allocate a full sampled predictor tensor.
- Tests demonstrate equivalence of predictor sampling, residual cubes, validity cubes, and final products for the
  chosen tile sizes, including masked/edge and `numberImages > 0` cases.
- Timing schema/version and user documentation are updated if persisted timing fields change.
- Thread/resource selection is opt-in or explicitly reported and never silently changes requested science settings.

## Deferred decisions

- GPU acceleration is not justified by the current measurements until a CPU locality/vectorization baseline is known.
- The best cache layout, tile dimension, and memory budget are hardware- and dataset-dependent.
- Reducing OR geometry or `p4.numberImages` can reduce work but changes the regression model; it is a science/config
  choice, not an implementation optimization.
