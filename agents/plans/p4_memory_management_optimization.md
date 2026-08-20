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
- The approximately 42.7 aggregate sampling-worker seconds per P4 elapsed second shows that many workers are already
  productive. Optimization should therefore prioritize data locality, memory bandwidth, and SIMD-friendly sampling
  before PCA backend, GPU, or temporal-selection changes.

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

### 2. Characterize same-image sampling

- Add focused, opt-in profiling only if needed to distinguish cubic-footprint loads, interpolation arithmetic, matrix
  writes, and finite-value promotion within same-image sampling.
- Inspect compiler vectorization and hardware counters to determine whether the dominant limit is memory bandwidth,
  cache misses, or scalar arithmetic.
- Retain the current timing granularity unless finer instrumentation changes measured throughput materially.

### 3. Evaluate bounded annulus-local batching

- Prototype processing a bounded block of search pixels at a time, organized so an input image and predictor-offset
  footprint are reused across multiple local fits before eviction.
- Consider a compact tile of sampled same-image predictor values, with its capacity selected from an explicit memory
  budget. The tile must be released before the next block and must not retain a full-annulus predictor tensor by
  default.
- Compare two layouts: search-pixel-major storage for direct P4PCA matrix assembly and image/predictor-major staging
  for contiguous sampling. Select from measurements, not intuition.
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
