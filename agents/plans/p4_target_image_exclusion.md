# P4 Target-Image Exclusion and GPU Plan

## Goal

Make detector-frame P4 genuinely out-of-sample in time: when predicting target pixel value `y[t]`, do not use the
same temporal sample in either the optimization-region basis or the regression target. Support an optional wider
temporal exclusion set without making fine-coadd reductions computationally unusable.

The first investigation will take the GPU path seriously. It will establish an exact batched-CUDA implementation as
a correctness and performance baseline, then determine whether production-scale execution should use direct batched
eigensolves, the Long and Males (2021) SVD downdate, or a hybrid of the two.

This work applies to detector-frame P4 only. Rotated-frame P4 remains a documented negative-result experiment.

## Required science operation

For one search pixel, the current code constructs

- `X`, a `T x K` matrix of optimization-region predictors; and
- `y`, a length-`T` vector containing the target-pixel time series.

It then performs one PCA regression using all rows and returns all `T` residuals. Consequently, prediction of `y[t]`
uses both the predictor row `X[t,:]` and the exact target value `y[t]` in the fit.

For target `t`, define an exclusion set `E_t`. The new regression is trained on

```text
X_train = X with rows E_t removed
y_train = y with entries E_t removed
```

and applied only to the held-out predictor row `X[t,:]`. The required product remains one residual per requested mode,

```text
r_t(m) = y[t] - prediction(X[t,:], X_train, y_train, m).
```

The minimum operation is `E_t={t}`. A later temporal radius may add nearby rows, but it must always contain `t`; there
is no wrap at either end of the sequence. Output-plane mode counts remain fixed for the annulus. A target whose
training set has insufficient structural or numerical rank is marked invalid for the affected output modes rather
than silently changing the meaning of that plane.

This definition initially treats an image as excluded when it is the central row of `X`. If `p4.numberImages > 0`, a
frame can also appear in another row's additional-image predictor columns. We must either exclude every row whose
selection references the target frame or explicitly document central-row exclusion. Since the additional-image
experiment is being reconsidered independently, the first exact implementation should use `p4.numberImages=0` and
make other values a deliberate validation error until that semantic decision is made.

## Current implementation boundary

The relevant path is compact and gives us a clean replacement point:

1. `P4Reduction::fitDetectorSearch()` samples `X` and `y` once for each search pixel.
2. `P4PCA::calculate()` forms either `X X^T` when `T <= K`, or `X^T X` and `X^T y` when `K < T`.
3. The current LAPACK `eigenSYEVR` call finds the requested largest eigenpairs in double precision.
4. `P4PCA` projects every training row and returns a `T x M` residual array.
5. The outer P4 loop is OpenMP-parallel over search pixels, with one large PCA workspace per worker.

The exclusion path needs only `T x M` held-out residual scalars and per-target/per-mode validity. It should not return
Gram matrices, eigenvectors, or predictor-space coefficients to the host.

## Useful temporal-Gram identity

The common current regime is `T <= K`. Define the full temporal Gram matrix

```text
C = X X^T.
```

For delete-one or delete-block exclusion, the training Gram matrix is simply the principal submatrix

```text
C_train = C[not E_t, not E_t].
```

The vector of inner products between the held-out row and all training rows is also already present:

```text
k_t = C[t, not E_t].
```

If `C_train u_j = lambda_j u_j`, the prediction after retaining `m` modes is

```text
prediction_t(m) = sum over j=1..m of
                  (k_t^T u_j) * (u_j^T y_train) / lambda_j,
```

using the largest modes and the same rank-tolerance contract as current P4. This avoids forming or downloading
predictor-space eigenvectors. It also means all delete-one matrices have the same dimension `T-1` and are directly
batchable. A wider exclusion region produces a small number of matrix dimensions near dataset boundaries; targets can
be grouped by retained dimension before batching.

When `K < T`, the smaller matrix is the predictor Gram matrix

```text
G = X^T X,
c = X^T y,
G_train = G - sum over i in E_t of X[i,:]^T X[i,:],
c_train = c - sum over i in E_t of X[i,:]^T y[i].
```

That is a rank-`|E_t|` downdate of fixed dimension `K`. It is algebraically convenient but becomes a very large dense
eigensolve for the fine-coadd case. It is the regime in which a true SVD downdate is most important.

## Hard GPU assessment

### What CUDA can do well

Current cuSOLVER provides `cusolverDnXsyevBatched`, including FP64 matrices and eigenvectors, for a contiguous batch of
same-sized symmetric matrices. This is a direct match for tiles of delete-one temporal Gram matrices. cuSOLVER calls
are asynchronous where possible and accept a stream, so sampling/assembly, eigensolve, projection, and transfer can be
pipelined. The library requires device and host workspace to be allocated by the caller, which also lets hciReduce
budget memory explicitly.

The initial CUDA baseline should therefore:

1. form or upload one full `C` per search pixel;
2. assemble a device-resident tile of principal submatrices and held-out cross-vectors;
3. solve the tile with `cusolverDnXsyevBatched` in double precision;
4. calculate rank and all requested held-out predictions on the device; and
5. copy back only residuals, validity/rank summaries, and timing counters.

For a wider exclusion set, tiles are grouped by matrix dimension. The device assembly kernel receives explicit retained
row indices, so endpoints need no special solver path.

### What streams do and do not solve

Use a bounded ring of pipeline slots, each owning a CUDA stream, events, device buffers, and any required pinned host
buffer. One GPU-owner thread schedules work; the existing OpenMP search-pixel loop must not create an independent set
of CUDA handles and large workspaces per CPU worker. A small CPU producer can sample the next search pixel while the
GPU processes the current one.

Streams can overlap independent small computations and data transfer when the device permits it. They cannot remove
the approximately `T` decompositions required by literal leave-one-out processing, and concurrent large eigensolves
may already saturate the GPU. The benchmark must separately report useful overlap and time spent waiting for a free
pipeline slot.

For production integration, prefer uploading the float input cube and compact interpolation geometry once, then
sampling `X` on the GPU. Sending a double `T x K` matrix for every search pixel would otherwise move tens to hundreds
of gigabytes in a complete reduction. The first solver benchmark may accept preformed matrices in order to isolate
cuSOLVER performance, but that is not the final data path.

### Workload reality check

The representative `finim_0048.fits` product records 621 retained images for most annuli, predictor counts up to 4,150,
and 5,024 total search pixels. Literal exclusion would perform 3,062,704 dense eigensolves across those recorded
annulus dimensions. One `620 x 620` FP64 matrix occupies about 2.93 MiB; a tile of 64 occupies about 188 MiB before
eigenvalues and solver workspace. The memory is manageable on a workstation GPU, but the aggregate cubic work remains
large.

At more than 8,000 images, `K` may become the smaller dimension. A single `4,000 x 4,000` FP64 matrix is about 122 MiB,
and a brute-force decomposition for every held-out frame and every search pixel is not credible as the production
algorithm. A favorable batched-GPU result at `T=621` must not be extrapolated as if it changed that scaling.

The Long and Males SVD downdate was developed for this exact repeated-reference-removal problem. It updates a
precomputed SVD after removing frames and reports near-linear runtime scaling, with published speedups increasing from
2.6x at 300 frames to 140x at 10,000 frames. The likely production direction is therefore:

```text
full SVD once per search pixel
    -> small delete-row downdate for each target/exclusion block
    -> held-out projection only
```

with GPU kernels and streams accelerating the full SVD, downdate batches, and projections. The direct batched
eigensolver remains valuable as an independent exact oracle and may be the faster implementation for modest `T`.

References:

- Long and Males (2021), [Unlocking starlight subtraction in full data rate exoplanet imaging by efficiently updating
  Karhunen-Loève eigenimages](https://arxiv.org/abs/2101.11634)
- NVIDIA, [cuSOLVER API reference](https://docs.nvidia.com/cuda/cusolver/)

## Proposed CUDA architecture

Keep CUDA optional and isolated behind a target-exclusion solver interface. CPU-only builds and the current in-sample
path must remain available.

```text
P4Reduction geometry and selections
                |
                v
bounded CPU producer or GPU sampler ----> pipeline slot N
                                           - stream + events
                                           - full Gram/SVD state
                                           - exclusion tile
                                           - solver workspace
                                           - residual/status output
                ^                                      |
                |                                      v
          slot N+1 prepared                 copy T x M results only
```

Build support should use an opt-in CMake option, `check_language(CUDA)`, `find_package(CUDAToolkit)`, and explicit
`CUDA::cudart`, `CUDA::cublas`, and `CUDA::cusolver` targets. Non-CUDA translation units see only a C++ interface and a
runtime availability query. An explicit CUDA request must fail clearly if the backend or device is unavailable; it
must not silently change algorithms.

The GPU memory budget should be a configurable fraction of `cudaMemGetInfo()` free memory. From the maximum annulus
dimensions and solver workspace query, choose the tile size at runtime. Reserve space for the input cube, geometry,
full Gram/SVD state, at least two pipeline slots when feasible, output buffers, and a safety margin. A tile size of one
must remain valid. The tile calculation must also enforce cuSOLVER's current
`n * leadingDimension * batchSize <= INT32_MAX` constraint rather than relying on memory capacity alone.

All initial science arithmetic remains FP64, matching current `P4PCA`. Mixed precision is a separate todo and must be
validated independently.

## Work sequence

### 1. Freeze semantics and build a CPU oracle — implemented

- Add an internal exact held-out solver that explicitly gathers `X_train` and `y_train`, uses the existing `P4PCA`
  eigensolver, and evaluates only `X[t,:]`.
- Cover delete-one first and represent exclusion as an explicit row-index set so a temporal window does not require
  another solver API.
- Preserve current eigenvalue ordering, rank tolerance, realized-mode counts, NaN/validity behavior, and all-double
  arithmetic.
- Validate the temporal-Gram held-out formula against explicit predictor-space coefficients.
- Keep this slow path available in tests and small diagnostic runs as the numerical oracle.

The production CPU oracle uses the existing KLIP configuration names and meanings: `adi.excludeMethod` selects
`none`, `pixel`, `angle`, or `imno`, and `adi.minDPx` supplies the inclusive threshold. Non-`none` methods therefore
exclude the target itself at zero. Each target receives an explicit ordered training-row set. Per-target/per-mode
validity is retained through compact output materialization. The initial supported surface is detector-frame,
`p4.numberImages=0`, ordinary full reductions without local or frozen-PSF products.

The first implementation gathered each target's training matrix, called the existing FP64 `P4PCA` solve, and
evaluated only the held-out row. A 48-worker remote run confirmed that this literal oracle kept the CPU occupied but
recomputed `X_train X_train^T` for every target and was only about 40% complete after 14 hours. The `T <= K`
production path now uses the temporal-Gram identity above: it forms `C = X X^T` once per search pixel, extracts each
training-set principal submatrix and held-out cross-vector, and evaluates the residual without constructing
predictor-space coefficients. Independent explicit refits remain the numerical test oracle, and the `K < T` path
retains explicit predictor-space refitting until its downdate implementation is available. This optimization removes
the repeated `T x K` copies and Gram products but still performs approximately `T` dense eigensolves per search pixel;
it is an intermediate direct baseline, not a replacement for the Long-Males downdate.

### 2. Create a representative standalone CUDA benchmark

- Generate and ingest dimension distributions matching `finim_0048.fits`, plus synthetic `T>K` cases representative
  of the 8,000-frame data.
- Benchmark full-matrix formation, exclusion-tile assembly, `cusolverDnXsyevBatched`, projection, and transfers
  separately with CUDA events.
- Sweep tile size, one versus multiple streams, and full batched eigenpairs versus a streamed selected-eigenpair
  cuSOLVER call if the latter wins despite lacking a batched API.
- Compare every residual and rank decision with the CPU oracle; do not use solver throughput alone as acceptance.
- Record GPU model, driver/runtime/cuSOLVER versions, dimensions, requested modes, tile size, peak allocation, and
  effective matrices per second.

The local CUDA 13.3 toolkit is sufficient to compile this benchmark, but the local NVIDIA driver is unavailable.
Execution and performance decisions therefore require the remote GPU machine.

#### Initial resident-decomposition harness (2026-08-26)

`p4DecompositionBenchmark` now provides the first deliberately narrow measurement. It generates a deterministic
Gram-like FP64 batch, preloads it on the selected GPU, warms every CPU/GPU solver workspace, and reports:

- serial CPU `eigenSYEVR` throughput with one OpenBLAS thread;
- OpenMP CPU batch throughput with one independent LAPACK workspace per simultaneous solve;
- one contiguous `cusolverDnXsyevBatched` call on one CUDA stream; and
- the same resident matrix batch divided among independently owned concurrent streams and cuSOLVER workspaces.

Host-to-device preload, device-to-device restoration of the destructively overwritten matrices, allocations,
workspace queries, and validation downloads are all outside the timed intervals. Timed GPU values include solver
launch overhead and synchronization because those are real costs of processing one resident batch. The report also
compares GPU eigenvalues with the CPU result and checks the GPU residual and orthogonality for the first matrix.

Build and run the AF Lep-sized default case with:

```bash
cmake -S . -B _build_cuda_benchmark \
    -DCMAKE_BUILD_TYPE=Release \
    -DHCIREDUCE_BUILD_CUDA_BENCHMARKS=ON
cmake --build _build_cuda_benchmark --target p4DecompositionBenchmark -j
./_build_cuda_benchmark/benchmarks/p4DecompositionBenchmark
```

Useful sweeps change `--matrices`, `--streams`, and `--cpu-threads`; `--dimension`, `--predictors`, `--repetitions`,
`--device`, and `--seed` are also explicit. `--cpu-only` permits CPU-path validation on a host without a working CUDA
device. The initial defaults use dimension 620, 2,480 synthetic predictors, 64 matrices, four streams, and three
repetitions. The source builds successfully with CUDA 13.3 on the local host and the small CPU-only path executes.
The first complete CPU/GPU measurement is recorded below.

#### First measured result: RTX 3050 Ti laptop GPU (2026-08-26)

The default resident benchmark was run with 64 `620 x 620` FP64 matrices, three repetitions, 2,480 synthetic
predictors, 20 simultaneous CPU decompositions, and four CUDA streams. OpenBLAS was fixed at one thread per CPU solve.
Transfers and destructive-input restoration remained outside every timed interval.

| Path | Mean batch time | Throughput | Relative to one-stream GPU |
|---|---:|---:|---:|
| Serial CPU | 2.7396 s | 23.36 matrices/s | 0.61x |
| 20-worker OpenMP CPU | 0.4698 s | 136.23 matrices/s | 3.58x |
| One cuSOLVER batch, one stream | 1.6805 s | 38.08 matrices/s | 1.00x |
| Four concurrent cuSOLVER batches/streams | 1.7377 s | 36.83 matrices/s | 0.97x |

The one-stream batch used 827.59 MiB of device workspace. Splitting the same work among four streams used essentially
the same total workspace but reduced throughput by 3.4%, indicating that this eigensolver batch already saturates the
useful GPU resources on the tested RTX 3050 Ti. Stream concurrency is not beneficial for this matrix size/device.

Numerical agreement is excellent: maximum relative CPU/GPU eigenvalue difference was `9.99e-16`, GPU normalized
residual was `7.05e-17`, and normalized orthogonality error was `1.44e-16`. FP64 cuSOLVER is therefore a suitable
numerical comparison oracle even though this direct algorithm is not competitive.

Applying the measured throughputs to the 3,062,704 decompositions implied by the `finim_0048.fits` geometry gives
solver-only lower bounds of approximately 6.24 hours for the 20-worker CPU batch, 22.34 hours for one GPU batch, and
23.10 hours for four GPU streams. These exclude Gram construction, exclusion-matrix assembly, projection, and all
data movement. Direct repeated eigendecomposition is rejected as the production target-exclusion algorithm on this
hardware. Continue with the Long-Males SVD downdate; retain this harness for numerical validation and measurements on
other CPU/GPU architectures.

### 3. Implement the temporal-Gram CUDA baseline — rejected before integration

The resident solver benchmark already gives a 22.34-hour lower bound for the reference geometry on the tested GPU,
before matrix assembly or projection. Do not implement this path in `P4Reduction`. The following design is retained
only in case a future GPU architecture changes the measured conclusion substantially:

- Add the optional CUDA backend and a single GPU-owner scheduler.
- Initially support detector-frame, `numberImages=0`, delete-one exclusion, and `T<=K`.
- Keep full temporal Gram data on the device, batch principal-submatrix construction, solve, project, and download only
  `T x M` residuals plus compact validity.
- Add wider exclusion sets by batching targets with equal retained dimension.
- Integrate GPU timing beneath P4 Sampling/Gram/EigenDecomposition/Projection without changing the existing high-level
  timing schema meaning.

### 4. Make an evidence-based algorithm decision — complete for the direct path

Before expanding the direct path, compare measured and extrapolated end-to-end time with these gates:

- a complete 621-frame reduction is operationally useful, not merely faster than a deliberately slow CPU oracle;
- memory stays within its declared GPU budget for every annulus;
- an 8,000-frame extrapolation is supported by scaling measurements rather than cubic wishful thinking; and
- transfers and CPU sampling do not erase solver gains.

The direct method missed even the 621-frame gate. Retain the standalone harness as an oracle and architecture probe,
but do not add a direct eigensolver backend to the reduction. Proceed to the SVD downdate.

### 5. Adapt the Long-Males SVD downdate

The detailed reusable mxlib design, exact factor-deletion equations, P4 factor-space regression, and revised CPU-first
work sequence are in [p4_svd_downdate.md](p4_svd_downdate.md). That plan supersedes the preliminary implementation
details below where they differ.

- Derive the held-out P4 regression directly from the downdated SVD, including numerical-rank truncation and multiple
  requested mode counts.
- Implement and test the row-removal update on CPU first against explicit recomputation.
- Batch the small downdate problems and held-out projections on the GPU; benchmark whether the initial full SVD is
  better on cuSOLVER or the existing CPU LAPACK path for each dimension regime.
- Support delete-block temporal exclusion without repeated single-row updates when a block update is numerically and
  computationally preferable.
- Define a periodic or residual-based recomputation policy if successive downdates accumulate unacceptable numerical
  error. Each target should normally be derived independently from the same full SVD, not from a long chain of
  destructive downdates.

### 6. Extend all P4 consumers

- Normal compact and retained-cube output paths receive held-out residuals directly.
- Pixel-local and negative-companion optimization reuse a persistent GPU context, but trial-dependent `X` and `y`
  invalidate any affected Gram/SVD state.
- PSF response calculation cannot retain one `K x M` coefficient matrix per search pixel anymore: coefficients vary
  with target `t`. Calculate each held-out PSF response on the device and retain only the requested response stamp or
  combined product.
- Jackknife refits must state whether target exclusion applies inside each block refit and avoid nesting an accidental
  brute-force exclusion loop inside the current optimizer loop.

## Configuration and provenance direction

Final names should be selected when the CPU oracle is integrated, but the configuration needs three separate ideas:

- exclusion mode: current in-sample behavior versus exclusion of the target row;
- optional temporal radius and its unit (frame index, elapsed time, or parallactic angle); and
- execution backend: `auto`, `cpu`, or `cuda`.

Elapsed time is the scientifically preferable temporal-window coordinate for irregular cadence, while a frame-count
radius is useful for testing and compatibility with existing image-neighbor concepts. If both are supported they must
be mutually exclusive. No endpoint wrapping is permitted.

FITS provenance must record at least the exclusion definition/unit/value, backend and algorithm (`explicit`,
`batched-eigen`, or `svd-downdate`), precision, device name, CUDA/cuSOLVER versions, peak GPU allocation, and actual
tile sizes. `HIERARCH P4 IN SAMPLE` should become zero for these products.

## Verification matrix

- Exact small random matrices for both `T<=K` and `K<T`.
- Delete-one and delete-block cases, including both endpoints and irregular temporal windows.
- Rank-deficient, repeated, zero, and nearly dependent predictor rows around the configured rank tolerance.
- Every requested mode supported, partially supported, or unsupported after exclusion.
- CPU oracle versus temporal-Gram identity, CUDA direct solve, and SVD downdate.
- Normal reduction, compact storage, retained cubes, PSF response, pixel-local optimizer, and jackknife paths.
- CUDA build with a device, CUDA build without a usable device, and a completely CPU-only build.
- Bounded GPU allocation and tile-size-one fallback under an artificially small memory budget.
- Repeated-run tolerance and explicit confirmation of whether the chosen cuSOLVER routines are deterministic on the
  tested software/hardware stack; require numerical equivalence, not bitwise CPU/GPU identity.

## Acceptance criteria

1. A held-out result agrees with explicit refitting within a documented FP64 tolerance and never uses `y[t]` or row
   `X[t,:]` in its training fit.
2. The configuration, FITS provenance, validity, and output-plane meanings are unambiguous.
3. CPU-only builds and legacy in-sample reductions remain supported and unchanged when exclusion is disabled.
4. GPU memory is bounded from an explicit estimate and discovered free memory; CPU thread count cannot multiply GPU
   workspaces accidentally.
5. The selected production algorithm has measured scaling on both the 621-frame reference case and a representative
   multi-thousand-frame case.
6. PSF response and optimizer products are computed with the same held-out model rather than quietly reverting to the
   in-sample coefficients.

## Current recommendation

Complete the exact CPU leave-out oracle in step 1, then move directly to the Long-Males SVD downdate in step 5. Use
the resident benchmark to validate future machines, but do not integrate direct repeated eigendecomposition into P4.
GPU work should now target batched downdates and held-out projections, with streams reconsidered only after measuring
those smaller operations; the full batched eigensolver did not benefit from concurrent streams on the first device.
