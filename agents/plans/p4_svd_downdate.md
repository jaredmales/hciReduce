# Reusable SVD Downdate and P4 Integration Plan

## Goal

Replace P4's repeated target-excluded eigendecompositions with one base factorization per search pixel followed by
small row-deletion downdates. Put the numerical row/column-deletion operation in mxlib's math library with a contract
that is useful outside hciReduce. A complete-rank base must preserve P4's exact exclusion semantics, mode planes,
rank threshold, validity, and held-out regression definition. A truncated base is explicitly approximate because the
global subspace was selected with the held-out row present; it must have different provenance and acceptance gates.

The first implementation is CPU/FP64. The existing exact held-out implementation remains the numerical oracle and
fallback during development. The direct GPU eigendecomposition benchmark and current exact temporal-Gram baseline
remain documented in [p4_target_image_exclusion.md](p4_target_image_exclusion.md).

## Decisions

- The reusable mxlib system will remove rows or columns from the matrix represented by supplied thin-SVD factors.
  The deletion itself is exact for those factors; truncating the base SVD is the caller's only model approximation.
- mxlib will distinguish a stable complement-preserving small-SVD realization from a covariance/leading-spectrum
  realization. The latter is the first P4 backend because P4 needs only predictor-side directions and squared leading
  singular values; it must not be presented as a full-relative-accuracy SVD routine.
- Only a complete-rank P4 base is labeled exact. A bounded-rank backend is labeled projected/approximate, and it cannot
  become the strict-exclusion default without an explicit science decision based on leakage and recovery tests.
- The published Long--Males reduced core will be implemented in the prototype or test harness as a separately named
  comparison backend. It must not be silently conflated with exact factor deletion.
- Every target exclusion is derived independently from one immutable base factorization. Downdates are never chained
  from target to target, so there is no periodic-recompute interval and no accumulated cross-target drift.
- A multi-row exclusion is one rank-`c` operation. It is not a sequence of destructive state updates.
- mxlib owns the linear-algebra contract, validation, workspaces, and solver status. P4 owns exclusion construction,
  science-mode selection, rank tolerance, invalid-output policy, timing, memory budgeting, and FITS provenance.
- CPU dense eigendecomposition is the correctness backend. A structured diagonal-minus-low-rank solver is a later
  performance backend, not a prerequisite for establishing the API and equations.
- The new path will not initially expand the current held-out surface: detector-frame P4, `p4.numberImages=0`, no
  pixel-local processing, and no frozen-PSF products.

## Why the reusable contract differs from the reduced Long--Males core

Long and Males (2021) write observations as columns. For a working-rank factorization

$$
A \approx U_q \Sigma_q V_q^T
$$

and a set `E` of deleted columns, let `V_E=V_q[E,:]`. After discarding the orthogonal-complement terms from Brand's
general low-rank SVD modification, their reduced core is

$$
K_{\mathrm{LM}}
  = \Sigma_q - \Sigma_q V_E^T V_E
  = \Sigma_q M_E,
\qquad
M_E = I - V_E^T V_E.
$$

The small SVD

$$
K_{\mathrm{LM}} = W\,\widetilde\Sigma\,Y^T
$$

gives the approximate updated left basis `U_q W`; KLIP does not need to form the updated right basis. The paper
recommends working at rank `q=k+c` when `c` columns are deleted, then retaining only `k` output modes. It also states
that the operation modifies the truncated representation, omits the changed subset mean, and is not numerically
identical to a fresh fit.

That reduced core includes an additional approximation that matters for P4. Whenever the deleted-side factor does
not span its complete coordinate space--including a truncated factor, and a full compact factor with fewer columns
than rows--`M_E` is generally not a projector, so

$$
K_{\mathrm{LM}}K_{\mathrm{LM}}^T
  = \Sigma_q M_E^2 \Sigma_q
  \ne \Sigma_q M_E \Sigma_q.
$$

The matrix on the right is the retained covariance represented by the supplied factors. Therefore the literal
Long--Males rotation does not, in general, diagonalize P4's held-out training covariance and cannot be inserted into
the existing independent-mode principal-component regression sum without changing its mathematics.

A one-row example must lock this distinction into the tests. For `X=[1;1]`, with `T=2` and `K=1`, the full compact
SVD has `sigma=sqrt(2)` and temporal factor `L=[1/sqrt(2);1/sqrt(2)]`. Removing either row leaves singular value one.
The symmetric core returns `lambda=1`, while the literal reduced core returns singular value `1/sqrt(2)`. Conversely,
when `T<=K` and the compact temporal factor is square orthogonal, the deleted-row outer product is a projector and
the two cores coincide. Both invariants belong in mxlib's shape tests.

The reusable mxlib operation should instead have this explicit contract:

> Given a thin SVD, return the requested one-sided singular system after physically deleting specified rows or
> columns from the matrix represented by those factors.

For full-rank factors this is the ordinary exact deletion problem. For truncated factors it is exact deletion from
the rank-`q` represented matrix and, equivalently, a Rayleigh--Ritz approximation to deletion from the original
matrix. A faithfully reduced Long--Males implementation can remain available under an explicitly projected name for
reproducibility and comparison.

## Generic row/column-deletion equations

Let

$$
A = U\Sigma V^T,
$$

where `U` is `m x q`, `V` is `n x q`, their columns are orthonormal, and the diagonal singular values in `Sigma` are
finite, nonnegative, and descending.

### Removing rows

Let `E` contain `c` row indices, and let `U_E=U[E,:]`. If `S` denotes the retained rows, then

$$
A_S^T A_S
  = V H_E V^T,
$$

with the small symmetric core

$$
Z_E = \Sigma U_E^T,
\qquad
H_E = \Sigma^2 - Z_E Z_E^T
    = \Sigma\left(I-U_E^TU_E\right)\Sigma.
$$

Solve

$$
H_E = W_E\Lambda_EW_E^T,
$$

with eigenvalues descending. The retained matrix's right singular directions are `V W_E`, and its singular values
are `sqrt(lambda_E)`. If a caller needs the retained left factor, its positive-singular-value columns are

$$
U'_S = U_S\Sigma W_E\Lambda_E^{-1/2}.
$$

### Removing columns

Column removal is the transpose-dual operation. With `V_E=V[E,:]`,

$$
H_E = \Sigma^2 - (\Sigma V_E^T)(\Sigma V_E^T)^T.
$$

The preserved left singular directions are `U W_E`. Positive-singular-value columns of the retained right factor
can be reconstructed only when requested.

The hot operation therefore needs only the singular values and rows of the singular-vector factor on the deleted
side. It returns a small preserved-side rotation and spectrum; it must not eagerly multiply large factors.

### Stable complement-preserving realization

The symmetric `H_E` core is ideal when only the leading preserved-side system is needed, but forming a covariance
squares singular-value conditioning. The general mxlib API also needs a small-SVD realization that retains Brand's
deleted-side orthogonal complement without materializing a large complement basis.

For row deletion, let `F=U_E` and construct a rank-revealing factor `B` such that

$$
B^T B = I_c-F F^T.
$$

The factor can be obtained from the small `c x c` matrix with a rank-revealing eigendecomposition or QR-compatible
square root. Define the at-most-`(q+c) x q` core

$$
K_E=
\begin{bmatrix}
  (I_q-F^TF)\Sigma \\
  B F\Sigma
\end{bmatrix}.
$$

It obeys

$$
K_E^T K_E
  = \Sigma(I_q-F^TF)\Sigma
  = H_E.
$$

If

$$
K_E=P_E\Sigma'_E Y_E^T,
$$

then the retained predictor-side rotation is `Y_E`, so the right singular directions are `V Y_E`. Column deletion
uses the transpose-dual construction and returns the preserved left-side rotation. This core is an exact Brand
deletion of the represented factors, avoids a full rectangular refit, and is the stable generic backend when small
singular values matter. The covariance backend and stable backend must agree on invariant leading results within
their documented numerical tolerances.

## P4 equation restructuring

### Definitions

For one detector search pixel, P4 has

$$
X \in \mathbb{R}^{T\times K},
\qquad
y \in \mathbb{R}^{T},
$$

where each row of `X` is the optimization-region predictor vector for one temporal sample and `y` is the target-pixel
time series. Use the working-rank SVD

$$
X \approx L\Sigma R^T,
$$

with

- `L`: `T x q` temporal singular factor;
- `Sigma`: `q x q` diagonal singular values; and
- `R`: `K x q` predictor-space singular factor.

For target `t`, `E_t` is the configured exclusion set and contains `t`; `S_t` is its complement. Let `c_t=|E_t|`.
The current detector-frame P4 fit is uncentered, so no subset-mean approximation is involved.

### One base factorization per search pixel

Let `m_max` be the largest realized output mode in the annulus and `c_max=max_t(c_t)`. The first truncated working
rank to test is

$$
q = \min\left(\operatorname{rank}_{\mathrm{struct}}(X),\ m_{\max}+c_{\max}\right).
$$

The `+c_max` rule is the Long--Males guard-rank heuristic, not an error bound. Complete-rank validation retains every
safely representable positive singular triplet. This gives `q=min(T,K)` only for a structurally full-rank design;
rank-deficient designs use their complete represented numerical rank under a base-safety criterion distinct from
P4's science `rankTolerance`. Convergence tests must also try extra guard modes before the truncated path is accepted.

Construct the base factors without a full rectangular SVD:

- If `T <= K`, form `C=X X^T`, find its largest `q` eigenpairs, set `L` to the temporal eigenvectors, and set
  `sigma_j=sqrt(lambda_j)`.
- If `K < T`, form `G=X^T X`, find its largest `q` predictor eigenpairs `R`, set
  `L=X R Sigma^{-1}`, and then discard `R` after the required temporal scores have been formed.

Division by an exactly zero or numerically unsafe singular value is forbidden. The full-data P4 science rank
threshold must not be used to discard base modes before a downdate: deletion can change eigenvalue ratios, and the
held-out rank rule applies to each deleted core. Retain the structurally requested base modes whenever they can be
formed safely; otherwise report the condition and use the defined invalid/fallback policy.

Precompute once

$$
g=L^Ty.
$$

P4's residual-only path needs `L`, `Sigma`, and `g`; it does not need to retain `R`.

### One exclusion downdate

For target `t`, gather

$$
L_E=L[E_t,:],
\qquad
y_E=y[E_t].
$$

Form

$$
Z_E=\Sigma L_E^T,
\qquad
H_E=\Sigma^2-Z_EZ_E^T,
$$

and find its requested largest eigenpairs

$$
H_E W_E = W_E\Lambda_E.
$$

The retained target projection in temporal-factor coordinates is

$$
a_E
  = L_{S_t}^T y_{S_t}
  = g-L_E^Ty_E.
$$

For the held-out predictor row, the predictor-basis coordinates follow directly from the SVD:

$$
x_t^T R = l_t\Sigma.
$$

Define the two rotated coordinate vectors

$$
d_E = W_E^T\Sigma a_E,
\qquad
h_t = W_E^T\Sigma l_t^T.
$$

The `m`-mode prediction and residual are

$$
\widehat y_t^{(m)}
  = \sum_{j=1}^{m}\frac{h_{t,j}d_{E,j}}{\lambda_{E,j}},
\qquad
r_t^{(m)}=y_t-\widehat y_t^{(m)}.
$$

This is the existing held-out principal-component regression written entirely in the base factor space. It does not
construct `X_train`, a principal temporal-Gram submatrix, a `K`-element coefficient vector, or updated predictor
directions. One cumulative sum produces every configured mode plane.

When `q` spans the complete safely represented rank of `X`, this is algebraically equivalent to an explicit fit of
`X[S_t,:]` in both `T<=K` and `K<T` regimes. With truncated `q`, it is the held-out regression constrained to the
global rank-`q` predictor subspace. Because that global subspace was selected using every row, including `t`, this is
not strictly out-of-sample even after deleting `t` inside the subspace; it must be labeled approximate and assessed
against the exact oracle.

### Structural and numerical support

For every target, structural support is bounded by

$$
r_{\mathrm{struct},t}=\min(K,T-c_t,q).
$$

The current strict numerical-rank rule is preserved exactly:

$$
\lambda_j > \texttt{rankTolerance}\,\lambda_{\max}.
$$

Do not compare singular values directly with `rankTolerance*sigma_max`; the equivalent singular-value threshold
would contain `sqrt(rankTolerance)`. A requested plane above either the structural or numerical rank remains NaN with
`sampleValidity=0`.

The initial dense backend computes and validates the complete `q`-value spectrum, even if it returns rotations only
for the largest requested modes. This is required both to diagnose a negative tail and to report the exact numerical
rank of the represented model. A later selected-spectrum optimization must perform a separate smallest-eigenvalue/PSD
check and must mark its rank as capped when it cannot count the entire supported spectrum.

Extend the result contract rather than overloading the current integer silently:

- `numericalRank` remains the rank of the model actually solved;
- a rank-model field distinguishes an exact complete-base fit from a projected/truncated-base fit;
- `numericalRankCapped` states that the integer is only the supported rank observed by a selected solve; and
- the annulus aggregate clears the capped flag if any target attaining the reported minimum has an exact count.

`P4 MINIMUM RANK` remains the minimum integer for backward readability, accompanied by explicit rank-model,
base-rank, full-base, and capped cards. Exact/direct products retain their existing meaning.

Rank-threshold comparisons also need an ambiguity policy. Estimate an eigenvalue uncertainty from solver residuals,
machine precision, dimension, and the core construction scale. If a requested boundary lies within that guard band
of `rankTolerance*lambda_max`, use the explicit retained-matrix oracle for that target rather than allowing harmless
rounding differences to flip a validity bit. Well-conditioned equivalence tests use condition-scaled prediction and
residual tolerances; a fixed decimal tolerance is not sufficient for every spectrum.

If `T-c_t=0`, P4 does not call mxlib: it preserves the current behavior of rank zero, NaN residuals, and
`sampleValidity=0` for every mode, and includes that target in minimum-rank/rank-invalid accounting. The generic
mxlib convenience API may reject deletion of its complete side without turning this valid P4 edge case into a
reduction exception.

## Compact exclusion representation

`P4Reduction::targetReferenceRows()` currently stores the complement of every exclusion set. This is `O(T^2)` index
storage even when only one or a few rows are excluded. At `T=8000`, it can exceed 500 MiB per annulus before vector
overhead.

Before integrating the downdate, replace the production representation with an `ExclusionPattern` that stores only
removed rows:

- `imno` uses inclusive/exclusive index spans where possible, with truncated endpoint spans and no wrap;
- `angle` and `pixel` use sorted unique explicit indices; and
- every representation exposes an arbitrary-index view accepted by the same mxlib API.

The exact oracle can traverse a lazy complement or gather retained rows on demand. Tests must continue to compare the
compact representation with current KLIP-compatible `adi.excludeMethod` and inclusive `adi.minDPx` semantics.
Count compact exclusion storage as annulus-persistent memory in the P4 budget; the current retained-row vectors are
not included there.

## Proposed mxlib API and ownership

### Files

The first mxlib change should be self-contained:

- `include/math/svdDowndate.hpp`
- `source/math/svdDowndate.cpp`
- `tests/include/math/svdDowndate_test.cpp`
- `include/math/CMakeLists.txt`
- `source/math/CMakeLists.txt`
- `include/math/math.hpp`
- the relevant Doxygen group definition and mxlib bibliography

The strict mxlib test inventory will require the same-stem test as soon as the public header/source is added.

### API layers

Names are provisional, but the separation of responsibilities is not:

1. `svdDeletionLeadingCore()` accepts descending singular values, the deleted-side factor rows, requested output rank,
   and a caller-owned workspace. It uses `H_E`, validates the full spectrum, and returns squared singular values,
   their safe square roots, and the requested preserved-side rotation.
2. `svdDeletionStableCore()` builds the complement-preserving `K_E`, performs its small SVD, and returns the
   preserved-side rotation and spectrum without a normal-equation accuracy claim.
3. `svdRemoveRows()` accepts `U` plus arbitrary sorted unique row indices, gathers `U_E`, and dispatches to an explicit
   backend. Its returned rotation applies to `V`.
4. `svdRemoveColumns()` accepts `V` plus arbitrary sorted unique column indices, gathers `V_E`, and dispatches to an
   explicit backend. Its returned rotation applies to `U`.
5. The high-level generic row/column API defaults to the stable core. P4 explicitly requests the leading covariance
   backend while its numerical contract is being validated.
6. Optional helpers may materialize updated factors. They are not used by P4 and should not be part of the hot path.
7. A separately named `svdProjectedModificationCore()` may implement the literal Long--Males `K_LM` operation for
   reproduction tests. It must document that different approximation contract.

The base-factor construction is intentionally separate. Callers may obtain factors from Gram eigenpairs, LAPACK,
randomized SVD, or a future CUDA backend without changing the deletion API.

Use concrete column-major dynamic `Eigen::Array<realT,...>` matrix/vector aliases, matching mxlib's existing LAPACK
interfaces, and `std::span<const Eigen::Index>` for deletion indices. Float and double are the only supported scalar
types. A result object owns its rotation and spectrum, and the caller may reuse it across calls; a workspace owns only
scratch. `result.prepare(q,maxOutputRank)` and `workspace.prepare(q,maxDeleted,backend)` reserve documented capacities
before the hot loop. Results are readable only after success and may be unspecified after a solver failure, matching
the existing mxlib numerical style. Accessors return unaligned borrowed views into result-owned storage; those views
remain valid until the result is prepared, assigned, moved, or destroyed and never refer to call-local storage.

### Workspace and threading

`svdDeletionWorkspace<realT>` is caller-owned, non-shared, and either safely move-only or explicitly nonmovable. It
owns the gathered factor rows, `Z`, symmetric core, complement factor, small SVD core, solver-side rotation/spectrum
scratch, and reusable LAPACK work memory. The result separately owns published outputs. Capacity preparation ensures
the target loop performs no routine allocation after warmup.

The existing `mx::math::syevrMem` owns raw buffers and is implicitly copyable. Keep the first change self-contained by
using new RAII `std::vector`/Eigen LAPACK work buffers and calling the low-level wrappers, rather than embedding or
copying `syevrMem`. Hardening that existing type is separate cleanup unless implementation evidence makes it necessary.

mxlib does not choose OpenMP or BLAS thread counts. P4 initially keeps its current outer parallelism of one search
pixel per OpenMP worker, but a serial inner solve is a requirement to enforce or verify, not an assumption. Benchmark
runs set and record `OPENBLAS_NUM_THREADS=1` (or the vendor equivalent), record the linked BLAS configuration and
effective process thread count, and detect oversubscription before interpreting timings. If hciReduce adds
vendor-conditional thread control, it does so once at the application level before the OpenMP region; mxlib does not
mutate global library state. Alternative outer/inner thread splits remain explicit hciReduce benchmarks.

### Validation and status

The public API validates before invoking LAPACK:

- finite float/double inputs and dimensions representable by `MXLAPACK_INT`;
- finite, nonnegative, descending base singular values;
- consistent factor width and requested output rank;
- sorted, unique, in-range deletion indices;
- a nonempty retained side for the row/column convenience APIs; and
- finite solver output in descending order.

Orthonormal columns on the supplied thin singular factor are a required precondition for the exact-deletion identity.
A separately callable base validator checks orthonormality once at a documented tolerance; the per-deletion hot core
does not repeat an `O(nq^2)` check. Optimized finite checks use `mx::math::isFinite()` so they remain valid under
mxlib's fast-math build.

The result/status distinguishes invalid input, allocation/preparation failure, LAPACK query/solve failure, materially
non-PSD core, and success with roundoff clamping. It retains the underlying LAPACK status. `prepare()` reports
allocation failure without leaving a falsely adequate capacity. No generic numerical-rank threshold is hidden in
mxlib; callers apply their own policy.

After explicitly symmetrizing `H_E`, clamp only negative eigenvalues consistent with a documented construction-error
scale such as `C*epsilon*q*(sigma_max^2+||Z_E Z_E^T||)`. Scaling only by the possibly cancellation-small norm of
`H_E` is insufficient. A more negative eigenvalue is a factor-consistency or cancellation failure and must be
reported, not unconditionally set to zero. High deleted-row leverage near one is reported as a rank-loss diagnostic,
not treated automatically as an error.

### Precision scope

Provide explicit float and double instantiations so the API is reusable, but integrate P4 with double only. The dense
`H_E` path uses normal equations and therefore squares singular-value conditioning. It is appropriate for P4's
thresholded leading spectrum and mirrors the current covariance-based method. A general full-spectrum consumer that
needs high relative accuracy in the smallest singular values uses the complement-preserving small-SVD backend; a
later Gu--Eisenstat structured backend improves its asymptotic cost. These backend-specific guarantees must be in the
public docs.

Do not use the existing `eigenGESDD()` in the target loop. It computes a full SVD, destructively modifies input,
allocates on every call, and currently lacks the workspace/error contract needed here.

## P4 integration boundary

### `P4PCA`

- Retain `P4PCA::calculateHeldOut()` as the exact direct oracle.
- Add a separate downdated implementation behind an explicit backend seam; do not branch deep inside the current
  direct loop.
- Factor base-spectrum preparation and factor-space response calculation into focused private helpers.
- Return the same `P4PCAResult`: `T x M` residuals, `T x M` sample validity, per-mode status, and rank diagnostics.
- Keep fallback deliberate and internally consistent. If one target has an ambiguous rank boundary, recompute every
  target row for that search pixel with the explicit oracle and record both the reason and the number of target rows.

### `P4Reduction`

- Build compact exclusion patterns once per annulus and calculate `c_max` directly.
- Allocate one mxlib downdate workspace per P4 worker.
- Include `L`, the base Gram/factorization scratch, the `q x q` core, rotations, deleted-row scratch, and LAPACK memory
  in `estimatedWorkerBytes()` and the `p4.memoryFraction` worker limit.
- In the `K<T` branch, budget the simultaneous peak while forming the base state (`X`, `G`, predictor eigenvectors
  `R`, temporal factor `L`, and solver workspace), not only the later steady `L+core` footprint.
- Preserve the current output cubes, compact/retained equivalence, validity propagation, mode-plane meanings, and
  finalization path.
- Keep the current validation errors for rotated-frame, `numberImages>0`, local, and frozen-PSF combinations.

### Configuration and provenance

During validation, expose the implemented exact choices

```text
p4.exclusionSolver=explicitRefit|factorDowndateExact
p4.deletionBackend=leadingCovariance|rankOneSecular
```

`factorDowndateExact` retains the complete safely represented rank. A future `factorDowndateProjected` choice would
use bounded `q` and be explicitly approximate; it is not currently accepted as a configuration value. The explicit
backend remains the default until the gates below pass. Do not add a user-visible guard-rank control until convergence
measurements show that it is scientifically necessary; sweeps can initially use a diagnostic override. If it becomes
public, define it as extra modes beyond `m_max+c_max`, not as a second interpretation of mode fractions.

`leadingCovariance` remains the compatibility deletion backend and supports one or more deleted rows.
`rankOneSecular` is an explicit one-row-only selection and is valid only with `factorDowndateExact`. Validate the
realized compact exclusion pattern before regression: every active target must delete exactly one row. Preserve an
all-rows-excluded target as the existing rank-zero invalid result, and reject any other multirow structured request
rather than silently switching backends.

Record at least:

- exclusion solver;
- base rank `q`, maximum requested mode, and `c_max`;
- whether `q` spans the full structural rank;
- dense or structured mxlib backend;
- any roundoff clamps, failures, or explicit-fallback count; and
- capped/minimum supported-rank information.

Record the selected mxlib backend for `factorDowndateExact` and `none` for the explicit oracle. Keep separate
search-pixel counts for rank-boundary, factor-validation, and structured deletion-solver fallback.

### Timing

Split held-out timing sufficiently to distinguish:

- base Gram construction;
- base eigendecomposition and temporal-factor formation;
- deletion-core assembly;
- reduced eigensolve;
- factor-space projection/residual accumulation; and
- explicit fallback, if enabled.

Preserve existing aggregate timing fields or advance the timing schema deliberately, updating console YAML,
`p4Timing.fits`, readers, and tests together.

The structured integration advances `p4Timing.fits` to schemas 6/7/8. The common columns insert reusable-base-factor,
row-deletion, and explicit-oracle-fallback worker seconds after the existing eigensolve aggregate. Base-factor and
row-deletion time subdivide the eigensolve aggregate. Explicit-fallback time spans the entire direct-oracle retry and
therefore overlaps the retry's Gram, eigensolve, and projection phase counters; it is diagnostic rather than an
additional component of total worker time. Schemas 7 and 8 append the existing local-PSF and reconstruction fields,
respectively.

## Current P4 integration and review checkpoint

The first P4 integration slice now has a concrete code boundary, but it is not yet a promoted performance backend:

- The first reusable mxlib SVD deletion implementation and focused test were committed as `f4aea38`; the ABI-v2
  correction described below supersedes that shared-library boundary and was committed as `b9bac77`.

- `P4TargetExclusions` represents `imno` deletions as per-target spans and arbitrary `angle`/`pixel` deletions as
  flattened index storage. The direct oracle materializes retained complements only on demand.
- `P4PCA::calculateHeldOut()` remains the direct retained-matrix oracle. The separate
  `P4PCA::calculateHeldOutDowndated()` path constructs a complete safely represented base factor and calls mxlib's
  configured double row-deletion backend independently for every target.
- `p4.exclusionSolver=explicitRefit` remains the default. `factorDowndateExact` is an explicit opt-in and has no
  unrecorded direct-refit fallback. If the base temporal factor fails only its orthogonality validation, or if any
  deleted spectrum is within a dimension- and scale-aware FP64 guard band of the configured rank threshold, the
  whole search pixel is recomputed with the explicit oracle so every residual, rank, and validity decision remains
  internally consistent. A recoverable numerical failure of the selected structured deletion solver uses the same
  whole-search-pixel fallback; invalid API contracts, dense deletion failures, and failures of the explicit oracle
  still reject the reduction. The initial unsupported-combination restrictions remain unchanged.
- P4 has per-worker mxlib deletion workspaces, deletion-aware memory estimates, and FITS cards for solver, rank model,
  full-base state, deletion backend, compact exclusion storage, base rank, maximum deletion count, clamp count, and
  explicit-fallback count. `P4 EXPLICIT FALLBACKS` is the number of target rows recomputed in each annulus, so one
  fallback search pixel contributes `T`; reason-specific search-pixel counts and the maximum rejected-factor defect
  with its paired tolerance are recorded separately. The deletion-solver fallback has its own search-pixel count.
  The configured solver and selected deletion backend are recorded separately; the direct oracle records deletion
  backend `none`.

The first external-consumer verification exposed an Eigen ABI hazard rather than a downdate-equation failure. mxlib
was compiled with `-march=native`, while the ordinary hciReduce build was generic. On the 48-core ROC run this first
appeared as `free(): invalid pointer` while worker workspaces were unwinding. All installed and build-tree libraries
contained the same mxlib source hash, ruling out stale headers or libraries. ASan/UBSan then located the earlier
undefined behavior: mxlib executed AVX aligned loads/stores against storage that satisfied only the generic
consumer's 16-byte Eigen alignment contract. Rebuilding hciReduce with matching native alignment made the focused
suite and reduced application clean, which isolated the failure from OpenMP and the eigensolver.

Opaque out-of-line result/workspace state alone was not a sufficient fix because the compiled shared-library API
still carried Eigen reference objects and alignment-sensitive header implementations across the boundary. The
second-generation ABI keeps the convenient Eigen-facing source API only in hidden inline consumer adapters. Those
adapters pass standard-layout, trivially-copyable pointer/dimension/stride/index descriptors to versioned compiled
entry points; mxlib maps all borrowed arrays explicitly unaligned and uses `DontAlign` for its owned Eigen storage.
The result and workspace types also carry an ABI-v2 tag, so mixing old and new headers/libraries fails at symbol
resolution before a differently laid-out object can be constructed. The implementation translation unit hides its
Eigen inline symbols as additional hardening. The argument-order convention remains explicit in `mxlib/AGENTS.md`:
output-only arguments precede semantic inputs, analogous to `y = f(x)`, while reusable input/output workspaces may
follow the inputs.

The permanent CMake fixture compiles the real shared-library producer translation unit with an explicit 32-byte Eigen
alignment contract and its consumer translation unit with a 16-byte cap. It exercises under-aligned inputs and
numerical reads from returned rotations and spectra, plus empty and malformed ABI descriptors. The ordinary focused
target and `MXLIB_ONE_TEST` path use the same mismatch, and a Linux symbol test requires the ABI-v2 exports while
rejecting obsolete Eigen-facing, untagged-handle, and consumer-adapter exports. A separate forced
64-byte-producer/16-byte-consumer test exercised both cores, both row-deletion backends, preparation, validation,
output views, and moves. The legacy GNU Make source list and target-scoped visibility flag are also covered by a
direct object build.

At this checkpoint the focused mxlib suite passes 399 assertions in 16 test cases with 841/841 executable
SVD-downdate lines covered. Matched and mismatched optimized hciReduce consumers each pass 855 assertions in 25
focused `P4PCA` cases, and the final staged mismatch passes all 56,598 assertions in 20 `P4Reduction` cases. The
mixed-alignment ASan/UBSan consumer also passes the 855 `P4PCA` assertions, and the reduced application completes all
four regressions and writes its final output without a sanitizer diagnostic.

A `factorNotOrthonormal` status after this memory-safety fix is a separate numerical question rather than another
allocator symptom. P4 measures `max|L^T L-I|` and retains the accepted tolerance. It does not use the rejected factor:
the existing exact held-out oracle recomputes every target row for that search pixel and the reduction continues. The
fallback reason and numerical measurements remain visible in console and FITS provenance. Compare `SYEVR`, `SYEVD`,
and a direct-SVD oracle as a later accuracy/performance experiment if the fallback rate or worker-time cost is
nontrivial; do not use a solver change merely to mask the earlier ABI failure.

The first post-ABI remote run reached annulus `[10,12)` before reporting `max|L^T L-I|=1.0274e-11` against the
mxlib automatic bound `64 epsilon 621=8.8249e-12`. This is the `T<=K` temporal-Gram branch, so the factor is copied
directly from the complete `SYEVR` eigensystem rather than formed through division by small singular values. The
observed defect is only `1.164` times the original bound. P4 therefore supplies the explicit solver-specific
acceptance bound `128 epsilon max(T,q)` while mxlib retains its stricter general default. A controlled eigensolver
test admits a defect between those two bounds. The matching local 621-image, 2988-predictor annulus completed all 132
search pixels under the P4 bound, and the next ROC run also completed that annulus. The following `[12,14)` annulus,
however, produced `max|L^T L-I|=4.8173e-11`, or `2.73` times the P4 bound and above even
`256 epsilon max(T,q)`. This second observation demonstrates that repeatedly relaxing a fixed constant would chase a
variable `SYEVR` orthogonality tail. P4 therefore retains the 128-scaled bound and sends any larger defect to the
recorded exact-oracle fallback. Do not reorthogonalize the temporal factor alone because that would make it
inconsistent with the unchanged singular values.

The exact-fallback slice passes 905 assertions in 25 focused `P4PCA` cases against the staged ABI-v2 mxlib build.
Those cases compare complete residual and validity products with `calculateHeldOut()`, cover accepted and rejected
factor defects, preserve the separate rank-boundary reason, and keep non-orthogonality-unrelated validation failures
fatal. The full `P4Reduction` executable passes 56,665 assertions in 20 cases, including a two-worker forced
factor-validation fallback and exact console, FITS-header, and `p4RegionSummary.fits` provenance checks. The
documentation target, formatting check, and whitespace check also pass; the documentation build retains only its
pre-existing warnings.

This is the requested stop-and-review point for P4 integration. Before beginning bounded-rank work or treating the
backend as an optimization, the open gates are an exact-path one-versus-multiple-worker product comparison with inner
BLAS threading controlled, a forced application-level rank-boundary fallback including its warning, forced
memory-limit behavior, production-dimension timing and peak-RSS measurement, and a remote real-data science
comparison. Passing those checks advances the full-rank correctness slice; it does not by itself satisfy the
structured-performance or science-promotion steps below.

The ROC timing baseline in `working/roc/finim_0059_timing.txt` used the exact factor path with the dense
`leadingCovariance` backend over 5,024 search pixels. The P4 algorithm took 3,041.519 elapsed seconds. Its accumulated
worker time was dominated by the existing combined eigensolve/downdate bucket: 132,504.042 worker seconds, or
98.745% of the measured worker total, versus 996.304 seconds for Gram construction and 475.500 seconds for
projection/residual work. This is essentially the same regression elapsed time as the preceding 3,032.06-second
run, so it is the comparison baseline for the structured integration rather than evidence of a dense-path speedup.

### Structured one-row mxlib checkpoint

Step 8 now has a working mxlib-only `rankOneSecular` backend. For one deleted singular-factor row it evaluates the
same normalized covariance core as `leadingCovariance`,

```text
H = diag(s^2) - z z^T,    z = s .* f,
```

by solving the positive rank-one problem `-H = diag(-s^2) + z z^T`. A LAPACK `LAED9` secular solve handles the active
system after a DLAED2-style pass deflates negligible update components and clustered poles. Recorded Givens rotations
embed the active eigenvectors back into the original factor coordinates in reverse order. The backend computes the
complete spectrum for PSD and rank diagnostics, reconstructs only the requested leading rotation columns, reuses all
workspace, and has `O(q^2)` time and storage per deletion after preparation.

The provider call is owned by the output-first `mx::math::laed9<float/double>()` wrappers in mxlib's
`templateLapack` layer rather than by the downdate implementation. Both typed wrappers have direct analytic tests and
complete executable-line coverage on the available OpenBLAS/LAPACK build.

The public contract is deliberately narrow at this checkpoint: an empty deletion returns the identity update,
exactly one row uses the structured solver, and more than one row returns `unsupportedDeletionCount`. It does not
silently select a dense backend. Finite/order/interlacing/PSD, eigenvector norm, and matrix-free residual checks turn
malformed or numerically unsafe solver output into explicit status values. Roundoff-scale negative covariance values
are clamped under the existing result contract; material indefiniteness is rejected.

The optimized focused executable passes 509 assertions in 20 test cases. The complete mxlib coverage target passes
all 178 executables, including the mixed-Eigen-alignment dynamic-symbol fixture, and reports 1,044/1,044 executable
lines and 80/80 functions covered in `source/math/svdDowndate.cpp`. Comparisons include direct row and column SVDs,
leading-projector equivalence, float and double, scalar and zero systems, repeated and clustered spectra, zero and
high leverage, roundoff clamping, overflow, invalid contracts, and injected LAPACK failures. The earlier mixed
32-byte-producer/16-byte-consumer ASan/UBSan run also completed without a diagnostic.

A local P4-shaped synthetic comparison used `q=621`, `outputRank=465`, a square orthogonal factor, and one deletion.
The dense `leadingCovariance` backend took 0.0392046 seconds per deletion and `rankOneSecular` took 0.00936694 seconds,
a 4.19-fold speedup. The relative squared-spectrum error was `5.20346e-15`, the leading-projector error was
`5.15436e-13`, and all 621 structured deletions completed in 7.82617 seconds with no failures. This is a useful local
kernel result, not yet a production P4 timing claim.

This was the requested stop-and-review point before P4 integration. The review accepted the reusable API and direct
SVD comparisons, so the following checkpoint records the deliberately opt-in P4 wiring. The current host verifies
the OpenBLAS/LAPACK ABI; a oneMKL compile/link smoke test remains open because oneMKL is not installed here.

### P4 structured one-row integration checkpoint

P4 now exposes `p4.deletionBackend=leadingCovariance|rankOneSecular` separately from
`p4.exclusionSolver=explicitRefit|factorDowndateExact`. Both defaults remain unchanged: direct explicit refits are the
science default, and `leadingCovariance` is the factor solver's compatibility backend. `rankOneSecular` is accepted
only with `factorDowndateExact`.

The application validates actual compact exclusions, not merely the configured `adi.minDPx`, before allocating
worker state or starting regression. Every active target must delete exactly one row. A target which deletes every
row remains the existing rank-zero invalid result and bypasses mxlib; every other multirow request is rejected up
front. P4 passes the same selected backend to workspace preparation and every independent deletion from the immutable
base factor, so backend selection cannot change partway through a search pixel.

A recoverable numerical `rankOneSecular` failure, including a positive LAPACK convergence status, abandons the
structured result for that search pixel and recomputes all of its target rows once with
`P4PCA::calculateHeldOut()`. It does not mix rows from the two solvers and does not silently invoke
`leadingCovariance`. A negative LAPACK status denotes an illegal argument and remains a hard contract failure, as do
invalid contracts, allocation failures, unsupported deletion counts, dense-backend failures, and failure of the
explicit oracle. FITS provenance records the selected backend and a distinct deletion-solver fallback search-pixel
count alongside the existing rank-boundary and factor-validation counts.

Timing retains the existing eigensolve/rank-selection aggregate and adds reusable-base-factor and row-deletion
subtotals plus a complete explicit-fallback retry timer. The latter overlaps phase-specific oracle counters and is
not additive. `p4Timing.fits` uses schemas 6, 7, and 8 for the base, local-PSF, and reconstructed-PSF layouts. Focused
local verification passes 1,121 assertions in the 28-case `P4PCA` suite and 65,382 assertions in the 20-case
`P4Reduction` suite. The latter forces a positive secular-solver failure through a complete reduction and verifies
the whole-pixel explicit result, warning, counters, and FITS provenance. A remote wall-time comparison remains the
promotion evidence for this opt-in backend; multiple-row structured deletion remains future work.

## Work sequence

### 1. Freeze the math with a standalone oracle

- Implement a small prototype of the symmetric core and the P4 factor-space prediction equation outside the
  production path.
- Compare full-rank results with explicit held-out PCR for deterministic tall, wide, and square matrices.
- Implement the literal Long--Males core alongside it and document the expected difference for truncated factors.
- Include the `X=[1;1]` wide-temporal counterexample and the square-temporal coincidence case so an orientation-limited
  derivation cannot pass.
- Confirm that the symmetric formulation matches the exact oracle in both `T<=K` and `K<T` regimes before designing
  around performance results.

### 2. Add the reusable mxlib deletion system

- Implement the symmetric leading-spectrum core and the complement-preserving small-SVD core with FP64 first, then
  the float instantiation.
- The dense leading backend computes the complete spectrum for PSD/rank diagnosis but publishes only requested
  rotations; the stable backend requests only the preserved-side vectors it needs from a reusable `GESVD` workspace.
- Implement the row/column convenience APIs and make their backend choice explicit.
- Add reusable allocation-free-after-prepare result/workspace objects and deterministic allocation/query/solve error
  seams.
- Add Doxygen equations, approximation/precision contract, references, and examples independent of P4/KLIP.
- Add CMake lists, umbrella include, bibliography, and strict same-stem tests.

### 3. Reach the mxlib correctness and coverage gate

- Obtain 100% executable-line coverage for every mxlib specialization each edited P4 function will call, including
  existing base-factor/eigensolver wrappers as well as the new double API, before integrating it.
- Cover allocation/query/solver failures as well as successful float/double paths.
- Add the counter-bearing implementation records to mxlib's expected inventory, then perform an explicit quantitative
  LCOV audit; presence in `expected_coverage_records.txt` does not itself enforce 100%.
- Run optimized tests plus ASan/UBSan where supported, and verify a clean installed mxlib through an external
  pkg-config consumer.
- Do not route P4 through untested `eigenGESDD`; its current float/double paths have no coverage.
- Put the required immediate Doxygen block before every Catch2 `TEST_CASE` and keep direct production-API calls (or
  Doxygen-only references when a failure seam obscures them).
- If any called mxlib API remains below 100%, add the concrete ownership follow-up to
  `agents/plans/mxlib_cleanup.md` as required before proceeding.

### 4. Replace retained-row lists with compact exclusions

- Introduce `ExclusionPattern` and reproduce every `none`, `imno`, `angle`, and `pixel` selection test.
- Verify inclusive boundaries, angle wrap, variable endpoint cardinality, and no dataset wrap.
- Give the exact oracle a lazy retained-row complement so it remains available without `O(T^2)` persistent indices.
- Update annulus memory estimates and failure-before-regression checks.

### 5. Integrate a full-rank P4 downdate backend

- Add the separate `P4PCA` backend and preserve the direct oracle untouched.
- Retain every safely representable positive singular triplet, which makes deletion algebraically exact for the
  represented design; do not divide through exact/unsafe zeros or use science `rankTolerance` as the base cutoff.
- Reject an unsafe or inconsistent complete base or deletion failure. Use deterministic whole-search-pixel explicit
  fallback only when a post-deletion eigenvalue is numerically indistinguishable from the science rank threshold.
- Compare every residual, rank decision, and validity bit with the direct implementation, using the rank-threshold
  ambiguity policy rather than demanding unstable bit equality at the boundary.
- Cover delete-one, endpoint windows, noncontiguous angle/pixel sets, structural rank loss, and both Gram orientations.
- Preserve the all-rows-excluded rank-zero/invalid result without calling mxlib.
- Verify one and multiple P4 workers produce equivalent products with inner BLAS threads controlled and effective
  process threading recorded.

### 6. Add bounded working rank and convergence checks

- Start at `q=min(rank_struct,m_max+c_max)` and sweep additional guard modes.
- Compare against the direct oracle at the per-target/per-mode residual level, then compare complete combined FITS
  products and companion SNR/merit on representative real data.
- Explicitly measure the remaining approximation that the global rank-`q` subspace was itself selected with the
  held-out target present. Repeat the positive-planet injection/recovery test as well as the real-source reduction.
- Set quantitative residual, target-leakage, injected-source recovery, and science-metric tolerances before considering
  the projected backend for routine use. It never receives exact-exclusion provenance.
- Distinguish truncation error from numerical failure. Increase `q` or explicitly fall back when a configured
  convergence policy fails; never hide a plane-meaning change.
- Decide whether guard rank needs public configuration only from these results.

### 7. Benchmark the dense backend at production dimensions

- Measure the existing exact baseline, full-rank symmetric downdate, and truncated symmetric downdate separately.
- Include the approximately 621-frame/high-`0.75`-mode AF Lep case and the greater-than-8,000-frame case with
  `K` around 4,000.
- Report base factorization, core assembly, reduced solve, projection, wall time, worker seconds, peak RSS, selected
  worker count, and output differences.
- Sweep outer worker count with inner BLAS/LAPACK fixed to one. A separate experiment may compare fewer outer workers
  with a multithreaded base/core solve.

Core assembly costs `O(c_t q^2)`, the dense per-target solve costs `O(q^3)`, and factor-space projection costs at
most `O(q^2)` when all leading rotations are returned. Storage is `O(Tq+q^2)` in addition to `X`. Near-linear scaling
in `T` assumes fixed `q`; P4's fractional modes can make `q` grow with `min(T,K)`. At `T=621` and a largest fraction
of `0.75`, `q` is about 466 before extra guards, so a dense downdate offers only a limited cubic reduction relative to
dimension 620. At `q` near 3,000, dense per-target solves may still be infeasible. These are acceptance measurements,
not details to defer until after making the backend default.

### 8. Add a structured deletion backend if dense solves remain limiting

- Use the diagonal-minus-rank-`c` form `H=diag(s^2)-Z Z^T` rather than treating `H` as an arbitrary dense matrix.
- The delete-one Gu--Eisenstat/secular eigendowndate with `O(q^2)` work is implemented in mxlib and wired into P4 as
  the explicit `rankOneSecular` backend described above.
- For small blocks, compare a stable block method with `c` within-target rank-one steps. Each target still starts from
  the immutable base spectrum; this does not introduce cross-target chaining.
- Keep the dense backend as a development oracle and as a fallback for supported, adequately conditioned block cases.
  An ill-conditioned/high-leverage structured failure falls back to the stable complement core or the direct retained
  refit, not automatically to cancellation-prone dense `H_E`.
- Revisit CPU batching and then GPU batching only after the structured CPU workload is known. The earlier GPU result
  rejected repeated full-size direct eigensolves; it did not benchmark this small structured kernel.

### 9. Promote only after science and operations acceptance

- Allow an `auto` choice to select the exact downdate only after full-rank equivalence, remote performance, memory
  accounting, and real-data science checks pass. It must not select the projected backend until the separately
  measured truncated-rank convergence and leakage results support an explicit science decision.
- Document the approximation boundary, guard rank, backend provenance, failure/fallback behavior, and interpretation
  of capped rank.
- Retain the explicit path for tests, diagnostics, and cases where `q` approaches full rank and the downdate has no
  useful performance advantage.

## Verification matrix

### mxlib

- Tall, wide, and square matrices.
- Row and column deletion; first, middle, last, contiguous block, and arbitrary noncontiguous indices.
- Empty deletion as identity if supported; duplicate/out-of-range indices and removal of the complete side rejected.
- Full-rank factors compared with a direct SVD of the physically retained matrix.
- Truncated factors compared with a direct SVD of the physically retained rank-`q` reconstruction.
- Exact rank loss, zero/nearly dependent rows, high-leverage deletions, zero singular values, and materially non-PSD
  cores.
- Repeated and clustered singular values, comparing spectra, reconstructions, projectors, and predictions rather than
  raw vector signs.
- Descending ordering, selected output rank, roundoff clamping, workspace reuse/growth/ownership, independent
  repeated calls, and injected LAPACK/allocation failures.
- Float and double API paths; optimized and sanitizer builds.

### hciReduce

- `T<=K`, `K<T`, and square designs.
- Delete-one and variable endpoint windows; `imno`, `angle`, and `pixel` exclusion sets.
- Full-rank downdate versus explicit oracle at tight FP64 tolerance.
- Truncated rank versus explicit rank-`q` projected oracle, plus convergence toward the complete refit as `q` grows.
- Excluded-target leakage: changing `y[E_t]` leaves the fitted prediction unchanged; changing `y[t]` changes only the
  final observed-minus-predicted residual. For the complete-rank backend, changing excluded predictor rows other than
  the held-out application row also leaves the trained model unchanged. Measure the corresponding leakage of the
  projected backend rather than assuming it is zero.
- All modes supported, partial numerical support, strict threshold equality, and structural rank loss.
- Exact residual, `sampleValidity`, per-mode status, and capped/minimum-rank behavior.
- Compact versus retained output, one versus multiple OpenMP workers, memory-limited worker selection, timing schema,
  and FITS provenance.
- Rejection of rotated frame, positive `numberImages`, local processing, and frozen-PSF products in the first slice.
- Complete AF Lep final products, science SNR/merit, wall time, and peak memory against the accepted exact baseline.

## Acceptance criteria

1. Full-rank factor deletion reproduces explicit held-out P4 predictions and validity in both matrix-shape regimes;
   degenerate subspaces are compared by invariant outputs rather than vector signs, and threshold-adjacent ranks use
   the explicit ambiguity fallback.
2. The mxlib double API called by P4 has 100% executable-line coverage and no copyable raw-workspace ownership.
3. Truncated working-rank results have a measured convergence story against explicit refits; `q=k+c` is not accepted
   solely because it is the published KLIP heuristic. Their residual/leakage/injection/science tolerances are
   quantitative, and they are never labeled exact.
4. The selected backend materially reduces remote wall time for the actual high-mode configurations, not only for a
   fixed-small-`q` synthetic benchmark.
5. Peak/per-worker memory is included in P4's automatic worker budget, and compact exclusions eliminate the retained
   `O(T^2)` index lists.
6. Solver, base rank, full/truncated state, rank/clamp/fallback diagnostics, and timing are recorded without changing
   science-plane meanings.
7. `adi.excludeMethod=none` and every unsupported held-out combination retain their current behavior.

## References

- Long and Males (2021), [Unlocking starlight subtraction in full data rate exoplanet imaging by efficiently updating
  Karhunen--Loeve eigenimages](https://arxiv.org/abs/2101.11634), especially Sections 4--5.
- Long and Males, [`minimal_downdate` reference implementation](https://github.com/joseph-long/svd-downdate-klip/blob/main/starbgone.py).
- Brand (2006), [Fast low-rank modifications of the thin singular value
  decomposition](https://doi.org/10.1016/j.laa.2005.07.021).
- Bunch and Nielsen (1978), [Updating the singular value decomposition](https://doi.org/10.1007/BF01397471).
- Gu and Eisenstat (1995), [Downdating the Singular Value
  Decomposition](https://doi.org/10.1137/S0895479893251472).

# Reviews

## After first mxlib implementation milestone

- See #1 in mxlib/AGENTS.md, where I have changed the file block slightly and changed the policy to exclude the \author tag

- In svdDowndate.cpp it looks like there are internal functions with doxygen docs.  I prefer to document those in the same style but non-doxygen, so, e.g. `///` should be `//`.  Update AGENTS.md with this preference.

- See #12 in mxlib/AGENTS.md, where I have updated the unit test doc guidelines

- svdDowndate_test.cpp should be in a \ingroup group for tests, but I'm not sure if we have established a standard or are following it.  I think we want a hierarchy math->math_files->math_test_files.

- How does Catch2 report REQUIRE failures in a wrapper function?  Does it report the line number of the REQUIRE but also the test case where it failed in a usable fashion?

- In relation to the above, I think I would prefer a style where the actual assertion is local to the TEST_CASE just for readability.  Using the wrapper makes sense for the same readability reason, but we could assert on its return value(s) in a transparent way.  If there is a solution here, update AGENTS.md

- svdDowndate_test.cpp same comment about doxygen in internal docs.  Though in this case there may be some merit in the test documentation being discoverable if there is a well structured grouping system to keep it out of the main API docs.

- For "SVD Row and Column Deletion" dox group, at the group level we should have expository text as well as example code for how to use the SVD downdate.

- Let's put the SVD docs and the Templatized BLAS and Lapack under a new "Linear Algebra" group
