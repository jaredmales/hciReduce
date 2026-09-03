# P4 and KLIP FP32 Evaluation and Decision Record

## Goal

Determine, independently for P4 and KLIP, whether moving principal-component arithmetic to IEEE single precision
provides a material end-to-end wall-time or memory benefit without unacceptable changes to numerical decisions,
validity, or science products.

The evaluation must distinguish three questions:

1. Does FP32 make the isolated Gram, eigensolve, and projection kernels faster on the CPUs used for reductions?
2. Does that improvement survive sampling, selection, threading, finalization, and I/O in a real reduction?
3. Are FP32 rank decisions and science products acceptable, including ill-conditioned and threshold-adjacent cases?

P4 and KLIP were evaluated independently. The accepted production policy is M32D64 for both algorithms: FP32 PCA
calculation with an FP64 direct eigensolve. The all-double P4 implementation remains a regression oracle, and a fully
FP32 eigensolve is deferred to future GPU-acceleration work.

## Status and decision framing

Planning began 2026-09-02 under [HCI-015](issues.md#hci-015--evaluate-fp32-arithmetic-for-p4-and-klip). On 2026-09-03,
the owner accepted `P4-M32D64` as the P4 production default and retained `KLIP-M32D64` as the KLIP production default.
The decision prioritizes the measured FP32 calculation/Gram speedup while retaining the lower-risk, established FP64
eigensolve. P4's exact factor-deletion path and its whole-search-pixel explicit fallback remain FP64. KLIP's exact
factor-deletion path likewise remains FP64 internally before its modes return to FP32 storage.

The owner acknowledged the thermally throttled timing host and the broader incomplete frozen gate set when making this
decision. The archived acceptance and run manifests remain unchanged historical evidence: their
`promotion_ready=false` values and report-only rows are not retroactively converted into formal gate passes. They
support the engineering judgment but do not claim that every originally proposed platform, path, injected-source, or
science threshold was completed. Fully FP32 eigensolves, FP32 factor deletion, guarded/runtime precision policies,
and additional CPU promotion work are tabled; reconsider them with future GPU acceleration.

Initial execution uses the owner-designated known-good files `working/kr.conf` and
`working/p4Reduce_afLepNaco.conf`. Their hashes, the complete input identity, controlled smoke profiles, environment,
and ignored artifact locations are frozen in
[`benchmarks/manifests/p4_klip_fp32_aflep.yaml`](../../benchmarks/manifests/p4_klip_fp32_aflep.yaml). The source files
themselves remain unchanged. The active blank `rankTolerance=` in the P4 file resolves to zero in the current parser, as confirmed
by the output FITS provenance; zero is therefore the configuration-faithful baseline. Runs at `1e-12` or a future
FP32-aware threshold are separate rank-policy experiments, not silent substitutions.

The key comparison was not simply `float` versus `double`. P4 needed a KLIP-style mixed candidate so the value of an
FP64 eigensolve could be separated from the cost of performing all surrounding dense arithmetic in FP64. KLIP already
provided that mixed baseline. The results supported the fixed mixed policy; no guarded runtime fallback was needed for
the accepted direct/refit scope.

## Current precision map

The precision inventory must classify science arithmetic separately from geometry, timing, configuration, counters,
and provenance. A `double` used for an angle, elapsed time, byte calculation, or FITS header is not part of the hot PCA
path and must not be changed merely to make the code textually all-float.

| Stage | P4 production path | KLIP production path | Decision consequence |
|---|---|---|---|
| Loaded/preprocessed images | FP32 | FP32 | Already common; preserve the existing observation lifecycle. |
| Spatial sampling | FP32 cubic/direct sampling | FP32 region extraction | Measure separately from PCA; do not change interpolation or geometry. |
| PCA input boundary | Sampled values retain their checked FP64 ingress, then the production PCA adapter converts them to FP32 worker scratch | Regional arrays remain FP32 | Preserve checked sampling and existing public P4 data boundaries while using the accepted calculation scalar internally. |
| Gram/covariance construction | FP32 Eigen products | FP32 `eigenSYRK` or Eigen products | This is the principal accepted arithmetic change. |
| Direct eigensolve | FP32 covariance promoted into an FP64 `SYEVR` workspace; successful eigenpairs return to FP32 before rank selection | FP32 covariance copied to an FP64 workspace for `SYEVR`; post-solve precision differs by Gram branch | Retain the established FP64 solver boundary in both algorithms. |
| Mode projection/residual | FP32 calculation, then checked conversion to the existing result boundary | FP32 | Preserve result layout, validity, and product contracts. |
| P4 centering, explicit held-out fits, coefficients, and frozen probes | FP32 calculation with FP64 eigensolve | Not applicable | Apply the accepted mixed policy consistently to normal, rotated-reference, local, and shared-PSF direct/refit paths. |
| Exact factor deletion | FP64 base factors, mxlib workspaces, updates, and whole-search-pixel fallback oracle | FP64 base factors and mxlib deletion path; modes cast to FP32 | Keep this opt-in path as an explicit FP64 exception. |
| Final residual products | FP32 | FP32 | Stored equality alone is insufficient because ranks, validity, and model responses can change before this boundary. |

The quoted 72% Gram/12% eigensolve split comes from an ordinary smaller-Gram P4 profile and is not universal. Recent
exact factor-deletion runs are instead dominated by eigensolve/deletion work. Each solver/path combination therefore
needs its own stage profile; neither historical result can stand in for the requested comparison.

KLIP's two adaptive direct branches are also not numerically identical at the solver boundary. The reference-Gram
branch normalizes and validates eigenvectors in FP64 before constructing FP32 KL modes. The pixel-Gram branch casts
eigenpairs back to FP32 first and performs its mode-validity checks there. Preserve and report those branch-specific
contracts rather than treating `KLIP-M32D64` as one opaque operation.

Relevant implementation boundaries are `P4PCA`, `P4Reduction::fitDetectorSearch()`, `P4Reduction::estimatedWorkerBytes()`,
`calcKLModesAdaptive()`, and `KLIPreduction::worker()`. The current P4 rank rule uses
`lambda > p4.rankTolerance*lambdaMax` with a default relative tolerance of `1e-12`; KLIP accepts positive finite
eigenvalues. Both policies need special scrutiny because `1e-12` is below FP32 machine precision and roundoff can move
small eigenvalues across zero or a requested-mode boundary.

The installed mxlib headers already expose `syevrMem<float>`, FP32 `SYEVR`, and float SVD-deletion entry points. That
establishes source/API availability, not numerical suitability, optimized-library dispatch, performance, or adequate
coverage. Those are explicit capability gates below.

The unversioned P4 prototype used FP32 for its Gram, cross-product, and CPU/GPU solve. That makes full FP32 a credible
experiment, but its benchmark provenance is not recoverable and it is not evidence for current P4 behavior.

The existing `p4DecompositionBenchmark` is an FP64 CPU/CUDA eigensolver benchmark. It excludes most transfers and
setup from GPU timing, validates only one matrix without pass/fail tolerances, and does not exercise Gram construction,
projection, residual formation, or application behavior. It may supply reusable option parsing and reporting, but it
cannot answer HCI-015 by itself.

## Candidate precision configurations

The experiments used compile-time/internal policies or benchmark-only instantiations. The accepted M32D64 policy is a
fixed implementation choice, not a user-facing configuration value.

| ID | Stored PCA operands and Gram | Eigensolve | Projection/residual | Purpose |
|---|---|---|---|---|
| `P4-D64` | FP64 | FP64 | FP64, then existing FP32 output cast | Retained regression and factor-deletion oracle. |
| `P4-M32D64` | FP32 | Covariance promoted only inside an FP64 solve; eigenpairs return to FP32 | FP32 | Accepted P4 default for direct, centered, explicit-held-out, coefficient, and frozen-probe calculation. |
| `P4-F32` | FP32 | FP32 `SSYEVR` | FP32 | Tabled; reconsider with future GPU acceleration. |
| `KLIP-M32D64` | FP32 | FP64 `DSYEVR` | FP32 | Retained KLIP direct-path default. |
| `KLIP-F32` | FP32 | FP32 `SSYEVR` | FP32 | Tabled; reconsider with future GPU acceleration. |
| `P4-F32-delete` / `KLIP-F32-delete` | FP32 operands, base factors, and deletion workspaces | FP32 direct and deletion solves | FP32 | Tabled with fully FP32 eigensolve work. |
| `*-F32-guarded` | Explicitly documented FP32/FP64 split | FP64 only at a defined ambiguity/failure boundary | Explicitly documented | Tabled; the accepted fixed M32D64 policy needs no runtime selector. |

The same input values and operation ordering must be used wherever the scalar types permit. A candidate may not gain
speed by changing mode counts, rank tolerance, Gram orientation, solver algorithm, reference selection, exclusion,
centering, or threading. If an FP32-aware rank floor or ambiguity fallback is needed, evaluate it as a separately
named follow-on policy against the unmodified FP32 policy; do not conceal it inside `P4-F32` or `KLIP-F32`.

For P4, run the direct in-sample and explicit-held-out paths before adding FP32 factor deletion. For KLIP, run the
shared and per-target direct paths before adding FP32 factor deletion. If deletion enters the promotion scope, first
evaluate it with float factors, workspaces, and the matching FP32 direct fallback. KLIP's direct fallback follows its
eigensolver template type; it is not FP64 in a pure `KLIP-F32` build. Any FP64 fallback is a separately named guarded
hybrid. A hybrid that silently promotes the whole fit on every target is not an all-FP32 result. Optional deletion
work must not block a conclusion about the direct paths.

## Constraints and non-goals

- Preserve preprocessing, masks, interpolation, region ownership, exclusion/selection, centering, mode ordering,
  residual finalization, and FITS product semantics throughout the comparison.
- Keep `P4-D64` and `KLIP-M32D64` callable in the experimental build so every candidate can be compared on identical
  in-memory inputs. Do not compare outputs from separately changed algorithms.
- Keep timing values, timestamps, geometry coordinates where double is part of an existing accuracy contract, byte
  accounting, optimizer parameters, and provenance types out of the blanket conversion.
- Do not compare raw eigenvector signs. Compare invariant subspaces/projectors, predictions, and residuals; a repeated
  or clustered eigenvalue admits multiple valid bases.
- Do not add CUDA, TF32, FP16/BF16, tensor-core kernels, fast-math changes, a different eigensolver, a new rank model,
  or general memory-management work to this experiment. CUDA precision is a later, separately controlled study.
- Do not change outer OpenMP or inner BLAS/LAPACK policy to make a precision candidate look faster. Fixed-concurrency
  measurements isolate arithmetic. For P4, a second automatic-memory-policy run measures operational gains from
  smaller workspaces. KLIP has no automatic worker-memory policy today; adding one remains owned by
  `klip_memory_management_optimization.md` rather than this evaluation.
- Do not treat lower RSS as sufficient for promotion unless it enables a previously infeasible reduction or produces
  a measured concurrency/wall-time benefit.
- Preserve the independently reached disposition for each algorithm. Their common M32D64 result does not imply that
  future GPU or lower-precision work must promote both together.

## Measurement protocol

Create a checked-in benchmark manifest and keep large/raw outputs under an ignored `working/fp32/` directory. Every
reported result must identify:

- hciReduce and mxlib commit hashes and whether installed headers/libraries match them;
- compiler and version, Release flags, Eigen version and vectorization macros, BLAS/LAPACK vendor and version, and
  linked libraries;
- CPU model, ISA, socket/NUMA topology, available memory, and affinity/pinning policy;
- outer OpenMP and inner BLAS/LAPACK thread counts;
- dataset/configuration identifiers, file checksums, dimensions, masks, region geometry, `T`, `K`, Gram orientation,
  requested modes, realized ranks, and exclusion/selection backend; and
- warm-up count, paired run order, repetition count, elapsed/worker timing, and peak RSS measurement method.

Alternate candidate order within one process for kernel benchmarks where workspace reuse is part of production
behavior. Use at least five measured paired repetitions for bounded runs and at least three for expensive full-data
runs. Report every sample plus median, range, and a robust spread; do not retain only the fastest run. Record CPU
frequency/thermal anomalies when available.

Measure these concurrency regimes:

1. for both algorithms, fixed equal worker counts with inner BLAS/LAPACK held at one thread, isolating scalar
   precision; and
2. for P4 only, the existing automatic worker policy after its byte estimator is made precision-aware, measuring any
   additional usable concurrency and actual peak RSS.

Production `ReductionTiming` already separates the important P4 phases and most KLIP phases. Add benchmark-local or
scoped KLIP Gram timing if needed rather than changing the persisted production timing schema during the initial
experiment. KLIP's current direct-fallback interval encloses its eigensolve and mode-construction intervals, so those
fields overlap and must not be summed as exclusive stage totals. Use exclusive benchmark counters or label fallback as
a nested diagnostic. Always report whole algorithm and whole application wall time alongside stage timings.

## Numerical and scientific comparison protocol

Before examining final candidate results, check in a tolerance table with units and rationale. It must distinguish
ordinary well-separated spectra from explicitly generated ambiguity-band cases, and it must include both internal
FP64 comparisons and final stored-FP32 comparisons. Tolerances may be scale-aware; a single global `Approx` margin is
not adequate.

For numerical kernels, record:

- Gram/covariance maximum absolute and relative Frobenius error against the FP64 construction, comparing KLIP's
  authoritative lower triangle or symmetrizing a copy outside the timed region before a full-matrix norm;
- eigenvalue error, normalized eigensolver residual `||GQ-QLambda||`, and orthogonality defect `||Q^TQ-I||`;
- projector distance or principal angles for each retained subspace, especially across eigenvalue clusters;
- prediction and residual maximum error, relative L2 error, RMS error, high-percentile error, and error normalized by
  the FP64 residual RMS;
- predictor coefficients and frozen-probe responses where those are production outputs of the computation;
- numerical/base rank, mode status, per-sample validity, clamp/fallback reason and count, and nonfinite handling; and
- determinism across repeated executions and one-versus-many outer workers.

For synthetic tests, use an independent FP64 direct SVD/QR oracle in addition to the current normal-equation path.
Cover both `T <= K` and `T > K`, including:

- well-conditioned full-rank matrices;
- exact and near rank deficiency;
- repeated and tightly clustered eigenvalues, including clusters crossing a requested-mode boundary;
- eigenvalues below, equal to, just above, and safely above `p4.rankTolerance*lambdaMax`;
- large dynamic range, centering cancellation, and FP32 overflow/underflow stress;
- realistic correlated noise and seeded production-sized Gram matrices; and
- held-out deletion patterns at endpoints, contiguous and noncontiguous multi-row sets, and complete structural rank
  loss.

Define a dimension- and precision-scaled rank-ambiguity band. FP32 and the authoritative baseline must agree outside
that band. Inside it, report the affected fraction of real fits and the resulting product sensitivity; either retain
FP64 for that decision, introduce an explicitly reviewed precision-aware rank policy, or reject the FP32 candidate.
Silent validity or plane-meaning changes are not ordinary floating-point error.

For complete products, compare:

- exact output shape, plane ordering, requested/realized mode metadata, selection flags, and finite/validity masks;
- per-plane maximum/RMS/high-percentile differences, correlation, fitted slope, and difference-to-baseline RMS;
- injected-source throughput and recovered contrast, aperture flux, centroid/peak astrometry, SNR, and detection rank;
- for P4, compact-versus-retained products, local evaluations, optimizer separation/position angle/contrast/merit,
  analytic PSF responses, reconstructed PSFs, and filtered products; and
- for KLIP, ADI and RDI products, shared and target-specific bases, selected/excluded libraries, right-reason products,
  and direct/factor-deletion fallback cases.

## Expected implementation surface

The precise design follows the capability spike, but the likely experimental touchpoints are:

- `src/common/P4PCA.hpp` and `src/common/P4PCA.cpp` for calculation/eigensolver scalar policies, precision-scaled
  numerical guards, every direct/centered/held-out/deletion path, results, and workspaces;
- `src/common/P4Reduction.hpp` and `src/common/P4Reduction.cpp` for sampled-value boundaries, worker storage, residual
  conversion, and precision-aware memory estimates;
- `src/common/P4PSFModel.hpp` and `src/common/P4PSFModel.cpp` where float templates currently feed double
  coefficient/probe matrices;
- `src/common/KLIPreduction.hpp`, `src/common/KLIPreduction.cpp`, and `src/apps/klipReduce.cpp` for the existing
  eigensolver template seam, explicit instantiations, and the separately hard-coded FP64 factor-deletion path;
- `tests/common/P4PCA_test.cpp`, `tests/common/P4Reduction_test.cpp`, `tests/common/KLIPreduction_test.cpp`, and
  application tests for FP64-oracle differential coverage and precision-specific bounds; and
- `benchmarks/CMakeLists.txt` plus a CPU precision benchmark available without CUDA. Production configuration and user
  documentation are touched only if a candidate is later accepted.

Use a compatibility-first P4 scalar seam. The first production-source patch adds a Doxygen-hidden, anonymous-namespace
`P4PCAKernel<CalculationT, EigensolverT>` for only the direct and centered paths, instantiates only
`<double,double>` in a normal build, and makes the existing concrete `P4PCA` methods forward to it. `P4PCA`,
`P4PCAResult`, normal-build exported symbols, the double failure-injection seam, held-out/probe paths, factor deletion,
`P4Reduction`, and `P4PSFModel` remain unchanged in that extraction. Its acceptance criterion is unchanged FP64
behavior.

After that extraction passes, explicitly instantiate `<float,double>` and `<float,float>` for the direct/centered
experiment while retaining double result storage. The mixed helper must preserve LAPACK status: copy the FP32 Gram
matrix into the FP64 workspace, call `eigenSYEVR` directly, and cast the selected eigenpairs back. Do not route this
through `calcEigenVecs`, which collapses a nonzero LAPACK status to `-1`. The explicit held-out/probe paths now use the
same typed core. Typed factor deletion and the shared-PSF coefficient path remain later stages because their tolerance
and fallback policies require separate design. The first reduction-level experiment uses only macro-gated internal
dispatch chosen outside OpenMP worker loops and does not change the existing `P4Reductionf` specialization ABI or
application configuration. It now routes detector-frame direct, rotated centered-in-place, explicit-held-out, and
target-held-out frozen-probe work through caller-owned per-worker typed workspaces. Non-double dispatch deliberately
rejects automatic memory limiting, exact factor deletion, shared-PSF coefficients, and local trials before reduction
state can be mutated.

## Initial execution record (2026-09-02)

The initial Release baseline passed `P4PCA` (2,087 assertions/30 cases), `P4Reduction` (75,452/22), `p4Reduce`
(571/14), `KLIPreduction` after adding the two FP32 capability cases (9,441,428/44), and `klipReduce` (36/7), with one
outer worker and one OpenBLAS thread for the focused reruns.

The two configs address the same contiguous 621-file, FP32, 256-by-256 AF Lep/NACO sequence. Its ordered basename
list hash is `498662d5d03eadb2492880fb74d8bf573cea58d6743e6832317d16e2545186ca`; the SHA-256 of the ordered
`sha256  basename` stream is `0d73f1c13877a225e17a69b78ae6b4ee219e782916421d24b7b8a156d4d6bcd7`. KLIP crops
these inputs to 128 by 128 for the configured `[6,60)` region. The repository manifest records the config and mask
hashes, input endpoints, dependency versions, exact overrides, output hashes, and concurrency hazards.

Four single-run smoke profiles established working baselines without supporting a speed claim:

| Case | Bound and numerical shape | Algorithm profile | Whole process |
|---|---|---|---|
| `P4-D64` | 96 frames, `[10,12)`, 132 fits, `K=2,988`, temporal Gram | 0.846 s: sampling 0.366 worker-s, Gram 0.384, eigensolve 0.081, projection 0.0009 | 1.72 s, 326,924 KiB peak RSS |
| `P4-D64` full geometry | 96 frames, all 19 annuli, 5,024 fits, `K=1,776`--`4,108`, temporal Gram, four workers | 12.021 s: sampling 20.969 worker-s, Gram 21.726, eigensolve 4.054, projection 0.031 | 14.18 s, 2,006,224 KiB peak RSS |
| `KLIP-M32D64` shared | 128 frames, `[6,60)`, 11,180 pixels, reference Gram | 0.453 s: eigensolve 0.0025 worker-s, mode construction 0.0022, subtraction 0.443 | 0.79 s, 272,276 KiB peak RSS |
| `KLIP-M32D64` target-specific | 32 frames, delete-self fits, `[6,60)`, reference Gram | 0.017 s: eigensolve 0.0039 worker-s, mode construction 0.0039, subtraction 0.0063 | 0.13 s, 266,264 KiB peak RSS |

The first capability test now runs `calcKLModesAdaptive<float>` and the current double-eigensolve path on identical
FP32 inputs in both Gram orientations. It compares individual-mode projectors, retained-subspace projectors, and
target residuals against one another and an independent FP64 direct-SVD oracle; its 1,852 assertions pass. A second
790-assertion test instantiates the actual FP32 `KLIPreduction::worker()` specialization and checks shared
reference-Gram and target-specific selected pixel-Gram residuals, finite output, mode clamping, and inclusion flags
against both the mixed worker and the FP64 SVD oracle. Binary inspection confirms that the installed
`mx::math::syevr<float>` wrapper calls `ssyevr_`, and the linked OpenBLAS exports both `ssyevr_` and `dsyevr_`. The
current 2026-08-30 mxlib LCOV report fully covers the corresponding double adapter paths and the lower-level float
LAPACK wrapper, but it does not emit the public `eigenSYRK<float>`/`eigenSYEVR<float>` adapter and float workspace
instantiations; the three executable lines in each of `syrk<float>` and `gemm<float>` are also at zero hits. The
hciReduce cases are valid linked integration exercises, but cannot satisfy mxlib's own coverage requirement. The
exact gap and required focused mxlib tests are recorded under `Known non-blocking ownership follow-ups` in
`mxlib_cleanup.md`.

The initial CPU-only precision-microbenchmark rows exercise the well-conditioned, uncentered, in-sample direct
microkernel. They hold the FP32 source values and operation order fixed, rotate the three candidates in balanced
triplets, and report conversion, Gram/cross-product, eigensolve, projection/residual, and total timing as raw CSV.
Separate pre-timing validation checks Gram, retained eigenvalues/projectors, residuals, solver residual, and
orthogonality. These are engineering smoke gates; this initial run did not cover numerical-rank selection, centering,
held-out fits, coefficients/probes, output materialization, or an independent accuracy oracle.

| Synthetic shape | `P4-D64` median | `P4-M32D64` median | `P4-F32` median | Initial observation |
|---|---:|---:|---:|---|
| `T=96`, `K=384`, 32 modes | 1.064 ms | 0.839 ms (1.269x) | 0.735 ms (1.448x) | Both FP32 Gram candidates pass the smoke gates. |
| `T=384`, `K=96`, 32 modes | 1.186 ms | 0.887 ms (1.337x) | 0.782 ms (1.516x) | The predictor-Gram branch shows the same ordering. |
| `T=96`, `K=2,988`, 48 modes | 3.371 ms | 1.998 ms (1.687x) | 1.943 ms (1.735x) | This matches the bounded AF Lep annulus shape. |
| `T=621`, `K=2,988`, 466 modes | 152.426 ms | 97.333 ms (1.566x) | 100.937 ms (1.510x) | FP32 Gram is about twice as fast, but `SSYEVR` is slightly slower than promoted `DSYEVR` here. |

At the 621-by-2,988 shape, the full-FP32 retained-projector relative error is `2.31e-5`, the residual relative error is
`4.91e-5`, and the orthogonality defect is `4.94e-4`; all pass only the declared smoke gates. This result makes the
mixed candidate especially important and justifies the next rank-sensitive experiment, but it is not a production
precision decision. Exact raw rows, file hashes, stage medians, and gates are in the manifest.

Exact baseline products and config copies are retained below `working/fp32/baseline/4869c2a/`, which is ignored by
Git. KLIP printed managed-environment stream-file-descriptor warnings while completing successfully; final timing
must be repeated outside that environment. KLIP cases must not run concurrently because current mask-cube generation
shares `/tmp/hciReduceAux` sidecars.

## Continued engineering record (2026-09-02)

The all-double direct/centered scalar seam is now extracted inside `P4PCA.cpp` without exposing the kernel template.
The full focused P4PCA suite passes 2,100 assertions in 30 cases, including the in-place centered failure seam. The
bounded AF Lep compatibility product is byte-identical to the pre-extraction product
(`9dff0820fe2f8390790a0c44499682b37f20f899c1d386013864cff439fd4774`), and the normal shared library retains the
same complete 664 dynamic symbol names and the same 11 P4PCA-related names. This establishes DD compatibility; it is
not evidence for a float candidate.

A default-off `HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION` build now instantiates that same anonymous kernel as
`<float,double>` and `<float,float>` and exposes only Doxygen-hidden experimental detail wrappers. The mixed adapter
promotes the FP32 Gram into the reusable FP64 workspace, preserves the native LAPACK status exactly, and casts
successful eigenpairs back before rank selection; bulk centering, projection, coefficients, residual arithmetic,
explicit target-held-out fits, and frozen-probe responses use the calculation scalar. With the option on, 4,564
assertions in 39 P4PCA cases pass, including both Gram orientations, direct and centered/in-place results,
coefficients, explicit held-out/probe validity and NaNs, noncontiguous requested modes against independent retained-row
SVD fits, a positive rank boundary, alternating workspace reuse, FP32-ingress overflow, alias rejection, invalid
policies, and native `37`/`-1000` failure propagation. ASan/UBSan passes the same suite. With the option off, 2,107
assertions in 30 cases pass, the complete 664-name normal
dynamic-symbol set is unchanged, and a second bounded AF Lep product after the held-out refactor is byte-identical.
At this checkpoint the typed code remained an internal evaluation surface: no centered-held-out production API
existed, factor deletion remained all-double, and the production worker-storage, memory-policy, configuration, and
default paths were unchanged. The option-on build was not a stable external ABI: its detail wrappers and FP32
dependency instantiations were visible dynamic symbols.

The same compile-time gate now reaches the reduction worker without adding a configuration value. An immutable
DD/M32D64/F32 function table is selected before each OpenMP loop, and every active worker receives a private typed
workspace. Focused synthetic reductions cover detector-frame direct fits, rotated centered-in-place fits, explicit
held-out rank-invalid validity/NaN propagation, and the target-held-out frozen-probe/PSF batch. DD is exact against the
ordinary production call, both float candidates preserve expected structure and declared numerical bounds, and F32
is exact between one and four outer workers. The option-on `P4Reduction` suite passes 294,152 assertions in 26 cases
(218,700 assertions in four experimental cases); the option-off suite passes 75,452 assertions in 22 cases. A fresh
option-off AF Lep bounded product after this dispatch integration is again byte-identical to the baseline, and the
original 664-name demangled dynamic-symbol set remains exact. These are capability and compatibility results, not
contract-bound real-product evidence for either float candidate.

The P4 microbenchmark now includes an independent FP64 thin-SVD structural oracle and two deliberately
rank-sensitive families in both Gram orientations. At synthetic `tau=1e-3`, the oracle and all candidates report
rank 58; one mixed and three all-FP32 eigenvalues fall inside their construction-specific ambiguity bands. At
`tau=1.1e-3`, the oracle and `P4-D64` report rank 59 while both FP32-Gram candidates report rank 58; those raw
differences are explicitly reported and confined to the diagnostic bands. At the configuration-faithful `tau=0`, an
exact-nullity-four construction has structural rank 60, while the raw solver ranks range from 61 to 63 depending on
orientation and candidate. All robust positive modes agree outside the bands, but the extra near-zero modes identified
a residual rank-policy risk that the owner later accepted for M32D64 and used as a reason to table F32 eigensolve.
These bands are bounded-construction diagnostics, not general forward-error or science tolerances.

The KLIP arithmetic benchmark compares the current `KLIP-M32D64` direct path with `KLIP-F32`, using warmed reusable
buffers, an independent FP64 thin-SVD oracle, and direct production-API cross-checks. Runs were pinned to CPU 4, with
one OpenBLAS and one OpenMP thread, two warmups, and six balanced repetitions:

| Scenario | `KLIP-M32D64` median | `KLIP-F32` median | Stage interpretation |
|---|---:|---:|---|
| Shared reference basis, `P=11,180`, `R=T=128`, 64 modes | 338.647 ms | 333.952 ms (1.014x) | `DSYEVR` to `SSYEVR` is 2.144x faster, but subtraction is about 98% of this arithmetic scope. |
| Per-target pixel Gram, `P=96`, `R=384`, `T=32`, 32 modes | 32.516 ms | 19.810 ms (1.641x) | The eigensolve is 1.901x faster and dominates this deliberately selected regime. |

Both candidates' observed metrics fall below the current engineering-smoke limits in both scenarios; staged results
match the production mapping exactly in the tested scope. These rows predate the frozen acceptance contract and lack
its required identity/gate columns, so they remain non-binding exploration rather than formal gate passes. The shared
case peaks at 323,756 KiB and confirms that validation no longer materializes a `P`-by-`P` projector. These results
make the target-pixel branch promising, while showing that an eigensolve speedup alone is nearly irrelevant for the
shared AF Lep shape. Selection, collapse, OpenMP scheduling, multi-mode product materialization, rank-boundary
spectra, and science-product gates remain outside this benchmark. Exact rows, hashes, medians, ranges, MADs, and scope
are recorded in the run manifest.

The acceptance contract is kept separately in
[`benchmarks/manifests/p4_klip_fp32_acceptance.yaml`](../../benchmarks/manifests/p4_klip_fp32_acceptance.yaml). Its
exact structural, engineering-smoke, stored-FP32, and synthetic ambiguity rules may be frozen during kernel work;
science and real-fit ambiguity limits were explicitly owner-TBD and therefore blocked a formal manifest promotion.
The later owner decision did not change those frozen results.

The emitters now also produce the uniform 94-column `hcireduce-pca-precision-metric/v2` record required by that frozen
contract. Each invocation verifies the exact contract, executable, and canonical source-value hashes and writes its
CSV atomically. Twelve CPU-4 runs are bound to contract v1: four P4 well-separated shapes, six P4 rank-sensitive
shapes, and two KLIP well-separated shapes. All 138 applicable frozen threshold rows pass; all 48 applicable
report-only rows contain values and intentionally carry no pass/fail result. Every row remains
`promotion_ready=false` because required case families and owner-approved science/performance limits are incomplete.

The contract-bound reruns reproduce the main performance signal while quantifying its variability. At the bounded
`T=96, K=2,988, M=48` P4 shape, total arithmetic medians are 3.674 ms (DD), 2.161 ms (M32D64, 1.700x), and 2.101 ms
(F32, 1.748x). At `T=621, K=2,988, M=466`, they are 164.839, 104.446 (1.578x), and 103.601 ms (1.591x). For KLIP's
shared AF Lep shape, F32 halves eigensolve time but is 1.35% slower overall in this six-repetition sample because
subtraction dominates; in the selected pixel-Gram shape it is 1.632x faster overall. Rank-boundary behavior is
unchanged from the earlier diagnostic: differences are confined to the declared construction-specific ambiguity
bands, including raw ranks 61--63 for a structural-rank-60 zero-threshold case. Exact rows and hashes are in the run
manifest; neither these microbenchmarks nor the synthetic reduction tests replace bounded float application products.

A benchmark-only application adapter now reuses the exact `p4Reduce` configuration, loading, checking, reduction, and
finalization lifecycle while selecting one internal precision policy before the reduction worker loop. It is built only
when both the CPU-benchmark and experimental-P4-precision options are enabled, has no production configuration or
default, rejects unsupported optimizer and non-DD paths before mutation, and prints the effective candidate in the run
log. Its archived executable and adjacent experimental library are hash-bound in two predeclared case descriptors:
[`p4_fp32_aflep_bounded_e2e_case.yaml`](../../benchmarks/manifests/p4_fp32_aflep_bounded_e2e_case.yaml) and
[`p4_fp32_aflep_science_roi_e2e_case.yaml`](../../benchmarks/manifests/p4_fp32_aflep_science_roi_e2e_case.yaml).

The first fixed-one-worker end-to-end matrix used the original 96-frame `[10,12)` iteration annulus, one warmup and six
balanced fresh-process repetitions per candidate. Every process completed, every `P4-D64` product remained byte-identical
to the established baseline, and every policy was bitwise repeatable. Median whole-process times were 1.840 s (DD),
1.675 s (M32D64), and 1.660 s (F32); median paired improvements were 9.21% and 9.76%. The corresponding algorithm medians
were 0.9022, 0.7271, and 0.7171 s. That geometry is not a science-product case: a two-pixel support annulus cannot contain
the strict 4-by-4 derotation interpolation footprint, so all 196,608 final values are expected NaNs even though all 132
regressions are valid. The result remains useful for lifecycle timing, compatibility, metadata, and determinism, but
stored-value and residual-science metrics are explicitly unavailable. The balanced timing evidence is also preliminary:
the uncontrolled run series reached 99.05 C on the TCPU sensor.

The second predeclared case widened the calculation annulus to `[7,15)` while scoring the original 132-pixel `[10,12)`
science ROI. It added an ordinary default-off production control and admitted every fresh process only after TCPU was at
most 65 C and package temperature at most 70 C. The calculation performed 560 valid fits with `K=2,988`, rank 96, and
modes 4, 19, and 48. The production control and every runner-DD product are byte-identical; all candidates have identical
headers and finite masks, 336 finite values per plane globally, all 132 ROI values finite per plane, no invalid fit or
fallback, and one deterministic product hash per policy across all repetitions. The pre-run 176-finite-pixel support
estimate was conservative; the measured value is 336 per plane.

| Candidate | Whole-process median and paired improvement | P4-algorithm median and paired improvement | Median Gram / eigensolve | Median peak RSS |
|---|---:|---:|---:|---:|
| `P4-D64` | 4.995 s | 3.9433 s | 1.8303 / 0.3501 worker-s | 447,922 KiB |
| `P4-M32D64` | 4.040 s; 19.26% | 2.9947 s; 24.12% | 0.8538 / 0.3476 worker-s | 448,192 KiB |
| `P4-F32` | 4.015 s; 19.62% | 2.9676 s; 24.87% | 0.8567 / 0.2967 worker-s | 448,102 KiB |

Thus FP32 Gram construction is about 2.14x faster in the real reduction. Keeping only the eigensolve in FP64 costs
little here: full FP32 adds a roughly 1.18x eigensolve speedup but only about 0.36 percentage point to the paired
whole-process improvement over M32D64. Neither candidate reduces measured RSS at the current FP64 application/storage
boundary. The thermal admission gate equalized launches, but every measured process still ended at the TCPU sensor's
99.05 C reading; these timings remain engineering evidence and cannot satisfy the owner-TBD performance gate without a
non-throttling environment or controlled-frequency repeat.

Stored FP32 differences are small in RMS terms but visibly mode-dependent. Across all 1,008 corresponding finite values,
candidate-minus-DD RMS is 0.0616% of baseline RMS for M32D64 and 0.0621% for F32. Across the 396 ROI values it is 0.1013%
and 0.1087%; on the most sensitive 48-mode ROI plane it is 0.7787% and 0.8352%. ROI correlations remain at least
0.9999696 and 0.9999651, while maximum absolute ROI differences are 0.00397 and 0.00358 input units. Maximum ULP distances
are millions because the worst ULP locations are close to zero, confirming that the report-only ULP distribution is not
a science tolerance. The run manifest records every raw timing/RSS/thermal sample, mask/hash, per-plane and pooled
percentile, RMS, correlation, slope, ULP convention, and worst coordinate. No science or performance pass is assigned:
limits and aggregation remain owner-TBD, per-fit real-rank ambiguity is not exported, and injected-source, full-621-frame,
memory-policy, other-path, cross-platform, and KLIP end-to-end evidence remain open.

A third predeclared cross-worker case reran each policy once with four outer OpenMP workers pinned to four distinct
performance cores and one BLAS/LAPACK thread per worker. For DD, M32D64, and F32 alike, the complete stored science
array and finite mask are bitwise identical to that policy's fixed-one-worker reference: zero payload mismatches and
maximum ULP distance zero. All 134 FITS header cards also agree except the two predeclared maximum/effective-worker
cards changing from 1 to 4. This closes same-backend repeatability only for the detector-frame direct AF Lep science-ROI
path. The single samples are not a performance comparison, and automatic memory-limited concurrency, other paths,
other schedules, and other platforms remain open.

## Disposition (2026-09-03)

The owner selected `P4-M32D64` for the production direct, centered, explicit-held-out, shared-coefficient, local-trial,
and frozen-probe paths. P4 retains its public all-double numerical API as the regression oracle and uses its all-double
implementation for exact factor deletion and the associated whole-search-pixel fallback. `KLIP-M32D64` remains the
production direct policy; KLIP exact factor deletion remains an FP64 internal exception. Production FITS provenance
must distinguish calculation, direct-eigensolve, and factor-deletion precision so a configured exception is visible.

This is an owner engineering decision made with the limitations above in view, not a declaration that the frozen
acceptance manifest passed. The marginal additional whole-process gain from F32 eigensolve in the measured P4 case did
not justify its rank-policy and validation risk, and KLIP's F32 benefit was workload dependent. Further FP32
eigensolve, rank-policy, deletion, guarded-fallback, and cross-platform exploration is tabled until the GPU-acceleration
work resumes.

## Work sequence

The sequence below is retained as the evaluation history. Items left open after the disposition are deferred evidence,
not blockers on the accepted M32D64 implementation unless the production change itself requires them.

### 1. Freeze the baseline and acceptance record

- [x] Capture current focused Release results for `P4PCA`, `P4Reduction`, `KLIPreduction`, and both applications.
- [~] Archive current P4 and KLIP output products, timing, peak RSS, configuration, dependency/build manifest, and input
      checksums for the selected bounded and production cases.
- [x] Check in the numerical/science tolerance table before inspecting final candidate products. Distinguish exact
      structural contracts, roundoff comparisons, rank-ambiguity handling, and science limits.
- [ ] Agree and check in performance-materiality thresholds with their scientific/operational rationale. Initial
      values to review are a 20% median regression-stage or 10% median whole-application improvement, paired-run
      uncertainty excluding a slowdown, and no greater than 5% unexplained regression on another primary workload;
      these are proposals, not established policy.
- [~] Record current per-stage proportions so the maximum plausible end-to-end speedup is explicit before writing an
      optimization.

### 2. Close the FP32 capability gate

- [x] Compile and run minimal FP32 `syevrMem`/`eigenSYEVR` and `calcEigenVecs` cases through the linked production
      BLAS/LAPACK, confirming that `SSYEVR` is actually selected and returns the requested spectrum in both Gram
      orientations.
- [ ] Compile and run float `svdDeletionResult`, `svdDeletionWorkspace`, factor validation, row deletion, and column
      deletion for every backend in a future FP32-deletion evaluation. This is tabled with FP32 eigensolve work.
- [x] For every edited function containing one or more mxlib calls, verify in mxlib's current LCOV report that every
      called mxlib API has 100% executable-line coverage, and exercise the exact scalar specialization directly, as
      required by the repository coverage gate. Record every gap under `Known non-blocking ownership follow-ups` in
      `mxlib_cleanup.md`.
- [~] Confirm Eigen vectorization/alignment and mxlib ABI compatibility for both scalar types. Record an unsupported or
      scalar fallback rather than assuming the nominal float API is optimized.

### 3. Add internal precision variants and a CPU benchmark

- [~] Generalize the numerical kernels by calculation scalar and eigensolver scalar without duplicating algorithms.
      The normal-build DD seam and mixed/FP32 direct, centered, held-out, probe, and reduction-worker dispatch are
      complete. The accepted production adapter applies M32D64 to shared-PSF coefficients and local trials while
      retaining the public DD oracle. Typed FP32 factor deletion remains tabled.
- [~] Preserve the test eigensolver seam, finite checks, requested-spectrum behavior, operation ordering, and explicit
      checked P4 input/output boundaries for every variant.
- [~] Add a CPU-oriented benchmark that times input conversion, centering, Gram construction, eigensolve, rank
      selection, mode/factor construction, projection/residual formation, and optional coefficients/probes separately.
      It must build and run without enabling the hciReduce CUDA benchmark. Repeat the portability check with a
      non-CUDA mxlib build when one is available.
- [x] Make P4 worker-memory estimates scalar-aware and test the formula. The production estimate includes simultaneous
      FP32 ingress, Gram/eigenvector, projection, coefficient, centering, and probe scratch while leaving exact FP64
      factor deletion free of mixed-policy storage. Candidate measurements used an equal fixed worker count.
- [~] Emit machine-readable raw rows plus a human-readable summary. The covered P4 and KLIP cases now emit validated,
      contract-bound v2 rows and reject output/identity failures; the full required case sets are not yet present.

### 4. Characterize kernel accuracy and speed

- [ ] Run all synthetic spectral families and both Gram orientations through `P4-D64`, `P4-M32D64`, and `P4-F32` in
      a future FP32-eigensolve study.
- [~] Exercise P4 in-sample, centered, explicit held-out, coefficient, frozen-probe, and direct-fallback paths. The
      typed PCA core and internal reduction dispatch cover the first five, including detector, rotated, held-out, and
      target-held-out probe workers. Contract-bound real-product evidence and any typed factor-deletion fallback remain.
      If factor deletion is in the candidate's intended promotion scope, add its float variant only after the
      corresponding direct path passes.
- [~] Run KLIP shared-basis and per-target direct paths through `KLIP-M32D64` and `KLIP-F32` for both adaptive Gram
      orientations. Use an FP64 oracle rather than relying only on the existing FP32 test SVD.
- [ ] In a future FP32-deletion study, characterize KLIP and P4 factor-deletion precision separately, including factor
      validation, rank-boundary, clamp, deletion-failure, and matching-precision direct-fallback behavior. The accepted
      production scope deliberately retains the existing FP64 deletion paths.
- [~] Report stage speedups, conversion cost, workspace sizes, allocations, throughput, and numerical diagnostics. The
      bounded real-reduction cases now include stage timing, whole-process timing, RSS, and stored-product diagnostics;
      allocation/workspace accounting and the remaining path/platform cases are still open. Do not infer application
      speedup from eigensolve throughput alone.

### 5. Run bounded end-to-end comparisons

- [x] Use a quick deterministic synthetic/FITS case in every normal edit/test cycle.
- [~] Use the documented 96-frame AF Lep P4 profile to iterate on stage timing, followed by the 621-frame NACO/AF Lep
      data with `K=2,988` and the configured `K=1,776`--`4,108` annulus range and mode fractions. Fixed-one-worker
      `[10,12)` timing and finite-product `[7,15)`/scored-`[10,12)` matrices are complete; repeat performance without
      thermal saturation, then run the complete configuration.
- [~] Compare fixed and parallel P4 workers with controlled inner threading. The direct detector-frame science-ROI
      cross-check is bitwise exact for all three policies at one versus four workers; automatic memory-policy
      concurrency and the remaining P4 paths/schedules/platforms are still open.
- [ ] Use the documented bounded KLIP NACO case (128 of 621 frames, up to 64 modes), then the complete configuration. Cover
      both reference-space and pixel-space Gram regimes with a synthetic or additional real region.
- [ ] Include a larger-`T` P4 case (the documented greater-than-8,000-frame scale if accessible) because compute,
      memory, Gram orientation, and rank behavior differ from the 621-frame case. Mark it unavailable rather than
      silently omitting the regime.
- [ ] Run detector-frame P4 first, then centered rotated P4; no exclusion and explicit held-out P4 first, then exact
      factor deletion if deletion is in the candidate scope; ordinary full reduction first, then
      local/optimizer/PSF/filter products.
- [ ] Run KLIP ADI and RDI, all references and selected/excluded references, shared and per-target bases, then optional
      factor deletion. Preserve reference-selection flags exactly.
- [ ] Repeat on the principal production CPU and at least one materially different CPU/BLAS platform when available;
      a single host cannot establish a universal default.

### 6. Make separate P4 and KLIP decisions

- [x] Promote `P4-M32D64` for ordinary, centered, explicit-held-out, local-trial, shared-coefficient, and frozen-probe
      work. Retain exact factor deletion, its whole-search-pixel fallback, and the public numerical oracle in FP64.
- [x] Retain `KLIP-M32D64` for direct KLIP. Keep exact factor deletion as its documented FP64 internal exception.
- [x] Table the observed F32 rank/validity questions rather than broadening production eigensolve precision. Revisit
      them as part of future GPU acceleration.
- [x] Update HCI-015 with the benchmark record and disposition. The production implementation must preserve a fixed
      policy, explicit provenance, and the retained FP64 oracles/exceptions.

### 7. Production follow-on, only for an accepted candidate

This is the production follow-on checklist for the accepted M32D64 policy.

- [x] Promote the fixed M32D64 policy without a runtime precision selector. Record effective
      calculation/eigensolver/deletion precision in FITS provenance, including experimental benchmark products.
- [x] Retain the authoritative all-double P4 baseline for regression tests and factor-deletion operation while the new
      path gains production experience.
- [~] Update numerical documentation, memory accounting, configuration docs, application help, and every affected
      benchmark/test oracle.
- [ ] Run focused Release, Debug/coverage, ASan/UBSan, full CTest, documentation, downstream shared/static consumer,
      and real-data product comparisons. Apply the mxlib coverage gate to every edited caller.

The 2026-09-03 review follow-up closed the production memory-accounting underestimate, made experimental diagnostic
FITS cards follow the active D64/M32D64/F32 dispatch, added an always-on reduction-level Gram-routing discriminator,
and applied `clang-format` to the mixed-precision blocks. The default Release build passed all 26 CTest targets; its
P4Reduction executable passed 75,478 assertions. The experimental Release P4Reduction executable passed 294,180
assertions across 28 cases, and the experimental benchmark target built successfully. The edited FITS-header caller's
mxlib dependencies remained at 100% executable-line coverage. Debug/coverage, sanitizers, documentation, downstream
consumers, and additional real-data comparisons remain open under the aggregate validation item above.

## Promotion gates

The promotion gates below are the original frozen contract. They remain useful requirements for any future broader
precision promotion, but the 2026-09-03 owner decision accepted the narrower fixed M32D64 policy without declaring all
of these gates passed:

1. Every non-ambiguous rank decision, mode-status entry, selection flag, finite/validity bit, plane meaning, and
   fallback contract agrees with the authoritative baseline. Every ambiguity-band case is identified and resolved by
   an approved policy.
2. Synthetic kernel residuals and projectors pass scale-aware bounds against the FP64 normal-equation and independent
   direct-SVD/QR oracles; no candidate introduces unexplained nonfinite values, overflow, or solver failure.
3. Real and injected products pass the predeclared residual, throughput/photometry, astrometry, SNR/ranking, optimizer,
   and PSF/filter bounds at all tested mode counts, not only a visually selected plane.
4. Serial and parallel results meet the same scientific limits with controlled inner threading, and reduced worker
   storage does not invalidate the memory-budget model.
5. Performance is material and repeatable at application scale under the owner-approved thresholds frozen in the
   acceptance record before final candidate measurements.
6. A memory-only candidate either enables a documented previously infeasible workload or translates its smaller
   per-worker storage into a measured end-to-end improvement. Use P4's existing automatic worker policy; coordinate
   any future KLIP automatic-worker experiment with its separate memory-management plan.

The frozen manifests continue to report their gate results exactly. Future FP32 eigensolve or deletion work must not
treat the M32D64 decision as evidence that these broader gates passed.

## Deliverables

- A source-level scalar/operation inventory separating hot science arithmetic from non-candidate doubles.
- A reproducible environment and dataset manifest with checksums and controlled thread settings.
- A CPU benchmark and machine-readable raw results for every precision configuration and matrix/path regime.
- A frozen tolerance table plus synthetic numerical, real-product, and injected-source comparison reports.
- [x] A P4 decision and a separate KLIP decision, including FP64 factor-deletion exceptions and platform limitations.
- [x] The accepted production implementation, precision provenance, numerical documentation, and validation record.

## Risks and deferred decisions

- Forming normal equations squares the design matrix's condition number. FP32 may be adequate for leading modes yet
  fail near rank or retained-mode boundaries; this cannot be inferred from final float storage alone.
- P4's current `1e-12` relative rank threshold does not by itself express an FP32-resolvable cutoff. Altering it would
  be a scientific/numerical-policy change and needs explicit provenance and validation.
- Clustered eigenspaces can change basis while preserving the projector. Conversely, a cluster split by a requested
  mode boundary can change the subtraction even when both eigensolvers are individually correct.
- FP32 centering and long dot products can accumulate error before the eigensolve; an FP64 eigensolve cannot recover
  information already lost in an FP32 Gram matrix.
- Float factor deletion changes epsilon-scaled orthogonality, PSD/clamp, and rank-boundary tolerances and may change
  fallback frequency. It must not inherit hard-coded double-epsilon bounds. Mechanically substituting FP32 epsilon in
  the current dimension-scaled formulas would produce relative scales of roughly `0.0047` to `0.019` at dimension
  621, which is far too large to accept without re-deriving and testing the policies.
- Halved P4 worker storage can increase its automatic OpenMP concurrency, which may be an operational gain but also a
  confounder. Future broader performance characterization should compare fixed-worker and automatic-worker results.
  KLIP remains a fixed-worker/RSS comparison unless its separate memory-management work supplies an automatic policy.
- BLAS/LAPACK FP32-to-FP64 speed ratios vary by CPU, library, and matrix dimension. Results from one GPU, laptop CPU,
  or eigensolve-only kernel are not a general production conclusion.
- A conditional fallback, runtime precision option, platform-dependent default, or FP32 eigensolve is tabled until
  GPU acceleration resumes. CUDA and lower-than-FP32 arithmetic remain separate issues.
- The wide explicit-held-out FP32 coefficient path currently exports calculation coefficients through FP64 and recasts
  them to FP32 after recorded PCA timing. Account for that conversion or retain typed scratch before interpreting
  held-out stage timing; the direct detector-frame AF Lep timings in this evaluation are unaffected.
- At the evaluation checkpoint, production probe-output alias rejection was exercised indirectly through the
  experimental wrapper, while exact DD reduction-dispatch equivalence covered only detector-frame direct reduction.
  Direct production alias coverage and centered, held-out, and frozen-probe regression cases belong in the production
  validation record.
- Extremely large synthetic dimensions can still overflow the held-out/probe result-column product or narrow an
  Eigen dimension/rank into the LAPACK `int` interface. Current real workloads are far smaller; add explicit checked
  dimension arithmetic before treating the experimental entry points as a supported external API.
