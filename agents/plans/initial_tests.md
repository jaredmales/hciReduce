# Tests needed for hciReduce

We need to create a test system for hciReduce.  We'll follow the general flow of the mxlib test system, including use of catch2.  See AGENTS.md for some guidance.

In parallel let's start a doxygen system, also following mxlib.  Docs should look just like the mxlib docs in appearance. We should include the same provisions for including test coverage output.

After getting organized and setting up the system, start with HCIobservation.hpp/cpp.  Propose a suite of tests.  We want to cover both logistics (i.e. does the file load the way it's supposed to, etc) and physics (do we apply the filter we think we do).

Take this in two parts: first do the organization, then document the proposed test suite below under "#plan".  Do not alter this prompt.

Note: a code review is in-scope, and this includes any identified issues with mxlib.

# plan

## Organization

- Use the mxlib-style Catch2 2.13.9 single-header runner, with one production test source per executable so failures
  remain isolated and template-heavy translation units do not accumulate in one binary.
- Provide on-demand `tests` / `hcireduceTests` targets, CTest registration through `HCIREDUCE_BUILD_TESTS`, and the
  `HCIREDUCE_ONE_TEST` / `hcireduceTestOneRun` workflow for focused development.
- Generate coverage in an isolated instrumented build.  Compile hciReduce after mxlib's exported optimization flags
  with `--coverage -O0 -fno-fast-math`, filter test/system sources, and embed the HTML report in Doxygen.
- Generate standalone Doxygen documentation with the same header, layout, Doxygen Awesome 2.3.4 theme, dark-mode
  toggle, and coverage page as mxlib.  Include test sources so production APIs acquire real `Referenced by` links.
- Keep fixtures deterministic and self-contained.  Tests will create tiny FITS and ASCII inputs at runtime in unique
  temporary directories and remove them through fixture teardown; no external observing data or network access is
  required.
- Use a small derived test harness to expose protected state and hooks.  Hide harness-only helpers with `\cond` and
  preserve links to real production calls with `#ifdef __DOXY_ONLY__` blocks where wrappers obscure them.
- Use Catch2 `TEST_CASE` for every top-level test; do not use `SCENARIO`.
- Add a brief Doxygen block before every `TEST_CASE`, naming the behavior and real production API.

## Dependency gate

- Complete the work in [mxlib_cleanup.md](mxlib_cleanup.md) before implementing the `HCIobservation` behavior tests
  below.  The existing hciReduce test, documentation, and coverage scaffold can remain in place while mxlib is cleaned.
- Fix dependency-owned behavior and add its primary regression test in mxlib.  After installing the corrected mxlib,
  retain only the hciReduce integration tests needed to verify configuration, orchestration, and scientific output.
- Do not copy or locally patch mxlib algorithms in hciReduce; this avoids maintaining the same correction in two
  repositories.

## Proposed `HCIobservation` test suite

Implement the phases in order.  A defect exposed by a focused test should be fixed before adding a broader test that
depends on the same behavior.

### Phase 0: harness and contracts [complete]

- [x] Add a `HCIobservation<float, mx::verbose::vv>` fixture exposing only the protected configuration/data needed by a
  given test, plus helpers for tiny FITS cubes, FITS headers, file lists, quality tables, and weight tables.
- [x] Pin OpenMP to one thread for numerical unit tests; add explicit one-thread-versus-many determinism coverage later.
- [x] Decide the following contracts before encoding assertions:
  - pixel time-series `rms` normalization uses root-mean-square (denominator `N`) or sample standard deviation
    (denominator `N-1`);
    ANSWER: use rms, N.
  - `minGoodFract` accepts exactly the requested fraction (recommended: inclusive `>=`);
    ANSWER: be inclusive, `>=`
  - median final combination either honors mask/minimum-good-fraction semantics or documents that it does not;
    ANSWER: it should honor the mask and min-good-fraction
  - azimuthal-USM configured widths are full widths or half widths (mxlib currently treats them as half widths);
    ANSWER: they should stay half-widths and have their variable names and config keys indicated this clearly.
  - even median-USM widths are rejected or normalized to an odd kernel size (recommended: reject explicitly).
    ANSWER: we should allow even widths, and handle the median by averaging the two central values.  mxlib should be updated if needed.

### Phase 1: configuration, lists, and pure operations [complete]

- [x] `setupConfig` / `loadConfig`: verify defaults, all keys, valid and invalid enum strings, numeric date units, boolean
  option types, mask/skip flags, thresholds, and output settings.  Verify duplicate or incorrectly typed option
  declarations are rejected by the test expectations.
- [x] `load_fileList` / `load_RDIfileList`: verify sorted directory scans by prefix/extension; list-file order and clearing;
  target and RDI relative-path prefixing; missing inputs; repeated loads and state reset.  Include a regression proving
  RDI loading neither mutates the target list nor leaves RDI paths relative.
- [x] `threshold`: cover below/equal/above values, stable order, basename matching, missing files/entries, malformed and
  duplicate rows, and the documented `qualityThreshold <= 0` disabled behavior.
- [x] `readWeights`: cover filename reordering, normalization (`[1,3] -> [0.25,0.75]`), missing/duplicate entries,
  insufficient rows, and zero/non-finite sums.
- [x] `preProcess_meanSub`: cover `none`, known mean-image and median-image cubes, mask reapplication, and invalid methods.
- [x] `preProcess_pixelTSNorm`: cover a known one-pixel series, zero variance, masked pixels, one-plane finite behavior,
  and the not-implemented sigma-clipped path.
- [x] `coaddImages`: cover disabled/no-op cases; mean and median; exact count and time boundaries; both limits together;
  final remainder groups; two-frame boundary cases; averaged MJD/keywords; start/end/delta/history cards; unrelated
  header provenance; no-date count-only coadds; invalid methods and mismatched metadata vectors; and RDI provenance.
- [x] `combineFinim`: cover `none`, mean, median, sigma mean, documented non-positive-sigma fallback, normalized weights,
  masks and the exact `minGoodFract` boundary, all-masked pixels, multiple reductions, and empty/mismatched cubes.

### Phase 2: FITS ingestion and masks [complete]

- [x] `readFiles`: verify filename order, centered crops for rectangular odd/even inputs, requested-size clamping, cube and
  dimension state, ISO-8601 and numeric dates, selected header propagation, NaN-to-zero handling, front/back deletion
  including oversized requests, skip-preprocess behavior, repeat reads, and hook order.
- [x] `readRDIFiles`: verify target-size precondition, undersized-reference rejection, independent RDI list/date
  keyword/format/unit handling, independent/input/no-mask modes, state flags, repeat reads, and RDI hook order.
- [x] `readMask` / `makeMaskCube`: verify exact-size masks, centered rectangular larger masks, inputs too small in either
  dimension, output plane count/content, and rebuilding after coaddition.

### Phase 3: preprocessing physics

- [x] Masking: verify excluded pixels remain excluded after every enabled stage and do not contaminate fitted/filter
  statistics.
- [x] Radial profile: use a synthetic radial background whose annular residual should be approximately zero plus a
  localized source; repeat with a heavily masked annulus to prove excluded zeros are omitted from the estimator.
- [x] Median USM: test constant and impulse responses against a direct small-kernel oracle, define edge behavior, require
  finite output, and exercise width validation.
- [x] Gaussian USM: test constant-to-zero response, an impulse against `delta - normalized Gaussian`, finite edges, and
  width/image-size limits.
- [x] Azimuthal USM: use radial and azimuthal synthetic backgrounds plus a localized source; verify configured width and
  `maxAz` semantics, even/odd image sizes, and a finite central pixel.
- [x] Add one small composite golden case that locks the documented stage order:
  mask -> radial profile -> median USM -> Gaussian USM -> azimuthal USM -> mean subtraction -> time-series
  normalization.
- Verification:
  - The focused suite passes 8 cases in normal and ASan/UBSan builds, including exact independent oracles for every
    filter stage and the complete composite sequence.
  - The full hciReduce suite passes 8/8, and the coverage build passes 8/8 with `HCIobservation.hpp` at 1028/1265
    executable lines (81.1%).
  - Doxygen renders all eight new cases under `HCIobservation_unit_tests` without adding a test/group warning.

### Phase 4: output and end-to-end behavior

- [x] `stdFitsHeader`: verify required cards, values, and provenance metadata.
- [x] `outputPreProcessed` / `outputRDIPreProcessed`: verify filenames, image count/content, headers, and nested-directory
  creation.
- [x] `writeFinim`: verify exact versus sequential naming, overwrite semantics, additional headers, and FITS round-trip.
- [x] `outputPSFSub`: verify every reduction/image, headers, and weight sidecars.
- [x] Run one tiny golden workflow through read -> coadd -> preprocess -> combine -> write and compare the final FITS data
  and essential headers to independently calculated values.
- Output contracts and corrections:
  - Preprocessed target images preserve rectangular row/column order; reference products use
    `<outputPrefix>RDI_######.fits`; both retain their corresponding input headers and append the standard metadata.
  - Final products recursively create nested directories, sequential names advance without overwriting, exact names
    overwrite as documented, and the corrected card is `MIN GOOD FRACTION`.
  - PSF-subtracted products validate cube/header/weight dimensions, include reduction and image indices, and create
    nested prefix directories plus the optional weight sidecar.
- Verification:
  - The focused seven-case output suite passes in normal and ASan/UBSan builds; the full and coverage suites pass 9/9.
  - `HCIobservation.hpp` is 1199/1335 executable lines (89.8%) in the final Phase 4 trace; the remaining uncovered
    output branches are low-level directory/file/stream failure translations.
  - Doxygen remains at the 101-warning legacy baseline and renders all output cases in `HCIobservation_unit_tests`.

### Phase 5: robustness and dependency regressions

- [x] Run the focused suite under AddressSanitizer and UndefinedBehaviorSanitizer; add MemorySanitizer where the toolchain
  permits.  Prioritize deletion bounds, empty/mismatched vectors and cubes, mask crops, and median kernel widths.
  - The complete 10-test suite passes under ASan+UBSan with leak detection disabled for the ptrace environment.
  - MemorySanitizer is unavailable on this host because no Clang toolchain is installed.
- [x] Compare deterministic outputs with one and multiple OpenMP threads.
  - A four-plane composite regression produces matching preprocessed cubes and final combinations with one and four
    OpenMP threads.
- [x] Confirm the mxlib regressions tracked in [mxlib_cleanup.md](mxlib_cleanup.md) against the installed dependency, then
  retain a small hciReduce integration assertion only where dependency behavior changes scientific output.
  - `/usr/local/include/mx/improc/imageFilters.hpp` is byte-identical to the mxlib `dev` header containing `9715752`
    and `b53c425`; no temporary include overlay remains.
  - Fresh normal, full executable, ASan/UBSan, documentation, and coverage builds all resolve the installed mxlib.
    Both runtime suites pass 10/10, and the dependency-sensitive preprocessing assertions retain their expected output.
- [x] Make `KLIPreduction::meanSubtract` deterministic for every accepted enum value.
  - `none` subtracts zero while retaining norm calculation; `imageMode` reports `notimpl`; unknown values report
    `invalidconfig`; all validation occurs before either cube or the norm vector is mutated.
  - Focused normal and ASan/UBSan tests pass. The exact mxlib `eigenCube` mean/median and image-median call paths remain
    at 100% executable-line coverage; KLIP true-RMS normalization is now calculated locally with denominator `N`.
- Verification:
  - Fresh installed-dependency normal, ASan/UBSan, and coverage suites pass 10/10; Doxygen stays at the 101-warning
    legacy baseline and adds the `KLIPreduction_unit_tests` group without group/test warnings.
  - The final installed-dependency coverage trace is 55.7% overall (1296/2328); `HCIobservation.hpp` is 1199/1335
    (89.8%) and the newly started `KLIPreduction.hpp` suite is 34/535 (6.4%).

### Phase 6: KLIP configuration, centering, and covariance selection

- [x] Register and load the complete current KLIP configuration surface.
  - `pixelTSNormMethod` is correctly described as a string, and the existing `m_pixelTSSigma` state now has a
    `klip.pixelTSSigma` configuration target.
  - Intermediate covariance, right-reason mask/projection, and inclusion-matrix FITS products are disabled by default.
    `klip.writeDiagnostics=true` enables them, and `klip.diagnosticDirectory` selects an automatically created output
    directory.
- [x] Complete scalar/image centering and normalization regressions.
  - Masked per-image mean and median centers exclude masked pixels and reapply the mask to references and targets.
  - `meanImage` and `medianImage` apply the reference-basis center to distinct target data rather than recomputing a
    target center.
  - Correlation normalization receives Euclidean norms rather than sums of squares. Pixel time-series normalization
    uses true RMS, applies the reference scale to targets, and maps zero/non-finite scales to finite zeros.
  - Unsupported methods, empty or mismatched cubes, and invalid masks are rejected before mutation.
- [x] Cover matrix extraction and covariance/reference selection.
  - Ordered and empty row/column selections are deterministic; inconsistent dimensions are rejected.
  - Image-number and angular exclusion boundaries are covered. Correlation selection keeps exactly the configured
    maximum with deterministic index tie-breaking and remains safe when exclusions leave fewer references.
  - A one-mode synthetic `worker()` reduction verifies the KL result and the disabled/enabled covariance diagnostic.
- [x] Implement and test every advertised reference-inclusion method.
  - `klip.includeMethod` accepts `all`, `corr`, `time`, `angle`, and `imno`; invalid strings and negative
    `includeRefNum` values are rejected during configuration.
  - Exclusion is applied before ranking. `all` retains every surviving reference regardless of a positive cap, while
    a zero cap retains every survivor for every method. Ranked ties resolve deterministically to the lower reference
    index.
  - Correlation uses the direct normalized target-reference dot product. Time uses absolute MJD separation, angle uses
    wrapped derotation-angle separation, and image-number selection uses absolute index separation.
  - RDI selection uses the independent reference count, timestamps, and derotation angles; a worker regression covers
    an RDI library with more references than target images. Missing or non-finite selection metadata is rejected before
    entering the OpenMP region.
- Verification:
  - The fresh full normal suite passes 10/10, and the focused ASan/UBSan suite passes 164 assertions in 19 test cases.
    Doxygen succeeds with 91 legacy warnings and no new test/group/helper diagnostics.
  - Coverage passes 10/10 at 74.2% overall (1842/2481). `KLIPreduction.hpp` rises to 414/676 (61.2%), `HCI.hpp`
    reaches 95/136 (69.9%), and `HCIobservation.hpp` remains 1199/1335 (89.8%).

### Phase 7: KLIP region and final-product orchestration [complete]

- [x] Exercise `regions()` with in-memory ADI and RDI cubes.
  - Multi-wedge ADI coverage verifies region geometry persistence, mask-aware extraction, multiple requested mode counts,
    inclusion-matrix dimensions, finite full-frame reconstruction, and preservation of masked pixels.
  - RDI coverage verifies an independent reference count and proves that RDI's effective no-exclusion behavior no
    longer permanently overwrites the configured ADI exclusion methods.
- [x] Validate region inputs before allocation or OpenMP work.
  - Empty/non-positive mode lists, mismatched geometry vectors, non-finite or reversed radii, inconsistent loaded cube
    dimensions, pixel exclusion at zero inner radius, and regions emptied by a mask now fail deterministically.
- [x] Exercise `finalProcess()` stage ordering and output gates.
  - Post-median subtraction, derotation, and mean combination are verified in order on a small exact cube.
  - Exact nested FITS output is round-tripped with both mode planes and the complete KLIP geometry, centering,
    normalization, right-reason, exclusion, and inclusion metadata.
- [x] Define and implement the saved-reduction resume contract for `processPSFSub()` / `readPSFSub()`.
  - The current `outputPSFSub()` format is authoritative: `<prefix>_<reduction>_<image><extension>`, with matching
    zero-based `REDUCTION` and `IMAGE` FITS cards. Loading infers both dimensions, requires a complete rectangular index
    grid, validates every header and image shape, and retains the first reduction's per-image headers for downstream
    derotation/provenance hooks.
  - `processPSFSub()` requires a positive `NMODES` list matching the inferred reduction count, restores that provenance,
    and applies the caller's current post-median, derotation, combination, weight-file, and output configuration. The
    saved weight sidecar is round-tripped through filename-based matching. Missing files, malformed names, mismatched
    header indices, missing mode metadata, and mismatched mode counts are covered.
- [x] Wire saved-reduction processing into the `klipReduce` application.
  - `mode=postprocess` loads `postprocess.directory`, `postprocess.prefix`, and `postprocess.extension` (default
    `.fits`), validates all three before execution, and dispatches directly to the restored `processPSFSub()` path.
    Historical `basic` and `normal` spellings remain accepted, while unknown modes fail deterministically.
  - An application-level test loads a real configuration, resumes a two-reduction FITS collection through `execute()`,
    and round-trips the configured final output. The generated command-line help exposes all three new options.
- Verification:
  - The optimized full suite and the coverage suite pass 11/11. Focused ASan/UBSan runs pass 19 assertions in both
    application cases, 144 assertions in eight output cases, and 248 assertions in 25 KLIP cases.
  - Coverage is 80.5% overall (2306/2863), with `klipReduce.cpp` at 121/222 (54.5%), `HCIobservation.hpp` at
    1287/1441 (89.3%), and `KLIPreduction.hpp` at 634/730 executable lines (86.8%). Every restored KLIP function and
    every application function relevant to postprocess configuration/dispatch is reached; the unrelated unfinished
    grid implementation and true process entry point account for most uncovered application lines.
  - Doxygen succeeds with 72 legacy warnings and no test/group diagnostics.

### Phase 8: ADI observation integration [complete]

- [x] Exercise target and RDI FITS-ingestion wrappers with independent derotators, propagated angle cards, and both
  configured and unconfigured failure paths.
- [x] Exercise target/RDI post-read and post-coadd angle hooks with valid and missing metadata.
- [x] Cover single-file and per-image-list fake injection, per-image scale lookup, padding/cropping, RDI scale factors,
  invalid vector/list/method inputs, missing files, and nested error propagation.
- [x] Cover zero- and nonzero-angle mask-cube construction, auxiliary angle/FITS products, PSF-subtracted image
  derotation, size validation, and ADI FITS provenance cards.
- [x] Correct fake-scale indexing to use the image index rather than the planet index, validate all parallel fake
  metadata before indexing, and preserve the input mask exactly when its derotation angle is zero.
- Verification:
  - The optimized full suite passes 12/12; the focused ASan/UBSan suite passes 63 assertions in six test cases.
  - Coverage passes 12/12 at 86.8% overall (2494/2873). `ADIobservation.hpp` is 295/295 executable lines (100%) and
    all 14 functions are reached. A malformed one-column scale row explicitly covers the unequal-column guard left
    reachable by `mx::ioutils::readColumns()`.
  - Doxygen renders all six cases under `ADIobservation_unit_tests` without new test/group diagnostics.
  - The mxlib coverage gate confirms the exact `get_curr_time`, `readColumns`, `pathFilename`, `fitsFile` read/write,
    `createDirectories`, and exception APIs used by the edited functions are at 100%. `rotateMask()` has no record in
    the current mxlib trace, so its focused coverage work is recorded in [mxlib_cleanup.md](mxlib_cleanup.md).

### Phase 9: HCI enumeration conversions [complete]

- [x] Exhaustively exercise every string-to-enum and enum-to-string conversion in `HCI.hpp`, including every
  documented spelling and every invalid string/enum error branch.
- [x] Correct the invalid-`includeMethod` diagnostic so it identifies the inclusion method rather than the exclusion
  method.
- Verification:
  - The focused optimized and ASan/UBSan targets pass 68 assertions in seven test cases, and the full coverage suite
    passes 13/13.
  - `HCI.hpp` is 136/136 executable lines and 14/14 functions (100%). Overall coverage rises to 87.8% (2522/2873),
    while `ADIobservation.hpp` and `ADIDerotator.hpp` remain at 100%.
  - Every mxlib exception API reached by the edited conversion functions is at 100% in the current mxlib trace.

### Phase 10: KLIP branch closure [complete]

- [x] Cover invalid current enum state during configuration, diagnostic directory/write failures, masked and unmasked
  centering variants, masked true-RMS normalization, invalid and non-finite covariance candidates, every minimum and
  maximum exclusion unit, reference-ranking metadata validation, right-reason diagnostics, region loading/state
  errors, preprocess-only completion, and malformed saved mode metadata.
- [x] Remove three structurally unreachable paths exposed by the coverage audit:
  - Target-mask dimensions were rechecked after target/reference dimensions and reference-mask dimensions had already
    established the same invariant.
  - `regions()` tested a negative return from `finalProcess()`, whose contract returns zero or throws.
  - `processPSFSub()` wrapped an exception from mxlib's explicitly exception-free `parseStringVector()` API.
- Verification:
  - The focused Release and ASan/UBSan targets pass 314 assertions in 32 test cases; the coverage suite passes 13/13.
  - `KLIPreduction.hpp` is 724/724 executable lines and 17/17 functions (100%). Overall coverage rises to 91.1%
    (2612/2867), with all four common headers other than `HCIobservation.hpp` now at 100%.
  - The mxlib coverage gate remains satisfied: `eigenCube`, `imageMedian`, `eigenSYRK`, `calcKLModes`, region geometry,
    `parseStringVector`, `ompLoopWatcher`, exception construction, and the diagnostic FITS/directory APIs used by the
    edited functions are all recorded at 100% for their HCI-used ranges in [mxlib_cleanup.md](mxlib_cleanup.md).

### Phase 11: HCIobservation branch closure [complete except platform fault injection]

- [x] Cover unequal-column, missing-file, missing-basename, empty-cube, non-finite/decreasing-date, zero-plane,
  masked/weighted combination, hook-failure, threshold-only child-process, mask-rebuild, and every output/readback
  error branch that can be induced deterministically with local fixtures.
- [x] Exercise target and RDI preprocessing/coaddition wrappers, every masked weighted/unweighted mean and sigma-mean
  overload, target/RDI/preprocessed/final/PSF-subtracted output failures, and malformed/duplicate/inconsistent saved
  PSF-subtracted collections.
- [x] Remove four paths proved structurally unreachable by their preceding invariants or mxlib contracts:
  - `threshold()` caught an exception from the exception-free `pathFilename()` helper.
  - Median and Gaussian preprocessing accumulated impossible negative returns after their arguments and kernel were
    already validated.
  - `readPSFSub()` rechecked complete indexed-grid membership after unique cardinality had established it, and
    rechecked zero image dimensions after `fitsFile::read()` had already rejected such a FITS input.
- [x] Flush the PSF-subtraction weight stream before checking it so a real write failure is observable and reported.
- Verification:
  - The full coverage suite passes 13/13. `HCIobservation.hpp` is 1424/1426 executable lines (99.9%) and 32/32
    functions; the only misses are the two-line translation of an operating-system error returned while incrementing
    a `std::filesystem::directory_iterator` after successful construction and entry inspection. Deterministically
    forcing that race requires filesystem fault injection, so the defensive production branch remains intact.
  - At this milestone the overall trace is 96.4% (2749/2852), with every other common header at 100%.
  - The mxlib coverage gate remains satisfied for the exact `readColumns`, path, FITS, header/card, directory,
    smoothing/filtering, combination, exception, and Git-provenance APIs called by the edited functions; their current
    100% records are enumerated in [mxlib_cleanup.md](mxlib_cleanup.md).

### Phase 12: klipReduce application branch closure [complete]

- [x] Cover command-line configuration-path extraction; unused, unknown, source-attributed, and positional-option
  diagnostics; non-postprocess geometry validation; and rejection of unsupported modes.
- [x] Run a real basic-mode FITS input through file discovery, ADI metadata extraction, preprocessing, and the
  documented preprocess-only stopping point.
- [x] Preserve the existing saved-product postprocess configuration and end-to-end dispatch coverage.
- [x] Remove the incomplete grid mode from `klipReduce`: delete its configuration targets, state, accepted mode,
  dispatch, unconditional process exits, and unreachable implementation. Preserve the scientific design work in the
  [future distributed grid-runner plan](grid_mode.md), which will orchestrate both KLIP and the forthcoming algorithm
  as independent executables.
- Verification:
  - The focused application suite passes 29 assertions in six test cases. The coverage suite passes 13/13;
    `klipReduce.cpp` is 102/102 executable lines and 7/7 functions (100%).
  - Overall coverage is 99.9% (2730/2732) with all 91 alias-filtered production functions reached. The only two
    remaining misses are the retained platform-fault translation in `HCIobservation::readPSFSub()` documented above.
  - The mxlib coverage gate remains satisfied: the edited application configuration and dispatch functions call the
    manifest-enforced `appConfigurator`, application lifecycle, exception, and saved-product APIs already recorded at
    100% in [mxlib_cleanup.md](mxlib_cleanup.md).

### Phase 13: adaptive smaller-Gram KLIP basis [complete]

- [x] Select the smaller KLIP Gram matrix automatically for each realized reference library.
  - For a vectorized region with `P` pixels and `R` retained references, preserve the legacy `X^T X` route when
    `R <= P`, including the square tie, and solve the pixel-space `X X^T` route only when `P < R`.
  - Preserve the configured output-plane order, largest-mode accumulation, duplicate/unordered requests, effective
    clamping to `min(P,R)`, exact-positive eigenvalue policy, and right-reason projection behavior without adding a
    configuration or output-provenance option.
  - Keep `cv.fits` as the complete symmetric reference-space covariance even when the pixel-space solver is selected.
    Normal use-all pixel-space reductions skip construction of the unused `R x R` matrix when diagnostics are off.
- [x] Repair worker safety and timing at the changed numerical seam.
  - Check every eigensolver return, reject an empty selected library, capture worker exceptions across OpenMP before
    rethrowing, share immutable master modes/projectors, and keep selected scratch and solver workspaces thread-local.
  - Replace raced timing accumulation with OpenMP reductions, reset it for each `regions()` call, and bracket the full
    multi-region calculation. The component totals are sums of per-call solve/mode/projection times and deliberately
    exclude Gram construction and reference selection.
- [x] Verify numerical equivalence and the expected performance crossover.
  - Direct-SVD and frozen legacy projectors/residuals agree across `P<R`, `R<P`, square, rank-deficient, zero,
    near-null, oversized-mode, unordered/duplicate-mode, workspace-reuse, RDI-selection, right-reason, and one-versus-
    three-thread cases. Solver allocation failure and malformed worker geometry propagate deterministically.
  - Controlled one-thread Release benchmarks using KLIP's production float arrays and double eigensolver improve
    `96 x 621` basis construction from 16.124 ms to 0.889 ms (18.13 times) and `32 x 621` from 15.086 ms to 0.202 ms
    (74.79 times). The retained reference/tie routes are unchanged within measurement noise, while pixel-route
    projector and residual relative differences stay below `2.5e-6`.
  - The `96 x 621` double eigensolver Gram storage falls from 2.942 MiB to 0.0703 MiB (41.84 times smaller). Per-target
    selection still constructs/collapses the full reference Gram before choosing its solve, so avoiding that preparation
    is a separate optimization opportunity.
- Verification:
  - The focused optimized and ASan/UBSan suites pass 1,133 assertions in 40 `TEST_CASE`s.
  - `KLIPreduction.hpp` is at 100% executable-line coverage; the only raw zero-count function record is the compiler's
    deleting-destructor alias, while every logical function is reached.
  - GCC ThreadSanitizer does not model libgomp's `omp critical` synchronization and reports the internally protected
    `ompLoopWatcher` counter as a race; the OpenMP contract, serial/three-thread equivalence tests, and independent
    source review confirm the changed worker state is target-disjoint or synchronized.
  - The mxlib gate remains satisfied for `eigenSYRK`, `calcEigenVecs`, `calcKLModes`, `isFinite`, `get_curr_time`,
    `ompLoopWatcher`, region insertion/masks, FITS diagnostics, and exception construction.

## Code-review findings that set test priority

Resolve these under focused regression tests before relying on the affected integration paths:

1. [resolved] `load_fileList(fileList, ...)` prefixes `m_fileList` instead of its `fileList` parameter, so an RDI list
   can corrupt the target list while leaving RDI entries unprefixed.
2. [resolved] RDI date and mask configuration is currently ignored in favor of target settings.
3. [resolved] Front/back deletions can erase outside the file-list bounds.
4. [resolved] Coadd grouping can violate count/time limits on the final pair, assumes dates exist for count-only
   coadds, uses target filenames for RDI history, and preserves headers from the wrong group members.
5. [resolved] The radial-profile stage includes zeroed masked pixels in its estimator, biasing the physical subtraction.
6. [resolved] Even median widths overrun mxlib's current `medianSmooth` workspace, while hciReduce also leaves filtered edge
   storage uninitialized.
7. [resolved] `combineFinim` dereferences empty input, does not mask median combinations, and does not implement its
   documented non-positive-sigma fallback.
8. [resolved] `readMask` uses the row dimension for both crop origins and can request an invalid block when only one
   dimension is oversized.
9. [resolved] `outputRDIPreProcessed` now writes distinct `RDI_`-suffixed files with RDI data and headers.
10. [resolved] `KLIPreduction::meanSubtract` implements `none` as zero subtraction and rejects `imageMode` before
    mutation because mode estimation is not implemented.

All mxlib-specific infrastructure, documentation, packaging, and scientific issues found during this review are
tracked in [mxlib_cleanup.md](mxlib_cleanup.md), which is the prerequisite work plan for this suite.
