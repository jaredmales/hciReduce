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
- [ ] Confirm the mxlib regressions tracked in [mxlib_cleanup.md](mxlib_cleanup.md) against the installed dependency, then
  retain a small hciReduce integration assertion only where dependency behavior changes scientific output.
  - The source-level dependency regressions and hciReduce integration assertions pass against the exact header overlay.
  - `/usr/local/include/mx/improc/imageFilters.hpp` still differs from mxlib `dev`; install commits `9715752` and
    `b53c425`, remove the overlay, and repeat the normal/full/coverage checks to close this item.
- [x] Make `KLIPreduction::meanSubtract` deterministic for every accepted enum value.
  - `none` subtracts zero while retaining norm calculation; `imageMode` reports `notimpl`; unknown values report
    `invalidconfig`; all validation occurs before either cube or the norm vector is mutated.
  - Focused normal and ASan/UBSan tests pass. The exact mxlib `eigenCube` mean/median, image-median, and vector-variance
    call paths remain at 100% executable-line coverage.
- Verification:
  - The normal, ASan/UBSan, and coverage suites pass 10/10; Doxygen stays at the 101-warning legacy baseline and adds
    the `KLIPreduction_unit_tests` group without group/test warnings.
  - The final source-level coverage trace is 51.5% overall; `HCIobservation.hpp` is 1199/1335 (89.8%) and the newly
    started `KLIPreduction.hpp` suite is 34/535 (6.4%).

## Code-review findings that set test priority

Resolve these under focused regression tests before relying on the affected integration paths:

1. `load_fileList(fileList, ...)` prefixes `m_fileList` instead of its `fileList` parameter, so an RDI list can corrupt
   the target list while leaving RDI entries unprefixed.
2. RDI date and mask configuration is currently ignored in favor of target settings.
3. Front/back deletions can erase outside the file-list bounds.
4. Coadd grouping can violate count/time limits on the final pair, assumes dates exist for count-only coadds, uses
   target filenames for RDI history, and preserves headers from the wrong group members.
5. The radial-profile stage includes zeroed masked pixels in its estimator, biasing the physical subtraction.
6. Even median widths overrun mxlib's current `medianSmooth` workspace, while hciReduce also leaves filtered edge
   storage uninitialized.
7. `combineFinim` dereferences empty input, does not mask median combinations, and does not implement its documented
   non-positive-sigma fallback.
8. `readMask` uses the row dimension for both crop origins and can request an invalid block when only one dimension is
   oversized.
9. [resolved] `outputRDIPreProcessed` now writes distinct `RDI_`-suffixed files with RDI data and headers.
10. [resolved] `KLIPreduction::meanSubtract` implements `none` as zero subtraction and rejects `imageMode` before
    mutation because mode estimation is not implemented.

All mxlib-specific infrastructure, documentation, packaging, and scientific issues found during this review are
tracked in [mxlib_cleanup.md](mxlib_cleanup.md), which is the prerequisite work plan for this suite.
