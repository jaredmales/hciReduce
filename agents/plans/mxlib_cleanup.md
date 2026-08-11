# mxlib cleanup before hciReduce tests

## Objective

Stabilize mxlib's test, coverage, documentation, packaging, and image-processing behavior before implementing the
hciReduce behavior suite.  Dependency-owned fixes and their primary regression tests belong in mxlib; hciReduce should
exercise them only through focused integration tests.

The hciReduce test/documentation scaffold may remain in place, but behavior-test implementation is gated on the exit
criteria in this plan.

## Test convention

- Use Catch2 `TEST_CASE` for every top-level test.  Do not add or retain `SCENARIO` / `SCENARIO_METHOD` declarations.
- Existing nested `SECTION`, `GIVEN`, `WHEN`, and `THEN` blocks may remain when they improve readability; the
  standardization applies to the top-level test declaration.
- Put a brief Doxygen block immediately before every `TEST_CASE`, state the behavior and real production API under
  test, and preserve real `Referenced by` links as required by `AGENTS.md`.
- Keep one production test source per executable and use the shared `tests/testMain.cpp`; individual test sources must
  not define their own Catch main.

## Ownership boundary

- Fix algorithms, numerical boundary semantics, reusable image operations, and mxlib build/install behavior in mxlib.
- Fix hciReduce configuration, FITS orchestration, target/RDI state, processing order, and application-specific output
  in hciReduce.
- When an mxlib defect changes hciReduce science output, put the exhaustive unit test in mxlib and one small end-to-end
  assertion in hciReduce.
- Reinstall mxlib and configure hciReduce from a fresh build tree after each dependency tranche; do not work around an
  uninstalled mxlib change with hciReduce include paths or copied code.

## Resolved scientific contracts

The `ANSWER:` decisions recorded in `initial_tests.md` establish these cleanup requirements:

- Pixel time-series RMS normalization uses the population root-mean-square denominator (`N`).
- `minGoodFract` is inclusive: exactly the requested good fraction is accepted (`>=`).
- Median final combination honors both the mask and `minGoodFract`.
- Azimuthal-USM radial and azimuthal widths remain half-widths; variable names and configuration keys must say so.
- Even median widths are valid.  The median is the arithmetic mean of the two central samples, so mxlib must allocate
  and process the actual even-sized window without overrunning its workspace.

## Phase 0: preserve state and establish baselines

- [x] Reconcile or explicitly preserve the current dirty mxlib worktree before editing; do not mix unrelated existing
      changes into cleanup commits.
- [x] Record the compiler, build options, optional CUDA/ISIO/OpenMP dependencies, installed prefix, and mxlib version
      used by hciReduce.
- [x] Capture baseline results for the default build, all current CTest entries, the docs target, install, and a small
      downstream `pkg-config` consumer.
- [x] Record the current test inventory and Doxygen warning count so each tranche has a measurable before/after result.
- [x] Keep the worktree organized into reviewable patches separated by test normalization, registration, coverage,
      packaging, scientific fixes, and documentation. Commit creation remains owner-controlled.

Baseline recorded 2026-08-09: mxlib 0.2.0 at `cfb45e5`, CMake 4.2.3, GCC 15.2.0, Doxygen 1.15.0,
CUDA toolkit 13.3.73, OpenMP 4.5, and `pkg-config` 2.5.1. The default configuration enabled CUDA, ISIO, OpenMP,
FFTW, and OpenBLAS. It built and installed from a source-local build directory; 25 of 26 registered tests passed, with
the sole failure caused by registering a CUDA test on a host with no CUDA device. Documentation completed with 675
warning records. The installed `pkg-config` consumer built and ran, while a true external build directory failed because
generated headers were written into the source tree.

## Phase 1: normalize and register the complete test suite

Initial audit: mxlib had 42 Catch2 sources plus one legacy standalone developer program. CMake registered 25 Catch2
sources unconditionally plus one CUDA-conditional source, leaving 16 Catch2 sources and the standalone program outside
the suite. There were 44 `SCENARIO` declarations across 23 test files. The first cleanup tranche brought that inventory
to 43 Catch2 sources plus the explicitly inventoried standalone program.

- [x] Mechanically convert all top-level `SCENARIO` / `SCENARIO_METHOD` declarations under `tests/include/` to
      `TEST_CASE` / `TEST_CASE_METHOD`, preserving tags, sections, assertions, and behavior.
- [x] Update mxlib's `AGENTS.md` to specify `TEST_CASE` only, matching hciReduce.
- [x] Add or repair the required Doxygen block before every `TEST_CASE`; remove legacy `\test` or prose-only links that
      manufacture API references.
- [x] Register every Catch2 source, retaining only genuinely feature-gated tests behind matching build options, and
      explicitly inventory the data-dependent standalone program outside CTest.
- [x] Remove the private Catch main from `appConfigurator_test.cpp` and run it through `tests/testMain.cpp` like every
      other source.
- [x] Verify `mxlibTests`, `mxlibTestRun`, `MXLIB_ONE_TEST`, `mxlibTestOneRun`, and direct CTest execution. The external
      OpenMP-enabled lean build passes all 40 tests through both `mxlibTestRun` and direct CTest; the focused
      `imageXCorrDiscrete` target passes 86 assertions.
- [x] Confirm each registered test has a unique target/name and that no test file is silently omitted. The lean
      configuration has 40 unique CTest names matching its 40 eligible Catch2 sources exactly.

### Complete production-file ownership inventory

A subsequent file-level audit found 193 project-owned `.hpp` files, one project-owned `.inc` data initializer, and 51
project-authored `.cpp` files after excluding generated files and bundled inih, optionparser, SOFA, and vendor sources.
The `.cpp` count includes four main-bearing error-enum build tools, which are inventoried even though their logic must
be extracted before Catch2 can call it directly. Headers and same-stem implementations share one canonical test source;
the inventory therefore requires 198 Catch2 sources, including the separately generated `error_t` interface test.

- [x] Add 155 canonical placeholder `TEST_CASE` sources so every project-owned header/implementation maps to a test.
      Each placeholder has an explicit behavioral-test `\todo` rather than pretending that inclusion is a behavioral
      assertion.
- [x] Replace the hand-maintained registration list with a strict, bidirectional CMake audit. Configuration fails for
      a missing production-to-test mapping, an orphan test, a source without `TEST_CASE`, a stale/overlapping feature
      gate, or an uninventoried standalone developer program.
- [x] Preserve one production test source per executable while compiling the shared Catch2 main translation unit only
      once. This keeps the complete inventory practical without changing test isolation.
- [x] Move the hard-coded `directPhaseSensor` developer program to `tests/manual/` and replace its canonical path with
      a documented Catch2 placeholder.
- [x] Explicitly inventory three obsolete duplicate AO headers (`directPhaseReconstructorD`,
      `directPhaseReconstructorOrtho`, and `pyramidSensorSepQuad`) plus the obsolete `psdFilterCuda` specialization
      without compiling their placeholders. Their guards, class names, dependencies, or specialized APIs conflict with
      maintained implementations; deletion or modernization is a separate owner decision, not something to conceal
      behind a feature gate.
- [x] Build and run the lean CUDA-off/ISIO-off matrix: 173/173 eligible tests pass. The remaining 21 sources are
      dependency-gated (17 CUDA, three ISIO, and one CUDA+ISIO), and the four obsolete headers are the explicit legacy
      exclusions above.
- [x] Put all 275 lexical `TEST_CASE` declarations into 197 leaf groups below User's Guide → Testing → Unit
      Tests. Doxygen renders all 274 active cases (one FFT case remains deliberately disabled by `#if 0`) and reports
      no test/group/duplicate/reference warnings.

## Phase 2: make coverage trustworthy

- [x] Instrument production object and library targets, not only test translation units.  Apply coverage options before
      production targets are created or through a dedicated interface target linked after exported optimization flags.
- [x] Ensure effective GCC compile ordering ends in `--coverage -O0 -fno-fast-math`, and link final test/library targets
      with `--coverage`.
- [x] Add `VERBATIM` to the coverage custom target so lcov receives literal recursive filters. First extract only the
      mxlib source root, then remove `tests/*`; this also avoids lcov's unanchored `/sys/*` matching mxlib's own
      `include/sys` and `source/sys` directories.
- [x] Verify `.gcno` files exist for every compiled production `.cpp`, execute the complete CTest suite, and confirm the
      report contains compiled implementation files as well as header/template instantiations. The final ISIO-enabled
      coverage build produced 226 `.gcno` files and 223 runtime `.gcda` files, including all 46 objlib units and the
      three CPU error-generator tools, and passed 176/176 CTests. The current filtered trace retains 146 mxlib-owned
      files at 48.8% line coverage (7,125/14,591); genhtml displays 51.1% function coverage after alias filtering
      (953/1,865).
      The three generator
      tools are intentionally present at zero (305 executable lines total) because their build-time executions are
      cleared before the unit-test run. No test, system, vendor, generated-header, or outside-source path survives the
      filters, while mxlib's own `include/sys` / `source/sys` paths remain in the report.
- [x] Decide whether the turnkey coverage target is GCC-only or add an explicit `llvm-cov gcov` path for Clang. The
      turnkey target is GCC/gcov-only and rejects other compiler families explicitly.
- [x] Verify `coverage_clean` removes only generated nested coverage state while preserving the parent build and
      unrelated files; the instrumented production object build is reproducible after cleaning.
- [x] Capture and merge an `lcov --initial` zero-counter baseline before the runtime trace so every compiled,
      instrumented translation unit remains in the denominator even if it is never loaded by a test.

The zero baseline cannot synthesize counters for uninstantiated inline/template header code: GCC only emits those
records after real API use or explicit instantiation. Include-only placeholders establish inventory ownership and
parsing in that test's include context, while their behavioral `\todo` items identify the work needed to make the
numeric header denominator more complete. The coverage documentation states this limitation explicitly rather than
presenting the report as a physical line census.

### Template instantiation ownership

- [x] Validate an explicit-instantiation probe with `aperturePhotometer<float>`. Compared with the include-only
      placeholder, the probe adds an honest zero-count trace for 95 executable header lines and nine source functions.
- [x] Move the common `aperturePhotometer<float>` and `aperturePhotometer<double>` instantiations into mxlib itself,
      with matching `extern template` declarations, so compile checking and the coverage denominator do not depend on
      a test-only instantiation. Keep arbitrary other template arguments available from the header.
- [x] Add test-translation-unit explicit instantiations for header-only class templates whose representative type is
      unambiguous, while keeping policy-, feature-, and user-type-dependent templates explicitly classified instead of
      inventing unsupported ABI promises. Forty-five test translation units now contain 76 Doxygen-hidden explicit
      instantiation definitions; all build in the lean matrix, while the feature-sensitive subset also builds with
      CUDA, ISIO, and OpenMP enabled.
- [x] Add a manifest-backed final-trace check so the coverage target fails when an expected instantiated template
      header is missing or has no emitted line/function records. The manifest enforces 48 executable headers with
      positive `LF` and `FNF` records; `units.hpp` and `tagT.hpp` are intentionally compile-owned but omitted because
      their constant/type-only instantiations emit no executable counters.
- [x] Record and repair deferred template-body compilation defects exposed by the instantiation sweep before promoting
      the affected headers into the enforced manifest. The sweep fixed missing includes, stale FITS/error APIs,
      undefined detector setters, invalid getter returns, non-void fallthroughs, and the gamma-distribution parameter
      typo; focused behavioral assertions now cover the repaired detector, leak, and aperture-photometry paths.

Coverage intentionally measures mxlib-owned `objlib` implementation, authored error-generator tools, and emitted
header/template code. The bundled levmar and SOFA object libraries are third-party code and remain outside the
instrumentation/reporting scope.

## Phase 3: correct downstream build and install behavior

- [x] Stop exporting `-O3`, `--fast-math`, and `-march=native` as mandatory consumer flags in `mxlib.pc`.  Keep required
      language/features/defines separate from mxlib's own optimization policy.
- [x] Verify Debug consumers are not silently re-optimized. Focused sanitizer builds and the fresh hciReduce coverage
      build both compile without inherited release optimization or fast-math flags.
- [x] Fix the `FFTWL_LIBDIRS_FLAGS` / `FFTW3L_LIBDIRS_FLAGS` typo in the FFTW long-double link-directory setup and add
      a configuration check that exposes it without another FFTW library masking the path. A sentinel `fftw3l.pc`
      with a unique library directory is now carried through to the generated `mxlib.pc`.
- [x] Provide a usable installed CMake package (`mxlibConfig.cmake` plus dependency discovery and usage requirements),
      or explicitly choose pkg-config as the supported consumer interface and remove the incomplete export.
- [x] Test installation to a non-`/usr/local` prefix and a clean downstream consumer without source-tree paths.
- [x] Ensure generated version headers reflect the installed revision and do not emit a stale "repository modified"
      pragma after a clean build/install. A synthetic clean snapshot was built, followed by an empty Git commit with no
      source timestamp change; the next library build refreshed both generated SHA values to the new commit while
      keeping the modified and untracked flags at zero. The installed library reports that same SHA. Unix Makefiles and
      Ninja were both exercised: Ninja rebuilt all then-current dependent production objects and relinked in the same
      invocation as the revision-only change, then preserved the generated-header, object, and library timestamps on
      the next unchanged build. `ninja -t missingdeps` reported no missing generated-file dependencies.
- [x] Generate error and revision headers exclusively under the build tree so a true external build directory leaves the
      source tree untouched.

## Phase 4: fix confirmed scientific and memory-safety defects

### `medianSmooth` kernel validation

- [x] Add focused `TEST_CASE`s for valid odd widths, zero/negative widths, even widths, widths larger than either image
      dimension, output edges, and constant/impulse inputs.
- [x] Prevent the current even-width buffer overrun while preserving even widths: allocate and fill exactly
      `width * width` samples and average the two central samples for an even-sized window.
- [x] Preserve and document the existing edge contract (caller initializes untouched edges), or intentionally replace
      it with a safer output-initialization contract and update callers/tests together.
- [x] Run the regression under AddressSanitizer and UndefinedBehaviorSanitizer. The focused image-filter executable
      passes 721 assertions in 10 test cases with both sanitizers enabled.

### `eigenCube` masked combination thresholds

- [x] Document and implement the inclusive `minGoodFract` boundary: a pixel with exactly the requested good fraction
      is accepted (`>=`), including `minGoodFract == 1` when all planes are good.
- [x] Apply the same boundary and range validation to masked mean, weighted mean, sigma mean, and weighted sigma mean.
- [x] Test fractions below, exactly at, and above the boundary; zero and one; all-masked data; mismatched masks/weights;
      and non-finite or out-of-range thresholds.
- [x] Use the cube's actual scalar type for invalid output rather than hard-coding `invalidNumber<float>()`.

### `azBoxKernel` and `precalcKernel` center/error handling

- [x] Add finite-kernel tests at the exact center and neighboring pixels for odd/even and non-square images.
- [x] Define the center orientation/shape without dividing by `rad0 == 0`; also handle the kernel's own `rad == 0`
      sample without producing NaNs.
- [x] Propagate `setKernel` failures during `precalcKernel` construction instead of silently retaining invalid or empty
      kernels.
- [x] Test zero/invalid widths, `maxAz` limits, normalization to unit sum, coordinate bounds, and error propagation.
      The complete image-filter executable passes 728 assertions in 10 test cases under both `-O3 --fast-math` and
      ASan/UBSan.

### `imageXCorrDiscrete` bounds and peak localization

- [x] Replace the asymmetric legacy regression fixture with centered reference crops and known zero, positive, and
      signed displacements.
- [x] Limit automatic and requested lag ranges to blocks that are valid for the target/reference dimensions.
- [x] Make centroid localization operate on finite, non-negative support around the correlation maximum rather than
      feeding a signed correlation plane directly to center-of-light.
- [x] Scale interpolated peak coordinates from the actual clipped row/column lag ranges and reject invalid tolerances
      or magnified dimensions.
- [x] Reject a reference mask that differs in either dimension, including the previously accepted one-dimension-only
      mismatch. The focused executable passes 86 assertions in one `TEST_CASE` under the project optimization flags.
- [x] Fall back to the discrete correlation maximum when geometry allows only a 3-sample lag axis, where the four-point
      cubic interpolator has no valid support. Zero-shift and boundary-shift regressions require finite exact results
      and an explicitly empty magnified image.
- [x] Use the same fallback when a 5-sample-or-larger correlation maximum lies on the outer boundary, whose cubic
      support is also truncated. Exact positive and negative boundary-shift regressions pass under the normal,
      instrumented, and ASan/UBSan builds.
- [x] Fall back before cubic magnification when the discrete correlation peak is non-positive or non-finite, because
      zero-filled unsupported magnification borders would otherwise win the peak search. A deterministic
      anti-correlated regression exercises a unique negative interior maximum; all 86 focused assertions are
      sanitizer-clean.

## Phase 5: repair documentation infrastructure and debt

- [x] Fix the malformed closing comment in `doc/groupdefs.dox` and other parser-breaking nested/unbalanced comments.
- [x] Remove the unsupported `dirs` navigation entry for current Doxygen.
- [x] Remove user-specific absolute `TAGFILES` paths and generate any mxlib tag file under the build tree, not the source
      tree.
- [x] Update stale test filenames, anchors, groups, and legacy test-link conventions after the complete suite is
      registered.
- [x] Fix the coverage placeholder's light-mode foreground variables and make the embedded report follow the
      parent page's dark-mode selection.
- [x] Remove or reconnect unused `coverage_prolog.html` / `coverage_epilog.html` assets.
- [x] Restore the class-page layout under Doxygen 1.15: update Doxygen Awesome from 2.3.4 to 2.4.2, disable the new
      right-hand page-outline panel, and constrain long member return types without breaking the mobile layout.
- [x] Drive the current Doxygen warning baseline down by category.  Parser/configuration/link warnings should reach
      zero; any deferred API-documentation warnings must be counted and recorded explicitly.

The fresh docs target reduced the baseline from 675 warnings plus one hard tag-file error to 512 warnings and no
errors. Parser/markup, taxonomy, unresolved-reference/citation, and layout/configuration categories are all zero. The
512 deferred warnings are API documentation debt: 415 unmatched or ambiguous declaration comments, 77 missing member
or parameter descriptions, and 20 invalid or duplicate parameter descriptions. The rebuilt `eigenCube` class page was
also checked in light and dark modes at desktop and narrow widths; the page outline, cramped member declarations, and
white breadcrumb/heading artifacts are gone.

## Additional audit items

- [x] Review all compiler warnings under `-Wall -Wextra -Wpedantic` and split confirmed defects from intentional
      template/feature-gated warnings.
- [x] Check the complete test suite with OpenMP enabled and disabled, and with optional CUDA/ISIO features in each
      supported configuration. The final lean CUDA-off/ISIO-off build passed 173/173 tests. The strict
      CUDA+ISIO+OpenMP build passed 194/194, with GPU cases reporting an explicit no-device skip on this host; the
      `pyramidSensor` test is correctly gated on both dependencies because its production header uses CUDA APIs.
- [x] Add ASan/UBSan CMake presets or documented build commands that do not inherit release optimization/fast-math
      flags.
- [x] Confirm the documentation and coverage workflows never write generated files into the source tree.

Known non-blocking ownership follow-ups:

- [~] Continue closing mxlib coverage gaps exposed by hciReduce integration tests.
  - Completed mxlib `dev` commits, each verified under its focused normal or coverage-instrumented target:
    - `57488e0`, `a081b61`, `7f37f22`, `81124f0`, and `0906081`: application lifecycle, configuration reload, and
      command-line/help coverage.
    - `5be0c4c`: command-line configurator state.
    - `32b583f`: path helpers.
    - `f6d4459` and `12a34fa`: HCI filename/value-column parsing.
    - `37a7fe6`: short filename filters.
    - `b8cef98`: duplicate command-line-only targets.
    - `974ce45`: missing column-file input.
    - `4309b8b`: unmasked and weighted `eigenCube` combinations.
    - `6c636cf`: FITS header-bearing image/cube write failures.
    - `d30df2e`: HCI FITS-header provenance.
    - `20ae905`: MJD-to-ISO header timestamps.
    - `f3b0c9c`: coadd FITS-card updates.
    - `23eb25a`: batch FITS image/header reads.
    - `a32ebba`: bulk FITS header-vector validation and error paths.
    - `29da072`: coadd start/end/delta metadata updates and invalid ISO-date input.
    - `8986d5e`: FITS Git-status HISTORY provenance for HCI final outputs.
    - `423f71f`: HCI raw-buffer batch FITS reads with selected headers.
    - `7a9d3fa`: HCI reduction-parameter FITS headers, including a boolean write/read round trip.
    - `70849e9`: HCI preprocessing radial-profile cleanup and azimuthal median filtering with a precomputed kernel.
    - `e83e3ae`: direct raw-buffer FITS image output with HCI-style headers.
    - `e16db6f`: persisted HCI COMMENT/HISTORY and typed parameter cards through FITS output.
    - `75c8be7`: nested, idempotent output-directory creation before HCI preprocessed-image writes.
    - `8db6073`: collision-safe sequential names for HCI final FITS outputs.
    - `01c851e`: complete `eigenCube` lifecycle and storage-mode coverage.
    - `e114dd5`: FITS file construction, state, dimensions, read-window reset, and pixel-coordinate allocation.
    - `85e35de`: raw-pointer, one-dimensional Eigen, configured-subset, and headerless batch FITS reads.
    - `6350b63`: FITS cube windows, direct writes, lazy opens, missing inputs, batch failures, and comment-free headers.
    - `da3cf8e`: FITS Eigen destination resize exception propagation.
    - `e838f22`: deterministic CFITSIO/allocation failure coverage, resource cleanup, correct header-write error
      propagation, and memory-safe one-dimensional subset handling.
    - `a49e3f3`: complete `fitsHeader`/`fitsHeaderCard` lifecycle, mutation, scalar-conversion, CONTINUE, CFITSIO-write,
      and deterministic allocation-failure coverage; transactional rollback for failed list/map mutations; safe
      missing-key lookup, exact duplicate-card erasure, self-assignment/self-append, and integer no-comment dispatch.
    - `decd9ca`: complete `timeUtils` calendar parsing, ISO8601/current/compact formatting, UTC-to-TAI conversion,
      timespec arithmetic/comparison/mean, clock, sleep, and external-failure coverage; safe broken-down-time failure
      propagation and support for the documented null output pointer in `timespecUTC2TAIMJD`.
    - `c30caa4`: complete `getFileNames` and path-conversion exception propagation/translation coverage, align the
      catch-all standard-exception result with the documented `std_exception` code, and add `fileUtils.hpp` to the
      permanent coverage manifest. The HCI-used `pathFilename` and `getFileNames` paths are both at 100%.
    - `a74ecca`: complete HCI string/double `readColumns` field parsing, skip-column, conversion, allocation, and
      stream-error coverage; prevent an empty comma-delimited field from indexing before the string, preserve read
      errors without appending partial records, and add `readColumns.hpp` to the permanent coverage manifest.
    - `353645e`: complete HCI scalar/vector `appConfigurator` access, defaults, indexed lookup, source logging, and
      empty-value coverage; make vector whitespace parsing bounds-safe for trailing empty fields; and add
      `appConfigurator.hpp` to the permanent coverage manifest.
    - `6d4f997`: complete defensive `azBoxKernel` coverage and masked radial-profile subtraction used by the HCI
      preprocessing sequence, including finite dimension overflow, inconsistent cached bounds, empty kernels, and
      off-center maximum-azimuth evaluation.
    - `0e88106`: replace the `eigenLapack` placeholder with behavioral `eigenSYRK`/`calcKLModes` coverage; accept valid
      zero-valued eigenvector components, remove an uninitialized LAPACK workspace-cache comparison, and release
      locally allocated workspaces on both `calcKLModes` error exits.
    - `f4fad51`: complete HCI radius/angle-image and annular-region coverage for centered rectangular images, radial
      and angular exclusions, wrapped sectors, clipped image bounds, optional masks, and negative angle normalization;
      add `imageMasks.hpp` to the permanent coverage manifest.
    - `e3d2e4d`: complete integer and real-valued `parseStringVector` coverage for the delimiter forms used by HCI
      reduction FITS headers, replace its deprecated internal `convertFromString` calls with `stoT`, and add
      `stringUtils.hpp` to the permanent coverage manifest.
    - `2accce5`: replace the `vectorUtils` placeholder with behavioral coverage for both `std::vector`
      `vectorVariance` overloads reached by HCI preprocessing, distinguish sample variance from root-mean-square in
      the regression values, and add `vectorUtils.hpp` to the permanent coverage manifest.
    - `a1693cb`: replace the `exception` placeholder with behavioral coverage for every constructor and accessor used
      by HCI, recursive standard and non-standard nested-exception extraction, and ladder-formatted reporting.
    - `3ff1eec`: complete the masked and unmasked `imageMedian` paths used by HCI, including both caller-owned and
      internally allocated workspaces.
    - `c7c72cd`: add behavioral cubic-convolution `imageRotate` coverage for HCI derotation, including output
      allocation, counterclockwise direction, flux preservation, interpolation, and rejected-edge pixels.
    - `44d4190`: complete `createDirectories` coverage for HCI output paths, including deterministic translation of
      an otherwise platform-dependent filesystem error into `error_t::filesystem` through the existing hidden test
      seam pattern.
  - Baseline verification:
    - CUDA-off/ISIO-on LCOV passed all 176 mxlib tests and generated
      `/tmp/mxlib-coverage-filelogistics/doc/html/coverage/index.html`.
    - Final LCOV after `a49e3f3` passed all 176 tests and regenerated
      `/tmp/mxlib-coverage-filelogistics/doc/html/coverage/index.html`: 48.2% overall executable-line coverage
      (7,023/14,565); HCI-relevant `eigenCube.hpp` is 100% (185/185 lines and 37/37 functions), `fitsFile.hpp`
      is 100% (491/491 lines and 46/46 functions), `fitsHeader.hpp` is 100% (259/259 lines and 30/30 functions),
      and `fitsHeaderCard.hpp` is 100% (563/563 lines and 80/80 functions).
    - Final LCOV after `decd9ca` passed all 176 tests and verified all 48 required records: 48.8% overall
      executable-line coverage (7,125/14,591), with `timeUtils.hpp` at 100% (30/30 lines and 4/4 alias-filtered
      functions) and `timeUtils.cpp` at 100% (140/140 lines and 24/24 functions).
    - Final LCOV after `c30caa4` passed all 176 tests and verified all 49 required records: 49.1% overall
      executable-line coverage (7,171/14,608) and 51.2% alias-filtered function coverage (957/1,868).
    - Final LCOV after `a74ecca` passed all 176 tests and verified all 50 required records: 49.5% overall
      executable-line coverage (7,240/14,629) and 51.4% alias-filtered function coverage (960/1,869), with
      `readColumns.hpp` at 100% (129/129 lines and 6/6 alias-filtered functions).
    - Final LCOV after `353645e` passed all 176 tests and verified all 51 required records: 49.8% overall
      executable-line coverage (7,282/14,631) and 51.8% alias-filtered function coverage (968/1,869), with
      `appConfigurator.hpp` at 100% (86/86 lines and 11/11 alias-filtered functions). Both `add` overloads are also
      completely covered in `appConfigurator.cpp`.
    - Final LCOV after `6d4f997` passed all 176 tests and verified all 51 required records: 49.9% overall
      executable-line coverage (7,307/14,631) and 51.8% alias-filtered function coverage (968/1,869).
      `imageFilters.hpp` has 100% function coverage (22/22), and every executable line in the HCI-used Gaussian,
      azimuthal, precalculated-kernel, filter, median-smoothing, and radial-profile call ranges is covered.
    - Final LCOV after `0e88106` passed all 176 tests and verified all 51 required records: 50.4% overall
      executable-line coverage (7,458/14,798) and 52.2% alias-filtered function coverage (977/1,873). The exact
      HCI-used `eigenSYRK` and `calcKLModes` ranges are at 100% (9/9 and 37/37 executable lines, respectively).
    - Final LCOV after `f4fad51` passed all 176 tests and verified all 52 required records: 50.6% overall
      executable-line coverage (7,512/14,851) and 52.3% alias-filtered function coverage (981/1,877). The exact
      HCI-used `radAngImage`, `annulusIndices`, `angleMod`, `angleDiff`, and `dtor` ranges and the annulus worker are
      all at 100% executable-line coverage.
    - Final LCOV after `e3d2e4d` passed all 176 tests and verified all 53 required records after a clean nested rebuild:
      50.6% overall executable-line coverage (7,532/14,871) and 52.3% alias-filtered function coverage (983/1,879).
      Both `parseStringVector` overloads are at 100% (10/10 executable lines each).
    - Final LCOV after `2accce5` passed all 176 tests and verified all 54 required records: 50.7% overall
      executable-line coverage (7,535/14,871) and 52.4% alias-filtered function coverage (984/1,879). The exact
      HCI-used `vectorVariance(const std::vector&)` overload is at 100% (3/3 executable lines), as is its
      supplied-mean helper (6/6 executable lines).
    - Final LCOV after `a1693cb` passed all 176 tests and verified all 54 required records: 50.8% overall
      executable-line coverage (7,563/14,887) and 52.7% alias-filtered function coverage (991/1,881).
      `exception.hpp` is at 100% (42/42 executable lines and 11/11 alias-filtered functions).
    - Final LCOV after `3ff1eec` passed all 176 tests and verified all 54 required records: 50.9% overall
      executable-line coverage (7,574/14,887) and 52.7% alias-filtered function coverage (991/1,881). The exact
      masked/unmasked `imageMedian` implementation and its public wrapper are at 100% (28/28 executable lines).
    - Final LCOV after `c7c72cd` passed all 176 tests and verified all 54 required records: 51.0% overall
      executable-line coverage (7,611/14,924) and 52.7% alias-filtered function coverage (992/1,882). The HCI-used
      `imageRotate` implementation is at 100% (37/37 executable lines).
    - Final LCOV after `44d4190` passed all 176 tests and verified all 54 required records: 51.0% overall
      executable-line coverage (7,617/14,929) and 52.8% alias-filtered function coverage (994/1,884). The HCI-used
      `createDirectories` implementation is at 100% (10/10 executable lines).
    - A fully instrumented ASan/UBSan `fitsFile` run with leak detection passed 345 assertions. It exposed and drove
      the fix for a pre-existing one-dimensional subimage coordinate overflow.
    - Fully instrumented ASan/UBSan `fitsHeader` and `fitsHeaderCard` runs with leak detection passed 477 assertions.
    - The fully instrumented ASan/UBSan `timeUtils` target passed all 18 test cases and 142 assertions with leak
      detection enabled.
    - The fully instrumented ASan/UBSan `fileUtils` target passed all 10 test cases and 117 assertions with leak
      detection enabled.
    - The fully instrumented ASan/UBSan `readColumns` target passed all 8 test cases and 98 assertions with leak
      detection enabled.
    - The fully instrumented ASan/UBSan `appConfigurator` target passed all 6 test cases and 151 assertions with leak
      detection enabled.
    - The fully instrumented ASan/UBSan `imageFilters` target passed all 13 test cases and 752 assertions with leak
      detection enabled.
    - A correctly instrumented ASan/UBSan `eigenLapack` target passed both test cases and 22 assertions with leak
      detection enabled, including both locally allocated workspace error exits.
    - The fully instrumented ASan/UBSan `imageMasks` target passed all 3 test cases and 29 assertions with leak
      detection enabled.
    - The fully instrumented ASan/UBSan `stringUtils` target passed both test cases and 129 assertions with leak
      detection enabled.
    - The fully instrumented ASan/UBSan `vectorUtils` target passed its sample-variance test case and both assertions
      with leak detection enabled.
    - The fully instrumented ASan/UBSan `exception` target passed all 3 test cases and 20 assertions with leak
      detection enabled.
    - The fully instrumented ASan/UBSan `imageUtils` target passed all 3 test cases and 15 assertions with leak
      detection enabled, including all `imageMedian` workspace paths.
    - The fully instrumented ASan/UBSan `imageTransforms` target passed both test cases and 11 assertions with leak
      detection enabled, including HCI-style cubic rotation.
  - Next verification:
    - Continue rerunning aggregate LCOV after subsequent focused commits before recording new percentages or claiming
      additional APIs are at 100%.
- [x] Bring the HCI FITS file I/O calls to 100% executable-line coverage.
  - APIs: `mx::fits::fitsFile<float, mx::verbose::vv>::write` Eigen image/cube overloads and the HCI bulk
    image/header `read` overload.
  - Completed coverage: successful float writes, unavailable-output errors, HCI-style Eigen-cube and raw-buffer
    multi-frame image/header reads, direct 2-D raw-buffer writes with persisted headers, empty inputs, mismatched
    header vectors, missing later files, and COMMENT/HISTORY plus typed parameter-card persistence.
  - Completed error coverage: deterministic open, dimension, size, close, allocation, subset-read, header-read,
    create-image, pixel-write, header-card, and final-close failures, including first and later batch elements.
  - Verification: `fitsFile.hpp` is 100% in the final filtered trace (491/491 executable lines and 46/46 functions),
    and the permanent coverage manifest now requires its counter record.
- [x] Bring HCI configuration and file-logistics calls to 100% executable-line coverage.
  - APIs: `appConfigurator::add`, typed `operator()`, float/string `readColumns`, `getFileNames`, and `pathFilename`.
  - Completed coverage: normal path helpers and filename filtering, plus deterministic allocation, filesystem,
    standard, and unknown-exception propagation/translation for path conversion and `getFileNames`; complete HCI
    filename/double `readColumns` parsing, skip-column, conversion, allocation, and stream-error paths; application
    lifecycle, reload, CLI filtering/unknown options, duplicate command-line-only registrations, and the HCI scalar
    and vector types through defaults, indexed lookup, empty values, and source logging. The HCI-used
    `pathFilename`, full `getFileNames` implementation, `readColumns.hpp`, and `appConfigurator.hpp` are at 100%;
    both `appConfigurator::add` overloads are completely covered as well. `createDirectories` is at 100%, including
    nested and already-existing output paths and deterministic translation of an unmapped filesystem error;
    `getSequentialFilename` covers existing-output collisions.
  - API follow-up: `appConfigurator::get` documents `-1` on conversion errors, but currently calls the exception-free
    `stoT` overload without requesting its error code, so malformed numeric configuration can produce a fallback
    value while returning success. Resolving whether to preserve the destination or store the fallback requires an
    explicit public error-contract decision.
- [x] Bring HCI `eigenCube<float>` combinations to 100% executable-line coverage.
  - APIs: mean, weighted mean, masked mean, median, masked median, and sigma mean, called by
    `HCIobservation::preProcess_meanSub`, `coaddImages`, and `combineFinim`.
  - Completed coverage: masked/weighted-masked combinations, threshold/shape validation, and unmasked
    mean/median/sigma-mean plus weighted mean/sigma-mean.
  - Verification: the final aggregate report records 100% for `eigenCube.hpp` (185/185 executable lines and 37/37
    functions).
- [x] Bring the HCI KL numerical calls to 100% executable-line coverage.
  - APIs: `mx::math::eigenSYRK` and `mx::math::calcKLModes`, called by `KLIPreduction::regions()` and
    `KLIPreduction::regions_radii()`.
  - Completed coverage: lower-triangle covariance construction, all- and limited-mode solutions, normalized and
    mutually orthogonal modes containing valid zero components, elapsed-time outputs, reusable caller-owned LAPACK
    workspace, incompatible reference geometry, and non-square covariance errors.
  - Correctness fixes: zero eigenvector components are no longer rejected as non-normal; the workspace cache no
    longer reads an uninitialized LAPACK output count; and locally allocated workspaces are released on both
    `calcKLModes` error returns.
  - Verification: the final aggregate trace records 100% for the exact HCI-used `eigenSYRK` (9/9) and `calcKLModes`
    (37/37) executable-line ranges. A fully instrumented ASan/UBSan run with leak detection passed 22 assertions.
- [x] Bring the HCI KL region-geometry calls to 100% executable-line coverage.
  - APIs: `mx::improc::radAngImage`, `mx::improc::annulusIndices`, `mx::math::angleMod`, `mx::math::angleDiff`, and
    `mx::math::dtor`, used to build and select KLIP annular sectors and enforce frame-separation angles.
  - Completed coverage: rectangular centered grids, cardinal degree angles, radius scaling, full and wrapped sectors,
    inner/outer radial rejection, clipping all four search bounds to the image, optional mask exclusion, and negative
    angle normalization.
  - Verification: the final aggregate trace records 100% for `radAngImage` (9/9), `annulusIndices` (2/2), its worker
    (40/40), `angleMod` (5/5), `angleDiff` (7/7), and `dtor` (2/2) executable lines. The image-mask header is now
    manifest-enforced, and an ASan/UBSan run with leak detection passed 29 assertions.
- [x] Bring HCI reduction-header vector parsing to 100% executable-line coverage.
  - API: the string-delimiter-set `mx::ioutils::parseStringVector` overload used by `KLIPreduction::loadReduction`,
    plus its default character-delimiter overload.
  - Completed coverage: integer mode lists, signed/fractional real-valued region bounds, replacement of existing vector
    contents, the default comma delimiter, multiple delimiter characters, and singleton input.
  - Cleanup: both overloads now call the supported `stoT` conversion API instead of mxlib's deprecated
    `convertFromString` wrapper.
  - Verification: both overloads are at 100% (10/10 executable lines each), `stringUtils.hpp` is manifest-enforced,
    and a fully instrumented ASan/UBSan run with leak detection passed 129 assertions.
  - HCI follow-up: `parseStringVector` returns `void` and always appends a final converted token, so
    `KLIPreduction::loadReduction()` cannot detect empty or malformed header content with its current `size() == 0`
    checks. HCI must validate header text/tokens explicitly or mxlib needs a separately designed status-returning API.
- [x] Bring the current HCI pixel-series variance dependency to 100% executable-line coverage.
  - API: `mx::math::vectorVariance(const std::vector&)`, currently called by
    `HCIobservation::preProcess_pixelTSNorm()`.
  - Completed coverage: the default overload computes the vector mean and delegates to the supplied-mean overload;
    the regression values distinguish sample variance about the computed mean from sample mean-square about zero.
  - Verification: the exact default overload is at 100% (3/3 executable lines), its supplied-mean helper is at 100%
    (6/6 executable lines), `vectorUtils.hpp` is manifest-enforced, and the sanitizer run passed with leak detection.
  - HCI follow-up: mxlib's documented implementation divides by `N - 1`, while the Phase-0 HCI contract requires a
    true root-mean-square with denominator `N`. HCI must replace this call and explicitly handle empty, one-plane,
    zero-valued, and non-finite series; the mxlib sample-variance API should not be redefined as RMS.
- [x] Bring HCI exception construction and nested-error reporting to 100% executable-line coverage.
  - APIs: all four `mx::exception` constructor forms, `what`, `message`, `file_name`, `line`, `code`,
    `unwind_exceptions`, and `print_exceptions`; HCI directly uses the code-only and code-plus-message constructors
    throughout validation and nested error propagation.
  - Completed coverage: default code/location, message-only, code-plus-message, code-only, all accessors, recursive
    standard nested exceptions, a non-standard nested payload, and exact ladder formatting.
  - Verification: `exception.hpp` is at 100% (42/42 executable lines and 11/11 alias-filtered functions), remains
    manifest-enforced, and the fully instrumented sanitizer run passed all 20 assertions with leak detection.
- [x] Bring HCI image-median calls to 100% executable-line coverage.
  - APIs: masked and unmasked `mx::improc::imageMedian` calls in `KLIPreduction::meanSubtract()`.
  - Completed coverage: even-sized unmasked data, a sparse valid-pixel mask that exercises both selected and skipped
    pixels, caller-retained workspace reuse, internal workspace allocation, and internal workspace release.
  - Verification: the shared implementation and unmasked public wrapper are at 100% (28/28 executable lines),
    `imageUtils.hpp` remains manifest-enforced, and the fully instrumented sanitizer run passed all 15 assertions.
- [x] Bring HCI image-rotation calls to 100% executable-line coverage.
  - API: `mx::improc::imageRotate` with `cubicConvolTransform<realT>`, used by
    `ADIobservation::outputPSFSub()` for derotation.
  - Completed coverage: automatic output allocation, exact counterclockwise impulse placement, cubic interpolation,
    flux preservation for an interior source, and the zero-valued edge branch where a complete kernel is unavailable.
  - Verification: `imageRotate` is at 100% (37/37 executable lines), `imageTransforms.hpp` remains
    manifest-enforced, and the fully instrumented sanitizer run passed all 11 assertions.
- [x] Bring HCI coadd metadata calls to 100% executable-line coverage.
  - APIs: `fitsHeader`, `fitsHeaderCard`, `ISO8601DateTimeStrMJD`, `get_curr_time`, and `pathFilename` used by
    `HCIobservation::coaddImages`.
  - Completed coverage: copied-header provenance, repeated HISTORY cards, MJD-to-ISO strings, singleton
    start/end/delta cards, invalid ISO-date input, final-output Git-status HISTORY provenance, and the scalar/string
    reduction-parameter card types emitted by `stdFitsHeader`.
  - Completed verification: `fitsHeader.hpp` is 100% (259/259 executable lines and 30/30 functions),
    `fitsHeaderCard.hpp` is 100% (563/563 executable lines and 80/80 functions); both are now required by the
    permanent coverage manifest. `timeUtils.hpp` is 100% (30/30 executable lines and 4/4 alias-filtered functions),
    `timeUtils.cpp` is 100% (140/140 executable lines and 24/24 functions), and the time header is also manifest
    enforced. The `pathFilename` implementation is also at 100% in the final trace.
- [x] Bring the remaining calls in recently edited HCI functions to 100% executable-line coverage.
  - `HCIobservation::preProcess`: completed `radprofim` plus `zeroNaNs` cleanup, Gaussian and median smoothing,
    azimuthal-kernel generation, kernel precalculation, and mean/median-filter calls in `70849e9` and `6d4f997`.
    The radial-profile regressions model the production cleanup sequence because profile interpolation can leave
    non-finite edge pixels before `zeroNaNs`, and now also verify that masked pixels are excluded and preserved.
  - `HCIobservation::stdFitsHeader`: the FITS-header append calls are covered through the `fitsHeader` and
    `fitsHeaderCard` suites, both of which are at 100%.
  - `HCIobservation::readPSFSub`: the FITS read and header-card calls are covered through the 100% `fitsFile`,
    `fitsHeader`, and `fitsHeaderCard` suites.
  - `ADIobservation` injection/output: the exact default cubic-convolution constructor, kernel, and `imageShift`
    specialization have counters on all 63 executable lines selected by that call path. The separately used
    `imageRotate` implementation is 37/37.
  - Remaining direct helpers: `fitsHeaderGitStatus` is 9/9, the default COMMENT/HISTORY tag constructors are 6/6,
    the two `error_report` overload bodies are 10/10, the HCI-used degree/radian and angle-normalization helpers are
    16/16, and `parentPath` plus `getSequentialFilename` are 21/21 executable lines in the current trace.
  - Verification: every executable line in the exact `gaussKernel`, `azBoxKernel`, `precalcKernel`, `filterImage`,
    `medianFilterImage`, `medianSmooth`, `radprofim`, and `zeroNaNs` overload ranges called by these HCI functions is
    covered. The remaining uncovered lines in `imageFilters.hpp` belong to unrelated smoothing APIs.
- [x] Define a real deep-copy constructor and move operations for `mx::improc::eigenCube`. Copying now allocates
      independent storage; moving and ownership transfer leave the source empty; and `clear()` resets ownership and
      the data pointer. Focused normal and ASan/UBSan regressions cover copy, assignment, move, and shallow transfer.
- [x] Fix `pywfsSlopeReconstructor::calcMask()`: a configured FITS mask is no longer overwritten by the circular-pupil
      path. A focused test writes a non-circular temporary FITS mask and verifies its content and derived measurement
      size after loading.

## Not related to hciReduce

- [ ] Decide whether to delete or modernize the four explicit legacy headers before enabling their placeholder targets.
- [ ] Resolve the incompatible duplicate `mx::AO::sim::wfMeasurement<T>` definitions in `generalIntegrator.hpp` and
      `leakyIntegrator.hpp`. They compile in isolated test translation units but cannot coexist safely in one program
      or one production explicit-instantiation translation unit.
- [ ] Consolidate the incompatible duplicate `mx::improc::xcorrPeakMethod` enums in `imageXCorrDiscrete.hpp` and
      `imageXCorrFFT.hpp`; both headers currently compile only when included in separate translation units.
- [ ] Continue the non-mechanical template sweep: free/member-function templates need exact representative signatures
      or real calls, while policy-, feature-, and user-type-dependent templates need an explicit supported-type
      decision before either test-local or production-library instantiation.
- [ ] Add a temporary FITS-cube fixture for behavioral coverage of `turbSequence::turbFnames()` success and clamping.
- [ ] Extract the four error-generator tools' logic from their `main()` functions so their ownership placeholders can
      become behavioral unit tests.
- [ ] Forward the documented `mean` argument through both `radprofim` overloads. The current implementations omit it
      when delegating to `radprof` and to the radius-image overload, so requesting mean profiles still selects the
      default median behavior; hciReduce currently uses the default median path.
- [ ] Forward the documented `pixbuf` argument from `annulusCoords` and `annulusIndices` to `annulusCoordsWorker`.
      Both wrappers currently drop nonzero buffers; hciReduce uses the unaffected default value of zero.

The focused sanitizer configuration uses `-DMXLIB_OPTIMIZE=-O1`,
`-DMXLIB_CXXFLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer"`, and matching
`CMAKE_EXE_LINKER_FLAGS` / `CMAKE_SHARED_LINKER_FLAGS`. Run the focused executable with
`ASAN_OPTIONS=detect_leaks=0 UBSAN_OPTIONS=halt_on_error=1`; LeakSanitizer cannot operate under this environment's
ptrace wrapper. The FFT regression fixtures now use fixed random seeds and float-relative Parseval tolerances; this
removed a nondeterministic absolute-tolerance failure observed only in the instrumented suite.

## Exit criteria before resuming hciReduce behavior tests

- [x] No project test source uses `SCENARIO` or defines a private Catch main.
- [x] All 198 canonical mxlib Catch2 sources are strictly inventoried. Every supported source is automatically
      registered in its eligible feature matrix; the one standalone developer program and four obsolete
      headers are explicitly inventoried outside CTest with concrete reasons.
- [x] Default, focused, complete CTest, install, and downstream-consumer builds pass from clean trees. The final clean
      OpenMP/Debug snapshot passed 40/40 CTests after the final xcorr and FFT changes, installed both library variants
      to a temporary prefix, and compiled and ran an umbrella-header consumer using only the installed `mxlib.pc`
      metadata.
- [x] Coverage includes compiled production `.cpp` files, uses literal filters, and is reproducible after cleaning.
- [x] The confirmed median, masked-combination, and azimuthal-kernel defects have regression tests and sanitizer-clean
      fixes.
- [x] Doxygen has no parser/configuration/broken-link warnings; the remaining 512 API-documentation warnings are
      enumerated above.
- [x] A freshly installed mxlib builds hciReduce without source-tree overrides or stale repository-status warnings.
      The Debug compile line contains no exported optimization flags, and `klipReduce` now has real dynamic dependencies
      on both `libhcireduce.so` and the temporary-prefix `libmxlib.so`.
- [x] Re-run the hciReduce infrastructure probe. A fresh coverage configuration places
      `--coverage -O0 -fno-fast-math` last, links the shared library and executable with coverage support, and produces
      `.gcno` files for all four hciReduce implementation units plus the CLI translation unit.
- [x] Complete Phase 0/1 of [initial_tests.md](initial_tests.md).
