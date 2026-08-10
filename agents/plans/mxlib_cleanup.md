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
      three CPU error-generator tools, and passed 176/176 CTests. The filtered trace retains 141 mxlib-owned files at
      38.0% line coverage (5,311/13,972); genhtml displays 40.0% function coverage after alias filtering (716/1,792).
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
      header is missing or has no emitted line/function records. The manifest enforces 44 executable headers with
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

- [ ] Bring the `mx::fits::fitsFile<float, mx::verbose::vv>::write` Eigen image/cube overloads used by the hciReduce
      test fixture to 100% executable-line coverage. The current mxlib report covers 240/470 executable lines in
      `fitsFile.hpp`; successful float writes are exercised, but the relevant overloads and their error paths are not
      complete.
- [ ] Bring the mxlib configuration and file-logistics APIs exercised by edited `HCIobservation` functions to 100%
      executable-line coverage: `appConfigurator::add` / typed `operator()`, the float/string `readColumns` overloads,
      `getFileNames`, and `pathFilename`. The current LCOV report shows incomplete executable lines in each owning
      header/source, so the hciReduce configuration, list, threshold, and weight regressions do not close the
      dependency-owned error paths.
- [ ] Bring the `eigenCube<float>` mean, weighted mean, masked mean, median, masked median, and sigma-mean overloads
      called by `HCIobservation::preProcess_meanSub`, `coaddImages`, and `combineFinim` to 100% executable-line
      coverage. The new hciReduce tests cover their integration contracts, but mxlib must own complete overload and
      validation-path coverage.
- [ ] Bring `fitsHeader`, `fitsHeaderCard`, `ISO8601DateTimeStrMJD`, `get_curr_time`, and `pathFilename` operations used
      by `HCIobservation::coaddImages` to 100% executable-line coverage, including copied-header provenance, repeated
      HISTORY cards, singleton start/end/delta cards, and date conversion errors.
- [ ] Bring the mxlib azimuthal-kernel/precalculation calls in `HCIobservation::preProcess`, FITS-header append calls in
      `stdFitsHeader`, and FITS read/header-card calls in `readPSFSub` to 100% executable-line coverage. Renaming the
      configured widths to explicit half-widths necessarily touched all three containing functions; their dependency
      overloads are not all at 100% in the current report.
- [ ] Define a real deep-copy constructor (and move operations) for `mx::improc::eigenCube`. Its implicit copy
      constructor shallow-copies the owned buffer while the custom copy assignment is deep; a Phase 1 test fixture
      initially exposed this as a double free during stack unwinding.
- [ ] Decide whether to delete or modernize the four explicit legacy headers before enabling their placeholder targets.
- [ ] Resolve the incompatible duplicate `mx::AO::sim::wfMeasurement<T>` definitions in `generalIntegrator.hpp` and
      `leakyIntegrator.hpp`. They compile in isolated test translation units but cannot coexist safely in one program
      or one production explicit-instantiation translation unit.
- [ ] Consolidate the incompatible duplicate `mx::improc::xcorrPeakMethod` enums in `imageXCorrDiscrete.hpp` and
      `imageXCorrFFT.hpp`; both headers currently compile only when included in separate translation units.
- [ ] Continue the non-mechanical template sweep: free/member-function templates need exact representative signatures
      or real calls, while policy-, feature-, and user-type-dependent templates need an explicit supported-type
      decision before either test-local or production-library instantiation.
- [ ] Fix `pywfsSlopeReconstructor::calcMask()`: a configured FITS mask is loaded and then unconditionally overwritten
      by the following circular-pupil block.
- [ ] Add a temporary FITS-cube fixture for behavioral coverage of `turbSequence::turbFnames()` success and clamping.
- [ ] Extract the four error-generator tools' logic from their `main()` functions so their ownership placeholders can
      become behavioral unit tests.

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
