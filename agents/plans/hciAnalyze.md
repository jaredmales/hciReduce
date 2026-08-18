# hciAnalyze

Please look at /home/jrmales/Source/mxApps/klipReduce/klipAnalyze.cpp

Here we want to adapt that to, at first, be only a basic SNR calculator.  We need to make it a full-fledged mx::application that accepts config files etc.  It needs to support the high and low pass filters.  Also want to add a planet pa (East of North, where north is up) and separation instead of x, y.  Also needs to support the fake planet keywords.  

This new application will be called `hciAnalyze`.

Please document a plan below, without altering this prompt.  We can being implementation after I have reviewed the plan.

# Plan

## Scope and first-release contract

Implement `hciAnalyze` as a standalone `mx::app::application` executable for measuring the SNR of one or more signals
in a final PSF-subtracted FITS image cube. The initial release excludes the legacy `grid.directory` aggregation mode,
centroid fitting, azimuthal-box filters, and diagnostic FITS products. Those features can be added later without
changing the basic measurement contract.

Adapt only the necessary parts of `/home/jrmales/Source/mxApps/klipReduce/klipAnalyze.cpp`: FITS-cube input,
invalid-pixel masking, optional Gaussian high/low-pass filtering, masked-annulus noise estimation, and aperture SNR
extraction. Implement this in a new source file; do not modify the legacy utility.

## Application and configuration

Create `src/apps/hciAnalyze.cpp` with the standard application lifecycle (`setupConfig`, `setConfigPathCL`,
`loadConfig`, `checkConfig`, and `execute`) and the exception-unwinding entry point used by the existing applications.
Set a local default config name of `hciAnalyze.conf`, support `--config`, and enable configuration-source reporting.

Register this proposed public interface:

| Group | Keys | Meaning |
| --- | --- | --- |
| root | `file` (`-f`) | Required input FITS image or cube. |
| `planet` | `sep`, `PA`, `R` | Explicit signal separations (pixels), PAs (degrees east of north), and exclusion-mask radii (pixels). |
| `snr` | `minRad`, `maxRad`, `apertureR` | Annulus overrides and aperture radius; non-positive annulus values derive from KLIP header metadata. |
| `filter` | `hpfGaussFW`, `lpfGaussFW` | Gaussian high-pass unsharp-mask and low-pass smoothing FWHMs in pixels; non-positive disables a filter. |
| `fake` | `sep`, `PA`, `contrast`, `R` | Explicit fake metadata, which supersedes `FAKESEP`, `FAKEPA`, and `FAKECONT` FITS keywords. |

The `planet`/`fake` coordinate convention will be canonical: the first image index is displayed FITS x and the second is
displayed FITS y. East is left on the displayed sky image. A PA of 0 degrees is north/up and a PA of +90 degrees is
east/left. Thus, matching `klipAnalyze`, `x = xCenter - sep * sin(PA)` and `y = yCenter + sep * cos(PA)`, with PA
converted from degrees. Use
the loaded image dimensions and centre `((xPixels - 1) / 2, (yPixels - 1) / 2)`. Do not carry forward
`klipAnalyze`'s ambiguous `planet.X`/`planet.Y` interface.

## Header metadata and validation

Resolve positions after reading the FITS header in this precedence order:

1. Complete explicit `planet.*` configuration.
2. Complete explicit `fake.*` configuration.
3. Complete `FAKESEP`, `FAKEPA`, and `FAKECONT` header metadata.

Parse the existing comma-separated fake-keyword encoding without mutating config state. Validate nonempty equal-length
`sep`/`PA` lists; require a present `contrast` list to match; and allow each `R` list to contain either one shared radius
or a radius per signal. Reject partial/mixed definitions, non-finite or negative separations/radii, a non-positive
aperture radius, and apertures/masks that cannot be placed meaningfully in the image. Explicit planets take precedence
over fake metadata; contrast is optional for a pure planet measurement.

Without an SNR-annulus override, require usable `REGMINR` and `REGMAXR` lists and derive bounds from their extrema.
Validate finite, nonnegative bounds with `minRad < maxRad`; report an actionable error rather than dereferencing an
empty legacy metadata vector. Parse `NMODES` or `FRACT NMODES` only to label per-plane results, never as a prerequisite
to compute SNR.

## Measurement implementation

1. Read the FITS image/cube with mxlib types and make a validity mask before replacing NaNs; invalid pixels remain
   excluded from noise and aperture measurement.
2. Apply high-pass Gaussian unsharp masking and then low-pass Gaussian smoothing per plane when configured, matching
   the legacy order and using mxlib filter helpers.
3. Build exclusion masks for every resolved signal and aperture masks for measurement, combine them with the validity
   mask, and call `stddevImageCube` in the validated radial annulus.
4. Form the legacy-style SNR cube while retaining validity masking, then report maximum aperture SNR for every signal
   and input plane. Emit stable, machine-readable rows: plane index, optional mode label, signal index, separation,
   PA, optional contrast, and SNR.

Do not write `planetMask.fits`, `apertureMask.fits`, or `snrc.fits` by default; an explicit diagnostic-output feature
can be designed later.

## Build, documentation, and tests

1. Add and install the `hciAnalyze` target in `src/CMakeLists.txt`, wire coverage options like the other apps, and add
   it to the aggregate test target and help test.
2. Add `tests/apps/hciAnalyze_test.cpp` using existing fixture helpers. Cover config registration/loading, north/east-left
   PA conversion, explicit-over-header precedence, `FAKESEP`/`FAKEPA`/`FAKECONT` parsing, vector/annulus failures,
   filter enable/disable behavior, NaN exclusion, and a deterministic cube with known per-plane SNRs. Follow the
   required Catch2 and Doxygen test-documentation conventions.
3. Add an `hciAnalyze` guide, configuration dictionary, group definitions, docs dependencies, and commented example
   config showing direct and header-driven fake-planet measurement.
4. Format only the changed source; build the executable and focused test, run focused tests and `hciAnalyze --help`,
   then broader tests if practical. Since the measurement calls mxlib APIs, inspect current mxlib LCOV coverage for
   each called API and record any shortfall in `agents/plans/mxlib_cleanup.md` under `Known non-blocking ownership
   follow-ups`.

## Completion criteria

`hciAnalyze --config example.conf` accepts normal mx application config files, locates a signal from the stated
east-of-north separation/PA convention or supported fake FITS keywords, applies requested filters, and prints validated
per-plane SNR measurements and always replaces the corresponding `_snr` FITS cube. The target, documentation, example,
and focused tests
build and pass.
