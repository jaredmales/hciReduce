# Future Distributed Grid Runner Plan

## Status

Deferred design draft, preserved for the future cross-algorithm runner.

Decision recorded 2026-08-14:

- Grid search is removed completely from `klipReduce`.
- `klipReduce` and the forthcoming new algorithm remain independent single-reduction executables.
- A future external runner will construct and execute grids across either algorithm.
- The likely implementation is Python so trial scheduling, restart, and distributed execution are not coupled to either
  C++ science executable.
- The scientific, coordinate, provenance, output, and failure-policy questions below remain relevant and are retained
  for that future design. The C++-local lifecycle proposal is historical context, not the selected implementation.

This document is intended to be iterative. Add decisions below as `ANSWER:` entries before implementation where the
scientific contract is not already clear.

## Historical klipReduce prototype

Before its removal, `klipReduce` advertised `mode=grid`, but its implementation was an incomplete prototype:

- It validates seven grid settings, then reaches an unconditional `exit(0)` at the first grid point.
- Its position-angle calculation is commented out, so the PA array is never populated.
- A second unconditional `exit(0)` marks the missing reduction-state reset.
- It mixes a separation offset in pixels with a position-angle offset in degrees as though they were Cartesian
  coordinates.
- The nominal output block is unreachable and does not uniquely associate a reduction product with its separation,
  position angle, contrast, and grid indices.
- It mutates one `KLIPreduction` object repeatedly even though file loading, fake injection, preprocessing, KLIP, and
  final processing all mutate that object's state.
- It uses process termination for internal control flow rather than reporting failure through the application API.

The test suite covered grid dispatch and every incomplete-configuration rejection, but deliberately did not encode
either unconditional exit as expected behavior. Removing this prototype restored `klipReduce` to its supported
single-reduction responsibilities.

## Future runner direction

The runner should treat each science executable as an external algorithm adapter with a stable single-run contract:

1. Generate the ordered trial tuples and a run identifier from a validated grid specification.
2. Materialize an algorithm-specific configuration for one tuple.
3. Invoke the selected executable as an isolated process, locally or through a distributed backend.
4. Capture its command, configuration fingerprint, standard output/error, exit status, timing, and product paths.
5. Commit the trial result to an atomic manifest so interrupted runs can be inspected and safely resumed.

The common runner owns grid generation, scheduling, concurrency, retry/resume policy, and the cross-algorithm manifest.
Each science executable owns only one reduction and its scientifically meaningful products. Algorithm adapters should
remain thin translations from a common trial description to executable-specific configuration.

The implementation language and distribution backends remain decisions, though Python is the current preference.

## Proposed future use case

The future runner should implement an independent fake-companion parameter sweep:

1. Define a polar grid of sky locations around a requested `(separation, position angle)` center.
2. At each grid location, run every requested contrast as a separate trial.
3. For each trial, begin with the same pristine configured observation.
4. Inject exactly one trial companion into the target sequence using the existing ADI convention.
5. Run the normal loading, preprocessing, KLIP subtraction, derotation, combination, and configured output stages.
6. Write a manifest and self-identifying products so every result can be mapped back to its trial parameters without
   relying on directory order or an implicit sequential filename.

Independent trials are the recommended baseline. Injecting all grid points into one cube would allow companions to
affect one another's KLIP basis and would not produce independently attributable recovery measurements.

The initial implementation should produce reduction products, not claim to measure recovered flux or select a
best-fitting companion. Automated photometry, a merit function, interpolation, and optimization can be added later
once their statistical contract is defined.

## Decisions to settle

### 1. Scientific objective

Recommended: treat grid mode as an injection/recovery product generator. One trial is one `(separation, PA, contrast)`
tuple and produces the same configured KLIP outputs as basic mode plus grid provenance.

Alternatives:

- A negative-companion fitting search that calculates a scalar residual metric and selects an optimum.
- A throughput study that injects positive companions and measures recovered flux.
- A bulk visualization mode that injects multiple well-separated companions into one reduction.

`ANSWER:`

### 2. Coordinate system

Recommended: define the grid directly in polar coordinates:

```text
separation(i) = centerSep + i * deltaSep
PA(j)         = normalize(centerPA + j * deltaPA)
```

where integer offsets span all samples within their configured half-widths and always include the center. Separation
is in pixels, PA is in degrees east of north, and PA normalization uses the same convention as the rest of hciReduce.
Cartesian coordinates may be derived for diagnostics but should not be used to combine pixel and degree offsets.

Questions:

- Should a requested separation below zero be rejected or omitted when a centered grid crosses the origin?
- Should duplicate PA samples after normalization be rejected or deduplicated?
- Is the desired PA range `[0, 360)`, `(-180, 180]`, or simply the existing `angleMod` convention?

`ANSWER:`

### 3. Grid endpoints

Recommended: use integer offsets from `-floor(halfWidth / delta)` through `+floor(halfWidth / delta)`. This preserves
the requested step exactly, includes the center, and never exceeds the half-width. A half-width of zero should mean
one sample on that axis. A positive half-width requires a finite, positive delta.

This means a half-width that is not an integer multiple of its delta does not place a point exactly on the boundary.

`ANSWER:`

### 4. Contrast semantics

Recommended: require finite, nonzero contrasts and retain support for both signs. Positive contrasts support throughput
experiments; negative contrasts support negative-companion residual studies. Preserve the user's contrast ordering in
the trial order and manifest.

Questions:

- Is a zero-contrast control trial useful enough to allow explicitly?
- Are contrasts always linear planet-to-star ratios, or do we need a magnitude-difference input option?

`ANSWER:`

### 5. Pipeline recomputation

Recommended initial behavior: every trial reloads pristine inputs, injects before preprocessing, and reruns the full
configured pipeline. This matches the existing `ADIobservation::postReadFiles()` injection point and correctly allows
the injected source to participate in coaddition, nonlinear preprocessing, covariance construction, reference
selection, and KLIP basis generation.

A future fixed-basis or cached-input mode could be much faster, but it is a scientifically different operation unless
we prove which intermediate state may be reused. It should therefore be a separately named option with comparison
tests, not an invisible optimization.

`ANSWER:`

### 6. RDI behavior

The current target post-read hook injects configured fakes into target images. The RDI post-read hook only extracts
angles, despite existing RDI flux- and separation-scale configuration members.

Recommended initial behavior: inject into target data only. Reference-library injection should be a separate explicit
option because contaminating the RDI basis changes the experiment and is not normally desired for a companion present
only in the target system.

`ANSWER:`

### 7. Output layout and provenance

Recommended layout below `output.directory`:

```text
grid/
  manifest.csv
  sep_000_pa_000_contrast_000/
    <the configured final and/or PSF-subtracted products>
  sep_000_pa_000_contrast_001/
    ...
```

Each product should contain grid-specific FITS cards for at least:

- radial index, PA index, and contrast index;
- absolute separation, PA, and contrast;
- total trial count and stable linear trial index;
- whether the run recomputed the full pipeline or used a future reuse strategy.

The manifest should contain those values, the relative product path, status, and any failure message. Filenames should
use indices rather than formatted floating-point values; the manifest and FITS headers retain the exact values.

Questions:

- Should the configured `output.fileName` become a basename inside each trial directory, or should grid mode enforce a
  fixed basename?
- Should PSF-subtracted intermediate products be permitted for every trial, given their potentially large volume?
- Should grid mode write final combined images by default even when ordinary output is disabled?

`ANSWER:`

### 8. Failure and restart policy

Recommended initial behavior: fail fast on invalid configuration and pipeline errors. Write the manifest atomically
after each successful trial so completed work is discoverable. Add an explicit resume/skip-complete option only after
the product identity and configuration fingerprint are defined; never infer compatibility from file existence alone.

`ANSWER:`

### 9. Parallelism

Recommended initial behavior: run trials serially and retain the existing OpenMP parallelism inside each reduction.
Trial-level parallelism would multiply memory usage and oversubscribe OpenMP/BLAS threads. It can be added later with
an explicit worker count and resource policy.

`ANSWER:`

## Historical C++-local repair architecture

The following architecture was drafted before deciding to remove grid mode from `klipReduce`. It is retained because
its state-isolation and provenance analysis informs the future process runner. It should not be implemented inside
`klipReduce`.

### A. Pure grid construction

Introduce a small grid specification and generated-point representation. Keep generation independent of FITS I/O and
KLIP so it can be exhaustively unit tested.

Conceptually:

```text
GridSpec
  center separation / PA
  separation / PA half-widths
  separation / PA deltas
  ordered contrasts

GridPoint
  radial / PA / contrast indices
  stable linear index
  absolute separation / normalized PA / contrast
```

The generator will validate finiteness, ranges, steps, overflow, duplicate points, and the maximum trial count before
allocating the complete list. We should also add a configurable trial-count safety limit so a unit typo cannot launch
an unexpectedly huge run.

Whether these types remain private application implementation details or become a small reusable common component can
be decided during implementation. They should not become a public mxlib API.

### B. Pristine observation lifecycle

Use a configured but unread `KLIPreduction` as the prototype. For each trial:

1. Create a fresh trial observation from that pristine configuration.
2. Set its one-element fake separation, PA, and contrast vectors.
3. Set trial-specific output paths and provenance.
4. Load target and RDI file lists from their configured sources.
5. Call `regions()` once.
6. Destroy the trial state before starting the next tuple.

The first implementation may copy the pristine configured observation if focused tests prove that all owned image,
header, list, mask, and reduction state is independent. If that ownership proof is uncomfortable, retain/reapply the
parsed configuration to a newly constructed observation instead. Do not add a broad `m_filesRead=false` escape hatch:
that flag alone does not restore file lists, coaddition metadata, image cubes, headers, masks, derotator angles, KLIP
outputs, or timing state.

### C. Output isolation

Give each trial a dedicated output directory. This allows existing final-image and PSF-subtracted writers to retain
their normal naming contracts without collisions and avoids dependence on sequential-filename discovery.

Add grid provenance through the normal KLIP final header path rather than reopening and editing FITS files after the
run. A small optional header/provenance input to the reduction is preferable to app-specific post-hoc mutation, but
this is an API decision and should remain scoped to hciReduce.

### D. Error handling

The obsolete `exit(0)` calls were removed with the `klipReduce` prototype. The future runner should treat each child
process's nonzero status as a failed trial, capture its diagnostics, and apply the configured failure/retry policy.

Every FITS, manifest, and directory write must be checked. A trial must not be marked complete until its configured
required outputs and manifest record have been flushed successfully.

### E. Performance follow-up

After the correctness-first implementation is validated on real data, profile time spent in:

- FITS reads and preprocessing;
- fake-image reading and shifting;
- covariance/reference selection and eigendecomposition;
- derotation/combination and output.

Only then consider caching. Candidate modes must state their scientific semantics:

- cache raw FITS data but clone it before injection;
- cache resized fake PSFs and per-image scale factors;
- cache preprocessing products only when injection is intentionally post-preprocessing;
- reuse a fixed KLIP basis only as an explicit fixed-basis mode.

## Future implementation phases

### Phase 0: settle contracts

- [ ] Answer the nine decision sections above.
- [ ] Provide one representative real-data invocation and describe the desired products.
- [ ] Decide the maximum expected grid dimensions, contrast count, and acceptable runtime/storage.

### Phase 1: runner model and validation

- [ ] Add the pure grid specification/point generator in the runner.
- [ ] Define the common trial schema and algorithm-adapter interface.
- [ ] Allow scientifically valid zero values where appropriate, especially PA zero and zero half-width.
- [ ] Reject non-finite values, invalid steps, negative separations according to the chosen contract, duplicate points,
  multiplication overflow, and excessive trial counts.
- [ ] Unit-test exact and non-integral boundaries, one-point axes, PA wrapping, contrast ordering, and invalid inputs.

### Phase 2: one complete external trial

- [ ] Invoke one selected science executable as an isolated process with a generated configuration.
- [ ] Verify the injection occurs once, before preprocessing, at the expected detector positions in every frame.
- [ ] Verify normal KLIP, derotation, combination, and output configuration remains effective.
- [ ] Add runner provenance, algorithm identity, configuration fingerprint, and deterministic output isolation.
- [ ] Exercise this path under ASan/UBSan and the mxlib coverage gate.

### Phase 3: complete sweep and manifest

- [ ] Execute the Cartesian product of radial, PA, and contrast samples in a documented stable order.
- [ ] Write and flush a manifest record for each completed trial.
- [ ] Test multiple points/contrasts for state isolation: no accumulated fake flux, reused header state, filename
  collision, or ordering drift.
- [ ] Test a mid-run error and the chosen fail-fast/partial-manifest behavior.

### Phase 4: scientific integration

- [ ] Compare a one-point grid result byte-for-tolerance with the equivalent basic-mode single-fake reduction.
- [ ] On a synthetic rotating sequence, confirm the injected companion lands at the requested derotated separation/PA.
- [ ] Verify positive and, if accepted, negative contrast behavior.
- [ ] Run a small real dataset and inspect final images, headers, manifest, runtime, and storage.

### Phase 5: coverage, documentation, and optional optimization

- [ ] Bring runner code to 100% executable-line and function coverage.
- [ ] Document the mode, coordinate convention, output layout, resource scaling, and restart behavior in the User's
  Guide and example configuration.
- [ ] Keep algorithm-side changes subject to the mxlib coverage gate; the external runner should not call mxlib.
- [ ] Profile the real workflow and propose caching or trial-level parallelism only with explicit equivalence/resource
  contracts and regression tests.

## Acceptance criteria

- Neither science executable contains grid scheduling or runner control flow.
- A one-point KLIP grid is scientifically equivalent to a direct `klipReduce` run with the same single fake.
- Every trial begins from identical pristine input and contains exactly its configured fake tuple.
- Polar grid values, wrapping, ordering, and output mapping are deterministic and documented.
- Every result is self-identifying through both the manifest and FITS provenance.
- Errors return a nonzero CLI status and preserve enough completed-state information to diagnose the run.
- Runner unit/integration coverage reaches 100%; both algorithm executables retain their own normal, sanitizer,
  documentation, and coverage verification.
