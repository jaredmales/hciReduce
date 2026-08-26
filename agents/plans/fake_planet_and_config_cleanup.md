# Fake Planet and Configuration Cleanup Plan

## Status

Implemented by `3b26865` and `c4dd646`. The corresponding completed tasks are tracked in
`agents/plans/developer_todo.md` under **Fake Planet Cleanup** and **Config Cleanup**.

## Goal

Make fake-source experiments usable directly on already-preprocessed input, distinguish known real companions from
injected fake companions throughout configuration and FITS provenance, and make every boolean command-line option
accept both an implicit true value and an explicit true/false value.

The cleanup must preserve existing INI files and bare boolean flags, keep P4 pixel-local trials non-accumulating, and
make the source represented by every output unambiguous to `hciAnalyze` and downstream tools.

## Current findings

### Fake injection is coupled to preprocessing

`ADIobservation::postReadFiles()` currently injects configured fakes only when `preProcess.skip=false`. This prevents
ordinary KLIP or full-field P4 reductions from injecting a companion into an input cube that is already at the
post-preprocessing algorithm stage. The workaround used for the positive-recovery experiment copies the configuration,
sets `preProcess.skip=false`, disables every actual preprocessing operation, enables `preProcess.only`, and writes a
second cube. That workaround is unnecessary I/O and makes the injection stage harder to reason about.

Removing the `skipPreProcess` condition without another guard would be incorrect. P4 pixel-local evaluation and the P4
negative-companion optimizer intentionally load one pristine post-preprocessing cube and sample each configured trial
on demand. Pre-injecting the configured trial during `readFiles()` would then inject it twice and would accumulate an
unwanted fixed source beneath every optimizer evaluation.

The full-frame and pixel-local contracts therefore need to be separated explicitly:

- ordinary full-field reductions inject a configured fake after input loading regardless of `preProcess.skip`;
- with `preProcess.skip=false`, the existing ordering remains injection followed by configured preprocessing;
- with `preProcess.skip=true`, injection occurs directly at the algorithm-input stage; and
- P4 pixel-local evaluation, including optimization and jackknife refits, suppresses inherited full-cube injection and
  applies its finite-amplitude trial and any requested fixed planet subtraction on demand.

### Preprocessing-only outputs lose fake provenance

`HCIobservation::outputPreProcessed()` writes the common HCI header but does not call the ADI fake-metadata writer.
Consequently, the positive-recovery `frame_*.fits` files contain the injected signal but not `FAKEFILE`, `FAKESEP`,
`FAKEPA`, or `FAKECONT`. Preprocessed products must carry the same source provenance as final products so a subsequent
reduction and `hciAnalyze` do not have to reconstruct it from an external script.

### Planet and fake sources are currently conflated

Reduction configuration and final headers have only the `fake` namespace and `FAKE*` keywords. `hciAnalyze` separately
has a `planet` namespace, but reduction applications cannot persist the known companion in an equivalent,
application-independent form. A final image containing both a known companion and injected fakes therefore cannot say
which tuple is astrophysical and which tuples were injected.

### Boolean CLI flags discard explicit values

Boolean options are registered with `mx::app::argType::True`. The underlying option parser treats the presence of such
an option as true and ignores an attached value, so `--flag=false` currently enables the flag. This caused the positive
recovery preparation command to enable `p4Optimize` accidentally.

Changing these targets to ordinary optional arguments is necessary but not sufficient. A bare optional `--flag` has no
stored value, so it must be interpreted as true through `appConfigurator::isSet()`. It must also override an earlier
`flag=false` from an INI file, while an explicit `--flag=false` must remain the last, authoritative value.

## Configuration contract

### Boolean command-line syntax

Every hciReduce boolean target will support:

```text
--flag             implicit true
--flag=true        explicit true
--flag=false       explicit false
```

INI syntax remains unchanged:

```ini
flag=true
flag=false
```

`--flag false` with a separate token will not be supported; the explicit value must be attached with `=` so it cannot
be mistaken for a positional argument. Command-line values retain precedence over configuration files. Repeating the
same boolean option on one command line is rejected as ambiguous; the implementation does not silently depend on
parser-internal ordering that is not exposed by mxlib.

Implement one shared hciReduce boolean loader rather than open-coding the `isSet()` behavior in each application. The
loader will pair `argType::Optional` registration with these cases:

| Sources | Result |
|---|---|
| target absent | retain the caller's default |
| INI `flag=true/false` | parse the explicit INI value |
| bare `--flag` | true, including over an INI `false` |
| `--flag=true/false` | parse the explicit command-line value |
| malformed attached value | actionable invalid-configuration error |

The first implementation step is a focused parser seam test using the installed mxlib version. If
`appConfigurator`'s public target metadata cannot distinguish an attached value from a bare option with correct source
precedence, make the smallest upstream mxlib change needed to expose that distinction; do not scan or rewrite raw
`argv` independently in each executable.

Apply the contract to all current boolean targets in `hciAnalyze`, `klipReduce`, `p4Reduce`, `HCIobservation`,
`ADIobservation`, `KLIPreduction`, and `P4Reduction`. This includes inverse controls such as
`combine.noDerotate`; only parsing changes, not their existing semantics or defaults. Correct boolean help metadata
while touching the registrations, including the current string-typed help entry for `adi.postMedSub`.

### Known planet metadata

Make the existing `planet` configuration section a shared contract for reduction applications and `hciAnalyze`, and
extend the existing analyzer section with contrast:

```ini
[planet]
sep=11.782
PA=262.051
contrast=0.004878
```

The three values are matching vectors, allowing more than one known companion without a later schema change. In the
reduction applications, if any one is specified, all three must be present, nonempty, finite, and equal length.
Separation and contrast must be nonnegative; PA is finite and retains the existing degrees-east-of-north convention. A
zero contrast is valid metadata but produces no flux when planet subtraction is requested. For backward compatibility,
explicit `hciAnalyze` planet positions may omit contrast; contrast is required in persisted `PLANET*` metadata.

The planet section is metadata only. It does not name a PSF template and does not inject anything by itself. Injection
continues to use `fake.fileName`, `fake.method`, and the optional fake scale file. The separate
`fake.subtractPlanet` flag described below can use that template to subtract the configured planet.

Persist planet tuples independently as comma-separated `PLANETSEP`, `PLANETPA`, and `PLANETCONT` FITS keywords.
Preserve them in ordinary final images, filtered images, local products, optimizer products, and preprocessing-only
images. Existing `FAKE*` keywords remain backward compatible and describe only the explicitly configured fake tuples.

### Subtracting the planet during fake injection

Add the boolean option:

```ini
[fake]
subtractPlanet=false
```

The default `false` preserves existing behavior. When true, the injection stage performs two independent operations:

- inject the configured fake tuples exactly as currently specified; and
- also inject every configured planet tuple with the negative of its `planet.contrast`.

Planet subtraction is a side operation, not composition of the fake list. It must not append, remove, reorder, or
otherwise modify the configured fake vectors. The `FAKESEP`, `FAKEPA`, and `FAKECONT` headers continue to describe only
the explicitly configured fake tuples, exactly as they do when `subtractPlanet=false`. The independent `PLANET*`
headers continue to describe the known planet with its positive physical contrast; they do not represent the negative
injection amplitude.

`subtractPlanet=true` requires complete planet metadata and a usable fake template. It uses the same template
normalization, rotation, interpolation, and optional per-frame scale conventions as ordinary fake injection, changing
only the sign of the configured planet contrast.

The operation must remain idempotent: ordinary full-cube injection subtracts the planet once, and every P4 local,
optimizer, or jackknife evaluation samples one fresh fixed negative-planet source in addition to its configured trial.
The fixed subtraction is not an optimizer parameter and does not count as another configured fake tuple for the
single-trial validation rule.

## `hciAnalyze` resolution contract

Extend the existing explicit `planet.sep` and `planet.PA` inputs with `planet.contrast`, and recognize the `PLANET*`
FITS keywords.

Resolve measurement targets in this order:

1. complete explicit `planet` configuration;
2. complete explicit `fake` configuration;
3. complete `PLANETSEP`, `PLANETPA`, and `PLANETCONT` header metadata; and
4. complete legacy `FAKESEP`, `FAKEPA`, and `FAKECONT` header metadata.

If both `PLANET*` and `FAKE*` occur in a header, the known planet is the default measurement target; users can select
injected sources through explicit `fake` configuration. Continue constructing the noise-exclusion mask from every
resolved target in the selected group.

## Work sequence

### 1. Establish the boolean parsing seam

- Add focused tests for bare, explicit true, explicit false, unset/default, malformed, config-plus-CLI precedence, and
  repeated-option behavior.
- Implement the shared loader and convert one low-risk target such as `showTiming` as the integration proof.
- Confirm generated help communicates the optional boolean value without implying that a following positional token
  is consumed.
- Convert the remaining boolean targets only after the precedence tests pass.

### 2. Add shared planet configuration and fixed subtraction

- Add documented planet state and accessors to `ADIobservation`, following the project's Data/accessor Doxygen
  section ordering.
- Add the default-false `fake.subtractPlanet` option through the shared boolean configuration path.
- Validate planet and fake vector completeness before reading or injecting images.
- Implement planet subtraction as an independent injection operation without mutating the configured fake vectors.
- Preserve existing behavior byte-for-byte when `subtractPlanet=false`.

### 3. Decouple injection from preprocessing

- Replace the `!m_skipPreProcess` injection condition with an explicit full-cube-injection policy.
- Inject after loading for both preprocessing modes; allow configured preprocessing to follow only when it is enabled.
- Add a virtual or otherwise explicit P4 policy hook so `localStampSize>0`, negative optimization, and jackknife paths
  defer injection to `P4TrialSource` and never modify the pristine cube.
- Exercise both `fake.method=single` and the existing per-frame list/scale paths without changing their normalization
  or interpolation conventions.

### 4. Unify FITS provenance

- Add the `PLANET*` keywords to the ADI provenance writer without adding the subtraction to `FAKE*`.
- Introduce one derived-observation provenance hook used by final and preprocessing-only output rather than duplicating
  header assembly.
- Ensure preprocessing-only products now contain `FAKEFILE`, optional scale-file provenance, configured `FAKE*`, and
  independent `PLANET*` metadata.
- Ensure downstream reductions preserve provenance without treating header metadata as a request to inject the same
  source again.

### 5. Extend `hciAnalyze`

- Extend and validate the existing planet namespace with contrast.
- Implement the explicit/header precedence above.
- Report the selected source category in diagnostics and continue printing contrast in the measurement table.
- Add coverage for planet-only, fake-only, combined planet/fake headers, explicit override, and partial vectors.

### 6. Documentation and examples

- Update the common configuration dictionary with the injection-stage rule, planet metadata, and
  `subtractPlanet` examples.
- Update the P4 pixel-local and optimizer guides to state explicitly that their configured trial is deferred even when
  ordinary `preProcess.skip=true` reductions now inject immediately.
- Update the `hciAnalyze` guide and example configuration with planet-source precedence.
- Audit every boolean example. Use `--flag` for implicit true and `--flag=false` only when explicitly disabling a value
  inherited from an INI file; distinguish CLI syntax from `flag=false` INI syntax.

### 7. Verification

- Add Catch2 tests with required Doxygen blocks for the shared boolean loader, ADI configuration/composition,
  skip-preprocessing injection, preprocessing-only provenance, P4 local non-double-injection, and `hciAnalyze`
  resolution.
- Add executable-level tests for at least one boolean in each application, including the exact regression cases
  `--p4Optimize.enabled=false`, `--preProcess.skip=false`, and `--p4Optimize.fitPosition=false`.
- Construct a small synthetic cube showing that full-field injection with `preProcess.skip=true` equals the original
  cube plus one expected rotated source, while a P4 local trial remains identical to the existing on-demand oracle.
- Verify complete planet/fake FITS round trips through preprocessing-only and final outputs.
- Run focused tests, the full CTest suite, application help tests, and `make docs`.
- For every edited function that calls mxlib, inspect the current mxlib LCOV report for each called API as required by
  `AGENTS.md`; record any shortfall in `agents/plans/mxlib_cleanup.md` under `Known non-blocking ownership follow-ups`.

## Completion criteria

- An ordinary KLIP or full-field P4 run can inject configured fakes directly into an already-preprocessed cube with
  `preProcess.skip=true`, without a preparation copy.
- P4 local evaluation and optimization still inject exactly one configured trial per evaluation into a pristine loaded
  cube, plus one fixed negative copy of each planet only when `subtractPlanet=true`.
- Planet and fake tuples remain distinguishable in configuration, every relevant FITS product, and `hciAnalyze`.
- `fake.subtractPlanet=true` subtracts each configured planet exactly once without changing fake vectors or `FAKE*`
  metadata.
- Every boolean CLI target accepts both bare implicit true and attached explicit true/false values with command-line
  precedence over INI files.
- Existing configurations and bare command invocations retain their prior behavior.

## Non-goals

- Do not change fake-template photometric calibration, cubic interpolation, derotation conventions, or per-frame scale
  semantics.
- Do not make one P4 optimizer evaluation support multiple simultaneous trial sources.
- Do not infer a planet from image content or silently copy fake tuples into planet metadata.
- Do not redesign `p4.numberImages` or other temporal-predictor algorithms as part of this cleanup.
