# KLIP Stage-D policy and Stage-E validation setup

## Status

The policy-freeze and validation runners are ready for ROC execution against
`working/roc/klip_stage_c_development_20260925`. A read-only preflight on
2026-09-26 confirmed that all 2,880 validation and held-out query-mode
combinations (360 distinct spatial queries evaluated in each of eight modes)
are present in the frozen generic Stage-C annular models. The experiment remains in `development_complete`: no held-out
score or validation positive has been opened.

The maintained runners are:

- [`freeze_klip_stage_d_policy.py`](../../scripts/freeze_klip_stage_d_policy.py); and
- [`run_klip_stage_e_validation.py`](../../scripts/run_klip_stage_e_validation.py).

## Immutable Stage-D policy

The freeze step selects the seven-method shortlist already recommended by the
completed Stage-C development report:

1. native image;
2. Gaussian FWHM 2.4;
3. Gaussian FWHM 3.6;
4. exact-response identity;
5. sparse-response identity;
6. raw rectangular Welch PSD with mixing 0.3; and
7. radial-standardized Hann Welch PSD with mixing 0.1 and hard precision
   truncation below 0.75 times the mean variance.

Mode 200 is primary. The other seven modes are reported separately and are not
maximized. All response and covariance methods use baseline-derived weights;
the development per-positive refit arm is retired.

For every selected method, radius, and mode, the receipt copies the existing
maximum-of-20 threshold without recalculation and records all five signed
scores at every calibration site. It also records the largest score, second
largest score, their gap, and angular geometry. Threshold sensitivity uses
4,096 deterministic paired bootstrap resamples shared across methods and
modes at a radius. A second diagnostic deletes four consecutive calibration
sites in angular order, wrapping around the annulus, for all 20 possible
blocks. These samples are spatially correlated and are not assigned a formal
false-alarm probability.

The policy fingerprints the complete calibration units and sites, all 108
development reductions and analyses, the reduction inputs and responses, both
new runners, all frozen commands, and all generated policy artifacts. Stage E
refuses to run if any fingerprint changes.

## Frozen validation gate

The two covariance candidates are evaluated against Gaussian FWHM 3.6 at mode
200. A candidate is accepted only when all of these conditions hold:

- its aggregate held-out-null exceedance count across the six radii is no
  larger than Gaussian 3.6;
- its recovery is at least Gaussian 3.6 at each of SNR 3, 5, and 7, and its
  total recovery is strictly larger;
- a 4,096-replicate paired, radius-stratified site bootstrap has a 95-percent
  lower bound above zero for mean maximum-SNR change at SNR 3 or 5, with no
  level whose 95-percent upper bound is below zero;
- pooled mean positive-minus-baseline throughput is in `[0.90, 1.10]` at every
  level and every radius-level mean is in `[0.85, 1.15]`; and
- candidate and comparator retain all 36 common sites at every level.

Gaussian FWHM 2.4 is reported as the stronger smoothing control. Exact and
sparse identity remain eligible scientific outcomes, but they do not count as
covariance successes.

## Execution order

Use the CPU affinity frozen by the parent experiment:

```bash
root=working/roc/klip_stage_c_development_20260925

taskset -c 12-27 python3 agents/plans/scripts/freeze_klip_stage_d_policy.py check
taskset -c 12-27 python3 agents/plans/scripts/run_klip_stage_e_validation.py check

taskset -c 12-27 python3 agents/plans/scripts/freeze_klip_stage_d_policy.py freeze "$root"
cat "$root/policy/README.md"
```

The freeze command reads development and calibration products only. After the
receipt has been reviewed, Stage E can be run as one resumable command:

```bash
taskset -c 12-27 python3 agents/plans/scripts/run_klip_stage_e_validation.py run "$root" --workers 4
```

For checkpoints, use the same frozen runner in this order:

```bash
taskset -c 12-27 python3 agents/plans/scripts/run_klip_stage_e_validation.py models "$root" --workers 4
taskset -c 12-27 python3 agents/plans/scripts/run_klip_stage_e_validation.py heldout "$root" --workers 4
taskset -c 12-27 python3 agents/plans/scripts/run_klip_stage_e_validation.py reduce "$root"
taskset -c 12-27 python3 agents/plans/scripts/run_klip_stage_e_validation.py analyze "$root" --workers 4
```

The repository entry point verifies and then executes the copy frozen in the
experiment. The first Stage-E launch exposed a driver-only unpacking error
before `models` began. The guarded repair is permitted only while the state is
`policy_frozen` and no model, held-out, or validation directory exists. It
archives the original runner and policy receipts under
`policy_repairs/stage_e_entry_unpack_20260926`, replaces only the runner,
updates its software fingerprint and repair receipt, and leaves every policy
artifact byte-for-byte unchanged.

`models` uses only the signal-free baseline and can run before any score is
opened. `heldout` exposes the 36 preassigned null sites. `reduce` then creates
the 108 validation positives, and `analyze` applies the frozen weights and
thresholds. Every phase verifies and reuses completed receipts, and archives
an incomplete task directory before retry.

## Validation performed before launch

- Both scripts compile under Python 3.13.
- Both deterministic self-checks pass.
- Git whitespace validation passes.
- The ROC preflight found zero missing validation or held-out generic queries.
- ROC still had no `heldout_analysis`, `validation_reductions`, or
  `validation_analysis` directory at the time of the preflight.

The full numerical policy generation and real estimator replay require the ROC
experiment tree and remain the next execution checkpoint.
