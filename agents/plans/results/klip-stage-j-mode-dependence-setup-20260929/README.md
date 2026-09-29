# KLIP Stage-J mode-dependence setup

## Purpose

Stage I fixed the KLIP mode at fraction 0.200. The planet and individual
injection sites can fluctuate differently as the retained KLIP fraction
changes, so one fixed fraction may either conceal or exaggerate the
planet-versus-injection comparison. Stage J measures the complete frozen mode
fraction grid:

`0.125, 0.150, 0.175, 0.200, 0.225, 0.250, 0.300, 0.350`.

No new KLIP reductions are required. The original planet cube and the twelve
nominal optimized-phase Stage-G injection cubes already contain every mode.

## Frozen analysis

For every mode, task, and method, Stage J reconstructs the Stage-I
optimized-nearest-pixel statistic at the common subpixel phase. It recalibrates
cubic-shifted identity, raw rectangular covariance, and radial-Hann covariance
weights using the corresponding mode's signal-free baseline and radius-12
calibration unit. Gaussian 2.4 and 3.6 and all integer-template controls are
retained.

Production `hciAnalyze` normalizes all 64 mode-by-method maps in one call per
science image. An independent annular oracle must reproduce every finite SNR
pixel. All unchanged Gaussian and integer-template injection values must
reproduce the completed Stage-G result at all eight modes.

The output retains every per-site curve in `mode_curves.csv` and
`results.json`. It reports three complementary comparisons:

1. **Modewise:** compare the planet with the twelve injections independently at
   every fixed fraction.
2. **Independent mode scan:** take the maximum over modes for every planet and
   injection curve before comparing them. This applies the same look-elsewhere
   operation at every location.
3. **Injection-selected mode:** select the planet's mode using the mean of all
   twelve injections. Each injection is evaluated at the mode selected by the
   other eleven sites, giving a leave-one-site-out estimate.

Method families choose their modes independently in the scanned and
injection-selected comparisons. The planet never selects the mode used for an
inferential comparison.

## Questions answered

The primary outputs are the raw- and radial-covariance minus shifted-identity
SNR distributions at every mode, after independent scanning, and after
injection-only mode selection. Secondary outputs include absolute SNR curves,
within-site mode standard deviations and ranges, shifted-minus-integer effects,
and histograms of the maximizing mode across injection sites.

This design distinguishes a general covariance penalty from a mode-specific
crossing and shows whether the planet's Stage-I extreme rank persists once
site-dependent mode fluctuations are included.

## ROC commands

After pulling the setup commit, run from the repository root:

```bash
root=working/roc/klip_stage_c_development_20260925
python3 agents/plans/scripts/run_klip_stage_j_mode_dependence.py check

taskset -c 12-27 python3 \
  agents/plans/scripts/run_klip_stage_j_mode_dependence.py prepare "$root"
```

Then launch the frozen command printed by `prepare`:

```bash
root=working/roc/klip_stage_c_development_20260925
taskset -c 12-27 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 python3 \
  "$root/stage_j_mode_dependence/software/run_klip_stage_j_mode_dependence.py" \
  run "$root" \
  > "$root/stage_j_mode_dependence/driver.log" 2>&1
```

The run is resumable. Monitor it with:

```bash
tail -f working/roc/klip_stage_c_development_20260925/stage_j_mode_dependence/driver.log
```

The primary report will be
`working/roc/klip_stage_c_development_20260925/stage_j_mode_dependence/results.md`.
