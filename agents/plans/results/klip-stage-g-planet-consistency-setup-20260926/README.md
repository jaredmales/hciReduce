# KLIP Stage-G planet-consistency injection setup

## Purpose

Stage G directly tests the mismatch identified after the known-planet closure.
The six comparable Stage-E validation positions do not make the planet an
absolute-SNR outlier, but the planet's covariance penalty relative to exact
identity is below all six injections. Stage G increases the angular sample,
uses the exact optimized contrast rather than interpolation, and applies the
same radius-three aperture statistic as Stage F.

The campaign also addresses the concern that the injection PSF may not model
the real planet. It separates subpixel registration from controlled PSF
broadening without tuning either control on the observed planet SNR.

## Frozen design

The runner selects 12 radius-12 centers by deterministic angular maximin from
the Stage-C calibration and held-out-null roles. These centers were never used
for a positive injection, and selection reads geometry only. Each center is
run in five arms, for 60 reductions total:

| Arm | Injection PSF | Subpixel phase |
| :--- | :--- | :--- |
| `nominal_integer` | Original `psf_reg_median.fits` | Pixel centered |
| `nominal_analysis_phase` | Original | Stage-F analysis phase: (+0.1688, -0.1294) pixels |
| `nominal_optimized_phase` | Original | Optimized negative-fit phase: (-0.2771, +0.4860) pixels |
| `blur0p9_optimized_phase` | Additional Gaussian FWHM 0.9 pixels | Optimized phase |
| `blur1p8_optimized_phase` | Additional Gaussian FWHM 1.8 pixels | Optimized phase |

The two broadenings are 0.25 and 0.5 lambda/D. They preserve the original PSF
sum and do not attempt to select a best planet model. They are controlled
proxies for chromatic effective-PSF changes, temporal AO variation, rotational
or registration smearing, and related morphology errors. They cannot identify
which physical effect is present.

Every injection uses contrast `0.004574362496414845`, the independently
optimized negative-companion result. The positive reductions retain all eight
KL modes, but mode 200 is the preregistered endpoint. The primary arm is
`nominal_optimized_phase`, and the primary paired quantities are:

- raw rectangular covariance minus exact identity SNR;
- radial-Hann covariance minus exact identity SNR; and
- both nearest-pixel and radius-three aperture-maximum forms.

Gaussian 3.6 versus exact identity, Gaussian 3.6 versus Gaussian 2.4, phase
effects, PSF-broadening effects, and finite-response fidelity are secondary
reported controls.

## Analysis contract

For each injection, the runner rebuilds the complete Stage-F aperture maps.
Weights come from the signal-free baseline, with the optimized-planet disk and
the union of every 11-by-11 aperture footprint excluded from training. The
source remains in candidate support. Generic frozen maps supply annular
normalization, while both the optimized planet and the injected source are
excluded from the annular noise profile.

The aperture is centered on the exact floating injection location and uses the
same Stage-F rule: native pixel centers within `snr.apertureR + 0.5`, with
`snr.apertureR=3`. Production `hciAnalyze` SNR must reproduce the independent
annular oracle for every map. The finite positive response at mode 200 is also
compared with the nominal exact response under identity, raw-covariance, and
radial-covariance metrics.

The Gaussian PSF controls use the same PSF file for fake-source injection and
optimized-planet subtraction because that is the existing `klipReduce`
interface. The known-planet region is excluded from training and annular noise,
but this remains a limitation when interpreting the broadened arms.

## Resumability and provenance

`prepare` verifies the completed Stage-E and Stage-F boundaries, all calibration
model receipts, the source binaries, configuration, response products, and
planet result. It creates and fingerprints the three PSF files, chooses sites,
preflights every distinct aperture at mode 200, freezes all 60 commands, and
copies the runner plus the exact Stage-F implementation.

`run` is resumable at both reduction and analysis granularity. Incomplete task
directories are archived before replay. A final receipt covers every nested
reduction and analysis receipt plus the aggregate JSON, CSV, and Markdown
reports.

## ROC commands

From the repository root after pulling the setup commit:

```bash
root=working/roc/klip_stage_c_development_20260925
python3 agents/plans/scripts/run_klip_stage_g_planet_consistency.py check \
  --config working/analyze.conf
taskset -c 12-27 python3 \
  agents/plans/scripts/run_klip_stage_g_planet_consistency.py prepare "$root" \
  --config working/analyze.conf
```

Then run the frozen copy printed by `prepare`:

```bash
root=working/roc/klip_stage_c_development_20260925
taskset -c 12-27 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 python3 \
  "$root/stage_g_planet_consistency/software/run_klip_stage_g_planet_consistency.py" \
  run "$root" --workers 4 \
  > "$root/stage_g_planet_consistency/driver.log" 2>&1
```

The command can be launched directly and left to completion. Progress is in:

```bash
tail -f working/roc/klip_stage_c_development_20260925/stage_g_planet_consistency/driver.log
```

The primary report will be
`working/roc/klip_stage_c_development_20260925/stage_g_planet_consistency/results.md`.

## Initial analysis repair (2026-09-27)

The first ROC run completed all 60 reductions, then failed during analysis.
`hciAnalyze` uses every `planet.*` entry both as an annular exclusion and as a
signal to score. The real-planet entry is needed as an exclusion, but the
Stage-G maps intentionally have no valid pixels in its three-pixel aperture.
The program therefore reported `signal 0 aperture contains no valid pixels`.

The repaired runner gives `hciAnalyze` a 60-pixel reporting aperture, as in the
validated Stage-E workflow. This setting only controls the program's summary
rows; it does not change the annular SNR map. Stage G continues to calculate
all injection endpoints over the reviewed three-pixel aperture. The repair
archives the original protocol, manifest, and runner, updates their
fingerprints, and requires all 60 reduction receipts before proceeding.

After pulling the repair commit on ROC, upgrade the prepared directory once:

```bash
root=working/roc/klip_stage_c_development_20260925
python3 agents/plans/scripts/run_klip_stage_g_planet_consistency.py \
  repair-analysis-aperture "$root"
```

Then rerun the frozen command. The runner reuses all completed reductions and
archives the incomplete analysis directories before rebuilding them:

```bash
taskset -c 12-27 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 python3 \
  "$root/stage_g_planet_consistency/software/run_klip_stage_g_planet_consistency.py" \
  run "$root" --workers 4 \
  > "$root/stage_g_planet_consistency/driver.log" 2>&1
```


## All-mode site repair (2026-09-27)

The repaired analysis completed 38 task receipts before exposing a second
preflight gap. The initial setup preflight exercised the complete aperture only
at mode 200. Five of the original twelve centers lose the required eight
training patches in one detector half at mode 125 for at least one phase:
`r12_cal03`, `r12_cal05`, `r12_cal07`, `r12_cal14`, and
`r12_cal15`. This is a geometry/support failure, not an injection-score
result.

A score-blind audit tested all 26 eligible radius-12 calibration and held-out
null centers at every phase, KL mode, and aperture pixel. Seventeen centers
have complete detector-half support. The repair retains the seven fully
supported original sites and fills the five open slots by the original
deterministic angular maximin rule, conditioned on those retained sites. The
replacements are `r12_hold02`, `r12_cal16`, `r12_cal13`,
`r12_cal02`, and `r12_cal18`. Each replacement then passed the
complete Stage-F amplitude-map construction for all three phases and all eight
KL modes.

The repair preserves 35 completed reductions and analyses from the seven valid
sites. It archives all products from the five replaced slots, including the
three completed optimized-phase analyses at `r12_cal07`, and queues 25
replacement reductions and analyses. No recovered SNR was used for support
screening or replacement selection.

After pulling the all-mode repair commit on ROC, upgrade the prepared directory
once:

```bash
root=working/roc/klip_stage_c_development_20260925
taskset -c 12-27 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 python3 \
  agents/plans/scripts/run_klip_stage_g_planet_consistency.py \
  repair-all-mode-sites "$root"
```

Then resume with the frozen runner:

```bash
taskset -c 12-27 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 python3 \
  "$root/stage_g_planet_consistency/software/run_klip_stage_g_planet_consistency.py" \
  run "$root" --workers 4 \
  > "$root/stage_g_planet_consistency/driver.log" 2>&1
```
