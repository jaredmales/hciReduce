# Step 5: inner-radius SNR 3/5/7 confirmation setup

**Completion note:** all 108 reductions completed, but the initial annular SNR
included each injection in its own noise sample. The
[completion and diagnosis](../p4-step5-inner-snr357-20260919/README.md) retain
that result and define a corrected saved-image analysis with no new reductions.

## Fixed source scale

This confirmatory study repeats the inner-radius comparison with nominal
source-only annular SNRs **3, 5, and 7**. The SNR-5 level is `1×`; the stored
brightness multipliers are therefore `0.6×`, `1×`, and `1.4×`.

For each site, let \(\sigma_\alpha\) be the identity matched-filter amplitude
map's interpolated baseline radial standard deviation and let \(f_{\rm small}\)
be the production small-sample multiplier. The contrast for nominal source SNR
\(S\) is frozen as

\[
c(S)=\frac{S\,\sigma_\alpha}{f_{\rm small}}.
\]

This is the source's expected increment in the same annular-SNR statistic used
by `hciAnalyze`. It does not use the candidate's baseline amplitude or any
positive image. The actual positive-image SNR can differ from the target
because of the local background, five-pixel maximum, changes to the annular
profile, and finite-source reduction nonlinearity. The report retains the mean
measured production search-maximum and fixed-center SNR at every radius and
level. The fixed-center mean directly checks the nominal 3/5/7 scale; the
five-pixel maximum is the detection statistic.

Identity is the fixed scale reference because its response-normalized amplitude
is in injected-contrast units. Every filter receives the same injected source;
the source is not rescaled separately for each method. The preceding paired
increment check found approximately unit response, so the baseline-derived
scale is suitable for freezing before positive reductions.

The earlier inner-radius studies used a null-threshold scale whose `1×` mean
SNR was generally near 2. They remain exploratory low-SNR comparisons and are
not the confirmatory SNR-5 result.

Applying the new formula to the completed study's frozen identity radial-noise
profiles predicts the following SNR-5 (`1×`) contrasts. Setup recalculates and
freezes them from the new root's verified baseline analyses before any positive
reduction.

| Radius (pixels) | Minimum `1×` contrast | Median `1×` contrast | Maximum `1×` contrast |
| ---: | ---: | ---: | ---: |
| 6 | `2.9251e-2` | `3.5031e-2` | `3.9159e-2` |
| 8 | `1.7395e-2` | `1.7395e-2` | `1.9830e-2` |
| 12 | `6.7705e-3` | `7.1493e-3` | `7.3329e-3` |
| 16 | `2.5187e-3` | `2.9848e-3` | `3.2378e-3` |
| 20 | `1.3360e-3` | `1.3781e-3` | `1.4144e-3` |
| 24 | `1.5275e-3` | `1.5420e-3` | `1.5500e-3` |

The SNR-3 and SNR-7 contrasts are 0.6 and 1.4 times these values. The large
inner contrasts are intentional: they test comparable nominal source strength
where the residual noise is much larger. They may also expose finite-source
nonlinearity, which is why measured SNR and contrast recovery remain required
outputs rather than assuming exact linear scaling.

## Fixed comparison

The study retains the completed inner experiment's geometry and policies:

- six sites at each radius 6, 8, 12, 16, 20, and 24 pixels;
- rectangular PSD radial half-widths ±5, ±10, and ±20 pixels;
- identity matched filtering and Gaussian FWHM 3.6 references;
- the existing known-planet guard, per-trial holdout, 11×11 stamps, five-pixel
  search, production annular normalization, and method-specific null
  thresholds;
- 20 geometry-selected baseline searches per active radius and method.

There are **108 new full P4 reductions**: 36 sites at three source levels. All
thresholds and SNR-derived contrasts are frozen from the baseline before the
first positive reduction. A second stage reuses those 108 images to compare
raw and post-mean patch-RMS rectangular PSD estimates at all three widths, with
zero additional P4 reductions.

## ROC setup and execution

After pulling the setup commit on ROC, prepare the full-reduction study from
the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/p4-step5-matplotlib \
  /opt/conda/envs/xpy3_13/bin/python3 \
  agents/plans/scripts/run_p4_step5_inner_rectangular.py setup \
  --parent working/roc/p4_psd_full_20260918 \
  --root working/roc/p4_inner_snr357_20260919 \
  --target-snrs 3 5 7
```

Start the resumable reductions and raw/reference analysis:

```sh
tmux new-session -d -s p4-inner-snr357 \
  "cd /home/jrmales/Source/mxApps/hciReduce && \
   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/p4-step5-matplotlib \
   /opt/conda/envs/xpy3_13/bin/python3 \
   working/roc/p4_inner_snr357_20260919/software/run_p4_step5_inner_rectangular.py run \
   --root working/roc/p4_inner_snr357_20260919 \
   > working/roc/p4_inner_snr357_20260919/driver.log 2>&1"
```

Monitor the first stage without changing it:

```sh
cat working/roc/p4_inner_snr357_20260919/state.json
tail -n 30 working/roc/p4_inner_snr357_20260919/driver.log
```

After that state is `complete`, prepare the saved-image patch-RMS comparison
with the current trial excluded through one effective lambda/D. This uses a new
root; the earlier 7.3-pixel-mask root is a retained failed design:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/p4-step5-matplotlib \
  /opt/conda/envs/xpy3_13/bin/python3 \
  agents/plans/scripts/compare_p4_step5_inner_patch_rms.py prepare \
  --study working/roc/p4_inner_snr357_20260919 \
  --root working/roc/p4_inner_snr357_patch_rms_lambdad_20260919 \
  --exclude-trial-from-annular \
  --cpus 0 1 2 3 4 5 6 7 8 9 10 11
```

Run that analysis in a second resumable session:

```sh
tmux new-session -d -s p4-inner-snr357-rms-lambdad \
  "cd /home/jrmales/Source/mxApps/hciReduce && \
   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/p4-step5-matplotlib \
   /opt/conda/envs/xpy3_13/bin/python3 \
   agents/plans/scripts/compare_p4_step5_inner_patch_rms.py run \
   --root working/roc/p4_inner_snr357_patch_rms_lambdad_20260919 \
   > working/roc/p4_inner_snr357_patch_rms_lambdad_20260919/driver.log 2>&1"
```

Monitor the saved-image analysis:

```sh
cat working/roc/p4_inner_snr357_patch_rms_lambdad_20260919/state.json
tail -n 30 working/roc/p4_inner_snr357_patch_rms_lambdad_20260919/driver.log
```

The failed `p4_inner_snr357_patch_rms_excluded_20260919` root records both the
strict-JSON diagnostic repair and the subsequent geometry failure at 140 of
164 baseline receipts. Its 7.3-pixel trial mask cannot support six injection
sites and is not reused. The replacement configures `planet.R=3.1`; production's
half-pixel boundary gives an effective 3.6-pixel or one-lambda/D trial mask.
The completed SNR-3/5/7 parent remains immutable.

## Post-summary runner repair

The first SNR-3/5/7 launch completed and preserved all 108 reductions and
measurements, then stopped while constructing the new fixed-center SNR summary.
The base runner records the five SNR values as `pixels`; the summary used the
`snr_pixels` name from the later patch-RMS schema. No final result product was
written.

After pulling the repair commit, update only the frozen runner:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/p4-step5-matplotlib \
  /opt/conda/envs/xpy3_13/bin/python3 \
  agents/plans/scripts/run_p4_step5_inner_rectangular.py repair \
  --root working/roc/p4_inner_snr357_20260919
```

This repair requires the exact `KeyError`, completed calibration, all 108 jobs
in frozen order, verified products for every measurement, unchanged inputs,
and no final result files. It replaces and rehashes only the frozen runner and
records every verified measurement in the repair receipt. Restart the same run
command; it reuses all completed reductions and analyses and regenerates only
the summaries.

## Local validation

The SNR-target mode passed Python syntax checks and a complete score-free audit
against the archived parent study. The audit reproduced 36 sites, 128 unique
calibration centers, and 108 positive jobs, and recorded target SNRs 3/5/7 as
brightness multipliers 0.6/1/1.4. The historical no-option mode retains its
original threshold-scaled protocol for reproducibility.
