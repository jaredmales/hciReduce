# Step 5: inner-radius rectangular-PSD recovery setup

## Fixed question

Does widening the pooled rectangular-PSD training-center range improve full P4
source recovery at small separations, and where does the larger sample count
stop compensating for radial covariance transfer?

The comparison fixes three covariance models: rectangular PSD with radial
half-widths ±5, ±10, and ±20 pixels, all using raw 11×11 patches, five-pixel
radial/angular center spacing, ensemble-mean subtraction, 21×21 FFT padding,
0.3 flat-spectrum mixing, and a minimum of eight patches. Identity matched
filtering and Gaussian smoothing with FWHM 3.6 pixels are references. Every
method uses the same production annular mean/stddev normalization and
small-sample correction.

## Injection grid

The score-free audit selects six sites at each radius and uses three fixed
brightness levels, for **108 full P4 positive reductions**.

| Radius (pixels) | Sites | ±5 valid | ±10 valid | ±20 valid | Minimum site spacing |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 6 | 6 | 0 | 0 | 6 | 1.00 px |
| 8 | 6 | 0 | 4 | 6 | 2.83 px |
| 12 | 6 | 5 | 6 | 6 | 6.00 px |
| 16 | 6 | 6 | 6 | 6 | 10.77 px |
| 20 | 6 | 6 | 6 | 6 | 15.81 px |
| 24 | 6 | 6 | 6 | 6 | 19.24 px |

Selection requires a complete ±20 five-pixel search and never reads a score or
positive image. Invalid narrower searches remain nondetections. The radius-6
sites are close because the known-planet guard leaves only ten eligible native
centers; their recovery counts are strongly correlated.

## Known planet and holdouts

The known planet is centered at `(139.169, 125.871)` in native `(x, y)` pixels.
The setup buffers its stored 7.3-pixel exclusion by the production half pixel.
For every injection, calibration search, and annular-noise pixel, the full
15×15 Gaussian kernel must be disjoint from that circle. Since the matched
filter uses an 11×11 stamp, this also keeps its candidate data clear of the
planet. Every covariance fit excludes the buffered planet circle and the
current trial's complete five-search 15×15 kernel union.

No pixel inside the configured planet exclusion therefore enters candidate
inputs, covariance training, or annular noise estimation. Residual planet wings
beyond that configured guard are not modeled. Other calibration and evaluation
footprints are not excluded as a global union. Each search receives its own
holdout so inner-ring training is not erased by the many overlapping trials.

## Calibration and source scale

Each active radius/method uses 20 geometry-selected planet-clean baseline
searches. The script chooses the narrowest radial band with complete support:

| Evaluation radius | Calibration radii | Notes |
| ---: | ---: | --- |
| 6 | 5–7 | ±20, identity, and Gaussian only; ±5/±10 are invalid |
| 8 | 7–9 | ±10, ±20, identity, and Gaussian; ±5 is invalid |
| 12, 16, 20, 24 | Same ring | All methods active where the evaluation grid has support |

The threshold is the maximum of the 20 five-pixel null search scores, with a
strict exceedance rule. All thresholds and contrasts are written and
fingerprinted before positive reductions start. The three injection contrasts
are 0.5, 0.75, and 1.0 times an identity-based reference: identity's threshold
times its interpolated amplitude-map radial standard deviation, divided by the
small-sample correction. Target-site pixel values do not set the contrast.
Injected-source neighborhoods remain in each positive image's annular mean and
standard-deviation estimate, matching the existing production path.

## ROC execution

After pulling this commit on ROC, run setup from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/p4-step5-matplotlib \
  /opt/conda/envs/xpy3_13/bin/python3 \
  agents/plans/scripts/run_p4_step5_inner_rectangular.py setup \
  --parent working/roc/p4_psd_full_20260918 \
  --root working/roc/p4_inner_rectangular_20260919
```

Setup verifies the completed parent study, all raw inputs, and frozen software;
then it copies the native executables and dependencies and writes the immutable
geometry. Start the resumable run in `tmux`:

```sh
tmux new-session -d -s p4-inner-rect \
  "cd /home/jrmales/Source/mxApps/hciReduce && \
   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/p4-step5-matplotlib \
   /opt/conda/envs/xpy3_13/bin/python3 \
   working/roc/p4_inner_rectangular_20260919/software/run_p4_step5_inner_rectangular.py run \
   --root working/roc/p4_inner_rectangular_20260919 \
   > working/roc/p4_inner_rectangular_20260919/driver.log 2>&1"
```

Monitor without changing the run:

```sh
cat working/roc/p4_inner_rectangular_20260919/state.json
tail -n 30 working/roc/p4_inner_rectangular_20260919/driver.log
```

The run first performs all per-search baseline analyses and freezes
`thresholds.json` and `jobs.json`. It then executes the 108 full reductions in
order and can resume completed jobs. Final products are `results.json`,
`results.md`, and `comparison.png`.

## Pre-calibration runner repair

The first launch stopped on the first baseline search before creating
thresholds, jobs, or positive reductions. Its sparse cube contained an
unsupported ±5 plane with too few finite radial-profile points for the
production GSL interpolator. The valid ±10, ±20, identity, and Gaussian planes
from that same saved cube completed successfully when replayed alone.

The runner now sends only methods with complete five-pixel support through
`hciAnalyze` and restores unsupported methods as `NaN` nondetections. Each
active plane is still checked against the independent annular-SNR oracle. To
repair the existing pre-calibration setup after pulling the fix, run:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/p4-step5-matplotlib \
  /opt/conda/envs/xpy3_13/bin/python3 \
  agents/plans/scripts/run_p4_step5_inner_rectangular.py repair \
  --root working/roc/p4_inner_rectangular_20260919
```

The repair verifies every unchanged frozen file and raw input, archives the
partial failed analysis under a numbered `pre_repair_failures/repair_NNNN/`
directory, replaces and rehashes only the frozen runner, and appends a
versioned repair record. Earlier repair records remain immutable and are
verified before each update. It refuses to run after calibration or any
positive reduction. Restart the same frozen run command, using a new log name
such as `driver-repaired.log`.

An analysis stopped before serialization completes may have output files but
no `complete.json` receipt. On resume, the runner preserves such a directory under
`interrupted_analysis/<analysis>/attempt_NNNN/` and recomputes only that
analysis. Directories with verified completion receipts are reused. This
allows calibration and later positive analysis to resume after a process
interruption without overwriting diagnostic artifacts or repeating completed
work.

For the existing 2026-09-19 root, the first repair remains recorded in the
legacy `repair.json`, and the resume correction remains in
`repair_0002.json`. The underlying failure at `baseline__null_x119_y124` was
then isolated to its unsupported ±5 method: validity, search score, and
nondetection handling were correct, but the five unavailable per-pixel
diagnostics were still represented as floating-point `NaN`. The strict JSON
writer rejected them. They are now recorded as JSON `null`, while the method
remains invalid with no score. Applying the next repair writes
`repair_0003.json`, archives only that unreceipted search under
`pre_repair_failures/repair_0003/`, and retains the 45 baseline analyses with
verified receipts.

## Preparation validation

The local `audit` action completed against the saved parent baseline and exact
response templates: 36 sites, 128 unique calibration centers, all six ±20
searches valid at every radius, and the narrower support counts shown above.
All selected valid center and neighbor PSD fits factored successfully. Python
syntax compilation and repository whitespace checks passed. No new reduction
was run, and no production or mxlib-calling function changed.
