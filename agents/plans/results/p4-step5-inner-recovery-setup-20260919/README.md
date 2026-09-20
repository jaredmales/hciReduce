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

The multiplier is specific to this protocol and site. It does not preserve the
preceding full ROC study's roughly `4e-4` scale. The inner-study median `1×`
contrasts at radii 6, 8, 12, 16, 20, and 24 pixels are `1.2937e-2`,
`6.4241e-3`, `2.4259e-3`, `1.5371e-3`, `7.7270e-4`, and `8.8063e-4`,
respectively. The accepted AF Lep b negative-injection coefficient is
`4.7639e-3`; at its 11.8-pixel separation, the radius-12 `1×` injections are
therefore about half its brightness rather than planet-matched injections.

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

The following launch reached 93 completed baseline searches before the ±10
method at `null_x123_y122` exposed an annular-normalization edge case. Its five
covariance fits were supported, but the innermost search pixel had no lower
bracketing radial-noise bin. The oracle returned an undefined value there,
while production converted the non-finite result to zero. The runner now
requires finite annular normalization at all five search pixels before a
method enters `hciAnalyze`; otherwise that method/search is an invalid
nondetection. The next repair is `repair_0004.json` and retains the 93 complete
searches.

Because this trial is in the frozen radius-8 ±10 calibration pool, the runner
replaces it before thresholding with a valid candidate from the already frozen
null union in the same radial band. It chooses the candidate maximizing the
minimum distance from retained trials and never reads a score during this
selection. The effective 20-search pools, rejected original trials, and
replacements are frozen in `effective_calibration_pools.json` before positive
reductions. The available replacement for this pool is `null_x121_y121`,
which already has a verified analysis receipt.

## Completed ROC result

The study completed all 164 baseline searches and 108 positive reductions.
The completion receipt records unchanged frozen inputs. A final audit found
108 jobs, 108 measurement records, 108 reduction receipts, and 272 analysis
receipts (164 baseline plus 108 positive), with identical unique job-name
sets. It reverified 49 frozen records, four repair records, four calibration
products, and six final products. The radius-8 ±10 replacement described
above was the only calibration-pool adjustment.

Recovery uses all six sites as the denominator. Invalid searches are
nondetections.

| Radius | Method | Valid sites | 0.5× | 0.75× | 1× | Held-out nulls |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 6 | Rectangular ±5 | 0/6 | 0/6 | 0/6 | 0/6 | 0/6 |
| 6 | Rectangular ±10 | 0/6 | 0/6 | 0/6 | 0/6 | 0/6 |
| 6 | Rectangular ±20 | 6/6 | 2/6 | 2/6 | 2/6 | 0/6 |
| 6 | Identity | 6/6 | 2/6 | 2/6 | 3/6 | 0/6 |
| 6 | Gaussian FWHM 3.6 | 6/6 | 1/6 | 3/6 | 3/6 | 0/6 |
| 8 | Rectangular ±5 | 0/6 | 0/6 | 0/6 | 0/6 | 0/6 |
| 8 | Rectangular ±10 | 3/6 | 1/6 | 2/6 | 2/6 | 0/6 |
| 8 | Rectangular ±20 | 6/6 | 1/6 | 2/6 | 3/6 | 0/6 |
| 8 | Identity | 6/6 | 1/6 | 3/6 | 3/6 | 0/6 |
| 8 | Gaussian FWHM 3.6 | 6/6 | 5/6 | 5/6 | 6/6 | 0/6 |
| 12 | Rectangular ±5 | 4/6 | 3/6 | 4/6 | 4/6 | 0/6 |
| 12 | Rectangular ±10 | 6/6 | 4/6 | 6/6 | 6/6 | 0/6 |
| 12 | Rectangular ±20 | 6/6 | 4/6 | 5/6 | 6/6 | 1/6 |
| 12 | Identity | 6/6 | 4/6 | 5/6 | 6/6 | 0/6 |
| 12 | Gaussian FWHM 3.6 | 6/6 | 2/6 | 3/6 | 5/6 | 0/6 |
| 16 | Rectangular ±5 | 6/6 | 0/6 | 0/6 | 1/6 | 0/6 |
| 16 | Rectangular ±10 | 6/6 | 0/6 | 0/6 | 2/6 | 0/6 |
| 16 | Rectangular ±20 | 6/6 | 0/6 | 0/6 | 1/6 | 0/6 |
| 16 | Identity | 6/6 | 0/6 | 1/6 | 3/6 | 0/6 |
| 16 | Gaussian FWHM 3.6 | 6/6 | 0/6 | 1/6 | 1/6 | 0/6 |
| 20 | Rectangular ±5 | 6/6 | 1/6 | 3/6 | 3/6 | 0/6 |
| 20 | Rectangular ±10 | 6/6 | 1/6 | 2/6 | 3/6 | 0/6 |
| 20 | Rectangular ±20 | 6/6 | 0/6 | 1/6 | 3/6 | 0/6 |
| 20 | Identity | 6/6 | 1/6 | 1/6 | 4/6 | 0/6 |
| 20 | Gaussian FWHM 3.6 | 6/6 | 1/6 | 2/6 | 3/6 | 0/6 |
| 24 | Rectangular ±5 | 6/6 | 0/6 | 1/6 | 3/6 | 0/6 |
| 24 | Rectangular ±10 | 6/6 | 0/6 | 0/6 | 3/6 | 0/6 |
| 24 | Rectangular ±20 | 6/6 | 0/6 | 1/6 | 3/6 | 0/6 |
| 24 | Identity | 6/6 | 0/6 | 1/6 | 3/6 | 0/6 |
| 24 | Gaussian FWHM 3.6 | 6/6 | 0/6 | 2/6 | 2/6 | 0/6 |

Across all 36 sites, the aggregate 0.5×/0.75×/1× recoveries are:

| Method | Valid sites | 0.5× | 0.75× | 1× | Held-out nulls /36 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Rectangular ±5 | 22/36 | 4 | 8 | 11 | 0 |
| Rectangular ±10 | 27/36 | 6 | 10 | 16 | 0 |
| Rectangular ±20 | 36/36 | 7 | 11 | 18 | 1 |
| Identity | 36/36 | 8 | 13 | 22 | 0 |
| Gaussian FWHM 3.6 | 36/36 | 9 | 16 | 20 | 0 |

Wider pooling achieves the intended coverage improvement, but ±20 does not
beat the fully supported references overall. Gaussian is strongest at radius
8, recovering 5/6, 5/6, and 6/6 versus ±20's 1/6, 2/6, and 3/6. Radius-12
±10 is the strongest localized covariance result at 4/6, 6/6, and 6/6, but
the six sites are correlated. The evidence does not justify pooled rectangular
PSD as the production small-separation default.

Final product fingerprints include:

| Product | Bytes | SHA-256 |
| --- | ---: | --- |
| `results.json` | 527308 | `26d0dc578cb7f3e73f00e9d5906a455151d4bcee710e9bf7b80f519571b10340` |
| `results.md` | 2054 | `1c7441298ede21ab6b8325950a370bbfee0d08d12b4b533a5a08bc12d91f17b2` |
| `comparison.png` | 221604 | `a1f2a026924d37cf797c6e4bc13b09b7bb87952dd8bfa98c77bda7b0d3f1687e` |
| `effective_calibration_pools.json` | 38137 | `595898481855463edd091f76debca8742512e91923aa8f5b4ac1d521a72ab1a6` |
| `thresholds.json` | 1393 | `8a96f48dfbb4a1920fe93bba9012b860270b320ca0b28c4b8e20b23aa330b645` |
| `jobs.json` | 32448 | `7f6a564ea83dbb95c94c0f5ae2314ebed576ab23ffe648a0035e2913baf9e1dc` |
| `baseline.json` | 459757 | `b44ab2401d74fca81a90d2a83b507bfc91c97ccb81b5b72b4cef9b9268ba20df` |

## Preparation validation

The local `audit` action completed against the saved parent baseline and exact
response templates: 36 sites, 128 unique calibration centers, all six ±20
searches valid at every radius, and the narrower support counts shown above.
All selected valid center and neighbor PSD fits factored successfully. Python
syntax compilation and repository whitespace checks passed. No new reduction
was run, and no production or mxlib-calling function changed.
