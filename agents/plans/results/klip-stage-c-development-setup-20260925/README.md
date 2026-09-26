# KLIP Stage-C immutable development preparer

## Purpose

Stage B selected two 11-pixel covariance arms: strict-radial
Hann/mixing-0.1 as the leading data-selected model and raw
rectangular/mixing-0.3 as the mandatory P4 prior. This preparer freezes the
development and validation experiment before any positive reduction is made.

The maintained runner is
[`prepare_klip_stage_c_development.py`](../../scripts/prepare_klip_stage_c_development.py).
It verifies the completed 47-pixel exact-response campaign and its Stage-A
lineage, including every sparse-response product.

## Frozen geometry

At each nominal radius 7.5, 10, 12, 16, 20, and 24 pixels, preparation starts
from native centers satisfying all of these conditions:

- complete 47-pixel exact responses in all eight KL modes at the center and
  four axial search pixels;
- finite 11-pixel final-image data in all eight modes at those five pixels;
- strict leave-site-out 3.6-pixel radial-profile support;
- at least eight accepted 11-pixel Welch patches in each detector half for
  every search pixel; and
- search centers outside the fixed seven-pixel known-planet disk.

For each radius, the preparer chooses the smallest predeclared training
half-width in 5, 10, 20, 40, and 60 pixels that leaves at least 38 common
candidates. It then assigns disjoint locations by deterministic angular
maximin selection, in this order:

| Role | Locations per radius | Total |
| :--- | ---: | ---: |
| Maximum-null calibration | 20 | 120 |
| Development injection | 6 | 36 |
| Validation injection | 6 | 36 |
| Held-out null | 6 | 36 |

This selection uses geometry and finite-support masks only. Baseline scores,
positive images, and recovery do not enter it.

## Source levels

The preparer computes site-specific physical contrasts for both development
and validation sites. At mode 200, it smooths the exact 47-pixel unit response
with Gaussian FWHM 3.6, takes the maximum source-only response over the fixed
five-pixel search, and divides by the interpolated one-pixel annular deviation
of the Gaussian-smoothed signal-free image. The annular profile excludes both
the known planet and the trial site with the production `R + 0.5` boundary.
The production small-sample factor then sets nominal source-only SNR 3, 5, and
7 contrasts.

This creates 108 development reduction commands and separately freezes 108
future validation commands. Validation positives remain unopened.

## Frozen method suite

The method manifest contains:

- native, Gaussian FWHM 3.6, exact identity, sparse identity, exact-response
  LPF 1.8 and 2.7, and exact fitted-mean isotropic permanent references;
- the development-only Gaussian widths 2.4, 3.0, and 4.2;
- raw 11-pixel rectangular/mixing-0.3;
- strict-radial 11-pixel Hann/mixing-0.1; and
- clipped and hard-truncated Hann precision at 0.5, 0.75, and 1.0 times mean
  covariance variance, plus the full inverse.

Positive covariance analysis will compare per-image refits with
baseline-frozen weights. Every method receives a separate maximum-of-20 null
threshold at each radius. Calibration must be receipted before a positive
image is analyzed.

## ROC preparation

```bash
git pull

root=working/roc/klip_stage_c_development_20260925

taskset -c 12-27 python3 \
  agents/plans/scripts/prepare_klip_stage_c_development.py check

taskset -c 12-27 python3 \
  agents/plans/scripts/prepare_klip_stage_c_development.py prepare "$root"
```

Preparation writes `protocol.json`, `geometry.json`, `contrasts.json`, all 108
development commands, copied software/configuration, and a strict manifest.
The verified receipt below enabled the separate calibration-first development
runner.

The first two verification attempts exposed packaging-only imports before any
reduction directory was created. The preparer initially omitted the KLIP PSD
extension module, and that module eagerly imported a P4 comparison chain used
only by its standalone self-check. The corrected package copies the extension
module and defers the P4 import to that self-check. Both untouched preparation
directories are archived and replaced rather than repaired in place.

## Verified ROC preparation

The final materialization from source commit 6b63d1a passed the frozen loader
and all hash, lineage, resource-affinity, and parent-completion checks. It
contains 18 methods, 228 disjoint sites, 108 development commands, and 108
unopened validation tasks. The state remains prepared with zero development
reductions; positive analysis has not started and no validation product has
been opened.

| Radius (px) | Selected half-width (px) | Eligible candidates | Minimum profile pixels | SNR 3 contrast min--mean--max | SNR 5 contrast min--mean--max | SNR 7 contrast min--mean--max |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 7.5 | 20 | 42 | 21 | 0.005266--0.008031--0.011586 | 0.008776--0.013385--0.019309 | 0.012286--0.018739--0.027033 |
| 10 | 20 | 44 | 38 | 0.002759--0.003817--0.004716 | 0.004598--0.006362--0.007859 | 0.006437--0.008906--0.011003 |
| 12 | 20 | 60 | 44 | 0.002183--0.003122--0.003893 | 0.003639--0.005204--0.006489 | 0.005095--0.007286--0.009084 |
| 16 | 10 | 85 | 44 | 0.001313--0.001467--0.001583 | 0.002188--0.002446--0.002639 | 0.003063--0.003424--0.003695 |
| 20 | 5 | 92 | 44 | 0.000950--0.001014--0.001094 | 0.001583--0.001691--0.001824 | 0.002216--0.002367--0.002554 |
| 24 | 5 | 160 | 44 | 0.000866--0.000949--0.001017 | 0.001444--0.001581--0.001695 | 0.002022--0.002214--0.002373 |

The inner three radii require a 20-pixel training half-width. Radius 16 narrows
to 10 pixels, and radii 20 and 24 support five pixels. This is the intended
small-separation behavior: the training band expands only where the fixed
38-site partition requires it.

| Prepared file | Bytes | SHA-256 |
| :--- | ---: | :--- |
| protocol.json | 156560 | ffbec258b9e1109217363f6cc4a6f9889f1cd09fc6fe1fa1e9ea94d3e3908c86 |
| geometry.json | 111804 | b86f53a2eddd41a440046636ef819d7ea7d732a9cd3e44637808c84e75fad50e |
| contrasts.json | 43372 | a8013cccef8f281c29a6510071175ac70424d507b9f3990356c7023ae755d58a |
| commands.json | 139733 | 403a05503e08e7bc557cc6b8b6e01a5c21a2ac4bdd8f24219f9bc37abc545655 |
| manifest.json | 23627 | 833679e372c7b5d1bc28a6825fbbe7baca89d7e297a604e709a17846c60a6893 |
| state.json | 135 | 950fc3a06cc635f78baf893ec3e5819e05f49067e6f111366f674914257f922b |


## Calibration-first development runner

The maintained runner is
[run_klip_stage_c_development.py](../../scripts/run_klip_stage_c_development.py).
Its run action enforces this order:

1. fit and receipt all baseline maps, development-site frozen weights, and the
   20-site maximum-null threshold for every radius, mode, and method;
2. run or verify all 108 development reductions; and
3. analyze the positive images with paired baseline-frozen and per-image-refit
   covariance at the five search pixels.

The calibration builds only the one-pixel annuli needed to bracket the frozen
searches. It passes each arbitrary amplitude cube through hciAnalyze for the
production mean-subtracted, small-sample-corrected annular SNR and requires an
independent oracle to agree. The radial Hann arm fits in strict leave-site-out
standardized coordinates. Its full, clipped, and hard-truncated weights are
all renormalized to unit physical response. The raw rectangular arm and every
permanent reference are retained alongside it.

The exact-response amplitude fields start at radius 6 pixels. Seven of the 38
radius-7.5 sites have one inward search pixel below radius 6.5, where production
interpolation requests the empty 5--6-pixel annulus and hciAnalyze converts the
undefined result to zero. They are `r7p5_cal02`, `r7p5_cal03`, `r7p5_cal08`,
`r7p5_cal15`, `r7p5_cal16`, `r7p5_dev03`, and `r7p5_hold02`; all six frozen
validation sites have complete bracketing. At only this inner boundary, the
runner requires the production value to be zero and substitutes the nearest
supported one-pixel annular mean and sample deviation while evaluating the
small-sample factor at the candidate pixel radius. It records every substitution
by mode and method in `annular_verification.json`. The policy does not
extrapolate a response template or an amplitude. Everywhere else, production
and oracle SNR must still agree within the strict tolerance.

The runner copies itself into the already prepared software directory on its
first invocation and freezes that copy plus hciAnalyze in a separate
development manifest. It never reads or runs a validation command. Interrupted
calibration units, reductions, and analysis tasks are archived; completed
receipts are recursively verified and reused. A pre-calibration runner
correction is permitted only while the state is still prepared, calibration is
incomplete, and no positive reduction directory exists. The previous runner
and development manifest are then fingerprinted and archived before the
corrected runner is frozen. Completed calibration units are retained, while
site records predating the boundary-audit schema are regenerated. An
analysis-only repair is separately permitted after calibration and all 108
reductions are complete, provided validation remains unopened and no final
development receipt exists. It preserves the base runner, development manifest,
calibration receipt, and reductions; archives all partial analysis products;
and freezes the corrected runner in a separate analysis manifest.

Run the complete development stage on ROC with:

```bash
git pull

root=working/roc/klip_stage_c_development_20260925

taskset -c 12-27 python3 \
  agents/plans/scripts/run_klip_stage_c_development.py check

taskset -c 12-27 python3 \
  agents/plans/scripts/run_klip_stage_c_development.py run "$root" --workers 4 \
  2>&1 | tee "$root/driver.log"
```

The phases can instead be resumed separately with calibrate, reduce, and
analyze. The four analysis workers share the frozen 16-core affinity while
BLAS and hciAnalyze remain single threaded. The KLIP reductions run one at a
time with all 16 OpenMP threads.

## ROC runner pilot

The committed runner passed two isolated pilots in /tmp without modifying the
canonical prepared directory. First, all eight radius-7.5 calibration units
completed against the real baseline, exact field, sparse field, and strict
radial profiles. A real fixed null site then passed production hciAnalyze
normalization with a maximum independent-oracle difference of 2.38e-7.

Second, a synthetic positive made by adding the exact 11-pixel response at one
frozen development site traversed the complete analysis path. At mode 200 its
unweighted response-fidelity cosine was 0.999999999999998, its projection was
0.9999999924, and its best-scaled relative residual was 6.03e-8. Frozen and
per-image-refit covariance both completed, and the largest hciAnalyze/oracle
SNR difference across the two arms was 9.54e-7. These pilots validate the FITS
orientation, 47-to-11 support crop, site-specific training exclusions,
precision replay, response-fidelity calculation, and production SNR bridge.

The first canonical launch then completed all 48 calibration units and two
null sites before `r7p5_cal02` exposed the unsupported inner annular bin. It
stopped before thresholds, reductions, or positive analysis. Isolated reruns
of all five affected calibration sites and the affected development site each
completed all 8 modes and 13 methods. Each recorded the expected 104
substitutions at its single affected search pixel; the largest
production/oracle difference on supported pixels was 4.77e-7. The resumable
correction therefore preserves the expensive 48 unit products and regenerates
the inexpensive site records under one audited schema.

The resumed canonical run then completed all 120 calibration sites, froze the
thresholds, and completed all 108 KLIP reductions. Parallel positive analysis
wrote 81 complete and 27 incomplete task directories before reporting a
response-fidelity diagnostic error. Its unweighted branch constructed an
identity covariance after masking the 121-pixel stamp, then applied the
original mask a second time; real response supports as small as 93 pixels
therefore indexed beyond that reduced matrix. Filter amplitudes, SNRs,
thresholds, and reductions were not involved in the failure.

The correction keeps the identity covariance in the full 121-pixel space and
applies the response-support mask once, matching the two PSD-metric branches.
A deterministic masked-support check now covers this case. An isolated replay
of the previously failing `r7p5_dev02_snr3` real reduction completed all eight
modes and both covariance arms. At mode 200 its unweighted fidelity cosine was
0.9999021, its projection was 0.9944667, and its relative residual was
0.0139931; both SNR arms agreed with the annular oracle within 4.77e-7. The
analysis-only repair will archive all 108 partial analysis directories and
regenerate them while retaining calibration and reductions. Validation remains
unopened.
