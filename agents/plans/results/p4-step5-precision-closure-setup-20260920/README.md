# Step 5: P4 precision-closure setup

## Questions fixed before the comparison

This closure study addresses three remaining questions from the completed
rectangular-PSD precision experiment:

1. Do the inner-radius gains persist when the covariance, fitted mean, and
   matched-filter weight are learned only from the uninjected baseline?
2. How sensitive are the maximum-null thresholds and recovery counts to the
   available common-support calibration locations?
3. What SNR do the selected policies measure for the known planet under
   `working/analyze.conf`?

The study reuses the completed 36 sites, 108 saved positive images, response
field, source masks, one-lambda/D trial exclusion, and production annular SNR.
It performs **no new P4 reduction**.

## Baseline-frozen covariance replay

Four raw rectangular ±20 PSD precision policies are fixed:

- the full inverse;
- eigenvalue clipping at 1.0 times the mean covariance variance;
- hard truncation below 0.5 times the mean covariance variance;
- hard truncation below 0.75 times the mean covariance variance as the
  aggressive diagnostic.

For each injection site and every candidate pixel needed by its five-pixel
search and annular normalization, the runner fits the PSD covariance and mean
once on the **uninjected baseline** using the same trial-specific training
exclusion as the completed comparison. It forms each unit-response weight
once, then applies that unchanged weight and baseline mean to the three saved
positive images at that site.

The paired control copies the completed maps that refit covariance separately
on each positive image. The frozen and refit members of each policy have
identical baseline maps by construction and therefore share the completed
maximum-null threshold. The runner independently recomputes every frozen
baseline amplitude and requires agreement with the completed parent map.

Identity with the 1.8-pixel response low pass and production Gaussian FWHM 3.6
remain copied references. They are not changed by covariance freezing.

This design isolates positive-image covariance adaptation. Any frozen-versus-
refit difference can come from changes in the covariance estimate or fitted
mean caused by the positive reduction. The injected source and P4 reduction
itself remain unchanged.

## Common-support threshold audit

For every nominal radius, all calibration trials inside the completed parent's
frozen radial band that are valid for all four precision policies form a paired
common-support set. Their expected counts are 20, 36, 24, 20, 20, and 20 at
radii 6, 8, 12, 16, 20, and 24 pixels.

The audit records each policy's median, upper-tail quantiles, second-largest
score, maximum, and gap between the two largest scores. It then makes 10,000
deterministic paired draws of 20 locations without replacement and reports the
resulting threshold, recovery-count, and held-out-null distributions for both
the refit and frozen positive maps. The random seed is `20260920`.

Only radii 8 and 12 have more than 20 common candidates, so only those radii
can vary under this fixed-band subset audit. Radius 6 and radii 16–24 explicitly
retain a single possible 20-location set. The audit diagnoses calibration-pool
sensitivity; it does not create independent null samples or justify a formal
false-alarm probability.

## Known-planet endpoint

The known planet is measured with the four precision policies, smoothed
identity, and Gaussian FWHM 3.6. The runner follows `working/analyze.conf`:

- lambda/D: 3.6 pixels;
- planet exclusion radius: 7 pixels plus the production half-pixel boundary;
- SNR minimum radius: 6 pixels;
- search aperture radius: 3 pixels;
- production annular mean subtraction and small-sample correction.

The precision maps cover every complete annulus required by all 39 native
pixels in the planet aperture. Covariance training excludes the known-source
circle and the 15-by-15 footprints of all aperture pixels. This is deliberately
conservative and matches the earlier planet diagnostic's treatment of source
contamination. The production `hciAnalyze` maps must agree with an independent
annular oracle.

Planet SNR is a descriptive endpoint only. It is not compared with the
injection thresholds and cannot select a method from one source.

## Validation completed before launch

The setup passed:

- local precision algebra and deterministic subset-sampling checks;
- Python compilation and strict-JSON summary generation;
- immutable preparation against the completed ROC parent, verifying 3,597
  frozen inputs;
- an end-to-end radius-6 site replay covering all three source levels;
- planet geometry with 604 fitted map pixels, 39/39 supported aperture pixels,
  and at least 55 covariance-training patches;
- an end-to-end planet `hciAnalyze` run and independent annular-oracle check.

## ROC commands

After pulling the setup commit, run the algebra check:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  MPLCONFIGDIR=/tmp/p4-step5-matplotlib \
  /opt/conda/envs/xpy3_13/bin/python3 \
  agents/plans/scripts/run_p4_step5_precision_closure.py check
```

Prepare the immutable closure study:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  MPLCONFIGDIR=/tmp/p4-step5-matplotlib \
  /opt/conda/envs/xpy3_13/bin/python3 \
  agents/plans/scripts/run_p4_step5_precision_closure.py prepare \
  --parent working/roc/p4_psd_precision_20260920 \
  --config working/analyze.conf \
  --root working/roc/p4_precision_closure_20260920 \
  --cpus 0 1 2 3 4 5 6 7 8 9 10 11
```

Launch the resumable replay:

```sh
tmux new-session -d -s p4-precision-closure \
  "cd /home/jrmales/Source/mxApps/hciReduce && \
   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
   MPLCONFIGDIR=/tmp/p4-step5-matplotlib \
   /opt/conda/envs/xpy3_13/bin/python3 \
   agents/plans/scripts/run_p4_step5_precision_closure.py run \
   --root working/roc/p4_precision_closure_20260920 \
   > working/roc/p4_precision_closure_20260920/driver.log 2>&1"
```

Monitor it with:

```sh
cat working/roc/p4_precision_closure_20260920/state.json
tail -n 30 working/roc/p4_precision_closure_20260920/driver.log
```

Final products will include `results.json`, `results.md`, and `comparison.png`,
plus the planet amplitude/SNR maps and their independent-oracle diagnostics.
