# Step 5 integration pilot — 2026-09-17

This is the first implementation checkpoint for [Step 5](../../Covariance-Aware-Matched-Filtering.md).
It verifies annular training and analysis products on the existing full AF Lep science image.
**It does not measure completeness, calibrated false-positive rates, or a covariance detection gain.**

## Fixed policy

The three methods use the same 621-frame, **sigma-mean-combined science image** and analytic 11-by-11,
**mean-combined response field**. This combination distinction was corrected from the original report on
2026-09-18 after checking both FITS headers; the pilot did not use mean-combined science.
Training stays at the candidate radius with a default 5-pixel arc step, whole-footprint source/candidate exclusions,
and complete finite bilinear samples. The known source has a 10-pixel exclusion radius. There is no extra candidate
guard beyond its response footprint. Covariance models require eight accepted patches; PCA retains at most three
modes, and the variance floor is 0.1 times the median training-pixel variance. These settings were fixed before
inspecting results. Identity uses the original zero-mean, unit-covariance filter.

## Coverage and runtime

| Model | Process time (seconds) | Valid candidates | Too few training patches | Invalid response support |
| --- | ---: | ---: | ---: | ---: |
| identity | 0.121 | 8,616 | 0 | 2,188 |
| diagonal | 1.622 | 7,679 | 937 | 2,188 |
| pca | 7.957 | 7,679 | 937 | 2,188 |

Times are single sequential process measurements, including input and diagnostic FITS I/O, using two OpenMP
workers, one BLAS thread, and CPU affinity 0/2. They are not repeated performance benchmarks. All original science
and response input hashes remained unchanged. No zero-variance outcomes occurred. Of 8,616 identity-valid positions,
7,679 remain eligible for covariance weighting; the additional training requirement excludes 937 positions.
The other 2,188 response locations already fail the original support policy.

| Radius (pixels) | Identity valid | Diagonal/PCA valid | Median accepted training patches |
| --- | ---: | ---: | ---: |
| 0–10 | 316 | 0 | — |
| 10–20 | 948 | 433 | 10 |
| 20–30 | 1,564 | 1,564 | 20 |
| 30–40 | 2,196 | 2,196 | 39 |
| 40–50 | 2,836 | 2,836 | 51 |
| 50–60 | 756 | 650 | 57 |

Training coverage is complete on existing valid response locations from 20–50 pixels. Inner radii have fewer
available patches after whole-neighborhood exclusions. Near the outer boundary, rotated interpolation footprints
lose support before some native response stamps do. Evaluation must use common eligible support and retain these
failures in the reported coverage. Patch counts include overlap and are not effective independent sample counts.

At the nearest pixel to the configured AF Lep location (first/second image indices 139, 126), only **six** training
patches survive. Identity is valid there, while diagonal and PCA return `tooFewSamples` under the fixed minimum of
eight. The pilot therefore supplies no covariance-filter measurement at that pixel. Inner-radius footprint,
guard, and regularization choices must be evaluated on held-out locations; a nearby valid aperture pixel must not
be presented as a covariance improvement at the known source.

## Products and validation

Each run exports amplitude, conditional sigma, signed score, support, absolute covariance floor, retained rank,
accepted/attempted/excluded/incomplete training counts, and explicit status, alongside the separate empirical
annular SNR map. Identity diagnostics are opt-in. Metadata marks conditional scores as uncalibrated.

The final build passed six regression suites: **107 test cases**, including 1,141 annular-training assertions and
445 application assertions. Tests cover geometry and vectorization, complete-footprint exclusions under source
mutation, missing data, dense covariance agreement, constant-image variance, identity preservation, CLI parsing,
FITS product values/metadata, and explicit insufficient-training status. A focused Doxygen build links the tests
with the production APIs; its included-fixture preprocessor warning is also reproduced with the pre-change test
file. The mxlib audit and concrete non-blocking upstream coverage follow-ups are in
[`mxlib_cleanup.md`](../../mxlib_cleanup.md).

## Archive and reproduction

Raw products, frozen analysis software and source files, input hashes, exact commands, regression logs, and the
mxlib coverage audit are in `working/roc/p4_noise_step5_20260917_final/`. The first pilot, before prioritizing
invalid response support over insufficient training in the status code, remains at
`working/roc/p4_noise_step5_20260917/`; its valid-candidate count is identical.

```sh
python3 agents/plans/scripts/run_p4_step5_pilot.py \
  --science working/roc/p4_analytic_step3_20260917/analytic_sparse_batch32/finim.fits \
  --manifest working/roc/p4_analytic_step3_20260917/analytic_sparse_batch32/finim_outputs/p4PSF_manifest.fits \
  --binary _build_fresh/src/hciAnalyze --library _build_fresh/src/libhcireduce.so \
  --lambda-d 2.5 --source-radius 10 \
  --output /tmp/p4-step5-pilot-replay
```

The output directory must be new. Next work is the predefined spatial split, threshold calibration, and fresh
raw positive-injection study described in the plan. The Step-3 local, baseline-subtracted injection products
remain response-bias diagnostics and do not substitute for those measurements.
