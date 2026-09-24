# KLIP Stage-B local mode-200 noise screen

This baseline-only screen uses exact response templates, raw 11-by-11 Welch training patches, candidate-specific five-footprint exclusions, and independent detector-half fits. Scores do not subtract a fitted candidate mean in the common 11/31/47 comparison.

## Controlling radii: 7.5, 10, and 12 pixels

| Response | Band | Method | Score variance | |Score mean| | Within ±1σ | Split weight cosine | Split score corr. | Median 5-pixel max | 11px opposite-half variance |
| ---: | :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 11 | narrow | identity | 6.164 | 0.400 | 0.233 | 1.000 | 1.000 | 1.719 | 3.968 |
| 11 | narrow | rectangular_m0.1 | 2.340 | 0.352 | 0.425 | 0.986 | 0.985 | 1.211 | 1.146 |
| 11 | narrow | rectangular_m0.3 | 2.633 | 0.363 | 0.383 | 0.989 | 0.989 | 1.310 | 1.335 |
| 11 | narrow | hann_m0.1 | 2.232 | 0.357 | 0.425 | 0.991 | 0.984 | 1.203 | 1.133 |
| 11 | narrow | hann_m0.3 | 2.508 | 0.366 | 0.392 | 0.993 | 0.988 | 1.291 | 1.306 |
| 11 | full | identity | 7.424 | 0.441 | 0.217 | 1.000 | 1.000 | 1.890 | 2.994 |
| 11 | full | rectangular_m0.1 | 2.981 | 0.397 | 0.333 | 0.994 | 0.996 | 1.534 | 1.082 |
| 11 | full | rectangular_m0.3 | 3.370 | 0.408 | 0.333 | 0.996 | 0.997 | 1.641 | 1.239 |
| 11 | full | hann_m0.1 | 2.821 | 0.382 | 0.350 | 0.998 | 0.998 | 1.425 | 1.049 |
| 11 | full | hann_m0.3 | 3.213 | 0.394 | 0.333 | 0.998 | 0.999 | 1.541 | 1.201 |
| 31 | narrow | identity | 6.827 | 0.473 | 0.258 | 1.000 | 1.000 | 1.977 | — |
| 31 | narrow | rectangular_m0.1 | 2.680 | 0.338 | 0.400 | 0.989 | 0.986 | 1.138 | — |
| 31 | narrow | rectangular_m0.3 | 3.044 | 0.369 | 0.367 | 0.991 | 0.990 | 1.186 | — |
| 31 | narrow | hann_m0.1 | 2.557 | 0.336 | 0.433 | 0.996 | 0.992 | 1.050 | — |
| 31 | narrow | hann_m0.3 | 2.936 | 0.363 | 0.433 | 0.997 | 0.994 | 1.108 | — |
| 31 | full | identity | 6.840 | 0.474 | 0.250 | 1.000 | 1.000 | 1.989 | — |
| 31 | full | rectangular_m0.1 | 2.645 | 0.338 | 0.408 | 0.990 | 0.989 | 1.066 | — |
| 31 | full | rectangular_m0.3 | 3.015 | 0.369 | 0.383 | 0.992 | 0.992 | 1.264 | — |
| 31 | full | hann_m0.1 | 2.544 | 0.336 | 0.450 | 0.996 | 0.995 | 1.027 | — |
| 31 | full | hann_m0.3 | 2.927 | 0.365 | 0.400 | 0.998 | 0.996 | 1.136 | — |
| 47 | narrow | identity | 6.002 | 0.467 | 0.250 | 1.000 | 1.000 | 1.883 | — |
| 47 | narrow | rectangular_m0.1 | 2.452 | 0.383 | 0.417 | 0.988 | 0.984 | 1.273 | — |
| 47 | narrow | rectangular_m0.3 | 2.805 | 0.398 | 0.417 | 0.991 | 0.988 | 1.378 | — |
| 47 | narrow | hann_m0.1 | 2.222 | 0.354 | 0.467 | 0.994 | 0.990 | 1.191 | — |
| 47 | narrow | hann_m0.3 | 2.580 | 0.372 | 0.433 | 0.996 | 0.992 | 1.193 | — |
| 47 | full | identity | 6.011 | 0.467 | 0.250 | 1.000 | 1.000 | 1.880 | — |
| 47 | full | rectangular_m0.1 | 2.399 | 0.383 | 0.417 | 0.991 | 0.987 | 1.273 | — |
| 47 | full | rectangular_m0.3 | 2.756 | 0.398 | 0.417 | 0.993 | 0.990 | 1.258 | — |
| 47 | full | hann_m0.1 | 2.207 | 0.354 | 0.467 | 0.996 | 0.992 | 1.125 | — |
| 47 | full | hann_m0.3 | 2.569 | 0.372 | 0.433 | 0.997 | 0.994 | 1.147 | — |

The score-variance columns use 120 directional scores per radius (60 fixed queries × two fits) and remain spatially correlated. The five-pixel maximum uses 24 directional searches per radius. The 11-pixel opposite-half column is the median, across query-specific fits, of the variance of genuinely disjoint held-out patches divided by the conditional prediction.

For 31- and 47-pixel responses, complete response-sized held-out training stamps do not exist under the frozen exclusion geometry. Their validation here is therefore limited to fixed candidate-null calibration, split-weight stability, and split-score agreement.

This first screen covers raw pixels and no candidate-mean subtraction. Radial-variance standardization, post-mean patch-RMS normalization, the direct 11-pixel PCA/diagonal grid, and the support-independent radial-mean control remain separate Stage-B arms.

## Result

The raw Welch model is a large improvement over isotropic weighting, but its
conditional candidate scores are not calibrated yet. Across the controlling
radii, narrow-band identity score variance has median values 6.16, 6.83, and
6.00 for the 11-, 31-, and 47-pixel responses. Hann/mixing-0.1 reduces those
values to 2.23, 2.56, and 2.22. Rectangular/mixing-0.3, retained as the P4 prior,
reaches 2.63, 3.04, and 2.80.

Hann/mixing-0.1 is the leading raw-pixel PSD family. Its split-weight cosine is
0.991--0.994 with the narrow bands and 0.996--0.998 with the full band over the
three response supports. The paired split scores correlate at 0.984--0.990 for
the narrow bands. The stronger model agreement from the full band does not
improve candidate calibration; it raises the 11-pixel primary score-variance
median from 2.23 to 2.82 and changes the larger-footprint result very little.

The 11-pixel opposite-half patch diagnostic isolates the discrepancy. For
Hann/mixing-0.1, its median variance divided by conditional prediction is 1.13
with the narrow bands and 1.05 with the full bands, while the corresponding
fixed candidate-null variances are 2.23 and 2.82. Thus the short-lag PSD predicts
generic disjoint patch projections fairly well, but the fixed candidate
locations retain additional spatial or radial heterogeneity. The largest
mismatch is at radius 12: narrow Hann/mixing-0.1 candidate variance is 3.84,
4.90, and 4.44 at response supports 11, 31, and 47. Center-only values remain
2.41--5.11 across the three inner radii, so the excess is not created by pooling
the four neighboring search pixels.

Fitted 11-pixel ensemble-mean subtraction does not resolve the mismatch. For
narrow Hann/mixing-0.1 it changes candidate variance by less than 0.05 at every
radius and slightly increases the median opposite-half mean-square error at the
three controlling radii. The 47-pixel response improves score variance at
radius 7.5 but worsens it at radius 12; the baseline alone does not support
selecting a response footprint.

## Decision

Do not promote a covariance candidate from this raw arm alone. Carry
Hann/mixing-0.1 as the leading data-selected PSD family and
rectangular/mixing-0.3 as the mandatory P4 prior while Stage B tests the
prespecified radial-variance standardization and patch-RMS control. The
radius-dependent candidate mismatch makes radial standardization the next
focused test. Direct 11-pixel PCA/diagonal covariance and the larger-support
radial-mean control remain required before freezing the injection shortlist.
No positive injection or known-planet value was inspected here.

## Reproduction and provenance

The local runner is
[`run_klip_stage_b_local_noise_screen.py`](../../scripts/run_klip_stage_b_local_noise_screen.py).
It consumed a 2,897,949-byte compact bundle with SHA-256
`ffc66934b7b8f535ac7e6bc45ec3430abb95eaad844e7b9fb5938d22688dc7b4`,
exported by
[`export_klip_stage_b_local_bundle.py`](../../scripts/export_klip_stage_b_local_bundle.py).
The bundle contains the mode-200 baseline, 348 exact 47-pixel response stamps,
validity, and the complete fixed preflight geometry. It retains the source
response and preflight receipt hashes without copying the 1.6-GB response
campaign.

The run completed 10,800 policy/query records and replayed all 1,080 query
geometries exactly in 285.8 seconds. The numerical extension check and every
saved-product receipt pass. Compact machine-readable values are in
[`primary.json`](primary.json), run receipts are in
[`verification.json`](verification.json), and the comparison figure is
[`comparison.png`](comparison.png). The complete 16.8-MB record table remains
under ignored `working/local` data.


## Follow-up: radial normalization

The subsequent
[strict radial-normalization screen](../klip-stage-b-radial-normalization-20260923/README.md)
reduces 11-pixel Hann/mixing-0.1 candidate score variance from 2.23 to 1.86 but
leaves radius 12 at 3.10. Its leave-site-out profile is unsupported on common
31- or 47-pixel response geometry. Radial scaling is therefore retained as the
leading 11-pixel coordinate choice, without promoting a final filter.
