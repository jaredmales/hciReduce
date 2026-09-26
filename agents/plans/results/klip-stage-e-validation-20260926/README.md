# KLIP Stage-E fresh validation result

## Status and provenance

Stage E completed on ROC under the immutable Stage-D policy. The run contains
48 baseline-only validation model units, 36 held-out null analyses, 108 fresh
validation reductions, and 108 validation analyses. The final state is
`validation_complete`; the known planet remains unopened.

The compact result set comes from
`working/roc/klip_stage_c_development_20260925` at hciReduce commit
`8f654f4e8e83940011192b8602922bac84d09a60`. A recursive audit on ROC verified
the policy receipt, all nested model, held-out, reduction, and analysis
products, and both guarded repair receipts before this copy was made. Absolute
ROC paths in the receipts are retained as provenance; the large FITS products
remain on ROC.

Two implementation errors interrupted execution. Neither changed the frozen
policy, thresholds, reductions, or filter scores:

1. The first entry point unpacked a two-value preparation result as three
   values and stopped before opening a Stage-E product. The guarded repair
   archived the original software receipt and changed only the driver.
2. The analysis initially imported the unrepaired Stage-C response-fidelity
   helper, which applied a support mask twice to an identity matrix. This
   affected only the finite-response diagnostic. The guarded analysis-only
   repair retained all models, held-out results, and reductions, archived the
   partial analysis tree, and regenerated all 108 analyses with the recorded
   full-space convention.

The repair receipts and archive manifest are in `repairs/`.

## Frozen experiment

Mode 200 is the preregistered primary endpoint. Each source level contains six
validation sites at each of radii 7.5, 10, 12, 16, 20, and 24 pixels, for 36
paired sites. A detection is the maximum SNR over the center and four cardinal
neighbors, compared with the unchanged method-, mode-, and radius-specific
maximum-of-20 calibration threshold. Covariance weights were estimated from
the signal-free baseline before held-out or validation products were opened.

The primary comparator is Gaussian smoothing with FWHM 3.6 pixels. Gaussian
FWHM 2.4 is the stronger smoothing control. `exact_identity` and
`sparse_identity` correlate the image with the exact or radially averaged
measured response without covariance weighting.

## Primary mode-200 result

| Target SNR | Method | Recovered / 36 | Mean maximum SNR | Mean throughput |
| ---: | :--- | ---: | ---: | ---: |
| 3 | Native | 23 | 2.9472 | 1.0053 |
| 3 | Gaussian 2.4 | 26 | 3.4244 | 1.0078 |
| 3 | Gaussian 3.6 | 23 | 3.1239 | 1.0089 |
| 3 | Exact identity | 28 | 3.7296 | 1.0064 |
| 3 | Sparse identity | **29** | 3.6716 | 1.0065 |
| 3 | Raw rectangular PSD | 26 | 3.6608 | 1.0057 |
| 3 | Radial Hann, truncation 0.75 | 27 | **3.7255** | 1.0059 |
| 5 | Native | 34 | 4.4075 | 1.0006 |
| 5 | Gaussian 2.4 | 35 | 5.4820 | 1.0025 |
| 5 | Gaussian 3.6 | 35 | 4.9758 | 1.0033 |
| 5 | Exact identity | 35 | **5.8870** | 1.0014 |
| 5 | Sparse identity | **36** | 5.7840 | 1.0016 |
| 5 | Raw rectangular PSD | 35 | 5.7607 | 1.0009 |
| 5 | Radial Hann, truncation 0.75 | **36** | 5.8651 | 1.0012 |
| 7 | Native | 36 | 5.8672 | 0.9918 |
| 7 | Gaussian 2.4 | 36 | 7.4345 | 0.9934 |
| 7 | Gaussian 3.6 | 36 | 6.7536 | 0.9937 |
| 7 | Exact identity | 36 | **7.9871** | 0.9932 |
| 7 | Sparse identity | 36 | 7.8734 | 0.9934 |
| 7 | Raw rectangular PSD | 36 | 7.8168 | 0.9931 |
| 7 | Radial Hann, truncation 0.75 | 36 | 7.9349 | 0.9934 |

## Preregistered covariance gate

Both covariance candidates pass every frozen acceptance condition against
Gaussian 3.6.

| Candidate | Accepted | Held-out, candidate / G3.6 | Recovery at SNR 3 / 5 / 7 | Total, candidate / G3.6 | Paired mean-SNR difference at 3, 95% interval |
| :--- | :---: | :---: | :---: | :---: | :---: |
| Raw rectangular PSD | Yes | 3 / 3 | 26 / 35 / 36 | 97 / 94 | +0.5369 [0.1429, 0.9134] |
| Radial Hann, truncation 0.75 | Yes | 3 / 3 | 27 / 36 / 36 | 99 / 94 | +0.6016 [0.2083, 0.9703] |

The raw candidate's paired mean-SNR differences at target SNR 5 and 7 are
+0.7849 [0.3614, 1.1802] and +1.0631 [0.6151, 1.5004]. The radial candidate's
corresponding differences are +0.8893 [0.4591, 1.2955] and +1.1812 [0.7010,
1.6449]. Common support and throughput conditions pass at every source level.

## Stronger Gaussian control

The Gaussian-2.4 comparisons were specified as a reported control rather than
the formal gate. The same radius-stratified paired bootstrap shows that radial
covariance also improves on this stronger Gaussian at every source level. Raw
rectangular has positive point estimates, but its intervals include zero and
its recovery total equals Gaussian 2.4.

| Candidate minus Gaussian 2.4 | Target SNR | Recovery, candidate / G2.4 | Mean SNR difference | 95% interval | Site wins / 36 |
| :--- | ---: | :---: | ---: | :---: | ---: |
| Raw rectangular PSD | 3 | 26 / 26 | +0.2364 | [-0.0593, 0.5177] | 21 |
| Raw rectangular PSD | 5 | 35 / 35 | +0.2786 | [-0.0626, 0.5967] | 20 |
| Raw rectangular PSD | 7 | 36 / 36 | +0.3823 | [-0.0104, 0.7485] | 20 |
| Radial Hann, truncation 0.75 | 3 | 27 / 26 | +0.3011 | [0.0080, 0.5776] | 21 |
| Radial Hann, truncation 0.75 | 5 | 36 / 35 | +0.3830 | [0.0292, 0.7099] | 20 |
| Radial Hann, truncation 0.75 | 7 | 36 / 36 | +0.5004 | [0.0781, 0.8919] | 20 |

## Covariance versus identity response filtering

Validation does not show an advantage from covariance weighting over applying
the measured response with identity covariance. Exact and sparse identity both
beat Gaussian 2.4 significantly at every level. Radial covariance is
statistically indistinguishable from either identity response filter, has two
fewer faint recoveries than sparse identity, and has three held-out exceedances
versus one for each identity method.

| Difference | Target SNR | Recovery, first / second | Mean SNR difference | 95% interval |
| :--- | ---: | :---: | ---: | :---: |
| Exact identity minus Gaussian 2.4 | 3 | 28 / 26 | +0.3052 | [0.0349, 0.5683] |
| Exact identity minus Gaussian 2.4 | 5 | 35 / 35 | +0.4049 | [0.0851, 0.7125] |
| Exact identity minus Gaussian 2.4 | 7 | 36 / 36 | +0.5526 | [0.1688, 0.9371] |
| Sparse identity minus Gaussian 2.4 | 3 | 29 / 26 | +0.2473 | [0.0291, 0.4478] |
| Sparse identity minus Gaussian 2.4 | 5 | 36 / 35 | +0.3019 | [0.0746, 0.5213] |
| Sparse identity minus Gaussian 2.4 | 7 | 36 / 36 | +0.4389 | [0.1941, 0.6774] |
| Radial covariance minus exact identity | 3 | 27 / 28 | -0.0041 | [-0.0968, 0.0904] |
| Radial covariance minus exact identity | 5 | 36 / 35 | -0.0219 | [-0.1234, 0.0832] |
| Radial covariance minus exact identity | 7 | 36 / 36 | -0.0522 | [-0.1929, 0.0960] |
| Radial covariance minus sparse identity | 3 | 27 / 29 | +0.0539 | [-0.0821, 0.1923] |
| Radial covariance minus sparse identity | 5 | 36 / 36 | +0.0811 | [-0.1058, 0.2654] |
| Radial covariance minus sparse identity | 7 | 36 / 36 | +0.0615 | [-0.1831, 0.2887] |

Sparse identity gives the best combined recovery: 101/108, compared with 99
for exact identity and radial covariance, 97 for raw rectangular and Gaussian
2.4, and 94 for Gaussian 3.6. At target SNR 3, sparse identity recovers 29/36
in every one of the eight KL modes; exact identity recovers 28--29 and radial
covariance 25--28. This all-mode consistency supports the identity result.

## Separation dependence at target SNR 3

Each cell gives recoveries out of six followed by mean maximum SNR in
parentheses.

| Radius (pixels) | Gaussian 2.4 | Gaussian 3.6 | Exact identity | Sparse identity | Raw rectangular | Radial covariance |
| ---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 7.5 | 5 (3.559) | 4 (3.028) | 6 (4.691) | 6 (4.565) | 6 (4.548) | 6 (4.724) |
| 10 | 4 (3.661) | 4 (3.579) | 6 (3.927) | 6 (3.594) | 6 (3.677) | 6 (3.733) |
| 12 | 4 (2.869) | 3 (2.661) | 4 (2.752) | 4 (2.715) | 4 (2.796) | 4 (2.764) |
| 16 | 5 (3.835) | 4 (3.555) | 6 (4.256) | 6 (4.346) | 5 (4.073) | 6 (4.212) |
| 20 | 5 (3.135) | 5 (2.894) | 4 (3.030) | 5 (3.022) | 3 (2.966) | 2 (2.989) |
| 24 | 3 (3.488) | 3 (3.027) | 2 (3.721) | 2 (3.787) | 2 (3.905) | 3 (3.931) |

The response methods' gain is strongest at the intended small separations.
Outer-radius recovery remains threshold-tail sensitive: higher mean SNR does
not always imply more threshold crossings because each method has its own
frozen maximum-null threshold.

## Held-out nulls and response fidelity

At mode 200, frozen-threshold held-out exceedances are 3 for each Gaussian and
covariance method, 1 for native, and 1 for each identity response filter. Under
the frozen calibration resamples, median exceedance counts are 3 for Gaussian
3.6 and both covariance methods and 2 for each identity filter; their 97.5th
percentiles are 6, 6, 6, 5, and 6 for Gaussian 3.6, raw covariance, radial
covariance, exact identity, and sparse identity, respectively.

The finite positive response agrees closely with the independent signal-free
response. Across all 108 mode-200 measurements, unweighted response cosine has
minimum 0.998292 and median 0.999850; projection scale has median 1.00691; and
best-scaled relative residual has median 0.01732 and maximum 0.05843. The
maximum production-versus-independent annular-SNR discrepancy is
`2.8611e-6`, within the fixed `rtol=2e-6`, `atol=2e-6` comparison.

## Interpretation and next step

The original KLIP result that response matched filtering was worse than
Gaussian smoothing does not reproduce with this validated response and
injection protocol. Both exact and sparse identity response filters beat both
Gaussian controls. Covariance mismatch therefore does not explain the old
deficit in this experiment, and covariance weighting does not improve the best
identity response result. Radial covariance remains a validated
Gaussian-beating method, but sparse identity is the strongest frozen method by
recovery and held-out behavior.

Stage F should apply every frozen reference and finalist to the original
planet-bearing cube using `working/analyze.conf`. That endpoint is descriptive:
it checks position, aperture-maximum SNR, amplitude, and support without
altering this validation conclusion or selecting a method from one planet.

## Preserved files

- `validation_results.json`, `.csv`, and `.md`: generated all-mode Stage-E
  summaries and the preregistered gate.
- `measurements.json`: compact records for all 108 validation analyses,
  including all seven methods, eight modes, annular-oracle errors, and response
  fidelity.
- `derived_comparisons.json`: stronger Gaussian, identity, separation, mode,
  held-out, and response-fidelity comparisons derived from the frozen records.
- `heldout_results.json`: all 36 held-out records and threshold sensitivity.
- `policy/`: immutable policy, threshold audit, selection sensitivity,
  commands, and fixed resamples.
- `experiment/`: frozen protocol, geometry, and contrasts.
- `*_complete.json`: copied completion receipts.
- `verification.json`: local hashes and the ROC recursive-verification record.
