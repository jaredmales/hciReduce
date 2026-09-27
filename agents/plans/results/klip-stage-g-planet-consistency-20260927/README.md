# KLIP Stage-G exact-contrast planet-consistency result

## Purpose

This archive preserves the completed direct comparison between the known
planet and exact-contrast injections at twelve score-blind radius-12 centers.
Every injection is measured with the same three-pixel aperture used for the
planet. Mode 200 and the nominal PSF at the optimized negative-companion phase
are the primary arm.

The five arms test integer, configured-analysis, and optimized subpixel phases,
plus flux-preserving additional Gaussian PSF broadening of 0.9 and 1.8 pixels
at the optimized phase. All filters, training rules, positions, and endpoints
were fixed before the final measurements were summarized.

## Completion and verification

The final ROC state contains 60 complete reductions and 60 complete analyses.
The reporting-aperture and all-mode site-support repairs are documented in the
[setup record](../klip-stage-g-planet-consistency-setup-20260926/README.md).
The second repair retained 35 complete reduction-analysis pairs and replaced
five geometrically unsupported centers without using recovered SNR.

Local verification checks:

- all three aggregate product hashes against `complete.json`;
- all 60 analysis-receipt hashes;
- each local per-task `results.json` and
  `annular_verification.json` against its receipt; and
- the final task names against the repaired protocol.

The largest disagreement between production `hciAnalyze` and the independent
annular oracle is `1.90735e-6`, within the frozen absolute and relative
tolerance of `2e-6`. The checks are saved in
[`verification.json`](verification.json).

## Primary aperture-maximum result

The table uses the nominal optimized-phase arm. Injection entries are means
and sample standard deviations over twelve sites. The standardized deviation
is relative to that injection distribution.

| Method | Injections | Planet | Standardized deviation |
| :--- | ---: | ---: | ---: |
| Gaussian 2.4 | 5.0228 +/- 1.2608 | 5.4186 | +0.31 |
| Gaussian 3.6 | 4.8351 +/- 1.5345 | 5.6496 | +0.53 |
| Exact identity | 4.8422 +/- 1.0442 | 4.2375 | -0.58 |
| Sparse identity | 4.8652 +/- 0.9730 | 4.2461 | -0.64 |
| Raw rectangular PSD | 4.5520 +/- 0.8887 | 3.6937 | -0.97 |
| Radial Hann, truncation 0.75 | 4.6116 +/- 0.8635 | 3.6461 | -1.12 |

Every planet value is within 1.12 injection sample standard deviations. There
is no absolute-SNR inconsistency at the predeclared aperture-maximum endpoint.

### Paired aperture-maximum differences

| Difference | Injections | Planet | Standardized deviation | Predictive p |
| :--- | ---: | ---: | ---: | ---: |
| Gaussian 3.6 - exact identity | -0.0071 +/- 1.3236 | +1.4121 | +1.07 | 0.325 |
| Gaussian 3.6 - Gaussian 2.4 | -0.1877 +/- 0.4341 | +0.2310 | +0.96 | 0.374 |
| Raw covariance - exact identity | -0.2902 +/- 0.2735 | -0.5438 | -0.93 | 0.392 |
| Radial covariance - exact identity | -0.2305 +/- 0.3021 | -0.5913 | -1.19 | 0.276 |

The planet lies inside the twelve-injection range for all four comparisons.
Its exchangeable lower-tail ranks are 4/13 for raw covariance minus identity
and 3/13 for radial covariance minus identity. The Gaussian differences have
upper-tail ranks of 3/13. The predictive p-values use the diagnostic
Student-t calculation frozen in the runner; the ranks make no normality
assumption.

The matched-aperture experiment therefore resolves the earlier apparent
covariance anomaly. The covariance penalties are somewhat larger for the
planet but statistically compatible with the injections.

This does not favor covariance weighting. Across the injections, raw and
radial covariance reduce aperture-maximum SNR relative to exact identity by
0.290 and 0.231 on average. The validated identity filters remain preferable.

## Localization-sensitive secondary result

At the single nearest pixel, the covariance discrepancy remains:

| Difference | Injections | Planet | Standardized deviation | Predictive p |
| :--- | ---: | ---: | ---: | ---: |
| Raw covariance - exact identity | -0.2143 +/- 0.2096 | -0.9521 | -3.52 | 0.0061 |
| Radial covariance - exact identity | -0.1992 +/- 0.2435 | -1.0040 | -3.31 | 0.0088 |
| Gaussian 3.6 - Gaussian 2.4 | -0.2631 +/- 0.2599 | +0.5854 | +3.26 | 0.0095 |

For all three comparisons, the planet lies beyond all twelve injection values;
the corresponding one-sided exchangeable rank is 1/13. This is a real
position-sensitive difference, but it is secondary because a one-pixel
displacement is acceptable for this analysis and the three-pixel maximum was
predeclared as primary.

## Phase and PSF controls

Relative to the nominal optimized-phase arm, mean aperture-maximum SNR changes
are:

| Control minus nominal optimized phase | Gaussian 3.6 | Exact identity | Raw covariance | Radial covariance |
| :--- | ---: | ---: | ---: | ---: |
| Integer phase | +0.3056 +/- 0.2684 | +0.3681 +/- 0.4633 | +0.3623 +/- 0.4222 | +0.3690 +/- 0.4684 |
| Configured analysis phase | +0.3369 +/- 0.3863 | +0.4435 +/- 0.5951 | +0.4353 +/- 0.4961 | +0.4514 +/- 0.5389 |
| Additional Gaussian FWHM 0.9 | -0.0897 +/- 0.0246 | -0.1210 +/- 0.0228 | -0.1255 +/- 0.0225 | -0.1199 +/- 0.0215 |
| Additional Gaussian FWHM 1.8 | -0.7925 +/- 0.1998 | -0.9933 +/- 0.1614 | -1.0400 +/- 0.1750 | -1.0403 +/- 0.1760 |

The tested broadenings consistently lower SNR and do not reproduce a selective
Gaussian-3.6 advantage like the planet's.

Response fidelity identifies a subpixel modeling limitation:

| Arm | Nominal-template cosine | Best-scaled residual |
| :--- | ---: | ---: |
| Nominal integer | 0.999853 +/- 0.000059 | 0.0168 +/- 0.0035 |
| Nominal analysis phase | 0.989886 +/- 0.004977 | 0.1374 +/- 0.0363 |
| Nominal optimized phase | 0.925365 +/- 0.034303 | 0.3683 +/- 0.0874 |
| Broadening 0.9, optimized phase | 0.925580 +/- 0.034181 | 0.3678 +/- 0.0871 |
| Broadening 1.8, optimized phase | 0.917054 +/- 0.033811 | 0.3903 +/- 0.0783 |

The integer-grid exact response is accurate for integer injections but does
not reproduce the finite response as well near the planet's almost half-pixel
column phase. The same limitation is present in the injections, and the
three-pixel aperture maximum largely calibrates it out. A subpixel-aware
response template is the direct next response-model test if improved
single-pixel localization or photometry is required.

## Interpretation

The Stage-G result supports three conclusions:

1. The real planet is statistically consistent with matched exact-contrast
   injections at the predeclared aperture-maximum endpoint.
2. Covariance weighting remains worse than identity response filtering for
   these injections.
3. The response model is sensitive to subpixel phase. That sensitivity affects
   fixed-pixel scores and response fidelity, but it does not create a
   statistically exceptional planet after maximizing over the accepted
   three-pixel aperture.

The twelve sites share one residual field, so diagnostic predictive
probabilities do not represent independent observing epochs. The smallest
one-sided exchangeable rank available here is 1/13.

## Archive contents

- `results.md`, `results.json`, and `results.csv`: aggregate ROC report.
- `analysis/*/results.json`: all 60 per-task measurements.
- `analysis/*/annular_verification.json`: production-oracle checks.
- `analysis/*/complete.json`: task completion receipts.
- `protocol.json`, `manifest.json`, `state.json`, and
  `complete.json`: final frozen provenance and completion state.
- `repairs/`: compact repair records.
- `verification.json`: local compact-archive verification and
  distribution-free ranks.
- `driver.log`: final ROC resume log.

