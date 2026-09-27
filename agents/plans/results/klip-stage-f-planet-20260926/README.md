# KLIP Stage-F known-planet closure result

## Status and provenance

Stage F completed on ROC after the Stage-E policy and validation result were
immutable. The completion receipt verifies 11 direct products, the state is
`closure_complete`, and the result records `known_planet_opened=true`,
`stage_e_policy_changed=false`, and `method_selection_performed=false`.

The result comes from
`working/roc/klip_stage_c_development_20260925/stage_f_planet` using the frozen
runner prepared from hciReduce commit `50fb08e`. The original planet-bearing
cube matched its Stage-A inventory before it was opened. The analysis uses the
source geometry in `working/analyze.conf`: separation 11.782 pixels, PA 262.051
degrees, a seven-pixel noise exclusion, and a three-pixel aperture.

A source-side audit verified every completion-product size and SHA-256 digest.
The large amplitude, response, policy-radius, and SNR FITS products remain on
ROC; this directory preserves the compact tables, diagnostics, commands,
receipts, logs, and their hashes.

## Mode-200 closure

| Method | Validation SNR-3 recovery | Held-out exceedances | Planet nearest-pixel SNR | Planet aperture-max SNR | Peak offset (pixels) | Peak contrast | Difference from optimized fit |
| :--- | :---: | :---: | ---: | ---: | ---: | ---: | ---: |
| Native | 23 / 36 | 1 | 3.4662 | 5.1536 | 0.841 | 0.005656 | +23.7% |
| Gaussian 2.4 | 26 / 36 | 3 | 5.0641 | 5.4186 | 0.841 | 0.005324 | +16.4% |
| Gaussian 3.6 | 23 / 36 | 3 | 5.6496 | **5.6496** | 0.213 | 0.005723 | +25.1% |
| Exact identity | 28 / 36 | 1 | 4.1813 | 4.2375 | 1.204 | 0.004283 | -6.4% |
| Sparse identity | **29 / 36** | **1** | 4.2461 | 4.2461 | 0.213 | 0.004171 | -8.8% |
| Exact response LPF 1.8 | Not in Stage E | Not in Stage E | 4.5358 | 4.5358 | 0.213 | 0.004343 | -5.1% |
| Exact response LPF 2.7 | Not in Stage E | Not in Stage E | 4.8662 | 4.8662 | 0.213 | 0.004668 | +2.0% |
| Exact identity, fitted mean | Not in Stage E | Not in Stage E | 4.1930 | 4.2794 | 1.204 | 0.004330 | -5.3% |
| Raw rectangular PSD | 26 / 36 | 3 | 3.2291 | 3.6937 | 1.204 | 0.003803 | -16.9% |
| Radial Hann, truncation 0.75 | 27 / 36 | 3 | 3.1773 | 3.6461 | 1.204 | 0.003781 | -17.3% |

The contrast reference is the independently optimized negative-companion fit,
0.00457436. Gaussian and native amplitudes are divided by their response to the
exact template for this descriptive comparison. Their broad kernels maximize
detection SNR but do not provide the closest contrast. Exact-response LPF 2.7
is closest to the optimized fit at mode 200.

## Behavior across KL modes

| Method | Mean aperture-max SNR | Minimum | Maximum |
| :--- | ---: | ---: | ---: |
| Native | 5.1010 | 4.9304 | 5.2919 |
| Gaussian 2.4 | 5.4666 | 5.3423 | 5.6071 |
| Gaussian 3.6 | **5.7076** | 5.5507 | 5.8313 |
| Exact identity | 4.3549 | 4.2375 | 4.5218 |
| Sparse identity | 4.3414 | 4.2395 | 4.4611 |
| Exact response LPF 1.8 | 4.6537 | 4.5358 | 4.7737 |
| Exact response LPF 2.7 | 4.9739 | 4.8662 | 5.0656 |
| Exact identity, fitted mean | 4.3915 | 4.2794 | 4.5592 |
| Raw rectangular PSD | 3.7976 | 3.6937 | 3.8713 |
| Radial Hann, truncation 0.75 | 3.7337 | 3.6461 | 3.8117 |

Gaussian 3.6 has the highest planet SNR in all eight modes. Smoothing the exact
response improves its SNR monotonically from identity through LPF 1.8 and LPF
2.7, but LPF 2.7 remains below Gaussian 3.6 in every mode. Sparse and exact
identity are nearly equal, and fitted-mean subtraction changes the result only
slightly. Neither response interpolation nor the fitted mean explains the
planet deficit.

Both covariance candidates reduce planet SNR in all eight modes. At mode 200,
raw rectangular loses 0.544 SNR relative to exact identity and radial Hann
loses 0.591. Their peak contrast estimates are about 17% below the optimized
fit, compared with 6--9% for the identity response filters. The peak covariance
fit has complete 121-pixel response support, 100 Welch samples, a minimum of 42
pixels in the radial profile, and 47 retained radial precision modes, so this
result is not caused by rejected aperture support or a failed fit.

## Closure interpretation

The fresh injections and the real planet answer different parts of the
question:

- On the immutable validation injections, exact and sparse identity response
  filters beat both Gaussian controls, and covariance does not improve on
  identity.
- On the one real planet, Gaussian 3.6 is consistently best, response smoothing
  progressively helps, and covariance makes the response filter worse.
- Exact-response filters give contrast estimates closer to the independent
  negative-companion fit than the highest-SNR Gaussian filter.

Covariance mismatch therefore does not explain the original real-planet
matched-filter deficit. At this location, the measured response's fine and
negative-lobe structure reduces detection SNR; progressively suppressing that
structure helps. The tested covariance models amplify the disagreement between
the real planet and the exact response rather than correcting it.

The injection-versus-planet difference remains scientifically important. The
finite injected response agrees very closely with the independently measured
template, so the completed experiment validates the estimator for sources made
with the injected PSF. It does not establish that those injections reproduce
every relevant property of the real companion and local residual. A second
real companion, epoch, or target is the appropriate next generalization test.
The single AF Lep planet cannot override the injection ranking.

## Verification

- All 39 native pixels in the configured aperture have common support.
- The production and independent annular SNR calculations agree to a maximum
  absolute error of `9.536743e-7`.
- The planet peak moves by at most about 1.2 pixels among the reported methods;
  this is within the previously accepted localization tolerance.
- The generic annular maps retain the nearest frozen Stage-C radial policy.
  Aperture weights were fitted from the signal-free baseline with the complete
  aperture excluded from training.

## Preserved files

- `planet/results.json`, `.csv`, and `.md`: generated planet and closure tables.
- `derived_comparisons.json`: reproducible all-mode SNR, contrast, ranking, and
  pairwise method comparisons.
- `planet/diagnostics.json`: aperture support, training, radial-profile, and
  retained-mode diagnostics.
- `planet/annular_verification.json`: per-method, per-mode oracle comparison.
- `manifest.json`, `state.json`, and `complete.json`: frozen inputs, final state,
  and source-product receipt.
- `analyze.conf`, commands, and logs: exact execution settings and records.
- `verification.json`: local compact-file hashes and source verification
  summary.
