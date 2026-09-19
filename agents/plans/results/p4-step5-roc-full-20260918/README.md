# Step 5: full PSD injection study on ROC

**The full study is prepared: one baseline and 90 positive full-image reductions, followed by automatic analysis.**
The user explicitly waived the separate small ROC numerical/throughput pilot. Native compilation, dependency
loading, input hashes, geometry, and local filter integration have been checked. This checkpoint fixes the
experiment before new-site scores or positive images are inspected; it contains no new recovery results.

## Fixed comparison

The [preceding PSD recovery experiment](../p4-step5-psd-recovery-20260918/README.md) motivates two fixed PSD
candidates and six controls. Every method shares the same reductions, sites, and calibration rule.

| Method key | Fixed estimator / statistic |
| --- | --- |
| `psd_hann_b0_m0.1` | Same-radius raw patches; Hann PSD window; 0.1 isotropic spectral mixture |
| `psd_rectangular_b5_m0.3` | Raw patches pooled over ±5 pixels; rectangular PSD window; 0.3 spectral mixture |
| `pca_b5_f1` | Raw ±5-pixel patches; three covariance modes; variance floor fraction 1 |
| `isotropic_b5` | Raw ±5-pixel patches; fitted mean and mean sample pixel variance |
| `identity` | Original response matched filter with `C=I`, without fitted mean or physical noise scale |
| `gaussian_snr` | Production `filter.lpfGaussFW=3.6`, followed by ordinary application annular SNR |
| `gaussian_raw` | The same production Gaussian-smoothed intensity |
| `identity_snr` | Response matched filter with `C=I`, followed by ordinary application annular SNR |

The PSD estimator is unchanged: 11×11 stamps, 21×21 FFT padding, finite lag covariance, positive spectral
mixture, and variance rescaling. Windows affect covariance estimation only; candidate pixels and templates
remain untapered. Radial/angular patch-center spacing is five pixels, with at least eight training patches.
No radial normalization or additional parameter sweep is included.

## Sites, source amplitudes, and holdouts

There are **30 new sites**, six at each of five nominal radii, with **0.5×, 0.75×, and 1×** brightnesses.
The levels concentrate on the transition between the earlier 0.5× and 1× measurements. The former 2× level,
which was already fully recovered, is omitted. These multipliers are not PSD SNRs.

| Radius (pixels) | New native centers `(x, y)` | Minimum same-radius patches across candidate searches |
| --- | --- | --- |
| 26 | (102, 128), (102, 124), (103, 119), (105, 114), (115, 104), (121, 102) | 15 |
| 32 | (121, 159), (96, 129), (96, 124), (97, 117), (115, 98), (121, 96) | 21 |
| 38 | (125, 166), (110, 161), (91, 140), (90, 127), (92, 114), (125, 89) | 28 |
| 44 | (127, 172), (108, 167), (92, 154), (84, 136), (85, 118), (127, 83) | 36 |
| 50 | (127, 178), (103, 171), (84, 153), (78, 128), (83, 105), (127, 77) | 25 |

Coordinates are zero-based native FITS-array pixels: `row` in the existing analysis records means `x`, and
`column` means `y`; the NumPy access is `[y, x]`. The field center is `(127.5, 127.5)`. Sites are selected from
integer centers within 0.65 pixel of each radius and azimuth 90–270 degrees, using geometry and finite support
only. Six evenly spaced angular ranks are chosen from the viable pool, with at least four pixels between
centers on a ring. The `b` suffix identifies angular rank, not an independent angular block.

Each new site's complete 15×15 Gaussian-kernel/five-pixel-search footprint is disjoint from the earlier
calibration/evaluation Gaussian measurement footprints, the three development measurement footprints, and
the known source circle. New sites can overlap one another. This reuses the same residual field, including
pixels that could previously have supplied training noise; it is not a blind new-field validation.

**Training and calibration use a separate holdout for each new site.** Masking all 30 new sites simultaneously
left only 0–2 same-radius patches in the initial geometry check. Keeping the former global historical holdout
also prevented adequate new inner-site sampling. The fixed replacement excludes:

- the known source circle;
- all 28 original eligible calibration searches, each with its complete Gaussian/search footprint;
- the current new site's complete Gaussian/search footprint.

Other new sites and historical non-calibration development/evaluation regions may supply training pixels.
The current target is excluded from every covariance fit used to calibrate or measure that target. All 28
calibration searches remain valid for every site. Means and covariances are fitted again on each positive
image with these same masks. Nonlocal changes caused by the full reduction can still alter training noise.

For each site and each method, the threshold is the `ceil((28+1)*0.95)` order statistic: the maximum of the
28 baseline calibration search scores. Searches maximize the signed statistic over the native center and
four axial one-pixel neighbors; detection requires strict exceedance. These are site-specific recalibrations,
not reused thresholds from the earlier global-mask experiment.

The site's reference contrast is its `isotropic_b5` calibration threshold times its baseline candidate
conditional sigma. This uses training covariance and the response template; target-stamp pixel values are
not used to set source brightness. **All 30 sets of thresholds and all 90 injection contrasts are written and
fingerprinted before any new-site baseline score or positive image is measured.** Baseline-exceeding sites
remain in the positive sample. Invalid searches count as nondetections; missing structural support aborts
visibly rather than silently changing the preassigned sample.

## ROC execution

The isolated run directory on `roc` (hostname `exao5`) is:

```text
/home/jrmales/Source/mxApps/hciReduce/working/roc/p4_psd_full_20260918
```

ROC's existing checkout is left intact. Native executables were built from archived C++ source commit
`882cb15`, with Release configuration, CPU benchmarks, and experimental P4 precision enabled. The frozen
executables are `p4ReductionPrecisionBenchmark`, `hciAnalyze`, and `hciGaussianReference`; their shared
libraries are copied into the isolated software directory and fingerprinted. No production C++ changed.

All **621 original input-frame hashes match on ROC**. The run preserves P4-M32D64, mode fraction 0.15,
256×256 full images, mean combination, the original reduction configuration and analytic response field,
and the cropped 12×12 injection PSF without renormalization. The baseline belongs to this full study; no
separate numerical-consistency or timing pilot is inserted.

Reductions use ROC's 24 physical cores (CPU IDs 0–23), OpenMP 24, and single-threaded BLAS. Analysis uses
CPU IDs 12–13, OpenMP 2, and single-threaded BLAS. Python is the existing `xpy3_13` environment:
NumPy 2.4.3, SciPy 1.16.3, Astropy 7.2.0, and Matplotlib 3.10.8. Native dependency, build, and Python
provenance is recorded in `manifest.json`; there were 740 GB available before launch.

The driver is intended to run in detached tmux session `p4-psd-full-20260918`. It performs the baseline,
calibration, 30 new-site baseline measurements, all 90 positive reductions and reference/filter analyses,
then writes `results.json`, `results.md`, `comparison.png`, and `complete.json`. Atomic `state.json`,
`driver.log`, per-trial completion records, commands, FITS products, and resource logs preserve progress.
It needs no connection or further input to finish. A lock prevents concurrent supervisors.

Completed reductions and reference maps can be reused after hash verification. Partial output directories
are deliberately not overwritten: preserve/quarantine a failed trial directory before restarting. A failure
sets `state.json` to `failed` and stops; completion sets it to `complete`. A second invocation after completion
verifies the summary products and returns.

Check progress without attaching to the supervisor:

```sh
ssh roc 'cat /home/jrmales/Source/mxApps/hciReduce/working/roc/p4_psd_full_20260918/state.json'
```

## Interpretation and verification

Recovery and baseline exceedances use the same fixed threshold procedure for all methods. Report raw
exact-center photometry separately: `positive amplitude / injected contrast − 1`, conditional-sigma coverage,
and the diagnostic paired increment relative to the baseline. Gaussian intensity is not a contrast estimate;
identity's `C=I` sigma is not a fitted physical uncertainty. Application-SNR references retain their ordinary
whole-annulus normalization, including trial neighborhoods.

The 30 positions share a residual field, calibration measurements, training pixels, and sometimes overlapping
measurement footprints. Counts are descriptive; they do not support independent-trial confidence intervals
or demonstrate equal underlying false-positive rates. Inner radii below 26 pixels remain a follow-up.

Local setup checks passed:

- independently reconstructed all 30 native masks;
- audited 150 candidate search stencils and 4,200 calibration search stencils;
- replacing target pixels with NaNs leaves 360 candidate/calibration covariance fits exactly unchanged;
- at a previously evaluated site, all four covariance filters and identity recover an exact template-only
  amplitude increment, without inspecting outcomes at the new sites;
- selected coordinates match a separate geometry prototype; original PSF and staged source hashes match.

The native build, library loading/help checks, and all input hashes pass on ROC. No mxlib-calling function
was edited, so this checkpoint adds no mxlib ownership follow-up.

### Artifacts

- [`protocol.json`](protocol.json): fixed methods, sites, amplitudes rule, exclusions, and compute settings.
- [`geometry.json`](geometry.json): every site's candidate and calibration training-count audit.
- [`native_manifest.json`](native_manifest.json): ROC-native executables, dependencies, source archive,
  build cache, Python environment, and all 621 source-frame fingerprints.
- [`setup_checks.json`](setup_checks.json), [`setup_check.py`](setup_check.py), and
  [`independent_geometry.json`](independent_geometry.json): local verification and separate selection audit.
- [`run_p4_step5_roc_full.py`](../../scripts/run_p4_step5_roc_full.py): maintained design/prepare/run driver.

The ignored run directory contains the frozen inputs and software. The tracked report preserves the protocol
and setup evidence; completed scientific results will be reviewed and archived at the next checkpoint.
