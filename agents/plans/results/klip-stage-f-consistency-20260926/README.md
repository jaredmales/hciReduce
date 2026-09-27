# KLIP Stage-F planet/injection consistency analysis

## Question and comparison

This checkpoint asks whether the mode-200 AF Lep planet is statistically
consistent with the fresh Stage-E injections. It uses the six held-out
validation positions at nominal radius 12 pixels and the independently
optimized planet contrast, `0.004574362496414845`. At each position the
reported value is linearly interpolated between the two bracketing SNR-3/5/7
contrasts. The interpolated source-only calibration spans SNR 3.707--6.285
across the six position angles.

The center comparison is exact: injection center SNR is compared with planet
nearest-pixel SNR. The search comparison is descriptive because Stage E uses
the maximum over the center and four cardinal neighbors, while Stage F uses
all 39 valid pixels in the radius-three aperture.

## Absolute SNR

| Method | Injection center SNR, mean ± sample std | Planet nearest-pixel SNR | Injection max-of-five SNR, mean ± sample std | Planet aperture maximum |
| :--- | ---: | ---: | ---: | ---: |
| Native | 3.783 ± 0.670 | 3.466 | 4.091 ± 0.781 | 5.154 |
| Gaussian 2.4 | 4.325 ± 0.931 | 5.064 | 4.387 ± 0.902 | 5.419 |
| Gaussian 3.6 | 4.008 ± 1.065 | 5.650 | 4.116 ± 1.058 | 5.650 |
| Exact identity | 3.902 ± 1.073 | 4.181 | 4.072 ± 1.010 | 4.237 |
| Sparse identity | 3.997 ± 1.004 | 4.246 | 4.103 ± 0.892 | 4.246 |
| Raw rectangular covariance | 3.954 ± 1.051 | 3.229 | 4.088 ± 0.961 | 3.694 |
| Radial-Hann covariance | 3.946 ± 1.064 | 3.177 | 4.078 ± 0.947 | 3.646 |

No individual centered planet SNR is a gross outlier: every method lies within
1.55 sample standard deviations of its comparable injection distribution.
The absolute scatter is dominated by position-dependent residual structure and
is strongly correlated among methods.

## Paired filter behavior

Pairing methods at each injection position removes much of that common spatial
scatter.

| Center-SNR difference | Injections, mean ± sample std | Planet | Standardized difference |
| :--- | ---: | ---: | ---: |
| Gaussian 3.6 − exact identity | +0.106 ± 1.190 | +1.468 | +1.15 |
| Gaussian 3.6 − sparse identity | +0.011 ± 0.788 | +1.404 | +1.77 |
| Gaussian 3.6 − Gaussian 2.4 | −0.317 ± 0.341 | +0.585 | +2.65 |
| Raw covariance − exact identity | +0.052 ± 0.180 | −0.952 | −5.57 |
| Radial covariance − exact identity | +0.044 ± 0.144 | −1.004 | −7.27 |

The planet's broad-Gaussian advantage is somewhat unusual but remains within
the observed absolute injection scatter. The covariance penalty relative to
the same exact-response identity filter is the specific inconsistency: the
planet value is below all six comparable injections. Under an independent,
exchangeable normal predictive model, the nominal two-sided probabilities are
0.0036 for raw covariance and 0.0011 for radial covariance. Those values are
diagnostic because spatial independence and normality are not established.
With six validation angles, the smallest one-sided exchangeable rank
probability is 1/7.

## Interpretation and next test

The existing data support a localized covariance-versus-identity anomaly, not
a population-level declaration that the planet is inconsistent with injected
sources. A direct closure test should inject the optimized contrast at
additional radius-12 position angles, use the same radius-three aperture as
the planet analysis, and preregister covariance-minus-exact-identity SNR as the
primary paired statistic.

The result also raises a source-model question. The finite injected response
matches the independently measured response very closely, which validates the
pipeline for the adopted injection PSF. It does not prove that this PSF is an
accurate model of the real planet. Plausible differences include subpixel
registration, temporal or rotational smearing, a static injection PSF versus
time-dependent AO quality, chromatic effective-PSF differences between the
star and companion, and off-axis or reduction-dependent morphology. The next
injection campaign should retain exact-PSF controls and add predeclared
broadened or perturbed PSF controls without selecting their parameters from the
planet SNR.

## Files

- `analysis.json` contains the six site values, summary statistics, paired
  differences, input hashes, and stated limitations.

## Direct follow-up

The [Stage-G setup](../klip-stage-g-planet-consistency-setup-20260926/README.md)
implements the exact-contrast, matched-aperture follow-up with additional
radius-12 positions, subpixel-phase controls, and predeclared PSF broadening.
