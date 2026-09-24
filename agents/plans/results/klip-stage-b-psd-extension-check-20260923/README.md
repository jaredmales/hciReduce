# KLIP Stage-B PSD extension numerical check

## Purpose

The decoupled geometry estimates noise from 11-by-11 patches but applies the
result behind response templates as large as 47 pixels. This check freezes a
finite-covariance extension that preserves the tested P4 estimator at 11 pixels
and avoids a dense 2,209-dimensional factorization.

## Numerical contract

For response width $s$, calculate every 11-pixel windowed periodogram directly
on the $(2s-1)$-square linear-lag grid: 21, 61, or 93 pixels. Average the powers,
rescale their spectral mean to the unwindowed sample variance, and apply the
fixed isotropic mixture. Because the input patches are only 11 pixels wide, the
inverse transform has zero covariance outside lags -10 through +10. This is an
explicit finite-correlation assumption rather than an extrapolation of
unmeasured long lags.

The nonnegative mixed spectrum defines a positive block-Toeplitz covariance on
the finite response grid. Apply it by zero-padded FFT convolution and solve
$Cw=t$ with preconditioned conjugate gradients. Normalize
$g=w/(t^Tw)$ to unit response. The candidate data and response remain untapered.

## Verification

The maintained
[`check_klip_stage_b_psd_extension.py`](../../scripts/check_klip_stage_b_psd_extension.py)
passes 18 numerical contracts:

- rectangular and Hann windows at isotropic mixtures 0.1, 0.3, and 1.0 exactly
  reproduce the tested P4 21-pixel spectrum, mean, variance, and finite
  11-by-11 covariance;
- iterative 11-pixel weights agree with the existing dense Cholesky solution;
- 31- and 47-pixel lag kernels are zero beyond the measured 10-pixel lag;
- the finite operators are symmetric and positive;
- mixture 1.0 gives the identity matched-filter weights and predicted variance
  at all three response sizes; and
- every iterative solve reaches the required residual.

For the deterministic test problem, the rectangular mixture-0.3 solve required
16, 18, and 18 iterations for response sizes 11, 31, and 47. Relative residuals
were $5.60\times10^{-11}$, $3.44\times10^{-11}$, and
$6.65\times10^{-11}$.

## Interpretation

The extension changes no measured 11-pixel covariance and introduces no fitted
long-range covariance. It lets the matched filter retain the measured KLIP
response tail while assigning its noise weights from the supported short-lag
PSD model. Comparing 11-, 31-, and 47-pixel response results will directly test
whether those outer response pixels add detection information under this
assumption.

The next noise-only screen must still test split stability, measured versus
predicted score variance, rectangular versus Hann estimation, radial
normalization, and the 40-pixel versus full-range inner training bands. This
numerical check establishes only the positive finite solve and its controls.

## Reproduction

```bash
python3 agents/plans/scripts/check_klip_stage_b_psd_extension.py check
```
