hciAnalyze needs to use the Mawet et al 2014 small-sample statistical correction to SNR.

1. first we will add the lambdaD parameter, which is the pixel scale in pixels per (lambda/D).  In configuration this should be a top-level (no section) setting.

2. The default value of lambdaD should be 2.5, and hciAnalyze will print a warning that it was not set. (It should not warn if it is set to 2.5 in config)

3. At each pixel of radius r in pixels, calculate n2 = 2*pi*r/lambdaD - 1.

4. Correct the S/N of that pixel with 1/sqrt(1 + 1/n2).

## Implementation status

- [x] Add the top-level `lambdaD` setting with a `2.5` default and an omission warning.
- [x] Convert pixel radius to lambda/D radius and calculate `n2 = 2*pi*r/lambdaD - 1`.
- [x] Apply the requested correction to every SNR-map pixel and record it in FITS provenance.
- [x] Add focused configuration, warning, numerical-correction, output-header, and documentation coverage.

Verification: all 26 test targets and `make docs` pass. Applying the correction to the external-PSF-filtered
`finim_0046.fits` example with `lambdaD=2.5` changes the best reported S/N from `7.13242` to `7.00900`, consistent
with the approximately `0.9830` correction factor at the companion separation.
