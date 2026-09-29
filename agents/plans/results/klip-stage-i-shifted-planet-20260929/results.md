# KLIP Stage-I shifted-template planet closure

All values below are mode-200 SNR differences at the nearest native
pixel to the common optimized-phase source position. Injection values
are means +/- sample standard deviations across twelve sites.

| Comparison | Injections | Range | Planet | Planet deviation |
| :--- | ---: | :--- | ---: | ---: |
| Shifted identity - integer identity | +0.0671 +/- 0.2753 | [-0.3246, +0.5612] | +0.2680 | +0.73 SD |
| Shifted raw covariance - integer raw covariance | +0.0580 +/- 0.2848 | [-0.4239, +0.5430] | +0.1293 | +0.25 SD |
| Shifted radial covariance - integer radial covariance | +0.0948 +/- 0.2869 | [-0.3775, +0.5743] | +0.1654 | +0.25 SD |
| Shifted raw covariance - shifted identity | -0.2234 +/- 0.2196 | [-0.5226, +0.1944] | -0.6921 | -2.13 SD |
| Shifted radial covariance - shifted identity | -0.1715 +/- 0.2648 | [-0.5766, +0.2412] | -0.7004 | -2.00 SD |
| Gaussian 3.6 - shifted identity | -0.0281 +/- 1.2056 | [-1.5064, +2.0770] | +0.5676 | +0.49 SD |
| Gaussian 3.6 - Gaussian 2.4 | -0.2631 +/- 0.2599 | [-0.6646, +0.1163] | -0.1827 | +0.31 SD |

## Absolute fixed-pixel SNR

| Method | Injection SNR | Planet SNR | Planet deviation |
| :--- | ---: | ---: | ---: |
| gaussian_fwhm2p4 | 4.8621 +/- 1.1576 | 5.2696 | +0.35 SD |
| gaussian_fwhm3p6 | 4.5990 +/- 1.3375 | 5.0869 | +0.36 SD |
| integer_exact_identity | 4.5600 +/- 1.0684 | 4.2514 | -0.29 SD |
| shifted_exact_identity | 4.6271 +/- 1.0872 | 4.5193 | -0.10 SD |
| integer_raw_rectangular_m0p3 | 4.3457 +/- 0.9195 | 3.6979 | -0.70 SD |
| shifted_raw_rectangular_m0p3 | 4.4037 +/- 0.9675 | 3.8273 | -0.60 SD |
| integer_radial_hann_m0p1_trunc0p75 | 4.3608 +/- 0.8930 | 3.6535 | -0.79 SD |
| shifted_radial_hann_m0p1_trunc0p75 | 4.4556 +/- 0.9405 | 3.8189 | -0.68 SD |

The report is descriptive. The predeclared paired differences and
exchangeable ranks determine whether response registration resolves the
fixed-pixel planet discrepancy.
