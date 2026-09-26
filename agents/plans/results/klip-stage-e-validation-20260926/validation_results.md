# KLIP Stage-E fresh validation

All values use the immutable Stage-D methods, thresholds, sites, contrasts, and baseline-derived weights. Mode 200 is primary; no maximization over methods or modes was performed.

## Mode-200 pooled comparison

| Target | Method | Recovered / 36 | Mean maximum SNR | Mean throughput |
| ---: | :--- | ---: | ---: | ---: |
| 3 | native | 23 | 2.9472 | 1.0053 |
| 3 | gaussian_fwhm2p4 | 26 | 3.4244 | 1.0078 |
| 3 | gaussian_fwhm3p6 | 23 | 3.1239 | 1.0089 |
| 3 | exact_identity | 28 | 3.7296 | 1.0064 |
| 3 | sparse_identity | 29 | 3.6716 | 1.0065 |
| 3 | raw_rectangular_m0p3 | 26 | 3.6608 | 1.0057 |
| 3 | radial_hann_m0p1_trunc0p75 | 27 | 3.7255 | 1.0059 |
| 5 | native | 34 | 4.4075 | 1.0006 |
| 5 | gaussian_fwhm2p4 | 35 | 5.4820 | 1.0025 |
| 5 | gaussian_fwhm3p6 | 35 | 4.9758 | 1.0033 |
| 5 | exact_identity | 35 | 5.8870 | 1.0014 |
| 5 | sparse_identity | 36 | 5.7840 | 1.0016 |
| 5 | raw_rectangular_m0p3 | 35 | 5.7607 | 1.0009 |
| 5 | radial_hann_m0p1_trunc0p75 | 36 | 5.8651 | 1.0012 |
| 7 | native | 36 | 5.8672 | 0.9918 |
| 7 | gaussian_fwhm2p4 | 36 | 7.4345 | 0.9934 |
| 7 | gaussian_fwhm3p6 | 36 | 6.7536 | 0.9937 |
| 7 | exact_identity | 36 | 7.9871 | 0.9932 |
| 7 | sparse_identity | 36 | 7.8734 | 0.9934 |
| 7 | raw_rectangular_m0p3 | 36 | 7.8168 | 0.9931 |
| 7 | radial_hann_m0p1_trunc0p75 | 36 | 7.9349 | 0.9934 |

## Frozen covariance acceptance gate

| Candidate | Accepted | Held-out (candidate / G3.6) | Recovery totals (candidate / G3.6) | SNR gate | Throughput | Support |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| raw_rectangular_m0p3 | yes | 3 / 3 | 97 / 94 | pass | pass | pass |
| radial_hann_m0p1_trunc0p75 | yes | 3 / 3 | 99 / 94 | pass | pass | pass |

The other seven KL modes are retained in the CSV and JSON as separate generalization checks. They are not combined into a maximum statistic.
