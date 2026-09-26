# KLIP Stage-D frozen validation policy

This receipt was written before any held-out-null score or validation positive was opened. It freezes mode 200 as primary, baseline-derived weights, the seven-method shortlist, the existing method-specific maximum-of-20 thresholds, and the Stage-E acceptance rule.

## Frozen methods

- `native`
- `gaussian_fwhm2p4`
- `gaussian_fwhm3p6`
- `exact_identity`
- `sparse_identity`
- `raw_rectangular_m0p3`
- `radial_hann_m0p1_trunc0p75`

The confirmatory covariance candidates are `raw_rectangular_m0p3` and `radial_hann_m0p1_trunc0p75`. The primary comparison is paired against `gaussian_fwhm3p6`; `gaussian_fwhm2p4` is reported as the stronger smoothing control.

## Mode-200 threshold audit

| Radius | Method | Threshold | Second largest | Gap | Bootstrap 2.5--97.5% | Block-delete 2.5--97.5% |
| ---: | :--- | ---: | ---: | ---: | :--- | :--- |
| 7.5 | native | 3.8744 | 3.6039 | 0.2705 | 2.3392--3.8744 | 2.8722--3.8744 |
| 7.5 | gaussian_fwhm2p4 | 3.0493 | 2.7668 | 0.2825 | 2.2118--3.0493 | 2.7668--3.0493 |
| 7.5 | gaussian_fwhm3p6 | 2.4984 | 2.4552 | 0.0432 | 2.0670--2.4984 | 2.2287--2.4984 |
| 7.5 | exact_identity | 2.3013 | 1.8685 | 0.4329 | 1.4536--2.3013 | 1.8685--2.3013 |
| 7.5 | sparse_identity | 2.1172 | 2.0565 | 0.0608 | 1.7850--2.1172 | 2.0565--2.1172 |
| 7.5 | raw_rectangular_m0p3 | 2.6729 | 2.0954 | 0.5775 | 1.9344--2.6729 | 2.0461--2.6729 |
| 7.5 | radial_hann_m0p1_trunc0p75 | 2.6740 | 2.1188 | 0.5552 | 1.4494--2.6740 | 2.1188--2.6740 |
| 10 | native | 2.1332 | 2.0603 | 0.0729 | 1.6841--2.1332 | 2.0603--2.1332 |
| 10 | gaussian_fwhm2p4 | 2.9003 | 2.8824 | 0.0179 | 1.3897--2.9003 | 2.0987--2.9003 |
| 10 | gaussian_fwhm3p6 | 3.1761 | 2.8372 | 0.3390 | 2.0315--3.1761 | 1.9883--3.1761 |
| 10 | exact_identity | 2.2241 | 2.0626 | 0.1615 | 1.4452--2.2241 | 1.8114--2.2241 |
| 10 | sparse_identity | 2.0022 | 1.9779 | 0.0243 | 1.4226--2.0022 | 1.9779--2.0022 |
| 10 | raw_rectangular_m0p3 | 2.4166 | 1.7672 | 0.6494 | 1.2430--2.4166 | 1.6808--2.4166 |
| 10 | radial_hann_m0p1_trunc0p75 | 2.2938 | 1.9274 | 0.3664 | 1.1889--2.2938 | 1.9274--2.2938 |
| 12 | native | 2.1029 | 1.8295 | 0.2734 | 1.5639--2.1029 | 1.8295--2.1029 |
| 12 | gaussian_fwhm2p4 | 2.4070 | 1.9113 | 0.4957 | 1.6656--2.4070 | 1.7108--2.4070 |
| 12 | gaussian_fwhm3p6 | 2.6219 | 2.5331 | 0.0888 | 1.7151--2.6219 | 1.7670--2.6219 |
| 12 | exact_identity | 2.5614 | 1.5193 | 1.0421 | 1.5053--2.5614 | 1.5193--2.5614 |
| 12 | sparse_identity | 2.3612 | 1.7851 | 0.5760 | 1.4262--2.3612 | 1.7851--2.3612 |
| 12 | raw_rectangular_m0p3 | 2.3355 | 1.7008 | 0.6346 | 1.5358--2.3355 | 1.6162--2.3355 |
| 12 | radial_hann_m0p1_trunc0p75 | 2.2572 | 1.6665 | 0.5907 | 1.4543--2.2572 | 1.6215--2.2572 |
| 16 | native | 2.6226 | 2.1725 | 0.4502 | 1.5956--2.6226 | 2.1725--2.6226 |
| 16 | gaussian_fwhm2p4 | 1.9763 | 1.8743 | 0.1019 | 1.2458--1.9763 | 1.8743--1.9763 |
| 16 | gaussian_fwhm3p6 | 2.4044 | 1.6066 | 0.7978 | 1.4064--2.4044 | 1.4675--2.4044 |
| 16 | exact_identity | 1.9432 | 1.4141 | 0.5291 | 1.2738--1.9432 | 1.3556--1.9432 |
| 16 | sparse_identity | 1.8589 | 1.5973 | 0.2616 | 1.3295--1.8589 | 1.3664--1.8589 |
| 16 | raw_rectangular_m0p3 | 2.3748 | 1.6103 | 0.7645 | 1.4647--2.3748 | 1.6103--2.3748 |
| 16 | radial_hann_m0p1_trunc0p75 | 2.2088 | 1.6105 | 0.5983 | 1.5461--2.2088 | 1.6105--2.2088 |
| 20 | native | 2.4598 | 2.3727 | 0.0871 | 1.2888--2.4598 | 2.3727--2.4598 |
| 20 | gaussian_fwhm2p4 | 2.3799 | 2.1159 | 0.2641 | 1.5222--2.3799 | 2.1159--2.3799 |
| 20 | gaussian_fwhm3p6 | 2.0300 | 1.7077 | 0.3223 | 1.2399--2.0300 | 1.7077--2.0300 |
| 20 | exact_identity | 2.5392 | 2.0613 | 0.4779 | 1.4291--2.5392 | 2.0613--2.5392 |
| 20 | sparse_identity | 2.3724 | 2.0468 | 0.3257 | 1.5186--2.3724 | 2.0468--2.3724 |
| 20 | raw_rectangular_m0p3 | 2.6974 | 1.9509 | 0.7464 | 1.4459--2.6974 | 1.9509--2.6974 |
| 20 | radial_hann_m0p1_trunc0p75 | 2.8395 | 2.2506 | 0.5889 | 1.2344--2.8395 | 2.2506--2.8395 |
| 24 | native | 3.4460 | 2.8858 | 0.5602 | 2.0577--3.4460 | 2.8858--3.4460 |
| 24 | gaussian_fwhm2p4 | 3.7124 | 2.2416 | 1.4709 | 1.4448--3.7124 | 2.2416--3.7124 |
| 24 | gaussian_fwhm3p6 | 3.3157 | 1.8428 | 1.4730 | 1.3157--3.3157 | 1.8428--3.3157 |
| 24 | exact_identity | 4.0673 | 2.1281 | 1.9392 | 1.1547--4.0673 | 2.1281--4.0673 |
| 24 | sparse_identity | 3.9431 | 2.3378 | 1.6053 | 1.1242--3.9431 | 2.3378--3.9431 |
| 24 | raw_rectangular_m0p3 | 4.2160 | 1.9829 | 2.2332 | 1.2201--4.2160 | 1.9829--4.2160 |
| 24 | radial_hann_m0p1_trunc0p75 | 4.0653 | 1.8441 | 2.2212 | 1.2074--4.0653 | 1.8441--4.0653 |

The 20 calibration locations are spatially correlated. Their maximum is retained for direct continuity with the P4 tests and is not assigned a formal false-alarm probability. Six separately assigned null locations per radius may now be exposed; they may reject a method but cannot alter this policy.

## Validation gate

A covariance candidate must have no more held-out exceedances than Gaussian 3.6, match or exceed its recovery at every source level with a strictly larger total, pass the frozen paired-SNR bootstrap rule, retain the frozen throughput ranges, and retain all common sites. Other KL modes are reported separately and are not maximized.
