# Full ROC PSD study

One baseline and 90 full-image injections; no separate ROC pilot.

| Method | 0.5× | 0.75× | 1× | Evaluation null exceedances |
| --- | --- | --- | --- | --- |
| psd_hann_b0_m0.1 | 14/30 | 16/30 | 20/30 | 3/30 |
| psd_rectangular_b5_m0.3 | 13/30 | 16/30 | 21/30 | 6/30 |
| pca_b5_f1 | 10/30 | 15/30 | 17/30 | 3/30 |
| isotropic_b5 | 10/30 | 14/30 | 17/30 | 3/30 |
| identity | 9/30 | 13/30 | 16/30 | 4/30 |
| gaussian_snr | 11/30 | 14/30 | 16/30 | 3/30 |
| gaussian_raw | 10/30 | 11/30 | 17/30 | 2/30 |
| identity_snr | 12/30 | 16/30 | 19/30 | 2/30 |

All thresholds and injection contrasts were fixed from calibration before new-site scores or positives. Invalid searches remain nondetections. The 30 sites are correlated; counts do not establish independent-trial confidence intervals.

See `results.json` for individual raw photometry and conditional errors, and `protocol.json` for the fixed holdouts and controls.
