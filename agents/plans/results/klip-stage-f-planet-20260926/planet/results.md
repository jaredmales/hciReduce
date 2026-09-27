# KLIP Stage-F known-planet closure

This is a descriptive single-planet endpoint performed after the Stage-E policy and validation result were immutable. It does not select or rescue a method.

## Mode-200 closure

| Method | Validation SNR-3 recovery | Held-out exceedances | Planet nearest-pixel SNR | Planet aperture-max SNR | Peak row, column | Peak offset (pixels) | Peak contrast estimate |
| :--- | :---: | :---: | ---: | ---: | :---: | ---: | ---: |
| Native | 23 / 36 | 1 | 3.4662 | 5.1536 | 76, 62 | 0.841 | 0.00565636 |
| Gaussian 2.4 | 26 / 36 | 3 | 5.0641 | 5.4186 | 76, 62 | 0.841 | 0.00532435 |
| Gaussian 3.6 | 23 / 36 | 3 | 5.6496 | 5.6496 | 75, 62 | 0.213 | 0.00572327 |
| Exact identity | 28 / 36 | 1 | 4.1813 | 4.2375 | 76, 61 | 1.204 | 0.00428321 |
| Sparse identity | 29 / 36 | 1 | 4.2461 | 4.2461 | 75, 62 | 0.213 | 0.0041712 |
| Exact response LPF 1.8 | not in Stage E | not in Stage E | 4.5358 | 4.5358 | 75, 62 | 0.213 | 0.00434282 |
| Exact response LPF 2.7 | not in Stage E | not in Stage E | 4.8662 | 4.8662 | 75, 62 | 0.213 | 0.00466787 |
| Exact identity, fitted mean | not in Stage E | not in Stage E | 4.1930 | 4.2794 | 76, 61 | 1.204 | 0.00432978 |
| Raw rectangular PSD | 26 / 36 | 3 | 3.2291 | 3.6937 | 76, 61 | 1.204 | 0.00380313 |
| Radial Hann, truncation 0.75 | 27 / 36 | 3 | 3.1773 | 3.6461 | 76, 61 | 1.204 | 0.00378083 |

## Planet aperture-maximum SNR by KL mode

| Method | 125 | 150 | 175 | 200 | 225 | 250 | 300 | 350 |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Native | 5.0183 | 4.9844 | 4.9304 | 5.1536 | 5.1840 | 5.2919 | 5.1496 | 5.0958 |
| Gaussian 2.4 | 5.4179 | 5.3423 | 5.3823 | 5.4186 | 5.5946 | 5.6071 | 5.5033 | 5.4666 |
| Gaussian 3.6 | 5.5507 | 5.6304 | 5.6556 | 5.6496 | 5.7615 | 5.7882 | 5.8313 | 5.7935 |
| Exact identity | 4.5218 | 4.4282 | 4.2993 | 4.2375 | 4.3926 | 4.3164 | 4.3166 | 4.3269 |
| Sparse identity | 4.4611 | 4.3074 | 4.2395 | 4.2461 | 4.3968 | 4.3578 | 4.3581 | 4.3640 |
| Exact response LPF 1.8 | 4.7737 | 4.6501 | 4.5599 | 4.5358 | 4.7186 | 4.6750 | 4.6558 | 4.6608 |
| Exact response LPF 2.7 | 5.0656 | 4.9353 | 4.9046 | 4.8662 | 5.0368 | 5.0142 | 4.9905 | 4.9776 |
| Exact identity, fitted mean | 4.5592 | 4.4670 | 4.3411 | 4.2794 | 4.4350 | 4.3540 | 4.3569 | 4.3391 |
| Raw rectangular PSD | 3.8352 | 3.8338 | 3.7409 | 3.6937 | 3.8713 | 3.8230 | 3.8066 | 3.7765 |
| Radial Hann, truncation 0.75 | 3.7757 | 3.7635 | 3.6791 | 3.6461 | 3.8117 | 3.7577 | 3.7489 | 3.6872 |

## Verification

- Analysis aperture: 39 native pixels within 3.5 pixels of the `working/analyze.conf` center.
- Maximum production-versus-independent annular-SNR error: `9.536743e-07`.
- Generic annular maps use the nearest frozen Stage-C radial policy. Aperture weights are fit only on the signal-free baseline and exclude the complete analysis aperture from training.
- Exact and sparse response identity, response-smoothed identity, fitted-mean identity, and both covariance filters retain the planet pixels in candidate support.

The injection and held-out columns remain the method-selection evidence. The planet columns are a consistency check on one real source.
