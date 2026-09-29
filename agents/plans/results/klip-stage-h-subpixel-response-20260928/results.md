# KLIP Stage-H subpixel-response result

Twelve planet-like phase responses at offset `(-0.2771, +0.4860)` pixels are compared at mode 200. Values are means +/- sample standard deviations.

The paired central response is `(positive-negative)/(2*contrast)`. The positive
target is the actual one-sided Stage-G injection response. Cosine is the retained
optimal matched-filter SNR fraction under the stated covariance metric.

## Paired central-response target

| Covariance metric | Model | Cosine | Projection scale | Best-scaled residual |
| :--- | :--- | ---: | ---: | ---: |
| Identity | Current integer | 0.925752 +/- 0.034646 | 0.91848 +/- 0.02944 | 0.3670 +/- 0.0888 |
| Identity | Cubic shifted | 0.993552 +/- 0.003880 | 0.99039 +/- 0.01215 | 0.1072 +/- 0.0382 |
| Identity | Paired exact | 1.000000 +/- 0.000000 | 1.00000 +/- 0.00000 | 0.0000 +/- 0.0000 |
| Raw rectangular | Current integer | 0.920476 +/- 0.035519 | 0.91130 +/- 0.03068 | 0.3803 +/- 0.0868 |
| Raw rectangular | Cubic shifted | 0.990793 +/- 0.006183 | 0.98881 +/- 0.01783 | 0.1269 +/- 0.0488 |
| Raw rectangular | Paired exact | 1.000000 +/- 0.000000 | 1.00000 +/- 0.00000 | 0.0000 +/- 0.0000 |
| Radial Hann | Current integer | 0.919875 +/- 0.034450 | 0.91108 +/- 0.03018 | 0.3825 +/- 0.0840 |
| Radial Hann | Cubic shifted | 0.989325 +/- 0.007182 | 0.98714 +/- 0.01931 | 0.1368 +/- 0.0520 |
| Radial Hann | Paired exact | 1.000000 +/- 0.000000 | 1.00000 +/- 0.00000 | 0.0000 +/- 0.0000 |

## Actual positive-response target

| Covariance metric | Model | Cosine | Projection scale | Best-scaled residual |
| :--- | :--- | ---: | ---: | ---: |
| Identity | Current integer | 0.925365 +/- 0.034303 | 0.91793 +/- 0.02922 | 0.3683 +/- 0.0874 |
| Identity | Cubic shifted | 0.993342 +/- 0.003988 | 0.99004 +/- 0.01545 | 0.1091 +/- 0.0384 |
| Identity | Paired exact | 0.999821 +/- 0.000143 | 0.99966 +/- 0.00695 | 0.0178 +/- 0.0067 |
| Raw rectangular | Current integer | 0.919788 +/- 0.034841 | 0.91127 +/- 0.02964 | 0.3826 +/- 0.0843 |
| Raw rectangular | Cubic shifted | 0.990301 +/- 0.006563 | 0.98906 +/- 0.01932 | 0.1301 +/- 0.0505 |
| Raw rectangular | Paired exact | 0.999480 +/- 0.000493 | 1.00023 +/- 0.00664 | 0.0296 +/- 0.0133 |
| Radial Hann | Current integer | 0.918823 +/- 0.033557 | 0.91098 +/- 0.02919 | 0.3857 +/- 0.0807 |
| Radial Hann | Cubic shifted | 0.988500 +/- 0.007881 | 0.98736 +/- 0.02114 | 0.1415 +/- 0.0550 |
| Radial Hann | Paired exact | 0.999084 +/- 0.000920 | 1.00013 +/- 0.00692 | 0.0388 +/- 0.0188 |

## Interpretation guide

- If the shifted model approaches the paired-exact model, the present mismatch is predominantly registration and a cheap phase-shift model is sufficient.
- If a substantial residual remains, the KLIP response shape itself depends on phase and the next implementation should measure a fractional-phase response grid.
- The paired-exact versus positive rows bound finite-amplitude asymmetry that neither integer nor shifted registration can remove.
