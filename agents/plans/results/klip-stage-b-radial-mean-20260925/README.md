# KLIP Stage-B support-independent radial-mean control

## Question and contract

Can a radial background model make the larger 31- or 47-pixel exact responses
better calibrated than the 11-pixel response? This mode-200 screen fits one
one-pixel annular mean profile to the signal-free final image. It excludes
native pixels within 7.5 pixels of the known planet, matching the `R + 0.5`
convention in `hciAnalyze`, and linearly interpolates the same profile across
each candidate stamp.

This is a mean-model control. The corrected fixed planet mask, exact response,
raw 11-pixel Welch PSD estimate, isotropic mixing, unit-response weights, and
conditional uncertainty remain identical to the corrected raw parent. The
runner checks every recomputed raw score against that parent before applying
the mean correction. It evaluates identity, rectangular/mixing-0.3, and
Hann/mixing-0.1 at all six radii, all three response supports, and both frozen
training bands.

## Controlling-radius result

The table gives the median candidate-score variance over radii 7.5, 10, and 12
for the primary narrow band. A useful correction would move variance toward
one and materially change the candidate statistic.

| Response | Method | Raw variance | Radial-mean variance | Ratio | Median absolute score change | Raw/radial correlation |
| ---: | :--- | ---: | ---: | ---: | ---: | ---: |
| 11 | identity | 6.512 | 6.444 | 0.989 | 0.038 | 0.99992 |
| 11 | rectangular/mixing-0.3 | 2.867 | 2.832 | 0.988 | 0.021 | 0.99984 |
| 11 | Hann/mixing-0.1 | 2.477 | 2.446 | 0.987 | 0.022 | 0.99980 |
| 31 | identity | 7.131 | 7.042 | 0.987 | 0.015 | 0.99992 |
| 31 | rectangular/mixing-0.3 | 3.324 | 3.276 | 0.986 | 0.016 | 0.99982 |
| 31 | Hann/mixing-0.1 | 2.875 | 2.831 | 0.985 | 0.013 | 0.99981 |
| 47 | identity | 6.879 | 6.755 | 0.982 | 0.022 | 0.99992 |
| 47 | rectangular/mixing-0.3 | 3.306 | 3.256 | 0.985 | 0.019 | 0.99985 |
| 47 | Hann/mixing-0.1 | 2.698 | 2.654 | 0.984 | 0.016 | 0.99984 |

The full-band ratios are similarly small, 0.982--0.989. For the 11-pixel Hann
arm, subtracting the same profile from disjoint opposite-half patches changes
the median measured/predicted variance from 1.133 to 1.116. The fitted profile
has 4--380 planet-excluded samples per annulus and means from -0.0286 to 0.0229
in final-image units.

## Radius dependence

The leading narrow-band Hann result remains strongly radius dependent.

| Radius | 11-pixel raw | 11-pixel radial mean | 31-pixel raw | 31-pixel radial mean | 47-pixel raw | 47-pixel radial mean |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 7.5 | 1.802 | 1.802 | 1.779 | 1.756 | 1.680 | 1.654 |
| 10 | 2.477 | 2.446 | 2.875 | 2.831 | 2.698 | 2.654 |
| 12 | 3.897 | 3.885 | 5.083 | 5.057 | 4.755 | 4.730 |
| 16 | 1.663 | 1.675 | 2.300 | 2.307 | 2.858 | 2.878 |
| 20 | 1.006 | 1.005 | 1.255 | 1.249 | 1.261 | 1.254 |
| 24 | 0.958 | 0.955 | 1.004 | 1.002 | 1.190 | 1.189 |

The correction does not repair radius 12, and the larger responses remain
worse than 11 pixels there. The mean projection itself has variance of only
about 0.001--0.002 conditional-sigma squared at the controlling radii.

## Decision

Do not add radial-mean subtraction to the shortlist and do not promote a larger
response support from this test. The change is too small to explain the
candidate calibration failure. Stage B is now closed with radial-standardized
11-pixel Hann/mixing-0.1 as the leading data-selected arm and raw 11-pixel
rectangular/mixing-0.3 as the mandatory P4 prior. Advance those two covariance
families, along with the permanent references, to the frozen development
injection stage.

## Reproduction and verification

The maintained runner is
[`run_klip_stage_b_local_radial_mean.py`](../../scripts/run_klip_stage_b_local_radial_mean.py).
It completed 6,480 policy/query records in 157.2 seconds. Deterministic checks
cover constant-profile recovery, support-independent evaluation, and exact
replay of every corrected raw parent score.

Compact primary and per-radius Hann results are in [`primary.json`](primary.json),
source-product hashes are in [`verification.json`](verification.json), and the
comparison is shown in [`comparison.png`](comparison.png). The 8.9-MB record
table remains under ignored `working/local` data.
