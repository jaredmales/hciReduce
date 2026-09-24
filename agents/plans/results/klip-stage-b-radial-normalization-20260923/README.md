# KLIP Stage-B strict radial-normalization screen

This mode-200 test changes only the input coordinates of the 11-pixel narrow-band PSD arm. For each fixed site, a 3.6-pixel-bin radial variance profile excludes the known planet and the union of all five complete candidate footprints. The same leave-site-out profile scales native training pixels, candidate data, and exact responses by standard deviation before interpolation and filtering.

## Controlling radii

| Method | Raw score variance | Normalized score variance | Normalized/raw | Center variance | Opposite-half variance | Split physical-weight cosine | Within ±1σ |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| identity | 6.164 | 5.024 | 0.818 | 5.135 | 3.923 | 1.000 | 0.267 |
| rectangular_m0.1 | 2.340 | 1.956 | 0.836 | 2.345 | 1.131 | 0.985 | 0.467 |
| rectangular_m0.3 | 2.633 | 2.191 | 0.832 | 2.611 | 1.304 | 0.989 | 0.425 |
| hann_m0.1 | 2.232 | 1.856 | 0.832 | 2.148 | 1.090 | 0.990 | 0.467 |
| hann_m0.3 | 2.508 | 2.072 | 0.826 | 2.410 | 1.275 | 0.993 | 0.442 |

## Normalization geometry

| Response support | Supported fixed sites | Minimum pixels in any radial bin |
| ---: | ---: | ---: |
| 11 | 72/72 | 21 |
| 31 | 22/72 | 0 |
| 47 | 0/72 | 0 |

Strict candidate-excluded profiles are supported for every 11-pixel site, with at least 21 pixels in every bin. They fail at 50 of 72 31-pixel sites and all 47-pixel sites because the larger held-out response footprints cover complete inner annuli. Those response supports were not normalized using a profile that reads the candidate or extrapolates across an unsupported radial interval.

The innermost 0--3.6-pixel KLIP bin is a structural zero: all 44 baseline pixels and all selected exact-response occurrences are bitwise zero. It receives the first positive-bin scale only after this gate; scaling a zero data/response coordinate cannot change its filter contribution.

## Result

Radial standardization is a useful partial correction, not a complete
calibration. For the leading Hann/mixing-0.1 arm, the median candidate score
variance over radii 7.5, 10, and 12 falls from 2.23 to 1.86, a median ratio of
0.832. The center-only median falls from 2.41 to 2.15, and the fraction of
candidate scores within one conditional sigma rises from 0.425 to 0.467.
Rectangular/mixing-0.3 improves similarly, from 2.63 to 2.19.

The correction does not remove the candidate-versus-generic-patch discrepancy.
Hann/mixing-0.1 opposite-half variance changes only from 1.13 to 1.09, while
its candidate variance remains 1.77, 1.86, and 3.10 at radii 7.5, 10, and 12.
The radius-12 value was 3.84 before normalization, so the radial scale explains
part, but not most, of that mismatch. At radius 20 the Hann candidate variance
slightly worsens from 1.07 to 1.19, and radius 24 is unchanged at 0.96.

Normalization does not buy its improvement through unstable weights. The
primary Hann/mixing-0.1 median physical-weight cosine changes from 0.991 to
0.990. The other PSD variants likewise retain their raw split stability.

## Geometry decision

Keep the strict candidate-excluded definition. It prevents the normalization
profile from reading the location whose score it calibrates, and it can be
frozen from the signal-free baseline for later positive-image work. This
profile is valid on common support for all 72 selected 11-pixel sites.

Do not extend this exact profile rule to the larger response footprints. Only
22 of 72 31-pixel sites have all radial bins supported, and none of the 72
47-pixel sites do. Filling those intervals from candidate pixels or endpoint
extrapolation would answer a different question. A larger-response
normalization arm would require an independently estimated baseline profile or
an explicitly cross-fitted angular policy.

The exact-zero KLIP core is handled without a variance floor. All 44 baseline
pixels inside 3.6 pixels and all 672 selected exact-response occurrences are
bitwise zero. Assigning the first positive-bin scale to those coordinates is
algebraically inert because both the data and response are zero there.

## Decision

Retain radial-standardized Hann/mixing-0.1 as the leading 11-pixel PSD arm and
raw rectangular/mixing-0.3 as the mandatory P4 prior. Conditional candidate
variance remains too large for promotion. The next Stage-B discriminator is
the prespecified post-ensemble-mean patch-RMS control on the supported
11-pixel geometry. Direct PCA/diagonal comparisons and the larger-support mean
question remain open.

## Reproduction and provenance

The maintained runner is
[`run_klip_stage_b_local_radial_normalization.py`](../../scripts/run_klip_stage_b_local_radial_normalization.py).
It verified the compact bundle and complete raw-screen receipt, evaluated 1,800
policy/query records in 8.84 seconds, and passed its masking, physical-amplitude,
structural-core, and product-hash gates. Compact values are in
[`primary.json`](primary.json), receipts are in
[`verification.json`](verification.json), and the comparison is shown in
[`comparison.png`](comparison.png). Complete profiles and records remain under
ignored `working/local` data.
