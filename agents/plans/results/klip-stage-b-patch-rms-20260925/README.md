# KLIP Stage-B post-mean patch-RMS control

## Question and contract

Does equalizing the amplitudes of individual covariance-training patches improve
the corrected 11-pixel KLIP PSD filters? For each detector-half fit, this test:

1. subtracts the ensemble pixelwise mean from every training patch;
2. divides each residual patch by its own 121-pixel RMS;
3. applies the selected rectangular or Hann spectral window without a second
   ensemble centering; and
4. rescales the spectrum to the original post-mean training variance.

Candidate data and exact responses are never divided by a candidate-dependent
patch RMS. They retain the parent coordinate system and the fixed seven-pixel
planet mask. The test is paired query by query with the corrected parent.

The frozen covariance families are Hann/mixing-0.1, the leading data-selected
PSD, and rectangular/mixing-0.3, the mandatory P4 prior. Raw coordinates use
both the narrow and full training bands. The leading strict radial coordinates
use their narrowest split-supported band.

## Raw-coordinate result

| Band | Method | Parent variance | Patch-RMS variance | Ratio | + Fitted mean | Opposite-half variance | Paired score correlation | Median abs. score change |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| narrow | rectangular/mixing-0.3 | 2.867 | 2.865 | 0.999 | 2.861 | 1.344 | 0.999833 | 0.018 |
| narrow | Hann/mixing-0.1 | 2.477 | 2.482 | 1.002 | 2.477 | 1.136 | 0.999745 | 0.023 |
| full | rectangular/mixing-0.3 | 3.623 | 3.655 | 1.009 | 3.646 | 1.252 | 0.999878 | 0.024 |
| full | Hann/mixing-0.1 | 3.083 | 3.104 | 1.007 | 3.100 | 1.060 | 0.999861 | 0.023 |

The primary narrow-band result is unchanged at the 0.2% level. Full-band
variance is slightly worse. Median maximum-to-minimum patch-RMS ratios are only
1.49 in the narrow bands and 1.79 in the full bands, and the fitted-mean control
does not improve calibration.

## Radial-coordinate interaction

| Method | Radial parent variance | + Patch-RMS variance | Ratio | + Fitted mean | Opposite-half variance | Paired score correlation | Median abs. score change |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| rectangular/mixing-0.3 | 2.412 | 2.403 | 0.996 | 2.402 | 1.298 | 0.999837 | 0.018 |
| Hann/mixing-0.1 | 2.092 | 2.089 | 0.999 | 2.088 | 1.091 | 0.999764 | 0.020 |

Patch RMS is likewise neutral after radial standardization. At radius 12, the
radial Hann variance changes from 3.153 to 3.168; the controlling-radius median
improvement to 2.089 is only 0.1%. At larger radii the Hann variance generally
increases. The paired scores remain almost identical.

## Narrow-band variance by radius

| Radius | Raw rectangular parent | Raw + patch RMS | Raw Hann parent | Raw + patch RMS | Radial Hann parent | Radial Hann + patch RMS |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 7.5 | 2.068 | 2.063 | 1.802 | 1.796 | 1.653 | 1.639 |
| 10 | 2.867 | 2.865 | 2.477 | 2.482 | 2.092 | 2.089 |
| 12 | 4.611 | 4.613 | 3.897 | 3.930 | 3.153 | 3.168 |
| 16 | 2.183 | 2.222 | 1.663 | 1.720 | 1.554 | 1.605 |
| 20 | 1.228 | 1.234 | 1.006 | 1.017 | 1.123 | 1.137 |
| 24 | 1.264 | 1.271 | 0.958 | 0.981 | 0.958 | 0.980 |

## Decision

Do not retain per-patch RMS normalization for KLIP. It neither improves
candidate calibration nor changes the fitted scores materially in raw or radial
coordinates. Keep radial-standardized 11-pixel Hann/mixing-0.1 as the leading
PSD arm and raw 11-pixel rectangular/mixing-0.3 as the mandatory P4 prior. The
next Stage-B discriminator is the direct 11-pixel diagonal/PCA covariance grid;
the larger-support radial-mean control remains open.

## Reproduction and verification

The raw-coordinate runner is
[`run_klip_stage_b_local_patch_rms.py`](../../scripts/run_klip_stage_b_local_patch_rms.py),
and the radial interaction is
[`run_klip_stage_b_local_radial_patch_rms.py`](../../scripts/run_klip_stage_b_local_radial_patch_rms.py).
They completed 1,440 and 720 policy/query records, respectively. Both replay
the verified bundle, parent receipts, frozen geometry, candidate masks, and
detector-half splits. Tests verify the normalization definition, preservation
of the raw fitted mean and variance scale, and unit physical source response
through simultaneous scaling and masked covariance solves.

Compact values are in [`primary.json`](primary.json), source-product receipts
are in [`verification.json`](verification.json), and the figures show the
[raw-coordinate](raw-coordinate-comparison.png) and
[radial-coordinate](radial-coordinate-comparison.png) comparisons. Full record
tables remain under ignored `working/local` data.

## Follow-up: direct diagonal/PCA covariance

The completed [direct-covariance screen](../klip-stage-b-direct-covariance-20260925/README.md)
finds no replacement for the PSD shortlist. The best direct model uses every
estimable PCA mode and floor 1.0, but its radial primary variance is 2.886 and
split-weight cosine is 0.730, versus 2.092 and 0.990 for radial Hann/mixing-0.1.
The remaining Stage-B arm is the larger-response radial-mean control.
