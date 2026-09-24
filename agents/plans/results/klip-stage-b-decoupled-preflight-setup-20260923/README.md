# KLIP Stage-B decoupled-footprint preflight result

## Question

The coupled-footprint audit found no usable 47-pixel covariance-training stamps
and inadequate 31-pixel coverage. This score-blind follow-up tested the geometry
required by a stationary PSD model: estimate the noise spectrum from 11-by-11
Welch patches while retaining an 11-, 31-, or 47-pixel response and source
exclusion.

## Result

The decoupled geometry passes at every primary radius.

| Radius | 11-pixel response | 31-pixel response | 47-pixel response |
| ---: | ---: | ---: | ---: |
| 7.5 | 20 | 40 | 40 |
| 10 | 20 | 40 | 40 |
| 12 | 20 | 40 | 40 |
| 16 | 10 | 20 | 40 |
| 20 | 10 | 20 | 20 |
| 24 | 5 | 10 | 20 |

Entries are the narrowest radial half-widths for which every selected site and
all five search pixels retain at least eight 11-pixel training patches wholly
within each detector half. Radius 6 remains unavailable because its five-pixel
candidate search is incomplete.

Even after excluding the union of five 47-pixel response footprints, the full
radial range supplies substantial training coverage:

| Radius | Minimum accepted patches | Minimum first half | Minimum second half | Rank ceiling |
| ---: | ---: | ---: | ---: | ---: |
| 7.5 | 206 | 76 | 76 | 121 |
| 10 | 183 | 61 | 60 | 121 |
| 12 | 200 | 58 | 58 | 121 |
| 16 | 193 | 47 | 48 | 121 |
| 20 | 185 | 39 | 40 | 121 |
| 24 | 186 | 39 | 40 | 121 |

The larger source exclusion costs samples but does not make 11-pixel PSD
training rank limited. At the controlling inner radii, a 40-pixel band is the
narrowest split-supported choice for both 31- and 47-pixel responses. The
next-wider full-range band remains a fixed stability control.

## Decision

Stage B will keep the 11-pixel half-overlap Welch estimator and compare its
weights behind central 11-, 31-, and 47-pixel exact responses. Direct empirical
PCA remains restricted to 11 pixels. The source exclusion always follows the
response footprint, so larger-template training cannot reuse pixels touched by
the modeled source.

The
[PSD extension check](../klip-stage-b-psd-extension-check-20260923/README.md)
freezes the remaining numerical question: how an 11-pixel periodogram defines
a positive finite covariance on the 31- and 47-pixel candidate grids.

## ROC receipt

The canonical result is
`working/roc/klip_stage_b_decoupled_preflight_20260923`. Its receipt records:

- `results.json`:
  `3d4f641763f711fb6f04740c9988776b87fbef62cf20f8e70b2ad85b5e15d2b6`;
- `results.md`:
  `0352c32bc749db520322b6802975fae2a7bd05b8c8dbec7263c43f905ec385ec`.

The detailed result contains 166,987 bytes of count distributions. Only the
compact report and receipt were transferred during the bandwidth-limited
review.
