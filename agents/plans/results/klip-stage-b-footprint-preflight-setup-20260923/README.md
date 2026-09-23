# KLIP Stage-B footprint preflight result

## Question

The accepted KLIP response contains 47-by-47 pixels, or 2,209 template
components. This score-blind preflight tested whether the response footprint
could also serve as the covariance-training patch size at the planned radii.
It measured response-energy capture, common search support, and half-overlap
training coverage without fitting a covariance or reading a baseline score.

## Response footprint

All three supports retain the same five-pixel search geometry at the six
primary radii. Radius 6 remains unsupported for the five-pixel search.

| Radius | Eligible 11-pixel centers | Eligible 31-pixel centers | Eligible 47-pixel centers |
| ---: | ---: | ---: | ---: |
| 7.5 | 42 | 42 | 42 |
| 10 | 44 | 44 | 44 |
| 12 | 60 | 60 | 60 |
| 16 | 91 | 91 | 91 |
| 20 | 113 | 113 | 113 |
| 24 | 160 | 160 | 160 |

Across the complete search field, 9,428 locations have complete 11-pixel
responses in all eight modes, 5,380 have complete 31-pixel responses, and
2,852 have complete 47-pixel responses.

Mode-200 median response energy retained relative to 47 pixels is:

| Radius | 11 pixels | 31 pixels |
| ---: | ---: | ---: |
| 7.5 | 0.7878 | 0.9477 |
| 10 | 0.8679 | 0.9198 |
| 12 | 0.9056 | 0.9528 |
| 16 | 0.9560 | 0.9935 |
| 20 | 0.9659 | 0.9961 |
| 24 | 0.9679 | 0.9970 |

At small separation, 11 pixels discards 9--21% median response energy. The
31-pixel crop retains about 92--95% there and more than 99% from radius 16
outward. These fractions vary little with KL mode.

## Same-footprint covariance coverage

Training patches used half-overlap steps of 5, 15, and 23 pixels for the three
supports. Every query excluded the known planet and the union of all five
candidate footprints. A valid split required at least eight accepted patches
wholly within each detector half-plane.

| Radius | Narrowest valid 11-pixel band | Valid 31-pixel band | Valid 47-pixel band |
| ---: | ---: | ---: | ---: |
| 7.5 | 20 | None | None |
| 10 | 20 | None | None |
| 12 | 20 | None | None |
| 16 | 10 | None | None |
| 20 | 10 | None | None |
| 24 | 5 | None | None |

The 11-pixel geometry is well supported. At radii 7.5, 10, and 12, the
20-pixel band has at least 61, 75, and 93 total patches at every query. The two
half-plane minima are 21/24, 24/27, and 32/33.

The larger direct covariances are infeasible under the declared exclusion:

- 31-pixel training has only 0--16 accepted patches even over the full radial
  range. Some queries have none, and every radius has zero patches in at least
  one held-out half.
- 47-pixel training has exactly zero accepted patches for every selected query
  at every radius, even over the complete `6 <= r < 60` center range.

This is a geometric failure, not merely a poorly conditioned high-dimensional
fit. The response footprint remains valid for candidate detection, but it
cannot also define the covariance-estimation patch in this 128-pixel image.

## Decision

Direct empirical PCA and finite covariance estimation remain confined to the
11-pixel support. Larger response templates require a covariance model whose
estimation footprint is independent of its application footprint.

The next
+[decoupled-footprint preflight](../klip-stage-b-decoupled-preflight-setup-20260923/README.md)
+keeps 11-pixel Welch patches and five-pixel sampling while expanding only the
+candidate/source exclusion to 11, 31, or 47 pixels. If that geometry passes,
+the resulting PSD can be evaluated on the larger template grid while preserving
+the complete measured response.

## ROC receipt

The canonical result is
+`working/roc/klip_stage_b_footprint_preflight_20260923`. Its receipt records:

- `results.json`:
+  `8e7786536304afeba58b56028f1d11ff2c96ea0420ccb975904a8ca71665cbe5`;
+- `results.md`:
+  `ca9539a81f271fda07e12368c7a35ebbd8a6637c4529f4a8d73a0ae3973d006b`.

The result contains 259,676 bytes of detailed per-support, per-radius, and
per-band count distributions.
