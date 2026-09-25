# KLIP Stage-C immutable development preparer

## Purpose

Stage B selected two 11-pixel covariance arms: strict-radial
Hann/mixing-0.1 as the leading data-selected model and raw
rectangular/mixing-0.3 as the mandatory P4 prior. This preparer freezes the
development and validation experiment before any positive reduction is made.

The maintained runner is
[`prepare_klip_stage_c_development.py`](../../scripts/prepare_klip_stage_c_development.py).
It verifies the completed 47-pixel exact-response campaign and its Stage-A
lineage, including every sparse-response product.

## Frozen geometry

At each nominal radius 7.5, 10, 12, 16, 20, and 24 pixels, preparation starts
from native centers satisfying all of these conditions:

- complete 47-pixel exact responses in all eight KL modes at the center and
  four axial search pixels;
- finite 11-pixel final-image data in all eight modes at those five pixels;
- strict leave-site-out 3.6-pixel radial-profile support;
- at least eight accepted 11-pixel Welch patches in each detector half for
  every search pixel; and
- search centers outside the fixed seven-pixel known-planet disk.

For each radius, the preparer chooses the smallest predeclared training
half-width in 5, 10, 20, 40, and 60 pixels that leaves at least 38 common
candidates. It then assigns disjoint locations by deterministic angular
maximin selection, in this order:

| Role | Locations per radius | Total |
| :--- | ---: | ---: |
| Maximum-null calibration | 20 | 120 |
| Development injection | 6 | 36 |
| Validation injection | 6 | 36 |
| Held-out null | 6 | 36 |

This selection uses geometry and finite-support masks only. Baseline scores,
positive images, and recovery do not enter it.

## Source levels

The preparer computes site-specific physical contrasts for both development
and validation sites. At mode 200, it smooths the exact 47-pixel unit response
with Gaussian FWHM 3.6, takes the maximum source-only response over the fixed
five-pixel search, and divides by the interpolated one-pixel annular deviation
of the Gaussian-smoothed signal-free image. The annular profile excludes both
the known planet and the trial site with the production `R + 0.5` boundary.
The production small-sample factor then sets nominal source-only SNR 3, 5, and
7 contrasts.

This creates 108 development reduction commands and separately freezes 108
future validation commands. Validation positives remain unopened.

## Frozen method suite

The method manifest contains:

- native, Gaussian FWHM 3.6, exact identity, sparse identity, exact-response
  LPF 1.8 and 2.7, and exact fitted-mean isotropic permanent references;
- the development-only Gaussian widths 2.4, 3.0, and 4.2;
- raw 11-pixel rectangular/mixing-0.3;
- strict-radial 11-pixel Hann/mixing-0.1; and
- clipped and hard-truncated Hann precision at 0.5, 0.75, and 1.0 times mean
  covariance variance, plus the full inverse.

Positive covariance analysis will compare per-image refits with
baseline-frozen weights. Every method receives a separate maximum-of-20 null
threshold at each radius. Calibration must be receipted before a positive
image is analyzed.

## ROC preparation

```bash
git pull

root=working/roc/klip_stage_c_development_20260925

taskset -c 12-27 python3 \
  agents/plans/scripts/prepare_klip_stage_c_development.py check

taskset -c 12-27 python3 \
  agents/plans/scripts/prepare_klip_stage_c_development.py prepare "$root"
```

Preparation writes `protocol.json`, `geometry.json`, `contrasts.json`, all 108
development commands, copied software/configuration, and a strict manifest.
The development reduction and analysis command will be added only after the
prepared geometry and contrast receipt pass review.

The first two verification attempts exposed packaging-only imports before any
reduction directory was created. The preparer initially omitted the KLIP PSD
extension module, and that module eagerly imported a P4 comparison chain used
only by its standalone self-check. The corrected package copies the extension
module and defers the P4 import to that self-check. Both untouched preparation
directories are archived and replaced rather than repaired in place.
