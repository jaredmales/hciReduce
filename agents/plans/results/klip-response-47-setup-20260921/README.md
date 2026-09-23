# KLIP 47-pixel exact-response campaign result

## Purpose

The stamp-convergence and contrast-linearity experiments showed that KLIP's
small-separation response is nonlocal and that 47 pixels is the first common
footprint satisfying the response-edge criteria. This campaign regenerated the
native exact response field with that promoted footprint before covariance
screening.

## Frozen calculation

The campaign preserved:

- the archived `klipReduce` executable and all 621 input frames;
- the bitwise-reproduced signal-free baseline;
- the fitted planet subtraction at separation 12.387750470790344 pixels,
  position angle 260.64315155951886 degrees, and contrast
  0.0045743624964148452;
- the same contrast as the positive/negative response half-amplitude;
- KL modes 125, 150, 175, 200, 225, 250, 300, and 350; and
- the complete `6 <= r < 60` search annulus.

One native `refitDifference` calculation measured all 11,192 integer search
locations, or 22,384 signed KLIP trials. It produced schema-2 `PIXEL_EXACT`
response and validity stamps with 47-by-47 support in every mode.

## Result

All preregistered gates pass.

| Check | Result | Diagnostic |
| --- | --- | --- |
| Signal-free baseline replay | Pass | Bitwise identical; maximum difference 0 |
| Search coordinates | Pass | Identical to all 11,192 archived locations |
| Central 11 pixels versus archive | Pass | Minimum cosine 1; maximum projection error 0 across every location and mode |
| Full 47 pixels versus independent derivatives | Pass | Minimum cosine 0.999932; maximum projection error $6.45\times10^{-5}$ |

The full-stamp comparison used the 18 independently reduced sites from the
contrast-linearity experiment. Per-mode worst cases were:

| Mode | Minimum cosine | Maximum projection error |
| ---: | ---: | ---: |
| 125 | 1.000000 | $1.31\times10^{-5}$ |
| 150 | 1.000000 | $1.51\times10^{-5}$ |
| 175 | 1.000000 | $1.86\times10^{-5}$ |
| 200 | 1.000000 | $1.38\times10^{-5}$ |
| 225 | 0.999932 | $6.45\times10^{-5}$ |
| 250 | 1.000000 | $1.75\times10^{-5}$ |
| 300 | 1.000000 | $2.02\times10^{-5}$ |
| 350 | 0.999999 | $2.26\times10^{-5}$ |

The native run took 107,517 seconds, or 29.87 hours, on CPUs 12--27. The
complete campaign occupies 1.6 GB. It finished without wrapper interruption or
product recovery.

## Interpretation

Increasing the stored footprint does not change the reduction or the response
core: every central 11-pixel stamp is exactly equal to its archived value. The
newly retained tail also agrees with external positive/negative reductions far
more closely than the promotion thresholds require. The 47-pixel field is
therefore the accepted exact-response oracle for covariance development.

The next
[Stage-B footprint preflight](../klip-stage-b-footprint-preflight-setup-20260923/README.md)
compares central 11-, 31-, and 47-pixel support. It measures response-energy
capture, common five-pixel search coverage, and half-overlap covariance-training
coverage before reading a baseline score or fitting a covariance.

## ROC receipt

The canonical campaign is
`working/roc/klip_response_47_20260921`. Its completion receipt records:

- `response/complete.json`:
  `6e532e909155ebd4f83435e5649c8f1ed047bb1d99d7821762e8d737252b5d95`;
- `results.json`:
  `0408bdf91866e589a932f58d8c60234ff4c9dae97b6bd3c82b8832ea5a6e52e2`;
- `results.md`:
  `ce692cb69cbeec1262fc8e019737d10d80b37ead7f2f8a230b86e82542e9bc76`.

The response receipt fingerprints the final image, manifest, coordinates, and
all 16 mode-specific response and validity products.

## Reproduction

```bash
git pull

root=working/roc/klip_response_47_20260921

taskset -c 12-27 python3 \
  agents/plans/scripts/run_klip_response_47.py check

taskset -c 12-27 python3 \
  agents/plans/scripts/run_klip_response_47.py prepare "$root"

taskset -c 12-27 python3 \
  "$root/software/run_klip_response_47.py" run "$root"
```

The native response loop is not resumable within a partially completed run. A
completed native product remains recoverable if only the validation wrapper
stops.
