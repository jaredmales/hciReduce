# KLIP 47-pixel exact-response campaign setup

## Purpose

The stamp-convergence and contrast-linearity experiments show that KLIP's
small-separation response is nonlocal and that 47 pixels is the first common
footprint satisfying the response-edge criteria. This campaign regenerates the
native exact response field with that promoted footprint before covariance
screening begins.

## Frozen calculation

The campaign uses the complete response-tail experiment, response-convergence
experiment, and Stage-A audit as a verified lineage. It preserves:

- the archived `klipReduce` executable and all 621 input frames;
- the bitwise-reproduced signal-free baseline;
- the fitted planet subtraction at separation 12.387750470790344 pixels,
  position angle 260.64315155951886 degrees, and contrast
  0.0045743624964148452;
- the same contrast as the positive/negative response half-amplitude;
- KL modes 125, 150, 175, 200, 225, 250, 300, and 350; and
- the complete `6 <= r < 60` search annulus.

One native `refitDifference` calculation measures all 11,192 integer search
locations, or 22,384 signed KLIP trials. It publishes a schema-2 `PIXEL_EXACT`
product with 47-by-47 response and validity stamps. The response accumulator is
expected to retain 988,925,120 bytes. The response and validity FITS products
will occupy about 1.6 GB before logs and receipts.

The earlier 11-pixel campaign took 16.63 hours while averaging the equivalent
of 42.45 CPU cores. The current reproducibility contract uses CPUs 12--27 and
recent eight-mode reductions take about 5 seconds each. Allow roughly 24--32
hours for this run; the single-process implementation can be faster because it
loads the 621 inputs only once.

The native response loop is not resumable within a partially completed run. A
failed response directory is preserved for diagnosis and must be moved aside
before a fresh prepared campaign is started. A completed native product is
recoverable if the wrapper stops during validation or hashing.

## Acceptance checks

After the native calculation, the runner requires:

1. schema-2 `PIXEL_EXACT`, `refitDifference`, and
   `PAIRED_FINAL_DIFFERENCE` metadata;
2. exactly 11,192 response locations, 22,384 signed trials, eight fixed modes,
   47-pixel stamps, and 988,925,120 retained bytes;
3. coordinates identical to the archived 11-pixel response field;
4. binary validity and finite supported values with a valid anchor everywhere;
5. a final signal-free image bitwise identical to the Stage-A baseline;
6. minimum cosine 0.999 and projection within 0.01 for the central 11 pixels at
   every archived location and mode; and
7. the same replay thresholds for the full 47-pixel response at the 18
   independent contrast-linearity sites.

The central replay detects any change caused by enlarging the extraction. The
external full-stamp replay verifies the newly stored tail against separately
executed positive/negative reductions.

## ROC commands

Run from the repository root in a durable ROC shell or tmux session:

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

The long command writes progress to `response/run.log`. Monitor it from another
shell with:

```bash
tail -f "$root/response/run.log"
```

On success the runner hashes all 19 response products, performs both replay
checks, and writes `results.json`, `results.md`, and `complete.json`.
