# KLIP Stage-F known-planet closure setup

## Purpose

This setup implements the last checkpoint in the frozen KLIP matched-filter
program. It measures the original planet-bearing KLIP cube only after the
Stage-E policy, held-out nulls, validation reductions, and validation result
are recursively verified. The planet remains a descriptive endpoint and cannot
select or rescue a method.

The runner is
[`run_klip_stage_f_planet.py`](../../scripts/run_klip_stage_f_planet.py).
It uses the canonical completed experiment at
`working/roc/klip_stage_c_development_20260925` and the analysis geometry in
`working/analyze.conf`.

## Methods

Stage F applies the seven immutable Stage-E methods:

- native image;
- Gaussian FWHM 2.4 and 3.6;
- exact- and sparse-response identity filtering;
- raw rectangular PSD with mixing 0.3; and
- radial-standardized Hann PSD with mixing 0.1 and hard truncation at 0.75
  mean variance.

It also reports the three permanent Stage-C controls that were fixed before
injection scores were read: exact-response low-pass filters at 1.8 and 2.7
pixels, and exact identity with the fitted PSD mean removed. These controls are
marked as absent from Stage E in the closure table; their planet values do not
supply validation evidence.

## Source-safe aperture calculation

The generic Stage-C maps mask the known-planet footprint and therefore cannot
be sampled at the planet. Stage F handles the two roles separately:

1. It applies the already frozen generic weights to the original science cube
   and stitches overlapping radial bands with the nearest predeclared Stage-C
   radius policy. Ties use the smaller radius. These maps supply the one-pixel
   annular noise statistics on common method support.
2. For every native pixel in the three-pixel analysis aperture, it rebuilds
   only the candidate weights from the signal-free baseline. Training excludes
   the optimized-planet disk and the union of every aperture pixel's complete
   11-by-11 footprint. Candidate support uses the exact or sparse response
   validity without the planet mask, so the source is retained in the
   amplitude estimate.

The `working/analyze.conf` center (separation 11.782 pixels, PA 262.051 degrees)
defines the three-pixel search aperture and seven-pixel SNR exclusion. The
optimized subtraction (separation 12.3877505 pixels, PA 260.643152 degrees)
continues to define the signal-free baseline and its training exclusion. The
runner records both centers.

Each aperture pixel uses the nearest frozen policy radius among 7.5, 10, 12,
16, 20, and 24 pixels. Thus the aperture can cross radial-policy boundaries
without choosing a policy from the observed planet. The stitched generic maps
provide continuous one-pixel annular coverage from the radius-7.5 through
radius-24 units.

## Frozen and verified behavior

`prepare` requires the response campaign's frozen CPU affinity (CPUs 12--27),
the unopened `validation_complete` state, and recursively verifies the immutable policy, 48 validation model units, 36 held-out analyses,
108 validation reductions, and 108 validation analyses. It follows the
promoted exact-response provenance to the original science cube and requires
that cube to match its Stage-A frozen inventory. It then copies this runner and
`working/analyze.conf` into `stage_f_planet/` and freezes all input, software,
and calibration-unit receipts.

`run` verifies that manifest again before opening the planet cube. It writes
nearest-pixel and aperture-maximum SNR, peak position, offset, amplitude,
response scale, contrast estimate, and support for every method and all eight
KL modes. The production `hciAnalyze` SNR maps must reproduce an independent
one-pixel annular oracle at `rtol=2e-6`, `atol=2e-6`.

A failed `run` is restartable. The next invocation moves the partial `planet/`
directory under `stage_f_planet/interrupted/` and replays from the immutable
manifest. A completed invocation verifies and reprints the existing report.

## ROC commands

From the repository root on ROC after pulling the commit that contains this setup:

```bash
root=working/roc/klip_stage_c_development_20260925
python3 agents/plans/scripts/run_klip_stage_f_planet.py check
taskset -c 12-27 python3 \
  agents/plans/scripts/run_klip_stage_f_planet.py prepare "$root" \
  --config working/analyze.conf
```

Run the frozen copy printed by `prepare`:

```bash
root=working/roc/klip_stage_c_development_20260925
taskset -c 12-27 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 python3 \
  "$root/stage_f_planet/software/run_klip_stage_f_planet.py" run "$root" \
  > "$root/stage_f_planet/driver.log" 2>&1
```

To monitor without `tmux`:

```bash
tail -f working/roc/klip_stage_c_development_20260925/stage_f_planet/driver.log
```

The primary report will be
`working/roc/klip_stage_c_development_20260925/stage_f_planet/planet/results.md`.
The same directory contains strict JSON and CSV tables, amplitude, response,
policy-radius, and SNR FITS products, the independent annular audit, and
candidate-fit diagnostics.

## Completed result

Stage F completed on ROC. The compact result, verification record, comparisons,
and scientific interpretation are preserved in the
[Stage-F known-planet closure result](../klip-stage-f-planet-20260926/README.md).
The large FITS products remain in the ROC run directory and are covered by the
preserved completion receipt.
