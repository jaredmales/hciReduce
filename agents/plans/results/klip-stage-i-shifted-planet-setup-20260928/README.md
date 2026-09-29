# KLIP Stage-I shifted-template planet setup

## Purpose

Stage H showed that cubic registration recovers most of the KLIP response
mismatch at the planet's fractional pixel phase. Stage I applies that model to
the real planet and the twelve matching Stage-G injections and asks whether it
resolves the fixed-pixel discrepancy.

The planet is evaluated at the independently optimized negative-companion
position, not the earlier configured analysis coordinate. Its nearest native
pixel and fractional phase therefore match the injection construction:
`(-0.27707,+0.48596)` pixels. This removes the approximately 0.67-pixel
position difference that was present when the Stage-F configured coordinate
was compared with optimized-phase injections.

## Frozen design

No new KLIP reductions are required. The runner reuses:

- the original planet science cube;
- the twelve completed nominal optimized-phase Stage-G injections;
- the mode-200 radius-12 signal-free calibration unit;
- the promoted 47-pixel exact-response field; and
- the completed Stage-F through Stage-H products.

Mode 200 and the nearest native pixel to the optimized source position are the
predeclared endpoint. The completed Stage-G aperture-maximum result remains
unchanged and is not retested.

Eight maps are measured:

| Family | Methods |
| :--- | :--- |
| Gaussian controls | FWHM 2.4 and 3.6 pixels |
| Identity response | Current integer and cubic shifted |
| Raw rectangular covariance | Current integer and cubic shifted |
| Radial-Hann covariance | Current integer and cubic shifted, truncation 0.75 |

The primary paired comparisons are shifted covariance minus shifted identity
and Gaussian 3.6 minus shifted identity. Shifted-minus-integer differences
measure the registration correction itself. Gaussian 3.6 minus Gaussian 2.4
is retained as the earlier fixed-pixel control.

## Annular-noise contract

Stage I does not change only the source pixel and reuse an integer-template
noise denominator. It reconstructs cubic-shifted identity and covariance
weights at all 245 radius-12 calibration positions in mode 200. The stored
integer weights must be reproduced before the shifted weights are accepted.

For each science image, the shifted weights are applied throughout the two
one-pixel annuli used to normalize the optimized nearest pixel. The source
pixel is then replaced with a candidate-specific shifted weight fitted from
the signal-free baseline with the complete radius-three source aperture
excluded from covariance training. The known planet and, for injections, the
trial source are excluded from annular noise.

Production `hciAnalyze` performs mean-centered, small-sample-corrected SNR
normalization. An independent oracle must reproduce every finite SNR pixel.
The unchanged Gaussian and integer-template injection SNRs must also reproduce
their completed Stage-G nearest-pixel values.

## Interpretation

The primary question is whether the planet's shifted raw- and radial-covariance
penalties relative to shifted identity fall inside the twelve-injection
distributions. The Gaussian-3.6 comparison is evaluated the same way.
Exchangeable ranks and the individual injection range accompany diagnostic
standard-deviation comparisons.

Because every task has the same phase and source statistic, a remaining planet
discrepancy cannot be attributed to the integer-grid response centering found
in Stage H. It would instead point to residual PSF morphology, local noise, or
another planet-specific effect.

## ROC commands

After pulling the setup commit, run from the repository root:

```bash
root=working/roc/klip_stage_c_development_20260925
python3 agents/plans/scripts/run_klip_stage_i_shifted_planet.py check
taskset -c 12-27 python3 \
  agents/plans/scripts/run_klip_stage_i_shifted_planet.py prepare "$root"
```

Then launch the frozen command printed by `prepare`:

```bash
root=working/roc/klip_stage_c_development_20260925
taskset -c 12-27 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 python3 \
  "$root/stage_i_shifted_planet/software/run_klip_stage_i_shifted_planet.py" \
  run "$root" \
  > "$root/stage_i_shifted_planet/driver.log" 2>&1
```

The run is resumable and can be left to completion. Monitor it with:

```bash
tail -f working/roc/klip_stage_c_development_20260925/stage_i_shifted_planet/driver.log
```

The primary report will be
`working/roc/klip_stage_c_development_20260925/stage_i_shifted_planet/results.md`.
