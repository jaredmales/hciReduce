# Step 5: measured-response smoothing setup

## Fixed question

Does removing fine spatial structure from the independently measured P4
response improve injection SNR under either identity weighting or the current
small-separation rectangular-PSD weighting?

This test follows the observation that production Gaussian FWHM-3.6 smoothing
has slightly higher mean injection SNR than the identity matched filter. If the
stored measured response were exact and the image noise were white, smoothing
that response could only reduce expected matched-filter SNR. An observed gain
would therefore show that the removed structure is not useful in these
injections. It could arise from response-estimation error, response variation,
or interaction with non-white residual noise.

## Fixed comparison

The response-template Gaussian low-pass FWHM grid is **0, 0.9, 1.8, 2.7, and
3.6 pixels**, or 0, 0.25, 0.5, 0.75, and 1.0 lambda/D. Every width is applied
to both:

- identity weights, using the smoothed response divided by its squared norm;
- raw rectangular PSD weights with ±20-pixel training-center pooling and the
  existing 0.3 flat-spectrum mixture.

The ±20 branch is used because it has complete support at all 36 inner sites.
The production Gaussian FWHM-3.6 image remains a reference. Its width describes
direct image smoothing; it is distinct from convolving an already processed
response with a 3.6-pixel kernel.

Only each native 11-by-11 response template changes. The implementation uses a
constant-zero exterior, a four-sigma Gaussian kernel cutoff, and no template
renormalization. The matched-filter amplitude normalization makes the detection
statistic invariant to a scalar rescaling of a template. The report separately
records paired amplitude throughput.

The completed SNR-3/5/7 reductions and the completed one-lambda/D-masked
analysis are immutable inputs. The test retains the same injections, known
planet and current-trial exclusions, candidate pixels, covariance samples,
fitted PSD mean, five-pixel search, and production `hciAnalyze` annular SNR.
It performs no new P4 reduction.

Within each weighting family, all smoothing widths reuse the same 20 frozen
calibration locations selected by the completed masked comparison. Each method
gets its own maximum-null threshold, frozen before positive images are read.
The zero-width identity, zero-width rectangular ±20, and production Gaussian
maps, SNRs, thresholds, and decisions must reproduce the completed parent.

## Interpretation diagnostics

For every smoothed matched filter, the output includes:

- mean five-pixel maximum and fixed-center production SNR;
- recovery and held-out null counts;
- positive-minus-baseline amplitude throughput divided by injected contrast;
- the predicted matched-SNR efficiency if the unsmoothed measured response
  were the exact source template.

For identity weighting, the last quantity is the cosine between the original
and smoothed templates. For PSD weighting it is the corresponding covariance
inner-product cosine. A measured SNR improvement despite a predicted
exact-template loss is the key signature that the removed response structure
does not reproduce as useful signal.

This first comparison does not truncate or otherwise change the PSD inverse.
That remains a separate follow-up after the template-only result.

## Validation and ROC commands

Run the local algebra checks after pulling the setup:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/p4-step5-matplotlib \
  /opt/conda/envs/xpy3_13/bin/python3 \
  agents/plans/scripts/compare_p4_step5_response_smoothing.py check
```

Prepare the immutable saved-image comparison on ROC:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/p4-step5-matplotlib \
  /opt/conda/envs/xpy3_13/bin/python3 \
  agents/plans/scripts/compare_p4_step5_response_smoothing.py prepare \
  --comparison working/roc/p4_inner_snr357_patch_rms_lambdad_20260919 \
  --root working/roc/p4_response_smoothing_20260920 \
  --cpus 0 1 2 3 4 5 6 7 8 9 10 11
```

Launch the resumable analysis:

```sh
tmux new-session -d -s p4-response-smoothing \
  "cd /home/jrmales/Source/mxApps/hciReduce && \
   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/p4-step5-matplotlib \
   /opt/conda/envs/xpy3_13/bin/python3 \
   agents/plans/scripts/compare_p4_step5_response_smoothing.py run \
   --root working/roc/p4_response_smoothing_20260920 \
   > working/roc/p4_response_smoothing_20260920/driver.log 2>&1"
```

Monitor it with:

```sh
cat working/roc/p4_response_smoothing_20260920/state.json
tail -n 30 working/roc/p4_response_smoothing_20260920/driver.log
```

The final products are `results.json`, `results.md`, and `comparison.png` under
the new root. This is a paired development reuse of one correlated field, not
an independent validation set.
