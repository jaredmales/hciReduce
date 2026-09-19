# Step 5: wider-pool post-mean patch-RMS comparison

## Fixed question

Does per-patch RMS normalization after ensemble mean-patch subtraction help when
the pooled rectangular PSD uses the wider radial ranges needed at small
separations?

This is a paired reanalysis of the completed inner-radius study. It compares
raw and post-mean patch-RMS rectangular PSD estimates at radial-center
half-widths **±5, ±10, and ±20 pixels**. The evaluation radii remain 6, 8, 12,
16, 20, and 24 pixels. Identity matched filtering and Gaussian smoothing with
FWHM 3.6 pixels remain reference methods.

The comparison reuses all 164 baseline images/searches and 108 saved positive
P4 images. It performs **zero new P4 reductions**. This keeps the reductions,
injection sites and brightnesses, known-planet guard, trial holdouts, response
templates, five-pixel searches, and production annular SNR fixed.

## Paired estimator

For each radial width and candidate pixel, extract the same raw aligned 11×11
training patches used by the completed study. Given patches \(x_j\), first form
the raw ensemble mean and residuals,

\[
\mu=\frac{1}{N}\sum_j x_j,\qquad r_j=x_j-\mu.
\]

The raw arm estimates its rectangular periodogram from \(r_j\). The normalized
arm computes

\[
s_j=\sqrt{\frac{1}{121}\sum_a r_{j,a}^2},\qquad z_j=\frac{r_j}{s_j},
\]

then averages the rectangular 21×21 zero-padded Fourier powers of \(z_j\)
without a second ensemble centering. Its spectrum is rescaled to the original
raw mean pixel variance before applying the same 0.3 flat-spectrum mixture.
The fitted mean \(\mu\), candidate stamp, and response template all remain in
raw contrast units. No spatial scalar mean is removed from an individual
patch, and the candidate is never normalized by its own RMS.

Any nonfinite or roundoff-scale patch RMS invalidates the normalized fit. The
driver does not discard a patch or refit the mean to avoid that condition.

## Calibration and validity

Each normalized method uses the same geometry-selected 20-search calibration
locations as its matching raw method. Validity-only replacement is allowed
from the already frozen candidate union, using the completed study's maximin
rule and no score values. The driver requires the resulting raw/normalized
locations to be identical. Each method then receives its own threshold: the
maximum of its 20 five-pixel baseline scores. Every threshold is frozen before
the saved positives are read.

All five search amplitudes and all five independently predicted annular-SNR
values must be finite. An incomplete search is a nondetection. The original
known-planet circle and the current trial's full search footprint remain
excluded from covariance training. Annular SNR retains the parent study's
known-planet exclusion and small-sample correction.

The raw ±5, ±10, and ±20 amplitude/SNR maps, identity and Gaussian maps,
effective calibration pools, thresholds, individual decisions, and summary
counts must reproduce the completed parent study. These checks make the new
patch weighting the only numerical change in the paired PSD arms.

## ROC execution

After pulling this commit on ROC, prepare the immutable comparison from the
repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/p4-step5-matplotlib \
  /opt/conda/envs/xpy3_13/bin/python3 \
  agents/plans/scripts/compare_p4_step5_inner_patch_rms.py prepare \
  --study working/roc/p4_inner_rectangular_20260919 \
  --root working/roc/p4_inner_patch_rms_20260919 \
  --cpus 0 1 2 3 4 5 6 7 8 9 10 11
```

Preparation requires the parent completion receipt, 272 completed parent
analysis receipts, and 108 completed positive reductions. It fingerprints
every parent product and script read by the comparison. Do not change those
files between `prepare` and completion.

Start the resumable analysis in `tmux`:

```sh
tmux new-session -d -s p4-inner-patch-rms \
  "cd /home/jrmales/Source/mxApps/hciReduce && \
   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/p4-step5-matplotlib \
   /opt/conda/envs/xpy3_13/bin/python3 \
   agents/plans/scripts/compare_p4_step5_inner_patch_rms.py run \
   --root working/roc/p4_inner_patch_rms_20260919 \
   > working/roc/p4_inner_patch_rms_20260919/driver.log 2>&1"
```

Monitor without changing the run:

```sh
cat working/roc/p4_inner_patch_rms_20260919/state.json
tail -n 30 working/roc/p4_inner_patch_rms_20260919/driver.log
```

The run first completes all 164 baseline reanalyses and writes a calibration
receipt. It then analyzes the 108 saved positive images. A task with a verified
completion receipt is reused; an interrupted unreceipted task is preserved
under `interrupted/` before recomputation. Final products are `results.json`,
`results.md`, and `comparison.png`, with aggregate and per-radius recovery plus
the individual raw-only and normalized-only decisions.

## Prelaunch validation

Local synthetic checks verify equal-RMS equivalence with the raw rectangular
estimator, fitted-mean and variance-scale preservation, unit response, rejection
of numerically constant training patches, paired validity-only pool replacement,
and strict JSON handling of unavailable search pixels. A mock completed parent
with 164 baseline and 108 positive analysis receipts passes preparation and
input-fingerprint verification. A synthetic full summary reproduces all raw
parent controls and writes the expected 48 per-radius groups, eight aggregate
rows, and 54 paired comparisons. Python syntax compilation and repository
whitespace checks also pass.

## Scope

This remains a development comparison on the same inspected, correlated
residual field. It tests whether the earlier null result for post-mean patch
RMS changes with wider pooling and inner-radius coverage. It is not an
independent false-positive or threshold-transfer validation.
