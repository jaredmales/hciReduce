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
the individual raw-only and normalized-only decisions. For every radius,
brightness, and method, the report also gives the arithmetic mean of the valid
five-pixel maximum search SNRs and the number of valid searches. Each
per-injection SNR is read directly from the frozen production `hciAnalyze`
output map, with its annular mean subtraction, interpolated annular standard
deviation, known-source exclusion, and small-sample correction. The independent
Python annular calculation only verifies that output. Invalid searches are
omitted from the mean rather than assigned zero. The JSON output additionally
preserves the mean center-pixel SNR, allowing a separate check of small
search-position shifts.

## Pre-calibration runner repair

The first launch stopped after retaining 119 completed baseline analyses. The
remaining 45 parent analysis receipts came from the earliest phase of the
inner-radius study and predate its `active_methods` metadata field. Their
amplitude maps, SNR maps, scores, validity flags, and completion receipts are
complete; only that later convenience field is absent. No threshold was
written and no positive analysis started.

A read-only audit of all 272 parent analyses derives method support directly
from finite values at the five search pixels in each frozen parent
`hciAnalyze` SNR map. The derived sets reproduce all 272 recorded model-validity
sets exactly and reproduce `active_methods` for all 227 receipts that contain
it. The runner now uses this map-derived check for every receipt. The full raw
SNR-map and result comparisons remain in place, so this changes compatibility
handling rather than a numerical estimator.

After pulling the repair commit on ROC, update only the frozen runner record:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/p4-step5-matplotlib \
  /opt/conda/envs/xpy3_13/bin/python3 \
  agents/plans/scripts/compare_p4_step5_inner_patch_rms.py repair \
  --root working/roc/p4_inner_patch_rms_20260919
```

The repair verifies every unchanged frozen input, requires the exact
`'active_methods'` failure, refuses to run after calibration or any positive
analysis, records old and new runner fingerprints in `repair_0001.json`, and
retains the 119 completed baseline receipts. Restart the same `tmux` run command
with a new log such as `driver-repaired.log`. Each of the 45 unreceipted task
directories is preserved under `interrupted/baseline/` before recomputation.

The resumed launch retained 160 baseline receipts and isolated a second legacy
difference in four radius-12 calibration searches: `null_x116_y125`,
`null_x116_y127`, `null_x116_y129`, and `null_x116_y131`. Their rectangular ±5
amplitude/oracle searches are incomplete, while the earliest production run
recorded finite SNR output after its nonfinite-pixel handling. Their ±10, ±20,
identity, and Gaussian support is complete. None of the four is in a rectangular
±5 calibration pool or an evaluation site; they enter radius-12 calibration
pools only for methods whose support is complete.

Raw/reference planes must replay the frozen parent production support so the
complete raw SNR maps remain exact controls. New normalized planes retain the
stricter requirement for five finite amplitudes and five finite oracle values.
Trial scoring independently requires all five amplitude and SNR pixels, so the
four rectangular ±5 searches remain invalid rather than treating production
zeros as measurements. This separation preserves the raw map replay without
weakening the normalized estimator or changing any rectangular ±5 calibration
location.

After pulling the second repair commit, run the same `repair` command. It
requires the exact rectangular ±5 validity failure, verifies the prior repair
and every unchanged input, writes `repair_0002.json`, and retains all 160
completed baseline analyses. Restart with a new log such as
`driver-repaired-2.log`; only the four unreceipted tasks are recomputed.

## Prelaunch validation

Local synthetic checks verify equal-RMS equivalence with the raw rectangular
estimator, fitted-mean and variance-scale preservation, unit response, rejection
of numerically constant training patches, paired validity-only pool replacement,
and strict JSON handling of unavailable search pixels. A mock completed parent
with 164 baseline and 108 positive analysis receipts passes preparation and
input-fingerprint verification. A synthetic full summary reproduces all raw
parent controls and writes the expected 48 per-radius groups, eight aggregate
rows, and 54 paired comparisons. It also verifies the per-radius and aggregate
means of supplied production-map search and center SNR values. A synthetic
two-stage failure-root check verifies both exact failure guards, cumulative
repair receipts, successive runner-fingerprint replacement, retained-baseline
counts, and the complete post-repair fingerprint set. Python syntax compilation
and repository whitespace checks also pass.

## Completed ROC result

The comparison completed on ROC at 2026-09-19 17:23:44 MDT. It contains all
164 baseline analyses and 108 saved-positive analyses, with no incomplete task
directories and no new P4 reductions. Calibration completed before the first
positive output. The [remote review](remote_review.json) verifies 3,415 distinct
paths, both repair receipts, all per-task products, calibration products, final
products, and unchanged frozen inputs. Raw/reference SNR maps reproduce the
parent exactly; the largest independent-oracle difference is
\(5.96\times10^{-8}\).

### Calibrated recovery

Patch RMS changes several individual decisions but gives no consistent
aggregate recovery improvement:

| Method | Valid sites | 0.5× | 0.75× | 1× | Held-out nulls |
| --- | ---: | ---: | ---: | ---: | ---: |
| Raw rectangular ±5 | 22/36 | 4 | 8 | 11 | 0/36 |
| Patch-RMS rectangular ±5 | 22/36 | 3 | 7 | 11 | 0/36 |
| Raw rectangular ±10 | 27/36 | 6 | 10 | 16 | 0/36 |
| Patch-RMS rectangular ±10 | 27/36 | 6 | 11 | 16 | 0/36 |
| Raw rectangular ±20 | 36/36 | 7 | 11 | 18 | 1/36 |
| Patch-RMS rectangular ±20 | 36/36 | 6 | 11 | 18 | 0/36 |
| Identity | 36/36 | 8 | 13 | 22 | 0/36 |
| Gaussian FWHM 3.6 | 36/36 | 9 | 16 | 20 | 0/36 |

Across all widths and levels, normalized-only versus raw-only decisions are
2/4 for ±5, 2/1 for ±10, and 3/4 for ±20. Thus normalization gains seven
individual detections and loses nine. Removing the raw ±20 held-out exceedance
is favorable, but one event in this correlated field is not enough to establish
a false-positive improvement.

### Mean production `hciAnalyze` SNR

The entries below are raw → patch-RMS arithmetic means of the valid five-pixel
maximum SNRs. Each SNR comes directly from the production `hciAnalyze` map.
Counts are the valid sites at that radius; unsupported searches are omitted.

| Radius | Width | Valid | 0.5× mean SNR | 0.75× mean SNR | 1× mean SNR |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 6 | ±5 | 0/6 | — | — | — |
| 6 | ±10 | 0/6 | — | — | — |
| 6 | ±20 | 6/6 | 1.3015 → 1.2982 | 1.4724 → 1.4704 | 1.6152 → 1.6131 |
| 8 | ±5 | 0/6 | — | — | — |
| 8 | ±10 | 3/6 | 1.5300 → 1.5324 | 1.6290 → 1.6458 | 1.6445 → 1.6575 |
| 8 | ±20 | 6/6 | 1.3943 → 1.4040 | 1.6984 → 1.7054 | 1.8839 → 1.8884 |
| 12 | ±5 | 4/6 | 1.5943 → 1.6374 | 1.7041 → 1.7524 | 1.7889 → 1.8214 |
| 12 | ±10 | 6/6 | 1.8649 → 1.8872 | 2.1019 → 2.1101 | 2.2934 → 2.2885 |
| 12 | ±20 | 6/6 | 1.8595 → 1.8557 | 2.0901 → 2.0764 | 2.2767 → 2.2547 |
| 16 | ±5 | 6/6 | 1.2081 → 1.2511 | 1.6925 → 1.7239 | 2.1351 → 2.1686 |
| 16 | ±10 | 6/6 | 1.2561 → 1.3124 | 1.8311 → 1.8379 | 2.3216 → 2.2986 |
| 16 | ±20 | 6/6 | 1.3013 → 1.3422 | 1.8363 → 1.8486 | 2.3295 → 2.3231 |
| 20 | ±5 | 6/6 | 1.8824 → 1.8946 | 2.3848 → 2.3936 | 2.8283 → 2.8340 |
| 20 | ±10 | 6/6 | 1.7756 → 1.8546 | 2.2814 → 2.3597 | 2.7282 → 2.8059 |
| 20 | ±20 | 6/6 | 1.6914 → 1.9248 | 2.2011 → 2.4350 | 2.6670 → 2.8822 |
| 24 | ±5 | 6/6 | 1.6314 → 1.6078 | 2.1702 → 2.1535 | 2.6982 → 2.6842 |
| 24 | ±10 | 6/6 | 1.5563 → 1.5770 | 2.1047 → 2.1213 | 2.6475 → 2.6548 |
| 24 | ±20 | 6/6 | 1.4907 → 1.5906 | 2.0372 → 2.1331 | 2.5767 → 2.6500 |

Reference mean SNRs are:

| Radius | Identity, 0.5× / 0.75× / 1× | Gaussian, 0.5× / 0.75× / 1× |
| ---: | --- | --- |
| 6 | 1.3024 / 1.4404 / 1.5441 | 0.9973 / 1.1631 / 1.3054 |
| 8 | 1.4276 / 1.7439 / 1.9444 | 1.7278 / 2.0231 / 2.1889 |
| 12 | 1.8019 / 2.0157 / 2.1918 | 1.7490 / 2.0152 / 2.2393 |
| 16 | 1.2964 / 1.7756 / 2.2249 | 1.2142 / 1.7000 / 2.1252 |
| 20 | 1.8245 / 2.3361 / 2.7835 | 1.6601 / 2.1920 / 2.6624 |
| 24 | 1.4619 / 2.0005 / 2.5033 | 1.0438 / 1.5751 / 2.0835 |

Aggregated over valid sites, the normalized-minus-raw search-SNR changes are
+0.0164/+0.0152/+0.0128 for ±5, +0.0399/+0.0263/+0.0141 for ±10, and
+0.0628/+0.0556/+0.0438 for ±20. The corresponding center-pixel changes are
−0.0180/−0.0170/−0.0159, −0.0026/+0.0006/+0.0050, and
+0.0241/+0.0291/+0.0335. Patch normalization therefore raises the five-pixel
search maximum more consistently than the fixed-center SNR. The ±20 center
does improve modestly, so its search-SNR gain is not only a one-pixel peak
change.

The largest local improvement is radius-20 ±20: search SNR increases by about
0.23 at all three brightnesses, and recovery changes from 0/1/3 to 1/2/3.
Radius-24 ±20 gains 0.10/0.096/0.073 in mean SNR without changing recovery.
At radius 12, normalized thresholds rise by 0.106, 0.086, and 0.120 for
±5/±10/±20; this offsets the small SNR changes and reduces recovery, most
clearly for ±20 from 4/5/6 to 2/3/6.

### Conclusion

Post-mean patch RMS has a real but modest effect on the filter statistic and
can improve wide-pool SNR locally, especially at radii 20–24. It does not
produce a stable calibrated recovery gain across radii. Identity and Gaussian
remain stronger aggregate references, and Gaussian remains markedly strongest
at radius 8. These correlated, already inspected images do not support making
patch RMS or pooled rectangular PSD a production default.

## Scope

This remains a development comparison on the same inspected, correlated
residual field. It tests whether the earlier null result for post-mean patch
RMS changes with wider pooling and inner-radius coverage. It is not an
independent false-positive or threshold-transfer validation.
