# KLIP covariance Stage-A result

## Scope

This checkpoint implements and exercises the response and geometry audit from
the [KLIP covariance test program](../klip-covariance-matched-filter-program-20260920/README.md).
The maintained runner is
[`run_klip_covariance_stage_a.py`](../../scripts/run_klip_covariance_stage_a.py).
It fingerprints all 621 raw inputs, both response product sets, the science
cubes, configuration, PSF, and executables before calculating any new score.

The canonical ROC run used CPUs 12--27. Its `klipReduce` executable is the
archived response-run binary; the newer `hciAnalyze` executable independently
reproduces the archived analysis results:

| Executable | Archived SHA-256 | ROC SHA-256 | Match |
| --- | --- | --- | --- |
| `klipReduce` | `6abd73a8a9e9e64df979923da3e16f2275e93f089df181e2b7e9dd042ec3a6e5` | `6abd73a8a9e9e64df979923da3e16f2275e93f089df181e2b7e9dd042ec3a6e5` | Yes |
| `hciAnalyze` | `f077812bf9390e05629d2af6b1993b7395fb1f48dd9766fb57600b3c6136d747` | `fb25e3123dd6f59c190ce7151a04e9d030368efeda9fb4d68adbfe043609dfec` | No |

## Audit results

| Check | Result | Diagnostic |
| --- | --- | --- |
| Exact product schema, modes, coordinates, and normalization | Pass | 11,192 positive-energy native templates in every mode |
| Five-pixel common-support geometry | Pass | Counts reproduce 42, 44, 60, 91, 113, and 160 at radii 7.5, 10, 12, 16, 20, and 24 pixels |
| Independent exact identity-filter reconstruction | Pass | Maximum difference from direct `hciAnalyze` annular SNR: $4.77\times10^{-7}$ |
| Archived planet controls | Pass | All 32 SNR values reproduce at the printed precision |
| Current versus archived signal-free baseline | Pass | Bitwise identical; maximum and RMS differences are zero |
| Preregistered 11-by-11 border-energy trigger | Triggered | Every primary mode-200 radial bin exceeds the one-percent median threshold |

The baseline compatibility gate passes exactly: the new signal-free cube and
the archived exact-response baseline have identical finite masks and identical
pixel values. The archived response may therefore be used with new ROC
reductions made by this frozen `klipReduce` binary. The different `hciAnalyze`
hash is acceptable for this stage because its independent identity-filter
replay agrees within $4.77\times10^{-7}$ in SNR and all 32 archived planet
values reproduce at their printed precision.

## Response footprint

The mode-200 response has substantial negative-lobe energy and measurable
energy at the edge of the 11-by-11 stamp:

| Radius (pixels) | Median negative energy fraction | Median border energy fraction |
| ---: | ---: | ---: |
| 6 | 0.4602 | 0.0631 |
| 7.5 | 0.4088 | 0.0363 |
| 10 | 0.3437 | 0.0135 |
| 12 | 0.3112 | 0.0154 |
| 16 | 0.2186 | 0.0150 |
| 20 | 0.1866 | 0.0157 |
| 24 | 0.1662 | 0.0185 |

The border result activates the preregistered larger-stamp response experiment
before an optimal matched-filter claim. At radius 7.5, one template also
exceeds the separate five-percent selected-location threshold, so later site
selection cannot avoid the larger-stamp question by relying on the median
criterion alone.

## Exact versus sparse response

The sparse radial response is close to the exact response outside the inner
boundary, but the mismatch grows at small separation:

| Radius (pixels) | Sampled exact locations | Median cosine | Median sparse-to-exact projection | Median relative residual |
| ---: | ---: | ---: | ---: | ---: |
| 6 | 12 | 0.8535 | 0.4155 | 0.5212 |
| 7.5 | 52 | 0.9752 | 1.0372 | 0.2215 |
| 10 | 56 | 0.9479 | 1.0486 | 0.3185 |
| 12 | 64 | 0.9864 | 0.9965 | 0.1647 |
| 16 | 64 | 0.9914 | 1.0922 | 0.1306 |
| 20 | 64 | 0.9958 | 1.0229 | 0.0921 |
| 24 | 64 | 0.9967 | 1.0217 | 0.0817 |

These are mode-200 values on common valid support. Sampling is deterministic:
up to 64 complete exact locations per radial bin, evenly indexed after sorting
by position angle.

## ROC receipt

The complete run is
`working/roc/klip_covariance_stage_a_20260921`. Its completion receipt records:

- `results.json`: `7c2708e6758bc42e156229b7765c1f50e5fe1a2cf5c30d659873120b0f4dac52`;
- `results.md`: `721568bfda8cde1559675b9231b80a2577693b2af9c8ee8bc5d9c50f50c40058`.

The receipt status is complete and lists the baseline, response audit, sparse
response audit, identity replay, planet controls, and report as completed. The
only activated follow-up is the larger-stamp response experiment.
That follow-up is frozen in the
[response-stamp convergence setup](../klip-response-stamp-convergence-setup-20260921/README.md).
