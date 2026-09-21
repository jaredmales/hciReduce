# KLIP covariance Stage-A implementation checkpoint

## Scope

This checkpoint implements and exercises the response and geometry audit from
the [KLIP covariance test program](../klip-covariance-matched-filter-program-20260920/README.md).
The maintained runner is
[`run_klip_covariance_stage_a.py`](../../scripts/run_klip_covariance_stage_a.py).
It fingerprints all 621 raw inputs, both response product sets, the science
cubes, configuration, PSF, and executables before calculating any new score.

The local end-to-end validation used 16 pinned CPUs and current build-tree
executables. It is an implementation check rather than the canonical ROC
receipt because both executable hashes differ from the archived response run:

| Executable | Archived SHA-256 | Validation SHA-256 |
| --- | --- | --- |
| `klipReduce` | `6abd73a8a9e9e64df979923da3e16f2275e93f089df181e2b7e9dd042ec3a6e5` | `5b50bde7af147f3eda9f6be349c3e381487c17ce54898dc06335d2aed040b695` |
| `hciAnalyze` | `f077812bf9390e05629d2af6b1993b7395fb1f48dd9766fb57600b3c6136d747` | `3ab9a69c5798512f9a35f7c96dda0c503d7e4274e819df11eb811779ed1d96e6` |

## Audit results

| Check | Result | Diagnostic |
| --- | --- | --- |
| Exact product schema, modes, coordinates, and normalization | Pass | 11,192 positive-energy native templates in every mode |
| Five-pixel common-support geometry | Pass | Counts reproduce 42, 44, 60, 91, 113, and 160 at radii 7.5, 10, 12, 16, 20, and 24 pixels |
| Independent exact identity-filter reconstruction | Pass | Maximum difference from direct `hciAnalyze` annular SNR: $4.77\times10^{-7}$ |
| Archived planet controls | Pass | All 32 SNR values reproduce at the printed precision |
| Current versus archived signal-free baseline | Fail | Same finite mask; maximum absolute difference 0.03318 and RMS difference $1.23\times10^{-4}$ |
| Preregistered 11-by-11 border-energy trigger | Triggered | Every primary mode-200 radial bin exceeds the one-percent median threshold |

The baseline failure is a deliberate blocking gate. The validation build has
broad floating-point differences at the strict `rtol=2e-6`, `atol=5e-7`
tolerance and isolated larger differences, led by 0.03318 at mode 125 and
0.01846 at mode 200. Injection reductions cannot be mixed with the archived
response until the canonical ROC run either reproduces the archived baseline
or freezes a compatible historical reduction binary. Regenerating the exact
response with the selected current build is the remaining fallback.

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

## ROC command

Run preparation and analysis under the same affinity. The copied runner in the
experiment directory is the frozen analysis entry point:

```bash
root=working/roc/klip_covariance_stage_a_20260921
taskset -c 12-27 python3 agents/plans/scripts/run_klip_covariance_stage_a.py check
taskset -c 12-27 python3 agents/plans/scripts/run_klip_covariance_stage_a.py prepare "$root"
taskset -c 12-27 python3 "$root/software/run_klip_covariance_stage_a.py" run "$root"
```

The runner writes strict `protocol.json`, `manifest.json`, `state.json`,
`results.json`, `results.md`, and completion receipts. It keeps a failed task
directory intact and verifies completed products before reusing them.
