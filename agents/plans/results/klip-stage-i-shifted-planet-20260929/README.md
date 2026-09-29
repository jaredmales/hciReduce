# KLIP Stage-I shifted-template planet result

## Result

Stage I completed the mode-200 reanalysis of the planet and twelve matched
optimized-phase injections. Cubic response registration raises the planet SNR
for identity and both covariance filters, and each planet shift lies inside
the corresponding injection distribution. The response-centering correction
therefore behaves normally for the planet.

Registration does not remove the covariance penalty relative to identity. The
planet's shifted raw-covariance penalty is `-0.6921` SNR and its shifted
radial-covariance penalty is `-0.7004` SNR. Both are below all twelve injection
values. They are respectively `-2.13` and `-2.00` injection sample standard
deviations from the injection means.

### Registration effect

| Shifted minus integer SNR | Injections | Planet | Planet deviation |
| :--- | ---: | ---: | ---: |
| Identity | +0.0671 +/- 0.2753 | +0.2680 | +0.73 SD |
| Raw rectangular covariance | +0.0580 +/- 0.2848 | +0.1293 | +0.25 SD |
| Radial-Hann covariance | +0.0948 +/- 0.2869 | +0.1654 | +0.25 SD |

The planet's identity response gains more than either covariance response.
Consequently, the raw covariance-minus-identity penalty changes from
`-0.5534` with integer templates to `-0.6921` with shifted templates; the
radial penalty changes from `-0.5979` to `-0.7004`.

### Shifted-filter closure

| Comparison | Injections | Range | Planet | Planet deviation | Predictive two-sided p |
| :--- | ---: | :--- | ---: | ---: | ---: |
| Raw covariance - identity | -0.2234 +/- 0.2196 | [-0.5226, +0.1944] | -0.6921 | -2.13 SD | 0.0649 |
| Radial covariance - identity | -0.1715 +/- 0.2648 | [-0.5766, +0.2412] | -0.7004 | -2.00 SD | 0.0813 |
| Gaussian 3.6 - identity | -0.0281 +/- 1.2056 | [-1.5064, +2.0770] | +0.5676 | +0.49 SD | 0.6443 |

For each covariance comparison, zero injections are below the planet. The
one-sided exchangeable lower-tail rank is therefore the minimum available
with twelve injections, `1/13 = 0.0769`. This is suggestive but cannot reach a
5% exchangeable threshold with the present sample.

### Absolute fixed-pixel SNR

| Method | Injection mean +/- sample SD | Planet |
| :--- | ---: | ---: |
| Gaussian 2.4 | 4.8621 +/- 1.1576 | 5.2696 |
| Gaussian 3.6 | 4.5990 +/- 1.3375 | 5.0869 |
| Integer identity | 4.5600 +/- 1.0684 | 4.2514 |
| Shifted identity | 4.6271 +/- 1.0872 | 4.5193 |
| Integer raw covariance | 4.3457 +/- 0.9195 | 3.6979 |
| Shifted raw covariance | 4.4037 +/- 0.9675 | 3.8273 |
| Integer radial covariance | 4.3608 +/- 0.8930 | 3.6535 |
| Shifted radial covariance | 4.4556 +/- 0.9405 | 3.8189 |

Gaussian 2.4 has the highest mean injection SNR and the highest planet SNR in
this predeclared comparison. Gaussian 3.6 remains statistically consistent
with shifted identity across the injections, while it exceeds shifted
identity by `0.568` SNR for the planet.

## Interpretation

Stage H established that cubic registration fixes most of the response-shape
mismatch. Stage I now shows that applying that correction consistently to the
source statistic and annular noise distribution does not explain the poorer
covariance performance. The remaining discrepancy is more consistent with
local planet noise, residual planet morphology, or another planet-specific
effect than with integer-grid template centering.

The result does not select a method from the single planet. Across the twelve
injections, Gaussian 2.4 is best on average, shifted identity is next, and the
two shifted covariance methods remain lower. A larger optimized-phase
injection sample would be required to resolve the planet's extreme rank below
the current `1/13` limit.

## Verification and archived products

The frozen completion verifier reran on ROC. All 245 shifted calibration fits
replayed their parent weights exactly, and 189 radius-12 policy positions were
applied per science image. The maximum production-versus-oracle SNR error was
`9.54e-7`; every unchanged injection statistic reproduced Stage G exactly.

The compact archive contains:

- [aggregate machine-readable result](results.json)
- [primary table](results.csv)
- [runner report](results.md)
- [verification record](verification.json)
- [frozen protocol](protocol.json)
- [frozen manifest](manifest.json)
- [completion receipt](complete.json)
- [shifted calibration weights](calibration/weights.npz)
- `analysis/`: per-task results, commands, logs, annular checks, and receipts

The 26 analysis FITS maps are omitted from the compact repository archive.
Their hashes and completion receipts were verified on ROC by the frozen
runner before the metadata products were copied.
