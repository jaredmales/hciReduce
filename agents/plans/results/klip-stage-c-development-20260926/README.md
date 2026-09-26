# KLIP Stage-C development result

## Status and provenance

The calibration-first Stage-C development campaign completed on ROC on
2026-09-26. The final state is `development_complete`; all 48 calibration
units, 120 calibration sites, 108 positive KLIP reductions, and 108 analysis
tasks have valid nested product fingerprints. The final receipt contains 113
products, including the analysis-repair manifest. No validation product or
held-out null was opened, and the generated result explicitly records
`selection_performed: false`.

The full preserved products are:

- [machine-readable result](results.json), with all 3,744 summary rows;
- [flat summary table](results.csv);
- [generated mode-200 and refit tables](results.md);
- [mode-200 recovery and SNR plot](mode200.png);
- [completion receipt](complete.json), [analysis-runner manifest](analysis_runner_manifest.json),
  and [analysis-repair record](analysis_repair.json); and
- [independent preservation check](verification.json).

The analysis repair archived 81 complete and 27 incomplete products from the
first analysis attempt, retained all calibration and reduction receipts, and
regenerated every analysis under one frozen runner. A recursive audit after
completion verified every calibration, reduction, analysis, and final-product
fingerprint. The largest production-versus-independent annular SNR difference
was 2.86e-6 under the fixed `rtol=2e-6`, `atol=2e-6` comparison.

## Experiment

Each of six radii (7.5, 10, 12, 16, 20, and 24 pixels) has six development
sites injected at physical contrasts targeting Gaussian-FWHM-3.6 source-only
SNR 3, 5, and 7 at KL mode 200. Detection uses the largest of the center and
four adjacent pixels. Every radius, mode, and method has its own threshold,
frozen as the largest score among 20 calibration sites before a positive image
was opened.

Mode 200 is the prespecified primary comparison. The other seven KL mode counts
are secondary outputs and cannot be searched jointly without a new joint-null
calibration. The table below pools the six radii, so each recovery count is out
of 36 development injections at a given brightness. `Mean SNR 3` averages the
six radius-level means. Paired wins compare the 36 individual SNR-3 searches
with Gaussian FWHM 3.6.

## Primary mode-200 result

| Method | Recovered at 3 | Recovered at 5 | Recovered at 7 | Mean SNR 3 | Paired wins vs G3.6 | Mean SNR change vs G3.6 |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| Sparse-response identity | **32** | 36 | 36 | 3.6216 | 25/36 | +0.4204 |
| Exact-response identity | 29 | 36 | 36 | **3.6764** | **26/36** | **+0.4752** |
| Radial Hann, truncate at 0.75 mean variance | 28 | 36 | 36 | 3.6519 | 24/36 | +0.4507 |
| Radial Hann, full precision | 28 | 36 | 36 | 3.5942 | 24/36 | +0.3931 |
| Raw rectangular, mixing 0.3 | 28 | 36 | 36 | 3.5714 | 23/36 | +0.3702 |
| Gaussian FWHM 2.4 | 27 | 35 | 36 | 3.5249 | 30/36 | +0.3237 |
| Gaussian FWHM 3.6 | 27 | 35 | 36 | 3.2012 | -- | -- |
| Native pixel | 23 | 34 | 35 | 3.0796 | -- | -0.1216 |

The exact-response matched filter is therefore measured correctly and improves
on Gaussian 3.6 in this development set. Gaussian 2.4 is a materially stronger
smoothing control than Gaussian 3.6. Relative to Gaussian 2.4, exact identity
wins 21/36 paired SNR-3 trials with mean difference +0.1515; radial truncation
at 0.75 also wins 21/36 with mean difference +0.1270.

Covariance weighting does not improve on the identity matched filters. The
best covariance policy gains one SNR-3 recovery over Gaussian 2.4 but loses one
to exact identity and four to sparse identity. Its small gain over Gaussian 2.4
is still suitable for the prespecified held-out validation, but the development
result does not support an optimal-covariance claim.

## Separation dependence at target SNR 3

Each cell is `recoveries / 6; mean maximum SNR`.

| Radius | Gaussian 3.6 | Gaussian 2.4 | Exact identity | Sparse identity | Raw rectangular | Radial truncation 0.75 |
| ---: | :--- | :--- | :--- | :--- | :--- | :--- |
| 7.5 | 4; 3.110 | 5; 3.547 | **6; 4.348** | 6; 4.120 | 6; 4.092 | 6; 4.255 |
| 10 | 2; 2.332 | 2; 2.696 | 4; 2.966 | **5; 2.818** | 4; 3.084 | 4; 3.166 |
| 12 | 6; 3.992 | **6; 4.110** | 6; 3.823 | 6; 3.810 | 6; 3.513 | 6; 3.496 |
| 16 | 5; 3.555 | 5; 3.818 | 6; 4.056 | **6; 4.114** | 6; 3.826 | 6; 3.875 |
| 20 | **6; 2.855** | 5; 3.174 | 4; 2.944 | 4; 2.894 | 3; 2.994 | 3; 2.991 |
| 24 | 4; 3.363 | 4; 3.805 | 3; 3.922 | **5; 3.974** | 3; 3.919 | 3; 4.129 |

The matched filters provide their clearest recovery advantage at 7.5 and 10
pixels, the separations that motivated this work. The picture is mixed at 20
and 24 pixels because the method-specific maximum-null thresholds are driven
by large calibration outliers. At radius 24, the largest-minus-second-largest
calibration-score gap is 1.94 for exact identity and 2.18--2.23 for the listed
covariance policies. The held-out nulls are needed to determine whether these
outer-radius thresholds are conservative or unstable; they cannot be revised
after those nulls are opened.

## Covariance refit and regularization

Across all eight modes, six radii, and three brightnesses, per-image covariance
refitting changes no recovery for seven of the eight covariance policies. It
loses one recovery for radial truncation at 0.5. Mean refit-minus-frozen SNR
changes range from -0.00070 to +0.00244, and the largest absolute cell change is
0.0492. Retain baseline-frozen covariance for validation; the refit arm adds
cost without a measurable detection benefit.

The radial precision policies are also nearly degenerate. At mode 200 and SNR
3 they all recover 28/36; their mean maximum SNR spans 3.5793--3.6519. Hard
truncation at 0.75 has the largest mean SNR and a calibration-gap audit comparable to
full precision, so it is the development-selected radial policy. Raw
rectangular/mixing-0.3 remains the mandatory independent P4-prior family.

Mode 300 gives slightly larger development SNR for most matched filters, but
this is a secondary result. For example, exact identity rises from 29 to 30
recoveries and radial truncation 0.75 rises from 28 to 29. Mode 200 remains the
confirmatory mode because it was fixed before development and the existing
thresholds do not calibrate a maximum over modes.

## Response fidelity and photometry

The measured exact response agrees closely with the source contribution in the
positive reductions. At mode 200, unweighted fidelity over all 108 images has
median cosine 0.999844, minimum cosine 0.998312, median projection 1.00686, and
maximum best-scaled residual 0.0581. The radius-7.5 median projection is
0.97324; the other radii are within about 1.3 percent of unity. The covariance
metrics emphasize small residual differences more strongly, reaching maximum
relative residual 0.104 for raw rectangular and 0.147 for radial Hann, but they
do not indicate a response-registration failure.

At SNR 3, the radius-pooled median relative contrast bias is +3.9 percent for
exact identity, +4.4 percent for sparse identity, -3.5 percent for radial
truncation 0.75, and -2.0 percent for Gaussian 3.6. Mean per-radius scatter is
19--21 percent for the matched filters and 32 percent for Gaussian 3.6, as
expected for these low-SNR measurements. Throughput is near unity except at
radius 7.5, where it falls from about 0.99 at the faint level to about 0.95 at
the bright level.

## Recommended frozen validation shortlist

Before opening held-out nulls or validation images, freeze this mode-200,
baseline-weight policy:

1. required controls: native, Gaussian FWHM 3.6, Gaussian FWHM 2.4, and exact-response identity;
2. leading response reference: sparse-response identity;
3. covariance families: raw rectangular/mixing-0.3 and radial-standardized Hann/mixing-0.1 with hard truncation at 0.75 mean variance; and
4. the already frozen method-specific maximum-of-20 calibration thresholds and five-pixel search.

The validation gate should retain the original paired comparison against
Gaussian 3.6 and additionally report Gaussian 2.4, now known to be the stronger
smoothing control. A covariance method should not be called successful unless
it improves paired validation SNR or recovery without increasing held-out-null
exceedances. Sparse identity and exact identity should be evaluated as leading
non-covariance outcomes rather than evidence for covariance weighting.

## Limits

This is a development selection result from six sites per radius in one
observing sequence. The 20 calibration maxima are correlated spatial samples,
and target SNR 5 and 7 are nearly saturated. No held-out false-positive rate,
validation completeness, known-planet statistic, second epoch, or second
target has entered the conclusions above. The three `r7p5_dev03` brightnesses
used the audited inner-annulus boundary substitution in both covariance arms;
all validation sites have complete annular bracketing.
