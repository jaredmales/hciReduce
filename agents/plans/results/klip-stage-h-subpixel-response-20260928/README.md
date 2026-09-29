# KLIP Stage-H subpixel-response result

## Result

Stage H completed all twelve negative reductions and paired them with the
existing nominal optimized-phase Stage-G positives. The frozen completion
verifier passed again on ROC, and the compact archive verifies all aggregate
products, analysis receipts, response stamps, and negative-reduction receipts.

The experiment cleanly identifies registration as the dominant response-model
error. At mode 200, cubic shifting the existing 47-pixel integer response to
the known `(-0.2771,+0.4860)`-pixel phase raises response cosine from about
`0.92` to `0.99` and reduces the best-scaled residual from about
`0.37-0.38` to `0.11-0.14`.

### Paired central-response target

Values are means plus or minus sample standard deviations over the twelve
sites.

| Covariance metric | Template | Cosine | Projection scale | Best-scaled residual |
| :--- | :--- | ---: | ---: | ---: |
| Identity | Current integer | 0.925752 +/- 0.034646 | 0.91848 +/- 0.02944 | 0.3670 +/- 0.0888 |
| Identity | Cubic shifted | 0.993552 +/- 0.003880 | 0.99039 +/- 0.01215 | 0.1072 +/- 0.0382 |
| Raw rectangular | Current integer | 0.920476 +/- 0.035519 | 0.91130 +/- 0.03068 | 0.3803 +/- 0.0868 |
| Raw rectangular | Cubic shifted | 0.990793 +/- 0.006183 | 0.98881 +/- 0.01783 | 0.1269 +/- 0.0488 |
| Radial Hann | Current integer | 0.919875 +/- 0.034450 | 0.91108 +/- 0.03018 | 0.3825 +/- 0.0840 |
| Radial Hann | Cubic shifted | 0.989325 +/- 0.007182 | 0.98714 +/- 0.01931 | 0.1368 +/- 0.0520 |

Every site improves under every primary comparison. Depending on covariance
metric, the mean cosine gain is `0.0678-0.0703`, residual reduction is
`0.2457-0.2597`, and absolute projection-error reduction is
`0.0685-0.0712`. All 108 site-by-metric primary improvements are positive.

The result is stable across KL mode. From modes 125 through 350, the shifted
mean cosine ranges are:

- identity: `0.99337-0.99382`;
- raw rectangular: `0.99026-0.99117`; and
- radial Hann: `0.98850-0.98975`.

## One-sided response control

The regenerated paired-exact response closely models the actual positive
injection. At mode 200 its positive-response cosine is `0.99982`,
`0.99948`, and `0.99908` under identity, raw rectangular, and radial-Hann
covariance. The corresponding residuals are `0.0178`, `0.0296`, and
`0.0388`, with projection scales consistent with one.

Finite-amplitude asymmetry is therefore small. It does not explain the
`0.107-0.137` residual left by cubic shifting; a smaller genuine
phase-dependent shape difference remains.

## Interpretation

Cubic registration is the appropriate first subpixel-aware template model. It
reduces the theoretical matched-filter SNR loss from roughly 7-8% to about
0.6-1.1% and reduces the noiseless contrast bias from roughly 8-9% to about
1-1.3%. A measured fractional-phase response grid could recover the remaining
shape difference, but the likely detection-SNR gain over cubic shifting is
small.

The response correction is similar under identity and covariance weighting.
It therefore does not explain why covariance filters underperform identity in
the injection study. It should primarily improve fixed-pixel localization and
photometry and reduce dependence on the three-pixel aperture maximum.

The next useful application is to use the cubic-shifted response for the known
planet and comparable optimized-phase injections, keeping the established
annular-noise and aperture rules fixed. That reanalysis can determine whether
the improved response centering removes the fixed-pixel planet discrepancy and
changes the Gaussian-3.6 comparison.

## Archived products

- [aggregate machine-readable result](results.json)
- [primary table](results.csv)
- [runner report](results.md)
- [verification record](verification.json)
- [frozen protocol](protocol.json)
- [frozen manifest](manifest.json)
- [completion receipt](complete.json)
- `analysis/`: per-site all-mode metrics and mode-200 response stamps
- `negative/`: paired-negative receipts and frozen commands
