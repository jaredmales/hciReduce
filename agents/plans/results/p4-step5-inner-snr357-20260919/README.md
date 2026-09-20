# Step 5: completed inner SNR-3/5/7 reductions and annular-mask diagnosis

## Completion

The study completed all **108 full P4 reductions and measurements**. Its final
receipt records unchanged frozen inputs and all 108 job names. The final
summary initially stopped on a field-name mismatch after the measurements were
complete; repair `repair_0001.json` verified every measurement product and
updated only the frozen runner. The resumed run reused all reductions and
analyses and completed normally.

The source contrasts were frozen from baseline identity amplitude-map noise so
that the nominal source-only SNRs were 3, 5, and 7. The SNR-5 (`1×`) median
contrasts from radius 6 through 24 pixels were `3.5031e-2`, `1.7395e-2`,
`7.1493e-3`, `2.9848e-3`, `1.3781e-3`, and `1.5420e-3`.

## The reported SNR did not reach the target

The completed analysis excluded AF Lep b from the annular-noise profile but
did not exclude the injected source. Identity's mean fixed-center SNRs were:

| Radius (pixels) | Nominal SNR 3 | Nominal SNR 5 | Nominal SNR 7 |
| ---: | ---: | ---: | ---: |
| 6 | 1.709 | 1.839 | 1.879 |
| 8 | 1.905 | 2.152 | 2.224 |
| 12 | 2.372 | 2.729 | 2.897 |
| 16 | 2.273 | 3.226 | 3.704 |
| 20 | 2.550 | 3.710 | 4.389 |
| 24 | 2.514 | 3.828 | 4.660 |

The failure to increase proportionally is caused primarily by the injection
entering its own one-pixel annular standard-deviation estimate. Relative to the
matching baseline-site profile, the mean identity noise-sigma ratios were:

| Radius (pixels) | Nominal SNR 3 | Nominal SNR 5 | Nominal SNR 7 |
| ---: | ---: | ---: | ---: |
| 6 | 1.898 | 2.779 | 3.604 |
| 8 | 1.590 | 2.346 | 3.108 |
| 12 | 1.591 | 2.113 | 2.658 |
| 16 | 1.164 | 1.456 | 1.815 |
| 20 | 1.117 | 1.319 | 1.582 |
| 24 | 1.084 | 1.245 | 1.462 |

This is not primarily lost matched-filter throughput. The mean paired identity
amplitude increment divided by injected contrast remains 0.91–1.05 over all
radii and levels. At radii 8–24 it is within about 5% of unity; the largest
departure is 0.91 at radius 6 and nominal SNR 7. The source response therefore
survives while its inclusion in the noise sample makes the reported SNR
saturate, most strongly where an injected PSF occupies a large fraction of the
small annulus.

The recovery counts from this unmasked normalization do not answer the intended
SNR-3/5/7 comparison and are not used to choose a filter.

## Corrected saved-image analysis

No new P4 reductions are required. The maintained raw/patch-RMS comparison can
reanalyze all 164 baseline searches and 108 positive images while excluding
both the known planet and the current trial from each production annular-noise
profile. Baseline searches receive the same current-trial exclusion as their
corresponding positive measurements, and new thresholds are frozen from those
masked baseline scores before positives are read.

`hciAnalyze` accepts vector-valued `planet.sep`, `planet.PA`, and `planet.R`, so
the corrected path supplies the known planet and current trial without changing
production code. An isolated radius-12 replay using two 7.3-pixel radii matched
an independent two-circle annular oracle to `2.72e-7`; its SNR-5 center
measurement was 4.126.

## Corrected-analysis launch repair

The first corrected launch stopped before calibration with 97 of 164 baseline
receipts complete and no positive analysis started. A verification-only raw
control diagnostic attempted a maximum difference even when the masked and
parent SNR maps had no common finite pixels, producing a NaN that strict JSON
rejected. The maintained runner now records that unavailable diagnostic as
JSON `null`. The exact-error repair retains all 97 completed receipts; restart
archives and recomputes the 67 unreceipted task directories.

The resumed 7.3-pixel analysis reached 140 baseline receipts, then showed that
the two large masks exhaust the usable annulus at four radius-6 sites, two
radius-8 sites, and 18 central calibration candidates. That failed root is
retained and will not be repaired into a different design.

## Lambda/D trial exclusion

A read-only sweep of all saved amplitude maps selected `planet.R=3.1` for the
trial. With production's half-pixel mask boundary, this excludes through 3.6
pixels, exactly one lambda/D. It is the largest fixed mask that gives complete
five-pixel identity and Gaussian support at all 36 injection sites. The known
planet retains its established 7.3-pixel exclusion.

Five of the 128 calibration candidates have no active method with these two
masks. They remain explicit invalid candidates. Each active method's 20-null
pool is reselected without scores from the narrowest radial band with enough
valid trials; raw and post-mean patch-RMS pairs require common support and use
the same maximin-selected locations.

The saved identity maps predict the following mean five-pixel-search / fixed-
center SNR after the lambda/D trial exclusion:

| Radius | Nominal 3 | Nominal 5 | Nominal 7 |
| ---: | ---: | ---: | ---: |
| 6 | 3.250 / 2.945 | 4.547 / 4.027 | 6.806 / 6.163 |
| 8 | 3.208 / 2.646 | 3.775 / 3.442 | 4.065 / 3.782 |
| 12 | 3.372 / 3.155 | 4.190 / 4.128 | 4.845 / 4.791 |
| 16 | 3.017 / 2.613 | 4.668 / 4.254 | 5.806 / 5.521 |
| 20 | 3.358 / 2.954 | 4.952 / 4.776 | 6.302 / 6.272 |
| 24 | 2.845 / 2.758 | 4.592 / 4.592 | 6.162 / 6.162 |

The radius-8 high levels still compress, and the radius-6 profiles retain only
3--6 native annular pixels after both masks. The final comparison therefore
reports measured SNR and sample support rather than treating 3/5/7 as achieved
values. At the most constrained supported radius-6 site, the production SNR
map and independent oracle agreed exactly on their common finite pixels. The
comparison reuses all 108 positive images and performs no new P4 reduction.

## Completed lambda/D comparison

The replacement comparison completed all **164 baseline analyses and 108
positive analyses** with no new P4 reductions. Its completion receipt verifies
that every frozen input is unchanged. Production annular SNR and the independent
oracle agree exactly on every common finite pixel.

The score-free pool reselection found 20 valid nulls for every active method,
but the smallest separations require wider calibration ranges. Radius-6 ±20
and the references use a ±2-pixel band with five replacements. Radius-8 ±10
uses ±4 with 13 replacements, while radius-8 ±20 and the references use ±1.
Radius-12 ±5 uses ±4 with 14 replacements. All raw/patch-RMS pairs use exactly
the same locations, so their paired comparison remains controlled; thresholds
between different widths can reflect these different radial null populations.

The methods with complete support at all 36 sites gave these aggregate recovery
counts:

| Method | Valid sites | Nominal 3 | Nominal 5 | Nominal 7 | Held-out nulls |
| --- | ---: | ---: | ---: | ---: | ---: |
| Raw rectangular ±20 | 36/36 | 19 | 33 | 34 | 0/36 |
| Patch-RMS rectangular ±20 | 36/36 | 20 | 33 | 34 | 0/36 |
| Identity matched filter | 36/36 | 23 | 35 | 35 | 0/36 |
| Gaussian FWHM 3.6 | 36/36 | 20 | 34 | 36 | 1/36 |

Post-mean patch-RMS normalization has no consistent recovery advantage. Across
all widths, radii, and levels it changes only four paired decisions: patch RMS
alone recovers one radius-6, nominal-3 injection at ±20 and one radius-16,
nominal-5 injection at ±10; raw alone recovers two radius-12, nominal-3
injections at ±5. At ±20 the mean patch-RMS-minus-raw search-SNR difference is
slightly negative at radii 6--12 and positive at radii 16--24, reaching about
`+0.30` at radius 20. The sign reversal and unchanged medium/high recovery do
not support adopting the normalization.

The amplitude response does not explain the recovery result. Mean paired
positive-minus-baseline throughput for raw and patch-RMS ±20 differs by less
than 0.01 at every radius and level. Both are near unity outside the innermost
sites. At radius 6, raw throughput declines from 0.984 to 0.908 over nominal
3--7, while patch RMS declines from 0.985 to 0.909; this common trend is
finite-source P4 behavior rather than covariance normalization.

Identity gives the best aggregate recovery at nominal 3 and 5 without a held-
out null exceedance. Gaussian reaches 36/36 at nominal 7 but has one null
exceedance. The raw ±20 covariance filter trails identity by four detections at
nominal 3 and two at nominal 5. With six correlated sites per radius, this is
development evidence rather than an independent filter ranking, but it gives
no evidence that the tested covariance weighting improves on identity.

Annular support limits the innermost result. At radius 6 the required profiles
contain as few as **2 native pixels** and have a median minimum of 3.5; at
radius 8 the corresponding values are 7 and 8. The nominal SNR sequence also
compresses most strongly at radius 8: identity's mean five-pixel search SNR is
3.208, 3.775, and 4.065. Radius-6 recovery and thresholds should therefore be
treated as sparse-sample diagnostics rather than stable completeness estimates.
