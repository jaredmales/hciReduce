# KLIP Stage-B known-planet footprint correction

## Finding

Radius 12 is effectively the fitted planet radius: the optimized planet lies at
11.783 pixels. The first Stage-B screens accounted for it in covariance
training, radial-profile estimation, and a seven-pixel center exclusion. They
did not mask outer pixels of candidate stamps and exact responses when those
pixels entered the planet disk. The 31- and 47-pixel candidate comparisons from
those screens therefore violated the intended conservative exclusion contract.

## Complete-footprint audit

Entries are complete five-query sites; parenthetical values are individual
query footprints. Site selection is frozen and score-blind.

| Radius | 11-pixel response | 31-pixel response | 47-pixel response |
| ---: | ---: | ---: | ---: |
| 7.5 | 7/12 (41/60) | 0/12 (0/60) | 0/12 (0/60) |
| 10 | 9/12 (45/60) | 0/12 (2/60) | 0/12 (0/60) |
| 12 | 9/12 (46/60) | 2/12 (12/60) | 0/12 (0/60) |
| 16 | 9/12 (53/60) | 4/12 (20/60) | 0/12 (0/60) |
| 20 | 10/12 (50/60) | 4/12 (23/60) | 2/12 (10/60) |
| 24 | 10/12 (54/60) | 6/12 (35/60) | 2/12 (11/60) |

A complete-footprint rejection would leave too little common angular coverage,
especially for the larger supports. The corrected calculation instead retains
every original site, assigns zero weight to native candidate and response
coordinates inside the fixed planet disk, restricts the PSD covariance to the
same valid rows and columns, and solves that positive principal submatrix. It
renormalizes every result to unit response on the retained coordinates.

## Corrected raw screen

These are median candidate-score variances over the controlling radii 7.5, 10,
and 12 for the narrow training bands.

| Method | 11 pixels | 31 pixels | 47 pixels |
| :--- | ---: | ---: | ---: |
| identity | 6.512 | 7.131 | 6.879 |
| rectangular, mixing 0.1 | 2.597 | 2.989 | 2.952 |
| rectangular, mixing 0.3 | 2.867 | 3.324 | 3.306 |
| Hann, mixing 0.1 | 2.477 | 2.875 | 2.698 |
| Hann, mixing 0.3 | 2.729 | 3.223 | 3.066 |

Masking raises rather than lowers every primary variance relative to the first
screen. The 11-pixel response is now best among the three supports for the
leading Hann/mixing-0.1 family. Its split-weight cosine remains 0.991. Median
retained template energy is 1.000, 0.989, and 0.983 for supports 11, 31, and
47, respectively, so the result is not driven by removing most of the response.

## Corrected radial normalization

The strict leave-site-out radial profile and the fixed candidate planet mask
are applied together on the supported 11-pixel geometry.

| Method | Masked raw variance | Masked normalized variance | Normalized/raw |
| :--- | ---: | ---: | ---: |
| identity | 6.512 | 5.316 | 0.818 |
| rectangular, mixing 0.1 | 2.597 | 2.204 | 0.849 |
| rectangular, mixing 0.3 | 2.867 | 2.412 | 0.841 |
| Hann, mixing 0.1 | 2.477 | 2.092 | 0.844 |
| Hann, mixing 0.3 | 2.729 | 2.280 | 0.835 |

For Hann/mixing-0.1, radius-12 variance changes from 3.897 raw to 3.153 after
normalization. It remains the largest radial mismatch. Restricting the original
unmasked calculation to the nine fully planet-clear 11-pixel sites makes the
radius-12 values larger still: 4.806 raw and 3.898 normalized. Planet overlap
therefore does not explain the radius-12 behavior, although it invalidated the
earlier larger-support comparison.

The subtraction residual RMS inside the seven-pixel disk is 0.1704, comparable
to 0.1762 in the planet-free portion of the 11.5--13.5-pixel annulus. This is a
descriptive check rather than grounds to relax the exclusion.

## Decision

Use only these planet-masked values for subsequent Stage-B decisions. Retain
radial-normalized 11-pixel Hann/mixing-0.1 as the leading data-selected PSD arm
and raw 11-pixel rectangular/mixing-0.3 as the mandatory P4 prior. Neither is
calibrated well enough for promotion: the Hann primary variance is 2.092 and
its radius-12 variance is 3.153. The next discriminator remains the prescribed
post-ensemble-mean patch-RMS control on this corrected 11-pixel geometry.

## Reproduction and verification

The audit is implemented by
[`audit_klip_stage_b_planet_footprints.py`](../../scripts/audit_klip_stage_b_planet_footprints.py),
the corrected raw screen by
[`run_klip_stage_b_local_planet_masked_screen.py`](../../scripts/run_klip_stage_b_local_planet_masked_screen.py),
and the corrected normalization by
[`run_klip_stage_b_local_planet_masked_normalization.py`](../../scripts/run_klip_stage_b_local_planet_masked_normalization.py).
The raw screen completed 10,800 records in 228.6 seconds; the normalized screen
completed 1,800 records. Dense principal-submatrix comparisons at 11 and 31
pixels, a masked isotropic endpoint, CG residual checks, geometry replay, and
saved-product hashes pass.

Compact machine-readable results are in [`primary.json`](primary.json), and
source-product receipts are in [`verification.json`](verification.json). The
figures show the [clear-site audit](clear-site-audit.png),
[masked raw comparison](masked-raw-comparison.png), and
[masked normalized comparison](masked-normalized-comparison.png). Full record
tables remain under ignored `working/local` data.
