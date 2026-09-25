# KLIP Stage-B direct diagonal/PCA covariance screen

## Question and contract

Can a finite empirical covariance estimated directly from the supported
11-pixel training patches calibrate KLIP candidate scores better than the
stationary PSD models? This mode-200 screen retains the corrected candidate
planet mask, all frozen score-blind sites, and independent detector-half fits.
It evaluates both raw and strict radial-standardized coordinates.

The fixed grid is:

- diagonal covariance using per-pixel sample variances and the production 0.1
  median-variance floor;
- PCA caps of 0, 3, 8, and every estimable centered sample mode;
- isotropic PCA floors of 0.1, 0.3, and 1.0 times the median pixel variance;
- the narrowest split-supported radial training band and the next wider fixed
  band; and
- candidate evaluation both with and without the fitted training mean.

Every candidate/template coordinate inside the fixed seven-pixel known-planet
disk has zero weight. The direct diagonal-plus-low-rank covariance is restricted
to the same valid principal submatrix before solving. The local solver copies
the production diagonal scaling and QR/Cholesky algebra and is checked against
dense principal-submatrix solves.

## Controlling-radius result

The table shows narrow-band median candidate-score variance over radii 7.5, 10,
and 12. Each PCA row uses the best of its predeclared floors, which is 1.0 in
every case; the complete floor grid remains in `primary.json`.

| Coordinates | Model | Variance | Fitted-mean variance | Opposite-half variance | Split weight cosine |
| :--- | :--- | ---: | ---: | ---: | ---: |
| raw | diagonal, floor 0.1 | 6.786 | 6.803 | 4.105 | 0.954 |
| raw | PCA, 0 modes, floor 1.0 | 6.615 | 6.618 | 3.999 | 1.000 |
| raw | PCA, 3 modes, floor 1.0 | 5.874 | 5.912 | 3.187 | 0.887 |
| raw | PCA, 8 modes, floor 1.0 | 4.455 | 4.434 | 2.555 | 0.777 |
| raw | PCA, all modes, floor 1.0 | 3.439 | 3.418 | 1.846 | 0.742 |
| radial | diagonal, floor 0.1 | 5.636 | 5.655 | 4.084 | 0.955 |
| radial | PCA, 0 modes, floor 1.0 | 5.406 | 5.410 | 3.938 | 1.000 |
| radial | PCA, 3 modes, floor 1.0 | 4.824 | 4.845 | 3.202 | 0.881 |
| radial | PCA, 8 modes, floor 1.0 | 3.800 | 3.781 | 2.559 | 0.770 |
| radial | PCA, all modes, floor 1.0 | 2.886 | 2.866 | 1.820 | 0.730 |

The corrected reference variances are 6.512 for raw identity, 2.867 for raw
rectangular/mixing-0.3, 2.477 for raw Hann/mixing-0.1, and 2.092 for radial
Hann/mixing-0.1. Diagonal covariance is slightly worse than identity. Adding
more PCA modes improves direct covariance, but the best direct model remains
38% above radial Hann and has much poorer split-weight stability.

Lower PCA floors are strongly miscalibrated. With every estimable mode, floor
0.3 gives variance 9.660 raw and 8.079 radial; floor 0.1 gives 29.096 and
24.268. Those floors assign too little conditional variance to the large
unresolved complement. The floor-1.0 all-mode model retains a median 27 modes
across the controlling radii, with a range of 20--30.

## Radius dependence

| Radius | Radial PCA, all modes, floor 1.0 | Radial Hann/mixing-0.1 parent |
| ---: | ---: | ---: |
| 7.5 | 2.742 | 1.653 |
| 10 | 2.886 | 2.092 |
| 12 | 5.004 | 3.153 |
| 16 | 4.072 | 1.554 |
| 20 | 2.341 | 1.123 |
| 24 | 2.661 | 0.958 |

Direct PCA is worse at every tested radius. It does not resolve the radius-12
mismatch and transfers particularly poorly at radii 16--24.

## Wider-band stability

For the all-mode floor-1.0 model, the wider band changes raw variance from
3.439 to 3.813 and radial variance from 2.886 to 2.909. The radial held-out
variance improves from 1.820 to 1.359 and split-weight cosine from 0.730 to
0.804, but fixed candidate calibration remains worse than radial Hann.

The strict radial profile is undefined at native radius 60 and beyond. At
radius 12, some wider-band stamps touch that boundary; the runner rejects those
geometry-defined patches rather than extrapolating the profile. Depending on
query and detector half, 0--31 patches are rejected while 110--160 training
patches remain. No primary narrow-band patch is rejected.

Fitted candidate-mean subtraction is neutral throughout the useful part of the
grid.

## Decision

Do not advance direct diagonal or PCA covariance. Direct PCA improves as the
floor and rank increase, but its best policy remains worse calibrated and much
less stable than the stationary PSD candidate. Keep radial-standardized
11-pixel Hann/mixing-0.1 as the leading data-selected arm and raw 11-pixel
rectangular/mixing-0.3 as the mandatory P4 prior. The remaining Stage-B loose
end is the larger-response support radial-mean control.

## Reproduction and verification

The maintained runner is
[`run_klip_stage_b_local_direct_covariance.py`](../../scripts/run_klip_stage_b_local_direct_covariance.py).
It completed 18,720 policy/query records. Deterministic checks cover every
model rank/floor, preservation of the centered sample trace, restricted
scaled-QR solves against dense covariance principal submatrices, the candidate
planet mask, and exact unit physical source response. All input and output
receipts pass.

Compact results are in [`primary.json`](primary.json), source-product hashes are
in [`verification.json`](verification.json), and the complete PCA grid is shown
in [`comparison.png`](comparison.png). The 58.7-MB record table remains under
ignored `working/local` data.
