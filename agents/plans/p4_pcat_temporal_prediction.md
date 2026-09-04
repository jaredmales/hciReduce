# P4 PCAT Temporal Prediction Plan

## Implementation status

The initial CPU detector-frame implementation is complete in the working tree. `P4TemporalPCA` learns a full
reference-series temporal covariance, caches one intercept-plus-mode least-squares factorization for each no-wrap
gap, and supplies exactly one PCAT prediction column to the ordinary P4 fit. The existing neighboring-frame
mean/median behavior remains the compatibility default and PCAT is selected explicitly with
`p4.temporalPredictor=pcat`.

Implemented scope:

- `or` and common-mask-valid detector-`annulus` PCAT reference regions;
- `pixelMean` and `none` reference-series centering;
- all-frame PCAT prediction with a configurable image-index gap, including endpoint handling;
- invalidation of the complete search-pixel P4 fit when a PCAT basis or outside-gap fit is unsupported;
- configuration and FITS provenance; and
- numerical leakage/gap tests plus reduction tests for both reference geometries.

Target-frame exclusion, pixel-local trials, frozen-PSF products, rotated regression, diagnostics, and scientific
performance validation remain deferred.

## Goal

Add an optional PCAT-style temporal predictor for each detector search pixel. Instead of using the current
`p4.numberImages>0` neighboring-frame mean/median summary, PCAT learns temporal eigenmodes from a configurable
spatial reference region and predicts the target pixel across a withheld interval in its own time series. That
prediction is appended as the one temporal-augmentation column in the ordinary P4 design matrix.

The central target sample must not influence its fitted PCAT coefficients. For target image `t`, PCAT fits only the
target-pixel samples outside a configurable temporal gap around `t`, then evaluates the fitted temporal model at `t`.
There is no wrap at the beginning or end of an observation.

For detector pixel `p`, ordinary P4 currently learns a spatial predictor vector `X_spatial[p,t,:]` and regresses the
target series `y_p[t]`. With PCAT enabled, its design row becomes

```text
X[p,t,:] = [ X_spatial[p,t,:], prediction_p(t) ].
```

PCAT therefore has the same role and one-column dimensional effect as the current positive-`p4.numberImages`
augmentation. P4 continues to learn and apply its ordinary spatial-plus-temporal regression coefficients; PCAT does
not directly subtract its prediction from the science image. A compact PCAT prediction/residual diagnostic may be
retained for validation, but is not a separate science product in the first implementation.

## Numerical model

For one detector search pixel `p`, let `R_p={r_1,...,r_N}` be its selected PCA reference pixels and let `T` be the
retained, ordered target-image count. The reference time-series matrix is

```text
Z[i,t] = I_t(r_i),                 shape N x T.
```

`R_p` must exclude the central target pixel and all geometry/mask-invalid samples. Its rows use the same
post-preprocessing detector-frame input cube as P4. Before PCA, each reference-pixel series is centered using only
the reference-region operation defined below; the centering rule and all retained reference-pixel identities are
recorded as provenance.

The time-domain covariance is

```text
C_time = Z^T Z,                    shape T x T.
```

Its leading orthonormal eigenvectors form temporal modes `U[:,0:m]`. For central pixel `p` with time series
`y_p`, and target image `t`, define a no-wrap temporal gap of half-width `g` images:

```text
G_t = { j : max(0,t-g) <= j <= min(T-1,t+g) }
S_t = {0,...,T-1} \ G_t.
```

The target fit includes an intercept plus the temporal modes. Its intercept is fitted only from the samples outside
the gap, so target-series centering does not leak the central sample or any other withheld sample. With
`B=[1,U[:,0:m]]`, solve

```text
b_p,t = argmin_b || B[S_t,:] b - y_p[S_t] ||_2,
prediction_p(t) = B[t,:] b_p,t.
```

The diagnostic PCAT residual is `y_p(t)-prediction_p(t)`. The predictor is an interpolation across the gap, not an
in-sample reconstruction. Endpoints naturally have one-sided training support. A target PCAT fit is invalid when the
retained outside-gap rows cannot support the requested basis numerically; it must never silently reduce the temporal
gap or change the requested mode count.

## Configurable science surface

The following implemented names are intentionally separate from
`p4.numberImages`, whose mean/median summary experiment did not improve the reduction.

| Key | Initial contract |
|---|---|
| `p4.temporalPredictor` | Selector: `none`, `neighborSummary` (default), or `pcat`. This makes the one-column augmentation explicit and prevents two competing temporal predictors in the first slice. |
| `p4.pcatGapImages` | Nonnegative integer temporal half-width `g`; the excluded interval contains `2g+1` images away from endpoints. `0` excludes only the central sample. |
| `p4.pcatReferenceRegion` | Reference-region selector. Initial supported values: `or` (the existing valid optimization region) and `annulus` (a configurable detector annulus with the current target's signal mask removed). |
| `p4.pcatReferenceMinRadius`, `p4.pcatReferenceMaxRadius` | Required ordered detector radii when `pcatReferenceRegion=annulus`; define a full annular reference region. |
| `p4.pcatModeFraction` | One finite fraction in `(0,1]`, realized against the PCAT numerical rank. It chooses the temporal basis size for the one appended PCAT predictor column and is separate from spatial-P4 `p4.modeFractions`. |
| `p4.pcatCentering` | Initial choices `pixelMean` and `none`; default `pixelMean`. This is the centering applied to `Z` only, not a second preprocessing pass on the science cube. The target fit always has its own outside-gap intercept. |

The first slice uses gap width in ordered retained images because it exactly describes the available samples and keeps
endpoint behavior unambiguous. A later physical-time or parallactic-angle gap may be added only after defining its
behavior for irregular cadence and observation gaps; it must not replace the image-count meaning silently. The first
PCAT implementation rejects `p4.numberImages>0`; `neighborSummary` remains the explicit existing alternative.

## Reference-region and masking contract

- The PCA reference region is configurable per detector search pixel, but it must be common to all temporal samples
  for that pixel. A pixel masked in any retained image is excluded from `R_p`; PCAT does not fill reference-region
  masks.
- `pcatReferenceRegion=or` uses the existing P4 optimization-region geometry after its signal-exclusion policy and
  common mask, omitting the central search pixel. This gives a direct first comparison to spatial P4.
- `pcatReferenceRegion=annulus` uses an independently configured full detector annulus. It must obey the same common
  mask and exclude the configured physical signal region around the current PCAT target pixel.
- The PCA basis never includes the target pixel being predicted. This avoids direct target-series leakage into the
  temporal modes, including inside the held-out gap.
- Detector-frame P4 is the initial supported coordinate frame. Rotated-frame PCAT has different sampling and mask
  semantics and is rejected until separately specified.

## Implementation seam

The new code should not overload `P4PCA`, whose current responsibility is spatial predictor regression. Introduce a
small numerical component, tentatively `P4TemporalPCA`, with an input-independent interface:

```text
configure(reference-series matrix, centering, requested mode)
predict(target-series, gap image count, requested mode)
    -> predictor column, diagnostic residuals, validity, numerical-rank diagnostics
```

`P4Reduction::fitDetectorSearch()` remains the authoritative P4 sampling/fit seam. `P4Reduction::regions()` owns
the PCAT geometry and annulus-wide configuration. For each annulus it will:

1. resolve and validate the target-specific reference-pixel set `R_p` for each search pixel;
2. sample `Z_p` and build/factor its temporal covariance once per search pixel;
3. sample the corresponding `y_p` series and calculate one gap-held-out PCAT predictor value for every target image;
4. append that predictor column after the same-image spatial P4 columns before calling the existing `P4PCA` path;
5. propagate PCAT invalidity to the ordinary P4 local-fit validity for every output plane; and
6. optionally retain compact prediction/residual diagnostics by search pixel and target image.

The initial CPU implementation may explicitly solve each `U[S_t,:]` least-squares problem. It must cache target-index
dependent restricted-mode factorizations because they are shared by every target image for one search pixel. Interior
targets have one common gap width, while only the two endpoint bands have smaller retained training dimensions.

No PCAT result may alter the ordinary P4 science product when `p4.temporalPredictor=none`, and no enabled PCAT result
may quietly alter P4's spatial mode counts, rank tolerance, spatial geometry, target-exclusion policy, or
final-combination method. It changes only the addition of one temporal predictor column, just as
`p4.numberImages>0` does today.

## Validation plan

### Numerical oracle

Implemented: known two-mode series reconstruct exactly across every gap; modifying all target values inside a gap does
not alter its prediction; and an oversized no-wrap gap marks its unsupported target invalid.

- Use synthetic reference pixels made from known temporal basis functions, with a target time series formed from the
  same modes. Verify exact gap-held-out reconstruction for every target time and endpoint.
- Add target components orthogonal to the retained modes and verify the expected residual rather than an artificial
  exact reconstruction.
- Compare the cached restricted-mode solve with an explicit least-squares solve independently for every target time,
  odd/even gap width, and every requested PCAT mode.
- Test mode-rank deficiency, insufficient outside-gap samples, non-finite input, common-mask holes, and invalid
  reference-region geometry. Each must produce the documented error or invalidity bit.

### Leakage and geometry tests

Implemented: reduction-level tests verify that PCAT appends one column and completes valid fits with both `or` and
`annulus` reference geometries. Remaining tests below are still useful for expanded coverage.

- Change `y_p[t]` and every sample inside `G_t` while retaining `y_p[S_t]`; `prediction_p(t)` must be unchanged.
- Change the target pixel's complete series while holding the reference region fixed; the temporal eigenmodes must be
  unchanged.
- Change a reference-region sample and confirm the modes and predictions can change.
- Verify that the `or` and `annulus` reference selections honor P4 signal exclusion, image edges, masks, and annulus
  ownership exactly.

### Scientific comparison

Run matched ordinary P4 reductions with `p4.temporalPredictor=none`, `neighborSummary`, and `pcat`, varying PCAT gap
widths, reference geometries, and mode fractions. Compare held-out residual variance and temporal autocorrelation at
empty locations; then compare injected-companion recovery and false-positive behavior. Do not claim an improvement
from in-sample residual reduction alone.

The earlier single-pixel mean/median temporal-summary result suggests radial-profile preprocessing has already removed
much of the per-pixel temporal mean. PCAT should therefore be judged on its ability to capture structured temporal
variation across a real withheld gap, not on whether it restores a mean level that preprocessing intentionally
removed.

## Deferred decisions

- Whether a later configuration should permit PCAT and a neighboring-frame summary as separate simultaneous predictor
  columns, rather than retaining the initial one-column mutually exclusive temporal-augmentation contract.
- A physical-time or parallactic-angle gap for irregularly sampled observations.
- Robust PCA/reference weighting, temporal noise weighting, and treatment of partially masked reference pixels.
- Frozen-PSF response, pixel-local fake injection, target-image exclusion, GPU acceleration, and a rotated-frame
  implementation. Each requires an explicit compatibility and validation plan before it is enabled.

## Reference

- Long et al. (2023), [PCAT](https://arxiv.org/abs/2303.05559).
