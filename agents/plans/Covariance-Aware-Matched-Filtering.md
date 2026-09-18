# Covariance Aware Matched Filtering

I think that the finite difference derivative of the KL basis that we are using in the PSF response estimation is reinventing the KLIP-FM formalism of https://arxiv.org/abs/1604.06097, and we are working towards re-implementing the matched-filtering follow-on work of https://arxiv.org/abs/1705.05477v1.

Please evaluate our current implementation of PSF response estimation (use the P4 version) and use of that response as a matched filter (in hciAnalyze) to compare and contrast to those papers. In particular, I'm interested if the derivations in Pueyo provide a way to speed up the calcualations.

The Ruffio matched filter result only uses a local standard deviation.  However I think we can use a full covariance.  See this discussion where we derived a PCA-based matched filter: [PCA Wiener filtering discussion](PCA_Wiener_filtering_literature.md).

Please document your findings below. Do not alter this prompt.

# Findings

## Summary

Reviewed 2026-09-14 against checkout `9bfe31e235d3da45eeebe0656acd7ca46f418d78`.

1. **Yes, an analytic derivative can replace the expensive paired P4 refits in the locally linear regime.**
   The current implementation differences complete local regression residuals, not individual KL eigenvectors.
   Its derivative includes changes in the basis and in the fitted regression coefficients. Both matter.
2. **The frozen P4 response is more restrictive than a frozen-basis projection.** It freezes the predictor
   coefficients as well. Differentiating only eigenvectors and adding that correction to the existing frozen
   response would therefore omit another term.
3. **Full residual-noise covariance weighting is a sound next step.** The difficult parts are estimating and
   regularizing that covariance, keeping the response and noise in the same coordinates, and calibrating detection
   significance. The PCA Gram matrix used to subtract speckles is not automatically this covariance.
4. **Response accuracy and noise weighting are separate improvements.** Implement and validate each against the
   present filter so that gains or regressions can be attributed to the right change.

The supplied [PCA/Wiener discussion](PCA_Wiener_filtering_literature.md) is preserved as a Markdown source.
Section 5 connects its derivation to the proposed implementation and distinguishes it from Wiener image estimation.

## 1. Connection to the papers

**Pueyo:** Appendix E, especially E14–E18, differentiates the covariance eigenpairs and propagates their changes
into the KL modes. E18 includes mode-normalization, mode-mixing, and source-in-reference contributions. Section II.4
and Appendix F turn this into a reusable first-order response, avoiding a new reduction for each trial flux or
spectrum. The expansion requires a suitable small-perturbation regime. This is the relevant analytic foundation
for replacing numerical refits, but P4 needs the derivative of its own regression operator.
[Pueyo 2016, v2](https://arxiv.org/html/1604.06097v2).

**Ruffio:** the premise needs one correction: Appendix A.2 already allows a general covariance matrix. A12 and A13
give covariance-weighted amplitude and S/N; A14 and A15 specialize to independent noise with a local variance per
exposure/wavelength. The implemented forward-model filter combines information from individual processed exposures
and wavelengths, using position-dependent responses. Section III.6 additionally calibrates the detection map
empirically. Thus implementing off-diagonal covariance here would extend the paper's practical noise model, within
its existing statistical framework.
[Ruffio et al. 2017, requested v1](https://arxiv.org/html/1705.05477v1).

The structural connection is a source-dependent PCA projection. The fitted objects differ: KLIP uses spatial
eigenimages of its reference library, while P4 predicts one detector pixel's time series from predictor-pixel time
series. This changes which eigensystem and coefficient dependence must be differentiated.

## 2. What the P4 implementation actually calculates

### Code map

Paths below are relative to the repository root; line numbers describe this review's checkout.

| Component | Production code | Relevant behavior |
| --- | --- | --- |
| Detector regression | `src/common/P4Reduction.cpp:2391`, `fitDetectorSearch()` | Samples the target and predictors, adds the trial source to both, and invokes direct PCA for in-sample fits. |
| Numerical regression | `src/common/P4PCA.cpp:472`, `P4PCAKernel::calculateValidated()` | Forms the smaller temporal/predictor Gram matrix, selects supported modes, and accumulates predictions. |
| Frozen response | `src/common/P4PSFModel.cpp:260`, `calculateLocalResponse()` | Applies the stored predictor coefficients to sampled source responses. |
| Paired refits | `src/common/P4Reduction.cpp:5678`, `calculateRefitDifferenceSamples()` | Performs positive and negative local detector fits and differences their residual time series. |
| Reconstruction | `src/common/P4Reduction.cpp:5988` and `src/common/P4PSFReconstructor.cpp:531` | Reconstructs source-centered sky stamps and combines valid frame contributions. |
| Sparse model | `src/common/RadialPSFModel.cpp:257`, `fit()` | Builds the angularly averaged radial response model used to populate the spatial field. |
| Filter | `src/common/P4PSFFilter.cpp:74`, `calculate()` | Returns signed correlation, template energy, amplitude, and support. |
| External analysis | `src/apps/hciAnalyze.cpp:719`, `filterCubePSFResponse()`; `:1502`, `analyzeCube()` | Loads P4 responses, filters final-image planes, then estimates the amplitude map's radial noise. |

For one detector search pixel, let `X` be the `T × P` predictor matrix and `y` its `T`-sample target time series.
The currently supported paired-refit path is **uncentered, in-sample detector regression**. It uses the direct
mixed-precision PCA path, not `calculateCenteredInPlace()`. With `X = U S Vᵀ` and `k` retained modes,

$$
\beta_k=V_k S_k^{-1}U_k^T y,\qquad
r_k=y-X\beta_k=(I-U_kU_k^T)y.
$$

For a unit-contrast trial source at sky location `q`, write its sampled predictor and target contributions as
`A_q` and `s_q`. Holding the fitted coefficients fixed gives

$$
h_{\mathrm{frozen},q}=s_q-A_q\beta_k.
$$

This describes the existing frozen-response contract. A refit instead differentiates

$$
r_k(\alpha)=y+\alpha s_q-(X+\alpha A_q)\beta_k(\alpha),
\qquad
r_k'(0)=s_q-A_q\beta_k-X\beta_k'(0).
$$

The missing `-X β′` includes the target's effect on the fitted coefficients, even when the basis does not move.
For example, if `A_q=0`, the frozen response is `s_q`, but the refitted response is
`(I-U_k U_kᵀ)s_q`. This is an especially useful regression test for an analytic implementation.

### Paired-refit scope and combination order

`calculateRefitDifferenceSamples()` evaluates `(r(+ε)-r(-ε))/(2ε)` with
`ε=psfResponse.refitContrast`. For each sampled source, it refits only the detector search pixels needed to
reconstruct the output stamp. Each fit still uses its complete regression time series. This is already much less
work than two full-field reductions per source.

The routine differences detector residuals **before** sky reconstruction and frame combination. Positive and
negative rank support must both be valid. The resulting stamps then enter angular averaging and region-aware
radial interpolation. Consequently the delivered template has three distinct approximations: finite perturbation
amplitude, spatial averaging/interpolation, and the response-combination rule.

At `processPSFProducts()` (`src/common/P4Reduction.cpp:6087`), science `sigmaMean` becomes an unclipped response
mean, preserving frame weights and minimum support. With fixed support, linear reconstruction and weighted means
commute with differencing. A median does not; independently clipping responses is also not the derivative of
science clipping. The current response is therefore not generally the derivative of every configurable final
science estimator. This distinction should stay explicit in provenance and validation.

Configuration currently restricts PSF responses to detector coordinates and excludes PCAT and post-median
subtraction. `refitDifference` additionally requires sparse radial sampling, `p4.numberImages=0`, and
`adi.excludeMethod=none`. These are appropriate initial boundaries for the analytic replacement too. The probe
uses a post-preprocessing PSF and unit per-frame scales; a source model with variable frame scaling must preserve
that normalization explicitly.

### Existing empirical evidence

The recorded [P4 response experiments](p4_psf_calculation_post_processing.md) support coefficient adaptation as
the main missing effect in the frozen response:

| Recorded experiment | Result |
| --- | --- |
| Frozen sparse response, 2026-09-07 | Fitted contrast was `0.2451` times the exact negative-companion result. |
| Dense response after source removal | Contrast ratio remained `0.2431`; source contamination and sparse sampling did not explain the discrepancy. |
| Mean versus sigma-clipped science | Finite-response projection ratios were `0.210678` and `0.210662`; clipping did not explain it either. |
| Paired refits after avoidance fix, 2026-09-12 | Contrast ratio `1.02539`, position difference `0.244` pixels, response cosine `0.97844`, elapsed time `6,574` seconds for `232` measurements. |

The earlier paired-refit run recorded `352,558` positive-plus-negative detector fits. These are prior measurements
reported in the plan, not experiments rerun for this review. The successful finite-amplitude result is evidence
that refitting matters; it does not establish infinitesimal convergence or performance on other datasets.

## 3. An analytic derivative suited to the current P4 path

The following is a derivation from the direct P4 regression above. It uses the same eigen-perturbation principle,
but differentiates the **retained temporal projector**, avoiding eigenvector sign alignment and unnecessary
within-subspace rotations.

Define

$$
G=XX^T,\quad \dot G=A_qX^T+XA_q^T,\quad
G u_i=\lambda_i u_i,\quad \Pi_k=\sum_{i=1}^{k}u_i u_i^T,
$$

with eigenvalues in descending order and a positive gap between retained and discarded groups. Then

$$
\dot\Pi_k=
\sum_{i\le k}\sum_{j>k}
\frac{u_j^T\dot G u_i}{\lambda_i-\lambda_j}
\left(u_j u_i^T+u_i u_j^T\right),
$$

and the complete local first-order response is

$$
\boxed{h_q=(I-\Pi_k)s_q-\dot\Pi_k y.}
$$

This is the replacement for the **whole paired residual difference**. It should not be added to the existing
frozen-coefficient template, which would mix two different decompositions and risk double counting.

Only retained–discarded interactions enter `Π̇`. Repeated eigenvalues wholly within either group are harmless to
the projector derivative. A repeated eigenvalue across the cutoff makes the chosen rank-`k` subspace ambiguous;
a small boundary gap makes the derivative sensitive. Rank-threshold crossings and discrete support changes also
need diagnostics and a finite-refit fallback. Do not silently clip eigenvalue-gap denominators and call the result
the same derivative.

The sum needs the discarded complement, including zero-eigenvalue directions when present. The current solver
returns `min(T,P)` eigenvectors for direct fits. When `P<T`, either account analytically for the remaining nullspace
or use a predictor-side/SVD formulation; simply summing the returned nonzero modes is incomplete.

### How to make it faster

1. **Reuse the unperturbed fit.** Retain an appropriate eigensystem and target projections while a detector pixel
   is active. Build the unit-source inputs with the existing trial-source interpolation conventions. Eliminate the
   two new perturbed Gram constructions and eigensolves per probe.
2. **Apply the derivative without constructing a dense `Π̇`.** If `U_a` contains retained modes, `U_b` discarded
   modes, and `B[j,i]=(u_jᵀ Ġ u_i)/(λ_i-λ_j)`, then

   $$
   \dot\Pi_k y=U_b B(U_a^Ty)+U_a B^T(U_b^Ty).
   $$

   Compute `Ġ U_a = A_q(Xᵀ U_a)+X(A_qᵀ U_a)` directly. For `T≤P`, its leading dense cost is
   `O(T P k)` rather than forming a new `O(T² P)` Gram matrix. The remaining basis products still cost work;
   this is not a constant-time response.
3. **Group overlapping probes by detector search pixel.** The existing source-first loop repeats the same base
   sampling and regression across overlapping stamp footprints. A detector-pixel-first pass can serve several
   source probes before releasing its workspace. Keep batches bounded and accumulate sky-stamp sums using the
   existing geometry; do not retain every pixel's full eigensystem or every frame-dependent response cube.
4. **Share work across requested mode counts.** Reuse the eigenbasis, source contractions, and target projections;
   each mode count defines its own retained/discarded split.

The likely saving is substantial when probes share fits, but no numerical speedup factor is established here.
Sampling, matrix products, reconstruction, and I/O remain. Benchmark the same `232`-sample configuration and report
wall time, peak RSS, unique detector fits, and response error against the paired-refit oracle.

### Extensions that require their own derivation

The centered rotated-frame implementation fits centered `X_c,y_c` but applies the coefficients to the uncentered
predictors. Its residual is `y-X β(X_c,y_c)`. Differentiating only a centered projector would omit the uncentered
application/mean terms. Target-held-out fits likewise require the derivative of each training-set coefficient fit
and its application to the held-out row. PCAT adds another data-dependent prediction stage. None is covered by the
simple direct-path formula without further work.

A derivative about data containing a bright real source is a local tangent there, not necessarily the finite
source contribution relative to source-free data. Preserve the negative-companion oracle; test amplitude sweeps
and, where necessary, relinearize about a source-subtracted estimate.

## 4. What `hciAnalyze` measures now

On valid stamp pixels, `P4PSFFilter::calculate()` forms

$$
N=t^T d,\qquad E=t^T t,\qquad \widehat\alpha=N/E.
$$

It retains signed negative response lobes, rejects insufficient support, and ignores nonfinite/out-of-bounds
science samples. This is the least-squares amplitude for uniform, uncorrelated noise. A constant local scalar
variance would cancel from this amplitude, but is needed for its uncertainty. `E` is template energy, not a noise
variance or a complete inverse amplitude variance.

The external P4 path loads coordinate-indexed final templates and applies this calculation independently to each
already-combined science mode plane. It does not combine per-exposure or per-wavelength likelihood contributions.
Ordinary image combination is generally not a sufficient statistic when frame responses and noise weights differ;
later covariance weighting of the final image cannot recover information discarded during that combination.

After filtering, `analyzeCube()` can apply configured Gaussian high/low-pass filtering to the amplitude map. It
then excludes configured source regions, calls `stddevImageCube(..., true)` to subtract the radial mean and divide
by the radial standard deviation, and applies `correctSmallSampleSNR()`. The correction uses
`n=2πr/(lambda/D)-1` and multiplies by `1/sqrt(1+1/n)`. The reported aperture value is a **maximum over pixels**,
not an aperture sum or a joint likelihood fit.

This is a useful empirically normalized detection map. It does not whiten correlated stamp pixels or supply an
exact false-alarm probability, especially after position/mode selection. The curvature covariance in
`agents/plans/scripts/fit_p4_matched_response.py` is a position-fit diagnostic derived from a local
`0.5*SNR²` surface; it is not a residual-pixel covariance already available to the filter.

## 5. A covariance-aware filter and its PCA implementation

### Connection to the supplied PCA/Wiener discussion

The conversation's matched-filter derivation is correct for an orthonormal eigenbasis of the **noise covariance
in the space being filtered**, with positive eigenvalues. Write `C=Q diag(ν_i) Qᵀ`, where the columns `q_i` are
noise eigenvectors, `ν_i` are noise variances, `t_i=q_iᵀt`, and `z_i=q_iᵀ(d-μ)`. Then

$$
\widehat\alpha=\frac{\sum_i t_i z_i/\nu_i}{\sum_i t_i^2/\nu_i},\qquad
S=\frac{\sum_i t_i z_i/\nu_i}{\sqrt{\sum_i t_i^2/\nu_i}}.
$$

This retains the coherent signed template projections and weights each mode by inverse noise variance. It does
not need a signal PSD. Removing the square root from the second denominator produces the amplitude estimator,
not normalized S/N.

Two qualifications matter when interpreting the statement that the PCA already supplies these quantities:

- Eigenvectors of a small image-to-image Gram matrix are not themselves spatial eigenimages. For a mean-centered
  `n × p` matrix of training stamps `R`, a positive eigenpair `RRᵀ v_i=γ_i v_i` gives spatial mode
  `q_i=Rᵀv_i/sqrt(γ_i)` and sample-covariance eigenvalue `ν_i=γ_i/(n-1)`. Preserve the actual centering and
  normalization convention. A global covariance scale cancels from amplitude, but changes S/N and uncertainty.
- The P4 eigensystem in Section 3 describes its local predictor regression. It does not directly provide the
  covariance eigenmodes of the final residual stamp. A residual-noise training step or explicit covariance
  propagation is still needed.

For Wiener **image estimation**, a zero-mean random signal uncorrelated with noise has linear minimum-MSE
estimator `s_hat=C_s(C_s+C_n)⁻¹d`. The scalar modal gains in the conversation are exact when both covariances
are diagonal in the chosen basis. A noise PCA basis alone generally does not diagonalize `C_s`.
For a fixed planet template and a zero-mean amplitude prior of variance `v_α`, `C_s=v_α t tᵀ`; its representation
in the noise basis includes cross-mode products `v_α t_i t_j`. Keeping only `t_i²` therefore changes the model.

For that rank-one model, the Wiener amplitude is `N_C/(E_C+1/v_α)`, using the quantities below. It approaches the
unregularized matched-filter amplitude as the prior variance grows. Thus the two ideas connect precisely, but
independent Wiener gains on each PCA component are not the same estimator as coherent template detection.

The discussion's PACO connection is also relevant: PACO learns local background means and spatial patch
covariances for detection; PACO ASDI extends the statistical modeling to spatio-spectral data. These are useful
precedents for covariance estimation and calibration, rather than evidence that the existing P4 Gram matrix is
already sufficient. See [Flasseur et al. 2018](https://olivier-flasseur.github.io/publication/2018-aa-paco/) and
[Flasseur et al. 2020](https://olivier-flasseur.github.io/publication/2020-aa-paco-asdi/).

### Statistical model and output quantities

For one candidate, use the same valid-pixel ordering for a data stamp `d`, a unit-contrast processed template `t`,
and a residual-noise covariance `C`. Subtract an estimated background mean `μ`, giving `z=d-μ`. Assume
`z=αt+n`, with fixed positive-definite `C=Cov(n)` independent of the trial amplitude. Minimizing
`(z-αt)ᵀ C⁻¹(z-αt)` gives

$$
N_C=t^T C^{-1}z,\qquad E_C=t^T C^{-1}t,
$$

$$
\boxed{\widehat\alpha=N_C/E_C,\quad
\sigma_\alpha=E_C^{-1/2},\quad
S=N_C/\sqrt{E_C}.}
$$

Here `σ_α` is conditional on the supplied covariance and template being correct. Covariance-estimation error and
template error require separate validation. Solve `Cw=t` and evaluate `wᵀz,wᵀt`; do not explicitly invert `C`.
Equivalently, whiten **both** `z` and `t` by a Cholesky factor or the covariance eigensystem, then use dot products.

For a signed unconstrained amplitude, the log-likelihood improvement over zero signal is `S²/2`. For a physical
nonnegative source, it is `max(0,S)²/2`; a large negative filter response should not become a positive detection
merely because it was squared. Fitting a local constant or gradient background jointly with the source is possible,
but requires projecting those nuisance components from both the data and template in the covariance metric.

### PCA as a covariance model, rather than another hard subtraction

Learn orthonormal noise modes `U_r` from comparable, source-free **residual** stamps. A useful regularized model is

$$
C=\tau^2I+U_r\,\operatorname{diag}(q_i)\,U_r^T,
\qquad \tau^2>0,\quad q_i\ge0.
$$

The `q_i` are excess variances above the positive noise floor. This gives

$$
C^{-1}v=\tau^{-2}\left[v-U_r\,
\operatorname{diag}\!\left(\frac{q_i}{\tau^2+q_i}\right)U_r^Tv\right].
$$

With `p` stamp pixels, applying this operator costs `O(pr)` after training. A diagonal floor `D` permits
spatially varying independent noise:

$$
C=D+LL^T,\qquad
C^{-1}=D^{-1}-D^{-1}L(I+L^TD^{-1}L)^{-1}L^TD^{-1},
$$

where `L=U_r diag(sqrt(q_i))`; only a small `r × r` system needs factorization in addition to diagonal scaling.
Shrinkage toward a diagonal covariance is another useful baseline.

In a complete covariance eigenbasis, whitening divides each mode coefficient by the square root of its variance.
Hard PCA subtraction instead discards selected modes. It corresponds to a limiting nuisance model with infinite
variance in those directions, and sacrifices any source information there. Truncating a sample covariance and
pseudoinverting it also discards the unmodeled complement; the positive floor above preserves that complement with
finite uncertainty. These choices should be explicit, rather than hidden behind a single PCA mode-count setting.

## 6. Estimating the right covariance

**Start with spatial covariance of final-image stamps.** This fits the existing `hciAnalyze` input contract and
provides a controlled comparison. Estimate noise from blank locations in comparable radius/region ranges,
excluding known sources and the tested location with a guard region. Tune rank and regularization on held-out
locations or blocks, then measure performance on separate injections/null samples.

### Overlapping annular patches: a Welch-style sampling scheme

The agreed starting geometry is a sequence of overlapping patches around an annulus, analogous to Welch's
overlapping-segment spectral estimator. Under stationarity, autocovariance and the power spectral density are
Fourier-transform pairs. Welch averages windowed segment powers; here we can average patch outer products while
retaining the full pixel covariance. The latter also retains correlations between Fourier modes and therefore
estimates more structure than a diagonal spectral-power model.

Let $\ell=\lambda/D$, let $a$ be the patch radius, and let $R_{\mathrm{ann}}$ be the radius of its center from the
star, all in image pixels. Half-width spacing means an azimuthal center-to-center arc step $\Delta s=a$ for a patch
of diameter $2a$. The approximate number of patches around a complete annulus is

$$
n\simeq\frac{2\pi R_{\mathrm{ann}}}{\Delta s}
=\frac{2\pi R_{\mathrm{ann}}}{a}.
$$

Thus a patch radius of one $\lambda/D$ at a separation of three $\lambda/D$ gives $n\simeq6\pi\simeq19$ samples
before masks and source exclusions. This describes 50% overlap along the azimuthal width, not 50% overlap in the
area of circular patches. Whether that patch radius captures the processed source's negative lobes remains a
response-support measurement; a larger required footprint reduces the count at the same separation.

Rotate each patch into the same radial/tangential coordinate system. For the resulting vectorized patches $x_j$,
use their mean $\bar x$ and the centered covariance estimate

$$
\widehat C=\frac{1}{n-1}\sum_{j=1}^{n}(x_j-\bar x)(x_j-\bar x)^T.
$$

Overlap improves coverage and can improve averaging. Its benefit depends on the correlations between the
quadratic estimates and any windowing; do not automatically replace 19 overlapping patches with either 19
independent samples or 9 useful samples. The usual $n-1$ normalization does not itself correct correlations
between training patches. Validate covariance scale and effective precision on held-out angular neighborhoods.
Welch's own variance/resolution tradeoff depends on overlap and tapering; see the
[MathWorks discussion of Welch estimation](https://www.mathworks.com/help/signal/ug/nonparametric-methods.html).

The raw centered covariance from 19 patches has rank at most 18, regardless of the number of stamp pixels.
Use the positive-floor PCA model or diagonal shrinkage from Section 5 so the unmeasured complement retains finite
uncertainty. Eigendecomposition of this covariance is a representation of the same matched filter; truncation or
regularization changes the covariance model. Fifty-percent spacing is an initial geometry to test, rather than a
claim of optimality for full covariance estimation. Choose tapering separately and propagate any measurement
taper through the data, template, and covariance consistently. A window used only inside a spectral estimator
is a different choice: the later PSD prototype reconstructs a covariance for untapered stamps, leaves the data
and template untapered, and explicitly retains the estimator's window-induced lag bias.

For candidate exclusion and validation, hold out whole angular neighborhoods, including every training patch
whose footprint intersects the candidate's exclusion region. Leaving out only the patch centered on the candidate
would still expose the covariance estimate to that source through overlapping neighbors. Compare the resulting
filter with identity weighting using held-out null scores, injected-source completeness, and amplitude uncertainty
calibration, as planned in Section 7.

### Pooling patches across radii and testing radial variance normalization

The next development comparison (agreed 2026-09-18) should explicitly pool training patches whose centers lie
at different stellar radii. The user's motivation is that the P4 optimization/predictor region (OR) is much wider
than the search region (SR), and making the OR smaller worsens the reduction. That motivates testing radial
correlations and shared residual structure over a wider range. It does not alone establish that final-image noise
covariances are interchangeable across that range; covariance shape and variance scale must both be measured.

Keep the current 11×11 patch and five-pixel angular arc spacing. Around candidate radius $R_{\mathrm{ann}}$, use
training-center radii from $R_{\mathrm{ann}}-b$ to $R_{\mathrm{ann}}+b$, with spacing $\Delta R=5$ pixels and initial
radial half-widths $b\in\{0,5,10,20\}$ pixels. Skip nonpositive center radii and retain the existing complete-support/exact-stencil
checks. Rotate each patch into the candidate's orientation without magnifying it radially. This changes the
locations supplying covariance samples, not the response footprint or the P4 reduction's OR/SR settings.

Use a four-way comparison to separate sampling from normalization:

| Training centers | Raw residual pixels | Radial-variance-normalized pixels |
| --- | --- | --- |
| Candidate radius only | Existing sampler/control | Normalization-only control |
| Radial band | Pooling-only change | Pooling plus normalization |

For the normalization arm, estimate a positive radial variance profile $\widehat v_{\mathrm{rad}}(\rho)$ from
unfiltered native image pixels, with known sources and all declared development/calibration/evaluation footprints
excluded. Divide by the **standard deviation**, $\widehat\sigma_{\mathrm{rad}}(\rho)=\sqrt{\widehat v_{\mathrm{rad}}(\rho)}$,
before rotating/interpolating training patches. Pixelwise normalization
handles a gradient across a stamp; dividing each entire stamp by its center's scale is a distinct approximation.
The profile estimator, binning/smoothing, minimum support, and endpoint behavior belong to the frozen policy.
Do not infer its scale from a held-out candidate or extrapolate through unsupported radii.

For candidate $q$, let $D_{\sigma,q}$ be the diagonal matrix of those standard deviations at its native stamp
pixels. If $\widehat C_{\mathrm{std}}$ and $\widehat\mu_{\mathrm{std}}$ are the regularized covariance and mean
estimated from standardized training patches, filter with

$$
t_{\mathrm{std}}=D_{\sigma,q}^{-1}t,\qquad
z_{\mathrm{std}}=D_{\sigma,q}^{-1}d-\widehat\mu_{\mathrm{std}}.
$$

Equivalently, in original image units use

$$
\widehat C_q=D_{\sigma,q}\widehat C_{\mathrm{std}}D_{\sigma,q},\qquad
\widehat\mu_q=D_{\sigma,q}\widehat\mu_{\mathrm{std}}.
$$

Transforming the response as well as the data preserves the contrast parameter. The standardized covariance
need not have a unit diagonal: its training samples undergo interpolation and overlap, and local covariance
shape may vary even after removal of a radial scale. Profile uncertainty also belongs to the eventual empirical
uncertainty-calibration problem.

First hold rank three and floor fraction 0.1 fixed to isolate the sampling/normalization changes. On the existing
development data, compare accepted patches by radius, covariance spectra and shape across bands, split-block
stability, held-out residual prediction/conditional coverage, and source leakage using the saved full-image
injections. Evaluate the profile and covariance with the same footprint exclusions. Then investigate rank/floor
separately rather than choosing every parameter from the same recovery curve. Equal weight per retained patch
is the initial convention; outer rings offer more centers, so record their contribution explicitly. More overlapping
patches do not imply an equal increase in independent information.

The [first radial-pooling audit](results/p4-step5-radial-pooling-20260918/README.md) checks geometry and variance
scale only. It confirms a useful increase in available patches, a strong inner variance gradient, and a large
outermost-bin variance rise. Assess processing-region boundaries and covariance shape before choosing a band;
neither the count increase nor variance normalization alone demonstrates a detection gain. Retain Gaussian and
identity references when the new filters are evaluated.

### Covariance constraints and later extensions

Key constraints:

- **Different matrices:** `XXᵀ` describes predictor time-series correlations used by P4. A final-stamp covariance
  describes correlations among output residual pixels after regression, interpolation, rotation, and combination.
  Their eigenvectors have different coordinates and generally different dimensions. They are not interchangeable.
- **Finite samples:** a covariance from `n` mean-centered training stamps has rank at most `n-1`. Overlapping
  stamps and neighboring frames further reduce effective independence. A nominally invertible sample covariance
  can still give unstable inverse weights. Use shrinkage/floors and hold-out validation.
- **Geometry:** if stamps are rotated into a radial/tangential frame for training, transform the data and template
  consistently too. Validate approximate azimuthal stationarity and handle region boundaries explicitly. Radial
  averaging of the template does not demonstrate stationarity of the noise.
- **Missing pixels:** restrict `C`, `t`, and `z` to the actual valid support before solving. The inverse of a
  restricted covariance is generally not the corresponding block of the full inverse. Zero-filling missing data
  and reusing unchanged weights produces a different statistic.
- **Data reuse:** training covariance on the source stamp can downweight the source itself. Training on noise
  processed differently from the science stamp produces the wrong covariance. Track the reduction, coordinates,
  mode, support, training mask, rank, regularization, and source exclusions in provenance.
- **Detection calibration:** retain signed amplitude, conditional uncertainty, raw score, and empirical
  calibration as distinct products. Use held-out nulls to check tails and injected sources to check completeness
  at a fixed false-positive rate. The existing small-sample multiplier alone does not account for covariance
  estimation, correlated trials, or searching over position and mode count.

For a fixed linear reduction `L_red`, propagation would give `C_out=L_red C_in L_redᵀ`. For adaptive P4, the
corresponding first-order object is its full data Jacobian `J`, giving `J C_in Jᵀ`. A few source-direction
responses `J s_q` do not determine that full Jacobian. Empirical residual covariance is therefore the practical
first route; reusing a speckle-fit Gram matrix would require an additional statistical derivation.

Later, retain time-resolved residual stamps and responses and combine likelihood contributions before image
collapse. Independent frame blocks allow sums of each block's `N_C` and `E_C`; correlated time blocks require a
temporal/block covariance model. The final cube's different PCA-mode planes are alternative reductions of the
same data, so they must not be counted as independent exposures.

## 7. Recommended implementation and validation sequence

This sequence was accepted on 2026-09-16. Begin with the derivative oracle; the overlapping annular-patch estimator
above supplies the initial sampling geometry for the later covariance-filter prototype.

1. **Establish the derivative oracle.** Sweep positive/negative amplitudes at representative radii, angles, mode
   counts, and eigengaps. Compare detector residuals before reconstruction, then final stamps. Separate finite
   amplitude error from mixed-precision cancellation, spatial interpolation, and frame-combination effects.
2. **Add the analytic direct-P4 response.** Implement the projector derivative for the current paired-refit scope,
   preserving `refitDifference` as the oracle/fallback. Test the target-only perturbation, predictor-only
   perturbation, both Gram orientations, rank deficiency, repeated internal eigenvalues, and cutoff crossings.
   Validate against the FP64 mathematical oracle first and the production mixed-precision path separately.
3. **Measure scientific and computational equivalence.** Reproduce the accepted post-avoidance AF Lep result,
   then broaden injections in position and brightness. Record contrast bias, astrometric error, template mismatch,
   amplitude dependence, wall time, and memory. Enabling response calculation must leave science products unchanged.
4. **Add noise weighting independently.** Extend the shared filter with identity, diagonal, and regularized PCA
   covariance choices. Preserve the current identity result. Verify normalization under template rescaling and
   compare low-rank solves with dense positive-definite solves, including changing valid support.
5. **Evaluate covariance gains on held-out data.** Compare detection completeness at a fixed false-positive rate
   and photometric uncertainty calibration. A larger value at the known planet alone is insufficient evidence.
   Include Gaussian smoothing (`filter.lpfGaussFW=3.6` for this dataset) as a reference before claiming a gain
   from either the response template or covariance weighting; distinguish filter shape from noise normalization.
6. **Consider time-resolved likelihoods after the final-image prototype.** Measure the information gained against
   storage and covariance-training costs before expanding the output format.

Current checkpoint (2026-09-18): the [initial fixed-policy Step-5 evaluation](results/p4-step5-evaluation-20260918/README.md)
is complete. It establishes neither a covariance detection gain nor calibrated conditional uncertainties;
the report records the remaining calibration and independent-validation work.

A useful response-approximation metric in the adopted noise model is

$$
\eta=\frac{\widetilde t^T C^{-1}t}
{\sqrt{(\widetilde t^T C^{-1}\widetilde t)(t^T C^{-1}t)}}.
$$

For a fixed correct covariance, `η` is the expected signed S/N relative to the exact-template optimum. Use it to
assess sparse/angular averaging alongside photometric scale error; ordinary image-space cosine alone can hide
mismatch in the most informative noise directions.

### Verification performed for this review

A standalone NumPy calculation checked the direct-P4 projector derivative against central differences at
`ε=1e-5`, using `(T,P)=(19,31)` and `(31,19)` and retained counts `1,5,12`. Maximum relative error was
`1.44e-9`. A case with repeated eigenvalues wholly inside retained/discarded groups gave `5.52e-10`.
The target-only perturbation limit was also checked. A regularized PCA inverse application agreed with a dense
solve to `8.80e-16` relative error, and its amplitude/S/N agreed with Cholesky-whitened dot products.

These initial-review checks verify the proposed algebra, not a production implementation or speedup. Existing
production tests were inspected, including the paired-local-reduction comparison in `tests/common/P4Reduction_test.cpp:2549` and
external P4 filtering in `tests/apps/hciAnalyze_test.cpp:349`. At that review stage no production source was changed
and no production test suite or ROC dataset was rerun. Only this findings document was edited; the mxlib
function-edit coverage gate was therefore not triggered. Remaining follow-ups were the analytic-response benchmark
and covariance estimation/calibration experiments described above. The supplied discussion is preserved without changes to its
content in [PCA_Wiener_filtering_literature.md](PCA_Wiener_filtering_literature.md).

### Implementation started: derivative-oracle groundwork (2026-09-16)

Step 1 now has a maintained numerical oracle in
[`tests/common/P4PCAResponse_test.cpp`](../../tests/common/P4PCAResponse_test.cpp). Seven Catch2 cases exercise the
actual FP64 `P4PCA::calculate()` and production mixed-precision `p4PCACalculateMixed()` entry points. The independent
reference constructs a known complete temporal eigensystem and evaluates the projector derivative from Section 3;
it does not reuse the production eigensolver or assume a thin SVD contains the discarded nullspace.

- Both Gram orientations, predictor-only and combined source perturbations, deficient baseline rank, and repeated
  eigenvalues within retained/discarded groups show second-order central-difference convergence. Halving the
  amplitude reduces the FP64 error by approximately four, with relative error below `5e-8` at `epsilon=1e-5`.
- The target-only case verifies coefficient adaptation: the derivative is `(I-Pi)s`, even when predictors are fixed.
- A separated but small cutoff eigengap requires smaller steps. A retained/discarded eigenvalue crossing instead
  produces differences growing as `1/epsilon`; rank-threshold crossings can make only one sign rank-supported.
- The mixed-precision sweep finds a usable finite-amplitude window in each tested case. Among the best sampled
  amplitudes for each case, the largest relative error was `0.143%` (the regression bound is `0.5%`). A diagonal example demonstrates complete
  cancellation at `epsilon=1e-10` in FP32 ingress while FP64 still resolves the response. These are fixture-specific
  results, not a universal contrast tolerance.

The final-stamp experiment driver is
[`run_p4_response_convergence.py`](scripts/run_p4_response_convergence.py). It requires preprocessed input images,
runs a zero-amplitude baseline and signed local trials for every requested radius/angle, and fixes direct
detector-frame P4 with mean combination and no temporal augmentation/exclusion. Each fresh output directory
preserves input hashes, commands, logs, timing, residuals, validity maps, and CSV/JSON comparisons. The driver checks
mode fractions, signed trial metadata, stamp origins, and valid support. It refuses existing output directories.

For adjacent amplitudes, it reports `||h_large-h_small||/||h_small||`, cosine, projection scale, and the normalized
even component `(r(+epsilon)+r(-epsilon)-2r(0))/(2epsilon)`. Comparisons use common valid pixels; support changes
are reported separately. Zero response norms produce missing diagnostics rather than an apparent exact match.
Neither member of an adjacent pair is designated the true derivative. Validity here describes final pixels;
detector-level eigengaps and changes in contributing temporal samples still need separate instrumentation.

#### Small AF Lep smoke experiment

Used the first 24 of the 621 local NACO `coadd5` frames, with radial-profile preprocessing performed once. Geometry
was the local profile configuration with search radii `[2,16)`, its existing predictor-wedge settings, and no
temporal augmentation. Trial radii were `7.8,10.8` pixels, position angles `0,90` degrees, and mode fractions
`0.05,0.15,0.3` (realized counts `1,3,7`). Stamps were `5x5`; the stored
`psf_reg_median.fits` template retained its original normalization. `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=2`.

The eight positive half-amplitudes below required 68 local reductions including the four baselines, producing
84 adjacent-amplitude comparisons. Every comparison retained all 25 pixels; neither final support nor signed
support changed. FITS diagnostics reported minimum rank 24 and zero rank-invalid counts. Summed subprocess wall
time was about 21 seconds; this small-run timing is not a production speedup measurement.

| Larger / smaller half-amplitude | Relative response change across the 12 position/mode cases |
| --- | --- |
| `0.03 / 0.01` | `0.0141%`–`0.2135%` |
| `0.01 / 0.003` | `0.0055%`–`0.4796%` |
| `0.003 / 0.001` | `0.0032%`–`0.7769%` |
| `0.001 / 0.0003` | `0.0100%`–`4.6201%` |
| `0.0003 / 0.0001` | `0.0397%`–`5.8602%` |
| `0.0001 / 0.00003` | `0.1350%`–`21.8909%` |
| `0.00003 / 0.00001` | `0.4338%`–`79.0559%` |

This is evidence against selecting the smallest amplitude automatically. The deterioration is consistent with
mixed-precision cancellation, but attributing the final-stamp error requires a higher-precision comparison and
detector-level diagnostics. Good agreement between two finite amplitudes alone does not bound their common bias.
The short sequence and small stamps do not replace the accepted full AF Lep validation or establish a universal
probe amplitude. The original smoke artifacts are under `/tmp/p4-response-convergence-20260916/broad_sweep`.

Example invocation after preparing `preprocessed/` and a geometry/angle configuration with absolute file paths:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 \
python3 agents/plans/scripts/run_p4_response_convergence.py \
  --binary _build_fresh/src/p4Reduce --config /path/to/smoke.conf \
  --inputs /path/to/preprocessed --psf /path/to/psf_reg_median.fits \
  --output /path/to/new-experiment --radii 7.8 10.8 --angles 0 90 \
  --modes 0.05 0.15 0.3 --stamp-size 5 \
  --amplitudes 0.03 0.01 0.003 0.001 0.0003 0.0001 0.00003 0.00001
```

Verification: all seven new cases passed (971 assertions), all 32 existing P4PCA cases passed (2204 assertions),
and the existing paired-refit/local-reduction equivalence case passed (40 assertions). The new source was
formatted with `clang-format`. No production implementation was changed, and the new functions call no upstream
mxlib APIs directly, so no mxlib coverage follow-up was introduced.

**Next in Step 1:** capture real detector regressions and cutoff gaps at the same trial locations, compare FP64
and mixed-precision amplitude sweeps before reconstruction, then propagate the responses through reconstruction
and mean combination. Establish the usable amplitude windows on longer sequences before treating paired refits
as a reference for the analytic production implementation in Step 2. The covariance-filter stages remain as
outlined above.

### Real detector precision comparison (2026-09-17)

The experimental precision build now offers a read-only `P4DetectorFitObserver` through
`p4ReductionReduceExperimental()`. It observes the actual sampled predictor matrix, target, retained counts,
rank threshold, and kernel residuals before conversion to image storage. The calling scope passes the observer
to each worker and restores the previous scope on both success and failure. This hook is compiled only with
`HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION`; it adds no observer dispatch to the default build.

`p4ReductionPrecisionBenchmark` adds `--capture-detectors` and `--capture-inputs`. Capture metadata identifies the
detector coordinate, dimensions, mode counts, numerical rank, supported modes, byte order, and array layout.
The companion binary contains column-major FP64 predictors and target when requested, followed by FP64 residual
storage. The latter holds promoted FP32 values for the mixed policy. Capture directories must be new, and I/O
failures propagate out of the reduction. Capture is restricted to local direct P4 without temporal augmentation
or target exclusion.

The sweep driver accepts `--precision P4-D64` or `P4-M32D64` with the benchmark binary. With
`--capture-detectors`, it records baseline and unit-source inputs and every signed trial's detector residuals.
[`analyze_p4_detector_convergence.py`](scripts/analyze_p4_detector_convergence.py) verifies matching experiment
inputs, exact equality of baseline/unit sampled matrices between policies, complete captures, geometry, and mode
counts. It calculates an independent SVD projector derivative using `A=X(1)-X(0)` and `s=y(1)-y(0)`. The complete
temporal complement is retained, including nullspace directions when `T>P`. Unresolved cutoff gaps and invalid
modes are reported separately; they are not treated as zero responses.

The reference is the derivative of continuous FP64 regression in that measured source direction. Unit-source
sampling is still FP32, so it is not an assertion that the entire floating-point pipeline is differentiable.
The analysis writes per-fit cutoff gaps, baseline agreement with SVD, per-fit derivative errors, aggregated
detector errors, and same-amplitude final-stamp differences between policies. Detector norms include all temporal
rows of the fits needed by the local stamp. Final-stamp differences use common valid pixels.

#### Results on 24 and 96 frames

Repeated the four positions and eight amplitudes above on both the original 24-frame subset and the first 96
frames, with preprocessing performed once per subset and identical input files supplied to both policies. The
96-frame realized mode counts were `4,14,28`. There were 305 captured detector fits across the four 24-frame
positions and 351 across the 96-frame positions. Every tested mode remained rank-supported and every final
comparison retained all 25 pixels. No cutoff gap was numerically unresolved by the analysis criterion.

Baseline FP64 residuals agreed with independent SVD residuals to `9.23e-14` (24 frames) and `8.57e-13` (96 frames),
normalized by the target norm. The smallest cutoff gap divided by the leading eigenvalue was `9.11e-5` and
`5.67e-7`, respectively; divided by the last retained eigenvalue, the minima were `0.0147` and `0.000817`.

Selected results below are maximum relative errors across the 12 position/mode cases, expressed as percentages.
Each detector error compares with the independent derivative. The final column compares mixed and FP64-PCA
final-stamp differences at the same amplitude; it does **not** treat the latter as an exact final-image derivative.
Column maxima need not occur at the same position/mode.

| Frames | Half-amplitude | FP64 detector error | Mixed detector error | Final mixed/FP64 difference |
| --- | --- | --- | --- | --- |
| 24 | `0.03` | `6.90%` | `6.90%` | `0.0275%` |
| 24 | `0.001` | `0.00998%` | `3.02%` | `0.942%` |
| 24 | `0.0001` | `0.000107%` | `22.5%` | `8.27%` |
| 24 | `0.00001` | `0.0000623%` | `193%` | `78.4%` |
| 96 | `0.03` | `26.2%` | `25.9%` | `0.228%` |
| 96 | `0.001` | `0.0759%` | `20.5%` | `11.4%` |
| 96 | `0.0001` | `0.000755%` | `105%` | `26.5%` |
| 96 | `0.00001` | `0.0000623%` | `716%` | `313%` |

The detector comparison separates finite-amplitude bias at large steps from mixed-precision error at small
steps. FP64 reaches a roughly `6e-7` relative-error floor in these sampled directions. Agreement between two
precision policies at a large amplitude does not establish negligible finite-amplitude bias.

**The final-image difference has an additional precision limit.** `P4-D64` selects FP64 PCA, but residual image
storage, spatial reconstruction, and frame combination remain FP32. At the smallest adjacent pair
(`3e-5` versus `1e-5`), the FP64-PCA final responses still differ by up to `17.0%` in the 24-frame subset and
`1.07%` in the 96-frame subset, even though pre-storage detector responses have converged. Thus subtracting two
stored final images is not an adequate tight-tolerance derivative oracle. Form the detector response before
image conversion, then reconstruct and combine that response. The existing `refitDifference` response path
already follows this ordering; the local-stamp sweep exposes why it matters.

Original artifacts are under `/tmp/p4-response-convergence-20260917`, with `24_D64`, `24_M32D64`, `96_D64`,
`96_M32D64`, `analysis24`, and `analysis96` subdirectories. Each precision sweep comprises 72 reductions: four
baselines, four unit-source captures, and 64 signed trials. These diagnostic timings include capture I/O and
are not response-speed benchmarks.

Build the optional benchmark with:

```sh
cmake -S . -B _build_response_precision -DHCIREDUCE_BUILD_TESTS=ON \
  -DHCIREDUCE_BUILD_CPU_BENCHMARKS=ON -DHCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION=ON
cmake --build _build_response_precision --target p4ReductionPrecisionBenchmark -j2
```

Use the earlier sweep invocation with `--binary _build_response_precision/benchmarks/p4ReductionPrecisionBenchmark`,
`--precision P4-D64 --capture-detectors`, then repeat with `P4-M32D64` and a different output directory. Analyze with:

```sh
OPENBLAS_NUM_THREADS=1 python3 agents/plans/scripts/analyze_p4_detector_convergence.py \
  --double /path/to/D64-sweep --mixed /path/to/M32D64-sweep --output /path/to/new-analysis
```

Verification: the experimental P4Reduction suite passed 33 cases / 349942 assertions, including exact replay of
captured FP64 and mixed regressions, unchanged science products, parallel capture, and observer restoration on
success/failure. The default-build P4Reduction suite passed 27 cases / 77657 assertions, and the seven response
oracle cases passed 971 assertions in the experimental build. All 68 comparable mixed-policy final images from
the 24-frame capture run matched the earlier production run bit for bit. Modified C++ sections were formatted.
The mxlib audit found no recorded instantiation of `eigenCube<float>::image(Index) const`; the concrete upstream
coverage follow-up is recorded in [mxlib_cleanup.md](mxlib_cleanup.md) under Known non-blocking ownership follow-ups.

### Reconstruction oracle validated (2026-09-17)

The existing `refitDifference` path was exercised with FP64 paired detector fits, subtraction before the float
image-storage boundary, local spatial reconstruction, mean combination, and the published radial PSF model.
No production behavior or default precision policy changed.

The maintained runner is [`run_p4_refit_convergence.py`](scripts/run_p4_refit_convergence.py). It reuses and verifies
the input, PSF, geometry, and precision provenance of a completed local sweep. A new output directory records the
commands, logs, timings, FITS products, and adjacent-amplitude comparisons. A response-disabled baseline checks
that enabling the products leaves the final science image exactly unchanged. An explicit `--binary` permits a
rebuilt experimental executable and records its new hash; other experiment inputs must retain their hashes.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 python3 agents/plans/scripts/run_p4_refit_convergence.py \
  --experiment /path/to/D64-sweep --output /path/to/new-refit-sweep --samples-per-radius 4
```

Both the 24- and 96-frame subsets used the same eight amplitudes and three mode fractions as the preceding
detector experiment. Each run measured four angular positions at each of the two requested radii (`7.8,10.8`),
then published `5x5` radial-model stamps at all 800 search coordinates. These come from eight distinct measurement
locations before angular averaging, not 800 independent response measurements; their noise is not assumed
independent. Selection uses available integer detector coordinates, so these measurements are distinct
from the preceding fractional-position local trials.

The table reports the maximum `||h_large-h_small||/||h_small||` over all published source/mode stamps. Values are
dimensionless fractions. Every comparison retained all 25 pixels per stamp, with no support changes.

| Larger / smaller half-amplitude | FP64, 24 frames | FP64, 96 frames | M32D64, 24 frames | M32D64, 96 frames |
| --- | --- | --- | --- | --- |
| `0.03 / 0.01` | `2.19e-3` | `3.31e-2` | `2.18e-3` | `3.33e-2` |
| `0.003 / 0.001` | `2.20e-5` | `1.19e-3` | `4.91e-3` | `9.57e-3` |
| `0.0003 / 0.0001` | `3.58e-7` | `9.76e-6` | `3.52e-2` | `5.58e-2` |
| `0.00003 / 0.00001` | `3.41e-7` | `2.49e-7` | `3.41e-1` | `5.77e-1` |

FP64 therefore converges through reconstruction at roughly the float response-storage floor in these cases.
Small-amplitude mixed refits remain unreliable even with the correct subtraction ordering, because that ordering
cannot recover information already lost inside the PCA kernel. At `epsilon=1e-5`, mixed versus FP64 templates at
the same amplitude differ by up to `29.0%` (24 frames) and `74.2%` (96 frames). Agreement at the largest amplitudes
still does not establish a negligible finite-amplitude bias. These results do not set a universal refit contrast.

A new experimental Catch2 case in [`P4Reduction_test.cpp`](../../tests/common/P4Reduction_test.cpp) checks the
reconstruction against an independent projector derivative, with and without derotation. It captures the exact
baseline and signed ingress matrices, uses Eigen's full temporal eigensystem to evaluate the derivative, then
applies the production geometry's interpolation weights and an independently accumulated FP64 mean. One on-axis
sample at its radial node avoids conflating this check with angular averaging. The source direction comes from
the signed ingress difference at `epsilon=1e-5`, so this specifically tests differentiation and reconstruction
after source sampling. Relative response errors were `1.8e-8`–`1.3e-7` across two modes and both derotation settings;
the regression bound is `2e-5`. The same case verifies unchanged final science pixels.

Artifacts: `/tmp/p4-response-convergence-20260917/refit{24,96}_{D64,M32D64}`. Each directory contains nine reductions:
one baseline plus eight response-enabled runs. All 32 response-enabled final science images equal their respective
baselines exactly. These short subsets and angular averages establish an oracle for this scope; they do not
replace the full scientific and performance validation in Step 3.

Verification: the experimental P4Reduction suite passed 34 cases / 354282 assertions (the new case accounts for
4340 assertions); the default build passed 27 cases / 77657 assertions. Both builds succeeded, modified C++ was
formatted, and the runner passed Python syntax checking and all four complete sweeps. The mxlib audit checked
the new test's float FITS constructor/read/write and mutable float-cube access calls: their exact instantiations
are present, and all executable lines in their function ranges are covered in the current LCOV report
(`fitsFile.hpp`: 491/491; `eigenCube.hpp`: 185/185). No new coverage gap was found.

### Analytic direct-P4 numerical kernel (2026-09-17)

Step 2 now has a production numerical API, `P4PCA::calculateResponse()`, in
[`P4PCA.hpp`](../../src/common/P4PCA.hpp) and [`P4PCA.cpp`](../../src/common/P4PCA.cpp). It evaluates the FP64 derivative
of the uncentered in-sample residual for supplied baseline arrays `X,y` and source directions `A,s`, including both
coefficient adaptation and motion of the retained subspace. It uses the complete smaller Gram eigensystem. On the
predictor-Gram path, unnormalized columns `X v_j` represent discarded predictor modes; an implicit complement term
includes the entire temporal nullspace without dividing by discarded zero/tiny singular values. Retaining all
predictor modes can therefore still give a nonzero response. Retaining all temporal modes gives zero.

Each requested count receives its derivative, baseline numerical rank, normalized cutoff gap, and one of
`differentiable`, `rankInsufficient`, `rankBoundary`, or `cutoffUnresolved`. Unavailable response columns are NaN.
The default numerical separation floor is `64 * epsilon_double * max(T,P)`, relative to the largest eigenvalue;
callers may increase it. Repeated eigenvalues entirely within a retained or discarded group are accepted. This
floor is a numerical-resolution guard, not a bound on response error or a physical source-amplitude tolerance.

The API is available in normal and experimental builds. It is not yet selectable through `psfResponse.method`:
the sparse sampling/configuration path still uses the existing paired refits. Automatic per-mode fallback and
product metadata belong to that next integration step. Numerical solver failures still throw; they are distinct
from unresolved per-mode boundaries. The derivative describes the FP64 mathematical fit, not differentiation of
the rounded mixed-precision program.

#### Analytic response verification

The known-eigensystem fixtures now compare the actual analytic API to their independently prescribed derivative
to relative tolerance `5e-12`, alongside the existing central-difference convergence tests. Coverage includes
both Gram orientations, target-only and predictor-only perturbations, deficient rank, internal repeated
eigenvalues, all retained predictor modes, all retained temporal modes, zero rank, narrow resolved gaps,
unresolved cutoffs, rank boundaries, output reuse, invalid source arrays, and malformed/failed eigensolvers.

The reconstruction integration case now also checks the analytic kernel against its independent Eigen derivative,
then reconstructs the analytic response using the production spatial weights and an independent FP64 masked mean.
Its comparison to the published FP64 refit template passes with and without derotation, while final science pixels
remain unchanged. This verifies the numerical kernel through the reconstruction oracle; it does not yet exercise
a user-selectable analytic-product path.

[`p4ResponseBenchmark.cpp`](../../benchmarks/p4ResponseBenchmark.cpp) replays captured baseline and unit-source
arrays through the actual analytic API and FP64 paired refits. The maintained
[`analyze_p4_analytic_response.py`](scripts/analyze_p4_analytic_response.py) compares both to the independent NumPy
SVD derivative. It records capture hashes, commands, native response arrays, diagnostics, executable and numerical
shared-library hashes, and summary tables. All three calculations use the same sampled unit-source direction;
the replayed paired fits perturb the captured FP64 arrays directly, so this comparison introduces no further
amplitude-dependent FP32 source interpolation.

| Frames | Captured fits / mode comparisons | Maximum analytic relative error per fit | Maximum analytic error after aggregation by position/mode | Maximum paired-refit error after aggregation, `epsilon=1e-5` |
| --- | --- | --- | --- | --- |
| 24 | 305 / 915 | `3.22e-12` | `6.31e-13` | `1.15e-8` |
| 96 | 351 / 1053 | `7.15e-11` | `2.11e-11` | `8.68e-8` |

All 1968 detector/mode comparisons were resolved and supported. Aggregation uses the norm over all captured time
samples at a position/mode, not an average of relative errors. Artifacts are under
`/tmp/p4-response-convergence-20260917/analytic24` and `analytic96`.

```sh
cmake --build _build_response_precision --target p4ResponseBenchmark -j2
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 python3 agents/plans/scripts/analyze_p4_analytic_response.py \
  --sweep /path/to/D64-sweep --binary _build_response_precision/benchmarks/p4ResponseBenchmark \
  --output /path/to/new-analytic-analysis
```

This correctness-first kernel is not yet a demonstrated speedup: summed single-call kernel times were about
`0.149 / 0.109` seconds (analytic / two refits) over the 24-frame captures and `1.64 / 1.38` seconds over the
96-frame captures. These timings exclude source sampling, image reconstruction, I/O, and process startup, and
do not measure reused baseline factorizations across source directions. Optimize and benchmark that reuse during
integration before claiming a computational gain.

Verification: both builds passed the response suite (10 cases / 1119 assertions). The experimental build also
passed P4PCA (42 / 4677) and P4Reduction (34 / 354318); the default build passed P4PCA (33 / 2220) and P4Reduction
(27 / 77657). Modified C++ sections were formatted, the replay script passed Python syntax checking, and both
real-data replay analyses completed with all requested comparisons available.

The mxlib coverage audit verified `isFinite<double>()` (7/7 executable lines), the FP64 `syevrMem` constructor and
destructor (3/3 each), and the reused `eigenSYEVR` array overload (60/60), along with the previously audited float
FITS/cube calls in the reconstruction test. Their exact instantiated function records are present in the current
LCOV report. No new upstream coverage gap was found.

The next integration milestone is recorded below. Baseline-factor reuse and the accepted AF Lep scientific result
belong to Step 3; covariance estimation and noise weighting remain separate later steps.

### Analytic sparse PSF products integrated (2026-09-17)

Step 2 now exposes `psfResponse.method=analytic` in production P4. Each sparse source measurement independently
samples baseline predictors/target and the unit-source direction, calculates the direct projector derivative in
FP64, and sends the detector time series through the same local reconstruction, mean combination, angular
averaging, and radial model as `refitDifference`. The ordinary science precision policy is unchanged. The source
sampler and published response arrays still use FP32; the analytic calculation does not recover precision lost
before its inputs are promoted. KLIP explicitly rejects this P4-only method.

The supported scope is detector-frame, uncentered, in-sample P4 without temporal augmentation or post-median
subtraction, with sparse radial sampling and model output or filtering. Science `mean` and `sigmaMean` are
supported; the response uses an unclipped mean and preserves configured weights. This is not a derivative of
sigma-clipping decisions. Known-companion trajectory avoidance uses the existing sparse-sampling rule.

Boundary policy is explicit:

- `psfResponse.refitContrast=0` leaves unresolved cutoff/rank-boundary modes unavailable.
- A positive half-amplitude enables paired FP64 refits of the same captured unit direction,
  `X +/- epsilon*A` and `y +/- epsilon*s`, for unresolved modes. Both signed fits must support the requested mode.
  This fallback is a finite-amplitude estimate; it does not certify differentiability at a degenerate boundary.
- Modes above the baseline numerical rank remain unavailable even if injection increases rank. Numerical failures
  propagate as errors. The standalone `refitDifference` method retains its existing precision-policy behavior.
- `psfResponse.analyticGapTolerance` is an optional larger relative boundary-resolution floor, default zero; the
  kernel's dimension-scaled FP64 numerical floor always applies.

Schema-8 products record response precision separately from science precision, fallback policy/amplitude, and
per-mode detector-fit counts. The new `measurement_diagnostics.fits` product records one row per measurement and
mode: `sourceIndex,row,column,modeIndex,analytic,rankInsufficient,rankBoundary,cutoffUnresolved,fallbackAttempted,
fallbackAccepted,unavailable`. A fallback retains its original boundary reason. Repeated reductions reset the
counts. A conservative scratch estimate limits analytic response workers under the configured memory budget.

The integration check also exposed an existing mismatch between detector-coordinate response products and
automatically cropped final images. Schema 8 publishes coordinates in the final-image frame and records the
detector-origin offset; its diagnostics use the same frame. Internal reconstruction still uses detector coordinates.
The cropped AF Lep products now load directly through the existing `hciAnalyze` external-filter interface.

#### Full-product equivalence

[`run_p4_analytic_products.py`](scripts/run_p4_analytic_products.py) reuses the completed FP64 refit experiment's
input/geometry provenance and checks current input hashes. It preserves the executed commands, binary/library
fingerprints, products, diagnostic counts, and per-stamp errors. Each run performs an independent science-only
baseline and an analytic-product reduction, then compares all 800 published response stamps at each of the three
mode fractions to the FP64 `epsilon=1e-5` refit products. The reference is the previously validated finite-amplitude
approximation, not an exact derivative.

| Frames | Analytic detector fits per mode | Compared published stamps across modes | Maximum relative error per stamp | Fallback / unavailable modes |
|---|---:|---:|---:|---:|
| 24 | 610 | 2400 | `3.15e-7` | `0 / 0` |
| 96 | 702 | 2400 | `1.96e-7` | `0 / 0` |

Both experiments passed with D64 and M32D64 science policies. Enabling response generation with filtering disabled
left the science images bit-for-bit unchanged, including NaN locations. Response support matched the refit
products, and analytic templates were bit-for-bit identical between science precision policies. Per-mode aggregate
relative errors were at most `1.63e-7`. All 3936 detector-mode outcomes across the two frame counts used the analytic
path. Headers agreed with the measurement diagnostics. These results validate the sampled/reconstructed response;
they do not establish photometric calibration or detection improvement.

Artifacts are under `/tmp/p4-response-convergence-20260917/analytic_products_v4_{24,96}_{D64,M32D64}`. Reproduction:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 python3 agents/plans/scripts/run_p4_analytic_products.py \
  --reference /path/to/completed/refit24_D64 \
  --binary _build_response_precision/benchmarks/p4ReductionPrecisionBenchmark \
  --precision P4-M32D64 --output /path/to/new-analytic-products
```

The integration test covers pure analytic responses, forced and mixed boundary fallback, disabled fallback,
baseline rank insufficiency, source avoidance, filtering, automatic crop coordinates, repeated reductions, and
per-mode provenance. The experimental P4Reduction suite passed 35 cases / 354742 assertions; the default build
passed 28 / 78081. KLIP passed 47 / 9441812, including rejection of `analytic`, and the default `hciAnalyze` suite
passed 13 / 166. The actual `hciAnalyze` executable also consumed a schema-8 cropped AF Lep response field.
The default production `p4Reduce` executable reproduced the 24-frame M32D64 benchmark's baseline science, analytic
science, and response planes bit-for-bit (`/tmp/p4-response-convergence-20260917/analytic_products_default24`).
Modified C++ ranges were formatted and the new driver passed Python syntax checking.

The mxlib coverage gate rechecked the current filtered LCOV report for the edited orchestration/configuration,
header, FITS, cube, geometry, timing, finite-check, and FP64 workspace calls. Called executable-line ranges are
covered, but exact float instantiations remain unverified for the const cube image view, cube assignments, and
KLIP's float annular-geometry/view-extraction calls. Concrete upstream tests are listed under
Known non-blocking ownership follow-ups in [mxlib_cleanup.md](mxlib_cleanup.md).

### Step 3: scientific and computational validation (2026-09-17)

**Step 3 is complete.** The full-data comparison, all 84 injection reductions, 72 conditional fits, and controlled
timing comparisons are recorded in the [final report](results/p4-step3-20260917/README.md). The report preserves
photometric and spatial limitations and a CPU-placement reproducibility finding. The earlier
[review checkpoint](results/p4-step3-20260917-review/README.md) remains as a record of the 54-trial pause.

#### Shared factors, source sampling, and image coordinates

`P4PCA::prepareResponse()` owns an FP64 baseline eigensystem and its rank/cutoff diagnostics. Successive
`calculateResponse()` calls apply independently sampled source directions to that baseline without another
eigensolve. The one-shot API delegates to these operations; the response formula and boundary policy are unchanged.

The reduction gathers detector requests from bounded source batches and factors each unique valid detector fit
once per batch. `psfResponse.analyticBatchSize` defaults to 32; 1 provides a comparison without cross-source reuse.
The existing memory budget limits both source batches and workers using conservative geometry, source-cache,
residual, factor, and scratch estimates. Worker-local factors are discarded after their consuming directions.
FITS provenance records the total `P4 PSF ANALYTIC FACTOR COUNT` and realized `P4 PSF ANALYTIC BATCH SIZE`.
The progress log reports a separate factor count for each batch.

Each source batch caches shifted detector pixels once per source/frame. Predictor interpolation reuses those
values with the original column-then-row accumulation order. For temporal-Gram directions with at most 20%
nonzero entries, the kernel forms `A*X^T` with an exact sparse representation and adds its transpose. Only exact
zeros are omitted. Dense directions retain the original product. Both 24- and 96-frame checks preserve all 4800
published stamps bit-for-bit; maximum per-stamp differences from the FP64 paired-refit oracle remain `3.15e-7`
and `1.96e-7`. A separate 621-by-3000 numerical check gives a dense/sparse analytic difference of `2.69e-16` and
an analytic/FP64 paired-refit difference of `2.05e-9` at half-amplitude `1e-5`.

Automatic cropping uses the same residual crop, derotation, and final crop with integrated filtering enabled or
disabled. All newly written P4 response methods publish final-image coordinates and detector-origin cards.
Detector dimensions are recorded separately from input-template dimensions, and the reader handles odd crop
parity. Legacy manifests retain their template-size fallback; old unequal detector/template products without
detector-dimension cards should be regenerated. A 24-frame D64 AF Lep run with integrated filtering
preserves the science pixels bit-for-bit and produces 900 finite filtered pixels. The synthetic regression covers
analytic, paired-refit, detector-local, and exact-sky coordinates, including a PSF smaller than the detector.

The experimental suites pass 42 cases / 4677 assertions (PCA), 12 / 1440 (response), and 35 / 354802 (reduction).
The default build passes 33 / 2220, 12 / 1440, and 28 / 78141, respectively; `hciAnalyze` passes 13 / 166.
The mxlib audit verifies the called FP64 eigensolver range at 60/60 executable lines, workspace construction and
destruction at 3/3 each, and cleanup at 13/13. Exact const-float cube views and float assignment instantiations
remain the non-blocking upstream coverage follow-ups recorded in [mxlib_cleanup.md](mxlib_cleanup.md).

#### Full AF Lep result

The 621-frame run reproduces the accepted geometry: 232 measurements at 58 radii, 221 excluded source centers,
and mode fraction 0.15. All 176566 detector responses are analytic, with zero rank, cutoff, fallback, or unavailable
outcomes. It uses 30166 baseline factorizations, 82.9% fewer than preparing one baseline per direction.
The response field has the same coordinates and finite support as the accepted paired-refit field.

A separate response-free reduction using the identical executable and recorded application/math libraries
produces **byte-for-byte identical science pixels**. Replaying the accepted response field on the current science
also reproduces its fit: the contrast changes by only `-4.08e-11` and separation by `-3.88e-8` pixel.

| Response field | Contrast | Separation (pixels) | PA (degrees) | SNR |
| --- | ---: | ---: | ---: | ---: |
| Accepted paired refit | 0.004884882108 | 11.49649046 | 262.32662164 | 4.32275222 |
| Analytic | 0.004856152652 | 11.48754721 | 262.32672117 | 4.33444033 |

The analytic fit differs from the accepted matched-refit result by `-0.588%` in contrast and about `0.009` pixel
in position. Relative to the accepted negative-planet optimizer, the contrast offsets are `+1.936%` (analytic)
and `+2.539%` (paired refit), and position offsets are `0.2525` and `0.2436` pixel. These are comparisons with the
accepted point estimate, rather than injected-source bias measurements.

Across all 11304 published stamps, analytic versus finite-amplitude refit relative differences have median
`5.61%`, 95th percentile `6.26%`, and maximum `80.56%`; the aggregate cosine is `0.999049`.
The largest relative difference is at radius `58.71` pixels, with only 19 of 121 stamp pixels available and a
small reference norm. The median differences are `0.89%` over radii 6--12, `1.94%` over 12--24, `5.09%` over
24--42, and `5.95%` over 42--60 pixels. The report retains these outer-edge outliers.

At the accepted integer peak, comparison with the archived negative-planet removal gives cosine `0.97805` and
best-scaled shape error `20.84%` for the analytic field, versus `0.97844` and `20.65%` for the accepted refit field.
This comparison combines finite amplitude, subpixel registration, clipping, and source-support differences:
the full-image removal shifts the full 256-pixel PSF; the accepted optimizer uses a 13-pixel local window and
14-pixel source crop; the 11-pixel response uses a 12-pixel source crop. The latter retains 99.08% of the full
PSF's squared energy before reduction. That fraction does not bound the post-reduction mismatch. The broader
injection study below uses matching source support to separate this issue from response calibration.

#### Injection protocol and computational measurements

The completed study uses six integer detector positions `(120,137)`, `(137,120)`, `(109,143)`, `(146,112)`,
`(98,98)`, and `(157,157)`, at radii approximately 12.1, 24.1, and 41.7 pixels. At each position it evaluates
zero contrast and both signs of 0.25, 1, and 4 times the accepted contrast `0.004763925929356391`, under arithmetic
mean and 5-sigma mean combination: 84 local reductions of all 621 frames. The PSF is zero-padded after its central
12-pixel crop, preserving the original contrast scale while allowing a 15-pixel local fitting window.
Both combination experiments use the same mean-combined response field. The sigma-mean comparison therefore also
measures this approximation for clipped science; the response model does not differentiate the final clipping rule.

Analysis measures fixed-position gain, fitted contrast and astrometric offsets, template cosine and shape error,
and one-sided versus symmetric response differences. The free-position estimate uses a bounded quadratic peak
within one pixel, with a complete 11-by-11 template at each candidate. Baseline subtraction measures conditional
response calibration; raw-image fits and unsuccessful fit statuses are retained separately. Detection completeness
and false-alarm calibration remain the held-out tests in Step 5.

All 84 reductions are complete. Across the six positions, fixed-position analytic contrast bias has median
`+3.376%`, `+2.232%`, and `-21.087%` for mean combination at 0.25, 1, and 4 times the reference brightness;
sigma-mean gives `+3.373%`, `+2.230%`, and `-21.087%`. Paired-refit medians are `+5.543%`, `+3.727%`, and
`-19.209%` for mean and `+5.489%`, `+3.731%`, and `-19.209%` for sigma-mean. The largest individual analytic
gain change between combinations is 0.125 percentage point. Brightest-case analytic bias spans `-38.54%` to
`-1.99%` for mean and `-38.54%` to `-1.93%` for sigma-mean.

Median analytic position errors are `0.0397`, `0.0391`, and `0.0365` pixel for mean, and `0.0397`, `0.0392`, and
`0.0365` for sigma-mean; individual errors reach about `0.204` pixel. Median template/response cosine decreases
from `0.9878` and `0.9870` to `0.9611` for mean, versus `0.9876`, `0.9869`, and `0.9612` for sigma-mean.
Median one-sided versus symmetric response differences grow from `0.95%` to `3.61%` across the brightness range
for mean and from `1.12%` to `3.61%` for sigma-mean; individual cases reach `8.31%`.

All 72 conditional fits converge. Four raw-image fits reach the search boundary: the faintest first-position
trial for both fields and combinations. Their statuses are retained. These results show a brightness limit for
both fixed response fields and include finite-amplitude effects, sparse interpolation, and native arithmetic.
Bias is measured against the **known added contrast**, not the negative-planet optimizer. Baseline subtraction
cancels much of the shared residual noise, so the scatter across positions is not a detection-SNR estimate.
The [final table and figures](results/p4-step3-20260917/README.md) retain every position and fit status.

| Full-data run | Wall time (seconds) | Peak resident memory (GiB) |
| --- | ---: | ---: |
| Same-build science-only control | 310.21 | 4.93 |
| Science plus analytic response field | 9894.06 | 11.92 |

These are observed whole-process measurements with 20 OpenMP workers and one BLAS thread on an i9-12900HK,
using the Release build. Small development checks overlapped portions of the analytic run. The archived refit
run took 6574 seconds with 48 OpenMP workers; that historical comparison does not establish a controlled speedup.
The 24/96-frame trials compare batches, source caching, sparse arithmetic, and production-precision/FP64
paired-refit controls serially, using two workers, one BLAS thread, and medians of three repeats. All processes
are pinned to performance-core CPUs 0 and 2. All 36 analytic trials retain identical templates across batches
and arithmetic variants; all 12 paired-refit controls preserve same-precision science. Two extra science-only
baselines provide the FP64 controls.

| Subset configuration | 24 frames: seconds / MiB | 96 frames: seconds / MiB |
| --- | ---: | ---: |
| Uncached dense, batch 8 | 1.75 / 298.2 | 9.38 / 354.1 |
| Cached dense, batch 8 | 0.97 / 340.8 | 5.36 / 516.3 |
| Cached sparse, batch 1 | 1.00 / 298.2 | 5.06 / 354.1 |
| Cached sparse, batch 8 | 0.95 / 340.7 | 4.93 / 515.0 |
| Paired refit, M32D64 | 3.33 / 298.1 | 17.34 / 354.1 |
| Paired refit, D64 | 3.41 / 297.8 | 16.06 / 354.1 |

Batch 8 is the realized size for a request of 32 in these subsets. Native refit uses half-contrast
`0.004763925929356391`; FP64 refit uses `1e-5`. Analytic trials use M32D64 science with FP64 responses, while
the FP64 refit also uses FP64 science. Native refit maximum stamp differences from the analytic field are
`0.114%` and `0.760%`; FP64 differences are at most `3.15e-7` and `1.96e-7` in relative norm.
The subsets use 5-pixel stamps, eight measurements, and three mode fractions; their roughly 3.5-fold speedup
over native paired refits does not establish a full-dataset speedup. Source caching supplies most of the measured
improvement; larger batches also increase memory. The full report retains all batch-one and batch-eight cases.

**CPU-placement limitation:** the initial unpinned timing attempt stopped after 12 successful trials when its
next science equality check failed. Native science differs by up to `0.0024900436401367188` between performance
cores 0/2 and efficiency cores 12/13. Two repeats of both science-only and response-enabled reductions on each
CPU set reproduce the corresponding result exactly; enabling responses changes neither within a set. Response
templates also remain identical. The internal arithmetic mechanism is not yet isolated. The corrected benchmark
pins one core class and retains strict equality checks; the original failure and diagnostics are archived.
Cross-core-class bitwise reproducibility is not asserted. The full-data and injection results retain their
original unpinned 20-worker execution policy and describe those recorded runs.

The maintained drivers are [`run_p4_step3_products.py`](scripts/run_p4_step3_products.py),
[`analyze_p4_step3.py`](scripts/analyze_p4_step3.py),
[`run_p4_step3_injections.py`](scripts/run_p4_step3_injections.py),
[`benchmark_p4_step3_reuse.py`](scripts/benchmark_p4_step3_reuse.py), and
[`summarize_p4_step3.py`](scripts/summarize_p4_step3.py). Final records are archived at
`working/roc/p4_analytic_step3_20260917`, with compact summaries and inspected figures under
[`results/p4-step3-20260917`](results/p4-step3-20260917/README.md). The archive preserves frozen software,
commands, input hashes, raw products, resources, analysis scripts, and the failed unpinned timing attempt.
All 621 input, configuration, and PSF hashes were reverified after the injection study. The original review
archive remains at `working/roc/p4_analytic_step3_20260917_review`.

**Step 3 limitations:** Bright-source response calibration, outer-edge support, and CPU-class reproducibility
remain explicit limitations of the measured scope; they are not claims of universal calibration or detection
completeness.

### Step 4: independent residual-noise weighting (2026-09-17)

**Step 4 is complete.** The shared filtering API now accepts identity, diagonal, and regularized PCA residual-noise
models. Response
estimation remains independent of noise weighting. Implementation is in
[`PSFNoiseModel.hpp`](../../src/common/PSFNoiseModel.hpp),
[`PSFNoiseModel.cpp`](../../src/common/PSFNoiseModel.cpp), and
[`P4PSFFilter`](../../src/common/P4PSFFilter.hpp).

#### Model and filter contract

- `PSFNoiseModel::identity(p)` uses unit variance and zero mean. `P4PSFFilter::calculateWeighted(..., noise)`
  delegates this choice to the existing `calculate(...)` path, preserving its accumulation order, amplitude,
  correlation, normalization, support fraction, and validity.
- `diagonal(variances, mean)` accepts positive per-pixel variances and an optional background mean.
  `lowRank(variances, factor, mean)` represents $C=D+LL^T$; the factor need not have orthogonal columns.
- `estimateDiagonal(samples, varianceFloor)` and `estimatePCA(samples, maximumModes, varianceFloor)` take a finite
  $n\times p$ matrix of already selected noise-only stamps, with $n\ge2$. They retain the sample mean and use
  $n-1$ covariance normalization. The diagonal estimator bounds each variance below by the supplied floor. The PCA
  estimator retains at most `min(maximumModes,n-1,p)` modes with eigenvalues above the floor, with factors
  $L_i=u_i\sqrt{\nu_i-\mathrm{varianceFloor}}$. Discarded directions retain the isotropic floor. Zero requested modes
  gives a floor-only covariance with the estimated mean.
- Pixel ordering is `row + column * stampRows`. Samples, response, science, and mean must share coordinates and
  image units. The floor is an absolute variance in squared image units and must be finite and strictly positive.
  The caller selects source-free training footprints; the estimator does not infer exclusions or correct for
  dependence among overlapping samples.
- Filtering subtracts the model mean from science only. It selects response-valid, in-bounds, finite science pixels,
  restricts both $D$ and the rows of $L$ to that support, and solves the resulting marginal covariance. It never
  selects a submatrix of the full inverse. The original odd-stamp geometry and full-area support-fraction policy
  still apply, including rectangular stamps.
- Results contain signed amplitude, weighted correlation $N_C$, normalization $E_C$, support, validity,
  `conditionalSigma` $=E_C^{-1/2}$, signed `score` $=N_C/\sqrt{E_C}$, and `nonnegativeLogLikelihood`
  $=\max(0,\mathrm{score})^2/2$. Invalid support or unrepresentable statistics give an invalid result; malformed
  inputs or failed covariance solves throw. Conditional uncertainties and scores assume the supplied mean,
  covariance, and response are fixed. They are not calibrated detection significances.

For example, after selecting and aligning the noise-only training stamps:

```cpp
const auto noise = mx::improc::PSFNoiseModel::estimatePCA(trainingSamples, maximumModes, varianceFloor);
const auto result = mx::improc::P4PSFFilter::calculateWeighted(
    science, response, validity, sourceRow, sourceColumn, minimumSupportFraction, noise);
```

#### Stable low-rank solve

The implementation uses an orthogonal-basis form of the same covariance solve as Woodbury. With
$B=D^{-1/2}L=Q[R;0]$, transform the right-hand side by $Q^TD^{-1/2}$, solve $I+RR^T$ in the leading coordinates,
leave the complementary coordinates unchanged, then transform back by $D^{-1/2}Q$. Householder QR and a small
Cholesky factorization avoid explicitly forming a $p\times p$ covariance or inverse.

A regression case exposed cancellation in the direct subtractive Woodbury formula: for $D=I$ and a factor
$L=(10^8,0,0)^T$, the precision along the first coordinate rounded to zero instead of approximately $10^{-16}$.
The orthogonal implementation retains this finite weight, including when only that one pixel remains in support.
Each call rebuilds the factorization for the retained pixels. For low rank $r\ll p$, setup/solve costs
$O(pr^2+r^3)$ and storage is $O(pr+r^2)$; repeated-support caching is a possible later optimization.

#### Verification and next work

The Release build passed all five selected regression suites (101 test cases):

| Suite | Test cases | Assertions |
| --- | ---: | ---: |
| `PSFNoiseModel` | 5 | 164 |
| `P4PSFFilter` | 8 | 291 |
| `P4Reduction` | 28 | 78,141 |
| `KLIPreduction` | 47 | 9,441,812 |
| `hciAnalyze` | 13 | 166 |

The noise/filter tests cover explicit identity equivalence to the original long-double accumulation, diagonal
and correlated dense-Cholesky references, reordered/masked/edge/NaN support, signed and rescaled templates,
covariance rescaling, background subtraction, sample centering and normalization, PCA truncation and its finite
complement, rank-deficient factors, strongly downweighted directions, and malformed or unrepresentable inputs.
`p4Reduce`, `klipReduce`, and `hciAnalyze` were rebuilt against the extended result structure. The application
regressions exercise their existing identity-filter behavior.

The final run used the existing `_build_fresh` Release configuration and homogeneous CPU affinity:

```sh
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 taskset -c 0,2 ctest --test-dir _build_fresh --output-on-failure -V \
    -R 'hcireduceTest_(common_(PSFNoiseModel|P4PSFFilter|P4Reduction|KLIPreduction)|apps_hciAnalyze)_test_cpp'
```

Changed C++ sections were formatted; new files pass `clang-format --dry-run --Werror`, and `git diff --check`
is clean. A focused Doxygen build produced no warnings and confirmed test references on both filter entry
points and the noise-model public APIs. No edited function calls an upstream mxlib API, so this step adds no
mxlib coverage-ownership follow-up.

Step 4 supplied the common numerical API. Application integration and annular training are tracked in Step 5
below. Rank and floor selection, uncertainty calibration, and completeness at a fixed false-positive rate require
held-out data. No covariance-related gain on AF Lep is claimed by the Step-4 numerical checks.

### Step 5 started: annular training and application integration (2026-09-17)

**Initial integration checkpoint; the completed fixed-policy evaluation is recorded below.** `PSFNoiseTraining`
now selects complete, finite annular patches and applies the Step-4
shared filter. `hciAnalyze` exposes this for both P4 and KLIP responses through `noise.model=diagonal|pca`.
Identity remains the default; `noise.outputDiagnostics=true` also exports its conditional products for comparisons.
The [configuration dictionary](../../doc/hciAnalyze_config.dox) documents the options and output roles.

#### Geometry and diagnostic contract

Training centers stay at the candidate's exact radius. The default arc step is the response half-width, with a
minimum of one pixel; the number of centers is rounded upward and their absolute angular phase is fixed. Each
training offset is rotated from the candidate angle to the training angle before bilinear sampling. Thus training
stamps are expressed in the candidate's native pixel coordinates, while science and response retain their original
sampling. All patches have the same radial footprint, including across concentric reduction-region boundaries.
Interpolation changes noise statistics; calibration must measure that consequence.

A training patch is removed if its enclosing circle, enlarged by `sqrt(2)` pixels for interpolation support,
intersects the candidate's enclosing circle plus `noise.guardRadius`, any resolved source-exclusion circle, or any
explicit `noise.excludeRows/Columns/Radii` circle. This conservative rule excludes whole overlapping neighborhoods.
Remaining patches require complete finite interpolation support. Defaults are eight accepted samples, three noise
PCA modes, and a floor of 0.1 times the median pixel variance of the centered training matrix. Those are explicit
initial settings, not settings chosen to improve the known planet. The full response footprint determines the patch
size: the existing 11-by-11 response uses a 5-pixel default arc step, rather than silently assuming a one-lambda/D
patch. The approximate 19-patch example in Section 6 still applies when the patch half-width is one lambda/D.

No-training and invalid-support outcomes have explicit status maps and never fall back to a different noise model.
Outputs separately retain amplitude, conditional sigma, signed score, support, accepted/attempted/excluded/incomplete
sample counts, actual rank, and absolute floor. Metadata preserves the covariance settings, response field, mode,
coordinate convention, exclusions, and uncalibrated status. Gaussian post-filtering is disallowed for this mode so
those conditional products remain tied to the stated response and covariance. The existing annular SNR calculation
is a separate output; its small-sample multiplier does not calibrate the conditional score.

#### Fixed-policy integration pilot

The maintained [`run_p4_step5_pilot.py`](scripts/run_p4_step5_pilot.py) freezes the analysis executable/library,
fingerprints the full Step-3 science image and response products, and runs identity, diagonal, and PCA sequentially
on homogeneous CPU cores. It measures runtime and valid training coverage by radius, without choosing settings
from the planet's score. This is an integration pilot, not a completeness or false-positive measurement.

The [pilot report](results/p4-step5-pilot-20260917/README.md) records the completed fixed-policy run. Identity,
diagonal, and PCA took 0.121, 1.622, and 7.957 seconds including diagnostic I/O. Of 8,616 identity-valid positions,
7,679 admit covariance training; 937 fail the eight-patch minimum. Another 2,188 response positions fail the original
support policy for all methods. Covariance training is complete from 20–50 pixels, with median accepted patch counts
of 20, 39, and 51 in the 20–30, 30–40, and 40–50 pixel bands. These overlapping counts are not independent sample counts.
At the nearest AF Lep pixel (139, 126), only six patches survive, so the fixed eight-patch policy marks both covariance
models invalid. Inner-radius footprint and regularization choices remain a development/held-out calibration task;
neighboring valid aperture pixels do not establish a gain at the source.

All six selected regression suites pass (107 test cases). The new geometry/filter suite has 1,141 assertions, and
`hciAnalyze` has 445. A focused Doxygen build verifies production API test links; its included-fixture preprocessor
warning is also present with the pre-change test file. The current mxlib LCOV audit found an unexercised exact
unsigned-long configuration overload, a 34/40-line CLI parser, and missing exact float cube lifecycle instantiations.
Concrete upstream tests are listed under `Known non-blocking ownership follow-ups` in
[`mxlib_cleanup.md`](mxlib_cleanup.md); passing downstream tests do not establish upstream coverage.

#### Held-out evaluation protocol and remaining work

1. Define disjoint angular neighborhoods before inspecting their scores. Reserve separate development,
   threshold-calibration, and evaluation neighborhoods, with full-footprint guards. Explicit exclusion circles can
   cover those neighborhoods; record their union and verify that no accepted training stencil reaches a held-out
   pixel. Use the same candidate support and exclusions for model comparisons.
2. Fix an initial operating point of 5% false positives per preassigned spatial trial. A trial uses a predeclared
   search aperture around its center, identically for nulls and injections. Thresholds come only from calibration
   neighborhoods, with achieved rates and uncertainty reported on evaluation neighborhoods. Count spatial trials
   explicitly and estimate uncertainty by angular blocks; do not count overlapped pixels as independent nulls or
   extrapolate this pilot to five-sigma tails.
3. Compare identity, diagonal, and PCA weighting, with a zero-mode PCA control to separate background centering and
   local variance normalization from correlated weighting. Select rank/floor on development neighborhoods only;
   freeze them before evaluation. Also report results for the fixed initial policy.
4. Run fresh positive injections through the reduction at held-out positions and brightnesses spanning the detection
   transition. Use raw recovered amplitudes and scores, with no baseline subtraction, for completeness and
   uncertainty coverage. Keep template/PSF normalization, response support, mode count, and aperture policy fixed.
   The six-position Step-3 baseline-subtracted study remains a response-bias diagnostic: its local 15-by-15 outputs
   neither supply annular training data nor establish the global effect of injections on the residual field.
5. Report recovery fraction at the frozen threshold, null exceedance rate, contrast bias, and empirical coverage of
   conditional intervals. Keep known-planet inspection separate. Record mask/sample-count failures and compare on
   common eligible support so a method cannot appear better merely by rejecting difficult locations.

The present checkpoint implements the sampler, application integration, and pilot machinery. The held-out
threshold and fresh-injection study are required before Step 5 can be marked complete or covariance gains claimed.

### Step 5 continued: exclusion audit and held-out null calibration (2026-09-18)

The discussion of the pilot's **10-pixel exclusion radius** exposed an unjustified inherited setting: it came
from `hciAnalyze`'s existing signal/SNR mask, rather than a response-support measurement or the Welch geometry.
The dataset's `working/analyze.conf` specifies **3.6 pixels per lambda/D**, making ten pixels about 2.8 lambda/D.
The first pilot explicitly used the application's 2.5-pixel default instead. Covariance sampling is specified in
pixels; `lambdaD` affects its physical interpretation and the separate legacy SNR calculation. New pilot commands
require an explicit scale and source radius. The old pilot remains archived with its actual settings.

The old footprint rule also added an 8.49-pixel enclosing/interpolation radius to that source mask, excluding
training centers within 18.49 pixels. Its separate candidate rule excluded centers within 15.56 pixels; these are
union tests, not additive guards. Consequently, six accepted patches at AF Lep was a property of that conservative
policy, not evidence of an intrinsic covariance-sampling limit.

`noise.exactExclusion=true` now checks every native pixel read with a nonzero bilinear weight. It withholds the
candidate's rectangular stamp, optionally expanded by a Euclidean guard, and native pixel centers inside the
specified exclusion circles. Any touching training patch is rejected in full. Enclosing circles are only an early
overlap screen; no interpolation margin is added a second time. The original enclosing-circle policy remains the
default. `noise.only=true` writes conditional maps without requiring an empirical SNR measurement at a listed source;
this supports sparse held-out response fields and explicit invalid-source outcomes.

The experiment retains the **11-by-11 response and five-pixel arc spacing**. At 10–20 pixels separation, a 3.6-pixel
circle contains median 80.1% of the stored response energy and 29.3% of its negative-lobe energy. Those fractions are
relative to the finite stored stamp, not the entire processed source. A provisional **7.3-pixel source mask** encloses
the stamp's 7.07-pixel half diagonal plus the known source's 0.21-pixel offset from the nearest response center.
Full-image injections measure the remaining wings and changes outside this declared footprint; the mask is not
asserted to make training source-free.

The [development and null-calibration report](results/p4-step5-development-20260918/README.md) records the geometry
ablation. Reducing the circle radius alone gives eight AF Lep patches; the exact stencil rule alone with radius ten
gives seven; using both gives eight and expands covariance-valid coverage from 7,679 to 7,919 positions. There are
15 proposed patches at the nearest AF Lep pixel. The approximately 19-patch example in Section 6 assumes a patch
half-width of one lambda/D.

The science FITS header audit also corrected the old report: the integration pilot used **sigma-mean science with
mean response templates**. The fresh full-image development/evaluation reductions explicitly use mean combination,
all 621 frames, fixed homogeneous CPU cores, the same unrenormalized 12-pixel source crop, and the existing 11-pixel
mean response field. The maintained full-image runner replaces neither entire images nor training annuli with local
injection cutouts.

For the held-out pilot, development, calibration, and evaluation pixel neighborhoods are disjoint. Their union is
excluded from every covariance fit. Sixty-four preassigned null searches use the same radius-one, five-pixel search
aperture; 28 calibration and 28 evaluation searches have common support, while the eight searches at radius 20 fail
the eight-patch minimum. Independent stencil checks match production counts at every search pixel. Inner-radius
training coverage in the earlier single-candidate pilot does not imply that these stricter simultaneous holdouts
can be supported there.

Rank three and floor fraction 0.1 remain fixed, with identity, diagonal, and zero-mode PCA controls; no model setting
was selected from evaluation scores. For a target 5% per-search false-positive rate, the prescribed conservative
order-statistic rule selects the largest of the 28 calibration scores. The separate evaluation exceedances are
2/28 (identity), 0/28 (diagonal), 1/28 (zero-mode PCA), and 1/28 (three-mode PCA). These overlapping spatial trials
and four angular blocks do not establish a precise tail probability or a detection gain.

Eighteen fresh full-image evaluations were launched at six preassigned evaluation sites, at 0.5, 1, and 2 times the
contrast scale from the frozen zero-mode threshold and training-only conditional sigma. Their runner freezes the
software, inputs, thresholds, masks, and response subset; it then measures all four filters automatically. Recovery
uses raw injected images, with no baseline subtraction. Invalid searches remain nondetections, and amplitude bias
and conditional-interval coverage use the exact injected position. The completed review follows.

### Step 5: initial held-out evaluation complete (2026-09-18)

All **18 full-image injections and 72 model measurements** finished successfully, with no invalid injection searches.
The [evaluation report](results/p4-step5-evaluation-20260918/README.md) preserves individual results, recovery and
photometry plots, and the independent product audit. All 621 input frames and recorded frozen/product hashes agree;
all 72 measurements reproduce exactly from their FITS maps. The reduction wall times total 1 h 58 min 52 s,
excluding the subsequent filter analysis.

Each brightness is a multiple of a site-specific reference contrast: the frozen zero-mode PCA threshold times
that site's baseline conditional sigma. It is not a measured SNR or each model's own threshold contrast. Three
brightnesses reuse six sites at nominal separations 26, 38, and 50 pixels, so these are not 18 independent noise
realizations. The rank-three/floor-0.1 policy and all thresholds remained fixed throughout evaluation.

| Method | Detections at 0.5× | Detections at 1× | Detections at 2× | Evaluation null exceedances | Positive amplitudes within ±1 conditional sigma |
| --- | ---: | ---: | ---: | ---: | ---: |
| Identity | 3/6 | 6/6 | 6/6 | 2/28 | 18/18 |
| Diagonal | 1/6 | 4/6 | 6/6 | 0/28 | 3/18 |
| Zero-mode PCA | 2/6 | 4/6 | 6/6 | 1/28 | 0/18 |
| Three-mode PCA | 2/6 | 5/6 | 6/6 | 1/28 | 1/18 |

**No covariance detection gain is established.** The target false-positive rate is common, but the observed null
rates differ and their uncertainty is large. Identity's greater recovery count is not a matched-rate superiority
result. Three-mode PCA recovers one more middle-brightness site than the zero-mode control; six sites cannot
establish a general gain. Some evaluation sites already exceed threshold before injection and remain in the
preassigned sample; the report identifies them explicitly.

The learned conditional sigmas fail to describe the raw amplitude errors in this sample. Identity's sigma assumes
unit pixel covariance, so its 18/18 coverage does not demonstrate successful noise calibration. At unmodified
evaluation null centers, the corresponding ±1-sigma coverage is 28/28, 6/28, 2/28, and 1/28, respectively; the
mismatch is therefore not solely finite-source nonlinearity. Scores remain uncalibrated as Gaussian significances.

Raw median contrast errors decrease with brightness: identity +33.6/+18.5/+10.8%, diagonal +20.4/+12.2/+8.0%,
zero-mode PCA +40.7/+22.0/+12.4%, and three-mode PCA +32.7/+17.1/+9.2%. These include substantial preexisting
background offsets. A **secondary paired-amplitude increment** subtracts the separately measured baseline amplitude:
its median errors are about +2–2.4%, with individual errors +0.36% to +6.14% across all methods and levels.
Those increments are neither negative-injection results nor raw completeness/coverage measurements. Learned weights
adapt between the baseline and positive image, so they do not isolate the response derivative alone.

This completes the initial fixed-policy Step-5 evaluation and supplies a review checkpoint. Remaining scientific
work is to investigate conditional-uncertainty calibration and rank/floor choices on development data, resolve
small-separation training support and source wings, and evaluate any revised policy on fresh held-out data with
more independent null/injection coverage. The inspected evaluation sample cannot serve as an untouched test of
subsequent tuning. No production code changed during this results review.

### Step 5 reference added: Gaussian smoothing (2026-09-18)

The user pointed out that the identity matched filter had not been shown to improve on the existing Gaussian
low-pass filter, and requested **`filter.lpfGaussFW=3.6`** as a reference. The
[Gaussian comparison](results/p4-step5-gaussian-20260918/README.md) now applies that fixed width to the same saved
baseline and 18 positive images, with the same 28 calibration searches, 28 evaluation searches, common eligibility,
and five-pixel search. New thresholds follow the original calibration-only rule and were frozen before processing
positive images. The extension was requested after inspecting the initial evaluation; it is not a fresh blind test.

| Statistic | Detections at 0.5× / 1× / 2× | Evaluation null exceedances |
| --- | --- | ---: |
| Gaussian 3.6 px, application annular SNR | 3/6, 6/6, 6/6 | 2/28 |
| Identity matched filter, original conditional score | 3/6, 6/6, 6/6 | 2/28 |
| Gaussian 3.6 px, smoothed intensity only | 3/6, 5/6, 6/6 | 2/28 |
| Identity matched filter, application annular SNR | 4/6, 6/6, 6/6 | 3/28 |

**No advantage over the Gaussian reference is established.** The usual Gaussian/SNR pipeline matches the original
identity recovery counts and observed null rate. Identity's one extra faint recovery when both filters use the
application SNR path comes with an additional null exceedance. The smoothing-only Gaussian misses one more
middle-brightness injection, but this six-site sample cannot establish a general difference.

The references separate the existing application's annular normalization from the smoothing operation. The new
optional `hciGaussianReference` benchmark exports the actual production Gaussian filter before SNR normalization.
For all 19 images, passing that output through the ordinary SNR path with smoothing disabled reproduces the direct
Gaussian/SNR output bit for bit. An independent mask-aware FP64 convolution agrees within 1.19e-6 of the peak
absolute image value. The production kernel is 15×15 at FWHM 3.6; its calibration/evaluation input-support unions
are disjoint. Both SNR controls retain the ordinary application's use of trial neighborhoods in annular
normalization, unlike the held-out covariance estimator; the report records this distinction and full-field
support differences. Raw Gaussian intensity is not reported as recovered contrast.

No production filter or application default changed, and no reduction was rerun. The Gaussian reference remains
part of subsequent validation, alongside the identity and covariance controls. Any revised policy needs fresh
held-out confirmation after its development choices are fixed.

### Step 5 development: radial pooling and subsequent ROC validation (2026-09-18)

The user requested testing patches at different radii, motivated by the larger OR required for good P4 reduction,
and suggested radial-variance normalization as a separate experiment. The four-way development comparison and
consistent data/template normalization are specified in Section 6 above. Production covariance training still
uses one radius; this checkpoint extends only the independent audit sampler and records a baseline feasibility check.

With the existing union of known-source and development/calibration/evaluation exclusions, accepted patches are:

| Candidate radius (pixels) | Same radius | ±5-pixel band | ±10-pixel band | ±20-pixel band |
| --- | ---: | ---: | ---: | ---: |
| 12.1, development | 3 | 8 | 15 | 37 |
| 24.1, development | 7 | 22 | 39 | 90 |
| 41.7, development | 24 | 69 | 82 | 99 |
| 20.4, previously insufficient null site | 6 | 19 | 34 | 74 |

These are overlapping counts, not independent samples or detection measurements. The native radial standard
deviation is roughly 0.13–0.17 in the 18–57.6-pixel bins, rises strongly inward, and reaches 0.84 in the outermost
57.6–60-pixel bin. A masked 3.6-pixel-bin variance profile supplies an initial normalization diagnostic, not a
selected final estimator. See the [audit report](results/p4-step5-radial-pooling-20260918/README.md) for its exact
profile, endpoint policy, sample covariance summaries, and figure. The default Python sampler retains its original
behavior and agrees with production geometry at 16 calibration/evaluation centers spanning all eight test radii.

The subsequent saved-image comparison is recorded below. Once a revised policy is frozen, the tentative fresh
confirmation study remains 30 new sites at three transition brightnesses
(90 full reductions shared by every filter). The user intends to run that study on **ROC** and expects at least
a 3–4× speedup. That would put the local approximately ten-hour estimate near 2.5–3.3 hours, before analysis;
this is a planning estimate, not a measured ROC benchmark. Validate numerical consistency and throughput there
before launching the batch. No new injections or ROC jobs were launched for this baseline audit.

### Step 5 development result: radial pooling and normalization (2026-09-18)

The [completed comparison](results/p4-step5-radial-comparison-20260918/README.md) applies the four sampling/normalization
combinations to the saved baseline, 18 evaluation positives, and three original development positives. Radial
half-widths are 0, 5, 10, and 20 pixels, with PCA rank three and variance-floor fraction 0.1 fixed. An independent
NumPy/SciPy implementation reproduces the original same-radius production control before comparing extensions;
production filtering remains unchanged. These previously inspected images are development data.

| Training / half-width | Detections at 0.5× / 1× / 2× | Original evaluation null exceedances | Positive ±1-sigma coverage |
| --- | --- | ---: | ---: |
| Raw / 0 px | 2/6, 5/6, 6/6 | 1/28 | 1/18 |
| Raw / 5 px | 3/6, 6/6, 6/6 | 2/28 | 2/18 |
| Raw / 10 px | 3/6, 5/6, 6/6 | 1/28 | 0/18 |
| Raw / 20 px | 4/6, 6/6, 6/6 | 2/28 | 2/18 |
| Normalized / 0 px | 3/6, 5/6, 6/6 | 1/28 | 1/18 |
| Normalized / 5 px | 3/6, 5/6, 6/6 | 2/28 | 3/18 |
| Normalized / 10 px | 2/6, 5/6, 6/6 | 1/28 | 0/18 |
| Normalized / 20 px | 4/6, 6/6, 6/6 | 2/28 | 3/18 |
| Identity reference | 3/6, 6/6, 6/6 | 2/28 | — |
| Gaussian 3.6 px + application SNR reference | 3/6, 6/6, 6/6 | 2/28 | — |

All pooled policies make the eight originally invalid radius-20 null searches trainable. Those searches remain
outside the original 28+28 common comparison set. At the radius-12 development position, ±5 pixels meets the
eight-patch minimum at the center but fails at one search neighbor; ±10 pixels is needed for the full five-pixel
search. Every normalized measurement uses consistent data/template scaling and agrees with an independent solve
using the equivalent physical-unit covariance.

Both ±20-pixel variants recover one extra faint site relative to identity and Gaussian/SNR at the same observed
null count. **No general detection gain is established:** eight policies and six reused sites provide a development
hint, and radial normalization gives no consistent benefit. Conditional sigmas still fail badly, including at
unmodified null centers, where only 0–4 of 28 estimates lie within ±1 sigma of zero.

Covariance fits to disjoint native y half-planes were compared on the same 16 eligible centers for every policy.
Median physical-weight cosine decreases from 0.895 to 0.640 for raw widths 0 to 20, and from 0.892 to 0.693 with
normalization. More patches have not made these weights agree better. Spatial covariance variation and finite-sample
mode estimation remain possible explanations; shared radial profiles and spatial correlations make this a descriptive
diagnostic rather than independent cross-validation. Source effects on training also remain measurable.

The regularized model represents about 45% of centered sample covariance trace at one radius, falling to 28–30%
with ±20-pixel pooling in the respective fitting coordinates. This is not a true-variance estimate or an amplitude-sigma
rescaling factor. **Follow-on test (completed below):** vary the variance-floor fraction through 0.1, 0.3, and 1.0 on these saved images,
holding rank three and sampling/normalization fixed. Check conditional coverage and split-weight stability alongside
recovery against identity/Gaussian. Modeling discarded-mode variance is a possible separate control. Retain the
outer radial-profile rise and covariance transfer across radii as open questions, then freeze a policy before the
fresh ROC study.

Verification covers 320 production baseline pixels (maximum relative numerical difference 5.62e-8), all 18 original
positive measurements and decisions, 25 independent interpolation-ring comparisons, and a 34-ring forbidden-pixel
NaN mutation check. Original data and thresholds are fingerprinted; no new reductions ran. The report contains
individual records, the common-support stability table, a figure, and the reproduction command.

### Step 5 development result: variance-floor comparison (2026-09-18)

The [completed floor comparison](results/p4-step5-variance-floor-20260918/README.md) tests fractions 0.1, 0.3, and 1.0
with all eight existing sampling/normalization policies, retaining at most three PCA modes and the same exclusions,
responses, profiles, and comparison support. All 24 calibration thresholds were frozen before filtering the same
21 saved positive images. The recovery comparison still contains 18 positives at six sites and three brightnesses;
the three original development positives remain separate. No new reductions or ROC jobs ran.

**Higher floors modestly improve conditional coverage and weight stability, but do not resolve the uncertainty
mismatch.** Ranges below span the eight sampling policies; they are not confidence intervals or pooled counts.

| Floor fraction | Positive ±1-sigma coverage, out of 18 | Null-center ±1-sigma coverage, out of 28 | Median modeled / sample covariance trace across policies |
| --- | --- | --- | --- |
| 0.1 | 0–3 | 0–4 | 0.282–0.453 |
| 0.3 | 0–5 | 3–5 | 0.476–0.630 |
| 1.0 | 0–6 | 6–8 | 1.153–1.259 |

Only two decisions change relative to each policy's floor-0.1 result. At floor 1.0, normalized ±10-pixel pooling
adds a faint recovery at radius 26, block 3, going from 2/6 to 3/6 with 1/28 null exceedances. Raw ±20-pixel pooling
adds a null exceedance at radius 34, block 2, going from 2/28 to 3/28 with faint recovery unchanged at 4/6. Floor
0.3 changes no detection decisions. All middle- and high-brightness recovery counts remain unchanged. Normalized
±20 retains 4/6, 6/6, 6/6 recoveries and 2/28 null exceedances at every floor, compared with 3/6, 6/6, 6/6 and 2/28
for identity and Gaussian/SNR. This reused development sample still cannot establish a general gain or select a
preferred policy.

Conditional sigma increases monotonically at every valid fit; the median positive sigma ratio between floors
1.0 and 0.1 is 2.99–3.14 across policies. Empirical calibration thresholds change with the scores, explaining why
larger sigmas need not materially change recovery. Split-weight agreement improves modestly on the same 16 eligible
centers: raw ±20 cosine increases from 0.640 to 0.708 and normalized ±20 from 0.693 to 0.751. Wider pooling still
has poorer agreement between angular halves than same-radius training. These remain descriptive, spatially related
splits, with a shared variance profile for normalized fits.

Floor 1.0 now represents more total variance than the empirical training sample in the fitting coordinates, yet
coverage remains poor. Total trace does not determine variance along a filter's weights, nor establish that the
training covariance represents the candidate location. Mean estimation, omitted directional correlations, spatial
transfer, and interpolation remain possible contributors. This test does not isolate one cause or justify a global
sigma correction. Floor 1.0 retains the leading covariance modes; it is not identity filtering.

**Follow-on diagnostic (completed below):** project disjoint held-out noise stamps through frozen weights fitted to another angular subset.
Compare the projected mean and variance with zero and the model's predicted amplitude variance; separately audit
interpolation versus native candidate statistics and radial transfer. This should distinguish mean offsets from
covariance mismatch before another rank/floor sweep. Production defaults remain unchanged, and a fresh ROC study
still requires a frozen policy and numerical/throughput check.

The optional floor parameter preserves the original helper default of 0.1. Verification reproduces 14,265 archived
scalar values across the previous null/positive measurements, profiles, and summaries, plus 320 direct production
baseline checks. All 24 thresholds/counts and the common 16-site stability comparison were independently checked.
Covariance/sigma ordering, the isotropic limit, consistent physical-unit normalization, and unchanged fingerprints
also pass. See the report for full tables, individual records, and the figure.

### Step 5 development result: projected noise and discarded modes (2026-09-18)

The [projected-noise diagnostic](results/p4-step5-projected-noise-20260918/README.md) uses the saved baseline at the
same 16 common split sites (radii 46 and 50), in both angular directions, for all 24 existing policies. Training
and validation patches read disjoint native y half-planes, including across pooled rings. Primary validation always
uses the opposite-half **same-radius ring**, so wider training does not change the validation sample. This gives
32 directional fits per policy, with overlapping/reused patches rather than independent trials. Normalized fits
condition on the existing shared radial profile. No positive images, reductions, or new detection thresholds ran.

For frozen amplitude weights $g=C^{-1}t/(t^TC^{-1}t)$, compare the variance of
$z_j^{\rm proj}=g^T(x_j-\widehat\mu)/\sigma_\alpha$ with one, where $\sigma_\alpha^2=g^TCg$.
Separate the centered variance from squared mean error. Ranges below are ranges of policy medians over the same
32 directional fits, not confidence intervals.

| Floor fraction | Training projected variance / model | Held-out same-radius projected variance / model |
| --- | --- | --- |
| 0.1 | 17.18–29.57 | 17.88–36.23 |
| 0.3 | 5.72–9.78 | 5.94–12.08 |
| 1.0 | 1.71–2.89 | 1.85–3.63 |

**The covariance approximation underestimates variance along the filter weights even in its own training data.**
At floor 1.0, squared mean offsets account for only 1.3–11.1% of held-out MSE in the policy medians. Explicit
eigenspace accounting places 98.01–99.66% of empirical training projection variance outside the retained three
modes. The isotropic floor underestimates this complement contribution by factors 1.72–2.97. The retained
empirical eigenvalues are represented exactly, so the in-training variance discrepancy is in the modeled complement.
This identifies a limitation of the three-mode-plus-isotropic approximation without establishing that an
unregularized full sample covariance would generalize.

Paired bilinear/nearest-neighbor validation projections have median variance ratios 0.814–0.904 at floor 1.0.
The exact bilinear white-noise gain relative to native-grid sampling is 0.867–0.902 in the policy medians, including
shared native inputs. These controls show a smaller interpolation effect than the covariance mismatch; they do
not justify a universal correction for correlated residuals. Nearest-neighbor resampling is not an exact native
candidate stamp. The original 28 native evaluation-center scores from full-training fits retain variance 4.3–6.8
at floor 1.0, but are a different spatial sample and fitting setup.

Radial transfer also remains imperfect. For ±20-pixel training at floor 1.0, held-out rings five pixels inward
have 2.20 times the same-radius projected variance for raw pixels and 1.87 times after normalization (median paired
ratios on 24 eligible directional fits). Ten pixels inward the ratios are 5.15 and 3.86 on only eight fits. These
counts and paired comparisons prevent attributing a change in validation support to a covariance improvement.
The experiment does not establish inner-separation performance or a preferred production policy.

**Follow-on comparison (completed below):** retain all empirical modes using trace-preserving shrinkage,
$C_\gamma=(1-\gamma_{\rm shrink})\widehat C+
\gamma_{\rm shrink}\operatorname{tr}(\widehat C)I/p$, for $\gamma_{\rm shrink}=0.1,0.3,1.0$.
Use the same sampling and fixed directional splits first, with the three-mode controls. This tests whether
preserving more measured correlations improves held-out variance prediction while regularizing the sample
nullspace. The isotropic endpoint retains the fitted mean, so it need not reproduce the original identity
reference's background handling. Keep identity/Gaussian references for subsequent recovery comparisons. No
shrinkage experiment or production change is claimed at this checkpoint.

Verification reproduces 5,828 archived control values, checks 144 ring geometries and 768 directional variance
identities, and independently recomputes 3,072 projection groups. Opposite-half NaN mutations leave 576 training
matrices unchanged with the profile fixed; analytic white-noise gains agree with dense covariance calculations
on 32 real stencils. The report contains the variance budgets, individual projections, radial-transfer counts,
figure, and reproduction command.

### Step 5 development result: covariance shrinkage (2026-09-18)

The [shrinkage comparison](results/p4-step5-shrinkage-20260918/README.md) tests the specified strengths 0.1, 0.3,
and 1.0 with all eight sampling/normalization choices. It retains the same 16 sites at radii 46 and 50, both angular
directions, and opposite-half same-radius validation patches. The three-mode, floor-1 control is reproduced on
identical samples. Training/validation native support is disjoint; normalized fits retain the shared fixed radial
profile. The 32 directional fits per policy remain correlated development measurements, not independent trials.

**These shrinkage settings do not improve the held-out projection diagnostic.** Ranges below span policy medians,
not confidence intervals. Actual variance ratios are paired by site and direction to the three-mode, floor-1
control, so increasing reported sigma cannot by itself appear as reduced noise.

| Shrinkage strength | Training variance / prediction | Held-out variance / prediction | Actual held-out amplitude variance / control | Split-weight cosine |
| --- | --- | --- | --- | --- |
| 0.1 | 0.0055–0.142 | 10.51–25.64 | 0.985–2.070 | 0.310–0.736 |
| 0.3 | 0.025–0.363 | 3.70–8.47 | 0.987–1.648 | 0.429–0.759 |
| 1.0 | 3.48–4.92 | 2.98–3.72 | 0.985–1.085 | 1.000 |

At ±20 pixels, strength 0.1 roughly doubles actual held-out amplitude variance: paired median ratios are 2.017
for raw and 2.070 for normalized pixels. Strength 0.3 gives 1.614 and 1.648. These weights also agree less well
between halves than the three-mode control. Near-unity same-radius variance ratios do not establish a gain.
Strength 1.0 is stable by construction because its covariance is isotropic, but its conditional variance remains
underestimated. Its fitted mean and radial scaling distinguish it from the original identity detection reference.

The training halves contain only 8–56 patches for 121 pixels, with measured sample ranks 7–55. At strength 0.1,
98.45–99.99% of the squared weight norm lies in the empirical sample nullspace in the policy medians; at 0.3,
90.52–99.90% does. These are weight-norm fractions in fitting coordinates, not true variance fractions. The
inverse favors directions assigned a small floor, achieving tiny training projections while held-out noise still
occupies those directions. Preserving total covariance trace and positive definiteness has not prevented this
failure to generalize.

The prior three-mode diagnostic identified omitted directional covariance; the full-mode test identifies poorly
supported empirical directions. Neither establishes that covariance weighting cannot help, and other shrinkage
strengths remain untested. No shrinkage policy is promoted to production or selected for the ROC study. This
checkpoint uses one saved baseline and generates no new source-recovery results, reductions, or thresholds.

**Follow-on PSD comparison (completed below):** a structured covariance in a fixed spatial-frequency representation,
following the user's Welch/PSD motivation. Average patch power with explicit window/edge handling and a positive floor, avoiding learned
eigenvectors from the small patch sample. Test unit response, consistent interpolation/normalization, and held-out
projected variance on the same splits before interpreting recovery. Local stationarity and radial transfer remain
assumptions to check. Retain Gaussian and identity references when saved-image recovery comparisons resume.

Verification reproduces 51,136 archived control scalars across 256 three-mode fits. All 768 shrinkage fits satisfy
trace, positive-variance, unit-response, and covariance-projection identities; physical-unit and standardized
solutions agree. Isotropic weights and band-independent actual centered validation variance are checked explicitly.
Independent SVD solves agree within 8.78e-15 relative precision-vector error, and 2,304 projection groups and paired
variance/MSE summaries are independently recomputed. The report contains all policies, sample-support diagnostics,
individual projections, stability comparisons, and the figure.

### Step 5 development result: Welch-style spectral covariance (2026-09-18)

The [PSD comparison](results/p4-step5-welch-psd-20260918/README.md) averages windowed patch powers in a fixed
spatial-frequency representation. It crosses rectangular/Hann windows and 0.1/0.3 isotropic spectral mixtures
with all eight sampling/normalization policies, retaining the same 16 outer sites and both angular directions.
The primary held-out ring, exact native exclusions, and shared radial profile are unchanged. Both the three-mode,
floor-1 and fitted-mean isotropic controls are reproduced on identical samples. No new reductions or positive-source
recovery measurements are part of this diagnostic.

**This structured model improves the held-out noise diagnostic.** Ranges below span the eight sampling-policy
medians for each window/mixture. Actual amplitude variance ratios are paired by site and training direction, so
inflating predicted sigma alone cannot account for the measured noise reduction. These ranges are descriptive,
not confidence intervals or independent trials.

| Estimation window / isotropic mixture | Held-out variance / prediction | Actual variance / three-mode control | Actual variance / isotropic control | Split-weight cosine |
| --- | --- | --- | --- | --- |
| Rectangular / 0.1 | 0.731–0.906 | 0.782–0.925 | 0.787–0.835 | 0.946–0.966 |
| Rectangular / 0.3 | 0.863–1.054 | 0.781–0.942 | 0.806–0.844 | 0.953–0.973 |
| Hann / 0.1 | 0.644–0.866 | 0.835–0.909 | 0.809–0.858 | 0.970–0.987 |
| Hann / 0.3 | 0.767–1.008 | 0.834–0.913 | 0.822–0.873 | 0.975–0.989 |

The three-mode control has held-out variance/prediction medians 1.85–3.63 and split-weight cosines 0.708–0.904.
PSD training variance/prediction medians span 0.84–1.35, avoiding the large training/validation gap of full empirical
shrinkage. The isotropic control retains the fitted mean and optional radial scaling, so its comparison does not
replace the original identity or `filter.lpfGaussFW=3.6` Gaussian detection references.

The estimator subtracts the unwindowed across-patch mean, averages windowed powers with normalization
$(n-1)\sum W^2$, and rescales their spectral mean to the unwindowed mean pixel variance $\bar v$. It then uses
$\widehat P_\beta=(1-\beta_{\rm PSD})\widehat P+\beta_{\rm PSD}\bar v$. A 21×21 zero-padded FFT keeps every
linear lag of an 11×11 patch distinct; the reconstructed 121×121 finite lag covariance does not wrap opposite
patch edges together. Its diagonal is $\bar v$ and its minimum eigenvalue is at least $\beta_{\rm PSD}\bar v$.
The finite covariance is solved directly; its eigenvectors need not be 11×11 Fourier modes. At mixture one,
either estimation window reproduces the isotropic control exactly.

The window is used only to estimate covariance; candidate data and response remain untapered and consistently
radial-scaled where requested. No correction divides out the window autocorrelation or finite-patch lag-overlap
factor. This biases long-lag correlations downward, including for the rectangular window, and is an explicit
structural assumption rather than an unbiased covariance reconstruction.

**Calibration and radial transfer still need work.** Normalized ±20 pixels with rectangular mixture 0.3 has
same-radius median variance/prediction 1.00, but its directional 10th–90th percentile range is 0.38–1.19. Its
−10-pixel ring has 4.24 times the same-radius projected variance on the eight eligible paired fits; raw pixels
give 5.17. Normalization has not removed the radial-transfer problem. Including mean offsets, PSD/PCA MSE
ratios span 0.769–0.997 in policy medians; the normalized ±20-pixel Hann arms barely improve MSE despite reduced
centered variance. These reused, correlated outer sites and the shared profile do not establish inner-separation
performance, native-candidate calibration, or a preferred production policy.

**Next comparison:** apply the structured PSD family to the existing saved null and positive-injection images,
retaining common support and the Gaussian/identity references. Fix the candidate grid and threshold rule before
examining recovery; use only designated calibration nulls to set thresholds. Check contrast error, native-pixel
noise prediction, false positives, and recovery. This can reuse the existing reductions. It remains development;
freeze any eventual policy before fresh ROC injections.

Verification reproduces 101,692 archived control scalars and checks all 1,024 PSD fits, 512 windowed isotropic
endpoints, unit response, positive variance, and consistent physical-unit/standardized solutions. Independent
direct linear-lag sums and dense solves agree with every PSD covariance and held-out same-radius projection:
maximum relative covariance and weight differences are 7.97e-16 and 2.53e-15. Independent calculations also check
8,448 projection groups, 6,144 paired variance/MSE/sigma ratios, and all 512 stability pairs. Exact white ensembles
and constant patches verify normalization and edge lags. The report records the complete design, equations,
policy tables, radial-transfer support, limitations, provenance, figure, and reproduction command.

## 8. Notation

Dimensions refer to one local regression or one vectorized stamp, as indicated. Reused symbols are listed
separately by context; the P4 regression eigensystem and the residual-noise eigensystem describe different spaces.

| Symbol | Meaning | Notes |
| --- | --- | --- |
| $T$ | Number of temporal samples in a local P4 fit | A scalar count; the superscript $T$ on a matrix or vector denotes transpose instead. |
| $P$ | Number of predictor columns in a local P4 fit | Uppercase; distinct from the stamp-pixel count $p$. |
| $X$ | P4 predictor matrix | $T\times P$; the supported paired-refit path uses uncentered predictors. |
| $y$ | Target pixel's time series | Length $T$; distinct from the final-image data stamp $d$. |
| $U,S,V$ | Singular-value decomposition factors of $X$ | $X=USV^T$; columns of $U$ are temporal modes, columns of $V$ are predictor-space modes, and matrix $S$ contains singular values. |
| $k$ | Number of retained P4 modes | Defines the regression truncation; distinct from the covariance-model rank $r$. |
| $U_k,S_k,V_k$ | Retained singular-value decomposition factors | Contain the first $k$ modes, ordered by decreasing singular value. |
| $\beta_k$ | Fitted predictor coefficients | Length $P$; $\beta_k=V_kS_k^{-1}U_k^Ty$. Depends on the source amplitude when the regression is refitted. |
| $r_k$ | P4 residual time series | Length $T$; $r_k=y-X\beta_k$. The subscript identifies the retained mode count. |
| $q$ | Trial source's sky location | A position label; distinct from either use of $q_i$ below. |
| $\alpha$ | Source amplitude or contrast | Multiplies a unit-amplitude source model; physical contrast requires the stated PSF normalization. |
| $\epsilon$ | Positive half-amplitude of a central-difference probe | Corresponds to `psfResponse.refitContrast`; trials use $+\epsilon$ and $-\epsilon$. |
| $A_q$ | Unit-source contribution to the predictor matrix | $T\times P$; $X(\alpha)=X+\alpha A_q$. Includes the actual predictor sampling/interpolation. |
| $s_q$ | Unit-source contribution to the target time series | Length $T$ locally; $y(\alpha)=y+\alpha s_q$. In the full-Jacobian discussion, the same notation denotes the complete input source direction. |
| $h_{\mathrm{frozen},q}$ | Response with predictor coefficients held fixed | $s_q-A_q\beta_k$; omits coefficient adaptation. |
| $h_q$ | Complete first-order local P4 response | $r_k'(0)=(I-\Pi_k)s_q-\dot\Pi_k y$; the dependence on $k$ is implicit. Must be reconstructed and combined to form a final template. |
| $G$ | Temporal Gram matrix of the predictors | $T\times T$; $G=XX^T$. It is not the final-stamp noise covariance $C$. |
| $u_i,\lambda_i$ | Temporal Gram eigenvector and eigenvalue | $Gu_i=\lambda_i u_i$; $u_i$ has length $T$, and eigenvalues are ordered decreasingly. Positive $\lambda_i$ are squared singular values of $X$. |
| $i,j$ | Mode indices | In the projector derivative, $i\le k$ is retained and $j>k$ is discarded. Elsewhere they index the relevant covariance modes. |
| $\Pi_k$ | Orthogonal projector onto retained temporal modes | $T\times T$; $\Pi_k=U_kU_k^T$. |
| $\dot G,\dot\Pi_k$ | Source-amplitude derivatives of the Gram matrix and projector | Evaluated at $\alpha=0$ for the specified source direction; $\dot G=A_qX^T+XA_q^T$. |
| $U_a,U_b$ | Retained and discarded temporal eigenvector blocks | $U_a=U_k$; $U_b$ includes the discarded complement, including nullspace directions when needed. |
| $B$ | Retained-to-discarded mode-coupling matrix | $B_{ji}=u_j^T\dot G u_i/(\lambda_i-\lambda_j)$; its rows index discarded modes and columns retained modes. Requires a separated cutoff. |
| $X_c,y_c$ | Temporally centered predictors and target | Used by the centered extension; its fitted coefficients are still applied to uncentered predictors. |
| $p$ | Number of pixels in a vectorized stamp | Lowercase; vectors and covariances must use matching valid support and pixel ordering. |
| $d$ | Observed final-image data stamp | Length $p$; one candidate location and one reduction-mode plane. |
| $t$ | Processed unit-amplitude source template | Length $p$; includes the reduction response and uses the same coordinates as $d$. |
| $\widetilde t$ | Approximate processed template | Used to assess spatial averaging/interpolation error relative to $t$. |
| $\mu$ | Estimated mean background stamp | Length $p$; subtraction defines $z=d-\mu$. |
| $z$ | Mean-subtracted data stamp | Length $p$; modeled as $z=\alpha t+n$. |
| $n$ (noise vector) | Residual noise realization | Length $p$, with zero mean in the statistical model; distinct from scalar sample counts below. |
| $C$ | Residual-noise covariance in the filtered data space | $p\times p$; $C=\operatorname{Cov}(n)$. Assumed fixed and positive definite for the stated filter equations. |
| $N,E$ | Unweighted correlation and template energy | $N=t^Td$, $E=t^Tt$; the current filter returns $\widehat\alpha=N/E$. |
| $N_C,E_C$ | Covariance-weighted correlation and template energy | $N_C=t^TC^{-1}z$, $E_C=t^TC^{-1}t$; $E_C$ is conditional inverse amplitude variance under the adopted model. |
| $\widehat\alpha$ | Estimated source amplitude | $N_C/E_C$ with covariance weighting; includes its sign. |
| $\sigma_\alpha$ | Conditional standard error of the amplitude | $E_C^{-1/2}$; does not include uncertainty in the supplied template or estimated covariance. |
| $S$ (scalar score) | Signed theoretical matched-filter S/N | $N_C/\sqrt{E_C}$; distinct from the singular-value matrix $S$ and from the empirically calibrated `hciAnalyze` map. |
| $w$ | Covariance-weighted template | Length $p$; solve $Cw=t$, then evaluate $N_C=w^Tz$ and $E_C=w^Tt$. |
| $g$ | Unit-response amplitude weights | $g=C^{-1}t/(t^TC^{-1}t)$, so $g^Tt=1$ and the conditional amplitude variance is $g^TCg$. Uses the chosen raw or standardized coordinates. |
| $z_j^{\rm proj}$ | Standardized frozen-weight projection of noise patch $j$ | $g^T(x_j-\widehat\mu)/\sigma_\alpha$; its centered sample variance diagnoses measured/model amplitude-variance mismatch. |
| $A_{\rm samp}$ | Native-pixel sampling operator for a patch | Bilinear or nearest-neighbor extraction; $\lVert A_{\rm samp}^Tg\rVert^2/\lVert g\rVert^2$ is the projected white-noise variance gain relative to native-grid sampling. |
| $Q$ | Complete orthonormal noise eigenbasis | $p\times p$; $C=Q\operatorname{diag}(\nu_i)Q^T$. Distinct from P4's temporal basis $U$. |
| $q_i$ (noise eigenvector) | Column $i$ of $Q$ | Length $p$; used in the complete covariance-eigenbasis formulation. Distinct from scalar excess variance $q_i$ below. |
| $\nu_i$ | Noise variance along eigenvector $q_i$ | Positive eigenvalue of $C$; modal matched-filter weights use $1/\nu_i$. |
| $t_i,z_i$ | Template and data coefficients in the noise eigenbasis | $t_i=q_i^Tt$, $z_i=q_i^Tz$; retain their signs for coherent template detection. |
| $R$ | Mean-centered matrix of noise-training stamps | $n\times p$; each row is a stamp, and the sample covariance is $R^TR/(n-1)$. |
| $n$ (training count) | Number of covariance-training stamps | Scalar; the centered sample covariance has rank at most $n-1$. Overlap can reduce effective independence. |
| $\ell$ | Diffraction scale used in the annular sampling geometry | $\ell=\lambda/D$, expressed in image pixels. |
| $a$ | Radius of a covariance-training patch | A footprint of diameter $2a$ has half-width center spacing $\Delta s=a$. |
| $R_{\mathrm{ann}}$ | Radius of the annulus of patch centers | Measured from the star in image pixels; distinct from the training matrix $R$. |
| $\Delta s$ | Azimuthal arc spacing between patch centers | The half-overlap starting geometry uses $\Delta s=a$, giving approximately $2\pi R_{\mathrm{ann}}/a$ patches. |
| $b,\Delta R$ | Radial pooling half-width and center-ring spacing | Development grid: $b=0,5,10,20$ pixels and $\Delta R=5$ pixels; independent of the response footprint and P4 OR/SR dimensions. |
| $\rho$ | Stellar radius of an individual native image pixel | Used to normalize pixels within a patch, not just its center. |
| $\widehat v_{\mathrm{rad}}(\rho),\widehat\sigma_{\mathrm{rad}}(\rho)$ | Estimated radial variance and its square-root standard deviation | Estimated outside all declared source/held-out footprints; normalize by standard deviation. |
| $D_{\sigma,q}$ | Diagonal matrix of radial standard deviations at candidate $q$ | Positive scale factors at native stamp pixels; distinct from covariance-floor matrix $D$. |
| $\widehat C_{\mathrm{std}},\widehat\mu_{\mathrm{std}}$ | Covariance and mean of standardized training patches | Estimate and regularize after native variance normalization and patch interpolation; the covariance need not have unit diagonal. |
| $t_{\mathrm{std}},z_{\mathrm{std}}$ | Consistently standardized template and background-subtracted candidate data | $t_{\mathrm{std}}=D_{\sigma,q}^{-1}t$, $z_{\mathrm{std}}=D_{\sigma,q}^{-1}d-\widehat\mu_{\mathrm{std}}$; contrast retains its original units. |
| $x_j,\bar x$ | Vectorized training patch and mean training patch | Length $p$ in a common radial/tangential coordinate system; $\bar x$ is the mean over retained training patches. |
| $\widehat C$ | Empirical centered patch covariance | $p\times p$; regularize before inversion, and calibrate the effect of overlap on its estimation. |
| $\gamma_{\rm shrink}, C_\gamma$ | Shrinkage fraction and covariance that continuously reduces all empirical modes toward isotropy | $C_\gamma=(1-\gamma_{\rm shrink})\widehat C+\gamma_{\rm shrink}\operatorname{tr}(\widehat C)I/p$; tested grid 0.1, 0.3, 1.0. The endpoint 1.0 is isotropic. Distinct from Gram eigenvalues $\gamma_i$. |
| $W,U_W$ | Spectral estimation window and its squared energy | $U_W=\sum_a W_a^2$; tested windows are rectangular and separable symmetric Hann on 11×11 patches. Used inside the covariance estimator, not applied to candidate data or templates. |
| $\bar v$ | Mean unwindowed empirical pixel variance | $\sum_j\lVert x_j-\bar x\rVert^2/[(n-1)p]$; the PSD model preserves covariance trace $p\bar v$. Distinct from the radial variance profile. |
| $\widehat P_{\rm raw}(k),\widehat P(k)$ | Averaged windowed patch power before and after rescaling | Forward FFT is unnormalized on a 21×21 padded grid. Divide summed powers by $(n-1)U_W$, then rescale their spectral mean to $\bar v$; $k$ here indexes spatial frequency, not retained P4 modes. |
| $\beta_{\rm PSD},\widehat P_\beta$ | Isotropic spectral mixing fraction and regularized power | $\widehat P_\beta=(1-\beta_{\rm PSD})\widehat P+\beta_{\rm PSD}\bar v$; tested mixtures 0.1 and 0.3, with 1.0 checked as the isotropic endpoint. Distinct from P4 regression coefficients $\beta_k$. |
| $\widehat c(u_a-u_b),C_{\rm PSD}$ | Inverse-transform lag kernel and reconstructed finite stamp covariance | $u_a$ is the two-dimensional coordinate of stamp pixel $a$; $(C_{\rm PSD})_{ab}=\widehat c(u_a-u_b)$. The inverse FFT divides by $21^2$; all lags −10 through +10 remain distinct and long-lag window bias is retained. |
| $v_i,\gamma_i$ | Eigenvector and eigenvalue of the training Gram matrix | $RR^Tv_i=\gamma_i v_i$; for $\gamma_i>0$, $q_i=R^Tv_i/\sqrt{\gamma_i}$ and $\nu_i=\gamma_i/(n-1)$. |
| $C_s,C_n$ | Signal and noise covariance matrices in Wiener estimation | Defined for zero-mean, uncorrelated random signal and noise. $C_n$ is $C$ when referring to the same filtered stamp space. |
| $\widehat s$ | Wiener estimate of the random signal image | $\widehat s=C_s(C_s+C_n)^{-1}d$ for zero-mean data; distinct from the unit-source input $s_q$. |
| $v_\alpha$ | Prior variance of the random source amplitude | The rank-one signal model has $C_s=v_\alpha tt^T$; distinct from the estimator variance $\sigma_\alpha^2$. |
| $r$ (covariance rank) | Number of retained noise-covariance modes | Scalar; sets the low-rank model size, independently of the P4 subtraction count $k$. |
| $U_r$ | Retained residual-noise eigenmodes | $p\times r$ with orthonormal columns; despite its letter, it is not a block of P4's temporal basis $U$. |
| $q_i$ (excess variance) | Correlated variance above the isotropic noise floor | Nonnegative scalar in $C=\tau^2I+U_r\operatorname{diag}(q_i)U_r^T$; total retained-mode variance is $\tau^2+q_i$. |
| $f_{\rm floor}$ | Fraction multiplying the median empirical pixel variance to set the isotropic floor | $\tau_f^2=f_{\rm floor}\operatorname{median}_j(\widehat C_{jj})$; development grid 0.1, 0.3, 1.0. Applied in the selected raw or standardized fitting coordinates. |
| $\tau^2$ | Isotropic residual-noise variance floor | Strictly positive; assigns finite uncertainty to the complement of $U_r$. |
| $D$ | Diagonal residual-noise covariance floor | $p\times p$ with positive diagonal entries; allows spatially varying independent noise in $C=D+LL^T$. |
| $L$ | Low-rank covariance factor | $p\times r$; $L=U_r\operatorname{diag}(\sqrt{q_i})$. Distinct from the reduction operator $L_{\mathrm{red}}$. |
| $v$ | Arbitrary vector to which the precision operator is applied | Length $p$ in the expression for $C^{-1}v$; distinct from indexed Gram eigenvectors $v_i$. |
| $L_{\mathrm{red}}$ | Fixed linear reduction operator | Maps the full input data vector into the selected output space. |
| $J$ | Full data Jacobian of an adaptive reduction | First-order input-to-output map at a chosen baseline; individual source responses give only its action on selected directions. |
| $C_{\mathrm{in}},C_{\mathrm{out}}$ | Input and output noise covariances | Fixed-linear propagation gives $C_{\mathrm{out}}=L_{\mathrm{red}}C_{\mathrm{in}}L_{\mathrm{red}}^T$; adaptive first-order propagation uses $J$ instead. |
| $\eta$ | Covariance-weighted template similarity | Expected signed S/N using $\widetilde t$, relative to the exact-template optimum using $t$, for the same correct fixed covariance. |
| $r$ (radius) | Angular separation expressed in image pixels | Used in `hciAnalyze`'s small-sample correction; distinct from residual $r_k$ and covariance rank $r$. |
| $\lambda/D$ | Diffraction angular scale, expressed in image pixels for the correction | Here $\lambda$ is wavelength and $D$ telescope aperture diameter; distinct from Gram eigenvalues $\lambda_i$ and covariance-floor matrix $D$. |
| $n$ (comparison count) | Approximate number of comparison resolution elements | $n=2\pi r/(\lambda/D)-1$ in the small-sample correction; not the training-stamp count or a noise vector. |
| $I$ | Identity matrix | Dimension follows the space of the expression, such as $T$, $p$, or $r$. |
| Superscript $T$ | Matrix or vector transpose | All displayed models use real-valued data. |
| Dot or prime | Derivative with respect to source amplitude | For example, $\dot\Pi_k$ and $r_k'(0)$ refer to the specified source direction at $\alpha=0$. |
| $\operatorname{diag}(a_i)$ | Diagonal matrix with entries $a_i$ | Used for singular values, covariance eigenvalues, and modal weights. |
| $O(\cdot)$ | Asymptotic operation-count scaling | Describes leading computational work, not a measured wall-time speedup. |
