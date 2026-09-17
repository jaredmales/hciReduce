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
claim of optimality for full covariance estimation. Choose tapering separately and propagate any taper through
the data, template, and covariance consistently.

For candidate exclusion and validation, hold out whole angular neighborhoods, including every training patch
whose footprint intersects the candidate's exclusion region. Leaving out only the patch centered on the candidate
would still expose the covariance estimate to that source through overlapping neighbors. Compare the resulting
filter with identity weighting using held-out null scores, injected-source completeness, and amplitude uncertainty
calibration, as planned in Section 7.

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
6. **Consider time-resolved likelihoods after the final-image prototype.** Measure the information gained against
   storage and covariance-training costs before expanding the output format.

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

### Step 3 implementation checkpoint: shared factors and consistent cropping (2026-09-17)

`P4PCA::prepareResponse()` now owns an FP64 baseline eigensystem and its rank/cutoff diagnostics;
`calculateResponse()` can apply successive independently sampled source directions without another eigensolve.
The one-shot API delegates to these operations. The response calculation and boundary policy are unchanged.

The reduction groups source measurements into batches, gathers their detector requests, and factors each unique
valid detector fit once per batch. `psfResponse.analyticBatchSize` defaults to 32; 1 provides a comparison without
cross-source reuse. The existing memory budget limits both the source batch and worker count using conservative
geometry, residual, factor, and scratch estimates. Factors remain worker-local and are discarded after their
consuming source directions. FITS provenance records `P4 PSF ANALYTIC FACTOR COUNT` and the realized
`P4 PSF ANALYTIC BATCH SIZE`, separately from response-direction and fallback counts.

Automatic cropping now uses the same residual crop, derotation, and final crop with integrated filtering enabled
or disabled. All newly written P4 response methods publish final-image coordinates plus detector-origin cards.
The analysis script honors those offsets, including odd detector-to-output size differences. A 24-frame D64
AF Lep check with integrated filtering produced science bit-for-bit identical to the unfiltered baseline and
900 finite filtered pixels. The synthetic regression covers analytic, paired-refit, detector-local, and exact-sky
product coordinates and unchanged science; its thin annulus provides only one or two usable pixels in some
3-by-3 stamps, so the analytic/refit finite-output check explicitly permits that support fraction.

The experimental P4Reduction suite passes 35 cases / 354799 assertions, including identical products and
per-measurement diagnostic counts for batch sizes 1 and 32 with fewer baseline factorizations in the shared run.
The numerical response suite passes 11 cases / 1417 assertions, including independent derivative checks for
successive directions, both Gram orientations, owned baseline inputs, and rejection after failed preparation.

The mxlib coverage gate rechecked the current filtered LCOV trace. The called FP64 `eigenSYEVR` range is 60/60
executable lines, workspace construction/destruction 3/3 each, and workspace cleanup 13/13. The called mutable
float cube views, allocation/access, FITS, finite-check, and configuration ranges are covered. Exact const-float
cube views and float assignment instantiations remain the previously documented non-blocking upstream coverage
follow-ups in [mxlib_cleanup.md](mxlib_cleanup.md).

The full-data run exposed repeated nested cubic interpolation in source sampling. Each source batch now caches
shifted detector pixels once per source/frame; predictor sampling reuses those values with the exact original
column-then-row accumulation order. The memory estimate includes the cache. All 2400 published stamps in the
24-frame check are bit-for-bit identical to the pre-cache analytic products, and science remains unchanged.
The realized eight-source batch uses 424 baseline factorizations for 610 source directions per mode.

An additional broad PCA regression caught a refactoring change in overflow handling for unresolved modes.
The original error behavior is restored. Experimental suites pass 42 / 4677 (PCA), 11 / 1417 (response), and
35 / 354799 (reduction). Default suites pass 33 / 2220, 11 / 1417, and 28 / 78138, respectively;
`hciAnalyze` passes 13 / 166. These include the cached source path. No new mxlib call is introduced by caching.

Stack sampling of the full response run then identified dense FP64 products as the dominant active work.
For temporal-Gram source directions with at most 20% nonzero entries, the kernel now forms `A*X^T` using an
exact sparse representation and adds its transpose. Only exact zeros are omitted; there is no source-value
threshold or change to the response formula. Dense directions retain the original multiplication path.
Both the 24- and 96-frame products (4800 published stamps total) remain bit-for-bit identical to the earlier
analytic outputs; their largest per-stamp errors against the FP64 paired-refit oracle remain `3.15e-7` and
`1.96e-7`. The new independent sparse-direction test brings the response suite to 12 cases / 1440 assertions.
The interrupted full runs retain logs and frozen software for controlled comparisons; their partial timings
are not complete-run benchmarks.

Full-data validation is in progress under `/tmp/p4-step3-aflep`. The 621-frame science baseline took 241.43 s
and 5,170,444 KiB peak RSS with 20 OpenMP workers and one BLAS thread. Its maximum difference from the archived
post-avoidance science image is `1.61e-6` per pixel. Applying the archived post-avoidance response field to this
current baseline reproduces the accepted fit: contrast `0.004884882067` (change `-4.08e-11`), separation
`11.49649042` pixels (change `-3.88e-8`), PA `262.32662166` degrees, and SNR `4.32275218`.
This checks reference replay; it is separate from the same-build bitwise science-invariance requirement.

**Still required for Step 3:** complete the full analytic scientific comparison, measure reuse time/memory,
and broaden injections in position and brightness. No new photometric-calibration or speedup conclusion is
claimed at this checkpoint. Covariance estimation and noise weighting remain Step 4 and later work.

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
| $x_j,\bar x$ | Vectorized training patch and mean training patch | Length $p$ in a common radial/tangential coordinate system; $\bar x$ is the mean over retained training patches. |
| $\widehat C$ | Empirical centered patch covariance | $p\times p$; regularize before inversion, and calibrate the effect of overlap on its estimation. |
| $v_i,\gamma_i$ | Eigenvector and eigenvalue of the training Gram matrix | $RR^Tv_i=\gamma_i v_i$; for $\gamma_i>0$, $q_i=R^Tv_i/\sqrt{\gamma_i}$ and $\nu_i=\gamma_i/(n-1)$. |
| $C_s,C_n$ | Signal and noise covariance matrices in Wiener estimation | Defined for zero-mean, uncorrelated random signal and noise. $C_n$ is $C$ when referring to the same filtered stamp space. |
| $\widehat s$ | Wiener estimate of the random signal image | $\widehat s=C_s(C_s+C_n)^{-1}d$ for zero-mean data; distinct from the unit-source input $s_q$. |
| $v_\alpha$ | Prior variance of the random source amplitude | The rank-one signal model has $C_s=v_\alpha tt^T$; distinct from the estimator variance $\sigma_\alpha^2$. |
| $r$ (covariance rank) | Number of retained noise-covariance modes | Scalar; sets the low-rank model size, independently of the P4 subtraction count $k$. |
| $U_r$ | Retained residual-noise eigenmodes | $p\times r$ with orthonormal columns; despite its letter, it is not a block of P4's temporal basis $U$. |
| $q_i$ (excess variance) | Correlated variance above the isotropic noise floor | Nonnegative scalar in $C=\tau^2I+U_r\operatorname{diag}(q_i)U_r^T$; total retained-mode variance is $\tau^2+q_i$. |
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
