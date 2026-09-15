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

These checks verify the proposed algebra, not a production implementation or speedup. Existing production tests
were inspected, including the paired-local-reduction comparison in `tests/common/P4Reduction_test.cpp:2549` and
external P4 filtering in `tests/apps/hciAnalyze_test.cpp:349`. No production source was changed and no production
test suite or ROC dataset was rerun. Only this findings document was edited; the mxlib function-edit coverage gate
is therefore not triggered. Remaining follow-ups are the analytic-response benchmark and covariance
estimation/calibration experiments described above. The supplied discussion is preserved without changes to its
content in [PCA_Wiener_filtering_literature.md](PCA_Wiener_filtering_literature.md).

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
