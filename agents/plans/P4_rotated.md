Task: implement solving eqn:sigma_m_simple below as a new selectable version of P4.  The math is derived below using Latex.  Please review the math and provide any comments / corrections.  Then formulate a plan for how to implement this form of P4 in the section marked "#Plan" below.  Do not alter this prompt or the latex below, rather update the plan section with comments and corrections.  Do not begin implementation until I have approved the plan.

 
\section{P4}

For context here is the derivation of the algorithm we have implemented as P4:

\begin{eqnarray}
\sigma^2_n & = & \sum_t \left( x_t[n] - \hat{x}_t[n] \right)^2 \nonumber \\
           & = & \sum_t \left( x_t[n] - \sum_k c_{k,n} x_t[n + \Delta n_k] \right)^2 \label{eqn:sigma_n}
\end{eqnarray}

\begin{equation}
\frac{\partial\sigma^2_n}{\partial c_{j,n}} = -2\sum_t \left( \left( x_t[n] - \sum_k c_{k,n}  x_t[n + \Delta n_k] \right)
\left( x_t[n + \Delta n_j]    \right)\right)
\end{equation}

\begin{equation}
\frac{\partial\sigma^2_n}{\partial c_{j,n}} = -2 \left( \sum_t x_t[n]  x_t[n + \Delta n_j] -  \sum_k c_{k,n}\sum_t  x_t[n + \Delta n_k]  x_t[n + \Delta n_j]   \right)
\end{equation}

\begin{equation}
\sum_k c_{k,n}  \sum_t x_t[n + \Delta n_k]  x_t[n + \Delta n_j] = \sum_t x_t[n]  x_t[n + \Delta n_j]
\end{equation}

This is a system of linear equations of the form:
\begin{equation}
\bm{A}_n \bm{c}_n = \bm{b}_n
\label{eqn:p4_lsq_norot}
\end{equation}
where the matrix elements are
\begin{equation}
A_{jk,n } =  \sum_t x_t[n + \Delta n_j]  x_t[n + \Delta n_k]
\label{eqn:matrix_elements_norot}
\end{equation}
and the result vector is
\begin{equation}
b_{j,n} = \sum_t x_t[n]  x_t[n + \Delta n_j].
\label{eqn:projection_norot}
\end{equation}
The coefficients which minimize the residual error are then
\begin{equation}
 \bm{c}_n = \bm{A}_n^{-1} \bm{b}_n
 \label{eqn:p4_lsq_norot_solution}
\end{equation}

We can solve $\bm{A}_n^{-1}$ by eigendecomposition.

#####################################3 
Next Version:

Now we want to implement a cost function in the post-rotation images. 

\section{Rotation as a Linear Operation}

The value of pixel $m$ in the rotated image is determined by the dot product of the un-rotated image with the interpolation kernel
\begin{eqnarray}
x_t^r[m] &=& \bm{r}_{t,m} \cdot \bm{x_t} \nonumber \\
         &=& \sum_n r_{t,m}[n] x_t[n]
\end{eqnarray}
and for the noise estimate
\begin{eqnarray}
\hat{x}_t^r[m] &=& \sum_n r_{t,m}[n] \hat{x}_t[n] \nonumber \\
&=& \sum_n r_{t,m}[n] \sum_k c_{k,n} x_t[n+\Delta n_k]
\end{eqnarray}

\section{Time-axis Mean-Subtracted Variance }
\begin{equation}
\sigma_m^2 = \sum_t \left[ \left( x_t^r[m] - \frac{1}{T} \sum_t x_t^r[m] \right) - \left( \hat{x}_t^r[m] - \frac{1}{T} \sum_t \hat{x}_t^r[m] \right) \right]^2
\end{equation}
Substituting the expressions for the rotations we have
\begin{align}
\begin{split}
\sigma_m^2  = \sum_t  & \left[ \left( \sum_n  r_{t,m}[n] x_t[n] - \frac{1}{T} \sum_t \sum_n r_{t,m}[n] x_t[n]\right) - \right.\\
    & \left. \left(  \sum_n r_{t,m}[n] \sum_k c_{k,n} x_t[n+\Delta n_k] - \frac{1}{T} \sum_t  \sum_n \sum_k c_{k,n} r_{t,m}[n] x_t[n+\Delta n_k] \right) \right]^2 \\
    \end{split} \\
\end{align}
Rearranging and changing the order of summation yields
\begin{align}
\begin{split}
\label{eqn:sigma_m_full}
\sigma_m^2  = \sum_t  & \left[ \left( \sum_n \left(r_{t,m}[n] x_t[n] - \frac{1}{T} \sum_t r_{t,m}[n] x_t[n] \right) \right) - \right.\\
   & \qquad \qquad \qquad \left.  \sum_n \sum_k c_{k,n}  \left( r_{t,m}[n] x_t[n+\Delta n_k] - \frac{1}{T} \sum_t r_{t,m}[n] x_t[n+\Delta n_k] \right) \right]^2  \\
    \end{split}
\end{align}
We now see a problem: the $c_{k,n}$ can not be brought outside the sum over $n$.  In order to arrive at a system of equations linear in the coefficients we must therefore make a modification.  Rather than have a projection vector unique to each pixel $n$ in the pupil-fixed images, we instead can use a projection vector unique to each pixel $m$ in the rotated images.  This is a rather significant change...

Next we introduce an expression for a general mean-subtracted rotated pixel time-series
\begin{equation}
\tilde{x}_t^r[m;\mathcal{\eta}] =  \sum_n \left(r_{t,m}[n]x[\mathcal{\eta}] - \frac{1}{T}\sum_t r_{t,m}[n]x[\mathcal{\eta}]\right)
\end{equation}
with parameter $\mathcal{\eta}$. Equation (\ref{eqn:sigma_m_full}) can now be simplified to
\begin{equation}
\sigma_m^2 = \sum_t \left( \tilde{x}_t^r[m;n] -  \sum_k c_{k,m} \tilde{x}_t^r[m;n+\Delta n_k]  \right)^2 
\label{eqn:sigma_m_simple}
\end{equation}


# Plan:

## Mathematical review and corrections

The derivation motivates a useful new P4 variant, but the following points must be made explicit before implementation.
The equations above are left unchanged as requested; the notation below is the proposed precise contract.

1. The model with the original coefficients \(c_{k,n}\) is still linear in those coefficients.  The difficulty after
   rotation is that it is no longer a small, independent \(K\)-coefficient fit at each pixel: all coefficients whose
   source pixels contribute to output pixel \(m\) become coupled in one much larger system.  Replacing \(c_{k,n}\)
   with \(c_{k,m}\) is therefore a new coefficient-tying assumption that restores an independent \(K\)-coefficient
   problem at each rotated pixel; it is not required to make the problem linear.

2. The time index in every temporal mean must be a new dummy index, for example \(\tau\), rather than reusing the
   outer index \(t\).  The proposed general series also needs the frame index and the source-pixel argument on
   \(x\).  Once the sum over \(n\) has been performed, \(n\) cannot remain as the parameter after the semicolon in
   Equation (\ref{eqn:sigma_m_simple}).

3. Let \(S_k\) denote the complete two-dimensional predictor-offset operation, including the current radial/tangential
   mapping and its interpolation, and let \(S_0\) be the identity.  A precise version of the quantities implied by
   the current derivation is
   \[
      u_{t,m}^{(k)} = (R_t S_k x_t)[m]
                     = \sum_n r_{t,m}[n]\,(S_k x_t)[n],
   \]
   with
   \[
      \bar u_m^{(k)} = \frac{1}{T}\sum_\tau u_{\tau,m}^{(k)}, \qquad
      \widetilde u_{t,m}^{(k)} = u_{t,m}^{(k)}-\bar u_m^{(k)}.
   \]
   The target is \(y_m(t)=\widetilde u_{t,m}^{(0)}\), and the design matrix is
   \(D_m(t,k)=\widetilde u_{t,m}^{(k)}\) for \(k=1,\ldots,K\).  The reviewed objective and normal equations are then
   \[
      \sigma_m^2 = \lVert y_m-D_m c_m\rVert_2^2,
      \qquad A_m=D_m^T D_m,
      \qquad b_m=D_m^T y_m.
   \]
   Element-wise,
   \[
      A_{jk,m}=\sum_t \widetilde u_{t,m}^{(j)}\widetilde u_{t,m}^{(k)},
      \qquad b_{j,m}=\sum_t \widetilde u_{t,m}^{(j)}\widetilde u_{t,m}^{(0)}.
   \]
   As in the existing P4 implementation, the numerical solution is a truncated pseudoinverse, not a literal
   unrestricted \(A_m^{-1}\).

4. Temporal centering imposes one structural null direction.  The maximum possible degrees of freedom are therefore
   \(\min(K,T-1)\), not \(\min(K,T)\).  The realized mode counts, numerical-rank checks, diagnostics, and FITS
   provenance must all use that limit.  The implementation must cap the reported numerical rank at this structural
   limit even when roundoff makes the nominally zero temporal eigenvalue slightly positive, including when
   `p4.rankTolerance=0`.

5. There are two scientifically distinct meanings of a "post-rotation predictor":

   - **shift then rotate**, as written above: \(u^{(k)}=(R_t S_k x_t)[m]\); or
   - **rotate then shift**: \(u^{(k)}=(S_k R_t x_t)[m]\), equivalently a fixed neighbor of \(m\) in each derotated
     image.

   These commute for ideal continuous rigid rotations only when the offset field is rotation-equivariant.  They are
   not identical for the actual discrete interpolation operators, masks, and image boundaries.

   **Decision:** use **rotate then shift** for the rotated variant.  Construct the OR wedge once in the post-rotation
   sky frame.  For output search pixel \(m\), let its fixed predictor coordinates be
   \(q_{m,k}=m\mathbin{\oplus}\Delta m_k\).  In frame \(t\), the detector coordinate supplying that predictor is
   \[
      p_{t,m,k}=C+Q_t^{-1}(q_{m,k}-C),
   \]
   where \(C\) is the rotation center and \(Q_t\) is the configured detector-to-sky derotation.  Thus
   \(u_{t,m}^{(k)}=(R_tx_t)[q_{m,k}]\), and the raw source of every fixed sky predictor is tracked through the inverse
   rotation.  In this variant the old \(\Delta n_k\) notation should be read as a sky-frame \(\Delta m_k\); the
   radial/tangential OR geometry rotates back into each detector frame automatically.

   **Implementation decision:** a lazy/direct sampler maps each \(q_{m,k}\) back to \(p_{t,m,k}\) and applies exactly
   one cubic interpolation in the original frame.  This is the most direct discrete realization of the approved
   operator and avoids the second smoothing operation introduced by first materializing an integer-grid rotated
   cube and then interpolating its fractional OR coordinates.  The materialized route remains a deliberately small
   test oracle, not a selectable production implementation.  Identity and exact-cardinal cases establish coordinate
   and sign agreement; fractional-angle tests against analytic or oversampled truth quantify the materialized
   route's additional smoothing, covariance change, validity footprint, throughput, runtime, and memory.  The
   rejected shift-then-rotate form remains useful as a small reference test proving that the operators differ.

6. The centered objective does not by itself specify what should be written as the science residual.  Three outputs
   are possible, and they do not preserve astrophysical flux in the same way:

   - \(u^{(0)}-D_m c_m\) subtracts the centered prediction from the uncentered target and preserves the target's
     temporal mean exactly;
   - \(u^{(0)}-U_m c_m\) subtracts an uncentered prediction and changes the mean by
     \(c_m^T\bar u_m\); and
   - \(y_m-D_m c_m\), or equivalently fitting and removing an intercept, forces a zero-mean residual and would remove
     a stationary companion's mean signal.

   **Decision (2026-08-18, superseding the mean-preserving rule):** estimate \(c_m\) from the centered objective and write
   \(u^{(0)}-U_m c_m\).  In the sky-aligned frame a stationary companion is absent from the variance being minimized,
   so it does not drive the coefficients.  The coefficient vector is nevertheless applied to the complete predictor
   values after inherited preprocessing and direct rotated sampling.  The residual mean is therefore
   \(\overline{u^{(0)}}-\bar u_m^T c_m\), not necessarily the target mean.  The adopted baseline convention is the
   raw post-preprocessing predictor DC together with the truncated minimum-norm coefficient solution.  Target-DC and
   predictor-DC tests lock the distinction, while injected-companion validation remains required to quantify
   throughput rather than assuming it from centering alone.

7. The same frame set must be used for the target and every predictor in one local fit.  The first implementation
   should require all \(T\) frames to have finite target values and complete target, predictor, and composite
   interpolation support; otherwise that output pixel is invalid in every mode plane.  Silently dropping a different
   set of frames per predictor would change the centering operator and the least-squares problem.  A future
   missing-frame implementation may use an explicit common valid set \(V_m\), with \(T_m=|V_m|\), but then mode
   realization and provenance must use \(\min(K,T_m-1)\).

8. The uncentered predictor application uses the predictor baseline remaining after inherited preprocessing.
   hciReduce provides two relevant operations that act before rotation, while a companion still moves through
   detector coordinates:

   - `preProcess.subradprof=true` subtracts a median radial profile independently from every detector-frame image,
     removing much of the axisymmetric stellar halo while generally leaving a localized companion; and
   - `preProcess.meanSubMethod=medianImage` forms the temporal median detector-frame image and subtracts it from every
     frame, providing a direct median-PSF subtraction.  A rotating companion is rejected from that median only to the
     extent allowed by the number of frames, field rotation, PSF width, and masking.

   These operations are distinct from the internal temporal centering of the rotated regression.  They modify the
   input cube from which both the uncentered science target and the centered fit are constructed; the internal
   centering then prevents the surviving stationary companion signal from driving the fit.  The two preprocessing
   operations may be combined in their existing order (per-frame radial-profile subtraction before temporal median
   image subtraction), but neither should be forced by rotated mode.  Their choices and order must remain in
   provenance, and injected-source tests must quantify any radial-profile bias or median-reference self-subtraction.

## Proposed user-visible contract

- Add one selector, proposed as `p4.regressionFrame=detector|rotated`, with `detector` as the default so all existing
  configurations and products remain unchanged.  Both `basic` and `normal` application modes continue to dispatch
  the selected P4 algorithm; this is not a new executable or application mode.
- Record the selected definition in configuration help, the configuration dictionary, diagnostics, progress output,
  and `HIERARCH P4 FRAME = detector|rotated` in FITS products.  The `rotated` definition includes application to the
  uncentered predictors.
- In `rotated` mode, temporal centering is part of the algorithm and is not an additional configurable preprocessing
  switch.  Record that fact in provenance.
- Continue to use the inherited detector-frame preprocessing configuration.  The User's Guide should recommend
  `preProcess.subradprof`, `preProcess.meanSubMethod=medianImage`, or both when the applied predictor baseline should
  first be reduced, while making clear that these are optional scientific choices with throughput consequences.  An
  already preprocessed sequence may continue to use `preProcess.skip=true` only when its preprocessing provenance is
  known.
- The rotated-mode science cubes are already in the final sky frame.  The final-processing lifecycle must combine
  them without applying ADI derotation a second time.  The implementation should carry an explicit data-frame state
  into shared final processing instead of relying on the user to set `combine.noDerotate` correctly.
- Initially reject `adi.postMedSub=true` in rotated mode.  It would impose a second sky-fixed baseline rule after the
  uncentered predictor subtraction and confound its throughput contract.
- Retain the existing annulus, optimization-region, exclusion-policy, rank-tolerance, diagnostic, combination, and
  output options.  `geom.minRadius`/`geom.maxRadius` describe the retained sky-frame output annuli in rotated mode.
- Preserve the existing detector-frame implementation as a separate tested branch.  Do not refactor its numerical
  path merely to share code unless exact regression tests show its products are unchanged.

## Implementation phases

Checklist status: `[x]` complete, `[~]` partially complete, `[ ]` not yet complete.

### Current implementation evidence (2026-08-18)

- The selectable CPU reference implementation is present: `detector` remains the default, while `rotated` uses
  fixed sky-frame predictors, direct inverse-mapped one-stage cubic sampling, all-frame support validation, centered
  PCA with the \(\min(K,T-1)\) limit, and uncentered predictor application in an explicitly sky-frame output
  lifecycle.
- Focused Release verification recorded 413 assertions in 14 P4PCA test cases, 1,337 assertions in 7 P4RotatedGrid
  test cases, and 43,559 assertions in 11 P4Reduction test cases.  The complete Release suite passed all 21 CTests,
  and the documentation target built with only the previously known unrelated warnings.
- A fresh `_build_fresh` `coverage` target passed all 21 CTests and generated the product-only report at 4,796/4,856
  executable lines (98.8%) and 256/257 filtered logical functions (99.6%).  The aggregate misses are confined to
  pre-existing `hciAnalyze.cpp` and `HCIobservation.hpp` paths.
- The current P4 LCOV snapshot reports 1,415/1,415 executable lines across `P4PCA.cpp`, `P4PixelGrid.cpp`,
  `P4RotatedGrid.cpp`, and `P4Reduction.cpp`; `P4PCA.cpp` itself is at 162/162 lines and 9/9 filtered logical
  functions.  Earlier sanitizer results predate the uncentered-application change and are not counted as current
  evidence.
- The direct sampler agrees with analytic/image-rotation references at tested coordinates.  On a fractional-angle,
  fractional-predictor oscillatory field it differs measurably from the materialized two-stage route and has lower
  aggregate squared interpolation error, supporting the one-stage production choice.
- The exact mxlib specializations newly used for cubic sampling, rotation, FITS image/header output, and finite-value
  checks are covered.  The remaining exact-specialization gap is `eigenCube<float>::median()` for the unmasked and
  validity-cube/masked `ADIobservation<float>::finalProcess()` paths; it is recorded under `Known non-blocking
  ownership follow-ups` in `agents/plans/mxlib_cleanup.md`.
- The AF Lep detector-versus-rotated comparison failed the scientific acceptance gate.  The implementation and its
  numerical tests remain useful as a reference for the specified centered-fit, uncentered-application rule, but the
  rotated variant must not be described or used as scientifically ready.  Injected-source throughput and astrometry,
  the rejected-operator oracle, and a controlled performance benchmark also remain incomplete.

### AF Lep scientific outcome (2026-08-18)

The 621-frame AF Lep/NACO sequence provides a controlled comparison because `finim_0013.fits` and
`finim_0014.fits` used the same radial-profile preprocessing, P4 geometry, mode fractions, and input frames; only
the regression frame changed.  `finim_0015.fits` tested rotated regression after detector-frame `meanImage`
subtraction.  In the science annulus \(8 \le r < 20\):

| Product | Regression/preprocessing | RMS at the 0.25 mode fraction | Best `hciAnalyze` SNR |
|---|---|---:|---:|
| `finim_0013.fits` | detector / radial profile | 0.220 | 4.79 |
| `finim_0014.fits` | rotated / radial profile | 18.121 | 2.28 |
| `finim_0015.fits` | rotated / `meanImage` | 6.323 | 2.26 |

The rotated output planes were nearly independent of mode count.  For the radial-profile case, the first and last
planes had correlation 0.999997 and differed by only 0.046 RMS in the annulus.  Independent target-only
reconstructions applied the same preprocessing, direct inverse-mapped cubic sampling, and iterative five-sigma mean
without P4.  They matched the last rotated planes as follows:

| Preprocessing | Target-only RMS | Rotated-P4 RMS | Correlation | Difference RMS |
|---|---:|---:|---:|---:|
| radial profile | 18.12296 | 18.12145 | 0.9999970 | 0.04399 |
| `meanImage` | 6.32310 | 6.32324 | 0.9999976 | 0.01406 |

This rules out rank failure, invalid support, angle sign, double derotation, preprocessing order, and final-combine
ordering as primary causes.  The outputs are effectively the no-P4 sky-aligned target stacks: the fitted raw
predictor correction has only about 0.2% of the target-stack RMS.  The implementation correctly minimizes

\[
    \left\|H\left(u_m^{(0)}-U_m c_m\right)\right\|_2^2,
\]

but that objective does not constrain the displayed residual mean

\[
    \overline{u_m^{(0)}}-\bar u_m^T c_m.
\]

Detector-frame `meanImage` subtraction does not make every fixed-sky time series zero mean because temporal
averaging and frame-dependent rotation do not commute.  A free or subsequently removed sky-frame intercept would
also remove a stationary companion, and the full coupled centered formulation has the same DC identifiability
problem.  A future scientifically useful formulation therefore needs a separately justified stellar-baseline model
or prior in addition to the centered dynamic predictor.  No replacement baseline rule or next implementation is
selected by this outcome.

The input sequence contains the rapid-rotation coaddition interval documented in `P4_rotated_full.md`, so these
products are not absolute golden data and the eventual accepted model must be retested on corrected coadds.  That
caveat does not invalidate the paired diagnosis: `finim_0013.fits` and `finim_0014.fits` use the same coadds and
preprocessing, so their direct comparison still isolates the regression-frame change.

### Phase 0: freeze the scientific definition

- [x] Use **rotate then shift** with a fixed sky-frame OR wedge.  For each frame, track every target/predictor sample
      back through the inverse derotation to its detector-frame source coordinate.
- [x] Use lazy/direct one-stage interpolation for production.  Map each fixed sky target/predictor coordinate through
      the inverse derotation and apply exactly one cubic interpolation in the original detector frame.  Retain a
      small materialized-two-stage oracle only for identity/cardinal/fractional comparison tests.
- [x] Use \(u^{(0)}-U_m c_m\): fit only centered temporal fluctuations, then apply the resulting coefficient vector
      to the original post-preprocessing predictor values.
- [x] Decide the FITS keyword names for regression frame, temporal centering, structural DOF, and residual-mean policy.
ANSWER: we can use long keywords and rely on HIERARCH.
 - P4 FRAME = detector  or rotated are good
 - Temporal centering: I think it's implicit if rotated
 - structural DOF: We report the mode count output, not the fraction input, so I think the K-1 issue is a detail here.  Document it, but it doesn't need to be in FITS heqder
 - residual mean: I think we answered this, nothing to specify

SUPERSEDING DECISION (2026-08-18): continue to write `HIERARCH P4 FRAME = detector|rotated`.  The centered-fit,
uncentered-application rule is intrinsic to the `rotated` definition and does not receive a separate compatibility
keyword.  The selected realized mode counts remain the numerical provenance.
Centering removes one temporal degree of freedom, so the limit is \(\min(K,T-1)\), not generally \(K-1\); if \(K<T\),
centering need not reduce the predictor-limited maximum.

- [x] Specify whether rotated mode initially requires every frame or may use a common per-pixel valid-frame subset.
      The recommendation is to require every frame for the first version.
ANSWER: agree

DECISION: every target and predictor in a retained local fit must have complete valid support in all \(T\) frames.
Any failure invalidates that output pixel for every mode plane; there is no frame dropping or per-predictor frame set
in the initial implementation.

The output-mean proposal below was approved and implemented on 2026-08-18.  The separate mask-validity proposal
remains unapproved and unimplemented.

- [x] **APPROVED — centered fit with an uncentered predictor subtraction.**  Estimate \(c_m\) from
      \[
          \underset{c_m}{\operatorname{minimize}}\;
          \left\|H\left(u_m^{(0)}-U_m c_m\right)\right\|_2^2,
      \]
      so a sky-stationary companion is absent from the temporal variance that determines the coefficients, but then
      interpret \(U_m c_m\) as the complete predictor-derived PSF estimate and write
      \[
          r_m = u_m^{(0)}-U_m c_m.
      \]
      The centered objective is equivalently an uncentered regression with a free intercept,
      \[
          \underset{a_m,c_m}{\operatorname{minimize}}\;
          \left\|u_m^{(0)}-a_m\mathbf{1}-U_m c_m\right\|_2^2,
          \qquad
          \widehat a_m=\overline{u_m^{(0)}}-\overline{u_m}^{\,T}c_m.
      \]
      The residual retains this fitted intercept; it is not the zero-mean residual
      \(H(u_m^{(0)}-U_m c_m)\).

      This proposal addresses an exact invariant of the implemented rule.  Since \(D_m=HU_m\),
      \[
          \overline{u_m^{(0)}-D_m c_m}=\overline{u_m^{(0)}}
      \]
      for every coefficient vector and every requested mode count.  An ordinary temporal mean combine therefore
      cannot remove or change the target's static PSF pedestal; mode planes can differ only through nonlinear steps
      such as clipping or changing validity.  The AF Lep result exposed this invariant as nearly identical output
      planes dominated by the derotated input mean.

      This rule also exposes a DC-level ambiguity that the centered objective cannot resolve.  Replacing a
      predictor by \(U_m+\mathbf{1}b^T\) leaves \(HU_m\), the objective, and the fitted \(c_m\) unchanged, but changes
      the residual mean by \(-b^Tc_m\).  The adopted convention is the raw post-preprocessing predictor zero point
      with the truncated minimum-norm solution, intrinsic to `HIERARCH P4 FRAME = rotated`.

      Required acceptance tests include a dense small-problem oracle; invariance to a constant added to the target;
      an explicit predictor-DC perturbation test that demonstrates and then locks the chosen baseline convention; a
      synthetic time-variable PSF with a known nonzero mean; and injected companions spanning separation, position
      angle, contrast, field rotation, and mode count.  Report residual variance, mode-to-mode differences,
      photometric throughput, astrometric bias, and behavior with no baseline subtraction, radial-profile
      subtraction, median-image subtraction, and their documented composition.  The dense oracle, target-DC,
      predictor-DC, and nonzero-mean synthetic tests now cover the numerical rule; the injected-source matrix remains
      a scientific acceptance gate.

- [ ] **REOPENED PROPOSAL — mask-aware whole-predictor pruning (NOT IMPLEMENTED).**  Continue to require the target
      sample at sky pixel \(m\) to have finite values and complete interpolation/mask support in every one of the
      common \(T\) frames.  Instead of invalidating the entire fit when any predictor fails support, retain predictor
      column \(k\) only when that predictor has finite values and complete support in all \(T\) frames, and otherwise
      remove that whole column before centering and solving.  This preserves one common frame set and one centering
      operator; it does not permit per-frame or per-predictor frame dropping.  The resulting local design has
      \(K_m\) retained columns and structural limit \(\min(K_m,T-1)\).  If no useful predictor columns remain, the
      output pixel is invalid.

      Before implementation, freeze how configured mode fractions/counts map onto a spatially varying \(K_m\).  The
      recommended contract is to give every output plane one global absolute requested mode count, use that count at
      a pixel only when it does not exceed the local numerical rank and \(\min(K_m,T-1)\), and otherwise mark that
      pixel invalid in that plane.  Silently re-realizing the same fraction against each local \(K_m\) would give one
      plane different model complexity at different pixels.  FITS/diagnostic provenance must expose the nominal
      predictor count, the distribution or map of \(K_m\), local rank failures, per-mode validity counts, and the
      eventual global-versus-local mode policy.  The exact aggregate keywords and whether a \(K_m\) diagnostic image
      is required remain Phase 0 decisions.

Phase 0 remains reopened only for mask-aware predictor pruning.  The uncentered predictor output rule is now part of
the implemented shortcut rotated contract.

### Phase 1: exact operator and numerical reference

- [~] Add a small pure reference component that constructs \(u_{t,m}^{(k)}\), centers target and predictors in double
      precision, forms the reviewed objective, and compares its result with a direct QR/SVD least-squares oracle.
      The centered PCA seam and an independent end-to-end thin-SVD oracle now cover this numerically, but there is no
      separate pure reference component that exposes the complete objective as a reusable test helper.
- [x] Represent offsets as checked two-dimensional coordinates or mapping operations.  Do not use flattened
      `n + delta` arithmetic that could wrap across image rows.
- [~] Add tests that show the analytic normal equations agree with finite differences of \(\sigma_m^2\), and that the
      truncated eigensolution agrees with direct SVD without comparing arbitrary eigenvector signs.
      Direct-SVD comparisons are present for both Gram-matrix branches; the finite-difference objective check remains.
- [x] Extend the pure PCA seam with an explicit centered/structural-rank path.  Preserve the current uncentered API and
      its behavior.  Test \(T=1\) rejection, \(T=2\), \(K<T\), \(K\ge T\), exact and near rank deficiency,
      non-contiguous requested modes, `rankTolerance=0`, and workspace reuse while switching paths.
- [x] Prove in tests that adding a constant to the target does not change fitted coefficients and passes that constant
      through exactly; prove separately that predictor DC leaves the centered fit unchanged but changes residual DC.
- [ ] Verify the centered fit and uncentered application after each supported baseline treatment: no preprocessing,
      per-frame radial-profile subtraction, detector-frame median-image subtraction, and their documented
      composition.

### Phase 2: rotation and geometry component

- [x] Extract or add a reusable rotation-sampling record that exactly matches `ADIobservation::derotate()` and
      `mx::improc::imageRotate` for sign, center, anchor, cubic weights, and complete nominal support.  Avoid a second
      handwritten rotation convention.
- [x] Use `P4PixelGrid`'s mapped predictor coordinates, or extract their geometry-only construction, to define the
      fixed \(q_{m,k}\) OR wedge in the sky frame without tying it to one image's stored interpolation kernel.
- [x] Implement the lazy/direct production path: precompute each frame's inverse rotation, map fixed \(q_{m,k}\) to
      detector coordinate \(p_{t,m,k}\), and sample the original frame once.  Generate bounded per-worker or per-tile
      interpolation records rather than materializing the prohibitive full \((m,t,k)\) record tensor.
- [~] Implement only a tiny materialized reference helper in the test harness: derotate a small image onto an integer
      sky grid and resample its OR points.  Compare it with the direct sampler at identity, cardinal, and fractional
      rotations; do not add a materialized production backend or configuration value.
      The fractional-angle/fractional-predictor comparison is present and the production backend remains direct-only;
      explicit identity and cardinal materialized-OR comparisons remain.
- [ ] Retain a tiny explicit **shift then rotate** \(R_tS_k\) oracle only to demonstrate the operator distinction; do
      not expose it as a third production configuration in this tranche.
- [x] Apply the common input mask and image-edge rules to every source pixel touched by the full composite operator.
      Negative or zero cubic coefficients do not relax nominal-footprint validity.
- [~] Test odd/even image sizes, explicit centers, zero/cardinal/fractional/wrapped angles, annulus boundaries, central
      masks, all four exclusion/mask states, exact image-edge acceptance/rejection, and transactional reconfiguration.
      Odd/even geometry, explicit centers, zero/cardinal/fractional angles, both exclusion policies, all-frame mask and
      edge invalidation, and transactional errors are covered; the remaining matrix includes explicit wrapped-angle
      and any still-unrepresented annulus/mask-state combinations.

### Phase 3: selectable P4 orchestration

- [x] Add the regression-frame enum, configuration registration/loading/validation, explicit default, CLI help, and
      parser tests.
- [x] Dispatch the existing detector path unchanged and add a rotated path that assembles one centered design matrix
      per retained sky pixel, calls the centered PCA operation, and writes already-rotated residual/validity cubes.
- [x] Realize rotated-mode fractions against \(\min(K,T-1)\); continue to reject zero or duplicate realized mode
      counts.  Rank-insufficient planes remain invalid/NaN with aggregated counts and warnings.
- [x] Generalize shared final processing with an explicit detector-frame versus sky-frame state.  Preserve the tested
      KLIP and detector-P4 lifecycle; sky-frame P4 skips only derotation and still uses validity-aware combination,
      FITS output, diagnostics, and error propagation.
- [x] Include the selected operator ordering, residual policy, structural DOF, realized modes/rank, validity counts, and
      data frame in FITS provenance and diagnostic summaries.
      Provenance uses `HIERARCH P4 FRAME`; centering, uncentered predictor application, and structural DOF are intrinsic
      to the documented `rotated` definition, while realized modes/rank and validity counts remain explicit.
- [~] Extend the existing annulus-local/overall progress output to distinguish rotation preparation, rotated feature
      assembly, eigensolve/projection, and final combination.
      Rotated geometry, annulus-local, overall, and final-completion progress are present; separate feature-assembly,
      eigensolve/projection, and final-combination stage records are not yet present.

### Phase 4: verification and performance gate

- [~] Unit-test identity-angle equivalence on zero-mean inputs, direct reference agreement, target- and predictor-DC
      behavior, noise-variance reduction, rank/mask invalidation, no double derotation, multiple annuli/mode planes,
      1-versus-multiple-thread determinism, and injected worker/solver failures.
      Direct-SVD agreement, constant preservation, rank/mask invalidation, no-double-derotation evidence, annulus/mode
      behavior, OpenMP determinism, and worker/solver failure propagation are covered.  The explicit detector-versus-
      rotated identity/zero-mean equivalence and dedicated noise-variance-reduction cases remain.
- [ ] Add injected companions across separation, position angle, contrast, rotation coverage, and mode count.  Report
      photometric throughput, astrometric bias, residual variance, and the validity footprint.  Regression driven by
      a sky-stationary companion or unacceptable subtraction through the adopted predictor baseline is a test failure.
- [ ] Repeat the injected-companion comparison with no baseline subtraction, radial-profile subtraction,
      median-image subtraction, and both.  Measure the retained static PSF pedestal separately from companion
      throughput so the preprocessing tradeoff is visible rather than absorbed into one residual-noise metric.
- [x] Compare current detector P4 and the direct rotated implementation on the AF Lep dataset.  The comparison above
      found healthy rank and support but rejected the raw-predictor DC convention: rotated planes were nearly
      mode-invariant, matched independently reconstructed no-P4 stacks, and had worse annular RMS and planet SNR than
      detector P4.  This completes the comparison as a failed scientific gate, not as acceptance.  Materialized
      two-stage and rejected shift-then-rotate realizations remain controlled-small-problem validation oracles rather
      than production candidates.
- [ ] Benchmark rotation preparation, direct feature sampling, centering, Gram construction, eigensolve, projection,
      and total wall time with controlled OpenMP/BLAS settings.  A small benchmark may report the materialized
      oracle's one-time cube rotation plus second-stage OR interpolation, but production acceptance is based on the
      direct path's scientific tests, full AF geometry runtime, and bounded memory use.
- [x] Keep the current smaller-Gram selection and thread-local workspaces.  Do not combine this work with the planned
      CUDA backend; establish the CPU scientific reference first so a future GPU implementation has an oracle.
- [~] Reach 100% executable-line and logical-function coverage for every new or changed production path, run Release,
      coverage, documentation, ASan/UBSan, and OpenMP determinism tests, and audit every newly called mxlib
      specialization against mxlib's 100% coverage gate before acceptance.
      The changed hciReduce production set is at 100% focused line/function coverage; Release, coverage, and
      documentation gates pass.  Acceptance remains partial because sanitizer and determinism checks have not been
      rerun for the superseding uncentered-application rule, the exact float `eigenCube::median()` mxlib coverage gap
      is recorded but unresolved, and the injected-source and remaining scientific-validation gates have not run.

### Phase 5: documentation and release criteria

- [x] Update the P4 User's Guide and configuration dictionary with the exact equations, predictor-frame choice,
      default, centered-fit and uncentered-application rules, validity behavior, performance implications, CLI form,
      config-file form, type/range/default, and examples.
- [x] Document the separate roles and exact execution order of detector-frame radial-profile/median-PSF preprocessing
      and rotated-frame temporal centering, including field-rotation and companion-throughput limitations.
- [x] Add a migration note stating that `detector` remains the default and reproduces the current algorithm, while
      `rotated` changes the coefficient field, structural DOF, validity support, and residual statistics.
- [ ] Treat the rotated variant as scientifically ready only after a new baseline/output contract is selected and
      passes synthetic throughput tests, AF Lep comparison, complete provenance, stable progress reporting, and full
      test/coverage/documentation gates.  The current raw-predictor DC convention failed the AF Lep gate.
