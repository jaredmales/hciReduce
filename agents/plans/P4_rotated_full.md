# Full rotated-frame P4 with the original coefficient field

## Status and scope

This document records a possible replacement for the current per-sky-pixel rotated P4 approximation described in
[`P4_rotated.md`](P4_rotated.md).  It returns to the original detector-frame coefficient field
\(c_{k,n}\) and minimizes a temporally centered residual after rotation.  The purpose of this document is to preserve
the mathematical discussion, compare feasible solvers, and define a staged investigation.  It does **not** approve
or begin an implementation.

The current rotated implementation ties one predictor vector \(c_{k,m}\) to each sky pixel \(m\).  That produces a
small independent problem at every \(m\), but it is a different model from the original P4 coefficient field.  The
original equations remain linear after rotation and centering; their design matrix is simply much larger and is no
longer block diagonal in detector pixel.

Two versions of that larger problem must be kept distinct:

1. A **local expanded problem** solves each sky pixel independently, allowing a separate local copy
   \(c_{k,n}^{(m)}\) of every detector coefficient touched by that sky pixel's rotation kernels.
2. A **globally shared problem** uses one literal detector-frame field \(c_{k,n}\) for every sky pixel.  Overlapping
   sky-pixel problems then share coefficients and must be solved together.

The local problem is a useful stepping stone and numerical oracle, but only the globally shared problem preserves the
original meaning of \(c_{k,n}\).

## Scientific objective

P4 is intended to predict the stellar PSF/noise background, not a companion.  Let \(x_t[n]\) be detector-frame pixel
\(n\) in frame \(t\), and let the original P4 model be

\[
   \widehat{x}_t[n] = \sum_{k\in\mathcal K_n} c_{k,n}\,g_{t,n,k},
   \qquad
   g_{t,n,k} = x_t[n+\Delta n_k],
\]

where \(\mathcal K_n\) is the accepted optimization-region predictor set for detector pixel \(n\).  The notation
\(g_{t,n,k}\) may include the checked two-dimensional sampling operation used by production code; flattened index
arithmetic is not implied.

For detector-to-sky rotation operator \(R_t\), define

\[
   (R_t x_t)[m] = \sum_n r_{t,m}[n]x_t[n].
\]

The desired objective measures temporal variation of the residual at a fixed sky pixel after rotation.  With

\[
   H = I-\frac{1}{T}\mathbf 1\mathbf 1^T,
\]

the objective is

\[
   J_m(c)=\frac{1}{2}\left\|H\left(y_m-Q_mc\right)\right\|_2^2,
   \qquad
   y_m[t]=(R_tx_t)[m].
\]

Centering in the sky frame removes a truly stationary companion from the variance being minimized.  It also removes
every time-constant component of the PSF prediction from the data term, which creates the identifiability problem
described below.

## Rotation-kernel support bounds the local problem

For a fixed sky pixel \(m\), let

\[
   F_{t,m}=\text{the complete nominal 4-by-4 footprint of }r_{t,m},
   \qquad
   S_m=\bigcup_t F_{t,m}.
\]

The nominal footprint is used even when a cardinal interpolation phase makes one of its coefficients algebraically
zero.  This matches the complete-support validity contract and the support counts below; it is deliberately more
conservative than defining \(F_{t,m}\) as only the nonzero algebraic support of \(r_{t,m}\).

Only detector pixels in \(S_m\) can affect the rotated value at \(m\).  The local column domain is therefore

\[
   \mathcal C_m=\{(n,k):n\in S_m,\ k\in\mathcal K_n\},
\]

not the complete image-wide coefficient field.  Its design matrix is

\[
   Q_m[t,(n,k)] = r_{t,m}[n]g_{t,n,k},
   \qquad (n,k)\in\mathcal C_m.
\]

The predictor *sample* domain can extend outside \(S_m\) through the offsets \(n+\Delta n_k\), but the unknown
coefficient domain cannot: a coefficient anchored at \(n\notin S_m\) is multiplied by zero by every rotation row
for this \(m\).

For a cubic rotation kernel, one row touches at most 16 detector pixels and hence roughly
\(16|\mathcal K_n|\) coefficient columns.  Across time, the union \(S_m\) follows the inverse-rotated path of \(m\)
through the detector and is typically a short arc or annular segment rather than the whole image.  Mask, edge, and
complete-interpolation-support rules may reduce the admissible domain further, but they must do so explicitly and
consistently for the entire fit.

A preliminary AF Lep count, before correcting the angle anomaly discussed below, found approximately
\(|S_m|=74\)--195 detector pixels.  The current per-sky-pixel annulus geometry provisionally gives 611--1,713
predictors per detector pixel, so a local expanded fit has roughly 45,000--334,000 columns but only \(T=621\) rows.
The full formulation still needs to define \(\mathcal K_n\) when a detector pixel in the rotation footprint crosses
a configured annulus boundary; the quoted predictor and column counts must be repeated after that decision.  A dense
double-precision \(Q_m\) would require about 0.22--1.66 GB for one sky pixel.  It is sparse by row, and its centered
rank is at most \(T-1=620\), but it is extremely underdetermined.  These counts must be repeated after the derotation
angles have been corrected.

## Gradient, Hessian, and normal equations

Define

\[
   \widetilde y_m=Hy_m,
   \qquad
   \widetilde Q_m=HQ_m.
\]

Then

\[
   J_m(c)=\frac12\|\widetilde y_m-\widetilde Q_mc\|_2^2,
\]

with

\[
   \nabla J_m(c)
     =Q_m^TH(Q_mc-y_m)
     =\widetilde Q_m^T(\widetilde Q_mc-\widetilde y_m),
\]

\[
   \nabla^2J_m
     =Q_m^THQ_m
     =\widetilde Q_m^T\widetilde Q_m.
\]

The normal equations are therefore

\[
   Q_m^THQ_m c=Q_m^THy_m.
\]

They are linear, but the primal matrix has tens to hundreds of thousands of rows and columns and rank no greater
than 620.  It must never be explicitly formed for production.

For the local expanded approximation, each \(m\) obtains its own solution \(c^{(m)}\).  If two rotation footprints
contain the same detector coefficient \((n,k)\), their independently fitted values need not agree.  This relaxation
can obtain a lower residual than the globally constrained model, but it does not produce one detector-frame PSF
model that can be applied consistently to the whole image.

## Globally shared coefficient formulation

To retain one original coefficient field, define the detector-frame prediction

\[
   P_t(c)[n]=\sum_{k\in\mathcal K_n}c_{k,n}g_{t,n,k}.
\]

Stacking all retained sky pixels and frames, let \(\mathcal H\) apply temporal centering independently at each sky
pixel,

\[
   (\mathcal Hz)_t[m]
      =z_t[m]-\frac{1}{T}\sum_\tau z_\tau[m].
\]

Define the matrix-free operator

\[
   (\mathcal A c)_t[m]
     =\mathcal H\{(R_tP_t(c))[m]\},
   \qquad
   b_t[m]=\mathcal H\{(R_tx_t)[m]\},
\]

where \(\mathcal H\) acts on the complete time stack in each expression.  The global objective is

\[
   J(c)=\frac12\|\mathcal Ac-b\|_2^2.
\]

Its gradient and Hessian action are

\[
   \nabla J(c)=\mathcal A^*(\mathcal Ac-b),
   \qquad
   \nabla^2J\,v=\mathcal A^*\mathcal Av.
\]

For a centered residual stack \(w_t[m]\), the adjoint accumulation is

\[
   (\mathcal A^*w)_{k,n}
      =\sum_t g_{t,n,k}\,[R_t^T(\mathcal Hw)_t][n].
\]

This gives an exact forward/adjoint pair without constructing the global design matrix:

1. form a detector prediction from the coefficient field and predictor samples;
2. rotate it into the sky frame;
3. center each sky-pixel time series;
4. for the adjoint, center the residual, apply the transpose of the exact rotation/interpolation operator, and
   accumulate predictor-weighted detector coefficients.

The adjoint must use the transpose of the implemented interpolation operator, not an ordinary inverse rotation.
A dot-product test, \(\langle\mathcal Ac,w\rangle=\langle c,\mathcal A^*w\rangle\), is a mandatory numerical gate.

## Temporal-dual formulation for one sky pixel

The local problem has many more columns than temporal samples.  Its natural dense object is therefore the temporal
Gram matrix

\[
   K_{0,m}=Q_mQ_m^T,
   \qquad
   K_m=HK_{0,m}H
      =\widetilde Q_m\widetilde Q_m^T,
\]

both of which are only \(T\times T\).  The raw Gram elements can be accumulated without materializing \(Q_m\):

\[
   K_{0,m}[t,s]
   =\sum_{n\in F_{t,m}\cap F_{s,m}}
      r_{t,m}[n]r_{s,m}[n]
      \sum_{k\in\mathcal K_n}g_{t,n,k}g_{s,n,k}.
\]

Frame pairs with disjoint rotation footprints contribute exactly zero to \(K_{0,m}\).  Temporal centering generally
makes \(K_m\) dense through row/column mean corrections, even when \(K_{0,m}\) is sparse.

Let the nonzero eigensystem be

\[
   K_mU=U\Lambda.
\]

For \(q\) retained eigenmodes, the truncated minimum-norm coefficient solution is

\[
   c_q=Q_m^THU_q\Lambda_q^{-1}U_q^THy_m.
\]

The centered fitted time series is especially simple:

\[
   HQ_mc_q=U_qU_q^THy_m.
\]

The raw, uncentered prediction can also be computed without constructing \(c_q\):

\[
   Q_mc_q
      =K_{0,m}HU_q\Lambda_q^{-1}U_q^THy_m
      =K_{0,m}U_q\Lambda_q^{-1}U_q^THy_m,
\]

where the second equality uses the fact that every retained nonzero eigenvector lies in the centered subspace.
This identity is important for testing proposed uncentered-output conventions against explicit tiny-problem
coefficients.

The temporal-dual formulation applies directly only to independent local \(m\) problems.  Once the same
\(c_{k,n}\) is constrained to be shared by overlapping sky pixels, their rows must be stacked and a single global
operator solved; independent \(T\times T\) eigensystems no longer enforce the shared field.

## Solver options

### Independent local expanded problems

| Priority | Method | Strengths | Limitations / role |
|---:|---|---|---|
| 1 | Temporal-dual eigensolve / TSVD | Uses a 621-by-621 matrix, naturally produces the existing sequence of mode-count products, and never needs the enormous coefficient vector | Requires a Gram construction and eigensolve for every sky pixel; gives local coefficient copies, not one shared field |
| 2 | Matrix-free LSMR or LSQR | Uses exact forward and adjoint operations, avoids squaring the condition number, and provides an independent check of the dual result | A separate solve or stopping rule is needed for each regularization choice; does not map as directly to retained eigenmode counts |
| 3 | Dual ridge solve | Solves \((K_m+\lambda I)\alpha=Hy_m\), with \(c=Q_m^TH\alpha\); Cholesky, CG, or the same eigensystem can serve many \(\lambda\) values | Changes the user-facing regularization from mode count to \(\lambda\), and still produces only a local field |
| 4 | CGLS / PCG on primal normal equations | Matrix-free and straightforward on a GPU | Squares the conditioning and needs a strong preconditioner; less useful as an independent numerical oracle |
| 5 | L-BFGS or first-order gradient methods | Need only the documented gradient | Discard the exact quadratic least-squares structure and are unlikely to outperform LSQR/LSMR |
| Oracle only | Explicit QR or SVD | Most direct solution for tiny synthetic cases | Prohibitive memory for AF Lep dimensions; never a production path |

The temporal-dual TSVD is the preferred first investigation for the local problem.  LSMR/LSQR should be retained as
an algorithmically independent reference.  Explicit primal normal equations and an unrestricted inverse are not
viable.

### One globally shared detector coefficient field

| Priority | Method | Assessment |
|---:|---|---|
| 1 | Matrix-free LSMR / LSQR | Best direct match to the global least-squares problem.  It uses only the exact global forward/adjoint pair and does not form either Gram matrix. |
| 2 | Damped LSMR / augmented LSQR | Preferred global ridge implementation because it keeps the better conditioning of the rectangular problem and admits an explicit prior operator. |
| 3 | Global ridge PCG / CGLS | Potentially GPU-friendly if an effective block, radial, or temporal preconditioner can be constructed; otherwise convergence may be poor. |
| 4 | ADMM / consensus decomposition | Allows per-sky-pixel or per-tile subproblems while enforcing agreement of overlapping \(c_{k,n}\).  Attractive for distributed execution, but introduces penalty tuning and a second convergence loop. |
| Research | Spatially regularized proximal methods | Appropriate only if sparsity, smoothness, or another physical coefficient prior is adopted deliberately. |

A global dense QR/SVD, a global primal Hessian, and a global direct inverse are outside the feasible design space.

## Identifiability and the static temporal mean

Temporal centering has two unavoidable consequences:

1. \(\operatorname{rank}(HQ_m)\le T-1\).  In the AF Lep local problem there are at least tens of thousands of
   coefficient null directions even before numerical rank loss.
2. The objective constrains only \(HQ_mc\).  If \(d\) satisfies \(HQ_md=0\), then \(c\) and \(c+d\) have the same
   objective while their raw predictions \(Q_mc\) can have different temporal means.

Consequently, the centered prediction is identifiable up to the chosen retained temporal subspace, but the temporal
mean of an uncentered PSF prediction is not determined by the data term.  TSVD, minimum norm, ridge, or an iterative
stopping rule selects one member of this solution family; that choice is a prior or convention, not information
recovered from the centered residual.

This explains the behavior of the current mean-preserving rotated implementation.  Subtracting
\(HQ_mc\) from the uncentered target leaves the target mean exactly unchanged, so ordinary temporal mean combination
cannot remove a static radial PSF pedestal.  Conversely, subtracting the raw \(Q_mc\) requires accepting and
validating the solver's implicit choice for the unconstrained mean.  Fitting a free sky-frame intercept would also
be unable to distinguish a stationary companion from a stationary PSF component.

Possible ways to make the raw mean scientifically meaningful include:

- detector-frame radial-profile or temporal median-PSF subtraction before the centered rotated regression;
- an explicit detector-frame baseline model estimated outside the centered objective;
- ridge or a physically justified spatial prior on \(c_{k,n}\);
- held-out-frame prediction that penalizes unstable null-space choices; and
- a joint probabilistic model with separately parameterized detector-fixed PSF and sky-fixed source terms.

None of these follows automatically from choosing a numerical solver.  The baseline/output contract must be frozen
before a full implementation is selected.

## Regularization choices to compare

- **Truncated SVD/eigendecomposition:** retains the current P4 mode-count interpretation and gives all requested
  truncations from one local eigensystem.  Small eigenvalues must be excluded using an explicit numerical threshold.
- **Isotropic ridge:** minimizes
  \(\|H(y-Qc)\|^2+\lambda\|c\|^2\).  It supplies a unique minimum-norm convention and has a simple local dual, but
  \(\lambda\) needs a new selection and provenance contract.
- **Spatial smoothness:** penalizes differences between neighboring \(c_{k,n}\).  This directly couples pixels and
  is most natural in the global formulation, but it changes the scientific model and generally requires an
  iterative solver.
- **Structured coefficient tying:** shares or smoothly varies coefficients within annuli, radial/tangential bins, or
  other geometry.  This interpolates between the current \(c_{k,m}\) approximation and a completely free global
  field, but any tying geometry is a scientific prior.
- **Sparsity or non-negativity:** possible through proximal or constrained methods, but neither is currently
  justified by the prototype.  P4 coefficients need not naturally be positive.
- **Early-stopped Krylov iteration:** acts as implicit spectral regularization.  It is useful for comparison, but an
  iteration tolerance is less transparent than an explicit rank or ridge parameter.

Regularization should be selected with held-out residual prediction and injected companions, not in-sample residual
variance alone.  Validation must measure photometric throughput, astrometric bias, residual noise, raw-mean
stability, and sensitivity to masks and angle coverage.

## Computational considerations

- A local \(K_{0,m}\) contains only about 386,000 values (about 3.1 MB in double precision), but naive construction
  is \(O(T^2|F_t\cap F_s||\mathcal K|)\), and a dense \(O(T^3)\) eigensolve would be repeated for every output pixel.
- Disjoint footprint pairs should be pruned before evaluating predictor inner products.  Geometry, predictor
  products, and repeated annulus structure should be cached or tiled only after a correct oracle exists.
- Centering should be applied as row/column mean operations; there is no need to construct \(H\).
- The local matrix-free forward cost is proportional to the active rotation-kernel samples and predictor count per
  iteration.  The adjoint reuses exactly the same records and must be deterministic under threading.
- The global forward/adjoint is well suited to GPU kernels: construct detector predictions, rotate/adjoint-rotate,
  center, and reduce coefficient gradients.  Transfer cost, storage of predictor samples, and preconditioning may
  dominate.
- Batched temporal eigensolves, GPU acceleration, and distributed consensus are performance phases, not substitutes
  for a CPU double-precision reference.

## AF Lep derotation-angle anomaly

The current AF Lep sequence is not yet a trustworthy golden dataset for this work.  The raw `POSANG` values in
`pp_000350` through `pp_000353` are \(-216.8021^\circ\), \(-123.2917^\circ\), \(-31.14813^\circ\), and
\(54.64996^\circ\).  After the configured \(+1.6^\circ\) offset and normalization, these become derotation angles
near \(144.7979^\circ\), \(238.3083^\circ\), \(330.45187^\circ\), and \(56.24996^\circ\), while `pp_000349`
and `pp_000354` give \(91.3809^\circ\) and \(90.41377^\circ\).  This is consistent with ordinary arithmetic
averaging across an angular wrap during coaddition, but that cause must be confirmed from the pre-coadd headers.

These angles change the rotation-kernel paths, inflate \(S_m\), alter overlap sparsity, and can corrupt both solver
conditioning and a visual/scientific comparison.  Before AF Lep establishes performance or correctness:

1. recover the contributing pre-coadd angles;
2. recompute each coadd angle with circular statistics or the appropriate unwrapped sequence;
3. verify the detector-to-sky sign and convention against neighboring frames and a known source;
4. regenerate or patch the test input with explicit provenance; and
5. repeat the \(|S_m|\), column-count, runtime, and memory estimates above.

## Staged investigation plan

Checklist status: `[ ]` not started.  All work below remains unapproved planning.

### Phase 0: freeze the mathematical and scientific contract

- [ ] Decide whether the intended endpoint is the independent local expanded problem or one globally shared
      detector coefficient field.  Do not call local copies \(c_{k,n}^{(m)}\) the original global coefficients.
- [ ] Define the detector-frame predictor geometry \(\mathcal K_n\), including which annulus owns a coefficient
      pixel when the nominal rotation footprint crosses a configured search-annulus boundary.
- [ ] Decide whether the written science residual subtracts the centered prediction, the raw minimum-norm/ridge
      prediction, or a separately modeled detector-frame baseline plus centered prediction.
- [ ] Specify the coefficient prior or regularization that makes a raw prediction meaningful, including its
      parameter-selection and FITS-provenance contract.
- [ ] Define masks, interpolation support, image edges, and whether every fit requires a common all-frame valid set.
- [ ] Correct and validate the AF Lep angle sequence before using its support sizes or products as acceptance data.

### Phase 1: tiny explicit mathematical oracle

- [ ] On a very small image and frame count, enumerate \(\mathcal C_m\), construct \(Q_m\) explicitly, and compare
      the analytic gradient with central finite differences.
- [ ] Compare dense SVD/QR, temporal-dual TSVD, ridge, and matrix-free LSMR/LSQR solutions on full-rank,
      rank-deficient, and deliberately ill-conditioned examples.
- [ ] Verify both \(HQ_mc_q=U_qU_q^THy_m\) and the raw-prediction formula against explicit coefficients for every
      retained rank.
- [ ] Add forward/adjoint dot-product tests for the local operator and the globally stacked rotation operator.
- [ ] Demonstrate explicitly that different coefficient solutions can have identical centered predictions but
      different raw temporal means.

### Phase 2: local temporal-dual research implementation

- [ ] Build \(K_{0,m}\) from rotation-footprint intersections without materializing \(Q_m\), initially in a
      standalone/reference component rather than the production reduction path.
- [ ] Compare every result with the tiny explicit oracle and with matrix-free LSMR/LSQR.
- [ ] Produce all requested TSVD ranks from one eigensystem and verify exact structural rank \(\le T-1\).
- [ ] Measure Gram construction, centering, eigensolve, raw/centered projection, peak memory, and scaling with
      \(|S_m|\), \(|\mathcal K_n|\), and \(T\).
- [ ] Determine whether per-pixel temporal eigensolves are viable before adding caching, batching, OpenMP, or GPU
      optimization.

### Phase 3: globally shared matrix-free prototype

- [ ] Implement the exact global forward and adjoint only after their small dense forms and dot-product identity are
      proven.
- [ ] Compare LSMR/LSQR, damped LSMR, and a preconditioned ridge solve on synthetic coefficient fields with known
      truth.
- [ ] Quantify the disagreement among independently fitted overlapping local coefficients and the residual penalty
      imposed by global sharing.
- [ ] Explore ADMM/consensus decomposition only if the monolithic matrix-free solve is not practical; record primal
      and dual convergence, communication volume, and penalty sensitivity.
- [ ] Keep CPU double precision as the scientific oracle before adding CUDA or distributed execution.

### Phase 4: scientific selection

- [ ] Use synthetic detector-fixed PSFs plus injected sky-fixed companions to compare centered subtraction, raw
      minimum-norm subtraction, ridge, and explicit baseline treatments.
- [ ] Sweep separation, position angle, contrast, field rotation, mode/ridge strength, masks, and kernel support.
- [ ] Report companion throughput, astrometric bias, residual variance, temporal-mean bias, finite support, and
      failure/conditioning diagnostics separately.
- [ ] Use blocked and interleaved held-out frames so regularization is not selected from in-sample residual variance.
- [ ] Repeat the accepted synthetic tests on the corrected AF Lep sequence and compare against the prototype and
      current detector/per-sky-pixel implementations.

### Phase 5: production decision and optimization

- [ ] Select local versus global coefficients, solver, regularization, baseline policy, and output contract from the
      scientific results rather than runtime alone.
- [ ] Define configuration names, defaults, progress reporting, diagnostics, failure policy, and complete FITS
      provenance only after that selection.
- [ ] Optimize footprint intersection, predictor-product reuse, tiling, batched eigensolves, and GPU kernels against
      the frozen CPU oracle.
- [ ] Require deterministic multi-thread tests, ASan/UBSan, full line/function coverage, mxlib coverage auditing,
      documentation, and real-data benchmarks before replacing or supplementing the current rotated path.

## Exit criteria for the planning stage

No production implementation should begin until Phase 0 resolves:

- whether coefficients are local copies or one shared detector field;
- how the unconstrained temporal mean of the prediction is modeled;
- which regularization defines a unique/stable solution;
- the valid-support and mask contract; and
- a corrected, provenance-preserving AF Lep angle sequence.

The first approved code should be the tiny explicit oracle, not an optimized solver.
