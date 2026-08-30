/** \file P4PCA.hpp
 * \brief Declares the pure numerical predictor used by Pixel Prediction Post-Processing.
 * \author Jared R. Males
 */

#ifndef P4PCA_hpp
#define P4PCA_hpp

#include <cstdint>
#include <span>
#include <vector>

#include <Eigen/Dense>

#include <mx/math/eigenLapack.hpp>
#include <mx/math/svdDowndate.hpp>

namespace mx
{
namespace improc
{

/// Rank support of one requested P4 residual plane.
/** All numerical failures other than insufficient rank throw, so these are the only two non-throwing plane states.
 *
 * \ingroup programming_library
 */
enum class P4PCAModeStatus : std::uint8_t
{
    rankSupported,   ///< The requested mode count is within the numerical rank.
    rankInsufficient ///< The requested mode count exceeds the numerical rank.
};

/// Numerical implementation used for target-frame exclusion.
/** \ingroup programming_library */
enum class P4ExclusionSolver : std::uint8_t
{
    explicitRefit,      ///< Refit the physically retained predictor rows independently for every target.
    factorDowndateExact ///< Delete target rows from one complete safely represented base factorization.
};

/// Reason an exact factor-downdated search pixel was recomputed with the explicit held-out oracle.
/** \ingroup programming_library */
enum class P4PCAFallbackReason : std::uint8_t
{
    none,             ///< The factor-downdated calculation completed without an explicit fallback.
    rankBoundary,     ///< A downdated eigenvalue was numerically indistinguishable from the rank threshold.
    factorValidation, ///< The base temporal factor failed its numerical orthogonality validation.
    deletionSolver    ///< The selected structured row-deletion solver reported a recoverable numerical failure.
};

/// Storage model used by a compact collection of target-specific deleted rows.
/** \ingroup programming_library */
enum class P4ExclusionStorage : std::uint8_t
{
    contiguousSpans, ///< Each target owns one half-open contiguous deleted-row span.
    explicitIndices  ///< Target rows index a flattened collection of arbitrary sorted deleted-row indices.
};

/// One half-open target-frame exclusion span.
/** \ingroup programming_library */
struct P4ExclusionSpan
{
    Eigen::Index begin{ 0 }; ///< First deleted temporal row.

    Eigen::Index end{ 0 };   ///< One past the last deleted temporal row.
};

/// Compact target-specific deleted-row patterns shared by explicit and factor-downdated P4 fits.
/** Contiguous image-number windows require only two indices per target. Angle- and pixel-based exclusions use one
 * flattened index array plus target offsets. Every nonempty pattern must delete its own target row.
 *
 * \ingroup programming_library
 */
class P4TargetExclusions
{
  public:
    /// Construct an empty collection representing in-sample fitting.
    P4TargetExclusions() = default;

    /// Construct validated half-open contiguous exclusion spans.
    static P4TargetExclusions
    fromSpans( Eigen::Index sampleCount, /**< [in] temporal row count */
               const std::vector<P4ExclusionSpan> &spans /**< [in] one target-containing span per temporal row */ );

    /// Construct validated arbitrary exclusion sets and flatten their storage.
    static P4TargetExclusions
    fromExplicit( Eigen::Index sampleCount, /**< [in] temporal row count */
                  const std::vector<std::vector<Eigen::Index>> &indices
                  /**< [in] one sorted, unique, target-containing index set per temporal row */ );

    /// Report whether this collection represents in-sample fitting.
    bool empty() const noexcept;

    /// Return the temporal row count represented by the collection.
    Eigen::Index sampleCount() const noexcept;

    /// Return the target-specific pattern count.
    Eigen::Index targetCount() const noexcept;

    /// Return the selected compact storage model.
    P4ExclusionStorage storage() const noexcept;

    /// Return the deleted-row count for one target.
    Eigen::Index deletedCount( Eigen::Index target /**< [in] temporal target row */ ) const;

    /// Return the largest target-specific deletion count.
    Eigen::Index maximumDeleted() const noexcept;

    /// Return a deleted-row view, materializing a contiguous span in caller-owned scratch when necessary.
    std::span<const Eigen::Index>
    deletedRows( std::vector<Eigen::Index> &spanScratch, /**< [in,out] reusable storage for a contiguous span */
                 Eigen::Index target /**< [in] temporal target row */ ) const;

    /// Return allocated vector-element storage in bytes, excluding object padding and allocator bookkeeping.
    std::size_t storageBytes() const noexcept;

  private:
    Eigen::Index m_sampleCount{ 0 }; ///< Temporal row count, or zero for the empty in-sample representation.

    P4ExclusionStorage m_storage{ P4ExclusionStorage::contiguousSpans }; ///< Selected compact representation.

    std::vector<P4ExclusionSpan> m_spans; ///< Target-indexed spans for contiguous storage.

    std::vector<Eigen::Index> m_offsets;  ///< Target offsets into `m_indices` for explicit storage.

    std::vector<Eigen::Index> m_indices;  ///< Flattened arbitrary deleted-row indices.

    Eigen::Index m_maximumDeleted{ 0 };   ///< Largest realized target-specific deletion count.
};

/// Results of one local P4 principal-component regression.
/** Residual columns and mode-status entries retain the order of the requested mode counts.
 *
 * \ingroup programming_library
 */
struct P4PCAResult
{
    /// Residual time series, with one row per sample and one column per requested mode count.
    Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> residuals;

    /// Rank-support state of each residual column.
    std::vector<P4PCAModeStatus> modeStatus;

    /// Number of eigenvalues above the configured relative rank threshold.
    int numericalRank{ 0 };

    /// Rank of the immutable base factorization used by a factor-downdated fit, or zero for a direct fit.
    int baseRank{ 0 };

    /// Whether the reported numerical rank is capped by a selected-spectrum solve.
    bool numericalRankCapped{ false };

    /// Number of roundoff-scale negative deletion-core eigenvalues clamped across target rows.
    std::size_t downdateClampCount{ 0 };

    /// Number of target rows recomputed by the explicit oracle after a recoverable downdate condition.
    std::size_t explicitFallbackCount{ 0 };

    /// Reason the complete search pixel was recomputed by the explicit oracle.
    P4PCAFallbackReason explicitFallbackReason{ P4PCAFallbackReason::none };

    /// Maximum absolute temporal-factor Gram defect that triggered factor-validation fallback, or zero.
    double factorOrthogonalityDefect{ 0 };

    /// Accepted temporal-factor Gram-defect tolerance when factor-validation fallback occurred, or zero.
    double factorOrthogonalityTolerance{ 0 };

    /// Optional per-sample, per-mode rank validity used by held-out fits; empty means modeStatus applies to every row.
    Eigen::Array<std::uint8_t, Eigen::Dynamic, Eigen::Dynamic> sampleValidity;

    /// Report whether one residual sample has sufficient numerical rank.
    bool sampleSupported( Eigen::Index sample, /**< [in] residual row */
                          std::size_t mode /**< [in] requested-mode index */ ) const;
};

/// Per-call aggregate-worker timing values for the numerical P4 PCA core.
/** \ingroup programming_library */
struct P4PCATiming
{
    double gramWorkerSeconds{ 0 };             ///< Time to form the selected normal equations.

    double eigensolveWorkerSeconds{ 0 };       ///< Time in the eigensolver and rank selection.

    double baseFactorWorkerSeconds{ 0 };       ///< Time to construct the complete reusable singular system.

    double deletionWorkerSeconds{ 0 };         ///< Time in target-specific row-deletion solves.

    double explicitFallbackWorkerSeconds{ 0 }; ///< Wall time in a complete explicit held-out fallback.

    double projectionWorkerSeconds{ 0 };       ///< Time to project modes and construct residuals.
};

/// Reusable worker-private storage for exact P4 factor deletion.
/** One instance may be reused for successive search pixels but must not be shared by concurrent workers.
 * \ingroup programming_library
 */
struct P4PCADowndateWorkspace
{
  private:
    friend struct P4PCA;

    Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> m_baseGram; ///< Base temporal or predictor Gram matrix.

    Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> m_baseEigenvectors; ///< Base Gram eigenvectors.

    Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> m_baseEigenvalues;  ///< Ascending base Gram eigenvalues.

    Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> m_temporalFactor;   ///< Descending temporal SVD factor `L`.

    Eigen::Array<double, Eigen::Dynamic, 1> m_singularValues;       ///< Descending complete safe base singular values.

    Eigen::Array<double, Eigen::Dynamic, 1> m_baseTargetProjection; ///< Full-data `L^T y` projection.

    Eigen::Array<double, Eigen::Dynamic, 1> m_retainedTargetProjection; ///< Target-excluded `L_S^T y_S`.

    Eigen::Array<double, Eigen::Dynamic, 1> m_scaledTargetProjection;   ///< `Sigma L_S^T y_S` scratch.

    Eigen::Array<double, Eigen::Dynamic, 1> m_scaledHeldOutRow;         ///< `Sigma l_t^T` scratch.

    std::vector<Eigen::Index> m_deletedRows;                            ///< Materialized contiguous deletion indices.

    mx::math::svdDeletionResult<double> m_deletionResult;       ///< Published spectrum and preserved-side rotation.

    mx::math::svdDeletionWorkspace<double> m_deletionWorkspace; ///< Reusable mxlib deletion scratch.
};

/// Pure all-double principal-component regression for one P4 search pixel.
/** This component forms the smaller of the temporal and predictor Gram matrices. calculate() uses the inputs directly,
 * while calculateCentered() centers the fit and applies the resulting coefficients to the uncentered predictors. A
 * caller-owned LAPACK workspace makes repeated calls efficient and permits one independent workspace per worker
 * thread.
 *
 * \ingroup programming_library
 */
struct P4PCA
{
    /// Dynamic all-double matrix used for predictors and residuals.
    using matrixT = Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>;

    /// Dynamic all-double column vector used for one target time series.
    using vectorT = Eigen::Array<double, Eigen::Dynamic, 1>;

    /// Reusable all-double LAPACK workspace. It must not be copied or shared by concurrent calls.
    using workspaceT = mx::math::syevrMem<double>;

    /// Calculate truncated-PCA residual time series for ordered retained-mode counts.
    /** On success, \p output is completely replaced. On exception it remains destructible, but its values are
     * unspecified. Requested counts above numerical rank produce complete quiet-NaN columns and do not invalidate
     * lower-count columns. When requested, \p coefficients receives one predictor-space coefficient vector per mode
     * count; rank-insufficient columns are quiet NaNs. The output objects must not own storage aliased by either input
     * array.
     */
    static void calculate( P4PCAResult &output,           /**< [out] residuals, mode states, and numerical rank */
                           const matrixT &predictors,     /**< [in] finite T-by-K local predictor matrix */
                           const vectorT &target,         /**< [in] finite T-sample target time series */
                           const std::vector<int> &modes, /**< [in] positive, strictly increasing retained counts */
                           double rankTolerance,  /**< [in] finite nonnegative threshold relative to lambdaMax */
                           workspaceT &workspace, /**< [in,out] caller-owned, non-shared LAPACK workspace */
                           P4PCATiming *timing = nullptr, /**< [out] optional per-call worker timing */
                           matrixT *coefficients = nullptr /**< [out] optional K-by-mode predictor coefficients */ );

    /// Calculate exact target-held-out residuals from explicit training-row sets.
    /** Each output row is predicted by a separate PCA fit containing only the complement of the corresponding
     * \p exclusions entry. Requested modes beyond an individual fit's structural or numerical rank remain NaN and are
     * marked invalid in P4PCAResult::sampleValidity. Coefficients are intentionally not returned because they differ
     * for every target row. When T is no greater than K, the implementation forms the full temporal Gram matrix once,
     * extracts each training-set principal submatrix, and applies the held-out Gram cross-vector. The K-less-than-T
     * path retains the explicit predictor-space refit.
     */
    static void calculateHeldOut(
        P4PCAResult &output,       /**< [out] held-out residuals, aggregate mode states, and minimum solved rank */
        const matrixT &predictors, /**< [in] finite full T-by-K predictor matrix */
        const vectorT &target,     /**< [in] finite full T-sample target time series */
        const P4TargetExclusions &exclusions, /**< [in] compact deleted-row pattern for every target row */
        const std::vector<int> &modes,        /**< [in] positive, strictly increasing retained counts */
        double rankTolerance,                 /**< [in] finite nonnegative threshold relative to lambdaMax */
        workspaceT &workspace,                /**< [in,out] caller-owned, non-shared LAPACK workspace */
        P4PCATiming *timing = nullptr /**< [out] optional aggregate worker timing for all refits */ );

    /// Calculate exact target-held-out residuals with the default dense deletion backend.
    /** The base keeps every safely representable singular triplet independently of \p rankTolerance. Each target's
     * compact exclusion set is passed to mxlib's full-spectrum leading-covariance deletion backend, after which the
     * existing held-out principal-component regression is evaluated entirely in factor coordinates. This path is
     * algebraically exact for the represented base design; the explicit calculateHeldOut() method remains the oracle.
     * If the temporal factor fails only its numerical orthogonality bound, or if any post-deletion eigenvalue is
     * numerically indistinguishable from the rank threshold, every target for that search pixel is recomputed with
     * calculateHeldOut(). P4PCAResult records the fallback count, reason, and factor defect when applicable. Other
     * unsafe base states and deletion-solver failures throw instead of falling back.
     */
    static void calculateHeldOutDowndated(
        P4PCAResult &output,                       /**< [out] held-out residuals and deletion diagnostics */
        const matrixT &predictors,                 /**< [in] finite full T-by-K predictor matrix */
        const vectorT &target,                     /**< [in] finite full T-sample target time series */
        const P4TargetExclusions &exclusions,      /**< [in] compact deleted-row pattern for every target row */
        const std::vector<int> &modes,             /**< [in] positive, strictly increasing retained counts */
        double rankTolerance,                      /**< [in] finite nonnegative threshold relative to lambdaMax */
        workspaceT &eigensolverWorkspace,          /**< [in,out] reusable base eigensolver workspace */
        P4PCADowndateWorkspace &downdateWorkspace, /**< [in,out] reusable factor-deletion storage */
        P4PCATiming *timing = nullptr /**< [out] optional aggregate base, deletion, and projection timing */ );

    /// Calculate exact target-held-out residuals with an explicitly selected row-deletion backend.
    /** `leadingCovariance` accepts arbitrary validated exclusion sizes. `rankOneSecular` requires exactly one deleted
     * row for every target that retains any training rows; all-row-excluded targets remain valid rank-zero skips. The
     * backend is selected once for the entire search pixel. Recoverable structured-solver failures discard any partial
     * target results and recompute the complete search pixel with calculateHeldOut(). Other invalid or resource-failure
     * states throw.
     */
    static void calculateHeldOutDowndated(
        P4PCAResult &output,                          /**< [out] held-out residuals and deletion diagnostics */
        const matrixT &predictors,                    /**< [in] finite full T-by-K predictor matrix */
        const vectorT &target,                        /**< [in] finite full T-sample target time series */
        const P4TargetExclusions &exclusions,         /**< [in] compact deleted-row pattern for every target row */
        const std::vector<int> &modes,                /**< [in] positive, strictly increasing retained counts */
        double rankTolerance,                         /**< [in] finite nonnegative threshold relative to lambdaMax */
        mx::math::svdDeletionBackend deletionBackend, /**< [in] supported dense or one-row structured backend */
        workspaceT &eigensolverWorkspace,             /**< [in,out] reusable base eigensolver workspace */
        P4PCADowndateWorkspace &downdateWorkspace,    /**< [in,out] reusable factor-deletion storage */
        P4PCATiming *timing = nullptr /**< [out] optional aggregate base, deletion, and projection timing */ );

    /// Calculate a temporally centered fit and apply its predictor coefficients to the uncentered data.
    /** Each predictor column and the target are centered over their T samples when estimating the truncated-PCA
     * coefficient vector. That vector is then applied to the original uncentered predictor matrix and subtracted from
     * the original uncentered target. The predictor means therefore affect the residual mean even though they do not
     * affect the fitted coefficients. Centering limits the structural degrees of freedom to min(K,T-1); requested
     * counts above that limit are rejected even when roundoff would make the centering null eigenvalue positive. T
     * must be at least two. The centered objective does not establish whether the applied predictor mean is a
     * physically valid baseline; that interpretation remains the caller's responsibility. On success, \p output is
     * completely replaced. On exception it remains destructible, but its values are unspecified. The output object
     * must not own storage aliased by either input array. When requested, \p coefficients receives the centered-fit
     * coefficient vectors applied to the original uncentered predictors.
     */
    static void
    calculateCentered( P4PCAResult &output,           /**< [out] residuals, mode states, and numerical rank */
                       const matrixT &predictors,     /**< [in] finite T-by-K uncentered predictor matrix */
                       const vectorT &target,         /**< [in] finite T-sample uncentered target time series */
                       const std::vector<int> &modes, /**< [in] positive, strictly increasing retained counts */
                       double rankTolerance,          /**< [in] finite nonnegative threshold relative to lambdaMax */
                       workspaceT &workspace,         /**< [in,out] caller-owned, non-shared LAPACK workspace */
                       P4PCATiming *timing = nullptr, /**< [out] optional per-call worker timing */
                       matrixT *coefficients = nullptr /**< [out] optional K-by-mode predictor coefficients */ );

    /// Calculate a centered fit while reusing the caller's predictor matrix as centering workspace.
    /** This is numerically equivalent to calculateCentered(), including applying the fitted coefficients to the
     * uncentered predictors. On return or after a post-validation exception, \p predictors no longer contains its
     * original values. The target remains unchanged. This destructive form avoids one complete T-by-K allocation for
     * callers that refill their predictor matrix before its next use.
     */
    static void calculateCenteredInPlace(
        P4PCAResult &output,           /**< [out] residuals, mode states, and numerical rank */
        matrixT &predictors,           /**< [in,out] finite T-by-K predictors, replaced by centered values */
        const vectorT &target,         /**< [in] finite T-sample uncentered target time series */
        const std::vector<int> &modes, /**< [in] positive, strictly increasing retained counts */
        double rankTolerance,          /**< [in] finite nonnegative threshold relative to lambdaMax */
        workspaceT &workspace,         /**< [in,out] caller-owned, non-shared LAPACK workspace */
        P4PCATiming *timing = nullptr, /**< [out] optional per-call worker timing */
        matrixT *coefficients = nullptr /**< [out] optional K-by-mode predictor coefficients */ );
};

/// \cond P4PCA_test_detail
namespace detail
{

/// Signature of the Doxygen-hidden eigensolver seam used by focused failure tests.
using P4PCAEigenSolverT =
    MXLAPACK_INT ( * )( P4PCA::matrixT &, P4PCA::matrixT &, P4PCA::matrixT &, int, P4PCA::workspaceT & );

/// Install a process-wide P4PCA test eigensolver.
void p4PCASetEigenSolverForTesting(
    P4PCAEigenSolverT solver /**< [in] test solver, or nullptr to select the production solver */ );

/// Restore the production P4PCA eigensolver.
void p4PCAResetEigenSolverForTesting();

} // namespace detail
/// \endcond

} // namespace improc
} // namespace mx

#endif // P4PCA_hpp
