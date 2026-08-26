/** \file P4PCA.hpp
 * \brief Declares the pure numerical predictor used by Pixel Prediction Post-Processing.
 * \author Jared R. Males
 */

#ifndef P4PCA_hpp
#define P4PCA_hpp

#include <cstdint>
#include <vector>

#include <Eigen/Dense>

#include <mx/math/eigenLapack.hpp>

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
    double gramWorkerSeconds{ 0 };       ///< Time to form the selected normal equations.

    double eigensolveWorkerSeconds{ 0 }; ///< Time in the eigensolver and rank selection.

    double projectionWorkerSeconds{ 0 }; ///< Time to project modes and construct residuals.
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
    /** Each output row is predicted by a separate PCA fit containing only the corresponding entries of
     * \p referenceRows. Requested modes beyond an individual fit's structural or numerical rank remain NaN and are
     * marked invalid in P4PCAResult::sampleValidity. Coefficients are intentionally not returned because they differ
     * for every target row.
     */
    static void calculateHeldOut(
        P4PCAResult &output,       /**< [out] held-out residuals, aggregate mode states, and minimum solved rank */
        const matrixT &predictors, /**< [in] finite full T-by-K predictor matrix */
        const vectorT &target,     /**< [in] finite full T-sample target time series */
        const std::vector<std::vector<std::size_t>> &referenceRows,
        /**< [in] retained training-row indices for each target row */
        const std::vector<int> &modes, /**< [in] positive, strictly increasing retained counts */
        double rankTolerance,          /**< [in] finite nonnegative threshold relative to lambdaMax */
        workspaceT &workspace,         /**< [in,out] caller-owned, non-shared LAPACK workspace */
        P4PCATiming *timing = nullptr /**< [out] optional aggregate worker timing for all refits */ );

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
