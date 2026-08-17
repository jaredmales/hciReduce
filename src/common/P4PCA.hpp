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
};

/// Pure all-double principal-component regression for one P4 search pixel.
/** This component forms the smaller of the temporal and predictor Gram matrices without centering or an intercept.
 * A caller-owned LAPACK workspace makes repeated calls efficient and permits one independent workspace per worker
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
     * lower-count columns. The output object must not own storage aliased by either input array.
     */
    static void calculate( P4PCAResult &output,           /**< [out] residuals, mode states, and numerical rank */
                           const matrixT &predictors,     /**< [in] finite T-by-K local predictor matrix */
                           const vectorT &target,         /**< [in] finite T-sample target time series */
                           const std::vector<int> &modes, /**< [in] positive, strictly increasing retained counts */
                           double rankTolerance, /**< [in] finite nonnegative threshold relative to lambdaMax */
                           workspaceT &workspace /**< [in,out] caller-owned, non-shared LAPACK workspace */ );

    /// Calculate temporally centered, mean-preserving truncated-PCA residuals.
    /** Each predictor column and the target are centered over their T samples before fitting. The fitted centered
     * prediction is subtracted from the original uncentered target, so every supported residual preserves the target
     * temporal mean. Centering limits the structural degrees of freedom to min(K,T-1); requested counts above that
     * limit are rejected even when roundoff would make the centering null eigenvalue positive. T must be at least
     * two. On success, \p output is completely replaced. On exception it remains destructible, but its values are
     * unspecified. The output object must not own storage aliased by either input array.
     */
    static void
    calculateCentered( P4PCAResult &output,           /**< [out] residuals, mode states, and numerical rank */
                       const matrixT &predictors,     /**< [in] finite T-by-K uncentered predictor matrix */
                       const vectorT &target,         /**< [in] finite T-sample uncentered target time series */
                       const std::vector<int> &modes, /**< [in] positive, strictly increasing retained counts */
                       double rankTolerance,          /**< [in] finite nonnegative threshold relative to lambdaMax */
                       workspaceT &workspace /**< [in,out] caller-owned, non-shared LAPACK workspace */ );
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
