/** \file P4PCA.cpp
 * \brief Implements the pure numerical predictor used by Pixel Prediction Post-Processing.
 * \author Jared R. Males
 */

#include "P4PCA.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include <mx/math/floatUtils.hpp>

namespace mx
{
namespace improc
{

namespace
{

/// Process-wide test override; production execution leaves this null.
std::atomic<detail::P4PCAEigenSolverT> p4PCAEigenSolverForTesting{ nullptr };

#ifdef HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION
/// Process-wide FP32 test override; production execution leaves this null.
std::atomic<detail::P4PCAFloatEigenSolverT> p4PCAFloatEigenSolverForTesting{ nullptr };
#endif

/// Determine whether every coefficient in an Eigen-like array is finite under fast-math.
template <typename arrayT>
bool p4PCAAllFinite( const arrayT &array /**< [in] array to validate */ )
{
    for( Eigen::Index column = 0; column < array.cols(); ++column )
    {
        for( Eigen::Index row = 0; row < array.rows(); ++row )
        {
            if( !mx::math::isFinite( array( row, column ) ) )
            {
                return false;
            }
        }
    }

    return true;
}

/// Return the P4 acceptance tolerance for a base eigensolver's temporal singular-vector factor.
/** Full-spectrum `SYEVR` orthogonality is dimension-scaled and can marginally exceed mxlib's generic automatic
 * factor-validation bound. This P4-specific bound remains tight enough to reject a materially inconsistent factor.
 */
double p4PCAFactorValidationTolerance( const P4PCA::matrixT &factor /**< [in] temporal singular-vector factor */ )
{
    return 128.0 * std::numeric_limits<double>::epsilon() *
           static_cast<double>( std::max( factor.rows(), factor.cols() ) );
}

/// Measure the maximum absolute error in a temporal factor's orthogonality Gram matrix.
double p4PCAFactorOrthogonalityDefect( const P4PCA::matrixT &factor /**< [in] temporal singular-vector factor */ )
{
    P4PCA::matrixT gram = ( factor.matrix().transpose() * factor.matrix() ).array();
    for( Eigen::Index mode = 0; mode < gram.rows(); ++mode )
    {
        gram( mode, mode ) -= 1.0;
    }
    return gram.abs().maxCoeff();
}

/// Format factor validation failures, including the orthogonality defect when that contract alone failed.
std::string p4PCAFactorValidationMessage( mx::math::svdDeletionStatus status, /**< [in] validation outcome */
                                          const P4PCA::matrixT &factor, /**< [in] temporal singular-vector factor */
                                          double tolerance /**< [in] accepted maximum Gram-matrix error */ )
{
    std::ostringstream message;
    message << "P4PCA downdate temporal factor validation failed with status "
            << mx::math::svdDeletionStatusName( status );
    if( status != mx::math::svdDeletionStatus::factorNotOrthonormal )
    {
        return message.str();
    }

    const double maximumError = p4PCAFactorOrthogonalityDefect( factor );
    message << std::setprecision( 17 ) << ", max|L^T L-I|=" << maximumError << ", acceptedTolerance=" << tolerance;
    return message.str();
}

/// Invoke the installed test solver or the production all-double mxlib eigensolver.
MXLAPACK_INT p4PCAEigenSolve( P4PCA::matrixT &eigenvectors, /**< [out] selected eigenvectors */
                              P4PCA::matrixT &eigenvalues,  /**< [out] selected eigenvalues */
                              P4PCA::matrixT &covariance,   /**< [in] lower-triangular covariance matrix */
                              int modeCount,                /**< [in] number of largest eigenpairs to calculate */
                              P4PCA::workspaceT &workspace /**< [in,out] reusable LAPACK workspace */ )
{
    const detail::P4PCAEigenSolverT testSolver = p4PCAEigenSolverForTesting.load( std::memory_order_acquire );
    if( testSolver )
    {
        return testSolver( eigenvectors, eigenvalues, covariance, modeCount, workspace );
    }

    const int dimension = static_cast<int>( covariance.rows() );
    return mx::math::eigenSYEVR( eigenvectors,
                                 eigenvalues,
                                 covariance,
                                 dimension - modeCount,
                                 dimension,
                                 'L',
                                 &workspace );
}

#ifdef HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION
/// Invoke the installed FP32 test solver or the production FP32 mxlib eigensolver.
MXLAPACK_INT p4PCAEigenSolve( detail::P4PCAFloatMatrixT &eigenvectors, /**< [out] selected eigenvectors */
                              detail::P4PCAFloatMatrixT &eigenvalues,  /**< [out] selected eigenvalues */
                              detail::P4PCAFloatMatrixT &covariance,   /**< [in] lower-triangular covariance matrix */
                              int modeCount, /**< [in] number of largest eigenpairs to calculate */
                              mx::math::syevrMem<float> &workspace /**< [in,out] reusable LAPACK workspace */ )
{
    const detail::P4PCAFloatEigenSolverT testSolver = p4PCAFloatEigenSolverForTesting.load( std::memory_order_acquire );
    if( testSolver )
    {
        return testSolver( eigenvectors, eigenvalues, covariance, modeCount, workspace );
    }

    const int dimension = static_cast<int>( covariance.rows() );
    return mx::math::eigenSYEVR( eigenvectors,
                                 eigenvalues,
                                 covariance,
                                 dimension - modeCount,
                                 dimension,
                                 'L',
                                 &workspace );
}
#endif

/// Dynamic matrix used at one calculation or eigensolver precision.
template <typename ScalarT>
using P4PCADynamicMatrixT = Eigen::Array<ScalarT, Eigen::Dynamic, Eigen::Dynamic>;

/// Dynamic column vector used at one calculation or eigensolver precision.
template <typename ScalarT>
using P4PCADynamicVectorT = Eigen::Array<ScalarT, Eigen::Dynamic, 1>;

/// Invoke the native calculation-precision solver or promote an FP32 Gram matrix for an FP64 eigensolve.
/** A mixed solve returns every nonzero native solver status without inspecting or converting its incomplete outputs.
 * Successful mixed eigenpairs are deliberately cast back to the calculation scalar before rank selection and
 * projection.
 */
template <typename CalculationT, typename EigensolverT>
MXLAPACK_INT p4PCAEigenSolveAdapted(
    P4PCADynamicMatrixT<CalculationT> &eigenvectors, /**< [out] selected calculation-precision eigenvectors */
    P4PCADynamicMatrixT<CalculationT> &eigenvalues,  /**< [out] selected calculation-precision eigenvalues */
    P4PCADynamicMatrixT<CalculationT> &covariance,   /**< [in] lower-triangular covariance matrix */
    int modeCount,                                   /**< [in] number of largest eigenpairs to calculate */
    mx::math::syevrMem<EigensolverT> &workspace /**< [in,out] reusable eigensolver and promotion storage */ )
{
    if constexpr( std::is_same_v<CalculationT, EigensolverT> )
    {
        return p4PCAEigenSolve( eigenvectors, eigenvalues, covariance, modeCount, workspace );
    }
    else
    {
        static_assert( std::is_same_v<CalculationT, float> && std::is_same_v<EigensolverT, double>,
                       "P4PCA supports only the FP32-calculation/FP64-eigensolver mixed policy" );
        workspace.cvd = covariance.template cast<EigensolverT>();
        const MXLAPACK_INT solverStatus =
            p4PCAEigenSolve( workspace.evecsd, workspace.evalsd, workspace.cvd, modeCount, workspace );
        if( solverStatus != 0 )
        {
            return solverStatus;
        }

        eigenvectors = workspace.evecsd.template cast<CalculationT>();
        eigenvalues = workspace.evalsd.template cast<CalculationT>();
        return solverStatus;
    }
}

/// Replace one recoverable factor-downdate attempt with the complete explicit held-out oracle result.
void p4PCAExplicitFallback(
    P4PCAResult &output,                             /**< [in,out] attempted diagnostics, then explicit result */
    const P4PCA::matrixT &predictors,                /**< [in] finite full predictor matrix */
    const P4PCA::vectorT &target,                    /**< [in] finite full target time series */
    const P4TargetExclusions &exclusions,            /**< [in] target-specific deleted rows */
    const std::vector<int> &modes,                   /**< [in] requested retained-mode counts */
    double rankTolerance,                            /**< [in] relative numerical-rank threshold */
    P4PCA::workspaceT &eigensolverWorkspace,         /**< [in,out] reusable explicit-oracle eigensolver storage */
    P4PCAFallbackReason reason,                      /**< [in] recoverable condition that selected the oracle */
    P4PCATiming *timing,                             /**< [in,out] optional attempted plus explicit timing */
    P4PCA::matrixT *probeResiduals = nullptr,        /**< [out] optional explicit frozen-probe result */
    const P4PCA::matrixT *probePredictors = nullptr, /**< [in] optional frozen predictor responses */
    const P4PCA::vectorT *probeTarget = nullptr /**< [in] optional direct frozen target response */ )
{
    const int attemptedBaseRank = output.baseRank;
    const std::size_t attemptedClampCount = output.downdateClampCount;
    const double factorOrthogonalityDefect = output.factorOrthogonalityDefect;
    const double factorOrthogonalityTolerance = output.factorOrthogonalityTolerance;
    P4PCAResult explicitResult;
    P4PCATiming explicitTiming;
    const auto explicitBegin = std::chrono::steady_clock::now();
    if( probeResiduals || probePredictors || probeTarget )
    {
        if( !probeResiduals || !probePredictors || !probeTarget )
        {
            throw std::invalid_argument( "P4PCA explicit fallback requires a complete probe request" );
        }
        P4PCA::calculateHeldOutProbe( explicitResult,
                                      *probeResiduals,
                                      predictors,
                                      target,
                                      *probePredictors,
                                      *probeTarget,
                                      exclusions,
                                      modes,
                                      rankTolerance,
                                      eigensolverWorkspace,
                                      timing ? &explicitTiming : nullptr );
    }
    else
    {
        P4PCA::calculateHeldOut( explicitResult,
                                 predictors,
                                 target,
                                 exclusions,
                                 modes,
                                 rankTolerance,
                                 eigensolverWorkspace,
                                 timing ? &explicitTiming : nullptr );
    }
    explicitResult.baseRank = attemptedBaseRank;
    explicitResult.downdateClampCount = attemptedClampCount;
    explicitResult.explicitFallbackCount = static_cast<std::size_t>( predictors.rows() );
    explicitResult.explicitFallbackReason = reason;
    explicitResult.factorOrthogonalityDefect = factorOrthogonalityDefect;
    explicitResult.factorOrthogonalityTolerance = factorOrthogonalityTolerance;
    output = std::move( explicitResult );
    if( timing )
    {
        timing->gramWorkerSeconds += explicitTiming.gramWorkerSeconds;
        timing->eigensolveWorkerSeconds += explicitTiming.eigensolveWorkerSeconds;
        timing->projectionWorkerSeconds += explicitTiming.projectionWorkerSeconds;
        timing->explicitFallbackWorkerSeconds +=
            std::chrono::duration<double>( std::chrono::steady_clock::now() - explicitBegin ).count();
    }
}

/// Return whether a structured deletion failure can be replaced safely by the complete explicit oracle.
bool p4PCARecoverableStructuredFailure(
    mx::math::svdDeletionStatus status, /**< [in] failed deletion status */
    MXLAPACK_INT lapackInfo /**< [in] underlying LAPACK status, when the solver was entered */ )
{
    switch( status )
    {
    case mx::math::svdDeletionStatus::solverFailure:
        return lapackInfo > 0;
    case mx::math::svdDeletionStatus::nonFiniteOutput:
    case mx::math::svdDeletionStatus::invalidSolverOutput:
    case mx::math::svdDeletionStatus::rescalingOverflow:
    case mx::math::svdDeletionStatus::nonPositiveSemidefinite:
        return true;
    case mx::math::svdDeletionStatus::notComputed:
    case mx::math::svdDeletionStatus::success:
    case mx::math::svdDeletionStatus::successWithClamping:
    case mx::math::svdDeletionStatus::invalidInput:
    case mx::math::svdDeletionStatus::allocationFailure:
    case mx::math::svdDeletionStatus::workspaceQueryFailure:
    case mx::math::svdDeletionStatus::factorNotOrthonormal:
    case mx::math::svdDeletionStatus::unsupportedDeletionCount:
        return false;
    }
    return false;
}

/// Validate one P4PCA request against its path-specific structural degree-of-freedom limit.
template <typename matrixT, typename vectorT>
void p4PCAValidateInputs( const matrixT &predictors,     /**< [in] predictor matrix to validate */
                          const vectorT &target,         /**< [in] target vector to validate */
                          const std::vector<int> &modes, /**< [in] requested retained-mode counts */
                          double rankTolerance,          /**< [in] relative numerical-rank threshold */
                          int maxDOF /**< [in] path-specific structural degree-of-freedom limit */ )
{
    const Eigen::Index sampleCount = predictors.rows();
    const Eigen::Index predictorCount = predictors.cols();

    if( sampleCount <= 0 || predictorCount <= 0 )
    {
        throw std::invalid_argument( "P4PCA requires a nonempty predictor matrix" );
    }

    if( target.rows() != sampleCount )
    {
        throw std::invalid_argument( "P4PCA target length must match the predictor row count" );
    }

    if( modes.empty() )
    {
        throw std::invalid_argument( "P4PCA requires at least one mode count" );
    }

    int previousMode{ 0 };
    for( const int mode : modes )
    {
        if( mode <= 0 )
        {
            throw std::invalid_argument( "P4PCA mode counts must be positive" );
        }

        if( mode <= previousMode )
        {
            throw std::invalid_argument( "P4PCA mode counts must be strictly increasing" );
        }

        if( mode > maxDOF )
        {
            throw std::invalid_argument( "P4PCA mode count exceeds maxDOF" );
        }

        previousMode = mode;
    }

    if( !mx::math::isFinite( rankTolerance ) || rankTolerance < 0 )
    {
        throw std::invalid_argument( "P4PCA rank tolerance must be finite and nonnegative" );
    }

    if( !p4PCAAllFinite( predictors ) || !p4PCAAllFinite( target ) )
    {
        throw std::invalid_argument( "P4PCA predictors and target must be finite" );
    }
}

/// Materialize the sorted complement of one validated deletion set.
void p4PCARetainedRows( std::vector<Eigen::Index> &retainedRows, /**< [out] retained temporal row indices */
                        Eigen::Index sampleCount,                /**< [in] complete temporal row count */
                        std::span<const Eigen::Index> deletedRows /**< [in] sorted rows to omit */ )
{
    retainedRows.clear();
    retainedRows.reserve( static_cast<std::size_t>( sampleCount - static_cast<Eigen::Index>( deletedRows.size() ) ) );
    auto deleted = deletedRows.begin();
    for( Eigen::Index row = 0; row < sampleCount; ++row )
    {
        if( deleted != deletedRows.end() && *deleted == row )
        {
            ++deleted;
            continue;
        }
        retainedRows.push_back( row );
    }
}

/// Subtract each column's calculation-precision temporal mean in place and optionally return those means.
template <typename arrayT>
void p4PCACenterColumns( arrayT &array, /**< [in,out] array whose columns are centered */
                         arrayT *columnMeans = nullptr /**< [out] optional one-row array of removed means */ )
{
    using scalarT = typename arrayT::Scalar;

    if( columnMeans )
    {
        columnMeans->resize( 1, array.cols() );
        columnMeans->setZero();
    }

    for( Eigen::Index column = 0; column < array.cols(); ++column )
    {
        scalarT scale{ 0 };
        for( Eigen::Index row = 0; row < array.rows(); ++row )
        {
            scale = std::max( scale, std::abs( array( row, column ) ) );
        }

        if( scale == 0 )
        {
            continue;
        }

        scalarT scaledSum{ 0 };
        for( Eigen::Index row = 0; row < array.rows(); ++row )
        {
            scaledSum += array( row, column ) / scale;
        }

        const scalarT scaledMean =
            std::clamp( scaledSum / static_cast<scalarT>( array.rows() ), scalarT{ -1 }, scalarT{ 1 } );
        const scalarT mean = scale * scaledMean;
        if( columnMeans )
        {
            ( *columnMeans )( 0, column ) = mean;
        }
        array.col( column ) -= mean;
    }
}

/// \cond P4PCA_detail

/// Translation-unit-local scalar-policy seam for the direct and centered P4 PCA paths.
/** Normal builds instantiate the all-double oracle and production mixed policy. The experimental precision option
 * additionally instantiates the all-FP32 policy without exposing this implementation type.
 *
 * \tparam CalculationT scalar used for normal equations, projections, and centering
 * \tparam EigensolverT scalar used by the reusable symmetric-eigensolver workspace
 */
template <typename CalculationT, typename EigensolverT>
struct P4PCAKernel
{
    /// Dynamic matrix using the calculation scalar.
    using matrixT = Eigen::Array<CalculationT, Eigen::Dynamic, Eigen::Dynamic>;

    /// Dynamic column vector using the calculation scalar.
    using vectorT = Eigen::Array<CalculationT, Eigen::Dynamic, 1>;

    /// Reusable LAPACK workspace using the eigensolver scalar.
    using workspaceT = mx::math::syevrMem<EigensolverT>;

    /// Calculate direct truncated-PCA residuals through this scalar policy.
    static void calculate( P4PCAResult &output,           /**< [out] residuals, mode states, and numerical rank */
                           const matrixT &predictors,     /**< [in] finite T-by-K local predictor matrix */
                           const vectorT &target,         /**< [in] finite T-sample target time series */
                           const std::vector<int> &modes, /**< [in] positive, strictly increasing retained counts */
                           double rankTolerance,  /**< [in] finite nonnegative threshold relative to lambdaMax */
                           workspaceT &workspace, /**< [in,out] caller-owned, non-shared LAPACK workspace */
                           P4PCATiming *timing,   /**< [out] optional per-call worker timing */
                           P4PCA::matrixT *coefficients /**< [out] optional K-by-mode double coefficients */ );

    /// Calculate a centered fit while preserving the caller's predictor matrix.
    static void
    calculateCentered( P4PCAResult &output,           /**< [out] residuals, mode states, and numerical rank */
                       const matrixT &predictors,     /**< [in] finite T-by-K uncentered predictor matrix */
                       const vectorT &target,         /**< [in] finite T-sample uncentered target time series */
                       const std::vector<int> &modes, /**< [in] positive, strictly increasing retained counts */
                       double rankTolerance,          /**< [in] finite nonnegative threshold relative to lambdaMax */
                       workspaceT &workspace,         /**< [in,out] caller-owned, non-shared LAPACK workspace */
                       P4PCATiming *timing,           /**< [out] optional per-call worker timing */
                       P4PCA::matrixT *coefficients /**< [out] optional K-by-mode double coefficients */ );

    /// Calculate a centered fit using the caller's predictor matrix as workspace.
    static void
    calculateCenteredInPlace( P4PCAResult &output,   /**< [out] residuals, mode states, and numerical rank */
                              matrixT &predictors,   /**< [in,out] finite predictors, replaced by centered values */
                              const vectorT &target, /**< [in] finite T-sample uncentered target time series */
                              const std::vector<int> &modes, /**< [in] positive, strictly increasing retained counts */
                              double rankTolerance,  /**< [in] finite nonnegative threshold relative to lambdaMax */
                              workspaceT &workspace, /**< [in,out] caller-owned, non-shared LAPACK workspace */
                              P4PCATiming *timing,   /**< [out] optional per-call worker timing */
                              P4PCA::matrixT *coefficients /**< [out] optional K-by-mode double coefficients */ );

  private:
    /// Calculate one already-validated request using selected fit and residual targets.
    static void
    calculateValidated( P4PCAResult &output,           /**< [out] completed regression result */
                        const matrixT &predictors,     /**< [in] predictors used to form the fit */
                        const vectorT &fitTarget,      /**< [in] target used to calculate projections */
                        const vectorT &residualTarget, /**< [in] initial target for output residuals */
                        const std::vector<int> &modes, /**< [in] requested retained-mode counts */
                        double rankTolerance,          /**< [in] relative numerical-rank threshold */
                        int maxDOF, /**< [in] structural degree-of-freedom limit and eigenpair count */
                        const matrixT *predictorMeans, /**< [in] optional means for uncentered application */
                        workspaceT &workspace,         /**< [in,out] reusable LAPACK workspace */
                        P4PCATiming *timing,           /**< [out] optional per-call worker timing */
                        P4PCA::matrixT *coefficients /**< [out] optional predictor-space coefficient vectors */ );
};

template <typename CalculationT, typename EigensolverT>
void P4PCAKernel<CalculationT, EigensolverT>::calculateValidated( P4PCAResult &output,
                                                                  const matrixT &predictors,
                                                                  const vectorT &fitTarget,
                                                                  const vectorT &residualTarget,
                                                                  const std::vector<int> &modes,
                                                                  double rankTolerance,
                                                                  int maxDOF,
                                                                  const matrixT *predictorMeans,
                                                                  workspaceT &workspace,
                                                                  P4PCATiming *timing,
                                                                  P4PCA::matrixT *coefficients )
{
    if( timing )
    {
        *timing = P4PCATiming{};
    }

    const auto secondsSince = []( const std::chrono::steady_clock::time_point &start )
    { return std::chrono::duration<double>( std::chrono::steady_clock::now() - start ).count(); };
    const Eigen::Index sampleCount = predictors.rows();
    const Eigen::Index predictorCount = predictors.cols();
    const bool useTemporalGram = sampleCount <= predictorCount;
    matrixT gram;
    vectorT crossProduct;
    const auto gramStart = std::chrono::steady_clock::now();
    if( useTemporalGram )
    {
        gram = ( predictors.matrix() * predictors.matrix().transpose() ).array();
    }
    else
    {
        gram = ( predictors.matrix().transpose() * predictors.matrix() ).array();
        crossProduct = ( predictors.matrix().transpose() * fitTarget.matrix() ).array();
    }
    if( timing )
    {
        timing->gramWorkerSeconds = secondsSince( gramStart );
    }

    if( !p4PCAAllFinite( gram ) || ( !useTemporalGram && !p4PCAAllFinite( crossProduct ) ) )
    {
        throw std::runtime_error( "P4PCA normal equations produced nonfinite values" );
    }

    vectorT predictorMeanProjection;
    if( predictorMeans )
    {
        if( useTemporalGram )
        {
            predictorMeanProjection = ( predictors.matrix() * predictorMeans->matrix().transpose() ).array();
        }
    }

    matrixT eigenvectors;
    matrixT eigenvalues;
    const auto eigensolveStart = std::chrono::steady_clock::now();
    const MXLAPACK_INT solverStatus =
        p4PCAEigenSolveAdapted<CalculationT, EigensolverT>( eigenvectors, eigenvalues, gram, maxDOF, workspace );
    if( solverStatus != 0 )
    {
        throw std::runtime_error( "P4PCA eigensolver failed with status " + std::to_string( solverStatus ) );
    }

    if( eigenvectors.rows() != gram.rows() || eigenvectors.cols() != maxDOF || eigenvalues.rows() < maxDOF ||
        eigenvalues.cols() < 1 )
    {
        throw std::runtime_error( "P4PCA eigensolver returned invalid dimensions" );
    }

    if( !p4PCAAllFinite( eigenvectors ) || !p4PCAAllFinite( eigenvalues.block( 0, 0, maxDOF, 1 ) ) )
    {
        throw std::runtime_error( "P4PCA eigensolver returned nonfinite results" );
    }

    for( int eigenIndex = 1; eigenIndex < maxDOF; ++eigenIndex )
    {
        if( eigenvalues( eigenIndex ) < eigenvalues( eigenIndex - 1 ) )
        {
            throw std::runtime_error( "P4PCA eigensolver returned unsorted eigenvalues" );
        }
    }

    const CalculationT largestEigenvalue = eigenvalues( maxDOF - 1 );
    int numericalRank{ 0 };
    if( largestEigenvalue > 0 )
    {
        const long double rankThreshold =
            static_cast<long double>( rankTolerance ) * static_cast<long double>( largestEigenvalue );
        for( int eigenIndex = 0; eigenIndex < maxDOF; ++eigenIndex )
        {
            if( static_cast<long double>( eigenvalues( eigenIndex ) ) > rankThreshold )
            {
                ++numericalRank;
            }
        }
    }
    if( timing )
    {
        timing->eigensolveWorkerSeconds = secondsSince( eigensolveStart );
    }

    output.residuals.resize( sampleCount, modes.size() );
    output.residuals.setConstant( std::numeric_limits<double>::quiet_NaN() );
    output.modeStatus.assign( modes.size(), P4PCAModeStatus::rankInsufficient );
    output.numericalRank = numericalRank;
    output.baseRank = 0;
    output.numericalRankCapped = false;
    output.downdateClampCount = 0;
    output.explicitFallbackCount = 0;
    output.explicitFallbackReason = P4PCAFallbackReason::none;
    output.factorOrthogonalityDefect = 0;
    output.factorOrthogonalityTolerance = 0;
    if( coefficients )
    {
        coefficients->resize( predictorCount, modes.size() );
        coefficients->setConstant( std::numeric_limits<double>::quiet_NaN() );
    }

    const auto unsupported = std::upper_bound( modes.begin(), modes.end(), numericalRank );
    if( unsupported == modes.begin() )
    {
        return;
    }

    const int largestSupportedMode = *( unsupported - 1 );
    const auto projectionStart = std::chrono::steady_clock::now();
    vectorT residual = residualTarget;
    vectorT accumulatedCoefficients;
    vectorT coefficientUpdate;
    if( coefficients )
    {
        accumulatedCoefficients.resize( predictorCount );
        accumulatedCoefficients.setZero();
        if( useTemporalGram )
        {
            coefficientUpdate.resize( predictorCount );
        }
    }
    std::size_t outputIndex{ 0 };

    for( int retainedMode = 1; retainedMode <= largestSupportedMode; ++retainedMode )
    {
        const int eigenIndex = maxDOF - retainedMode;
        const CalculationT eigenvalue = eigenvalues( eigenIndex );
        const vectorT eigenvector = eigenvectors.col( eigenIndex );
        CalculationT coefficient{ 0 };
        CalculationT predictionMean{ 0 };
        vectorT projectedMode;
        if( useTemporalGram )
        {
            coefficient = eigenvector.matrix().dot( fitTarget.matrix() );
            projectedMode = eigenvector * coefficient;
            if( coefficients )
            {
                coefficientUpdate = ( predictors.matrix().transpose() * eigenvector.matrix() ).array();
                coefficientUpdate *= coefficient / eigenvalue;
                accumulatedCoefficients += coefficientUpdate;
            }
            if( predictorMeans )
            {
                predictionMean =
                    eigenvector.matrix().dot( predictorMeanProjection.matrix() ) * coefficient / eigenvalue;
            }
        }
        else
        {
            coefficient = eigenvector.matrix().dot( crossProduct.matrix() ) / eigenvalue;
            projectedMode = ( predictors.matrix() * eigenvector.matrix() ).array();
            projectedMode *= coefficient;
            if( coefficients )
            {
                accumulatedCoefficients += eigenvector * coefficient;
            }
            if( predictorMeans )
            {
                predictionMean = ( predictorMeans->matrix() * eigenvector.matrix() )( 0, 0 ) * coefficient;
            }
        }

        if( !mx::math::isFinite( coefficient ) || !mx::math::isFinite( predictionMean ) ||
            !p4PCAAllFinite( projectedMode ) || ( coefficients && !p4PCAAllFinite( accumulatedCoefficients ) ) )
        {
            throw std::runtime_error( "P4PCA mode projection produced nonfinite values" );
        }

        if( predictorMeans )
        {
            projectedMode += predictionMean;
        }

        residual -= projectedMode;
        if( !p4PCAAllFinite( residual ) )
        {
            throw std::runtime_error( "P4PCA residual calculation produced nonfinite values" );
        }

        if( outputIndex < modes.size() && modes[outputIndex] == retainedMode )
        {
            if constexpr( std::is_same_v<CalculationT, double> )
            {
                output.residuals.col( outputIndex ) = residual;
            }
            else
            {
                output.residuals.col( outputIndex ) = residual.template cast<double>();
            }
            output.modeStatus[outputIndex] = P4PCAModeStatus::rankSupported;
            if( coefficients )
            {
                if constexpr( std::is_same_v<CalculationT, double> )
                {
                    coefficients->col( outputIndex ) = accumulatedCoefficients;
                }
                else
                {
                    coefficients->col( outputIndex ) = accumulatedCoefficients.template cast<double>();
                }
            }
            ++outputIndex;
        }
    }
    if( timing )
    {
        timing->projectionWorkerSeconds = secondsSince( projectionStart );
    }
}

template <typename CalculationT, typename EigensolverT>
void P4PCAKernel<CalculationT, EigensolverT>::calculate( P4PCAResult &output,
                                                         const matrixT &predictors,
                                                         const vectorT &target,
                                                         const std::vector<int> &modes,
                                                         double rankTolerance,
                                                         workspaceT &workspace,
                                                         P4PCATiming *timing,
                                                         P4PCA::matrixT *coefficients )
{
    output.sampleValidity.resize( 0, 0 );
    const Eigen::Index sampleCount = predictors.rows();
    const Eigen::Index predictorCount = predictors.cols();
    const int maxDOF = static_cast<int>( std::min( sampleCount, predictorCount ) );
    p4PCAValidateInputs( predictors, target, modes, rankTolerance, maxDOF );
    calculateValidated( output,
                        predictors,
                        target,
                        target,
                        modes,
                        rankTolerance,
                        maxDOF,
                        nullptr,
                        workspace,
                        timing,
                        coefficients );
}

template <typename CalculationT, typename EigensolverT>
void P4PCAKernel<CalculationT, EigensolverT>::calculateCentered( P4PCAResult &output,
                                                                 const matrixT &predictors,
                                                                 const vectorT &target,
                                                                 const std::vector<int> &modes,
                                                                 double rankTolerance,
                                                                 workspaceT &workspace,
                                                                 P4PCATiming *timing,
                                                                 P4PCA::matrixT *coefficients )
{
    matrixT centeredPredictors = predictors;
    calculateCenteredInPlace( output,
                              centeredPredictors,
                              target,
                              modes,
                              rankTolerance,
                              workspace,
                              timing,
                              coefficients );
}

template <typename CalculationT, typename EigensolverT>
void P4PCAKernel<CalculationT, EigensolverT>::calculateCenteredInPlace( P4PCAResult &output,
                                                                        matrixT &predictors,
                                                                        const vectorT &target,
                                                                        const std::vector<int> &modes,
                                                                        double rankTolerance,
                                                                        workspaceT &workspace,
                                                                        P4PCATiming *timing,
                                                                        P4PCA::matrixT *coefficients )
{
    output.sampleValidity.resize( 0, 0 );
    const Eigen::Index sampleCount = predictors.rows();
    const Eigen::Index predictorCount = predictors.cols();
    if( sampleCount == 1 )
    {
        throw std::invalid_argument( "P4PCA centered regression requires at least two samples" );
    }

    const int maxDOF = static_cast<int>( std::min( predictorCount, sampleCount - 1 ) );
    p4PCAValidateInputs( predictors, target, modes, rankTolerance, maxDOF );

    vectorT centeredTarget = target;
    matrixT predictorMeans;
    p4PCACenterColumns( predictors, &predictorMeans );
    p4PCACenterColumns( centeredTarget );
    if( !p4PCAAllFinite( predictors ) || !p4PCAAllFinite( centeredTarget ) )
    {
        throw std::runtime_error( "P4PCA temporal centering produced nonfinite values" );
    }

    calculateValidated( output,
                        predictors,
                        centeredTarget,
                        target,
                        modes,
                        rankTolerance,
                        maxDOF,
                        &predictorMeans,
                        workspace,
                        timing,
                        coefficients );
}

template struct P4PCAKernel<double, double>;
template struct P4PCAKernel<float, double>;

#ifdef HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION
template struct P4PCAKernel<float, float>;
#endif

/// \endcond
} // namespace

P4TargetExclusions P4TargetExclusions::fromSpans( Eigen::Index sampleCount, const std::vector<P4ExclusionSpan> &spans )
{
    if( sampleCount <= 0 || spans.size() != static_cast<std::size_t>( sampleCount ) )
    {
        throw std::invalid_argument( "P4 exclusion spans must provide one pattern per temporal row" );
    }

    P4TargetExclusions result;
    result.m_sampleCount = sampleCount;
    result.m_storage = P4ExclusionStorage::contiguousSpans;
    result.m_spans = spans;
    for( Eigen::Index target = 0; target < sampleCount; ++target )
    {
        const P4ExclusionSpan &span = spans[static_cast<std::size_t>( target )];
        if( span.begin < 0 || span.begin > target || span.end <= target || span.end > sampleCount )
        {
            throw std::invalid_argument( "P4 exclusion span must be in range, nonempty, and contain its target" );
        }
        result.m_maximumDeleted = std::max( result.m_maximumDeleted, span.end - span.begin );
    }
    return result;
}

P4TargetExclusions P4TargetExclusions::fromExplicit( Eigen::Index sampleCount,
                                                     const std::vector<std::vector<Eigen::Index>> &indices )
{
    if( sampleCount <= 0 || indices.size() != static_cast<std::size_t>( sampleCount ) )
    {
        throw std::invalid_argument( "P4 explicit exclusions must provide one pattern per temporal row" );
    }

    std::size_t totalIndices{ 0 };
    for( const auto &targetIndices : indices )
    {
        if( targetIndices.size() > std::numeric_limits<std::size_t>::max() - totalIndices )
        {
            throw std::length_error( "P4 explicit exclusion storage exceeds size_t range" );
        }
        totalIndices += targetIndices.size();
    }
    if( totalIndices > static_cast<std::size_t>( std::numeric_limits<Eigen::Index>::max() ) )
    {
        throw std::length_error( "P4 explicit exclusion storage exceeds Eigen index range" );
    }

    P4TargetExclusions result;
    result.m_sampleCount = sampleCount;
    result.m_storage = P4ExclusionStorage::explicitIndices;
    result.m_offsets.reserve( static_cast<std::size_t>( sampleCount ) + 1 );
    result.m_indices.reserve( totalIndices );
    result.m_offsets.push_back( 0 );
    for( Eigen::Index target = 0; target < sampleCount; ++target )
    {
        const auto &targetIndices = indices[static_cast<std::size_t>( target )];
        bool containsTarget{ false };
        Eigen::Index previous{ -1 };
        for( const Eigen::Index row : targetIndices )
        {
            if( row < 0 || row >= sampleCount || row <= previous )
            {
                throw std::invalid_argument( "P4 explicit exclusion rows must be sorted, unique, and in range" );
            }
            containsTarget = containsTarget || row == target;
            previous = row;
        }
        if( !containsTarget )
        {
            throw std::invalid_argument( "P4 explicit exclusion pattern must contain its target row" );
        }
        result.m_indices.insert( result.m_indices.end(), targetIndices.begin(), targetIndices.end() );
        result.m_offsets.push_back( static_cast<Eigen::Index>( result.m_indices.size() ) );
        result.m_maximumDeleted =
            std::max( result.m_maximumDeleted, static_cast<Eigen::Index>( targetIndices.size() ) );
    }
    return result;
}

bool P4TargetExclusions::empty() const noexcept
{
    return m_sampleCount == 0;
}

Eigen::Index P4TargetExclusions::sampleCount() const noexcept
{
    return m_sampleCount;
}

Eigen::Index P4TargetExclusions::targetCount() const noexcept
{
    return m_sampleCount;
}

P4ExclusionStorage P4TargetExclusions::storage() const noexcept
{
    return m_storage;
}

Eigen::Index P4TargetExclusions::deletedCount( Eigen::Index target ) const
{
    if( target < 0 || target >= m_sampleCount )
    {
        throw std::out_of_range( "P4 exclusion target is out of range" );
    }
    if( m_storage == P4ExclusionStorage::contiguousSpans )
    {
        const P4ExclusionSpan &span = m_spans[static_cast<std::size_t>( target )];
        return span.end - span.begin;
    }
    return m_offsets[static_cast<std::size_t>( target ) + 1] - m_offsets[static_cast<std::size_t>( target )];
}

Eigen::Index P4TargetExclusions::maximumDeleted() const noexcept
{
    return m_maximumDeleted;
}

std::span<const Eigen::Index> P4TargetExclusions::deletedRows( std::vector<Eigen::Index> &spanScratch,
                                                               Eigen::Index target ) const
{
    if( target < 0 || target >= m_sampleCount )
    {
        throw std::out_of_range( "P4 exclusion target is out of range" );
    }
    if( m_storage == P4ExclusionStorage::contiguousSpans )
    {
        const P4ExclusionSpan &span = m_spans[static_cast<std::size_t>( target )];
        spanScratch.resize( static_cast<std::size_t>( span.end - span.begin ) );
        std::iota( spanScratch.begin(), spanScratch.end(), span.begin );
        return spanScratch;
    }

    const Eigen::Index begin = m_offsets[static_cast<std::size_t>( target )];
    const Eigen::Index end = m_offsets[static_cast<std::size_t>( target ) + 1];
    return std::span<const Eigen::Index>( m_indices.data() + begin, static_cast<std::size_t>( end - begin ) );
}

std::size_t P4TargetExclusions::storageBytes() const noexcept
{
    const long double bytes = static_cast<long double>( m_spans.capacity() ) * sizeof( P4ExclusionSpan ) +
                              static_cast<long double>( m_offsets.capacity() ) * sizeof( Eigen::Index ) +
                              static_cast<long double>( m_indices.capacity() ) * sizeof( Eigen::Index );
    if( bytes > static_cast<long double>( std::numeric_limits<std::size_t>::max() ) )
    {
        return std::numeric_limits<std::size_t>::max();
    }
    return static_cast<std::size_t>( bytes );
}

namespace detail
{

/// \cond P4PCA_test_detail
void p4PCASetEigenSolverForTesting( P4PCAEigenSolverT solver )
{
    p4PCAEigenSolverForTesting.store( solver, std::memory_order_release );
}

void p4PCAResetEigenSolverForTesting()
{
    p4PCAEigenSolverForTesting.store( nullptr, std::memory_order_release );
}

void p4PCACalculateMixed( P4PCAResult &output,
                          const P4PCA::matrixT &predictors,
                          const P4PCA::vectorT &target,
                          const std::vector<int> &modes,
                          double rankTolerance,
                          P4PCAMixedWorkspace &workspace,
                          P4PCATiming *timing,
                          P4PCA::matrixT *coefficients )
{
    workspace.floatPredictors = predictors.cast<float>();
    workspace.floatTarget = target.cast<float>();
    P4PCAKernel<float, double>::calculate( output,
                                           workspace.floatPredictors,
                                           workspace.floatTarget,
                                           modes,
                                           rankTolerance,
                                           workspace.doubleEigensolver,
                                           timing,
                                           coefficients );
}

void p4PCACalculateCenteredMixed( P4PCAResult &output,
                                  const P4PCA::matrixT &predictors,
                                  const P4PCA::vectorT &target,
                                  const std::vector<int> &modes,
                                  double rankTolerance,
                                  P4PCAMixedWorkspace &workspace,
                                  P4PCATiming *timing,
                                  P4PCA::matrixT *coefficients )
{
    workspace.floatPredictors = predictors.cast<float>();
    workspace.floatTarget = target.cast<float>();
    P4PCAKernel<float, double>::calculateCentered( output,
                                                   workspace.floatPredictors,
                                                   workspace.floatTarget,
                                                   modes,
                                                   rankTolerance,
                                                   workspace.doubleEigensolver,
                                                   timing,
                                                   coefficients );
}

void p4PCACalculateCenteredInPlaceMixed( P4PCAResult &output,
                                         P4PCA::matrixT &predictors,
                                         const P4PCA::vectorT &target,
                                         const std::vector<int> &modes,
                                         double rankTolerance,
                                         P4PCAMixedWorkspace &workspace,
                                         P4PCATiming *timing,
                                         P4PCA::matrixT *coefficients )
{
    workspace.floatPredictors = predictors.cast<float>();
    workspace.floatTarget = target.cast<float>();
    if( !p4PCAAllFinite( workspace.floatPredictors ) || !p4PCAAllFinite( workspace.floatTarget ) )
    {
        throw std::invalid_argument( "P4PCA values and rank tolerance must be finite and valid" );
    }
    try
    {
        P4PCAKernel<float, double>::calculateCenteredInPlace( output,
                                                              workspace.floatPredictors,
                                                              workspace.floatTarget,
                                                              modes,
                                                              rankTolerance,
                                                              workspace.doubleEigensolver,
                                                              timing,
                                                              coefficients );
    }
    catch( ... )
    {
        predictors = workspace.floatPredictors.cast<double>();
        throw;
    }
    predictors = workspace.floatPredictors.cast<double>();
}

#ifdef HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION
void p4PCASetFloatEigenSolverForTesting( P4PCAFloatEigenSolverT solver )
{
    p4PCAFloatEigenSolverForTesting.store( solver, std::memory_order_release );
}

void p4PCAResetFloatEigenSolverForTesting()
{
    p4PCAFloatEigenSolverForTesting.store( nullptr, std::memory_order_release );
}

void p4PCACalculateExperimental( P4PCAResult &output,
                                 const P4PCA::matrixT &predictors,
                                 const P4PCA::vectorT &target,
                                 const std::vector<int> &modes,
                                 double rankTolerance,
                                 P4PCAPrecisionPolicy precisionPolicy,
                                 P4PCAExperimentalWorkspace &workspace,
                                 P4PCATiming *timing,
                                 P4PCA::matrixT *coefficients )
{
    switch( precisionPolicy )
    {
    case P4PCAPrecisionPolicy::doubleDouble:
        P4PCAKernel<double, double>::calculate( output,
                                                predictors,
                                                target,
                                                modes,
                                                rankTolerance,
                                                workspace.doubleEigensolver,
                                                timing,
                                                coefficients );
        return;
    case P4PCAPrecisionPolicy::floatDouble:
        p4PCACalculateMixed( output, predictors, target, modes, rankTolerance, workspace, timing, coefficients );
        return;
    case P4PCAPrecisionPolicy::floatFloat:
        workspace.floatPredictors = predictors.cast<float>();
        workspace.floatTarget = target.cast<float>();
        P4PCAKernel<float, float>::calculate( output,
                                              workspace.floatPredictors,
                                              workspace.floatTarget,
                                              modes,
                                              rankTolerance,
                                              workspace.floatEigensolver,
                                              timing,
                                              coefficients );
        return;
    }

    throw std::invalid_argument( "P4PCA experimental precision policy is invalid" );
}

void p4PCACalculateCenteredExperimental( P4PCAResult &output,
                                         const P4PCA::matrixT &predictors,
                                         const P4PCA::vectorT &target,
                                         const std::vector<int> &modes,
                                         double rankTolerance,
                                         P4PCAPrecisionPolicy precisionPolicy,
                                         P4PCAExperimentalWorkspace &workspace,
                                         P4PCATiming *timing,
                                         P4PCA::matrixT *coefficients )
{
    switch( precisionPolicy )
    {
    case P4PCAPrecisionPolicy::doubleDouble:
        P4PCAKernel<double, double>::calculateCentered( output,
                                                        predictors,
                                                        target,
                                                        modes,
                                                        rankTolerance,
                                                        workspace.doubleEigensolver,
                                                        timing,
                                                        coefficients );
        return;
    case P4PCAPrecisionPolicy::floatDouble:
        p4PCACalculateCenteredMixed( output,
                                     predictors,
                                     target,
                                     modes,
                                     rankTolerance,
                                     workspace,
                                     timing,
                                     coefficients );
        return;
    case P4PCAPrecisionPolicy::floatFloat:
        workspace.floatPredictors = predictors.cast<float>();
        workspace.floatTarget = target.cast<float>();
        P4PCAKernel<float, float>::calculateCentered( output,
                                                      workspace.floatPredictors,
                                                      workspace.floatTarget,
                                                      modes,
                                                      rankTolerance,
                                                      workspace.floatEigensolver,
                                                      timing,
                                                      coefficients );
        return;
    }

    throw std::invalid_argument( "P4PCA experimental precision policy is invalid" );
}

void p4PCACalculateCenteredInPlaceExperimental( P4PCAResult &output,
                                                P4PCA::matrixT &predictors,
                                                const P4PCA::vectorT &target,
                                                const std::vector<int> &modes,
                                                double rankTolerance,
                                                P4PCAPrecisionPolicy precisionPolicy,
                                                P4PCAExperimentalWorkspace &workspace,
                                                P4PCATiming *timing,
                                                P4PCA::matrixT *coefficients )
{
    switch( precisionPolicy )
    {
    case P4PCAPrecisionPolicy::doubleDouble:
        P4PCAKernel<double, double>::calculateCenteredInPlace( output,
                                                               predictors,
                                                               target,
                                                               modes,
                                                               rankTolerance,
                                                               workspace.doubleEigensolver,
                                                               timing,
                                                               coefficients );
        return;
    case P4PCAPrecisionPolicy::floatDouble:
        p4PCACalculateCenteredInPlaceMixed( output,
                                            predictors,
                                            target,
                                            modes,
                                            rankTolerance,
                                            workspace,
                                            timing,
                                            coefficients );
        return;
    case P4PCAPrecisionPolicy::floatFloat:
        workspace.floatPredictors = predictors.cast<float>();
        workspace.floatTarget = target.cast<float>();
        if( !p4PCAAllFinite( workspace.floatPredictors ) || !p4PCAAllFinite( workspace.floatTarget ) )
        {
            throw std::invalid_argument( "P4PCA values and rank tolerance must be finite and valid" );
        }
        try
        {
            P4PCAKernel<float, float>::calculateCenteredInPlace( output,
                                                                 workspace.floatPredictors,
                                                                 workspace.floatTarget,
                                                                 modes,
                                                                 rankTolerance,
                                                                 workspace.floatEigensolver,
                                                                 timing,
                                                                 coefficients );
        }
        catch( ... )
        {
            predictors = workspace.floatPredictors.cast<double>();
            throw;
        }
        predictors = workspace.floatPredictors.cast<double>();
        return;
    }

    throw std::invalid_argument( "P4PCA experimental precision policy is invalid" );
}
#endif
/// \endcond

} // namespace detail

bool P4PCAResult::sampleSupported( Eigen::Index sample, std::size_t mode ) const
{
    if( sample < 0 || sample >= residuals.rows() || mode >= modeStatus.size() )
    {
        throw std::out_of_range( "P4PCA result sample or mode is out of range" );
    }
    if( sampleValidity.size() == 0 )
    {
        return modeStatus[mode] == P4PCAModeStatus::rankSupported;
    }
    if( sampleValidity.rows() != residuals.rows() || sampleValidity.cols() != residuals.cols() )
    {
        throw std::logic_error( "P4PCA sample-validity dimensions do not match residuals" );
    }
    return sampleValidity( sample, static_cast<Eigen::Index>( mode ) ) != 0;
}

void P4PCA::calculate( P4PCAResult &output,
                       const matrixT &predictors,
                       const vectorT &target,
                       const std::vector<int> &modes,
                       double rankTolerance,
                       workspaceT &workspace,
                       P4PCATiming *timing,
                       matrixT *coefficients )
{
    P4PCAKernel<double,
                double>::calculate( output, predictors, target, modes, rankTolerance, workspace, timing, coefficients );
}

namespace
{

/// Reject a frozen-probe output that is also one of its immutable matrix inputs.
void p4PCAValidateProbeOutputAliasing( P4PCA::matrixT &probeResiduals,   /**< [out] proposed probe output */
                                       const P4PCA::matrixT &predictors, /**< [in] science predictors */
                                       const P4PCA::matrixT &probePredictors /**< [in] frozen-probe predictors */ )
{
    if( &probeResiduals == &predictors || &probeResiduals == &probePredictors )
    {
        throw std::invalid_argument( "P4PCA held-out probe output must not alias predictor inputs" );
    }
}

/// Calculate explicit target-held-out science residuals and an optional frozen linear-probe response.
template <typename CalculationT, typename EigensolverT>
void p4PCACalculateHeldOut(
    P4PCAResult &output,                                      /**< [out] held-out science result */
    P4PCA::matrixT *probeResiduals,                           /**< [out] optional FP64 held-out probe responses */
    const P4PCADynamicMatrixT<CalculationT> &predictors,      /**< [in] full predictor matrix */
    const P4PCADynamicVectorT<CalculationT> &target,          /**< [in] full regression target */
    const P4PCADynamicMatrixT<CalculationT> *probePredictors, /**< [in] optional probe predictor responses */
    const P4PCADynamicVectorT<CalculationT> *probeTarget,     /**< [in] optional direct probe response */
    const P4TargetExclusions &exclusions,                     /**< [in] target-specific deleted rows */
    const std::vector<int> &modes,                            /**< [in] requested retained-mode counts */
    double rankTolerance,                                     /**< [in] relative numerical-rank threshold */
    mx::math::syevrMem<EigensolverT> &workspace,              /**< [in,out] reusable eigensolver workspace */
    P4PCATiming *timing /**< [out] optional aggregate timing */ )
{
    using matrixT = P4PCADynamicMatrixT<CalculationT>;
    using vectorT = P4PCADynamicVectorT<CalculationT>;

    const Eigen::Index sampleCount = predictors.rows();
    const Eigen::Index predictorCount = predictors.cols();
    if( sampleCount <= 0 || predictorCount <= 0 || target.rows() != sampleCount ||
        exclusions.sampleCount() != sampleCount || exclusions.targetCount() != sampleCount )
    {
        throw std::invalid_argument( "P4PCA held-out inputs must have consistent nonempty dimensions" );
    }
    if( modes.empty() )
    {
        throw std::invalid_argument( "P4PCA held-out regression requires at least one mode count" );
    }
    int previousMode{ 0 };
    for( const int mode : modes )
    {
        if( mode <= previousMode )
        {
            throw std::invalid_argument( "P4PCA held-out mode counts must be positive and strictly increasing" );
        }
        previousMode = mode;
    }
    if( !mx::math::isFinite( rankTolerance ) || rankTolerance < 0 || !p4PCAAllFinite( predictors ) ||
        !p4PCAAllFinite( target ) )
    {
        throw std::invalid_argument( "P4PCA held-out values and rank tolerance must be finite and valid" );
    }
    const bool calculateProbe = probeResiduals || probePredictors || probeTarget;
    if( calculateProbe &&
        ( !probeResiduals || !probePredictors || !probeTarget || probePredictors->rows() <= 0 ||
          probePredictors->cols() != predictorCount || probeTarget->rows() != probePredictors->rows() ||
          !p4PCAAllFinite( *probePredictors ) || !p4PCAAllFinite( *probeTarget ) ) )
    {
        throw std::invalid_argument( "P4PCA held-out probe inputs must have consistent finite dimensions" );
    }

    output.residuals.resize( sampleCount, static_cast<Eigen::Index>( modes.size() ) );
    output.residuals.setConstant( std::numeric_limits<double>::quiet_NaN() );
    output.sampleValidity.resize( sampleCount, static_cast<Eigen::Index>( modes.size() ) );
    output.sampleValidity.setZero();
    output.modeStatus.assign( modes.size(), P4PCAModeStatus::rankInsufficient );
    output.numericalRank = std::numeric_limits<int>::max();
    output.baseRank = 0;
    output.numericalRankCapped = false;
    output.downdateClampCount = 0;
    output.explicitFallbackCount = 0;
    output.explicitFallbackReason = P4PCAFallbackReason::none;
    output.factorOrthogonalityDefect = 0;
    output.factorOrthogonalityTolerance = 0;
    if( calculateProbe )
    {
        probeResiduals->resize( probePredictors->rows(), sampleCount * static_cast<Eigen::Index>( modes.size() ) );
        probeResiduals->setConstant( std::numeric_limits<double>::quiet_NaN() );
    }
    if( timing )
    {
        *timing = P4PCATiming{};
    }

    std::vector<Eigen::Index> deletedRows;
    std::vector<Eigen::Index> retainedRows;
    if( sampleCount <= predictorCount )
    {
        const auto secondsSince = []( const std::chrono::steady_clock::time_point &start )
        { return std::chrono::duration<double>( std::chrono::steady_clock::now() - start ).count(); };
        const auto fullGramStart = std::chrono::steady_clock::now();
        const matrixT fullGram = ( predictors.matrix() * predictors.matrix().transpose() ).array();
        matrixT fullProbeCrossProduct;
        if( calculateProbe )
        {
            fullProbeCrossProduct = ( probePredictors->matrix() * predictors.matrix().transpose() ).array();
        }
        if( !p4PCAAllFinite( fullGram ) || ( calculateProbe && !p4PCAAllFinite( fullProbeCrossProduct ) ) )
        {
            throw std::runtime_error( "P4PCA normal equations produced nonfinite values" );
        }
        if( timing )
        {
            timing->gramWorkerSeconds = secondsSince( fullGramStart );
        }

        matrixT trainingGram;
        vectorT trainingTarget;
        vectorT heldOutCrossProduct;
        matrixT trainingProbeCrossProduct;
        matrixT eigenvectors;
        matrixT eigenvalues;
        for( Eigen::Index heldOut = 0; heldOut < sampleCount; ++heldOut )
        {
            const std::span<const Eigen::Index> excluded = exclusions.deletedRows( deletedRows, heldOut );
            p4PCARetainedRows( retainedRows, sampleCount, excluded );
            if( retainedRows.empty() )
            {
                output.numericalRank = 0;
                continue;
            }

            const int structuralRank = static_cast<int>( retainedRows.size() );
            output.numericalRank = std::min( output.numericalRank, structuralRank );
            const auto supportedEnd = std::upper_bound( modes.begin(), modes.end(), structuralRank );
            if( supportedEnd == modes.begin() )
            {
                continue;
            }

            const auto extractStart = std::chrono::steady_clock::now();
            trainingGram.resize( structuralRank, structuralRank );
            trainingTarget.resize( structuralRank );
            heldOutCrossProduct.resize( structuralRank );
            if( calculateProbe )
            {
                trainingProbeCrossProduct.resize( probePredictors->rows(), structuralRank );
            }
            for( Eigen::Index trainingColumn = 0; trainingColumn < structuralRank; ++trainingColumn )
            {
                const Eigen::Index sourceColumn = retainedRows[static_cast<std::size_t>( trainingColumn )];
                trainingTarget( trainingColumn ) = target( sourceColumn );
                heldOutCrossProduct( trainingColumn ) = fullGram( heldOut, sourceColumn );
                if( calculateProbe )
                {
                    trainingProbeCrossProduct.col( trainingColumn ) = fullProbeCrossProduct.col( sourceColumn );
                }
                for( Eigen::Index trainingRow = 0; trainingRow < structuralRank; ++trainingRow )
                {
                    const Eigen::Index sourceRow = retainedRows[static_cast<std::size_t>( trainingRow )];
                    trainingGram( trainingRow, trainingColumn ) = fullGram( sourceRow, sourceColumn );
                }
            }
            if( timing )
            {
                timing->gramWorkerSeconds += secondsSince( extractStart );
            }

            const auto eigensolveStart = std::chrono::steady_clock::now();
            const MXLAPACK_INT solverStatus = p4PCAEigenSolveAdapted<CalculationT, EigensolverT>( eigenvectors,
                                                                                                  eigenvalues,
                                                                                                  trainingGram,
                                                                                                  structuralRank,
                                                                                                  workspace );
            if( solverStatus != 0 )
            {
                throw std::runtime_error( "P4PCA eigensolver failed with status " + std::to_string( solverStatus ) );
            }
            if( eigenvectors.rows() != structuralRank || eigenvectors.cols() != structuralRank ||
                eigenvalues.rows() < structuralRank || eigenvalues.cols() < 1 )
            {
                throw std::runtime_error( "P4PCA eigensolver returned invalid dimensions" );
            }
            if( !p4PCAAllFinite( eigenvectors ) || !p4PCAAllFinite( eigenvalues.block( 0, 0, structuralRank, 1 ) ) )
            {
                throw std::runtime_error( "P4PCA eigensolver returned nonfinite results" );
            }
            for( int eigenIndex = 1; eigenIndex < structuralRank; ++eigenIndex )
            {
                if( eigenvalues( eigenIndex ) < eigenvalues( eigenIndex - 1 ) )
                {
                    throw std::runtime_error( "P4PCA eigensolver returned unsorted eigenvalues" );
                }
            }

            const CalculationT largestEigenvalue = eigenvalues( structuralRank - 1 );
            int numericalRank{ 0 };
            if( largestEigenvalue > 0 )
            {
                const long double rankThreshold =
                    static_cast<long double>( rankTolerance ) * static_cast<long double>( largestEigenvalue );
                for( int eigenIndex = 0; eigenIndex < structuralRank; ++eigenIndex )
                {
                    if( static_cast<long double>( eigenvalues( eigenIndex ) ) > rankThreshold )
                    {
                        ++numericalRank;
                    }
                }
            }
            output.numericalRank = std::min( output.numericalRank, numericalRank );
            if( timing )
            {
                timing->eigensolveWorkerSeconds += secondsSince( eigensolveStart );
            }

            const auto numericalEnd = std::upper_bound( modes.begin(), supportedEnd, numericalRank );
            if( numericalEnd == modes.begin() )
            {
                continue;
            }
            const int largestSupportedMode = *( numericalEnd - 1 );
            const auto projectionStart = std::chrono::steady_clock::now();
            CalculationT prediction{ 0 };
            vectorT probePrediction;
            vectorT probeMode;
            if( calculateProbe )
            {
                probePrediction.resize( probePredictors->rows() );
                probePrediction.setZero();
            }
            std::size_t outputIndex{ 0 };
            for( int retainedMode = 1; retainedMode <= largestSupportedMode; ++retainedMode )
            {
                const int eigenIndex = structuralRank - retainedMode;
                const auto eigenvector = eigenvectors.col( eigenIndex );
                const CalculationT heldOutProjection = heldOutCrossProduct.matrix().dot( eigenvector.matrix() );
                const CalculationT targetProjection = eigenvector.matrix().dot( trainingTarget.matrix() );
                prediction += heldOutProjection * ( targetProjection / eigenvalues( eigenIndex ) );
                if( calculateProbe )
                {
                    probeMode = ( trainingProbeCrossProduct.matrix() * eigenvector.matrix() ).array();
                    probePrediction += probeMode * ( targetProjection / eigenvalues( eigenIndex ) );
                }
                if( !mx::math::isFinite( prediction ) || ( calculateProbe && !p4PCAAllFinite( probePrediction ) ) )
                {
                    throw std::runtime_error( "P4PCA held-out prediction produced a nonfinite residual" );
                }
                if( outputIndex < modes.size() && modes[outputIndex] == retainedMode )
                {
                    const CalculationT residual = target( heldOut ) - prediction;
                    if( !mx::math::isFinite( residual ) )
                    {
                        throw std::runtime_error( "P4PCA held-out prediction produced a nonfinite residual" );
                    }
                    output.residuals( heldOut, static_cast<Eigen::Index>( outputIndex ) ) =
                        static_cast<double>( residual );
                    output.sampleValidity( heldOut, static_cast<Eigen::Index>( outputIndex ) ) = 1;
                    output.modeStatus[outputIndex] = P4PCAModeStatus::rankSupported;
                    if( calculateProbe )
                    {
                        const Eigen::Index probeColumn = heldOut * static_cast<Eigen::Index>( modes.size() ) +
                                                         static_cast<Eigen::Index>( outputIndex );
                        const vectorT probeResidual = *probeTarget - probePrediction;
                        if( !p4PCAAllFinite( probeResidual ) )
                        {
                            throw std::runtime_error( "P4PCA held-out probe produced a nonfinite response" );
                        }
                        if constexpr( std::is_same_v<CalculationT, double> )
                        {
                            probeResiduals->col( probeColumn ) = probeResidual;
                        }
                        else
                        {
                            probeResiduals->col( probeColumn ) = probeResidual.template cast<double>();
                        }
                    }
                    ++outputIndex;
                }
            }
            if( timing )
            {
                timing->projectionWorkerSeconds += secondsSince( projectionStart );
            }
        }
        return;
    }

    matrixT trainingPredictors;
    vectorT trainingTarget;
    P4PCA::matrixT exportedCoefficients;
    matrixT calculationCoefficients;
    for( Eigen::Index heldOut = 0; heldOut < sampleCount; ++heldOut )
    {
        const std::span<const Eigen::Index> excluded = exclusions.deletedRows( deletedRows, heldOut );
        p4PCARetainedRows( retainedRows, sampleCount, excluded );
        if( retainedRows.empty() )
        {
            output.numericalRank = 0;
            continue;
        }

        const int structuralRank = static_cast<int>(
            std::min<std::size_t>( retainedRows.size(), static_cast<std::size_t>( predictorCount ) ) );
        output.numericalRank = std::min( output.numericalRank, structuralRank );
        const auto supportedEnd = std::upper_bound( modes.begin(), modes.end(), structuralRank );
        if( supportedEnd == modes.begin() )
        {
            continue;
        }
        const std::vector<int> supportedModes( modes.begin(), supportedEnd );
        trainingPredictors.resize( static_cast<Eigen::Index>( retainedRows.size() ), predictorCount );
        trainingTarget.resize( static_cast<Eigen::Index>( retainedRows.size() ) );
        for( std::size_t training = 0; training < retainedRows.size(); ++training )
        {
            const Eigen::Index source = retainedRows[training];
            trainingPredictors.row( static_cast<Eigen::Index>( training ) ) = predictors.row( source );
            trainingTarget( static_cast<Eigen::Index>( training ) ) = target( source );
        }

        P4PCAResult fit;
        P4PCATiming fitTiming;
        P4PCAKernel<CalculationT, EigensolverT>::calculate( fit,
                                                            trainingPredictors,
                                                            trainingTarget,
                                                            supportedModes,
                                                            rankTolerance,
                                                            workspace,
                                                            &fitTiming,
                                                            &exportedCoefficients );
        const matrixT *coefficients = nullptr;
        if constexpr( std::is_same_v<CalculationT, double> )
        {
            coefficients = &exportedCoefficients;
        }
        else
        {
            calculationCoefficients = exportedCoefficients.template cast<CalculationT>();
            coefficients = &calculationCoefficients;
        }
        output.numericalRank = std::min( output.numericalRank, fit.numericalRank );
        if( timing )
        {
            timing->gramWorkerSeconds += fitTiming.gramWorkerSeconds;
            timing->eigensolveWorkerSeconds += fitTiming.eigensolveWorkerSeconds;
            timing->projectionWorkerSeconds += fitTiming.projectionWorkerSeconds;
        }
        const auto heldOutProjectionBegin = std::chrono::steady_clock::now();
        for( std::size_t mode = 0; mode < supportedModes.size(); ++mode )
        {
            if( fit.modeStatus[mode] == P4PCAModeStatus::rankInsufficient )
            {
                continue;
            }
            const CalculationT prediction = predictors.row( heldOut ).matrix().dot(
                coefficients->col( static_cast<Eigen::Index>( mode ) ).matrix() );
            const CalculationT residual = target( heldOut ) - prediction;
            if( !mx::math::isFinite( residual ) )
            {
                throw std::runtime_error( "P4PCA held-out prediction produced a nonfinite residual" );
            }
            output.residuals( heldOut, static_cast<Eigen::Index>( mode ) ) = static_cast<double>( residual );
            output.sampleValidity( heldOut, static_cast<Eigen::Index>( mode ) ) = 1;
            output.modeStatus[mode] = P4PCAModeStatus::rankSupported;
            if( calculateProbe )
            {
                const Eigen::Index probeColumn =
                    heldOut * static_cast<Eigen::Index>( modes.size() ) + static_cast<Eigen::Index>( mode );
                const vectorT probeResidual =
                    *probeTarget -
                    ( probePredictors->matrix() * coefficients->col( static_cast<Eigen::Index>( mode ) ).matrix() )
                        .array();
                if( !p4PCAAllFinite( probeResidual ) )
                {
                    throw std::runtime_error( "P4PCA held-out probe produced a nonfinite response" );
                }
                if constexpr( std::is_same_v<CalculationT, double> )
                {
                    probeResiduals->col( probeColumn ) = probeResidual;
                }
                else
                {
                    probeResiduals->col( probeColumn ) = probeResidual.template cast<double>();
                }
            }
        }
        if( timing )
        {
            timing->projectionWorkerSeconds +=
                std::chrono::duration<double>( std::chrono::steady_clock::now() - heldOutProjectionBegin ).count();
        }
    }
}

/// Dispatch one explicit held-out request through the production mixed-precision scalar policy.
void p4PCACalculateHeldOutMixedImpl(
    P4PCAResult &output,                    /**< [out] held-out science result */
    P4PCA::matrixT *probeResiduals,         /**< [out] optional FP64 frozen-probe responses */
    const P4PCA::matrixT &predictors,       /**< [in] common FP64 predictor ingress */
    const P4PCA::vectorT &target,           /**< [in] common FP64 target ingress */
    const P4PCA::matrixT *probePredictors,  /**< [in] optional common FP64 probe predictors */
    const P4PCA::vectorT *probeTarget,      /**< [in] optional common FP64 direct probe response */
    const P4TargetExclusions &exclusions,   /**< [in] target-specific deleted rows */
    const std::vector<int> &modes,          /**< [in] requested retained-mode counts */
    double rankTolerance,                   /**< [in] relative numerical-rank threshold */
    detail::P4PCAMixedWorkspace &workspace, /**< [in,out] reusable mixed-policy workspace */
    P4PCATiming *timing /**< [out] optional aggregate timing */ )
{
    workspace.floatPredictors = predictors.cast<float>();
    workspace.floatTarget = target.cast<float>();
    if( probePredictors && probeTarget )
    {
        workspace.floatProbePredictors = probePredictors->cast<float>();
        workspace.floatProbeTarget = probeTarget->cast<float>();
    }
    p4PCACalculateHeldOut<float, double>( output,
                                          probeResiduals,
                                          workspace.floatPredictors,
                                          workspace.floatTarget,
                                          probePredictors ? &workspace.floatProbePredictors : nullptr,
                                          probeTarget ? &workspace.floatProbeTarget : nullptr,
                                          exclusions,
                                          modes,
                                          rankTolerance,
                                          workspace.doubleEigensolver,
                                          timing );
}

#ifdef HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION
/// Dispatch one explicit held-out request through the selected experimental scalar policy.
void p4PCACalculateHeldOutExperimentalImpl(
    P4PCAResult &output,                           /**< [out] held-out science result */
    P4PCA::matrixT *probeResiduals,                /**< [out] optional FP64 frozen-probe responses */
    const P4PCA::matrixT &predictors,              /**< [in] common FP64 predictor ingress */
    const P4PCA::vectorT &target,                  /**< [in] common FP64 target ingress */
    const P4PCA::matrixT *probePredictors,         /**< [in] optional common FP64 probe predictors */
    const P4PCA::vectorT *probeTarget,             /**< [in] optional common FP64 direct probe response */
    const P4TargetExclusions &exclusions,          /**< [in] target-specific deleted rows */
    const std::vector<int> &modes,                 /**< [in] requested retained-mode counts */
    double rankTolerance,                          /**< [in] relative numerical-rank threshold */
    detail::P4PCAPrecisionPolicy precisionPolicy,  /**< [in] calculation/eigensolver scalar combination */
    detail::P4PCAExperimentalWorkspace &workspace, /**< [in,out] reusable policy workspace */
    P4PCATiming *timing /**< [out] optional aggregate timing */ )
{
    switch( precisionPolicy )
    {
    case detail::P4PCAPrecisionPolicy::doubleDouble:
        p4PCACalculateHeldOut<double, double>( output,
                                               probeResiduals,
                                               predictors,
                                               target,
                                               probePredictors,
                                               probeTarget,
                                               exclusions,
                                               modes,
                                               rankTolerance,
                                               workspace.doubleEigensolver,
                                               timing );
        return;
    case detail::P4PCAPrecisionPolicy::floatDouble:
        p4PCACalculateHeldOutMixedImpl( output,
                                        probeResiduals,
                                        predictors,
                                        target,
                                        probePredictors,
                                        probeTarget,
                                        exclusions,
                                        modes,
                                        rankTolerance,
                                        workspace,
                                        timing );
        return;
    case detail::P4PCAPrecisionPolicy::floatFloat:
        workspace.floatPredictors = predictors.cast<float>();
        workspace.floatTarget = target.cast<float>();
        if( probePredictors && probeTarget )
        {
            workspace.floatProbePredictors = probePredictors->cast<float>();
            workspace.floatProbeTarget = probeTarget->cast<float>();
        }
        p4PCACalculateHeldOut<float, float>( output,
                                             probeResiduals,
                                             workspace.floatPredictors,
                                             workspace.floatTarget,
                                             probePredictors ? &workspace.floatProbePredictors : nullptr,
                                             probeTarget ? &workspace.floatProbeTarget : nullptr,
                                             exclusions,
                                             modes,
                                             rankTolerance,
                                             workspace.floatEigensolver,
                                             timing );
        return;
    }

    throw std::invalid_argument( "P4PCA experimental precision policy is invalid" );
}
#endif

} // namespace

void P4PCA::calculateHeldOut( P4PCAResult &output,
                              const matrixT &predictors,
                              const vectorT &target,
                              const P4TargetExclusions &exclusions,
                              const std::vector<int> &modes,
                              double rankTolerance,
                              workspaceT &workspace,
                              P4PCATiming *timing )
{
    p4PCACalculateHeldOut<double, double>( output,
                                           nullptr,
                                           predictors,
                                           target,
                                           nullptr,
                                           nullptr,
                                           exclusions,
                                           modes,
                                           rankTolerance,
                                           workspace,
                                           timing );
}

void P4PCA::calculateHeldOutProbe( P4PCAResult &output,
                                   matrixT &probeResiduals,
                                   const matrixT &predictors,
                                   const vectorT &target,
                                   const matrixT &probePredictors,
                                   const vectorT &probeTarget,
                                   const P4TargetExclusions &exclusions,
                                   const std::vector<int> &modes,
                                   double rankTolerance,
                                   workspaceT &workspace,
                                   P4PCATiming *timing )
{
    p4PCAValidateProbeOutputAliasing( probeResiduals, predictors, probePredictors );
    p4PCACalculateHeldOut<double, double>( output,
                                           &probeResiduals,
                                           predictors,
                                           target,
                                           &probePredictors,
                                           &probeTarget,
                                           exclusions,
                                           modes,
                                           rankTolerance,
                                           workspace,
                                           timing );
}

namespace detail
{

/// \cond P4PCA_test_detail
void p4PCACalculateHeldOutMixed( P4PCAResult &output,
                                 const P4PCA::matrixT &predictors,
                                 const P4PCA::vectorT &target,
                                 const P4TargetExclusions &exclusions,
                                 const std::vector<int> &modes,
                                 double rankTolerance,
                                 P4PCAMixedWorkspace &workspace,
                                 P4PCATiming *timing )
{
    p4PCACalculateHeldOutMixedImpl( output,
                                    nullptr,
                                    predictors,
                                    target,
                                    nullptr,
                                    nullptr,
                                    exclusions,
                                    modes,
                                    rankTolerance,
                                    workspace,
                                    timing );
}

void p4PCACalculateHeldOutProbeMixed( P4PCAResult &output,
                                      P4PCA::matrixT &probeResiduals,
                                      const P4PCA::matrixT &predictors,
                                      const P4PCA::vectorT &target,
                                      const P4PCA::matrixT &probePredictors,
                                      const P4PCA::vectorT &probeTarget,
                                      const P4TargetExclusions &exclusions,
                                      const std::vector<int> &modes,
                                      double rankTolerance,
                                      P4PCAMixedWorkspace &workspace,
                                      P4PCATiming *timing )
{
    p4PCAValidateProbeOutputAliasing( probeResiduals, predictors, probePredictors );
    p4PCACalculateHeldOutMixedImpl( output,
                                    &probeResiduals,
                                    predictors,
                                    target,
                                    &probePredictors,
                                    &probeTarget,
                                    exclusions,
                                    modes,
                                    rankTolerance,
                                    workspace,
                                    timing );
}
/// \endcond

} // namespace detail

#ifdef HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION
namespace detail
{

/// \cond P4PCA_test_detail
void p4PCACalculateHeldOutExperimental( P4PCAResult &output,
                                        const P4PCA::matrixT &predictors,
                                        const P4PCA::vectorT &target,
                                        const P4TargetExclusions &exclusions,
                                        const std::vector<int> &modes,
                                        double rankTolerance,
                                        P4PCAPrecisionPolicy precisionPolicy,
                                        P4PCAExperimentalWorkspace &workspace,
                                        P4PCATiming *timing )
{
    p4PCACalculateHeldOutExperimentalImpl( output,
                                           nullptr,
                                           predictors,
                                           target,
                                           nullptr,
                                           nullptr,
                                           exclusions,
                                           modes,
                                           rankTolerance,
                                           precisionPolicy,
                                           workspace,
                                           timing );
}

void p4PCACalculateHeldOutProbeExperimental( P4PCAResult &output,
                                             P4PCA::matrixT &probeResiduals,
                                             const P4PCA::matrixT &predictors,
                                             const P4PCA::vectorT &target,
                                             const P4PCA::matrixT &probePredictors,
                                             const P4PCA::vectorT &probeTarget,
                                             const P4TargetExclusions &exclusions,
                                             const std::vector<int> &modes,
                                             double rankTolerance,
                                             P4PCAPrecisionPolicy precisionPolicy,
                                             P4PCAExperimentalWorkspace &workspace,
                                             P4PCATiming *timing )
{
    p4PCAValidateProbeOutputAliasing( probeResiduals, predictors, probePredictors );
    p4PCACalculateHeldOutExperimentalImpl( output,
                                           &probeResiduals,
                                           predictors,
                                           target,
                                           &probePredictors,
                                           &probeTarget,
                                           exclusions,
                                           modes,
                                           rankTolerance,
                                           precisionPolicy,
                                           workspace,
                                           timing );
}
/// \endcond

} // namespace detail
#endif

void P4PCA::calculateHeldOutDowndated( P4PCAResult &output,
                                       const matrixT &predictors,
                                       const vectorT &target,
                                       const P4TargetExclusions &exclusions,
                                       const std::vector<int> &modes,
                                       double rankTolerance,
                                       workspaceT &eigensolverWorkspace,
                                       P4PCADowndateWorkspace &downdateWorkspace,
                                       P4PCATiming *timing )
{
    calculateHeldOutDowndated( output,
                               predictors,
                               target,
                               exclusions,
                               modes,
                               rankTolerance,
                               mx::math::svdDeletionBackend::leadingCovariance,
                               eigensolverWorkspace,
                               downdateWorkspace,
                               timing );
}

void P4PCA::calculateHeldOutDowndated( P4PCAResult &output,
                                       const matrixT &predictors,
                                       const vectorT &target,
                                       const P4TargetExclusions &exclusions,
                                       const std::vector<int> &modes,
                                       double rankTolerance,
                                       mx::math::svdDeletionBackend deletionBackend,
                                       workspaceT &eigensolverWorkspace,
                                       P4PCADowndateWorkspace &downdateWorkspace,
                                       P4PCATiming *timing )
{
    calculateHeldOutDowndatedImpl( output,
                                   nullptr,
                                   predictors,
                                   target,
                                   nullptr,
                                   nullptr,
                                   exclusions,
                                   modes,
                                   rankTolerance,
                                   deletionBackend,
                                   eigensolverWorkspace,
                                   downdateWorkspace,
                                   timing );
}

void P4PCA::calculateHeldOutDowndatedImpl( P4PCAResult &output,
                                           matrixT *probeResiduals,
                                           const matrixT &predictors,
                                           const vectorT &target,
                                           const matrixT *probePredictors,
                                           const vectorT *probeTarget,
                                           const P4TargetExclusions &exclusions,
                                           const std::vector<int> &modes,
                                           double rankTolerance,
                                           mx::math::svdDeletionBackend deletionBackend,
                                           workspaceT &eigensolverWorkspace,
                                           P4PCADowndateWorkspace &downdateWorkspace,
                                           P4PCATiming *timing )
{
    const Eigen::Index sampleCount = predictors.rows();
    const Eigen::Index predictorCount = predictors.cols();
    const Eigen::Index structuralBaseRank = std::min( sampleCount, predictorCount );
    if( structuralBaseRank > static_cast<Eigen::Index>( std::numeric_limits<int>::max() ) )
    {
        throw std::invalid_argument( "P4PCA downdate base rank exceeds the supported integer range" );
    }
    p4PCAValidateInputs( predictors, target, modes, rankTolerance, static_cast<int>( structuralBaseRank ) );
    if( exclusions.sampleCount() != sampleCount || exclusions.targetCount() != sampleCount )
    {
        throw std::invalid_argument( "P4PCA downdate exclusions must match the predictor row count" );
    }
    if( deletionBackend != mx::math::svdDeletionBackend::leadingCovariance &&
        deletionBackend != mx::math::svdDeletionBackend::rankOneSecular )
    {
        throw std::invalid_argument( "P4PCA supports only leadingCovariance and rankOneSecular deletion backends" );
    }
    const bool calculateProbe = probeResiduals || probePredictors || probeTarget;
    if( calculateProbe &&
        ( !probeResiduals || !probePredictors || !probeTarget || probePredictors->rows() <= 0 ||
          probePredictors->cols() != predictorCount || probeTarget->rows() != probePredictors->rows() ||
          !p4PCAAllFinite( *probePredictors ) || !p4PCAAllFinite( *probeTarget ) ) )
    {
        throw std::invalid_argument( "P4PCA downdated probe inputs must have consistent finite dimensions" );
    }

    output.residuals.resize( sampleCount, static_cast<Eigen::Index>( modes.size() ) );
    output.residuals.setConstant( std::numeric_limits<double>::quiet_NaN() );
    output.sampleValidity.resize( sampleCount, static_cast<Eigen::Index>( modes.size() ) );
    output.sampleValidity.setZero();
    output.modeStatus.assign( modes.size(), P4PCAModeStatus::rankInsufficient );
    output.numericalRank = std::numeric_limits<int>::max();
    output.baseRank = 0;
    output.numericalRankCapped = false;
    output.downdateClampCount = 0;
    output.explicitFallbackCount = 0;
    output.explicitFallbackReason = P4PCAFallbackReason::none;
    output.factorOrthogonalityDefect = 0;
    output.factorOrthogonalityTolerance = 0;
    if( calculateProbe )
    {
        probeResiduals->resize( probePredictors->rows(), sampleCount * static_cast<Eigen::Index>( modes.size() ) );
        probeResiduals->setConstant( std::numeric_limits<double>::quiet_NaN() );
    }
    if( timing )
    {
        *timing = P4PCATiming{};
    }

    Eigen::Index maximumActiveDeleted{ 0 };
    bool hasRetainedTarget{ false };
    for( Eigen::Index targetIndex = 0; targetIndex < sampleCount; ++targetIndex )
    {
        const Eigen::Index deletedCount = exclusions.deletedCount( targetIndex );
        if( deletedCount < sampleCount )
        {
            hasRetainedTarget = true;
            maximumActiveDeleted = std::max( maximumActiveDeleted, deletedCount );
        }
    }
    if( deletionBackend == mx::math::svdDeletionBackend::rankOneSecular && hasRetainedTarget )
    {
        for( Eigen::Index targetIndex = 0; targetIndex < sampleCount; ++targetIndex )
        {
            const Eigen::Index deletedCount = exclusions.deletedCount( targetIndex );
            if( deletedCount < sampleCount && deletedCount != 1 )
            {
                throw std::invalid_argument(
                    "P4PCA rankOneSecular deletion requires exactly one deleted row for every retained target" );
            }
        }
    }
    if( !hasRetainedTarget )
    {
        output.numericalRank = 0;
        return;
    }

    const auto secondsSince = []( const std::chrono::steady_clock::time_point &start )
    { return std::chrono::duration<double>( std::chrono::steady_clock::now() - start ).count(); };
    const bool temporalGram = sampleCount <= predictorCount;
    const auto gramStart = std::chrono::steady_clock::now();
    if( temporalGram )
    {
        downdateWorkspace.m_baseGram = ( predictors.matrix() * predictors.matrix().transpose() ).array();
    }
    else
    {
        downdateWorkspace.m_baseGram = ( predictors.matrix().transpose() * predictors.matrix() ).array();
    }
    if( !p4PCAAllFinite( downdateWorkspace.m_baseGram ) )
    {
        throw std::runtime_error( "P4PCA downdate base Gram matrix contains nonfinite values" );
    }
    if( timing )
    {
        timing->gramWorkerSeconds = secondsSince( gramStart );
    }

    const auto baseSolveStart = std::chrono::steady_clock::now();
    const int gramDimension = static_cast<int>( downdateWorkspace.m_baseGram.rows() );
    const MXLAPACK_INT baseStatus = p4PCAEigenSolve( downdateWorkspace.m_baseEigenvectors,
                                                     downdateWorkspace.m_baseEigenvalues,
                                                     downdateWorkspace.m_baseGram,
                                                     gramDimension,
                                                     eigensolverWorkspace );
    if( baseStatus != 0 )
    {
        throw std::runtime_error( "P4PCA downdate base eigensolver failed with status " +
                                  std::to_string( baseStatus ) );
    }
    if( downdateWorkspace.m_baseEigenvectors.rows() != gramDimension ||
        downdateWorkspace.m_baseEigenvectors.cols() != gramDimension ||
        downdateWorkspace.m_baseEigenvalues.rows() < gramDimension || downdateWorkspace.m_baseEigenvalues.cols() < 1 ||
        !p4PCAAllFinite( downdateWorkspace.m_baseEigenvectors ) ||
        !p4PCAAllFinite( downdateWorkspace.m_baseEigenvalues.block( 0, 0, gramDimension, 1 ) ) )
    {
        throw std::runtime_error( "P4PCA downdate base eigensolver returned an invalid singular system" );
    }
    for( int index = 1; index < gramDimension; ++index )
    {
        if( downdateWorkspace.m_baseEigenvalues( index ) < downdateWorkspace.m_baseEigenvalues( index - 1 ) )
        {
            throw std::runtime_error( "P4PCA downdate base eigensolver returned unsorted eigenvalues" );
        }
    }

    const double spectrumScale = std::max( std::abs( downdateWorkspace.m_baseEigenvalues( 0 ) ),
                                           std::abs( downdateWorkspace.m_baseEigenvalues( gramDimension - 1 ) ) );
    const double baseSafetyTolerance = 64.0 * std::numeric_limits<double>::epsilon() *
                                       static_cast<double>( std::max( 1, gramDimension ) ) * spectrumScale;
    if( downdateWorkspace.m_baseEigenvalues( 0 ) < -baseSafetyTolerance )
    {
        throw std::runtime_error( "P4PCA downdate base Gram matrix is materially non-positive-semidefinite" );
    }

    Eigen::Index baseRank{ 0 };
    for( Eigen::Index source = 0; source < structuralBaseRank; ++source )
    {
        const double eigenvalue = downdateWorkspace.m_baseEigenvalues( source );
        if( eigenvalue < -baseSafetyTolerance )
        {
            throw std::runtime_error( "P4PCA downdate base Gram matrix has a negative eigenvalue" );
        }
        if( eigenvalue > baseSafetyTolerance )
        {
            ++baseRank;
        }
    }

    if( temporalGram )
    {
        downdateWorkspace.m_temporalFactor.resize( sampleCount, baseRank );
        downdateWorkspace.m_singularValues.resize( baseRank );
        for( Eigen::Index mode = 0; mode < baseRank; ++mode )
        {
            const Eigen::Index source = structuralBaseRank - 1 - mode;
            const double eigenvalue = downdateWorkspace.m_baseEigenvalues( source );
            downdateWorkspace.m_singularValues( mode ) = std::sqrt( eigenvalue );
            downdateWorkspace.m_temporalFactor.col( mode ) = downdateWorkspace.m_baseEigenvectors.col( source );
        }
    }
    else
    {
        downdateWorkspace.m_temporalFactor.resize( sampleCount, baseRank );
        downdateWorkspace.m_singularValues.resize( baseRank );
        for( Eigen::Index mode = 0; mode < baseRank; ++mode )
        {
            const Eigen::Index source = structuralBaseRank - 1 - mode;
            const double singularValue = std::sqrt( downdateWorkspace.m_baseEigenvalues( source ) );
            downdateWorkspace.m_singularValues( mode ) = singularValue;
            downdateWorkspace.m_temporalFactor.col( mode ).matrix().noalias() =
                predictors.matrix() * downdateWorkspace.m_baseEigenvectors.col( source ).matrix();
            downdateWorkspace.m_temporalFactor.col( mode ) /= singularValue;
        }
    }
    output.baseRank = static_cast<int>( baseRank );
    if( baseRank == 0 )
    {
        output.numericalRank = 0;
        if( timing )
        {
            timing->eigensolveWorkerSeconds = secondsSince( baseSolveStart );
            timing->baseFactorWorkerSeconds = timing->eigensolveWorkerSeconds;
        }
        return;
    }

    const double factorTolerance = p4PCAFactorValidationTolerance( downdateWorkspace.m_temporalFactor );
    const mx::math::svdDeletionStatus factorStatus =
        mx::math::validateSvdDeletionFactor( downdateWorkspace.m_temporalFactor, factorTolerance );
    if( !mx::math::svdDeletionSucceeded( factorStatus ) )
    {
        if( factorStatus != mx::math::svdDeletionStatus::factorNotOrthonormal )
        {
            throw std::runtime_error(
                p4PCAFactorValidationMessage( factorStatus, downdateWorkspace.m_temporalFactor, factorTolerance ) );
        }
        output.factorOrthogonalityDefect = p4PCAFactorOrthogonalityDefect( downdateWorkspace.m_temporalFactor );
        output.factorOrthogonalityTolerance = factorTolerance;
        if( timing )
        {
            timing->eigensolveWorkerSeconds = secondsSince( baseSolveStart );
            timing->baseFactorWorkerSeconds = timing->eigensolveWorkerSeconds;
        }
        p4PCAExplicitFallback( output,
                               predictors,
                               target,
                               exclusions,
                               modes,
                               rankTolerance,
                               eigensolverWorkspace,
                               P4PCAFallbackReason::factorValidation,
                               timing,
                               probeResiduals,
                               probePredictors,
                               probeTarget );
        return;
    }

    if( calculateProbe )
    {
        if( temporalGram )
        {
            downdateWorkspace.m_probeBaseFactor = ( probePredictors->matrix() * predictors.matrix().transpose() *
                                                    downdateWorkspace.m_temporalFactor.matrix() )
                                                      .array();
            for( Eigen::Index mode = 0; mode < baseRank; ++mode )
            {
                downdateWorkspace.m_probeBaseFactor.col( mode ) /= downdateWorkspace.m_singularValues( mode );
            }
        }
        else
        {
            downdateWorkspace.m_probeBaseFactor.resize( probePredictors->rows(), baseRank );
            for( Eigen::Index mode = 0; mode < baseRank; ++mode )
            {
                const Eigen::Index source = structuralBaseRank - 1 - mode;
                downdateWorkspace.m_probeBaseFactor.col( mode ).matrix().noalias() =
                    probePredictors->matrix() * downdateWorkspace.m_baseEigenvectors.col( source ).matrix();
            }
        }
        if( !p4PCAAllFinite( downdateWorkspace.m_probeBaseFactor ) )
        {
            throw std::runtime_error( "P4PCA downdate probe base factor contains nonfinite values" );
        }
    }

    const Eigen::Index maximumOutputRank = std::min<Eigen::Index>( modes.back(), baseRank );
    const mx::math::svdDeletionStatus resultPrepareStatus =
        downdateWorkspace.m_deletionResult.prepare( baseRank, maximumOutputRank );
    const mx::math::svdDeletionStatus workspacePrepareStatus =
        downdateWorkspace.m_deletionWorkspace.prepare( baseRank, maximumActiveDeleted, deletionBackend );
    if( !mx::math::svdDeletionSucceeded( resultPrepareStatus ) ||
        !mx::math::svdDeletionSucceeded( workspacePrepareStatus ) )
    {
        const mx::math::svdDeletionStatus status =
            !mx::math::svdDeletionSucceeded( resultPrepareStatus ) ? resultPrepareStatus : workspacePrepareStatus;
        throw std::runtime_error( "P4PCA downdate workspace preparation failed with status " +
                                  std::string( mx::math::svdDeletionStatusName( status ) ) );
    }

    downdateWorkspace.m_baseTargetProjection =
        ( downdateWorkspace.m_temporalFactor.matrix().transpose() * target.matrix() ).array();
    if( !p4PCAAllFinite( downdateWorkspace.m_baseTargetProjection ) )
    {
        throw std::runtime_error( "P4PCA downdate base target projection contains nonfinite values" );
    }
    if( timing )
    {
        timing->eigensolveWorkerSeconds = secondsSince( baseSolveStart );
        timing->baseFactorWorkerSeconds = timing->eigensolveWorkerSeconds;
    }

    bool requiresExplicitOracle{ false };
    for( Eigen::Index heldOut = 0; heldOut < sampleCount; ++heldOut )
    {
        const std::span<const Eigen::Index> deleted =
            exclusions.deletedRows( downdateWorkspace.m_deletedRows, heldOut );
        const Eigen::Index retainedCount = sampleCount - static_cast<Eigen::Index>( deleted.size() );
        if( retainedCount == 0 )
        {
            output.numericalRank = 0;
            continue;
        }

        const int structuralRank = static_cast<int>( std::min<Eigen::Index>( baseRank, retainedCount ) );
        output.numericalRank = std::min( output.numericalRank, structuralRank );
        const auto structuralEnd = std::upper_bound( modes.begin(), modes.end(), structuralRank );
        if( structuralEnd == modes.begin() )
        {
            continue;
        }

        const Eigen::Index requestedOutputRank = *( structuralEnd - 1 );
        const auto deletionStart = std::chrono::steady_clock::now();
        const mx::math::svdDeletionStatus deletionStatus =
            mx::math::svdRemoveRows<double>( downdateWorkspace.m_deletionResult,
                                             downdateWorkspace.m_singularValues,
                                             downdateWorkspace.m_temporalFactor,
                                             deleted,
                                             requestedOutputRank,
                                             downdateWorkspace.m_deletionWorkspace,
                                             deletionBackend );
        const double deletionSeconds = secondsSince( deletionStart );
        if( timing )
        {
            timing->eigensolveWorkerSeconds += deletionSeconds;
            timing->deletionWorkerSeconds += deletionSeconds;
        }
        if( !mx::math::svdDeletionSucceeded( deletionStatus ) )
        {
            if( deletionBackend == mx::math::svdDeletionBackend::rankOneSecular &&
                p4PCARecoverableStructuredFailure( deletionStatus, downdateWorkspace.m_deletionResult.lapackInfo() ) )
            {
                p4PCAExplicitFallback( output,
                                       predictors,
                                       target,
                                       exclusions,
                                       modes,
                                       rankTolerance,
                                       eigensolverWorkspace,
                                       P4PCAFallbackReason::deletionSolver,
                                       timing,
                                       probeResiduals,
                                       probePredictors,
                                       probeTarget );
                return;
            }
            throw std::runtime_error( "P4PCA row deletion failed for target " + std::to_string( heldOut ) +
                                      " with status " + mx::math::svdDeletionStatusName( deletionStatus ) +
                                      " and LAPACK status " +
                                      std::to_string( downdateWorkspace.m_deletionResult.lapackInfo() ) );
        }
        output.downdateClampCount +=
            static_cast<std::size_t>( downdateWorkspace.m_deletionResult.clampedEigenvalues() );
        const auto squaredSingularValues = downdateWorkspace.m_deletionResult.squaredSingularValues();
        int numericalRank{ 0 };
        const double largestEigenvalue = squaredSingularValues( 0 );
        const double rankThreshold = rankTolerance * largestEigenvalue;
        const double ambiguityScale =
            std::max( std::abs( downdateWorkspace.m_baseEigenvalues( structuralBaseRank - 1 ) ),
                      std::abs( largestEigenvalue ) );
        const double ambiguityTolerance = 256.0 * std::numeric_limits<double>::epsilon() *
                                          static_cast<double>( std::max<Eigen::Index>( 1, baseRank ) ) * ambiguityScale;
        if( largestEigenvalue > 0 )
        {
            for( int mode = 0; mode < structuralRank; ++mode )
            {
                if( std::abs( squaredSingularValues( mode ) - rankThreshold ) <= ambiguityTolerance )
                {
                    requiresExplicitOracle = true;
                    break;
                }
                if( squaredSingularValues( mode ) > rankThreshold )
                {
                    ++numericalRank;
                }
            }
        }
        if( requiresExplicitOracle )
        {
            break;
        }
        output.numericalRank = std::min( output.numericalRank, numericalRank );
        const auto numericalEnd = std::upper_bound( modes.begin(), structuralEnd, numericalRank );
        if( numericalEnd == modes.begin() )
        {
            continue;
        }

        const auto projectionStart = std::chrono::steady_clock::now();
        downdateWorkspace.m_retainedTargetProjection = downdateWorkspace.m_baseTargetProjection;
        for( const Eigen::Index row : deleted )
        {
            downdateWorkspace.m_retainedTargetProjection -=
                downdateWorkspace.m_temporalFactor.row( row ).transpose() * target( row );
        }
        downdateWorkspace.m_scaledTargetProjection =
            downdateWorkspace.m_singularValues * downdateWorkspace.m_retainedTargetProjection;
        downdateWorkspace.m_scaledHeldOutRow =
            downdateWorkspace.m_singularValues * downdateWorkspace.m_temporalFactor.row( heldOut ).transpose();
        const auto rotation = downdateWorkspace.m_deletionResult.rotation();
        double prediction{ 0 };
        if( calculateProbe )
        {
            downdateWorkspace.m_probePrediction.resize( probePredictors->rows() );
            downdateWorkspace.m_probePrediction.setZero();
        }
        std::size_t outputIndex{ 0 };
        const int largestSupportedMode = *( numericalEnd - 1 );
        for( int retainedMode = 1; retainedMode <= largestSupportedMode; ++retainedMode )
        {
            const Eigen::Index mode = retainedMode - 1;
            const double targetCoordinate =
                rotation.col( mode ).matrix().dot( downdateWorkspace.m_scaledTargetProjection.matrix() );
            const double heldOutCoordinate =
                rotation.col( mode ).matrix().dot( downdateWorkspace.m_scaledHeldOutRow.matrix() );
            prediction += heldOutCoordinate * targetCoordinate / squaredSingularValues( mode );
            if( calculateProbe )
            {
                downdateWorkspace.m_probeMode =
                    ( downdateWorkspace.m_probeBaseFactor.matrix() * rotation.col( mode ).matrix() ).array();
                downdateWorkspace.m_probePrediction +=
                    downdateWorkspace.m_probeMode * ( targetCoordinate / squaredSingularValues( mode ) );
            }
            if( !mx::math::isFinite( prediction ) ||
                ( calculateProbe && !p4PCAAllFinite( downdateWorkspace.m_probePrediction ) ) )
            {
                throw std::runtime_error( "P4PCA factor-space prediction produced a nonfinite value" );
            }
            if( outputIndex < modes.size() && modes[outputIndex] == retainedMode )
            {
                const double residual = target( heldOut ) - prediction;
                if( !mx::math::isFinite( residual ) )
                {
                    throw std::runtime_error( "P4PCA factor-space residual produced a nonfinite value" );
                }
                output.residuals( heldOut, static_cast<Eigen::Index>( outputIndex ) ) = residual;
                output.sampleValidity( heldOut, static_cast<Eigen::Index>( outputIndex ) ) = 1;
                output.modeStatus[outputIndex] = P4PCAModeStatus::rankSupported;
                if( calculateProbe )
                {
                    const Eigen::Index probeColumn =
                        heldOut * static_cast<Eigen::Index>( modes.size() ) + static_cast<Eigen::Index>( outputIndex );
                    probeResiduals->col( probeColumn ) = *probeTarget - downdateWorkspace.m_probePrediction;
                    if( !p4PCAAllFinite( probeResiduals->col( probeColumn ) ) )
                    {
                        throw std::runtime_error( "P4PCA downdated probe produced a nonfinite response" );
                    }
                }
                ++outputIndex;
            }
        }
        if( timing )
        {
            timing->projectionWorkerSeconds += secondsSince( projectionStart );
        }
    }

    if( requiresExplicitOracle )
    {
        p4PCAExplicitFallback( output,
                               predictors,
                               target,
                               exclusions,
                               modes,
                               rankTolerance,
                               eigensolverWorkspace,
                               P4PCAFallbackReason::rankBoundary,
                               timing,
                               probeResiduals,
                               probePredictors,
                               probeTarget );
    }
}

void P4PCA::calculateHeldOutProbeDowndated( P4PCAResult &output,
                                            matrixT &probeResiduals,
                                            const matrixT &predictors,
                                            const vectorT &target,
                                            const matrixT &probePredictors,
                                            const vectorT &probeTarget,
                                            const P4TargetExclusions &exclusions,
                                            const std::vector<int> &modes,
                                            double rankTolerance,
                                            mx::math::svdDeletionBackend deletionBackend,
                                            workspaceT &eigensolverWorkspace,
                                            P4PCADowndateWorkspace &downdateWorkspace,
                                            P4PCATiming *timing )
{
    calculateHeldOutDowndatedImpl( output,
                                   &probeResiduals,
                                   predictors,
                                   target,
                                   &probePredictors,
                                   &probeTarget,
                                   exclusions,
                                   modes,
                                   rankTolerance,
                                   deletionBackend,
                                   eigensolverWorkspace,
                                   downdateWorkspace,
                                   timing );
}

void P4PCA::calculateCentered( P4PCAResult &output,
                               const matrixT &predictors,
                               const vectorT &target,
                               const std::vector<int> &modes,
                               double rankTolerance,
                               workspaceT &workspace,
                               P4PCATiming *timing,
                               matrixT *coefficients )
{
    P4PCAKernel<double, double>::calculateCentered( output,
                                                    predictors,
                                                    target,
                                                    modes,
                                                    rankTolerance,
                                                    workspace,
                                                    timing,
                                                    coefficients );
}

void P4PCA::calculateCenteredInPlace( P4PCAResult &output,
                                      matrixT &predictors,
                                      const vectorT &target,
                                      const std::vector<int> &modes,
                                      double rankTolerance,
                                      workspaceT &workspace,
                                      P4PCATiming *timing,
                                      matrixT *coefficients )
{
    P4PCAKernel<double, double>::calculateCenteredInPlace( output,
                                                           predictors,
                                                           target,
                                                           modes,
                                                           rankTolerance,
                                                           workspace,
                                                           timing,
                                                           coefficients );
}

} // namespace improc
} // namespace mx
