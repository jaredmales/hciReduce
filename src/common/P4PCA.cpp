/** \file P4PCA.cpp
 * \brief Implements the pure numerical predictor used by Pixel Prediction Post-Processing.
 * \author Jared R. Males
 */

#include "P4PCA.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <chrono>
#include <limits>
#include <stdexcept>
#include <string>

#include <mx/math/floatUtils.hpp>

namespace mx
{
namespace improc
{

namespace
{

/// Process-wide test override; production execution leaves this null.
std::atomic<detail::P4PCAEigenSolverT> p4PCAEigenSolverForTesting{ nullptr };

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

/// Validate one P4PCA request against its path-specific structural degree-of-freedom limit.
void p4PCAValidateInputs( const P4PCA::matrixT &predictors, /**< [in] predictor matrix to validate */
                          const P4PCA::vectorT &target,     /**< [in] target vector to validate */
                          const std::vector<int> &modes,    /**< [in] requested retained-mode counts */
                          double rankTolerance,             /**< [in] relative numerical-rank threshold */
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

/// Subtract each column's double-precision temporal mean in place and optionally return those means.
template <typename arrayT>
void p4PCACenterColumns( arrayT &array, /**< [in,out] array whose columns are centered */
                         P4PCA::matrixT *columnMeans = nullptr /**< [out] optional one-row array of removed means */ )
{
    if( columnMeans )
    {
        columnMeans->resize( 1, array.cols() );
        columnMeans->setZero();
    }

    for( Eigen::Index column = 0; column < array.cols(); ++column )
    {
        double scale{ 0 };
        for( Eigen::Index row = 0; row < array.rows(); ++row )
        {
            scale = std::max( scale, std::abs( array( row, column ) ) );
        }

        if( scale == 0 )
        {
            continue;
        }

        double scaledSum{ 0 };
        for( Eigen::Index row = 0; row < array.rows(); ++row )
        {
            scaledSum += array( row, column ) / scale;
        }

        const double scaledMean = std::clamp( scaledSum / static_cast<double>( array.rows() ), -1.0, 1.0 );
        const double mean = scale * scaledMean;
        if( columnMeans )
        {
            ( *columnMeans )( 0, column ) = mean;
        }
        array.col( column ) -= mean;
    }
}

/// Calculate one already-validated P4PCA request using its selected fit and residual targets.
void p4PCACalculateValidated(
    P4PCAResult &output,                  /**< [out] completed regression result */
    const P4PCA::matrixT &predictors,     /**< [in] predictors used to form the fit */
    const P4PCA::vectorT &fitTarget,      /**< [in] target used to calculate projections */
    const P4PCA::vectorT &residualTarget, /**< [in] initial target for output residuals */
    const std::vector<int> &modes,        /**< [in] requested retained-mode counts */
    double rankTolerance,                 /**< [in] relative numerical-rank threshold */
    int maxDOF,                           /**< [in] structural degree-of-freedom limit and eigenpair count */
    const P4PCA::matrixT *predictorMeans, /**< [in] optional means for uncentered application */
    P4PCA::workspaceT &workspace,         /**< [in,out] reusable LAPACK workspace */
    P4PCATiming *timing,                  /**< [out] optional per-call worker timing */
    P4PCA::matrixT *coefficients /**< [out] optional predictor-space coefficient vectors */ )
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
    P4PCA::matrixT gram;
    P4PCA::vectorT crossProduct;
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

    P4PCA::vectorT predictorMeanProjection;
    if( predictorMeans )
    {
        if( useTemporalGram )
        {
            predictorMeanProjection = ( predictors.matrix() * predictorMeans->matrix().transpose() ).array();
        }
    }

    P4PCA::matrixT eigenvectors;
    P4PCA::matrixT eigenvalues;
    const auto eigensolveStart = std::chrono::steady_clock::now();
    const MXLAPACK_INT solverStatus = p4PCAEigenSolve( eigenvectors, eigenvalues, gram, maxDOF, workspace );
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

    const double largestEigenvalue = eigenvalues( maxDOF - 1 );
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
    P4PCA::vectorT residual = residualTarget;
    P4PCA::vectorT accumulatedCoefficients;
    P4PCA::vectorT coefficientUpdate;
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
        const double eigenvalue = eigenvalues( eigenIndex );
        const P4PCA::vectorT eigenvector = eigenvectors.col( eigenIndex );
        double coefficient{ 0 };
        double predictionMean{ 0 };
        P4PCA::vectorT projectedMode;
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
            output.residuals.col( outputIndex ) = residual;
            output.modeStatus[outputIndex] = P4PCAModeStatus::rankSupported;
            if( coefficients )
            {
                coefficients->col( outputIndex ) = accumulatedCoefficients;
            }
            ++outputIndex;
        }
    }
    if( timing )
    {
        timing->projectionWorkerSeconds = secondsSince( projectionStart );
    }
}

} // namespace

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
    output.sampleValidity.resize( 0, 0 );
    const Eigen::Index sampleCount = predictors.rows();
    const Eigen::Index predictorCount = predictors.cols();
    const int maxDOF = static_cast<int>( std::min( sampleCount, predictorCount ) );
    p4PCAValidateInputs( predictors, target, modes, rankTolerance, maxDOF );
    p4PCACalculateValidated( output,
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

void P4PCA::calculateHeldOut( P4PCAResult &output,
                              const matrixT &predictors,
                              const vectorT &target,
                              const std::vector<std::vector<std::size_t>> &referenceRows,
                              const std::vector<int> &modes,
                              double rankTolerance,
                              workspaceT &workspace,
                              P4PCATiming *timing )
{
    const Eigen::Index sampleCount = predictors.rows();
    const Eigen::Index predictorCount = predictors.cols();
    if( sampleCount <= 0 || predictorCount <= 0 || target.rows() != sampleCount ||
        referenceRows.size() != static_cast<std::size_t>( sampleCount ) )
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

    output.residuals.resize( sampleCount, static_cast<Eigen::Index>( modes.size() ) );
    output.residuals.setConstant( std::numeric_limits<double>::quiet_NaN() );
    output.sampleValidity.resize( sampleCount, static_cast<Eigen::Index>( modes.size() ) );
    output.sampleValidity.setZero();
    output.modeStatus.assign( modes.size(), P4PCAModeStatus::rankInsufficient );
    output.numericalRank = std::numeric_limits<int>::max();
    if( timing )
    {
        *timing = P4PCATiming{};
    }

    matrixT trainingPredictors;
    vectorT trainingTarget;
    matrixT coefficients;
    for( Eigen::Index heldOut = 0; heldOut < sampleCount; ++heldOut )
    {
        const std::vector<std::size_t> &rows = referenceRows[static_cast<std::size_t>( heldOut )];
        if( rows.empty() )
        {
            output.numericalRank = 0;
            continue;
        }
        std::size_t previousRow{ 0 };
        bool first{ true };
        for( const std::size_t row : rows )
        {
            if( row >= static_cast<std::size_t>( sampleCount ) || row == static_cast<std::size_t>( heldOut ) ||
                ( !first && row <= previousRow ) )
            {
                throw std::invalid_argument(
                    "P4PCA held-out reference rows must be ordered, unique, in range, and omit their target" );
            }
            previousRow = row;
            first = false;
        }

        const int structuralRank =
            static_cast<int>( std::min<std::size_t>( rows.size(), static_cast<std::size_t>( predictorCount ) ) );
        output.numericalRank = std::min( output.numericalRank, structuralRank );
        const auto supportedEnd = std::upper_bound( modes.begin(), modes.end(), structuralRank );
        if( supportedEnd == modes.begin() )
        {
            continue;
        }
        const std::vector<int> supportedModes( modes.begin(), supportedEnd );
        trainingPredictors.resize( static_cast<Eigen::Index>( rows.size() ), predictorCount );
        trainingTarget.resize( static_cast<Eigen::Index>( rows.size() ) );
        for( std::size_t training = 0; training < rows.size(); ++training )
        {
            const Eigen::Index source = static_cast<Eigen::Index>( rows[training] );
            trainingPredictors.row( static_cast<Eigen::Index>( training ) ) = predictors.row( source );
            trainingTarget( static_cast<Eigen::Index>( training ) ) = target( source );
        }

        P4PCAResult fit;
        P4PCATiming fitTiming;
        calculate( fit,
                   trainingPredictors,
                   trainingTarget,
                   supportedModes,
                   rankTolerance,
                   workspace,
                   &fitTiming,
                   &coefficients );
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
            const double prediction = predictors.row( heldOut ).matrix().dot(
                coefficients.col( static_cast<Eigen::Index>( mode ) ).matrix() );
            const double residual = target( heldOut ) - prediction;
            if( !mx::math::isFinite( residual ) )
            {
                throw std::runtime_error( "P4PCA held-out prediction produced a nonfinite residual" );
            }
            output.residuals( heldOut, static_cast<Eigen::Index>( mode ) ) = residual;
            output.sampleValidity( heldOut, static_cast<Eigen::Index>( mode ) ) = 1;
            output.modeStatus[mode] = P4PCAModeStatus::rankSupported;
        }
        if( timing )
        {
            timing->projectionWorkerSeconds +=
                std::chrono::duration<double>( std::chrono::steady_clock::now() - heldOutProjectionBegin ).count();
        }
    }
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

void P4PCA::calculateCenteredInPlace( P4PCAResult &output,
                                      matrixT &predictors,
                                      const vectorT &target,
                                      const std::vector<int> &modes,
                                      double rankTolerance,
                                      workspaceT &workspace,
                                      P4PCATiming *timing,
                                      matrixT *coefficients )
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

    p4PCACalculateValidated( output,
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

} // namespace improc
} // namespace mx
