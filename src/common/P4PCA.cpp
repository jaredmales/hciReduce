/** \file P4PCA.cpp
 * \brief Implements the pure numerical predictor used by Pixel Prediction Post-Processing.
 * \author Jared R. Males
 */

#include "P4PCA.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
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

    return mx::math::calcEigenVecs<double>( eigenvectors,
                                            eigenvalues,
                                            covariance,
                                            modeCount,
                                            false,
                                            false,
                                            &workspace,
                                            nullptr );
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

/// Subtract each column's double-precision temporal mean in place.
template <typename arrayT>
void p4PCACenterColumns( arrayT &array /**< [in,out] array whose columns are centered */ )
{
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
        array.col( column ) -= mean;
    }
}

/// Calculate one already-validated P4PCA request using its selected fit and residual targets.
void p4PCACalculateValidated( P4PCAResult &output,                  /**< [out] completed regression result */
                              const P4PCA::matrixT &predictors,     /**< [in] predictors used to form the fit */
                              const P4PCA::vectorT &fitTarget,      /**< [in] target used to calculate projections */
                              const P4PCA::vectorT &residualTarget, /**< [in] initial target for output residuals */
                              const std::vector<int> &modes,        /**< [in] requested retained-mode counts */
                              double rankTolerance,                 /**< [in] relative numerical-rank threshold */
                              int maxDOF,            /**< [in] structural degree-of-freedom limit and eigenpair count */
                              bool centerPrediction, /**< [in] enforce a zero-mean fitted prediction */
                              P4PCA::workspaceT &workspace /**< [in,out] reusable LAPACK workspace */ )
{
    const Eigen::Index sampleCount = predictors.rows();
    const Eigen::Index predictorCount = predictors.cols();
    const bool useTemporalGram = sampleCount <= predictorCount;
    P4PCA::matrixT gram;
    P4PCA::vectorT crossProduct;
    if( useTemporalGram )
    {
        gram = ( predictors.matrix() * predictors.matrix().transpose() ).array();
    }
    else
    {
        gram = ( predictors.matrix().transpose() * predictors.matrix() ).array();
        crossProduct = ( predictors.matrix().transpose() * fitTarget.matrix() ).array();
    }

    if( !p4PCAAllFinite( gram ) || ( !useTemporalGram && !p4PCAAllFinite( crossProduct ) ) )
    {
        throw std::runtime_error( "P4PCA normal equations produced nonfinite values" );
    }

    P4PCA::matrixT eigenvectors;
    P4PCA::matrixT eigenvalues;
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

    output.residuals.resize( sampleCount, modes.size() );
    output.residuals.setConstant( std::numeric_limits<double>::quiet_NaN() );
    output.modeStatus.assign( modes.size(), P4PCAModeStatus::rankInsufficient );
    output.numericalRank = numericalRank;

    const auto unsupported = std::upper_bound( modes.begin(), modes.end(), numericalRank );
    if( unsupported == modes.begin() )
    {
        return;
    }

    const int largestSupportedMode = *( unsupported - 1 );
    P4PCA::vectorT residual = residualTarget;
    std::size_t outputIndex{ 0 };

    for( int retainedMode = 1; retainedMode <= largestSupportedMode; ++retainedMode )
    {
        const int eigenIndex = maxDOF - retainedMode;
        const double eigenvalue = eigenvalues( eigenIndex );
        const P4PCA::vectorT eigenvector = eigenvectors.col( eigenIndex );
        double coefficient{ 0 };
        P4PCA::vectorT projectedMode;
        if( useTemporalGram )
        {
            coefficient = eigenvector.matrix().dot( fitTarget.matrix() );
            projectedMode = eigenvector * coefficient;
        }
        else
        {
            coefficient = eigenvector.matrix().dot( crossProduct.matrix() ) / eigenvalue;
            projectedMode = ( predictors.matrix() * eigenvector.matrix() ).array();
            projectedMode *= coefficient;
        }

        if( !mx::math::isFinite( coefficient ) || !p4PCAAllFinite( projectedMode ) )
        {
            throw std::runtime_error( "P4PCA mode projection produced nonfinite values" );
        }

        if( centerPrediction )
        {
            p4PCACenterColumns( projectedMode );
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
            ++outputIndex;
        }
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

void P4PCA::calculate( P4PCAResult &output,
                       const matrixT &predictors,
                       const vectorT &target,
                       const std::vector<int> &modes,
                       double rankTolerance,
                       workspaceT &workspace )
{
    const Eigen::Index sampleCount = predictors.rows();
    const Eigen::Index predictorCount = predictors.cols();
    const int maxDOF = static_cast<int>( std::min( sampleCount, predictorCount ) );
    p4PCAValidateInputs( predictors, target, modes, rankTolerance, maxDOF );
    p4PCACalculateValidated( output, predictors, target, target, modes, rankTolerance, maxDOF, false, workspace );
}

void P4PCA::calculateCentered( P4PCAResult &output,
                               const matrixT &predictors,
                               const vectorT &target,
                               const std::vector<int> &modes,
                               double rankTolerance,
                               workspaceT &workspace )
{
    const Eigen::Index sampleCount = predictors.rows();
    const Eigen::Index predictorCount = predictors.cols();
    if( sampleCount == 1 )
    {
        throw std::invalid_argument( "P4PCA centered regression requires at least two samples" );
    }

    const int maxDOF = static_cast<int>( std::min( predictorCount, sampleCount - 1 ) );
    p4PCAValidateInputs( predictors, target, modes, rankTolerance, maxDOF );

    matrixT centeredPredictors = predictors;
    vectorT centeredTarget = target;
    p4PCACenterColumns( centeredPredictors );
    p4PCACenterColumns( centeredTarget );
    if( !p4PCAAllFinite( centeredPredictors ) || !p4PCAAllFinite( centeredTarget ) )
    {
        throw std::runtime_error( "P4PCA temporal centering produced nonfinite values" );
    }

    p4PCACalculateValidated( output,
                             centeredPredictors,
                             centeredTarget,
                             target,
                             modes,
                             rankTolerance,
                             maxDOF,
                             true,
                             workspace );
}

} // namespace improc
} // namespace mx
