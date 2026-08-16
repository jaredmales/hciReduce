/** \file P4PCA.cpp
 * \brief Implements the pure numerical predictor used by Pixel Prediction Post-Processing.
 * \author Jared R. Males
 */

#include "P4PCA.hpp"

#include <algorithm>
#include <atomic>
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

    const int maxDOF = static_cast<int>( std::min( sampleCount, predictorCount ) );
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

    const bool useTemporalGram = sampleCount <= predictorCount;
    matrixT gram;
    vectorT crossProduct;
    if( useTemporalGram )
    {
        gram = ( predictors.matrix() * predictors.matrix().transpose() ).array();
    }
    else
    {
        gram = ( predictors.matrix().transpose() * predictors.matrix() ).array();
        crossProduct = ( predictors.matrix().transpose() * target.matrix() ).array();
    }

    if( !p4PCAAllFinite( gram ) || ( !useTemporalGram && !p4PCAAllFinite( crossProduct ) ) )
    {
        throw std::runtime_error( "P4PCA normal equations produced nonfinite values" );
    }

    matrixT eigenvectors;
    matrixT eigenvalues;
    const MXLAPACK_INT solverStatus = p4PCAEigenSolve( eigenvectors, eigenvalues, gram, maxDOF, workspace );
    if( solverStatus != 0 )
    {
        throw std::runtime_error( "P4PCA eigensolver failed with status " + std::to_string( solverStatus ) );
    }

    if( eigenvectors.rows() != maxDOF || eigenvectors.cols() != maxDOF || eigenvalues.rows() < maxDOF ||
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
    vectorT residual = target;
    std::size_t outputIndex{ 0 };

    for( int retainedMode = 1; retainedMode <= largestSupportedMode; ++retainedMode )
    {
        const int eigenIndex = maxDOF - retainedMode;
        const double eigenvalue = eigenvalues( eigenIndex );
        const vectorT eigenvector = eigenvectors.col( eigenIndex );
        double coefficient{ 0 };
        vectorT projectedMode;
        if( useTemporalGram )
        {
            coefficient = eigenvector.matrix().dot( target.matrix() );
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

} // namespace improc
} // namespace mx
