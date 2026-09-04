/** \file P4TemporalPCA.cpp
 * \brief Implements time-domain PCA prediction for Pixel Prediction Post-Processing.
 * \author Jared R. Males
 */

#include "P4TemporalPCA.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#include <mx/math/floatUtils.hpp>

namespace mx
{
namespace improc
{

template <typename arrayT>
bool P4TemporalPCA::allFinite( const arrayT &array )
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

void P4TemporalPCA::configure( const matrixT &referenceSeries,
                               int requestedModes,
                               P4TemporalPCACentering centering,
                               double rankTolerance,
                               int gapImages,
                               workspaceT &workspace )
{
    if( referenceSeries.rows() <= 0 || referenceSeries.cols() <= 0 || requestedModes <= 0 || gapImages < 0 ||
        !mx::math::isFinite( rankTolerance ) || rankTolerance < 0 || rankTolerance >= 1 )
    {
        throw std::invalid_argument( "P4 temporal PCA configuration inputs are invalid" );
    }
    if( requestedModes > std::min( referenceSeries.rows(), referenceSeries.cols() ) ||
        referenceSeries.cols() > static_cast<Eigen::Index>( std::numeric_limits<int>::max() ) )
    {
        throw std::invalid_argument( "P4 temporal PCA requested mode count is structurally invalid" );
    }
    if( centering != P4TemporalPCACentering::pixelMean && centering != P4TemporalPCACentering::none )
    {
        throw std::invalid_argument( "P4 temporal PCA centering policy is invalid" );
    }
    if( !allFinite( referenceSeries ) )
    {
        throw std::invalid_argument( "P4 temporal PCA reference series must be finite" );
    }

    matrixT centered = referenceSeries;
    if( centering == P4TemporalPCACentering::pixelMean )
    {
        for( Eigen::Index reference = 0; reference < centered.rows(); ++reference )
        {
            double mean = centered( reference, 0 );
            for( Eigen::Index sample = 1; sample < centered.cols(); ++sample )
            {
                mean += ( centered( reference, sample ) - mean ) / static_cast<double>( sample + 1 );
            }
            for( Eigen::Index sample = 0; sample < centered.cols(); ++sample )
            {
                centered( reference, sample ) -= mean;
            }
        }
    }
    if( !allFinite( centered ) )
    {
        throw std::runtime_error( "P4 temporal PCA centering produced a nonfinite reference series" );
    }

    matrixT covariance = ( centered.matrix().transpose() * centered.matrix() ).array();
    if( !allFinite( covariance ) )
    {
        throw std::runtime_error( "P4 temporal PCA covariance is nonfinite" );
    }

    matrixT eigenvectors;
    matrixT eigenvalues;
    const int dimension = static_cast<int>( covariance.rows() );
    const MXLAPACK_INT solverStatus = mx::math::eigenSYEVR( eigenvectors,
                                                            eigenvalues,
                                                            covariance,
                                                            dimension - requestedModes,
                                                            dimension,
                                                            'L',
                                                            &workspace );
    if( solverStatus != 0 )
    {
        throw std::runtime_error( "P4 temporal PCA eigensolver failed with status " + std::to_string( solverStatus ) );
    }
    if( eigenvectors.rows() != dimension || eigenvectors.cols() != requestedModes ||
        eigenvalues.rows() < requestedModes || eigenvalues.cols() < 1 || !allFinite( eigenvectors ) ||
        !allFinite( eigenvalues.block( 0, 0, requestedModes, 1 ) ) )
    {
        throw std::runtime_error( "P4 temporal PCA eigensolver returned invalid output" );
    }
    for( int mode = 1; mode < requestedModes; ++mode )
    {
        if( eigenvalues( mode ) < eigenvalues( mode - 1 ) )
        {
            throw std::runtime_error( "P4 temporal PCA eigensolver returned unsorted eigenvalues" );
        }
    }
    const double largestEigenvalue = eigenvalues( requestedModes - 1 );
    m_numericalRank = 0;
    if( largestEigenvalue > 0 )
    {
        const long double threshold = static_cast<long double>( rankTolerance ) * largestEigenvalue;
        for( int mode = 0; mode < requestedModes; ++mode )
        {
            if( static_cast<long double>( eigenvalues( mode ) ) > threshold )
            {
                ++m_numericalRank;
            }
        }
    }
    m_modes = std::move( eigenvectors );
    m_requestedModes = requestedModes;
    m_gapFits.clear();
    m_retainedRows.clear();
    m_gapFits.resize( static_cast<std::size_t>( dimension ) );
    m_retainedRows.resize( static_cast<std::size_t>( dimension ) );
    for( int target = 0; target < dimension; ++target )
    {
        const int gapBegin = gapImages > target ? 0 : target - gapImages;
        const int followingSamples = dimension - 1 - target;
        const int gapEnd = gapImages > followingSamples ? dimension - 1 : target + gapImages;
        std::vector<Eigen::Index> &retained = m_retainedRows[static_cast<std::size_t>( target )];
        retained.reserve( static_cast<std::size_t>( dimension - ( gapEnd - gapBegin + 1 ) ) );
        for( int sample = 0; sample < dimension; ++sample )
        {
            if( sample < gapBegin || sample > gapEnd )
            {
                retained.push_back( sample );
            }
        }
        if( retained.size() <= static_cast<std::size_t>( requestedModes ) )
        {
            continue;
        }
        Eigen::MatrixXd design( static_cast<Eigen::Index>( retained.size() ), requestedModes + 1 );
        for( std::size_t row = 0; row < retained.size(); ++row )
        {
            design( static_cast<Eigen::Index>( row ), 0 ) = 1;
            design.row( static_cast<Eigen::Index>( row ) ).tail( requestedModes ) =
                m_modes.row( retained[row] ).matrix();
        }
        m_gapFits[static_cast<std::size_t>( target )].compute( design );
        m_gapFits[static_cast<std::size_t>( target )].setThreshold( rankTolerance );
    }
}

void P4TemporalPCA::predict( P4TemporalPCAResult &result, const vectorT &targetSeries ) const
{
    if( m_requestedModes <= 0 || m_modes.rows() <= 0 ||
        m_gapFits.size() != static_cast<std::size_t>( m_modes.rows() ) || m_retainedRows.size() != m_gapFits.size() ||
        targetSeries.rows() != m_modes.rows() || !allFinite( targetSeries ) )
    {
        throw std::invalid_argument( "P4 temporal PCA prediction inputs are invalid" );
    }

    result.predictions.resize( m_modes.rows() );
    result.predictions.setConstant( std::numeric_limits<double>::quiet_NaN() );
    result.residuals.resize( m_modes.rows() );
    result.residuals.setConstant( std::numeric_limits<double>::quiet_NaN() );
    result.validity.resize( m_modes.rows() );
    result.validity.setZero();
    result.numericalRank = m_numericalRank;
    if( m_numericalRank != m_requestedModes )
    {
        return;
    }

    for( Eigen::Index target = 0; target < m_modes.rows(); ++target )
    {
        const std::size_t index = static_cast<std::size_t>( target );
        const std::vector<Eigen::Index> &retained = m_retainedRows[index];
        const Eigen::ColPivHouseholderQR<Eigen::MatrixXd> &fit = m_gapFits[index];
        if( retained.size() <= static_cast<std::size_t>( m_requestedModes ) || fit.rank() != m_requestedModes + 1 )
        {
            continue;
        }
        Eigen::VectorXd samples( static_cast<Eigen::Index>( retained.size() ) );
        for( std::size_t row = 0; row < retained.size(); ++row )
        {
            samples( static_cast<Eigen::Index>( row ) ) = targetSeries( retained[row] );
        }
        const Eigen::VectorXd coefficients = fit.solve( samples );
        const double prediction =
            coefficients( 0 ) + m_modes.row( target ).matrix().dot( coefficients.tail( m_requestedModes ) );
        if( !mx::math::isFinite( prediction ) )
        {
            throw std::runtime_error( "P4 temporal PCA prediction is nonfinite" );
        }
        result.predictions( target ) = prediction;
        result.residuals( target ) = targetSeries( target ) - prediction;
        result.validity( target ) = 1;
    }
}

Eigen::Index P4TemporalPCA::sampleCount() const noexcept
{
    return m_modes.rows();
}

int P4TemporalPCA::requestedModes() const noexcept
{
    return m_requestedModes;
}

} // namespace improc
} // namespace mx
