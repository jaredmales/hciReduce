/** \file PSFNoiseModel.cpp
 * \brief Implements regularized residual-noise covariance estimation and restricted solves.
 */

#include "PSFNoiseModel.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace mx
{
namespace improc
{
namespace
{

/// Center finite noise samples after validating the sample count and absolute variance floor.
PSFNoiseModel::matrixT centeredNoiseSamples( PSFNoiseModel::vectorT &mean, /**< [out] sample mean in pixel order */
                                             const PSFNoiseModel::matrixT &samples, /**< [in] noise-only sample rows */
                                             double varianceFloor /**< [in] positive regularization variance */ )
{
    if( samples.rows() < 2 || samples.cols() <= 0 || !samples.allFinite() || !std::isfinite( varianceFloor ) ||
        varianceFloor <= 0 )
    {
        throw std::invalid_argument(
            "noise training requires at least two finite samples and a positive variance floor" );
    }
    mean = samples.colwise().mean().transpose();
    PSFNoiseModel::matrixT centered = samples.rowwise() - mean.transpose();
    if( !mean.allFinite() || !centered.allFinite() )
    {
        throw std::overflow_error( "noise sample centering is nonfinite" );
    }
    return centered;
}

} // namespace

PSFNoiseModel::PSFNoiseModel( PSFNoiseKind kind, vectorT variances, matrixT factor, vectorT mean )
    : m_kind( kind ), m_variances( std::move( variances ) ), m_factor( std::move( factor ) ),
      m_mean( std::move( mean ) )
{
    if( m_variances.size() <= 0 || !m_variances.allFinite() || ( m_variances.array() <= 0 ).any() ||
        m_factor.rows() != m_variances.size() || !m_factor.allFinite() )
    {
        throw std::invalid_argument(
            "noise covariance requires positive finite variances and a compatible finite factor" );
    }
    if( m_mean.size() == 0 )
    {
        m_mean = vectorT::Zero( m_variances.size() );
    }
    if( m_mean.size() != m_variances.size() || !m_mean.allFinite() )
    {
        throw std::invalid_argument( "noise background mean must be finite and match the stamp-pixel count" );
    }
}

PSFNoiseModel PSFNoiseModel::identity( Eigen::Index pixels )
{
    if( pixels <= 0 )
    {
        throw std::invalid_argument( "identity noise requires a positive stamp-pixel count" );
    }
    return PSFNoiseModel( PSFNoiseKind::identity, vectorT::Ones( pixels ), matrixT( pixels, 0 ), vectorT{} );
}

PSFNoiseModel PSFNoiseModel::diagonal( const vectorT &variances, const vectorT &mean )
{
    return PSFNoiseModel( PSFNoiseKind::diagonal, variances, matrixT( variances.size(), 0 ), mean );
}

PSFNoiseModel PSFNoiseModel::lowRank( const vectorT &variances, const matrixT &factor, const vectorT &mean )
{
    return PSFNoiseModel( PSFNoiseKind::pca, variances, factor, mean );
}

PSFNoiseModel PSFNoiseModel::estimateDiagonal( const matrixT &samples, double varianceFloor )
{
    vectorT mean;
    const matrixT centered = centeredNoiseSamples( mean, samples, varianceFloor );
    const vectorT variances = centered.colwise().squaredNorm().transpose() / static_cast<double>( samples.rows() - 1 );
    if( !variances.allFinite() )
    {
        throw std::overflow_error( "noise sample variances are nonfinite" );
    }
    return diagonal( variances.cwiseMax( varianceFloor ), mean );
}

PSFNoiseModel PSFNoiseModel::estimatePCA( const matrixT &samples, std::size_t maximumModes, double varianceFloor )
{
    vectorT mean;
    const matrixT centered = centeredNoiseSamples( mean, samples, varianceFloor );
    const Eigen::Index count = static_cast<Eigen::Index>(
        std::min( maximumModes, static_cast<std::size_t>( std::min( samples.rows() - 1, samples.cols() ) ) ) );
    if( count == 0 )
    {
        return lowRank( vectorT::Constant( samples.cols(), varianceFloor ), matrixT( samples.cols(), 0 ), mean );
    }
    Eigen::JacobiSVD<matrixT> decomposition( centered, Eigen::ComputeThinV );
    if( decomposition.info() != Eigen::Success || !decomposition.singularValues().allFinite() )
    {
        throw std::runtime_error( "noise sample PCA failed" );
    }
    const vectorT variances =
        decomposition.singularValues().head( count ).array().square() / static_cast<double>( samples.rows() - 1 );
    if( !variances.allFinite() )
    {
        throw std::overflow_error( "noise PCA variances are nonfinite" );
    }
    Eigen::Index retained{ 0 };
    while( retained < count && variances( retained ) > varianceFloor )
    {
        ++retained;
    }
    const matrixT factor = decomposition.matrixV().leftCols( retained ) *
                           ( variances.head( retained ).array() - varianceFloor ).sqrt().matrix().asDiagonal();
    return lowRank( vectorT::Constant( samples.cols(), varianceFloor ), factor, mean );
}

PSFNoiseKind PSFNoiseModel::kind() const
{
    return m_kind;
}

Eigen::Index PSFNoiseModel::pixels() const
{
    return m_variances.size();
}

Eigen::Index PSFNoiseModel::modes() const
{
    return m_factor.cols();
}

const PSFNoiseModel::vectorT &PSFNoiseModel::mean() const
{
    return m_mean;
}

PSFNoiseModel::vectorT PSFNoiseModel::solve( const vectorT &value, const std::vector<Eigen::Index> &support ) const
{
    if( support.size() != static_cast<std::size_t>( value.size() ) || !value.allFinite() )
    {
        throw std::invalid_argument( "noise solve requires one finite vector entry per retained pixel" );
    }
    std::vector<bool> seen( static_cast<std::size_t>( pixels() ), false );
    vectorT scale( value.size() );
    matrixT factor( value.size(), modes() );
    for( Eigen::Index row = 0; row < value.size(); ++row )
    {
        const Eigen::Index pixel = support[static_cast<std::size_t>( row )];
        if( pixel < 0 || pixel >= pixels() || seen[static_cast<std::size_t>( pixel )] )
        {
            throw std::invalid_argument( "noise solve support must contain distinct in-range stamp indices" );
        }
        seen[static_cast<std::size_t>( pixel )] = true;
        scale( row ) = 1 / std::sqrt( m_variances( pixel ) );
        factor.row( row ) = scale( row ) * m_factor.row( pixel );
    }
    vectorT scaled = scale.array() * value.array();
    if( modes() > 0 && value.size() > 0 )
    {
        if( !factor.allFinite() || !scaled.allFinite() )
        {
            throw std::overflow_error( "scaled noise covariance is nonfinite" );
        }
        // For B = D^-1/2 L = Q [R;0], solve I + R R^T in the leading Q coordinates.
        // This is equivalent to Woodbury without subtracting nearly equal vectors when B is large.
        const Eigen::Index count = std::min( value.size(), modes() );
        const Eigen::HouseholderQR<matrixT> rotation( factor );
        const matrixT triangular = rotation.matrixQR().topRows( count ).template triangularView<Eigen::Upper>();
        const matrixT reduced = matrixT::Identity( count, count ) + triangular * triangular.transpose();
        if( !reduced.allFinite() )
        {
            throw std::overflow_error( "reduced noise covariance is nonfinite" );
        }
        Eigen::LLT<matrixT> decomposition( reduced );
        if( decomposition.info() != Eigen::Success )
        {
            throw std::runtime_error( "regularized noise covariance factorization failed" );
        }
        vectorT rotated = rotation.householderQ().adjoint() * scaled;
        rotated.head( count ) = decomposition.solve( rotated.head( count ) ).eval();
        scaled = rotation.householderQ() * rotated;
    }
    const vectorT result = scale.array() * scaled.array();
    if( !result.allFinite() )
    {
        throw std::overflow_error( "noise covariance solve is nonfinite" );
    }
    return result;
}

} // namespace improc
} // namespace mx
