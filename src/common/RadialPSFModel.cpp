/** \file RadialPSFModel.cpp
 * \brief Implements azimuthally averaged PSF response models sampled at discrete radii.
 * \author Jared R. Males
 */

#include "RadialPSFModel.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <unordered_set>

#include "P4PixelGrid.hpp"

namespace mx
{
namespace improc
{

namespace
{

/// Validate a radial grid shared by selection and model construction.
void validateRadialPSFRadii( const std::vector<double> &radii /**< [in] candidate radial grid */ )
{
    if( radii.empty() )
    {
        throw std::invalid_argument( "radial PSF model requires at least one radius" );
    }
    for( std::size_t index = 0; index < radii.size(); ++index )
    {
        if( !std::isfinite( radii[index] ) || radii[index] < 0 || ( index != 0 && radii[index] <= radii[index - 1] ) )
        {
            throw std::invalid_argument( "radial PSF radii must be finite, nonnegative, and strictly increasing" );
        }
    }
}

/// Return whether a response contains finite values at every valid sample.
bool validRadialPSFValues( const RadialPSFModel::imageT &response, /**< [in] response values */
                           const RadialPSFModel::validityT &validity /**< [in] response validity */ )
{
    for( Eigen::Index column = 0; column < response.cols(); ++column )
    {
        for( Eigen::Index row = 0; row < response.rows(); ++row )
        {
            if( validity( row, column ) != 0 && !std::isfinite( response( row, column ) ) )
            {
                return false;
            }
        }
    }
    return true;
}

} // namespace

RadialPSFModel::RadialPSFModel( std::vector<double> radii, int stampRows, int stampColumns )
    : m_radii( std::move( radii ) ), m_stampRows( stampRows ), m_stampColumns( stampColumns )
{
    validateRadialPSFRadii( m_radii );
    if( stampRows <= 0 || stampColumns <= 0 )
    {
        throw std::invalid_argument( "radial PSF stamp dimensions must be positive" );
    }
    m_responses.resize( m_radii.size() );
    m_validities.resize( m_radii.size() );
    m_sampleCounts.assign( m_radii.size(), 0 );
}

std::vector<RadialPSFSample> RadialPSFModel::selectSamples( const std::vector<RadialPSFSource> &sources,
                                                            double centerRow,
                                                            double centerColumn,
                                                            const std::vector<double> &radii,
                                                            std::size_t samplesPerRadius )
{
    validateRadialPSFRadii( radii );
    if( sources.empty() )
    {
        throw std::invalid_argument( "radial PSF sampling requires at least one available source coordinate" );
    }
    if( !std::isfinite( centerRow ) || !std::isfinite( centerColumn ) )
    {
        throw std::invalid_argument( "radial PSF sampling center must be finite" );
    }
    if( samplesPerRadius == 0 )
    {
        throw std::invalid_argument( "radial PSF sampling requires a positive angular sample count" );
    }
    std::unordered_set<std::size_t> sourceIndices;
    for( const RadialPSFSource &source : sources )
    {
        if( !std::isfinite( source.row ) || !std::isfinite( source.column ) ||
            !sourceIndices.insert( source.sourceIndex ).second )
        {
            throw std::invalid_argument( "radial PSF source coordinates must be finite with unique indices" );
        }
    }

    std::vector<RadialPSFSample> samples;
    if( radii.size() > std::numeric_limits<std::size_t>::max() / samplesPerRadius )
    {
        throw std::length_error( "radial PSF sample count overflows size_t" );
    }
    samples.reserve( radii.size() * samplesPerRadius );
    for( std::size_t radiusIndex = 0; radiusIndex < radii.size(); ++radiusIndex )
    {
        std::unordered_set<std::size_t> selectedAtRadius;
        for( std::size_t angularIndex = 0; angularIndex < samplesPerRadius; ++angularIndex )
        {
            const double requestedAngle =
                2.0 * std::numbers::pi * static_cast<double>( angularIndex ) / static_cast<double>( samplesPerRadius );
            const double requestedRow = centerRow + radii[radiusIndex] * std::sin( requestedAngle );
            const double requestedColumn = centerColumn + radii[radiusIndex] * std::cos( requestedAngle );
            const RadialPSFSource *nearest = nullptr;
            double nearestDistanceSquared{ 0 };
            for( const RadialPSFSource &source : sources )
            {
                const double deltaRow = source.row - requestedRow;
                const double deltaColumn = source.column - requestedColumn;
                const double distanceSquared = deltaRow * deltaRow + deltaColumn * deltaColumn;
                if( nearest == nullptr || distanceSquared < nearestDistanceSquared ||
                    ( distanceSquared == nearestDistanceSquared && source.sourceIndex < nearest->sourceIndex ) )
                {
                    nearest = &source;
                    nearestDistanceSquared = distanceSquared;
                }
            }
            if( selectedAtRadius.insert( nearest->sourceIndex ).second )
            {
                const double deltaRow = nearest->row - centerRow;
                const double deltaColumn = nearest->column - centerColumn;
                samples.push_back( { nearest->sourceIndex,
                                     radiusIndex,
                                     std::hypot( deltaRow, deltaColumn ),
                                     std::atan2( deltaRow, deltaColumn ) } );
            }
        }
        if( selectedAtRadius.empty() )
        {
            throw std::runtime_error( "radial PSF sampling did not select a source coordinate" );
        }
    }
    return samples;
}

void RadialPSFModel::fit( const std::vector<imageT> &responses,
                          const std::vector<validityT> &validities,
                          const std::vector<RadialPSFSample> &samples )
{
    if( responses.size() != samples.size() || validities.size() != samples.size() || samples.empty() )
    {
        throw std::invalid_argument( "radial PSF measurements and metadata must be nonempty and have equal lengths" );
    }

    std::vector<Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>> sums( m_radii.size() );
    std::vector<Eigen::Array<std::size_t, Eigen::Dynamic, Eigen::Dynamic>> counts( m_radii.size() );
    m_sampleCounts.assign( m_radii.size(), 0 );
    for( std::size_t radiusIndex = 0; radiusIndex < m_radii.size(); ++radiusIndex )
    {
        sums[radiusIndex] = Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>::Zero( m_stampRows, m_stampColumns );
        counts[radiusIndex] =
            Eigen::Array<std::size_t, Eigen::Dynamic, Eigen::Dynamic>::Zero( m_stampRows, m_stampColumns );
    }

    for( std::size_t sampleIndex = 0; sampleIndex < samples.size(); ++sampleIndex )
    {
        const imageT &measured = responses[sampleIndex];
        const validityT &measuredValidity = validities[sampleIndex];
        const RadialPSFSample &sample = samples[sampleIndex];
        if( sample.radiusIndex >= m_radii.size() || measured.rows() != m_stampRows ||
            measured.cols() != m_stampColumns || measuredValidity.rows() != m_stampRows ||
            measuredValidity.cols() != m_stampColumns || !std::isfinite( sample.radius ) || sample.radius < 0 ||
            !std::isfinite( sample.angle ) || !validRadialPSFValues( measured, measuredValidity ) )
        {
            throw std::invalid_argument( "radial PSF measurement is inconsistent with configured geometry" );
        }
        imageT canonical;
        validityT canonicalValidity;
        rotate( canonical, canonicalValidity, measured, measuredValidity, sample.angle );
        for( int column = 0; column < m_stampColumns; ++column )
        {
            for( int row = 0; row < m_stampRows; ++row )
            {
                if( canonicalValidity( row, column ) != 0 )
                {
                    sums[sample.radiusIndex]( row, column ) += canonical( row, column );
                    ++counts[sample.radiusIndex]( row, column );
                }
            }
        }
        ++m_sampleCounts[sample.radiusIndex];
    }

    for( std::size_t radiusIndex = 0; radiusIndex < m_radii.size(); ++radiusIndex )
    {
        if( m_sampleCounts[radiusIndex] == 0 )
        {
            throw std::invalid_argument( "every radial PSF bin requires at least one measurement" );
        }
        m_responses[radiusIndex].resize( m_stampRows, m_stampColumns );
        m_responses[radiusIndex].setZero();
        m_validities[radiusIndex].resize( m_stampRows, m_stampColumns );
        m_validities[radiusIndex].setZero();
        for( int column = 0; column < m_stampColumns; ++column )
        {
            for( int row = 0; row < m_stampRows; ++row )
            {
                const std::size_t count = counts[radiusIndex]( row, column );
                if( count != 0 )
                {
                    m_responses[radiusIndex]( row, column ) =
                        static_cast<float>( sums[radiusIndex]( row, column ) / static_cast<double>( count ) );
                    m_validities[radiusIndex]( row, column ) = 1;
                }
            }
        }
    }
}

void RadialPSFModel::rotate(
    imageT &output, validityT &outputValidity, const imageT &input, const validityT &inputValidity, double angle )
{
    if( input.rows() <= 0 || input.cols() <= 0 || input.rows() != inputValidity.rows() ||
        input.cols() != inputValidity.cols() )
    {
        throw std::invalid_argument( "radial PSF rotation requires matching nonempty response and validity arrays" );
    }
    if( !std::isfinite( angle ) || !validRadialPSFValues( input, inputValidity ) )
    {
        throw std::invalid_argument( "radial PSF rotation requires finite geometry and valid response samples" );
    }

    output.resize( input.rows(), input.cols() );
    output.setZero();
    outputValidity.resize( input.rows(), input.cols() );
    outputValidity.setZero();
    const double centerRow = 0.5 * static_cast<double>( input.rows() - 1 );
    const double centerColumn = 0.5 * static_cast<double>( input.cols() - 1 );
    const double cosine = std::cos( angle );
    const double sine = std::sin( angle );
    using gridT = P4PixelGridf;
    gridT::transformT transform;
    for( int outputColumn = 0; outputColumn < output.cols(); ++outputColumn )
    {
        const double deltaColumn = static_cast<double>( outputColumn ) - centerColumn;
        for( int outputRow = 0; outputRow < output.rows(); ++outputRow )
        {
            const double deltaRow = static_cast<double>( outputRow ) - centerRow;
            const double inputRow = centerRow + deltaRow * cosine + deltaColumn * sine;
            const double inputColumn = centerColumn - deltaRow * sine + deltaColumn * cosine;
            const double floorRow = std::floor( inputRow );
            const double floorColumn = std::floor( inputColumn );
            if( floorRow < static_cast<double>( std::numeric_limits<int>::min() + gridT::leftBuffer ) ||
                floorRow > static_cast<double>( std::numeric_limits<int>::max() - gridT::width ) ||
                floorColumn < static_cast<double>( std::numeric_limits<int>::min() + gridT::leftBuffer ) ||
                floorColumn > static_cast<double>( std::numeric_limits<int>::max() - gridT::width ) )
            {
                continue;
            }
            const int footprintRow = static_cast<int>( floorRow ) - gridT::leftBuffer;
            const int footprintColumn = static_cast<int>( floorColumn ) - gridT::leftBuffer;
            gridT::kernelT kernel;
            transform( kernel,
                       static_cast<float>( inputRow - floorRow ),
                       static_cast<float>( inputColumn - floorColumn ) );
            bool valid{ true };
            float value{ 0 };
            for( int columnOffset = 0; columnOffset < gridT::width && valid; ++columnOffset )
            {
                for( int rowOffset = 0; rowOffset < gridT::width; ++rowOffset )
                {
                    const float weight = kernel( rowOffset, columnOffset );
                    if( weight == 0 )
                    {
                        continue;
                    }
                    const int inputSampleRow = footprintRow + rowOffset;
                    const int inputSampleColumn = footprintColumn + columnOffset;
                    if( inputSampleRow < 0 || inputSampleRow >= input.rows() || inputSampleColumn < 0 ||
                        inputSampleColumn >= input.cols() )
                    {
                        continue;
                    }
                    if( inputValidity( inputSampleRow, inputSampleColumn ) == 0 )
                    {
                        valid = false;
                        break;
                    }
                    value += input( inputSampleRow, inputSampleColumn ) * weight;
                }
            }
            if( valid && std::isfinite( value ) )
            {
                output( outputRow, outputColumn ) = value;
                outputValidity( outputRow, outputColumn ) = 1;
            }
        }
    }
}

void RadialPSFModel::response( imageT &output, validityT &outputValidity, double requestedRadius, double angle ) const
{
    if( !std::isfinite( angle ) )
    {
        throw std::invalid_argument( "radial PSF response angle must be finite" );
    }
    const std::size_t radiusIndex = nearestRadiusIndex( requestedRadius );
    if( m_responses[radiusIndex].size() == 0 )
    {
        throw std::logic_error( "radial PSF model must be fitted before evaluating a response" );
    }
    rotate( output, outputValidity, m_responses[radiusIndex], m_validities[radiusIndex], -angle );
}

std::size_t RadialPSFModel::radiusCount() const noexcept
{
    return m_radii.size();
}

double RadialPSFModel::radius( std::size_t radiusIndex ) const
{
    return m_radii.at( radiusIndex );
}

const RadialPSFModel::imageT &RadialPSFModel::canonicalResponse( std::size_t radiusIndex ) const
{
    return m_responses.at( radiusIndex );
}

const RadialPSFModel::validityT &RadialPSFModel::canonicalValidity( std::size_t radiusIndex ) const
{
    return m_validities.at( radiusIndex );
}

std::size_t RadialPSFModel::sampleCount( std::size_t radiusIndex ) const
{
    return m_sampleCounts.at( radiusIndex );
}

std::size_t RadialPSFModel::nearestRadiusIndex( double requestedRadius ) const
{
    if( !std::isfinite( requestedRadius ) || requestedRadius < 0 )
    {
        throw std::invalid_argument( "radial PSF response radius must be finite and nonnegative" );
    }
    const auto upper = std::lower_bound( m_radii.begin(), m_radii.end(), requestedRadius );
    if( upper == m_radii.begin() )
    {
        return 0;
    }
    if( upper == m_radii.end() )
    {
        return m_radii.size() - 1;
    }
    const std::size_t upperIndex = static_cast<std::size_t>( upper - m_radii.begin() );
    const std::size_t lowerIndex = upperIndex - 1;
    return requestedRadius - m_radii[lowerIndex] <= m_radii[upperIndex] - requestedRadius ? lowerIndex : upperIndex;
}

} // namespace improc
} // namespace mx
