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
    : RadialPSFModel( std::move( radii ), {}, stampRows, stampColumns )
{
}

RadialPSFModel::RadialPSFModel( std::vector<double> radii,
                                std::vector<std::size_t> radiusRegions,
                                int stampRows,
                                int stampColumns )
    : m_radii( std::move( radii ) ), m_radiusRegions( std::move( radiusRegions ) ), m_stampRows( stampRows ),
      m_stampColumns( stampColumns )
{
    validateRadialPSFRadii( m_radii );
    if( m_radiusRegions.empty() )
    {
        m_radiusRegions.assign( m_radii.size(), 0 );
    }
    if( m_radiusRegions.size() != m_radii.size() || !std::is_sorted( m_radiusRegions.begin(), m_radiusRegions.end() ) )
    {
        throw std::invalid_argument( "radial PSF radius regions must be nondecreasing and correspond to radii" );
    }
    if( stampRows <= 0 || stampColumns <= 0 )
    {
        throw std::invalid_argument( "radial PSF stamp dimensions must be positive" );
    }
    m_responses.resize( m_radii.size() );
    m_validities.resize( m_radii.size() );
    m_sampleCounts.assign( m_radii.size(), 0 );
}

std::vector<std::size_t> RadialPSFModel::angularSampleCounts( const std::vector<double> &radii,
                                                              std::size_t samplesPerRadius,
                                                              double sampleArcStep )
{
    validateRadialPSFRadii( radii );
    if( !std::isfinite( sampleArcStep ) || sampleArcStep < 0 || ( samplesPerRadius == 0 ) == ( sampleArcStep == 0 ) )
    {
        throw std::invalid_argument( "radial PSF sampling requires exactly one positive angular control" );
    }

    std::vector<std::size_t> counts( radii.size(), samplesPerRadius );
    if( samplesPerRadius != 0 )
    {
        return counts;
    }

    for( std::size_t radiusIndex = 0; radiusIndex < radii.size(); ++radiusIndex )
    {
        const long double circumference = 2.0L * std::numbers::pi_v<long double> * radii[radiusIndex];
        const long double count = std::max( 1.0L, std::ceil( circumference / sampleArcStep ) );
        if( count > static_cast<long double>( std::numeric_limits<std::size_t>::max() ) )
        {
            throw std::length_error( "radial PSF angular sample count overflows size_t" );
        }
        counts[radiusIndex] = static_cast<std::size_t>( count );
    }
    return counts;
}

double RadialPSFModel::detectorPixelRadiusTolerance() noexcept
{
    return std::sqrt( 0.5 ) + 8.0 * std::numeric_limits<double>::epsilon();
}

std::vector<RadialPSFSample> RadialPSFModel::selectSamples( const std::vector<RadialPSFSource> &sources,
                                                            double centerRow,
                                                            double centerColumn,
                                                            const std::vector<double> &radii,
                                                            std::size_t samplesPerRadius )
{
    return selectSamples( sources, centerRow, centerColumn, radii, angularSampleCounts( radii, samplesPerRadius, 0 ) );
}

std::vector<RadialPSFSample> RadialPSFModel::selectSamples( const std::vector<RadialPSFSource> &sources,
                                                            double centerRow,
                                                            double centerColumn,
                                                            const std::vector<double> &radii,
                                                            const std::vector<std::size_t> &samplesPerRadius )
{
    std::vector<RadialPSFSource> globalSources = sources;
    for( RadialPSFSource &source : globalSources )
    {
        source.regionIndex = 0;
    }
    return selectSamples( globalSources,
                          centerRow,
                          centerColumn,
                          radii,
                          samplesPerRadius,
                          std::vector<std::size_t>( radii.size(), 0 ) );
}

std::vector<RadialPSFSample> RadialPSFModel::selectSamples( const std::vector<RadialPSFSource> &sources,
                                                            double centerRow,
                                                            double centerColumn,
                                                            const std::vector<double> &radii,
                                                            const std::vector<std::size_t> &samplesPerRadius,
                                                            const std::vector<std::size_t> &radiusRegions )
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
    if( samplesPerRadius.size() != radii.size() || std::any_of( samplesPerRadius.begin(),
                                                                samplesPerRadius.end(),
                                                                []( std::size_t count ) { return count == 0; } ) )
    {
        throw std::invalid_argument( "radial PSF radii and positive angular sample counts must correspond" );
    }
    if( radiusRegions.size() != radii.size() || !std::is_sorted( radiusRegions.begin(), radiusRegions.end() ) )
    {
        throw std::invalid_argument( "radial PSF radius regions must be nondecreasing and correspond to radii" );
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
    std::size_t requestedCount{ 0 };
    for( const std::size_t count : samplesPerRadius )
    {
        if( requestedCount > std::numeric_limits<std::size_t>::max() - count )
        {
            throw std::length_error( "radial PSF sample count overflows size_t" );
        }
        requestedCount += count;
    }
    samples.reserve( requestedCount );
    const double radialTolerance = detectorPixelRadiusTolerance();
    for( std::size_t radiusIndex = 0; radiusIndex < radii.size(); ++radiusIndex )
    {
        std::vector<const RadialPSFSource *> radialSources;
        radialSources.reserve( sources.size() );
        for( const RadialPSFSource &source : sources )
        {
            const double sourceRadius = std::hypot( source.row - centerRow, source.column - centerColumn );
            if( source.regionIndex == radiusRegions[radiusIndex] &&
                std::abs( sourceRadius - radii[radiusIndex] ) <= radialTolerance )
            {
                radialSources.push_back( &source );
            }
        }
        if( radialSources.empty() )
        {
            throw std::runtime_error( "radial PSF sampling found no detector pixel near a requested radius" );
        }

        std::unordered_set<std::size_t> selectedAtRadius;
        for( std::size_t angularIndex = 0; angularIndex < samplesPerRadius[radiusIndex]; ++angularIndex )
        {
            const double requestedAngle = 2.0 * std::numbers::pi * static_cast<double>( angularIndex ) /
                                          static_cast<double>( samplesPerRadius[radiusIndex] );
            const double requestedRow = centerRow + radii[radiusIndex] * std::sin( requestedAngle );
            const double requestedColumn = centerColumn + radii[radiusIndex] * std::cos( requestedAngle );
            const RadialPSFSource *nearest = nullptr;
            double nearestDistanceSquared{ 0 };
            for( const RadialPSFSource *source : radialSources )
            {
                const double deltaRow = source->row - requestedRow;
                const double deltaColumn = source->column - requestedColumn;
                const double distanceSquared = deltaRow * deltaRow + deltaColumn * deltaColumn;
                if( nearest == nullptr || distanceSquared < nearestDistanceSquared ||
                    ( distanceSquared == nearestDistanceSquared && source->sourceIndex < nearest->sourceIndex ) )
                {
                    nearest = source;
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
                                     std::atan2( deltaRow, deltaColumn ),
                                     radiusRegions[radiusIndex] } );
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
        if( sample.radiusIndex >= m_radii.size() || sample.regionIndex != m_radiusRegions[sample.radiusIndex] ||
            measured.rows() != m_stampRows || measured.cols() != m_stampColumns ||
            measuredValidity.rows() != m_stampRows || measuredValidity.cols() != m_stampColumns ||
            !std::isfinite( sample.radius ) || sample.radius < 0 || !std::isfinite( sample.angle ) ||
            !validRadialPSFValues( measured, measuredValidity ) )
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
    if( m_radiusRegions.front() != m_radiusRegions.back() )
    {
        throw std::logic_error( "region-aware radial PSF response requires an explicit source region" );
    }
    response( output, outputValidity, requestedRadius, angle, 0, m_radii.size() );
}

void RadialPSFModel::response(
    imageT &output, validityT &outputValidity, double requestedRadius, double angle, std::size_t regionIndex ) const
{
    const auto begin = std::lower_bound( m_radiusRegions.begin(), m_radiusRegions.end(), regionIndex );
    const auto end = std::upper_bound( begin, m_radiusRegions.end(), regionIndex );
    if( begin == end )
    {
        throw std::invalid_argument( "radial PSF response region has no configured radial nodes" );
    }
    response( output,
              outputValidity,
              requestedRadius,
              angle,
              static_cast<std::size_t>( begin - m_radiusRegions.begin() ),
              static_cast<std::size_t>( end - m_radiusRegions.begin() ) );
}

void RadialPSFModel::response( imageT &output,
                               validityT &outputValidity,
                               double requestedRadius,
                               double angle,
                               std::size_t beginIndex,
                               std::size_t endIndex ) const
{
    if( !std::isfinite( requestedRadius ) || requestedRadius < 0 || !std::isfinite( angle ) )
    {
        throw std::invalid_argument( "radial PSF response radius and angle must be finite and nonnegative" );
    }
    if( beginIndex >= endIndex || endIndex > m_radii.size() )
    {
        throw std::invalid_argument( "radial PSF response requires a nonempty valid radial-node range" );
    }
    const auto rangeBegin = m_radii.begin() + static_cast<std::ptrdiff_t>( beginIndex );
    const auto rangeEnd = m_radii.begin() + static_cast<std::ptrdiff_t>( endIndex );
    const auto upper = std::lower_bound( rangeBegin, rangeEnd, requestedRadius );
    std::size_t lowerIndex{ 0 };
    std::size_t upperIndex{ 0 };
    double upperFraction{ 0 };
    if( upper == rangeEnd )
    {
        lowerIndex = upperIndex = endIndex - 1;
    }
    else if( upper == rangeBegin || *upper == requestedRadius )
    {
        lowerIndex = upperIndex = static_cast<std::size_t>( upper - m_radii.begin() );
    }
    else
    {
        upperIndex = static_cast<std::size_t>( upper - m_radii.begin() );
        lowerIndex = upperIndex - 1;
        upperFraction = ( requestedRadius - m_radii[lowerIndex] ) / ( m_radii[upperIndex] - m_radii[lowerIndex] );
    }
    if( m_responses[lowerIndex].size() == 0 || m_responses[upperIndex].size() == 0 )
    {
        throw std::logic_error( "radial PSF model must be fitted before evaluating a response" );
    }

    if( lowerIndex == upperIndex )
    {
        rotate( output, outputValidity, m_responses[lowerIndex], m_validities[lowerIndex], -angle );
        return;
    }

    imageT interpolated = imageT::Zero( m_stampRows, m_stampColumns );
    validityT interpolatedValidity = validityT::Zero( m_stampRows, m_stampColumns );
    const double lowerFraction = 1.0 - upperFraction;
    for( int column = 0; column < m_stampColumns; ++column )
    {
        for( int row = 0; row < m_stampRows; ++row )
        {
            if( m_validities[lowerIndex]( row, column ) == 0 || m_validities[upperIndex]( row, column ) == 0 )
            {
                continue;
            }
            const double value = lowerFraction * m_responses[lowerIndex]( row, column ) +
                                 upperFraction * m_responses[upperIndex]( row, column );
            if( std::isfinite( value ) )
            {
                interpolated( row, column ) = static_cast<float>( value );
                interpolatedValidity( row, column ) = 1;
            }
        }
    }
    rotate( output, outputValidity, interpolated, interpolatedValidity, -angle );
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

} // namespace improc
} // namespace mx
