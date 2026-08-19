/** \file P4PixelGrid.cpp
 * \brief Implements the local pixel geometry used by Pixel Prediction Post-Processing.
 * \author Jared R. Males
 */

#include "P4PixelGrid.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <utility>

#include <mx/improc/imageMasks.hpp>
#include <mx/math/floatUtils.hpp>
#include <mx/math/geo.hpp>

namespace mx
{
namespace improc
{

namespace
{

/// Determine whether a required geometry scalar is finite under fast-math.
bool p4PixelGridFinite( double value /**< [in] scalar to validate */ )
{
    return mx::math::isFinite( value );
}

/// Validate that an exclusion policy enumerator is supported.
bool p4PixelGridValidPolicy( P4ExclusionPolicy policy /**< [in] policy to validate */ )
{
    return policy == P4ExclusionPolicy::sampleCenter || policy == P4ExclusionPolicy::kernelSupport;
}

/// Combine independently detected target- and predictor-mask states.
P4PixelInvalidReason p4PixelGridInvalidReason( bool targetMasked, /**< [in] whether the target is masked */
                                               bool predictorMasked /**< [in] whether a predictor touches the mask */ )
{
    if( targetMasked && predictorMasked )
    {
        return P4PixelInvalidReason::targetAndPredictorMasked;
    }

    if( targetMasked )
    {
        return P4PixelInvalidReason::targetMasked;
    }

    if( predictorMasked )
    {
        return P4PixelInvalidReason::predictorMasked;
    }

    return P4PixelInvalidReason::none;
}

/// Return whether a nominal interpolation footprint overlaps a disk about a target pixel.
template <typename gridT>
bool p4PixelGridFootprintOverlapsDisk( int footprintRow,                /**< [in] first nominal-footprint row */
                                       int footprintColumn,             /**< [in] first nominal-footprint column */
                                       const P4PixelCoordinate &target, /**< [in] exclusion-disk center */
                                       double effectiveRadius /**< [in] disk radius */ )
{
    const double radiusSquared = effectiveRadius * effectiveRadius;
    for( int rowOffset = 0; rowOffset < gridT::width; ++rowOffset )
    {
        const double deltaRow = static_cast<double>( footprintRow ) + static_cast<double>( rowOffset ) -
                                static_cast<double>( target.row() );
        for( int columnOffset = 0; columnOffset < gridT::width; ++columnOffset )
        {
            const double deltaColumn = static_cast<double>( footprintColumn ) + static_cast<double>( columnOffset ) -
                                       static_cast<double>( target.column() );
            if( deltaRow * deltaRow + deltaColumn * deltaColumn <= radiusSquared )
            {
                return true;
            }
        }
    }

    return false;
}

/// Map one canonical predictor offset to one search pixel using direct scalar rotation.
template <typename gridT>
std::pair<double, double> p4PixelGridMap( const P4PixelCoordinate &searchPixel, /**< [in] target coordinate */
                                          const P4PixelCoordinate &offset,      /**< [in] canonical predictor offset */
                                          double xCenter,                       /**< [in] image-center row */
                                          double yCenter /**< [in] image-center column */ )
{
    const double deltaRow = static_cast<double>( searchPixel.row() ) - xCenter;
    const double deltaColumn = static_cast<double>( searchPixel.column() ) - yCenter;
    const double rotation = std::atan2( deltaColumn, deltaRow ) - 0.5 * std::numbers::pi_v<double>;
    const double cosine = std::cos( rotation );
    const double sine = std::sin( rotation );
    const double mappedRow = static_cast<double>( searchPixel.row() ) + static_cast<double>( offset.row() ) * cosine -
                             static_cast<double>( offset.column() ) * sine;
    const double mappedColumn = static_cast<double>( searchPixel.column() ) +
                                static_cast<double>( offset.row() ) * sine +
                                static_cast<double>( offset.column() ) * cosine;
    return { mappedRow, mappedColumn };
}

} // namespace

P4PixelGridRegion::P4PixelGridRegion( double searchInnerRadius,
                                      double searchOuterRadius,
                                      double optimizationDeltaRadiusInner,
                                      double optimizationDeltaRadiusOuter,
                                      double optimizationArcHalfWidth,
                                      double optimizationMaxHalfAngle,
                                      double psfRadius,
                                      P4ExclusionPolicy exclusionPolicy,
                                      double exclusionRadiusBuffer )
    : searchInnerRadius( searchInnerRadius ), searchOuterRadius( searchOuterRadius ),
      optimizationDeltaRadiusInner( optimizationDeltaRadiusInner ),
      optimizationDeltaRadiusOuter( optimizationDeltaRadiusOuter ),
      optimizationArcHalfWidth( optimizationArcHalfWidth ), optimizationMaxHalfAngle( optimizationMaxHalfAngle ),
      psfRadius( psfRadius ), exclusionPolicy( exclusionPolicy ), exclusionRadiusBuffer( exclusionRadiusBuffer )
{
}

P4PixelCoordinate::P4PixelCoordinate( int row, int column ) : m_row( row ), m_column( column )
{
}

int P4PixelCoordinate::row() const noexcept
{
    return m_row;
}

int P4PixelCoordinate::column() const noexcept
{
    return m_column;
}

P4SearchPixelRecord::P4SearchPixelRecord( const P4PixelCoordinate &coordinate, P4PixelInvalidReason invalidReason )
    : m_coordinate( coordinate ), m_invalidReason( invalidReason )
{
}

const P4PixelCoordinate &P4SearchPixelRecord::coordinate() const noexcept
{
    return m_coordinate;
}

P4PixelInvalidReason P4SearchPixelRecord::invalidReason() const noexcept
{
    return m_invalidReason;
}

bool P4SearchPixelRecord::valid() const noexcept
{
    return m_invalidReason == P4PixelInvalidReason::none;
}

bool P4SearchPixelRecord::targetMasked() const noexcept
{
    return m_invalidReason == P4PixelInvalidReason::targetMasked ||
           m_invalidReason == P4PixelInvalidReason::targetAndPredictorMasked;
}

bool P4SearchPixelRecord::predictorMasked() const noexcept
{
    return m_invalidReason == P4PixelInvalidReason::predictorMasked ||
           m_invalidReason == P4PixelInvalidReason::targetAndPredictorMasked;
}

template <typename transformT>
P4InterpolationRecord<transformT>::P4InterpolationRecord(
    double row, double column, int footprintRow, int footprintColumn, const kernelT &kernel )
    : m_row( row ), m_column( column ), m_footprintRow( footprintRow ), m_footprintColumn( footprintColumn ),
      m_kernel( kernel )
{
}

template <typename transformT>
double P4InterpolationRecord<transformT>::row() const noexcept
{
    return m_row;
}

template <typename transformT>
double P4InterpolationRecord<transformT>::column() const noexcept
{
    return m_column;
}

template <typename transformT>
int P4InterpolationRecord<transformT>::footprintRow() const noexcept
{
    return m_footprintRow;
}

template <typename transformT>
int P4InterpolationRecord<transformT>::footprintColumn() const noexcept
{
    return m_footprintColumn;
}

template <typename transformT>
const typename P4InterpolationRecord<transformT>::kernelT &P4InterpolationRecord<transformT>::kernel() const noexcept
{
    return m_kernel;
}

template <typename transformT>
P4PixelGrid<transformT>::P4PixelGrid() = default;

template <typename transformT>
int P4PixelGrid<transformT>::checkedFootprintOrigin( double mappedCoordinate )
{
    if( !p4PixelGridFinite( mappedCoordinate ) )
    {
        throw std::invalid_argument( "P4PixelGrid mapped interpolation coordinate must be finite" );
    }

    const double footprintOrigin = std::floor( mappedCoordinate ) - static_cast<double>( leftBuffer );
    if( footprintOrigin < static_cast<double>( std::numeric_limits<int>::min() ) ||
        footprintOrigin > static_cast<double>( std::numeric_limits<int>::max() ) )
    {
        throw std::invalid_argument( "P4PixelGrid mapped interpolation coordinate exceeds integer range" );
    }

    return static_cast<int>( footprintOrigin );
}

template <typename transformT>
std::size_t P4PixelGrid<transformT>::checkedInterpolationCount( std::size_t searchPixelCount,
                                                                std::size_t predictorCount )
{
    if( searchPixelCount != 0 && predictorCount > std::numeric_limits<std::size_t>::max() / searchPixelCount )
    {
        throw std::length_error( "P4PixelGrid interpolation-record count exceeds size_t range" );
    }

    return searchPixelCount * predictorCount;
}

template <typename transformT>
void P4PixelGrid<transformT>::resize( int rows, int columns )
{
    resize( rows, columns, 0.5 * static_cast<double>( rows - 1 ), 0.5 * static_cast<double>( columns - 1 ) );
}

template <typename transformT>
void P4PixelGrid<transformT>::resize( int rows, int columns, double xCenter, double yCenter )
{
    if( rows <= 0 || columns <= 0 )
    {
        throw std::invalid_argument( "P4PixelGrid image dimensions must be positive" );
    }

    if( !p4PixelGridFinite( xCenter ) || !p4PixelGridFinite( yCenter ) || xCenter < 0 || yCenter < 0 ||
        xCenter > static_cast<double>( rows - 1 ) || yCenter > static_cast<double>( columns - 1 ) )
    {
        throw std::invalid_argument( "P4PixelGrid center must be finite and inside the image" );
    }

    geometryImageT radiusImage( rows, columns );
    geometryImageT angleImage;
    mx::improc::radAngImage<mx::math::degreesT<double>>( radiusImage, angleImage, xCenter, yCenter );

    m_rows = rows;
    m_columns = columns;
    m_xCenter = xCenter;
    m_yCenter = yCenter;
    m_radiusImage = std::move( radiusImage );
    m_angleImage = std::move( angleImage );
    m_region.reset();
    m_regionComplete = false;
    m_searchPixels.clear();
    m_candidateOffsets.clear();
    m_predictorOffsets.clear();
    m_interpolations.clear();
}

template <typename transformT>
void P4PixelGrid<transformT>::region( const P4PixelGridRegion &configuration, const imageT *commonMask )
{
    regionImpl( configuration, commonMask, false );
}

template <typename transformT>
void P4PixelGrid<transformT>::candidateRegion( const P4PixelGridRegion &configuration )
{
    regionImpl( configuration, nullptr, true );
}

template <typename transformT>
void P4PixelGrid<transformT>::regionImpl( const P4PixelGridRegion &configuration,
                                          const imageT *commonMask,
                                          bool candidatesOnly )
{
    if( !resized() )
    {
        throw std::logic_error( "P4PixelGrid must be resized before configuring a region" );
    }

    const double scalars[]{ configuration.searchInnerRadius,
                            configuration.searchOuterRadius,
                            configuration.optimizationDeltaRadiusInner,
                            configuration.optimizationDeltaRadiusOuter,
                            configuration.optimizationArcHalfWidth,
                            configuration.optimizationMaxHalfAngle,
                            configuration.psfRadius,
                            configuration.exclusionRadiusBuffer };
    if( !std::all_of( std::begin( scalars ), std::end( scalars ), p4PixelGridFinite ) )
    {
        throw std::invalid_argument( "P4PixelGrid region values must be finite" );
    }

    if( configuration.searchInnerRadius < 0 || configuration.searchOuterRadius <= configuration.searchInnerRadius )
    {
        throw std::invalid_argument( "P4PixelGrid search radii must define a nonempty annulus" );
    }

    if( configuration.optimizationDeltaRadiusInner <= 0 || configuration.optimizationDeltaRadiusOuter <= 0 ||
        configuration.optimizationArcHalfWidth < 0 || configuration.optimizationMaxHalfAngle <= 0 ||
        configuration.optimizationMaxHalfAngle > 180 || configuration.psfRadius <= 0 ||
        configuration.exclusionRadiusBuffer < 0 )
    {
        throw std::invalid_argument( "P4PixelGrid optimization and exclusion values are outside their valid ranges" );
    }

    if( !p4PixelGridValidPolicy( configuration.exclusionPolicy ) )
    {
        throw std::invalid_argument( "P4PixelGrid exclusion policy is unsupported" );
    }

    if( commonMask && ( commonMask->rows() != m_rows || commonMask->cols() != m_columns ) )
    {
        throw std::invalid_argument( "P4PixelGrid common mask dimensions must match the image" );
    }

    if( commonMask )
    {
        for( int column = 0; column < commonMask->cols(); ++column )
        {
            for( int row = 0; row < commonMask->rows(); ++row )
            {
                if( !mx::math::isFinite( ( *commonMask )( row, column ) ) )
                {
                    throw std::invalid_argument( "P4PixelGrid common mask must contain only finite values" );
                }
            }
        }
    }

    const double regionMiddleRadius =
        configuration.searchInnerRadius + 0.5 * ( configuration.searchOuterRadius - configuration.searchInnerRadius );
    const double effectiveRadius = configuration.psfRadius + configuration.exclusionRadiusBuffer;
    if( !p4PixelGridFinite( regionMiddleRadius ) || regionMiddleRadius <= 0 || !p4PixelGridFinite( effectiveRadius ) )
    {
        throw std::invalid_argument( "P4PixelGrid derived region values must be finite" );
    }

    geometryImageT *noGeometryMask{ nullptr };
    const auto searchCoordinates =
        mx::improc::annulusCoords<mx::math::degreesT<double>>( m_radiusImage,
                                                               m_angleImage,
                                                               m_xCenter,
                                                               m_yCenter,
                                                               configuration.searchInnerRadius,
                                                               configuration.searchOuterRadius,
                                                               0.0,
                                                               360.0,
                                                               noGeometryMask,
                                                               0.0 );
    if( searchCoordinates.empty() )
    {
        throw std::invalid_argument( "P4PixelGrid search annulus contains no image pixels" );
    }

    std::vector<P4PixelCoordinate> searchPixels;
    searchPixels.reserve( searchCoordinates.size() );
    for( const auto &coordinate : searchCoordinates )
    {
        searchPixels.emplace_back( coordinate[0], coordinate[1] );
    }

    const double referenceRow = m_xCenter + 0.5;
    const double referenceColumn = m_yCenter + 0.5;
    const double canonicalRow = std::floor( referenceRow );
    const double canonicalColumn = std::floor( referenceColumn + regionMiddleRadius );
    const double canonicalRadius = canonicalColumn - referenceColumn;
    const double canonicalAngle = std::atan2( canonicalColumn - referenceColumn, canonicalRow - referenceRow ) * 180.0 /
                                  std::numbers::pi_v<double>;
    const double optimizationInnerRadius = canonicalRadius - configuration.optimizationDeltaRadiusInner;
    const double optimizationOuterRadius = canonicalRadius + configuration.optimizationDeltaRadiusOuter;
    if( optimizationOuterRadius + 0.5 <= 0 )
    {
        throw std::invalid_argument( "P4PixelGrid optimization annulus contains no possible pixels" );
    }

    double optimizationHalfAngle = configuration.optimizationMaxHalfAngle;
    if( configuration.optimizationArcHalfWidth > 0 )
    {
        optimizationHalfAngle =
            configuration.optimizationArcHalfWidth / ( 2.0 * std::numbers::pi_v<double> * regionMiddleRadius ) * 360.0;
        optimizationHalfAngle = std::min( optimizationHalfAngle, configuration.optimizationMaxHalfAngle );
    }
    const double angularTolerance = 16.0 * std::numeric_limits<double>::epsilon() *
                                    std::max( { 1.0, std::abs( canonicalAngle ), optimizationHalfAngle } );

    const double maximumBufferedRadius = optimizationOuterRadius + 0.5;
    const double canonicalMarginValue = std::ceil( std::abs( canonicalRadius ) + maximumBufferedRadius + 2.0 );
    const double maximumMargin = static_cast<double>( ( std::numeric_limits<int>::max() - 1 ) / 2 );
    if( !p4PixelGridFinite( canonicalMarginValue ) || canonicalMarginValue <= 0 ||
        canonicalMarginValue > maximumMargin )
    {
        throw std::invalid_argument( "P4PixelGrid canonical optimization geometry is too large" );
    }
    const int canonicalMargin = static_cast<int>( canonicalMarginValue );

    const double maximumUsefulRadius =
        std::hypot( static_cast<double>( m_rows ), static_cast<double>( m_columns ) ) + static_cast<double>( width );
    if( maximumBufferedRadius > maximumUsefulRadius )
    {
        throw std::invalid_argument( "P4PixelGrid optimization geometry cannot fit inside the image" );
    }

    const int canonicalSize = 2 * canonicalMargin + 1;
    const double localReferenceRow = static_cast<double>( canonicalMargin ) - ( canonicalRow - referenceRow );
    const double localReferenceColumn = static_cast<double>( canonicalMargin ) - ( canonicalColumn - referenceColumn );
    geometryImageT optimizationRadiusImage( canonicalSize, canonicalSize );
    geometryImageT optimizationAngleImage;
    mx::improc::radAngImage<mx::math::degreesT<double>>( optimizationRadiusImage,
                                                         optimizationAngleImage,
                                                         localReferenceRow,
                                                         localReferenceColumn );

    const auto optimizationCoordinates = mx::improc::annulusCoords<mx::math::degreesT<double>>( optimizationRadiusImage,
                                                                                                optimizationAngleImage,
                                                                                                localReferenceRow,
                                                                                                localReferenceColumn,
                                                                                                optimizationInnerRadius,
                                                                                                optimizationOuterRadius,
                                                                                                0.0,
                                                                                                360.0,
                                                                                                noGeometryMask,
                                                                                                0.5 );

    std::vector<P4PixelCoordinate> candidateOffsets;
    candidateOffsets.reserve( optimizationCoordinates.size() );
    for( const auto &coordinate : optimizationCoordinates )
    {
        const double angularDifference =
            mx::math::angleDiff<mx::math::degreesT<double>>( canonicalAngle,
                                                             optimizationAngleImage( coordinate[0], coordinate[1] ) );
        if( std::abs( angularDifference ) > optimizationHalfAngle + angularTolerance )
        {
            continue;
        }

        const P4PixelCoordinate offset( coordinate[0] - canonicalMargin, coordinate[1] - canonicalMargin );
        candidateOffsets.push_back( offset );
    }

    if( candidateOffsets.empty() )
    {
        throw std::invalid_argument( "P4PixelGrid optimization wedge contains no predictor candidates" );
    }

    if( candidatesOnly )
    {
        std::vector<P4SearchPixelRecord> candidateSearchRecords;
        candidateSearchRecords.reserve( searchPixels.size() );
        for( const P4PixelCoordinate &searchPixel : searchPixels )
        {
            candidateSearchRecords.emplace_back( searchPixel, P4PixelInvalidReason::none );
        }

        m_region = configuration;
        m_regionComplete = false;
        m_searchPixels = std::move( candidateSearchRecords );
        m_candidateOffsets = std::move( candidateOffsets );
        m_predictorOffsets.clear();
        m_interpolations.clear();
        return;
    }

    std::vector<P4PixelCoordinate> predictorOffsets;
    predictorOffsets.reserve( candidateOffsets.size() );
    for( const P4PixelCoordinate &offset : candidateOffsets )
    {
        if( configuration.exclusionPolicy == P4ExclusionPolicy::sampleCenter &&
            std::hypot( static_cast<double>( offset.row() ), static_cast<double>( offset.column() ) ) <=
                effectiveRadius )
        {
            continue;
        }

        bool excluded{ false };
        if( configuration.exclusionPolicy == P4ExclusionPolicy::kernelSupport )
        {
            for( const P4PixelCoordinate &searchPixel : searchPixels )
            {
                const auto mapped =
                    p4PixelGridMap<P4PixelGrid<transformT>>( searchPixel, offset, m_xCenter, m_yCenter );
                const int footprintRow = checkedFootprintOrigin( mapped.first );
                const int footprintColumn = checkedFootprintOrigin( mapped.second );
                if( p4PixelGridFootprintOverlapsDisk<P4PixelGrid<transformT>>( footprintRow,
                                                                               footprintColumn,
                                                                               searchPixel,
                                                                               effectiveRadius ) )
                {
                    excluded = true;
                    break;
                }
            }
        }

        if( !excluded )
        {
            predictorOffsets.push_back( offset );
        }
    }

    if( predictorOffsets.empty() )
    {
        throw std::invalid_argument( "P4PixelGrid support exclusion removed every predictor column" );
    }

    transformT transform;
    std::vector<P4SearchPixelRecord> searchRecords;
    searchRecords.reserve( searchPixels.size() );
    std::vector<interpolationRecordT> interpolations;
    interpolations.reserve( checkedInterpolationCount( searchPixels.size(), predictorOffsets.size() ) );

    for( const P4PixelCoordinate &searchPixel : searchPixels )
    {
        const bool targetMasked = commonMask && ( *commonMask )( searchPixel.row(), searchPixel.column() ) == 0;
        bool predictorMasked{ false };

        for( const P4PixelCoordinate &offset : predictorOffsets )
        {
            const auto mapped = p4PixelGridMap<P4PixelGrid<transformT>>( searchPixel, offset, m_xCenter, m_yCenter );
            const int footprintRow = checkedFootprintOrigin( mapped.first );
            const int footprintColumn = checkedFootprintOrigin( mapped.second );

            if( footprintRow < 0 || footprintColumn < 0 || static_cast<std::int64_t>( footprintRow ) + width > m_rows ||
                static_cast<std::int64_t>( footprintColumn ) + width > m_columns )
            {
                throw std::invalid_argument( "P4PixelGrid mapped interpolation footprint crosses an image edge" );
            }

            if( commonMask )
            {
                for( int rowOffset = 0; rowOffset < width; ++rowOffset )
                {
                    for( int columnOffset = 0; columnOffset < width; ++columnOffset )
                    {
                        if( ( *commonMask )( footprintRow + rowOffset, footprintColumn + columnOffset ) == 0 )
                        {
                            predictorMasked = true;
                        }
                    }
                }
            }

            kernelT kernel;
            transform( kernel,
                       static_cast<realT>( mapped.first - std::floor( mapped.first ) ),
                       static_cast<realT>( mapped.second - std::floor( mapped.second ) ) );
            interpolations.push_back(
                interpolationRecordT( mapped.first, mapped.second, footprintRow, footprintColumn, kernel ) );
        }

        searchRecords.emplace_back( searchPixel, p4PixelGridInvalidReason( targetMasked, predictorMasked ) );
    }

    m_region = configuration;
    m_regionComplete = true;
    m_searchPixels = std::move( searchRecords );
    m_candidateOffsets = std::move( candidateOffsets );
    m_predictorOffsets = std::move( predictorOffsets );
    m_interpolations = std::move( interpolations );
}

template <typename transformT>
int P4PixelGrid<transformT>::rows() const noexcept
{
    return m_rows;
}

template <typename transformT>
int P4PixelGrid<transformT>::columns() const noexcept
{
    return m_columns;
}

template <typename transformT>
double P4PixelGrid<transformT>::xCenter() const noexcept
{
    return m_xCenter;
}

template <typename transformT>
double P4PixelGrid<transformT>::yCenter() const noexcept
{
    return m_yCenter;
}

template <typename transformT>
bool P4PixelGrid<transformT>::resized() const noexcept
{
    return m_rows > 0 && m_columns > 0;
}

template <typename transformT>
bool P4PixelGrid<transformT>::regionConfigured() const noexcept
{
    return m_region.has_value() && m_regionComplete;
}

template <typename transformT>
bool P4PixelGrid<transformT>::candidateRegionConfigured() const noexcept
{
    return m_region.has_value() && !m_searchPixels.empty() && !m_candidateOffsets.empty();
}

template <typename transformT>
const P4PixelGridRegion &P4PixelGrid<transformT>::regionConfiguration() const
{
    if( !m_region )
    {
        throw std::logic_error( "P4PixelGrid has no configured region" );
    }
    return *m_region;
}

template <typename transformT>
double P4PixelGrid<transformT>::effectiveExclusionRadius() const
{
    const P4PixelGridRegion &configuration = regionConfiguration();
    return configuration.psfRadius + configuration.exclusionRadiusBuffer;
}

template <typename transformT>
std::size_t P4PixelGrid<transformT>::searchPixelCount() const noexcept
{
    return m_searchPixels.size();
}

template <typename transformT>
std::size_t P4PixelGrid<transformT>::predictorCount() const noexcept
{
    return m_predictorOffsets.size();
}

template <typename transformT>
std::size_t P4PixelGrid<transformT>::candidatePredictorCount() const noexcept
{
    return m_candidateOffsets.size();
}

template <typename transformT>
const P4SearchPixelRecord &P4PixelGrid<transformT>::searchPixel( std::size_t searchIndex ) const
{
    if( searchIndex >= m_searchPixels.size() )
    {
        throw std::out_of_range( "P4PixelGrid search-pixel index is out of range" );
    }
    return m_searchPixels[searchIndex];
}

template <typename transformT>
const P4PixelCoordinate &P4PixelGrid<transformT>::predictorOffset( std::size_t predictorIndex ) const
{
    if( predictorIndex >= m_predictorOffsets.size() )
    {
        throw std::out_of_range( "P4PixelGrid predictor index is out of range" );
    }
    return m_predictorOffsets[predictorIndex];
}

template <typename transformT>
const P4PixelCoordinate &P4PixelGrid<transformT>::candidatePredictorOffset( std::size_t candidateIndex ) const
{
    if( candidateIndex >= m_candidateOffsets.size() )
    {
        throw std::out_of_range( "P4PixelGrid candidate-predictor index is out of range" );
    }
    return m_candidateOffsets[candidateIndex];
}

template <typename transformT>
std::pair<double, double> P4PixelGrid<transformT>::candidateCoordinate( std::size_t searchIndex,
                                                                        std::size_t candidateIndex ) const
{
    const P4SearchPixelRecord &search = searchPixel( searchIndex );
    const P4PixelCoordinate &candidate = candidatePredictorOffset( candidateIndex );
    return p4PixelGridMap<P4PixelGrid<transformT>>( search.coordinate(), candidate, m_xCenter, m_yCenter );
}

template <typename transformT>
const typename P4PixelGrid<transformT>::interpolationRecordT &
P4PixelGrid<transformT>::interpolation( std::size_t searchIndex, std::size_t predictorIndex ) const
{
    if( searchIndex >= m_searchPixels.size() || predictorIndex >= m_predictorOffsets.size() )
    {
        throw std::out_of_range( "P4PixelGrid interpolation index is out of range" );
    }
    return m_interpolations[searchIndex * m_predictorOffsets.size() + predictorIndex];
}

template <typename transformT>
typename P4PixelGrid<transformT>::realT
P4PixelGrid<transformT>::sample( const imageT &image, std::size_t searchIndex, std::size_t predictorIndex ) const
{
    if( image.rows() != m_rows || image.cols() != m_columns )
    {
        throw std::invalid_argument( "P4PixelGrid sampled image dimensions must match the grid" );
    }

    if( !searchPixel( searchIndex ).valid() )
    {
        throw std::logic_error( "P4PixelGrid cannot sample a common-mask-invalid local fit" );
    }

    const interpolationRecordT &record = interpolation( searchIndex, predictorIndex );
    realT value{ 0 };
    for( int rowOffset = 0; rowOffset < width; ++rowOffset )
    {
        for( int columnOffset = 0; columnOffset < width; ++columnOffset )
        {
            value += image( record.footprintRow() + rowOffset, record.footprintColumn() + columnOffset ) *
                     record.kernel()( rowOffset, columnOffset );
        }
    }
    return value;
}

template class P4InterpolationRecord<mx::improc::cubicConvolTransform<float>>;
template class P4PixelGrid<mx::improc::cubicConvolTransform<float>>;

} // namespace improc
} // namespace mx
