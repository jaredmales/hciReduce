/** \file P4RotatedGrid.cpp
 * \brief Implements direct sky-frame sampling geometry for rotated P4 regression.
 * \author Jared R. Males
 */

#include "P4RotatedGrid.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#include <mx/math/floatUtils.hpp>

namespace mx
{
namespace improc
{

namespace
{

/// Return the integer representation of a rotated-grid invalidity mask.
constexpr std::uint8_t p4RotatedReasonValue( P4RotatedInvalidReason reason /**< [in] invalidity mask */ ) noexcept
{
    return static_cast<std::uint8_t>( reason );
}

/// Add one invalidity bit to a combined rotated-grid reason.
void p4RotatedAddReason( P4RotatedInvalidReason &combined, /**< [in,out] combined invalidity mask */
                         P4RotatedInvalidReason reason /**< [in] invalidity bit to add */ ) noexcept
{
    combined = static_cast<P4RotatedInvalidReason>( p4RotatedReasonValue( combined ) | p4RotatedReasonValue( reason ) );
}

/// Report whether a scalar is finite under fast-math.
bool p4RotatedFinite( double value /**< [in] scalar to validate */ ) noexcept
{
    return mx::math::isFinite( value );
}

} // namespace

P4RotatedSearchRecord::P4RotatedSearchRecord( const P4PixelCoordinate &coordinate,
                                              P4RotatedInvalidReason invalidReason )
    : m_coordinate( coordinate ), m_invalidReason( invalidReason )
{
    constexpr std::uint8_t supportedReasons =
        p4RotatedReasonValue( P4RotatedInvalidReason::targetFootprintOutside ) |
        p4RotatedReasonValue( P4RotatedInvalidReason::predictorFootprintOutside ) |
        p4RotatedReasonValue( P4RotatedInvalidReason::targetFootprintMasked ) |
        p4RotatedReasonValue( P4RotatedInvalidReason::predictorFootprintMasked );
    if( ( p4RotatedReasonValue( invalidReason ) & ~supportedReasons ) != 0 )
    {
        throw std::invalid_argument( "P4RotatedSearchRecord invalidity mask is unsupported" );
    }
}

const P4PixelCoordinate &P4RotatedSearchRecord::coordinate() const noexcept
{
    return m_coordinate;
}

P4RotatedInvalidReason P4RotatedSearchRecord::invalidReason() const noexcept
{
    return m_invalidReason;
}

bool P4RotatedSearchRecord::valid() const noexcept
{
    return m_invalidReason == P4RotatedInvalidReason::none;
}

bool P4RotatedSearchRecord::hasReason( P4RotatedInvalidReason reason ) const noexcept
{
    return reason != P4RotatedInvalidReason::none &&
           ( p4RotatedReasonValue( m_invalidReason ) & p4RotatedReasonValue( reason ) ) ==
               p4RotatedReasonValue( reason );
}

P4RotatedInterpolationRecord::P4RotatedInterpolationRecord(
    double row, double column, int footprintRow, int footprintColumn, const kernelT &kernel, bool footprintInside )
    : m_row( row ), m_column( column ), m_footprintRow( footprintRow ), m_footprintColumn( footprintColumn ),
      m_kernel( kernel ), m_footprintInside( footprintInside )
{
}

double P4RotatedInterpolationRecord::row() const noexcept
{
    return m_row;
}

double P4RotatedInterpolationRecord::column() const noexcept
{
    return m_column;
}

int P4RotatedInterpolationRecord::footprintRow() const noexcept
{
    return m_footprintRow;
}

int P4RotatedInterpolationRecord::footprintColumn() const noexcept
{
    return m_footprintColumn;
}

const P4RotatedInterpolationRecord::kernelT &P4RotatedInterpolationRecord::kernel() const noexcept
{
    return m_kernel;
}

bool P4RotatedInterpolationRecord::footprintInside() const noexcept
{
    return m_footprintInside;
}

P4RotatedGrid::P4RotatedGrid() = default;

void P4RotatedGrid::configure( const P4PixelGridf &skyGrid,
                               const std::vector<double> &derotationAngles,
                               const imageT *rawCommonMask )
{
    if( !skyGrid.candidateRegionConfigured() )
    {
        throw std::logic_error( "P4RotatedGrid requires candidate-configured sky geometry" );
    }
    if( derotationAngles.empty() )
    {
        throw std::invalid_argument( "P4RotatedGrid requires at least one derotation angle" );
    }
    if( rawCommonMask && ( rawCommonMask->rows() != skyGrid.rows() || rawCommonMask->cols() != skyGrid.columns() ) )
    {
        throw std::invalid_argument( "P4RotatedGrid detector common-mask dimensions must match the sky grid" );
    }
    if( rawCommonMask )
    {
        for( int column = 0; column < rawCommonMask->cols(); ++column )
        {
            for( int row = 0; row < rawCommonMask->rows(); ++row )
            {
                if( !mx::math::isFinite( ( *rawCommonMask )( row, column ) ) )
                {
                    throw std::invalid_argument( "P4RotatedGrid detector common mask must contain only finite values" );
                }
            }
        }
    }

    std::vector<FrameRotation> frameRotations;
    frameRotations.reserve( derotationAngles.size() );
    for( double angle : derotationAngles )
    {
        if( !p4RotatedFinite( angle ) )
        {
            throw std::invalid_argument( "P4RotatedGrid derotation angles must be finite" );
        }
        const FrameRotation rotation{ std::cos( angle ), std::sin( angle ) };
        frameRotations.push_back( rotation );
    }

    const int rows = skyGrid.rows();
    const int columns = skyGrid.columns();
    const double xCenter = skyGrid.xCenter();
    const double yCenter = skyGrid.yCenter();
    const std::size_t searchCount = skyGrid.searchPixelCount();
    const std::size_t candidateCount = skyGrid.candidatePredictorCount();
    const P4PixelGridRegion &configuration = skyGrid.regionConfiguration();
    const double effectiveRadius = skyGrid.effectiveExclusionRadius();
    std::vector<std::size_t> retainedCandidates;
    retainedCandidates.reserve( candidateCount );
    for( std::size_t candidate = 0; candidate < candidateCount; ++candidate )
    {
        const P4PixelCoordinate &offset = skyGrid.candidatePredictorOffset( candidate );
        bool excluded = false;
        if( configuration.exclusionPolicy == P4ExclusionPolicy::sampleCenter )
        {
            excluded = std::hypot( static_cast<double>( offset.row() ), static_cast<double>( offset.column() ) ) <=
                       effectiveRadius;
        }
        else if( configuration.exclusionPolicy == P4ExclusionPolicy::kernelSupport )
        {
            for( std::size_t search = 0; search < searchCount && !excluded; ++search )
            {
                const P4PixelCoordinate &targetCoordinate = skyGrid.searchPixel( search ).coordinate();
                const SkyCoordinate targetSky{ static_cast<double>( targetCoordinate.row() ),
                                               static_cast<double>( targetCoordinate.column() ) };
                const std::pair<double, double> mappedCandidate = skyGrid.candidateCoordinate( search, candidate );
                const SkyCoordinate predictorSky{ mappedCandidate.first, mappedCandidate.second };
                for( const FrameRotation &rotation : frameRotations )
                {
                    const MappedFootprint target = mapFootprint( rows, columns, xCenter, yCenter, rotation, targetSky );
                    const MappedFootprint predictor =
                        mapFootprint( rows, columns, xCenter, yCenter, rotation, predictorSky );
                    if( footprintOverlapsExclusion( predictor, target, effectiveRadius ) )
                    {
                        excluded = true;
                        break;
                    }
                }
            }
        }
        if( !excluded )
        {
            retainedCandidates.push_back( candidate );
        }
    }
    if( retainedCandidates.empty() )
    {
        throw std::invalid_argument( "P4RotatedGrid central exclusion removed every direct predictor column" );
    }

    std::vector<P4PixelCoordinate> predictorOffsets;
    predictorOffsets.reserve( retainedCandidates.size() );
    for( std::size_t candidate : retainedCandidates )
    {
        predictorOffsets.push_back( skyGrid.candidatePredictorOffset( candidate ) );
    }

    std::vector<SkyCoordinate> skyPredictors;
    skyPredictors.reserve( checkedCoordinateCount( searchCount, retainedCandidates.size() ) );
    for( std::size_t search = 0; search < searchCount; ++search )
    {
        for( std::size_t candidate : retainedCandidates )
        {
            const std::pair<double, double> coordinate = skyGrid.candidateCoordinate( search, candidate );
            skyPredictors.push_back( SkyCoordinate{ coordinate.first, coordinate.second } );
        }
    }

    std::vector<P4RotatedSearchRecord> searchPixels;
    searchPixels.reserve( searchCount );
    for( std::size_t search = 0; search < searchCount; ++search )
    {
        const P4PixelCoordinate &targetCoordinate = skyGrid.searchPixel( search ).coordinate();
        const SkyCoordinate targetSky{ static_cast<double>( targetCoordinate.row() ),
                                       static_cast<double>( targetCoordinate.column() ) };
        P4RotatedInvalidReason invalidReason = P4RotatedInvalidReason::none;
        for( std::size_t frame = 0; frame < frameRotations.size(); ++frame )
        {
            const MappedFootprint target =
                mapFootprint( rows, columns, xCenter, yCenter, frameRotations[frame], targetSky );
            if( !target.inside )
            {
                p4RotatedAddReason( invalidReason, P4RotatedInvalidReason::targetFootprintOutside );
            }
            else if( !maskFootprintValid( rawCommonMask, target ) )
            {
                p4RotatedAddReason( invalidReason, P4RotatedInvalidReason::targetFootprintMasked );
            }

            for( std::size_t predictor = 0; predictor < retainedCandidates.size(); ++predictor )
            {
                const SkyCoordinate &predictorSky = skyPredictors[search * retainedCandidates.size() + predictor];
                const MappedFootprint record =
                    mapFootprint( rows, columns, xCenter, yCenter, frameRotations[frame], predictorSky );
                if( !record.inside )
                {
                    p4RotatedAddReason( invalidReason, P4RotatedInvalidReason::predictorFootprintOutside );
                }
                else if( !maskFootprintValid( rawCommonMask, record ) )
                {
                    p4RotatedAddReason( invalidReason, P4RotatedInvalidReason::predictorFootprintMasked );
                }
            }
        }
        searchPixels.emplace_back( targetCoordinate, invalidReason );
    }

    m_rows = rows;
    m_columns = columns;
    m_xCenter = xCenter;
    m_yCenter = yCenter;
    m_frameRotations = std::move( frameRotations );
    m_searchPixels = std::move( searchPixels );
    m_predictorOffsets = std::move( predictorOffsets );
    m_candidatePredictorIndices = std::move( retainedCandidates );
    m_skyPredictors = std::move( skyPredictors );
}

bool P4RotatedGrid::configured() const noexcept
{
    return m_rows > 0 && m_columns > 0 && !m_frameRotations.empty() && !m_searchPixels.empty() &&
           !m_predictorOffsets.empty();
}

int P4RotatedGrid::rows() const noexcept
{
    return m_rows;
}

int P4RotatedGrid::columns() const noexcept
{
    return m_columns;
}

double P4RotatedGrid::xCenter() const noexcept
{
    return m_xCenter;
}

double P4RotatedGrid::yCenter() const noexcept
{
    return m_yCenter;
}

std::size_t P4RotatedGrid::frameCount() const noexcept
{
    return m_frameRotations.size();
}

std::size_t P4RotatedGrid::searchPixelCount() const noexcept
{
    return m_searchPixels.size();
}

std::size_t P4RotatedGrid::predictorCount() const noexcept
{
    return m_predictorOffsets.size();
}

const P4RotatedSearchRecord &P4RotatedGrid::searchPixel( std::size_t searchIndex ) const
{
    if( searchIndex >= m_searchPixels.size() )
    {
        throw std::out_of_range( "P4RotatedGrid search-pixel index is out of range" );
    }
    return m_searchPixels[searchIndex];
}

const P4PixelCoordinate &P4RotatedGrid::predictorOffset( std::size_t predictorIndex ) const
{
    if( predictorIndex >= m_predictorOffsets.size() )
    {
        throw std::out_of_range( "P4RotatedGrid predictor index is out of range" );
    }
    return m_predictorOffsets[predictorIndex];
}

std::size_t P4RotatedGrid::candidatePredictorIndex( std::size_t predictorIndex ) const
{
    if( predictorIndex >= m_candidatePredictorIndices.size() )
    {
        throw std::out_of_range( "P4RotatedGrid predictor index is out of range" );
    }
    return m_candidatePredictorIndices[predictorIndex];
}

P4RotatedGrid::interpolationRecordT P4RotatedGrid::targetInterpolation( std::size_t frameIndex,
                                                                        std::size_t searchIndex ) const
{
    if( !configured() )
    {
        throw std::logic_error( "P4RotatedGrid is not configured" );
    }
    const P4PixelCoordinate &target = searchPixel( searchIndex ).coordinate();
    return interpolation(
        frameIndex,
        SkyCoordinate{ static_cast<double>( target.row() ), static_cast<double>( target.column() ) } );
}

P4RotatedGrid::interpolationRecordT P4RotatedGrid::predictorInterpolation( std::size_t frameIndex,
                                                                           std::size_t searchIndex,
                                                                           std::size_t predictorIndex ) const
{
    if( !configured() )
    {
        throw std::logic_error( "P4RotatedGrid is not configured" );
    }
    return interpolation( frameIndex, skyPredictor( searchIndex, predictorIndex ) );
}

P4RotatedGrid::realT
P4RotatedGrid::sampleTarget( const imageT &rawImage, std::size_t frameIndex, std::size_t searchIndex ) const
{
    requireSampleable( rawImage, frameIndex, searchIndex );
    return sampleRecord( rawImage, targetInterpolation( frameIndex, searchIndex ) );
}

P4RotatedGrid::realT P4RotatedGrid::samplePredictor( const imageT &rawImage,
                                                     std::size_t frameIndex,
                                                     std::size_t searchIndex,
                                                     std::size_t predictorIndex ) const
{
    requireSampleable( rawImage, frameIndex, searchIndex );
    return sampleRecord( rawImage, predictorInterpolation( frameIndex, searchIndex, predictorIndex ) );
}

void P4RotatedGrid::sampleFrame( const imageT &rawImage,
                                 std::size_t frameIndex,
                                 std::size_t searchIndex,
                                 realT &target,
                                 std::vector<realT> &predictors ) const
{
    requireSampleable( rawImage, frameIndex, searchIndex );
    if( predictors.size() != predictorCount() )
    {
        throw std::invalid_argument( "P4RotatedGrid predictor output size must equal the configured K" );
    }

    const realT nextTarget = sampleRecord( rawImage, targetInterpolation( frameIndex, searchIndex ) );
    std::vector<realT> nextPredictors( predictorCount() );
    for( std::size_t predictor = 0; predictor < predictorCount(); ++predictor )
    {
        nextPredictors[predictor] =
            sampleRecord( rawImage, predictorInterpolation( frameIndex, searchIndex, predictor ) );
    }

    target = nextTarget;
    predictors = std::move( nextPredictors );
}

P4RotatedGrid::interpolationRecordT P4RotatedGrid::interpolation( std::size_t frameIndex,
                                                                  const SkyCoordinate &skyCoordinate ) const
{
    if( frameIndex >= m_frameRotations.size() )
    {
        throw std::out_of_range( "P4RotatedGrid frame index is out of range" );
    }
    return makeInterpolation( m_rows, m_columns, m_xCenter, m_yCenter, m_frameRotations[frameIndex], skyCoordinate );
}

P4RotatedGrid::interpolationRecordT P4RotatedGrid::makeInterpolation( int rows,
                                                                      int columns,
                                                                      double xCenter,
                                                                      double yCenter,
                                                                      const FrameRotation &rotation,
                                                                      const SkyCoordinate &skyCoordinate )
{
    const MappedFootprint footprint = mapFootprint( rows, columns, xCenter, yCenter, rotation, skyCoordinate );
    interpolationRecordT::kernelT kernel;
    P4PixelGridf::transformT transform;
    transform( kernel,
               static_cast<realT>( footprint.row - std::floor( footprint.row ) ),
               static_cast<realT>( footprint.column - std::floor( footprint.column ) ) );
    return interpolationRecordT( footprint.row,
                                 footprint.column,
                                 footprint.footprintRow,
                                 footprint.footprintColumn,
                                 kernel,
                                 footprint.inside );
}

P4RotatedGrid::MappedFootprint P4RotatedGrid::mapFootprint( int rows,
                                                            int columns,
                                                            double xCenter,
                                                            double yCenter,
                                                            const FrameRotation &rotation,
                                                            const SkyCoordinate &skyCoordinate )
{
    const double deltaRow = skyCoordinate.row - xCenter;
    const double deltaColumn = skyCoordinate.column - yCenter;
    const double rawRow = deltaRow * rotation.cosine + deltaColumn * rotation.sine + xCenter;
    const double rawColumn = -deltaRow * rotation.sine + deltaColumn * rotation.cosine + yCenter;
    if( !p4RotatedFinite( rawRow ) || !p4RotatedFinite( rawColumn ) )
    {
        throw std::invalid_argument( "P4RotatedGrid mapped detector coordinates must be finite" );
    }

    const double footprintRowValue = std::floor( rawRow ) - static_cast<double>( leftBuffer );
    const double footprintColumnValue = std::floor( rawColumn ) - static_cast<double>( leftBuffer );
    if( footprintRowValue < static_cast<double>( std::numeric_limits<int>::min() ) ||
        footprintRowValue > static_cast<double>( std::numeric_limits<int>::max() ) ||
        footprintColumnValue < static_cast<double>( std::numeric_limits<int>::min() ) ||
        footprintColumnValue > static_cast<double>( std::numeric_limits<int>::max() ) )
    {
        throw std::invalid_argument( "P4RotatedGrid mapped footprint exceeds integer range" );
    }

    const int footprintRow = static_cast<int>( footprintRowValue );
    const int footprintColumn = static_cast<int>( footprintColumnValue );
    const bool inside = footprintRow >= 0 && footprintColumn >= 0 &&
                        static_cast<std::int64_t>( footprintRow ) + width <= rows &&
                        static_cast<std::int64_t>( footprintColumn ) + width <= columns;

    return MappedFootprint{ rawRow, rawColumn, footprintRow, footprintColumn, inside };
}

std::size_t P4RotatedGrid::checkedCoordinateCount( std::size_t searchCount, std::size_t predictorCount )
{
    if( searchCount != 0 && predictorCount > std::numeric_limits<std::size_t>::max() / searchCount )
    {
        throw std::length_error( "P4RotatedGrid fixed sky-coordinate count exceeds size_t range" );
    }
    return searchCount * predictorCount;
}

bool P4RotatedGrid::maskFootprintValid( const imageT *rawCommonMask, const MappedFootprint &footprint )
{
    if( !rawCommonMask )
    {
        return true;
    }
    for( int rowOffset = 0; rowOffset < width; ++rowOffset )
    {
        for( int columnOffset = 0; columnOffset < width; ++columnOffset )
        {
            if( ( *rawCommonMask )( footprint.footprintRow + rowOffset, footprint.footprintColumn + columnOffset ) ==
                0 )
            {
                return false;
            }
        }
    }
    return true;
}

bool P4RotatedGrid::footprintOverlapsExclusion( const MappedFootprint &predictor,
                                                const MappedFootprint &target,
                                                double effectiveRadius )
{
    const double radiusSquared = effectiveRadius * effectiveRadius;
    for( int rowOffset = 0; rowOffset < width; ++rowOffset )
    {
        const double deltaRow = static_cast<double>( predictor.footprintRow + rowOffset ) - target.row;
        for( int columnOffset = 0; columnOffset < width; ++columnOffset )
        {
            const double deltaColumn = static_cast<double>( predictor.footprintColumn + columnOffset ) - target.column;
            if( deltaRow * deltaRow + deltaColumn * deltaColumn <= radiusSquared )
            {
                return true;
            }
        }
    }
    return false;
}

P4RotatedGrid::realT P4RotatedGrid::sampleRecord( const imageT &rawImage, const interpolationRecordT &record )
{
    if( !record.footprintInside() )
    {
        throw std::logic_error( "P4RotatedGrid cannot sample an incomplete detector footprint" );
    }

    realT value{ 0 };
    for( int rowOffset = 0; rowOffset < width; ++rowOffset )
    {
        for( int columnOffset = 0; columnOffset < width; ++columnOffset )
        {
            const realT source = rawImage( record.footprintRow() + rowOffset, record.footprintColumn() + columnOffset );
            if( !mx::math::isFinite( source ) )
            {
                throw std::domain_error( "P4RotatedGrid detector interpolation footprint contains a non-finite value" );
            }
            value += source * record.kernel()( rowOffset, columnOffset );
        }
    }
    if( !mx::math::isFinite( value ) )
    {
        throw std::domain_error( "P4RotatedGrid direct interpolation result is non-finite" );
    }
    return value;
}

void P4RotatedGrid::requireSampleable( const imageT &rawImage, std::size_t frameIndex, std::size_t searchIndex ) const
{
    if( !configured() )
    {
        throw std::logic_error( "P4RotatedGrid is not configured" );
    }
    if( rawImage.rows() != m_rows || rawImage.cols() != m_columns )
    {
        throw std::invalid_argument( "P4RotatedGrid detector image dimensions must match the configured grid" );
    }
    if( frameIndex >= frameCount() )
    {
        throw std::out_of_range( "P4RotatedGrid frame index is out of range" );
    }
    if( !searchPixel( searchIndex ).valid() )
    {
        throw std::logic_error( "P4RotatedGrid cannot sample an all-frame-invalid local fit" );
    }
}

const P4RotatedGrid::SkyCoordinate &P4RotatedGrid::skyPredictor( std::size_t searchIndex,
                                                                 std::size_t predictorIndex ) const
{
    if( searchIndex >= searchPixelCount() || predictorIndex >= predictorCount() )
    {
        throw std::out_of_range( "P4RotatedGrid fixed sky predictor index is out of range" );
    }
    return m_skyPredictors[searchIndex * predictorCount() + predictorIndex];
}

} // namespace improc
} // namespace mx
