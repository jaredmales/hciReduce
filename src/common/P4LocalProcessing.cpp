/** \file P4LocalProcessing.cpp
 * \brief Implements sparse geometry and trial-source sampling for pixel-local P4 processing.
 * \author Jared R. Males
 */

#include "P4LocalProcessing.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

#include <mx/math/floatUtils.hpp>
#include <mx/math/geo.hpp>

namespace mx
{
namespace improc
{

namespace
{

/// Report whether a double remains finite under fast-math.
bool p4LocalFinite( double value /**< [in] scalar to validate */ ) noexcept
{
    return mx::math::isFinite( value );
}

/// Return a phase-preserving crop dimension at least as large as the requested width.
int p4LocalCropDimension( int requested, /**< [in] positive nominal crop dimension */
                          int source,    /**< [in] positive source-template dimension */
                          int detector /**< [in] positive detector-image dimension */ )
{
    if( requested <= 0 || source <= 0 || detector <= 0 )
    {
        throw std::invalid_argument( "P4 local source dimensions must be positive" );
    }
    if( source % 2 != detector % 2 )
    {
        throw std::invalid_argument( "P4 local source template and detector dimensions must have matching parity" );
    }

    int actual = requested;
    if( actual % 2 != detector % 2 )
    {
        if( actual == std::numeric_limits<int>::max() )
        {
            throw std::invalid_argument( "P4 local source crop dimension exceeds integer range" );
        }
        ++actual;
    }
    if( actual > source || actual > detector )
    {
        throw std::invalid_argument( "P4 local source crop exceeds the template or detector dimensions" );
    }
    return actual;
}

/// Temporary detector coordinate used while constructing a sparse dependency map.
using p4LocalCoordinateT = std::pair<int, int>;

/// Temporary ordered detector footprint for one local output sample.
struct P4LocalPendingOutput
{
    bool valid{ false }; ///< Whether every required detector pixel has P4 ownership.

    int frame{ -1 };     ///< Physical target frame supplying every dependency in this output record.

    std::vector<std::pair<p4LocalCoordinateT, float>> samples; ///< Detector coordinates and cubic weights.
};

} // namespace

P4LocalSearchRequest::P4LocalSearchRequest( const P4PixelCoordinate &coordinate,
                                            int region,
                                            std::size_t searchIndex,
                                            std::vector<int> frames )
    : m_coordinate( coordinate ), m_region( region ), m_searchIndex( searchIndex ), m_frames( std::move( frames ) )
{
    if( region < 0 || m_frames.empty() || !std::is_sorted( m_frames.begin(), m_frames.end() ) ||
        std::adjacent_find( m_frames.begin(), m_frames.end() ) != m_frames.end() || m_frames.front() < 0 )
    {
        throw std::invalid_argument( "P4 local search request is invalid" );
    }
}

const P4PixelCoordinate &P4LocalSearchRequest::coordinate() const noexcept
{
    return m_coordinate;
}

int P4LocalSearchRequest::region() const noexcept
{
    return m_region;
}

std::size_t P4LocalSearchRequest::searchIndex() const noexcept
{
    return m_searchIndex;
}

const std::vector<int> &P4LocalSearchRequest::frames() const noexcept
{
    return m_frames;
}

P4LocalResidualSample::P4LocalResidualSample( std::size_t requestIndex, std::size_t frameOffset, float weight )
    : m_requestIndex( requestIndex ), m_frameOffset( frameOffset ), m_weight( weight )
{
    if( !mx::math::isFinite( weight ) )
    {
        throw std::invalid_argument( "P4 local residual weight must be finite" );
    }
}

std::size_t P4LocalResidualSample::requestIndex() const noexcept
{
    return m_requestIndex;
}

std::size_t P4LocalResidualSample::frameOffset() const noexcept
{
    return m_frameOffset;
}

float P4LocalResidualSample::weight() const noexcept
{
    return m_weight;
}

P4LocalOutputSample::P4LocalOutputSample() = default;

P4LocalOutputSample::P4LocalOutputSample( std::vector<P4LocalResidualSample> samples )
    : m_valid( true ), m_samples( std::move( samples ) )
{
    if( m_samples.empty() )
    {
        throw std::invalid_argument( "P4 local valid output sample requires residual dependencies" );
    }
}

bool P4LocalOutputSample::valid() const noexcept
{
    return m_valid;
}

const std::vector<P4LocalResidualSample> &P4LocalOutputSample::samples() const noexcept
{
    return m_samples;
}

P4LocalGeometry::P4LocalGeometry() = default;

void P4LocalGeometry::configure( int rows,
                                 int columns,
                                 int stampSize,
                                 double sourceRow,
                                 double sourceColumn,
                                 const std::vector<double> &angles,
                                 bool derotate,
                                 const lookupImageT &ownership,
                                 const lookupImageT &searchIndexLookup )
{
    if( rows <= 0 || columns <= 0 || stampSize <= 0 || stampSize % 2 == 0 || angles.empty() )
    {
        throw std::invalid_argument( "P4 local geometry requires positive dimensions, an odd stamp, and frames" );
    }
    if( ownership.rows() != rows || ownership.cols() != columns || searchIndexLookup.rows() != rows ||
        searchIndexLookup.cols() != columns )
    {
        throw std::invalid_argument( "P4 local ownership lookups must match the detector dimensions" );
    }
    if( !p4LocalFinite( sourceRow ) || !p4LocalFinite( sourceColumn ) )
    {
        throw std::invalid_argument( "P4 local source coordinates must be finite" );
    }
    if( !std::all_of( angles.begin(), angles.end(), p4LocalFinite ) )
    {
        throw std::invalid_argument( "P4 local derotation angles must be finite" );
    }

    const int centerRow = static_cast<int>( std::floor( sourceRow + 0.5 ) );
    const int centerColumn = static_cast<int>( std::floor( sourceColumn + 0.5 ) );
    const int originRow = centerRow - stampSize / 2;
    const int originColumn = centerColumn - stampSize / 2;
    if( originRow < 0 || originColumn < 0 || originRow > rows - stampSize || originColumn > columns - stampSize )
    {
        throw std::invalid_argument( "P4 local output stamp crosses the full-image edge" );
    }

    const long double recordCountValue = static_cast<long double>( stampSize ) * static_cast<long double>( stampSize ) *
                                         static_cast<long double>( angles.size() );
    if( recordCountValue > static_cast<long double>( std::numeric_limits<std::size_t>::max() ) )
    {
        throw std::length_error( "P4 local output dependency count exceeds size_t range" );
    }
    const std::size_t recordCount = static_cast<std::size_t>( recordCountValue );
    std::vector<P4LocalPendingOutput> pending( recordCount );
    std::map<p4LocalCoordinateT, std::set<int>> requestedFrames;
    const float xCenter = 0.5F * static_cast<float>( rows - 1 );
    const float yCenter = 0.5F * static_cast<float>( columns - 1 );
    const int width = P4PixelGridf::width;
    const int leftBuffer = P4PixelGridf::leftBuffer;
    const int upperRow = rows - width + leftBuffer;
    const int upperColumn = columns - width + leftBuffer;
    cubicConvolTransform<float> transform;

    const auto pendingOffset = [stampSize, &angles]( int stampRow, int stampColumn, std::size_t frame )
    {
        return frame +
               angles.size() * ( static_cast<std::size_t>( stampRow ) +
                                 static_cast<std::size_t>( stampSize ) * static_cast<std::size_t>( stampColumn ) );
    };

    for( int stampColumn = 0; stampColumn < stampSize; ++stampColumn )
    {
        for( int stampRow = 0; stampRow < stampSize; ++stampRow )
        {
            const int outputRow = originRow + stampRow;
            const int outputColumn = originColumn + stampColumn;
            for( std::size_t frame = 0; frame < angles.size(); ++frame )
            {
                P4LocalPendingOutput &output = pending[pendingOffset( stampRow, stampColumn, frame )];
                if( !derotate || angles[frame] == 0 )
                {
                    if( ownership( outputRow, outputColumn ) < 0 || searchIndexLookup( outputRow, outputColumn ) < 0 )
                    {
                        continue;
                    }
                    output.valid = true;
                    output.frame = static_cast<int>( frame );
                    output.samples.push_back( { { outputRow, outputColumn }, 1 } );
                    requestedFrames[{ outputRow, outputColumn }].insert( static_cast<int>( frame ) );
                    continue;
                }

                const float angle = static_cast<float>( angles[frame] );
                const float cosine = std::cos( angle );
                const float sine = std::sin( angle );
                float centerRowCosine = xCenter * cosine;
                float centerRowSine = xCenter * sine;
                const float centerColumnCosine = yCenter * cosine;
                const float centerColumnSine = yCenter * sine;
                centerRowCosine += centerColumnSine;
                centerRowSine -= centerColumnCosine;
                const float rowCosine = static_cast<float>( outputRow ) * cosine - centerRowCosine;
                const float rowSine = -( static_cast<float>( outputRow ) * sine - centerRowSine );
                const float detectorRow = rowCosine + static_cast<float>( outputColumn ) * sine + xCenter;
                const float detectorColumn = rowSine + static_cast<float>( outputColumn ) * cosine + yCenter;
                const int anchorRow = static_cast<int>( detectorRow );
                const int anchorColumn = static_cast<int>( detectorColumn );
                if( anchorRow <= leftBuffer || anchorRow >= upperRow || anchorColumn <= leftBuffer ||
                    anchorColumn >= upperColumn )
                {
                    continue;
                }

                P4PixelGridf::kernelT kernel;
                transform( kernel,
                           detectorRow - static_cast<float>( anchorRow ),
                           detectorColumn - static_cast<float>( anchorColumn ) );
                const int footprintRow = anchorRow - leftBuffer;
                const int footprintColumn = anchorColumn - leftBuffer;
                bool complete{ true };
                for( int columnOffset = 0; columnOffset < width && complete; ++columnOffset )
                {
                    for( int rowOffset = 0; rowOffset < width; ++rowOffset )
                    {
                        const int detectorPixelRow = footprintRow + rowOffset;
                        const int detectorPixelColumn = footprintColumn + columnOffset;
                        if( ownership( detectorPixelRow, detectorPixelColumn ) < 0 ||
                            searchIndexLookup( detectorPixelRow, detectorPixelColumn ) < 0 )
                        {
                            complete = false;
                            break;
                        }
                    }
                }
                if( !complete )
                {
                    continue;
                }

                output.valid = true;
                output.frame = static_cast<int>( frame );
                output.samples.reserve( static_cast<std::size_t>( width * width ) );
                for( int columnOffset = 0; columnOffset < width; ++columnOffset )
                {
                    for( int rowOffset = 0; rowOffset < width; ++rowOffset )
                    {
                        const int detectorPixelRow = footprintRow + rowOffset;
                        const int detectorPixelColumn = footprintColumn + columnOffset;
                        output.samples.push_back(
                            { { detectorPixelRow, detectorPixelColumn }, kernel( rowOffset, columnOffset ) } );
                        requestedFrames[{ detectorPixelRow, detectorPixelColumn }].insert( static_cast<int>( frame ) );
                    }
                }
            }
        }
    }

    std::vector<P4LocalSearchRequest> searchRequests;
    searchRequests.reserve( requestedFrames.size() );
    std::map<p4LocalCoordinateT, std::size_t> requestIndices;
    std::size_t sparseSampleCount{ 0 };
    for( const auto &[coordinate, frames] : requestedFrames )
    {
        const std::int64_t region = ownership( coordinate.first, coordinate.second );
        const std::int64_t searchIndex = searchIndexLookup( coordinate.first, coordinate.second );
        if( region < 0 || region > std::numeric_limits<int>::max() || searchIndex < 0 )
        {
            throw std::logic_error( "P4 local requested coordinate lost its ownership lookup" );
        }
        std::vector<int> orderedFrames( frames.begin(), frames.end() );
        if( sparseSampleCount > std::numeric_limits<std::size_t>::max() - orderedFrames.size() )
        {
            throw std::length_error( "P4 local sparse sample count exceeds size_t range" );
        }
        sparseSampleCount += orderedFrames.size();
        const std::size_t requestIndex = searchRequests.size();
        requestIndices.emplace( coordinate, requestIndex );
        searchRequests.emplace_back( P4PixelCoordinate( coordinate.first, coordinate.second ),
                                     static_cast<int>( region ),
                                     static_cast<std::size_t>( searchIndex ),
                                     std::move( orderedFrames ) );
    }

    std::vector<P4LocalOutputSample> outputSamples;
    outputSamples.reserve( pending.size() );
    for( P4LocalPendingOutput &output : pending )
    {
        if( !output.valid )
        {
            outputSamples.emplace_back();
            continue;
        }
        std::vector<P4LocalResidualSample> samples;
        samples.reserve( output.samples.size() );
        for( const auto &[coordinate, weight] : output.samples )
        {
            const std::size_t requestIndex = requestIndices.at( coordinate );
            const std::vector<int> &frames = searchRequests[requestIndex].frames();
            const auto found = std::lower_bound( frames.begin(), frames.end(), output.frame );
            if( found == frames.end() || *found != output.frame )
            {
                throw std::logic_error( "P4 local output frame lost its sparse request" );
            }
            samples.emplace_back( requestIndex,
                                  static_cast<std::size_t>( std::distance( frames.begin(), found ) ),
                                  weight );
        }
        outputSamples.emplace_back( std::move( samples ) );
    }

    m_stampSize = stampSize;
    m_originRow = originRow;
    m_originColumn = originColumn;
    m_sourceRow = sourceRow;
    m_sourceColumn = sourceColumn;
    m_frameCount = angles.size();
    m_sparseSampleCount = sparseSampleCount;
    m_searchRequests = std::move( searchRequests );
    m_outputSamples = std::move( outputSamples );
}

bool P4LocalGeometry::configured() const noexcept
{
    return m_stampSize > 0 && m_frameCount > 0 && !m_outputSamples.empty();
}

int P4LocalGeometry::stampSize() const noexcept
{
    return m_stampSize;
}

int P4LocalGeometry::originRow() const noexcept
{
    return m_originRow;
}

int P4LocalGeometry::originColumn() const noexcept
{
    return m_originColumn;
}

double P4LocalGeometry::sourceRow() const noexcept
{
    return m_sourceRow;
}

double P4LocalGeometry::sourceColumn() const noexcept
{
    return m_sourceColumn;
}

std::size_t P4LocalGeometry::frameCount() const noexcept
{
    return m_frameCount;
}

const std::vector<P4LocalSearchRequest> &P4LocalGeometry::searchRequests() const noexcept
{
    return m_searchRequests;
}

std::size_t P4LocalGeometry::outputOffset( int stampRow, int stampColumn, std::size_t frame ) const
{
    if( stampRow < 0 || stampRow >= m_stampSize || stampColumn < 0 || stampColumn >= m_stampSize ||
        frame >= m_frameCount )
    {
        throw std::out_of_range( "P4 local output sample index is out of range" );
    }
    return frame + m_frameCount * ( static_cast<std::size_t>( stampRow ) +
                                    static_cast<std::size_t>( m_stampSize ) * static_cast<std::size_t>( stampColumn ) );
}

const P4LocalOutputSample &P4LocalGeometry::outputSample( int stampRow, int stampColumn, std::size_t frame ) const
{
    if( !configured() )
    {
        throw std::logic_error( "P4 local geometry is not configured" );
    }
    return m_outputSamples[outputOffset( stampRow, stampColumn, frame )];
}

std::size_t P4LocalGeometry::sparseSampleCount() const noexcept
{
    return m_sparseSampleCount;
}

std::size_t P4LocalGeometry::storageBytes() const
{
    const auto checkedAdd = []( std::size_t left, std::size_t right )
    {
        if( left > std::numeric_limits<std::size_t>::max() - right )
        {
            throw std::overflow_error( "P4 local geometry storage exceeds size_t range" );
        }
        return left + right;
    };
    const auto checkedMultiply = []( std::size_t count, std::size_t elementSize )
    {
        if( elementSize != 0 && count > std::numeric_limits<std::size_t>::max() / elementSize )
        {
            throw std::overflow_error( "P4 local geometry storage exceeds size_t range" );
        }
        return count * elementSize;
    };

    std::size_t bytes = checkedMultiply( m_searchRequests.capacity(), sizeof( P4LocalSearchRequest ) );
    for( const P4LocalSearchRequest &request : m_searchRequests )
    {
        bytes = checkedAdd( bytes, checkedMultiply( request.frames().capacity(), sizeof( int ) ) );
    }
    bytes = checkedAdd( bytes, checkedMultiply( m_outputSamples.capacity(), sizeof( P4LocalOutputSample ) ) );
    for( const P4LocalOutputSample &sample : m_outputSamples )
    {
        bytes = checkedAdd( bytes, checkedMultiply( sample.samples().capacity(), sizeof( P4LocalResidualSample ) ) );
    }
    return bytes;
}

P4TrialSource::P4TrialSource() = default;

void P4TrialSource::configure( const imageT &sourceTemplate,
                               int detectorRows,
                               int detectorColumns,
                               int requestedStampSize,
                               const std::vector<double> &angles,
                               double separation,
                               double positionAngle,
                               double contrast,
                               const std::vector<float> &scales )
{
    if( sourceTemplate.rows() <= 0 || sourceTemplate.cols() <= 0 || detectorRows <= 0 || detectorColumns <= 0 ||
        requestedStampSize <= 0 || angles.empty() || angles.size() != scales.size() )
    {
        throw std::invalid_argument( "P4 local trial source dimensions and frame vectors must be nonempty" );
    }
    if( !p4LocalFinite( separation ) || separation < 0 || !p4LocalFinite( positionAngle ) ||
        !p4LocalFinite( contrast ) || contrast == 0 || !std::all_of( angles.begin(), angles.end(), p4LocalFinite ) )
    {
        throw std::invalid_argument( "P4 local trial position and contrast must be finite and valid" );
    }
    for( int column = 0; column < sourceTemplate.cols(); ++column )
    {
        for( int row = 0; row < sourceTemplate.rows(); ++row )
        {
            if( !mx::math::isFinite( sourceTemplate( row, column ) ) )
            {
                throw std::invalid_argument( "P4 local source template must contain only finite values" );
            }
        }
    }
    if( !std::all_of( scales.begin(), scales.end(), []( float scale ) { return mx::math::isFinite( scale ); } ) )
    {
        throw std::invalid_argument( "P4 local source scales must be finite" );
    }

    const int cropRows = p4LocalCropDimension( requestedStampSize, sourceTemplate.rows(), detectorRows );
    const int cropColumns = p4LocalCropDimension( requestedStampSize, sourceTemplate.cols(), detectorColumns );
    const int sourceRow = ( sourceTemplate.rows() - cropRows ) / 2;
    const int sourceColumn = ( sourceTemplate.cols() - cropColumns ) / 2;
    const int detectorRow = ( detectorRows - cropRows ) / 2;
    const int detectorColumn = ( detectorColumns - cropColumns ) / 2;
    imageT prepared = imageT::Zero( detectorRows, detectorColumns );
    prepared.block( detectorRow, detectorColumn, cropRows, cropColumns ) =
        sourceTemplate.block( sourceRow, sourceColumn, cropRows, cropColumns );

    std::vector<FrameShift> shifts;
    shifts.reserve( angles.size() );
    cubicConvolTransform<float> transform;
    const float separationValue = static_cast<float>( separation );
    const float positionAngleRadians = mx::math::dtor( -static_cast<float>( positionAngle ) );
    for( std::size_t frame = 0; frame < angles.size(); ++frame )
    {
        const float angle = positionAngleRadians + static_cast<float>( angles[frame] );
        FrameShift shift;
        shift.rowShift = separationValue * std::sin( angle );
        shift.columnShift = separationValue * std::cos( angle );
        shift.scale = static_cast<float>( contrast * static_cast<double>( scales[frame] ) );
        if( !p4LocalFinite( shift.rowShift ) || !p4LocalFinite( shift.columnShift ) ||
            !mx::math::isFinite( shift.scale ) )
        {
            throw std::invalid_argument( "P4 local source shift or scaled contrast is non-finite" );
        }
        shift.integral =
            shift.rowShift == std::floor( shift.rowShift ) && shift.columnShift == std::floor( shift.columnShift );
        if( !shift.integral )
        {
            const float rowPhase = 1.0F - ( shift.rowShift - std::floor( shift.rowShift ) );
            const float columnPhase = 1.0F - ( shift.columnShift - std::floor( shift.columnShift ) );
            transform( shift.kernel, rowPhase, columnPhase );
        }
        shifts.push_back( std::move( shift ) );
    }

    m_template = std::move( prepared );
    m_cropRows = cropRows;
    m_cropColumns = cropColumns;
    m_shifts = std::move( shifts );
}

bool P4TrialSource::configured() const noexcept
{
    return m_template.rows() > 0 && m_template.cols() > 0 && !m_shifts.empty();
}

int P4TrialSource::cropRows() const noexcept
{
    return m_cropRows;
}

int P4TrialSource::cropColumns() const noexcept
{
    return m_cropColumns;
}

float P4TrialSource::value( std::size_t frame, int row, int column ) const
{
    if( !configured() )
    {
        throw std::logic_error( "P4 local trial source is not configured" );
    }
    if( frame >= m_shifts.size() || row < 0 || row >= m_template.rows() || column < 0 || column >= m_template.cols() )
    {
        throw std::out_of_range( "P4 local trial-source sample index is out of range" );
    }

    const FrameShift &shift = m_shifts[frame];
    if( shift.integral )
    {
        const int sourceRow = row - static_cast<int>( shift.rowShift );
        const int sourceColumn = column - static_cast<int>( shift.columnShift );
        if( sourceRow < 0 || sourceRow >= m_template.rows() || sourceColumn < 0 || sourceColumn >= m_template.cols() )
        {
            return 0;
        }
        return shift.scale * m_template( sourceRow, sourceColumn );
    }

    const int anchorRow = static_cast<int>( static_cast<float>( row ) - shift.rowShift );
    const int anchorColumn = static_cast<int>( static_cast<float>( column ) - shift.columnShift );
    const int leftBuffer = P4PixelGridf::leftBuffer;
    const int upperRow = m_template.rows() - P4PixelGridf::width + leftBuffer;
    const int upperColumn = m_template.cols() - P4PixelGridf::width + leftBuffer;
    if( anchorRow <= leftBuffer || anchorRow >= upperRow || anchorColumn <= leftBuffer || anchorColumn >= upperColumn )
    {
        return 0;
    }

    float value{ 0 };
    const int footprintRow = anchorRow - leftBuffer;
    const int footprintColumn = anchorColumn - leftBuffer;
    for( int columnOffset = 0; columnOffset < P4PixelGridf::width; ++columnOffset )
    {
        for( int rowOffset = 0; rowOffset < P4PixelGridf::width; ++rowOffset )
        {
            value += m_template( footprintRow + rowOffset, footprintColumn + columnOffset ) *
                     shift.kernel( rowOffset, columnOffset );
        }
    }
    return shift.scale * value;
}

float P4TrialSource::sample( std::size_t frame, const P4PixelGridf::interpolationRecordT &record ) const
{
    float result{ 0 };
    for( int columnOffset = 0; columnOffset < P4PixelGridf::width; ++columnOffset )
    {
        for( int rowOffset = 0; rowOffset < P4PixelGridf::width; ++rowOffset )
        {
            result += value( frame, record.footprintRow() + rowOffset, record.footprintColumn() + columnOffset ) *
                      record.kernel()( rowOffset, columnOffset );
        }
    }
    return result;
}

} // namespace improc
} // namespace mx
