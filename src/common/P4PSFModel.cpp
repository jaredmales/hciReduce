/** \file P4PSFModel.cpp
 * \brief Implements compact frozen-model PSF calculations for Pixel Prediction Post-Processing.
 * \author Jared R. Males
 */

#include "P4PSFModel.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include <mx/math/floatUtils.hpp>

namespace mx
{
namespace improc
{

namespace
{

/// Determine whether every value in an Eigen-like array is finite under fast-math.
template <typename arrayT>
bool p4PSFAllFinite( const arrayT &array /**< [in] values to validate */ )
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

/// Sample a finite float PSF template with the production cubic convention and zero exterior padding.
float p4PSFSample( const P4PSFModel::imageT &psfTemplate, /**< [in] finite PSF template */
                   double row,                            /**< [in] template row coordinate */
                   double column /**< [in] template column coordinate */ )
{
    if( !mx::math::isFinite( row ) || !mx::math::isFinite( column ) )
    {
        throw std::invalid_argument( "P4 PSF sample coordinates must be finite" );
    }

    const double floorRow = std::floor( row );
    const double floorColumn = std::floor( column );
    constexpr double integerMinimum = static_cast<double>( std::numeric_limits<int>::min() + 1 );
    constexpr double integerMaximum = static_cast<double>( std::numeric_limits<int>::max() - 2 );
    if( floorRow < integerMinimum || floorRow > integerMaximum || floorColumn < integerMinimum ||
        floorColumn > integerMaximum )
    {
        return 0;
    }

    using transformT = P4PSFModel::gridT::transformT;
    transformT transform;
    P4PSFModel::gridT::kernelT kernel;
    transform( kernel, static_cast<float>( row - floorRow ), static_cast<float>( column - floorColumn ) );

    const int footprintRow = static_cast<int>( floorRow ) - P4PSFModel::gridT::leftBuffer;
    const int footprintColumn = static_cast<int>( floorColumn ) - P4PSFModel::gridT::leftBuffer;
    float value{ 0 };
    for( int rowOffset = 0; rowOffset < P4PSFModel::gridT::width; ++rowOffset )
    {
        const int imageRow = footprintRow + rowOffset;
        if( imageRow < 0 || imageRow >= psfTemplate.rows() )
        {
            continue;
        }
        for( int columnOffset = 0; columnOffset < P4PSFModel::gridT::width; ++columnOffset )
        {
            const int imageColumn = footprintColumn + columnOffset;
            if( imageColumn >= 0 && imageColumn < psfTemplate.cols() )
            {
                value += psfTemplate( imageRow, imageColumn ) * kernel( rowOffset, columnOffset );
            }
        }
    }
    return value;
}

/// Convert a checked signed lookup extent to an Eigen dimension.
Eigen::Index p4PSFLookupExtent( std::int64_t minimum, /**< [in] inclusive minimum lookup coordinate */
                                std::int64_t maximum /**< [in] inclusive maximum lookup coordinate */ )
{
    if( maximum < minimum )
    {
        throw std::logic_error( "P4 PSF lookup bounds are reversed" );
    }
    const std::uint64_t extent = static_cast<std::uint64_t>( maximum - minimum ) + 1;
    if( extent > static_cast<std::uint64_t>( std::numeric_limits<Eigen::Index>::max() ) )
    {
        throw std::length_error( "P4 PSF lookup dimension exceeds Eigen index range" );
    }
    return static_cast<Eigen::Index>( extent );
}

} // namespace

P4PSFModel::P4PSFModel( const imageT &psfTemplate, int stampSize ) : P4PSFModel( psfTemplate, stampSize, stampSize )
{
}

P4PSFModel::P4PSFModel( const imageT &psfTemplate, int stampRows, int stampColumns )
{
    if( psfTemplate.rows() <= 0 || psfTemplate.cols() <= 0 )
    {
        throw std::invalid_argument( "P4 PSF template must be nonempty" );
    }
    if( !p4PSFAllFinite( psfTemplate ) )
    {
        throw std::invalid_argument( "P4 PSF template must contain only finite values" );
    }
    if( stampRows <= 0 || stampColumns <= 0 )
    {
        throw std::invalid_argument( "P4 PSF stamp dimensions must be positive" );
    }

    const double templateCenterRow = 0.5 * static_cast<double>( psfTemplate.rows() - 1 );
    const double templateCenterColumn = 0.5 * static_cast<double>( psfTemplate.cols() - 1 );
    const double stampCenterRow = 0.5 * static_cast<double>( stampRows - 1 );
    const double stampCenterColumn = 0.5 * static_cast<double>( stampColumns - 1 );
    m_minimumRowIndex = static_cast<std::int64_t>( std::floor( stampCenterRow - templateCenterRow ) ) - 3;
    m_minimumColumnIndex = static_cast<std::int64_t>( std::floor( stampCenterColumn - templateCenterColumn ) ) - 3;
    const std::int64_t maximumRowIndex =
        static_cast<std::int64_t>(
            std::ceil( stampCenterRow + static_cast<double>( psfTemplate.rows() - 1 ) - templateCenterRow ) ) +
        3;
    const std::int64_t maximumColumnIndex =
        static_cast<std::int64_t>(
            std::ceil( stampCenterColumn + static_cast<double>( psfTemplate.cols() - 1 ) - templateCenterColumn ) ) +
        3;
    const Eigen::Index shiftedRows = p4PSFLookupExtent( m_minimumRowIndex, maximumRowIndex );
    const Eigen::Index shiftedColumns = p4PSFLookupExtent( m_minimumColumnIndex, maximumColumnIndex );
    if( static_cast<long double>( shiftedRows ) * static_cast<long double>( shiftedColumns ) >
        static_cast<long double>( std::numeric_limits<std::size_t>::max() ) / sizeof( float ) )
    {
        throw std::length_error( "P4 shifted PSF template exceeds size_t range" );
    }

    m_shiftedTemplate.resize( shiftedRows, shiftedColumns );
    for( Eigen::Index column = 0; column < shiftedColumns; ++column )
    {
        const std::int64_t columnIndex = m_minimumColumnIndex + column;
        const double templateColumn = templateCenterColumn - stampCenterColumn + static_cast<double>( columnIndex );
        for( Eigen::Index row = 0; row < shiftedRows; ++row )
        {
            const std::int64_t rowIndex = m_minimumRowIndex + row;
            const double templateRow = templateCenterRow - stampCenterRow + static_cast<double>( rowIndex );
            m_shiftedTemplate( row, column ) = p4PSFSample( psfTemplate, templateRow, templateColumn );
        }
    }
    m_stampRows = stampRows;
    m_stampColumns = stampColumns;
}

int P4PSFModel::stampRows() const noexcept
{
    return m_stampRows;
}

int P4PSFModel::stampColumns() const noexcept
{
    return m_stampColumns;
}

std::size_t P4PSFModel::storageBytes() const noexcept
{
    return static_cast<std::size_t>( m_shiftedTemplate.size() ) * sizeof( float );
}

float P4PSFModel::sampleTemplate( double deltaRow, double deltaColumn ) const
{
    if( !mx::math::isFinite( deltaRow ) || !mx::math::isFinite( deltaColumn ) )
    {
        throw std::invalid_argument( "P4 PSF detector offsets must be finite" );
    }
    const double stampCenterRow = 0.5 * static_cast<double>( m_stampRows - 1 );
    const double stampCenterColumn = 0.5 * static_cast<double>( m_stampColumns - 1 );
    return p4PSFSample( m_shiftedTemplate,
                        stampCenterRow + deltaRow - static_cast<double>( m_minimumRowIndex ),
                        stampCenterColumn + deltaColumn - static_cast<double>( m_minimumColumnIndex ) );
}

float P4PSFModel::shiftedTemplateValue( std::int64_t rowIndex, std::int64_t columnIndex ) const noexcept
{
    const std::int64_t storedRow = rowIndex - m_minimumRowIndex;
    const std::int64_t storedColumn = columnIndex - m_minimumColumnIndex;
    if( storedRow < 0 || storedColumn < 0 || storedRow >= m_shiftedTemplate.rows() ||
        storedColumn >= m_shiftedTemplate.cols() )
    {
        return 0;
    }
    return m_shiftedTemplate( static_cast<Eigen::Index>( storedRow ), static_cast<Eigen::Index>( storedColumn ) );
}

void P4PSFModel::calculateLocalResponse( imageT &output,
                                         const gridT &grid,
                                         std::size_t searchIndex,
                                         const Eigen::Ref<const coefficientT> &coefficients ) const
{
    imageT components;
    calculateLocalResponseComponents( components, grid, searchIndex, {}, 0, coefficients );
    output.resize( m_stampRows, m_stampColumns );
    for( int stampColumn = 0; stampColumn < m_stampColumns; ++stampColumn )
    {
        for( int stampRow = 0; stampRow < m_stampRows; ++stampRow )
        {
            const Eigen::Index stampPixel =
                static_cast<Eigen::Index>( stampRow ) + static_cast<Eigen::Index>( m_stampRows ) * stampColumn;
            output( stampRow, stampColumn ) = components( stampPixel, 0 );
        }
    }
}

void P4PSFModel::calculateLocalResponseComponents( imageT &output,
                                                   const gridT &grid,
                                                   std::size_t searchIndex,
                                                   const std::vector<P4PixelCoordinate> &temporalOffsets,
                                                   std::size_t temporalImageCount,
                                                   const Eigen::Ref<const coefficientT> &coefficients ) const
{
    if( !grid.regionConfigured() )
    {
        throw std::invalid_argument( "P4 PSF calculation requires a complete detector-frame pixel grid" );
    }
    const P4SearchPixelRecord &search = grid.searchPixel( searchIndex );
    if( !search.valid() )
    {
        throw std::invalid_argument( "P4 PSF calculation requires a valid local fit" );
    }
    if( temporalImageCount != 0 && temporalOffsets.empty() )
    {
        throw std::invalid_argument( "P4 temporal PSF calculation requires predictor offsets" );
    }
    if( temporalImageCount != 0 &&
        temporalImageCount > std::numeric_limits<std::size_t>::max() / temporalOffsets.size() )
    {
        throw std::length_error( "P4 temporal PSF predictor count overflows size_t" );
    }
    const std::size_t temporalPredictorCount = temporalImageCount * temporalOffsets.size();
    if( grid.predictorCount() > std::numeric_limits<std::size_t>::max() - temporalPredictorCount )
    {
        throw std::length_error( "P4 PSF predictor count overflows size_t" );
    }
    const std::size_t predictorCount = grid.predictorCount() + temporalPredictorCount;
    if( predictorCount > static_cast<std::size_t>( std::numeric_limits<Eigen::Index>::max() ) ||
        coefficients.rows() != static_cast<Eigen::Index>( predictorCount ) )
    {
        throw std::invalid_argument( "P4 PSF coefficient count must match same-image and temporal predictors" );
    }
    if( !p4PSFAllFinite( coefficients ) )
    {
        throw std::invalid_argument( "P4 PSF coefficients must contain only finite values" );
    }

    const P4PixelCoordinate &target = search.coordinate();
    const Eigen::Index stampPixels =
        static_cast<Eigen::Index>( m_stampRows ) * static_cast<Eigen::Index>( m_stampColumns );
    if( temporalImageCount == std::numeric_limits<std::size_t>::max() ||
        temporalImageCount + 1 > static_cast<std::size_t>( std::numeric_limits<Eigen::Index>::max() ) )
    {
        throw std::length_error( "P4 PSF temporal component count exceeds Eigen range" );
    }
    output.resize( stampPixels, static_cast<Eigen::Index>( temporalImageCount + 1 ) );
    for( int stampColumn = 0; stampColumn < m_stampColumns; ++stampColumn )
    {
        for( int stampRow = 0; stampRow < m_stampRows; ++stampRow )
        {
            double residual = shiftedTemplateValue( stampRow, stampColumn );

            for( std::size_t predictor = 0; predictor < grid.predictorCount(); ++predictor )
            {
                const gridT::interpolationRecordT &record = grid.interpolation( searchIndex, predictor );
                float predictorSample{ 0 };
                for( int rowOffset = 0; rowOffset < gridT::width; ++rowOffset )
                {
                    const std::int64_t rowIndex =
                        static_cast<std::int64_t>( stampRow ) + record.footprintRow() + rowOffset - target.row();
                    for( int columnOffset = 0; columnOffset < gridT::width; ++columnOffset )
                    {
                        const std::int64_t columnIndex = static_cast<std::int64_t>( stampColumn ) +
                                                         record.footprintColumn() + columnOffset - target.column();
                        predictorSample +=
                            shiftedTemplateValue( rowIndex, columnIndex ) * record.kernel()( rowOffset, columnOffset );
                    }
                }
                residual -=
                    coefficients( static_cast<Eigen::Index>( predictor ) ) * static_cast<double>( predictorSample );
            }

            if( !mx::math::isFinite( residual ) )
            {
                throw std::runtime_error( "P4 local PSF calculation produced a nonfinite value" );
            }
            const float stored = static_cast<float>( residual );
            if( !mx::math::isFinite( stored ) )
            {
                throw std::overflow_error( "P4 local PSF value exceeds float storage range" );
            }
            const Eigen::Index stampPixel =
                static_cast<Eigen::Index>( stampRow ) + static_cast<Eigen::Index>( m_stampRows ) * stampColumn;
            output( stampPixel, 0 ) = stored;

            for( std::size_t temporalImage = 0; temporalImage < temporalImageCount; ++temporalImage )
            {
                double temporalResponse{ 0 };
                for( std::size_t predictor = 0; predictor < temporalOffsets.size(); ++predictor )
                {
                    const P4PixelCoordinate &offset = temporalOffsets[predictor];
                    const std::size_t coefficientIndex =
                        grid.predictorCount() + temporalImage * temporalOffsets.size() + predictor;
                    temporalResponse -= coefficients( static_cast<Eigen::Index>( coefficientIndex ) ) *
                                        static_cast<double>( shiftedTemplateValue( stampRow + offset.row(),
                                                                                   stampColumn + offset.column() ) );
                }
                if( !mx::math::isFinite( temporalResponse ) )
                {
                    throw std::runtime_error( "P4 temporal PSF component calculation produced a nonfinite value" );
                }
                const float temporalStored = static_cast<float>( temporalResponse );
                if( !mx::math::isFinite( temporalStored ) )
                {
                    throw std::overflow_error( "P4 temporal PSF component exceeds float storage range" );
                }
                output( stampPixel, static_cast<Eigen::Index>( temporalImage + 1 ) ) = temporalStored;
            }
        }
    }
}

} // namespace improc
} // namespace mx
