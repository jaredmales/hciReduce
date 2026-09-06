/** \file KLIPPSFModel.cpp
 * \brief Implements frozen-basis KLIP response calculations for sparse PSF probes.
 * \author Jared R. Males
 */

#include "KLIPPSFModel.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace mx
{
namespace improc
{

namespace
{

/// Return whether every value in an Eigen-like array is finite.
template <typename arrayT>
bool klipPSFAllFinite( const arrayT &array /**< [in] values to validate */ )
{
    for( Eigen::Index column = 0; column < array.cols(); ++column )
    {
        for( Eigen::Index row = 0; row < array.rows(); ++row )
        {
            if( !std::isfinite( array( row, column ) ) )
            {
                return false;
            }
        }
    }
    return true;
}

} // namespace

KLIPPSFLinearAccumulator::KLIPPSFLinearAccumulator( std::size_t sampleCount,
                                                    std::size_t modeCount,
                                                    int stampSize,
                                                    std::size_t workerCount,
                                                    std::size_t frameCount,
                                                    double minimumGoodFraction )
    : m_sampleCount( sampleCount ), m_modeCount( modeCount ), m_stampSize( stampSize ), m_workerCount( workerCount ),
      m_frameCount( frameCount ), m_minimumGoodFraction( minimumGoodFraction )
{
    if( sampleCount == 0 || modeCount == 0 || stampSize <= 0 || stampSize % 2 == 0 || workerCount == 0 ||
        frameCount == 0 || frameCount > std::numeric_limits<std::uint32_t>::max() ||
        !std::isfinite( minimumGoodFraction ) || minimumGoodFraction < 0 || minimumGoodFraction > 1 )
    {
        throw std::invalid_argument( "KLIP linear PSF accumulator dimensions and support fraction must be valid" );
    }
    if( sampleCount > std::numeric_limits<std::size_t>::max() / modeCount ||
        sampleCount * modeCount > std::numeric_limits<std::size_t>::max() / workerCount ||
        sampleCount > std::numeric_limits<std::size_t>::max() / workerCount )
    {
        throw std::length_error( "KLIP linear PSF accumulator column count overflows size_t" );
    }
    const std::size_t responseColumns = sampleCount * modeCount * workerCount;
    const std::size_t supportColumns = sampleCount * workerCount;
    if( responseColumns > static_cast<std::size_t>( std::numeric_limits<Eigen::Index>::max() ) ||
        supportColumns > static_cast<std::size_t>( std::numeric_limits<Eigen::Index>::max() ) )
    {
        throw std::length_error( "KLIP linear PSF accumulator column count exceeds Eigen range" );
    }
    const Eigen::Index stampPixels = static_cast<Eigen::Index>( stampSize ) * stampSize;
    m_responseSums = imageT::Zero( stampPixels, static_cast<Eigen::Index>( responseColumns ) );
    m_weightSums = imageT::Zero( stampPixels, static_cast<Eigen::Index>( supportColumns ) );
    m_validCounts = Eigen::Array<std::uint32_t, Eigen::Dynamic, Eigen::Dynamic>::Zero(
        stampPixels,
        static_cast<Eigen::Index>( supportColumns ) );
}

void KLIPPSFLinearAccumulator::addResponse( std::size_t workerIndex,
                                            std::size_t sampleIndex,
                                            std::size_t modeIndex,
                                            float frameWeight,
                                            const imageT &contribution )
{
    if( workerIndex >= m_workerCount || sampleIndex >= m_sampleCount || modeIndex >= m_modeCount ||
        contribution.rows() != m_stampSize || contribution.cols() != m_stampSize || !std::isfinite( frameWeight ) ||
        !klipPSFAllFinite( contribution ) )
    {
        throw std::invalid_argument( "KLIP linear PSF response contribution is inconsistent with accumulator" );
    }
    const Eigen::Index column =
        static_cast<Eigen::Index>( workerIndex * m_sampleCount * m_modeCount + sampleIndex * m_modeCount + modeIndex );
    for( int stampColumn = 0; stampColumn < m_stampSize; ++stampColumn )
    {
        for( int stampRow = 0; stampRow < m_stampSize; ++stampRow )
        {
            const Eigen::Index pixel = stampRow + static_cast<Eigen::Index>( m_stampSize ) * stampColumn;
            m_responseSums( pixel, column ) += frameWeight * contribution( stampRow, stampColumn );
        }
    }
}

void KLIPPSFLinearAccumulator::addFrameSupport( std::size_t workerIndex,
                                                std::size_t sampleIndex,
                                                float frameWeight,
                                                const validityT &validity )
{
    if( workerIndex >= m_workerCount || sampleIndex >= m_sampleCount || validity.rows() != m_stampSize ||
        validity.cols() != m_stampSize || !std::isfinite( frameWeight ) )
    {
        throw std::invalid_argument( "KLIP linear PSF frame support is inconsistent with accumulator" );
    }
    const Eigen::Index column = static_cast<Eigen::Index>( workerIndex * m_sampleCount + sampleIndex );
    for( int stampColumn = 0; stampColumn < m_stampSize; ++stampColumn )
    {
        for( int stampRow = 0; stampRow < m_stampSize; ++stampRow )
        {
            if( validity( stampRow, stampColumn ) == 0 )
            {
                continue;
            }
            const Eigen::Index pixel = stampRow + static_cast<Eigen::Index>( m_stampSize ) * stampColumn;
            if( m_validCounts( pixel, column ) == std::numeric_limits<std::uint32_t>::max() )
            {
                throw std::overflow_error( "KLIP linear PSF valid-frame count overflow" );
            }
            ++m_validCounts( pixel, column );
            m_weightSums( pixel, column ) += frameWeight;
        }
    }
}

void KLIPPSFLinearAccumulator::finalize( std::vector<imageT> &responses, std::vector<validityT> &validities ) const
{
    const std::size_t responseCount = m_sampleCount * m_modeCount;
    responses.assign( responseCount, imageT::Zero( m_stampSize, m_stampSize ) );
    validities.assign( responseCount, validityT::Zero( m_stampSize, m_stampSize ) );
    for( std::size_t sample = 0; sample < m_sampleCount; ++sample )
    {
        for( int stampColumn = 0; stampColumn < m_stampSize; ++stampColumn )
        {
            for( int stampRow = 0; stampRow < m_stampSize; ++stampRow )
            {
                const Eigen::Index pixel = stampRow + static_cast<Eigen::Index>( m_stampSize ) * stampColumn;
                std::uint64_t validCount{ 0 };
                float weightSum{ 0 };
                for( std::size_t worker = 0; worker < m_workerCount; ++worker )
                {
                    const Eigen::Index supportColumn = static_cast<Eigen::Index>( worker * m_sampleCount + sample );
                    validCount += m_validCounts( pixel, supportColumn );
                    weightSum += m_weightSums( pixel, supportColumn );
                }
                const bool enoughSupport =
                    validCount > 0 &&
                    static_cast<double>( validCount ) >= m_minimumGoodFraction * static_cast<double>( m_frameCount );
                if( !enoughSupport || !std::isfinite( weightSum ) || weightSum == 0 )
                {
                    continue;
                }
                for( std::size_t mode = 0; mode < m_modeCount; ++mode )
                {
                    const std::size_t responseIndex = sample * m_modeCount + mode;
                    float responseSum{ 0 };
                    for( std::size_t worker = 0; worker < m_workerCount; ++worker )
                    {
                        const Eigen::Index responseColumn =
                            static_cast<Eigen::Index>( worker * responseCount + sample * m_modeCount + mode );
                        responseSum += m_responseSums( pixel, responseColumn );
                    }
                    const float value = responseSum / weightSum;
                    if( std::isfinite( value ) )
                    {
                        responses[responseIndex]( stampRow, stampColumn ) = value;
                        validities[responseIndex]( stampRow, stampColumn ) = 1;
                    }
                }
            }
        }
    }
}

std::size_t KLIPPSFLinearAccumulator::storageBytes() const noexcept
{
    return static_cast<std::size_t>( m_responseSums.size() ) * sizeof( imageT::Scalar ) +
           static_cast<std::size_t>( m_weightSums.size() ) * sizeof( imageT::Scalar ) +
           static_cast<std::size_t>( m_validCounts.size() ) * sizeof( std::uint32_t );
}

KLIPPSFModel::KLIPPSFModel( const imageT &psfTemplate,
                            int detectorRows,
                            int detectorColumns,
                            double detectorCenterRow,
                            double detectorCenterColumn,
                            int stampSize )
    : m_template( psfTemplate, stampSize ), m_detectorRows( detectorRows ), m_detectorColumns( detectorColumns ),
      m_detectorCenterRow( detectorCenterRow ), m_detectorCenterColumn( detectorCenterColumn ), m_stampSize( stampSize )
{
    if( detectorRows <= 0 || detectorColumns <= 0 || !std::isfinite( detectorCenterRow ) ||
        !std::isfinite( detectorCenterColumn ) )
    {
        throw std::invalid_argument( "KLIP PSF detector geometry must be positive and finite" );
    }
    if( stampSize <= 0 || stampSize % 2 == 0 )
    {
        throw std::invalid_argument( "KLIP PSF response stamp size must be positive and odd" );
    }
}

KLIPPSFModel::regionIndexT
KLIPPSFModel::regionIndex( const std::vector<std::size_t> &indices, int detectorRows, int detectorColumns )
{
    if( indices.empty() || detectorRows <= 0 || detectorColumns <= 0 )
    {
        throw std::invalid_argument( "KLIP PSF region geometry must be nonempty and positive" );
    }
    const std::size_t detectorPixels = static_cast<std::size_t>( detectorRows ) * detectorColumns;
    if( indices.size() > static_cast<std::size_t>( std::numeric_limits<int>::max() ) )
    {
        throw std::length_error( "KLIP PSF region exceeds integer lookup range" );
    }
    regionIndexT output = regionIndexT::Constant( detectorRows, detectorColumns, -1 );
    for( std::size_t regionPixel = 0; regionPixel < indices.size(); ++regionPixel )
    {
        if( indices[regionPixel] >= detectorPixels )
        {
            throw std::out_of_range( "KLIP PSF region index is outside the detector" );
        }
        const int row = static_cast<int>( indices[regionPixel] % static_cast<std::size_t>( detectorRows ) );
        const int column = static_cast<int>( indices[regionPixel] / static_cast<std::size_t>( detectorRows ) );
        if( output( row, column ) >= 0 )
        {
            throw std::invalid_argument( "KLIP PSF region indices must be unique" );
        }
        output( row, column ) = static_cast<int>( regionPixel );
    }
    return output;
}

void KLIPPSFModel::probe( vectorT &output,
                          const std::vector<std::size_t> &indices,
                          double sourceSkyRow,
                          double sourceSkyColumn,
                          double derotationAngle,
                          HCI::meanSub centering,
                          const imageT &cutoutMask ) const
{
    if( indices.empty() || indices.size() > static_cast<std::size_t>( std::numeric_limits<Eigen::Index>::max() ) ||
        !std::isfinite( sourceSkyRow ) || !std::isfinite( sourceSkyColumn ) || !std::isfinite( derotationAngle ) )
    {
        throw std::invalid_argument( "KLIP PSF probe geometry must be nonempty and finite" );
    }
    if( centering != HCI::meanSub::none && centering != HCI::meanSub::imageMean )
    {
        throw std::invalid_argument( "KLIP PSF probes support only none or imageMean centering" );
    }
    const bool masked = cutoutMask.size() != 0;
    if( masked && ( cutoutMask.rows() != static_cast<Eigen::Index>( indices.size() ) || cutoutMask.cols() != 1 ||
                    !klipPSFAllFinite( cutoutMask ) ) )
    {
        throw std::invalid_argument( "KLIP PSF cutout mask must be finite and match the region" );
    }

    const std::pair<double, double> sourceDetector = inverseRotate( sourceSkyRow, sourceSkyColumn, derotationAngle );
    const std::size_t detectorPixels = static_cast<std::size_t>( m_detectorRows ) * m_detectorColumns;
    output.resize( static_cast<Eigen::Index>( indices.size() ) );
    for( std::size_t regionPixel = 0; regionPixel < indices.size(); ++regionPixel )
    {
        if( indices[regionPixel] >= detectorPixels )
        {
            throw std::out_of_range( "KLIP PSF probe region index is outside the detector" );
        }
        const int row = static_cast<int>( indices[regionPixel] % static_cast<std::size_t>( m_detectorRows ) );
        const int column = static_cast<int>( indices[regionPixel] / static_cast<std::size_t>( m_detectorRows ) );
        output( static_cast<Eigen::Index>( regionPixel ) ) =
            m_template.sampleTemplate( static_cast<double>( row ) - sourceDetector.first,
                                       static_cast<double>( column ) - sourceDetector.second );
    }

    if( centering == HCI::meanSub::imageMean )
    {
        float mean{ 0 };
        if( masked )
        {
            const float maskSum = cutoutMask.sum();
            if( !std::isfinite( maskSum ) || maskSum <= 0 )
            {
                throw std::invalid_argument( "KLIP PSF cutout mask must select at least one pixel" );
            }
            mean = ( output * cutoutMask.col( 0 ) ).sum() / maskSum;
        }
        else
        {
            mean = output.mean();
        }
        output -= mean;
    }
    if( masked )
    {
        output *= cutoutMask.col( 0 );
    }
    if( !klipPSFAllFinite( output ) )
    {
        throw std::runtime_error( "KLIP PSF probe calculation produced nonfinite values" );
    }
}

void KLIPPSFModel::residuals( std::vector<vectorT> &outputs,
                              const vectorT &probe,
                              const imageT &klModes,
                              const std::vector<int> &modeCounts )
{
    if( probe.size() <= 0 || klModes.rows() <= 0 || klModes.cols() != probe.rows() || modeCounts.empty() ||
        !klipPSFAllFinite( probe ) || !klipPSFAllFinite( klModes ) )
    {
        throw std::invalid_argument( "KLIP PSF residual inputs have inconsistent or nonfinite dimensions" );
    }
    outputs.resize( modeCounts.size() );
    for( std::size_t outputIndex = 0; outputIndex < modeCounts.size(); ++outputIndex )
    {
        if( modeCounts[outputIndex] <= 0 )
        {
            throw std::invalid_argument( "KLIP PSF mode counts must be positive" );
        }
        outputs[outputIndex] = probe;
        const int retainedModes = std::min<int>( modeCounts[outputIndex], klModes.rows() );
        for( int offset = 0; offset < retainedModes; ++offset )
        {
            const int mode = klModes.rows() - 1 - offset;
            const float coefficient = klModes.row( mode ).matrix().dot( probe.matrix() );
            outputs[outputIndex] -= coefficient * klModes.row( mode ).transpose();
        }
        if( !klipPSFAllFinite( outputs[outputIndex] ) )
        {
            throw std::runtime_error( "KLIP PSF residual calculation produced nonfinite values" );
        }
    }
}

void KLIPPSFModel::projectedResidual( vectorT &output, const vectorT &probe, const imageT &projectionMatrix )
{
    if( probe.size() <= 0 || projectionMatrix.rows() != probe.rows() || projectionMatrix.cols() != probe.rows() ||
        !klipPSFAllFinite( probe ) || !klipPSFAllFinite( projectionMatrix ) )
    {
        throw std::invalid_argument( "KLIP PSF projection inputs have inconsistent or nonfinite dimensions" );
    }
    output = probe - ( projectionMatrix.matrix() * probe.matrix() ).array();
    if( !klipPSFAllFinite( output ) )
    {
        throw std::runtime_error( "KLIP PSF projected residual contains nonfinite values" );
    }
}

void KLIPPSFModel::accumulate( imageT &output,
                               validityT &outputValidity,
                               const vectorT &regionResponse,
                               const regionIndexT &regionIndices,
                               double sourceSkyRow,
                               double sourceSkyColumn,
                               double derotationAngle ) const
{
    if( regionResponse.size() <= 0 || regionIndices.rows() != m_detectorRows ||
        regionIndices.cols() != m_detectorColumns || !std::isfinite( sourceSkyRow ) ||
        !std::isfinite( sourceSkyColumn ) || !std::isfinite( derotationAngle ) || !klipPSFAllFinite( regionResponse ) )
    {
        throw std::invalid_argument( "KLIP PSF accumulation inputs are inconsistent or nonfinite" );
    }
    if( output.size() == 0 && outputValidity.size() == 0 )
    {
        output = imageT::Zero( m_stampSize, m_stampSize );
        outputValidity = validityT::Zero( m_stampSize, m_stampSize );
    }
    if( output.rows() != m_stampSize || output.cols() != m_stampSize || outputValidity.rows() != m_stampSize ||
        outputValidity.cols() != m_stampSize )
    {
        throw std::invalid_argument( "KLIP PSF accumulated stamp dimensions do not match the model" );
    }

    using gridT = P4PSFModel::gridT;
    gridT::transformT transform;
    const double stampCenter = 0.5 * static_cast<double>( m_stampSize - 1 );
    for( int stampColumn = 0; stampColumn < m_stampSize; ++stampColumn )
    {
        const double skyColumn = sourceSkyColumn + static_cast<double>( stampColumn ) - stampCenter;
        for( int stampRow = 0; stampRow < m_stampSize; ++stampRow )
        {
            const double skyRow = sourceSkyRow + static_cast<double>( stampRow ) - stampCenter;
            const std::pair<double, double> detector = inverseRotate( skyRow, skyColumn, derotationAngle );
            const int floorRow = static_cast<int>( std::floor( detector.first ) );
            const int floorColumn = static_cast<int>( std::floor( detector.second ) );
            const int upperRow = m_detectorRows - gridT::width + gridT::leftBuffer;
            const int upperColumn = m_detectorColumns - gridT::width + gridT::leftBuffer;
            if( floorRow <= gridT::leftBuffer || floorRow >= upperRow || floorColumn <= gridT::leftBuffer ||
                floorColumn >= upperColumn )
            {
                continue;
            }
            const int footprintRow = floorRow - gridT::leftBuffer;
            const int footprintColumn = floorColumn - gridT::leftBuffer;

            gridT::kernelT kernel;
            transform( kernel,
                       static_cast<float>( detector.first - std::floor( detector.first ) ),
                       static_cast<float>( detector.second - std::floor( detector.second ) ) );
            float contribution{ 0 };
            for( int columnOffset = 0; columnOffset < gridT::width; ++columnOffset )
            {
                for( int rowOffset = 0; rowOffset < gridT::width; ++rowOffset )
                {
                    const int regionPixel = regionIndices( footprintRow + rowOffset, footprintColumn + columnOffset );
                    if( regionPixel >= 0 )
                    {
                        if( regionPixel >= regionResponse.rows() )
                        {
                            throw std::invalid_argument( "KLIP PSF region lookup exceeds the response vector" );
                        }
                        contribution += regionResponse( regionPixel ) * kernel( rowOffset, columnOffset );
                    }
                }
            }
            output( stampRow, stampColumn ) += contribution;
            outputValidity( stampRow, stampColumn ) = 1;
        }
    }
}

int KLIPPSFModel::stampSize() const noexcept
{
    return m_stampSize;
}

std::pair<double, double> KLIPPSFModel::inverseRotate( double row, double column, double angle ) const
{
    const double deltaRow = row - m_detectorCenterRow;
    const double deltaColumn = column - m_detectorCenterColumn;
    return { m_detectorCenterRow + deltaRow * std::cos( angle ) + deltaColumn * std::sin( angle ),
             m_detectorCenterColumn - deltaRow * std::sin( angle ) + deltaColumn * std::cos( angle ) };
}

} // namespace improc
} // namespace mx
