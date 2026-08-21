/** \file P4PSFReconstructor.cpp
 * \brief Implements compact detector-to-sky reconstruction for P4 frozen-model PSFs.
 * \author Jared R. Males
 */

#include "P4PSFReconstructor.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include <mx/math/floatUtils.hpp>

namespace mx
{
namespace improc
{

P4PSFReconstructor::P4PSFReconstructor( int detectorRows,
                                        int detectorColumns,
                                        double detectorCenterRow,
                                        double detectorCenterColumn,
                                        int outputStampSize,
                                        int localStampRows,
                                        int localStampColumns )
    : m_detectorRows( detectorRows ), m_detectorColumns( detectorColumns ), m_detectorCenterRow( detectorCenterRow ),
      m_detectorCenterColumn( detectorCenterColumn ), m_outputStampSize( outputStampSize ),
      m_localStampRows( localStampRows ), m_localStampColumns( localStampColumns )
{
    if( detectorRows <= 0 || detectorColumns <= 0 )
    {
        throw std::invalid_argument( "P4 PSF reconstruction detector dimensions must be positive" );
    }
    if( !mx::math::isFinite( detectorCenterRow ) || !mx::math::isFinite( detectorCenterColumn ) )
    {
        throw std::invalid_argument( "P4 PSF reconstruction detector center must be finite" );
    }
    if( outputStampSize <= 0 || localStampRows <= 0 || localStampColumns <= 0 )
    {
        throw std::invalid_argument( "P4 PSF reconstruction stamp dimensions must be positive" );
    }
    const long double localPixels =
        static_cast<long double>( localStampRows ) * static_cast<long double>( localStampColumns );
    if( localPixels > static_cast<long double>( std::numeric_limits<Eigen::Index>::max() ) )
    {
        throw std::length_error( "P4 local PSF stamp exceeds Eigen index range" );
    }
}

std::pair<double, double> P4PSFReconstructor::inverseRotate( double row, double column, double angle ) const
{
    const double deltaRow = row - m_detectorCenterRow;
    const double deltaColumn = column - m_detectorCenterColumn;
    return { m_detectorCenterRow + deltaRow * std::cos( angle ) + deltaColumn * std::sin( angle ),
             m_detectorCenterColumn - deltaRow * std::sin( angle ) + deltaColumn * std::cos( angle ) };
}

bool P4PSFReconstructor::sampleLocalResponse(
    float &value, const imageT &localModels, Eigen::Index modelColumn, double deltaRow, double deltaColumn ) const
{
    const double sampleRow = 0.5 * static_cast<double>( m_localStampRows - 1 ) + deltaRow;
    const double sampleColumn = 0.5 * static_cast<double>( m_localStampColumns - 1 ) + deltaColumn;
    if( !mx::math::isFinite( sampleRow ) || !mx::math::isFinite( sampleColumn ) )
    {
        return false;
    }

    const double floorRow = std::floor( sampleRow );
    const double floorColumn = std::floor( sampleColumn );
    constexpr int rightBuffer = gridT::width - gridT::leftBuffer - 1;
    if( floorRow < static_cast<double>( std::numeric_limits<int>::min() + gridT::leftBuffer ) ||
        floorRow > static_cast<double>( std::numeric_limits<int>::max() - rightBuffer ) ||
        floorColumn < static_cast<double>( std::numeric_limits<int>::min() + gridT::leftBuffer ) ||
        floorColumn > static_cast<double>( std::numeric_limits<int>::max() - rightBuffer ) )
    {
        return false;
    }

    const int footprintRow = static_cast<int>( floorRow ) - gridT::leftBuffer;
    const int footprintColumn = static_cast<int>( floorColumn ) - gridT::leftBuffer;
    if( footprintRow < 0 || footprintColumn < 0 || footprintRow + gridT::width > m_localStampRows ||
        footprintColumn + gridT::width > m_localStampColumns )
    {
        return false;
    }

    gridT::transformT transform;
    gridT::kernelT kernel;
    transform( kernel, static_cast<float>( sampleRow - floorRow ), static_cast<float>( sampleColumn - floorColumn ) );

    float sample{ 0 };
    for( int columnOffset = 0; columnOffset < gridT::width; ++columnOffset )
    {
        for( int rowOffset = 0; rowOffset < gridT::width; ++rowOffset )
        {
            const Eigen::Index stampPixel = static_cast<Eigen::Index>( footprintRow + rowOffset ) +
                                            static_cast<Eigen::Index>( m_localStampRows ) *
                                                static_cast<Eigen::Index>( footprintColumn + columnOffset );
            sample += localModels( stampPixel, modelColumn ) * kernel( rowOffset, columnOffset );
        }
    }
    value = sample;
    return true;
}

void P4PSFReconstructor::reconstructFrame( imageT &output,
                                           validityT &outputValidity,
                                           const imageT &localModels,
                                           const validityT &localValidity,
                                           const searchIndexT &searchIndex,
                                           std::size_t modeIndex,
                                           double sourceSkyRow,
                                           double sourceSkyColumn,
                                           double derotationAngle ) const
{
    const Eigen::Index localPixels =
        static_cast<Eigen::Index>( m_localStampRows ) * static_cast<Eigen::Index>( m_localStampColumns );
    if( localModels.rows() != localPixels || localValidity.rows() <= 0 || localValidity.cols() <= 0 )
    {
        throw std::invalid_argument( "P4 compact local PSF dimensions are inconsistent with reconstruction geometry" );
    }
    if( modeIndex >= static_cast<std::size_t>( localValidity.cols() ) )
    {
        throw std::out_of_range( "P4 PSF reconstruction mode index is out of range" );
    }
    const long double expectedColumns =
        static_cast<long double>( localValidity.rows() ) * static_cast<long double>( localValidity.cols() );
    if( expectedColumns != static_cast<long double>( localModels.cols() ) )
    {
        throw std::invalid_argument( "P4 compact local PSF column count does not match validity dimensions" );
    }
    if( searchIndex.rows() != m_detectorRows || searchIndex.cols() != m_detectorColumns )
    {
        throw std::invalid_argument( "P4 PSF search-index image does not match detector dimensions" );
    }
    if( !mx::math::isFinite( sourceSkyRow ) || !mx::math::isFinite( sourceSkyColumn ) ||
        !mx::math::isFinite( derotationAngle ) )
    {
        throw std::invalid_argument( "P4 PSF reconstruction coordinates and angle must be finite" );
    }

    output.resize( m_outputStampSize, m_outputStampSize );
    output.setZero();
    outputValidity.resize( m_outputStampSize, m_outputStampSize );
    outputValidity.setZero();

    const std::pair<double, double> sourceDetector = inverseRotate( sourceSkyRow, sourceSkyColumn, derotationAngle );
    const double outputCenter = 0.5 * static_cast<double>( m_outputStampSize - 1 );
    int minimumFootprintRow = m_detectorRows;
    int minimumFootprintColumn = m_detectorColumns;
    int maximumFootprintRow = -1;
    int maximumFootprintColumn = -1;
    for( int stampColumn = 0; stampColumn < m_outputStampSize; ++stampColumn )
    {
        const double skyColumn = sourceSkyColumn + static_cast<double>( stampColumn ) - outputCenter;
        for( int stampRow = 0; stampRow < m_outputStampSize; ++stampRow )
        {
            const double skyRow = sourceSkyRow + static_cast<double>( stampRow ) - outputCenter;
            const std::pair<double, double> detector = inverseRotate( skyRow, skyColumn, derotationAngle );
            const double floorRow = std::floor( detector.first );
            const double floorColumn = std::floor( detector.second );
            if( floorRow < 0 || floorColumn < 0 || floorRow > m_detectorRows - 1 ||
                floorColumn > m_detectorColumns - 1 )
            {
                continue;
            }
            const int footprintRow = static_cast<int>( floorRow ) - gridT::leftBuffer;
            const int footprintColumn = static_cast<int>( floorColumn ) - gridT::leftBuffer;
            if( footprintRow < 0 || footprintColumn < 0 || footprintRow + gridT::width > m_detectorRows ||
                footprintColumn + gridT::width > m_detectorColumns )
            {
                continue;
            }
            minimumFootprintRow = std::min( minimumFootprintRow, footprintRow );
            minimumFootprintColumn = std::min( minimumFootprintColumn, footprintColumn );
            maximumFootprintRow = std::max( maximumFootprintRow, footprintRow + gridT::width - 1 );
            maximumFootprintColumn = std::max( maximumFootprintColumn, footprintColumn + gridT::width - 1 );
        }
    }

    if( maximumFootprintRow < minimumFootprintRow || maximumFootprintColumn < minimumFootprintColumn )
    {
        return;
    }

    const int cacheRows = maximumFootprintRow - minimumFootprintRow + 1;
    const int cacheColumns = maximumFootprintColumn - minimumFootprintColumn + 1;
    imageT detectorResponse = imageT::Zero( cacheRows, cacheColumns );
    validityT detectorValidity = validityT::Zero( cacheRows, cacheColumns );
    for( int detectorColumn = minimumFootprintColumn; detectorColumn <= maximumFootprintColumn; ++detectorColumn )
    {
        for( int detectorRow = minimumFootprintRow; detectorRow <= maximumFootprintRow; ++detectorRow )
        {
            const int search = searchIndex( detectorRow, detectorColumn );
            if( search < 0 || search >= localValidity.rows() ||
                localValidity( search, static_cast<Eigen::Index>( modeIndex ) ) == 0 )
            {
                continue;
            }
            const Eigen::Index modelColumn =
                static_cast<Eigen::Index>( search ) * localValidity.cols() + static_cast<Eigen::Index>( modeIndex );
            float localSample{ 0 };
            if( sampleLocalResponse( localSample,
                                     localModels,
                                     modelColumn,
                                     static_cast<double>( detectorRow ) - sourceDetector.first,
                                     static_cast<double>( detectorColumn ) - sourceDetector.second ) )
            {
                detectorResponse( detectorRow - minimumFootprintRow, detectorColumn - minimumFootprintColumn ) =
                    localSample;
                detectorValidity( detectorRow - minimumFootprintRow, detectorColumn - minimumFootprintColumn ) = 1;
            }
        }
    }

    gridT::transformT transform;
    for( int stampColumn = 0; stampColumn < m_outputStampSize; ++stampColumn )
    {
        const double skyColumn = sourceSkyColumn + static_cast<double>( stampColumn ) - outputCenter;
        for( int stampRow = 0; stampRow < m_outputStampSize; ++stampRow )
        {
            const double skyRow = sourceSkyRow + static_cast<double>( stampRow ) - outputCenter;
            const std::pair<double, double> detector = inverseRotate( skyRow, skyColumn, derotationAngle );
            const double floorRow = std::floor( detector.first );
            const double floorColumn = std::floor( detector.second );
            if( floorRow < 0 || floorColumn < 0 || floorRow > m_detectorRows - 1 ||
                floorColumn > m_detectorColumns - 1 )
            {
                continue;
            }
            const int footprintRow = static_cast<int>( floorRow ) - gridT::leftBuffer;
            const int footprintColumn = static_cast<int>( floorColumn ) - gridT::leftBuffer;
            if( footprintRow < minimumFootprintRow || footprintColumn < minimumFootprintColumn ||
                footprintRow + gridT::width - 1 > maximumFootprintRow ||
                footprintColumn + gridT::width - 1 > maximumFootprintColumn )
            {
                continue;
            }

            gridT::kernelT kernel;
            transform( kernel,
                       static_cast<float>( detector.first - floorRow ),
                       static_cast<float>( detector.second - floorColumn ) );
            bool valid{ true };
            float response{ 0 };
            for( int columnOffset = 0; columnOffset < gridT::width && valid; ++columnOffset )
            {
                for( int rowOffset = 0; rowOffset < gridT::width; ++rowOffset )
                {
                    const int cacheRow = footprintRow + rowOffset - minimumFootprintRow;
                    const int cacheColumn = footprintColumn + columnOffset - minimumFootprintColumn;
                    if( detectorValidity( cacheRow, cacheColumn ) == 0 )
                    {
                        valid = false;
                        break;
                    }
                    response += detectorResponse( cacheRow, cacheColumn ) * kernel( rowOffset, columnOffset );
                }
            }
            if( valid )
            {
                output( stampRow, stampColumn ) = response;
                outputValidity( stampRow, stampColumn ) = 1;
            }
        }
    }
}

void P4PSFReconstructor::reconstructFrames( cubeT &output,
                                            cubeT &outputValidity,
                                            const imageT &localModels,
                                            const validityT &localValidity,
                                            const searchIndexT &searchIndex,
                                            std::size_t modeIndex,
                                            double sourceSkyRow,
                                            double sourceSkyColumn,
                                            const std::vector<double> &derotationAngles ) const
{
    if( derotationAngles.empty() )
    {
        throw std::invalid_argument( "P4 PSF reconstruction requires at least one derotation angle" );
    }
    if( derotationAngles.size() > static_cast<std::size_t>( std::numeric_limits<int>::max() ) )
    {
        throw std::length_error( "P4 PSF reconstruction frame count exceeds cube range" );
    }

    output.resize( m_outputStampSize, m_outputStampSize, static_cast<int>( derotationAngles.size() ) );
    output.setZero();
    outputValidity.resize( m_outputStampSize, m_outputStampSize, static_cast<int>( derotationAngles.size() ) );
    outputValidity.setZero();
    imageT frame;
    validityT validity;
    for( std::size_t image = 0; image < derotationAngles.size(); ++image )
    {
        reconstructFrame( frame,
                          validity,
                          localModels,
                          localValidity,
                          searchIndex,
                          modeIndex,
                          sourceSkyRow,
                          sourceSkyColumn,
                          derotationAngles[image] );
        output.image( static_cast<int>( image ) ) = frame;
        outputValidity.image( static_cast<int>( image ) ) = validity.cast<float>();
    }
}

void P4PSFReconstructor::combineFrames( imageT &output,
                                        validityT &outputValidity,
                                        cubeT &frames,
                                        cubeT &frameValidity,
                                        HCI::combine method,
                                        std::vector<float> weights,
                                        float sigmaThreshold,
                                        float minimumGoodFraction )
{
    if( frames.rows() <= 0 || frames.cols() <= 0 || frames.planes() <= 0 || frameValidity.rows() != frames.rows() ||
        frameValidity.cols() != frames.cols() || frameValidity.planes() != frames.planes() )
    {
        throw std::invalid_argument( "P4 PSF response and validity cubes must have matching nonempty dimensions" );
    }
    if( !weights.empty() && weights.size() != static_cast<std::size_t>( frames.planes() ) )
    {
        throw std::invalid_argument( "P4 PSF combination weight count must match the frame count" );
    }
    if( !mx::math::isFinite( sigmaThreshold ) || !mx::math::isFinite( minimumGoodFraction ) ||
        minimumGoodFraction < 0 || minimumGoodFraction > 1 )
    {
        throw std::invalid_argument( "P4 PSF combination thresholds must be finite and minimum-good in [0,1]" );
    }
    for( const float weight : weights )
    {
        if( !mx::math::isFinite( weight ) )
        {
            throw std::invalid_argument( "P4 PSF combination weights must be finite" );
        }
    }

    if( method == HCI::combine::median )
    {
        frames.median( output, frameValidity, minimumGoodFraction );
    }
    else if( method == HCI::combine::mean || ( method == HCI::combine::sigmaMean && sigmaThreshold <= 0 ) )
    {
        if( weights.empty() )
        {
            frames.mean( output, frameValidity, minimumGoodFraction );
        }
        else
        {
            frames.mean( output, weights, frameValidity, minimumGoodFraction );
        }
    }
    else if( method == HCI::combine::sigmaMean )
    {
        if( weights.empty() )
        {
            frames.sigmaMean( output, frameValidity, sigmaThreshold, minimumGoodFraction );
        }
        else
        {
            frames.sigmaMean( output, weights, frameValidity, sigmaThreshold, minimumGoodFraction );
        }
    }
    else if( method == HCI::combine::none )
    {
        throw std::invalid_argument( "P4 combined PSF reconstruction requires a final combination method" );
    }
    else
    {
        throw std::invalid_argument( "P4 PSF reconstruction received an unsupported combination method" );
    }

    outputValidity.resize( output.rows(), output.cols() );
    for( Eigen::Index column = 0; column < output.cols(); ++column )
    {
        for( Eigen::Index row = 0; row < output.rows(); ++row )
        {
            if( mx::math::isFinite( output( row, column ) ) )
            {
                outputValidity( row, column ) = 1;
            }
            else
            {
                output( row, column ) = 0;
                outputValidity( row, column ) = 0;
            }
        }
    }
}

void P4PSFReconstructor::reconstructCombined( imageT &output,
                                              validityT &outputValidity,
                                              const imageT &localModels,
                                              const validityT &localValidity,
                                              const searchIndexT &searchIndex,
                                              std::size_t modeIndex,
                                              double sourceSkyRow,
                                              double sourceSkyColumn,
                                              const std::vector<double> &derotationAngles,
                                              HCI::combine method,
                                              const std::vector<float> &weights,
                                              float sigmaThreshold,
                                              float minimumGoodFraction ) const
{
    cubeT frames;
    cubeT frameValidity;
    reconstructFrames( frames,
                       frameValidity,
                       localModels,
                       localValidity,
                       searchIndex,
                       modeIndex,
                       sourceSkyRow,
                       sourceSkyColumn,
                       derotationAngles );
    combineFrames( output,
                   outputValidity,
                   frames,
                   frameValidity,
                   method,
                   weights,
                   sigmaThreshold,
                   minimumGoodFraction );
}

} // namespace improc
} // namespace mx
