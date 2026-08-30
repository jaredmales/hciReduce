/** \file P4PSFReconstructor_test.cpp
 * \brief Tests compact detector-to-sky reconstruction for P4 frozen-model PSFs.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "src/common/P4PSFModel.hpp"
#include "src/common/P4PSFReconstructor.hpp"

#include <cmath>
#include <limits>
#include <vector>

namespace unitTest
{
namespace P4PSFReconstructor_test
{

/** \cond P4PSFReconstructor_test_harness */
using modelT = mx::improc::P4PSFModel;
using reconstructorT = mx::improc::P4PSFReconstructor;
using gridT = modelT::gridT;
using imageT = modelT::imageT;
using coefficientT = modelT::coefficientT;
using validityT = reconstructorT::validityT;
using searchIndexT = reconstructorT::searchIndexT;

/// Return a broad annulus whose search pixels exercise spatially variable reconstruction.
mx::improc::P4PixelGridRegion testRegion()
{
    return mx::improc::P4PixelGridRegion( 2.0,
                                          14.0,
                                          2.0,
                                          2.0,
                                          3.0,
                                          50.0,
                                          0.5,
                                          mx::improc::P4ExclusionPolicy::sampleCenter,
                                          0.0 );
}

/// Construct a finite asymmetric centered test template.
imageT asymmetricTemplate( int rows, /**< [in] template rows */
                           int columns /**< [in] template columns */ )
{
    imageT psf( rows, columns );
    const double centerRow = 0.5 * static_cast<double>( rows - 1 );
    const double centerColumn = 0.5 * static_cast<double>( columns - 1 );
    for( int column = 0; column < columns; ++column )
    {
        for( int row = 0; row < rows; ++row )
        {
            const double deltaRow = static_cast<double>( row ) - centerRow;
            const double deltaColumn = static_cast<double>( column ) - centerColumn;
            psf( row, column ) =
                static_cast<float>( std::exp( -0.11 * deltaRow * deltaRow - 0.06 * deltaColumn * deltaColumn ) *
                                    ( 1.0 + 0.025 * deltaRow - 0.018 * deltaColumn + 0.003 * deltaRow * deltaColumn ) );
        }
    }
    return psf;
}

/// Cubically sample an image with zero exterior padding.
float sampleImageZero( const imageT &image, /**< [in] finite image */
                       double row,          /**< [in] sample row */
                       double column /**< [in] sample column */ )
{
    const double floorRow = std::floor( row );
    const double floorColumn = std::floor( column );
    gridT::transformT transform;
    gridT::kernelT kernel;
    transform( kernel, static_cast<float>( row - floorRow ), static_cast<float>( column - floorColumn ) );
    const int footprintRow = static_cast<int>( floorRow ) - gridT::leftBuffer;
    const int footprintColumn = static_cast<int>( floorColumn ) - gridT::leftBuffer;
    float value{ 0 };
    for( int columnOffset = 0; columnOffset < gridT::width; ++columnOffset )
    {
        for( int rowOffset = 0; rowOffset < gridT::width; ++rowOffset )
        {
            const int imageRow = footprintRow + rowOffset;
            const int imageColumn = footprintColumn + columnOffset;
            if( imageRow >= 0 && imageRow < image.rows() && imageColumn >= 0 && imageColumn < image.cols() )
            {
                value += image( imageRow, imageColumn ) * kernel( rowOffset, columnOffset );
            }
        }
    }
    return value;
}

/// Map one sky coordinate to the detector coordinate used by production image rotation.
std::pair<double, double> inverseRotate( double row,          /**< [in] sky row */
                                         double column,       /**< [in] sky column */
                                         double centerRow,    /**< [in] rotation-center row */
                                         double centerColumn, /**< [in] rotation-center column */
                                         double angle /**< [in] counterclockwise output rotation in radians */ )
{
    const double deltaRow = row - centerRow;
    const double deltaColumn = column - centerColumn;
    return { centerRow + deltaRow * std::cos( angle ) + deltaColumn * std::sin( angle ),
             centerColumn - deltaRow * std::sin( angle ) + deltaColumn * std::cos( angle ) };
}

/// Construct a detector-sampled image containing one shifted template.
imageT injectedDetectorImage( const imageT &psfTemplate, /**< [in] centered template */
                              int detectorRows,          /**< [in] detector rows */
                              int detectorColumns,       /**< [in] detector columns */
                              double sourceRow,          /**< [in] detector source row */
                              double sourceColumn /**< [in] detector source column */ )
{
    imageT detector( detectorRows, detectorColumns );
    const double templateCenterRow = 0.5 * static_cast<double>( psfTemplate.rows() - 1 );
    const double templateCenterColumn = 0.5 * static_cast<double>( psfTemplate.cols() - 1 );
    for( int column = 0; column < detectorColumns; ++column )
    {
        for( int row = 0; row < detectorRows; ++row )
        {
            detector( row, column ) =
                sampleImageZero( psfTemplate,
                                 templateCenterRow + static_cast<double>( row ) - sourceRow,
                                 templateCenterColumn + static_cast<double>( column ) - sourceColumn );
        }
    }
    return detector;
}

/// Apply explicitly known spatially varying P4 coefficients to a full injected detector image.
imageT directDetectorResponse( const imageT &psfTemplate,                     /**< [in] centered template */
                               const gridT &grid,                             /**< [in] complete P4 geometry */
                               const std::vector<coefficientT> &coefficients, /**< [in] coefficients by search pixel */
                               double sourceRow,                              /**< [in] detector source row */
                               double sourceColumn /**< [in] detector source column */ )
{
    const imageT detector = injectedDetectorImage( psfTemplate, grid.rows(), grid.columns(), sourceRow, sourceColumn );
    imageT response = imageT::Zero( grid.rows(), grid.columns() );
    for( std::size_t search = 0; search < grid.searchPixelCount(); ++search )
    {
        if( !grid.searchPixel( search ).valid() )
        {
            continue;
        }
        const mx::improc::P4PixelCoordinate &coordinate = grid.searchPixel( search ).coordinate();
        double value = detector( coordinate.row(), coordinate.column() );
        for( std::size_t predictor = 0; predictor < grid.predictorCount(); ++predictor )
        {
            value -= coefficients[search]( static_cast<Eigen::Index>( predictor ) ) *
                     grid.sample( detector, search, predictor );
        }
        response( coordinate.row(), coordinate.column() ) = static_cast<float>( value );
    }
    return response;
}

/// Prepare compact local models and lookup state for one test mode.
void prepareCompactModels( imageT &localModels,                     /**< [out] flattened compact response columns */
                           validityT &localValidity,                /**< [out] search/mode validity */
                           searchIndexT &searchIndex,               /**< [out] detector search lookup */
                           std::vector<coefficientT> &coefficients, /**< [out] coefficients by search pixel */
                           const modelT &model,                     /**< [in] prepared local-model calculator */
                           const gridT &grid /**< [in] complete P4 geometry */ )
{
    localModels.resize( static_cast<Eigen::Index>( model.stampRows() ) * model.stampColumns(),
                        static_cast<Eigen::Index>( grid.searchPixelCount() ) );
    localModels.setZero();
    localValidity.resize( static_cast<Eigen::Index>( grid.searchPixelCount() ), 1 );
    localValidity.setZero();
    searchIndex = searchIndexT::Constant( grid.rows(), grid.columns(), -1 );
    coefficients.clear();
    coefficients.reserve( grid.searchPixelCount() );
    imageT localResponse;
    for( std::size_t search = 0; search < grid.searchPixelCount(); ++search )
    {
        coefficientT coefficient( static_cast<Eigen::Index>( grid.predictorCount() ) );
        for( Eigen::Index predictor = 0; predictor < coefficient.rows(); ++predictor )
        {
            coefficient( predictor ) = 0.0025 * std::sin( 0.09 * static_cast<double>( predictor + 2 * search ) ) +
                                       0.0002 * std::cos( 0.13 * static_cast<double>( 3 * predictor + search ) );
        }
        coefficients.push_back( coefficient );
        if( !grid.searchPixel( search ).valid() )
        {
            continue;
        }
        model.calculateLocalResponse( localResponse, grid, search, coefficients.back() );
        const Eigen::Index column = static_cast<Eigen::Index>( search );
        for( int stampColumn = 0; stampColumn < model.stampColumns(); ++stampColumn )
        {
            for( int stampRow = 0; stampRow < model.stampRows(); ++stampRow )
            {
                localModels( stampRow + model.stampRows() * stampColumn, column ) =
                    localResponse( stampRow, stampColumn );
            }
        }
        localValidity( column, 0 ) = 1;
        const mx::improc::P4PixelCoordinate &coordinate = grid.searchPixel( search ).coordinate();
        searchIndex( coordinate.row(), coordinate.column() ) = static_cast<int>( search );
    }
}

/** \endcond */

/// Verify compact reconstruction matches a direct frozen detector injection.
/** This exercises mx::improc::P4PSFReconstructor::reconstructFrame() for zero and nonzero angles, fractional inverse
 * coordinates, spatially varying coefficients, asymmetric images/templates, and odd/even final stamps.
 * \ingroup P4PSFReconstructor_unit_tests
 */
TEST_CASE( "P4 compact sky reconstruction matches direct frozen injection", "[P4PSFReconstructor][oracle][derotation]" )
{
    constexpr int detectorRows = 43;
    constexpr int detectorColumns = 45;
    constexpr double centerRow = 21.0;
    constexpr double centerColumn = 22.0;
    gridT grid;
    grid.resize( detectorRows, detectorColumns, centerRow, centerColumn );
    grid.region( testRegion(), nullptr );

    const int outputSize = GENERATE( 3, 4 );
    const double angle = GENERATE( 0.0, 0.23, -0.31 );
    const imageT psfTemplate = asymmetricTemplate( 15, 16 );
    constexpr int localRows = 19;
    constexpr int localColumns = 20;
    const modelT model( psfTemplate, localRows, localColumns );
    imageT localModels;
    validityT localValidity;
    searchIndexT searchIndex;
    std::vector<coefficientT> coefficients;
    prepareCompactModels( localModels, localValidity, searchIndex, coefficients, model, grid );

    constexpr double sourceSkyRow = 21.0;
    constexpr double sourceSkyColumn = 29.0;
    const std::pair<double, double> sourceDetector =
        inverseRotate( sourceSkyRow, sourceSkyColumn, centerRow, centerColumn, angle );
    const imageT directDetector =
        directDetectorResponse( psfTemplate, grid, coefficients, sourceDetector.first, sourceDetector.second );

    reconstructorT
        reconstructor( detectorRows, detectorColumns, centerRow, centerColumn, outputSize, localRows, localColumns );
    imageT compact;
    validityT compactValidity;
    reconstructor.reconstructFrame( compact,
                                    compactValidity,
                                    localModels,
                                    localValidity,
                                    searchIndex,
                                    0,
                                    sourceSkyRow,
                                    sourceSkyColumn,
                                    angle );

    const double outputCenter = 0.5 * static_cast<double>( outputSize - 1 );
    std::size_t validCount{ 0 };
    for( int stampColumn = 0; stampColumn < outputSize; ++stampColumn )
    {
        for( int stampRow = 0; stampRow < outputSize; ++stampRow )
        {
            const std::pair<double, double> detector = inverseRotate( sourceSkyRow + stampRow - outputCenter,
                                                                      sourceSkyColumn + stampColumn - outputCenter,
                                                                      centerRow,
                                                                      centerColumn,
                                                                      angle );
            const int footprintRow = static_cast<int>( std::floor( detector.first ) ) - gridT::leftBuffer;
            const int footprintColumn = static_cast<int>( std::floor( detector.second ) ) - gridT::leftBuffer;
            bool expectedValid = footprintRow >= 0 && footprintColumn >= 0 &&
                                 footprintRow + gridT::width <= detectorRows &&
                                 footprintColumn + gridT::width <= detectorColumns;
            for( int columnOffset = 0; columnOffset < gridT::width && expectedValid; ++columnOffset )
            {
                for( int rowOffset = 0; rowOffset < gridT::width; ++rowOffset )
                {
                    if( searchIndex( footprintRow + rowOffset, footprintColumn + columnOffset ) < 0 )
                    {
                        expectedValid = false;
                        break;
                    }
                }
            }
            REQUIRE( compactValidity( stampRow, stampColumn ) == static_cast<std::uint8_t>( expectedValid ) );
            if( expectedValid )
            {
                ++validCount;
                REQUIRE( compact( stampRow, stampColumn ) ==
                         Approx( sampleImageZero( directDetector, detector.first, detector.second ) ).margin( 1e-6 ) );
            }
            else
            {
                REQUIRE( compact( stampRow, stampColumn ) == 0 );
            }
        }
    }
    REQUIRE( validCount == static_cast<std::size_t>( outputSize * outputSize ) );
}

/// Verify temporal response components follow the selected physical image in every central frame.
/** This exercises mx::improc::P4PSFModel::calculateLocalResponse(), mx::improc::P4PSFModel::sampleTemplate(), and
 * mx::improc::P4PSFReconstructor::reconstructCombinedTemporal() against direct frozen-model detector cubes,
 * including one-sided endpoint selections and rotations beyond the compact same-image stamp.
 * \ingroup P4PSFReconstructor_unit_tests
 */
TEST_CASE( "P4 temporal compact reconstruction matches direct selected-image injection",
           "[P4PSFReconstructor][temporal][oracle]" )
{
    constexpr int detectorRows = 43;
    constexpr int detectorColumns = 45;
    constexpr double centerRow = 21.0;
    constexpr double centerColumn = 22.0;
    constexpr int outputSize = 3;
    constexpr int localRows = 19;
    constexpr int localColumns = 20;
    constexpr double sourceSkyRow = 21.0;
    constexpr double sourceSkyColumn = 29.0;
    const std::vector<double> angles{ 0.0, 1.6, -1.2, 2.4 };
    const reconstructorT::temporalSelectionT selections{ {
        { 0, 1, 2 },
        { 1, 0, 2 },
        { 2, 1, 3 },
        { 3, 2, 1 },
    } };
    const std::vector<mx::improc::P4PixelCoordinate> temporalOffsets{ mx::improc::P4PixelCoordinate( -1, 0 ),
                                                                      mx::improc::P4PixelCoordinate( 0, 0 ),
                                                                      mx::improc::P4PixelCoordinate( 1, 0 ) };
    constexpr std::size_t componentCount = 3;

    gridT grid;
    grid.resize( detectorRows, detectorColumns, centerRow, centerColumn );
    grid.region( testRegion(), nullptr );
    const imageT psfTemplate = asymmetricTemplate( 25, 26 );
    const modelT model( psfTemplate, localRows, localColumns );
    const Eigen::Index localPixels = static_cast<Eigen::Index>( localRows * localColumns );
    imageT localModels( localPixels, static_cast<Eigen::Index>( grid.searchPixelCount() ) );
    localModels.setZero();
    imageT temporalCoefficients( static_cast<Eigen::Index>( ( componentCount - 1 ) * temporalOffsets.size() ),
                                 static_cast<Eigen::Index>( grid.searchPixelCount() ) );
    temporalCoefficients.setZero();
    validityT localValidity( static_cast<Eigen::Index>( grid.searchPixelCount() ), 1 );
    localValidity.setZero();
    searchIndexT searchIndex = searchIndexT::Constant( detectorRows, detectorColumns, -1 );
    std::vector<coefficientT> coefficients;
    coefficients.reserve( grid.searchPixelCount() );
    imageT localResponse;
    for( std::size_t search = 0; search < grid.searchPixelCount(); ++search )
    {
        coefficientT coefficient(
            static_cast<Eigen::Index>( grid.predictorCount() + ( componentCount - 1 ) * temporalOffsets.size() ) );
        for( Eigen::Index predictor = 0; predictor < coefficient.rows(); ++predictor )
        {
            coefficient( predictor ) = 0.0017 * std::sin( 0.11 * static_cast<double>( predictor + 3 * search ) );
        }
        coefficients.push_back( coefficient );
        if( !grid.searchPixel( search ).valid() )
        {
            continue;
        }
        model.calculateLocalResponse( localResponse,
                                      grid,
                                      search,
                                      coefficients.back().head( static_cast<Eigen::Index>( grid.predictorCount() ) ) );
        for( int column = 0; column < localColumns; ++column )
        {
            for( int row = 0; row < localRows; ++row )
            {
                localModels( static_cast<Eigen::Index>( row + localRows * column ),
                             static_cast<Eigen::Index>( search ) ) = localResponse( row, column );
            }
        }
        for( Eigen::Index coefficient = static_cast<Eigen::Index>( grid.predictorCount() );
             coefficient < coefficients.back().rows();
             ++coefficient )
        {
            temporalCoefficients( coefficient - static_cast<Eigen::Index>( grid.predictorCount() ),
                                  static_cast<Eigen::Index>( search ) ) =
                static_cast<float>( -coefficients.back()( coefficient ) );
        }
        localValidity( static_cast<Eigen::Index>( search ), 0 ) = 1;
        const mx::improc::P4PixelCoordinate &coordinate = grid.searchPixel( search ).coordinate();
        searchIndex( coordinate.row(), coordinate.column() ) = static_cast<int>( search );
    }

    reconstructorT::cubeT directFrames;
    reconstructorT::cubeT directValidity;
    directFrames.resize( outputSize, outputSize, static_cast<int>( angles.size() ) );
    directValidity.resize( outputSize, outputSize, static_cast<int>( angles.size() ) );
    for( int image = 0; image < directValidity.planes(); ++image )
    {
        directValidity.image( image ).setOnes();
    }
    for( std::size_t central = 0; central < angles.size(); ++central )
    {
        std::vector<imageT> injected;
        for( const int selected : selections[0][central] )
        {
            const std::pair<double, double> source =
                inverseRotate( sourceSkyRow, sourceSkyColumn, centerRow, centerColumn, angles[selected] );
            injected.push_back(
                injectedDetectorImage( psfTemplate, detectorRows, detectorColumns, source.first, source.second ) );
        }
        imageT detectorResponse = imageT::Zero( detectorRows, detectorColumns );
        for( std::size_t search = 0; search < grid.searchPixelCount(); ++search )
        {
            if( !grid.searchPixel( search ).valid() )
            {
                continue;
            }
            const mx::improc::P4PixelCoordinate &coordinate = grid.searchPixel( search ).coordinate();
            double response = injected[0]( coordinate.row(), coordinate.column() );
            for( std::size_t predictor = 0; predictor < grid.predictorCount(); ++predictor )
            {
                response -= coefficients[search]( static_cast<Eigen::Index>( predictor ) ) *
                            grid.sample( injected[0], search, predictor );
            }
            for( std::size_t component = 1; component < componentCount; ++component )
            {
                for( std::size_t predictor = 0; predictor < temporalOffsets.size(); ++predictor )
                {
                    const mx::improc::P4PixelCoordinate &offset = temporalOffsets[predictor];
                    const Eigen::Index coefficientIndex = static_cast<Eigen::Index>(
                        grid.predictorCount() + ( component - 1 ) * temporalOffsets.size() + predictor );
                    response -=
                        coefficients[search]( coefficientIndex ) *
                        injected[component]( coordinate.row() + offset.row(), coordinate.column() + offset.column() );
                }
            }
            detectorResponse( coordinate.row(), coordinate.column() ) = static_cast<float>( response );
        }
        for( int stampColumn = 0; stampColumn < outputSize; ++stampColumn )
        {
            for( int stampRow = 0; stampRow < outputSize; ++stampRow )
            {
                const std::pair<double, double> detector = inverseRotate( sourceSkyRow + stampRow - 1,
                                                                          sourceSkyColumn + stampColumn - 1,
                                                                          centerRow,
                                                                          centerColumn,
                                                                          angles[central] );
                directFrames.image( static_cast<int>( central ) )( stampRow, stampColumn ) =
                    sampleImageZero( detectorResponse, detector.first, detector.second );
            }
        }
    }
    imageT expected;
    directFrames.mean( expected, directValidity, 1.0F );

    reconstructorT
        reconstructor( detectorRows, detectorColumns, centerRow, centerColumn, outputSize, localRows, localColumns );
    imageT compact;
    validityT compactValidity;
    const std::vector<int> searchRegions( grid.searchPixelCount(), 0 );
    reconstructor.reconstructCombinedTemporal( compact,
                                               compactValidity,
                                               localModels,
                                               temporalCoefficients,
                                               temporalOffsets,
                                               model,
                                               localValidity,
                                               searchIndex,
                                               searchRegions,
                                               { componentCount },
                                               selections,
                                               sourceSkyRow,
                                               sourceSkyColumn,
                                               angles,
                                               mx::improc::HCI::combine::mean,
                                               {},
                                               0,
                                               1 );
    REQUIRE( compactValidity.sum() == outputSize * outputSize );
    for( int column = 0; column < outputSize; ++column )
    {
        for( int row = 0; row < outputSize; ++row )
        {
            REQUIRE( compact( row, column ) == Approx( expected( row, column ) ).margin( 2e-6 ) );
        }
    }
}

/// Verify bounded multi-frame reconstruction follows every supported final estimator.
/** This exercises mx::improc::P4PSFReconstructor::reconstructCombined() and
 * mx::improc::P4PSFReconstructor::combineFrames() against direct frozen detector injections for unweighted mean,
 * weighted mean, median, sigma fallback, and active sigma clipping.
 * \ingroup P4PSFReconstructor_unit_tests
 */
TEST_CASE( "P4 compact response stack matches direct combination", "[P4PSFReconstructor][oracle][combine]" )
{
    constexpr int detectorRows = 43;
    constexpr int detectorColumns = 45;
    constexpr double centerRow = 21.0;
    constexpr double centerColumn = 22.0;
    constexpr int outputSize = 3;
    constexpr int localRows = 19;
    constexpr int localColumns = 20;
    constexpr double sourceSkyRow = 21.0;
    constexpr double sourceSkyColumn = 29.0;
    const std::vector<double> angles{ 0.0, 0.23, -0.31, 0.11 };

    gridT grid;
    grid.resize( detectorRows, detectorColumns, centerRow, centerColumn );
    grid.region( testRegion(), nullptr );
    const imageT psfTemplate = asymmetricTemplate( 15, 16 );
    const modelT model( psfTemplate, localRows, localColumns );
    imageT localModels;
    validityT localValidity;
    searchIndexT searchIndex;
    std::vector<coefficientT> coefficients;
    prepareCompactModels( localModels, localValidity, searchIndex, coefficients, model, grid );

    reconstructorT::cubeT directFrames;
    reconstructorT::cubeT directValidity;
    directFrames.resize( outputSize, outputSize, static_cast<int>( angles.size() ) );
    directValidity.resize( outputSize, outputSize, static_cast<int>( angles.size() ) );
    for( int image = 0; image < directValidity.planes(); ++image )
    {
        directValidity.image( image ).setOnes();
    }
    for( std::size_t image = 0; image < angles.size(); ++image )
    {
        const std::pair<double, double> sourceDetector =
            inverseRotate( sourceSkyRow, sourceSkyColumn, centerRow, centerColumn, angles[image] );
        const imageT detectorResponse =
            directDetectorResponse( psfTemplate, grid, coefficients, sourceDetector.first, sourceDetector.second );
        for( int stampColumn = 0; stampColumn < outputSize; ++stampColumn )
        {
            for( int stampRow = 0; stampRow < outputSize; ++stampRow )
            {
                const std::pair<double, double> detector = inverseRotate( sourceSkyRow + stampRow - 1,
                                                                          sourceSkyColumn + stampColumn - 1,
                                                                          centerRow,
                                                                          centerColumn,
                                                                          angles[image] );
                directFrames.image( static_cast<int>( image ) )( stampRow, stampColumn ) =
                    sampleImageZero( detectorResponse, detector.first, detector.second );
            }
        }
    }

    const mx::improc::HCI::combine method = GENERATE( mx::improc::HCI::combine::mean,
                                                      mx::improc::HCI::combine::median,
                                                      mx::improc::HCI::combine::sigmaMean );
    const bool weighted = GENERATE( false, true );
    if( method == mx::improc::HCI::combine::median && weighted )
    {
        SUCCEED( "median does not use configured weights" );
        return;
    }
    const float sigmaThreshold = method == mx::improc::HCI::combine::sigmaMean ? GENERATE( -1.0F, 1.25F ) : 0.0F;
    std::vector<float> weights = weighted ? std::vector<float>{ 0.1F, 0.2F, 0.3F, 0.4F } : std::vector<float>{};

    reconstructorT
        reconstructor( detectorRows, detectorColumns, centerRow, centerColumn, outputSize, localRows, localColumns );
    imageT compact;
    validityT compactValidity;
    reconstructor.reconstructCombined( compact,
                                       compactValidity,
                                       localModels,
                                       localValidity,
                                       searchIndex,
                                       0,
                                       sourceSkyRow,
                                       sourceSkyColumn,
                                       angles,
                                       method,
                                       weights,
                                       sigmaThreshold,
                                       1.0F );

    imageT expected;
    if( method == mx::improc::HCI::combine::median )
    {
        directFrames.median( expected, directValidity, 1.0F );
    }
    else if( method == mx::improc::HCI::combine::mean || sigmaThreshold <= 0 )
    {
        if( weights.empty() )
        {
            directFrames.mean( expected, directValidity, 1.0F );
        }
        else
        {
            directFrames.mean( expected, weights, directValidity, 1.0F );
        }
    }
    else if( weights.empty() )
    {
        directFrames.sigmaMean( expected, directValidity, sigmaThreshold, 1.0F );
    }
    else
    {
        directFrames.sigmaMean( expected, weights, directValidity, sigmaThreshold, 1.0F );
    }
    REQUIRE( compactValidity.sum() == outputSize * outputSize );
    for( int column = 0; column < outputSize; ++column )
    {
        for( int row = 0; row < outputSize; ++row )
        {
            REQUIRE( compact( row, column ) == Approx( expected( row, column ) ).margin( 1e-6 ) );
        }
    }
}

/// Verify target-specific compact columns reproduce explicit per-frame reconstruction and combination.
/** This exercises mx::improc::P4PSFReconstructor::reconstructCombinedTargeted() against an oracle assembled from
 * mx::improc::P4PSFReconstructor::reconstructFrame() and mx::improc::P4PSFReconstructor::combineFrames(), including
 * search-major/target-major model layout, distinct response columns, and target-specific invalid search pixels.
 * \ingroup P4PSFReconstructor_unit_tests
 */
TEST_CASE( "P4 targeted response matches per-frame reconstruction oracle",
           "[P4PSFReconstructor][targeted][oracle][combine]" )
{
    constexpr int detectorRows = 43;
    constexpr int detectorColumns = 45;
    constexpr double centerRow = 21.0;
    constexpr double centerColumn = 22.0;
    constexpr int outputSize = 3;
    constexpr int localRows = 19;
    constexpr int localColumns = 20;
    constexpr double sourceSkyRow = 21.0;
    constexpr double sourceSkyColumn = 29.0;
    const std::vector<double> angles{ 0.0, 0.23, -0.31, 0.11 };

    gridT grid;
    grid.resize( detectorRows, detectorColumns, centerRow, centerColumn );
    grid.region( testRegion(), nullptr );
    const imageT psfTemplate = asymmetricTemplate( 15, 16 );
    const modelT model( psfTemplate, localRows, localColumns );
    imageT sharedModels;
    validityT sharedValidity;
    searchIndexT searchIndex;
    std::vector<coefficientT> coefficients;
    prepareCompactModels( sharedModels, sharedValidity, searchIndex, coefficients, model, grid );

    const Eigen::Index targetCount = static_cast<Eigen::Index>( angles.size() );
    imageT targetedModels( sharedModels.rows(), sharedModels.cols() * targetCount );
    validityT targetedValidity( sharedValidity.rows(), targetCount );
    targetedValidity.setZero();
    for( Eigen::Index search = 0; search < sharedModels.cols(); ++search )
    {
        for( Eigen::Index target = 0; target < targetCount; ++target )
        {
            targetedModels.col( search * targetCount + target ) =
                sharedModels.col( search ) * static_cast<float>( 1.0 + 0.2 * target );
            targetedValidity( search, target ) = sharedValidity( search, 0 );
        }
    }
    const int firstInvalidSearch = searchIndex( 21, 29 );
    const int secondInvalidSearch = searchIndex( 22, 28 );
    REQUIRE( firstInvalidSearch >= 0 );
    REQUIRE( secondInvalidSearch >= 0 );
    targetedValidity( firstInvalidSearch, 1 ) = 0;
    targetedValidity( secondInvalidSearch, 2 ) = 0;

    reconstructorT
        reconstructor( detectorRows, detectorColumns, centerRow, centerColumn, outputSize, localRows, localColumns );
    reconstructorT::cubeT frames;
    reconstructorT::cubeT frameValidity;
    frames.resize( outputSize, outputSize, static_cast<int>( angles.size() ) );
    frameValidity.resize( outputSize, outputSize, static_cast<int>( angles.size() ) );
    frames.setZero();
    frameValidity.setZero();
    for( std::size_t target = 0; target < angles.size(); ++target )
    {
        imageT frame;
        validityT validity;
        reconstructor.reconstructFrame( frame,
                                        validity,
                                        targetedModels,
                                        targetedValidity,
                                        searchIndex,
                                        target,
                                        sourceSkyRow,
                                        sourceSkyColumn,
                                        angles[target] );
        frames.image( static_cast<int>( target ) ) = frame;
        frameValidity.image( static_cast<int>( target ) ) = validity.cast<float>();
    }
    REQUIRE( frameValidity.image( 1 ).sum() < frameValidity.image( 0 ).sum() );

    const mx::improc::HCI::combine method = GENERATE( mx::improc::HCI::combine::mean,
                                                      mx::improc::HCI::combine::median,
                                                      mx::improc::HCI::combine::sigmaMean );
    const float sigmaThreshold = method == mx::improc::HCI::combine::sigmaMean ? 1.25F : 0.0F;
    const std::vector<float> weights = method == mx::improc::HCI::combine::median
                                           ? std::vector<float>{}
                                           : std::vector<float>{ 0.1F, 0.2F, 0.3F, 0.4F };
    constexpr float minimumGoodFraction = 0.5F;
    imageT expected;
    validityT expectedValidity;
    reconstructorT::combineFrames( expected,
                                   expectedValidity,
                                   frames,
                                   frameValidity,
                                   method,
                                   weights,
                                   sigmaThreshold,
                                   minimumGoodFraction );

    imageT actual;
    validityT actualValidity;
    reconstructor.reconstructCombinedTargeted( actual,
                                               actualValidity,
                                               targetedModels,
                                               targetedValidity,
                                               searchIndex,
                                               sourceSkyRow,
                                               sourceSkyColumn,
                                               angles,
                                               method,
                                               weights,
                                               sigmaThreshold,
                                               minimumGoodFraction );
    REQUIRE( actualValidity.rows() == expectedValidity.rows() );
    REQUIRE( actualValidity.cols() == expectedValidity.cols() );
    for( int column = 0; column < outputSize; ++column )
    {
        for( int row = 0; row < outputSize; ++row )
        {
            REQUIRE( actualValidity( row, column ) == expectedValidity( row, column ) );
            REQUIRE( actual( row, column ) == Approx( expected( row, column ) ).margin( 1e-6 ) );
        }
    }
}

/// Verify reconstruction validity follows complete detector and local-model support.
/** This exercises mx::improc::P4PSFReconstructor::reconstructFrame() invalid propagation and contract validation.
 * \ingroup P4PSFReconstructor_unit_tests
 */
TEST_CASE( "P4 compact sky reconstruction validates support and dimensions", "[P4PSFReconstructor][validation]" )
{
    REQUIRE_THROWS( reconstructorT( 0, 5, 2, 2, 3, 9, 9 ) );
    REQUIRE_THROWS( reconstructorT( 5, 5, std::numeric_limits<double>::quiet_NaN(), 2, 3, 9, 9 ) );
    REQUIRE_THROWS( reconstructorT( 5, 5, 2, 2, 0, 9, 9 ) );

    constexpr int detectorRows = 43;
    constexpr int detectorColumns = 45;
    gridT grid;
    grid.resize( detectorRows, detectorColumns, 21.0, 22.0 );
    grid.region( testRegion(), nullptr );
    const imageT psfTemplate = asymmetricTemplate( 15, 16 );
    const modelT model( psfTemplate, 19, 20 );
    imageT localModels;
    validityT localValidity;
    searchIndexT searchIndex;
    std::vector<coefficientT> coefficients;
    prepareCompactModels( localModels, localValidity, searchIndex, coefficients, model, grid );
    reconstructorT reconstructor( detectorRows, detectorColumns, 21.0, 22.0, 3, 19, 20 );
    imageT output;
    validityT outputValidity;

    REQUIRE_THROWS( reconstructor.reconstructFrame( output,
                                                    outputValidity,
                                                    localModels.topRows( 1 ),
                                                    localValidity,
                                                    searchIndex,
                                                    0,
                                                    21,
                                                    29,
                                                    0 ) );
    REQUIRE_THROWS(
        reconstructor
            .reconstructFrame( output, outputValidity, localModels, localValidity, searchIndex, 1, 21, 29, 0 ) );
    REQUIRE_THROWS( reconstructor.reconstructFrame( output,
                                                    outputValidity,
                                                    localModels,
                                                    localValidity,
                                                    searchIndex,
                                                    0,
                                                    std::numeric_limits<double>::infinity(),
                                                    29,
                                                    0 ) );
    reconstructorT::cubeT frames;
    reconstructorT::cubeT frameValidity;
    REQUIRE_THROWS(
        reconstructor
            .reconstructFrames( frames, frameValidity, localModels, localValidity, searchIndex, 0, 21, 29, {} ) );
    frames.resize( 1, 1, 2 );
    frames.setZero();
    frameValidity.resize( 1, 1, 2 );
    for( int image = 0; image < frameValidity.planes(); ++image )
    {
        frameValidity.image( image ).setOnes();
    }
    REQUIRE_THROWS( reconstructorT::combineFrames( output,
                                                   outputValidity,
                                                   frames,
                                                   frameValidity,
                                                   mx::improc::HCI::combine::none,
                                                   {},
                                                   0,
                                                   0 ) );
    REQUIRE_THROWS( reconstructorT::combineFrames( output,
                                                   outputValidity,
                                                   frames,
                                                   frameValidity,
                                                   mx::improc::HCI::combine::mean,
                                                   { 1 },
                                                   0,
                                                   0 ) );

    reconstructor.reconstructFrame( output, outputValidity, localModels, localValidity, searchIndex, 0, 21, 29, 0 );
    REQUIRE( outputValidity.sum() == 9 );
    const int invalidSearch = searchIndex( 21, 29 );
    REQUIRE( invalidSearch >= 0 );
    localValidity( invalidSearch, 0 ) = 0;
    reconstructor.reconstructFrame( output, outputValidity, localModels, localValidity, searchIndex, 0, 21, 29, 0 );
    REQUIRE( outputValidity.sum() < 9 );

#ifdef __DOXY_ONLY__
    // clang-format off
    mx::improc::P4PSFReconstructor::reconstructFrame(
        output, outputValidity, localModels, localValidity, searchIndex, 0, 21, 29, 0 );
    // clang-format on
#endif
}

} // namespace P4PSFReconstructor_test
} // namespace unitTest
