/** \file P4PSFModel_test.cpp
 * \brief Tests compact frozen-model PSF calculations for Pixel Prediction Post-Processing.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "src/common/P4PSFModel.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace unitTest
{
namespace P4PSFModel_test
{

/** \cond P4PSFModel_test_harness */
using modelT = mx::improc::P4PSFModel;
using gridT = modelT::gridT;
using imageT = modelT::imageT;
using coefficientT = modelT::coefficientT;
using probeMatrixT = modelT::probeMatrixT;
using probeVectorT = modelT::probeVectorT;

/// Return a compact region safely contained by the test detector image.
mx::improc::P4PixelGridRegion testRegion()
{
    return mx::improc::P4PixelGridRegion( 5.0,
                                          6.0,
                                          2.0,
                                          2.0,
                                          3.0,
                                          45.0,
                                          0.5,
                                          mx::improc::P4ExclusionPolicy::sampleCenter,
                                          0.0 );
}

/// Find one test search coordinate, requiring that it occurs exactly once.
std::size_t findSearch( const gridT &grid, /**< [in] configured test grid */
                        int row,           /**< [in] desired detector row */
                        int column /**< [in] desired detector column */ )
{
    std::size_t found = grid.searchPixelCount();
    for( std::size_t index = 0; index < grid.searchPixelCount(); ++index )
    {
        const mx::improc::P4PixelCoordinate &coordinate = grid.searchPixel( index ).coordinate();
        if( coordinate.row() == row && coordinate.column() == column )
        {
            REQUIRE( found == grid.searchPixelCount() );
            found = index;
        }
    }
    REQUIRE( found < grid.searchPixelCount() );
    return found;
}

/// Construct a finite asymmetric centered template with no reflection or transpose symmetry.
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
                static_cast<float>( std::exp( -0.13 * deltaRow * deltaRow - 0.07 * deltaColumn * deltaColumn ) *
                                    ( 1.0 + 0.03 * deltaRow - 0.02 * deltaColumn + 0.004 * deltaRow * deltaColumn ) );
        }
    }
    return psf;
}

/// Independently sample a test template with production cubic interpolation and zero exterior padding.
float sampleTemplate( const imageT &psfTemplate, /**< [in] centered PSF template */
                      double row,                /**< [in] template row coordinate */
                      double column /**< [in] template column coordinate */ )
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
            if( imageRow >= 0 && imageRow < psfTemplate.rows() && imageColumn >= 0 && imageColumn < psfTemplate.cols() )
            {
                value += psfTemplate( imageRow, imageColumn ) * kernel( rowOffset, columnOffset );
            }
        }
    }
    return value;
}

/// Construct a detector image containing the shifted, detector-sampled test PSF.
imageT injectedDetectorImage( const imageT &psfTemplate, /**< [in] centered PSF template */
                              const gridT &grid,         /**< [in] detector geometry */
                              double sourceRow,          /**< [in] detector source-center row */
                              double sourceColumn /**< [in] detector source-center column */ )
{
    imageT detector( grid.rows(), grid.columns() );
    const double templateCenterRow = 0.5 * static_cast<double>( psfTemplate.rows() - 1 );
    const double templateCenterColumn = 0.5 * static_cast<double>( psfTemplate.cols() - 1 );
    for( int column = 0; column < grid.columns(); ++column )
    {
        for( int row = 0; row < grid.rows(); ++row )
        {
            detector( row, column ) =
                sampleTemplate( psfTemplate,
                                templateCenterRow + static_cast<double>( row ) - sourceRow,
                                templateCenterColumn + static_cast<double>( column ) - sourceColumn );
        }
    }
    return detector;
}

/// Calculate one local stamp by explicitly injecting and sampling a complete detector image per element.
imageT bruteForceResponse( const imageT &psfTemplate,        /**< [in] centered PSF template */
                           const gridT &grid,                /**< [in] configured detector grid */
                           std::size_t searchIndex,          /**< [in] search-pixel index */
                           const coefficientT &coefficients, /**< [in] predictor coefficients */
                           int stampSize /**< [in] output stamp size */ )
{
    imageT response( stampSize, stampSize );
    const mx::improc::P4PixelCoordinate &target = grid.searchPixel( searchIndex ).coordinate();
    const double stampCenter = 0.5 * static_cast<double>( stampSize - 1 );
    for( int stampColumn = 0; stampColumn < stampSize; ++stampColumn )
    {
        const double deltaColumn = static_cast<double>( stampColumn ) - stampCenter;
        for( int stampRow = 0; stampRow < stampSize; ++stampRow )
        {
            const double deltaRow = static_cast<double>( stampRow ) - stampCenter;
            const imageT detector = injectedDetectorImage( psfTemplate,
                                                           grid,
                                                           static_cast<double>( target.row() ) - deltaRow,
                                                           static_cast<double>( target.column() ) - deltaColumn );
            double residual = detector( target.row(), target.column() );
            for( std::size_t predictor = 0; predictor < grid.predictorCount(); ++predictor )
            {
                residual -= coefficients( predictor ) * grid.sample( detector, searchIndex, predictor );
            }
            response( stampRow, stampColumn ) = static_cast<float>( residual );
        }
    }
    return response;
}

/// Return a wider test region whose complete derotation footprints remain locally supported.
mx::improc::P4PixelGridRegion wideTestRegion()
{
    return mx::improc::P4PixelGridRegion( 3.0,
                                          12.0,
                                          2.0,
                                          2.0,
                                          3.0,
                                          45.0,
                                          0.5,
                                          mx::improc::P4ExclusionPolicy::sampleCenter,
                                          0.0 );
}

/// Map one output coordinate to the input coordinate used by production image rotation.
std::pair<double, double> inverseRotate( double row,          /**< [in] output row */
                                         double column,       /**< [in] output column */
                                         double centerRow,    /**< [in] geometric image-center row */
                                         double centerColumn, /**< [in] geometric image-center column */
                                         double angle /**< [in] counterclockwise output rotation in radians */ )
{
    const double deltaRow = row - centerRow;
    const double deltaColumn = column - centerColumn;
    return { centerRow + deltaRow * std::cos( angle ) + deltaColumn * std::sin( angle ),
             centerColumn - deltaRow * std::sin( angle ) + deltaColumn * std::cos( angle ) };
}

/// Sample one local-response stamp with the production cubic convention and complete support.
float sampleLocalResponse( const imageT &response, /**< [in] local response sampled about its geometric center */
                           double deltaRow,        /**< [in] detector target-minus-source row offset */
                           double deltaColumn /**< [in] detector target-minus-source column offset */ )
{
    const double centerRow = 0.5 * static_cast<double>( response.rows() - 1 );
    const double centerColumn = 0.5 * static_cast<double>( response.cols() - 1 );
    const double sampleRow = centerRow + deltaRow;
    const double sampleColumn = centerColumn + deltaColumn;
    const int anchorRow = static_cast<int>( std::floor( sampleRow ) );
    const int anchorColumn = static_cast<int>( std::floor( sampleColumn ) );
    const int footprintRow = anchorRow - gridT::leftBuffer;
    const int footprintColumn = anchorColumn - gridT::leftBuffer;
    REQUIRE( footprintRow >= 0 );
    REQUIRE( footprintColumn >= 0 );
    REQUIRE( footprintRow + gridT::width <= response.rows() );
    REQUIRE( footprintColumn + gridT::width <= response.cols() );

    gridT::transformT transform;
    gridT::kernelT kernel;
    transform( kernel,
               static_cast<float>( sampleRow - std::floor( sampleRow ) ),
               static_cast<float>( sampleColumn - std::floor( sampleColumn ) ) );
    return ( response.block( footprintRow, footprintColumn, gridT::width, gridT::width ) * kernel ).sum();
}

/// Calculate a direct frozen-model detector response for one source position.
imageT directDetectorResponse( const imageT &psfTemplate,                     /**< [in] centered template */
                               const gridT &grid,                             /**< [in] complete detector geometry */
                               const std::vector<coefficientT> &coefficients, /**< [in] coefficients by search pixel */
                               double sourceRow,                              /**< [in] detector source row */
                               double sourceColumn /**< [in] detector source column */ )
{
    const imageT detector = injectedDetectorImage( psfTemplate, grid, sourceRow, sourceColumn );
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

/// Require two float images to agree element by element.
void requireApprox( const imageT &actual,   /**< [in] calculated values */
                    const imageT &expected, /**< [in] oracle values */
                    double tolerance = 2e-5 /**< [in] absolute tolerance */ )
{
    REQUIRE( actual.rows() == expected.rows() );
    REQUIRE( actual.cols() == expected.cols() );
    for( Eigen::Index column = 0; column < actual.cols(); ++column )
    {
        for( Eigen::Index row = 0; row < actual.rows(); ++row )
        {
            REQUIRE( actual( row, column ) == Approx( expected( row, column ) ).margin( tolerance ) );
        }
    }
}
/** \endcond */

/// Verify the compact local operator matches explicit detector injection and P4 sampling.
/** This exercises mx::improc::P4PSFModel::calculateLocalResponse() with asymmetric odd/even template and stamp
 * geometry, an off-axis search pixel, fractional predictor coordinates, and zero exterior template padding.
 * \ingroup P4PSFModel_unit_tests
 */
TEST_CASE( "P4 local PSF response matches brute-force detector sampling", "[P4PSFModel][oracle][geometry]" )
{
    gridT grid;
    grid.resize( 41, 43, 20.0, 21.0 );
    grid.region( testRegion(), nullptr );
    const std::size_t searchIndex = findSearch( grid, 23, 25 );

    coefficientT coefficients( static_cast<Eigen::Index>( grid.predictorCount() ) );
    for( Eigen::Index predictor = 0; predictor < coefficients.rows(); ++predictor )
    {
        coefficients( predictor ) = ( predictor % 2 == 0 ? 1.0 : -1.0 ) * 0.02 / static_cast<double>( predictor + 1 );
    }

    SECTION( "odd template and stamp" )
    {
        const imageT psfTemplate = asymmetricTemplate( 11, 13 );
        const modelT model( psfTemplate, 5 );
        imageT compact;
        model.calculateLocalResponse( compact, grid, searchIndex, coefficients );
        REQUIRE( model.stampRows() == 5 );
        REQUIRE( model.stampColumns() == 5 );
        REQUIRE( model.storageBytes() > static_cast<std::size_t>( psfTemplate.size() ) * sizeof( float ) );
        requireApprox( compact, bruteForceResponse( psfTemplate, grid, searchIndex, coefficients, 5 ) );
    }

    SECTION( "even template and stamp" )
    {
        const imageT psfTemplate = asymmetricTemplate( 10, 12 );
        const modelT model( psfTemplate, 4 );
        imageT compact;
        model.calculateLocalResponse( compact, grid, searchIndex, coefficients );
        requireApprox( compact, bruteForceResponse( psfTemplate, grid, searchIndex, coefficients, 4 ) );
    }

    SECTION( "small template exercises zero padding" )
    {
        const imageT psfTemplate = asymmetricTemplate( 3, 4 );
        const modelT model( psfTemplate, 7 );
        imageT compact;
        model.calculateLocalResponse( compact, grid, searchIndex, coefficients );
        requireApprox( compact, bruteForceResponse( psfTemplate, grid, searchIndex, coefficients, 7 ) );
    }
}

/// Verify frozen probe inputs reproduce the ordinary compact residual operator.
/** This exercises mx::improc::P4PSFModel::responseInputs() and
 * mx::improc::P4PSFModel::calculateLocalResponse() with finite asymmetric coefficients, fractional predictor
 * phases, rectangular odd/even stamps, and template-edge zero padding.
 * \ingroup P4PSFModel_unit_tests
 */
TEST_CASE( "P4 PSF probe inputs reproduce local response", "[P4PSFModel][probe][oracle][geometry]" )
{
    gridT grid;
    grid.resize( 41, 43, 20.0, 21.0 );
    grid.region( testRegion(), nullptr );
    const std::size_t searchIndex = findSearch( grid, 23, 25 );

    coefficientT coefficients( static_cast<Eigen::Index>( grid.predictorCount() ) );
    for( Eigen::Index predictor = 0; predictor < coefficients.rows(); ++predictor )
    {
        coefficients( predictor ) = 0.013 * std::sin( 0.17 * static_cast<double>( predictor + 1 ) ) -
                                    0.007 * std::cos( 0.11 * static_cast<double>( 2 * predictor + 3 ) );
    }

    const int stampRows = GENERATE( 4, 5, 7 );
    const int stampColumns = GENERATE( 5, 6 );
    const int templateRows = stampRows == 7 ? 3 : 10;
    const int templateColumns = stampRows == 7 ? 4 : 13;
    const imageT psfTemplate = asymmetricTemplate( templateRows, templateColumns );
    const modelT model( psfTemplate, stampRows, stampColumns );

    probeVectorT probeTarget;
    probeMatrixT probePredictors;
    model.responseInputs( probeTarget, probePredictors, grid, searchIndex );
    REQUIRE( probeTarget.rows() == stampRows * stampColumns );
    REQUIRE( probePredictors.rows() == stampRows * stampColumns );
    REQUIRE( probePredictors.cols() == coefficients.rows() );

    const probeVectorT probeResidual =
        ( probeTarget.matrix() - probePredictors.matrix() * coefficients.matrix() ).array();
    imageT localResponse;
    model.calculateLocalResponse( localResponse, grid, searchIndex, coefficients );
    for( int stampColumn = 0; stampColumn < stampColumns; ++stampColumn )
    {
        for( int stampRow = 0; stampRow < stampRows; ++stampRow )
        {
            const Eigen::Index stampPixel = stampRow + stampRows * stampColumn;
            REQUIRE( probeResidual( stampPixel ) == Approx( localResponse( stampRow, stampColumn ) ).margin( 2e-6 ) );
        }
    }
}

/// Verify temporal PSF components preserve predictor-slot ordering and direct detector sampling.
/** This exercises mx::improc::P4PSFModel::calculateLocalResponseComponents() against explicit shifted-template
 * images for two temporal slots with direct PSF-disk predictor offsets.
 * \ingroup P4PSFModel_unit_tests
 */
TEST_CASE( "P4 temporal local PSF components match direct sampling", "[P4PSFModel][temporal][oracle]" )
{
    gridT grid;
    grid.resize( 41, 43, 20.0, 21.0 );
    grid.region( testRegion(), nullptr );
    const std::size_t searchIndex = findSearch( grid, 23, 25 );
    const mx::improc::P4PixelCoordinate &target = grid.searchPixel( searchIndex ).coordinate();
    const std::vector<mx::improc::P4PixelCoordinate> temporalOffsets{ mx::improc::P4PixelCoordinate( -1, 0 ),
                                                                      mx::improc::P4PixelCoordinate( 0, 1 ) };
    constexpr std::size_t temporalImageCount = 2;
    coefficientT coefficients(
        static_cast<Eigen::Index>( grid.predictorCount() + temporalImageCount * temporalOffsets.size() ) );
    for( Eigen::Index predictor = 0; predictor < coefficients.rows(); ++predictor )
    {
        coefficients( predictor ) = 0.001 * static_cast<double>( predictor + 1 );
    }

    constexpr int stampSize = 5;
    const imageT psfTemplate = asymmetricTemplate( 11, 13 );
    const modelT model( psfTemplate, stampSize );
    imageT components;
    model.calculateLocalResponseComponents( components,
                                            grid,
                                            searchIndex,
                                            temporalOffsets,
                                            temporalImageCount,
                                            coefficients );
    REQUIRE( components.rows() == stampSize * stampSize );
    REQUIRE( components.cols() == 1 + static_cast<Eigen::Index>( temporalImageCount ) );

    const double stampCenter = 0.5 * static_cast<double>( stampSize - 1 );
    for( int stampColumn = 0; stampColumn < stampSize; ++stampColumn )
    {
        for( int stampRow = 0; stampRow < stampSize; ++stampRow )
        {
            const double sourceRow = static_cast<double>( target.row() ) - ( stampRow - stampCenter );
            const double sourceColumn = static_cast<double>( target.column() ) - ( stampColumn - stampCenter );
            const imageT detector = injectedDetectorImage( psfTemplate, grid, sourceRow, sourceColumn );
            double baseResponse = detector( target.row(), target.column() );
            for( std::size_t predictor = 0; predictor < grid.predictorCount(); ++predictor )
            {
                baseResponse -= coefficients( static_cast<Eigen::Index>( predictor ) ) *
                                grid.sample( detector, searchIndex, predictor );
            }
            const Eigen::Index stampPixel = stampRow + stampSize * stampColumn;
            REQUIRE( components( stampPixel, 0 ) == Approx( baseResponse ).margin( 2e-5 ) );
            for( std::size_t temporalImage = 0; temporalImage < temporalImageCount; ++temporalImage )
            {
                double temporalResponse{ 0 };
                for( std::size_t predictor = 0; predictor < temporalOffsets.size(); ++predictor )
                {
                    const mx::improc::P4PixelCoordinate &offset = temporalOffsets[predictor];
                    const Eigen::Index coefficientIndex = static_cast<Eigen::Index>(
                        grid.predictorCount() + temporalImage * temporalOffsets.size() + predictor );
                    temporalResponse -= coefficients( coefficientIndex ) *
                                        detector( target.row() + offset.row(), target.column() + offset.column() );
                }
                REQUIRE( components( stampPixel, static_cast<Eigen::Index>( temporalImage + 1 ) ) ==
                         Approx( temporalResponse ).margin( 2e-5 ) );
            }
        }
    }
}

/// Measure compact-stamp reconstruction against direct frozen injection after nonzero derotation.
/** This exercises mx::improc::P4PSFModel::calculateLocalResponse() as the compact input to the exact
 * mx::improc::imageRotate() detector-to-sky sampling convention.
 * \ingroup P4PSFModel_unit_tests
 */
TEST_CASE( "P4 compact local stamps reconstruct a rotated frozen response", "[P4PSFModel][oracle][derotation]" )
{
    constexpr int detectorRows = 41;
    constexpr int detectorColumns = 43;
    constexpr int localStampRows = 17;
    constexpr int localStampColumns = 18;
    constexpr double angle = 0.23;
    gridT grid;
    grid.resize( detectorRows, detectorColumns, 20.0, 21.0 );
    grid.region( wideTestRegion(), nullptr );
    const imageT psfTemplate = asymmetricTemplate( 15, 16 );
    const modelT model( psfTemplate, localStampRows, localStampColumns );

    std::vector<coefficientT> coefficients;
    std::vector<imageT> localResponses;
    coefficients.reserve( grid.searchPixelCount() );
    localResponses.resize( grid.searchPixelCount() );
    Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic> searchIndex =
        Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic>::Constant( detectorRows, detectorColumns, -1 );
    for( std::size_t search = 0; search < grid.searchPixelCount(); ++search )
    {
        coefficientT coefficient( static_cast<Eigen::Index>( grid.predictorCount() ) );
        for( Eigen::Index predictor = 0; predictor < coefficient.rows(); ++predictor )
        {
            coefficient( predictor ) = 0.003 * std::sin( 0.07 * static_cast<double>( predictor + 3 * search ) );
        }
        coefficients.push_back( coefficient );
        if( grid.searchPixel( search ).valid() )
        {
            model.calculateLocalResponse( localResponses[search], grid, search, coefficients.back() );
            const mx::improc::P4PixelCoordinate &coordinate = grid.searchPixel( search ).coordinate();
            searchIndex( coordinate.row(), coordinate.column() ) = static_cast<int>( search );
        }
    }

    const int sourceSkyRow = 20;
    const int sourceSkyColumn = 28;
    const auto sourceDetector = inverseRotate( sourceSkyRow, sourceSkyColumn, 20.0, 21.0, angle );
    const imageT detectorResponse =
        directDetectorResponse( psfTemplate, grid, coefficients, sourceDetector.first, sourceDetector.second );
    imageT directSkyResponse;
    mx::improc::imageRotate( directSkyResponse, detectorResponse, angle, mx::improc::cubicConvolTransform<float>() );

    gridT::transformT transform;
    double maximumError{ 0 };
    for( int skyColumn = sourceSkyColumn - 1; skyColumn <= sourceSkyColumn + 1; ++skyColumn )
    {
        for( int skyRow = sourceSkyRow - 1; skyRow <= sourceSkyRow + 1; ++skyRow )
        {
            const auto detectorCoordinate = inverseRotate( skyRow, skyColumn, 20.0, 21.0, angle );
            const int anchorRow = static_cast<int>( std::floor( detectorCoordinate.first ) );
            const int anchorColumn = static_cast<int>( std::floor( detectorCoordinate.second ) );
            const int footprintRow = anchorRow - gridT::leftBuffer;
            const int footprintColumn = anchorColumn - gridT::leftBuffer;
            gridT::kernelT kernel;
            transform( kernel,
                       static_cast<float>( detectorCoordinate.first - std::floor( detectorCoordinate.first ) ),
                       static_cast<float>( detectorCoordinate.second - std::floor( detectorCoordinate.second ) ) );

            float compactValue{ 0 };
            for( int columnOffset = 0; columnOffset < gridT::width; ++columnOffset )
            {
                for( int rowOffset = 0; rowOffset < gridT::width; ++rowOffset )
                {
                    const int detectorRow = footprintRow + rowOffset;
                    const int detectorColumn = footprintColumn + columnOffset;
                    const int search = searchIndex( detectorRow, detectorColumn );
                    REQUIRE( search >= 0 );
                    compactValue +=
                        sampleLocalResponse( localResponses[static_cast<std::size_t>( search )],
                                             static_cast<double>( detectorRow ) - sourceDetector.first,
                                             static_cast<double>( detectorColumn ) - sourceDetector.second ) *
                        kernel( rowOffset, columnOffset );
                }
            }
            maximumError =
                std::max( maximumError,
                          std::abs( static_cast<double>( compactValue - directSkyResponse( skyRow, skyColumn ) ) ) );
        }
    }
    INFO( "maximum compact reconstruction error: " << maximumError );
    REQUIRE( maximumError <= 1e-6 );
}

/// Verify the local PSF seam rejects incomplete geometry and invalid numerical inputs.
/** This exercises validation in mx::improc::P4PSFModel::calculateLocalResponse().
 * \ingroup P4PSFModel_unit_tests
 */
TEST_CASE( "P4 local PSF response validates its complete contract", "[P4PSFModel][validation]" )
{
    gridT grid;
    grid.resize( 41, 43, 20.0, 21.0 );
    grid.region( testRegion(), nullptr );
    const std::size_t searchIndex = findSearch( grid, 23, 25 );
    coefficientT coefficients( static_cast<Eigen::Index>( grid.predictorCount() ) );
    coefficients.setZero();
    const imageT validTemplate = asymmetricTemplate( 9, 9 );
    const modelT model( validTemplate, 3 );
    imageT output;

    imageT emptyTemplate;
    REQUIRE_THROWS_AS( modelT( emptyTemplate, 3 ), std::invalid_argument );

    imageT nonfiniteTemplate = validTemplate;
    nonfiniteTemplate( 0, 0 ) = std::numeric_limits<float>::quiet_NaN();
    REQUIRE_THROWS_AS( modelT( nonfiniteTemplate, 3 ), std::invalid_argument );
    REQUIRE_THROWS_AS( modelT( validTemplate, 0 ), std::invalid_argument );
    REQUIRE_THROWS_AS( modelT( validTemplate, 3, 0 ), std::invalid_argument );

    gridT incompleteGrid;
    incompleteGrid.resize( 41, 43, 20.0, 21.0 );
    REQUIRE_THROWS_AS( model.calculateLocalResponse( output, incompleteGrid, 0, coefficients ), std::invalid_argument );
    REQUIRE_THROWS_AS( model.calculateLocalResponse( output, grid, grid.searchPixelCount(), coefficients ),
                       std::out_of_range );

    coefficientT wrongCount( coefficients.rows() - 1 );
    wrongCount.setZero();
    REQUIRE_THROWS_AS( model.calculateLocalResponse( output, grid, searchIndex, wrongCount ), std::invalid_argument );

    coefficients( 0 ) = std::numeric_limits<double>::infinity();
    REQUIRE_THROWS_AS( model.calculateLocalResponse( output, grid, searchIndex, coefficients ), std::invalid_argument );
    coefficients.setZero();

    imageT commonMask = imageT::Ones( grid.rows(), grid.columns() );
    const mx::improc::P4PixelCoordinate &target = grid.searchPixel( searchIndex ).coordinate();
    commonMask( target.row(), target.column() ) = 0;
    gridT maskedGrid;
    maskedGrid.resize( 41, 43, 20.0, 21.0 );
    maskedGrid.region( testRegion(), &commonMask );
    const std::size_t maskedIndex = findSearch( maskedGrid, target.row(), target.column() );
    REQUIRE_THROWS_AS( model.calculateLocalResponse( output, maskedGrid, maskedIndex, coefficients ),
                       std::invalid_argument );
}

} // namespace P4PSFModel_test
} // namespace unitTest
