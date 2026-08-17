/** \file P4RotatedGrid_test.cpp
 * \brief Tests direct sky-frame sampling geometry for rotated P4 regression.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "src/common/P4RotatedGrid.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include <mx/improc/imageTransforms.hpp>

namespace mx
{
namespace improc
{

/** \cond P4RotatedGrid_test_harness */
class P4RotatedGridTestAccess
{
  public:
    /// Invoke the production fixed-coordinate overflow check.
    static std::size_t checkedCoordinateCount( std::size_t searchCount, /**< [in] target count */
                                               std::size_t predictorCount /**< [in] predictor count */ )
    {
        return P4RotatedGrid::checkedCoordinateCount( searchCount, predictorCount );
    }

    /// Invoke the production raw-coordinate and footprint range checks.
    static void mapFootprint( int rows,       /**< [in] raw-image rows */
                              int columns,    /**< [in] raw-image columns */
                              double xCenter, /**< [in] center row */
                              double yCenter, /**< [in] center column */
                              double cosine,  /**< [in] frame cosine */
                              double sine,    /**< [in] frame sine */
                              double skyRow,  /**< [in] fixed sky row */
                              double skyColumn /**< [in] fixed sky column */ )
    {
        static_cast<void>( P4RotatedGrid::mapFootprint( rows,
                                                        columns,
                                                        xCenter,
                                                        yCenter,
                                                        P4RotatedGrid::FrameRotation{ cosine, sine },
                                                        P4RotatedGrid::SkyCoordinate{ skyRow, skyColumn } ) );
    }

    /// Invoke the production direct sampler with an intentionally incomplete record.
    static float sampleIncomplete( const P4RotatedGrid::imageT &image /**< [in] raw test image */ )
    {
        const auto record = P4RotatedGrid::makeInterpolation( image.rows(),
                                                              image.cols(),
                                                              2.0,
                                                              2.0,
                                                              P4RotatedGrid::FrameRotation{ 1.0, 0.0 },
                                                              P4RotatedGrid::SkyCoordinate{ 0.0, 0.0 } );
        return P4RotatedGrid::sampleRecord( image, record );
    }
};
/** \endcond */

} // namespace improc
} // namespace mx

namespace unitTest
{
namespace P4RotatedGrid_test
{

/** \cond P4RotatedGrid_test_harness */
using skyGridT = mx::improc::P4PixelGridf;
using rotatedGridT = mx::improc::P4RotatedGrid;
using imageT = rotatedGridT::imageT;
using policyT = mx::improc::P4ExclusionPolicy;
using reasonT = mx::improc::P4RotatedInvalidReason;
using regionT = mx::improc::P4PixelGridRegion;

/// Return a compact region with enough candidates for exclusion-policy tests.
regionT testRegion( policyT policy = policyT::sampleCenter, /**< [in] central-exclusion policy */
                    double psfRadius = 0.5,                 /**< [in] physical exclusion radius */
                    double radiusBuffer = 0.0 /**< [in] explicit exclusion comparison buffer */ )
{
    return regionT( 5.0, 6.0, 2.0, 2.0, 4.0, 60.0, psfRadius, policy, radiusBuffer );
}

/// Construct candidate-only fixed sky geometry.
skyGridT candidateGrid( int rows = 31,         /**< [in] image rows */
                        int columns = 33,      /**< [in] image columns */
                        double xCenter = 15.0, /**< [in] center row */
                        double yCenter = 16.0, /**< [in] center column */
                        const regionT &region = testRegion() /**< [in] candidate-region configuration */ )
{
    skyGridT grid;
    grid.resize( rows, columns, xCenter, yCenter );
    grid.candidateRegion( region );
    return grid;
}

/// Find one search coordinate, requiring exactly one match.
std::size_t findSearch( const skyGridT &grid, /**< [in] configured candidate grid */
                        int row,              /**< [in] desired fixed sky row */
                        int column /**< [in] desired fixed sky column */ )
{
    std::size_t result = grid.searchPixelCount();
    for( std::size_t search = 0; search < grid.searchPixelCount(); ++search )
    {
        const auto &coordinate = grid.searchPixel( search ).coordinate();
        if( coordinate.row() == row && coordinate.column() == column )
        {
            REQUIRE( result == grid.searchPixelCount() );
            result = search;
        }
    }
    REQUIRE( result < grid.searchPixelCount() );
    return result;
}

/// Find one all-frame-valid rotated search pixel.
std::size_t findValidSearch( const rotatedGridT &grid /**< [in] configured rotated grid */ )
{
    for( std::size_t search = 0; search < grid.searchPixelCount(); ++search )
    {
        if( grid.searchPixel( search ).valid() )
        {
            return search;
        }
    }
    FAIL( "expected at least one valid rotated P4 search pixel" );
    return grid.searchPixelCount();
}

/// Map one fixed sky coordinate to its raw source using the documented inverse derotation.
std::pair<double, double> rawCoordinate( double skyRow,    /**< [in] fixed sky row */
                                         double skyColumn, /**< [in] fixed sky column */
                                         double xCenter,   /**< [in] center row */
                                         double yCenter,   /**< [in] center column */
                                         double angle /**< [in] counter-clockwise derotation in radians */ )
{
    const double cosine = std::cos( angle );
    const double sine = std::sin( angle );
    const double deltaRow = skyRow - xCenter;
    const double deltaColumn = skyColumn - yCenter;
    return { deltaRow * cosine + deltaColumn * sine + xCenter, -deltaRow * sine + deltaColumn * cosine + yCenter };
}

/// Return whether one raw nominal predictor footprint overlaps the target-centered exclusion disk.
bool rawFootprintOverlaps( const std::pair<double, double> &predictor, /**< [in] raw predictor coordinate */
                           const std::pair<double, double> &target,    /**< [in] raw target coordinate */
                           double radius /**< [in] effective exclusion radius */ )
{
    const int firstRow = static_cast<int>( std::floor( predictor.first ) ) - rotatedGridT::leftBuffer;
    const int firstColumn = static_cast<int>( std::floor( predictor.second ) ) - rotatedGridT::leftBuffer;
    for( int rowOffset = 0; rowOffset < rotatedGridT::width; ++rowOffset )
    {
        for( int columnOffset = 0; columnOffset < rotatedGridT::width; ++columnOffset )
        {
            const double deltaRow = static_cast<double>( firstRow + rowOffset ) - target.first;
            const double deltaColumn = static_cast<double>( firstColumn + columnOffset ) - target.second;
            if( deltaRow * deltaRow + deltaColumn * deltaColumn <= radius * radius )
            {
                return true;
            }
        }
    }
    return false;
}

/// Independently calculate retained pre-exclusion candidate indices for the direct policy.
std::vector<std::size_t> expectedCandidates( const skyGridT &skyGrid, /**< [in] candidate-only sky grid */
                                             const std::vector<double> &angles /**< [in] derotation angles */ )
{
    const regionT &configuration = skyGrid.regionConfiguration();
    const double effectiveRadius = skyGrid.effectiveExclusionRadius();
    std::vector<std::size_t> retained;
    for( std::size_t candidate = 0; candidate < skyGrid.candidatePredictorCount(); ++candidate )
    {
        const auto &offset = skyGrid.candidatePredictorOffset( candidate );
        bool excluded = false;
        if( configuration.exclusionPolicy == policyT::sampleCenter )
        {
            excluded = std::hypot( static_cast<double>( offset.row() ), static_cast<double>( offset.column() ) ) <=
                       effectiveRadius;
        }
        else
        {
            for( std::size_t search = 0; search < skyGrid.searchPixelCount() && !excluded; ++search )
            {
                const auto &targetSky = skyGrid.searchPixel( search ).coordinate();
                const auto predictorSky = skyGrid.candidateCoordinate( search, candidate );
                for( double angle : angles )
                {
                    const auto target = rawCoordinate( targetSky.row(),
                                                       targetSky.column(),
                                                       skyGrid.xCenter(),
                                                       skyGrid.yCenter(),
                                                       angle );
                    const auto predictor = rawCoordinate( predictorSky.first,
                                                          predictorSky.second,
                                                          skyGrid.xCenter(),
                                                          skyGrid.yCenter(),
                                                          angle );
                    if( rawFootprintOverlaps( predictor, target, effectiveRadius ) )
                    {
                        excluded = true;
                        break;
                    }
                }
            }
        }
        if( !excluded )
        {
            retained.push_back( candidate );
        }
    }
    return retained;
}

/// Fill an image with an affine field that cubic convolution reproduces exactly away from edges.
void fillAffine( imageT &image /**< [out] image to fill */ )
{
    for( int column = 0; column < image.cols(); ++column )
    {
        for( int row = 0; row < image.rows(); ++row )
        {
            image( row, column ) = 2.0F * static_cast<float>( row ) - 3.0F * static_cast<float>( column ) + 7.0F;
        }
    }
}

/// Evaluate the smooth high-spatial-frequency field used to expose double interpolation.
float oscillatoryField( double row, /**< [in] continuous row coordinate */
                        double column /**< [in] continuous column coordinate */ )
{
    return static_cast<float>( std::sin( 0.73 * row + 0.19 * column ) + 0.6 * std::cos( 0.31 * row - 0.83 * column ) );
}

/// Sample one complete cubic footprint in an already materialized image.
float cubicSample( const imageT &image, /**< [in] image to sample */
                   double row,          /**< [in] continuous row coordinate */
                   double column /**< [in] continuous column coordinate */ )
{
    skyGridT::kernelT kernel;
    skyGridT::transformT transform;
    transform( kernel,
               static_cast<float>( row - std::floor( row ) ),
               static_cast<float>( column - std::floor( column ) ) );
    const int firstRow = static_cast<int>( std::floor( row ) ) - skyGridT::leftBuffer;
    const int firstColumn = static_cast<int>( std::floor( column ) ) - skyGridT::leftBuffer;
    REQUIRE( firstRow >= 0 );
    REQUIRE( firstColumn >= 0 );
    REQUIRE( firstRow + skyGridT::width <= image.rows() );
    REQUIRE( firstColumn + skyGridT::width <= image.cols() );

    float value{ 0 };
    for( int rowOffset = 0; rowOffset < skyGridT::width; ++rowOffset )
    {
        for( int columnOffset = 0; columnOffset < skyGridT::width; ++columnOffset )
        {
            value += image( firstRow + rowOffset, firstColumn + columnOffset ) * kernel( rowOffset, columnOffset );
        }
    }
    return value;
}
/** \endcond */

/** \brief Verifies candidate-only P4PixelGrid construction exposes raw wedge geometry without complete-region state.
 *
 * Exercises `P4PixelGrid::candidateRegion()`, candidate accessors, and normal `P4PixelGrid::region()` parity.
 * \ingroup P4RotatedGrid_unit_tests
 */
TEST_CASE( "P4PixelGrid exposes pre-exclusion candidate geometry", "[P4RotatedGrid][P4PixelGrid][geometry]" )
{
#ifdef __DOXY_ONLY__
    // clang-format off
    mx::improc::P4PixelGridf doxygenGrid;
    mx::improc::P4PixelGridRegion doxygenRegion( 1, 2, 1, 1, 1, 30, 0.5,
                                                 mx::improc::P4ExclusionPolicy::sampleCenter, 0 );
    doxygenGrid.candidateRegion( doxygenRegion );
    doxygenGrid.candidateRegionConfigured();
    doxygenGrid.candidatePredictorCount();
    doxygenGrid.candidatePredictorOffset( 0 );
    doxygenGrid.candidateCoordinate( 0, 0 );
    // clang-format on
#endif

    const regionT configuration = testRegion( policyT::sampleCenter, 1.5 );
    skyGridT candidates = candidateGrid( 31, 33, 15.0, 16.0, configuration );
    REQUIRE( candidates.candidateRegionConfigured() );
    REQUIRE_FALSE( candidates.regionConfigured() );
    REQUIRE( candidates.searchPixelCount() > 0 );
    REQUIRE( candidates.candidatePredictorCount() > 0 );
    REQUIRE( candidates.predictorCount() == 0 );
    REQUIRE_THROWS_AS( candidates.candidatePredictorOffset( candidates.candidatePredictorCount() ), std::out_of_range );
    REQUIRE_THROWS_AS( candidates.candidateCoordinate( candidates.searchPixelCount(), 0 ), std::out_of_range );

    const std::size_t originalSearchCount = candidates.searchPixelCount();
    const std::size_t originalCandidateCount = candidates.candidatePredictorCount();
    const regionT invalidConfiguration( 5.0, 5.0, 2.0, 2.0, 4.0, 60.0, 1.5, policyT::sampleCenter, 0.0 );
    REQUIRE_THROWS_AS( candidates.candidateRegion( invalidConfiguration ), std::invalid_argument );
    REQUIRE( candidates.searchPixelCount() == originalSearchCount );
    REQUIRE( candidates.candidatePredictorCount() == originalCandidateCount );
    REQUIRE( candidates.regionConfiguration().searchInnerRadius == configuration.searchInnerRadius );

    skyGridT complete;
    complete.resize( 31, 33, 15.0, 16.0 );
    complete.region( configuration, nullptr );
    REQUIRE( complete.regionConfigured() );
    REQUIRE( complete.candidateRegionConfigured() );
    REQUIRE( complete.searchPixelCount() == candidates.searchPixelCount() );
    REQUIRE( complete.candidatePredictorCount() == candidates.candidatePredictorCount() );
    REQUIRE( complete.predictorCount() < complete.candidatePredictorCount() );
    for( std::size_t candidate = 0; candidate < candidates.candidatePredictorCount(); ++candidate )
    {
        REQUIRE( complete.candidatePredictorOffset( candidate ).row() ==
                 candidates.candidatePredictorOffset( candidate ).row() );
        REQUIRE( complete.candidatePredictorOffset( candidate ).column() ==
                 candidates.candidatePredictorOffset( candidate ).column() );
    }

    const regionT edgeConfiguration( 3.5, 5.5, 1.0, 1.0, 2.0, 45.0, 0.25, policyT::sampleCenter, 0.0 );
    skyGridT edgeCandidates;
    edgeCandidates.resize( 11, 11 );
    REQUIRE_NOTHROW( edgeCandidates.candidateRegion( edgeConfiguration ) );
    REQUIRE( edgeCandidates.candidateRegionConfigured() );
    skyGridT edgeComplete;
    edgeComplete.resize( 11, 11 );
    REQUIRE_THROWS_AS( edgeComplete.region( edgeConfiguration, nullptr ), std::invalid_argument );
    REQUIRE_FALSE( edgeComplete.regionConfigured() );

    const regionT emptyCandidateConfiguration( 0.0, 0.2, 0.01, 0.01, 0.01, 1.0, 0.01, policyT::sampleCenter, 0.0 );
    skyGridT emptyCandidates;
    emptyCandidates.resize( 11, 11 );
    REQUIRE_THROWS_AS( emptyCandidates.candidateRegion( emptyCandidateConfiguration ), std::invalid_argument );
    REQUIRE_FALSE( emptyCandidates.candidateRegionConfigured() );

    candidates.resize( 31, 33 );
    REQUIRE_FALSE( candidates.candidateRegionConfigured() );
    REQUIRE( candidates.candidatePredictorCount() == 0 );
}

/** \brief Verifies direct inverse-rotation mapping for identity, cardinal, and fractional angles and explicit centers.
 *
 * Exercises `P4RotatedGrid::configure()`, target interpolation, and predictor interpolation geometry.
 * \ingroup P4RotatedGrid_unit_tests
 */
TEST_CASE( "P4RotatedGrid maps fixed sky samples directly into raw frames", "[P4RotatedGrid][geometry][rotation]" )
{
#ifdef __DOXY_ONLY__
    mx::improc::P4RotatedGrid doxygenGrid;
    mx::improc::P4PixelGridf doxygenSkyGrid;
    std::vector<double> doxygenAngles;
    doxygenGrid.configure( doxygenSkyGrid, doxygenAngles, nullptr );
    doxygenGrid.configured();
    doxygenGrid.rows();
    doxygenGrid.columns();
    doxygenGrid.xCenter();
    doxygenGrid.yCenter();
    doxygenGrid.frameCount();
    doxygenGrid.searchPixelCount();
    doxygenGrid.predictorCount();
    doxygenGrid.searchPixel( 0 );
    doxygenGrid.predictorOffset( 0 );
    doxygenGrid.candidatePredictorIndex( 0 );
    doxygenGrid.targetInterpolation( 0, 0 );
    doxygenGrid.predictorInterpolation( 0, 0, 0 );
    mx::improc::P4RotatedInterpolationRecord doxygenInterpolation;
    doxygenInterpolation.row();
    doxygenInterpolation.column();
    doxygenInterpolation.footprintRow();
    doxygenInterpolation.footprintColumn();
    doxygenInterpolation.kernel();
    doxygenInterpolation.footprintInside();
    mx::improc::P4RotatedSearchRecord doxygenSearchRecord( mx::improc::P4PixelCoordinate( 0, 0 ),
                                                           mx::improc::P4RotatedInvalidReason::none );
    doxygenSearchRecord.coordinate();
    doxygenSearchRecord.invalidReason();
    doxygenSearchRecord.valid();
    doxygenSearchRecord.hasReason( mx::improc::P4RotatedInvalidReason::targetFootprintOutside );
#endif

    const double xCenter = 15.25;
    const double yCenter = 16.75;
    skyGridT skyGrid = candidateGrid( 32, 34, xCenter, yCenter );
    const std::vector<double> angles{ 0.0, 0.5 * std::numbers::pi_v<double>, 0.293 };
    rotatedGridT grid;
    grid.configure( skyGrid, angles, nullptr );

    REQUIRE( grid.configured() );
    REQUIRE( grid.rows() == 32 );
    REQUIRE( grid.columns() == 34 );
    REQUIRE( grid.xCenter() == xCenter );
    REQUIRE( grid.yCenter() == yCenter );
    REQUIRE( grid.frameCount() == angles.size() );
    REQUIRE( grid.searchPixelCount() == skyGrid.searchPixelCount() );
    REQUIRE( grid.predictorCount() > 0 );

    const std::size_t search = findValidSearch( grid );
    const auto &target = grid.searchPixel( search ).coordinate();
    for( std::size_t frame = 0; frame < angles.size(); ++frame )
    {
        const auto expected = rawCoordinate( target.row(), target.column(), xCenter, yCenter, angles[frame] );
        const auto actual = grid.targetInterpolation( frame, search );
        REQUIRE( actual.row() == Approx( expected.first ).margin( 2e-12 ) );
        REQUIRE( actual.column() == Approx( expected.second ).margin( 2e-12 ) );
    }

    REQUIRE( grid.targetInterpolation( 0, search ).row() == Approx( target.row() ) );
    REQUIRE( grid.targetInterpolation( 0, search ).column() == Approx( target.column() ) );
    const auto cardinal = grid.targetInterpolation( 1, search );
    REQUIRE( cardinal.row() == Approx( xCenter + static_cast<double>( target.column() ) - yCenter ).margin( 2e-12 ) );
    REQUIRE( cardinal.column() == Approx( yCenter - static_cast<double>( target.row() ) + xCenter ).margin( 2e-12 ) );

    const std::size_t predictor = 0;
    const std::size_t candidate = grid.candidatePredictorIndex( predictor );
    const auto fixedSky = skyGrid.candidateCoordinate( search, candidate );
    const auto expectedPredictor = rawCoordinate( fixedSky.first, fixedSky.second, xCenter, yCenter, angles[2] );
    const auto actualPredictor = grid.predictorInterpolation( 2, search, predictor );
    REQUIRE( actualPredictor.row() == Approx( expectedPredictor.first ).margin( 2e-12 ) );
    REQUIRE( actualPredictor.column() == Approx( expectedPredictor.second ).margin( 2e-12 ) );
    REQUIRE( grid.predictorOffset( predictor ).row() == skyGrid.candidatePredictorOffset( candidate ).row() );
    REQUIRE( grid.predictorOffset( predictor ).column() == skyGrid.candidatePredictorOffset( candidate ).column() );

    REQUIRE_THROWS_AS( grid.targetInterpolation( grid.frameCount(), search ), std::out_of_range );
    REQUIRE_THROWS_AS( grid.targetInterpolation( 0, grid.searchPixelCount() ), std::out_of_range );
    REQUIRE_THROWS_AS( grid.predictorInterpolation( 0, search, grid.predictorCount() ), std::out_of_range );
    REQUIRE_THROWS_AS( grid.predictorOffset( grid.predictorCount() ), std::out_of_range );
    REQUIRE_THROWS_AS( grid.candidatePredictorIndex( grid.predictorCount() ), std::out_of_range );
}

/** \brief Verifies one direct cubic sample avoids the measurable extra error of materialized two-stage sampling.
 *
 * Exercises `P4RotatedGrid::sampleTarget()`, `samplePredictor()`, and `sampleFrame()` against analytic truth, the
 * interior integer-target `imageRotate` oracle, and fractional predictors sampled again from a materialized rotated
 * image.
 * \ingroup P4RotatedGrid_unit_tests
 */
TEST_CASE( "P4RotatedGrid applies one direct cubic interpolation", "[P4RotatedGrid][sample][analytic]" )
{
#ifdef __DOXY_ONLY__
    mx::improc::P4RotatedGrid doxygenGrid;
    mx::improc::P4RotatedGrid::imageT doxygenImage;
    float doxygenTarget;
    std::vector<float> doxygenPredictors;
    doxygenGrid.sampleTarget( doxygenImage, 0, 0 );
    doxygenGrid.samplePredictor( doxygenImage, 0, 0, 0 );
    doxygenGrid.sampleFrame( doxygenImage, 0, 0, doxygenTarget, doxygenPredictors );
#endif

    skyGridT skyGrid = candidateGrid();
    constexpr double angle = 0.293;
    rotatedGridT grid;
    grid.configure( skyGrid, { angle }, nullptr );
    const std::size_t search = findValidSearch( grid );

    imageT image( grid.rows(), grid.columns() );
    fillAffine( image );
    const auto targetRecord = grid.targetInterpolation( 0, search );
    const float expectedTarget = static_cast<float>( 2.0 * targetRecord.row() - 3.0 * targetRecord.column() + 7.0 );
    REQUIRE( grid.sampleTarget( image, 0, search ) == Approx( expectedTarget ).margin( 2e-4 ) );

    for( std::size_t predictor = 0; predictor < grid.predictorCount(); ++predictor )
    {
        const auto record = grid.predictorInterpolation( 0, search, predictor );
        const float expected = static_cast<float>( 2.0 * record.row() - 3.0 * record.column() + 7.0 );
        REQUIRE( grid.samplePredictor( image, 0, search, predictor ) == Approx( expected ).margin( 3e-4 ) );
    }

    float target = -999;
    std::vector<float> predictors( grid.predictorCount(), -999 );
    grid.sampleFrame( image, 0, search, target, predictors );
    REQUIRE( target == Approx( expectedTarget ).margin( 2e-4 ) );
    for( std::size_t predictor = 0; predictor < grid.predictorCount(); ++predictor )
    {
        REQUIRE( predictors[predictor] == Approx( grid.samplePredictor( image, 0, search, predictor ) ) );
    }

    imageT rotated;
    mx::improc::imageRotate( rotated, image, angle, mx::improc::cubicConvolTransform<float>() );
    const auto &targetCoordinate = grid.searchPixel( search ).coordinate();
    REQUIRE( target == Approx( rotated( targetCoordinate.row(), targetCoordinate.column() ) ).margin( 3e-4 ) );

    imageT oscillatory( grid.rows(), grid.columns() );
    for( int column = 0; column < oscillatory.cols(); ++column )
    {
        for( int row = 0; row < oscillatory.rows(); ++row )
        {
            oscillatory( row, column ) = oscillatoryField( row, column );
        }
    }
    mx::improc::imageRotate( rotated, oscillatory, angle, mx::improc::cubicConvolTransform<float>() );

    std::size_t fractionalPredictorCount{ 0 };
    double directSquaredError{ 0 };
    double materializedSquaredError{ 0 };
    double maximumMethodDifference{ 0 };
    for( std::size_t predictor = 0; predictor < grid.predictorCount(); ++predictor )
    {
        const auto sky = skyGrid.candidateCoordinate( search, grid.candidatePredictorIndex( predictor ) );
        const int firstRow = static_cast<int>( std::floor( sky.first ) ) - skyGridT::leftBuffer;
        const int firstColumn = static_cast<int>( std::floor( sky.second ) ) - skyGridT::leftBuffer;
        const bool fractional = std::abs( sky.first - std::round( sky.first ) ) > 0.1 ||
                                std::abs( sky.second - std::round( sky.second ) ) > 0.1;
        const bool complete = firstRow >= 0 && firstColumn >= 0 && firstRow + skyGridT::width <= rotated.rows() &&
                              firstColumn + skyGridT::width <= rotated.cols();
        if( fractional && complete )
        {
            const auto directRecord = grid.predictorInterpolation( 0, search, predictor );
            const float analytic = oscillatoryField( directRecord.row(), directRecord.column() );
            const float direct = grid.samplePredictor( oscillatory, 0, search, predictor );
            const float materialized = cubicSample( rotated, sky.first, sky.second );
            directSquaredError += static_cast<double>( direct - analytic ) * ( direct - analytic );
            materializedSquaredError += static_cast<double>( materialized - analytic ) * ( materialized - analytic );
            maximumMethodDifference =
                std::max( maximumMethodDifference, static_cast<double>( std::abs( direct - materialized ) ) );
            ++fractionalPredictorCount;
        }
    }
    REQUIRE( fractionalPredictorCount > 1 );
    REQUIRE( maximumMethodDifference > 1e-3 );
    REQUIRE( directSquaredError < materializedSquaredError );
}

/** \brief Verifies exact annulus-wide direct exclusion for sample-center and raw kernel-support policies.
 *
 * Exercises `P4RotatedGrid::configure()` candidate filtering and fixed retained `K` provenance.
 * \ingroup P4RotatedGrid_unit_tests
 */
TEST_CASE( "P4RotatedGrid applies exclusion once in direct raw geometry", "[P4RotatedGrid][kernelSupport]" )
{
    const std::vector<double> angles{ 0.0, 0.25 * std::numbers::pi_v<double>, 0.67 };
    skyGridT centerSky = candidateGrid( 31, 33, 15.0, 16.0, testRegion( policyT::sampleCenter, 1.25 ) );
    rotatedGridT centerGrid;
    centerGrid.configure( centerSky, angles, nullptr );
    const auto expectedCenter = expectedCandidates( centerSky, angles );
    REQUIRE( centerGrid.predictorCount() == expectedCenter.size() );
    for( std::size_t predictor = 0; predictor < expectedCenter.size(); ++predictor )
    {
        REQUIRE( centerGrid.candidatePredictorIndex( predictor ) == expectedCenter[predictor] );
    }

    skyGridT supportSky = candidateGrid( 31, 33, 15.0, 16.0, testRegion( policyT::kernelSupport, 1.25 ) );
    rotatedGridT supportGrid;
    supportGrid.configure( supportSky, angles, nullptr );
    const auto expectedSupport = expectedCandidates( supportSky, angles );
    REQUIRE( supportGrid.predictorCount() == expectedSupport.size() );
    REQUIRE( supportGrid.predictorCount() < centerGrid.predictorCount() );
    for( std::size_t predictor = 0; predictor < expectedSupport.size(); ++predictor )
    {
        REQUIRE( supportGrid.candidatePredictorIndex( predictor ) == expectedSupport[predictor] );
    }

    for( std::size_t search = 0; search < supportGrid.searchPixelCount(); ++search )
    {
        for( std::size_t frame = 0; frame < supportGrid.frameCount(); ++frame )
        {
            const auto target = supportGrid.targetInterpolation( frame, search );
            for( std::size_t predictor = 0; predictor < supportGrid.predictorCount(); ++predictor )
            {
                const auto sample = supportGrid.predictorInterpolation( frame, search, predictor );
                REQUIRE_FALSE( rawFootprintOverlaps( { sample.row(), sample.column() },
                                                     { target.row(), target.column() },
                                                     supportSky.effectiveExclusionRadius() ) );
            }
        }
    }
}

/** \brief Verifies all-frame target and predictor footprints enforce raw common-mask validity.
 *
 * Exercises `P4RotatedGrid::configure()` all-frame `P4RotatedSearchRecord` mask aggregation.
 * \ingroup P4RotatedGrid_unit_tests
 */
TEST_CASE( "P4RotatedGrid validates every raw footprint against the common mask", "[P4RotatedGrid][mask][validity]" )
{
    skyGridT skyGrid = candidateGrid();
    const std::vector<double> angles{ 0.0, 0.293 };
    rotatedGridT baseline;
    baseline.configure( skyGrid, angles, nullptr );
    const std::size_t search = findValidSearch( baseline );

    imageT mask( baseline.rows(), baseline.columns() );
    mask.setOnes();
    const auto targetSecondFrame = baseline.targetInterpolation( 1, search );
    mask( targetSecondFrame.footprintRow(), targetSecondFrame.footprintColumn() ) = 0;
    rotatedGridT targetMasked;
    targetMasked.configure( skyGrid, angles, &mask );
    REQUIRE_FALSE( targetMasked.searchPixel( search ).valid() );
    REQUIRE( targetMasked.searchPixel( search ).invalidReason() != reasonT::none );
    REQUIRE( targetMasked.searchPixel( search ).hasReason( reasonT::targetFootprintMasked ) );
    REQUIRE_FALSE( targetMasked.searchPixel( search ).hasReason( reasonT::targetFootprintOutside ) );

    mask.setOnes();
    const auto predictorSecondFrame = baseline.predictorInterpolation( 1, search, 0 );
    mask( predictorSecondFrame.footprintRow(), predictorSecondFrame.footprintColumn() ) = 0;
    rotatedGridT predictorMasked;
    predictorMasked.configure( skyGrid, angles, &mask );
    REQUIRE_FALSE( predictorMasked.searchPixel( search ).valid() );
    REQUIRE( predictorMasked.searchPixel( search ).hasReason( reasonT::predictorFootprintMasked ) );

    imageT wrongMask( baseline.rows() - 1, baseline.columns() );
    wrongMask.setOnes();
    REQUIRE_THROWS_AS( predictorMasked.configure( skyGrid, angles, &wrongMask ), std::invalid_argument );
    REQUIRE( predictorMasked.frameCount() == angles.size() );
    REQUIRE( predictorMasked.searchPixel( search ).hasReason( reasonT::predictorFootprintMasked ) );

    mask.setOnes();
    mask( 0, 0 ) = std::numeric_limits<float>::quiet_NaN();
    REQUIRE_THROWS_AS( predictorMasked.configure( skyGrid, angles, &mask ), std::invalid_argument );
    REQUIRE( predictorMasked.searchPixel( search ).hasReason( reasonT::predictorFootprintMasked ) );
}

/** \brief Verifies true nominal cubic edge support for odd/even geometry and all-frame invalidation.
 *
 * Exercises target/predictor `P4RotatedInterpolationRecord::footprintInside()` and edge invalidity bits.
 * \ingroup P4RotatedGrid_unit_tests
 */
TEST_CASE( "P4RotatedGrid uses true complete cubic support at raw edges", "[P4RotatedGrid][edge][validity]" )
{
    const regionT edgeRegion( 3.5, 5.5, 1.0, 1.0, 2.0, 45.0, 0.25, policyT::sampleCenter, 0.0 );

    skyGridT oddSky = candidateGrid( 11, 11, 5.0, 5.0, edgeRegion );
    rotatedGridT odd;
    odd.configure( oddSky, { 0.0 }, nullptr );
    const std::size_t completeBorderSearch = findSearch( oddSky, 1, 5 );
    const auto completeBorder = odd.targetInterpolation( 0, completeBorderSearch );
    REQUIRE( completeBorder.footprintRow() == 0 );
    REQUIRE( completeBorder.footprintInside() );
    REQUIRE_FALSE( odd.searchPixel( completeBorderSearch ).hasReason( reasonT::targetFootprintOutside ) );

    const std::size_t incompleteBorderSearch = findSearch( oddSky, 0, 5 );
    const auto incompleteBorder = odd.targetInterpolation( 0, incompleteBorderSearch );
    REQUIRE( incompleteBorder.footprintRow() == -1 );
    REQUIRE_FALSE( incompleteBorder.footprintInside() );
    REQUIRE( odd.searchPixel( incompleteBorderSearch ).hasReason( reasonT::targetFootprintOutside ) );

    skyGridT evenSky = candidateGrid( 12, 14, 5.5, 6.5, edgeRegion );
    rotatedGridT even;
    even.configure( evenSky, { 0.0, 0.5 * std::numbers::pi_v<double> }, nullptr );
    REQUIRE( even.rows() == 12 );
    REQUIRE( even.columns() == 14 );
    REQUIRE( even.frameCount() == 2 );
    REQUIRE( even.searchPixelCount() == evenSky.searchPixelCount() );
    for( std::size_t search = 0; search < even.searchPixelCount(); ++search )
    {
        if( !even.targetInterpolation( 0, search ).footprintInside() ||
            !even.targetInterpolation( 1, search ).footprintInside() )
        {
            REQUIRE( even.searchPixel( search ).hasReason( reasonT::targetFootprintOutside ) );
        }
    }
}

/** \brief Verifies configuration and sampling errors preserve valid state and caller outputs transactionally.
 *
 * Exercises `P4RotatedGrid::configure()` validation and finite-value checks in all sampling APIs.
 * \ingroup P4RotatedGrid_unit_tests
 */
TEST_CASE( "P4RotatedGrid rejects inconsistent and non-finite inputs transactionally",
           "[P4RotatedGrid][error][transaction]" )
{
    rotatedGridT empty;
    imageT image( 31, 33 );
    image.setOnes();
    REQUIRE_FALSE( empty.configured() );
    REQUIRE_THROWS_AS( empty.sampleTarget( image, 0, 0 ), std::logic_error );
    REQUIRE_THROWS_AS( empty.targetInterpolation( 0, 0 ), std::logic_error );
    REQUIRE_THROWS_AS( empty.predictorInterpolation( 0, 0, 0 ), std::logic_error );

    skyGridT resizedOnly;
    resizedOnly.resize( 31, 33 );
    REQUIRE_THROWS_AS( empty.configure( resizedOnly, { 0.0 }, nullptr ), std::logic_error );

    skyGridT skyGrid = candidateGrid();
    REQUIRE_THROWS_AS( empty.configure( skyGrid, {}, nullptr ), std::invalid_argument );
    REQUIRE_THROWS_AS( empty.configure( skyGrid, { std::numeric_limits<double>::quiet_NaN() }, nullptr ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( empty.configure( skyGrid, { std::numeric_limits<double>::infinity() }, nullptr ),
                       std::invalid_argument );

    skyGridT fullyExcludedSky = candidateGrid( 31, 33, 15.0, 16.0, testRegion( policyT::sampleCenter, 100.0 ) );
    REQUIRE_THROWS_AS( empty.configure( fullyExcludedSky, { 0.0 }, nullptr ), std::invalid_argument );

    REQUIRE( mx::improc::P4RotatedGridTestAccess::checkedCoordinateCount( 3, 4 ) == 12 );
    REQUIRE_THROWS_AS(
        mx::improc::P4RotatedGridTestAccess::checkedCoordinateCount( 2, std::numeric_limits<std::size_t>::max() ),
        std::length_error );
    REQUIRE_THROWS_AS( mx::improc::P4RotatedGridTestAccess::mapFootprint( 5,
                                                                          5,
                                                                          std::numeric_limits<double>::infinity(),
                                                                          2.0,
                                                                          1.0,
                                                                          0.0,
                                                                          0.0,
                                                                          0.0 ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( mx::improc::P4RotatedGridTestAccess::mapFootprint(
                           5,
                           5,
                           0.0,
                           0.0,
                           1.0,
                           0.0,
                           static_cast<double>( std::numeric_limits<int>::max() ) + 16.0,
                           0.0 ),
                       std::invalid_argument );
    imageT incompleteImage( 5, 5 );
    incompleteImage.setOnes();
    REQUIRE_THROWS_AS( mx::improc::P4RotatedGridTestAccess::sampleIncomplete( incompleteImage ), std::logic_error );

    rotatedGridT grid;
    grid.configure( skyGrid, { 0.0, 0.293 }, nullptr );
    const std::size_t search = findValidSearch( grid );
    const std::size_t originalFrames = grid.frameCount();
    const std::size_t originalPredictors = grid.predictorCount();
    REQUIRE_THROWS_AS( grid.configure( skyGrid, {}, nullptr ), std::invalid_argument );
    REQUIRE( grid.frameCount() == originalFrames );
    REQUIRE( grid.predictorCount() == originalPredictors );

    imageT wrongImage( grid.rows(), grid.columns() - 1 );
    wrongImage.setOnes();
    REQUIRE_THROWS_AS( grid.sampleTarget( wrongImage, 0, search ), std::invalid_argument );
    REQUIRE_THROWS_AS( grid.sampleTarget( image, grid.frameCount(), search ), std::out_of_range );
    REQUIRE_THROWS_AS( grid.samplePredictor( image, 0, search, grid.predictorCount() ), std::out_of_range );

    float target = 123;
    std::vector<float> wrongPredictors( grid.predictorCount() - 1, 456 );
    REQUIRE_THROWS_AS( grid.sampleFrame( image, 0, search, target, wrongPredictors ), std::invalid_argument );
    REQUIRE( target == 123 );
    REQUIRE(
        std::all_of( wrongPredictors.begin(), wrongPredictors.end(), []( float value ) { return value == 456; } ) );

    std::vector<float> predictors( grid.predictorCount(), 456 );
    const auto targetRecord = grid.targetInterpolation( 0, search );
    image( targetRecord.footprintRow(), targetRecord.footprintColumn() ) = std::numeric_limits<float>::quiet_NaN();
    REQUIRE_THROWS_AS( grid.sampleFrame( image, 0, search, target, predictors ), std::domain_error );
    REQUIRE( target == 123 );
    REQUIRE( std::all_of( predictors.begin(), predictors.end(), []( float value ) { return value == 456; } ) );

    image.setZero();
    const auto fractionalRecord = grid.targetInterpolation( 1, search );
    for( int rowOffset = 0; rowOffset < rotatedGridT::width; ++rowOffset )
    {
        for( int columnOffset = 0; columnOffset < rotatedGridT::width; ++columnOffset )
        {
            image( fractionalRecord.footprintRow() + rowOffset, fractionalRecord.footprintColumn() + columnOffset ) =
                fractionalRecord.kernel()( rowOffset, columnOffset ) < 0 ? -std::numeric_limits<float>::max()
                                                                         : std::numeric_limits<float>::max();
        }
    }
    REQUIRE_THROWS_AS( grid.sampleTarget( image, 1, search ), std::domain_error );

    image.setOnes();
    std::size_t invalidSearch = grid.searchPixelCount();
    const regionT edgeRegion( 3.5, 5.5, 1.0, 1.0, 2.0, 45.0, 0.25, policyT::sampleCenter, 0.0 );
    skyGridT edgeSky = candidateGrid( 11, 11, 5.0, 5.0, edgeRegion );
    rotatedGridT edgeGrid;
    edgeGrid.configure( edgeSky, { 0.0 }, nullptr );
    for( std::size_t candidate = 0; candidate < edgeGrid.searchPixelCount(); ++candidate )
    {
        if( !edgeGrid.searchPixel( candidate ).valid() )
        {
            invalidSearch = candidate;
            break;
        }
    }
    REQUIRE( invalidSearch < edgeGrid.searchPixelCount() );
    imageT edgeImage( 11, 11 );
    edgeImage.setOnes();
    REQUIRE_THROWS_AS( edgeGrid.sampleTarget( edgeImage, 0, invalidSearch ), std::logic_error );

    REQUIRE_THROWS_AS(
        mx::improc::P4RotatedSearchRecord( mx::improc::P4PixelCoordinate( 1, 2 ), static_cast<reasonT>( 16 ) ),
        std::invalid_argument );
}

} // namespace P4RotatedGrid_test
} // namespace unitTest
