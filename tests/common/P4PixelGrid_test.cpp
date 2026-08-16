/** \file P4PixelGrid_test.cpp
 * \brief Tests the local pixel geometry used by Pixel Prediction Post-Processing.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "src/common/P4PixelGrid.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include <omp.h>

namespace mx
{
namespace improc
{

/** \cond P4PixelGrid_test_harness */
class P4PixelGridTestAccess
{
  public:
    /// Invoke the production mapped-coordinate range check.
    static int checkedFootprintOrigin( double mappedCoordinate /**< [in] mapped coordinate to check */ )
    {
        return P4PixelGridf::checkedFootprintOrigin( mappedCoordinate );
    }

    /// Invoke the production interpolation-count overflow check.
    static std::size_t checkedInterpolationCount( std::size_t searchPixelCount, /**< [in] search-pixel count */
                                                  std::size_t predictorCount /**< [in] predictor count */ )
    {
        return P4PixelGridf::checkedInterpolationCount( searchPixelCount, predictorCount );
    }
};
/** \endcond */

} // namespace improc
} // namespace mx

namespace unitTest
{
namespace P4PixelGrid_test
{

/** \cond P4PixelGrid_test_harness */
using gridT = mx::improc::P4PixelGridf;
using regionT = mx::improc::P4PixelGridRegion;
using policyT = mx::improc::P4ExclusionPolicy;
using reasonT = mx::improc::P4PixelInvalidReason;

/// Return a compact region that remains safely inside a 41-by-43 image.
regionT safeRegion( policyT policy = policyT::sampleCenter, /**< [in] central-exclusion policy */
                    double psfRadius = 0.5,                 /**< [in] physical exclusion radius */
                    double radiusBuffer = 0.0 /**< [in] explicit exclusion comparison buffer */ )
{
    return regionT( 5.0, 6.0, 2.0, 2.0, 3.0, 45.0, psfRadius, policy, radiusBuffer );
}

/// Find one search coordinate, requiring that it occurs exactly once.
std::size_t findSearch( const gridT &grid, /**< [in] configured grid */
                        int row,           /**< [in] desired row */
                        int column /**< [in] desired column */ )
{
    std::size_t found = grid.searchPixelCount();
    for( std::size_t index = 0; index < grid.searchPixelCount(); ++index )
    {
        if( grid.searchPixel( index ).coordinate().row() == row &&
            grid.searchPixel( index ).coordinate().column() == column )
        {
            REQUIRE( found == grid.searchPixelCount() );
            found = index;
        }
    }
    REQUIRE( found < grid.searchPixelCount() );
    return found;
}

/// Return whether two grids have identical public geometry and kernels.
void requireSameGrid( const gridT &left, /**< [in] first grid */
                      const gridT &right /**< [in] second grid */ )
{
    REQUIRE( left.searchPixelCount() == right.searchPixelCount() );
    REQUIRE( left.predictorCount() == right.predictorCount() );
    for( std::size_t searchIndex = 0; searchIndex < left.searchPixelCount(); ++searchIndex )
    {
        const auto &leftSearch = left.searchPixel( searchIndex );
        const auto &rightSearch = right.searchPixel( searchIndex );
        REQUIRE( leftSearch.coordinate().row() == rightSearch.coordinate().row() );
        REQUIRE( leftSearch.coordinate().column() == rightSearch.coordinate().column() );
        REQUIRE( leftSearch.invalidReason() == rightSearch.invalidReason() );

        for( std::size_t predictorIndex = 0; predictorIndex < left.predictorCount(); ++predictorIndex )
        {
            const auto &leftInterpolation = left.interpolation( searchIndex, predictorIndex );
            const auto &rightInterpolation = right.interpolation( searchIndex, predictorIndex );
            REQUIRE( leftInterpolation.row() == rightInterpolation.row() );
            REQUIRE( leftInterpolation.column() == rightInterpolation.column() );
            REQUIRE( leftInterpolation.footprintRow() == rightInterpolation.footprintRow() );
            REQUIRE( leftInterpolation.footprintColumn() == rightInterpolation.footprintColumn() );
            REQUIRE( ( leftInterpolation.kernel() == rightInterpolation.kernel() ).all() );
        }
    }
}

/// Independently enumerate canonical sample-center offsets for one region.
std::set<std::pair<int, int>>
expectedSampleCenterOffsets( double xCenter, /**< [in] grid center row */
                             double yCenter, /**< [in] grid center column */
                             const regionT &configuration /**< [in] region to enumerate */ )
{
    const double middleRadius =
        configuration.searchInnerRadius + 0.5 * ( configuration.searchOuterRadius - configuration.searchInnerRadius );
    const double referenceRow = xCenter + 0.5;
    const double referenceColumn = yCenter + 0.5;
    const double canonicalRow = std::floor( referenceRow );
    const double canonicalColumn = std::floor( referenceColumn + middleRadius );
    const double deltaRow = canonicalRow - referenceRow;
    const double deltaColumn = canonicalColumn - referenceColumn;
    const double canonicalAngle = std::atan2( deltaColumn, deltaRow ) * 180.0 / std::numbers::pi_v<double>;
    const double innerRadius = deltaColumn - configuration.optimizationDeltaRadiusInner;
    const double outerRadius = deltaColumn + configuration.optimizationDeltaRadiusOuter;
    const double halfAngle =
        std::min( configuration.optimizationArcHalfWidth / ( 2.0 * std::numbers::pi_v<double> * middleRadius ) * 360.0,
                  configuration.optimizationMaxHalfAngle );
    const double effectiveRadius = configuration.psfRadius + configuration.exclusionRadiusBuffer;
    const double angularTolerance =
        16.0 * std::numeric_limits<double>::epsilon() * std::max( { 1.0, std::abs( canonicalAngle ), halfAngle } );
    const int offsetLimit = static_cast<int>( std::ceil( std::abs( deltaColumn ) + outerRadius + 2.0 ) );

    std::set<std::pair<int, int>> offsets;
    for( int columnOffset = -offsetLimit; columnOffset <= offsetLimit; ++columnOffset )
    {
        for( int rowOffset = -offsetLimit; rowOffset <= offsetLimit; ++rowOffset )
        {
            const double radius = std::hypot( rowOffset + deltaRow, columnOffset + deltaColumn );
            if( radius < innerRadius - 0.5 || radius >= outerRadius + 0.5 )
            {
                continue;
            }

            double angle =
                std::atan2( columnOffset + deltaColumn, rowOffset + deltaRow ) * 180.0 / std::numbers::pi_v<double>;
            if( angle < 0 )
            {
                angle += 360.0;
            }
            double difference = angle - canonicalAngle;
            while( difference < -180.0 )
            {
                difference += 360.0;
            }
            while( difference >= 180.0 )
            {
                difference -= 360.0;
            }
            if( std::abs( difference ) > halfAngle + angularTolerance )
            {
                continue;
            }

            if( std::hypot( static_cast<double>( rowOffset ), static_cast<double>( columnOffset ) ) <= effectiveRadius )
            {
                continue;
            }
            offsets.emplace( rowOffset, columnOffset );
        }
    }
    return offsets;
}

/// Restore the process-wide OpenMP thread count after one focused test.
class OmpThreadCountGuard
{
  public:
    /// Capture the current OpenMP maximum thread count.
    OmpThreadCountGuard() : m_threadCount( omp_get_max_threads() )
    {
    }

    /// Restore the captured OpenMP thread count.
    ~OmpThreadCountGuard()
    {
        omp_set_num_threads( m_threadCount );
    }

    OmpThreadCountGuard( const OmpThreadCountGuard & ) = delete;
    OmpThreadCountGuard &operator=( const OmpThreadCountGuard & ) = delete;

  private:
    int m_threadCount; ///< Captured maximum thread count.
};
/** \endcond */

/** \brief Verifies P4PixelGrid resize validation, default/explicit centers, and transactional state replacement.
 *
 * Exercises the real `P4PixelGrid::resize()` API.
 * \ingroup P4PixelGrid_unit_tests
 */
TEST_CASE( "P4PixelGrid validates and replaces image geometry", "[P4PixelGrid][resize]" )
{
#ifdef __DOXY_ONLY__
    // clang-format off
    using DoxygenGrid = mx::improc::P4PixelGrid<mx::improc::cubicConvolTransform<float>>;
    DoxygenGrid doxygenGrid;
    void (DoxygenGrid::*doxygenExplicitResize)(int, int, double, double) = &DoxygenGrid::resize;
    doxygenGrid.resize( 4, 5, 1.0, 2.0 );
    // clang-format on
#endif

    mx::improc::P4PixelGridf grid;
    REQUIRE_FALSE( grid.resized() );
    REQUIRE_FALSE( grid.regionConfigured() );
    REQUIRE_THROWS_AS( grid.regionConfiguration(), std::logic_error );
    REQUIRE_THROWS_AS( grid.effectiveExclusionRadius(), std::logic_error );

    REQUIRE_THROWS_AS( grid.resize( 0, 4 ), std::invalid_argument );
    REQUIRE_THROWS_AS( grid.resize( 4, 0 ), std::invalid_argument );
    REQUIRE_THROWS_AS( grid.resize( 4, 5, std::numeric_limits<double>::quiet_NaN(), 2.0 ), std::invalid_argument );
    REQUIRE_THROWS_AS( grid.resize( 4, 5, 1.0, std::numeric_limits<double>::infinity() ), std::invalid_argument );
    REQUIRE_THROWS_AS( grid.resize( 4, 5, -0.1, 2.0 ), std::invalid_argument );
    REQUIRE_THROWS_AS( grid.resize( 4, 5, 1.0, 5.0 ), std::invalid_argument );

    grid.resize( 40, 42 );
    REQUIRE( grid.resized() );
    REQUIRE( grid.rows() == 40 );
    REQUIRE( grid.columns() == 42 );
    REQUIRE( grid.xCenter() == 19.5 );
    REQUIRE( grid.yCenter() == 20.5 );
    REQUIRE( grid.searchPixelCount() == 0 );
    REQUIRE( grid.predictorCount() == 0 );

    grid.resize( 41, 43, 20.25, 21.75 );
    REQUIRE( grid.rows() == 41 );
    REQUIRE( grid.columns() == 43 );
    REQUIRE( grid.xCenter() == 20.25 );
    REQUIRE( grid.yCenter() == 21.75 );

    REQUIRE_THROWS_AS( grid.resize( -1, 43 ), std::invalid_argument );
    REQUIRE( grid.rows() == 41 );
    REQUIRE( grid.columns() == 43 );
    REQUIRE( grid.xCenter() == 20.25 );
    REQUIRE( grid.yCenter() == 21.75 );
}

/** \brief Verifies P4PixelGrid search membership, half-pixel canonical offsets, and cardinal mapping.
 *
 * Exercises the real `P4PixelGrid::region()`, search-pixel, predictor-offset, and interpolation APIs.
 * \ingroup P4PixelGrid_unit_tests
 */
TEST_CASE( "P4PixelGrid preserves annulus and half-pixel mapping geometry", "[P4PixelGrid][geometry]" )
{
#ifdef __DOXY_ONLY__
    // clang-format off
    using DoxygenTransform = mx::improc::cubicConvolTransform<float>;
    mx::improc::P4PixelCoordinate doxygenCoordinate( 0, 0 );
    doxygenCoordinate.row();
    doxygenCoordinate.column();
    mx::improc::P4PixelGrid<DoxygenTransform> doxygenGrid;
    const mx::improc::P4SearchPixelRecord &doxygenSearch = doxygenGrid.searchPixel( 0 );
    doxygenSearch.coordinate();
    doxygenSearch.invalidReason();
    doxygenSearch.valid();
    doxygenSearch.targetMasked();
    doxygenSearch.predictorMasked();
    const mx::improc::P4InterpolationRecord<DoxygenTransform> &doxygenInterpolation = doxygenGrid.interpolation( 0, 0 );
    doxygenInterpolation.row();
    doxygenInterpolation.column();
    doxygenInterpolation.footprintRow();
    doxygenInterpolation.footprintColumn();
    doxygenInterpolation.kernel();
    // clang-format on
#endif

    mx::improc::P4PixelGridf grid;
    grid.resize( 41, 43, 20.0, 21.0 );
    const mx::improc::P4PixelGridRegion configuration = safeRegion();
    grid.region( configuration, nullptr );

    REQUIRE( grid.regionConfigured() );
    REQUIRE( &grid.regionConfiguration() != &configuration );
    REQUIRE( grid.regionConfiguration().searchInnerRadius == 5.0 );
    REQUIRE( grid.regionConfiguration().searchOuterRadius == 6.0 );
    REQUIRE( grid.regionConfiguration().optimizationDeltaRadiusInner == 2.0 );
    REQUIRE( grid.regionConfiguration().optimizationDeltaRadiusOuter == 2.0 );
    REQUIRE( grid.regionConfiguration().optimizationArcHalfWidth == 3.0 );
    REQUIRE( grid.regionConfiguration().optimizationMaxHalfAngle == 45.0 );
    REQUIRE( grid.regionConfiguration().psfRadius == 0.5 );
    REQUIRE( grid.regionConfiguration().exclusionPolicy == policyT::sampleCenter );
    REQUIRE( grid.regionConfiguration().exclusionRadiusBuffer == 0.0 );
    REQUIRE( grid.effectiveExclusionRadius() == 0.5 );
    REQUIRE( grid.searchPixelCount() > 4 );
    REQUIRE( grid.predictorCount() > 0 );

    for( std::size_t index = 0; index < grid.searchPixelCount(); ++index )
    {
        const auto &coordinate = grid.searchPixel( index ).coordinate();
        const double radius = std::hypot( static_cast<double>( coordinate.row() ) - 20.0,
                                          static_cast<double>( coordinate.column() ) - 21.0 );
        REQUIRE( radius >= 5.0 );
        REQUIRE( radius < 6.0 );
        REQUIRE( grid.searchPixel( index ).valid() );
        REQUIRE( grid.searchPixel( index ).invalidReason() == reasonT::none );
        REQUIRE_FALSE( grid.searchPixel( index ).targetMasked() );
        REQUIRE_FALSE( grid.searchPixel( index ).predictorMasked() );
    }

    for( std::size_t predictorIndex = 0; predictorIndex < grid.predictorCount(); ++predictorIndex )
    {
        const auto &offset = grid.predictorOffset( predictorIndex );
        REQUIRE( std::hypot( static_cast<double>( offset.row() ), static_cast<double>( offset.column() ) ) > 0.5 );
    }

    const std::size_t north = findSearch( grid, 20, 26 );
    const std::size_t east = findSearch( grid, 25, 21 );
    const std::size_t south = findSearch( grid, 20, 16 );
    const std::size_t west = findSearch( grid, 15, 21 );
    const std::size_t offAxis = findSearch( grid, 23, 25 );
    static_cast<void>( offAxis );

    for( std::size_t predictorIndex = 0; predictorIndex < grid.predictorCount(); ++predictorIndex )
    {
        const auto &offset = grid.predictorOffset( predictorIndex );
        const auto &northMap = grid.interpolation( north, predictorIndex );
        const auto &eastMap = grid.interpolation( east, predictorIndex );
        const auto &southMap = grid.interpolation( south, predictorIndex );
        const auto &westMap = grid.interpolation( west, predictorIndex );
        const auto &offAxisMap = grid.interpolation( offAxis, predictorIndex );

        REQUIRE( northMap.row() == Approx( 20.0 + offset.row() ).margin( 1e-12 ) );
        REQUIRE( northMap.column() == Approx( 26.0 + offset.column() ).margin( 1e-12 ) );
        REQUIRE( eastMap.row() == Approx( 25.0 + offset.column() ).margin( 1e-12 ) );
        REQUIRE( eastMap.column() == Approx( 21.0 - offset.row() ).margin( 1e-12 ) );
        REQUIRE( southMap.row() == Approx( 20.0 - offset.row() ).margin( 1e-12 ) );
        REQUIRE( southMap.column() == Approx( 16.0 - offset.column() ).margin( 1e-12 ) );
        REQUIRE( westMap.row() == Approx( 15.0 - offset.column() ).margin( 1e-12 ) );
        REQUIRE( westMap.column() == Approx( 21.0 + offset.row() ).margin( 1e-12 ) );
        REQUIRE( offAxisMap.row() == Approx( 23.0 + 0.8 * offset.row() + 0.6 * offset.column() ).margin( 1e-12 ) );
        REQUIRE( offAxisMap.column() == Approx( 25.0 - 0.6 * offset.row() + 0.8 * offset.column() ).margin( 1e-12 ) );
    }

    REQUIRE_THROWS_AS( grid.searchPixel( grid.searchPixelCount() ), std::out_of_range );
    REQUIRE_THROWS_AS( grid.predictorOffset( grid.predictorCount() ), std::out_of_range );
    REQUIRE_THROWS_AS( grid.interpolation( grid.searchPixelCount(), 0 ), std::out_of_range );
    REQUIRE_THROWS_AS( grid.interpolation( 0, grid.predictorCount() ), std::out_of_range );
}

/** \brief Verifies P4PixelGrid supports odd/even rectangles and clears stale records after reconfiguration.
 *
 * Exercises successful calls to the real P4PixelGrid resize and region APIs.
 * \ingroup P4PixelGrid_unit_tests
 */
TEST_CASE( "P4PixelGrid supports rectangular geometry and reconfiguration", "[P4PixelGrid][geometry]" )
{
    mx::improc::P4PixelGridf oddGrid;
    oddGrid.resize( 41, 43 );
    oddGrid.region( safeRegion(), nullptr );
    REQUIRE( oddGrid.xCenter() == 20.0 );
    REQUIRE( oddGrid.yCenter() == 21.0 );
    REQUIRE( oddGrid.searchPixelCount() > 0 );
    REQUIRE( oddGrid.predictorCount() > 0 );

    mx::improc::P4PixelGridf evenGrid;
    evenGrid.resize( 40, 42 );
    evenGrid.region( safeRegion(), nullptr );
    REQUIRE( evenGrid.xCenter() == 19.5 );
    REQUIRE( evenGrid.yCenter() == 20.5 );
    REQUIRE( evenGrid.searchPixelCount() > 0 );
    REQUIRE( evenGrid.predictorCount() > 0 );

    const std::size_t firstSearchCount = oddGrid.searchPixelCount();
    const std::size_t firstPredictorCount = oddGrid.predictorCount();
    const double firstMappedRow = oddGrid.interpolation( 0, 0 ).row();
    const mx::improc::P4PixelGridRegion changed( 7.0, 8.0, 3.0, 3.0, 1.0, 20.0, 0.5, policyT::sampleCenter, 0.0 );
    oddGrid.region( changed, nullptr );
    REQUIRE( oddGrid.regionConfiguration().searchInnerRadius == 7.0 );
    REQUIRE( oddGrid.regionConfiguration().searchOuterRadius == 8.0 );
    REQUIRE( oddGrid.searchPixelCount() != firstSearchCount );
    REQUIRE( oddGrid.predictorCount() != firstPredictorCount );
    REQUIRE( oddGrid.interpolation( 0, 0 ).row() != firstMappedRow );
    for( std::size_t searchIndex = 0; searchIndex < oddGrid.searchPixelCount(); ++searchIndex )
    {
        const auto &coordinate = oddGrid.searchPixel( searchIndex ).coordinate();
        const double radius = std::hypot( static_cast<double>( coordinate.row() ) - oddGrid.xCenter(),
                                          static_cast<double>( coordinate.column() ) - oddGrid.yCenter() );
        REQUIRE( radius >= 7.0 );
        REQUIRE( radius < 8.0 );
    }

    oddGrid.resize( 45, 47, 22.25, 23.5 );
    REQUIRE_FALSE( oddGrid.regionConfigured() );
    REQUIRE( oddGrid.searchPixelCount() == 0 );
    REQUIRE( oddGrid.predictorCount() == 0 );
}

/** \brief Verifies P4PixelGrid optimization membership for a capped half-angle that wraps through zero degrees.
 *
 * Exercises exact canonical radial/azimuthal boundaries in the real P4PixelGrid region API.
 * \ingroup P4PixelGrid_unit_tests
 */
TEST_CASE( "P4PixelGrid caps and wraps optimization wedges", "[P4PixelGrid][geometry][wrap]" )
{
    const mx::improc::P4PixelGridRegion capped( 4.0, 6.0, 6.0, 4.0, 100.0, 110.0, 0.5, policyT::sampleCenter, 0.5 );
    mx::improc::P4PixelGridf grid;
    grid.resize( 61, 65, 30.25, 31.75 );
    grid.region( capped, nullptr );

    const std::set<std::pair<int, int>> expected =
        expectedSampleCenterOffsets( grid.xCenter(), grid.yCenter(), capped );
    std::set<std::pair<int, int>> actual;
    for( std::size_t predictorIndex = 0; predictorIndex < grid.predictorCount(); ++predictorIndex )
    {
        const auto &offset = grid.predictorOffset( predictorIndex );
        actual.emplace( offset.row(), offset.column() );
    }
    REQUIRE( actual == expected );

    const double referenceRow = grid.xCenter() + 0.5;
    const double referenceColumn = grid.yCenter() + 0.5;
    const double canonicalRow = std::floor( referenceRow );
    const double canonicalColumn = std::floor( referenceColumn + 5.0 );
    bool belowTenDegrees{ false };
    bool aboveThreeHundredFiftyDegrees{ false };
    for( const auto &offset : actual )
    {
        double angle = std::atan2( offset.second + canonicalColumn - referenceColumn,
                                   offset.first + canonicalRow - referenceRow ) *
                       180.0 / std::numbers::pi_v<double>;
        if( angle < 0 )
        {
            angle += 360.0;
        }
        belowTenDegrees = belowTenDegrees || angle < 10.0;
        aboveThreeHundredFiftyDegrees = aboveThreeHundredFiftyDegrees || angle > 350.0;
    }
    REQUIRE( belowTenDegrees );
    REQUIRE( aboveThreeHundredFiftyDegrees );

    const mx::improc::P4PixelGridRegion narrow( 4.0, 6.0, 6.0, 4.0, 2.0, 110.0, 0.5, policyT::sampleCenter, 0.5 );
    mx::improc::P4PixelGridf narrowGrid;
    narrowGrid.resize( 61, 65, 30.25, 31.75 );
    narrowGrid.region( narrow, nullptr );
    REQUIRE( narrowGrid.predictorCount() < grid.predictorCount() );
}

/** \brief Verifies P4PixelGrid preserves all AF Lep predictor counts and both first-annulus angular boundary rays.
 *
 * Exercises the real P4PixelGrid region and predictor-offset APIs with the prototype comparison geometry.
 * \ingroup P4PixelGrid_unit_tests
 */
TEST_CASE( "P4PixelGrid preserves AF Lep predictor counts and wedge boundaries", "[P4PixelGrid][geometry][aflep]" )
{
    mx::improc::P4PixelGridf grid;
    grid.resize( 128, 128 );

    constexpr std::array<double, 7> innerRadii{ 6, 8, 10, 12, 14, 16, 18 };
    constexpr std::array<double, 7> outerRadii{ 8, 10, 12, 14, 16, 18, 20 };
    constexpr std::array<std::size_t, 7> expectedPredictorCounts{ 2235, 2478, 2724, 2970, 3212, 3458, 3704 };

    for( std::size_t regionIndex = 0; regionIndex < innerRadii.size(); ++regionIndex )
    {
        grid.region(
            regionT( innerRadii[regionIndex], outerRadii[regionIndex], 8, 30, 32, 90, 1.5, policyT::sampleCenter, 0.5 ),
            nullptr );
        REQUIRE( grid.predictorCount() == expectedPredictorCounts[regionIndex] );
    }

    grid.region( regionT( 6, 8, 8, 30, 32, 90, 1.5, policyT::sampleCenter, 0.5 ), nullptr );

    std::set<std::pair<int, int>> offsets;
    for( std::size_t predictorIndex = 0; predictorIndex < grid.predictorCount(); ++predictorIndex )
    {
        const auto &offset = grid.predictorOffset( predictorIndex );
        offsets.emplace( offset.row(), offset.column() );
    }

    for( int radius = 1; radius <= 37; ++radius )
    {
        REQUIRE( offsets.count( { radius, -7 } ) == 1 );
        REQUIRE( offsets.count( { -radius, -7 } ) == 1 );
    }
}

/** \brief Verifies P4PixelGrid constant and linear-ramp sampling at cardinal and subpixel positions.
 *
 * Exercises the real `P4PixelGrid::sample()` API and its precomputed cubic kernels.
 * \ingroup P4PixelGrid_unit_tests
 */
TEST_CASE( "P4PixelGrid samples constant and linear images", "[P4PixelGrid][sample]" )
{
    mx::improc::P4PixelGridf grid;
    grid.resize( 41, 43, 20.0, 21.0 );
    grid.region( safeRegion(), nullptr );

    gridT::imageT constant = gridT::imageT::Constant( 41, 43, 3.25F );
    gridT::imageT ramp( 41, 43 );
    for( int row = 0; row < ramp.rows(); ++row )
    {
        for( int column = 0; column < ramp.cols(); ++column )
        {
            ramp( row, column ) = 2.0F * static_cast<float>( row ) + 3.0F * static_cast<float>( column ) + 7.0F;
        }
    }

    const std::vector<std::size_t> searchIndices{ findSearch( grid, 20, 26 ),
                                                  findSearch( grid, 25, 21 ),
                                                  findSearch( grid, 20, 16 ),
                                                  findSearch( grid, 15, 21 ),
                                                  findSearch( grid, 23, 25 ) };
    for( const std::size_t searchIndex : searchIndices )
    {
        for( std::size_t predictorIndex = 0; predictorIndex < grid.predictorCount(); ++predictorIndex )
        {
            const auto &record = grid.interpolation( searchIndex, predictorIndex );
            REQUIRE( record.kernel().sum() == Approx( 1.0F ).margin( 2e-6 ) );
            REQUIRE( grid.sample( constant, searchIndex, predictorIndex ) == Approx( 3.25F ).margin( 2e-5 ) );
            const float expected = static_cast<float>( 2.0 * record.row() + 3.0 * record.column() + 7.0 );
            REQUIRE( grid.sample( ramp, searchIndex, predictorIndex ) == Approx( expected ).margin( 2e-4 ) );
        }
    }

    gridT::imageT wrongSize( 40, 43 );
    REQUIRE_THROWS_AS( grid.sample( wrongSize, 0, 0 ), std::invalid_argument );
    REQUIRE_THROWS_AS( grid.sample( constant, grid.searchPixelCount(), 0 ), std::out_of_range );
    REQUIRE_THROWS_AS( grid.sample( constant, 0, grid.predictorCount() ), std::out_of_range );
}

/** \brief Verifies both P4PixelGrid central-exclusion policies and their shared explicit radius buffer.
 *
 * Exercises the real `P4PixelGrid::region()` API with `sampleCenter` and `kernelSupport`.
 * \ingroup P4PixelGrid_unit_tests
 */
TEST_CASE( "P4PixelGrid applies explicit center and whole-kernel exclusion", "[P4PixelGrid][exclusion]" )
{
    mx::improc::P4PixelGridf sampleGrid;
    sampleGrid.resize( 41, 43, 20.0, 21.0 );
    sampleGrid.region( safeRegion( policyT::sampleCenter, 0.5, 0.0 ), nullptr );

    bool sampleFootprintOverlapsTarget{ false };
    for( std::size_t searchIndex = 0; searchIndex < sampleGrid.searchPixelCount(); ++searchIndex )
    {
        const auto &target = sampleGrid.searchPixel( searchIndex ).coordinate();
        for( std::size_t predictorIndex = 0; predictorIndex < sampleGrid.predictorCount(); ++predictorIndex )
        {
            const auto &record = sampleGrid.interpolation( searchIndex, predictorIndex );
            REQUIRE_FALSE( ( record.row() == target.row() && record.column() == target.column() ) );
            sampleFootprintOverlapsTarget =
                sampleFootprintOverlapsTarget ||
                ( record.footprintRow() <= target.row() && target.row() < record.footprintRow() + gridT::width &&
                  record.footprintColumn() <= target.column() &&
                  target.column() < record.footprintColumn() + gridT::width );
        }
    }
    REQUIRE( sampleFootprintOverlapsTarget );

    mx::improc::P4PixelGridf bufferedSampleGrid;
    bufferedSampleGrid.resize( 41, 43, 20.0, 21.0 );
    bufferedSampleGrid.region( safeRegion( policyT::sampleCenter, 0.5, 0.5 ), nullptr );
    REQUIRE( bufferedSampleGrid.effectiveExclusionRadius() == 1.0 );
    REQUIRE( bufferedSampleGrid.predictorCount() <= sampleGrid.predictorCount() );
    for( std::size_t predictorIndex = 0; predictorIndex < bufferedSampleGrid.predictorCount(); ++predictorIndex )
    {
        const auto &offset = bufferedSampleGrid.predictorOffset( predictorIndex );
        REQUIRE( std::hypot( static_cast<double>( offset.row() ), static_cast<double>( offset.column() ) ) > 1.0 );
    }
    for( std::size_t searchIndex = 0; searchIndex < bufferedSampleGrid.searchPixelCount(); ++searchIndex )
    {
        const auto &target = bufferedSampleGrid.searchPixel( searchIndex ).coordinate();
        for( std::size_t predictorIndex = 0; predictorIndex < bufferedSampleGrid.predictorCount(); ++predictorIndex )
        {
            const auto &record = bufferedSampleGrid.interpolation( searchIndex, predictorIndex );
            REQUIRE( std::hypot( record.row() - target.row(), record.column() - target.column() ) >
                     bufferedSampleGrid.effectiveExclusionRadius() );
        }
    }

    mx::improc::P4PixelGridf supportGrid;
    supportGrid.resize( 41, 43, 20.0, 21.0 );
    supportGrid.region( safeRegion( policyT::kernelSupport, 0.5, 0.0 ), nullptr );
    REQUIRE( supportGrid.predictorCount() < sampleGrid.predictorCount() );
    REQUIRE( supportGrid.searchPixelCount() == sampleGrid.searchPixelCount() );

    const double radiusSquared = supportGrid.effectiveExclusionRadius() * supportGrid.effectiveExclusionRadius();
    for( std::size_t searchIndex = 0; searchIndex < supportGrid.searchPixelCount(); ++searchIndex )
    {
        const auto &target = supportGrid.searchPixel( searchIndex ).coordinate();
        for( std::size_t predictorIndex = 0; predictorIndex < supportGrid.predictorCount(); ++predictorIndex )
        {
            const auto &record = supportGrid.interpolation( searchIndex, predictorIndex );
            for( int rowOffset = 0; rowOffset < gridT::width; ++rowOffset )
            {
                for( int columnOffset = 0; columnOffset < gridT::width; ++columnOffset )
                {
                    const double deltaRow = record.footprintRow() + rowOffset - target.row();
                    const double deltaColumn = record.footprintColumn() + columnOffset - target.column();
                    REQUIRE( deltaRow * deltaRow + deltaColumn * deltaColumn > radiusSquared );
                }
            }
        }
    }

    mx::improc::P4PixelGridf prototypeRadiusGrid;
    prototypeRadiusGrid.resize( 41, 43, 20.0, 21.0 );
    prototypeRadiusGrid.region( safeRegion( policyT::sampleCenter, 1.5, 0.5 ), nullptr );
    REQUIRE( prototypeRadiusGrid.effectiveExclusionRadius() == 2.0 );
}

/** \brief Verifies P4PixelGrid preserves annulus-wide K while reporting all common-mask invalidity reasons.
 *
 * Exercises the real `P4PixelGrid::region()` and `P4SearchPixelRecord` validity APIs.
 * \ingroup P4PixelGrid_unit_tests
 */
TEST_CASE( "P4PixelGrid reports target and predictor mask invalidation", "[P4PixelGrid][mask]" )
{
    mx::improc::P4PixelGridf reference;
    reference.resize( 41, 43, 20.0, 21.0 );
    reference.region( safeRegion( policyT::kernelSupport, 0.5, 0.0 ), nullptr );
    const std::size_t searchIndex = findSearch( reference, 20, 26 );
    const auto &target = reference.searchPixel( searchIndex ).coordinate();
    const auto &footprint = reference.interpolation( searchIndex, 0 );
    const int predictorMaskRow = footprint.footprintRow();
    const int predictorMaskColumn = footprint.footprintColumn();
    REQUIRE_FALSE( ( predictorMaskRow == target.row() && predictorMaskColumn == target.column() ) );

    gridT::imageT targetMask = gridT::imageT::Ones( 41, 43 );
    targetMask( target.row(), target.column() ) = 0;
    mx::improc::P4PixelGridf targetMaskedGrid;
    targetMaskedGrid.resize( 41, 43, 20.0, 21.0 );
    targetMaskedGrid.region( safeRegion( policyT::kernelSupport, 0.5, 0.0 ), &targetMask );
    REQUIRE( targetMaskedGrid.predictorCount() == reference.predictorCount() );
    REQUIRE( targetMaskedGrid.searchPixel( searchIndex ).invalidReason() == reasonT::targetMasked );
    REQUIRE_FALSE( targetMaskedGrid.searchPixel( searchIndex ).valid() );
    REQUIRE( targetMaskedGrid.searchPixel( searchIndex ).targetMasked() );
    REQUIRE_FALSE( targetMaskedGrid.searchPixel( searchIndex ).predictorMasked() );

    gridT::imageT predictorMask = gridT::imageT::Ones( 41, 43 );
    predictorMask( predictorMaskRow, predictorMaskColumn ) = 0;
    mx::improc::P4PixelGridf predictorMaskedGrid;
    predictorMaskedGrid.resize( 41, 43, 20.0, 21.0 );
    predictorMaskedGrid.region( safeRegion( policyT::kernelSupport, 0.5, 0.0 ), &predictorMask );
    REQUIRE( predictorMaskedGrid.predictorCount() == reference.predictorCount() );
    REQUIRE( predictorMaskedGrid.searchPixel( searchIndex ).invalidReason() == reasonT::predictorMasked );
    REQUIRE_FALSE( predictorMaskedGrid.searchPixel( searchIndex ).targetMasked() );
    REQUIRE( predictorMaskedGrid.searchPixel( searchIndex ).predictorMasked() );

    predictorMask( target.row(), target.column() ) = 0;
    mx::improc::P4PixelGridf combinedMaskedGrid;
    combinedMaskedGrid.resize( 41, 43, 20.0, 21.0 );
    combinedMaskedGrid.region( safeRegion( policyT::kernelSupport, 0.5, 0.0 ), &predictorMask );
    REQUIRE( combinedMaskedGrid.predictorCount() == reference.predictorCount() );
    REQUIRE( combinedMaskedGrid.searchPixel( searchIndex ).invalidReason() == reasonT::targetAndPredictorMasked );
    REQUIRE( combinedMaskedGrid.searchPixel( searchIndex ).targetMasked() );
    REQUIRE( combinedMaskedGrid.searchPixel( searchIndex ).predictorMasked() );

    gridT::imageT image = gridT::imageT::Ones( 41, 43 );
    REQUIRE_THROWS_AS( combinedMaskedGrid.sample( image, searchIndex, 0 ), std::logic_error );

    gridT::imageT wrongMask( 40, 43 );
    REQUIRE_THROWS_AS( reference.region( safeRegion(), &wrongMask ), std::invalid_argument );
    REQUIRE( reference.predictorCount() == combinedMaskedGrid.predictorCount() );

    gridT::imageT nonfiniteMask = gridT::imageT::Ones( 41, 43 );
    nonfiniteMask( 0, 0 ) = std::numeric_limits<float>::quiet_NaN();
    REQUIRE_THROWS_AS( reference.region( safeRegion(), &nonfiniteMask ), std::invalid_argument );
}

/** \brief Verifies P4PixelGrid rejects invalid, empty, edge-crossing, and all-excluded regions transactionally.
 *
 * Exercises validation in the real `P4PixelGrid::region()` API.
 * \ingroup P4PixelGrid_unit_tests
 */
TEST_CASE( "P4PixelGrid rejects invalid region configurations", "[P4PixelGrid][validation]" )
{
    mx::improc::P4PixelGridf emptyGrid;
    REQUIRE_THROWS_AS( emptyGrid.region( safeRegion(), nullptr ), std::logic_error );

    mx::improc::P4PixelGridf grid;
    grid.resize( 41, 43, 20.0, 21.0 );
    grid.region( safeRegion(), nullptr );
    const std::size_t oldSearchCount = grid.searchPixelCount();
    const std::size_t oldPredictorCount = grid.predictorCount();

    REQUIRE_THROWS_AS( grid.region( regionT( -1, 2, 1, 1, 1, 30, 0.5, policyT::sampleCenter, 0 ), nullptr ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( grid.region( regionT( 2, 2, 1, 1, 1, 30, 0.5, policyT::sampleCenter, 0 ), nullptr ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( grid.region( regionT( 2, 3, 0, 1, 1, 30, 0.5, policyT::sampleCenter, 0 ), nullptr ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( grid.region( regionT( 2, 3, 1, 0, 1, 30, 0.5, policyT::sampleCenter, 0 ), nullptr ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( grid.region( regionT( 2, 3, 1, 1, 0, 30, 0.5, policyT::sampleCenter, 0 ), nullptr ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( grid.region( regionT( 2, 3, 1, 1, 1, 0, 0.5, policyT::sampleCenter, 0 ), nullptr ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( grid.region( regionT( 2, 3, 1, 1, 1, 181, 0.5, policyT::sampleCenter, 0 ), nullptr ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( grid.region( regionT( 2, 3, 1, 1, 1, 180, 0.5, policyT::sampleCenter, 0 ), nullptr ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( grid.region( regionT( 2, 3, 1, 1, 1, 30, 0, policyT::sampleCenter, 0 ), nullptr ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( grid.region( regionT( 2, 3, 1, 1, 1, 30, 0.5, policyT::sampleCenter, -0.1 ), nullptr ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS(
        grid.region( regionT( 2, 3, 1, 1, std::numeric_limits<double>::infinity(), 30, 0.5, policyT::sampleCenter, 0 ),
                     nullptr ),
        std::invalid_argument );
    REQUIRE_THROWS_AS( grid.region( regionT( 2, 3, 1, 1, 1, 30, 0.5, static_cast<policyT>( 99 ), 0 ), nullptr ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( grid.region( regionT( 0.1, 0.2, 1, 1, 1, 30, 0.5, policyT::sampleCenter, 0 ), nullptr ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( grid.region( regionT( 5, 6, 2, 2, 3, 45, 100, policyT::sampleCenter, 0 ), nullptr ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( grid.region( regionT( 5, 6, 2, 2, 3, 45, 100, policyT::kernelSupport, 0 ), nullptr ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( grid.region( regionT( 5,
                                             6,
                                             2,
                                             2,
                                             3,
                                             45,
                                             std::numeric_limits<double>::max(),
                                             policyT::sampleCenter,
                                             std::numeric_limits<double>::max() ),
                                    nullptr ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS(
        grid.region( regionT( 1e308, 1.1e308, 1, 1e308, 3, 45, 0.5, policyT::sampleCenter, 0 ), nullptr ),
        std::invalid_argument );
    REQUIRE_THROWS_AS( grid.region( regionT( 5, 6, 2, 1e20, 3, 45, 0.5, policyT::sampleCenter, 0 ), nullptr ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( grid.region( regionT( 5, 6, 2, 100, 3, 45, 0.5, policyT::sampleCenter, 0 ), nullptr ),
                       std::invalid_argument );

    REQUIRE( mx::improc::P4PixelGridTestAccess::checkedFootprintOrigin( 4.75 ) == 3 );
    const double minimumMapped =
        static_cast<double>( std::numeric_limits<int>::min() ) + static_cast<double>( gridT::leftBuffer );
    const double maximumMapped =
        static_cast<double>( std::numeric_limits<int>::max() ) + static_cast<double>( gridT::leftBuffer );
    REQUIRE( mx::improc::P4PixelGridTestAccess::checkedFootprintOrigin( minimumMapped ) ==
             std::numeric_limits<int>::min() );
    REQUIRE( mx::improc::P4PixelGridTestAccess::checkedFootprintOrigin( maximumMapped ) ==
             std::numeric_limits<int>::max() );
    REQUIRE_THROWS_AS( mx::improc::P4PixelGridTestAccess::checkedFootprintOrigin( minimumMapped - 1.0 ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( mx::improc::P4PixelGridTestAccess::checkedFootprintOrigin( maximumMapped + 1.0 ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS(
        mx::improc::P4PixelGridTestAccess::checkedFootprintOrigin( std::numeric_limits<double>::quiet_NaN() ),
        std::invalid_argument );
    REQUIRE(
        mx::improc::P4PixelGridTestAccess::checkedInterpolationCount( 0, std::numeric_limits<std::size_t>::max() ) ==
        0 );
    REQUIRE( mx::improc::P4PixelGridTestAccess::checkedInterpolationCount( 7, 11 ) == 77 );
    REQUIRE_THROWS_AS(
        mx::improc::P4PixelGridTestAccess::checkedInterpolationCount( 2, std::numeric_limits<std::size_t>::max() ),
        std::length_error );

    REQUIRE( grid.searchPixelCount() == oldSearchCount );
    REQUIRE( grid.predictorCount() == oldPredictorCount );
    REQUIRE( grid.regionConfiguration().searchInnerRadius == 5.0 );

    mx::improc::P4PixelGridf small;
    small.resize( 15, 15 );
    REQUIRE_THROWS_AS( small.region( regionT( 4, 5, 1, 5, 3, 45, 0.5, policyT::sampleCenter, 0 ), nullptr ),
                       std::invalid_argument );

    mx::improc::P4PixelGridf exactEdge;
    exactEdge.resize( 9, 9 );
    exactEdge.region( regionT( 1, 2, 1, 1, 3, 45, 0.5, policyT::sampleCenter, 0 ), nullptr );
    bool touchesExactEdge{ false };
    for( std::size_t searchIndex = 0; searchIndex < exactEdge.searchPixelCount(); ++searchIndex )
    {
        for( std::size_t predictorIndex = 0; predictorIndex < exactEdge.predictorCount(); ++predictorIndex )
        {
            const auto &record = exactEdge.interpolation( searchIndex, predictorIndex );
            REQUIRE( record.footprintRow() >= 0 );
            REQUIRE( record.footprintColumn() >= 0 );
            REQUIRE( record.footprintRow() + gridT::width <= exactEdge.rows() );
            REQUIRE( record.footprintColumn() + gridT::width <= exactEdge.columns() );
            touchesExactEdge = touchesExactEdge || record.footprintRow() == 0 || record.footprintColumn() == 0 ||
                               record.footprintRow() + gridT::width == exactEdge.rows() ||
                               record.footprintColumn() + gridT::width == exactEdge.columns();
        }
    }
    REQUIRE( touchesExactEdge );
    const std::size_t exactEdgePredictorCount = exactEdge.predictorCount();
    REQUIRE_THROWS_AS( exactEdge.region( regionT( 1, 2, 2, 2, 3, 45, 0.5, policyT::sampleCenter, 0 ), nullptr ),
                       std::invalid_argument );
    REQUIRE( exactEdge.predictorCount() == exactEdgePredictorCount );

    mx::improc::P4PixelGridf central;
    central.resize( 9, 9, 4.0, 4.1 );
    REQUIRE_THROWS_AS( central.region( regionT( 0, 0.11, 0.1, 0.01, 0.1, 30, 0.1, policyT::sampleCenter, 0 ), nullptr ),
                       std::invalid_argument );

    mx::improc::P4PixelGridf underflow;
    underflow.resize( 9, 9, 4.0, 4.0 );
    REQUIRE_THROWS_AS(
        underflow.region(
            regionT( 0, std::numeric_limits<double>::denorm_min(), 0.1, 1, 0.1, 30, 0.1, policyT::sampleCenter, 0 ),
            nullptr ),
        std::invalid_argument );
}

/** \brief Verifies P4PixelGrid geometry and read-only sampling are deterministic across OpenMP worker counts.
 *
 * Exercises the real `P4PixelGrid::region()` and `P4PixelGrid::sample()` APIs.
 * \ingroup P4PixelGrid_unit_tests
 */
TEST_CASE( "P4PixelGrid is deterministic across OpenMP worker counts", "[P4PixelGrid][openmp]" )
{
#ifdef __DOXY_ONLY__
    // clang-format off
    mx::improc::P4PixelGrid<mx::improc::cubicConvolTransform<float>> doxygenGrid;
    mx::improc::P4PixelGridRegion doxygenRegion( 5, 6, 2, 2, 3, 45, 0.5,
                                                mx::improc::P4ExclusionPolicy::kernelSupport, 0 );
    doxygenGrid.resize( 41, 43, 20, 21 );
    doxygenGrid.region( doxygenRegion, nullptr );
    mx::improc::P4PixelGrid<mx::improc::cubicConvolTransform<float>>::imageT doxygenImage;
    doxygenGrid.sample( doxygenImage, 0, 0 );
    // clang-format on
#endif

    const OmpThreadCountGuard threadCountGuard;
    std::vector<gridT> grids( 4 );
    omp_set_num_threads( 1 );
    grids[0].resize( 41, 43, 20.0, 21.0 );
    grids[0].region( safeRegion( policyT::kernelSupport, 0.5, 0.0 ), nullptr );

    omp_set_num_threads( 3 );
#pragma omp parallel for
    for( int index = 1; index < static_cast<int>( grids.size() ); ++index )
    {
        grids[index].resize( 41, 43, 20.0, 21.0 );
        grids[index].region( safeRegion( policyT::kernelSupport, 0.5, 0.0 ), nullptr );
    }

    for( std::size_t index = 1; index < grids.size(); ++index )
    {
        requireSameGrid( grids[0], grids[index] );
    }

    gridT::imageT image( 41, 43 );
    for( int row = 0; row < image.rows(); ++row )
    {
        for( int column = 0; column < image.cols(); ++column )
        {
            image( row, column ) = static_cast<float>( row * image.cols() + column );
        }
    }

    const std::size_t sampleCount = grids[0].searchPixelCount() * grids[0].predictorCount();
    std::vector<float> serial( sampleCount );
    std::vector<float> parallel( sampleCount );
    omp_set_num_threads( 1 );
#pragma omp parallel for
    for( std::size_t index = 0; index < sampleCount; ++index )
    {
        serial[index] = grids[0].sample( image, index / grids[0].predictorCount(), index % grids[0].predictorCount() );
    }
    omp_set_num_threads( 3 );
#pragma omp parallel for
    for( std::size_t index = 0; index < sampleCount; ++index )
    {
        parallel[index] =
            grids[0].sample( image, index / grids[0].predictorCount(), index % grids[0].predictorCount() );
    }
    REQUIRE( serial == parallel );
}

} // namespace P4PixelGrid_test
} // namespace unitTest
