/** \file HCIobservation_coadd_test.cpp
 * \brief Tests HCIobservation image coaddition and metadata provenance.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "HCIobservation_test_fixture.hpp"

#include <limits>

namespace unitTest
{
namespace HCIobservation_coadd_test
{

/// \cond HCIobservation_test_harness
namespace
{
void makeCoaddInputs( HCIobservationTestHarness::cubeT &cube,
                      std::vector<std::string> &files,
                      std::vector<double> &dates,
                      std::vector<HCIobservationTestHarness::fitsHeaderT> &heads,
                      const std::vector<float> &values,
                      const std::vector<double> &seconds )
{
    cube.resize( 1, 1, values.size() );
    files.clear();
    dates.clear();
    heads.clear();
    files.resize( values.size() );
    dates.resize( values.size() );
    heads.resize( values.size() );

    constexpr double startMJD = 60000;
    for( size_t index = 0; index < values.size(); ++index )
    {
        cube.image( index )( 0, 0 ) = values[index];
        files[index] = "/inputs/frame" + std::to_string( index ) + ".fits";
        dates[index] = startMJD + seconds[index] / 86400.0;
        heads[index].append( "ANGLE", 10.0 * index, "test angle" );
        heads[index].append( "ORIGIN", static_cast<int>( index ), "source header index" );
    }
}
} // namespace
/// \endcond

/// Verify HCIobservation::coaddImages enforces image-count groups, combines remainders, and preserves group headers.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation count-limited mean coaddition", "[HCIobservation][coadd]" )
{
    HCIobservationTestHarness observation;
    HCIobservationTestHarness::cubeT cube;
    std::vector<std::string> files;
    std::vector<double> dates;
    std::vector<HCIobservationTestHarness::fitsHeaderT> heads;
    makeCoaddInputs( cube, files, dates, heads, { 1, 2, 3, 4, 5 }, { 0, 1, 2, 3, 4 } );

    observation.coaddImages( mx::improc::HCI::coadd::mean, 2, 0, { "ANGLE" }, files, "DATE-OBS", dates, heads, cube );

    REQUIRE( cube.planes() == 3 );
    REQUIRE( cube.image( 0 )( 0, 0 ) == Approx( 1.5 ) );
    REQUIRE( cube.image( 1 )( 0, 0 ) == Approx( 3.5 ) );
    REQUIRE( cube.image( 2 )( 0, 0 ) == Approx( 5 ) );
    REQUIRE( dates.size() == 3 );
    REQUIRE( ( dates[0] - 60000 ) * 86400 == Approx( 0.5 ).margin( 1e-6 ) );
    REQUIRE( heads.size() == 3 );
    REQUIRE( heads[0]["ORIGIN"].value<int>() == 0 );
    REQUIRE( heads[1]["ORIGIN"].value<int>() == 2 );
    REQUIRE( heads[2]["ORIGIN"].value<int>() == 4 );
    REQUIRE( heads[0]["ANGLE"].value<double>() == Approx( 5 ) );
    REQUIRE( heads[0]["START ANGLE"].value<double>() == Approx( 0 ) );
    REQUIRE( heads[0]["END ANGLE"].value<double>() == Approx( 10 ) );
    REQUIRE( heads[0]["DELTA ANGLE"].value<double>() == Approx( 10 ) );
    REQUIRE( heads[0]["IMAGES COADDED"].value<int>() == 2 );
    REQUIRE( heads[0].count( "HISTORY" ) == 2 );
    REQUIRE( heads[2]["START ANGLE"].value<double>() == Approx( 40 ) );
    REQUIRE( heads[2]["END ANGLE"].value<double>() == Approx( 40 ) );
    REQUIRE( heads[2]["DELTA ANGLE"].value<double>() == Approx( 0 ) );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::coaddImages();
#endif
    // clang-format on
}

/// Verify HCIobservation::coaddImages unwraps only the configured angular metadata before averaging.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation wrap-aware angle coaddition", "[HCIobservation][coadd][angle]" )
{
    const auto runCase =
        []( const std::vector<double> &angles, double angleScale, double expectedAngle, double expectedAngleDelta )
    {
        HCIobservationTestHarness observation;
        HCIobservationTestHarness::cubeT cube;
        std::vector<std::string> files;
        std::vector<double> dates;
        std::vector<HCIobservationTestHarness::fitsHeaderT> heads;
        std::vector<float> values;
        std::vector<double> seconds;
        values.reserve( angles.size() );
        seconds.reserve( angles.size() );
        for( size_t index = 0; index < angles.size(); ++index )
        {
            values.push_back( static_cast<float>( index + 1 ) );
            seconds.push_back( static_cast<double>( index ) );
        }
        makeCoaddInputs( cube, files, dates, heads, values, seconds );

        double scalarSum = 0;
        for( size_t index = 0; index < angles.size(); ++index )
        {
            heads[index]["ANGLE"].value( angles[index] );
            heads[index].append( "SCALAR", angles[index], "ordinary scalar metadata" );
            scalarSum += angles[index];
        }

        observation.m_coaddAngleKeyword = "ANGLE";
        observation.m_coaddAngleScale = static_cast<float>( angleScale );
        observation.coaddImages( mx::improc::HCI::coadd::mean,
                                 static_cast<int>( angles.size() ),
                                 0,
                                 { "ANGLE", "SCALAR" },
                                 files,
                                 "DATE-OBS",
                                 dates,
                                 heads,
                                 cube );

        REQUIRE( cube.planes() == 1 );
        REQUIRE( heads[0]["ANGLE"].value<double>() == Approx( expectedAngle ) );
        REQUIRE( heads[0]["START ANGLE"].value<double>() == Approx( angles.front() ) );
        REQUIRE( heads[0]["END ANGLE"].value<double>() == Approx( angles.back() ) );
        REQUIRE( heads[0]["DELTA ANGLE"].value<double>() == Approx( expectedAngleDelta ) );
        REQUIRE( heads[0]["SCALAR"].value<double>() == Approx( scalarSum / angles.size() ) );
        REQUIRE( heads[0]["DELTA SCALAR"].value<double>() == Approx( angles.back() - angles.front() ) );
    };

    runCase( { 358, 359, 1, 2 }, 1, 360, 4 );
    runCase( { 2, 1, 359, 358 }, 1, 0, -4 );
    runCase( { -270, -270, 89, 89 }, 1, -270.5, -1 );
    runCase( { 10, 20, 40 }, 1, 70.0 / 3.0, 30 );
    runCase( { 179, 179.5, 0.5, 1 }, 2, 180, 2 );
    runCase( { -179, -179.5, -0.5, -1 }, -2, -180, -2 );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::coaddImages();
#endif
    // clang-format on
}

/// Verify HCIobservation::coaddImages limits each group by its wrap-aware field-angle span.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation angle-limited coaddition", "[HCIobservation][coadd][angle][limit]" )
{
    HCIobservationTestHarness observation;
    HCIobservationTestHarness::cubeT cube;
    std::vector<std::string> files;
    std::vector<double> dates;
    std::vector<HCIobservationTestHarness::fitsHeaderT> heads;

    makeCoaddInputs( cube, files, dates, heads, { 1, 2, 3, 4 }, { 0, 1, 2, 3 } );
    heads[0]["ANGLE"].value( 358.0 );
    heads[1]["ANGLE"].value( 359.0 );
    heads[2]["ANGLE"].value( 0.0 );
    heads[3]["ANGLE"].value( 1.0 );
    observation.m_coaddAngleKeyword = "ANGLE";
    observation.m_coaddAngleScale = 1;
    observation.m_coaddMaxAngle = 2;

    observation.coaddImages( mx::improc::HCI::coadd::mean, 0, 0, { "ANGLE" }, files, "DATE-OBS", dates, heads, cube );

    REQUIRE( cube.planes() == 2 );
    REQUIRE( cube.image( 0 )( 0, 0 ) == Approx( 2 ) );
    REQUIRE( cube.image( 1 )( 0, 0 ) == Approx( 4 ) );
    REQUIRE( heads[0]["ANGLE"].value<double>() == Approx( 359 ) );
    REQUIRE( heads[0]["DELTA ANGLE"].value<double>() == Approx( 2 ) );
    REQUIRE( heads[0]["IMAGES COADDED"].value<int>() == 3 );
    REQUIRE( heads[1]["ANGLE"].value<double>() == Approx( 1 ) );
    REQUIRE( heads[1]["IMAGES COADDED"].value<int>() == 1 );

    makeCoaddInputs( cube, files, dates, heads, { 1, 2, 3 }, { 0, 0.2, 0.4 } );
    heads[0]["ANGLE"].value( 0.0 );
    heads[1]["ANGLE"].value( 3.6 );
    heads[2]["ANGLE"].value( 7.2 );
    observation.m_coaddMaxAngle = 5;
    observation.coaddImages( mx::improc::HCI::coadd::mean, 0, 0, { "ANGLE" }, files, "DATE-OBS", dates, heads, cube );
    REQUIRE( cube.planes() == 2 );
    REQUIRE( heads[0]["IMAGES COADDED"].value<int>() == 2 );
    REQUIRE( heads[0]["DELTA ANGLE"].value<double>() == Approx( 3.6 ) );

    makeCoaddInputs( cube, files, dates, heads, { 1, 2, 3 }, { 0, 1, 2 } );
    heads[0]["ANGLE"].value( 0.0 );
    heads[1]["ANGLE"].value( 4.0 );
    heads[2]["ANGLE"].value( -4.0 );
    observation.m_coaddMaxAngle = 5;
    observation.coaddImages( mx::improc::HCI::coadd::mean, 0, 0, { "ANGLE" }, files, "DATE-OBS", dates, heads, cube );
    REQUIRE( cube.planes() == 2 );
    REQUIRE( heads[0]["IMAGES COADDED"].value<int>() == 2 );
    REQUIRE( heads[1]["IMAGES COADDED"].value<int>() == 1 );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::coaddImages();
#endif
    // clang-format on
}

/// Verify HCIobservation::coaddImages includes exact time boundaries and separates larger gaps.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation time-limited coaddition", "[HCIobservation][coadd][time]" )
{
    HCIobservationTestHarness observation;
    HCIobservationTestHarness::cubeT cube;
    std::vector<std::string> files;
    std::vector<double> dates;
    std::vector<HCIobservationTestHarness::fitsHeaderT> heads;

    makeCoaddInputs( cube, files, dates, heads, { 1, 3 }, { 0, 2 } );
    observation.coaddImages( mx::improc::HCI::coadd::mean, 0, 2, {}, files, "DATE-OBS", dates, heads, cube );
    REQUIRE( cube.planes() == 1 );
    REQUIRE( cube.image( 0 )( 0, 0 ) == Approx( 2 ) );

    makeCoaddInputs( cube, files, dates, heads, { 1, 3 }, { 0, 2.01 } );
    observation.coaddImages( mx::improc::HCI::coadd::mean, 0, 2, {}, files, "DATE-OBS", dates, heads, cube );
    REQUIRE( cube.planes() == 2 );

    makeCoaddInputs( cube, files, dates, heads, { 1, 2, 10, 11, 12 }, { 0, 1, 10, 11, 12 } );
    observation.coaddImages( mx::improc::HCI::coadd::mean, 3, 2, {}, files, "DATE-OBS", dates, heads, cube );
    REQUIRE( cube.planes() == 2 );
    REQUIRE( cube.image( 0 )( 0, 0 ) == Approx( 1.5 ) );
    REQUIRE( cube.image( 1 )( 0, 0 ) == Approx( 11 ) );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::coaddImages();
#endif
    // clang-format on
}

/// Verify HCIobservation::coaddImages supports count-only median coadds without dates and uses supplied RDI provenance.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation date-free median and RDI coaddition", "[HCIobservation][coadd][RDI]" )
{
    HCIobservationTestHarness observation;
    HCIobservationTestHarness::cubeT cube;
    std::vector<std::string> files{ "/rdi/ref-a.fits", "/rdi/ref-b.fits", "/rdi/ref-c.fits", "/rdi/ref-d.fits" };
    std::vector<double> dates;
    std::vector<HCIobservationTestHarness::fitsHeaderT> heads( 4 );
    cube.resize( 1, 1, 4 );
    cube.image( 0 )( 0, 0 ) = 1;
    cube.image( 1 )( 0, 0 ) = 4;
    cube.image( 2 )( 0, 0 ) = 7;
    cube.image( 3 )( 0, 0 ) = 10;

    observation.coaddImages( mx::improc::HCI::coadd::median, 4, 0, {}, files, "", dates, heads, cube );
    REQUIRE( cube.planes() == 1 );
    REQUIRE( cube.image( 0 )( 0, 0 ) == Approx( 5.5 ) );
    REQUIRE( dates.empty() );
    REQUIRE( heads[0].count( "HISTORY" ) == 4 );
    REQUIRE( heads[0]["IMAGES COADDED"].value<int>() == 4 );

    std::vector<std::string> history;
    for( auto &card : heads[0] )
    {
        if( card.keyword() == "HISTORY" )
        {
            history.push_back( card.comment() );
        }
    }
    REQUIRE( history == std::vector<std::string>{ "coadded file: ref-a.fits",
                                                  "coadded file: ref-b.fits",
                                                  "coadded file: ref-c.fits",
                                                  "coadded file: ref-d.fits" } );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::coaddImages();
#endif
    // clang-format on
}

/// Verify HCIobservation::coaddImages no-op and invalid metadata contracts.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation coadd validation", "[HCIobservation][coadd]" )
{
    HCIobservationTestHarness observation;
    HCIobservationTestHarness::cubeT cube;
    std::vector<std::string> files;
    std::vector<double> dates;
    std::vector<HCIobservationTestHarness::fitsHeaderT> heads;
    makeCoaddInputs( cube, files, dates, heads, { 1, 2 }, { 0, 1 } );

    observation.coaddImages( mx::improc::HCI::coadd::none, 1, 0, {}, files, "DATE-OBS", dates, heads, cube );
    REQUIRE( cube.planes() == 2 );
    observation.coaddImages( mx::improc::HCI::coadd::mean, 0, 0, {}, files, "DATE-OBS", dates, heads, cube );
    REQUIRE( cube.planes() == 2 );

    dates.clear();
    REQUIRE_THROWS(
        observation.coaddImages( mx::improc::HCI::coadd::mean, 0, 1, {}, files, "DATE-OBS", dates, heads, cube ) );

    dates = { 60000 };
    REQUIRE_THROWS(
        observation.coaddImages( mx::improc::HCI::coadd::mean, 1, 0, {}, files, "DATE-OBS", dates, heads, cube ) );

    dates = { 60000, 60000.1 };
    files.pop_back();
    REQUIRE_THROWS(
        observation.coaddImages( mx::improc::HCI::coadd::mean, 1, 0, {}, files, "DATE-OBS", dates, heads, cube ) );

    files.push_back( "/inputs/frame1.fits" );
    heads.pop_back();
    REQUIRE_THROWS(
        observation.coaddImages( mx::improc::HCI::coadd::mean, 1, 0, {}, files, "DATE-OBS", dates, heads, cube ) );

    heads.resize( 2 );
    REQUIRE_THROWS( observation.coaddImages( static_cast<mx::improc::HCI::coadd>( 99 ),
                                             1,
                                             0,
                                             {},
                                             files,
                                             "DATE-OBS",
                                             dates,
                                             heads,
                                             cube ) );

    HCIobservationTestHarness::cubeT emptyCube;
    std::vector<std::string> emptyFiles;
    std::vector<double> emptyDates;
    std::vector<HCIobservationTestHarness::fitsHeaderT> emptyHeads;
    REQUIRE_THROWS( observation.coaddImages( mx::improc::HCI::coadd::mean,
                                             1,
                                             0,
                                             {},
                                             emptyFiles,
                                             "DATE-OBS",
                                             emptyDates,
                                             emptyHeads,
                                             emptyCube ) );

    makeCoaddInputs( cube, files, dates, heads, { 1, 2 }, { 0, 1 } );
    dates[1] = std::numeric_limits<double>::quiet_NaN();
    REQUIRE_THROWS(
        observation.coaddImages( mx::improc::HCI::coadd::mean, 1, 0, {}, files, "DATE-OBS", dates, heads, cube ) );

    dates = { 60001, 60000 };
    REQUIRE_THROWS(
        observation.coaddImages( mx::improc::HCI::coadd::mean, 1, 0, {}, files, "DATE-OBS", dates, heads, cube ) );

    makeCoaddInputs( cube, files, dates, heads, { 1, 2 }, { 0, 1 } );
    observation.m_coaddMaxAngle = std::numeric_limits<float>::quiet_NaN();
    REQUIRE_THROWS(
        observation.coaddImages( mx::improc::HCI::coadd::mean, 1, 0, {}, files, "DATE-OBS", dates, heads, cube ) );

    observation.m_coaddMaxAngle = 1;
    REQUIRE_THROWS(
        observation.coaddImages( mx::improc::HCI::coadd::mean, 1, 0, {}, files, "DATE-OBS", dates, heads, cube ) );

    observation.m_coaddAngleKeyword = "ANGLE";
    observation.m_coaddAngleScale = std::numeric_limits<float>::quiet_NaN();
    REQUIRE_THROWS(
        observation
            .coaddImages( mx::improc::HCI::coadd::mean, 1, 0, { "ANGLE" }, files, "DATE-OBS", dates, heads, cube ) );

    observation.m_coaddAngleScale = 1;
    REQUIRE_THROWS(
        observation.coaddImages( mx::improc::HCI::coadd::mean, 1, 0, {}, files, "DATE-OBS", dates, heads, cube ) );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::coaddImages();
#endif
    // clang-format on
}

} // namespace HCIobservation_coadd_test
} // namespace unitTest
