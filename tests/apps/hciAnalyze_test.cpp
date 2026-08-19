/** \file hciAnalyze_test.cpp
 * \brief Tests hciAnalyze configuration, position resolution, filtering, and SNR measurement.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "../common/HCIobservation_test_fixture.hpp"

#include <cmath>
#include <limits>

#define HCIREDUCE_HCIANALYZE_NO_MAIN
#include "src/apps/hciAnalyze.cpp"

namespace unitTest
{
namespace hciAnalyze_test
{

/// \cond hciAnalyze_test_harness
struct appHarness : public hciAnalyze
{
    using hciAnalyze::analyzeCube;
    using hciAnalyze::checkConfig;
    using hciAnalyze::config;
    using hciAnalyze::filterCube;
    using hciAnalyze::loadConfig;
    using hciAnalyze::m_fakeContrast;
    using hciAnalyze::m_fakePositionAngle;
    using hciAnalyze::m_fakeSeparation;
    using hciAnalyze::m_fakeSpecified;
    using hciAnalyze::m_file;
    using hciAnalyze::m_planetPositionAngle;
    using hciAnalyze::m_planetRadius;
    using hciAnalyze::m_planetSeparation;
    using hciAnalyze::m_planetSpecified;
    using hciAnalyze::m_results;
    using hciAnalyze::m_signals;
    using hciAnalyze::m_snrApertureRadius;
    using hciAnalyze::m_snrMaxRadius;
    using hciAnalyze::m_snrMinRadius;
    using hciAnalyze::resolveSignals;
    using hciAnalyze::setupConfig;
    using hciAnalyze::snrAnnulus;
};
/// \endcond

/// Verify hciAnalyze maps the public east-of-north PA convention into image coordinates.
/** \ingroup hciAnalyze_unit_tests */
TEST_CASE( "hciAnalyze PA coordinate convention", "[hciAnalyze][coordinates]" )
{
    const auto [northX, northY] = hciAnalyze::signalCoordinates( 2, 0, 11, 11 );
    REQUIRE( northX == Approx( 5 ) );
    REQUIRE( northY == Approx( 7 ) );

    const auto [eastX, eastY] = hciAnalyze::signalCoordinates( 2, 90, 11, 11 );
    REQUIRE( eastX == Approx( 3 ) );
    REQUIRE( eastY == Approx( 5 ) );
}

/// Verify hciAnalyze validates explicit planet and fake coordinate configuration before input processing.
/** \ingroup hciAnalyze_unit_tests */
TEST_CASE( "hciAnalyze configuration validation", "[hciAnalyze][config][validation]" )
{
    appHarness application;
    application.m_file = "input.fits";
    application.m_planetSeparation = { 3 };
    application.m_planetPositionAngle = { 0, 90 };
    application.m_planetSpecified = true;
    REQUIRE_THROWS( application.checkConfig() );

    application.m_planetPositionAngle = { 0 };
    application.m_fakeSeparation = { 3 };
    application.m_fakePositionAngle = { 0 };
    application.m_fakeSpecified = true;
    REQUIRE_THROWS( application.checkConfig() );
}

/// Verify hciAnalyze reads a normal configuration file and measures its configured FITS cube.
/** \ingroup hciAnalyze_unit_tests */
TEST_CASE( "hciAnalyze config-file SNR dispatch", "[hciAnalyze][config][execute]" )
{
    TestDirectory directory;
    hciAnalyze::cubeT cube( 15, 15, 1 );
    for( int column = 0; column < cube.cols(); ++column )
    {
        for( int row = 0; row < cube.rows(); ++row )
        {
            cube.image( 0 )( row, column ) = static_cast<float>( ( row + column ) % 5 );
        }
    }
    cube.image( 0 )( 5, 7 ) = 20;
    writeFitsCube( directory.file( "product.fits" ), cube );

    const auto configPath = directory.file( "hciAnalyze.conf" );
    writeTextFile( configPath,
                   "file=" + directory.file( "product.fits" ).string() +
                       "\n"
                       "[planet]\n"
                       "sep=2\n"
                       "PA=0\n"
                       "[snr]\n"
                       "minRad=1\n"
                       "maxRad=6\n"
                       "apertureR=1\n" );

    appHarness application;
    std::string invokedName = "hciAnalyze-test";
    std::string configOption = "--config";
    std::string configName = configPath.string();
    char *arguments[]{ invokedName.data(), configOption.data(), configName.data() };
    REQUIRE( application.main( 3, arguments ) == 0 );
    REQUIRE( application.m_results.size() == 1 );
    REQUIRE( std::isfinite( application.m_results[0].m_snr ) );
    REQUIRE( std::filesystem::exists( directory.file( "product_snr.fits" ) ) );

    hciAnalyze::cubeT snrCube;
    hciAnalyze::fitsHeaderT snrHeader;
    mx::fits::fitsFile<float, mx::verbose::vv> reader;
    REQUIRE( reader.read( snrCube, snrHeader, directory.file( "product_snr.fits" ).string() ) == mx::error_t::noerror );
    REQUIRE( snrHeader["HCIANAL"].String() == "hciAnalyze" );
    REQUIRE( snrHeader["SNRMINR"].value<float>() == Approx( 1 ) );
    REQUIRE( snrHeader["SNRMAXR"].value<float>() == Approx( 6 ) );
    REQUIRE( snrHeader["SNRAPER"].value<float>() == Approx( 1 ) );
    REQUIRE( snrHeader["HPFGFW"].value<float>() == Approx( 0 ) );
    REQUIRE( snrHeader["LPFGFW"].value<float>() == Approx( 0 ) );
    REQUIRE( snrHeader["HCISEPS"].String().starts_with( "2" ) );
    REQUIRE( snrHeader["HCIPAS"].String().starts_with( "0" ) );
    REQUIRE( snrHeader["HCIRADS"].String().starts_with( "10" ) );

    appHarness overwriteApplication;
    REQUIRE( overwriteApplication.main( 3, arguments ) == 0 );
    REQUIRE( std::filesystem::exists( directory.file( "product_snr.fits" ) ) );
}

/// Verify hciAnalyze parses FITS fake-planet keywords and lets explicit planet coordinates take precedence.
/** \ingroup hciAnalyze_unit_tests */
TEST_CASE( "hciAnalyze fake keyword resolution", "[hciAnalyze][header][fake]" )
{
    appHarness application;
    hciAnalyze::fitsHeaderT header;
    header.append<std::string>( "FAKESEP", "2", "fake separation" );
    header.append<std::string>( "FAKEPA", "90", "fake position angle" );
    header.append<std::string>( "FAKECONT", "1e-4", "fake contrast" );

    application.resolveSignals( header, 11, 11 );
    REQUIRE( application.m_signals.size() == 1 );
    REQUIRE( application.m_signals[0].m_x == Approx( 3 ) );
    REQUIRE( application.m_signals[0].m_contrast == Approx( 1e-4F ) );

    application.m_planetSeparation = { 1 };
    application.m_planetPositionAngle = { 0 };
    application.m_planetSpecified = true;
    application.resolveSignals( header, 11, 11 );
    REQUIRE( application.m_signals[0].m_y == Approx( 6 ) );
    REQUIRE( std::isnan( application.m_signals[0].m_contrast ) );
}

/// Verify hciAnalyze rejects missing header annulus metadata and accepts explicit annulus bounds.
/** \ingroup hciAnalyze_unit_tests */
TEST_CASE( "hciAnalyze SNR annulus validation", "[hciAnalyze][header][annulus]" )
{
    appHarness application;
    hciAnalyze::fitsHeaderT header;
    REQUIRE_THROWS( application.snrAnnulus( header ) );

    header.append<std::string>( "P4MINR", "2,4", "P4 region minima" );
    header.append<std::string>( "P4MAXR", "5,7", "P4 region maxima" );
    const auto [p4Minimum, p4Maximum] = application.snrAnnulus( header );
    REQUIRE( p4Minimum == Approx( 2 ) );
    REQUIRE( p4Maximum == Approx( 7 ) );

    application.m_snrMinRadius = 1;
    application.m_snrMaxRadius = 4;
    const auto [minimum, maximum] = application.snrAnnulus( header );
    REQUIRE( minimum == Approx( 1 ) );
    REQUIRE( maximum == Approx( 4 ) );
}

/// Verify hciAnalyze applies Gaussian filters to every cube plane, including partial edge kernels.
/** \ingroup hciAnalyze_unit_tests */
TEST_CASE( "hciAnalyze Gaussian filtering", "[hciAnalyze][filter]" )
{
    appHarness application;
    hciAnalyze::cubeT cube( 15, 15, 1 );
    hciAnalyze::cubeT invalidMask( 15, 15, 1 );
    invalidMask.setZero();
    cube.setZero();
    cube.image( 0 )( 7, 7 ) = 1;

    application.filterCube( cube, invalidMask, 2, 0 );
    REQUIRE( cube.image( 0 )( 7, 7 ) < 1 );
    REQUIRE( cube.image( 0 )( 7, 7 ) > 0 );

    cube.setZero();
    cube.image( 0 )( 7, 7 ) = 1;
    application.filterCube( cube, invalidMask, 0, 2 );
    REQUIRE( cube.image( 0 )( 7, 7 ) < 1 );
    REQUIRE( cube.image( 0 )( 7, 7 ) > 0 );

    cube.image( 0 ).setConstant( 1 );
    application.filterCube( cube, invalidMask, 0, 2 );
    REQUIRE( cube.image( 0 )( 0, 0 ) == Approx( 1 ) );
    REQUIRE( cube.image( 0 )( 0, 7 ) == Approx( 1 ) );
    REQUIRE( cube.image( 0 )( 7, 0 ) == Approx( 1 ) );
    REQUIRE( cube.image( 0 )( 14, 14 ) == Approx( 1 ) );

    cube.image( 0 ).setConstant( 1 );
    cube.image( 0 )( 7, 7 ) = 0;
    invalidMask.image( 0 )( 7, 7 ) = 1;
    application.filterCube( cube, invalidMask, 0, 2 );
    REQUIRE( cube.image( 0 )( 7, 6 ) == Approx( 1 ) );
    REQUIRE( cube.image( 0 )( 7, 7 ) == 0 );
}

/// Verify hciAnalyze reports the maximum inside snr.apertureR while excluding invalid input pixels.
/** \ingroup hciAnalyze_unit_tests */
TEST_CASE( "hciAnalyze cube SNR measurement", "[hciAnalyze][snr][nan]" )
{
    appHarness application;
    application.m_snrMinRadius = 1;
    application.m_snrMaxRadius = 6;
    application.m_snrApertureRadius = 1;
    application.m_planetSeparation = { 5 };
    application.m_planetPositionAngle = { 0 };
    application.m_planetRadius = { 3 };
    application.m_planetSpecified = true;

    hciAnalyze::fitsHeaderT header;
    hciAnalyze::cubeT cube( 31, 31, 2 );
    for( int plane = 0; plane < cube.planes(); ++plane )
    {
        for( int column = 0; column < cube.cols(); ++column )
        {
            for( int row = 0; row < cube.rows(); ++row )
            {
                cube.image( plane )( row, column ) = static_cast<float>( ( row + 2 * column + plane ) % 7 );
            }
        }
    }
    cube.image( 0 )( 0, 0 ) = std::numeric_limits<float>::quiet_NaN();
    cube.image( 0 )( 15, 10 ) = 20;
    cube.image( 0 )( 15, 13 ) = 200;

    application.resolveSignals( header, cube.rows(), cube.cols() );
    application.analyzeCube( cube, header );

    REQUIRE( application.m_results.size() == 2 );
    REQUIRE( std::isfinite( application.m_results[0].m_snr ) );
    REQUIRE( application.m_results[0].m_snr < 20 );
    REQUIRE( std::isfinite( application.m_results[1].m_snr ) );

    // clang-format off
#ifdef __DOXY_ONLY__
    hciAnalyze doxygenApplication;
    doxygenApplication.signalCoordinates( 0, 0, 1, 1 );
#endif
    // clang-format on
}

} // namespace hciAnalyze_test
} // namespace unitTest
