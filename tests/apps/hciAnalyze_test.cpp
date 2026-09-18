/** \file hciAnalyze_test.cpp
 * \brief Tests hciAnalyze configuration, position resolution, filtering, and SNR measurement.
 */

#include "../catch2/catch.hpp"

#include "../common/HCIobservation_test_fixture.hpp"

#include <cmath>
#include <limits>
#include <numbers>
#include <sstream>

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
    using hciAnalyze::correctSmallSampleSNR;
    using hciAnalyze::execute;
    using hciAnalyze::filterCube;
    using hciAnalyze::filterCubePSFResponse;
    using hciAnalyze::loadConfig;
    using hciAnalyze::m_diagnostics;
    using hciAnalyze::m_fakeContrast;
    using hciAnalyze::m_fakePositionAngle;
    using hciAnalyze::m_fakeSeparation;
    using hciAnalyze::m_fakeSpecified;
    using hciAnalyze::m_file;
    using hciAnalyze::m_highPassFwhm;
    using hciAnalyze::m_lambdaD;
    using hciAnalyze::m_lambdaDSpecified;
    using hciAnalyze::m_noiseConfig;
    using hciAnalyze::m_noiseDiagnostics;
    using hciAnalyze::m_noiseExcludeColumns;
    using hciAnalyze::m_noiseExcludeRadii;
    using hciAnalyze::m_noiseExcludeRows;
    using hciAnalyze::m_noiseModel;
    using hciAnalyze::m_noiseOnly;
    using hciAnalyze::m_planetContrast;
    using hciAnalyze::m_planetPositionAngle;
    using hciAnalyze::m_planetRadius;
    using hciAnalyze::m_planetSeparation;
    using hciAnalyze::m_planetSpecified;
    using hciAnalyze::m_psfResponse;
    using hciAnalyze::m_psfResponseMinimumSupport;
    using hciAnalyze::m_results;
    using hciAnalyze::m_signals;
    using hciAnalyze::m_snrApertureRadius;
    using hciAnalyze::m_snrMaxRadius;
    using hciAnalyze::m_snrMinRadius;
    using hciAnalyze::resolveSignals;
    using hciAnalyze::setupConfig;
    using hciAnalyze::snrAnnulus;
};

/// Capture standard-error output and restore its stream buffer at scope exit.
class ErrorCapture
{
  public:
    /// Redirect standard error into the owned string buffer.
    ErrorCapture() : m_buffer(), m_original( std::cerr.rdbuf( m_buffer.rdbuf() ) )
    {
    }

    /// Restore the original standard-error stream buffer.
    ~ErrorCapture()
    {
        std::cerr.rdbuf( m_original );
    }

    /// Disallow copying a live standard-error capture.
    ErrorCapture( const ErrorCapture & ) = delete;

    /// Disallow assigning a live standard-error capture.
    ErrorCapture &operator=( const ErrorCapture & ) = delete;

    /// Return all text captured so far.
    std::string str() const
    {
        return m_buffer.str();
    }

  private:
    std::ostringstream m_buffer; ///< Captured standard-error text.
    std::streambuf *m_original;  ///< Original standard-error buffer restored by the destructor.
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
    application.setupConfig();
    REQUIRE( application.config.m_targets.count( "filter.psfResponse" ) == 1 );
    REQUIRE( application.config.m_targets.count( "lambdaD" ) == 1 );
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

    application.m_fakeSpecified = false;
    application.m_lambdaD = 0;
    REQUIRE_THROWS( application.checkConfig() );
}

/// Verify hciAnalyze honors an attached false value for its diagnostics boolean option.
/** This exercises hciAnalyze::loadConfig() through the real command-line parser.
 * \ingroup hciAnalyze_unit_tests
 */
TEST_CASE( "hciAnalyze explicit false command-line option", "[hciAnalyze][config][bool]" )
{
    appHarness application;
    application.setupConfig();

    std::string invoked = "hciAnalyze";
    std::string diagnostics = "--diagnostics=false";
    char *arguments[]{ invoked.data(), diagnostics.data() };
    application.config.parseCommandLine( 2, arguments );
    REQUIRE_NOTHROW( application.loadConfig() );
    REQUIRE_FALSE( application.m_diagnostics );

    // clang-format off
#ifdef __DOXY_ONLY__
    hciAnalyze doxygenApplication;
    doxygenApplication.loadConfig();
#endif
    // clang-format on
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
                       "lambdaD=2.5" +
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
    std::string explicitWarnings;
    {
        ErrorCapture errors;
        REQUIRE( application.main( 3, arguments ) == 0 );
        explicitWarnings = errors.str();
    }
    REQUIRE( explicitWarnings.find( "lambdaD was not set" ) == std::string::npos );
    REQUIRE( application.m_lambdaDSpecified );
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
    REQUIRE( snrHeader["LAMBDAD"].value<float>() == Approx( 2.5 ) );
    REQUIRE( snrHeader["SNRMEAN"].value<int>() == 1 );
    REQUIRE( snrHeader["SNRSMALL"].value<int>() == 1 );
    REQUIRE( snrHeader["HPFGFW"].value<float>() == Approx( 0 ) );
    REQUIRE( snrHeader["LPFGFW"].value<float>() == Approx( 0 ) );
    REQUIRE( snrHeader["HCISEPS"].String().starts_with( "2" ) );
    REQUIRE( snrHeader["HCIPAS"].String().starts_with( "0" ) );
    REQUIRE( snrHeader["HCIRADS"].String().starts_with( "10" ) );

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
    appHarness overwriteApplication;
    std::string defaultWarnings;
    {
        ErrorCapture errors;
        REQUIRE( overwriteApplication.main( 3, arguments ) == 0 );
        defaultWarnings = errors.str();
    }
    REQUIRE( defaultWarnings.find( "lambdaD was not set; using 2.5" ) != std::string::npos );
    REQUIRE_FALSE( overwriteApplication.m_lambdaDSpecified );
    REQUIRE( std::filesystem::exists( directory.file( "product_snr.fits" ) ) );
}

/// Verify hciAnalyze applies the requested Mawet small-sample correction as a function of lambda/D radius.
/** \ingroup hciAnalyze_unit_tests */
TEST_CASE( "hciAnalyze small-sample SNR correction", "[hciAnalyze][snr][small-sample]" )
{
    appHarness application;
    application.m_lambdaD = 2.5F;
    hciAnalyze::cubeT snrCube( 11, 11, 2 );
    snrCube.cube().setOnes();

    application.correctSmallSampleSNR( snrCube );

    const double comparisonSamples = 4.0 * std::numbers::pi - 1.0;
    const double expected = 1.0 / std::sqrt( 1.0 + 1.0 / comparisonSamples );
    REQUIRE( snrCube.image( 0 )( 5, 10 ) == Approx( expected ) );
    REQUIRE( snrCube.image( 1 )( 5, 10 ) == Approx( expected ) );
    REQUIRE( snrCube.image( 0 )( 5, 5 ) == 0 );
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

    application.m_planetSpecified = false;
    header.append<std::string>( "PLANETSEP", "1,2", "planet separations" );
    header.append<std::string>( "PLANETPA", "0,180", "planet position angles" );
    header.append<std::string>( "PLANETCONT", "0.001,0.002", "planet contrasts" );
    application.resolveSignals( header, 11, 11 );
    REQUIRE( application.m_signals.size() == 2 );
    REQUIRE( application.m_signals[0].m_y == Approx( 6 ) );
    REQUIRE( application.m_signals[0].m_contrast == Approx( 0.001F ) );
    REQUIRE( application.m_signals[1].m_y == Approx( 3 ) );
    REQUIRE( application.m_signals[1].m_contrast == Approx( 0.002F ) );
}

/// Verify hciAnalyze rejects missing header annulus metadata and accepts explicit annulus bounds.
/** \ingroup hciAnalyze_unit_tests */
TEST_CASE( "hciAnalyze SNR annulus validation", "[hciAnalyze][header][annulus]" )
{
    appHarness application;
    hciAnalyze::fitsHeaderT header;
    REQUIRE_THROWS( application.snrAnnulus( header ) );

    header.append<std::string>( "P4 MIN RADIUS", "2,4", "P4 region minima" );
    header.append<std::string>( "P4 MAX RADIUS", "5,7", "P4 region maxima" );
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

/// Verify hciAnalyze loads a complete P4 response field and applies its normalized per-pixel filters.
/** \ingroup hciAnalyze_unit_tests */
TEST_CASE( "hciAnalyze external per-pixel P4 PSF filtering", "[hciAnalyze][filter][p4PSF]" )
{
    TestDirectory directory;
    const std::filesystem::path manifestPath = directory.file( "field_manifest.fits" );

    hciAnalyze::fitsHeaderT manifestHeader;
    REQUIRE( manifestHeader.append<int>( "P4 PSF PRODUCT SCHEMA", 2, "schema" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<std::string>( "P4 PSF PRODUCT", "MANIFEST", "role" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<int>( "P4 PSF COMPLETE", 1, "complete" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<int>( "P4 PSF MODE COUNT", 1, "modes" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<std::string>( "P4 PSF SOURCE COUNT", "1", "sources" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<int>( "P4 PSF STAMP SIZE", 3, "stamp" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<float>( "P4 PSF FILTER MIN GOOD FRACTION", 1, "support" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<std::string>( "P4 MODE FRACTIONS", "0.5", "modes" ) == mx::error_t::noerror );
    hciAnalyze::imageT manifest( 1, 1 );
    manifest( 0, 0 ) = 1;
    writeFitsImage( manifestPath, manifest, &manifestHeader );

    hciAnalyze::fitsHeaderT coordinateHeader;
    REQUIRE( coordinateHeader.append<std::string>( "P4 PSF PRODUCT", "COORDINATES", "role" ) == mx::error_t::noerror );
    hciAnalyze::imageT coordinates( 1, 4 );
    coordinates << 4, 4, 0, 0;
    writeFitsImage( directory.file( "field_coordinates.fits" ), coordinates, &coordinateHeader );

    hciAnalyze::fitsHeaderT modelHeader;
    REQUIRE( modelHeader.append<std::string>( "P4 PSF PRODUCT", "MODEL", "role" ) == mx::error_t::noerror );
    REQUIRE( modelHeader.append<int>( "P4 PSF MODE INDEX", 0, "mode" ) == mx::error_t::noerror );
    hciAnalyze::cubeT models( 3, 3, 1 );
    models.setZero();
    models.image( 0 )( 1, 1 ) = 2;
    writeFitsCube( directory.file( "field_model_0000.fits" ), models, &modelHeader );

    hciAnalyze::fitsHeaderT validityHeader;
    REQUIRE( validityHeader.append<std::string>( "P4 PSF PRODUCT", "VALIDITY", "role" ) == mx::error_t::noerror );
    REQUIRE( validityHeader.append<int>( "P4 PSF MODE INDEX", 0, "mode" ) == mx::error_t::noerror );
    hciAnalyze::imageT validity( 1, 1 );
    validity( 0, 0 ) = 1;
    writeFitsImage( directory.file( "field_validity_0000.fits" ), validity, &validityHeader );

    appHarness application;
    application.m_psfResponse = manifestPath.string();
    hciAnalyze::cubeT science( 9, 9, 1 );
    science.setZero();
    science.image( 0 )( 4, 4 ) = 6;
    hciAnalyze::fitsHeaderT scienceHeader;
    REQUIRE( scienceHeader.append<std::string>( "P4 MODE FRACTIONS", "0.5", "modes" ) == mx::error_t::noerror );

    application.filterCubePSFResponse( science, scienceHeader );
    REQUIRE( science.image( 0 )( 4, 4 ) == Approx( 3 ) );
    REQUIRE( std::isnan( science.image( 0 )( 0, 0 ) ) );
    REQUIRE( application.m_psfResponseMinimumSupport == Approx( 1 ) );

    hciAnalyze::fitsHeaderT mismatchedHeader;
    REQUIRE( mismatchedHeader.append<std::string>( "P4 MODE FRACTIONS", "0.6", "modes" ) == mx::error_t::noerror );
    hciAnalyze::cubeT mismatchedScience( 9, 9, 1 );
    mismatchedScience.setZero();
    REQUIRE_THROWS( application.filterCubePSFResponse( mismatchedScience, mismatchedHeader ) );

    manifest( 0, 0 ) = 0;
    writeFitsImage( manifestPath, manifest, &manifestHeader );
    hciAnalyze::cubeT rejectedScience( 9, 9, 1 );
    rejectedScience.setZero();
    REQUIRE_THROWS( application.filterCubePSFResponse( rejectedScience, scienceHeader ) );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::P4PSFFilter::calculate( {}, {}, {}, 0, 0, 1 );
#endif
    // clang-format on
}

/// Verify hciAnalyze reconstructs sparse radial KLIP responses and applies them at every eligible science pixel.
/** This exercises hciAnalyze::filterCubePSFResponse(), mx::improc::RadialPSFModel::fit(),
 * mx::improc::RadialPSFModel::response(), and mx::improc::P4PSFFilter::calculate().
 * \ingroup hciAnalyze_unit_tests
 */
TEST_CASE( "hciAnalyze external sparse radial KLIP PSF filtering", "[hciAnalyze][filter][klipPSF]" )
{
    TestDirectory directory;
    const std::filesystem::path manifestPath = directory.file( "radial_manifest.fits" );

    hciAnalyze::fitsHeaderT manifestHeader;
    REQUIRE( manifestHeader.append<int>( "KLIP PSF PRODUCT SCHEMA", 1, "schema" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<std::string>( "KLIP PSF PRODUCT", "MANIFEST", "role" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<int>( "KLIP PSF COMPLETE", 1, "complete" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<int>( "KLIP PSF STAMP SIZE", 3, "stamp" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<float>( "KLIP PSF FILTER MIN GOOD FRACTION", 1, "support" ) ==
             mx::error_t::noerror );
    REQUIRE( manifestHeader.append<std::string>( "KLIP PSF SAMPLE RADII", "1,2,4", "radii" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<std::string>( "KLIP PSF SPATIAL MODEL", "REGION_RADIAL_LINEAR", "model" ) ==
             mx::error_t::noerror );
    REQUIRE( manifestHeader.append<std::string>( "NMODES", "3", "modes" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<std::string>( "REGMINR", "0,2", "region minima" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<std::string>( "REGMAXR", "2,5", "region maxima" ) == mx::error_t::noerror );
    hciAnalyze::imageT manifest( 1, 1 );
    manifest( 0, 0 ) = 1;
    writeFitsImage( manifestPath, manifest, &manifestHeader );

    hciAnalyze::fitsHeaderT responseHeader;
    REQUIRE( responseHeader.append<int>( "KLIP PSF PRODUCT SCHEMA", 1, "schema" ) == mx::error_t::noerror );
    REQUIRE( responseHeader.append<std::string>( "KLIP PSF PRODUCT", "RADIAL_RESPONSE", "role" ) ==
             mx::error_t::noerror );
    REQUIRE( responseHeader.append<int>( "KLIP PSF MODE COUNT", 3, "mode" ) == mx::error_t::noerror );
    hciAnalyze::cubeT responses( 3, 3, 3 );
    responses.setZero();
    responses.image( 0 )( 1, 1 ) = 2;
    responses.image( 1 )( 1, 1 ) = 4;
    responses.image( 2 )( 1, 1 ) = 8;
    writeFitsCube( directory.file( "radial_mode000_radial_response.fits" ), responses, &responseHeader );

    hciAnalyze::fitsHeaderT validityHeader;
    REQUIRE( validityHeader.append<int>( "KLIP PSF PRODUCT SCHEMA", 1, "schema" ) == mx::error_t::noerror );
    REQUIRE( validityHeader.append<std::string>( "KLIP PSF PRODUCT", "RADIAL_VALIDITY", "role" ) ==
             mx::error_t::noerror );
    REQUIRE( validityHeader.append<int>( "KLIP PSF MODE COUNT", 3, "mode" ) == mx::error_t::noerror );
    hciAnalyze::cubeT validities( 3, 3, 3 );
    validities.cube().setOnes();
    writeFitsCube( directory.file( "radial_mode000_radial_validity.fits" ), validities, &validityHeader );

    appHarness application;
    application.m_psfResponse = manifestPath.string();
    hciAnalyze::cubeT science( 11, 11, 1 );
    science.setZero();
    science.image( 0 )( 5, 6 ) = 6;
    science.image( 0 )( 5, 8 ) = 18;
    science.image( 0 )( 5, 9 ) = 16;
    hciAnalyze::fitsHeaderT scienceHeader;
    REQUIRE( scienceHeader.append<std::string>( "NMODES", "3", "modes" ) == mx::error_t::noerror );
    REQUIRE( scienceHeader.append<std::string>( "REGMINR", "0,2", "region minima" ) == mx::error_t::noerror );
    REQUIRE( scienceHeader.append<std::string>( "REGMAXR", "2,5", "region maxima" ) == mx::error_t::noerror );

    application.filterCubePSFResponse( science, scienceHeader );
    REQUIRE( science.image( 0 )( 5, 6 ) == Approx( 3 ) );
    REQUIRE( science.image( 0 )( 5, 8 ) == Approx( 3 ) );
    REQUIRE( science.image( 0 )( 5, 9 ) == Approx( 2 ) );
    REQUIRE( std::isnan( science.image( 0 )( 0, 0 ) ) );
    REQUIRE( application.m_psfResponseMinimumSupport == Approx( 1 ) );

    hciAnalyze::fitsHeaderT mismatchedHeader;
    REQUIRE( mismatchedHeader.append<std::string>( "NMODES", "4", "modes" ) == mx::error_t::noerror );
    REQUIRE( mismatchedHeader.append<std::string>( "REGMINR", "0,2", "region minima" ) == mx::error_t::noerror );
    REQUIRE( mismatchedHeader.append<std::string>( "REGMAXR", "2,5", "region maxima" ) == mx::error_t::noerror );
    hciAnalyze::cubeT mismatchedScience( 11, 11, 1 );
    mismatchedScience.setZero();
    REQUIRE_THROWS( application.filterCubePSFResponse( mismatchedScience, mismatchedHeader ) );

    manifest( 0, 0 ) = 0;
    writeFitsImage( manifestPath, manifest, &manifestHeader );
    hciAnalyze::cubeT rejectedScience( 11, 11, 1 );
    rejectedScience.setZero();
    REQUIRE_THROWS( application.filterCubePSFResponse( rejectedScience, scienceHeader ) );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::RadialPSFModel doxygenRadialModel( { 1 }, 1, 1 );
    doxygenRadialModel.fit( {}, {}, {} );
    doxygenRadialModel.response( {}, {}, 1, 0 );
    mx::improc::P4PSFFilter::calculate( {}, {}, {}, 0, 0, 1 );
#endif
    // clang-format on
}

/// Verify hciAnalyze preserves an exact KLIP response at its anchor angle and rotates it at other azimuths.
/** This exercises hciAnalyze::filterCubePSFResponse(), mx::improc::RadialPSFModel::rotate(), and
 * mx::improc::P4PSFFilter::calculate().
 * \ingroup hciAnalyze_unit_tests
 */
TEST_CASE( "hciAnalyze external exact azimuthal KLIP PSF filtering", "[hciAnalyze][filter][klipPSF][exact]" )
{
    TestDirectory directory;
    const std::filesystem::path manifestPath = directory.file( "exact_manifest.fits" );

    hciAnalyze::fitsHeaderT manifestHeader;
    REQUIRE( manifestHeader.append<int>( "KLIP PSF PRODUCT SCHEMA", 1, "schema" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<std::string>( "KLIP PSF PRODUCT", "MANIFEST", "role" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<int>( "KLIP PSF COMPLETE", 1, "complete" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<int>( "KLIP PSF STAMP SIZE", 3, "stamp" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<float>( "KLIP PSF FILTER MIN GOOD FRACTION", 1, "support" ) ==
             mx::error_t::noerror );
    REQUIRE( manifestHeader.append<std::string>( "KLIP PSF SAMPLE RADII", "1", "radius" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<std::string>( "KLIP PSF SPATIAL MODEL", "EXACT_AZIMUTHAL", "model" ) ==
             mx::error_t::noerror );
    REQUIRE( manifestHeader.append<double>( "KLIP PSF SAMPLE ANGLE", 0, "angle" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<int>( "KLIP PSF EXACT ROW", 4, "row" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<int>( "KLIP PSF EXACT COLUMN", 5, "column" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<std::string>( "NMODES", "3", "modes" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<std::string>( "REGMINR", "0", "region minimum" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<std::string>( "REGMAXR", "3", "region maximum" ) == mx::error_t::noerror );
    hciAnalyze::imageT manifest( 1, 1 );
    manifest( 0, 0 ) = 1;
    writeFitsImage( manifestPath, manifest, &manifestHeader );

    hciAnalyze::fitsHeaderT responseHeader;
    REQUIRE( responseHeader.append<int>( "KLIP PSF PRODUCT SCHEMA", 1, "schema" ) == mx::error_t::noerror );
    REQUIRE( responseHeader.append<std::string>( "KLIP PSF PRODUCT", "RADIAL_RESPONSE", "role" ) ==
             mx::error_t::noerror );
    REQUIRE( responseHeader.append<int>( "KLIP PSF MODE COUNT", 3, "mode" ) == mx::error_t::noerror );
    hciAnalyze::cubeT responses( 3, 3, 1 );
    responses.setZero();
    responses.image( 0 )( 1, 1 ) = 2;
    responses.image( 0 )( 1, 2 ) = 1;
    writeFitsCube( directory.file( "exact_mode000_radial_response.fits" ), responses, &responseHeader );

    hciAnalyze::fitsHeaderT validityHeader;
    REQUIRE( validityHeader.append<int>( "KLIP PSF PRODUCT SCHEMA", 1, "schema" ) == mx::error_t::noerror );
    REQUIRE( validityHeader.append<std::string>( "KLIP PSF PRODUCT", "RADIAL_VALIDITY", "role" ) ==
             mx::error_t::noerror );
    REQUIRE( validityHeader.append<int>( "KLIP PSF MODE COUNT", 3, "mode" ) == mx::error_t::noerror );
    hciAnalyze::cubeT validities( 3, 3, 1 );
    validities.cube().setOnes();
    writeFitsCube( directory.file( "exact_mode000_radial_validity.fits" ), validities, &validityHeader );

    appHarness application;
    application.m_psfResponse = manifestPath.string();
    hciAnalyze::cubeT science( 9, 9, 1 );
    science.setZero();
    science.image( 0 )( 4, 5 ) = 6;
    science.image( 0 )( 4, 6 ) = 3;
    science.image( 0 )( 5, 4 ) = 6;
    science.image( 0 )( 6, 4 ) = 3;
    hciAnalyze::fitsHeaderT scienceHeader;
    REQUIRE( scienceHeader.append<std::string>( "NMODES", "3", "modes" ) == mx::error_t::noerror );
    REQUIRE( scienceHeader.append<std::string>( "REGMINR", "0", "region minimum" ) == mx::error_t::noerror );
    REQUIRE( scienceHeader.append<std::string>( "REGMAXR", "3", "region maximum" ) == mx::error_t::noerror );

    application.filterCubePSFResponse( science, scienceHeader );
    REQUIRE( science.image( 0 )( 4, 5 ) == Approx( 3 ) );
    REQUIRE( science.image( 0 )( 5, 4 ) == Approx( 3 ) );

    manifestHeader.erase( "KLIP PSF SAMPLE ANGLE" );
    writeFitsImage( manifestPath, manifest, &manifestHeader );
    hciAnalyze::cubeT rejectedScience( 9, 9, 1 );
    rejectedScience.setZero();
    REQUIRE_THROWS( application.filterCubePSFResponse( rejectedScience, scienceHeader ) );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::RadialPSFModel::rotate( {}, {}, {}, {}, 0 );
    mx::improc::P4PSFFilter::calculate( {}, {}, {}, 0, 0, 1 );
#endif
    // clang-format on
}

/// Verify hciAnalyze applies each persisted exact KLIP response only at its integer source coordinate.
/** This exercises hciAnalyze::filterCubePSFResponse() and mx::improc::P4PSFFilter::calculate().
 * \ingroup hciAnalyze_unit_tests
 */
TEST_CASE( "hciAnalyze external per-pixel KLIP PSF filtering", "[hciAnalyze][filter][klipPSF][pixel]" )
{
    TestDirectory directory;
    const std::filesystem::path manifestPath = directory.file( "pixel_manifest.fits" );

    hciAnalyze::fitsHeaderT manifestHeader;
    REQUIRE( manifestHeader.append<int>( "KLIP PSF PRODUCT SCHEMA", 2, "schema" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<std::string>( "KLIP PSF PRODUCT", "MANIFEST", "role" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<int>( "KLIP PSF COMPLETE", 1, "complete" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<int>( "KLIP PSF STAMP SIZE", 3, "stamp" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<float>( "KLIP PSF FILTER MIN GOOD FRACTION", 1, "support" ) ==
             mx::error_t::noerror );
    REQUIRE( manifestHeader.append<std::string>( "KLIP PSF SAMPLE RADII", "", "radii" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<std::string>( "KLIP PSF SPATIAL MODEL", "PIXEL_EXACT", "model" ) ==
             mx::error_t::noerror );
    REQUIRE( manifestHeader.append<std::string>( "KLIP PSF MEASUREMENT COUNT", "2", "count" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<std::string>( "NMODES", "3", "modes" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<std::string>( "REGMINR", "0", "region minimum" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<std::string>( "REGMAXR", "5", "region maximum" ) == mx::error_t::noerror );
    hciAnalyze::imageT manifest( 1, 1 );
    manifest( 0, 0 ) = 1;
    writeFitsImage( manifestPath, manifest, &manifestHeader );

    hciAnalyze::fitsHeaderT coordinateHeader;
    REQUIRE( coordinateHeader.append<int>( "KLIP PSF PRODUCT SCHEMA", 2, "schema" ) == mx::error_t::noerror );
    REQUIRE( coordinateHeader.append<std::string>( "KLIP PSF PRODUCT", "PIXEL_COORDINATES", "role" ) ==
             mx::error_t::noerror );
    hciAnalyze::imageT coordinates( 2, 4 );
    coordinates << 5, 7, 0, 0, 7, 5, 0, 1;
    writeFitsImage( directory.file( "pixel_coordinates.fits" ), coordinates, &coordinateHeader );

    hciAnalyze::fitsHeaderT responseHeader;
    REQUIRE( responseHeader.append<int>( "KLIP PSF PRODUCT SCHEMA", 2, "schema" ) == mx::error_t::noerror );
    REQUIRE( responseHeader.append<std::string>( "KLIP PSF PRODUCT", "PIXEL_RESPONSE", "role" ) ==
             mx::error_t::noerror );
    REQUIRE( responseHeader.append<int>( "KLIP PSF MODE COUNT", 3, "mode" ) == mx::error_t::noerror );
    hciAnalyze::cubeT responses( 3, 3, 2 );
    responses.setZero();
    responses.image( 0 )( 1, 1 ) = 2;
    responses.image( 1 )( 1, 1 ) = 4;
    writeFitsCube( directory.file( "pixel_mode000_pixel_response.fits" ), responses, &responseHeader );

    hciAnalyze::fitsHeaderT validityHeader;
    REQUIRE( validityHeader.append<int>( "KLIP PSF PRODUCT SCHEMA", 2, "schema" ) == mx::error_t::noerror );
    REQUIRE( validityHeader.append<std::string>( "KLIP PSF PRODUCT", "PIXEL_VALIDITY", "role" ) ==
             mx::error_t::noerror );
    REQUIRE( validityHeader.append<int>( "KLIP PSF MODE COUNT", 3, "mode" ) == mx::error_t::noerror );
    hciAnalyze::cubeT validities( 3, 3, 2 );
    validities.cube().setOnes();
    writeFitsCube( directory.file( "pixel_mode000_pixel_validity.fits" ), validities, &validityHeader );

    appHarness application;
    application.m_psfResponse = manifestPath.string();
    hciAnalyze::cubeT science( 11, 11, 1 );
    science.setZero();
    science.image( 0 )( 5, 7 ) = 6;
    science.image( 0 )( 7, 5 ) = 20;
    hciAnalyze::fitsHeaderT scienceHeader;
    REQUIRE( scienceHeader.append<std::string>( "NMODES", "3", "modes" ) == mx::error_t::noerror );
    REQUIRE( scienceHeader.append<std::string>( "REGMINR", "0", "region minimum" ) == mx::error_t::noerror );
    REQUIRE( scienceHeader.append<std::string>( "REGMAXR", "5", "region maximum" ) == mx::error_t::noerror );

    application.filterCubePSFResponse( science, scienceHeader );
    REQUIRE( science.image( 0 )( 5, 7 ) == Approx( 3 ) );
    REQUIRE( science.image( 0 )( 7, 5 ) == Approx( 5 ) );
    REQUIRE( std::isnan( science.image( 0 )( 5, 6 ) ) );

    coordinates( 1, 3 ) = 0;
    writeFitsImage( directory.file( "pixel_coordinates.fits" ), coordinates, &coordinateHeader );
    hciAnalyze::cubeT rejectedScience( 11, 11, 1 );
    rejectedScience.setZero();
    REQUIRE_THROWS( application.filterCubePSFResponse( rejectedScience, scienceHeader ) );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::P4PSFFilter::calculate( {}, {}, {}, 0, 0, 1 );
#endif
    // clang-format on
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

/** \brief Verify hciAnalyze validates and loads covariance geometry and regularization options.
 * \ingroup hciAnalyze_unit_tests
 */
TEST_CASE( "hciAnalyze covariance configuration", "[hciAnalyze][noise][config]" )
{
    appHarness application;
    application.setupConfig();
    std::vector<std::string> options{ "hciAnalyze",
                                      "--noise.model=pca",
                                      "--noise.minimumSamples=12",
                                      "--noise.maximumModes=2",
                                      "--noise.floorFraction=0.25",
                                      "--noise.arcStep=1.25",
                                      "--noise.guardRadius=2",
                                      "--noise.exactExclusion=true",
                                      "--noise.only=true",
                                      "--lambdaD=2.5" };
    std::vector<char *> arguments;
    for( auto &option : options )
    {
        arguments.push_back( option.data() );
    }
    application.config.parseCommandLine( static_cast<int>( arguments.size() ), arguments.data() );
    application.loadConfig();
    REQUIRE( application.m_noiseModel == "pca" );
    REQUIRE( application.m_noiseConfig.m_minimumSamples == 12 );
    REQUIRE( application.m_noiseConfig.m_maximumModes == 2 );
    REQUIRE( application.m_noiseConfig.m_floorFraction == 0.25 );
    REQUIRE( application.m_noiseConfig.m_arcStep == 1.25 );
    REQUIRE( application.m_noiseConfig.m_guardRadius == 2 );
    REQUIRE( application.m_noiseConfig.m_exactExclusion );
    REQUIRE( application.m_noiseOnly );
    application.m_file = "input.fits";
    REQUIRE( application.config.m_targets.count( "noise.model" ) == 1 );
    REQUIRE( application.config.m_targets.count( "noise.floorFraction" ) == 1 );
    application.m_noiseModel = "pca";
    REQUIRE_THROWS( application.checkConfig() );
    application.m_psfResponse = "field_manifest.fits";
    REQUIRE_NOTHROW( application.checkConfig() );
    REQUIRE( application.m_noiseConfig.m_kind == mx::improc::PSFNoiseKind::pca );
    application.m_highPassFwhm = 2;
    REQUIRE_THROWS( application.checkConfig() );
    application.m_highPassFwhm = 0;
    application.m_noiseConfig.m_floorFraction = 0;
    REQUIRE_THROWS( application.checkConfig() );
    application.m_noiseConfig.m_floorFraction = 0.1;
    application.m_noiseConfig.m_minimumSamples = 1;
    REQUIRE_THROWS( application.checkConfig() );
    application.m_noiseConfig.m_minimumSamples = 8;
    application.m_noiseConfig.m_maximumModes = std::numeric_limits<std::size_t>::max();
    REQUIRE_THROWS( application.checkConfig() );
    application.m_noiseConfig.m_maximumModes = 3;
    application.m_noiseExcludeRows = { 20 };
    REQUIRE_THROWS( application.checkConfig() );
    application.m_noiseExcludeColumns = { 30 };
    application.m_noiseExcludeRadii = { -1 };
    REQUIRE_THROWS( application.checkConfig() );
    application.m_noiseExcludeRadii = { 5 };
    REQUIRE_NOTHROW( application.checkConfig() );
    application.m_noiseModel = "unknown";
    REQUIRE_THROWS( application.checkConfig() );
    // clang-format off
#ifdef __DOXY_ONLY__
    hciAnalyze::checkConfig();
    hciAnalyze::checkNoiseConfig();
#endif
    // clang-format on
}

/** \brief Verify hciAnalyze routes P4 response filtering through PSFNoiseTraining and writes conditional diagnostics.
 * \ingroup hciAnalyze_unit_tests
 */
TEST_CASE( "hciAnalyze writes covariance filter products", "[hciAnalyze][noise][products]" )
{
    // clang-format off
#ifdef __DOXY_ONLY__
    hciAnalyze::execute();
    hciAnalyze::filterCubePSFResponse( {}, {} );
    hciAnalyze::prepareNoiseProducts( {} );
    hciAnalyze::applyPSFFilter( {}, {}, {}, 0, 0, 0 );
    hciAnalyze::writeNoiseProducts( {} );
#endif
    // clang-format on
    TestDirectory directory;
    const std::filesystem::path manifestPath = directory.file( "field_manifest.fits" );

    hciAnalyze::fitsHeaderT manifestHeader;
    REQUIRE( manifestHeader.append<int>( "P4 PSF PRODUCT SCHEMA", 2, "schema" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<std::string>( "P4 PSF PRODUCT", "MANIFEST", "role" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<int>( "P4 PSF COMPLETE", 1, "complete" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<int>( "P4 PSF MODE COUNT", 1, "modes" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<std::string>( "P4 PSF SOURCE COUNT", "1", "sources" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<int>( "P4 PSF STAMP SIZE", 3, "stamp" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<float>( "P4 PSF FILTER MIN GOOD FRACTION", 1, "support" ) == mx::error_t::noerror );
    REQUIRE( manifestHeader.append<std::string>( "P4 MODE FRACTIONS", "0.5", "modes" ) == mx::error_t::noerror );
    hciAnalyze::imageT manifest( 1, 1 );
    manifest( 0, 0 ) = 1;
    writeFitsImage( manifestPath, manifest, &manifestHeader );

    hciAnalyze::fitsHeaderT coordinateHeader;
    REQUIRE( coordinateHeader.append<std::string>( "P4 PSF PRODUCT", "COORDINATES", "role" ) == mx::error_t::noerror );
    hciAnalyze::imageT coordinates( 1, 4 );
    coordinates << 48, 32, 16, -90;
    writeFitsImage( directory.file( "field_coordinates.fits" ), coordinates, &coordinateHeader );

    hciAnalyze::fitsHeaderT modelHeader;
    REQUIRE( modelHeader.append<std::string>( "P4 PSF PRODUCT", "MODEL", "role" ) == mx::error_t::noerror );
    REQUIRE( modelHeader.append<int>( "P4 PSF MODE INDEX", 0, "mode" ) == mx::error_t::noerror );
    hciAnalyze::cubeT models( 3, 3, 1 );
    models.setZero();
    models.image( 0 )( 1, 1 ) = 2;
    writeFitsCube( directory.file( "field_model_0000.fits" ), models, &modelHeader );

    hciAnalyze::fitsHeaderT validityHeader;
    REQUIRE( validityHeader.append<std::string>( "P4 PSF PRODUCT", "VALIDITY", "role" ) == mx::error_t::noerror );
    REQUIRE( validityHeader.append<int>( "P4 PSF MODE INDEX", 0, "mode" ) == mx::error_t::noerror );
    hciAnalyze::imageT validity( 1, 1 );
    validity( 0, 0 ) = 1;
    writeFitsImage( directory.file( "field_validity_0000.fits" ), validity, &validityHeader );

    hciAnalyze::cubeT original( 65, 65, 1 );
    for( int column = 0; column < 65; ++column )
    {
        for( int row = 0; row < 65; ++row )
        {
            original.image( 0 )( row, column ) = std::sin( 0.3 * row ) + std::cos( 0.2 * column );
        }
    }
    hciAnalyze::fitsHeaderT scienceHeader;
    REQUIRE( scienceHeader.append<std::string>( "P4 MODE FRACTIONS", "0.5", "modes" ) == mx::error_t::noerror );
    const auto responseValidity = mx::improc::P4PSFFilter::validityT::Ones( 3, 3 ).eval();
    for( const std::string model : { "identity", "diagonal", "pca" } )
    {
        appHarness application;
        application.m_psfResponse = manifestPath.string();
        application.m_file = directory.file( model + ".fits" ).string();
        application.m_noiseModel = model;
        application.m_noiseDiagnostics = true;
        application.m_noiseConfig.m_exactExclusion = model == "pca";
        application.m_noiseExcludeRows = { 32 };
        application.m_noiseExcludeColumns = { 48 };
        application.m_noiseExcludeRadii = { 3 };
        application.checkConfig();
        const mx::improc::PSFNoiseTrainingResult expected =
            mx::improc::PSFNoiseTraining::calculate( original.image( 0 ),
                                                     models.image( 0 ),
                                                     responseValidity,
                                                     48,
                                                     32,
                                                     1,
                                                     application.m_noiseConfig,
                                                     { { 32, 48, 3 } } );
        REQUIRE( expected.m_filter.valid );
        hciAnalyze::cubeT science = original;
        application.filterCubePSFResponse( science, scienceHeader );
        REQUIRE( science.image( 0 )( 48, 32 ) == Approx( expected.m_filter.amplitude ).epsilon( 1e-6 ) );
        REQUIRE( std::isnan( science.image( 0 )( 0, 0 ) ) );
        const std::array<std::pair<const char *, double>, 11> checks{
            std::pair{ "psf_amplitude", expected.m_filter.amplitude },
            { "psf_sigma", expected.m_filter.conditionalSigma },
            { "psf_score", expected.m_filter.score },
            { "noise_samples", static_cast<double>( expected.m_samples ) },
            { "noise_modes", static_cast<double>( expected.m_modes ) },
            { "noise_floor", expected.m_varianceFloor },
            { "noise_status", 0 },
            { "psf_support", 1 },
            { "noise_attempted", static_cast<double>( expected.m_attempted ) },
            { "noise_excluded", static_cast<double>( expected.m_excluded ) },
            { "noise_incomplete", static_cast<double>( expected.m_incomplete ) } };
        mx::fits::fitsFile<float, mx::verbose::vv> reader;
        for( const auto &[role, value] : checks )
        {
            hciAnalyze::cubeT product;
            hciAnalyze::fitsHeaderT header;
            REQUIRE( reader.read( product, header, directory.file( model + "_" + role + ".fits" ).string() ) ==
                     mx::error_t::noerror );
            REQUIRE( product.image( 0 )( 48, 32 ) == Approx( value ).epsilon( 1e-6 ) );
            REQUIRE( std::isnan( product.image( 0 )( 0, 0 ) ) );
            std::string recordedModel = header["HCIA NOISE MODEL"].String();
            recordedModel.erase( recordedModel.find_last_not_of( ' ' ) + 1 );
            REQUIRE( recordedModel == model );
            REQUIRE( header["HCIA NOISE EXCLUSION"].String().find(
                         model == "pca" ? "NONZERO_STENCIL_PIXEL_CENTERS" : "CIRCUMCIRCLE_PLUS_SQRT2" ) == 0 );
            REQUIRE( header["HCIA NOISE CALIBRATED"].value<int>() == 0 );
            REQUIRE( header["SNRSMALL"].value<int>() == 0 );
            std::string recordedExclusions = header["HCIA NOISE EXCLUSIONS"].String();
            recordedExclusions.erase( recordedExclusions.find_last_not_of( ' ' ) + 1 );
            REQUIRE( recordedExclusions == "32,48,3" );
        }
        // A declared source can have no response/training support in a sparse field. Conditional-only execution
        // still writes its explicit statuses without requiring an unrelated annular SNR aperture measurement.
        writeFitsCube( directory.file( model + ".fits" ), original, &scienceHeader );
        application.m_noiseOnly = true;
        application.m_noiseDiagnostics = false;
        application.m_planetSpecified = true;
        application.m_planetSeparation = { 16 };
        application.m_planetPositionAngle = { 0 };
        application.m_planetRadius = { 3 };
        REQUIRE( application.execute() == EXIT_SUCCESS );
        REQUIRE( application.m_noiseDiagnostics );
        REQUIRE( application.m_results.empty() );
        REQUIRE_FALSE( std::filesystem::exists( directory.file( model + "_snr.fits" ) ) );
        if( model == "identity" )
        {
            continue;
        }
        // Too little training is explicit and does not silently produce an identity-filter amplitude.
        application.m_noiseConfig.m_minimumSamples = 1000;
        science = original;
        application.filterCubePSFResponse( science, scienceHeader );
        REQUIRE( std::isnan( science.image( 0 )( 48, 32 ) ) );
        hciAnalyze::cubeT status;
        hciAnalyze::fitsHeaderT statusHeader;
        REQUIRE( reader.read( status, statusHeader, directory.file( model + "_noise_status.fits" ).string() ) ==
                 mx::error_t::noerror );
        REQUIRE( status.image( 0 )( 48, 32 ) == 1 );
    }
}

} // namespace hciAnalyze_test
} // namespace unitTest
