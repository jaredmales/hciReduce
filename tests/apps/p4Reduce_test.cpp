/** \file p4Reduce_test.cpp
 * \brief Tests p4Reduce configuration, dispatch, FITS integration, and process errors.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "../common/HCIobservation_test_fixture.hpp"

#include <cstdlib>
#include <filesystem>
#include <sstream>
#include <string>

#define HCIREDUCE_P4REDUCE_NO_MAIN
#include "src/apps/p4Reduce.cpp"

namespace unitTest
{
namespace p4Reduce_test
{

/// \cond p4Reduce_test_harness
using appT = p4Reduce;

struct appHarness : public appT
{
    using appT::checkConfig;
    using appT::config;
    using appT::execute;
    using appT::loadConfig;
    using appT::m_mode;
    using appT::m_obs;
    using appT::m_showTiming;
    using appT::setupConfig;
};

/// Restore one output stream after capturing text written during a test scope.
class StreamCapture
{
  public:
    /// Redirect the supplied output stream into an owned string buffer.
    explicit StreamCapture( std::ostream &stream /**< [in,out] stream to redirect */ )
        : m_stream( stream ), m_buffer(), m_original( stream.rdbuf( m_buffer.rdbuf() ) )
    {
    }

    /// Restore the original output buffer.
    ~StreamCapture()
    {
        m_stream.rdbuf( m_original );
    }

    StreamCapture( const StreamCapture & ) = delete;
    StreamCapture &operator=( const StreamCapture & ) = delete;

    /// Return all text captured so far.
    std::string str() const
    {
        return m_buffer.str();
    }

  private:
    std::ostream &m_stream;      ///< Redirected output stream.
    std::ostringstream m_buffer; ///< Captured text.
    std::streambuf *m_original;  ///< Output buffer restored at destruction.
};

/// Construct a complete compact p4Reduce configuration.
std::string p4Configuration( const std::filesystem::path &inputDirectory,  /**< [in] target FITS directory */
                             const std::filesystem::path &outputDirectory, /**< [in] final-product directory */
                             const std::string &mode,                      /**< [in] basic or normal */
                             bool preprocessOnly /**< [in] whether to stop after preprocessing */ )
{
    std::ostringstream configuration;
    configuration << "mode=" << mode << '\n'
                  << "[input]\n"
                  << "directory=" << inputDirectory.string() << '\n'
                  << "prefix=target_\n"
                  << "extension=.fits\n"
                  << "dateKeyword=\n"
                  << "angleKeyword=ANGLE\n"
                  << "angleScale=1\n"
                  << "angleConstant=0\n"
                  << "[geom]\n"
                  << "minRadius=5\n"
                  << "maxRadius=6\n"
                  << "[p4]\n"
                  << "modeFractions=0.5\n"
                  << "regressionFrame=detector\n"
                  << "numberImages=0\n"
                  << "orDeltaRadiusInner=2\n"
                  << "orDeltaRadiusOuter=2\n"
                  << "orArcHalfWidth=4\n"
                  << "orMaxHalfAngle=60\n"
                  << "psfRadius=0.5\n"
                  << "exclusionPolicy=kernelSupport\n"
                  << "exclusionRadiusBuffer=0\n"
                  << "rankTolerance=1e-12\n"
                  << "writeDiagnostics=false\n"
                  << "[preProcess]\n"
                  << "only=" << ( preprocessOnly ? "true" : "false" ) << '\n'
                  << "[combine]\n"
                  << "method=mean\n"
                  << "minGoodFract=0\n"
                  << "[output]\n"
                  << "directory=" << outputDirectory.string() << '\n'
                  << "fileName=p4-final.fits\n"
                  << "exactFName=true\n";
    return configuration.str();
}

/// Construct the AF Lep comparison geometry configuration without depending on an untracked example file.
std::string p4AfLepConfiguration( const std::filesystem::path &inputDirectory /**< [in] target FITS directory */ )
{
    std::ostringstream configuration;
    configuration << "mode=normal\n"
                  << "[input]\n"
                  << "directory=" << inputDirectory.string() << '\n'
                  << "prefix=pp\n"
                  << "extension=.fits\n"
                  << "fileList=\n"
                  << "dateKeyword=\n"
                  << "maskFile=\n"
                  << "angleKeyword=POSANG\n"
                  << "angleScale=1\n"
                  << "angleConstant=1.6\n"
                  << "[rdi]\n"
                  << "directory=\n"
                  << "prefix=\n"
                  << "fileList=\n"
                  << "[coadd]\n"
                  << "method=none\n"
                  << "[preProcess]\n"
                  << "skip=true\n"
                  << "meanSubMethod=none\n"
                  << "pixelTSNormMethod=none\n"
                  << "only=false\n"
                  << "[adi]\n"
                  << "postMedSub=false\n"
                  << "[geom]\n"
                  << "minRadius=3.5,6,8,10,12,14,16,18\n"
                  << "maxRadius=6,8,10,12,14,16,18,20\n"
                  << "[p4]\n"
                  << "modeFractions=0.01,0.03,0.05,0.07,0.09,0.11,0.13,0.15,0.17,0.19,0.21,0.23,0.25,0.27,0.29\n"
                  << "regressionFrame=detector\n"
                  << "orDeltaRadiusInner=8\n"
                  << "orDeltaRadiusOuter=15\n"
                  << "orArcHalfWidth=32\n"
                  << "orMaxHalfAngle=90\n"
                  << "psfRadius=1.5\n"
                  << "exclusionPolicy=sampleCenter\n"
                  << "exclusionRadiusBuffer=0\n"
                  << "rankTolerance=0\n"
                  << "writeDiagnostics=false\n"
                  << "[combine]\n"
                  << "method=sigmaMean\n"
                  << "weightFile=\n"
                  << "sigmaThreshold=5\n"
                  << "minGoodFract=0\n"
                  << "noDerotate=false\n"
                  << "[output]\n"
                  << "directory=.\n"
                  << "fileName=finim_\n"
                  << "exactFName=false\n"
                  << "outputPSFSub=false\n";
    return configuration.str();
}

/// Write one constant P4 target with a configured ADI angle.
void writeTarget( const std::filesystem::path &path, /**< [in] output FITS path */
                  float value,                       /**< [in] constant science value */
                  float angle /**< [in] ADI angle in degrees */ )
{
    HCIobservationTestHarness::imageT image( 31, 31 );
    image.setConstant( value );
    HCIobservationTestHarness::fitsHeaderT header;
    header.append<float>( "ANGLE", angle, "test parallactic angle" );
    writeFitsImage( path, image, &header );
}
/// \endcond

/// Verify p4Reduce registers and loads application and reduction configuration without grid or postprocess targets.
/** \ingroup p4Reduce_unit_tests */
TEST_CASE( "p4Reduce configuration registration", "[p4Reduce][config]" )
{
    TestDirectory directory;
    appHarness application;
    application.setupConfig();

    REQUIRE( application.m_mode == "basic" );
    REQUIRE( application.config.m_targets.at( "mode" ).clType == mx::app::argType::Required );
    REQUIRE( application.config.m_targets.at( "showTiming" ).helpType == "bool" );
    REQUIRE( application.config.m_targets.at( "p4.modeFractions" ).helpType == "vector<realT>" );
    REQUIRE( application.config.m_targets.at( "p4.regressionFrame" ).helpType == "string" );
    REQUIRE( application.config.m_targets.at( "p4.numberImages" ).helpType == "int" );
    REQUIRE( application.config.m_targets.at( "p4.memoryFraction" ).helpType == "double" );
    REQUIRE( application.config.m_targets.at( "p4.exclusionPolicy" ).clType == mx::app::argType::Required );

    for( const auto &[name, target] : application.config.m_targets )
    {
        static_cast<void>( target );
        REQUIRE_FALSE( name.starts_with( "grid." ) );
        REQUIRE_FALSE( name.starts_with( "postprocess." ) );
    }

    const auto configPath = directory.file( "normal.conf" );
    writeTextFile( configPath,
                   "showTiming=true\n" +
                       p4Configuration( directory.path(), directory.file( "output" ), "normal", true ) );
    REQUIRE( application.config.readConfig( configPath.string() ) == 0 );
    REQUIRE_NOTHROW( application.loadConfig() );
    REQUIRE( application.m_mode == "normal" );
    REQUIRE( application.m_showTiming );
    REQUIRE( application.m_obs.m_minRadius == std::vector<float>{ 5 } );
    REQUIRE( application.m_obs.m_maxRadius == std::vector<float>{ 6 } );
    REQUIRE( application.m_obs.m_modeFractions == std::vector<float>{ 0.5F } );
    REQUIRE( application.m_obs.m_regressionFrame == mx::improc::P4RegressionFrame::detector );
    REQUIRE( application.m_obs.m_numberImages == 0 );
    REQUIRE( application.m_obs.m_memoryFraction == Approx( 0.8 ) );
    REQUIRE( application.m_obs.m_exclusionPolicy == mx::improc::P4ExclusionPolicy::kernelSupport );

    // clang-format off
#ifdef __DOXY_ONLY__
    p4Reduce doxygenApplication;
    doxygenApplication.setupConfig();
    doxygenApplication.loadConfig();
#endif
    // clang-format on
}

/// Verify the AF Lep/NACO comparison configuration loads its support geometry and controlled settings.
/** \ingroup p4Reduce_unit_tests */
TEST_CASE( "p4Reduce AF Lep prototype configuration", "[p4Reduce][config][prototype][AF-Lep]" )
{
    TestDirectory directory;
    const auto configPath = directory.file( "aflep-comparison.conf" );
    writeTextFile( configPath, p4AfLepConfiguration( directory.path() ) );

    appHarness application;
    application.setupConfig();
    REQUIRE( application.config.readConfig( configPath.string() ) == 0 );
    REQUIRE( application.config.m_unusedConfigs.empty() );
    REQUIRE_NOTHROW( application.loadConfig() );
    REQUIRE_NOTHROW( application.checkConfig() );

    const auto configuredValue = [&application]( const std::string &name ) -> const std::string &
    {
        const auto &target = application.config.m_targets.at( name );
        REQUIRE( target.set );
        REQUIRE_FALSE( target.values.empty() );
        return target.values.back();
    };

    REQUIRE( application.m_mode == "normal" );
    REQUIRE( configuredValue( "input.directory" ) == directory.path().string() );
    REQUIRE( configuredValue( "input.prefix" ) == "pp" );
    REQUIRE( configuredValue( "input.fileList" ).empty() );
    REQUIRE( configuredValue( "input.dateKeyword" ).empty() );
    REQUIRE( configuredValue( "input.maskFile" ).empty() );
    REQUIRE( configuredValue( "rdi.directory" ).empty() );
    REQUIRE( configuredValue( "rdi.prefix" ).empty() );
    REQUIRE( configuredValue( "rdi.fileList" ).empty() );
    REQUIRE( configuredValue( "coadd.method" ) == "none" );
    REQUIRE( configuredValue( "preProcess.skip" ) == "true" );
    REQUIRE( configuredValue( "preProcess.meanSubMethod" ) == "none" );
    REQUIRE( configuredValue( "preProcess.pixelTSNormMethod" ) == "none" );
    REQUIRE( configuredValue( "preProcess.only" ) == "false" );
    REQUIRE( configuredValue( "p4.regressionFrame" ) == "detector" );
    REQUIRE( configuredValue( "output.fileName" ) == "finim_" );
    REQUIRE( configuredValue( "output.exactFName" ) == "false" );

    REQUIRE( application.m_obs.m_derotF.m_angleKeyword == "POSANG" );
    REQUIRE( application.m_obs.m_derotF.m_angleScale == 1 );
    REQUIRE( application.m_obs.m_derotF.m_angleConstant == 1.6F );
    REQUIRE_FALSE( application.m_obs.m_postMedSub );
    REQUIRE( application.m_obs.m_minRadius == std::vector<float>{ 3.5F, 6, 8, 10, 12, 14, 16, 18 } );
    REQUIRE( application.m_obs.m_maxRadius == std::vector<float>{ 6, 8, 10, 12, 14, 16, 18, 20 } );
    REQUIRE( application.m_obs.m_modeFractions == std::vector<float>{ 0.01F,
                                                                      0.03F,
                                                                      0.05F,
                                                                      0.07F,
                                                                      0.09F,
                                                                      0.11F,
                                                                      0.13F,
                                                                      0.15F,
                                                                      0.17F,
                                                                      0.19F,
                                                                      0.21F,
                                                                      0.23F,
                                                                      0.25F,
                                                                      0.27F,
                                                                      0.29F } );
    REQUIRE( application.m_obs.m_regressionFrame == mx::improc::P4RegressionFrame::detector );
    REQUIRE( application.m_obs.m_orDeltaRadiusInner == 8 );
    REQUIRE( application.m_obs.m_orDeltaRadiusOuter == 15 );
    REQUIRE( application.m_obs.m_orArcHalfWidth == 32 );
    REQUIRE( application.m_obs.m_orMaxHalfAngle == 90 );
    REQUIRE( application.m_obs.m_psfRadius == 1.5F );
    REQUIRE( application.m_obs.m_exclusionPolicy == mx::improc::P4ExclusionPolicy::sampleCenter );
    REQUIRE( application.m_obs.m_exclusionRadiusBuffer == 0 );
    REQUIRE( application.m_obs.m_rankTolerance == 0 );
    REQUIRE_FALSE( application.m_obs.m_writeDiagnostics );
    REQUIRE( application.m_obs.m_combineMethod == mx::improc::HCI::combine::sigmaMean );
    REQUIRE( application.m_obs.m_weightFile.empty() );
    REQUIRE( application.m_obs.m_sigmaThreshold == 5 );
    REQUIRE( application.m_obs.m_minGoodFract == 0 );
    REQUIRE( application.m_obs.m_doDerotate );
    REQUIRE( application.m_obs.m_outputDir == "." );
    REQUIRE( application.m_obs.m_finimName == "finim_" );
    REQUIRE_FALSE( application.m_obs.m_exactFinimName );
    REQUIRE_FALSE( application.m_obs.m_doOutputPSFSub );
}

/// Verify p4Reduce reports unused targets and rejects unknown values and positional command-line inputs.
/** \ingroup p4Reduce_unit_tests */
TEST_CASE( "p4Reduce configuration diagnostics", "[p4Reduce][config][diagnostics]" )
{
    SECTION( "source locations enabled" )
    {
        appHarness application;
        application.setupConfig();

        mx::app::configTarget unusedTarget;
        unusedTarget.name = "unused.target";
        application.config.m_targets[unusedTarget.name] = unusedTarget;

        mx::app::configTarget sourcedUnknown;
        sourcedUnknown.name = "unknown.sourced";
        sourcedUnknown.sources = { "unknown.conf" };
        application.config.m_unusedConfigs[sourcedUnknown.name] = sourcedUnknown;

        mx::app::configTarget unsourcedUnknown;
        unsourcedUnknown.name = "unknown.unsourced";
        application.config.m_unusedConfigs[unsourcedUnknown.name] = unsourcedUnknown;
        application.config.nonOptions = { "unexpected" };

        StreamCapture errors( std::cerr );
        REQUIRE_THROWS( application.loadConfig() );
        const std::string output = errors.str();
        REQUIRE( output.find( "unused.target" ) != std::string::npos );
        REQUIRE( output.find( "unknown.sourced [unknown.conf]" ) != std::string::npos );
        REQUIRE( output.find( "unknown.unsourced" ) != std::string::npos );
        REQUIRE( output.find( "unrecognized command line arguments" ) != std::string::npos );
    }

    SECTION( "source locations disabled" )
    {
        appHarness application;
        application.setupConfig();
        application.config.m_sources = false;

        mx::app::configTarget unknown;
        unknown.name = "unknown.target";
        unknown.sources = { "must-not-print.conf" };
        application.config.m_unusedConfigs[unknown.name] = unknown;

        StreamCapture errors( std::cerr );
        REQUIRE_THROWS( application.loadConfig() );
        const std::string output = errors.str();
        REQUIRE( output.find( "unknown.target" ) != std::string::npos );
        REQUIRE( output.find( "must-not-print.conf" ) == std::string::npos );
    }

    SECTION( "dotted config-file key" )
    {
        TestDirectory directory;
        const auto configPath = directory.file( "dotted.conf" );
        writeTextFile( configPath, "p4.psfFile=/tmp/ignored.fits\n" );

        appHarness application;
        application.setupConfig();
        REQUIRE( application.config.readConfig( configPath.string() ) == 0 );
        REQUIRE_FALSE( application.config.m_unusedConfigs.empty() );

        StreamCapture errors( std::cerr );
        REQUIRE_THROWS( application.loadConfig() );
        const std::string output = errors.str();
        REQUIRE( output.find( "p4.psfFile" ) != std::string::npos );
        REQUIRE( output.find( configPath.string() ) != std::string::npos );
    }
}

/// Verify p4Reduce accepts only the two equivalent configured-run modes in the initial release.
/** \ingroup p4Reduce_unit_tests */
TEST_CASE( "p4Reduce mode validation", "[p4Reduce][config][mode]" )
{
    appHarness application;

    application.m_mode = "basic";
    REQUIRE_NOTHROW( application.checkConfig() );
    application.m_mode = "normal";
    REQUIRE_NOTHROW( application.checkConfig() );

    for( const std::string &invalid : { "postprocess", "grid", "unknown" } )
    {
        application.m_mode = invalid;
        REQUIRE_THROWS( application.checkConfig() );
    }

    // clang-format off
#ifdef __DOXY_ONLY__
    p4Reduce doxygenApplication;
    doxygenApplication.checkConfig();
#endif
    // clang-format on
}

/// Verify command-line P4 values override values loaded from the selected configuration file.
/** \ingroup p4Reduce_unit_tests */
TEST_CASE( "p4Reduce command-line precedence", "[p4Reduce][config][command-line][precedence]" )
{
    OpenMPThreadGuard threads( 1 );
    TestDirectory directory;
    writeTarget( directory.file( "target_000.fits" ), 1, 4 );

    const auto configPath = directory.file( "precedence.conf" );
    writeTextFile( configPath, p4Configuration( directory.path(), directory.file( "output" ), "basic", true ) );

    appHarness application;
    std::string invokedName = "p4Reduce-test";
    std::string configOption = "--config";
    std::string configName = configPath.string();
    std::string modeOption = "--mode";
    std::string modeValue = "normal";
    std::string frameOption = "--p4.regressionFrame";
    std::string frameValue = "rotated";
    char *arguments[]{ invokedName.data(),
                       configOption.data(),
                       configName.data(),
                       modeOption.data(),
                       modeValue.data(),
                       frameOption.data(),
                       frameValue.data() };

    REQUIRE( application.main( 7, arguments ) == 0 );
    REQUIRE( application.m_mode == "normal" );
    REQUIRE( application.m_obs.m_regressionFrame == mx::improc::P4RegressionFrame::rotated );
    REQUIRE( application.m_obs.m_derotF.m_angles == std::vector<float>{ 4 } );
}

/// Verify p4Reduce basic mode discovers real FITS input and stops at the inherited preprocessing-only gate.
/** \ingroup p4Reduce_unit_tests */
TEST_CASE( "p4Reduce basic preprocess-only dispatch", "[p4Reduce][execute][basic][preprocess][FITS]" )
{
    OpenMPThreadGuard threads( 1 );
    TestDirectory directory;
    writeTarget( directory.file( "target_000.fits" ), 1, 12 );

    const auto configPath = directory.file( "basic.conf" );
    writeTextFile( configPath, p4Configuration( directory.path(), directory.file( "output" ), "basic", true ) );

    appHarness application;
    std::string invokedName = "p4Reduce-test";
    std::string configOption = "--config";
    std::string configName = configPath.string();
    char *arguments[]{ invokedName.data(), configOption.data(), configName.data() };

    REQUIRE( application.main( 3, arguments ) == 0 );
    REQUIRE( application.m_mode == "basic" );
    REQUIRE( application.m_obs.m_derotF.m_angles == std::vector<float>{ 12 } );
    REQUIRE( application.m_obs.m_psfsub.empty() );

    // clang-format off
#ifdef __DOXY_ONLY__
    p4Reduce doxygenApplication;
    doxygenApplication.execute();
    mx::improc::P4Reductionf doxygenReduction;
    doxygenReduction.reduce();
#endif
    // clang-format on
}

/// Verify p4Reduce preprocessing-only execution does not require P4 reduction geometry.
/** ingroup p4Reduce_unit_tests */
TEST_CASE( "p4Reduce preprocess-only without P4 geometry", "[p4Reduce][execute][preprocess][geometry]" )
{
    OpenMPThreadGuard threads( 1 );
    TestDirectory directory;
    writeTarget( directory.file( "target_000.fits" ), 1, 12 );

    const auto configPath = directory.file( "preprocess-only.conf" );
    writeTextFile( configPath,
                   "[input]\n"
                   "directory=" +
                       directory.path().string() +
                       "\n"
                       "prefix=target_\n"
                       "extension=.fits\n"
                       "dateKeyword=\n"
                       "angleKeyword=ANGLE\n"
                       "angleScale=1\n"
                       "angleConstant=0\n"
                       "[preProcess]\n"
                       "only=true\n" );

    appHarness application;
    std::string invokedName = "p4Reduce-test";
    std::string configOption = "--config";
    std::string configName = configPath.string();
    char *arguments[]{ invokedName.data(), configOption.data(), configName.data() };

    REQUIRE( application.main( 3, arguments ) == 0 );
    REQUIRE( application.m_obs.m_psfsub.empty() );
}

/// Verify p4Reduce normal mode performs a real-FITS P4 reduction, derotation, combination, header, and output
/// round-trip.
/** \ingroup p4Reduce_unit_tests */
TEST_CASE( "p4Reduce normal FITS end-to-end", "[p4Reduce][execute][normal][FITS][derotate][combine][output]" )
{
    OpenMPThreadGuard threads( 1 );
    TestDirectory directory;
    writeTarget( directory.file( "target_002.fits" ), 3, -11 );
    writeTarget( directory.file( "target_000.fits" ), 1, 0 );
    writeTarget( directory.file( "target_001.fits" ), 2, 17 );

    const auto outputDirectory = directory.file( "nested/output" );
    const auto configPath = directory.file( "normal-e2e.conf" );
    writeTextFile( configPath, p4Configuration( directory.path(), outputDirectory, "normal", false ) );

    appHarness application;
    std::string invokedName = "p4Reduce-test";
    std::string configOption = "--config";
    std::string configName = configPath.string();
    char *arguments[]{ invokedName.data(), configOption.data(), configName.data() };

    REQUIRE( application.main( 3, arguments ) == 0 );
    REQUIRE( application.m_mode == "normal" );
    REQUIRE( application.m_obs.m_derotF.m_angles == std::vector<float>{ 0, 17, -11 } );
    REQUIRE( application.m_obs.m_realizedModes == std::vector<std::vector<int>>{ { 1 } } );
    REQUIRE( application.m_obs.m_psfsub.empty() );
    REQUIRE( application.m_obs.m_psfsubValidity.empty() );

    mx::improc::eigenCube<float> output;
    mx::fits::fitsHeader<mx::verbose::vv> header;
    mx::fits::fitsFile<float, mx::verbose::vv> reader;
    REQUIRE( reader.read( output, header, ( outputDirectory / "p4-final.fits" ).string() ) == mx::error_t::noerror );
    REQUIRE( output.rows() == 31 );
    REQUIRE( output.cols() == 31 );
    REQUIRE( output.planes() == 1 );

    std::size_t finitePixels{ 0 };
    for( int row = 0; row < output.rows(); ++row )
    {
        for( int column = 0; column < output.cols(); ++column )
        {
            if( mx::math::isFinite( output.image( 0 )( row, column ) ) )
            {
                ++finitePixels;
                REQUIRE( output.image( 0 )( row, column ) == Approx( 0 ).margin( 3e-5 ) );
            }
        }
    }
    REQUIRE( finitePixels > 0 );
    REQUIRE( finitePixels < static_cast<std::size_t>( output.rows() * output.cols() ) );

    REQUIRE( header["P4 ALGORITHM"].String().starts_with( "P4-PCA" ) );
    REQUIRE( header["P4 FRAME"].String().starts_with( "detector" ) );
    REQUIRE( header["P4 IN SAMPLE"].value<char>() == 1 );
    REQUIRE( header["P4 RDI"].value<char>() == 0 );
    REQUIRE( header["P4 MODE FRACTIONS"].String().starts_with( "0.5" ) );
    REQUIRE_FALSE( header["P4 PREDICTOR COUNT"].String().empty() );
    REQUIRE( header["P4 REALIZED MODES 000"].String().starts_with( "1" ) );
    REQUIRE( header["COMBINATION METHOD"].String().starts_with( "mean" ) );
}

/// Verify p4ReduceMain reports a nested configured-run failure and treats command-line help as success.
/** \ingroup p4Reduce_unit_tests */
TEST_CASE( "p4Reduce process entry behavior", "[p4Reduce][main][help][exception]" )
{
    SECTION( "help" )
    {
        std::string invokedName = "p4Reduce-test";
        std::string helpOption = "--help";
        char *arguments[]{ invokedName.data(), helpOption.data() };

        StreamCapture standardOutput( std::cout );
        StreamCapture errorOutput( std::cerr );
        REQUIRE( p4ReduceMain( 2, arguments ) == EXIT_SUCCESS );
        const std::string help = standardOutput.str() + errorOutput.str();
        REQUIRE( help.find( "p4.modeFractions" ) != std::string::npos );
        REQUIRE( help.find( "p4.regressionFrame" ) != std::string::npos );
        REQUIRE( help.find( "grid." ) == std::string::npos );
        REQUIRE( help.find( "postprocess." ) == std::string::npos );
    }

    SECTION( "nested target-loading failure" )
    {
        TestDirectory directory;
        const auto configPath = directory.file( "empty-input.conf" );
        writeTextFile( configPath, p4Configuration( directory.path(), directory.file( "output" ), "normal", false ) );

        std::string invokedName = "p4Reduce-test";
        std::string configOption = "--config";
        std::string configName = configPath.string();
        char *arguments[]{ invokedName.data(), configOption.data(), configName.data() };

        StreamCapture errors( std::cerr );
        REQUIRE( p4ReduceMain( 3, arguments ) == EXIT_FAILURE );
        const std::string output = errors.str();
        REQUIRE( output.find( "exception(s) encountered during execution" ) != std::string::npos );
        REQUIRE( output.find( "P4 target loading" ) != std::string::npos );
        REQUIRE( output.find( "0 length" ) != std::string::npos );
        REQUIRE( output.find( "To get help try: p4Reduce-test -h" ) != std::string::npos );
    }

    // clang-format off
#ifdef __DOXY_ONLY__
    p4ReduceMain( 0, nullptr );
#endif
    // clang-format on
}

} // namespace p4Reduce_test
} // namespace unitTest
