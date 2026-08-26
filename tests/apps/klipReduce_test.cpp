/** \file klipReduce_test.cpp
 * \brief Tests klipReduce saved-product configuration and dispatch.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "../common/HCIobservation_test_fixture.hpp"

#define HCIREDUCE_KLIPREDUCE_NO_MAIN
#include "src/apps/klipReduce.cpp"

namespace unitTest
{
namespace klipReduce_test
{

/// \cond klipReduce_test_harness
using appT = klipReduce<float, double, mx::verbose::vv>;

struct appHarness : public appT
{
    using appT::checkConfig;
    using appT::config;
    using appT::execute;
    using appT::loadConfig;
    using appT::m_configPathCL;
    using appT::m_mode;
    using appT::m_obs;
    using appT::m_postprocessDirectory;
    using appT::m_postprocessExtension;
    using appT::m_postprocessPrefix;
    using appT::m_showTiming;
    using appT::setConfigPathCL;
    using appT::setupConfig;
};

void writeSavedProduct(
    const std::filesystem::path &path, int reduction, int imageIndex, float value, const std::string &modeList )
{
    mx::improc::eigenImage<float> image( 2, 2 );
    image.setConstant( value );
    mx::fits::fitsHeader<mx::verbose::vv> header;
    header.append<int>( "REDUCTION", reduction, "reduction index" );
    header.append<int>( "IMAGE", imageIndex, "input image index" );
    header.append<std::string>( "NMODES", modeList, "KLIP mode counts" );
    header.append<float>( "ANGLE", static_cast<float>( imageIndex * 10 ), "parallactic angle" );

    mx::fits::fitsFile<float, mx::verbose::vv> writer;
    if( writer.write( path.string(), image, header ) != mx::error_t::noerror )
    {
        throw std::runtime_error( "could not write saved KLIP fixture" );
    }
}
/// \endcond

/// Verify klipReduce extracts the command-line configuration target into the application path.
/** \ingroup klipReduce_unit_tests */
TEST_CASE( "klipReduce command-line configuration path", "[klipReduce][config][command-line]" )
{
    appHarness application;
    mx::app::configTarget target;
    target.name = "config";
    target.set = true;
    target.values = { "/tmp/klipReduce-test.conf" };
    application.config.m_targets[target.name] = target;

    application.setConfigPathCL();

    REQUIRE( application.m_configPathCL == "/tmp/klipReduce-test.conf" );
}

/// Verify klipReduce reports unused targets and rejects unknown configuration entries and positional arguments.
/** \ingroup klipReduce_unit_tests */
TEST_CASE( "klipReduce configuration diagnostics", "[klipReduce][config][diagnostics]" )
{
    SECTION( "source locations enabled" )
    {
        appHarness application;
        application.setupConfig();

        mx::app::configTarget unusedTarget;
        unusedTarget.name = "unused.target";
        application.config.m_targets[unusedTarget.name] = unusedTarget;

        mx::app::configTarget unknownTarget;
        unknownTarget.name = "unknown.target";
        unknownTarget.sources = { "unknown.conf" };
        application.config.m_unusedConfigs[unknownTarget.name] = unknownTarget;
        application.config.nonOptions = { "unexpected" };

        REQUIRE_THROWS( application.loadConfig() );
        REQUIRE( application.config.m_targets.at( "unused.target" ).used == false );
    }

    SECTION( "source locations disabled" )
    {
        appHarness application;
        application.setupConfig();
        application.config.m_sources = false;

        mx::app::configTarget unknownTarget;
        unknownTarget.name = "unknown.target";
        application.config.m_unusedConfigs[unknownTarget.name] = unknownTarget;

        REQUIRE_THROWS( application.loadConfig() );
    }

    SECTION( "dotted config-file key" )
    {
        TestDirectory directory;
        const auto configPath = directory.file( "dotted.conf" );
        writeTextFile( configPath, "klip.Nmodes=1,2\n" );

        appHarness application;
        application.setupConfig();
        REQUIRE( application.config.readConfig( configPath.string() ) == 0 );
        REQUIRE_FALSE( application.config.m_unusedConfigs.empty() );
        REQUIRE_THROWS( application.loadConfig() );
    }
}

/// Verify klipReduce reports every missing or inconsistent non-postprocess geometry requirement.
/** \ingroup klipReduce_unit_tests */
TEST_CASE( "klipReduce reduction configuration validation", "[klipReduce][config][validation]" )
{
    appHarness application;
    application.m_mode = "normal";

    REQUIRE_NOTHROW( application.checkConfig() );

    application.m_obs.m_Nmodes = { 1 };
    application.m_obs.m_minRadius = { 1, 2 };
    application.m_obs.m_maxRadius = { 3 };
    application.m_obs.m_minAngle = { 0, 90 };
    application.m_obs.m_maxAngle = { 180 };

    REQUIRE_NOTHROW( application.checkConfig() );
}

/// Verify klipReduce executes a complete basic-mode dispatch through the documented preprocess-only stopping point.
/** \ingroup klipReduce_unit_tests */
TEST_CASE( "klipReduce basic preprocess-only dispatch", "[klipReduce][execute][basic][preprocess]" )
{
    TestDirectory directory;
    HCIobservationTestHarness::imageT image( 5, 5 );
    image.setConstant( 1 );
    HCIobservationTestHarness::fitsHeaderT header;
    header.append<float>( "ANGLE", 0, "parallactic angle" );
    writeFitsImage( directory.file( "target.fits" ), image, &header );

    const auto configPath = directory.file( "basic.conf" );
    writeTextFile( configPath,
                   "mode=basic\n"
                   "[input]\n"
                   "directory=" +
                       directory.path().string() +
                       "\n"
                       "prefix=target\n"
                       "extension=.fits\n"
                       "dateKeyword=\n"
                       "angleKeyword=ANGLE\n"
                       "angleScale=1\n"
                       "[preProcess]\n"
                       "only=true\n" );

    appHarness application;
    std::string invokedName = "klipReduce-test";
    std::string configOption = "--config";
    std::string configName = configPath.string();
    char *arguments[]{ invokedName.data(), configOption.data(), configName.data() };

    REQUIRE( application.main( 3, arguments ) == 0 );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv> doxygenObservation;
    doxygenObservation.preprocessingOnly();
#endif
    // clang-format on
}

/// Verify klipReduce registers, loads, and validates saved-product configuration.
/** \ingroup klipReduce_unit_tests */
TEST_CASE( "klipReduce postprocess configuration", "[klipReduce][config][postprocess]" )
{
    TestDirectory directory;
    appHarness application;
    application.setupConfig();

    REQUIRE( application.config.m_targets.at( "postprocess.directory" ).clType == mx::app::argType::Required );
    REQUIRE( application.config.m_targets.at( "postprocess.prefix" ).helpType == "string" );
    REQUIRE( application.config.m_targets.at( "postprocess.extension" ).helpType == "string" );
    REQUIRE( application.config.m_targets.at( "showTiming" ).helpType == "bool" );
    REQUIRE( application.config.m_targets.at( "mode" ).helpExplanation.find( "grid" ) == std::string::npos );
    bool hasGridTarget = false;
    for( const auto &[name, target] : application.config.m_targets )
    {
        static_cast<void>( target );
        hasGridTarget = hasGridTarget || name.starts_with( "grid." );
    }
    REQUIRE_FALSE( hasGridTarget );

    const auto configPath = directory.file( "postprocess.conf" );
    writeTextFile( configPath,
                   "mode=postprocess\n"
                   "showTiming=true\n"
                   "[postprocess]\n"
                   "directory=" +
                       directory.path().string() +
                       "\n"
                       "prefix=saved\n"
                       "extension=.fit\n" );
    REQUIRE( application.config.readConfig( configPath.string() ) == 0 );
    application.loadConfig();

    REQUIRE( application.m_mode == "postprocess" );
    REQUIRE( application.m_postprocessDirectory == directory.path().string() );
    REQUIRE( application.m_postprocessPrefix == "saved" );
    REQUIRE( application.m_postprocessExtension == ".fit" );
    REQUIRE( application.m_showTiming );
    REQUIRE_NOTHROW( application.checkConfig() );

    application.m_postprocessPrefix.clear();
    REQUIRE_THROWS( application.checkConfig() );
    application.m_postprocessPrefix = "saved";
    application.m_mode = "grid";
    REQUIRE_THROWS( application.checkConfig() );
    application.m_mode = "unknown";
    REQUIRE_THROWS( application.checkConfig() );
}

/// Verify klipReduce honors an attached false value for an application boolean option.
/** This exercises klipReduce::loadConfig() through the real command-line parser.
 * \ingroup klipReduce_unit_tests
 */
TEST_CASE( "klipReduce explicit false command-line option", "[klipReduce][config][bool]" )
{
    appHarness application;
    application.setupConfig();

    std::string invoked = "klipReduce";
    std::string showTiming = "--showTiming=false";
    char *arguments[]{ invoked.data(), showTiming.data() };
    application.config.parseCommandLine( 2, arguments );
    REQUIRE_NOTHROW( application.loadConfig() );
    REQUIRE_FALSE( application.m_showTiming );

    // clang-format off
#ifdef __DOXY_ONLY__
    klipReduce doxygenApplication;
    doxygenApplication.loadConfig();
#endif
    // clang-format on
}

/// Verify klipReduce dispatches postprocess mode through KLIPreduction::processPSFSub and writes the configured result.
/** \ingroup klipReduce_unit_tests */
TEST_CASE( "klipReduce postprocess dispatch", "[klipReduce][execute][postprocess]" )
{
    TestDirectory directory;
    const auto savedDirectory = directory.file( "saved" );
    const auto outputDirectory = directory.file( "output" );
    std::filesystem::create_directories( savedDirectory );

    writeSavedProduct( savedDirectory / "product_000_000.fit", 0, 0, 1, "1,3" );
    writeSavedProduct( savedDirectory / "product_000_001.fit", 0, 1, 3, "1,3" );
    writeSavedProduct( savedDirectory / "product_001_000.fit", 1, 0, 3, "1,3" );
    writeSavedProduct( savedDirectory / "product_001_001.fit", 1, 1, 5, "1,3" );

    appHarness application;
    application.setupConfig();
    const auto configPath = directory.file( "dispatch.conf" );
    writeTextFile( configPath,
                   "mode=postprocess\n"
                   "[postprocess]\n"
                   "directory=" +
                       savedDirectory.string() +
                       "\n"
                       "prefix=product\n"
                       "extension=.fit\n"
                       "[input]\n"
                       "angleKeyword=ANGLE\n"
                       "angleScale=1\n"
                       "[combine]\n"
                       "method=mean\n"
                       "noDerotate=true\n"
                       "[output]\n"
                       "directory=" +
                       outputDirectory.string() +
                       "\n"
                       "fileName=resumed.fits\n"
                       "exactFName=true\n" );
    std::string invokedName = "klipReduce-test";
    std::string configOption = "--config";
    std::string configName = configPath.string();
    char *arguments[]{ invokedName.data(), configOption.data(), configName.data() };
    REQUIRE( application.main( 3, arguments ) == 0 );

    mx::improc::eigenCube<float> output;
    mx::fits::fitsHeader<mx::verbose::vv> header;
    mx::fits::fitsFile<float, mx::verbose::vv> reader;
    REQUIRE( reader.read( output, header, ( outputDirectory / "resumed.fits" ).string() ) == mx::error_t::noerror );
    REQUIRE( output.rows() == 2 );
    REQUIRE( output.cols() == 2 );
    REQUIRE( output.planes() == 2 );
    REQUIRE( output.image( 0 ).isConstant( 2 ) );
    REQUIRE( output.image( 1 ).isConstant( 4 ) );
    REQUIRE( header["NMODES"].String().starts_with( "1,3" ) );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::KLIPreduction<float,
                               mx::improc::ADIDerotator<float, mx::verbose::vv>,
                               double,
                               mx::verbose::vv>::processPSFSub( "", "" );
#endif
    // clang-format on
}

} // namespace klipReduce_test
} // namespace unitTest
