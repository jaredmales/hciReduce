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
    using appT::m_mode;
    using appT::m_postprocessDirectory;
    using appT::m_postprocessExtension;
    using appT::m_postprocessPrefix;
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

    const auto configPath = directory.file( "postprocess.conf" );
    writeTextFile( configPath,
                   "mode=postprocess\n"
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
    REQUIRE_NOTHROW( application.checkConfig() );

    application.m_postprocessPrefix.clear();
    REQUIRE_THROWS( application.checkConfig() );
    application.m_postprocessPrefix = "saved";
    application.m_mode = "unknown";
    REQUIRE_THROWS( application.checkConfig() );
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
