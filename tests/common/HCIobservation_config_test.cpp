/** \file HCIobservation_config_test.cpp
 * \brief Tests HCIobservation configuration registration and loading.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "HCIobservation_test_fixture.hpp"

#include <stdexcept>

namespace unitTest
{
namespace HCIobservation_config_test
{

/// \cond HCIobservation_test_harness
namespace
{
void readObservationConfig( HCIobservationTestHarness &observation,
                            const std::filesystem::path &path,
                            const std::string &contents )
{
    writeTextFile( path, contents );
    mx::app::appConfigurator config;
    if( observation.setupConfig( config ) != 0 )
    {
        throw std::runtime_error( "could not register observation configuration" );
    }
    if( config.readConfig( path.string() ) != 0 )
    {
        throw std::runtime_error( "could not read observation configuration" );
    }
    if( observation.loadConfig( config ) != 0 )
    {
        throw std::runtime_error( "could not load observation configuration" );
    }
}
} // namespace
/// \endcond

/// Verify HCIobservation::setupConfig registers correctly typed, uniquely named configuration targets.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation configuration metadata", "[HCIobservation][config]" )
{
    HCIobservationTestHarness observation;
    mx::app::appConfigurator config;
    REQUIRE( observation.setupConfig( config ) == 0 );

    REQUIRE( config.m_targets.at( "input.dateUnit" ).clType == mx::app::argType::Required );
    REQUIRE( config.m_targets.at( "input.dateUnit" ).helpType == "float" );
    REQUIRE( config.m_targets.at( "rdi.dateUnit" ).clType == mx::app::argType::Required );
    REQUIRE( config.m_targets.at( "rdi.dateUnit" ).helpType == "float" );
    REQUIRE( config.m_targets.at( "input.dateIsISO8601" ).clType == mx::app::argType::True );
    REQUIRE( config.m_targets.at( "rdi.dateIsISO8601" ).clType == mx::app::argType::True );
    REQUIRE( config.m_targets.at( "input.thresholdOnly" ).clType == mx::app::argType::True );
    REQUIRE( config.m_targets.at( "input.thresholdOnly" ).helpType == "bool" );
    REQUIRE( config.m_targets.at( "coadd.maxAngle" ).helpType == "float" );
    REQUIRE( config.m_targets.at( "preProcess.mask" ).helpType == "bool" );
    REQUIRE( config.m_targets.at( "input.maskFile" ).helpType == "string" );
    REQUIRE( config.m_targets.at( "rdi.useInputMask" ).clType == mx::app::argType::True );
    REQUIRE( config.m_targets.at( "preProcess.medianUSM_fwhm" ).helpType == "int" );
    REQUIRE( config.m_targets.count( "preProcess.azUSM_azHalfWidth" ) == 1 );
    REQUIRE( config.m_targets.count( "preProcess.azUSM_radHalfWidth" ) == 1 );
    REQUIRE( config.m_targets.count( "preProcess.pixelTSSigma" ) == 1 );

    REQUIRE( observation.loadConfig( config ) == 0 );
    REQUIRE( observation.m_extension == ".fits" );
    REQUIRE( observation.m_qualityThreshold == 0 );
    REQUIRE( observation.m_dateKeyword == "DATE-OBS" );
    REQUIRE( observation.m_dateIsISO8601 );
    REQUIRE( observation.m_dateUnit == 1 );
    REQUIRE( observation.m_RDIextension == ".fits" );
    REQUIRE( observation.m_RDIdateKeyword == "DATE-OBS" );
    REQUIRE( observation.m_RDIdateIsISO8601 );
    REQUIRE( observation.m_RDIdateUnit == 1 );
    REQUIRE( observation.m_coaddMethod == mx::improc::HCI::coadd::none );
    REQUIRE( observation.m_coaddMaxAngle == 0 );
    REQUIRE( observation.m_preProcess_mask );
    REQUIRE( observation.m_preProcess_azUSM_maxAz == 45 );
    REQUIRE( observation.m_preProcess_meanSubMethod == mx::improc::HCI::meanSub::none );
    REQUIRE( observation.m_preProcess_pixelTSNormMethod == mx::improc::HCI::pixelTSNorm::none );
    REQUIRE( observation.m_pixelTSSigma == 3 );
    REQUIRE( observation.m_combineMethod == mx::improc::HCI::combine::mean );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::setupConfig( config );
    mx::improc::HCIobservation<float, mx::verbose::vv>::loadConfig( config );
#endif
    // clang-format on
}

/// Verify HCIobservation::loadConfig loads numeric, boolean, enum, preprocessing, and output settings.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation configuration loading", "[HCIobservation][config]" )
{
    TestDirectory directory;
    HCIobservationTestHarness observation;

    readObservationConfig( observation,
                           directory.file( "observation.conf" ),
                           "[input]\n"
                           "directory=/targets\n"
                           "prefix=tgt_\n"
                           "extension=fit\n"
                           "fileList=targets.list\n"
                           "deleteFront=1\n"
                           "deleteBack=2\n"
                           "qualityFile=quality.txt\n"
                           "qualityThreshold=0.75\n"
                           "thresholdOnly=true\n"
                           "dateKeyword=MJD\n"
                           "dateIsISO8601=false\n"
                           "dateUnit=0.000011574074\n"
                           "imSize=17\n"
                           "maskFile=mask.fits\n"
                           "[rdi]\n"
                           "directory=/references\n"
                           "prefix=ref_\n"
                           "extension=fits\n"
                           "fileList=references.list\n"
                           "deleteFront=3\n"
                           "deleteBack=4\n"
                           "qualityFile=rdi-quality.txt\n"
                           "qualityThreshold=0.5\n"
                           "dateKeyword=RDI-MJD\n"
                           "dateIsISO8601=false\n"
                           "dateUnit=0.5\n"
                           "maskFile=rdi-mask.fits\n"
                           "useInputMask=false\n"
                           "[coadd]\n"
                           "method=median\n"
                           "maxImno=3\n"
                           "maxTime=1.25\n"
                           "maxAngle=0.75\n"
                           "keywords=ANGLE,EXPTIME\n"
                           "[preProcess]\n"
                           "skip=true\n"
                           "beforeCoadd=true\n"
                           "mask=false\n"
                           "subradprof=true\n"
                           "azUSM_azHalfWidth=7\n"
                           "azUSM_maxAz=30\n"
                           "azUSM_radHalfWidth=2\n"
                           "medianUSM_fwhm=4\n"
                           "gaussUSM_fwhm=1.5\n"
                           "meanSubMethod=medianImage\n"
                           "pixelTSNormMethod=rms\n"
                           "pixelTSSigma=2.75\n"
                           "outputPrefix=/tmp/pre-\n"
                           "only=true\n"
                           "[combine]\n"
                           "method=sigmaMean\n"
                           "weightFile=weights.txt\n"
                           "sigmaThreshold=2.5\n"
                           "minGoodFract=0.8\n"
                           "[output]\n"
                           "fileName=result_\n"
                           "exactFName=true\n"
                           "directory=/output\n"
                           "outputPSFSub=true\n"
                           "psfSubPrefix=psf_\n" );

    REQUIRE( observation.m_directory == "/targets" );
    REQUIRE( observation.m_prefix == "tgt_" );
    REQUIRE( observation.m_extension == "fit" );
    REQUIRE( observation.m_fileListFile == "targets.list" );
    REQUIRE( observation.m_deleteFront == 1 );
    REQUIRE( observation.m_deleteBack == 2 );
    REQUIRE( observation.m_qualityFile == "quality.txt" );
    REQUIRE( observation.m_qualityThreshold == Approx( 0.75 ) );
    REQUIRE( observation.m_thresholdOnly );
    REQUIRE( observation.m_dateKeyword == "MJD" );
    REQUIRE_FALSE( observation.m_dateIsISO8601 );
    REQUIRE( observation.m_dateUnit == Approx( 0.000011574074 ) );
    REQUIRE( observation.m_imSize == 17 );
    REQUIRE( observation.m_maskFile == "mask.fits" );

    REQUIRE( observation.m_RDIdirectory == "/references" );
    REQUIRE( observation.m_RDIprefix == "ref_" );
    REQUIRE( observation.m_RDIextension == "fits" );
    REQUIRE( observation.m_RDIfileListFile == "references.list" );
    REQUIRE( observation.m_RDIdeleteFront == 3 );
    REQUIRE( observation.m_RDIdeleteBack == 4 );
    REQUIRE( observation.m_RDIqualityFile == "rdi-quality.txt" );
    REQUIRE( observation.m_RDIqualityThreshold == Approx( 0.5 ) );
    REQUIRE( observation.m_RDIdateKeyword == "RDI-MJD" );
    REQUIRE_FALSE( observation.m_RDIdateIsISO8601 );
    REQUIRE( observation.m_RDIdateUnit == Approx( 0.5 ) );
    REQUIRE( observation.m_RDImaskFile == "rdi-mask.fits" );
    REQUIRE_FALSE( observation.m_RDImaskUseInput );

    REQUIRE( observation.m_coaddMethod == mx::improc::HCI::coadd::median );
    REQUIRE( observation.m_coaddMaxImno == 3 );
    REQUIRE( observation.m_coaddMaxTime == Approx( 1.25 ) );
    REQUIRE( observation.m_coaddMaxAngle == Approx( 0.75 ) );
    REQUIRE( observation.m_coaddKeywords == std::vector<std::string>{ "ANGLE", "EXPTIME" } );
    REQUIRE( observation.m_skipPreProcess );
    REQUIRE( observation.m_preProcess_beforeCoadd );
    REQUIRE_FALSE( observation.m_preProcess_mask );
    REQUIRE( observation.m_preProcess_subradprof );
    REQUIRE( observation.m_preProcess_azUSM_azHalfWidth == Approx( 7 ) );
    REQUIRE( observation.m_preProcess_azUSM_maxAz == Approx( 30 ) );
    REQUIRE( observation.m_preProcess_azUSM_radHalfWidth == Approx( 2 ) );
    REQUIRE( observation.m_preProcess_medianUSM_fwhm == 4 );
    REQUIRE( observation.m_preProcess_gaussUSM_fwhm == Approx( 1.5 ) );
    REQUIRE( observation.m_preProcess_meanSubMethod == mx::improc::HCI::meanSub::medianImage );
    REQUIRE( observation.m_preProcess_pixelTSNormMethod == mx::improc::HCI::pixelTSNorm::rms );
    REQUIRE( observation.m_pixelTSSigma == Approx( 2.75 ) );
    REQUIRE( observation.m_preProcess_outputPrefix == "/tmp/pre-" );
    REQUIRE( observation.m_preProcess_only );

    REQUIRE( observation.m_combineMethod == mx::improc::HCI::combine::sigmaMean );
    REQUIRE( observation.m_weightFile == "weights.txt" );
    REQUIRE( observation.m_sigmaThreshold == Approx( 2.5 ) );
    REQUIRE( observation.m_minGoodFract == Approx( 0.8 ) );
    REQUIRE( observation.m_finimName == "result_" );
    REQUIRE( observation.m_exactFinimName );
    REQUIRE( observation.m_outputDir == "/output" );
    REQUIRE( observation.m_doOutputPSFSub );
    REQUIRE( observation.m_PSFSubPrefix == "psf_" );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::setupConfig( config );
    mx::improc::HCIobservation<float, mx::verbose::vv>::loadConfig( config );
#endif
    // clang-format on
}

/// Verify HCIobservation::loadConfig rejects unsupported configuration enum strings.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation invalid configuration methods", "[HCIobservation][config]" )
{
    TestDirectory directory;

    HCIobservationTestHarness observation;
    REQUIRE_THROWS(
        readObservationConfig( observation, directory.file( "bad-coadd.conf" ), "[coadd]\nmethod=not-a-method\n" ) );

    HCIobservationTestHarness observation2;
    REQUIRE_THROWS( readObservationConfig( observation2,
                                           directory.file( "bad-mean.conf" ),
                                           "[preProcess]\nmeanSubMethod=imageMean\n" ) );

    HCIobservationTestHarness observation3;
    REQUIRE_THROWS( readObservationConfig( observation3,
                                           directory.file( "bad-norm.conf" ),
                                           "[preProcess]\npixelTSNormMethod=not-a-method\n" ) );

    HCIobservationTestHarness observation4;
    REQUIRE_THROWS( readObservationConfig( observation4,
                                           directory.file( "bad-combine.conf" ),
                                           "[combine]\nmethod=not-a-method\n" ) );

    {
        HCIobservationTestHarness invalidCurrentMean;
        mx::app::appConfigurator config;
        invalidCurrentMean.setupConfig( config );
        invalidCurrentMean.m_preProcess_meanSubMethod = static_cast<mx::improc::HCI::meanSub>( 99 );
        REQUIRE_THROWS( invalidCurrentMean.loadConfig( config ) );
    }

    {
        HCIobservationTestHarness invalidCurrentNorm;
        mx::app::appConfigurator config;
        invalidCurrentNorm.setupConfig( config );
        invalidCurrentNorm.m_preProcess_pixelTSNormMethod = static_cast<mx::improc::HCI::pixelTSNorm>( 99 );
        REQUIRE_THROWS( invalidCurrentNorm.loadConfig( config ) );
    }

    {
        HCIobservationTestHarness invalidCurrentCombination;
        mx::app::appConfigurator config;
        invalidCurrentCombination.setupConfig( config );
        invalidCurrentCombination.m_combineMethod = static_cast<mx::improc::HCI::combine>( 99 );
        REQUIRE_THROWS( invalidCurrentCombination.loadConfig( config ) );
    }

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::loadConfig( config );
#endif
    // clang-format on
}

} // namespace HCIobservation_config_test
} // namespace unitTest
