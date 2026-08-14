/** \file ADIobservation_test.cpp
 * \brief Tests angular-differential observation ingestion, fake injection, masks, and derotation.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "HCIobservation_test_fixture.hpp"
#include "src/common/ADIDerotator.hpp"
#include "src/common/ADIobservation.hpp"

namespace unitTest
{
namespace ADIobservation_test
{

/// \cond ADIobservation_test_harness
using observationT =
    mx::improc::ADIobservation<float, mx::improc::ADIDerotator<float, mx::verbose::vv>, mx::verbose::vv>;

struct observationHarness : public observationT
{
    using cubeT = mx::improc::eigenCube<float>;

    using observationT::m_auxDataDir;
    using observationT::m_dateKeyword;
    using observationT::m_fileList;
    using observationT::m_filesDeleted;
    using observationT::m_heads;
    using observationT::m_imSize;
    using observationT::m_mask;
    using observationT::m_maskCube;
    using observationT::m_Ncols;
    using observationT::m_Nims;
    using observationT::m_Nrows;
    using observationT::m_psfsub;
    using observationT::m_RDIdateKeyword;
    using observationT::m_RDIfileList;
    using observationT::m_RDIfilesDeleted;
    using observationT::m_RDIheads;
    using observationT::m_refIms;
    using observationT::m_skipPreProcess;
    using observationT::m_tgtIms;
};

void configureDerotator( mx::improc::ADIDerotator<float, mx::verbose::vv> &derotator )
{
    derotator.angleKeyword( "ANGLE" );
    derotator.m_angleScale = 1;
    derotator.m_angleConstant = 0;
}

observationHarness::fitsHeaderT angleHeader( float angle )
{
    observationHarness::fitsHeaderT header;
    header.append<float>( "ANGLE", angle, "parallactic angle" );
    return header;
}
/// \endcond

/// Verify ADIobservation rejects invalid default and configured fake-planet methods.
/** \ingroup ADIobservation_unit_tests */
TEST_CASE( "ADIobservation configuration errors", "[ADIobservation][config]" )
{
    TestDirectory directory;

    observationHarness invalidDefault;
    mx::app::appConfigurator defaultConfig;
    invalidDefault.setupConfig( defaultConfig );
    invalidDefault.m_fakeMethod = static_cast<mx::improc::HCI::fake>( 99 );
    REQUIRE_THROWS( invalidDefault.loadConfig( defaultConfig ) );

    observationHarness invalidConfigured;
    mx::app::appConfigurator configured;
    invalidConfigured.setupConfig( configured );
    const auto path = directory.file( "invalid-fake.conf" );
    writeTextFile( path, "[fake]\nmethod=invalid\n" );
    REQUIRE( configured.readConfig( path.string() ) == 0 );
    REQUIRE_THROWS( invalidConfigured.loadConfig( configured ) );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::ADIobservation<float,
                                mx::improc::ADIDerotator<float, mx::verbose::vv>,
                                mx::verbose::vv>::setupConfig( configured );
    mx::improc::ADIobservation<float,
                                mx::improc::ADIDerotator<float, mx::verbose::vv>,
                                mx::verbose::vv>::loadConfig( configured );
#endif
    // clang-format on
}

/// Verify ADIobservation target and RDI wrappers propagate angle keywords and populate their derotators.
/** \ingroup ADIobservation_unit_tests */
TEST_CASE( "ADIobservation FITS ingestion", "[ADIobservation][readFiles][readRDIFiles]" )
{
    TestDirectory directory;
    observationHarness observation;
    configureDerotator( observation.m_derotF );
    configureDerotator( observation.m_RDIderotF );
    observation.m_dateKeyword.clear();
    observation.m_RDIdateKeyword.clear();
    observation.m_skipPreProcess = true;
    observation.m_imSize = 3;

    observationHarness::imageT image( 3, 3 );
    image.setOnes();
    auto targetHead0 = angleHeader( 10 );
    auto targetHead1 = angleHeader( 20 );
    const auto target0 = directory.file( "target0.fits" );
    const auto target1 = directory.file( "target1.fits" );
    writeFitsImage( target0, image, &targetHead0 );
    image.setConstant( 2 );
    writeFitsImage( target1, image, &targetHead1 );
    observation.m_fileList = { target0.string(), target1.string() };
    observation.readFiles();

    REQUIRE( observation.m_tgtIms.planes() == 2 );
    REQUIRE( observation.m_derotF.m_angles == std::vector<float>{ 10, 20 } );
    REQUIRE( observation.m_heads[0].count( "ANGLE" ) == 1 );

    auto referenceHead0 = angleHeader( 30 );
    auto referenceHead1 = angleHeader( 40 );
    const auto reference0 = directory.file( "reference0.fits" );
    const auto reference1 = directory.file( "reference1.fits" );
    image.setConstant( 3 );
    writeFitsImage( reference0, image, &referenceHead0 );
    image.setConstant( 4 );
    writeFitsImage( reference1, image, &referenceHead1 );
    observation.m_RDIfileList = { reference0.string(), reference1.string() };
    observation.readRDIFiles();

    REQUIRE( observation.m_refIms.planes() == 2 );
    REQUIRE( observation.m_RDIderotF.m_angles == std::vector<float>{ 30, 40 } );
    REQUIRE( observation.m_RDIheads[1].count( "ANGLE" ) == 1 );

    observationHarness unconfigured;
    unconfigured.m_fileList = { target0.string() };
    REQUIRE_THROWS( unconfigured.readFiles() );
    unconfigured.m_Nrows = 3;
    unconfigured.m_Ncols = 3;
    unconfigured.m_imSize = 3;
    unconfigured.m_RDIfileList = { reference0.string() };
    REQUIRE_THROWS( unconfigured.readRDIFiles() );

    observationHarness missingTarget;
    configureDerotator( missingTarget.m_derotF );
    missingTarget.m_dateKeyword.clear();
    missingTarget.m_fileList = { directory.file( "missing-target.fits" ).string() };
    REQUIRE_THROWS( missingTarget.readFiles() );

    observationHarness missingReference;
    configureDerotator( missingReference.m_RDIderotF );
    missingReference.m_RDIdateKeyword.clear();
    missingReference.m_Nrows = 3;
    missingReference.m_Ncols = 3;
    missingReference.m_imSize = 3;
    missingReference.m_RDIfileList = { directory.file( "missing-reference.fits" ).string() };
    REQUIRE_THROWS( missingReference.readRDIFiles() );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::ADIobservation<float,
                                mx::improc::ADIDerotator<float, mx::verbose::vv>,
                                mx::verbose::vv>::readFiles();
    mx::improc::ADIobservation<float,
                                mx::improc::ADIDerotator<float, mx::verbose::vv>,
                                mx::verbose::vv>::readRDIFiles();
#endif
    // clang-format on
}

/// Verify ADIobservation fake injection honors per-image scales and its exact single-image size contract.
/** \ingroup ADIobservation_unit_tests */
TEST_CASE( "ADIobservation fake injection", "[ADIobservation][injectFake]" )
{
    TestDirectory directory;
    observationHarness observation;
    configureDerotator( observation.m_derotF );
    observation.m_derotF.m_angles = { 0, 0 };
    observation.m_fakeMethod = mx::improc::HCI::fake::single;
    observation.m_fakeSep = { 0 };
    observation.m_fakePA = { 0 };
    observation.m_fakeContrast = { 1 };

    observationHarness::imageT fakePSF( 5, 5 );
    fakePSF.setZero();
    fakePSF( 2, 2 ) = 1;
    observation.m_fakeFileName = directory.file( "fake.fits" ).string();
    writeFitsImage( observation.m_fakeFileName, fakePSF );

    observationHarness::cubeT images( 5, 5, 2 );
    images.setZero();
    std::vector<std::string> fileList{ directory.file( "target0.fits" ).string(),
                                       directory.file( "target1.fits" ).string() };
    observation.m_fakeScaleFileName = directory.file( "scales.dat" ).string();
    writeTextFile( observation.m_fakeScaleFileName, "target0.fits 2\ntarget1.fits 3\n" );

    observation.injectFake( images, fileList, observation.m_derotF, 1, 1 );
    REQUIRE( images.image( 0 )( 2, 2 ) == Approx( 2 ) );
    REQUIRE( images.image( 1 )( 2, 2 ) == Approx( 3 ) );

    observation.m_fakePA.clear();
    REQUIRE_THROWS( observation.injectFake( images, fileList, observation.m_derotF, 1, 1 ) );
    observation.m_fakePA = { 0 };
    fileList.pop_back();
    REQUIRE_THROWS( observation.injectFake( images, fileList, observation.m_derotF, 1, 1 ) );
    fileList.push_back( directory.file( "target1.fits" ).string() );

    writeTextFile( observation.m_fakeScaleFileName, "target0.fits 2\n" );
    REQUIRE_THROWS( observation.injectFake( images, fileList, observation.m_derotF, 1, 1 ) );
    writeTextFile( observation.m_fakeScaleFileName, "target0.fits \n" );
    REQUIRE_THROWS( observation.injectFake( images, fileList, observation.m_derotF, 1, 1 ) );
    observation.m_fakeScaleFileName = directory.file( "missing-scales.dat" ).string();
    REQUIRE_THROWS( observation.injectFake( images, fileList, observation.m_derotF, 1, 1 ) );

    observationHarness::cubeT directImages( 5, 5, 1 );
    directImages.setZero();
    observationHarness::imageT smallFake( 3, 3 );
    smallFake.setZero();
    smallFake( 1, 1 ) = 2;
    observation.injectFake( smallFake, directImages, 0, 0, 0, 0, 0.5F, 2, 1, 1 );
    REQUIRE( directImages.image( 0 )( 2, 2 ) == Approx( 2 ) );

    observationHarness::imageT mismatchedFake( 3, 5 );
    mismatchedFake.setOnes();
    REQUIRE_THROWS( observation.injectFake( mismatchedFake, directImages, 0, 0, 0, 0, 1, 1, 1, 1 ) );
    mismatchedFake.resize( 7, 5 );
    mismatchedFake.setOnes();
    REQUIRE_THROWS( observation.injectFake( mismatchedFake, directImages, 0, 0, 0, 0, 1, 1, 1, 1 ) );

    observationHarness::imageT largeFake( 7, 7 );
    largeFake.setZero();
    largeFake( 3, 3 ) = 1;
    directImages.setZero();
    REQUIRE_NOTHROW( observation.injectFake( largeFake, directImages, 0, 0, 0, 0, 1, 1, 1, 1 ) );
    REQUIRE( directImages.image( 0 )( 2, 2 ) == Approx( 1 ) );

    observation.m_fakeScaleFileName.clear();
    observation.m_fakeMethod = mx::improc::HCI::fake::list;
    observation.m_fakeFileName = directory.file( "fake-list.dat" ).string();
    writeTextFile( observation.m_fakeFileName, directory.file( "fake.fits" ).string() + "\n" );
    REQUIRE_THROWS( observation.injectFake( images, fileList, observation.m_derotF, 1, 1 ) );

    writeTextFile( observation.m_fakeFileName,
                   directory.file( "fake.fits" ).string() + "\n" + directory.file( "fake.fits" ).string() + "\n" );
    images.setZero();
    REQUIRE_NOTHROW( observation.injectFake( images, fileList, observation.m_derotF, 1, 1 ) );
    REQUIRE( images.image( 0 )( 2, 2 ) == Approx( 1 ) );
    REQUIRE( images.image( 1 )( 2, 2 ) == Approx( 1 ) );

    writeTextFile( observation.m_fakeFileName,
                   directory.file( "fake.fits" ).string() + "\n" + directory.file( "missing-fake.fits" ).string() +
                       "\n" );
    REQUIRE_THROWS( observation.injectFake( images, fileList, observation.m_derotF, 1, 1 ) );
    observation.m_fakeFileName = directory.file( "missing-list.dat" ).string();
    REQUIRE_THROWS( observation.injectFake( images, fileList, observation.m_derotF, 1, 1 ) );

    observation.m_fakeMethod = mx::improc::HCI::fake::single;
    observation.m_fakeFileName = directory.file( "missing-fake.fits" ).string();
    REQUIRE_THROWS( observation.injectFake( images, fileList, observation.m_derotF, 1, 1 ) );

    observationHarness::imageT invalidFake( 3, 5 );
    invalidFake.setOnes();
    observation.m_fakeFileName = directory.file( "invalid-fake.fits" ).string();
    writeFitsImage( observation.m_fakeFileName, invalidFake );
    REQUIRE_THROWS( observation.injectFake( images, fileList, observation.m_derotF, 1, 1 ) );

    observation.m_fakeMethod = static_cast<mx::improc::HCI::fake>( 99 );
    REQUIRE_THROWS( observation.injectFake( images, fileList, observation.m_derotF, 1, 1 ) );

    observationHarness postRead;
    configureDerotator( postRead.m_derotF );
    postRead.m_heads = { angleHeader( 0 ) };
    postRead.m_fileList = { directory.file( "target0.fits" ).string() };
    postRead.m_tgtIms.resize( 5, 5, 1 );
    postRead.m_tgtIms.setZero();
    postRead.m_fakeMethod = mx::improc::HCI::fake::single;
    postRead.m_fakeFileName = directory.file( "fake.fits" ).string();
    postRead.m_fakeSep = { 0 };
    postRead.m_fakePA = { 0 };
    postRead.m_fakeContrast = { 1 };
    REQUIRE_NOTHROW( postRead.postReadFiles() );
    REQUIRE( postRead.m_tgtIms.image( 0 )( 2, 2 ) == Approx( 1 ) );
    postRead.m_fakeFileName = directory.file( "missing-fake.fits" ).string();
    REQUIRE_THROWS( postRead.postReadFiles() );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::ADIobservation<float,
                                mx::improc::ADIDerotator<float, mx::verbose::vv>,
                                mx::verbose::vv>::injectFake( images, fileList, observation.m_derotF, 1, 1 );
    mx::improc::ADIobservation<float,
                                mx::improc::ADIDerotator<float, mx::verbose::vv>,
                                mx::verbose::vv>::injectFake( fakePSF, directImages, 0, 0, 0, 0, 1, 1, 1, 1 );
    mx::improc::ADIobservation<float,
                                mx::improc::ADIDerotator<float, mx::verbose::vv>,
                                mx::verbose::vv>::postReadFiles();
#endif
    // clang-format on
}

/// Verify ADIobservation refreshes target and RDI angles after coaddition and rejects missing angle metadata.
/** \ingroup ADIobservation_unit_tests */
TEST_CASE( "ADIobservation angle hooks", "[ADIobservation][postReadFiles][postCoadd][postRDIReadFiles][postRDICoadd]" )
{
    observationHarness observation;
    configureDerotator( observation.m_derotF );
    configureDerotator( observation.m_RDIderotF );
    observation.m_fileList = { "target.fits" };
    observation.m_RDIfileList = { "reference.fits" };
    observation.m_heads = { angleHeader( 12 ) };
    observation.m_RDIheads = { angleHeader( 34 ) };

    REQUIRE_NOTHROW( observation.postCoadd() );
    REQUIRE( observation.m_derotF.m_angles == std::vector<float>{ 12 } );
    REQUIRE_NOTHROW( observation.postRDIReadFiles() );
    REQUIRE_NOTHROW( observation.postRDICoadd() );
    REQUIRE( observation.m_RDIderotF.m_angles == std::vector<float>{ 34 } );

    observation.m_heads = { observationHarness::fitsHeaderT{} };
    observation.m_RDIheads = { observationHarness::fitsHeaderT{} };
    REQUIRE_THROWS( observation.postReadFiles() );
    REQUIRE_THROWS( observation.postCoadd() );
    REQUIRE_THROWS( observation.postRDIReadFiles() );
    REQUIRE_THROWS( observation.postRDICoadd() );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::ADIobservation<float,
                                mx::improc::ADIDerotator<float, mx::verbose::vv>,
                                mx::verbose::vv>::postCoadd();
    mx::improc::ADIobservation<float,
                                mx::improc::ADIDerotator<float, mx::verbose::vv>,
                                mx::verbose::vv>::postRDIReadFiles();
    mx::improc::ADIobservation<float,
                                mx::improc::ADIDerotator<float, mx::verbose::vv>,
                                mx::verbose::vv>::postRDICoadd();
#endif
    // clang-format on
}

/// Verify ADIobservation rotates masks and PSF-subtracted planes using the configured angle sequence.
/** \ingroup ADIobservation_unit_tests */
TEST_CASE( "ADIobservation mask and image derotation", "[ADIobservation][makeMaskCube][derotate]" )
{
    OpenMPThreadGuard threads( 1 );
    TestDirectory directory;
    observationHarness observation;
    configureDerotator( observation.m_derotF );
    observation.m_derotF.m_angles = { 0, 90 };
    observation.m_Nrows = 7;
    observation.m_Ncols = 7;
    observation.m_Nims = 2;
    observation.m_auxDataDir = directory.file( "aux" ).string() + "/";
    observation.m_mask.resize( 7, 7 );
    observation.m_mask.setZero();
    observation.m_mask( 2, 3 ) = 1;

    observation.makeMaskCube();
    REQUIRE( observation.m_maskCube.planes() == 2 );
    REQUIRE( observation.m_maskCube.image( 0 ).isApprox( observation.m_mask ) );
    REQUIRE( observation.m_maskCube.image( 1 ).sum() == Approx( 1 ) );
    REQUIRE( std::filesystem::exists( directory.file( "aux/angles.dat" ) ) );
    REQUIRE( std::filesystem::exists( directory.file( "aux/maskCube.fits" ) ) );

    observation.m_psfsub.resize( 1 );
    observation.m_psfsub[0].resize( 7, 7, 2 );
    observation.m_psfsub[0].setZero();
    observation.m_psfsub[0].image( 0 )( 2, 3 ) = 1;
    observation.m_psfsub[0].image( 1 )( 2, 3 ) = 1;
    const auto unrotated = observation.m_psfsub[0].image( 0 );
    observation.derotate();
    REQUIRE( observation.m_psfsub[0].image( 0 ).isApprox( unrotated ) );
    REQUIRE_FALSE( observation.m_psfsub[0].image( 1 ).isApprox( unrotated ) );
    REQUIRE( observation.m_psfsub[0].image( 1 ).isFinite().all() );

    observation.m_mask.resize( 6, 7 );
    REQUIRE_THROWS( observation.makeMaskCube() );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::ADIobservation<float,
                                mx::improc::ADIDerotator<float, mx::verbose::vv>,
                                mx::verbose::vv>::makeMaskCube();
    mx::improc::ADIobservation<float,
                                mx::improc::ADIDerotator<float, mx::verbose::vv>,
                                mx::verbose::vv>::derotate();
#endif
    // clang-format on
}

/// Verify ADIobservation records fake-planet and post-processing provenance in FITS headers.
/** \ingroup ADIobservation_unit_tests */
TEST_CASE( "ADIobservation FITS metadata", "[ADIobservation][stdFitsHeader]" )
{
    observationHarness observation;
    observation.m_postMedSub = true;
    observation.m_fakeFileName = "fake.fits";
    observation.m_fakeScaleFileName = "scales.dat";
    observation.m_fakeSep = { 2, 4 };
    observation.m_fakePA = { 10, 20 };
    observation.m_fakeContrast = { 1e-3F, 2e-3F };

    observationHarness::fitsHeaderT header;
    observation.stdFitsHeader( &header );
    REQUIRE( header["POSTMEDS"].type() == mx::fits::fitsType<bool>() );
    REQUIRE( header["POSTMEDS"].valueGood() );
    REQUIRE( header["FAKEFILE"].String().starts_with( "fake.fits" ) );
    REQUIRE( header["FAKESCFL"].String().starts_with( "scales.dat" ) );
    REQUIRE( header["FAKESEP"].String().starts_with( "2,4" ) );
    REQUIRE( header["FAKEPA"].String().starts_with( "10,20" ) );
    REQUIRE( header["FAKECONT"].String().starts_with( "0.001,0.002" ) );
    REQUIRE_NOTHROW( observation.stdFitsHeader( nullptr ) );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::ADIobservation<float,
                                mx::improc::ADIDerotator<float, mx::verbose::vv>,
                                mx::verbose::vv>::stdFitsHeader( &header );
#endif
    // clang-format on
}

} // namespace ADIobservation_test
} // namespace unitTest
