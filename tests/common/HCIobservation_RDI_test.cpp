/** \file HCIobservation_RDI_test.cpp
 * \brief Tests HCIobservation reference FITS ingestion and independent RDI masks.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "HCIobservation_test_fixture.hpp"

#include <limits>
#include <stdexcept>

namespace unitTest
{
namespace HCIobservation_RDI_test
{

/// \cond HCIobservation_test_harness
namespace
{
HCIobservationTestHarness::imageT indexedRDIImage( int rows, int cols, float offset = 0 )
{
    HCIobservationTestHarness::imageT image( rows, cols );
    for( int row = 0; row < rows; ++row )
    {
        for( int col = 0; col < cols; ++col )
        {
            image( row, col ) = offset + 10 * row + col;
        }
    }
    return image;
}

void configureTargetSize( HCIobservationTestHarness &observation, int size = 3 )
{
    observation.m_imSize = size;
    observation.m_Nrows = size;
    observation.m_Ncols = size;
    observation.m_Npix = size * size;
}

struct throwingRDIHookHarness : public HCIobservationTestHarness
{
    bool m_throwPostRead{ false };
    bool m_throwPostCoadd{ false };

    void postRDIReadFiles() override
    {
        if( m_throwPostRead )
        {
            throw std::runtime_error( "RDI post-read test failure" );
        }
        HCIobservationTestHarness::postRDIReadFiles();
    }

    void postRDICoadd() override
    {
        if( m_throwPostCoadd )
        {
            throw std::runtime_error( "RDI post-coadd test failure" );
        }
        HCIobservationTestHarness::postRDICoadd();
    }
};
} // namespace
/// \endcond

/// Verify HCIobservation::readRDIFiles crops references and uses independent numeric date/header configuration.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation RDI FITS ingestion", "[HCIobservation][readRDIFiles][FITS]" )
{
    TestDirectory directory;
    HCIobservationTestHarness observation;
    configureTargetSize( observation );

    auto image = indexedRDIImage( 5, 6 );
    image( 2, 2 ) = std::numeric_limits<float>::quiet_NaN();
    HCIobservationTestHarness::fitsHeaderT header;
    header.append<float>( "RSTAMP", 4, "reference time" );
    header.append<int>( "KEEP", 31, "selected value" );
    const auto path = directory.file( "reference.fits" );
    writeFitsImage( path, image, &header );

    observation.m_RDIfileList = { path.string() };
    observation.m_dateKeyword = "WRONG";
    observation.m_dateIsISO8601 = true;
    observation.m_RDIdateKeyword = "RSTAMP";
    observation.m_RDIdateIsISO8601 = false;
    observation.m_RDIdateUnit = 2;
    observation.m_RDIkeywords = { "KEEP" };
    observation.m_skipPreProcess = true;
    observation.readRDIFiles();

    REQUIRE( observation.m_RDIfilesRead );
    REQUIRE( observation.m_RDIfilesDeleted );
    REQUIRE( observation.m_refIms.rows() == 3 );
    REQUIRE( observation.m_refIms.cols() == 3 );
    REQUIRE( observation.m_refIms.planes() == 1 );
    REQUIRE( observation.m_refIms.image( 0 )( 0, 0 ) == Approx( 11 ) );
    REQUIRE( observation.m_refIms.image( 0 )( 1, 1 ) == 0 );
    REQUIRE( observation.m_RDIheads[0]["KEEP"].Int() == 31 );
    REQUIRE( observation.m_RDIimageMJD == std::vector<double>{ 8 } );
    REQUIRE( observation.m_hookEvents == std::vector<std::string>{ "rdi-read" } );

    observation.m_RDIdateKeyword.clear();
    observation.m_RDIfilesDeleted = false;
    observation.readRDIFiles();
    REQUIRE( observation.m_RDIimageMJD.empty() );
    REQUIRE( observation.m_hookEvents == std::vector<std::string>{ "rdi-read", "rdi-read" } );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::readRDIFiles();
#endif
    // clang-format on
}

/// Verify HCIobservation::readRDIFiles enforces target size, reference dimensions, and deletion bounds.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation RDI validation", "[HCIobservation][readRDIFiles][validation]" )
{
    TestDirectory directory;
    HCIobservationTestHarness observation;
    observation.m_RDIfileList = { directory.file( "reference.fits" ).string() };
    REQUIRE_THROWS( observation.readRDIFiles() );

    configureTargetSize( observation );
    auto tooSmall = indexedRDIImage( 2, 5 );
    writeFitsImage( observation.m_RDIfileList[0], tooSmall );
    observation.m_RDIdateKeyword.clear();
    observation.m_skipPreProcess = true;
    REQUIRE_THROWS( observation.readRDIFiles() );
    REQUIRE_FALSE( observation.m_RDIfilesRead );

    const std::vector<std::string> original{ "first.fits", "second.fits" };
    observation.m_RDIfileList = original;
    observation.m_RDIdeleteFront = 1;
    observation.m_RDIdeleteBack = 2;
    observation.m_RDIfilesDeleted = false;
    REQUIRE_THROWS( observation.readRDIFiles() );
    REQUIRE( observation.m_RDIfileList == original );
    REQUIRE_FALSE( observation.m_RDIfilesDeleted );

    observation.m_RDIdeleteFront = -1;
    observation.m_RDIdeleteBack = 0;
    REQUIRE_THROWS( observation.readRDIFiles() );
    REQUIRE( observation.m_RDIfileList == original );

    HCIobservationTestHarness empty;
    configureTargetSize( empty );
    REQUIRE_THROWS( empty.readRDIFiles() );

    HCIobservationTestHarness emptiedByDeletion;
    configureTargetSize( emptiedByDeletion );
    emptiedByDeletion.m_RDIfileList = { "only.fits" };
    emptiedByDeletion.m_RDIdeleteFront = 1;
    REQUIRE_THROWS( emptiedByDeletion.readRDIFiles() );
    REQUIRE( emptiedByDeletion.m_RDIfileList.empty() );
}

/// Verify HCIobservation::readRDIFiles applies valid deletion, quality, ISO-date, and centered-crop settings.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation RDI selection and ISO date", "[HCIobservation][readRDIFiles][quality][date]" )
{
    TestDirectory directory;
    const auto firstPath = directory.file( "first.fits" );
    const auto middlePath = directory.file( "middle.fits" );
    const auto lastPath = directory.file( "last.fits" );
    HCIobservationTestHarness::fitsHeaderT header;
    header.append<std::string>( "DATE-OBS", "2024-02-03T04:05:06", "observation time" );
    writeFitsImage( firstPath, indexedRDIImage( 5, 3 ), &header );
    writeFitsImage( middlePath, indexedRDIImage( 5, 3, 100 ), &header );
    writeFitsImage( lastPath, indexedRDIImage( 5, 3, 200 ), &header );

    const auto qualityPath = directory.file( "quality.txt" );
    writeTextFile( qualityPath, "middle.fits 2\n" );

    HCIobservationTestHarness observation;
    configureTargetSize( observation );
    observation.m_RDIfileList = { firstPath.string(), middlePath.string(), lastPath.string() };
    observation.m_RDIdeleteFront = 1;
    observation.m_RDIdeleteBack = 1;
    observation.m_RDIqualityFile = qualityPath.string();
    observation.m_RDIqualityThreshold = 2;
    observation.m_RDIdateKeyword = "DATE-OBS";
    observation.m_RDIdateIsISO8601 = true;
    observation.m_skipPreProcess = true;

    observation.readRDIFiles();

    REQUIRE( observation.m_RDIfileList == std::vector<std::string>{ middlePath.string() } );
    REQUIRE( observation.m_refIms.rows() == 3 );
    REQUIRE( observation.m_refIms.cols() == 3 );
    REQUIRE( observation.m_refIms.image( 0 )( 0, 0 ) == Approx( 110 ) );
    REQUIRE( observation.m_RDIimageMJD.size() == 1 );
    REQUIRE( observation.m_RDIimageMJD[0] == Approx( mx::sys::ISO8601date2mjd( "2024-02-03T04:05:06" ) ) );
    REQUIRE( observation.m_hookEvents == std::vector<std::string>{ "rdi-read" } );
}

/// Verify HCIobservation::readRDIFiles propagates first and later reference FITS failures.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation RDI FITS failures", "[HCIobservation][readRDIFiles][FITS][validation]" )
{
    TestDirectory directory;
    HCIobservationTestHarness observation;
    configureTargetSize( observation, 2 );
    observation.m_RDIdateKeyword.clear();
    observation.m_skipPreProcess = true;
    observation.m_RDIfileList = { directory.file( "missing-first.fits" ).string() };

    REQUIRE_THROWS( observation.readRDIFiles() );
    REQUIRE_FALSE( observation.m_RDIfilesRead );

    const auto validPath = directory.file( "valid.fits" );
    auto valid = indexedRDIImage( 2, 2 );
    writeFitsImage( validPath, valid );
    observation.m_RDIfileList = { validPath.string(), directory.file( "missing-later.fits" ).string() };
    observation.m_RDIfilesDeleted = false;

    REQUIRE_THROWS( observation.readRDIFiles() );
    REQUIRE_FALSE( observation.m_RDIfilesRead );
}

/// Verify HCIobservation::readRDIFiles adds operation context to failures from each delegated processing stage.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation RDI delegated failures", "[HCIobservation][readRDIFiles][validation][delegation]" )
{
    TestDirectory directory;
    const auto imagePath = directory.file( "reference.fits" );
    writeFitsImage( imagePath, indexedRDIImage( 3, 3 ) );

    HCIobservationTestHarness thresholdFailure;
    configureTargetSize( thresholdFailure );
    thresholdFailure.m_RDIfileList = { imagePath.string() };
    thresholdFailure.m_RDIqualityFile = directory.file( "missing-quality.txt" ).string();
    thresholdFailure.m_RDIqualityThreshold = 1;
    REQUIRE_THROWS( thresholdFailure.readRDIFiles() );

    const auto qualityPath = directory.file( "quality.txt" );
    writeTextFile( qualityPath, "reference.fits 1\n" );
    HCIobservationTestHarness thresholdEmpty;
    configureTargetSize( thresholdEmpty );
    thresholdEmpty.m_RDIfileList = { imagePath.string() };
    thresholdEmpty.m_RDIqualityFile = qualityPath.string();
    thresholdEmpty.m_RDIqualityThreshold = 2;
    REQUIRE_THROWS( thresholdEmpty.readRDIFiles() );

    throwingRDIHookHarness postReadFailure;
    configureTargetSize( postReadFailure );
    postReadFailure.m_RDIfileList = { imagePath.string() };
    postReadFailure.m_RDIdateKeyword.clear();
    postReadFailure.m_skipPreProcess = true;
    postReadFailure.m_throwPostRead = true;
    REQUIRE_THROWS( postReadFailure.readRDIFiles() );

    HCIobservationTestHarness beforeCoaddFailure;
    configureTargetSize( beforeCoaddFailure );
    beforeCoaddFailure.m_RDIfileList = { imagePath.string() };
    beforeCoaddFailure.m_RDIdateKeyword.clear();
    beforeCoaddFailure.m_preProcess_beforeCoadd = true;
    beforeCoaddFailure.m_preProcess_meanSubMethod = static_cast<mx::improc::HCI::meanSub>( 99 );
    REQUIRE_THROWS( beforeCoaddFailure.readRDIFiles() );

    HCIobservationTestHarness coaddFailure;
    configureTargetSize( coaddFailure );
    coaddFailure.m_RDIfileList = { imagePath.string() };
    coaddFailure.m_RDIdateKeyword.clear();
    coaddFailure.m_coaddMethod = static_cast<mx::improc::HCI::coadd>( 99 );
    coaddFailure.m_coaddMaxImno = 1;
    REQUIRE_THROWS( coaddFailure.readRDIFiles() );

    throwingRDIHookHarness postCoaddFailure;
    configureTargetSize( postCoaddFailure );
    postCoaddFailure.m_RDIfileList = { imagePath.string() };
    postCoaddFailure.m_RDIdateKeyword.clear();
    postCoaddFailure.m_coaddMethod = mx::improc::HCI::coadd::mean;
    postCoaddFailure.m_coaddMaxImno = 1;
    postCoaddFailure.m_throwPostCoadd = true;
    REQUIRE_THROWS( postCoaddFailure.readRDIFiles() );

    HCIobservationTestHarness afterCoaddFailure;
    configureTargetSize( afterCoaddFailure );
    afterCoaddFailure.m_RDIfileList = { imagePath.string() };
    afterCoaddFailure.m_RDIdateKeyword.clear();
    afterCoaddFailure.m_preProcess_meanSubMethod = static_cast<mx::improc::HCI::meanSub>( 99 );
    REQUIRE_THROWS( afterCoaddFailure.readRDIFiles() );
}

/// Verify HCIobservation::readRDIFiles applies and rebuilds a distinct RDI mask through coaddition.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation distinct RDI mask", "[HCIobservation][readRDIFiles][mask][coadd]" )
{
    TestDirectory directory;
    HCIobservationTestHarness observation;
    configureTargetSize( observation );

    HCIobservationTestHarness::imageT image( 3, 3 );
    image.setConstant( 2 );
    const auto firstPath = directory.file( "first.fits" );
    const auto secondPath = directory.file( "second.fits" );
    writeFitsImage( firstPath, image );
    writeFitsImage( secondPath, image );

    HCIobservationTestHarness::imageT mask( 3, 3 );
    mask.setOnes();
    mask( 1, 1 ) = 0;
    const auto maskPath = directory.file( "rdi-mask.fits" );
    writeFitsImage( maskPath, mask );

    observation.m_RDIfileList = { firstPath.string(), secondPath.string() };
    observation.m_RDIdateKeyword.clear();
    observation.m_RDImaskFile = maskPath.string();
    observation.m_maskFile = "target-mask-is-different.fits";
    observation.m_mask.resize( 3, 3 );
    observation.m_mask.setZero();
    observation.m_preProcess_mask = true;
    observation.m_coaddMethod = mx::improc::HCI::coadd::mean;
    observation.m_coaddMaxImno = 2;
    observation.readRDIFiles();

    REQUIRE( observation.m_refIms.planes() == 1 );
    REQUIRE( observation.m_refIms.image( 0 )( 0, 0 ) == Approx( 2 ) );
    REQUIRE( observation.m_refIms.image( 0 )( 1, 1 ) == 0 );
    REQUIRE( observation.m_RDImask.isApprox( mask ) );
    REQUIRE( observation.m_RDImaskCube.planes() == 1 );
    REQUIRE( observation.m_RDImaskCube.image( 0 ).isApprox( mask ) );
    REQUIRE( observation.m_hookEvents == std::vector<std::string>{ "rdi-read", "rdi-coadd" } );
}

/// Verify HCIobservation::readRDIFiles can reuse the input mask or process references with no mask independently.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation input and absent RDI masks", "[HCIobservation][readRDIFiles][mask]" )
{
    TestDirectory directory;
    HCIobservationTestHarness::imageT image( 3, 3 );
    image.setConstant( 3 );
    const auto imagePath = directory.file( "reference.fits" );
    writeFitsImage( imagePath, image );

    HCIobservationTestHarness reused;
    configureTargetSize( reused );
    reused.m_RDIfileList = { imagePath.string() };
    reused.m_RDIdateKeyword.clear();
    reused.m_maskFile = "configured-input-mask.fits";
    reused.m_mask.resize( 3, 3 );
    reused.m_mask.setOnes();
    reused.m_mask( 0, 0 ) = 0;
    reused.m_RDImaskUseInput = true;
    reused.m_preProcess_mask = true;
    reused.readRDIFiles();
    REQUIRE( reused.m_refIms.image( 0 )( 0, 0 ) == 0 );
    REQUIRE( reused.m_refIms.image( 0 )( 1, 1 ) == Approx( 3 ) );
    REQUIRE( reused.m_RDImask.isApprox( reused.m_mask ) );

    HCIobservationTestHarness unmasked;
    configureTargetSize( unmasked );
    unmasked.m_RDIfileList = { imagePath.string() };
    unmasked.m_RDIdateKeyword.clear();
    unmasked.m_maskFile = "target-only-mask.fits";
    unmasked.m_mask.resize( 3, 3 );
    unmasked.m_mask.setZero();
    unmasked.m_preProcess_mask = true;
    unmasked.readRDIFiles();
    REQUIRE( unmasked.m_refIms.image( 0 ).isApprox( image ) );
    REQUIRE( unmasked.m_RDImask.size() == 0 );
    REQUIRE( unmasked.m_RDImaskCube.planes() == 0 );

    HCIobservationTestHarness missingInputMask;
    configureTargetSize( missingInputMask );
    missingInputMask.m_RDIfileList = { imagePath.string() };
    missingInputMask.m_RDIdateKeyword.clear();
    missingInputMask.m_RDImaskUseInput = true;
    missingInputMask.m_skipPreProcess = true;
    REQUIRE_THROWS( missingInputMask.readRDIFiles() );
    REQUIRE_FALSE( missingInputMask.m_RDIfilesRead );
}

} // namespace HCIobservation_RDI_test
} // namespace unitTest
