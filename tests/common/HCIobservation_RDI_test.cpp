/** \file HCIobservation_RDI_test.cpp
 * \brief Tests HCIobservation reference FITS ingestion and independent RDI masks.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "HCIobservation_test_fixture.hpp"

#include <limits>

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
}

} // namespace HCIobservation_RDI_test
} // namespace unitTest
