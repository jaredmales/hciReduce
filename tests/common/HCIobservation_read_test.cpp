/** \file HCIobservation_read_test.cpp
 * \brief Tests HCIobservation target FITS ingestion and mask construction.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "HCIobservation_test_fixture.hpp"

#include <cmath>
#include <limits>

namespace unitTest
{
namespace HCIobservation_read_test
{

/// \cond HCIobservation_test_harness
namespace
{
HCIobservationTestHarness::imageT indexedImage( int rows, int cols, float offset = 0 )
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
} // namespace
/// \endcond

/// Verify HCIobservation::readFiles preserves order, crops rectangular images, propagates headers, and removes NaNs.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation target FITS ingestion", "[HCIobservation][readFiles][FITS]" )
{
    TestDirectory directory;
    HCIobservationTestHarness observation;

    auto first = indexedImage( 5, 6 );
    auto second = indexedImage( 5, 6, 100 );
    first( 2, 2 ) = std::numeric_limits<float>::quiet_NaN();

    HCIobservationTestHarness::fitsHeaderT firstHeader;
    firstHeader.append( "DATE-OBS", "2024-01-02T03:04:05", "observation time" );
    firstHeader.append<int>( "KEEP", 17, "selected value" );
    HCIobservationTestHarness::fitsHeaderT secondHeader;
    secondHeader.append( "DATE-OBS", "2024-01-02T03:04:06", "observation time" );
    secondHeader.append<int>( "KEEP", 29, "selected value" );

    const auto firstPath = directory.file( "first.fits" );
    const auto secondPath = directory.file( "second.fits" );
    writeFitsImage( firstPath, first, &firstHeader );
    writeFitsImage( secondPath, second, &secondHeader );

    observation.m_fileList = { secondPath.string(), firstPath.string() };
    observation.m_imSize = 3;
    observation.m_dateKeyword = "DATE-OBS";
    observation.m_dateIsISO8601 = true;
    observation.m_keywords = { "KEEP" };
    observation.m_skipPreProcess = true;
    observation.readFiles();

    REQUIRE( observation.m_filesRead );
    REQUIRE( observation.m_filesDeleted );
    REQUIRE( observation.m_Nims == 2 );
    REQUIRE( observation.m_Nrows == 3 );
    REQUIRE( observation.m_Ncols == 3 );
    REQUIRE( observation.m_Npix == 9 );
    REQUIRE( observation.m_tgtIms.image( 0 )( 0, 0 ) == Approx( 111 ) );
    REQUIRE( observation.m_tgtIms.image( 0 )( 2, 2 ) == Approx( 133 ) );
    REQUIRE( observation.m_tgtIms.image( 1 )( 0, 0 ) == Approx( 11 ) );
    REQUIRE( observation.m_tgtIms.image( 1 )( 1, 1 ) == 0 );
    REQUIRE( observation.m_heads[0]["KEEP"].Int() == 29 );
    REQUIRE( observation.m_heads[1]["KEEP"].Int() == 17 );
    REQUIRE( observation.m_imageMJD.size() == 2 );
    REQUIRE( observation.m_imageMJD[0] == Approx( mx::sys::ISO8601date2mjd( "2024-01-02T03:04:06" ) ) );
    REQUIRE( observation.m_hookEvents == std::vector<std::string>{ "target-read" } );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::readFiles();
#endif
    // clang-format on
}

/// Verify HCIobservation::readFiles clamps oversized requests, scales numeric dates, and clears stale date state.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation repeated target FITS ingestion", "[HCIobservation][readFiles][state]" )
{
    TestDirectory directory;
    HCIobservationTestHarness observation;

    const auto path = directory.file( "numeric-date.fits" );
    auto image = indexedImage( 2, 4 );
    HCIobservationTestHarness::fitsHeaderT header;
    header.append<float>( "STAMP", 3, "numeric time" );
    writeFitsImage( path, image, &header );

    observation.m_fileList = { path.string() };
    observation.m_imSize = 9;
    observation.m_dateKeyword = "STAMP";
    observation.m_dateIsISO8601 = false;
    observation.m_dateUnit = 2;
    observation.m_skipPreProcess = true;
    observation.readFiles();

    REQUIRE( observation.m_imSize == 2 );
    REQUIRE( observation.m_tgtIms.image( 0 )( 0, 0 ) == Approx( 1 ) );
    REQUIRE( observation.m_tgtIms.image( 0 )( 1, 1 ) == Approx( 12 ) );
    REQUIRE( observation.m_imageMJD == std::vector<double>{ 6 } );

    observation.m_dateKeyword.clear();
    observation.m_filesDeleted = false;
    observation.readFiles();
    REQUIRE( observation.m_imageMJD.empty() );
    REQUIRE( observation.m_hookEvents == std::vector<std::string>{ "target-read", "target-read" } );
}

/// Verify HCIobservation::readFiles rejects invalid deletion ranges without partially mutating the list.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation target deletion bounds", "[HCIobservation][readFiles][validation]" )
{
    HCIobservationTestHarness observation;
    const std::vector<std::string> original{ "first.fits", "second.fits" };
    observation.m_fileList = original;
    observation.m_deleteFront = 1;
    observation.m_deleteBack = 2;

    REQUIRE_THROWS( observation.readFiles() );
    REQUIRE( observation.m_fileList == original );
    REQUIRE_FALSE( observation.m_filesDeleted );
    REQUIRE_FALSE( observation.m_filesRead );

    observation.m_deleteFront = -1;
    observation.m_deleteBack = 0;
    REQUIRE_THROWS( observation.readFiles() );
    REQUIRE( observation.m_fileList == original );
}

/// Verify HCIobservation::readFiles propagates failures from the first and later FITS inputs.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation target FITS failures", "[HCIobservation][readFiles][FITS][validation]" )
{
    TestDirectory directory;
    HCIobservationTestHarness observation;
    observation.m_dateKeyword.clear();
    observation.m_skipPreProcess = true;
    observation.m_fileList = { directory.file( "missing-first.fits" ).string() };

    REQUIRE_THROWS( observation.readFiles() );
    REQUIRE_FALSE( observation.m_filesRead );

    const auto validPath = directory.file( "valid.fits" );
    auto valid = indexedImage( 2, 2 );
    writeFitsImage( validPath, valid );
    observation.m_fileList = { validPath.string(), directory.file( "missing-later.fits" ).string() };
    observation.m_filesDeleted = false;

    REQUIRE_THROWS( observation.readFiles() );
    REQUIRE_FALSE( observation.m_filesRead );
}

/// Verify HCIobservation::readFiles honors preprocessing skip and invokes post-coadd only when coaddition runs.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation target read hook order", "[HCIobservation][readFiles][hooks]" )
{
    TestDirectory directory;
    const auto firstPath = directory.file( "first.fits" );
    const auto secondPath = directory.file( "second.fits" );
    auto first = indexedImage( 2, 2 );
    auto second = indexedImage( 2, 2, 10 );
    writeFitsImage( firstPath, first );
    writeFitsImage( secondPath, second );

    HCIobservationTestHarness skipped;
    skipped.m_fileList = { firstPath.string(), secondPath.string() };
    skipped.m_dateKeyword.clear();
    skipped.m_skipPreProcess = true;
    skipped.m_coaddMethod = mx::improc::HCI::coadd::mean;
    skipped.m_coaddMaxImno = 2;
    skipped.readFiles();
    REQUIRE( skipped.m_tgtIms.planes() == 2 );
    REQUIRE( skipped.m_hookEvents == std::vector<std::string>{ "target-read" } );

    HCIobservationTestHarness processed;
    processed.m_fileList = { firstPath.string(), secondPath.string() };
    processed.m_dateKeyword.clear();
    processed.m_coaddMethod = mx::improc::HCI::coadd::mean;
    processed.m_coaddMaxImno = 2;
    processed.readFiles();
    REQUIRE( processed.m_tgtIms.planes() == 1 );
    REQUIRE( processed.m_hookEvents == std::vector<std::string>{ "target-read", "target-coadd" } );
}

/// Verify HCIobservation::readMask centers rectangular masks and HCIobservation::makeMaskCube replicates them.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation mask ingestion", "[HCIobservation][readMask][makeMaskCube]" )
{
    TestDirectory directory;
    HCIobservationTestHarness observation;
    observation.m_imSize = 3;
    observation.m_Nrows = 3;
    observation.m_Ncols = 3;
    observation.m_Nims = 2;
    observation.m_maskFile = directory.file( "mask.fits" ).string();

    auto mask = indexedImage( 5, 7 );
    writeFitsImage( observation.m_maskFile, mask );
    observation.readMask();

    REQUIRE( observation.m_mask.rows() == 3 );
    REQUIRE( observation.m_mask.cols() == 3 );
    REQUIRE( observation.m_mask( 0, 0 ) == Approx( 12 ) );
    REQUIRE( observation.m_mask( 2, 2 ) == Approx( 34 ) );
    REQUIRE( observation.m_maskCube.planes() == 2 );
    REQUIRE( observation.m_maskCube.image( 0 ).isApprox( observation.m_mask ) );
    REQUIRE( observation.m_maskCube.image( 1 ).isApprox( observation.m_mask ) );

    observation.m_Nims = 1;
    observation.makeMaskCube();
    REQUIRE( observation.m_maskCube.planes() == 1 );

    auto tooShort = indexedImage( 2, 5 );
    writeFitsImage( directory.file( "short-mask.fits" ), tooShort );
    observation.m_maskFile = directory.file( "short-mask.fits" ).string();
    REQUIRE_THROWS( observation.readMask() );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::readMask();
    mx::improc::HCIobservation<float, mx::verbose::vv>::makeMaskCube();
#endif
    // clang-format on
}

} // namespace HCIobservation_read_test
} // namespace unitTest
