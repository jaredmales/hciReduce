/** \file HCIobservation_read_test.cpp
 * \brief Tests HCIobservation target FITS ingestion and mask construction.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "HCIobservation_test_fixture.hpp"

#include <cmath>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <sys/wait.h>
#include <unistd.h>

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

struct throwingTargetHookHarness : public HCIobservationTestHarness
{
    bool m_throwPostRead{ false };
    bool m_throwPostCoadd{ false };

    void postReadFiles() override
    {
        if( m_throwPostRead )
        {
            throw std::runtime_error( "post-read test failure" );
        }
        HCIobservationTestHarness::postReadFiles();
    }

    void postCoadd() override
    {
        if( m_throwPostCoadd )
        {
            throw std::runtime_error( "post-coadd test failure" );
        }
        HCIobservationTestHarness::postCoadd();
    }
};

struct secondMaskFailureHarness : public HCIobservationTestHarness
{
    int m_maskCalls{ 0 };

    void makeMaskCube() override
    {
        ++m_maskCalls;
        if( m_maskCalls == 2 )
        {
            throw std::runtime_error( "second mask-cube test failure" );
        }
        HCIobservationTestHarness::makeMaskCube();
    }
};

struct firstMaskFailureHarness : public HCIobservationTestHarness
{
    void makeMaskCube() override
    {
        throw std::runtime_error( "first mask-cube test failure" );
    }
};
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

    const auto narrowPath = directory.file( "narrow.fits" );
    writeFitsImage( narrowPath, indexedImage( 5, 3 ) );
    HCIobservationTestHarness narrow;
    narrow.m_fileList = { narrowPath.string() };
    narrow.m_imSize = 4;
    narrow.m_dateKeyword.clear();
    narrow.m_skipPreProcess = true;
    narrow.readFiles();
    REQUIRE( narrow.m_imSize == 3 );
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

    HCIobservationTestHarness empty;
    REQUIRE_THROWS( empty.readFiles() );

    HCIobservationTestHarness emptiedByDeletion;
    emptiedByDeletion.m_fileList = { "only.fits" };
    emptiedByDeletion.m_deleteFront = 1;
    REQUIRE_THROWS( emptiedByDeletion.readFiles() );
    REQUIRE( emptiedByDeletion.m_fileList.empty() );
}

/// Verify HCIobservation::readFiles applies valid list deletions and quality selection before automatic cropping.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation target selection and automatic size", "[HCIobservation][readFiles][quality][crop]" )
{
    TestDirectory directory;
    const auto firstPath = directory.file( "first.fits" );
    const auto middlePath = directory.file( "middle.fits" );
    const auto lastPath = directory.file( "last.fits" );
    writeFitsImage( firstPath, indexedImage( 5, 3 ) );
    writeFitsImage( middlePath, indexedImage( 5, 3, 100 ) );
    writeFitsImage( lastPath, indexedImage( 5, 3, 200 ) );

    const auto qualityPath = directory.file( "quality.txt" );
    writeTextFile( qualityPath, "middle.fits 2\n" );

    HCIobservationTestHarness observation;
    observation.m_fileList = { firstPath.string(), middlePath.string(), lastPath.string() };
    observation.m_deleteFront = 1;
    observation.m_deleteBack = 1;
    observation.m_qualityFile = qualityPath.string();
    observation.m_qualityThreshold = 2;
    observation.m_dateKeyword.clear();
    observation.m_skipPreProcess = true;

    observation.readFiles();

    REQUIRE( observation.m_fileList == std::vector<std::string>{ middlePath.string() } );
    REQUIRE( observation.m_imSize == 3 );
    REQUIRE( observation.m_tgtIms.rows() == 3 );
    REQUIRE( observation.m_tgtIms.cols() == 3 );
    REQUIRE( observation.m_tgtIms.image( 0 )( 0, 0 ) == Approx( 110 ) );
    REQUIRE( observation.m_hookEvents == std::vector<std::string>{ "target-read" } );
}

/// Verify HCIobservation::readFiles threshold-only mode prints the selected list and exits successfully.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation threshold-only mode", "[HCIobservation][readFiles][quality][thresholdOnly]" )
{
    TestDirectory directory;
    const auto imagePath = directory.file( "selected.fits" );
    const auto qualityPath = directory.file( "quality.txt" );
    writeFitsImage( imagePath, indexedImage( 1, 1 ) );
    writeTextFile( qualityPath, "selected.fits 2\n" );

    const pid_t child = fork();
    REQUIRE( child >= 0 );
    if( child == 0 )
    {
        HCIobservationTestHarness observation;
        observation.m_fileList = { imagePath.string() };
        observation.m_qualityFile = qualityPath.string();
        observation.m_qualityThreshold = 1;
        observation.m_thresholdOnly = true;
        try
        {
            observation.readFiles();
        }
        catch( ... )
        {
            std::_Exit( 2 );
        }
        std::_Exit( 3 );
    }

    int status = 0;
    REQUIRE( waitpid( child, &status, 0 ) == child );
    REQUIRE( WIFEXITED( status ) );
    REQUIRE( WEXITSTATUS( status ) == 0 );
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

    HCIobservationTestHarness badWeights;
    badWeights.m_fileList = { validPath.string() };
    badWeights.m_weightFile = directory.file( "missing-weights.txt" ).string();
    badWeights.m_skipPreProcess = true;
    REQUIRE_THROWS( badWeights.readFiles() );
    REQUIRE_FALSE( badWeights.m_filesRead );
}

/// Verify HCIobservation::readFiles adds operation context to failures from each delegated processing stage.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation target delegated failures", "[HCIobservation][readFiles][validation][delegation]" )
{
    TestDirectory directory;
    const auto imagePath = directory.file( "valid.fits" );
    writeFitsImage( imagePath, indexedImage( 3, 3 ) );

    HCIobservationTestHarness thresholdFailure;
    thresholdFailure.m_fileList = { imagePath.string() };
    thresholdFailure.m_qualityFile = directory.file( "missing-quality.txt" ).string();
    thresholdFailure.m_qualityThreshold = 1;
    REQUIRE_THROWS( thresholdFailure.readFiles() );

    const auto qualityPath = directory.file( "quality.txt" );
    writeTextFile( qualityPath, "valid.fits 1\n" );
    HCIobservationTestHarness thresholdEmpty;
    thresholdEmpty.m_fileList = { imagePath.string() };
    thresholdEmpty.m_qualityFile = qualityPath.string();
    thresholdEmpty.m_qualityThreshold = 2;
    REQUIRE_THROWS( thresholdEmpty.readFiles() );

    throwingTargetHookHarness postReadFailure;
    postReadFailure.m_fileList = { imagePath.string() };
    postReadFailure.m_dateKeyword.clear();
    postReadFailure.m_skipPreProcess = true;
    postReadFailure.m_throwPostRead = true;
    REQUIRE_THROWS( postReadFailure.readFiles() );

    HCIobservationTestHarness maskFailure;
    maskFailure.m_fileList = { imagePath.string() };
    maskFailure.m_dateKeyword.clear();
    maskFailure.m_maskFile = directory.file( "missing-mask.fits" ).string();
    REQUIRE_THROWS( maskFailure.readFiles() );

    HCIobservationTestHarness beforeCoaddFailure;
    beforeCoaddFailure.m_fileList = { imagePath.string() };
    beforeCoaddFailure.m_dateKeyword.clear();
    beforeCoaddFailure.m_preProcess_beforeCoadd = true;
    beforeCoaddFailure.m_preProcess_meanSubMethod = static_cast<mx::improc::HCI::meanSub>( 99 );
    REQUIRE_THROWS( beforeCoaddFailure.readFiles() );

    HCIobservationTestHarness coaddFailure;
    coaddFailure.m_fileList = { imagePath.string() };
    coaddFailure.m_dateKeyword.clear();
    coaddFailure.m_coaddMethod = static_cast<mx::improc::HCI::coadd>( 99 );
    coaddFailure.m_coaddMaxImno = 1;
    REQUIRE_THROWS( coaddFailure.readFiles() );

    throwingTargetHookHarness postCoaddFailure;
    postCoaddFailure.m_fileList = { imagePath.string() };
    postCoaddFailure.m_dateKeyword.clear();
    postCoaddFailure.m_coaddMethod = mx::improc::HCI::coadd::mean;
    postCoaddFailure.m_coaddMaxImno = 1;
    postCoaddFailure.m_throwPostCoadd = true;
    REQUIRE_THROWS( postCoaddFailure.readFiles() );

    HCIobservationTestHarness afterCoaddFailure;
    afterCoaddFailure.m_fileList = { imagePath.string() };
    afterCoaddFailure.m_dateKeyword.clear();
    afterCoaddFailure.m_preProcess_meanSubMethod = static_cast<mx::improc::HCI::meanSub>( 99 );
    REQUIRE_THROWS( afterCoaddFailure.readFiles() );

    const auto outputBlocker = directory.file( "output-blocker" );
    writeTextFile( outputBlocker, "regular file" );
    HCIobservationTestHarness outputFailure;
    outputFailure.m_fileList = { imagePath.string() };
    outputFailure.m_dateKeyword.clear();
    outputFailure.m_preProcess_outputPrefix = ( outputBlocker / "target_" ).string();
    REQUIRE_THROWS( outputFailure.readFiles() );

    const auto maskPath = directory.file( "mask.fits" );
    HCIobservationTestHarness::imageT mask( 3, 3 );
    mask.setOnes();
    writeFitsImage( maskPath, mask );
    secondMaskFailureHarness maskRebuildFailure;
    maskRebuildFailure.m_fileList = { imagePath.string() };
    maskRebuildFailure.m_dateKeyword.clear();
    maskRebuildFailure.m_maskFile = maskPath.string();
    maskRebuildFailure.m_coaddMethod = mx::improc::HCI::coadd::mean;
    maskRebuildFailure.m_coaddMaxImno = 1;
    REQUIRE_THROWS( maskRebuildFailure.readFiles() );
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

    observation.m_maskFile = directory.file( "missing-mask.fits" ).string();
    REQUIRE_THROWS( observation.readMask() );

    firstMaskFailureHarness maskCubeFailure;
    maskCubeFailure.m_imSize = 3;
    maskCubeFailure.m_maskFile = directory.file( "mask.fits" ).string();
    REQUIRE_THROWS( maskCubeFailure.readMask() );

    observation.m_mask.resize( 2, 3 );
    observation.m_Nrows = 3;
    observation.m_Ncols = 3;
    observation.m_Nims = 1;
    REQUIRE_THROWS( observation.makeMaskCube() );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::readMask();
    mx::improc::HCIobservation<float, mx::verbose::vv>::makeMaskCube();
#endif
    // clang-format on
}

} // namespace HCIobservation_read_test
} // namespace unitTest
