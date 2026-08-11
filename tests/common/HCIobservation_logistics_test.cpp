/** \file HCIobservation_logistics_test.cpp
 * \brief Tests HCIobservation file-list, threshold, and weight logistics.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "HCIobservation_test_fixture.hpp"

#include <limits>

namespace unitTest
{
namespace HCIobservation_logistics_test
{

/// Verify that HCIobservation::load_fileList scans, sorts, prefixes, and resets lists correctly.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation file-list loading", "[HCIobservation][files]" )
{
    TestDirectory directory;
    HCIobservationTestHarness observation;

    writeTextFile( directory.file( "target_b.fits" ), "" );
    writeTextFile( directory.file( "target_a.fits" ), "" );
    writeTextFile( directory.file( "target_ignore.txt" ), "" );
    writeTextFile( directory.file( "other.fits" ), "" );

    std::vector<std::string> files{ "stale" };
    REQUIRE( observation.load_fileList( files, "", directory.path().string(), "target_", "fits" ) ==
             mx::error_t::noerror );
    REQUIRE( files == std::vector<std::string>{ directory.file( "target_a.fits" ).string(),
                                                directory.file( "target_b.fits" ).string() } );

    const auto listPath = directory.file( "targets.list" );
    writeTextFile( listPath, "target_b.fits\ntarget_a.fits\n" );
    REQUIRE( observation.load_fileList( files, listPath.string(), directory.path().string(), "", "" ) ==
             mx::error_t::noerror );
    REQUIRE( files == std::vector<std::string>{ directory.file( "target_b.fits" ).string(),
                                                directory.file( "target_a.fits" ).string() } );

    writeTextFile( listPath, "target_a.fits\n" );
    REQUIRE( observation.load_fileList( files, listPath.string(), "", "", "" ) == mx::error_t::noerror );
    REQUIRE( files == std::vector<std::string>{ "target_a.fits" } );

    REQUIRE( observation.load_fileList( files, "", "", "", "" ) == mx::error_t::noerror );
    REQUIRE( files.empty() );

    REQUIRE( observation.load_fileList( files, "", directory.file( "missing" ).string(), "", "fits" ) !=
             mx::error_t::noerror );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::load_fileList( files, "", "", "", "" );
#endif
    // clang-format on
}

/// Verify that HCIobservation::load_RDIfileList keeps target and RDI list state independent.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation RDI file-list independence", "[HCIobservation][files][RDI]" )
{
    TestDirectory directory;
    HCIobservationTestHarness observation;

    const auto targetList = directory.file( "targets.list" );
    const auto rdiList = directory.file( "rdi.list" );
    writeTextFile( targetList, "target.fits\n" );
    writeTextFile( rdiList, "reference.fits\n" );

    observation.m_fileListFile = targetList.string();
    observation.m_directory = directory.file( "target-dir" ).string();
    observation.m_filesDeleted = true;
    observation.m_filesRead = true;
    REQUIRE( observation.load_fileList() == mx::error_t::noerror );
    const auto targets = observation.m_fileList;
    REQUIRE_FALSE( observation.m_filesDeleted );
    REQUIRE_FALSE( observation.m_filesRead );

    observation.m_RDIfileListFile = rdiList.string();
    observation.m_RDIdirectory = directory.file( "rdi-dir" ).string();
    observation.m_RDIfilesDeleted = true;
    observation.m_RDIfilesRead = true;
    REQUIRE( observation.load_RDIfileList() == mx::error_t::noerror );
    REQUIRE_FALSE( observation.m_RDIfilesDeleted );
    REQUIRE_FALSE( observation.m_RDIfilesRead );
    REQUIRE( observation.m_fileList == targets );
    REQUIRE( observation.m_RDIfileList ==
             std::vector<std::string>{ directory.file( "rdi-dir/reference.fits" ).string() } );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::load_fileList();
    mx::improc::HCIobservation<float, mx::verbose::vv>::load_RDIfileList();
#endif
    // clang-format on
}

/// Verify threshold boundary, ordering, basename matching, validation, and disabled behavior in
/// HCIobservation::threshold.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation quality thresholding", "[HCIobservation][threshold]" )
{
    TestDirectory directory;
    HCIobservationTestHarness observation;
    std::vector<std::string> files{ "/data/a.fits", "/data/b.fits", "/data/c.fits" };

    const auto quality = directory.file( "quality.txt" );
    writeTextFile( quality, "elsewhere/c.fits 3\na.fits 1\nb.fits 2\n" );
    observation.threshold( files, quality.string(), 2 );
    REQUIRE( files == std::vector<std::string>{ "/data/b.fits", "/data/c.fits" } );

    files = { "/data/a.fits" };
    observation.threshold( files, "", 0 );
    REQUIRE( files == std::vector<std::string>{ "/data/a.fits" } );

    REQUIRE_THROWS( observation.threshold( files, "", 1 ) );

    writeTextFile( quality, "other.fits 2\n" );
    REQUIRE_THROWS( observation.threshold( files, quality.string(), 1 ) );

    writeTextFile( quality, "a.fits 1\na.fits 2\n" );
    REQUIRE_THROWS( observation.threshold( files, quality.string(), 1 ) );

    writeTextFile( quality, "a.fits nan\n" );
    REQUIRE_THROWS( observation.threshold( files, quality.string(), 1 ) );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::threshold( files, quality.string(), 1 );
#endif
    // clang-format on
}

/// Verify HCIobservation::readWeights maps by basename, normalizes, and rejects ambiguous or invalid tables.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation weight loading", "[HCIobservation][weights]" )
{
    TestDirectory directory;
    HCIobservationTestHarness observation;
    observation.m_fileList = { "/data/a.fits", "/data/b.fits" };
    observation.m_weightFile = directory.file( "weights.txt" ).string();

    writeTextFile( observation.m_weightFile, "b.fits 3\na.fits 1\n" );
    observation.readWeights();
    REQUIRE( observation.m_comboWeights.size() == 2 );
    REQUIRE( observation.m_comboWeights[0] == Approx( 0.25 ) );
    REQUIRE( observation.m_comboWeights[1] == Approx( 0.75 ) );

    writeTextFile( observation.m_weightFile, "a.fits 1\n" );
    REQUIRE_THROWS( observation.readWeights() );

    writeTextFile( observation.m_weightFile, "a.fits 1\na.fits 2\nb.fits 3\n" );
    REQUIRE_THROWS( observation.readWeights() );

    writeTextFile( observation.m_weightFile, "a.fits 1\nb.fits -1\n" );
    REQUIRE_THROWS( observation.readWeights() );

    writeTextFile( observation.m_weightFile, "a.fits nan\nb.fits 1\n" );
    REQUIRE_THROWS( observation.readWeights() );

    observation.m_weightFile.clear();
    REQUIRE_THROWS( observation.readWeights() );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::readWeights();
#endif
    // clang-format on
}

} // namespace HCIobservation_logistics_test
} // namespace unitTest
