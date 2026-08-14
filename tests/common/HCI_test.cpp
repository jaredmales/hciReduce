/** \file HCI_test.cpp
 * \brief Tests high-contrast imaging configuration enum conversions.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "src/common/HCI.hpp"

namespace unitTest
{
namespace HCI_test
{

/// Verify every coaddition method converts bidirectionally and invalid inputs are rejected.
/** \ingroup HCI_unit_tests */
TEST_CASE( "HCI coadd method conversions", "[HCI][coadd]" )
{
    using namespace mx::improc::HCI;
    const std::vector<std::pair<coadd, std::string>> methods{ { coadd::none, "none" },
                                                              { coadd::median, "median" },
                                                              { coadd::mean, "mean" } };

    for( const auto &[method, name] : methods )
    {
        REQUIRE( coaddToStr<mx::verbose::vv>( method ) == name );
        REQUIRE( coaddFmStr<mx::verbose::vv>( name ) == method );
    }

    REQUIRE_THROWS( coaddToStr<mx::verbose::vv>( coadd::invalid ) );
    REQUIRE_THROWS( coaddFmStr<mx::verbose::vv>( "invalid" ) );
}

/// Verify every mean-subtraction method converts bidirectionally and invalid inputs are rejected.
/** \ingroup HCI_unit_tests */
TEST_CASE( "HCI mean-subtraction method conversions", "[HCI][meanSub]" )
{
    using namespace mx::improc::HCI;
    const std::vector<std::pair<meanSub, std::string>> methods{ { meanSub::none, "none" },
                                                                { meanSub::meanImage, "meanImage" },
                                                                { meanSub::medianImage, "medianImage" },
                                                                { meanSub::imageMean, "imageMean" },
                                                                { meanSub::imageMedian, "imageMedian" },
                                                                { meanSub::imageMode, "imageMode" } };

    for( const auto &[method, name] : methods )
    {
        REQUIRE( meanSubToStr<mx::verbose::vv>( method ) == name );
        REQUIRE( meanSubFmStr<mx::verbose::vv>( name ) == method );
    }

    REQUIRE_THROWS( meanSubToStr<mx::verbose::vv>( static_cast<meanSub>( 99 ) ) );
    REQUIRE_THROWS( meanSubFmStr<mx::verbose::vv>( "invalid" ) );
}

/// Verify every pixel-time-series normalization converts bidirectionally and invalid inputs are rejected.
/** \ingroup HCI_unit_tests */
TEST_CASE( "HCI pixel normalization conversions", "[HCI][pixelTSNorm]" )
{
    using namespace mx::improc::HCI;
    const std::vector<std::pair<pixelTSNorm, std::string>> methods{
        { pixelTSNorm::none, "none" },
        { pixelTSNorm::rms, "rms" },
        { pixelTSNorm::rmsSigmaClipped, "rmsSigmaClipped" } };

    for( const auto &[method, name] : methods )
    {
        REQUIRE( pixelTSNormToStr<mx::verbose::vv>( method ) == name );
        REQUIRE( pixelTSNormFmStr<mx::verbose::vv>( name ) == method );
    }

    REQUIRE_THROWS( pixelTSNormToStr<mx::verbose::vv>( static_cast<pixelTSNorm>( 99 ) ) );
    REQUIRE_THROWS( pixelTSNormFmStr<mx::verbose::vv>( "invalid" ) );
}

/// Verify every final-combination method converts bidirectionally and invalid inputs are rejected.
/** \ingroup HCI_unit_tests */
TEST_CASE( "HCI combination method conversions", "[HCI][combine]" )
{
    using namespace mx::improc::HCI;
    const std::vector<std::pair<combine, std::string>> methods{ { combine::none, "none" },
                                                                { combine::median, "median" },
                                                                { combine::mean, "mean" },
                                                                { combine::sigmaMean, "sigmaMean" } };

    for( const auto &[method, name] : methods )
    {
        REQUIRE( combineToStr<mx::verbose::vv>( method ) == name );
        REQUIRE( combineFmStr<mx::verbose::vv>( name ) == method );
    }

    REQUIRE_THROWS( combineToStr<mx::verbose::vv>( static_cast<combine>( 99 ) ) );
    REQUIRE_THROWS( combineFmStr<mx::verbose::vv>( "invalid" ) );
}

/// Verify every fake-PSF method converts bidirectionally and invalid inputs are rejected.
/** \ingroup HCI_unit_tests */
TEST_CASE( "HCI fake method conversions", "[HCI][fake]" )
{
    using namespace mx::improc::HCI;
    const std::vector<std::pair<fake, std::string>> methods{ { fake::single, "single" }, { fake::list, "list" } };

    for( const auto &[method, name] : methods )
    {
        REQUIRE( fakeToStr<mx::verbose::vv>( method ) == name );
        REQUIRE( fakeFmStr<mx::verbose::vv>( name ) == method );
    }

    REQUIRE_THROWS( fakeToStr<mx::verbose::vv>( static_cast<fake>( 99 ) ) );
    REQUIRE_THROWS( fakeFmStr<mx::verbose::vv>( "invalid" ) );
}

/// Verify every reference-exclusion method converts bidirectionally and invalid inputs are rejected.
/** \ingroup HCI_unit_tests */
TEST_CASE( "HCI exclusion method conversions", "[HCI][exclude]" )
{
    using namespace mx::improc::HCI;
    const std::vector<std::pair<exclude, std::string>> methods{ { exclude::none, "none" },
                                                                { exclude::pixel, "pixel" },
                                                                { exclude::angle, "angle" },
                                                                { exclude::imno, "imno" } };

    for( const auto &[method, name] : methods )
    {
        REQUIRE( excludeToStr<mx::verbose::vv>( method ) == name );
        REQUIRE( excludeFmStr<mx::verbose::vv>( name ) == method );
    }

    REQUIRE_THROWS( excludeToStr<mx::verbose::vv>( static_cast<exclude>( 99 ) ) );
    REQUIRE_THROWS( excludeFmStr<mx::verbose::vv>( "invalid" ) );
}

/// Verify every reference-inclusion method converts bidirectionally and invalid inputs are rejected.
/** \ingroup HCI_unit_tests */
TEST_CASE( "HCI inclusion method conversions", "[HCI][include]" )
{
    using namespace mx::improc::HCI;
    const std::vector<std::pair<include, std::string>> methods{ { include::all, "all" },
                                                                { include::corr, "corr" },
                                                                { include::time, "time" },
                                                                { include::angle, "angle" },
                                                                { include::imno, "imno" } };

    for( const auto &[method, name] : methods )
    {
        REQUIRE( includeToStr<mx::verbose::vv>( method ) == name );
        REQUIRE( includeFmStr<mx::verbose::vv>( name ) == method );
    }

    try
    {
        includeToStr<mx::verbose::vv>( static_cast<include>( 99 ) );
        FAIL( "invalid inclusion enum was accepted" );
    }
    catch( const mx::exception<mx::verbose::vv> &error )
    {
        REQUIRE( error.message() == "got an invalid include method (bug)" );
    }
    REQUIRE_THROWS( includeFmStr<mx::verbose::vv>( "invalid" ) );
}

} // namespace HCI_test
} // namespace unitTest
