/** \file P4PSFFilter_test.cpp
 * \brief Tests normalized local filtering with spatially variable P4 PSFs.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "src/common/P4PSFFilter.hpp"

#include <limits>

namespace unitTest
{
namespace P4PSFFilter_test
{

using filterT = mx::improc::P4PSFFilter;
using imageT = filterT::imageT;
using validityT = filterT::validityT;

/// Verify the production filter returns the normalized amplitude of a signed asymmetric response.
/** This exercises mx::improc::P4PSFFilter::calculate() with complete support.
 * \ingroup P4PSFFilter_unit_tests
 */
TEST_CASE( "P4 local PSF filter estimates normalized amplitude", "[P4PSFFilter][amplitude][signed]" )
{
    imageT response( 3, 3 );
    response << 0.1F, -0.2F, 0.05F, -0.3F, 1.0F, 0.2F, 0.08F, -0.1F, 0.04F;
    const validityT validity = validityT::Ones( 3, 3 );
    imageT science = imageT::Zero( 7, 7 );
    science.block( 2, 3, 3, 3 ) = 2.5F * response;

    const mx::improc::P4PSFFilterResult result = filterT::calculate( science, response, validity, 3, 4, 1.0 );
    REQUIRE( result.valid );
    REQUIRE( result.amplitude == Approx( 2.5 ).margin( 1e-6 ) );
    REQUIRE( result.correlation == Approx( 2.5 * result.normalization ).margin( 1e-6 ) );
    REQUIRE( result.normalization > 0 );
    REQUIRE( result.supportFraction == Approx( 1.0 ) );
}

/// Verify edge, PSF-validity, and finite-science samples control local filter support.
/** This exercises mx::improc::P4PSFFilter::calculate() with partial support and its configured threshold.
 * \ingroup P4PSFFilter_unit_tests
 */
TEST_CASE( "P4 local PSF filter reports partial support", "[P4PSFFilter][support][edge][validity]" )
{
    const imageT response = imageT::Ones( 3, 3 );
    validityT validity = validityT::Ones( 3, 3 );
    imageT science = imageT::Constant( 4, 4, 3.0F );

    mx::improc::P4PSFFilterResult result = filterT::calculate( science, response, validity, 0, 0, 0.4 );
    REQUIRE( result.valid );
    REQUIRE( result.amplitude == Approx( 3.0 ) );
    REQUIRE( result.normalization == Approx( 4.0 ) );
    REQUIRE( result.supportFraction == Approx( 4.0 / 9.0 ) );

    result = filterT::calculate( science, response, validity, 0, 0, 0.5 );
    REQUIRE_FALSE( result.valid );

    validity( 1, 1 ) = 0;
    science( 0, 1 ) = std::numeric_limits<float>::quiet_NaN();
    result = filterT::calculate( science, response, validity, 0, 0, 0.2 );
    REQUIRE( result.valid );
    REQUIRE( result.supportFraction == Approx( 2.0 / 9.0 ) );
}

/// Verify degenerate and malformed local filtering inputs are rejected or explicitly invalid.
/** This exercises mx::improc::P4PSFFilter::calculate() validation and zero-normalization behavior.
 * \ingroup P4PSFFilter_unit_tests
 */
TEST_CASE( "P4 local PSF filter validates its contract", "[P4PSFFilter][validation]" )
{
    imageT science = imageT::Ones( 5, 5 );
    imageT response = imageT::Zero( 3, 3 );
    validityT validity = validityT::Ones( 3, 3 );
    REQUIRE_FALSE( filterT::calculate( science, response, validity, 2, 2, 1.0 ).valid );

    imageT empty;
    REQUIRE_THROWS_AS( filterT::calculate( empty, response, validity, 0, 0, 1.0 ), std::invalid_argument );
    REQUIRE_THROWS_AS( filterT::calculate( science, imageT::Zero( 2, 3 ), validityT::Ones( 2, 3 ), 2, 2, 1.0 ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( filterT::calculate( science, response, validityT::Ones( 2, 3 ), 2, 2, 1.0 ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( filterT::calculate( science, response, validity, -1, 2, 1.0 ), std::out_of_range );
    REQUIRE_THROWS_AS( filterT::calculate( science, response, validity, 2, 2, -0.1 ), std::invalid_argument );
    REQUIRE_THROWS_AS(
        filterT::calculate( science, response, validity, 2, 2, std::numeric_limits<double>::quiet_NaN() ),
        std::invalid_argument );

    response( 1, 1 ) = std::numeric_limits<float>::infinity();
    REQUIRE_THROWS_AS( filterT::calculate( science, response, validity, 2, 2, 1.0 ), std::invalid_argument );
}

} // namespace P4PSFFilter_test
} // namespace unitTest
