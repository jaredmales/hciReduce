/** \file HCIobservation_preprocess_test.cpp
 * \brief Tests HCIobservation pure preprocessing operations.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "HCIobservation_test_fixture.hpp"

#include <cmath>

namespace unitTest
{
namespace HCIobservation_preprocess_test
{

/// \cond HCIobservation_test_harness
namespace
{
void makeMeanSubCube( HCIobservationTestHarness::cubeT &cube )
{
    cube.resize( 1, 2, 3 );
    cube.image( 0 ) << 1, 10;
    cube.image( 1 ) << 2, 20;
    cube.image( 2 ) << 6, 30;
}
} // namespace
/// \endcond

/// Verify HCIobservation::preProcess_meanSub mean, median, no-op, mask, and validation behavior.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation mean-image preprocessing", "[HCIobservation][preprocess][mean]" )
{
    OpenMPThreadGuard threads;
    HCIobservationTestHarness observation;

    HCIobservationTestHarness::cubeT cube;
    makeMeanSubCube( cube );
    HCIobservationTestHarness::cubeT original;
    original = cube;
    observation.m_preProcess_meanSubMethod = mx::improc::HCI::meanSub::none;
    observation.preProcess_meanSub( cube );
    REQUIRE( cube.image( 0 ).isApprox( original.image( 0 ) ) );
    REQUIRE( cube.image( 1 ).isApprox( original.image( 1 ) ) );
    REQUIRE( cube.image( 2 ).isApprox( original.image( 2 ) ) );

    makeMeanSubCube( cube );
    observation.m_preProcess_meanSubMethod = mx::improc::HCI::meanSub::meanImage;
    observation.preProcess_meanSub( cube );
    REQUIRE( cube.image( 0 )( 0, 0 ) == Approx( -2 ) );
    REQUIRE( cube.image( 1 )( 0, 0 ) == Approx( -1 ) );
    REQUIRE( cube.image( 2 )( 0, 0 ) == Approx( 3 ) );
    REQUIRE( cube.image( 0 )( 0, 1 ) == Approx( -10 ) );
    REQUIRE( cube.image( 1 )( 0, 1 ) == Approx( 0 ) );
    REQUIRE( cube.image( 2 )( 0, 1 ) == Approx( 10 ) );

    makeMeanSubCube( cube );
    observation.m_preProcess_meanSubMethod = mx::improc::HCI::meanSub::medianImage;
    observation.preProcess_meanSub( cube );
    REQUIRE( cube.image( 0 )( 0, 0 ) == Approx( -1 ) );
    REQUIRE( cube.image( 1 )( 0, 0 ) == Approx( 0 ) );
    REQUIRE( cube.image( 2 )( 0, 0 ) == Approx( 4 ) );

    makeMeanSubCube( cube );
    observation.m_maskFile = "configured-mask.fits";
    observation.m_preProcess_mask = true;
    observation.m_mask.resize( 1, 2 );
    observation.m_mask << 1, 0;
    observation.preProcess_meanSub( cube );
    REQUIRE( cube.image( 0 )( 0, 1 ) == 0 );
    REQUIRE( cube.image( 1 )( 0, 1 ) == 0 );
    REQUIRE( cube.image( 2 )( 0, 1 ) == 0 );

    observation.m_preProcess_meanSubMethod = mx::improc::HCI::meanSub::imageMean;
    REQUIRE_THROWS( observation.preProcess_meanSub( cube ) );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::preProcess_meanSub( cube );
#endif
    // clang-format on
}

/// Verify HCIobservation::preProcess_pixelTSNorm uses the population RMS and remains finite at edge cases.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation pixel time-series RMS normalization", "[HCIobservation][preprocess][rms]" )
{
    OpenMPThreadGuard threads;
    HCIobservationTestHarness observation;
    observation.m_preProcess_pixelTSNormMethod = mx::improc::HCI::pixelTSNorm::rms;

    HCIobservationTestHarness::cubeT cube( 1, 1, 2 );
    cube.image( 0 )( 0, 0 ) = 0;
    cube.image( 1 )( 0, 0 ) = 2;
    observation.preProcess_pixelTSNorm( cube );
    REQUIRE( cube.image( 0 )( 0, 0 ) == Approx( 0 ) );
    REQUIRE( cube.image( 1 )( 0, 0 ) == Approx( std::sqrt( 2.0 ) ) );

    cube.resize( 1, 1, 1 );
    cube.image( 0 )( 0, 0 ) = 4;
    observation.preProcess_pixelTSNorm( cube );
    REQUIRE( cube.image( 0 )( 0, 0 ) == Approx( 1 ) );
    REQUIRE( std::isfinite( cube.image( 0 )( 0, 0 ) ) );

    cube.resize( 1, 1, 2 );
    cube.setZero();
    observation.preProcess_pixelTSNorm( cube );
    REQUIRE( cube.image( 0 )( 0, 0 ) == 0 );
    REQUIRE( cube.image( 1 )( 0, 0 ) == 0 );

    cube.image( 0 )( 0, 0 ) = 3;
    cube.image( 1 )( 0, 0 ) = 4;
    observation.m_maskFile = "configured-mask.fits";
    observation.m_preProcess_mask = true;
    observation.m_mask.resize( 1, 1 );
    observation.m_mask( 0, 0 ) = 0;
    observation.preProcess_pixelTSNorm( cube );
    REQUIRE( cube.image( 0 )( 0, 0 ) == 3 );
    REQUIRE( cube.image( 1 )( 0, 0 ) == 4 );

    observation.m_preProcess_pixelTSNormMethod = mx::improc::HCI::pixelTSNorm::rmsSigmaClipped;
    REQUIRE_THROWS( observation.preProcess_pixelTSNorm( cube ) );

    observation.m_preProcess_pixelTSNormMethod = static_cast<mx::improc::HCI::pixelTSNorm>( 99 );
    REQUIRE_THROWS( observation.preProcess_pixelTSNorm( cube ) );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::preProcess_pixelTSNorm( cube );
#endif
    // clang-format on
}

} // namespace HCIobservation_preprocess_test
} // namespace unitTest
