/** \file HCIobservation_combine_test.cpp
 * \brief Tests HCIobservation final image combination.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "HCIobservation_test_fixture.hpp"

#include <mx/improc/imageUtils.hpp>

#include <cmath>

namespace unitTest
{
namespace HCIobservation_combine_test
{

/// Verify HCIobservation::combineFinim no-op, mean, weighted-mean, sigma fallback, and validation behavior.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation final mean combinations", "[HCIobservation][combine]" )
{
    OpenMPThreadGuard threads;
    HCIobservationTestHarness observation;

    observation.m_combineMethod = mx::improc::HCI::combine::none;
    observation.combineFinim();
    REQUIRE( observation.m_finim.planes() == 0 );

    observation.m_psfsub.resize( 2 );
    observation.m_psfsub[0].resize( 1, 1, 3 );
    observation.m_psfsub[1].resize( 1, 1, 3 );
    observation.m_psfsub[0].image( 0 )( 0, 0 ) = 1;
    observation.m_psfsub[0].image( 1 )( 0, 0 ) = 3;
    observation.m_psfsub[0].image( 2 )( 0, 0 ) = 8;
    observation.m_psfsub[1].image( 0 )( 0, 0 ) = 2;
    observation.m_psfsub[1].image( 1 )( 0, 0 ) = 4;
    observation.m_psfsub[1].image( 2 )( 0, 0 ) = 6;

    observation.m_combineMethod = mx::improc::HCI::combine::mean;
    observation.combineFinim();
    REQUIRE( observation.m_finim.planes() == 2 );
    REQUIRE( observation.m_finim.image( 0 )( 0, 0 ) == Approx( 4 ) );
    REQUIRE( observation.m_finim.image( 1 )( 0, 0 ) == Approx( 4 ) );

    observation.m_comboWeights = { 0.25, 0.25, 0.5 };
    observation.combineFinim();
    REQUIRE( observation.m_finim.image( 0 )( 0, 0 ) == Approx( 5 ) );
    REQUIRE( observation.m_finim.image( 1 )( 0, 0 ) == Approx( 4.5 ) );

    observation.m_combineMethod = mx::improc::HCI::combine::sigmaMean;
    observation.m_sigmaThreshold = 0;
    observation.combineFinim();
    REQUIRE( observation.m_finim.image( 0 )( 0, 0 ) == Approx( 5 ) );

    observation.m_sigmaThreshold = 3;
    observation.m_comboWeights.clear();
    observation.m_psfsub[0].image( 0 ).setConstant( 7 );
    observation.m_psfsub[0].image( 1 ).setConstant( 7 );
    observation.m_psfsub[0].image( 2 ).setConstant( 7 );
    observation.combineFinim();
    REQUIRE( observation.m_finim.image( 0 )( 0, 0 ) == Approx( 7 ) );

    observation.m_comboWeights = { 1, 2 };
    REQUIRE_THROWS( observation.combineFinim() );

    observation.m_comboWeights.clear();
    observation.m_combineMethod = static_cast<mx::improc::HCI::combine>( 99 );
    REQUIRE_THROWS( observation.combineFinim() );

    observation.m_combineMethod = mx::improc::HCI::combine::mean;
    observation.m_psfsub[1].resize( 2, 1, 3 );
    REQUIRE_THROWS( observation.combineFinim() );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::combineFinim();
#endif
    // clang-format on
}

/// Verify HCIobservation::combineFinim applies masks and inclusive minimum-good-fraction semantics to medians.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation masked median combination", "[HCIobservation][combine][mask]" )
{
    OpenMPThreadGuard threads;
    HCIobservationTestHarness observation;
    observation.m_psfsub.resize( 1 );
    observation.m_psfsub[0].resize( 1, 2, 4 );
    observation.m_psfsub[0].image( 0 ) << 1, 10;
    observation.m_psfsub[0].image( 1 ) << 5, 20;
    observation.m_psfsub[0].image( 2 ) << 9, 30;
    observation.m_psfsub[0].image( 3 ) << 11, 40;

    observation.m_maskFile = "configured-mask.fits";
    observation.m_maskCube.resize( 1, 2, 4 );
    observation.m_maskCube.image( 0 ) << 1, 0;
    observation.m_maskCube.image( 1 ) << 1, 0;
    observation.m_maskCube.image( 2 ) << 0, 0;
    observation.m_maskCube.image( 3 ) << 0, 0;
    observation.m_minGoodFract = 0.5;
    observation.m_combineMethod = mx::improc::HCI::combine::median;

    observation.combineFinim();
    REQUIRE( observation.m_finim.image( 0 )( 0, 0 ) == Approx( 3 ) );
    REQUIRE( mx::improc::isInvalidPixel( observation.m_finim.image( 0 )( 0, 1 ) ) );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::combineFinim();
#endif
    // clang-format on
}

/// Verify HCIobservation::combineFinim rejects missing and empty PSF-subtracted inputs.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation empty combination inputs", "[HCIobservation][combine]" )
{
    HCIobservationTestHarness observation;
    observation.m_combineMethod = mx::improc::HCI::combine::mean;
    REQUIRE_THROWS( observation.combineFinim() );

    observation.m_psfsub.resize( 1 );
    REQUIRE_THROWS( observation.combineFinim() );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::combineFinim();
#endif
    // clang-format on
}

} // namespace HCIobservation_combine_test
} // namespace unitTest
