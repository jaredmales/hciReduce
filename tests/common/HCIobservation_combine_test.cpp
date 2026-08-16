/** \file HCIobservation_combine_test.cpp
 * \brief Tests HCIobservation final image combination.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "HCIobservation_test_fixture.hpp"

#include <mx/improc/imageUtils.hpp>

#include <cmath>
#include <limits>

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

/// Verify HCIobservation::combineFinim exercises every weighted and masked mean and sigma-mean combination.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation weighted masked combinations", "[HCIobservation][combine][mask][weights]" )
{
    OpenMPThreadGuard threads;
    HCIobservationTestHarness observation;
    observation.m_psfsub.resize( 1 );
    observation.m_psfsub[0].resize( 1, 1, 3 );
    observation.m_psfsub[0].image( 0 )( 0, 0 ) = 1;
    observation.m_psfsub[0].image( 1 )( 0, 0 ) = 3;
    observation.m_psfsub[0].image( 2 )( 0, 0 ) = 5;
    observation.m_maskCube.resize( 1, 1, 3 );
    observation.m_maskCube.cube().setOnes();
    observation.m_minGoodFract = 1;

    observation.m_maskFile = "configured-mask.fits";
    observation.m_combineMethod = mx::improc::HCI::combine::mean;
    observation.m_comboWeights = { 0.25F, 0.25F, 0.5F };
    observation.combineFinim();
    REQUIRE( observation.m_finim.image( 0 )( 0, 0 ) == Approx( 3.5F ) );

    observation.m_comboWeights.clear();
    observation.combineFinim();
    REQUIRE( observation.m_finim.image( 0 )( 0, 0 ) == Approx( 3 ) );

    observation.m_combineMethod = mx::improc::HCI::combine::sigmaMean;
    observation.m_sigmaThreshold = 0;
    observation.m_comboWeights = { 0.25F, 0.25F, 0.5F };
    observation.combineFinim();
    REQUIRE( observation.m_finim.image( 0 )( 0, 0 ) == Approx( 3.5F ) );

    observation.m_comboWeights.clear();
    observation.combineFinim();
    REQUIRE( observation.m_finim.image( 0 )( 0, 0 ) == Approx( 3 ) );

    observation.m_maskFile.clear();
    observation.combineFinim();
    REQUIRE( observation.m_finim.image( 0 )( 0, 0 ) == Approx( 3 ) );

    observation.m_sigmaThreshold = 10;
    observation.m_comboWeights = { 0.25F, 0.25F, 0.5F };
    observation.m_maskFile = "configured-mask.fits";
    observation.combineFinim();
    REQUIRE( observation.m_finim.image( 0 )( 0, 0 ) == Approx( 3.5F ) );

    observation.m_maskFile.clear();
    observation.combineFinim();
    REQUIRE( observation.m_finim.image( 0 )( 0, 0 ) == Approx( 3.5F ) );

    observation.m_comboWeights.clear();
    observation.m_maskFile = "configured-mask.fits";
    observation.combineFinim();
    REQUIRE( observation.m_finim.image( 0 )( 0, 0 ) == Approx( 3 ) );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::combineFinim();
#endif
    // clang-format on
}

/// Verify HCIobservation::combineFinim uses each reduction's authoritative validity cube for every combination mode.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation per-reduction validity combination", "[HCIobservation][combine][validity][mask][weights]" )
{
    OpenMPThreadGuard threads;
    HCIobservationTestHarness observation;
    observation.m_psfsub.resize( 2 );
    observation.m_psfsubValidity.resize( 2 );
    for( size_t reduction = 0; reduction < 2; ++reduction )
    {
        observation.m_psfsub[reduction].resize( 1, 2, 3 );
        observation.m_psfsubValidity[reduction].resize( 1, 2, 3 );
    }

    observation.m_psfsub[0].image( 0 ) << 1, 2;
    observation.m_psfsub[0].image( 1 ) << std::numeric_limits<float>::quiet_NaN(), 4;
    observation.m_psfsub[0].image( 2 ) << 5, 6;
    observation.m_psfsubValidity[0].image( 0 ) << 1, 0;
    observation.m_psfsubValidity[0].image( 1 ) << 0, 0;
    observation.m_psfsubValidity[0].image( 2 ) << 1, 0;

    observation.m_psfsub[1].image( 0 ) << 2, 10;
    observation.m_psfsub[1].image( 1 ) << 4, 20;
    observation.m_psfsub[1].image( 2 ) << 8, 30;
    observation.m_psfsubValidity[1].image( 0 ) << 0, 1;
    observation.m_psfsubValidity[1].image( 1 ) << 1, 1;
    observation.m_psfsubValidity[1].image( 2 ) << 1, 0;

    observation.m_maskFile = "legacy-mask-must-not-win.fits";
    observation.m_maskCube.resize( 1, 2, 3 );
    observation.m_maskCube.cube().setZero();
    observation.m_minGoodFract = 0;

    observation.m_combineMethod = mx::improc::HCI::combine::median;
    observation.combineFinim();
    REQUIRE( observation.m_finim.image( 0 )( 0, 0 ) == Approx( 3 ) );
    REQUIRE( mx::improc::isInvalidPixel( observation.m_finim.image( 0 )( 0, 1 ) ) );
    REQUIRE( observation.m_finim.image( 1 )( 0, 0 ) == Approx( 6 ) );
    REQUIRE( observation.m_finim.image( 1 )( 0, 1 ) == Approx( 15 ) );
    REQUIRE( observation.m_psfsub[0].image( 1 )( 0, 0 ) == 0 );

    observation.m_combineMethod = mx::improc::HCI::combine::mean;
    observation.combineFinim();
    REQUIRE( observation.m_finim.image( 0 )( 0, 0 ) == Approx( 3 ) );
    REQUIRE( observation.m_finim.image( 1 )( 0, 0 ) == Approx( 6 ) );

    observation.m_comboWeights = { 0.25F, 0.25F, 0.5F };
    observation.combineFinim();
    REQUIRE( observation.m_finim.image( 0 )( 0, 0 ) == Approx( 11.0F / 3.0F ) );
    REQUIRE( observation.m_finim.image( 1 )( 0, 0 ) == Approx( 20.0F / 3.0F ) );

    observation.m_combineMethod = mx::improc::HCI::combine::sigmaMean;
    observation.m_sigmaThreshold = 0;
    observation.combineFinim();
    REQUIRE( observation.m_finim.image( 0 )( 0, 0 ) == Approx( 11.0F / 3.0F ) );

    observation.m_comboWeights.clear();
    observation.combineFinim();
    REQUIRE( observation.m_finim.image( 1 )( 0, 1 ) == Approx( 15 ) );

    observation.m_sigmaThreshold = 10;
    observation.m_comboWeights = { 0.25F, 0.25F, 0.5F };
    observation.combineFinim();
    REQUIRE( observation.m_finim.image( 1 )( 0, 0 ) == Approx( 20.0F / 3.0F ) );

    observation.m_comboWeights.clear();
    observation.combineFinim();
    REQUIRE( observation.m_finim.image( 0 )( 0, 0 ) == Approx( 3 ) );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::combineFinim();
#endif
    // clang-format on
}

/// Verify HCIobservation::validatePSFSubValidity and combineFinim reject malformed state before mutating science.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation validity validation", "[HCIobservation][combine][validity][validation]" )
{
    HCIobservationTestHarness observation;
    REQUIRE_NOTHROW( observation.validatePSFSubValidity() );

    observation.m_combineMethod = mx::improc::HCI::combine::mean;
    observation.m_psfsub.resize( 1 );
    observation.m_psfsub[0].resize( 1, 1, 1 );
    observation.m_psfsub[0].cube().setOnes();

    observation.m_psfsubValidity.resize( 2 );
    REQUIRE_THROWS( observation.combineFinim() );

    observation.m_psfsubValidity.resize( 1 );
    REQUIRE_THROWS( observation.combineFinim() );

    observation.m_psfsubValidity[0].resize( 1, 1, 1 );
    observation.m_psfsubValidity[0].cube().setConstant( 0.5F );
    REQUIRE_THROWS( observation.combineFinim() );
    REQUIRE( observation.m_psfsub[0].image( 0 )( 0, 0 ) == 1 );

    observation.m_psfsubValidity[0].cube().setConstant( std::numeric_limits<float>::quiet_NaN() );
    REQUIRE_THROWS( observation.combineFinim() );

    observation.m_psfsubValidity[0].cube().setOnes();
    observation.m_psfsub[0].cube().setConstant( std::numeric_limits<float>::quiet_NaN() );
    REQUIRE_THROWS( observation.combineFinim() );

    observation.m_psfsubValidity[0].cube().setZero();
    REQUIRE_NOTHROW( observation.combineFinim() );
    REQUIRE( observation.m_psfsub[0].image( 0 )( 0, 0 ) == 0 );
    REQUIRE( mx::improc::isInvalidPixel( observation.m_finim.image( 0 )( 0, 0 ) ) );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::validatePSFSubValidity();
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
