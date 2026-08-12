/** \file KLIPreduction_test.cpp
 * \brief Tests KLIP reduction centering and normalization behavior.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "src/common/ADIDerotator.hpp"
#include "src/common/KLIPreduction.hpp"

namespace unitTest
{
namespace KLIPreduction_test
{

/// \cond KLIPreduction_test_harness
using reductionT =
    mx::improc::KLIPreduction<float, mx::improc::ADIDerotator<float, mx::verbose::vv>, double, mx::verbose::vv>;
/// \endcond

/// Verify KLIPreduction::meanSubtract leaves both cubes unchanged for `none` while calculating reference norms.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP mean subtraction none", "[KLIPreduction][meanSubtract][none]" )
{
    reductionT reduction;
    reduction.m_meanSubMethod = mx::improc::HCI::meanSub::none;
    reduction.m_pixelTSNormMethod = mx::improc::HCI::pixelTSNorm::none;

    mx::improc::eigenCube<float> references( 2, 2, 2 );
    references.image( 0 ) << 1, 3, 2, 4;
    references.image( 1 ) << 10, 14, 12, 16;
    mx::improc::eigenCube<float> targets( 2, 2, 1 );
    targets.image( 0 ) << 20, 22, 21, 23;
    const mx::improc::eigenCube<float> originalReferences = references;
    const mx::improc::eigenCube<float> originalTargets = targets;
    reductionT::imageT mask;
    std::vector<float> norms;

    reduction.meanSubtract( references, targets, mask, norms );

    REQUIRE( references.image( 0 ).isApprox( originalReferences.image( 0 ) ) );
    REQUIRE( references.image( 1 ).isApprox( originalReferences.image( 1 ) ) );
    REQUIRE( targets.image( 0 ).isApprox( originalTargets.image( 0 ) ) );
    REQUIRE( norms.size() == 2 );
    REQUIRE( norms[0] == Approx( 5 ) );
    REQUIRE( norms[1] == Approx( 20 ) );
}

/// Verify KLIPreduction::meanSubtract rejects unimplemented and invalid methods before mutating either cube.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP mean subtraction method validation", "[KLIPreduction][meanSubtract][validation]" )
{
    reductionT reduction;
    mx::improc::eigenCube<float> references( 1, 2, 1 );
    references.image( 0 ) << 1, 2;
    mx::improc::eigenCube<float> targets = references;
    const mx::improc::eigenCube<float> original = references;
    reductionT::imageT mask;
    std::vector<float> norms;

    reduction.m_meanSubMethod = mx::improc::HCI::meanSub::imageMode;
    REQUIRE_THROWS( reduction.meanSubtract( references, targets, mask, norms ) );
    REQUIRE( references.image( 0 ).isApprox( original.image( 0 ) ) );
    REQUIRE( targets.image( 0 ).isApprox( original.image( 0 ) ) );
    REQUIRE( norms.empty() );

    reduction.m_meanSubMethod = static_cast<mx::improc::HCI::meanSub>( 99 );
    REQUIRE_THROWS( reduction.meanSubtract( references, targets, mask, norms ) );
    REQUIRE( references.image( 0 ).isApprox( original.image( 0 ) ) );
    REQUIRE( targets.image( 0 ).isApprox( original.image( 0 ) ) );
    REQUIRE( norms.empty() );
}

} // namespace KLIPreduction_test
} // namespace unitTest
