/** \file P4NegativeOptimizer_test.cpp
 * \brief Tests contrast-only negative-companion optimization for pixel-local P4 products.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "src/common/P4NegativeOptimizer.hpp"

#include <cmath>
#include <limits>

namespace unitTest
{
namespace P4NegativeOptimizer_test
{

/// Construct one fully valid local evaluation with a constant residual plane.
mx::improc::P4LocalEvaluation<float> constantEvaluation( double value, /**< [in] residual value */
                                                         int size = 11 /**< [in] square stamp size */ )
{
    mx::improc::P4LocalEvaluation<float> evaluation;
    evaluation.residual.resize( size, size, 1 );
    evaluation.residual.image( 0 ).setConstant( static_cast<float>( value ) );
    evaluation.validity.resize( size, size, 1 );
    evaluation.validity.image( 0 ).setOnes();
    evaluation.sourceRow = 0.5 * static_cast<double>( size - 1 );
    evaluation.sourceColumn = 0.5 * static_cast<double>( size - 1 );
    evaluation.elapsedSeconds = 0.25;
    evaluation.timing.geometryElapsedSeconds = 0.05;
    evaluation.timing.regressionElapsedSeconds = 0.2;
    return evaluation;
}

/// Verify p4LocalL2Merit uses the continuous aperture, requires complete validity, and rejects edge crossings.
/** \ingroup P4NegativeOptimizer_unit_tests */
TEST_CASE( "P4 local optimizer aperture merit", "[P4NegativeOptimizer][merit][validity][edge]" )
{
    auto evaluation = constantEvaluation( 2 );
    evaluation.sourceRow = 5.3;
    evaluation.sourceColumn = 5.2;
    REQUIRE( mx::improc::p4LocalL2Merit( evaluation, 0, 5 ) == Approx( 4 ) );

    evaluation.validity.image( 0 )( 5, 5 ) = 0;
    REQUIRE_THROWS( mx::improc::p4LocalL2Merit( evaluation, 0, 5 ) );
    evaluation.validity.image( 0 )( 5, 5 ) = 1;
    evaluation.residual.image( 0 )( 5, 5 ) = std::numeric_limits<float>::quiet_NaN();
    REQUIRE_THROWS( mx::improc::p4LocalL2Merit( evaluation, 0, 5 ) );

    evaluation.residual.image( 0 )( 5, 5 ) = 2;
    REQUIRE_THROWS( mx::improc::p4LocalL2Merit( evaluation, 1, 5 ) );
    REQUIRE_THROWS( mx::improc::p4LocalL2Merit( evaluation, 0, 5.3 ) );
}

/// Verify optimizeP4Contrast finds a bounded quadratic minimum and preserves dense-scan diagnostics.
/** \ingroup P4NegativeOptimizer_unit_tests */
TEST_CASE( "P4 bounded contrast optimizer", "[P4NegativeOptimizer][Brent][dense][timing]" )
{
    constexpr double expectedContrast = -0.012345;
    mx::improc::P4ContrastOptimizerConfig configuration;
    configuration.contrastLower = -0.05;
    configuration.contrastUpper = 0;
    configuration.validationSamples = 11;
    configuration.maxEvaluations = 48;
    configuration.parameterTolerance = 1e-8;
    configuration.meritTolerance = 1e-10;

    const auto evaluate = [=]( double contrast ) { return constantEvaluation( contrast - expectedContrast ); };
    const auto result = mx::improc::optimizeP4Contrast( configuration, evaluate, 0, 1 );
    REQUIRE( result.converged );
    REQUIRE( result.denseAgreement );
    REQUIRE( result.status == "converged" );
    REQUIRE( result.bestContrast == Approx( expectedContrast ).margin( 2e-8 ) );
    REQUIRE( result.bestMerit == Approx( 0 ).margin( 1e-14 ) );
    REQUIRE( result.evaluationCount == result.samples.size() );
    REQUIRE( result.evaluationCount > configuration.validationSamples );
    REQUIRE( result.evaluationCount <= configuration.maxEvaluations );
    REQUIRE( result.evaluationElapsedSeconds == Approx( 0.25 * result.evaluationCount ) );
    REQUIRE( result.timing.geometryElapsedSeconds == Approx( 0.05 * result.evaluationCount ) );
    REQUIRE( result.timing.regressionElapsedSeconds == Approx( 0.2 * result.evaluationCount ) );
    for( std::size_t index = 0; index < configuration.validationSamples; ++index )
    {
        REQUIRE( result.samples[index].denseValidation );
    }
}

/// Verify optimizeP4Contrast rejects inconsistent bounds, scan budgets, tolerances, and callbacks.
/** \ingroup P4NegativeOptimizer_unit_tests */
TEST_CASE( "P4 contrast optimizer validation", "[P4NegativeOptimizer][validation][budget]" )
{
    const auto evaluate = []( double contrast ) { return constantEvaluation( contrast ); };
    mx::improc::P4ContrastOptimizerConfig configuration;

    configuration.contrastUpper = 0.1;
    REQUIRE_THROWS( mx::improc::optimizeP4Contrast( configuration, evaluate, 0, 1 ) );
    configuration = {};
    configuration.validationSamples = 2;
    REQUIRE_THROWS( mx::improc::optimizeP4Contrast( configuration, evaluate, 0, 1 ) );
    configuration = {};
    configuration.maxEvaluations = configuration.validationSamples + 2;
    REQUIRE_THROWS( mx::improc::optimizeP4Contrast( configuration, evaluate, 0, 1 ) );
    configuration = {};
    configuration.parameterTolerance = 0;
    REQUIRE_THROWS( mx::improc::optimizeP4Contrast( configuration, evaluate, 0, 1 ) );
    configuration = {};
    REQUIRE_THROWS( mx::improc::optimizeP4Contrast( configuration, {}, 0, 1 ) );
}

} // namespace P4NegativeOptimizer_test
} // namespace unitTest
