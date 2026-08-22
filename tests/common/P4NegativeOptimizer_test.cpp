/** \file P4NegativeOptimizer_test.cpp
 * \brief Tests negative-companion optimization for pixel-local P4 products.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "src/common/P4NegativeOptimizer.hpp"

#include <algorithm>
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

/// Verify an explicit full-image aperture retains identical sky pixels as the trial source and stamp origin move.
/** \ingroup P4NegativeOptimizer_unit_tests */
TEST_CASE( "P4 fixed sky aperture merit", "[P4NegativeOptimizer][merit][position][aperture]" )
{
    auto evaluation = constantEvaluation( 2, 13 );
    evaluation.originRow = 100;
    evaluation.originColumn = 100;
    evaluation.sourceRow = 106.1;
    evaluation.sourceColumn = 106;
    evaluation.residual.image( 0 )( 1, 6 ) = 10;

    constexpr double expectedFixedMerit = 420.0 / 81.0;
    REQUIRE( mx::improc::p4LocalL2Merit( evaluation, 0, 5 ) == Approx( 4 ) );
    REQUIRE( mx::improc::p4LocalL2Merit( evaluation, 0, 5, 106, 106 ) == Approx( expectedFixedMerit ) );

    auto shifted = constantEvaluation( 2, 13 );
    shifted.originRow = 101;
    shifted.originColumn = 101;
    shifted.sourceRow = 105.8;
    shifted.sourceColumn = 106.2;
    shifted.residual.image( 0 )( 0, 5 ) = 10;
    REQUIRE( mx::improc::p4LocalL2Merit( shifted, 0, 5, 106, 106 ) == Approx( expectedFixedMerit ) );
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

/// Verify the polar/Cartesian conversions preserve the P4 center and position-angle convention.
/** \ingroup P4NegativeOptimizer_unit_tests */
TEST_CASE( "P4 optimizer Cartesian coordinates", "[P4NegativeOptimizer][coordinates][center]" )
{
    const mx::improc::P4LocalTrial initial{ 12.7, 260.1, -0.004 };
    const auto offset = mx::improc::p4TrialCartesianOffset( initial );
    const auto roundTrip = mx::improc::p4CartesianOffsetTrial( offset, initial.contrast );
    REQUIRE( roundTrip.separation == Approx( initial.separation ).epsilon( 1e-12 ) );
    REQUIRE( roundTrip.positionAngle == Approx( initial.positionAngle ).epsilon( 1e-12 ) );
    REQUIRE( roundTrip.contrast == initial.contrast );

    const auto north = mx::improc::p4TrialCartesianOffset( { 5, 0, 0 } );
    REQUIRE( north.row == Approx( 0 ).margin( 1e-15 ) );
    REQUIRE( north.column == Approx( 5 ) );
    const auto east = mx::improc::p4TrialCartesianOffset( { 5, 90, 0 } );
    REQUIRE( east.row == Approx( -5 ) );
    REQUIRE( east.column == Approx( 0 ).margin( 1e-15 ) );
}

/// Verify the full optimizer recovers a bounded Cartesian position and contrast on a smooth synthetic merit.
/** \ingroup P4NegativeOptimizer_unit_tests */
TEST_CASE( "P4 joint position contrast optimizer", "[P4NegativeOptimizer][position][contrast][simplex][dense]" )
{
    const mx::improc::P4LocalTrial initial{ 10, 0, -0.02 };
    const auto initialOffset = mx::improc::p4TrialCartesianOffset( initial );
    constexpr double expectedRowDelta = 0.22;
    constexpr double expectedColumnDelta = 0.31;
    constexpr double expectedContrast = -0.0123;
    constexpr double apertureRow = 100;
    constexpr double apertureColumn = 100;

    mx::improc::P4ContrastOptimizerConfig configuration;
    configuration.contrastLower = -0.05;
    configuration.contrastUpper = 0;
    configuration.validationSamples = 7;
    configuration.maxEvaluations = 192;
    configuration.parameterTolerance = 1e-4;
    configuration.positionTolerance = 1e-4;
    configuration.meritTolerance = 1e-5;

    const auto evaluate = [=]( const mx::improc::P4LocalTrial &trial )
    {
        const auto offset = mx::improc::p4TrialCartesianOffset( trial );
        const double rowError = offset.row - initialOffset.row - expectedRowDelta;
        const double columnError = offset.column - initialOffset.column - expectedColumnDelta;
        const double contrastError = ( trial.contrast - expectedContrast ) / 0.01;
        const double objective =
            0.05 + rowError * rowError + 1.5 * columnError * columnError + contrastError * contrastError;
        auto evaluation = constantEvaluation( std::sqrt( objective ), 13 );
        evaluation.sourceRow = apertureRow + offset.row - initialOffset.row;
        evaluation.sourceColumn = apertureColumn + offset.column - initialOffset.column;
        evaluation.originRow = static_cast<int>( std::floor( evaluation.sourceRow + 0.5 ) ) - 6;
        evaluation.originColumn = static_cast<int>( std::floor( evaluation.sourceColumn + 0.5 ) ) - 6;
        const int boundaryRow = static_cast<int>( apertureRow - 5 ) - evaluation.originRow;
        const int boundaryColumn = static_cast<int>( apertureColumn ) - evaluation.originColumn;
        evaluation.residual.image( 0 )( boundaryRow, boundaryColumn ) =
            static_cast<float>( std::sqrt( objective + 4 ) );
        return evaluation;
    };
    const auto result = mx::improc::optimizeP4PositionContrast( configuration, initial, 1, evaluate, 0, 5 );
    REQUIRE( result.converged );
    REQUIRE( result.positionConverged );
    REQUIRE( result.contrastConverged );
    REQUIRE( result.denseAgreement );
    REQUIRE( result.status == "converged" );
    REQUIRE( result.bestRowDelta == Approx( expectedRowDelta ).margin( 5e-4 ) );
    REQUIRE( result.bestColumnDelta == Approx( expectedColumnDelta ).margin( 5e-4 ) );
    REQUIRE( result.bestTrial.contrast == Approx( expectedContrast ).margin( 2e-5 ) );
    REQUIRE( result.apertureRow == apertureRow );
    REQUIRE( result.apertureColumn == apertureColumn );
    REQUIRE( result.bestMerit == Approx( 0.05 + 4.0 / 81.0 ).margin( 1e-6 ) );
    REQUIRE( result.evaluationCount == result.samples.size() );
    REQUIRE( result.evaluationCount <= configuration.maxEvaluations );
    REQUIRE( result.evaluationCount >= mx::improc::p4PositionContrastMinimumEvaluations( 7 ) - 4 );
    REQUIRE( std::any_of( result.samples.begin(),
                          result.samples.end(),
                          []( const auto &sample ) { return sample.stage == "simplex"; } ) );
    REQUIRE( std::any_of( result.samples.begin(),
                          result.samples.end(),
                          []( const auto &sample ) { return sample.stage == "final-dense"; } ) );
}

/// Verify ordinary jackknife covariance uses the delete-one-block multiplier and Cartesian parameter order.
/** \ingroup P4NegativeOptimizer_unit_tests */
TEST_CASE( "P4 jackknife covariance", "[P4NegativeOptimizer][jackknife][covariance][uncertainty]" )
{
    const std::vector<mx::improc::P4JackknifeEstimate> estimates{ { 1, 2, -1 }, { 2, 4, -2 }, { 3, 6, -3 } };
    const auto joint = mx::improc::p4JackknifeStatistics( estimates, true );
    REQUIRE( joint.blockCount == 3 );
    REQUIRE( joint.mean.rowDelta == Approx( 2 ) );
    REQUIRE( joint.mean.columnDelta == Approx( 4 ) );
    REQUIRE( joint.mean.contrast == Approx( -2 ) );
    REQUIRE( joint.covariance[0] == Approx( 4.0 / 3.0 ) );
    REQUIRE( joint.covariance[1] == Approx( 8.0 / 3.0 ) );
    REQUIRE( joint.covariance[2] == Approx( -4.0 / 3.0 ) );
    REQUIRE( joint.covariance[4] == Approx( 16.0 / 3.0 ) );
    REQUIRE( joint.covariance[5] == Approx( -8.0 / 3.0 ) );
    REQUIRE( joint.covariance[8] == Approx( 4.0 / 3.0 ) );
    REQUIRE( joint.rowStandardError == Approx( std::sqrt( 4.0 / 3.0 ) ) );
    REQUIRE( joint.columnStandardError == Approx( std::sqrt( 16.0 / 3.0 ) ) );
    REQUIRE( joint.contrastStandardError == Approx( std::sqrt( 4.0 / 3.0 ) ) );

    const auto contrastOnly = mx::improc::p4JackknifeStatistics( estimates, false );
    REQUIRE( contrastOnly.mean.rowDelta == 0 );
    REQUIRE( contrastOnly.mean.columnDelta == 0 );
    REQUIRE( contrastOnly.covariance[0] == 0 );
    REQUIRE( contrastOnly.covariance[4] == 0 );
    REQUIRE( contrastOnly.covariance[8] == Approx( 4.0 / 3.0 ) );

    REQUIRE_THROWS( mx::improc::p4JackknifeStatistics( { estimates.front() }, true ) );
    auto nonFinite = estimates;
    nonFinite.back().contrast = std::numeric_limits<double>::quiet_NaN();
    REQUIRE_THROWS( mx::improc::p4JackknifeStatistics( nonFinite, false ) );
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

    const mx::improc::P4LocalTrial initial{ 5, 20, -0.01 };
    const auto jointEvaluate = []( const mx::improc::P4LocalTrial &trial )
    { return constantEvaluation( trial.contrast ); };
    configuration = {};
    configuration.maxEvaluations =
        mx::improc::p4PositionContrastMinimumEvaluations( configuration.validationSamples ) - 1;
    REQUIRE_THROWS( mx::improc::optimizeP4PositionContrast( configuration, initial, 1, jointEvaluate, 0, 1 ) );
    configuration.maxEvaluations = mx::improc::p4PositionContrastMinimumEvaluations( configuration.validationSamples );
    REQUIRE_THROWS( mx::improc::optimizeP4PositionContrast( configuration, initial, 0, jointEvaluate, 0, 1 ) );
    configuration.positionTolerance = 0;
    REQUIRE_THROWS( mx::improc::optimizeP4PositionContrast( configuration, initial, 1, jointEvaluate, 0, 1 ) );
}

} // namespace P4NegativeOptimizer_test
} // namespace unitTest
