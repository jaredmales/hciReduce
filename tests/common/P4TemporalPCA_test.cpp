/** \file P4TemporalPCA_test.cpp
 * \brief Tests gap-held-out time-domain PCA prediction for Pixel Prediction Post-Processing.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "src/common/P4TemporalPCA.hpp"

#include <cmath>

namespace unitTest
{
namespace P4TemporalPCA_test
{

/** \cond P4TemporalPCA_test_harness */
using temporalPCAT = mx::improc::P4TemporalPCA;
using centeringT = mx::improc::P4TemporalPCACentering;
using resultT = mx::improc::P4TemporalPCAResult;

/// Construct three reference-pixel series with a two-dimensional centered temporal subspace.
temporalPCAT::matrixT referenceSeries( Eigen::Index samples /**< [in] positive time-series sample count */ )
{
    temporalPCAT::matrixT references( 3, samples );
    for( Eigen::Index sample = 0; sample < samples; ++sample )
    {
        const double time = static_cast<double>( sample ) - 4.0;
        references( 0, sample ) = 5.0 + time;
        references( 1, sample ) = -2.0 + time * time;
        references( 2, sample ) = 3.0 + 2.0 * time - 0.5 * time * time;
    }
    return references;
}

/// Construct one target series exactly represented by an intercept plus the reference temporal subspace.
temporalPCAT::vectorT targetSeries( Eigen::Index samples /**< [in] positive time-series sample count */ )
{
    temporalPCAT::vectorT target( samples );
    for( Eigen::Index sample = 0; sample < samples; ++sample )
    {
        const double time = static_cast<double>( sample ) - 4.0;
        target( sample ) = 11.0 + 3.0 * time - 0.5 * time * time;
    }
    return target;
}
/** \endcond */

/// Verify P4TemporalPCA predicts a central series from a two-mode temporal basis outside every configured gap.
/** This exercises mx::improc::P4TemporalPCA::configure() and mx::improc::P4TemporalPCA::predict().
 * \ingroup P4TemporalPCA_unit_tests
 */
TEST_CASE( "P4 temporal PCA exactly predicts a gap-held-out low-rank central series", "[P4TemporalPCA][prediction]" )
{
    constexpr Eigen::Index samples = 9;
    temporalPCAT predictor;
    temporalPCAT::workspaceT workspace;
    predictor.configure( referenceSeries( samples ), 2, centeringT::pixelMean, 1e-12, 1, workspace );

    const temporalPCAT::vectorT target = targetSeries( samples );
    resultT result;
    predictor.predict( result, target );

    REQUIRE( predictor.sampleCount() == samples );
    REQUIRE( predictor.requestedModes() == 2 );
    REQUIRE( result.numericalRank == 2 );
    for( Eigen::Index sample = 0; sample < samples; ++sample )
    {
        REQUIRE( result.validity( sample ) == 1 );
        REQUIRE( result.predictions( sample ) == Approx( target( sample ) ).margin( 1e-10 ) );
        REQUIRE( result.residuals( sample ) == Approx( 0.0 ).margin( 1e-10 ) );
    }
}

/// Verify P4TemporalPCA never fits a target prediction with samples inside its no-wrap temporal gap.
/** This exercises mx::improc::P4TemporalPCA::configure() and mx::improc::P4TemporalPCA::predict().
 * \ingroup P4TemporalPCA_unit_tests
 */
TEST_CASE( "P4 temporal PCA excludes the complete configured prediction gap", "[P4TemporalPCA][gap]" )
{
    constexpr Eigen::Index samples = 9;
    constexpr Eigen::Index target = 4;
    temporalPCAT predictor;
    temporalPCAT::workspaceT workspace;
    predictor.configure( referenceSeries( samples ), 2, centeringT::pixelMean, 1e-12, 1, workspace );

    const temporalPCAT::vectorT unmodified = targetSeries( samples );
    resultT baseline;
    predictor.predict( baseline, unmodified );

    temporalPCAT::vectorT modified = unmodified;
    modified( target - 1 ) += 1000;
    modified( target ) -= 500;
    modified( target + 1 ) += 250;
    resultT result;
    predictor.predict( result, modified );

    REQUIRE( result.validity( target ) == 1 );
    REQUIRE( result.predictions( target ) == Approx( baseline.predictions( target ) ).margin( 1e-10 ) );
    REQUIRE( result.predictions( target ) == Approx( unmodified( target ) ).margin( 1e-10 ) );
}

/// Verify P4TemporalPCA marks samples unsupported when their no-wrap temporal gaps leave too little fit support.
/** This exercises mx::improc::P4TemporalPCA::configure() and mx::improc::P4TemporalPCA::predict().
 * \ingroup P4TemporalPCA_unit_tests
 */
TEST_CASE( "P4 temporal PCA marks oversized temporal gaps invalid", "[P4TemporalPCA][gap][validity]" )
{
    constexpr Eigen::Index samples = 9;
    temporalPCAT predictor;
    temporalPCAT::workspaceT workspace;
    predictor.configure( referenceSeries( samples ), 1, centeringT::pixelMean, 1e-12, 4, workspace );

    resultT result;
    predictor.predict( result, targetSeries( samples ) );

    REQUIRE( result.validity( 0 ) == 1 );
    REQUIRE( result.validity( samples / 2 ) == 0 );
    REQUIRE( std::isnan( result.predictions( samples / 2 ) ) );
    REQUIRE( std::isnan( result.residuals( samples / 2 ) ) );
}

/// Verify P4TemporalPCA invalidates every prediction when the requested temporal basis is rank-deficient.
/** This exercises mx::improc::P4TemporalPCA::configure() and mx::improc::P4TemporalPCA::predict().
 * \ingroup P4TemporalPCA_unit_tests
 */
TEST_CASE( "P4 temporal PCA rejects a rank-deficient requested basis", "[P4TemporalPCA][rank][validity]" )
{
    constexpr Eigen::Index samples = 9;
    temporalPCAT predictor;
    temporalPCAT::workspaceT workspace;
    predictor.configure( referenceSeries( samples ), 3, centeringT::pixelMean, 1e-12, 0, workspace );

    resultT result;
    predictor.predict( result, targetSeries( samples ) );

    REQUIRE( result.numericalRank == 2 );
    REQUIRE( ( result.validity == 0 ).all() );
    for( Eigen::Index sample = 0; sample < samples; ++sample )
    {
        REQUIRE( std::isnan( result.predictions( sample ) ) );
        REQUIRE( std::isnan( result.residuals( sample ) ) );
    }
}

} // namespace P4TemporalPCA_test
} // namespace unitTest
