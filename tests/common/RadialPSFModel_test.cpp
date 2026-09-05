/** \file RadialPSFModel_test.cpp
 * \brief Tests sparse polar PSF response averaging and nearest-radius interpolation.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "src/common/RadialPSFModel.hpp"

#include <cmath>
#include <limits>
#include <numbers>
#include <vector>

namespace unitTest
{
namespace RadialPSFModel_test
{

/** \cond RadialPSFModel_test_harness */
using modelT = mx::improc::RadialPSFModel;
using imageT = modelT::imageT;
using validityT = modelT::validityT;

/// Construct an asymmetric canonical response whose orientation is unambiguous.
imageT asymmetricResponse( int size /**< [in] positive odd stamp size */ )
{
    imageT response( size, size );
    const double center = 0.5 * static_cast<double>( size - 1 );
    for( int column = 0; column < size; ++column )
    {
        for( int row = 0; row < size; ++row )
        {
            const double deltaRow = static_cast<double>( row ) - center;
            const double deltaColumn = static_cast<double>( column ) - center;
            response( row, column ) =
                static_cast<float>( std::exp( -0.16 * ( deltaRow * deltaRow + 0.7 * deltaColumn * deltaColumn ) ) *
                                    ( 1.0 + 0.08 * deltaRow - 0.05 * deltaColumn + 0.015 * deltaRow * deltaColumn ) );
        }
    }
    return response;
}

/** \endcond */

/// Verify polar sample selection is deterministic and removes repeated integer-pixel locations within one radius.
/** This exercises mx::improc::RadialPSFModel::selectSamples() with shuffled candidates, exact ties, and more requested
 * angles than unique pixels on the inner circle.
 * \ingroup RadialPSFModel_unit_tests
 */
TEST_CASE( "Radial PSF sample selection is deterministic", "[RadialPSFModel][selection]" )
{
    const std::vector<mx::improc::RadialPSFSource> sources{ { 8, 5, 7 },
                                                            { 2, 5, 3 },
                                                            { 7, 7, 5 },
                                                            { 1, 3, 5 },
                                                            { 4, 5, 5 },
                                                            { 3, 5, 6 } };
    const std::vector<mx::improc::RadialPSFSample> samples = modelT::selectSamples( sources, 5, 5, { 0, 2 }, 8 );

    REQUIRE( samples.size() == 6 );
    REQUIRE( samples[0].radiusIndex == 0 );
    REQUIRE( samples[0].sourceIndex == 4 );
    REQUIRE( samples[0].radius == Approx( 0 ) );
    for( std::size_t index = 1; index < samples.size() - 1; ++index )
    {
        REQUIRE( samples[index].radiusIndex == 1 );
    }
    REQUIRE( samples[1].sourceIndex == 8 );
    REQUIRE( samples[2].sourceIndex == 3 );
    REQUIRE( samples[3].sourceIndex == 7 );
    REQUIRE( samples[4].sourceIndex == 2 );
    REQUIRE( samples[5].sourceIndex == 1 );
    REQUIRE( samples[2].radius == Approx( 1 ) );

    REQUIRE_THROWS_AS( modelT::selectSamples( {}, 5, 5, { 2 }, 4 ), std::invalid_argument );
    REQUIRE_THROWS_AS( modelT::selectSamples( sources, 5, 5, { 2 }, 0 ), std::invalid_argument );
}

/// Verify common-angle averaging recovers an azimuth-independent response and carries invalid support explicitly.
/** This exercises mx::improc::RadialPSFModel::rotate(), mx::improc::RadialPSFModel::fit(), and
 * mx::improc::RadialPSFModel::response() for multiple source angles.
 * \ingroup RadialPSFModel_unit_tests
 */
TEST_CASE( "Radial PSF responses align and average at a common angle", "[RadialPSFModel][rotation][average]" )
{
    constexpr int size = 13;
    const imageT canonical = asymmetricResponse( size );
    const validityT complete = validityT::Ones( size, size );
    imageT directionInput = imageT::Zero( size, size );
    directionInput( size / 2, size / 2 + 2 ) = 1;
    imageT directionOutput;
    validityT directionValidity;
    modelT::rotate( directionOutput, directionValidity, directionInput, complete, -0.5 * std::numbers::pi );
    REQUIRE( directionOutput( size / 2 + 2, size / 2 ) == Approx( 1 ) );

    const std::vector<double> angles{ 0.0, 0.5 * std::numbers::pi, std::numbers::pi, -0.5 * std::numbers::pi };
    std::vector<imageT> measured( angles.size() );
    std::vector<validityT> measuredValidity( angles.size() );
    std::vector<mx::improc::RadialPSFSample> samples;
    for( std::size_t index = 0; index < angles.size(); ++index )
    {
        modelT::rotate( measured[index], measuredValidity[index], canonical, complete, -angles[index] );
        samples.push_back( { index, 0, 6, angles[index] } );
    }
    measuredValidity[1]( size / 2, size / 2 + 2 ) = 0;

    modelT model( { 6 }, size, size );
    model.fit( measured, measuredValidity, samples );
    REQUIRE( model.sampleCount( 0 ) == angles.size() );
    REQUIRE( model.canonicalValidity( 0 )( size / 2, size / 2 ) == 1 );
    for( int column = 2; column < size - 2; ++column )
    {
        for( int row = 2; row < size - 2; ++row )
        {
            if( model.canonicalValidity( 0 )( row, column ) != 0 )
            {
                REQUIRE( model.canonicalResponse( 0 )( row, column ) ==
                         Approx( canonical( row, column ) ).margin( 2e-6 ) );
            }
        }
    }

    imageT restored;
    validityT restoredValidity;
    model.response( restored, restoredValidity, 5.9, 0.5 * std::numbers::pi );
    for( int column = 3; column < size - 3; ++column )
    {
        for( int row = 3; row < size - 3; ++row )
        {
            if( restoredValidity( row, column ) != 0 && measuredValidity[1]( row, column ) != 0 )
            {
                REQUIRE( restored( row, column ) == Approx( measured[1]( row, column ) ).margin( 3e-6 ) );
            }
        }
    }
}

/// Verify nearest-radius interpolation uses the closest canonical radial response and resolves exact ties inward.
/** This exercises mx::improc::RadialPSFModel::fit() and mx::improc::RadialPSFModel::response() with three radial bins.
 * \ingroup RadialPSFModel_unit_tests
 */
TEST_CASE( "Radial PSF interpolation selects the nearest radius", "[RadialPSFModel][nearest]" )
{
    constexpr int size = 7;
    const validityT validity = validityT::Ones( size, size );
    std::vector<imageT> responses;
    std::vector<validityT> validities;
    std::vector<mx::improc::RadialPSFSample> samples;
    for( std::size_t radiusIndex = 0; radiusIndex < 3; ++radiusIndex )
    {
        responses.push_back( imageT::Constant( size, size, static_cast<float>( radiusIndex + 1 ) ) );
        validities.push_back( validity );
        samples.push_back( { radiusIndex, radiusIndex, 4.0 + 4.0 * radiusIndex, 0 } );
    }
    modelT model( { 4, 8, 12 }, size, size );
    model.fit( responses, validities, samples );

    imageT output;
    validityT outputValidity;
    model.response( output, outputValidity, 6, 0 );
    REQUIRE( output( size / 2, size / 2 ) == Approx( 1 ) );
    model.response( output, outputValidity, 6.01, 0 );
    REQUIRE( output( size / 2, size / 2 ) == Approx( 2 ) );
    model.response( output, outputValidity, 20, 0 );
    REQUIRE( output( size / 2, size / 2 ) == Approx( 3 ) );

    REQUIRE_THROWS_AS( model.response( output, outputValidity, -1, 0 ), std::invalid_argument );
    REQUIRE_THROWS_AS( modelT( { 4, 4 }, size, size ), std::invalid_argument );
}

} // namespace RadialPSFModel_test
} // namespace unitTest
