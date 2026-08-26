/** \file P4LocalProcessing_test.cpp
 * \brief Tests sparse geometry and trial-source sampling for pixel-local P4 processing.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "src/common/P4LocalProcessing.hpp"

#include <cmath>
#include <vector>

#include <mx/improc/imageTransforms.hpp>
#include <mx/math/geo.hpp>

namespace unitTest
{
namespace P4LocalProcessing_test
{

/// Sparse geometry type under test.
using geometryT = mx::improc::P4LocalGeometry;

/// Trial-source image type under test.
using imageT = mx::improc::P4TrialSource::imageT;

/// Verify the sparse trial-source sampler preserves geometric-center phase and matches production image shifting.
/** This exercises mx::improc::P4TrialSource::configure() and mx::improc::P4TrialSource::value() against
 * mx::improc::imageShift().
 * \ingroup P4LocalProcessing_unit_tests
 */
TEST_CASE( "P4 local trial source matches full image shift", "[P4LocalProcessing][source][shift][phase]" )
{
    imageT sourceTemplate( 16, 16 );
    for( int column = 0; column < sourceTemplate.cols(); ++column )
    {
        for( int row = 0; row < sourceTemplate.rows(); ++row )
        {
            const float deltaRow = static_cast<float>( row ) - 7.5F;
            const float deltaColumn = static_cast<float>( column ) - 7.5F;
            sourceTemplate( row, column ) =
                std::exp( -0.17F * deltaRow * deltaRow - 0.11F * deltaColumn * deltaColumn ) *
                ( 1.0F + 0.04F * deltaRow - 0.03F * deltaColumn );
        }
    }

    const std::vector<double> angles{ 0.17, -0.31 };
    const std::vector<float> scales{ 1, 1.7F };
    constexpr double separation = 2.3;
    constexpr double positionAngle = 31.0;
    constexpr double contrast = -0.4;
    mx::improc::P4TrialSource source;
    source.configure( sourceTemplate, 16, 16, 5, angles, separation, positionAngle, contrast, scales );
    REQUIRE( source.configured() );
    REQUIRE( source.cropRows() == 6 );
    REQUIRE( source.cropColumns() == 6 );

    imageT prepared = imageT::Zero( 16, 16 );
    prepared.block( 5, 5, 6, 6 ) = sourceTemplate.block( 5, 5, 6, 6 );
    for( std::size_t frame = 0; frame < angles.size(); ++frame )
    {
        const float angle =
            mx::math::dtor( -static_cast<float>( positionAngle ) ) + static_cast<float>( angles[frame] );
        const float rowShift = static_cast<float>( separation ) * std::sin( angle );
        const float columnShift = static_cast<float>( separation ) * std::cos( angle );
        imageT shifted;
        mx::improc::imageShift( shifted, prepared, rowShift, columnShift, mx::improc::cubicConvolTransform<float>() );
        shifted *= static_cast<float>( contrast ) * scales[frame];
        for( int column = 0; column < shifted.cols(); ++column )
        {
            for( int row = 0; row < shifted.rows(); ++row )
            {
                REQUIRE( source.value( frame, row, column ) == Approx( shifted( row, column ) ).margin( 3e-7 ) );
            }
        }
    }

    REQUIRE_NOTHROW( source.addSource( angles, separation, positionAngle, -contrast, scales ) );
    for( std::size_t frame = 0; frame < angles.size(); ++frame )
    {
        for( int column = 0; column < sourceTemplate.cols(); ++column )
        {
            for( int row = 0; row < sourceTemplate.rows(); ++row )
            {
                REQUIRE( source.value( frame, row, column ) == Approx( 0 ).margin( 3e-7 ) );
            }
        }
    }

    mx::improc::P4TrialSource integral;
    integral.configure( sourceTemplate, 16, 16, 5, { 0 }, 2, 0, 0.5, { 2 } );
    imageT integralExpected;
    mx::improc::imageShift( integralExpected, prepared, 0, 2, mx::improc::cubicConvolTransform<float>() );
    for( int column = 0; column < integralExpected.cols(); ++column )
    {
        for( int row = 0; row < integralExpected.rows(); ++row )
        {
            REQUIRE( integral.value( 0, row, column ) == Approx( integralExpected( row, column ) ) );
        }
    }

    REQUIRE_THROWS(
        source.configure( sourceTemplate, 15, 16, 5, angles, separation, positionAngle, contrast, scales ) );
    REQUIRE_THROWS( source.value( angles.size(), 0, 0 ) );

    mx::improc::P4TrialSource zeroSource;
    REQUIRE_NOTHROW( zeroSource.configure( sourceTemplate, 16, 16, 5, angles, separation, positionAngle, 0, scales ) );
    for( std::size_t frame = 0; frame < angles.size(); ++frame )
    {
        for( int column = 0; column < sourceTemplate.cols(); ++column )
        {
            for( int row = 0; row < sourceTemplate.rows(); ++row )
            {
                REQUIRE( zeroSource.value( frame, row, column ) == 0 );
            }
        }
    }

    imageT oddTemplate( 15, 15 );
    for( int column = 0; column < oddTemplate.cols(); ++column )
    {
        for( int row = 0; row < oddTemplate.rows(); ++row )
        {
            oddTemplate( row, column ) = static_cast<float>( row + 2 * column ) / 100.0F;
        }
    }
    constexpr float wrappedPositionAngle = 391.0F;
    constexpr float positiveContrast = 0.3F;
    mx::improc::P4TrialSource oddSource;
    oddSource.configure( oddTemplate, 15, 15, 5, { -0.27 }, 3.1, wrappedPositionAngle, positiveContrast, { 0.8F } );
    REQUIRE( oddSource.cropRows() == 5 );
    REQUIRE( oddSource.cropColumns() == 5 );
    imageT oddPrepared = imageT::Zero( 15, 15 );
    oddPrepared.block( 5, 5, 5, 5 ) = oddTemplate.block( 5, 5, 5, 5 );
    const float oddAngle = mx::math::dtor( -wrappedPositionAngle ) - 0.27F;
    imageT oddExpected;
    mx::improc::imageShift( oddExpected,
                            oddPrepared,
                            3.1F * std::sin( oddAngle ),
                            3.1F * std::cos( oddAngle ),
                            mx::improc::cubicConvolTransform<float>() );
    oddExpected *= positiveContrast * 0.8F;
    for( int column = 0; column < oddExpected.cols(); ++column )
    {
        for( int row = 0; row < oddExpected.rows(); ++row )
        {
            REQUIRE( oddSource.value( 0, row, column ) == Approx( oddExpected( row, column ) ).margin( 3e-7 ) );
        }
    }
}

/// Verify sparse local derotation exactly reconstructs the matching crop of production full-image rotation.
/** This exercises mx::improc::P4LocalGeometry::configure() and
 * mx::improc::P4LocalGeometry::outputSample() against mx::improc::imageRotate().
 * \ingroup P4LocalProcessing_unit_tests
 */
TEST_CASE( "P4 local sparse derotation matches full image rotation",
           "[P4LocalProcessing][geometry][derotation][sparse]" )
{
    constexpr int rows = 31;
    constexpr int columns = 31;
    constexpr int stampSize = 5;
    const std::vector<double> angles{ 0, 0.23, -0.31 };
    geometryT::lookupImageT ownership( rows, columns );
    geometryT::lookupImageT searchIndices( rows, columns );
    ownership.setZero();
    for( int column = 0; column < columns; ++column )
    {
        for( int row = 0; row < rows; ++row )
        {
            searchIndices( row, column ) = row + rows * column;
        }
    }

    geometryT geometry;
    geometry.configure( rows, columns, stampSize, 15.2, 14.7, angles, true, ownership, searchIndices );
    REQUIRE( geometry.configured() );
    REQUIRE( geometry.stampSize() == stampSize );
    REQUIRE( geometry.originRow() == 13 );
    REQUIRE( geometry.originColumn() == 13 );
    REQUIRE( geometry.frameCount() == angles.size() );
    REQUIRE_FALSE( geometry.searchRequests().empty() );
    REQUIRE( geometry.sparseSampleCount() >= geometry.searchRequests().size() );
    REQUIRE( geometry.storageBytes() > 0 );

    geometryT halfPixelGeometry;
    halfPixelGeometry.configure( rows, columns, stampSize, 14.5, 15.5, { 0 }, true, ownership, searchIndices );
    REQUIRE( halfPixelGeometry.sourceRow() == 14.5 );
    REQUIRE( halfPixelGeometry.sourceColumn() == 15.5 );
    REQUIRE( halfPixelGeometry.originRow() == 13 );
    REQUIRE( halfPixelGeometry.originColumn() == 14 );

    std::vector<imageT> detectorFrames;
    std::vector<imageT> rotatedFrames;
    detectorFrames.reserve( angles.size() );
    rotatedFrames.reserve( angles.size() );
    for( std::size_t frame = 0; frame < angles.size(); ++frame )
    {
        imageT detector( rows, columns );
        for( int column = 0; column < columns; ++column )
        {
            for( int row = 0; row < rows; ++row )
            {
                detector( row, column ) = static_cast<float>( 0.2 * row - 0.13 * column + 0.01 * row * column + frame );
            }
        }
        detectorFrames.push_back( detector );
        if( angles[frame] == 0 )
        {
            rotatedFrames.push_back( detector );
        }
        else
        {
            imageT rotated;
            mx::improc::imageRotate( rotated,
                                     detector,
                                     static_cast<float>( angles[frame] ),
                                     mx::improc::cubicConvolTransform<float>() );
            rotatedFrames.push_back( std::move( rotated ) );
        }
    }

    for( int stampColumn = 0; stampColumn < stampSize; ++stampColumn )
    {
        for( int stampRow = 0; stampRow < stampSize; ++stampRow )
        {
            for( std::size_t frame = 0; frame < angles.size(); ++frame )
            {
                const mx::improc::P4LocalOutputSample &output = geometry.outputSample( stampRow, stampColumn, frame );
                REQUIRE( output.valid() );
                float sparse{ 0 };
                for( const mx::improc::P4LocalResidualSample &sample : output.samples() )
                {
                    const mx::improc::P4LocalSearchRequest &request = geometry.searchRequests()[sample.requestIndex()];
                    REQUIRE( request.frames()[sample.frameOffset()] == static_cast<int>( frame ) );
                    sparse += detectorFrames[frame]( request.coordinate().row(), request.coordinate().column() ) *
                              sample.weight();
                }
                REQUIRE( sparse == Approx( rotatedFrames[frame]( geometry.originRow() + stampRow,
                                                                 geometry.originColumn() + stampColumn ) )
                                       .margin( 3e-6 ) );
            }
        }
    }

    ownership( 15, 15 ) = -1;
    searchIndices( 15, 15 ) = -1;
    geometryT incomplete;
    incomplete.configure( rows, columns, stampSize, 15, 15, { 0 }, true, ownership, searchIndices );
    REQUIRE_FALSE( incomplete.outputSample( 2, 2, 0 ).valid() );
    REQUIRE_THROWS( incomplete.outputSample( stampSize, 0, 0 ) );
}

} // namespace P4LocalProcessing_test
} // namespace unitTest
