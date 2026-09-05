/** \file KLIPPSFModel_test.cpp
 * \brief Tests sparse frozen-basis KLIP PSF response calculations.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "src/common/KLIPPSFModel.hpp"

#include <numbers>
#include <vector>

#include <mx/improc/imageTransforms.hpp>

namespace unitTest
{
namespace KLIPPSFModel_test
{

/** \cond KLIPPSFModel_test_harness */
using modelT = mx::improc::KLIPPSFModel;
using imageT = modelT::imageT;
using vectorT = modelT::vectorT;

/// Return every full-image index in column-major image order.
std::vector<std::size_t> fullRegion( int rows /**< [in] image rows */, int columns /**< [in] image columns */ )
{
    std::vector<std::size_t> indices;
    indices.reserve( static_cast<std::size_t>( rows ) * columns );
    for( int column = 0; column < columns; ++column )
    {
        for( int row = 0; row < rows; ++row )
        {
            indices.push_back( static_cast<std::size_t>( column ) * rows + row );
        }
    }
    return indices;
}

/** \endcond */

/// Verify a sampled KLIP probe follows detector geometry, image-mean centering, and the regional mask.
/** This exercises mx::improc::KLIPPSFModel::probe() with the same linear preprocessing used by the production KLIP
 * worker.
 * \ingroup KLIPPSFModel_unit_tests
 */
TEST_CASE( "KLIP PSF probes use production centering", "[KLIPPSFModel][probe]" )
{
    imageT psf = imageT::Zero( 5, 5 );
    psf( 2, 2 ) = 1;
    modelT model( psf, 9, 9, 4, 4, 5 );
    const std::vector<std::size_t> indices{ 4 * 9 + 4, 5 * 9 + 4, 4 * 9 + 5, 5 * 9 + 5 };

    vectorT uncentered;
    model.probe( uncentered, indices, 4, 4, 0, mx::improc::HCI::meanSub::none );
    REQUIRE( uncentered( 0 ) == Approx( 1 ) );
    REQUIRE( uncentered.tail( 3 ).isZero( 1e-6 ) );

    imageT mask( 4, 1 );
    mask << 1, 1, 0, 0;
    vectorT centered;
    model.probe( centered, indices, 4, 4, 0, mx::improc::HCI::meanSub::imageMean, mask );
    REQUIRE( centered( 0 ) == Approx( 0.5 ) );
    REQUIRE( centered( 1 ) == Approx( -0.5 ) );
    REQUIRE( centered( 2 ) == Approx( 0 ) );
    REQUIRE( centered( 3 ) == Approx( 0 ) );

    REQUIRE_THROWS_AS( model.probe( centered, indices, 4, 4, 0, mx::improc::HCI::meanSub::imageMedian ),
                       std::invalid_argument );
}

/// Verify frozen KL-mode projection returns direct-minus-projected responses for each configured mode count.
/** This exercises mx::improc::KLIPPSFModel::residuals() and
 * mx::improc::KLIPPSFModel::projectedResidual() against explicit orthogonal projections.
 * \ingroup KLIPPSFModel_unit_tests
 */
TEST_CASE( "KLIP PSF responses use the frozen basis", "[KLIPPSFModel][projection]" )
{
    vectorT probe( 3 );
    probe << 2, 3, 4;
    imageT modes = imageT::Zero( 2, 3 );
    modes( 0, 1 ) = 1;
    modes( 1, 0 ) = 1;

    std::vector<vectorT> residuals;
    modelT::residuals( residuals, probe, modes, { 1, 2, 8 } );
    REQUIRE( residuals.size() == 3 );
    REQUIRE( residuals[0].isApprox( ( vectorT( 3 ) << 0, 3, 4 ).finished() ) );
    REQUIRE( residuals[1].isApprox( ( vectorT( 3 ) << 0, 0, 4 ).finished() ) );
    REQUIRE( residuals[2].isApprox( residuals[1] ) );

    imageT projection = modes.matrix().transpose() * modes.matrix();
    vectorT projected;
    modelT::projectedResidual( projected, probe, projection );
    REQUIRE( projected.isApprox( residuals[1] ) );
}

/// Verify regional residuals accumulate into the same compact stamp as a full detector image rotation.
/** This exercises mx::improc::KLIPPSFModel::regionIndex() and mx::improc::KLIPPSFModel::accumulate() for complete
 * and split non-overlapping regions.
 * \ingroup KLIPPSFModel_unit_tests
 */
TEST_CASE( "KLIP PSF regional responses accumulate in the sky frame", "[KLIPPSFModel][derotation]" )
{
    imageT psf = imageT::Zero( 5, 5 );
    psf( 2, 2 ) = 1;
    modelT model( psf, 11, 11, 5, 5, 5 );

    std::vector<std::size_t> firstIndices;
    std::vector<std::size_t> secondIndices;
    for( int column = 0; column < 11; ++column )
    {
        for( int row = 0; row < 11; ++row )
        {
            ( column < 5 ? firstIndices : secondIndices ).push_back( static_cast<std::size_t>( column * 11 + row ) );
        }
    }
    const auto firstLookup = modelT::regionIndex( firstIndices, 11, 11 );
    const auto secondLookup = modelT::regionIndex( secondIndices, 11, 11 );
    vectorT firstResponse( static_cast<Eigen::Index>( firstIndices.size() ) );
    vectorT secondResponse( static_cast<Eigen::Index>( secondIndices.size() ) );
    for( Eigen::Index index = 0; index < firstResponse.rows(); ++index )
    {
        firstResponse( index ) = static_cast<float>( firstIndices[static_cast<std::size_t>( index )] );
    }
    for( Eigen::Index index = 0; index < secondResponse.rows(); ++index )
    {
        secondResponse( index ) = static_cast<float>( secondIndices[static_cast<std::size_t>( index )] );
    }

    imageT stamp;
    modelT::validityT validity;
    model.accumulate( stamp, validity, firstResponse, firstLookup, 5, 5, 0 );
    model.accumulate( stamp, validity, secondResponse, secondLookup, 5, 5, 0 );
    REQUIRE( validity.isConstant( 1 ) );
    for( int column = 0; column < 5; ++column )
    {
        for( int row = 0; row < 5; ++row )
        {
            REQUIRE( stamp( row, column ) == Approx( static_cast<float>( ( column + 3 ) * 11 + row + 3 ) ) );
        }
    }

    imageT rotatedStamp;
    modelT::validityT rotatedValidity;
    const auto allIndices = fullRegion( 11, 11 );
    vectorT allResponse( static_cast<Eigen::Index>( allIndices.size() ) );
    for( Eigen::Index index = 0; index < allResponse.rows(); ++index )
    {
        allResponse( index ) = static_cast<float>( allIndices[static_cast<std::size_t>( index )] );
    }
    model.accumulate( rotatedStamp,
                      rotatedValidity,
                      allResponse,
                      modelT::regionIndex( allIndices, 11, 11 ),
                      5,
                      6,
                      0.5 * std::numbers::pi );
    REQUIRE( rotatedValidity.sum() == 20 );
    imageT detectorResponse( 11, 11 );
    for( int column = 0; column < detectorResponse.cols(); ++column )
    {
        for( int row = 0; row < detectorResponse.rows(); ++row )
        {
            detectorResponse( row, column ) = static_cast<float>( column * 11 + row );
        }
    }
    imageT directSkyResponse;
    mx::improc::imageRotate( directSkyResponse,
                             detectorResponse,
                             0.5 * std::numbers::pi,
                             mx::improc::cubicConvolTransform<float>() );
    for( int column = 0; column < rotatedStamp.cols(); ++column )
    {
        for( int row = 0; row < rotatedStamp.rows(); ++row )
        {
            REQUIRE( rotatedStamp( row, column ) == Approx( directSkyResponse( row + 3, column + 4 ) ).margin( 1e-5 ) );
        }
    }
}

} // namespace KLIPPSFModel_test
} // namespace unitTest
