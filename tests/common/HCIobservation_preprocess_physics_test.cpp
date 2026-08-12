/** \file HCIobservation_preprocess_physics_test.cpp
 * \brief Tests HCIobservation image-domain preprocessing behavior.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "HCIobservation_test_fixture.hpp"

#include <cmath>
#include <limits>

namespace unitTest
{
namespace HCIobservation_preprocess_physics_test
{

/// \cond HCIobservation_test_harness
namespace
{
bool allFinite( const HCIobservationTestHarness::imageT &image )
{
    for( int col = 0; col < image.cols(); ++col )
    {
        for( int row = 0; row < image.rows(); ++row )
        {
            if( !std::isfinite( image( row, col ) ) )
            {
                return false;
            }
        }
    }
    return true;
}

HCIobservationTestHarness::imageT azimuthalOracle( const HCIobservationTestHarness::imageT &image,
                                                   float radialHalfWidth,
                                                   float azimuthalHalfWidth,
                                                   float maximumAzimuth )
{
    using imageT = HCIobservationTestHarness::imageT;
    using kernelT = mx::improc::azBoxKernel<imageT>;

    kernelT kernel( radialHalfWidth, azimuthalHalfWidth, maximumAzimuth );
    mx::improc::precalcKernel<kernelT> cache( kernel,
                                              image.rows(),
                                              image.cols(),
                                              0.5F * ( image.rows() - 1 ),
                                              0.5F * ( image.cols() - 1 ) );
    imageT filtered;
    mx::improc::medianFilterImage( filtered, image, cache );
    return image - filtered;
}

HCIobservationTestHarness::imageT syntheticPolarImage( int rows, int cols )
{
    HCIobservationTestHarness::imageT image( rows, cols );
    const float centerRow = 0.5F * ( rows - 1 );
    const float centerCol = 0.5F * ( cols - 1 );
    for( int col = 0; col < cols; ++col )
    {
        for( int row = 0; row < rows; ++row )
        {
            const float y = col - centerCol;
            const float x = row - centerRow;
            image( row, col ) = 3 * std::hypot( x, y ) + std::atan2( y, x );
        }
    }
    return image;
}

HCIobservationTestHarness::cubeT compositeOracle( const HCIobservationTestHarness::cubeT &input,
                                                  const HCIobservationTestHarness::imageT &mask )
{
    using imageT = HCIobservationTestHarness::imageT;
    HCIobservationTestHarness::cubeT result = input;

    imageT radius( result.rows(), result.cols() );
    mx::improc::radiusImage( radius );
    for( int plane = 0; plane < result.planes(); ++plane )
    {
        result.image( plane ) *= mask;
        imageT profile;
        auto image = result.image( plane );
        mx::improc::radprofim( profile, image, radius, &mask, true );
        mx::improc::zeroNaNs( image );
        result.image( plane ) = image * mask;
    }

    imageT medianEdgeMask( result.rows(), result.cols() );
    medianEdgeMask.setOnes();
    medianEdgeMask.row( 0 ).setZero();
    medianEdgeMask.row( medianEdgeMask.rows() - 1 ).setZero();
    medianEdgeMask.col( 0 ).setZero();
    medianEdgeMask.col( medianEdgeMask.cols() - 1 ).setZero();
    for( int plane = 0; plane < result.planes(); ++plane )
    {
        imageT filtered( result.rows(), result.cols() );
        filtered.setZero();
        const imageT image = result.image( plane );
        REQUIRE( mx::improc::medianSmooth( filtered, image, 3 ) == 0 );
        result.image( plane ) = ( image - filtered ) * medianEdgeMask * mask;
    }

    const mx::improc::gaussKernel<imageT, 2> gaussian( 1 );
    for( int plane = 0; plane < result.planes(); ++plane )
    {
        imageT filtered;
        const imageT image = result.image( plane );
        REQUIRE( mx::improc::filterImage( filtered, image, gaussian, 0 ) == mx::error_t::noerror );
        result.image( plane ) = ( image - filtered ) * mask;
    }

    for( int plane = 0; plane < result.planes(); ++plane )
    {
        const imageT image = result.image( plane );
        result.image( plane ) = azimuthalOracle( image, 1, 2, 45 ) * mask;
    }

    for( int col = 0; col < result.cols(); ++col )
    {
        for( int row = 0; row < result.rows(); ++row )
        {
            if( mask( row, col ) == 0 )
            {
                continue;
            }

            float mean = 0;
            for( int plane = 0; plane < result.planes(); ++plane )
            {
                mean += result.image( plane )( row, col );
            }
            mean /= result.planes();

            float sumSquares = 0;
            for( int plane = 0; plane < result.planes(); ++plane )
            {
                result.image( plane )( row, col ) -= mean;
                sumSquares += result.image( plane )( row, col ) * result.image( plane )( row, col );
            }

            const float rms = std::sqrt( sumSquares / result.planes() );
            if( rms > 0 )
            {
                for( int plane = 0; plane < result.planes(); ++plane )
                {
                    result.image( plane )( row, col ) /= rms;
                }
            }
        }
    }

    return result;
}
} // namespace
/// \endcond

/// Verify HCIobservation::preProcess subtracts radial backgrounds while preserving a localized source.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation radial-profile preprocessing", "[HCIobservation][preprocess][radial-profile]" )
{
    OpenMPThreadGuard threads;
    HCIobservationTestHarness observation;
    observation.m_preProcess_subradprof = true;

    HCIobservationTestHarness::cubeT cube( 11, 11, 1 );
    cube.image( 0 ).setConstant( 5 );
    cube.image( 0 )( 5, 7 ) = 105;
    observation.preProcess( cube );

    REQUIRE( allFinite( cube.image( 0 ) ) );
    REQUIRE( cube.image( 0 )( 5, 5 ) == Approx( 0 ).margin( 1e-5 ) );
    REQUIRE( cube.image( 0 )( 5, 7 ) > 90 );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::preProcess( cube );
#endif
    // clang-format on
}

/// Verify HCIobservation::preProcess excludes heavily masked pixels from radial-profile statistics.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation masked radial-profile preprocessing", "[HCIobservation][preprocess][radial-profile][mask]" )
{
    OpenMPThreadGuard threads;
    HCIobservationTestHarness observation;
    observation.m_preProcess_subradprof = true;
    observation.m_maskFile = "configured-mask.fits";
    observation.m_mask.resize( 11, 11 );
    observation.m_mask.setZero();
    for( int col = 5; col < observation.m_mask.cols(); ++col )
    {
        observation.m_mask( 5, col ) = 1;
    }

    HCIobservationTestHarness::cubeT cube( 11, 11, 1 );
    cube.image( 0 ).setConstant( 5 );
    observation.preProcess( cube );

    REQUIRE( allFinite( cube.image( 0 ) ) );
    REQUIRE( cube.image( 0 ).abs().maxCoeff() == Approx( 0 ).margin( 1e-5 ) );
}

/// Verify HCIobservation::preProcess median USM uses exact odd/even windows and deterministic zero edges.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation median unsharp masking", "[HCIobservation][preprocess][median-USM]" )
{
    OpenMPThreadGuard threads;
    HCIobservationTestHarness observation;

    HCIobservationTestHarness::cubeT cube( 5, 5, 1 );
    cube.image( 0 ).setConstant( 7 );
    observation.m_preProcess_medianUSM_fwhm = 3;
    observation.preProcess( cube );
    REQUIRE( allFinite( cube.image( 0 ) ) );
    REQUIRE( cube.image( 0 ).abs().maxCoeff() == 0 );

    for( int col = 0; col < cube.cols(); ++col )
    {
        for( int row = 0; row < cube.rows(); ++row )
        {
            cube.image( 0 )( row, col ) = row + 5 * col;
        }
    }
    observation.m_preProcess_medianUSM_fwhm = 4;
    observation.preProcess( cube );

    for( int col = 0; col < cube.cols(); ++col )
    {
        for( int row = 0; row < cube.rows(); ++row )
        {
            const bool hasCompleteWindow = row >= 2 && row < 4 && col >= 2 && col < 4;
            if( hasCompleteWindow )
            {
                REQUIRE( cube.image( 0 )( row, col ) == Approx( 3 ) );
            }
            else
            {
                REQUIRE( cube.image( 0 )( row, col ) == 0 );
            }
        }
    }
}

/// Verify HCIobservation::preProcess rejects a median USM window larger than either image dimension.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation median unsharp-mask validation", "[HCIobservation][preprocess][median-USM][validation]" )
{
    OpenMPThreadGuard threads;
    HCIobservationTestHarness observation;
    observation.m_preProcess_medianUSM_fwhm = 4;

    HCIobservationTestHarness::cubeT cube( 3, 5, 1 );
    cube.image( 0 ).setOnes();
    REQUIRE_THROWS( observation.preProcess( cube ) );
    REQUIRE( cube.image( 0 ).minCoeff() == 1 );
    REQUIRE( cube.image( 0 ).maxCoeff() == 1 );
}

/// Verify HCIobservation::preProcess Gaussian USM maps a constant to zero and matches its normalized impulse kernel.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation Gaussian unsharp masking", "[HCIobservation][preprocess][Gaussian-USM]" )
{
    OpenMPThreadGuard threads;
    HCIobservationTestHarness observation;
    observation.m_preProcess_gaussUSM_fwhm = 1;

    HCIobservationTestHarness::cubeT cube( 17, 17, 1 );
    cube.image( 0 ).setConstant( 7 );
    observation.preProcess( cube );
    REQUIRE( allFinite( cube.image( 0 ) ) );
    REQUIRE( cube.image( 0 ).abs().maxCoeff() == Approx( 0 ).margin( 1e-5 ) );

    cube.image( 0 ).setZero();
    cube.image( 0 )( 8, 8 ) = 1;
    observation.preProcess( cube );

    mx::improc::gaussKernel<HCIobservationTestHarness::imageT, 2> kernel( 1 );
    for( int col = 0; col < cube.cols(); ++col )
    {
        for( int row = 0; row < cube.rows(); ++row )
        {
            float expected = row == 8 && col == 8 ? 1 : 0;
            if( std::abs( row - 8 ) <= 1 && std::abs( col - 8 ) <= 1 )
            {
                expected -= kernel.kernel( row - 7, col - 7 );
            }
            REQUIRE( cube.image( 0 )( row, col ) == Approx( expected ).margin( 1e-6 ) );
        }
    }
}

/// Verify HCIobservation::preProcess validates Gaussian USM widths before constructing an unsafe kernel.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation Gaussian unsharp-mask validation", "[HCIobservation][preprocess][Gaussian-USM][validation]" )
{
    OpenMPThreadGuard threads;
    HCIobservationTestHarness observation;
    HCIobservationTestHarness::cubeT cube( 9, 9, 1 );
    cube.image( 0 ).setOnes();

    observation.m_preProcess_gaussUSM_fwhm = 5;
    REQUIRE_THROWS( observation.preProcess( cube ) );

    observation.m_preProcess_gaussUSM_fwhm = std::numeric_limits<float>::max();
    REQUIRE_THROWS( observation.preProcess( cube ) );

    observation.m_preProcess_gaussUSM_fwhm = std::numeric_limits<float>::infinity();
    REQUIRE_THROWS( observation.preProcess( cube ) );

    observation.m_preProcess_gaussUSM_fwhm = -1;
    REQUIRE_THROWS( observation.preProcess( cube ) );
}

/// Verify HCIobservation::preProcess azimuthal USM honors half-widths and maxAz on odd and even images.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation azimuthal unsharp masking", "[HCIobservation][preprocess][azimuthal-USM]" )
{
    OpenMPThreadGuard threads;

    for( const int size : { 9, 10 } )
    {
        HCIobservationTestHarness observation;
        observation.m_preProcess_azUSM_radHalfWidth = 1;
        observation.m_preProcess_azUSM_azHalfWidth = 2;
        observation.m_preProcess_azUSM_maxAz = 30;

        const auto input = syntheticPolarImage( size, size );
        const auto expected = azimuthalOracle( input, 1, 2, 30 );
        HCIobservationTestHarness::cubeT cube( size, size, 1 );
        cube.image( 0 ) = input;
        observation.preProcess( cube );

        REQUIRE( allFinite( cube.image( 0 ) ) );
        REQUIRE( cube.image( 0 ).isApprox( expected, 1e-6 ) );
        REQUIRE( std::isfinite( cube.image( 0 )( size / 2, size / 2 ) ) );
    }

    auto input = syntheticPolarImage( 11, 11 );
    input( 5, 8 ) += 50;

    HCIobservationTestHarness limited;
    limited.m_preProcess_azUSM_radHalfWidth = 1;
    limited.m_preProcess_azUSM_azHalfWidth = 4;
    limited.m_preProcess_azUSM_maxAz = 15;
    HCIobservationTestHarness::cubeT limitedCube( 11, 11, 1 );
    limitedCube.image( 0 ) = input;
    limited.preProcess( limitedCube );

    HCIobservationTestHarness unlimited;
    unlimited.m_preProcess_azUSM_radHalfWidth = 1;
    unlimited.m_preProcess_azUSM_azHalfWidth = 4;
    unlimited.m_preProcess_azUSM_maxAz = 180;
    HCIobservationTestHarness::cubeT unlimitedCube( 11, 11, 1 );
    unlimitedCube.image( 0 ) = input;
    unlimited.preProcess( unlimitedCube );

    REQUIRE( limitedCube.image( 0 ).isApprox( azimuthalOracle( input, 1, 4, 15 ), 1e-6 ) );
    REQUIRE( unlimitedCube.image( 0 ).isApprox( azimuthalOracle( input, 1, 4, 180 ), 1e-6 ) );
    REQUIRE( ( limitedCube.image( 0 ) - unlimitedCube.image( 0 ) ).abs().maxCoeff() > 1e-3 );
    REQUIRE( limitedCube.image( 0 )( 5, 8 ) > 40 );
}

/// Verify HCIobservation::preProcess preserves the documented complete preprocessing stage order.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation composite preprocessing order", "[HCIobservation][preprocess][stage-order]" )
{
    OpenMPThreadGuard threads;
    HCIobservationTestHarness observation;
    observation.m_preProcess_subradprof = true;
    observation.m_preProcess_medianUSM_fwhm = 3;
    observation.m_preProcess_gaussUSM_fwhm = 1;
    observation.m_preProcess_azUSM_radHalfWidth = 1;
    observation.m_preProcess_azUSM_azHalfWidth = 2;
    observation.m_preProcess_azUSM_maxAz = 45;
    observation.m_preProcess_meanSubMethod = mx::improc::HCI::meanSub::meanImage;
    observation.m_preProcess_pixelTSNormMethod = mx::improc::HCI::pixelTSNorm::rms;
    observation.m_maskFile = "configured-mask.fits";
    observation.m_mask.resize( 9, 9 );
    observation.m_mask.setOnes();
    observation.m_mask( 1, 1 ) = 0;

    HCIobservationTestHarness::cubeT cube( 9, 9, 2 );
    for( int col = 0; col < cube.cols(); ++col )
    {
        for( int row = 0; row < cube.rows(); ++row )
        {
            const float radius = std::hypot( row - 4.0F, col - 4.0F );
            cube.image( 0 )( row, col ) = 5 + radius;
            cube.image( 1 )( row, col ) = 9 + 0.5F * radius;
        }
    }
    cube.image( 0 )( 4, 6 ) += 10;
    cube.image( 1 )( 5, 6 ) += 6;

    const auto expected = compositeOracle( cube, observation.m_mask );
    observation.preProcess( cube );

    REQUIRE( allFinite( cube.image( 0 ) ) );
    REQUIRE( allFinite( cube.image( 1 ) ) );
    REQUIRE( cube.image( 0 ).isApprox( expected.image( 0 ), 1e-4 ) );
    REQUIRE( cube.image( 1 ).isApprox( expected.image( 1 ), 1e-4 ) );
    REQUIRE( cube.image( 0 )( 1, 1 ) == 0 );
    REQUIRE( cube.image( 1 )( 1, 1 ) == 0 );
}

/// Verify HCIobservation preprocessing and combination are deterministic across OpenMP team sizes.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation OpenMP preprocessing determinism", "[HCIobservation][preprocess][OpenMP]" )
{
    const auto process = []( int threads )
    {
        OpenMPThreadGuard threadGuard( threads );
        HCIobservationTestHarness observation;
        observation.m_preProcess_subradprof = true;
        observation.m_preProcess_medianUSM_fwhm = 3;
        observation.m_preProcess_gaussUSM_fwhm = 1;
        observation.m_preProcess_azUSM_radHalfWidth = 1;
        observation.m_preProcess_azUSM_azHalfWidth = 2;
        observation.m_preProcess_azUSM_maxAz = 45;
        observation.m_preProcess_meanSubMethod = mx::improc::HCI::meanSub::meanImage;
        observation.m_preProcess_pixelTSNormMethod = mx::improc::HCI::pixelTSNorm::rms;
        observation.m_maskFile = "configured-mask.fits";
        observation.m_mask.resize( 11, 11 );
        observation.m_mask.setOnes();
        observation.m_mask( 2, 3 ) = 0;
        observation.m_mask( 8, 7 ) = 0;

        HCIobservationTestHarness::cubeT cube( 11, 11, 4 );
        for( int plane = 0; plane < cube.planes(); ++plane )
        {
            for( int col = 0; col < cube.cols(); ++col )
            {
                for( int row = 0; row < cube.rows(); ++row )
                {
                    const float radius = std::hypot( row - 5.0F, col - 5.0F );
                    cube.image( plane )( row, col ) = 3 * radius + 0.2F * plane * row - 0.1F * ( plane + 1 ) * col;
                }
            }
            cube.image( plane )( 5, 7 ) += 10 + plane;
        }

        observation.preProcess( cube );
        observation.m_maskFile.clear();
        observation.m_psfsub = { cube };
        observation.m_combineMethod = mx::improc::HCI::combine::mean;
        observation.combineFinim();
        return std::pair{ cube, observation.m_finim };
    };

    const auto oneThread = process( 1 );
    const auto fourThreads = process( 4 );
    for( int plane = 0; plane < oneThread.first.planes(); ++plane )
    {
        REQUIRE( oneThread.first.image( plane ).isApprox( fourThreads.first.image( plane ), 1e-6 ) );
    }
    REQUIRE( oneThread.second.image( 0 ).isApprox( fourThreads.second.image( 0 ), 1e-6 ) );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::preProcess( oneThread.first );
    mx::improc::HCIobservation<float, mx::verbose::vv>::combineFinim();
#endif
    // clang-format on
}

} // namespace HCIobservation_preprocess_physics_test
} // namespace unitTest
