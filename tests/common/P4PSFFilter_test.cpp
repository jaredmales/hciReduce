/** \file P4PSFFilter_test.cpp
 * \brief Tests normalized local filtering with spatially variable P4 PSFs.
 */

#include "../catch2/catch.hpp"

#include "src/common/P4PSFFilter.hpp"

#include <cmath>
#include <limits>

namespace unitTest
{
namespace P4PSFFilter_test
{

/// Production local filter used by the identity regression tests.
using filterT = mx::improc::P4PSFFilter;

/// Float science and response stamp type.
using imageT = filterT::imageT;

/// Per-pixel response validity type.
using validityT = filterT::validityT;

/// Verify the production filter returns the normalized amplitude of a signed asymmetric response.
/** This exercises mx::improc::P4PSFFilter::calculate() with complete support.
 * \ingroup P4PSFFilter_unit_tests
 */
TEST_CASE( "P4 local PSF filter estimates normalized amplitude", "[P4PSFFilter][amplitude][signed]" )
{
    imageT response( 3, 3 );
    response << 0.1F, -0.2F, 0.05F, -0.3F, 1.0F, 0.2F, 0.08F, -0.1F, 0.04F;
    const validityT validity = validityT::Ones( 3, 3 );
    imageT science = imageT::Zero( 7, 7 );
    science.block( 2, 3, 3, 3 ) = 2.5F * response;

    const mx::improc::P4PSFFilterResult result = filterT::calculate( science, response, validity, 3, 4, 1.0 );
    REQUIRE( result.valid );
    REQUIRE( result.amplitude == Approx( 2.5 ).margin( 1e-6 ) );
    REQUIRE( result.correlation == Approx( 2.5 * result.normalization ).margin( 1e-6 ) );
    REQUIRE( result.normalization > 0 );
    REQUIRE( result.supportFraction == Approx( 1.0 ) );
}

/// Verify edge, PSF-validity, and finite-science samples control local filter support.
/** This exercises mx::improc::P4PSFFilter::calculate() with partial support and its configured threshold.
 * \ingroup P4PSFFilter_unit_tests
 */
TEST_CASE( "P4 local PSF filter reports partial support", "[P4PSFFilter][support][edge][validity]" )
{
    const imageT response = imageT::Ones( 3, 3 );
    validityT validity = validityT::Ones( 3, 3 );
    imageT science = imageT::Constant( 4, 4, 3.0F );

    mx::improc::P4PSFFilterResult result = filterT::calculate( science, response, validity, 0, 0, 0.4 );
    REQUIRE( result.valid );
    REQUIRE( result.amplitude == Approx( 3.0 ) );
    REQUIRE( result.normalization == Approx( 4.0 ) );
    REQUIRE( result.supportFraction == Approx( 4.0 / 9.0 ) );

    result = filterT::calculate( science, response, validity, 0, 0, 0.5 );
    REQUIRE_FALSE( result.valid );

    validity( 1, 1 ) = 0;
    science( 0, 1 ) = std::numeric_limits<float>::quiet_NaN();
    result = filterT::calculate( science, response, validity, 0, 0, 0.2 );
    REQUIRE( result.valid );
    REQUIRE( result.supportFraction == Approx( 2.0 / 9.0 ) );
}

/// Verify degenerate and malformed local filtering inputs are rejected or explicitly invalid.
/** This exercises mx::improc::P4PSFFilter::calculate() validation and zero-normalization behavior.
 * \ingroup P4PSFFilter_unit_tests
 */
TEST_CASE( "P4 local PSF filter validates its contract", "[P4PSFFilter][validation]" )
{
    imageT science = imageT::Ones( 5, 5 );
    imageT response = imageT::Zero( 3, 3 );
    validityT validity = validityT::Ones( 3, 3 );
    REQUIRE_FALSE( filterT::calculate( science, response, validity, 2, 2, 1.0 ).valid );

    imageT empty;
    REQUIRE_THROWS_AS( filterT::calculate( empty, response, validity, 0, 0, 1.0 ), std::invalid_argument );
    REQUIRE_THROWS_AS( filterT::calculate( science, imageT::Zero( 2, 3 ), validityT::Ones( 2, 3 ), 2, 2, 1.0 ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( filterT::calculate( science, response, validityT::Ones( 2, 3 ), 2, 2, 1.0 ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( filterT::calculate( science, response, validity, -1, 2, 1.0 ), std::out_of_range );
    REQUIRE_THROWS_AS( filterT::calculate( science, response, validity, 2, 2, -0.1 ), std::invalid_argument );
    REQUIRE_THROWS_AS(
        filterT::calculate( science, response, validity, 2, 2, std::numeric_limits<double>::quiet_NaN() ),
        std::invalid_argument );

    response( 1, 1 ) = std::numeric_limits<float>::infinity();
    REQUIRE_THROWS_AS( filterT::calculate( science, response, validity, 2, 2, 1.0 ), std::invalid_argument );
}

/// Verify shared PSF filter products retain the exact and sequential P4 naming contract.
/** This exercises mx::improc::psfFilterProductPath() and mx::improc::psfFilterDiagnosticPath().
 * \ingroup P4PSFFilter_unit_tests
 */
TEST_CASE( "PSF filter product paths are deterministic", "[P4PSFFilter][products][naming]" )
{
    REQUIRE( mx::improc::psfFilterProductPath( "/tmp/finim.fits", "filtered", false ) == "/tmp/finim_filtered.fits" );
    REQUIRE( mx::improc::psfFilterDiagnosticPath( "/tmp/finim.fits", "filter_support", false ) ==
             "/tmp/finim_outputs/finim_filter_support.fits" );
    REQUIRE( mx::improc::psfFilterProductPath( "/tmp/finim0000.fits", "filtered", true ) ==
             "/tmp/finim_filtered_0000.fits" );
    REQUIRE( mx::improc::psfFilterDiagnosticPath( "/tmp/finim_0000.fits", "filter_validity", true ) ==
             "/tmp/finim_0000_outputs/finim_filter_validity_0000.fits" );
    REQUIRE_THROWS_AS( mx::improc::psfFilterProductPath( "/tmp/finim.fits", "", false ), std::invalid_argument );
    REQUIRE_THROWS_AS( mx::improc::psfFilterProductPath( "/tmp/finim.fits", "filtered", true ), std::invalid_argument );
}

/** \brief Verify explicit identity noise preserves P4PSFFilter::calculate results and the original accumulation.
 * \ingroup P4PSFFilter_unit_tests
 */
TEST_CASE( "Explicit identity preserves the original PSF filter", "[P4PSFFilter][noise][identity]" )
{
    imageT response( 3, 5 );
    for( Eigen::Index pixel = 0; pixel < response.size(); ++pixel )
    {
        response.data()[pixel] = static_cast<float>( pixel - 4 ) / 7;
    }
    imageT science = imageT::Constant( 7, 7, 0.37F );
    science( 1, 1 ) = std::numeric_limits<float>::quiet_NaN();
    validityT validity = validityT::Ones( 3, 5 );
    validity( 0, 2 ) = 0;
    const auto noise = mx::improc::PSFNoiseModel::identity( response.size() );
    for( const int center : { 0, 1, 3, 6 } )
    {
        for( const double threshold : { 0.0, 0.5, 1.0 } )
        {
            const auto legacy =
                mx::improc::P4PSFFilter::calculate( science, response, validity, center, center, threshold );
            const auto explicitIdentity = mx::improc::P4PSFFilter::calculateWeighted( science,
                                                                                      response,
                                                                                      validity,
                                                                                      center,
                                                                                      center,
                                                                                      threshold,
                                                                                      noise );
            REQUIRE( explicitIdentity.amplitude == legacy.amplitude );
            REQUIRE( explicitIdentity.correlation == legacy.correlation );
            REQUIRE( explicitIdentity.normalization == legacy.normalization );
            REQUIRE( explicitIdentity.supportFraction == legacy.supportFraction );
            REQUIRE( explicitIdentity.valid == legacy.valid );
            REQUIRE( explicitIdentity.conditionalSigma == legacy.conditionalSigma );
            REQUIRE( explicitIdentity.score == legacy.score );
            REQUIRE( explicitIdentity.nonnegativeLogLikelihood == legacy.nonnegativeLogLikelihood );
            long double correlation = 0;
            long double normalization = 0;
            int count = 0;
            for( int column = 0; column < response.cols(); ++column )
            {
                for( int row = 0; row < response.rows(); ++row )
                {
                    const int ir = center + row - 1;
                    const int ic = center + column - 2;
                    if( validity( row, column ) && ir >= 0 && ir < science.rows() && ic >= 0 && ic < science.cols() &&
                        std::isfinite( science( ir, ic ) ) )
                    {
                        correlation += static_cast<long double>( response( row, column ) ) * science( ir, ic );
                        normalization += static_cast<long double>( response( row, column ) ) * response( row, column );
                        ++count;
                    }
                }
            }
            REQUIRE( legacy.correlation == static_cast<double>( correlation ) );
            REQUIRE( legacy.normalization == static_cast<double>( normalization ) );
            REQUIRE( legacy.supportFraction == static_cast<double>( count ) / response.size() );
            if( legacy.valid )
            {
                REQUIRE( legacy.amplitude == static_cast<double>( correlation / normalization ) );
            }
        }
    }
}

/** \brief Verify P4PSFFilter::calculateWeighted against a dense likelihood on each retained support.
 * \ingroup P4PSFFilter_unit_tests
 */
TEST_CASE( "Weighted PSF filtering matches dense likelihoods with changing support", "[P4PSFFilter][noise][dense]" )
{
    imageT response( 3, 5 );
    imageT science( 5, 7 );
    Eigen::VectorXd mean( response.size() );
    Eigen::MatrixXd factor( response.size(), 2 );
    const Eigen::VectorXd variances = Eigen::VectorXd::LinSpaced( response.size(), 0.2, 1.4 );
    for( Eigen::Index pixel = 0; pixel < response.size(); ++pixel )
    {
        response.data()[pixel] = static_cast<float>( pixel % 5 - 2 ) * 0.25F;
        mean( pixel ) = 0.1 * pixel;
        factor( pixel, 0 ) = 0.1 * pixel;
        factor( pixel, 1 ) = std::sin( static_cast<double>( pixel ) );
    }
    for( Eigen::Index pixel = 0; pixel < science.size(); ++pixel )
    {
        science.data()[pixel] = static_cast<float>( pixel % 11 - 4 ) * 0.4F;
    }
    for( const bool correlated : { false, true } )
    {
        const auto noise = correlated ? mx::improc::PSFNoiseModel::lowRank( variances, factor, mean )
                                      : mx::improc::PSFNoiseModel::diagonal( variances, mean );
        Eigen::MatrixXd covariance = variances.asDiagonal();
        if( correlated )
        {
            covariance += factor * factor.transpose();
        }
        for( const int center : { 2, 0, 4, 2 } )
        {
            validityT validity = validityT::Ones( 3, 5 );
            if( center != 2 )
            {
                validity( 1, 1 ) = 0;
                science( 1, 1 ) = std::numeric_limits<float>::quiet_NaN();
            }
            std::vector<Eigen::Index> support;
            std::vector<double> data;
            for( int column = 0; column < response.cols(); ++column )
            {
                for( int row = 0; row < response.rows(); ++row )
                {
                    const int ir = center + row - 1;
                    const int ic = center + column - 2;
                    if( validity( row, column ) && ir >= 0 && ir < science.rows() && ic >= 0 && ic < science.cols() &&
                        std::isfinite( science( ir, ic ) ) )
                    {
                        support.push_back( row + column * response.rows() );
                        data.push_back( science( ir, ic ) );
                    }
                }
            }
            const Eigen::Index count = static_cast<Eigen::Index>( support.size() );
            Eigen::MatrixXd marginal( count, count );
            Eigen::VectorXd signal( count ), centered( count );
            for( Eigen::Index row = 0; row < count; ++row )
            {
                signal( row ) = response.data()[support[row]];
                centered( row ) = data[row] - mean( support[row] );
                for( Eigen::Index column = 0; column < count; ++column )
                {
                    marginal( row, column ) = covariance( support[row], support[column] );
                }
            }
            const Eigen::LLT<Eigen::MatrixXd> dense( marginal );
            REQUIRE( dense.info() == Eigen::Success );
            const double numerator = signal.dot( dense.solve( centered ) );
            const double denominator = signal.dot( dense.solve( signal ) );
            const auto actual =
                mx::improc::P4PSFFilter::calculateWeighted( science, response, validity, center, center, 0, noise );
            REQUIRE( actual.valid );
            REQUIRE( actual.correlation == Approx( numerator ).margin( 1e-11 ) );
            REQUIRE( actual.normalization == Approx( denominator ).epsilon( 1e-12 ) );
            REQUIRE( actual.amplitude == Approx( numerator / denominator ).margin( 1e-11 ) );
            REQUIRE( actual.conditionalSigma == Approx( 1 / std::sqrt( denominator ) ).epsilon( 1e-12 ) );
            REQUIRE( actual.score == Approx( numerator / std::sqrt( denominator ) ).margin( 1e-11 ) );
            REQUIRE( actual.nonnegativeLogLikelihood ==
                     Approx( 0.5 * std::pow( std::max( 0.0, actual.score ), 2 ) ).margin( 1e-11 ) );
            REQUIRE( actual.supportFraction == static_cast<double>( count ) / response.size() );
            if( count < response.size() )
            {
                REQUIRE_FALSE(
                    mx::improc::P4PSFFilter::calculateWeighted( science, response, validity, center, center, 1, noise )
                        .valid );
            }
        }
    }
}

/** \brief Verify P4PSFFilter::calculateWeighted amplitude, conditional uncertainty, and score under template/covariance
 * scaling.
 * \ingroup P4PSFFilter_unit_tests
 */
TEST_CASE( "Weighted PSF normalization respects template and noise units", "[P4PSFFilter][noise][normalization]" )
{
    imageT response( 3, 3 );
    response << -0.25F, 0, 0.125F, -0.5F, 1, -0.25F, 0.125F, 0.25F, -0.125F;
    const imageT science = 2.5F * response + 3.0F;
    const validityT validity = validityT::Ones( 3, 3 );
    const Eigen::VectorXd variance = Eigen::VectorXd::LinSpaced( 9, 0.25, 2 );
    const Eigen::VectorXd mean = Eigen::VectorXd::Constant( 9, 3 );
    const Eigen::MatrixXd factor = Eigen::VectorXd::LinSpaced( 9, -1, 1 );
    for( const bool correlated : { false, true } )
    {
        const auto noise = correlated ? mx::improc::PSFNoiseModel::lowRank( variance, factor, mean )
                                      : mx::improc::PSFNoiseModel::diagonal( variance, mean );
        const auto reference =
            mx::improc::P4PSFFilter::calculateWeighted( science, response, validity, 1, 1, 1, noise );
        REQUIRE( reference.valid );
        REQUIRE( reference.amplitude == Approx( 2.5 ).epsilon( 1e-13 ) );
        for( const float scale : { 0.25F, 4.0F, -2.0F } )
        {
            const imageT scaledResponse = scale * response;
            const auto scaled =
                mx::improc::P4PSFFilter::calculateWeighted( science, scaledResponse, validity, 1, 1, 1, noise );
            REQUIRE( scaled.valid );
            REQUIRE( scaled.amplitude == Approx( reference.amplitude / scale ).epsilon( 1e-13 ) );
            REQUIRE( scaled.conditionalSigma ==
                     Approx( reference.conditionalSigma / std::abs( scale ) ).epsilon( 1e-13 ) );
            REQUIRE( scaled.score == Approx( ( scale > 0 ? 1 : -1 ) * reference.score ).epsilon( 1e-13 ) );
            if( scale < 0 )
            {
                REQUIRE( scaled.nonnegativeLogLikelihood == 0 );
            }
        }
        const auto scaledNoise = correlated ? mx::improc::PSFNoiseModel::lowRank( 9 * variance, 3 * factor, mean )
                                            : mx::improc::PSFNoiseModel::diagonal( 9 * variance, mean );
        const auto scaled =
            mx::improc::P4PSFFilter::calculateWeighted( science, response, validity, 1, 1, 1, scaledNoise );
        REQUIRE( scaled.valid );
        REQUIRE( scaled.amplitude == Approx( reference.amplitude ).epsilon( 1e-13 ) );
        REQUIRE( scaled.conditionalSigma == Approx( 3 * reference.conditionalSigma ).epsilon( 1e-13 ) );
        REQUIRE( scaled.score == Approx( reference.score / 3 ).epsilon( 1e-13 ) );
    }
}

/** \brief Verify P4PSFFilter::calculateWeighted enforces support and finite-input contracts.
 * \ingroup P4PSFFilter_unit_tests
 */
TEST_CASE( "Weighted PSF filtering rejects malformed or unusable stamps", "[P4PSFFilter][noise][validation]" )
{
    imageT science = imageT::Ones( 3, 3 );
    imageT response = imageT::Ones( 3, 3 );
    validityT validity = validityT::Ones( 3, 3 );
    const auto noise = mx::improc::PSFNoiseModel::diagonal( Eigen::VectorXd::Ones( 9 ) );
    const auto hugeVariance = mx::improc::PSFNoiseModel::diagonal( Eigen::VectorXd::Constant( 9, 1e270 ) );
    const auto underflow = mx::improc::P4PSFFilter::calculateWeighted( science,
                                                                       imageT::Constant( 3, 3, 1e-30F ),
                                                                       validity,
                                                                       1,
                                                                       1,
                                                                       0,
                                                                       hugeVariance );
    REQUIRE_FALSE( underflow.valid );
    REQUIRE( underflow.normalization == 0 );
    REQUIRE_THROWS_AS( mx::improc::P4PSFFilter::calculateWeighted( science,
                                                                   response,
                                                                   validity,
                                                                   1,
                                                                   1,
                                                                   1,
                                                                   mx::improc::PSFNoiseModel::identity( 8 ) ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( mx::improc::P4PSFFilter::calculateWeighted( science, response, validity, -1, 1, 1, noise ),
                       std::out_of_range );
    REQUIRE_THROWS_AS( mx::improc::P4PSFFilter::calculateWeighted( science, response, validity, 1, 1, 1.1, noise ),
                       std::invalid_argument );
    REQUIRE_FALSE(
        mx::improc::P4PSFFilter::calculateWeighted( science, response, validityT::Zero( 3, 3 ), 1, 1, 0, noise )
            .valid );
    REQUIRE_FALSE(
        mx::improc::P4PSFFilter::calculateWeighted( science, imageT::Zero( 3, 3 ), validity, 1, 1, 0, noise ).valid );
    response( 0, 0 ) = std::numeric_limits<float>::infinity();
    // Response validation still occurs before the image-edge exclusion.
    REQUIRE_THROWS_AS( mx::improc::P4PSFFilter::calculateWeighted( science, response, validity, 0, 0, 0, noise ),
                       std::invalid_argument );
    validity( 0, 0 ) = 0;
    REQUIRE( mx::improc::P4PSFFilter::calculateWeighted( science, response, validity, 0, 0, 0, noise ).valid );
    science.setConstant( std::numeric_limits<float>::quiet_NaN() );
    REQUIRE_FALSE( mx::improc::P4PSFFilter::calculateWeighted( science, response, validity, 1, 1, 0, noise ).valid );
}

} // namespace P4PSFFilter_test
} // namespace unitTest
