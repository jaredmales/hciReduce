/** \file PSFNoiseTraining_test.cpp
 * \brief Tests annular patch geometry, source exclusion, and production noise-weighted filtering.
 */

#include "../catch2/catch.hpp"

#include "src/common/PSFNoiseTraining.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

namespace unitTest
{
namespace PSFNoiseTraining_test
{

/// Float science and response image type used by the production filter.
using imageT = mx::improc::P4PSFFilter::imageT;

/** \brief Verify PSFNoiseTraining::sample uses half-width spacing and aligns every offset to the candidate lattice.
 * \ingroup PSFNoiseTraining_unit_tests
 */
TEST_CASE( "Annular noise patches have consistent radial geometry", "[PSFNoiseTraining][geometry]" )
{
    imageT science( 81, 81 );
    for( int column = 0; column < science.cols(); ++column )
    {
        for( int row = 0; row < science.rows(); ++row )
        {
            science( row, column ) = 3 * row - 2 * column + 5;
        }
    }
    mx::improc::PSFNoiseTrainingConfig config;
    const auto samples = mx::improc::PSFNoiseTraining::sample( science, 9, 9, 52, 40, config, {} );
    REQUIRE( samples.m_attempted == 19 );
    REQUIRE( samples.m_arcStep == Approx( 24 * std::numbers::pi / 19 ) );
    REQUIRE( samples.m_samples.cols() == 81 );
    REQUIRE( samples.m_samples.rows() > 2 );
    REQUIRE( samples.m_attempted == samples.m_excluded + samples.m_incomplete + samples.m_samples.rows() );
    for( Eigen::Index sample = 0; sample < samples.m_samples.rows(); ++sample )
    {
        const double row = samples.m_centers( sample, 0 );
        const double column = samples.m_centers( sample, 1 );
        REQUIRE( std::hypot( row - 40, column - 40 ) == Approx( 12 ) );
        const double angle = std::atan2( column - 40, row - 40 );
        for( int sc = 0; sc < 9; ++sc )
        {
            for( int sr = 0; sr < 9; ++sr )
            {
                const double rr = row + std::cos( angle ) * ( sr - 4 ) - std::sin( angle ) * ( sc - 4 );
                const double cc = column + std::sin( angle ) * ( sr - 4 ) + std::cos( angle ) * ( sc - 4 );
                REQUIRE( samples.m_samples( sample, sr + sc * 9 ) == Approx( 3 * rr - 2 * cc + 5 ).margin( 1e-12 ) );
            }
        }
    }
    // A different candidate orientation and a rectangular stamp exercise both rotation signs and vectorization.
    const auto rotated = mx::improc::PSFNoiseTraining::sample( science, 3, 5, 40, 52, config, {} );
    for( Eigen::Index sample = 0; sample < rotated.m_samples.rows(); ++sample )
    {
        const double row = rotated.m_centers( sample, 0 );
        const double column = rotated.m_centers( sample, 1 );
        const double angle = std::atan2( column - 40, row - 40 ) - std::numbers::pi / 2;
        const double rr = row + std::cos( angle ) - 2 * std::sin( angle );
        const double cc = column + std::sin( angle ) + 2 * std::cos( angle );
        REQUIRE( rotated.m_samples( sample, 14 ) == Approx( 3 * rr - 2 * cc + 5 ).margin( 1e-12 ) );
    }
}

/** \brief Verify PSFNoiseTraining::sample excludes whole interpolated footprints and never fills missing data.
 * \ingroup PSFNoiseTraining_unit_tests
 */
TEST_CASE( "Annular noise training withholds overlapping source neighborhoods",
           "[PSFNoiseTraining][exclusion][support]" )
{
    imageT science = imageT::Ones( 81, 81 );
    mx::improc::PSFNoiseTrainingConfig config;
    config.m_guardRadius = 2;
    const std::vector<mx::improc::PSFNoiseExclusion> exclusions{ { 40, 20, 4 }, { 20, 40, 6 } };
    const auto before = mx::improc::PSFNoiseTraining::sample( science, 5, 5, 60, 40, config, exclusions );
    for( int column = 0; column < science.cols(); ++column )
    {
        for( int row = 0; row < science.rows(); ++row )
        {
            if( std::hypot( row - 60, column - 40 ) <= std::sqrt( 8.0 ) + 2 ||
                std::hypot( row - 40, column - 20 ) <= 4 || std::hypot( row - 20, column - 40 ) <= 6 )
            {
                science( row, column ) = std::numeric_limits<float>::quiet_NaN();
            }
        }
    }
    const auto after = mx::improc::PSFNoiseTraining::sample( science, 5, 5, 60, 40, config, exclusions );
    REQUIRE( after.m_samples == before.m_samples );
    REQUIRE( after.m_centers == before.m_centers );
    REQUIRE( after.m_excluded == before.m_excluded );
    REQUIRE( after.m_incomplete == 0 );
    for( Eigen::Index sample = 0; sample < after.m_samples.rows(); ++sample )
    {
        for( const auto &exclusion : exclusions )
        {
            REQUIRE( std::hypot( after.m_centers( sample, 0 ) - exclusion.m_row,
                                 after.m_centers( sample, 1 ) - exclusion.m_column ) >
                     after.m_footprintRadius + exclusion.m_radius );
        }
    }
    REQUIRE( after.m_samples.rows() > 0 );
    science( static_cast<int>( std::round( after.m_centers( 0, 0 ) ) ),
             static_cast<int>( std::round( after.m_centers( 0, 1 ) ) ) ) = std::numeric_limits<float>::infinity();
    const auto missing = mx::improc::PSFNoiseTraining::sample( science, 5, 5, 60, 40, config, exclusions );
    REQUIRE( missing.m_incomplete > 0 );
    REQUIRE( missing.m_samples.rows() < after.m_samples.rows() );
    REQUIRE( missing.m_samples.allFinite() );
    const auto edge = mx::improc::PSFNoiseTraining::sample( science, 5, 5, 80, 40, config, {} );
    REQUIRE( edge.m_incomplete > 0 );
    REQUIRE( edge.m_samples.allFinite() );
}

/** \brief Verify PSFNoiseTraining::sample retains disjoint stencils without reading any withheld native pixel.
 * \ingroup PSFNoiseTraining_unit_tests
 */
TEST_CASE( "Exact annular exclusions protect pixel support without enclosing-circle losses",
           "[PSFNoiseTraining][exclusion][stencil]" )
{
    imageT science = imageT::Ones( 81, 81 );
    mx::improc::PSFNoiseTrainingConfig config;
    config.m_guardRadius = 1.5;
    const std::vector<mx::improc::PSFNoiseExclusion> exclusions{ { 36.3, 21.7, 4.2 }, { 21, 43, 0 } };
    const auto conservative = mx::improc::PSFNoiseTraining::sample( science, 7, 11, 59, 44, config, exclusions );
    config.m_exactExclusion = true;
    const auto exact = mx::improc::PSFNoiseTraining::sample( science, 7, 11, 59, 44, config, exclusions );
    REQUIRE( exact.m_attempted == conservative.m_attempted );
    REQUIRE( exact.m_samples.rows() > conservative.m_samples.rows() );
    REQUIRE( exact.m_excluded < conservative.m_excluded );
    REQUIRE( exact.m_incomplete == 0 );
    // Poison the entire declared holdout, not just the candidate center or the rotated sample points.
    for( int column = 0; column < science.cols(); ++column )
    {
        for( int row = 0; row < science.rows(); ++row )
        {
            const int closestRow = std::clamp( row, 56, 62 );
            const int closestColumn = std::clamp( column, 39, 49 );
            if( std::hypot( row - closestRow, column - closestColumn ) <= 1.5 ||
                std::hypot( row - 36.3, column - 21.7 ) <= 4.2 || ( row == 21 && column == 43 ) )
            {
                science( row, column ) = std::numeric_limits<float>::quiet_NaN();
            }
        }
    }
    const auto poisoned = mx::improc::PSFNoiseTraining::sample( science, 7, 11, 59, 44, config, exclusions );
    REQUIRE( poisoned.m_samples == exact.m_samples );
    REQUIRE( poisoned.m_centers == exact.m_centers );
    REQUIRE( poisoned.m_excluded == exact.m_excluded );
    REQUIRE( poisoned.m_incomplete == 0 );
    REQUIRE( poisoned.m_samples.minCoeff() == 1 );
    REQUIRE( poisoned.m_samples.maxCoeff() == 1 );
    // A masked point outside the holdout still makes every stencil that reads it incomplete.
    science( static_cast<int>( std::floor( exact.m_centers( 0, 0 ) ) ),
             static_cast<int>( std::floor( exact.m_centers( 0, 1 ) ) ) ) = std::numeric_limits<float>::quiet_NaN();
    const auto incomplete = mx::improc::PSFNoiseTraining::sample( science, 7, 11, 59, 44, config, exclusions );
    REQUIRE( incomplete.m_excluded == exact.m_excluded );
    REQUIRE( incomplete.m_incomplete > 0 );
    REQUIRE( incomplete.m_samples.rows() < exact.m_samples.rows() );
    REQUIRE( incomplete.m_samples.allFinite() );

    // With integer-aligned samples, an unused zero-weight neighbor must not exclude the sample itself.
    science.setOnes();
    config.m_guardRadius = 0;
    config.m_arcStep = 10 * std::numbers::pi;
    const auto adjacent = mx::improc::PSFNoiseTraining::sample( science, 1, 1, 60, 40, config, { { 41, 60, 0 } } );
    const auto touching = mx::improc::PSFNoiseTraining::sample( science, 1, 1, 60, 40, config, { { 40, 60, 0 } } );
    REQUIRE( adjacent.m_attempted == 4 );
    REQUIRE( adjacent.m_samples.rows() == 3 );
    REQUIRE( touching.m_samples.rows() == 2 );
}

/** \brief Verify PSFNoiseTraining::calculate against a dense PCA covariance and source-isolated amplitude increments.
 * \ingroup PSFNoiseTraining_unit_tests
 */
TEST_CASE( "Annular covariance filtering matches an independent dense solve", "[PSFNoiseTraining][filter][dense]" )
{
    imageT science( 65, 65 );
    for( int column = 0; column < science.cols(); ++column )
    {
        for( int row = 0; row < science.rows(); ++row )
        {
            science( row, column ) = std::sin( 0.23 * row ) + std::cos( 0.31 * column ) + 0.03 * ( row % 7 );
        }
    }
    imageT response( 3, 3 );
    response << -0.25F, 0, 0.125F, -0.5F, 1, -0.25F, 0.125F, 0.25F, -0.125F;
    const auto validity = mx::improc::P4PSFFilter::validityT::Ones( 3, 3 ).eval();
    mx::improc::PSFNoiseTrainingConfig config;
    config.m_maximumModes = 2;
    const auto training = mx::improc::PSFNoiseTraining::sample( science, 3, 3, 48, 32, config, {} );
    const Eigen::VectorXd mean = training.m_samples.colwise().mean().transpose();
    const Eigen::MatrixXd centered = training.m_samples.rowwise() - mean.transpose();
    const Eigen::MatrixXd sampleCovariance = centered.transpose() * centered / ( centered.rows() - 1 );
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> decomposition( sampleCovariance );
    REQUIRE( decomposition.info() == Eigen::Success );
    for( const auto kind :
         { mx::improc::PSFNoiseKind::identity, mx::improc::PSFNoiseKind::diagonal, mx::improc::PSFNoiseKind::pca } )
    {
        config.m_kind = kind;
        const auto result =
            mx::improc::PSFNoiseTraining::calculate( science, response, validity, 48, 32, 1, config, {} );
        REQUIRE( result.m_status == mx::improc::PSFNoiseTrainingStatus::valid );
        REQUIRE( result.m_filter.valid );
        Eigen::MatrixXd covariance = Eigen::MatrixXd::Identity( 9, 9 );
        if( kind == mx::improc::PSFNoiseKind::diagonal )
        {
            covariance = sampleCovariance.diagonal().cwiseMax( result.m_varianceFloor ).asDiagonal();
        }
        if( kind == mx::improc::PSFNoiseKind::pca )
        {
            covariance *= result.m_varianceFloor;
            for( int mode = 7; mode < 9; ++mode )
            {
                const Eigen::VectorXd u = decomposition.eigenvectors().col( mode );
                covariance +=
                    std::max( 0.0, decomposition.eigenvalues()( mode ) - result.m_varianceFloor ) * u * u.transpose();
            }
        }
        Eigen::VectorXd signal( 9 ), data( 9 );
        for( int column = 0; column < 3; ++column )
        {
            for( int row = 0; row < 3; ++row )
            {
                const int pixel = row + 3 * column;
                signal( pixel ) = response( row, column );
                data( pixel ) = science( 47 + row, 31 + column ) -
                                ( kind == mx::improc::PSFNoiseKind::identity ? 0 : mean( pixel ) );
            }
        }
        const Eigen::VectorXd weight = covariance.llt().solve( signal );
        REQUIRE( result.m_filter.amplitude == Approx( weight.dot( data ) / weight.dot( signal ) ).margin( 1e-11 ) );
        imageT injected = science;
        injected.block( 47, 31, 3, 3 ) += 2 * response;
        const auto added =
            mx::improc::PSFNoiseTraining::calculate( injected, response, validity, 48, 32, 1, config, {} );
        REQUIRE( added.m_filter.valid );
        REQUIRE( added.m_varianceFloor == result.m_varianceFloor );
        REQUIRE( added.m_samples == result.m_samples );
        REQUIRE( added.m_filter.conditionalSigma == result.m_filter.conditionalSigma );
        REQUIRE( added.m_filter.amplitude - result.m_filter.amplitude == Approx( 2 ).margin( 2e-7 ) );
    }
}

/** \brief Verify PSFNoiseTraining rejects invalid policies and reports insufficient or degenerate training explicitly.
 * \ingroup PSFNoiseTraining_unit_tests
 */
TEST_CASE( "Annular covariance training validates its contract", "[PSFNoiseTraining][validation]" )
{
    const imageT science = imageT::Zero( 31, 31 );
    const imageT response = imageT::Ones( 3, 3 );
    const auto validity = mx::improc::P4PSFFilter::validityT::Ones( 3, 3 ).eval();
    mx::improc::PSFNoiseTrainingConfig config;
    REQUIRE( mx::improc::PSFNoiseTraining::calculate( science, response, validity, 25, 15, 1, config, {} ).m_status ==
             mx::improc::PSFNoiseTrainingStatus::zeroVariance );
    REQUIRE( mx::improc::PSFNoiseTraining::calculate( science, response, validity, 15, 15, 1, config, {} ).m_status ==
             mx::improc::PSFNoiseTrainingStatus::tooFewSamples );
    REQUIRE(
        mx::improc::PSFNoiseTraining::calculate( imageT::Ones( 31, 31 ), response, validity, 25, 15, 1, config, {} )
            .m_status == mx::improc::PSFNoiseTrainingStatus::zeroVariance );
    const auto noSupport = mx::improc::PSFNoiseTraining::calculate( science,
                                                                    response,
                                                                    mx::improc::P4PSFFilter::validityT::Zero( 3, 3 ),
                                                                    25,
                                                                    15,
                                                                    1,
                                                                    config,
                                                                    {} );
    REQUIRE( noSupport.m_status == mx::improc::PSFNoiseTrainingStatus::invalidFilter );
    REQUIRE( noSupport.m_attempted == 0 );
    config.m_minimumSamples = 10000;
    REQUIRE( mx::improc::PSFNoiseTraining::calculate( science, response, validity, 25, 15, 1, config, {} ).m_status ==
             mx::improc::PSFNoiseTrainingStatus::tooFewSamples );
    config.m_minimumSamples = 1;
    REQUIRE_THROWS_AS( mx::improc::PSFNoiseTraining::sample( science, 3, 3, 25, 15, config, {} ),
                       std::invalid_argument );
    config.m_minimumSamples = 8;
    config.m_floorFraction = 0;
    REQUIRE_THROWS_AS( mx::improc::PSFNoiseTraining::sample( science, 3, 3, 25, 15, config, {} ),
                       std::invalid_argument );
    config.m_floorFraction = 0.1;
    REQUIRE_THROWS_AS( mx::improc::PSFNoiseTraining::sample( science, 2, 3, 25, 15, config, {} ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( mx::improc::PSFNoiseTraining::sample( science, 3, 3, -1, 15, config, {} ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( mx::improc::PSFNoiseTraining::sample( science, 3, 3, 25, 15, config, { { 0, 0, -1 } } ),
                       std::invalid_argument );
    config.m_arcStep = std::numeric_limits<double>::min();
    REQUIRE_THROWS_AS( mx::improc::PSFNoiseTraining::sample( science, 3, 3, 25, 15, config, {} ),
                       std::invalid_argument );
    config.m_arcStep = 1;
    config.m_guardRadius = std::numeric_limits<double>::quiet_NaN();
    REQUIRE_THROWS_AS( mx::improc::PSFNoiseTraining::sample( science, 3, 3, 25, 15, config, {} ),
                       std::invalid_argument );
}

} // namespace PSFNoiseTraining_test
} // namespace unitTest
