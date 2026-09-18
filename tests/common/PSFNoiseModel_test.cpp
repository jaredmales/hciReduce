/** \file PSFNoiseModel_test.cpp
 * \brief Tests residual-noise estimation and marginal covariance solves against dense references.
 */

#include "../catch2/catch.hpp"

#include "src/common/PSFNoiseModel.hpp"

#include <limits>
#include <numeric>
#include <random>

namespace unitTest
{
namespace PSFNoiseModel_test
{

/// Double-precision pixel vector used by the dense oracles.
using vectorT = Eigen::VectorXd;

/// Double-precision covariance and training matrix used by the dense oracles.
using matrixT = Eigen::MatrixXd;

/** \brief Verify PSFNoiseModel identity, diagonal, and lowRank solves against dense Cholesky for changing support.
 * \ingroup PSFNoiseModel_unit_tests
 */
TEST_CASE( "Noise covariance solves agree with dense marginal covariances", "[PSFNoiseModel][solve][support]" )
{
    std::mt19937 random( 417 );
    std::normal_distribution<double> normal;
    constexpr Eigen::Index pixels = 15;
    const vectorT variances = vectorT::LinSpaced( pixels, 0.2, 2.0 );
    vectorT value( pixels );
    matrixT factor( pixels, 4 );
    for( Eigen::Index row = 0; row < pixels; ++row )
    {
        value( row ) = normal( random );
        for( Eigen::Index column = 0; column < factor.cols(); ++column )
        {
            factor( row, column ) = normal( random );
        }
    }
    // Include dependent, nonorthogonal modes. Orthogonality cannot be assumed after masking.
    factor.col( 3 ) = 0.7 * factor.col( 0 ) + 0.2 * factor.col( 1 );
    std::vector<Eigen::Index> full( pixels );
    std::iota( full.begin(), full.end(), 0 );
    const std::vector<std::vector<Eigen::Index>> supports{ full, { 14, 0, 5, 8, 2, 11 }, { 7 }, {}, full };
    for( const double factorScale : { 0.0, 1.0, 100.0 } )
    {
        const matrixT scaledFactor = factorScale * factor;
        const mx::improc::PSFNoiseModel model = mx::improc::PSFNoiseModel::lowRank( variances, scaledFactor );
        const mx::improc::PSFNoiseModel diagonal = mx::improc::PSFNoiseModel::diagonal( variances );
        const mx::improc::PSFNoiseModel identity = mx::improc::PSFNoiseModel::identity( pixels );
        matrixT covariance = variances.asDiagonal();
        covariance += scaledFactor * scaledFactor.transpose();
        REQUIRE( model.kind() == mx::improc::PSFNoiseKind::pca );
        REQUIRE( model.pixels() == pixels );
        REQUIRE( model.modes() == 4 );
        for( const auto &support : supports )
        {
            CAPTURE( factorScale, support.size() );
            const Eigen::Index count = static_cast<Eigen::Index>( support.size() );
            matrixT marginal( count, count );
            vectorT rhs( count );
            vectorT diagonalExpected( count );
            for( Eigen::Index row = 0; row < count; ++row )
            {
                rhs( row ) = value( support[row] );
                diagonalExpected( row ) = rhs( row ) / variances( support[row] );
                for( Eigen::Index column = 0; column < count; ++column )
                {
                    marginal( row, column ) = covariance( support[row], support[column] );
                }
            }
            REQUIRE( identity.solve( rhs, support ).isApprox( rhs, 1e-14 ) );
            REQUIRE( diagonal.solve( rhs, support ).isApprox( diagonalExpected, 1e-14 ) );
            const vectorT actual = model.solve( rhs, support );
            if( count == 0 )
            {
                REQUIRE( actual.size() == 0 );
                continue;
            }
            const Eigen::LLT<matrixT> dense( marginal );
            REQUIRE( dense.info() == Eigen::Success );
            const vectorT expected = dense.solve( rhs );
            REQUIRE( ( actual - expected ).norm() <= 2e-9 * expected.norm() );
            REQUIRE( ( marginal * actual - rhs ).norm() <= 2e-9 * rhs.norm() );
        }
    }
}

/** \brief Verify PSFNoiseModel::solve marginalizes missing pixels rather than cropping the full precision matrix.
 * \ingroup PSFNoiseModel_unit_tests
 */
TEST_CASE( "Missing noise pixels are marginalized", "[PSFNoiseModel][marginalization]" )
{
    const vectorT variances = vectorT::Ones( 2 );
    const matrixT factor = matrixT::Ones( 2, 1 );
    const mx::improc::PSFNoiseModel model = mx::improc::PSFNoiseModel::lowRank( variances, factor );
    const vectorT actual = model.solve( vectorT::Ones( 1 ), { 0 } );
    // C = [[2,1],[1,2]], so inverse(C_subset)=1/2, while subset(inverse(C))=2/3.
    REQUIRE( actual( 0 ) == Approx( 0.5 ).epsilon( 1e-13 ) );
    REQUIRE( std::abs( actual( 0 ) - 2.0 / 3.0 ) > 0.1 );
}

/** \brief Verify PSFNoiseModel sample centering, n-1 normalization, PCA truncation, and finite complement variance.
 * \ingroup PSFNoiseModel_unit_tests
 */
TEST_CASE( "Noise training retains the mean and regularizes the PCA complement", "[PSFNoiseModel][training][PCA]" )
{
    matrixT samples( 4, 3 );
    samples << 13, 21, 30, 13, 19, 30, 7, 21, 30, 7, 19, 30;
    const vectorT expectedMean = ( vectorT( 3 ) << 10, 20, 30 ).finished();
    const matrixT rotation = ( matrixT( 3, 3 ) << 0.8, -0.6, 0, 0.6, 0.8, 0, 0, 0, 1 ).finished();
    const matrixT rotatedSamples = samples * rotation.transpose();
    const std::vector<Eigen::Index> support{ 0, 1, 2 };
    const vectorT value = ( vectorT( 3 ) << 2, -3, 4 ).finished();
    const mx::improc::PSFNoiseModel diagonal = mx::improc::PSFNoiseModel::estimateDiagonal( samples, 0.25 );
    const vectorT diagonalExpected = ( vectorT( 3 ) << 2.0 / 12, -3.0 / ( 4.0 / 3 ), 4 / 0.25 ).finished();
    REQUIRE( diagonal.kind() == mx::improc::PSFNoiseKind::diagonal );
    REQUIRE( diagonal.mean().isApprox( expectedMean, 1e-14 ) );
    REQUIRE( diagonal.solve( value, support ).isApprox( diagonalExpected, 1e-14 ) );
    for( const double floor : { 0.25, 2.0, 20.0 } )
    {
        for( const std::size_t modes : { 0U, 1U, 2U, 100U } )
        {
            CAPTURE( floor, modes );
            const mx::improc::PSFNoiseModel model =
                mx::improc::PSFNoiseModel::estimatePCA( rotatedSamples, modes, floor );
            vectorT spectrum = vectorT::Constant( 3, floor );
            Eigen::Index retained = 0;
            if( modes >= 1 && 12 > floor )
            {
                spectrum( 0 ) = 12;
                ++retained;
            }
            if( modes >= 2 && 4.0 / 3 > floor )
            {
                spectrum( 1 ) = 4.0 / 3;
                ++retained;
            }
            const matrixT covariance = rotation * spectrum.asDiagonal() * rotation.transpose();
            const vectorT expected = covariance.llt().solve( value );
            REQUIRE( model.modes() == retained );
            REQUIRE( model.mean().isApprox( rotation * expectedMean, 1e-14 ) );
            REQUIRE( model.solve( value, support ).isApprox( expected, 1e-12 ) );
        }
    }
    const matrixT constant = matrixT::Constant( 3, 5, 8 );
    const mx::improc::PSFNoiseModel constantModel = mx::improc::PSFNoiseModel::estimatePCA( constant, 100, 0.5 );
    REQUIRE( constantModel.modes() == 0 );
    REQUIRE( constantModel.mean().isApprox( vectorT::Constant( 5, 8 ), 1e-14 ) );
    REQUIRE(
        constantModel.solve( vectorT::Ones( 5 ), { 0, 1, 2, 3, 4 } ).isApprox( vectorT::Constant( 5, 2 ), 1e-14 ) );
    const matrixT twoSamples = ( matrixT( 2, 4 ) << 1, 2, 3, 4, -1, -2, -3, -4 ).finished();
    const mx::improc::PSFNoiseModel rankOne = mx::improc::PSFNoiseModel::estimatePCA( twoSamples, 100, 0.25 );
    REQUIRE( rankOne.modes() == 1 );
    const vectorT direction = twoSamples.row( 0 ).transpose();
    REQUIRE( rankOne.solve( direction, { 0, 1, 2, 3 } ).isApprox( direction / 60, 1e-12 ) );
}

/** \brief Verify PSFNoiseModel::solve preserves finite precision along strongly downweighted modes.
 * \ingroup PSFNoiseModel_unit_tests
 */
TEST_CASE( "Low-rank solves retain precision in high-variance directions", "[PSFNoiseModel][precision]" )
{
    const Eigen::VectorXd variance = Eigen::VectorXd::Ones( 3 );
    Eigen::MatrixXd factor = Eigen::MatrixXd::Zero( 3, 1 );
    factor( 0, 0 ) = 1e8;
    const mx::improc::PSFNoiseModel model = mx::improc::PSFNoiseModel::lowRank( variance, factor );
    const Eigen::VectorXd parallel = ( Eigen::VectorXd( 3 ) << 1, 0, 0 ).finished();
    const Eigen::VectorXd actual = model.solve( parallel, { 0, 1, 2 } );
    REQUIRE( actual( 0 ) == Approx( 1e-16 ).epsilon( 1e-13 ) );
    REQUIRE( actual.tail( 2 ).isZero() );
    REQUIRE( model.solve( Eigen::VectorXd::Ones( 1 ), { 0 } )( 0 ) == Approx( 1e-16 ).epsilon( 1e-13 ) );
}

/** \brief Verify PSFNoiseModel rejects malformed covariance, training, support, and nonfinite arithmetic.
 * \ingroup PSFNoiseModel_unit_tests
 */
TEST_CASE( "Noise model validates covariance and training inputs", "[PSFNoiseModel][validation]" )
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    const vectorT ones = vectorT::Ones( 3 );
    REQUIRE_THROWS_AS( mx::improc::PSFNoiseModel::identity( 0 ), std::invalid_argument );
    REQUIRE_THROWS_AS( mx::improc::PSFNoiseModel::identity( -1 ), std::invalid_argument );
    REQUIRE_THROWS_AS( mx::improc::PSFNoiseModel::diagonal( vectorT{} ), std::invalid_argument );
    for( const double bad : { 0.0, -1.0, nan, inf } )
    {
        REQUIRE_THROWS_AS( mx::improc::PSFNoiseModel::diagonal( vectorT::Constant( 3, bad ) ), std::invalid_argument );
        REQUIRE_THROWS_AS( mx::improc::PSFNoiseModel::estimateDiagonal( matrixT::Ones( 3, 3 ), bad ),
                           std::invalid_argument );
        REQUIRE_THROWS_AS( mx::improc::PSFNoiseModel::estimatePCA( matrixT::Ones( 3, 3 ), 2, bad ),
                           std::invalid_argument );
    }
    REQUIRE_THROWS_AS( mx::improc::PSFNoiseModel::lowRank( ones, matrixT::Ones( 2, 2 ) ), std::invalid_argument );
    REQUIRE_THROWS_AS( mx::improc::PSFNoiseModel::lowRank( ones, matrixT::Constant( 3, 1, inf ) ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( mx::improc::PSFNoiseModel::diagonal( ones, vectorT::Ones( 2 ) ), std::invalid_argument );
    REQUIRE_THROWS_AS( mx::improc::PSFNoiseModel::diagonal( ones, vectorT::Constant( 3, nan ) ),
                       std::invalid_argument );
    for( const matrixT &bad : std::vector<matrixT>{ matrixT( 1, 3 ), matrixT( 3, 0 ), matrixT::Constant( 3, 3, nan ) } )
    {
        REQUIRE_THROWS_AS( mx::improc::PSFNoiseModel::estimateDiagonal( bad, 1 ), std::invalid_argument );
        REQUIRE_THROWS_AS( mx::improc::PSFNoiseModel::estimatePCA( bad, 1, 1 ), std::invalid_argument );
    }
    const mx::improc::PSFNoiseModel model = mx::improc::PSFNoiseModel::identity( 3 );
    REQUIRE( model.kind() == mx::improc::PSFNoiseKind::identity );
    REQUIRE( model.modes() == 0 );
    REQUIRE( model.mean().isZero() );
    REQUIRE_THROWS_AS( model.solve( ones, { 0, 1 } ), std::invalid_argument );
    REQUIRE_THROWS_AS( model.solve( ones, { 0, 1, 1 } ), std::invalid_argument );
    REQUIRE_THROWS_AS( model.solve( ones, { -1, 1, 2 } ), std::invalid_argument );
    REQUIRE_THROWS_AS( model.solve( ones, { 0, 1, 3 } ), std::invalid_argument );
    REQUIRE_THROWS_AS( model.solve( vectorT::Constant( 3, nan ), { 0, 1, 2 } ), std::invalid_argument );
    const matrixT huge = ( matrixT( 2, 1 ) << 1e200, -1e200 ).finished();
    REQUIRE_THROWS_AS( mx::improc::PSFNoiseModel::estimateDiagonal( huge, 1 ), std::overflow_error );
    REQUIRE_THROWS_AS( mx::improc::PSFNoiseModel::estimatePCA( huge, 1, 1 ), std::overflow_error );
    const mx::improc::PSFNoiseModel overflowing =
        mx::improc::PSFNoiseModel::lowRank( ones, matrixT::Constant( 3, 1, 1e200 ) );
    REQUIRE_THROWS_AS( overflowing.solve( ones, { 0, 1, 2 } ), std::overflow_error );
    const mx::improc::PSFNoiseModel tinyVariance =
        mx::improc::PSFNoiseModel::diagonal( vectorT::Constant( 3, 1e-300 ) );
    REQUIRE_THROWS_AS( tinyVariance.solve( vectorT::Constant( 3, 1e100 ), { 0, 1, 2 } ), std::overflow_error );
}

} // namespace PSFNoiseModel_test
} // namespace unitTest
