/** \file P4PCA_test.cpp
 * \brief Tests the pure numerical predictor used by Pixel Prediction Post-Processing.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "src/common/P4PCA.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace unitTest
{
namespace P4PCA_test
{

/** \cond P4PCA_test_harness */
using pcaT = mx::improc::P4PCA;
using resultT = mx::improc::P4PCAResult;
using statusT = mx::improc::P4PCAModeStatus;

/// Compare two Eigen-like vectors coefficient by coefficient.
template <typename actualT, typename expectedT>
void requireApprox( const actualT &actual,     /**< [in] values produced by P4PCA */
                    const expectedT &expected, /**< [in] expected reference values */
                    double tolerance = 1e-10 /**< [in] absolute comparison tolerance */ )
{
    REQUIRE( actual.size() == expected.size() );
    for( Eigen::Index index = 0; index < actual.size(); ++index )
    {
        REQUIRE( actual( index ) == Approx( expected( index ) ).margin( tolerance ) );
    }
}

/// Return a direct thin-SVD residual using the requested largest singular modes.
pcaT::vectorT svdResidual( const pcaT::matrixT &predictors, /**< [in] predictor matrix */
                           const pcaT::vectorT &target,     /**< [in] target time series */
                           int modeCount /**< [in] number of largest singular modes */ )
{
    Eigen::JacobiSVD<Eigen::MatrixXd> decomposition( predictors.matrix(), Eigen::ComputeThinU | Eigen::ComputeThinV );
    const Eigen::MatrixXd basis = decomposition.matrixU().leftCols( modeCount );
    return ( target.matrix() - basis * ( basis.transpose() * target.matrix() ) ).array();
}

/// Return a copy with each column's temporal mean removed.
template <typename arrayT>
arrayT centeredColumns( const arrayT &input /**< [in] values to center */ )
{
    arrayT centered = input;
    for( Eigen::Index column = 0; column < centered.cols(); ++column )
    {
        centered.col( column ) -= centered.col( column ).mean();
    }
    return centered;
}

/// Return direct thin-SVD coefficients for a centered fit using the requested largest singular modes.
pcaT::vectorT centeredSvdCoefficients( const pcaT::matrixT &predictors, /**< [in] uncentered predictor matrix */
                                       const pcaT::vectorT &target,     /**< [in] uncentered target time series */
                                       int modeCount /**< [in] number of largest centered singular modes */ )
{
    const pcaT::matrixT centeredPredictors = centeredColumns( predictors );
    const pcaT::vectorT centeredTarget = centeredColumns( target );
    Eigen::JacobiSVD<Eigen::MatrixXd> decomposition( centeredPredictors.matrix(),
                                                     Eigen::ComputeThinU | Eigen::ComputeThinV );
    const Eigen::MatrixXd temporalBasis = decomposition.matrixU().leftCols( modeCount );
    const Eigen::MatrixXd predictorBasis = decomposition.matrixV().leftCols( modeCount );
    return ( predictorBasis * decomposition.singularValues().head( modeCount ).cwiseInverse().asDiagonal() *
             temporalBasis.transpose() * centeredTarget.matrix() )
        .array();
}

/// Return the direct thin-SVD residual for a centered fit applied to the uncentered predictors.
pcaT::vectorT centeredFitUncenteredResidual( const pcaT::matrixT &predictors, /**< [in] predictor matrix */
                                             const pcaT::vectorT &target,     /**< [in] target time series */
                                             int modeCount /**< [in] number of largest centered singular modes */ )
{
    const pcaT::vectorT coefficients = centeredSvdCoefficients( predictors, target, modeCount );
    return ( target.matrix() - predictors.matrix() * coefficients.matrix() ).array();
}

/// Determine whether an Eigen-like array contains only finite values under fast-math.
template <typename arrayT>
bool allFinite( const arrayT &array /**< [in] values to inspect */ )
{
    for( Eigen::Index column = 0; column < array.cols(); ++column )
    {
        for( Eigen::Index row = 0; row < array.rows(); ++row )
        {
            if( !mx::math::isFinite( array( row, column ) ) )
            {
                return false;
            }
        }
    }
    return true;
}

/// Determine whether an Eigen-like array contains only NaNs under fast-math.
template <typename arrayT>
bool allNan( const arrayT &array /**< [in] values to inspect */ )
{
    for( Eigen::Index column = 0; column < array.cols(); ++column )
    {
        for( Eigen::Index row = 0; row < array.rows(); ++row )
        {
            if( !mx::math::isNan( array( row, column ) ) )
            {
                return false;
            }
        }
    }
    return true;
}

/// Controlled eigensolver outcomes used to exercise P4PCA's post-solver validation.
enum class fakeSolverBehavior
{
    valid,
    failure,
    invalidDimensions,
    invalidEigenvalueShape,
    nonfiniteEigenvalue,
    nonfiniteEigenvector,
    unsortedEigenvalues,
    nonpositiveEigenvalues,
    tinyEigenvalue,
    residualOverflow
};

/// Outcome supplied by the controlled eigensolver.
fakeSolverBehavior solverBehavior{ fakeSolverBehavior::failure };

/// Number of invocations of the controlled eigensolver.
int solverCalls{ 0 };

/// Gram matrix received by the controlled eigensolver.
pcaT::matrixT solverGram;

/// Eigenpair count received by the controlled eigensolver.
int solverModeCount{ 0 };

/// Supply deterministic eigensolver failures and malformed results to P4PCA.
MXLAPACK_INT fakeEigenSolver( pcaT::matrixT &eigenvectors, /**< [out] controlled eigenvectors */
                              pcaT::matrixT &eigenvalues,  /**< [out] controlled eigenvalues */
                              pcaT::matrixT &covariance,   /**< [in] covariance matrix */
                              int modeCount,               /**< [in] requested largest eigenpair count */
                              pcaT::workspaceT &workspace /**< [in,out] unused caller workspace */ )
{
    static_cast<void>( workspace );
    ++solverCalls;
    solverGram = covariance;
    solverModeCount = modeCount;

    if( solverBehavior == fakeSolverBehavior::failure )
    {
        return 37;
    }

    if( solverBehavior == fakeSolverBehavior::invalidDimensions )
    {
        eigenvectors.resize( 0, 0 );
        eigenvalues.resize( 0, 0 );
        return 0;
    }

    const int predictorCount = covariance.rows();
    const int firstEigenvector = predictorCount - modeCount;
    eigenvectors.resize( predictorCount, modeCount );
    eigenvectors.setZero();
    eigenvalues.resize( predictorCount, 1 );
    eigenvalues.setZero();

    for( int mode = 0; mode < modeCount; ++mode )
    {
        eigenvectors( firstEigenvector + mode, mode ) = 1;
        eigenvalues( mode ) = firstEigenvector + mode + 1;
    }

    if( solverBehavior == fakeSolverBehavior::invalidEigenvalueShape )
    {
        eigenvalues.transposeInPlace();
        return 0;
    }

    if( solverBehavior == fakeSolverBehavior::nonfiniteEigenvalue )
    {
        eigenvalues( 0 ) = std::numeric_limits<double>::quiet_NaN();
    }
    else if( solverBehavior == fakeSolverBehavior::nonfiniteEigenvector )
    {
        eigenvectors( 0, 0 ) = std::numeric_limits<double>::infinity();
    }
    else if( solverBehavior == fakeSolverBehavior::unsortedEigenvalues )
    {
        eigenvalues( 0 ) = 2;
        eigenvalues( 1 ) = 1;
    }
    else if( solverBehavior == fakeSolverBehavior::nonpositiveEigenvalues )
    {
        eigenvalues( 0 ) = -2;
        eigenvalues( 1 ) = -1;
    }
    else if( solverBehavior == fakeSolverBehavior::tinyEigenvalue )
    {
        eigenvalues( 0 ) = std::numeric_limits<double>::denorm_min();
        eigenvalues( 1 ) = 1;
    }
    else if( solverBehavior == fakeSolverBehavior::residualOverflow )
    {
        eigenvectors.col( 1 ) << 0.5, -1;
    }

    return 0;
}

/// Restore the production eigensolver after a controlled test scope.
struct solverReset
{
    /// Install the controlled eigensolver and reset its call count.
    explicit solverReset( fakeSolverBehavior behavior /**< [in] controlled outcome */ )
    {
        solverBehavior = behavior;
        solverCalls = 0;
        solverGram.resize( 0, 0 );
        solverModeCount = 0;
        mx::improc::detail::p4PCASetEigenSolverForTesting( &fakeEigenSolver );
    }

    /// Restore the production eigensolver.
    ~solverReset()
    {
        mx::improc::detail::p4PCAResetEigenSolverForTesting();
    }
};
/** \endcond */

/// Verify P4PCA's real normal-equation solve produces exact uncentered truncated residuals.
/** \ingroup P4PCA_unit_tests */
TEST_CASE( "P4PCA computes exact truncated residuals", "[P4PCA][regression]" )
{
    pcaT::matrixT predictors( 3, 2 );
    predictors << 1, 0, 0, 2, 0, 0;
    pcaT::vectorT target( 3 );
    target << 3, 5, 7;

    resultT result;
    pcaT::workspaceT workspace;
    mx::improc::P4PCA::calculate( result, predictors, target, { 1, 2 }, 1e-12, workspace );

    REQUIRE( result.residuals.rows() == 3 );
    REQUIRE( result.residuals.cols() == 2 );
    REQUIRE( result.numericalRank == 2 );
    REQUIRE( result.modeStatus == std::vector<statusT>{ statusT::rankSupported, statusT::rankSupported } );

    pcaT::vectorT oneMode( 3 );
    oneMode << 3, 0, 7;
    pcaT::vectorT twoModes( 3 );
    twoModes << 0, 0, 7;
    requireApprox( result.residuals.col( 0 ), oneMode );
    requireApprox( result.residuals.col( 1 ), twoModes );
    REQUIRE( result.residuals.col( 1 ).matrix().norm() <= result.residuals.col( 0 ).matrix().norm() );
}

/// Verify P4PCA sends the smaller shape-dependent Gram matrix to its eigensolver seam.
/** \ingroup P4PCA_unit_tests */
TEST_CASE( "P4PCA selects the smaller Gram matrix", "[P4PCA][regression][solver]" )
{
    resultT result;
    pcaT::workspaceT workspace;

    SECTION( "T is less than K" )
    {
        pcaT::matrixT predictors( 2, 3 );
        predictors << 1, 2, 3, 4, 5, 6;
        pcaT::vectorT target( 2 );
        target << 7, 8;
        const pcaT::matrixT expected = ( predictors.matrix() * predictors.matrix().transpose() ).array();

        solverReset reset( fakeSolverBehavior::failure );
        REQUIRE_THROWS_WITH( mx::improc::P4PCA::calculate( result, predictors, target, { 1 }, 0, workspace ),
                             "P4PCA eigensolver failed with status 37" );
        REQUIRE( solverModeCount == 2 );
        requireApprox( solverGram, expected );
    }

    SECTION( "T equals K" )
    {
        pcaT::matrixT predictors( 2, 2 );
        predictors << 1, 2, 0, 3;
        pcaT::vectorT target( 2 );
        target << 4, 5;
        const pcaT::matrixT expected = ( predictors.matrix() * predictors.matrix().transpose() ).array();
        const pcaT::matrixT predictorGram = ( predictors.matrix().transpose() * predictors.matrix() ).array();
        REQUIRE( ( expected - predictorGram ).matrix().norm() > 1 );

        solverReset reset( fakeSolverBehavior::failure );
        REQUIRE_THROWS_WITH( mx::improc::P4PCA::calculate( result, predictors, target, { 1 }, 0, workspace ),
                             "P4PCA eigensolver failed with status 37" );
        REQUIRE( solverModeCount == 2 );
        requireApprox( solverGram, expected );
    }

    SECTION( "T exceeds K" )
    {
        pcaT::matrixT predictors( 3, 2 );
        predictors << 1, 2, 3, 4, 5, 6;
        pcaT::vectorT target( 3 );
        target << 7, 8, 9;
        const pcaT::matrixT expected = ( predictors.matrix().transpose() * predictors.matrix() ).array();

        solverReset reset( fakeSolverBehavior::failure );
        REQUIRE_THROWS_WITH( mx::improc::P4PCA::calculate( result, predictors, target, { 1 }, 0, workspace ),
                             "P4PCA eigensolver failed with status 37" );
        REQUIRE( solverModeCount == 2 );
        requireApprox( solverGram, expected );
    }
}

/// Verify P4PCA agrees with direct SVD for tall, wide, and orthogonally transformed predictors.
/** \ingroup P4PCA_unit_tests */
TEST_CASE( "P4PCA agrees with direct SVD across matrix shapes", "[P4PCA][reference]" )
{
    pcaT::workspaceT workspace;

    SECTION( "T exceeds K" )
    {
        pcaT::matrixT predictors( 4, 2 );
        predictors << 1, 0, 0, 2, 1, 0, 0, 0;
        pcaT::vectorT target( 4 );
        target << 3, 5, -1, 7;

        resultT result;
        mx::improc::P4PCA::calculate( result, predictors, target, { 1, 2 }, 1e-12, workspace );
        requireApprox( result.residuals.col( 0 ), svdResidual( predictors, target, 1 ) );
        requireApprox( result.residuals.col( 1 ), svdResidual( predictors, target, 2 ) );
    }

    SECTION( "T is less than K" )
    {
        pcaT::matrixT predictors( 3, 5 );
        predictors << 3, 0.2, 0, 1, -0.5, 0.1, 2, 0.4, -0.3, 0.8, 0, -0.2, 1, 0.7, 2;
        pcaT::vectorT target( 3 );
        target << 5, -1, 2;

        resultT result;
        mx::improc::P4PCA::calculate( result, predictors, target, { 1, 2, 3 }, 1e-12, workspace );
        REQUIRE( result.numericalRank == 3 );
        requireApprox( result.residuals.col( 0 ), svdResidual( predictors, target, 1 ) );
        requireApprox( result.residuals.col( 1 ), svdResidual( predictors, target, 2 ) );
        requireApprox( result.residuals.col( 2 ), svdResidual( predictors, target, 3 ) );
    }

    SECTION( "T equals K" )
    {
        pcaT::matrixT predictors( 3, 3 );
        predictors << 3, 0.2, 1, 0.1, 2, -0.3, 0.7, -0.2, 1;
        pcaT::vectorT target( 3 );
        target << 5, -1, 2;

        resultT result;
        mx::improc::P4PCA::calculate( result, predictors, target, { 1, 2, 3 }, 1e-12, workspace );
        REQUIRE( result.numericalRank == 3 );
        requireApprox( result.residuals.col( 0 ), svdResidual( predictors, target, 1 ) );
        requireApprox( result.residuals.col( 1 ), svdResidual( predictors, target, 2 ) );
        requireApprox( result.residuals.col( 2 ), svdResidual( predictors, target, 3 ) );
    }

    SECTION( "orthogonal predictor basis" )
    {
        pcaT::matrixT predictors( 4, 2 );
        predictors << 1, 0, 0, 2, 1, 0, 0, 0;
        pcaT::vectorT target( 4 );
        target << 3, 5, -1, 7;

        const double angle = 0.6;
        Eigen::Matrix2d rotation;
        rotation << std::cos( angle ), -std::sin( angle ), std::sin( angle ), std::cos( angle );
        const pcaT::matrixT rotatedPredictors = ( predictors.matrix() * rotation ).array();

        resultT original;
        resultT rotated;
        pcaT::workspaceT rotatedWorkspace;
        mx::improc::P4PCA::calculate( original, predictors, target, { 1, 2 }, 1e-12, workspace );
        mx::improc::P4PCA::calculate( rotated, rotatedPredictors, target, { 1, 2 }, 1e-12, rotatedWorkspace );
        REQUIRE( rotated.numericalRank == original.numericalRank );
        requireApprox( rotated.residuals, original.residuals );
    }

    SECTION( "orthogonal wide-predictor basis" )
    {
        pcaT::matrixT predictors( 3, 4 );
        predictors << 3, 0.2, 1, -0.5, 0.1, 2, -0.3, 0.8, 0.7, -0.2, 1, 2;
        pcaT::vectorT target( 3 );
        target << 5, -1, 2;

        const double angle = 0.6;
        Eigen::Matrix4d rotation = Eigen::Matrix4d::Identity();
        rotation( 0, 0 ) = std::cos( angle );
        rotation( 0, 3 ) = -std::sin( angle );
        rotation( 3, 0 ) = std::sin( angle );
        rotation( 3, 3 ) = std::cos( angle );
        const pcaT::matrixT rotatedPredictors = ( predictors.matrix() * rotation ).array();

        resultT original;
        resultT rotated;
        pcaT::workspaceT rotatedWorkspace;
        mx::improc::P4PCA::calculate( original, predictors, target, { 1, 2, 3 }, 1e-12, workspace );
        mx::improc::P4PCA::calculate( rotated, rotatedPredictors, target, { 1, 2, 3 }, 1e-12, rotatedWorkspace );
        REQUIRE( rotated.numericalRank == original.numericalRank );
        requireApprox( rotated.residuals, original.residuals );
    }

    SECTION( "non-contiguous requested counts" )
    {
        pcaT::matrixT predictors( 3, 5 );
        predictors << 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 3, 0, 0;
        pcaT::vectorT target( 3 );
        target << 4, 5, 6;

        resultT result;
        mx::improc::P4PCA::calculate( result, predictors, target, { 1, 3 }, 1e-12, workspace );
        REQUIRE( result.numericalRank == 3 );
        requireApprox( result.residuals.col( 0 ), svdResidual( predictors, target, 1 ) );
        requireApprox( result.residuals.col( 1 ), svdResidual( predictors, target, 3 ) );
    }

    SECTION( "repeated nonzero singular values" )
    {
        pcaT::matrixT predictors( 3, 5 );
        predictors << 3, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0;
        pcaT::vectorT target( 3 );
        target << 4, 5, 6;

        resultT result;
        mx::improc::P4PCA::calculate( result, predictors, target, { 1, 3 }, 1e-12, workspace );
        REQUIRE( result.numericalRank == 3 );
        requireApprox( result.residuals.col( 0 ), svdResidual( predictors, target, 1 ) );
        requireApprox( result.residuals.col( 1 ), svdResidual( predictors, target, 3 ) );
    }
}

/// Verify P4PCA preserves supported lower planes and marks only rank-insufficient planes as NaN.
/** \ingroup P4PCA_unit_tests */
TEST_CASE( "P4PCA marks rank-insufficient planes", "[P4PCA][rank]" )
{
    pcaT::workspaceT workspace;

    SECTION( "exact rank deficiency" )
    {
        pcaT::matrixT predictors( 3, 2 );
        predictors << 1, 2, 2, 4, 3, 6;
        pcaT::vectorT target( 3 );
        target << 1, 0, 1;

        resultT result;
        mx::improc::P4PCA::calculate( result, predictors, target, { 1, 2 }, 1e-10, workspace );
        REQUIRE( result.numericalRank == 1 );
        REQUIRE( result.modeStatus[0] == statusT::rankSupported );
        REQUIRE( result.modeStatus[1] == statusT::rankInsufficient );
        REQUIRE( allFinite( result.residuals.col( 0 ) ) );
        REQUIRE( allNan( result.residuals.col( 1 ) ) );
    }

    SECTION( "wide exact rank deficiency" )
    {
        pcaT::matrixT predictors( 2, 3 );
        predictors << 1, 2, 3, 2, 4, 6;
        pcaT::vectorT target( 2 );
        target << 1, -1;

        resultT result;
        mx::improc::P4PCA::calculate( result, predictors, target, { 1, 2 }, 1e-10, workspace );
        REQUIRE( result.numericalRank == 1 );
        REQUIRE( result.modeStatus[0] == statusT::rankSupported );
        REQUIRE( result.modeStatus[1] == statusT::rankInsufficient );
        REQUIRE( allFinite( result.residuals.col( 0 ) ) );
        REQUIRE( allNan( result.residuals.col( 1 ) ) );
    }

    SECTION( "wide exact rank deficiency with zero tolerance" )
    {
        pcaT::matrixT predictors( 2, 3 );
        predictors << 1, 2, 3, 2, 4, 6;
        pcaT::vectorT target( 2 );
        target << 1, -1;

        resultT result;
        mx::improc::P4PCA::calculate( result, predictors, target, { 1, 2 }, 0, workspace );
        REQUIRE( result.numericalRank == 1 );
        REQUIRE( result.modeStatus[0] == statusT::rankSupported );
        REQUIRE( result.modeStatus[1] == statusT::rankInsufficient );
        REQUIRE( allFinite( result.residuals.col( 0 ) ) );
        REQUIRE( allNan( result.residuals.col( 1 ) ) );
    }

    SECTION( "near-zero eigenvalue" )
    {
        pcaT::matrixT predictors( 3, 2 );
        predictors << 1, 0, 0, 1e-4, 0, 0;
        pcaT::vectorT target( 3 );
        target << 2, 3, 4;

        resultT result;
        mx::improc::P4PCA::calculate( result, predictors, target, { 1, 2 }, 1e-6, workspace );
        REQUIRE( result.numericalRank == 1 );
        REQUIRE( allFinite( result.residuals.col( 0 ) ) );
        REQUIRE( allNan( result.residuals.col( 1 ) ) );
    }

    SECTION( "wide near-zero eigenvalue" )
    {
        pcaT::matrixT predictors( 2, 3 );
        predictors << 1, 0, 0, 0, 1e-4, 0;
        pcaT::vectorT target( 2 );
        target << 2, 3;

        resultT result;
        mx::improc::P4PCA::calculate( result, predictors, target, { 1, 2 }, 1e-6, workspace );
        REQUIRE( result.numericalRank == 1 );
        REQUIRE( allFinite( result.residuals.col( 0 ) ) );
        REQUIRE( allNan( result.residuals.col( 1 ) ) );
    }

    SECTION( "zero predictor matrix" )
    {
        pcaT::matrixT predictors( 3, 2 );
        predictors.setZero();
        pcaT::vectorT target( 3 );
        target << 2, 3, 4;

        resultT result;
        mx::improc::P4PCA::calculate( result, predictors, target, { 1, 2 }, 0, workspace );
        REQUIRE( result.numericalRank == 0 );
        REQUIRE( result.modeStatus == std::vector<statusT>{ statusT::rankInsufficient, statusT::rankInsufficient } );
        REQUIRE( allNan( result.residuals ) );
    }

    SECTION( "wide zero predictor matrix" )
    {
        pcaT::matrixT predictors( 2, 3 );
        predictors.setZero();
        pcaT::vectorT target( 2 );
        target << 2, 3;

        resultT result;
        mx::improc::P4PCA::calculate( result, predictors, target, { 1, 2 }, 0, workspace );
        REQUIRE( result.numericalRank == 0 );
        REQUIRE( result.modeStatus == std::vector<statusT>{ statusT::rankInsufficient, statusT::rankInsufficient } );
        REQUIRE( allNan( result.residuals ) );
    }

    SECTION( "strict threshold equality" )
    {
        pcaT::matrixT predictors( 2, 2 );
        predictors.matrix().setIdentity();
        pcaT::vectorT target( 2 );
        target << 2, 3;

        resultT result;
        solverReset reset( fakeSolverBehavior::valid );
        mx::improc::P4PCA::calculate( result, predictors, target, { 1, 2 }, 0.5, workspace );
        REQUIRE( result.numericalRank == 1 );
        REQUIRE( result.modeStatus[0] == statusT::rankSupported );
        REQUIRE( result.modeStatus[1] == statusT::rankInsufficient );
        REQUIRE( allNan( result.residuals.col( 1 ) ) );
    }
}

/// Verify P4PCA rejects empty or incompatible arrays and invalid retained-mode sequences.
/** \ingroup P4PCA_unit_tests */
TEST_CASE( "P4PCA rejects invalid shapes and mode counts", "[P4PCA][validation]" )
{
    resultT result;
    pcaT::workspaceT workspace;
    pcaT::vectorT target( 2 );
    target << 1, 2;
    pcaT::matrixT predictors( 2, 2 );
    predictors.matrix().setIdentity();

    pcaT::matrixT noSamples( 0, 2 );
    pcaT::vectorT noTarget( 0 );
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculate( result, noSamples, noTarget, { 1 }, 0, workspace ),
                       std::invalid_argument );

    pcaT::matrixT noPredictors( 2, 0 );
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculate( result, noPredictors, target, { 1 }, 0, workspace ),
                       std::invalid_argument );

    pcaT::vectorT wrongTarget( 3 );
    wrongTarget.setZero();
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculate( result, predictors, wrongTarget, { 1 }, 0, workspace ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculate( result, predictors, target, {}, 0, workspace ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculate( result, predictors, target, { 0 }, 0, workspace ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculate( result, predictors, target, { -1 }, 0, workspace ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculate( result, predictors, target, { 1, 1 }, 0, workspace ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculate( result, predictors, target, { 2, 1 }, 0, workspace ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculate( result, predictors, target, { 3 }, 0, workspace ),
                       std::invalid_argument );
}

/// Verify P4PCA rejects nonfinite inputs and normal-equation overflow under inherited fast-math flags.
/** \ingroup P4PCA_unit_tests */
TEST_CASE( "P4PCA rejects invalid numeric input", "[P4PCA][validation][finite]" )
{
    resultT result;
    pcaT::workspaceT workspace;
    pcaT::matrixT predictors( 2, 2 );
    predictors.matrix().setIdentity();
    pcaT::vectorT target( 2 );
    target << 1, 2;

    predictors( 0, 0 ) = std::numeric_limits<double>::quiet_NaN();
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculate( result, predictors, target, { 1 }, 0, workspace ),
                       std::invalid_argument );
    predictors.matrix().setIdentity();
    predictors( 0, 0 ) = std::numeric_limits<double>::infinity();
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculate( result, predictors, target, { 1 }, 0, workspace ),
                       std::invalid_argument );
    predictors.matrix().setIdentity();

    target( 0 ) = std::numeric_limits<double>::quiet_NaN();
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculate( result, predictors, target, { 1 }, 0, workspace ),
                       std::invalid_argument );
    target( 0 ) = std::numeric_limits<double>::infinity();
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculate( result, predictors, target, { 1 }, 0, workspace ),
                       std::invalid_argument );
    target << 1, 2;

    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculate( result, predictors, target, { 1 }, -1, workspace ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculate( result,
                                                     predictors,
                                                     target,
                                                     { 1 },
                                                     std::numeric_limits<double>::quiet_NaN(),
                                                     workspace ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculate( result,
                                                     predictors,
                                                     target,
                                                     { 1 },
                                                     std::numeric_limits<double>::infinity(),
                                                     workspace ),
                       std::invalid_argument );

    mx::improc::P4PCA::calculate( result, predictors, target, { 1, 2 }, 0, workspace );
    REQUIRE( result.numericalRank == 2 );
    mx::improc::P4PCA::calculate( result, predictors, target, { 1 }, 1, workspace );
    REQUIRE( result.numericalRank == 0 );
    REQUIRE( allNan( result.residuals ) );

    pcaT::matrixT overflowingCovariance( 1, 1 );
    overflowingCovariance( 0, 0 ) = std::numeric_limits<double>::max();
    pcaT::vectorT oneTarget( 1 );
    oneTarget << 1;
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculate( result, overflowingCovariance, oneTarget, { 1 }, 0, workspace ),
                       std::runtime_error );

    pcaT::matrixT overflowingCrossProduct( 2, 1 );
    overflowingCrossProduct.setOnes();
    pcaT::vectorT hugeTarget( 2 );
    hugeTarget.setConstant( std::numeric_limits<double>::max() );
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculate( result, overflowingCrossProduct, hugeTarget, { 1 }, 0, workspace ),
                       std::runtime_error );
}

/// Verify P4PCA fully replaces caller output and safely reuses one workspace by reference.
/** \ingroup P4PCA_unit_tests */
TEST_CASE( "P4PCA replaces output and reuses caller workspace", "[P4PCA][workspace]" )
{
    resultT result;
    result.residuals.resize( 7, 7 );
    result.residuals.setConstant( 42 );
    result.modeStatus.assign( 7, statusT::rankInsufficient );
    result.numericalRank = 99;

    pcaT::matrixT firstPredictors( 3, 2 );
    firstPredictors << 1, 0, 0, 2, 0, 0;
    pcaT::vectorT firstTarget( 3 );
    firstTarget << 3, 5, 7;
    pcaT::workspaceT workspace;
    mx::improc::P4PCA::calculate( result, firstPredictors, firstTarget, { 1, 2 }, 1e-12, workspace );

    const auto *work = workspace.work;
    const auto *integerWork = workspace.iWork;

    pcaT::matrixT secondPredictors( 4, 2 );
    secondPredictors << 1, 0, 0, 2, 1, 0, 0, 0;
    pcaT::vectorT secondTarget( 4 );
    secondTarget << 3, 5, -1, 7;
    mx::improc::P4PCA::calculate( result, secondPredictors, secondTarget, { 1 }, 1e-12, workspace );

    REQUIRE( result.residuals.rows() == 4 );
    REQUIRE( result.residuals.cols() == 1 );
    REQUIRE( result.modeStatus == std::vector<statusT>{ statusT::rankSupported } );
    REQUIRE( result.numericalRank == 2 );
    REQUIRE( workspace.work == work );
    REQUIRE( workspace.iWork == integerWork );

    pcaT::matrixT widePredictors( 3, 4 );
    widePredictors << 3, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 0;
    pcaT::vectorT wideTarget( 3 );
    wideTarget << 5, -1, 2;
    mx::improc::P4PCA::calculate( result, widePredictors, wideTarget, { 1, 3 }, 1e-12, workspace );

    REQUIRE( result.residuals.rows() == 3 );
    REQUIRE( result.residuals.cols() == 2 );
    REQUIRE( result.numericalRank == 3 );
    requireApprox( result.residuals.col( 0 ), svdResidual( widePredictors, wideTarget, 1 ) );
    requireApprox( result.residuals.col( 1 ), svdResidual( widePredictors, wideTarget, 3 ) );

    mx::improc::P4PCA::calculate( result, secondPredictors, secondTarget, { 1 }, 1e-12, workspace );
    REQUIRE( result.residuals.rows() == 4 );
    REQUIRE( result.residuals.cols() == 1 );
    requireApprox( result.residuals.col( 0 ), svdResidual( secondPredictors, secondTarget, 1 ) );
}

/// Verify P4PCA propagates solver failures and rejects malformed finite/nonfinite solver output.
/** \ingroup P4PCA_unit_tests */
TEST_CASE( "P4PCA propagates eigensolver and invalid solver output", "[P4PCA][solver][validation]" )
{
    pcaT::matrixT predictors( 2, 2 );
    predictors.matrix().setIdentity();
    pcaT::vectorT target( 2 );
    target.setOnes();
    resultT result;
    pcaT::workspaceT workspace;

    SECTION( "solver status" )
    {
        solverReset reset( fakeSolverBehavior::failure );
        REQUIRE_THROWS_WITH( mx::improc::P4PCA::calculate( result, predictors, target, { 1 }, 0, workspace ),
                             "P4PCA eigensolver failed with status 37" );
        REQUIRE( solverCalls == 1 );
    }

    SECTION( "invalid dimensions" )
    {
        solverReset reset( fakeSolverBehavior::invalidDimensions );
        REQUIRE_THROWS_AS( mx::improc::P4PCA::calculate( result, predictors, target, { 1 }, 0, workspace ),
                           std::runtime_error );
    }

    SECTION( "nonfinite eigenvalue" )
    {
        solverReset reset( fakeSolverBehavior::nonfiniteEigenvalue );
        REQUIRE_THROWS_AS( mx::improc::P4PCA::calculate( result, predictors, target, { 1, 2 }, 0, workspace ),
                           std::runtime_error );
    }

    SECTION( "invalid eigenvalue shape" )
    {
        solverReset reset( fakeSolverBehavior::invalidEigenvalueShape );
        REQUIRE_THROWS_AS( mx::improc::P4PCA::calculate( result, predictors, target, { 1, 2 }, 0, workspace ),
                           std::runtime_error );
    }

    SECTION( "nonfinite eigenvector" )
    {
        solverReset reset( fakeSolverBehavior::nonfiniteEigenvector );
        REQUIRE_THROWS_AS( mx::improc::P4PCA::calculate( result, predictors, target, { 1, 2 }, 0, workspace ),
                           std::runtime_error );
    }

    SECTION( "unsorted eigenvalues" )
    {
        solverReset reset( fakeSolverBehavior::unsortedEigenvalues );
        REQUIRE_THROWS_WITH( mx::improc::P4PCA::calculate( result, predictors, target, { 1, 2 }, 0, workspace ),
                             "P4PCA eigensolver returned unsorted eigenvalues" );
    }

    SECTION( "all nonpositive eigenvalues" )
    {
        solverReset reset( fakeSolverBehavior::nonpositiveEigenvalues );
        mx::improc::P4PCA::calculate( result, predictors, target, { 1, 2 }, 0, workspace );
        REQUIRE( result.numericalRank == 0 );
        REQUIRE( allNan( result.residuals ) );
        REQUIRE( result.modeStatus == std::vector<statusT>{ statusT::rankInsufficient, statusT::rankInsufficient } );
    }

    SECTION( "nonfinite projection" )
    {
        pcaT::matrixT tallPredictors( 3, 2 );
        tallPredictors << 1, 0, 0, 1, 0, 0;
        pcaT::vectorT tallTarget( 3 );
        tallTarget.setOnes();
        solverReset reset( fakeSolverBehavior::tinyEigenvalue );
        REQUIRE_THROWS_AS( mx::improc::P4PCA::calculate( result, tallPredictors, tallTarget, { 2 }, 0, workspace ),
                           std::runtime_error );
    }

    SECTION( "nonfinite residual" )
    {
        solverReset reset( fakeSolverBehavior::residualOverflow );
        target.setConstant( std::numeric_limits<double>::max() );
        REQUIRE_THROWS_AS( mx::improc::P4PCA::calculate( result, predictors, target, { 1 }, 0, workspace ),
                           std::runtime_error );
    }
}

/// Verify centered-fit P4PCA applies its coefficients to uncentered predictor data in both Gram branches.
/** This exercises mx::improc::P4PCA::calculateCentered() across both adaptive Gram-matrix branches. */
TEST_CASE( "P4PCA centered fit with uncentered application agrees with direct SVD",
           "[P4PCA][centered][uncentered][reference]" )
{
    pcaT::workspaceT workspace;

    SECTION( "K is less than T with non-contiguous modes" )
    {
        pcaT::matrixT predictors( 5, 3 );
        predictors << 13, -4, 3, 8, 0, 2, 11, -5, 4, 10, -6, -2, 8, -5, 3;
        pcaT::vectorT target( 5 );
        target << 9, 4, 8, -1, 5;

        resultT result;
        mx::improc::P4PCA::calculateCentered( result, predictors, target, { 1, 3 }, 1e-12, workspace );

        REQUIRE( result.numericalRank == 3 );
        REQUIRE( result.modeStatus == std::vector<statusT>{ statusT::rankSupported, statusT::rankSupported } );
        requireApprox( result.residuals.col( 0 ), centeredFitUncenteredResidual( predictors, target, 1 ) );
        requireApprox( result.residuals.col( 1 ), centeredFitUncenteredResidual( predictors, target, 3 ) );
    }

    SECTION( "K equals or exceeds T" )
    {
        pcaT::matrixT predictors( 4, 5 );
        predictors << 1, 2, 0, 4, -1, 3, -1, 2, 0, 5, -2, 4, 1, -3, 2, 0, 1, -4, 2, 3;
        pcaT::vectorT target( 4 );
        target << 11, 3, 8, -2;

        resultT result;
        mx::improc::P4PCA::calculateCentered( result, predictors, target, { 1, 2, 3 }, 0, workspace );

        REQUIRE( result.numericalRank == 3 );
        REQUIRE( result.modeStatus ==
                 std::vector<statusT>{ statusT::rankSupported, statusT::rankSupported, statusT::rankSupported } );
        requireApprox( result.residuals.col( 0 ), centeredFitUncenteredResidual( predictors, target, 1 ) );
        requireApprox( result.residuals.col( 1 ), centeredFitUncenteredResidual( predictors, target, 2 ) );
        requireApprox( result.residuals.col( 2 ), centeredFitUncenteredResidual( predictors, target, 3 ) );
    }

    SECTION( "T equals two" )
    {
        pcaT::matrixT predictors( 2, 3 );
        predictors << 1, 4, -2, 3, 1, 6;
        pcaT::vectorT target( 2 );
        target << 7, 11;

        resultT result;
        mx::improc::P4PCA::calculateCentered( result, predictors, target, { 1 }, 0, workspace );

        REQUIRE( result.numericalRank == 1 );
        REQUIRE( result.modeStatus == std::vector<statusT>{ statusT::rankSupported } );
        requireApprox( result.residuals.col( 0 ), centeredFitUncenteredResidual( predictors, target, 1 ) );
    }
}

/// Verify centered P4PCA uses min(K,T-1) eigenpairs and the corresponding smaller Gram matrix.
/** This directly exercises mx::improc::P4PCA::calculateCentered() structural-rank dispatch. */
TEST_CASE( "P4PCA centered regression enforces structural degrees of freedom", "[P4PCA][centered][rank][solver]" )
{
    resultT result;
    pcaT::workspaceT workspace;

    SECTION( "predictor limited" )
    {
        pcaT::matrixT predictors( 5, 2 );
        predictors << 1, 8, 2, 4, 4, 3, -1, 2, 3, -2;
        pcaT::vectorT target( 5 );
        target << 7, 1, 3, 8, -2;
        const pcaT::matrixT centered = centeredColumns( predictors );
        const pcaT::matrixT expected = ( centered.matrix().transpose() * centered.matrix() ).array();

        solverReset reset( fakeSolverBehavior::failure );
        REQUIRE_THROWS_WITH( mx::improc::P4PCA::calculateCentered( result, predictors, target, { 1, 2 }, 0, workspace ),
                             "P4PCA eigensolver failed with status 37" );
        REQUIRE( solverModeCount == 2 );
        requireApprox( solverGram, expected );
    }

    SECTION( "sample limited" )
    {
        pcaT::matrixT predictors( 3, 5 );
        predictors << 1, 2, 0, 4, -1, 3, -1, 2, 0, 5, -2, 4, 1, -3, 2;
        pcaT::vectorT target( 3 );
        target << 7, 1, 3;
        const pcaT::matrixT centered = centeredColumns( predictors );
        const pcaT::matrixT expected = ( centered.matrix() * centered.matrix().transpose() ).array();

        solverReset reset( fakeSolverBehavior::failure );
        REQUIRE_THROWS_WITH( mx::improc::P4PCA::calculateCentered( result, predictors, target, { 1, 2 }, 0, workspace ),
                             "P4PCA eigensolver failed with status 37" );
        REQUIRE( solverModeCount == 2 );
        REQUIRE( solverGram.rows() == 3 );
        requireApprox( solverGram, expected );
    }
}

/// Verify centered P4PCA ignores target DC but applies the post-preprocessing predictor DC.
/** This locks the centered-fit, uncentered-application contract of mx::improc::P4PCA::calculateCentered(). */
TEST_CASE( "P4PCA centered fit applies uncentered predictor baselines", "[P4PCA][centered][uncentered][mean]" )
{
    pcaT::matrixT predictors( 5, 3 );
    predictors << 13, -4, 3, 8, 0, 2, 11, -5, 4, 10, -6, -2, 8, -5, 3;
    pcaT::vectorT target( 5 );
    target << 9, 4, 8, -1, 5;
    pcaT::vectorT shiftedTarget = target + 37;

    resultT original;
    resultT shifted;
    pcaT::workspaceT workspace;
    mx::improc::P4PCA::calculateCentered( original, predictors, target, { 1, 2, 3 }, 1e-12, workspace );
    mx::improc::P4PCA::calculateCentered( shifted, predictors, shiftedTarget, { 1, 2, 3 }, 1e-12, workspace );

    REQUIRE( shifted.numericalRank == original.numericalRank );
    REQUIRE( shifted.modeStatus == original.modeStatus );
    for( Eigen::Index column = 0; column < original.residuals.cols(); ++column )
    {
        pcaT::vectorT offset = shifted.residuals.col( column ) - original.residuals.col( column );
        pcaT::vectorT expectedOffset( target.rows() );
        expectedOffset.setConstant( 37 );
        requireApprox( offset, expectedOffset );
        requireApprox( original.residuals.col( column ),
                       centeredFitUncenteredResidual( predictors, target, column + 1 ) );
        requireApprox( shifted.residuals.col( column ),
                       centeredFitUncenteredResidual( predictors, shiftedTarget, column + 1 ) );
    }

    pcaT::matrixT predictorShift = predictors;
    pcaT::vectorT baseline( predictors.cols() );
    baseline << 8, -3, 5;
    for( Eigen::Index row = 0; row < predictorShift.rows(); ++row )
    {
        predictorShift.row( row ) += baseline.transpose();
    }

    resultT shiftedPredictors;
    mx::improc::P4PCA::calculateCentered( shiftedPredictors, predictorShift, target, { 1, 2, 3 }, 1e-12, workspace );
    REQUIRE( shiftedPredictors.numericalRank == original.numericalRank );
    REQUIRE( shiftedPredictors.modeStatus == original.modeStatus );
    for( Eigen::Index column = 0; column < original.residuals.cols(); ++column )
    {
        const pcaT::vectorT coefficients = centeredSvdCoefficients( predictors, target, column + 1 );
        pcaT::vectorT expectedOffset( target.rows() );
        expectedOffset.setConstant( -baseline.matrix().dot( coefficients.matrix() ) );
        requireApprox( shiftedPredictors.residuals.col( column ) - original.residuals.col( column ), expectedOffset );
        requireApprox( shiftedPredictors.residuals.col( column ),
                       centeredFitUncenteredResidual( predictorShift, target, column + 1 ) );
    }

    pcaT::vectorT constantTarget( 5 );
    constantTarget.setConstant( 12.5 );
    mx::improc::P4PCA::calculateCentered( shifted, predictors, constantTarget, { 1, 3 }, 1e-12, workspace );
    requireApprox( shifted.residuals.col( 0 ), constantTarget );
    requireApprox( shifted.residuals.col( 1 ), constantTarget );

    constantTarget.setZero();
    mx::improc::P4PCA::calculateCentered( shifted, predictors, constantTarget, { 1, 3 }, 1e-12, workspace );
    requireApprox( shifted.residuals.col( 0 ), constantTarget );
    requireApprox( shifted.residuals.col( 1 ), constantTarget );
}

/// Verify centered P4PCA handles numerical rank deficiency without admitting its structural null direction.
/** This exercises mx::improc::P4PCA::calculateCentered() at zero and nonzero rank tolerances. */
TEST_CASE( "P4PCA centered regression reports only supported numerical rank", "[P4PCA][centered][rank]" )
{
    pcaT::workspaceT workspace;

    SECTION( "structural null is excluded at zero tolerance" )
    {
        pcaT::matrixT predictors( 3, 3 );
        predictors << 1, 0, 2, 0, 1, 4, -1, -1, 6;
        pcaT::vectorT target( 3 );
        target << 4, -2, 7;

        resultT result;
        mx::improc::P4PCA::calculateCentered( result, predictors, target, { 1, 2 }, 0, workspace );
        REQUIRE( result.numericalRank == 2 );
        REQUIRE( result.modeStatus == std::vector<statusT>{ statusT::rankSupported, statusT::rankSupported } );
    }

    SECTION( "exact centered rank deficiency" )
    {
        pcaT::matrixT predictors( 4, 3 );
        predictors << 1, 2, 8, 2, 4, 8, 3, 6, 8, 4, 8, 8;
        pcaT::vectorT target( 4 );
        target << 4, -2, 7, 1;

        resultT result;
        mx::improc::P4PCA::calculateCentered( result, predictors, target, { 1, 2, 3 }, 1e-10, workspace );
        REQUIRE( result.numericalRank == 1 );
        REQUIRE( result.modeStatus ==
                 std::vector<statusT>{ statusT::rankSupported, statusT::rankInsufficient, statusT::rankInsufficient } );
        REQUIRE( allFinite( result.residuals.col( 0 ) ) );
        REQUIRE( allNan( result.residuals.col( 1 ) ) );
        REQUIRE( allNan( result.residuals.col( 2 ) ) );
    }

    SECTION( "constant predictors have zero centered rank" )
    {
        pcaT::matrixT predictors( 3, 2 );
        predictors << 2, -4, 2, -4, 2, -4;
        pcaT::vectorT target( 3 );
        target << 4, -2, 7;

        resultT result;
        mx::improc::P4PCA::calculateCentered( result, predictors, target, { 1, 2 }, 0, workspace );
        REQUIRE( result.numericalRank == 0 );
        REQUIRE( result.modeStatus == std::vector<statusT>{ statusT::rankInsufficient, statusT::rankInsufficient } );
        REQUIRE( allNan( result.residuals ) );
    }
}

/// Verify centered P4PCA rejects invalid sample counts, structural mode counts, and numeric inputs.
/** This exercises validation performed by mx::improc::P4PCA::calculateCentered(). */
TEST_CASE( "P4PCA centered regression rejects invalid requests", "[P4PCA][centered][validation]" )
{
    resultT result;
    pcaT::workspaceT workspace;

    pcaT::matrixT oneSample( 1, 2 );
    oneSample << 1, 2;
    pcaT::vectorT oneTarget( 1 );
    oneTarget << 3;
    REQUIRE_THROWS_WITH( mx::improc::P4PCA::calculateCentered( result, oneSample, oneTarget, { 1 }, 0, workspace ),
                         "P4PCA centered regression requires at least two samples" );

    pcaT::matrixT sampleLimited( 3, 4 );
    sampleLimited.setRandom();
    pcaT::vectorT target( 3 );
    target << 1, 2, 3;
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculateCentered( result, sampleLimited, target, { 3 }, 0, workspace ),
                       std::invalid_argument );

    pcaT::matrixT predictorLimited( 4, 2 );
    predictorLimited.setRandom();
    pcaT::vectorT tallTarget( 4 );
    tallTarget << 1, 2, 3, 4;
    REQUIRE_THROWS_AS(
        mx::improc::P4PCA::calculateCentered( result, predictorLimited, tallTarget, { 3 }, 0, workspace ),
        std::invalid_argument );

    pcaT::vectorT wrongTarget( 4 );
    wrongTarget.setZero();
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculateCentered( result, sampleLimited, wrongTarget, { 1 }, 0, workspace ),
                       std::invalid_argument );

    sampleLimited( 0, 0 ) = std::numeric_limits<double>::quiet_NaN();
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculateCentered( result, sampleLimited, target, { 1 }, 0, workspace ),
                       std::invalid_argument );

    pcaT::matrixT overflowingCenter( 3, 1 );
    overflowingCenter << std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(),
        -std::numeric_limits<double>::max();
    REQUIRE_THROWS_WITH( mx::improc::P4PCA::calculateCentered( result, overflowingCenter, target, { 1 }, 0, workspace ),
                         "P4PCA temporal centering produced nonfinite values" );
}

/// Verify one caller workspace can alternate between legacy and centered P4PCA paths.
/** This exercises both mx::improc::P4PCA::calculate() and mx::improc::P4PCA::calculateCentered(). */
TEST_CASE( "P4PCA reuses one workspace across centered and uncentered paths", "[P4PCA][centered][workspace]" )
{
    pcaT::matrixT predictors( 4, 3 );
    predictors << 1, 0, 2, 3, -1, 4, -2, 4, 1, 0, 1, -3;
    pcaT::vectorT target( 4 );
    target << 7, 1, 5, -2;

    resultT result;
    pcaT::workspaceT workspace;
    mx::improc::P4PCA::calculate( result, predictors, target, { 1, 3 }, 1e-12, workspace );
    requireApprox( result.residuals.col( 0 ), svdResidual( predictors, target, 1 ) );
    requireApprox( result.residuals.col( 1 ), svdResidual( predictors, target, 3 ) );

    mx::improc::P4PCA::calculateCentered( result, predictors, target, { 1, 3 }, 1e-12, workspace );
    requireApprox( result.residuals.col( 0 ), centeredFitUncenteredResidual( predictors, target, 1 ) );
    requireApprox( result.residuals.col( 1 ), centeredFitUncenteredResidual( predictors, target, 3 ) );

    mx::improc::P4PCA::calculate( result, predictors, target, { 2 }, 1e-12, workspace );
    requireApprox( result.residuals.col( 0 ), svdResidual( predictors, target, 2 ) );
}

} // namespace P4PCA_test
} // namespace unitTest
