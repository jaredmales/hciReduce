/** \file P4PCA_test.cpp
 * \brief Tests the pure numerical predictor used by Pixel Prediction Post-Processing.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "src/common/P4PCA.hpp"
#include "src/common/ReductionTiming.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace unitTest
{
namespace P4PCA_test
{

/** \cond P4PCA_test_harness */
using pcaT = mx::improc::P4PCA;
using reductionTimingT = mx::improc::ReductionTiming;
using resultT = mx::improc::P4PCAResult;
using statusT = mx::improc::P4PCAModeStatus;
using timingT = mx::improc::P4PCATiming;

/// Convert target-specific retained rows into the compact deleted-row contract used by production P4.
mx::improc::P4TargetExclusions exclusionsFromRetained(
    const std::vector<std::vector<std::size_t>> &retained /**< [in] target-specific retained temporal rows */ )
{
    std::vector<std::vector<Eigen::Index>> deleted( retained.size() );
    for( std::size_t target = 0; target < retained.size(); ++target )
    {
        auto retainedRow = retained[target].begin();
        for( std::size_t row = 0; row < retained.size(); ++row )
        {
            if( retainedRow != retained[target].end() && *retainedRow == row )
            {
                ++retainedRow;
            }
            else
            {
                deleted[target].push_back( static_cast<Eigen::Index>( row ) );
            }
        }
    }
    return mx::improc::P4TargetExclusions::fromExplicit( static_cast<Eigen::Index>( retained.size() ), deleted );
}

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

/// Return direct thin-SVD predictor coefficients using the requested largest singular modes.
pcaT::vectorT svdCoefficients( const pcaT::matrixT &predictors, /**< [in] predictor matrix */
                               const pcaT::vectorT &target,     /**< [in] target time series */
                               int modeCount /**< [in] number of largest singular modes */ )
{
    Eigen::JacobiSVD<Eigen::MatrixXd> decomposition( predictors.matrix(), Eigen::ComputeThinU | Eigen::ComputeThinV );
    return ( decomposition.matrixV().leftCols( modeCount ) *
             decomposition.singularValues().head( modeCount ).cwiseInverse().asDiagonal() *
             decomposition.matrixU().leftCols( modeCount ).transpose() * target.matrix() )
        .array();
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

/// Verify held-out science and probe responses against independently assembled retained-row fits.
void requireHeldOutProbeMatchesIndependent(
    const resultT &result,                                 /**< [in] held-out science result to verify */
    const pcaT::matrixT &probeResiduals,                   /**< [in] target-major frozen-probe responses to verify */
    const pcaT::matrixT &predictors,                       /**< [in] full predictor matrix */
    const pcaT::vectorT &target,                           /**< [in] full regression target */
    const pcaT::matrixT &probePredictors,                  /**< [in] frozen predictor-probe responses */
    const pcaT::vectorT &probeTarget,                      /**< [in] direct frozen target-probe response */
    const std::vector<std::vector<std::size_t>> &retained, /**< [in] target-specific retained rows */
    const std::vector<int> &modes,                         /**< [in] retained-mode counts */
    double tolerance /**< [in] absolute comparison tolerance */ )
{
    REQUIRE( result.residuals.rows() == predictors.rows() );
    REQUIRE( result.residuals.cols() == static_cast<Eigen::Index>( modes.size() ) );
    REQUIRE( result.sampleValidity.rows() == predictors.rows() );
    REQUIRE( result.sampleValidity.cols() == static_cast<Eigen::Index>( modes.size() ) );
    REQUIRE( probeResiduals.rows() == probePredictors.rows() );
    REQUIRE( probeResiduals.cols() == predictors.rows() * static_cast<Eigen::Index>( modes.size() ) );
    REQUIRE( retained.size() == static_cast<std::size_t>( predictors.rows() ) );

    std::vector<bool> anySupported( modes.size(), false );
    for( Eigen::Index heldOut = 0; heldOut < predictors.rows(); ++heldOut )
    {
        const std::vector<std::size_t> &rows = retained[static_cast<std::size_t>( heldOut )];
        const int structuralRank =
            static_cast<int>( std::min<std::size_t>( rows.size(), static_cast<std::size_t>( predictors.cols() ) ) );
        const auto supportedEnd = std::upper_bound( modes.begin(), modes.end(), structuralRank );
        const std::vector<int> supportedModes( modes.begin(), supportedEnd );

        resultT independent;
        pcaT::matrixT coefficients;
        if( !supportedModes.empty() )
        {
            pcaT::matrixT trainingPredictors( static_cast<Eigen::Index>( rows.size() ), predictors.cols() );
            pcaT::vectorT trainingTarget( static_cast<Eigen::Index>( rows.size() ) );
            for( std::size_t training = 0; training < rows.size(); ++training )
            {
                const Eigen::Index source = static_cast<Eigen::Index>( rows[training] );
                trainingPredictors.row( static_cast<Eigen::Index>( training ) ) = predictors.row( source );
                trainingTarget( static_cast<Eigen::Index>( training ) ) = target( source );
            }

            pcaT::workspaceT workspace;
            mx::improc::P4PCA::calculate( independent,
                                          trainingPredictors,
                                          trainingTarget,
                                          supportedModes,
                                          1e-12,
                                          workspace,
                                          nullptr,
                                          &coefficients );
        }

        for( std::size_t mode = 0; mode < modes.size(); ++mode )
        {
            const bool independentSupported =
                mode < supportedModes.size() && independent.modeStatus[mode] == statusT::rankSupported;
            const Eigen::Index probeColumn =
                heldOut * static_cast<Eigen::Index>( modes.size() ) + static_cast<Eigen::Index>( mode );
            CAPTURE( heldOut, mode, predictors.rows(), predictors.cols() );
            REQUIRE( result.sampleSupported( heldOut, mode ) == independentSupported );
            if( independentSupported )
            {
                anySupported[mode] = true;
                const double expectedResidual =
                    target( heldOut ) - predictors.row( heldOut ).matrix().dot( coefficients.col( mode ).matrix() );
                const pcaT::vectorT expectedProbe =
                    ( probeTarget.matrix() - probePredictors.matrix() * coefficients.col( mode ).matrix() ).array();
                REQUIRE( result.residuals( heldOut, static_cast<Eigen::Index>( mode ) ) ==
                         Approx( expectedResidual ).margin( tolerance ) );
                requireApprox( probeResiduals.col( probeColumn ), expectedProbe, tolerance );
                REQUIRE( allFinite( probeResiduals.col( probeColumn ) ) );
            }
            else
            {
                REQUIRE( mx::math::isNan( result.residuals( heldOut, static_cast<Eigen::Index>( mode ) ) ) );
                REQUIRE( allNan( probeResiduals.col( probeColumn ) ) );
            }
        }
    }

    for( std::size_t mode = 0; mode < modes.size(); ++mode )
    {
        REQUIRE( result.modeStatus[mode] ==
                 ( anySupported[mode] ? statusT::rankSupported : statusT::rankInsufficient ) );
    }
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
    marginallyNonorthonormalEigenvectors,
    excessivelyNonorthonormalEigenvectors,
    nonorthonormalEigenvectors,
    residualOverflow
};

/// Outcome supplied by the controlled eigensolver.
fakeSolverBehavior solverBehavior{ fakeSolverBehavior::failure };

/// Number of invocations of the controlled eigensolver.
int solverCalls{ 0 };

/// Status returned by the controlled FP64 eigensolver's failure behavior.
MXLAPACK_INT solverFailureStatus{ 37 };

/// Whether controlled tests solve a diagonal oracle exactly after the first injected base result.
bool solverUsesDiagonalOracleAfterFirst{ false };

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

    if( solverUsesDiagonalOracleAfterFirst && solverCalls > 1 )
    {
        const int dimension = static_cast<int>( covariance.rows() );
        if( modeCount != dimension )
        {
            return 38;
        }
        eigenvectors.resize( dimension, dimension );
        eigenvectors.setZero();
        eigenvalues.resize( dimension, 1 );
        for( int index = 0; index < dimension; ++index )
        {
            eigenvectors( index, index ) = 1;
            eigenvalues( index ) = covariance( index, index );
        }
        return 0;
    }

    if( solverBehavior == fakeSolverBehavior::failure )
    {
        return solverFailureStatus;
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
    else if( solverBehavior == fakeSolverBehavior::marginallyNonorthonormalEigenvectors )
    {
        const double defect = 96.0 * std::numeric_limits<double>::epsilon() * predictorCount;
        eigenvectors.col( 0 ) *= std::sqrt( 1.0 + defect );
    }
    else if( solverBehavior == fakeSolverBehavior::excessivelyNonorthonormalEigenvectors )
    {
        const double defect = 160.0 * std::numeric_limits<double>::epsilon() * predictorCount;
        eigenvectors.col( 0 ) *= std::sqrt( 1.0 + defect );
    }
    else if( solverBehavior == fakeSolverBehavior::nonorthonormalEigenvectors )
    {
        eigenvectors.col( 0 ) *= 2.0;
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
    explicit solverReset(
        fakeSolverBehavior behavior, /**< [in] controlled first or repeated outcome */
        bool diagonalOracleAfterFirst = false /**< [in] whether later calls solve the diagonal test Gram exactly */ )
    {
        solverBehavior = behavior;
        solverCalls = 0;
        solverFailureStatus = 37;
        solverUsesDiagonalOracleAfterFirst = diagonalOracleAfterFirst;
        solverGram.resize( 0, 0 );
        solverModeCount = 0;
        mx::improc::detail::p4PCASetEigenSolverForTesting( &fakeEigenSolver );
    }

    /// Restore the production eigensolver.
    ~solverReset()
    {
        mx::improc::detail::p4PCAResetEigenSolverForTesting();
        solverUsesDiagonalOracleAfterFirst = false;
        solverFailureStatus = 37;
    }
};

#ifdef HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION

/// Number of invocations of the controlled experimental FP32 eigensolver.
int floatSolverCalls{ 0 };

/// Status returned by the controlled experimental FP32 eigensolver.
MXLAPACK_INT floatSolverStatus{ 37 };

/// Gram matrix received by the controlled experimental FP32 eigensolver.
mx::improc::detail::P4PCAFloatMatrixT floatSolverGram;

/// Eigenpair count received by the controlled experimental FP32 eigensolver.
int floatSolverModeCount{ 0 };

/// Return a configured status from the experimental FP32 eigensolver seam.
MXLAPACK_INT
fakeFloatEigenSolver( mx::improc::detail::P4PCAFloatMatrixT &eigenvectors, /**< [out] unused controlled eigenvectors */
                      mx::improc::detail::P4PCAFloatMatrixT &eigenvalues,  /**< [out] unused controlled eigenvalues */
                      mx::improc::detail::P4PCAFloatMatrixT &covariance,   /**< [in] covariance matrix */
                      int modeCount, /**< [in] requested largest eigenpair count */
                      mx::math::syevrMem<float> &workspace /**< [in,out] unused caller workspace */ )
{
    static_cast<void>( eigenvectors );
    static_cast<void>( eigenvalues );
    static_cast<void>( workspace );
    ++floatSolverCalls;
    floatSolverGram = covariance;
    floatSolverModeCount = modeCount;
    return floatSolverStatus;
}

/// Restore both native eigensolvers after an experimental controlled-failure scope.
struct experimentalSolverReset
{
    /// Install controlled native eigensolvers that return one exact status.
    explicit experimentalSolverReset( MXLAPACK_INT status /**< [in] status returned by both native solvers */ )
    {
        solverBehavior = fakeSolverBehavior::failure;
        solverCalls = 0;
        solverFailureStatus = status;
        solverGram.resize( 0, 0 );
        solverModeCount = 0;
        floatSolverCalls = 0;
        floatSolverStatus = status;
        floatSolverGram.resize( 0, 0 );
        floatSolverModeCount = 0;
        mx::improc::detail::p4PCASetEigenSolverForTesting( &fakeEigenSolver );
        mx::improc::detail::p4PCASetFloatEigenSolverForTesting( &fakeFloatEigenSolver );
    }

    /// Restore both production native eigensolvers.
    ~experimentalSolverReset()
    {
        mx::improc::detail::p4PCAResetEigenSolverForTesting();
        mx::improc::detail::p4PCAResetFloatEigenSolverForTesting();
        solverFailureStatus = 37;
        floatSolverStatus = 37;
    }
};

/// Compare held-out scalar-policy results while preserving NaN and sample-validity semantics.
void requireExperimentalHeldOutApprox(
    const resultT &actual,              /**< [in] scalar-policy held-out science result */
    const pcaT::matrixT &actualProbe,   /**< [in] scalar-policy target-major probe response */
    const resultT &expected,            /**< [in] all-double held-out science result */
    const pcaT::matrixT &expectedProbe, /**< [in] all-double target-major probe response */
    double tolerance /**< [in] absolute tolerance for supported values */ )
{
    REQUIRE( actual.residuals.rows() == expected.residuals.rows() );
    REQUIRE( actual.residuals.cols() == expected.residuals.cols() );
    REQUIRE( actual.sampleValidity.rows() == expected.sampleValidity.rows() );
    REQUIRE( actual.sampleValidity.cols() == expected.sampleValidity.cols() );
    REQUIRE( actualProbe.rows() == expectedProbe.rows() );
    REQUIRE( actualProbe.cols() == expectedProbe.cols() );
    REQUIRE( actual.modeStatus == expected.modeStatus );
    REQUIRE( actual.numericalRank == expected.numericalRank );
    REQUIRE( actual.baseRank == expected.baseRank );
    REQUIRE( actual.numericalRankCapped == expected.numericalRankCapped );

    for( Eigen::Index target = 0; target < expected.residuals.rows(); ++target )
    {
        for( Eigen::Index mode = 0; mode < expected.residuals.cols(); ++mode )
        {
            CAPTURE( target, mode );
            REQUIRE( actual.sampleValidity( target, mode ) == expected.sampleValidity( target, mode ) );
            const Eigen::Index probeColumn = target * expected.residuals.cols() + mode;
            if( expected.sampleValidity( target, mode ) )
            {
                REQUIRE( actual.residuals( target, mode ) ==
                         Approx( expected.residuals( target, mode ) ).margin( tolerance ) );
                requireApprox( actualProbe.col( probeColumn ), expectedProbe.col( probeColumn ), tolerance );
            }
            else
            {
                REQUIRE( mx::math::isNan( actual.residuals( target, mode ) ) );
                REQUIRE( mx::math::isNan( expected.residuals( target, mode ) ) );
                REQUIRE( allNan( actualProbe.col( probeColumn ) ) );
                REQUIRE( allNan( expectedProbe.col( probeColumn ) ) );
            }
        }
    }
}

#endif // HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION

/// Inject an allocation failure only while mxlib validates the supplied temporal factor.
void failFactorValidationAllocation(
    mx::math::detail::svdDeletionTestOperation operation /**< [in] mxlib operation being attempted */ )
{
    if( operation == mx::math::detail::svdDeletionTestOperation::validateFactor )
    {
        throw std::bad_alloc();
    }
}

/// Number of deterministic structured deletion failures supplied by the test hook.
int structuredDeletionSolveCalls{ 0 };

/// LAPACK status returned by the deterministic structured deletion test hook.
MXLAPACK_INT structuredDeletionSolveInfo{ 73 };

/// Return a deterministic LAPACK failure from the structured deletion eigensolver seam.
MXLAPACK_INT
failStructuredDeletionSolve( double *eigenvalues,      /**< [out] unused updated eigenvalues */
                             double *eigenvectors,     /**< [out] unused updated eigenvectors */
                             double *secularWorkspace, /**< [out] unused secular-equation workspace */
                             MXLAPACK_INT dimension,   /**< [in] unused secular-equation dimension */
                             MXLAPACK_INT outputRank,  /**< [in] unused requested output rank */
                             double updateMagnitude,   /**< [in] unused rank-one update magnitude */
                             double *diagonal,         /**< [in] unused diagonal spectrum */
                             double *deflationAdjustedUpdate /**< [in] unused deflation-adjusted update vector */ )
{
    ++structuredDeletionSolveCalls;
    static_cast<void>( eigenvalues );
    static_cast<void>( eigenvectors );
    static_cast<void>( secularWorkspace );
    static_cast<void>( dimension );
    static_cast<void>( outputRank );
    static_cast<void>( updateMagnitude );
    static_cast<void>( diagonal );
    static_cast<void>( deflationAdjustedUpdate );
    return structuredDeletionSolveInfo;
}

/// Restore mxlib's process-wide SVD-deletion hooks after a controlled failure test.
struct svdDeletionHookReset
{
    /// Clear preexisting hooks before installing one controlled outcome.
    svdDeletionHookReset()
    {
        mx::math::detail::svdDeletionHooks<double>() = {};
    }

    /// Restore ordinary mxlib execution.
    ~svdDeletionHookReset()
    {
        mx::math::detail::svdDeletionHooks<double>() = {};
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

/// Verify P4PCA optionally returns the exact predictor coefficients in both adaptive Gram branches.
/** This exercises mx::improc::P4PCA::calculate() with and without its optional coefficient output. */
TEST_CASE( "P4PCA returns optional predictor coefficients", "[P4PCA][regression][coefficients]" )
{
    SECTION( "predictor Gram" )
    {
        pcaT::matrixT predictors( 5, 3 );
        predictors << 1, 0, 2, 0, 2, -1, 3, 1, 0, -2, 0, 1, 0.5, -1, 4;
        pcaT::vectorT target( 5 );
        target << 3, 5, -1, 7, 2;

        resultT result;
        resultT resultWithoutCoefficients;
        pcaT::matrixT coefficients;
        pcaT::workspaceT workspace;
        pcaT::workspaceT workspaceWithoutCoefficients;
        pcaT::calculate( result, predictors, target, { 1, 3 }, 1e-12, workspace, nullptr, &coefficients );
        pcaT::calculate( resultWithoutCoefficients, predictors, target, { 1, 3 }, 1e-12, workspaceWithoutCoefficients );

        REQUIRE( coefficients.rows() == predictors.cols() );
        REQUIRE( coefficients.cols() == 2 );
        REQUIRE( ( result.residuals == resultWithoutCoefficients.residuals ).all() );
        requireApprox( coefficients.col( 0 ), svdCoefficients( predictors, target, 1 ) );
        requireApprox( coefficients.col( 1 ), svdCoefficients( predictors, target, 3 ) );
        requireApprox( ( target.matrix() - predictors.matrix() * coefficients.col( 0 ).matrix() ).array(),
                       result.residuals.col( 0 ) );
        requireApprox( ( target.matrix() - predictors.matrix() * coefficients.col( 1 ).matrix() ).array(),
                       result.residuals.col( 1 ) );
    }

    SECTION( "temporal Gram" )
    {
        pcaT::matrixT predictors( 3, 5 );
        predictors << 3, 0.2, 0, 1, -0.5, 0.1, 2, 0.4, -0.3, 0.8, 0, -0.2, 1, 0.7, 2;
        pcaT::vectorT target( 3 );
        target << 5, -1, 2;

        resultT result;
        resultT resultWithoutCoefficients;
        pcaT::matrixT coefficients;
        pcaT::workspaceT workspace;
        pcaT::workspaceT workspaceWithoutCoefficients;
        pcaT::calculate( result, predictors, target, { 1, 2, 3 }, 1e-12, workspace, nullptr, &coefficients );
        pcaT::calculate( resultWithoutCoefficients,
                         predictors,
                         target,
                         { 1, 2, 3 },
                         1e-12,
                         workspaceWithoutCoefficients );

        REQUIRE( coefficients.rows() == predictors.cols() );
        REQUIRE( coefficients.cols() == 3 );
        REQUIRE( ( result.residuals == resultWithoutCoefficients.residuals ).all() );
        for( Eigen::Index output = 0; output < coefficients.cols(); ++output )
        {
            requireApprox( coefficients.col( output ), svdCoefficients( predictors, target, output + 1 ) );
            requireApprox( ( target.matrix() - predictors.matrix() * coefficients.col( output ).matrix() ).array(),
                           result.residuals.col( output ) );
        }
    }

    SECTION( "disabled output remains untouched" )
    {
        pcaT::matrixT predictors( 3, 2 );
        predictors << 1, 0, 0, 2, 0, 0;
        pcaT::vectorT target( 3 );
        target << 3, 5, 7;
        pcaT::matrixT coefficients( 1, 1 );
        coefficients( 0, 0 ) = 42;

        resultT result;
        pcaT::workspaceT workspace;
        pcaT::calculate( result, predictors, target, { 1, 2 }, 1e-12, workspace );

        REQUIRE( coefficients.rows() == 1 );
        REQUIRE( coefficients.cols() == 1 );
        REQUIRE( coefficients( 0, 0 ) == 42 );
    }
}

/// Verify P4PCA reports finite nonnegative timing components without changing its numerical result.
/** \ingroup P4PCA_unit_tests */
TEST_CASE( "P4PCA reports numerical-stage timing", "[P4PCA][timing]" )
{
    pcaT::matrixT predictors( 3, 2 );
    predictors << 1, 0, 0, 2, 0, 0;
    pcaT::vectorT target( 3 );
    target << 3, 5, 7;

    resultT result;
    pcaT::workspaceT workspace;
    timingT timing;
    timing.gramWorkerSeconds = -1;
    timing.eigensolveWorkerSeconds = -1;
    timing.baseFactorWorkerSeconds = -1;
    timing.deletionWorkerSeconds = -1;
    timing.explicitFallbackWorkerSeconds = -1;
    timing.projectionWorkerSeconds = -1;
    pcaT::calculate( result, predictors, target, { 1, 2 }, 1e-12, workspace, &timing );

    REQUIRE( result.numericalRank == 2 );
    REQUIRE( std::isfinite( timing.gramWorkerSeconds ) );
    REQUIRE( std::isfinite( timing.eigensolveWorkerSeconds ) );
    REQUIRE( std::isfinite( timing.baseFactorWorkerSeconds ) );
    REQUIRE( std::isfinite( timing.deletionWorkerSeconds ) );
    REQUIRE( std::isfinite( timing.explicitFallbackWorkerSeconds ) );
    REQUIRE( std::isfinite( timing.projectionWorkerSeconds ) );
    REQUIRE( timing.gramWorkerSeconds >= 0 );
    REQUIRE( timing.eigensolveWorkerSeconds >= 0 );
    REQUIRE( timing.baseFactorWorkerSeconds == 0 );
    REQUIRE( timing.deletionWorkerSeconds == 0 );
    REQUIRE( timing.explicitFallbackWorkerSeconds == 0 );
    REQUIRE( timing.projectionWorkerSeconds >= 0 );
}

/// Verify ReductionTiming clears both elapsed and aggregate-worker measurements between reductions.
/** \ingroup P4PCA_unit_tests */
TEST_CASE( "ReductionTiming resets all measurements", "[P4PCA][timing]" )
{
    reductionTimingT timing;
    timing.geometryElapsedSeconds = 1;
    timing.regressionElapsedSeconds = 2;
    timing.samplingWorkerSeconds = 3;
    timing.sameImageSamplingWorkerSeconds = 4;
    timing.temporalSamplingWorkerSeconds = 5;
    timing.gramWorkerSeconds = 6;
    timing.eigensolveWorkerSeconds = 7;
    timing.baseFactorWorkerSeconds = 8;
    timing.deletionWorkerSeconds = 9;
    timing.explicitFallbackWorkerSeconds = 10;
    timing.modeWorkerSeconds = 11;
    timing.projectionWorkerSeconds = 12;
    timing.psfWorkerSeconds = 13;
    timing.psfReconstructionElapsedSeconds = 14;

    timing.reset();

    REQUIRE( timing.geometryElapsedSeconds == 0 );
    REQUIRE( timing.regressionElapsedSeconds == 0 );
    REQUIRE( timing.samplingWorkerSeconds == 0 );
    REQUIRE( timing.sameImageSamplingWorkerSeconds == 0 );
    REQUIRE( timing.temporalSamplingWorkerSeconds == 0 );
    REQUIRE( timing.gramWorkerSeconds == 0 );
    REQUIRE( timing.eigensolveWorkerSeconds == 0 );
    REQUIRE( timing.baseFactorWorkerSeconds == 0 );
    REQUIRE( timing.deletionWorkerSeconds == 0 );
    REQUIRE( timing.explicitFallbackWorkerSeconds == 0 );
    REQUIRE( timing.modeWorkerSeconds == 0 );
    REQUIRE( timing.projectionWorkerSeconds == 0 );
    REQUIRE( timing.psfWorkerSeconds == 0 );
    REQUIRE( timing.psfReconstructionElapsedSeconds == 0 );
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
        pcaT::matrixT coefficients;
        solverReset reset( fakeSolverBehavior::failure );
        REQUIRE_THROWS_WITH(
            mx::improc::P4PCA::calculate( result, predictors, target, { 1 }, 0, workspace, nullptr, &coefficients ),
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

/// Verify exact held-out P4 refits match independently assembled training problems and retain sample-level validity.
/** This directly exercises mx::improc::P4PCA::calculateHeldOut() and
 * mx::improc::P4PCAResult::sampleSupported().
 */
TEST_CASE( "P4PCA exact held-out regression", "[P4PCA][held-out]" )
{
    pcaT::matrixT predictors( 4, 2 );
    predictors << 1, 0, 0, 1, 1, 1, 2, -1;
    pcaT::vectorT target( 4 );
    target << 2, -1, 1, 5;
    const std::vector<std::vector<std::size_t>> references{ { 1, 2, 3 }, { 0, 2, 3 }, { 0, 1, 3 }, { 0, 1, 2 } };
    const mx::improc::P4TargetExclusions exclusions = exclusionsFromRetained( references );
    resultT heldOut;
    pcaT::workspaceT workspace;
    mx::improc::P4PCATiming timing;
    mx::improc::P4PCA::calculateHeldOut( heldOut, predictors, target, exclusions, { 1, 2 }, 1e-12, workspace, &timing );

    REQUIRE( heldOut.residuals.rows() == 4 );
    REQUIRE( heldOut.residuals.cols() == 2 );
    REQUIRE( heldOut.sampleValidity.rows() == 4 );
    REQUIRE( heldOut.sampleValidity.cols() == 2 );
    REQUIRE( timing.gramWorkerSeconds >= 0 );
    REQUIRE( timing.eigensolveWorkerSeconds >= 0 );
    REQUIRE( timing.projectionWorkerSeconds >= 0 );
    for( Eigen::Index heldOutRow = 0; heldOutRow < 4; ++heldOutRow )
    {
        pcaT::matrixT trainingPredictors( 3, 2 );
        pcaT::vectorT trainingTarget( 3 );
        for( Eigen::Index training = 0; training < 3; ++training )
        {
            const Eigen::Index source = static_cast<Eigen::Index>( references[heldOutRow][training] );
            trainingPredictors.row( training ) = predictors.row( source );
            trainingTarget( training ) = target( source );
        }
        resultT independent;
        pcaT::matrixT coefficients;
        pcaT::workspaceT independentWorkspace;
        mx::improc::P4PCA::calculate( independent,
                                      trainingPredictors,
                                      trainingTarget,
                                      { 1, 2 },
                                      1e-12,
                                      independentWorkspace,
                                      nullptr,
                                      &coefficients );
        for( std::size_t mode = 0; mode < 2; ++mode )
        {
            REQUIRE( heldOut.sampleSupported( heldOutRow, mode ) );
            const double expected =
                target( heldOutRow ) - predictors.row( heldOutRow ).matrix().dot( coefficients.col( mode ).matrix() );
            REQUIRE( heldOut.residuals( heldOutRow, static_cast<Eigen::Index>( mode ) ) == Approx( expected ) );
        }
    }

    const std::vector<std::vector<std::size_t>> limitedReferences{ { 1 }, { 0, 2, 3 }, { 0, 1, 3 }, { 2 } };
    mx::improc::P4PCA::calculateHeldOut( heldOut,
                                         predictors,
                                         target,
                                         exclusionsFromRetained( limitedReferences ),
                                         { 1, 2 },
                                         1e-12,
                                         workspace );
    REQUIRE( heldOut.sampleSupported( 0, 0 ) );
    REQUIRE_FALSE( heldOut.sampleSupported( 0, 1 ) );
    REQUIRE( heldOut.sampleSupported( 1, 1 ) );
    REQUIRE_FALSE( heldOut.sampleSupported( 3, 1 ) );
    REQUIRE_THROWS_AS( heldOut.sampleSupported( 4, 0 ), std::out_of_range );
}

/// Verify the reused temporal Gram path agrees with independent explicit refits for every held-out target.
/** This directly exercises mx::improc::P4PCA::calculateHeldOut() when T is no greater than K and compares its output
 * with mx::improc::P4PCA::calculate() on explicitly assembled training matrices.
 */
TEST_CASE( "P4PCA temporal Gram held-out regression matches explicit refits", "[P4PCA][held-out][temporal-gram]" )
{
    pcaT::matrixT predictors( 4, 6 );
    predictors << 1, 0, 0, 0, 1, 2, 0, 1, 0, 0, 2, -1, 0, 0, 1, 0, -1, 1, 0, 0, 0, 1, 1, 0.5;
    pcaT::vectorT target( 4 );
    target << 2, -1, 3, 0.5;
    const std::vector<std::vector<std::size_t>> references{ { 1, 2, 3 }, { 2, 3 }, { 0, 1, 3 }, { 0, 1 } };
    const mx::improc::P4TargetExclusions exclusions = exclusionsFromRetained( references );
    const std::vector<int> modes{ 1, 2, 3 };

    resultT heldOut;
    pcaT::workspaceT workspace;
    mx::improc::P4PCA::calculateHeldOut( heldOut, predictors, target, exclusions, modes, 1e-12, workspace );

    for( Eigen::Index heldOutRow = 0; heldOutRow < predictors.rows(); ++heldOutRow )
    {
        const std::vector<std::size_t> &rows = references[static_cast<std::size_t>( heldOutRow )];
        pcaT::matrixT trainingPredictors( static_cast<Eigen::Index>( rows.size() ), predictors.cols() );
        pcaT::vectorT trainingTarget( static_cast<Eigen::Index>( rows.size() ) );
        for( std::size_t training = 0; training < rows.size(); ++training )
        {
            const Eigen::Index source = static_cast<Eigen::Index>( rows[training] );
            trainingPredictors.row( static_cast<Eigen::Index>( training ) ) = predictors.row( source );
            trainingTarget( static_cast<Eigen::Index>( training ) ) = target( source );
        }

        const auto supportedEnd = std::upper_bound( modes.begin(), modes.end(), static_cast<int>( rows.size() ) );
        const std::vector<int> supportedModes( modes.begin(), supportedEnd );
        resultT independent;
        pcaT::matrixT coefficients;
        pcaT::workspaceT independentWorkspace;
        mx::improc::P4PCA::calculate( independent,
                                      trainingPredictors,
                                      trainingTarget,
                                      supportedModes,
                                      1e-12,
                                      independentWorkspace,
                                      nullptr,
                                      &coefficients );
        for( std::size_t mode = 0; mode < modes.size(); ++mode )
        {
            if( mode >= supportedModes.size() )
            {
                REQUIRE_FALSE( heldOut.sampleSupported( heldOutRow, mode ) );
                continue;
            }
            REQUIRE( heldOut.sampleSupported( heldOutRow, mode ) );
            const double expected =
                target( heldOutRow ) - predictors.row( heldOutRow ).matrix().dot( coefficients.col( mode ).matrix() );
            REQUIRE( heldOut.residuals( heldOutRow, static_cast<Eigen::Index>( mode ) ) ==
                     Approx( expected ).margin( 1e-10 ) );
        }
    }
}

/// Verify target-specific frozen-probe responses reproduce independent retained-row coefficient fits.
/** This directly exercises mx::improc::P4PCA::calculateHeldOutProbe() in both adaptive Gram orientations with
 * target-specific multi-row exclusions and rank-insufficient target/mode combinations.
 *
 * \ingroup P4PCA_unit_tests
 */
TEST_CASE( "P4PCA held-out frozen probe matches independent retained-row fits", "[P4PCA][held-out][probe][reference]" )
{
    const auto compare = []( const pcaT::matrixT &predictors,
                             const pcaT::vectorT &target,
                             const pcaT::matrixT &probePredictors,
                             const pcaT::vectorT &probeTarget,
                             const std::vector<std::vector<std::size_t>> &retained,
                             const std::vector<int> &modes )
    {
        resultT result;
        pcaT::matrixT probeResiduals;
        pcaT::workspaceT workspace;
        const mx::improc::P4TargetExclusions exclusions = exclusionsFromRetained( retained );
        mx::improc::P4PCA::calculateHeldOutProbe( result,
                                                  probeResiduals,
                                                  predictors,
                                                  target,
                                                  probePredictors,
                                                  probeTarget,
                                                  exclusions,
                                                  modes,
                                                  1e-12,
                                                  workspace );

        requireHeldOutProbeMatchesIndependent( result,
                                               probeResiduals,
                                               predictors,
                                               target,
                                               probePredictors,
                                               probeTarget,
                                               retained,
                                               modes,
                                               2e-10 );
    };

    SECTION( "T is no greater than K" )
    {
        pcaT::matrixT predictors( 4, 6 );
        predictors << 4, 0.2, 0.1, 0, 0.3, -0.2, 0.1, 3, 0.4, 0.2, 0, 0.1, 0.2, -0.1, 2, 0.3, 0.1, 0.2, 0.1, 0.2, 0.3,
            1.5, -0.2, 0.1;
        pcaT::vectorT target( 4 );
        target << 2, -1, 3, 0.5;
        pcaT::matrixT probePredictors( 3, 6 );
        probePredictors << 0.3, -0.2, 0.7, 0.1, -0.4, 0.5, -0.1, 0.8, 0.2, -0.3, 0.6, 0.4, 0.5, 0.1, -0.2, 0.9, 0.3,
            -0.6;
        pcaT::vectorT probeTarget( 3 );
        probeTarget << 1.2, -0.7, 0.3;
        const std::vector<std::vector<std::size_t>> retained{ { 1, 2, 3 }, { 0, 2 }, { 3 }, {} };

        compare( predictors, target, probePredictors, probeTarget, retained, { 1, 2, 3, 4 } );
    }

    SECTION( "K is less than T" )
    {
        pcaT::matrixT predictors( 6, 3 );
        predictors << 4, 0.2, 0.1, 0.1, 3, 0.4, 0.2, -0.1, 2, 0.1, 0.2, 0.3, 1, -0.4, 0.2, -0.3, 0.5, 1.1;
        pcaT::vectorT target( 6 );
        target << 2, -1, 3, 0.5, 4, -2;
        pcaT::matrixT probePredictors( 4, 3 );
        probePredictors << 0.3, -0.2, 0.7, -0.1, 0.8, 0.2, 0.5, 0.1, -0.2, -0.4, 0.6, 0.9;
        pcaT::vectorT probeTarget( 4 );
        probeTarget << 1.2, -0.7, 0.3, 0.8;
        const std::vector<std::vector<std::size_t>> retained{ { 1, 2, 3, 4, 5 },
                                                              { 0, 2, 3 },
                                                              { 0, 1 },
                                                              { 4 },
                                                              {},
                                                              { 0, 1, 2 } };

        compare( predictors, target, probePredictors, probeTarget, retained, { 1, 2, 3 } );
    }
}

/// Verify explicit, dense-downdated, and secular-downdated frozen-probe responses agree target by target.
/** This directly compares mx::improc::P4PCA::calculateHeldOutProbe() and
 * mx::improc::P4PCA::calculateHeldOutProbeDowndated() with both supported row-deletion backends against independent
 * retained-row coefficient fits. It covers both Gram orientations, non-contiguous modes, and structurally
 * unsupported response columns.
 *
 * \ingroup P4PCA_unit_tests
 */
TEST_CASE( "P4PCA downdated held-out frozen probe matches explicit retained-row fits",
           "[P4PCA][held-out][probe][downdate][structured][exact]" )
{
    const auto compare = []( const pcaT::matrixT &predictors,
                             const pcaT::vectorT &target,
                             const pcaT::matrixT &probePredictors,
                             const pcaT::vectorT &probeTarget,
                             const std::vector<int> &modes )
    {
        std::vector<std::vector<Eigen::Index>> deleted( static_cast<std::size_t>( predictors.rows() ) );
        std::vector<std::vector<std::size_t>> retained( static_cast<std::size_t>( predictors.rows() ) );
        for( Eigen::Index heldOut = 0; heldOut < predictors.rows(); ++heldOut )
        {
            deleted[static_cast<std::size_t>( heldOut )].push_back( heldOut );
            for( Eigen::Index row = 0; row < predictors.rows(); ++row )
            {
                if( row != heldOut )
                {
                    retained[static_cast<std::size_t>( heldOut )].push_back( static_cast<std::size_t>( row ) );
                }
            }
        }
        const mx::improc::P4TargetExclusions exclusions =
            mx::improc::P4TargetExclusions::fromExplicit( predictors.rows(), deleted );

        resultT explicitResult;
        resultT denseResult;
        resultT structuredResult;
        pcaT::matrixT explicitProbe;
        pcaT::matrixT denseProbe;
        pcaT::matrixT structuredProbe;
        pcaT::workspaceT explicitWorkspace;
        pcaT::workspaceT denseWorkspace;
        pcaT::workspaceT structuredWorkspace;
        mx::improc::P4PCADowndateWorkspace denseDowndateWorkspace;
        mx::improc::P4PCADowndateWorkspace structuredDowndateWorkspace;

        mx::improc::P4PCA::calculateHeldOutProbe( explicitResult,
                                                  explicitProbe,
                                                  predictors,
                                                  target,
                                                  probePredictors,
                                                  probeTarget,
                                                  exclusions,
                                                  modes,
                                                  1e-12,
                                                  explicitWorkspace );
        mx::improc::P4PCA::calculateHeldOutProbeDowndated( denseResult,
                                                           denseProbe,
                                                           predictors,
                                                           target,
                                                           probePredictors,
                                                           probeTarget,
                                                           exclusions,
                                                           modes,
                                                           1e-12,
                                                           mx::math::svdDeletionBackend::leadingCovariance,
                                                           denseWorkspace,
                                                           denseDowndateWorkspace );
        mx::improc::P4PCA::calculateHeldOutProbeDowndated( structuredResult,
                                                           structuredProbe,
                                                           predictors,
                                                           target,
                                                           probePredictors,
                                                           probeTarget,
                                                           exclusions,
                                                           modes,
                                                           1e-12,
                                                           mx::math::svdDeletionBackend::rankOneSecular,
                                                           structuredWorkspace,
                                                           structuredDowndateWorkspace );

        REQUIRE( denseResult.explicitFallbackCount == 0 );
        REQUIRE( denseResult.explicitFallbackReason == mx::improc::P4PCAFallbackReason::none );
        REQUIRE( structuredResult.explicitFallbackCount == 0 );
        REQUIRE( structuredResult.explicitFallbackReason == mx::improc::P4PCAFallbackReason::none );
        REQUIRE( denseResult.modeStatus == explicitResult.modeStatus );
        REQUIRE( structuredResult.modeStatus == explicitResult.modeStatus );
        REQUIRE( ( denseResult.sampleValidity == explicitResult.sampleValidity ).all() );
        REQUIRE( ( structuredResult.sampleValidity == explicitResult.sampleValidity ).all() );

        requireHeldOutProbeMatchesIndependent( explicitResult,
                                               explicitProbe,
                                               predictors,
                                               target,
                                               probePredictors,
                                               probeTarget,
                                               retained,
                                               modes,
                                               2e-10 );
        requireHeldOutProbeMatchesIndependent( denseResult,
                                               denseProbe,
                                               predictors,
                                               target,
                                               probePredictors,
                                               probeTarget,
                                               retained,
                                               modes,
                                               2e-8 );
        requireHeldOutProbeMatchesIndependent( structuredResult,
                                               structuredProbe,
                                               predictors,
                                               target,
                                               probePredictors,
                                               probeTarget,
                                               retained,
                                               modes,
                                               2e-8 );
    };

    SECTION( "T is no greater than K" )
    {
        pcaT::matrixT predictors( 4, 6 );
        predictors << 4, 0.2, 0.1, 0, 0.3, -0.2, 0.1, 3, 0.4, 0.2, 0, 0.1, 0.2, -0.1, 2, 0.3, 0.1, 0.2, 0.1, 0.2, 0.3,
            1.5, -0.2, 0.1;
        pcaT::vectorT target( 4 );
        target << 2, -1, 3, 0.5;
        pcaT::matrixT probePredictors( 3, 6 );
        probePredictors << 0.3, -0.2, 0.7, 0.1, -0.4, 0.5, -0.1, 0.8, 0.2, -0.3, 0.6, 0.4, 0.5, 0.1, -0.2, 0.9, 0.3,
            -0.6;
        pcaT::vectorT probeTarget( 3 );
        probeTarget << 1.2, -0.7, 0.3;

        compare( predictors, target, probePredictors, probeTarget, { 1, 2, 3, 4 } );
    }

    SECTION( "K is less than T" )
    {
        pcaT::matrixT predictors( 6, 4 );
        predictors << 4, 0.2, 0.1, 0, 0.1, 3, 0.4, 0.2, 0.2, -0.1, 2, 0.3, 0.1, 0.2, 0.3, 1.5, 1, -0.4, 0.2, 0.7, -0.3,
            0.5, 1.1, -0.2;
        pcaT::vectorT target( 6 );
        target << 2, -1, 3, 0.5, 4, -2;
        pcaT::matrixT probePredictors( 3, 4 );
        probePredictors << 0.3, -0.2, 0.7, 0.1, -0.1, 0.8, 0.2, -0.3, 0.5, 0.1, -0.2, 0.9;
        pcaT::vectorT probeTarget( 3 );
        probeTarget << 1.2, -0.7, 0.3;

        compare( predictors, target, probePredictors, probeTarget, { 1, 2, 4 } );
    }
}

/// Verify full-rank factor deletion reproduces explicit held-out P4 regression across matrix orientations.
/** This directly compares mx::improc::P4PCA::calculateHeldOutDowndated() with
 * mx::improc::P4PCA::calculateHeldOut() for wide, square, and tall predictor matrices, variable deletion sets,
 * workspace reuse, structural rank loss, and the all-rows-deleted edge case.
 * \ingroup P4PCA_unit_tests
 */
TEST_CASE( "P4PCA exact factor downdate matches explicit held-out refits", "[P4PCA][held-out][downdate][exact]" )
{
    pcaT::workspaceT explicitWorkspace;
    pcaT::workspaceT baseWorkspace;
    mx::improc::P4PCADowndateWorkspace downdateWorkspace;

    const auto compare = [&]( const pcaT::matrixT &predictors,
                              const pcaT::vectorT &target,
                              const std::vector<std::vector<Eigen::Index>> &deleted,
                              const std::vector<int> &modes,
                              int expectedBaseRank )
    {
        const auto exclusions = mx::improc::P4TargetExclusions::fromExplicit( predictors.rows(), deleted );
        resultT explicitResult;
        resultT downdatedResult;
        mx::improc::P4PCA::calculateHeldOut( explicitResult,
                                             predictors,
                                             target,
                                             exclusions,
                                             modes,
                                             1e-12,
                                             explicitWorkspace );
        mx::improc::P4PCA::calculateHeldOutDowndated( downdatedResult,
                                                      predictors,
                                                      target,
                                                      exclusions,
                                                      modes,
                                                      1e-12,
                                                      baseWorkspace,
                                                      downdateWorkspace );

        REQUIRE( downdatedResult.baseRank == expectedBaseRank );
        REQUIRE_FALSE( downdatedResult.numericalRankCapped );
        REQUIRE( downdatedResult.explicitFallbackCount == 0 );
        REQUIRE( downdatedResult.explicitFallbackReason == mx::improc::P4PCAFallbackReason::none );
        REQUIRE( downdatedResult.numericalRank == explicitResult.numericalRank );
        REQUIRE( downdatedResult.modeStatus == explicitResult.modeStatus );
        REQUIRE( downdatedResult.sampleValidity.rows() == explicitResult.sampleValidity.rows() );
        REQUIRE( downdatedResult.sampleValidity.cols() == explicitResult.sampleValidity.cols() );
        for( Eigen::Index mode = 0; mode < downdatedResult.residuals.cols(); ++mode )
        {
            for( Eigen::Index sample = 0; sample < downdatedResult.residuals.rows(); ++sample )
            {
                CAPTURE( sample, mode, predictors.rows(), predictors.cols() );
                REQUIRE( downdatedResult.sampleSupported( sample, static_cast<std::size_t>( mode ) ) ==
                         explicitResult.sampleSupported( sample, static_cast<std::size_t>( mode ) ) );
                if( downdatedResult.sampleSupported( sample, static_cast<std::size_t>( mode ) ) )
                {
                    REQUIRE( downdatedResult.residuals( sample, mode ) ==
                             Approx( explicitResult.residuals( sample, mode ) ).margin( 2e-9 ) );
                }
                else
                {
                    REQUIRE( mx::math::isNan( downdatedResult.residuals( sample, mode ) ) );
                }
            }
        }
    };

    SECTION( "T is less than K" )
    {
        pcaT::matrixT predictors( 4, 6 );
        predictors << 1, 0, 0, 0, 1, 2, 0, 1, 0, 0, 2, -1, 0, 0, 1, 0, -1, 1, 0, 0, 0, 1, 1, 0.5;
        pcaT::vectorT target( 4 );
        target << 2, -1, 3, 0.5;
        compare( predictors, target, { { 0, 2 }, { 1 }, { 0, 2 }, { 1, 3 } }, { 1, 2, 3 }, 4 );
    }

    SECTION( "T equals K" )
    {
        pcaT::matrixT predictors( 4, 4 );
        predictors << 2, 0, 1, -1, 1, 3, 0, 1, 0, 1, 2, 2, 1, -1, 1, 3;
        pcaT::vectorT target( 4 );
        target << 1, 4, -2, 3;
        compare( predictors, target, { { 0 }, { 1, 3 }, { 0, 2 }, { 3 } }, { 1, 2, 3 }, 4 );
    }

    SECTION( "K is less than T" )
    {
        pcaT::matrixT predictors( 6, 4 );
        predictors << 1, 0, 2, -1, 0, 1, -1, 2, 2, 1, 0, 1, -1, 2, 1, 0, 3, -1, 2, 2, 1, 2, -2, 1;
        pcaT::vectorT target( 6 );
        target << 2, -1, 3, 0.5, 4, -2;
        compare( predictors, target, { { 0 }, { 1, 4 }, { 0, 2 }, { 3, 5 }, { 1, 4 }, { 2, 5 } }, { 1, 2, 3 }, 4 );
    }

    SECTION( "rank-deficient T is less than K" )
    {
        pcaT::matrixT predictors( 4, 6 );
        predictors << 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 2, -1, 0, 0, 0, 0;
        pcaT::vectorT target( 4 );
        target << 2, -1, 3, 0.5;
        compare( predictors, target, { { 0 }, { 1 }, { 2 }, { 3 } }, { 1, 2 }, 2 );
    }

    SECTION( "rank-deficient T equals K" )
    {
        pcaT::matrixT predictors( 4, 4 );
        predictors << 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 2, 0, 2, -1, 1, 0;
        pcaT::vectorT target( 4 );
        target << 1, 4, -2, 3;
        compare( predictors, target, { { 0 }, { 1 }, { 2 }, { 3 } }, { 1, 2 }, 2 );
    }

    SECTION( "rank-deficient K is less than T" )
    {
        pcaT::matrixT predictors( 6, 3 );
        predictors << 1, 0, 1, 0, 1, 1, 1, 1, 2, 2, -1, 1, -1, 2, 1, 3, 1, 4;
        pcaT::vectorT target( 6 );
        target << 2, -1, 3, 0.5, 4, -2;
        compare( predictors, target, { { 0 }, { 1 }, { 2 }, { 3 }, { 4 }, { 5 } }, { 1, 2 }, 2 );
    }

    SECTION( "every target has no retained row" )
    {
        pcaT::matrixT predictors( 3, 4 );
        predictors << 1, 0, 2, 1, 0, 1, -1, 2, 2, 1, 0, -1;
        pcaT::vectorT target( 3 );
        target << 2, -1, 3;
        const auto exclusions =
            mx::improc::P4TargetExclusions::fromExplicit( 3, { { 0, 1, 2 }, { 0, 1, 2 }, { 0, 1, 2 } } );
        resultT result;
        solverReset reset( fakeSolverBehavior::failure );
        mx::improc::P4PCA::calculateHeldOutDowndated( result,
                                                      predictors,
                                                      target,
                                                      exclusions,
                                                      { 1, 2 },
                                                      1e-12,
                                                      baseWorkspace,
                                                      downdateWorkspace );
        REQUIRE( solverCalls == 0 );
        REQUIRE( result.baseRank == 0 );
        REQUIRE( result.numericalRank == 0 );
        REQUIRE( allNan( result.residuals ) );
        REQUIRE( result.sampleValidity.count() == 0 );
    }
}

/// Verify structured one-row deletion agrees with dense deletion and explicit refits in every Gram orientation.
/** This directly exercises mx::improc::P4PCA::calculateHeldOutDowndated() with
 * mx::math::svdDeletionBackend::rankOneSecular for every target-row position when T is less than, equal to, and
 * greater than K.
 *
 * \ingroup P4PCA_unit_tests
 */
TEST_CASE( "P4PCA structured one-row deletion matches dense and explicit held-out refits",
           "[P4PCA][held-out][downdate][structured][exact]" )
{
    const auto compare = []( const pcaT::matrixT &predictors, const pcaT::vectorT &target )
    {
        std::vector<std::vector<Eigen::Index>> deleted( static_cast<std::size_t>( predictors.rows() ) );
        for( Eigen::Index targetIndex = 0; targetIndex < predictors.rows(); ++targetIndex )
        {
            deleted[static_cast<std::size_t>( targetIndex )].push_back( targetIndex );
        }
        const auto exclusions = mx::improc::P4TargetExclusions::fromExplicit( predictors.rows(), deleted );
        const int retainedRank = static_cast<int>( std::min( predictors.rows() - 1, predictors.cols() ) );
        std::vector<int> modes;
        for( int mode = 1; mode <= retainedRank; ++mode )
        {
            modes.push_back( mode );
        }

        resultT explicitResult;
        resultT denseResult;
        resultT structuredResult;
        pcaT::workspaceT explicitWorkspace;
        pcaT::workspaceT denseBaseWorkspace;
        pcaT::workspaceT structuredBaseWorkspace;
        mx::improc::P4PCADowndateWorkspace denseDowndateWorkspace;
        mx::improc::P4PCADowndateWorkspace structuredDowndateWorkspace;
        mx::improc::P4PCA::calculateHeldOut( explicitResult,
                                             predictors,
                                             target,
                                             exclusions,
                                             modes,
                                             1e-12,
                                             explicitWorkspace );
        mx::improc::P4PCA::calculateHeldOutDowndated( denseResult,
                                                      predictors,
                                                      target,
                                                      exclusions,
                                                      modes,
                                                      1e-12,
                                                      mx::math::svdDeletionBackend::leadingCovariance,
                                                      denseBaseWorkspace,
                                                      denseDowndateWorkspace );
        mx::improc::P4PCA::calculateHeldOutDowndated( structuredResult,
                                                      predictors,
                                                      target,
                                                      exclusions,
                                                      modes,
                                                      1e-12,
                                                      mx::math::svdDeletionBackend::rankOneSecular,
                                                      structuredBaseWorkspace,
                                                      structuredDowndateWorkspace );

        REQUIRE( denseResult.explicitFallbackCount == 0 );
        REQUIRE( denseResult.explicitFallbackReason == mx::improc::P4PCAFallbackReason::none );
        REQUIRE( structuredResult.explicitFallbackCount == 0 );
        REQUIRE( structuredResult.explicitFallbackReason == mx::improc::P4PCAFallbackReason::none );
        REQUIRE( structuredResult.baseRank == denseResult.baseRank );
        REQUIRE( structuredResult.numericalRank == explicitResult.numericalRank );
        REQUIRE( structuredResult.numericalRank == denseResult.numericalRank );
        REQUIRE( structuredResult.modeStatus == explicitResult.modeStatus );
        REQUIRE( structuredResult.modeStatus == denseResult.modeStatus );
        REQUIRE( ( structuredResult.sampleValidity == explicitResult.sampleValidity ).all() );
        REQUIRE( ( structuredResult.sampleValidity == denseResult.sampleValidity ).all() );
        for( Eigen::Index mode = 0; mode < structuredResult.residuals.cols(); ++mode )
        {
            for( Eigen::Index sample = 0; sample < structuredResult.residuals.rows(); ++sample )
            {
                CAPTURE( sample, mode, predictors.rows(), predictors.cols() );
                REQUIRE( structuredResult.sampleSupported( sample, static_cast<std::size_t>( mode ) ) );
                REQUIRE( structuredResult.residuals( sample, mode ) ==
                         Approx( explicitResult.residuals( sample, mode ) ).margin( 1e-8 ) );
                REQUIRE( structuredResult.residuals( sample, mode ) ==
                         Approx( denseResult.residuals( sample, mode ) ).margin( 1e-8 ) );
            }
        }
    };

    SECTION( "T is less than K" )
    {
        pcaT::matrixT predictors( 4, 6 );
        predictors << 4, 0.2, 0.1, 0, 0.3, -0.2, 0.1, 3, 0.4, 0.2, 0, 0.1, 0.2, -0.1, 2, 0.3, 0.1, 0.2, 0.1, 0.2, 0.3,
            1.5, -0.2, 0.1;
        pcaT::vectorT target( 4 );
        target << 2, -1, 3, 0.5;
        compare( predictors, target );
    }

    SECTION( "T equals K" )
    {
        pcaT::matrixT predictors( 4, 4 );
        predictors << 4, 0.2, 0.1, 0, 0.1, 3, 0.4, 0.2, 0.2, -0.1, 2, 0.3, 0.1, 0.2, 0.3, 1.5;
        pcaT::vectorT target( 4 );
        target << 1, 4, -2, 3;
        compare( predictors, target );
    }

    SECTION( "K is less than T" )
    {
        pcaT::matrixT predictors( 6, 4 );
        predictors << 4, 0.2, 0.1, 0, 0.1, 3, 0.4, 0.2, 0.2, -0.1, 2, 0.3, 0.1, 0.2, 0.3, 1.5, 1, -0.4, 0.2, 0.7, -0.3,
            0.5, 1.1, -0.2;
        pcaT::vectorT target( 6 );
        target << 2, -1, 3, 0.5, 4, -2;
        compare( predictors, target );
    }
}

/// Verify the structured backend rejects active multi-row deletion while preserving an all-rows-deleted target.
/** This directly exercises mx::improc::P4PCA::calculateHeldOutDowndated() with the one-row deletion contract and
 * compares the supported mixed edge case with mx::improc::P4PCA::calculateHeldOut().
 * \ingroup P4PCA_unit_tests
 */
TEST_CASE( "P4PCA structured deletion validates active deletion counts",
           "[P4PCA][held-out][downdate][structured][validation]" )
{
    pcaT::matrixT predictors( 3, 4 );
    predictors << 4, 0.2, 0.1, 0, 0.1, 3, 0.4, 0.2, 0.2, -0.1, 2, 0.3;
    pcaT::vectorT target( 3 );
    target << 2, -1, 3;

    SECTION( "active multi-row deletion is rejected" )
    {
        const auto exclusions = mx::improc::P4TargetExclusions::fromExplicit( 3, { { 0 }, { 0, 1 }, { 2 } } );
        resultT result;
        pcaT::workspaceT baseWorkspace;
        mx::improc::P4PCADowndateWorkspace downdateWorkspace;
        REQUIRE_THROWS_AS( mx::improc::P4PCA::calculateHeldOutDowndated( result,
                                                                         predictors,
                                                                         target,
                                                                         exclusions,
                                                                         { 1, 2 },
                                                                         1e-12,
                                                                         mx::math::svdDeletionBackend::rankOneSecular,
                                                                         baseWorkspace,
                                                                         downdateWorkspace ),
                           std::invalid_argument );
    }

    SECTION( "an all-rows-deleted target is skipped" )
    {
        const auto exclusions = mx::improc::P4TargetExclusions::fromExplicit( 3, { { 0 }, { 0, 1, 2 }, { 2 } } );
        resultT explicitResult;
        resultT structuredResult;
        pcaT::workspaceT explicitWorkspace;
        pcaT::workspaceT baseWorkspace;
        mx::improc::P4PCADowndateWorkspace downdateWorkspace;
        mx::improc::P4PCA::calculateHeldOut( explicitResult,
                                             predictors,
                                             target,
                                             exclusions,
                                             { 1, 2 },
                                             1e-12,
                                             explicitWorkspace );
        mx::improc::P4PCA::calculateHeldOutDowndated( structuredResult,
                                                      predictors,
                                                      target,
                                                      exclusions,
                                                      { 1, 2 },
                                                      1e-12,
                                                      mx::math::svdDeletionBackend::rankOneSecular,
                                                      baseWorkspace,
                                                      downdateWorkspace );

        REQUIRE( structuredResult.explicitFallbackCount == 0 );
        REQUIRE( structuredResult.explicitFallbackReason == mx::improc::P4PCAFallbackReason::none );
        REQUIRE( ( structuredResult.sampleValidity == explicitResult.sampleValidity ).all() );
        REQUIRE_FALSE( structuredResult.sampleSupported( 1, 0 ) );
        REQUIRE_FALSE( structuredResult.sampleSupported( 1, 1 ) );
        REQUIRE( allNan( structuredResult.residuals.row( 1 ) ) );
        for( const Eigen::Index sample : { Eigen::Index{ 0 }, Eigen::Index{ 2 } } )
        {
            for( Eigen::Index mode = 0; mode < structuredResult.residuals.cols(); ++mode )
            {
                REQUIRE( structuredResult.residuals( sample, mode ) ==
                         Approx( explicitResult.residuals( sample, mode ) ).margin( 1e-8 ) );
            }
        }
    }
}

/// Verify a structured secular-solver failure recomputes the complete search pixel with the explicit oracle.
/** This directly exercises mx::improc::P4PCA::calculateHeldOutProbeDowndated() through mxlib's LAED9 test seam and
 * compares both fallback products with mx::improc::P4PCA::calculateHeldOutProbe().
 * \ingroup P4PCA_unit_tests
 */
TEST_CASE( "P4PCA structured deletion solver failure uses the explicit held-out oracle",
           "[P4PCA][held-out][downdate][structured][fallback]" )
{
    pcaT::matrixT predictors( 3, 4 );
    predictors << 4, 0.2, 0.1, 0, 0.1, 3, 0.4, 0.2, 0.2, -0.1, 2, 0.3;
    pcaT::vectorT target( 3 );
    target << 2, -1, 3;
    pcaT::matrixT probePredictors( 2, 4 );
    probePredictors << 0.3, -0.2, 0.7, 0.1, -0.1, 0.8, 0.2, -0.3;
    pcaT::vectorT probeTarget( 2 );
    probeTarget << 1.2, -0.7;
    const auto exclusions = mx::improc::P4TargetExclusions::fromExplicit( 3, { { 0 }, { 1 }, { 2 } } );

    resultT explicitResult;
    pcaT::matrixT explicitProbe;
    pcaT::workspaceT explicitWorkspace;
    mx::improc::P4PCA::calculateHeldOutProbe( explicitResult,
                                              explicitProbe,
                                              predictors,
                                              target,
                                              probePredictors,
                                              probeTarget,
                                              exclusions,
                                              { 1, 2 },
                                              1e-12,
                                              explicitWorkspace );

    resultT structuredResult;
    pcaT::matrixT structuredProbe;
    pcaT::workspaceT baseWorkspace;
    mx::improc::P4PCADowndateWorkspace downdateWorkspace;
    timingT timing;
    svdDeletionHookReset reset;
    structuredDeletionSolveCalls = 0;
    structuredDeletionSolveInfo = 73;
    mx::math::detail::svdDeletionHooks<double>().laed9 = &failStructuredDeletionSolve;
    REQUIRE_NOTHROW( mx::improc::P4PCA::calculateHeldOutProbeDowndated( structuredResult,
                                                                        structuredProbe,
                                                                        predictors,
                                                                        target,
                                                                        probePredictors,
                                                                        probeTarget,
                                                                        exclusions,
                                                                        { 1, 2 },
                                                                        1e-12,
                                                                        mx::math::svdDeletionBackend::rankOneSecular,
                                                                        baseWorkspace,
                                                                        downdateWorkspace,
                                                                        &timing ) );

    REQUIRE( structuredDeletionSolveCalls == 1 );
    REQUIRE( structuredResult.explicitFallbackCount == static_cast<std::size_t>( predictors.rows() ) );
    REQUIRE( structuredResult.explicitFallbackReason == mx::improc::P4PCAFallbackReason::deletionSolver );
    REQUIRE( structuredResult.modeStatus == explicitResult.modeStatus );
    REQUIRE( ( structuredResult.sampleValidity == explicitResult.sampleValidity ).all() );
    for( Eigen::Index mode = 0; mode < structuredResult.residuals.cols(); ++mode )
    {
        for( Eigen::Index sample = 0; sample < structuredResult.residuals.rows(); ++sample )
        {
            REQUIRE( structuredResult.residuals( sample, mode ) ==
                     Approx( explicitResult.residuals( sample, mode ) ).margin( 1e-11 ) );
        }
    }
    REQUIRE( structuredProbe.rows() == explicitProbe.rows() );
    REQUIRE( structuredProbe.cols() == explicitProbe.cols() );
    for( Eigen::Index column = 0; column < structuredProbe.cols(); ++column )
    {
        for( Eigen::Index row = 0; row < structuredProbe.rows(); ++row )
        {
            REQUIRE( structuredProbe( row, column ) == Approx( explicitProbe( row, column ) ).margin( 1e-11 ) );
        }
    }
    REQUIRE( timing.baseFactorWorkerSeconds >= 0 );
    REQUIRE( timing.deletionWorkerSeconds >= 0 );
    REQUIRE( timing.explicitFallbackWorkerSeconds >= 0 );

    SECTION( "an illegal LAPACK argument remains fatal" )
    {
        resultT contractFailureResult;
        pcaT::workspaceT contractFailureBaseWorkspace;
        mx::improc::P4PCADowndateWorkspace contractFailureDowndateWorkspace;
        structuredDeletionSolveCalls = 0;
        structuredDeletionSolveInfo = -7;
        REQUIRE_THROWS_WITH( mx::improc::P4PCA::calculateHeldOutDowndated( contractFailureResult,
                                                                           predictors,
                                                                           target,
                                                                           exclusions,
                                                                           { 1, 2 },
                                                                           1e-12,
                                                                           mx::math::svdDeletionBackend::rankOneSecular,
                                                                           contractFailureBaseWorkspace,
                                                                           contractFailureDowndateWorkspace ),
                             Catch::Matchers::Contains( "status solverFailure" ) &&
                                 Catch::Matchers::Contains( "LAPACK status -7" ) );
        REQUIRE( structuredDeletionSolveCalls == 1 );
        REQUIRE( contractFailureResult.explicitFallbackReason == mx::improc::P4PCAFallbackReason::none );
    }
}

/// Verify P4 accepts ordinary factor error, explicitly refits orthogonality failures, and rejects other statuses.
/** This directly exercises mx::improc::P4PCA::calculateHeldOutDowndated() through its controlled eigensolver seam.
 * \ingroup P4PCA_unit_tests
 */
TEST_CASE( "P4PCA exact factor downdate handles factor validation outcomes",
           "[P4PCA][held-out][downdate][validation][fallback]" )
{
    pcaT::matrixT predictors( 3, 4 );
    predictors << 1, 0, 0, 0, 0, std::sqrt( 2.0 ), 0, 0, 0, 0, std::sqrt( 3.0 ), 0;
    pcaT::vectorT target( 3 );
    target << 2, -1, 3;
    const auto exclusions = mx::improc::P4TargetExclusions::fromExplicit( 3, { { 0 }, { 1 }, { 2 } } );

    resultT explicitResult;
    pcaT::workspaceT explicitWorkspace;
    mx::improc::P4PCA::calculateHeldOut( explicitResult,
                                         predictors,
                                         target,
                                         exclusions,
                                         { 1, 2 },
                                         1e-12,
                                         explicitWorkspace );
    const auto requireExplicitMatch = [&]( const resultT &result )
    {
        REQUIRE( result.modeStatus == explicitResult.modeStatus );
        REQUIRE( ( result.sampleValidity == explicitResult.sampleValidity ).all() );
        for( Eigen::Index mode = 0; mode < result.residuals.cols(); ++mode )
        {
            for( Eigen::Index sample = 0; sample < result.residuals.rows(); ++sample )
            {
                REQUIRE( result.residuals( sample, mode ) ==
                         Approx( explicitResult.residuals( sample, mode ) ).margin( 1e-11 ) );
            }
        }
    };

    SECTION( "DSYEVR-scale defect is accepted" )
    {
        resultT downdatedResult;
        pcaT::workspaceT baseWorkspace;
        mx::improc::P4PCADowndateWorkspace downdateWorkspace;
        solverReset reset( fakeSolverBehavior::marginallyNonorthonormalEigenvectors );

        REQUIRE_NOTHROW( mx::improc::P4PCA::calculateHeldOutDowndated( downdatedResult,
                                                                       predictors,
                                                                       target,
                                                                       exclusions,
                                                                       { 1, 2 },
                                                                       1e-12,
                                                                       baseWorkspace,
                                                                       downdateWorkspace ) );
        REQUIRE( downdatedResult.explicitFallbackCount == 0 );
        REQUIRE( downdatedResult.explicitFallbackReason == mx::improc::P4PCAFallbackReason::none );
        REQUIRE( downdatedResult.factorOrthogonalityDefect == 0 );
        REQUIRE( downdatedResult.factorOrthogonalityTolerance == 0 );
        requireExplicitMatch( downdatedResult );
    }

    SECTION( "defect above the P4 acceptance bound uses the explicit oracle" )
    {
        resultT result;
        pcaT::workspaceT eigensolverWorkspace;
        mx::improc::P4PCADowndateWorkspace downdateWorkspace;
        solverReset reset( fakeSolverBehavior::excessivelyNonorthonormalEigenvectors, true );

        REQUIRE_NOTHROW( mx::improc::P4PCA::calculateHeldOutDowndated( result,
                                                                       predictors,
                                                                       target,
                                                                       exclusions,
                                                                       { 1, 2 },
                                                                       1e-12,
                                                                       eigensolverWorkspace,
                                                                       downdateWorkspace ) );
        REQUIRE( solverCalls > 1 );
        REQUIRE( result.explicitFallbackCount == static_cast<std::size_t>( predictors.rows() ) );
        REQUIRE( result.explicitFallbackReason == mx::improc::P4PCAFallbackReason::factorValidation );
        REQUIRE( result.factorOrthogonalityDefect > result.factorOrthogonalityTolerance );
        REQUIRE( result.factorOrthogonalityTolerance > 0 );
        requireExplicitMatch( result );
    }

    SECTION( "material defect also uses the explicit oracle without relaxing tolerance" )
    {
        resultT result;
        pcaT::workspaceT eigensolverWorkspace;
        mx::improc::P4PCADowndateWorkspace downdateWorkspace;
        solverReset reset( fakeSolverBehavior::nonorthonormalEigenvectors, true );

        REQUIRE_NOTHROW( mx::improc::P4PCA::calculateHeldOutDowndated( result,
                                                                       predictors,
                                                                       target,
                                                                       exclusions,
                                                                       { 1, 2 },
                                                                       1e-12,
                                                                       eigensolverWorkspace,
                                                                       downdateWorkspace ) );
        REQUIRE( result.explicitFallbackCount == static_cast<std::size_t>( predictors.rows() ) );
        REQUIRE( result.explicitFallbackReason == mx::improc::P4PCAFallbackReason::factorValidation );
        REQUIRE( result.factorOrthogonalityDefect == Approx( 3.0 ) );
        REQUIRE( result.factorOrthogonalityDefect > result.factorOrthogonalityTolerance );
        requireExplicitMatch( result );
    }

    SECTION( "non-factor validation status remains a hard error" )
    {
        resultT result;
        pcaT::workspaceT eigensolverWorkspace;
        mx::improc::P4PCADowndateWorkspace downdateWorkspace;
        svdDeletionHookReset reset;
        mx::math::detail::svdDeletionHooks<double>().operation = &failFactorValidationAllocation;

        try
        {
            mx::improc::P4PCA::calculateHeldOutDowndated( result,
                                                          predictors,
                                                          target,
                                                          exclusions,
                                                          { 1, 2 },
                                                          1e-12,
                                                          eigensolverWorkspace,
                                                          downdateWorkspace );
            FAIL( "expected non-factor validation failure" );
        }
        catch( const std::runtime_error &error )
        {
            const std::string message = error.what();
            REQUIRE( message.find( "status allocationFailure" ) != std::string::npos );
        }
    }
}

/// Verify threshold-adjacent downdated spectra are recomputed by the explicit held-out oracle.
/** This exercises mx::improc::P4PCA::calculateHeldOutDowndated() and mx::improc::P4PCA::calculateHeldOut() in both
 * Gram orientations at an exactly equal strict rank boundary.
 * \ingroup P4PCA_unit_tests
 */
TEST_CASE( "P4PCA exact factor downdate records rank-boundary fallback", "[P4PCA][held-out][downdate][rank][fallback]" )
{
    const auto compareFallback = []( const pcaT::matrixT &predictors, const pcaT::vectorT &target )
    {
        std::vector<std::vector<Eigen::Index>> deleted( static_cast<std::size_t>( predictors.rows() ) );
        for( Eigen::Index targetIndex = 0; targetIndex < predictors.rows(); ++targetIndex )
        {
            deleted[static_cast<std::size_t>( targetIndex )].push_back( targetIndex );
        }
        const auto exclusions = mx::improc::P4TargetExclusions::fromExplicit( predictors.rows(), deleted );
        pcaT::workspaceT explicitWorkspace;
        pcaT::workspaceT baseWorkspace;
        mx::improc::P4PCADowndateWorkspace downdateWorkspace;
        resultT explicitResult;
        resultT downdatedResult;
        mx::improc::P4PCA::calculateHeldOut( explicitResult,
                                             predictors,
                                             target,
                                             exclusions,
                                             { 1, 2 },
                                             0.5,
                                             explicitWorkspace );
        mx::improc::P4PCA::calculateHeldOutDowndated( downdatedResult,
                                                      predictors,
                                                      target,
                                                      exclusions,
                                                      { 1, 2 },
                                                      0.5,
                                                      baseWorkspace,
                                                      downdateWorkspace );

        REQUIRE( downdatedResult.explicitFallbackCount == static_cast<std::size_t>( predictors.rows() ) );
        REQUIRE( downdatedResult.explicitFallbackReason == mx::improc::P4PCAFallbackReason::rankBoundary );
        REQUIRE( downdatedResult.factorOrthogonalityDefect == 0 );
        REQUIRE( downdatedResult.factorOrthogonalityTolerance == 0 );
        REQUIRE( downdatedResult.baseRank == std::min( predictors.rows(), predictors.cols() ) );
        REQUIRE( downdatedResult.numericalRank == explicitResult.numericalRank );
        REQUIRE( downdatedResult.modeStatus == explicitResult.modeStatus );
        REQUIRE( ( downdatedResult.sampleValidity == explicitResult.sampleValidity ).all() );
        for( Eigen::Index mode = 0; mode < explicitResult.residuals.cols(); ++mode )
        {
            for( Eigen::Index sample = 0; sample < explicitResult.residuals.rows(); ++sample )
            {
                if( explicitResult.sampleSupported( sample, static_cast<std::size_t>( mode ) ) )
                {
                    REQUIRE( downdatedResult.residuals( sample, mode ) ==
                             Approx( explicitResult.residuals( sample, mode ) ).margin( 1e-12 ) );
                }
                else
                {
                    REQUIRE( mx::math::isNan( downdatedResult.residuals( sample, mode ) ) );
                }
            }
        }
    };

    SECTION( "temporal Gram" )
    {
        pcaT::matrixT predictors = pcaT::matrixT::Zero( 3, 3 );
        predictors( 0, 0 ) = 2;
        predictors( 1, 1 ) = std::sqrt( 2.0 );
        predictors( 2, 2 ) = 1;
        pcaT::vectorT target( 3 );
        target << 1, -2, 3;
        compareFallback( predictors, target );
    }

    SECTION( "predictor Gram" )
    {
        pcaT::matrixT predictors = pcaT::matrixT::Zero( 4, 3 );
        predictors( 0, 0 ) = 2;
        predictors( 1, 1 ) = std::sqrt( 2.0 );
        predictors( 2, 2 ) = 1;
        pcaT::vectorT target( 4 );
        target << 1, -2, 3, 0.5;
        compareFallback( predictors, target );
    }
}

/// Verify complete-base factor deletion removes excluded target and predictor information from each fitted model.
/** This directly exercises mx::improc::P4PCA::calculateHeldOutDowndated() with noncontiguous multi-row exclusions.
 * \ingroup P4PCA_unit_tests
 */
TEST_CASE( "P4PCA exact factor downdate has no excluded-row leakage", "[P4PCA][held-out][downdate][leakage]" )
{
    pcaT::matrixT predictors( 5, 6 );
    predictors << 1, 0, 0, 2, -1, 1, 0, 1, 2, 0, 1, -1, 2, -1, 1, 0, 2, 1, -1, 2, 0, 1, 1, 0, 1, 1, -2, 1, 0, 2;
    pcaT::vectorT target( 5 );
    target << 2, -1, 3, 0.5, 4;
    const auto exclusions = mx::improc::P4TargetExclusions::fromExplicit( 5, { { 0 }, { 1 }, { 0, 2 }, { 3 }, { 4 } } );
    pcaT::workspaceT workspace;
    mx::improc::P4PCADowndateWorkspace downdateWorkspace;
    resultT baseline;
    mx::improc::P4PCA::calculateHeldOutDowndated( baseline,
                                                  predictors,
                                                  target,
                                                  exclusions,
                                                  { 1, 2 },
                                                  1e-12,
                                                  workspace,
                                                  downdateWorkspace );

    pcaT::vectorT changedTarget = target;
    changedTarget( 0 ) += 100;
    resultT targetChanged;
    mx::improc::P4PCA::calculateHeldOutDowndated( targetChanged,
                                                  predictors,
                                                  changedTarget,
                                                  exclusions,
                                                  { 1, 2 },
                                                  1e-12,
                                                  workspace,
                                                  downdateWorkspace );
    REQUIRE( targetChanged.residuals.row( 2 ).isApprox( baseline.residuals.row( 2 ), 2e-9 ) );

    pcaT::matrixT changedPredictors = predictors;
    changedPredictors.row( 0 ) *= 20;
    resultT predictorsChanged;
    mx::improc::P4PCA::calculateHeldOutDowndated( predictorsChanged,
                                                  changedPredictors,
                                                  target,
                                                  exclusions,
                                                  { 1, 2 },
                                                  1e-12,
                                                  workspace,
                                                  downdateWorkspace );
    REQUIRE( predictorsChanged.residuals.row( 2 ).isApprox( baseline.residuals.row( 2 ), 2e-8 ) );

    changedTarget = target;
    changedTarget( 2 ) += 7;
    resultT observedChanged;
    mx::improc::P4PCA::calculateHeldOutDowndated( observedChanged,
                                                  predictors,
                                                  changedTarget,
                                                  exclusions,
                                                  { 1, 2 },
                                                  1e-12,
                                                  workspace,
                                                  downdateWorkspace );
    for( Eigen::Index mode = 0; mode < baseline.residuals.cols(); ++mode )
    {
        const double baselinePrediction = target( 2 ) - baseline.residuals( 2, mode );
        const double changedPrediction = changedTarget( 2 ) - observedChanged.residuals( 2, mode );
        REQUIRE( changedPrediction == Approx( baselinePrediction ).margin( 2e-9 ) );
    }
}

/// Verify held-out deletion validation rejects target omission, duplicates, and inconsistent dimensions.
/** This directly exercises mx::improc::P4TargetExclusions, mx::improc::P4PCA::calculateHeldOut(), and
 * mx::improc::P4PCA::calculateHeldOutProbe().
 */
TEST_CASE( "P4PCA held-out input validation", "[P4PCA][held-out][validation]" )
{
    pcaT::matrixT predictors( 2, 1 );
    predictors << 1, 2;
    pcaT::vectorT target( 2 );
    target << 3, 4;
    resultT result;
    pcaT::workspaceT workspace;

    REQUIRE_THROWS_AS( mx::improc::P4TargetExclusions::fromExplicit( 2, { { 0 }, { 0 } } ), std::invalid_argument );
    REQUIRE_THROWS_AS( mx::improc::P4TargetExclusions::fromExplicit( 2, { { 0, 0 }, { 1 } } ), std::invalid_argument );
    REQUIRE_THROWS_AS( mx::improc::P4TargetExclusions::fromExplicit( 2, { { 0 } } ), std::invalid_argument );

    pcaT::matrixT widePredictors( 2, 3 );
    widePredictors << 1, 0, 2, 0, 1, -1;
    const auto mismatched = mx::improc::P4TargetExclusions::fromSpans( 3, { { 0, 1 }, { 1, 2 }, { 2, 3 } } );
    REQUIRE_THROWS_AS(
        mx::improc::P4PCA::calculateHeldOut( result, widePredictors, target, mismatched, { 1 }, 0, workspace ),
        std::invalid_argument );

    const auto exclusions = mx::improc::P4TargetExclusions::fromSpans( 2, { { 0, 1 }, { 1, 2 } } );
    pcaT::matrixT probePredictors( 1, 1 );
    probePredictors << 0.5;
    pcaT::vectorT probeTarget( 1 );
    probeTarget << 1;
    const pcaT::matrixT originalPredictors = predictors;
    const pcaT::matrixT originalProbePredictors = probePredictors;
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculateHeldOutProbe( result,
                                                                 predictors,
                                                                 predictors,
                                                                 target,
                                                                 probePredictors,
                                                                 probeTarget,
                                                                 exclusions,
                                                                 { 1 },
                                                                 0,
                                                                 workspace ),
                       std::invalid_argument );
    requireApprox( predictors, originalPredictors, 0 );
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculateHeldOutProbe( result,
                                                                 probePredictors,
                                                                 predictors,
                                                                 target,
                                                                 probePredictors,
                                                                 probeTarget,
                                                                 exclusions,
                                                                 { 1 },
                                                                 0,
                                                                 workspace ),
                       std::invalid_argument );
    requireApprox( probePredictors, originalProbePredictors, 0 );
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
        pcaT::matrixT coefficients;
        mx::improc::P4PCA::calculateCentered( result,
                                              predictors,
                                              target,
                                              { 1, 3 },
                                              1e-12,
                                              workspace,
                                              nullptr,
                                              &coefficients );

        REQUIRE( result.numericalRank == 3 );
        REQUIRE( result.modeStatus == std::vector<statusT>{ statusT::rankSupported, statusT::rankSupported } );
        requireApprox( result.residuals.col( 0 ), centeredFitUncenteredResidual( predictors, target, 1 ) );
        requireApprox( result.residuals.col( 1 ), centeredFitUncenteredResidual( predictors, target, 3 ) );
        requireApprox( coefficients.col( 0 ), centeredSvdCoefficients( predictors, target, 1 ) );
        requireApprox( coefficients.col( 1 ), centeredSvdCoefficients( predictors, target, 3 ) );
        requireApprox( ( target.matrix() - predictors.matrix() * coefficients.col( 0 ).matrix() ).array(),
                       result.residuals.col( 0 ) );
        requireApprox( ( target.matrix() - predictors.matrix() * coefficients.col( 1 ).matrix() ).array(),
                       result.residuals.col( 1 ) );
    }

    SECTION( "K equals or exceeds T" )
    {
        pcaT::matrixT predictors( 4, 5 );
        predictors << 1, 2, 0, 4, -1, 3, -1, 2, 0, 5, -2, 4, 1, -3, 2, 0, 1, -4, 2, 3;
        pcaT::vectorT target( 4 );
        target << 11, 3, 8, -2;

        resultT result;
        pcaT::matrixT coefficients;
        mx::improc::P4PCA::calculateCentered( result,
                                              predictors,
                                              target,
                                              { 1, 2, 3 },
                                              0,
                                              workspace,
                                              nullptr,
                                              &coefficients );

        REQUIRE( result.numericalRank == 3 );
        REQUIRE( result.modeStatus ==
                 std::vector<statusT>{ statusT::rankSupported, statusT::rankSupported, statusT::rankSupported } );
        requireApprox( result.residuals.col( 0 ), centeredFitUncenteredResidual( predictors, target, 1 ) );
        requireApprox( result.residuals.col( 1 ), centeredFitUncenteredResidual( predictors, target, 2 ) );
        requireApprox( result.residuals.col( 2 ), centeredFitUncenteredResidual( predictors, target, 3 ) );
        for( Eigen::Index output = 0; output < coefficients.cols(); ++output )
        {
            requireApprox( coefficients.col( output ), centeredSvdCoefficients( predictors, target, output + 1 ) );
            requireApprox( ( target.matrix() - predictors.matrix() * coefficients.col( output ).matrix() ).array(),
                           result.residuals.col( output ) );
        }
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

/// Verify destructive centered P4PCA matches the preserving API and reuses predictor storage.
/** This directly compares mx::improc::P4PCA::calculateCenteredInPlace() with
 * mx::improc::P4PCA::calculateCentered().
 */
TEST_CASE( "P4PCA in-place centered fit preserves results", "[P4PCA][centered][memory]" )
{
    pcaT::matrixT predictors( 5, 3 );
    predictors << 1, 8, -2, 2, 4, 1, 4, 3, 7, -1, 2, 5, 3, -2, 0;
    const pcaT::matrixT expectedCentered = centeredColumns( predictors );
    pcaT::matrixT inPlacePredictors = predictors;
    pcaT::vectorT target( 5 );
    target << 9, 4, 8, -1, 5;
    const pcaT::vectorT originalTarget = target;

    resultT preservingResult;
    resultT inPlaceResult;
    pcaT::workspaceT preservingWorkspace;
    pcaT::workspaceT inPlaceWorkspace;
    mx::improc::P4PCA::calculateCentered( preservingResult, predictors, target, { 1, 3 }, 1e-12, preservingWorkspace );
    mx::improc::P4PCA::calculateCenteredInPlace( inPlaceResult,
                                                 inPlacePredictors,
                                                 target,
                                                 { 1, 3 },
                                                 1e-12,
                                                 inPlaceWorkspace );

    REQUIRE( inPlaceResult.numericalRank == preservingResult.numericalRank );
    REQUIRE( inPlaceResult.modeStatus == preservingResult.modeStatus );
    requireApprox( inPlaceResult.residuals, preservingResult.residuals );
    requireApprox( inPlacePredictors, expectedCentered );
    requireApprox( target, originalTarget );
}

/// Verify centered P4PCA uses min(K,T-1) eigenpairs and the corresponding smaller Gram matrix.
/** This directly exercises mx::improc::P4PCA::calculateCentered() and
 * mx::improc::P4PCA::calculateCenteredInPlace() structural-rank dispatch through the controlled eigensolver seam.
 */
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

    SECTION( "in-place entry point" )
    {
        pcaT::matrixT predictors( 3, 5 );
        predictors << 1, 2, 0, 4, -1, 3, -1, 2, 0, 5, -2, 4, 1, -3, 2;
        pcaT::vectorT target( 3 );
        target << 7, 1, 3;
        const pcaT::matrixT centered = centeredColumns( predictors );
        const pcaT::matrixT expected = ( centered.matrix() * centered.matrix().transpose() ).array();

        solverReset reset( fakeSolverBehavior::failure );
        REQUIRE_THROWS_WITH(
            mx::improc::P4PCA::calculateCenteredInPlace( result, predictors, target, { 1, 2 }, 0, workspace ),
            "P4PCA eigensolver failed with status 37" );
        REQUIRE( solverCalls == 1 );
        REQUIRE( solverModeCount == 2 );
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
        pcaT::matrixT coefficients;
        mx::improc::P4PCA::calculateCentered( result,
                                              predictors,
                                              target,
                                              { 1, 2, 3 },
                                              1e-10,
                                              workspace,
                                              nullptr,
                                              &coefficients );
        REQUIRE( result.numericalRank == 1 );
        REQUIRE( result.modeStatus ==
                 std::vector<statusT>{ statusT::rankSupported, statusT::rankInsufficient, statusT::rankInsufficient } );
        REQUIRE( allFinite( result.residuals.col( 0 ) ) );
        REQUIRE( allNan( result.residuals.col( 1 ) ) );
        REQUIRE( allNan( result.residuals.col( 2 ) ) );
        REQUIRE( allFinite( coefficients.col( 0 ) ) );
        REQUIRE( allNan( coefficients.col( 1 ) ) );
        REQUIRE( allNan( coefficients.col( 2 ) ) );
    }

    SECTION( "constant predictors have zero centered rank" )
    {
        pcaT::matrixT predictors( 3, 2 );
        predictors << 2, -4, 2, -4, 2, -4;
        pcaT::vectorT target( 3 );
        target << 4, -2, 7;

        resultT result;
        pcaT::matrixT coefficients;
        mx::improc::P4PCA::calculateCentered( result,
                                              predictors,
                                              target,
                                              { 1, 2 },
                                              0,
                                              workspace,
                                              nullptr,
                                              &coefficients );
        REQUIRE( result.numericalRank == 0 );
        REQUIRE( result.modeStatus == std::vector<statusT>{ statusT::rankInsufficient, statusT::rankInsufficient } );
        REQUIRE( allNan( result.residuals ) );
        REQUIRE( coefficients.rows() == predictors.cols() );
        REQUIRE( coefficients.cols() == 2 );
        REQUIRE( allNan( coefficients ) );
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

/** Verify the production mixed-precision adapters cover direct, centered, held-out, and frozen-probe fits.
 * This exercises mx::improc::detail::p4PCACalculateMixed(),
 * mx::improc::detail::p4PCACalculateCenteredInPlaceMixed(),
 * mx::improc::detail::p4PCACalculateHeldOutMixed(), and
 * mx::improc::detail::p4PCACalculateHeldOutProbeMixed().
 */
TEST_CASE( "P4PCA production policy is FP32 calculation with FP64 eigensolve", "[P4PCA][precision][production][mixed]" )
{
    static_assert( std::is_same_v<typename mx::improc::detail::P4PCAFloatMatrixT::Scalar, float> );

    mx::improc::detail::P4PCAFloatMatrixT nativePredictors( 5, 3 );
    nativePredictors << 13, -4, 3, 8, 0, 2, 11, -5, 4, 10, -6, -2, 8, -5, 3;
    mx::improc::detail::P4PCAFloatVectorT nativeTarget( 5 );
    nativeTarget << 9, 4, 8, -1, 5;
    const pcaT::matrixT predictors = nativePredictors.cast<double>();
    const pcaT::vectorT target = nativeTarget.cast<double>();
    const std::vector<int> modes{ 1, 2 };

    mx::improc::detail::P4PCAMixedWorkspace workspace;
    static_assert( std::is_same_v<typename decltype( workspace.doubleEigensolver.cvd )::Scalar, double> );
    resultT direct;
    pcaT::matrixT coefficients;
    mx::improc::detail::p4PCACalculateMixed( direct,
                                             predictors,
                                             target,
                                             modes,
                                             1e-7,
                                             workspace,
                                             nullptr,
                                             &coefficients );
    REQUIRE( direct.numericalRank == 3 );
    REQUIRE( direct.modeStatus == std::vector<statusT>{ statusT::rankSupported, statusT::rankSupported } );
    for( Eigen::Index output = 0; output < static_cast<Eigen::Index>( modes.size() ); ++output )
    {
        requireApprox( direct.residuals.col( output ), svdResidual( predictors, target, modes[output] ), 2e-4 );
        requireApprox( coefficients.col( output ), svdCoefficients( predictors, target, modes[output] ), 2e-4 );
    }

    pcaT::matrixT centeredPredictors = predictors;
    resultT centered;
    mx::improc::detail::p4PCACalculateCenteredInPlaceMixed( centered,
                                                            centeredPredictors,
                                                            target,
                                                            modes,
                                                            1e-7,
                                                            workspace );
    requireApprox( centeredPredictors, centeredColumns( nativePredictors ).cast<double>(), 2e-6 );
    for( Eigen::Index output = 0; output < static_cast<Eigen::Index>( modes.size() ); ++output )
    {
        requireApprox( centered.residuals.col( output ),
                       centeredFitUncenteredResidual( predictors, target, modes[output] ),
                       2e-4 );
    }

    std::vector<mx::improc::P4ExclusionSpan> spans;
    for( Eigen::Index row = 0; row < predictors.rows(); ++row )
    {
        spans.push_back( { row, row + 1 } );
    }
    const mx::improc::P4TargetExclusions exclusions =
        mx::improc::P4TargetExclusions::fromSpans( predictors.rows(), spans );
    resultT heldOut;
    mx::improc::detail::p4PCACalculateHeldOutMixed( heldOut, predictors, target, exclusions, modes, 1e-7, workspace );
    resultT heldOutOracle;
    pcaT::workspaceT oracleWorkspace;
    pcaT::calculateHeldOut( heldOutOracle, predictors, target, exclusions, modes, 1e-7, oracleWorkspace );
    REQUIRE( ( heldOut.sampleValidity == heldOutOracle.sampleValidity ).all() );
    requireApprox( heldOut.residuals, heldOutOracle.residuals, 3e-4 );

    pcaT::matrixT probePredictors( 2, predictors.cols() );
    probePredictors << 1, -0.5, 0.25, -0.25, 0.75, 1.5;
    pcaT::vectorT probeTarget( 2 );
    probeTarget << 2, -1;
    pcaT::matrixT probeResiduals;
    resultT probe;
    mx::improc::detail::p4PCACalculateHeldOutProbeMixed( probe,
                                                         probeResiduals,
                                                         predictors,
                                                         target,
                                                         probePredictors,
                                                         probeTarget,
                                                         exclusions,
                                                         modes,
                                                         1e-7,
                                                         workspace );
    pcaT::matrixT probeOracleResiduals;
    resultT probeOracle;
    pcaT::calculateHeldOutProbe( probeOracle,
                                 probeOracleResiduals,
                                 predictors,
                                 target,
                                 probePredictors,
                                 probeTarget,
                                 exclusions,
                                 modes,
                                 1e-7,
                                 oracleWorkspace );
    REQUIRE( ( probe.sampleValidity == probeOracle.sampleValidity ).all() );
    requireApprox( probe.residuals, probeOracle.residuals, 3e-4 );
    requireApprox( probeResiduals, probeOracleResiduals, 3e-4 );
}

/** Verify the production mixed adapter preserves native FP64 eigensolver failures.
 * This exercises mx::improc::detail::p4PCACalculateMixed() through the production DSYEVR seam.
 */
TEST_CASE( "P4PCA production mixed policy preserves eigensolver failure status",
           "[P4PCA][precision][production][mixed][solver][failure]" )
{
    pcaT::matrixT predictors( 4, 3 );
    predictors << 1, 2, 0, 3, -1, 2, -2, 4, 1, 0, 1, -3;
    pcaT::vectorT target( 4 );
    target << 7, 1, 5, -2;
    mx::improc::detail::P4PCAMixedWorkspace workspace;
    resultT output;
    solverReset reset( fakeSolverBehavior::failure );
    REQUIRE_THROWS_WITH( mx::improc::detail::p4PCACalculateMixed( output, predictors, target, { 1, 2 }, 0, workspace ),
                         "P4PCA eigensolver failed with status 37" );
    REQUIRE( solverCalls == 1 );
}

#ifdef HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION

/** Verify every experimental scalar policy preserves direct-fit results and coefficients in both Gram orientations.
 * This exercises mx::improc::detail::p4PCACalculateExperimental() on identical FP32-quantized ingress values.
 */
TEST_CASE( "P4PCA experimental precision policies agree for direct fits", "[P4PCA][precision][experimental][direct]" )
{
    using policyT = mx::improc::detail::P4PCAPrecisionPolicy;
    constexpr std::array policies{ policyT::doubleDouble, policyT::floatDouble, policyT::floatFloat };

    SECTION( "predictor Gram" )
    {
        mx::improc::detail::P4PCAFloatMatrixT floatPredictors( 5, 3 );
        floatPredictors << 13, -4, 3, 8, 0, 2, 11, -5, 4, 10, -6, -2, 8, -5, 3;
        mx::improc::detail::P4PCAFloatVectorT floatTarget( 5 );
        floatTarget << 9, 4, 8, -1, 5;
        const pcaT::matrixT predictors = floatPredictors.cast<double>();
        const pcaT::vectorT target = floatTarget.cast<double>();
        const std::vector<int> modes{ 1, 2, 3 };

        resultT baseline;
        pcaT::matrixT baselineCoefficients;
        mx::improc::detail::P4PCAExperimentalWorkspace workspace;
        mx::improc::detail::p4PCACalculateExperimental( baseline,
                                                        predictors,
                                                        target,
                                                        modes,
                                                        1e-7,
                                                        policyT::doubleDouble,
                                                        workspace,
                                                        nullptr,
                                                        &baselineCoefficients );

        REQUIRE( baseline.numericalRank == 3 );
        for( Eigen::Index output = 0; output < static_cast<Eigen::Index>( modes.size() ); ++output )
        {
            requireApprox( baseline.residuals.col( output ), svdResidual( predictors, target, modes[output] ), 1e-10 );
            requireApprox( baselineCoefficients.col( output ),
                           svdCoefficients( predictors, target, modes[output] ),
                           1e-10 );
        }

        for( const policyT policy : policies )
        {
            INFO( "policy=" << static_cast<int>( policy ) );
            resultT candidate;
            pcaT::matrixT coefficients;
            mx::improc::detail::p4PCACalculateExperimental( candidate,
                                                            predictors,
                                                            target,
                                                            modes,
                                                            1e-7,
                                                            policy,
                                                            workspace,
                                                            nullptr,
                                                            &coefficients );

            REQUIRE( candidate.numericalRank == baseline.numericalRank );
            REQUIRE( candidate.modeStatus == baseline.modeStatus );
            requireApprox( candidate.residuals, baseline.residuals, 2e-4 );
            requireApprox( coefficients, baselineCoefficients, 2e-4 );
        }
    }

    SECTION( "temporal Gram" )
    {
        mx::improc::detail::P4PCAFloatMatrixT floatPredictors( 4, 5 );
        floatPredictors << 1, 2, 0, 4, -1, 3, -1, 2, 0, 5, -2, 4, 1, -3, 2, 0, 1, -4, 2, 3;
        mx::improc::detail::P4PCAFloatVectorT floatTarget( 4 );
        floatTarget << 11, 3, 8, -2;
        const pcaT::matrixT predictors = floatPredictors.cast<double>();
        const pcaT::vectorT target = floatTarget.cast<double>();
        const std::vector<int> modes{ 1, 2, 4 };

        resultT baseline;
        pcaT::matrixT baselineCoefficients;
        mx::improc::detail::P4PCAExperimentalWorkspace workspace;
        mx::improc::detail::p4PCACalculateExperimental( baseline,
                                                        predictors,
                                                        target,
                                                        modes,
                                                        0,
                                                        policyT::doubleDouble,
                                                        workspace,
                                                        nullptr,
                                                        &baselineCoefficients );

        REQUIRE( baseline.numericalRank == 4 );
        for( Eigen::Index output = 0; output < static_cast<Eigen::Index>( modes.size() ); ++output )
        {
            requireApprox( baseline.residuals.col( output ), svdResidual( predictors, target, modes[output] ), 1e-10 );
            requireApprox( baselineCoefficients.col( output ),
                           svdCoefficients( predictors, target, modes[output] ),
                           1e-10 );
        }

        for( const policyT policy : policies )
        {
            INFO( "policy=" << static_cast<int>( policy ) );
            resultT candidate;
            pcaT::matrixT coefficients;
            mx::improc::detail::p4PCACalculateExperimental( candidate,
                                                            predictors,
                                                            target,
                                                            modes,
                                                            0,
                                                            policy,
                                                            workspace,
                                                            nullptr,
                                                            &coefficients );

            REQUIRE( candidate.numericalRank == baseline.numericalRank );
            REQUIRE( candidate.modeStatus == baseline.modeStatus );
            requireApprox( candidate.residuals, baseline.residuals, 2e-4 );
            requireApprox( coefficients, baselineCoefficients, 2e-4 );
        }
    }
}

/** Verify every experimental scalar policy preserves exact rank loss and unsupported-plane semantics.
 * This exercises mx::improc::detail::p4PCACalculateExperimental() at zero rank tolerance in both Gram orientations.
 */
TEST_CASE( "P4PCA experimental precision policies agree for exact rank deficiency",
           "[P4PCA][precision][experimental][direct][rank]" )
{
    using policyT = mx::improc::detail::P4PCAPrecisionPolicy;
    constexpr std::array policies{ policyT::doubleDouble, policyT::floatDouble, policyT::floatFloat };

    const auto exerciseRankDeficiency = [&]( const mx::improc::detail::P4PCAFloatMatrixT &floatPredictors,
                                             const mx::improc::detail::P4PCAFloatVectorT &floatTarget )
    {
        const pcaT::matrixT predictors = floatPredictors.cast<double>();
        const pcaT::vectorT target = floatTarget.cast<double>();
        const std::vector<int> modes{ 1, 2, 3 };
        resultT baseline;
        pcaT::matrixT baselineCoefficients;
        mx::improc::detail::P4PCAExperimentalWorkspace workspace;
        mx::improc::detail::p4PCACalculateExperimental( baseline,
                                                        predictors,
                                                        target,
                                                        modes,
                                                        0,
                                                        policyT::doubleDouble,
                                                        workspace,
                                                        nullptr,
                                                        &baselineCoefficients );

        REQUIRE( baseline.numericalRank == 2 );
        REQUIRE( baseline.modeStatus ==
                 std::vector<statusT>{ statusT::rankSupported, statusT::rankSupported, statusT::rankInsufficient } );
        REQUIRE( allNan( baseline.residuals.col( 2 ) ) );
        REQUIRE( allNan( baselineCoefficients.col( 2 ) ) );

        for( const policyT policy : policies )
        {
            INFO( "policy=" << static_cast<int>( policy ) );
            resultT candidate;
            pcaT::matrixT coefficients;
            mx::improc::detail::p4PCACalculateExperimental( candidate,
                                                            predictors,
                                                            target,
                                                            modes,
                                                            0,
                                                            policy,
                                                            workspace,
                                                            nullptr,
                                                            &coefficients );

            REQUIRE( candidate.numericalRank == baseline.numericalRank );
            REQUIRE( candidate.modeStatus == baseline.modeStatus );
            requireApprox( candidate.residuals.leftCols( 2 ), baseline.residuals.leftCols( 2 ), 2e-6 );
            requireApprox( coefficients.leftCols( 2 ), baselineCoefficients.leftCols( 2 ), 2e-6 );
            REQUIRE( allNan( candidate.residuals.col( 2 ) ) );
            REQUIRE( allNan( coefficients.col( 2 ) ) );
        }
    };

    SECTION( "predictor Gram" )
    {
        mx::improc::detail::P4PCAFloatMatrixT predictors( 5, 3 );
        predictors << 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        mx::improc::detail::P4PCAFloatVectorT target( 5 );
        target << 3, 5, 7, -2, 1;
        exerciseRankDeficiency( predictors, target );
    }

    SECTION( "temporal Gram" )
    {
        mx::improc::detail::P4PCAFloatMatrixT predictors( 3, 5 );
        predictors << 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0;
        mx::improc::detail::P4PCAFloatVectorT target( 3 );
        target << 3, 5, 7;
        exerciseRankDeficiency( predictors, target );
    }
}

/** Verify every experimental scalar policy preserves centered and in-place results in both Gram orientations.
 * This exercises mx::improc::detail::p4PCACalculateCenteredExperimental() and
 * mx::improc::detail::p4PCACalculateCenteredInPlaceExperimental() on identical FP32-quantized ingress values.
 */
TEST_CASE( "P4PCA experimental precision policies agree for centered fits",
           "[P4PCA][precision][experimental][centered]" )
{
    using policyT = mx::improc::detail::P4PCAPrecisionPolicy;
    constexpr std::array policies{ policyT::doubleDouble, policyT::floatDouble, policyT::floatFloat };

    const auto exerciseShape = [&]( const mx::improc::detail::P4PCAFloatMatrixT &floatPredictors,
                                    const mx::improc::detail::P4PCAFloatVectorT &floatTarget,
                                    const std::vector<int> &modes,
                                    double rankTolerance )
    {
        const pcaT::matrixT predictors = floatPredictors.cast<double>();
        const pcaT::vectorT target = floatTarget.cast<double>();
        const pcaT::matrixT doubleCentered = centeredColumns( predictors );
        const pcaT::matrixT floatCentered = centeredColumns( floatPredictors ).cast<double>();

        resultT baseline;
        pcaT::matrixT baselineCoefficients;
        mx::improc::detail::P4PCAExperimentalWorkspace workspace;
        mx::improc::detail::p4PCACalculateCenteredExperimental( baseline,
                                                                predictors,
                                                                target,
                                                                modes,
                                                                rankTolerance,
                                                                policyT::doubleDouble,
                                                                workspace,
                                                                nullptr,
                                                                &baselineCoefficients );

        for( Eigen::Index output = 0; output < static_cast<Eigen::Index>( modes.size() ); ++output )
        {
            requireApprox( baseline.residuals.col( output ),
                           centeredFitUncenteredResidual( predictors, target, modes[output] ),
                           1e-10 );
            requireApprox( baselineCoefficients.col( output ),
                           centeredSvdCoefficients( predictors, target, modes[output] ),
                           1e-10 );
        }

        for( const policyT policy : policies )
        {
            INFO( "policy=" << static_cast<int>( policy ) );
            resultT preserving;
            pcaT::matrixT preservingCoefficients;
            mx::improc::detail::p4PCACalculateCenteredExperimental( preserving,
                                                                    predictors,
                                                                    target,
                                                                    modes,
                                                                    rankTolerance,
                                                                    policy,
                                                                    workspace,
                                                                    nullptr,
                                                                    &preservingCoefficients );

            pcaT::matrixT inPlacePredictors = predictors;
            resultT inPlace;
            pcaT::matrixT inPlaceCoefficients;
            mx::improc::detail::p4PCACalculateCenteredInPlaceExperimental( inPlace,
                                                                           inPlacePredictors,
                                                                           target,
                                                                           modes,
                                                                           rankTolerance,
                                                                           policy,
                                                                           workspace,
                                                                           nullptr,
                                                                           &inPlaceCoefficients );

            REQUIRE( preserving.numericalRank == baseline.numericalRank );
            REQUIRE( preserving.modeStatus == baseline.modeStatus );
            REQUIRE( inPlace.numericalRank == preserving.numericalRank );
            REQUIRE( inPlace.modeStatus == preserving.modeStatus );
            requireApprox( preserving.residuals, baseline.residuals, 3e-4 );
            requireApprox( preservingCoefficients, baselineCoefficients, 3e-4 );
            requireApprox( inPlace.residuals, preserving.residuals, 1e-10 );
            requireApprox( inPlaceCoefficients, preservingCoefficients, 1e-10 );
            requireApprox( inPlacePredictors, policy == policyT::doubleDouble ? doubleCentered : floatCentered, 2e-6 );
        }
    };

    SECTION( "predictor Gram" )
    {
        mx::improc::detail::P4PCAFloatMatrixT predictors( 5, 3 );
        predictors << 13, -4, 3, 8, 0, 2, 11, -5, 4, 10, -6, -2, 8, -5, 3;
        mx::improc::detail::P4PCAFloatVectorT target( 5 );
        target << 9, 4, 8, -1, 5;
        exerciseShape( predictors, target, { 1, 3 }, 1e-7 );
    }

    SECTION( "temporal Gram" )
    {
        mx::improc::detail::P4PCAFloatMatrixT predictors( 4, 5 );
        predictors << 1, 2, 0, 4, -1, 3, -1, 2, 0, 5, -2, 4, 1, -3, 2, 0, 1, -4, 2, 3;
        mx::improc::detail::P4PCAFloatVectorT target( 4 );
        target << 11, 3, 8, -2;
        exerciseShape( predictors, target, { 1, 2, 3 }, 0 );
    }
}

/** Verify every experimental scalar policy preserves explicit held-out science and frozen-probe responses.
 * This exercises mx::improc::detail::p4PCACalculateHeldOutExperimental() and
 * mx::improc::detail::p4PCACalculateHeldOutProbeExperimental() in both adaptive Gram orientations, including
 * non-contiguous mode requests, target-specific structural rank loss, and universally unsupported requested modes.
 */
TEST_CASE( "P4PCA experimental precision policies agree for held-out frozen probes",
           "[P4PCA][precision][experimental][held-out][probe]" )
{
    using policyT = mx::improc::detail::P4PCAPrecisionPolicy;
    constexpr std::array policies{ policyT::doubleDouble, policyT::floatDouble, policyT::floatFloat };

    const auto exerciseShape = [&]( const mx::improc::detail::P4PCAFloatMatrixT &floatPredictors,
                                    const mx::improc::detail::P4PCAFloatVectorT &floatTarget,
                                    const mx::improc::detail::P4PCAFloatMatrixT &floatProbePredictors,
                                    const mx::improc::detail::P4PCAFloatVectorT &floatProbeTarget,
                                    const std::vector<std::vector<std::size_t>> &retained )
    {
        const pcaT::matrixT predictors = floatPredictors.cast<double>();
        const pcaT::vectorT target = floatTarget.cast<double>();
        const pcaT::matrixT probePredictors = floatProbePredictors.cast<double>();
        const pcaT::vectorT probeTarget = floatProbeTarget.cast<double>();
        const mx::improc::P4TargetExclusions exclusions = exclusionsFromRetained( retained );
        const std::vector<int> modes{ 1, 3, 4 };

        resultT baseline;
        pcaT::matrixT baselineProbe;
        pcaT::workspaceT baselineWorkspace;
        mx::improc::P4PCA::calculateHeldOutProbe( baseline,
                                                  baselineProbe,
                                                  predictors,
                                                  target,
                                                  probePredictors,
                                                  probeTarget,
                                                  exclusions,
                                                  modes,
                                                  1e-6,
                                                  baselineWorkspace );

        REQUIRE( baseline.modeStatus.back() == statusT::rankInsufficient );
        REQUIRE( allNan( baseline.residuals.col( static_cast<Eigen::Index>( modes.size() - 1 ) ) ) );
        for( Eigen::Index targetIndex = 0; targetIndex < predictors.rows(); ++targetIndex )
        {
            REQUIRE( allNan( baselineProbe.col( targetIndex * static_cast<Eigen::Index>( modes.size() ) +
                                                static_cast<Eigen::Index>( modes.size() - 1 ) ) ) );
        }

        mx::improc::detail::P4PCAExperimentalWorkspace workspace;
        for( const policyT policy : policies )
        {
            INFO( "policy=" << static_cast<int>( policy ) );
            resultT probeCandidate;
            pcaT::matrixT probeCandidateResiduals;
            mx::improc::detail::p4PCACalculateHeldOutProbeExperimental( probeCandidate,
                                                                        probeCandidateResiduals,
                                                                        predictors,
                                                                        target,
                                                                        probePredictors,
                                                                        probeTarget,
                                                                        exclusions,
                                                                        modes,
                                                                        1e-6,
                                                                        policy,
                                                                        workspace );
            requireExperimentalHeldOutApprox( probeCandidate, probeCandidateResiduals, baseline, baselineProbe, 7e-4 );
            requireHeldOutProbeMatchesIndependent( probeCandidate,
                                                   probeCandidateResiduals,
                                                   predictors,
                                                   target,
                                                   probePredictors,
                                                   probeTarget,
                                                   retained,
                                                   modes,
                                                   7e-4 );

            resultT scienceCandidate;
            mx::improc::detail::p4PCACalculateHeldOutExperimental( scienceCandidate,
                                                                   predictors,
                                                                   target,
                                                                   exclusions,
                                                                   modes,
                                                                   1e-6,
                                                                   policy,
                                                                   workspace );
            REQUIRE( scienceCandidate.modeStatus == probeCandidate.modeStatus );
            REQUIRE( scienceCandidate.numericalRank == probeCandidate.numericalRank );
            REQUIRE( ( scienceCandidate.sampleValidity == probeCandidate.sampleValidity ).all() );
            for( Eigen::Index column = 0; column < scienceCandidate.residuals.cols(); ++column )
            {
                if( scienceCandidate.modeStatus[static_cast<std::size_t>( column )] == statusT::rankSupported )
                {
                    for( Eigen::Index row = 0; row < scienceCandidate.residuals.rows(); ++row )
                    {
                        if( scienceCandidate.sampleValidity( row, column ) )
                        {
                            REQUIRE( scienceCandidate.residuals( row, column ) ==
                                     Approx( probeCandidate.residuals( row, column ) ).margin( 1e-10 ) );
                        }
                        else
                        {
                            REQUIRE( mx::math::isNan( scienceCandidate.residuals( row, column ) ) );
                        }
                    }
                }
                else
                {
                    REQUIRE( allNan( scienceCandidate.residuals.col( column ) ) );
                }
            }
        }
    };

    SECTION( "predictor Gram" )
    {
        mx::improc::detail::P4PCAFloatMatrixT predictors( 6, 3 );
        predictors << 4, 0.2F, 0.1F, 0.1F, 3, 0.4F, 0.2F, -0.1F, 2, 0.1F, 0.2F, 0.3F, 1, -0.4F, 0.2F, -0.3F, 0.5F, 1.1F;
        mx::improc::detail::P4PCAFloatVectorT target( 6 );
        target << 2, -1, 3, 0.5F, 4, -2;
        mx::improc::detail::P4PCAFloatMatrixT probePredictors( 4, 3 );
        probePredictors << 0.3F, -0.2F, 0.7F, -0.1F, 0.8F, 0.2F, 0.5F, 0.1F, -0.2F, -0.4F, 0.6F, 0.9F;
        mx::improc::detail::P4PCAFloatVectorT probeTarget( 4 );
        probeTarget << 1.2F, -0.7F, 0.3F, 0.8F;
        const std::vector<std::vector<std::size_t>> retained{ { 1, 2, 3, 4, 5 },
                                                              { 0, 2, 3 },
                                                              { 0, 1 },
                                                              { 4 },
                                                              {},
                                                              { 0, 1, 2 } };
        exerciseShape( predictors, target, probePredictors, probeTarget, retained );
    }

    SECTION( "temporal Gram" )
    {
        mx::improc::detail::P4PCAFloatMatrixT predictors( 4, 6 );
        predictors << 4, 0.2F, 0.1F, 0, 0.3F, -0.2F, 0.1F, 3, 0.4F, 0.2F, 0, 0.1F, 0.2F, -0.1F, 2, 0.3F, 0.1F, 0.2F,
            0.1F, 0.2F, 0.3F, 1.5F, -0.2F, 0.1F;
        mx::improc::detail::P4PCAFloatVectorT target( 4 );
        target << 2, -1, 3, 0.5F;
        mx::improc::detail::P4PCAFloatMatrixT probePredictors( 3, 6 );
        probePredictors << 0.3F, -0.2F, 0.7F, 0.1F, -0.4F, 0.5F, -0.1F, 0.8F, 0.2F, -0.3F, 0.6F, 0.4F, 0.5F, 0.1F,
            -0.2F, 0.9F, 0.3F, -0.6F;
        mx::improc::detail::P4PCAFloatVectorT probeTarget( 3 );
        probeTarget << 1.2F, -0.7F, 0.3F;
        const std::vector<std::vector<std::size_t>> retained{ { 1, 2, 3 }, { 0, 2 }, { 3 }, {} };
        exerciseShape( predictors, target, probePredictors, probeTarget, retained );
    }
}

/** Verify positive rank thresholds retain only eigenvalues strictly above the boundary for every scalar policy.
 * This exercises mx::improc::detail::p4PCACalculateExperimental() with exactly representable spectra in both Gram
 * orientations.
 */
TEST_CASE( "P4PCA experimental precision policies preserve positive rank boundaries",
           "[P4PCA][precision][experimental][rank][boundary]" )
{
    using policyT = mx::improc::detail::P4PCAPrecisionPolicy;
    constexpr std::array policies{ policyT::doubleDouble, policyT::floatDouble, policyT::floatFloat };

    const auto exerciseShape = [&]( const mx::improc::detail::P4PCAFloatMatrixT &floatPredictors,
                                    const mx::improc::detail::P4PCAFloatVectorT &floatTarget )
    {
        const pcaT::matrixT predictors = floatPredictors.cast<double>();
        const pcaT::vectorT target = floatTarget.cast<double>();
        mx::improc::detail::P4PCAExperimentalWorkspace workspace;
        for( const policyT policy : policies )
        {
            INFO( "policy=" << static_cast<int>( policy ) );
            resultT output;
            pcaT::matrixT coefficients;
            mx::improc::detail::p4PCACalculateExperimental( output,
                                                            predictors,
                                                            target,
                                                            { 1, 2 },
                                                            0.25,
                                                            policy,
                                                            workspace,
                                                            nullptr,
                                                            &coefficients );
            REQUIRE( output.numericalRank == 1 );
            REQUIRE( output.modeStatus == std::vector<statusT>{ statusT::rankSupported, statusT::rankInsufficient } );
            REQUIRE( allFinite( output.residuals.col( 0 ) ) );
            REQUIRE( allNan( output.residuals.col( 1 ) ) );
            REQUIRE( allFinite( coefficients.col( 0 ) ) );
            REQUIRE( allNan( coefficients.col( 1 ) ) );
        }
    };

    SECTION( "predictor Gram" )
    {
        mx::improc::detail::P4PCAFloatMatrixT predictors( 3, 2 );
        predictors << 2, 0, 0, 1, 0, 0;
        mx::improc::detail::P4PCAFloatVectorT target( 3 );
        target << 3, -2, 1;
        exerciseShape( predictors, target );
    }

    SECTION( "temporal Gram" )
    {
        mx::improc::detail::P4PCAFloatMatrixT predictors( 2, 3 );
        predictors << 2, 0, 0, 0, 1, 0;
        mx::improc::detail::P4PCAFloatVectorT target( 2 );
        target << 3, -2;
        exerciseShape( predictors, target );
    }
}

/** Verify one experimental workspace can alternate scalar policies, dimensions, and Gram orientations.
 * This repeatedly exercises mx::improc::detail::p4PCACalculateExperimental() and compares each reused-workspace
 * result with the same policy in a fresh workspace.
 */
TEST_CASE( "P4PCA experimental precision workspace supports alternating policies and dimensions",
           "[P4PCA][precision][experimental][workspace]" )
{
    using policyT = mx::improc::detail::P4PCAPrecisionPolicy;
    mx::improc::detail::P4PCAExperimentalWorkspace reusedWorkspace;
    const auto exercise = [&]( const pcaT::matrixT &predictors,
                               const pcaT::vectorT &target,
                               const std::vector<int> &modes,
                               policyT policy )
    {
        resultT reused;
        pcaT::matrixT reusedCoefficients;
        mx::improc::detail::p4PCACalculateExperimental( reused,
                                                        predictors,
                                                        target,
                                                        modes,
                                                        1e-6,
                                                        policy,
                                                        reusedWorkspace,
                                                        nullptr,
                                                        &reusedCoefficients );

        mx::improc::detail::P4PCAExperimentalWorkspace freshWorkspace;
        resultT fresh;
        pcaT::matrixT freshCoefficients;
        mx::improc::detail::p4PCACalculateExperimental( fresh,
                                                        predictors,
                                                        target,
                                                        modes,
                                                        1e-6,
                                                        policy,
                                                        freshWorkspace,
                                                        nullptr,
                                                        &freshCoefficients );

        REQUIRE( reused.numericalRank == fresh.numericalRank );
        REQUIRE( reused.modeStatus == fresh.modeStatus );
        requireApprox( reused.residuals, fresh.residuals, 0 );
        requireApprox( reusedCoefficients, freshCoefficients, 0 );
    };

    pcaT::matrixT tall( 5, 3 );
    tall << 13, -4, 3, 8, 0, 2, 11, -5, 4, 10, -6, -2, 8, -5, 3;
    pcaT::vectorT tallTarget( 5 );
    tallTarget << 9, 4, 8, -1, 5;
    pcaT::matrixT wide( 4, 5 );
    wide << 1, 2, 0, 4, -1, 3, -1, 2, 0, 5, -2, 4, 1, -3, 2, 0, 1, -4, 2, 3;
    pcaT::vectorT wideTarget( 4 );
    wideTarget << 11, 3, 8, -2;
    pcaT::matrixT shortTall( 6, 2 );
    shortTall << 4, 1, 3, -2, 2, 0.5, -1, 3, 0.25, -0.5, 2.5, 1.5;
    pcaT::vectorT shortTallTarget( 6 );
    shortTallTarget << 2, -1, 3, 0.5, 4, -2;

    exercise( tall, tallTarget, { 1, 3 }, policyT::doubleDouble );
    exercise( wide, wideTarget, { 1, 2, 4 }, policyT::floatFloat );
    exercise( shortTall, shortTallTarget, { 1, 2 }, policyT::floatDouble );
    exercise( wide, wideTarget, { 1, 3 }, policyT::doubleDouble );
    exercise( tall, tallTarget, { 1, 2 }, policyT::floatDouble );
    exercise( shortTall, shortTallTarget, { 1, 2 }, policyT::floatFloat );
}

/** Verify every experimental entry point rejects invalid policies and FP64 values outside FP32 range.
 * This exercises the direct, centered, in-place, held-out, and frozen-probe adapters before unsafe FP32 arithmetic.
 */
TEST_CASE( "P4PCA experimental precision adapters reject invalid policies and FP32 overflow",
           "[P4PCA][precision][experimental][validation]" )
{
    using policyT = mx::improc::detail::P4PCAPrecisionPolicy;
    const policyT invalidPolicy = static_cast<policyT>( 255 );
    pcaT::matrixT predictors( 3, 2 );
    predictors << 2, 0, 0, 1, 1, -1;
    pcaT::vectorT target( 3 );
    target << 3, -2, 1;
    pcaT::matrixT probePredictors( 2, 2 );
    probePredictors << 0.5, -0.25, 1, 0.75;
    pcaT::vectorT probeTarget( 2 );
    probeTarget << 1, -1;
    const mx::improc::P4TargetExclusions exclusions =
        mx::improc::P4TargetExclusions::fromSpans( 3, { { 0, 1 }, { 1, 2 }, { 2, 3 } } );
    mx::improc::detail::P4PCAExperimentalWorkspace workspace;
    resultT output;
    pcaT::matrixT probeResiduals;

    REQUIRE_THROWS_AS( mx::improc::detail::p4PCACalculateExperimental( output,
                                                                       predictors,
                                                                       target,
                                                                       { 1 },
                                                                       0,
                                                                       invalidPolicy,
                                                                       workspace ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( mx::improc::detail::p4PCACalculateCenteredExperimental( output,
                                                                               predictors,
                                                                               target,
                                                                               { 1 },
                                                                               0,
                                                                               invalidPolicy,
                                                                               workspace ),
                       std::invalid_argument );
    pcaT::matrixT inPlacePredictors = predictors;
    REQUIRE_THROWS_AS( mx::improc::detail::p4PCACalculateCenteredInPlaceExperimental( output,
                                                                                      inPlacePredictors,
                                                                                      target,
                                                                                      { 1 },
                                                                                      0,
                                                                                      invalidPolicy,
                                                                                      workspace ),
                       std::invalid_argument );
    requireApprox( inPlacePredictors, predictors, 0 );
    REQUIRE_THROWS_AS( mx::improc::detail::p4PCACalculateHeldOutExperimental( output,
                                                                              predictors,
                                                                              target,
                                                                              exclusions,
                                                                              { 1 },
                                                                              0,
                                                                              invalidPolicy,
                                                                              workspace ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( mx::improc::detail::p4PCACalculateHeldOutProbeExperimental( output,
                                                                                   probeResiduals,
                                                                                   predictors,
                                                                                   target,
                                                                                   probePredictors,
                                                                                   probeTarget,
                                                                                   exclusions,
                                                                                   { 1 },
                                                                                   0,
                                                                                   invalidPolicy,
                                                                                   workspace ),
                       std::invalid_argument );

    pcaT::matrixT overflowingPredictors = predictors;
    overflowingPredictors( 0, 0 ) = std::numeric_limits<double>::max();
    pcaT::vectorT overflowingTarget = target;
    overflowingTarget( 0 ) = std::numeric_limits<double>::max();
    pcaT::matrixT overflowingProbePredictors = probePredictors;
    overflowingProbePredictors( 0, 0 ) = std::numeric_limits<double>::max();
    pcaT::vectorT overflowingProbeTarget = probeTarget;
    overflowingProbeTarget( 0 ) = std::numeric_limits<double>::max();

    for( const policyT policy : std::array{ policyT::doubleDouble, policyT::floatDouble, policyT::floatFloat } )
    {
        INFO( "alias policy=" << static_cast<int>( policy ) );
        pcaT::matrixT aliasedPredictors = predictors;
        const pcaT::matrixT originalAliasedPredictors = aliasedPredictors;
        REQUIRE_THROWS_AS( mx::improc::detail::p4PCACalculateHeldOutProbeExperimental( output,
                                                                                       aliasedPredictors,
                                                                                       aliasedPredictors,
                                                                                       target,
                                                                                       probePredictors,
                                                                                       probeTarget,
                                                                                       exclusions,
                                                                                       { 1 },
                                                                                       0,
                                                                                       policy,
                                                                                       workspace ),
                           std::invalid_argument );
        requireApprox( aliasedPredictors, originalAliasedPredictors, 0 );

        pcaT::matrixT aliasedProbePredictors = probePredictors;
        const pcaT::matrixT originalAliasedProbePredictors = aliasedProbePredictors;
        REQUIRE_THROWS_AS( mx::improc::detail::p4PCACalculateHeldOutProbeExperimental( output,
                                                                                       aliasedProbePredictors,
                                                                                       predictors,
                                                                                       target,
                                                                                       aliasedProbePredictors,
                                                                                       probeTarget,
                                                                                       exclusions,
                                                                                       { 1 },
                                                                                       0,
                                                                                       policy,
                                                                                       workspace ),
                           std::invalid_argument );
        requireApprox( aliasedProbePredictors, originalAliasedProbePredictors, 0 );
    }

    for( const policyT policy : std::array{ policyT::floatDouble, policyT::floatFloat } )
    {
        INFO( "policy=" << static_cast<int>( policy ) );
        REQUIRE_THROWS_AS( mx::improc::detail::p4PCACalculateExperimental( output,
                                                                           overflowingPredictors,
                                                                           target,
                                                                           { 1 },
                                                                           0,
                                                                           policy,
                                                                           workspace ),
                           std::invalid_argument );
        REQUIRE_THROWS_AS( mx::improc::detail::p4PCACalculateExperimental( output,
                                                                           predictors,
                                                                           overflowingTarget,
                                                                           { 1 },
                                                                           0,
                                                                           policy,
                                                                           workspace ),
                           std::invalid_argument );
        REQUIRE_THROWS_AS( mx::improc::detail::p4PCACalculateCenteredExperimental( output,
                                                                                   overflowingPredictors,
                                                                                   target,
                                                                                   { 1 },
                                                                                   0,
                                                                                   policy,
                                                                                   workspace ),
                           std::invalid_argument );
        REQUIRE_THROWS_AS( mx::improc::detail::p4PCACalculateCenteredExperimental( output,
                                                                                   predictors,
                                                                                   overflowingTarget,
                                                                                   { 1 },
                                                                                   0,
                                                                                   policy,
                                                                                   workspace ),
                           std::invalid_argument );

        pcaT::matrixT inPlaceOverflowingPredictors = overflowingPredictors;
        const pcaT::matrixT originalOverflowingPredictors = inPlaceOverflowingPredictors;
        REQUIRE_THROWS_AS( mx::improc::detail::p4PCACalculateCenteredInPlaceExperimental( output,
                                                                                          inPlaceOverflowingPredictors,
                                                                                          target,
                                                                                          { 1 },
                                                                                          0,
                                                                                          policy,
                                                                                          workspace ),
                           std::invalid_argument );
        requireApprox( inPlaceOverflowingPredictors, originalOverflowingPredictors, 0 );

        pcaT::matrixT inPlaceTargetOverflowPredictors = predictors;
        REQUIRE_THROWS_AS(
            mx::improc::detail::p4PCACalculateCenteredInPlaceExperimental( output,
                                                                           inPlaceTargetOverflowPredictors,
                                                                           overflowingTarget,
                                                                           { 1 },
                                                                           0,
                                                                           policy,
                                                                           workspace ),
            std::invalid_argument );
        requireApprox( inPlaceTargetOverflowPredictors, predictors, 0 );

        REQUIRE_THROWS_AS( mx::improc::detail::p4PCACalculateHeldOutExperimental( output,
                                                                                  overflowingPredictors,
                                                                                  target,
                                                                                  exclusions,
                                                                                  { 1 },
                                                                                  0,
                                                                                  policy,
                                                                                  workspace ),
                           std::invalid_argument );
        REQUIRE_THROWS_AS( mx::improc::detail::p4PCACalculateHeldOutExperimental( output,
                                                                                  predictors,
                                                                                  overflowingTarget,
                                                                                  exclusions,
                                                                                  { 1 },
                                                                                  0,
                                                                                  policy,
                                                                                  workspace ),
                           std::invalid_argument );
        REQUIRE_THROWS_AS( mx::improc::detail::p4PCACalculateHeldOutProbeExperimental( output,
                                                                                       probeResiduals,
                                                                                       predictors,
                                                                                       target,
                                                                                       overflowingProbePredictors,
                                                                                       probeTarget,
                                                                                       exclusions,
                                                                                       { 1 },
                                                                                       0,
                                                                                       policy,
                                                                                       workspace ),
                           std::invalid_argument );
        REQUIRE_THROWS_AS( mx::improc::detail::p4PCACalculateHeldOutProbeExperimental( output,
                                                                                       probeResiduals,
                                                                                       predictors,
                                                                                       target,
                                                                                       probePredictors,
                                                                                       overflowingProbeTarget,
                                                                                       exclusions,
                                                                                       { 1 },
                                                                                       0,
                                                                                       policy,
                                                                                       workspace ),
                           std::invalid_argument );
    }
}

/** Verify held-out and frozen-probe experimental adapters preserve exact native eigensolver statuses.
 * This exercises mx::improc::detail::p4PCACalculateHeldOutExperimental() and
 * mx::improc::detail::p4PCACalculateHeldOutProbeExperimental() in both Gram orientations.
 */
TEST_CASE( "P4PCA experimental held-out adapters preserve eigensolver failures",
           "[P4PCA][precision][experimental][held-out][probe][solver][failure]" )
{
    using policyT = mx::improc::detail::P4PCAPrecisionPolicy;
    constexpr std::array policies{ policyT::doubleDouble, policyT::floatDouble, policyT::floatFloat };

    const auto exerciseShape = [&]( const pcaT::matrixT &predictors, const pcaT::vectorT &target )
    {
        std::vector<mx::improc::P4ExclusionSpan> spans;
        for( Eigen::Index row = 0; row < predictors.rows(); ++row )
        {
            spans.push_back( { row, row + 1 } );
        }
        const mx::improc::P4TargetExclusions exclusions =
            mx::improc::P4TargetExclusions::fromSpans( predictors.rows(), spans );
        pcaT::matrixT probePredictors( 2, predictors.cols() );
        for( Eigen::Index column = 0; column < probePredictors.cols(); ++column )
        {
            probePredictors( 0, column ) = 0.25 * static_cast<double>( column + 1 );
            probePredictors( 1, column ) = -0.1 * static_cast<double>( column + 2 );
        }
        pcaT::vectorT probeTarget( 2 );
        probeTarget << 1, -1;

        for( const MXLAPACK_INT status : std::array<MXLAPACK_INT, 2>{ 37, -1000 } )
        {
            INFO( "status=" << status );
            experimentalSolverReset reset( status );
            mx::improc::detail::P4PCAExperimentalWorkspace workspace;
            const std::string expected = "P4PCA eigensolver failed with status " + std::to_string( status );
            for( const policyT policy : policies )
            {
                INFO( "policy=" << static_cast<int>( policy ) );
                resultT output;
                REQUIRE_THROWS_WITH( mx::improc::detail::p4PCACalculateHeldOutExperimental( output,
                                                                                            predictors,
                                                                                            target,
                                                                                            exclusions,
                                                                                            { 1, 2 },
                                                                                            0,
                                                                                            policy,
                                                                                            workspace ),
                                     expected );
                pcaT::matrixT probeResiduals;
                REQUIRE_THROWS_WITH( mx::improc::detail::p4PCACalculateHeldOutProbeExperimental( output,
                                                                                                 probeResiduals,
                                                                                                 predictors,
                                                                                                 target,
                                                                                                 probePredictors,
                                                                                                 probeTarget,
                                                                                                 exclusions,
                                                                                                 { 1, 2 },
                                                                                                 0,
                                                                                                 policy,
                                                                                                 workspace ),
                                     expected );
            }
        }
    };

    SECTION( "predictor Gram" )
    {
        pcaT::matrixT predictors( 4, 2 );
        predictors << 2, 0, 0, 1, 1, -1, 0.5, 2;
        pcaT::vectorT target( 4 );
        target << 3, -2, 1, 0.5;
        exerciseShape( predictors, target );
    }

    SECTION( "temporal Gram" )
    {
        pcaT::matrixT predictors( 3, 4 );
        predictors << 2, 0, 1, -1, 0, 1, 2, 0.5, 1, -1, 0, 2;
        pcaT::vectorT target( 3 );
        target << 3, -2, 1;
        exerciseShape( predictors, target );
    }
}

/** Verify experimental DD, mixed, and FP32 adapters preserve native eigensolver statuses without remapping.
 * This exercises all three experimental P4 entry points through both controlled native solver seams.
 */
TEST_CASE( "P4PCA experimental precision adapters preserve eigensolver failures",
           "[P4PCA][precision][experimental][solver][failure]" )
{
    using policyT = mx::improc::detail::P4PCAPrecisionPolicy;
    constexpr std::array policies{ policyT::doubleDouble, policyT::floatDouble, policyT::floatFloat };
    pcaT::matrixT predictors( 4, 3 );
    predictors << 1, 2, 0, 3, -1, 2, -2, 4, 1, 0, 1, -3;
    pcaT::vectorT target( 4 );
    target << 7, 1, 5, -2;

    for( const MXLAPACK_INT status : std::array<MXLAPACK_INT, 2>{ 37, -1000 } )
    {
        INFO( "status=" << status );
        experimentalSolverReset reset( status );
        mx::improc::detail::P4PCAExperimentalWorkspace workspace;
        const std::string expected = "P4PCA eigensolver failed with status " + std::to_string( status );

        for( const policyT policy : policies )
        {
            INFO( "direct policy=" << static_cast<int>( policy ) );
            resultT output;
            REQUIRE_THROWS_WITH( mx::improc::detail::p4PCACalculateExperimental( output,
                                                                                 predictors,
                                                                                 target,
                                                                                 { 1, 2 },
                                                                                 0,
                                                                                 policy,
                                                                                 workspace ),
                                 expected );
        }

        for( const policyT policy : policies )
        {
            INFO( "centered policy=" << static_cast<int>( policy ) );
            resultT output;
            REQUIRE_THROWS_WITH( mx::improc::detail::p4PCACalculateCenteredExperimental( output,
                                                                                         predictors,
                                                                                         target,
                                                                                         { 1, 2 },
                                                                                         0,
                                                                                         policy,
                                                                                         workspace ),
                                 expected );

            pcaT::matrixT inPlacePredictors = predictors;
            REQUIRE_THROWS_WITH( mx::improc::detail::p4PCACalculateCenteredInPlaceExperimental( output,
                                                                                                inPlacePredictors,
                                                                                                target,
                                                                                                { 1, 2 },
                                                                                                0,
                                                                                                policy,
                                                                                                workspace ),
                                 expected );
            pcaT::matrixT expectedCentered;
            if( policy == policyT::doubleDouble )
            {
                expectedCentered = centeredColumns( predictors );
            }
            else
            {
                const mx::improc::detail::P4PCAFloatMatrixT floatPredictors = predictors.cast<float>();
                expectedCentered = centeredColumns( floatPredictors ).cast<double>();
            }
            requireApprox( inPlacePredictors, expectedCentered, 2e-6 );
        }

        REQUIRE( solverCalls == 6 );
        REQUIRE( floatSolverCalls == 3 );
        REQUIRE( solverModeCount == 3 );
        REQUIRE( floatSolverModeCount == 3 );

        const mx::improc::detail::P4PCAFloatMatrixT floatPredictors = predictors.cast<float>();
        const mx::improc::detail::P4PCAFloatMatrixT floatCentered = centeredColumns( floatPredictors );
        const pcaT::matrixT expectedMixedGram =
            ( floatCentered.matrix().transpose() * floatCentered.matrix() ).template cast<double>().array();
        requireApprox( solverGram, expectedMixedGram, 0 );
        requireApprox( floatSolverGram, ( floatCentered.matrix().transpose() * floatCentered.matrix() ).array(), 0 );
    }
}

#endif // HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION

} // namespace P4PCA_test
} // namespace unitTest
