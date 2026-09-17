/** \file P4PCAResponse_test.cpp
 * \brief Establishes amplitude-convergence oracles for the direct P4 residual response.
 */

#include "../catch2/catch.hpp"

#include "src/common/P4PCA.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <vector>

namespace unitTest
{
namespace P4PCAResponse_test
{

/** \cond P4PCAResponse_test_harness */
using pcaT = mx::improc::P4PCA;
using resultT = mx::improc::P4PCAResult;
using statusT = mx::improc::P4PCAModeStatus;

/// Construct a complete orthonormal cosine basis without using the production eigensolver.
Eigen::MatrixXd cosineBasis( Eigen::Index size /**< [in] dimension of the square basis */ )
{
    Eigen::MatrixXd basis( size, size );
    for( Eigen::Index column = 0; column < size; ++column )
    {
        for( Eigen::Index row = 0; row < size; ++row )
        {
            basis( row, column ) = std::sqrt( ( column == 0 ? 1.0 : 2.0 ) / size ) *
                                   std::cos( std::numbers::pi * ( row + 0.5 ) * column / size );
        }
    }
    return basis;
}

/// Inputs whose exact temporal eigensystem is known, including the discarded nullspace.
struct ResponseFixture
{
    /// Baseline predictor time series.
    pcaT::matrixT m_predictors;

    /// Direction in which the predictors change with source amplitude.
    pcaT::matrixT m_predictorSource;

    /// Baseline target time series.
    pcaT::vectorT m_target;

    /// Direction in which the target changes with source amplitude.
    pcaT::vectorT m_targetSource;

    /// Complete temporal eigenbasis constructed independently of P4PCA.
    Eigen::MatrixXd m_basis;

    /// Descending temporal eigenvalues, including zeros outside the predictor rank.
    Eigen::VectorXd m_eigenvalues;

    /// Construct a deterministic fixture with prescribed nonzero and zero singular values.
    ResponseFixture( Eigen::Index rows,    /**< [in] temporal sample count */
                     Eigen::Index columns, /**< [in] predictor count */
                     const std::vector<double> &spectrum /**< [in] descending Gram eigenvalues */ );

    /// Evaluate the projector derivative using the independently prescribed eigensystem.
    pcaT::vectorT derivative( int modes /**< [in] retained count with a separated cutoff */ ) const;
};

ResponseFixture::ResponseFixture( Eigen::Index rows, Eigen::Index columns, const std::vector<double> &spectrum )
    : m_predictorSource( rows, columns ), m_target( rows ), m_targetSource( rows ), m_basis( cosineBasis( rows ) ),
      m_eigenvalues( Eigen::VectorXd::Zero( rows ) )
{
    const Eigen::Index count = static_cast<Eigen::Index>( spectrum.size() );
    REQUIRE( count <= std::min( rows, columns ) );
    Eigen::VectorXd singularValues( count );
    for( Eigen::Index index = 0; index < count; ++index )
    {
        m_eigenvalues( index ) = spectrum[index];
        singularValues( index ) = std::sqrt( spectrum[index] );
    }
    m_predictors = ( m_basis.leftCols( count ) * singularValues.asDiagonal() *
                     cosineBasis( columns ).leftCols( count ).transpose() )
                       .array();
    for( Eigen::Index row = 0; row < rows; ++row )
    {
        m_target( row ) = 1.2 + std::sin( 0.7 * ( row + 1 ) );
        m_targetSource( row ) = std::cos( 0.4 * ( row + 2 ) );
        for( Eigen::Index column = 0; column < columns; ++column )
        {
            m_predictorSource( row, column ) =
                0.4 * std::sin( 0.3 * ( row + 1 ) * ( column + 2 ) ) + 0.1 * std::cos( row + 2 * column );
        }
    }
}

pcaT::vectorT ResponseFixture::derivative( int modes ) const
{
    const Eigen::MatrixXd projector = m_basis.leftCols( modes ) * m_basis.leftCols( modes ).transpose();
    const Eigen::MatrixXd gramDerivative = m_predictorSource.matrix() * m_predictors.matrix().transpose() +
                                           m_predictors.matrix() * m_predictorSource.matrix().transpose();
    Eigen::MatrixXd projectorDerivative = Eigen::MatrixXd::Zero( m_basis.rows(), m_basis.rows() );
    for( Eigen::Index retained = 0; retained < modes; ++retained )
    {
        for( Eigen::Index discarded = modes; discarded < m_basis.cols(); ++discarded )
        {
            const double gap = m_eigenvalues( retained ) - m_eigenvalues( discarded );
            REQUIRE( gap > 0 );
            const double coupling = m_basis.col( discarded ).dot( gramDerivative * m_basis.col( retained ) ) / gap;
            projectorDerivative += coupling * ( m_basis.col( discarded ) * m_basis.col( retained ).transpose() +
                                                m_basis.col( retained ) * m_basis.col( discarded ).transpose() );
        }
    }
    return ( m_targetSource.matrix() - projector * m_targetSource.matrix() - projectorDerivative * m_target.matrix() )
        .array();
}
/** \endcond */

/** Verify mx::improc::P4PCA::calculate() converges quadratically to the complete projector derivative.
 * P4PCA::calculateResponse() must agree with the independent derivative. Both Gram orientations,
 * predictor-only perturbations, rank deficiency, and repeated eigenvalues wholly
 * inside retained or discarded groups are covered using a known complete temporal eigensystem.
 */
TEST_CASE( "Direct P4 paired refits converge to the complete residual derivative", "[P4PCA][response][convergence]" )
{
    const std::array<std::vector<double>, 3> spectra{ std::vector<double>{ 36, 25, 16, 9, 4, 1 },
                                                      { 25, 25, 9, 9, 1, 1 },
                                                      { 36, 25, 16, 9, 0, 0 } };
    const std::array<std::vector<int>, 3> modeSets{ std::vector<int>{ 1, 3, 5 }, { 2, 4 }, { 1, 3, 4 } };
    for( const bool temporalGram : { true, false } )
    {
        for( std::size_t spectrum = 0; spectrum < spectra.size(); ++spectrum )
        {
            for( const bool predictorOnly : { false, true } )
            {
                CAPTURE( temporalGram, spectrum, predictorOnly );
                ResponseFixture fixture( temporalGram ? 6 : 9, temporalGram ? 9 : 6, spectra[spectrum] );
                if( predictorOnly )
                {
                    fixture.m_targetSource.setZero();
                }
                const auto &modes = modeSets[spectrum];
                Eigen::MatrixXd expected( fixture.m_target.size(), modes.size() );
                for( std::size_t mode = 0; mode < modes.size(); ++mode )
                {
                    expected.col( mode ) = fixture.derivative( modes[mode] ).matrix();
                }
                pcaT::workspaceT workspace;
                mx::improc::P4PCAResponseResult analytic;
                mx::improc::P4PCA::calculateResponse( analytic,
                                                      fixture.m_predictors,
                                                      fixture.m_target,
                                                      fixture.m_predictorSource,
                                                      fixture.m_targetSource,
                                                      modes,
                                                      1e-12,
                                                      workspace );
                REQUIRE( analytic.modeStatus == std::vector<mx::improc::P4PCAResponseStatus>(
                                                    modes.size(),
                                                    mx::improc::P4PCAResponseStatus::differentiable ) );
                REQUIRE( ( analytic.responses.matrix() - expected ).norm() / expected.norm() < 5e-12 );
                Eigen::VectorXd previousError = Eigen::VectorXd::Zero( modes.size() );
                for( const double amplitude : { 0.04, 0.02, 0.01, 1e-5 } )
                {
                    resultT positive, negative;
                    mx::improc::P4PCA::calculate( positive,
                                                  fixture.m_predictors + amplitude * fixture.m_predictorSource,
                                                  fixture.m_target + amplitude * fixture.m_targetSource,
                                                  modes,
                                                  1e-12,
                                                  workspace );
                    mx::improc::P4PCA::calculate( negative,
                                                  fixture.m_predictors - amplitude * fixture.m_predictorSource,
                                                  fixture.m_target - amplitude * fixture.m_targetSource,
                                                  modes,
                                                  1e-12,
                                                  workspace );
                    REQUIRE( positive.modeStatus == std::vector<statusT>( modes.size(), statusT::rankSupported ) );
                    REQUIRE( negative.modeStatus == positive.modeStatus );
                    const Eigen::MatrixXd measured =
                        ( ( positive.residuals - negative.residuals ) / ( 2 * amplitude ) ).matrix();
                    for( std::size_t mode = 0; mode < modes.size(); ++mode )
                    {
                        const double error =
                            ( measured.col( mode ) - expected.col( mode ) ).norm() / expected.col( mode ).norm();
                        CAPTURE( modes[mode], amplitude, error );
                        if( amplitude == 1e-5 )
                        {
                            REQUIRE( error < 5e-8 );
                        }
                        else if( amplitude < 0.04 )
                        {
                            REQUIRE( error < 0.27 * previousError( mode ) );
                            REQUIRE( error > 0.23 * previousError( mode ) );
                        }
                        previousError( mode ) = error;
                    }
                }
            }
        }
    }
}

/** Verify mx::improc::P4PCA::calculate() refits the target coefficients even when predictors are unchanged.
 * P4PCA::calculateResponse() returns (I-Pi)s rather than the frozen-coefficient response s.
 */
TEST_CASE( "Direct P4 target-only response includes coefficient adaptation", "[P4PCA][response][target-only]" )
{
    ResponseFixture fixture( 9, 6, { 36, 25, 16, 9, 4, 1 } );
    fixture.m_predictorSource.setZero();
    pcaT::workspaceT workspace;
    for( const int modes : { 1, 3, 5 } )
    {
        const pcaT::vectorT expected = fixture.derivative( modes );
        REQUIRE( ( expected - fixture.m_targetSource ).matrix().norm() > 0.1 );
        mx::improc::P4PCAResponseResult analytic;
        mx::improc::P4PCA::calculateResponse( analytic,
                                              fixture.m_predictors,
                                              fixture.m_target,
                                              fixture.m_predictorSource,
                                              fixture.m_targetSource,
                                              { modes },
                                              1e-12,
                                              workspace );
        REQUIRE( ( analytic.responses.col( 0 ) - expected ).matrix().norm() < 1e-12 );
        for( const double amplitude : { 0.1, 0.001, 1e-5 } )
        {
            resultT positive, negative;
            mx::improc::P4PCA::calculate( positive,
                                          fixture.m_predictors,
                                          fixture.m_target + amplitude * fixture.m_targetSource,
                                          { modes },
                                          1e-12,
                                          workspace );
            mx::improc::P4PCA::calculate( negative,
                                          fixture.m_predictors,
                                          fixture.m_target - amplitude * fixture.m_targetSource,
                                          { modes },
                                          1e-12,
                                          workspace );
            const pcaT::vectorT measured =
                ( positive.residuals.col( 0 ) - negative.residuals.col( 0 ) ) / ( 2 * amplitude );
            CAPTURE( modes, amplitude );
            REQUIRE( ( measured - expected ).matrix().norm() < 1e-9 );
        }
    }
}

/** Verify mx::improc::detail::p4PCACalculateMixed() has a usable finite-amplitude window in both Gram orientations.
 * The reference uses the prescribed FP64 eigensystem; no assumption is made that the smallest step is best.
 */
TEST_CASE( "Mixed P4 paired refits have a finite-amplitude convergence window", "[P4PCA][response][mixed]" )
{
    for( const bool temporalGram : { true, false } )
    {
        ResponseFixture fixture( temporalGram ? 6 : 9, temporalGram ? 9 : 6, { 36, 25, 16, 9, 4, 1 } );
        mx::improc::detail::P4PCAMixedWorkspace workspace;
        for( const int modes : { 1, 3, 5 } )
        {
            const pcaT::vectorT expected = fixture.derivative( modes );
            double bestError = std::numeric_limits<double>::infinity();
            std::vector<double> errors;
            for( const double amplitude : { 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0004, 0.00008, 0.000016 } )
            {
                resultT positive, negative;
                mx::improc::detail::p4PCACalculateMixed( positive,
                                                         fixture.m_predictors + amplitude * fixture.m_predictorSource,
                                                         fixture.m_target + amplitude * fixture.m_targetSource,
                                                         { modes },
                                                         1e-6,
                                                         workspace );
                mx::improc::detail::p4PCACalculateMixed( negative,
                                                         fixture.m_predictors - amplitude * fixture.m_predictorSource,
                                                         fixture.m_target - amplitude * fixture.m_targetSource,
                                                         { modes },
                                                         1e-6,
                                                         workspace );
                REQUIRE( positive.modeStatus == std::vector<statusT>{ statusT::rankSupported } );
                REQUIRE( negative.modeStatus == positive.modeStatus );
                const pcaT::vectorT measured =
                    ( positive.residuals.col( 0 ) - negative.residuals.col( 0 ) ) / ( 2 * amplitude );
                const double error = ( measured - expected ).matrix().norm() / expected.matrix().norm();
                REQUIRE( std::isfinite( error ) );
                errors.push_back( error );
                bestError = std::min( bestError, error );
            }
            CAPTURE( temporalGram, modes, errors, bestError );
            // A finite-difference oracle in the production precision needs a looser bound than the FP64 limit.
            REQUIRE( bestError < 5e-3 );
        }
    }
}

/** Verify FP32 ingress can erase a source that mx::improc::P4PCA::calculate() still resolves in FP64.
 * This exercises mx::improc::detail::p4PCACalculateMixed() with an exactly diagonal, well-separated Gram matrix.
 */
TEST_CASE( "Mixed P4 tiny paired refits can cancel completely", "[P4PCA][response][mixed][cancellation]" )
{
    const pcaT::matrixT predictors = ( Eigen::Vector3d( 3, 2, 1 ).asDiagonal() ).toDenseMatrix().array();
    const pcaT::vectorT target = Eigen::Vector3d( 2, 2, 2 ).array();
    const pcaT::vectorT source = Eigen::Vector3d( 0.5, 1, 1 ).array();
    const pcaT::vectorT expected = Eigen::Vector3d( 0, 1, 1 ).array();
    constexpr double amplitude = 1e-10;
    resultT positive, negative;
    pcaT::workspaceT doubleWorkspace;
    mx::improc::P4PCA::calculate( positive, predictors, target + amplitude * source, { 1 }, 0, doubleWorkspace );
    mx::improc::P4PCA::calculate( negative, predictors, target - amplitude * source, { 1 }, 0, doubleWorkspace );
    const pcaT::vectorT measured = ( positive.residuals.col( 0 ) - negative.residuals.col( 0 ) ) / ( 2 * amplitude );
    REQUIRE( ( measured - expected ).matrix().norm() < 2e-6 );

    mx::improc::detail::P4PCAMixedWorkspace mixedWorkspace;
    mx::improc::detail::p4PCACalculateMixed( positive,
                                             predictors,
                                             target + amplitude * source,
                                             { 1 },
                                             0,
                                             mixedWorkspace );
    mx::improc::detail::p4PCACalculateMixed( negative,
                                             predictors,
                                             target - amplitude * source,
                                             { 1 },
                                             0,
                                             mixedWorkspace );
    REQUIRE( ( positive.residuals - negative.residuals ).matrix().norm() == 0 );
}

/** Verify mx::improc::P4PCA::calculate() needs smaller source steps near a separated truncation boundary. */
TEST_CASE( "Direct P4 response convergence depends on the cutoff eigengap", "[P4PCA][response][eigengap]" )
{
    constexpr double leading = 2.001;
    const pcaT::matrixT predictors = Eigen::Vector3d( leading, 2, 1 ).asDiagonal().toDenseMatrix().array();
    pcaT::matrixT source = pcaT::matrixT::Zero( 3, 3 );
    source( 1, 0 ) = 1;
    const pcaT::vectorT target = Eigen::Vector3d( 1, 1, 1 ).array();
    const double coupling = leading / ( leading * leading - 4 );
    const pcaT::vectorT expected = Eigen::Vector3d( -coupling, -coupling, 0 ).array();
    pcaT::workspaceT workspace;
    double previousError = 0;
    for( const double amplitude : { 0.01, 4e-5, 2e-5, 1e-5 } )
    {
        resultT positive, negative;
        mx::improc::P4PCA::calculate( positive, predictors + amplitude * source, target, { 1 }, 0, workspace );
        mx::improc::P4PCA::calculate( negative, predictors - amplitude * source, target, { 1 }, 0, workspace );
        const pcaT::vectorT measured =
            ( positive.residuals.col( 0 ) - negative.residuals.col( 0 ) ) / ( 2 * amplitude );
        const double error = ( measured - expected ).matrix().norm() / expected.matrix().norm();
        CAPTURE( amplitude, error );
        if( amplitude == 0.01 )
        {
            REQUIRE( error > 0.5 );
        }
        else
        {
            REQUIRE( error < 1e-3 );
            if( amplitude < 4e-5 )
            {
                REQUIRE( error < 0.26 * previousError );
                REQUIRE( error > 0.24 * previousError );
            }
        }
        previousError = error;
    }
}

/** Verify a retained/discarded degeneracy makes mx::improc::P4PCA::calculate() discontinuous along a crossing. */
TEST_CASE( "Direct P4 cutoff crossing has no finite derivative", "[P4PCA][response][cutoff]" )
{
    const pcaT::matrixT predictors = 2 * Eigen::Matrix2d::Identity().array();
    const pcaT::matrixT source = Eigen::Vector2d( 1, -1 ).asDiagonal().toDenseMatrix().array();
    const pcaT::vectorT target = Eigen::Vector2d( 1, 1 ).array();
    pcaT::workspaceT workspace;
    double previousNorm = 0;
    for( const double amplitude : { 0.01, 0.005, 0.0025 } )
    {
        resultT positive, negative;
        mx::improc::P4PCA::calculate( positive, predictors + amplitude * source, target, { 1 }, 0, workspace );
        mx::improc::P4PCA::calculate( negative, predictors - amplitude * source, target, { 1 }, 0, workspace );
        REQUIRE( positive.modeStatus == std::vector<statusT>{ statusT::rankSupported } );
        REQUIRE( negative.modeStatus == positive.modeStatus );
        const double norm = ( positive.residuals - negative.residuals ).matrix().norm() / ( 2 * amplitude );
        REQUIRE( norm == Approx( 1 / ( std::sqrt( 2.0 ) * amplitude ) ).epsilon( 1e-12 ) );
        if( previousNorm > 0 )
        {
            REQUIRE( norm == Approx( 2 * previousNorm ).epsilon( 1e-12 ) );
        }
        previousNorm = norm;
    }
}

/** Verify mx::improc::P4PCA::calculate() exposes asymmetric support when a probe crosses the rank threshold. */
TEST_CASE( "Direct P4 paired refits must check rank support on both sides", "[P4PCA][response][rank]" )
{
    const pcaT::matrixT predictors = Eigen::Vector2d( 2, 0.2 ).asDiagonal().toDenseMatrix().array();
    const pcaT::matrixT source = Eigen::Vector2d( 0, 1 ).asDiagonal().toDenseMatrix().array();
    const pcaT::vectorT target = Eigen::Vector2d( 1, 1 ).array();
    pcaT::workspaceT workspace;
    resultT positive, negative;
    mx::improc::P4PCA::calculate( positive, predictors + 0.001 * source, target, { 1, 2 }, 0.01, workspace );
    mx::improc::P4PCA::calculate( negative, predictors - 0.001 * source, target, { 1, 2 }, 0.01, workspace );
    REQUIRE( positive.numericalRank == 2 );
    REQUIRE( negative.numericalRank == 1 );
    REQUIRE( positive.modeStatus == std::vector<statusT>{ statusT::rankSupported, statusT::rankSupported } );
    REQUIRE( negative.modeStatus == std::vector<statusT>{ statusT::rankSupported, statusT::rankInsufficient } );
}

/** Verify P4PCA::calculateResponse() retains temporal nullspace motion when all predictor modes are retained.
 * A full temporal projector instead has identically zero residual derivative. Reusing the output must replace
 * prior dimensions, statuses and diagnostics, including when the next baseline has zero rank.
 */
TEST_CASE( "Analytic P4 response distinguishes full temporal and predictor rank", "[P4PCA][response][analytic]" )
{
    pcaT::workspaceT workspace;
    mx::improc::P4PCAResponseResult response;
    for( const bool temporalGram : { false, true } )
    {
        ResponseFixture fixture( temporalGram ? 6 : 9, temporalGram ? 9 : 6, { 36, 25, 16, 9, 4, 1 } );
        mx::improc::P4PCA::calculateResponse( response,
                                              fixture.m_predictors,
                                              fixture.m_target,
                                              fixture.m_predictorSource,
                                              fixture.m_targetSource,
                                              { 6 },
                                              1e-12,
                                              workspace );
        REQUIRE( response.modeStatus == std::vector{ mx::improc::P4PCAResponseStatus::differentiable } );
        REQUIRE( response.numericalRank == 6 );
        REQUIRE( ( response.responses.col( 0 ) - fixture.derivative( 6 ) ).matrix().norm() < 1e-12 );
        if( temporalGram )
        {
            REQUIRE( response.responses.abs().maxCoeff() == 0 );
            REQUIRE( std::isinf( response.relativeCutoffGaps[0] ) );
        }
        else
        {
            REQUIRE( response.responses.matrix().norm() > 0.1 );
            REQUIRE( response.relativeCutoffGaps[0] == Approx( 1.0 / 36 ) );
        }
    }
    mx::improc::P4PCA::calculateResponse( response,
                                          pcaT::matrixT::Zero( 3, 2 ),
                                          pcaT::vectorT::Ones( 3 ),
                                          pcaT::matrixT::Ones( 3, 2 ),
                                          pcaT::vectorT::Ones( 3 ),
                                          { 1, 2 },
                                          0,
                                          workspace );
    REQUIRE( response.numericalRank == 0 );
    REQUIRE( response.responses.rows() == 3 );
    REQUIRE( response.responses.cols() == 2 );
    REQUIRE( response.responses.isNaN().all() );
    REQUIRE( response.modeStatus == std::vector( 2, mx::improc::P4PCAResponseStatus::rankInsufficient ) );
}

/** Verify P4PCA::calculateResponse() separates rank failure, rank boundaries and unresolved cutoff gaps.
 * Degeneracy inside a retained/discarded group does not invalidate other requested counts. A caller may increase
 * the resolution floor; unresolved columns must remain NaN while usable columns are returned.
 */
TEST_CASE( "Analytic P4 response reports unresolved boundaries per mode", "[P4PCA][response][analytic][boundary]" )
{
    using responseStatus = mx::improc::P4PCAResponseStatus;
    pcaT::workspaceT workspace;
    mx::improc::P4PCAResponseResult response;
    for( const bool temporalGram : { true, false } )
    {
        ResponseFixture fixture( temporalGram ? 6 : 9, temporalGram ? 9 : 6, { 25, 25, 9, 9, 0, 0 } );
        mx::improc::P4PCA::calculateResponse( response,
                                              fixture.m_predictors,
                                              fixture.m_target,
                                              fixture.m_predictorSource,
                                              fixture.m_targetSource,
                                              { 1, 2, 3, 4, 5 },
                                              1e-12,
                                              workspace );
        REQUIRE( response.modeStatus == std::vector{ responseStatus::cutoffUnresolved,
                                                     responseStatus::differentiable,
                                                     responseStatus::cutoffUnresolved,
                                                     responseStatus::differentiable,
                                                     responseStatus::rankInsufficient } );
        for( const int index : { 0, 2, 4 } )
        {
            REQUIRE( response.responses.col( index ).isNaN().all() );
        }
        for( const int index : { 1, 3 } )
        {
            REQUIRE( ( response.responses.col( index ) - fixture.derivative( index + 1 ) ).matrix().norm() < 1e-12 );
        }
    }
    const pcaT::matrixT predictors = Eigen::Vector3d( 2, 0.5, 0.25 ).asDiagonal().toDenseMatrix().array();
    mx::improc::P4PCA::calculateResponse( response,
                                          predictors,
                                          pcaT::vectorT::Ones( 3 ),
                                          predictors,
                                          pcaT::vectorT::Ones( 3 ),
                                          { 1, 2, 3 },
                                          0.0625 - 1e-15,
                                          workspace );
    REQUIRE( response.modeStatus == std::vector{ responseStatus::differentiable,
                                                 responseStatus::rankBoundary,
                                                 responseStatus::rankInsufficient } );
    REQUIRE( response.responses.rightCols( 2 ).isNaN().all() );

    const pcaT::matrixT narrow = Eigen::Vector3d( 2.001, 2, 1 ).asDiagonal().toDenseMatrix().array();
    pcaT::matrixT source = pcaT::matrixT::Zero( 3, 3 );
    source( 1, 0 ) = 1;
    mx::improc::P4PCA::calculateResponse( response,
                                          narrow,
                                          pcaT::vectorT::Ones( 3 ),
                                          source,
                                          pcaT::vectorT::Zero( 3 ),
                                          { 1, 2 },
                                          0,
                                          workspace );
    const double coupling = 2.001 / ( 2.001 * 2.001 - 4 );
    REQUIRE( response.responses( 0, 0 ) == Approx( -coupling ).epsilon( 1e-12 ) );
    REQUIRE( response.responses( 1, 0 ) == Approx( -coupling ).epsilon( 1e-12 ) );
    REQUIRE( response.modeStatus == std::vector( 2, responseStatus::differentiable ) );
    mx::improc::P4PCA::calculateResponse( response,
                                          narrow,
                                          pcaT::vectorT::Ones( 3 ),
                                          source,
                                          pcaT::vectorT::Zero( 3 ),
                                          { 1, 2 },
                                          0,
                                          workspace,
                                          0.01 );
    REQUIRE( response.relativeResolution == 0.01 );
    REQUIRE( response.modeStatus == std::vector{ responseStatus::cutoffUnresolved, responseStatus::differentiable } );
    REQUIRE( response.responses.col( 0 ).isNaN().all() );
}

/** Verify P4PCA::calculateResponse() rejects malformed/nonfinite source inputs and invalid numerical controls. */
TEST_CASE( "Analytic P4 response validates source directions and controls", "[P4PCA][response][analytic][validation]" )
{
    ResponseFixture fixture( 6, 9, { 36, 25, 16, 9, 4, 1 } );
    pcaT::workspaceT workspace;
    mx::improc::P4PCAResponseResult response;
    for( const double invalid :
         { -1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity() } )
    {
        REQUIRE_THROWS_AS( mx::improc::P4PCA::calculateResponse( response,
                                                                 fixture.m_predictors,
                                                                 fixture.m_target,
                                                                 fixture.m_predictorSource,
                                                                 fixture.m_targetSource,
                                                                 { 1 },
                                                                 0,
                                                                 workspace,
                                                                 invalid ),
                           std::invalid_argument );
    }
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculateResponse( response,
                                                             fixture.m_predictors,
                                                             fixture.m_target,
                                                             pcaT::matrixT::Zero( 5, 9 ),
                                                             fixture.m_targetSource,
                                                             { 1 },
                                                             0,
                                                             workspace ),
                       std::invalid_argument );
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculateResponse( response,
                                                             fixture.m_predictors,
                                                             fixture.m_target,
                                                             fixture.m_predictorSource,
                                                             pcaT::vectorT::Zero( 5 ),
                                                             { 1 },
                                                             0,
                                                             workspace ),
                       std::invalid_argument );
    fixture.m_predictorSource( 0, 0 ) = std::numeric_limits<double>::quiet_NaN();
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculateResponse( response,
                                                             fixture.m_predictors,
                                                             fixture.m_target,
                                                             fixture.m_predictorSource,
                                                             fixture.m_targetSource,
                                                             { 1 },
                                                             0,
                                                             workspace ),
                       std::invalid_argument );
    fixture.m_predictorSource.setZero();
    fixture.m_targetSource( 0 ) = std::numeric_limits<double>::infinity();
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculateResponse( response,
                                                             fixture.m_predictors,
                                                             fixture.m_target,
                                                             fixture.m_predictorSource,
                                                             fixture.m_targetSource,
                                                             { 1 },
                                                             0,
                                                             workspace ),
                       std::invalid_argument );
    fixture.m_targetSource.setZero();
    for( const auto &modes : { std::vector<int>{}, { 0 }, { 2, 1 }, { 1, 1 }, { 7 } } )
    {
        REQUIRE_THROWS_AS( mx::improc::P4PCA::calculateResponse( response,
                                                                 fixture.m_predictors,
                                                                 fixture.m_target,
                                                                 fixture.m_predictorSource,
                                                                 fixture.m_targetSource,
                                                                 modes,
                                                                 0,
                                                                 workspace ),
                           std::invalid_argument );
    }
    response.responses = fixture.m_predictors;
    REQUIRE_THROWS_AS( mx::improc::P4PCA::calculateResponse( response,
                                                             response.responses,
                                                             fixture.m_target,
                                                             fixture.m_predictorSource,
                                                             fixture.m_targetSource,
                                                             { 1 },
                                                             0,
                                                             workspace ),
                       std::invalid_argument );
}

} // namespace P4PCAResponse_test
} // namespace unitTest
