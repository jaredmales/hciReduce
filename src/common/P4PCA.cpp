/** \file P4PCA.cpp
 * \brief Implements the pure numerical predictor used by Pixel Prediction Post-Processing.
 * \author Jared R. Males
 */

#include "P4PCA.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <chrono>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

#include <mx/math/floatUtils.hpp>

namespace mx
{
namespace improc
{

namespace
{

/// Process-wide test override; production execution leaves this null.
std::atomic<detail::P4PCAEigenSolverT> p4PCAEigenSolverForTesting{ nullptr };

/// Determine whether every coefficient in an Eigen-like array is finite under fast-math.
template <typename arrayT>
bool p4PCAAllFinite( const arrayT &array /**< [in] array to validate */ )
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

/// Invoke the installed test solver or the production all-double mxlib eigensolver.
MXLAPACK_INT p4PCAEigenSolve( P4PCA::matrixT &eigenvectors, /**< [out] selected eigenvectors */
                              P4PCA::matrixT &eigenvalues,  /**< [out] selected eigenvalues */
                              P4PCA::matrixT &covariance,   /**< [in] lower-triangular covariance matrix */
                              int modeCount,                /**< [in] number of largest eigenpairs to calculate */
                              P4PCA::workspaceT &workspace /**< [in,out] reusable LAPACK workspace */ )
{
    const detail::P4PCAEigenSolverT testSolver = p4PCAEigenSolverForTesting.load( std::memory_order_acquire );
    if( testSolver )
    {
        return testSolver( eigenvectors, eigenvalues, covariance, modeCount, workspace );
    }

    const int dimension = static_cast<int>( covariance.rows() );
    return mx::math::eigenSYEVR( eigenvectors,
                                 eigenvalues,
                                 covariance,
                                 dimension - modeCount,
                                 dimension,
                                 'L',
                                 &workspace );
}

/// Validate one P4PCA request against its path-specific structural degree-of-freedom limit.
void p4PCAValidateInputs( const P4PCA::matrixT &predictors, /**< [in] predictor matrix to validate */
                          const P4PCA::vectorT &target,     /**< [in] target vector to validate */
                          const std::vector<int> &modes,    /**< [in] requested retained-mode counts */
                          double rankTolerance,             /**< [in] relative numerical-rank threshold */
                          int maxDOF /**< [in] path-specific structural degree-of-freedom limit */ )
{
    const Eigen::Index sampleCount = predictors.rows();
    const Eigen::Index predictorCount = predictors.cols();

    if( sampleCount <= 0 || predictorCount <= 0 )
    {
        throw std::invalid_argument( "P4PCA requires a nonempty predictor matrix" );
    }

    if( target.rows() != sampleCount )
    {
        throw std::invalid_argument( "P4PCA target length must match the predictor row count" );
    }

    if( modes.empty() )
    {
        throw std::invalid_argument( "P4PCA requires at least one mode count" );
    }

    int previousMode{ 0 };
    for( const int mode : modes )
    {
        if( mode <= 0 )
        {
            throw std::invalid_argument( "P4PCA mode counts must be positive" );
        }

        if( mode <= previousMode )
        {
            throw std::invalid_argument( "P4PCA mode counts must be strictly increasing" );
        }

        if( mode > maxDOF )
        {
            throw std::invalid_argument( "P4PCA mode count exceeds maxDOF" );
        }

        previousMode = mode;
    }

    if( !mx::math::isFinite( rankTolerance ) || rankTolerance < 0 )
    {
        throw std::invalid_argument( "P4PCA rank tolerance must be finite and nonnegative" );
    }

    if( !p4PCAAllFinite( predictors ) || !p4PCAAllFinite( target ) )
    {
        throw std::invalid_argument( "P4PCA predictors and target must be finite" );
    }
}

/// Materialize the sorted complement of one validated deletion set.
void p4PCARetainedRows( std::vector<Eigen::Index> &retainedRows, /**< [out] retained temporal row indices */
                        Eigen::Index sampleCount,                /**< [in] complete temporal row count */
                        std::span<const Eigen::Index> deletedRows /**< [in] sorted rows to omit */ )
{
    retainedRows.clear();
    retainedRows.reserve( static_cast<std::size_t>( sampleCount - static_cast<Eigen::Index>( deletedRows.size() ) ) );
    auto deleted = deletedRows.begin();
    for( Eigen::Index row = 0; row < sampleCount; ++row )
    {
        if( deleted != deletedRows.end() && *deleted == row )
        {
            ++deleted;
            continue;
        }
        retainedRows.push_back( row );
    }
}

/// Subtract each column's double-precision temporal mean in place and optionally return those means.
template <typename arrayT>
void p4PCACenterColumns( arrayT &array, /**< [in,out] array whose columns are centered */
                         P4PCA::matrixT *columnMeans = nullptr /**< [out] optional one-row array of removed means */ )
{
    if( columnMeans )
    {
        columnMeans->resize( 1, array.cols() );
        columnMeans->setZero();
    }

    for( Eigen::Index column = 0; column < array.cols(); ++column )
    {
        double scale{ 0 };
        for( Eigen::Index row = 0; row < array.rows(); ++row )
        {
            scale = std::max( scale, std::abs( array( row, column ) ) );
        }

        if( scale == 0 )
        {
            continue;
        }

        double scaledSum{ 0 };
        for( Eigen::Index row = 0; row < array.rows(); ++row )
        {
            scaledSum += array( row, column ) / scale;
        }

        const double scaledMean = std::clamp( scaledSum / static_cast<double>( array.rows() ), -1.0, 1.0 );
        const double mean = scale * scaledMean;
        if( columnMeans )
        {
            ( *columnMeans )( 0, column ) = mean;
        }
        array.col( column ) -= mean;
    }
}

/// Calculate one already-validated P4PCA request using its selected fit and residual targets.
void p4PCACalculateValidated(
    P4PCAResult &output,                  /**< [out] completed regression result */
    const P4PCA::matrixT &predictors,     /**< [in] predictors used to form the fit */
    const P4PCA::vectorT &fitTarget,      /**< [in] target used to calculate projections */
    const P4PCA::vectorT &residualTarget, /**< [in] initial target for output residuals */
    const std::vector<int> &modes,        /**< [in] requested retained-mode counts */
    double rankTolerance,                 /**< [in] relative numerical-rank threshold */
    int maxDOF,                           /**< [in] structural degree-of-freedom limit and eigenpair count */
    const P4PCA::matrixT *predictorMeans, /**< [in] optional means for uncentered application */
    P4PCA::workspaceT &workspace,         /**< [in,out] reusable LAPACK workspace */
    P4PCATiming *timing,                  /**< [out] optional per-call worker timing */
    P4PCA::matrixT *coefficients /**< [out] optional predictor-space coefficient vectors */ )
{
    if( timing )
    {
        *timing = P4PCATiming{};
    }

    const auto secondsSince = []( const std::chrono::steady_clock::time_point &start )
    { return std::chrono::duration<double>( std::chrono::steady_clock::now() - start ).count(); };
    const Eigen::Index sampleCount = predictors.rows();
    const Eigen::Index predictorCount = predictors.cols();
    const bool useTemporalGram = sampleCount <= predictorCount;
    P4PCA::matrixT gram;
    P4PCA::vectorT crossProduct;
    const auto gramStart = std::chrono::steady_clock::now();
    if( useTemporalGram )
    {
        gram = ( predictors.matrix() * predictors.matrix().transpose() ).array();
    }
    else
    {
        gram = ( predictors.matrix().transpose() * predictors.matrix() ).array();
        crossProduct = ( predictors.matrix().transpose() * fitTarget.matrix() ).array();
    }
    if( timing )
    {
        timing->gramWorkerSeconds = secondsSince( gramStart );
    }

    if( !p4PCAAllFinite( gram ) || ( !useTemporalGram && !p4PCAAllFinite( crossProduct ) ) )
    {
        throw std::runtime_error( "P4PCA normal equations produced nonfinite values" );
    }

    P4PCA::vectorT predictorMeanProjection;
    if( predictorMeans )
    {
        if( useTemporalGram )
        {
            predictorMeanProjection = ( predictors.matrix() * predictorMeans->matrix().transpose() ).array();
        }
    }

    P4PCA::matrixT eigenvectors;
    P4PCA::matrixT eigenvalues;
    const auto eigensolveStart = std::chrono::steady_clock::now();
    const MXLAPACK_INT solverStatus = p4PCAEigenSolve( eigenvectors, eigenvalues, gram, maxDOF, workspace );
    if( solverStatus != 0 )
    {
        throw std::runtime_error( "P4PCA eigensolver failed with status " + std::to_string( solverStatus ) );
    }

    if( eigenvectors.rows() != gram.rows() || eigenvectors.cols() != maxDOF || eigenvalues.rows() < maxDOF ||
        eigenvalues.cols() < 1 )
    {
        throw std::runtime_error( "P4PCA eigensolver returned invalid dimensions" );
    }

    if( !p4PCAAllFinite( eigenvectors ) || !p4PCAAllFinite( eigenvalues.block( 0, 0, maxDOF, 1 ) ) )
    {
        throw std::runtime_error( "P4PCA eigensolver returned nonfinite results" );
    }

    for( int eigenIndex = 1; eigenIndex < maxDOF; ++eigenIndex )
    {
        if( eigenvalues( eigenIndex ) < eigenvalues( eigenIndex - 1 ) )
        {
            throw std::runtime_error( "P4PCA eigensolver returned unsorted eigenvalues" );
        }
    }

    const double largestEigenvalue = eigenvalues( maxDOF - 1 );
    int numericalRank{ 0 };
    if( largestEigenvalue > 0 )
    {
        const long double rankThreshold =
            static_cast<long double>( rankTolerance ) * static_cast<long double>( largestEigenvalue );
        for( int eigenIndex = 0; eigenIndex < maxDOF; ++eigenIndex )
        {
            if( static_cast<long double>( eigenvalues( eigenIndex ) ) > rankThreshold )
            {
                ++numericalRank;
            }
        }
    }
    if( timing )
    {
        timing->eigensolveWorkerSeconds = secondsSince( eigensolveStart );
    }

    output.residuals.resize( sampleCount, modes.size() );
    output.residuals.setConstant( std::numeric_limits<double>::quiet_NaN() );
    output.modeStatus.assign( modes.size(), P4PCAModeStatus::rankInsufficient );
    output.numericalRank = numericalRank;
    output.baseRank = 0;
    output.numericalRankCapped = false;
    output.downdateClampCount = 0;
    output.explicitFallbackCount = 0;
    if( coefficients )
    {
        coefficients->resize( predictorCount, modes.size() );
        coefficients->setConstant( std::numeric_limits<double>::quiet_NaN() );
    }

    const auto unsupported = std::upper_bound( modes.begin(), modes.end(), numericalRank );
    if( unsupported == modes.begin() )
    {
        return;
    }

    const int largestSupportedMode = *( unsupported - 1 );
    const auto projectionStart = std::chrono::steady_clock::now();
    P4PCA::vectorT residual = residualTarget;
    P4PCA::vectorT accumulatedCoefficients;
    P4PCA::vectorT coefficientUpdate;
    if( coefficients )
    {
        accumulatedCoefficients.resize( predictorCount );
        accumulatedCoefficients.setZero();
        if( useTemporalGram )
        {
            coefficientUpdate.resize( predictorCount );
        }
    }
    std::size_t outputIndex{ 0 };

    for( int retainedMode = 1; retainedMode <= largestSupportedMode; ++retainedMode )
    {
        const int eigenIndex = maxDOF - retainedMode;
        const double eigenvalue = eigenvalues( eigenIndex );
        const P4PCA::vectorT eigenvector = eigenvectors.col( eigenIndex );
        double coefficient{ 0 };
        double predictionMean{ 0 };
        P4PCA::vectorT projectedMode;
        if( useTemporalGram )
        {
            coefficient = eigenvector.matrix().dot( fitTarget.matrix() );
            projectedMode = eigenvector * coefficient;
            if( coefficients )
            {
                coefficientUpdate = ( predictors.matrix().transpose() * eigenvector.matrix() ).array();
                coefficientUpdate *= coefficient / eigenvalue;
                accumulatedCoefficients += coefficientUpdate;
            }
            if( predictorMeans )
            {
                predictionMean =
                    eigenvector.matrix().dot( predictorMeanProjection.matrix() ) * coefficient / eigenvalue;
            }
        }
        else
        {
            coefficient = eigenvector.matrix().dot( crossProduct.matrix() ) / eigenvalue;
            projectedMode = ( predictors.matrix() * eigenvector.matrix() ).array();
            projectedMode *= coefficient;
            if( coefficients )
            {
                accumulatedCoefficients += eigenvector * coefficient;
            }
            if( predictorMeans )
            {
                predictionMean = ( predictorMeans->matrix() * eigenvector.matrix() )( 0, 0 ) * coefficient;
            }
        }

        if( !mx::math::isFinite( coefficient ) || !mx::math::isFinite( predictionMean ) ||
            !p4PCAAllFinite( projectedMode ) || ( coefficients && !p4PCAAllFinite( accumulatedCoefficients ) ) )
        {
            throw std::runtime_error( "P4PCA mode projection produced nonfinite values" );
        }

        if( predictorMeans )
        {
            projectedMode += predictionMean;
        }

        residual -= projectedMode;
        if( !p4PCAAllFinite( residual ) )
        {
            throw std::runtime_error( "P4PCA residual calculation produced nonfinite values" );
        }

        if( outputIndex < modes.size() && modes[outputIndex] == retainedMode )
        {
            output.residuals.col( outputIndex ) = residual;
            output.modeStatus[outputIndex] = P4PCAModeStatus::rankSupported;
            if( coefficients )
            {
                coefficients->col( outputIndex ) = accumulatedCoefficients;
            }
            ++outputIndex;
        }
    }
    if( timing )
    {
        timing->projectionWorkerSeconds = secondsSince( projectionStart );
    }
}

} // namespace

P4TargetExclusions P4TargetExclusions::fromSpans( Eigen::Index sampleCount, const std::vector<P4ExclusionSpan> &spans )
{
    if( sampleCount <= 0 || spans.size() != static_cast<std::size_t>( sampleCount ) )
    {
        throw std::invalid_argument( "P4 exclusion spans must provide one pattern per temporal row" );
    }

    P4TargetExclusions result;
    result.m_sampleCount = sampleCount;
    result.m_storage = P4ExclusionStorage::contiguousSpans;
    result.m_spans = spans;
    for( Eigen::Index target = 0; target < sampleCount; ++target )
    {
        const P4ExclusionSpan &span = spans[static_cast<std::size_t>( target )];
        if( span.begin < 0 || span.begin > target || span.end <= target || span.end > sampleCount )
        {
            throw std::invalid_argument( "P4 exclusion span must be in range, nonempty, and contain its target" );
        }
        result.m_maximumDeleted = std::max( result.m_maximumDeleted, span.end - span.begin );
    }
    return result;
}

P4TargetExclusions P4TargetExclusions::fromExplicit( Eigen::Index sampleCount,
                                                     const std::vector<std::vector<Eigen::Index>> &indices )
{
    if( sampleCount <= 0 || indices.size() != static_cast<std::size_t>( sampleCount ) )
    {
        throw std::invalid_argument( "P4 explicit exclusions must provide one pattern per temporal row" );
    }

    std::size_t totalIndices{ 0 };
    for( const auto &targetIndices : indices )
    {
        if( targetIndices.size() > std::numeric_limits<std::size_t>::max() - totalIndices )
        {
            throw std::length_error( "P4 explicit exclusion storage exceeds size_t range" );
        }
        totalIndices += targetIndices.size();
    }
    if( totalIndices > static_cast<std::size_t>( std::numeric_limits<Eigen::Index>::max() ) )
    {
        throw std::length_error( "P4 explicit exclusion storage exceeds Eigen index range" );
    }

    P4TargetExclusions result;
    result.m_sampleCount = sampleCount;
    result.m_storage = P4ExclusionStorage::explicitIndices;
    result.m_offsets.reserve( static_cast<std::size_t>( sampleCount ) + 1 );
    result.m_indices.reserve( totalIndices );
    result.m_offsets.push_back( 0 );
    for( Eigen::Index target = 0; target < sampleCount; ++target )
    {
        const auto &targetIndices = indices[static_cast<std::size_t>( target )];
        bool containsTarget{ false };
        Eigen::Index previous{ -1 };
        for( const Eigen::Index row : targetIndices )
        {
            if( row < 0 || row >= sampleCount || row <= previous )
            {
                throw std::invalid_argument( "P4 explicit exclusion rows must be sorted, unique, and in range" );
            }
            containsTarget = containsTarget || row == target;
            previous = row;
        }
        if( !containsTarget )
        {
            throw std::invalid_argument( "P4 explicit exclusion pattern must contain its target row" );
        }
        result.m_indices.insert( result.m_indices.end(), targetIndices.begin(), targetIndices.end() );
        result.m_offsets.push_back( static_cast<Eigen::Index>( result.m_indices.size() ) );
        result.m_maximumDeleted =
            std::max( result.m_maximumDeleted, static_cast<Eigen::Index>( targetIndices.size() ) );
    }
    return result;
}

bool P4TargetExclusions::empty() const noexcept
{
    return m_sampleCount == 0;
}

Eigen::Index P4TargetExclusions::sampleCount() const noexcept
{
    return m_sampleCount;
}

Eigen::Index P4TargetExclusions::targetCount() const noexcept
{
    return m_sampleCount;
}

P4ExclusionStorage P4TargetExclusions::storage() const noexcept
{
    return m_storage;
}

Eigen::Index P4TargetExclusions::deletedCount( Eigen::Index target ) const
{
    if( target < 0 || target >= m_sampleCount )
    {
        throw std::out_of_range( "P4 exclusion target is out of range" );
    }
    if( m_storage == P4ExclusionStorage::contiguousSpans )
    {
        const P4ExclusionSpan &span = m_spans[static_cast<std::size_t>( target )];
        return span.end - span.begin;
    }
    return m_offsets[static_cast<std::size_t>( target ) + 1] - m_offsets[static_cast<std::size_t>( target )];
}

Eigen::Index P4TargetExclusions::maximumDeleted() const noexcept
{
    return m_maximumDeleted;
}

std::span<const Eigen::Index> P4TargetExclusions::deletedRows( std::vector<Eigen::Index> &spanScratch,
                                                               Eigen::Index target ) const
{
    if( target < 0 || target >= m_sampleCount )
    {
        throw std::out_of_range( "P4 exclusion target is out of range" );
    }
    if( m_storage == P4ExclusionStorage::contiguousSpans )
    {
        const P4ExclusionSpan &span = m_spans[static_cast<std::size_t>( target )];
        spanScratch.resize( static_cast<std::size_t>( span.end - span.begin ) );
        std::iota( spanScratch.begin(), spanScratch.end(), span.begin );
        return spanScratch;
    }

    const Eigen::Index begin = m_offsets[static_cast<std::size_t>( target )];
    const Eigen::Index end = m_offsets[static_cast<std::size_t>( target ) + 1];
    return std::span<const Eigen::Index>( m_indices.data() + begin, static_cast<std::size_t>( end - begin ) );
}

std::size_t P4TargetExclusions::storageBytes() const noexcept
{
    const long double bytes = static_cast<long double>( m_spans.capacity() ) * sizeof( P4ExclusionSpan ) +
                              static_cast<long double>( m_offsets.capacity() ) * sizeof( Eigen::Index ) +
                              static_cast<long double>( m_indices.capacity() ) * sizeof( Eigen::Index );
    if( bytes > static_cast<long double>( std::numeric_limits<std::size_t>::max() ) )
    {
        return std::numeric_limits<std::size_t>::max();
    }
    return static_cast<std::size_t>( bytes );
}

namespace detail
{

/// \cond P4PCA_test_detail
void p4PCASetEigenSolverForTesting( P4PCAEigenSolverT solver )
{
    p4PCAEigenSolverForTesting.store( solver, std::memory_order_release );
}

void p4PCAResetEigenSolverForTesting()
{
    p4PCAEigenSolverForTesting.store( nullptr, std::memory_order_release );
}
/// \endcond

} // namespace detail

bool P4PCAResult::sampleSupported( Eigen::Index sample, std::size_t mode ) const
{
    if( sample < 0 || sample >= residuals.rows() || mode >= modeStatus.size() )
    {
        throw std::out_of_range( "P4PCA result sample or mode is out of range" );
    }
    if( sampleValidity.size() == 0 )
    {
        return modeStatus[mode] == P4PCAModeStatus::rankSupported;
    }
    if( sampleValidity.rows() != residuals.rows() || sampleValidity.cols() != residuals.cols() )
    {
        throw std::logic_error( "P4PCA sample-validity dimensions do not match residuals" );
    }
    return sampleValidity( sample, static_cast<Eigen::Index>( mode ) ) != 0;
}

void P4PCA::calculate( P4PCAResult &output,
                       const matrixT &predictors,
                       const vectorT &target,
                       const std::vector<int> &modes,
                       double rankTolerance,
                       workspaceT &workspace,
                       P4PCATiming *timing,
                       matrixT *coefficients )
{
    output.sampleValidity.resize( 0, 0 );
    const Eigen::Index sampleCount = predictors.rows();
    const Eigen::Index predictorCount = predictors.cols();
    const int maxDOF = static_cast<int>( std::min( sampleCount, predictorCount ) );
    p4PCAValidateInputs( predictors, target, modes, rankTolerance, maxDOF );
    p4PCACalculateValidated( output,
                             predictors,
                             target,
                             target,
                             modes,
                             rankTolerance,
                             maxDOF,
                             nullptr,
                             workspace,
                             timing,
                             coefficients );
}

void P4PCA::calculateHeldOut( P4PCAResult &output,
                              const matrixT &predictors,
                              const vectorT &target,
                              const P4TargetExclusions &exclusions,
                              const std::vector<int> &modes,
                              double rankTolerance,
                              workspaceT &workspace,
                              P4PCATiming *timing )
{
    const Eigen::Index sampleCount = predictors.rows();
    const Eigen::Index predictorCount = predictors.cols();
    if( sampleCount <= 0 || predictorCount <= 0 || target.rows() != sampleCount ||
        exclusions.sampleCount() != sampleCount || exclusions.targetCount() != sampleCount )
    {
        throw std::invalid_argument( "P4PCA held-out inputs must have consistent nonempty dimensions" );
    }
    if( modes.empty() )
    {
        throw std::invalid_argument( "P4PCA held-out regression requires at least one mode count" );
    }
    int previousMode{ 0 };
    for( const int mode : modes )
    {
        if( mode <= previousMode )
        {
            throw std::invalid_argument( "P4PCA held-out mode counts must be positive and strictly increasing" );
        }
        previousMode = mode;
    }
    if( !mx::math::isFinite( rankTolerance ) || rankTolerance < 0 || !p4PCAAllFinite( predictors ) ||
        !p4PCAAllFinite( target ) )
    {
        throw std::invalid_argument( "P4PCA held-out values and rank tolerance must be finite and valid" );
    }

    output.residuals.resize( sampleCount, static_cast<Eigen::Index>( modes.size() ) );
    output.residuals.setConstant( std::numeric_limits<double>::quiet_NaN() );
    output.sampleValidity.resize( sampleCount, static_cast<Eigen::Index>( modes.size() ) );
    output.sampleValidity.setZero();
    output.modeStatus.assign( modes.size(), P4PCAModeStatus::rankInsufficient );
    output.numericalRank = std::numeric_limits<int>::max();
    output.baseRank = 0;
    output.numericalRankCapped = false;
    output.downdateClampCount = 0;
    output.explicitFallbackCount = 0;
    if( timing )
    {
        *timing = P4PCATiming{};
    }

    std::vector<Eigen::Index> deletedRows;
    std::vector<Eigen::Index> retainedRows;
    if( sampleCount <= predictorCount )
    {
        const auto secondsSince = []( const std::chrono::steady_clock::time_point &start )
        { return std::chrono::duration<double>( std::chrono::steady_clock::now() - start ).count(); };
        const auto fullGramStart = std::chrono::steady_clock::now();
        const matrixT fullGram = ( predictors.matrix() * predictors.matrix().transpose() ).array();
        if( !p4PCAAllFinite( fullGram ) )
        {
            throw std::runtime_error( "P4PCA normal equations produced nonfinite values" );
        }
        if( timing )
        {
            timing->gramWorkerSeconds = secondsSince( fullGramStart );
        }

        matrixT trainingGram;
        vectorT trainingTarget;
        vectorT heldOutCrossProduct;
        matrixT eigenvectors;
        matrixT eigenvalues;
        for( Eigen::Index heldOut = 0; heldOut < sampleCount; ++heldOut )
        {
            const std::span<const Eigen::Index> excluded = exclusions.deletedRows( deletedRows, heldOut );
            p4PCARetainedRows( retainedRows, sampleCount, excluded );
            if( retainedRows.empty() )
            {
                output.numericalRank = 0;
                continue;
            }

            const int structuralRank = static_cast<int>( retainedRows.size() );
            output.numericalRank = std::min( output.numericalRank, structuralRank );
            const auto supportedEnd = std::upper_bound( modes.begin(), modes.end(), structuralRank );
            if( supportedEnd == modes.begin() )
            {
                continue;
            }

            const auto extractStart = std::chrono::steady_clock::now();
            trainingGram.resize( structuralRank, structuralRank );
            trainingTarget.resize( structuralRank );
            heldOutCrossProduct.resize( structuralRank );
            for( Eigen::Index trainingColumn = 0; trainingColumn < structuralRank; ++trainingColumn )
            {
                const Eigen::Index sourceColumn = retainedRows[static_cast<std::size_t>( trainingColumn )];
                trainingTarget( trainingColumn ) = target( sourceColumn );
                heldOutCrossProduct( trainingColumn ) = fullGram( heldOut, sourceColumn );
                for( Eigen::Index trainingRow = 0; trainingRow < structuralRank; ++trainingRow )
                {
                    const Eigen::Index sourceRow = retainedRows[static_cast<std::size_t>( trainingRow )];
                    trainingGram( trainingRow, trainingColumn ) = fullGram( sourceRow, sourceColumn );
                }
            }
            if( timing )
            {
                timing->gramWorkerSeconds += secondsSince( extractStart );
            }

            const auto eigensolveStart = std::chrono::steady_clock::now();
            const MXLAPACK_INT solverStatus =
                p4PCAEigenSolve( eigenvectors, eigenvalues, trainingGram, structuralRank, workspace );
            if( solverStatus != 0 )
            {
                throw std::runtime_error( "P4PCA eigensolver failed with status " + std::to_string( solverStatus ) );
            }
            if( eigenvectors.rows() != structuralRank || eigenvectors.cols() != structuralRank ||
                eigenvalues.rows() < structuralRank || eigenvalues.cols() < 1 )
            {
                throw std::runtime_error( "P4PCA eigensolver returned invalid dimensions" );
            }
            if( !p4PCAAllFinite( eigenvectors ) || !p4PCAAllFinite( eigenvalues.block( 0, 0, structuralRank, 1 ) ) )
            {
                throw std::runtime_error( "P4PCA eigensolver returned nonfinite results" );
            }
            for( int eigenIndex = 1; eigenIndex < structuralRank; ++eigenIndex )
            {
                if( eigenvalues( eigenIndex ) < eigenvalues( eigenIndex - 1 ) )
                {
                    throw std::runtime_error( "P4PCA eigensolver returned unsorted eigenvalues" );
                }
            }

            const double largestEigenvalue = eigenvalues( structuralRank - 1 );
            int numericalRank{ 0 };
            if( largestEigenvalue > 0 )
            {
                const long double rankThreshold =
                    static_cast<long double>( rankTolerance ) * static_cast<long double>( largestEigenvalue );
                for( int eigenIndex = 0; eigenIndex < structuralRank; ++eigenIndex )
                {
                    if( static_cast<long double>( eigenvalues( eigenIndex ) ) > rankThreshold )
                    {
                        ++numericalRank;
                    }
                }
            }
            output.numericalRank = std::min( output.numericalRank, numericalRank );
            if( timing )
            {
                timing->eigensolveWorkerSeconds += secondsSince( eigensolveStart );
            }

            const auto numericalEnd = std::upper_bound( modes.begin(), supportedEnd, numericalRank );
            if( numericalEnd == modes.begin() )
            {
                continue;
            }
            const int largestSupportedMode = *( numericalEnd - 1 );
            const auto projectionStart = std::chrono::steady_clock::now();
            double prediction{ 0 };
            std::size_t outputIndex{ 0 };
            for( int retainedMode = 1; retainedMode <= largestSupportedMode; ++retainedMode )
            {
                const int eigenIndex = structuralRank - retainedMode;
                const auto eigenvector = eigenvectors.col( eigenIndex );
                const double heldOutProjection = heldOutCrossProduct.matrix().dot( eigenvector.matrix() );
                const double targetProjection = eigenvector.matrix().dot( trainingTarget.matrix() );
                prediction += heldOutProjection * ( targetProjection / eigenvalues( eigenIndex ) );
                if( !mx::math::isFinite( prediction ) )
                {
                    throw std::runtime_error( "P4PCA held-out prediction produced a nonfinite residual" );
                }
                if( outputIndex < modes.size() && modes[outputIndex] == retainedMode )
                {
                    const double residual = target( heldOut ) - prediction;
                    if( !mx::math::isFinite( residual ) )
                    {
                        throw std::runtime_error( "P4PCA held-out prediction produced a nonfinite residual" );
                    }
                    output.residuals( heldOut, static_cast<Eigen::Index>( outputIndex ) ) = residual;
                    output.sampleValidity( heldOut, static_cast<Eigen::Index>( outputIndex ) ) = 1;
                    output.modeStatus[outputIndex] = P4PCAModeStatus::rankSupported;
                    ++outputIndex;
                }
            }
            if( timing )
            {
                timing->projectionWorkerSeconds += secondsSince( projectionStart );
            }
        }
        return;
    }

    matrixT trainingPredictors;
    vectorT trainingTarget;
    matrixT coefficients;
    for( Eigen::Index heldOut = 0; heldOut < sampleCount; ++heldOut )
    {
        const std::span<const Eigen::Index> excluded = exclusions.deletedRows( deletedRows, heldOut );
        p4PCARetainedRows( retainedRows, sampleCount, excluded );
        if( retainedRows.empty() )
        {
            output.numericalRank = 0;
            continue;
        }

        const int structuralRank = static_cast<int>(
            std::min<std::size_t>( retainedRows.size(), static_cast<std::size_t>( predictorCount ) ) );
        output.numericalRank = std::min( output.numericalRank, structuralRank );
        const auto supportedEnd = std::upper_bound( modes.begin(), modes.end(), structuralRank );
        if( supportedEnd == modes.begin() )
        {
            continue;
        }
        const std::vector<int> supportedModes( modes.begin(), supportedEnd );
        trainingPredictors.resize( static_cast<Eigen::Index>( retainedRows.size() ), predictorCount );
        trainingTarget.resize( static_cast<Eigen::Index>( retainedRows.size() ) );
        for( std::size_t training = 0; training < retainedRows.size(); ++training )
        {
            const Eigen::Index source = retainedRows[training];
            trainingPredictors.row( static_cast<Eigen::Index>( training ) ) = predictors.row( source );
            trainingTarget( static_cast<Eigen::Index>( training ) ) = target( source );
        }

        P4PCAResult fit;
        P4PCATiming fitTiming;
        calculate( fit,
                   trainingPredictors,
                   trainingTarget,
                   supportedModes,
                   rankTolerance,
                   workspace,
                   &fitTiming,
                   &coefficients );
        output.numericalRank = std::min( output.numericalRank, fit.numericalRank );
        if( timing )
        {
            timing->gramWorkerSeconds += fitTiming.gramWorkerSeconds;
            timing->eigensolveWorkerSeconds += fitTiming.eigensolveWorkerSeconds;
            timing->projectionWorkerSeconds += fitTiming.projectionWorkerSeconds;
        }
        const auto heldOutProjectionBegin = std::chrono::steady_clock::now();
        for( std::size_t mode = 0; mode < supportedModes.size(); ++mode )
        {
            if( fit.modeStatus[mode] == P4PCAModeStatus::rankInsufficient )
            {
                continue;
            }
            const double prediction = predictors.row( heldOut ).matrix().dot(
                coefficients.col( static_cast<Eigen::Index>( mode ) ).matrix() );
            const double residual = target( heldOut ) - prediction;
            if( !mx::math::isFinite( residual ) )
            {
                throw std::runtime_error( "P4PCA held-out prediction produced a nonfinite residual" );
            }
            output.residuals( heldOut, static_cast<Eigen::Index>( mode ) ) = residual;
            output.sampleValidity( heldOut, static_cast<Eigen::Index>( mode ) ) = 1;
            output.modeStatus[mode] = P4PCAModeStatus::rankSupported;
        }
        if( timing )
        {
            timing->projectionWorkerSeconds +=
                std::chrono::duration<double>( std::chrono::steady_clock::now() - heldOutProjectionBegin ).count();
        }
    }
}

void P4PCA::calculateHeldOutDowndated( P4PCAResult &output,
                                       const matrixT &predictors,
                                       const vectorT &target,
                                       const P4TargetExclusions &exclusions,
                                       const std::vector<int> &modes,
                                       double rankTolerance,
                                       workspaceT &eigensolverWorkspace,
                                       P4PCADowndateWorkspace &downdateWorkspace,
                                       P4PCATiming *timing )
{
    const Eigen::Index sampleCount = predictors.rows();
    const Eigen::Index predictorCount = predictors.cols();
    const Eigen::Index structuralBaseRank = std::min( sampleCount, predictorCount );
    if( structuralBaseRank > static_cast<Eigen::Index>( std::numeric_limits<int>::max() ) )
    {
        throw std::invalid_argument( "P4PCA downdate base rank exceeds the supported integer range" );
    }
    p4PCAValidateInputs( predictors, target, modes, rankTolerance, static_cast<int>( structuralBaseRank ) );
    if( exclusions.sampleCount() != sampleCount || exclusions.targetCount() != sampleCount )
    {
        throw std::invalid_argument( "P4PCA downdate exclusions must match the predictor row count" );
    }

    output.residuals.resize( sampleCount, static_cast<Eigen::Index>( modes.size() ) );
    output.residuals.setConstant( std::numeric_limits<double>::quiet_NaN() );
    output.sampleValidity.resize( sampleCount, static_cast<Eigen::Index>( modes.size() ) );
    output.sampleValidity.setZero();
    output.modeStatus.assign( modes.size(), P4PCAModeStatus::rankInsufficient );
    output.numericalRank = std::numeric_limits<int>::max();
    output.baseRank = 0;
    output.numericalRankCapped = false;
    output.downdateClampCount = 0;
    output.explicitFallbackCount = 0;
    if( timing )
    {
        *timing = P4PCATiming{};
    }

    Eigen::Index maximumActiveDeleted{ 0 };
    bool hasRetainedTarget{ false };
    for( Eigen::Index targetIndex = 0; targetIndex < sampleCount; ++targetIndex )
    {
        const Eigen::Index deletedCount = exclusions.deletedCount( targetIndex );
        if( deletedCount < sampleCount )
        {
            hasRetainedTarget = true;
            maximumActiveDeleted = std::max( maximumActiveDeleted, deletedCount );
        }
    }
    if( !hasRetainedTarget )
    {
        output.numericalRank = 0;
        return;
    }

    const auto secondsSince = []( const std::chrono::steady_clock::time_point &start )
    { return std::chrono::duration<double>( std::chrono::steady_clock::now() - start ).count(); };
    const bool temporalGram = sampleCount <= predictorCount;
    const auto gramStart = std::chrono::steady_clock::now();
    if( temporalGram )
    {
        downdateWorkspace.m_baseGram = ( predictors.matrix() * predictors.matrix().transpose() ).array();
    }
    else
    {
        downdateWorkspace.m_baseGram = ( predictors.matrix().transpose() * predictors.matrix() ).array();
    }
    if( !p4PCAAllFinite( downdateWorkspace.m_baseGram ) )
    {
        throw std::runtime_error( "P4PCA downdate base Gram matrix contains nonfinite values" );
    }
    if( timing )
    {
        timing->gramWorkerSeconds = secondsSince( gramStart );
    }

    const auto baseSolveStart = std::chrono::steady_clock::now();
    const int gramDimension = static_cast<int>( downdateWorkspace.m_baseGram.rows() );
    const MXLAPACK_INT baseStatus = p4PCAEigenSolve( downdateWorkspace.m_baseEigenvectors,
                                                     downdateWorkspace.m_baseEigenvalues,
                                                     downdateWorkspace.m_baseGram,
                                                     gramDimension,
                                                     eigensolverWorkspace );
    if( baseStatus != 0 )
    {
        throw std::runtime_error( "P4PCA downdate base eigensolver failed with status " +
                                  std::to_string( baseStatus ) );
    }
    if( downdateWorkspace.m_baseEigenvectors.rows() != gramDimension ||
        downdateWorkspace.m_baseEigenvectors.cols() != gramDimension ||
        downdateWorkspace.m_baseEigenvalues.rows() < gramDimension || downdateWorkspace.m_baseEigenvalues.cols() < 1 ||
        !p4PCAAllFinite( downdateWorkspace.m_baseEigenvectors ) ||
        !p4PCAAllFinite( downdateWorkspace.m_baseEigenvalues.block( 0, 0, gramDimension, 1 ) ) )
    {
        throw std::runtime_error( "P4PCA downdate base eigensolver returned an invalid singular system" );
    }
    for( int index = 1; index < gramDimension; ++index )
    {
        if( downdateWorkspace.m_baseEigenvalues( index ) < downdateWorkspace.m_baseEigenvalues( index - 1 ) )
        {
            throw std::runtime_error( "P4PCA downdate base eigensolver returned unsorted eigenvalues" );
        }
    }

    const double spectrumScale = std::max( std::abs( downdateWorkspace.m_baseEigenvalues( 0 ) ),
                                           std::abs( downdateWorkspace.m_baseEigenvalues( gramDimension - 1 ) ) );
    const double baseSafetyTolerance = 64.0 * std::numeric_limits<double>::epsilon() *
                                       static_cast<double>( std::max( 1, gramDimension ) ) * spectrumScale;
    if( downdateWorkspace.m_baseEigenvalues( 0 ) < -baseSafetyTolerance )
    {
        throw std::runtime_error( "P4PCA downdate base Gram matrix is materially non-positive-semidefinite" );
    }

    Eigen::Index baseRank{ 0 };
    for( Eigen::Index source = 0; source < structuralBaseRank; ++source )
    {
        const double eigenvalue = downdateWorkspace.m_baseEigenvalues( source );
        if( eigenvalue < -baseSafetyTolerance )
        {
            throw std::runtime_error( "P4PCA downdate base Gram matrix has a negative eigenvalue" );
        }
        if( eigenvalue > baseSafetyTolerance )
        {
            ++baseRank;
        }
    }

    if( temporalGram )
    {
        downdateWorkspace.m_temporalFactor.resize( sampleCount, baseRank );
        downdateWorkspace.m_singularValues.resize( baseRank );
        for( Eigen::Index mode = 0; mode < baseRank; ++mode )
        {
            const Eigen::Index source = structuralBaseRank - 1 - mode;
            const double eigenvalue = downdateWorkspace.m_baseEigenvalues( source );
            downdateWorkspace.m_singularValues( mode ) = std::sqrt( eigenvalue );
            downdateWorkspace.m_temporalFactor.col( mode ) = downdateWorkspace.m_baseEigenvectors.col( source );
        }
    }
    else
    {
        downdateWorkspace.m_temporalFactor.resize( sampleCount, baseRank );
        downdateWorkspace.m_singularValues.resize( baseRank );
        for( Eigen::Index mode = 0; mode < baseRank; ++mode )
        {
            const Eigen::Index source = structuralBaseRank - 1 - mode;
            const double singularValue = std::sqrt( downdateWorkspace.m_baseEigenvalues( source ) );
            downdateWorkspace.m_singularValues( mode ) = singularValue;
            downdateWorkspace.m_temporalFactor.col( mode ).matrix().noalias() =
                predictors.matrix() * downdateWorkspace.m_baseEigenvectors.col( source ).matrix();
            downdateWorkspace.m_temporalFactor.col( mode ) /= singularValue;
        }
    }
    output.baseRank = static_cast<int>( baseRank );
    if( baseRank == 0 )
    {
        output.numericalRank = 0;
        if( timing )
        {
            timing->eigensolveWorkerSeconds = secondsSince( baseSolveStart );
        }
        return;
    }

    const mx::math::svdDeletionStatus factorStatus =
        mx::math::validateSvdDeletionFactor( downdateWorkspace.m_temporalFactor );
    if( !mx::math::svdDeletionSucceeded( factorStatus ) )
    {
        throw std::runtime_error( "P4PCA downdate temporal factor validation failed with status " +
                                  std::string( mx::math::svdDeletionStatusName( factorStatus ) ) );
    }

    const Eigen::Index maximumOutputRank = std::min<Eigen::Index>( modes.back(), baseRank );
    const mx::math::svdDeletionStatus resultPrepareStatus =
        downdateWorkspace.m_deletionResult.prepare( baseRank, maximumOutputRank );
    const mx::math::svdDeletionStatus workspacePrepareStatus =
        downdateWorkspace.m_deletionWorkspace.prepare( baseRank,
                                                       maximumActiveDeleted,
                                                       mx::math::svdDeletionBackend::leadingCovariance );
    if( !mx::math::svdDeletionSucceeded( resultPrepareStatus ) ||
        !mx::math::svdDeletionSucceeded( workspacePrepareStatus ) )
    {
        const mx::math::svdDeletionStatus status =
            !mx::math::svdDeletionSucceeded( resultPrepareStatus ) ? resultPrepareStatus : workspacePrepareStatus;
        throw std::runtime_error( "P4PCA downdate workspace preparation failed with status " +
                                  std::string( mx::math::svdDeletionStatusName( status ) ) );
    }

    downdateWorkspace.m_baseTargetProjection =
        ( downdateWorkspace.m_temporalFactor.matrix().transpose() * target.matrix() ).array();
    if( !p4PCAAllFinite( downdateWorkspace.m_baseTargetProjection ) )
    {
        throw std::runtime_error( "P4PCA downdate base target projection contains nonfinite values" );
    }
    if( timing )
    {
        timing->eigensolveWorkerSeconds = secondsSince( baseSolveStart );
    }

    bool requiresExplicitOracle{ false };
    for( Eigen::Index heldOut = 0; heldOut < sampleCount; ++heldOut )
    {
        const std::span<const Eigen::Index> deleted =
            exclusions.deletedRows( downdateWorkspace.m_deletedRows, heldOut );
        const Eigen::Index retainedCount = sampleCount - static_cast<Eigen::Index>( deleted.size() );
        if( retainedCount == 0 )
        {
            output.numericalRank = 0;
            continue;
        }

        const int structuralRank = static_cast<int>( std::min<Eigen::Index>( baseRank, retainedCount ) );
        output.numericalRank = std::min( output.numericalRank, structuralRank );
        const auto structuralEnd = std::upper_bound( modes.begin(), modes.end(), structuralRank );
        if( structuralEnd == modes.begin() )
        {
            continue;
        }

        const Eigen::Index requestedOutputRank = *( structuralEnd - 1 );
        const auto deletionStart = std::chrono::steady_clock::now();
        const mx::math::svdDeletionStatus deletionStatus =
            mx::math::svdRemoveRows<double>( downdateWorkspace.m_deletionResult,
                                             downdateWorkspace.m_singularValues,
                                             downdateWorkspace.m_temporalFactor,
                                             deleted,
                                             requestedOutputRank,
                                             downdateWorkspace.m_deletionWorkspace,
                                             mx::math::svdDeletionBackend::leadingCovariance );
        if( !mx::math::svdDeletionSucceeded( deletionStatus ) )
        {
            throw std::runtime_error( "P4PCA row deletion failed for target " + std::to_string( heldOut ) +
                                      " with status " + mx::math::svdDeletionStatusName( deletionStatus ) +
                                      " and LAPACK status " +
                                      std::to_string( downdateWorkspace.m_deletionResult.lapackInfo() ) );
        }
        output.downdateClampCount +=
            static_cast<std::size_t>( downdateWorkspace.m_deletionResult.clampedEigenvalues() );
        if( timing )
        {
            timing->eigensolveWorkerSeconds += secondsSince( deletionStart );
        }

        const auto squaredSingularValues = downdateWorkspace.m_deletionResult.squaredSingularValues();
        int numericalRank{ 0 };
        const double largestEigenvalue = squaredSingularValues( 0 );
        const double rankThreshold = rankTolerance * largestEigenvalue;
        const double ambiguityScale =
            std::max( std::abs( downdateWorkspace.m_baseEigenvalues( structuralBaseRank - 1 ) ),
                      std::abs( largestEigenvalue ) );
        const double ambiguityTolerance = 256.0 * std::numeric_limits<double>::epsilon() *
                                          static_cast<double>( std::max<Eigen::Index>( 1, baseRank ) ) * ambiguityScale;
        if( largestEigenvalue > 0 )
        {
            for( int mode = 0; mode < structuralRank; ++mode )
            {
                if( std::abs( squaredSingularValues( mode ) - rankThreshold ) <= ambiguityTolerance )
                {
                    requiresExplicitOracle = true;
                    break;
                }
                if( squaredSingularValues( mode ) > rankThreshold )
                {
                    ++numericalRank;
                }
            }
        }
        if( requiresExplicitOracle )
        {
            break;
        }
        output.numericalRank = std::min( output.numericalRank, numericalRank );
        const auto numericalEnd = std::upper_bound( modes.begin(), structuralEnd, numericalRank );
        if( numericalEnd == modes.begin() )
        {
            continue;
        }

        const auto projectionStart = std::chrono::steady_clock::now();
        downdateWorkspace.m_retainedTargetProjection = downdateWorkspace.m_baseTargetProjection;
        for( const Eigen::Index row : deleted )
        {
            downdateWorkspace.m_retainedTargetProjection -=
                downdateWorkspace.m_temporalFactor.row( row ).transpose() * target( row );
        }
        downdateWorkspace.m_scaledTargetProjection =
            downdateWorkspace.m_singularValues * downdateWorkspace.m_retainedTargetProjection;
        downdateWorkspace.m_scaledHeldOutRow =
            downdateWorkspace.m_singularValues * downdateWorkspace.m_temporalFactor.row( heldOut ).transpose();
        const auto rotation = downdateWorkspace.m_deletionResult.rotation();
        double prediction{ 0 };
        std::size_t outputIndex{ 0 };
        const int largestSupportedMode = *( numericalEnd - 1 );
        for( int retainedMode = 1; retainedMode <= largestSupportedMode; ++retainedMode )
        {
            const Eigen::Index mode = retainedMode - 1;
            const double targetCoordinate =
                rotation.col( mode ).matrix().dot( downdateWorkspace.m_scaledTargetProjection.matrix() );
            const double heldOutCoordinate =
                rotation.col( mode ).matrix().dot( downdateWorkspace.m_scaledHeldOutRow.matrix() );
            prediction += heldOutCoordinate * targetCoordinate / squaredSingularValues( mode );
            if( !mx::math::isFinite( prediction ) )
            {
                throw std::runtime_error( "P4PCA factor-space prediction produced a nonfinite value" );
            }
            if( outputIndex < modes.size() && modes[outputIndex] == retainedMode )
            {
                const double residual = target( heldOut ) - prediction;
                if( !mx::math::isFinite( residual ) )
                {
                    throw std::runtime_error( "P4PCA factor-space residual produced a nonfinite value" );
                }
                output.residuals( heldOut, static_cast<Eigen::Index>( outputIndex ) ) = residual;
                output.sampleValidity( heldOut, static_cast<Eigen::Index>( outputIndex ) ) = 1;
                output.modeStatus[outputIndex] = P4PCAModeStatus::rankSupported;
                ++outputIndex;
            }
        }
        if( timing )
        {
            timing->projectionWorkerSeconds += secondsSince( projectionStart );
        }
    }

    if( requiresExplicitOracle )
    {
        const int attemptedBaseRank = output.baseRank;
        const std::size_t attemptedClampCount = output.downdateClampCount;
        P4PCAResult explicitResult;
        P4PCATiming explicitTiming;
        calculateHeldOut( explicitResult,
                          predictors,
                          target,
                          exclusions,
                          modes,
                          rankTolerance,
                          eigensolverWorkspace,
                          &explicitTiming );
        explicitResult.baseRank = attemptedBaseRank;
        explicitResult.downdateClampCount = attemptedClampCount;
        explicitResult.explicitFallbackCount = static_cast<std::size_t>( sampleCount );
        output = std::move( explicitResult );
        if( timing )
        {
            timing->gramWorkerSeconds += explicitTiming.gramWorkerSeconds;
            timing->eigensolveWorkerSeconds += explicitTiming.eigensolveWorkerSeconds;
            timing->projectionWorkerSeconds += explicitTiming.projectionWorkerSeconds;
        }
    }
}

void P4PCA::calculateCentered( P4PCAResult &output,
                               const matrixT &predictors,
                               const vectorT &target,
                               const std::vector<int> &modes,
                               double rankTolerance,
                               workspaceT &workspace,
                               P4PCATiming *timing,
                               matrixT *coefficients )
{
    matrixT centeredPredictors = predictors;
    calculateCenteredInPlace( output,
                              centeredPredictors,
                              target,
                              modes,
                              rankTolerance,
                              workspace,
                              timing,
                              coefficients );
}

void P4PCA::calculateCenteredInPlace( P4PCAResult &output,
                                      matrixT &predictors,
                                      const vectorT &target,
                                      const std::vector<int> &modes,
                                      double rankTolerance,
                                      workspaceT &workspace,
                                      P4PCATiming *timing,
                                      matrixT *coefficients )
{
    output.sampleValidity.resize( 0, 0 );
    const Eigen::Index sampleCount = predictors.rows();
    const Eigen::Index predictorCount = predictors.cols();
    if( sampleCount == 1 )
    {
        throw std::invalid_argument( "P4PCA centered regression requires at least two samples" );
    }

    const int maxDOF = static_cast<int>( std::min( predictorCount, sampleCount - 1 ) );
    p4PCAValidateInputs( predictors, target, modes, rankTolerance, maxDOF );

    vectorT centeredTarget = target;
    matrixT predictorMeans;
    p4PCACenterColumns( predictors, &predictorMeans );
    p4PCACenterColumns( centeredTarget );
    if( !p4PCAAllFinite( predictors ) || !p4PCAAllFinite( centeredTarget ) )
    {
        throw std::runtime_error( "P4PCA temporal centering produced nonfinite values" );
    }

    p4PCACalculateValidated( output,
                             predictors,
                             centeredTarget,
                             target,
                             modes,
                             rankTolerance,
                             maxDOF,
                             &predictorMeans,
                             workspace,
                             timing,
                             coefficients );
}

} // namespace improc
} // namespace mx
