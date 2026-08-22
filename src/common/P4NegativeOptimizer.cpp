/** \file P4NegativeOptimizer.cpp
 * \brief Implements contrast-only negative-companion optimization for pixel-local P4 reductions.
 * \author Jared R. Males
 */

#include "P4NegativeOptimizer.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>

#include <mx/math/floatUtils.hpp>

namespace mx
{
namespace improc
{

namespace
{

/// One cached local evaluation and its scalar merit.
struct P4CachedContrastEvaluation
{
    double merit{ 0 };                   ///< Calculated aperture merit.

    P4LocalEvaluation<float> evaluation; ///< Owning local P4 result.
};

/// Report whether a double remains finite under fast-math.
bool p4OptimizerFinite( double value /**< [in] scalar to validate */ ) noexcept
{
    return mx::math::isFinite( value );
}

/// Add one evaluation's detailed timing to an optimizer-session total.
void p4AccumulateTiming( ReductionTiming &total, /**< [in,out] accumulated session timing */
                         const ReductionTiming &evaluation /**< [in] one local evaluation timing */ )
{
    total.geometryElapsedSeconds += evaluation.geometryElapsedSeconds;
    total.regressionElapsedSeconds += evaluation.regressionElapsedSeconds;
    total.samplingWorkerSeconds += evaluation.samplingWorkerSeconds;
    total.sameImageSamplingWorkerSeconds += evaluation.sameImageSamplingWorkerSeconds;
    total.temporalSamplingWorkerSeconds += evaluation.temporalSamplingWorkerSeconds;
    total.gramWorkerSeconds += evaluation.gramWorkerSeconds;
    total.eigensolveWorkerSeconds += evaluation.eigensolveWorkerSeconds;
    total.modeWorkerSeconds += evaluation.modeWorkerSeconds;
    total.projectionWorkerSeconds += evaluation.projectionWorkerSeconds;
    total.psfWorkerSeconds += evaluation.psfWorkerSeconds;
    total.psfReconstructionElapsedSeconds += evaluation.psfReconstructionElapsedSeconds;
}

} // namespace

double p4LocalL2Merit( const P4LocalEvaluation<float> &evaluation, std::size_t modeIndex, double apertureRadius )
{
    if( evaluation.residual.rows() <= 0 || evaluation.residual.cols() <= 0 || evaluation.residual.planes() <= 0 ||
        evaluation.validity.rows() != evaluation.residual.rows() ||
        evaluation.validity.cols() != evaluation.residual.cols() ||
        evaluation.validity.planes() != evaluation.residual.planes() )
    {
        throw std::invalid_argument( "P4 local merit requires matching nonempty residual and validity cubes" );
    }
    if( modeIndex >= static_cast<std::size_t>( evaluation.residual.planes() ) )
    {
        throw std::invalid_argument( "P4 local merit mode index is outside the result cube" );
    }
    if( !p4OptimizerFinite( apertureRadius ) || apertureRadius <= 0 )
    {
        throw std::invalid_argument( "P4 local merit aperture radius must be finite and positive" );
    }

    const double sourceRow = evaluation.sourceRow - static_cast<double>( evaluation.originRow );
    const double sourceColumn = evaluation.sourceColumn - static_cast<double>( evaluation.originColumn );
    if( !p4OptimizerFinite( sourceRow ) || !p4OptimizerFinite( sourceColumn ) || sourceRow - apertureRadius < -0.5 ||
        sourceRow + apertureRadius > static_cast<double>( evaluation.residual.rows() ) - 0.5 ||
        sourceColumn - apertureRadius < -0.5 ||
        sourceColumn + apertureRadius > static_cast<double>( evaluation.residual.cols() ) - 0.5 )
    {
        throw std::invalid_argument( "P4 local merit aperture crosses the local stamp edge" );
    }

    const auto &residual = evaluation.residual.image( static_cast<int>( modeIndex ) );
    const auto &validity = evaluation.validity.image( static_cast<int>( modeIndex ) );
    const double radiusSquared = apertureRadius * apertureRadius;
    long double squareSum{ 0 };
    std::size_t sampleCount{ 0 };
    for( int column = 0; column < residual.cols(); ++column )
    {
        for( int row = 0; row < residual.rows(); ++row )
        {
            const double rowOffset = static_cast<double>( row ) - sourceRow;
            const double columnOffset = static_cast<double>( column ) - sourceColumn;
            if( rowOffset * rowOffset + columnOffset * columnOffset > radiusSquared )
            {
                continue;
            }
            const float value = residual( row, column );
            if( validity( row, column ) == 0 || !mx::math::isFinite( value ) )
            {
                throw std::invalid_argument( "P4 local merit aperture contains an invalid residual sample" );
            }
            squareSum += static_cast<long double>( value ) * static_cast<long double>( value );
            ++sampleCount;
        }
    }
    if( sampleCount == 0 )
    {
        throw std::invalid_argument( "P4 local merit aperture contains no pixel centers" );
    }
    const double merit = static_cast<double>( squareSum / static_cast<long double>( sampleCount ) );
    if( !p4OptimizerFinite( merit ) )
    {
        throw std::invalid_argument( "P4 local merit is non-finite" );
    }
    return merit;
}

P4ContrastOptimizationResult optimizeP4Contrast( const P4ContrastOptimizerConfig &configuration,
                                                 const std::function<P4LocalEvaluation<float>( double )> &evaluate,
                                                 std::size_t modeIndex,
                                                 double apertureRadius )
{
    if( !evaluate || !p4OptimizerFinite( configuration.contrastLower ) ||
        !p4OptimizerFinite( configuration.contrastUpper ) ||
        configuration.contrastLower >= configuration.contrastUpper || configuration.contrastUpper > 0 )
    {
        throw std::invalid_argument( "P4 contrast optimizer requires finite ordered nonpositive bounds" );
    }
    if( configuration.validationSamples < 3 || configuration.maxEvaluations < configuration.validationSamples + 3 )
    {
        throw std::invalid_argument( "P4 contrast optimizer budget cannot contain the dense scan and refinement" );
    }
    if( !p4OptimizerFinite( configuration.parameterTolerance ) || configuration.parameterTolerance <= 0 ||
        !p4OptimizerFinite( configuration.meritTolerance ) || configuration.meritTolerance < 0 ||
        !p4OptimizerFinite( apertureRadius ) || apertureRadius <= 0 )
    {
        throw std::invalid_argument( "P4 contrast optimizer tolerances and aperture are invalid" );
    }

    P4ContrastOptimizationResult result;
    result.bestMerit = std::numeric_limits<double>::infinity();
    std::map<double, P4CachedContrastEvaluation> cache;
    const auto sampleContrast = [&]( double contrast, bool dense ) -> const P4CachedContrastEvaluation &
    {
        const auto cached = cache.find( contrast );
        if( cached != cache.end() )
        {
            return cached->second;
        }
        if( cache.size() >= configuration.maxEvaluations )
        {
            throw std::runtime_error( "P4 contrast optimizer exhausted its evaluation budget" );
        }
        P4CachedContrastEvaluation entry;
        entry.evaluation = evaluate( contrast );
        entry.merit = p4LocalL2Merit( entry.evaluation, modeIndex, apertureRadius );
        auto inserted = cache.emplace( contrast, std::move( entry ) );
        result.samples.push_back(
            { contrast, inserted.first->second.merit, dense, inserted.first->second.evaluation.elapsedSeconds } );
        result.evaluationElapsedSeconds += inserted.first->second.evaluation.elapsedSeconds;
        p4AccumulateTiming( result.timing, inserted.first->second.evaluation.timing );
        if( inserted.first->second.merit < result.bestMerit )
        {
            result.bestContrast = contrast;
            result.bestMerit = inserted.first->second.merit;
            result.bestEvaluation = inserted.first->second.evaluation;
        }
        return inserted.first->second;
    };

    const double denseStep = ( configuration.contrastUpper - configuration.contrastLower ) /
                             static_cast<double>( configuration.validationSamples - 1 );
    std::size_t denseMinimum{ 0 };
    double denseMinimumMerit = std::numeric_limits<double>::infinity();
    std::vector<double> denseContrasts;
    denseContrasts.reserve( configuration.validationSamples );
    for( std::size_t index = 0; index < configuration.validationSamples; ++index )
    {
        const double contrast = index + 1 == configuration.validationSamples
                                    ? configuration.contrastUpper
                                    : configuration.contrastLower + static_cast<double>( index ) * denseStep;
        denseContrasts.push_back( contrast );
        const double merit = sampleContrast( contrast, true ).merit;
        if( merit < denseMinimumMerit )
        {
            denseMinimum = index;
            denseMinimumMerit = merit;
        }
    }

    double lower = denseMinimum == 0 ? denseContrasts[0] : denseContrasts[denseMinimum - 1];
    double upper = denseMinimum + 1 == denseContrasts.size() ? denseContrasts.back() : denseContrasts[denseMinimum + 1];
    double x = denseContrasts[denseMinimum];
    double w = x;
    double v = x;
    double fx = sampleContrast( x, false ).merit;
    double fw = fx;
    double fv = fx;
    double movement{ 0 };
    double previousMovement{ 0 };
    constexpr double goldenSection = 0.3819660112501051;
    constexpr double machineTolerance = 1e-12;

    while( cache.size() < configuration.maxEvaluations )
    {
        const double midpoint = 0.5 * ( lower + upper );
        const double tolerance = configuration.parameterTolerance + machineTolerance * std::abs( x );
        const double doubleTolerance = 2 * tolerance;
        if( std::abs( x - midpoint ) <= doubleTolerance - 0.5 * ( upper - lower ) )
        {
            result.converged = true;
            break;
        }

        bool parabolicStep{ false };
        if( std::abs( previousMovement ) > tolerance )
        {
            const double r = ( x - w ) * ( fx - fv );
            double q = ( x - v ) * ( fx - fw );
            double p = ( x - v ) * q - ( x - w ) * r;
            q = 2 * ( q - r );
            if( q > 0 )
            {
                p = -p;
            }
            q = std::abs( q );
            const double savedMovement = previousMovement;
            previousMovement = movement;
            if( q > 0 && std::abs( p ) < std::abs( 0.5 * q * savedMovement ) && p > q * ( lower - x ) &&
                p < q * ( upper - x ) )
            {
                movement = p / q;
                const double proposal = x + movement;
                if( proposal - lower < doubleTolerance || upper - proposal < doubleTolerance )
                {
                    movement = std::copysign( tolerance, midpoint - x );
                }
                parabolicStep = true;
            }
        }
        if( !parabolicStep )
        {
            previousMovement = x < midpoint ? upper - x : lower - x;
            movement = goldenSection * previousMovement;
        }

        const double candidate = x + ( std::abs( movement ) >= tolerance
                                           ? movement
                                           : std::copysign( tolerance, movement == 0 ? midpoint - x : movement ) );
        const double candidateMerit = sampleContrast( candidate, false ).merit;
        if( candidateMerit <= fx )
        {
            if( candidate >= x )
            {
                lower = x;
            }
            else
            {
                upper = x;
            }
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = candidate;
            fx = candidateMerit;
        }
        else
        {
            if( candidate < x )
            {
                lower = candidate;
            }
            else
            {
                upper = candidate;
            }
            if( candidateMerit <= fw || w == x )
            {
                v = w;
                fv = fw;
                w = candidate;
                fw = candidateMerit;
            }
            else if( candidateMerit <= fv || v == x || v == w )
            {
                v = candidate;
                fv = candidateMerit;
            }
        }
    }

    result.evaluationCount = cache.size();
    const double basinLower = denseContrasts[denseMinimum == 0 ? 0 : denseMinimum - 1];
    const double basinUpper =
        denseContrasts[denseMinimum + 1 == denseContrasts.size() ? denseMinimum : denseMinimum + 1];
    const double meritScale = std::max( std::abs( denseMinimumMerit ), std::numeric_limits<double>::min() );
    result.denseAgreement = result.bestContrast >= basinLower - configuration.parameterTolerance &&
                            result.bestContrast <= basinUpper + configuration.parameterTolerance &&
                            result.bestMerit <= denseMinimumMerit + configuration.meritTolerance * meritScale;
    if( result.converged && result.denseAgreement )
    {
        result.status = "converged";
    }
    else if( !result.converged )
    {
        result.status = "evaluation budget exhausted before contrast convergence";
    }
    else
    {
        result.status = "bounded minimum does not agree with dense validation basin";
    }
    return result;
}

} // namespace improc
} // namespace mx
