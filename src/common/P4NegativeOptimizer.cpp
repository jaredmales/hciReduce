/** \file P4NegativeOptimizer.cpp
 * \brief Implements negative-companion optimization for pixel-local P4 reductions.
 * \author Jared R. Males
 */

#include "P4NegativeOptimizer.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numbers>
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

/// One cached joint-optimizer vertex and its owning local result.
struct P4JointVertex
{
    std::array<double, 3> point{};       ///< Row delta, column delta, and contrast.

    P4LocalTrial trial;                  ///< Polar trial passed to the P4 evaluator.

    double merit{ 0 };                   ///< Calculated aperture merit.

    P4LocalEvaluation<float> evaluation; ///< Owning local P4 result.
};

/// Internal signal that preserves the best joint result when its reserved budget is exhausted.
class P4JointBudgetExhausted
{
};

constexpr std::size_t p4ContrastRefinementReserve = 16; ///< Brent evaluations reserved beyond each dense scan.

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

/// Reflect a proposed bounded-simplex coordinate back into its inclusive interval.
double p4ReflectIntoBounds( double value, /**< [in] proposed coordinate */
                            double lower, /**< [in] inclusive lower bound */
                            double upper /**< [in] inclusive upper bound */ )
{
    while( value < lower || value > upper )
    {
        if( value < lower )
        {
            value = lower + ( lower - value );
        }
        if( value > upper )
        {
            value = upper - ( value - upper );
        }
    }
    return std::clamp( value, lower, upper );
}

} // namespace

P4CartesianOffset p4TrialCartesianOffset( const P4LocalTrial &trial )
{
    if( !p4OptimizerFinite( trial.separation ) || trial.separation < 0 || !p4OptimizerFinite( trial.positionAngle ) )
    {
        throw std::invalid_argument( "P4 trial position must have finite nonnegative separation and finite PA" );
    }
    const double radians = -trial.positionAngle * std::numbers::pi / 180.0;
    return { trial.separation * std::sin( radians ), trial.separation * std::cos( radians ) };
}

P4LocalTrial p4CartesianOffsetTrial( const P4CartesianOffset &offset, double contrast )
{
    if( !p4OptimizerFinite( offset.row ) || !p4OptimizerFinite( offset.column ) || !p4OptimizerFinite( contrast ) )
    {
        throw std::invalid_argument( "P4 Cartesian trial offset and contrast must be finite" );
    }
    const double separation = std::hypot( offset.row, offset.column );
    if( separation == 0 )
    {
        return { 0, 0, contrast };
    }
    double positionAngle = -std::atan2( offset.row, offset.column ) * 180.0 / std::numbers::pi;
    positionAngle = std::fmod( positionAngle, 360.0 );
    if( positionAngle < 0 )
    {
        positionAngle += 360.0;
    }
    return { separation, positionAngle, contrast };
}

std::size_t p4PositionContrastMinimumEvaluations( std::size_t validationSamples )
{
    constexpr std::size_t simplexVertices = 4;
    if( validationSamples >
        ( std::numeric_limits<std::size_t>::max() - simplexVertices ) / 2 - p4ContrastRefinementReserve )
    {
        throw std::overflow_error( "P4 joint optimizer validation grid exceeds the evaluation-count range" );
    }
    return 2 * ( validationSamples + p4ContrastRefinementReserve ) + simplexVertices;
}

double p4LocalL2Merit( const P4LocalEvaluation<float> &evaluation,
                       std::size_t modeIndex,
                       double apertureRadius,
                       double apertureRow,
                       double apertureColumn )
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

    const double localApertureRow = apertureRow - static_cast<double>( evaluation.originRow );
    const double localApertureColumn = apertureColumn - static_cast<double>( evaluation.originColumn );
    if( !p4OptimizerFinite( localApertureRow ) || !p4OptimizerFinite( localApertureColumn ) ||
        localApertureRow - apertureRadius < -0.5 ||
        localApertureRow + apertureRadius > static_cast<double>( evaluation.residual.rows() ) - 0.5 ||
        localApertureColumn - apertureRadius < -0.5 ||
        localApertureColumn + apertureRadius > static_cast<double>( evaluation.residual.cols() ) - 0.5 )
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
            const double rowOffset = static_cast<double>( evaluation.originRow + row ) - apertureRow;
            const double columnOffset = static_cast<double>( evaluation.originColumn + column ) - apertureColumn;
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

double p4LocalL2Merit( const P4LocalEvaluation<float> &evaluation, std::size_t modeIndex, double apertureRadius )
{
    return p4LocalL2Merit( evaluation, modeIndex, apertureRadius, evaluation.sourceRow, evaluation.sourceColumn );
}

namespace
{

/// Implement contrast minimization with either evaluation-centered or fixed full-image aperture support.
P4ContrastOptimizationResult optimizeP4ContrastImpl( const P4ContrastOptimizerConfig &configuration,
                                                     const std::function<P4LocalEvaluation<float>( double )> &evaluate,
                                                     std::size_t modeIndex,
                                                     double apertureRadius,
                                                     bool fixedAperture,
                                                     double apertureRow,
                                                     double apertureColumn )
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
        entry.merit = fixedAperture
                          ? p4LocalL2Merit( entry.evaluation, modeIndex, apertureRadius, apertureRow, apertureColumn )
                          : p4LocalL2Merit( entry.evaluation, modeIndex, apertureRadius );
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
    const double meritScale = std::max( std::abs( denseMinimumMerit ), 1.0 );
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

} // namespace

P4ContrastOptimizationResult optimizeP4Contrast( const P4ContrastOptimizerConfig &configuration,
                                                 const std::function<P4LocalEvaluation<float>( double )> &evaluate,
                                                 std::size_t modeIndex,
                                                 double apertureRadius )
{
    return optimizeP4ContrastImpl( configuration, evaluate, modeIndex, apertureRadius, false, 0, 0 );
}

P4ContrastOptimizationResult optimizeP4Contrast( const P4ContrastOptimizerConfig &configuration,
                                                 const std::function<P4LocalEvaluation<float>( double )> &evaluate,
                                                 std::size_t modeIndex,
                                                 double apertureRadius,
                                                 double apertureRow,
                                                 double apertureColumn )
{
    if( !p4OptimizerFinite( apertureRow ) || !p4OptimizerFinite( apertureColumn ) )
    {
        throw std::invalid_argument( "P4 fixed merit-aperture coordinates must be finite" );
    }
    return optimizeP4ContrastImpl( configuration,
                                   evaluate,
                                   modeIndex,
                                   apertureRadius,
                                   true,
                                   apertureRow,
                                   apertureColumn );
}

P4PositionContrastOptimizationResult
optimizeP4PositionContrast( const P4ContrastOptimizerConfig &configuration,
                            const P4LocalTrial &initialTrial,
                            double positionBound,
                            const std::function<P4LocalEvaluation<float>( const P4LocalTrial & )> &evaluate,
                            std::size_t modeIndex,
                            double apertureRadius )
{
    const P4CartesianOffset initialOffset = p4TrialCartesianOffset( initialTrial );
    if( !evaluate || !p4OptimizerFinite( initialTrial.contrast ) ||
        initialTrial.contrast < configuration.contrastLower || initialTrial.contrast > configuration.contrastUpper ||
        !p4OptimizerFinite( positionBound ) || positionBound <= 0 ||
        !p4OptimizerFinite( configuration.positionTolerance ) || configuration.positionTolerance <= 0 ||
        configuration.maxEvaluations < p4PositionContrastMinimumEvaluations( configuration.validationSamples ) )
    {
        throw std::invalid_argument( "P4 joint optimizer requires a valid trial, positive position controls, and "
                                     "sufficient total evaluation budget" );
    }

    const std::size_t contrastStageBudget = configuration.validationSamples + p4ContrastRefinementReserve;
    P4ContrastOptimizerConfig contrastConfiguration = configuration;
    contrastConfiguration.maxEvaluations = contrastStageBudget;
    const auto evaluateInitialContrast = [&evaluate, &initialTrial]( double contrast )
    {
        P4LocalTrial trial = initialTrial;
        trial.contrast = contrast;
        return evaluate( trial );
    };
    const P4ContrastOptimizationResult seed =
        optimizeP4Contrast( contrastConfiguration, evaluateInitialContrast, modeIndex, apertureRadius );

    P4PositionContrastOptimizationResult result;
    result.bestTrial = initialTrial;
    result.bestTrial.contrast = seed.bestContrast;
    result.apertureRow = seed.bestEvaluation.sourceRow;
    result.apertureColumn = seed.bestEvaluation.sourceColumn;
    result.bestMerit = seed.bestMerit;
    result.evaluationCount = seed.evaluationCount;
    result.evaluationElapsedSeconds = seed.evaluationElapsedSeconds;
    result.timing = seed.timing;
    result.bestEvaluation = seed.bestEvaluation;
    result.samples.reserve( configuration.maxEvaluations );
    for( const P4ContrastMeritSample &sample : seed.samples )
    {
        P4LocalTrial trial = initialTrial;
        trial.contrast = sample.contrast;
        result.samples.push_back( { trial,
                                    0,
                                    0,
                                    sample.merit,
                                    sample.denseValidation ? "seed-dense" : "seed-brent",
                                    sample.elapsedSeconds } );
    }

    using pointT = std::array<double, 3>;
    const pointT lowerBounds{ -positionBound, -positionBound, configuration.contrastLower };
    const pointT upperBounds{ positionBound, positionBound, configuration.contrastUpper };
    const std::size_t jointEvaluationLimit = configuration.maxEvaluations - contrastStageBudget;
    std::map<pointT, P4JointVertex> cache;
    P4JointVertex seedVertex;
    seedVertex.point = { 0, 0, seed.bestContrast };
    seedVertex.trial = result.bestTrial;
    seedVertex.merit = seed.bestMerit;
    seedVertex.evaluation = seed.bestEvaluation;
    cache.emplace( seedVertex.point, seedVertex );

    const auto boundedPoint = [&lowerBounds, &upperBounds]( pointT point )
    {
        for( std::size_t coordinate = 0; coordinate < point.size(); ++coordinate )
        {
            point[coordinate] =
                p4ReflectIntoBounds( point[coordinate], lowerBounds[coordinate], upperBounds[coordinate] );
        }
        return point;
    };
    const auto samplePoint = [&]( pointT point ) -> const P4JointVertex &
    {
        point = boundedPoint( point );
        const auto cached = cache.find( point );
        if( cached != cache.end() )
        {
            return cached->second;
        }
        if( result.evaluationCount >= jointEvaluationLimit )
        {
            throw P4JointBudgetExhausted{};
        }

        P4JointVertex vertex;
        vertex.point = point;
        vertex.trial =
            p4CartesianOffsetTrial( { initialOffset.row + point[0], initialOffset.column + point[1] }, point[2] );
        vertex.evaluation = evaluate( vertex.trial );
        vertex.merit =
            p4LocalL2Merit( vertex.evaluation, modeIndex, apertureRadius, result.apertureRow, result.apertureColumn );
        auto inserted = cache.emplace( point, std::move( vertex ) );
        const P4JointVertex &stored = inserted.first->second;
        result.samples.push_back(
            { stored.trial, point[0], point[1], stored.merit, "simplex", stored.evaluation.elapsedSeconds } );
        ++result.evaluationCount;
        result.evaluationElapsedSeconds += stored.evaluation.elapsedSeconds;
        p4AccumulateTiming( result.timing, stored.evaluation.timing );
        if( stored.merit < result.bestMerit )
        {
            result.bestTrial = stored.trial;
            result.bestRowDelta = point[0];
            result.bestColumnDelta = point[1];
            result.bestMerit = stored.merit;
            result.bestEvaluation = stored.evaluation;
        }
        return stored;
    };

    const double positionStep =
        std::min( positionBound, std::max( 0.25 * positionBound, 16 * configuration.positionTolerance ) );
    const double contrastRange = configuration.contrastUpper - configuration.contrastLower;
    const double contrastStep =
        std::min( contrastRange, std::max( 0.1 * contrastRange, 16 * configuration.parameterTolerance ) );
    double contrastVertex = seed.bestContrast;
    if( seed.bestContrast - configuration.contrastLower >= configuration.contrastUpper - seed.bestContrast )
    {
        contrastVertex = std::max( configuration.contrastLower, seed.bestContrast - contrastStep );
    }
    else
    {
        contrastVertex = std::min( configuration.contrastUpper, seed.bestContrast + contrastStep );
    }

    std::array<P4JointVertex, 4> simplex;
    try
    {
        simplex[0] = seedVertex;
        simplex[1] = samplePoint( { positionStep, 0, seed.bestContrast } );
        simplex[2] = samplePoint( { 0, positionStep, seed.bestContrast } );
        simplex[3] = samplePoint( { 0, 0, contrastVertex } );

        const auto meritOrder = []( const P4JointVertex &left, const P4JointVertex &right )
        { return left.merit < right.merit; };
        while( result.evaluationCount < jointEvaluationLimit )
        {
            std::sort( simplex.begin(), simplex.end(), meritOrder );
            double positionDiameter{ 0 };
            double contrastDiameter{ 0 };
            for( std::size_t vertex = 1; vertex < simplex.size(); ++vertex )
            {
                positionDiameter = std::max( positionDiameter,
                                             std::hypot( simplex[vertex].point[0] - simplex[0].point[0],
                                                         simplex[vertex].point[1] - simplex[0].point[1] ) );
                contrastDiameter =
                    std::max( contrastDiameter, std::abs( simplex[vertex].point[2] - simplex[0].point[2] ) );
            }
            const double meritScale = std::max( std::abs( simplex[0].merit ), 1.0 );
            const double meritSpread = simplex.back().merit - simplex.front().merit;
            if( positionDiameter <= configuration.positionTolerance &&
                contrastDiameter <= configuration.parameterTolerance &&
                meritSpread <= configuration.meritTolerance * meritScale )
            {
                result.positionConverged = true;
                break;
            }

            pointT centroid{};
            for( std::size_t vertex = 0; vertex < 3; ++vertex )
            {
                for( std::size_t coordinate = 0; coordinate < centroid.size(); ++coordinate )
                {
                    centroid[coordinate] += simplex[vertex].point[coordinate] / 3.0;
                }
            }
            pointT reflected{};
            for( std::size_t coordinate = 0; coordinate < reflected.size(); ++coordinate )
            {
                reflected[coordinate] = centroid[coordinate] + centroid[coordinate] - simplex[3].point[coordinate];
            }
            const P4JointVertex reflectedVertex = samplePoint( reflected );
            if( reflectedVertex.merit < simplex[0].merit )
            {
                pointT expanded{};
                for( std::size_t coordinate = 0; coordinate < expanded.size(); ++coordinate )
                {
                    expanded[coordinate] =
                        centroid[coordinate] + 2 * ( reflectedVertex.point[coordinate] - centroid[coordinate] );
                }
                const P4JointVertex expandedVertex = samplePoint( expanded );
                simplex[3] = expandedVertex.merit < reflectedVertex.merit ? expandedVertex : reflectedVertex;
                continue;
            }
            if( reflectedVertex.merit < simplex[2].merit )
            {
                simplex[3] = reflectedVertex;
                continue;
            }

            const bool outsideContraction = reflectedVertex.merit < simplex[3].merit;
            pointT contracted{};
            for( std::size_t coordinate = 0; coordinate < contracted.size(); ++coordinate )
            {
                const double endpoint =
                    outsideContraction ? reflectedVertex.point[coordinate] : simplex[3].point[coordinate];
                contracted[coordinate] = centroid[coordinate] + 0.5 * ( endpoint - centroid[coordinate] );
            }
            const P4JointVertex contractedVertex = samplePoint( contracted );
            if( contractedVertex.merit < ( outsideContraction ? reflectedVertex.merit : simplex[3].merit ) )
            {
                simplex[3] = contractedVertex;
                continue;
            }

            for( std::size_t vertex = 1; vertex < simplex.size(); ++vertex )
            {
                pointT shrunk{};
                for( std::size_t coordinate = 0; coordinate < shrunk.size(); ++coordinate )
                {
                    shrunk[coordinate] = simplex[0].point[coordinate] +
                                         0.5 * ( simplex[vertex].point[coordinate] - simplex[0].point[coordinate] );
                }
                simplex[vertex] = samplePoint( shrunk );
            }
        }
    }
    catch( const P4JointBudgetExhausted & )
    {
    }

    const double fittedRowDelta = result.bestRowDelta;
    const double fittedColumnDelta = result.bestColumnDelta;
    const P4CartesianOffset fittedOffset{ initialOffset.row + fittedRowDelta,
                                          initialOffset.column + fittedColumnDelta };
    P4ContrastOptimizerConfig finalContrastConfiguration = configuration;
    finalContrastConfiguration.maxEvaluations = configuration.maxEvaluations - result.evaluationCount;
    const auto evaluateFinalContrast = [&evaluate, &fittedOffset]( double contrast )
    { return evaluate( p4CartesianOffsetTrial( fittedOffset, contrast ) ); };
    const P4ContrastOptimizationResult finalContrast = optimizeP4Contrast( finalContrastConfiguration,
                                                                           evaluateFinalContrast,
                                                                           modeIndex,
                                                                           apertureRadius,
                                                                           result.apertureRow,
                                                                           result.apertureColumn );
    for( const P4ContrastMeritSample &sample : finalContrast.samples )
    {
        const P4LocalTrial trial = p4CartesianOffsetTrial( fittedOffset, sample.contrast );
        result.samples.push_back( { trial,
                                    fittedRowDelta,
                                    fittedColumnDelta,
                                    sample.merit,
                                    sample.denseValidation ? "final-dense" : "final-brent",
                                    sample.elapsedSeconds } );
    }
    result.evaluationCount += finalContrast.evaluationCount;
    result.evaluationElapsedSeconds += finalContrast.evaluationElapsedSeconds;
    p4AccumulateTiming( result.timing, finalContrast.timing );
    result.contrastConverged = finalContrast.converged;
    result.denseAgreement = finalContrast.denseAgreement;
    result.bestTrial = p4CartesianOffsetTrial( fittedOffset, finalContrast.bestContrast );
    result.bestMerit = finalContrast.bestMerit;
    result.bestEvaluation = finalContrast.bestEvaluation;
    result.converged = result.positionConverged && result.contrastConverged && result.denseAgreement;
    if( result.converged )
    {
        result.status = "converged";
    }
    else if( !result.positionConverged )
    {
        result.status = "evaluation budget exhausted before position convergence";
    }
    else if( !result.contrastConverged )
    {
        result.status = "evaluation budget exhausted before final contrast convergence";
    }
    else
    {
        result.status = "final contrast does not agree with dense validation basin";
    }
    return result;
}

P4JackknifeStatistics p4JackknifeStatistics( const std::vector<P4JackknifeEstimate> &estimates, bool fitPosition )
{
    if( estimates.size() < 2 )
    {
        throw std::invalid_argument( "P4 jackknife covariance requires at least two block estimates" );
    }

    P4JackknifeStatistics statistics;
    statistics.blockCount = estimates.size();
    for( const P4JackknifeEstimate &estimate : estimates )
    {
        if( !p4OptimizerFinite( estimate.contrast ) ||
            ( fitPosition &&
              ( !p4OptimizerFinite( estimate.rowDelta ) || !p4OptimizerFinite( estimate.columnDelta ) ) ) )
        {
            throw std::invalid_argument( "P4 jackknife estimates must be finite" );
        }
        if( fitPosition )
        {
            statistics.mean.rowDelta += estimate.rowDelta;
            statistics.mean.columnDelta += estimate.columnDelta;
        }
        statistics.mean.contrast += estimate.contrast;
    }
    const double inverseBlockCount = 1.0 / static_cast<double>( estimates.size() );
    statistics.mean.rowDelta *= inverseBlockCount;
    statistics.mean.columnDelta *= inverseBlockCount;
    statistics.mean.contrast *= inverseBlockCount;

    for( const P4JackknifeEstimate &estimate : estimates )
    {
        const std::array<double, 3> delta{ fitPosition ? estimate.rowDelta - statistics.mean.rowDelta : 0,
                                           fitPosition ? estimate.columnDelta - statistics.mean.columnDelta : 0,
                                           estimate.contrast - statistics.mean.contrast };
        for( std::size_t row = 0; row < delta.size(); ++row )
        {
            for( std::size_t column = 0; column < delta.size(); ++column )
            {
                statistics.covariance[row * delta.size() + column] += delta[row] * delta[column];
            }
        }
    }
    const double multiplier = static_cast<double>( estimates.size() - 1 ) * inverseBlockCount;
    for( double &value : statistics.covariance )
    {
        value *= multiplier;
    }
    statistics.rowStandardError = std::sqrt( std::max( 0.0, statistics.covariance[0] ) );
    statistics.columnStandardError = std::sqrt( std::max( 0.0, statistics.covariance[4] ) );
    statistics.contrastStandardError = std::sqrt( std::max( 0.0, statistics.covariance[8] ) );
    return statistics;
}

} // namespace improc
} // namespace mx
