/** \file P4NegativeOptimizer.hpp
 * \brief Declares contrast-only negative-companion optimization for pixel-local P4 reductions.
 * \author Jared R. Males
 */

#ifndef P4NegativeOptimizer_hpp
#define P4NegativeOptimizer_hpp

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

#include "P4Reduction.hpp"

namespace mx
{
namespace improc
{

/// Configuration for one deterministic bounded contrast-only P4 optimization.
/** \ingroup programming_library */
struct P4ContrastOptimizerConfig
{
    double contrastLower{ -0.05 };       ///< Inclusive lower signed contrast bound.

    double contrastUpper{ 0 };           ///< Inclusive upper signed contrast bound.

    std::size_t maxEvaluations{ 64 };    ///< Maximum number of distinct local P4 evaluations.

    std::size_t validationSamples{ 21 }; ///< Uniform dense-grid samples spanning both contrast bounds.

    double parameterTolerance{ 1e-5 };   ///< Absolute contrast tolerance for bounded minimization.

    double meritTolerance{ 1e-6 };       ///< Relative merit tolerance used for dense-basin agreement.
};

/// One evaluated contrast and its aperture merit.
/** \ingroup programming_library */
struct P4ContrastMeritSample
{
    double contrast{ 0 };          ///< Signed trial contrast.

    double merit{ 0 };             ///< Uniform aperture mean-square residual.

    bool denseValidation{ false }; ///< Whether the sample belongs to the required dense validation grid.

    double elapsedSeconds{ 0 };    ///< Total elapsed time for the local P4 evaluation.
};

/// Complete result of one contrast-only P4 optimization.
/** \ingroup programming_library */
struct P4ContrastOptimizationResult
{
    bool converged{ false };              ///< Whether the bounded search met its contrast-width tolerance.

    bool denseAgreement{ false };         ///< Whether the result agrees with the sampled dense-grid basin.

    double bestContrast{ 0 };             ///< Best finite evaluated signed contrast.

    double bestMerit{ 0 };                ///< Merit at `bestContrast`.

    std::size_t evaluationCount{ 0 };     ///< Number of distinct local P4 calls.

    double evaluationElapsedSeconds{ 0 }; ///< Sum of elapsed time across distinct local evaluations.

    ReductionTiming timing;               ///< Sum of detailed algorithm timings across distinct local evaluations.

    std::vector<P4ContrastMeritSample> samples; ///< Distinct evaluations in execution order.

    P4LocalEvaluation<float> bestEvaluation;    ///< Owning local products for the best evaluated contrast.

    std::string status;                         ///< Stable human-readable convergence status.
};

/// Calculate uniform mean-square residual merit inside a continuous-coordinate circular aperture.
/** Every pixel center inside the aperture must have a finite residual and nonzero validity. The aperture itself must
 * remain inside the local stamp's pixel-edge footprint.
 *
 * \returns the finite mean-square residual merit
 * \throws std::invalid_argument for invalid dimensions, mode, aperture, support, validity, or numeric values
 */
double p4LocalL2Merit( const P4LocalEvaluation<float> &evaluation, /**< [in] owning local P4 result */
                       std::size_t modeIndex,                      /**< [in] selected output-plane index */
                       double apertureRadius /**< [in] positive aperture radius in pixels */ );

/// Minimize one bounded signed contrast and validate it against a dense scan.
/** The dense grid is evaluated first. A Brent search then refines the sampled basin while retaining the best finite
 * evaluation. Cached contrasts never consume the evaluation budget twice.
 *
 * \returns the best evaluation and convergence diagnostics
 * \throws std::invalid_argument for an invalid optimizer contract
 * \throws std::runtime_error when the evaluation budget cannot complete the required dense scan
 */
P4ContrastOptimizationResult optimizeP4Contrast(
    const P4ContrastOptimizerConfig &configuration, /**< [in] validated search bounds and stopping controls */
    const std::function<P4LocalEvaluation<float>( double )> &evaluate,
    /**< [in] callback evaluating one signed contrast */
    std::size_t modeIndex, /**< [in] selected output-plane index */
    double apertureRadius /**< [in] positive merit-aperture radius in pixels */ );

} // namespace improc
} // namespace mx

#endif // P4NegativeOptimizer_hpp
