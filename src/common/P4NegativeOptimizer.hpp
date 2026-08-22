/** \file P4NegativeOptimizer.hpp
 * \brief Declares negative-companion optimization for pixel-local P4 reductions.
 * \author Jared R. Males
 */

#ifndef P4NegativeOptimizer_hpp
#define P4NegativeOptimizer_hpp

#include <array>
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

    double positionTolerance{ 1e-3 };    ///< Absolute Cartesian position tolerance in pixels for joint fitting.

    double meritTolerance{ 1e-6 };       ///< Combined absolute/relative merit convergence tolerance.
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

/// Cartesian sky offset from the `0.5*(N-1)` image center.
/** \ingroup programming_library */
struct P4CartesianOffset
{
    double row{ 0 };    ///< Signed source-row offset in pixels.

    double column{ 0 }; ///< Signed source-column offset in pixels.
};

/// One evaluated joint position-and-contrast trial and its aperture merit.
/** \ingroup programming_library */
struct P4PositionContrastMeritSample
{
    P4LocalTrial trial;         ///< Evaluated separation, position angle, and signed contrast.

    double rowDelta{ 0 };       ///< Cartesian row displacement from the configured initial trial.

    double columnDelta{ 0 };    ///< Cartesian column displacement from the configured initial trial.

    double merit{ 0 };          ///< Uniform aperture mean-square residual.

    std::string stage;          ///< Stable evaluation stage: seed, simplex, or final validation.

    double elapsedSeconds{ 0 }; ///< Total elapsed time for the local P4 evaluation.
};

/// Complete result of one bounded joint position-and-contrast P4 optimization.
/** \ingroup programming_library */
struct P4PositionContrastOptimizationResult
{
    bool converged{ false };              ///< Whether position and final contrast validation both converged.

    bool positionConverged{ false };      ///< Whether the bounded Cartesian simplex met its stopping criteria.

    bool contrastConverged{ false };      ///< Whether the final fixed-position contrast refinement converged.

    bool denseAgreement{ false };         ///< Whether the final contrast agrees with its dense validation basin.

    P4LocalTrial bestTrial;               ///< Best final separation, position angle, and signed contrast.

    double bestRowDelta{ 0 };             ///< Fitted Cartesian row displacement from the initial trial.

    double bestColumnDelta{ 0 };          ///< Fitted Cartesian column displacement from the initial trial.

    double apertureRow{ 0 };              ///< Fixed full-image merit-aperture row at the initial trial.

    double apertureColumn{ 0 };           ///< Fixed full-image merit-aperture column at the initial trial.

    double bestMerit{ 0 };                ///< Merit at `bestTrial`.

    std::size_t evaluationCount{ 0 };     ///< Number of distinct local P4 calls across all stages.

    double evaluationElapsedSeconds{ 0 }; ///< Sum of elapsed time across all local P4 calls.

    ReductionTiming timing;               ///< Sum of detailed algorithm timings across all local evaluations.

    std::vector<P4PositionContrastMeritSample> samples; ///< Evaluations in execution order.

    P4LocalEvaluation<float> bestEvaluation;            ///< Owning local products for the fitted trial.

    std::string status;                                 ///< Stable human-readable convergence status.
};

/// One delete-one-block fitted estimate in Cartesian position and signed contrast.
/** \ingroup programming_library */
struct P4JackknifeEstimate
{
    double rowDelta{ 0 };    ///< Fitted Cartesian row displacement from the common initial trial.

    double columnDelta{ 0 }; ///< Fitted Cartesian column displacement from the common initial trial.

    double contrast{ 0 };    ///< Fitted signed contrast.
};

/// Ordinary delete-one-block jackknife mean, covariance, and standard errors.
/** Covariance storage uses row-major `(rowDelta,columnDelta,contrast)` parameter order.
 * \ingroup programming_library
 */
struct P4JackknifeStatistics
{
    std::size_t blockCount{ 0 };        ///< Number of delete-one-block estimates in the calculation.

    P4JackknifeEstimate mean;           ///< Arithmetic mean of the delete-one-block estimates.

    std::array<double, 9> covariance{}; ///< Ordinary jackknife covariance in row-major parameter order.

    double rowStandardError{ 0 };       ///< Square root of the Cartesian-row covariance diagonal.

    double columnStandardError{ 0 };    ///< Square root of the Cartesian-column covariance diagonal.

    double contrastStandardError{ 0 };  ///< Square root of the contrast covariance diagonal.
};

/// Convert a separation and position angle to its Cartesian sky offset.
/** This is the exact inverse convention used by P4 local processing: row is `separation*sin(-PA)` and column is
 * `separation*cos(-PA)`.
 *
 * \returns the signed Cartesian source offset from image center
 * \throws std::invalid_argument for a non-finite or negative trial position
 */
P4CartesianOffset p4TrialCartesianOffset( const P4LocalTrial &trial /**< [in] polar local trial */ );

/// Convert one Cartesian sky offset and contrast to the P4 separation/PA convention.
/** Position angle is normalized to `[0,360)`; the zero-separation position angle is zero.
 *
 * \returns the equivalent polar local trial
 * \throws std::invalid_argument for a non-finite coordinate or contrast
 */
P4LocalTrial p4CartesianOffsetTrial( const P4CartesianOffset &offset, /**< [in] Cartesian source offset */
                                     double contrast /**< [in] signed trial contrast */ );

/// Return the minimum total evaluation budget accepted for joint position-and-contrast optimization.
/** The budget contains two dense+Brent contrast stages and four initial simplex vertices.
 *
 * \returns the minimum number of local P4 evaluations
 * \throws std::overflow_error when the requested validation grid cannot be represented safely
 */
std::size_t p4PositionContrastMinimumEvaluations(
    std::size_t validationSamples /**< [in] dense-grid samples used by each contrast stage */ );

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

/// Calculate uniform mean-square residual merit inside a fixed full-image aperture.
/** The explicit center keeps the sampled sky pixels invariant when the trial source moves. Every included pixel
 * center must be present in the local stamp and have a finite residual and nonzero validity.
 *
 * \returns the finite mean-square residual merit
 * \throws std::invalid_argument for invalid dimensions, mode, aperture, support, validity, or numeric values
 */
double p4LocalL2Merit( const P4LocalEvaluation<float> &evaluation, /**< [in] owning local P4 result */
                       std::size_t modeIndex,                      /**< [in] selected output-plane index */
                       double apertureRadius,                      /**< [in] positive aperture radius in pixels */
                       double apertureRow,                         /**< [in] fixed full-image aperture-center row */
                       double apertureColumn /**< [in] fixed full-image aperture-center column */ );

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

/// Minimize one bounded signed contrast using a fixed full-image merit aperture.
/** This variant preserves identical sky-pixel support across a surrounding joint position search.
 *
 * \returns the best evaluation and convergence diagnostics
 * \throws std::invalid_argument for an invalid optimizer or aperture contract
 * \throws std::runtime_error when the evaluation budget cannot complete the required dense scan
 */
P4ContrastOptimizationResult optimizeP4Contrast(
    const P4ContrastOptimizerConfig &configuration, /**< [in] validated search bounds and stopping controls */
    const std::function<P4LocalEvaluation<float>( double )> &evaluate,
    /**< [in] callback evaluating one signed contrast */
    std::size_t modeIndex, /**< [in] selected output-plane index */
    double apertureRadius, /**< [in] positive merit-aperture radius in pixels */
    double apertureRow,    /**< [in] fixed full-image aperture-center row */
    double apertureColumn /**< [in] fixed full-image aperture-center column */ );

/// Jointly optimize bounded Cartesian position and signed contrast.
/** A fixed-position contrast fit seeds a bounded three-dimensional Nelder-Mead search. Every stage measures the same
 * full-image aperture centered on the initial trial. A second dense+Brent contrast fit at the selected position
 * provides the final photometric result and dense-basin validation. Cartesian position bounds form a square centered
 * on the configured initial trial.
 *
 * \returns the best joint trial, owning local products, and convergence diagnostics
 * \throws std::invalid_argument for an invalid joint optimizer contract
 */
P4PositionContrastOptimizationResult optimizeP4PositionContrast(
    const P4ContrastOptimizerConfig &configuration, /**< [in] total budget, contrast bounds, and tolerances */
    const P4LocalTrial &initialTrial,               /**< [in] center of the Cartesian search box */
    double positionBound,                           /**< [in] symmetric row/column bound in pixels */
    const std::function<P4LocalEvaluation<float>( const P4LocalTrial & )> &evaluate,
    /**< [in] callback evaluating one complete local trial */
    std::size_t modeIndex, /**< [in] selected output-plane index */
    double apertureRadius /**< [in] positive merit-aperture radius in pixels */ );

/// Calculate ordinary delete-one-block jackknife covariance for fitted P4 parameters.
/** Contrast is always included. Cartesian entries are calculated only when `fitPosition` is true and otherwise
 * remain zero. The covariance multiplier is `(B-1)/B` for `B` delete-one-block estimates.
 *
 * \returns the estimate mean, covariance, and diagonal standard errors
 * \throws std::invalid_argument for fewer than two estimates or any non-finite included estimate
 */
P4JackknifeStatistics p4JackknifeStatistics(
    const std::vector<P4JackknifeEstimate> &estimates, /**< [in] one fitted estimate per omitted time block */
    bool fitPosition /**< [in] whether to include the two Cartesian position coordinates */ );

} // namespace improc
} // namespace mx

#endif // P4NegativeOptimizer_hpp
