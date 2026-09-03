/** \file P4Reduction.hpp
 * \brief Declares the observation orchestrator for Pixel Prediction Post-Processing.
 * \author Jared R. Males
 */

#ifndef P4Reduction_hpp
#define P4Reduction_hpp

#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <type_traits>
#include <vector>

#include <mx/improc/imageTransforms.hpp>

#include "ADIDerotator.hpp"
#include "ADIobservation.hpp"
#include "P4LocalProcessing.hpp"
#include "P4PCA.hpp"
#include "P4PixelGrid.hpp"
#include "P4PSFFilter.hpp"
#include "P4PSFModel.hpp"
#include "P4PSFReconstructor.hpp"
#include "P4RotatedGrid.hpp"
#include "ReductionTiming.hpp"

namespace mx
{
namespace improc
{

/// Coordinate frame in which P4 learns its temporal regression.
/** \ingroup programming_library */
enum class P4RegressionFrame : std::uint8_t
{
    detector, ///< Learn independent coefficient vectors at fixed detector pixels, preserving the original algorithm.
    rotated   ///< Learn at fixed sky pixels after direct inverse-mapped sampling and temporal centering.
};

/** \cond P4Reduction_test_harness */
class P4ReductionTestAccess;
/** \endcond */

/// Aggregate outcome of one P4 search annulus.
/** Counts are accumulated deterministically after the parallel search-pixel loop.
 *
 * \ingroup programming_library
 */
struct P4RegionStatistics
{
    std::size_t searchPixelCount{ 0 }; ///< Number of owned search pixels in the annulus.

    std::size_t targetImageCount{ 0 }; ///< Number of central target-image rows retained for the annulus.

    int temporalNumberImages{ 0 };     ///< Qualifying earlier and later predictor images retained per central target.

    double temporalPsfRadius{ 0 };     ///< Effective temporal exclusion radius in pixels at the mean annulus radius.

    std::size_t predictorCount{ 0 };   ///< Annulus-wide predictor-column count.

    std::size_t maximumExcludedRows{ 0 };   ///< Largest target-specific deleted-row count in the annulus.

    std::size_t exclusionStorageBytes{ 0 }; ///< Persistent bytes owned by compact target exclusions.

    int maximumDegreesOfFreedom{ 0 };       ///< Structural mode limit for the selected regression frame and annulus.

    std::size_t estimatedWorkerBytes{ 0 };  ///< Conservative peak private allocation for one regression worker.

    int maximumWorkerCount{ 0 };            ///< OpenMP/search-pixel worker maximum before memory limiting.

    int effectiveWorkerCount{ 0 };          ///< Worker count selected after applying the configured memory policy.

    int minimumNumericalRank{ 0 }; ///< Minimum rank among common-mask-valid local fits, or zero when none are valid.

    int minimumBaseRank{ 0 };      ///< Minimum factor base rank, or zero when the explicit solver is used.

    std::size_t downdateClampCount{ 0 };    ///< Roundoff-scale negative core eigenvalues clamped in this annulus.

    std::size_t explicitFallbackCount{ 0 }; ///< Target rows recomputed by the explicit oracle for any reason.

    std::size_t rankBoundaryFallbackPixelCount{ 0 }; ///< Search pixels explicitly refit at a rank-boundary ambiguity.

    std::size_t factorValidationFallbackPixelCount{ 0 }; ///< Search pixels explicitly refit after factor validation.

    /// Search pixels explicitly refit after deletion-solver failure.
    std::size_t deletionSolverFallbackPixelCount{ 0 };

    double maximumFactorOrthogonalityDefect{ 0 }; ///< Largest factor defect among validation-fallback search pixels.

    double factorOrthogonalityToleranceAtMaximumDefect{ 0 }; ///< Tolerance paired with the largest factor defect.

    std::size_t validLocalFitCount{ 0 };                     ///< Number of local fits accepted by the common mask.

    std::size_t maskedLocalFitCount{ 0 };         ///< Number of complete local fits rejected by the common mask.

    std::size_t supportInvalidLocalFitCount{ 0 }; ///< Direct rotated fits lacking all-frame edge or mask support.

    std::vector<std::size_t> rankInvalidCounts; ///< Search pixels with one or more rank-invalid target fits per plane.
};

/// One finite-amplitude pixel-local trial evaluated against the pristine target cube.
/** \ingroup programming_library */
struct P4LocalTrial
{
    double separation{ 0 };    ///< Nonnegative sky separation in pixels.

    double positionAngle{ 0 }; ///< Position angle in degrees east of north.

    double contrast{ 0 };      ///< Signed source contrast, including zero for the unperturbed baseline.
};

/// Owning result of one finite-amplitude pixel-local P4 evaluation.
/** \tparam realT local residual and validity storage type
 * \ingroup programming_library
 */
template <typename realT>
struct P4LocalEvaluation
{
    eigenCube<realT> residual;          ///< Combined local residual cube in configured mode order.

    eigenCube<realT> validity;          ///< Nonzero combined local validity cube matching `residual`.

    int originRow{ 0 };                 ///< Full-image row occupied by local element `(0,0)`.

    int originColumn{ 0 };              ///< Full-image column occupied by local element `(0,0)`.

    double sourceRow{ 0 };              ///< Continuous full-image trial-source row.

    double sourceColumn{ 0 };           ///< Continuous full-image trial-source column.

    std::size_t searchCount{ 0 };       ///< Unique detector regressions evaluated by the local path.

    std::size_t sparseSampleCount{ 0 }; ///< Requested detector-pixel/frame residual pairs.

    double elapsedSeconds{ 0 };         ///< Total elapsed time for this evaluation.

    ReductionTiming timing;             ///< Detailed algorithm timing for this evaluation.
};

/// Target-only Pixel Prediction Post-Processing reduction orchestrator.
/** P4 learns one temporal regression per search pixel. Predictor geometry is fixed within each annulus,
 * while numerical workspaces are private to OpenMP workers. The initial supported implementation uses float image
 * storage, all-double normal equations, and float cubic-convolution interpolation. Detector-frame regression is the
 * compatibility default; rotated-frame regression directly samples fixed sky coordinates from each preprocessed,
 * unrotated detector frame.
 *
 * \tparam _realT image-storage arithmetic type
 * \tparam _derotFunctObj ADI derotation policy
 * \tparam verboseT mxlib exception verbosity policy
 * \ingroup programming_library
 */
template <typename _realT, class _derotFunctObj, class verboseT>
struct P4Reduction : public ADIobservation<_realT, _derotFunctObj, verboseT>
{
    static_assert( std::is_same_v<_realT, float>, "initial P4Reduction supports float image storage only" );

    /// Image-storage arithmetic type.
    using realT = _realT;

    /// Inherited FITS-header type.
    using fitsHeaderT = typename ADIobservation<realT, _derotFunctObj, verboseT>::fitsHeaderT;

    /// Inherited two-dimensional image type.
    using imageT = typename ADIobservation<realT, _derotFunctObj, verboseT>::imageT;

    /// Compact local-PSF validity storage type.
    using psfValidityT = Eigen::Array<std::uint8_t, Eigen::Dynamic, Eigen::Dynamic>;

    /// Fixed float cubic-convolution P4 pixel-grid used by the initial production implementation.
    using pixelGridT = P4PixelGridf;

    /// Fixed float direct sampler used by rotated-frame regression.
    using rotatedGridT = P4RotatedGrid;

    /** \name P4 Configuration - Data
     * @{
     */

    std::vector<realT> m_minRadius;     ///< Inclusive inner radii of the ordered search annuli.

    std::vector<realT> m_maxRadius;     ///< Exclusive outer radii of the ordered search annuli.

    std::vector<realT> m_modeFractions; ///< Strictly increasing PCA fractions defining output planes.

    P4RegressionFrame m_regressionFrame{ P4RegressionFrame::detector }; ///< Frame of the learned regression.

    int m_numberImages{ 0 }; ///< Qualifying earlier and later images appended to each detector-frame predictor row.

    realT m_minDPx{ 0 };     ///< Minimum target/reference displacement interpreted by `m_excludeMethod`.

    HCI::exclude m_excludeMethod{ HCI::exclude::none }; ///< KLIP-compatible target-frame exclusion method.

    P4ExclusionSolver m_exclusionSolver{ P4ExclusionSolver::explicitRefit }; ///< Target-frame exclusion solver.

    mx::math::svdDeletionBackend m_deletionBackend{ mx::math::svdDeletionBackend::leadingCovariance };
    ///< mxlib row-deletion backend used by the exact factor solver.

    realT m_orDeltaRadiusInner{ std::numeric_limits<realT>::quiet_NaN() }; ///< Inward OR radial extent in pixels.

    realT m_orDeltaRadiusOuter{ std::numeric_limits<realT>::quiet_NaN() }; ///< Outward OR radial extent in pixels.

    realT m_orArcHalfWidth{ std::numeric_limits<realT>::quiet_NaN() }; ///< OR azimuthal half-width in pixels; zero uses
                                                                       ///< only the angular cap.

    realT m_orMaxHalfAngle{
        std::numeric_limits<realT>::quiet_NaN() };                ///< OR angular half-width cap in degrees, up to 180.

    realT m_psfRadius{ std::numeric_limits<realT>::quiet_NaN() }; ///< Physical signal-exclusion radius in pixels.

    std::string m_psfFile;              ///< Optional post-preprocessing PSF template enabling frozen-model calculation.

    int m_psfStampSize{ 0 };            ///< Square frozen-model PSF stamp size; required when `m_psfFile` is set.

    bool m_outputPSFModels{ false };    ///< Whether to reconstruct and write compact final-frame PSF fields.

    bool m_psfFilter{ false };          ///< Whether to apply the spatially variable normalized PSF filter.

    realT m_psfFilterMinGoodFract{ 1 }; ///< Minimum usable local-stamp fraction required by PSF filtering.

    std::string m_psfOutputPrefix{ "p4PSF_" }; ///< Prefix for compact products in the final image's output directory.

    int m_localStampSize{ 0 }; ///< Square pixel-local result and nominal source-crop width; zero disables the path.

    std::optional<P4ExclusionPolicy> m_exclusionPolicy; ///< Explicit central signal-exclusion policy.

    realT m_exclusionRadiusBuffer{ std::numeric_limits<realT>::quiet_NaN() }; ///< Added exclusion-radius buffer.

    double m_rankTolerance{ 1e-12 }; ///< Relative numerical-rank threshold in `[0,1)`.

    double m_memoryFraction{ 0.8 };  ///< Fraction of available memory P4 may reserve; zero disables automatic limiting.

    bool m_writeDiagnostics{ false };         ///< Whether checked P4 diagnostic FITS products are written.

    std::string m_diagnosticDirectory{ "." }; ///< Optional explicit destination; `.` selects the auxiliary directory.

    /// @}

    /** \name P4 Reduction Results - Data
     * @{
     */

    std::vector<std::vector<int>> m_realizedModes;      ///< Retained mode counts indexed by annulus then output plane.

    std::vector<P4RegionStatistics> m_regionStatistics; ///< Per-annulus geometry, rank, and invalidity summaries.

    std::vector<std::vector<std::vector<int>>> m_temporalSelections;
    ///< Per-annulus central and neighboring target-image indices used in detector-frame predictor rows.

    std::vector<imageT> m_localPSFModels;
    ///< Compact same-image local responses by annulus; columns are search-major then mode-minor.

    std::vector<imageT> m_localPSFTemporalCoefficients;
    ///< Temporal predictor coefficients by annulus; rows are slot-major offsets and columns are search-major modes.

    std::vector<psfValidityT> m_localPSFValidity;
    ///< Rank and geometry validity by annulus, search pixel, and requested output mode.

    std::vector<std::size_t> m_localPSFComponentCounts;
    ///< Same-image plus realized temporal response-component count for every annulus.

    int m_localPSFRows{ 0 };       ///< Phase-matched, support-padded local response row count.

    int m_localPSFColumns{ 0 };    ///< Phase-matched, support-padded local response column count.

    int m_psfTemplateRows{ 0 };    ///< Configured PSF-template row count used by the current reduction.

    int m_psfTemplateColumns{ 0 }; ///< Configured PSF-template column count used by the current reduction.

    Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic> m_ownership; ///< Owning annulus index, or -1 outside support.

    std::size_t m_availableMemoryBytes{ 0 };  ///< Linux available-memory snapshot, or zero when unavailable/disabled.

    std::size_t m_memoryBudgetBytes{ 0 };     ///< Bytes selected from available memory for future P4 allocations.

    std::size_t m_compactResidualBytes{ 0 };  ///< Estimated bytes retained by compact residual and validity arrays.

    std::size_t m_targetExclusionBytes{ 0 };  ///< Bytes retained by compact target-specific deletion patterns.

    std::size_t m_materializationBytes{ 0 };  ///< Estimated bytes in one full residual/validity materialization pair.

    std::size_t m_localPSFBytes{ 0 };         ///< Bytes retained by opt-in compact local PSF stamps and validity.

    std::size_t m_psfModeBatchSize{ 0 };      ///< Target-held-out PSF modes retained in one memory-bounded batch.

    std::size_t m_psfDowndateClampCount{ 0 }; ///< Eigenvalue clamps accumulated by held-out PSF response passes.

    std::size_t m_psfExplicitFallbackCount{ 0 }; ///< Row-refit invocations across held-out PSF response batches.

    std::size_t m_psfRankBoundaryFallbackPixelCount{ 0 }; ///< PSF search/batch fits falling back at a rank boundary.

    std::size_t m_psfFactorValidationFallbackPixelCount{ 0 }; ///< PSF search/batch fits with invalid base factors.

    std::size_t m_psfDeletionSolverFallbackPixelCount{ 0 }; ///< PSF search/batch fits with structured-solver fallback.

    double m_psfMaximumFactorOrthogonalityDefect{ 0 };      ///< Largest PSF-pass base-factor orthogonality defect.

    double m_psfFactorOrthogonalityToleranceAtMaximumDefect{ 0 }; ///< Tolerance paired with the largest PSF defect.

    std::size_t m_psfModelBytes{ 0 };          ///< Bytes retained by the shared precomputed PSF-template model.

    std::size_t m_psfReconstructionBytes{ 0 }; ///< Conservative peak scratch for one final PSF mode.

    std::size_t m_psfFilterBytes{ 0 };     ///< Bytes retained by filtered, normalization, support, and validity cubes.

    eigenCube<realT> m_localFinalValidity; ///< Combined local-result validity cube in output-mode order.

    std::optional<std::vector<int>> m_localIncludedFrames; ///< Physical target frames retained by one resampled local
                                                           ///< evaluation; unset retains the complete sequence.

    int m_localOriginRow{ 0 };                             ///< Full-image row occupied by local result element `(0,0)`.

    int m_localOriginColumn{ 0 };              ///< Full-image column occupied by local result element `(0,0)`.

    double m_localSourceRow{ 0 };              ///< Continuous full-image row of the configured trial source.

    double m_localSourceColumn{ 0 };           ///< Continuous full-image column of the configured trial source.

    int m_localTemplateRows{ 0 };              ///< Actual phase-preserving internal trial-template crop rows.

    int m_localTemplateColumns{ 0 };           ///< Actual phase-preserving internal trial-template crop columns.

    std::size_t m_localSearchCount{ 0 };       ///< Unique detector regressions evaluated by pixel-local processing.

    std::size_t m_localSparseSampleCount{ 0 }; ///< Requested detector-pixel/frame residual pairs.

    std::size_t m_localGeometryBytes{ 0 };     ///< Retained sparse request and output-dependency storage.

    ReductionTiming m_timing; ///< Instance-owned elapsed and aggregate-worker timing record for the current reduction.

    /// @}

    /// Construct a P4 reduction and enable strict pristine/working target finite-value checks.
    P4Reduction();

    /// Register inherited observation and P4-specific configuration targets.
    void setupConfig( mx::app::appConfigurator &config /**< [in,out] application configuration */ );

    /// Load inherited observation and P4-specific configuration values.
    void loadConfig( mx::app::appConfigurator &config /**< [in] parsed application configuration */ );

    /// Defer configured target-source injection when pixel-local processing samples sources on demand.
    bool deferTargetFakeInjection() const override;

    /// Run the configured ordered search annuli through P4 and the shared ADI final lifecycle.
    /** \returns 0 on success, including a configured preprocessing-only stop.
     * \throws mx::exception for invalid configuration, data, geometry, numerical, diagnostic, or output state.
     */
    int reduce();

    /// Evaluate one local trial without writing intermediate products or changing the configured trial tuple.
    /** The returned cubes own their storage and remain valid after later evaluations. The pristine loaded target cube,
     * file selection, headers, angles, configured fake tuple, and output controls are preserved.
     *
     * \returns an owning local residual, validity, coordinate, count, and timing result
     * \throws mx::exception for invalid local configuration, input, geometry, or numerical state
     */
    P4LocalEvaluation<realT> evaluateLocal( const P4LocalTrial &trial /**< [in] replacement single local trial */ );

    /// Evaluate one local trial using an explicit subset of physical target frames.
    /** The selected frames are used as both central targets and temporal predictors. The selection must be nonempty,
     * strictly increasing, unique, and in range after the observation has loaded. All state guarantees of the
     * complete-sequence overload apply.
     *
     * \returns an owning local result combined from only the selected frames
     * \throws mx::exception for an invalid frame selection or local reduction state
     */
    P4LocalEvaluation<realT> evaluateLocal(
        const P4LocalTrial &trial, /**< [in] replacement single local trial */
        const std::vector<int> &includedFrames /**< [in] ordered physical target-frame indices to retain */ );

    /// Return the current loaded target-frame count, or zero before target loading.
    std::size_t targetFrameCount() const;

    /// Return the auxiliary-product directory derived from the currently resolved final-image path.
    std::string auxiliaryOutputDirectory() const;

    /// Run explicit ordered search-annulus radii through P4 and the shared ADI final lifecycle.
    /** The supplied vectors replace the configured `geom.minRadius` and `geom.maxRadius` values.
     *
     * \returns 0 on success, including a configured preprocessing-only stop.
     */
    int regions( const std::vector<realT> &minimumRadii, /**< [in] inclusive ordered inner radii */
                 const std::vector<realT> &maximumRadii /**< [in] exclusive ordered outer radii */ );

    /// Print the current reduction-stage timing summary.
    void dump_times() const;

    /// Append complete P4 configuration and realized reduction outcomes to a FITS header.
    void appendReductionHeader( fitsHeaderT &head /**< [in,out] header receiving P4 provenance */ ) const;

    /// Apply the shared ADI final-processing lifecycle with P4 provenance.
    /** \returns 0 on success.
     * \throws mx::exception if the regression frame is invalid, rotated regression is paired with post-median
     * subtraction, or final processing or output fails.
     */
    int finalProcess( bool reportProgress = true /**< [in] whether to report final-processing stages to stderr */ );

  private:
    /** \cond P4Reduction_test_harness */
    friend class P4ReductionTestAccess;
    /** \endcond */

    /// Convert a supported exclusion policy to its stable configuration spelling.
    static std::string exclusionPolicyString( P4ExclusionPolicy policy /**< [in] supported exclusion policy */ );

    /// Parse an exact exclusion-policy spelling.
    static P4ExclusionPolicy parseExclusionPolicy( const std::string &value /**< [in] configuration spelling */ );

    /// Convert a supported target-frame exclusion solver to its stable configuration spelling.
    static std::string exclusionSolverString( P4ExclusionSolver solver /**< [in] supported exclusion solver */ );

    /// Parse an exact target-frame exclusion-solver spelling.
    static P4ExclusionSolver parseExclusionSolver( const std::string &value /**< [in] configuration spelling */ );

    /// Convert a supported mxlib row-deletion backend to its stable configuration spelling.
    static std::string
    deletionBackendString( mx::math::svdDeletionBackend backend /**< [in] supported row-deletion backend */ );

    /// Parse an exact supported mxlib row-deletion-backend spelling.
    static mx::math::svdDeletionBackend
    parseDeletionBackend( const std::string &value /**< [in] configuration spelling */ );

    /// Convert a supported regression frame to its stable configuration spelling.
    static std::string regressionFrameString( P4RegressionFrame frame /**< [in] supported regression frame */ );

    /// Parse an exact regression-frame spelling.
    static P4RegressionFrame parseRegressionFrame( const std::string &value /**< [in] configuration spelling */ );

    /// Validate configuration that is independent of loaded image dimensions.
    void validateConfiguration() const;

    /// Report whether the complete target cube contains only finite values under fast-math.
    bool targetCubeFinite() const;

    /// Promote one sampled predictor only when float interpolation remained finite.
    static double checkedPredictorPromotion( realT value /**< [in] sampled predictor value */ );

    /// Convert one all-double residual only when it remains finite in image storage.
    static realT checkedResidualCast( double value /**< [in] all-double residual value */ );

    /// Return the structurally available degrees of freedom after checking the P4PCA integer boundary.
    static int checkedMaximumDegreesOfFreedom( int imageCount, /**< [in] positive temporal sample count */
                                               std::size_t predictorCount, /**< [in] predictor-column count */
                                               bool temporallyCentered = false /**< [in] whether centering removes one
                                                                                         temporal degree of freedom */ );

    /// Build compact KLIP-compatible deleted-row patterns for every central target row.
    static P4TargetExclusions targetExclusions(
        const std::vector<std::vector<int>> &selections, /**< [in] central physical image in each row's first slot */
        HCI::exclude method,                             /**< [in] none, pixel, angle, or image-number exclusion */
        realT minDPx,                                    /**< [in] nonnegative exclusion threshold */
        realT minimumRadius,                             /**< [in] annulus inner radius for pixel conversion */
        const std::vector<double> &derotationAngles /**< [in] radians angles indexed by physical target image */ );

    /// Conservatively estimate the peak private allocation for one P4 regression worker.
    static std::size_t estimatedWorkerBytes( std::size_t targetImageCount, /**< [in] temporal sample count */
                                             std::size_t predictorCount,   /**< [in] predictor-column count */
                                             std::size_t modeCount,        /**< [in] requested residual-column count */
                                             bool includeCoefficients = false, /**< [in] include optional K-by-mode
                                                                                          coefficient output */
                                             std::size_t psfStampPixels = 0,   /**< [in] optional float PSF scratch */
                                             bool exactHeldOut = false, /**< [in] include target-exclusion scratch */
                                             P4ExclusionSolver exclusionSolver = P4ExclusionSolver::explicitRefit,
                                             /**< [in] target-exclusion numerical implementation */
                                             std::size_t maximumDeletedRows = 0,
                                             /**< [in] maximum deletion count for factor workspace */
                                             std::size_t maximumRetainedMode = 0
                                             /**< [in] largest retained mode, or zero for a conservative maximum */,
                                             std::size_t probeSampleCount = 0
                                             /**< [in] frozen-probe response samples, or zero when disabled */ );

    /// Return the phase-matched local-model dimension needed for final-stamp reconstruction.
    static int localPSFModelDimension( int outputStampSize, /**< [in] positive square final-stamp size */
                                       int templateDimension /**< [in] positive template rows or columns */ );

    /// Return the exact retained byte count for one compact local-PSF annulus.
    static std::size_t localPSFBytes( std::size_t searchPixelCount, /**< [in] owned search-pixel count */
                                      std::size_t modeCount,        /**< [in] requested output-mode count */
                                      std::size_t temporalPredictorCount,
                                      /**< [in] retained temporal coefficients per pixel and mode */
                                      int stampRows, /**< [in] positive local-stamp row count */
                                      int stampColumns /**< [in] positive local-stamp column count */ );

    /// Conservatively estimate peak scratch used to reconstruct and persist one final PSF mode.
    static std::size_t psfReconstructionBytes( std::size_t searchPixelCount, /**< [in] owned output positions */
                                               std::size_t targetImageCount, /**< [in] target-frame count */
                                               std::size_t temporalPredictorCount,
                                               /**< [in] retained temporal coefficients per search pixel */
                                               int outputStampSize, /**< [in] square final-stamp size */
                                               int localStampRows,  /**< [in] support-padded local rows */
                                               int localStampColumns /**< [in] support-padded local columns */ );

    /// Return the exact retained byte count for the four full-frame PSF-filter products.
    static std::size_t psfFilterBytes( int rows,    /**< [in] positive final-image rows */
                                       int columns, /**< [in] positive final-image columns */
                                       std::size_t modeCount /**< [in] positive final-image mode count */ );

    /// Select the number of workers that fit after persistent future allocations.
    static int memoryLimitedWorkerCount( int requestedWorkers,    /**< [in] positive OpenMP/search-pixel maximum */
                                         std::size_t budgetBytes, /**< [in] bytes available for future P4 allocations */
                                         std::size_t persistentBytes, /**< [in] bytes retained while workers run */
                                         std::size_t workerBytes /**< [in] conservative private bytes per worker */ );

    /// Assign one search pixel to exactly one annulus.
    static void claimOwnership( Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic> &ownership,
                                /**< [in,out] ownership image initialized to -1 */
                                const P4PixelCoordinate &coordinate, /**< [in] owned search coordinate */
                                int region /**< [in] nonnegative annulus index */ );

    /// Assemble and solve one detector-frame search-pixel regression through the authoritative production path.
    void fitDetectorSearch(
        P4PCAResult &result,                             /**< [out] residuals, rank, and mode status */
        P4PCA::matrixT &predictors,                      /**< [in,out] reusable target-by-predictor matrix */
        P4PCA::vectorT &target,                          /**< [in,out] reusable target time-series vector */
        P4PCA::matrixT *coefficients,                    /**< [out] optional predictor coefficients, or nullptr */
        P4PCA::matrixT *probeResiduals,                  /**< [out] optional held-out frozen-probe responses */
        const pixelGridT &grid,                          /**< [in] configured detector-frame P4 geometry */
        std::size_t search,                              /**< [in] annulus-local search-pixel index */
        const std::vector<std::vector<int>> &selections, /**< [in] central and neighboring images per PCA row */
        const std::vector<P4PixelCoordinate> &temporalOffsets,
        /**< [in] direct additional-image predictor offsets */
        const std::vector<int> &modes, /**< [in] requested realized integer modes */
        P4PCA::workspaceT &workspace,  /**< [in,out] reusable all-FP64 or factor-deletion workspace */
        detail::P4PCAMixedWorkspace &mixedWorkspace,
        /**< [in,out] reusable FP32-calculation and FP64-eigensolver workspace */
        P4PCADowndateWorkspace *downdateWorkspace,
        /**< [in,out] optional reusable factor-deletion workspace */
        P4PCATiming &timing,              /**< [out] PCA phase timing */
        double &sameImageSamplingSeconds, /**< [out] target and same-image OR worker seconds */
        double &temporalSamplingSeconds,  /**< [out] additional-image sampling worker seconds */
        const P4TargetExclusions *exclusions,
        /**< [in] optional per-target deleted rows; nullptr selects in-sample fitting */
        const P4PCA::matrixT *probePredictors, /**< [in] optional frozen predictor responses */
        const P4PCA::vectorT *probeTarget,     /**< [in] optional direct frozen target response */
        const P4TrialSource *trialSource = nullptr /**< [in] optional finite-amplitude trial perturbation */ ) const;

    /// Load one finite per-frame fake scale vector using the inherited filename-matching convention.
    std::vector<realT> localTrialScales() const;

    /// Execute one configured pixel-local finite-amplitude refit and final combination.
    int processLocalTrial(
        const std::vector<pixelGridT> &grids, /**< [in] configured detector-frame annular geometry */
        const std::vector<P4PixelCoordinate> &temporalOffsets,
        /**< [in] direct additional-image predictor offsets */
        const std::vector<double> &derotationAngles, /**< [in] one finite radians angle per target frame */
        const std::vector<P4TargetExclusions> &regionExclusions
        /**< [in] optional target-specific deleted rows by annulus */ );

    /// Atomically write the combined local validity cube next to its final residual product.
    void writeLocalValidity( const std::string &finalImagePath, /**< [in] resolved local residual FITS path */
                             const fitsHeaderT &finalHeader /**< [in] P4 provenance appended to the final header */ );

    /// Reconstruct, optionally persist, and optionally filter the final-frame PSF field one mode at a time.
    void processPSFProducts(
        const std::vector<pixelGridT> &grids, /**< [in] retained detector-frame annulus geometry */
        const P4PSFModel &psfModel,           /**< [in] prepared full-support source template */
        const std::vector<P4TargetExclusions> &regionExclusions,
        /**< [in] optional target-specific deleted rows by annulus */
        const std::string &finalImagePath, /**< [in] resolved path supplying filter naming and provenance */
        const fitsHeaderT &finalHeader /**< [in] ADI and P4 cards mirrored from the ordinary final image */ );

    /// Calculate one memory-bounded batch of target-held-out frozen PSF responses.
    void calculateHeldOutPSFBatch(
        std::vector<imageT> &localModels,         /**< [out] one stamp matrix per requested output mode */
        std::vector<psfValidityT> &localValidity, /**< [out] one search-by-target validity matrix per mode */
        const std::vector<pixelGridT> &grids,     /**< [in] retained detector-frame annulus geometry */
        const P4PSFModel &psfModel,               /**< [in] prepared full-support source template */
        const std::vector<P4TargetExclusions> &regionExclusions,
        /**< [in] nonempty target-specific deleted rows by annulus */
        const std::vector<std::size_t> &regionOffsets, /**< [in] global first-search index by annulus */
        std::size_t searchPixelCount,                  /**< [in] total global search-pixel count */
        std::size_t firstOutput,                       /**< [in] first requested output-mode index */
        std::size_t outputCount /**< [in] positive number of consecutive modes in this batch */ );

    /// Write one enabled image-like diagnostic with P4 provenance and checked directory and FITS errors.
    template <typename dataT>
    void writeDiagnostic(
        const std::string &fileName, /**< [in] diagnostic basename */
        const dataT &data,           /**< [in] image-like diagnostic values */
        const fitsHeaderT *additionalHeader = nullptr /**< [in] optional diagnostic-specific cards */ ) const;
};

/// Supported production P4 reduction specialization.
using P4Reductionf = P4Reduction<float, ADIDerotator<float, verbose::vv>, verbose::vv>;

#ifdef HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION
/** \cond P4Reduction_precision_experiment */
namespace detail
{

/// Run one P4 reduction through an explicitly selected experimental PCA scalar policy.
/** Non-DD policies currently require disabled automatic memory limiting and reject factor deletion, shared PSF
 * coefficients, and local-trial processing before the reduction object is mutated.
 */
int p4ReductionReduceExperimental(
    P4Reductionf &reduction, /**< [in,out] configured production-layout reduction object */
    P4PCAPrecisionPolicy precisionPolicy /**< [in] PCA calculation/eigensolver scalar combination */ );

} // namespace detail
/** \endcond */
#endif

extern template struct P4Reduction<float, ADIDerotator<float, verbose::vv>, verbose::vv>;

} // namespace improc
} // namespace mx

#endif // P4Reduction_hpp
