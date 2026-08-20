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
#include "P4PCA.hpp"
#include "ReductionTiming.hpp"
#include "P4PixelGrid.hpp"
#include "P4RotatedGrid.hpp"

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

    int maximumDegreesOfFreedom{ 0 };  ///< Structural mode limit for the selected regression frame and annulus.

    std::size_t estimatedWorkerBytes{ 0 }; ///< Conservative peak private allocation for one regression worker.

    int maximumWorkerCount{ 0 };           ///< OpenMP/search-pixel worker maximum before memory limiting.

    int effectiveWorkerCount{ 0 };         ///< Worker count selected after applying the configured memory policy.

    int minimumNumericalRank{ 0 }; ///< Minimum rank among common-mask-valid local fits, or zero when none are valid.

    std::size_t validLocalFitCount{ 0 };          ///< Number of local fits accepted by the common mask.

    std::size_t maskedLocalFitCount{ 0 };         ///< Number of complete local fits rejected by the common mask.

    std::size_t supportInvalidLocalFitCount{ 0 }; ///< Direct rotated fits lacking all-frame edge or mask support.

    std::vector<std::size_t> rankInvalidCounts;   ///< Rank-insufficient search pixels for each requested output plane.
};

/// Target-only Pixel Prediction Post-Processing reduction orchestrator.
/** P4 learns one in-sample temporal regression per search pixel. Predictor geometry is fixed within each annulus,
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

    realT m_orDeltaRadiusInner{ std::numeric_limits<realT>::quiet_NaN() }; ///< Inward OR radial extent in pixels.

    realT m_orDeltaRadiusOuter{ std::numeric_limits<realT>::quiet_NaN() }; ///< Outward OR radial extent in pixels.

    realT m_orArcHalfWidth{ std::numeric_limits<realT>::quiet_NaN() }; ///< OR azimuthal half-width in pixels; zero uses
                                                                       ///< only the angular cap.

    realT m_orMaxHalfAngle{
        std::numeric_limits<realT>::quiet_NaN() };                ///< OR angular half-width cap in degrees, up to 180.

    realT m_psfRadius{ std::numeric_limits<realT>::quiet_NaN() }; ///< Physical signal-exclusion radius in pixels.

    std::optional<P4ExclusionPolicy> m_exclusionPolicy;           ///< Explicit central signal-exclusion policy.

    realT m_exclusionRadiusBuffer{ std::numeric_limits<realT>::quiet_NaN() }; ///< Added exclusion-radius buffer.

    double m_rankTolerance{ 1e-12 }; ///< Relative numerical-rank threshold in `[0,1)`.

    double m_memoryFraction{ 0.8 };  ///< Fraction of available memory P4 may reserve; zero disables automatic limiting.

    bool m_writeDiagnostics{ false };         ///< Whether checked P4 diagnostic FITS products are written.

    std::string m_diagnosticDirectory{ "." }; ///< Destination directory for enabled diagnostics.

    /// @}

    /** \name P4 Reduction Results - Data
     * @{
     */

    std::vector<std::vector<int>> m_realizedModes;      ///< Retained mode counts indexed by annulus then output plane.

    std::vector<P4RegionStatistics> m_regionStatistics; ///< Per-annulus geometry, rank, and invalidity summaries.

    std::vector<std::vector<std::vector<int>>> m_temporalSelections;
    ///< Per-annulus central and neighboring target-image indices used in detector-frame predictor rows.

    Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic> m_ownership; ///< Owning annulus index, or -1 outside support.

    std::size_t m_availableMemoryBytes{ 0 }; ///< Linux available-memory snapshot, or zero when unavailable/disabled.

    std::size_t m_memoryBudgetBytes{ 0 };    ///< Bytes selected from available memory for future P4 allocations.

    std::size_t m_compactResidualBytes{ 0 }; ///< Estimated bytes retained by compact residual and validity arrays.

    std::size_t m_materializationBytes{ 0 }; ///< Estimated bytes in one full residual/validity materialization pair.

    ReductionTiming m_timing; ///< Instance-owned elapsed and aggregate-worker timing record for the current reduction.

    /// @}

    /// Construct a P4 reduction and enable strict pristine/working target finite-value checks.
    P4Reduction();

    /// Register inherited observation and P4-specific configuration targets.
    void setupConfig( mx::app::appConfigurator &config /**< [in,out] application configuration */ );

    /// Load inherited observation and P4-specific configuration values.
    void loadConfig( mx::app::appConfigurator &config /**< [in] parsed application configuration */ );

    /// Run the configured ordered search annuli through P4 and the shared ADI final lifecycle.
    /** \returns 0 on success, including a configured preprocessing-only stop.
     * \throws mx::exception for invalid configuration, data, geometry, numerical, diagnostic, or output state.
     */
    int reduce();

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
    int finalProcess();

  private:
    /** \cond P4Reduction_test_harness */
    friend class P4ReductionTestAccess;
    /** \endcond */

    /// Convert a supported exclusion policy to its stable configuration spelling.
    static std::string exclusionPolicyString( P4ExclusionPolicy policy /**< [in] supported exclusion policy */ );

    /// Parse an exact exclusion-policy spelling.
    static P4ExclusionPolicy parseExclusionPolicy( const std::string &value /**< [in] configuration spelling */ );

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

    /// Conservatively estimate the peak private allocation for one P4 regression worker.
    static std::size_t estimatedWorkerBytes( std::size_t targetImageCount, /**< [in] temporal sample count */
                                             std::size_t predictorCount,   /**< [in] predictor-column count */
                                             std::size_t modeCount /**< [in] requested residual-column count */ );

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

    /// Write one enabled image-like diagnostic with P4 provenance and checked directory and FITS errors.
    template <typename dataT>
    void writeDiagnostic(
        const std::string &fileName, /**< [in] diagnostic basename */
        const dataT &data,           /**< [in] image-like diagnostic values */
        const fitsHeaderT *additionalHeader = nullptr /**< [in] optional diagnostic-specific cards */ ) const;
};

/// Supported production P4 reduction specialization.
using P4Reductionf = P4Reduction<float, ADIDerotator<float, verbose::vv>, verbose::vv>;

extern template struct P4Reduction<float, ADIDerotator<float, verbose::vv>, verbose::vv>;

} // namespace improc
} // namespace mx

#endif // P4Reduction_hpp
