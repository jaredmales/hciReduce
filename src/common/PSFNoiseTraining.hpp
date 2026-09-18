/** \file PSFNoiseTraining.hpp
 * \brief Declares candidate-excluded annular noise training for local PSF filters.
 */

#ifndef PSFNoiseTraining_hpp
#define PSFNoiseTraining_hpp

#include "P4PSFFilter.hpp"

namespace mx
{
namespace improc
{

/// A circular source or held-out neighborhood in science-image coordinates.
struct PSFNoiseExclusion
{
    double m_row{ 0 };    ///< Center along the first image index.

    double m_column{ 0 }; ///< Center along the second image index.

    double m_radius{ 0 }; ///< Excluded radius in pixels, before dilation by the training footprint.
};

/// Explicit geometry and regularization choices for annular residual-noise training.
struct PSFNoiseTrainingConfig
{
    PSFNoiseKind m_kind{ PSFNoiseKind::pca }; ///< Noise representation; identity bypasses training.

    double m_arcStep{
        0 }; ///< Center arc spacing; zero uses (max(stampRows,stampColumns)-1)/2, bounded below by one pixel.

    double m_guardRadius{ 0 };      ///< Extra holdout beyond the candidate circle, or rectangle with exact exclusions.

    bool m_exactExclusion{ false }; ///< Check actual nonzero interpolation stencils instead of enclosing circles.

    std::size_t m_minimumSamples{ 8 }; ///< Minimum complete training patches; does not imply independent samples.

    std::size_t m_maximumModes{ 3 };   ///< Maximum retained covariance modes; zero uses only the variance floor.

    double m_floorFraction{ 0.1 };     ///< Positive fraction of the median training-pixel variance used as a floor.
};

/// Complete finite training stamps with counts describing all attempted annular locations.
struct PSFNoiseSamples
{
    PSFNoiseModel::matrixT m_samples; ///< One aligned noise stamp per row, in column-major stamp-pixel order.

    PSFNoiseModel::matrixT m_centers; ///< Accepted science-image coordinates, with row and column in two columns.

    std::size_t m_attempted{ 0 };     ///< Equally spaced angular locations before footprint exclusions.

    std::size_t m_excluded{ 0 };      ///< Locations removed by candidate or supplied exclusion neighborhoods.

    std::size_t m_incomplete{ 0 };    ///< Remaining locations rejected for an out-of-bounds or nonfinite interpolant.

    double m_arcStep{ 0 };            ///< Realized arc spacing after rounding the number of angular locations upward.

    double m_footprintRadius{ 0 };    ///< Conservative enclosing radius including the bilinear interpolation stencil.
};

/// Outcome of training and filtering one candidate.
enum class PSFNoiseTrainingStatus
{
    valid,         ///< Training and filter statistics are usable.
    tooFewSamples, ///< Footprint exclusions or invalid pixels left fewer than the required training patches.
    zeroVariance,  ///< Training variance cannot supply a positive finite covariance floor.
    invalidFilter  ///< Response support or filter normalization is unusable.
};

/// Conditional matched-filter result and its local noise-training diagnostics.
struct PSFNoiseTrainingResult
{
    P4PSFFilterResult m_filter; ///< Amplitude and conditional statistics on the candidate's valid support.

    PSFNoiseTrainingStatus m_status{ PSFNoiseTrainingStatus::tooFewSamples }; ///< Explicit reason for unusable results.

    std::size_t m_samples{ 0 };    ///< Accepted patches, including dependence caused by overlap.

    std::size_t m_attempted{ 0 };  ///< Angular locations before exclusions and finite-support checks.

    std::size_t m_excluded{ 0 };   ///< Locations removed by candidate or supplied exclusion neighborhoods.

    std::size_t m_incomplete{ 0 }; ///< Locations rejected for incomplete interpolation support.

    Eigen::Index m_modes{ 0 };     ///< Covariance modes actually retained above the floor.

    double m_varianceFloor{ 0 };   ///< Absolute variance floor in squared science-image units.

    double m_arcStep{ 0 };         ///< Realized training-center arc spacing in pixels.
};

/// Train at fixed radius while withholding the entire candidate and known-source neighborhoods.
/** A training offset is rotated from the candidate angle to the training-center angle before bilinear sampling.
 * Thus every sample is expressed in the candidate's native stamp coordinates; science and response stay on their
 * original pixel lattice. All patches have the same center radius and radial footprint, including across concentric
 * reduction-region boundaries. Interpolation changes the sampled noise and requires held-out calibration.
 * The default exclusions conservatively dilate circles by the whole rotated stamp and interpolation stencil.
 * Exact exclusions instead withhold the candidate's rectangular pixel support (with an optional Euclidean guard)
 * and the pixel centers inside supplied circles, rejecting every patch that reads any of those pixels. Training never
 * zero-fills missing pixels. No taper, effective-sample correction, or significance calibration is implied.
 * \ingroup programming_library
 */
class PSFNoiseTraining
{
  public:
    /// Extract complete, aligned, candidate-excluded annular training patches.
    static PSFNoiseSamples
    sample( P4PSFFilter::imageConstRefT science,  /**< [in] final residual image before matched filtering */
            int stampRows,                        /**< [in] positive odd response-stamp row count */
            int stampColumns,                     /**< [in] positive odd response-stamp column count */
            int sourceRow,                        /**< [in] integer candidate center in science rows */
            int sourceColumn,                     /**< [in] integer candidate center in science columns */
            const PSFNoiseTrainingConfig &config, /**< [in] spacing, guard, and training policy */
            const std::vector<PSFNoiseExclusion>
                &exclusions /**< [in] source and held-out circles in science coordinates */ );

    /// Estimate a local covariance and apply the shared weighted PSF filter.
    /** Insufficient or zero-variance training returns a status without silently falling back to identity weighting.
     * Malformed inputs and failed numerical covariance solves throw. Identity preserves the original filter path.
     */
    static PSFNoiseTrainingResult calculate(
        P4PSFFilter::imageConstRefT science,     /**< [in] final residual image before matched filtering */
        P4PSFFilter::imageConstRefT response,    /**< [in] signed response on the candidate's native pixel lattice */
        P4PSFFilter::validityConstRefT validity, /**< [in] usable response pixels */
        int sourceRow,                           /**< [in] integer candidate row */
        int sourceColumn,                        /**< [in] integer candidate column */
        double minimumSupportFraction,           /**< [in] required usable response fraction */
        const PSFNoiseTrainingConfig &config,    /**< [in] explicit geometry and covariance policy */
        const std::vector<PSFNoiseExclusion>
            &exclusions /**< [in] source and held-out circles in science coordinates */ );
};

} // namespace improc
} // namespace mx

#endif // PSFNoiseTraining_hpp
