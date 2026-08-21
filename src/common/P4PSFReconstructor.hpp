/** \file P4PSFReconstructor.hpp
 * \brief Declares compact detector-to-sky reconstruction for P4 frozen-model PSFs.
 * \author Jared R. Males
 */

#ifndef P4PSFReconstructor_hpp
#define P4PSFReconstructor_hpp

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <mx/improc/eigenCube.hpp>

#include "HCI.hpp"
#include "P4PixelGrid.hpp"

namespace mx
{
namespace improc
{

/// Reconstruct one frame of a sky-centered P4 PSF from compact detector-local responses.
/** The detector-local model matrix stores one column per search-pixel/mode pair, with column-major local-stamp
 * pixels down each column. Reconstruction exactly composes the cubic sampling of each local response with the cubic
 * detector-to-sky interpolation used by `imageRotate`. A sample is valid only when both complete footprints are
 * available and every contributing detector search pixel has a valid local model.
 *
 * \ingroup programming_library
 */
class P4PSFReconstructor
{
  public:
    /// Float image and compact local-model storage.
    using imageT = Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic>;

    /// Byte-valued response validity storage.
    using validityT = Eigen::Array<std::uint8_t, Eigen::Dynamic, Eigen::Dynamic>;

    /// Detector-coordinate to local-model search-index lookup; negative entries are unsupported.
    using searchIndexT = Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic>;

    /// Production cubic interpolation geometry.
    using gridT = P4PixelGridf;

    /// Bounded frame-stack storage used by the configured final estimator.
    using cubeT = eigenCube<float>;

    /// Construct and validate fixed detector, output-stamp, and local-stamp geometry.
    P4PSFReconstructor( int detectorRows,            /**< [in] positive detector row count */
                        int detectorColumns,         /**< [in] positive detector column count */
                        double detectorCenterRow,    /**< [in] finite detector rotation-center row */
                        double detectorCenterColumn, /**< [in] finite detector rotation-center column */
                        int outputStampSize,         /**< [in] positive square final-stamp size */
                        int localStampRows,          /**< [in] positive local-response row count */
                        int localStampColumns /**< [in] positive local-response column count */ );

    /// Reconstruct one derotated PSF response frame about one sky coordinate.
    /** `localModels` must have `localStampRows * localStampColumns` rows and search-major/mode-minor columns.
     * `localValidity` must contain one row per search pixel and one column per mode. On success, both outputs are
     * replaced with square `outputStampSize` arrays. Invalid response samples are stored as zero with validity zero.
     */
    void
    reconstructFrame( imageT &output,                  /**< [out] one reconstructed sky response stamp */
                      validityT &outputValidity,       /**< [out] per-stamp-element validity */
                      const imageT &localModels,       /**< [in] compact detector-local response columns */
                      const validityT &localValidity,  /**< [in] search-pixel/mode validity */
                      const searchIndexT &searchIndex, /**< [in] detector-to-search-index lookup */
                      std::size_t modeIndex,           /**< [in] zero-based requested output-mode index */
                      double sourceSkyRow,             /**< [in] finite final-frame source-center row */
                      double sourceSkyColumn,          /**< [in] finite final-frame source-center column */
                      double derotationAngle /**< [in] finite counterclockwise image rotation in radians */ ) const;

    /// Reconstruct a bounded stack of sky-response stamps for the supplied derotation angles.
    /** Output planes preserve the input angle order. The validity cube contains exact zero/one values and invalid
     * response samples remain zero for safe masked combination.
     */
    void reconstructFrames( cubeT &output,                   /**< [out] response stamp cube by target frame */
                            cubeT &outputValidity,           /**< [out] float zero/one stamp validity cube */
                            const imageT &localModels,       /**< [in] compact detector-local response columns */
                            const validityT &localValidity,  /**< [in] search-pixel/mode validity */
                            const searchIndexT &searchIndex, /**< [in] detector-to-search-index lookup */
                            std::size_t modeIndex,           /**< [in] zero-based requested output-mode index */
                            double sourceSkyRow,             /**< [in] finite final-frame source-center row */
                            double sourceSkyColumn,          /**< [in] finite final-frame source-center column */
                            const std::vector<double> &derotationAngles
                            /**< [in] finite counterclockwise image rotations in radians */ ) const;

    /// Combine one reconstructed frame stack with the same estimators used for the science cube.
    /** `combine::none` is rejected because it does not define one final PSF. Invalid output samples are stored as
     * zero with validity zero; valid output values preserve their signed response and template normalization.
     */
    static void combineFrames( imageT &output,             /**< [out] combined response stamp */
                               validityT &outputValidity,  /**< [out] finite combined-sample validity */
                               cubeT &frames,              /**< [in,out] response frames; invalids may be sanitized */
                               cubeT &frameValidity,       /**< [in] float zero/one response validity */
                               HCI::combine method,        /**< [in] supported final estimator other than none */
                               std::vector<float> weights, /**< [in] empty or one normalized weight per frame */
                               float sigmaThreshold,       /**< [in] sigma threshold; nonpositive falls back to mean */
                               float minimumGoodFraction /**< [in] required valid-frame fraction in `[0,1]` */ );

    /// Reconstruct and combine one final sky-position response using bounded frame scratch.
    void reconstructCombined( imageT &output,                  /**< [out] combined response stamp */
                              validityT &outputValidity,       /**< [out] finite combined-sample validity */
                              const imageT &localModels,       /**< [in] compact detector-local response columns */
                              const validityT &localValidity,  /**< [in] search-pixel/mode validity */
                              const searchIndexT &searchIndex, /**< [in] detector-to-search-index lookup */
                              std::size_t modeIndex,           /**< [in] zero-based requested output-mode index */
                              double sourceSkyRow,             /**< [in] finite final-frame source-center row */
                              double sourceSkyColumn,          /**< [in] finite final-frame source-center column */
                              const std::vector<double> &derotationAngles,
                              /**< [in] finite counterclockwise image rotations in radians */
                              HCI::combine method,               /**< [in] supported final estimator other than none */
                              const std::vector<float> &weights, /**< [in] empty or one normalized weight per frame */
                              float sigmaThreshold,              /**< [in] configured sigma threshold */
                              float minimumGoodFraction /**< [in] required valid-frame fraction in `[0,1]` */ ) const;

  private:
    /// Map one final-frame coordinate to its detector-frame coordinate.
    std::pair<double, double> inverseRotate( double row,    /**< [in] final-frame row */
                                             double column, /**< [in] final-frame column */
                                             double angle /**< [in] counterclockwise rotation in radians */ ) const;

    /// Sample one flattened local-response column, requiring its complete cubic footprint.
    bool sampleLocalResponse( float &value,              /**< [out] interpolated response when valid */
                              const imageT &localModels, /**< [in] flattened local-response columns */
                              Eigen::Index modelColumn,  /**< [in] selected search-pixel/mode column */
                              double deltaRow,           /**< [in] detector target-minus-source row offset */
                              double deltaColumn /**< [in] detector target-minus-source column offset */ ) const;

    int m_detectorRows{ 0 };            ///< Positive detector row count.

    int m_detectorColumns{ 0 };         ///< Positive detector column count.

    double m_detectorCenterRow{ 0 };    ///< Detector rotation-center row.

    double m_detectorCenterColumn{ 0 }; ///< Detector rotation-center column.

    int m_outputStampSize{ 0 };         ///< Square final response-stamp size.

    int m_localStampRows{ 0 };          ///< Support-padded local response row count.

    int m_localStampColumns{ 0 };       ///< Support-padded local response column count.
};

} // namespace improc
} // namespace mx

#endif // P4PSFReconstructor_hpp
