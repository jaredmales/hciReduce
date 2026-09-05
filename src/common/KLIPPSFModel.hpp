/** \file KLIPPSFModel.hpp
 * \brief Declares frozen-basis KLIP response calculations for sparse PSF probes.
 * \author Jared R. Males
 */

#ifndef KLIPPSFModel_hpp
#define KLIPPSFModel_hpp

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "HCI.hpp"
#include "P4PSFModel.hpp"

namespace mx
{
namespace improc
{

/// Propagate a compact PSF probe through a frozen KLIP basis and into a derotated response stamp.
/** The source coordinate is specified in the final sky frame. For each target frame, the source is inverse-rotated
 * into detector coordinates, sampled on one flattened KLIP region, processed by the supported linear centering
 * rule, and projected through the already-computed KL modes. Region responses can then be accumulated directly into
 * a compact sky-frame stamp without materializing a full detector response image.
 *
 * \ingroup programming_library
 */
class KLIPPSFModel
{
  public:
    /// Float image, matrix, and response-stamp storage.
    using imageT = Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic>;

    /// Float flattened regional probe or response storage.
    using vectorT = Eigen::Array<float, Eigen::Dynamic, 1>;

    /// Full-detector lookup from a pixel to its flattened regional index; negative values are outside the region.
    using regionIndexT = Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic>;

    /// Byte-valued stamp-validity storage.
    using validityT = Eigen::Array<std::uint8_t, Eigen::Dynamic, Eigen::Dynamic>;

    /// Construct a frozen-response calculator for fixed detector and output-stamp geometry.
    KLIPPSFModel( const imageT &psfTemplate,   /**< [in] finite centered detector-sampled PSF before KLIP centering */
                  int detectorRows,            /**< [in] positive detector row count */
                  int detectorColumns,         /**< [in] positive detector column count */
                  double detectorCenterRow,    /**< [in] finite detector rotation-center row */
                  double detectorCenterColumn, /**< [in] finite detector rotation-center column */
                  int stampSize /**< [in] positive odd square response-stamp size */ );

    /// Build a checked full-detector lookup for one flattened KLIP region.
    static regionIndexT regionIndex( const std::vector<std::size_t> &indices, /**< [in] unique full-image indices */
                                     int detectorRows,                        /**< [in] positive detector row count */
                                     int detectorColumns /**< [in] positive detector column count */ );

    /// Sample and linearly preprocess one target-frame PSF probe on a flattened KLIP region.
    /** `none` and `imageMean` are supported. The optional cutout mask must have `indices.size()` rows and one column;
     * it is applied with the same ordering as KLIP's regional preprocessing.
     */
    void probe( vectorT &output,                         /**< [out] processed flattened regional probe */
                const std::vector<std::size_t> &indices, /**< [in] flattened full-image region coordinates */
                double sourceSkyRow,                     /**< [in] source-center row in the final sky frame */
                double sourceSkyColumn,                  /**< [in] source-center column in the final sky frame */
                double derotationAngle,                  /**< [in] target-frame counterclockwise derotation angle */
                HCI::meanSub centering,                  /**< [in] supported KLIP centering rule */
                const imageT &cutoutMask = imageT() /**< [in] optional region mask */ ) const;

    /// Calculate residual probes for an ordered list of retained KL-mode counts.
    static void residuals( std::vector<vectorT> &outputs, /**< [out] one residual vector per requested mode count */
                           const vectorT &probe,          /**< [in] processed target-frame probe */
                           const imageT &klModes,         /**< [in] KL modes stored by row in ascending order */
                           const std::vector<int> &modeCounts /**< [in] positive requested retained-mode counts */ );

    /// Calculate a residual probe using an explicitly supplied projection matrix.
    static void projectedResidual( vectorT &output,      /**< [out] residual probe */
                                   const vectorT &probe, /**< [in] processed target-frame probe */
                                   const imageT &projectionMatrix /**< [in] frozen pixel-space projection matrix */ );

    /// Add one flattened regional response to its directly derotated compact sky-frame stamp.
    /** Output storage is initialized on the first call. A stamp element is valid exactly when its cubic detector
     * interpolation anchor is accepted by the production image-rotation boundary policy. Contributions from separate
     * non-overlapping regions may be accumulated by repeated calls.
     */
    void accumulate( imageT &output,                    /**< [in,out] accumulated sky-frame response stamp */
                     validityT &outputValidity,         /**< [in,out] geometric stamp validity */
                     const vectorT &regionResponse,     /**< [in] flattened response for this region */
                     const regionIndexT &regionIndices, /**< [in] detector-to-region lookup */
                     double sourceSkyRow,               /**< [in] source-center row in the final sky frame */
                     double sourceSkyColumn,            /**< [in] source-center column in the final sky frame */
                     double derotationAngle /**< [in] target-frame counterclockwise derotation angle */ ) const;

    /// Return the configured response-stamp size.
    int stampSize() const noexcept;

  private:
    /// Map one sky coordinate into the detector frame for the supplied derotation angle.
    std::pair<double, double> inverseRotate( double row,    /**< [in] sky-frame row */
                                             double column, /**< [in] sky-frame column */
                                             double angle /**< [in] counterclockwise derotation angle */ ) const;

    P4PSFModel m_template;              ///< Prepared compact source template with zero exterior support.

    int m_detectorRows{ 0 };            ///< Positive detector row count.

    int m_detectorColumns{ 0 };         ///< Positive detector column count.

    double m_detectorCenterRow{ 0 };    ///< Detector rotation-center row.

    double m_detectorCenterColumn{ 0 }; ///< Detector rotation-center column.

    int m_stampSize{ 0 };               ///< Positive odd response-stamp width and height.
};

} // namespace improc
} // namespace mx

#endif // KLIPPSFModel_hpp
