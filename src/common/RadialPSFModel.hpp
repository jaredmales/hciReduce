/** \file RadialPSFModel.hpp
 * \brief Declares azimuthally averaged PSF response models sampled at discrete radii.
 * \author Jared R. Males
 */

#ifndef RadialPSFModel_hpp
#define RadialPSFModel_hpp

#include <cstddef>
#include <cstdint>
#include <vector>

#include <Eigen/Dense>

namespace mx
{
namespace improc
{

/// One available final-image coordinate at which a PSF response can be measured.
struct RadialPSFSource
{
    std::size_t sourceIndex{ 0 }; ///< Stable caller-owned index of this source coordinate.

    double row{ 0 };              ///< Final-image source-center row.

    double column{ 0 };           ///< Final-image source-center column.
};

/// One selected response measurement assigned to a configured radial bin.
struct RadialPSFSample
{
    std::size_t sourceIndex{ 0 }; ///< Stable caller-owned index of the selected source coordinate.

    std::size_t radiusIndex{ 0 }; ///< Zero-based configured radial-bin index.

    double radius{ 0 };           ///< Actual selected source radius in pixels.

    double angle{ 0 };            ///< Actual source angle in radians from positive column toward positive row.
};

/// Azimuthally average sparse PSF measurements and linearly interpolate the result in radius.
/** Each measured response is rotated so its source separation points along the positive-column axis. Samples in the
 * same configured radial bin are then averaged independently at every valid stamp element. A response requested at
 * an arbitrary source coordinate linearly interpolates the bracketing canonical radial averages and rotates the
 * result back to the source angle. Requests outside the sampled interval use the nearest endpoint. Rotation uses the
 * production cubic-convolution kernel, treats samples beyond the compact stamp as zero, and propagates validity
 * through every nonzero in-bounds interpolation weight.
 *
 * \ingroup programming_library
 */
class RadialPSFModel
{
  public:
    /// Float response-stamp storage.
    using imageT = Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic>;

    /// Byte-valued response validity storage.
    using validityT = Eigen::Array<std::uint8_t, Eigen::Dynamic, Eigen::Dynamic>;

    /// Construct an empty model for strictly increasing discrete radii and fixed stamp geometry.
    RadialPSFModel( std::vector<double> radii, /**< [in] finite strictly increasing nonnegative sample radii */
                    int stampRows,             /**< [in] positive response-stamp row count */
                    int stampColumns /**< [in] positive response-stamp column count */ );

    /// Calculate the requested angular sample count at every configured radius.
    /** Exactly one angular control must be positive. A positive fixed count is repeated at every radius. Otherwise,
     * the count is the minimum needed to keep the circular arc between adjacent samples no larger than the requested
     * step, with one sample at radius zero.
     */
    static std::vector<std::size_t> angularSampleCounts(
        const std::vector<double> &radii, /**< [in] finite strictly increasing nonnegative radii */
        std::size_t samplesPerRadius,     /**< [in] fixed positive count, or zero for arc spacing */
        double sampleArcStep /**< [in] positive maximum azimuthal arc step, or zero for fixed count */ );

    /// Return the radial tolerance used to associate integer detector pixels with requested radial nodes.
    static double detectorPixelRadiusTolerance() noexcept;

    /// Select the nearest distinct, radius-constrained coordinate at every requested angle and configured radius.
    /** Candidate detector pixels must be within detectorPixelRadiusTolerance() of the requested radius. Duplicate
     * coordinates caused by integer-pixel sampling or excluded angular sectors are retained only once within each
     * radial bin. Ties are resolved by the stable source index, so selection is deterministic regardless of candidate
     * ordering.
     */
    static std::vector<RadialPSFSample>
    selectSamples( const std::vector<RadialPSFSource> &sources, /**< [in] available final-image source coordinates */
                   double centerRow,                            /**< [in] finite image-center row */
                   double centerColumn,                         /**< [in] finite image-center column */
                   const std::vector<double> &radii, /**< [in] finite strictly increasing nonnegative radii */
                   std::size_t samplesPerRadius /**< [in] positive number of uniformly spaced requested angles */ );

    /// Select radius-constrained samples using an independently requested angular count at every radius.
    static std::vector<RadialPSFSample>
    selectSamples( const std::vector<RadialPSFSource> &sources, /**< [in] available final-image source coordinates */
                   double centerRow,                            /**< [in] finite image-center row */
                   double centerColumn,                         /**< [in] finite image-center column */
                   const std::vector<double> &radii, /**< [in] finite strictly increasing nonnegative radii */
                   const std::vector<std::size_t> &samplesPerRadius
                   /**< [in] positive requested angular count corresponding one-to-one with radii */ );

    /// Fit canonical radial averages from response measurements corresponding one-to-one with selected samples.
    void fit( const std::vector<imageT> &responses,     /**< [in] measured final-frame response stamps */
              const std::vector<validityT> &validities, /**< [in] per-element measurement validity */
              const std::vector<RadialPSFSample> &samples /**< [in] radial-bin and source-angle metadata */ );

    /// Rotate one response and its validity by a counterclockwise angle about their geometric center.
    static void rotate( imageT &output,                 /**< [out] rotated response stamp */
                        validityT &outputValidity,      /**< [out] rotated per-element validity */
                        const imageT &input,            /**< [in] finite values wherever input validity is nonzero */
                        const validityT &inputValidity, /**< [in] exact zero/nonzero input validity */
                        double angle /**< [in] finite counterclockwise rotation in radians */ );

    /// Evaluate the linearly interpolated canonical radial response at one requested source angle.
    void response( imageT &output,            /**< [out] source-oriented response stamp */
                   validityT &outputValidity, /**< [out] source-oriented per-element validity */
                   double radius,             /**< [in] finite nonnegative requested radius */
                   double angle /**< [in] finite source angle from positive column toward positive row */ ) const;

    /// Return the number of configured radial templates.
    std::size_t radiusCount() const noexcept;

    /// Return one configured radial sample location.
    double radius( std::size_t radiusIndex /**< [in] zero-based radial index */ ) const;

    /// Return one canonical averaged response stamp.
    const imageT &canonicalResponse( std::size_t radiusIndex /**< [in] zero-based radial index */ ) const;

    /// Return the validity of one canonical averaged response stamp.
    const validityT &canonicalValidity( std::size_t radiusIndex /**< [in] zero-based radial index */ ) const;

    /// Return the number of distinct measurements assigned to one radial bin.
    std::size_t sampleCount( std::size_t radiusIndex /**< [in] zero-based radial index */ ) const;

  private:
    std::vector<double> m_radii;             ///< Strictly increasing configured radial sample locations.

    int m_stampRows{ 0 };                    ///< Positive response-stamp row count.

    int m_stampColumns{ 0 };                 ///< Positive response-stamp column count.

    std::vector<imageT> m_responses;         ///< Canonically oriented averaged response by radial bin.

    std::vector<validityT> m_validities;     ///< Per-element validity of each canonical radial response.

    std::vector<std::size_t> m_sampleCounts; ///< Distinct measurement count assigned to each radial bin.
};

} // namespace improc
} // namespace mx

#endif // RadialPSFModel_hpp
