/** \file P4PixelGrid.hpp
 * \brief Declares the local pixel geometry used by Pixel Prediction Post-Processing.
 * \author Jared R. Males
 */

#ifndef P4PixelGrid_hpp
#define P4PixelGrid_hpp

#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include <mx/improc/imageTransforms.hpp>

namespace mx
{
namespace improc
{

/// Central signal-exclusion rule applied while constructing a P4 predictor grid.
/** Both rules compare against the same explicitly configured effective exclusion radius. Under `sampleCenter`, only
 * the interpolation sample center is excluded. Under `kernelSupport`, a canonical predictor column is excluded from
 * the whole annulus when any mapped nominal interpolation footprint overlaps the exclusion disk.
 *
 * \ingroup programming_library
 */
enum class P4ExclusionPolicy : std::uint8_t
{
    sampleCenter, ///< Exclude predictor sample centers inside the effective disk.
    kernelSupport ///< Exclude complete nominal interpolation footprints that overlap the effective disk.
};

/// Common-mask reason that makes one P4 search-pixel fit invalid.
/** The reason applies to a complete local fit without changing the annulus-wide predictor count.
 *
 * \ingroup programming_library
 */
enum class P4PixelInvalidReason : std::uint8_t
{
    none,                    ///< The target and every predictor footprint are valid.
    targetMasked,            ///< The target pixel is masked.
    predictorMasked,         ///< At least one predictor footprint touches a masked pixel.
    targetAndPredictorMasked ///< Both the target and at least one predictor footprint are masked.
};

/// Complete, explicit geometry configuration for one P4 search annulus.
/** There are deliberately no default geometry or exclusion-policy values.
 *
 * \ingroup programming_library
 */
struct P4PixelGridRegion
{
    /// Construct a complete P4 region description.
    P4PixelGridRegion( double searchInnerRadius,            /**< [in] inclusive search-annulus inner radius */
                       double searchOuterRadius,            /**< [in] exclusive search-annulus outer radius */
                       double optimizationDeltaRadiusInner, /**< [in] inward optimization-wedge radial extent */
                       double optimizationDeltaRadiusOuter, /**< [in] outward optimization-wedge radial extent */
                       double optimizationArcHalfWidth,     /**< [in] azimuthal half-width as arc length in pixels; zero
                                                               uses only the angular cap */
                       double optimizationMaxHalfAngle, /**< [in] maximum azimuthal half-angle in degrees in (0,180] */
                       double psfRadius,                /**< [in] physical central exclusion radius */
                       P4ExclusionPolicy exclusionPolicy, /**< [in] required central-exclusion rule */
                       double exclusionRadiusBuffer /**< [in] explicit comparison buffer added to \p psfRadius */ );

    double searchInnerRadius;            ///< Inclusive search-annulus inner radius in pixels.

    double searchOuterRadius;            ///< Exclusive search-annulus outer radius in pixels.

    double optimizationDeltaRadiusInner; ///< Inward optimization-wedge radial extent in pixels.

    double optimizationDeltaRadiusOuter; ///< Outward optimization-wedge radial extent in pixels.

    double optimizationArcHalfWidth;     ///< Optimization-wedge azimuthal half-width as arc length in pixels; zero uses
                                         ///< only the angular cap.

    double optimizationMaxHalfAngle;     ///< Maximum optimization-wedge half-angle in degrees, up to 180.

    double psfRadius;                    ///< Physical central signal-exclusion radius in pixels.

    P4ExclusionPolicy exclusionPolicy;   ///< Required central-exclusion rule.

    double exclusionRadiusBuffer;        ///< Explicit radius-comparison buffer in pixels.
};

/// Immutable integer image coordinate used by P4 geometry records.
/** \ingroup programming_library */
class P4PixelCoordinate
{
  public:
    /// Construct an integer row/column coordinate.
    P4PixelCoordinate( int row, /**< [in] zero-based image row */
                       int column /**< [in] zero-based image column */ );

    /// Return the zero-based image row.
    int row() const noexcept;

    /// Return the zero-based image column.
    int column() const noexcept;

  private:
    int m_row;    ///< Zero-based image row.

    int m_column; ///< Zero-based image column.
};

/// Immutable common-mask validity record for one P4 search pixel.
/** \ingroup programming_library */
class P4SearchPixelRecord
{
  public:
    /// Construct a search-pixel record.
    P4SearchPixelRecord( const P4PixelCoordinate &coordinate, /**< [in] target image coordinate */
                         P4PixelInvalidReason invalidReason /**< [in] complete local-fit validity reason */ );

    /// Return the target image coordinate.
    const P4PixelCoordinate &coordinate() const noexcept;

    /// Return the complete local-fit validity reason.
    P4PixelInvalidReason invalidReason() const noexcept;

    /// Report whether the complete local fit is valid.
    bool valid() const noexcept;

    /// Report whether the target pixel is masked.
    bool targetMasked() const noexcept;

    /// Report whether any mapped predictor footprint touches the common mask.
    bool predictorMasked() const noexcept;

  private:
    P4PixelCoordinate m_coordinate;       ///< Target image coordinate.

    P4PixelInvalidReason m_invalidReason; ///< Complete local-fit validity reason.
};

template <typename transformT>
class P4PixelGrid;

/** \cond P4PixelGrid_test_harness */
class P4PixelGridTestAccess;
/** \endcond */

/// Immutable interpolation coordinate and kernel for one mapped P4 predictor sample.
/** \ingroup programming_library */
template <typename transformT>
class P4InterpolationRecord
{
  public:
    /// Arithmetic type used by the interpolation transform.
    using realT = typename transformT::arithT;

    /// Fixed-size interpolation kernel type.
    using kernelT = Eigen::Array<realT, transformT::width, transformT::width>;

    /// Return the mapped subpixel row coordinate.
    double row() const noexcept;

    /// Return the mapped subpixel column coordinate.
    double column() const noexcept;

    /// Return the first image row consumed by the nominal interpolation footprint.
    int footprintRow() const noexcept;

    /// Return the first image column consumed by the nominal interpolation footprint.
    int footprintColumn() const noexcept;

    /// Return the precomputed interpolation kernel.
    const kernelT &kernel() const noexcept;

  private:
    friend class P4PixelGrid<transformT>;

    /// Construct a complete interpolation record.
    P4InterpolationRecord( double row,          /**< [in] mapped subpixel row */
                           double column,       /**< [in] mapped subpixel column */
                           int footprintRow,    /**< [in] first nominal-footprint row */
                           int footprintColumn, /**< [in] first nominal-footprint column */
                           const kernelT &kernel /**< [in] precomputed interpolation kernel */ );

    double m_row;          ///< Mapped subpixel row coordinate.

    double m_column;       ///< Mapped subpixel column coordinate.

    int m_footprintRow;    ///< First nominal-footprint image row.

    int m_footprintColumn; ///< First nominal-footprint image column.

    kernelT m_kernel;      ///< Precomputed interpolation kernel.
};

/// Validated geometry and interpolation grid for one P4 search annulus.
/** Geometry is calculated in double precision. The transform controls only the stored kernel and sampled-image
 * arithmetic. `resize()` and `region()` are transactional: a rejected operation leaves the preceding valid state
 * unchanged. A successful resize clears any previously configured region.
 *
 * The prototype's half-pixel convention is preserved by constructing the canonical optimization wedge around
 * `(xCenter + 0.5, yCenter + 0.5)`, applying a 0.5-pixel radial enumeration buffer, and recording integer offsets
 * from its rounded-down wedge center. The buffered radial interval remains inner-inclusive and outer-exclusive, while
 * both nominal angular boundary rays are included symmetrically. All mapped samples must have their complete nominal
 * transform footprint inside the image.
 *
 * \ingroup programming_library
 */
template <typename _transformT>
class P4PixelGrid
{
  public:
    /// Interpolation transform type.
    using transformT = _transformT;

    /// Arithmetic type used by the transform and sampled images.
    using realT = typename transformT::arithT;

    /// Sampled image and common-mask storage type.
    using imageT = Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic>;

    /// One precomputed mapped interpolation record.
    using interpolationRecordT = P4InterpolationRecord<transformT>;

    /// Fixed-size interpolation kernel type.
    using kernelT = typename interpolationRecordT::kernelT;

    /// Nominal interpolation-kernel width.
    static constexpr int width = static_cast<int>( transformT::width );

    /// Number of pixels before the transform anchor in each dimension.
    static constexpr int leftBuffer = static_cast<int>( transformT::lbuff );

    /// Construct an empty grid.
    P4PixelGrid();

    /// Resize using the canonical image center `(0.5*(rows-1), 0.5*(columns-1))`.
    void resize( int rows, /**< [in] positive image row count */
                 int columns /**< [in] positive image column count */ );

    /// Resize using an explicit finite center inside the image.
    void resize( int rows,       /**< [in] positive image row count */
                 int columns,    /**< [in] positive image column count */
                 double xCenter, /**< [in] center row coordinate */
                 double yCenter /**< [in] center column coordinate */ );

    /// Build one search-annulus geometry and all interpolation kernels.
    /** A zero-valued common-mask element is invalid. Masked targets and predictor-footprint intersections mark the
     * complete local search-pixel fit invalid without changing the shared predictor count. Pass `nullptr` explicitly
     * when no common mask applies.
     */
    void region( const P4PixelGridRegion &configuration, /**< [in] complete region and exclusion configuration */
                 const imageT *commonMask /**< [in] optional image-wide common mask, or nullptr */ );

    /// Build only search pixels and pre-exclusion canonical predictor candidates.
    /** This path stops before central exclusion, mapped-coordinate edge checks, common-mask handling, and kernel
     * construction. It supports direct rotated sampling whose detector-frame footprints must be evaluated only after
     * inverse derotation. A successful candidate build is not reported as a complete `regionConfigured()` result.
     */
    void candidateRegion( const P4PixelGridRegion &configuration /**< [in] complete region configuration */ );

    /// Return the configured image row count.
    int rows() const noexcept;

    /// Return the configured image column count.
    int columns() const noexcept;

    /// Return the configured center row coordinate.
    double xCenter() const noexcept;

    /// Return the configured center column coordinate.
    double yCenter() const noexcept;

    /// Report whether a successful resize has configured image geometry.
    bool resized() const noexcept;

    /// Report whether a successful region build is available.
    bool regionConfigured() const noexcept;

    /// Report whether search pixels and pre-exclusion canonical candidates are available.
    bool candidateRegionConfigured() const noexcept;

    /// Return the most recently configured region.
    const P4PixelGridRegion &regionConfiguration() const;

    /// Return the effective exclusion radius `psfRadius + exclusionRadiusBuffer`.
    double effectiveExclusionRadius() const;

    /// Return the number of search pixels in the annulus.
    std::size_t searchPixelCount() const noexcept;

    /// Return the annulus-wide number of predictor columns.
    std::size_t predictorCount() const noexcept;

    /// Return the number of canonical wedge candidates before central exclusion is applied.
    std::size_t candidatePredictorCount() const noexcept;

    /// Return one immutable search-pixel coordinate and common-mask validity record.
    const P4SearchPixelRecord &searchPixel( std::size_t searchIndex /**< [in] zero-based search-pixel index */ ) const;

    /// Return one immutable canonical predictor offset.
    const P4PixelCoordinate &
    predictorOffset( std::size_t predictorIndex /**< [in] zero-based annulus-wide predictor index */ ) const;

    /// Return one canonical wedge candidate before central exclusion is applied.
    const P4PixelCoordinate &
    candidatePredictorOffset( std::size_t candidateIndex /**< [in] zero-based pre-exclusion candidate index */ ) const;

    /// Map one pre-exclusion canonical candidate into the fixed image frame at a search pixel.
    std::pair<double, double>
    candidateCoordinate( std::size_t searchIndex, /**< [in] zero-based search-pixel index */
                         std::size_t candidateIndex /**< [in] zero-based pre-exclusion candidate index */ ) const;

    /// Return one immutable mapped interpolation record.
    const interpolationRecordT &
    interpolation( std::size_t searchIndex, /**< [in] zero-based search-pixel index */
                   std::size_t predictorIndex /**< [in] zero-based annulus-wide predictor index */ ) const;

    /// Sample one mapped predictor from an image.
    /** Sampling an invalid local fit is rejected so a masked detector value cannot be consumed accidentally.
     */
    realT sample( const imageT &image,     /**< [in] image matching the configured grid dimensions */
                  std::size_t searchIndex, /**< [in] zero-based valid search-pixel index */
                  std::size_t predictorIndex /**< [in] zero-based predictor index */ ) const;

  private:
    /** \cond P4PixelGrid_test_harness */
    friend class P4PixelGridTestAccess;
    /** \endcond */

    /// Double-precision radius/angle working-image type.
    using geometryImageT = Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>;

    /// Convert a mapped coordinate into a checked nominal-footprint origin.
    static int checkedFootprintOrigin( double mappedCoordinate /**< [in] finite mapped coordinate */ );

    /// Return a checked search-major interpolation-record count.
    static std::size_t checkedInterpolationCount( std::size_t searchPixelCount, /**< [in] number of search pixels */
                                                  std::size_t predictorCount /**< [in] shared predictor count */ );

    /// Implement complete and candidate-only region construction transactionally.
    void regionImpl( const P4PixelGridRegion &configuration, /**< [in] complete region configuration */
                     const imageT *commonMask,               /**< [in] optional complete-region common mask */
                     bool candidatesOnly /**< [in] whether to stop before exclusion and interpolation */ );

    int m_rows{ 0 };              ///< Configured input-image row count.

    int m_columns{ 0 };           ///< Configured input-image column count.

    double m_xCenter{ 0 };        ///< Configured image-center row coordinate.

    double m_yCenter{ 0 };        ///< Configured image-center column coordinate.

    geometryImageT m_radiusImage; ///< Double-precision radius of every target-image pixel.

    geometryImageT m_angleImage;  ///< Double-precision position angle of every target-image pixel in degrees.

    std::optional<P4PixelGridRegion> m_region;          ///< Successfully constructed region, if any.

    bool m_regionComplete{ false };                     ///< Whether the current region includes interpolation data.

    std::vector<P4SearchPixelRecord> m_searchPixels;    ///< Search-pixel coordinates and validity records.

    std::vector<P4PixelCoordinate> m_candidateOffsets;  ///< Canonical wedge before central exclusion.

    std::vector<P4PixelCoordinate> m_predictorOffsets;  ///< Shared canonical predictor offsets.

    std::vector<interpolationRecordT> m_interpolations; ///< Search-major mapped interpolation records.
};

/// Supported production P4 pixel-grid specialization.
using P4PixelGridf = P4PixelGrid<mx::improc::cubicConvolTransform<float>>;

extern template class P4InterpolationRecord<mx::improc::cubicConvolTransform<float>>;
extern template class P4PixelGrid<mx::improc::cubicConvolTransform<float>>;

} // namespace improc
} // namespace mx

#endif // P4PixelGrid_hpp
