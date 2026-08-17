/** \file P4RotatedGrid.hpp
 * \brief Declares direct sky-frame sampling geometry for rotated P4 regression.
 * \author Jared R. Males
 */

#ifndef P4RotatedGrid_hpp
#define P4RotatedGrid_hpp

#include "P4PixelGrid.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace mx
{
namespace improc
{

/** \cond P4RotatedGrid_test_harness */
class P4RotatedGridTestAccess;
/** \endcond */

/// Independent causes that invalidate one rotated P4 local fit.
/** Values form a bit mask because target/predictor and edge/mask failures may occur together.
 *
 * \ingroup programming_library
 */
enum class P4RotatedInvalidReason : std::uint8_t
{
    none = 0,                      ///< Every required detector-frame interpolation footprint is valid.
    targetFootprintOutside = 1,    ///< A target footprint crosses an unrotated detector-image edge.
    predictorFootprintOutside = 2, ///< A predictor footprint crosses an unrotated detector-image edge.
    targetFootprintMasked = 4,     ///< A target footprint touches a zero-valued detector common-mask pixel.
    predictorFootprintMasked = 8   ///< A predictor footprint touches a zero-valued detector common-mask pixel.
};

/// Immutable all-frame validity record for one rotated P4 search pixel.
/** \ingroup programming_library */
class P4RotatedSearchRecord
{
  public:
    /// Construct one complete all-frame validity record.
    P4RotatedSearchRecord( const P4PixelCoordinate &coordinate, /**< [in] fixed sky-frame target coordinate */
                           P4RotatedInvalidReason invalidReason /**< [in] combined invalidity bit mask */ );

    /// Return the fixed sky-frame target coordinate.
    const P4PixelCoordinate &coordinate() const noexcept;

    /// Return the combined invalidity bit mask.
    P4RotatedInvalidReason invalidReason() const noexcept;

    /// Report whether the local fit has complete support in every configured frame.
    bool valid() const noexcept;

    /// Report whether one specific invalidity bit is present.
    bool hasReason( P4RotatedInvalidReason reason /**< [in] one non-combined invalidity bit */ ) const noexcept;

  private:
    P4PixelCoordinate m_coordinate;         ///< Fixed sky-frame target coordinate.

    P4RotatedInvalidReason m_invalidReason; ///< Combined all-frame invalidity bit mask.
};

/// One on-demand direct interpolation record in a preprocessed, unrotated detector frame.
/** \ingroup programming_library */
class P4RotatedInterpolationRecord
{
  public:
    /// Fixed float cubic-convolution kernel type.
    using kernelT = P4PixelGridf::kernelT;

    /// Return the mapped unrotated detector-image row coordinate.
    double row() const noexcept;

    /// Return the mapped unrotated detector-image column coordinate.
    double column() const noexcept;

    /// Return the first row consumed by the nominal cubic footprint.
    int footprintRow() const noexcept;

    /// Return the first column consumed by the nominal cubic footprint.
    int footprintColumn() const noexcept;

    /// Return the direct cubic-convolution kernel.
    const kernelT &kernel() const noexcept;

    /// Report whether the complete nominal footprint is inside the configured detector image.
    bool footprintInside() const noexcept;

  private:
    friend class P4RotatedGrid;

    /// Construct one complete on-demand interpolation record.
    P4RotatedInterpolationRecord( double row,            /**< [in] unrotated detector-image row coordinate */
                                  double column,         /**< [in] unrotated detector-image column coordinate */
                                  int footprintRow,      /**< [in] nominal-footprint first row */
                                  int footprintColumn,   /**< [in] nominal-footprint first column */
                                  const kernelT &kernel, /**< [in] direct cubic kernel */
                                  bool footprintInside /**< [in] whether full nominal support is inside */ );

    double m_row;           ///< Mapped unrotated detector-image row coordinate.

    double m_column;        ///< Mapped unrotated detector-image column coordinate.

    int m_footprintRow;     ///< First nominal-footprint row.

    int m_footprintColumn;  ///< First nominal-footprint column.

    kernelT m_kernel;       ///< Direct cubic-convolution kernel.

    bool m_footprintInside; ///< Whether full nominal support is inside the detector image.
};

/// Direct one-stage detector-image sampler for a fixed sky-frame P4 predictor wedge.
/** `configure()` copies fixed target and predictor coordinates from a candidate-configured `P4PixelGridf`, maps those
 * sky coordinates through every inverse derotation, and records one all-frame validity result per search pixel. The
 * detector-frame common mask is applied to every pixel of every nominal target and predictor footprint; mask validity
 * previously attached to a complete source sky grid is deliberately ignored. Use `P4PixelGridf::candidateRegion()`
 * so pre-exclusion candidates reach this direct policy without prior sky-frame interpolation checks.
 *
 * Detector-frame interpolation records and kernels are generated on demand. Storage is proportional to `M*K + T`,
 * not `M*T*K`. Under `kernelSupport`, a source predictor column is retained only when its exact inverse-rotated
 * detector footprint avoids the configured exclusion disk for every search pixel and frame. This filtering remains
 * annulus-wide, so every retained local fit has the same `K`.
 *
 * The input image for sampling is the preprocessed but unrotated detector image used by P4 regression, not the
 * pristine observation before inherited preprocessing.
 *
 * `configure()` is transactional: rejected input leaves the preceding successful configuration unchanged.
 *
 * \ingroup programming_library
 */
class P4RotatedGrid
{
  public:
    /// Image storage accepted by the production float implementation.
    using imageT = P4PixelGridf::imageT;

    /// Arithmetic type used by the direct cubic transform.
    using realT = P4PixelGridf::realT;

    /// On-demand unrotated detector-frame interpolation record type.
    using interpolationRecordT = P4RotatedInterpolationRecord;

    /// Nominal direct interpolation-kernel width.
    static constexpr int width = P4PixelGridf::width;

    /// Number of pixels before the direct transform anchor.
    static constexpr int leftBuffer = P4PixelGridf::leftBuffer;

    /// Construct an empty rotated grid.
    P4RotatedGrid();

    /// Configure fixed sky geometry, inverse derotations, and all-frame detector-footprint validity.
    /** Derotation angles use the same counter-clockwise radians convention as `mx::improc::imageRotate`. At least
     * one finite angle is required. A zero-valued common-mask element invalidates every local fit whose target or
     * predictor nominal footprint touches it. Pass `nullptr` explicitly when no detector common mask applies.
     */
    void configure( const P4PixelGridf &skyGrid,                 /**< [in] configured fixed sky-frame P4 geometry */
                    const std::vector<double> &derotationAngles, /**< [in] one radians angle per unrotated frame */
                    const imageT *rawCommonMask /**< [in] optional detector-frame common mask, or nullptr */ );

    /// Report whether a successful configuration is available.
    bool configured() const noexcept;

    /// Return the configured unrotated detector-image row count.
    int rows() const noexcept;

    /// Return the configured unrotated detector-image column count.
    int columns() const noexcept;

    /// Return the configured rotation-center row coordinate.
    double xCenter() const noexcept;

    /// Return the configured rotation-center column coordinate.
    double yCenter() const noexcept;

    /// Return the configured preprocessed, unrotated detector-frame count `T`.
    std::size_t frameCount() const noexcept;

    /// Return the fixed sky-frame search-pixel count `M`.
    std::size_t searchPixelCount() const noexcept;

    /// Return the retained annulus-wide predictor count `K`.
    std::size_t predictorCount() const noexcept;

    /// Return one immutable all-frame search-pixel validity record.
    const P4RotatedSearchRecord &
    searchPixel( std::size_t searchIndex /**< [in] zero-based search-pixel index */ ) const;

    /// Return one retained canonical sky-frame predictor offset.
    const P4PixelCoordinate &
    predictorOffset( std::size_t predictorIndex /**< [in] zero-based retained predictor index */ ) const;

    /// Return the pre-exclusion candidate index for one retained column.
    std::size_t
    candidatePredictorIndex( std::size_t predictorIndex /**< [in] zero-based retained predictor index */ ) const;

    /// Construct the direct detector-frame interpolation record for one target sample.
    interpolationRecordT
    targetInterpolation( std::size_t frameIndex, /**< [in] zero-based unrotated detector-frame index */
                         std::size_t searchIndex /**< [in] zero-based search-pixel index */ ) const;

    /// Construct the direct detector-frame interpolation record for one predictor sample.
    interpolationRecordT
    predictorInterpolation( std::size_t frameIndex,  /**< [in] zero-based unrotated detector-frame index */
                            std::size_t searchIndex, /**< [in] zero-based search-pixel index */
                            std::size_t predictorIndex /**< [in] zero-based retained predictor index */ ) const;

    /// Directly sample one target from its preprocessed, unrotated detector image.
    realT sampleTarget( const imageT &rawImage, /**< [in] preprocessed detector image matching the grid */
                        std::size_t frameIndex, /**< [in] zero-based unrotated detector-frame index */
                        std::size_t searchIndex /**< [in] zero-based all-frame-valid search index */ ) const;

    /// Directly sample one predictor from its preprocessed, unrotated detector image.
    realT samplePredictor( const imageT &rawImage,  /**< [in] preprocessed detector image matching the grid */
                           std::size_t frameIndex,  /**< [in] zero-based unrotated detector-frame index */
                           std::size_t searchIndex, /**< [in] zero-based all-frame-valid search index */
                           std::size_t predictorIndex /**< [in] zero-based retained predictor index */ ) const;

    /// Transactionally sample one target and the complete retained predictor row for a frame.
    /** `predictors` must already have exactly `predictorCount()` elements. Neither output is changed if any required
     * detector-image input value or interpolated result is non-finite.
     */
    void sampleFrame( const imageT &rawImage,  /**< [in] preprocessed detector image matching the grid */
                      std::size_t frameIndex,  /**< [in] zero-based unrotated detector-frame index */
                      std::size_t searchIndex, /**< [in] zero-based all-frame-valid search index */
                      realT &target,           /**< [out] direct target value */
                      std::vector<realT> &predictors /**< [out] complete retained predictor row */ ) const;

  private:
    /** \cond P4RotatedGrid_test_harness */
    friend class P4RotatedGridTestAccess;
    /** \endcond */

    /// Cached inverse-rotation trigonometric pair for one frame.
    struct FrameRotation
    {
        double cosine; ///< Cosine of the configured derotation angle.

        double sine;   ///< Sine of the configured derotation angle.
    };

    /// Fixed mapped sky coordinate for one predictor sample.
    struct SkyCoordinate
    {
        double row;    ///< Fixed sky-frame row coordinate.

        double column; ///< Fixed sky-frame column coordinate.
    };

    /// Lightweight detector-frame mapping used while validating all required footprints.
    struct MappedFootprint
    {
        double row;          ///< Mapped unrotated detector-image row coordinate.

        double column;       ///< Mapped unrotated detector-image column coordinate.

        int footprintRow;    ///< First nominal-footprint row.

        int footprintColumn; ///< First nominal-footprint column.

        bool inside;         ///< Whether full nominal support is inside the detector image.
    };

    /// Construct an on-demand record for one fixed sky coordinate.
    interpolationRecordT interpolation( std::size_t frameIndex, /**< [in] zero-based frame index */
                                        const SkyCoordinate &skyCoordinate /**< [in] fixed sky coordinate */ ) const;

    /// Construct one direct interpolation record from explicit transactional geometry.
    static interpolationRecordT makeInterpolation( int rows,                      /**< [in] detector-image rows */
                                                   int columns,                   /**< [in] detector-image columns */
                                                   double xCenter,                /**< [in] center row */
                                                   double yCenter,                /**< [in] center column */
                                                   const FrameRotation &rotation, /**< [in] frame trigonometry */
                                                   const SkyCoordinate &skyCoordinate /**< [in] fixed sky point */ );

    /// Map one fixed sky coordinate without constructing an interpolation kernel.
    static MappedFootprint mapFootprint( int rows,                      /**< [in] detector-image rows */
                                         int columns,                   /**< [in] detector-image columns */
                                         double xCenter,                /**< [in] center row */
                                         double yCenter,                /**< [in] center column */
                                         const FrameRotation &rotation, /**< [in] frame trigonometry */
                                         const SkyCoordinate &skyCoordinate /**< [in] fixed sky point */ );

    /// Return a checked search-major fixed sky-coordinate count.
    static std::size_t checkedCoordinateCount( std::size_t searchCount, /**< [in] fixed target count */
                                               std::size_t predictorCount /**< [in] retained predictor count */ );

    /// Return whether one complete detector footprint is nonzero in the optional common mask.
    static bool maskFootprintValid( const imageT *rawCommonMask, /**< [in] optional validated common mask */
                                    const MappedFootprint &footprint /**< [in] inside raw footprint */ );

    /// Return whether a detector predictor footprint overlaps the configured exclusion disk.
    static bool footprintOverlapsExclusion( const MappedFootprint &predictor, /**< [in] predictor footprint */
                                            const MappedFootprint &target,    /**< [in] target sample center */
                                            double effectiveRadius /**< [in] exclusion-disk radius */ );

    /// Apply one direct record while requiring finite source and result values.
    static realT sampleRecord( const imageT &rawImage, /**< [in] dimension-validated detector image */
                               const interpolationRecordT &record /**< [in] complete inside footprint */ );

    /// Validate common sampling preconditions for one local fit and frame.
    void requireSampleable( const imageT &rawImage, /**< [in] detector image to validate */
                            std::size_t frameIndex, /**< [in] zero-based unrotated detector-frame index */
                            std::size_t searchIndex /**< [in] zero-based search-pixel index */ ) const;

    /// Return one search-major fixed sky predictor coordinate.
    const SkyCoordinate &skyPredictor( std::size_t searchIndex, /**< [in] zero-based search-pixel index */
                                       std::size_t predictorIndex /**< [in] zero-based predictor index */ ) const;

    int m_rows{ 0 };                                      ///< Configured detector-image row count.

    int m_columns{ 0 };                                   ///< Configured detector-image column count.

    double m_xCenter{ 0 };                                ///< Configured common rotation-center row.

    double m_yCenter{ 0 };                                ///< Configured common rotation-center column.

    std::vector<FrameRotation> m_frameRotations;          ///< Per-frame inverse-rotation trigonometry.

    std::vector<P4RotatedSearchRecord> m_searchPixels;    ///< Fixed targets and all-frame validity.

    std::vector<P4PixelCoordinate> m_predictorOffsets;    ///< Retained canonical predictor offsets.

    std::vector<std::size_t> m_candidatePredictorIndices; ///< Retained pre-exclusion candidate indices.

    std::vector<SkyCoordinate> m_skyPredictors;           ///< Search-major fixed sky predictor coordinates.
};

} // namespace improc
} // namespace mx

#endif // P4RotatedGrid_hpp
