/** \file P4LocalProcessing.hpp
 * \brief Declares sparse geometry and trial-source sampling for pixel-local P4 processing.
 * \author Jared R. Males
 */

#ifndef P4LocalProcessing_hpp
#define P4LocalProcessing_hpp

#include <cstddef>
#include <cstdint>
#include <vector>

#include <Eigen/Core>

#include "P4PixelGrid.hpp"

namespace mx
{
namespace improc
{

/** \cond P4LocalProcessing_test_harness */
class P4LocalProcessingTestAccess;
/** \endcond */

/// One detector search pixel and the target frames needed by a local sky stamp.
/** \ingroup programming_library */
class P4LocalSearchRequest
{
  public:
    /// Construct one immutable sparse search request.
    P4LocalSearchRequest( const P4PixelCoordinate &coordinate, /**< [in] detector search-pixel coordinate */
                          int region,                          /**< [in] owning annulus index */
                          std::size_t searchIndex,             /**< [in] annulus-local search index */
                          std::vector<int> frames /**< [in] sorted unique target-frame indices */ );

    /// Return the detector search-pixel coordinate.
    const P4PixelCoordinate &coordinate() const noexcept;

    /// Return the owning annulus index.
    int region() const noexcept;

    /// Return the annulus-local search index.
    std::size_t searchIndex() const noexcept;

    /// Return the sorted unique target-frame indices needed from this fit.
    const std::vector<int> &frames() const noexcept;

  private:
    P4PixelCoordinate m_coordinate; ///< Detector search-pixel coordinate.

    int m_region{ -1 };             ///< Owning annulus index.

    std::size_t m_searchIndex{ 0 }; ///< Annulus-local search index.

    std::vector<int> m_frames;      ///< Sorted unique requested target frames.
};

/// One weighted sparse residual sample used by local derotation.
/** \ingroup programming_library */
class P4LocalResidualSample
{
  public:
    /// Construct one immutable sparse residual sample.
    P4LocalResidualSample( std::size_t requestIndex, /**< [in] index into the local search-request vector */
                           std::size_t frameOffset,  /**< [in] requested frame's offset within that search request */
                           float weight /**< [in] cubic derotation weight */ );

    /// Return the local search-request index.
    std::size_t requestIndex() const noexcept;

    /// Return the requested frame's offset within its search request.
    std::size_t frameOffset() const noexcept;

    /// Return the cubic derotation weight.
    float weight() const noexcept;

  private:
    std::size_t m_requestIndex{ 0 }; ///< Local search-request index.

    std::size_t m_frameOffset{ 0 };  ///< Frame offset within the search request.

    float m_weight{ 0 };             ///< Cubic derotation weight.
};

/// Sparse detector footprint for one output-stamp pixel in one target frame.
/** \ingroup programming_library */
class P4LocalOutputSample
{
  public:
    /// Construct one invalid output sample with no residual dependencies.
    P4LocalOutputSample();

    /// Construct one valid output sample from its ordered residual dependencies.
    explicit P4LocalOutputSample( std::vector<P4LocalResidualSample> samples
                                  /**< [in] ordered one- or sixteen-pixel residual footprint */ );

    /// Report whether the complete detector footprint belongs to configured P4 search annuli.
    bool valid() const noexcept;

    /// Return the ordered sparse residual dependencies.
    const std::vector<P4LocalResidualSample> &samples() const noexcept;

  private:
    bool m_valid{ false };                        ///< Whether complete detector geometry is available.

    std::vector<P4LocalResidualSample> m_samples; ///< Ordered sparse detector residual dependencies.
};

/// Exact sparse derotation dependency map for one square local sky stamp.
/** The map reproduces the detector coordinates, nominal cubic footprint, and interpolation kernel used by
 * `mx::improc::imageRotate`. Each unique detector search pixel is represented once, with only the target frames
 * requested by the local output stamp.
 *
 * \ingroup programming_library
 */
class P4LocalGeometry
{
  public:
    /// Signed integer image used for annulus and search-index lookup.
    using lookupImageT = Eigen::Array<std::int64_t, Eigen::Dynamic, Eigen::Dynamic>;

    /// Construct an unconfigured sparse geometry.
    P4LocalGeometry();

    /// Configure one local output stamp transactionally.
    void configure( int rows,                          /**< [in] full detector-image row count */
                    int columns,                       /**< [in] full detector-image column count */
                    int stampSize,                     /**< [in] positive odd local sky-stamp width */
                    double sourceRow,                  /**< [in] finite continuous final-sky source row */
                    double sourceColumn,               /**< [in] finite continuous final-sky source column */
                    const std::vector<double> &angles, /**< [in] one finite derotation angle per target frame */
                    bool derotate,                     /**< [in] whether nonzero angles rotate residual frames */
                    const lookupImageT &ownership,     /**< [in] detector coordinate to owning annulus */
                    const lookupImageT &searchIndexLookup /**< [in] detector coordinate to annulus-local index */ );

    /// Report whether a successful configuration is available.
    bool configured() const noexcept;

    /// Return the requested square output-stamp width.
    int stampSize() const noexcept;

    /// Return the integer full-image row of local stamp element `(0,0)`.
    int originRow() const noexcept;

    /// Return the integer full-image column of local stamp element `(0,0)`.
    int originColumn() const noexcept;

    /// Return the continuous final-sky trial-source row.
    double sourceRow() const noexcept;

    /// Return the continuous final-sky trial-source column.
    double sourceColumn() const noexcept;

    /// Return the configured target-frame count.
    std::size_t frameCount() const noexcept;

    /// Return every unique detector search request in deterministic coordinate order.
    const std::vector<P4LocalSearchRequest> &searchRequests() const noexcept;

    /// Return one output pixel/frame sparse dependency record.
    const P4LocalOutputSample &outputSample( int stampRow,    /**< [in] zero-based local stamp row */
                                             int stampColumn, /**< [in] zero-based local stamp column */
                                             std::size_t frame /**< [in] zero-based target frame */ ) const;

    /// Return the number of requested sparse `(detector pixel, target frame)` residual samples.
    std::size_t sparseSampleCount() const noexcept;

    /// Return the retained dynamic storage used by sparse requests and output dependencies.
    std::size_t storageBytes() const;

  private:
    /** \cond P4LocalProcessing_test_harness */
    friend class P4LocalProcessingTestAccess;
    /** \endcond */

    /// Return the checked column-major output-record offset.
    std::size_t outputOffset( int stampRow,    /**< [in] zero-based local stamp row */
                              int stampColumn, /**< [in] zero-based local stamp column */
                              std::size_t frame /**< [in] zero-based target frame */ ) const;

    int m_stampSize{ 0 };                               ///< Requested output-stamp width, or zero before configuration.

    int m_originRow{ 0 };                               ///< Full-image row occupied by local stamp element `(0,0)`.

    int m_originColumn{ 0 };                            ///< Full-image column occupied by local stamp element `(0,0)`.

    double m_sourceRow{ 0 };                            ///< Continuous final-sky trial-source row.

    double m_sourceColumn{ 0 };                         ///< Continuous final-sky trial-source column.

    std::size_t m_frameCount{ 0 };                      ///< Configured target-frame count.

    std::size_t m_sparseSampleCount{ 0 };               ///< Unique requested detector-pixel/frame pairs.

    std::vector<P4LocalSearchRequest> m_searchRequests; ///< Unique detector fits needed by the stamp.

    std::vector<P4LocalOutputSample> m_outputSamples;   ///< Column-major stamp/frame dependency records.
};

/// Prepared centered trial source sampled without materializing shifted full-frame images.
/** The configured template is center-cropped and zero-padded to the detector dimensions. When the requested odd
 * stamp and detector/template sampling phases differ, the internal crop grows by one pixel to preserve the geometric
 * center. `value()` reproduces the ordinary `imageShift` whole-pixel or cubic-convolution value at one detector
 * pixel, then applies the configured per-frame scale and contrast.
 *
 * \ingroup programming_library
 */
class P4TrialSource
{
  public:
    /// Float image type used by production P4.
    using imageT = P4PixelGridf::imageT;

    /// Construct an unconfigured trial source.
    P4TrialSource();

    /// Configure one immutable trial source transactionally.
    void configure( const imageT &sourceTemplate,      /**< [in] finite centered source template */
                    int detectorRows,                  /**< [in] target detector-image row count */
                    int detectorColumns,               /**< [in] target detector-image column count */
                    int requestedStampSize,            /**< [in] positive nominal centered crop width */
                    const std::vector<double> &angles, /**< [in] one finite derotation angle per target frame */
                    double separation,                 /**< [in] finite nonnegative sky separation in pixels */
                    double positionAngle,              /**< [in] finite PA in degrees east of north */
                    double contrast,                   /**< [in] finite trial contrast, including zero */
                    const std::vector<float> &scales /**< [in] one finite per-frame flux scale */ );

    /// Add another source using the already-prepared template and detector geometry.
    void addSource( const std::vector<double> &angles, /**< [in] one finite derotation angle per target frame */
                    double separation,                 /**< [in] finite nonnegative sky separation in pixels */
                    double positionAngle,              /**< [in] finite PA in degrees east of north */
                    double contrast,                   /**< [in] finite signed contrast */
                    const std::vector<float> &scales /**< [in] one finite per-frame flux scale */ );

    /// Report whether a successful configuration is available.
    bool configured() const noexcept;

    /// Return the actual phase-preserving internal crop row count.
    int cropRows() const noexcept;

    /// Return the actual phase-preserving internal crop column count.
    int cropColumns() const noexcept;

    /// Return the prepared trial perturbation at one integer detector pixel.
    float value( std::size_t frame, /**< [in] zero-based target frame */
                 int row,           /**< [in] integer detector row */
                 int column /**< [in] integer detector column */ ) const;

    /// Sample the prepared perturbation through one production P4 interpolation record.
    float sample( std::size_t frame, /**< [in] zero-based target frame */
                  const P4PixelGridf::interpolationRecordT &record
                  /**< [in] production predictor interpolation record */ ) const;

  private:
    /** \cond P4LocalProcessing_test_harness */
    friend class P4LocalProcessingTestAccess;
    /** \endcond */

    /// One precomputed physical-frame source shift and cubic kernel.
    struct FrameShift
    {
        float rowShift{ 0 };          ///< Detector-row translation in production float arithmetic.

        float columnShift{ 0 };       ///< Detector-column translation in production float arithmetic.

        float scale{ 0 };             ///< Product of trial contrast and per-frame scale.

        bool integral{ false };       ///< Whether both translations are exact integers.

        P4PixelGridf::kernelT kernel; ///< Cubic shift kernel when translation is fractional.
    };

    imageT m_template;                ///< Center-cropped template zero-padded to detector dimensions.

    int m_cropRows{ 0 };              ///< Actual phase-preserving internal crop rows.

    int m_cropColumns{ 0 };           ///< Actual phase-preserving internal crop columns.

    std::vector<std::vector<FrameShift>> m_sourceShifts; ///< Per-source prepared shifts and scales for each frame.
};

} // namespace improc
} // namespace mx

#endif // P4LocalProcessing_hpp
