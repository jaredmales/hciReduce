/** \file P4PSFModel.hpp
 * \brief Declares compact frozen-model PSF calculations for Pixel Prediction Post-Processing.
 * \author Jared R. Males
 */

#ifndef P4PSFModel_hpp
#define P4PSFModel_hpp

#include <cstddef>
#include <cstdint>

#include <Eigen/Dense>

#include "P4PixelGrid.hpp"

namespace mx
{
namespace improc
{

/// Pure local frozen-model PSF calculation for detector-frame P4.
/** The supplied template is interpreted as a detector-sampled, post-preprocessing PSF centered at its geometric
 * `(rows-1)/2,(columns-1)/2` coordinate. Values outside the template are zero. The returned stamp represents the
 * residual response as the source moves about one fixed detector search pixel; its geometric center is zero source
 * offset. Predictor sampling composes template shifting with the exact stored P4 cubic interpolation kernels.
 *
 * \ingroup programming_library
 */
class P4PSFModel
{
  public:
    /// Float image and PSF-stamp storage.
    using imageT = Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic>;

    /// Double-precision predictor coefficient vector.
    using coefficientT = Eigen::Array<double, Eigen::Dynamic, 1>;

    /// Detector-frame P4 geometry using the production float interpolation transform.
    using gridT = P4PixelGridf;

    /// Construct and validate a centered template for one square response-stamp size.
    /** The finite template is sampled once onto the output stamp's fixed subpixel phase. Values beyond its supplied
     * extent remain zero. The resulting object is immutable and safe to share among concurrent local calculations.
     */
    P4PSFModel( const imageT &psfTemplate, /**< [in] finite centered detector-sampled PSF */
                int stampSize /**< [in] positive output width and height */ );

    /// Construct and validate a centered template for a rectangular response-stamp grid.
    /** Matching each stamp dimension's parity to the corresponding template dimension preserves the template's
     * native integer-sample phase when the response is later interpolated at fractional offsets.
     */
    P4PSFModel( const imageT &psfTemplate, /**< [in] finite centered detector-sampled PSF */
                int stampRows,             /**< [in] positive output row count */
                int stampColumns /**< [in] positive output column count */ );

    /// Return the configured response-stamp row count.
    int stampRows() const noexcept;

    /// Return the configured response-stamp column count.
    int stampColumns() const noexcept;

    /// Return the exact byte count of the retained precomputed float template samples.
    std::size_t storageBytes() const noexcept;

    /// Calculate one compact local residual PSF stamp from a frozen coefficient vector.
    /** The grid must contain a complete valid detector-frame region, \p searchIndex must select a valid local fit,
     * and the coefficient count must equal the grid predictor count. On success, \p output is completely replaced;
     * on exception its contents are unspecified.
     */
    void calculateLocalResponse( imageT &output,          /**< [out] configured float residual PSF stamp */
                                 const gridT &grid,       /**< [in] complete detector-frame P4 geometry */
                                 std::size_t searchIndex, /**< [in] zero-based search-pixel index */
                                 const Eigen::Ref<const coefficientT> &coefficients
                                 /**< [in] finite predictor coefficients */ ) const;

  private:
    /// Return one phase-shifted template sample or zero outside the stored support.
    float shiftedTemplateValue( std::int64_t rowIndex, /**< [in] integer response-grid row */
                                std::int64_t columnIndex /**< [in] integer response-grid column */ ) const noexcept;

    imageT m_shiftedTemplate;               ///< Template sampled at the output stamp's fixed subpixel phase.

    std::int64_t m_minimumRowIndex{ 0 };    ///< Response-grid row represented by stored row zero.

    std::int64_t m_minimumColumnIndex{ 0 }; ///< Response-grid column represented by stored column zero.

    int m_stampRows{ 0 };                   ///< Positive response-stamp row count.

    int m_stampColumns{ 0 };                ///< Positive response-stamp column count.
};

} // namespace improc
} // namespace mx

#endif // P4PSFModel_hpp
