/** \file P4PSFFilter.hpp
 * \brief Declares normalized local filtering and product naming for spatially variable PSF responses.
 */

#ifndef P4PSFFilter_hpp
#define P4PSFFilter_hpp

#include <cstdint>
#include <string>

#include <Eigen/Dense>

#include "PSFNoiseModel.hpp"

namespace mx
{
namespace improc
{

/// One normalized local matched-filter result.
/** \ingroup programming_library */
struct P4PSFFilterResult
{
    double amplitude{ 0 };        ///< Estimated source amplitude relative to the stored PSF-template normalization.

    double correlation{ 0 };      ///< Signed weighted correlation t^T C^-1 (science - mean).

    double normalization{ 0 };    ///< Retained weighted response norm t^T C^-1 t.

    double supportFraction{ 0 };  ///< Fraction of stamp samples with usable PSF and science values.

    bool valid{ false };          ///< Whether support and normalization satisfy the filtering contract.

    double conditionalSigma{ 0 }; ///< Amplitude uncertainty conditional on the supplied covariance and template.

    double score{ 0 }; ///< Signed correlation divided by sqrt(normalization); requires empirical detection calibration.

    double nonnegativeLogLikelihood{ 0 }; ///< Log-likelihood improvement for a nonnegative source, max(0,score)^2 / 2.
};

/// Derive a filtered-image product path from a resolved final-image path.
std::string psfFilterProductPath( const std::string &finalImagePath, /**< [in] resolved final-image output path */
                                  const std::string &role, /**< [in] filename role inserted before the sequence */
                                  bool sequential /**< [in] whether the final image uses a four-digit sequence */ );

/// Derive a filter-diagnostic path in the final image's auxiliary-product directory.
std::string psfFilterDiagnosticPath( const std::string &finalImagePath, /**< [in] resolved final-image output path */
                                     const std::string &role, /**< [in] filename role inserted before the sequence */
                                     bool sequential /**< [in] whether the final image uses a four-digit sequence */ );

/// Apply one signed local PSF response to its corresponding final-image neighborhood.
/** Filtering uses `sum(H * I) / sum(H * H)` over response-valid, in-bounds, finite science samples. Stamp dimensions
 * must be odd so the geometric response center maps exactly onto the integer source coordinate. The support fraction
 * is the number of retained samples divided by the full stamp area. The weighted calculation instead solves
 * C w = t on those retained samples and returns w^T (science - mean) / w^T t. Its conditional uncertainty and score
 * do not account for covariance/template estimation error or dependence between overlapping training stamps.
 *
 * \ingroup programming_library
 */
class P4PSFFilter
{
  public:
    /// Float image type used by science and response products.
    using imageT = Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic>;

    /// Read-only view of a science or response image without copying cube planes.
    using imageConstRefT = Eigen::Ref<const imageT>;

    /// Byte-valued per-response-sample validity type.
    using validityT = Eigen::Array<std::uint8_t, Eigen::Dynamic, Eigen::Dynamic>;

    /// Read-only view of response validity without copying mapped arrays.
    using validityConstRefT = Eigen::Ref<const validityT>;

    /// Calculate one normalized local matched-filter result.
    static P4PSFFilterResult
    calculate( imageConstRefT science,     /**< [in] final science image for one mode */
               imageConstRefT response,    /**< [in] signed local response stamp */
               validityConstRefT validity, /**< [in] valid response samples */
               int sourceRow,              /**< [in] integer source-center image row */
               int sourceColumn,           /**< [in] integer source-center image column */
               double minimumSupportFraction /**< [in] required usable fraction in `[0,1]` */ );

    /// Calculate a covariance-weighted result using an independent residual-noise model.
    /** An identity model delegates to the original calculation. Otherwise the marginal covariance is restricted
     * to the current support before solving; the mean is subtracted from science only. The caller is responsible
     * for matching the model's column-major pixel coordinates to the response stamp and excluding source footprints
     * from training. Malformed inputs or a failed covariance solve throw; insufficient support or nonfinite filter
     * statistics produce an invalid result. Diagnostic values are meaningful only when valid is true.
     */
    static P4PSFFilterResult calculateWeighted(
        imageConstRefT science,        /**< [in] final science image for one mode */
        imageConstRefT response,       /**< [in] signed local response stamp */
        validityConstRefT validity,    /**< [in] valid response samples */
        int sourceRow,                 /**< [in] integer source-center image row */
        int sourceColumn,              /**< [in] integer source-center image column */
        double minimumSupportFraction, /**< [in] required usable fraction in `[0,1]` */
        const PSFNoiseModel &noise /**< [in] covariance and background mean on the full response lattice */ );
};

} // namespace improc
} // namespace mx

#endif // P4PSFFilter_hpp
