/** \file P4PSFFilter.hpp
 * \brief Declares normalized local filtering and product naming for spatially variable PSF responses.
 * \author Jared R. Males
 */

#ifndef P4PSFFilter_hpp
#define P4PSFFilter_hpp

#include <cstdint>
#include <string>

#include <Eigen/Dense>

namespace mx
{
namespace improc
{

/// One normalized local matched-filter result.
/** \ingroup programming_library */
struct P4PSFFilterResult
{
    double amplitude{ 0 };       ///< Estimated source amplitude relative to the stored PSF-template normalization.

    double correlation{ 0 };     ///< Signed response-science inner product before normalization.

    double normalization{ 0 };   ///< Retained signed-response squared norm.

    double supportFraction{ 0 }; ///< Fraction of stamp samples with usable PSF and science values.

    bool valid{ false };         ///< Whether support and normalization satisfy the filtering contract.
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
 * is the number of retained samples divided by the full stamp area.
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
};

} // namespace improc
} // namespace mx

#endif // P4PSFFilter_hpp
