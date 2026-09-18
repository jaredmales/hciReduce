/** \file PSFNoiseModel.hpp
 * \brief Declares positive-definite residual-noise models for local PSF filtering.
 */

#ifndef PSFNoiseModel_hpp
#define PSFNoiseModel_hpp

#include <cstddef>
#include <vector>

#include <Eigen/Dense>

namespace mx
{
namespace improc
{

/// Covariance representation used to weight a local response stamp.
enum class PSFNoiseKind
{
    identity, ///< Unit independent noise and zero background mean.
    diagonal, ///< Independent noise with a positive variance for each stamp pixel.
    pca       ///< Positive diagonal floor plus a finite low-rank covariance factor.
};

/// Immutable covariance and background mean on one complete response-stamp lattice.
/** Pixels use column-major stamp ordering: index = row + column * stampRows. A low-rank model represents
 * C = D + L L^T, with strictly positive D. Missing pixels are marginalized by restricting D and the rows of L
 * before solving; selecting entries from the full precision matrix would give a different statistic.
 * Training samples must already share the science/template coordinates and exclude candidate/source footprints.
 * Sample-count normalization does not correct dependence between overlapping training patches.
 *
 * \ingroup programming_library
 */
class PSFNoiseModel
{
  public:
    /// Double-precision vector in stamp-pixel order.
    using vectorT = Eigen::VectorXd;

    /// Double-precision training samples or low-rank covariance factor.
    using matrixT = Eigen::MatrixXd;

  private:
    /// \name Covariance Model - Data
    ///@{
    PSFNoiseKind m_kind; ///< Selected representation, retained even when its low-rank factor is empty.

    vectorT m_variances; ///< Strictly positive diagonal covariance floor in full stamp order.

    matrixT m_factor;    ///< Owned factor L; rows follow stamp pixels, columns represent noise modes.

    vectorT m_mean;      ///< Owned background mean to subtract from science, never from the response template.
    ///@}

    /// Validate and own one covariance representation.
    PSFNoiseModel( PSFNoiseKind kind, /**< [in] selected covariance representation */
                   vectorT variances, /**< [in] positive diagonal covariance floor */
                   matrixT factor,    /**< [in] low-rank factor with one row per stamp pixel */
                   vectorT mean /**< [in] background mean; empty selects zero */ );

  public:
    /// \name Covariance Model
    ///@{
    /// Construct unit independent noise with zero background mean.
    static PSFNoiseModel identity( Eigen::Index pixels /**< [in] positive full stamp-pixel count */ );

    /// Construct independent noise from per-pixel variances and an optional known mean.
    static PSFNoiseModel
    diagonal( const vectorT &variances, /**< [in] finite strictly positive pixel variances */
              const vectorT &mean = vectorT{} /**< [in] finite background mean; empty selects zero */ );

    /// Construct C = D + L L^T without requiring orthogonal columns of L.
    static PSFNoiseModel
    lowRank( const vectorT &variances, /**< [in] finite strictly positive diagonal floor D */
             const matrixT &factor,    /**< [in] finite factor L with one row per stamp pixel */
             const vectorT &mean = vectorT{} /**< [in] finite background mean; empty selects zero */ );

    /// Estimate a mean and sample variances, bounded below by an absolute variance floor.
    static PSFNoiseModel
    estimateDiagonal( const matrixT &samples, /**< [in] at least two finite noise-only rows */
                      double varianceFloor /**< [in] finite positive variance floor in image units squared */ );

    /// Estimate centered PCA modes and their excess variances above an isotropic floor.
    /** For sample eigenvalues nu_i, retain at most maximumModes factors sqrt(max(nu_i-floor,0))*u_i.
     * Covariance uses n-1 normalization. Unretained directions keep the positive floor; they are not projected out.
     */
    static PSFNoiseModel
    estimatePCA( const matrixT &samples,   /**< [in] finite noise-only rows in common stamp coordinates */
                 std::size_t maximumModes, /**< [in] maximum retained modes; zero selects the floor only */
                 double varianceFloor /**< [in] finite positive isotropic floor in image units squared */ );

    /// Return the selected covariance representation.
    PSFNoiseKind kind() const;

    /// Return the number of pixels in the complete stamp lattice.
    Eigen::Index pixels() const;

    /// Return the number of retained low-rank factor columns.
    Eigen::Index modes() const;

    /// Return the background mean in full stamp-pixel order.
    const vectorT &mean() const;

    /// Solve the marginal covariance on an ordered subset of distinct valid stamp pixels.
    /** Uses diagonal scaling, Householder QR, and a Cholesky solve within the factor's span. This is algebraically
     * equivalent to Woodbury without its subtractive cancellation in strongly downweighted directions. The
     * factorization is rebuilt for each support; changing masks never reuses an incompatible full-stamp inverse.
     */
    vectorT
    solve( const vectorT &value, /**< [in] finite vector in the same order as support */
           const std::vector<Eigen::Index> &support /**< [in] distinct full-stamp indices; may be empty */ ) const;
    ///@}
};

} // namespace improc
} // namespace mx

#endif // PSFNoiseModel_hpp
