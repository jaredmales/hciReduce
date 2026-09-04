/** \file P4TemporalPCA.hpp
 * \brief Declares time-domain PCA prediction for Pixel Prediction Post-Processing.
 * \author Jared R. Males
 */

#ifndef P4TemporalPCA_hpp
#define P4TemporalPCA_hpp

#include <cstdint>
#include <vector>

#include <Eigen/Dense>

#include <mx/math/eigenLapack.hpp>

namespace mx
{
namespace improc
{

/// Centering applied to every reference-pixel time series before temporal PCA.
/** \ingroup programming_library */
enum class P4TemporalPCACentering : std::uint8_t
{
    pixelMean, ///< Subtract the full retained-series mean independently from every reference pixel.
    none       ///< Preserve the uncentered reference time series.
};

/// One held-out temporal PCA predictor calculation.
/** \ingroup programming_library */
struct P4TemporalPCAResult
{
    Eigen::Array<double, Eigen::Dynamic, 1> predictions;    ///< One gap-held-out prediction per target image.

    Eigen::Array<double, Eigen::Dynamic, 1> residuals;      ///< Target series minus the corresponding prediction.

    Eigen::Array<std::uint8_t, Eigen::Dynamic, 1> validity; ///< One when the target's restricted fit is supported.

    int numericalRank{ 0 }; ///< Numerical rank of the selected temporal reference modes.
};

/// Reusable time-domain PCA basis and gap-held-out central-pixel predictor.
/** The input matrix stores one reference-pixel time series per row. configure() learns a temporal basis once from
 * those reference series. predict() then fits the target series outside a no-wrap image-index gap independently for
 * every target image, returning one prediction suitable for one appended ordinary-P4 predictor column.
 *
 * \ingroup programming_library
 */
class P4TemporalPCA
{
  public:
    /// Dynamic all-double reference-series matrix type.
    using matrixT = Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>;

    /// Dynamic all-double time-series vector type.
    using vectorT = Eigen::Array<double, Eigen::Dynamic, 1>;

    /// Reusable all-double LAPACK workspace.
    using workspaceT = mx::math::syevrMem<double>;

    /// Learn a requested number of temporal PCA modes from finite reference time series.
    void configure( const matrixT &referenceSeries,   /**< [in] finite reference-pixel by target-image matrix */
                    int requestedModes,               /**< [in] positive temporal basis size */
                    P4TemporalPCACentering centering, /**< [in] reference-series centering policy */
                    double rankTolerance,             /**< [in] finite relative temporal-rank threshold */
                    int gapImages,                    /**< [in] nonnegative no-wrap temporal-gap half-width */
                    workspaceT &workspace /**< [in,out] non-shared LAPACK workspace */ );

    /// Predict every target-series sample from data outside its configured temporal gap.
    void predict( P4TemporalPCAResult &result, /**< [out] target-indexed predictions, residuals, and validity */
                  const vectorT &targetSeries /**< [in] finite central-pixel target time series */ ) const;

    /// Return the configured target-image count.
    Eigen::Index sampleCount() const noexcept;

    /// Return the requested temporal basis size.
    int requestedModes() const noexcept;

  private:
    /// Return whether every value of an Eigen-like input remains finite under fast-math.
    template <typename arrayT>
    static bool allFinite( const arrayT &array /**< [in] finite-candidate array */ );

    matrixT m_modes; ///< Target-image by requested-mode temporal eigenvector matrix.

    std::vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd>> m_gapFits;
    ///< Target-indexed least-squares factorizations of intercept-plus-mode rows outside each gap.

    std::vector<std::vector<Eigen::Index>> m_retainedRows; ///< Target-indexed retained target-series row indices.

    int m_requestedModes{ 0 };                             ///< Positive selected temporal basis size.

    int m_numericalRank{ 0 }; ///< Selected-basis numerical rank after relative thresholding.
};

} // namespace improc
} // namespace mx

#endif // P4TemporalPCA_hpp
