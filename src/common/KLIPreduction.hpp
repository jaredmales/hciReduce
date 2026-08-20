/** \file KLIPreduction.hpp
 * \author Jared R. Males
 * \brief Declarations and definitions for an implementation of the
 * Karhunen-Loeve Image Processing (KLIP) algorithm.
 *
 */

#ifndef __KLIPreduction_hpp__
#define __KLIPreduction_hpp__

#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <omp.h>

#include <mx/ipc/ompLoopWatcher.hpp>
#include <mx/ioutils/stringUtils.hpp>
#include <mx/math/eigenLapack.hpp>
#include <mx/math/floatUtils.hpp>
#include <mx/math/geo.hpp>
#include <mx/sigproc/gramSchmidt.hpp>
using namespace mx::sigproc;

#include "ADIobservation.hpp"
#include "ReductionTiming.hpp"

namespace mx
{
namespace improc
{

// double t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, tf;
// double dcut, dscv, dklims, dgemm, dsyevr, dcfs, drot, dcombo, dread;

/// \cond KLIPreduction_detail
/// Calculate KL modes with the smaller of the reference- and pixel-space Gram matrices.
template <typename _evCalcT = double, typename eigenT, typename eigenT1>
MXLAPACK_INT calcKLModesAdaptive(
    eigenT &klModes,                /**< [out] KL modes stored by row in ascending eigenvalue order */
    eigenT &referenceCovariance,    /**< [in] reference-space Gram matrix; unused when pixels < references */
    const eigenT1 &referenceImages, /**< [in] vectorized reference images, with images stored by column */
    int maximumModeCount = 0,       /**< [in] maximum number of largest-eigenvalue modes, or zero for all */
    math::syevrMem<_evCalcT> *workspace = nullptr, /**< [in,out] optional reusable eigensolver workspace */
    double *eigenSeconds = nullptr,                /**< [out] optional eigensolver elapsed time */
    double *modeSeconds = nullptr /**< [out] optional spatial-mode construction elapsed time */ )
{
    const int pixelCount = referenceImages.rows();
    const int referenceCount = referenceImages.cols();

    if( eigenSeconds != nullptr )
    {
        *eigenSeconds = 0;
    }
    if( modeSeconds != nullptr )
    {
        *modeSeconds = 0;
    }
    if( pixelCount <= 0 || referenceCount <= 0 )
    {
        return -1;
    }

    if( referenceCount <= pixelCount )
    {
        return math::calcKLModes<_evCalcT>( klModes,
                                            referenceCovariance,
                                            referenceImages,
                                            maximumModeCount,
                                            workspace,
                                            eigenSeconds,
                                            modeSeconds );
    }

    const int modeCount = maximumModeCount <= 0 ? pixelCount : std::min( maximumModeCount, pixelCount );
    eigenT pixelCovariance;
    eigenT eigenvectors;
    eigenT eigenvalues;
    pixelCovariance = ( referenceImages.matrix() * referenceImages.matrix().transpose() ).array();

    const MXLAPACK_INT solverStatus = math::calcEigenVecs<_evCalcT>( eigenvectors,
                                                                     eigenvalues,
                                                                     pixelCovariance,
                                                                     modeCount,
                                                                     false,
                                                                     false,
                                                                     workspace,
                                                                     eigenSeconds );
    if( solverStatus != 0 )
    {
        return solverStatus;
    }
    const double modeStart = modeSeconds == nullptr ? 0 : sys::get_curr_time();
    klModes = eigenvectors.matrix().transpose().array();
    for( int mode = 0; mode < modeCount; ++mode )
    {
        bool validMode = math::isFinite( eigenvalues( mode ) ) && eigenvalues( mode ) > 0;
        for( int pixel = 0; validMode && pixel < pixelCount; ++pixel )
        {
            validMode = math::isFinite( klModes( mode, pixel ) );
        }
        if( !validMode )
        {
            klModes.row( mode ).setZero();
        }
    }
    if( modeSeconds != nullptr )
    {
        *modeSeconds = sys::get_curr_time() - modeStart;
    }

    return 0;
}
/// \endcond

/// An implementation of the Karhunen-Loeve Image Processing (KLIP) algorithm.
/** KLIP\cite soummer_2012 is a principle components analysis (PCA) based
 * technique for PSF estimation.
 *
 * For a vectorized reference library \f$X\f$ with \f$P\f$ region pixels and \f$R\f$ reference images, the mode
 * calculation automatically uses the smaller Gram matrix. It retains the reference-space \f$X^T X\f$ path when
 * \f$R \leq P\f$ (including the equal-dimension tie) and uses the pixel-space \f$X X^T\f$ path when \f$P < R\f$.
 * Returned mode rows are the largest-eigenvalue modes in ascending order, and a request larger than the available
 * dimension is clamped to \f$\min(P,R)\f$. Nonpositive or non-finite modes are represented by zero rows.
 *
 * When diagnostics are enabled, `cv.fits` remains the complete symmetric reference-space \f$X^T X\f$ covariance
 * even when the pixel-space solve is selected. At a truncation boundary that crosses an exactly or nearly degenerate
 * eigenvalue cluster, finite-precision eigensolvers may choose different valid bases (and potentially different
 * truncated subspaces); include the complete cluster when that distinction matters.
 *
 * \tparam _realT  is the floating point type in which to do calculations
 * \tparam _derotFunctObj the ADIobservation derotator class.
 * \tparam _evCalcT the real type in which to do eigen-decomposition.  Should
 * generally be double for stable results. \ingroup programming_library
 */
template <typename _realT, class _derotFunctObj, typename _evCalcT, class verboseT>
struct KLIPreduction : public ADIobservation<_realT, _derotFunctObj, verboseT>
{
    typedef _realT realT;
    typedef _evCalcT evCalcT;

    typedef Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> imageT;

    /// FITS-header type shared with ADIobservation final processing.
    typedef typename ADIobservation<realT, _derotFunctObj, verboseT>::fitsHeaderT fitsHeaderT;

    /// Default c'tor
    KLIPreduction();

    virtual ~KLIPreduction();

    /// Register KLIP configuration targets.
    void setupConfig( mx::app::appConfigurator &config /**< [in.out] application configuration */ );

    /// Load KLIP settings from an application configuration.
    void loadConfig( mx::app::appConfigurator &config /**< [in] parsed application configuration */ );

    /// Specify how the data are centered for PCA within each search region
    /** Can have the following values:
     * - <b>HCI::meanSub::none</b> = no value is subtracted
     * - <b>HCI::meanSub::imageMean</b> = the mean of each image (within the search
     * region) is subtracted from itself
     * - <b>HCI::meanSub::imageMedian</b> = the median of each image (within the search
     * region) is subtracted from itself
     * - <b>HCI::meanSub::imageMode</b>  = the mode of each image (within the search
     * region) is subtracted from itself
     * - <b>HCI::meanSub::meanImage</b> = the mean image of the data is subtracted from
     * each image
     * - <b>HCI::meanSub::medianImage</b> = the median image of the data is subtracted
     * from each image
     *
     * \note `HCI::meanSub::imageMode` is recognized by the configuration parser but is not implemented.
     */
    HCI::meanSub m_meanSubMethod{ HCI::meanSub::imageMean };

    /// Specify if each pixel time-series is normalized
    /** This normalizaton is applied after centering. Can have the following values:
     * - <b>HCI::pixelTSNorm::none</b>: no normalization (the default)
     * - <b>HCI::pixelTSNorm::rms</b>: divide by the time-series rms
     * - <b>HCI::pixelTSNorm::rmsSigmaClipped</b>: divide by the sigma-slipped time-series rms.
     *                                                   The sigma is provided by m_pixelTSSigma.
     */
    HCI::pixelTSNorm m_pixelTSNormMethod{ HCI::pixelTSNorm::none };

    realT m_pixelTSSigma{ 3 }; ///< Sigma-clipping parameter for pixel time-series normalization

    int m_padSize{ 4 };

    /// Specifies the number of modes to include in the PSF.
    /** The output image is a cube with a plane for each entry in m_Nmodes.
     * Only the number of eigenvalues required for the maximum value of m_Nmodes
     * are calculated, so this can have an important impact on speed.
     *
     * Can be initialized as:
     * \code
     * red.m_Nmodes={5,10,15,20};
     * \endcode
     *
     */
    std::vector<int> m_Nmodes;

    std::vector<realT> m_minRadius;
    std::vector<realT> m_maxRadius;
    std::vector<realT> m_minAngle;
    std::vector<realT> m_maxAngle;
    int m_nWedges{ 0 };

    int m_maxNmodes{ 0 };

    /// Specify the minimum pixel difference at the inner edge of the search
    /// region
    realT m_minDPx{ 0 };

    /// Specify the maximum pixel difference at the inner edge of the search
    /// region
    realT m_maxDPx{ 0 };

    /// Controls how reference images are excluded, if at all, from the
    /// covariance matrix for each target image based on a minimum criterion.
    /** Can have the following values:
     *  - <b>HCI::exclude::none</b> = no exclusion, all images included [default]
     *  - <b>HCI::exclude::pixel</b> = exclude based on pixels of rotation at the
     * inner edge of the region
     *  - <b>HCI::exclude::angle</b> = exclude based on degrees of rotation at the
     * inner edge of the region
     *  - <b>HCI::exclude::imno</b> = exclude based on number of images
     */
    HCI::exclude m_excludeMethod{ HCI::exclude::none };

    /// Controls how reference images are excluded, if at all, from the
    /// covariance matrix for each target image based on a maximum criterion.
    /** Can have the following values:
     *  - <b>HCI::exclude::none</b> = no exclusion, all images included [default]
     *  - <b>HCI::exclude::pixel</b> = exclude based on pixels of rotation at the
     * inner edge of the region
     *  - <b>HCI::exclude::angle</b> = exclude based on degrees of rotation at the
     * inner edge of the region
     *  - <b>HCI::exclude::imno</b> = exclude based on number of images
     */
    HCI::exclude m_excludeMethodMax{ HCI::exclude::none };

    /// Maximum number of reference images to include in the covariance matrix.
    /** If greater than zero, at most this many images are retained according to
     * m_includeMethod after exclusion. If zero, all surviving reference images
     * are included.
     */
    int m_includeRefNum{ 0 };

    /// Controls how number of included images is calculated.
    /** The number of included images is calculated after exclusion is complete.
     * Can have the following values:
     * - <b>HCI::include::all</b> = all remaining images are included, ignoring m_includeRefNum [default]
     * - <b>HCI::include::corr</b> = the m_includeRefNum of the remaining images
     * which are most correlated with the target are included
     * - <b>HCI::include::time</b> = the m_includeRefNum of the remaining images
     * which are closest in time to the target are included
     * - <b>HCI::include::angle</b> = the m_includeRefNum of the remaining images
     * which are closest in angle to the target are included
     * - <b>HCI::include::imno</b> = the m_includeRefNum of the remaining images
     * which are closest in image number to the target are included
     */
    HCI::include m_includeMethod{ HCI::include::all };

    eigenImage<int> m_imsIncluded;

    bool m_rightReason = false;

    realT m_rightReasonRadius = 2.5;

    /// Whether intermediate KLIP diagnostic FITS products are written.
    bool m_writeDiagnostics{ false };

    /// Directory in which intermediate KLIP diagnostic FITS products are written.
    std::string m_diagnosticDirectory{ "." };

    /// Write an intermediate KLIP diagnostic FITS product when diagnostics are enabled.
    template <typename dataT>
    void writeDiagnostic( const std::string &fileName, /**< [in] basename of the diagnostic product */
                          const dataT &data /**< [in] image-like data to write */ ) const;

    /// Subtract the basis mean from each of the images
    /** The mean is subtracted according to m_meanSubMethod.
     */
    void meanSubtract( eigenCube<realT> &rims, ///< [in.out] The reference images.  These are
                                               ///< mean subtracted on output.
                       eigenCube<realT> &tims, ///< [in.out] The target images, which can be the same
                                               ///< cube as rims (tested by pointer comparison), in which
                                               ///< case they will be ignored.  Mean subtractedon output.
                       imageT &cmask,          ///< [in] the cutout mask.  Ignored if empty.
                       std::vector<realT> &sds ///< [out] The standard deviation of the mean
                                               ///< subtracted reference images.
    );

    std::vector<realT> m_minr;
    std::vector<realT> m_maxr;
    std::vector<realT> m_minq;
    std::vector<realT> m_maxq;

    /// Run KLIP in a set of geometric search regions.
    /** The arguments are 4 vectors, where each entry defines one component of
     * the search region.
     *
     * \returns 0 on success
     * \returns -1 on error
     *
     */
    int regions( const std::vector<realT> &minr, ///< [in]
                 const std::vector<realT> &maxr, ///< [in]
                 const std::vector<realT> &minq, ///< [in]
                 const std::vector<realT> &maxq  ///< [in]
    );

    /// Run KLIP in a geometric search region.
    /** \overload
     *
     * \returns 0 on success
     * \returns -1 on error
     *
     */
    int regions( realT minr, ///< [in]
                 realT maxr, ///< [in]
                 realT minq, ///< [in]
                 realT maxq  ///< [in]
    )
    {
        std::vector<realT> vminr( 1, minr );
        std::vector<realT> vmaxr( 1, maxr );
        std::vector<realT> vminq( 1, minq );
        std::vector<realT> vmaxq( 1, maxq );

        return regions( vminr, vmaxr, vminq, vmaxq );
    }

    /// Perform KLIP subtraction for one flattened image region.
    void worker( eigenCube<realT> &rims,   /**< [in.out] flattened reference images */
                 eigenCube<realT> &tims,   /**< [in.out] flattened target images */
                 imageT &cmask,            /**< [in] flattened optional mask */
                 std::vector<size_t> &idx, /**< [in] full-image indices for the region */
                 realT dang,               /**< [in] minimum exclusion threshold */
                 realT dangMax /**< [in] maximum exclusion threshold */ );

    /// Append KLIP-specific reduction parameters to a FITS header.
    void appendReductionHeader( fitsHeaderT &head /**< [in.out] header receiving the KLIP cards */ );

    /// Apply the shared ADI final-processing lifecycle with KLIP provenance.
    /** \returns 0 on success.
     * \throws mx::exception if final processing or output fails.
     */
    int finalProcess();

    /// Load a saved KLIP reduction and apply the currently configured final-processing stages.
    /** The saved mode list is restored from the required `NMODES` header card; post-median subtraction, derotation,
     * combination, and output behavior remain controlled by the current configuration.
     *
     * \returns 0 on success.
     * \throws mx::exception if the saved collection or mode metadata is invalid.
     */
    int processPSFSub( const std::string &directory, /**< [in] directory containing the saved products */
                       const std::string &prefix,    /**< [in] literal filename prefix */
                       const std::string &extension = ".fits" /**< [in] literal filename extension */ );

  private:
    /// Instance-owned KLIP algorithm timing record for the current reduction.
    ReductionTiming m_algorithmTiming;

  public:
    /// Return the current instance-owned KLIP algorithm timing snapshot.
    const ReductionTiming &algorithmTiming() const
    {
        return m_algorithmTiming;
    }

    /// Print the current reduction-stage timing summary.
    void dump_times()
    {
        printf( "KLIP reduction times: \n" );
        printf( "  Total time: %f sec\n", this->t_end - this->t_begin );
        printf( "    Loading: %f sec\n", this->t_load_end - this->t_load_begin );
        printf( "    Fake Injection: %f sec\n", this->t_fake_end - this->t_fake_begin );
        printf( "    Coadding: %f sec\n", this->t_coadd_end - this->t_coadd_begin );
        printf( "    Preprocessing: %f sec\n", this->t_preproc_end - this->t_preproc_begin );
        printf( "      Az USM: %f sec\n", this->t_azusm_end - this->t_azusm_begin );
        printf( "      Gauss USM: %f sec\n", this->t_gaussusm_end - this->t_gaussusm_begin );
        printf( "    KLIP algorithm: %f elapsed real sec\n", m_algorithmTiming.regressionElapsedSeconds );
        const double klipWorker = m_algorithmTiming.eigensolveWorkerSeconds + m_algorithmTiming.modeWorkerSeconds +
                                  m_algorithmTiming.projectionWorkerSeconds;
        const auto percentage = [klipWorker]( double seconds )
        { return klipWorker > 0 ? seconds / klipWorker * 100 : 0; };
        printf( "      EigenDecomposition %f cpu sec (%f%%)\n",
                m_algorithmTiming.eigensolveWorkerSeconds,
                percentage( m_algorithmTiming.eigensolveWorkerSeconds ) );
        printf( "      KL image calc %f cpu sec (%f%%)\n",
                m_algorithmTiming.modeWorkerSeconds,
                percentage( m_algorithmTiming.modeWorkerSeconds ) );
        printf( "      PSF calc/sub %f cpu sec (%f%%)\n",
                m_algorithmTiming.projectionWorkerSeconds,
                percentage( m_algorithmTiming.projectionWorkerSeconds ) );
        printf( "    Derotation: %f sec\n", this->t_derotate_end - this->t_derotate_begin );
        printf( "    Combination: %f sec\n", this->t_combo_end - this->t_combo_begin );
    }
};

template <typename realT, class derotFunctObj, typename evCalcT, class verboseT>
KLIPreduction<realT, derotFunctObj, evCalcT, verboseT>::KLIPreduction()
{
}

template <typename realT, class derotFunctObj, typename evCalcT, class verboseT>
KLIPreduction<realT, derotFunctObj, evCalcT, verboseT>::~KLIPreduction()
{
}

template <typename realT, class derotFunctObj, typename evCalcT, class verboseT>
void KLIPreduction<realT, derotFunctObj, evCalcT, verboseT>::setupConfig( mx::app::appConfigurator &config )
{
    HCIobservation<realT, verboseT>::setupConfig( config );

    config.add( "adi.minDPx",
                "",
                "adi.minDPx",
                mx::app::argType::Required,
                "adi",
                "minDPx",
                false,
                "float",
                "Specify the minimum angle or pixel difference at the inner edge of the search region" );

    config.add( "adi.maxDPx",
                "",
                "adi.maxDPx",
                mx::app::argType::Required,
                "adi",
                "maxDPx",
                false,
                "float",
                "Specify the maximum angle or pixel difference at the inner edge of the search region" );

    config.add( "adi.excludeMethod",
                "",
                "adi.excludeMethod",
                mx::app::argType::Required,
                "adi",
                "excludeMethod",
                false,
                "string",
                "Method for minimum exclusion.  Values are none (default), pixel, angle, imno." );

    config.add( "adi.excludeMethodMax",
                "",
                "adi.excludeMethodMax",
                mx::app::argType::Required,
                "adi",
                "excludeMethodMax",
                false,
                "string",
                "Method for maximum exclusion.  Values are none (default), pixel, angle, imno." );

    ADIobservation<realT, derotFunctObj, verboseT>::setupConfig( config );

    /*>>>> geom */
    config.add( "geom.minRadius",
                "",
                "geom.minRadius",
                mx::app::argType::Required,
                "geom",
                "minRadius",
                false,
                "vector<realT>",
                "The minimum radius of the search regions" );

    config.add( "geom.maxRadius",
                "",
                "geom.maxRadius",
                mx::app::argType::Required,
                "geom",
                "maxRadius",
                false,
                "vector<realT>",
                "The maximum radius of the search regions" );

    config.add( "geom.minAngle",
                "",
                "geom.minAngle",
                mx::app::argType::Required,
                "geom",
                "minAngle",
                false,
                "vector<realT>",
                "The minimum angle of the search regions" );

    config.add( "geom.maxAngle",
                "",
                "geom.maxAngle",
                mx::app::argType::Required,
                "geom",
                "maxAngle",
                false,
                "vector<realT>",
                "The maximum angle of the search regions" );

    config.add( "geom.nWedges",
                "",
                "geom.nWedges",
                mx::app::argType::Required,
                "geom",
                "nWedges",
                false,
                "",
                "The number of angular wedges.  Overrides minAngle and maxAngle, "
                "and expands minRadius and maxRadius" );

    /*>>>> klip */
    config.add( "klip.meanSubMethod",
                "",
                "klip.meanSubMethod",
                mx::app::argType::Required,
                "klip",
                "meanSubMethod",
                false,
                "string",
                "The method of mean subtraction for PCA: imageMean, imageMedian, meanImage, or medianImage." );

    config.add( "klip.pixelTSNormMethod",
                "",
                "klip.pixelTSNormMethod",
                mx::app::argType::Required,
                "klip",
                "pixelTSNormMethod",
                false,
                "string",
                "The method of pixel time-series normalization for PCA: none or rms." );

    config.add( "klip.pixelTSSigma",
                "",
                "klip.pixelTSSigma",
                mx::app::argType::Required,
                "klip",
                "pixelTSSigma",
                false,
                "float",
                "Sigma-clipping threshold for pixel time-series normalization" );

    config.add( "klip.includeRefNum",
                "",
                "klip.includeRefNum",
                mx::app::argType::Required,
                "klip",
                "includeRefNum",
                false,
                "int",
                "The maximum number of references to include after exclusion." );

    config.add( "klip.includeMethod",
                "",
                "klip.includeMethod",
                mx::app::argType::Required,
                "klip",
                "includeMethod",
                false,
                "string",
                "Reference inclusion method: all, corr, time, angle, or imno." );

    config.add( "klip.Nmodes",
                "",
                "klip.Nmodes",
                mx::app::argType::Required,
                "klip",
                "Nmodes",
                false,
                "vector<int>",
                "The number of modes to included in the PSF estimate." );

    config.add( "klip.rightReason",
                "",
                "klip.rightReason",
                mx::app::argType::Required,
                "klip",
                "rightReason",
                false,
                "bool",
                "Whether or not the right reason mask is applied" );

    config.add( "klip.rrRadius",
                "",
                "klip.rrRadius",
                mx::app::argType::Required,
                "klip",
                "rrRadius",
                false,
                "float",
                "The radius of the right reason mask" );

    config.add( "klip.writeDiagnostics",
                "",
                "klip.writeDiagnostics",
                mx::app::argType::True,
                "klip",
                "writeDiagnostics",
                false,
                "bool",
                "Whether intermediate KLIP diagnostic FITS products are written" );

    config.add( "klip.diagnosticDirectory",
                "",
                "klip.diagnosticDirectory",
                mx::app::argType::Required,
                "klip",
                "diagnosticDirectory",
                false,
                "string",
                "Directory in which intermediate KLIP diagnostic FITS products are written" );

    /*<<<< klip */
}

template <typename realT, class derotFunctObj, typename evCalcT, class verboseT>
void KLIPreduction<realT, derotFunctObj, evCalcT, verboseT>::loadConfig( mx::app::appConfigurator &config )
{
    HCIobservation<realT, verboseT>::loadConfig( config );
    config( m_minDPx, "adi.minDPx" );
    config( m_maxDPx, "adi.maxDPx" );

    std::string em = HCI::excludeToStr<verboseT>( m_excludeMethod );
    config( em, "adi.excludeMethod" );
    m_excludeMethod = HCI::excludeFmStr<verboseT>( em );

    em = HCI::excludeToStr<verboseT>( m_excludeMethodMax );
    config( em, "adi.excludeMethodMax" );
    m_excludeMethodMax = HCI::excludeFmStr<verboseT>( em );

    ADIobservation<realT, derotFunctObj, verboseT>::loadConfig( config );

    /*>>>> geom */
    config( m_minRadius, "geom.minRadius" );
    config( m_maxRadius, "geom.maxRadius" );
    config( m_minAngle, "geom.minAngle" );
    config( m_maxAngle, "geom.maxAngle" );
    config( m_nWedges, "geom.nWedges" );
    /*<<<< geom */

    /*>>>> klip */
    std::string ms;
    try
    {
        ms = HCI::meanSubToStr<verboseT>( m_meanSubMethod );
    }
    catch( ... )
    {
        std::throw_with_nested( mx::exception<verboseT>( mx::error_t::exception ) );
    }

    config( ms, "klip.meanSubMethod" );

    try
    {
        m_meanSubMethod = HCI::meanSubFmStr<verboseT>( ms );
    }
    catch( ... )
    {
        std::throw_with_nested(
            mx::exception<verboseT>( mx::error_t::invalidconfig, "klip.meanSubMethod is not valid" ) );
    }

    std::string ptsnm;
    try
    {
        ptsnm = HCI::pixelTSNormToStr<verboseT>( m_pixelTSNormMethod );
    }
    catch( ... )
    {
        std::throw_with_nested( mx::exception<verboseT>( mx::error_t::exception ) );
    }

    config( ptsnm, "klip.pixelTSNormMethod" );

    try
    {
        m_pixelTSNormMethod = HCI::pixelTSNormFmStr<verboseT>( ptsnm );
    }
    catch( ... )
    {
        std::throw_with_nested( mx::exception<verboseT>( mx::error_t::invalidconfig,
                                                         "invalid pixel time-series "
                                                         "normalization method" ) );
    }

    config( m_includeRefNum, "klip.includeRefNum" );
    if( m_includeRefNum < 0 )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig, "klip.includeRefNum must be non-negative" );
    }

    std::string includeMethod;
    try
    {
        includeMethod = HCI::includeToStr<verboseT>( m_includeMethod );
    }
    catch( ... )
    {
        std::throw_with_nested( mx::exception<verboseT>( mx::error_t::exception ) );
    }
    config( includeMethod, "klip.includeMethod" );
    try
    {
        m_includeMethod = HCI::includeFmStr<verboseT>( includeMethod );
    }
    catch( ... )
    {
        std::throw_with_nested(
            mx::exception<verboseT>( mx::error_t::invalidconfig, "klip.includeMethod is not valid" ) );
    }
    config( m_Nmodes, "klip.Nmodes" );
    config( m_rightReason, "klip.rightReason" );
    config( m_rightReasonRadius, "klip.rrRadius" );
    config( m_pixelTSSigma, "klip.pixelTSSigma" );
    config( m_writeDiagnostics, "klip.writeDiagnostics" );
    config( m_diagnosticDirectory, "klip.diagnosticDirectory" );
}

template <typename realT, class derotFunctObj, typename evCalcT, class verboseT>
template <typename dataT>
void KLIPreduction<realT, derotFunctObj, evCalcT, verboseT>::writeDiagnostic( const std::string &fileName,
                                                                              const dataT &data ) const
{
    if( !m_writeDiagnostics )
    {
        return;
    }

    std::string path = fileName;
    if( !m_diagnosticDirectory.empty() && m_diagnosticDirectory != "." )
    {
        path = m_diagnosticDirectory + "/" + fileName;
    }

    const std::string outputParent = ioutils::parentPath( path );
    if( !outputParent.empty() )
    {
        const mx::error_t result = ioutils::createDirectories( outputParent );
        if( result != mx::error_t::noerror )
        {
            throw mx::exception<verboseT>( result, "could not create KLIP diagnostic output directory" );
        }
    }

    fits::fitsFile<typename dataT::Scalar, verboseT> writer;
    const mx::error_t result = writer.write( path, data );
    if( result != mx::error_t::noerror )
    {
        throw mx::exception<verboseT>( result, "could not write KLIP diagnostic product " + path );
    }
}

template <typename realT, class derotFunctObj, typename evCalcT, class verboseT>
void KLIPreduction<realT, derotFunctObj, evCalcT, verboseT>::meanSubtract( eigenCube<realT> &rims,
                                                                           eigenCube<realT> &tims,
                                                                           imageT &cmask,
                                                                           std::vector<realT> &norms )
{
    if( m_meanSubMethod == HCI::meanSub::imageMode )
    {
        throw mx::exception<verboseT>( mx::error_t::notimpl, "image-mode subtraction is not implemented" );
    }
    if( m_meanSubMethod != HCI::meanSub::none && m_meanSubMethod != HCI::meanSub::meanImage &&
        m_meanSubMethod != HCI::meanSub::medianImage && m_meanSubMethod != HCI::meanSub::imageMean &&
        m_meanSubMethod != HCI::meanSub::imageMedian )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig, "invalid KLIP mean-subtraction method" );
    }

    if( m_pixelTSNormMethod == HCI::pixelTSNorm::rmsSigmaClipped )
    {
        throw mx::exception<verboseT>( mx::error_t::notimpl,
                                       "pixelTSNormMethod is rmsSigmaClipped, which is not implemented" );
    }
    if( m_pixelTSNormMethod != HCI::pixelTSNorm::none && m_pixelTSNormMethod != HCI::pixelTSNorm::rms )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                       "invalid KLIP pixel time-series normalization method" );
    }

    if( rims.rows() <= 0 || rims.cols() <= 0 || rims.planes() <= 0 || tims.rows() != rims.rows() ||
        tims.cols() != rims.cols() || tims.planes() <= 0 )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidarg, "invalid reference or target cube geometry" );
    }

    int maskPix = 0;

    if( ( cmask.rows() > 0 ) != ( cmask.cols() > 0 ) )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidarg, "KLIP mask has incomplete geometry" );
    }

    if( cmask.rows() > 0 && cmask.cols() > 0 )
    {
        if( cmask.rows() != rims.rows() || cmask.cols() != rims.cols() )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidarg,
                                           "mask does not have same size as reference images" );
        }

        maskPix = cmask.sum();
        if( maskPix <= 0 )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidarg, "KLIP mask does not select any pixels" );
        }
    }

    norms.resize( rims.planes() );

    if( m_meanSubMethod == HCI::meanSub::meanImage || m_meanSubMethod == HCI::meanSub::medianImage )
    {
        imageT mean;

        if( m_meanSubMethod == HCI::meanSub::meanImage )
        {
            rims.mean( mean );
        }
        else if( m_meanSubMethod == HCI::meanSub::medianImage )
        {
            rims.median( mean );
        }

        realT immean;
        for( int n = 0; n < rims.planes(); ++n )
        {
            rims.image( n ) -= mean;

            if( maskPix > 0 )
            {
                rims.image( n ) *= cmask;

                immean = rims.image( n ).sum() / maskPix;

                norms[n] = std::sqrt( ( ( rims.image( n ) - immean ) * cmask ).square().sum() );
            }
            else
            {
                immean = rims.image( n ).mean();

                norms[n] = std::sqrt( ( rims.image( n ) - immean ).square().sum() );
            }
        }

        if( &tims != &rims )
        {
            for( int n = 0; n < tims.planes(); ++n )
            {
                tims.image( n ) -= mean;

                if( maskPix > 0 )
                {
                    tims.image( n ) *= cmask;
                }
            }
        }
    }
    else
    {
        realT mean{ 0 };
        std::vector<realT> work; // Working memmory for median calc

        realT immean;
        for( int n = 0; n < rims.planes(); ++n )
        {
            if( m_meanSubMethod == HCI::meanSub::imageMean )
            {
                mean = maskPix > 0 ? ( rims.image( n ) * cmask ).sum() / maskPix : rims.image( n ).mean();
            }
            else if( m_meanSubMethod == HCI::meanSub::imageMedian )
            {
                mean =
                    maskPix > 0 ? imageMedian( rims.image( n ), &cmask, &work ) : imageMedian( rims.image( n ), &work );
            }

            rims.image( n ) -= mean;

            if( maskPix > 0 )
            {
                rims.image( n ) *= cmask;
                immean = rims.image( n ).sum() / maskPix;
                norms[n] = std::sqrt( ( ( rims.image( n ) - immean ) * cmask ).square().sum() );
            }
            else
            {
                // Because we might not have used the mean, we need to re-mean to
                //  make this the standard deviation
                immean = rims.image( n ).mean();
                norms[n] = std::sqrt( ( rims.image( n ) - immean ).square().sum() );
            }
        }

        if( &tims != &rims )
        {
            for( int n = 0; n < tims.planes(); ++n )
            {
                if( maskPix > 0 )
                {
                    if( m_meanSubMethod == HCI::meanSub::imageMean )
                    {
                        mean = ( tims.image( n ) * cmask ).sum() / maskPix;
                    }
                    else if( m_meanSubMethod == HCI::meanSub::imageMedian )
                    {
                        mean = imageMedian( tims.image( n ), &cmask, &work );
                    }

                    tims.image( n ) -= mean;
                    tims.image( n ) *= cmask;
                }
                else
                {
                    if( m_meanSubMethod == HCI::meanSub::imageMean )
                    {
                        mean = tims.image( n ).mean();
                    }
                    else if( m_meanSubMethod == HCI::meanSub::imageMedian )
                    {
                        mean = imageMedian( tims.image( n ), &work );
                    }

                    tims.image( n ) -= mean;
                }
            }
        }
    }

    if( m_pixelTSNormMethod == HCI::pixelTSNorm::rms )
    {
        std::cerr << "normalizing pixels\n";

        for( int cc = 0; cc < rims.cols(); ++cc )
        {
            for( int rr = 0; rr < rims.rows(); ++rr )
            {
                if( maskPix > 0 )
                {
                    if( cmask( rr, cc ) == 0 )
                    {
                        continue;
                    }
                }

                realT sumSquares{ 0 };
                for( int pp = 0; pp < rims.planes(); ++pp )
                {
                    sumSquares += rims.image( pp )( rr, cc ) * rims.image( pp )( rr, cc );
                }

                const realT rms = std::sqrt( sumSquares / rims.planes() );

                if( std::isfinite( rms ) && rms > std::numeric_limits<realT>::epsilon() )
                {
                    for( int pp = 0; pp < rims.planes(); ++pp )
                    {
                        rims.image( pp )( rr, cc ) /= rms;
                    }
                    if( &tims != &rims )
                    {
                        for( int pp = 0; pp < tims.planes(); ++pp )
                        {
                            tims.image( pp )( rr, cc ) /= rms;
                        }
                    }
                }
                else
                {
                    for( int pp = 0; pp < rims.planes(); ++pp )
                    {
                        rims.image( pp )( rr, cc ) = 0;
                    }
                    if( &tims != &rims )
                    {
                        for( int pp = 0; pp < tims.planes(); ++pp )
                        {
                            tims.image( pp )( rr, cc ) = 0;
                        }
                    }
                }
            }
        }
    }
}

template <typename realT, class derotFunctObj, typename evCalcT, class verboseT>
int KLIPreduction<realT, derotFunctObj, evCalcT, verboseT>::regions( const std::vector<realT> &minr,
                                                                     const std::vector<realT> &maxr,
                                                                     const std::vector<realT> &minq,
                                                                     const std::vector<realT> &maxq )
{
    m_algorithmTiming.reset();

    const bool preprocessingOnly = this->preprocessingOnly();

    this->t_begin = sys::get_curr_time();

    if( !preprocessingOnly )
    {
        if( m_Nmodes.empty() )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidarg, "KLIP requires at least one mode count" );
        }
        if( minr.empty() || minr.size() != maxr.size() || minr.size() != minq.size() || minr.size() != maxq.size() )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidarg,
                                           "KLIP region vectors must be nonempty and have equal lengths" );
        }
        for( const int modeCount : m_Nmodes )
        {
            if( modeCount <= 0 )
            {
                throw mx::exception<verboseT>( mx::error_t::invalidarg, "KLIP mode counts must be positive" );
            }
        }
        for( size_t region = 0; region < minr.size(); ++region )
        {
            if( !math::isFinite( minr[region] ) || !math::isFinite( maxr[region] ) || !math::isFinite( minq[region] ) ||
                !math::isFinite( maxq[region] ) || minr[region] < 0 || maxr[region] <= minr[region] )
            {
                throw mx::exception<verboseT>( mx::error_t::invalidarg, "invalid KLIP region geometry" );
            }
        }
        if( m_padSize < 0 )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidarg, "KLIP image padding must be nonnegative" );
        }

        m_minr = minr;
        m_maxr = maxr;
        m_minq = minq;
        m_maxq = maxq;

        m_maxNmodes = m_Nmodes[0];
        for( size_t i = 1; i < m_Nmodes.size(); ++i )
        {
            if( m_Nmodes[i] > m_maxNmodes )
                m_maxNmodes = m_Nmodes[i];
        }
    }

    std::cerr << "Beginning\n";

    if( !preprocessingOnly && this->m_imSize == 0 )
    {
        this->m_imSize = 2 * ( *std::max_element( m_maxr.begin(), m_maxr.end() ) + m_padSize );

        std::cerr << "set image size based on regions to " << this->m_imSize << "\n";
    }

    if( !this->m_filesRead )
    {
        try
        {
            this->readFiles();
        }
        catch( const std::exception &e )
        {
            std::throw_with_nested( mx::exception<verboseT>( mx::error_t::std_exception, "from readFiles" ) );
        }
    }

    // CHECK IF RDI HERE
    if( !this->m_RDIfilesRead && this->m_RDIfileList.size() != 0 )
    {
        try
        {
            this->readRDIFiles();
        }
        catch( ... )
        {
            std::throw_with_nested( mx::exception<verboseT>( mx::error_t::exception, "from readRDIFiles" ) );
        }
    }

    if( this->m_Nims <= 0 || this->m_Nrows <= 0 || this->m_Ncols <= 0 || this->m_tgtIms.planes() != this->m_Nims ||
        this->m_tgtIms.rows() != this->m_Nrows || this->m_tgtIms.cols() != this->m_Ncols )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "invalid target cube state for KLIP regions" );
    }
    if( this->m_refIms.planes() > 0 &&
        ( this->m_refIms.rows() != this->m_Nrows || this->m_refIms.cols() != this->m_Ncols ) )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "reference and target image dimensions must match" );
    }

    if( preprocessingOnly )
    {
        std::cerr << "Pre-processing complete, stopping.\n";
        return 0;
    }

    std::cerr << "allocating psf subtracted cubes\n";
    this->m_psfsub.resize( m_Nmodes.size() );
    for( size_t n = 0; n < m_Nmodes.size(); ++n )
    {
        this->m_psfsub[n].resize( this->m_Nrows, this->m_Ncols, this->m_Nims );
        this->m_psfsub[n].cube().setZero();
    }

    // Make radius and angle images
    imageT rIm( this->m_Nrows, this->m_Ncols );
    imageT qIm( this->m_Nrows, this->m_Ncols );

    radAngImage<math::degreesT<realT>>( rIm, qIm, .5 * ( this->m_Nrows - 1 ), .5 * ( this->m_Ncols - 1 ) );

    const int referenceCount = this->m_refIms.planes() > 0 ? this->m_refIms.planes() : this->m_Nims;
    m_imsIncluded.resize( this->m_Nims, referenceCount );
    m_imsIncluded.setConstant( 1 );

    if( this->m_refIms.planes() > 0 ) // RDI
    {
        std::cerr << "******* RDI MODE **********\n";
    }
    else // ADI
    {
        std::cerr << "******* ADI MODE **********\n";
    }

    std::cerr << "processing " << minr.size() << " regions\n";

    //******** For each region do this:
    const double workerBegin = sys::get_curr_time();
    for( size_t regno = 0; regno < minr.size(); ++regno )
    {
        std::cerr << "  region " << regno + 1 << ": " << m_minr[regno] << "-" << m_maxr[regno] << " pixels, ";
        std::cerr << m_minq[regno] << "-" << m_maxq[regno] << " degrees.          \n";

        imageT *maskPtr = nullptr;

        if( this->m_mask.rows() == this->m_Nrows && this->m_mask.cols() == this->m_Ncols )
        {
            maskPtr = &this->m_mask;
        }

        std::vector<size_t> idx = annulusIndices<math::degreesT<realT>>( rIm,
                                                                         qIm,
                                                                         .5 * ( this->m_Nrows - 1 ),
                                                                         .5 * ( this->m_Ncols - 1 ),
                                                                         m_minr[regno],
                                                                         m_maxr[regno],
                                                                         m_minq[regno],
                                                                         m_maxq[regno],
                                                                         maskPtr );

        if( idx.empty() )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidarg, "KLIP region contains no usable pixels" );
        }

        // Create storage for the R-ims and psf-subbed Ims
        eigenCube<realT> tims( idx.size(), 1, this->m_Nims );

        //------If doing RDI, create bims
        eigenCube<realT> rims;

        //------Get the mask cutout too
        imageT cmask;

        if( this->m_refIms.planes() > 0 )
        {
            rims.resize( idx.size(), 1, this->m_refIms.planes() );
        }

        for( int i = 0; i < this->m_Nims; ++i )
        {
            auto tim = tims.image( i );
            cutImageRegion( tim, this->m_tgtIms.image( i ), idx, false );
        }

        for( int p = 0; p < this->m_refIms.planes(); ++p )
        {
            auto rim = rims.image( p );
            cutImageRegion( rim, this->m_refIms.image( p ), idx, false );
        }

        if( this->m_maskFile != "" )
        {
            cutImageRegion( cmask, this->m_mask, idx, true );
        }

#if 0

            std::vector<realT> sds;

            //*** First mean subtract ***//
            try{
            meanSubtract( tims, tims, cmask, sds );}
            catch(...)
            {
                std::throw_with_nested(mx::exception<verboseT>(mx::error_t::exception, "from meanSubtract"));
            }


            mx::fits::fitsFile<realT> ffff;
            ffff.write("tims.fits", tims);

            std::ofstream fout("idx.dat");
            std::ofstream aout("derot.dat");
            for(auto i : idx)
            {
                fout << i << "\n";
            }
            fout.close();

            for(int pp =0; pp < tims.planes(); ++pp)
            {
                aout << this->m_derotF.derotAngle( pp ) << "\n";
            }
            aout.close();

            fout.open("dims.dat");
            fout << this->m_Nrows << "\n";
            fout << this->m_Ncols << "\n";
            fout.close();


            exit(0);
#endif

        realT dang = 0;
        realT dangMax = 0;

        const HCI::exclude configuredExcludeMethod = m_excludeMethod;
        const HCI::exclude configuredExcludeMethodMax = m_excludeMethodMax;
        HCI::exclude effectiveExcludeMethod = configuredExcludeMethod;
        HCI::exclude effectiveExcludeMethodMax = configuredExcludeMethodMax;

        if( m_minDPx < 0 || this->m_refIms.planes() > 0 )
        {
            effectiveExcludeMethod = HCI::exclude::none;
        }
        if( m_maxDPx < 0 || this->m_refIms.planes() > 0 )
        {
            effectiveExcludeMethodMax = HCI::exclude::none;
        }

        if( ( effectiveExcludeMethod == HCI::exclude::pixel || effectiveExcludeMethodMax == HCI::exclude::pixel ) &&
            minr[regno] <= 0 )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidarg,
                                           "pixel-based KLIP exclusion requires a positive inner radius" );
        }

        if( effectiveExcludeMethod == HCI::exclude::pixel )
        {
            dang = fabs( atan( m_minDPx / minr[regno] ) );
        }
        else if( effectiveExcludeMethod == HCI::exclude::angle )
        {
            dang = math::dtor( m_minDPx );
        }
        else if( effectiveExcludeMethod == HCI::exclude::imno )
        {
            dang = m_minDPx;
        }

        if( effectiveExcludeMethodMax == HCI::exclude::pixel )
        {
            dangMax = fabs( atan( m_maxDPx / minr[regno] ) );
        }
        else if( effectiveExcludeMethodMax == HCI::exclude::angle )
        {
            dangMax = math::dtor( m_maxDPx );
        }
        else if( effectiveExcludeMethodMax == HCI::exclude::imno )
        {
            dangMax = m_maxDPx;
        }

        m_excludeMethod = effectiveExcludeMethod;
        m_excludeMethodMax = effectiveExcludeMethodMax;
        try
        {
            if( this->m_refIms.planes() > 0 )
            {
                worker( rims, tims, cmask, idx, dang, dangMax );
            }
            else
            {
                worker( tims, tims, cmask, idx, dang, dangMax );
            }
        }
        catch( ... )
        {
            m_excludeMethod = configuredExcludeMethod;
            m_excludeMethodMax = configuredExcludeMethodMax;
            throw;
        }
        m_excludeMethod = configuredExcludeMethod;
        m_excludeMethodMax = configuredExcludeMethodMax;
    }
    m_algorithmTiming.regressionElapsedSeconds = sys::get_curr_time() - workerBegin;

    writeDiagnostic( "imsIncluded.fits", m_imsIncluded );

    finalProcess();

    this->t_end = sys::get_curr_time();

    dump_times();

    return 0;
}

/// Reference-image candidate state used while selecting a target's KLIP library.
struct cvEntry
{
    /// Correlation score or distance used to rank the candidate.
    double cvVal{ 0 };

    /// Whether the candidate remains admissible after exclusion and ranking.
    bool included{ true };
};

/// Extract a square matrix containing the selected rows and columns.
template <typename eigenT, typename eigenTin>
void extractRowsAndCols( eigenT &out,        /**< [out] selected square matrix */
                         const eigenTin &in, /**< [in] source matrix */
                         const std::vector<size_t> &idx /**< [in] ordered source indices */ )
{
    out.resize( idx.size(), idx.size() );

    for( size_t i = 0; i < idx.size(); ++i )
    {
        for( size_t j = 0; j < idx.size(); ++j )
        {
            out( i, j ) = in( idx[i], idx[j] );
        }
    }
}

/// Extract selected columns from a matrix.
template <typename eigenT, typename eigenTin>
void extractCols( eigenT &out,        /**< [out] selected columns */
                  const eigenTin &in, /**< [in] source matrix */
                  const std::vector<size_t> &idx /**< [in] ordered source column indices */ )
{
    out.resize( in.rows(), idx.size() );

    for( size_t i = 0; i < idx.size(); ++i )
    {
        out.col( i ) = in.col( idx[i] ); // it1->index);
    }
}

/// Select a target image's admissible reference covariance rows, columns, and images.
template <typename realT, typename eigenT, typename eigenTv, class referenceDerotT, class targetDerotT>
void collapseCovar( eigenT &cutCV,                           /**< [out] selected covariance matrix */
                    const eigenT &CV,                        /**< [in] lower-triangular covariance matrix */
                    const std::vector<realT> &sds,           /**< [in] reference Euclidean norms */
                    eigenT &rimsCut,                         /**< [out] selected vectorized references */
                    const eigenTv &rims,                     /**< [in] vectorized references */
                    const eigenTv &tims,                     /**< [in] vectorized target images */
                    int imno,                                /**< [in] target image index */
                    double dang,                             /**< [in] minimum exclusion threshold */
                    double dangMax,                          /**< [in] maximum exclusion threshold */
                    HCI::exclude excludeMethod,              /**< [in] minimum exclusion method */
                    HCI::exclude excludeMethodMax,           /**< [in] maximum exclusion method */
                    HCI::include includeMethod,              /**< [in] reference inclusion ranking method */
                    int includeRefNum,                       /**< [in] maximum ranked references retained */
                    const referenceDerotT &referenceDerotF,  /**< [in] reference derotator */
                    const targetDerotT &targetDerotF,        /**< [in] target derotator */
                    const std::vector<double> &referenceMJD, /**< [in] reference timestamps */
                    const std::vector<double> &targetMJD,    /**< [in] target timestamps */
                    eigenImage<int> &imsIncluded /**< [out] per-target inclusion flags */ )
{
    const int referenceCount = rims.cols();
    const int targetCount = tims.cols();
    if( referenceCount <= 0 || targetCount <= 0 || imno < 0 || imno >= targetCount || CV.rows() < referenceCount ||
        CV.cols() < referenceCount || sds.size() < static_cast<size_t>( referenceCount ) ||
        rims.rows() != tims.rows() || imsIncluded.rows() < targetCount || imsIncluded.cols() < referenceCount ||
        includeRefNum < 0 )
    {
        throw std::invalid_argument( "invalid KLIP covariance-selection dimensions" );
    }
    if( includeMethod != HCI::include::all && includeMethod != HCI::include::corr &&
        includeMethod != HCI::include::time && includeMethod != HCI::include::angle &&
        includeMethod != HCI::include::imno )
    {
        throw std::invalid_argument( "invalid KLIP reference inclusion method" );
    }
    if( includeMethod == HCI::include::time && includeRefNum > 0 &&
        ( referenceMJD.size() < static_cast<size_t>( referenceCount ) ||
          targetMJD.size() < static_cast<size_t>( targetCount ) ) )
    {
        throw std::invalid_argument( "time-based KLIP inclusion requires target and reference timestamps" );
    }

    std::vector<cvEntry> allidx( referenceCount );

    realT targetNorm{ 0 };
    if( includeMethod == HCI::include::corr && includeRefNum > 0 )
    {
        targetNorm = std::sqrt( tims.col( imno ).square().sum() );
    }
    for( int i = 0; i < referenceCount; ++i )
    {
        if( includeMethod == HCI::include::corr && includeRefNum > 0 )
        {
            const realT denominator = sds[i] * targetNorm;
            if( std::isfinite( denominator ) && denominator > std::numeric_limits<realT>::epsilon() )
            {
                allidx[i].cvVal = ( rims.col( i ) * tims.col( imno ) ).sum() / denominator;
                if( !std::isfinite( allidx[i].cvVal ) )
                {
                    allidx[i].cvVal = -std::numeric_limits<double>::infinity();
                }
            }
            else
            {
                allidx[i].cvVal = -std::numeric_limits<double>::infinity();
            }
        }
        else if( includeMethod == HCI::include::time && includeRefNum > 0 )
        {
            if( !std::isfinite( referenceMJD[i] ) || !std::isfinite( targetMJD[imno] ) )
            {
                throw std::invalid_argument( "time-based KLIP inclusion requires finite timestamps" );
            }
            allidx[i].cvVal = std::abs( referenceMJD[i] - targetMJD[imno] );
        }
        else if( includeMethod == HCI::include::angle && includeRefNum > 0 )
        {
            allidx[i].cvVal = std::abs( math::angleDiff<math::radiansT<realT>>( referenceDerotF.derotAngle( i ),
                                                                                targetDerotF.derotAngle( imno ) ) );
        }
        else if( includeMethod == HCI::include::imno && includeRefNum > 0 )
        {
            allidx[i].cvVal = std::abs( static_cast<long>( i ) - static_cast<long>( imno ) );
        }
    }

    if( excludeMethod == HCI::exclude::pixel || excludeMethod == HCI::exclude::angle )
    {
        for( int j = 0; j < referenceCount; ++j )
        {
            if( std::abs( math::angleDiff<math::radiansT<realT>>( referenceDerotF.derotAngle( j ),
                                                                  targetDerotF.derotAngle( imno ) ) ) <= dang )
            {
                allidx[j].included = false;
            }
        }
    }
    else if( excludeMethod == HCI::exclude::imno )
    {
        for( int j = 0; j < referenceCount; ++j )
        {
            if( fabs( (long)j - imno ) <= dang )
            {
                allidx[j].included = false;
            }
        }
    }

    if( excludeMethodMax == HCI::exclude::pixel || excludeMethodMax == HCI::exclude::angle )
    {
        for( int j = 0; j < referenceCount; ++j )
        {
            if( std::abs( math::angleDiff<math::radiansT<realT>>( referenceDerotF.derotAngle( j ),
                                                                  targetDerotF.derotAngle( imno ) ) ) > dangMax )
                allidx[j].included = false;
        }
    }
    else if( excludeMethodMax == HCI::exclude::imno )
    {
        for( int j = 0; j < referenceCount; ++j )
        {
            if( fabs( (long)j - imno ) > dangMax )
                allidx[j].included = false;
        }
    }

    if( includeMethod != HCI::include::all && includeRefNum > 0 )
    {
        std::vector<size_t> included;
        for( size_t j = 0; j < allidx.size(); ++j )
        {
            if( allidx[j].included )
            {
                included.push_back( j );
            }
        }

        if( static_cast<size_t>( includeRefNum ) < included.size() )
        {
            std::sort( included.begin(),
                       included.end(),
                       [&allidx, includeMethod]( size_t lhs, size_t rhs )
                       {
                           if( allidx[lhs].cvVal == allidx[rhs].cvVal )
                           {
                               return lhs < rhs;
                           }
                           if( includeMethod == HCI::include::corr )
                           {
                               return allidx[lhs].cvVal > allidx[rhs].cvVal;
                           }
                           return allidx[lhs].cvVal < allidx[rhs].cvVal;
                       } );
            for( size_t j = static_cast<size_t>( includeRefNum ); j < included.size(); ++j )
            {
                allidx[included[j]].included = false;
            }
        }
    }

    std::vector<size_t> keepidx;
    for( int j = 0; j < referenceCount; ++j )
    {
        imsIncluded( imno, j ) = allidx[j].included;

        if( allidx[j].included )
            keepidx.push_back( j );
    }

    // std::cerr << "  Keeping " << keepidx.size() << " reference images out of
    // "
    // << Nims << " (" << Nims-keepidx.size() << " rejected)\n";

    extractRowsAndCols( cutCV, CV, keepidx );
    extractCols( rimsCut, rims, keepidx );
}

template <typename realT, class derotFunctObj, typename evCalcT, class verboseT>
void KLIPreduction<realT, derotFunctObj, evCalcT, verboseT>::worker(
    eigenCube<realT> &rims, eigenCube<realT> &tims, imageT &cmask, std::vector<size_t> &idx, realT dang, realT dangMax )
{

    const bool isRDI = &rims != &tims;
    const derotFunctObj &referenceDerotF = isRDI ? this->m_RDIderotF : this->m_derotF;
    const std::vector<double> &referenceMJD = isRDI ? this->m_RDIimageMJD : this->m_imageMJD;
    const bool selectionActive = m_includeMethod != HCI::include::all && m_includeRefNum > 0;
    const bool exclusionActive = m_excludeMethod != HCI::exclude::none || m_excludeMethodMax != HCI::exclude::none;
    const bool angleExclusionActive =
        m_excludeMethod == HCI::exclude::pixel || m_excludeMethod == HCI::exclude::angle ||
        m_excludeMethodMax == HCI::exclude::pixel || m_excludeMethodMax == HCI::exclude::angle;
    const bool useAllReferences = !exclusionActive && !selectionActive;

    if( rims.planes() <= 0 || tims.planes() <= 0 )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidarg, "KLIP requires reference and target images" );
    }
    if( rims.rows() <= 0 || rims.cols() <= 0 || rims.rows() != tims.rows() || rims.cols() != tims.cols() ||
        static_cast<size_t>( rims.rows() ) * static_cast<size_t>( rims.cols() ) != idx.size() )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "invalid KLIP worker image dimensions" );
    }

    if( m_includeRefNum < 0 || ( m_includeMethod != HCI::include::all && m_includeMethod != HCI::include::corr &&
                                 m_includeMethod != HCI::include::time && m_includeMethod != HCI::include::angle &&
                                 m_includeMethod != HCI::include::imno ) )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidarg, "invalid KLIP reference inclusion configuration" );
    }
    if( selectionActive && m_includeMethod == HCI::include::time &&
        ( referenceMJD.size() < static_cast<size_t>( rims.planes() ) ||
          this->m_imageMJD.size() < static_cast<size_t>( tims.planes() ) ) )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidarg,
                                       "time-based KLIP inclusion requires target and reference timestamps" );
    }
    if( selectionActive && m_includeMethod == HCI::include::time )
    {
        for( int i = 0; i < rims.planes(); ++i )
        {
            if( !std::isfinite( referenceMJD[i] ) )
            {
                throw mx::exception<verboseT>( mx::error_t::invalidarg,
                                               "time-based KLIP inclusion requires finite reference timestamps" );
            }
        }
        for( int i = 0; i < tims.planes(); ++i )
        {
            if( !std::isfinite( this->m_imageMJD[i] ) )
            {
                throw mx::exception<verboseT>( mx::error_t::invalidarg,
                                               "time-based KLIP inclusion requires finite target timestamps" );
            }
        }
    }
    if( ( angleExclusionActive || ( selectionActive && m_includeMethod == HCI::include::angle ) ) &&
        ( referenceDerotF.m_angles.size() < static_cast<size_t>( rims.planes() ) ||
          this->m_derotF.m_angles.size() < static_cast<size_t>( tims.planes() ) ) )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidarg,
                                       "angle-based KLIP selection requires target and reference angles" );
    }
    if( angleExclusionActive || ( selectionActive && m_includeMethod == HCI::include::angle ) )
    {
        for( int i = 0; i < rims.planes(); ++i )
        {
            if( !std::isfinite( referenceDerotF.derotAngle( i ) ) )
            {
                throw mx::exception<verboseT>( mx::error_t::invalidarg,
                                               "angle-based KLIP selection requires finite reference angles" );
            }
        }
        for( int i = 0; i < tims.planes(); ++i )
        {
            if( !std::isfinite( this->m_derotF.derotAngle( i ) ) )
            {
                throw mx::exception<verboseT>( mx::error_t::invalidarg,
                                               "angle-based KLIP selection requires finite target angles" );
            }
        }
    }

    std::vector<realT> sds;

    imageT meanim;

    //*** First mean subtract ***//
    try
    {
        meanSubtract( rims, tims, cmask, sds );
    }
    catch( ... )
    {
        std::throw_with_nested( mx::exception<verboseT>( mx::error_t::exception, "from meanSubtract" ) );
    }

    //*** Form the reference covariance when required by the solve, selection, or diagnostics.
    imageT cv;
    const int pixelCount = rims.cube().rows();
    const int referenceCount = rims.cube().cols();
    const bool formReferenceCovariance = !useAllReferences || referenceCount <= pixelCount || m_writeDiagnostics;
    if( formReferenceCovariance )
    {
        math::eigenSYRK( cv, rims.cube() );
    }
    // eigenSYRK fills only the lower triangle; selection and full-matrix diagnostics inspect both triangles.
    if( !useAllReferences || m_writeDiagnostics )
    {
        for( int column = 1; column < cv.cols(); ++column )
        {
            for( int row = 0; row < column; ++row )
            {
                cv( row, column ) = cv( column, row );
            }
        }
    }
    if( m_writeDiagnostics )
    {
        writeDiagnostic( "cv.fits", cv );
    }
    ipc::ompLoopWatcher<> status( tims.planes(), std::cerr );

    // Pre-calculate KL images once if we are exclude none OR IF RDI
    imageT master_klims;
    imageT master_projMat;

    eigenImage<realT> rrMask;
    if( m_rightReason )
    {
        rrMask.resize( tims.rows(), tims.rows() );
        rrMask.setConstant( 1 );

        imageT rrmim;
        rrmim.resize( this->m_Nrows, this->m_Ncols );
        for( int rr = 0; rr < rrMask.rows(); ++rr )
        {
            rrmim.setConstant( 1 );

            // Calculate the 2D coords of this pixel
            int jj = idx[rr] / this->m_Ncols;
            int ii = idx[rr] - this->m_Ncols * jj;

            maskCircle( rrmim, ii, jj, m_rightReasonRadius, 0, 0 );

            // Extract the 2D r.r. mask into the row-vector image.
            for( int cc = 0; cc < rrMask.cols(); ++cc )
            {
                rrMask( rr, cc ) = rrmim.data()[idx[cc]];
            }
        }

        writeDiagnostic( "rrMask.fits", rrMask );
    }

    if( useAllReferences )
    {
        double teigenv{ 0 };
        double tklim{ 0 };
        math::syevrMem<evCalcT> workspace;

        const MXLAPACK_INT solverStatus =
            calcKLModesAdaptive<evCalcT>( master_klims, cv, rims.cube(), m_maxNmodes, &workspace, &teigenv, &tklim );
        if( solverStatus != 0 )
        {
            throw mx::exception<verboseT>( mx::error_t::lapackerr,
                                           "KLIP eigensolver failed with status " + std::to_string( solverStatus ) );
        }

        if( m_rightReason )
        {
            master_projMat = ( master_klims.matrix().transpose() * master_klims.matrix() ).array();
            writeDiagnostic( "projMat.fits", master_projMat );
            master_projMat *= rrMask;
            writeDiagnostic( "projMatrr.fits", master_projMat );
        }

        m_algorithmTiming.eigensolveWorkerSeconds += teigenv;
        m_algorithmTiming.modeWorkerSeconds += tklim;
    }

    const imageT &sharedMasterKlims = master_klims;
    const imageT &sharedMasterProjMat = master_projMat;
    std::exception_ptr workerException;
    std::atomic<bool> workerFailed{ false };
    double parallelEigenSeconds{ 0 };
    double parallelModeSeconds{ 0 };
    double parallelPsfSeconds{ 0 };

    // clang-format off
    #pragma omp parallel reduction( + : parallelEigenSeconds, parallelModeSeconds, parallelPsfSeconds )
    // clang-format on
    {
        imageT cfs; // The coefficients
        imageT psf;
        imageT rims_cut;
        imageT cv_cut;
        imageT localKlims;
        imageT localProjMat;

        math::syevrMem<evCalcT> mem;

        // clang-format off
        #pragma omp for
        // clang-format on
        for( int imno = 0; imno < tims.planes(); ++imno )
        {
            if( workerFailed.load( std::memory_order_relaxed ) )
            {
                continue;
            }

            try
            {
                const imageT *activeKlims = &sharedMasterKlims;
                const imageT *activeProjMat = &sharedMasterProjMat;
                if( !useAllReferences )
                {
                    collapseCovar<realT>( cv_cut,
                                          cv,
                                          sds,
                                          rims_cut,
                                          rims.asVectors(),
                                          tims.asVectors(),
                                          imno,
                                          dang,
                                          dangMax,
                                          this->m_excludeMethod,
                                          this->m_excludeMethodMax,
                                          this->m_includeMethod,
                                          this->m_includeRefNum,
                                          referenceDerotF,
                                          this->m_derotF,
                                          referenceMJD,
                                          this->m_imageMJD,
                                          m_imsIncluded );
                    if( rims_cut.cols() <= 0 )
                    {
                        throw mx::exception<verboseT>( mx::error_t::invalidarg,
                                                       "KLIP target " + std::to_string( imno ) +
                                                           " has no admissible references" );
                    }

                    double teigenv{ 0 };
                    double tklim{ 0 };
                    const MXLAPACK_INT solverStatus = calcKLModesAdaptive<evCalcT>( localKlims,
                                                                                    cv_cut,
                                                                                    rims_cut,
                                                                                    m_maxNmodes,
                                                                                    &mem,
                                                                                    &teigenv,
                                                                                    &tklim );
                    if( solverStatus != 0 )
                    {
                        throw mx::exception<verboseT>( mx::error_t::lapackerr,
                                                       "KLIP eigensolver failed for target " + std::to_string( imno ) +
                                                           " with status " + std::to_string( solverStatus ) );
                    }

                    if( m_rightReason )
                    {
                        localProjMat = ( localKlims.matrix().transpose() * localKlims.matrix() ).array();
                        localProjMat *= rrMask;
                    }

                    parallelEigenSeconds += teigenv;
                    parallelModeSeconds += tklim;
                    activeKlims = &localKlims;
                    activeProjMat = &localProjMat;
                }

                cfs.resize( 1, activeKlims->rows() );
                if( cfs.size() <= 0 )
                {
                    throw mx::exception<verboseT>( mx::error_t::sizeerr, "KLIP eigensolver returned no modes" );
                }

                const double psfStart = sys::get_curr_time();

                if( !m_rightReason )
                {
                    for( int j = 0; j < cfs.size(); ++j )
                    {
                        cfs( j ) = activeKlims->row( j ).matrix().dot( tims.cube().col( imno ).matrix() );
                    }

                    for( size_t mode_i = 0; mode_i < m_Nmodes.size(); ++mode_i )
                    {
                        psf = cfs( cfs.size() - 1 ) * activeKlims->row( cfs.size() - 1 );

                        // Count down, since eigenvalues are returned in increasing
                        // order
                        //   handle case where cfs.size() < m_Nmodes[mode_i], i.e.
                        //   when more modes than images.
                        for( int j = cfs.size() - 2; j >= cfs.size() - m_Nmodes[mode_i] && j >= 0; --j )
                        {
                            psf += cfs( j ) * activeKlims->row( j );
                        }

                        insertImageRegion( this->m_psfsub[mode_i].cube().col( imno ),
                                           tims.cube().col( imno ) - psf.transpose(),
                                           idx );
                    }
                }
                else
                {
                    psf = activeProjMat->matrix() * tims.cube().col( imno ).matrix();
                    insertImageRegion( this->m_psfsub[0].cube().col( imno ), tims.cube().col( imno ) - psf, idx );
                }

                parallelPsfSeconds += sys::get_curr_time() - psfStart;
                status.incrementAndOutputStatus();
            }
            catch( ... )
            {
                workerFailed.store( true, std::memory_order_relaxed );
                // clang-format off
                #pragma omp critical( KLIPreduction_worker_exception )
                // clang-format on
                {
                    if( workerException == nullptr )
                    {
                        workerException = std::current_exception();
                    }
                }
            }

        } // for imno
    } // openmp parrallel

    m_algorithmTiming.eigensolveWorkerSeconds += parallelEigenSeconds;
    m_algorithmTiming.modeWorkerSeconds += parallelModeSeconds;
    m_algorithmTiming.projectionWorkerSeconds += parallelPsfSeconds;
    if( workerException != nullptr )
    {
        std::rethrow_exception( workerException );
    }
}

template <typename realT, class derotFunctObj, typename evCalcT, class verboseT>
void KLIPreduction<realT, derotFunctObj, evCalcT, verboseT>::appendReductionHeader( fitsHeaderT &head )
{
    head.append( "", fits::fitsCommentType(), "----------------------------------------" );
    head.append( "", fits::fitsCommentType(), "mx::KLIPreduction parameters:" );
    head.append( "", fits::fitsCommentType(), "----------------------------------------" );

    head.append( "MEAN SUB METHOD", HCI::meanSubToStr<verboseT>( m_meanSubMethod ), "PCA mean subtraction method" );
    head.append( "PIXTS NORM METHOD", HCI::pixelTSNormToStr<verboseT>( m_pixelTSNormMethod ), "Pixel TS norm method" );

    std::stringstream str;

    if( m_Nmodes.size() > 0 )
    {
        for( size_t nm = 0; nm < m_Nmodes.size() - 1; ++nm )
            str << m_Nmodes[nm] << ",";
        str << m_Nmodes[m_Nmodes.size() - 1];
        head.template append<char *>( "NMODES", (char *)str.str().c_str(), "number of modes" );
    }

    head.template append<bool>( "RIGHT REASON", m_rightReason, "whether or not the right reason mask is used" );
    if( m_rightReason )
    {
        head.template append<realT>( "RIGHT REASON RADIUS", m_rightReasonRadius, "radius of the right reason mask" );
    }

    if( m_minr.size() > 0 )
    {
        str.str( "" );
        for( size_t nm = 0; nm < m_minr.size() - 1; ++nm )
            str << m_minr[nm] << ",";
        str << m_minr[m_minr.size() - 1];
        head.template append<char *>( "REGMINR", (char *)str.str().c_str(), "region inner edge(s)" );
    }

    if( m_maxr.size() > 0 )
    {
        str.str( "" );
        for( size_t nm = 0; nm < m_maxr.size() - 1; ++nm )
            str << m_maxr[nm] << ",";
        str << m_maxr[m_maxr.size() - 1];
        head.template append<char *>( "REGMAXR", (char *)str.str().c_str(), "region outer edge(s)" );
    }

    if( m_minq.size() > 0 )
    {
        str.str( "" );
        for( size_t nm = 0; nm < m_minq.size() - 1; ++nm )
            str << m_minq[nm] << ",";
        str << m_minq[m_minq.size() - 1];
        head.template append<char *>( "REGMINQ", (char *)str.str().c_str(), "region minimum angle(s)" );
    }

    if( m_maxq.size() > 0 )
    {
        str.str( "" );
        for( size_t nm = 0; nm < m_maxq.size() - 1; ++nm )
            str << m_maxq[nm] << ",";
        str << m_maxq[m_maxq.size() - 1];
        head.template append<char *>( "REGMAXQ", (char *)str.str().c_str(), "region maximum angle(s)" );
    }

    head.template append<std::string>( "EXMTHDMN",
                                       HCI::excludeToStr<verboseT>( m_excludeMethod ),
                                       "exclusion method (min)" );

    head.template append<realT>( "MINDPX", m_minDPx, "minimum delta (units based on EXMTHDMN)" );

    head.template append<std::string>( "EXMTHDMX",
                                       HCI::excludeToStr<verboseT>( m_excludeMethodMax ),
                                       "exclusion method (max)" );
    head.template append<realT>( "MAXDPX", m_maxDPx, "maximum delta (units based on EXMTHDMX)" );

    head.template append<std::string>( "INMTHDMX", HCI::includeToStr<verboseT>( m_includeMethod ), "inclusion method" );
    head.template append<int>( "INCLREFN", m_includeRefNum, "number of images included by INMTHDMX" );
}

template <typename realT, class derotFunctObj, typename evCalcT, class verboseT>
int KLIPreduction<realT, derotFunctObj, evCalcT, verboseT>::finalProcess()
{
    fitsHeaderT algorithmHeader;
    const fitsHeaderT *algorithmHeaderPointer = nullptr;

    if( this->m_doWriteFinim == true || this->m_doOutputPSFSub == true )
    {
        appendReductionHeader( algorithmHeader );
        algorithmHeaderPointer = &algorithmHeader;
    }

    return this->ADIobservation<realT, derotFunctObj, verboseT>::finalProcess( algorithmHeaderPointer );
}

template <typename realT, class derotFunctObj, typename evCalcT, class verboseT>
int KLIPreduction<realT, derotFunctObj, evCalcT, verboseT>::processPSFSub( const std::string &directory,
                                                                           const std::string &prefix,
                                                                           const std::string &extension )
{
    this->readPSFSub( directory, prefix, extension );

    if( this->m_heads.empty() || this->m_heads.front().count( "NMODES" ) == 0 )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidarg, "saved KLIP reductions require an NMODES header card" );
    }

    std::vector<int> modeCounts;
    ioutils::parseStringVector( modeCounts, this->m_heads.front()["NMODES"].String(), ',' );
    if( modeCounts.size() != this->m_psfsub.size() ||
        std::any_of( modeCounts.begin(), modeCounts.end(), []( int count ) { return count <= 0; } ) )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidarg,
                                       "saved KLIP NMODES metadata does not match the reduction collection" );
    }
    m_Nmodes = std::move( modeCounts );

    if( !this->m_weightFile.empty() )
    {
        this->readWeights();
    }

    return finalProcess();
}

/*
template <typename realT, class derotFunctObj, typename evCalcT, class verboseT>
int KLIPreduction<realT, derotFunctObj, evCalcT, verboseT>::processPSFSub( const std::string &dir,
                                                                    const std::string &prefix,
                                                                    const std::string &ext )

{
    std::cerr << "Beginning PSF Subtracted Image Processing\n";

    // Load first file to condigure based on its header.
    std::vector<std::string> flist = ioutils::getFileNames( dir, prefix, "000", ext );

    fits::fitsHeader fh;
    eigenImage<realT> im;
    fits::fitsFile<realT> ff;

    ff.read( im, fh, flist[0] );

    if( fh.count( "MEANSUBM" ) == 0 )
    {
        mxError( "KLIPReduction", MXE_PARAMNOTSET, "MEANSUBM not found in FITS header." );
        return -1;
    }
    m_meanSubMethod = HCI::meanSubMethodFmStr( fh["MEANSUBM"].String() );
    std::cerr << "meanSubMethod: " << HCI::meanSubMethodStr( m_meanSubMethod ) << "\n";

    if( fh.count( "NMODES" ) == 0 )
    {
        mxError( "KLIPReduction", MXE_PARAMNOTSET, "NMODES not found in FITS header." );
        return -1;
    }
    ioutils::parseStringVector( m_Nmodes, fh["NMODES"].String(), "," );
    if( m_Nmodes.size() == 0 )
    {
        mxError( "KLIPReduction", MXE_PARSEERR, "NMODES vector did not parse correctly." );
        return -1;
    }
    std::cerr << "nModes: " << fh["NMODES"].String() << "\n";

    /* -- REGMINR --
    if( fh.count( "REGMINR" ) == 0 )
    {
        mxError( "KLIPReduction", MXE_PARAMNOTSET, "REGMINR not found in FITS header." );
        return -1;
    }
    ioutils::parseStringVector( m_minr, fh["REGMINR"].String(), "," );
    if( m_minr.size() == 0 )
    {
        mxError( "KLIPReduction", MXE_PARSEERR, "REGMINR vector did not parse correctly." );
        return -1;
    }
    std::cerr << "minr: " << fh["REGMINR"].String() << "\n";

    /* -- REGMAXR --
    if( fh.count( "REGMAXR" ) == 0 )
    {
        mxError( "KLIPReduction", MXE_PARAMNOTSET, "REGMAXR not found in FITS header." );
        return -1;
    }
    ioutils::parseStringVector( m_maxr, fh["REGMAXR"].String(), "," );
    if( m_maxr.size() == 0 )
    {
        mxError( "KLIPReduction", MXE_PARSEERR, "REGMAXR vector did not parse correctly." );
        return -1;
    }
    std::cerr << "minr: " << fh["REGMAXR"].String() << "\n";

    /* -- REGMINQ --
    if( fh.count( "REGMINQ" ) == 0 )
    {
        mxError( "KLIPReduction", MXE_PARAMNOTSET, "REGMINQ not found in FITS header." );
        return -1;
    }
    ioutils::parseStringVector( m_minq, fh["REGMINQ"].String(), "," );
    if( m_minq.size() == 0 )
    {
        mxError( "KLIPReduction", MXE_PARSEERR, "REGMINQ vector did not parse correctly." );
        return -1;
    }
    std::cerr << "minq: " << fh["REGMINQ"].String() << "\n";

    /* -- REGMAXR --
    if( fh.count( "REGMAXR" ) == 0 )
    {
        mxError( "KLIPReduction", MXE_PARAMNOTSET, "REGMAXR not found in FITS header." );
        return -1;
    }
    ioutils::parseStringVector( m_maxq, fh["REGMAXR"].String(), "," );
    if( m_maxq.size() == 0 )
    {
        mxError( "KLIPReduction", MXE_PARSEERR, "REGMAXR vector did not parse correctly." );
        return -1;
    }
    std::cerr << "minr: " << fh["REGMAXR"].String() << "\n";

    if( fh.count( "EXMTHDMN" ) == 0 )
    {
        mxError( "KLIPReduction", MXE_PARAMNOTSET, "EXMTHDMN not found in FITS header." );
        return -1;
    }
    m_excludeMethod = HCI::excludeMethodFmStr( fh["EXMTHDMN"].String() );
    std::cerr << "excludeMethod: " << HCI::excludeMethodStr( m_excludeMethod ) << "\n";

    if( fh.count( "MINDPX" ) == 0 )
    {
        mxError( "KLIPReduction", MXE_PARAMNOTSET, "MINDPX not found in FITS header." );
        return -1;
    }
    m_minDPx = fh["MINDPX"].value<realT>();
    std::cerr << "minDPx: " << m_minDPx << "\n";

    if( fh.count( "EXMTHDMX" ) == 0 )
    {
        mxError( "KLIPReduction", MXE_PARAMNOTSET, "EXMTHDMX not found in FITS header." );
        return -1;
    }
    m_excludeMethodMax = HCI::excludeMethodFmStr( fh["EXMTHDMX"].String() );
    std::cerr << "excludeMethodMax: " << HCI::excludeMethodStr( m_excludeMethodMax ) << "\n";

    if( fh.count( "MAXDPX" ) == 0 )
    {
        mxError( "KLIPReduction", MXE_PARAMNOTSET, "MAXDPX not found in FITS header." );
        return -1;
    }
    m_maxDPx = fh["MAXDPX"].value<realT>();
    std::cerr << "maxDPx: " << m_maxDPx << "\n";

    if( fh.count( "INMTHDMX" ) == 0 )
    {
        mxError( "KLIPReduction", MXE_PARAMNOTSET, "INMTHDMX not found in FITS header." );
        return -1;
    }
    m_includeMethod = HCI::includeMethodFmStr( fh["INMTHDMX"].String() );
    std::cerr << "includeMethod: " << HCI::includeMethodStr( m_includeMethod ) << "\n";

    if( fh.count( "INCLREFN" ) == 0 )
    {
        mxError( "KLIPReduction", MXE_PARAMNOTSET, "INCLREFN not found in FITS header." );
        return -1;
    }
    m_includeRefNum = fh["INCLREFN"].value<int>();
    std::cerr << "includedRefNum: " << m_includeRefNum << "\n";

    this->m_skipPreProcess = true;

    this->m_keywords.clear();

    //this->readPSFSub( dir, prefix, ext, m_Nmodes.size() );

    finalProcess();

    return 0;
}*/

///@}

template <typename realT, class verboseT>
class ADIDerotator;

extern template struct KLIPreduction<float, ADIDerotator<float, verbose::vv>, double, verbose::vv>;

} // namespace improc
} // namespace mx

#endif //__KLIPreduction_hpp__
