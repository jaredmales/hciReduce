/** \file KLIPreduction.hpp
 * \author Jared R. Males
 * \brief Declarations and definitions for an implementation of the
 * Karhunen-Loeve Image Processing (KLIP) algorithm. \ingroup hc_imaging_files
 * \ingroup image_processing_files
 *
 */

#ifndef __KLIPreduction_hpp__
#define __KLIPreduction_hpp__

#include <map>
#include <vector>

#include <omp.h>

#include <mx/ipc/ompLoopWatcher.hpp>
#include <mx/math/eigenLapack.hpp>
#include <mx/math/geo.hpp>
#include <mx/sigproc/gramSchmidt.hpp>
using namespace mx::sigproc;

#include "ADIobservation.hpp"

namespace mx
{
namespace improc
{

// double t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, tf;
// double dcut, dscv, dklims, dgemm, dsyevr, dcfs, drot, dcombo, dread;

/// An implementation of the Karhunen-Loeve Image Processing (KLIP) algorithm.
/** KLIP\cite soummer_2012 is a principle components analysis (PCA) based
 * technique for PSF estimation.
 *
 *
 * \tparam _realT  is the floating point type in which to do calculations
 * \tparam _derotFunctObj the ADIobservation derotator class.
 * \tparam _evCalcT the real type in which to do eigen-decomposition.  Should
 * generally be double for stable results. \ingroup hc_imaging
 */
template <typename _realT, class _derotFunctObj, typename _evCalcT, class verboseT>
struct KLIPreduction : public ADIobservation<_realT, _derotFunctObj, verboseT>
{
    typedef _realT realT;
    typedef _evCalcT evCalcT;

    typedef Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> imageT;

    /// Default c'tor
    KLIPreduction();

    virtual ~KLIPreduction();

    void setupConfig( mx::app::appConfigurator &config );

    void loadConfig( mx::app::appConfigurator &config );

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

    /// Number of reference images to include in the covariance matrix
    /** If > 0, then at most this many images, determined by highest
     * cross-correlation, are included. This is determined after
     * rotational/image-number exclusion. If == 0, then all reference images are
     * included.
     */
    int m_includeRefNum{ 0 };

    /// Controls how number of included images is calculated.
    /** The number of included images is calculated after exclusion is complete.
     * Can have the following values:
     * - <b>HCI::includeAll</b> = all remaining images are included [default]
     * - <b>HCI::includeCorr</b> = the m_includeRefNum of the remaining images
     * which are most correlated with the target are included
     * - <b>HCI::includeTime</b> = the m_includeRefNum of the remaining images
     * which are closest in time to the target are included
     * - <b>HCI::includeAngle</b> = the m_includeRefNum of the remaining images
     * which are closest in angle to the target are included
     * - <b>HCI::includeImno</b> = the m_includeRefNum of the remaining images
     * which are closest in image number to the target are included
     */
    HCI::include m_includeMethod { HCI::include::all };

    eigenImage<int> m_imsIncluded;

    bool m_rightReason = false;

    realT m_rightReasonRadius = 2.5;

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

    void worker( eigenCube<realT> &rims,
                 eigenCube<realT> &tims,
                 imageT &cmask,
                 std::vector<size_t> &idx,
                 realT dang,
                 realT dangMax );

    int finalProcess();

    // int processPSFSub( const std::string &dir, const std::string &prefix, const std::string &ext );

    double t_worker_begin{ 0 };
    double t_worker_end{ 0 };

    double t_eigenv{ 0 };
    double t_klim{ 0 };
    double t_psf{ 0 };

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
        printf( "    KLIP algorithm: %f elapsed real sec\n", this->t_worker_end - this->t_worker_begin );
        double klip_cpu = this->t_eigenv + this->t_klim + this->t_psf;
        printf( "      EigenDecomposition %f cpu sec (%f%%)\n", this->t_eigenv, this->t_eigenv / klip_cpu * 100 );
        printf( "      KL image calc %f cpu sec (%f%%)\n", this->t_klim, this->t_klim / klip_cpu * 100 );
        printf( "      PSF calc/sub %f cpu sec (%f%%)\n", this->t_psf, this->t_psf / klip_cpu * 100 );
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
                "int",
                "The method of pixel time-series normalization for PCA: none or rms." );

    config.add( "klip.includeRefNum",
                "",
                "klip.includeRefNum",
                mx::app::argType::Required,
                "klip",
                "includeRefNum",
                false,
                "int",
                "The number of references to include, based on correlation." );

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
        std::throw_with_nested(mx::exception<verboseT>(mx::error_t::exception));
    }

    config( ms, "klip.meanSubMethod" );

    try
    {
        m_meanSubMethod = HCI::meanSubFmStr<verboseT>( ms );
    }
    catch( ... )
    {
        std::throw_with_nested(mx::exception<verboseT>(mx::error_t::invalidconfig, "klip.meanSubMethod is not valid"));
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
    config( m_Nmodes, "klip.Nmodes" );
    config( m_rightReason, "klip.rightReason" );
    config( m_rightReasonRadius, "klip.rrRadius" );
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

    norms.resize( rims.planes() );

    int maskPix = 0;

    if( cmask.rows() > 0 && cmask.cols() > 0 )
    {
        if( cmask.rows() != rims.rows() || cmask.cols() != rims.cols() )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidarg,
                                           "mask does not have same size as reference images" );
        }

        if( cmask.rows() != tims.rows() || cmask.cols() != tims.cols() )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidarg, "mask does not have same size as target images" );
        }

        maskPix = cmask.sum();
    }

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

                norms[n] = ( ( rims.image( n ) - immean ) * cmask ).square().sum();
            }
            else
            {
                immean = rims.image( n ).mean();

                norms[n] = ( rims.image( n ) - immean ).square().sum();
            }
        }

        if( &tims != &rims )
        {
            if( m_meanSubMethod == HCI::meanSub::meanImage )
            {
                tims.mean( mean );
            }
            else if( m_meanSubMethod == HCI::meanSub::medianImage )
            {
                tims.median( mean );
            }

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
                mean = rims.image( n ).mean();
            }
            else if( m_meanSubMethod == HCI::meanSub::imageMedian )
            {
                mean = imageMedian( rims.image( n ), &work );
            }

            rims.image( n ) -= mean;

            if( maskPix > 0 )
            {
                rims.image( n ) *= cmask;
                immean = rims.image( n ).sum() / maskPix;
                norms[n] = ( ( rims.image( n ) - immean ) * cmask ).square().sum();
            }
            else
            {
                // Because we might not have used the mean, we need to re-mean to
                //  make this the standard deviation
                immean = rims.image( n ).mean();
                norms[n] = ( rims.image( n ) - immean ).square().sum();
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

    if( m_pixelTSNormMethod == HCI::pixelTSNorm::rmsSigmaClipped )
    {
        throw mx::exception<verboseT>( mx::error_t::notimpl,
                                       "pixelTSNormMethod is rmsSigmaClipped, "
                                       "which is not implemented" );
    }

    if( m_pixelTSNormMethod != HCI::pixelTSNorm::none )
    {
        std::cerr << "normalizing pixels\n";
        std::vector<realT> pixs( rims.planes() );

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

                // We bother to load a vector in prep to add sigma clipping later.
                for( int pp = 0; pp < rims.planes(); ++pp )
                {
                    pixs[pp] = rims.image( pp )( rr, cc );
                }

                realT sd = sqrt( math::vectorVariance( pixs ) );

                for( int pp = 0; pp < rims.planes(); ++pp )
                {
                    rims.image( pp )( rr, cc ) /= sd;
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
    this->t_begin = sys::get_curr_time();

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

    std::cerr << "Beginning\n";

    if( this->m_imSize == 0 )
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

    if( this->m_preProcess_only && !this->m_skipPreProcess )
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

    m_imsIncluded.resize( this->m_Nims, this->m_Nims );
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

        // Create storage for the R-ims and psf-subbed Ims
        eigenCube<realT> tims( idx.size(), 1, this->m_Nims );

        //------If doing RDI, create bims
        eigenCube<realT> rims;

        //------Get the mask cutout too
        imageT cmask;

        if( this->m_refIms.planes() > 0 )
        {
            rims.resize( idx.size(), 1, this->m_Nims );
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

        if( m_minDPx < 0 )
        {
            m_excludeMethod = HCI::exclude::none;
        }
        if( m_maxDPx < 0 )
        {
            m_excludeMethodMax = HCI::exclude::none;
        }

        //------- If doing RDI, excludeMethod and excludeMethodMax must be none!
        if( this->m_refIms.planes() > 0 )
        {
            m_excludeMethod = HCI::exclude::none;
            m_excludeMethodMax = HCI::exclude::none;
        }

        if( m_excludeMethod == HCI::exclude::pixel )
        {
            dang = fabs( atan( m_minDPx / minr[regno] ) );
        }
        else if( m_excludeMethod == HCI::exclude::angle )
        {
            dang = math::dtor( m_minDPx );
        }
        else if( m_excludeMethod == HCI::exclude::imno )
        {
            dang = m_minDPx;
        }

        if( m_excludeMethodMax == HCI::exclude::pixel )
        {
            dangMax = fabs( atan( m_maxDPx / minr[regno] ) );
        }
        else if( m_excludeMethodMax == HCI::exclude::angle )
        {
            dangMax = math::dtor( m_maxDPx );
        }
        else if( m_excludeMethodMax == HCI::exclude::imno )
        {
            dangMax = m_maxDPx;
        }

        //------- If doing RDI, call this with rims, bims
        //*** Dispatch the work
        if( this->m_refIms.planes() > 0 ) // RDI
        {
            worker( rims, tims, cmask, idx, dang, dangMax );
        }
        else // ADI
        {
            worker( tims, tims, cmask, idx, dang, dangMax );
        }
    }

    fits::fitsFile<int> ffii;
    ffii.write( "imsIncluded.fits", m_imsIncluded );

    if( finalProcess() < 0 )
    {
        std::cerr << "Error in final processing\n";
    }

    this->t_end = sys::get_curr_time();

    dump_times();

    return 0;
}

struct cvEntry
{
    int index;
    double cvVal;
    double angle;
    bool included{ true };
};

template <typename eigenT, typename eigenTin>
void extractRowsAndCols( eigenT &out, const eigenTin &in, const std::vector<size_t> &idx )
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

template <typename eigenT, typename eigenTin>
void extractCols( eigenT &out, const eigenTin &in, const std::vector<size_t> &idx )
{
    out.resize( in.rows(), idx.size() );

    for( size_t i = 0; i < idx.size(); ++i )
    {
        out.col( i ) = in.col( idx[i] ); // it1->index);
    }
}

template <typename realT, typename eigenT, typename eigenTv, class derotFunctObj>
void collapseCovar( eigenT &cutCV,
                    const eigenT &CV,
                    const std::vector<realT> &sds,
                    eigenT &rimsCut,
                    const eigenTv &rims,
                    int imno,
                    double dang,
                    double dangMax,
                    int Nims,
                    HCI::exclude excludeMethod,
                    HCI::exclude excludeMethodMax,
                    int includeRefNum,
                    const derotFunctObj &derotF,
                    eigenImage<int> &imsIncluded )
{
    std::vector<cvEntry> allidx( Nims );

    // Initialize the vector of cvEntries
    for( int i = 0; i < Nims; ++i )
    {
        allidx[i].index = i;
        allidx[i].angle = derotF.derotAngle( i );

        // CV is lower-triangular
        if( i <= imno )
        {
            allidx[i].cvVal = CV( imno, i ) / ( sds[i] * sds[imno] );
        }
        else
        {
            allidx[i].cvVal = CV( i, imno ) / ( sds[i] * sds[imno] );
        }
    }

    if( excludeMethod == HCI::exclude::pixel || excludeMethod == HCI::exclude::angle )
    {
        for( size_t j = 0; j < Nims; ++j )
        {
            if( fabs( math::angleDiff<math::radiansT<realT>>( derotF.derotAngle( j ), derotF.derotAngle( imno ) ) ) <=
                dang )
            {
                allidx[j].included = false;
            }
        }
    }
    else if( excludeMethod == HCI::exclude::imno )
    {
        for( size_t j = 0; j < Nims; ++j )
        {
            if( fabs( (long)j - imno ) <= dang )
            {
                allidx[j].included = false;
            }
        }
    }

    if( excludeMethodMax == HCI::exclude::pixel || excludeMethodMax == HCI::exclude::angle )
    {
        for( size_t j = 0; j < Nims; ++j )
        {
            if( fabs( math::angleDiff<math::radiansT<realT>>( derotF.derotAngle( j ), derotF.derotAngle( imno ) ) ) >
                dangMax )
                allidx[j].included = false;
        }
    }
    else if( excludeMethodMax == HCI::exclude::imno )
    {
        for( size_t j = 0; j < Nims; ++j )
        {
            if( fabs( (long)j - imno ) > dangMax )
                allidx[j].included = false;
        }
    }

    if( includeRefNum > 0 && (size_t)includeRefNum < allidx.size() )
    {
        long kept = 0;
        for( size_t j = 0; j < Nims; ++j )
        {
            if( allidx[j].included == true )
            {
                ++kept;
            }
        }

        // Get a vector for sorting
        std::vector<realT> cvVal;
        cvVal.resize( kept );
        size_t k = 0;
        for( size_t j = 0; j < Nims; ++j )
        {
            if( allidx[j].included == true )
            {
                cvVal[k] = allidx[j].cvVal;
                ++k;
            }
        }

        // Partially sort the correlation values
        std::nth_element( cvVal.begin(), cvVal.begin() + ( kept - includeRefNum ), cvVal.end() );

        realT mincorr = cvVal[kept - includeRefNum];
        // std::cerr << "    Minimum correlation: " << mincorr << "\n";

        for( size_t j = 0; j < Nims; ++j )
        {
            if( allidx[j].cvVal < mincorr )
            {
                allidx[j].included = false;
            }
        }
    }

    std::vector<size_t> keepidx;
    for( size_t j = 0; j < Nims; ++j )
    {
        imsIncluded( imno, j ) = allidx[j].included;

        if( allidx[j].included )
            keepidx.push_back( j );
    }

    // std::cerr << "  Keeping " << keepidx.size() << " reference images out of
    // "
    // << Nims << " (" << Nims-keepidx.size() << " rejected)\n";

    if( keepidx.size() == 0 )
    {
        std::cerr << "\n\n" << imno << "\n\n";
    }

    extractRowsAndCols( cutCV, CV, keepidx );
    extractCols( rimsCut, rims, keepidx );
}

template <typename realT, class derotFunctObj, typename evCalcT, class verboseT>
void KLIPreduction<realT, derotFunctObj, evCalcT, verboseT>::worker(
    eigenCube<realT> &rims, eigenCube<realT> &tims, imageT &cmask, std::vector<size_t> &idx, realT dang, realT dangMax )
{

    t_worker_begin = sys::get_curr_time();

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

    //*** Form lower-triangle covariance matrix
    imageT cv;

    math::eigenSYRK( cv, rims.cube() );

    fits::fitsFile<realT> ff;
    ff.write( "cv.fits", cv );
    ipc::ompLoopWatcher<> status( this->m_Nims, std::cerr );

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

        ff.write( "rrMask.fits", rrMask );
    }

    if( m_excludeMethod == HCI::exclude::none && m_excludeMethodMax == HCI::exclude::none && m_includeRefNum == 0 )
    {
        double teigenv;
        double tklim;

        math::calcKLModes<double>( master_klims, cv, rims.cube(), m_maxNmodes, nullptr, &teigenv, &tklim );

        if( m_rightReason )
        {
            master_projMat = ( master_klims.matrix().transpose() * master_klims.matrix() ).array();
            ff.write( "projMat.fits", master_projMat );
            master_projMat *= rrMask;
            ff.write( "projMatrr.fits", master_projMat );
        }

        t_eigenv += teigenv;
        t_klim += tklim;
    }

    // clang-format off
    #pragma omp parallel // num_threads(20) // clang-format on
    {
        // We need local copies for each thread.  Only way this works, for
        // whatever reason.
        imageT cfs; // The coefficients
        imageT psf;
        imageT rims_cut;
        imageT cv_cut;
        imageT klims;
        imageT projMat;

        math::syevrMem<evCalcT> mem;

        if( m_excludeMethod == HCI::exclude::none && m_excludeMethodMax == HCI::exclude::none &&
            m_includeRefNum == 0 ) // OR RDI
        {
            klims = master_klims;
            if( m_rightReason )
            {
                projMat = master_projMat;
            }
        }

        // clang-format off
        #pragma omp for // clang-format on
        for( int imno = 0; imno < this->m_Nims; ++imno )
        {
            if( m_excludeMethod != HCI::exclude::none || m_excludeMethodMax != HCI::exclude::none ||
                m_includeRefNum != 0 )
            {
                collapseCovar<realT>( cv_cut,
                                      cv,
                                      sds,
                                      rims_cut,
                                      rims.asVectors(),
                                      imno,
                                      dang,
                                      dangMax,
                                      this->m_Nims,
                                      this->m_excludeMethod,
                                      this->m_excludeMethodMax,
                                      this->m_includeRefNum,
                                      this->m_derotF,
                                      m_imsIncluded );

                /**** Now calculate the K-L Images ****/
                double teigenv, tklim;
                math::calcKLModes( klims, cv_cut, rims_cut, m_maxNmodes, &mem, &teigenv, &tklim );

                if( m_rightReason )
                {
                    projMat = ( klims.matrix().transpose() * klims.matrix() ).array();
                    projMat *= rrMask;
                }

                t_eigenv += teigenv;
                t_klim += tklim;
            }
            cfs.resize( 1, klims.rows() );

            double t0 = sys::get_curr_time();

            if( !m_rightReason )
            {
                for( int j = 0; j < cfs.size(); ++j )
                {
                    cfs( j ) = klims.row( j ).matrix().dot( tims.cube().col( imno ).matrix() );
                }

                for( size_t mode_i = 0; mode_i < m_Nmodes.size(); ++mode_i )
                {
                    psf = cfs( cfs.size() - 1 ) * klims.row( cfs.size() - 1 );

                    // Count down, since eigenvalues are returned in increasing
                    // order
                    //   handle case where cfs.size() < m_Nmodes[mode_i], i.e.
                    //   when more modes than images.
                    for( int j = cfs.size() - 2; j >= cfs.size() - m_Nmodes[mode_i] && j >= 0; --j )
                    {
                        psf += cfs( j ) * klims.row( j );
                    }

                    // #pragma omp critical
                    insertImageRegion( this->m_psfsub[mode_i].cube().col( imno ),
                                       tims.cube().col( imno ) - psf.transpose(),
                                       idx );
                }
            }
            else
            {
                psf = projMat.matrix() * tims.cube().col( imno ).matrix();
                insertImageRegion( this->m_psfsub[0].cube().col( imno ), tims.cube().col( imno ) - psf, idx );
            }

            t_psf += ( sys::get_curr_time() - t0 ); /// omp_get_num_threads();

            status.incrementAndOutputStatus();

        } // for imno
    } // openmp parrallel

    t_worker_end = sys::get_curr_time();
}

template <typename realT, class derotFunctObj, typename evCalcT, class verboseT>
int KLIPreduction<realT, derotFunctObj, evCalcT, verboseT>::finalProcess()
{
    if( this->m_postMedSub )
    {
        std::cerr << "Subtracting medians in post\n";

        for( size_t n = 0; n < this->m_psfsub.size(); ++n )
        {
            // clang-format off
            #pragma omp parallel
            // clang-format on
            {
                eigenImage<realT> medim;

                this->m_psfsub[n].median( medim );

                // clang-format off
                #pragma omp for
                // clang-format on
                for( int i = 0; i < this->m_psfsub[n].planes(); ++i )
                {
                    this->m_psfsub[n].image( i ) -= medim;
                }
            }
        }
    }

    if( this->m_doDerotate )
    {
        std::cerr << "derotating\n";
        this->derotate();
    }

    if( this->m_combineMethod != HCI:: combine::none )
    {
        std::cerr << "combining\n";
        this->combineFinim();
    }

    if( this->m_doWriteFinim == true || this->m_doOutputPSFSub == true )
    {
        std::cerr << "writing\n";

        fits::fitsHeader<verboseT> head;

        this->ADIobservation<realT, derotFunctObj, verboseT>::stdFitsHeader( &head );

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
            head.template append<realT>( "RIGHT REASON RADIUS",
                                         m_rightReasonRadius,
                                         "radius of the right reason mask" );
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

        if( this->m_doWriteFinim == true && this->m_combineMethod != HCI:: combine::none )
        {
            this->writeFinim( &head );
        }

        if( this->m_doOutputPSFSub )
        {
            this->outputPSFSub( &head );
        }
    }

    return 0;
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
