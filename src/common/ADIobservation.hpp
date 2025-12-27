/** \file ADIobservation.hpp
 * \author Jared R. Males
 * \brief Defines the ADI high contrast imaging data type.
 * \ingroup hc_imaging_files
 * \ingroup image_processing_files
 *
 */

#ifndef __ADIobservation_hpp__
#define __ADIobservation_hpp__
#include <mx/ioutils/fileUtils.hpp>

#include "HCIobservation.hpp"
#include <mx/ioutils/fits/fitsHeader.hpp>

#include <mx/improc/imagePads.hpp>

namespace mx
{
namespace improc
{

/// Process an angular differential imaging (ADI) observation
/** Angular differential imaging (ADI) uses sky rotation to differentiate real objects from
 * speckles.
 *
 * \tparam realT is the floating point type in which to do calculations
 *
 * \tparam _derotFunctObj
 * \parblock
 * is the derotation object with the following minimum interface:
 * \code
 * template<typename _realT>
 * struct derotF
 * {
 *    typedef  _realT realT;
 *
 *    //Vector of keywords to extract from the fits headers
 *    std::vector<std::string> keywords;
 *
 *    //Vector(s) to hold the keyword values
 *    std::vector<realT> keyValue1;
 *
 *    ///To allow ADIobservation to check for errors.
 *    bool isSetup()
 *    {
 *      if( <any condition indicating not set up>) return false;
 *      return true;
 *    }
 *
 *    //Method called by HCIobservation to get keyword-values
 *    void extractKeywords(vector<fitsHeader> & heads)
 *    {
 *       keyValue1 = headersToValues<float>(heads, "KEYWORD1");
 *    }
 *
 *    //Calculate the derotation angle for a given image number
 *    realT derotAngle(size_t imno) const
 *    {
 *       return some_function_of(keyValue1[imno]); //This function uses keyValue1[imno] to produce the derotation angle
 * in radians.
 *    }
 * };
 * \endcode
 * \endparblock
 *
 * \ingroup hc_imaging
 */
template <typename _realT, class _derotFunctObj, class verboseT>
struct ADIobservation : public HCIobservation<_realT, verboseT>
{
    typedef _realT realT;

    typedef _derotFunctObj derotFunctObj;

    typedef HCIobservation<_realT, verboseT>::imageT imageT;

    typedef HCIobservation<_realT, verboseT>::fitsFileT fitsFileT;
    typedef HCIobservation<_realT, verboseT>::fitsHeaderT fitsHeaderT;
    typedef HCIobservation<_realT, verboseT>::fitsHeaderCardT fitsHeaderCardT;

    derotFunctObj m_derotF;

    derotFunctObj m_RDIderotF;

    bool m_doDerotate{ true };

    bool m_postMedSub{ false };

    ADIobservation();

    void setupConfig( mx::app::appConfigurator &config );

    void loadConfig( mx::app::appConfigurator &config );

    /// Read in the target files
    /** First sets up the keywords, then calls HCIobservation readFiles
     */
    void readFiles();

    /// Actions to take after the files are first read in by HCIobservation
    /** Exracts the angle keywords from the FITS headers, and checks that a valid value is found for each file.
     * Performs fake planet injectio if configured.
     *
     * \throws
     */
    virtual void postReadFiles();

    /// Post target coadd actions.
    /** Here updates derotation for new average values.
     */
    virtual void postCoadd();

    /// Read in the RDI files
    /** First sets up the keywords, then calls HCIobservation readRDIFiles
     */
    void readRDIFiles();

    /// Post reference read actions, including fake injection
    virtual void postRDIReadFiles();

    /// Post reference coadd actions.
    /** Here updates derotation for new average values.
     */
    virtual void postRDICoadd();

    /// Read in already PSF-subtracted files
    /** Used to take up final processing after applying some non-klipReduce processing steps to
     * PSF-subtracted images.
     */
    // int readPSFSub( const std::string &dir, const std::string &prefix, const std::string &ext, size_t nReductions );

    /** \name Fake Planets
     * @{
     */
    HCI::fake m_fakeMethod{ HCI::fake::single }; ///< Method for reading fake files, either HCI::single or HCI::list.

    std::string m_fakeFileName;        ///< FITS file containing the fake planet PSF to inject or a list of fake images

    std::string m_fakeScaleFileName;   ///< One-column text file containing a scale factor for each point in time.

    std::vector<realT> m_fakeSep;      ///< Separation(s) of the fake planet(s)
    std::vector<realT> m_fakePA;       ///< Position angles(s) of the fake planet(s)
    std::vector<realT> m_fakeContrast; ///< Contrast(s) of the fake planet(s)

    realT m_fakeRDIFluxScale{ 1 };     /**< Flux scaling to apply to fake planets injected in RDI.
                                            Depend on the assumed spectrum in SDI.*/
    realT m_fakeRDISepScale{ 1 };      /**< Scaling to apply to fake planet separation in RDI.
                                            This should be the ratio of wavelengths for SDI.*/

    /// Inject the fake plants
    /**
     *  \todo should pad the fake before calling the single image version
     *  \todo throw exceptions for all errors, and switch to void
     */
    void injectFake( eigenCube<realT> &ims,              ///< [in.out] the image cube in which to inject the fakes.
                     std::vector<std::string> &fileList, /**< [in] a list of file paths used for per-image fake PSFs. If
                                                                   empty, then m_fakeFileName is used.*/
                     derotFunctObj &derotF,              ///< [in] the derotation object
                     realT RDIfluxScale,                 /**< [in] the flux scaling for RDI.  In SDI,
                                                                   this is from the planet spectrum.*/
                     realT RDISepScale                   /**< [in] the separation scale for RDI.
                                                                   In SDI, this is the ratio of wavelengths
                                                                   after lambda/D scaling. */
    );

    /// Inject the fake planet into a single image
    /**
     *  \todo should pad the fake before this point
     *  \todo throw exceptions for all errors, and switch to void
     */
    void injectFake( imageT &fakePSF,
                     eigenCube<realT> &ims,
                     int image_i,
                     realT derotAngle,
                     realT PA,
                     realT sep,
                     realT contrast,
                     realT scale,
                     realT RDIfluxScale,
                     realT RDISepScale );

    /// @}

    void stdFitsHeader( fitsHeaderT *head );

    virtual void makeMaskCube();

    /// De-rotate the PSF subtracted images
    void derotate();

    double t_fake_begin{ 0 };
    double t_fake_end{ 0 };

    double t_derotate_begin{ 0 };
    double t_derotate_end{ 0 };
};

template <typename realT, class derotFunctObj, class verboseT>
ADIobservation<realT, derotFunctObj, verboseT>::ADIobservation()
{
}

template <typename realT, class derotFunctObj, class verboseT>
void ADIobservation<realT, derotFunctObj, verboseT>::setupConfig( mx::app::appConfigurator &config )
{

    m_derotF.setupConfig( config );

    config.add( "adi.postMedSub",
                "",
                "adi.postMedSub",
                mx::app::argType::True,
                "adi",
                "postMedSub",
                false,
                "string",
                "If true, the median image is subtracted after post-processing, before de-rotation" );

    config.add( "fake.method",
                "",
                "fake.method",
                mx::app::argType::Required,
                "fake",
                "method",
                false,
                "string",
                "How the fake PSF is specified by fileName: single, if a single PSF is used (default); or list, if "
                "1 PSF per miage is used." );

    config.add( "fake.fileName",
                "",
                "fake.fileName",
                mx::app::argType::Required,
                "fake",
                "fileName",
                false,
                "string",
                "Full path to FITS file containing the fake planet PSF to inject, or a file with a list of FITS "
                "file paths." );

    config.add( "fake.scaleFileName",
                "",
                "fake.scaleFileName",
                mx::app::argType::Required,
                "fake",
                "scaleFileName",
                false,
                "string",
                "Path to one-column text file containing a scale factor for each point in time." );

    config.add( "fake.sep",
                "",
                "fake.sep",
                mx::app::argType::Required,
                "fake",
                "sep",
                false,
                "vector<float>",
                "Separation(s) of the fake planet(s) in pixels." );

    config.add( "fake.PA",
                "",
                "fake.PA",
                mx::app::argType::Required,
                "fake",
                "PA",
                false,
                "vector<float>",
                "Position angles(s) of the fake planet(s)" );

    config.add( "fake.contrast",
                "",
                "fake.contrast",
                mx::app::argType::Required,
                "fake",
                "contrast",
                false,
                "vector<float>",
                "Contrast(s) of the fake planet(s)" );

    config.add( "fake.RDIFluxScale",
                "",
                "fake.RDIFluxScale",
                mx::app::argType::Required,
                "fake",
                "RDIFluxScale",
                false,
                "vector<float>",
                "Flux scaling for the planets injected into the RDI images" );

    config.add( "fake.RDISepScale",
                "",
                "fake.RDISepScale",
                mx::app::argType::Required,
                "fake",
                "RDISepScale",
                false,
                "vector<float>",
                "Separation scaling for the planets injected into the RDI images" );

    config.add( "combine.noDerotate",
                "",
                "combine.noDerotate",
                mx::app::argType::True,
                "combine",
                "noDerotate",
                false,
                "bool",
                "Do not derotate before combining." );
}

template <typename realT, class derotFunctObj, class verboseT>
void ADIobservation<realT, derotFunctObj, verboseT>::loadConfig( mx::app::appConfigurator &config )
{
    m_derotF.loadConfig( config );

    m_RDIderotF.angleKeyword( m_derotF.m_angleKeyword );
    m_RDIderotF.m_angleScale = m_derotF.m_angleScale;
    m_RDIderotF.m_angleConstant = m_derotF.m_angleConstant;

    config( m_postMedSub, "adi.postMedSub" );

    std::string fakestr;
    try
    {
        fakestr = fakeToStr<verboseT>( m_fakeMethod ); // get default
    }
    catch( ... )
    {
        std::throw_with_nested( mx::exception<verboseT>( mx::error_t::exception ) );
    }

    config( fakestr, "fake.method" );

    try
    {
        m_fakeMethod = HCI::fakeFmStr<verboseT>( fakestr );
    }
    catch( ... )
    {
        std::throw_with_nested( mx::exception<verboseT>( mx::error_t::invalidconfig, "invalid fake method" ) );
    }

    config( m_fakeFileName, "fake.fileName" );
    config( m_fakeScaleFileName, "fake.scaleFileName" );
    config( m_fakeSep, "fake.sep" );
    config( m_fakePA, "fake.PA" );
    config( m_fakeContrast, "fake.contrast" );
    config( m_fakeRDIFluxScale, "fake.RDIFluxScale" );
    config( m_fakeRDISepScale, "fake.RDISepScale" );

    bool noDer = !m_doDerotate;
    config( noDer, "combine.noDerotate" );
    m_doDerotate = !noDer;
}

template <typename realT, class derotFunctObj, class verboseT>
void ADIobservation<realT, derotFunctObj, verboseT>::readFiles()
{
    this->m_keywords.clear();

    if( !m_derotF.isSetup() )
    {
        throw mx::exception<verboseT>( mx::error_t::paramnotset, "Derotator is not configured." );
    }

    /*----- Append the ADI keywords to propagate them if needed -----*/

    for( size_t i = 0; i < m_derotF.m_keywords.size(); ++i )
    {
        this->m_keywords.push_back( m_derotF.m_keywords[i] );
    }

    try
    {
        HCIobservation<realT, verboseT>::readFiles();
    }
    catch( ... )
    {
        std::throw_with_nested( mx::exception<verboseT>( mx::error_t::exception, "from readFiles" ) );
    }
}

template <typename realT, class derotFunctObj, class verboseT>
void ADIobservation<realT, derotFunctObj, verboseT>::postReadFiles()
{
    std::vector<size_t> bad;
    mx::error_t errc = m_derotF.extractKeywords( this->m_heads, bad );

    if( !!errc )
    {
        for( size_t n = 0; n < bad.size(); ++n )
        {
            std::cerr << this->m_fileList[bad[n]] << " conversion failed for " << m_derotF.m_angleKeyword << "\n";
        }

        throw mx::exception<verboseT>( mx::error_t::invalidarg, "bad derotation angles in FITS header" );
    }

    if( m_fakeFileName != "" && !this->m_skipPreProcess )
    {
        std::cerr << "Injecting fakes in target images...\n";

        try
        {
            injectFake( this->m_tgtIms, this->m_fileList, m_derotF, 1, 1 );
        }
        catch( ... )
        {
            std::throw_with_nested( mx::exception<verboseT>( mx::error_t::exception, "injecting fake" ) );
        }
    }
}

template <typename realT, class derotFunctObj, class verboseT>
void ADIobservation<realT, derotFunctObj, verboseT>::postCoadd()
{
    std::vector<size_t> bad;
    mx::error_t errc = m_derotF.extractKeywords( this->m_heads, bad );

    if( !!errc )
    {
        for( size_t n = 0; n < bad.size(); ++n )
        {
            std::cerr << this->m_fileList[bad[n]] << " conversion failed for " << m_derotF.m_angleKeyword << "\n";
        }

        throw mx::exception<verboseT>( mx::error_t::invalidarg, "bad derotation angles in FITS header" );
    }
}

template <typename realT, class derotFunctObj, class verboseT>
void ADIobservation<realT, derotFunctObj, verboseT>::readRDIFiles()
{
    this->m_RDIkeywords.clear();

    if( !m_RDIderotF.isSetup() )
    {
        throw mx::exception<verboseT>( mx::error_t::paramnotset, "Derotator is not configured." );
    }

    /*----- Append the ADI keywords to propagate them if needed -----*/

    for( size_t i = 0; i < m_RDIderotF.m_keywords.size(); ++i )
    {
        this->m_RDIkeywords.push_back( m_RDIderotF.m_keywords[i] );
    }

    try
    {
        HCIobservation<realT, verboseT>::readRDIFiles();
    }
    catch( ... )
    {
        std::throw_with_nested( mx::exception<verboseT>( mx::error_t::exception, "from readRDIFiles" ) );
    }
}

template <typename realT, class derotFunctObj, class verboseT>
void ADIobservation<realT, derotFunctObj, verboseT>::postRDIReadFiles()
{
    std::vector<size_t> bad;
    mx::error_t errc = m_RDIderotF.extractKeywords( this->m_RDIheads, bad );

    if( !!errc )
    {
        for( size_t n = 0; n < bad.size(); ++n )
        {
            std::cerr << this->m_RDIfileList[bad[n]] << " conversion failed for " << m_RDIderotF.m_angleKeyword << "\n";
        }

        throw mx::exception<verboseT>( mx::error_t::invalidarg, "bad derotation angles in FITS header" );
    }
}

template <typename realT, class derotFunctObj, class verboseT>
void ADIobservation<realT, derotFunctObj, verboseT>::postRDICoadd()
{
    std::vector<size_t> bad;
    mx::error_t errc = m_RDIderotF.extractKeywords( this->m_RDIheads, bad );

    if( !!errc )
    {
        for( size_t n = 0; n < bad.size(); ++n )
        {
            std::cerr << this->m_RDIfileList[bad[n]] << " conversion failed for " << m_RDIderotF.m_angleKeyword << "\n";
        }

        throw mx::exception<verboseT>( mx::error_t::invalidarg, "bad derotation angles in FITS header" );
    }
}

/*
template <typename realT, class derotFunctObj, class verboseT>
int ADIobservation<realT, derotFunctObj, verboseT>::readPSFSub( const std::string &dir,
                                                        const std::string &prefix,
                                                        const std::string &ext,
                                                        size_t nReductions )
{
    // Load first file to condigure based on its header.
    std::vector<std::string> flist = ioutils::getFileNames( dir, prefix, "000", ext );

    fits::fitsHeader fh;
    eigenImage<realT> im;
    fits::fitsFile<realT> ff;

    ff.read( im, fh, flist[0] );

    if( !m_derotF.isSetup() )
    {
        mxError( "ADIobservation::readFiles", MXE_PARAMNOTSET, "Derotator is not configured." );
        return -1;
    }

    if( fh.count( "POSTMEDS" ) != 0 )
    {
        m_postMedSub = fh["POSTMEDS"].Int();
        std::cerr << "postMedSub: " << m_postMedSub << "\n";
    }

    if( fh.count( "FAKEFILE" ) != 0 )
    {
        m_fakeFileName = fh["FAKEFILE"].String();
        std::cerr << "fakeFileName: " << m_fakeFileName << "\n";
    }

    if( fh.count( "FAKESCFL" ) != 0 )
    {
        m_fakeScaleFileName = fh["FAKESCFL"].String();
        std::cerr << "fakeScaleFileName: " << m_fakeScaleFileName << "\n";
    }

    if( fh.count( "FAKESEP" ) != 0 )
    {
        ioutils::parseStringVector( m_fakeSep, fh["FAKESEP"].String(), "," );

        if( m_fakeSep.size() == 0 )
        {
            mxError( "KLIPReduction", MXE_PARSEERR, "FAKESEP vector did not parse correctly." );
            return -1;
        }
        std::cerr << "fakeSep: " << fh["FAKESEP"].String() << "\n";
    }

    if( fh.count( "FAKEPA" ) != 0 )
    {
        ioutils::parseStringVector( m_fakePA, fh["FAKEPA"].String(), "," );

        if( m_fakePA.size() == 0 )
        {
            mxError( "KLIPReduction", MXE_PARSEERR, "FAKEPA vector did not parse correctly." );
            return -1;
        }
        std::cerr << "fakePA: " << fh["FAKEPA"].String() << "\n";
    }

    if( fh.count( "FAKECONT" ) != 0 )
    {
        ioutils::parseStringVector( m_fakeContrast, fh["FAKECONT"].String(), "," );

        if( m_fakeContrast.size() == 0 )
        {
            mxError( "KLIPReduction", MXE_PARSEERR, "FAKECONT vector did not parse correctly." );
            return -1;
        }
        std::cerr << "fakeContrast: " << fh["FAKECONT"].String() << "\n";
    }

    /*----- Append the ADI keywords to propagate them if needed -----

    for( size_t i = 0; i < m_derotF.m_keywords.size(); ++i )
    {
        this->m_keywords.push_back( m_derotF.m_keywords[i] );
    }

    if( HCIobservation<realT,verboseT>::readPSFSub( dir, prefix, ext, nReductions ) < 0 )
        return -1;

    m_derotF.extractKeywords( this->m_heads );

    return 0;
}*/

template <typename realT, class derotFunctObj, class verboseT>
void ADIobservation<realT, derotFunctObj, verboseT>::injectFake( eigenCube<realT> &ims,
                                                                 std::vector<std::string> &fileList,
                                                                 derotFunctObj &derotF,
                                                                 realT RDIFluxScale,
                                                                 realT RDISepScale )
{
    t_fake_begin = sys::get_curr_time();

    // typedef Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> imT;
    imageT fakePSF;
    std::vector<std::string> fakeFiles; // used if m_fakeMethod == HCI::list

    fits::fitsFile<realT> ff;
    std::ifstream scaleFin; // for reading the scale file.

    // Fake Scale -- default to 1, read from file otherwise
    std::vector<realT> fakeScale( ims.planes(), 1.0 );
    if( m_fakeScaleFileName != "" )
    {
        std::vector<std::string> sfileNames;
        std::vector<realT> imS;

        // Read the scale file and load it into a map
        mx::error_t errc = ioutils::readColumns( m_fakeScaleFileName, sfileNames, imS );

        if( errc != mx::error_t::noerror )
        {
            throw mx::exception<verboseT>( errc, "reading fake scale file" );
        }

        if( sfileNames.size() != imS.size() )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "fake scale file must be two columns of: fileName scale" );
        }

        std::map<std::string, realT> scales;
        for( size_t i = 0; i < sfileNames.size(); ++i )
        {
            scales[ioutils::pathFilename( sfileNames[i].c_str() )] = imS[i];
        }

        for( size_t i = 0; i < fileList.size(); ++i )
        {
            if( scales.count( ioutils::pathFilename( fileList[i].c_str() ) ) > 0 )
            {
                fakeScale[i] = scales[ioutils::pathFilename( fileList[i].c_str() )];
            }
            else
            {
                throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                               std::format( "File name not found in fakeScaleFile: {}",
                                                            ioutils::pathFilename( fileList[i].c_str() ) ) );
            }
        }
    } // if(fakeScaleFileName != "")

    if( m_fakeMethod == HCI::fake::single )
    {
        mx::error_t errc = ff.read( fakePSF, m_fakeFileName );
        if( errc != error_t::noerror )
        {
            throw mx::exception<verboseT>( errc, "reading fake PSF" );
        }
    }
    else if( m_fakeMethod == HCI::fake::list )
    {
        mx::error_t errc = ioutils::readColumns( m_fakeFileName, fakeFiles );
        if( errc != mx::error_t::noerror )
        {
            throw mx::exception<verboseT>( errc, "reading fake PSF filenames" );
        }
    }

    for( int i = 0; i < ims.planes(); ++i )
    {
        if( m_fakeMethod == HCI::fake::list )
        {
            mx::error_t errc = ff.read( fakePSF, fakeFiles[i] );
            if( errc != error_t::noerror )
            {
                throw mx::exception<verboseT>( errc, std::format( "reading fake PSF {}", i ) );
            }
        }

        for( size_t j = 0; j < m_fakeSep.size(); ++j )
        {
            try
            {
                injectFake( fakePSF,
                            ims,
                            i,
                            derotF.derotAngle( i ),
                            m_fakePA[j],
                            m_fakeSep[j],
                            m_fakeContrast[j],
                            fakeScale[j],
                            RDIFluxScale,
                            RDISepScale );
            }
            catch( ... )
            {
                std::throw_with_nested( mx::exception<verboseT>( mx::error_t::exception, "from injectFake" ) );
            }
        }
    }

    t_fake_end = sys::get_curr_time();
}

template <typename realT, class derotFunctObj, class verboseT>
void ADIobservation<realT, derotFunctObj, verboseT>::injectFake( imageT &fakePSF,
                                                                 eigenCube<realT> &ims,
                                                                 int image_i,
                                                                 realT derotAngle,
                                                                 realT PA,
                                                                 realT sep,
                                                                 realT contrast,
                                                                 realT scale,
                                                                 realT RDIFluxScale,
                                                                 realT RDISepScale )
{
    // Check for correct sizing
    if( ( fakePSF.rows() < ims.rows() && fakePSF.cols() >= ims.cols() ) ||
        ( fakePSF.rows() >= ims.rows() && fakePSF.cols() < ims.cols() ) )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr,
                                       "fake PSF has different dimensions and can't be sized properly" );
    }

    // Check if fake needs to be padded out
    if( fakePSF.rows() < ims.rows() && fakePSF.cols() < ims.cols() )
    {
        imageT pfake( ims.rows(), ims.cols() );
        padImage( pfake, fakePSF, 0.5 * ( ims.rows() - fakePSF.rows() ), 0 );
        fakePSF = pfake;
    }

    // Check if fake needs to be cut down
    if( fakePSF.rows() > ims.rows() && fakePSF.cols() > ims.cols() )
    {
        imageT cfake( ims.rows(), ims.cols() );
        cutPaddedImage( cfake, fakePSF, 0.5 * ( fakePSF.rows() - ims.rows() ) );
        fakePSF = cfake;
    }

    if( fakePSF.rows() != ims.rows() || fakePSF.cols() != ims.cols() )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr,
                                       "fake PSF has different dimensions and can't be "
                                       "sized properly (is it even in rows and cols?)" );
    }

    /*** Now shift to the separation and PA, scale, apply contrast, and inject ***/
    // allocate shifted fake psf
    imageT shiftFake( fakePSF.rows(), fakePSF.cols() );

    realT ang, dx, dy;

    ang = math::dtor( -1 * PA ) + derotAngle;

    dx = sep * RDISepScale * sin( ang );
    dy = sep * RDISepScale * cos( ang );

    imageShift( shiftFake, fakePSF, dx, dy, cubicConvolTransform<realT>() );

    ims.image( image_i ) = ims.image( image_i ) + shiftFake * scale * RDIFluxScale * contrast;
}

template <typename realT, class derotFunctObj, class verboseT>
void ADIobservation<realT, derotFunctObj, verboseT>::makeMaskCube()
{
    if( this->m_mask.rows() != this->m_Nrows || this->m_mask.cols() != this->m_Ncols )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                       std::format( "Mask is not the same size as images.\n"
                                                    "    Mask:   rows={}\n"
                                                    "            cols={}\n"
                                                    "    Images: rows={}\n"
                                                    "            cols={}\n",
                                                    this->m_mask.rows(),
                                                    this->m_mask.cols(),
                                                    this->m_Nrows,
                                                    this->m_Ncols ) );
    }

    this->m_maskCube.resize( this->m_Nrows, this->m_Ncols, this->m_Nims );

    // clang-format off
    #pragma omp parallel // clang-format on
    {
        imageT rm;

        // clang-format off
        #pragma omp for // clang-format on
        for( int i = 0; i < this->m_Nims; ++i )
        {
            rotateMask( rm, this->m_mask, m_derotF.derotAngle( i ) );
            this->m_maskCube.image( i ) = rm;
        }
    }

    ioutils::createDirectories( this->m_auxDataDir );
    std::ofstream fout( this->m_auxDataDir + "angles.dat" );
    for( int i = 0; i < this->m_Nims; ++i )
    {
        fout << m_derotF.derotAngle( i ) << "\n";
    }
    fout.close();

    fits::fitsFile<realT> ff;
    ff.write( this->m_auxDataDir + "maskCube.fits", this->m_maskCube );
}

template <typename realT, class derotFunctObj, class verboseT>
void ADIobservation<realT, derotFunctObj, verboseT>::derotate()
{
    t_derotate_begin = sys::get_curr_time();

    for( size_t n = 0; n < this->m_psfsub.size(); ++n )
    {
#pragma omp parallel
        {
            imageT rotim;
            realT derot;

#pragma omp for
            for( int i = 0; i < this->m_psfsub[n].planes(); ++i )
            {
                derot = m_derotF.derotAngle( i );
                if( derot != 0 )
                {
                    imageRotate( rotim, this->m_psfsub[n].image( i ), derot, cubicConvolTransform<realT>() );
                    this->m_psfsub[n].image( i ) = rotim;
                }
            }
        }
    }

    t_derotate_end = sys::get_curr_time();
}

// If fakeFileName == "" or skipPreProcess == true then use the structure of propagated values

template <typename realT, class derotFunctObj, class verboseT>
void ADIobservation<realT, derotFunctObj, verboseT>::stdFitsHeader( fitsHeaderT *head )
{
    if( head == 0 )
        return;

    head->append( "", fits::fitsCommentType(), "----------------------------------------" );
    head->append( "", fits::fitsCommentType(), "mx::ADIobservation parameters:" );
    head->append( "", fits::fitsCommentType(), "----------------------------------------" );

    head->append( "POSTMEDS", m_postMedSub, "median subtraction after processing" );

    if( m_fakeFileName != "" )
    {
        head->append( "FAKEFILE", m_fakeFileName, "name of fake planet PSF file" );

        if( m_fakeScaleFileName != "" )
        {
            head->append( "FAKESCFL", m_fakeScaleFileName, "name of fake planet scale file name" );
        }

        std::stringstream str;

        if( m_fakeSep.size() > 0 )
        {
            for( size_t nm = 0; nm < m_fakeSep.size() - 1; ++nm )
            {
                str << m_fakeSep[nm] << ",";
            }
            str << m_fakeSep[m_fakeSep.size() - 1];
            head->append( "FAKESEP", (const char *)str.str().c_str(), "separation of fake planets" );
        }

        if( m_fakePA.size() > 0 )
        {
            str.str( "" );
            for( size_t nm = 0; nm < m_fakePA.size() - 1; ++nm )
            {
                str << m_fakePA[nm] << ",";
            }
            str << m_fakePA[m_fakePA.size() - 1];
            head->append( "FAKEPA", (char *)str.str().c_str(), "PA of fake planets" );
        }

        if( m_fakeContrast.size() > 0 )
        {
            str.str( "" );
            for( size_t nm = 0; nm < m_fakeContrast.size() - 1; ++nm )
            {
                str << m_fakeContrast[nm] << ",";
            }
            str << m_fakeContrast[m_fakeContrast.size() - 1];
            head->template append<char *>( "FAKECONT", (char *)str.str().c_str(), "Contrast of fake planets" );
        }
    }
}

template <typename realT, class verboseT>
class ADIDerotator;

extern template class ADIobservation<float, ADIDerotator<float, verbose::vv>, verbose::vv>;

} // namespace improc
} // namespace mx

#endif //__ADIobservation_hpp__
