
#include <iostream>

#include <mx/app/application.hpp>
using namespace mx::app;

#include "../common/ADIDerotator.hpp"
#include "../common/KLIPreduction.hpp"
using namespace mx::improc;

#include <libgen.h>

/// A program to run the KLIP pipeline
/**
 *
 */
template <typename _realT, typename _evCalcT = double>
class klipReduce : public application
{
  public:
    typedef _realT realT;
    typedef _evCalcT evCalcT;

  protected:

    // Rotation setup [ADIobservation]
    std::string angleKeyword;
    float angleScale;
    float angleConstant;

    // Fake Planet Injection [ADIobservation]
    std::string fakeMethod{ "single" };
    std::string fakeFileName;        ///< FITS file containing the fake planet PSF to inject
    std::string fakeScaleFileName;   ///< One-column text file containing a scale factor for each point in time.

    std::vector<realT> fakeSep;      ///< Separation(s) of the fake planet(s)
    std::vector<realT> fakePA;       ///< Position angles(s) of the fake planet(s)
    std::vector<realT> fakeContrast; ///< Contrast(s) of the fake planet(s)
    realT m_fakeRDIFluxScale{ 1 };
    realT m_fakeRDISepScale{ 1 };


    // ADI
    bool m_postMedSub{ false };

    // KLIP parameters

    std::string m_meanSubMethod;
    std::string m_pixelTSNormMethod {"none"};

    realT minDPx{ 0 };
    realT maxDPx{ 0 };
    std::string excludeMethod;
    std::string excludeMethodMax;
    int includeRefNum;
    std::vector<int> Nmodes;
    std::vector<realT> minRadius;
    std::vector<realT> maxRadius;
    std::vector<realT> minAngle;
    std::vector<realT> maxAngle;
    int nWedges{ 0 };

    bool m_rightReason{ false };
    realT m_rightReasonRadius{ 2.5 };

    // Combination and Output

    bool noDerotate;

    std::string combineMethod;
    std::string weightFile;
    _realT sigmaThreshold;
    _realT minGoodFract{ 0.0 };

    std::string outputFile;
    bool exactFName;
    std::string outputDir;

    // Output the individual PSF subtracted images
    bool outputPSFSub;
    std::string psfSubPrefix;

    //************************************//
    // Mode of execution                  //
    std::string mode; // Choices are basic [default] and grid [grid manages a fake planet grid]

    // Executes a grid of fake planets.
    realT gridCenterSep;              ///< The separation of the grid center [pixels].
    realT gridCenterPA;               ///< The PA of the grid center [deg E of N].
    realT gridHalfWidthSep;           ///< The grid half-width in radius [pixels]
    realT gridDeltaSep;               ///< The grid spacing in radius [pixels]
    realT gridHalfWidthPA;            ///< The grid half-wdith in PA [pixels]
    realT gridDeltaPA;                ///< The grid spacing in PA [pixels]
    std::vector<realT> gridContrasts; ///< The grid contrasts, possibly negative.

    int doGrid();

    KLIPreduction<realT, ADIDerotator<realT>, evCalcT> m_obs;

  public:
    klipReduce()
    {
        m_configPathGlobal_env = "KLIPREDUCE_GLOBAL_CONFIG";
        m_configPathLocal = "klipReduce.conf";
        m_requireConfigPathLocal = false;

        angleScale = 1;
        angleConstant = 0;

        includeRefNum = 0;

        noDerotate = false;

        sigmaThreshold = 0;

        exactFName = false;

        outputPSFSub = false;

        // Grid setup
        gridCenterSep = 0;
        gridCenterPA = 0;
        gridHalfWidthSep = 0;
        gridDeltaSep = 0;
        gridHalfWidthPA = 0;
        gridDeltaPA = 0;

        mode = "basic";

        config.m_sources = true;
    }

    ~klipReduce()
    {
    }

    // This sets up the configuration
    void setupConfig()
    {
        m_obs.setupConfig( config );

        config.add( "input.angleKeyword",
                    "",
                    "input.angleKeyword",
                    argType::Required,
                    "input",
                    "angleKeyword",
                    false,
                    "string",
                    "The header keyword to use for the rotation angle of each image." );
        config.add( "input.angleScale",
                    "",
                    "input.angleScale",
                    argType::Required,
                    "input",
                    "angleScale",
                    false,
                    "float",
                    "The scale to apply to the angle, default is 1." );
        config.add( "input.angleConstant",
                    "",
                    "input.angleConstant",
                    argType::Required,
                    "input",
                    "angleConstant",
                    false,
                    "float",
                    "The offset to apply to each angle (e.g. the North angle), default is 0." );


        /*>>>> preProcess */


        /*<<<< preProcess */

        /*>>>> adi */

        config.add( "adi.minDPx",
                    "",
                    "adi.minDPx",
                    argType::Required,
                    "adi",
                    "minDPx",
                    false,
                    "float",
                    "Specify the minimum angle or pixel difference at the inner edge of the search region" );
        config.add( "adi.maxDPx",
                    "",
                    "adi.maxDPx",
                    argType::Required,
                    "adi",
                    "maxDPx",
                    false,
                    "float",
                    "Specify the maximum angle or pixel difference at the inner edge of the search region" );
        config.add( "adi.excludeMethod",
                    "",
                    "adi.excludeMethod",
                    argType::Required,
                    "adi",
                    "excludeMethod",
                    false,
                    "string",
                    "Method for minimum exclusion.  Values are none (default), pixel, angle, imno." );
        config.add( "adi.excludeMethodMax",
                    "",
                    "adi.excludeMethodMax",
                    argType::Required,
                    "adi",
                    "excludeMethodMax"
    /** \name Pre-Processing
     * @{

    /// Do the pre-processing
    void preProcess( eigenCube<realT> &ims /**< [in] the image cube, should be either m_tgtIms or m_refIms  );

    /// Do mean subtraction as part of pre-processing
    void preProcess_meanSub( eigenCube<realT> &ims /**< [in] the image cube, should be either m_tgtIms or m_refIms  );

    /// Do pixel time-series normalization as part of pre-processing
    void
    preProcess_pixelTSNorm( eigenCube<realT> &ims /**< [in] the image cube, should be either m_tgtIms or m_refIms  );

    ///@}

    /** \name Image Combination
     * These options control how the final image combination is performed.
     * @{


    /// Determine how to combine the PSF subtracted images
    /** Possibilities are
     * - HCI::noCombine -- do not combine
     * - HCI::combineMethod::median -- [default] final image is the median
     * - HCI::combineMethod::mean -- final image is the simple mean
     * - HCI::weightedMeanCombine -- final image is the weighted mean.  m_weightFile must be provided.
     * - HCI::combineMethod::sigmaMean -- final image is sigma clipped mean.  If m_sigmaThreshold \<= 0, then it reverts
     * to meanCombine.

    HCI::combineMethod m_combineMethod{ HCI::combineMethod::mean };

    /// Specifies a file containing the image weights, for combining with weighted mean.
    /** This 2-column space-delimited ASCII file containing  filenames and weights. It must be specified before
     * readFiles() is executed.  In readFiles this is loaded after the m_fileList is cutdown and matched to the
     * remaining files.

    std::string m_weightFile;

    /// Vector to hold the image weights read from the m_weightFile.
    /** After readWeights is executed by readFiles, this will contain the normalized weights.
     * \todo check how comboWeights are handled in coadding

    std::vector<realT> m_comboWeights;

    /// The standard deviation threshold used if combineMethod == HCI::combineMethod::sigmaMean.
    realT m_sigmaThreshold{ 0 };

    /// The minimum fraction of good (un-masked) pixels to include in the final combination (0.0 to 1.0). If not met,
    /// then the pixel will be NaN-ed.
    realT m_minGoodFract{ 0.0 };

    /// Read the image weights from m_weightFile
    /** This is called by readFiles().
     *
     * \returns 0 on success, -1 on  error.

    int readWeights();

    /// Combine the images into a single final image.
    /** Images are combined by the method specified in \ref combineMethod

    void combineFinim();

    ///@}*/

    /** \name Output
     * These options control the ouput of the final combined images and the individual PSF subtracted images.
     * @{


    /// Location for temporary auxilliary output files (e.g. masks)
    std::string m_auxDataDir{ "/tmp/hciReduceAux/" };

    /// Whether or not to move the temp. aux files.



    ///@}*/,
                    false,
                    "string",
                    "Method for maximum exclusion.  Values are none (default), pixel, angle, imno." );
        config.add( "adi.postMedSub",
                    "",
                    "adi.postMedSub",
                    argType::True,
                    "adi",
                    "postMedSub",
                    false,
                    "string",
                    "If true, the median image is subtracted after post-processing, before de-rotation" );

        /*<<<< adi */

        /*>>>> geom */
        config.add( "geom.minRadius",
                    "",
                    "geom.minRadius",
                    argType::Required,
                    "geom",
                    "minRadius",
                    false,
                    "vector<realT>",
                    "The minimum radius of the search regions" );
        config.add( "geom.maxRadius",
                    "",
                    "geom.maxRadius",
                    argType::Required,
                    "geom",
                    "maxRadius",
                    false,
                    "vector<realT>",
                    "The maximum radius of the search regions" );
        config.add( "geom.minAngle",
                    "",
                    "geom.minAngle",
                    argType::Required,
                    "geom",
                    "minAngle",
                    false,
                    "vector<realT>",
                    "The minimum angle of the search regions" );
        config.add( "geom.maxAngle",
                    "",
                    "geom.maxAngle",
                    argType::Required,
                    "geom",
                    "maxAngle",
                    false,
                    "vector<realT>",
                    "The maximum angle of the search regions" );
        config.add(
            "geom.nWedges",
            "",
            "geom.nWedges",
            argType::Required,
            "geom",
            "nWedges",
            false,
            "",
            "The number of angular wedges.  Overrides minAngle and maxAngle, and expands minRadius and maxRadius" );

        /*<<<< geom */

        /*>>>> klip */
        config.add( "klip.meanSubMethod",
                    "",
                    "klip.meanSubMethod",
                    argType::Required,
                    "klip",
                    "meanSubMethod",
                    false,
                    "string",
                    "The method of mean subtraction for PCA: imageMean, imageMedian, meanImage, or medianImage." );
        config.add( "klip.pixelTSNormMethod",
                    "",
                    "klip.pixelTSNormMethod",
                    argType::Required,
                    "klip",
                    "pixelTSNormMethod",
                    false,
                    "int",
                    "The method of pixel time-series normalization for PCA: none or rms." );
        config.add( "klip.includeRefNum",
                    "",
                    "klip.includeRefNum",
                    argType::Required,
                    "klip",
                    "includeRefNum",
                    false,
                    "int",
                    "The number of references to include, based on correlation." );
        config.add( "klip.Nmodes",
                    "",
                    "klip.Nmodes",
                    argType::Required,
                    "klip",
                    "Nmodes",
                    false,
                    "vector<int>",
                    "The number of modes to included in the PSF estimate." );
        config.add( "klip.rightReason",
                    "",
                    "klip.rightReason",
                    argType::Required,
                    "klip",
                    "rightReason",
                    false,
                    "bool",
                    "Whether or not the right reason mask is applied" );
        config.add( "klip.rrRadius",
                    "",
                    "klip.rrRadius",
                    argType::Required,
                    "klip",
                    "rrRadius",
                    false,
                    "float",
                    "The radius of the right reason mask" );

        /*<<<< klip */

        /*>>>> combine */
        config.add( "combine.noDerotate",
                    "",
                    "combine.noDerotate",
                    argType::True,
                    "combine",
                    "noDerotate",
                    false,
                    "bool",
                    "Do not derotate before combining." );

        config.add( "combine.method",
                    "",
                    "combine.method",
                    argType::Required,
                    "combine",
                    "method",
                    false,
                    "string",
                    "Averaging method for final combination: mean, median, weighted, sigma" );
        config.add( "combine.weightFile",
                    "",
                    "combine.weightFile",
                    argType::Required,
                    "combine",
                    "weightFile",
                    false,
                    "string",
                    "File containing weights for the weighted combo.  Two column format: filename weight" );
        config.add( "combine.sigmaThreshold",
                    "",
                    "combine.sigmaThreshold",
                    argType::Required,
                    "combine",
                    "sigmaThreshold",
                    false,
                    "float",
                    "Clipping threshold for sigma clipped mean combination." );
        config.add( "combine.minGoodFract",
                    "",
                    "combine.minGoodFract",
                    argType::Required,
                    "combine",
                    "minGoodFract",
                    false,
                    "float",
                    "Minimum fraction of good/un-masked pixels to include in final image, otherwise pixel is NaN-ed." );

        /*<<<< combine */

        /*>>>> output */
        config.add( "output.fileName",
                    "",
                    "output.fileName",
                    argType::Required,
                    "output",
                    "fileName",
                    false,
                    "string",
                    "Prefix for output file name.  A 4 digit 0-padded number is appended." );
        config.add( "output.exactFName",
                    "",
                    "output.exactFName",
                    argType::True,
                    "output",
                    "exactFName",
                    false,
                    "bool",
                    "Used outputFile exactly as specified, without appending a number or .fits" );
        config.add( "output.directory",
                    "",
                    "output.directory",
                    argType::Required,
                    "output",
                    "directory",
                    false,
                    "string",
                    "The directory where to output files." );
        config.add( "output.outputPSFSub",
                    "",
                    "output.outputPSFSub",
                    argType::True,
                    "output",
                    "outputPSFSub",
                    false,
                    "bool",
                    "Output the PSF subtracted images (default false)" );
        config.add( "output.psfSubPrefix",
                    "",
                    "output.psfSubPrefix",
                    argType::Required,
                    "output",
                    "psfSubPrefix",
                    false,
                    "string",
                    "Prefix of the PSF subtracted output files." );
        /*<<<< output */

        /*>>>> fake */
        config.add( "fake.method",
                    "",
                    "fake.method",
                    argType::Required,
                    "fake",
                    "method",
                    false,
                    "string",
                    "How the fake PSF is specified by fileName: single, if a single PSF is used (default); or list, if "
                    "1 PSF per miage is used." );
        config.add( "fake.fileName",
                    "",
                    "fake.fileName",
                    argType::Required,
                    "fake",
                    "fileName",
                    false,
                    "string",
                    "Full path to FITS file containing the fake planet PSF to inject, or a file with a list of FITS "
                    "file paths." );
        config.add( "fake.scaleFileName",
                    "",
                    "fake.scaleFileName",
                    argType::Required,
                    "fake",
                    "scaleFileName",
                    false,
                    "string",
                    "Path to one-column text file containing a scale factor for each point in time." );
        config.add( "fake.sep",
                    "",
                    "fake.sep",
                    argType::Required,
                    "fake",
                    "sep",
                    false,
                    "vector<float>",
                    "Separation(s) of the fake planet(s) in pixels." );
        config.add( "fake.PA",
                    "",
                    "fake.PA",
                    argType::Required,
                    "fake",
                    "PA",
                    false,
                    "vector<float>",
                    "Position angles(s) of the fake planet(s)" );
        config.add( "fake.contrast",
                    "",
                    "fake.contrast",
                    argType::Required,
                    "fake",
                    "contrast",
                    false,
                    "vector<float>",
                    "Contrast(s) of the fake planet(s)" );
        config.add( "fake.RDIFluxScale",
                    "",
                    "fake.RDIFluxScale",
                    argType::Required,
                    "fake",
                    "RDIFluxScale",
                    false,
                    "vector<float>",
                    "Flux scaling for the planets injected into the RDI images" );
        config.add( "fake.RDISepScale",
                    "",
                    "fake.RDISepScale",
                    argType::Required,
                    "fake",
                    "RDISepScale",
                    false,
                    "vector<float>",
                    "Separation scaling for the planets injected into the RDI images" );

        /*<<<< fake */

        config.add( "mode",
                    "",
                    "mode",
                    argType::Required,
                    "",
                    "mode",
                    false,
                    "string",
                    "The mode of operation: either \"grid\" or \"normal\" (the default)" );
        config.add( "grid.centerSep",
                    "",
                    "grid.centerSep",
                    argType::Required,
                    "grid",
                    "centerSep",
                    false,
                    "float",
                    "The grid center in separation [pixels]" );
        config.add( "grid.centerPA",
                    "",
                    "grid.centerPA",
                    argType::Required,
                    "grid",
                    "centerPA",
                    false,
                    "float",
                    "The grid center in position angle [degrees]" );
        config.add( "grid.halfWidthSep",
                    "",
                    "grid.halfWidthSep",
                    argType::Required,
                    "grid",
                    "halfWidthSep",
                    false,
                    "float",
                    "The half width of the grid in spearation [pixels]" );
        config.add( "grid.halfWidthPA",
                    "",
                    "grid.halfWidthPA",
                    argType::Required,
                    "grid",
                    "halfWidthPA",
                    false,
                    "float",
                    "The half width of the grid in PA [degrees]" );
        config.add( "grid.deltaSep",
                    "",
                    "grid.deltaSep",
                    argType::Required,
                    "grid",
                    "deltaSep",
                    false,
                    "float",
                    "The grid step size in separation [pixels]" );
        config.add( "grid.deltaPA",
                    "",
                    "grid.deltaPA",
                    argType::Required,
                    "grid",
                    "deltaPA",
                    false,
                    "float",
                    "The grid step size in PA [degrees]" );
        config.add( "grid.contrasts",
                    "",
                    "grid.contrasts",
                    argType::Required,
                    "grid",
                    "contrasts",
                    false,
                    "vector<float>",
                    "The contrast grid [planet:star]." );
        // config.add("","", "",argType::Required, "", ""));
    }

    virtual void setConfigPathCL()
    {
        config.get<std::string>( m_configPathCL, "config" );
    }

    void loadConfig()
    {
        m_obs.loadConfig( config );

        config( angleKeyword, "input.angleKeyword" );
        config( angleScale, "input.angleScale" );
        config( angleConstant, "input.angleConstant" );


        /*>>>> adi */
        config( minDPx, "adi.minDPx" );
        config( maxDPx, "adi.maxDPx" );
        config( excludeMethod, "adi.excludeMethod" );
        config( excludeMethodMax, "adi.excludeMethodMax" );
        config( m_postMedSub, "adi.postMedSub" );

        /*<<<< adi */

        /*>>>> geom */
        config( minRadius, "geom.minRadius" );
        config( maxRadius, "geom.maxRadius" );
        config( minAngle, "geom.minAngle" );
        config( maxAngle, "geom.maxAngle" );
        config( nWedges, "geom.nWedges" );
        /*<<<< geom */

        /*>>>> klip */
        config( m_meanSubMethod, "klip.meanSubMethod" );
        config( m_pixelTSNormMethod, "klip.pixelTSNormMethod" );

        config( includeRefNum, "klip.includeRefNum" );
        config( Nmodes, "klip.Nmodes" );
        config( m_rightReason, "klip.rightReason" );
        config( m_rightReasonRadius, "klip.rrRadius" );

        /*<<<< klip */

        /*>>>> combine */

        config( noDerotate, "combine.noDerotate" );

        config( combineMethod, "combine.method" );

        config( weightFile, "combine.weightFile" );

        config( sigmaThreshold, "combine.sigmaThreshold" );
        config( minGoodFract, "combine.minGoodFract" );

        /*<<<< combine */

        /*>>>> output */
        config( outputFile, "output.fileName" );
        config( exactFName, "output.exactFName" );
        config( outputDir, "output.directory" );
        config( outputPSFSub, "output.outputPSFSub" );
        config( psfSubPrefix, "output.psfSubPrefix" );
        /*<<<< output */

        /*>>>> fake */
        config( fakeMethod, "fake.method" );
        config( fakeFileName, "fake.fileName" );
        config( fakeScaleFileName, "fake.scaleFileName" );
        config( fakeSep, "fake.sep" );
        config( fakePA, "fake.PA" );
        config( fakeContrast, "fake.contrast" );
        config( m_fakeRDIFluxScale, "fake.RDIFluxScale" );
        config( m_fakeRDISepScale, "fake.RDISepScale" );
        /*<<<< fake */

        config( gridCenterSep, "grid.centerSep" );

        config( gridCenterPA, "grid.centerPA" );
        config( gridHalfWidthSep, "grid.halfWidthSep" );
        config( gridDeltaSep, "grid.deltaSep" );
        config( gridHalfWidthPA, "grid.halfWidthPA" );
        config( gridDeltaPA, "grid.deltaPA" );
        config( gridContrasts, "grid.contrasts" );

        config( mode, "mode" );

        // This checks for unused config options, printing the banner only once no matter how many there are.
        // This will catch both bad options, and options we aren't actually using (debugging).
        bool unusedPrinted = false;
        for( auto it = config.m_targets.begin(); it != config.m_targets.end(); ++it )
        {
            if( it->second.used == false )
            {
                if( !unusedPrinted )
                {
                    std::cerr << "****************************************************\n";
                    std::cerr << "WARNING: unused config options (this is a programmer error):\n";
                    unusedPrinted = true;
                }

                std::cerr << "   " << it->second.name << '\n';
            }
        }

        unusedPrinted = false;
        if( config.m_unusedConfigs.size() > 0 )
        {
            if( !unusedPrinted )
            {
                std::cerr << "****************************************************\n";
                std::cerr << "WARNING: unrecognized config options:\n";
                unusedPrinted = true;
            }

            for( auto it = config.m_unusedConfigs.begin(); it != config.m_unusedConfigs.end(); ++it )
            {
                std::cerr << "   " << it->second.name;
                if( config.m_sources )
                    std::cerr << " [" << it->second.sources[0] << "]\n";
                else
                    std::cerr << "\n";
            }

            std::cerr << "****************************************************\n";
        }

        if( config.nonOptions.size() > 0 )
        {
            std::cerr << "****************************************************\n";
            std::cerr << "WARNING: unrecognized command line arguments\n";
        }
    }

    void printUsage()
    {
        fprintf( stderr, "%s: Perform a KLIP reduction using the mxlib pipeline.\n\n", invokedName.c_str() );
        fprintf( stderr,
                 "   usage: %s -D directory -P prefix [-E extension] -n \"x,y,z\" -r X -R x \n\n",
                 invokedName.c_str() );
        fprintf( stderr, "   For usage and full documentation see somewhere \n\n" );
    }

    bool m_configError{ false };

    void checkConfig()
    {

        if( doHelp )
        {
            m_configError = true;
            return;
        }

        m_obs.m_derotF.angleKeyword( angleKeyword );
        m_obs.m_derotF.m_angleScale = angleScale;
        m_obs.m_derotF.m_angleConstant = angleConstant;

        m_obs.m_RDIderotF.angleKeyword( angleKeyword );
        m_obs.m_RDIderotF.m_angleScale = angleScale;
        m_obs.m_RDIderotF.m_angleConstant = angleConstant;

        if( fakeMethod == "list" )
        {
            m_obs.m_fakeMethod = HCI::list;
        }
        else
        {
            m_obs.m_fakeMethod = HCI::single;
        }

        m_obs.m_fakeFileName = fakeFileName;
        m_obs.m_fakeScaleFileName = fakeScaleFileName;
        m_obs.m_fakeSep = fakeSep;
        m_obs.m_fakePA = fakePA;
        m_obs.m_fakeContrast = fakeContrast;
        m_obs.m_RDIFluxScale = m_fakeRDIFluxScale;
        m_obs.m_RDISepScale = m_fakeRDISepScale;



        m_obs.m_minDPx = minDPx;
        m_obs.m_maxDPx = maxDPx;

        if( excludeMethod != "" )
        {
            if( excludeMethod == "none" )
            {
                m_obs.m_excludeMethod = HCI::excludeNone;
            }
            else if( excludeMethod == "pixel" )
            {
                m_obs.m_excludeMethod = HCI::excludePixel;
            }
            else if( excludeMethod == "angle" )
            {
                m_obs.m_excludeMethod = HCI::excludeAngle;
            }
            else if( excludeMethod == "imno" )
            {
                m_obs.m_excludeMethod = HCI::excludeImno;
            }
            else
            {
                std::cerr << invokedName << ": invalid excludeMethod.\n";
                m_configError = true;
            }
        }

        if( excludeMethodMax != "" )
        {
            if( excludeMethodMax == "none" )
            {
                m_obs.m_excludeMethodMax = HCI::excludeNone;
            }
            else if( excludeMethodMax == "pixel" )
            {
                m_obs.m_excludeMethodMax = HCI::excludePixel;
            }
            else if( excludeMethodMax == "angle" )
            {
                m_obs.m_excludeMethodMax = HCI::excludeAngle;
            }
            else if( excludeMethodMax == "imno" )
            {
                m_obs.m_excludeMethodMax = HCI::excludeImno;
            }
            else
            {
                std::cerr << invokedName << ": invalid excludeMethodMax.\n";
                m_configError = true;
            }
        }

        m_obs.m_postMedSub = m_postMedSub;

        //KLIP:
        m_obs.m_meanSubMethod = HCI::meanSubMethodStr( m_meanSubMethod );

        //KLIP:
        m_obs.m_pixelTSNormMethod == HCI::pixelTSNormMethodStr( m_pixelTSNormMethod );

        m_obs.m_includeRefNum = includeRefNum;

        if( Nmodes.size() == 0 && mode != "postprocess" )
        {
            std::cerr << invokedName << ": must specify number of modes (Nmodes)\n";
            m_configError = true;
        }
        else
            m_obs.m_Nmodes = Nmodes;

        m_obs.m_rightReason = m_rightReason;
        m_obs.m_rightReasonRadius = m_rightReasonRadius;

        if( minRadius.size() == 0 && mode != "postprocess" )
        {
            std::cerr << invokedName << ": must specify minimum radii of KLIP regions (minRadius)\n";
            m_configError = true;
        }

        if( maxRadius.size() == 0 && mode != "postprocess" )
        {
            std::cerr << invokedName << ": must specify maximum radii of KLIP regions (maxRadius)\n";
            m_configError = true;
        }

        if( minRadius.size() != maxRadius.size() && mode != "postprocess" )
        {
            std::cerr << invokedName << ": number of minimum and maximum radii must be equal\n";
            m_configError = true;
        }

        if( minAngle.size() != maxAngle.size() && mode != "postprocess" )
        {
            std::cerr << invokedName << ": number of minimum and maximum angles must be equal\n";
            m_configError = true;
        }

        m_obs.m_doDerotate = !noDerotate;

        if( combineMethod != "" )
        {
            if( combineMethod == "none" )
            {
                m_obs.m_combineMethod = HCI::combineMethod::none;
            }
            else if( combineMethod == "median" )
            {
                m_obs.m_combineMethod = HCI::combineMethod::median;
            }
            else if( combineMethod == "mean" )
            {
                m_obs.m_combineMethod = HCI::combineMethod::mean;
            }
            else if( combineMethod == "sigma" )
            {
                m_obs.m_combineMethod = HCI::combineMethod::sigmaMean;
            }
            else
            {
                std::cerr << invokedName << ": invalid combine method.\n";
                m_configError = true;
            }
        }

        m_obs.m_weightFile = weightFile;
        m_obs.m_sigmaThreshold = sigmaThreshold;
        m_obs.m_minGoodFract = minGoodFract;

        if( outputFile != "" )
        {
            m_obs.m_finimName = outputFile;
        }

        m_obs.m_exactFinimName = exactFName;

        m_obs.m_outputDir = outputDir;

        m_obs.m_doOutputPSFSub = outputPSFSub;
        if( psfSubPrefix != "" )
        {
            m_obs.m_PSFSubPrefix = psfSubPrefix;
        }

        return;
    }

    virtual int execute()
    {
        if( m_configError )
        {
            printUsage();
            return -1;
        }

        if( mode == "grid" )
        {
            return doGrid();
        }
        /*else if( mode == "postprocess" )
        {
            return m_obs.processPSFSub( directory, prefix, extension );
        }*/
        else
        {
            if( nWedges > 0 )
            {
                if( minRadius.size() > 1 || maxRadius.size() > 1 )
                {
                    std::cerr << "Error: nWedges set but min/maxRadius have more than one entry" << "\n";
                    return -1;
                }

                if( minAngle.size() > 0 || maxAngle.size() > 0 )
                {
                    std::cerr << "Warning: nWedges set but min/maxAngle have more than one entry" << "\n";
                    minAngle.clear();
                    maxAngle.clear();
                }

                if( ( 360 % nWedges ) != 0 )
                {
                    std::cerr << "Error: nWedges must be a divisor of 360\n";
                    return -1;
                }

                realT mnr = minRadius[0];
                realT mxr = maxRadius[0];

                minRadius.resize( nWedges, mnr );
                maxRadius.resize( nWedges, mxr );

                int dang = 360 / nWedges;
                minAngle.resize( nWedges );
                maxAngle.resize( nWedges );

                for( size_t n = 0; n < minAngle.size(); ++n )
                {
                    minAngle[n] = n * dang;
                    maxAngle[n] = n * dang + dang;
                }
            }
            else
            {
                if( minAngle.size() == 0 || maxAngle.size() == 0 )
                {
                    minAngle.resize( minRadius.size(), 0 );
                    maxAngle.resize( maxRadius.size(), 360 );
                }
            }

            m_obs.load_fileList();
            return m_obs.regions( minRadius, maxRadius, minAngle, maxAngle );
        }
    }
};

template <typename realT, typename evCalcT>
int klipReduce<realT, evCalcT>::doGrid()
{
    if( gridCenterSep == 0 )
    {
        mxError( "klipReduce", MXE_PARAMNOTSET, "Grid center separation not set (grid.centerSep)" );
        return -1;
    }

    if( gridCenterPA == 0 )
    {
        mxError( "klipReduce", MXE_PARAMNOTSET, "Grid center PA not set (grid.centerPA)" );
        return -1;
    }

    if( gridHalfWidthSep == 0 )
    {
        mxError( "klipReduce", MXE_PARAMNOTSET, "Grid half-width in radius not set (grid.halfWidthSep)" );
        return -1;
    }

    if( gridHalfWidthPA == 0 )
    {
        mxError( "klipReduce", MXE_PARAMNOTSET, "Grid half-width in PA not set (grid.halfWidthPA)" );
        return -1;
    }

    if( gridDeltaSep == 0 )
    {
        mxError( "klipReduce", MXE_PARAMNOTSET, "Grid spacing in radius not set (gridDeltaSep)" );
        return -1;
    }

    if( gridDeltaPA == 0 )
    {
        mxError( "klipReduce", MXE_PARAMNOTSET, "Grid spacing in PA not set (grid.deltaPA)" );
        return -1;
    }

    if( gridContrasts.size() == 0 )
    {
        mxError( "klipReduce", MXE_PARAMNOTSET, "Grid contrasts not set (grid.contrasts)" );
        return -1;
    }

    realT x0, y0;

    x0 = -1 * gridCenterSep * sin( mx::math::dtor( gridCenterPA ) );
    y0 = gridCenterSep * cos( mx::math::dtor( gridCenterPA ) );

    //   std::cerr << gridCenterSep << " " << x0 << " " << y0 << "\n";

    int Nrad = 2 * floor( gridHalfWidthSep / gridDeltaSep ) + 1;
    int Npa = 2 * floor( gridHalfWidthPA / gridDeltaPA ) + 1;

    Eigen::Array<realT, -1, -1> sep, pa;

    sep.resize( Nrad, Npa );
    pa.resize( Nrad, Npa );

    realT xp, yp, q, x, y;

    std::vector<realT> xv, yv;

    for( int i = 0; i < Nrad; ++i )
    {
        xp = ( -0.5 * ( Nrad - 1 ) + i ) * gridDeltaSep;

        for( int j = 0; j < Npa; ++j )
        {
            yp = ( -0.5 * ( Npa - 1 ) + j ) * gridDeltaPA;

            q = mx::math::dtor( 90 - gridCenterPA );

            x = ( x0 + xp * cos( q ) + yp * sin( q ) );
            y = ( y0 - xp * sin( q ) + yp * cos( q ) );

            xv.push_back( x );
            yv.push_back( y );

            sep( i, j ) = sqrt( pow( x, 2 ) + pow( y, 2 ) );

            std::cerr << "THIS WON'T WORK UNTIL YOU FIX ANGLEMOD\n";
            exit( 0 );
            // pa(i,j) = mx::math::angleMod<mx::math::degreesT<realT>>(mx::math::rtod( atan2(y, x))  - 90.0);

            // std::cerr << sep(i,j) << " " << pa(i,j) << "\n";

            for( size_t k = 0; k < gridContrasts.size(); ++k )
            {
                // m_obs.m_filesRead = false;
                std::cerr << "THIS WON'T WORK UNTIL YOU FIX FILESREAD IS PROTECTED\n";
                exit( 0 );

                m_obs.m_fakeSep = { sep( i, j ) };
                m_obs.m_fakePA = { pa( i, j ) };
                m_obs.m_fakeContrast = { gridContrasts[k] };

                // std::cerr << sep(i,j) << " " << pa(i,j) << " " << gridContrasts[k] << "\n";
                std::vector<realT> minMaxQ( minRadius.size(), 0 );
                m_obs.regions( minRadius, maxRadius, minMaxQ, minMaxQ );
            }
        }
    }

    mx::fits::fitsFile<realT> ff;

    std::string fn;
    fn = "gridSep.fits";
    if( m_obs.m_outputDir != "" )
        fn = outputDir + "/" + fn;

    ff.write( fn, sep );

    fn = "gridPA.fits";
    if( m_obs.m_outputDir != "" )
        fn = outputDir + "/" + fn;
    ff.write( fn, pa );

    fn = "gridContrasts.dat";
    if( m_obs.m_outputDir != "" )
        fn = outputDir + "/" + fn;
    std::ofstream fout;
    fout.open( fn );
    for( size_t i = 0; i < gridContrasts.size(); ++i )
        fout << gridContrasts[i] << "\n";
    fout.close();

    return 0;
}

int main( int argc, char **argv )
{

#ifdef _OPENMP
    std::cerr << "I think maximum number of threads is: " << omp_get_max_threads() << "\n";
#endif

    klipReduce<float, double> kr;

    try
    {
        kr.main( argc, argv );
    }
    catch( const std::exception &e )
    {
        std::cerr << argv[0] << " exception caught:\n  " << e.what() << "\n";
        return -1;
    }

    return 0;
}
