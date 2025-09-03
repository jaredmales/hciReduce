/** \file HCIobservation.hpp
 * \author Jared R. Males
 * \brief Defines the basic high contrast imaging data type.
 * \ingroup hc_imaging_files
 * \ingroup image_processing_files
 *
 */

#ifndef __HCIobservation_hpp__
#define __HCIobservation_hpp__

#include <vector>
#include <map>
#include <string>
#include <fstream>

#include <sys/stat.h>

#include <mx/mxlib.hpp>

#include <mx/app/appConfigurator.hpp>

#include <mx/mxException.hpp>

#include <mx/math/templateBLAS.hpp>
#include <mx/sys/timeUtils.hpp>
#include <mx/ioutils/fileUtils.hpp>
#include <mx/ioutils/readColumns.hpp>
#include <mx/ioutils/fits/fitsFile.hpp>
#include <mx/ipc/ompLoopWatcher.hpp>

#include <mx/improc/eigenImage.hpp>
#include <mx/improc/eigenCube.hpp>
#include <mx/improc/imageFilters.hpp>
#include <mx/improc/imageMasks.hpp>
#include <mx/improc/imageTransforms.hpp>
#include <mx/improc/imageUtils.hpp>

namespace mx
{

namespace improc
{

/// Namespace for high contrast imaging enums.
/** \ingroup hc_imaging_enums
 */
namespace HCI
{

/// Possible coadding methods
/** \ingroup hc_imaging_enums
 */
enum class coaddMethod
{
    none,   ///< Do not combine the images.
    median, ///< Combine with the median.
    mean,   ///< Combine with the mean.
    invalid ///< An invalid method
};

/// Get the string name of the coaddMethod
/**
 * \returns a string with the name of the coaddMethod
 */
std::string coaddMethodStr( coaddMethod method /**< [in] one of the coaddMethod enum members */ );

/// Get the coaddMethod from the corresponding string name
/**
 * \returns the coaddMethod enum member corresponding to the string name.
 */
coaddMethod coaddMethodStr( const std::string &method /**< [in] the string name of the coaddMethod */ );

/// Mean subtraction methods
/** These control how the data in each search region is centered to meet the PCA
 * requirement. \ingroup hc_imaging_enums
 */
enum class meanSubMethod
{
    none,        ///< No mean subtraction
    meanImage,   ///< The mean image of the data is subtracted from each image
    medianImage, ///< The median image of the data is subtracted from each image
    imageMean,   /**< The mean of each image (within the search region) is
                      subtracted from itself*/
    imageMedian, /**< The median of each image (within the search region) is
                      subtracted from itself*/
    imageMode    /**< The mode of each image (within the search region) is
                      subtracted from itself*/
};

std::string meanSubMethodStr( meanSubMethod method );

meanSubMethod meanSubMethodStr( const std::string &method );

enum class pixelTSNormMethod
{
    none,           ///< no pixel time series norm
    rms,            ///< the rms of the pixel time series
    rmsSigmaClipped ///< the sigma clipped rms of the pixel time series
};

std::string pixelTSNormMethodStr( pixelTSNormMethod method );

pixelTSNormMethod pixelTSNormMethodStr( const std::string &method );

/// Possible combination methods
/** \ingroup hc_imaging_enums
 */
enum class combineMethod
{
    none,     ///< Do not combine the images.
    median,   ///< Combine with the median.
    mean,     ///< Combine with the mean.
    sigmaMean ///< Combine with the sigma clipped mean.
};

/// Get the string name of the combineMethod
/**
 * \returns a string with the name of the combineMethod
 */
std::string combineMethodStr( combineMethod method /**< [in] one of the combineMethod enum members */ );

/// Get the combineMethod from the corresponding string name
/**
 * \returns the combineMethods enum member corresponding to the string name.
 */
combineMethod combineMethodFmStr( const std::string &method /**< [in] the string name of the combineMethod */ );

} // namespace HCI

/// The basic high contrast imaging data type
/** This class manages file reading, resizing, co-adding, pre-processing (masking and filtering),
 * and final image combination.
 *
 * \tparam _realT is the floating point type in which to do all arithmetic.
 * \tparam verboseT sets the verbosity of error reporting
 *
 * \ingroup hc_imaging
 */
template <typename _realT, class verboseT = mx::verbose::vvv>
struct HCIobservation
{

  public:
    /// The arithmetic type used for calculations.  Does not have to match the type in images on disk.
    typedef _realT realT;

    /// The Eigen image array type basted on realT
    typedef Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> imageT;

    typedef mx::fits::fitsFile<realT, verboseT> fitsFileT;
    typedef mx::fits::fitsHeader<verboseT> fitsHeaderT;
    typedef mx::fits::fitsHeaderCard<verboseT> fitsHeaderCardT;

  protected:
    /** \name Input Target Images Configuration
     * Options to control which files are read, how they are read, what meta data is extracted
     * from FITS headers, and image size.
     * @{
     */

    /// Directory to search for input files
    std::string m_directory;

    /// Prefix of the input files. Can be empty if all files in the directory are to be used.
    std::string m_prefix;

    /// Extension of the input files.  Default is .fits.
    std::string m_extension{ ".fits" };

    /// Path to a file containing a list of input files to read.
    /** If directory is also specified as \ref m_directory this must contain paths relative to that directory.
     *  This can be set on construction, configuration, or by calling \ref load_fileList
     */
    std::string m_fileListFile;

    /// The number of files to delete from the front of the list.  Default is 0.
    int m_deleteFront{ 0 };

    /// The number of files to delete from the back of the list.  Default is 0.
    int m_deleteBack{ 0 };

    /// The path to the file containing a list of image quality in 'filename quality' pairs, one entry per image.
    /** If this is not empty and \ref m_qualityThreshold is > 0, then only images where
     * qualityValue >= qualityThreshold are read and processed.
     *
     * The only restriction on m_qualityThreshold is that it is > 0 and is in the same units as the m_qualityThreshold.
     */
    std::string m_qualityFile;

    /// Threshold to apply to quality values read from \ref m_qualityFile.
    /** If <= 0, then thresholding is not performed.
     */
    realT m_qualityThreshold{ 0 };

    /// Perform thresholding only, and print the list of files which pass.
    /** Prints the names and qualities of the files which pass threshold, and then stops.
     *
     */
    bool m_thresholdOnly{ false };

    /// The FITS keyword to use for the image date.
    /** Specifies the keyword corresponding to the date.  This is
     *  the "DATE" keyword for file write time, and usually
     *  "DATE-OBS" for actual observation time.
     *
     * Default is "DATE-OBS".
     *
     * If empty "", then image date is not read.
     */
    std::string m_dateKeyword{ "DATE-OBS" };

    /// Whether or not the date is in ISO 8601 format
    /**
     * Default is true.
     */
    bool m_dateIsISO8601{ true };

    /// If the date is not ISO 8601, this specifies the conversion to Julian Days (e.g. seconds to days)
    /**
     * Default is 1.0.
     */
    realT m_dateUnit{ 1.0 };

    /// Vector of FITS header keywords to read from the files in m_fileList.
    std::vector<std::string> m_keywords;

    /// The max size to read in from the images.
    /** Set to 0 to use images uncut.
     *
     * Image sizes are not increased if this is larger than their size on disk.
     */
    int m_imSize{ 0 };

    ///@}

    ///\name The Input Target Images Data
    /** @{
     */

    /// The list of input files to read.
    /** Constructed by either reading m_fileListFile or searching m_directory for the
     * files specified by m_prefix and m_extension.
     */
    std::vector<std::string> m_fileList;

    /// The target image cube
    eigenCube<realT> m_tgtIms;

    int m_Nims{ 0 };  ///< Number of images in m_tgtIms
    int m_Nrows{ 0 }; ///< Number of rows of the images in m_tgtIms
    int m_Ncols{ 0 }; ///< Number of columns of the images in m_tgtIms
    int m_Npix{ 0 };  ///< Pixels per image, that is Nrows*Ncols

    /// Vector of target image times, in MJD.
    std::vector<double> m_imageMJD;

    /// Vector of FITS headers, one per file, populated with the values for the keywords.
    std::vector<fitsHeaderT> m_heads;

    /// Whether or not the m_fileList has been read.
    bool m_filesRead{ false };

    /// Whether or not the specified files have been deleted from m_fileList
    bool m_filesDeleted{ false };

    ///@}

    /** \name RDI Input Reference Images Configuration
     * Options for Reference Differential Imaging (RDI) input files to control which files are read for
     * the references, how they are read, and what meta data is extracted from FITS headers.
     * @{
     */

    /// Directory to search for input reference files
    std::string m_RDIdirectory;

    /// Prefix of the input reference files. Can be empty if all files in the directory are to be used.
    std::string m_RDIprefix;

    /// Extension of the input reference files.  Default is to match the target input files.
    std::string m_RDIextension{ ".fits" };

    /// Path to a file containing a list of input reference files to read.
    /** If directory is also specified as \ref m_RDIdirectory this must contain paths relative to that directory.
     *  This can be set on construction, configuration, or by calling \ref load_RDIfileList
     */
    std::string m_RDIfileListFile;

    /// The number of reference files to delete from the front of the list.  Default is 0.
    int m_RDIdeleteFront{ 0 };

    /// The number of reference files to delete from the back of the list.  Default is 0.
    int m_RDIdeleteBack{ 0 };

    /// The path to the file containing a list of reference image quality in 'filename quality' pairs, one entry per
    /// image.
    /** If this is not empty and \ref m_RDIqualityThreshold is > 0, then only reference images where
     * RDIqualityValue >= RDIqualityThreshold are read and processed.
     *
     * The only restriction on m_RDIqualityThreshold is that it is > 0 and is in the same units as
     * m_RDIqualityThreshold.
     */
    std::string m_RDIqualityFile;

    /// Threshold to apply to qualityValues read from \ref m_RDIqualityFile.
    /** If <= 0, then thresholding is not performed on reference data.
     */
    realT m_RDIqualityThreshold{ 0 };

    /// The FITS keyword to use for the reference image date.
    /** Specifies the keyword corresponding to the date.  This is
     *  the "DATE" keyword for file write time, and usually
     *  "DATE-OBS" for actual observation time.
     *
     * Default is to follow the main input \ref m_dateKeyword.
     *
     * If empty "", then image date is not read.
     */
    std::string m_RDIdateKeyword{ "DATE-OBS" };

    /// Whether or not the date in reference images is in ISO 8601 format
    /**
     * Default is to follow the main input \ref m_dateIsISO8601
     */
    bool m_RDIdateIsISO8601{ true };

    /// If the reference image date is not ISO 8601, this specifies the conversion to Julian Days (e.g. seconds to days)
    /**
     * The default is follow hte main image \ref m_dateUnit;
     */
    realT m_RDIdateUnit{ 1.0 };

    /// Vector of FITS header keywords to read from the files in m_RDIfileList.
    std::vector<std::string> m_RDIkeywords;

    ///@}

    ///\name RDI Input Reference Images Data
    /** @{
     */

    /// The list of input reference files to read.
    /** Constructed by either reading m_RDIfileListFile or searching m_directory for the
     * files specified by m_RDIprefix and m_RDIextension.
     */
    std::vector<std::string> m_RDIfileList;

    /// The optional reference image cube
    eigenCube<realT> m_refIms;

    /// Vector of reference image times, in MJD.
    std::vector<double> m_RDIimageMJD;

    /// Vector of FITS headers, one per reference file, populated with the values for the keywords.
    std::vector<fitsHeaderT> m_RDIheads;

    /// Whether or not the reference files have been read.
    bool m_RDIfilesRead{ false };

    /// Whether or not the specified files have been deleted from m_RDIfileList
    bool m_RDIfilesDeleted{ false };

    ///@}

    /** \name Mask Configuration
     * A 1/0 mask can be supplied, which is used in pre-processing and in image combination.
     * @{
     */

    /// Specify a mask file to apply to the input images
    /**No mask is applied if this is empty.
     */
    std::string m_maskFile;

    /// Specify a mask file to apply to the reference images
    /**No mask is applied if this is empty, unless m_RDImaskUseInput is true.
     */
    std::string m_RDImaskFile;

    /// Specify that the input mask file should be used for the reference images
    /** This will override m_RDImaskFile.
     */
    bool m_RDImaskUseInput{ false };

    ///@}

    /** \name Mask Data
     *
     * @{
     */

    imageT m_mask;                  ///< The mask

    eigenCube<realT> m_maskCube;    /**< A cube of masks, one for each input image, which may be modified
                                         versions (e.g. rotated) of mask. */

    imageT m_RDImask;               ///< The mask for RDI images

    eigenCube<realT> m_RDImaskCube; /**< A cube of masks, one for each reference image, which may be modified
                                      versions (e.g. rotated) of mask. */

    ///@}

    /** \name Coadding Configuration
     * These parameters control whether and how the images are coadded after being read.  Coadding can
     * be done up to a given number of images, and/or a given elapsed time.
     *
     * Averages the values of given Keywords as well.
     * @{
     */

    /// The method to use for coadding the input images.
    /** Possibilities are
     * - HCI::coaddMethod::none -- [default] do not combine.  This turns off coadding.
     * - HCI::coaddMethod::median --  coadded image is the median
     * - HCI::coaddMethod::mean -- coadded image is the simple mean
     *
     * No other types of combination are currently supported for coadding.
     */
    HCI::coaddMethod m_coaddMethod{ HCI::coaddMethod::none };

    /// Maximum number of images to coadd at a time.
    int m_coaddMaxImno{ 0 };

    /// Maximum elapsed time over which to coadd the images, in seconds.
    realT m_coaddMaxTime{ 0 };

    /// Specify FITS keywords from the input images whose values will be averaged and replaced.
    std::vector<std::string> m_coaddKeywords;

    ///@}
    //-- Coadding Configuration

    /** \name Pre-Processing Configuraton
     * These options control the pre-processing masking and filtering.
     * They are performed in the following order:
     * -# mask applied (enabled by m_preProcess_mask)
     * -# radial profile subtraction (enabled by m_preProcess_subradprof)
     * -# mask applied (enabled by m_preProcess_mask)
     * -# symmetric median unsharp mask (m_preProcess_gaussUSM_fwhm)
     * -# symmetric Gaussian unsharp mask (m_preProcess_gaussUSM_fwhm)
     * -# mask applied (enabled by m_preProcess_mask)
     * -# azimuthal unsharp mask (m_preProcess_azUSM_azW, and m_preProcess_azUSM_radW)
     * -# mask applied (enabled by m_preProcess_mask)
     * -# mean subtraction (enabled by m_preProcess_meanSubMethod)
     * -# mask applied (enabled by m_preProcess_mask)
     * -# pixel time-series normalization (enabled by m_preProcess_pixelTSNormMethod)
     * @{
     */

    bool m_skipPreProcess{ false };         ///< Don't do any of the pre-processing steps (including coadding).

    bool m_preProcess_beforeCoadd{ false }; ///< controls whether pre-processing takes place before or after coadding

    bool m_preProcess_mask{ true };         ///< If true, the mask is applied during each pre-processing step.

    bool m_preProcess_subradprof{ false };  ///< If true, a radial profile is subtracted from each image.

    /// Azimuthal boxcar width for azimuthal unsharp mask [pixels]
    /** If this is 0 then azimuthal-USM is not performed.
     */
    realT m_preProcess_azUSM_azW{ 0 };

    /// Maximum azimuthal boxcar width for azimuthal unsharp mask [degrees]
    /** Limits width close to center, preventing wrap-around.  Default is 45 degrees.  Set to 0 for no maximum.
     */
    realT m_preProcess_azUSM_maxAz{ 45 };

    /// Radial boxcar width for azimuthal unsharp mask [pixels]
    /** If this is 0 then azimuthal-USM is not performed.
     */
    realT m_preProcess_azUSM_radW{ 0 };

    /// Kernel full-width for symmetric box median unsharp mask (USM)
    /** USM is not performed if this is 0.
     */
    int m_preProcess_medianUSM_fwhm{ 0 };

    /// Kernel FWHM for symmetric Gaussian unsharp mask (USM)
    /** USM is not performed if this is 0.
     */
    realT m_preProcess_gaussUSM_fwhm{ 0 };

    /// The mean subtraction method during pre-processing
    /** Can only be none, meanImage, or medianImage
     */
    HCI::meanSubMethod m_preProcess_meanSubMethod{ HCI::meanSubMethod::none };

    /// Specify if each pixel time-series is normalized
    /** This normalizaton is applied after centering. Can have the following values:
     * - <b>HCI::pixelTSNormMethod::none</b>: no normalization (the default)
     * - <b>HCI::pixelTSNormMethod::rms</b>: divide by the time-series rms
     * - <b>HCI::pixelTSNormMethod::rmsSigmaClipped</b>: divide by the sigma-slipped time-series rms.
     *                                                   The sigma is provided by m_preProcess_pixelTSSigma.
     */
    HCI::pixelTSNormMethod m_preProcess_pixelTSNormMethod{ HCI::pixelTSNormMethod::none };

    realT m_pixelTSSigma{ 3 }; ///< Sigma-clipping parameter for pixel time-series normalization

    /// Set path and file prefix to output the pre-processed images.
    /** If empty, then pre-processed images are not output.
     */
    std::string m_preProcess_outputPrefix;

    /// If true, then we stop after pre-processing.
    bool m_preProcess_only{ false };

    ///@}
    //--pre-processing configuration

  public:
    ///\name Construction and Initialization
    /** @{
     */
    /// Default c'tor
    HCIobservation();

    int setupConfig( mx::app::appConfigurator &config );

    int loadConfig( mx::app::appConfigurator &config );

  protected:
    /// Load the file list (internal worker)
    /** Populates the fileList vector by either reading fileListFile (if it is not "") or by
     * searching on disk for files which match the given parameters "directory/prefix*.extension".
     *
     */
    mx::error_t load_fileList( std::vector<std::string> &fileList,
                               const std::string &fileListFile,
                               const std::string &directory,
                               const std::string &prefix,
                               const std::string &extension );

  public:
    /// Load the file list
    /** Populates the \ref m_fileList vector by either reading m_fileListFile (if it is not "") or by
     * searching on disk for files which match the given parameters "m_directory/m_prefix*.m_extension".
     *
     */
    mx::error_t load_fileList();

    /// Load the RDI reference file list
    /** Populates the \ref m_RDIfileList vector by either reading m_RDIfileListFile (if it is not "") or by
     * searching on disk for files which match the given parameters "m_RDIdirectory/m_RDIprefix*.m_RDIextension".
     *
     */
    mx::error_t load_RDIfileList();

    ///@}

    ///\name The Reduced Data
    /** @{
     */
    /// The PSF subtracted images
    /** This is a vector of cubes so that it can contain results from different reductions,
     * e.g. different modes when using KLIP.
     */
    std::vector<eigenCube<realT>> m_psfsub;

    /// The final combined images, one for each cube in psfsub.
    eigenCube<realT> m_finim;

    ///@}

    /// Input
    /** @{ */
    /// Read the list of files, cut to size, and preprocess.
    /**
     * \returns 0 on success, -1 on  error.
     */
    int readFiles();

    /// Perform post-read actions for the target images, for use by derived classes
    virtual int postReadFiles();

    /// Perform post-coadd actions for the target images, for use by derived classes.
    /**
     * \returns 0 on success
     * \returns \<0 on error.
     */
    virtual int postCoadd();

    ///@}

    /// Read the list of reference files, cut to size, and preprocess.
    /** The target files must be read with \ref readFiles() before calling this method.
     *
     * \returns 0 on success, -1 on  error.
     */
    int readRDIFiles();

    /// Perform post-read actions for the RDI images, for use by derived classes
    virtual int postRDIReadFiles();

    /// Perform post-coadd actions, for use by derived classes.
    /** A key example is to update keywords after the averaging occurs in coaddImages().
     *
     * \returns 0 on success
     * \returns \<0 on error.
     */
    virtual int postRDICoadd();

    ///@}
    //--RDI

    /** \name Thresholding
     * Thresholds are applied to a list of files before it is read, based on the image qualities supplied.
     * @{
     */

    /// Read the image qualities from a qualityFile and apply the threshold to a fileList
    /** This is called by readFiles().
     *
     * \returns 0 on success, -1 on  error.
     */
    int threshold( std::vector<std::string> &fileList, ///< [in.out] the fileList to threshold
                   const std::string &qualityFile, ///< [in] the path to the file containing qualities, one per file.
                   realT qualityThreshold          ///< [in] the quality threshold to apply
    );

    ///@}

    /** \name Masking
     * A 1/0 mask can be supplied, which is used in pre-processing and in image combination.
     * @{
     */

    /// Read the mask file, resizing to imSize if needed.
    void readMask();

    /// Populate the mask cube which is used for post-processing.
    /** Derived classes can do this as appropriate, e.g. by rotating the mask.
     * \throws mx::err::invalidconfig if mask is not the same size as the images
     *
     * \todo this should probably be makePostMaskCube
     */
    virtual void makeMaskCube();

    ///@}

    /** \name Coadding
     *
     * @{
     */

    /// Coadd the images
    void coaddImages( HCI::coaddMethod coaddMethod,
                      int coaddMaxImno,
                      int coaddMaxTime,
                      std::vector<std::string> &coaddKeywords,
                      std::vector<double> &imageMJD,
                      std::vector<fitsHeaderT> &heads,
                      eigenCube<realT> &ims );

    ///@} -- coadding

  public:
    /** \name Pre-Processing
     * @{
     */
    /// Do the pre-processing
    void preProcess( eigenCube<realT> &ims /**< [in] the image cube, should be either m_tgtIms or m_refIms */ );

    /// Do mean subtraction as part of pre-processing
    void preProcess_meanSub( eigenCube<realT> &ims /**< [in] the image cube, should be either m_tgtIms or m_refIms */ );

    /// Do pixel time-series normalization as part of pre-processing
    void
    preProcess_pixelTSNorm( eigenCube<realT> &ims /**< [in] the image cube, should be either m_tgtIms or m_refIms */ );

    ///@}

    /** \name Image Combination
     * These options control how the final image combination is performed.
     * @{
     */

    /// Determine how to combine the PSF subtracted images
    /** Possibilities are
     * - HCI::noCombine -- do not combine
     * - HCI::combineMethod::median -- [default] final image is the median
     * - HCI::combineMethod::mean -- final image is the simple mean
     * - HCI::weightedMeanCombine -- final image is the weighted mean.  m_weightFile must be provided.
     * - HCI::combineMethod::sigmaMean -- final image is sigma clipped mean.  If m_sigmaThreshold \<= 0, then it reverts
     * to meanCombine.
     */
    HCI::combineMethod m_combineMethod{ HCI::combineMethod::mean };

    /// Specifies a file containing the image weights, for combining with weighted mean.
    /** This 2-column space-delimited ASCII file containing  filenames and weights. It must be specified before
     * readFiles() is executed.  In readFiles this is loaded after the m_fileList is cutdown and matched to the
     * remaining files.
     */
    std::string m_weightFile;

    /// Vector to hold the image weights read from the m_weightFile.
    /** After readWeights is executed by readFiles, this will contain the normalized weights.
     * \todo check how comboWeights are handled in coadding
     */
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
     */
    int readWeights();

    /// Combine the images into a single final image.
    /** Images are combined by the method specified in \ref combineMethod
     */
    void combineFinim();

    ///@}

    /** \name Output
     * These options control the ouput of the final combined images and the individual PSF subtracted images.
     * @{
     */

    /// Location for temporary auxilliary output files (e.g. masks)
    std::string m_auxDataDir{ "/tmp/hciReduceAux/" };

    /// Whether or not to move the temp. aux files.
    bool m_moveAuxDataDir{ true };

    /// Set whether the final combined image is written to disk
    int m_doWriteFinim{ 1 };

    /// The directory where to write output files.
    std::string m_outputDir;

    /// The base file name of the output final image
    /** The complete name is formed by combining with a sequential number and the ".fits" extension.
     * that is: m_finimName0000.fits.  This behavior can be modified with m_exactFinimName.
     */
    std::string m_finimName{ "finim_" };

    /// Use m_finimName exactly as specified, without appending a number or an extension.
    /** Output is still FITS format, regardless of extension.  This will overwrite
     * an existing file without asking.
     */
    bool m_exactFinimName{ false };

    /// Controls whether or not the individual PSF subtracted images are written to disk.
    /** - true -- write to disk
     * - false -- [default] don't write to disk
     */
    bool m_doOutputPSFSub{ false };

    /// Prefix of the FITS file names used to write individual PSF subtracted images to disk if m_doOutputPSFSub is
    /// true.
    std::string m_PSFSubPrefix;

    /// Output the pre-processed target images
    void outputPreProcessed();

    /// Output the pre-processed reference images
    void outputRDIPreProcessed();

    /// Fill in the HCIobservation standard FITS header
    /**
     */
    void stdFitsHeader( fitsHeaderT &head /**< [in.out] the fistHeader structure which will
                                                             have cards appended to it. */
    );

    /// Write the final combined image to disk
    /**
     */
    void writeFinim( fitsHeaderT *addHead = 0 );

    /// Write the PSF subtracted images to disk
    /**
     */
    void outputPSFSub( fitsHeaderT *addHead = 0 );

    ///@}

    /// Read in already PSF-subtracted files
    /** Used to take up final processing after applying some non-klipReduce processing steps to
     * PSF-subtracted images.
     */
    // int readPSFSub( const std::string &dir, const std::string &prefix, const std::string &ext, size_t nReductions );

    double t_begin{ 0 };
    double t_end{ 0 };

    double t_load_begin{ 0 };
    double t_load_end{ 0 };

    double t_coadd_begin{ 0 };
    double t_coadd_end{ 0 };

    double t_preproc_begin{ 0 };
    double t_preproc_end{ 0 };

    double t_azusm_begin{ 0 };
    double t_azusm_end{ 0 };

    double t_gaussusm_begin{ 0 };
    double t_gaussusm_end{ 0 };

    double t_combo_begin{ 0 };
    double t_combo_end{ 0 };
};

// -- construction and initialization

template <typename _realT, class verboseT>
HCIobservation<_realT, verboseT>::HCIobservation()
{
}

template <typename _realT, class verboseT>
int HCIobservation<_realT, verboseT>::setupConfig( mx::app::appConfigurator &config )
{
    config.add( "input.directory",
                "D",
                "input.directory",
                mx::app::argType::Required,
                "input",
                "directory",
                false,
                "string",
                "Directory to search for input files" );

    config.add( "input.prefix",
                "P",
                "input.prefix",
                mx::app::argType::Required,
                "input",
                "prefix",
                false,
                "string",
                "Prefix of the input files. Can be empty if all files in the directory are to be used." );

    config.add( "input.extension",
                "E",
                "input.extension",
                mx::app::argType::Required,
                "input",
                "extension",
                false,
                "string",
                "Extension of the input files. Default is `.fits`" );

    config.add( "input.fileList",
                "F",
                "input.fileList",
                mx::app::argType::Required,
                "input",
                "fileList",
                false,
                "string",
                "Path to a list of input files to read. If directory is also specified this must contain paths "
                "relative to that directory" );

    config.add( "input.deleteFront",
                "",
                "input.deleteFront",
                mx::app::argType::Required,
                "input",
                "deleteFront",
                false,
                "int",
                "The number of files to delete from the front of the list.  Default is 0." );

    config.add( "input.deleteBack",
                "",
                "input.deleteBack",
                mx::app::argType::Required,
                "input",
                "deleteBack",
                false,
                "int",
                "The number of files to delete from the back of the list.  Default is 0." );

    config.add( "input.qualityFile",
                "",
                "input.qualityFile",
                mx::app::argType::Required,
                "input",
                "qualityFile",
                false,
                "string",
                "The path to the file containing a list of image quality in 'filename number' pairs, one entry per "
                "each image." );

    config.add( "input.qualityThreshold",
                "",
                "input.qualityThreshold",
                mx::app::argType::Required,
                "input",
                "qualityThreshold",
                false,
                "",
                "Threshold to apply to quality values read from the qualityFile. If <= 0, then thresholding is not "
                "performed." );

    config.add( "input.thresholdOnly",
                "",
                "input.thresholdOnly",
                mx::app::argType::True,
                "input",
                "thresholdOnly",
                false,
                "bool",
                "Perform thresholding only, and print the list of files which pass." );

    config.add( "input.dateKeyword",
                "",
                "input.dateKeyword",
                mx::app::argType::Required,
                "input",
                "dateKeyword",
                false,
                "string",
                "Name of the FITS keyword to use for the image date.  Default is `DATE-OBS`" );

    config.add( "input.dateIsISO8601",
                "",
                "input.dateIsISO8601",
                mx::app::argType::True,
                "input",
                "dateIsISO8601",
                false,
                "bool",
                "Whether or not the date is in ISO 8601 format." );

    config.add( "input.dateUnit",
                "",
                "input.dateUnit",
                mx::app::argType::True,
                "input",
                "dateUnit",
                false,
                "bool",
                "If the date is not ISO 8601, this specifies the conversion to Julian Days (e.g. seconds to days)" );

    config.add( "input.imSize",
                "S",
                "input.imSize",
                mx::app::argType::Required,
                "input",
                "imSize",
                false,
                "int",
                "The max size to read in from the images.  Set to 0 to use images uncut." );

    config.add( "input.maskFile",
                "",
                "input.maskFile",
                mx::app::argType::Required,
                "input",
                "maskFile",
                false,
                "string",
                "Path to a file containing a 1/0 mask.  0 pixels are excluded from analysis." );

    config.add( "rdi.directory",
                "",
                "rdi.directory",
                mx::app::argType::Required,
                "rdi",
                "directory",
                false,
                "string",
                "Directory to search for RDI reference files" );

    config.add( "rdi.prefix",
                "",
                "rdi.prefix",
                mx::app::argType::Required,
                "rdi",
                "prefix",
                false,
                "string",
                "Prefix of the RDI reference files.  Can be empty if all files in the RDI directory are to be used." );

    config.add( "rdi.extension",
                "",
                "rdi.extension",
                mx::app::argType::Required,
                "rdi",
                "extension",
                false,
                "string",
                "Extension of the input reference files.  Default is to match the main input.extension." );

    config.add( "rdi.fileList",
                "",
                "rdi.fileList",
                mx::app::argType::Required,
                "rdi",
                "fileList",
                false,
                "string",
                "Path to a list of input reference files to read. If RDI directory is also specified this must contain "
                "paths relative to that directory" );

    config.add( "rdi.deleteFront",
                "",
                "rdi.deleteFront",
                mx::app::argType::Required,
                "rdi",
                "deleteFront",
                false,
                "int",
                "The number of files to delete from the front of the RDI file list.  Default is 0." );

    config.add( "rdi.deleteBack",
                "",
                "rdi.deleteBack",
                mx::app::argType::Required,
                "rdi",
                "deleteBack",
                false,
                "int",
                "The number of files to delete from the back of the RDI file list.  Default is 0." );

    config.add( "rdi.qualityFile",
                "",
                "rdi.qualityFile",
                mx::app::argType::Required,
                "rdi",
                "qualityFile",
                false,
                "string",
                "The path to the file containing a list of image quality for the RDI images in 'filename number' "
                "pairs, one entry per each image." );

    config.add( "rdi.qualityThreshold",
                "",
                "rdi.qualityThreshold",
                mx::app::argType::Required,
                "rdi",
                "qualityThreshold",
                false,
                "",
                "Threshold to apply to quality values read from the RDIqualityFile. If <= 0, then thresholding is not "
                "performed." );

    config.add( "rdi.dateKeyword",
                "",
                "rdi.dateKeyword",
                mx::app::argType::Required,
                "rdi",
                "dateKeyword",
                false,
                "string",
                "Name of the FITS keyword to use for the reference image date.  Default is to follow the main "
                "input.dateKeyword." );

    config.add( "rdi.dateIsISO8601",
                "",
                "rdi.dateIsISO8601",
                mx::app::argType::True,
                "rdi",
                "dateIsISO8601",
                false,
                "bool",
                "Whether or not the reference image date is in ISO 8601 format. Default is to follow the main "
                "input.dateIsISO8601." );

    config.add( "rdi.dateUnit",
                "",
                "rdi.dateUnit",
                mx::app::argType::True,
                "rdi",
                "dateUnit",
                false,
                "bool",
                "If the reference image date is not ISO 8601, this specifies the conversion to Julian Days (e.g. "
                "seconds to days).  Default is to follow the main input.dateUnits." );

    config.add( "rdi.maskFile",
                "",
                "rdi.maskFile",
                mx::app::argType::Required,
                "rdi",
                "maskFile",
                false,
                "string",
                "Path to a file containing a 1/0 mask for the reference images. "
                "0 pixels are excluded from analysis.  Defaults to using input.maskFile.  To mask the "
                "main input images but not the references you must set this to empty " );

    config.add( "rdi.useInputMask",
                "",
                "rdi.useInputMask",
                mx::app::argType::Required,
                "rdi",
                "useInputMask",
                false,
                "bool",
                "If true then the main input.maskFile is used for masking "
                "the reference images. This overrides rdi.maskFile" );

    config.add( "coadd.method",
                "",
                "coadd.method",
                mx::app::argType::Required,
                "coadd",
                "method",
                false,
                "string",
                "The method to use for coadding the input images.  Options are none (default), median, and mean.  If "
                "none no coadding is performed." );

    config.add( "coadd.maxImno",
                "",
                "coadd.maxImno",
                mx::app::argType::Required,
                "coadd",
                "maxImno",
                false,
                "int",
                "Maximum number of images to coadd at a time." );

    config.add( "coadd.maxTime",
                "",
                "coadd.maxTime",
                mx::app::argType::Required,
                "coadd",
                "maxTime",
                false,
                "float",
                "Maximum elapsed time over which to coadd the images, in seconds." );

    config.add( "coadd.keywords",
                "",
                "coadd.keywords",
                mx::app::argType::Required,
                "coadd",
                "keywords",
                false,
                "vector<string>",
                "Specify FITS keywords from the input images whose values will be averaged and replaced." );

    config.add( "preProcess.skip",
                "",
                "preProcess.skip",
                mx::app::argType::True,
                "preProcess",
                "skip",
                false,
                "bool",
                "If true, then pre-processing is skipped.  Default is false." );

    config.add( "preProcess.beforeCoadd",
                "",
                "preProcess.beforeCoadd",
                mx::app::argType::True,
                "preProcess",
                "beforeCoadd",
                false,
                "bool",
                "Controls whether pre-processing takes place before (true) or after (false, default) coadding." );

    config.add( "preProcess.mask",
                "",
                "preProcess.mask",
                mx::app::argType::True,
                "preProcess",
                "mask",
                false,
                "string",
                "Determines if mask is applied for pre-processing." );

    config.add( "preProcess.subradprof",
                "",
                "preProcess.subradprof",
                mx::app::argType::True,
                "preProcess",
                "subradprof",
                false,
                "bool",
                "If true, the radial profile is subtracted from each image." );

    config.add( "preProcess.azUSM_azW",
                "",
                "preProcess.azUSM_azW",
                mx::app::argType::Required,
                "preProcess",
                "azUSM_azW",
                false,
                "float",
                "The azimuth USM boxcar azimuthal width in pixels.  Enabled if both azW and radW are non-zero." );

    config.add( "preProcess.azUSM_azW",
                "",
                "preProcess.azUSM_azW",
                mx::app::argType::Required,
                "preProcess",
                "azUSM_azW",
                false,
                "float",
                "The azimuth USM boxcar azimuthal width in pixels.  Enabled if both azW and radW are non-zero." );

    config.add( "preProcess.azUSM_maxAz",
                "",
                "preProcess.azUSM_maxAz",
                mx::app::argType::Required,
                "preProcess",
                "azUSM_maxAz",
                false,
                "float",
                "Maximum azimuthal boxcar width for azimuthal unsharp mask in degrees. Limits width close to center, "
                "preventing wrap-around.  Default is 45 degrees.  Set to 0 for no maximum." );

    config.add( "preProcess.azUSM_radW",
                "",
                "preProcess.azUSM_radW",
                mx::app::argType::Required,
                "preProcess",
                "azUSM_radW",
                false,
                "float",
                "The azimuth USM boxcar radial width in pixels.  Enabled if both azW and radW are non-zero." );

    config.add( "preProcess.medianUSM_fwhm",
                "",
                "preProcess.medianUSM_fwhm",
                mx::app::argType::Required,
                "preProcess",
                "medianUSM_fwhm",
                false,
                "float",
                "The median USM kernel full-width at half max.  Enabled if non-zero." );

    config.add( "preProcess.gaussUSM_fwhm",
                "",
                "preProcess.gaussUSM_fwhm",
                mx::app::argType::Required,
                "preProcess",
                "gaussUSM_fwhm",
                false,
                "float",
                "The gaussian USM kernel full-width at half max.  Enabled if non-zero." );

    config.add( "preProcess.meanSubMethod",
                "",
                "preProcess.meanSubMethod",
                mx::app::argType::Required,
                "preProcess",
                "meanSubMethod",
                false,
                "string",
                "The mean subtraction method during pre-processing. Options are none, meanImage, and medianImage." );

    config.add( "preProcess.pixelTSNormMethod",
                "",
                "preProcess.pixelTSNormMethod",
                mx::app::argType::Required,
                "preProcess",
                "pixelTSNormMethod",
                false,
                "string",
                "The pixel time-series normalization method during pre-processing. "
                "Options are none, rms, rsmSigmaClipped." );

    config.add( "preProcess.outputPrefix",
                "",
                "preProcess.outputPrefix",
                mx::app::argType::Required,
                "preProcess",
                "outputPrefix",
                false,
                "string",
                "If not empty, then this prefix (which should be a full path) is added to file names and the "
                "pre-processed images are output" );

    config.add( "preProcess.only",
                "",
                "preProcess.only",
                mx::app::argType::True,
                "preProcess",
                "only",
                false,
                "bool",
                "If true, stop after pre-processing.  Default is false." );

    return 0;
}

template <typename _realT, class verboseT>
int HCIobservation<_realT, verboseT>::loadConfig( mx::app::appConfigurator &config )
{
    config( m_directory, "input.directory" );
    config( m_prefix, "input.prefix" );
    config( m_extension, "input.extension" );

    config( m_fileListFile, "input.fileList" );

    config( m_deleteFront, "input.deleteFront" );
    config( m_deleteBack, "input.deleteBack" );

    config( m_qualityFile, "input.qualityFile" );
    config( m_qualityThreshold, "input.qualityThreshold" );
    config( m_thresholdOnly, "input.thresholdOnly" );

    config( m_dateKeyword, "input.dateKeyword" );
    config( m_dateIsISO8601, "input.dateIsISO8601" );
    config( m_dateUnit, "input.dateUnit" );

    config( m_imSize, "input.imSize" );

    config( m_maskFile, "input.maskFile" );

    config( m_RDIdirectory, "rdi.directory" );
    config( m_RDIprefix, "rdi.prefix" );
    config( m_RDIextension, "rdi.extension" );

    config( m_RDIfileListFile, "rdi.fileList" );

    config( m_RDIdeleteFront, "rdi.deleteFront" );
    config( m_RDIdeleteBack, "rdi.deleteBack" );

    config( m_RDIqualityFile, "rdi.qualityFile" );
    config( m_RDIqualityThreshold, "rdi.qualityThreshold" );

    config( m_RDIdateKeyword, "rdi.dateKeyword" );
    config( m_RDIdateIsISO8601, "rdi.dateIsISO8601" );
    config( m_RDIdateUnit, "rdi.dateUnit" );

    config( m_RDImaskFile, "rdi.maskFile" );

    config( m_RDImaskUseInput, "rdi.useInputMask" );

    std::string coaddMethodStr = HCI::coaddMethodStr( m_coaddMethod ); // get default
    config( coaddMethodStr, "coadd.method" );
    try
    {
        m_coaddMethod = HCI::coaddMethodStr( coaddMethodStr );
    }
    catch( ... )
    {
        std::throw_with_nested(
            mx::err::invalidconfig( "HCIobservation::loadConfig", __FILE__, __LINE__, "invalid coadd method" ) );
    }

    config( m_coaddMaxImno, "coadd.maxImno" );
    config( m_coaddMaxTime, "coadd.maxTime" );
    config( m_coaddKeywords, "coadd.keywords" );

    config( m_preProcess_beforeCoadd, "preProcess.beforeCoadd" );
    config( m_preProcess_mask, "preProcess.mask" );
    config( m_preProcess_subradprof, "preProcess.subradprof" );
    config( m_preProcess_azUSM_azW, "preProcess.azUSM_azW" );
    config( m_preProcess_azUSM_maxAz, "preProcess.azUSM_maxAz" );
    config( m_preProcess_azUSM_radW, "preProcess.azUSM_radW" );
    config( m_preProcess_medianUSM_fwhm, "preProcess.medianUSM_fwhm" );
    config( m_preProcess_gaussUSM_fwhm, "preProcess.gaussUSM_fwhm" );

    std::string ppmsm = HCI::meanSubMethodStr( m_preProcess_meanSubMethod );
    config( ppmsm, "preProcess.meanSubMethod" );
    try
    {
        m_preProcess_meanSubMethod = HCI::meanSubMethodStr( ppmsm );

        if( m_preProcess_meanSubMethod != HCI::meanSubMethod::none &&
            m_preProcess_meanSubMethod != HCI::meanSubMethod::meanImage &&
            m_preProcess_meanSubMethod != HCI::meanSubMethod::medianImage )
        {
            std::string msg = "Mean subtraction by " + HCI::meanSubMethodStr( m_preProcess_meanSubMethod );
            msg += " can't be done in pre-processing. Only meanImage or medianImage can be used in pre.";
            mxThrowException( mx::err::invalidconfig, "HCIobservation::loadConfig", msg );
        }
    }
    catch( ... )
    {
        std::throw_with_nested( mx::err::invalidconfig( "HCIobservation::loadConfig",
                                                        __FILE__,
                                                        __LINE__,
                                                        "invalid pre-processing mean subtraction method" ) );
    }

    std::string ptsnm = HCI::pixelTSNormMethodStr( m_preProcess_pixelTSNormMethod );
    config( ptsnm, "preProcess.pixelTSNormMethod" );
    try
    {
        m_preProcess_pixelTSNormMethod = HCI::pixelTSNormMethodStr( ptsnm );
    }
    catch( ... )
    {
        std::throw_with_nested( mx::err::invalidconfig( "HCIobservation::loadConfig",
                                                        __FILE__,
                                                        __LINE__,
                                                        "invalid pixel time-series normalization method" ) );
    }

    config( m_preProcess_outputPrefix, "preProcess.outputPrefix" );

    config( m_preProcess_only, "preProcess.only" );
    config( m_skipPreProcess, "preProcess.skip" );

    return 0;
}

template <typename _realT, class verboseT>
mx::error_t HCIobservation<_realT, verboseT>::load_fileList( std::vector<std::string> &fileList,
                                                             const std::string &fileListFile,
                                                             const std::string &directory,
                                                             const std::string &prefix,
                                                             const std::string &extension )
{
    if( fileListFile != "" )
    {
        fileList.clear(); // otherwise readColumns appends

        mx::error_t errc = ioutils::readColumns( fileListFile, fileList );
        if( errc != mx::error_t::noerror )
        {
            return mx::error_report<verboseT>( errc, "error reading " + fileListFile );
        }

        if( directory != "" )
        {
            std::string dir = directory;
            if( dir.back() != '/' )
            {
                dir += '/';
            }

            for( auto &&fpath : m_fileList )
            {
                fpath = dir + fpath;
            }
        }
    }
    else
    {
        mx::error_t errc = ioutils::getFileNames( fileList, directory, prefix, "", extension );
        if( errc != mx::error_t::noerror )
        {
            return mx::error_report<verboseT>(
                errc,
                std::format( "error getting file names for {}/{}*.{}", directory, prefix, extension ) );
        }
    }

    return mx::error_t::noerror;
}

template <typename _realT, class verboseT>
mx::error_t HCIobservation<_realT, verboseT>::load_fileList()
{
    mx::error_t errc = load_fileList( m_fileList, m_fileListFile, m_directory, m_prefix, m_extension );
    m_filesDeleted = false;

    if( errc != mx::error_t::noerror )
    {
        return mx::error_report<verboseT>( errc, "error loading file list" );
    }
    return mx::error_t::noerror;
}

template <typename _realT, class verboseT>
mx::error_t HCIobservation<_realT, verboseT>::load_RDIfileList()
{
    mx::error_t errc = load_fileList( m_RDIfileList, m_RDIfileListFile, m_RDIdirectory, m_RDIprefix, m_RDIextension );
    m_RDIfilesDeleted = false;

    if( errc != mx::error_t::noerror )
    {
        return mx::error_report<verboseT>( errc, "error loading file list" );
    }
    return mx::error_t::noerror;
}

// --< construction and initialization

template <typename realT, class verboseT>
int HCIobservation<realT, verboseT>::readFiles()
{
    if( m_fileList.size() == 0 )
    {
        mxThrowException( err::invalidconfig,
                          "HCIobservation<realT, verboseT>::readFiles",
                          "The target fileList has 0 length, there are no files to be read." );
    }

    // First make the list deletions
    if( !m_filesDeleted )
    {
        if( m_deleteFront > 0 )
        {
            m_fileList.erase( m_fileList.begin(), m_fileList.begin() + m_deleteFront );
        }

        if( m_deleteBack > 0 )
        {
            m_fileList.erase( m_fileList.end() - m_deleteBack, m_fileList.end() );
        }
        m_filesDeleted = true;
    }

    if( m_fileList.size() == 0 )
    {
        mxThrowException( err::invalidconfig,
                          "HCIobservation<realT, verboseT>::readFiles",
                          "The target fileList has 0 length, there are no files to be read after deletions." );
    }

    if( m_qualityFile != "" )
    {
        std::cerr << "Thresholding target images...";
        size_t origsize = m_fileList.size();

        if( threshold( m_fileList, m_qualityFile, m_qualityThreshold ) < 0 )
            return -1;

        if( m_fileList.size() == 0 )
        {
            mxThrowException( err::invalidconfig,
                              "HCIobservation<realT, verboseT>::readFiles",
                              "The fileList has 0 length, there are no files to be read after thresholding." );
        }

        std::cerr << "Done.  Selected " << m_fileList.size() << " out of " << origsize << "\n";

        if( m_thresholdOnly )
        {
            std::cout << "#Files which passed thresholding:\n";
            for( size_t i = 0; i < m_fileList.size(); ++i )
            {
                std::cout << m_fileList[i] << "\n";
            }

            exit( 0 );
        }
    }

    if( m_weightFile != "" )
    {
        if( readWeights() < 0 )
            return -1;
    }

    /*----- Append the HCI keywords to propagate them if needed -----*/

    fitsHeaderT head;

    if( m_dateKeyword != "" )
        head.append( m_dateKeyword );

    for( size_t i = 0; i < m_keywords.size(); ++i )
    {
        if( head.count( m_keywords[i] ) == 0 )
            head.append( m_keywords[i] );
    }

    m_heads.clear(); // This is necessary to make sure heads.resize() copies head on a 2nd call
    m_heads.resize( m_fileList.size(), head );

    /*----- Read in first image to test size -----*/

    Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> im;

    fitsFileT f( m_fileList[0] );

    f.read( im );

    // We set imSize to match the first image, but we make it a square.
    if( m_imSize == 0 )
    {
        m_imSize = im.rows();
        if( m_imSize > im.cols() )
        {
            m_imSize = im.cols();
        }
    }
    else
    {
        // Now make sure we don't read too much.
        if( m_imSize > im.rows() )
        {
            m_imSize = im.rows();
        }
        if( m_imSize > im.cols() )
        {
            m_imSize = im.cols();
        }
    }

    // And now set the read size so we only read what we want.
    // the +0.1 is just to make sure we don't have a problem with precision (we shouldn't)/
    f.setReadSize( floor( 0.5 * ( im.rows() - 1 ) - 0.5 * ( m_imSize - 1 ) + 0.1 ),
                   floor( 0.5 * ( im.cols() - 1.0 ) - 0.5 * ( m_imSize - 1.0 ) + 0.1 ),
                   m_imSize,
                   m_imSize );

    im.resize( m_imSize, m_imSize );

    std::cerr << "image size is " << m_imSize << "\n";

    /**** Right here is we go to coadd.
     */

    // But if not we just do it
    m_tgtIms.resize( im.rows(), im.cols(), m_fileList.size() );

    t_load_begin = sys::get_curr_time();

    f.read( m_tgtIms.data(), m_heads, m_fileList );

    f.setReadSize();

    /* read in the image timestamps, depending on how MJD is stored in the header */
    if( m_dateKeyword != "" )
    {
        m_imageMJD.resize( m_heads.size() );

        if( m_dateIsISO8601 )
        {
            for( size_t i = 0; i < m_imageMJD.size(); ++i )
            {
                m_imageMJD[i] = sys::ISO8601date2mjd( m_heads[i][m_dateKeyword].String() );
            }
        }
        else
        {
            for( size_t i = 0; i < m_imageMJD.size(); ++i )
            {
                m_imageMJD[i] = m_heads[i][m_dateKeyword].template value<realT>() * m_dateUnit;
            }
        }
    }

    t_load_end = sys::get_curr_time();

    m_Nims = m_tgtIms.planes();
    m_Nrows = m_tgtIms.rows();
    m_Ncols = m_tgtIms.cols();
    m_Npix = m_tgtIms.rows() * m_tgtIms.cols();

    std::cerr << "loading complete\n";

    std::cerr << "zero-ing NaNs\n";
    zeroNaNCube( m_tgtIms );
    std::cerr << "done\n";

    /*** Now do the post-read actions ***/
    if( postReadFiles() < 0 )
    {
        return -1;
    }

    /*** Read in the mask if present ***/
    readMask();

    /*** Now begin processing ***/
    if( !m_skipPreProcess )
    {
        /*** Now do any pre-processing ***/
        if( m_preProcess_beforeCoadd )
        {
            preProcess( m_tgtIms );
        }

        if( m_coaddMethod != HCI::coaddMethod::none )
        {
            std::cerr << "Coadding target images...\n";
            coaddImages( m_coaddMethod,
                         m_coaddMaxImno,
                         m_coaddMaxTime,
                         m_coaddKeywords,
                         m_imageMJD,
                         m_heads,
                         m_tgtIms );

            m_Nims = m_tgtIms.planes();
            m_Nrows = m_tgtIms.rows();
            m_Ncols = m_tgtIms.cols();
            m_Npix = m_tgtIms.rows() * m_tgtIms.cols();

            std::cerr << "number of target images after coadding: " << m_Nims << "\n";

            if( postCoadd() < 0 )
            {
                std::cerr << "Post coadd error " << __FILE__ << " " << __LINE__ << "\n";
                return -1;
            }
            std::cerr << "Done.\n";

            // Re-make the mask cube if we coadded...
            if( m_maskFile != "" )
            {
                makeMaskCube();
            }
        }

        /*** Now do any pre-processing if not done already***/
        if( !m_preProcess_beforeCoadd )
        {
            preProcess( m_tgtIms );
        }

        outputPreProcessed();
    }

    m_filesRead = true;

    return 0;
} // readFiles()

template <typename _realT, class verboseT>
int HCIobservation<_realT, verboseT>::postReadFiles()
{
    return 0;
}

template <typename _realT, class verboseT>
int HCIobservation<_realT, verboseT>::postCoadd()
{
    return 0;
}

//------------------- readRDIFiles
template <typename _realT, class verboseT>
int HCIobservation<_realT, verboseT>::readRDIFiles()
{

    /* First check if the target files have been read */
    if( m_Nrows == 0 || m_Ncols == 0 )
    {
        mxError( "HCIobservation", MXE_PARAMNOTSET, "The target image size must be set before reading RDI files." );
        return -1;
    }

    if( m_RDIfileList.size() == 0 )
    {
        mxError( "HCIobservation", MXE_FILENOTFOUND, "The RDI fileList has 0 length, there are no files to be read." );
        return -1;
    }

    // First make the list deletions
    if( !m_RDIfilesDeleted )
    {
        if( m_RDIdeleteFront > 0 )
        {
            m_RDIfileList.erase( m_RDIfileList.begin(), m_RDIfileList.begin() + m_RDIdeleteFront );
        }

        if( m_RDIdeleteBack > 0 )
        {
            m_RDIfileList.erase( m_RDIfileList.end() - m_RDIdeleteBack, m_RDIfileList.end() );
        }
        m_RDIfilesDeleted = true;
    }

    if( m_RDIfileList.size() == 0 )
    {
        mxError( "HCIobservation",
                 MXE_FILENOTFOUND,
                 "The RDI fileList has 0 length, there are no files to be read after deletions." );
        return -1;
    }

    if( m_RDIqualityFile != "" )
    {
        std::cerr << "Thresholding RDI images...";
        size_t origsize = m_RDIfileList.size();

        if( threshold( m_RDIfileList, m_RDIqualityFile, m_RDIqualityThreshold ) < 0 )
            return -1;

        if( m_RDIfileList.size() == 0 )
        {
            mxError( "HCIobservation",
                     MXE_FILENOTFOUND,
                     "The fileList has 0 length, there are no files to be read after thresholding." );
            return -1;
        }

        std::cerr << "Done.  Selected " << m_RDIfileList.size() << " out of " << origsize << "\n";
    }

    /*----- Append the HCI keywords to propagate them if needed -----*/

    fitsHeaderT head;

    if( m_dateKeyword != "" )
        head.append( m_dateKeyword ); // Currently assuming the MJD keyword will be the same

    for( size_t i = 0; i < m_RDIkeywords.size(); ++i )
    {
        head.append( m_RDIkeywords[i] );
    }

    m_RDIheads.clear(); // This is necessary to make sure heads.resize() copies head on a 2nd call
    m_RDIheads.resize( m_RDIfileList.size(), head );

    /*----- Read in first image to adjust size ----*/
    Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> im;

    fitsFileT f( m_RDIfileList[0] );

    f.read( im );

    if( im.rows() < m_imSize || im.cols() < m_imSize )
    {
        mxError( "HCIobservation", MXE_SIZEERR, "The reference images are too small, do not match the target images." );
        return -1;
    }

    // And now set the read size so we only read what we want.
    // the +0.1 is just to make sure we don't have a problem with precision (we shouldn't)/
    f.setReadSize( floor( 0.5 * ( im.rows() - 1 ) - 0.5 * ( m_imSize - 1 ) + 0.1 ),
                   floor( 0.5 * ( im.cols() - 1.0 ) - 0.5 * ( m_imSize - 1.0 ) + 0.1 ),
                   m_imSize,
                   m_imSize );

    m_refIms.resize( m_imSize, m_imSize, m_RDIfileList.size() );

    t_load_begin = sys::get_curr_time();

    f.read( m_refIms.data(), m_RDIheads, m_RDIfileList );

    f.setReadSize();

    if( m_dateKeyword != "" )
    {
        m_RDIimageMJD.resize( m_RDIheads.size() );

        if( m_dateIsISO8601 )
        {
            for( size_t i = 0; i < m_RDIimageMJD.size(); ++i )
            {
                m_RDIimageMJD[i] = sys::ISO8601date2mjd( m_RDIheads[i][m_dateKeyword].String() );
            }
        }
        else
        {
            for( size_t i = 0; i < m_RDIimageMJD.size(); ++i )
            {
                m_RDIimageMJD[i] = m_RDIheads[i][m_dateKeyword].template value<realT>() * m_dateUnit;
            }
        }
    }

    t_load_end = sys::get_curr_time();

    std::cerr << "loading complete\n";

    std::cerr << "zero-ing NaNs\n";
    zeroNaNCube( m_refIms );
    std::cerr << "Done.\n";

    /*** Now do the post-read actions ***/
    if( postRDIReadFiles() < 0 )
        return -1;

    /*** Now begin processing ***/
    if( !m_skipPreProcess )
    {
        /*** Now do any pre-processing ***/
        if( m_preProcess_beforeCoadd )
        {
            preProcess( m_refIms );
        }

        if( m_coaddMethod != HCI::coaddMethod::none )
        {
            std::cerr << "Coadding reference images...\n";
            coaddImages( m_coaddMethod,
                         m_coaddMaxImno,
                         m_coaddMaxTime,
                         m_coaddKeywords,
                         m_RDIimageMJD,
                         m_RDIheads,
                         m_refIms );

            std::cerr << "number of reference images after coadding: " << m_refIms.planes() << "\n";

            if( postRDICoadd() < 0 )
            {
                std::cerr << "Post coadd error " << __FILE__ << " " << __LINE__ << "\n";
                return -1;
            }
            std::cerr << "Done.\n";
        }

        /*** Now do any pre-processing if not done already***/
        if( !m_preProcess_beforeCoadd )
        {
            preProcess( m_refIms );
        }

        // outputRDIPreProcessed();
    }

    m_RDIfilesRead = true;

    return 0;
} // readRDIFiles()

template <typename _realT, class verboseT>
int HCIobservation<_realT, verboseT>::postRDIReadFiles()
{
    return 0;
}

template <typename _realT, class verboseT>
int HCIobservation<_realT, verboseT>::postRDICoadd()
{
    return 0;
}

template <typename _realT, class verboseT>
int HCIobservation<_realT, verboseT>::threshold( std::vector<std::string> &fileList,
                                                 const std::string &qualityFile,
                                                 realT qualityThreshold )
{
    if( qualityFile == "" )
    {
        mxError( "HCIobservation::threshold", MXE_PARAMNOTSET, "qualityFile not set" );
        return -1;
    }

    int origsize = fileList.size();

    std::vector<std::string> qfileNames;
    std::vector<realT> imQ;

    // Read the quality file and load it into a map
    ioutils::readColumns( qualityFile, qfileNames, imQ );

    std::map<std::string, realT> quality;
    for( size_t i = 0; i < qfileNames.size(); ++i )
        quality[ioutils::pathFilename( qfileNames[i].c_str() )] = imQ[i];

    realT q;

    for( size_t i = 0; i < fileList.size(); ++i )
    {
        try
        {
            q = quality.at( ioutils::pathFilename( fileList[i].c_str() ) );
        }
        catch( ... )
        {
            q = qualityThreshold - 1; // Cause it to be erased
        }

        if( q < qualityThreshold )
        {
            fileList.erase( fileList.begin() + i );
            --i;
        }
    }

    return 0;
}

template <typename _realT, class verboseT>
void HCIobservation<_realT, verboseT>::coaddImages( HCI::coaddMethod coaddMethod,
                                                    int coaddMaxImno,
                                                    int coaddMaxTime,
                                                    std::vector<std::string> &coaddKeywords,
                                                    std::vector<double> &imageMJD,
                                                    std::vector<fitsHeaderT> &heads,
                                                    eigenCube<realT> &ims )
{
    std::cerr << "***************************************************************\n";
    std::cerr << "                       *** WARNING ***                         \n";
    std::cerr << "       coadding is poorly tested.  proceed with caution.       \n";
    std::cerr << "***************************************************************\n";

    // Validate setup
    if( coaddMaxImno <= 0 && coaddMaxTime <= 0 )
        return;

    // Validate combine method
    if( coaddMethod == HCI::coaddMethod::none )
    {
        return;
    }

    int Nims = ims.planes();
    int Nrows = ims.rows();
    int Ncols = ims.cols();

    t_coadd_begin = sys::get_curr_time();

    std::vector<imageT> coadds;
    std::vector<std::vector<std::string>> coaddFileNames;

    // We do all math here in double, to avoid losing precision
    std::vector<double> avgMJD;
    std::vector<double> startMJD;
    std::vector<double> endMJD;

    std::vector<std::vector<double>> avgVals;
    std::vector<std::vector<double>> startVals;
    std::vector<std::vector<double>> endVals;

    // Index range of images for next coadd
    int im0, imF;
    im0 = 0;
    imF = 1;

    // Accumulate images to coadd into a cube
    eigenCube<realT> imsToCoadd;

    // The filenames of the images included in the coadd
    std::vector<std::string> imsCoadded;

    // Temporary for combination.
    imageT coadd;

    // Accumulate values
    double initMJD;

    std::vector<double> initVals;
    initVals.resize( coaddKeywords.size() );

    std::vector<double> startVal;
    startVal.resize( coaddKeywords.size() );

    std::vector<double> endVal;
    endVal.resize( coaddKeywords.size() );

    while( im0 < Nims )
    {
        // Initialize the accumulators
        initMJD = imageMJD[im0];
        startMJD.push_back( initMJD );
        endMJD.push_back( initMJD );

        for( size_t i = 0; i < coaddKeywords.size(); ++i )
        {
            initVals[i] = heads[im0][coaddKeywords[i]].template value<double>();
            startVal[i] = initVals[i];
            endVal[i] = initVals[i];
        }

        // Now increment imF, then test whether each variable is now outside the range
        bool increment = true;
        while( increment == true )
        {
            ++imF;

            if( imF >= Nims ) // out of images
            {
                imF = Nims;
                increment = false;
                break;
            }

            if( imF - im0 > coaddMaxImno && coaddMaxImno > 0 ) // too many images in coadd
            {
                --imF;
                increment = false;
                break;
            }

            ///\todo this should really include end of exposure too.
            if( ( imageMJD[imF] - imageMJD[im0] ) * 86400.0 > coaddMaxTime && coaddMaxTime > 0 )
            {
                --imF;
                increment = false;
                break;
            }

        } // while(increment == true)
        // At this point, imF is the first image NOT to include in the next coadd.

        // Now actually do the accumulation
        ///\todo this needs to handle averaging of angles
        for( int imno = im0 + 1; imno < imF; ++imno )
        {
            initMJD += imageMJD[imno];
            endMJD.back() = imageMJD[imno]; // After the last one, this will be the last one

            for( size_t i = 0; i < coaddKeywords.size(); ++i )
            {
                endVal[i] = heads[imno][coaddKeywords[i]]
                                .template value<double>(); // After the last one, this will be the last one
                initVals[i] += endVal[i];
            }
        }

        // And then turn them into an average
        initMJD /= ( imF - im0 );
        for( size_t i = 0; i < coaddKeywords.size(); ++i )
        {
            initVals[i] /= ( imF - im0 );
        }

        // Extract the images into the temporary
        imsToCoadd.resize( Nrows, Ncols, imF - im0 );
        imsCoadded.resize( imF - im0 );
        for( int i = 0; i < ( imF - im0 ); ++i )
        {
            imsToCoadd.image( i ) = ims.image( im0 + i );
            imsCoadded[i] = ioutils::pathFilename( m_fileList[im0 + i] );
        }

        // Here do the combine and insert into the vector
        if( coaddMethod == HCI::coaddMethod::median )
        {
            imsToCoadd.median( coadd );
        }
        if( coaddMethod == HCI::coaddMethod::mean )
        {
            imsToCoadd.mean( coadd );
        }
        coadds.push_back( coadd );
        coaddFileNames.push_back( imsCoadded );

        // Insert the new averages
        avgMJD.push_back( initMJD );
        avgVals.push_back( initVals );
        startVals.push_back( startVal );
        endVals.push_back( endVal );

        im0 = imF;
        imF = im0 + 1;
    } // while(im0 < Nims)

    // Now resize ims and copy the coadds to the new cube
    ims.resize( Nrows, Ncols, coadds.size() );
    Nims = coadds.size();

    for( int i = 0; i < Nims; ++i )
    {
        ims.image( i ) = coadds[i];
    }

    // Now deal with imageMJD and headers
    imageMJD.erase( imageMJD.begin() + Nims, imageMJD.end() );
    heads.erase( heads.begin() + Nims, heads.end() );

    for( int i = 0; i < Nims; ++i )
    {
        imageMJD[i] = avgMJD[i];
        heads[i][m_dateKeyword].value( mx::sys::ISO8601DateTimeStrMJD( imageMJD[i], 1 ) );
        heads[i]["START " + m_dateKeyword].value( mx::sys::ISO8601DateTimeStrMJD( startMJD[i], 1 ) );
        heads[i]["END " + m_dateKeyword].value( mx::sys::ISO8601DateTimeStrMJD( endMJD[i], 1 ) );
        heads[i].append( "DELTA " + m_dateKeyword,
                         ( endMJD[i] - startMJD[i] ) * 86400,
                         "change in " + m_dateKeyword + " in seconds." );

        for( size_t j = 0; j < coaddKeywords.size(); ++j )
        {
            heads[i][coaddKeywords[j]].value( avgVals[i][j] );
            heads[i]["START " + coaddKeywords[j]].value( startVals[i][j] );
            heads[i]["END " + coaddKeywords[j]].value( endVals[i][j] );
            heads[i]["DELTA " + coaddKeywords[j]].value( endVals[i][j] - startVals[i][j] );
        }

        heads[i].append( "IMAGES COADDED", coaddFileNames[i].size(), "number of images included in this coadd" );
        for( size_t j = 0; j < coaddFileNames[i].size(); ++j )
        {
            // fits::fitsHistoryType t;
            // fits::fitsHeaderCard c("HISTORY", t , "coadded file: " + coaddFileNames[i][j]);
            heads[i].append( "HISTORY", fits::fitsHistoryType(), "coadded file: " + coaddFileNames[i][j] );
        }
    }

    t_coadd_end = sys::get_curr_time();

} // void HCIobservation<_realT,verboseT>::coaddImages()

template <typename _realT, class verboseT>
void HCIobservation<_realT, verboseT>::readMask()
{
    /*** Load the mask ***/
    if( m_maskFile != "" )
    {
        std::cerr << "creating mask cube\n";
        fitsFileT ff;
        ff.read( m_mask, m_maskFile );

        ///\todo here re-size mask if needed to match imSize
        if( m_mask.rows() > m_imSize || m_mask.cols() > m_imSize )
        {
            imageT tmask = m_mask.block( (int)( 0.5 * ( m_mask.rows() - 1 ) - 0.5 * ( m_imSize - 1 ) ),
                                         (int)( 0.5 * ( m_mask.rows() - 1 ) - 0.5 * ( m_imSize - 1 ) ),
                                         m_imSize,
                                         m_imSize );
            m_mask = tmask;
        }

        makeMaskCube();
    }
}

template <typename realT, class verboseT>
void HCIobservation<realT, verboseT>::makeMaskCube()
{
    if( m_mask.rows() != m_Nrows || m_mask.cols() != m_Ncols )
    {
        // clang-format off
        std::string message = "Mask is not the same size as images.\n";
                   message += "    Mask:   rows=" + std::to_string(m_mask.rows()) + "\n";
                   message += "            cols=" + std::to_string(m_mask.cols()) + "\n";
                   message += "    Images: rows=" + std::to_string(m_Nrows) + "\n";
                   message += "            cols=" + std::to_string(m_Ncols) + "\n";
        // clang-format on

        mxThrowException( err::invalidconfig, "HCIobservation<realT, verboseT>::makeMaskCube", message );
    }

    m_maskCube.resize( m_Nrows, m_Ncols, m_Nims );

    for( int i = 0; i < m_Nims; ++i )
    {
        m_maskCube.image( i ) = m_mask;
    }
}

template <typename _realT, class verboseT>
void HCIobservation<_realT, verboseT>::preProcess( eigenCube<realT> &ims )
{
    t_preproc_begin = sys::get_curr_time();

    // The mask is applied first, and then after each subsequent P.P. step.
    if( m_maskFile != "" && m_preProcess_mask )
    {
        std::cerr << "Masking . . .\n";
#pragma omp parallel for
        for( int i = 0; i < ims.planes(); ++i )
        {
            ims.image( i ) *= m_mask;
        }
        std::cerr << "done\n";
    }

    if( m_preProcess_subradprof )
    {
        std::cerr << "subtracting radial profile . . .\n";

#pragma omp parallel
        {
            imageT rp;

#pragma omp for
            for( int i = 0; i < ims.planes(); ++i )
            {
                Eigen::Map<Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic>> imRef( ims.image( i ) );
                radprofim( rp, imRef, true );
                zeroNaNs( imRef );
            }
        }

        std::cerr << "done\n";

        if( m_maskFile != "" && m_preProcess_mask )
        {
            std::cerr << "Masking . . .\n";
#pragma omp parallel for
            for( int i = 0; i < ims.planes(); ++i )
            {
                ims.image( i ) *= m_mask;
            }
            std::cerr << "done\n";
        }
    }

    if( m_preProcess_medianUSM_fwhm > 0 )
    {
        std::cerr << "applying median USM . . .\n";

        imageT medmask;
        medmask.resize( ims.rows(), ims.cols() );
        medmask.setConstant( 1 );

        for( int z = 0; z < (int)( 0.5 * m_preProcess_medianUSM_fwhm ); ++z )
        {
            medmask.row( z ) = 0;
            medmask.row( medmask.rows() - 1 - z ) = 0;
            medmask.col( z ) = 0;
            medmask.col( medmask.cols() - 1 - z ) = 0;
        }

#pragma omp parallel
        {
            imageT fim, im;
            fim.resize( ims.rows(), ims.cols() );
            im.resize( ims.rows(), ims.cols() );

#pragma omp for
            for( int i = 0; i < ims.planes(); ++i )
            {
                im = ims.image( i );
                medianSmooth( fim, im, m_preProcess_medianUSM_fwhm );
                ims.image( i ) = ( im - fim ) * medmask;
            }
        }

        if( m_maskFile != "" && m_preProcess_mask )
        {
            std::cerr << "Masking . . .\n";
#pragma omp parallel for
            for( int i = 0; i < ims.planes(); ++i )
            {
                ims.image( i ) *= m_mask;
            }
        }
        std::cerr << "done\n";
    }

    if( m_preProcess_gaussUSM_fwhm > 0 )
    {
        std::cerr << "applying Gauss USM . . .\n";
        t_gaussusm_begin = sys::get_curr_time();

#pragma omp parallel for
        for( int i = 0; i < ims.planes(); ++i )
        {
            imageT fim, im;
            im = ims.image( i );
            filterImage( fim,
                         im,
                         gaussKernel<eigenImage<_realT>, 2>( m_preProcess_gaussUSM_fwhm ),
                         0.5 * ( ims.cols() - 1 ) - m_preProcess_gaussUSM_fwhm * 4 );
            im = ( im - fim );
            ims.image( i ) = im;
        }

        if( m_maskFile != "" && m_preProcess_mask )
        {
            std::cerr << "Masking . . .\n";
            // clang-format off
            #pragma omp parallel for // clang-format on
            for( int i = 0; i < ims.planes(); ++i )
            {
                ims.image( i ) *= m_mask;
            }
        }
        t_gaussusm_end = sys::get_curr_time();
        std::cerr << "done\n";
    }

    if( m_preProcess_azUSM_azW && m_preProcess_azUSM_radW )
    {
        ipc::ompLoopWatcher<> status( ims.planes(), std::cerr );

        std::cerr << "applying azimuthal USM . . .\n";
        t_azusm_begin = sys::get_curr_time();

        azBoxKernel<eigenImage<realT>> azbK( m_preProcess_azUSM_radW,
                                             m_preProcess_azUSM_azW,
                                             m_preProcess_azUSM_maxAz );

        precalcKernel pcK( azbK, ims.rows(), ims.cols(), 0.5 * ( ims.rows() - 1 ), 0.5 * ( ims.cols() - 1 ) );

        // clang-format off
        #pragma omp parallel for // clang-format on
        for( int i = 0; i < ims.planes(); ++i )
        {
            imageT fim, im;
            im = ims.image( i );
            medianFilterImage( fim, im, pcK );

            im = ( im - fim );
            ims.image( i ) = im;
            status.incrementAndOutputStatus();
        }

        status.clearOutput();
        status.outputFinalStatus();
        std::cerr << '\n';

        if( m_maskFile != "" && m_preProcess_mask )
        {
            std::cerr << "Masking . . .\n";
#pragma omp parallel for
            for( int i = 0; i < ims.planes(); ++i )
            {
                ims.image( i ) *= m_mask;
            }
        }

        t_azusm_end = sys::get_curr_time();
        std::cerr << "done (" << t_azusm_end - t_azusm_begin << " sec)                                \n";
    }

    preProcess_meanSub( ims );

    preProcess_pixelTSNorm( ims );

    t_preproc_end = sys::get_curr_time();

} // void HCIobservation<_realT,verboseT>::preProcess()

template <typename _realT, class verboseT>
void HCIobservation<_realT, verboseT>::preProcess_meanSub( eigenCube<realT> &ims )
{
    if( m_preProcess_meanSubMethod == HCI::meanSubMethod::none )
    {
        return;
    }
    else if( m_preProcess_meanSubMethod != HCI::meanSubMethod::meanImage &&
             m_preProcess_meanSubMethod != HCI::meanSubMethod::medianImage )
    {
        std::string msg = "Mean subtraction by " + HCI::meanSubMethodStr( m_preProcess_meanSubMethod );
        msg += " can't be done in pre-processing. Only meanImage or medianImage can be used in pre.";
        mxThrowException( err::invalidconfig, "HCIobservation::preProcess_meanSub", msg );
    }

    imageT mean;

    if( m_preProcess_meanSubMethod == HCI::meanSubMethod::meanImage )
    {
        ims.mean( mean );
    }
    else if( m_preProcess_meanSubMethod == HCI::meanSubMethod::medianImage )
    {
        ims.median( mean );
    }

#pragma omp parallel for
    for( int n = 0; n < ims.planes(); ++n )
    {
        ims.image( n ) -= mean;

        if( m_maskFile != "" && m_preProcess_mask )
        {
            ims.image( n ) *= m_mask;
        }
    }
}

template <typename _realT, class verboseT>
void HCIobservation<_realT, verboseT>::preProcess_pixelTSNorm( eigenCube<realT> &ims )
{
    if( m_preProcess_pixelTSNormMethod == HCI::pixelTSNormMethod::none )
    {
        return;
    }

    if( m_preProcess_pixelTSNormMethod == HCI::pixelTSNormMethod::rmsSigmaClipped )
    {
        mxThrowException( err::notimpl,
                          "KlipReduction::preProcess_pixelTSNorm",
                          "pixelTSNormMethod is rmsSigmaClipped, which is not implemented" );
    }

    std::cerr << "normalizing pixels\n";

#pragma omp parallel
    {
        std::vector<realT> pixs( ims.planes() );

#pragma omp for
        for( int cc = 0; cc < ims.cols(); ++cc )
        {
            for( int rr = 0; rr < ims.rows(); ++rr )
            {
                if( m_maskFile != "" && m_preProcess_mask )
                {
                    if( m_mask( rr, cc ) == 0 )
                    {
                        continue;
                    }
                }

                // We bother to load a vector in prep to add sigma clipping later.
                for( int pp = 0; pp < ims.planes(); ++pp )
                {
                    pixs[pp] = ims.image( pp )( rr, cc );
                }

                realT sd = sqrt( math::vectorVariance( pixs, 0 ) );

                if( sd == 0 )
                    continue;

                for( int pp = 0; pp < ims.planes(); ++pp )
                {
                    ims.image( pp )( rr, cc ) /= sd;
                }
            }
        }
    }
}

template <typename _realT, class verboseT>
int HCIobservation<_realT, verboseT>::readWeights()
{
    std::ifstream fin;
    std::string str;

    if( m_weightFile == "" )
    {
        mxError( "HCIobservation::readWeights", MXE_PARAMNOTSET, "m_weightFile not set" );
        return -1;
    }

    // Read the weight file and load it into a map
    std::vector<std::string> wfileNames;
    std::vector<realT> imW;

    if( ioutils::readColumns( m_weightFile, wfileNames, imW ) != mx::error_t::noerror )
    {
        return -1;
    }

    if( imW.size() < m_fileList.size() )
    {
        mxError( "HCIobservation::readWeights", MXE_SIZEERR, "not enough weights specified" );
        return -1;
    }

    std::map<std::string, realT> weights;
    for( size_t i = 0; i < wfileNames.size(); ++i )
        weights[ioutils::pathFilename( wfileNames[i].c_str() )] = imW[i];

    m_comboWeights.resize( m_fileList.size() );

    realT wi;
    realT weightSum = 0;
    for( size_t i = 0; i < m_fileList.size(); ++i )
    {
        try
        {
            wi = weights.at( ioutils::pathFilename( m_fileList[i].c_str() ) );
        }
        catch( ... )
        {
            mxError( "HCIobservation::readWeights", MXE_NOTFOUND, "Weight for a file in m_fileList not found." );
            return -1;
        }
        m_comboWeights[i] = wi;
        weightSum += wi;
    }

    // Finally normalize the weights
    for( size_t i = 0; i < m_comboWeights.size(); ++i )
    {
        m_comboWeights[i] /= weightSum;
    }

    return 0;
} // int HCIobservation<_realT,verboseT>::readWeights()

template <typename _realT, class verboseT>
void HCIobservation<_realT, verboseT>::combineFinim()
{
    if( m_combineMethod == HCI::combineMethod::none )
    {
        return;
    }

    t_combo_begin = sys::get_curr_time();

    // Create and size temporary image for averaging
    imageT tfinim;

    m_finim.resize( m_psfsub[0].rows(), m_psfsub[0].cols(), m_psfsub.size() );

    // Now cycle through each set of psf subtractions
    for( size_t n = 0; n < m_psfsub.size(); ++n )
    {
        if( m_combineMethod == HCI::combineMethod::median )
        {
            m_psfsub[n].median( tfinim );
            m_finim.image( n ) = tfinim;
        }
        else if( m_combineMethod == HCI::combineMethod::mean )
        {
            if( m_comboWeights.size() == (size_t)m_Nims )
            {
                if( m_maskFile != "" )
                {
                    m_psfsub[n].mean( tfinim, m_comboWeights, m_maskCube, m_minGoodFract );
                }
                else
                {
                    m_psfsub[n].mean( tfinim, m_comboWeights );
                }
            }
            else
            {
                if( m_maskFile != "" )
                {
                    m_psfsub[n].mean( tfinim, m_maskCube, m_minGoodFract );
                }
                else
                {
                    m_psfsub[n].mean( tfinim );
                }
            }
            m_finim.image( n ) = tfinim;
        }
        else if( m_combineMethod == HCI::combineMethod::sigmaMean )
        {
            if( m_comboWeights.size() == (size_t)m_Nims )
            {
                if( m_maskFile != "" )
                {
                    m_psfsub[n].sigmaMean( tfinim, m_comboWeights, m_maskCube, m_sigmaThreshold, m_minGoodFract );
                }
                else
                {
                    m_psfsub[n].sigmaMean( tfinim, m_comboWeights, m_sigmaThreshold );
                }
            }
            else
            {
                if( m_maskFile != "" )
                {
                    m_psfsub[n].sigmaMean( tfinim, m_maskCube, m_sigmaThreshold, m_minGoodFract );
                }
                else
                {
                    m_psfsub[n].sigmaMean( tfinim, m_sigmaThreshold );
                }
            }
            m_finim.image( n ) = tfinim;
        }
    }

    t_combo_end = sys::get_curr_time();
} // void HCIobservation<_realT,verboseT>::combineFinim()

template <typename _realT, class verboseT>
void HCIobservation<_realT, verboseT>::outputPreProcessed()
{
    if( m_preProcess_outputPrefix == "" )
    {
        return;
    }

    std::string dir = ioutils::parentPath( m_preProcess_outputPrefix );
    ioutils::createDirectories( dir );

    std::string fname;

    fitsFileT ff;

    /** \todo Should add a HISTORY card here */
    char nstr[16];
    for( int i = 0; i < m_Nims; ++i )
    {
        snprintf( nstr, sizeof( nstr ), "%06d", i );
        fname = m_preProcess_outputPrefix + nstr + ".fits";

        fitsHeaderT fh = m_heads[i];
        stdFitsHeader( fh );
        ff.write( fname, m_tgtIms.image( i ).data(), m_Ncols, m_Nrows, 1, fh );
    }
} // void HCIobservation<_realT,verboseT>::outputPreProcessed()

template <typename _realT, class verboseT>
void HCIobservation<_realT, verboseT>::stdFitsHeader( fitsHeaderT &head )
{
    head.append( "", fits::fitsCommentType(), "----------------------------------------" );
    head.append( "", fits::fitsCommentType(), "mx::HCIobservation parameters:" );
    head.append( "", fits::fitsCommentType(), "----------------------------------------" );

    head.template append<int>( "FDELFRNT", m_deleteFront, "images deleted from front of file list" );
    head.template append<int>( "FDELBACK", m_deleteBack, "images deleted from back of file list" );

    head.append( "QFILE", ioutils::pathFilename( m_qualityFile.c_str() ), "quality file for thresholding" );
    head.template append<realT>( "QTHRESH", m_qualityThreshold, "quality threshold" );
    head.template append<int>( "NUMIMS", m_Nims, "number of images processed" );

    head.template append<int>( "IMSIZE", m_imSize, "image size after reading" );

    head.template append<std::string>( "COADMTHD", HCI::coaddMethodStr( m_coaddMethod ), "coadd combination method" );
    if( m_coaddMethod != HCI::coaddMethod::none )
    {
        head.template append<int>( "COADIMNO", m_coaddMaxImno, "max number of images in each coadd" );
        head.template append<realT>( "COADTIME", m_coaddMaxTime, "max time in each coadd" );
    }
    else
    {
        head.template append<int>( "COADIMNO", 0, "max number of images in each coadd" );
        head.template append<realT>( "COADTIME", 0, "max time in each coadd" );
    }

    head.append( "MASKFILE", m_maskFile, "mask file" );

    head.template append<int>( "PREPROC BEFORE", m_preProcess_beforeCoadd, "pre-process before coadd flag" );
    head.template append<int>( "PREPROC MASK", m_preProcess_mask, "pre-process mask flag" );
    head.template append<int>( "PREPROC SUBRADPROF",
                               m_preProcess_subradprof,
                               "pre-process subtract radial profile flag" );
    head.template append<realT>( "PREPROC AZUSM AZWIDTH",
                                 m_preProcess_azUSM_azW,
                                 "pre-process azimuthal USM azimuthal width [pixels]" );
    head.template append<realT>( "PREPROC AZUSM MAXAZ",
                                 m_preProcess_azUSM_maxAz,
                                 "pre-process azimuthal USM maximum azimuthal width [degrees]" );
    head.template append<realT>( "PREPROC AZUSM RADWIDTH",
                                 m_preProcess_azUSM_radW,
                                 "pre-process azimuthal USM radial width [pixels]" );
    head.template append<realT>( "PREPROC MEDIANUSM FWHM",
                                 m_preProcess_medianUSM_fwhm,
                                 "pre-process median USM fwhm [pixels]" );
    head.template append<realT>( "PREPROC GAUSSUSM FWHM",
                                 m_preProcess_gaussUSM_fwhm,
                                 "pre-process Gaussian USM fwhm [pixels]" );
    head.template append<std::string>( "PREPROC MEANSUB METHOD",
                                       HCI::meanSubMethodStr( m_preProcess_meanSubMethod ),
                                       "pre-process mean subtraction method" );
    head.template append<std::string>( "PREPROC PIXELTSNORM METHOD",
                                       HCI::pixelTSNormMethodStr( m_preProcess_pixelTSNormMethod ),
                                       "pre-process pixel time-series norm method" );
}

template <typename _realT, class verboseT>
void HCIobservation<_realT, verboseT>::writeFinim( fitsHeaderT *addHead )
{
    std::string fname = m_finimName;

    if( m_outputDir != "" )
    {
        mkdir( m_outputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );

        fname = m_outputDir + "/" + fname;
    }

    if( !m_exactFinimName )
    {
        fname = ioutils::getSequentialFilename( fname, ".fits" );
    }

    fitsHeaderT head;

    // Add HCIobservation standard header:
    stdFitsHeader( head );

    // Now add the final combination details:
    head.template append<std::string>( "COMBINATION METHOD",
                                       HCI::combineMethodStr( m_combineMethod ),
                                       "combination method" );

    if( m_weightFile != "" )
        head.append( "WEIGHT FILE", m_weightFile, "file containing weights for combination" );

    if( m_combineMethod == HCI::combineMethod::sigmaMean )
        head.template append<realT>( "SIGMA THRESHOLD", m_sigmaThreshold, "threshold for sigma clipping" );

    head.template append<realT>( "MIN FOOD FRACTION", m_minGoodFract, "minimum good fraction for inclusion" );
    if( addHead )
    {
        head.append( *addHead );
    }

    fits::fitsHeaderGitStatus( head, "mxlib_uncomp", MXLIB_UNCOMP_CURRENT_SHA1, MXLIB_UNCOMP_REPO_MODIFIED );

    fitsFileT f;

    f.write( fname, m_finim, head );

    std::cerr << "Final image written to: " << fname << "\n";
} // void HCIobservation<_realT,verboseT>::writeFinim(fitsHeaderT * addHead)

template <typename _realT, class verboseT>
void HCIobservation<_realT, verboseT>::outputPSFSub( fitsHeaderT *addHead )
{

    std::string fname;

    fitsHeaderT head;

    // Add the HCIobservation standard fits header
    stdFitsHeader( head );

    if( addHead )
    {
        head.append( *addHead );
    }

    fits::fitsHeaderGitStatus( head, "mxlib_uncomp", MXLIB_UNCOMP_CURRENT_SHA1, MXLIB_UNCOMP_REPO_MODIFIED );

    fitsFileT f;

    std::ofstream wout;

    if( m_comboWeights.size() > 0 )
    {
        fname = m_PSFSubPrefix + "weights.dat";
        if( m_outputDir != "" )
        {
            mkdir( m_outputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
            fname = m_outputDir + "/" + fname;
        }
        wout.open( fname );
        std::cerr << "writing comboWeights: " << fname << "\n";
    }

    char num[256];
    for( size_t n = 0; n < m_psfsub.size(); ++n )
    {
        for( int p = 0; p < m_psfsub[n].planes(); ++p )
        {
            snprintf( num, 256, "_%03zu_%05d.fits", n, p );
            fname = m_PSFSubPrefix + num;

            if( m_outputDir != "" )
            {
                mkdir( m_outputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
                fname = m_outputDir + "/" + fname;
            }

            fitsHeaderT h = head;

            h.append( m_heads[p] );

            f.write( fname, m_psfsub[n].image( p ).data(), m_psfsub[n].rows(), m_psfsub[n].cols(), 1, h );

            if( m_comboWeights.size() > 0 && n == 0 )
            {
                wout << fname << " " << m_comboWeights[p] << "\n";
            }
        }
    }

    if( m_comboWeights.size() > 0 )
    {
        wout.close();
    }
} // void HCIobservation<_realT,verboseT>::outputPSFSub(fitsHeaderT * addHead)

/*
template <typename _realT, class verboseT>
int HCIobservation<_realT,verboseT>::readPSFSub( const std::string &dir,
                                        const std::string &prefix,
                                        const std::string &ext,
                                        size_t nReductions )
{

    m_psfsub.resize( nReductions );

    // Load first file to condigure based on its header.
    std::vector<std::string> flist = ioutils::getFileNames( dir, prefix, "000", ext );
    fitsHeaderT fh;
    eigenImage<realT> im;
    fitsFileT ff;

    ff.read( im, fh, flist[0] );

    if( fh.count( "FDELFRNT" ) == 0 )
    {
        mxError( "KLIPReduction", MXE_PARAMNOTSET, "FDELFRNT not found in FITS header." );
        return -1;
    }
    m_deleteFront = fh["FDELFRNT"].template value<int>();
    std::cerr << "deleteFront: " << m_deleteFront << "\n";

    if( fh.count( "FDELBACK" ) == 0 )
    {
        mxError( "KLIPReduction", MXE_PARAMNOTSET, "FDELBACK not found in FITS header." );
        return -1;
    }
    m_deleteBack = fh["FDELBACK"].template value<int>();
    std::cerr << "deleteBack: " << m_deleteBack << "\n";

    if( fh.count( "QFILE" ) == 0 )
    {
        mxError( "KLIPReduction", MXE_PARAMNOTSET, "QFILE not found in FITS header." );
        return -1;
    }
    m_qualityFile = fh["QFILE"].String();
    std::cerr << "qualityFile: " << m_qualityFile << "\n";

    if( fh.count( "QTHRESH" ) == 0 )
    {
        mxError( "KLIPReduction", MXE_PARAMNOTSET, "QTHRESH not found in FITS header." );
        return -1;
    }
    m_qualityThreshold = fh["QTHRESH"].template value<realT>();
    std::cerr << "qualityThreshold: " << m_qualityThreshold << "\n";

    if( fh.count( "COADMTHD" ) == 0 )
    {
        mxError( "KLIPReduction", MXE_PARAMNOTSET, "COADMTHD not found in FITS header." );
        return -1;
    }
    m_coaddMethod = HCI::combineMethodFmStr( fh["COADMTHD"].String() );
    std::cerr << "coaddMethod: " << m_coaddMethod << "\n";

    if( fh.count( "COADIMNO" ) != 0 )
    {
        m_coaddMaxImno = fh["COADIMNO"].template value<int>();
        std::cerr << "coaddMaxImno: " << m_coaddMaxImno << "\n";
    }

    if( fh.count( "COADTIME" ) != 0 )
    {
        m_coaddMaxImno = fh["COADTIME"].template value<realT>();
        std::cerr << "coaddMaxtime: " << m_coaddMaxTime << "\n";
    }

    if( m_maskFile == "" )
    {
        if( fh.count( "MASKFILE" ) == 0 )
        {
            mxError( "KLIPReduction",
                     MXE_PARAMNOTSET,
                     "MASKFILE not found in FITS header and not set in configuration." );
            return -1;
        }
        m_maskFile = fh["MASKFILE"].String();
    }
    std::cerr << "maskFile: " << m_maskFile << "\n";

    if( fh.count( "PPBEFORE" ) == 0 )
    {
        mxError( "KLIPReduction", MXE_PARAMNOTSET, "PPBEFORE not found in FITS header." );
        return -1;
    }
    m_preProcess_beforeCoadd = fh["PPBEFORE"].template value<int>();
    std::cerr << "preProcess_beforeCoadd: " << m_preProcess_beforeCoadd << "\n";

    if( fh.count( "PPMASK" ) == 0 )
    {
        mxError( "KLIPReduction", MXE_PARAMNOTSET, "PPMASK not found in FITS header." );
        return -1;
    }
    m_preProcess_mask = fh["PPMASK"].template value<int>();
    std::cerr << "preProcess_mask: " << m_preProcess_mask << "\n";

    if( fh.count( "PPSUBRAD" ) == 0 )
    {
        mxError( "KLIPReduction", MXE_PARAMNOTSET, "PPSUBRAD not found in FITS header." );
        return -1;
    }
    m_preProcess_subradprof = fh["PPSUBRAD"].template value<int>();
    std::cerr << "preProcess_subradprof: " << m_preProcess_subradprof << "\n";

    if( fh.count( "PPAUSMAW" ) == 0 )
    {
        mxError( "KLIPReduction", MXE_PARAMNOTSET, "PPAUSMAW not found in FITS header." );
        return -1;
    }
    m_preProcess_azUSM_azW = fh["PPAUSMAW"].template value<realT>();
    std::cerr << "preProcess_azUSM_azW: " << m_preProcess_azUSM_azW << "\n";

    if( fh.count( "PPAUSMRW" ) == 0 )
    {
        mxError( "KLIPReduction", MXE_PARAMNOTSET, "PPAUSMRW not found in FITS header." );
        return -1;
    }
    m_preProcess_azUSM_radW = fh["PPAUSMRW"].template value<realT>();
    std::cerr << "preProcess_azUSM_radW: " << m_preProcess_azUSM_radW << "\n";

    if( fh.count( "PPGUSMFW" ) == 0 )
    {
        mxError( "KLIPReduction", MXE_PARAMNOTSET, "PPGUSMFW not found in FITS header." );
        return -1;
    }
    m_preProcess_gaussUSM_fwhm = fh["PPGUSMFW"].template value<realT>();
    std::cerr << "preProcess_gaussUSM_fwhm: " << m_preProcess_gaussUSM_fwhm << "\n";

    fitsHeaderT head;

    if( m_dateKeyword != "" )
        head.append( m_dateKeyword );

    for( size_t i = 0; i < m_keywords.size(); ++i )
    {
        head.append( m_keywords[i] );
    }

    for( size_t n = 0; n < nReductions; ++n )
    {
        char nstr[5];
        int nwr = snprintf( nstr, sizeof( nstr ), "%03zu", n );
        if( nwr < 0 || n >= sizeof( nstr ) )
        {
            std::cerr << "possibly bad formatting in filename\n";
        }

        std::string nprefix = prefix + "_" + nstr + "_";
        load_fileList( dir, nprefix, ext );

        if( m_fileList.size() == 0 )
        {
            mxError( "HCIobservation",
                     MXE_FILENOTFOUND,
                     "The m_fileList has 0 length, there are no files to be read." );
            return -1;
        }

        Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> im;

        fitsFileT f( m_fileList[0] );

        fitsHeaderT fh = head;
        f.read( im, fh );

        // We set imSize to match the first image, but we make it a square.
        if( m_imSize == 0 )
        {
            m_imSize = im.rows();
            if( m_imSize > im.cols() )
            {
                m_imSize = im.cols();
            }
        }
        else
        {
            // Now make sure we don't read too much.
            if( m_imSize > im.rows() )
            {
                m_imSize = im.rows();
            }
            if( m_imSize > im.cols() )
            {
                m_imSize = im.cols();
            }
        }

        // the +0.1 is just to make sure we don't have a problem with precision (we shouldn't)/
        f.setReadSize( floor( 0.5 * ( im.rows() - 1 ) - 0.5 * ( m_imSize - 1 ) + 0.1 ),
                       floor( 0.5 * ( im.cols() - 1.0 ) - 0.5 * ( m_imSize - 1.0 ) + 0.1 ),
                       m_imSize,
                       m_imSize );
        im.resize( m_imSize, m_imSize );

        std::cerr << "image size is " << m_imSize << "\n";

        if( n > 0 )
        {
            if( m_fileList.size() != (size_t)m_Nims )
            {
                mxError( "HCIobservation", MXE_INVALIDARG, "Different number of images in reductions." );
                return -1;
            }
            if( m_Nrows != im.rows() )
            {
                mxError( "HCIobservation", MXE_INVALIDARG, "Different number of rows in reductions." );
                return -1;
            }
            if( m_Ncols != im.cols() )
            {
                mxError( "HCIobservation", MXE_INVALIDARG, "Different number of cols in reductions." );
                return -1;
            }
        }
        else
        {
            std::cerr << "found " << nReductions << " sets of " << m_fileList.size() << " " << im.rows() << " x "
                      << im.cols() << " files\n";
        }
        m_Nims = m_fileList.size();
        m_Nrows = im.rows();
        m_Ncols = im.cols();
        m_Npix = im.rows() * im.cols();

        m_psfsub[n].resize( m_Nrows, m_Ncols, m_Nims );

        m_heads.clear(); // This is necessary to make sure heads.resize() copies head on a 2nd call
        m_heads.resize( m_fileList.size(), head );

        t_load_begin = sys::get_curr_time();

        f.read( m_psfsub[n].data(), m_heads, m_fileList );

        f.setReadSize();

        if( m_dateKeyword != "" )
        {
            m_imageMJD.resize( m_heads.size() );

            if( m_dateIsISO8601 )
            {
                for( size_t i = 0; i < m_imageMJD.size(); ++i )
                {
                    m_imageMJD[i] = sys::ISO8601date2mjd( m_heads[i][m_dateKeyword].String() );
                }
            }
            else
            {
                for( size_t i = 0; i < m_imageMJD.size(); ++i )
                {
                    m_imageMJD[i] = m_heads[i][m_dateKeyword].template value<realT>() * m_dateUnit;
                }
            }
        }

        t_load_end = sys::get_curr_time();

        for( size_t n = 0; n < m_psfsub.size(); ++n )
        {
            zeroNaNCube( m_psfsub[n] );
        }
    }

    if( m_weightFile != "" )
    {
        std::vector<std::string> fn;
        ioutils::readColumns( m_weightFile, fn, m_comboWeights );

        std::cerr << "read: " << m_weightFile << " (" << m_comboWeights.size() << ")\n";
    }

    /*** Now do the post-read actions ***
    if( postReadFiles() < 0 )
    {
        return -1;
    }

    readMask();

    m_filesRead = true;

    return 0;
}
*/
///@}

extern template class HCIobservation<float>;
extern template class HCIobservation<double>;

} // namespace improc
} // namespace mx

#endif //__HCIobservation_hpp__
