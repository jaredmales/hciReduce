/** \file HCIobservation.hpp
 * \author Jared R. Males
 * \brief Defines the basic high contrast imaging data type.
 * \ingroup hc_imaging_files
 * \ingroup image_processing_files
 *
 */

#ifndef __HCIobservation_hpp__
#define __HCIobservation_hpp__

#include <algorithm>
#include <charconv>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <mx/mxlib.hpp>

#include <mx/app/appConfigurator.hpp>

// #include <mx/mxException.hpp>

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

#include "HCI.hpp"

namespace mx
{

namespace improc
{

/// The basic high contrast imaging data type
/** This class manages file reading, resizing, co-adding, pre-processing (masking and filtering),
 * final image combination, and output.
 *
 * \tparam _realT is the floating point type in which to do all arithmetic.
 * \tparam verboseT sets the verbosity of error reporting
 *
 * \ingroup hc_imaging
 */
template <typename _realT, class verboseT>
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
     * - HCI::coadd::none -- [default] do not combine.  This turns off coadding.
     * - HCI::coadd::median --  coadded image is the median
     * - HCI::coadd::mean -- coadded image is the simple mean
     *
     * No other types of combination are currently supported for coadding.
     */
    HCI::coadd m_coaddMethod{ HCI::coadd::none };

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
     * -# symmetric median unsharp mask (m_preProcess_medianUSM_fwhm)
     * -# symmetric Gaussian unsharp mask (m_preProcess_gaussUSM_fwhm)
     * -# mask applied (enabled by m_preProcess_mask)
     * -# azimuthal unsharp mask (m_preProcess_azUSM_azHalfWidth and m_preProcess_azUSM_radHalfWidth)
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

    /// Azimuthal boxcar half-width for azimuthal unsharp mask [pixels]
    /** If this is 0 then azimuthal-USM is not performed.
     */
    realT m_preProcess_azUSM_azHalfWidth{ 0 };

    /// Maximum azimuthal boxcar width for azimuthal unsharp mask [degrees]
    /** Limits width close to center, preventing wrap-around.  Default is 45 degrees.  Set to 0 for no maximum.
     */
    realT m_preProcess_azUSM_maxAz{ 45 };

    /// Radial boxcar half-width for azimuthal unsharp mask [pixels]
    /** If this is 0 then azimuthal-USM is not performed.
     */
    realT m_preProcess_azUSM_radHalfWidth{ 0 };

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
    HCI::meanSub m_preProcess_meanSubMethod{ HCI::meanSub::none };

    /// Specify if each pixel time-series is normalized
    /** This normalizaton is applied after centering. Can have the following values:
     * - <b>HCI::pixelTSNorm::none</b>: no normalization (the default)
     * - <b>HCI::pixelTSNorm::rms</b>: divide by the time-series rms
     * - <b>HCI::pixelTSNorm::rmsSigmaClipped</b>: divide by the sigma-slipped time-series rms.
     *                                                   The sigma is provided by m_preProcess_pixelTSSigma.
     */
    HCI::pixelTSNorm m_preProcess_pixelTSNormMethod{ HCI::pixelTSNorm::none };

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

    /// Register this observation's configuration options.
    int setupConfig( mx::app::appConfigurator &config /**< [in.out] application configuration manager */ );

    /// Load this observation's state from registered configuration values.
    int loadConfig( mx::app::appConfigurator &config /**< [in] populated application configuration manager */ );

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

    /// Read and center-crop one configured mask file to the input image size.
    void readMaskFile( imageT &mask,                    /**< [out] loaded and cropped mask */
                       const std::string &maskFileName, /**< [in] FITS mask path */
                       const std::string &description /**< [in] mask role for diagnostics */ );

    /// Replicate a validated image mask into a cube.
    void populateMaskCube( eigenCube<realT> &maskCube, /**< [out] replicated mask cube */
                           const imageT &mask,         /**< [in] two-dimensional mask */
                           int rows,                   /**< [in] required row count */
                           int cols,                   /**< [in] required column count */
                           int planes /**< [in] output plane count */ );

    /// Select the target or RDI mask associated with an image cube.
    const imageT *preProcessMask( const eigenCube<realT> &ims /**< [in] cube being preprocessed */ ) const;

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
    void readFiles();

    /// Perform post-read actions for the target images, for use by derived classes
    virtual void postReadFiles();

    /// Perform post-coadd actions for the target images, for use by derived classes.
    /**
     * \returns 0 on success
     * \returns \<0 on error.
     */
    virtual void postCoadd();

    ///@}

    /// Read the list of reference files, cut to size, and preprocess.
    /** The target files must be read with \ref readFiles() before calling this method.
     *
     * \throws on  error.
     */
    void readRDIFiles();

    /// Perform post-read actions for the RDI images, for use by derived classes
    virtual void postRDIReadFiles();

    /// Perform post-coadd actions, for use by derived classes.
    /** A key example is to update keywords after the averaging occurs in coaddImages().
     *
     * \returns 0 on success
     * \returns \<0 on error.
     */
    virtual void postRDICoadd();

    ///@}
    //--RDI

    /** \name Thresholding
     * Thresholds are applied to a list of files before it is read, based on the image qualities supplied.
     * @{
     */

    /// Read the image qualities from a qualityFile and apply the threshold to a fileList
    /** This is called by readFiles().
     *
     * \throws on  error.
     */
    void threshold( std::vector<std::string> &fileList, ///< [in.out] the fileList to threshold
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
    void coaddImages( HCI::coadd coaddMethod,                        /**< [in] image combination method */
                      int coaddMaxImno,                              /**< [in] maximum images per output group */
                      realT coaddMaxTime,                            /**< [in] maximum group duration in seconds */
                      const std::vector<std::string> &coaddKeywords, /**< [in] numeric header cards to average */
                      const std::vector<std::string> &fileList,      /**< [in] filenames for HISTORY provenance */
                      const std::string &dateKeyword,                /**< [in] date card name, or empty for no dates */
                      std::vector<double> &imageMJD,                 /**< [in.out] image dates in MJD */
                      std::vector<fitsHeaderT> &heads,               /**< [in.out] per-image FITS headers */
                      eigenCube<realT> &ims /**< [in.out] image cube */ );

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
     * - HCI::combine::median -- [default] final image is the median
     * - HCI::combine::mean -- final image is the simple mean
     * - HCI::weightedMeanCombine -- final image is the weighted mean.  m_weightFile must be provided.
     * - HCI::combine::sigmaMean -- final image is sigma clipped mean.  If m_sigmaThreshold \<= 0, then it reverts
     * to meanCombine.
     */
    HCI::combine m_combineMethod{ HCI::combine::mean };

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

    /// The standard deviation threshold used if m_combineMethod == HCI::combine::sigmaMean.
    realT m_sigmaThreshold{ 0 };

    /// The minimum fraction of good (un-masked) pixels to include in the final combination (0.0 to 1.0). If not met,
    /// then the pixel will be NaN-ed.
    realT m_minGoodFract{ 0.0 };

    /// Read the image weights from m_weightFile
    /** This is called by readFiles().
     *
     * \throws on error
     */
    void readWeights();

    /// Combine the images into a single final image.
    /** Images are combined by the method specified in \ref m_combineMethod
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

    /// Output the pre-processed target images using the configured output prefix and six-digit sequence numbers.
    /** \throws mx::exception if the cube/header sizes are invalid or an output operation fails. */
    void outputPreProcessed();

    /// Output the pre-processed reference images using an `RDI_` suffix after the configured output prefix.
    /** \throws mx::exception if the cube/header sizes are invalid or an output operation fails. */
    void outputRDIPreProcessed();

    /// Append the standard HCIobservation reduction metadata to a FITS header.
    void stdFitsHeader( fitsHeaderT &head /**< [in.out] the fitsHeader structure which will
                                                             have cards appended to it. */
    );

    /// Write the final combined image cube to its exact or next sequential output path.
    /** \throws mx::exception if the cube is empty or an output operation fails. */
    void writeFinim( fitsHeaderT *addHead = 0 );

    /// Write every PSF-subtracted reduction/image and the optional combination-weight sidecar.
    /** \throws mx::exception if the input sizes are inconsistent or an output operation fails. */
    void outputPSFSub( fitsHeaderT *addHead = 0 );

    ///@}

    /// Read a complete saved collection of PSF-subtracted FITS images.
    /** Files must use the naming and header contract produced by outputPSFSub(). The reduction and image counts are
     * inferred from a complete, zero-based rectangular index grid.
     *
     * \throws mx::exception if the directory, names, headers, index grid, or image dimensions are invalid.
     */
    void readPSFSub( const std::string &directory, /**< [in] directory containing the saved products */
                     const std::string &prefix,    /**< [in] literal filename prefix used by outputPSFSub() */
                     const std::string &extension = ".fits" /**< [in] literal filename extension */ );

    double t_begin{ 0 }; ///< Start time of the complete reduction.
    double t_end{ 0 };   ///< End time of the complete reduction.

    double t_load_begin{ 0 }; ///< Start time of input loading.
    double t_load_end{ 0 };   ///< End time of input loading.

    double t_coadd_begin{ 0 }; ///< Start time of image coaddition.
    double t_coadd_end{ 0 };   ///< End time of image coaddition.

    double t_preproc_begin{ 0 }; ///< Start time of preprocessing.
    double t_preproc_end{ 0 };   ///< End time of preprocessing.

    double t_azusm_begin{ 0 }; ///< Start time of azimuthal unsharp masking.
    double t_azusm_end{ 0 };   ///< End time of azimuthal unsharp masking.

    double t_gaussusm_begin{ 0 }; ///< Start time of Gaussian unsharp masking.
    double t_gaussusm_end{ 0 };   ///< End time of Gaussian unsharp masking.

    double t_combo_begin{ 0 }; ///< Start time of final image combination.
    double t_combo_end{ 0 };   ///< End time of final image combination.
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
                mx::app::argType::Required,
                "input",
                "dateUnit",
                false,
                "float",
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
                mx::app::argType::Required,
                "rdi",
                "dateUnit",
                false,
                "float",
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
                mx::app::argType::True,
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
                "bool",
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

    config.add( "preProcess.azUSM_azHalfWidth",
                "",
                "preProcess.azUSM_azHalfWidth",
                mx::app::argType::Required,
                "preProcess",
                "azUSM_azHalfWidth",
                false,
                "float",
                "The azimuth USM boxcar azimuthal half-width in pixels. Enabled when both half-widths are non-zero." );

    // Deprecated compatibility alias. The explicitly named half-width key wins when both are set.
    config.add( "preProcess.azUSM_azW",
                "",
                "preProcess.azUSM_azW",
                mx::app::argType::Required,
                "preProcess",
                "azUSM_azW",
                false,
                "float",
                "Deprecated alias for preProcess.azUSM_azHalfWidth." );

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

    config.add( "preProcess.azUSM_radHalfWidth",
                "",
                "preProcess.azUSM_radHalfWidth",
                mx::app::argType::Required,
                "preProcess",
                "azUSM_radHalfWidth",
                false,
                "float",
                "The azimuth USM boxcar radial half-width in pixels. Enabled when both half-widths are non-zero." );

    config.add( "preProcess.azUSM_radW",
                "",
                "preProcess.azUSM_radW",
                mx::app::argType::Required,
                "preProcess",
                "azUSM_radW",
                false,
                "float",
                "Deprecated alias for preProcess.azUSM_radHalfWidth." );

    config.add( "preProcess.medianUSM_fwhm",
                "",
                "preProcess.medianUSM_fwhm",
                mx::app::argType::Required,
                "preProcess",
                "medianUSM_fwhm",
                false,
                "int",
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
                "Options are none, rms, rmsSigmaClipped." );

    config.add( "preProcess.pixelTSSigma",
                "",
                "preProcess.pixelTSSigma",
                mx::app::argType::Required,
                "preProcess",
                "pixelTSSigma",
                false,
                "float",
                "Sigma-clipping threshold for pixel time-series RMS normalization." );

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

    config.add( "combine.method",
                "",
                "combine.method",
                mx::app::argType::Required,
                "combine",
                "method",
                false,
                "string",
                "Averaging method for final combination: mean, median, weighted, sigma" );

    config.add( "combine.weightFile",
                "",
                "combine.weightFile",
                mx::app::argType::Required,
                "combine",
                "weightFile",
                false,
                "string",
                "File containing weights for the weighted combo.  Two column format: filename weight" );

    config.add( "combine.sigmaThreshold",
                "",
                "combine.sigmaThreshold",
                mx::app::argType::Required,
                "combine",
                "sigmaThreshold",
                false,
                "float",
                "Clipping threshold for sigma clipped mean combination." );

    config.add( "combine.minGoodFract",
                "",
                "combine.minGoodFract",
                mx::app::argType::Required,
                "combine",
                "minGoodFract",
                false,
                "float",
                "Minimum fraction of good/un-masked pixels to include in final image, otherwise pixel is NaN-ed." );

    config.add( "output.fileName",
                "",
                "output.fileName",
                mx::app::argType::Required,
                "output",
                "fileName",
                false,
                "string",
                "Prefix for output file name.  A 4 digit 0-padded number is appended." );

    config.add( "output.exactFName",
                "",
                "output.exactFName",
                mx::app::argType::True,
                "output",
                "exactFName",
                false,
                "bool",
                "Used outputFile exactly as specified, without appending a number or .fits" );

    config.add( "output.directory",
                "",
                "output.directory",
                mx::app::argType::Required,
                "output",
                "directory",
                false,
                "string",
                "The directory where to output files." );

    config.add( "output.outputPSFSub",
                "",
                "output.outputPSFSub",
                mx::app::argType::True,
                "output",
                "outputPSFSub",
                false,
                "bool",
                "Output the PSF subtracted images (default false)" );

    config.add( "output.psfSubPrefix",
                "",
                "output.psfSubPrefix",
                mx::app::argType::Required,
                "output",
                "psfSubPrefix",
                false,
                "string",
                "Prefix of the PSF subtracted output files." );

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

    std::string coaddToStr = HCI::coaddToStr<verboseT>( m_coaddMethod ); // get default
    config( coaddToStr, "coadd.method" );
    try
    {
        m_coaddMethod = HCI::coaddFmStr<verboseT>( coaddToStr );
    }
    catch( ... )
    {
        std::throw_with_nested( mx::exception<verboseT>( mx::error_t::invalidconfig, "invalid coadd method" ) );
    }

    config( m_coaddMaxImno, "coadd.maxImno" );
    config( m_coaddMaxTime, "coadd.maxTime" );
    config( m_coaddKeywords, "coadd.keywords" );

    config( m_preProcess_beforeCoadd, "preProcess.beforeCoadd" );
    config( m_preProcess_mask, "preProcess.mask" );
    config( m_preProcess_subradprof, "preProcess.subradprof" );
    config( m_preProcess_azUSM_azHalfWidth, "preProcess.azUSM_azW" );
    config( m_preProcess_azUSM_azHalfWidth, "preProcess.azUSM_azHalfWidth" );
    config( m_preProcess_azUSM_maxAz, "preProcess.azUSM_maxAz" );
    config( m_preProcess_azUSM_radHalfWidth, "preProcess.azUSM_radW" );
    config( m_preProcess_azUSM_radHalfWidth, "preProcess.azUSM_radHalfWidth" );
    config( m_preProcess_medianUSM_fwhm, "preProcess.medianUSM_fwhm" );
    config( m_preProcess_gaussUSM_fwhm, "preProcess.gaussUSM_fwhm" );

    std::string ppmsm;
    try
    {
        ppmsm = HCI::meanSubToStr<verboseT>( m_preProcess_meanSubMethod );
    }
    catch( ... )
    {
        std::throw_with_nested( mx::exception<verboseT>( mx::error_t::exception ) );
    }

    config( ppmsm, "preProcess.meanSubMethod" );
    try
    {
        m_preProcess_meanSubMethod = HCI::meanSubFmStr<verboseT>( ppmsm );

        if( m_preProcess_meanSubMethod != HCI::meanSub::none && m_preProcess_meanSubMethod != HCI::meanSub::meanImage &&
            m_preProcess_meanSubMethod != HCI::meanSub::medianImage )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           std::format( "Mean subtraction by {} "
                                                        "can't be done in pre-processing. "
                                                        "Only meanImage or medianImage can be used in pre.",
                                                        HCI::meanSubToStr<verboseT>( m_preProcess_meanSubMethod ) ) );
        }
    }
    catch( ... )
    {
        std::throw_with_nested( mx::exception<verboseT>( mx::error_t::invalidconfig,
                                                         "invalid pre-processing mean "
                                                         "subtraction method" ) );
    }

    std::string ptsnm;
    try
    {
        ptsnm = HCI::pixelTSNormToStr<verboseT>( m_preProcess_pixelTSNormMethod );
    }
    catch( ... )
    {
        std::throw_with_nested( mx::exception<verboseT>( mx::error_t::exception ) );
    }

    config( ptsnm, "preProcess.pixelTSNormMethod" );
    config( m_pixelTSSigma, "preProcess.pixelTSSigma" );

    try
    {
        m_preProcess_pixelTSNormMethod = HCI::pixelTSNormFmStr<verboseT>( ptsnm );
    }
    catch( ... )
    {
        std::throw_with_nested( mx::exception<verboseT>( mx::error_t::invalidconfig,
                                                         "invalid pixel time-series "
                                                         "normalization method" ) );
    }

    config( m_preProcess_outputPrefix, "preProcess.outputPrefix" );

    config( m_preProcess_only, "preProcess.only" );
    config( m_skipPreProcess, "preProcess.skip" );

    std::string cmbm;
    try
    {
        cmbm = HCI::combineToStr<verboseT>( m_combineMethod );
    }
    catch( ... )
    {
        std::throw_with_nested( mx::exception<verboseT>( mx::error_t::exception ) );
    }

    config( cmbm, "combine.method" );

    try
    {
        m_combineMethod = HCI::combineFmStr<verboseT>( cmbm );
    }
    catch( ... )
    {
        std::throw_with_nested( mx::exception<verboseT>( mx::error_t::invalidconfig, "invalid combination method" ) );
    }

    config( m_weightFile, "combine.weightFile" );

    config( m_sigmaThreshold, "combine.sigmaThreshold" );
    config( m_minGoodFract, "combine.minGoodFract" );

    config( m_finimName, "output.fileName" );
    config( m_exactFinimName, "output.exactFName" );
    config( m_outputDir, "output.directory" );
    config( m_doOutputPSFSub, "output.outputPSFSub" );
    config( m_PSFSubPrefix, "output.psfSubPrefix" );

    return 0;
}

template <typename _realT, class verboseT>
mx::error_t HCIobservation<_realT, verboseT>::load_fileList( std::vector<std::string> &fileList,
                                                             const std::string &fileListFile,
                                                             const std::string &directory,
                                                             const std::string &prefix,
                                                             const std::string &extension )
{
    fileList.clear();

    if( fileListFile != "" )
    {
        mx::error_t errc = ioutils::readColumns<mx::ioutils::readColSpaceDelim, verboseT>( fileListFile, fileList );
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

            for( auto &&fpath : fileList )
            {
                fpath = dir + fpath;
            }
        }
    }
    else if( directory != "" )
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
    m_filesRead = false;

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
    m_RDIfilesRead = false;

    if( errc != mx::error_t::noerror )
    {
        return mx::error_report<verboseT>( errc, "error loading file list" );
    }
    return mx::error_t::noerror;
}

// --< construction and initialization

template <typename realT, class verboseT>
void HCIobservation<realT, verboseT>::readFiles()
{
    m_filesRead = false;

    if( m_fileList.size() == 0 )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                       "The target fileList has 0 length, there are no files to be read." );
    }

    // First make the list deletions
    if( !m_filesDeleted )
    {
        if( m_deleteFront < 0 || m_deleteBack < 0 || static_cast<size_t>( m_deleteFront ) > m_fileList.size() ||
            static_cast<size_t>( m_deleteBack ) > m_fileList.size() - static_cast<size_t>( m_deleteFront ) )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "Target front/back deletions exceed the file-list bounds." );
        }

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
        throw mx::exception<verboseT>(
            mx::error_t::invalidconfig,
            "The target fileList has 0 length, there are no files to be read after deletions." );
    }

    if( m_qualityFile != "" )
    {
        std::cerr << "Thresholding target images...";
        size_t origsize = m_fileList.size();

        try
        {
            threshold( m_fileList, m_qualityFile, m_qualityThreshold );
        }
        catch( const std::exception &e )
        {
            std::throw_with_nested( mx::exception<verboseT>( mx::error_t::std_exception, "from threshold" ) );
        }

        if( m_fileList.size() == 0 )
        {
            throw mx::exception<verboseT>(
                mx::error_t::invalidconfig,
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
        try
        {
            readWeights();
        }
        catch( const std::exception &e )
        {
            std::throw_with_nested( mx::exception<verboseT>( mx::error_t::std_exception, "reading weights" ) );
        }
    }

    /*----- Append the HCI keywords to propagate them if needed -----*/

    fitsHeaderT head;

    if( m_dateKeyword != "" )
    {
        head.append( m_dateKeyword );
    }

    for( size_t i = 0; i < m_keywords.size(); ++i )
    {
        if( head.count( m_keywords[i] ) == 0 )
        {
            head.append( m_keywords[i] );
        }
    }
    m_heads.clear(); // This is necessary to make sure heads.resize() copies head on a 2nd call
    m_heads.resize( m_fileList.size(), head );

    /*----- Read in first image to test size -----*/

    Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> im;

    fitsFileT f( m_fileList[0] );

    mx::error_t errc = f.read( im );
    if( errc != mx::error_t::noerror )
    {
        throw mx::exception<verboseT>( errc, "Could not read the first target image." );
    }

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

    errc = f.read( m_tgtIms.data(), m_heads, m_fileList );
    if( errc != mx::error_t::noerror )
    {
        throw mx::exception<verboseT>( errc, "Could not read the target image list." );
    }

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
    else
    {
        m_imageMJD.clear();
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
    try
    {
        postReadFiles();
    }
    catch( ... )
    {
        std::throw_with_nested( mx::exception<verboseT>( mx::error_t::std_exception, "from postReadFiles" ) );
    }

    /*** Read in the mask if present ***/
    try
    {
        readMask();
    }
    catch( const std::exception &e )
    {
        std::throw_with_nested( mx::exception<verboseT>( mx::error_t::std_exception, "from readMask" ) );
    }

    /*** Now begin processing ***/
    if( !m_skipPreProcess )
    {
        /*** Now do any pre-processing ***/
        if( m_preProcess_beforeCoadd )
        {
            try
            {
                preProcess( m_tgtIms );
            }
            catch( const std::exception &e )
            {
                std::throw_with_nested( mx::exception<verboseT>( mx::error_t::std_exception, "from preProcess" ) );
            }
        }

        if( m_coaddMethod != HCI::coadd::none )
        {
            std::cerr << "Coadding target images...\n";
            try
            {
                coaddImages( m_coaddMethod,
                             m_coaddMaxImno,
                             m_coaddMaxTime,
                             m_coaddKeywords,
                             m_fileList,
                             m_dateKeyword,
                             m_imageMJD,
                             m_heads,
                             m_tgtIms );
            }
            catch( const std::exception &e )
            {
                std::throw_with_nested( mx::exception<verboseT>( mx::error_t::std_exception, "from coaddImages" ) );
            }

            m_Nims = m_tgtIms.planes();
            m_Nrows = m_tgtIms.rows();
            m_Ncols = m_tgtIms.cols();
            m_Npix = m_tgtIms.rows() * m_tgtIms.cols();

            std::cerr << "number of target images after coadding: " << m_Nims << "\n";

            try
            {
                postCoadd();
            }
            catch( ... )
            {
                std::throw_with_nested( mx::exception<verboseT>( mx::error_t::std_exception, "from postCoadd" ) );
            }
            std::cerr << "Done.\n";

            // Re-make the mask cube if we coadded...
            if( m_maskFile != "" )
            {
                try
                {
                    makeMaskCube();
                }
                catch( const std::exception &e )
                {
                    std::throw_with_nested(
                        mx::exception<verboseT>( mx::error_t::std_exception, "from makeMaskCube" ) );
                }
            }
        }

        /*** Now do any pre-processing if not done already***/
        if( !m_preProcess_beforeCoadd )
        {
            try
            {
                preProcess( m_tgtIms );
            }
            catch( const std::exception &e )
            {
                std::throw_with_nested( mx::exception<verboseT>( mx::error_t::std_exception, "from preProcess" ) );
            }
        }

        try
        {
            outputPreProcessed();
        }
        catch( const std::exception &e )
        {
            std::throw_with_nested( mx::exception<verboseT>( mx::error_t::std_exception, "from outputPreProcessed" ) );
        }
    }

    m_filesRead = true;

} // readFiles()

template <typename _realT, class verboseT>
void HCIobservation<_realT, verboseT>::postReadFiles()
{
}

template <typename _realT, class verboseT>
void HCIobservation<_realT, verboseT>::postCoadd()
{
}

//------------------- readRDIFiles
template <typename _realT, class verboseT>
void HCIobservation<_realT, verboseT>::readRDIFiles()
{

    m_RDIfilesRead = false;

    /* First check if the target files have been read */
    if( m_Nrows == 0 || m_Ncols == 0 )
    {
        throw mx::exception<verboseT>( mx::error_t::paramnotset,
                                       "The target image size must be set before reading RDI files." );
    }

    if( m_RDIfileList.size() == 0 )
    {
        throw mx::exception<verboseT>( mx::error_t::filenotfound,
                                       "The RDI fileList has 0 length, there are no files to be read." );
    }

    // First make the list deletions
    if( !m_RDIfilesDeleted )
    {
        if( m_RDIdeleteFront < 0 || m_RDIdeleteBack < 0 ||
            static_cast<size_t>( m_RDIdeleteFront ) > m_RDIfileList.size() ||
            static_cast<size_t>( m_RDIdeleteBack ) > m_RDIfileList.size() - static_cast<size_t>( m_RDIdeleteFront ) )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "RDI front/back deletions exceed the file-list bounds." );
        }

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
        throw mx::exception<verboseT>( mx::error_t::filenotfound,
                                       "The RDI fileList has 0 length, there are no "
                                       "files to be read after deletions." );
    }

    if( m_RDIqualityFile != "" )
    {
        std::cerr << "Thresholding RDI images...";
        size_t origsize = m_RDIfileList.size();

        try
        {
            threshold( m_RDIfileList, m_RDIqualityFile, m_RDIqualityThreshold );
        }
        catch( const std::exception &e )
        {
            std::throw_with_nested( mx::exception<verboseT>( mx::error_t::std_exception, "while thresholding" ) );
        }

        if( m_RDIfileList.size() == 0 )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "The fileList has 0 length, there are no "
                                           "files to be read after thresholding." );
        }

        std::cerr << "Done.  Selected " << m_RDIfileList.size() << " out of " << origsize << "\n";
    }

    /*----- Append the HCI keywords to propagate them if needed -----*/

    fitsHeaderT head;

    if( m_RDIdateKeyword != "" )
    {
        head.append( m_RDIdateKeyword );
    }

    for( size_t i = 0; i < m_RDIkeywords.size(); ++i )
    {
        head.append( m_RDIkeywords[i] );
    }

    m_RDIheads.clear(); // This is necessary to make sure heads.resize() copies head on a 2nd call
    m_RDIheads.resize( m_RDIfileList.size(), head );

    /*----- Read in first image to adjust size ----*/
    Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> im;

    fitsFileT f( m_RDIfileList[0] );

    mx::error_t errc = f.read( im );
    if( errc != mx::error_t::noerror )
    {
        throw mx::exception<verboseT>( errc, "Could not read the first RDI image." );
    }

    if( im.rows() < m_imSize || im.cols() < m_imSize )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr,
                                       "The reference images are too small, do not match the target images." );
    }

    // And now set the read size so we only read what we want.
    // the +0.1 is just to make sure we don't have a problem with precision (we shouldn't)/
    f.setReadSize( floor( 0.5 * ( im.rows() - 1 ) - 0.5 * ( m_imSize - 1 ) + 0.1 ),
                   floor( 0.5 * ( im.cols() - 1.0 ) - 0.5 * ( m_imSize - 1.0 ) + 0.1 ),
                   m_imSize,
                   m_imSize );

    m_refIms.resize( m_imSize, m_imSize, m_RDIfileList.size() );

    t_load_begin = sys::get_curr_time();

    errc = f.read( m_refIms.data(), m_RDIheads, m_RDIfileList );
    if( errc != mx::error_t::noerror )
    {
        throw mx::exception<verboseT>( errc, "Could not read the RDI image list." );
    }

    f.setReadSize();

    if( m_RDIdateKeyword != "" )
    {
        m_RDIimageMJD.resize( m_RDIheads.size() );

        if( m_RDIdateIsISO8601 )
        {
            for( size_t i = 0; i < m_RDIimageMJD.size(); ++i )
            {
                m_RDIimageMJD[i] = sys::ISO8601date2mjd( m_RDIheads[i][m_RDIdateKeyword].String() );
            }
        }
        else
        {
            for( size_t i = 0; i < m_RDIimageMJD.size(); ++i )
            {
                m_RDIimageMJD[i] = m_RDIheads[i][m_RDIdateKeyword].template value<realT>() * m_RDIdateUnit;
            }
        }
    }
    else
    {
        m_RDIimageMJD.clear();
    }

    t_load_end = sys::get_curr_time();

    std::cerr << "loading complete\n";

    std::cerr << "zero-ing NaNs\n";
    zeroNaNCube( m_refIms );
    std::cerr << "Done.\n";

    /*** Now do the post-read actions ***/
    try
    {
        postRDIReadFiles();
    }
    catch( ... )
    {
        std::throw_with_nested( mx::exception<verboseT>( mx::error_t::exception, "from postRDIReadFiles" ) );
    }

    /*** Read the independently configured RDI mask, if present. ***/
    try
    {
        if( m_RDImaskUseInput )
        {
            if( m_maskFile.empty() || m_mask.size() == 0 )
            {
                throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                               "The input mask must be loaded before it can be reused for RDI." );
            }
            m_RDImask = m_mask;
            populateMaskCube( m_RDImaskCube, m_RDImask, m_refIms.rows(), m_refIms.cols(), m_refIms.planes() );
        }
        else if( !m_RDImaskFile.empty() )
        {
            readMaskFile( m_RDImask, m_RDImaskFile, "RDI" );
            populateMaskCube( m_RDImaskCube, m_RDImask, m_refIms.rows(), m_refIms.cols(), m_refIms.planes() );
        }
        else
        {
            m_RDImask.resize( 0, 0 );
            m_RDImaskCube.clear();
        }
    }
    catch( ... )
    {
        std::throw_with_nested( mx::exception<verboseT>( mx::error_t::exception, "from RDI mask setup" ) );
    }

    /*** Now begin processing ***/
    if( !m_skipPreProcess )
    {
        /*** Now do any pre-processing ***/
        if( m_preProcess_beforeCoadd )
        {
            try
            {
                preProcess( m_refIms );
            }
            catch( ... )
            {
                std::throw_with_nested( mx::exception<verboseT>( mx::error_t::exception, "from preProcess" ) );
            }
        }

        if( m_coaddMethod != HCI::coadd::none )
        {
            std::cerr << "Coadding reference images...\n";
            try
            {
                coaddImages( m_coaddMethod,
                             m_coaddMaxImno,
                             m_coaddMaxTime,
                             m_coaddKeywords,
                             m_RDIfileList,
                             m_RDIdateKeyword,
                             m_RDIimageMJD,
                             m_RDIheads,
                             m_refIms );
            }
            catch( ... )
            {
                std::throw_with_nested( mx::exception<verboseT>( mx::error_t::exception, "from coaddImages" ) );
            }

            std::cerr << "number of reference images after coadding: " << m_refIms.planes() << "\n";

            try
            {
                postRDICoadd();
            }
            catch( ... )
            {
                std::throw_with_nested( mx::exception<verboseT>( mx::error_t::exception, "from postRDICoadd" ) );
            }
            std::cerr << "Done.\n";

            if( m_RDImask.size() != 0 )
            {
                populateMaskCube( m_RDImaskCube, m_RDImask, m_refIms.rows(), m_refIms.cols(), m_refIms.planes() );
            }
        }

        /*** Now do any pre-processing if not done already***/
        if( !m_preProcess_beforeCoadd )
        {
            try
            {
                preProcess( m_refIms );
            }
            catch( ... )
            {
                std::throw_with_nested( mx::exception<verboseT>( mx::error_t::exception, "from preProcess" ) );
            }
        }

        // outputRDIPreProcessed();
    }

    m_RDIfilesRead = true;

} // readRDIFiles()

template <typename _realT, class verboseT>
void HCIobservation<_realT, verboseT>::postRDIReadFiles()
{
}

template <typename _realT, class verboseT>
void HCIobservation<_realT, verboseT>::postRDICoadd()
{
}

template <typename _realT, class verboseT>
void HCIobservation<_realT, verboseT>::threshold( std::vector<std::string> &fileList,
                                                  const std::string &qualityFile,
                                                  realT qualityThreshold )
{
    if( qualityThreshold <= 0 )
    {
        return;
    }

    if( qualityFile == "" )
    {
        throw mx::exception<verboseT>( mx::error_t::paramnotset, "qualityFile not set" );
    }

    std::vector<std::string> qfileNames;
    std::vector<realT> imQ;

    // Read the quality file and load it into a map
    mx::error_t errc = ioutils::readColumns<mx::ioutils::readColSpaceDelim, verboseT>( qualityFile, qfileNames, imQ );
    if( errc != mx::error_t::noerror )
    {
        throw mx::exception<verboseT>( errc, "error reading quality file " + qualityFile );
    }

    if( qfileNames.size() != imQ.size() )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "quality file columns have different lengths" );
    }

    std::map<std::string, realT> quality;
    for( size_t i = 0; i < qfileNames.size(); ++i )
    {
        std::string fname = ioutils::pathFilename( qfileNames[i].c_str() );
        if( !std::isfinite( imQ[i] ) )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidarg, "non-finite quality for filename " + fname );
        }
        if( !quality.emplace( fname, imQ[i] ).second )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidarg, "duplicate quality for filename " + fname );
        }
    }

    realT q;

    for( size_t i = 0; i < fileList.size(); ++i )
    {
        const std::string fname = ioutils::pathFilename( fileList[i].c_str() );

        try
        {
            q = quality.at( fname );
        }
        catch( const std::exception &e )
        {
            std::throw_with_nested(
                mx::exception<verboseT>( mx::error_t::std_exception,
                                         std::format( "getting quality for filename {}", fname ) ) );
        }

        if( q < qualityThreshold )
        {
            fileList.erase( fileList.begin() + i );
            --i;
        }
    }

    // return 0;
}

template <typename _realT, class verboseT>
void HCIobservation<_realT, verboseT>::coaddImages( HCI::coadd coaddMethod,
                                                    int coaddMaxImno,
                                                    realT coaddMaxTime,
                                                    const std::vector<std::string> &coaddKeywords,
                                                    const std::vector<std::string> &fileList,
                                                    const std::string &dateKeyword,
                                                    std::vector<double> &imageMJD,
                                                    std::vector<fitsHeaderT> &heads,
                                                    eigenCube<realT> &ims )
{
    if( coaddMethod == HCI::coadd::none || ( coaddMaxImno <= 0 && coaddMaxTime <= 0 ) )
    {
        return;
    }

    if( coaddMethod != HCI::coadd::mean && coaddMethod != HCI::coadd::median )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig, "invalid image coadd method" );
    }

    const int inputImages = ims.planes();
    const int inputRows = ims.rows();
    const int inputCols = ims.cols();
    if( inputImages <= 0 )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "cannot coadd an empty image cube" );
    }
    if( fileList.size() != static_cast<size_t>( inputImages ) || heads.size() != static_cast<size_t>( inputImages ) )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "coadd filenames and headers must match the image count" );
    }

    const bool haveDates = !imageMJD.empty();
    if( haveDates && imageMJD.size() != static_cast<size_t>( inputImages ) )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "coadd date count must match the image count" );
    }
    if( coaddMaxTime > 0 && !haveDates )
    {
        throw mx::exception<verboseT>( mx::error_t::paramnotset, "time-limited coadding requires image dates" );
    }
    if( haveDates )
    {
        for( size_t index = 0; index < imageMJD.size(); ++index )
        {
            if( !std::isfinite( imageMJD[index] ) || ( index > 0 && imageMJD[index] < imageMJD[index - 1] ) )
            {
                throw mx::exception<verboseT>( mx::error_t::invalidarg,
                                               "coadd dates must be finite and nondecreasing" );
            }
        }
    }

    t_coadd_begin = sys::get_curr_time();

    std::vector<std::pair<int, int>> ranges;
    for( int begin = 0; begin < inputImages; )
    {
        int end = begin + 1;
        while( end < inputImages )
        {
            if( coaddMaxImno > 0 && end - begin >= coaddMaxImno )
            {
                break;
            }
            if( coaddMaxTime > 0 )
            {
                const double elapsedSeconds = ( imageMJD[end] - imageMJD[begin] ) * 86400.0;
                const double timeTolerance = 1e-6 * std::max( 1.0, static_cast<double>( coaddMaxTime ) );
                if( elapsedSeconds - coaddMaxTime > timeTolerance )
                {
                    break;
                }
            }
            ++end;
        }
        ranges.emplace_back( begin, end );
        begin = end;
    }

    std::vector<imageT> coadds;
    std::vector<fitsHeaderT> outputHeads;
    std::vector<double> outputMJD;
    coadds.reserve( ranges.size() );
    outputHeads.reserve( ranges.size() );
    outputMJD.reserve( ranges.size() );

    eigenCube<realT> groupCube;
    imageT coadd;

    for( const auto &[begin, end] : ranges )
    {
        const int groupSize = end - begin;
        groupCube.resize( ims.rows(), ims.cols(), groupSize );
        for( int offset = 0; offset < groupSize; ++offset )
        {
            groupCube.image( offset ) = ims.image( begin + offset );
        }

        if( coaddMethod == HCI::coadd::mean )
        {
            groupCube.mean( coadd );
        }
        else
        {
            groupCube.median( coadd );
        }
        coadds.push_back( coadd );

        fitsHeaderT header = heads[begin];
        if( haveDates )
        {
            double sum = 0;
            for( int index = begin; index < end; ++index )
            {
                sum += imageMJD[index];
            }
            const double average = sum / groupSize;
            outputMJD.push_back( average );

            if( !dateKeyword.empty() )
            {
                header[dateKeyword].value( mx::sys::ISO8601DateTimeStrMJD( average, 1 ) );
                header["START " + dateKeyword].value( mx::sys::ISO8601DateTimeStrMJD( imageMJD[begin], 1 ) );
                header["END " + dateKeyword].value( mx::sys::ISO8601DateTimeStrMJD( imageMJD[end - 1], 1 ) );
                header.append( "DELTA " + dateKeyword,
                               ( imageMJD[end - 1] - imageMJD[begin] ) * 86400,
                               "change in " + dateKeyword + " in seconds." );
            }
        }

        for( const auto &keyword : coaddKeywords )
        {
            const double start = heads[begin][keyword].template value<double>();
            const double finish = heads[end - 1][keyword].template value<double>();
            double sum = 0;
            for( int index = begin; index < end; ++index )
            {
                sum += heads[index][keyword].template value<double>();
            }
            header[keyword].value( sum / groupSize );
            header["START " + keyword].value( start );
            header["END " + keyword].value( finish );
            header["DELTA " + keyword].value( finish - start );
        }

        header.append( "IMAGES COADDED", groupSize, "number of images included in this coadd" );
        for( int index = begin; index < end; ++index )
        {
            header.append( "HISTORY",
                           fits::fitsHistoryType(),
                           "coadded file: " + ioutils::pathFilename( fileList[index] ) );
        }
        outputHeads.push_back( std::move( header ) );
    }

    ims.resize( inputRows, inputCols, coadds.size() );
    for( size_t index = 0; index < coadds.size(); ++index )
    {
        ims.image( index ) = coadds[index];
    }
    heads = std::move( outputHeads );
    if( haveDates )
    {
        imageMJD = std::move( outputMJD );
    }

    t_coadd_end = sys::get_curr_time();

} // void HCIobservation<_realT,verboseT>::coaddImages()

template <typename _realT, class verboseT>
void HCIobservation<_realT, verboseT>::readMaskFile( imageT &mask,
                                                     const std::string &maskFileName,
                                                     const std::string &description )
{
    fitsFileT ff;
    mx::error_t errc = ff.read( mask, maskFileName );
    if( errc != mx::error_t::noerror )
    {
        throw mx::exception<verboseT>( errc, std::format( "Could not read the {} mask.", description ) );
    }

    if( m_imSize <= 0 || mask.rows() < m_imSize || mask.cols() < m_imSize )
    {
        throw mx::exception<verboseT>(
            mx::error_t::sizeerr,
            std::format( "The {} mask is smaller than the input image size.", description ) );
    }

    if( mask.rows() > m_imSize || mask.cols() > m_imSize )
    {
        imageT cropped = mask.block( static_cast<int>( 0.5 * ( mask.rows() - 1 ) - 0.5 * ( m_imSize - 1 ) ),
                                     static_cast<int>( 0.5 * ( mask.cols() - 1 ) - 0.5 * ( m_imSize - 1 ) ),
                                     m_imSize,
                                     m_imSize );
        mask = cropped;
    }
}

template <typename _realT, class verboseT>
void HCIobservation<_realT, verboseT>::readMask()
{
    /*** Load the mask ***/
    if( m_maskFile != "" )
    {
        std::cerr << "creating mask cube\n";
        readMaskFile( m_mask, m_maskFile, "input" );
        try
        {
            makeMaskCube();
        }
        catch( ... )
        {
            std::throw_with_nested( mx::exception<verboseT>( mx::error_t::exception ) );
        }
    }
}

template <typename realT, class verboseT>
void HCIobservation<realT, verboseT>::populateMaskCube(
    eigenCube<realT> &maskCube, const imageT &mask, int rows, int cols, int planes )
{
    if( mask.rows() != rows || mask.cols() != cols )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                       std::format( "Mask is not the same size as images.\n"
                                                    "    Mask:   rows={}\n"
                                                    "            cols={}\n"
                                                    "    Images: rows={}\n"
                                                    "            cols={}\n",
                                                    mask.rows(),
                                                    mask.cols(),
                                                    rows,
                                                    cols ) );
    }

    maskCube.resize( rows, cols, planes );

    for( int i = 0; i < planes; ++i )
    {
        maskCube.image( i ) = mask;
    }
}

template <typename realT, class verboseT>
void HCIobservation<realT, verboseT>::makeMaskCube()
{
    populateMaskCube( m_maskCube, m_mask, m_Nrows, m_Ncols, m_Nims );
}

template <typename realT, class verboseT>
const typename HCIobservation<realT, verboseT>::imageT *
HCIobservation<realT, verboseT>::preProcessMask( const eigenCube<realT> &ims ) const
{
    if( &ims == &m_refIms )
    {
        if( ( m_RDImaskUseInput || !m_RDImaskFile.empty() ) && m_RDImask.size() != 0 )
        {
            return &m_RDImask;
        }
        return nullptr;
    }

    if( !m_maskFile.empty() && m_mask.size() != 0 )
    {
        return &m_mask;
    }
    return nullptr;
}

template <typename _realT, class verboseT>
void HCIobservation<_realT, verboseT>::preProcess( eigenCube<realT> &ims )
{
    t_preproc_begin = sys::get_curr_time();
    const imageT *mask = preProcessMask( ims );

    // The mask is applied first, and then after each subsequent P.P. step.
    if( mask != nullptr && m_preProcess_mask )
    {
        std::cerr << "Masking . . .\n";
#pragma omp parallel for
        for( int i = 0; i < ims.planes(); ++i )
        {
            ims.image( i ) *= *mask;
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
                if( mask != nullptr && m_preProcess_mask )
                {
                    imageT radius;
                    radius.resize( imRef.rows(), imRef.cols() );
                    radiusImage( radius );
                    radprofim( rp, imRef, radius, mask, true );
                }
                else
                {
                    radprofim( rp, imRef, true );
                }
                zeroNaNs( imRef );
            }
        }

        std::cerr << "done\n";

        if( mask != nullptr && m_preProcess_mask )
        {
            std::cerr << "Masking . . .\n";
#pragma omp parallel for
            for( int i = 0; i < ims.planes(); ++i )
            {
                ims.image( i ) *= *mask;
            }
            std::cerr << "done\n";
        }
    }

    if( m_preProcess_medianUSM_fwhm > 0 )
    {
        if( m_preProcess_medianUSM_fwhm > ims.rows() || m_preProcess_medianUSM_fwhm > ims.cols() )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "The median USM width cannot exceed the image dimensions." );
        }

        std::cerr << "applying median USM . . .\n";

        imageT medmask;
        medmask.resize( ims.rows(), ims.cols() );
        medmask.setConstant( 1 );

        const int before = m_preProcess_medianUSM_fwhm / 2;
        const int after = m_preProcess_medianUSM_fwhm - before - 1;
        for( int z = 0; z < before; ++z )
        {
            medmask.row( z ) = 0;
            medmask.col( z ) = 0;
        }
        for( int z = 0; z < after; ++z )
        {
            medmask.row( medmask.rows() - 1 - z ) = 0;
            medmask.col( medmask.cols() - 1 - z ) = 0;
        }

#pragma omp parallel
        {
            imageT fim, im;
            fim.resize( ims.rows(), ims.cols() );
            fim.setZero();
            im.resize( ims.rows(), ims.cols() );

#pragma omp for
            for( int i = 0; i < ims.planes(); ++i )
            {
                im = ims.image( i );
                // Width and output dimensions were validated above, which is medianSmooth's complete error contract.
                static_cast<void>( medianSmooth( fim, im, m_preProcess_medianUSM_fwhm ) );
                ims.image( i ) = ( im - fim ) * medmask;
            }
        }

        if( mask != nullptr && m_preProcess_mask )
        {
            std::cerr << "Masking . . .\n";
#pragma omp parallel for
            for( int i = 0; i < ims.planes(); ++i )
            {
                ims.image( i ) *= *mask;
            }
        }
        std::cerr << "done\n";
    }

    if( !std::isfinite( m_preProcess_gaussUSM_fwhm ) || m_preProcess_gaussUSM_fwhm < 0 )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                       "The Gaussian USM FWHM must be finite and nonnegative." );
    }

    if( m_preProcess_gaussUSM_fwhm > 0 )
    {
        const realT scaledWidth = 2 * m_preProcess_gaussUSM_fwhm;
        if( scaledWidth > std::numeric_limits<int>::max() )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "The Gaussian USM kernel exceeds the supported dimensions." );
        }

        int kernelWidth = static_cast<int>( scaledWidth );
        if( kernelWidth % 2 == 0 )
        {
            ++kernelWidth;
        }
        if( kernelWidth > ims.rows() || kernelWidth > ims.cols() )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "The Gaussian USM kernel cannot exceed the image dimensions." );
        }

        std::cerr << "applying Gauss USM . . .\n";
        t_gaussusm_begin = sys::get_curr_time();
        const gaussKernel<eigenImage<_realT>, 2> kernel( m_preProcess_gaussUSM_fwhm );
#pragma omp parallel for
        for( int i = 0; i < ims.planes(); ++i )
        {
            imageT fim, im;
            im = ims.image( i );
            // gaussKernel::setKernel always succeeds, so this filter specialization cannot return an error.
            static_cast<void>(
                filterImage( fim,
                             im,
                             kernel,
                             0.5 * ( std::min( ims.rows(), ims.cols() ) - 1 ) - m_preProcess_gaussUSM_fwhm * 4 ) );
            im = ( im - fim );
            ims.image( i ) = im;
        }

        if( mask != nullptr && m_preProcess_mask )
        {
            std::cerr << "Masking . . .\n";
            // clang-format off
            #pragma omp parallel for // clang-format on
            for( int i = 0; i < ims.planes(); ++i )
            {
                ims.image( i ) *= *mask;
            }
        }
        t_gaussusm_end = sys::get_curr_time();
        std::cerr << "done\n";
    }

    if( m_preProcess_azUSM_azHalfWidth && m_preProcess_azUSM_radHalfWidth )
    {
        ipc::ompLoopWatcher<> status( ims.planes(), std::cerr );

        std::cerr << "applying azimuthal USM . . .\n";
        t_azusm_begin = sys::get_curr_time();

        azBoxKernel<eigenImage<realT>> azbK( m_preProcess_azUSM_radHalfWidth,
                                             m_preProcess_azUSM_azHalfWidth,
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

        if( mask != nullptr && m_preProcess_mask )
        {
            std::cerr << "Masking . . .\n";
#pragma omp parallel for
            for( int i = 0; i < ims.planes(); ++i )
            {
                ims.image( i ) *= *mask;
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
    if( m_preProcess_meanSubMethod == HCI::meanSub::none )
    {
        return;
    }
    else if( m_preProcess_meanSubMethod != HCI::meanSub::meanImage &&
             m_preProcess_meanSubMethod != HCI::meanSub::medianImage )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                       std::format( "Mean subtraction by {} "
                                                    "can't be done in pre-processing. "
                                                    "Only meanImage or medianImage can be used in pre.",
                                                    HCI::meanSubToStr<verboseT>( m_preProcess_meanSubMethod ) ) );
    }

    imageT mean;
    const imageT *mask = preProcessMask( ims );

    if( m_preProcess_meanSubMethod == HCI::meanSub::meanImage )
    {
        ims.mean( mean );
    }
    else if( m_preProcess_meanSubMethod == HCI::meanSub::medianImage )
    {
        ims.median( mean );
    }

#pragma omp parallel for
    for( int n = 0; n < ims.planes(); ++n )
    {
        ims.image( n ) -= mean;

        if( mask != nullptr && m_preProcess_mask )
        {
            ims.image( n ) *= *mask;
        }
    }
}

template <typename _realT, class verboseT>
void HCIobservation<_realT, verboseT>::preProcess_pixelTSNorm( eigenCube<realT> &ims )
{
    if( m_preProcess_pixelTSNormMethod == HCI::pixelTSNorm::none )
    {
        return;
    }

    if( m_preProcess_pixelTSNormMethod == HCI::pixelTSNorm::rmsSigmaClipped )
    {
        throw mx::exception<verboseT>( mx::error_t::notimpl,
                                       "pixelTSNormMethod is rmsSigmaClipped, which is not implemented" );
    }

    if( m_preProcess_pixelTSNormMethod != HCI::pixelTSNorm::rms )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                       "pixelTSNormMethod is not a supported normalization method" );
    }

    if( ims.planes() == 0 )
    {
        return;
    }

    std::cerr << "normalizing pixels\n";
    const imageT *mask = preProcessMask( ims );

    // clang-format off
    #pragma omp parallel // clang-format on
    {
        std::vector<realT> pixs( ims.planes() );

        // clang-format off
        #pragma omp for // clang-format on
        for( int cc = 0; cc < ims.cols(); ++cc )
        {
            for( int rr = 0; rr < ims.rows(); ++rr )
            {
                if( mask != nullptr && m_preProcess_mask )
                {
                    if( ( *mask )( rr, cc ) == 0 )
                    {
                        continue;
                    }
                }

                // We bother to load a vector in prep to add sigma clipping later.
                for( int pp = 0; pp < ims.planes(); ++pp )
                {
                    pixs[pp] = ims.image( pp )( rr, cc );
                }

                realT sumSquares = 0;
                for( int pp = 0; pp < ims.planes(); ++pp )
                {
                    sumSquares += pixs[pp] * pixs[pp];
                }

                realT sd = std::sqrt( sumSquares / static_cast<realT>( ims.planes() ) );

                if( sd == 0 || !std::isfinite( sd ) )
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
void HCIobservation<_realT, verboseT>::readWeights()
{
    if( m_weightFile == "" )
    {
        throw mx::exception<verboseT>( mx::error_t::paramnotset, "m_weightFile not set" );
    }

    // Read the weight file and load it into a map
    std::vector<std::string> wfileNames;
    std::vector<realT> imW;

    mx::error_t errc = ioutils::readColumns<mx::ioutils::readColSpaceDelim, verboseT>( m_weightFile, wfileNames, imW );
    if( errc != mx::error_t::noerror )
    {
        throw mx::exception<verboseT>( errc );
    }

    if( wfileNames.size() != imW.size() )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "weight file columns have different lengths" );
    }

    if( imW.size() < m_fileList.size() )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "not enough weights specified" );
    }

    std::map<std::string, realT> weights;

    for( size_t i = 0; i < wfileNames.size(); ++i )
    {
        std::string fname = ioutils::pathFilename( wfileNames[i].c_str() );
        if( !std::isfinite( imW[i] ) )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidarg, "non-finite weight for filename " + fname );
        }
        if( !weights.emplace( fname, imW[i] ).second )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidarg, "duplicate weight for filename " + fname );
        }
    }

    m_comboWeights.resize( m_fileList.size() );

    realT wi;
    realT weightSum = 0;
    for( size_t i = 0; i < m_fileList.size(); ++i )
    {
        try
        {
            wi = weights.at( ioutils::pathFilename( m_fileList[i].c_str() ) );
        }
        catch( const std::exception &e )
        {
            std::throw_with_nested( mx::exception<verboseT>( mx::error_t::std_exception, "finding weights for file" ) );
        }

        m_comboWeights[i] = wi;
        weightSum += wi;
    }

    if( !std::isfinite( weightSum ) || weightSum == 0 )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidarg, "weight sum must be finite and non-zero" );
    }

    // Finally normalize the weights
    for( size_t i = 0; i < m_comboWeights.size(); ++i )
    {
        m_comboWeights[i] /= weightSum;
    }

    // return 0;
} // void HCIobservation<_realT,verboseT>::readWeights()

template <typename _realT, class verboseT>
void HCIobservation<_realT, verboseT>::combineFinim()
{
    if( m_combineMethod == HCI::combine::none )
    {
        return;
    }

    if( m_psfsub.empty() )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "no PSF-subtracted image cubes to combine" );
    }

    const int rows = m_psfsub.front().rows();
    const int cols = m_psfsub.front().cols();
    const int planes = m_psfsub.front().planes();
    if( rows <= 0 || cols <= 0 || planes <= 0 )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "PSF-subtracted image cubes must not be empty" );
    }

    for( const auto &cube : m_psfsub )
    {
        if( cube.rows() != rows || cube.cols() != cols || cube.planes() != planes )
        {
            throw mx::exception<verboseT>( mx::error_t::sizeerr,
                                           "PSF-subtracted image cubes must have matching dimensions" );
        }
    }

    if( !m_comboWeights.empty() && m_comboWeights.size() != static_cast<size_t>( planes ) )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "combination weight count must match the image count" );
    }

    if( m_combineMethod != HCI::combine::mean && m_combineMethod != HCI::combine::median &&
        m_combineMethod != HCI::combine::sigmaMean )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig, "invalid final image combination method" );
    }

    t_combo_begin = sys::get_curr_time();

    // Create and size temporary image for averaging
    imageT tfinim;

    m_finim.resize( rows, cols, m_psfsub.size() );

    // Now cycle through each set of psf subtractions
    for( size_t n = 0; n < m_psfsub.size(); ++n )
    {
        if( m_combineMethod == HCI::combine::median )
        {
            if( m_maskFile != "" )
            {
                m_psfsub[n].median( tfinim, m_maskCube, m_minGoodFract );
            }
            else
            {
                m_psfsub[n].median( tfinim );
            }
            m_finim.image( n ) = tfinim;
        }
        else if( m_combineMethod == HCI::combine::mean )
        {
            if( !m_comboWeights.empty() )
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
        else if( m_combineMethod == HCI::combine::sigmaMean )
        {
            if( m_sigmaThreshold <= 0 )
            {
                if( !m_comboWeights.empty() )
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
                else if( m_maskFile != "" )
                {
                    m_psfsub[n].mean( tfinim, m_maskCube, m_minGoodFract );
                }
                else
                {
                    m_psfsub[n].mean( tfinim );
                }
            }
            else if( !m_comboWeights.empty() )
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
    if( !dir.empty() )
    {
        const mx::error_t result = ioutils::createDirectories( dir );
        if( result != mx::error_t::noerror )
        {
            throw mx::exception<verboseT>( result, "creating the preprocessed-image output directory" );
        }
    }

    if( m_tgtIms.planes() <= 0 || m_tgtIms.rows() <= 0 || m_tgtIms.cols() <= 0 )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "the target image cube is empty" );
    }
    if( m_heads.size() != static_cast<size_t>( m_tgtIms.planes() ) )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "the target header count does not match the image count" );
    }

    std::string fname;

    fitsFileT ff;

    /** \todo Should add a HISTORY card here */
    for( int i = 0; i < m_tgtIms.planes(); ++i )
    {
        fname = std::format( "{}{:06}.fits", m_preProcess_outputPrefix, i );

        fitsHeaderT fh = m_heads[i];
        stdFitsHeader( fh );
        const mx::error_t result =
            ff.write( fname, m_tgtIms.image( i ).data(), m_tgtIms.rows(), m_tgtIms.cols(), 1, fh );
        if( result != mx::error_t::noerror )
        {
            throw mx::exception<verboseT>( result, "writing preprocessed target image " + fname );
        }
    }
} // void HCIobservation<_realT,verboseT>::outputPreProcessed()

template <typename _realT, class verboseT>
void HCIobservation<_realT, verboseT>::outputRDIPreProcessed()
{
    if( m_preProcess_outputPrefix.empty() )
    {
        return;
    }

    const std::string dir = ioutils::parentPath( m_preProcess_outputPrefix );
    if( !dir.empty() )
    {
        const mx::error_t result = ioutils::createDirectories( dir );
        if( result != mx::error_t::noerror )
        {
            throw mx::exception<verboseT>( result, "creating the preprocessed-RDI output directory" );
        }
    }

    if( m_refIms.planes() <= 0 || m_refIms.rows() <= 0 || m_refIms.cols() <= 0 )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "the RDI image cube is empty" );
    }
    if( m_RDIheads.size() != static_cast<size_t>( m_refIms.planes() ) )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "the RDI header count does not match the image count" );
    }

    fitsFileT ff;
    for( int i = 0; i < m_refIms.planes(); ++i )
    {
        const std::string fname = std::format( "{}RDI_{:06}.fits", m_preProcess_outputPrefix, i );
        fitsHeaderT head = m_RDIheads[i];
        stdFitsHeader( head );
        const mx::error_t result =
            ff.write( fname, m_refIms.image( i ).data(), m_refIms.rows(), m_refIms.cols(), 1, head );
        if( result != mx::error_t::noerror )
        {
            throw mx::exception<verboseT>( result, "writing preprocessed RDI image " + fname );
        }
    }
} // void HCIobservation<_realT,verboseT>::outputRDIPreProcessed()

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

    head.template append<std::string>( "COADMTHD",
                                       HCI::coaddToStr<verboseT>( m_coaddMethod ),
                                       "coadd combination method" );
    if( m_coaddMethod != HCI::coadd::none )
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
                                 m_preProcess_azUSM_azHalfWidth,
                                 "pre-process azimuthal USM azimuthal half-width [pixels]" );
    head.template append<realT>( "PREPROC AZUSM MAXAZ",
                                 m_preProcess_azUSM_maxAz,
                                 "pre-process azimuthal USM maximum azimuthal width [degrees]" );
    head.template append<realT>( "PREPROC AZUSM RADWIDTH",
                                 m_preProcess_azUSM_radHalfWidth,
                                 "pre-process azimuthal USM radial half-width [pixels]" );
    head.template append<realT>( "PREPROC MEDIANUSM FWHM",
                                 m_preProcess_medianUSM_fwhm,
                                 "pre-process median USM fwhm [pixels]" );
    head.template append<realT>( "PREPROC GAUSSUSM FWHM",
                                 m_preProcess_gaussUSM_fwhm,
                                 "pre-process Gaussian USM fwhm [pixels]" );
    head.template append<std::string>( "PREPROC MEANSUB METHOD",
                                       HCI::meanSubToStr<verboseT>( m_preProcess_meanSubMethod ),
                                       "pre-process mean subtraction method" );
    head.template append<std::string>( "PREPROC PIXELTSNORM METHOD",
                                       HCI::pixelTSNormToStr<verboseT>( m_preProcess_pixelTSNormMethod ),
                                       "pre-process pixel time-series norm method" );
}

template <typename _realT, class verboseT>
void HCIobservation<_realT, verboseT>::writeFinim( fitsHeaderT *addHead )
{
    if( m_finim.rows() <= 0 || m_finim.cols() <= 0 || m_finim.planes() <= 0 )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "the final image cube is empty" );
    }

    std::string fname = m_finimName;

    if( !m_outputDir.empty() )
    {
        fname = m_outputDir + "/" + fname;
    }

    const std::string outputParent = ioutils::parentPath( fname );
    if( !outputParent.empty() )
    {
        const mx::error_t result = ioutils::createDirectories( outputParent );
        if( result != mx::error_t::noerror )
        {
            throw mx::exception<verboseT>( result, "creating the final-image output directory" );
        }
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
                                       HCI::combineToStr<verboseT>( m_combineMethod ),
                                       "combination method" );

    if( m_weightFile != "" )
        head.append( "WEIGHT FILE", m_weightFile, "file containing weights for combination" );

    if( m_combineMethod == HCI::combine::sigmaMean )
        head.template append<realT>( "SIGMA THRESHOLD", m_sigmaThreshold, "threshold for sigma clipping" );

    head.template append<realT>( "MIN GOOD FRACTION", m_minGoodFract, "minimum good fraction for inclusion" );
    if( addHead )
    {
        head.append( *addHead );
    }

    fits::fitsHeaderGitStatus( head, "mxlib_uncomp", MXLIB_UNCOMP_CURRENT_SHA1, MXLIB_UNCOMP_REPO_MODIFIED );

    fitsFileT f;

    const mx::error_t result = f.write( fname, m_finim, head );
    if( result != mx::error_t::noerror )
    {
        throw mx::exception<verboseT>( result, "writing final image " + fname );
    }

    std::cerr << "Final image written to: " << fname << "\n";
} // void HCIobservation<_realT,verboseT>::writeFinim(fitsHeaderT * addHead)

template <typename _realT, class verboseT>
void HCIobservation<_realT, verboseT>::outputPSFSub( fitsHeaderT *addHead )
{
    if( m_psfsub.empty() )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "the PSF-subtracted image collection is empty" );
    }

    const int rows = m_psfsub.front().rows();
    const int cols = m_psfsub.front().cols();
    const int planes = m_psfsub.front().planes();
    if( rows <= 0 || cols <= 0 || planes <= 0 )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "the PSF-subtracted image cubes are empty" );
    }
    for( const auto &cube : m_psfsub )
    {
        if( cube.rows() != rows || cube.cols() != cols || cube.planes() != planes )
        {
            throw mx::exception<verboseT>( mx::error_t::sizeerr,
                                           "the PSF-subtracted image cubes have inconsistent dimensions" );
        }
    }
    if( m_heads.size() != static_cast<size_t>( planes ) )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr,
                                       "the target header count does not match the PSF-subtracted image count" );
    }
    if( !m_comboWeights.empty() && m_comboWeights.size() != static_cast<size_t>( planes ) )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr,
                                       "the combination weight count does not match the PSF-subtracted image count" );
    }

    const auto outputPath = [this]( const std::string &name )
    {
        if( m_outputDir.empty() )
        {
            return name;
        }
        return m_outputDir + "/" + name;
    };

    const auto ensureParent = []( const std::string &path )
    {
        const std::string parent = ioutils::parentPath( path );
        if( parent.empty() )
        {
            return;
        }
        const mx::error_t result = ioutils::createDirectories( parent );
        if( result != mx::error_t::noerror )
        {
            throw mx::exception<verboseT>( result, "creating a PSF-subtracted output directory" );
        }
    };

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

    if( !m_comboWeights.empty() )
    {
        const std::string weightName = outputPath( m_PSFSubPrefix + "weights.dat" );
        ensureParent( weightName );
        wout.open( weightName );
        if( !wout )
        {
            throw mx::exception<verboseT>( mx::error_t::fileoerr, "opening combination-weight output " + weightName );
        }
        std::cerr << "writing comboWeights: " << weightName << "\n";
    }

    for( size_t n = 0; n < m_psfsub.size(); ++n )
    {
        for( int p = 0; p < m_psfsub[n].planes(); ++p )
        {
            const std::string fname = outputPath( std::format( "{}_{:03}_{:05}.fits", m_PSFSubPrefix, n, p ) );
            ensureParent( fname );

            fitsHeaderT h = head;

            h.append( m_heads[p] );
            h.template append<int>( "REDUCTION", static_cast<int>( n ), "reduction index" );
            h.template append<int>( "IMAGE", p, "input image index" );

            const mx::error_t result =
                f.write( fname, m_psfsub[n].image( p ).data(), m_psfsub[n].rows(), m_psfsub[n].cols(), 1, h );
            if( result != mx::error_t::noerror )
            {
                throw mx::exception<verboseT>( result, "writing PSF-subtracted image " + fname );
            }

            if( !m_comboWeights.empty() && n == 0 )
            {
                wout << fname << " " << m_comboWeights[p] << "\n";
                wout.flush();
                if( !wout )
                {
                    throw mx::exception<verboseT>( mx::error_t::fileoerr, "writing combination weights" );
                }
            }
        }
    }

    if( !m_comboWeights.empty() )
    {
        wout.close();
    }
} // void HCIobservation<_realT,verboseT>::outputPSFSub(fitsHeaderT * addHead)

template <typename _realT, class verboseT>
void HCIobservation<_realT, verboseT>::readPSFSub( const std::string &directory,
                                                   const std::string &prefix,
                                                   const std::string &extension )
{
    namespace fs = std::filesystem;

    if( prefix.empty() || extension.empty() )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidarg,
                                       "saved PSF-subtracted prefix and extension must not be empty" );
    }

    std::error_code filesystemError;
    if( !fs::is_directory( directory, filesystemError ) || filesystemError )
    {
        throw mx::exception<verboseT>( mx::error_t::dirnotfound, "saved PSF-subtracted directory is not readable" );
    }

    using indexT = std::pair<size_t, size_t>;
    std::map<indexT, std::string> indexedFiles;
    const std::string filenamePrefix = prefix + "_";

    fs::directory_iterator entry( directory, filesystemError );
    const fs::directory_iterator end;
    if( filesystemError )
    {
        throw mx::exception<verboseT>( mx::error_t::dirnotfound,
                                       "could not inspect the saved PSF-subtracted directory" );
    }

    while( entry != end )
    {
        const bool isRegularFile = entry->is_regular_file( filesystemError );
        if( filesystemError )
        {
            throw mx::exception<verboseT>( mx::error_t::fileoerr,
                                           "could not inspect an entry in the saved PSF-subtracted directory" );
        }

        if( isRegularFile )
        {
            const std::string filename = entry->path().filename().string();
            if( filename.starts_with( filenamePrefix ) && filename.ends_with( extension ) )
            {
                const size_t bodyLength = filename.size() - filenamePrefix.size() - extension.size();
                const std::string_view body( filename.data() + filenamePrefix.size(), bodyLength );
                const size_t separator = body.find( '_' );
                if( separator == std::string_view::npos || separator == 0 || separator + 1 >= body.size() ||
                    body.find( '_', separator + 1 ) != std::string_view::npos )
                {
                    throw mx::exception<verboseT>( mx::error_t::parseerr,
                                                   "malformed saved PSF-subtracted filename " + filename );
                }

                size_t reduction = 0;
                size_t image = 0;
                const auto reductionResult = std::from_chars( body.data(), body.data() + separator, reduction );
                const auto imageResult =
                    std::from_chars( body.data() + separator + 1, body.data() + body.size(), image );
                if( reductionResult.ec != std::errc{} || reductionResult.ptr != body.data() + separator ||
                    imageResult.ec != std::errc{} || imageResult.ptr != body.data() + body.size() )
                {
                    throw mx::exception<verboseT>( mx::error_t::parseerr,
                                                   "malformed saved PSF-subtracted filename " + filename );
                }

                const auto [inserted, unique] =
                    indexedFiles.emplace( indexT{ reduction, image }, entry->path().string() );
                static_cast<void>( inserted );
                if( !unique )
                {
                    throw mx::exception<verboseT>( mx::error_t::invalidarg,
                                                   "multiple saved PSF-subtracted files have the same indices" );
                }
            }
        }

        entry.increment( filesystemError );
        if( filesystemError )
        {
            throw mx::exception<verboseT>( mx::error_t::fileoerr,
                                           "could not traverse the saved PSF-subtracted directory" );
        }
    }

    if( indexedFiles.empty() )
    {
        throw mx::exception<verboseT>( mx::error_t::filenotfound, "no saved PSF-subtracted images were found" );
    }

    const size_t reductionCount = indexedFiles.rbegin()->first.first + 1;
    size_t imageCount = 0;
    for( const auto &[index, path] : indexedFiles )
    {
        static_cast<void>( path );
        imageCount = std::max( imageCount, index.second + 1 );
    }
    if( indexedFiles.size() != reductionCount * imageCount )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr,
                                       "saved PSF-subtracted images do not form a complete index grid" );
    }

    std::vector<eigenCube<realT>> psfsub( reductionCount );
    std::vector<fitsHeaderT> heads( imageCount );
    std::vector<std::string> fileList( imageCount );
    fitsFileT fitsFile;
    int rows = 0;
    int cols = 0;

    for( size_t reduction = 0; reduction < reductionCount; ++reduction )
    {
        for( size_t imageIndex = 0; imageIndex < imageCount; ++imageIndex )
        {
            const std::string &filename = indexedFiles.at( indexT{ reduction, imageIndex } );

            imageT image;
            fitsHeaderT header;
            const mx::error_t result = fitsFile.read( image, header, filename );
            if( result != mx::error_t::noerror )
            {
                throw mx::exception<verboseT>( result, "reading saved PSF-subtracted image " + filename );
            }
            if( header.count( "REDUCTION" ) == 0 || header.count( "IMAGE" ) == 0 ||
                header["REDUCTION"].template value<int>() != static_cast<int>( reduction ) ||
                header["IMAGE"].template value<int>() != static_cast<int>( imageIndex ) )
            {
                throw mx::exception<verboseT>( mx::error_t::invalidarg,
                                               "saved PSF-subtracted header indices do not match the filename" );
            }

            if( reduction == 0 && imageIndex == 0 )
            {
                rows = image.rows();
                cols = image.cols();
                for( auto &cube : psfsub )
                {
                    cube.resize( rows, cols, imageCount );
                }
            }
            else if( image.rows() != rows || image.cols() != cols )
            {
                throw mx::exception<verboseT>( mx::error_t::sizeerr,
                                               "saved PSF-subtracted images have inconsistent dimensions" );
            }

            psfsub[reduction].image( imageIndex ) = image;
            if( reduction == 0 )
            {
                heads[imageIndex] = header;
                fileList[imageIndex] = filename;
            }
        }
    }

    m_psfsub = std::move( psfsub );
    m_heads = std::move( heads );
    m_fileList = std::move( fileList );
    m_Nims = static_cast<int>( imageCount );
    m_Nrows = rows;
    m_Ncols = cols;
    m_Npix = rows * cols;
    m_imSize = std::min( rows, cols );
    m_filesRead = true;
    postReadFiles();
}

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
    m_preProcess_azUSM_azHalfWidth = fh["PPAUSMAW"].template value<realT>();
    std::cerr << "preProcess_azUSM_azHalfWidth: " << m_preProcess_azUSM_azHalfWidth << "\n";

    if( fh.count( "PPAUSMRW" ) == 0 )
    {
        mxError( "KLIPReduction", MXE_PARAMNOTSET, "PPAUSMRW not found in FITS header." );
        return -1;
    }
    m_preProcess_azUSM_radHalfWidth = fh["PPAUSMRW"].template value<realT>();
    std::cerr << "preProcess_azUSM_radHalfWidth: " << m_preProcess_azUSM_radHalfWidth << "\n";

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
        ioutils::readColumns<mx::ioutils::readColSpaceDelim, verboseT>( m_weightFile, fn, m_comboWeights );

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

extern template class HCIobservation<float, mx::verbose::vv>;

} // namespace improc
} // namespace mx

#endif //__HCIobservation_hpp__
