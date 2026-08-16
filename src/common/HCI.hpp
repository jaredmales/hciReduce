/** \file HCI.hpp
 * \author Jared R. Males
 * \brief Defines the configuration types for high contrast imaging.
 *
 */

#ifndef __HCI_hpp__
#define __HCI_hpp__

#include <string>
#include <mx/mxlib.hpp>

namespace mx
{
namespace improc
{

/// Namespace for high contrast imaging configuration enums.
/** \ingroup programming_library
 */
namespace HCI
{

/// Possible coadding methods
/** \ingroup programming_library
 */
enum class coadd
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
template <class verboseT>
std::string coaddToStr( coadd method /**< [in] one of the \ref coadd method enum members */ )
{
    if( method == coadd::none )
    {
        return "none";
    }
    else if( method == coadd::median )
    {
        return "median";
    }
    else if( method == coadd::mean )
    {
        return "mean";
    }
    else
    {
        throw mx::exception<verboseT>( mx::error_t::invalidarg, "got an invalid coadd method (bug)" );
    }
}

/// Get the coaddMethod from the corresponding string name
/**
 * \returns the coaddMethod enum member corresponding to the string name.
 */
template <class verboseT>
coadd coaddFmStr( const std::string &method /**< [in] the string name of the \ref coadd method */ )
{
    if( method == "none" )
    {
        return coadd::none;
    }
    else if( method == "median" )
    {
        return coadd::median;
    }
    else if( method == "mean" )
    {
        return coadd::mean;
    }
    else
    {
        throw mx::exception<verboseT>( mx::error_t::invalidarg, method + " is not a valid coadd method" );
    }
}

/// Mean subtraction methods
/** These control how the data in each search region is centered to meet the PCA
 * requirement. \ingroup programming_library
 */
enum class meanSub
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

template <class verboseT>
std::string meanSubToStr( const meanSub &method /**< [in] The \ref meanSub method to convert to its string name*/ )
{
    if( method == meanSub::none )
    {
        return "none";
    }
    else if( method == meanSub::meanImage )
    {
        return "meanImage";
    }
    else if( method == meanSub::medianImage )
    {
        return "medianImage";
    }
    else if( method == meanSub::imageMean )
    {
        return "imageMean";
    }
    else if( method == meanSub::imageMedian )
    {
        return "imageMedian";
    }
    else if( method == meanSub::imageMode )
    {
        return "imageMode";
    }
    else
    {
        throw mx::exception( mx::error_t::invalidarg, "got an invalid mean sub method (bug)" );
    }
}

template <class verboseT>
meanSub meanSubFmStr( const std::string &method /**< [in] The string name of a \ref meanSub method to convert */ )
{
    if( method == "none" )
    {
        return meanSub::none;
    }
    else if( method == "meanImage" )
    {
        return meanSub::meanImage;
    }
    else if( method == "medianImage" )
    {
        return meanSub::medianImage;
    }
    else if( method == "imageMean" )
    {
        return meanSub::imageMean;
    }
    else if( method == "imageMedian" )
    {
        return meanSub::imageMedian;
    }
    else if( method == "imageMode" )
    {
        return meanSub::imageMode;
    }
    else
    {
        throw mx::exception( mx::error_t::invalidarg, method + " is not a valid mean sub method" );
    }
}

enum class pixelTSNorm
{
    none,           ///< no pixel time series norm
    rms,            ///< the rms of the pixel time series
    rmsSigmaClipped ///< the sigma clipped rms of the pixel time series
};

template <class verboseT>
std::string
pixelTSNormToStr( const pixelTSNorm &method /**< [in] The \ref pixelTSNorm method to convert to its string name*/ )
{
    if( method == pixelTSNorm::none )
    {
        return "none";
    }
    else if( method == pixelTSNorm::rms )
    {
        return "rms";
    }
    else if( method == pixelTSNorm::rmsSigmaClipped )
    {
        return "rmsSigmaClipped";
    }
    else
    {
        throw mx::exception<verboseT>( mx::error_t::invalidarg, "got an invalid pixelTSNorm method (bug)" );
    }
}

template <class verboseT>
pixelTSNorm pixelTSNormFmStr( const std::string &method /**< [in] The string name of a \ref
                                                                  pixelTSNorm method to convert */
)
{
    if( method == "none" )
    {
        return pixelTSNorm::none;
    }
    else if( method == "rms" )
    {
        return pixelTSNorm::rms;
    }
    else if( method == "rmsSigmaClipped" )
    {
        return pixelTSNorm::rmsSigmaClipped;
    }
    else
    {
        throw mx::exception<verboseT>( mx::error_t::invalidarg, method + " is not a valid pixelTSNorm method" );
    }
}

/// Possible combination methods
/** \ingroup programming_library
 */
enum class combine
{
    none,     ///< Do not combine the images.
    median,   ///< Combine with the median.
    mean,     ///< Combine with the mean.
    sigmaMean ///< Combine with the sigma clipped mean.
};

/// Get the string name of the \ref combine method
/**
 * \returns a string with the name of the \ref combine method
 */
template <class verboseT>
std::string combineToStr( combine method /**< [in] one of the \ref combine enum members */ )
{
    if( method == combine::none )
    {
        return "none";
    }
    else if( method == combine::median )
    {
        return "median";
    }
    else if( method == combine::mean )
    {
        return "mean";
    }
    else if( method == combine::sigmaMean )
    {
        return "sigmaMean";
    }
    else
    {
        throw mx::exception<verboseT>( mx::error_t::invalidarg, "got an invalid combine method (bug)" );
    }
}

/// Get the \ref combine method from the corresponding string name
/**
 * \returns the \ref combine methods enum member corresponding to the string name.
 */
template <class verboseT>
combine combineFmStr( const std::string &method /**< [in] the string name of the \ref combine method */ )
{
    if( method == "none" )
    {
        return combine::none;
    }
    else if( method == "median" )
    {
        return combine::median;
    }
    else if( method == "mean" )
    {
        return combine::mean;
    }
    else if( method == "sigmaMean" )
    {
        return combine::sigmaMean;
    }
    else
    {
        throw mx::exception<verboseT>( mx::error_t::invalidarg, method + " is not a valid combine method" );
    }
}

/// Fake injection PSF file specification methods
/** \ingroup programming_library
 */
enum class fake
{
    single, ///< A single PSF is used
    list    ///< A list of PSF files, one per input image, is used.
};

/// Get the string name of a \ref fake injection method
/**
 * \returns the string name corresponding to the \ref fake injection method
 */
template <class verboseT>
std::string fakeToStr( const fake &method /**< [in] the \ref fake injection method */
)
{
    if( method == fake::single )
    {
        return "single";
    }
    else if( method == fake::list )
    {
        return "list";
    }
    else
    {
        throw mx::exception<verboseT>( mx::error_t::invalidarg, "got an invalid fake method (bug)" );
    }
}

/// Get the \ref fake injection method from its string name
/**
 * \returns the corresponding member of the \ref fake methods enum
 */
template <class verboseT>
fake fakeFmStr( const std::string &method /**< [in] the \ref fake injection method name*/ )
{
    if( method == "single" )
    {
        return fake::single;
    }
    else if( method == "list" )
    {
        return fake::list;
    }
    else
    {
        throw mx::exception<verboseT>( mx::error_t::invalidarg, method + " is not a valid fake method" );
    }
}

/// Image exclusion methods
/** \ingroup programming_library
 */
enum class exclude
{
    none,  ///< Exclude no images
    pixel, ///< Exclude by pixels of rotation
    angle, ///< Exclude by angle of rotation
    imno   ///< Exclude by number of images
};

/// Get the string name of an \ref exclude exclusion method
/**
 * \returns the string name corresponding to the \ref exclude
 *
 * \throws mx::exception with mx::error_t::invalidarg if the method is not valid (this is a bug)
 */
template <class verboseT>
std::string excludeToStr( const exclude &method /**< [in] the \ref exclude method to convert to string */ )
{
    if( method == exclude::none )
    {
        return "none";
    }
    else if( method == exclude::pixel )
    {
        return "pixel";
    }
    else if( method == exclude::angle )
    {
        return "angle";
    }
    else if( method == exclude::imno )
    {
        return "imno";
    }
    else
    {
        throw mx::exception<verboseT>( mx::error_t::invalidarg, "got an invalid exclude method (bug)" );
    }
}

/// Get the \ref exclude exclusion method from its string name
/**
 * \returns the \ref exclude corresponding to the name
 *
 * \throws mx::exception with mx::error_t::invalidarg if string name is not valid
 */
template <class verboseT>
exclude excludeFmStr( const std::string &method /**< [in] the string name of the \ref exclude method */ )
{
    if( method == "none" )
    {
        return exclude::none;
    }
    else if( method == "pixel" )
    {
        return exclude::pixel;
    }
    else if( method == "angle" )
    {
        return exclude::angle;
    }
    else if( method == "imno" )
    {
        return exclude::imno;
    }
    else
    {
        throw mx::exception<verboseT>( mx::error_t::invalidarg, method + " is not a valid exclude method" );
    }
}

/// Image inclusion methods
/** \ingroup programming_library
 */
enum class include
{
    all,   ///< include all images
    corr,  ///< include images which are most correlated with the target
    time,  ///< include images which are closest in time to the target
    angle, ///< include images which are closest in angle to the target
    imno   ///< include images which are closest in image number to the target
};

/// Get the string name of an \ref include method
/**
 * \returns the string name corresponding to the \ref include
 *
 * \throws mx::exception with mx::error_t::invalidarg if the method is not valid (this is a bug)
 */
template <class verboseT>
std::string includeToStr( const include &method /**< [in] the \ref include method to convert to string */ )
{
    if( method == include::all )
    {
        return "all";
    }
    else if( method == include::corr )
    {
        return "corr";
    }
    else if( method == include::time )
    {
        return "time";
    }
    else if( method == include::angle )
    {
        return "angle";
    }
    else if( method == include::imno )
    {
        return "imno";
    }
    else
    {
        throw mx::exception<verboseT>( mx::error_t::invalidarg, "got an invalid include method (bug)" );
    }
}

/// Get the \ref include method from its string name
/**
 * \returns the \ref include corresponding to the name
 *
 * \throws mx::exception with mx::error_t::invalidarg if string name is not valid
 */

template <class verboseT>
include includeFmStr( const std::string &method /**< [in] the string name of the \ref include method */ )
{
    if( method == "all" )
    {
        return include::all;
    }
    else if( method == "corr" )
    {
        return include::corr;
    }
    else if( method == "time" )
    {
        return include::time;
    }
    else if( method == "angle" )
    {
        return include::angle;
    }
    else if( method == "imno" )
    {
        return include::imno;
    }
    else
    {
        throw mx::exception<verboseT>( mx::error_t::invalidarg, method + " is not a valid include method" );
    }
}

} // namespace HCI

} // namespace improc
} // namespace mx

#endif //__HCI_hpp__
