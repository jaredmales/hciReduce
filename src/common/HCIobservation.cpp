/** \file HCIobservation.cpp
 * \author Jared R. Males
 * \brief Instantiation of the basic high contrast imaging data type.
 * \ingroup hc_imaging_files
 * \ingroup image_processing_files
 *
 */

#include "HCIobservation.hpp"

namespace mx
{

namespace improc
{

namespace HCI
{

std::string coaddMethodStr( coaddMethod method )
{
    if( method == coaddMethod::none )
    {
        return "none";
    }
    else if( method == coaddMethod::median )
    {
        return "median";
    }
    else if( method == coaddMethod::mean )
    {
        return "mean";
    }
    else
    {
        mxThrowException(mx::err::invalidarg,"coaddMethodStr", "got an invalid coadd method (bug)");
    }

}

coaddMethod coaddMethodStr( const std::string &method )
{
    if( method == "none" )
    {
        return coaddMethod::none;
    }
    else if( method == "median" )
    {
        return coaddMethod::median;
    }
    else if( method == "mean" )
    {
        return coaddMethod::mean;
    }
    else
    {
        mxThrowException(mx::err::invalidarg,"coaddMethodStr", method + " is not a valid coadd method");
    }
}

std::string meanSubMethodStr( meanSubMethod method )
{
    if( method == meanSubMethod::none )
    {
        return "none";
    }
    else if( method == meanSubMethod::meanImage )
    {
        return "meanImage";
    }
    else if( method == meanSubMethod::medianImage )
    {
        return "medianImage";
    }
    else if( method == meanSubMethod::imageMean )
    {
        return "imageMean";
    }
    else if( method == meanSubMethod::imageMedian )
    {
        return "imageMedian";
    }
    else if( method == meanSubMethod::imageMode )
    {
        return "imageMode";
    }
    else
    {
        mxThrowException(mx::err::invalidarg, "meanSubMethodStr", "got an invalid mean sub method (bug)");
    }
}

meanSubMethod meanSubMethodStr( const std::string &method )
{
    if( method == "none" )
    {
        return meanSubMethod::none;
    }
    else if( method == "meanImage" )
    {
        return meanSubMethod::meanImage;
    }
    else if( method == "medianImage" )
    {
        return meanSubMethod::medianImage;
    }
    else if( method == "imageMean" )
    {
        return meanSubMethod::imageMean;
    }
    else if( method == "imageMedian" )
    {
        return meanSubMethod::imageMedian;
    }
    else if( method == "imageMode" )
    {
        return meanSubMethod::imageMode;
    }
    else
    {
        mxThrowException(mx::err::invalidarg, "meanSubMethodStr", method + " is not a valid mean sub method");
    }
}

std::string pixelTSNormMethodStr( pixelTSNormMethod method )
{
    if( method == pixelTSNormMethod::none )
    {
        return "none";
    }
    else if( method == pixelTSNormMethod::rms )
    {
        return "rms";
    }
    else if( method == pixelTSNormMethod::rmsSigmaClipped )
    {
        return "rmsSigmaClipped";
    }
    else
    {
        mxThrowException(mx::err::invalidarg,"pixelTSNormMethodStr", "got an invalid pixelTSNorm method (bug)");
    }
}

pixelTSNormMethod pixelTSNormMethodStr( const std::string &method )
{
    if( method == "none" )
    {
        return pixelTSNormMethod::none;
    }
    else if( method == "rms" )
    {
        return pixelTSNormMethod::rms;
    }
    else if( method == "rmsSigmaClipped" )
    {
        return pixelTSNormMethod::rmsSigmaClipped;
    }
    else
    {
        mxThrowException(mx::err::invalidarg,"pixelTSNormMethodStr", method + " is not a valid pixelTSNorm method");
    }
}

std::string combineMethodStr( combineMethod method )
{
    if( method == combineMethod::none )
    {
        return "none";
    }
    else if( method == combineMethod::median )
    {
        return "median";
    }
    else if( method == combineMethod::mean )
    {
        return "mean";
    }
    else if( method == combineMethod::sigmaMean )
    {
        return "sigmaMean";
    }
    else
    {
        mxThrowException(mx::err::invalidarg, "combineMethodStr", "got an invalid combine method (bug)");
    }

}

combineMethod combineMethodFmStr( const std::string &method )
{
    if( method == "none" )
    {
        return combineMethod::none;
    }
    else if( method == "median" )
    {
        return combineMethod::median;
    }
    else if( method == "mean" )
    {
        return combineMethod::mean;
    }
    else if( method == "sigmaMean" )
    {
        return combineMethod::sigmaMean;
    }
    else
    {
        mxThrowException(mx::err::invalidarg, "combineMethodFmStr", method + " is not a valid combine method");
    }
}

} // namespace HCI

template class HCIobservation<float>;
template class HCIobservation<double>;

} // namespace improc
} // namespace mx
