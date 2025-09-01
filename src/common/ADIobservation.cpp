/** \file ADIobservation.cpp
 * \author Jared R. Males
 * \brief Instantiates the ADI high contrast imaging data type.
 * \ingroup hc_imaging_files
 * \ingroup image_processing_files
 *
 */

#include "ADIDerotator.hpp"
#include "ADIobservation.hpp"

namespace mx
{
namespace improc
{
namespace HCI
{
std::string fakeMethodsStr( int method )
{
    if( method == single )
    {
        return "single";
    }
    else if( method == list )
    {
        return "list";
    }
    else
    {
        return "unknown";
    }
}

int fakeMethodFmStr( const std::string &method )
{
    if( method == "single" )
    {
        return single;
    }
    else if( method == "list" )
    {
        return list;
    }
    else
    {
        return -1;
    }
}
} // namespace HCI

template class ADIobservation<float, ADIDerotator<float, verbose::o>, verbose::o>;
template class ADIobservation<float, ADIDerotator<float, verbose::v>, verbose::v>;
template class ADIobservation<float, ADIDerotator<float, verbose::vv>, verbose::vv>;
template class ADIobservation<float, ADIDerotator<float, verbose::vvv>, verbose::vvv>;

template class ADIobservation<double, ADIDerotator<double, verbose::o>, verbose::o>;
template class ADIobservation<double, ADIDerotator<double, verbose::v>, verbose::v>;
template class ADIobservation<double, ADIDerotator<double, verbose::vv>, verbose::vv>;
template class ADIobservation<double, ADIDerotator<double, verbose::vvv>, verbose::vvv>;

} // namespace improc
} // namespace mx
