/** \file ADIobservation.cpp
 * \brief Instantiates the ADI high contrast imaging data type.
 *
 */

#include "ADIDerotator.hpp"
#include "ADIobservation.hpp"

namespace mx
{
namespace improc
{

template class ADIobservation<float, ADIDerotator<float, verbose::vv>, verbose::vv>;

} // namespace improc
} // namespace mx
