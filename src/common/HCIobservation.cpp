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

template class HCIobservation<float, mx::verbose::vv>;

} // namespace improc
} // namespace mx
