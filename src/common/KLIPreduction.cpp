/** \file KLIPreduction.cpp
 * \author Jared R. Males
 * \brief Instantiations of an implementation of the Karhunen-Loeve Image Processing (KLIP) algorithm.
 * \ingroup hc_imaging_files
 * \ingroup image_processing_files
 *
 */

#include "ADIDerotator.hpp"
#include "KLIPreduction.hpp"

namespace mx
{
namespace improc
{

template struct KLIPreduction<float, ADIDerotator<float, verbose::vv>, double, verbose::vv>;

} // namespace improc
} // namespace mx
