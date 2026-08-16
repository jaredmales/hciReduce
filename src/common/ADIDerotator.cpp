/** \file ADIDerotator.cpp
 * \author Jared R. Males
 * \brief Implements a generic ADI derotator class.
 *
 */

#include "ADIDerotator.hpp"

namespace mx
{
namespace improc
{

template struct ADIDerotator<float, mx::verbose::o>;
template struct ADIDerotator<float, mx::verbose::v>;
template struct ADIDerotator<float, mx::verbose::vv>;
template struct ADIDerotator<float, mx::verbose::vvv>;

template struct ADIDerotator<double, mx::verbose::o>;
template struct ADIDerotator<double, mx::verbose::v>;
template struct ADIDerotator<double, mx::verbose::vv>;
template struct ADIDerotator<double, mx::verbose::vvv>;

} // namespace improc
} // namespace mx
