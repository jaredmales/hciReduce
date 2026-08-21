/** \file P4PSFFilter.cpp
 * \brief Implements normalized local filtering with spatially variable P4 PSFs.
 * \author Jared R. Males
 */

#include "P4PSFFilter.hpp"

#include <cmath>
#include <stdexcept>

namespace mx
{
namespace improc
{

P4PSFFilterResult P4PSFFilter::calculate( const imageT &science,
                                          const imageT &response,
                                          const validityT &validity,
                                          int sourceRow,
                                          int sourceColumn,
                                          double minimumSupportFraction )
{
    if( science.rows() <= 0 || science.cols() <= 0 || response.rows() <= 0 || response.cols() <= 0 )
    {
        throw std::invalid_argument( "P4 PSF filtering requires nonempty science and response images" );
    }
    if( response.rows() != validity.rows() || response.cols() != validity.cols() )
    {
        throw std::invalid_argument( "P4 PSF response and validity dimensions must match" );
    }
    if( response.rows() % 2 == 0 || response.cols() % 2 == 0 )
    {
        throw std::invalid_argument( "P4 PSF filtering requires odd response dimensions" );
    }
    if( sourceRow < 0 || sourceRow >= science.rows() || sourceColumn < 0 || sourceColumn >= science.cols() )
    {
        throw std::out_of_range( "P4 PSF filter source coordinate is outside the science image" );
    }
    if( !std::isfinite( minimumSupportFraction ) || minimumSupportFraction < 0 || minimumSupportFraction > 1 )
    {
        throw std::invalid_argument( "P4 PSF filter minimum support fraction must be finite and in [0,1]" );
    }

    const int centerRow = response.rows() / 2;
    const int centerColumn = response.cols() / 2;
    std::size_t usableSamples{ 0 };
    long double correlation{ 0 };
    long double normalization{ 0 };
    for( int stampColumn = 0; stampColumn < response.cols(); ++stampColumn )
    {
        const int imageColumn = sourceColumn + stampColumn - centerColumn;
        for( int stampRow = 0; stampRow < response.rows(); ++stampRow )
        {
            if( validity( stampRow, stampColumn ) == 0 )
            {
                continue;
            }
            const float responseValue = response( stampRow, stampColumn );
            if( !std::isfinite( responseValue ) )
            {
                throw std::invalid_argument( "valid P4 PSF response samples must be finite" );
            }
            const int imageRow = sourceRow + stampRow - centerRow;
            if( imageRow < 0 || imageRow >= science.rows() || imageColumn < 0 || imageColumn >= science.cols() )
            {
                continue;
            }
            const float scienceValue = science( imageRow, imageColumn );
            if( !std::isfinite( scienceValue ) )
            {
                continue;
            }
            ++usableSamples;
            correlation += static_cast<long double>( responseValue ) * scienceValue;
            normalization += static_cast<long double>( responseValue ) * responseValue;
        }
    }

    P4PSFFilterResult result;
    result.supportFraction = static_cast<double>( usableSamples ) / static_cast<double>( response.size() );
    result.correlation = static_cast<double>( correlation );
    result.normalization = static_cast<double>( normalization );
    if( result.supportFraction < minimumSupportFraction || normalization <= 0 || !std::isfinite( result.correlation ) ||
        !std::isfinite( result.normalization ) )
    {
        return result;
    }
    result.amplitude = static_cast<double>( correlation / normalization );
    result.valid = std::isfinite( result.amplitude );
    return result;
}

} // namespace improc
} // namespace mx
