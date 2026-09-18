/** \file P4PSFFilter.cpp
 * \brief Implements normalized local filtering and product naming for spatially variable PSF responses.
 */

#include "P4PSFFilter.hpp"

#include <cctype>
#include <cmath>
#include <filesystem>
#include <stdexcept>
#include <string_view>

namespace mx
{
namespace improc
{

namespace
{

/// Validate the common image geometry and support policy before either filter calculation.
void validatePSFFilterInputs( P4PSFFilter::imageConstRefT science,     /**< [in] final science image */
                              P4PSFFilter::imageConstRefT response,    /**< [in] signed response stamp */
                              P4PSFFilter::validityConstRefT validity, /**< [in] response validity */
                              int sourceRow,                           /**< [in] image row of the stamp center */
                              int sourceColumn,                        /**< [in] image column of the stamp center */
                              double minimumSupportFraction /**< [in] required usable fraction */ )
{
    if( science.rows() <= 0 || science.cols() <= 0 || response.rows() <= 0 || response.cols() <= 0 )
    {
        throw std::invalid_argument( "PSF filtering requires nonempty science and response images" );
    }
    if( response.rows() != validity.rows() || response.cols() != validity.cols() )
    {
        throw std::invalid_argument( "PSF response and validity dimensions must match" );
    }
    if( response.rows() % 2 == 0 || response.cols() % 2 == 0 )
    {
        throw std::invalid_argument( "PSF filtering requires odd response dimensions" );
    }
    if( sourceRow < 0 || sourceRow >= science.rows() || sourceColumn < 0 || sourceColumn >= science.cols() )
    {
        throw std::out_of_range( "PSF filter source coordinate is outside the science image" );
    }
    if( !std::isfinite( minimumSupportFraction ) || minimumSupportFraction < 0 || minimumSupportFraction > 1 )
    {
        throw std::invalid_argument( "PSF filter minimum support fraction must be finite and in [0,1]" );
    }
}

/// Normalize the filter and calculate statistics conditional on a fixed covariance and template.
P4PSFFilterResult normalizedPSFFilterResult( long double correlation,   /**< [in] weighted data-template correlation */
                                             long double normalization, /**< [in] weighted template norm */
                                             double supportFraction,    /**< [in] usable fraction of the full stamp */
                                             double minimumSupportFraction /**< [in] required usable fraction */ )
{
    P4PSFFilterResult result;
    result.supportFraction = supportFraction;
    result.correlation = static_cast<double>( correlation );
    result.normalization = static_cast<double>( normalization );
    if( result.supportFraction < minimumSupportFraction || result.normalization <= 0 ||
        !std::isfinite( result.correlation ) || !std::isfinite( result.normalization ) )
    {
        return result;
    }
    result.amplitude = static_cast<double>( correlation / normalization );
    const long double root = std::sqrt( normalization );
    result.conditionalSigma = static_cast<double>( 1 / root );
    result.score = static_cast<double>( correlation / root );
    result.nonnegativeLogLikelihood =
        correlation > 0 ? static_cast<double>( 0.5L * correlation / normalization * correlation ) : 0;
    result.valid = std::isfinite( result.amplitude ) && std::isfinite( result.conditionalSigma ) &&
                   std::isfinite( result.score ) && std::isfinite( result.nonnegativeLogLikelihood );
    return result;
}

} // namespace

std::string psfFilterProductPath( const std::string &finalImagePath, const std::string &role, bool sequential )
{
    const std::filesystem::path finalPath( finalImagePath );
    std::string fileName = finalPath.filename().string();
    if( fileName.empty() || role.empty() )
    {
        throw std::invalid_argument( "PSF filter-product naming requires a final-image filename and role" );
    }

    if( sequential )
    {
        constexpr std::size_t sequenceDigits = 4;
        constexpr std::string_view extension = ".fits";
        if( fileName.size() < sequenceDigits + extension.size() || !fileName.ends_with( extension ) )
        {
            throw std::invalid_argument( "sequential final-image path does not match the expected FITS naming" );
        }
        const std::size_t sequenceBegin = fileName.size() - extension.size() - sequenceDigits;
        for( std::size_t index = sequenceBegin; index < sequenceBegin + sequenceDigits; ++index )
        {
            if( !std::isdigit( static_cast<unsigned char>( fileName[index] ) ) )
            {
                throw std::invalid_argument( "sequential final-image path does not end in four digits" );
            }
        }
        if( sequenceBegin > 0 && !std::isalnum( static_cast<unsigned char>( fileName[sequenceBegin - 1] ) ) )
        {
            fileName.insert( sequenceBegin, role + fileName.substr( sequenceBegin - 1, 1 ) );
        }
        else
        {
            fileName.insert( sequenceBegin, "_" + role + "_" );
        }
    }
    else
    {
        const std::filesystem::path filePath( fileName );
        const std::string extension = filePath.extension().string();
        fileName = filePath.stem().string() + "_" + role + extension;
    }

    return ( finalPath.parent_path() / fileName ).string();
}

std::string psfFilterDiagnosticPath( const std::string &finalImagePath, const std::string &role, bool sequential )
{
    const std::filesystem::path finalPath( finalImagePath );
    const std::string finalStem = finalPath.stem().string();
    if( finalStem.empty() )
    {
        throw std::invalid_argument( "PSF filter-diagnostic naming requires a final-image filename" );
    }
    const std::filesystem::path productPath( psfFilterProductPath( finalImagePath, role, sequential ) );
    return ( finalPath.parent_path() / ( finalStem + "_outputs" ) / productPath.filename() ).string();
}

P4PSFFilterResult P4PSFFilter::calculate( imageConstRefT science,
                                          imageConstRefT response,
                                          validityConstRefT validity,
                                          int sourceRow,
                                          int sourceColumn,
                                          double minimumSupportFraction )
{
    validatePSFFilterInputs( science, response, validity, sourceRow, sourceColumn, minimumSupportFraction );

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
                throw std::invalid_argument( "valid PSF response samples must be finite" );
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

    return normalizedPSFFilterResult( correlation,
                                      normalization,
                                      static_cast<double>( usableSamples ) / static_cast<double>( response.size() ),
                                      minimumSupportFraction );
}

P4PSFFilterResult P4PSFFilter::calculateWeighted( imageConstRefT science,
                                                  imageConstRefT response,
                                                  validityConstRefT validity,
                                                  int sourceRow,
                                                  int sourceColumn,
                                                  double minimumSupportFraction,
                                                  const PSFNoiseModel &noise )
{
    validatePSFFilterInputs( science, response, validity, sourceRow, sourceColumn, minimumSupportFraction );
    if( noise.pixels() != response.size() )
    {
        throw std::invalid_argument( "noise covariance must match the full response-stamp pixel count" );
    }
    if( noise.kind() == PSFNoiseKind::identity )
    {
        return calculate( science, response, validity, sourceRow, sourceColumn, minimumSupportFraction );
    }

    const int centerRow = response.rows() / 2;
    const int centerColumn = response.cols() / 2;
    std::vector<Eigen::Index> support;
    support.reserve( static_cast<std::size_t>( response.size() ) );
    PSFNoiseModel::vectorT signal( response.size() );
    PSFNoiseModel::vectorT centered( response.size() );
    Eigen::Index retained{ 0 };
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
                throw std::invalid_argument( "valid PSF response samples must be finite" );
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
            const Eigen::Index pixel = stampRow + stampColumn * response.rows();
            support.push_back( pixel );
            signal( retained ) = responseValue;
            centered( retained ) = static_cast<double>( scienceValue ) - noise.mean()( pixel );
            ++retained;
        }
    }
    signal.conservativeResize( retained );
    centered.conservativeResize( retained );
    const PSFNoiseModel::vectorT weight = noise.solve( signal, support );
    long double correlation{ 0 };
    long double normalization{ 0 };
    for( Eigen::Index sample = 0; sample < retained; ++sample )
    {
        correlation += static_cast<long double>( weight( sample ) ) * centered( sample );
        normalization += static_cast<long double>( weight( sample ) ) * signal( sample );
    }
    return normalizedPSFFilterResult( correlation,
                                      normalization,
                                      static_cast<double>( retained ) / static_cast<double>( response.size() ),
                                      minimumSupportFraction );
}

} // namespace improc
} // namespace mx
