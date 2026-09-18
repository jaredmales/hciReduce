/** \file PSFNoiseTraining.cpp
 * \brief Implements overlapping annular sampling and conditional covariance-weighted filtering.
 */

#include "PSFNoiseTraining.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace mx
{
namespace improc
{
namespace
{

/// Reject invalid training policies before accessing science or allocating sample matrices.
void validateNoiseTraining( const PSFNoiseTrainingConfig &config /**< [in] geometry and covariance settings */ )
{
    if( !std::isfinite( config.m_arcStep ) || config.m_arcStep < 0 || !std::isfinite( config.m_guardRadius ) ||
        config.m_guardRadius < 0 || config.m_minimumSamples < 2 || !std::isfinite( config.m_floorFraction ) ||
        config.m_floorFraction <= 0 || config.m_floorFraction > 1 ||
        ( config.m_kind != PSFNoiseKind::identity && config.m_kind != PSFNoiseKind::diagonal &&
          config.m_kind != PSFNoiseKind::pca ) )
    {
        throw std::invalid_argument(
            "noise training requires nonnegative spacing/guard, at least two samples, and a floor fraction in (0,1]" );
    }
}

/// Interpolate only from finite in-bounds pixels with nonzero bilinear weights.
bool sampleNoisePixel( double &value,                       /**< [out] interpolated value when usable */
                       P4PSFFilter::imageConstRefT science, /**< [in] unfiltered residual image */
                       double row,                          /**< [in] first-index image coordinate */
                       double column /**< [in] second-index image coordinate */ )
{
    if( row < 0 || row > science.rows() - 1 || column < 0 || column > science.cols() - 1 )
    {
        return false;
    }
    const int firstRow = static_cast<int>( std::floor( row ) );
    const int firstColumn = static_cast<int>( std::floor( column ) );
    const double rowFraction = row - firstRow;
    const double columnFraction = column - firstColumn;
    double pixels[2][2]{};
    for( int dc = 0; dc < 2; ++dc )
    {
        for( int dr = 0; dr < 2; ++dr )
        {
            const double weight =
                ( dr == 0 ? 1 - rowFraction : rowFraction ) * ( dc == 0 ? 1 - columnFraction : columnFraction );
            if( weight == 0 )
            {
                continue;
            }
            const float pixel = science( firstRow + dr, firstColumn + dc );
            if( !std::isfinite( pixel ) )
            {
                return false;
            }
            pixels[dr][dc] = pixel;
        }
    }
    value = std::lerp( std::lerp( pixels[0][0], pixels[1][0], rowFraction ),
                       std::lerp( pixels[0][1], pixels[1][1], rowFraction ),
                       columnFraction );
    return true;
}

} // namespace

PSFNoiseSamples PSFNoiseTraining::sample( P4PSFFilter::imageConstRefT science,
                                          int stampRows,
                                          int stampColumns,
                                          int sourceRow,
                                          int sourceColumn,
                                          const PSFNoiseTrainingConfig &config,
                                          const std::vector<PSFNoiseExclusion> &exclusions )
{
    validateNoiseTraining( config );
    if( science.rows() <= 0 || science.cols() <= 0 || stampRows <= 0 || stampColumns <= 0 || stampRows % 2 == 0 ||
        stampColumns % 2 == 0 || sourceRow < 0 || sourceRow >= science.rows() || sourceColumn < 0 ||
        sourceColumn >= science.cols() )
    {
        throw std::invalid_argument( "annular noise sampling requires odd stamps and an in-bounds candidate" );
    }
    for( const auto &exclusion : exclusions )
    {
        if( !std::isfinite( exclusion.m_row ) || !std::isfinite( exclusion.m_column ) ||
            !std::isfinite( exclusion.m_radius ) || exclusion.m_radius < 0 )
        {
            throw std::invalid_argument( "noise exclusion centers and nonnegative radii must be finite" );
        }
    }
    const double centerRow = 0.5 * ( science.rows() - 1 );
    const double centerColumn = 0.5 * ( science.cols() - 1 );
    const double radius = std::hypot( sourceRow - centerRow, sourceColumn - centerColumn );
    const double angle = std::atan2( sourceColumn - centerColumn, sourceRow - centerRow );
    const int halfRows = stampRows / 2;
    const int halfColumns = stampColumns / 2;
    const double stampRadius = std::hypot( halfRows, halfColumns );
    const double step = config.m_arcStep > 0 ? config.m_arcStep : std::max( 1, std::max( halfRows, halfColumns ) );
    const double count = std::ceil( 2 * std::numbers::pi * radius / step );
    const Eigen::Index pixels = static_cast<Eigen::Index>( stampRows ) * stampColumns;
    if( !std::isfinite( count ) || count >= static_cast<double>( std::numeric_limits<int>::max() ) )
    {
        throw std::invalid_argument( "requested annular noise sampling exceeds the supported sample count" );
    }
    PSFNoiseSamples result;
    result.m_attempted = static_cast<std::size_t>( count );
    result.m_arcStep = count > 0 ? 2 * std::numbers::pi * radius / count : 0;
    result.m_footprintRadius = stampRadius + std::sqrt( 2.0 );
    result.m_samples.resize( static_cast<Eigen::Index>( count ), pixels );
    result.m_centers.resize( static_cast<Eigen::Index>( count ), 2 );
    Eigen::Index retained = 0;
    for( std::size_t index = 0; index < result.m_attempted; ++index )
    {
        const double trainingAngle = 2 * std::numbers::pi * index / count;
        const double row = centerRow + radius * std::cos( trainingAngle );
        const double column = centerColumn + radius * std::sin( trainingAngle );
        bool excluded = std::hypot( row - sourceRow, column - sourceColumn ) <=
                        result.m_footprintRadius + stampRadius + config.m_guardRadius;
        for( const auto &exclusion : exclusions )
        {
            excluded = excluded || std::hypot( row - exclusion.m_row, column - exclusion.m_column ) <=
                                       result.m_footprintRadius + exclusion.m_radius;
        }
        if( excluded )
        {
            ++result.m_excluded;
            continue;
        }
        const double cosine = std::cos( trainingAngle - angle );
        const double sine = std::sin( trainingAngle - angle );
        bool complete = true;
        for( int sc = 0; sc < stampColumns && complete; ++sc )
        {
            for( int sr = 0; sr < stampRows; ++sr )
            {
                const double dr = sr - halfRows;
                const double dc = sc - halfColumns;
                double value;
                if( !sampleNoisePixel( value,
                                       science,
                                       row + cosine * dr - sine * dc,
                                       column + sine * dr + cosine * dc ) )
                {
                    complete = false;
                    break;
                }
                result.m_samples( retained, sr + static_cast<Eigen::Index>( sc ) * stampRows ) = value;
            }
        }
        if( !complete )
        {
            ++result.m_incomplete;
            continue;
        }
        result.m_centers( retained, 0 ) = row;
        result.m_centers( retained, 1 ) = column;
        ++retained;
    }
    result.m_samples.conservativeResize( retained, pixels );
    result.m_centers.conservativeResize( retained, 2 );
    return result;
}

PSFNoiseTrainingResult PSFNoiseTraining::calculate( P4PSFFilter::imageConstRefT science,
                                                    P4PSFFilter::imageConstRefT response,
                                                    P4PSFFilter::validityConstRefT validity,
                                                    int sourceRow,
                                                    int sourceColumn,
                                                    double minimumSupportFraction,
                                                    const PSFNoiseTrainingConfig &config,
                                                    const std::vector<PSFNoiseExclusion> &exclusions )
{
    validateNoiseTraining( config );
    // Validate the shared filter contract even when training will be insufficient.
    PSFNoiseTrainingResult result;
    result.m_filter =
        P4PSFFilter::calculate( science, response, validity, sourceRow, sourceColumn, minimumSupportFraction );
    if( !result.m_filter.valid )
    {
        result.m_status = PSFNoiseTrainingStatus::invalidFilter;
        return result;
    }
    if( config.m_kind == PSFNoiseKind::identity )
    {
        result.m_varianceFloor = 1;
        result.m_status = result.m_filter.valid ? PSFNoiseTrainingStatus::valid : PSFNoiseTrainingStatus::invalidFilter;
        return result;
    }
    result.m_filter = P4PSFFilterResult{};
    const PSFNoiseSamples training =
        sample( science, response.rows(), response.cols(), sourceRow, sourceColumn, config, exclusions );
    result.m_samples = static_cast<std::size_t>( training.m_samples.rows() );
    result.m_attempted = training.m_attempted;
    result.m_excluded = training.m_excluded;
    result.m_incomplete = training.m_incomplete;
    result.m_arcStep = training.m_arcStep;
    if( result.m_samples < config.m_minimumSamples )
    {
        return result;
    }
    const PSFNoiseModel::matrixT centered = training.m_samples.rowwise() - training.m_samples.colwise().mean();
    PSFNoiseModel::vectorT variances = centered.colwise().squaredNorm().transpose() / ( result.m_samples - 1 );
    if( !variances.allFinite() )
    {
        throw std::overflow_error( "annular noise sample variance is nonfinite" );
    }
    std::sort( variances.data(), variances.data() + variances.size() );
    // Odd stamp dimensions give an odd pixel count and an unambiguous median.
    result.m_varianceFloor = config.m_floorFraction * variances( variances.size() / 2 );
    if( !std::isfinite( result.m_varianceFloor ) || result.m_varianceFloor <= 0 )
    {
        result.m_status = PSFNoiseTrainingStatus::zeroVariance;
        return result;
    }
    const PSFNoiseModel model =
        config.m_kind == PSFNoiseKind::diagonal
            ? PSFNoiseModel::estimateDiagonal( training.m_samples, result.m_varianceFloor )
            : PSFNoiseModel::estimatePCA( training.m_samples, config.m_maximumModes, result.m_varianceFloor );
    result.m_modes = model.modes();
    result.m_filter = P4PSFFilter::calculateWeighted( science,
                                                      response,
                                                      validity,
                                                      sourceRow,
                                                      sourceColumn,
                                                      minimumSupportFraction,
                                                      model );
    result.m_status = result.m_filter.valid ? PSFNoiseTrainingStatus::valid : PSFNoiseTrainingStatus::invalidFilter;
    return result;
}

} // namespace improc
} // namespace mx
