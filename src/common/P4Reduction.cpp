/** \file P4Reduction.cpp
 * \brief Implements the observation orchestrator for Pixel Prediction Post-Processing.
 * \author Jared R. Males
 */

#include "ADIDerotator.hpp"
#include "P4Reduction.hpp"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>

#include <omp.h>

#include <mx/ioutils/fileUtils.hpp>
#include <mx/ioutils/fits/fitsFile.hpp>
#include <mx/ipc/ompLoopWatcher.hpp>
#include <mx/math/floatUtils.hpp>
#include <mx/math/geo.hpp>
#include <mx/sys/timeUtils.hpp>

namespace mx
{
namespace improc
{

namespace
{

constexpr int p4AutomaticCropPadding{ 4 };

/// Derive a filter-product path from the resolved final-image path.
std::string p4FilterProductPath( const std::string &finalImagePath, /**< [in] resolved final-image output path */
                                 const std::string &role, /**< [in] filename role inserted before the sequence */
                                 bool sequential /**< [in] whether the final image uses a four-digit sequence */ )
{
    const std::filesystem::path finalPath( finalImagePath );
    std::string fileName = finalPath.filename().string();
    if( fileName.empty() || role.empty() )
    {
        throw std::invalid_argument( "P4 filter-product naming requires a final-image filename and role" );
    }

    if( sequential )
    {
        constexpr std::size_t sequenceDigits = 4;
        constexpr std::string_view extension = ".fits";
        if( fileName.size() < sequenceDigits + extension.size() || !fileName.ends_with( extension ) )
        {
            throw std::invalid_argument( "P4 sequential final-image path does not match the expected FITS naming" );
        }
        const std::size_t sequenceBegin = fileName.size() - extension.size() - sequenceDigits;
        for( std::size_t index = sequenceBegin; index < sequenceBegin + sequenceDigits; ++index )
        {
            if( !std::isdigit( static_cast<unsigned char>( fileName[index] ) ) )
            {
                throw std::invalid_argument( "P4 sequential final-image path does not end in four digits" );
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

/// Derive the P4 auxiliary-product directory from the resolved final-image path.
std::filesystem::path
p4AuxiliaryProductDirectory( const std::string &finalImagePath /**< [in] resolved final-image output path */ )
{
    const std::filesystem::path finalPath( finalImagePath );
    const std::string finalStem = finalPath.stem().string();
    if( finalStem.empty() )
    {
        throw std::invalid_argument( "P4 auxiliary-product naming requires a final-image filename" );
    }

    return finalPath.parent_path() / ( finalStem + "_outputs" );
}

/// Derive a filter-diagnostic path inside the final image's P4 PSF product directory.
std::string p4FilterDiagnosticPath( const std::string &finalImagePath, /**< [in] resolved final-image output path */
                                    const std::string &role, /**< [in] filename role inserted before the sequence */
                                    bool sequential /**< [in] whether the final image uses a four-digit sequence */ )
{
    const std::filesystem::path productPath( p4FilterProductPath( finalImagePath, role, sequential ) );
    return ( p4AuxiliaryProductDirectory( finalImagePath ) / productPath.filename() ).string();
}

/// Join a vector into a stable comma-delimited FITS value.
template <typename valueT>
std::string p4Join( const std::vector<valueT> &values /**< [in] ordered values */ )
{
    std::ostringstream stream;
    stream << std::setprecision( std::numeric_limits<long double>::digits10 + 1 );
    for( std::size_t index = 0; index < values.size(); ++index )
    {
        if( index != 0 )
        {
            stream << ',';
        }
        stream << values[index];
    }
    return stream.str();
}

/// Construct a stable zero-padded diagnostic or FITS-card index.
std::string p4Index( std::size_t index, /**< [in] zero-based index */
                     int width = 3 /**< [in] minimum decimal field width */ )
{
    std::ostringstream stream;
    stream << std::setfill( '0' ) << std::setw( width ) << index;
    return stream.str();
}

/// Select central target images and their qualifying temporal predictor images for one detector-frame annulus.
std::vector<std::vector<int>> p4TemporalSelections( const std::vector<double> &derotationAngles,
                                                    /**< [in] derotation angles in input-image order, in radians */
                                                    const std::vector<int> &frameIndices,
                                                    /**< [in] ordered physical frame indices eligible for selection */
                                                    double minimumRadius, /**< [in] annulus inner radius in pixels */
                                                    double psfRadius,     /**< [in] physical PSF radius in pixels */
                                                    int numberImages /**< [in] required qualifying images per side */ )
{
    std::vector<std::vector<int>> selections;
    for( std::size_t centralOffset = 0; centralOffset < frameIndices.size(); ++centralOffset )
    {
        const int central = frameIndices[centralOffset];
        std::vector<int> earlier;
        std::vector<int> later;
        for( std::size_t candidateOffset = centralOffset; candidateOffset-- > 0; )
        {
            const int candidate = frameIndices[candidateOffset];
            const double displacement =
                std::abs( mx::math::angleDiff<mx::math::radiansT<double>>( derotationAngles[candidate],
                                                                           derotationAngles[central] ) ) *
                minimumRadius;
            if( displacement >= psfRadius )
            {
                earlier.push_back( candidate );
            }
        }
        for( std::size_t candidateOffset = centralOffset + 1; candidateOffset < frameIndices.size(); ++candidateOffset )
        {
            const int candidate = frameIndices[candidateOffset];
            const double displacement =
                std::abs( mx::math::angleDiff<mx::math::radiansT<double>>( derotationAngles[candidate],
                                                                           derotationAngles[central] ) ) *
                minimumRadius;
            if( displacement >= psfRadius )
            {
                later.push_back( candidate );
            }
        }
        const std::size_t requiredPerSide = static_cast<std::size_t>( numberImages );
        const std::size_t earlierCount = std::min( earlier.size(), requiredPerSide );
        const std::size_t laterCount = std::min( later.size(), requiredPerSide );
        const std::size_t requiredTotal = 2 * requiredPerSide;
        const std::size_t additionalEarlier =
            std::min( earlier.size() - earlierCount, requiredTotal - earlierCount - laterCount );
        const std::size_t additionalLater = requiredTotal - earlierCount - laterCount - additionalEarlier;
        if( additionalLater > later.size() - laterCount )
        {
            continue;
        }
        std::vector<int> selection{ central };
        selection.insert( selection.end(), earlier.begin(), earlier.begin() + earlierCount + additionalEarlier );
        selection.insert( selection.end(), later.begin(), later.begin() + laterCount + additionalLater );
        selections.push_back( std::move( selection ) );
    }
    return selections;
}

/// Report whether one temporal selection can realize every configured uncentered detector-frame PCA mode.
bool p4TemporalSelectionSupportsModes( const std::vector<std::vector<int>> &selections,
                                       /**< [in] retained central targets and their temporal predictor images */
                                       std::size_t basePredictorCount, /**< [in] same-image OR predictor count */
                                       std::size_t temporalPredictorCount,
                                       /**< [in] physical-PSF predictor pixels per neighboring image */
                                       const std::vector<float> &modeFractions /**< [in] requested PCA fractions */ )
{
    if( selections.empty() || selections.front().empty() ||
        selections.front().size() - 1 > std::numeric_limits<std::size_t>::max() / temporalPredictorCount )
    {
        return false;
    }
    const std::size_t neighborPredictorCount = temporalPredictorCount * ( selections.front().size() - 1 );
    if( basePredictorCount > std::numeric_limits<std::size_t>::max() - neighborPredictorCount )
    {
        return false;
    }
    const std::size_t predictorCount = basePredictorCount + neighborPredictorCount;
    if( predictorCount > static_cast<std::size_t>( std::numeric_limits<int>::max() ) ||
        selections.size() > static_cast<std::size_t>( std::numeric_limits<int>::max() ) )
    {
        return false;
    }
    const int maximumDegreesOfFreedom =
        std::min( static_cast<int>( selections.size() ), static_cast<int>( predictorCount ) );
    int previousMode{ 0 };
    for( const float fraction : modeFractions )
    {
        const int mode = static_cast<int>( std::floor( static_cast<double>( fraction ) * maximumDegreesOfFreedom ) );
        if( mode <= 0 || mode <= previousMode )
        {
            return false;
        }
        previousMode = mode;
    }
    return true;
}

/// Select the greatest positive temporal exclusion radius that retains structurally usable detector-frame PCA rows.
std::pair<double, std::vector<std::vector<int>>>
p4TemporalSelectionWithFallback( const std::vector<double> &derotationAngles,
                                 /**< [in] derotation angles in input-image order, in radians */
                                 const std::vector<int> &frameIndices,
                                 /**< [in] ordered physical frame indices eligible for selection */
                                 double meanRadius,              /**< [in] mean annulus radius in pixels */
                                 double requestedPsfRadius,      /**< [in] configured physical PSF radius in pixels */
                                 int numberImages,               /**< [in] requested qualifying images per side */
                                 std::size_t basePredictorCount, /**< [in] same-image OR predictor count */
                                 std::size_t temporalPredictorCount,
                                 /**< [in] physical-PSF predictor pixels per neighboring image */
                                 const std::vector<float> &modeFractions /**< [in] requested PCA fractions */ )
{
    std::vector<double> candidateRadii{ requestedPsfRadius };
    for( std::size_t firstOffset = 0; firstOffset < frameIndices.size(); ++firstOffset )
    {
        const int first = frameIndices[firstOffset];
        for( std::size_t secondOffset = firstOffset + 1; secondOffset < frameIndices.size(); ++secondOffset )
        {
            const int second = frameIndices[secondOffset];
            const double radius =
                std::abs( mx::math::angleDiff<mx::math::radiansT<double>>( derotationAngles[first],
                                                                           derotationAngles[second] ) ) *
                meanRadius;
            if( radius > 0 && radius < requestedPsfRadius )
            {
                candidateRadii.push_back( radius );
            }
        }
    }
    std::sort( candidateRadii.begin(), candidateRadii.end(), std::greater<double>() );
    candidateRadii.erase( std::unique( candidateRadii.begin(), candidateRadii.end() ), candidateRadii.end() );
    for( const double radius : candidateRadii )
    {
        std::vector<std::vector<int>> selections =
            p4TemporalSelections( derotationAngles, frameIndices, meanRadius, radius, numberImages );
        if( p4TemporalSelectionSupportsModes( selections, basePredictorCount, temporalPredictorCount, modeFractions ) )
        {
            return { radius, std::move( selections ) };
        }
    }
    return { 0, {} };
}

/// Enumerate direct detector-pixel offsets inside the physical temporal PSF predictor disk.
std::vector<P4PixelCoordinate> p4TemporalPredictorOffsets( double psfRadius /**< [in] physical PSF radius in pixels */ )
{
    const int maximumOffset = static_cast<int>( std::ceil( psfRadius ) );
    std::vector<P4PixelCoordinate> offsets;
    for( int row = -maximumOffset; row <= maximumOffset; ++row )
    {
        for( int column = -maximumOffset; column <= maximumOffset; ++column )
        {
            if( std::hypot( static_cast<double>( row ), static_cast<double>( column ) ) <= psfRadius )
            {
                offsets.emplace_back( row, column );
            }
        }
    }
    return offsets;
}

/// Report whether every direct temporal PSF predictor pixel is inside the image and common mask.
bool p4TemporalPredictorsValid( const P4PixelCoordinate &searchCoordinate,
                                /**< [in] central search-pixel coordinate */
                                const std::vector<P4PixelCoordinate> &offsets,
                                /**< [in] physical temporal PSF predictor offsets */
                                int rows,    /**< [in] detector-image row count */
                                int columns, /**< [in] detector-image column count */
                                const P4PixelGridf::imageT *mask /**< [in] optional common mask */ )
{
    for( const P4PixelCoordinate &offset : offsets )
    {
        const int row = searchCoordinate.row() + offset.row();
        const int column = searchCoordinate.column() + offset.column();
        if( row < 0 || row >= rows || column < 0 || column >= columns || ( mask && ( *mask )( row, column ) == 0 ) )
        {
            return false;
        }
    }
    return true;
}

/** \cond P4ProgressOutput */

/// Decorate one annulus-local OpenMP loop watcher with overall P4 status.
/** The object is used only through the serialized output operations performed by `mx::ipc::ompLoopWatcher` while a
 * region is running. Region transitions occur outside the OpenMP parallel region.
 */
class P4ProgressOutput
{
  private:
    std::ostream &m_output;                ///< Destination stream owned by the caller.

    std::ostringstream m_status;           ///< One status assembled by `mx::ipc::ompLoopWatcher`.

    std::size_t m_totalRegions{ 0 };       ///< Total configured annulus count.

    std::size_t m_totalSearchPixels{ 0 };  ///< Total search-pixel count across all annuli.

    bool m_rotated{ false };               ///< Whether status describes rotated-frame regression.

    std::size_t m_region{ 0 };             ///< Zero-based active annulus index.

    std::size_t m_regionSearchPixels{ 0 }; ///< Search-pixel count in the active annulus.

    std::size_t m_predictorCount{ 0 };     ///< Predictor count in the active annulus.

    std::size_t m_regionStart{ 0 };        ///< Overall completed count at the active annulus start.

    std::size_t m_regionCompleted{ 0 };    ///< Active-annulus counter captured from the watcher status.

    double m_minimumRadius{ 0 };           ///< Inclusive active-annulus inner radius.

    double m_maximumRadius{ 0 };           ///< Exclusive active-annulus outer radius.

    bool m_captureRegionCounter{ true };   ///< Whether the next size value is the watcher counter.

    /// Return a finite completion percentage for a nonempty reduction.
    double percentage( std::size_t completed /**< [in] completed overall search pixels */ ) const
    {
        return 100.0 * static_cast<double>( completed ) / static_cast<double>( m_totalSearchPixels );
    }

    /// Write the stable active-annulus label.
    void writePrefix()
    {
        m_output << "P4 ";
        if( m_rotated )
        {
            m_output << "rotated ";
        }
        m_output << "annulus " << m_region + 1 << " / " << m_totalRegions << " [" << m_minimumRadius << ','
                 << m_maximumRadius << "), K=" << m_predictorCount << ": ";
    }

    /// Publish one complete watcher status and reset the assembly buffer.
    void publish()
    {
        std::string status = m_status.str();
        if( !status.empty() && status.back() == '\n' )
        {
            status.pop_back();
        }
        const std::size_t overallCompleted = m_regionStart + m_regionCompleted;
        writePrefix();
        m_output << status << " | overall " << overallCompleted << " / " << m_totalSearchPixels << " ("
                 << percentage( overallCompleted ) << "%)           \r" << std::flush;
        m_status.str( std::string() );
        m_status.clear();
        m_captureRegionCounter = true;
    }

    /// Clear a conservative terminal-width status line before writing a shorter record.
    void clearLine()
    {
        m_output << '\r' << std::string( 160, ' ' ) << '\r';
    }

  public:
    /// Construct a decorator for one complete P4 regression.
    P4ProgressOutput( std::ostream &output,          /**< [in,out] destination stream */
                      std::size_t totalRegions,      /**< [in] total configured annuli */
                      std::size_t totalSearchPixels, /**< [in] total search pixels */
                      bool rotated )                 /**< [in] whether regression is in the rotated frame */
        : m_output( output ), m_totalRegions( totalRegions ), m_totalSearchPixels( totalSearchPixels ),
          m_rotated( rotated )
    {
    }

    /// Begin reporting one active annulus and emit its initial zero-progress status.
    void beginRegion( std::size_t region,             /**< [in] zero-based annulus index */
                      std::size_t regionSearchPixels, /**< [in] active-annulus search pixels */
                      std::size_t predictorCount,     /**< [in] active-annulus predictors */
                      std::size_t overallCompleted,   /**< [in] completed prior-annulus pixels */
                      double minimumRadius,           /**< [in] inclusive annulus inner radius */
                      double maximumRadius )          /**< [in] exclusive annulus outer radius */
    {
        m_region = region;
        m_regionSearchPixels = regionSearchPixels;
        m_predictorCount = predictorCount;
        m_regionStart = overallCompleted;
        m_regionCompleted = 0;
        m_minimumRadius = minimumRadius;
        m_maximumRadius = maximumRadius;
        writePrefix();
        m_output << "0 / " << m_regionSearchPixels << " (0%) | overall " << overallCompleted << " / "
                 << m_totalSearchPixels << " (" << percentage( overallCompleted ) << "%)           \r" << std::flush;
    }

    /// Finish one annulus with exact local and overall counts and a newline.
    void completeRegion( std::size_t overallCompleted /**< [in] completed overall search pixels */ )
    {
        clearLine();
        writePrefix();
        m_output << m_regionSearchPixels << " / " << m_regionSearchPixels << " (100%) complete | overall "
                 << overallCompleted << " / " << m_totalSearchPixels << " (" << percentage( overallCompleted ) << "%)\n"
                 << std::flush;
    }

    /// Clear the in-place status and report interruption before an exception is rethrown.
    void interruptRegion()
    {
        clearLine();
        writePrefix();
        m_output << "interrupted\n" << std::flush;
    }

    /// Report the final overall regression time.
    void completeReduction( double seconds /**< [in] elapsed regression seconds */ )
    {
        m_output << "P4 regression complete: " << m_totalSearchPixels << " / " << m_totalSearchPixels
                 << " search pixels in " << seconds << " s\n"
                 << std::flush;
    }

    /// Capture a watcher size value, treating the first one as its annulus-local counter.
    P4ProgressOutput &operator<<( std::size_t value /**< [in] watcher counter or total */ )
    {
        if( m_captureRegionCounter )
        {
            m_regionCompleted = value;
            m_captureRegionCounter = false;
        }
        m_status << value;
        return *this;
    }

    /// Capture any non-size watcher status component.
    template <typename valueT>
    P4ProgressOutput &operator<<( const valueT &value /**< [in] status component */ )
    {
        m_status << value;
        return *this;
    }

    /// Publish the assembled status when the watcher flushes its output.
    P4ProgressOutput &operator<<( std::ostream &( *manipulator )(std::ostream &)/**< [in] stream manipulator */ )
    {
        static_cast<void>( manipulator );
        publish();
        return *this;
    }
};

/** \endcond */

#ifdef HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION
/** \cond P4Reduction_precision_experiment */

/// Experimental direct-PCA dispatch signature shared by every reduction worker.
using P4ReductionDirectFunction =
    void ( * )( P4PCAResult &output,                                 /**< [out] completed regression result */
                const P4PCA::matrixT &predictors,                    /**< [in] sampled predictor matrix */
                const P4PCA::vectorT &target,                        /**< [in] sampled target series */
                const std::vector<int> &modes,                       /**< [in] requested retained-mode counts */
                double rankTolerance,                                /**< [in] relative numerical-rank threshold */
                P4PCA::workspaceT &doubleWorkspace,                  /**< [in,out] ordinary FP64 workspace */
                detail::P4PCAExperimentalWorkspace &policyWorkspace, /**< [in,out] typed experimental workspace */
                P4PCATiming *timing,                                 /**< [out] optional stage timing */
                P4PCA::matrixT *coefficients /**< [out] optional predictor coefficients */ );

/// Experimental centered-PCA dispatch signature shared by every reduction worker.
using P4ReductionCenteredFunction =
    void ( * )( P4PCAResult &output,                /**< [out] completed centered regression result */
                P4PCA::matrixT &predictors,         /**< [in,out] sampled predictors, centered in place */
                const P4PCA::vectorT &target,       /**< [in] sampled target series */
                const std::vector<int> &modes,      /**< [in] requested retained-mode counts */
                double rankTolerance,               /**< [in] relative numerical-rank threshold */
                P4PCA::workspaceT &doubleWorkspace, /**< [in,out] ordinary FP64 workspace */
                detail::P4PCAExperimentalWorkspace &policyWorkspace, /**< [in,out] typed experimental workspace */
                P4PCATiming *timing,                                 /**< [out] optional stage timing */
                P4PCA::matrixT *coefficients /**< [out] optional predictor coefficients */ );

/// Experimental held-out PCA dispatch signature shared by every reduction worker.
using P4ReductionHeldOutFunction =
    void ( * )( P4PCAResult &output,                                 /**< [out] completed held-out regression result */
                const P4PCA::matrixT &predictors,                    /**< [in] sampled predictor matrix */
                const P4PCA::vectorT &target,                        /**< [in] sampled target series */
                const P4TargetExclusions &exclusions,                /**< [in] target-specific deleted rows */
                const std::vector<int> &modes,                       /**< [in] requested retained-mode counts */
                double rankTolerance,                                /**< [in] relative numerical-rank threshold */
                P4PCA::workspaceT &doubleWorkspace,                  /**< [in,out] ordinary FP64 workspace */
                detail::P4PCAExperimentalWorkspace &policyWorkspace, /**< [in,out] typed experimental workspace */
                P4PCATiming *timing /**< [out] optional stage timing */ );

/// Experimental held-out frozen-probe dispatch signature shared by every reduction worker.
using P4ReductionHeldOutProbeFunction =
    void ( * )( P4PCAResult &output,                                 /**< [out] completed held-out regression result */
                P4PCA::matrixT &probeResiduals,                      /**< [out] frozen-probe residual responses */
                const P4PCA::matrixT &predictors,                    /**< [in] sampled predictor matrix */
                const P4PCA::vectorT &target,                        /**< [in] sampled target series */
                const P4PCA::matrixT &probePredictors,               /**< [in] frozen-probe predictor responses */
                const P4PCA::vectorT &probeTarget,                   /**< [in] frozen-probe direct response */
                const P4TargetExclusions &exclusions,                /**< [in] target-specific deleted rows */
                const std::vector<int> &modes,                       /**< [in] requested retained-mode counts */
                double rankTolerance,                                /**< [in] relative numerical-rank threshold */
                P4PCA::workspaceT &doubleWorkspace,                  /**< [in,out] ordinary FP64 workspace */
                detail::P4PCAExperimentalWorkspace &policyWorkspace, /**< [in,out] typed experimental workspace */
                P4PCATiming *timing /**< [out] optional stage timing */ );

/// Immutable PCA call targets selected before entering a reduction worker loop.
struct P4ReductionPCADispatch
{
    P4ReductionDirectFunction direct;             ///< Direct detector-frame evaluator.

    P4ReductionCenteredFunction centeredInPlace;  ///< Destructive centered rotated-frame evaluator.

    P4ReductionHeldOutFunction heldOut;           ///< Explicit held-out evaluator.

    P4ReductionHeldOutProbeFunction heldOutProbe; ///< Explicit held-out frozen-probe evaluator.
};

/// Invoke one direct PCA policy without worker-loop policy selection.
template <detail::P4PCAPrecisionPolicy precisionPolicy>
void p4ReductionCalculateDirect(
    P4PCAResult &output,                                       /**< [out] completed regression result */
    const P4PCA::matrixT &predictors,                          /**< [in] sampled predictor matrix */
    const P4PCA::vectorT &target,                              /**< [in] sampled target series */
    const std::vector<int> &modes,                             /**< [in] requested retained-mode counts */
    double rankTolerance,                                      /**< [in] relative numerical-rank threshold */
    P4PCA::workspaceT &doubleWorkspace,                        /**< [in,out] ordinary FP64 workspace */
    detail::P4PCAExperimentalWorkspace &experimentalWorkspace, /**< [in,out] typed experimental workspace */
    P4PCATiming *timing,                                       /**< [out] optional stage timing */
    P4PCA::matrixT *coefficients /**< [out] optional predictor coefficients */ )
{
    if constexpr( precisionPolicy == detail::P4PCAPrecisionPolicy::doubleDouble )
    {
        static_cast<void>( experimentalWorkspace );
        P4PCA::calculate( output, predictors, target, modes, rankTolerance, doubleWorkspace, timing, coefficients );
    }
    else
    {
        static_cast<void>( doubleWorkspace );
        detail::p4PCACalculateExperimental( output,
                                            predictors,
                                            target,
                                            modes,
                                            rankTolerance,
                                            precisionPolicy,
                                            experimentalWorkspace,
                                            timing,
                                            coefficients );
    }
}

/// Invoke one centered PCA policy without worker-loop policy selection.
template <detail::P4PCAPrecisionPolicy precisionPolicy>
void p4ReductionCalculateCenteredInPlace(
    P4PCAResult &output,                                       /**< [out] completed centered regression result */
    P4PCA::matrixT &predictors,                                /**< [in,out] predictors, centered in place */
    const P4PCA::vectorT &target,                              /**< [in] sampled target series */
    const std::vector<int> &modes,                             /**< [in] requested retained-mode counts */
    double rankTolerance,                                      /**< [in] relative numerical-rank threshold */
    P4PCA::workspaceT &doubleWorkspace,                        /**< [in,out] ordinary FP64 workspace */
    detail::P4PCAExperimentalWorkspace &experimentalWorkspace, /**< [in,out] typed experimental workspace */
    P4PCATiming *timing,                                       /**< [out] optional stage timing */
    P4PCA::matrixT *coefficients /**< [out] optional predictor coefficients */ )
{
    if constexpr( precisionPolicy == detail::P4PCAPrecisionPolicy::doubleDouble )
    {
        static_cast<void>( experimentalWorkspace );
        P4PCA::calculateCenteredInPlace( output,
                                         predictors,
                                         target,
                                         modes,
                                         rankTolerance,
                                         doubleWorkspace,
                                         timing,
                                         coefficients );
    }
    else
    {
        static_cast<void>( doubleWorkspace );
        detail::p4PCACalculateCenteredInPlaceExperimental( output,
                                                           predictors,
                                                           target,
                                                           modes,
                                                           rankTolerance,
                                                           precisionPolicy,
                                                           experimentalWorkspace,
                                                           timing,
                                                           coefficients );
    }
}

/// Invoke one explicit held-out PCA policy without worker-loop policy selection.
template <detail::P4PCAPrecisionPolicy precisionPolicy>
void p4ReductionCalculateHeldOut(
    P4PCAResult &output,                                       /**< [out] completed held-out regression result */
    const P4PCA::matrixT &predictors,                          /**< [in] sampled predictor matrix */
    const P4PCA::vectorT &target,                              /**< [in] sampled target series */
    const P4TargetExclusions &exclusions,                      /**< [in] target-specific deleted rows */
    const std::vector<int> &modes,                             /**< [in] requested retained-mode counts */
    double rankTolerance,                                      /**< [in] relative numerical-rank threshold */
    P4PCA::workspaceT &doubleWorkspace,                        /**< [in,out] ordinary FP64 workspace */
    detail::P4PCAExperimentalWorkspace &experimentalWorkspace, /**< [in,out] typed experimental workspace */
    P4PCATiming *timing /**< [out] optional stage timing */ )
{
    if constexpr( precisionPolicy == detail::P4PCAPrecisionPolicy::doubleDouble )
    {
        static_cast<void>( experimentalWorkspace );
        P4PCA::calculateHeldOut( output,
                                 predictors,
                                 target,
                                 exclusions,
                                 modes,
                                 rankTolerance,
                                 doubleWorkspace,
                                 timing );
    }
    else
    {
        static_cast<void>( doubleWorkspace );
        detail::p4PCACalculateHeldOutExperimental( output,
                                                   predictors,
                                                   target,
                                                   exclusions,
                                                   modes,
                                                   rankTolerance,
                                                   precisionPolicy,
                                                   experimentalWorkspace,
                                                   timing );
    }
}

/// Invoke one explicit held-out frozen-probe policy without worker-loop policy selection.
template <detail::P4PCAPrecisionPolicy precisionPolicy>
void p4ReductionCalculateHeldOutProbe(
    P4PCAResult &output,                                       /**< [out] completed held-out regression result */
    P4PCA::matrixT &probeResiduals,                            /**< [out] frozen-probe residual responses */
    const P4PCA::matrixT &predictors,                          /**< [in] sampled predictor matrix */
    const P4PCA::vectorT &target,                              /**< [in] sampled target series */
    const P4PCA::matrixT &probePredictors,                     /**< [in] frozen-probe predictor responses */
    const P4PCA::vectorT &probeTarget,                         /**< [in] frozen-probe direct response */
    const P4TargetExclusions &exclusions,                      /**< [in] target-specific deleted rows */
    const std::vector<int> &modes,                             /**< [in] requested retained-mode counts */
    double rankTolerance,                                      /**< [in] relative numerical-rank threshold */
    P4PCA::workspaceT &doubleWorkspace,                        /**< [in,out] ordinary FP64 workspace */
    detail::P4PCAExperimentalWorkspace &experimentalWorkspace, /**< [in,out] typed experimental workspace */
    P4PCATiming *timing /**< [out] optional stage timing */ )
{
    if constexpr( precisionPolicy == detail::P4PCAPrecisionPolicy::doubleDouble )
    {
        static_cast<void>( experimentalWorkspace );
        P4PCA::calculateHeldOutProbe( output,
                                      probeResiduals,
                                      predictors,
                                      target,
                                      probePredictors,
                                      probeTarget,
                                      exclusions,
                                      modes,
                                      rankTolerance,
                                      doubleWorkspace,
                                      timing );
    }
    else
    {
        static_cast<void>( doubleWorkspace );
        detail::p4PCACalculateHeldOutProbeExperimental( output,
                                                        probeResiduals,
                                                        predictors,
                                                        target,
                                                        probePredictors,
                                                        probeTarget,
                                                        exclusions,
                                                        modes,
                                                        rankTolerance,
                                                        precisionPolicy,
                                                        experimentalWorkspace,
                                                        timing );
    }
}

/// Materialize immutable PCA call targets for one previously validated policy.
P4ReductionPCADispatch
p4ReductionPCADispatch( detail::P4PCAPrecisionPolicy precisionPolicy /**< [in] validated scalar policy */ )
{
    switch( precisionPolicy )
    {
    case detail::P4PCAPrecisionPolicy::doubleDouble:
        return { &p4ReductionCalculateDirect<detail::P4PCAPrecisionPolicy::doubleDouble>,
                 &p4ReductionCalculateCenteredInPlace<detail::P4PCAPrecisionPolicy::doubleDouble>,
                 &p4ReductionCalculateHeldOut<detail::P4PCAPrecisionPolicy::doubleDouble>,
                 &p4ReductionCalculateHeldOutProbe<detail::P4PCAPrecisionPolicy::doubleDouble> };
    case detail::P4PCAPrecisionPolicy::floatDouble:
        return { &p4ReductionCalculateDirect<detail::P4PCAPrecisionPolicy::floatDouble>,
                 &p4ReductionCalculateCenteredInPlace<detail::P4PCAPrecisionPolicy::floatDouble>,
                 &p4ReductionCalculateHeldOut<detail::P4PCAPrecisionPolicy::floatDouble>,
                 &p4ReductionCalculateHeldOutProbe<detail::P4PCAPrecisionPolicy::floatDouble> };
    case detail::P4PCAPrecisionPolicy::floatFloat:
        return { &p4ReductionCalculateDirect<detail::P4PCAPrecisionPolicy::floatFloat>,
                 &p4ReductionCalculateCenteredInPlace<detail::P4PCAPrecisionPolicy::floatFloat>,
                 &p4ReductionCalculateHeldOut<detail::P4PCAPrecisionPolicy::floatFloat>,
                 &p4ReductionCalculateHeldOutProbe<detail::P4PCAPrecisionPolicy::floatFloat> };
    }
    throw std::invalid_argument( "unknown experimental P4 precision policy" );
}

/// Thread-local precision selection active only while the internal entry point calls `reduce()`.
struct P4ReductionPrecisionSelection
{
    bool active{ false }; ///< Whether the calling thread requested an experimental dispatch.

    detail::P4PCAPrecisionPolicy policy{ detail::P4PCAPrecisionPolicy::doubleDouble }; ///< Selected policy.
};

/// Calling-thread selection inherited explicitly by each reduction worker binding.
thread_local P4ReductionPrecisionSelection p4ReductionPrecisionSelection;

/// Restore a calling thread's previous precision selection at scope exit.
class P4ReductionPrecisionSelectionScope
{
  public:
    /// Install one validated reduction-level precision selection.
    explicit P4ReductionPrecisionSelectionScope(
        detail::P4PCAPrecisionPolicy policy /**< [in] validated scalar policy to install */ )
        : m_previous( p4ReductionPrecisionSelection )
    {
        p4ReductionPrecisionSelection = { true, policy };
    }

    /// Restore the selection that preceded this scope.
    ~P4ReductionPrecisionSelectionScope()
    {
        p4ReductionPrecisionSelection = m_previous;
    }

    /// Prevent copying a scope that uniquely owns restoration of the calling-thread selection.
    P4ReductionPrecisionSelectionScope(
        const P4ReductionPrecisionSelectionScope &other /**< [in] disabled copy source */ ) = delete;

    /// Prevent overwriting a scope that uniquely owns restoration of the calling-thread selection.
    P4ReductionPrecisionSelectionScope &
    operator=( const P4ReductionPrecisionSelectionScope &other /**< [in] disabled assignment source */ ) = delete;

  private:
    P4ReductionPrecisionSelection m_previous; ///< Calling-thread state restored at destruction.
};

/// Resolve the calling thread's selection once, before a worker loop begins.
std::optional<P4ReductionPCADispatch> p4ReductionSelectedPCADispatch()
{
    if( !p4ReductionPrecisionSelection.active )
    {
        return std::nullopt;
    }
    return p4ReductionPCADispatch( p4ReductionPrecisionSelection.policy );
}

/// Worker-local binding of one immutable dispatch to one caller-owned typed workspace.
struct P4ReductionWorkerPrecisionContext
{
    const P4ReductionPCADispatch *dispatch{ nullptr };        ///< Shared immutable call targets.

    detail::P4PCAExperimentalWorkspace *workspace{ nullptr }; ///< Private typed policy storage.
};

/// Active worker-local experimental dispatch binding, or null for the ordinary production path.
thread_local const P4ReductionWorkerPrecisionContext *p4ReductionWorkerPrecisionContext{ nullptr };

/// Restore a worker thread's previous dispatch binding at scope exit.
class P4ReductionWorkerPrecisionScope
{
  public:
    /// Bind one worker-private experimental workspace to the preselected dispatch.
    P4ReductionWorkerPrecisionScope(
        const P4ReductionPCADispatch &dispatch, /**< [in] immutable selected call targets */
        detail::P4PCAExperimentalWorkspace &workspace /**< [in,out] worker-private typed storage */ )
        : m_previous( p4ReductionWorkerPrecisionContext ), m_context{ &dispatch, &workspace }
    {
        p4ReductionWorkerPrecisionContext = &m_context;
    }

    /// Restore the worker binding that preceded this scope.
    ~P4ReductionWorkerPrecisionScope()
    {
        p4ReductionWorkerPrecisionContext = m_previous;
    }

    /// Prevent copying a scope that uniquely owns restoration of the worker binding.
    P4ReductionWorkerPrecisionScope( const P4ReductionWorkerPrecisionScope &other /**< [in] disabled copy source */ ) =
        delete;

    /// Prevent overwriting a scope that uniquely owns restoration of the worker binding.
    P4ReductionWorkerPrecisionScope &
    operator=( const P4ReductionWorkerPrecisionScope &other /**< [in] disabled assignment source */ ) = delete;

  private:
    const P4ReductionWorkerPrecisionContext *m_previous; ///< Binding restored at destruction.

    P4ReductionWorkerPrecisionContext m_context;         ///< Binding active for this scope.
};

/// Invoke the active direct evaluator, falling back to the production path outside an experiment.
void p4ReductionDispatchDirect( P4PCAResult &output,              /**< [out] completed regression result */
                                const P4PCA::matrixT &predictors, /**< [in] sampled predictor matrix */
                                const P4PCA::vectorT &target,     /**< [in] sampled target series */
                                const std::vector<int> &modes,    /**< [in] requested retained-mode counts */
                                double rankTolerance,             /**< [in] relative numerical-rank threshold */
                                P4PCA::workspaceT &workspace,     /**< [in,out] ordinary FP64 workspace */
                                detail::P4PCAMixedWorkspace &mixedWorkspace,
                                /**< [in,out] production mixed-precision workspace */
                                P4PCATiming *timing,              /**< [out] optional stage timing */
                                P4PCA::matrixT *coefficients /**< [out] optional predictor coefficients */ )
{
    if( !p4ReductionWorkerPrecisionContext )
    {
        static_cast<void>( workspace );
        detail::p4PCACalculateMixed(
            output, predictors, target, modes, rankTolerance, mixedWorkspace, timing, coefficients );
        return;
    }
    p4ReductionWorkerPrecisionContext->dispatch->direct( output,
                                                         predictors,
                                                         target,
                                                         modes,
                                                         rankTolerance,
                                                         workspace,
                                                         *p4ReductionWorkerPrecisionContext->workspace,
                                                         timing,
                                                         coefficients );
}

/// Invoke the active centered evaluator, falling back to the production path outside an experiment.
void p4ReductionDispatchCenteredInPlace( P4PCAResult &output,        /**< [out] completed centered regression result */
                                         P4PCA::matrixT &predictors, /**< [in,out] predictors, centered in place */
                                         const P4PCA::vectorT &target,  /**< [in] sampled target series */
                                         const std::vector<int> &modes, /**< [in] requested retained-mode counts */
                                         double rankTolerance,          /**< [in] relative numerical-rank threshold */
                                         P4PCA::workspaceT &workspace,  /**< [in,out] ordinary FP64 workspace */
                                         detail::P4PCAMixedWorkspace &mixedWorkspace,
                                         /**< [in,out] production mixed-precision workspace */
                                         P4PCATiming *timing,           /**< [out] optional stage timing */
                                         P4PCA::matrixT *coefficients /**< [out] optional predictor coefficients */ )
{
    if( !p4ReductionWorkerPrecisionContext )
    {
        static_cast<void>( workspace );
        detail::p4PCACalculateCenteredInPlaceMixed(
            output, predictors, target, modes, rankTolerance, mixedWorkspace, timing, coefficients );
        return;
    }
    p4ReductionWorkerPrecisionContext->dispatch->centeredInPlace( output,
                                                                  predictors,
                                                                  target,
                                                                  modes,
                                                                  rankTolerance,
                                                                  workspace,
                                                                  *p4ReductionWorkerPrecisionContext->workspace,
                                                                  timing,
                                                                  coefficients );
}

/// Invoke the active held-out evaluator, falling back to the production path outside an experiment.
void p4ReductionDispatchHeldOut( P4PCAResult &output,              /**< [out] completed held-out regression result */
                                 const P4PCA::matrixT &predictors, /**< [in] sampled predictor matrix */
                                 const P4PCA::vectorT &target,     /**< [in] sampled target series */
                                 const P4TargetExclusions &exclusions, /**< [in] target-specific deleted rows */
                                 const std::vector<int> &modes,        /**< [in] requested retained-mode counts */
                                 double rankTolerance,                 /**< [in] relative numerical-rank threshold */
                                 P4PCA::workspaceT &workspace,         /**< [in,out] ordinary FP64 workspace */
                                 detail::P4PCAMixedWorkspace &mixedWorkspace,
                                 /**< [in,out] production mixed-precision workspace */
                                 P4PCATiming *timing /**< [out] optional stage timing */ )
{
    if( !p4ReductionWorkerPrecisionContext )
    {
        static_cast<void>( workspace );
        detail::p4PCACalculateHeldOutMixed(
            output, predictors, target, exclusions, modes, rankTolerance, mixedWorkspace, timing );
        return;
    }
    p4ReductionWorkerPrecisionContext->dispatch->heldOut( output,
                                                          predictors,
                                                          target,
                                                          exclusions,
                                                          modes,
                                                          rankTolerance,
                                                          workspace,
                                                          *p4ReductionWorkerPrecisionContext->workspace,
                                                          timing );
}

/// Invoke the active held-out probe evaluator, falling back to the production path outside an experiment.
void p4ReductionDispatchHeldOutProbe(
    P4PCAResult &output,                   /**< [out] completed held-out regression result */
    P4PCA::matrixT &probeResiduals,        /**< [out] frozen-probe residual responses */
    const P4PCA::matrixT &predictors,      /**< [in] sampled predictor matrix */
    const P4PCA::vectorT &target,          /**< [in] sampled target series */
    const P4PCA::matrixT &probePredictors, /**< [in] frozen-probe predictor responses */
    const P4PCA::vectorT &probeTarget,     /**< [in] frozen-probe direct response */
    const P4TargetExclusions &exclusions,  /**< [in] target-specific deleted rows */
    const std::vector<int> &modes,         /**< [in] requested retained-mode counts */
    double rankTolerance,                  /**< [in] relative numerical-rank threshold */
    P4PCA::workspaceT &workspace,          /**< [in,out] ordinary FP64 workspace */
    detail::P4PCAMixedWorkspace &mixedWorkspace, /**< [in,out] production mixed-precision workspace */
    P4PCATiming *timing /**< [out] optional stage timing */ )
{
    if( !p4ReductionWorkerPrecisionContext )
    {
        static_cast<void>( workspace );
        detail::p4PCACalculateHeldOutProbeMixed( output,
                                                 probeResiduals,
                                                 predictors,
                                                 target,
                                                 probePredictors,
                                                 probeTarget,
                                                 exclusions,
                                                 modes,
                                                 rankTolerance,
                                                 mixedWorkspace,
                                                 timing );
        return;
    }
    p4ReductionWorkerPrecisionContext->dispatch->heldOutProbe( output,
                                                               probeResiduals,
                                                               predictors,
                                                               target,
                                                               probePredictors,
                                                               probeTarget,
                                                               exclusions,
                                                               modes,
                                                               rankTolerance,
                                                               workspace,
                                                               *p4ReductionWorkerPrecisionContext->workspace,
                                                               timing );
}

/** \endcond */
#endif // HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION

/// Read Linux MemAvailable in bytes, or report that automatic discovery is unavailable.
std::optional<std::size_t> p4AvailableMemoryBytes()
{
    std::ifstream memoryInformation( "/proc/meminfo" );
    std::string key;
    std::uint64_t valueKiB{ 0 };
    std::string unit;
    while( memoryInformation >> key >> valueKiB >> unit )
    {
        if( key != "MemAvailable:" )
        {
            continue;
        }
        if( unit != "kB" || valueKiB > std::numeric_limits<std::size_t>::max() / 1024 )
        {
            return std::nullopt;
        }
        return static_cast<std::size_t>( valueKiB ) * 1024;
    }
    return std::nullopt;
}

} // namespace

template <typename realT, class derotFunctObj, class verboseT>
P4Reduction<realT, derotFunctObj, verboseT>::P4Reduction()
{
    static_assert( std::is_same_v<realT, float>,
                   "the initial P4Reduction specialization requires float image storage and cubic kernels" );
    this->m_rejectNonFiniteTargetInput = true;
}

template <typename realT, class derotFunctObj, class verboseT>
void P4Reduction<realT, derotFunctObj, verboseT>::setupConfig( mx::app::appConfigurator &config )
{
    HCIobservation<realT, verboseT>::setupConfig( config );
    ADIobservation<realT, derotFunctObj, verboseT>::setupConfig( config );

    config.add( "geom.minRadius",
                "",
                "geom.minRadius",
                mx::app::argType::Required,
                "geom",
                "minRadius",
                false,
                "vector<realT>",
                "Inclusive inner radii of the ordered, non-overlapping P4 search annuli" );
    config.add( "geom.maxRadius",
                "",
                "geom.maxRadius",
                mx::app::argType::Required,
                "geom",
                "maxRadius",
                false,
                "vector<realT>",
                "Exclusive outer radii of the ordered, non-overlapping P4 search annuli" );
    config.add( "p4.modeFractions",
                "",
                "p4.modeFractions",
                mx::app::argType::Required,
                "p4",
                "modeFractions",
                false,
                "vector<realT>",
                "Strictly increasing PCA fractions in (0,1]" );
    config.add( "p4.regressionFrame",
                "",
                "p4.regressionFrame",
                mx::app::argType::Required,
                "p4",
                "regressionFrame",
                false,
                "string",
                "Regression coordinate frame: detector or rotated; default detector" );
    config.add( "p4.numberImages",
                "",
                "p4.numberImages",
                mx::app::argType::Required,
                "p4",
                "numberImages",
                false,
                "int",
                "Qualifying earlier and later detector images appended to each predictor row; default 0" );
    config.add( "adi.minDPx",
                "",
                "adi.minDPx",
                mx::app::argType::Required,
                "adi",
                "minDPx",
                false,
                "float",
                "Minimum target/reference displacement interpreted by adi.excludeMethod; default 0" );
    config.add( "adi.excludeMethod",
                "",
                "adi.excludeMethod",
                mx::app::argType::Required,
                "adi",
                "excludeMethod",
                false,
                "string",
                "Target-frame exclusion method: none, pixel, angle, or imno; default none" );
    config.add( "solver.exclusionSolver",
                "",
                "solver.exclusionSolver",
                mx::app::argType::Required,
                "solver",
                "exclusionSolver",
                false,
                "string",
                "Target-exclusion solver: explicitRefit or factorDowndateExact; default explicitRefit" );
    config.add( "solver.deletionBackend",
                "",
                "solver.deletionBackend",
                mx::app::argType::Required,
                "solver",
                "deletionBackend",
                false,
                "string",
                "Factor row-deletion backend: leadingCovariance or rankOneSecular; default leadingCovariance" );
    config.add( "p4.orDeltaRadiusInner",
                "",
                "p4.orDeltaRadiusInner",
                mx::app::argType::Required,
                "p4",
                "orDeltaRadiusInner",
                false,
                "float",
                "Positive inward radial extent of the local predictor wedge" );
    config.add( "p4.orDeltaRadiusOuter",
                "",
                "p4.orDeltaRadiusOuter",
                mx::app::argType::Required,
                "p4",
                "orDeltaRadiusOuter",
                false,
                "float",
                "Positive outward radial extent of the local predictor wedge" );
    config.add(
        "p4.orArcHalfWidth",
        "",
        "p4.orArcHalfWidth",
        mx::app::argType::Required,
        "p4",
        "orArcHalfWidth",
        false,
        "float",
        "Nonnegative predictor-wedge azimuthal half-width as arc length in pixels; zero uses only orMaxHalfAngle" );
    config.add( "p4.orMaxHalfAngle",
                "",
                "p4.orMaxHalfAngle",
                mx::app::argType::Required,
                "p4",
                "orMaxHalfAngle",
                false,
                "float",
                "Predictor-wedge angular half-width cap in the range (0,180] degrees" );
    config.add( "p4.psfRadius",
                "",
                "p4.psfRadius",
                mx::app::argType::Required,
                "p4",
                "psfRadius",
                false,
                "float",
                "Positive physical central signal-exclusion radius in pixels" );
    config.add( "p4.psfFile",
                "",
                "p4.psfFile",
                mx::app::argType::Required,
                "p4",
                "psfFile",
                false,
                "string",
                "Optional post-preprocessing centered FITS PSF template enabling frozen-model calculation" );
    config.add( "p4.psfStampSize",
                "",
                "p4.psfStampSize",
                mx::app::argType::Required,
                "p4",
                "psfStampSize",
                false,
                "int",
                "Positive square frozen-model PSF stamp size; required when p4.psfFile is set" );
    config.add( "p4.outputPSFModels",
                "",
                "p4.outputPSFModels",
                mx::app::argType::Optional,
                "p4",
                "outputPSFModels",
                false,
                "bool",
                "Reconstruct and write compact final-frame frozen-model PSF products" );
    config.add( "p4.psfFilter",
                "",
                "p4.psfFilter",
                mx::app::argType::Optional,
                "p4",
                "psfFilter",
                false,
                "bool",
                "Apply the spatially variable normalized PSF filter and write separate full-image products" );
    config.add( "p4.psfFilterMinGoodFract",
                "",
                "p4.psfFilterMinGoodFract",
                mx::app::argType::Required,
                "p4",
                "psfFilterMinGoodFract",
                false,
                "float",
                "Minimum usable local-stamp fraction for PSF filtering in [0,1]; default 1" );
    config.add( "p4.psfOutputPrefix",
                "",
                "p4.psfOutputPrefix",
                mx::app::argType::Required,
                "p4",
                "psfOutputPrefix",
                false,
                "string",
                "Prefix for compact PSF products inside the final image's _outputs directory; default p4PSF_" );
    config.add( "p4.localStampSize",
                "",
                "p4.localStampSize",
                mx::app::argType::Required,
                "p4",
                "localStampSize",
                false,
                "int",
                "Positive odd pixel-local result width; zero disables local finite-amplitude refitting" );
    config.add( "p4.exclusionPolicy",
                "",
                "p4.exclusionPolicy",
                mx::app::argType::Required,
                "p4",
                "exclusionPolicy",
                false,
                "string",
                "Required central exclusion rule: sampleCenter or kernelSupport" );
    config.add( "p4.exclusionRadiusBuffer",
                "",
                "p4.exclusionRadiusBuffer",
                mx::app::argType::Required,
                "p4",
                "exclusionRadiusBuffer",
                false,
                "float",
                "Explicit nonnegative radius buffer added to p4.psfRadius" );
    config.add( "p4.rankTolerance",
                "",
                "p4.rankTolerance",
                mx::app::argType::Required,
                "p4",
                "rankTolerance",
                false,
                "double",
                "Relative numerical-rank threshold in [0,1), default 1e-12" );
    config.add( "p4.memoryFraction",
                "",
                "p4.memoryFraction",
                mx::app::argType::Required,
                "p4",
                "memoryFraction",
                false,
                "double",
                "Fraction of Linux MemAvailable usable by future P4 allocations in [0,1]; zero disables limiting" );
    config.add( "p4.writeDiagnostics",
                "",
                "p4.writeDiagnostics",
                mx::app::argType::Optional,
                "p4",
                "writeDiagnostics",
                false,
                "bool",
                "Whether P4 ownership, validity, geometry, rank, and timing diagnostics are written" );
    config.add( "p4.diagnosticDirectory",
                "",
                "p4.diagnosticDirectory",
                mx::app::argType::Required,
                "p4",
                "diagnosticDirectory",
                false,
                "string",
                "Explicit P4 diagnostic destination; default . groups diagnostics in the final image's _outputs "
                "directory" );
}

template <typename realT, class derotFunctObj, class verboseT>
void P4Reduction<realT, derotFunctObj, verboseT>::loadConfig( mx::app::appConfigurator &config )
{
    HCIobservation<realT, verboseT>::loadConfig( config );
    ADIobservation<realT, derotFunctObj, verboseT>::loadConfig( config );

    config( m_minRadius, "geom.minRadius" );
    config( m_maxRadius, "geom.maxRadius" );
    config( m_modeFractions, "p4.modeFractions" );

    std::string regressionFrame = regressionFrameString( m_regressionFrame );
    config( regressionFrame, "p4.regressionFrame" );
    try
    {
        m_regressionFrame = parseRegressionFrame( regressionFrame );
    }
    catch( ... )
    {
        std::throw_with_nested(
            mx::exception<verboseT>( mx::error_t::invalidconfig, "p4.regressionFrame is not valid" ) );
    }

    std::string deletionBackend = deletionBackendString( m_deletionBackend );
    config( deletionBackend, "solver.deletionBackend" );
    try
    {
        m_deletionBackend = parseDeletionBackend( deletionBackend );
    }
    catch( ... )
    {
        std::throw_with_nested(
            mx::exception<verboseT>( mx::error_t::invalidconfig, "solver.deletionBackend is not valid" ) );
    }

    config( m_numberImages, "p4.numberImages" );
    config( m_minDPx, "adi.minDPx" );
    std::string excludeMethod = HCI::excludeToStr<verboseT>( m_excludeMethod );
    config( excludeMethod, "adi.excludeMethod" );
    m_excludeMethod = HCI::excludeFmStr<verboseT>( excludeMethod );

    std::string exclusionSolver = exclusionSolverString( m_exclusionSolver );
    config( exclusionSolver, "solver.exclusionSolver" );
    try
    {
        m_exclusionSolver = parseExclusionSolver( exclusionSolver );
    }
    catch( ... )
    {
        std::throw_with_nested(
            mx::exception<verboseT>( mx::error_t::invalidconfig, "solver.exclusionSolver is not valid" ) );
    }

    config( m_orDeltaRadiusInner, "p4.orDeltaRadiusInner" );
    config( m_orDeltaRadiusOuter, "p4.orDeltaRadiusOuter" );
    config( m_orArcHalfWidth, "p4.orArcHalfWidth" );
    config( m_orMaxHalfAngle, "p4.orMaxHalfAngle" );
    config( m_psfRadius, "p4.psfRadius" );
    config( m_psfFile, "p4.psfFile" );
    config( m_psfStampSize, "p4.psfStampSize" );
    loadBoolConfig<verboseT>( config, m_outputPSFModels, "p4.outputPSFModels" );
    loadBoolConfig<verboseT>( config, m_psfFilter, "p4.psfFilter" );
    config( m_psfFilterMinGoodFract, "p4.psfFilterMinGoodFract" );
    config( m_psfOutputPrefix, "p4.psfOutputPrefix" );
    config( m_localStampSize, "p4.localStampSize" );

    std::string policy;
    if( m_exclusionPolicy.has_value() )
    {
        policy = exclusionPolicyString( *m_exclusionPolicy );
    }
    config( policy, "p4.exclusionPolicy" );
    if( policy.empty() )
    {
        m_exclusionPolicy.reset();
    }
    else
    {
        try
        {
            m_exclusionPolicy = parseExclusionPolicy( policy );
        }
        catch( ... )
        {
            std::throw_with_nested(
                mx::exception<verboseT>( mx::error_t::invalidconfig, "p4.exclusionPolicy is not valid" ) );
        }
    }

    config( m_exclusionRadiusBuffer, "p4.exclusionRadiusBuffer" );
    config( m_rankTolerance, "p4.rankTolerance" );
    config( m_memoryFraction, "p4.memoryFraction" );
    loadBoolConfig<verboseT>( config, m_writeDiagnostics, "p4.writeDiagnostics" );
    config( m_diagnosticDirectory, "p4.diagnosticDirectory" );
}

template <typename realT, class derotFunctObj, class verboseT>
bool P4Reduction<realT, derotFunctObj, verboseT>::deferTargetFakeInjection() const
{
    return m_localStampSize > 0;
}

template <typename realT, class derotFunctObj, class verboseT>
std::string P4Reduction<realT, derotFunctObj, verboseT>::exclusionPolicyString( P4ExclusionPolicy policy )
{
    if( policy == P4ExclusionPolicy::sampleCenter )
    {
        return "sampleCenter";
    }
    if( policy == P4ExclusionPolicy::kernelSupport )
    {
        return "kernelSupport";
    }
    throw std::invalid_argument( "unsupported P4 exclusion policy" );
}

template <typename realT, class derotFunctObj, class verboseT>
P4ExclusionPolicy P4Reduction<realT, derotFunctObj, verboseT>::parseExclusionPolicy( const std::string &value )
{
    if( value == "sampleCenter" )
    {
        return P4ExclusionPolicy::sampleCenter;
    }
    if( value == "kernelSupport" )
    {
        return P4ExclusionPolicy::kernelSupport;
    }
    throw std::invalid_argument( "unsupported P4 exclusion policy: " + value );
}

template <typename realT, class derotFunctObj, class verboseT>
std::string P4Reduction<realT, derotFunctObj, verboseT>::exclusionSolverString( P4ExclusionSolver solver )
{
    if( solver == P4ExclusionSolver::explicitRefit )
    {
        return "explicitRefit";
    }
    if( solver == P4ExclusionSolver::factorDowndateExact )
    {
        return "factorDowndateExact";
    }
    throw std::invalid_argument( "unsupported P4 target-exclusion solver" );
}

template <typename realT, class derotFunctObj, class verboseT>
P4ExclusionSolver P4Reduction<realT, derotFunctObj, verboseT>::parseExclusionSolver( const std::string &value )
{
    if( value == "explicitRefit" )
    {
        return P4ExclusionSolver::explicitRefit;
    }
    if( value == "factorDowndateExact" )
    {
        return P4ExclusionSolver::factorDowndateExact;
    }
    throw std::invalid_argument( "unsupported P4 target-exclusion solver: " + value );
}

template <typename realT, class derotFunctObj, class verboseT>
std::string P4Reduction<realT, derotFunctObj, verboseT>::deletionBackendString( mx::math::svdDeletionBackend backend )
{
    if( backend == mx::math::svdDeletionBackend::leadingCovariance ||
        backend == mx::math::svdDeletionBackend::rankOneSecular )
    {
        return mx::math::svdDeletionBackendName( backend );
    }
    throw std::invalid_argument( "unsupported P4 row-deletion backend" );
}

template <typename realT, class derotFunctObj, class verboseT>
mx::math::svdDeletionBackend
P4Reduction<realT, derotFunctObj, verboseT>::parseDeletionBackend( const std::string &value )
{
    if( value == "leadingCovariance" )
    {
        return mx::math::svdDeletionBackend::leadingCovariance;
    }
    if( value == "rankOneSecular" )
    {
        return mx::math::svdDeletionBackend::rankOneSecular;
    }
    throw std::invalid_argument( "unsupported P4 row-deletion backend: " + value );
}

template <typename realT, class derotFunctObj, class verboseT>
std::string P4Reduction<realT, derotFunctObj, verboseT>::regressionFrameString( P4RegressionFrame frame )
{
    if( frame == P4RegressionFrame::detector )
    {
        return "detector";
    }
    if( frame == P4RegressionFrame::rotated )
    {
        return "rotated";
    }
    throw std::invalid_argument( "unsupported P4 regression frame" );
}

template <typename realT, class derotFunctObj, class verboseT>
P4RegressionFrame P4Reduction<realT, derotFunctObj, verboseT>::parseRegressionFrame( const std::string &value )
{
    if( value == "detector" )
    {
        return P4RegressionFrame::detector;
    }
    if( value == "rotated" )
    {
        return P4RegressionFrame::rotated;
    }
    throw std::invalid_argument( "unsupported P4 regression frame: " + value );
}

template <typename realT, class derotFunctObj, class verboseT>
void P4Reduction<realT, derotFunctObj, verboseT>::validateConfiguration() const
{
    static_cast<void>( regressionFrameString( m_regressionFrame ) );
    static_cast<void>( HCI::excludeToStr<verboseT>( m_excludeMethod ) );
    static_cast<void>( exclusionSolverString( m_exclusionSolver ) );
    static_cast<void>( deletionBackendString( m_deletionBackend ) );
    if( !mx::math::isFinite( m_minDPx ) || m_minDPx < 0 )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig, "adi.minDPx must be finite and nonnegative" );
    }
    if( m_excludeMethod != HCI::exclude::none )
    {
        if( m_regressionFrame != P4RegressionFrame::detector )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "adi.excludeMethod is initially supported only for detector-frame P4" );
        }
        if( m_numberImages > 0 )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "P4 target-frame exclusion initially requires p4.numberImages=0" );
        }
    }
    else if( m_exclusionSolver != P4ExclusionSolver::explicitRefit )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                       "solver.exclusionSolver requires adi.excludeMethod other than none" );
    }
    if( m_deletionBackend == mx::math::svdDeletionBackend::rankOneSecular &&
        m_exclusionSolver != P4ExclusionSolver::factorDowndateExact )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                       "solver.deletionBackend=rankOneSecular requires "
                                       "solver.exclusionSolver=factorDowndateExact" );
    }
    if( m_localStampSize < 0 )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig, "p4.localStampSize must be nonnegative" );
    }
    if( m_localStampSize > 0 )
    {
        if( m_localStampSize % 2 == 0 )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "p4.localStampSize must be odd when local processing is enabled" );
        }
        if( m_regressionFrame != P4RegressionFrame::detector )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "p4.localStampSize is supported only for detector-frame P4" );
        }
        if( !this->m_skipPreProcess )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "p4.localStampSize initially requires preProcess.skip=true" );
        }
        if( !this->m_doDerotate )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig, "p4.localStampSize requires ADI derotation" );
        }
        if( this->m_postMedSub )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "p4.localStampSize does not yet support adi.postMedSub=true" );
        }
        if( this->m_combineMethod == HCI::combine::none )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "p4.localStampSize requires a final combination method" );
        }
        if( this->m_doOutputPSFSub )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "p4.localStampSize does not support output.outputPSFSub=true" );
        }
        if( this->m_fakeMethod != HCI::fake::single || this->m_fakeFileName.empty() )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "p4.localStampSize requires fake.method=single and fake.fileName" );
        }
        if( this->m_fakeSep.size() != 1 || this->m_fakePA.size() != 1 || this->m_fakeContrast.size() != 1 )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "p4.localStampSize requires exactly one fake.sep, fake.PA, and "
                                           "fake.contrast value" );
        }
        if( !mx::math::isFinite( this->m_fakeSep[0] ) || this->m_fakeSep[0] < 0 ||
            !mx::math::isFinite( this->m_fakePA[0] ) || !mx::math::isFinite( this->m_fakeContrast[0] ) )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "local fake separation, PA, and contrast must be finite" );
        }
        if( !m_psfFile.empty() || m_outputPSFModels || m_psfFilter )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "p4.localStampSize does not yet support P4 frozen-PSF products" );
        }
    }
    if( m_regressionFrame == P4RegressionFrame::rotated && this->m_postMedSub )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                       "adi.postMedSub is not supported for rotated-frame P4 regression" );
    }
    if( m_numberImages < 0 || m_numberImages > ( std::numeric_limits<int>::max() - 1 ) / 2 )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                       "p4.numberImages must be nonnegative and within the supported range" );
    }
    if( m_regressionFrame == P4RegressionFrame::rotated && m_numberImages > 0 )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                       "p4.numberImages is supported only for detector-frame P4 regression" );
    }
    if( !m_psfFile.empty() )
    {
        if( m_psfStampSize <= 0 )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "p4.psfStampSize must be positive when p4.psfFile is set" );
        }
        if( m_regressionFrame != P4RegressionFrame::detector )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "p4.psfFile is initially supported only for detector-frame P4" );
        }
        if( this->m_postMedSub )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "p4.psfFile is not yet supported with adi.postMedSub=true" );
        }
    }
    if( m_outputPSFModels || m_psfFilter )
    {
        if( m_psfFile.empty() )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "P4 PSF output or filtering requires p4.psfFile" );
        }
        if( m_psfOutputPrefix.empty() )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "p4.psfOutputPrefix must not be empty when PSF products are enabled" );
        }
        if( this->m_combineMethod == HCI::combine::none )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "P4 PSF output or filtering requires a final combination method" );
        }
    }
    if( !mx::math::isFinite( m_psfFilterMinGoodFract ) || m_psfFilterMinGoodFract < 0 || m_psfFilterMinGoodFract > 1 )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                       "p4.psfFilterMinGoodFract must be finite and in [0,1]" );
    }
    if( m_psfFilter )
    {
        if( m_psfStampSize % 2 == 0 )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig, "p4.psfFilter requires an odd p4.psfStampSize" );
        }
    }

    if( m_minRadius.empty() || m_minRadius.size() != m_maxRadius.size() )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                       "P4 search-radius vectors must be nonempty and have equal lengths" );
    }
    if( m_modeFractions.empty() )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig, "P4 requires at least one mode fraction" );
    }
    if( !m_exclusionPolicy.has_value() )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig, "P4 requires an explicit exclusion policy" );
    }
    static_cast<void>( exclusionPolicyString( *m_exclusionPolicy ) );

    const realT geometry[]{ m_orDeltaRadiusInner,
                            m_orDeltaRadiusOuter,
                            m_orArcHalfWidth,
                            m_orMaxHalfAngle,
                            m_psfRadius,
                            m_exclusionRadiusBuffer };
    for( const realT value : geometry )
    {
        if( !mx::math::isFinite( value ) )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig, "P4 geometry values must be finite" );
        }
    }
    if( m_orDeltaRadiusInner <= 0 || m_orDeltaRadiusOuter <= 0 || m_orArcHalfWidth < 0 || m_orMaxHalfAngle <= 0 ||
        m_orMaxHalfAngle > 180 || m_psfRadius <= 0 || m_exclusionRadiusBuffer < 0 )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig, "P4 geometry values are outside valid ranges" );
    }
    if( !mx::math::isFinite( m_rankTolerance ) || m_rankTolerance < 0 || m_rankTolerance >= 1 )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig, "p4.rankTolerance must be finite and in [0,1)" );
    }
    if( !mx::math::isFinite( m_memoryFraction ) || m_memoryFraction < 0 || m_memoryFraction > 1 )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig, "p4.memoryFraction must be finite and in [0,1]" );
    }
    if( m_writeDiagnostics && m_diagnosticDirectory.empty() )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                       "p4.diagnosticDirectory must not be empty when diagnostics are enabled" );
    }

    realT previousOuter{ -1 };
    for( std::size_t region = 0; region < m_minRadius.size(); ++region )
    {
        if( !mx::math::isFinite( m_minRadius[region] ) || !mx::math::isFinite( m_maxRadius[region] ) ||
            m_minRadius[region] < 0 || m_maxRadius[region] <= m_minRadius[region] ||
            ( region != 0 && m_minRadius[region] < previousOuter ) )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "P4 search annuli must be finite, ordered, and non-overlapping" );
        }
        previousOuter = m_maxRadius[region];
    }

    realT previousFraction{ 0 };
    for( const realT fraction : m_modeFractions )
    {
        if( !mx::math::isFinite( fraction ) || fraction <= 0 || fraction > 1 || fraction <= previousFraction )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "P4 mode fractions must be finite, strictly increasing, and in (0,1]" );
        }
        previousFraction = fraction;
    }
}

template <typename realT, class derotFunctObj, class verboseT>
bool P4Reduction<realT, derotFunctObj, verboseT>::targetCubeFinite() const
{
    for( int image = 0; image < this->m_tgtIms.planes(); ++image )
    {
        for( int column = 0; column < this->m_tgtIms.cols(); ++column )
        {
            for( int row = 0; row < this->m_tgtIms.rows(); ++row )
            {
                if( !mx::math::isFinite( this->m_tgtIms.image( image )( row, column ) ) )
                {
                    return false;
                }
            }
        }
    }
    return true;
}

template <typename realT, class derotFunctObj, class verboseT>
double P4Reduction<realT, derotFunctObj, verboseT>::checkedPredictorPromotion( realT value )
{
    if( !mx::math::isFinite( value ) )
    {
        throw std::runtime_error( "P4 interpolation produced a nonfinite predictor value" );
    }
    return static_cast<double>( value );
}

template <typename realT, class derotFunctObj, class verboseT>
realT P4Reduction<realT, derotFunctObj, verboseT>::checkedResidualCast( double value )
{
    if( !mx::math::isFinite( value ) )
    {
        throw std::runtime_error( "P4PCA returned a nonfinite rank-supported residual" );
    }
    const realT converted = static_cast<realT>( value );
    if( !mx::math::isFinite( converted ) )
    {
        throw std::overflow_error( "P4 residual exceeds image-storage range" );
    }
    return converted;
}

template <typename realT, class derotFunctObj, class verboseT>
int P4Reduction<realT, derotFunctObj, verboseT>::checkedMaximumDegreesOfFreedom( int imageCount,
                                                                                 std::size_t predictorCount,
                                                                                 bool temporallyCentered )
{
    if( predictorCount > static_cast<std::size_t>( std::numeric_limits<int>::max() ) )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "P4 predictor count exceeds supported range" );
    }
    if( temporallyCentered && imageCount < 2 )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                       "rotated-frame P4 requires at least two temporal samples" );
    }
    const int temporalDegreesOfFreedom = imageCount - ( temporallyCentered ? 1 : 0 );
    return std::min( temporalDegreesOfFreedom, static_cast<int>( predictorCount ) );
}

template <typename realT, class derotFunctObj, class verboseT>
P4TargetExclusions
P4Reduction<realT, derotFunctObj, verboseT>::targetExclusions( const std::vector<std::vector<int>> &selections,
                                                               HCI::exclude method,
                                                               realT minDPx,
                                                               realT minimumRadius,
                                                               const std::vector<double> &derotationAngles )
{
    if( method == HCI::exclude::none )
    {
        return P4TargetExclusions();
    }
    if( selections.empty() || !mx::math::isFinite( minDPx ) || minDPx < 0 )
    {
        throw std::invalid_argument( "P4 target exclusion requires temporal rows and a nonnegative finite threshold" );
    }

    double threshold = static_cast<double>( minDPx );
    if( method == HCI::exclude::pixel )
    {
        if( !mx::math::isFinite( minimumRadius ) || minimumRadius <= 0 )
        {
            throw std::invalid_argument( "pixel target exclusion requires a positive annulus inner radius" );
        }
        threshold = std::abs( std::atan( threshold / static_cast<double>( minimumRadius ) ) );
    }
    else if( method == HCI::exclude::angle )
    {
        threshold = mx::math::dtor( threshold );
    }
    else if( method != HCI::exclude::imno )
    {
        throw std::invalid_argument( "unsupported P4 target-frame exclusion method" );
    }

    for( const auto &selection : selections )
    {
        if( selection.empty() )
        {
            throw std::invalid_argument( "P4 target exclusion encountered an empty temporal selection" );
        }
    }

    const Eigen::Index targetCount = static_cast<Eigen::Index>( selections.size() );
    if( method == HCI::exclude::imno )
    {
        std::vector<P4ExclusionSpan> spans( selections.size() );
        for( Eigen::Index target = 0; target < targetCount; ++target )
        {
            const double first = std::ceil( std::max( 0.0, static_cast<double>( target ) - threshold ) );
            const double last = std::floor(
                std::min( static_cast<double>( targetCount - 1 ), static_cast<double>( target ) + threshold ) );
            spans[static_cast<std::size_t>( target )] = { static_cast<Eigen::Index>( first ),
                                                          static_cast<Eigen::Index>( last ) + 1 };
        }
        return P4TargetExclusions::fromSpans( targetCount, spans );
    }

    std::vector<std::vector<Eigen::Index>> rows( selections.size() );
    for( std::size_t targetRow = 0; targetRow < selections.size(); ++targetRow )
    {
        const int targetImage = selections[targetRow][0];
        if( targetImage < 0 || static_cast<std::size_t>( targetImage ) >= derotationAngles.size() ||
            !mx::math::isFinite( derotationAngles[static_cast<std::size_t>( targetImage )] ) )
        {
            throw std::invalid_argument( "P4 target exclusion requires a finite angle for every central image" );
        }
        for( std::size_t candidateRow = 0; candidateRow < selections.size(); ++candidateRow )
        {
            const int candidateImage = selections[candidateRow][0];
            if( candidateImage < 0 || static_cast<std::size_t>( candidateImage ) >= derotationAngles.size() ||
                !mx::math::isFinite( derotationAngles[static_cast<std::size_t>( candidateImage )] ) )
            {
                throw std::invalid_argument( "P4 target exclusion requires a finite angle for every central image" );
            }
            const double difference = std::abs( mx::math::angleDiff<mx::math::radiansT<double>>(
                derotationAngles[static_cast<std::size_t>( candidateImage )],
                derotationAngles[static_cast<std::size_t>( targetImage )] ) );
            if( difference <= threshold )
            {
                rows[targetRow].push_back( static_cast<Eigen::Index>( candidateRow ) );
            }
        }
    }
    return P4TargetExclusions::fromExplicit( targetCount, rows );
}

template <typename realT, class derotFunctObj, class verboseT>
std::size_t P4Reduction<realT, derotFunctObj, verboseT>::estimatedWorkerBytes( std::size_t targetImageCount,
                                                                               std::size_t predictorCount,
                                                                               std::size_t modeCount,
                                                                               bool includeCoefficients,
                                                                               std::size_t psfStampPixels,
                                                                               bool exactHeldOut,
                                                                               P4ExclusionSolver exclusionSolver,
                                                                               std::size_t maximumDeletedRows,
                                                                               std::size_t maximumRetainedMode,
                                                                               std::size_t probeSampleCount )
{
    if( targetImageCount == 0 || predictorCount == 0 || modeCount == 0 )
    {
        throw std::invalid_argument( "P4 worker dimensions must be positive" );
    }
    const long double targets = static_cast<long double>( targetImageCount );
    const long double predictors = static_cast<long double>( predictorCount );
    const long double modes = static_cast<long double>( modeCount );
    const long double dimension = std::min( targets, predictors );
    const long double probeSamples = static_cast<long double>( probeSampleCount );
    const long double maximumDeleted = static_cast<long double>( maximumDeletedRows );
    long double heldOutValues{ 0 };
    if( exactHeldOut )
    {
        if( exclusionSolver == P4ExclusionSolver::factorDowndateExact )
        {
            const long double retainedModes =
                maximumRetainedMode == 0 ? dimension
                                         : std::min( dimension, static_cast<long double>( maximumRetainedMode ) );
            heldOutValues = targets * dimension + 2 * dimension * dimension + 2 * maximumDeleted * dimension +
                            dimension * retainedModes + 10 * dimension;
        }
        else
        {
            heldOutValues = targets <= predictors ? dimension * dimension + targets * modes
                                                  : targets * predictors + predictors * modes + targets * modes;
        }
    }
    const long double probeFactorValues =
        exclusionSolver == P4ExclusionSolver::explicitRefit && targets <= predictors ? 2 * targets : dimension;
    const long double probeValues =
        probeSamples == 0 ? 0 : probeSamples * ( predictors + targets * modes + probeFactorValues + 3 );
    const long double doubleValues =
        targets * predictors + targets * modes + 2 * dimension * dimension + 5 * targets + 2 * predictors +
        70 * dimension + ( includeCoefficients ? predictors * ( modes + 2 ) : 0 ) + heldOutValues + probeValues;
    const bool mixedCalculation = !exactHeldOut || exclusionSolver != P4ExclusionSolver::factorDowndateExact;
    const long double floatValues =
        mixedCalculation
            ? targets * predictors + targets +
                  ( probeSamples == 0 ? 0 : probeSamples * predictors + probeSamples )
            : 0;
    constexpr long double safetyFactor = 1.25L;
    constexpr std::size_t fixedMargin = 1024 * 1024;
    const long double estimated = safetyFactor * ( sizeof( double ) * doubleValues +
                                                   sizeof( float ) * floatValues +
                                                   sizeof( realT ) * static_cast<long double>( psfStampPixels ) ) +
                                  fixedMargin;
    if( estimated > static_cast<long double>( std::numeric_limits<std::size_t>::max() ) )
    {
        throw std::overflow_error( "P4 worker memory estimate exceeds size_t range" );
    }
    return static_cast<std::size_t>( std::ceil( estimated ) );
}

template <typename realT, class derotFunctObj, class verboseT>
int P4Reduction<realT, derotFunctObj, verboseT>::localPSFModelDimension( int outputStampSize, int templateDimension )
{
    if( outputStampSize <= 0 || templateDimension <= 0 )
    {
        throw std::invalid_argument( "P4 PSF output and template dimensions must be positive" );
    }
    const long double outputHalfExtent = 0.5L * static_cast<long double>( outputStampSize - 1 );
    const long double rotatedHalfExtent = std::sqrt( 2.0L ) * outputHalfExtent;

    // Reconstruction composes one four-sample PSF-shift footprint with one four-sample derotation footprint. The
    // phase-matched exact bound differs by half a sample for integer- and half-integer-centered local grids.
    const long double dimension = templateDimension % 2 == 0 ? 2.0L * std::floor( rotatedHalfExtent + 0.5L ) + 8.0L
                                                             : 2.0L * std::floor( rotatedHalfExtent ) + 9.0L;
    if( dimension > static_cast<long double>( std::numeric_limits<int>::max() ) )
    {
        throw std::overflow_error( "P4 local PSF dimension exceeds int range" );
    }
    return static_cast<int>( dimension );
}

template <typename realT, class derotFunctObj, class verboseT>
std::size_t P4Reduction<realT, derotFunctObj, verboseT>::localPSFBytes( std::size_t searchPixelCount,
                                                                        std::size_t modeCount,
                                                                        std::size_t temporalPredictorCount,
                                                                        int stampRows,
                                                                        int stampColumns )
{
    if( searchPixelCount == 0 || modeCount == 0 || stampRows <= 0 || stampColumns <= 0 )
    {
        throw std::invalid_argument( "P4 local PSF dimensions must be positive" );
    }
    const long double elementCount = static_cast<long double>( searchPixelCount ) * modeCount;
    const long double stampPixels = static_cast<long double>( stampRows ) * stampColumns;
    const long double bytes =
        elementCount * ( sizeof( realT ) * ( stampPixels + temporalPredictorCount ) + sizeof( std::uint8_t ) );
    if( bytes > static_cast<long double>( std::numeric_limits<std::size_t>::max() ) )
    {
        throw std::overflow_error( "P4 local PSF byte count exceeds size_t range" );
    }
    return static_cast<std::size_t>( bytes );
}

template <typename realT, class derotFunctObj, class verboseT>
std::size_t P4Reduction<realT, derotFunctObj, verboseT>::psfReconstructionBytes( std::size_t searchPixelCount,
                                                                                 std::size_t targetImageCount,
                                                                                 std::size_t temporalPredictorCount,
                                                                                 int outputStampSize,
                                                                                 int localStampRows,
                                                                                 int localStampColumns )
{
    if( searchPixelCount == 0 || targetImageCount == 0 || outputStampSize <= 0 || localStampRows <= 0 ||
        localStampColumns <= 0 )
    {
        throw std::invalid_argument( "P4 PSF reconstruction dimensions must be positive" );
    }
    const long double searches = static_cast<long double>( searchPixelCount );
    const long double frames = static_cast<long double>( targetImageCount );
    const long double outputPixels = static_cast<long double>( outputStampSize ) * outputStampSize;
    const long double localPixels = static_cast<long double>( localStampRows ) * localStampColumns;
    const long double bytes = searches * ( localPixels + temporalPredictorCount ) * sizeof( realT ) +
                              searches * sizeof( std::uint8_t ) + searches * outputPixels * sizeof( realT ) +
                              searches * sizeof( realT ) + searches * 4 * sizeof( realT ) +
                              2 * outputPixels * frames * sizeof( realT ) + 2 * frames * sizeof( realT ) +
                              outputPixels * ( sizeof( realT ) + sizeof( std::uint8_t ) ) +
                              localPixels * ( sizeof( realT ) + sizeof( std::uint8_t ) );
    if( bytes > static_cast<long double>( std::numeric_limits<std::size_t>::max() ) )
    {
        throw std::overflow_error( "P4 PSF reconstruction byte count exceeds size_t range" );
    }
    return static_cast<std::size_t>( bytes );
}

template <typename realT, class derotFunctObj, class verboseT>
std::size_t P4Reduction<realT, derotFunctObj, verboseT>::psfFilterBytes( int rows, int columns, std::size_t modeCount )
{
    if( rows <= 0 || columns <= 0 || modeCount == 0 )
    {
        throw std::invalid_argument( "P4 PSF filter dimensions must be positive" );
    }
    const long double bytes = 4.0L * static_cast<long double>( rows ) * static_cast<long double>( columns ) *
                              static_cast<long double>( modeCount ) * sizeof( realT );
    if( bytes > static_cast<long double>( std::numeric_limits<std::size_t>::max() ) )
    {
        throw std::overflow_error( "P4 PSF filter byte count exceeds size_t range" );
    }
    return static_cast<std::size_t>( bytes );
}

template <typename realT, class derotFunctObj, class verboseT>
int P4Reduction<realT, derotFunctObj, verboseT>::memoryLimitedWorkerCount( int requestedWorkers,
                                                                           std::size_t budgetBytes,
                                                                           std::size_t persistentBytes,
                                                                           std::size_t workerBytes )
{
    if( requestedWorkers <= 0 || workerBytes == 0 )
    {
        throw std::invalid_argument( "P4 worker-count inputs must be positive" );
    }
    if( budgetBytes <= persistentBytes )
    {
        return 0;
    }
    const std::size_t supportedWorkers = ( budgetBytes - persistentBytes ) / workerBytes;
    return static_cast<int>( std::min<std::size_t>( static_cast<std::size_t>( requestedWorkers ), supportedWorkers ) );
}

template <typename realT, class derotFunctObj, class verboseT>
void P4Reduction<realT, derotFunctObj, verboseT>::claimOwnership(
    Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic> &ownership, const P4PixelCoordinate &coordinate, int region )
{
    if( ownership( coordinate.row(), coordinate.column() ) != -1 )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig, "P4 search annuli assign a pixel more than once" );
    }
    ownership( coordinate.row(), coordinate.column() ) = region;
}

template <typename realT, class derotFunctObj, class verboseT>
void P4Reduction<realT, derotFunctObj, verboseT>::fitDetectorSearch(
    P4PCAResult &result,
    P4PCA::matrixT &predictors,
    P4PCA::vectorT &target,
    P4PCA::matrixT *coefficients,
    P4PCA::matrixT *probeResiduals,
    const pixelGridT &grid,
    std::size_t search,
    const std::vector<std::vector<int>> &selections,
    const std::vector<P4PixelCoordinate> &temporalOffsets,
    const std::vector<int> &modes,
    P4PCA::workspaceT &workspace,
    detail::P4PCAMixedWorkspace &mixedWorkspace,
    P4PCADowndateWorkspace *downdateWorkspace,
    P4PCATiming &timing,
    double &sameImageSamplingSeconds,
    double &temporalSamplingSeconds,
    const P4TargetExclusions *exclusions,
    const P4PCA::matrixT *probePredictors,
    const P4PCA::vectorT *probeTarget,
    const P4TrialSource *trialSource ) const
{
    if( selections.empty() || modes.empty() || search >= grid.searchPixelCount() ||
        !grid.searchPixel( search ).valid() )
    {
        throw std::invalid_argument( "P4 detector fit requires valid geometry, temporal rows, and modes" );
    }
    const std::size_t basePredictorCount = grid.predictorCount();
    const std::size_t temporalImageCount = selections.front().size() - 1;
    const std::size_t predictorCount = basePredictorCount + temporalImageCount * temporalOffsets.size();
    predictors.resize( static_cast<Eigen::Index>( selections.size() ), static_cast<Eigen::Index>( predictorCount ) );
    target.resize( static_cast<Eigen::Index>( selections.size() ) );

    const P4PixelCoordinate &coordinate = grid.searchPixel( search ).coordinate();
    sameImageSamplingSeconds = 0;
    temporalSamplingSeconds = 0;
    for( std::size_t targetIndex = 0; targetIndex < selections.size(); ++targetIndex )
    {
        const std::vector<int> &selection = selections[targetIndex];
        if( selection.size() != temporalImageCount + 1 )
        {
            throw std::invalid_argument( "P4 detector temporal selections must have a uniform width" );
        }
        const int centralImage = selection[0];
        const double sameImageSamplingBegin = omp_get_wtime();
        double targetValue =
            static_cast<double>( this->m_tgtIms.image( centralImage )( coordinate.row(), coordinate.column() ) );
        if( trialSource )
        {
            targetValue += checkedPredictorPromotion(
                trialSource->value( static_cast<std::size_t>( centralImage ), coordinate.row(), coordinate.column() ) );
        }
        if( !mx::math::isFinite( targetValue ) )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidarg,
                                           "P4 trial-perturbed detector target is non-finite" );
        }
        target( static_cast<Eigen::Index>( targetIndex ) ) = targetValue;
        for( std::size_t predictor = 0; predictor < basePredictorCount; ++predictor )
        {
            double predictorValue =
                checkedPredictorPromotion( grid.sample( this->m_tgtIms.image( centralImage ), search, predictor ) );
            if( trialSource )
            {
                predictorValue +=
                    checkedPredictorPromotion( trialSource->sample( static_cast<std::size_t>( centralImage ),
                                                                    grid.interpolation( search, predictor ) ) );
            }
            if( !mx::math::isFinite( predictorValue ) )
            {
                throw mx::exception<verboseT>( mx::error_t::invalidarg,
                                               "P4 trial-perturbed same-image predictor is non-finite" );
            }
            predictors( static_cast<Eigen::Index>( targetIndex ), static_cast<Eigen::Index>( predictor ) ) =
                predictorValue;
        }
        sameImageSamplingSeconds += omp_get_wtime() - sameImageSamplingBegin;

        const double temporalSamplingBegin = omp_get_wtime();
        for( std::size_t source = 1; source < selection.size(); ++source )
        {
            const int image = selection[source];
            for( std::size_t predictor = 0; predictor < temporalOffsets.size(); ++predictor )
            {
                const int row = coordinate.row() + temporalOffsets[predictor].row();
                const int column = coordinate.column() + temporalOffsets[predictor].column();
                double predictorValue = checkedPredictorPromotion( this->m_tgtIms.image( image )( row, column ) );
                if( trialSource )
                {
                    predictorValue += checkedPredictorPromotion(
                        trialSource->value( static_cast<std::size_t>( image ), row, column ) );
                }
                if( !mx::math::isFinite( predictorValue ) )
                {
                    throw mx::exception<verboseT>( mx::error_t::invalidarg,
                                                   "P4 trial-perturbed additional-image predictor is non-finite" );
                }
                const std::size_t columnIndex =
                    basePredictorCount + ( source - 1 ) * temporalOffsets.size() + predictor;
                predictors( static_cast<Eigen::Index>( targetIndex ), static_cast<Eigen::Index>( columnIndex ) ) =
                    predictorValue;
            }
        }
        temporalSamplingSeconds += omp_get_wtime() - temporalSamplingBegin;
    }

    if( exclusions )
    {
        if( coefficients )
        {
            throw std::invalid_argument( "P4 held-out fitting does not provide shared predictor coefficients" );
        }
        if( m_exclusionSolver == P4ExclusionSolver::factorDowndateExact )
        {
            if( !downdateWorkspace )
            {
                throw std::invalid_argument( "P4 exact factor deletion requires a worker-private workspace" );
            }
            if( probeResiduals || probePredictors || probeTarget )
            {
                if( !probeResiduals || !probePredictors || !probeTarget )
                {
                    throw std::invalid_argument( "P4 held-out fitting requires a complete frozen-probe request" );
                }
                P4PCA::calculateHeldOutProbeDowndated( result,
                                                       *probeResiduals,
                                                       predictors,
                                                       target,
                                                       *probePredictors,
                                                       *probeTarget,
                                                       *exclusions,
                                                       modes,
                                                       m_rankTolerance,
                                                       m_deletionBackend,
                                                       workspace,
                                                       *downdateWorkspace,
                                                       &timing );
            }
            else
            {
                P4PCA::calculateHeldOutDowndated( result,
                                                  predictors,
                                                  target,
                                                  *exclusions,
                                                  modes,
                                                  m_rankTolerance,
                                                  m_deletionBackend,
                                                  workspace,
                                                  *downdateWorkspace,
                                                  &timing );
            }
        }
        else
        {
            if( probeResiduals || probePredictors || probeTarget )
            {
                if( !probeResiduals || !probePredictors || !probeTarget )
                {
                    throw std::invalid_argument( "P4 held-out fitting requires a complete frozen-probe request" );
                }
#ifdef HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION
                p4ReductionDispatchHeldOutProbe( result,
                                                 *probeResiduals,
                                                 predictors,
                                                 target,
                                                 *probePredictors,
                                                 *probeTarget,
                                                 *exclusions,
                                                 modes,
                                                 m_rankTolerance,
                                                 workspace,
                                                 mixedWorkspace,
                                                 &timing );
#else
                static_cast<void>( workspace );
                detail::p4PCACalculateHeldOutProbeMixed( result,
                                                         *probeResiduals,
                                                         predictors,
                                                         target,
                                                         *probePredictors,
                                                         *probeTarget,
                                                         *exclusions,
                                                         modes,
                                                         m_rankTolerance,
                                                         mixedWorkspace,
                                                         &timing );
#endif
            }
            else
            {
#ifdef HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION
                p4ReductionDispatchHeldOut( result,
                                            predictors,
                                            target,
                                            *exclusions,
                                            modes,
                                            m_rankTolerance,
                                            workspace,
                                            mixedWorkspace,
                                            &timing );
#else
                static_cast<void>( workspace );
                detail::p4PCACalculateHeldOutMixed(
                    result, predictors, target, *exclusions, modes, m_rankTolerance, mixedWorkspace, &timing );
#endif
            }
        }
    }
    else
    {
        if( probeResiduals || probePredictors || probeTarget )
        {
            throw std::invalid_argument( "P4 frozen-probe fitting requires target exclusions" );
        }
#ifdef HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION
        p4ReductionDispatchDirect( result,
                                   predictors,
                                   target,
                                   modes,
                                   m_rankTolerance,
                                   workspace,
                                   mixedWorkspace,
                                   &timing,
                                   coefficients );
#else
        static_cast<void>( workspace );
        detail::p4PCACalculateMixed(
            result, predictors, target, modes, m_rankTolerance, mixedWorkspace, &timing, coefficients );
#endif
    }
}

template <typename realT, class derotFunctObj, class verboseT>
std::vector<realT> P4Reduction<realT, derotFunctObj, verboseT>::localTrialScales() const
{
    std::vector<realT> scales( static_cast<std::size_t>( this->m_Nims ), static_cast<realT>( 1 ) );
    if( this->m_fakeScaleFileName.empty() )
    {
        return scales;
    }
    if( this->m_fileList.size() != static_cast<std::size_t>( this->m_Nims ) )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                       "local fake scale lookup requires one filename per target image" );
    }

    std::vector<std::string> scaleFileNames;
    std::vector<realT> scaleValues;
    const mx::error_t readResult = mx::ioutils::readColumns( this->m_fakeScaleFileName, scaleFileNames, scaleValues );
    if( readResult != mx::error_t::noerror )
    {
        throw mx::exception<verboseT>( readResult, "reading local fake scale file" );
    }
    if( scaleFileNames.size() != scaleValues.size() )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                       "local fake scale file must contain filename and scale columns" );
    }

    std::map<std::string, realT> scaleByFile;
    for( std::size_t index = 0; index < scaleFileNames.size(); ++index )
    {
        const std::string name = mx::ioutils::pathFilename( scaleFileNames[index].c_str() );
        if( !mx::math::isFinite( scaleValues[index] ) || !scaleByFile.emplace( name, scaleValues[index] ).second )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "local fake scale entries must be finite and have unique filenames" );
        }
    }
    for( std::size_t image = 0; image < scales.size(); ++image )
    {
        const std::string name = mx::ioutils::pathFilename( this->m_fileList[image].c_str() );
        const auto found = scaleByFile.find( name );
        if( found == scaleByFile.end() )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                           "target filename is missing from the local fake scale file: " + name );
        }
        scales[image] = found->second;
    }
    return scales;
}

template <typename realT, class derotFunctObj, class verboseT>
int P4Reduction<realT, derotFunctObj, verboseT>::processLocalTrial(
    const std::vector<pixelGridT> &grids,
    const std::vector<P4PixelCoordinate> &temporalOffsets,
    const std::vector<double> &derotationAngles,
    const std::vector<P4TargetExclusions> &regionExclusions )
{
#ifdef HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION
    const std::optional<P4ReductionPCADispatch> experimentalDispatch = p4ReductionSelectedPCADispatch();
#endif
    if( grids.size() != m_regionStatistics.size() || grids.size() != m_temporalSelections.size() ||
        grids.size() != m_realizedModes.size() || grids.size() != regionExclusions.size() ||
        derotationAngles.size() != static_cast<std::size_t>( this->m_Nims ) )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr,
                                       "P4 local processing state does not match the configured reduction" );
    }

    const double localGeometryBegin = omp_get_wtime();
    P4LocalGeometry::lookupImageT ownership = m_ownership.template cast<std::int64_t>();
    P4LocalGeometry::lookupImageT searchIndexLookup =
        P4LocalGeometry::lookupImageT::Constant( this->m_Nrows, this->m_Ncols, -1 );
    for( std::size_t region = 0; region < grids.size(); ++region )
    {
        for( std::size_t search = 0; search < grids[region].searchPixelCount(); ++search )
        {
            const P4PixelCoordinate &coordinate = grids[region].searchPixel( search ).coordinate();
            if( search > static_cast<std::size_t>( std::numeric_limits<std::int64_t>::max() ) )
            {
                throw mx::exception<verboseT>( mx::error_t::sizeerr, "P4 local search index exceeds lookup range" );
            }
            searchIndexLookup( coordinate.row(), coordinate.column() ) = static_cast<std::int64_t>( search );
        }
    }

    const double imageCenterRow = 0.5 * static_cast<double>( this->m_Nrows - 1 );
    const double imageCenterColumn = 0.5 * static_cast<double>( this->m_Ncols - 1 );
    const realT positionAngleRadians = mx::math::dtor( -this->m_fakePA[0] );
    const realT separation = this->m_fakeSep[0];
    const double sourceRow = imageCenterRow + static_cast<double>( separation * std::sin( positionAngleRadians ) );
    const double sourceColumn =
        imageCenterColumn + static_cast<double>( separation * std::cos( positionAngleRadians ) );

    P4LocalGeometry geometry;
    try
    {
        geometry.configure( this->m_Nrows,
                            this->m_Ncols,
                            m_localStampSize,
                            sourceRow,
                            sourceColumn,
                            derotationAngles,
                            true,
                            ownership,
                            searchIndexLookup );
    }
    catch( ... )
    {
        std::throw_with_nested(
            mx::exception<verboseT>( mx::error_t::invalidconfig, "could not construct P4 local geometry" ) );
    }

    imageT sourceTemplate;
    mx::fits::fitsFile<realT, verboseT> reader;
    const mx::error_t readResult = reader.read( sourceTemplate, this->m_fakeFileName );
    if( readResult != mx::error_t::noerror )
    {
        throw mx::exception<verboseT>( readResult, "could not read local fake template " + this->m_fakeFileName );
    }
    P4TrialSource trialSource;
    try
    {
        const std::vector<realT> trialScales = localTrialScales();
        trialSource.configure( sourceTemplate,
                               this->m_Nrows,
                               this->m_Ncols,
                               m_localStampSize,
                               derotationAngles,
                               separation,
                               static_cast<double>( this->m_fakePA[0] ),
                               static_cast<double>( this->m_fakeContrast[0] ),
                               trialScales );
        if( this->m_subtractPlanet )
        {
            for( std::size_t planet = 0; planet < this->m_planetSep.size(); ++planet )
            {
                trialSource.addSource( derotationAngles,
                                       static_cast<double>( this->m_planetSep[planet] ),
                                       static_cast<double>( this->m_planetPA[planet] ),
                                       -static_cast<double>( this->m_planetContrast[planet] ),
                                       trialScales );
            }
        }
    }
    catch( ... )
    {
        std::throw_with_nested(
            mx::exception<verboseT>( mx::error_t::invalidconfig, "could not prepare P4 local trial source" ) );
    }

    m_localOriginRow = geometry.originRow();
    m_localOriginColumn = geometry.originColumn();
    m_localSourceRow = geometry.sourceRow();
    m_localSourceColumn = geometry.sourceColumn();
    m_localTemplateRows = trialSource.cropRows();
    m_localTemplateColumns = trialSource.cropColumns();
    m_localSearchCount = geometry.searchRequests().size();
    m_localSparseSampleCount = geometry.sparseSampleCount();
    m_localGeometryBytes = geometry.storageBytes();
    m_timing.geometryElapsedSeconds += omp_get_wtime() - localGeometryBegin;

    using localResidualT = Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic>;
    using localValidityT = Eigen::Array<std::uint8_t, Eigen::Dynamic, Eigen::Dynamic>;
    std::vector<localResidualT> localResiduals( m_localSearchCount );
    std::vector<localValidityT> localValidity( m_localSearchCount );
    for( std::size_t request = 0; request < m_localSearchCount; ++request )
    {
        const std::size_t frameCount = geometry.searchRequests()[request].frames().size();
        if( frameCount > static_cast<std::size_t>( std::numeric_limits<Eigen::Index>::max() ) ||
            m_modeFractions.size() > static_cast<std::size_t>( std::numeric_limits<Eigen::Index>::max() ) )
        {
            throw mx::exception<verboseT>( mx::error_t::sizeerr, "P4 local residual dimensions exceed Eigen range" );
        }
        localResiduals[request].resize( static_cast<Eigen::Index>( frameCount ),
                                        static_cast<Eigen::Index>( m_modeFractions.size() ) );
        localResiduals[request].setZero();
        localValidity[request].resize( static_cast<Eigen::Index>( frameCount ),
                                       static_cast<Eigen::Index>( m_modeFractions.size() ) );
        localValidity[request].setZero();
    }

    const long double compactBytes = static_cast<long double>( m_localSparseSampleCount ) *
                                     static_cast<long double>( m_modeFractions.size() ) *
                                     ( sizeof( realT ) + sizeof( std::uint8_t ) );
    const long double templateBytes =
        static_cast<long double>( this->m_Nrows ) * static_cast<long double>( this->m_Ncols ) * sizeof( realT );
    const long double materializationBytes =
        2.0L * static_cast<long double>( m_localStampSize ) * static_cast<long double>( m_localStampSize ) *
        static_cast<long double>( this->m_Nims ) * static_cast<long double>( m_modeFractions.size() ) * sizeof( realT );
    if( compactBytes > static_cast<long double>( std::numeric_limits<std::size_t>::max() ) ||
        templateBytes > static_cast<long double>( std::numeric_limits<std::size_t>::max() ) ||
        materializationBytes > static_cast<long double>( std::numeric_limits<std::size_t>::max() ) )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "P4 local memory estimate exceeds size_t range" );
    }
    m_compactResidualBytes = static_cast<std::size_t>( compactBytes );
    m_psfModelBytes = static_cast<std::size_t>( templateBytes );
    m_materializationBytes = static_cast<std::size_t>( materializationBytes );

    if( m_memoryFraction != 0 )
    {
        const std::optional<std::size_t> availableMemory = p4AvailableMemoryBytes();
        if( availableMemory )
        {
            m_availableMemoryBytes = *availableMemory;
            m_memoryBudgetBytes = static_cast<std::size_t>(
                std::floor( static_cast<long double>( m_memoryFraction ) * m_availableMemoryBytes ) );
            if( m_memoryBudgetBytes == 0 )
            {
                throw mx::exception<verboseT>( mx::error_t::allocerr,
                                               "p4.memoryFraction selects a zero-byte automatic memory budget" );
            }
        }
        else
        {
            std::cerr << "WARNING: P4 could not read Linux MemAvailable; automatic worker limiting is disabled\n";
        }
    }

    if( m_compactResidualBytes > std::numeric_limits<std::size_t>::max() - m_psfModelBytes ||
        m_compactResidualBytes + m_psfModelBytes > std::numeric_limits<std::size_t>::max() - m_localGeometryBytes ||
        m_compactResidualBytes + m_psfModelBytes + m_localGeometryBytes >
            std::numeric_limits<std::size_t>::max() - m_targetExclusionBytes )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "P4 local persistent byte count overflow" );
    }
    const std::size_t persistentBytes =
        m_compactResidualBytes + m_psfModelBytes + m_localGeometryBytes + m_targetExclusionBytes;
    if( m_memoryBudgetBytes != 0 &&
        ( persistentBytes > m_memoryBudgetBytes || m_materializationBytes > m_memoryBudgetBytes - persistentBytes ) )
    {
        throw mx::exception<verboseT>(
            mx::error_t::allocerr,
            "P4 local residuals and output materialization exceed the automatic memory budget; increase "
            "p4.memoryFraction, reduce p4.localStampSize, or set p4.memoryFraction=0 to disable the budget" );
    }

    std::vector<std::vector<int>> targetRowByFrame( grids.size(),
                                                    std::vector<int>( static_cast<std::size_t>( this->m_Nims ), -1 ) );
    for( std::size_t region = 0; region < grids.size(); ++region )
    {
        for( std::size_t targetIndex = 0; targetIndex < m_temporalSelections[region].size(); ++targetIndex )
        {
            const int frame = m_temporalSelections[region][targetIndex][0];
            targetRowByFrame[region][static_cast<std::size_t>( frame )] = static_cast<int>( targetIndex );
        }
    }

    std::size_t maximumWorkerBytes{ 0 };
    for( std::size_t region = 0; region < grids.size(); ++region )
    {
        P4RegionStatistics &statistics = m_regionStatistics[region];
        const std::size_t estimatedMaximumExcludedRows =
            m_exclusionSolver == P4ExclusionSolver::factorDowndateExact &&
                    m_deletionBackend == mx::math::svdDeletionBackend::rankOneSecular
                ? 1
                : statistics.maximumExcludedRows;
        statistics.estimatedWorkerBytes =
            estimatedWorkerBytes( statistics.targetImageCount,
                                  statistics.predictorCount,
                                  m_modeFractions.size(),
                                  false,
                                  0,
                                  m_excludeMethod != HCI::exclude::none,
                                  m_exclusionSolver,
                                  estimatedMaximumExcludedRows,
                                  static_cast<std::size_t>( m_realizedModes[region].back() ) );
        maximumWorkerBytes = std::max( maximumWorkerBytes, statistics.estimatedWorkerBytes );
    }
    const int requestedWorkers =
        std::max( 1,
                  std::min( omp_get_max_threads(),
                            static_cast<int>( std::min<std::size_t>(
                                std::max<std::size_t>( 1, m_localSearchCount ),
                                static_cast<std::size_t>( std::numeric_limits<int>::max() ) ) ) ) );
    int effectiveWorkers = requestedWorkers;
    if( m_memoryBudgetBytes != 0 )
    {
        effectiveWorkers =
            memoryLimitedWorkerCount( requestedWorkers, m_memoryBudgetBytes, persistentBytes, maximumWorkerBytes );
        if( effectiveWorkers == 0 )
        {
            throw mx::exception<verboseT>(
                mx::error_t::allocerr,
                "one P4 local worker does not fit the automatic memory budget; increase p4.memoryFraction or set "
                "p4.memoryFraction=0 to disable the budget" );
        }
    }
    for( std::size_t region = 0; region < grids.size(); ++region )
    {
        P4RegionStatistics &statistics = m_regionStatistics[region];
        statistics.maximumWorkerCount = requestedWorkers;
        statistics.effectiveWorkerCount = effectiveWorkers;
        statistics.validLocalFitCount = 0;
        statistics.maskedLocalFitCount = 0;
        statistics.supportInvalidLocalFitCount = 0;
        statistics.minimumNumericalRank = 0;
        statistics.minimumBaseRank = 0;
        statistics.downdateClampCount = 0;
        statistics.explicitFallbackCount = 0;
        statistics.rankBoundaryFallbackPixelCount = 0;
        statistics.factorValidationFallbackPixelCount = 0;
        statistics.deletionSolverFallbackPixelCount = 0;
        statistics.maximumFactorOrthogonalityDefect = 0;
        statistics.factorOrthogonalityToleranceAtMaximumDefect = 0;
        statistics.rankInvalidCounts.assign( m_modeFractions.size(), 0 );
    }

    constexpr double bytesPerMiB = 1024.0 * 1024.0;
    std::cerr << "P4 local geometry: " << m_localStampSize << 'x' << m_localStampSize << " sky pixels, "
              << m_localSearchCount << " unique detector fits, " << m_localSparseSampleCount
              << " sparse residual samples\n";
    std::cerr << "P4 local memory: " << static_cast<double>( m_compactResidualBytes ) / bytesPerMiB
              << " MiB residuals, " << static_cast<double>( m_localGeometryBytes ) / bytesPerMiB
              << " MiB sparse geometry, " << static_cast<double>( m_targetExclusionBytes ) / bytesPerMiB
              << " MiB target exclusions, " << static_cast<double>( m_materializationBytes ) / bytesPerMiB
              << " MiB output materialization, workers " << effectiveWorkers << " / " << requestedWorkers << '\n';

    std::vector<int> requestRanks( m_localSearchCount, 0 );
    std::vector<int> requestBaseRanks( m_localSearchCount, 0 );
    std::vector<std::size_t> requestDowndateClampCounts( m_localSearchCount, 0 );
    std::vector<std::size_t> requestExplicitFallbackCounts( m_localSearchCount, 0 );
    std::vector<P4PCAFallbackReason> requestFallbackReasons( m_localSearchCount, P4PCAFallbackReason::none );
    std::vector<double> requestFactorDefects( m_localSearchCount, 0 );
    std::vector<double> requestFactorTolerances( m_localSearchCount, 0 );
    std::vector<std::uint8_t> requestAttempted( m_localSearchCount, 0 );
    std::exception_ptr workerException;
    std::atomic<bool> failed{ false };
    const imageT *mask = this->m_mask.size() == 0 ? nullptr : &this->m_mask;
    const double regressionBegin = omp_get_wtime();

    // clang-format off
#pragma omp parallel num_threads(effectiveWorkers)
    // clang-format on
    {
        P4PCA::workspaceT workspace;
        detail::P4PCAMixedWorkspace mixedWorkspace;
#ifdef HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION
        std::optional<detail::P4PCAExperimentalWorkspace> experimentalWorkspace;
        std::optional<P4ReductionWorkerPrecisionScope> experimentalPrecisionScope;
        if( experimentalDispatch )
        {
            experimentalWorkspace.emplace();
            experimentalPrecisionScope.emplace( *experimentalDispatch, *experimentalWorkspace );
        }
#endif
        P4PCADowndateWorkspace downdateWorkspace;
        P4PCA::matrixT predictors;
        P4PCA::vectorT target;
        P4PCAResult result;
        double threadSamplingSeconds{ 0 };
        double threadSameImageSamplingSeconds{ 0 };
        double threadTemporalSamplingSeconds{ 0 };
        double threadGramSeconds{ 0 };
        double threadEigensolveSeconds{ 0 };
        double threadBaseFactorSeconds{ 0 };
        double threadDeletionSeconds{ 0 };
        double threadExplicitFallbackSeconds{ 0 };
        double threadProjectionSeconds{ 0 };

        // clang-format off
#pragma omp for schedule(static)
        // clang-format on
        for( std::size_t requestIndex = 0; requestIndex < m_localSearchCount; ++requestIndex )
        {
            if( failed.load( std::memory_order_acquire ) )
            {
                continue;
            }
            try
            {
                const P4LocalSearchRequest &request = geometry.searchRequests()[requestIndex];
                const std::size_t region = static_cast<std::size_t>( request.region() );
                const pixelGridT &grid = grids[region];
                const bool temporalValid = m_regionStatistics[region].temporalNumberImages == 0 ||
                                           p4TemporalPredictorsValid( request.coordinate(),
                                                                      temporalOffsets,
                                                                      this->m_Nrows,
                                                                      this->m_Ncols,
                                                                      mask );
                if( !grid.searchPixel( request.searchIndex() ).valid() || !temporalValid )
                {
                    continue;
                }

                double sameImageSamplingSeconds{ 0 };
                double temporalSamplingSeconds{ 0 };
                P4PCATiming pcaTiming;
                const P4TargetExclusions *exclusions =
                    regionExclusions[region].empty() ? nullptr : &regionExclusions[region];
                fitDetectorSearch( result,
                                   predictors,
                                   target,
                                   nullptr,
                                   nullptr,
                                   grid,
                                   request.searchIndex(),
                                   m_temporalSelections[region],
                                   temporalOffsets,
                                   m_realizedModes[region],
                                   workspace,
                                   mixedWorkspace,
                                   &downdateWorkspace,
                                   pcaTiming,
                                   sameImageSamplingSeconds,
                                   temporalSamplingSeconds,
                                   exclusions,
                                   nullptr,
                                   nullptr,
                                   &trialSource );
                threadSamplingSeconds += sameImageSamplingSeconds + temporalSamplingSeconds;
                threadSameImageSamplingSeconds += sameImageSamplingSeconds;
                threadTemporalSamplingSeconds += temporalSamplingSeconds;
                threadGramSeconds += pcaTiming.gramWorkerSeconds;
                threadEigensolveSeconds += pcaTiming.eigensolveWorkerSeconds;
                threadBaseFactorSeconds += pcaTiming.baseFactorWorkerSeconds;
                threadDeletionSeconds += pcaTiming.deletionWorkerSeconds;
                threadExplicitFallbackSeconds += pcaTiming.explicitFallbackWorkerSeconds;
                threadProjectionSeconds += pcaTiming.projectionWorkerSeconds;
                requestAttempted[requestIndex] = 1;
                requestRanks[requestIndex] = result.numericalRank;
                requestBaseRanks[requestIndex] = result.baseRank;
                requestDowndateClampCounts[requestIndex] = result.downdateClampCount;
                requestExplicitFallbackCounts[requestIndex] = result.explicitFallbackCount;
                requestFallbackReasons[requestIndex] = result.explicitFallbackReason;
                requestFactorDefects[requestIndex] = result.factorOrthogonalityDefect;
                requestFactorTolerances[requestIndex] = result.factorOrthogonalityTolerance;

                const double residualApplyBegin = omp_get_wtime();
                for( std::size_t output = 0; output < m_modeFractions.size(); ++output )
                {
                    for( std::size_t frameOffset = 0; frameOffset < request.frames().size(); ++frameOffset )
                    {
                        const int frame = request.frames()[frameOffset];
                        const int targetRow = targetRowByFrame[region][static_cast<std::size_t>( frame )];
                        if( targetRow < 0 || !result.sampleSupported( static_cast<Eigen::Index>( targetRow ), output ) )
                        {
                            continue;
                        }
                        localResiduals[requestIndex]( static_cast<Eigen::Index>( frameOffset ),
                                                      static_cast<Eigen::Index>( output ) ) =
                            checkedResidualCast( result.residuals( targetRow, static_cast<Eigen::Index>( output ) ) );
                        localValidity[requestIndex]( static_cast<Eigen::Index>( frameOffset ),
                                                     static_cast<Eigen::Index>( output ) ) = 1;
                    }
                }
                threadProjectionSeconds += omp_get_wtime() - residualApplyBegin;
            }
            catch( ... )
            {
                // clang-format off
#pragma omp critical(P4ReductionLocalException)
                // clang-format on
                {
                    if( !workerException )
                    {
                        workerException = std::current_exception();
                    }
                }
                failed.store( true, std::memory_order_release );
            }
        }

        // clang-format off
#pragma omp critical(P4ReductionLocalTiming)
        // clang-format on
        {
            m_timing.samplingWorkerSeconds += threadSamplingSeconds;
            m_timing.sameImageSamplingWorkerSeconds += threadSameImageSamplingSeconds;
            m_timing.temporalSamplingWorkerSeconds += threadTemporalSamplingSeconds;
            m_timing.gramWorkerSeconds += threadGramSeconds;
            m_timing.eigensolveWorkerSeconds += threadEigensolveSeconds;
            m_timing.baseFactorWorkerSeconds += threadBaseFactorSeconds;
            m_timing.deletionWorkerSeconds += threadDeletionSeconds;
            m_timing.explicitFallbackWorkerSeconds += threadExplicitFallbackSeconds;
            m_timing.projectionWorkerSeconds += threadProjectionSeconds;
        }
    }

    if( workerException )
    {
        try
        {
            std::rethrow_exception( workerException );
        }
        catch( ... )
        {
            std::throw_with_nested( mx::exception<verboseT>( mx::error_t::exception, "P4 local worker failed" ) );
        }
    }
    m_timing.regressionElapsedSeconds = omp_get_wtime() - regressionBegin;

    for( std::size_t requestIndex = 0; requestIndex < m_localSearchCount; ++requestIndex )
    {
        const std::size_t region = static_cast<std::size_t>( geometry.searchRequests()[requestIndex].region() );
        P4RegionStatistics &statistics = m_regionStatistics[region];
        if( requestAttempted[requestIndex] == 0 )
        {
            ++statistics.maskedLocalFitCount;
            continue;
        }
        ++statistics.validLocalFitCount;
        if( statistics.validLocalFitCount == 1 )
        {
            statistics.minimumNumericalRank = requestRanks[requestIndex];
            if( m_exclusionSolver == P4ExclusionSolver::factorDowndateExact )
            {
                statistics.minimumBaseRank = requestBaseRanks[requestIndex];
            }
        }
        else
        {
            statistics.minimumNumericalRank = std::min( statistics.minimumNumericalRank, requestRanks[requestIndex] );
            if( m_exclusionSolver == P4ExclusionSolver::factorDowndateExact )
            {
                statistics.minimumBaseRank = std::min( statistics.minimumBaseRank, requestBaseRanks[requestIndex] );
            }
        }
        if( m_exclusionSolver == P4ExclusionSolver::factorDowndateExact )
        {
            statistics.downdateClampCount += requestDowndateClampCounts[requestIndex];
            statistics.explicitFallbackCount += requestExplicitFallbackCounts[requestIndex];
            switch( requestFallbackReasons[requestIndex] )
            {
            case P4PCAFallbackReason::none:
                break;
            case P4PCAFallbackReason::rankBoundary:
                ++statistics.rankBoundaryFallbackPixelCount;
                break;
            case P4PCAFallbackReason::factorValidation:
                ++statistics.factorValidationFallbackPixelCount;
                if( statistics.factorValidationFallbackPixelCount == 1 ||
                    requestFactorDefects[requestIndex] > statistics.maximumFactorOrthogonalityDefect )
                {
                    statistics.maximumFactorOrthogonalityDefect = requestFactorDefects[requestIndex];
                    statistics.factorOrthogonalityToleranceAtMaximumDefect = requestFactorTolerances[requestIndex];
                }
                break;
            case P4PCAFallbackReason::deletionSolver:
                ++statistics.deletionSolverFallbackPixelCount;
                break;
            }
        }
        for( std::size_t output = 0; output < m_modeFractions.size(); ++output )
        {
            bool modeValid{ false };
            for( Eigen::Index frameOffset = 0; frameOffset < localValidity[requestIndex].rows(); ++frameOffset )
            {
                modeValid =
                    modeValid || localValidity[requestIndex]( frameOffset, static_cast<Eigen::Index>( output ) );
            }
            if( !modeValid )
            {
                ++statistics.rankInvalidCounts[output];
            }
        }
    }

    for( std::size_t region = 0; region < grids.size(); ++region )
    {
        const P4RegionStatistics &statistics = m_regionStatistics[region];
        if( statistics.explicitFallbackCount != 0 )
        {
            std::ostringstream fallbackWarning;
            fallbackWarning << "WARNING: P4 local region " << region << " recomputed "
                            << statistics.explicitFallbackCount << " target rows with the explicit oracle at "
                            << statistics.rankBoundaryFallbackPixelCount << " rank-boundary search pixels and "
                            << statistics.factorValidationFallbackPixelCount << " factor-validation search pixels and "
                            << statistics.deletionSolverFallbackPixelCount << " deletion-solver search pixels";
            if( statistics.factorValidationFallbackPixelCount != 0 )
            {
                fallbackWarning << "; maximum factor defect " << std::setprecision( 17 )
                                << statistics.maximumFactorOrthogonalityDefect << " at tolerance "
                                << statistics.factorOrthogonalityToleranceAtMaximumDefect;
            }
            std::cerr << fallbackWarning.str() << '\n';
        }
        for( std::size_t output = 0; output < m_modeFractions.size(); ++output )
        {
            if( statistics.rankInvalidCounts[output] != 0 )
            {
                std::cerr << "WARNING: P4 local region " << region << " output " << output << " requested "
                          << m_realizedModes[region][output]
                          << " modes above numerical rank for every requested target at "
                          << statistics.rankInvalidCounts[output] << " search pixels\n";
            }
        }
    }

    this->t_derotate_begin = mx::sys::get_curr_time();
    this->m_psfsub.resize( m_modeFractions.size() );
    this->m_psfsubValidity.resize( m_modeFractions.size() );
    for( std::size_t output = 0; output < m_modeFractions.size(); ++output )
    {
        this->m_psfsub[output].resize( m_localStampSize, m_localStampSize, this->m_Nims );
        this->m_psfsub[output].setZero();
        this->m_psfsubValidity[output].resize( m_localStampSize, m_localStampSize, this->m_Nims );
        this->m_psfsubValidity[output].setZero();
    }

    for( int stampColumn = 0; stampColumn < m_localStampSize; ++stampColumn )
    {
        for( int stampRow = 0; stampRow < m_localStampSize; ++stampRow )
        {
            for( std::size_t frame = 0; frame < static_cast<std::size_t>( this->m_Nims ); ++frame )
            {
                const P4LocalOutputSample &sample = geometry.outputSample( stampRow, stampColumn, frame );
                if( !sample.valid() )
                {
                    continue;
                }
                for( std::size_t output = 0; output < m_modeFractions.size(); ++output )
                {
                    realT value{ 0 };
                    bool valid{ true };
                    for( const P4LocalResidualSample &dependency : sample.samples() )
                    {
                        if( localValidity[dependency.requestIndex()](
                                static_cast<Eigen::Index>( dependency.frameOffset() ),
                                static_cast<Eigen::Index>( output ) ) == 0 )
                        {
                            valid = false;
                            break;
                        }
                        value += localResiduals[dependency.requestIndex()](
                                     static_cast<Eigen::Index>( dependency.frameOffset() ),
                                     static_cast<Eigen::Index>( output ) ) *
                                 dependency.weight();
                    }
                    if( valid )
                    {
                        this->m_psfsub[output].image( static_cast<int>( frame ) )( stampRow, stampColumn ) = value;
                        this->m_psfsubValidity[output].image( static_cast<int>( frame ) )( stampRow, stampColumn ) = 1;
                    }
                }
            }
        }
    }
    this->t_derotate_end = mx::sys::get_curr_time();

    std::cerr << "combining local sky stamp\n";
    this->combineFinim();
    m_localFinalValidity.resize( m_localStampSize, m_localStampSize, static_cast<int>( m_modeFractions.size() ) );
    for( int output = 0; output < m_localFinalValidity.planes(); ++output )
    {
        for( int stampColumn = 0; stampColumn < m_localStampSize; ++stampColumn )
        {
            for( int stampRow = 0; stampRow < m_localStampSize; ++stampRow )
            {
                m_localFinalValidity.image( output )( stampRow, stampColumn ) =
                    mx::math::isFinite( this->m_finim.image( output )( stampRow, stampColumn ) ) ? 1 : 0;
            }
        }
    }
    this->materializePSFSubInvalid();

    fitsHeaderT algorithmHeader;
    this->stdFitsHeader( &algorithmHeader );
    appendReductionHeader( algorithmHeader );
    if( this->m_doWriteFinim )
    {
        const std::string finalImagePath = this->finalImageOutputPath();
        std::cerr << "writing\n";
        fitsHeaderT residualHeader = algorithmHeader;
        residualHeader.template append<std::string>( "P4 PRODUCT ROLE", "LOCAL_RESIDUAL", "P4 product role" );
        this->writeFinimAtPath( finalImagePath, &residualHeader );
        writeLocalValidity( finalImagePath, algorithmHeader );
    }
    this->t_end = mx::sys::get_curr_time();
    return 0;
}

template <typename realT, class derotFunctObj, class verboseT>
void P4Reduction<realT, derotFunctObj, verboseT>::writeLocalValidity( const std::string &finalImagePath,
                                                                      const fitsHeaderT &finalHeader )
{
    const std::string path = p4FilterDiagnosticPath( finalImagePath, "local_validity", !this->m_exactFinimName );
    const std::string parent = mx::ioutils::parentPath( path );
    if( !parent.empty() )
    {
        const mx::error_t directoryResult = mx::ioutils::createDirectories( parent );
        if( directoryResult != mx::error_t::noerror )
        {
            throw mx::exception<verboseT>( directoryResult, "could not create P4 local validity directory" );
        }
    }

    fitsHeaderT header;
    fitsHeaderT additionalHeader = finalHeader;
    this->finalImageHeader( header, &additionalHeader );
    header.template append<std::string>( "P4 PRODUCT ROLE", "LOCAL_VALIDITY", "P4 product role" );
    static std::atomic<unsigned long long> sequence{ 0 };
    const std::string temporaryPath =
        path + ".tmp." + std::to_string( sequence.fetch_add( 1, std::memory_order_relaxed ) );
    mx::fits::fitsFile<realT, verboseT> writer;
    const mx::error_t writeResult = writer.write( temporaryPath, m_localFinalValidity, header );
    if( writeResult != mx::error_t::noerror )
    {
        std::error_code cleanupError;
        std::filesystem::remove( temporaryPath, cleanupError );
        throw mx::exception<verboseT>( writeResult, "could not write P4 local validity product " + path );
    }
    std::error_code renameError;
    std::filesystem::rename( temporaryPath, path, renameError );
    if( renameError )
    {
        std::error_code cleanupError;
        std::filesystem::remove( temporaryPath, cleanupError );
        throw mx::exception<verboseT>( mx::error_t::fileoerr,
                                       "could not publish P4 local validity product " + path + ": " +
                                           renameError.message() );
    }
    std::cerr << "P4 local validity written to: " << path << '\n';
}

template <typename realT, class derotFunctObj, class verboseT>
int P4Reduction<realT, derotFunctObj, verboseT>::reduce()
{
    if( !( this->m_preProcess_only && !this->m_skipPreProcess ) )
    {
        validateConfiguration();
    }
    if( !this->m_RDIfileList.empty() || !this->m_RDIfileListFile.empty() || !this->m_RDIdirectory.empty() ||
        !this->m_RDIprefix.empty() || this->m_refIms.planes() != 0 )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                       "P4 initially supports target-only ADI and rejects RDI input" );
    }
    if( !this->m_filesRead && this->m_fileList.empty() )
    {
        const mx::error_t listResult = this->load_fileList();
        if( listResult != mx::error_t::noerror )
        {
            throw mx::exception<verboseT>( listResult, "could not discover P4 target input files" );
        }
    }
    return regions( m_minRadius, m_maxRadius );
}

template <typename realT, class derotFunctObj, class verboseT>
P4LocalEvaluation<realT> P4Reduction<realT, derotFunctObj, verboseT>::evaluateLocal( const P4LocalTrial &trial )
{
    if( m_localStampSize <= 0 )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                       "P4 local evaluation requires positive p4.localStampSize" );
    }
    if( !mx::math::isFinite( trial.separation ) || trial.separation < 0 || !mx::math::isFinite( trial.positionAngle ) ||
        !mx::math::isFinite( trial.contrast ) )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidarg,
                                       "P4 local trial separation, PA, and contrast must be finite and valid" );
    }

    const std::vector<realT> configuredSeparation = this->m_fakeSep;
    const std::vector<realT> configuredPositionAngle = this->m_fakePA;
    const std::vector<realT> configuredContrast = this->m_fakeContrast;
    const int configuredWriteFinal = this->m_doWriteFinim;
    const bool configuredWriteDiagnostics = m_writeDiagnostics;
    this->m_fakeSep = { static_cast<realT>( trial.separation ) };
    this->m_fakePA = { static_cast<realT>( trial.positionAngle ) };
    this->m_fakeContrast = { static_cast<realT>( trial.contrast ) };
    this->m_doWriteFinim = false;
    m_writeDiagnostics = false;

    try
    {
        const int result = reduce();
        if( result != 0 )
        {
            throw mx::exception<verboseT>( mx::error_t::exception,
                                           "P4 local evaluation returned a nonzero reduction status" );
        }
        P4LocalEvaluation<realT> evaluation;
        evaluation.residual = this->m_finim;
        evaluation.validity = m_localFinalValidity;
        evaluation.originRow = m_localOriginRow;
        evaluation.originColumn = m_localOriginColumn;
        evaluation.sourceRow = m_localSourceRow;
        evaluation.sourceColumn = m_localSourceColumn;
        evaluation.searchCount = m_localSearchCount;
        evaluation.sparseSampleCount = m_localSparseSampleCount;
        evaluation.elapsedSeconds = this->t_end - this->t_begin;
        evaluation.timing = m_timing;

        this->m_fakeSep = configuredSeparation;
        this->m_fakePA = configuredPositionAngle;
        this->m_fakeContrast = configuredContrast;
        this->m_doWriteFinim = configuredWriteFinal;
        m_writeDiagnostics = configuredWriteDiagnostics;
        return evaluation;
    }
    catch( ... )
    {
        this->m_fakeSep = configuredSeparation;
        this->m_fakePA = configuredPositionAngle;
        this->m_fakeContrast = configuredContrast;
        this->m_doWriteFinim = configuredWriteFinal;
        m_writeDiagnostics = configuredWriteDiagnostics;
        throw;
    }
}

template <typename realT, class derotFunctObj, class verboseT>
P4LocalEvaluation<realT>
P4Reduction<realT, derotFunctObj, verboseT>::evaluateLocal( const P4LocalTrial &trial,
                                                            const std::vector<int> &includedFrames )
{
    if( includedFrames.empty() || !std::is_sorted( includedFrames.begin(), includedFrames.end(), std::less<int>() ) ||
        std::adjacent_find( includedFrames.begin(), includedFrames.end() ) != includedFrames.end() ||
        includedFrames.front() < 0 )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidarg,
                                       "P4 local included frames must be nonempty, unique, and strictly increasing" );
    }
    const std::optional<std::vector<int>> configuredIncludedFrames = m_localIncludedFrames;
    m_localIncludedFrames = includedFrames;
    try
    {
        P4LocalEvaluation<realT> evaluation = evaluateLocal( trial );
        m_localIncludedFrames = configuredIncludedFrames;
        return evaluation;
    }
    catch( ... )
    {
        m_localIncludedFrames = configuredIncludedFrames;
        throw;
    }
}

template <typename realT, class derotFunctObj, class verboseT>
std::size_t P4Reduction<realT, derotFunctObj, verboseT>::targetFrameCount() const
{
    return this->m_Nims > 0 ? static_cast<std::size_t>( this->m_Nims ) : 0;
}

template <typename realT, class derotFunctObj, class verboseT>
std::string P4Reduction<realT, derotFunctObj, verboseT>::auxiliaryOutputDirectory() const
{
    return p4AuxiliaryProductDirectory( this->finalImageOutputPath() ).string();
}

template <typename realT, class derotFunctObj, class verboseT>
int P4Reduction<realT, derotFunctObj, verboseT>::regions( const std::vector<realT> &minimumRadii,
                                                          const std::vector<realT> &maximumRadii )
{
#ifdef HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION
    const std::optional<P4ReductionPCADispatch> experimentalDispatch = p4ReductionSelectedPCADispatch();
#endif
    m_minRadius = minimumRadii;
    m_maxRadius = maximumRadii;
    if( !( this->m_preProcess_only && !this->m_skipPreProcess ) )
    {
        validateConfiguration();
    }

    if( !this->m_RDIfileList.empty() || !this->m_RDIfileListFile.empty() || !this->m_RDIdirectory.empty() ||
        !this->m_RDIprefix.empty() || this->m_refIms.planes() != 0 )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                       "P4 initially supports target-only ADI and rejects RDI input" );
    }

    m_automaticOutputSize = 0;
    m_automaticDerotationSize = 0;
    if( !( this->m_preProcess_only && !this->m_skipPreProcess ) && this->m_imSize == 0 )
    {
        m_automaticOutputSize =
            2 * ( *std::max_element( m_maxRadius.begin(), m_maxRadius.end() ) + p4AutomaticCropPadding );

        std::cerr << "set final output size based on regions to " << m_automaticOutputSize << '\n';
    }

    this->t_begin = mx::sys::get_curr_time();
    this->m_psfsub.clear();
    this->m_psfsubValidity.clear();
    m_realizedModes.clear();
    m_regionStatistics.clear();
    m_temporalSelections.clear();
    m_localPSFModels.clear();
    m_localPSFTemporalCoefficients.clear();
    m_localPSFValidity.clear();
    m_localPSFComponentCounts.clear();
    m_localPSFRows = 0;
    m_localPSFColumns = 0;
    m_psfTemplateRows = 0;
    m_psfTemplateColumns = 0;
    m_ownership.resize( 0, 0 );
    m_availableMemoryBytes = 0;
    m_memoryBudgetBytes = 0;
    m_compactResidualBytes = 0;
    m_targetExclusionBytes = 0;
    m_materializationBytes = 0;
    m_localPSFBytes = 0;
    m_psfModeBatchSize = 0;
    m_psfDowndateClampCount = 0;
    m_psfExplicitFallbackCount = 0;
    m_psfRankBoundaryFallbackPixelCount = 0;
    m_psfFactorValidationFallbackPixelCount = 0;
    m_psfDeletionSolverFallbackPixelCount = 0;
    m_psfMaximumFactorOrthogonalityDefect = 0;
    m_psfFactorOrthogonalityToleranceAtMaximumDefect = 0;
    m_psfModelBytes = 0;
    m_psfReconstructionBytes = 0;
    m_psfFilterBytes = 0;
    m_localFinalValidity.resize( 0, 0, 0 );
    m_localOriginRow = 0;
    m_localOriginColumn = 0;
    m_localSourceRow = 0;
    m_localSourceColumn = 0;
    m_localTemplateRows = 0;
    m_localTemplateColumns = 0;
    m_localSearchCount = 0;
    m_localSparseSampleCount = 0;
    m_localGeometryBytes = 0;
    m_timing.reset();

    if( !this->m_filesRead )
    {
        try
        {
            this->readFiles();
        }
        catch( ... )
        {
            std::throw_with_nested( mx::exception<verboseT>( mx::error_t::exception, "from P4 target loading" ) );
        }
    }

    if( this->m_Nims <= 0 || this->m_Nrows <= 0 || this->m_Ncols <= 0 || this->m_tgtIms.planes() != this->m_Nims ||
        this->m_tgtIms.rows() != this->m_Nrows || this->m_tgtIms.cols() != this->m_Ncols )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "invalid target cube state for P4 reduction" );
    }
    if( m_automaticOutputSize != 0 )
    {
        m_automaticDerotationSize = m_automaticOutputSize + ( ( this->m_Nrows - m_automaticOutputSize ) & 1 );
    }
    std::vector<int> activeFrames;
    if( m_localIncludedFrames )
    {
        activeFrames = *m_localIncludedFrames;
        if( activeFrames.empty() || activeFrames.back() >= this->m_Nims )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidarg,
                                           "P4 local included frame index is outside the loaded target sequence" );
        }
    }
    else
    {
        activeFrames.reserve( static_cast<std::size_t>( this->m_Nims ) );
        for( int image = 0; image < this->m_Nims; ++image )
        {
            activeFrames.push_back( image );
        }
    }
    if( !targetCubeFinite() )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidarg,
                                       "P4 target data must remain finite after preprocessing" );
    }
    if( this->m_mask.size() != 0 && ( this->m_mask.rows() != this->m_Nrows || this->m_mask.cols() != this->m_Ncols ) )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "P4 common-mask dimensions must match target images" );
    }

    if( this->m_preProcess_only && !this->m_skipPreProcess )
    {
        std::cerr << "Pre-processing complete, stopping.\n";
        this->t_end = mx::sys::get_curr_time();
        return 0;
    }

    const bool calculatePSF = !m_psfFile.empty();
    const bool processPSF = m_outputPSFModels || m_psfFilter;
    const bool targetHeldOutPSF = calculatePSF && m_excludeMethod != HCI::exclude::none;
    const bool sharedPSF = calculatePSF && !targetHeldOutPSF;
    std::optional<P4PSFModel> psfModel;
    if( calculatePSF )
    {
        try
        {
            imageT psfTemplate;
            mx::fits::fitsFile<realT, verboseT> reader;
            const mx::error_t readResult = reader.read( psfTemplate, m_psfFile );
            if( readResult != mx::error_t::noerror )
            {
                throw mx::exception<verboseT>( readResult, "could not read p4.psfFile " + m_psfFile );
            }
            m_psfTemplateRows = psfTemplate.rows();
            m_psfTemplateColumns = psfTemplate.cols();
            m_localPSFRows = localPSFModelDimension( m_psfStampSize, psfTemplate.rows() );
            m_localPSFColumns = localPSFModelDimension( m_psfStampSize, psfTemplate.cols() );
            psfModel.emplace( psfTemplate, m_localPSFRows, m_localPSFColumns );
            m_psfModelBytes = psfModel->storageBytes();
        }
        catch( ... )
        {
            std::throw_with_nested(
                mx::exception<verboseT>( mx::error_t::invalidconfig, "could not initialize the P4 PSF template" ) );
        }
    }

    const double geometryBegin = omp_get_wtime();
    std::vector<pixelGridT> grids( m_minRadius.size() );
    std::vector<rotatedGridT> rotatedGrids( m_minRadius.size() );
    m_realizedModes.resize( m_minRadius.size() );
    m_regionStatistics.resize( m_minRadius.size() );
    m_temporalSelections.resize( m_minRadius.size() );
    m_localPSFComponentCounts.resize( m_minRadius.size(), 1 );
    m_ownership.resize( this->m_Nrows, this->m_Ncols );
    m_ownership.setConstant( -1 );

    const imageT *mask = this->m_mask.size() == 0 ? nullptr : &this->m_mask;
    const std::vector<P4PixelCoordinate> temporalPredictorOffsets = p4TemporalPredictorOffsets( m_psfRadius );
    std::vector<double> derotationAngles;
    if( m_regressionFrame == P4RegressionFrame::rotated || m_numberImages > 0 || m_localStampSize > 0 ||
        m_excludeMethod == HCI::exclude::pixel || m_excludeMethod == HCI::exclude::angle )
    {
        derotationAngles.reserve( static_cast<std::size_t>( this->m_Nims ) );
        for( int image = 0; image < this->m_Nims; ++image )
        {
            const double angle = static_cast<double>( this->m_derotF.derotAngle( static_cast<std::size_t>( image ) ) );
            if( !mx::math::isFinite( angle ) )
            {
                throw mx::exception<verboseT>( mx::error_t::invalidarg, "P4 derotation angles must remain finite" );
            }
            derotationAngles.push_back( angle );
        }
    }

    for( std::size_t region = 0; region < m_minRadius.size(); ++region )
    {
        std::cerr << "P4 ";
        if( m_regressionFrame == P4RegressionFrame::rotated )
        {
            std::cerr << "rotated ";
        }
        std::cerr << "geometry annulus " << region + 1 << " / " << m_minRadius.size() << " [" << m_minRadius[region]
                  << ',' << m_maxRadius[region] << ")... " << std::flush;
        const P4PixelGridRegion configuration( m_minRadius[region],
                                               m_maxRadius[region],
                                               m_orDeltaRadiusInner,
                                               m_orDeltaRadiusOuter,
                                               m_orArcHalfWidth,
                                               m_orMaxHalfAngle,
                                               m_psfRadius,
                                               *m_exclusionPolicy,
                                               m_exclusionRadiusBuffer );
        try
        {
            grids[region].resize( this->m_Nrows, this->m_Ncols );
            if( m_regressionFrame == P4RegressionFrame::detector )
            {
                grids[region].region( configuration, mask );
            }
            else
            {
                grids[region].candidateRegion( configuration );
                rotatedGrids[region].configure( grids[region], derotationAngles, mask );
                grids[region] = pixelGridT();
            }
        }
        catch( ... )
        {
            std::cerr << "failed\n";
            std::throw_with_nested(
                mx::exception<verboseT>( mx::error_t::invalidconfig,
                                         "could not construct P4 region " + std::to_string( region ) ) );
        }

        const std::size_t searchPixelCount = m_regressionFrame == P4RegressionFrame::detector
                                                 ? grids[region].searchPixelCount()
                                                 : rotatedGrids[region].searchPixelCount();
        const std::size_t basePredictorCount = m_regressionFrame == P4RegressionFrame::detector
                                                   ? grids[region].predictorCount()
                                                   : rotatedGrids[region].predictorCount();
        double temporalPsfRadius{ 0 };
        if( m_regressionFrame == P4RegressionFrame::detector )
        {
            if( m_numberImages == 0 )
            {
                m_temporalSelections[region].reserve( activeFrames.size() );
                for( const int image : activeFrames )
                {
                    m_temporalSelections[region].push_back( { image } );
                }
            }
            else
            {
                const double meanRadius =
                    0.5 * ( static_cast<double>( m_minRadius[region] ) + static_cast<double>( m_maxRadius[region] ) );
                auto selection = p4TemporalSelectionWithFallback( derotationAngles,
                                                                  activeFrames,
                                                                  meanRadius,
                                                                  2 * static_cast<double>( m_psfRadius ),
                                                                  m_numberImages,
                                                                  basePredictorCount,
                                                                  temporalPredictorOffsets.size(),
                                                                  m_modeFractions );
                temporalPsfRadius = selection.first;
                m_temporalSelections[region] = std::move( selection.second );
                if( m_temporalSelections[region].empty() )
                {
                    m_temporalSelections[region].reserve( activeFrames.size() );
                    for( const int image : activeFrames )
                    {
                        m_temporalSelections[region].push_back( { image } );
                    }
                }
            }
        }
        else
        {
            m_temporalSelections[region].reserve( activeFrames.size() );
            for( const int image : activeFrames )
            {
                m_temporalSelections[region].push_back( { image } );
            }
        }
        const std::size_t targetImageCount = m_temporalSelections[region].size();
        const std::size_t temporalImageCount = m_temporalSelections[region].front().size() - 1;
        m_localPSFComponentCounts[region] = temporalImageCount + 1;
        if( temporalImageCount > std::numeric_limits<std::size_t>::max() / temporalPredictorOffsets.size() )
        {
            throw mx::exception<verboseT>( mx::error_t::sizeerr,
                                           "P4 multi-image predictor count exceeds size_t range" );
        }
        const std::size_t temporalPredictorCount = temporalImageCount * temporalPredictorOffsets.size();
        if( basePredictorCount > std::numeric_limits<std::size_t>::max() - temporalPredictorCount )
        {
            throw mx::exception<verboseT>( mx::error_t::sizeerr,
                                           "P4 multi-image predictor count exceeds size_t range" );
        }
        const std::size_t predictorCount = basePredictorCount + temporalPredictorCount;
        std::cerr << "done: " << searchPixelCount << " search pixels, K=" << predictorCount << '\n';

        const int maximumDegreesOfFreedom =
            checkedMaximumDegreesOfFreedom( static_cast<int>( targetImageCount ),
                                            predictorCount,
                                            m_regressionFrame == P4RegressionFrame::rotated );
        int previousMode{ 0 };
        for( const realT fraction : m_modeFractions )
        {
            const int mode = static_cast<int>(
                std::floor( static_cast<double>( fraction ) * static_cast<double>( maximumDegreesOfFreedom ) ) );
            if( mode <= 0 )
            {
                throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                               "P4 mode fraction realizes to zero modes in region " +
                                                   std::to_string( region ) );
            }
            if( mode <= previousMode )
            {
                throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                               "P4 mode fractions realize to duplicate modes in region " +
                                                   std::to_string( region ) );
            }
            m_realizedModes[region].push_back( mode );
            previousMode = mode;
        }

        for( std::size_t search = 0; search < searchPixelCount; ++search )
        {
            const P4PixelCoordinate &coordinate = m_regressionFrame == P4RegressionFrame::detector
                                                      ? grids[region].searchPixel( search ).coordinate()
                                                      : rotatedGrids[region].searchPixel( search ).coordinate();
            claimOwnership( m_ownership, coordinate, static_cast<int>( region ) );
        }

        P4RegionStatistics &statistics = m_regionStatistics[region];
        statistics.searchPixelCount = searchPixelCount;
        statistics.targetImageCount = targetImageCount;
        statistics.temporalNumberImages = static_cast<int>( temporalImageCount / 2 );
        statistics.temporalPsfRadius = temporalPsfRadius;
        statistics.predictorCount = predictorCount;
        statistics.maximumDegreesOfFreedom = maximumDegreesOfFreedom;
        statistics.rankInvalidCounts.assign( m_modeFractions.size(), 0 );
    }
    m_timing.geometryElapsedSeconds = omp_get_wtime() - geometryBegin;

    std::vector<P4TargetExclusions> regionExclusions( m_regionStatistics.size() );
    if( m_excludeMethod != HCI::exclude::none )
    {
        for( std::size_t region = 0; region < m_regionStatistics.size(); ++region )
        {
            regionExclusions[region] = targetExclusions( m_temporalSelections[region],
                                                         m_excludeMethod,
                                                         m_minDPx,
                                                         m_minRadius[region],
                                                         derotationAngles );
            P4RegionStatistics &statistics = m_regionStatistics[region];
            statistics.maximumExcludedRows = static_cast<std::size_t>( regionExclusions[region].maximumDeleted() );
            statistics.exclusionStorageBytes = regionExclusions[region].storageBytes();
            if( m_exclusionSolver == P4ExclusionSolver::factorDowndateExact &&
                m_deletionBackend == mx::math::svdDeletionBackend::rankOneSecular )
            {
                const Eigen::Index targetCount = regionExclusions[region].targetCount();
                for( Eigen::Index target = 0; target < targetCount; ++target )
                {
                    const Eigen::Index deletedCount = regionExclusions[region].deletedCount( target );
                    if( deletedCount < targetCount && deletedCount != 1 )
                    {
                        throw mx::exception<verboseT>(
                            mx::error_t::invalidconfig,
                            "solver.deletionBackend=rankOneSecular requires exactly one excluded row for every "
                            "retained target; annulus " +
                                std::to_string( region ) + ", target row " + std::to_string( target ) + " excludes " +
                                std::to_string( deletedCount ) );
                    }
                }
            }
            if( m_targetExclusionBytes > std::numeric_limits<std::size_t>::max() - statistics.exclusionStorageBytes )
            {
                throw mx::exception<verboseT>( mx::error_t::sizeerr,
                                               "P4 compact target-exclusion byte count overflow" );
            }
            m_targetExclusionBytes += statistics.exclusionStorageBytes;
        }
    }

    if( m_localStampSize > 0 )
    {
        return processLocalTrial( grids, temporalPredictorOffsets, derotationAngles, regionExclusions );
    }

    std::size_t totalSearchPixels{ 0 };
    for( const P4RegionStatistics &statistics : m_regionStatistics )
    {
        totalSearchPixels += statistics.searchPixelCount;
    }

    const bool compactFinalization = this->m_combineMethod != HCI::combine::none && !this->m_doOutputPSFSub;
    using compactResidualT = Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic>;
    using compactValidityT = Eigen::Array<std::uint8_t, Eigen::Dynamic, Eigen::Dynamic>;
    std::vector<compactResidualT> compactResiduals;
    std::vector<compactValidityT> compactValidity;

    std::size_t psfStampPixels{ 0 };
    std::size_t heldOutPSFModeBytes{ 0 };
    if( calculatePSF )
    {
        const long double localStampPixels =
            static_cast<long double>( m_localPSFRows ) * static_cast<long double>( m_localPSFColumns );
        if( localStampPixels > static_cast<long double>( std::numeric_limits<std::size_t>::max() ) )
        {
            throw mx::exception<verboseT>( mx::error_t::sizeerr, "P4 local PSF dimensions exceed size_t range" );
        }
        psfStampPixels = static_cast<std::size_t>( localStampPixels );
        for( std::size_t region = 0; region < m_regionStatistics.size(); ++region )
        {
            const P4RegionStatistics &statistics = m_regionStatistics[region];
            std::size_t retainedSearchCount = statistics.searchPixelCount;
            std::size_t retainedModeCount = m_modeFractions.size();
            std::size_t temporalCoefficientCount =
                ( m_localPSFComponentCounts[region] - 1 ) * temporalPredictorOffsets.size();
            if( targetHeldOutPSF )
            {
                if( statistics.targetImageCount != static_cast<std::size_t>( this->m_Nims ) ||
                    retainedSearchCount > std::numeric_limits<std::size_t>::max() / statistics.targetImageCount )
                {
                    throw mx::exception<verboseT>(
                        mx::error_t::sizeerr,
                        "P4 target-held-out PSF dimensions require one target row per loaded image" );
                }
                retainedSearchCount *= statistics.targetImageCount;
                retainedModeCount = 1;
                temporalCoefficientCount = 0;
            }
            const std::size_t regionBytes = localPSFBytes( retainedSearchCount,
                                                           retainedModeCount,
                                                           temporalCoefficientCount,
                                                           m_localPSFRows,
                                                           m_localPSFColumns );
            std::size_t &totalBytes = targetHeldOutPSF ? heldOutPSFModeBytes : m_localPSFBytes;
            if( totalBytes > std::numeric_limits<std::size_t>::max() - regionBytes )
            {
                throw mx::exception<verboseT>( mx::error_t::sizeerr, "P4 local PSF byte count overflow" );
            }
            totalBytes += regionBytes;
        }
        if( targetHeldOutPSF )
        {
            m_localPSFBytes = heldOutPSFModeBytes;
            m_psfModeBatchSize = 1;
        }
        if( processPSF )
        {
            m_psfReconstructionBytes = psfReconstructionBytes(
                totalSearchPixels,
                static_cast<std::size_t>( this->m_Nims ),
                ( *std::max_element( m_localPSFComponentCounts.begin(), m_localPSFComponentCounts.end() ) - 1 ) *
                    temporalPredictorOffsets.size(),
                m_psfStampSize,
                m_localPSFRows,
                m_localPSFColumns );
            if( m_psfFilter )
            {
                m_psfFilterBytes = psfFilterBytes( this->m_Nrows, this->m_Ncols, m_modeFractions.size() );
                if( m_psfReconstructionBytes > std::numeric_limits<std::size_t>::max() - m_psfFilterBytes )
                {
                    throw mx::exception<verboseT>( mx::error_t::sizeerr, "P4 PSF output byte count overflow" );
                }
                m_psfReconstructionBytes += m_psfFilterBytes;
            }
        }
    }

    const auto configureMemoryBudget = [&]()
    {
        if( m_memoryFraction == 0 )
        {
            return;
        }
        const std::optional<std::size_t> availableMemory = p4AvailableMemoryBytes();
        if( !availableMemory )
        {
            std::cerr << "WARNING: P4 could not read Linux MemAvailable; automatic worker limiting is disabled\n";
            return;
        }
        m_availableMemoryBytes = *availableMemory;
        m_memoryBudgetBytes = static_cast<std::size_t>(
            std::floor( static_cast<long double>( m_memoryFraction ) * m_availableMemoryBytes ) );
        if( m_memoryBudgetBytes == 0 )
        {
            throw mx::exception<verboseT>( mx::error_t::allocerr,
                                           "p4.memoryFraction selects a zero-byte automatic memory budget" );
        }
    };

    if( compactFinalization )
    {
        for( std::size_t region = 0; region < m_regionStatistics.size(); ++region )
        {
            const P4RegionStatistics &statistics = m_regionStatistics[region];
            if( statistics.searchPixelCount > std::numeric_limits<std::size_t>::max() / m_modeFractions.size() )
            {
                throw mx::exception<verboseT>( mx::error_t::sizeerr, "P4 compact residual dimensions overflow" );
            }
            const std::size_t compactColumnCount = statistics.searchPixelCount * m_modeFractions.size();
            if( compactColumnCount > static_cast<std::size_t>( std::numeric_limits<Eigen::Index>::max() ) ||
                statistics.targetImageCount > static_cast<std::size_t>( std::numeric_limits<Eigen::Index>::max() ) ||
                ( compactColumnCount != 0 && statistics.targetImageCount > std::numeric_limits<std::size_t>::max() /
                                                                               compactColumnCount / sizeof( realT ) ) )
            {
                throw mx::exception<verboseT>( mx::error_t::sizeerr, "P4 compact residual dimensions overflow" );
            }
            const std::size_t residualBytes = statistics.targetImageCount * compactColumnCount * sizeof( realT );
            const std::size_t validityBytes = statistics.targetImageCount * compactColumnCount * sizeof( std::uint8_t );
            if( residualBytes > std::numeric_limits<std::size_t>::max() - validityBytes ||
                m_compactResidualBytes > std::numeric_limits<std::size_t>::max() - residualBytes - validityBytes )
            {
                throw mx::exception<verboseT>( mx::error_t::sizeerr, "P4 compact residual byte count overflow" );
            }
            m_compactResidualBytes += residualBytes + validityBytes;
        }

        const long double materializationBytes = 2.0L * static_cast<long double>( this->m_Nrows ) *
                                                 static_cast<long double>( this->m_Ncols ) *
                                                 static_cast<long double>( this->m_Nims ) * sizeof( realT );
        if( materializationBytes > static_cast<long double>( std::numeric_limits<std::size_t>::max() ) )
        {
            throw mx::exception<verboseT>( mx::error_t::sizeerr, "P4 materialization byte count overflow" );
        }
        m_materializationBytes = static_cast<std::size_t>( materializationBytes );

        configureMemoryBudget();
        if( m_compactResidualBytes > std::numeric_limits<std::size_t>::max() - m_localPSFBytes ||
            m_compactResidualBytes + m_localPSFBytes > std::numeric_limits<std::size_t>::max() - m_psfModelBytes ||
            m_compactResidualBytes + m_localPSFBytes + m_psfModelBytes >
                std::numeric_limits<std::size_t>::max() - m_targetExclusionBytes )
        {
            throw mx::exception<verboseT>( mx::error_t::sizeerr, "P4 persistent byte count overflow" );
        }
        const std::size_t persistentBytes =
            m_compactResidualBytes + m_localPSFBytes + m_psfModelBytes + m_targetExclusionBytes;
        const std::size_t postRegressionScratch = std::max( m_materializationBytes, m_psfReconstructionBytes );
        if( m_memoryBudgetBytes != 0 &&
            ( persistentBytes > m_memoryBudgetBytes || postRegressionScratch > m_memoryBudgetBytes - persistentBytes ) )
        {
            throw mx::exception<verboseT>(
                mx::error_t::allocerr,
                "P4 compact residuals, target exclusions, local PSFs, and post-regression output scratch exceed the "
                "automatic memory "
                "budget; "
                "increase p4.memoryFraction, reduce the frame count, or set p4.memoryFraction=0 to disable the "
                "automatic budget" );
        }

        compactResiduals.resize( m_regionStatistics.size() );
        compactValidity.resize( m_regionStatistics.size() );
        for( std::size_t region = 0; region < m_regionStatistics.size(); ++region )
        {
            const P4RegionStatistics &statistics = m_regionStatistics[region];
            const std::size_t compactColumnCount = statistics.searchPixelCount * m_modeFractions.size();
            compactResiduals[region].resize( static_cast<Eigen::Index>( statistics.targetImageCount ),
                                             static_cast<Eigen::Index>( compactColumnCount ) );
            compactValidity[region].resize( static_cast<Eigen::Index>( statistics.targetImageCount ),
                                            static_cast<Eigen::Index>( compactColumnCount ) );
            compactValidity[region].setZero();
        }
    }
    else
    {
        this->m_psfsub.resize( m_modeFractions.size() );
        this->m_psfsubValidity.resize( m_modeFractions.size() );
        for( std::size_t output = 0; output < m_modeFractions.size(); ++output )
        {
            this->m_psfsub[output].resize( this->m_Nrows, this->m_Ncols, this->m_Nims );
            this->m_psfsub[output].setZero();
            this->m_psfsubValidity[output].resize( this->m_Nrows, this->m_Ncols, this->m_Nims );
            this->m_psfsubValidity[output].setZero();
        }
        configureMemoryBudget();
        if( m_localPSFBytes > std::numeric_limits<std::size_t>::max() - m_psfModelBytes ||
            m_localPSFBytes + m_psfModelBytes > std::numeric_limits<std::size_t>::max() - m_targetExclusionBytes )
        {
            throw mx::exception<verboseT>( mx::error_t::sizeerr, "P4 persistent byte count overflow" );
        }
        const std::size_t persistentBytes = m_localPSFBytes + m_psfModelBytes + m_targetExclusionBytes;
        if( m_memoryBudgetBytes != 0 && ( persistentBytes > m_memoryBudgetBytes ||
                                          m_psfReconstructionBytes > m_memoryBudgetBytes - persistentBytes ) )
        {
            throw mx::exception<verboseT>(
                mx::error_t::allocerr,
                "P4 target exclusions and local PSF models exceed the automatic memory budget; increase "
                "p4.memoryFraction, reduce the search or stamp size, or set p4.memoryFraction=0 to disable the "
                "automatic budget" );
        }
    }

    if( targetHeldOutPSF )
    {
        const std::size_t outputCount = m_modeFractions.size();
        std::size_t selectedBatch = outputCount;
        if( m_memoryBudgetBytes != 0 )
        {
            const long double fixedBytes = static_cast<long double>( m_compactResidualBytes ) +
                                           static_cast<long double>( m_targetExclusionBytes ) + m_psfModelBytes +
                                           m_psfReconstructionBytes;
            long double bestWork = std::numeric_limits<long double>::infinity();
            selectedBatch = 0;
            for( std::size_t batch = 1; batch <= outputCount; ++batch )
            {
                const long double batchBytes = static_cast<long double>( heldOutPSFModeBytes ) * batch;
                if( fixedBytes + batchBytes >= static_cast<long double>( m_memoryBudgetBytes ) )
                {
                    continue;
                }
                const std::size_t persistentBytes = static_cast<std::size_t>( fixedBytes + batchBytes );
                long double regionWork{ 0 };
                bool supported{ true };
                for( std::size_t region = 0; region < m_regionStatistics.size(); ++region )
                {
                    const P4RegionStatistics &statistics = m_regionStatistics[region];
                    const std::size_t maximumDeletedRows =
                        m_exclusionSolver == P4ExclusionSolver::factorDowndateExact &&
                                m_deletionBackend == mx::math::svdDeletionBackend::rankOneSecular
                            ? 1
                            : statistics.maximumExcludedRows;
                    const std::size_t workerBytes =
                        estimatedWorkerBytes( statistics.targetImageCount,
                                              statistics.predictorCount,
                                              batch,
                                              false,
                                              0,
                                              true,
                                              m_exclusionSolver,
                                              maximumDeletedRows,
                                              static_cast<std::size_t>( m_realizedModes[region].back() ),
                                              psfStampPixels );
                    const int requestedWorkers =
                        std::max( 1,
                                  std::min( omp_get_max_threads(),
                                            static_cast<int>( std::min<std::size_t>(
                                                statistics.searchPixelCount,
                                                static_cast<std::size_t>( std::numeric_limits<int>::max() ) ) ) ) );
                    const int workers =
                        memoryLimitedWorkerCount( requestedWorkers, m_memoryBudgetBytes, persistentBytes, workerBytes );
                    if( workers == 0 )
                    {
                        supported = false;
                        break;
                    }
                    regionWork += static_cast<long double>( statistics.searchPixelCount ) / workers;
                }
                if( !supported )
                {
                    continue;
                }
                const std::size_t passes = ( outputCount + batch - 1 ) / batch;
                const long double work = static_cast<long double>( passes ) * regionWork;
                if( work < bestWork || ( work == bestWork && batch > selectedBatch ) )
                {
                    bestWork = work;
                    selectedBatch = batch;
                }
            }
            if( selectedBatch == 0 )
            {
                throw mx::exception<verboseT>(
                    mx::error_t::allocerr,
                    "one target-held-out P4 PSF response mode and worker do not fit the automatic memory budget; "
                    "increase p4.memoryFraction, reduce the frame/search/stamp size, or set p4.memoryFraction=0 to "
                    "disable the automatic budget" );
            }
        }
        if( heldOutPSFModeBytes > std::numeric_limits<std::size_t>::max() / selectedBatch )
        {
            throw mx::exception<verboseT>( mx::error_t::sizeerr, "P4 target-held-out PSF batch byte count overflow" );
        }
        m_psfModeBatchSize = selectedBatch;
        m_localPSFBytes = heldOutPSFModeBytes * selectedBatch;
    }

    if( sharedPSF )
    {
        m_localPSFModels.resize( m_regionStatistics.size() );
        m_localPSFTemporalCoefficients.resize( m_regionStatistics.size() );
        m_localPSFValidity.resize( m_regionStatistics.size() );
        for( std::size_t region = 0; region < m_regionStatistics.size(); ++region )
        {
            const P4RegionStatistics &statistics = m_regionStatistics[region];
            if( statistics.searchPixelCount > std::numeric_limits<std::size_t>::max() / m_modeFractions.size() )
            {
                throw mx::exception<verboseT>( mx::error_t::sizeerr, "P4 local PSF dimensions overflow" );
            }
            const std::size_t compactColumnCount = statistics.searchPixelCount * m_modeFractions.size();
            if( psfStampPixels > static_cast<std::size_t>( std::numeric_limits<Eigen::Index>::max() ) ||
                compactColumnCount > static_cast<std::size_t>( std::numeric_limits<Eigen::Index>::max() ) )
            {
                throw mx::exception<verboseT>( mx::error_t::sizeerr, "P4 local PSF dimensions exceed Eigen range" );
            }
            m_localPSFModels[region].resize( static_cast<Eigen::Index>( psfStampPixels ),
                                             static_cast<Eigen::Index>( compactColumnCount ) );
            m_localPSFModels[region].setZero();
            const std::size_t temporalCoefficientCount =
                ( m_localPSFComponentCounts[region] - 1 ) * temporalPredictorOffsets.size();
            if( temporalCoefficientCount > static_cast<std::size_t>( std::numeric_limits<Eigen::Index>::max() ) )
            {
                throw mx::exception<verboseT>( mx::error_t::sizeerr,
                                               "P4 temporal PSF coefficient dimensions exceed Eigen range" );
            }
            m_localPSFTemporalCoefficients[region].resize( static_cast<Eigen::Index>( temporalCoefficientCount ),
                                                           static_cast<Eigen::Index>( compactColumnCount ) );
            m_localPSFTemporalCoefficients[region].setZero();
            m_localPSFValidity[region].resize( static_cast<Eigen::Index>( statistics.searchPixelCount ),
                                               static_cast<Eigen::Index>( m_modeFractions.size() ) );
            m_localPSFValidity[region].setZero();
        }
    }

    constexpr double bytesPerGiB = 1024.0 * 1024.0 * 1024.0;
    if( m_memoryBudgetBytes != 0 )
    {
        std::cerr << "P4 memory policy: " << static_cast<double>( m_availableMemoryBytes ) / bytesPerGiB
                  << " GiB available, " << static_cast<double>( m_memoryBudgetBytes ) / bytesPerGiB << " GiB budget, "
                  << static_cast<double>( m_compactResidualBytes ) / bytesPerGiB << " GiB compact residuals, "
                  << static_cast<double>( m_targetExclusionBytes ) / bytesPerGiB << " GiB target exclusions, "
                  << static_cast<double>( m_materializationBytes ) / bytesPerGiB << " GiB output materialization";
        if( calculatePSF )
        {
            std::cerr << ", " << static_cast<double>( m_localPSFBytes ) / bytesPerGiB
                      << ( targetHeldOutPSF ? " GiB held-out PSF batch, " : " GiB local PSFs, " )
                      << static_cast<double>( m_psfModelBytes ) / bytesPerGiB << " GiB PSF template model";
            if( targetHeldOutPSF )
            {
                std::cerr << ", " << m_psfModeBatchSize << " PSF modes per batch";
            }
            if( processPSF )
            {
                std::cerr << ", " << static_cast<double>( m_psfReconstructionBytes ) / bytesPerGiB
                          << " GiB PSF reconstruction scratch";
            }
        }
        std::cerr << '\n';
    }
    else if( m_memoryFraction == 0 )
    {
        std::cerr << "P4 memory policy: automatic worker limiting disabled by p4.memoryFraction=0\n";
    }

    for( std::size_t region = 0; region < m_regionStatistics.size(); ++region )
    {
        P4RegionStatistics &statistics = m_regionStatistics[region];
        const std::size_t estimatedMaximumExcludedRows =
            m_exclusionSolver == P4ExclusionSolver::factorDowndateExact &&
                    m_deletionBackend == mx::math::svdDeletionBackend::rankOneSecular
                ? 1
                : statistics.maximumExcludedRows;
        statistics.estimatedWorkerBytes =
            estimatedWorkerBytes( statistics.targetImageCount,
                                  statistics.predictorCount,
                                  m_modeFractions.size(),
                                  sharedPSF,
                                  sharedPSF ? psfStampPixels : 0,
                                  m_excludeMethod != HCI::exclude::none,
                                  m_exclusionSolver,
                                  estimatedMaximumExcludedRows,
                                  static_cast<std::size_t>( m_realizedModes[region].back() ) );
        statistics.maximumWorkerCount =
            std::max( 1,
                      std::min( omp_get_max_threads(),
                                static_cast<int>( std::min<std::size_t>(
                                    statistics.searchPixelCount,
                                    static_cast<std::size_t>( std::numeric_limits<int>::max() ) ) ) ) );
        statistics.effectiveWorkerCount = statistics.maximumWorkerCount;
        if( m_memoryBudgetBytes != 0 )
        {
            statistics.effectiveWorkerCount = memoryLimitedWorkerCount(
                statistics.maximumWorkerCount,
                m_memoryBudgetBytes,
                m_compactResidualBytes + m_targetExclusionBytes + ( sharedPSF ? m_localPSFBytes : 0 ) + m_psfModelBytes,
                statistics.estimatedWorkerBytes );
            if( statistics.effectiveWorkerCount == 0 )
            {
                throw mx::exception<verboseT>(
                    mx::error_t::allocerr,
                    "one P4 worker does not fit the automatic memory budget; increase p4.memoryFraction, reduce the "
                    "frame or predictor count, or set p4.memoryFraction=0 to disable the automatic budget" );
            }
        }
        std::cerr << "P4 memory annulus " << region + 1 << " / " << m_regionStatistics.size() << ": worker estimate "
                  << static_cast<double>( statistics.estimatedWorkerBytes ) / ( 1024.0 * 1024.0 ) << " MiB, workers "
                  << statistics.effectiveWorkerCount << " / " << statistics.maximumWorkerCount << '\n';
    }

    const double regressionBegin = omp_get_wtime();
    P4ProgressOutput progressOutput( std::cerr,
                                     grids.size(),
                                     totalSearchPixels,
                                     m_regressionFrame == P4RegressionFrame::rotated );
    std::size_t completedSearchPixels{ 0 };
    for( std::size_t region = 0; region < grids.size(); ++region )
    {
        const pixelGridT &grid = grids[region];
        const rotatedGridT &rotatedGrid = rotatedGrids[region];
        const bool rotated = m_regressionFrame == P4RegressionFrame::rotated;
        const std::size_t searchPixelCount = m_regionStatistics[region].searchPixelCount;
        const std::vector<std::vector<int>> &temporalSelections = m_temporalSelections[region];
        const std::size_t targetImageCount = temporalSelections.size();
        const P4TargetExclusions *exclusions = regionExclusions[region].empty() ? nullptr : &regionExclusions[region];
        const bool usesTemporalPredictors = m_regionStatistics[region].temporalNumberImages > 0;
        const std::size_t predictorCount = m_regionStatistics[region].predictorCount;
        const std::vector<int> &modes = m_realizedModes[region];
        const int effectiveWorkerCount = m_regionStatistics[region].effectiveWorkerCount;
        progressOutput.beginRegion( region,
                                    searchPixelCount,
                                    predictorCount,
                                    completedSearchPixels,
                                    static_cast<double>( m_minRadius[region] ),
                                    static_cast<double>( m_maxRadius[region] ) );
        mx::ipc::ompLoopWatcher<P4ProgressOutput, true, true, true, true, true> progressWatcher( searchPixelCount,
                                                                                                 progressOutput );
        std::exception_ptr workerException;
        std::atomic<bool> failed{ false };
        std::size_t validLocalFitCount{ 0 };
        std::size_t maskedLocalFitCount{ 0 };
        std::size_t supportInvalidLocalFitCount{ 0 };
        int minimumNumericalRank{ std::numeric_limits<int>::max() };
        int minimumBaseRank{ std::numeric_limits<int>::max() };
        std::size_t downdateClampCount{ 0 };
        std::size_t explicitFallbackCount{ 0 };
        std::size_t rankBoundaryFallbackPixelCount{ 0 };
        std::size_t factorValidationFallbackPixelCount{ 0 };
        std::size_t deletionSolverFallbackPixelCount{ 0 };
        double maximumFactorOrthogonalityDefect{ 0 };
        double factorOrthogonalityToleranceAtMaximumDefect{ 0 };
        std::vector<std::size_t> rankInvalidCounts( modes.size(), 0 );

        // clang-format off
#pragma omp parallel num_threads(effectiveWorkerCount)
        // clang-format on
        {
            P4PCA::workspaceT workspace;
            detail::P4PCAMixedWorkspace mixedWorkspace;
#ifdef HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION
            std::optional<detail::P4PCAExperimentalWorkspace> experimentalWorkspace;
            std::optional<P4ReductionWorkerPrecisionScope> experimentalPrecisionScope;
            if( experimentalDispatch )
            {
                experimentalWorkspace.emplace();
                experimentalPrecisionScope.emplace( *experimentalDispatch, *experimentalWorkspace );
            }
#endif
            P4PCADowndateWorkspace downdateWorkspace;
            P4PCA::matrixT predictors( static_cast<Eigen::Index>( targetImageCount ),
                                       static_cast<Eigen::Index>( predictorCount ) );
            P4PCA::vectorT target( static_cast<Eigen::Index>( targetImageCount ) );
            P4PCAResult result;
            P4PCA::matrixT coefficients;
            imageT localPSFResponse;
            std::size_t threadValid{ 0 };
            std::size_t threadMasked{ 0 };
            std::size_t threadSupportInvalid{ 0 };
            int threadMinimumRank{ std::numeric_limits<int>::max() };
            int threadMinimumBaseRank{ std::numeric_limits<int>::max() };
            std::size_t threadDowndateClampCount{ 0 };
            std::size_t threadExplicitFallbackCount{ 0 };
            std::size_t threadRankBoundaryFallbackPixelCount{ 0 };
            std::size_t threadFactorValidationFallbackPixelCount{ 0 };
            std::size_t threadDeletionSolverFallbackPixelCount{ 0 };
            double threadMaximumFactorOrthogonalityDefect{ 0 };
            double threadFactorOrthogonalityToleranceAtMaximumDefect{ 0 };
            std::vector<std::size_t> threadRankInvalid( modes.size(), 0 );
            double threadSamplingSeconds{ 0 };
            double threadSameImageSamplingSeconds{ 0 };
            double threadTemporalSamplingSeconds{ 0 };
            double threadGramSeconds{ 0 };
            double threadEigensolveSeconds{ 0 };
            double threadBaseFactorSeconds{ 0 };
            double threadDeletionSeconds{ 0 };
            double threadExplicitFallbackSeconds{ 0 };
            double threadProjectionSeconds{ 0 };
            double threadPSFSeconds{ 0 };

            // clang-format off
#pragma omp for schedule(static)
            // clang-format on
            for( std::size_t search = 0; search < searchPixelCount; ++search )
            {
                const bool attempt = !failed.load( std::memory_order_acquire );
                if( attempt )
                {
                    try
                    {
                        const bool searchValid =
                            rotated ? rotatedGrid.searchPixel( search ).valid() : grid.searchPixel( search ).valid();
                        const P4PixelCoordinate &coordinate = rotated ? rotatedGrid.searchPixel( search ).coordinate()
                                                                      : grid.searchPixel( search ).coordinate();
                        const bool temporalPredictorsValid =
                            !usesTemporalPredictors || p4TemporalPredictorsValid( coordinate,
                                                                                  temporalPredictorOffsets,
                                                                                  this->m_Nrows,
                                                                                  this->m_Ncols,
                                                                                  mask );
                        if( !searchValid || !temporalPredictorsValid )
                        {
                            if( rotated )
                            {
                                ++threadSupportInvalid;
                            }
                            else
                            {
                                ++threadMasked;
                            }
                        }
                        else
                        {
                            double sameImageSamplingSeconds{ 0 };
                            double temporalSamplingSeconds{ 0 };
                            P4PCATiming pcaTiming;
                            const int row = coordinate.row();
                            const int column = coordinate.column();
                            if( rotated )
                            {
                                const double sameImageSamplingBegin = omp_get_wtime();
                                for( std::size_t targetIndex = 0; targetIndex < targetImageCount; ++targetIndex )
                                {
                                    const int image = temporalSelections[targetIndex][0];
                                    target( static_cast<Eigen::Index>( targetIndex ) ) = checkedPredictorPromotion(
                                        rotatedGrid.sampleTarget( this->m_tgtIms.image( image ),
                                                                  static_cast<std::size_t>( image ),
                                                                  search ) );
                                    for( std::size_t predictor = 0; predictor < predictorCount; ++predictor )
                                    {
                                        predictors( static_cast<Eigen::Index>( targetIndex ), predictor ) =
                                            checkedPredictorPromotion(
                                                rotatedGrid.samplePredictor( this->m_tgtIms.image( image ),
                                                                             static_cast<std::size_t>( image ),
                                                                             search,
                                                                             predictor ) );
                                    }
                                }
                                sameImageSamplingSeconds += omp_get_wtime() - sameImageSamplingBegin;
                            }
                            else
                            {
                                fitDetectorSearch( result,
                                                   predictors,
                                                   target,
                                                   sharedPSF ? &coefficients : nullptr,
                                                   nullptr,
                                                   grid,
                                                   search,
                                                   temporalSelections,
                                                   temporalPredictorOffsets,
                                                   modes,
                                                   workspace,
                                                   mixedWorkspace,
                                                   &downdateWorkspace,
                                                   pcaTiming,
                                                   sameImageSamplingSeconds,
                                                   temporalSamplingSeconds,
                                                   exclusions,
                                                   nullptr,
                                                   nullptr );
                            }

                            threadSamplingSeconds += sameImageSamplingSeconds + temporalSamplingSeconds;
                            threadSameImageSamplingSeconds += sameImageSamplingSeconds;
                            threadTemporalSamplingSeconds += temporalSamplingSeconds;

                            if( rotated )
                            {
#ifdef HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION
                                p4ReductionDispatchCenteredInPlace( result,
                                                                    predictors,
                                                                    target,
                                                                    modes,
                                                                    m_rankTolerance,
                                                                    workspace,
                                                                    mixedWorkspace,
                                                                    &pcaTiming,
                                                                    sharedPSF ? &coefficients : nullptr );
#else
                                static_cast<void>( workspace );
                                detail::p4PCACalculateCenteredInPlaceMixed( result,
                                                                            predictors,
                                                                            target,
                                                                            modes,
                                                                            m_rankTolerance,
                                                                            mixedWorkspace,
                                                                            &pcaTiming,
                                                                            sharedPSF ? &coefficients : nullptr );
#endif
                            }
                            threadGramSeconds += pcaTiming.gramWorkerSeconds;
                            threadEigensolveSeconds += pcaTiming.eigensolveWorkerSeconds;
                            threadBaseFactorSeconds += pcaTiming.baseFactorWorkerSeconds;
                            threadDeletionSeconds += pcaTiming.deletionWorkerSeconds;
                            threadExplicitFallbackSeconds += pcaTiming.explicitFallbackWorkerSeconds;
                            threadProjectionSeconds += pcaTiming.projectionWorkerSeconds;
                            ++threadValid;
                            threadMinimumRank = std::min( threadMinimumRank, result.numericalRank );
                            if( m_exclusionSolver == P4ExclusionSolver::factorDowndateExact )
                            {
                                threadMinimumBaseRank = std::min( threadMinimumBaseRank, result.baseRank );
                                threadDowndateClampCount += result.downdateClampCount;
                                threadExplicitFallbackCount += result.explicitFallbackCount;
                                switch( result.explicitFallbackReason )
                                {
                                case P4PCAFallbackReason::none:
                                    break;
                                case P4PCAFallbackReason::rankBoundary:
                                    ++threadRankBoundaryFallbackPixelCount;
                                    break;
                                case P4PCAFallbackReason::factorValidation:
                                    ++threadFactorValidationFallbackPixelCount;
                                    if( threadFactorValidationFallbackPixelCount == 1 ||
                                        result.factorOrthogonalityDefect > threadMaximumFactorOrthogonalityDefect )
                                    {
                                        threadMaximumFactorOrthogonalityDefect = result.factorOrthogonalityDefect;
                                        threadFactorOrthogonalityToleranceAtMaximumDefect =
                                            result.factorOrthogonalityTolerance;
                                    }
                                    break;
                                case P4PCAFallbackReason::deletionSolver:
                                    ++threadDeletionSolverFallbackPixelCount;
                                    break;
                                }
                            }

                            const double residualApplyBegin = omp_get_wtime();
                            for( std::size_t output = 0; output < modes.size(); ++output )
                            {
                                const Eigen::Index compactColumn =
                                    static_cast<Eigen::Index>( search * modes.size() + output );
                                bool allTargetsValid{ true };
                                for( std::size_t targetIndex = 0; targetIndex < targetImageCount; ++targetIndex )
                                {
                                    if( !result.sampleSupported( static_cast<Eigen::Index>( targetIndex ), output ) )
                                    {
                                        allTargetsValid = false;
                                        continue;
                                    }
                                    const int image = temporalSelections[targetIndex][0];
                                    const realT residual = checkedResidualCast(
                                        result.residuals( static_cast<Eigen::Index>( targetIndex ), output ) );
                                    if( compactFinalization )
                                    {
                                        compactResiduals[region]( static_cast<Eigen::Index>( targetIndex ),
                                                                  compactColumn ) = residual;
                                        compactValidity[region]( static_cast<Eigen::Index>( targetIndex ),
                                                                 compactColumn ) = 1;
                                    }
                                    else
                                    {
                                        this->m_psfsub[output].image( image )( row, column ) = residual;
                                        this->m_psfsubValidity[output].image( image )( row, column ) = 1;
                                    }
                                }
                                if( !allTargetsValid )
                                {
                                    ++threadRankInvalid[output];
                                }
                            }
                            threadProjectionSeconds += omp_get_wtime() - residualApplyBegin;

                            if( sharedPSF )
                            {
                                const double psfBegin = omp_get_wtime();
                                for( std::size_t output = 0; output < modes.size(); ++output )
                                {
                                    if( result.modeStatus[output] == P4PCAModeStatus::rankInsufficient )
                                    {
                                        continue;
                                    }
                                    const Eigen::Index compactColumn =
                                        static_cast<Eigen::Index>( search * modes.size() + output );
                                    psfModel->calculateLocalResponse(
                                        localPSFResponse,
                                        grid,
                                        search,
                                        coefficients.col( static_cast<Eigen::Index>( output ) )
                                            .head( static_cast<Eigen::Index>( grid.predictorCount() ) ) );
                                    for( int stampColumn = 0; stampColumn < m_localPSFColumns; ++stampColumn )
                                    {
                                        for( int stampRow = 0; stampRow < m_localPSFRows; ++stampRow )
                                        {
                                            const Eigen::Index stampPixel =
                                                static_cast<Eigen::Index>( stampRow + m_localPSFRows * stampColumn );
                                            m_localPSFModels[region]( stampPixel, compactColumn ) =
                                                localPSFResponse( stampRow, stampColumn );
                                        }
                                    }
                                    for( Eigen::Index coefficient = static_cast<Eigen::Index>( grid.predictorCount() );
                                         coefficient < coefficients.rows();
                                         ++coefficient )
                                    {
                                        const double value =
                                            -coefficients( coefficient, static_cast<Eigen::Index>( output ) );
                                        const float stored = static_cast<float>( value );
                                        if( !mx::math::isFinite( stored ) )
                                        {
                                            throw std::overflow_error(
                                                "P4 temporal PSF coefficient exceeds float storage range" );
                                        }
                                        m_localPSFTemporalCoefficients[region](
                                            coefficient - static_cast<Eigen::Index>( grid.predictorCount() ),
                                            compactColumn ) = stored;
                                    }
                                    m_localPSFValidity[region]( static_cast<Eigen::Index>( search ),
                                                                static_cast<Eigen::Index>( output ) ) = 1;
                                }
                                threadPSFSeconds += omp_get_wtime() - psfBegin;
                            }
                        }
                        progressWatcher.incrementAndOutputStatus();
                    }
                    catch( ... )
                    {
                        // clang-format off
#pragma omp critical(P4ReductionException)
                        // clang-format on
                        {
                            if( !workerException )
                            {
                                workerException = std::current_exception();
                            }
                        }
                        failed.store( true, std::memory_order_release );
                    }
                }
            }

            // clang-format off
#pragma omp critical(P4ReductionStatistics)
            // clang-format on
            {
                validLocalFitCount += threadValid;
                maskedLocalFitCount += threadMasked;
                supportInvalidLocalFitCount += threadSupportInvalid;
                minimumNumericalRank = std::min( minimumNumericalRank, threadMinimumRank );
                minimumBaseRank = std::min( minimumBaseRank, threadMinimumBaseRank );
                downdateClampCount += threadDowndateClampCount;
                explicitFallbackCount += threadExplicitFallbackCount;
                rankBoundaryFallbackPixelCount += threadRankBoundaryFallbackPixelCount;
                if( threadFactorValidationFallbackPixelCount != 0 &&
                    ( factorValidationFallbackPixelCount == 0 ||
                      threadMaximumFactorOrthogonalityDefect > maximumFactorOrthogonalityDefect ) )
                {
                    maximumFactorOrthogonalityDefect = threadMaximumFactorOrthogonalityDefect;
                    factorOrthogonalityToleranceAtMaximumDefect = threadFactorOrthogonalityToleranceAtMaximumDefect;
                }
                factorValidationFallbackPixelCount += threadFactorValidationFallbackPixelCount;
                deletionSolverFallbackPixelCount += threadDeletionSolverFallbackPixelCount;
                for( std::size_t output = 0; output < modes.size(); ++output )
                {
                    rankInvalidCounts[output] += threadRankInvalid[output];
                }
                m_timing.samplingWorkerSeconds += threadSamplingSeconds;
                m_timing.sameImageSamplingWorkerSeconds += threadSameImageSamplingSeconds;
                m_timing.temporalSamplingWorkerSeconds += threadTemporalSamplingSeconds;
                m_timing.gramWorkerSeconds += threadGramSeconds;
                m_timing.eigensolveWorkerSeconds += threadEigensolveSeconds;
                m_timing.baseFactorWorkerSeconds += threadBaseFactorSeconds;
                m_timing.deletionWorkerSeconds += threadDeletionSeconds;
                m_timing.explicitFallbackWorkerSeconds += threadExplicitFallbackSeconds;
                m_timing.projectionWorkerSeconds += threadProjectionSeconds;
                m_timing.psfWorkerSeconds += threadPSFSeconds;
            }
        }

        if( workerException )
        {
            progressOutput.interruptRegion();
            try
            {
                std::rethrow_exception( workerException );
            }
            catch( ... )
            {
                std::throw_with_nested(
                    mx::exception<verboseT>( mx::error_t::exception,
                                             "P4 worker failed in region " + std::to_string( region ) ) );
            }
        }

        completedSearchPixels += searchPixelCount;
        progressOutput.completeRegion( completedSearchPixels );

        P4RegionStatistics &statistics = m_regionStatistics[region];
        statistics.validLocalFitCount = validLocalFitCount;
        statistics.maskedLocalFitCount = maskedLocalFitCount;
        statistics.supportInvalidLocalFitCount = supportInvalidLocalFitCount;
        statistics.minimumNumericalRank =
            minimumNumericalRank == std::numeric_limits<int>::max() ? 0 : minimumNumericalRank;
        statistics.minimumBaseRank = minimumBaseRank == std::numeric_limits<int>::max() ? 0 : minimumBaseRank;
        statistics.downdateClampCount = downdateClampCount;
        statistics.explicitFallbackCount = explicitFallbackCount;
        statistics.rankBoundaryFallbackPixelCount = rankBoundaryFallbackPixelCount;
        statistics.factorValidationFallbackPixelCount = factorValidationFallbackPixelCount;
        statistics.deletionSolverFallbackPixelCount = deletionSolverFallbackPixelCount;
        statistics.maximumFactorOrthogonalityDefect = maximumFactorOrthogonalityDefect;
        statistics.factorOrthogonalityToleranceAtMaximumDefect = factorOrthogonalityToleranceAtMaximumDefect;
        statistics.rankInvalidCounts = rankInvalidCounts;
        if( explicitFallbackCount != 0 )
        {
            std::ostringstream fallbackWarning;
            fallbackWarning << "WARNING: P4 region " << region << " recomputed " << explicitFallbackCount
                            << " target rows with the explicit oracle at " << rankBoundaryFallbackPixelCount
                            << " rank-boundary search pixels and " << factorValidationFallbackPixelCount
                            << " factor-validation search pixels and " << deletionSolverFallbackPixelCount
                            << " deletion-solver search pixels";
            if( factorValidationFallbackPixelCount != 0 )
            {
                fallbackWarning << "; maximum factor defect " << std::setprecision( 17 )
                                << maximumFactorOrthogonalityDefect << " at tolerance "
                                << factorOrthogonalityToleranceAtMaximumDefect;
            }
            std::cerr << fallbackWarning.str() << '\n';
        }
        for( std::size_t output = 0; output < modes.size(); ++output )
        {
            if( rankInvalidCounts[output] != 0 )
            {
                std::cerr << "WARNING: P4 region " << region << " output " << output << " requested " << modes[output]
                          << " modes above numerical rank for one or more targets at " << rankInvalidCounts[output]
                          << " search pixels\n";
            }
        }
    }
    m_timing.regressionElapsedSeconds = omp_get_wtime() - regressionBegin;
    progressOutput.completeReduction( m_timing.regressionElapsedSeconds );

    if( m_writeDiagnostics )
    {
        const imageT ownershipDiagnostic = m_ownership.template cast<realT>();
        writeDiagnostic( "p4Ownership.fits", ownershipDiagnostic );
        for( std::size_t output = 0; output < this->m_psfsubValidity.size(); ++output )
        {
            writeDiagnostic( "p4Validity_" + p4Index( output ) + ".fits", this->m_psfsubValidity[output] );
        }
        for( std::size_t region = 0; region < m_temporalSelections.size(); ++region )
        {
            const std::vector<std::vector<int>> &selection = m_temporalSelections[region];
            imageT selectionDiagnostic( static_cast<Eigen::Index>( selection.size() ),
                                        static_cast<Eigen::Index>( selection.front().size() ) );
            for( std::size_t target = 0; target < selection.size(); ++target )
            {
                for( std::size_t image = 0; image < selection[target].size(); ++image )
                {
                    selectionDiagnostic( static_cast<Eigen::Index>( target ), static_cast<Eigen::Index>( image ) ) =
                        static_cast<realT>( selection[target][image] );
                }
            }
            writeDiagnostic( "p4TemporalSelection_" + p4Index( region ) + ".fits", selectionDiagnostic );
        }
        for( std::size_t region = 0; region < grids.size(); ++region )
        {
            int offsetRadius{ 0 };
            const std::size_t predictorCount = m_regressionFrame == P4RegressionFrame::detector
                                                   ? grids[region].predictorCount()
                                                   : rotatedGrids[region].predictorCount();
            const auto predictorOffset = [&]( std::size_t predictor ) -> const P4PixelCoordinate &
            {
                return m_regressionFrame == P4RegressionFrame::detector
                           ? grids[region].predictorOffset( predictor )
                           : rotatedGrids[region].predictorOffset( predictor );
            };
            for( std::size_t predictor = 0; predictor < predictorCount; ++predictor )
            {
                const P4PixelCoordinate &offset = predictorOffset( predictor );
                offsetRadius = std::max( { offsetRadius, std::abs( offset.row() ), std::abs( offset.column() ) } );
            }
            imageT canonical = imageT::Zero( 2 * offsetRadius + 1, 2 * offsetRadius + 1 );
            for( std::size_t predictor = 0; predictor < predictorCount; ++predictor )
            {
                const P4PixelCoordinate &offset = predictorOffset( predictor );
                canonical( offset.row() + offsetRadius, offset.column() + offsetRadius ) = 1;
            }
            writeDiagnostic( "p4CanonicalOR_" + p4Index( region ) + ".fits", canonical );
        }

        imageT summary( static_cast<Eigen::Index>( m_regionStatistics.size() ),
                        static_cast<Eigen::Index>( 20 + 2 * m_modeFractions.size() ) );
        for( std::size_t region = 0; region < m_regionStatistics.size(); ++region )
        {
            const P4RegionStatistics &statistics = m_regionStatistics[region];
            summary( region, 0 ) = m_minRadius[region];
            summary( region, 1 ) = m_maxRadius[region];
            summary( region, 2 ) = static_cast<double>( statistics.searchPixelCount );
            summary( region, 3 ) = static_cast<double>( statistics.targetImageCount );
            summary( region, 4 ) = statistics.temporalNumberImages;
            summary( region, 5 ) = statistics.temporalPsfRadius;
            summary( region, 6 ) = static_cast<double>( statistics.predictorCount );
            summary( region, 7 ) = statistics.maximumDegreesOfFreedom;
            summary( region, 8 ) = statistics.minimumNumericalRank;
            summary( region, 9 ) = static_cast<realT>( statistics.validLocalFitCount );
            summary( region, 10 ) = static_cast<realT>( statistics.maskedLocalFitCount );
            summary( region, 11 ) = static_cast<realT>( statistics.supportInvalidLocalFitCount );
            for( std::size_t output = 0; output < m_modeFractions.size(); ++output )
            {
                summary( region, 12 + output ) = static_cast<realT>( m_realizedModes[region][output] );
                summary( region, 12 + m_modeFractions.size() + output ) =
                    static_cast<realT>( statistics.rankInvalidCounts[output] );
            }
            const std::size_t workerColumn = 12 + 2 * m_modeFractions.size();
            summary( region, workerColumn ) = static_cast<realT>( statistics.estimatedWorkerBytes );
            summary( region, workerColumn + 1 ) = static_cast<realT>( statistics.maximumWorkerCount );
            summary( region, workerColumn + 2 ) = static_cast<realT>( statistics.effectiveWorkerCount );
            summary( region, workerColumn + 3 ) = static_cast<realT>( statistics.rankBoundaryFallbackPixelCount );
            summary( region, workerColumn + 4 ) = static_cast<realT>( statistics.factorValidationFallbackPixelCount );
            summary( region, workerColumn + 5 ) = static_cast<realT>( statistics.deletionSolverFallbackPixelCount );
            summary( region, workerColumn + 6 ) = static_cast<realT>( statistics.maximumFactorOrthogonalityDefect );
            summary( region, workerColumn + 7 ) =
                static_cast<realT>( statistics.factorOrthogonalityToleranceAtMaximumDefect );
        }
        fitsHeaderT summaryHeader;
        summaryHeader.template append<int>( "P4 REGION SUMMARY SCHEMA",
                                            2,
                                            "P4 region-summary diagnostic schema version" );
        summaryHeader.template append<int>( "P4 REGION SUMMARY MODE COUNT",
                                            static_cast<int>( m_modeFractions.size() ),
                                            "output-plane entries in each indexed column group" );
        std::ostringstream summaryColumns;
        summaryColumns << "minRadius,maxRadius,searchPixelCount,targetImageCount,temporalNumberImages,"
                          "temporalPsfRadius,predictorCount,maximumDegreesOfFreedom,minimumNumericalRank,"
                          "validLocalFitCount,maskedLocalFitCount,supportInvalidLocalFitCount";
        for( std::size_t output = 0; output < m_modeFractions.size(); ++output )
        {
            summaryColumns << ",realizedModes[" << p4Index( output ) << ']';
        }
        for( std::size_t output = 0; output < m_modeFractions.size(); ++output )
        {
            summaryColumns << ",rankInvalidCounts[" << p4Index( output ) << ']';
        }
        summaryColumns << ",estimatedWorkerBytes,maximumWorkerCount,effectiveWorkerCount,"
                          "rankBoundaryFallbackPixelCount,factorValidationFallbackPixelCount,"
                          "deletionSolverFallbackPixelCount,maximumFactorOrthogonalityDefect,"
                          "factorOrthogonalityToleranceAtMaximumDefect";
        summaryHeader.template append<std::string>( "P4 REGION SUMMARY COLUMNS",
                                                    summaryColumns.str(),
                                                    "exact columns; [NNN] indexes output-plane groups" );
        writeDiagnostic( "p4RegionSummary.fits", summary, &summaryHeader );
    }

    const auto writeTimingDiagnostic = [&]()
    {
        if( !m_writeDiagnostics )
        {
            return;
        }
        imageT timing( 1, processPSF ? 13 : ( calculatePSF ? 12 : 11 ) );
        timing( 0, 0 ) = static_cast<realT>( m_timing.geometryElapsedSeconds );
        timing( 0, 1 ) = static_cast<realT>( m_timing.regressionElapsedSeconds );
        timing( 0, 2 ) = static_cast<realT>( m_timing.samplingWorkerSeconds );
        timing( 0, 3 ) = static_cast<realT>( m_timing.sameImageSamplingWorkerSeconds );
        timing( 0, 4 ) = static_cast<realT>( m_timing.temporalSamplingWorkerSeconds );
        timing( 0, 5 ) = static_cast<realT>( m_timing.gramWorkerSeconds );
        timing( 0, 6 ) = static_cast<realT>( m_timing.eigensolveWorkerSeconds );
        timing( 0, 7 ) = static_cast<realT>( m_timing.baseFactorWorkerSeconds );
        timing( 0, 8 ) = static_cast<realT>( m_timing.deletionWorkerSeconds );
        timing( 0, 9 ) = static_cast<realT>( m_timing.explicitFallbackWorkerSeconds );
        timing( 0, 10 ) = static_cast<realT>( m_timing.projectionWorkerSeconds );
        if( calculatePSF )
        {
            timing( 0, 11 ) = static_cast<realT>( m_timing.psfWorkerSeconds );
        }
        if( processPSF )
        {
            timing( 0, 12 ) = static_cast<realT>( m_timing.psfReconstructionElapsedSeconds );
        }
        fitsHeaderT timingHeader;
        timingHeader.template append<int>( "P4 TIMING SCHEMA",
                                           processPSF ? 8 : ( calculatePSF ? 7 : 6 ),
                                           "P4 timing diagnostic schema version" );
        timingHeader.template append<std::string>( "P4 TIMING COLUMNS",
                                                   "geometryElapsed,regressionElapsed,samplingWorker,"
                                                   "sameImageSamplingWorker,temporalSamplingWorker,gramWorker,"
                                                   "eigensolveWorker,baseFactorWorker,deletionWorker,"
                                                   "explicitFallbackWorker,projectionWorker" +
                                                       std::string( calculatePSF ? ",psfWorker" : "" ) +
                                                       std::string( processPSF ? ",psfReconstructionElapsed" : "" ),
                                                   "P4 timing columns in seconds" );
        writeDiagnostic( "p4Timing.fits", timing, &timingHeader );
    };

    if( targetHeldOutPSF && !processPSF )
    {
        std::vector<std::size_t> regionOffsets( grids.size(), 0 );
        std::size_t regionOffset{ 0 };
        for( std::size_t region = 0; region < grids.size(); ++region )
        {
            regionOffsets[region] = regionOffset;
            regionOffset += grids[region].searchPixelCount();
        }
        for( std::size_t firstOutput = 0; firstOutput < m_modeFractions.size(); firstOutput += m_psfModeBatchSize )
        {
            const std::size_t outputCount = std::min( m_psfModeBatchSize, m_modeFractions.size() - firstOutput );
            std::vector<imageT> unpublishedModels;
            std::vector<psfValidityT> unpublishedValidity;
            calculateHeldOutPSFBatch( unpublishedModels,
                                      unpublishedValidity,
                                      grids,
                                      *psfModel,
                                      regionExclusions,
                                      regionOffsets,
                                      totalSearchPixels,
                                      firstOutput,
                                      outputCount );
        }
    }

    fitsHeaderT finalHeader;
    this->stdFitsHeader( &finalHeader );
    appendReductionHeader( finalHeader );
    const std::string finalImagePath = this->finalImageOutputPath();

    if( !compactFinalization )
    {
        const int result = finalProcess();
        if( result == 0 && processPSF )
        {
            const double reconstructionBegin = omp_get_wtime();
            processPSFProducts( grids, *psfModel, regionExclusions, finalImagePath, finalHeader );
            m_timing.psfReconstructionElapsedSeconds = omp_get_wtime() - reconstructionBegin;
        }
        if( result == 0 )
        {
            cropAutomaticFinalImage();
            if( m_automaticOutputSize != 0 && this->m_doWriteFinim && this->m_combineMethod != HCI::combine::none )
            {
                std::cerr << "writing\n";
                this->writeFinimAtPath( finalImagePath, &finalHeader );
            }
        }
        writeTimingDiagnostic();
        this->t_end = mx::sys::get_curr_time();
        return result;
    }

    const int combinedImageSize =
        m_automaticDerotationSize != 0 && !m_psfFilter ? m_automaticDerotationSize : this->m_Nrows;
    eigenCube<realT> combinedImages;
    combinedImages.resize( combinedImageSize, combinedImageSize, static_cast<int>( m_modeFractions.size() ) );
    const int configuredWriteFinim = this->m_doWriteFinim;
    const bool configuredOutputPSFSub = this->m_doOutputPSFSub;
    this->m_doWriteFinim = false;
    this->m_doOutputPSFSub = false;
    double derotationSeconds{ 0 };
    double combinationSeconds{ 0 };
    int result{ 0 };
    if( this->m_postMedSub )
    {
        std::cerr << "Subtracting medians in post\n";
    }
    if( m_regressionFrame == P4RegressionFrame::detector && this->m_doDerotate )
    {
        std::cerr << "derotating\n";
    }
    std::cerr << "combining\n";
    try
    {
        for( std::size_t output = 0; output < m_modeFractions.size(); ++output )
        {
            this->m_psfsub.resize( 1 );
            this->m_psfsubValidity.resize( 1 );
            this->m_psfsub[0].resize( this->m_Nrows, this->m_Ncols, this->m_Nims );
            this->m_psfsub[0].setZero();
            this->m_psfsubValidity[0].resize( this->m_Nrows, this->m_Ncols, this->m_Nims );
            this->m_psfsubValidity[0].setZero();

            for( std::size_t region = 0; region < m_regionStatistics.size(); ++region )
            {
                const std::size_t searchPixelCount = m_regionStatistics[region].searchPixelCount;
                const std::vector<std::vector<int>> &temporalSelections = m_temporalSelections[region];
                for( std::size_t search = 0; search < searchPixelCount; ++search )
                {
                    const P4PixelCoordinate &coordinate = m_regressionFrame == P4RegressionFrame::detector
                                                              ? grids[region].searchPixel( search ).coordinate()
                                                              : rotatedGrids[region].searchPixel( search ).coordinate();
                    const Eigen::Index compactColumn =
                        static_cast<Eigen::Index>( search * m_modeFractions.size() + output );
                    for( std::size_t targetIndex = 0; targetIndex < temporalSelections.size(); ++targetIndex )
                    {
                        if( compactValidity[region]( static_cast<Eigen::Index>( targetIndex ), compactColumn ) == 0 )
                        {
                            continue;
                        }
                        const int image = temporalSelections[targetIndex][0];
                        this->m_psfsub[0].image( image )( coordinate.row(), coordinate.column() ) =
                            compactResiduals[region]( static_cast<Eigen::Index>( targetIndex ), compactColumn );
                        this->m_psfsubValidity[0].image( image )( coordinate.row(), coordinate.column() ) = 1;
                    }
                }
            }

            if( m_writeDiagnostics )
            {
                writeDiagnostic( "p4Validity_" + p4Index( output ) + ".fits", this->m_psfsubValidity[0] );
            }

            result = finalProcess( false );
            if( result != 0 )
            {
                break;
            }
            derotationSeconds += this->t_derotate_end - this->t_derotate_begin;
            combinationSeconds += this->t_combo_end - this->t_combo_begin;
            combinedImages.image( static_cast<int>( output ) ) = this->m_finim.image( 0 );
            this->m_psfsub.clear();
            this->m_psfsubValidity.clear();
        }
    }
    catch( ... )
    {
        this->m_doWriteFinim = configuredWriteFinim;
        this->m_doOutputPSFSub = configuredOutputPSFSub;
        throw;
    }
    this->m_doWriteFinim = configuredWriteFinim;
    this->m_doOutputPSFSub = configuredOutputPSFSub;
    if( result == 0 )
    {
        this->m_finim = std::move( combinedImages );
        const double timingReference = mx::sys::get_curr_time();
        this->t_derotate_begin = timingReference;
        this->t_derotate_end = timingReference + derotationSeconds;
        this->t_combo_begin = timingReference;
        this->t_combo_end = timingReference + combinationSeconds;
        if( processPSF )
        {
            const double reconstructionBegin = omp_get_wtime();
            processPSFProducts( grids, *psfModel, regionExclusions, finalImagePath, finalHeader );
            m_timing.psfReconstructionElapsedSeconds = omp_get_wtime() - reconstructionBegin;
        }
        cropAutomaticFinalImage();
        writeTimingDiagnostic();
        if( configuredWriteFinim )
        {
            std::cerr << "writing\n";
            this->writeFinimAtPath( finalImagePath, &finalHeader );
        }
    }
    else
    {
        writeTimingDiagnostic();
    }
    this->t_end = mx::sys::get_curr_time();
    return result;
}

template <typename realT, class derotFunctObj, class verboseT>
void P4Reduction<realT, derotFunctObj, verboseT>::calculateHeldOutPSFBatch(
    std::vector<imageT> &localModels,
    std::vector<psfValidityT> &localValidity,
    const std::vector<pixelGridT> &grids,
    const P4PSFModel &psfModel,
    const std::vector<P4TargetExclusions> &regionExclusions,
    const std::vector<std::size_t> &regionOffsets,
    std::size_t searchPixelCount,
    std::size_t firstOutput,
    std::size_t outputCount )
{
#ifdef HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION
    const std::optional<P4ReductionPCADispatch> experimentalDispatch = p4ReductionSelectedPCADispatch();
#endif
    if( outputCount == 0 || firstOutput >= m_modeFractions.size() ||
        outputCount > m_modeFractions.size() - firstOutput || grids.size() != m_regionStatistics.size() ||
        regionExclusions.size() != grids.size() || regionOffsets.size() != grids.size() || searchPixelCount == 0 )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr,
                                       "P4 held-out PSF batch dimensions do not match the reduction state" );
    }
    const std::size_t targetCount = static_cast<std::size_t>( this->m_Nims );
    const Eigen::Index localPixels =
        static_cast<Eigen::Index>( m_localPSFRows ) * static_cast<Eigen::Index>( m_localPSFColumns );
    if( targetCount == 0 || outputCount > std::numeric_limits<std::size_t>::max() / targetCount ||
        targetCount * outputCount > static_cast<std::size_t>( std::numeric_limits<Eigen::Index>::max() ) ||
        searchPixelCount > std::numeric_limits<std::size_t>::max() / targetCount ||
        searchPixelCount * targetCount > static_cast<std::size_t>( std::numeric_limits<Eigen::Index>::max() ) )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr,
                                       "P4 held-out PSF batch dimensions exceed the Eigen index range" );
    }
    const Eigen::Index modelColumns = static_cast<Eigen::Index>( searchPixelCount * targetCount );
    localModels.resize( outputCount );
    localValidity.resize( outputCount );
    for( std::size_t output = 0; output < outputCount; ++output )
    {
        localModels[output].resize( localPixels, modelColumns );
        localModels[output].setZero();
        localValidity[output].resize( static_cast<Eigen::Index>( searchPixelCount ),
                                      static_cast<Eigen::Index>( targetCount ) );
        localValidity[output].setZero();
    }

    const imageT *mask = this->m_mask.size() == 0 ? nullptr : &this->m_mask;
    const std::vector<P4PixelCoordinate> temporalOffsets = p4TemporalPredictorOffsets( m_psfRadius );
    for( std::size_t region = 0; region < grids.size(); ++region )
    {
        const pixelGridT &grid = grids[region];
        const P4RegionStatistics &statistics = m_regionStatistics[region];
        const std::vector<std::vector<int>> &selections = m_temporalSelections[region];
        if( selections.size() != targetCount || regionOffsets[region] > searchPixelCount ||
            grid.searchPixelCount() > searchPixelCount - regionOffsets[region] ||
            regionExclusions[region].targetCount() != static_cast<Eigen::Index>( targetCount ) )
        {
            throw mx::exception<verboseT>( mx::error_t::sizeerr,
                                           "P4 held-out PSF annulus does not retain one row per target frame" );
        }
        for( std::size_t target = 0; target < targetCount; ++target )
        {
            if( selections[target].size() != 1 || selections[target][0] != static_cast<int>( target ) )
            {
                throw mx::exception<verboseT>(
                    mx::error_t::invalidconfig,
                    "P4 held-out PSF reconstruction requires p4.numberImages=0 and ordered target rows" );
            }
        }
        std::vector<int> modes;
        modes.reserve( outputCount );
        for( std::size_t output = 0; output < outputCount; ++output )
        {
            modes.push_back( m_realizedModes[region][firstOutput + output] );
        }

        const std::size_t maximumDeletedRows = m_exclusionSolver == P4ExclusionSolver::factorDowndateExact &&
                                                       m_deletionBackend == mx::math::svdDeletionBackend::rankOneSecular
                                                   ? 1
                                                   : statistics.maximumExcludedRows;
        const std::size_t workerBytes = estimatedWorkerBytes( targetCount,
                                                              statistics.predictorCount,
                                                              outputCount,
                                                              false,
                                                              0,
                                                              true,
                                                              m_exclusionSolver,
                                                              maximumDeletedRows,
                                                              static_cast<std::size_t>( modes.back() ),
                                                              static_cast<std::size_t>( localPixels ) );
        const int requestedWorkers =
            std::max( 1,
                      std::min( omp_get_max_threads(),
                                static_cast<int>( std::min<std::size_t>(
                                    grid.searchPixelCount(),
                                    static_cast<std::size_t>( std::numeric_limits<int>::max() ) ) ) ) );
        int effectiveWorkers = requestedWorkers;
        if( m_memoryBudgetBytes != 0 )
        {
            const long double persistent = static_cast<long double>( m_compactResidualBytes ) + m_targetExclusionBytes +
                                           m_psfModelBytes + m_localPSFBytes + m_psfReconstructionBytes;
            if( persistent > static_cast<long double>( std::numeric_limits<std::size_t>::max() ) )
            {
                throw mx::exception<verboseT>( mx::error_t::sizeerr, "P4 held-out PSF persistent byte count overflow" );
            }
            effectiveWorkers = memoryLimitedWorkerCount( requestedWorkers,
                                                         m_memoryBudgetBytes,
                                                         static_cast<std::size_t>( persistent ),
                                                         workerBytes );
            if( effectiveWorkers == 0 )
            {
                throw mx::exception<verboseT>( mx::error_t::allocerr,
                                               "one P4 held-out PSF worker does not fit the automatic memory budget" );
            }
        }
        std::cerr << "P4 held-out PSF batch " << firstOutput + 1 << '-' << firstOutput + outputCount << " / "
                  << m_modeFractions.size() << ", annulus " << region + 1 << " / " << grids.size() << ": workers "
                  << effectiveWorkers << " / " << requestedWorkers << '\n';

        std::exception_ptr workerException;
        std::atomic<bool> failed{ false };
        double responseWorkerSeconds{ 0 };
        std::size_t batchDowndateClampCount{ 0 };
        std::size_t batchExplicitFallbackCount{ 0 };
        std::size_t batchRankBoundaryFallbackPixelCount{ 0 };
        std::size_t batchFactorValidationFallbackPixelCount{ 0 };
        std::size_t batchDeletionSolverFallbackPixelCount{ 0 };
        double batchMaximumFactorDefect{ 0 };
        double batchFactorToleranceAtMaximumDefect{ 0 };
        // clang-format off
#pragma omp parallel num_threads(effectiveWorkers)
        // clang-format on
        {
            P4PCA::workspaceT workspace;
            detail::P4PCAMixedWorkspace mixedWorkspace;
#ifdef HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION
            std::optional<detail::P4PCAExperimentalWorkspace> experimentalWorkspace;
            std::optional<P4ReductionWorkerPrecisionScope> experimentalPrecisionScope;
            if( experimentalDispatch )
            {
                experimentalWorkspace.emplace();
                experimentalPrecisionScope.emplace( *experimentalDispatch, *experimentalWorkspace );
            }
#endif
            P4PCADowndateWorkspace downdateWorkspace;
            P4PCA::matrixT predictors;
            P4PCA::vectorT target;
            P4PCA::matrixT probePredictors;
            P4PCA::vectorT probeTarget;
            P4PCA::matrixT probeResiduals;
            P4PCAResult result;
            double threadResponseSeconds{ 0 };
            std::size_t threadDowndateClampCount{ 0 };
            std::size_t threadExplicitFallbackCount{ 0 };
            std::size_t threadRankBoundaryFallbackPixelCount{ 0 };
            std::size_t threadFactorValidationFallbackPixelCount{ 0 };
            std::size_t threadDeletionSolverFallbackPixelCount{ 0 };
            double threadMaximumFactorDefect{ 0 };
            double threadFactorToleranceAtMaximumDefect{ 0 };

            // clang-format off
#pragma omp for schedule(static)
            // clang-format on
            for( std::size_t search = 0; search < grid.searchPixelCount(); ++search )
            {
                if( failed.load( std::memory_order_acquire ) )
                {
                    continue;
                }
                try
                {
                    const P4PixelCoordinate &coordinate = grid.searchPixel( search ).coordinate();
                    const bool temporalValid =
                        statistics.temporalNumberImages == 0 ||
                        p4TemporalPredictorsValid( coordinate, temporalOffsets, this->m_Nrows, this->m_Ncols, mask );
                    if( !grid.searchPixel( search ).valid() || !temporalValid )
                    {
                        continue;
                    }
                    const double responseBegin = omp_get_wtime();
                    psfModel.responseInputs( probeTarget, probePredictors, grid, search );
                    double sameImageSamplingSeconds{ 0 };
                    double temporalSamplingSeconds{ 0 };
                    P4PCATiming timing;
                    fitDetectorSearch( result,
                                       predictors,
                                       target,
                                       nullptr,
                                       &probeResiduals,
                                       grid,
                                       search,
                                       selections,
                                       temporalOffsets,
                                       modes,
                                       workspace,
                                       mixedWorkspace,
                                       &downdateWorkspace,
                                       timing,
                                       sameImageSamplingSeconds,
                                       temporalSamplingSeconds,
                                       &regionExclusions[region],
                                       &probePredictors,
                                       &probeTarget );
                    threadDowndateClampCount += result.downdateClampCount;
                    threadExplicitFallbackCount += result.explicitFallbackCount;
                    switch( result.explicitFallbackReason )
                    {
                    case P4PCAFallbackReason::none:
                        break;
                    case P4PCAFallbackReason::rankBoundary:
                        ++threadRankBoundaryFallbackPixelCount;
                        break;
                    case P4PCAFallbackReason::factorValidation:
                        ++threadFactorValidationFallbackPixelCount;
                        if( threadFactorValidationFallbackPixelCount == 1 ||
                            result.factorOrthogonalityDefect > threadMaximumFactorDefect )
                        {
                            threadMaximumFactorDefect = result.factorOrthogonalityDefect;
                            threadFactorToleranceAtMaximumDefect = result.factorOrthogonalityTolerance;
                        }
                        break;
                    case P4PCAFallbackReason::deletionSolver:
                        ++threadDeletionSolverFallbackPixelCount;
                        break;
                    }
                    const std::size_t globalSearch = regionOffsets[region] + search;
                    for( std::size_t targetIndex = 0; targetIndex < targetCount; ++targetIndex )
                    {
                        const Eigen::Index modelColumn =
                            static_cast<Eigen::Index>( globalSearch * targetCount + targetIndex );
                        for( std::size_t output = 0; output < outputCount; ++output )
                        {
                            if( !result.sampleSupported( static_cast<Eigen::Index>( targetIndex ), output ) )
                            {
                                continue;
                            }
                            const Eigen::Index probeColumn =
                                static_cast<Eigen::Index>( targetIndex * outputCount + output );
                            for( Eigen::Index pixel = 0; pixel < localPixels; ++pixel )
                            {
                                localModels[output]( pixel, modelColumn ) =
                                    checkedResidualCast( probeResiduals( pixel, probeColumn ) );
                            }
                            localValidity[output]( static_cast<Eigen::Index>( globalSearch ),
                                                   static_cast<Eigen::Index>( targetIndex ) ) = 1;
                        }
                    }
                    threadResponseSeconds += omp_get_wtime() - responseBegin;
                }
                catch( ... )
                {
                    // clang-format off
#pragma omp critical(P4HeldOutPSFException)
                    // clang-format on
                    {
                        if( !workerException )
                        {
                            workerException = std::current_exception();
                        }
                    }
                    failed.store( true, std::memory_order_release );
                }
            }

            // clang-format off
#pragma omp critical(P4HeldOutPSFDiagnostics)
            // clang-format on
            {
                responseWorkerSeconds += threadResponseSeconds;
                batchDowndateClampCount += threadDowndateClampCount;
                batchExplicitFallbackCount += threadExplicitFallbackCount;
                batchRankBoundaryFallbackPixelCount += threadRankBoundaryFallbackPixelCount;
                batchFactorValidationFallbackPixelCount += threadFactorValidationFallbackPixelCount;
                batchDeletionSolverFallbackPixelCount += threadDeletionSolverFallbackPixelCount;
                if( threadFactorValidationFallbackPixelCount != 0 &&
                    ( batchFactorValidationFallbackPixelCount == threadFactorValidationFallbackPixelCount ||
                      threadMaximumFactorDefect > batchMaximumFactorDefect ) )
                {
                    batchMaximumFactorDefect = threadMaximumFactorDefect;
                    batchFactorToleranceAtMaximumDefect = threadFactorToleranceAtMaximumDefect;
                }
            }
        }
        if( workerException )
        {
            try
            {
                std::rethrow_exception( workerException );
            }
            catch( const std::exception &error )
            {
                std::throw_with_nested(
                    mx::exception<verboseT>( mx::error_t::exception,
                                             "P4 held-out PSF response calculation failed in annulus " +
                                                 std::to_string( region ) + ": " + error.what() ) );
            }
            catch( ... )
            {
                std::throw_with_nested( mx::exception<verboseT>(
                    mx::error_t::exception,
                    "P4 held-out PSF response calculation failed in annulus " + std::to_string( region ) ) );
            }
        }
        m_psfDowndateClampCount += batchDowndateClampCount;
        m_psfExplicitFallbackCount += batchExplicitFallbackCount;
        m_psfRankBoundaryFallbackPixelCount += batchRankBoundaryFallbackPixelCount;
        m_psfFactorValidationFallbackPixelCount += batchFactorValidationFallbackPixelCount;
        m_psfDeletionSolverFallbackPixelCount += batchDeletionSolverFallbackPixelCount;
        if( batchFactorValidationFallbackPixelCount != 0 &&
            ( m_psfFactorValidationFallbackPixelCount == batchFactorValidationFallbackPixelCount ||
              batchMaximumFactorDefect > m_psfMaximumFactorOrthogonalityDefect ) )
        {
            m_psfMaximumFactorOrthogonalityDefect = batchMaximumFactorDefect;
            m_psfFactorOrthogonalityToleranceAtMaximumDefect = batchFactorToleranceAtMaximumDefect;
        }
        if( batchExplicitFallbackCount != 0 )
        {
            std::cerr << "WARNING: P4 held-out PSF batch " << firstOutput + 1 << '-' << firstOutput + outputCount
                      << " annulus " << region + 1 << " / " << grids.size() << " recomputed "
                      << batchExplicitFallbackCount << " target rows with the explicit oracle at "
                      << batchRankBoundaryFallbackPixelCount << " rank-boundary search pixels and "
                      << batchFactorValidationFallbackPixelCount << " factor-validation search pixels and "
                      << batchDeletionSolverFallbackPixelCount << " deletion-solver search pixels";
            if( batchFactorValidationFallbackPixelCount != 0 )
            {
                std::cerr << "; maximum factor defect " << std::setprecision( 17 ) << batchMaximumFactorDefect
                          << " at tolerance " << batchFactorToleranceAtMaximumDefect;
            }
            std::cerr << '\n';
        }
        m_timing.psfWorkerSeconds += responseWorkerSeconds;
    }
}

template <typename realT, class derotFunctObj, class verboseT>
void P4Reduction<realT, derotFunctObj, verboseT>::processPSFProducts(
    const std::vector<pixelGridT> &grids,
    const P4PSFModel &psfModel,
    const std::vector<P4TargetExclusions> &regionExclusions,
    const std::string &finalImagePath,
    const fitsHeaderT &finalHeader )
{
    const bool targetHeldOutPSF = m_excludeMethod != HCI::exclude::none;
    const bool sharedStateValid = m_localPSFModels.size() == grids.size() &&
                                  m_localPSFTemporalCoefficients.size() == grids.size() &&
                                  m_localPSFValidity.size() == grids.size();
    if( grids.size() != m_regionStatistics.size() || m_localPSFComponentCounts.size() != grids.size() ||
        regionExclusions.size() != grids.size() || ( targetHeldOutPSF ? m_psfModeBatchSize == 0 : !sharedStateValid ) )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr,
                                       "P4 local PSF state does not match detector-region geometry" );
    }
    if( targetHeldOutPSF )
    {
        for( const P4TargetExclusions &exclusions : regionExclusions )
        {
            if( exclusions.empty() )
            {
                throw mx::exception<verboseT>( mx::error_t::sizeerr,
                                               "P4 target-held-out PSF state is missing annular exclusions" );
            }
        }
    }
    if( m_psfFilter && ( this->m_finim.rows() != this->m_Nrows || this->m_finim.cols() != this->m_Ncols ||
                         this->m_finim.planes() != static_cast<int>( m_modeFractions.size() ) ) )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr,
                                       "P4 PSF filtering requires the complete combined final-image cube" );
    }
    if( m_psfFilter && finalImagePath.empty() )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidarg,
                                       "P4 PSF filtering requires a resolved final-image output path" );
    }

    std::size_t searchPixelCount{ 0 };
    for( const P4RegionStatistics &statistics : m_regionStatistics )
    {
        if( searchPixelCount > std::numeric_limits<std::size_t>::max() - statistics.searchPixelCount )
        {
            throw mx::exception<verboseT>( mx::error_t::sizeerr, "P4 final PSF search count overflow" );
        }
        searchPixelCount += statistics.searchPixelCount;
    }
    if( searchPixelCount == 0 || searchPixelCount > static_cast<std::size_t>( std::numeric_limits<int>::max() ) ||
        searchPixelCount > static_cast<std::size_t>( std::numeric_limits<Eigen::Index>::max() ) )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "P4 final PSF search count is outside output range" );
    }

    using reconstructorT = P4PSFReconstructor;
    const std::vector<P4PixelCoordinate> temporalOffsets = p4TemporalPredictorOffsets( m_psfRadius );
    reconstructorT::searchIndexT searchIndex =
        reconstructorT::searchIndexT::Constant( this->m_Nrows, this->m_Ncols, -1 );
    imageT coordinates( static_cast<Eigen::Index>( searchPixelCount ), 4 );
    std::vector<std::size_t> regionOffsets( grids.size(), 0 );
    std::vector<int> searchRegions( searchPixelCount, -1 );
    std::size_t globalSearch{ 0 };
    for( std::size_t region = 0; region < grids.size(); ++region )
    {
        regionOffsets[region] = globalSearch;
        for( std::size_t search = 0; search < grids[region].searchPixelCount(); ++search )
        {
            const P4PixelCoordinate &coordinate = grids[region].searchPixel( search ).coordinate();
            if( searchIndex( coordinate.row(), coordinate.column() ) >= 0 )
            {
                throw mx::exception<verboseT>( mx::error_t::exception,
                                               "P4 final PSF coordinate has multiple local-model owners" );
            }
            searchIndex( coordinate.row(), coordinate.column() ) = static_cast<int>( globalSearch );
            coordinates( static_cast<Eigen::Index>( globalSearch ), 0 ) = static_cast<realT>( coordinate.row() );
            coordinates( static_cast<Eigen::Index>( globalSearch ), 1 ) = static_cast<realT>( coordinate.column() );
            coordinates( static_cast<Eigen::Index>( globalSearch ), 2 ) = static_cast<realT>( region );
            coordinates( static_cast<Eigen::Index>( globalSearch ), 3 ) = static_cast<realT>( search );
            searchRegions[globalSearch] = static_cast<int>( region );
            ++globalSearch;
        }
    }

    std::vector<double> derotationAngles( static_cast<std::size_t>( this->m_Nims ), 0 );
    if( this->m_doDerotate )
    {
        for( int image = 0; image < this->m_Nims; ++image )
        {
            const double angle = static_cast<double>( this->m_derotF.derotAngle( static_cast<std::size_t>( image ) ) );
            if( !mx::math::isFinite( angle ) )
            {
                throw mx::exception<verboseT>( mx::error_t::invalidarg,
                                               "P4 PSF output derotation angles must be finite" );
            }
            derotationAngles[static_cast<std::size_t>( image )] = angle;
        }
    }

    const reconstructorT reconstructor( this->m_Nrows,
                                        this->m_Ncols,
                                        grids.front().xCenter(),
                                        grids.front().yCenter(),
                                        m_psfStampSize,
                                        m_localPSFRows,
                                        m_localPSFColumns );
    const reconstructorT centerReconstructor( this->m_Nrows,
                                              this->m_Ncols,
                                              grids.front().xCenter(),
                                              grids.front().yCenter(),
                                              1,
                                              m_localPSFRows,
                                              m_localPSFColumns );

    const std::size_t componentStride =
        *std::max_element( m_localPSFComponentCounts.begin(), m_localPSFComponentCounts.end() );
    if( temporalOffsets.empty() || componentStride > static_cast<std::size_t>( std::numeric_limits<int>::max() ) ||
        ( componentStride - 1 ) > std::numeric_limits<std::size_t>::max() / temporalOffsets.size() )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "P4 final PSF temporal coefficient dimensions overflow" );
    }
    const std::size_t temporalCoefficientCount = ( componentStride - 1 ) * temporalOffsets.size();
    if( temporalCoefficientCount > static_cast<std::size_t>( std::numeric_limits<int>::max() ) ||
        temporalCoefficientCount > static_cast<std::size_t>( std::numeric_limits<Eigen::Index>::max() ) )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr,
                                       "P4 final PSF temporal coefficient dimensions exceed output range" );
    }
    const std::size_t oneWorkerBytes = psfReconstructionBytes( 1,
                                                               static_cast<std::size_t>( this->m_Nims ),
                                                               temporalCoefficientCount,
                                                               m_psfStampSize,
                                                               m_localPSFRows,
                                                               m_localPSFColumns );
    const std::size_t globalReconstructionBytes =
        m_psfReconstructionBytes > oneWorkerBytes ? m_psfReconstructionBytes - oneWorkerBytes : 0;
    const int requestedWorkers =
        std::max( 1,
                  std::min( omp_get_max_threads(),
                            static_cast<int>( std::min<std::size_t>(
                                searchPixelCount,
                                static_cast<std::size_t>( std::numeric_limits<int>::max() ) ) ) ) );
    int effectiveWorkers = requestedWorkers;
    if( m_memoryBudgetBytes != 0 )
    {
        const long double baseBytes = static_cast<long double>( m_compactResidualBytes ) + m_targetExclusionBytes +
                                      m_localPSFBytes + m_psfModelBytes + globalReconstructionBytes;
        if( baseBytes > static_cast<long double>( m_memoryBudgetBytes ) )
        {
            throw mx::exception<verboseT>( mx::error_t::allocerr,
                                           "P4 PSF reconstruction global state exceeds the automatic memory budget" );
        }
        const std::size_t availableForWorkers = m_memoryBudgetBytes - static_cast<std::size_t>( baseBytes );
        effectiveWorkers = std::min( requestedWorkers,
                                     std::max( 1,
                                               static_cast<int>( std::min<std::size_t>(
                                                   availableForWorkers / oneWorkerBytes,
                                                   static_cast<std::size_t>( std::numeric_limits<int>::max() ) ) ) ) );
    }
    const long double selectedReconstructionBytes = static_cast<long double>( globalReconstructionBytes ) +
                                                    static_cast<long double>( effectiveWorkers ) * oneWorkerBytes;
    if( selectedReconstructionBytes > static_cast<long double>( std::numeric_limits<std::size_t>::max() ) )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "P4 PSF reconstruction worker byte count overflow" );
    }
    const std::size_t selectedReconstructionByteCount = static_cast<std::size_t>( selectedReconstructionBytes );
    std::cerr << "P4 PSF reconstruction workers: " << effectiveWorkers << " / " << requestedWorkers << '\n';

    const std::filesystem::path productDirectory = p4AuxiliaryProductDirectory( finalImagePath );
    const auto productPath = [&]( const std::string &suffix )
    { return ( productDirectory / ( m_psfOutputPrefix + suffix ) ).string(); };
    const auto writeProduct = [&]<typename dataT>( const std::string &path, const dataT &data, fitsHeaderT &header )
    {
        const std::string parent = mx::ioutils::parentPath( path );
        if( !parent.empty() )
        {
            const mx::error_t directoryResult = mx::ioutils::createDirectories( parent );
            if( directoryResult != mx::error_t::noerror )
            {
                throw mx::exception<verboseT>( directoryResult, "could not create P4 PSF product output directory" );
            }
        }
        static std::atomic<unsigned long long> sequence{ 0 };
        const std::string temporaryPath =
            path + ".tmp." + std::to_string( sequence.fetch_add( 1, std::memory_order_relaxed ) );
        mx::fits::fitsFile<realT, verboseT> writer;
        const mx::error_t writeResult = writer.write( temporaryPath, data, header );
        if( writeResult != mx::error_t::noerror )
        {
            std::error_code cleanupError;
            std::filesystem::remove( temporaryPath, cleanupError );
            throw mx::exception<verboseT>( writeResult, "could not write temporary P4 PSF product " + path );
        }
        std::error_code renameError;
        std::filesystem::rename( temporaryPath, path, renameError );
        if( renameError )
        {
            std::error_code cleanupError;
            std::filesystem::remove( temporaryPath, cleanupError );
            throw mx::exception<verboseT>( mx::error_t::fileoerr,
                                           "could not publish P4 PSF product " + path + ": " + renameError.message() );
        }
    };
    const auto productHeader = [&]( const std::string &product, std::size_t output )
    {
        fitsHeaderT header;
        fitsHeaderT finalHeaderCopy( finalHeader );
        this->finalImageHeader( header, &finalHeaderCopy );
        header.template append<int>( "P4 PSF PRODUCT SCHEMA", m_psfFilter ? 2 : 1, "frozen-model PSF product schema" );
        header.template append<std::string>( "P4 PSF PRODUCT", product, "compact PSF product role" );
        header.template append<std::string>( "P4 PSF TEMPLATE", m_psfFile, "post-preprocessing centered template" );
        header.template append<std::string>( "P4 PSF TEMPLATE STAGE", "P4_INPUT", "template processing stage" );
        header.template append<std::string>( "P4 PSF NORMALIZATION", "STORED", "template normalization convention" );
        header.template append<std::string>( "P4 PSF RESPONSE", "FROZEN_SIGNED", "forward-model convention" );
        header.template append<std::string>( "P4 PSF COEFFICIENT SCOPE",
                                             targetHeldOutPSF ? "TARGET_HELD_OUT" : "SHARED_IN_SAMPLE",
                                             "training-row scope of the fitted PSF response" );
        header.template append<int>( "P4 PSF MODE BATCH",
                                     static_cast<int>( targetHeldOutPSF ? m_psfModeBatchSize : 0 ),
                                     "held-out response modes retained per bounded batch" );
        if( product == "MANIFEST" && targetHeldOutPSF )
        {
            header.template append<std::string>( "P4 PSF DOWNDATE CLAMPS",
                                                 std::to_string( m_psfDowndateClampCount ),
                                                 "held-out response deletion eigenvalue clamps" );
            header.template append<std::string>( "P4 PSF EXPLICIT FALLBACK ROWS",
                                                 std::to_string( m_psfExplicitFallbackCount ),
                                                 "held-out response row-refit invocations" );
            header.template append<std::string>( "P4 PSF RANK FALLBACK FITS",
                                                 std::to_string( m_psfRankBoundaryFallbackPixelCount ),
                                                 "held-out response search/batch rank fallbacks" );
            header.template append<std::string>( "P4 PSF FACTOR FALLBACK FITS",
                                                 std::to_string( m_psfFactorValidationFallbackPixelCount ),
                                                 "held-out response search/batch factor fallbacks" );
            header.template append<std::string>( "P4 PSF SOLVER FALLBACK FITS",
                                                 std::to_string( m_psfDeletionSolverFallbackPixelCount ),
                                                 "held-out response search/batch solver fallbacks" );
            header.template append<double>( "P4 PSF MAX FACTOR DEFECT",
                                            m_psfMaximumFactorOrthogonalityDefect,
                                            "largest held-out response factor defect" );
            header.template append<double>( "P4 PSF FACTOR DEFECT TOLERANCE",
                                            m_psfFactorOrthogonalityToleranceAtMaximumDefect,
                                            "tolerance paired with largest PSF factor defect" );
        }
        header.template append<int>( "P4 PSF TEMPLATE ROWS", m_psfTemplateRows, "input template row count" );
        header.template append<int>( "P4 PSF TEMPLATE COLUMNS", m_psfTemplateColumns, "input template column count" );
        header.template append<double>( "P4 PSF TEMPLATE CENTER ROW",
                                        0.5 * static_cast<double>( m_psfTemplateRows - 1 ),
                                        "geometric template-center row" );
        header.template append<double>( "P4 PSF TEMPLATE CENTER COLUMN",
                                        0.5 * static_cast<double>( m_psfTemplateColumns - 1 ),
                                        "geometric template-center column" );
        header.template append<std::string>( "P4 PSF COMBINATION",
                                             HCI::combineToStr<verboseT>( this->m_combineMethod ),
                                             "filter-template combination approximation" );
        header.template append<realT>( "P4 PSF MIN GOOD FRACTION",
                                       this->m_minGoodFract,
                                       "minimum valid-frame fraction" );
        header.template append<realT>( "P4 PSF SIGMA THRESHOLD",
                                       this->m_sigmaThreshold,
                                       "configured response clipping threshold" );
        header.template append<int>( "P4 PSF WEIGHT COUNT",
                                     static_cast<int>( this->m_comboWeights.size() ),
                                     "configured final-combination weights" );
        header.template append<int>( "P4 PSF TEMPORAL NUMBER IMAGES",
                                     m_numberImages,
                                     "additional predictor images per direction" );
        header.template append<int>( "P4 PSF COMPONENT STRIDE",
                                     static_cast<int>( componentStride ),
                                     "maximum retained temporal response components" );
        header.template append<int>( "P4 PSF TEMPORAL COEFFICIENT COUNT",
                                     static_cast<int>( temporalCoefficientCount ),
                                     "retained temporal coefficients per search pixel and mode" );
        header.template append<int>( "P4 PSF STAMP SIZE", m_psfStampSize, "square final PSF stamp size" );
        header.template append<int>( "P4 LOCAL PSF ROWS", m_localPSFRows, "support-padded local response rows" );
        header.template append<int>( "P4 LOCAL PSF COLUMNS",
                                     m_localPSFColumns,
                                     "support-padded local response columns" );
        header.template append<std::string>( "P4 LOCAL PSF BYTES",
                                             std::to_string( m_localPSFBytes ),
                                             "retained compact local PSF bytes" );
        header.template append<std::string>( "P4 PSF RECONSTRUCTION BYTES",
                                             std::to_string( selectedReconstructionByteCount ),
                                             "estimated PSF-product peak scratch" );
        header.template append<int>( "P4 PSF MODEL OUTPUT", m_outputPSFModels ? 1 : 0, "compact models enabled" );
        header.template append<int>( "P4 PSF FILTER", m_psfFilter ? 1 : 0, "normalized filtering enabled" );
        header.template append<realT>( "P4 PSF FILTER MIN GOOD FRACTION",
                                       m_psfFilterMinGoodFract,
                                       "minimum usable stamp fraction" );
        header.template append<std::string>( "P4 PSF FILTER BYTES",
                                             std::to_string( m_psfFilterBytes ),
                                             "retained filter-product bytes" );
        header.template append<std::string>( "P4 PSF FILTER EQUATION",
                                             "SUM(H*I)/SUM(H*H)",
                                             "signed normalized local filter" );
        header.template append<std::string>( "P4 PSF FILTER SOURCE",
                                             "CURRENT_FINAL_IMAGE",
                                             "source final-image identity" );
        header.template append<std::string>( "P4 PSF FILTER SOURCE OUTPUT NAME",
                                             this->m_finimName,
                                             "configured final-image output name" );
        header.template append<std::string>( "P4 PSF FILTER SOURCE PATH",
                                             finalImagePath,
                                             "resolved final-image output path" );
        header.template append<int>( "P4 PSF FILTER SOURCE EXACT NAME",
                                     this->m_exactFinimName ? 1 : 0,
                                     "final-image exact-name policy" );
        header.template append<std::string>( "P4 PSF SOURCE COUNT",
                                             std::to_string( searchPixelCount ),
                                             "coordinate-indexed source positions" );
        header.template append<int>( "P4 PSF MODE COUNT",
                                     static_cast<int>( m_modeFractions.size() ),
                                     "configured compact response products" );
        if( output < m_modeFractions.size() )
        {
            header.template append<int>( "P4 PSF MODE INDEX", static_cast<int>( output ), "zero-based output index" );
            header.template append<realT>( "P4 PSF MODE FRACTION",
                                           m_modeFractions[output],
                                           "configured PCA mode fraction" );
        }
        return header;
    };

    imageT completion( 1, 1 );
    completion( 0, 0 ) = 0;
    fitsHeaderT incompleteHeader = productHeader( "MANIFEST", m_modeFractions.size() );
    incompleteHeader.template append<int>( "P4 PSF COMPLETE", 0, "complete product set available" );
    writeProduct( productPath( "manifest.fits" ), completion, incompleteHeader );

    if( m_outputPSFModels )
    {
        fitsHeaderT coordinateHeader = productHeader( "COORDINATES", m_modeFractions.size() );
        coordinateHeader.template append<std::string>( "P4 PSF COORDINATE COLUMNS",
                                                       "row,column,region,regionSearch",
                                                       "coordinate image columns" );
        writeProduct( productPath( "coordinates.fits" ), coordinates, coordinateHeader );
    }

    eigenCube<realT> filtered;
    eigenCube<realT> filterNormalization;
    eigenCube<realT> filterSupport;
    eigenCube<realT> filterValidity;
    if( m_psfFilter )
    {
        const int outputCount = static_cast<int>( m_modeFractions.size() );
        filtered.resize( this->m_Nrows, this->m_Ncols, outputCount );
        filterNormalization.resize( this->m_Nrows, this->m_Ncols, outputCount );
        filterSupport.resize( this->m_Nrows, this->m_Ncols, outputCount );
        filterValidity.resize( this->m_Nrows, this->m_Ncols, outputCount );
        filtered.cube().setConstant( invalidNumber<realT>() );
        filterNormalization.cube().setConstant( invalidNumber<realT>() );
        filterSupport.cube().setConstant( invalidNumber<realT>() );
        filterValidity.setZero();
    }

    const Eigen::Index localPixels =
        static_cast<Eigen::Index>( m_localPSFRows ) * static_cast<Eigen::Index>( m_localPSFColumns );
    if( searchPixelCount > static_cast<std::size_t>( std::numeric_limits<Eigen::Index>::max() ) )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr, "P4 final PSF component dimensions exceed Eigen range" );
    }
    const std::size_t modeBatchSize = targetHeldOutPSF ? m_psfModeBatchSize : 1;
    for( std::size_t firstOutput = 0; firstOutput < m_modeFractions.size(); firstOutput += modeBatchSize )
    {
        const std::size_t currentBatchSize = std::min( modeBatchSize, m_modeFractions.size() - firstOutput );
        std::vector<imageT> heldOutModels;
        std::vector<psfValidityT> heldOutValidity;
        if( targetHeldOutPSF )
        {
            calculateHeldOutPSFBatch( heldOutModels,
                                      heldOutValidity,
                                      grids,
                                      psfModel,
                                      regionExclusions,
                                      regionOffsets,
                                      searchPixelCount,
                                      firstOutput,
                                      currentBatchSize );
        }
        for( std::size_t batchOutput = 0; batchOutput < currentBatchSize; ++batchOutput )
        {
            const std::size_t output = firstOutput + batchOutput;
            std::cerr << "reconstructing P4 PSF model " << output + 1 << " / " << m_modeFractions.size() << '\n';
            imageT localModels;
            imageT temporalCoefficients( static_cast<Eigen::Index>( temporalCoefficientCount ),
                                         static_cast<Eigen::Index>( searchPixelCount ) );
            temporalCoefficients.setZero();
            reconstructorT::validityT localValidity;
            if( targetHeldOutPSF )
            {
                localModels = std::move( heldOutModels[batchOutput] );
                localValidity = std::move( heldOutValidity[batchOutput] );
            }
            else
            {
                localModels.resize( localPixels, static_cast<Eigen::Index>( searchPixelCount ) );
                localModels.setZero();
                localValidity.resize( static_cast<Eigen::Index>( searchPixelCount ), 1 );
                localValidity.setZero();
                for( std::size_t region = 0; region < grids.size(); ++region )
                {
                    for( std::size_t search = 0; search < grids[region].searchPixelCount(); ++search )
                    {
                        const Eigen::Index globalSearchIndex =
                            static_cast<Eigen::Index>( regionOffsets[region] + search );
                        const Eigen::Index sourceColumn =
                            static_cast<Eigen::Index>( search * m_modeFractions.size() + output );
                        localModels.col( globalSearchIndex ) = m_localPSFModels[region].col( sourceColumn );
                        if( m_localPSFTemporalCoefficients[region].rows() != 0 )
                        {
                            temporalCoefficients.col( globalSearchIndex )
                                .head( m_localPSFTemporalCoefficients[region].rows() ) =
                                m_localPSFTemporalCoefficients[region].col( sourceColumn );
                        }
                        localValidity( globalSearchIndex, 0 ) =
                            m_localPSFValidity[region]( static_cast<Eigen::Index>( search ),
                                                        static_cast<Eigen::Index>( output ) );
                    }
                }
            }

            eigenCube<realT> finalModels;
            imageT finalValidity;
            if( m_outputPSFModels )
            {
                finalModels.resize( m_psfStampSize, m_psfStampSize, static_cast<int>( searchPixelCount ) );
                finalValidity.resize( static_cast<Eigen::Index>( searchPixelCount ), 1 );
                finalValidity.setZero();
            }
            std::atomic<bool> reconstructionFailed{ false };
            std::exception_ptr reconstructionException;
            // clang-format off
#pragma omp parallel num_threads(effectiveWorkers)
            // clang-format on
            {
                imageT combined;
                reconstructorT::validityT combinedValidity;
                imageT centerCombined;
                reconstructorT::validityT centerValidity;

                // clang-format off
#pragma omp for schedule(static)
                // clang-format on
                for( std::size_t source = 0; source < searchPixelCount; ++source )
                {
                    if( reconstructionFailed.load( std::memory_order_acquire ) )
                    {
                        continue;
                    }
                    try
                    {
                        const double sourceRow =
                            static_cast<double>( coordinates( static_cast<Eigen::Index>( source ), 0 ) );
                        const double sourceColumn =
                            static_cast<double>( coordinates( static_cast<Eigen::Index>( source ), 1 ) );
                        if( targetHeldOutPSF )
                        {
                            reconstructor.reconstructCombinedTargeted( combined,
                                                                       combinedValidity,
                                                                       localModels,
                                                                       localValidity,
                                                                       searchIndex,
                                                                       sourceRow,
                                                                       sourceColumn,
                                                                       derotationAngles,
                                                                       this->m_combineMethod,
                                                                       this->m_comboWeights,
                                                                       this->m_sigmaThreshold,
                                                                       this->m_minGoodFract );
                        }
                        else
                        {
                            reconstructor.reconstructCombinedTemporal( combined,
                                                                       combinedValidity,
                                                                       localModels,
                                                                       temporalCoefficients,
                                                                       temporalOffsets,
                                                                       psfModel,
                                                                       localValidity,
                                                                       searchIndex,
                                                                       searchRegions,
                                                                       m_localPSFComponentCounts,
                                                                       m_temporalSelections,
                                                                       sourceRow,
                                                                       sourceColumn,
                                                                       derotationAngles,
                                                                       this->m_combineMethod,
                                                                       this->m_comboWeights,
                                                                       this->m_sigmaThreshold,
                                                                       this->m_minGoodFract );
                        }
                        const int imageRow = static_cast<int>( sourceRow );
                        const int imageColumn = static_cast<int>( sourceColumn );
                        if( m_psfFilter )
                        {
                            const P4PSFFilterResult filterResult =
                                P4PSFFilter::calculate( this->m_finim.image( output ),
                                                        combined,
                                                        combinedValidity,
                                                        imageRow,
                                                        imageColumn,
                                                        m_psfFilterMinGoodFract );
                            filterSupport.image( static_cast<int>( output ) )( imageRow, imageColumn ) =
                                static_cast<realT>( filterResult.supportFraction );
                            if( mx::math::isFinite( filterResult.normalization ) && filterResult.normalization >= 0 &&
                                filterResult.normalization <= std::numeric_limits<realT>::max() )
                            {
                                filterNormalization.image( static_cast<int>( output ) )( imageRow, imageColumn ) =
                                    static_cast<realT>( filterResult.normalization );
                            }
                            if( filterResult.valid && mx::math::isFinite( filterResult.amplitude ) &&
                                std::abs( filterResult.amplitude ) <= std::numeric_limits<realT>::max() )
                            {
                                filtered.image( static_cast<int>( output ) )( imageRow, imageColumn ) =
                                    static_cast<realT>( filterResult.amplitude );
                                filterValidity.image( static_cast<int>( output ) )( imageRow, imageColumn ) = 1;
                            }
                        }
                        if( m_outputPSFModels )
                        {
                            if( targetHeldOutPSF )
                            {
                                centerReconstructor.reconstructCombinedTargeted( centerCombined,
                                                                                 centerValidity,
                                                                                 localModels,
                                                                                 localValidity,
                                                                                 searchIndex,
                                                                                 sourceRow,
                                                                                 sourceColumn,
                                                                                 derotationAngles,
                                                                                 this->m_combineMethod,
                                                                                 this->m_comboWeights,
                                                                                 this->m_sigmaThreshold,
                                                                                 this->m_minGoodFract );
                            }
                            else
                            {
                                centerReconstructor.reconstructCombinedTemporal( centerCombined,
                                                                                 centerValidity,
                                                                                 localModels,
                                                                                 temporalCoefficients,
                                                                                 temporalOffsets,
                                                                                 psfModel,
                                                                                 localValidity,
                                                                                 searchIndex,
                                                                                 searchRegions,
                                                                                 m_localPSFComponentCounts,
                                                                                 m_temporalSelections,
                                                                                 sourceRow,
                                                                                 sourceColumn,
                                                                                 derotationAngles,
                                                                                 this->m_combineMethod,
                                                                                 this->m_comboWeights,
                                                                                 this->m_sigmaThreshold,
                                                                                 this->m_minGoodFract );
                            }
                            finalValidity( static_cast<Eigen::Index>( source ), 0 ) = centerValidity( 0, 0 );
                            for( int column = 0; column < combined.cols(); ++column )
                            {
                                for( int row = 0; row < combined.rows(); ++row )
                                {
                                    if( combinedValidity( row, column ) == 0 )
                                    {
                                        combined( row, column ) = invalidNumber<realT>();
                                    }
                                }
                            }
                            finalModels.image( static_cast<int>( source ) ) = combined;
                        }
                    }
                    catch( ... )
                    {
                        // clang-format off
#pragma omp critical(P4PSFReconstructionException)
                        // clang-format on
                        {
                            if( !reconstructionException )
                            {
                                reconstructionException = std::current_exception();
                            }
                        }
                        reconstructionFailed.store( true, std::memory_order_release );
                    }
                }
            }
            if( reconstructionException )
            {
                try
                {
                    std::rethrow_exception( reconstructionException );
                }
                catch( const std::exception &error )
                {
                    std::throw_with_nested( mx::exception<verboseT>( mx::error_t::exception,
                                                                     "P4 final PSF reconstruction failed for output " +
                                                                         std::to_string( output ) + ": " +
                                                                         error.what() ) );
                }
                catch( ... )
                {
                    std::throw_with_nested( mx::exception<verboseT>( mx::error_t::exception,
                                                                     "P4 final PSF reconstruction failed for output " +
                                                                         std::to_string( output ) ) );
                }
            }

            if( m_outputPSFModels )
            {
                fitsHeaderT modelHeader = productHeader( "MODEL", output );
                modelHeader.template append<std::string>( "P4 PSF PLANE ORDER",
                                                          "COORDINATES",
                                                          "plane mapping product" );
                writeProduct( productPath( "model_" + p4Index( output, 4 ) + ".fits" ), finalModels, modelHeader );
                fitsHeaderT validityHeader = productHeader( "VALIDITY", output );
                writeProduct( productPath( "validity_" + p4Index( output, 4 ) + ".fits" ),
                              finalValidity,
                              validityHeader );
            }
        }
    }

    if( m_psfFilter )
    {
        fitsHeaderT filteredHeader = productHeader( "FILTERED", m_modeFractions.size() );
        writeProduct( p4FilterProductPath( finalImagePath, "filtered", !this->m_exactFinimName ),
                      filtered,
                      filteredHeader );
        fitsHeaderT normalizationHeader = productHeader( "FILTER_NORMALIZATION", m_modeFractions.size() );
        writeProduct( p4FilterDiagnosticPath( finalImagePath, "filter_normalization", !this->m_exactFinimName ),
                      filterNormalization,
                      normalizationHeader );
        fitsHeaderT supportHeader = productHeader( "FILTER_SUPPORT", m_modeFractions.size() );
        writeProduct( p4FilterDiagnosticPath( finalImagePath, "filter_support", !this->m_exactFinimName ),
                      filterSupport,
                      supportHeader );
        fitsHeaderT filterValidityHeader = productHeader( "FILTER_VALIDITY", m_modeFractions.size() );
        writeProduct( p4FilterDiagnosticPath( finalImagePath, "filter_validity", !this->m_exactFinimName ),
                      filterValidity,
                      filterValidityHeader );
    }

    completion( 0, 0 ) = 1;
    fitsHeaderT completeHeader = productHeader( "MANIFEST", m_modeFractions.size() );
    completeHeader.template append<int>( "P4 PSF COMPLETE", 1, "complete product set available" );
    writeProduct( productPath( "manifest.fits" ), completion, completeHeader );
    m_psfReconstructionBytes = selectedReconstructionByteCount;
}

template <typename realT, class derotFunctObj, class verboseT>
void P4Reduction<realT, derotFunctObj, verboseT>::dump_times() const
{
    const double workerSeconds = m_timing.samplingWorkerSeconds + m_timing.gramWorkerSeconds +
                                 m_timing.eigensolveWorkerSeconds + m_timing.projectionWorkerSeconds +
                                 m_timing.psfWorkerSeconds;
    const auto percentage = [workerSeconds]( double seconds )
    { return workerSeconds > 0 ? seconds / workerSeconds * 100 : 0; };

    printf( "P4 reduction times: \n" );
    printf( "  Total time: %f sec\n", this->t_end - this->t_begin );
    printf( "    Loading: %f sec\n", this->t_load_end - this->t_load_begin );
    printf( "    Fake Injection: %f sec\n", this->t_fake_end - this->t_fake_begin );
    printf( "    Coadding: %f sec\n", this->t_coadd_end - this->t_coadd_begin );
    printf( "    Preprocessing: %f sec\n", this->t_preproc_end - this->t_preproc_begin );
    printf( "      Az USM: %f sec\n", this->t_azusm_end - this->t_azusm_begin );
    printf( "      Gauss USM: %f sec\n", this->t_gaussusm_end - this->t_gaussusm_begin );
    printf( "    P4 algorithm: %f elapsed real sec\n", m_timing.regressionElapsedSeconds );
    printf( "      Geometry %f elapsed real sec\n", m_timing.geometryElapsedSeconds );
    printf( "      Sampling %f worker sec (%f%%)\n",
            m_timing.samplingWorkerSeconds,
            percentage( m_timing.samplingWorkerSeconds ) );
    const auto samplingPercentage = [this]( double seconds )
    { return m_timing.samplingWorkerSeconds > 0 ? seconds / m_timing.samplingWorkerSeconds * 100 : 0; };
    printf( "        Same-image target/OR %f worker sec (%f%% of sampling)\n",
            m_timing.sameImageSamplingWorkerSeconds,
            samplingPercentage( m_timing.sameImageSamplingWorkerSeconds ) );
    printf( "        Additional-image PSF disk %f worker sec (%f%% of sampling)\n",
            m_timing.temporalSamplingWorkerSeconds,
            samplingPercentage( m_timing.temporalSamplingWorkerSeconds ) );
    printf( "      Gram construction %f worker sec (%f%%)\n",
            m_timing.gramWorkerSeconds,
            percentage( m_timing.gramWorkerSeconds ) );
    printf( "      EigenDecomposition %f worker sec (%f%%)\n",
            m_timing.eigensolveWorkerSeconds,
            percentage( m_timing.eigensolveWorkerSeconds ) );
    if( m_exclusionSolver == P4ExclusionSolver::factorDowndateExact )
    {
        const auto eigensolvePercentage = [this]( double seconds )
        { return m_timing.eigensolveWorkerSeconds > 0 ? seconds / m_timing.eigensolveWorkerSeconds * 100 : 0; };
        printf( "        Base factorization %f worker sec (%f%% of eigensolve)\n",
                m_timing.baseFactorWorkerSeconds,
                eigensolvePercentage( m_timing.baseFactorWorkerSeconds ) );
        printf( "        Row deletion %f worker sec (%f%% of eigensolve)\n",
                m_timing.deletionWorkerSeconds,
                eigensolvePercentage( m_timing.deletionWorkerSeconds ) );
        if( m_timing.explicitFallbackWorkerSeconds > 0 )
        {
            printf( "      Explicit fallback %f worker sec (overlaps phase totals)\n",
                    m_timing.explicitFallbackWorkerSeconds );
        }
    }
    printf( "      Projection/residual %f worker sec (%f%%)\n",
            m_timing.projectionWorkerSeconds,
            percentage( m_timing.projectionWorkerSeconds ) );
    if( !m_psfFile.empty() )
    {
        printf( "      Local PSF calculation %f worker sec (%f%%)\n",
                m_timing.psfWorkerSeconds,
                percentage( m_timing.psfWorkerSeconds ) );
    }
    if( m_outputPSFModels || m_psfFilter )
    {
        printf( "      PSF field reconstruction/filtering %f elapsed real sec\n",
                m_timing.psfReconstructionElapsedSeconds );
    }
    printf( "    Derotation: %f sec\n", this->t_derotate_end - this->t_derotate_begin );
    printf( "    Combination: %f sec\n", this->t_combo_end - this->t_combo_begin );
}

template <typename realT, class derotFunctObj, class verboseT>
template <typename dataT>
void P4Reduction<realT, derotFunctObj, verboseT>::writeDiagnostic( const std::string &fileName,
                                                                   const dataT &data,
                                                                   const fitsHeaderT *additionalHeader ) const
{
    if( !m_writeDiagnostics )
    {
        return;
    }

    std::string path;
    if( m_diagnosticDirectory.empty() || m_diagnosticDirectory == "." )
    {
        path = ( p4AuxiliaryProductDirectory( this->finalImageOutputPath() ) / fileName ).string();
    }
    else
    {
        path = m_diagnosticDirectory + "/" + fileName;
    }

    const std::string outputParent = mx::ioutils::parentPath( path );
    if( !outputParent.empty() )
    {
        const mx::error_t directoryResult = mx::ioutils::createDirectories( outputParent );
        if( directoryResult != mx::error_t::noerror )
        {
            throw mx::exception<verboseT>( directoryResult, "could not create P4 diagnostic output directory" );
        }
    }

    static std::atomic<unsigned long long> diagnosticSequence{ 0 };
    const std::string temporaryPath =
        path + ".tmp." + std::to_string( diagnosticSequence.fetch_add( 1, std::memory_order_relaxed ) );
    fitsHeaderT diagnosticHeader;
    appendReductionHeader( diagnosticHeader );
    if( additionalHeader )
    {
        fitsHeaderT additionalCopy = *additionalHeader;
        diagnosticHeader.append( additionalCopy );
    }
    mx::fits::fitsFile<typename dataT::Scalar, verboseT> writer;
    const mx::error_t writeResult = writer.write( temporaryPath, data, diagnosticHeader );
    if( writeResult != mx::error_t::noerror )
    {
        std::error_code cleanupError;
        std::filesystem::remove( temporaryPath, cleanupError );
        throw mx::exception<verboseT>( writeResult, "could not write temporary P4 diagnostic product " + path );
    }

    std::error_code renameError;
    std::filesystem::rename( temporaryPath, path, renameError );
    if( renameError )
    {
        std::error_code cleanupError;
        std::filesystem::remove( temporaryPath, cleanupError );
        throw mx::exception<verboseT>( mx::error_t::fileoerr,
                                       "could not publish P4 diagnostic product " + path + ": " +
                                           renameError.message() );
    }
}

template <typename realT, class derotFunctObj, class verboseT>
void P4Reduction<realT, derotFunctObj, verboseT>::appendReductionHeader( fitsHeaderT &head ) const
{
    head.append( "", fits::fitsCommentType(), "----------------------------------------" );
    head.append( "", fits::fitsCommentType(), "P4 reduction parameters:" );
    head.append( "", fits::fitsCommentType(), "----------------------------------------" );
    const bool factorDeletionActive =
        m_excludeMethod != HCI::exclude::none && m_exclusionSolver == P4ExclusionSolver::factorDowndateExact;
    head.template append<std::string>( "P4 ALGORITHM", "P4-PCA", "pixel prediction algorithm" );
    head.template append<std::string>( "P4POLCY",
                                       factorDeletionActive ? "D64-FACTOR-DOWNDATE" : "M32D64",
                                       "effective numerical precision policy" );
    head.template append<std::string>(
        "P4CALCPR", factorDeletionActive ? "FP64" : "FP32", "normal-equation and projection precision" );
    head.template append<std::string>( "P4EIGPR", "FP64", "symmetric eigensolve precision" );
    head.template append<std::string>(
        "P4RANKPR", factorDeletionActive ? "FP64" : "FP32", "eigenpair and rank-decision precision" );
    head.template append<std::string>( "P4OUTPR",
                                       std::is_same_v<realT, float> ? "FP32" : "FP64",
                                       "reduction product storage precision" );
    head.template append<std::string>(
        "P4FDELPR", factorDeletionActive ? "FP64" : "N/A", "exact factor-deletion precision" );
    head.template append<std::string>( "P4 FRAME", regressionFrameString( m_regressionFrame ), "regression frame" );
    head.template append<int>( "P4 IN SAMPLE",
                               m_excludeMethod == HCI::exclude::none ? 1 : 0,
                               "whether target rows enter their own temporal fits" );
    head.template append<std::string>( "P4 ADI EXCLUDE METHOD",
                                       HCI::excludeToStr<verboseT>( m_excludeMethod ),
                                       "target-frame exclusion method" );
    head.template append<std::string>( "P4 EXCLUSION SOLVER",
                                       exclusionSolverString( m_exclusionSolver ),
                                       "target-frame exclusion solver" );
    head.template append<std::string>( "P4 RANK MODEL",
                                       factorDeletionActive ? "completeBaseExact" : "direct",
                                       "held-out rank model" );
    head.template append<int>( "P4 FULL BASE",
                               factorDeletionActive ? 1 : 0,
                               "whether factor solver retained a complete safe base" );
    head.template append<std::string>( "P4 DELETION BACKEND",
                                       factorDeletionActive ? deletionBackendString( m_deletionBackend ) : "none",
                                       "mxlib SVD deletion backend" );
    head.template append<realT>( "P4 ADI MIN DPX", m_minDPx, "minimum target/reference displacement" );
    head.template append<int>( "P4 RDI", 0, "target-only ADI implementation" );
    head.template append<int>( "P4 NUMBER IMAGES", m_numberImages, "qualifying earlier and later predictor images" );
    head.template append<int>( "P4 LOCAL STAMP SIZE", m_localStampSize, "pixel-local sky result width; zero disables" );
    if( m_localStampSize > 0 )
    {
        head.template append<std::string>( "P4 LOCAL INJECTION STAGE",
                                           "POST_PREPROCESS_P4_INPUT",
                                           "trial injection stage" );
        head.template append<std::string>( "P4 LOCAL CENTER CONVENTION",
                                           "0.5*(N-1)",
                                           "full-image geometric center convention" );
        head.template append<std::string>( "P4 LOCAL LATTICE ANCHOR",
                                           "floor(source+0.5)",
                                           "integer output-stamp anchor convention" );
        head.template append<int>( "P4 LOCAL ORIGIN ROW", m_localOriginRow, "full-image row of local element 0,0" );
        head.template append<int>( "P4 LOCAL ORIGIN COLUMN",
                                   m_localOriginColumn,
                                   "full-image column of local element 0,0" );
        head.template append<double>( "P4 LOCAL SOURCE ROW", m_localSourceRow, "continuous full-image source row" );
        head.template append<double>( "P4 LOCAL SOURCE COLUMN",
                                      m_localSourceColumn,
                                      "continuous full-image source column" );
        head.template append<int>( "P4 LOCAL TEMPLATE ROWS", m_localTemplateRows, "phase-preserving template rows" );
        head.template append<int>( "P4 LOCAL TEMPLATE COLUMNS",
                                   m_localTemplateColumns,
                                   "phase-preserving template columns" );
        head.template append<std::string>( "P4 LOCAL SEARCH COUNT",
                                           std::to_string( m_localSearchCount ),
                                           "unique detector regressions" );
        head.template append<std::string>( "P4 LOCAL SPARSE SAMPLE COUNT",
                                           std::to_string( m_localSparseSampleCount ),
                                           "retained detector residual samples" );
        head.template append<std::string>( "P4 LOCAL GEOMETRY BYTES",
                                           std::to_string( m_localGeometryBytes ),
                                           "retained sparse geometry storage" );
    }
    head.template append<std::string>( "P4 MODE FRACTIONS",
                                       p4Join( m_modeFractions ),
                                       "ordered output mode fractions" );
    head.template append<std::string>( "P4 MIN RADIUS", p4Join( m_minRadius ), "search-annulus inner radii" );
    head.template append<std::string>( "P4 MAX RADIUS", p4Join( m_maxRadius ), "search-annulus outer radii" );
    head.template append<realT>( "P4 OR DELTA RADIUS INNER", m_orDeltaRadiusInner, "OR inward radial extent" );
    head.template append<realT>( "P4 OR DELTA RADIUS OUTER", m_orDeltaRadiusOuter, "OR outward radial extent" );
    head.template append<realT>( "P4 OR ARC HALF WIDTH", m_orArcHalfWidth, "OR arc half-width" );
    head.template append<realT>( "P4 OR MAX HALF ANGLE", m_orMaxHalfAngle, "OR maximum half-angle" );
    head.template append<realT>( "P4 PSF RADIUS", m_psfRadius, "physical PSF exclusion radius" );
    head.template append<std::string>( "P4 EXCLUSION POLICY",
                                       m_exclusionPolicy ? exclusionPolicyString( *m_exclusionPolicy ) : "invalid",
                                       "exclusion policy" );
    head.template append<realT>( "P4 EXCLUSION RADIUS BUFFER", m_exclusionRadiusBuffer, "exclusion-radius buffer" );
    head.template append<double>( "P4 RANK TOLERANCE", m_rankTolerance, "relative rank threshold" );
    head.template append<double>( "P4 MEMORY FRACTION", m_memoryFraction, "available-memory budget fraction" );
    head.template append<std::string>( "P4 AVAILABLE MEMORY BYTES",
                                       std::to_string( m_availableMemoryBytes ),
                                       "available-memory snapshot" );
    head.template append<std::string>( "P4 MEMORY BUDGET BYTES",
                                       std::to_string( m_memoryBudgetBytes ),
                                       "future-allocation budget" );
    head.template append<std::string>( "P4 COMPACT RESIDUAL BYTES",
                                       std::to_string( m_compactResidualBytes ),
                                       "compact residual estimate" );
    head.template append<std::string>( "P4 TARGET EXCLUSION BYTES",
                                       std::to_string( m_targetExclusionBytes ),
                                       "compact deleted-row storage" );
    head.template append<std::string>( "P4 MATERIALIZATION BYTES",
                                       std::to_string( m_materializationBytes ),
                                       "one output-pair estimate" );
    head.template append<int>( "P4 WRITE DIAGNOSTICS", m_writeDiagnostics ? 1 : 0, "diagnostics enabled" );

    std::vector<std::size_t> predictorCounts;
    std::vector<int> degreesOfFreedom;
    std::vector<int> minimumRanks;
    std::vector<int> minimumBaseRanks;
    std::vector<std::size_t> maximumExcludedRows;
    std::vector<std::size_t> exclusionStorageBytes;
    std::vector<std::size_t> downdateClampCounts;
    std::vector<std::size_t> explicitFallbackCounts;
    std::vector<std::size_t> rankBoundaryFallbackPixelCounts;
    std::vector<std::size_t> factorValidationFallbackPixelCounts;
    std::vector<std::size_t> deletionSolverFallbackPixelCounts;
    std::vector<double> maximumFactorOrthogonalityDefects;
    std::vector<double> factorOrthogonalityTolerancesAtMaximumDefect;
    std::vector<std::size_t> searchCounts;
    std::vector<std::size_t> targetImageCounts;
    std::vector<int> temporalImageCounts;
    std::vector<double> temporalPsfRadii;
    std::vector<std::size_t> workerBytes;
    std::vector<int> maximumWorkerCounts;
    std::vector<int> effectiveWorkerCounts;
    std::vector<std::size_t> validCounts;
    std::vector<std::size_t> maskedCounts;
    std::vector<std::size_t> supportInvalidCounts;
    predictorCounts.reserve( m_regionStatistics.size() );
    degreesOfFreedom.reserve( m_regionStatistics.size() );
    minimumRanks.reserve( m_regionStatistics.size() );
    minimumBaseRanks.reserve( m_regionStatistics.size() );
    maximumExcludedRows.reserve( m_regionStatistics.size() );
    exclusionStorageBytes.reserve( m_regionStatistics.size() );
    downdateClampCounts.reserve( m_regionStatistics.size() );
    explicitFallbackCounts.reserve( m_regionStatistics.size() );
    rankBoundaryFallbackPixelCounts.reserve( m_regionStatistics.size() );
    factorValidationFallbackPixelCounts.reserve( m_regionStatistics.size() );
    deletionSolverFallbackPixelCounts.reserve( m_regionStatistics.size() );
    maximumFactorOrthogonalityDefects.reserve( m_regionStatistics.size() );
    factorOrthogonalityTolerancesAtMaximumDefect.reserve( m_regionStatistics.size() );
    searchCounts.reserve( m_regionStatistics.size() );
    targetImageCounts.reserve( m_regionStatistics.size() );
    temporalImageCounts.reserve( m_regionStatistics.size() );
    temporalPsfRadii.reserve( m_regionStatistics.size() );
    workerBytes.reserve( m_regionStatistics.size() );
    maximumWorkerCounts.reserve( m_regionStatistics.size() );
    effectiveWorkerCounts.reserve( m_regionStatistics.size() );
    validCounts.reserve( m_regionStatistics.size() );
    maskedCounts.reserve( m_regionStatistics.size() );
    supportInvalidCounts.reserve( m_regionStatistics.size() );
    for( const P4RegionStatistics &statistics : m_regionStatistics )
    {
        predictorCounts.push_back( statistics.predictorCount );
        degreesOfFreedom.push_back( statistics.maximumDegreesOfFreedom );
        minimumRanks.push_back( statistics.minimumNumericalRank );
        minimumBaseRanks.push_back( statistics.minimumBaseRank );
        maximumExcludedRows.push_back( statistics.maximumExcludedRows );
        exclusionStorageBytes.push_back( statistics.exclusionStorageBytes );
        downdateClampCounts.push_back( statistics.downdateClampCount );
        explicitFallbackCounts.push_back( statistics.explicitFallbackCount );
        rankBoundaryFallbackPixelCounts.push_back( statistics.rankBoundaryFallbackPixelCount );
        factorValidationFallbackPixelCounts.push_back( statistics.factorValidationFallbackPixelCount );
        deletionSolverFallbackPixelCounts.push_back( statistics.deletionSolverFallbackPixelCount );
        maximumFactorOrthogonalityDefects.push_back( statistics.maximumFactorOrthogonalityDefect );
        factorOrthogonalityTolerancesAtMaximumDefect.push_back(
            statistics.factorOrthogonalityToleranceAtMaximumDefect );
        searchCounts.push_back( statistics.searchPixelCount );
        targetImageCounts.push_back( statistics.targetImageCount );
        temporalImageCounts.push_back( statistics.temporalNumberImages );
        temporalPsfRadii.push_back( statistics.temporalPsfRadius );
        workerBytes.push_back( statistics.estimatedWorkerBytes );
        maximumWorkerCounts.push_back( statistics.maximumWorkerCount );
        effectiveWorkerCounts.push_back( statistics.effectiveWorkerCount );
        validCounts.push_back( statistics.validLocalFitCount );
        maskedCounts.push_back( statistics.maskedLocalFitCount );
        supportInvalidCounts.push_back( statistics.supportInvalidLocalFitCount );
    }
    head.template append<std::string>( "P4 PREDICTOR COUNT", p4Join( predictorCounts ), "predictor counts by annulus" );
    head.template append<std::string>( "P4 MAX DOF", p4Join( degreesOfFreedom ), "maximum DOF by annulus" );
    head.template append<std::string>( "P4 MINIMUM RANK", p4Join( minimumRanks ), "minimum numerical rank by annulus" );
    head.template append<std::string>( "P4 BASE RANK", p4Join( minimumBaseRanks ), "minimum base rank by annulus" );
    head.template append<std::string>( "P4 MAX DELETED ROWS",
                                       p4Join( maximumExcludedRows ),
                                       "maximum target deletion count by annulus" );
    head.template append<std::string>( "P4 EXCLUSION STORAGE BYTES",
                                       p4Join( exclusionStorageBytes ),
                                       "compact exclusion bytes by annulus" );
    head.template append<std::string>( "P4 DELETION CLAMPS",
                                       p4Join( downdateClampCounts ),
                                       "roundoff clamps by annulus" );
    head.template append<std::string>( "P4 EXPLICIT FALLBACKS",
                                       p4Join( explicitFallbackCounts ),
                                       "recomputed target rows by annulus" );
    head.template append<std::string>( "P4 RANK BOUNDARY FALLBACK PIXELS",
                                       p4Join( rankBoundaryFallbackPixelCounts ),
                                       "rank-boundary fallback search pixels by annulus" );
    head.template append<std::string>( "P4 FACTOR VALIDATION FALLBACK PIXELS",
                                       p4Join( factorValidationFallbackPixelCounts ),
                                       "factor-validation fallback search pixels by annulus" );
    head.template append<std::string>( "P4 DELETION SOLVER FALLBACK PIXELS",
                                       p4Join( deletionSolverFallbackPixelCounts ),
                                       "deletion-solver fallback search pixels by annulus" );
    head.template append<std::string>( "P4 MAX FACTOR ORTHOGONALITY DEFECT",
                                       p4Join( maximumFactorOrthogonalityDefects ),
                                       "maximum rejected factor defect by annulus" );
    head.template append<std::string>( "P4 FACTOR TOLERANCE AT MAX DEFECT",
                                       p4Join( factorOrthogonalityTolerancesAtMaximumDefect ),
                                       "factor tolerance paired with maximum defect" );
    head.template append<std::string>( "P4 SEARCH PIXEL COUNT", p4Join( searchCounts ), "search pixels by annulus" );
    head.template append<std::string>( "P4 TARGET IMAGE COUNT",
                                       p4Join( targetImageCounts ),
                                       "retained target images by annulus" );
    head.template append<std::string>( "P4 TEMPORAL NUMBER IMAGES",
                                       p4Join( temporalImageCounts ),
                                       "effective temporal images by annulus" );
    head.template append<std::string>( "P4 TEMPORAL PSF RADIUS",
                                       p4Join( temporalPsfRadii ),
                                       "effective temporal PSF radii by annulus" );
    head.template append<std::string>( "P4 WORKER BYTES", p4Join( workerBytes ), "worker estimates by annulus" );
    head.template append<std::string>( "P4 MAXIMUM WORKERS",
                                       p4Join( maximumWorkerCounts ),
                                       "worker maxima before memory limiting" );
    head.template append<std::string>( "P4 EFFECTIVE WORKERS",
                                       p4Join( effectiveWorkerCounts ),
                                       "selected workers by annulus" );
    head.template append<std::string>( "P4 VALID FIT COUNT", p4Join( validCounts ), "valid local fits by annulus" );
    head.template append<std::string>( "P4 MASK INVALID FIT COUNT",
                                       p4Join( maskedCounts ),
                                       "detector-grid mask-invalid fits" );
    head.template append<std::string>( "P4 SUPPORT INVALID FIT COUNT",
                                       p4Join( supportInvalidCounts ),
                                       "direct support-invalid fits" );
    for( std::size_t region = 0; region < m_realizedModes.size(); ++region )
    {
        head.template append<std::string>( "P4 REALIZED MODES " + p4Index( region ),
                                           p4Join( m_realizedModes[region] ),
                                           "realized modes for one annulus" );
        head.template append<std::string>( "P4 RANK INVALID COUNT " + p4Index( region ),
                                           p4Join( m_regionStatistics[region].rankInvalidCounts ),
                                           "pixels with rank-invalid target fits by output plane" );
    }
}

template <typename realT, class derotFunctObj, class verboseT>
void P4Reduction<realT, derotFunctObj, verboseT>::cropAutomaticFinalImage()
{
    if( m_automaticOutputSize == 0 || this->m_finim.rows() == 0 || this->m_finim.cols() == 0 ||
        this->m_finim.planes() == 0 )
    {
        return;
    }
    if( m_automaticOutputSize > this->m_finim.rows() || m_automaticOutputSize > this->m_finim.cols() )
    {
        throw mx::exception<verboseT>( mx::error_t::sizeerr,
                                       "automatic P4 final-image crop exceeds the combined residual dimensions" );
    }

    const int firstRow = static_cast<int>(
        std::floor( 0.5 * ( this->m_finim.rows() - 1 ) - 0.5 * ( m_automaticOutputSize - 1 ) + 0.1 ) );
    const int firstColumn = static_cast<int>(
        std::floor( 0.5 * ( this->m_finim.cols() - 1 ) - 0.5 * ( m_automaticOutputSize - 1 ) + 0.1 ) );
    eigenCube<realT> cropped;
    cropped.resize( m_automaticOutputSize, m_automaticOutputSize, this->m_finim.planes() );
    for( int plane = 0; plane < this->m_finim.planes(); ++plane )
    {
        cropped.image( plane ) =
            this->m_finim.image( plane ).block( firstRow, firstColumn, m_automaticOutputSize, m_automaticOutputSize );
    }
    this->m_finim = std::move( cropped );
}

template <typename realT, class derotFunctObj, class verboseT>
void P4Reduction<realT, derotFunctObj, verboseT>::cropAutomaticFinalResiduals()
{
    if( m_automaticDerotationSize == 0 || m_psfFilter )
    {
        return;
    }

    const auto crop = [this]( auto &cube )
    {
        if( m_automaticDerotationSize > cube.rows() || m_automaticDerotationSize > cube.cols() )
        {
            throw mx::exception<verboseT>( mx::error_t::sizeerr,
                                           "automatic P4 residual crop exceeds the residual-cube dimensions" );
        }

        const int firstRow =
            static_cast<int>( std::floor( 0.5 * ( cube.rows() - 1 ) - 0.5 * ( m_automaticDerotationSize - 1 ) + 0.1 ) );
        const int firstColumn =
            static_cast<int>( std::floor( 0.5 * ( cube.cols() - 1 ) - 0.5 * ( m_automaticDerotationSize - 1 ) + 0.1 ) );
        std::remove_cvref_t<decltype( cube )> cropped;
        cropped.resize( m_automaticDerotationSize, m_automaticDerotationSize, cube.planes() );
        for( int plane = 0; plane < cube.planes(); ++plane )
        {
            cropped.image( plane ) = cube.image( plane ).block( firstRow,
                                                                firstColumn,
                                                                m_automaticDerotationSize,
                                                                m_automaticDerotationSize );
        }
        cube = std::move( cropped );
    };

    for( auto &cube : this->m_psfsub )
    {
        crop( cube );
    }
    for( auto &cube : this->m_psfsubValidity )
    {
        crop( cube );
    }
}

template <typename realT, class derotFunctObj, class verboseT>
int P4Reduction<realT, derotFunctObj, verboseT>::finalProcess( bool reportProgress )
{
    static_cast<void>( regressionFrameString( m_regressionFrame ) );
    if( m_regressionFrame == P4RegressionFrame::rotated && this->m_postMedSub )
    {
        throw mx::exception<verboseT>( mx::error_t::invalidconfig,
                                       "adi.postMedSub is not supported for rotated-frame P4 regression" );
    }

    fitsHeaderT algorithmHeader;
    const fitsHeaderT *algorithmHeaderPointer{ nullptr };
    if( this->m_doWriteFinim || this->m_doOutputPSFSub )
    {
        appendReductionHeader( algorithmHeader );
        algorithmHeaderPointer = &algorithmHeader;
    }
    const ADIDataFrame dataFrame =
        m_regressionFrame == P4RegressionFrame::rotated ? ADIDataFrame::sky : ADIDataFrame::detector;
    const std::string configuredPSFSubPrefix = this->m_PSFSubPrefix;
    const int configuredWriteFinim = this->m_doWriteFinim;
    if( m_automaticOutputSize != 0 )
    {
        this->m_doWriteFinim = false;
    }
    if( this->m_doOutputPSFSub )
    {
        const std::filesystem::path auxiliaryDirectory( auxiliaryOutputDirectory() );
        const std::filesystem::path relativeDirectory =
            this->m_outputDir.empty()
                ? auxiliaryDirectory
                : auxiliaryDirectory.lexically_relative( std::filesystem::path( this->m_outputDir ) );
        this->m_PSFSubPrefix = relativeDirectory.string();
        if( !this->m_PSFSubPrefix.empty() && !configuredPSFSubPrefix.empty() )
        {
            this->m_PSFSubPrefix += "/";
        }
        this->m_PSFSubPrefix += configuredPSFSubPrefix;
    }

    int result{ 0 };
    try
    {
        cropAutomaticFinalResiduals();
        result = this->ADIobservation<realT, derotFunctObj, verboseT>::finalProcess( algorithmHeaderPointer,
                                                                                     dataFrame,
                                                                                     reportProgress );
    }
    catch( ... )
    {
        this->m_PSFSubPrefix = configuredPSFSubPrefix;
        this->m_doWriteFinim = configuredWriteFinim;
        throw;
    }
    this->m_PSFSubPrefix = configuredPSFSubPrefix;
    this->m_doWriteFinim = configuredWriteFinim;
    return result;
}

#ifdef HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION
namespace detail
{

/** \cond P4Reduction_precision_experiment */
int p4ReductionReduceExperimental( P4Reductionf &reduction, P4PCAPrecisionPolicy precisionPolicy )
{
    switch( precisionPolicy )
    {
    case P4PCAPrecisionPolicy::doubleDouble:
    case P4PCAPrecisionPolicy::floatDouble:
    case P4PCAPrecisionPolicy::floatFloat:
        break;
    default:
        throw std::invalid_argument( "unknown experimental P4 precision policy" );
    }

    if( precisionPolicy == P4PCAPrecisionPolicy::floatFloat &&
        reduction.m_exclusionSolver == P4ExclusionSolver::factorDowndateExact )
    {
        throw std::invalid_argument( "experimental FP32 eigensolve does not support exact factor deletion" );
    }

    P4ReductionPrecisionSelectionScope selectionScope( precisionPolicy );
    return reduction.reduce();
}
/** \endcond */

} // namespace detail
#endif

template struct P4Reduction<float, ADIDerotator<float, verbose::vv>, verbose::vv>;

} // namespace improc
} // namespace mx
