/** \file P4Reduction.cpp
 * \brief Implements the observation orchestrator for Pixel Prediction Post-Processing.
 * \author Jared R. Males
 */

#include "ADIDerotator.hpp"
#include "P4Reduction.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include <omp.h>

#include <mx/ioutils/fileUtils.hpp>
#include <mx/ioutils/fits/fitsFile.hpp>
#include <mx/ipc/ompLoopWatcher.hpp>
#include <mx/math/floatUtils.hpp>
#include <mx/sys/timeUtils.hpp>

namespace mx
{
namespace improc
{

namespace
{

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
                                                    double minimumRadius, /**< [in] annulus inner radius in pixels */
                                                    double psfRadius,     /**< [in] physical PSF radius in pixels */
                                                    int numberImages /**< [in] required qualifying images per side */ )
{
    std::vector<std::vector<int>> selections;
    for( std::size_t central = 0; central < derotationAngles.size(); ++central )
    {
        std::vector<int> earlier;
        std::vector<int> later;
        for( std::size_t candidate = central; candidate-- > 0; )
        {
            const double displacement =
                std::abs( mx::math::angleDiff<mx::math::radiansT<double>>( derotationAngles[candidate],
                                                                           derotationAngles[central] ) ) *
                minimumRadius;
            if( displacement >= psfRadius )
            {
                earlier.push_back( static_cast<int>( candidate ) );
            }
        }
        for( std::size_t candidate = central + 1; candidate < derotationAngles.size(); ++candidate )
        {
            const double displacement =
                std::abs( mx::math::angleDiff<mx::math::radiansT<double>>( derotationAngles[candidate],
                                                                           derotationAngles[central] ) ) *
                minimumRadius;
            if( displacement >= psfRadius )
            {
                later.push_back( static_cast<int>( candidate ) );
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
        std::vector<int> selection{ static_cast<int>( central ) };
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
                                 double meanRadius,              /**< [in] mean annulus radius in pixels */
                                 double requestedPsfRadius,      /**< [in] configured physical PSF radius in pixels */
                                 int numberImages,               /**< [in] requested qualifying images per side */
                                 std::size_t basePredictorCount, /**< [in] same-image OR predictor count */
                                 std::size_t temporalPredictorCount,
                                 /**< [in] physical-PSF predictor pixels per neighboring image */
                                 const std::vector<float> &modeFractions /**< [in] requested PCA fractions */ )
{
    std::vector<double> candidateRadii{ requestedPsfRadius };
    for( std::size_t first = 0; first < derotationAngles.size(); ++first )
    {
        for( std::size_t second = first + 1; second < derotationAngles.size(); ++second )
        {
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
            p4TemporalSelections( derotationAngles, meanRadius, radius, numberImages );
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
    config.add( "p4.writeDiagnostics",
                "",
                "p4.writeDiagnostics",
                mx::app::argType::True,
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
                "Destination directory for enabled P4 diagnostics" );
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

    config( m_numberImages, "p4.numberImages" );

    config( m_orDeltaRadiusInner, "p4.orDeltaRadiusInner" );
    config( m_orDeltaRadiusOuter, "p4.orDeltaRadiusOuter" );
    config( m_orArcHalfWidth, "p4.orArcHalfWidth" );
    config( m_orMaxHalfAngle, "p4.orMaxHalfAngle" );
    config( m_psfRadius, "p4.psfRadius" );

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
    config( m_writeDiagnostics, "p4.writeDiagnostics" );
    config( m_diagnosticDirectory, "p4.diagnosticDirectory" );
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
int P4Reduction<realT, derotFunctObj, verboseT>::regions( const std::vector<realT> &minimumRadii,
                                                          const std::vector<realT> &maximumRadii )
{
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

    this->t_begin = mx::sys::get_curr_time();
    this->m_psfsub.clear();
    this->m_psfsubValidity.clear();
    m_realizedModes.clear();
    m_regionStatistics.clear();
    m_temporalSelections.clear();
    m_ownership.resize( 0, 0 );
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

    const double geometryBegin = omp_get_wtime();
    std::vector<pixelGridT> grids( m_minRadius.size() );
    std::vector<rotatedGridT> rotatedGrids( m_minRadius.size() );
    m_realizedModes.resize( m_minRadius.size() );
    m_regionStatistics.resize( m_minRadius.size() );
    m_temporalSelections.resize( m_minRadius.size() );
    m_ownership.resize( this->m_Nrows, this->m_Ncols );
    m_ownership.setConstant( -1 );

    const imageT *mask = this->m_mask.size() == 0 ? nullptr : &this->m_mask;
    const std::vector<P4PixelCoordinate> temporalPredictorOffsets = p4TemporalPredictorOffsets( m_psfRadius );
    std::vector<double> derotationAngles;
    if( m_regressionFrame == P4RegressionFrame::rotated || m_numberImages > 0 )
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
                m_temporalSelections[region].reserve( static_cast<std::size_t>( this->m_Nims ) );
                for( int image = 0; image < this->m_Nims; ++image )
                {
                    m_temporalSelections[region].push_back( { image } );
                }
            }
            else
            {
                const double meanRadius =
                    0.5 * ( static_cast<double>( m_minRadius[region] ) + static_cast<double>( m_maxRadius[region] ) );
                auto selection = p4TemporalSelectionWithFallback( derotationAngles,
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
                    m_temporalSelections[region].reserve( static_cast<std::size_t>( this->m_Nims ) );
                    for( int image = 0; image < this->m_Nims; ++image )
                    {
                        m_temporalSelections[region].push_back( { image } );
                    }
                }
            }
        }
        else
        {
            m_temporalSelections[region].reserve( static_cast<std::size_t>( this->m_Nims ) );
            for( int image = 0; image < this->m_Nims; ++image )
            {
                m_temporalSelections[region].push_back( { image } );
            }
        }
        const std::size_t targetImageCount = m_temporalSelections[region].size();
        const std::size_t temporalImageCount = m_temporalSelections[region].front().size() - 1;
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

    std::size_t totalSearchPixels{ 0 };
    for( const P4RegionStatistics &statistics : m_regionStatistics )
    {
        totalSearchPixels += statistics.searchPixelCount;
    }

    this->m_psfsub.resize( m_modeFractions.size() );
    this->m_psfsubValidity.resize( m_modeFractions.size() );
    for( std::size_t output = 0; output < m_modeFractions.size(); ++output )
    {
        this->m_psfsub[output].resize( this->m_Nrows, this->m_Ncols, this->m_Nims );
        this->m_psfsub[output].setZero();
        this->m_psfsubValidity[output].resize( this->m_Nrows, this->m_Ncols, this->m_Nims );
        this->m_psfsubValidity[output].setZero();
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
        const bool usesTemporalPredictors = m_regionStatistics[region].temporalNumberImages > 0;
        const std::size_t predictorCount = m_regionStatistics[region].predictorCount;
        const std::vector<int> &modes = m_realizedModes[region];
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
        std::vector<std::size_t> rankInvalidCounts( modes.size(), 0 );

        // clang-format off
#pragma omp parallel
        // clang-format on
        {
            P4PCA::workspaceT workspace;
            P4PCA::matrixT predictors( static_cast<Eigen::Index>( targetImageCount ),
                                       static_cast<Eigen::Index>( predictorCount ) );
            P4PCA::vectorT target( static_cast<Eigen::Index>( targetImageCount ) );
            P4PCAResult result;
            std::size_t threadValid{ 0 };
            std::size_t threadMasked{ 0 };
            std::size_t threadSupportInvalid{ 0 };
            int threadMinimumRank{ std::numeric_limits<int>::max() };
            std::vector<std::size_t> threadRankInvalid( modes.size(), 0 );
            double threadSamplingSeconds{ 0 };
            double threadSameImageSamplingSeconds{ 0 };
            double threadTemporalSamplingSeconds{ 0 };
            double threadGramSeconds{ 0 };
            double threadEigensolveSeconds{ 0 };
            double threadProjectionSeconds{ 0 };

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
                                const std::size_t basePredictorCount = grid.predictorCount();
                                for( std::size_t targetIndex = 0; targetIndex < targetImageCount; ++targetIndex )
                                {
                                    const std::vector<int> &selection = temporalSelections[targetIndex];
                                    const int centralImage = selection[0];
                                    const double sameImageSamplingBegin = omp_get_wtime();
                                    target( static_cast<Eigen::Index>( targetIndex ) ) =
                                        static_cast<double>( this->m_tgtIms.image( centralImage )( row, column ) );
                                    for( std::size_t predictor = 0; predictor < basePredictorCount; ++predictor )
                                    {
                                        predictors( static_cast<Eigen::Index>( targetIndex ), predictor ) =
                                            checkedPredictorPromotion(
                                                grid.sample( this->m_tgtIms.image( centralImage ),
                                                             search,
                                                             predictor ) );
                                    }
                                    sameImageSamplingSeconds += omp_get_wtime() - sameImageSamplingBegin;

                                    const double temporalSamplingBegin = omp_get_wtime();
                                    for( std::size_t source = 1; source < selection.size(); ++source )
                                    {
                                        const int image = selection[source];
                                        for( std::size_t predictor = 0; predictor < temporalPredictorOffsets.size();
                                             ++predictor )
                                        {
                                            predictors( static_cast<Eigen::Index>( targetIndex ),
                                                        basePredictorCount +
                                                            ( source - 1 ) * temporalPredictorOffsets.size() +
                                                            predictor ) =
                                                checkedPredictorPromotion( this->m_tgtIms.image( image )(
                                                    row + temporalPredictorOffsets[predictor].row(),
                                                    column + temporalPredictorOffsets[predictor].column() ) );
                                        }
                                    }
                                    temporalSamplingSeconds += omp_get_wtime() - temporalSamplingBegin;
                                }
                            }

                            threadSamplingSeconds += sameImageSamplingSeconds + temporalSamplingSeconds;
                            threadSameImageSamplingSeconds += sameImageSamplingSeconds;
                            threadTemporalSamplingSeconds += temporalSamplingSeconds;
                            P4PCATiming pcaTiming;

                            if( rotated )
                            {
                                P4PCA::calculateCentered( result,
                                                          predictors,
                                                          target,
                                                          modes,
                                                          m_rankTolerance,
                                                          workspace,
                                                          &pcaTiming );
                            }
                            else
                            {
                                P4PCA::calculate( result,
                                                  predictors,
                                                  target,
                                                  modes,
                                                  m_rankTolerance,
                                                  workspace,
                                                  &pcaTiming );
                            }
                            threadGramSeconds += pcaTiming.gramWorkerSeconds;
                            threadEigensolveSeconds += pcaTiming.eigensolveWorkerSeconds;
                            threadProjectionSeconds += pcaTiming.projectionWorkerSeconds;
                            ++threadValid;
                            threadMinimumRank = std::min( threadMinimumRank, result.numericalRank );

                            const double residualApplyBegin = omp_get_wtime();
                            for( std::size_t output = 0; output < modes.size(); ++output )
                            {
                                if( result.modeStatus[output] == P4PCAModeStatus::rankInsufficient )
                                {
                                    ++threadRankInvalid[output];
                                    continue;
                                }
                                for( std::size_t targetIndex = 0; targetIndex < targetImageCount; ++targetIndex )
                                {
                                    const int image = temporalSelections[targetIndex][0];
                                    const realT residual = checkedResidualCast(
                                        result.residuals( static_cast<Eigen::Index>( targetIndex ), output ) );
                                    this->m_psfsub[output].image( image )( row, column ) = residual;
                                    this->m_psfsubValidity[output].image( image )( row, column ) = 1;
                                }
                            }
                            threadProjectionSeconds += omp_get_wtime() - residualApplyBegin;
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
                for( std::size_t output = 0; output < modes.size(); ++output )
                {
                    rankInvalidCounts[output] += threadRankInvalid[output];
                }
                m_timing.samplingWorkerSeconds += threadSamplingSeconds;
                m_timing.sameImageSamplingWorkerSeconds += threadSameImageSamplingSeconds;
                m_timing.temporalSamplingWorkerSeconds += threadTemporalSamplingSeconds;
                m_timing.gramWorkerSeconds += threadGramSeconds;
                m_timing.eigensolveWorkerSeconds += threadEigensolveSeconds;
                m_timing.projectionWorkerSeconds += threadProjectionSeconds;
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
        statistics.rankInvalidCounts = rankInvalidCounts;
        for( std::size_t output = 0; output < modes.size(); ++output )
        {
            if( rankInvalidCounts[output] != 0 )
            {
                std::cerr << "WARNING: P4 region " << region << " output " << output << " requested " << modes[output]
                          << " modes above numerical rank at " << rankInvalidCounts[output] << " search pixels\n";
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
            const std::size_t predictorCount = m_regionStatistics[region].predictorCount;
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
                        static_cast<Eigen::Index>( 12 + 2 * m_modeFractions.size() ) );
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
        }
        writeDiagnostic( "p4RegionSummary.fits", summary );

        imageT timing( 1, 8 );
        timing << static_cast<realT>( m_timing.geometryElapsedSeconds ),
            static_cast<realT>( m_timing.regressionElapsedSeconds ),
            static_cast<realT>( m_timing.samplingWorkerSeconds ),
            static_cast<realT>( m_timing.sameImageSamplingWorkerSeconds ),
            static_cast<realT>( m_timing.temporalSamplingWorkerSeconds ),
            static_cast<realT>( m_timing.gramWorkerSeconds ), static_cast<realT>( m_timing.eigensolveWorkerSeconds ),
            static_cast<realT>( m_timing.projectionWorkerSeconds );
        fitsHeaderT timingHeader;
        timingHeader.template append<int>( "P4 TIMING SCHEMA", 3, "P4 timing diagnostic schema version" );
        timingHeader.template append<std::string>( "P4 TIMING COLUMNS",
                                                   "geometryElapsed,regressionElapsed,samplingWorker,"
                                                   "sameImageSamplingWorker,temporalSamplingWorker,gramWorker,"
                                                   "eigensolveWorker,projectionWorker",
                                                   "P4 timing columns in seconds" );
        writeDiagnostic( "p4Timing.fits", timing, &timingHeader );
    }

    const int result = finalProcess();
    this->t_end = mx::sys::get_curr_time();
    return result;
}

template <typename realT, class derotFunctObj, class verboseT>
void P4Reduction<realT, derotFunctObj, verboseT>::dump_times() const
{
    const double workerSeconds = m_timing.samplingWorkerSeconds + m_timing.gramWorkerSeconds +
                                 m_timing.eigensolveWorkerSeconds + m_timing.projectionWorkerSeconds;
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
    printf( "      Projection/residual %f worker sec (%f%%)\n",
            m_timing.projectionWorkerSeconds,
            percentage( m_timing.projectionWorkerSeconds ) );
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

    std::string path = fileName;
    if( !m_diagnosticDirectory.empty() && m_diagnosticDirectory != "." )
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
    head.template append<std::string>( "P4 ALGORITHM", "P4-PCA", "pixel prediction algorithm" );
    head.template append<std::string>( "P4 FRAME", regressionFrameString( m_regressionFrame ), "regression frame" );
    head.template append<int>( "P4 IN SAMPLE", 1, "in-sample temporal regression" );
    head.template append<int>( "P4 RDI", 0, "target-only ADI implementation" );
    head.template append<int>( "P4 NUMBER IMAGES", m_numberImages, "qualifying earlier and later predictor images" );
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
    head.template append<int>( "P4 WRITE DIAGNOSTICS", m_writeDiagnostics ? 1 : 0, "diagnostics enabled" );

    std::vector<std::size_t> predictorCounts;
    std::vector<int> degreesOfFreedom;
    std::vector<int> minimumRanks;
    std::vector<std::size_t> searchCounts;
    std::vector<std::size_t> targetImageCounts;
    std::vector<int> temporalImageCounts;
    std::vector<double> temporalPsfRadii;
    std::vector<std::size_t> validCounts;
    std::vector<std::size_t> maskedCounts;
    std::vector<std::size_t> supportInvalidCounts;
    predictorCounts.reserve( m_regionStatistics.size() );
    degreesOfFreedom.reserve( m_regionStatistics.size() );
    minimumRanks.reserve( m_regionStatistics.size() );
    searchCounts.reserve( m_regionStatistics.size() );
    targetImageCounts.reserve( m_regionStatistics.size() );
    temporalImageCounts.reserve( m_regionStatistics.size() );
    temporalPsfRadii.reserve( m_regionStatistics.size() );
    validCounts.reserve( m_regionStatistics.size() );
    maskedCounts.reserve( m_regionStatistics.size() );
    supportInvalidCounts.reserve( m_regionStatistics.size() );
    for( const P4RegionStatistics &statistics : m_regionStatistics )
    {
        predictorCounts.push_back( statistics.predictorCount );
        degreesOfFreedom.push_back( statistics.maximumDegreesOfFreedom );
        minimumRanks.push_back( statistics.minimumNumericalRank );
        searchCounts.push_back( statistics.searchPixelCount );
        targetImageCounts.push_back( statistics.targetImageCount );
        temporalImageCounts.push_back( statistics.temporalNumberImages );
        temporalPsfRadii.push_back( statistics.temporalPsfRadius );
        validCounts.push_back( statistics.validLocalFitCount );
        maskedCounts.push_back( statistics.maskedLocalFitCount );
        supportInvalidCounts.push_back( statistics.supportInvalidLocalFitCount );
    }
    head.template append<std::string>( "P4 PREDICTOR COUNT", p4Join( predictorCounts ), "predictor counts by annulus" );
    head.template append<std::string>( "P4 MAX DOF", p4Join( degreesOfFreedom ), "maximum DOF by annulus" );
    head.template append<std::string>( "P4 MINIMUM RANK", p4Join( minimumRanks ), "minimum numerical rank by annulus" );
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
                                           "rank-invalid pixels by output plane" );
    }
}

template <typename realT, class derotFunctObj, class verboseT>
int P4Reduction<realT, derotFunctObj, verboseT>::finalProcess()
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
    return this->ADIobservation<realT, derotFunctObj, verboseT>::finalProcess( algorithmHeaderPointer, dataFrame );
}

template struct P4Reduction<float, ADIDerotator<float, verbose::vv>, verbose::vv>;

} // namespace improc
} // namespace mx
