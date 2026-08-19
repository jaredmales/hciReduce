/** \file P4Reduction.cpp
 * \brief Implements the observation orchestrator for Pixel Prediction Post-Processing.
 * \author Jared R. Males
 */

#include "ADIDerotator.hpp"
#include "P4Reduction.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <filesystem>
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
    m_ownership.resize( 0, 0 );
    m_geometrySeconds = 0;
    m_regressionSeconds = 0;

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
    m_ownership.resize( this->m_Nrows, this->m_Ncols );
    m_ownership.setConstant( -1 );

    const imageT *mask = this->m_mask.size() == 0 ? nullptr : &this->m_mask;
    std::vector<double> derotationAngles;
    if( m_regressionFrame == P4RegressionFrame::rotated )
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
        const std::size_t predictorCount = m_regressionFrame == P4RegressionFrame::detector
                                               ? grids[region].predictorCount()
                                               : rotatedGrids[region].predictorCount();
        std::cerr << "done: " << searchPixelCount << " search pixels, K=" << predictorCount << '\n';

        const int maximumDegreesOfFreedom =
            checkedMaximumDegreesOfFreedom( this->m_Nims,
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
        statistics.predictorCount = predictorCount;
        statistics.maximumDegreesOfFreedom = maximumDegreesOfFreedom;
        statistics.rankInvalidCounts.assign( m_modeFractions.size(), 0 );
    }
    m_geometrySeconds = omp_get_wtime() - geometryBegin;

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
            P4PCA::matrixT predictors( this->m_Nims, static_cast<Eigen::Index>( predictorCount ) );
            P4PCA::vectorT target( this->m_Nims );
            P4PCAResult result;
            std::size_t threadValid{ 0 };
            std::size_t threadMasked{ 0 };
            std::size_t threadSupportInvalid{ 0 };
            int threadMinimumRank{ std::numeric_limits<int>::max() };
            std::vector<std::size_t> threadRankInvalid( modes.size(), 0 );

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
                        if( !searchValid )
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
                            const P4PixelCoordinate &coordinate = rotated
                                                                      ? rotatedGrid.searchPixel( search ).coordinate()
                                                                      : grid.searchPixel( search ).coordinate();
                            const int row = coordinate.row();
                            const int column = coordinate.column();
                            if( rotated )
                            {
                                for( int image = 0; image < this->m_Nims; ++image )
                                {
                                    target( image ) = checkedPredictorPromotion(
                                        rotatedGrid.sampleTarget( this->m_tgtIms.image( image ),
                                                                  static_cast<std::size_t>( image ),
                                                                  search ) );
                                    for( std::size_t predictor = 0; predictor < predictorCount; ++predictor )
                                    {
                                        predictors( image, predictor ) = checkedPredictorPromotion(
                                            rotatedGrid.samplePredictor( this->m_tgtIms.image( image ),
                                                                         static_cast<std::size_t>( image ),
                                                                         search,
                                                                         predictor ) );
                                    }
                                }
                            }
                            else
                            {
                                for( int image = 0; image < this->m_Nims; ++image )
                                {
                                    target( image ) =
                                        static_cast<double>( this->m_tgtIms.image( image )( row, column ) );
                                    for( std::size_t predictor = 0; predictor < predictorCount; ++predictor )
                                    {
                                        predictors( image, predictor ) = checkedPredictorPromotion(
                                            grid.sample( this->m_tgtIms.image( image ), search, predictor ) );
                                    }
                                }
                            }

                            if( rotated )
                            {
                                P4PCA::calculateCentered( result,
                                                          predictors,
                                                          target,
                                                          modes,
                                                          m_rankTolerance,
                                                          workspace );
                            }
                            else
                            {
                                P4PCA::calculate( result, predictors, target, modes, m_rankTolerance, workspace );
                            }
                            ++threadValid;
                            threadMinimumRank = std::min( threadMinimumRank, result.numericalRank );

                            for( std::size_t output = 0; output < modes.size(); ++output )
                            {
                                if( result.modeStatus[output] == P4PCAModeStatus::rankInsufficient )
                                {
                                    ++threadRankInvalid[output];
                                    continue;
                                }
                                for( int image = 0; image < this->m_Nims; ++image )
                                {
                                    const realT residual = checkedResidualCast( result.residuals( image, output ) );
                                    this->m_psfsub[output].image( image )( row, column ) = residual;
                                    this->m_psfsubValidity[output].image( image )( row, column ) = 1;
                                }
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
                for( std::size_t output = 0; output < modes.size(); ++output )
                {
                    rankInvalidCounts[output] += threadRankInvalid[output];
                }
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
    m_regressionSeconds = omp_get_wtime() - regressionBegin;
    progressOutput.completeReduction( m_regressionSeconds );

    if( m_writeDiagnostics )
    {
        const imageT ownershipDiagnostic = m_ownership.template cast<realT>();
        writeDiagnostic( "p4Ownership.fits", ownershipDiagnostic );
        for( std::size_t output = 0; output < this->m_psfsubValidity.size(); ++output )
        {
            writeDiagnostic( "p4Validity_" + p4Index( output ) + ".fits", this->m_psfsubValidity[output] );
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
                        static_cast<Eigen::Index>( 9 + 2 * m_modeFractions.size() ) );
        for( std::size_t region = 0; region < m_regionStatistics.size(); ++region )
        {
            const P4RegionStatistics &statistics = m_regionStatistics[region];
            summary( region, 0 ) = m_minRadius[region];
            summary( region, 1 ) = m_maxRadius[region];
            summary( region, 2 ) = static_cast<double>( statistics.searchPixelCount );
            summary( region, 3 ) = static_cast<double>( statistics.predictorCount );
            summary( region, 4 ) = statistics.maximumDegreesOfFreedom;
            summary( region, 5 ) = statistics.minimumNumericalRank;
            summary( region, 6 ) = static_cast<realT>( statistics.validLocalFitCount );
            summary( region, 7 ) = static_cast<realT>( statistics.maskedLocalFitCount );
            summary( region, 8 ) = static_cast<realT>( statistics.supportInvalidLocalFitCount );
            for( std::size_t output = 0; output < m_modeFractions.size(); ++output )
            {
                summary( region, 9 + output ) = static_cast<realT>( m_realizedModes[region][output] );
                summary( region, 9 + m_modeFractions.size() + output ) =
                    static_cast<realT>( statistics.rankInvalidCounts[output] );
            }
        }
        writeDiagnostic( "p4RegionSummary.fits", summary );

        imageT timing( 1, 2 );
        timing << static_cast<realT>( m_geometrySeconds ), static_cast<realT>( m_regressionSeconds );
        writeDiagnostic( "p4Timing.fits", timing );
    }

    const int result = finalProcess();
    this->t_end = mx::sys::get_curr_time();
    return result;
}

template <typename realT, class derotFunctObj, class verboseT>
template <typename dataT>
void P4Reduction<realT, derotFunctObj, verboseT>::writeDiagnostic( const std::string &fileName,
                                                                   const dataT &data ) const
{
    if( !m_writeDiagnostics )
    {
        return;
    }

    std::string path = fileName;
    if( m_diagnosticDirectory != "." )
    {
        const mx::error_t directoryResult = mx::ioutils::createDirectories( m_diagnosticDirectory );
        if( directoryResult != mx::error_t::noerror )
        {
            throw mx::exception<verboseT>( directoryResult, "could not create P4 diagnostic directory" );
        }
        path = m_diagnosticDirectory + "/" + fileName;
    }

    static std::atomic<unsigned long long> diagnosticSequence{ 0 };
    const std::string temporaryPath =
        path + ".tmp." + std::to_string( diagnosticSequence.fetch_add( 1, std::memory_order_relaxed ) );
    fitsHeaderT diagnosticHeader;
    appendReductionHeader( diagnosticHeader );
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
    head.template append<std::string>( "P4ALGOR", "P4-PCA", "pixel prediction algorithm" );
    head.template append<std::string>( "P4 FRAME", regressionFrameString( m_regressionFrame ), "regression frame" );
    head.template append<int>( "P4INSAMP", 1, "in-sample temporal regression" );
    head.template append<int>( "P4RDI", 0, "target-only ADI implementation" );
    head.template append<std::string>( "P4MODFR", p4Join( m_modeFractions ), "ordered output mode fractions" );
    head.template append<std::string>( "P4MINR", p4Join( m_minRadius ), "search-annulus inner radii" );
    head.template append<std::string>( "P4MAXR", p4Join( m_maxRadius ), "search-annulus outer radii" );
    head.template append<realT>( "P4DRIN", m_orDeltaRadiusInner, "OR inward radial extent" );
    head.template append<realT>( "P4DROUT", m_orDeltaRadiusOuter, "OR outward radial extent" );
    head.template append<realT>( "P4ARCHW", m_orArcHalfWidth, "OR arc half-width" );
    head.template append<realT>( "P4MAXHA", m_orMaxHalfAngle, "OR maximum half-angle" );
    head.template append<realT>( "P4PSFR", m_psfRadius, "physical PSF exclusion radius" );
    head.template append<std::string>( "P4EXCL",
                                       m_exclusionPolicy ? exclusionPolicyString( *m_exclusionPolicy ) : "invalid",
                                       "exclusion policy" );
    head.template append<realT>( "P4EXBUF", m_exclusionRadiusBuffer, "exclusion-radius buffer" );
    head.template append<double>( "P4RKTOL", m_rankTolerance, "relative rank threshold" );
    head.template append<int>( "P4DIAG", m_writeDiagnostics ? 1 : 0, "diagnostics enabled" );

    std::vector<std::size_t> predictorCounts;
    std::vector<int> degreesOfFreedom;
    std::vector<int> minimumRanks;
    std::vector<std::size_t> searchCounts;
    std::vector<std::size_t> validCounts;
    std::vector<std::size_t> maskedCounts;
    std::vector<std::size_t> supportInvalidCounts;
    predictorCounts.reserve( m_regionStatistics.size() );
    degreesOfFreedom.reserve( m_regionStatistics.size() );
    minimumRanks.reserve( m_regionStatistics.size() );
    searchCounts.reserve( m_regionStatistics.size() );
    validCounts.reserve( m_regionStatistics.size() );
    maskedCounts.reserve( m_regionStatistics.size() );
    supportInvalidCounts.reserve( m_regionStatistics.size() );
    for( const P4RegionStatistics &statistics : m_regionStatistics )
    {
        predictorCounts.push_back( statistics.predictorCount );
        degreesOfFreedom.push_back( statistics.maximumDegreesOfFreedom );
        minimumRanks.push_back( statistics.minimumNumericalRank );
        searchCounts.push_back( statistics.searchPixelCount );
        validCounts.push_back( statistics.validLocalFitCount );
        maskedCounts.push_back( statistics.maskedLocalFitCount );
        supportInvalidCounts.push_back( statistics.supportInvalidLocalFitCount );
    }
    head.template append<std::string>( "P4K", p4Join( predictorCounts ), "predictor counts by annulus" );
    head.template append<std::string>( "P4DOF", p4Join( degreesOfFreedom ), "maximum DOF by annulus" );
    head.template append<std::string>( "P4RANK", p4Join( minimumRanks ), "minimum numerical rank by annulus" );
    head.template append<std::string>( "P4SRCH", p4Join( searchCounts ), "search pixels by annulus" );
    head.template append<std::string>( "P4VALID", p4Join( validCounts ), "valid local fits by annulus" );
    head.template append<std::string>( "P4MASK", p4Join( maskedCounts ), "detector-grid mask-invalid fits" );
    head.template append<std::string>( "P4SUP", p4Join( supportInvalidCounts ), "direct support-invalid fits" );
    for( std::size_t region = 0; region < m_realizedModes.size(); ++region )
    {
        head.template append<std::string>( "P4M" + p4Index( region ),
                                           p4Join( m_realizedModes[region] ),
                                           "realized modes for one annulus" );
        head.template append<std::string>( "P4I" + p4Index( region ),
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
