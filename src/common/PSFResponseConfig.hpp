/** \file PSFResponseConfig.hpp
 * \brief Declares shared sparse PSF-response configuration and radial-grid policy.
 */

#ifndef PSFResponseConfig_hpp
#define PSFResponseConfig_hpp

#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <mx/app/appConfigurator.hpp>

#include "ConfigUtils.hpp"

namespace mx
{
namespace improc
{

/// Measurement operator used to build a sparse radial PSF response model.
/** \ingroup programming_library */
enum class PSFResponseMethod : std::uint8_t
{
    skyExact,       ///< Propagate each sampled sky response through the fitted reduction operator.
    detectorLocal,  ///< Approximate each response with one detector-local frozen operator.
    refitDifference ///< Measure each response with paired finite-amplitude refits.
};

/// Resolved radial nodes and interpolation-region ownership for sparse PSF-response sampling.
/** \ingroup programming_library */
struct PSFResponseSamplingGrid
{
    std::vector<double> radii;        ///< Strictly increasing radii selected for response measurement.

    std::vector<std::size_t> regions; ///< Region index assigned one-to-one with `radii`, or all zero globally.

    bool regionAware{ false };        ///< Whether interpolation must remain within each reduction region.
};

/// Shared P4 and KLIP controls for sparse radial PSF-response estimation.
/** All controls are registered in the common `[psfResponse]` section. Algorithm implementations validate whether a
 * selected measurement method is available.
 *
 * \tparam realT response-configuration floating-point type
 * \ingroup programming_library
 */
template <typename realT>
struct PSFResponseConfig
{
    /** \name Sparse PSF Response Configuration - Data
     * @{
     */

    std::string m_psfFile;               ///< Optional centered PSF template enabling response measurement.

    int m_psfStampSize{ 0 };             ///< Positive square response-stamp size; sparse and filter paths require odd.

    std::vector<realT> m_psfSampleRadii; ///< Explicit strictly increasing response radii; empty allows region nodes.

    int m_psfRadiiPerRegion{ 0 };        ///< Interior radial nodes per reduction region; zero selects explicit radii.

    int m_psfSamplesPerRadius{ 0 };      ///< Uniform angular response measurements requested at every radial node.

    realT m_psfSampleArcStep{ 0 };       ///< Maximum azimuthal arc spacing; zero selects the fixed angular count.

    PSFResponseMethod m_psfSamplingMode{ PSFResponseMethod::skyExact }; ///< Sparse response measurement operator.

    realT m_psfSampleAvoidRadius{ 0 };  ///< Radius kept clear of configured known-source locations or trajectories.

    realT m_psfRefitContrast{ 0 };      ///< Positive half-amplitude for paired refit-difference measurements.

    bool m_outputPSFModels{ false };    ///< Whether canonical radial response and validity products are written.

    bool m_psfFilter{ false };          ///< Whether the spatially varying normalized PSF filter is applied.

    realT m_psfFilterMinGoodFract{ 1 }; ///< Minimum usable local-stamp fraction required by filtering.

    std::string m_psfOutputPrefix;      ///< Prefix for response products in the final image's auxiliary directory.

    /** @} */

    /** \name Sparse PSF Response Configuration
     * @{
     */

    /// Construct shared response controls with an algorithm-specific product prefix.
    explicit PSFResponseConfig( std::string outputPrefix /**< [in] default auxiliary-product prefix */ )
        : m_psfOutputPrefix( std::move( outputPrefix ) )
    {
    }

    /// Convert a supported response method to its stable configuration spelling.
    static std::string methodString( PSFResponseMethod method /**< [in] supported response method */ )
    {
        if( method == PSFResponseMethod::skyExact )
        {
            return "skyExact";
        }
        if( method == PSFResponseMethod::detectorLocal )
        {
            return "detectorLocal";
        }
        if( method == PSFResponseMethod::refitDifference )
        {
            return "refitDifference";
        }
        throw std::invalid_argument( "unsupported PSF response method" );
    }

    /// Parse an exact response-method configuration spelling.
    static PSFResponseMethod methodFromString( const std::string &value /**< [in] configuration spelling */ )
    {
        if( value == "skyExact" )
        {
            return PSFResponseMethod::skyExact;
        }
        if( value == "detectorLocal" )
        {
            return PSFResponseMethod::detectorLocal;
        }
        if( value == "refitDifference" )
        {
            return PSFResponseMethod::refitDifference;
        }
        throw std::invalid_argument( "unsupported PSF response method: " + value );
    }

    /// Register the common PSF-response configuration section.
    void setupPSFResponseConfig( mx::app::appConfigurator &config /**< [in,out] application configuration */ )
    {
        config.add( "psfResponse.file",
                    "",
                    "psfResponse.file",
                    mx::app::argType::Required,
                    "psfResponse",
                    "file",
                    false,
                    "string",
                    "Centered FITS PSF template used for sparse response measurement" );

        config.add( "psfResponse.stampSize",
                    "",
                    "psfResponse.stampSize",
                    mx::app::argType::Required,
                    "psfResponse",
                    "stampSize",
                    false,
                    "int",
                    "Positive square response-stamp size; sparse and filter paths require odd" );

        config.add( "psfResponse.sampleRadii",
                    "",
                    "psfResponse.sampleRadii",
                    mx::app::argType::Required,
                    "psfResponse",
                    "sampleRadii",
                    false,
                    "float vector",
                    "Explicit strictly increasing sparse response radii" );

        config.add( "psfResponse.radiiPerRegion",
                    "",
                    "psfResponse.radiiPerRegion",
                    mx::app::argType::Required,
                    "psfResponse",
                    "radiiPerRegion",
                    false,
                    "int",
                    "Interior radial response nodes per reduction region; mutually exclusive with explicit radii" );

        config.add( "psfResponse.samplesPerRadius",
                    "",
                    "psfResponse.samplesPerRadius",
                    mx::app::argType::Required,
                    "psfResponse",
                    "samplesPerRadius",
                    false,
                    "int",
                    "Uniform angular response measurement count at every radial node" );

        config.add( "psfResponse.sampleArcStep",
                    "",
                    "psfResponse.sampleArcStep",
                    mx::app::argType::Required,
                    "psfResponse",
                    "sampleArcStep",
                    false,
                    "float",
                    "Maximum azimuthal response-sample arc spacing; mutually exclusive with fixed count" );

        config.add( "psfResponse.method",
                    "",
                    "psfResponse.method",
                    mx::app::argType::Required,
                    "psfResponse",
                    "method",
                    false,
                    "string",
                    "Response method: skyExact, detectorLocal, or refitDifference" );

        config.add( "psfResponse.sampleAvoidRadius",
                    "",
                    "psfResponse.sampleAvoidRadius",
                    mx::app::argType::Required,
                    "psfResponse",
                    "sampleAvoidRadius",
                    false,
                    "float",
                    "Nonnegative distance kept clear of configured known sources" );

        config.add( "psfResponse.refitContrast",
                    "",
                    "psfResponse.refitContrast",
                    mx::app::argType::Required,
                    "psfResponse",
                    "refitContrast",
                    false,
                    "float",
                    "Positive half-amplitude for paired refit-difference response measurements" );

        config.add( "psfResponse.outputModels",
                    "",
                    "psfResponse.outputModels",
                    mx::app::argType::Optional,
                    "psfResponse",
                    "outputModels",
                    false,
                    "bool",
                    "Write canonical radial response and validity products" );

        config.add( "psfResponse.filter",
                    "",
                    "psfResponse.filter",
                    mx::app::argType::Optional,
                    "psfResponse",
                    "filter",
                    false,
                    "bool",
                    "Apply the spatially variable normalized PSF response filter" );

        config.add( "psfResponse.filterMinGoodFract",
                    "",
                    "psfResponse.filterMinGoodFract",
                    mx::app::argType::Required,
                    "psfResponse",
                    "filterMinGoodFract",
                    false,
                    "float",
                    "Minimum usable response-stamp fraction for filtering in [0,1]" );

        config.add( "psfResponse.outputPrefix",
                    "",
                    "psfResponse.outputPrefix",
                    mx::app::argType::Required,
                    "psfResponse",
                    "outputPrefix",
                    false,
                    "string",
                    "Prefix for sparse radial PSF response products" );
    }

    /// Load the common PSF-response configuration section.
    template <class verboseT>
    void loadPSFResponseConfig( mx::app::appConfigurator &config /**< [in] parsed application configuration */ )
    {
        config( m_psfFile, "psfResponse.file" );
        config( m_psfStampSize, "psfResponse.stampSize" );
        config( m_psfSampleRadii, "psfResponse.sampleRadii" );
        config( m_psfRadiiPerRegion, "psfResponse.radiiPerRegion" );
        config( m_psfSamplesPerRadius, "psfResponse.samplesPerRadius" );
        config( m_psfSampleArcStep, "psfResponse.sampleArcStep" );
        std::string method = methodString( m_psfSamplingMode );
        config( method, "psfResponse.method" );
        try
        {
            m_psfSamplingMode = methodFromString( method );
        }
        catch( ... )
        {
            std::throw_with_nested(
                mx::exception<verboseT>( mx::error_t::invalidconfig, "psfResponse.method is not supported" ) );
        }
        config( m_psfSampleAvoidRadius, "psfResponse.sampleAvoidRadius" );
        config( m_psfRefitContrast, "psfResponse.refitContrast" );
        loadBoolConfig<verboseT>( config, m_outputPSFModels, "psfResponse.outputModels" );
        loadBoolConfig<verboseT>( config, m_psfFilter, "psfResponse.filter" );
        config( m_psfFilterMinGoodFract, "psfResponse.filterMinGoodFract" );
        config( m_psfOutputPrefix, "psfResponse.outputPrefix" );
    }

    /// Report whether any non-default response control requests validation or measurement.
    bool psfResponseRequested() const
    {
        return !m_psfFile.empty() || m_psfStampSize != 0 || !m_psfSampleRadii.empty() || m_psfRadiiPerRegion != 0 ||
               m_psfSamplesPerRadius != 0 || m_psfSampleArcStep != 0 ||
               m_psfSamplingMode != PSFResponseMethod::skyExact || m_psfSampleAvoidRadius != 0 ||
               m_psfRefitContrast != 0 || m_outputPSFModels || m_psfFilter ||
               m_psfFilterMinGoodFract != static_cast<realT>( 1 );
    }

    /// Resolve explicit radii or uniformly spaced interior nodes within every reduction region.
    PSFResponseSamplingGrid resolvedPSFResponseGrid(
        const std::vector<realT> &minimumRadii, /**< [in] inclusive ordered region inner radii */
        const std::vector<realT> &maximumRadii /**< [in] exclusive ordered region outer radii */ ) const
    {
        PSFResponseSamplingGrid grid;
        if( !m_psfSampleRadii.empty() )
        {
            grid.radii.reserve( m_psfSampleRadii.size() );
            grid.regions.assign( m_psfSampleRadii.size(), 0 );
            for( const realT radius : m_psfSampleRadii )
            {
                grid.radii.push_back( static_cast<double>( radius ) );
            }
            return grid;
        }
        if( m_psfRadiiPerRegion <= 0 )
        {
            return grid;
        }
        if( minimumRadii.empty() || minimumRadii.size() != maximumRadii.size() )
        {
            throw std::invalid_argument( "PSF response regions must be nonempty and paired" );
        }

        const std::size_t nodesPerRegion = static_cast<std::size_t>( m_psfRadiiPerRegion );
        if( minimumRadii.size() > std::numeric_limits<std::size_t>::max() / nodesPerRegion )
        {
            throw std::length_error( "PSF response radial node count overflows size_t" );
        }
        grid.radii.reserve( minimumRadii.size() * nodesPerRegion );
        grid.regions.reserve( minimumRadii.size() * nodesPerRegion );
        grid.regionAware = true;
        for( std::size_t region = 0; region < minimumRadii.size(); ++region )
        {
            const double minimum = static_cast<double>( minimumRadii[region] );
            const double maximum = static_cast<double>( maximumRadii[region] );
            for( int node = 0; node < m_psfRadiiPerRegion; ++node )
            {
                const double fraction =
                    static_cast<double>( node + 1 ) / ( static_cast<double>( m_psfRadiiPerRegion ) + 1.0 );
                grid.radii.push_back( minimum + fraction * ( maximum - minimum ) );
                grid.regions.push_back( region );
            }
        }
        return grid;
    }

    /** @} */
};

} // namespace improc
} // namespace mx

#endif // PSFResponseConfig_hpp
