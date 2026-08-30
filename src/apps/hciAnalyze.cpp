/** \file hciAnalyze.cpp
 * \brief Defines the hciAnalyze SNR-measurement command-line application.
 * \author Jared R. Males
 */

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <iostream>
#include <limits>
#include <numbers>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include <mx/app/application.hpp>
#include <mx/improc/eigenCube.hpp>
#include <mx/improc/imageFilters.hpp>
#include <mx/improc/imageMasks.hpp>
#include <mx/improc/imageUtils.hpp>
#include <mx/ioutils/fits/fitsFile.hpp>
#include <mx/ioutils/stringUtils.hpp>

#include "src/common/ConfigUtils.hpp"
#include "src/common/P4PSFFilter.hpp"

/// One configured or header-derived signal to measure.
struct hciAnalyzeSignal
{
    float m_separation{ 0 };       ///< Angular separation from image centre in pixels.
    float m_positionAngle{ 0 };    ///< Position angle in degrees east of north.
    float m_exclusionRadius{ 10 }; ///< Radius excluded from the noise estimate in pixels.
    float m_x{ 0 };                ///< Resolved FITS/display x coordinate, the first image index.
    float m_y{ 0 };                ///< Resolved FITS/display y coordinate, the second image index.
    float m_contrast{ std::numeric_limits<float>::quiet_NaN() }; ///< Optional fake-planet contrast.
};

/// A per-plane SNR measurement.
struct hciAnalyzeResult
{
    int m_plane{ 0 };                                            ///< Zero-based input cube plane index.
    float m_mode{ std::numeric_limits<float>::quiet_NaN() };     ///< Optional KLIP mode label from the FITS header.
    int m_signal{ 0 };                                           ///< Zero-based measured signal index.
    float m_separation{ 0 };                                     ///< Signal separation in pixels.
    float m_positionAngle{ 0 };                                  ///< Signal position angle in degrees east of north.
    float m_contrast{ std::numeric_limits<float>::quiet_NaN() }; ///< Optional fake-planet contrast.
    float m_snr{ std::numeric_limits<float>::quiet_NaN() };      ///< Maximum SNR in the configured aperture.
};

/// Command-line application for measuring signal-to-noise ratios in a PSF-subtracted image cube.
class hciAnalyze : public mx::app::application
{
  public:
    using realT = float;                                       ///< Scalar type used for FITS image calculations.
    using cubeT = mx::improc::eigenCube<realT>;                ///< FITS image-cube type.
    using imageT = mx::improc::eigenImage<realT>;              ///< FITS image type.
    using fitsHeaderT = mx::fits::fitsHeader<mx::verbose::vv>; ///< FITS header type.

  protected:
    std::string m_file;                       ///< Input PSF-subtracted FITS image or cube.

    std::vector<realT> m_planetSeparation;    ///< Explicit measured-signal separations in pixels.
    std::vector<realT> m_planetPositionAngle; ///< Explicit measured-signal PAs, degrees east of north.
    std::vector<realT> m_planetContrast;      ///< Explicit positive physical planet contrasts.
    std::vector<realT> m_planetRadius;        ///< Optional explicit exclusion radii in pixels.

    std::vector<realT> m_fakeSeparation;      ///< Explicit fake-planet separations in pixels.
    std::vector<realT> m_fakePositionAngle;   ///< Explicit fake-planet PAs, degrees east of north.
    std::vector<realT> m_fakeContrast;        ///< Optional explicit fake-planet contrasts.
    std::vector<realT> m_fakeRadius;          ///< Optional fake-planet exclusion radii in pixels.

    realT m_snrMinRadius{ 0 };                ///< Inner SNR annulus radius; non-positive derives from FITS metadata.
    realT m_snrMaxRadius{ 0 };                ///< Outer SNR annulus radius; non-positive derives from FITS metadata.
    realT m_snrApertureRadius{ 2 };           ///< Radius of each SNR aperture in pixels.
    realT m_lambdaD{ 2.5 };                   ///< Image scale in pixels per diffraction resolution element.
    bool m_lambdaDSpecified{ false };         ///< Whether lambdaD was supplied explicitly instead of defaulted.

    realT m_highPassFwhm{ 0 };                ///< Gaussian high-pass unsharp-mask FWHM in pixels.
    realT m_lowPassFwhm{ 0 };                 ///< Gaussian low-pass smoothing FWHM in pixels.
    std::string m_perPixelPSF;                ///< P4 PSF manifest supplying the spatially variable filter field.
    realT m_perPixelPSFMinimumSupport{ 1 };   ///< Minimum usable response-stamp fraction read from the manifest.
    bool m_diagnostics{ false };              ///< Whether to print resolved measurement diagnostics to standard error.

    bool m_planetSpecified{ false };          ///< True when any explicit planet target was supplied.
    bool m_fakeSpecified{ false };            ///< True when any explicit fake target was supplied.

    std::vector<hciAnalyzeSignal> m_signals;  ///< Resolved signals for the current input cube.
    std::vector<hciAnalyzeResult> m_results;  ///< Per-plane, per-signal SNR results from the current input cube.

  public:
    /// Construct hciAnalyze with its standard configuration search paths.
    hciAnalyze();

    /// Register hciAnalyze configuration targets.
    void setupConfig() override;

    /// Load hciAnalyze configuration values and record explicit coordinate sources.
    void loadConfig() override;

    /// Validate configuration values that do not depend on the input FITS header.
    void checkConfig() override;

    /// Read, filter, and measure the configured FITS image cube.
    int execute() override;

    /// Convert a separation and PA into FITS/display x and y coordinates.
    /** East is left in the displayed sky image: PA 0 degrees is up and PA +90 degrees is left.
     */
    static std::pair<realT, realT> signalCoordinates( realT separation,    /**< [in] separation in pixels */
                                                      realT positionAngle, /**< [in] PA in degrees east of north */
                                                      int xPixels,         /**< [in] number of first-index x pixels */
                                                      int yPixels /**< [in] number of second-index y pixels */ );

  protected:
    /// Return whether a configuration target was explicitly set.
    bool targetSpecified( const std::string &name /**< [in] registered configuration target name */ ) const;

    /// Validate a configured separation/PA/radius/contrast vector set.
    void validateSignalVectors( const std::vector<realT> &separation,    /**< [in] separations in pixels */
                                const std::vector<realT> &positionAngle, /**< [in] PAs in degrees */
                                const std::vector<realT> &radius,        /**< [in] optional exclusion radii */
                                const std::vector<realT> &contrast,      /**< [in] optional contrasts */
                                const std::string &source,               /**< [in] source description for errors */
                                bool requireContrast, /**< [in] whether a contrast vector is required */
                                bool nonnegativeContrast = false
                                /**< [in] whether every supplied contrast must be nonnegative */ ) const;

    /// Parse a comma-separated floating-point FITS keyword vector.
    std::vector<realT> headerVector( fitsHeaderT &header, /**< [in] input FITS header */
                                     const std::string &keyword /**< [in] FITS keyword name */ ) const;

    /// Resolve configured or header fake-planet metadata into image coordinates.
    void resolveSignals( fitsHeaderT &header, /**< [in] input FITS header */
                         int xPixels,         /**< [in] number of first-index x pixels */
                         int yPixels /**< [in] number of second-index y pixels */ );

    /// Determine the SNR annulus from configuration and FITS header metadata.
    std::pair<realT, realT> snrAnnulus( fitsHeaderT &header /**< [in] input FITS header */ ) const;

    /// Apply mask-aware Gaussian high-pass and low-pass filters to every cube plane.
    void filterCube( cubeT &cube,              /**< [in,out] image cube to filter */
                     const cubeT &invalidMask, /**< [in] non-zero marks samples excluded from every kernel */
                     realT highPassFwhm,       /**< [in] high-pass FWHM, non-positive disables */
                     realT lowPassFwhm /**< [in] low-pass FWHM, non-positive disables */ ) const;

    /// Apply an externally generated spatially variable P4 PSF field to every matching cube plane.
    void filterCubePerPixelPSF( cubeT &cube, /**< [in,out] science cube replaced by matched-filter amplitudes */
                                fitsHeaderT &scienceHeader /**< [in] science header supplying mode labels */ );

    /// Apply the Mawet et al. small-sample multiplicative correction to every SNR-cube pixel.
    void correctSmallSampleSNR( cubeT &snrCube /**< [in,out] uncorrected SNR cube */ ) const;

    /// Measure all resolved signals in an already-read FITS cube.
    void analyzeCube( cubeT &cube, /**< [in,out] input image cube, replaced with zeroed-invalid pixels */
                      fitsHeaderT &header /**< [in] input FITS header */ );

    /// Print stable machine-readable result rows.
    void printResults() const;

    /// Print resolved configuration, geometry, and per-plane measurement diagnostics.
    void printDiagnostics( realT minRadius, /**< [in] resolved inner SNR annulus radius */
                           realT maxRadius, /**< [in] resolved outer SNR annulus radius */
                           const cubeT &invalidMask /**< [in] invalid-pixel mask from the input cube */ ) const;

    /// Write or replace the full SNR cube next to the configured input FITS file.
    std::filesystem::path
    writeSNRMap( const cubeT &snrCube, /**< [in] calculated SNR image cube */
                 fitsHeaderT &header,  /**< [in,out] input FITS header augmented with SNR settings */
                 realT minRadius,      /**< [in] resolved inner SNR annulus radius */
                 realT maxRadius /**< [in] resolved outer SNR annulus radius */ ) const;
};

hciAnalyze::hciAnalyze()
{
    m_configPathGlobal_env = "HCIANALYZE_GLOBAL_CONFIG";
    m_configPathLocal = "hciAnalyze.conf";
    m_requireConfigPathLocal = false;
    config.m_sources = true;
}

void hciAnalyze::setupConfig()
{
    config
        .add( "file", "f", "file", mx::app::argType::Required, "", "file", true, "string", "input FITS image or cube" );
    config.add( "lambdaD",
                "",
                "lambdaD",
                mx::app::argType::Required,
                "",
                "lambdaD",
                false,
                "float",
                "image scale in pixels per lambda/D; default 2.5 with a warning when omitted" );

    config.add( "planet.sep",
                "",
                "planet.sep",
                mx::app::argType::Required,
                "planet",
                "sep",
                false,
                "vector<float>",
                "signal separations in pixels" );
    config.add( "planet.PA",
                "",
                "planet.PA",
                mx::app::argType::Required,
                "planet",
                "PA",
                false,
                "vector<float>",
                "signal position angles in degrees east of north; east is left" );
    config.add( "planet.contrast",
                "",
                "planet.contrast",
                mx::app::argType::Required,
                "planet",
                "contrast",
                false,
                "vector<float>",
                "positive physical planet contrasts" );
    config.add( "planet.R",
                "",
                "planet.R",
                mx::app::argType::Required,
                "planet",
                "R",
                false,
                "vector<float>",
                "signal exclusion-mask radii in pixels" );

    config.add( "fake.sep",
                "",
                "fake.sep",
                mx::app::argType::Required,
                "fake",
                "sep",
                false,
                "vector<float>",
                "fake-planet separations in pixels" );
    config.add( "fake.PA",
                "",
                "fake.PA",
                mx::app::argType::Required,
                "fake",
                "PA",
                false,
                "vector<float>",
                "fake-planet position angles in degrees east of north; east is left" );
    config.add( "fake.contrast",
                "",
                "fake.contrast",
                mx::app::argType::Required,
                "fake",
                "contrast",
                false,
                "vector<float>",
                "fake-planet contrasts" );
    config.add( "fake.R",
                "",
                "fake.R",
                mx::app::argType::Required,
                "fake",
                "R",
                false,
                "vector<float>",
                "fake-planet exclusion-mask radii in pixels" );

    config.add( "snr.minRad",
                "",
                "snr.minRad",
                mx::app::argType::Required,
                "snr",
                "minRad",
                false,
                "float",
                "inner SNR annulus radius in pixels; non-positive uses REGMINR or P4 MIN RADIUS" );
    config.add( "snr.maxRad",
                "",
                "snr.maxRad",
                mx::app::argType::Required,
                "snr",
                "maxRad",
                false,
                "float",
                "outer SNR annulus radius in pixels; non-positive uses REGMAXR or P4 MAX RADIUS" );
    config.add( "snr.apertureR",
                "",
                "snr.apertureR",
                mx::app::argType::Required,
                "snr",
                "apertureR",
                false,
                "float",
                "SNR aperture radius in pixels" );

    config.add( "filter.hpfGaussFW",
                "",
                "filter.hpfGaussFW",
                mx::app::argType::Required,
                "filter",
                "hpfGaussFW",
                false,
                "float",
                "Gaussian high-pass unsharp-mask FWHM in pixels; non-positive disables" );
    config.add( "filter.lpfGaussFW",
                "",
                "filter.lpfGaussFW",
                mx::app::argType::Required,
                "filter",
                "lpfGaussFW",
                false,
                "float",
                "Gaussian low-pass smoothing FWHM in pixels; non-positive disables" );
    config.add( "filter.perPixelPSF",
                "",
                "filter.perPixelPSF",
                mx::app::argType::Required,
                "filter",
                "perPixelPSF",
                false,
                "string",
                "complete P4 per-pixel PSF manifest used for spatially variable matched filtering" );
    config.add( "diagnostics",
                "d",
                "diagnostics",
                mx::app::argType::Optional,
                "",
                "diagnostics",
                false,
                "bool",
                "print resolved configuration, geometry, and measurement diagnostics to standard error" );
}

void hciAnalyze::loadConfig()
{
    config( m_file, "file" );
    config( m_lambdaD, "lambdaD" );
    config( m_planetSeparation, "planet.sep" );
    config( m_planetPositionAngle, "planet.PA" );
    config( m_planetContrast, "planet.contrast" );
    config( m_planetRadius, "planet.R" );
    config( m_fakeSeparation, "fake.sep" );
    config( m_fakePositionAngle, "fake.PA" );
    config( m_fakeContrast, "fake.contrast" );
    config( m_fakeRadius, "fake.R" );
    config( m_snrMinRadius, "snr.minRad" );
    config( m_snrMaxRadius, "snr.maxRad" );
    config( m_snrApertureRadius, "snr.apertureR" );
    config( m_highPassFwhm, "filter.hpfGaussFW" );
    config( m_lowPassFwhm, "filter.lpfGaussFW" );
    config( m_perPixelPSF, "filter.perPixelPSF" );
    mx::improc::loadBoolConfig<mx::verbose::vv>( config, m_diagnostics, "diagnostics" );

    m_planetSpecified = targetSpecified( "planet.sep" ) || targetSpecified( "planet.PA" ) ||
                        targetSpecified( "planet.contrast" ) || targetSpecified( "planet.R" );
    m_fakeSpecified = targetSpecified( "fake.sep" ) || targetSpecified( "fake.PA" ) ||
                      targetSpecified( "fake.contrast" ) || targetSpecified( "fake.R" );
    m_lambdaDSpecified = targetSpecified( "lambdaD" );
    if( !m_lambdaDSpecified )
    {
        std::cerr << "hciAnalyze warning: lambdaD was not set; using 2.5 pixels per lambda/D\n";
    }
}

void hciAnalyze::checkConfig()
{
    if( m_file.empty() )
    {
        throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig, "hciAnalyze requires file" );
    }
    if( m_planetSpecified && m_fakeSpecified )
    {
        throw mx::exception<mx::verbose::vv>(
            mx::error_t::invalidconfig,
            "specify either explicit planet coordinates or explicit fake metadata, not both" );
    }
    if( m_planetSpecified )
    {
        validateSignalVectors( m_planetSeparation,
                               m_planetPositionAngle,
                               m_planetRadius,
                               m_planetContrast,
                               "planet",
                               false,
                               true );
    }
    if( m_fakeSpecified )
    {
        validateSignalVectors( m_fakeSeparation, m_fakePositionAngle, m_fakeRadius, m_fakeContrast, "fake", false );
    }
    if( !std::isfinite( m_snrApertureRadius ) || m_snrApertureRadius <= 0 )
    {
        throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig, "snr.apertureR must be finite and positive" );
    }
    if( !std::isfinite( m_lambdaD ) || m_lambdaD <= 0 )
    {
        throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig, "lambdaD must be finite and positive" );
    }
    if( !std::isfinite( m_snrMinRadius ) || !std::isfinite( m_snrMaxRadius ) || !std::isfinite( m_highPassFwhm ) ||
        !std::isfinite( m_lowPassFwhm ) )
    {
        throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig, "SNR and filter radii must be finite" );
    }
    if( m_snrMinRadius > 0 && m_snrMaxRadius > 0 && m_snrMinRadius >= m_snrMaxRadius )
    {
        throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                              "snr.minRad must be smaller than snr.maxRad" );
    }
}

int hciAnalyze::execute()
{
    cubeT cube;
    fitsHeaderT header;
    mx::fits::fitsFile<realT, mx::verbose::vv> fitsFile;
    const mx::error_t result = fitsFile.read( cube, header, m_file );
    if( result != mx::error_t::noerror )
    {
        throw mx::exception<mx::verbose::vv>( result, "reading hciAnalyze input " + m_file );
    }
    if( cube.rows() == 0 || cube.cols() == 0 || cube.planes() == 0 )
    {
        throw mx::exception<mx::verbose::vv>( mx::error_t::sizeerr, "input FITS cube has no image planes" );
    }

    resolveSignals( header, cube.rows(), cube.cols() );
    analyzeCube( cube, header );
    printResults();
    return EXIT_SUCCESS;
}

std::pair<hciAnalyze::realT, hciAnalyze::realT>
hciAnalyze::signalCoordinates( realT separation, realT positionAngle, int xPixels, int yPixels )
{
    const realT angleRadians = positionAngle * std::numbers::pi_v<realT> / 180;
    const realT xCenter = 0.5F * static_cast<realT>( xPixels - 1 );
    const realT yCenter = 0.5F * static_cast<realT>( yPixels - 1 );
    return { xCenter - separation * std::sin( angleRadians ), yCenter + separation * std::cos( angleRadians ) };
}

bool hciAnalyze::targetSpecified( const std::string &name ) const
{
    return config.m_targets.at( name ).set;
}

void hciAnalyze::validateSignalVectors( const std::vector<realT> &separation,
                                        const std::vector<realT> &positionAngle,
                                        const std::vector<realT> &radius,
                                        const std::vector<realT> &contrast,
                                        const std::string &source,
                                        bool requireContrast,
                                        bool nonnegativeContrast ) const
{
    if( separation.empty() || positionAngle.empty() || separation.size() != positionAngle.size() )
    {
        throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                              source + ".sep and " + source + ".PA must be nonempty and equal-sized" );
    }
    if( !radius.empty() && radius.size() != 1 && radius.size() != separation.size() )
    {
        throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                              source + ".R must have one value or one value per signal" );
    }
    if( ( requireContrast && contrast.empty() ) || ( !contrast.empty() && contrast.size() != separation.size() ) )
    {
        throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                              source + ".contrast must have one value per signal" );
    }
    for( const realT value : separation )
    {
        if( !std::isfinite( value ) || value < 0 )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                  source + ".sep must be finite and nonnegative" );
        }
    }
    for( const realT value : positionAngle )
    {
        if( !std::isfinite( value ) )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig, source + ".PA must be finite" );
        }
    }
    for( const realT value : radius )
    {
        if( !std::isfinite( value ) || value < 0 )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                  source + ".R must be finite and nonnegative" );
        }
    }
    for( const realT value : contrast )
    {
        if( !std::isfinite( value ) || ( nonnegativeContrast && value < 0 ) )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                  source + ".contrast must be finite" +
                                                      ( nonnegativeContrast ? " and nonnegative" : "" ) );
        }
    }
}

std::vector<hciAnalyze::realT> hciAnalyze::headerVector( fitsHeaderT &header, const std::string &keyword ) const
{
    if( header.count( keyword ) == 0 )
    {
        return {};
    }
    std::vector<realT> values;
    mx::ioutils::parseStringVector( values, header[keyword].String() );
    return values;
}

void hciAnalyze::resolveSignals( fitsHeaderT &header, int xPixels, int yPixels )
{
    std::vector<realT> separation;
    std::vector<realT> positionAngle;
    std::vector<realT> radius;
    std::vector<realT> contrast;
    std::string source;
    bool requireContrast{ false };

    if( m_planetSpecified )
    {
        separation = m_planetSeparation;
        positionAngle = m_planetPositionAngle;
        radius = m_planetRadius;
        contrast = m_planetContrast;
        source = "planet";
    }
    else if( m_fakeSpecified )
    {
        separation = m_fakeSeparation;
        positionAngle = m_fakePositionAngle;
        radius = m_fakeRadius;
        contrast = m_fakeContrast;
        source = "fake";
    }
    else if( header.count( "PLANETSEP" ) != 0 || header.count( "PLANETPA" ) != 0 || header.count( "PLANETCONT" ) != 0 )
    {
        separation = headerVector( header, "PLANETSEP" );
        positionAngle = headerVector( header, "PLANETPA" );
        contrast = headerVector( header, "PLANETCONT" );
        source = "PLANET FITS keywords";
        requireContrast = true;
    }
    else
    {
        separation = headerVector( header, "FAKESEP" );
        positionAngle = headerVector( header, "FAKEPA" );
        contrast = headerVector( header, "FAKECONT" );
        source = "FAKE FITS keywords";
        requireContrast = true;
    }

    validateSignalVectors( separation,
                           positionAngle,
                           radius,
                           contrast,
                           source,
                           requireContrast,
                           source == "planet" || source == "PLANET FITS keywords" );
    m_signals.clear();
    m_signals.reserve( separation.size() );
    for( size_t index = 0; index < separation.size(); ++index )
    {
        const realT exclusionRadius = radius.empty() ? 10 : radius.size() == 1 ? radius.front() : radius[index];
        const auto [x, y] = signalCoordinates( separation[index], positionAngle[index], xPixels, yPixels );
        const bool apertureIntersects = x + m_snrApertureRadius >= 0 && x - m_snrApertureRadius <= xPixels - 1 &&
                                        y + m_snrApertureRadius >= 0 && y - m_snrApertureRadius <= yPixels - 1;
        const bool exclusionIntersects = x + exclusionRadius >= 0 && x - exclusionRadius <= xPixels - 1 &&
                                         y + exclusionRadius >= 0 && y - exclusionRadius <= yPixels - 1;
        if( !apertureIntersects || !exclusionIntersects )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                  "signal " + std::to_string( index ) +
                                                      " aperture or exclusion mask misses the image" );
        }

        hciAnalyzeSignal signal;
        signal.m_separation = separation[index];
        signal.m_positionAngle = positionAngle[index];
        signal.m_exclusionRadius = exclusionRadius;
        signal.m_x = x;
        signal.m_y = y;
        if( !contrast.empty() )
        {
            signal.m_contrast = contrast[index];
        }
        m_signals.push_back( signal );
    }
}

std::pair<hciAnalyze::realT, hciAnalyze::realT> hciAnalyze::snrAnnulus( fitsHeaderT &header ) const
{
    realT minRadius = m_snrMinRadius;
    realT maxRadius = m_snrMaxRadius;
    if( minRadius <= 0 )
    {
        const std::vector<realT> regionMinima = header.count( "REGMINR" ) ? headerVector( header, "REGMINR" )
                                                : header.count( "P4 MIN RADIUS" )
                                                    ? headerVector( header, "P4 MIN RADIUS" )
                                                    : headerVector( header, "P4MINR" );
        if( regionMinima.empty() )
        {
            throw mx::exception<mx::verbose::vv>(
                mx::error_t::invalidconfig,
                "snr.minRad is required when neither REGMINR nor P4 MIN RADIUS is in the FITS header" );
        }
        minRadius = *std::min_element( regionMinima.begin(), regionMinima.end() );
    }
    if( maxRadius <= 0 )
    {
        const std::vector<realT> regionMaxima = header.count( "REGMAXR" ) ? headerVector( header, "REGMAXR" )
                                                : header.count( "P4 MAX RADIUS" )
                                                    ? headerVector( header, "P4 MAX RADIUS" )
                                                    : headerVector( header, "P4MAXR" );
        if( regionMaxima.empty() )
        {
            throw mx::exception<mx::verbose::vv>(
                mx::error_t::invalidconfig,
                "snr.maxRad is required when neither REGMAXR nor P4 MAX RADIUS is in the FITS header" );
        }
        maxRadius = *std::max_element( regionMaxima.begin(), regionMaxima.end() );
    }
    if( !std::isfinite( minRadius ) || !std::isfinite( maxRadius ) || minRadius < 0 || maxRadius <= minRadius )
    {
        throw mx::exception<mx::verbose::vv>(
            mx::error_t::invalidconfig,
            "resolved SNR annulus must be finite, nonnegative, and have minRad < maxRad" );
    }
    return { minRadius, maxRadius };
}

void hciAnalyze::filterCube( cubeT &cube, const cubeT &invalidMask, realT highPassFwhm, realT lowPassFwhm ) const
{
    auto filterMasked = [&invalidMask]( imageT &filtered, const imageT &image, const auto &kernel, int plane )
    {
        typename std::remove_cvref_t<decltype( kernel )>::arrayT kernelArray;
        // gaussKernel::setKernel is position-independent and always succeeds.
        kernel.setKernel( 0, 0, kernelArray );

        const int rowRadius = ( kernelArray.rows() - 1 ) / 2;
        const int columnRadius = ( kernelArray.cols() - 1 ) / 2;
        filtered.resize( image.rows(), image.cols() );
        filtered.setZero();

        // clang-format off
        #pragma omp parallel for // clang-format on
        for( int row = 0; row < image.rows(); ++row )
        {
            for( int column = 0; column < image.cols(); ++column )
            {
                if( invalidMask.image( plane )( row, column ) != 0 )
                {
                    continue;
                }

                realT weightedSum{ 0 };
                realT normalization{ 0 };
                for( int kernelRow = 0; kernelRow < kernelArray.rows(); ++kernelRow )
                {
                    const int imageRow = row + kernelRow - rowRadius;
                    if( imageRow < 0 || imageRow >= image.rows() )
                    {
                        continue;
                    }
                    for( int kernelColumn = 0; kernelColumn < kernelArray.cols(); ++kernelColumn )
                    {
                        const int imageColumn = column + kernelColumn - columnRadius;
                        if( imageColumn < 0 || imageColumn >= image.cols() ||
                            invalidMask.image( plane )( imageRow, imageColumn ) != 0 )
                        {
                            continue;
                        }
                        const realT weight = kernelArray( kernelRow, kernelColumn );
                        weightedSum += weight * image( imageRow, imageColumn );
                        normalization += weight;
                    }
                }
                filtered( row, column ) = normalization > 0 ? weightedSum / normalization : 0;
            }
        }
    };

    for( int plane = 0; plane < cube.planes(); ++plane )
    {
        imageT filtered;
        if( highPassFwhm > 0 )
        {
            const imageT input = cube.image( plane );
            filterMasked( filtered, input, mx::improc::gaussKernel<imageT, 2>( highPassFwhm ), plane );
            cube.image( plane ) = input - filtered;
        }
        if( lowPassFwhm > 0 )
        {
            filterMasked( filtered, cube.image( plane ), mx::improc::gaussKernel<imageT, 4>( lowPassFwhm ), plane );
            cube.image( plane ) = filtered;
        }
    }
}

void hciAnalyze::filterCubePerPixelPSF( cubeT &cube, fitsHeaderT &scienceHeader )
{
    const std::filesystem::path manifestPath{ m_perPixelPSF };
    const std::string manifestName = manifestPath.filename().string();
    constexpr std::string_view manifestSuffix{ "manifest.fits" };
    if( !manifestName.ends_with( manifestSuffix ) )
    {
        throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                              "filter.perPixelPSF must name a P4 PSF manifest.fits product" );
    }
    const std::string productPrefix =
        ( manifestPath.parent_path() / manifestName.substr( 0, manifestName.size() - manifestSuffix.size() ) ).string();

    mx::fits::fitsFile<realT, mx::verbose::vv> reader;
    imageT manifest;
    fitsHeaderT manifestHeader;
    mx::error_t readResult = reader.read( manifest, manifestHeader, manifestPath.string() );
    if( readResult != mx::error_t::noerror )
    {
        throw mx::exception<mx::verbose::vv>( readResult, "reading per-pixel PSF manifest " + m_perPixelPSF );
    }
    const auto requireHeader = [&manifestHeader]( const std::string &keyword )
    {
        if( manifestHeader.count( keyword ) == 0 )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                  "per-pixel PSF manifest is missing " + keyword );
        }
    };
    for( const std::string &keyword : { "P4 PSF PRODUCT SCHEMA",
                                        "P4 PSF PRODUCT",
                                        "P4 PSF COMPLETE",
                                        "P4 PSF MODE COUNT",
                                        "P4 PSF SOURCE COUNT",
                                        "P4 PSF STAMP SIZE",
                                        "P4 PSF FILTER MIN GOOD FRACTION" } )
    {
        requireHeader( keyword );
    }
    if( manifest.rows() != 1 || manifest.cols() != 1 || manifest( 0, 0 ) != 1 ||
        !manifestHeader["P4 PSF PRODUCT"].String().starts_with( "MANIFEST" ) ||
        manifestHeader["P4 PSF COMPLETE"].value<int>() != 1 )
    {
        throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                              "filter.perPixelPSF does not name a complete P4 PSF manifest" );
    }
    const int schema = manifestHeader["P4 PSF PRODUCT SCHEMA"].value<int>();
    const int modeCount = manifestHeader["P4 PSF MODE COUNT"].value<int>();
    const int stampSize = manifestHeader["P4 PSF STAMP SIZE"].value<int>();
    m_perPixelPSFMinimumSupport = manifestHeader["P4 PSF FILTER MIN GOOD FRACTION"].value<realT>();
    if( schema < 1 || modeCount <= 0 || modeCount != cube.planes() || stampSize <= 0 || stampSize % 2 == 0 ||
        !std::isfinite( m_perPixelPSFMinimumSupport ) || m_perPixelPSFMinimumSupport < 0 ||
        m_perPixelPSFMinimumSupport > 1 )
    {
        throw mx::exception<mx::verbose::vv>(
            mx::error_t::invalidconfig,
            "per-pixel PSF schema, mode count, odd stamp size, or minimum support is incompatible with the input" );
    }

    const std::vector<realT> scienceModes =
        scienceHeader.count( "NMODES" )              ? headerVector( scienceHeader, "NMODES" )
        : scienceHeader.count( "FRACT NMODES" )      ? headerVector( scienceHeader, "FRACT NMODES" )
        : scienceHeader.count( "P4 MODE FRACTIONS" ) ? headerVector( scienceHeader, "P4 MODE FRACTIONS" )
                                                     : headerVector( scienceHeader, "P4MODFR" );
    const std::vector<realT> responseModes = manifestHeader.count( "P4 MODE FRACTIONS" )
                                                 ? headerVector( manifestHeader, "P4 MODE FRACTIONS" )
                                                 : headerVector( manifestHeader, "P4MODFR" );
    if( scienceModes.size() != static_cast<std::size_t>( modeCount ) || responseModes.size() != scienceModes.size() )
    {
        throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                              "science and per-pixel PSF products require matching mode labels" );
    }
    for( int mode = 0; mode < modeCount; ++mode )
    {
        const realT scale = std::max( { realT{ 1 }, std::abs( scienceModes[mode] ), std::abs( responseModes[mode] ) } );
        if( !std::isfinite( scienceModes[mode] ) || !std::isfinite( responseModes[mode] ) ||
            std::abs( scienceModes[mode] - responseModes[mode] ) > 8 * std::numeric_limits<realT>::epsilon() * scale )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                  "science and per-pixel PSF mode labels differ at plane " +
                                                      std::to_string( mode ) );
        }
    }

    imageT coordinates;
    fitsHeaderT coordinateHeader;
    const std::string coordinatePath = productPrefix + "coordinates.fits";
    readResult = reader.read( coordinates, coordinateHeader, coordinatePath );
    if( readResult != mx::error_t::noerror )
    {
        throw mx::exception<mx::verbose::vv>( readResult, "reading per-pixel PSF coordinates " + coordinatePath );
    }
    if( coordinates.rows() <= 0 || coordinates.cols() != 4 || coordinateHeader.count( "P4 PSF PRODUCT" ) == 0 ||
        !coordinateHeader["P4 PSF PRODUCT"].String().starts_with( "COORDINATES" ) )
    {
        throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                              "per-pixel PSF coordinates have an invalid role or shape" );
    }
    const std::string sourceCountString = manifestHeader["P4 PSF SOURCE COUNT"].String();
    std::size_t parsedCharacters{ 0 };
    std::size_t sourceCount{ 0 };
    try
    {
        sourceCount = std::stoull( sourceCountString, &parsedCharacters );
    }
    catch( const std::exception & )
    {
        throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                              "per-pixel PSF manifest source count is not an integer" );
    }
    const bool trailingCharactersAreWhitespace =
        std::all_of( sourceCountString.begin() + static_cast<std::ptrdiff_t>( parsedCharacters ),
                     sourceCountString.end(),
                     []( unsigned char character ) { return std::isspace( character ) != 0; } );
    if( !trailingCharactersAreWhitespace || sourceCount != static_cast<std::size_t>( coordinates.rows() ) )
    {
        throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                              "per-pixel PSF source count '" + sourceCountString + "' does not match " +
                                                  std::to_string( coordinates.rows() ) + " coordinate rows" );
    }

    std::vector<std::pair<int, int>> sourceCoordinates( sourceCount );
    std::vector<std::uint8_t> coordinateOwners( static_cast<std::size_t>( cube.rows() ) * cube.cols(), 0 );
    for( std::size_t source = 0; source < sourceCount; ++source )
    {
        const realT rowValue = coordinates( static_cast<Eigen::Index>( source ), 0 );
        const realT columnValue = coordinates( static_cast<Eigen::Index>( source ), 1 );
        if( !std::isfinite( rowValue ) || !std::isfinite( columnValue ) || rowValue != std::trunc( rowValue ) ||
            columnValue != std::trunc( columnValue ) || rowValue < 0 || rowValue >= cube.rows() || columnValue < 0 ||
            columnValue >= cube.cols() )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                  "per-pixel PSF source coordinate is not a valid image pixel" );
        }
        const int row = static_cast<int>( rowValue );
        const int column = static_cast<int>( columnValue );
        const std::size_t owner = static_cast<std::size_t>( row ) + static_cast<std::size_t>( cube.rows() ) * column;
        if( coordinateOwners[owner] != 0 )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                  "per-pixel PSF coordinate product contains duplicate pixels" );
        }
        coordinateOwners[owner] = 1;
        sourceCoordinates[source] = { row, column };
    }

    cubeT filtered( cube.rows(), cube.cols(), cube.planes() );
    filtered.cube().setConstant( std::numeric_limits<realT>::quiet_NaN() );
    for( int mode = 0; mode < modeCount; ++mode )
    {
        const std::string modeIndex = std::format( "{:04d}", mode );
        const std::string modelPath = productPrefix + "model_" + modeIndex + ".fits";
        const std::string validityPath = productPrefix + "validity_" + modeIndex + ".fits";
        cubeT models;
        imageT sourceValidity;
        fitsHeaderT modelHeader;
        fitsHeaderT validityHeader;
        readResult = reader.read( models, modelHeader, modelPath );
        if( readResult != mx::error_t::noerror )
        {
            throw mx::exception<mx::verbose::vv>( readResult, "reading per-pixel PSF model " + modelPath );
        }
        readResult = reader.read( sourceValidity, validityHeader, validityPath );
        if( readResult != mx::error_t::noerror )
        {
            throw mx::exception<mx::verbose::vv>( readResult, "reading per-pixel PSF validity " + validityPath );
        }
        if( models.rows() != stampSize || models.cols() != stampSize ||
            models.planes() != static_cast<int>( sourceCount ) || sourceValidity.rows() != coordinates.rows() ||
            sourceValidity.cols() != 1 || modelHeader.count( "P4 PSF PRODUCT" ) == 0 ||
            validityHeader.count( "P4 PSF PRODUCT" ) == 0 ||
            !modelHeader["P4 PSF PRODUCT"].String().starts_with( "MODEL" ) ||
            !validityHeader["P4 PSF PRODUCT"].String().starts_with( "VALIDITY" ) ||
            modelHeader.count( "P4 PSF MODE INDEX" ) == 0 || validityHeader.count( "P4 PSF MODE INDEX" ) == 0 ||
            modelHeader["P4 PSF MODE INDEX"].value<int>() != mode ||
            validityHeader["P4 PSF MODE INDEX"].value<int>() != mode )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                  "per-pixel PSF model or validity product is inconsistent" );
        }

        for( std::size_t source = 0; source < sourceCount; ++source )
        {
            if( !std::isfinite( sourceValidity( static_cast<Eigen::Index>( source ), 0 ) ) ||
                sourceValidity( static_cast<Eigen::Index>( source ), 0 ) <= 0 )
            {
                continue;
            }
            mx::improc::P4PSFFilter::validityT responseValidity( stampSize, stampSize );
            for( int column = 0; column < stampSize; ++column )
            {
                for( int row = 0; row < stampSize; ++row )
                {
                    responseValidity( row, column ) =
                        std::isfinite( models.image( static_cast<int>( source ) )( row, column ) ) ? 1 : 0;
                }
            }
            const auto [sourceRow, sourceColumn] = sourceCoordinates[source];
            const mx::improc::P4PSFFilterResult result =
                mx::improc::P4PSFFilter::calculate( cube.image( mode ),
                                                    models.image( static_cast<int>( source ) ),
                                                    responseValidity,
                                                    sourceRow,
                                                    sourceColumn,
                                                    m_perPixelPSFMinimumSupport );
            if( result.valid && std::isfinite( result.amplitude ) &&
                std::abs( result.amplitude ) <= std::numeric_limits<realT>::max() )
            {
                filtered.image( mode )( sourceRow, sourceColumn ) = static_cast<realT>( result.amplitude );
            }
        }
    }
    cube = std::move( filtered );
}

void hciAnalyze::correctSmallSampleSNR( cubeT &snrCube ) const
{
    const double centerRow = 0.5 * static_cast<double>( snrCube.rows() - 1 );
    const double centerColumn = 0.5 * static_cast<double>( snrCube.cols() - 1 );
    for( int column = 0; column < snrCube.cols(); ++column )
    {
        for( int row = 0; row < snrCube.rows(); ++row )
        {
            const double radiusPixels =
                std::hypot( static_cast<double>( row ) - centerRow, static_cast<double>( column ) - centerColumn );
            const double resolutionRadius = radiusPixels / static_cast<double>( m_lambdaD );
            const double comparisonSamples = 2.0 * std::numbers::pi * resolutionRadius - 1.0;
            const realT correction = comparisonSamples > 0
                                         ? static_cast<realT>( 1.0 / std::sqrt( 1.0 + 1.0 / comparisonSamples ) )
                                         : realT{ 0 };
            for( int plane = 0; plane < snrCube.planes(); ++plane )
            {
                snrCube.image( plane )( row, column ) *= correction;
            }
        }
    }
}

void hciAnalyze::analyzeCube( cubeT &cube, fitsHeaderT &header )
{
    const auto [minRadius, maxRadius] = snrAnnulus( header );
    if( !m_perPixelPSF.empty() )
    {
        filterCubePerPixelPSF( cube, header );
    }
    cubeT invalidMask;
    mx::improc::zeroNaNCube( cube, &invalidMask );
    filterCube( cube, invalidMask, m_highPassFwhm, m_lowPassFwhm );

    imageT signalMask( cube.rows(), cube.cols() );
    signalMask.setOnes();
    std::vector<imageT> apertureMasks;
    apertureMasks.reserve( m_signals.size() );
    for( const hciAnalyzeSignal &signal : m_signals )
    {
        mx::improc::maskCircle( signalMask, signal.m_x, signal.m_y, signal.m_exclusionRadius, 0.0F );
        imageT apertureMask( cube.rows(), cube.cols() );
        apertureMask.setZero();
        mx::improc::maskCircle( apertureMask, signal.m_x, signal.m_y, m_snrApertureRadius, 1.0F );
        apertureMasks.push_back( std::move( apertureMask ) );
    }

    cubeT noiseMask( cube.rows(), cube.cols(), cube.planes() );
    for( int plane = 0; plane < cube.planes(); ++plane )
    {
        noiseMask.image( plane ) = ( 1.0F - invalidMask.image( plane ) ) * signalMask;
    }

    cubeT snrCube;
    mx::improc::stddevImageCube( snrCube, cube, noiseMask, minRadius, maxRadius, true );
    correctSmallSampleSNR( snrCube );
    mx::improc::zeroNaNCube( snrCube );
    const std::filesystem::path snrMapPath = writeSNRMap( snrCube, header, minRadius, maxRadius );

    const std::vector<realT> modes = header.count( "NMODES" )              ? headerVector( header, "NMODES" )
                                     : header.count( "FRACT NMODES" )      ? headerVector( header, "FRACT NMODES" )
                                     : header.count( "P4 MODE FRACTIONS" ) ? headerVector( header, "P4 MODE FRACTIONS" )
                                                                           : headerVector( header, "P4MODFR" );
    m_results.clear();
    m_results.reserve( static_cast<size_t>( cube.planes() ) * m_signals.size() );
    for( int plane = 0; plane < cube.planes(); ++plane )
    {
        for( size_t signalIndex = 0; signalIndex < m_signals.size(); ++signalIndex )
        {
            realT apertureMaximum = -std::numeric_limits<realT>::infinity();
            bool hasValidAperturePixel{ false };
            for( int column = 0; column < cube.cols(); ++column )
            {
                for( int row = 0; row < cube.rows(); ++row )
                {
                    if( apertureMasks[signalIndex]( row, column ) == 0 ||
                        invalidMask.image( plane )( row, column ) != 0 )
                    {
                        continue;
                    }
                    apertureMaximum = std::max( apertureMaximum, snrCube.image( plane )( row, column ) );
                    hasValidAperturePixel = true;
                }
            }
            if( !hasValidAperturePixel )
            {
                throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                      "signal " + std::to_string( signalIndex ) +
                                                          " aperture contains no valid pixels" );
            }
            hciAnalyzeResult result;
            result.m_plane = plane;
            result.m_mode = modes.size() == static_cast<size_t>( cube.planes() )
                                ? modes[plane]
                                : std::numeric_limits<realT>::quiet_NaN();
            result.m_signal = static_cast<int>( signalIndex );
            result.m_separation = m_signals[signalIndex].m_separation;
            result.m_positionAngle = m_signals[signalIndex].m_positionAngle;
            result.m_contrast = m_signals[signalIndex].m_contrast;
            result.m_snr = apertureMaximum;
            m_results.push_back( result );
        }
    }
    if( m_diagnostics )
    {
        printDiagnostics( minRadius, maxRadius, invalidMask );
        std::cerr << "  SNR map: " << snrMapPath.string() << '\n';
    }
}

void hciAnalyze::printResults() const
{
    std::cout << "plane mode signal separation position_angle contrast snr\n";
    for( const hciAnalyzeResult &result : m_results )
    {
        std::cout << result.m_plane << ' ' << result.m_mode << ' ' << result.m_signal << ' ' << result.m_separation
                  << ' ' << result.m_positionAngle << ' ' << result.m_contrast << ' ' << result.m_snr << '\n';
    }
}

void hciAnalyze::printDiagnostics( realT minRadius, realT maxRadius, const cubeT &invalidMask ) const
{
    std::cerr << "hciAnalyze diagnostics:\n"
              << "  input: " << m_file << '\n'
              << "  SNR annulus: " << minRadius << " to " << maxRadius << " pixels\n"
              << "  aperture radius: " << m_snrApertureRadius << " pixels\n"
              << "  lambda/D scale: " << m_lambdaD << " pixels per lambda/D\n"
              << "  high-pass FWHM: " << m_highPassFwhm << " pixels\n"
              << "  low-pass FWHM: " << m_lowPassFwhm << " pixels\n"
              << "  per-pixel PSF manifest: " << ( m_perPixelPSF.empty() ? "disabled" : m_perPixelPSF ) << '\n';
    for( size_t signalIndex = 0; signalIndex < m_signals.size(); ++signalIndex )
    {
        const hciAnalyzeSignal &signal = m_signals[signalIndex];
        std::cerr << "  signal " << signalIndex << ": sep=" << signal.m_separation << " PA=" << signal.m_positionAngle
                  << " x=" << signal.m_x << " y=" << signal.m_y << " exclusion-radius=" << signal.m_exclusionRadius
                  << " contrast=" << signal.m_contrast << '\n';
    }
    for( int plane = 0; plane < invalidMask.planes(); ++plane )
    {
        const int invalidPixels = static_cast<int>( invalidMask.image( plane ).sum() );
        std::cerr << "  plane " << plane << ": invalid-pixels=" << invalidPixels << '\n';
    }
    for( const hciAnalyzeResult &result : m_results )
    {
        std::cerr << "  result: plane=" << result.m_plane << " mode=" << result.m_mode << " signal=" << result.m_signal
                  << " aperture-SNR=" << result.m_snr << '\n';
    }
}

std::filesystem::path
hciAnalyze::writeSNRMap( const cubeT &snrCube, fitsHeaderT &header, realT minRadius, realT maxRadius ) const
{
    for( const std::string &keyword : { "SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "NAXIS3", "EXTEND" } )
    {
        if( header.count( keyword ) == 0 )
        {
            continue;
        }
        const mx::error_t eraseResult = header.erase( keyword );
        if( eraseResult != mx::error_t::noerror )
        {
            throw mx::exception<mx::verbose::vv>( eraseResult, "removing structural FITS header " + keyword );
        }
    }
    const auto addHeader = [&header]( const std::string &keyword, const auto &value, const std::string &comment )
    {
        if( header.count( keyword ) != 0 )
        {
            const mx::error_t eraseResult = header.erase( keyword );
            if( eraseResult != mx::error_t::noerror )
            {
                throw mx::exception<mx::verbose::vv>( eraseResult, "replacing hciAnalyze SNR header " + keyword );
            }
        }
        const mx::error_t appendResult = header.append( keyword, value, comment );
        if( appendResult != mx::error_t::noerror )
        {
            throw mx::exception<mx::verbose::vv>( appendResult, "adding hciAnalyze SNR header " + keyword );
        }
    };
    const auto signalValues = [this]( const auto member )
    {
        std::string values;
        for( const hciAnalyzeSignal &signal : m_signals )
        {
            if( !values.empty() )
            {
                values += ',';
            }
            values += std::format( "{:.8g}", signal.*member );
        }
        return values;
    };

    addHeader( "HCIANAL", std::string{ "hciAnalyze" }, "SNR cube produced by" );
    addHeader( "SNRMINR", minRadius, "SNR annulus inner radius [pix]" );
    addHeader( "SNRMAXR", maxRadius, "SNR annulus outer radius [pix]" );
    addHeader( "SNRAPER", m_snrApertureRadius, "SNR aperture radius [pix]" );
    addHeader( "LAMBDAD", m_lambdaD, "image scale [pix per lambda/D]" );
    addHeader( "SNRMEAN", 1, "annular mean subtracted" );
    addHeader( "SNRSMALL", 1, "Mawet small-sample correction applied" );
    addHeader( "HPFGFW", m_highPassFwhm, "high-pass Gaussian FWHM [pix]" );
    addHeader( "LPFGFW", m_lowPassFwhm, "low-pass Gaussian FWHM [pix]" );
    addHeader( "HCIAPSF", m_perPixelPSF, "per-pixel PSF manifest" );
    addHeader( "HCIAPSM", m_perPixelPSFMinimumSupport, "per-pixel PSF minimum support" );
    addHeader( "HCISEPS", signalValues( &hciAnalyzeSignal::m_separation ), "signal separations [pix]" );
    addHeader( "HCIPAS", signalValues( &hciAnalyzeSignal::m_positionAngle ), "signal PAs [deg E of N]" );
    addHeader( "HCIRADS", signalValues( &hciAnalyzeSignal::m_exclusionRadius ), "signal exclusion radii [pix]" );
    addHeader( "HCIXS", signalValues( &hciAnalyzeSignal::m_x ), "signal first-index x coordinates [pix]" );
    addHeader( "HCIYS", signalValues( &hciAnalyzeSignal::m_y ), "signal second-index y coordinates [pix]" );

    const std::filesystem::path inputPath{ m_file };
    const std::filesystem::path outputPath =
        inputPath.parent_path() / ( inputPath.stem().string() + "_snr" + inputPath.extension().string() );
    mx::fits::fitsFile<realT, mx::verbose::vv> fitsFile;
    const mx::error_t result = fitsFile.write( outputPath.string(), snrCube, header );
    if( result != mx::error_t::noerror )
    {
        throw mx::exception<mx::verbose::vv>( result, "writing hciAnalyze SNR map " + outputPath.string() );
    }
    return outputPath;
}

/// Run hciAnalyze with standard exception reporting.
/** \returns `EXIT_SUCCESS` on help or successful execution; otherwise `EXIT_FAILURE`.
 */
int hciAnalyzeMain( int argc, /**< [in] standard command-line argument count */
                    char **argv /**< [in] standard command-line argument vector */ )
{
    const std::string argv0 = argc > 0 && argv != nullptr && argv[0] != nullptr ? argv[0] : "hciAnalyze";
    hciAnalyze analysis;
    try
    {
        const int result = analysis.main( argc, argv );
        return result == 1 ? EXIT_SUCCESS : result;
    }
    catch( const std::exception &exception )
    {
        std::vector<std::string> messages;
        mx::unwind_exceptions( messages, exception );
        mx::print_exceptions( messages, std::format( "{}: exception(s) encountered during execution", argv0 ) );
        std::cerr << std::format( "\nTo get help try: {} -h\n\n", argv0 );
        return EXIT_FAILURE;
    }
}

#ifndef HCIREDUCE_HCIANALYZE_NO_MAIN
/// hciAnalyze process entry point.
int main( int argc, /**< [in] standard command-line argument count */
          char **argv /**< [in] standard command-line argument vector */ )
{
    return hciAnalyzeMain( argc, argv );
}
#endif
