/** \file hciAnalyze.cpp
 * \brief Defines the hciAnalyze SNR-measurement command-line application.
 */

#include <algorithm>
#include <array>
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
#include <optional>
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
#include "src/common/PSFNoiseTraining.hpp"
#include "src/common/RadialPSFModel.hpp"

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
    std::string m_psfResponse;                ///< P4 or KLIP response manifest supplying the spatially variable filter.
    realT m_psfResponseMinimumSupport{ 1 };   ///< Minimum usable response-stamp fraction read from the manifest.
    std::string m_noiseModel{ "identity" };   ///< Residual-noise weighting: identity, diagonal, or pca.
    bool m_noiseDiagnostics{ false };         ///< Also write conditional diagnostic products for identity weighting.
    mx::improc::PSFNoiseTrainingConfig m_noiseConfig; ///< Annular geometry and explicit covariance regularization.
    std::vector<double> m_noiseExcludeRows;    ///< Additional source or held-out circle centers, first image index.
    std::vector<double> m_noiseExcludeColumns; ///< Additional source or held-out circle centers, second image index.
    std::vector<double> m_noiseExcludeRadii;   ///< Additional exclusion radii, one per supplied center.
    std::vector<mx::improc::PSFNoiseExclusion> m_noiseExclusions; ///< Resolved signals and additional held-out circles.
    std::array<cubeT, 11> m_noiseProducts;    ///< Amplitude, sigma, score, sample/mode counts, floor, status, support,
                                              ///< attempted/excluded/incomplete counts.

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

    /// Apply an externally generated P4 or KLIP PSF response field to every matching cube plane.
    void filterCubePSFResponse( cubeT &cube, /**< [in,out] science cube replaced by matched-filter amplitudes */
                                fitsHeaderT &scienceHeader /**< [in] science header supplying mode labels */ );

    /// Validate and resolve the optional covariance-weighting configuration.
    void checkNoiseConfig();

    /// Allocate conditional-statistic products and resolve all training exclusions.
    void prepareNoiseProducts( const cubeT &cube /**< [in] unfiltered science cube supplying output dimensions */ );

    /// Apply the selected noise model to one native response stamp and retain training diagnostics.
    mx::improc::P4PSFFilterResult
    applyPSFFilter( mx::improc::P4PSFFilter::imageConstRefT science,     /**< [in] unfiltered science image */
                    mx::improc::P4PSFFilter::imageConstRefT response,    /**< [in] native signed response stamp */
                    mx::improc::P4PSFFilter::validityConstRefT validity, /**< [in] valid response pixels */
                    int sourceRow,                                       /**< [in] candidate row */
                    int sourceColumn,                                    /**< [in] candidate column */
                    int mode /**< [in] output mode-plane index */ );

    /// Write covariance-filter diagnostics separately from the empirical annular SNR product.
    void writeNoiseProducts( const fitsHeaderT &scienceHeader /**< [in] original science metadata */ ) const;

    /// Reconstruct and apply one sparse radial KLIP PSF response field.
    void filterCubeKLIPPSFResponse( cubeT &cube, /**< [in,out] science cube replaced by matched-filter amplitudes */
                                    fitsHeaderT &scienceHeader,       /**< [in] science header supplying mode labels */
                                    const std::string &productPrefix, /**< [in] sibling-product path prefix */
                                    const imageT &manifest,           /**< [in] loaded completion image */
                                    fitsHeaderT &manifestHeader /**< [in] loaded sparse-response manifest header */ );

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
    config.add( "filter.psfResponse",
                "",
                "filter.psfResponse",
                mx::app::argType::Required,
                "filter",
                "psfResponse",
                false,
                "string",
                "complete P4 or KLIP PSF response manifest used for spatially variable matched filtering" );
    config.add( "noise.model",
                "",
                "noise.model",
                mx::app::argType::Required,
                "noise",
                "model",
                false,
                "string",
                "identity, diagonal, or pca residual-noise weighting; requires filter.psfResponse" );
    config.add( "noise.outputDiagnostics",
                "",
                "noise.outputDiagnostics",
                mx::app::argType::Optional,
                "noise",
                "outputDiagnostics",
                false,
                "bool",
                "also write conditional products for identity weighting" );
    config.add( "noise.arcStep",
                "",
                "noise.arcStep",
                mx::app::argType::Required,
                "noise",
                "arcStep",
                false,
                "double",
                "training-center arc spacing in pixels; zero uses half the response width" );
    config.add( "noise.guardRadius",
                "",
                "noise.guardRadius",
                mx::app::argType::Required,
                "noise",
                "guardRadius",
                false,
                "double",
                "additional candidate exclusion radius in pixels beyond the full response footprint" );
    config.add( "noise.minimumSamples",
                "",
                "noise.minimumSamples",
                mx::app::argType::Required,
                "noise",
                "minimumSamples",
                false,
                "size_t",
                "minimum complete noise patches; default 8; overlap does not imply independence" );
    config.add( "noise.maximumModes",
                "",
                "noise.maximumModes",
                mx::app::argType::Required,
                "noise",
                "maximumModes",
                false,
                "size_t",
                "maximum retained noise PCA modes; default 3; zero retains only the variance floor" );
    config.add( "noise.floorFraction",
                "",
                "noise.floorFraction",
                mx::app::argType::Required,
                "noise",
                "floorFraction",
                false,
                "double",
                "variance floor as a fraction of the median training-pixel variance in (0,1]; default 0.1" );
    config.add( "noise.excludeRows",
                "",
                "noise.excludeRows",
                mx::app::argType::Required,
                "noise",
                "excludeRows",
                false,
                "vector<double>",
                "additional source or held-out centers along the first image index" );
    config.add( "noise.excludeColumns",
                "",
                "noise.excludeColumns",
                mx::app::argType::Required,
                "noise",
                "excludeColumns",
                false,
                "vector<double>",
                "additional source or held-out centers along the second image index" );
    config.add( "noise.excludeRadii",
                "",
                "noise.excludeRadii",
                mx::app::argType::Required,
                "noise",
                "excludeRadii",
                false,
                "vector<double>",
                "one exclusion radius in pixels per additional center" );
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
    config( m_psfResponse, "filter.psfResponse" );
    config( m_noiseModel, "noise.model" );
    mx::improc::loadBoolConfig<mx::verbose::vv>( config, m_noiseDiagnostics, "noise.outputDiagnostics" );
    config( m_noiseConfig.m_arcStep, "noise.arcStep" );
    config( m_noiseConfig.m_guardRadius, "noise.guardRadius" );
    config( m_noiseConfig.m_minimumSamples, "noise.minimumSamples" );
    config( m_noiseConfig.m_maximumModes, "noise.maximumModes" );
    config( m_noiseConfig.m_floorFraction, "noise.floorFraction" );
    config( m_noiseExcludeRows, "noise.excludeRows" );
    config( m_noiseExcludeColumns, "noise.excludeColumns" );
    config( m_noiseExcludeRadii, "noise.excludeRadii" );
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
    checkNoiseConfig();
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

void hciAnalyze::checkNoiseConfig()
{
    if( m_noiseModel == "identity" )
    {
        m_noiseConfig.m_kind = mx::improc::PSFNoiseKind::identity;
    }
    else if( m_noiseModel == "diagonal" )
    {
        m_noiseConfig.m_kind = mx::improc::PSFNoiseKind::diagonal;
    }
    else if( m_noiseModel == "pca" )
    {
        m_noiseConfig.m_kind = mx::improc::PSFNoiseKind::pca;
    }
    else
    {
        throw std::invalid_argument( "noise.model must be identity, diagonal, or pca" );
    }
    if( m_noiseModel != "identity" && ( m_psfResponse.empty() || m_highPassFwhm > 0 || m_lowPassFwhm > 0 ) )
    {
        throw std::invalid_argument( "covariance weighting requires filter.psfResponse and disabled Gaussian filters" );
    }
    if( !std::isfinite( m_noiseConfig.m_arcStep ) || m_noiseConfig.m_arcStep < 0 ||
        !std::isfinite( m_noiseConfig.m_guardRadius ) || m_noiseConfig.m_guardRadius < 0 ||
        m_noiseConfig.m_minimumSamples < 2 || m_noiseConfig.m_minimumSamples > std::numeric_limits<int>::max() ||
        m_noiseConfig.m_maximumModes > std::numeric_limits<int>::max() ||
        !std::isfinite( m_noiseConfig.m_floorFraction ) || m_noiseConfig.m_floorFraction <= 0 ||
        m_noiseConfig.m_floorFraction > 1 )
    {
        throw std::invalid_argument( "noise settings require nonnegative finite geometry, at least two samples, "
                                     "nonnegative modes, and floorFraction in (0,1]" );
    }
    if( m_noiseExcludeRows.size() != m_noiseExcludeColumns.size() ||
        m_noiseExcludeRows.size() != m_noiseExcludeRadii.size() )
    {
        throw std::invalid_argument( "noise.excludeRows, excludeColumns, and excludeRadii must have equal lengths" );
    }
    for( std::size_t index = 0; index < m_noiseExcludeRows.size(); ++index )
    {
        if( !std::isfinite( m_noiseExcludeRows[index] ) || !std::isfinite( m_noiseExcludeColumns[index] ) ||
            !std::isfinite( m_noiseExcludeRadii[index] ) || m_noiseExcludeRadii[index] < 0 )
        {
            throw std::invalid_argument( "additional noise exclusions require finite centers and nonnegative radii" );
        }
    }
}

void hciAnalyze::prepareNoiseProducts( const cubeT &cube )
{
    checkNoiseConfig();
    if( m_noiseModel == "identity" && !m_noiseDiagnostics )
    {
        return;
    }
    m_noiseExclusions.clear();
    for( const auto &signal : m_signals )
    {
        m_noiseExclusions.push_back( { signal.m_x, signal.m_y, signal.m_exclusionRadius } );
    }
    for( std::size_t index = 0; index < m_noiseExcludeRows.size(); ++index )
    {
        m_noiseExclusions.push_back(
            { m_noiseExcludeRows[index], m_noiseExcludeColumns[index], m_noiseExcludeRadii[index] } );
    }
    for( auto &product : m_noiseProducts )
    {
        product.resize( cube.rows(), cube.cols(), cube.planes() );
        product.cube().setConstant( std::numeric_limits<realT>::quiet_NaN() );
    }
}

mx::improc::P4PSFFilterResult hciAnalyze::applyPSFFilter( mx::improc::P4PSFFilter::imageConstRefT science,
                                                          mx::improc::P4PSFFilter::imageConstRefT response,
                                                          mx::improc::P4PSFFilter::validityConstRefT validity,
                                                          int sourceRow,
                                                          int sourceColumn,
                                                          int mode )
{
    if( m_noiseModel == "identity" && !m_noiseDiagnostics )
    {
        return mx::improc::P4PSFFilter::calculate( science,
                                                   response,
                                                   validity,
                                                   sourceRow,
                                                   sourceColumn,
                                                   m_psfResponseMinimumSupport );
    }
    const mx::improc::PSFNoiseTrainingResult result =
        mx::improc::PSFNoiseTraining::calculate( science,
                                                 response,
                                                 validity,
                                                 sourceRow,
                                                 sourceColumn,
                                                 m_psfResponseMinimumSupport,
                                                 m_noiseConfig,
                                                 m_noiseExclusions );
    const double missing = std::numeric_limits<double>::quiet_NaN();
    const std::array<double, 11> values{ result.m_filter.valid ? result.m_filter.amplitude : missing,
                                         result.m_filter.valid ? result.m_filter.conditionalSigma : missing,
                                         result.m_filter.valid ? result.m_filter.score : missing,
                                         static_cast<double>( result.m_samples ),
                                         static_cast<double>( result.m_modes ),
                                         result.m_varianceFloor > 0 ? result.m_varianceFloor : missing,
                                         static_cast<double>( result.m_status ),
                                         result.m_status == mx::improc::PSFNoiseTrainingStatus::tooFewSamples ||
                                                 result.m_status == mx::improc::PSFNoiseTrainingStatus::zeroVariance
                                             ? missing
                                             : result.m_filter.supportFraction,
                                         static_cast<double>( result.m_attempted ),
                                         static_cast<double>( result.m_excluded ),
                                         static_cast<double>( result.m_incomplete ) };
    for( std::size_t index = 0; index < values.size(); ++index )
    {
        m_noiseProducts[index].image( mode )( sourceRow, sourceColumn ) =
            std::isfinite( values[index] ) && std::abs( values[index] ) <= std::numeric_limits<realT>::max()
                ? static_cast<realT>( values[index] )
                : std::numeric_limits<realT>::quiet_NaN();
    }
    return result.m_filter;
}

void hciAnalyze::writeNoiseProducts( const fitsHeaderT &scienceHeader ) const
{
    if( m_noiseModel == "identity" && !m_noiseDiagnostics )
    {
        return;
    }
    if( m_file.empty() )
    {
        throw std::invalid_argument( "noise diagnostic products require an input file path" );
    }
    fitsHeaderT header = scienceHeader;
    const auto add = [&header]( const std::string &keyword, const auto &value, const std::string &comment )
    {
        if( header.count( keyword ) && header.erase( keyword ) != mx::error_t::noerror )
        {
            throw std::runtime_error( "replacing noise diagnostic header " + keyword );
        }
        if( header.append( keyword, value, comment ) != mx::error_t::noerror )
        {
            throw std::runtime_error( "adding noise diagnostic header " + keyword );
        }
    };
    for( const std::string &keyword : { "SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "NAXIS3", "EXTEND" } )
    {
        if( header.count( keyword ) && header.erase( keyword ) != mx::error_t::noerror )
        {
            throw std::runtime_error( "removing structural noise-product header " + keyword );
        }
    }
    add( "HCIA NOISE SCHEMA", 1, "conditional covariance-filter diagnostics" );
    add( "HCIA NOISE MODEL", m_noiseModel, "residual-noise covariance representation" );
    add( "HCIA NOISE TRAINED", m_noiseModel == "identity" ? 0 : 1, "whether a noise covariance was estimated" );
    add( "HCIA NOISE ARC STEP", m_noiseConfig.m_arcStep, "requested spacing [pix]; 0=half stamp width" );
    add( "HCIA NOISE GUARD", m_noiseConfig.m_guardRadius, "extra candidate guard beyond stamp radius [pix]" );
    add( "HCIA NOISE MIN SAMPLES",
         static_cast<int>( m_noiseConfig.m_minimumSamples ),
         "minimum patches, not independent samples" );
    add( "HCIA NOISE MAX MODES", static_cast<int>( m_noiseConfig.m_maximumModes ), "maximum retained noise PCA rank" );
    add( "HCIA NOISE FLOOR FRACTION", m_noiseConfig.m_floorFraction, "floor / median training pixel variance" );
    add( "HCIA NOISE SAMPLING", std::string( "BILINEAR_NATIVE_ALIGNED" ), "fixed radius, absolute angular phase zero" );
    add( "HCIA NOISE EXCLUSION", std::string( "CIRCUMCIRCLE_PLUS_SQRT2" ), "whole training footprint dilation [pix]" );
    add( "HCIA NOISE STATUS",
         std::string( "0=valid,1=few_samples,2=zero_variance,3=invalid_filter" ),
         "NaN=no usable response" );
    add( "HCIA NOISE CALIBRATED", 0, "conditional quantities; no empirical significance calibration" );
    add( "SNRSMALL", 0, "no small-sample multiplier on conditional products" );
    add( "HCIAPSF", m_psfResponse, "PSF response manifest" );
    add( "HCIAPSM", m_psfResponseMinimumSupport, "minimum response support fraction" );
    std::string exclusions;
    for( const auto &exclusion : m_noiseExclusions )
    {
        if( !exclusions.empty() )
        {
            exclusions += ";";
        }
        exclusions += std::format( "{:.17g},{:.17g},{:.17g}", exclusion.m_row, exclusion.m_column, exclusion.m_radius );
    }
    add( "HCIA NOISE EXCLUSIONS", exclusions, "resolved row,column,radius circles; semicolon separated" );
    constexpr std::array<const char *, 11> roles{ "psf_amplitude",
                                                  "psf_sigma",
                                                  "psf_score",
                                                  "noise_samples",
                                                  "noise_modes",
                                                  "noise_floor",
                                                  "noise_status",
                                                  "psf_support",
                                                  "noise_attempted",
                                                  "noise_excluded",
                                                  "noise_incomplete" };
    constexpr std::array<const char *, 11> units{ "template amplitude",
                                                  "template amplitude",
                                                  "conditional score",
                                                  "patches",
                                                  "modes",
                                                  "image units squared",
                                                  "status code",
                                                  "fraction",
                                                  "patches",
                                                  "patches",
                                                  "patches" };
    mx::fits::fitsFile<realT, mx::verbose::vv> writer;
    for( std::size_t index = 0; index < roles.size(); ++index )
    {
        add( "HCIA NOISE PRODUCT", std::string( roles[index] ), "diagnostic quantity" );
        add( "BUNIT", std::string( units[index] ), "diagnostic units" );
        const std::string path = mx::improc::psfFilterProductPath( m_file, roles[index], false );
        const mx::error_t result = writer.write( path, m_noiseProducts[index], header );
        if( result != mx::error_t::noerror )
        {
            throw mx::exception<mx::verbose::vv>( result, "writing covariance-filter diagnostic " + path );
        }
    }
}

void hciAnalyze::filterCubePSFResponse( cubeT &cube, fitsHeaderT &scienceHeader )
{
    prepareNoiseProducts( cube );
    const std::filesystem::path manifestPath{ m_psfResponse };
    const std::string manifestName = manifestPath.filename().string();
    constexpr std::string_view manifestSuffix{ "manifest.fits" };
    if( !manifestName.ends_with( manifestSuffix ) )
    {
        throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                              "filter.psfResponse must name a PSF response manifest.fits product" );
    }
    const std::string productPrefix =
        ( manifestPath.parent_path() / manifestName.substr( 0, manifestName.size() - manifestSuffix.size() ) ).string();

    mx::fits::fitsFile<realT, mx::verbose::vv> reader;
    imageT manifest;
    fitsHeaderT manifestHeader;
    mx::error_t readResult = reader.read( manifest, manifestHeader, manifestPath.string() );
    if( readResult != mx::error_t::noerror )
    {
        throw mx::exception<mx::verbose::vv>( readResult, "reading PSF response manifest " + m_psfResponse );
    }
    if( manifestHeader.count( "KLIP PSF PRODUCT" ) != 0 )
    {
        filterCubeKLIPPSFResponse( cube, scienceHeader, productPrefix, manifest, manifestHeader );
        writeNoiseProducts( scienceHeader );
        return;
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
                                              "filter.psfResponse does not name a complete P4 PSF manifest" );
    }
    const int schema = manifestHeader["P4 PSF PRODUCT SCHEMA"].value<int>();
    const int modeCount = manifestHeader["P4 PSF MODE COUNT"].value<int>();
    const int stampSize = manifestHeader["P4 PSF STAMP SIZE"].value<int>();
    m_psfResponseMinimumSupport = manifestHeader["P4 PSF FILTER MIN GOOD FRACTION"].value<realT>();
    if( schema < 1 || modeCount <= 0 || modeCount != cube.planes() || stampSize <= 0 || stampSize % 2 == 0 ||
        !std::isfinite( m_psfResponseMinimumSupport ) || m_psfResponseMinimumSupport < 0 ||
        m_psfResponseMinimumSupport > 1 )
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
            const mx::improc::P4PSFFilterResult result = applyPSFFilter( cube.image( mode ),
                                                                         models.image( static_cast<int>( source ) ),
                                                                         responseValidity,
                                                                         sourceRow,
                                                                         sourceColumn,
                                                                         mode );
            if( result.valid && std::isfinite( result.amplitude ) &&
                std::abs( result.amplitude ) <= std::numeric_limits<realT>::max() )
            {
                filtered.image( mode )( sourceRow, sourceColumn ) = static_cast<realT>( result.amplitude );
            }
        }
    }
    cube = std::move( filtered );
    writeNoiseProducts( scienceHeader );
}

void hciAnalyze::filterCubeKLIPPSFResponse( cubeT &cube,
                                            fitsHeaderT &scienceHeader,
                                            const std::string &productPrefix,
                                            const imageT &manifest,
                                            fitsHeaderT &manifestHeader )
{
    const auto requireHeader = [&manifestHeader]( const std::string &keyword )
    {
        if( manifestHeader.count( keyword ) == 0 )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                  "KLIP PSF response manifest is missing " + keyword );
        }
    };
    for( const std::string &keyword : { "KLIP PSF PRODUCT SCHEMA",
                                        "KLIP PSF PRODUCT",
                                        "KLIP PSF COMPLETE",
                                        "KLIP PSF STAMP SIZE",
                                        "KLIP PSF FILTER MIN GOOD FRACTION",
                                        "KLIP PSF SAMPLE RADII",
                                        "KLIP PSF SPATIAL MODEL",
                                        "NMODES",
                                        "REGMINR",
                                        "REGMAXR" } )
    {
        requireHeader( keyword );
    }
    if( manifest.rows() != 1 || manifest.cols() != 1 || manifest( 0, 0 ) != 1 ||
        !manifestHeader["KLIP PSF PRODUCT"].String().starts_with( "MANIFEST" ) ||
        manifestHeader["KLIP PSF COMPLETE"].value<int>() != 1 )
    {
        throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                              "filter.psfResponse does not name a complete KLIP PSF manifest" );
    }

    const int schema = manifestHeader["KLIP PSF PRODUCT SCHEMA"].value<int>();
    const int stampSize = manifestHeader["KLIP PSF STAMP SIZE"].value<int>();
    m_psfResponseMinimumSupport = manifestHeader["KLIP PSF FILTER MIN GOOD FRACTION"].value<realT>();
    const std::vector<realT> radiusValues = headerVector( manifestHeader, "KLIP PSF SAMPLE RADII" );
    const std::vector<realT> responseModes = headerVector( manifestHeader, "NMODES" );
    const std::vector<realT> minimumRadii = headerVector( manifestHeader, "REGMINR" );
    const std::vector<realT> maximumRadii = headerVector( manifestHeader, "REGMAXR" );
    const std::string spatialModel = manifestHeader["KLIP PSF SPATIAL MODEL"].String();
    const bool pixelExact = spatialModel.starts_with( "PIXEL_EXACT" );
    const bool radialSchema = schema == 1 && !pixelExact && !radiusValues.empty();
    const bool pixelSchema = schema == 2 && pixelExact;
    if( ( !radialSchema && !pixelSchema ) || stampSize <= 0 || stampSize % 2 == 0 ||
        !std::isfinite( m_psfResponseMinimumSupport ) || m_psfResponseMinimumSupport < 0 ||
        m_psfResponseMinimumSupport > 1 ||
        radiusValues.size() > static_cast<std::size_t>( std::numeric_limits<int>::max() ) ||
        responseModes.size() != static_cast<std::size_t>( cube.planes() ) || minimumRadii.empty() ||
        minimumRadii.size() != maximumRadii.size() )
    {
        throw mx::exception<mx::verbose::vv>(
            mx::error_t::invalidconfig,
            "KLIP PSF schema, modes, odd stamp size, radial grid, regions, or minimum support is incompatible with "
            "the input" );
    }

    std::vector<double> radii;
    radii.reserve( radiusValues.size() );
    for( const realT radius : radiusValues )
    {
        if( !std::isfinite( radius ) || radius < 0 || ( !radii.empty() && radius <= radii.back() ) )
        {
            throw mx::exception<mx::verbose::vv>(
                mx::error_t::invalidconfig,
                "KLIP PSF response radii must be finite, nonnegative, and strictly increasing" );
        }
        radii.push_back( static_cast<double>( radius ) );
    }

    for( std::size_t region = 0; region < minimumRadii.size(); ++region )
    {
        if( !std::isfinite( minimumRadii[region] ) || !std::isfinite( maximumRadii[region] ) ||
            minimumRadii[region] < 0 || maximumRadii[region] <= minimumRadii[region] ||
            ( region != 0 && minimumRadii[region] < maximumRadii[region - 1] ) )
        {
            throw mx::exception<mx::verbose::vv>(
                mx::error_t::invalidconfig,
                "KLIP PSF response regions must be finite, positive-width, ordered, and non-overlapping" );
        }
    }

    const bool regionAware = spatialModel.starts_with( "REGION_RADIAL_LINEAR" );
    const bool exactAzimuthal = spatialModel.starts_with( "EXACT_AZIMUTHAL" );
    if( !pixelExact && !regionAware && !exactAzimuthal && !spatialModel.starts_with( "RADIAL_LINEAR" ) )
    {
        throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                              "KLIP PSF manifest has an unsupported spatial model" );
    }
    double exactSampleAngle{ 0 };
    int exactSampleRow{ -1 };
    int exactSampleColumn{ -1 };
    if( exactAzimuthal )
    {
        requireHeader( "KLIP PSF SAMPLE ANGLE" );
        requireHeader( "KLIP PSF EXACT ROW" );
        requireHeader( "KLIP PSF EXACT COLUMN" );
        exactSampleAngle = manifestHeader["KLIP PSF SAMPLE ANGLE"].value<double>();
        exactSampleRow = manifestHeader["KLIP PSF EXACT ROW"].value<int>();
        exactSampleColumn = manifestHeader["KLIP PSF EXACT COLUMN"].value<int>();
        if( radii.size() != 1 || !std::isfinite( exactSampleAngle ) )
        {
            throw mx::exception<mx::verbose::vv>(
                mx::error_t::invalidconfig,
                "an exact azimuthal KLIP PSF response requires one radius and a finite sample angle" );
        }
    }
    std::vector<std::size_t> radiusRegions( radii.size(), 0 );
    if( regionAware )
    {
        std::vector<std::size_t> nodesPerRegion( minimumRadii.size(), 0 );
        for( std::size_t radiusIndex = 0; radiusIndex < radii.size(); ++radiusIndex )
        {
            std::optional<std::size_t> owner;
            for( std::size_t region = 0; region < minimumRadii.size(); ++region )
            {
                if( radii[radiusIndex] >= static_cast<double>( minimumRadii[region] ) &&
                    radii[radiusIndex] < static_cast<double>( maximumRadii[region] ) )
                {
                    if( owner.has_value() )
                    {
                        throw mx::exception<mx::verbose::vv>(
                            mx::error_t::invalidconfig,
                            "KLIP PSF radial node belongs to more than one interpolation region" );
                    }
                    owner = region;
                }
            }
            if( !owner.has_value() )
            {
                throw mx::exception<mx::verbose::vv>(
                    mx::error_t::invalidconfig,
                    "KLIP PSF radial node does not belong to an interpolation region" );
            }
            radiusRegions[radiusIndex] = *owner;
            ++nodesPerRegion[*owner];
        }
        if( std::any_of( nodesPerRegion.begin(),
                         nodesPerRegion.end(),
                         []( std::size_t count ) { return count == 0; } ) )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                  "every KLIP PSF interpolation region requires a radial node" );
        }
    }

    const std::vector<realT> scienceModes = scienceHeader.count( "NMODES" )
                                                ? headerVector( scienceHeader, "NMODES" )
                                                : headerVector( scienceHeader, "FRACT NMODES" );
    const std::vector<realT> scienceMinimumRadii = headerVector( scienceHeader, "REGMINR" );
    const std::vector<realT> scienceMaximumRadii = headerVector( scienceHeader, "REGMAXR" );
    if( scienceModes.size() != responseModes.size() || scienceMinimumRadii.size() != minimumRadii.size() ||
        scienceMaximumRadii.size() != maximumRadii.size() )
    {
        throw mx::exception<mx::verbose::vv>(
            mx::error_t::invalidconfig,
            "science and KLIP PSF products require matching mode labels and regions" );
    }
    for( std::size_t region = 0; region < minimumRadii.size(); ++region )
    {
        const realT scale = std::max( { realT{ 1 },
                                        std::abs( minimumRadii[region] ),
                                        std::abs( maximumRadii[region] ),
                                        std::abs( scienceMinimumRadii[region] ),
                                        std::abs( scienceMaximumRadii[region] ) } );
        if( !std::isfinite( scienceMinimumRadii[region] ) || !std::isfinite( scienceMaximumRadii[region] ) ||
            std::abs( scienceMinimumRadii[region] - minimumRadii[region] ) >
                8 * std::numeric_limits<realT>::epsilon() * scale ||
            std::abs( scienceMaximumRadii[region] - maximumRadii[region] ) >
                8 * std::numeric_limits<realT>::epsilon() * scale )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                  "science and KLIP PSF region bounds differ at region " +
                                                      std::to_string( region ) );
        }
    }
    std::vector<int> modeCounts( responseModes.size() );
    for( std::size_t mode = 0; mode < responseModes.size(); ++mode )
    {
        const realT scale = std::max( { realT{ 1 }, std::abs( scienceModes[mode] ), std::abs( responseModes[mode] ) } );
        if( !std::isfinite( scienceModes[mode] ) || !std::isfinite( responseModes[mode] ) || responseModes[mode] <= 0 ||
            responseModes[mode] != std::trunc( responseModes[mode] ) ||
            responseModes[mode] > static_cast<realT>( std::numeric_limits<int>::max() ) ||
            std::abs( scienceModes[mode] - responseModes[mode] ) > 8 * std::numeric_limits<realT>::epsilon() * scale )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                  "science and KLIP PSF mode labels differ at plane " +
                                                      std::to_string( mode ) );
        }
        modeCounts[mode] = static_cast<int>( responseModes[mode] );
    }

    if( pixelExact )
    {
        requireHeader( "KLIP PSF MEASUREMENT COUNT" );
        const std::string sourceCountString = manifestHeader["KLIP PSF MEASUREMENT COUNT"].String();
        std::size_t parsedCharacters{ 0 };
        std::size_t sourceCount{ 0 };
        try
        {
            sourceCount = std::stoull( sourceCountString, &parsedCharacters );
        }
        catch( const std::exception & )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                  "KLIP exact response measurement count is not an integer" );
        }
        const bool trailingCharactersAreWhitespace =
            std::all_of( sourceCountString.begin() + static_cast<std::ptrdiff_t>( parsedCharacters ),
                         sourceCountString.end(),
                         []( unsigned char character ) { return std::isspace( character ) != 0; } );
        if( !trailingCharactersAreWhitespace || sourceCount == 0 ||
            sourceCount > static_cast<std::size_t>( std::numeric_limits<int>::max() ) )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                  "KLIP exact response measurement count is outside output range" );
        }

        mx::fits::fitsFile<realT, mx::verbose::vv> reader;
        imageT coordinates;
        fitsHeaderT coordinateHeader;
        const std::string coordinatePath = productPrefix + "coordinates.fits";
        mx::error_t readResult = reader.read( coordinates, coordinateHeader, coordinatePath );
        if( readResult != mx::error_t::noerror )
        {
            throw mx::exception<mx::verbose::vv>( readResult, "reading KLIP exact coordinates " + coordinatePath );
        }
        if( coordinates.rows() != static_cast<Eigen::Index>( sourceCount ) || coordinates.cols() != 4 ||
            coordinateHeader.count( "KLIP PSF PRODUCT SCHEMA" ) == 0 ||
            coordinateHeader["KLIP PSF PRODUCT SCHEMA"].value<int>() != schema ||
            coordinateHeader.count( "KLIP PSF PRODUCT" ) == 0 ||
            !coordinateHeader["KLIP PSF PRODUCT"].String().starts_with( "PIXEL_COORDINATES" ) )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                  "KLIP exact response coordinates are inconsistent" );
        }

        std::vector<std::pair<int, int>> sourceCoordinates( sourceCount );
        std::vector<std::uint8_t> coordinateOwners( static_cast<std::size_t>( cube.rows() ) * cube.cols(), 0 );
        const double centerRow = 0.5 * static_cast<double>( cube.rows() - 1 );
        const double centerColumn = 0.5 * static_cast<double>( cube.cols() - 1 );
        for( std::size_t source = 0; source < sourceCount; ++source )
        {
            const realT rowValue = coordinates( static_cast<Eigen::Index>( source ), 0 );
            const realT columnValue = coordinates( static_cast<Eigen::Index>( source ), 1 );
            const realT regionValue = coordinates( static_cast<Eigen::Index>( source ), 2 );
            const realT indexValue = coordinates( static_cast<Eigen::Index>( source ), 3 );
            if( !std::isfinite( rowValue ) || !std::isfinite( columnValue ) || !std::isfinite( regionValue ) ||
                !std::isfinite( indexValue ) || rowValue != std::trunc( rowValue ) ||
                columnValue != std::trunc( columnValue ) || regionValue != std::trunc( regionValue ) ||
                indexValue != static_cast<realT>( source ) || rowValue < 0 || rowValue >= cube.rows() ||
                columnValue < 0 || columnValue >= cube.cols() || regionValue < 0 ||
                regionValue >= static_cast<realT>( minimumRadii.size() ) )
            {
                throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                      "KLIP exact response coordinate row is invalid" );
            }
            const int row = static_cast<int>( rowValue );
            const int column = static_cast<int>( columnValue );
            const std::size_t region = static_cast<std::size_t>( regionValue );
            const double radius =
                std::hypot( static_cast<double>( row ) - centerRow, static_cast<double>( column ) - centerColumn );
            const std::size_t owner =
                static_cast<std::size_t>( row ) + static_cast<std::size_t>( cube.rows() ) * column;
            if( coordinateOwners[owner] != 0 || radius < static_cast<double>( minimumRadii[region] ) ||
                radius >= static_cast<double>( maximumRadii[region] ) )
            {
                throw mx::exception<mx::verbose::vv>(
                    mx::error_t::invalidconfig,
                    "KLIP exact response coordinate is duplicated or outside its declared region" );
            }
            coordinateOwners[owner] = 1;
            sourceCoordinates[source] = { row, column };
        }

        cubeT filtered( cube.rows(), cube.cols(), cube.planes() );
        filtered.cube().setConstant( std::numeric_limits<realT>::quiet_NaN() );
        const int responseCenter = stampSize / 2;
        for( int mode = 0; mode < cube.planes(); ++mode )
        {
            const std::string modeIndex = std::format( "{:03d}", mode );
            const std::string responsePath = productPrefix + "mode" + modeIndex + "_pixel_response.fits";
            const std::string validityPath = productPrefix + "mode" + modeIndex + "_pixel_validity.fits";
            cubeT responseCube;
            cubeT validityCube;
            fitsHeaderT responseHeader;
            fitsHeaderT validityHeader;
            readResult = reader.read( responseCube, responseHeader, responsePath );
            if( readResult != mx::error_t::noerror )
            {
                throw mx::exception<mx::verbose::vv>( readResult, "reading KLIP exact response " + responsePath );
            }
            readResult = reader.read( validityCube, validityHeader, validityPath );
            if( readResult != mx::error_t::noerror )
            {
                throw mx::exception<mx::verbose::vv>( readResult, "reading KLIP exact validity " + validityPath );
            }
            if( responseCube.rows() != stampSize || responseCube.cols() != stampSize ||
                responseCube.planes() != static_cast<int>( sourceCount ) || validityCube.rows() != stampSize ||
                validityCube.cols() != stampSize || validityCube.planes() != responseCube.planes() ||
                responseHeader.count( "KLIP PSF PRODUCT SCHEMA" ) == 0 ||
                validityHeader.count( "KLIP PSF PRODUCT SCHEMA" ) == 0 ||
                responseHeader["KLIP PSF PRODUCT SCHEMA"].value<int>() != schema ||
                validityHeader["KLIP PSF PRODUCT SCHEMA"].value<int>() != schema ||
                responseHeader.count( "KLIP PSF PRODUCT" ) == 0 || validityHeader.count( "KLIP PSF PRODUCT" ) == 0 ||
                !responseHeader["KLIP PSF PRODUCT"].String().starts_with( "PIXEL_RESPONSE" ) ||
                !validityHeader["KLIP PSF PRODUCT"].String().starts_with( "PIXEL_VALIDITY" ) ||
                responseHeader.count( "KLIP PSF MODE COUNT" ) == 0 ||
                validityHeader.count( "KLIP PSF MODE COUNT" ) == 0 ||
                responseHeader["KLIP PSF MODE COUNT"].value<int>() != modeCounts[static_cast<std::size_t>( mode )] ||
                validityHeader["KLIP PSF MODE COUNT"].value<int>() != modeCounts[static_cast<std::size_t>( mode )] )
            {
                throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                      "KLIP exact response or validity product is inconsistent" );
            }

            for( std::size_t source = 0; source < sourceCount; ++source )
            {
                mx::improc::P4PSFFilter::validityT responseValidity( stampSize, stampSize );
                for( int column = 0; column < stampSize; ++column )
                {
                    for( int row = 0; row < stampSize; ++row )
                    {
                        const realT validity = validityCube.image( static_cast<int>( source ) )( row, column );
                        const realT response = responseCube.image( static_cast<int>( source ) )( row, column );
                        if( !std::isfinite( validity ) || ( validity != 0 && validity != 1 ) ||
                            ( validity != 0 && !std::isfinite( response ) ) )
                        {
                            throw mx::exception<mx::verbose::vv>(
                                mx::error_t::invalidconfig,
                                "KLIP exact response contains invalid values or non-binary validity" );
                        }
                        responseValidity( row, column ) = validity != 0 ? 1 : 0;
                    }
                }
                if( responseValidity( responseCenter, responseCenter ) == 0 )
                {
                    continue;
                }
                const auto [sourceRow, sourceColumn] = sourceCoordinates[source];
                const mx::improc::P4PSFFilterResult result =
                    applyPSFFilter( cube.image( mode ),
                                    responseCube.image( static_cast<int>( source ) ),
                                    responseValidity,
                                    sourceRow,
                                    sourceColumn,
                                    mode );
                if( result.valid && std::isfinite( result.amplitude ) &&
                    std::abs( result.amplitude ) <= std::numeric_limits<realT>::max() )
                {
                    filtered.image( mode )( sourceRow, sourceColumn ) = static_cast<realT>( result.amplitude );
                }
            }
        }
        cube = std::move( filtered );
        return;
    }

    mx::fits::fitsFile<realT, mx::verbose::vv> reader;
    cubeT filtered( cube.rows(), cube.cols(), cube.planes() );
    filtered.cube().setConstant( std::numeric_limits<realT>::quiet_NaN() );
    const double centerRow = 0.5 * static_cast<double>( cube.rows() - 1 );
    const double centerColumn = 0.5 * static_cast<double>( cube.cols() - 1 );
    if( exactAzimuthal )
    {
        const double deltaRow = static_cast<double>( exactSampleRow ) - centerRow;
        const double deltaColumn = static_cast<double>( exactSampleColumn ) - centerColumn;
        const double anchorRadius = std::hypot( deltaRow, deltaColumn );
        const double anchorAngle = std::atan2( deltaRow, deltaColumn );
        const double angleDifference =
            std::atan2( std::sin( anchorAngle - exactSampleAngle ), std::cos( anchorAngle - exactSampleAngle ) );
        const double radiusScale = std::max( 1.0, anchorRadius );
        const double geometryTolerance = 8 * static_cast<double>( std::numeric_limits<realT>::epsilon() );
        bool anchorInRegion{ false };
        for( std::size_t region = 0; region < minimumRadii.size(); ++region )
        {
            if( anchorRadius >= static_cast<double>( minimumRadii[region] ) &&
                anchorRadius < static_cast<double>( maximumRadii[region] ) )
            {
                anchorInRegion = true;
                break;
            }
        }
        if( exactSampleRow < 0 || exactSampleRow >= cube.rows() || exactSampleColumn < 0 ||
            exactSampleColumn >= cube.cols() ||
            std::abs( anchorRadius - radii.front() ) > geometryTolerance * radiusScale ||
            std::abs( angleDifference ) > geometryTolerance || !anchorInRegion )
        {
            throw mx::exception<mx::verbose::vv>(
                mx::error_t::invalidconfig,
                "exact azimuthal KLIP PSF anchor is inconsistent with the science image and response geometry" );
        }
    }
    const int responseCenter = stampSize / 2;
    for( int mode = 0; mode < cube.planes(); ++mode )
    {
        const std::string modeIndex = std::format( "{:03d}", mode );
        const std::string responsePath = productPrefix + "mode" + modeIndex + "_radial_response.fits";
        const std::string validityPath = productPrefix + "mode" + modeIndex + "_radial_validity.fits";
        cubeT responseCube;
        cubeT validityCube;
        fitsHeaderT responseHeader;
        fitsHeaderT validityHeader;
        mx::error_t readResult = reader.read( responseCube, responseHeader, responsePath );
        if( readResult != mx::error_t::noerror )
        {
            throw mx::exception<mx::verbose::vv>( readResult, "reading KLIP PSF response " + responsePath );
        }
        readResult = reader.read( validityCube, validityHeader, validityPath );
        if( readResult != mx::error_t::noerror )
        {
            throw mx::exception<mx::verbose::vv>( readResult, "reading KLIP PSF validity " + validityPath );
        }
        if( responseCube.rows() != stampSize || responseCube.cols() != stampSize ||
            responseCube.planes() != static_cast<int>( radii.size() ) || validityCube.rows() != stampSize ||
            validityCube.cols() != stampSize || validityCube.planes() != responseCube.planes() ||
            responseHeader.count( "KLIP PSF PRODUCT SCHEMA" ) == 0 ||
            validityHeader.count( "KLIP PSF PRODUCT SCHEMA" ) == 0 ||
            responseHeader["KLIP PSF PRODUCT SCHEMA"].value<int>() != schema ||
            validityHeader["KLIP PSF PRODUCT SCHEMA"].value<int>() != schema ||
            responseHeader.count( "KLIP PSF PRODUCT" ) == 0 || validityHeader.count( "KLIP PSF PRODUCT" ) == 0 ||
            !responseHeader["KLIP PSF PRODUCT"].String().starts_with( "RADIAL_RESPONSE" ) ||
            !validityHeader["KLIP PSF PRODUCT"].String().starts_with( "RADIAL_VALIDITY" ) ||
            responseHeader.count( "KLIP PSF MODE COUNT" ) == 0 || validityHeader.count( "KLIP PSF MODE COUNT" ) == 0 ||
            responseHeader["KLIP PSF MODE COUNT"].value<int>() != modeCounts[static_cast<std::size_t>( mode )] ||
            validityHeader["KLIP PSF MODE COUNT"].value<int>() != modeCounts[static_cast<std::size_t>( mode )] )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                  "KLIP PSF radial response or validity product is inconsistent" );
        }

        std::vector<mx::improc::RadialPSFModel::imageT> responses( radii.size() );
        std::vector<mx::improc::RadialPSFModel::validityT> validities( radii.size() );
        std::vector<mx::improc::RadialPSFSample> samples( radii.size() );
        for( std::size_t radiusIndex = 0; radiusIndex < radii.size(); ++radiusIndex )
        {
            responses[radiusIndex].resize( stampSize, stampSize );
            validities[radiusIndex].resize( stampSize, stampSize );
            for( int column = 0; column < stampSize; ++column )
            {
                for( int row = 0; row < stampSize; ++row )
                {
                    const realT validity = validityCube.image( static_cast<int>( radiusIndex ) )( row, column );
                    const realT response = responseCube.image( static_cast<int>( radiusIndex ) )( row, column );
                    if( !std::isfinite( validity ) || ( validity != 0 && validity != 1 ) ||
                        ( validity != 0 && !std::isfinite( response ) ) )
                    {
                        throw mx::exception<mx::verbose::vv>(
                            mx::error_t::invalidconfig,
                            "KLIP PSF radial response contains invalid values or non-binary validity" );
                    }
                    responses[radiusIndex]( row, column ) = validity != 0 ? response : 0;
                    validities[radiusIndex]( row, column ) = validity != 0 ? 1 : 0;
                }
            }
            samples[radiusIndex] = { radiusIndex, radiusIndex, radii[radiusIndex], 0, radiusRegions[radiusIndex] };
        }

        std::optional<mx::improc::RadialPSFModel> radialModel;
        if( exactAzimuthal )
        {
            radialModel.reset();
        }
        else if( regionAware )
        {
            radialModel.emplace( radii, radiusRegions, stampSize, stampSize );
        }
        else
        {
            radialModel.emplace( radii, stampSize, stampSize );
        }
        if( radialModel.has_value() )
        {
            radialModel->fit( responses, validities, samples );
        }

        mx::improc::RadialPSFModel::imageT response;
        mx::improc::RadialPSFModel::validityT responseValidity;
        for( int column = 0; column < cube.cols(); ++column )
        {
            const double deltaColumn = static_cast<double>( column ) - centerColumn;
            for( int row = 0; row < cube.rows(); ++row )
            {
                const double deltaRow = static_cast<double>( row ) - centerRow;
                const double radius = std::hypot( deltaRow, deltaColumn );
                std::optional<std::size_t> activeRegion;
                for( std::size_t region = 0; region < minimumRadii.size(); ++region )
                {
                    if( radius >= static_cast<double>( minimumRadii[region] ) &&
                        radius < static_cast<double>( maximumRadii[region] ) )
                    {
                        activeRegion = region;
                        break;
                    }
                }
                if( !activeRegion.has_value() )
                {
                    continue;
                }
                const double requestedAngle = std::atan2( deltaRow, deltaColumn );
                if( exactAzimuthal )
                {
                    mx::improc::RadialPSFModel::rotate( response,
                                                        responseValidity,
                                                        responses.front(),
                                                        validities.front(),
                                                        exactSampleAngle - requestedAngle );
                }
                else if( regionAware )
                {
                    radialModel->response( response, responseValidity, radius, requestedAngle, *activeRegion );
                }
                else
                {
                    radialModel->response( response, responseValidity, radius, requestedAngle );
                }
                if( responseValidity( responseCenter, responseCenter ) == 0 )
                {
                    continue;
                }
                const mx::improc::P4PSFFilterResult result =
                    applyPSFFilter( cube.image( mode ), response, responseValidity, row, column, mode );
                if( result.valid && std::isfinite( result.amplitude ) &&
                    std::abs( result.amplitude ) <= std::numeric_limits<realT>::max() )
                {
                    filtered.image( mode )( row, column ) = static_cast<realT>( result.amplitude );
                }
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
    if( !m_psfResponse.empty() )
    {
        filterCubePSFResponse( cube, header );
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
              << "  PSF response manifest: " << ( m_psfResponse.empty() ? "disabled" : m_psfResponse ) << '\n';
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
    addHeader( "HCIAPSF", m_psfResponse, "PSF response manifest" );
    addHeader( "HCIAPSM", m_psfResponseMinimumSupport, "PSF response minimum support" );
    if( m_noiseModel != "identity" )
    {
        addHeader( "HCIA NOISE MODEL", m_noiseModel, "noise weighting applied before empirical annular SNR" );
        addHeader( "HCIA NOISE MAX MODES", static_cast<int>( m_noiseConfig.m_maximumModes ), "maximum noise PCA rank" );
        addHeader( "HCIA NOISE FLOOR FRACTION",
                   m_noiseConfig.m_floorFraction,
                   "floor / median training pixel variance" );
        addHeader( "HCIA NOISE CALIBRATED", 0, "covariance-estimation uncertainty is not empirically calibrated" );
    }
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
