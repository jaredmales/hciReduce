/** \file p4ReductionPrecisionBenchmark.cpp
 * \brief Runs a configured production-layout P4 reduction with an experimental precision policy.
 */

#include <bit>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

#define HCIREDUCE_P4REDUCE_NO_MAIN
#include "src/apps/p4Reduce.cpp"

/// Persist detector inputs and pre-storage residuals in a simple benchmark-only binary format.
class P4DetectorCapture final : public mx::improc::detail::P4DetectorFitObserver
{
    std::filesystem::path m_directory; ///< Fresh output directory, owned by this benchmark run.

    bool m_inputs;                     ///< Whether to include the sampled predictor matrix and target before residuals.

    std::mutex m_mutex;                ///< Serializes file creation and duplicate-coordinate checks across workers.

  public:
    /// Create a fresh capture directory.
    P4DetectorCapture( const std::string &directory, /**< [in] nonexistent capture path */
                       bool inputs /**< [in] include ingress inputs as well as residuals */ );

    /// Write one coordinate's metadata and column-major FP64 arrays with checked I/O.
    void
    observe( const mx::improc::P4PixelCoordinate &coordinate, /**< [in] detector target pixel */
             const mx::improc::P4PCA::matrixT &predictors,    /**< [in] sampled FP64 predictors */
             const mx::improc::P4PCA::vectorT &target,        /**< [in] sampled FP64 target */
             const std::vector<int> &modes,                   /**< [in] retained counts */
             double rankTolerance,                            /**< [in] relative rank threshold */
             const mx::improc::P4PCAResult &result /**< [in] kernel residuals before image conversion */ ) override;
};

P4DetectorCapture::P4DetectorCapture( const std::string &directory, bool inputs )
    : m_directory( directory ), m_inputs( inputs )
{
    if( !std::filesystem::create_directories( m_directory ) )
    {
        throw std::invalid_argument( "detector capture requires a new directory" );
    }
}

void P4DetectorCapture::observe( const mx::improc::P4PixelCoordinate &coordinate,
                                 const mx::improc::P4PCA::matrixT &predictors,
                                 const mx::improc::P4PCA::vectorT &target,
                                 const std::vector<int> &modes,
                                 double rankTolerance,
                                 const mx::improc::P4PCAResult &result )
{
    const std::lock_guard<std::mutex> lock( m_mutex );
    const std::string stem = std::format( "fit_{}_{}", coordinate.row(), coordinate.column() );
    const auto metadataPath = m_directory / ( stem + ".json" );
    if( std::filesystem::exists( metadataPath ) )
    {
        throw std::runtime_error( "duplicate detector coordinate in capture" );
    }
    std::ofstream binary;
    binary.exceptions( std::ios::failbit | std::ios::badbit );
    binary.open( m_directory / ( stem + ".bin" ), std::ios::binary );
    if( m_inputs )
    {
        binary.write( reinterpret_cast<const char *>( predictors.data() ), predictors.size() * sizeof( double ) );
        binary.write( reinterpret_cast<const char *>( target.data() ), target.size() * sizeof( double ) );
    }
    binary.write( reinterpret_cast<const char *>( result.residuals.data() ),
                  result.residuals.size() * sizeof( double ) );
    binary.close();
    std::ofstream metadata;
    metadata.exceptions( std::ios::failbit | std::ios::badbit );
    metadata.open( metadataPath );
    metadata << std::setprecision( std::numeric_limits<double>::max_digits10 )
             << "{\"schema\":1,\"row\":" << coordinate.row() << ",\"column\":" << coordinate.column()
             << ",\"rows\":" << predictors.rows() << ",\"columns\":" << predictors.cols()
             << ",\"inputs\":" << ( m_inputs ? "true" : "false" ) << ",\"dtype\":\""
             << ( std::endian::native == std::endian::little ? "<f8" : ">f8" )
             << "\",\"order\":\"F\",\"rank_tolerance\":" << rankTolerance << ",\"rank\":" << result.numericalRank
             << ",\"modes\":[";
    for( std::size_t index = 0; index < modes.size(); ++index )
    {
        metadata << ( index == 0 ? "" : "," ) << modes[index];
    }
    metadata << "],\"supported\":[";
    for( std::size_t index = 0; index < modes.size(); ++index )
    {
        metadata << ( index == 0 ? "" : "," )
                 << ( result.modeStatus[index] == mx::improc::P4PCAModeStatus::rankSupported ? "true" : "false" );
    }
    metadata << "]}\n";
    metadata.close();
}

/// Benchmark-only p4Reduce adapter selecting one internal P4 PCA precision policy.
class p4ReductionPrecisionBenchmark final : public p4Reduce
{
  protected:
    std::string m_precisionName; ///< Required frozen candidate identifier supplied by the benchmark caller.

    std::string m_captureDirectory; ///< Optional new detector-capture directory.

    bool m_captureInputs{ false };  ///< Include ingress matrices when detector capture is enabled.

    mx::improc::detail::P4PCAPrecisionPolicy m_precisionPolicy{
        mx::improc::detail::P4PCAPrecisionPolicy::doubleDouble }; ///< Validated internal precision selection.

    /// Parse one exact benchmark policy spelling.
    static mx::improc::detail::P4PCAPrecisionPolicy
    parsePrecisionPolicy( const std::string &value /**< [in] exact frozen P4 candidate identifier */ );

  public:
    /// Construct a runner without implicit global or local p4Reduce configuration files.
    p4ReductionPrecisionBenchmark();

    /// Register the production p4Reduce controls and the required benchmark precision selector.
    void setupConfig() override;

    /// Consume the precision selector before running the production unused-target audit.
    void loadConfig() override;

    /// Validate production dispatch, the benchmark policy, and unsupported optimizer use.
    void checkConfig() override;

    /// Execute and time one production-layout reduction through the selected experimental policy.
    int execute() override;
};

mx::improc::detail::P4PCAPrecisionPolicy p4ReductionPrecisionBenchmark::parsePrecisionPolicy( const std::string &value )
{
    if( value == "P4-D64" )
    {
        return mx::improc::detail::P4PCAPrecisionPolicy::doubleDouble;
    }
    if( value == "P4-M32D64" )
    {
        return mx::improc::detail::P4PCAPrecisionPolicy::floatDouble;
    }
    if( value == "P4-F32" )
    {
        return mx::improc::detail::P4PCAPrecisionPolicy::floatFloat;
    }
    throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                          "precision must be exactly one of P4-D64, P4-M32D64, or P4-F32" );
}

p4ReductionPrecisionBenchmark::p4ReductionPrecisionBenchmark()
{
    m_configPathGlobal.clear();
    m_configPathGlobal_env.clear();
    m_configPathLocal.clear();
    m_configPathLocal_env.clear();
}

void p4ReductionPrecisionBenchmark::setupConfig()
{
    p4Reduce::setupConfig();
    config.add( "benchmark.precision",
                "",
                "precision",
                mx::app::argType::Required,
                "benchmark",
                "precision",
                true,
                "string",
                "required frozen P4 candidate: exactly P4-D64, P4-M32D64, or P4-F32; non-D64 candidates also require "
                "p4.memoryFraction=0" );
    config.add( "benchmark.captureDetectors",
                "",
                "capture-detectors",
                mx::app::argType::Required,
                "benchmark",
                "captureDetectors",
                false,
                "string",
                "optional new directory for direct local detector inputs and residuals" );
    config.add( "benchmark.captureInputs",
                "",
                "capture-inputs",
                mx::app::argType::Optional,
                "benchmark",
                "captureInputs",
                false,
                "bool",
                "include ingress predictors and targets in detector captures; default false" );
}

void p4ReductionPrecisionBenchmark::loadConfig()
{
    config( m_precisionName, "benchmark.precision" );
    config( m_captureDirectory, "benchmark.captureDetectors" );
    mx::improc::loadBoolConfig<mx::verbose::vv>( config, m_captureInputs, "benchmark.captureInputs" );
    p4Reduce::loadConfig();
}

void p4ReductionPrecisionBenchmark::checkConfig()
{
    p4Reduce::checkConfig();
    m_precisionPolicy = parsePrecisionPolicy( m_precisionName );
    if( m_optimizeEnabled )
    {
        throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                              "the P4 precision benchmark does not support p4Optimize.enabled" );
    }
    if( !m_captureDirectory.empty() && ( m_obs.m_localStampSize <= 0 || m_obs.m_numberImages != 0 ||
                                         m_obs.m_temporalPredictor != mx::improc::P4TemporalPredictor::none ||
                                         m_obs.m_excludeMethod != mx::improc::HCI::exclude::none ) )
    {
        throw std::invalid_argument(
            "detector capture requires local direct P4 without temporal augmentation/exclusion" );
    }
    if( m_captureInputs && m_captureDirectory.empty() )
    {
        throw std::invalid_argument( "capture-inputs requires capture-detectors" );
    }
}

int p4ReductionPrecisionBenchmark::execute()
{
    const auto begin = std::chrono::steady_clock::now();
    std::unique_ptr<P4DetectorCapture> capture;
    if( !m_captureDirectory.empty() )
    {
        capture = std::make_unique<P4DetectorCapture>( m_captureDirectory, m_captureInputs );
    }
    const int result = mx::improc::detail::p4ReductionReduceExperimental( m_obs, m_precisionPolicy, capture.get() );
    const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - begin;

    std::cout << std::setprecision( std::numeric_limits<double>::max_digits10 )
              << "p4_reduction_precision_benchmark policy=" << m_precisionName << " elapsed_seconds=" << elapsed.count()
              << " status=" << result << '\n';
    if( result == 0 && m_showTiming )
    {
        m_obs.dump_times();
    }
    return result;
}

/// Run the benchmark application with p4Reduce-compatible exception reporting.
/** The application framework's help status is normalized to `EXIT_SUCCESS`; any reduction status is otherwise
 * propagated.
 *
 * \returns `EXIT_SUCCESS` on help or successful execution, otherwise `EXIT_FAILURE`.
 */
int p4ReductionPrecisionBenchmarkMain( int argc, /**< [in] standard command-line argument count */
                                       char **argv /**< [in] standard command-line argument vector */ )
{
    const std::string argv0 =
        argc > 0 && argv != nullptr && argv[0] != nullptr ? argv[0] : "p4ReductionPrecisionBenchmark";
    p4ReductionPrecisionBenchmark benchmark;

    try
    {
        const int result = benchmark.main( argc, argv );
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

/// P4 reduction precision benchmark process entry point.
int main( int argc, /**< [in] standard command-line argument count */
          char **argv /**< [in] standard command-line argument vector */ )
{
    return p4ReductionPrecisionBenchmarkMain( argc, argv );
}
