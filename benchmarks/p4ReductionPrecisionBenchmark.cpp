/** \file p4ReductionPrecisionBenchmark.cpp
 * \brief Runs a configured production-layout P4 reduction with an experimental precision policy.
 * \author Jared R. Males
 */

#include <chrono>
#include <cstdlib>
#include <format>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#define HCIREDUCE_P4REDUCE_NO_MAIN
#include "src/apps/p4Reduce.cpp"

/// Benchmark-only p4Reduce adapter selecting one internal P4 PCA precision policy.
class p4ReductionPrecisionBenchmark final : public p4Reduce
{
  protected:
    std::string m_precisionName; ///< Required frozen candidate identifier supplied by the benchmark caller.

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
}

void p4ReductionPrecisionBenchmark::loadConfig()
{
    config( m_precisionName, "benchmark.precision" );
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
}

int p4ReductionPrecisionBenchmark::execute()
{
    const auto begin = std::chrono::steady_clock::now();
    const int result = mx::improc::detail::p4ReductionReduceExperimental( m_obs, m_precisionPolicy );
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
