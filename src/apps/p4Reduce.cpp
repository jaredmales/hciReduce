/** \file p4Reduce.cpp
 * \brief Defines the p4Reduce command-line application.
 * \author Jared R. Males
 */

#include <cstdlib>
#include <format>
#include <iostream>
#include <string>
#include <vector>

#include <mx/app/application.hpp>

#include "../common/ADIDerotator.hpp"
#include "../common/P4Reduction.hpp"

/// Command-line application for Pixel Prediction Post-Processing.
/** The initial application supports one configured target-only P4 reduction. `basic` and `normal` are exact aliases;
 * saved-product postprocessing remains unavailable until the P4 validity-checkpoint format is defined.
 */
class p4Reduce : public mx::app::application
{
  protected:
    std::string m_mode{ "basic" };  ///< Execution mode: basic or normal.

    mx::improc::P4Reductionf m_obs; ///< Configured P4 observation and reduction state.

  public:
    /// Construct the application with its standard configuration paths.
    p4Reduce()
    {
        m_configPathGlobal_env = "P4REDUCE_GLOBAL_CONFIG";
        m_configPathLocal = "p4Reduce.conf";
        m_requireConfigPathLocal = false;
        config.m_sources = true;
    }

    /// Register application dispatch and P4 reduction configuration targets.
    void setupConfig() override
    {
        config.add( "mode",
                    "",
                    "mode",
                    mx::app::argType::Required,
                    "",
                    "mode",
                    false,
                    "string",
                    "The mode of operation: basic or normal (exact aliases; default basic)" );

        m_obs.setupConfig( config );
    }

    /// Load application and P4 reduction configuration values and report unused inputs.
    void loadConfig() override
    {
        m_obs.loadConfig( config );
        config( m_mode, "mode" );

        bool unusedPrinted{ false };
        for( const auto &[name, target] : config.m_targets )
        {
            if( target.used )
            {
                continue;
            }

            if( !unusedPrinted )
            {
                std::cerr << "****************************************************\n";
                std::cerr << "WARNING: unused config options (this is a programmer error):\n";
                unusedPrinted = true;
            }
            std::cerr << "   " << name << '\n';
        }

        if( !config.m_unusedConfigs.empty() )
        {
            std::cerr << "****************************************************\n";
            std::cerr << "WARNING: unrecognized config options:\n";
            for( const auto &[name, target] : config.m_unusedConfigs )
            {
                std::cerr << "   " << name;
                if( config.m_sources && !target.sources.empty() )
                {
                    std::cerr << " [" << target.sources.front() << ']';
                }
                std::cerr << '\n';
            }
            std::cerr << "****************************************************\n";
        }

        if( !config.nonOptions.empty() )
        {
            std::cerr << "****************************************************\n";
            std::cerr << "WARNING: unrecognized command line arguments\n";
        }
    }

    /// Validate application-level dispatch without duplicating P4 scientific validation.
    void checkConfig() override
    {
        if( m_mode != "basic" && m_mode != "normal" )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                  "invalid p4Reduce mode: " + m_mode +
                                                      "; supported modes are basic and normal" );
        }
    }

    /// Execute one fully configured P4 reduction.
    /** \returns zero on success, including a configured preprocessing-only stop.
     * \throws mx::exception when configuration, input, reduction, or output fails.
     */
    int execute() override
    {
        return m_obs.reduce();
    }
};

/// Run p4Reduce with command-line exception reporting.
/** The application framework's help status is normalized to `EXIT_SUCCESS`; any reduction status is otherwise
 * propagated. A standard exception is unwound, printed with a help hint, and returned as `EXIT_FAILURE`.
 *
 * \returns `EXIT_SUCCESS` on help or successful execution, otherwise `EXIT_FAILURE`.
 */
int p4ReduceMain( int argc, /**< [in] standard command-line argument count */
                  char **argv /**< [in] standard command-line argument vector */ )
{
    const std::string argv0 = argc > 0 && argv != nullptr && argv[0] != nullptr ? argv[0] : "p4Reduce";
    p4Reduce reduction;

    try
    {
        const int result = reduction.main( argc, argv );
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

#ifndef HCIREDUCE_P4REDUCE_NO_MAIN
/// p4Reduce process entry point.
int main( int argc, /**< [in] standard command-line argument count */
          char **argv /**< [in] standard command-line argument vector */ )
{
    return p4ReduceMain( argc, argv );
}
#endif
