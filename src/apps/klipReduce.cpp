
/** \file klipReduce.cpp
 * \brief Defines the klipReduce command-line application.
 * \author Jared R. Males
 */

#include <iostream>

#include <mx/app/application.hpp>
using namespace mx::app;

#include "../common/ADIDerotator.hpp"
#include "../common/KLIPreduction.hpp"
using namespace mx::improc;

//#include <libgen.h>

/// A program to run the KLIP pipeline
/**
 *
 */
template <typename _realT, typename _evCalcT, class _verboseT>
class klipReduce : public application
{
  public:
    typedef _realT realT;
    typedef _evCalcT evCalcT;
    typedef _verboseT verboseT;

  protected:
    //************************************//
    // Mode of execution                  //
    std::string m_mode{ "basic" };                 ///< Execution mode: basic, normal, or postprocess.

    std::string m_postprocessDirectory;            ///< Directory containing saved PSF-subtracted products.
    std::string m_postprocessPrefix;               ///< Literal prefix of saved PSF-subtracted products.
    std::string m_postprocessExtension{ ".fits" }; ///< Literal extension of saved PSF-subtracted products.

    KLIPreduction<realT, ADIDerotator<realT, verboseT>, evCalcT, verboseT> m_obs;

  public:
    klipReduce()
    {
        m_configPathGlobal_env = "KLIPREDUCE_GLOBAL_CONFIG";
        m_configPathLocal = "klipReduce.conf";
        m_requireConfigPathLocal = false;

        config.m_sources = true;
    }

    ~klipReduce()
    {
    }

    /// Register application and reduction configuration targets.
    void setupConfig()
    {
        config.add( "mode",
                    "",
                    "mode",
                    argType::Required,
                    "",
                    "mode",
                    false,
                    "string",
                    "The mode of operation: basic/normal (the default) or postprocess" );

        config.add( "postprocess.directory",
                    "",
                    "postprocess.directory",
                    argType::Required,
                    "postprocess",
                    "directory",
                    false,
                    "string",
                    "Directory containing saved PSF-subtracted FITS products" );

        config.add( "postprocess.prefix",
                    "",
                    "postprocess.prefix",
                    argType::Required,
                    "postprocess",
                    "prefix",
                    false,
                    "string",
                    "Literal prefix of saved PSF-subtracted FITS products" );

        config.add( "postprocess.extension",
                    "",
                    "postprocess.extension",
                    argType::Required,
                    "postprocess",
                    "extension",
                    false,
                    "string",
                    "Literal extension of saved PSF-subtracted FITS products (default: .fits)" );

        m_obs.setupConfig( config );
    }

    /// Load the command-line configuration path.
    virtual void setConfigPathCL()
    {
        config.get<std::string>( m_configPathCL, "config" );
    }

    /// Load application and reduction configuration values.
    void loadConfig()
    {
        m_obs.loadConfig( config );

        config( m_mode, "mode" );
        config( m_postprocessDirectory, "postprocess.directory" );
        config( m_postprocessPrefix, "postprocess.prefix" );
        config( m_postprocessExtension, "postprocess.extension" );

        // This checks for unused config options, printing the banner only once no matter how many there are.
        // This will catch both bad options, and options we aren't actually using (debugging).
        bool unusedPrinted = false;
        for( auto it = config.m_targets.begin(); it != config.m_targets.end(); ++it )
        {
            if( it->second.used == false )
            {
                if( !unusedPrinted )
                {
                    std::cerr << "****************************************************\n";
                    std::cerr << "WARNING: unused config options (this is a programmer error):\n";
                    unusedPrinted = true;
                }

                std::cerr << "   " << it->second.name << '\n';
            }
        }

        unusedPrinted = false;
        if( config.m_unusedConfigs.size() > 0 )
        {
            if( !unusedPrinted )
            {
                std::cerr << "****************************************************\n";
                std::cerr << "WARNING: unrecognized config options:\n";
                unusedPrinted = true;
            }

            for( auto it = config.m_unusedConfigs.begin(); it != config.m_unusedConfigs.end(); ++it )
            {
                std::cerr << "   " << it->second.name;
                if( config.m_sources )
                    std::cerr << " [" << it->second.sources[0] << "]\n";
                else
                    std::cerr << "\n";
            }

            std::cerr << "****************************************************\n";
        }

        if( config.nonOptions.size() > 0 )
        {
            std::cerr << "****************************************************\n";
            std::cerr << "WARNING: unrecognized command line arguments\n";
        }
    }

    /// Validate configuration required by the selected execution mode.
    void checkConfig()
    {
        if( m_mode != "basic" && m_mode != "normal" && m_mode != "postprocess" )
        {
            throw mx::exception<verboseT>( mx::error_t::invalidconfig, "invalid klipReduce mode: " + m_mode );
        }

        if( m_mode == "postprocess" )
        {
            if( m_postprocessDirectory.empty() || m_postprocessPrefix.empty() || m_postprocessExtension.empty() )
            {
                throw mx::exception<verboseT>(
                    mx::error_t::invalidconfig,
                    "postprocess mode requires postprocess.directory, postprocess.prefix, and postprocess.extension" );
            }
            return;
        }

        if( m_obs.preprocessingOnly() )
        {
            return;
        }

        // KLIP:

        if( m_mode != "postprocess" )
        {
            if( m_obs.m_Nmodes.size() == 0 )
            {
                std::cerr << invokedName << ": must specify number of modes (Nmodes)\n";
            }

            if( m_obs.m_minRadius.size() == 0 )
            {
                std::cerr << invokedName << ": must specify minimum radii of KLIP regions (minRadius)\n";
            }

            if( m_obs.m_maxRadius.size() == 0 )
            {
                std::cerr << invokedName << ": must specify maximum radii of KLIP regions (maxRadius)\n";
            }

            if( m_obs.m_minRadius.size() != m_obs.m_maxRadius.size() )
            {
                std::cerr << invokedName << ": number of minimum and maximum radii must be equal\n";
            }

            if( m_obs.m_minAngle.size() != m_obs.m_maxAngle.size() )
            {
                std::cerr << invokedName << ": number of minimum and maximum angles must be equal\n";
            }
        }

        return;
    }

    /// Execute the selected reduction mode.
    virtual int execute()
    {

        if( m_mode == "postprocess" )
        {
            return m_obs.processPSFSub( m_postprocessDirectory, m_postprocessPrefix, m_postprocessExtension );
        }
        else
        {
            /*
            if( nWedges > 0 )
            {
                if( minRadius.size() > 1 || maxRadius.size() > 1 )
                {
                    std::cerr << "Error: nWedges set but min/maxRadius have more than one entry" << "\n";
                    return -1;
                }

                if( minAngle.size() > 0 || maxAngle.size() > 0 )
                {
                    std::cerr << "Warning: nWedges set but min/maxAngle have more than one entry" << "\n";
                    minAngle.clear();
                    maxAngle.clear();
                }

                if( ( 360 % nWedges ) != 0 )
                {
                    std::cerr << "Error: nWedges must be a divisor of 360\n";
                    return -1;
                }

                realT mnr = minRadius[0];
                realT mxr = maxRadius[0];

                minRadius.resize( nWedges, mnr );
                maxRadius.resize( nWedges, mxr );

                int dang = 360 / nWedges;
                minAngle.resize( nWedges );
                maxAngle.resize( nWedges );

                for( size_t n = 0; n < minAngle.size(); ++n )
                {
                    minAngle[n] = n * dang;
                    maxAngle[n] = n * dang + dang;
                }
            }
            else*/
            {
                if( m_obs.m_minAngle.size() == 0 || m_obs.m_maxAngle.size() == 0 )
                {
                    m_obs.m_minAngle.resize( m_obs.m_minRadius.size(), 0 );
                    m_obs.m_maxAngle.resize( m_obs.m_maxRadius.size(), 360 );
                }
            }

            m_obs.load_fileList();

            // if(m_obs.m_RDIfileListFile != "" ||  m_obs.m_RDIdirectory != "")
            //{
            m_obs.load_RDIfileList();
            //}

            return m_obs.regions( m_obs.m_minRadius, m_obs.m_maxRadius, m_obs.m_minAngle, m_obs.m_maxAngle );
        }
    }
};

#ifndef HCIREDUCE_KLIPREDUCE_NO_MAIN
int main( int argc, char **argv )
{

    std::string argv0 = argv[0];

    klipReduce<float, double, mx::verbose::vv> kr;

    try
    {
        kr.main( argc, argv );
    }
    catch( const std::exception &e )
    {
        std::vector<std::string> whats;
        mx::unwind_exceptions(whats, e);

        mx::print_exceptions(whats, std::format("{}: exception(s) encountered during execution", argv0));
        std::cerr << std::format("\nTo get help try: {} -h\n\n", argv0);
        return -1;
    }

    return 0;
}
#endif
