
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
    std::string mode {"basic"}; // Choices are basic [default] and grid [grid manages a fake planet grid]

    // Executes a grid of fake planets.
    realT gridCenterSep {0};              ///< The separation of the grid center [pixels].
    realT gridCenterPA {0};               ///< The PA of the grid center [deg E of N].
    realT gridHalfWidthSep {0};           ///< The grid half-width in radius [pixels]
    realT gridDeltaSep {0};               ///< The grid spacing in radius [pixels]
    realT gridHalfWidthPA {0};            ///< The grid half-wdith in PA [pixels]
    realT gridDeltaPA {0};                ///< The grid spacing in PA [pixels]
    std::vector<realT> gridContrasts; ///< The grid contrasts, possibly negative.

    int doGrid();

    KLIPreduction<realT, ADIDerotator<realT,verboseT>, evCalcT, verboseT> m_obs;

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

    // This sets up the configuration
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
                    "The mode of operation: either \"grid\" or \"normal\" (the default)" );

        m_obs.setupConfig( config );

        config.add( "grid.centerSep",
                    "",
                    "grid.centerSep",
                    argType::Required,
                    "grid",
                    "centerSep",
                    false,
                    "float",
                    "The grid center in separation [pixels]" );

        config.add( "grid.centerPA",
                    "",
                    "grid.centerPA",
                    argType::Required,
                    "grid",
                    "centerPA",
                    false,
                    "float",
                    "The grid center in position angle [degrees]" );

        config.add( "grid.halfWidthSep",
                    "",
                    "grid.halfWidthSep",
                    argType::Required,
                    "grid",
                    "halfWidthSep",
                    false,
                    "float",
                    "The half width of the grid in spearation [pixels]" );

        config.add( "grid.halfWidthPA",
                    "",
                    "grid.halfWidthPA",
                    argType::Required,
                    "grid",
                    "halfWidthPA",
                    false,
                    "float",
                    "The half width of the grid in PA [degrees]" );

        config.add( "grid.deltaSep",
                    "",
                    "grid.deltaSep",
                    argType::Required,
                    "grid",
                    "deltaSep",
                    false,
                    "float",
                    "The grid step size in separation [pixels]" );

        config.add( "grid.deltaPA",
                    "",
                    "grid.deltaPA",
                    argType::Required,
                    "grid",
                    "deltaPA",
                    false,
                    "float",
                    "The grid step size in PA [degrees]" );

        config.add( "grid.contrasts",
                    "",
                    "grid.contrasts",
                    argType::Required,
                    "grid",
                    "contrasts",
                    false,
                    "vector<float>",
                    "The contrast grid [planet:star]." );

    }

    virtual void setConfigPathCL()
    {
        config.get<std::string>( m_configPathCL, "config" );
    }

    void loadConfig()
    {
        m_obs.loadConfig( config );

        config( gridCenterSep, "grid.centerSep" );

        config( gridCenterPA, "grid.centerPA" );
        config( gridHalfWidthSep, "grid.halfWidthSep" );
        config( gridDeltaSep, "grid.deltaSep" );
        config( gridHalfWidthPA, "grid.halfWidthPA" );
        config( gridDeltaPA, "grid.deltaPA" );
        config( gridContrasts, "grid.contrasts" );

        config( mode, "mode" );

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

    void checkConfig()
    {

        // KLIP:

        if( mode != "postprocess" )
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

    virtual int execute()
    {

        if( mode == "grid" )
        {
            return doGrid();
        }
        /*else if( mode == "postprocess" )
        {
            return m_obs.processPSFSub( directory, prefix, extension );
        }*/
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

template <typename realT, typename evCalcT, class verboseT>
int klipReduce<realT, evCalcT, verboseT>::doGrid()
{
    if( gridCenterSep == 0 )
    {
        mx::error_report<verboseT>( mx::error_t::paramnotset, "Grid center separation not set (grid.centerSep)" );
        return -1;
    }

    if( gridCenterPA == 0 )
    {
        mx::error_report<verboseT>( mx::error_t::paramnotset, "Grid center PA not set (grid.centerPA)" );
        return -1;
    }

    if( gridHalfWidthSep == 0 )
    {
        mx::error_report<verboseT>( mx::error_t::paramnotset, "Grid half-width in radius not set (grid.halfWidthSep)" );
        return -1;
    }

    if( gridHalfWidthPA == 0 )
    {
        mx::error_report<verboseT>( mx::error_t::paramnotset, "Grid half-width in PA not set (grid.halfWidthPA)" );
        return -1;
    }

    if( gridDeltaSep == 0 )
    {
        mx::error_report<verboseT>( mx::error_t::paramnotset, "Grid spacing in radius not set (gridDeltaSep)" );
        return -1;
    }

    if( gridDeltaPA == 0 )
    {
        mx::error_report<verboseT>( mx::error_t::paramnotset, "Grid spacing in PA not set (grid.deltaPA)" );
        return -1;
    }

    if( gridContrasts.size() == 0 )
    {
        mx::error_report<verboseT>( mx::error_t::paramnotset, "Grid contrasts not set (grid.contrasts)" );
        return -1;
    }

    realT x0, y0;

    x0 = -1 * gridCenterSep * sin( mx::math::dtor( gridCenterPA ) );
    y0 = gridCenterSep * cos( mx::math::dtor( gridCenterPA ) );

    //   std::cerr << gridCenterSep << " " << x0 << " " << y0 << "\n";

    int Nrad = 2 * floor( gridHalfWidthSep / gridDeltaSep ) + 1;
    int Npa = 2 * floor( gridHalfWidthPA / gridDeltaPA ) + 1;

    Eigen::Array<realT, -1, -1> sep, pa;

    sep.resize( Nrad, Npa );
    pa.resize( Nrad, Npa );

    realT xp, yp, q, x, y;

    std::vector<realT> xv, yv;

    for( int i = 0; i < Nrad; ++i )
    {
        xp = ( -0.5 * ( Nrad - 1 ) + i ) * gridDeltaSep;

        for( int j = 0; j < Npa; ++j )
        {
            yp = ( -0.5 * ( Npa - 1 ) + j ) * gridDeltaPA;

            q = mx::math::dtor( 90 - gridCenterPA );

            x = ( x0 + xp * cos( q ) + yp * sin( q ) );
            y = ( y0 - xp * sin( q ) + yp * cos( q ) );

            xv.push_back( x );
            yv.push_back( y );

            sep( i, j ) = sqrt( pow( x, 2 ) + pow( y, 2 ) );

            std::cerr << "THIS WON'T WORK UNTIL YOU FIX ANGLEMOD\n";
            exit( 0 );
            // pa(i,j) = mx::math::angleMod<mx::math::degreesT<realT>>(mx::math::rtod( atan2(y, x))  - 90.0);

            // std::cerr << sep(i,j) << " " << pa(i,j) << "\n";

            for( size_t k = 0; k < gridContrasts.size(); ++k )
            {
                // m_obs.m_filesRead = false;
                std::cerr << "THIS WON'T WORK UNTIL YOU FIX FILESREAD IS PROTECTED\n";
                exit( 0 );

                m_obs.m_fakeSep = { sep( i, j ) };
                m_obs.m_fakePA = { pa( i, j ) };
                m_obs.m_fakeContrast = { gridContrasts[k] };

                // std::cerr << sep(i,j) << " " << pa(i,j) << " " << gridContrasts[k] << "\n";
                std::vector<realT> minMaxQ( m_obs.m_minRadius.size(), 0 );
                m_obs.regions( m_obs.m_minRadius, m_obs.m_maxRadius, minMaxQ, minMaxQ );
            }
        }
    }

    mx::fits::fitsFile<realT> ff;

    std::string fn;
    fn = "gridSep.fits";
    if( m_obs.m_outputDir != "" )
    {
        fn = m_obs.m_outputDir + "/" + fn;
    }

    ff.write( fn, sep );

    fn = "gridPA.fits";
    if( m_obs.m_outputDir != "" )
    {
        fn = m_obs.m_outputDir + "/" + fn;
    }
    ff.write( fn, pa );

    fn = "gridContrasts.dat";
    if( m_obs.m_outputDir != "" )
    {
        fn = m_obs.m_outputDir + "/" + fn;
    }
    std::ofstream fout;
    fout.open( fn );
    for( size_t i = 0; i < gridContrasts.size(); ++i )
    {
        fout << gridContrasts[i] << "\n";
    }
    fout.close();

    return 0;
}

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
