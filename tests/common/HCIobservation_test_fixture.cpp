/** \file HCIobservation_test_fixture.cpp
 * \brief Implements shared temporary-input helpers for HCIobservation tests.
 * \author Jared R. Males
 */

#include "HCIobservation_test_fixture.hpp"

#include <atomic>
#include <fstream>
#include <stdexcept>
#include <system_error>

#include <omp.h>
#include <unistd.h>

/// \cond HCIobservation_test_harness
namespace unitTest
{
namespace
{
std::atomic<unsigned long> testDirectorySequence{ 0 };
}

void HCIobservationTestHarness::postReadFiles()
{
    m_hookEvents.push_back( "target-read" );
}

void HCIobservationTestHarness::postCoadd()
{
    m_hookEvents.push_back( "target-coadd" );
}

void HCIobservationTestHarness::postRDIReadFiles()
{
    m_hookEvents.push_back( "rdi-read" );
}

void HCIobservationTestHarness::postRDICoadd()
{
    m_hookEvents.push_back( "rdi-coadd" );
}

TestDirectory::TestDirectory()
{
    std::error_code error;

    for( unsigned int attempt = 0; attempt < 100; ++attempt )
    {
        const unsigned long sequence = testDirectorySequence.fetch_add( 1 );
        m_path = std::filesystem::temp_directory_path() /
                 ( "hcireduce-test-" + std::to_string( getpid() ) + "-" + std::to_string( sequence ) );

        if( std::filesystem::create_directory( m_path, error ) )
        {
            return;
        }

        if( error )
        {
            throw std::runtime_error( "could not create test directory: " + error.message() );
        }
    }

    throw std::runtime_error( "could not allocate a unique test directory" );
}

TestDirectory::~TestDirectory()
{
    std::error_code error;
    std::filesystem::remove_all( m_path, error );
}

const std::filesystem::path &TestDirectory::path() const
{
    return m_path;
}

std::filesystem::path TestDirectory::file( const std::string &name ) const
{
    return m_path / name;
}

OpenMPThreadGuard::OpenMPThreadGuard( int threads ) : m_previousThreads( omp_get_max_threads() )
{
    if( threads <= 0 )
    {
        throw std::invalid_argument( "OpenMP thread count must be positive" );
    }

    omp_set_num_threads( threads );
}

OpenMPThreadGuard::~OpenMPThreadGuard()
{
    omp_set_num_threads( m_previousThreads );
}

void writeTextFile( const std::filesystem::path &path, const std::string &contents )
{
    std::ofstream output( path );
    if( !output )
    {
        throw std::runtime_error( "could not open test file for writing: " + path.string() );
    }

    output << contents;
    if( !output )
    {
        throw std::runtime_error( "could not write test file: " + path.string() );
    }
}

void writeFitsImage( const std::filesystem::path &path,
                     const HCIobservationTestHarness::imageT &image,
                     HCIobservationTestHarness::fitsHeaderT *header )
{
    HCIobservationTestHarness::fitsFileT fitsFile;
    mx::error_t error = mx::error_t::noerror;

    if( header == nullptr )
    {
        error = fitsFile.write( path.string(), image );
    }
    else
    {
        error = fitsFile.write( path.string(), image, *header );
    }

    if( error != mx::error_t::noerror )
    {
        throw mx::exception<mx::verbose::vv>( error, "writing test FITS image" );
    }
}

void writeFitsCube( const std::filesystem::path &path,
                    const HCIobservationTestHarness::cubeT &cube,
                    HCIobservationTestHarness::fitsHeaderT *header )
{
    HCIobservationTestHarness::fitsFileT fitsFile;
    mx::error_t error = mx::error_t::noerror;

    if( header == nullptr )
    {
        error = fitsFile.write( path.string(), cube );
    }
    else
    {
        error = fitsFile.write( path.string(), cube, *header );
    }

    if( error != mx::error_t::noerror )
    {
        throw mx::exception<mx::verbose::vv>( error, "writing test FITS cube" );
    }
}

} // namespace unitTest
/// \endcond
