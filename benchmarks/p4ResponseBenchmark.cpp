/** \file p4ResponseBenchmark.cpp
 * \brief Replays captured direct-P4 inputs through analytic and FP64 paired-refit responses.
 */

#include "src/common/P4PCA.hpp"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace
{

/// Parse a complete numeric argument, rejecting trailing characters.
template <typename valueT>
valueT number( const std::string &text /**< [in] numeric argument */ )
{
    std::istringstream stream( text );
    valueT value;
    if( !( stream >> value ) || stream.peek() != std::char_traits<char>::eof() )
    {
        throw std::invalid_argument( "invalid numeric argument: " + text );
    }
    return value;
}

/// Read native-endian FP64 ingress arrays from a capture with a validated byte count.
void readInputs( mx::improc::P4PCA::matrixT &predictors, /**< [out] pre-sized predictor array */
                 mx::improc::P4PCA::vectorT &target,     /**< [out] pre-sized target vector */
                 const std::string &path,                /**< [in] capture binary path */
                 std::size_t expectedBytes /**< [in] size including the unused captured residuals */ )
{
    if( std::filesystem::file_size( path ) != expectedBytes )
    {
        throw std::runtime_error( "incorrect capture size: " + path );
    }
    std::ifstream stream( path, std::ios::binary );
    stream.exceptions( std::ios::failbit | std::ios::badbit );
    stream.read( reinterpret_cast<char *>( predictors.data() ), predictors.size() * sizeof( double ) );
    stream.read( reinterpret_cast<char *>( target.data() ), target.size() * sizeof( double ) );
}

/// Serialize the public per-mode response status without relying on enum ordinals.
const char *statusName( mx::improc::P4PCAResponseStatus status /**< [in] analytic availability */ )
{
    using statusT = mx::improc::P4PCAResponseStatus;
    switch( status )
    {
        case statusT::differentiable:
            return "differentiable";
        case statusT::rankInsufficient:
            return "rankInsufficient";
        case statusT::rankBoundary:
            return "rankBoundary";
        case statusT::cutoffUnresolved:
            return "cutoffUnresolved";
    }
    throw std::logic_error( "unknown response status" );
}

} // namespace

/// Compare one captured baseline/unit-source pair, writing column-major native FP64 responses and JSON diagnostics.
int main( int argc, /**< [in] argument count */
          char **argv /**< [in] dimensions, tolerances, modes, and capture/output paths */ )
{
    try
    {
        if( argc != 9 )
        {
            throw std::invalid_argument(
                "usage: p4ResponseBenchmark ROWS COLUMNS RANK_TOL MODES_CSV EPS BASELINE.bin UNIT.bin OUTPUT.bin" );
        }
        const int rows = number<int>( argv[1] ), columns = number<int>( argv[2] );
        const double rankTolerance = number<double>( argv[3] ), amplitude = number<double>( argv[5] );
        std::vector<int> modes;
        std::istringstream list( argv[4] );
        std::string token;
        while( std::getline( list, token, ',' ) )
        {
            modes.push_back( number<int>( token ) );
        }
        if( rows <= 0 || columns <= 0 || modes.empty() || !std::isfinite( amplitude ) || amplitude <= 0 ||
            std::filesystem::exists( argv[8] ) )
        {
            throw std::invalid_argument( "require positive dimensions/amplitude, modes, and a new output file" );
        }
        const std::size_t width = static_cast<std::size_t>( columns ) + 1 + modes.size();
        if( width > std::numeric_limits<std::size_t>::max() / sizeof( double ) / rows )
        {
            throw std::length_error( "capture dimensions overflow byte count" );
        }
        const std::size_t bytes = width * rows * sizeof( double );
        using pcaT = mx::improc::P4PCA;
        pcaT::matrixT predictors( rows, columns ), source( rows, columns );
        pcaT::vectorT target( rows ), targetSource( rows );
        readInputs( predictors, target, argv[6], bytes );
        readInputs( source, targetSource, argv[7], bytes );
        source -= predictors;
        targetSource -= target;

        pcaT::workspaceT analyticWorkspace, refitWorkspace;
        mx::improc::P4PCAResponseResult analytic;
        mx::improc::P4PCAResult positive, negative;
        const auto begin = std::chrono::steady_clock::now();
        pcaT::calculateResponse( analytic,
                                 predictors,
                                 target,
                                 source,
                                 targetSource,
                                 modes,
                                 rankTolerance,
                                 analyticWorkspace );
        const auto middle = std::chrono::steady_clock::now();
        pcaT::calculate( positive,
                         predictors + amplitude * source,
                         target + amplitude * targetSource,
                         modes,
                         rankTolerance,
                         refitWorkspace );
        pcaT::calculate( negative,
                         predictors - amplitude * source,
                         target - amplitude * targetSource,
                         modes,
                         rankTolerance,
                         refitWorkspace );
        const pcaT::matrixT paired = ( positive.residuals - negative.residuals ) / ( 2 * amplitude );
        const auto end = std::chrono::steady_clock::now();
        std::ofstream output( argv[8], std::ios::binary );
        output.exceptions( std::ios::failbit | std::ios::badbit );
        output.write( reinterpret_cast<const char *>( analytic.responses.data() ),
                      analytic.responses.size() * sizeof( double ) );
        output.write( reinterpret_cast<const char *>( paired.data() ), paired.size() * sizeof( double ) );
        output.close();

        std::cout << std::setprecision( 17 ) << "{\"schema\":1,\"rank\":" << analytic.numericalRank
                  << ",\"relative_resolution\":" << analytic.relativeResolution
                  << ",\"analytic_seconds\":" << std::chrono::duration<double>( middle - begin ).count()
                  << ",\"paired_seconds\":" << std::chrono::duration<double>( end - middle ).count() << ",\"modes\":[";
        for( std::size_t mode = 0; mode < modes.size(); ++mode )
        {
            if( mode != 0 )
            {
                std::cout << ',';
            }
            std::cout << "{\"count\":" << modes[mode] << ",\"status\":\"" << statusName( analytic.modeStatus[mode] )
                      << "\",\"paired_supported\":"
                      << ( positive.sampleSupported( 0, mode ) && negative.sampleSupported( 0, mode ) ? "true"
                                                                                                      : "false" )
                      << ",\"relative_cutoff_gap\":";
            const double gap = analytic.relativeCutoffGaps[mode];
            if( std::isfinite( gap ) )
            {
                std::cout << gap;
            }
            else
            {
                std::cout << "null";
            }
            std::cout << '}';
        }
        std::cout << "]}\n";
    }
    catch( const std::exception &error )
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
