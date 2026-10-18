/** \file contrast_curve_test.cpp
 * \brief Tests the contrast-curve Python command line and its reduction workflow.
 */

#include "../catch2/catch.hpp"

#include <filesystem>
#include <sys/wait.h>
#include <unistd.h>

namespace unitTest
{
namespace contrast_curve_test
{

/// \cond contrast_curve_test_harness
/// Run one behavioral check against the real contrast_curve.py module and command line.
int runCheck( const char *name /**< [in] check selected by the Python driver */ )
{
    const auto driver =
        std::filesystem::path( __FILE__ ).parent_path().parent_path() / "scripts/contrast_curve_checks.py";
    const pid_t child = fork();
    if( child < 0 )
    {
        return -1;
    }
    if( child == 0 )
    {
        execlp( "python3", "python3", driver.c_str(), name, static_cast<char *>( nullptr ) );
        _exit( 127 );
    }
    int status = 0;
    if( waitpid( child, &status, 0 ) != child || !WIFEXITED( status ) )
    {
        return -1;
    }
    return WEXITSTATUS( status );
}
/// \endcond

/** \brief Verify contrast_curve.py sampling, contrast interpolation, configuration precedence, and SNR parsing. */
TEST_CASE( "Contrast curve sampling and parsing", "[contrast_curve]" )
{
    REQUIRE( runCheck( "grid" ) == 0 );

    // clang-format off
#ifdef __DOXY_ONLY__
    contrast_curve::samples( minimum, maximum, step );
    contrast_curve::contrasts( args, radii );
    contrast_curve::make_plan( args );
    contrast_curve::parse_snr( text, trial );
#endif
    // clang-format on
}

/** \brief Verify contrast_curve.py runs isolated injections, averages signed SNRs by mode, and resumes safely. */
TEST_CASE( "Contrast curve reduction workflow", "[contrast_curve]" )
{
    REQUIRE( runCheck( "workflow" ) == 0 );

    // clang-format off
#ifdef __DOXY_ONLY__
    contrast_curve::main( argv );
    contrast_curve::run_trial( plan, root, trial, resume );
    contrast_curve::summarize( plan, results );
#endif
    // clang-format on
}

/** \brief Verify contrast_curve.py records failures, reuses completed reductions, and counts NaN aperture SNRs
 * separately by mode when reporting incomplete PA coverage.
 */
TEST_CASE( "Contrast curve failed and incomplete measurements", "[contrast_curve]" )
{
    REQUIRE( runCheck( "failure" ) == 0 );

    // clang-format off
#ifdef __DOXY_ONLY__
    contrast_curve::run_trial( plan, root, trial, resume );
    contrast_curve::summarize( plan, results );
#endif
    // clang-format on
}

/** \brief Verify contrast_curve.py run_batch() uses one P4 process, preserves completed products after failure,
 * and resumes missing reductions and analyses independently.
 */
TEST_CASE( "Contrast curve P4 baseline reuse and batch recovery", "[contrast_curve][batch]" )
{
    REQUIRE( runCheck( "batch" ) == 0 );

    // clang-format off
#ifdef __DOXY_ONLY__
    contrast_curve::batch_command( plan, root, keep_going );
    contrast_curve::batch_row( root, trial );
    contrast_curve::run_batch( plan, root, trials, resume, keep_going );
    contrast_curve::run_trial( plan, root, trial, resume, batch_error );
#endif
    // clang-format on
}

} // namespace contrast_curve_test
} // namespace unitTest
