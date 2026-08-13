/** \file KLIPreduction_test.cpp
 * \brief Tests KLIP reduction centering and normalization behavior.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "HCIobservation_test_fixture.hpp"
#include "src/common/ADIDerotator.hpp"
#include "src/common/KLIPreduction.hpp"

namespace unitTest
{
namespace KLIPreduction_test
{

/// \cond KLIPreduction_test_harness
using reductionT =
    mx::improc::KLIPreduction<float, mx::improc::ADIDerotator<float, mx::verbose::vv>, double, mx::verbose::vv>;

void readReductionConfig( reductionT &reduction, const std::filesystem::path &path, const std::string &contents )
{
    writeTextFile( path, contents );
    mx::app::appConfigurator config;
    reduction.setupConfig( config );
    if( config.readConfig( path.string() ) != 0 )
    {
        throw std::runtime_error( "could not read KLIP configuration" );
    }
    reduction.loadConfig( config );
}

struct testDerotator
{
    std::vector<double> m_angles;

    double derotAngle( size_t index ) const
    {
        return m_angles.at( index );
    }
};

struct reductionHarness : public reductionT
{
    using reductionT::m_exactFinimName;
    using reductionT::m_filesRead;
    using reductionT::m_finim;
    using reductionT::m_finimName;
    using reductionT::m_heads;
    using reductionT::m_imageMJD;
    using reductionT::m_imSize;
    using reductionT::m_mask;
    using reductionT::m_maskFile;
    using reductionT::m_Ncols;
    using reductionT::m_Nims;
    using reductionT::m_Npix;
    using reductionT::m_Nrows;
    using reductionT::m_outputDir;
    using reductionT::m_psfsub;
    using reductionT::m_PSFSubPrefix;
    using reductionT::m_RDIfilesRead;
    using reductionT::m_RDIimageMJD;
    using reductionT::m_refIms;
    using reductionT::m_tgtIms;

    void postReadFiles() override
    {
    }
};

void prepareRegionReduction( reductionHarness &reduction )
{
    reduction.m_filesRead = true;
    reduction.m_RDIfilesRead = true;
    reduction.m_imSize = 3;
    reduction.m_Nrows = 3;
    reduction.m_Ncols = 3;
    reduction.m_Nims = 3;
    reduction.m_Npix = 9;
    reduction.m_tgtIms.resize( 3, 3, 3 );
    reduction.m_tgtIms.image( 0 ) << 1, 0, 2, -1, 3, 1, 2, -2, 0;
    reduction.m_tgtIms.image( 1 ) << 0, 2, -1, 1, -2, 3, -1, 1, 2;
    reduction.m_tgtIms.image( 2 ) << 2, -1, 0, 3, 1, -2, 0, 2, -1;
    reduction.m_Nmodes = { 1, 2 };
    reduction.m_meanSubMethod = mx::improc::HCI::meanSub::none;
    reduction.m_pixelTSNormMethod = mx::improc::HCI::pixelTSNorm::none;
    reduction.m_excludeMethod = mx::improc::HCI::exclude::none;
    reduction.m_excludeMethodMax = mx::improc::HCI::exclude::none;
    reduction.m_includeMethod = mx::improc::HCI::include::all;
    reduction.m_includeRefNum = 0;
    reduction.m_doDerotate = false;
    reduction.m_combineMethod = mx::improc::HCI::combine::none;
    reduction.m_doWriteFinim = false;
    reduction.m_doOutputPSFSub = false;
}
/// \endcond

/// Verify KLIPreduction registers and loads its diagnostic-output configuration with safe defaults.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP diagnostic configuration", "[KLIPreduction][config][diagnostics]" )
{
    reductionT defaults;
    mx::app::appConfigurator config;
    defaults.setupConfig( config );

    REQUIRE( config.m_targets.at( "klip.writeDiagnostics" ).clType == mx::app::argType::True );
    REQUIRE( config.m_targets.at( "klip.writeDiagnostics" ).helpType == "bool" );
    REQUIRE( config.m_targets.at( "klip.diagnosticDirectory" ).clType == mx::app::argType::Required );
    REQUIRE( config.m_targets.at( "klip.diagnosticDirectory" ).helpType == "string" );
    REQUIRE( config.m_targets.at( "klip.pixelTSNormMethod" ).helpType == "string" );
    REQUIRE( config.m_targets.at( "klip.pixelTSSigma" ).helpType == "float" );
    REQUIRE( config.m_targets.at( "klip.includeMethod" ).clType == mx::app::argType::Required );
    REQUIRE( config.m_targets.at( "klip.includeMethod" ).helpType == "string" );

    defaults.loadConfig( config );
    REQUIRE_FALSE( defaults.m_writeDiagnostics );
    REQUIRE( defaults.m_diagnosticDirectory == "." );
    REQUIRE( defaults.m_includeMethod == mx::improc::HCI::include::all );
    REQUIRE( defaults.m_includeRefNum == 0 );

    TestDirectory directory;
    reductionT configured;
    const std::filesystem::path diagnosticDirectory = directory.file( "nested/diagnostics" );
    readReductionConfig( configured,
                         directory.file( "klip.conf" ),
                         "[klip]\nwriteDiagnostics=true\ndiagnosticDirectory=" + diagnosticDirectory.string() + "\n" );

    REQUIRE( configured.m_writeDiagnostics );
    REQUIRE( configured.m_diagnosticDirectory == diagnosticDirectory.string() );
}

/// Verify KLIPreduction::loadConfig loads geometric, selection, centering, and normalization settings.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP configuration loading", "[KLIPreduction][config]" )
{
    TestDirectory directory;
    reductionT reduction;
    readReductionConfig( reduction,
                         directory.file( "klip.conf" ),
                         "[adi]\n"
                         "minDPx=1.5\n"
                         "maxDPx=7.5\n"
                         "excludeMethod=pixel\n"
                         "excludeMethodMax=angle\n"
                         "[geom]\n"
                         "minRadius=2,8\n"
                         "maxRadius=7,15\n"
                         "minAngle=-45,30\n"
                         "maxAngle=15,75\n"
                         "nWedges=6\n"
                         "[klip]\n"
                         "meanSubMethod=imageMedian\n"
                         "pixelTSNormMethod=rms\n"
                         "pixelTSSigma=2.25\n"
                         "includeRefNum=4\n"
                         "includeMethod=angle\n"
                         "Nmodes=1,3,5\n"
                         "rightReason=true\n"
                         "rrRadius=3.5\n" );

    REQUIRE( reduction.m_minDPx == Approx( 1.5 ) );
    REQUIRE( reduction.m_maxDPx == Approx( 7.5 ) );
    REQUIRE( reduction.m_excludeMethod == mx::improc::HCI::exclude::pixel );
    REQUIRE( reduction.m_excludeMethodMax == mx::improc::HCI::exclude::angle );
    REQUIRE( reduction.m_minRadius == std::vector<float>{ 2, 8 } );
    REQUIRE( reduction.m_maxRadius == std::vector<float>{ 7, 15 } );
    REQUIRE( reduction.m_minAngle == std::vector<float>{ -45, 30 } );
    REQUIRE( reduction.m_maxAngle == std::vector<float>{ 15, 75 } );
    REQUIRE( reduction.m_nWedges == 6 );
    REQUIRE( reduction.m_meanSubMethod == mx::improc::HCI::meanSub::imageMedian );
    REQUIRE( reduction.m_pixelTSNormMethod == mx::improc::HCI::pixelTSNorm::rms );
    REQUIRE( reduction.m_pixelTSSigma == Approx( 2.25 ) );
    REQUIRE( reduction.m_includeRefNum == 4 );
    REQUIRE( reduction.m_includeMethod == mx::improc::HCI::include::angle );
    REQUIRE( reduction.m_Nmodes == std::vector<int>{ 1, 3, 5 } );
    REQUIRE( reduction.m_rightReason );
    REQUIRE( reduction.m_rightReasonRadius == Approx( 3.5 ) );
}

/// Verify KLIPreduction::loadConfig accepts every documented reference-inclusion method.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP inclusion-method configuration", "[KLIPreduction][config][include]" )
{
    TestDirectory directory;
    const std::vector<std::pair<std::string, mx::improc::HCI::include>> methods{
        { "all", mx::improc::HCI::include::all },
        { "corr", mx::improc::HCI::include::corr },
        { "time", mx::improc::HCI::include::time },
        { "angle", mx::improc::HCI::include::angle },
        { "imno", mx::improc::HCI::include::imno },
    };

    for( const auto &[name, method] : methods )
    {
        const std::filesystem::path configPath = directory.file( "include-" + name + ".conf" );
        writeTextFile( configPath, "[klip]\nincludeMethod=" + name + "\n" );
        reductionT reduction;
        mx::app::appConfigurator config;
        reduction.setupConfig( config );
        REQUIRE( config.readConfig( configPath.string() ) == 0 );
        reduction.loadConfig( config );
        REQUIRE( reduction.m_includeMethod == method );
    }
}

/// Verify KLIPreduction::loadConfig rejects unsupported enum strings.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP invalid configuration methods", "[KLIPreduction][config][validation]" )
{
    TestDirectory directory;

    reductionT badMean;
    REQUIRE_THROWS(
        readReductionConfig( badMean, directory.file( "bad-mean.conf" ), "[klip]\nmeanSubMethod=not-a-method\n" ) );

    reductionT badNorm;
    REQUIRE_THROWS(
        readReductionConfig( badNorm, directory.file( "bad-norm.conf" ), "[klip]\npixelTSNormMethod=not-a-method\n" ) );

    reductionT badExclusion;
    REQUIRE_THROWS( readReductionConfig( badExclusion,
                                         directory.file( "bad-exclusion.conf" ),
                                         "[adi]\nexcludeMethod=not-a-method\n" ) );

    reductionT badInclusion;
    REQUIRE_THROWS( readReductionConfig( badInclusion,
                                         directory.file( "bad-inclusion.conf" ),
                                         "[klip]\nincludeMethod=not-a-method\n" ) );

    reductionT negativeReferenceCount;
    REQUIRE_THROWS( readReductionConfig( negativeReferenceCount,
                                         directory.file( "negative-reference-count.conf" ),
                                         "[klip]\nincludeRefNum=-1\n" ) );
}

/// Verify KLIPreduction::writeDiagnostic is inert by default and writes exact FITS data when enabled.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP diagnostic output gate", "[KLIPreduction][diagnostics]" )
{
    TestDirectory directory;
    reductionT reduction;
    reduction.m_diagnosticDirectory = directory.file( "nested/diagnostics" ).string();
    reductionT::imageT image( 2, 2 );
    image << 1, 2, 3, 4;

    reduction.writeDiagnostic( "cv.fits", image );
    REQUIRE_FALSE( std::filesystem::exists( directory.file( "nested/diagnostics/cv.fits" ) ) );

    reduction.m_writeDiagnostics = true;
    reduction.writeDiagnostic( "cv.fits", image );
    REQUIRE( std::filesystem::exists( directory.file( "nested/diagnostics/cv.fits" ) ) );

    reductionT::imageT loaded;
    mx::fits::fitsFile<float, mx::verbose::vv> reader;
    REQUIRE( reader.read( loaded, directory.file( "nested/diagnostics/cv.fits" ).string() ) == mx::error_t::noerror );
    REQUIRE( loaded.isApprox( image ) );
}

/// Verify KLIPreduction::meanSubtract leaves both cubes unchanged for `none` while calculating reference norms.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP mean subtraction none", "[KLIPreduction][meanSubtract][none]" )
{
    reductionT reduction;
    reduction.m_meanSubMethod = mx::improc::HCI::meanSub::none;
    reduction.m_pixelTSNormMethod = mx::improc::HCI::pixelTSNorm::none;

    mx::improc::eigenCube<float> references( 2, 2, 2 );
    references.image( 0 ) << 1, 3, 2, 4;
    references.image( 1 ) << 10, 14, 12, 16;
    mx::improc::eigenCube<float> targets( 2, 2, 1 );
    targets.image( 0 ) << 20, 22, 21, 23;
    const mx::improc::eigenCube<float> originalReferences = references;
    const mx::improc::eigenCube<float> originalTargets = targets;
    reductionT::imageT mask;
    std::vector<float> norms;

    reduction.meanSubtract( references, targets, mask, norms );

    REQUIRE( references.image( 0 ).isApprox( originalReferences.image( 0 ) ) );
    REQUIRE( references.image( 1 ).isApprox( originalReferences.image( 1 ) ) );
    REQUIRE( targets.image( 0 ).isApprox( originalTargets.image( 0 ) ) );
    REQUIRE( norms.size() == 2 );
    REQUIRE( norms[0] == Approx( std::sqrt( 5.0 ) ) );
    REQUIRE( norms[1] == Approx( std::sqrt( 20.0 ) ) );
}

/// Verify KLIPreduction::meanSubtract applies masked per-image mean and median centering to both libraries.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP masked per-image centering", "[KLIPreduction][meanSubtract][mask]" )
{
    for( const mx::improc::HCI::meanSub method :
         { mx::improc::HCI::meanSub::imageMean, mx::improc::HCI::meanSub::imageMedian } )
    {
        reductionT reduction;
        reduction.m_meanSubMethod = method;
        reduction.m_pixelTSNormMethod = mx::improc::HCI::pixelTSNorm::none;

        mx::improc::eigenCube<float> references( 1, 4, 2 );
        references.image( 0 ) << 1, 2, 100, 4;
        references.image( 1 ) << 5, 6, 200, 8;
        mx::improc::eigenCube<float> targets( 1, 4, 1 );
        targets.image( 0 ) << 9, 10, 300, 12;
        reductionT::imageT mask( 1, 4 );
        mask << 1, 1, 0, 1;
        std::vector<float> norms;

        reduction.meanSubtract( references, targets, mask, norms );

        if( method == mx::improc::HCI::meanSub::imageMean )
        {
            reductionT::imageT expected( 1, 4 );
            expected << -4.0 / 3.0, -1.0 / 3.0, 0, 5.0 / 3.0;
            REQUIRE( references.image( 0 ).isApprox( expected ) );
            REQUIRE( references.image( 1 ).isApprox( expected ) );
            REQUIRE( targets.image( 0 ).isApprox( expected ) );
        }
        else
        {
            reductionT::imageT expected( 1, 4 );
            expected << -1, 0, 0, 2;
            REQUIRE( references.image( 0 ).isApprox( expected ) );
            REQUIRE( references.image( 1 ).isApprox( expected ) );
            REQUIRE( targets.image( 0 ).isApprox( expected ) );
        }

        REQUIRE( norms.size() == 2 );
        REQUIRE( norms[0] == Approx( std::sqrt( 14.0 / 3.0 ) ) );
        REQUIRE( norms[1] == Approx( std::sqrt( 14.0 / 3.0 ) ) );
    }
}

/// Verify KLIPreduction::meanSubtract applies the reference mean or median image to distinct target data.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP reference-image centering", "[KLIPreduction][meanSubtract][reference]" )
{
    for( const mx::improc::HCI::meanSub method :
         { mx::improc::HCI::meanSub::meanImage, mx::improc::HCI::meanSub::medianImage } )
    {
        reductionT reduction;
        reduction.m_meanSubMethod = method;
        reduction.m_pixelTSNormMethod = mx::improc::HCI::pixelTSNorm::none;

        mx::improc::eigenCube<float> references( 1, 2, 3 );
        references.image( 0 ) << 1, 10;
        references.image( 1 ) << 3, 14;
        references.image( 2 ) << 9, 30;
        mx::improc::eigenCube<float> targets( 1, 2, 1 );
        targets.image( 0 ) << 13, 38;
        reductionT::imageT mask;
        std::vector<float> norms;

        reduction.meanSubtract( references, targets, mask, norms );

        reductionT::imageT expectedTarget( 1, 2 );
        if( method == mx::improc::HCI::meanSub::meanImage )
        {
            expectedTarget << 13.0 - 13.0 / 3.0, 20;
        }
        else
        {
            expectedTarget << 10, 24;
        }
        REQUIRE( targets.image( 0 ).isApprox( expectedTarget ) );

        reductionT::imageT referenceCenter;
        if( method == mx::improc::HCI::meanSub::meanImage )
        {
            references.mean( referenceCenter );
        }
        else
        {
            references.median( referenceCenter );
        }
        REQUIRE( referenceCenter.isZero( 1e-5 ) );
        REQUIRE( norms.size() == 3 );
    }
}

/// Verify KLIPreduction::meanSubtract uses true RMS pixel scaling, applies it to targets, and handles zero series.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP pixel time-series RMS normalization", "[KLIPreduction][meanSubtract][rms]" )
{
    reductionT reduction;
    reduction.m_meanSubMethod = mx::improc::HCI::meanSub::meanImage;
    reduction.m_pixelTSNormMethod = mx::improc::HCI::pixelTSNorm::rms;

    mx::improc::eigenCube<float> references( 1, 2, 2 );
    references.image( 0 ) << 1, 4;
    references.image( 1 ) << 3, 4;
    mx::improc::eigenCube<float> targets( 1, 2, 1 );
    targets.image( 0 ) << 6, 10;
    reductionT::imageT mask;
    std::vector<float> norms;

    reduction.meanSubtract( references, targets, mask, norms );

    reductionT::imageT expected0( 1, 2 );
    expected0 << -1, 0;
    reductionT::imageT expected1( 1, 2 );
    expected1 << 1, 0;
    reductionT::imageT expectedTarget( 1, 2 );
    expectedTarget << 4, 0;
    REQUIRE( references.image( 0 ).isApprox( expected0 ) );
    REQUIRE( references.image( 1 ).isApprox( expected1 ) );
    REQUIRE( targets.image( 0 ).isApprox( expectedTarget ) );
    REQUIRE( references.cube().isFinite().all() );
    REQUIRE( targets.cube().isFinite().all() );
}

/// Verify extractRowsAndCols and extractCols preserve index order and support empty selections.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP matrix extraction helpers", "[KLIPreduction][covariance][extract]" )
{
    Eigen::ArrayXXf matrix( 4, 4 );
    matrix << 0, 1, 2, 3, 10, 11, 12, 13, 20, 21, 22, 23, 30, 31, 32, 33;
    Eigen::ArrayXXf square;
    mx::improc::extractRowsAndCols( square, matrix, std::vector<size_t>{ 3, 1 } );
    Eigen::ArrayXXf expectedSquare( 2, 2 );
    expectedSquare << 33, 31, 13, 11;
    REQUIRE( square.isApprox( expectedSquare ) );

    Eigen::ArrayXXf vectors( 2, 4 );
    vectors << 0, 1, 2, 3, 10, 11, 12, 13;
    Eigen::ArrayXXf columns;
    mx::improc::extractCols( columns, vectors, std::vector<size_t>{ 2, 0 } );
    Eigen::ArrayXXf expectedColumns( 2, 2 );
    expectedColumns << 2, 0, 12, 10;
    REQUIRE( columns.isApprox( expectedColumns ) );

    mx::improc::extractRowsAndCols( square, matrix, {} );
    mx::improc::extractCols( columns, vectors, {} );
    REQUIRE( square.size() == 0 );
    REQUIRE( columns.rows() == 2 );
    REQUIRE( columns.cols() == 0 );
}

/// Verify collapseCovar applies image-number and angular exclusion boundaries.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP covariance exclusion", "[KLIPreduction][covariance][exclude]" )
{
    Eigen::ArrayXXf covariance = Eigen::MatrixXf::Identity( 4, 4 ).array();
    Eigen::ArrayXXf references( 2, 4 );
    references << 0, 1, 2, 3, 10, 11, 12, 13;
    const std::vector<float> norms( 4, 1 );
    testDerotator derotator{ { 0.0, 0.1, 0.5, 1.0 } };
    Eigen::ArrayXXf selectedCovariance;
    Eigen::ArrayXXf selectedReferences;
    mx::improc::eigenImage<int> included( 4, 4 );
    included.setZero();

    mx::improc::collapseCovar<float>( selectedCovariance,
                                      covariance,
                                      norms,
                                      selectedReferences,
                                      references,
                                      references,
                                      1,
                                      1,
                                      0,
                                      mx::improc::HCI::exclude::imno,
                                      mx::improc::HCI::exclude::none,
                                      mx::improc::HCI::include::all,
                                      0,
                                      derotator,
                                      derotator,
                                      {},
                                      {},
                                      included );
    REQUIRE( included( 1, 0 ) == 0 );
    REQUIRE( included( 1, 1 ) == 0 );
    REQUIRE( included( 1, 2 ) == 0 );
    REQUIRE( included( 1, 3 ) == 1 );
    REQUIRE( selectedCovariance.rows() == 1 );
    REQUIRE( selectedReferences.isApprox( references.col( 3 ) ) );

    included.setZero();
    mx::improc::collapseCovar<float>( selectedCovariance,
                                      covariance,
                                      norms,
                                      selectedReferences,
                                      references,
                                      references,
                                      1,
                                      0.15,
                                      0,
                                      mx::improc::HCI::exclude::angle,
                                      mx::improc::HCI::exclude::none,
                                      mx::improc::HCI::include::all,
                                      0,
                                      derotator,
                                      derotator,
                                      {},
                                      {},
                                      included );
    REQUIRE( included( 1, 0 ) == 0 );
    REQUIRE( included( 1, 1 ) == 0 );
    REQUIRE( included( 1, 2 ) == 1 );
    REQUIRE( included( 1, 3 ) == 1 );
    REQUIRE( selectedReferences.cols() == 2 );
    REQUIRE( selectedReferences.col( 0 ).isApprox( references.col( 2 ) ) );
    REQUIRE( selectedReferences.col( 1 ).isApprox( references.col( 3 ) ) );
}

/// Verify collapseCovar keeps an exact deterministic correlation-ranked maximum after exclusions.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP covariance correlation selection", "[KLIPreduction][covariance][include]" )
{
    Eigen::ArrayXXf covariance = Eigen::MatrixXf::Identity( 4, 4 ).array();
    Eigen::ArrayXXf references( 2, 4 );
    references << 1, 1, 0, -1, 0, 0, 1, 0;
    Eigen::ArrayXXf targets( 2, 1 );
    targets << 1, 0;
    const std::vector<float> norms( 4, 1 );
    testDerotator referenceDerotator{ { 0.0, 0.1, 0.5, 1.0 } };
    testDerotator targetDerotator{ { 0.1 } };
    Eigen::ArrayXXf selectedCovariance;
    Eigen::ArrayXXf selectedReferences;
    mx::improc::eigenImage<int> included( 1, 4 );
    included.setZero();

    mx::improc::collapseCovar<float>( selectedCovariance,
                                      covariance,
                                      norms,
                                      selectedReferences,
                                      references,
                                      targets,
                                      0,
                                      0,
                                      0,
                                      mx::improc::HCI::exclude::none,
                                      mx::improc::HCI::exclude::none,
                                      mx::improc::HCI::include::corr,
                                      1,
                                      referenceDerotator,
                                      targetDerotator,
                                      {},
                                      {},
                                      included );
    REQUIRE( included( 0, 0 ) == 1 );
    REQUIRE( included( 0, 1 ) == 0 );
    REQUIRE( included( 0, 2 ) == 0 );
    REQUIRE( included( 0, 3 ) == 0 );
    REQUIRE( selectedReferences.cols() == 1 );
    REQUIRE( selectedReferences.col( 0 ).isApprox( references.col( 0 ) ) );

    included.setZero();
    mx::improc::collapseCovar<float>( selectedCovariance,
                                      covariance,
                                      norms,
                                      selectedReferences,
                                      references,
                                      targets,
                                      0,
                                      0,
                                      0,
                                      mx::improc::HCI::exclude::imno,
                                      mx::improc::HCI::exclude::none,
                                      mx::improc::HCI::include::corr,
                                      1,
                                      referenceDerotator,
                                      targetDerotator,
                                      {},
                                      {},
                                      included );
    REQUIRE( included( 0, 0 ) == 0 );
    REQUIRE( included( 0, 1 ) == 1 );
    REQUIRE( included.row( 0 ).sum() == 1 );
    REQUIRE( selectedReferences.col( 0 ).isApprox( references.col( 1 ) ) );
}

/// Verify collapseCovar implements all, time, wrapped-angle, and image-number inclusion ordering.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP covariance inclusion methods", "[KLIPreduction][covariance][include]" )
{
    Eigen::ArrayXXf covariance = Eigen::MatrixXf::Identity( 4, 4 ).array();
    Eigen::ArrayXXf references( 2, 4 );
    references << 1, 2, 3, 4, 5, 6, 7, 8;
    Eigen::ArrayXXf targets = Eigen::ArrayXXf::Ones( 2, 2 );
    const std::vector<float> norms{ std::sqrt( 26.0F ), std::sqrt( 40.0F ), std::sqrt( 58.0F ), std::sqrt( 80.0F ) };
    testDerotator referenceDerotator{ { -3.1, 0.2, 2.0, -1.0 } };
    testDerotator targetDerotator{ { 3.1, 0.5 } };
    const std::vector<double> referenceMJD{ 10, 5, 12, 7 };
    const std::vector<double> targetMJD{ 8, 20 };
    Eigen::ArrayXXf selectedCovariance;
    Eigen::ArrayXXf selectedReferences;
    mx::improc::eigenImage<int> included( 2, 4 );

    included.setZero();
    mx::improc::collapseCovar<float>( selectedCovariance,
                                      covariance,
                                      norms,
                                      selectedReferences,
                                      references,
                                      targets,
                                      0,
                                      0,
                                      0,
                                      mx::improc::HCI::exclude::none,
                                      mx::improc::HCI::exclude::none,
                                      mx::improc::HCI::include::all,
                                      1,
                                      referenceDerotator,
                                      targetDerotator,
                                      referenceMJD,
                                      targetMJD,
                                      included );
    REQUIRE( included.row( 0 ).isOnes() );
    REQUIRE( selectedReferences.cols() == 4 );

    included.setZero();
    mx::improc::collapseCovar<float>( selectedCovariance,
                                      covariance,
                                      norms,
                                      selectedReferences,
                                      references,
                                      targets,
                                      0,
                                      0,
                                      0,
                                      mx::improc::HCI::exclude::none,
                                      mx::improc::HCI::exclude::none,
                                      mx::improc::HCI::include::time,
                                      0,
                                      referenceDerotator,
                                      targetDerotator,
                                      {},
                                      {},
                                      included );
    REQUIRE( included.row( 0 ).isOnes() );
    REQUIRE( selectedReferences.cols() == 4 );

    included.setZero();
    mx::improc::collapseCovar<float>( selectedCovariance,
                                      covariance,
                                      norms,
                                      selectedReferences,
                                      references,
                                      targets,
                                      0,
                                      0,
                                      0,
                                      mx::improc::HCI::exclude::none,
                                      mx::improc::HCI::exclude::none,
                                      mx::improc::HCI::include::time,
                                      2,
                                      referenceDerotator,
                                      targetDerotator,
                                      referenceMJD,
                                      targetMJD,
                                      included );
    REQUIRE( included( 0, 0 ) == 1 );
    REQUIRE( included( 0, 1 ) == 0 );
    REQUIRE( included( 0, 2 ) == 0 );
    REQUIRE( included( 0, 3 ) == 1 );

    included.setZero();
    mx::improc::collapseCovar<float>( selectedCovariance,
                                      covariance,
                                      norms,
                                      selectedReferences,
                                      references,
                                      targets,
                                      0,
                                      0,
                                      0,
                                      mx::improc::HCI::exclude::none,
                                      mx::improc::HCI::exclude::none,
                                      mx::improc::HCI::include::angle,
                                      1,
                                      referenceDerotator,
                                      targetDerotator,
                                      referenceMJD,
                                      targetMJD,
                                      included );
    REQUIRE( included( 0, 0 ) == 1 );
    REQUIRE( included.row( 0 ).sum() == 1 );

    included.setZero();
    mx::improc::collapseCovar<float>( selectedCovariance,
                                      covariance,
                                      norms,
                                      selectedReferences,
                                      references,
                                      targets,
                                      1,
                                      0,
                                      0,
                                      mx::improc::HCI::exclude::none,
                                      mx::improc::HCI::exclude::none,
                                      mx::improc::HCI::include::imno,
                                      2,
                                      referenceDerotator,
                                      targetDerotator,
                                      referenceMJD,
                                      targetMJD,
                                      included );
    REQUIRE( included( 1, 0 ) == 1 );
    REQUIRE( included( 1, 1 ) == 1 );
    REQUIRE( included( 1, 2 ) == 0 );
    REQUIRE( included( 1, 3 ) == 0 );
}

/// Verify collapseCovar rejects inconsistent covariance-selection dimensions.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP covariance input validation", "[KLIPreduction][covariance][validation]" )
{
    Eigen::ArrayXXf covariance = Eigen::MatrixXf::Identity( 2, 2 ).array();
    Eigen::ArrayXXf references = Eigen::ArrayXXf::Ones( 2, 2 );
    Eigen::ArrayXXf selectedCovariance;
    Eigen::ArrayXXf selectedReferences;
    mx::improc::eigenImage<int> included( 2, 2 );
    testDerotator derotator{ { 0.0, 0.1 } };

    REQUIRE_THROWS_AS( mx::improc::collapseCovar<float>( selectedCovariance,
                                                         covariance,
                                                         std::vector<float>{ 1 },
                                                         selectedReferences,
                                                         references,
                                                         references,
                                                         0,
                                                         0,
                                                         0,
                                                         mx::improc::HCI::exclude::none,
                                                         mx::improc::HCI::exclude::none,
                                                         mx::improc::HCI::include::all,
                                                         0,
                                                         derotator,
                                                         derotator,
                                                         {},
                                                         {},
                                                         included ),
                       std::invalid_argument );

    REQUIRE_THROWS_AS( mx::improc::collapseCovar<float>( selectedCovariance,
                                                         covariance,
                                                         std::vector<float>{ 1, 1 },
                                                         selectedReferences,
                                                         references,
                                                         references,
                                                         0,
                                                         0,
                                                         0,
                                                         mx::improc::HCI::exclude::none,
                                                         mx::improc::HCI::exclude::none,
                                                         mx::improc::HCI::include::time,
                                                         1,
                                                         derotator,
                                                         derotator,
                                                         {},
                                                         {},
                                                         included ),
                       std::invalid_argument );
}

/// Verify KLIPreduction::worker performs a deterministic one-mode subtraction and gates covariance diagnostics.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP worker one-mode reduction", "[KLIPreduction][worker][diagnostics]" )
{
    OpenMPThreadGuard threads( 1 );
    TestDirectory directory;
    reductionHarness reduction;
    reduction.m_meanSubMethod = mx::improc::HCI::meanSub::imageMean;
    reduction.m_pixelTSNormMethod = mx::improc::HCI::pixelTSNorm::none;
    reduction.m_includeMethod = mx::improc::HCI::include::all;
    reduction.m_includeRefNum = 1;
    reduction.m_Nrows = 1;
    reduction.m_Ncols = 3;
    reduction.m_Nims = 2;
    reduction.m_Nmodes = { 1 };
    reduction.m_maxNmodes = 1;
    reduction.m_psfsub.resize( 1 );
    reduction.m_psfsub[0].resize( 1, 3, 2 );
    reduction.m_psfsub[0].cube().setZero();
    reduction.m_imsIncluded.resize( 2, 2 );
    reduction.m_imsIncluded.setConstant( 1 );
    reduction.m_diagnosticDirectory = directory.file( "disabled" ).string();

    mx::improc::eigenCube<float> images( 3, 1, 2 );
    images.image( 0 ) << 1, -1, 0;
    images.image( 1 ) << 1, 1, -2;
    reductionT::imageT mask;
    std::vector<size_t> indices{ 0, 1, 2 };

    reduction.worker( images, images, mask, indices, 0, 0 );

    REQUIRE_FALSE( std::filesystem::exists( directory.file( "disabled/cv.fits" ) ) );
    REQUIRE( reduction.m_psfsub[0].image( 0 ).isApprox( ( reductionT::imageT( 1, 3 ) << 1, -1, 0 ).finished(), 1e-5 ) );
    REQUIRE( reduction.m_psfsub[0].image( 1 ).isZero( 1e-5 ) );

    reduction.m_writeDiagnostics = true;
    reduction.m_diagnosticDirectory = directory.file( "enabled" ).string();
    reduction.m_psfsub[0].cube().setZero();
    images.image( 0 ) << 1, -1, 0;
    images.image( 1 ) << 1, 1, -2;
    reduction.worker( images, images, mask, indices, 0, 0 );

    REQUIRE( std::filesystem::exists( directory.file( "enabled/cv.fits" ) ) );
    REQUIRE_FALSE( std::filesystem::exists( directory.file( "enabled/rrMask.fits" ) ) );
    reductionT::imageT covariance;
    mx::fits::fitsFile<float, mx::verbose::vv> reader;
    REQUIRE( reader.read( covariance, directory.file( "enabled/cv.fits" ).string() ) == mx::error_t::noerror );
    REQUIRE( covariance.rows() == 2 );
    REQUIRE( covariance.cols() == 2 );
    REQUIRE( covariance( 0, 0 ) == Approx( 2 ) );
    REQUIRE( covariance( 1, 0 ) == Approx( 0 ) );
    REQUIRE( covariance( 1, 1 ) == Approx( 6 ) );
}

/// Verify KLIPreduction::worker ranks an RDI library independently of the target image count.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP worker RDI correlation inclusion", "[KLIPreduction][worker][RDI][include]" )
{
    OpenMPThreadGuard threads( 1 );
    reductionHarness reduction;
    reduction.m_meanSubMethod = mx::improc::HCI::meanSub::none;
    reduction.m_pixelTSNormMethod = mx::improc::HCI::pixelTSNorm::none;
    reduction.m_includeMethod = mx::improc::HCI::include::corr;
    reduction.m_includeRefNum = 2;
    reduction.m_Nrows = 1;
    reduction.m_Ncols = 3;
    reduction.m_Nims = 1;
    reduction.m_Nmodes = { 1 };
    reduction.m_maxNmodes = 1;
    reduction.m_psfsub.resize( 1 );
    reduction.m_psfsub[0].resize( 1, 3, 1 );
    reduction.m_psfsub[0].cube().setZero();
    reduction.m_imsIncluded.resize( 1, 3 );
    reduction.m_imsIncluded.setZero();

    mx::improc::eigenCube<float> references( 3, 1, 3 );
    references.image( 0 ) << 1, 0, -1;
    references.image( 1 ) << 0, 1, -1;
    references.image( 2 ) << -1, 0, 1;
    mx::improc::eigenCube<float> targets( 3, 1, 1 );
    targets.image( 0 ) << 1, 0, -1;
    reductionT::imageT mask;
    std::vector<size_t> indices{ 0, 1, 2 };

    reduction.worker( references, targets, mask, indices, 0, 0 );

    REQUIRE( reduction.m_imsIncluded.rows() == 1 );
    REQUIRE( reduction.m_imsIncluded.cols() == 3 );
    REQUIRE( reduction.m_imsIncluded( 0, 0 ) == 1 );
    REQUIRE( reduction.m_imsIncluded( 0, 1 ) == 1 );
    REQUIRE( reduction.m_imsIncluded( 0, 2 ) == 0 );
    REQUIRE( reduction.m_psfsub[0].image( 0 ).isFinite().all() );
}

/// Verify KLIPreduction::worker validates inclusion metadata before entering its parallel region.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP worker inclusion metadata validation", "[KLIPreduction][worker][include][validation]" )
{
    reductionHarness reduction;
    reduction.m_meanSubMethod = mx::improc::HCI::meanSub::none;
    reduction.m_pixelTSNormMethod = mx::improc::HCI::pixelTSNorm::none;
    reduction.m_includeMethod = mx::improc::HCI::include::time;
    reduction.m_includeRefNum = 1;

    mx::improc::eigenCube<float> references( 2, 1, 2 );
    references.cube().setOnes();
    mx::improc::eigenCube<float> targets( 2, 1, 1 );
    targets.cube().setOnes();
    reductionT::imageT mask;
    std::vector<size_t> indices{ 0, 1 };

    REQUIRE_THROWS( reduction.worker( references, targets, mask, indices, 0, 0 ) );

    reduction.m_RDIimageMJD = { 1, std::numeric_limits<double>::quiet_NaN() };
    reduction.m_imageMJD = { 2 };
    REQUIRE_THROWS( reduction.worker( references, targets, mask, indices, 0, 0 ) );
}

/// Verify KLIPreduction::regions extracts a masked annulus, runs each requested mode count, and preserves geometry.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP region orchestration", "[KLIPreduction][regions][ADI][mask]" )
{
    OpenMPThreadGuard threads( 1 );
    reductionHarness reduction;
    prepareRegionReduction( reduction );
    reduction.m_maskFile = "in-memory-mask";
    reduction.m_mask.resize( 3, 3 );
    reduction.m_mask.setOnes();
    reduction.m_mask( 1, 0 ) = 0;

    REQUIRE( reduction.regions( std::vector<float>{ 0, 0 },
                                std::vector<float>{ 2, 2 },
                                std::vector<float>{ 0, 180 },
                                std::vector<float>{ 180, 360 } ) == 0 );

    REQUIRE( reduction.m_minr == std::vector<float>{ 0, 0 } );
    REQUIRE( reduction.m_maxr == std::vector<float>{ 2, 2 } );
    REQUIRE( reduction.m_minq == std::vector<float>{ 0, 180 } );
    REQUIRE( reduction.m_maxq == std::vector<float>{ 180, 360 } );
    REQUIRE( reduction.m_maxNmodes == 2 );
    REQUIRE( reduction.m_psfsub.size() == 2 );
    REQUIRE( reduction.m_imsIncluded.rows() == 3 );
    REQUIRE( reduction.m_imsIncluded.cols() == 3 );
    REQUIRE( reduction.m_imsIncluded.isOnes() );

    for( auto &cube : reduction.m_psfsub )
    {
        REQUIRE( cube.rows() == 3 );
        REQUIRE( cube.cols() == 3 );
        REQUIRE( cube.planes() == 3 );
        REQUIRE( cube.cube().isFinite().all() );
        REQUIRE( cube.image( 0 )( 1, 0 ) == 0 );
    }
    REQUIRE_FALSE( reduction.m_psfsub[0].cube().isZero() );
}

/// Verify KLIPreduction::regions uses an independent RDI library without permanently changing exclusion settings.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP RDI region orchestration", "[KLIPreduction][regions][RDI]" )
{
    OpenMPThreadGuard threads( 1 );
    reductionHarness reduction;
    prepareRegionReduction( reduction );
    reduction.m_Nmodes = { 1 };
    reduction.m_refIms.resize( 3, 3, 4 );
    reduction.m_refIms.image( 0 ) << 2, 1, 0, -1, 1, 2, 0, -2, 1;
    reduction.m_refIms.image( 1 ) << 0, -1, 2, 2, 1, 0, -2, 1, 1;
    reduction.m_refIms.image( 2 ) << 1, 2, -1, 0, -2, 1, 2, 0, 1;
    reduction.m_refIms.image( 3 ) << -1, 0, 1, 1, 2, -2, 0, 1, 2;
    reduction.m_excludeMethod = mx::improc::HCI::exclude::angle;
    reduction.m_excludeMethodMax = mx::improc::HCI::exclude::imno;
    reduction.m_minDPx = 2;
    reduction.m_maxDPx = 5;

    REQUIRE( reduction.regions( 0, 2, 0, 360 ) == 0 );

    REQUIRE( reduction.m_excludeMethod == mx::improc::HCI::exclude::angle );
    REQUIRE( reduction.m_excludeMethodMax == mx::improc::HCI::exclude::imno );
    REQUIRE( reduction.m_imsIncluded.rows() == 3 );
    REQUIRE( reduction.m_imsIncluded.cols() == 4 );
    REQUIRE( reduction.m_imsIncluded.isOnes() );
    REQUIRE( reduction.m_psfsub.size() == 1 );
    REQUIRE( reduction.m_psfsub[0].cube().isFinite().all() );
}

/// Verify KLIPreduction::regions rejects malformed modes, geometry, and preloaded image state before unsafe access.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP region validation", "[KLIPreduction][regions][validation]" )
{
    reductionHarness reduction;

    REQUIRE_THROWS( reduction.regions( 0, 2, 0, 360 ) );
    REQUIRE( reduction.m_minr.empty() );
    REQUIRE( reduction.m_psfsub.empty() );

    reduction.m_Nmodes = { 1 };
    REQUIRE_THROWS( reduction.regions( std::vector<float>{ 0 },
                                       std::vector<float>{ 2, 3 },
                                       std::vector<float>{ 0 },
                                       std::vector<float>{ 360 } ) );
    REQUIRE( reduction.m_minr.empty() );

    reduction.m_Nmodes = { 0 };
    REQUIRE_THROWS( reduction.regions( 0, 2, 0, 360 ) );

    reduction.m_Nmodes = { 1 };
    REQUIRE_THROWS( reduction.regions( 2, 2, 0, 360 ) );
    REQUIRE_THROWS( reduction.regions( -1, 2, 0, 360 ) );
    REQUIRE_THROWS( reduction.regions( 0, std::numeric_limits<float>::quiet_NaN(), 0, 360 ) );

    prepareRegionReduction( reduction );
    reduction.m_Nrows = 2;
    REQUIRE_THROWS( reduction.regions( 0, 2, 0, 360 ) );

    prepareRegionReduction( reduction );
    reduction.m_excludeMethod = mx::improc::HCI::exclude::pixel;
    reduction.m_minDPx = 1;
    REQUIRE_THROWS( reduction.regions( 0, 2, 0, 360 ) );
    REQUIRE( reduction.m_excludeMethod == mx::improc::HCI::exclude::pixel );

    prepareRegionReduction( reduction );
    reduction.m_maskFile = "empty-mask";
    reduction.m_mask.resize( 3, 3 );
    reduction.m_mask.setZero();
    REQUIRE_THROWS( reduction.regions( 0, 2, 0, 360 ) );
}

/// Verify KLIPreduction::finalProcess applies post-median subtraction, derotation, and final combination in order.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP final processing", "[KLIPreduction][finalProcess][combine]" )
{
    OpenMPThreadGuard threads( 1 );
    reductionHarness reduction;
    reduction.m_psfsub.resize( 1 );
    reduction.m_psfsub[0].resize( 2, 2, 3 );
    reduction.m_psfsub[0].image( 0 ) << 1, 2, 3, 4;
    reduction.m_psfsub[0].image( 1 ) << 3, 4, 5, 6;
    reduction.m_psfsub[0].image( 2 ) << 5, 6, 7, 8;
    reduction.m_postMedSub = true;
    reduction.m_doDerotate = true;
    reduction.m_derotF.m_angleScale = 1;
    reduction.m_derotF.m_angles = { 0, 0, 0 };
    reduction.m_combineMethod = mx::improc::HCI::combine::mean;
    reduction.m_doWriteFinim = false;
    reduction.m_doOutputPSFSub = false;

    REQUIRE( reduction.finalProcess() == 0 );

    const reductionT::imageT negative = reductionT::imageT::Constant( 2, 2, -2 );
    const reductionT::imageT positive = reductionT::imageT::Constant( 2, 2, 2 );
    REQUIRE( reduction.m_psfsub[0].image( 0 ).isApprox( negative ) );
    REQUIRE( reduction.m_psfsub[0].image( 1 ).isZero() );
    REQUIRE( reduction.m_psfsub[0].image( 2 ).isApprox( positive ) );
    REQUIRE( reduction.m_finim.rows() == 2 );
    REQUIRE( reduction.m_finim.cols() == 2 );
    REQUIRE( reduction.m_finim.planes() == 1 );
    REQUIRE( reduction.m_finim.image( 0 ).isZero() );
}

/// Verify KLIPreduction::finalProcess writes the complete KLIP configuration into a final FITS header.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP final output metadata", "[KLIPreduction][finalProcess][output][header]" )
{
    OpenMPThreadGuard threads( 1 );
    TestDirectory directory;
    reductionHarness reduction;
    reduction.m_psfsub.resize( 2 );
    for( size_t mode = 0; mode < reduction.m_psfsub.size(); ++mode )
    {
        reduction.m_psfsub[mode].resize( 2, 2, 2 );
        reduction.m_psfsub[mode].image( 0 ).setConstant( static_cast<float>( mode + 1 ) );
        reduction.m_psfsub[mode].image( 1 ).setConstant( static_cast<float>( mode + 3 ) );
    }
    reduction.m_Nims = 2;
    reduction.m_Nrows = 2;
    reduction.m_Ncols = 2;
    reduction.m_imSize = 2;
    reduction.m_Nmodes = { 1, 2 };
    reduction.m_meanSubMethod = mx::improc::HCI::meanSub::imageMedian;
    reduction.m_pixelTSNormMethod = mx::improc::HCI::pixelTSNorm::rms;
    reduction.m_rightReason = true;
    reduction.m_rightReasonRadius = 3.5;
    reduction.m_minr = { 1, 2 };
    reduction.m_maxr = { 2, 3 };
    reduction.m_minq = { 0, 180 };
    reduction.m_maxq = { 180, 360 };
    reduction.m_excludeMethod = mx::improc::HCI::exclude::angle;
    reduction.m_excludeMethodMax = mx::improc::HCI::exclude::imno;
    reduction.m_minDPx = 1.5;
    reduction.m_maxDPx = 4;
    reduction.m_includeMethod = mx::improc::HCI::include::corr;
    reduction.m_includeRefNum = 7;
    reduction.m_doDerotate = false;
    reduction.m_combineMethod = mx::improc::HCI::combine::mean;
    reduction.m_doWriteFinim = true;
    reduction.m_doOutputPSFSub = false;
    reduction.m_outputDir = directory.file( "nested" ).string();
    reduction.m_finimName = "klip-final.fits";
    reduction.m_exactFinimName = true;

    REQUIRE( reduction.finalProcess() == 0 );

    mx::improc::eigenCube<float> output;
    mx::fits::fitsHeader<mx::verbose::vv> header;
    mx::fits::fitsFile<float, mx::verbose::vv> reader;
    REQUIRE( reader.read( output, header, directory.file( "nested/klip-final.fits" ).string() ) ==
             mx::error_t::noerror );
    REQUIRE( output.rows() == 2 );
    REQUIRE( output.cols() == 2 );
    REQUIRE( output.planes() == 2 );
    REQUIRE( output.image( 0 )( 0, 0 ) == Approx( 2 ) );
    REQUIRE( output.image( 1 )( 0, 0 ) == Approx( 3 ) );
    REQUIRE( header["MEAN SUB METHOD"].String().starts_with( "imageMedian" ) );
    REQUIRE( header["PIXTS NORM METHOD"].String().starts_with( "rms" ) );
    REQUIRE( header["NMODES"].String().starts_with( "1,2" ) );
    REQUIRE( header["RIGHT REASON"].value<char>() == 1 );
    REQUIRE( header["RIGHT REASON RADIUS"].value<float>() == Approx( 3.5 ) );
    REQUIRE( header["REGMINR"].String().starts_with( "1,2" ) );
    REQUIRE( header["REGMAXR"].String().starts_with( "2,3" ) );
    REQUIRE( header["REGMINQ"].String().starts_with( "0,180" ) );
    REQUIRE( header["REGMAXQ"].String().starts_with( "180,360" ) );
    REQUIRE( header["EXMTHDMN"].String().starts_with( "angle" ) );
    REQUIRE( header["MINDPX"].value<float>() == Approx( 1.5 ) );
    REQUIRE( header["EXMTHDMX"].String().starts_with( "imno" ) );
    REQUIRE( header["MAXDPX"].value<float>() == Approx( 4 ) );
    REQUIRE( header["INMTHDMX"].String().starts_with( "corr" ) );
    REQUIRE( header["INCLREFN"].value<int>() == 7 );
}

/// Verify KLIPreduction::processPSFSub resumes the current final-processing configuration from saved KLIP products.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP saved-reduction processing", "[KLIPreduction][processPSFSub][output][combine]" )
{
    OpenMPThreadGuard threads( 1 );
    TestDirectory directory;
    reductionHarness producer;
    producer.m_psfsub.resize( 2 );
    producer.m_heads.resize( 2 );
    producer.m_heads[0].append<int>( "SOURCE", 10, "source image" );
    producer.m_heads[1].append<int>( "SOURCE", 20, "source image" );
    producer.m_Nims = 2;
    producer.m_Nrows = 2;
    producer.m_Ncols = 2;
    producer.m_imSize = 2;
    producer.m_Nmodes = { 1, 3 };
    for( size_t reduction = 0; reduction < producer.m_psfsub.size(); ++reduction )
    {
        producer.m_psfsub[reduction].resize( 2, 2, 2 );
        producer.m_psfsub[reduction].image( 0 ).setConstant( 2 * reduction + 1 );
        producer.m_psfsub[reduction].image( 1 ).setConstant( 2 * reduction + 3 );
    }
    producer.m_doDerotate = false;
    producer.m_combineMethod = mx::improc::HCI::combine::none;
    producer.m_doWriteFinim = false;
    producer.m_doOutputPSFSub = true;
    producer.m_outputDir = directory.path().string();
    producer.m_PSFSubPrefix = "saved";
    producer.m_comboWeights = { 0.25F, 0.75F };
    REQUIRE( producer.finalProcess() == 0 );

    reductionHarness resumed;
    resumed.m_doDerotate = false;
    resumed.m_combineMethod = mx::improc::HCI::combine::mean;
    resumed.m_doWriteFinim = false;
    resumed.m_doOutputPSFSub = false;
    resumed.m_weightFile = directory.file( "savedweights.dat" ).string();
    REQUIRE( resumed.processPSFSub( directory.path().string(), "saved" ) == 0 );
    REQUIRE( resumed.m_Nmodes == std::vector<int>{ 1, 3 } );
    REQUIRE( resumed.m_comboWeights == std::vector<float>{ 0.25F, 0.75F } );
    REQUIRE( resumed.m_psfsub.size() == 2 );
    REQUIRE( resumed.m_finim.planes() == 2 );
    REQUIRE( resumed.m_finim.image( 0 ).isConstant( 2.5F ) );
    REQUIRE( resumed.m_finim.image( 1 ).isConstant( 4.5F ) );

    producer.m_Nmodes.clear();
    producer.m_outputDir = directory.file( "missing-metadata" ).string();
    REQUIRE( producer.finalProcess() == 0 );
    reductionHarness missingMetadata;
    missingMetadata.m_doDerotate = false;
    missingMetadata.m_combineMethod = mx::improc::HCI::combine::mean;
    missingMetadata.m_doWriteFinim = false;
    REQUIRE_THROWS( missingMetadata.processPSFSub( producer.m_outputDir, "saved" ) );

    producer.m_Nmodes = { 1 };
    producer.m_outputDir = directory.file( "mismatched-metadata" ).string();
    REQUIRE( producer.finalProcess() == 0 );
    reductionHarness mismatchedMetadata;
    mismatchedMetadata.m_doDerotate = false;
    mismatchedMetadata.m_combineMethod = mx::improc::HCI::combine::mean;
    mismatchedMetadata.m_doWriteFinim = false;
    REQUIRE_THROWS( mismatchedMetadata.processPSFSub( producer.m_outputDir, "saved" ) );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::KLIPreduction<float,
                               mx::improc::ADIDerotator<float, mx::verbose::vv>,
                               double,
                               mx::verbose::vv>::processPSFSub( "", "" );
#endif
    // clang-format on
}

/// Verify KLIPreduction::meanSubtract rejects unimplemented and invalid methods before mutating either cube.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP mean subtraction method validation", "[KLIPreduction][meanSubtract][validation]" )
{
    reductionT reduction;
    mx::improc::eigenCube<float> references( 1, 2, 1 );
    references.image( 0 ) << 1, 2;
    mx::improc::eigenCube<float> targets = references;
    const mx::improc::eigenCube<float> original = references;
    reductionT::imageT mask;
    std::vector<float> norms;

    reduction.m_meanSubMethod = mx::improc::HCI::meanSub::imageMode;
    REQUIRE_THROWS( reduction.meanSubtract( references, targets, mask, norms ) );
    REQUIRE( references.image( 0 ).isApprox( original.image( 0 ) ) );
    REQUIRE( targets.image( 0 ).isApprox( original.image( 0 ) ) );
    REQUIRE( norms.empty() );

    reduction.m_meanSubMethod = static_cast<mx::improc::HCI::meanSub>( 99 );
    REQUIRE_THROWS( reduction.meanSubtract( references, targets, mask, norms ) );
    REQUIRE( references.image( 0 ).isApprox( original.image( 0 ) ) );
    REQUIRE( targets.image( 0 ).isApprox( original.image( 0 ) ) );
    REQUIRE( norms.empty() );
}

/// Verify KLIPreduction::meanSubtract rejects invalid normalization, cube, and mask state before mutation.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP mean subtraction input validation", "[KLIPreduction][meanSubtract][validation]" )
{
    reductionT reduction;
    reduction.m_meanSubMethod = mx::improc::HCI::meanSub::none;
    mx::improc::eigenCube<float> references( 1, 2, 1 );
    references.image( 0 ) << 1, 2;
    mx::improc::eigenCube<float> targets = references;
    const mx::improc::eigenCube<float> original = references;
    reductionT::imageT mask;
    std::vector<float> norms{ 17 };

    reduction.m_pixelTSNormMethod = mx::improc::HCI::pixelTSNorm::rmsSigmaClipped;
    REQUIRE_THROWS( reduction.meanSubtract( references, targets, mask, norms ) );
    REQUIRE( references.image( 0 ).isApprox( original.image( 0 ) ) );
    REQUIRE( targets.image( 0 ).isApprox( original.image( 0 ) ) );
    REQUIRE( norms == std::vector<float>{ 17 } );

    reduction.m_pixelTSNormMethod = static_cast<mx::improc::HCI::pixelTSNorm>( 99 );
    REQUIRE_THROWS( reduction.meanSubtract( references, targets, mask, norms ) );
    REQUIRE( references.image( 0 ).isApprox( original.image( 0 ) ) );
    REQUIRE( norms == std::vector<float>{ 17 } );

    reduction.m_pixelTSNormMethod = mx::improc::HCI::pixelTSNorm::none;
    mx::improc::eigenCube<float> wrongTargets( 2, 1, 1 );
    wrongTargets.image( 0 ) << 1, 2;
    REQUIRE_THROWS( reduction.meanSubtract( references, wrongTargets, mask, norms ) );
    REQUIRE( references.image( 0 ).isApprox( original.image( 0 ) ) );
    REQUIRE( norms == std::vector<float>{ 17 } );

    mask.resize( 1, 2 );
    mask.setZero();
    REQUIRE_THROWS( reduction.meanSubtract( references, targets, mask, norms ) );
    REQUIRE( references.image( 0 ).isApprox( original.image( 0 ) ) );
    REQUIRE( targets.image( 0 ).isApprox( original.image( 0 ) ) );
    REQUIRE( norms == std::vector<float>{ 17 } );
}

} // namespace KLIPreduction_test
} // namespace unitTest
