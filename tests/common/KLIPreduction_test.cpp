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
    using reductionT::m_Ncols;
    using reductionT::m_Nims;
    using reductionT::m_Nrows;
    using reductionT::m_psfsub;
};
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

    defaults.loadConfig( config );
    REQUIRE_FALSE( defaults.m_writeDiagnostics );
    REQUIRE( defaults.m_diagnosticDirectory == "." );

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
    REQUIRE( reduction.m_Nmodes == std::vector<int>{ 1, 3, 5 } );
    REQUIRE( reduction.m_rightReason );
    REQUIRE( reduction.m_rightReasonRadius == Approx( 3.5 ) );
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
                                      1,
                                      1,
                                      0,
                                      4,
                                      mx::improc::HCI::exclude::imno,
                                      mx::improc::HCI::exclude::none,
                                      0,
                                      derotator,
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
                                      1,
                                      0.15,
                                      0,
                                      4,
                                      mx::improc::HCI::exclude::angle,
                                      mx::improc::HCI::exclude::none,
                                      0,
                                      derotator,
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
    Eigen::ArrayXXf covariance = Eigen::ArrayXXf::Zero( 4, 4 );
    covariance( 1, 0 ) = 0.8;
    covariance( 1, 1 ) = 1.0;
    covariance( 2, 1 ) = 0.8;
    covariance( 3, 1 ) = 0.5;
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
                                      1,
                                      0,
                                      0,
                                      4,
                                      mx::improc::HCI::exclude::none,
                                      mx::improc::HCI::exclude::none,
                                      2,
                                      derotator,
                                      included );
    REQUIRE( included( 1, 0 ) == 1 );
    REQUIRE( included( 1, 1 ) == 1 );
    REQUIRE( included( 1, 2 ) == 0 );
    REQUIRE( included( 1, 3 ) == 0 );
    REQUIRE( selectedReferences.cols() == 2 );

    included.setZero();
    mx::improc::collapseCovar<float>( selectedCovariance,
                                      covariance,
                                      norms,
                                      selectedReferences,
                                      references,
                                      1,
                                      1,
                                      0,
                                      4,
                                      mx::improc::HCI::exclude::imno,
                                      mx::improc::HCI::exclude::none,
                                      2,
                                      derotator,
                                      included );
    REQUIRE( included( 1, 0 ) == 0 );
    REQUIRE( included( 1, 1 ) == 0 );
    REQUIRE( included( 1, 2 ) == 0 );
    REQUIRE( included( 1, 3 ) == 1 );
    REQUIRE( selectedReferences.cols() == 1 );
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
                                                         0,
                                                         0,
                                                         0,
                                                         2,
                                                         mx::improc::HCI::exclude::none,
                                                         mx::improc::HCI::exclude::none,
                                                         0,
                                                         derotator,
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
