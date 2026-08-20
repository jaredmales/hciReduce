/** \file P4Reduction_test.cpp
 * \brief Tests observation orchestration for Pixel Prediction Post-Processing.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "HCIobservation_test_fixture.hpp"
#include "src/common/ADIDerotator.hpp"
#include "src/common/P4Reduction.hpp"

#include <cmath>
#include <filesystem>
#include <limits>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

namespace mx
{
namespace improc
{

/** \cond P4Reduction_test_harness */
class P4ReductionTestAccess
{
  public:
    /// Invoke the production sampled-predictor finite boundary.
    static double checkedPredictorPromotion( float value /**< [in] predictor value */ )
    {
        return P4Reductionf::checkedPredictorPromotion( value );
    }

    /// Invoke the production all-double to float residual boundary.
    static float checkedResidualCast( double value /**< [in] residual value */ )
    {
        return P4Reductionf::checkedResidualCast( value );
    }

    /// Invoke the production P4PCA integer-boundary check.
    static int checkedMaximumDegreesOfFreedom( int imageCount,             /**< [in] sample count */
                                               std::size_t predictorCount, /**< [in] predictor count */
                                               bool temporallyCentered = false /**< [in] whether centering is used */ )
    {
        return P4Reductionf::checkedMaximumDegreesOfFreedom( imageCount, predictorCount, temporallyCentered );
    }

    /// Invoke the production unique-pixel ownership assignment.
    static void claimOwnership( Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic> &ownership,
                                /**< [in,out] ownership image */
                                const P4PixelCoordinate &coordinate, /**< [in] search coordinate */
                                int region /**< [in] annulus index */ )
    {
        P4Reductionf::claimOwnership( ownership, coordinate, region );
    }

    /// Invoke checked atomic publication of one enabled float diagnostic.
    static void writeDiagnostic( const P4Reductionf &reduction, /**< [in] configured reduction */
                                 const std::string &fileName,   /**< [in] diagnostic basename */
                                 const P4Reductionf::imageT &image /**< [in] diagnostic image */ )
    {
        reduction.writeDiagnostic( fileName, image );
    }
};
/** \endcond */

} // namespace improc
} // namespace mx

namespace unitTest
{
namespace P4Reduction_test
{

/** \cond P4Reduction_test_harness */
using reductionT = mx::improc::P4Reductionf;
using policyT = mx::improc::P4ExclusionPolicy;

struct reductionHarness : public reductionT
{
    using reductionT::m_combineMethod;
    using reductionT::m_dateKeyword;
    using reductionT::m_diagnosticDirectory;
    using reductionT::m_directory;
    using reductionT::m_doDerotate;
    using reductionT::m_doOutputPSFSub;
    using reductionT::m_doWriteFinim;
    using reductionT::m_fileList;
    using reductionT::m_filesRead;
    using reductionT::m_finim;
    using reductionT::m_heads;
    using reductionT::m_imSize;
    using reductionT::m_mask;
    using reductionT::m_maskFile;
    using reductionT::m_Ncols;
    using reductionT::m_Nims;
    using reductionT::m_Npix;
    using reductionT::m_Nrows;
    using reductionT::m_outputDir;
    using reductionT::m_prefix;
    using reductionT::m_preProcess_only;
    using reductionT::m_psfsub;
    using reductionT::m_PSFSubPrefix;
    using reductionT::m_psfsubValidity;
    using reductionT::m_RDIdirectory;
    using reductionT::m_RDIfileList;
    using reductionT::m_RDIfileListFile;
    using reductionT::m_RDIfilesRead;
    using reductionT::m_RDIprefix;
    using reductionT::m_refIms;
    using reductionT::m_rejectNonFiniteTargetInput;
    using reductionT::m_skipPreProcess;
    using reductionT::m_tgtIms;

    bool m_postReadCalled{ false }; ///< Whether inherited input loading reached the post-read hook.

    /// Suppress dataset-specific post-read behavior in in-memory tests.
    void postReadFiles() override
    {
        m_postReadCalled = true;
    }
};

/// Capture `std::cerr` for one serialized progress-reporting test scope.
class CerrCapture
{
  private:
    std::ostringstream m_buffer; ///< Captured output storage constructed before stream redirection.

    std::streambuf *m_original;  ///< Original `std::cerr` buffer restored at scope exit.

  public:
    /// Redirect `std::cerr` into the owned buffer.
    CerrCapture() : m_buffer(), m_original( std::cerr.rdbuf( m_buffer.rdbuf() ) )
    {
    }

    /// Restore the original `std::cerr` buffer.
    ~CerrCapture()
    {
        std::cerr.rdbuf( m_original );
    }

    /// Return the output captured so far.
    std::string str() const
    {
        return m_buffer.str();
    }
};

/// Configure compact geometry and a preloaded cube suitable for P4 tests.
void prepareReduction( reductionHarness &reduction, /**< [out] configured reduction */
                       int imageCount = 3,          /**< [in] number of temporal samples */
                       int rows = 31,               /**< [in] image row count */
                       int columns = 31 /**< [in] image column count */ )
{
    reduction.m_filesRead = true;
    reduction.m_RDIfilesRead = true;
    reduction.m_Nims = imageCount;
    reduction.m_Nrows = rows;
    reduction.m_Ncols = columns;
    reduction.m_Npix = rows * columns;
    reduction.m_tgtIms.resize( rows, columns, imageCount );
    for( int image = 0; image < imageCount; ++image )
    {
        reduction.m_tgtIms.image( image ).setConstant( static_cast<float>( image + 1 ) );
    }
    reduction.m_minRadius = { 5 };
    reduction.m_maxRadius = { 6 };
    reduction.m_modeFractions = { 0.5f };
    reduction.m_orDeltaRadiusInner = 2;
    reduction.m_orDeltaRadiusOuter = 2;
    reduction.m_orArcHalfWidth = 4;
    reduction.m_orMaxHalfAngle = 60;
    reduction.m_psfRadius = 0.5f;
    reduction.m_exclusionPolicy = policyT::kernelSupport;
    reduction.m_exclusionRadiusBuffer = 0;
    reduction.m_rankTolerance = 1e-6;
    reduction.m_derotF.m_angleScale = 1;
    reduction.m_derotF.m_angles.assign( static_cast<std::size_t>( imageCount ), 0 );
    reduction.m_writeDiagnostics = false;
    reduction.m_doDerotate = false;
    reduction.m_combineMethod = mx::improc::HCI::combine::none;
    reduction.m_doWriteFinim = false;
    reduction.m_doOutputPSFSub = false;
}

/// Load one P4 configuration file through the production setup/load API.
void readReductionConfig( reductionT &reduction,             /**< [out] configured reduction */
                          const std::filesystem::path &path, /**< [in] temporary configuration path */
                          const std::string &configurationText /**< [in] configuration contents */ )
{
    writeTextFile( path, configurationText );
    mx::app::appConfigurator config;
    reduction.setupConfig( config );
    if( config.readConfig( path.string() ) != 0 )
    {
        throw std::runtime_error( "could not read P4 configuration" );
    }
    reduction.loadConfig( config );
}

/// Compare the science, validity, ownership, modes, and count outcomes of two completed reductions.
void requireSameReduction( const reductionHarness &left, /**< [in] first completed reduction */
                           const reductionHarness &right /**< [in] second completed reduction */ )
{
    REQUIRE( left.m_realizedModes == right.m_realizedModes );
    REQUIRE( ( left.m_ownership == right.m_ownership ).all() );
    REQUIRE( left.m_psfsub.size() == right.m_psfsub.size() );
    REQUIRE( left.m_psfsubValidity.size() == right.m_psfsubValidity.size() );
    REQUIRE( left.m_regionStatistics.size() == right.m_regionStatistics.size() );
    for( std::size_t region = 0; region < left.m_regionStatistics.size(); ++region )
    {
        const auto &leftStatistics = left.m_regionStatistics[region];
        const auto &rightStatistics = right.m_regionStatistics[region];
        REQUIRE( leftStatistics.searchPixelCount == rightStatistics.searchPixelCount );
        REQUIRE( leftStatistics.predictorCount == rightStatistics.predictorCount );
        REQUIRE( leftStatistics.maximumDegreesOfFreedom == rightStatistics.maximumDegreesOfFreedom );
        REQUIRE( leftStatistics.minimumNumericalRank == rightStatistics.minimumNumericalRank );
        REQUIRE( leftStatistics.validLocalFitCount == rightStatistics.validLocalFitCount );
        REQUIRE( leftStatistics.maskedLocalFitCount == rightStatistics.maskedLocalFitCount );
        REQUIRE( leftStatistics.supportInvalidLocalFitCount == rightStatistics.supportInvalidLocalFitCount );
        REQUIRE( leftStatistics.rankInvalidCounts == rightStatistics.rankInvalidCounts );
    }
    for( std::size_t output = 0; output < left.m_psfsub.size(); ++output )
    {
        for( int image = 0; image < left.m_psfsub[output].planes(); ++image )
        {
            REQUIRE( ( left.m_psfsubValidity[output].image( image ) == right.m_psfsubValidity[output].image( image ) )
                         .all() );
            for( int column = 0; column < left.m_psfsub[output].cols(); ++column )
            {
                for( int row = 0; row < left.m_psfsub[output].rows(); ++row )
                {
                    if( left.m_psfsubValidity[output].image( image )( row, column ) == 1 )
                    {
                        REQUIRE( left.m_psfsub[output].image( image )( row, column ) ==
                                 Approx( right.m_psfsub[output].image( image )( row, column ) ).margin( 2e-5 ) );
                    }
                    else
                    {
                        REQUIRE( mx::math::isNan( left.m_psfsub[output].image( image )( row, column ) ) );
                        REQUIRE( mx::math::isNan( right.m_psfsub[output].image( image )( row, column ) ) );
                    }
                }
            }
        }
    }
}

/// Controlled P4PCA solver failure used to verify serial rethrow after an OpenMP loop.
MXLAPACK_INT failingEigenSolver( mx::improc::P4PCA::matrixT &eigenvectors, /**< [out] unused eigenvectors */
                                 mx::improc::P4PCA::matrixT &eigenvalues,  /**< [out] unused eigenvalues */
                                 mx::improc::P4PCA::matrixT &covariance,   /**< [in] unused covariance */
                                 int modeCount,                            /**< [in] unused mode count */
                                 mx::improc::P4PCA::workspaceT &workspace /**< [in,out] unused workspace */ )
{
    static_cast<void>( eigenvectors );
    static_cast<void>( eigenvalues );
    static_cast<void>( covariance );
    static_cast<void>( modeCount );
    static_cast<void>( workspace );
    return 73;
}

/// Restore the production P4PCA eigensolver after one controlled failure scope.
struct eigenSolverReset
{
    /// Install the failing eigensolver.
    eigenSolverReset()
    {
        mx::improc::detail::p4PCASetEigenSolverForTesting( &failingEigenSolver );
    }

    /// Restore the production eigensolver.
    ~eigenSolverReset()
    {
        mx::improc::detail::p4PCAResetEigenSolverForTesting();
    }
};

/// Sample one fixed sky coordinate directly from an unrotated detector image.
float directRotatedSample( const reductionT::imageT &image, /**< [in] preprocessed detector-frame image */
                           double skyRow,                   /**< [in] fixed sky-frame row coordinate */
                           double skyColumn,                /**< [in] fixed sky-frame column coordinate */
                           double angle /**< [in] counter-clockwise derotation angle in radians */ )
{
    const double xCenter = 0.5 * static_cast<double>( image.rows() - 1 );
    const double yCenter = 0.5 * static_cast<double>( image.cols() - 1 );
    const double cosine = std::cos( angle );
    const double sine = std::sin( angle );
    const double deltaRow = skyRow - xCenter;
    const double deltaColumn = skyColumn - yCenter;
    const double rawRow = deltaRow * cosine + deltaColumn * sine + xCenter;
    const double rawColumn = -deltaRow * sine + deltaColumn * cosine + yCenter;
    const int footprintRow =
        static_cast<int>( std::floor( rawRow ) ) - static_cast<int>( mx::improc::P4PixelGridf::leftBuffer );
    const int footprintColumn =
        static_cast<int>( std::floor( rawColumn ) ) - static_cast<int>( mx::improc::P4PixelGridf::leftBuffer );
    if( footprintRow < 0 || footprintColumn < 0 ||
        footprintRow + static_cast<int>( mx::improc::P4PixelGridf::width ) > image.rows() ||
        footprintColumn + static_cast<int>( mx::improc::P4PixelGridf::width ) > image.cols() )
    {
        throw std::out_of_range( "independent rotated sample crosses the detector edge" );
    }

    mx::improc::P4PixelGridf::kernelT kernel;
    mx::improc::cubicConvolTransform<float> transform;
    transform( kernel,
               static_cast<float>( rawRow - std::floor( rawRow ) ),
               static_cast<float>( rawColumn - std::floor( rawColumn ) ) );
    float result{ 0 };
    for( int rowOffset = 0; rowOffset < static_cast<int>( mx::improc::P4PixelGridf::width ); ++rowOffset )
    {
        for( int columnOffset = 0; columnOffset < static_cast<int>( mx::improc::P4PixelGridf::width ); ++columnOffset )
        {
            result +=
                image( footprintRow + rowOffset, footprintColumn + columnOffset ) * kernel( rowOffset, columnOffset );
        }
    }
    return result;
}
/** \endcond */

/// Verify P4Reduction setupConfig/loadConfig expose required invalid defaults and load every supported P4 option.
/** \ingroup P4Reduction_unit_tests */
TEST_CASE( "P4 reduction configuration", "[P4Reduction][config]" )
{
    reductionHarness defaults;
    REQUIRE( defaults.m_rejectNonFiniteTargetInput );
    REQUIRE( defaults.m_minRadius.empty() );
    REQUIRE( defaults.m_maxRadius.empty() );
    REQUIRE( defaults.m_modeFractions.empty() );
    REQUIRE( defaults.m_regressionFrame == mx::improc::P4RegressionFrame::detector );
    REQUIRE( defaults.m_numberImages == 0 );
    REQUIRE_FALSE( defaults.m_exclusionPolicy.has_value() );
    REQUIRE( mx::math::isNan( defaults.m_orDeltaRadiusInner ) );
    REQUIRE( mx::math::isNan( defaults.m_exclusionRadiusBuffer ) );
    REQUIRE( defaults.m_rankTolerance == Approx( 1e-12 ) );
    REQUIRE_FALSE( defaults.m_writeDiagnostics );
    REQUIRE( defaults.m_diagnosticDirectory == "." );

    mx::app::appConfigurator registered;
    defaults.setupConfig( registered );
    REQUIRE( registered.m_targets.at( "p4.modeFractions" ).clType == mx::app::argType::Required );
    REQUIRE( registered.m_targets.at( "p4.regressionFrame" ).clType == mx::app::argType::Required );
    REQUIRE( registered.m_targets.at( "p4.numberImages" ).helpType == "int" );
    REQUIRE( registered.m_targets.at( "p4.writeDiagnostics" ).clType == mx::app::argType::True );
    REQUIRE( registered.m_targets.at( "p4.orMaxHalfAngle" ).helpType == "float" );
    REQUIRE_NOTHROW( defaults.loadConfig( registered ) );
    REQUIRE_FALSE( defaults.m_exclusionPolicy.has_value() );

    reductionHarness existingPolicy;
    existingPolicy.m_exclusionPolicy = policyT::sampleCenter;
    mx::app::appConfigurator existingPolicyConfig;
    existingPolicy.setupConfig( existingPolicyConfig );
    REQUIRE_NOTHROW( existingPolicy.loadConfig( existingPolicyConfig ) );
    REQUIRE( existingPolicy.m_exclusionPolicy == policyT::sampleCenter );

    TestDirectory directory;
    reductionHarness configured;
    readReductionConfig( configured,
                         directory.file( "p4.conf" ),
                         "[geom]\nminRadius=5,8\nmaxRadius=6,9\n"
                         "[p4]\nmodeFractions=0.25,0.5\nregressionFrame=rotated\n"
                         "orDeltaRadiusInner=2\norDeltaRadiusOuter=3\n"
                         "orArcHalfWidth=4\norMaxHalfAngle=90\npsfRadius=1.5\n"
                         "exclusionPolicy=sampleCenter\nexclusionRadiusBuffer=0.5\nrankTolerance=1e-8\n"
                         "writeDiagnostics=true\ndiagnosticDirectory=" +
                             directory.file( "diagnostics" ).string() + "\n" );
    REQUIRE( configured.m_minRadius == std::vector<float>{ 5, 8 } );
    REQUIRE( configured.m_maxRadius == std::vector<float>{ 6, 9 } );
    REQUIRE( configured.m_modeFractions == std::vector<float>{ 0.25f, 0.5f } );
    REQUIRE( configured.m_regressionFrame == mx::improc::P4RegressionFrame::rotated );
    REQUIRE( configured.m_numberImages == 0 );
    REQUIRE( configured.m_orDeltaRadiusInner == 2 );
    REQUIRE( configured.m_orDeltaRadiusOuter == 3 );
    REQUIRE( configured.m_orArcHalfWidth == 4 );
    REQUIRE( configured.m_orMaxHalfAngle == 90 );
    REQUIRE( configured.m_psfRadius == 1.5f );
    REQUIRE( configured.m_exclusionPolicy == policyT::sampleCenter );
    REQUIRE( configured.m_exclusionRadiusBuffer == 0.5f );
    REQUIRE( configured.m_rankTolerance == Approx( 1e-8 ) );
    REQUIRE( configured.m_writeDiagnostics );
    reductionT::fitsHeaderT configuredHeader;
    configured.appendReductionHeader( configuredHeader );
    REQUIRE( configuredHeader["P4 FRAME"].String().starts_with( "rotated" ) );
    REQUIRE( configuredHeader["P4 EXCLUSION POLICY"].String().starts_with( "sampleCenter" ) );
    REQUIRE( configuredHeader["P4 MIN RADIUS"].String().find( ',' ) != std::string::npos );

    reductionHarness kernelPolicy;
    readReductionConfig( kernelPolicy,
                         directory.file( "kernel-policy.conf" ),
                         "[p4]\nexclusionPolicy=kernelSupport\n" );
    REQUIRE( kernelPolicy.m_exclusionPolicy == policyT::kernelSupport );

    reductionHarness invalidPolicy;
    REQUIRE_THROWS( readReductionConfig( invalidPolicy,
                                         directory.file( "invalid-policy.conf" ),
                                         "[p4]\nexclusionPolicy=automatic\n" ) );

    reductionHarness invalidFrame;
    REQUIRE_THROWS( readReductionConfig( invalidFrame,
                                         directory.file( "invalid-frame.conf" ),
                                         "[p4]\nregressionFrame=materialized\n" ) );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::P4Reduction<float, mx::improc::ADIDerotator<float, mx::verbose::vv>, mx::verbose::vv>
        doxygenReduction;
    mx::app::appConfigurator doxygenConfig;
    doxygenReduction.setupConfig( doxygenConfig );
    doxygenReduction.loadConfig( doxygenConfig );
#endif
    // clang-format on
}

/// Verify P4Reduction enforces both finite interpolation and all-double-to-float storage boundaries.
/** \ingroup P4Reduction_unit_tests */
TEST_CASE( "P4 reduction arithmetic boundaries", "[P4Reduction][finite][conversion]" )
{
    REQUIRE( mx::improc::P4ReductionTestAccess::checkedPredictorPromotion( 2.5f ) == Approx( 2.5 ) );
    REQUIRE_THROWS(
        mx::improc::P4ReductionTestAccess::checkedPredictorPromotion( std::numeric_limits<float>::infinity() ) );
    REQUIRE( mx::improc::P4ReductionTestAccess::checkedResidualCast( 3.25 ) == Approx( 3.25f ) );
    REQUIRE_THROWS(
        mx::improc::P4ReductionTestAccess::checkedResidualCast( std::numeric_limits<double>::quiet_NaN() ) );
    REQUIRE_THROWS( mx::improc::P4ReductionTestAccess::checkedResidualCast(
        static_cast<double>( std::numeric_limits<float>::max() ) * 2.0 ) );
    REQUIRE( mx::improc::P4ReductionTestAccess::checkedMaximumDegreesOfFreedom( 3, 8 ) == 3 );
    REQUIRE( mx::improc::P4ReductionTestAccess::checkedMaximumDegreesOfFreedom( 8, 3 ) == 3 );
    REQUIRE( mx::improc::P4ReductionTestAccess::checkedMaximumDegreesOfFreedom( 3, 8, true ) == 2 );
    REQUIRE( mx::improc::P4ReductionTestAccess::checkedMaximumDegreesOfFreedom( 8, 3, true ) == 3 );
    REQUIRE_THROWS( mx::improc::P4ReductionTestAccess::checkedMaximumDegreesOfFreedom( 1, 3, true ) );
    REQUIRE_THROWS( mx::improc::P4ReductionTestAccess::checkedMaximumDegreesOfFreedom(
        3,
        static_cast<std::size_t>( std::numeric_limits<int>::max() ) + 1 ) );

    Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic> ownership =
        Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic>::Constant( 2, 2, -1 );
    const mx::improc::P4PixelCoordinate coordinate( 1, 1 );
    REQUIRE_NOTHROW( mx::improc::P4ReductionTestAccess::claimOwnership( ownership, coordinate, 3 ) );
    REQUIRE( ownership( 1, 1 ) == 3 );
    REQUIRE_THROWS( mx::improc::P4ReductionTestAccess::claimOwnership( ownership, coordinate, 4 ) );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::P4Reduction<float,
                            mx::improc::ADIDerotator<float, mx::verbose::vv>,
                            mx::verbose::vv>::checkedPredictorPromotion( 1 );
    mx::improc::P4Reduction<float,
                            mx::improc::ADIDerotator<float, mx::verbose::vv>,
                            mx::verbose::vv>::checkedResidualCast( 1 );
    mx::improc::P4Reduction<float,
                            mx::improc::ADIDerotator<float, mx::verbose::vv>,
                            mx::verbose::vv>::checkedMaximumDegreesOfFreedom( 1, 1 );
    mx::improc::P4Reduction<float,
                            mx::improc::ADIDerotator<float, mx::verbose::vv>,
                            mx::verbose::vv>::claimOwnership( ownership, coordinate, 0 );
#endif
    // clang-format on
}

/// Verify P4Reduction predicts an exactly rank-one synthetic cube and populates authoritative validity cubes.
/** \ingroup P4Reduction_unit_tests */
TEST_CASE( "P4 reduction exact synthetic prediction", "[P4Reduction][reduce][prediction][validity]" )
{
    OpenMPThreadGuard threads( 1 );
    reductionHarness reduction;
    prepareReduction( reduction );

    REQUIRE( reduction.reduce() == 0 );
    REQUIRE( reduction.m_realizedModes == std::vector<std::vector<int>>{ { 1 } } );
    REQUIRE( reduction.m_regionStatistics.size() == 1 );
    REQUIRE( reduction.m_regionStatistics[0].predictorCount > 0 );
    REQUIRE( reduction.m_regionStatistics[0].minimumNumericalRank == 1 );
    REQUIRE( reduction.m_psfsub.size() == 1 );
    REQUIRE( reduction.m_psfsubValidity.size() == 1 );

    std::size_t validSamples{ 0 };
    for( int image = 0; image < reduction.m_Nims; ++image )
    {
        for( int column = 0; column < reduction.m_Ncols; ++column )
        {
            for( int row = 0; row < reduction.m_Nrows; ++row )
            {
                if( reduction.m_psfsubValidity[0].image( image )( row, column ) == 1 )
                {
                    ++validSamples;
                    REQUIRE( reduction.m_psfsub[0].image( image )( row, column ) == Approx( 0 ).margin( 2e-5 ) );
                }
            }
        }
    }
    REQUIRE( validSamples == reduction.m_regionStatistics[0].searchPixelCount * 3 );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::P4Reduction<float, mx::improc::ADIDerotator<float, mx::verbose::vv>, mx::verbose::vv>
        doxygenReduction;
    doxygenReduction.reduce();
    doxygenReduction.regions( { 5 }, { 6 } );
#endif
    // clang-format on
}

/// Verify detector-frame P4 uses nearest qualifying temporal images with one-sided endpoint replacement.
/** This exercises P4Reduction::reduce() with multi-image detector predictors. */
/** \ingroup P4Reduction_unit_tests */
TEST_CASE( "P4 detector multi-image temporal selection", "[P4Reduction][reduce][detector][temporal]" )
{
    SECTION( "one neighbor per direction supports decreasing rotation" )
    {
        OpenMPThreadGuard threads( 1 );
        reductionHarness reduction;
        prepareReduction( reduction, 5 );
        reduction.m_numberImages = 1;
        reduction.m_derotF.m_angles = { 80, 60, 40, 20, 0 };

        REQUIRE( reduction.reduce() == 0 );
        REQUIRE( reduction.m_temporalSelections ==
                 std::vector<std::vector<std::vector<int>>>{
                     { { 0, 1, 2 }, { 1, 0, 2 }, { 2, 1, 3 }, { 3, 2, 4 }, { 4, 3, 2 } } } );
        REQUIRE( reduction.m_regionStatistics[0].targetImageCount == 5 );
        REQUIRE( reduction.m_regionStatistics[0].temporalNumberImages == 1 );
        REQUIRE( reduction.m_regionStatistics[0].temporalPsfRadius == Approx( 2 * reduction.m_psfRadius ) );
        REQUIRE( reduction.m_regionStatistics[0].predictorCount == 19 );

        for( int image = 0; image < reduction.m_Nims; ++image )
        {
            REQUIRE( reduction.m_psfsubValidity[0].image( image ).maxCoeff() == 1 );
        }
    }

    SECTION( "multiple neighbors per direction worsen endpoint truncation" )
    {
        OpenMPThreadGuard threads( 1 );
        reductionHarness reduction;
        prepareReduction( reduction, 7 );
        reduction.m_numberImages = 2;
        reduction.m_derotF.m_angles = { 0, 20, 40, 60, 80, 100, 120 };

        REQUIRE( reduction.reduce() == 0 );
        REQUIRE( reduction.m_temporalSelections[0].front() == std::vector<int>{ 0, 1, 2, 3, 4 } );
        REQUIRE( reduction.m_temporalSelections[0].at( 2 ) == std::vector<int>{ 2, 1, 0, 3, 4 } );
        REQUIRE( reduction.m_temporalSelections[0].back() == std::vector<int>{ 6, 5, 4, 3, 2 } );
        REQUIRE( reduction.m_regionStatistics[0].targetImageCount == 7 );
        REQUIRE( reduction.m_regionStatistics[0].temporalNumberImages == 2 );
        REQUIRE( reduction.m_regionStatistics[0].predictorCount == 21 );
    }

    SECTION( "unattainable radius degrades to the largest structurally usable value" )
    {
        OpenMPThreadGuard threads( 1 );
        reductionHarness reduction;
        prepareReduction( reduction, 4 );
        reduction.m_numberImages = 1;
        reduction.m_derotF.m_angles = { 0, 2, 4, 6 };

        REQUIRE( reduction.reduce() == 0 );
        REQUIRE( reduction.m_regionStatistics[0].temporalNumberImages == 1 );
        REQUIRE( reduction.m_regionStatistics[0].temporalPsfRadius > 0 );
        REQUIRE( reduction.m_regionStatistics[0].temporalPsfRadius < 2 * reduction.m_psfRadius );
        REQUIRE( reduction.m_regionStatistics[0].targetImageCount >= 2 );
    }

    SECTION( "zero usable rotation reverts the annulus to same-image predictors" )
    {
        OpenMPThreadGuard threads( 1 );
        reductionHarness reduction;
        prepareReduction( reduction, 5 );
        reduction.m_numberImages = 1;

        REQUIRE( reduction.reduce() == 0 );
        REQUIRE( reduction.m_regionStatistics[0].temporalNumberImages == 0 );
        REQUIRE( reduction.m_regionStatistics[0].temporalPsfRadius == 0 );
        REQUIRE( reduction.m_regionStatistics[0].targetImageCount == 5 );
    }
}

/// Verify rotated P4 centers its fit, applies the predictor to uncentered data, and reports sky support.
/** \ingroup P4Reduction_unit_tests */
TEST_CASE( "P4 rotated reduction applies uncentered predictors",
           "[P4Reduction][rotated][centered][uncentered][validity]" )
{
    OpenMPThreadGuard threads( 1 );
    reductionHarness reduction;
    prepareReduction( reduction );
    reduction.m_regressionFrame = mx::improc::P4RegressionFrame::rotated;
    reduction.m_doDerotate = true;
    reduction.m_derotF.m_angles = { 0, 17, -11 };

    REQUIRE( reduction.reduce() == 0 );
    REQUIRE( reduction.m_realizedModes == std::vector<std::vector<int>>{ { 1 } } );
    REQUIRE( reduction.m_regionStatistics.size() == 1 );
    const auto &statistics = reduction.m_regionStatistics.front();
    REQUIRE( statistics.maximumDegreesOfFreedom == 2 );
    REQUIRE( statistics.minimumNumericalRank == 1 );
    REQUIRE( statistics.maskedLocalFitCount == 0 );
    REQUIRE( statistics.validLocalFitCount + statistics.supportInvalidLocalFitCount == statistics.searchPixelCount );

    std::size_t validPixels{ 0 };
    for( int column = 0; column < reduction.m_Ncols; ++column )
    {
        for( int row = 0; row < reduction.m_Nrows; ++row )
        {
            if( reduction.m_psfsubValidity[0].image( 0 )( row, column ) == 0 )
            {
                continue;
            }
            ++validPixels;
            for( int image = 0; image < reduction.m_Nims; ++image )
            {
                REQUIRE( reduction.m_psfsubValidity[0].image( image )( row, column ) == 1 );
                REQUIRE( reduction.m_psfsub[0].image( image )( row, column ) == Approx( 0 ).margin( 3e-5 ) );
            }
        }
    }
    REQUIRE( validPixels == statistics.validLocalFitCount );
    for( int image = 0; image < reduction.m_Nims; ++image )
    {
        REQUIRE( reduction.m_psfsubValidity[0].image( image ).sum() == Approx( statistics.validLocalFitCount ) );
    }

    reductionT::fitsHeaderT header;
    reduction.appendReductionHeader( header );
    REQUIRE( header["P4 FRAME"].String().starts_with( "rotated" ) );
    REQUIRE( header["P4 SUPPORT INVALID FIT COUNT"].String() ==
             std::to_string( statistics.supportInvalidLocalFitCount ) );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::P4Reduction<float, mx::improc::ADIDerotator<float, mx::verbose::vv>, mx::verbose::vv>
        doxygenReduction;
    doxygenReduction.reduce();
#endif
    // clang-format on
}

/// Verify the complete rotated worker against an independent direct-sampling and thin-SVD numerical oracle.
/** \ingroup P4Reduction_unit_tests */
TEST_CASE( "P4 rotated reduction matches an independent numerical oracle", "[P4Reduction][rotated][oracle][frame]" )
{
    OpenMPThreadGuard threads( 1 );
    reductionHarness reduction;
    prepareReduction( reduction, 5 );
    reduction.m_regressionFrame = mx::improc::P4RegressionFrame::rotated;
    reduction.m_exclusionPolicy = policyT::sampleCenter;
    reduction.m_modeFractions = { 0.25F };
    reduction.m_doDerotate = true;
    reduction.m_derotF.m_angles = { 3.25F, 17.5F, -11.75F, 33.125F, -27.625F };

    const double xCenter = 0.5 * static_cast<double>( reduction.m_Nrows - 1 );
    const double yCenter = 0.5 * static_cast<double>( reduction.m_Ncols - 1 );
    for( int image = 0; image < reduction.m_Nims; ++image )
    {
        const double time = static_cast<double>( image - 2 );
        const double quadraticTime = time * time - 2.0;
        const double alternatingTime = image % 2 == 0 ? -1.0 : 1.0;
        for( int column = 0; column < reduction.m_Ncols; ++column )
        {
            const double centeredColumn = static_cast<double>( column ) - yCenter;
            for( int row = 0; row < reduction.m_Nrows; ++row )
            {
                const double centeredRow = static_cast<double>( row ) - xCenter;
                const double spatial = 20.0 + 0.07 * centeredRow + 0.11 * centeredColumn +
                                       0.013 * centeredRow * centeredColumn + 0.004 * centeredRow * centeredRow -
                                       0.006 * centeredColumn * centeredColumn;
                const double linearCoefficient = 1.2 + 0.025 * centeredRow - 0.017 * centeredColumn;
                const double quadraticCoefficient = 0.22 + 0.003 * centeredRow * centeredColumn;
                const double alternatingCoefficient = 0.09 + 0.002 * centeredRow - 0.003 * centeredColumn;
                reduction.m_tgtIms.image( image )( row, column ) =
                    static_cast<float>( spatial + linearCoefficient * time + quadraticCoefficient * quadraticTime +
                                        alternatingCoefficient * alternatingTime );
            }
        }
    }

    mx::improc::P4PixelGridf candidateGrid;
    candidateGrid.resize( reduction.m_Nrows, reduction.m_Ncols );
    candidateGrid.candidateRegion(
        mx::improc::P4PixelGridRegion( 5.0, 6.0, 2.0, 2.0, 4.0, 60.0, 0.5, policyT::sampleCenter, 0.0 ) );

    constexpr int targetRow = 12;
    constexpr int targetColumn = 19;
    std::size_t targetSearch = candidateGrid.searchPixelCount();
    for( std::size_t search = 0; search < candidateGrid.searchPixelCount(); ++search )
    {
        const mx::improc::P4PixelCoordinate &coordinate = candidateGrid.searchPixel( search ).coordinate();
        if( coordinate.row() == targetRow && coordinate.column() == targetColumn )
        {
            REQUIRE( targetSearch == candidateGrid.searchPixelCount() );
            targetSearch = search;
        }
    }
    REQUIRE( targetSearch < candidateGrid.searchPixelCount() );

    std::vector<std::size_t> retainedCandidates;
    for( std::size_t candidate = 0; candidate < candidateGrid.candidatePredictorCount(); ++candidate )
    {
        const mx::improc::P4PixelCoordinate &offset = candidateGrid.candidatePredictorOffset( candidate );
        if( std::hypot( static_cast<double>( offset.row() ), static_cast<double>( offset.column() ) ) > 0.5 )
        {
            retainedCandidates.push_back( candidate );
        }
    }
    REQUIRE_FALSE( retainedCandidates.empty() );

    Eigen::VectorXd target( reduction.m_Nims );
    Eigen::MatrixXd predictors( reduction.m_Nims, static_cast<Eigen::Index>( retainedCandidates.size() ) );
    const double targetAngle =
        std::atan2( static_cast<double>( targetColumn ) - yCenter, static_cast<double>( targetRow ) - xCenter ) -
        0.5 * std::numbers::pi_v<double>;
    const double targetCosine = std::cos( targetAngle );
    const double targetSine = std::sin( targetAngle );
    for( int image = 0; image < reduction.m_Nims; ++image )
    {
        const double derotationAngle =
            static_cast<double>( reduction.m_derotF.m_angles[image] ) * std::numbers::pi_v<double> / 180.0;
        target( image ) =
            directRotatedSample( reduction.m_tgtIms.image( image ), targetRow, targetColumn, derotationAngle );
        for( std::size_t predictor = 0; predictor < retainedCandidates.size(); ++predictor )
        {
            const mx::improc::P4PixelCoordinate &offset =
                candidateGrid.candidatePredictorOffset( retainedCandidates[predictor] );
            const double predictorSkyRow = static_cast<double>( targetRow ) +
                                           static_cast<double>( offset.row() ) * targetCosine -
                                           static_cast<double>( offset.column() ) * targetSine;
            const double predictorSkyColumn = static_cast<double>( targetColumn ) +
                                              static_cast<double>( offset.row() ) * targetSine +
                                              static_cast<double>( offset.column() ) * targetCosine;
            predictors( image, static_cast<Eigen::Index>( predictor ) ) =
                directRotatedSample( reduction.m_tgtIms.image( image ),
                                     predictorSkyRow,
                                     predictorSkyColumn,
                                     derotationAngle );
        }
    }

    Eigen::VectorXd centeredTarget = target;
    centeredTarget.array() -= centeredTarget.mean();
    Eigen::MatrixXd centeredPredictors = predictors;
    for( Eigen::Index predictor = 0; predictor < centeredPredictors.cols(); ++predictor )
    {
        centeredPredictors.col( predictor ).array() -= centeredPredictors.col( predictor ).mean();
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> singularValueDecomposition( centeredPredictors,
                                                                  Eigen::ComputeThinU | Eigen::ComputeThinV );
    REQUIRE( singularValueDecomposition.singularValues()( 0 ) >
             5.0 * singularValueDecomposition.singularValues()( 1 ) );
    const Eigen::VectorXd leadingCoefficients = singularValueDecomposition.matrixV().col( 0 ) *
                                                ( singularValueDecomposition.matrixU().col( 0 ).dot( centeredTarget ) /
                                                  singularValueDecomposition.singularValues()( 0 ) );
    const Eigen::VectorXd expectedResidual = target - predictors * leadingCoefficients;

    REQUIRE( reduction.reduce() == 0 );
    REQUIRE( reduction.m_ownership( targetRow, targetColumn ) == 0 );
    REQUIRE( reduction.m_regionStatistics.front().predictorCount == retainedCandidates.size() );
    REQUIRE( reduction.m_realizedModes == std::vector<std::vector<int>>{ { 1 } } );

    bool accidentalDerotationChangesScience{ false };
    bool accidentalDerotationChangesValidity{ false };
    for( int image = 0; image < reduction.m_Nims; ++image )
    {
        REQUIRE( reduction.m_psfsubValidity[0].image( image )( targetRow, targetColumn ) == 1 );
        const float actualResidual = reduction.m_psfsub[0].image( image )( targetRow, targetColumn );
        REQUIRE( actualResidual == Approx( expectedResidual( image ) ).margin( 3e-4 ) );

        reductionT::imageT sanitized = reduction.m_psfsub[0].image( image );
        for( int column = 0; column < reduction.m_Ncols; ++column )
        {
            for( int row = 0; row < reduction.m_Nrows; ++row )
            {
                if( reduction.m_psfsubValidity[0].image( image )( row, column ) == 0 )
                {
                    sanitized( row, column ) = 0;
                }
            }
        }

        const float derotationAngle = reduction.m_derotF.derotAngle( static_cast<std::size_t>( image ) );
        reductionT::imageT doubleRotatedScience;
        reductionT::imageT doubleRotatedValidity;
        mx::improc::imageRotate( doubleRotatedScience,
                                 sanitized,
                                 derotationAngle,
                                 mx::improc::cubicConvolTransform<float>() );
        mx::improc::imageRotate( doubleRotatedValidity,
                                 reduction.m_psfsubValidity[0].image( image ),
                                 derotationAngle,
                                 mx::improc::detail::completeCubicFootprintTransform<float>() );
        const bool hypotheticalValid = doubleRotatedValidity( targetRow, targetColumn ) == 1;
        if( !hypotheticalValid )
        {
            doubleRotatedScience( targetRow, targetColumn ) = 0;
        }
        accidentalDerotationChangesValidity |= !hypotheticalValid;
        accidentalDerotationChangesScience |=
            std::abs( doubleRotatedScience( targetRow, targetColumn ) - actualResidual ) > 1e-3F;
    }
    REQUIRE( accidentalDerotationChangesValidity );
    REQUIRE( accidentalDerotationChangesScience );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::P4Reduction<float, mx::improc::ADIDerotator<float, mx::verbose::vv>, mx::verbose::vv>
        doxygenReduction;
    doxygenReduction.reduce();
#endif
    // clang-format on
}

/// Verify P4Reduction preserves an excluded central signal and rejects complete local fits touched by the common mask.
/** \ingroup P4Reduction_unit_tests */
TEST_CASE( "P4 reduction central exclusion and common mask", "[P4Reduction][exclusion][mask][validity]" )
{
    OpenMPThreadGuard threads( 1 );

    SECTION( "kernel-support exclusion preserves an orthogonal injected signal" )
    {
        reductionHarness reduction;
        prepareReduction( reduction );
        const int targetRow = 15;
        const int targetColumn = 20;
        const float signal[]{ 1, -2, 1 };
        for( int image = 0; image < reduction.m_Nims; ++image )
        {
            reduction.m_tgtIms.image( image )( targetRow, targetColumn ) += signal[image];
        }

        REQUIRE( reduction.reduce() == 0 );
        for( int image = 0; image < reduction.m_Nims; ++image )
        {
            REQUIRE( reduction.m_psfsubValidity[0].image( image )( targetRow, targetColumn ) == 1 );
            REQUIRE( reduction.m_psfsub[0].image( image )( targetRow, targetColumn ) ==
                     Approx( signal[image] ).margin( 2e-4 ) );
        }
    }

    SECTION( "common mask rejects entire target and predictor-touched local fits" )
    {
        reductionHarness reduction;
        prepareReduction( reduction );
        reduction.m_mask.resize( reduction.m_Nrows, reduction.m_Ncols );
        reduction.m_mask.setOnes();
        reduction.m_mask( 15, 20 ) = 0;

        REQUIRE( reduction.reduce() == 0 );
        const auto &statistics = reduction.m_regionStatistics.front();
        REQUIRE( statistics.validLocalFitCount + statistics.maskedLocalFitCount == statistics.searchPixelCount );
        REQUIRE( statistics.maskedLocalFitCount > 1 );
        for( int image = 0; image < reduction.m_Nims; ++image )
        {
            REQUIRE( reduction.m_psfsubValidity[0].image( image )( 15, 20 ) == 0 );
            REQUIRE( mx::math::isNan( reduction.m_psfsub[0].image( image )( 15, 20 ) ) );
        }
    }

    SECTION( "rotated common mask rejects all-frame direct support without changing K" )
    {
        reductionHarness reduction;
        prepareReduction( reduction );
        reduction.m_regressionFrame = mx::improc::P4RegressionFrame::rotated;
        reduction.m_mask.resize( reduction.m_Nrows, reduction.m_Ncols );
        reduction.m_mask.setOnes();
        reduction.m_mask( 15, 20 ) = 0;

        REQUIRE( reduction.reduce() == 0 );
        const auto &statistics = reduction.m_regionStatistics.front();
        REQUIRE( statistics.maskedLocalFitCount == 0 );
        REQUIRE( statistics.supportInvalidLocalFitCount > 1 );
        REQUIRE( statistics.validLocalFitCount + statistics.supportInvalidLocalFitCount ==
                 statistics.searchPixelCount );
        for( int image = 0; image < reduction.m_Nims; ++image )
        {
            REQUIRE( reduction.m_psfsubValidity[0].image( image )( 15, 20 ) == 0 );
            REQUIRE( mx::math::isNan( reduction.m_psfsub[0].image( image )( 15, 20 ) ) );
        }
    }

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::P4Reduction<float, mx::improc::ADIDerotator<float, mx::verbose::vv>, mx::verbose::vv>
        doxygenReduction;
    doxygenReduction.reduce();
#endif
    // clang-format on
}

/// Verify P4Reduction maps fractions independently per annulus and marks only rank-unsupported output planes invalid.
/** \ingroup P4Reduction_unit_tests */
TEST_CASE( "P4 reduction annulus modes and rank invalidity", "[P4Reduction][annuli][modes][rank]" )
{
    OpenMPThreadGuard threads( 1 );
    reductionHarness reduction;
    prepareReduction( reduction, 5 );
    reduction.m_minRadius = { 5, 7 };
    reduction.m_maxRadius = { 6, 8 };
    reduction.m_modeFractions = { 0.2f, 0.4f };

    std::ostringstream warning;
    std::streambuf *const oldBuffer = std::cerr.rdbuf( warning.rdbuf() );
    const int result = reduction.reduce();
    std::cerr.rdbuf( oldBuffer );

    REQUIRE( result == 0 );
    REQUIRE( reduction.m_realizedModes == std::vector<std::vector<int>>{ { 1, 2 }, { 1, 2 } } );
    REQUIRE( reduction.m_regionStatistics.size() == 2 );
    std::size_t ownedPixels{ 0 };
    for( int column = 0; column < reduction.m_ownership.cols(); ++column )
    {
        for( int row = 0; row < reduction.m_ownership.rows(); ++row )
        {
            if( reduction.m_ownership( row, column ) >= 0 )
            {
                ++ownedPixels;
            }
        }
    }

    std::size_t expectedOwnedPixels{ 0 };
    for( const auto &statistics : reduction.m_regionStatistics )
    {
        expectedOwnedPixels += statistics.searchPixelCount;
        REQUIRE( statistics.validLocalFitCount + statistics.maskedLocalFitCount == statistics.searchPixelCount );
        REQUIRE( statistics.rankInvalidCounts[0] == 0 );
        REQUIRE( statistics.rankInvalidCounts[1] == statistics.validLocalFitCount );
    }
    REQUIRE( ownedPixels == expectedOwnedPixels );
    REQUIRE( warning.str().find( "above numerical rank" ) != std::string::npos );

    for( int image = 0; image < reduction.m_Nims; ++image )
    {
        REQUIRE( reduction.m_psfsubValidity[0].image( image ).sum() == Approx( ownedPixels ) );
        REQUIRE( reduction.m_psfsubValidity[1].image( image ).sum() == Approx( 0 ) );
    }

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::P4Reduction<float, mx::improc::ADIDerotator<float, mx::verbose::vv>, mx::verbose::vv>
        doxygenReduction;
    doxygenReduction.reduce();
#endif
    // clang-format on
}

/// Verify P4Reduction rejects invalid geometry, fraction mappings, cube/mask states, and image-edge footprints.
/** \ingroup P4Reduction_unit_tests */
TEST_CASE( "P4 reduction validation", "[P4Reduction][validation][edge]" )
{
    OpenMPThreadGuard threads( 1 );

    SECTION( "missing and mismatched required vectors" )
    {
        reductionHarness missing;
        REQUIRE_THROWS( missing.reduce() );

        reductionHarness mismatch;
        prepareReduction( mismatch );
        mismatch.m_maxRadius.push_back( 7 );
        REQUIRE_THROWS( mismatch.reduce() );

        reductionHarness missingFractions;
        prepareReduction( missingFractions );
        missingFractions.m_modeFractions.clear();
        REQUIRE_THROWS( missingFractions.reduce() );

        reductionHarness missingPolicy;
        prepareReduction( missingPolicy );
        missingPolicy.m_exclusionPolicy.reset();
        REQUIRE_THROWS( missingPolicy.reduce() );
    }

    SECTION( "overlapping or invalid annuli" )
    {
        reductionHarness overlap;
        prepareReduction( overlap );
        overlap.m_minRadius = { 5, 5.5f };
        overlap.m_maxRadius = { 6, 7 };
        REQUIRE_THROWS( overlap.reduce() );

        reductionHarness invalidRadius;
        prepareReduction( invalidRadius );
        invalidRadius.m_minRadius = { -1 };
        REQUIRE_THROWS( invalidRadius.reduce() );

        reductionHarness emptyAnnulus;
        prepareReduction( emptyAnnulus );
        emptyAnnulus.m_minRadius = { 0.1f };
        emptyAnnulus.m_maxRadius = { 0.2f };
        REQUIRE_THROWS( emptyAnnulus.reduce() );
    }

    SECTION( "invalid optimization geometry and policy" )
    {
        reductionHarness fullAnnulus;
        prepareReduction( fullAnnulus );
        fullAnnulus.m_orArcHalfWidth = 0;
        fullAnnulus.m_orMaxHalfAngle = 180;
        REQUIRE_NOTHROW( fullAnnulus.reduce() );

        reductionHarness invalidBuffer;
        prepareReduction( invalidBuffer );
        invalidBuffer.m_exclusionRadiusBuffer = -1;
        REQUIRE_THROWS( invalidBuffer.reduce() );

        reductionHarness nonfiniteGeometry;
        prepareReduction( nonfiniteGeometry );
        nonfiniteGeometry.m_orArcHalfWidth = std::numeric_limits<float>::quiet_NaN();
        REQUIRE_THROWS( nonfiniteGeometry.reduce() );

        reductionHarness invalidPolicy;
        prepareReduction( invalidPolicy );
        invalidPolicy.m_exclusionPolicy = static_cast<policyT>( 99 );
        REQUIRE_THROWS( invalidPolicy.reduce() );

        reductionHarness invalidFrame;
        prepareReduction( invalidFrame );
        invalidFrame.m_regressionFrame = static_cast<mx::improc::P4RegressionFrame>( 99 );
        REQUIRE_THROWS( invalidFrame.reduce() );

        reductionHarness negativeNumberImages;
        prepareReduction( negativeNumberImages );
        negativeNumberImages.m_numberImages = -1;
        REQUIRE_THROWS( negativeNumberImages.reduce() );

        reductionHarness rotatedPostMedian;
        prepareReduction( rotatedPostMedian );
        rotatedPostMedian.m_regressionFrame = mx::improc::P4RegressionFrame::rotated;
        rotatedPostMedian.m_postMedSub = true;
        REQUIRE_THROWS( rotatedPostMedian.reduce() );

        reductionHarness rotatedMultiImage;
        prepareReduction( rotatedMultiImage );
        rotatedMultiImage.m_regressionFrame = mx::improc::P4RegressionFrame::rotated;
        rotatedMultiImage.m_numberImages = 1;
        REQUIRE_THROWS_WITH(
            rotatedMultiImage.reduce(),
            Catch::Matchers::Contains( "p4.numberImages is supported only for detector-frame P4 regression" ) );

        reductionHarness nonfiniteDerotation;
        prepareReduction( nonfiniteDerotation );
        nonfiniteDerotation.m_regressionFrame = mx::improc::P4RegressionFrame::rotated;
        nonfiniteDerotation.m_derotF.m_angles[1] = std::numeric_limits<float>::quiet_NaN();
        REQUIRE_THROWS( nonfiniteDerotation.reduce() );

        reductionHarness invalidTolerance;
        prepareReduction( invalidTolerance );
        invalidTolerance.m_rankTolerance = 1;
        REQUIRE_THROWS( invalidTolerance.reduce() );
    }

    SECTION( "invalid fractions and realized mode collisions" )
    {
        reductionHarness invalidFraction;
        prepareReduction( invalidFraction );
        invalidFraction.m_modeFractions = { 0.5f, 0.25f };
        REQUIRE_THROWS( invalidFraction.reduce() );

        reductionHarness zeroMode;
        prepareReduction( zeroMode );
        zeroMode.m_modeFractions = { 0.2f };
        REQUIRE_THROWS( zeroMode.reduce() );

        reductionHarness duplicateMode;
        prepareReduction( duplicateMode );
        duplicateMode.m_modeFractions = { 0.5f, 0.6f };
        REQUIRE_THROWS( duplicateMode.reduce() );

        reductionHarness oneFrameRotated;
        prepareReduction( oneFrameRotated, 1 );
        oneFrameRotated.m_regressionFrame = mx::improc::P4RegressionFrame::rotated;
        REQUIRE_THROWS( oneFrameRotated.reduce() );
    }

    SECTION( "invalid cube and mask state" )
    {
        reductionHarness invalidCube;
        prepareReduction( invalidCube );
        invalidCube.m_Nrows += 1;
        REQUIRE_THROWS( invalidCube.reduce() );

        reductionHarness invalidMask;
        prepareReduction( invalidMask );
        invalidMask.m_mask.resize( 2, 2 );
        REQUIRE_THROWS( invalidMask.reduce() );
    }

    SECTION( "complete OR interpolation footprint crosses an edge" )
    {
        reductionHarness edge;
        prepareReduction( edge );
        edge.m_minRadius = { 13 };
        edge.m_maxRadius = { 14 };
        REQUIRE_THROWS( edge.reduce() );
    }

    SECTION( "diagnostic destination is required only when enabled" )
    {
        reductionHarness diagnostics;
        prepareReduction( diagnostics );
        diagnostics.m_writeDiagnostics = true;
        diagnostics.m_diagnosticDirectory.clear();
        REQUIRE_THROWS( diagnostics.reduce() );
    }

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::P4Reduction<float, mx::improc::ADIDerotator<float, mx::verbose::vv>, mx::verbose::vv>
        doxygenReduction;
    doxygenReduction.reduce();
#endif
    // clang-format on
}

/// Verify P4Reduction owns target discovery, rejects every RDI source, enforces finite data, and honors
/// preprocess-only.
/** \ingroup P4Reduction_unit_tests */
TEST_CASE( "P4 reduction input lifecycle", "[P4Reduction][input][RDI][finite][preprocess]" )
{
    OpenMPThreadGuard threads( 1 );

    SECTION( "all configured or preloaded RDI sources are rejected before target discovery" )
    {
        auto requireRDIRejected = []( auto configureRDI )
        {
            reductionHarness reduction;
            prepareReduction( reduction );
            reduction.m_filesRead = false;
            reduction.m_fileList.clear();
            reduction.m_directory = "/definitely/not/a/P4/target/directory";
            configureRDI( reduction );
            REQUIRE_THROWS( reduction.reduce() );
            REQUIRE( reduction.m_fileList.empty() );
        };

        requireRDIRejected( []( reductionHarness &reduction ) { reduction.m_RDIfileList = { "reference.fits" }; } );
        requireRDIRejected( []( reductionHarness &reduction ) { reduction.m_RDIfileListFile = "references.txt"; } );
        requireRDIRejected( []( reductionHarness &reduction ) { reduction.m_RDIdirectory = "references"; } );
        requireRDIRejected( []( reductionHarness &reduction ) { reduction.m_RDIprefix = "rdi_"; } );
        requireRDIRejected( []( reductionHarness &reduction ) { reduction.m_refIms.resize( 2, 2, 1 ); } );

        reductionHarness explicitRegions;
        prepareReduction( explicitRegions );
        explicitRegions.m_RDIprefix = "rdi_";
        REQUIRE_THROWS( explicitRegions.regions( { 5 }, { 6 } ) );
    }

    SECTION( "configured target discovery failure is reported" )
    {
        reductionHarness reduction;
        prepareReduction( reduction );
        reduction.m_filesRead = false;
        reduction.m_fileList.clear();
        reduction.m_dateKeyword.clear();
        reduction.m_directory = "/definitely/not/a/P4/target/directory";
        REQUIRE_THROWS( reduction.reduce() );
    }

    SECTION( "configured target discovery owns the complete finite-input run" )
    {
        TestDirectory directory;
        mx::fits::fitsFile<float, mx::verbose::vv> writer;
        reductionT::fitsHeaderT inputHeader;
        inputHeader.append<float>( "ANGLE", 0, "test parallactic angle" );
        for( int imageIndex = 0; imageIndex < 3; ++imageIndex )
        {
            reductionT::imageT image = reductionT::imageT::Constant( 31, 31, static_cast<float>( imageIndex + 1 ) );
            REQUIRE( writer.write( directory.file( "target_" + std::to_string( imageIndex ) + ".fits" ).string(),
                                   image,
                                   inputHeader ) == mx::error_t::noerror );
        }

        reductionHarness reduction;
        prepareReduction( reduction );
        reduction.m_filesRead = false;
        reduction.m_fileList.clear();
        reduction.m_directory = directory.path().string();
        reduction.m_prefix = "target_";
        reduction.m_dateKeyword.clear();
        reduction.m_imSize = 0;
        reduction.m_derotF.angleKeyword( "ANGLE" );
        reduction.m_derotF.m_angleScale = 1;
        REQUIRE( reduction.reduce() == 0 );
        REQUIRE( reduction.m_fileList.size() == 3 );
        REQUIRE( reduction.m_postReadCalled );
        REQUIRE( reduction.m_psfsub.size() == 1 );
    }

    SECTION( "strict input policy rejects pristine nonfinite FITS data" )
    {
        TestDirectory directory;
        reductionHarness reduction;
        prepareReduction( reduction );
        reduction.m_filesRead = false;
        reduction.m_fileList.clear();
        reduction.m_dateKeyword.clear();
        reductionT::imageT image( 31, 31 );
        image.setOnes();
        image( 15, 20 ) = std::numeric_limits<float>::quiet_NaN();
        mx::fits::fitsFile<float, mx::verbose::vv> writer;
        reductionT::fitsHeaderT inputHeader;
        inputHeader.append<float>( "ANGLE", 0, "test parallactic angle" );
        REQUIRE( writer.write( directory.file( "nonfinite.fits" ).string(), image, inputHeader ) ==
                 mx::error_t::noerror );
        reduction.m_fileList = { directory.file( "nonfinite.fits" ).string() };
        reduction.m_derotF.angleKeyword( "ANGLE" );
        reduction.m_derotF.m_angleScale = 1;
        REQUIRE_THROWS( reduction.reduce() );
        REQUIRE_FALSE( reduction.m_postReadCalled );
    }

    SECTION( "post-preprocess target state must remain finite" )
    {
        reductionHarness reduction;
        prepareReduction( reduction );
        reduction.m_tgtIms.image( 1 )( 15, 20 ) = std::numeric_limits<float>::infinity();
        REQUIRE_THROWS( reduction.reduce() );
    }

    SECTION( "preprocess-only stop leaves no stale reduction products" )
    {
        reductionHarness reduction;
        prepareReduction( reduction );
        reduction.m_psfsub.resize( 1 );
        reduction.m_psfsubValidity.resize( 1 );
        reduction.m_preProcess_only = true;
        reduction.m_skipPreProcess = false;
        reduction.m_minRadius.clear();
        reduction.m_maxRadius.clear();
        reduction.m_modeFractions.clear();
        reduction.m_orDeltaRadiusInner = std::numeric_limits<float>::quiet_NaN();
        reduction.m_orDeltaRadiusOuter = std::numeric_limits<float>::quiet_NaN();
        reduction.m_orArcHalfWidth = std::numeric_limits<float>::quiet_NaN();
        reduction.m_orMaxHalfAngle = std::numeric_limits<float>::quiet_NaN();
        reduction.m_psfRadius = std::numeric_limits<float>::quiet_NaN();
        reduction.m_exclusionPolicy.reset();
        REQUIRE( reduction.reduce() == 0 );
        REQUIRE( reduction.m_psfsub.empty() );
        REQUIRE( reduction.m_psfsubValidity.empty() );
        REQUIRE( reduction.m_regionStatistics.empty() );
    }

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::P4Reduction<float, mx::improc::ADIDerotator<float, mx::verbose::vv>, mx::verbose::vv>
        doxygenReduction;
    doxygenReduction.reduce();
    doxygenReduction.regions( { 5 }, { 6 } );
#endif
    // clang-format on
}

/// Verify P4Reduction progress, OpenMP determinism, and serial rethrow of worker failures.
/** \ingroup P4Reduction_unit_tests */
TEST_CASE( "P4 reduction OpenMP determinism and exception capture", "[P4Reduction][OpenMP][exception]" )
{
    SECTION( "one and many workers are deterministic" )
    {
        reductionHarness serial;
        prepareReduction( serial, 5 );
        serial.m_minRadius = { 5, 7 };
        serial.m_maxRadius = { 6, 8 };
        serial.m_modeFractions = { 0.2f, 0.4f };
        std::string serialProgress;
        {
            OpenMPThreadGuard threads( 1 );
            CerrCapture capture;
            REQUIRE( serial.reduce() == 0 );
            serialProgress = capture.str();
        }

        reductionHarness parallel;
        prepareReduction( parallel, 5 );
        parallel.m_minRadius = { 5, 7 };
        parallel.m_maxRadius = { 6, 8 };
        parallel.m_modeFractions = { 0.2f, 0.4f };
        std::string parallelProgress;
        {
            OpenMPThreadGuard threads( 3 );
            CerrCapture capture;
            REQUIRE( parallel.reduce() == 0 );
            parallelProgress = capture.str();
        }
        requireSameReduction( serial, parallel );

        const std::size_t totalSearchPixels =
            serial.m_regionStatistics[0].searchPixelCount + serial.m_regionStatistics[1].searchPixelCount;
        const std::size_t firstRegionSearchPixels = serial.m_regionStatistics[0].searchPixelCount;
        for( const std::string *progress : { &serialProgress, &parallelProgress } )
        {
            REQUIRE( progress->find( "P4 geometry annulus 1 / 2 [5,6)... done:" ) != std::string::npos );
            REQUIRE( progress->find( "P4 geometry annulus 2 / 2 [7,8)... done:" ) != std::string::npos );
            REQUIRE( progress->find( "P4 annulus 1 / 2 [5,6), K=" ) != std::string::npos );
            REQUIRE( progress->find( "P4 annulus 2 / 2 [7,8), K=" ) != std::string::npos );
            REQUIRE( progress->find( "| overall " ) != std::string::npos );
            REQUIRE( progress->find( "s/loop" ) != std::string::npos );
            REQUIRE( progress->find( "s left" ) != std::string::npos );
            REQUIRE( progress->find( "complete | overall " + std::to_string( firstRegionSearchPixels ) + " / " +
                                     std::to_string( totalSearchPixels ) ) != std::string::npos );
            REQUIRE( progress->find( "complete | overall " + std::to_string( totalSearchPixels ) + " / " +
                                     std::to_string( totalSearchPixels ) + " (100%)" ) != std::string::npos );
            REQUIRE( progress->find( "P4 regression complete: " + std::to_string( totalSearchPixels ) + " / " +
                                     std::to_string( totalSearchPixels ) + " search pixels in " ) !=
                     std::string::npos );
            REQUIRE( progress->back() == '\n' );

            std::size_t progressRecordCount{ 0 };
            std::size_t completeRecordCount{ 0 };
            std::size_t recordBegin{ 0 };
            while( recordBegin < progress->size() )
            {
                const std::size_t recordEnd = progress->find_first_of( "\r\n", recordBegin );
                const std::string record = progress->substr( recordBegin, recordEnd - recordBegin );
                if( record.starts_with( "P4 annulus" ) )
                {
                    ++progressRecordCount;
                    const std::size_t overall = record.find( "| overall " );
                    REQUIRE( overall != std::string::npos );
                    REQUIRE( record.find( "| overall ", overall + 1 ) == std::string::npos );
                    REQUIRE( record.find( "P4 annulus", 1 ) == std::string::npos );
                    if( record.find( "complete | overall" ) != std::string::npos )
                    {
                        ++completeRecordCount;
                    }
                }
                if( recordEnd == std::string::npos )
                {
                    break;
                }
                recordBegin = recordEnd + 1;
            }
            REQUIRE( progressRecordCount >= 6 );
            REQUIRE( completeRecordCount == 2 );
        }
    }

    SECTION( "rotated regression is deterministic across worker counts" )
    {
        const auto prepareRotated = []( reductionHarness &reduction )
        {
            prepareReduction( reduction, 5 );
            reduction.m_minRadius = { 5, 7 };
            reduction.m_maxRadius = { 6, 8 };
            reduction.m_modeFractions = { 0.25F, 0.5F };
            reduction.m_regressionFrame = mx::improc::P4RegressionFrame::rotated;
            reduction.m_derotF.m_angles = { 0, 17, -11, 33, -27 };
            reduction.m_doDerotate = true;
            for( int image = 0; image < reduction.m_Nims; ++image )
            {
                for( int column = 0; column < reduction.m_Ncols; ++column )
                {
                    for( int row = 0; row < reduction.m_Nrows; ++row )
                    {
                        reduction.m_tgtIms.image( image )( row, column ) =
                            0.02F * static_cast<float>( row * row ) + 0.03F * static_cast<float>( column ) +
                            static_cast<float>( image + 1 ) *
                                ( 0.1F * static_cast<float>( row ) - 0.07F * static_cast<float>( column ) ) +
                            static_cast<float>( image % 2 ) * 0.0005F * static_cast<float>( row * column );
                    }
                }
            }
        };

        reductionHarness serial;
        prepareRotated( serial );
        std::string serialProgress;
        {
            OpenMPThreadGuard threads( 1 );
            CerrCapture capture;
            REQUIRE( serial.reduce() == 0 );
            serialProgress = capture.str();
        }

        reductionHarness parallel;
        prepareRotated( parallel );
        std::string parallelProgress;
        {
            OpenMPThreadGuard threads( 3 );
            CerrCapture capture;
            REQUIRE( parallel.reduce() == 0 );
            parallelProgress = capture.str();
        }
        requireSameReduction( serial, parallel );

        for( const std::string *progress : { &serialProgress, &parallelProgress } )
        {
            REQUIRE( progress->find( "P4 rotated geometry annulus 1 / 2 [5,6)... done:" ) != std::string::npos );
            REQUIRE( progress->find( "P4 rotated geometry annulus 2 / 2 [7,8)... done:" ) != std::string::npos );
            REQUIRE( progress->find( "P4 rotated annulus 1 / 2 [5,6), K=" ) != std::string::npos );
            REQUIRE( progress->find( "P4 rotated annulus 2 / 2 [7,8), K=" ) != std::string::npos );
            REQUIRE( progress->find( "P4 regression complete:" ) != std::string::npos );
            REQUIRE( progress->back() == '\n' );
        }
    }

    SECTION( "worker eigensolver exception is captured and rethrown outside OpenMP" )
    {
        OpenMPThreadGuard threads( 3 );
        reductionHarness reduction;
        prepareReduction( reduction );
        eigenSolverReset solver;
        std::string failedProgress;
        {
            CerrCapture capture;
            REQUIRE_THROWS( reduction.reduce() );
            failedProgress = capture.str();
        }
        REQUIRE( failedProgress.find( "P4 annulus 1 / 1 [5,6), K=" ) != std::string::npos );
        REQUIRE( failedProgress.find( "interrupted" ) != std::string::npos );
        REQUIRE( failedProgress.find( "P4 regression complete:" ) == std::string::npos );
        REQUIRE( failedProgress.back() == '\n' );
    }

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::P4Reduction<float, mx::improc::ADIDerotator<float, mx::verbose::vv>, mx::verbose::vv>
        doxygenReduction;
    doxygenReduction.reduce();
#endif
    // clang-format on
}

/// Verify P4Reduction diagnostics are opt-in, atomically published, science-neutral, and complete in FITS provenance.
/** \ingroup P4Reduction_unit_tests */
TEST_CASE( "P4 reduction diagnostics and provenance", "[P4Reduction][diagnostics][header][finalProcess]" )
{
    OpenMPThreadGuard threads( 1 );
    TestDirectory directory;

    reductionHarness baseline;
    prepareReduction( baseline );
    baseline.m_diagnosticDirectory = directory.file( "disabled" ).string();
    REQUIRE( baseline.reduce() == 0 );
    REQUIRE_FALSE( std::filesystem::exists( directory.file( "disabled" ) ) );

    reductionHarness diagnostic;
    prepareReduction( diagnostic );
    diagnostic.m_writeDiagnostics = true;
    diagnostic.m_diagnosticDirectory = directory.file( "diagnostics" ).string();
    REQUIRE( diagnostic.reduce() == 0 );
    requireSameReduction( baseline, diagnostic );
    REQUIRE( std::filesystem::exists( directory.file( "diagnostics/p4Ownership.fits" ) ) );
    REQUIRE( std::filesystem::exists( directory.file( "diagnostics/p4Validity_000.fits" ) ) );
    REQUIRE( std::filesystem::exists( directory.file( "diagnostics/p4TemporalSelection_000.fits" ) ) );
    REQUIRE( std::filesystem::exists( directory.file( "diagnostics/p4CanonicalOR_000.fits" ) ) );
    REQUIRE( std::filesystem::exists( directory.file( "diagnostics/p4RegionSummary.fits" ) ) );
    REQUIRE( std::filesystem::exists( directory.file( "diagnostics/p4Timing.fits" ) ) );

    reductionT::imageT summary;
    reductionT::fitsHeaderT summaryHeader;
    mx::fits::fitsFile<float, mx::verbose::vv> reader;
    REQUIRE( reader.read( summary, summaryHeader, directory.file( "diagnostics/p4RegionSummary.fits" ).string() ) ==
             mx::error_t::noerror );
    REQUIRE( summaryHeader["P4 FRAME"].String().starts_with( "detector" ) );

    mx::improc::eigenCube<float> diagnosticValidity;
    reductionT::fitsHeaderT diagnosticValidityHeader;
    REQUIRE( reader.read( diagnosticValidity,
                          diagnosticValidityHeader,
                          directory.file( "diagnostics/p4Validity_000.fits" ).string() ) == mx::error_t::noerror );
    REQUIRE( diagnosticValidityHeader["P4 FRAME"].String().starts_with( "detector" ) );
    REQUIRE( diagnosticValidity.planes() == diagnostic.m_Nims );
    REQUIRE( diagnosticValidity.image( 0 ).isApprox( diagnostic.m_psfsubValidity[0].image( 0 ), 0 ) );

    REQUIRE( summary.rows() == 1 );
    REQUIRE( summary.cols() == 14 );
    REQUIRE( summary( 0, 9 ) + summary( 0, 10 ) + summary( 0, 11 ) == Approx( summary( 0, 2 ) ) );
    REQUIRE( summary( 0, 12 ) == Approx( diagnostic.m_realizedModes[0][0] ) );
    REQUIRE( summary( 0, 13 ) == Approx( diagnostic.m_regionStatistics[0].rankInvalidCounts[0] ) );

    reductionT::imageT timing;
    reductionT::fitsHeaderT timingHeader;
    REQUIRE( reader.read( timing, timingHeader, directory.file( "diagnostics/p4Timing.fits" ).string() ) ==
             mx::error_t::noerror );
    REQUIRE( timing.rows() == 1 );
    REQUIRE( timing.cols() == 8 );
    REQUIRE( timingHeader["P4 TIMING SCHEMA"].Int() == 3 );
    REQUIRE( timingHeader["P4 TIMING COLUMNS"].String().starts_with( "geometryElapsed" ) );
    REQUIRE( timingHeader["P4 TIMING COLUMNS"].String().find( "sameImageSamplingWorker" ) != std::string::npos );
    REQUIRE( timingHeader["P4 TIMING COLUMNS"].String().find( "temporalSamplingWorker" ) != std::string::npos );
    for( Eigen::Index column = 0; column < timing.cols(); ++column )
    {
        REQUIRE( std::isfinite( timing( 0, column ) ) );
        REQUIRE( timing( 0, column ) >= 0 );
    }
    REQUIRE( timing( 0, 2 ) == Approx( timing( 0, 3 ) + timing( 0, 4 ) ).margin( 1e-6 ) );

    reductionT::fitsHeaderT header;
    diagnostic.appendReductionHeader( header );
    REQUIRE( header["P4 ALGORITHM"].String().starts_with( "P4-PCA" ) );
    REQUIRE( header["P4 FRAME"].String().starts_with( "detector" ) );
    REQUIRE( header["P4 IN SAMPLE"].Int() == 1 );
    REQUIRE( header["P4 RDI"].Int() == 0 );
    REQUIRE( header["P4 NUMBER IMAGES"].Int() == 0 );
    REQUIRE( header["P4 EXCLUSION POLICY"].String().starts_with( "kernelSupport" ) );
    REQUIRE( header["P4 MODE FRACTIONS"].String().starts_with( "0.5" ) );
    REQUIRE( header["P4 REALIZED MODES 000"].String().starts_with( "1" ) );
    REQUIRE( header["P4 VALID FIT COUNT"].String() ==
             std::to_string( diagnostic.m_regionStatistics[0].validLocalFitCount ) );
    REQUIRE( header["P4 MASK INVALID FIT COUNT"].String() ==
             std::to_string( diagnostic.m_regionStatistics[0].maskedLocalFitCount ) );
    REQUIRE( header["P4 SUPPORT INVALID FIT COUNT"].String() ==
             std::to_string( diagnostic.m_regionStatistics[0].supportInvalidLocalFitCount ) );
    diagnostic.m_heads.resize( diagnostic.m_Nims );
    diagnostic.m_doOutputPSFSub = true;
    diagnostic.m_outputDir = directory.path().string();
    diagnostic.m_PSFSubPrefix = "p4-psf";
    REQUIRE( diagnostic.finalProcess() == 0 );
    REQUIRE( std::filesystem::exists( directory.file( "p4-psf_000_00000.fits" ) ) );
    reductionT::imageT persistedScience;
    reductionT::fitsHeaderT persistedHeader;
    REQUIRE( reader.read( persistedScience, persistedHeader, directory.file( "p4-psf_000_00000.fits" ).string() ) ==
             mx::error_t::noerror );
    REQUIRE( persistedHeader["P4 ALGORITHM"].String().starts_with( "P4-PCA" ) );
    REQUIRE( persistedHeader["P4 FRAME"].String().starts_with( "detector" ) );

    reductionHarness atomic;
    prepareReduction( atomic );
    atomic.m_regressionFrame = mx::improc::P4RegressionFrame::rotated;
    atomic.m_writeDiagnostics = true;
    atomic.m_diagnosticDirectory = directory.file( "atomic" ).string();
    reductionT::imageT product = reductionT::imageT::Constant( 2, 2, 1 );
    mx::improc::P4ReductionTestAccess::writeDiagnostic( atomic, "replace.fits", product );
    product.setConstant( 2 );
    mx::improc::P4ReductionTestAccess::writeDiagnostic( atomic, "replace.fits", product );
    reductionT::imageT replaced;
    reductionT::fitsHeaderT replacedHeader;
    REQUIRE( reader.read( replaced, replacedHeader, directory.file( "atomic/replace.fits" ).string() ) ==
             mx::error_t::noerror );
    REQUIRE( replaced.isConstant( 2 ) );
    REQUIRE( replacedHeader["P4 FRAME"].String().starts_with( "rotated" ) );

    reductionHarness rotatedPostMedian;
    prepareReduction( rotatedPostMedian );
    rotatedPostMedian.m_regressionFrame = mx::improc::P4RegressionFrame::rotated;
    rotatedPostMedian.m_postMedSub = true;
    REQUIRE_THROWS( rotatedPostMedian.finalProcess() );

    reductionHarness invalidFrame;
    prepareReduction( invalidFrame );
    invalidFrame.m_regressionFrame = static_cast<mx::improc::P4RegressionFrame>( 255 );
    REQUIRE_THROWS( invalidFrame.finalProcess() );

    std::filesystem::create_directory( directory.file( "atomic/occupied.fits" ) );
    REQUIRE_THROWS( mx::improc::P4ReductionTestAccess::writeDiagnostic( atomic, "occupied.fits", product ) );
    for( const auto &entry : std::filesystem::directory_iterator( directory.file( "atomic" ) ) )
    {
        REQUIRE( entry.path().filename().string().find( "occupied.fits.tmp" ) == std::string::npos );
    }

    reductionHarness unavailableDirectory;
    prepareReduction( unavailableDirectory );
    unavailableDirectory.m_writeDiagnostics = true;
    unavailableDirectory.m_diagnosticDirectory = directory.file( "regular-file" ).string();
    writeTextFile( directory.file( "regular-file" ), "not a directory" );
    REQUIRE_THROWS(
        mx::improc::P4ReductionTestAccess::writeDiagnostic( unavailableDirectory, "product.fits", product ) );

    reductionHarness nestedParent;
    prepareReduction( nestedParent );
    nestedParent.m_writeDiagnostics = true;
    nestedParent.m_diagnosticDirectory = ".";
    const std::string nestedProduct = directory.file( "missing/parent/product.fits" ).string();
    REQUIRE_NOTHROW( mx::improc::P4ReductionTestAccess::writeDiagnostic( nestedParent, nestedProduct, product ) );
    REQUIRE( std::filesystem::exists( nestedProduct ) );

    reductionHarness disabledWriter;
    prepareReduction( disabledWriter );
    disabledWriter.m_diagnosticDirectory = directory.file( "disabled-writer" ).string();
    REQUIRE_NOTHROW(
        mx::improc::P4ReductionTestAccess::writeDiagnostic( disabledWriter, "not-written.fits", product ) );
    REQUIRE_FALSE( std::filesystem::exists( directory.file( "disabled-writer" ) ) );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::P4Reduction<float, mx::improc::ADIDerotator<float, mx::verbose::vv>, mx::verbose::vv>
        doxygenReduction;
    mx::improc::P4Reduction<float,
                            mx::improc::ADIDerotator<float, mx::verbose::vv>,
                            mx::verbose::vv>::fitsHeaderT doxygenHeader;
    mx::improc::P4Reduction<float,
                            mx::improc::ADIDerotator<float, mx::verbose::vv>,
                            mx::verbose::vv>::imageT doxygenDiagnostic;
    doxygenReduction.appendReductionHeader( doxygenHeader );
    doxygenReduction.writeDiagnostic( "diagnostic.fits", doxygenDiagnostic );
    doxygenReduction.finalProcess();
#endif
    // clang-format on
}

} // namespace P4Reduction_test
} // namespace unitTest
