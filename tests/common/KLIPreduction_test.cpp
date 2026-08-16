/** \file KLIPreduction_test.cpp
 * \brief Tests KLIP reduction centering and normalization behavior.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "HCIobservation_test_fixture.hpp"
#include "src/common/ADIDerotator.hpp"
#include "src/common/KLIPreduction.hpp"

#include <algorithm>
#include <array>
#include <cmath>

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
    using reductionT::m_preProcess_only;
    using reductionT::m_psfsub;
    using reductionT::m_PSFSubPrefix;
    using reductionT::m_RDIfileList;
    using reductionT::m_RDIfilesRead;
    using reductionT::m_RDIimageMJD;
    using reductionT::m_refIms;
    using reductionT::m_skipPreProcess;
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

/// Compare two Eigen-like arrays coefficient by coefficient.
template <typename actualT, typename expectedT>
void requireApprox( const actualT &actual,     /**< [in] values produced by the adaptive KL basis */
                    const expectedT &expected, /**< [in] reference values */
                    double tolerance = 2e-5 /**< [in] absolute comparison tolerance */ )
{
    REQUIRE( actual.rows() == expected.rows() );
    REQUIRE( actual.cols() == expected.cols() );
    for( Eigen::Index column = 0; column < actual.cols(); ++column )
    {
        for( Eigen::Index row = 0; row < actual.rows(); ++row )
        {
            REQUIRE( actual( row, column ) == Approx( expected( row, column ) ).margin( tolerance ) );
        }
    }
}

/// Form a spatial projector from the largest rows of an ascending KL-mode matrix.
reductionT::imageT modeProjector( const reductionT::imageT &modes, /**< [in] ascending KL modes stored by row */
                                  int modeCount /**< [in] number of largest modes to retain */ )
{
    if( modeCount <= 0 || modeCount > modes.rows() )
    {
        throw std::invalid_argument( "invalid test projector mode count" );
    }

    const Eigen::MatrixXf selected = modes.matrix().bottomRows( modeCount );
    return ( selected.transpose() * selected ).array();
}

/// Form the direct thin-SVD spatial projector for a vectorized reference library.
reductionT::imageT svdProjector( const reductionT::imageT &references, /**< [in] pixels-by-references matrix */
                                 int modeCount /**< [in] number of largest singular modes to retain */ )
{
    Eigen::JacobiSVD<Eigen::MatrixXf> decomposition( references.matrix(), Eigen::ComputeThinU );
    const Eigen::MatrixXf selected = decomposition.matrixU().leftCols( modeCount );
    return ( selected * selected.transpose() ).array();
}

/// Form a target residual from a spatial projector.
reductionT::imageT projectedResidual( const reductionT::imageT &projector, /**< [in] spatial projector */
                                      const reductionT::imageT &target /**< [in] target column */ )
{
    return ( target.matrix() - projector.matrix() * target.matrix() ).array();
}

/// Initialize the direct-worker state shared by adaptive-basis tests.
void prepareWorkerReduction( reductionHarness &reduction, /**< [out] configured reduction harness */
                             int pixelCount,              /**< [in] flattened region size */
                             int targetCount,             /**< [in] number of target images */
                             int referenceCount,          /**< [in] number of reference images */
                             const std::vector<int> &modeCounts /**< [in] configured KL mode counts */ )
{
    reduction.m_meanSubMethod = mx::improc::HCI::meanSub::none;
    reduction.m_pixelTSNormMethod = mx::improc::HCI::pixelTSNorm::none;
    reduction.m_includeMethod = mx::improc::HCI::include::all;
    reduction.m_includeRefNum = 0;
    reduction.m_excludeMethod = mx::improc::HCI::exclude::none;
    reduction.m_excludeMethodMax = mx::improc::HCI::exclude::none;
    reduction.m_Nrows = 1;
    reduction.m_Ncols = pixelCount;
    reduction.m_Npix = pixelCount;
    reduction.m_Nims = targetCount;
    reduction.m_Nmodes = modeCounts;
    reduction.m_maxNmodes = *std::max_element( modeCounts.begin(), modeCounts.end() );
    reduction.m_psfsub.resize( modeCounts.size() );
    for( auto &output : reduction.m_psfsub )
    {
        output.resize( pixelCount, 1, targetCount );
        output.cube().setZero();
    }
    reduction.m_imsIncluded.resize( targetCount, referenceCount );
    reduction.m_imsIncluded.setOnes();
}

/// Reject an mxlib eigensolver workspace allocation for deterministic failure testing.
void *rejectEigenAllocation( std::size_t /**< [in] ignored requested byte count */ )
{
    return nullptr;
}

/// Restore mxlib's eigensolver allocation hook after a scoped failure test.
class EigenAllocatorGuard
{
  public:
    /// Install the deterministic allocation failure hook.
    EigenAllocatorGuard() : m_previous( mx::math::detail::eigenLapackTestHooks<double>::allocator )
    {
        mx::math::detail::eigenLapackTestHooks<double>::allocator = rejectEigenAllocation;
    }

    /// Restore the allocator that was active before construction.
    ~EigenAllocatorGuard()
    {
        mx::math::detail::eigenLapackTestHooks<double>::allocator = m_previous;
    }

    /// Scoped allocator guards cannot be copied.
    EigenAllocatorGuard( const EigenAllocatorGuard & ) = delete;

    /// Scoped allocator guards cannot be assigned.
    EigenAllocatorGuard &operator=( const EigenAllocatorGuard & ) = delete;

  private:
    /// Allocator hook restored at scope exit.
    mx::math::detail::eigenLapackTestHooks<double>::allocatorT m_previous;
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

    mx::app::appConfigurator config;
    reductionT badCurrentMean;
    badCurrentMean.setupConfig( config );
    badCurrentMean.m_meanSubMethod = static_cast<mx::improc::HCI::meanSub>( 99 );
    REQUIRE_THROWS( badCurrentMean.loadConfig( config ) );

    config = mx::app::appConfigurator{};
    reductionT badCurrentNorm;
    badCurrentNorm.setupConfig( config );
    badCurrentNorm.m_pixelTSNormMethod = static_cast<mx::improc::HCI::pixelTSNorm>( 99 );
    REQUIRE_THROWS( badCurrentNorm.loadConfig( config ) );

    config = mx::app::appConfigurator{};
    reductionT badCurrentInclusion;
    badCurrentInclusion.setupConfig( config );
    badCurrentInclusion.m_includeMethod = static_cast<mx::improc::HCI::include>( 99 );
    REQUIRE_THROWS( badCurrentInclusion.loadConfig( config ) );
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

    writeTextFile( directory.file( "not-a-directory" ), "blocking file\n" );
    reduction.m_diagnosticDirectory = directory.file( "not-a-directory/child" ).string();
    REQUIRE_THROWS( reduction.writeDiagnostic( "cv.fits", image ) );

    reduction.m_diagnosticDirectory = ".";
    REQUIRE_THROWS( reduction.writeDiagnostic( directory.path().string(), image ) );
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

/// Verify KLIPreduction::meanSubtract reapplies a mask after reference-image centering.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP masked reference-image centering", "[KLIPreduction][meanSubtract][reference][mask]" )
{
    reductionT reduction;
    reduction.m_meanSubMethod = mx::improc::HCI::meanSub::meanImage;
    reduction.m_pixelTSNormMethod = mx::improc::HCI::pixelTSNorm::none;

    mx::improc::eigenCube<float> references( 1, 3, 2 );
    references.image( 0 ) << 1, 100, 3;
    references.image( 1 ) << 3, 200, 7;
    mx::improc::eigenCube<float> targets( 1, 3, 1 );
    targets.image( 0 ) << 12, 300, 15;
    reductionT::imageT mask( 1, 3 );
    mask << 1, 0, 1;
    std::vector<float> norms;

    reduction.meanSubtract( references, targets, mask, norms );

    REQUIRE( references.image( 0 ).isApprox( ( reductionT::imageT( 1, 3 ) << -1, 0, -2 ).finished() ) );
    REQUIRE( references.image( 1 ).isApprox( ( reductionT::imageT( 1, 3 ) << 1, 0, 2 ).finished() ) );
    REQUIRE( targets.image( 0 ).isApprox( ( reductionT::imageT( 1, 3 ) << 10, 0, 10 ).finished() ) );
    REQUIRE( norms.size() == 2 );
    REQUIRE( norms[0] == Approx( std::sqrt( 0.5F ) ) );
    REQUIRE( norms[1] == Approx( std::sqrt( 0.5F ) ) );
}

/// Verify KLIPreduction::meanSubtract independently centers unmasked target images for both per-image methods.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP unmasked target-image centering", "[KLIPreduction][meanSubtract][target]" )
{
    for( const mx::improc::HCI::meanSub method :
         { mx::improc::HCI::meanSub::imageMean, mx::improc::HCI::meanSub::imageMedian } )
    {
        reductionT reduction;
        reduction.m_meanSubMethod = method;
        reduction.m_pixelTSNormMethod = mx::improc::HCI::pixelTSNorm::none;

        mx::improc::eigenCube<float> references( 1, 3, 1 );
        references.image( 0 ) << 1, 2, 4;
        mx::improc::eigenCube<float> targets( 1, 3, 1 );
        targets.image( 0 ) << 9, 10, 12;
        reductionT::imageT mask;
        std::vector<float> norms;

        reduction.meanSubtract( references, targets, mask, norms );

        if( method == mx::improc::HCI::meanSub::imageMean )
        {
            const reductionT::imageT expected =
                ( reductionT::imageT( 1, 3 ) << -4.0F / 3.0F, -1.0F / 3.0F, 5.0F / 3.0F ).finished();
            REQUIRE( references.image( 0 ).isApprox( expected ) );
            REQUIRE( targets.image( 0 ).isApprox( expected ) );
        }
        else
        {
            const reductionT::imageT expected = ( reductionT::imageT( 1, 3 ) << -1, 0, 2 ).finished();
            REQUIRE( references.image( 0 ).isApprox( expected ) );
            REQUIRE( targets.image( 0 ).isApprox( expected ) );
        }
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

/// Verify KLIPreduction::meanSubtract excludes masked pixels from RMS normalization.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP masked pixel RMS normalization", "[KLIPreduction][meanSubtract][rms][mask]" )
{
    reductionT reduction;
    reduction.m_meanSubMethod = mx::improc::HCI::meanSub::none;
    reduction.m_pixelTSNormMethod = mx::improc::HCI::pixelTSNorm::rms;

    mx::improc::eigenCube<float> references( 1, 2, 2 );
    references.image( 0 ) << 1, 99;
    references.image( 1 ) << 3, 88;
    mx::improc::eigenCube<float> targets( 1, 2, 1 );
    targets.image( 0 ) << 6, 77;
    reductionT::imageT mask( 1, 2 );
    mask << 1, 0;
    std::vector<float> norms;

    reduction.meanSubtract( references, targets, mask, norms );

    const float rms = std::sqrt( 5.0F );
    REQUIRE( references.image( 0 )( 0, 0 ) == Approx( 1.0F / rms ) );
    REQUIRE( references.image( 1 )( 0, 0 ) == Approx( 3.0F / rms ) );
    REQUIRE( targets.image( 0 )( 0, 0 ) == Approx( 6.0F / rms ) );
    REQUIRE( references.image( 0 )( 0, 1 ) == 0 );
    REQUIRE( references.image( 1 )( 0, 1 ) == 0 );
    REQUIRE( targets.image( 0 )( 0, 1 ) == 0 );
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

    included.setZero();
    mx::improc::collapseCovar<float>( selectedCovariance,
                                      covariance,
                                      norms,
                                      selectedReferences,
                                      references,
                                      references,
                                      1,
                                      0,
                                      0.45,
                                      mx::improc::HCI::exclude::none,
                                      mx::improc::HCI::exclude::angle,
                                      mx::improc::HCI::include::all,
                                      0,
                                      derotator,
                                      derotator,
                                      {},
                                      {},
                                      included );
    REQUIRE( included.row( 1 ).sum() == 3 );
    REQUIRE( included( 1, 3 ) == 0 );

    included.setZero();
    mx::improc::collapseCovar<float>( selectedCovariance,
                                      covariance,
                                      norms,
                                      selectedReferences,
                                      references,
                                      references,
                                      1,
                                      0,
                                      1,
                                      mx::improc::HCI::exclude::none,
                                      mx::improc::HCI::exclude::imno,
                                      mx::improc::HCI::include::all,
                                      0,
                                      derotator,
                                      derotator,
                                      {},
                                      {},
                                      included );
    REQUIRE( included.row( 1 ).sum() == 3 );
    REQUIRE( included( 1, 3 ) == 0 );

    included.setZero();
    mx::improc::collapseCovar<float>( selectedCovariance,
                                      covariance,
                                      norms,
                                      selectedReferences,
                                      references,
                                      references,
                                      1,
                                      10,
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
    REQUIRE( included.row( 1 ).isZero() );
    REQUIRE( selectedCovariance.size() == 0 );
    REQUIRE( selectedReferences.cols() == 0 );
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

/// Verify collapseCovar demotes zero-denominator and non-finite correlation candidates.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP covariance invalid correlations", "[KLIPreduction][covariance][include][nonfinite]" )
{
    Eigen::ArrayXXf covariance = Eigen::MatrixXf::Identity( 3, 3 ).array();
    Eigen::ArrayXXf references( 2, 3 );
    references << 1, std::numeric_limits<float>::infinity(), 0, 0, 0, 1;
    Eigen::ArrayXXf targets( 2, 1 );
    targets << 0, 1;
    const std::vector<float> norms{ 0, 1, 1 };
    testDerotator referenceDerotator{ { 0, 0, 0 } };
    testDerotator targetDerotator{ { 0 } };
    Eigen::ArrayXXf selectedCovariance;
    Eigen::ArrayXXf selectedReferences;
    mx::improc::eigenImage<int> included( 1, 3 );
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

    REQUIRE( included( 0, 0 ) == 0 );
    REQUIRE( included( 0, 1 ) == 0 );
    REQUIRE( included( 0, 2 ) == 1 );
    REQUIRE( selectedReferences.col( 0 ).isApprox( references.col( 2 ) ) );
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
                                                         static_cast<mx::improc::HCI::include>( 99 ),
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
                                                         { 0, std::numeric_limits<double>::quiet_NaN() },
                                                         { 0, 1 },
                                                         included ),
                       std::invalid_argument );
}

/// Verify calcKLModesAdaptive agrees with direct SVD and the legacy reference-Gram projector for every matrix shape.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP adaptive basis shape equivalence", "[KLIPreduction][basis][reference]" )
{
    reductionT::imageT target;
    reductionT::imageT covariance;
    reductionT::imageT modes;
    mx::math::syevrMem<double> workspace;

    SECTION( "pixels less than references" )
    {
        reductionT::imageT references( 3, 5 );
        references << 3, 0.2, 0, 1, -0.5, 0.1, 2, 0.4, -0.3, 0.8, 0, -0.2, 1, 0.7, 2;
        target.resize( 3, 1 );
        target << 5, -1, 2;

        double eigenSeconds = -1;
        double modeSeconds = -1;
        REQUIRE( mx::improc::calcKLModesAdaptive<double>( modes,
                                                          covariance,
                                                          references,
                                                          3,
                                                          &workspace,
                                                          &eigenSeconds,
                                                          &modeSeconds ) == 0 );
        REQUIRE( modes.rows() == 3 );
        REQUIRE( modes.cols() == 3 );
        REQUIRE( eigenSeconds >= 0 );
        REQUIRE( modeSeconds >= 0 );

        reductionT::imageT legacyCovariance;
        reductionT::imageT legacyModes;
        mx::math::eigenSYRK( legacyCovariance, references );
        REQUIRE( mx::math::calcKLModes<double>( legacyModes, legacyCovariance, references, 3 ) == 0 );

        for( int modeCount = 1; modeCount <= 3; ++modeCount )
        {
            const reductionT::imageT actualProjector = modeProjector( modes, modeCount );
            requireApprox( actualProjector, svdProjector( references, modeCount ) );
            requireApprox( actualProjector, modeProjector( legacyModes, modeCount ) );
            requireApprox( projectedResidual( actualProjector, target ),
                           projectedResidual( svdProjector( references, modeCount ), target ) );
        }
    }

    SECTION( "references less than pixels" )
    {
        reductionT::imageT references( 5, 3 );
        references << 3, 0.2, 1, 0.1, 2, -0.3, 0.7, -0.2, 1, -0.5, 0.8, 2, 0.4, -1.1, 0.6;
        target.resize( 5, 1 );
        target << 5, -1, 2, 0.5, -2;
        mx::math::eigenSYRK( covariance, references );

        REQUIRE( mx::improc::calcKLModesAdaptive<double>( modes, covariance, references, 3, &workspace ) == 0 );
        REQUIRE( modes.rows() == 3 );
        REQUIRE( modes.cols() == 5 );

        for( int modeCount = 1; modeCount <= 3; ++modeCount )
        {
            const reductionT::imageT actualProjector = modeProjector( modes, modeCount );
            requireApprox( actualProjector, svdProjector( references, modeCount ) );
            requireApprox( projectedResidual( actualProjector, target ),
                           projectedResidual( svdProjector( references, modeCount ), target ) );
        }
    }

    SECTION( "pixels equal references" )
    {
        reductionT::imageT references( 3, 3 );
        references << 3, 0.2, 1, 0.1, 2, -0.3, 0.7, -0.2, 1;
        target.resize( 3, 1 );
        target << 5, -1, 2;

        mx::math::eigenSYRK( covariance, references );
        REQUIRE( mx::improc::calcKLModesAdaptive<double>( modes, covariance, references, 3, &workspace ) == 0 );
        reductionT::imageT legacyModes;
        REQUIRE( mx::math::calcKLModes<double>( legacyModes, covariance, references, 3 ) == 0 );
        for( int modeCount = 1; modeCount <= 3; ++modeCount )
        {
            const reductionT::imageT actualProjector = modeProjector( modes, modeCount );
            requireApprox( actualProjector, svdProjector( references, modeCount ) );
            requireApprox( actualProjector, modeProjector( legacyModes, modeCount ) );
            requireApprox( projectedResidual( actualProjector, target ),
                           projectedResidual( svdProjector( references, modeCount ), target ) );
        }
    }
}

/// Verify calcKLModesAdaptive preserves legacy mode clamping and exact-positive rank behavior.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP adaptive basis mode and rank behavior", "[KLIPreduction][basis][rank]" )
{
    reductionT::imageT covariance;
    reductionT::imageT modes;
    mx::math::syevrMem<double> workspace;

    SECTION( "empty reference dimensions are rejected and clear timing outputs" )
    {
        double eigenSeconds = -1;
        double modeSeconds = -1;
        reductionT::imageT noPixels( 0, 3 );
        REQUIRE( mx::improc::calcKLModesAdaptive<double>( modes,
                                                          covariance,
                                                          noPixels,
                                                          1,
                                                          &workspace,
                                                          &eigenSeconds,
                                                          &modeSeconds ) == -1 );
        REQUIRE( eigenSeconds == 0 );
        REQUIRE( modeSeconds == 0 );

        eigenSeconds = -1;
        modeSeconds = -1;
        reductionT::imageT noReferences( 3, 0 );
        REQUIRE( mx::improc::calcKLModesAdaptive<double>( modes,
                                                          covariance,
                                                          noReferences,
                                                          1,
                                                          &workspace,
                                                          &eigenSeconds,
                                                          &modeSeconds ) == -1 );
        REQUIRE( eigenSeconds == 0 );
        REQUIRE( modeSeconds == 0 );
    }

    SECTION( "mode requests clamp to the available dimension" )
    {
        reductionT::imageT wideReferences( 3, 5 );
        wideReferences << 3, 0.2, 0, 1, -0.5, 0.1, 2, 0.4, -0.3, 0.8, 0, -0.2, 1, 0.7, 2;

        REQUIRE( mx::improc::calcKLModesAdaptive<double>( modes, covariance, wideReferences, 99, &workspace ) == 0 );
        REQUIRE( modes.rows() == 3 );
        REQUIRE( modes.cols() == 3 );
        requireApprox( modeProjector( modes, 3 ), svdProjector( wideReferences, 3 ) );

        REQUIRE( mx::improc::calcKLModesAdaptive<double>( modes, covariance, wideReferences, 2, &workspace ) == 0 );
        REQUIRE( modes.rows() == 2 );
        requireApprox( modeProjector( modes, 1 ), svdProjector( wideReferences, 1 ) );
        requireApprox( modeProjector( modes, 2 ), svdProjector( wideReferences, 2 ) );

        REQUIRE( mx::improc::calcKLModesAdaptive<double>( modes, covariance, wideReferences, 0, &workspace ) == 0 );
        REQUIRE( modes.rows() == 3 );
        REQUIRE( mx::improc::calcKLModesAdaptive<double>( modes, covariance, wideReferences, -7, &workspace ) == 0 );
        REQUIRE( modes.rows() == 3 );

        reductionT::imageT tallReferences( 5, 3 );
        tallReferences << 3, 0.2, 1, 0.1, 2, -0.3, 0.7, -0.2, 1, -0.5, 0.8, 2, 0.4, -1.1, 0.6;
        mx::math::eigenSYRK( covariance, tallReferences );
        REQUIRE( mx::improc::calcKLModesAdaptive<double>( modes, covariance, tallReferences, 99, &workspace ) == 0 );
        REQUIRE( modes.rows() == 3 );
        REQUIRE( modes.cols() == 5 );
    }

    SECTION( "wide exact rank deficiency" )
    {
        reductionT::imageT references = reductionT::imageT::Zero( 3, 5 );
        references( 0, 0 ) = 3;
        references( 1, 1 ) = 2;
        REQUIRE( mx::improc::calcKLModesAdaptive<double>( modes, covariance, references, 3, &workspace ) == 0 );
        REQUIRE( modes.row( 0 ).isZero() );
        requireApprox( modeProjector( modes, 2 ), modeProjector( modes, 3 ) );

        reductionT::imageT target( 3, 1 );
        target << 1, 2, 3;
        reductionT::imageT expected( 3, 1 );
        expected << 0, 0, 3;
        requireApprox( projectedResidual( modeProjector( modes, 2 ), target ), expected );
        requireApprox( projectedResidual( modeProjector( modes, 3 ), target ), expected );
    }

    SECTION( "tall exact rank deficiency" )
    {
        reductionT::imageT references = reductionT::imageT::Zero( 5, 3 );
        references( 0, 0 ) = 3;
        references( 1, 1 ) = 2;
        mx::math::eigenSYRK( covariance, references );
        REQUIRE( mx::improc::calcKLModesAdaptive<double>( modes, covariance, references, 3, &workspace ) == 0 );
        REQUIRE( modes.row( 0 ).isZero() );
        requireApprox( modeProjector( modes, 2 ), modeProjector( modes, 3 ) );
    }

    SECTION( "zero and positive near-zero spatial eigenvalues" )
    {
        reductionT::imageT references = reductionT::imageT::Zero( 3, 5 );
        REQUIRE( mx::improc::calcKLModesAdaptive<double>( modes, covariance, references, 3, &workspace ) == 0 );
        REQUIRE( modes.isZero() );

        references( 0, 0 ) = 3;
        references( 1, 1 ) = 2;
        references( 2, 2 ) = 1e-4F;
        REQUIRE( mx::improc::calcKLModesAdaptive<double>( modes, covariance, references, 3, &workspace ) == 0 );
        REQUIRE_FALSE( modes.row( 0 ).isZero() );
        reductionT::imageT identity( 3, 3 );
        identity.matrix().setIdentity();
        requireApprox( modeProjector( modes, 3 ), identity );
    }
}

/// Verify one calcKLModesAdaptive workspace can be reused across both Gram branches and changing dimensions.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP adaptive basis workspace reuse", "[KLIPreduction][basis][workspace]" )
{
    mx::math::syevrMem<double> workspace;
    reductionT::imageT covariance;
    reductionT::imageT modes;

    reductionT::imageT wideReferences( 3, 5 );
    wideReferences << 3, 0.2, 0, 1, -0.5, 0.1, 2, 0.4, -0.3, 0.8, 0, -0.2, 1, 0.7, 2;
    REQUIRE( mx::improc::calcKLModesAdaptive<double>( modes, covariance, wideReferences, 2, &workspace ) == 0 );
    requireApprox( modeProjector( modes, 2 ), svdProjector( wideReferences, 2 ) );

    reductionT::imageT tallReferences( 5, 2 );
    tallReferences << 3, 0.2, 0.1, 2, 0.7, -0.2, -0.5, 0.8, 1, -0.3;
    mx::math::eigenSYRK( covariance, tallReferences );
    REQUIRE( mx::improc::calcKLModesAdaptive<double>( modes, covariance, tallReferences, 2, &workspace ) == 0 );
    requireApprox( modeProjector( modes, 2 ), svdProjector( tallReferences, 2 ) );

    reductionT::imageT largerWideReferences( 4, 6 );
    largerWideReferences << 3, 0.2, 0, 1, -0.5, 0.4, 0.1, 2, 0.4, -0.3, 0.8, -0.2, 0, -0.2, 1, 0.7, 2, 0.3, 0.5, -0.1,
        0.2, 1.3, -0.7, 1.8;
    REQUIRE( mx::improc::calcKLModesAdaptive<double>( modes, covariance, largerWideReferences, 3, &workspace ) == 0 );
    requireApprox( modeProjector( modes, 3 ), svdProjector( largerWideReferences, 3 ) );
}

/// Verify KLIPreduction::worker uses the spatial Gram basis once for an unfiltered image library.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP worker adaptive master basis", "[KLIPreduction][worker][basis][OpenMP]" )
{
    TestDirectory directory;
    const std::vector<int> modeCounts{ 2, 1, 9, 1 };
    mx::improc::eigenCube<float> input( 2, 1, 4 );
    input.image( 0 ) << 3, 0.2;
    input.image( 1 ) << 0.1, 2;
    input.image( 2 ) << 1, -0.3;
    input.image( 3 ) << -0.5, 0.8;
    const reductionT::imageT referenceMatrix = input.cube();
    reductionT::imageT mask;
    std::vector<size_t> indices{ 0, 1 };

    reductionHarness serial;
    prepareWorkerReduction( serial, 2, 4, 4, modeCounts );
    mx::improc::eigenCube<float> serialInput = input;
    {
        OpenMPThreadGuard threads( 1 );
        serial.worker( serialInput, serialInput, mask, indices, 0, 0 );
    }

    reductionHarness parallel;
    prepareWorkerReduction( parallel, 2, 4, 4, modeCounts );
    mx::improc::eigenCube<float> parallelInput = input;
    {
        OpenMPThreadGuard threads( 3 );
        parallel.worker( parallelInput, parallelInput, mask, indices, 0, 0 );
    }

    REQUIRE( serial.m_Nmodes == modeCounts );
    REQUIRE( parallel.m_Nmodes == modeCounts );
    REQUIRE( serial.m_imsIncluded.isOnes() );
    REQUIRE( parallel.m_imsIncluded.isOnes() );
    for( size_t modeIndex = 0; modeIndex < modeCounts.size(); ++modeIndex )
    {
        const int effectiveModeCount = std::min( modeCounts[modeIndex], 2 );
        const reductionT::imageT projector = svdProjector( referenceMatrix, effectiveModeCount );
        requireApprox( serial.m_psfsub[modeIndex].cube(), parallel.m_psfsub[modeIndex].cube() );
        for( int targetIndex = 0; targetIndex < 4; ++targetIndex )
        {
            requireApprox( serial.m_psfsub[modeIndex].cube().col( targetIndex ),
                           projectedResidual( projector, referenceMatrix.col( targetIndex ) ) );
        }
    }
    requireApprox( serial.m_psfsub[1].cube(), serial.m_psfsub[3].cube() );

    for( const reductionHarness *result : { &serial, &parallel } )
    {
        REQUIRE( std::isfinite( result->t_eigenv ) );
        REQUIRE( std::isfinite( result->t_klim ) );
        REQUIRE( std::isfinite( result->t_psf ) );
        REQUIRE( result->t_eigenv >= 0 );
        REQUIRE( result->t_klim >= 0 );
        REQUIRE( result->t_psf >= 0 );
    }

    reductionHarness diagnostic;
    prepareWorkerReduction( diagnostic, 2, 4, 4, modeCounts );
    diagnostic.m_writeDiagnostics = true;
    diagnostic.m_diagnosticDirectory = directory.file( "adaptive-diagnostics" ).string();
    mx::improc::eigenCube<float> diagnosticInput = input;
    {
        OpenMPThreadGuard threads( 1 );
        diagnostic.worker( diagnosticInput, diagnosticInput, mask, indices, 0, 0 );
    }

    reductionT::imageT savedCovariance;
    mx::fits::fitsFile<float, mx::verbose::vv> reader;
    REQUIRE( reader.read( savedCovariance, directory.file( "adaptive-diagnostics/cv.fits" ).string() ) ==
             mx::error_t::noerror );
    REQUIRE( savedCovariance.rows() == 4 );
    REQUIRE( savedCovariance.cols() == 4 );
    const reductionT::imageT expectedCovariance =
        ( referenceMatrix.matrix().transpose() * referenceMatrix.matrix() ).array();
    for( int column = 0; column < 4; ++column )
    {
        for( int row = 0; row < 4; ++row )
        {
            REQUIRE( savedCovariance( row, column ) == Approx( expectedCovariance( row, column ) ).margin( 1e-5 ) );
        }
    }
}

/// Verify KLIPreduction::worker applies the spatial Gram basis after per-target RDI selection.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP worker adaptive selected RDI basis", "[KLIPreduction][worker][basis][RDI][include][OpenMP]" )
{
    const std::vector<int> modeCounts{ 1, 2, 7 };
    mx::improc::eigenCube<float> references( 2, 1, 4 );
    references.image( 0 ) << 3, 0.2;
    references.image( 1 ) << 0.1, 2;
    references.image( 2 ) << 1, -0.3;
    references.image( 3 ) << -0.5, 0.8;
    mx::improc::eigenCube<float> targets( 2, 1, 4 );
    targets.image( 0 ) << 0.8, -0.2;
    targets.image( 1 ) << -1, 1.5;
    targets.image( 2 ) << 2, 0.4;
    targets.image( 3 ) << 0.3, -0.7;
    const reductionT::imageT referenceMatrix = references.cube();
    const reductionT::imageT targetMatrix = targets.cube();
    reductionT::imageT mask;
    std::vector<size_t> indices{ 0, 1 };

    reductionHarness serial;
    prepareWorkerReduction( serial, 2, 4, 4, modeCounts );
    serial.m_includeMethod = mx::improc::HCI::include::imno;
    serial.m_includeRefNum = 3;
    serial.m_imsIncluded.setZero();
    mx::improc::eigenCube<float> serialReferences = references;
    mx::improc::eigenCube<float> serialTargets = targets;
    {
        OpenMPThreadGuard threads( 1 );
        serial.worker( serialReferences, serialTargets, mask, indices, 0, 0 );
    }

    reductionHarness parallel;
    prepareWorkerReduction( parallel, 2, 4, 4, modeCounts );
    parallel.m_includeMethod = mx::improc::HCI::include::imno;
    parallel.m_includeRefNum = 3;
    parallel.m_imsIncluded.setZero();
    mx::improc::eigenCube<float> parallelReferences = references;
    mx::improc::eigenCube<float> parallelTargets = targets;
    {
        OpenMPThreadGuard threads( 3 );
        parallel.worker( parallelReferences, parallelTargets, mask, indices, 0, 0 );
    }

    Eigen::ArrayXXi expectedIncluded( 4, 4 );
    expectedIncluded << 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1;
    REQUIRE( ( serial.m_imsIncluded == expectedIncluded ).all() );
    REQUIRE( ( parallel.m_imsIncluded == expectedIncluded ).all() );

    for( int targetIndex = 0; targetIndex < 4; ++targetIndex )
    {
        const std::array<int, 3> selectedIndices =
            targetIndex < 2 ? std::array<int, 3>{ 0, 1, 2 } : std::array<int, 3>{ 1, 2, 3 };
        reductionT::imageT selectedReferences( 2, 3 );
        for( int selectedIndex = 0; selectedIndex < 3; ++selectedIndex )
        {
            selectedReferences.col( selectedIndex ) = referenceMatrix.col( selectedIndices[selectedIndex] );
        }

        for( size_t modeIndex = 0; modeIndex < modeCounts.size(); ++modeIndex )
        {
            const int effectiveModeCount = std::min( modeCounts[modeIndex], 2 );
            const reductionT::imageT expected =
                projectedResidual( svdProjector( selectedReferences, effectiveModeCount ),
                                   targetMatrix.col( targetIndex ) );
            requireApprox( serial.m_psfsub[modeIndex].cube().col( targetIndex ), expected );
            requireApprox( parallel.m_psfsub[modeIndex].cube().col( targetIndex ), expected );
        }
    }

    for( size_t modeIndex = 0; modeIndex < modeCounts.size(); ++modeIndex )
    {
        requireApprox( serial.m_psfsub[modeIndex].cube(), parallel.m_psfsub[modeIndex].cube() );
    }
    for( const reductionHarness *result : { &serial, &parallel } )
    {
        REQUIRE( std::isfinite( result->t_eigenv ) );
        REQUIRE( std::isfinite( result->t_klim ) );
        REQUIRE( std::isfinite( result->t_psf ) );
        REQUIRE( result->t_eigenv >= 0 );
        REQUIRE( result->t_klim >= 0 );
        REQUIRE( result->t_psf >= 0 );
    }
}

/// Verify KLIPreduction::worker applies right-reason support masking to an adaptive spatial basis.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP worker adaptive right-reason basis", "[KLIPreduction][worker][basis][rightReason][OpenMP]" )
{
    mx::improc::eigenCube<float> input( 4, 1, 5 );
    input.image( 0 ) << 3, 0.2, 1, -0.5;
    input.image( 1 ) << 0.1, 2, -0.3, 0.8;
    input.image( 2 ) << 1, -0.3, 1.5, 0.4;
    input.image( 3 ) << -0.5, 0.8, 0.2, 2;
    input.image( 4 ) << 0.7, 1.1, -0.6, 0.3;
    const reductionT::imageT inputMatrix = input.cube();
    reductionT::imageT mask;
    std::vector<size_t> indices{ 0, 1, 2, 3 };

    reductionHarness serial;
    prepareWorkerReduction( serial, 4, 5, 5, { 2 } );
    serial.m_rightReason = true;
    serial.m_rightReasonRadius = 100;
    mx::improc::eigenCube<float> serialInput = input;
    {
        OpenMPThreadGuard threads( 1 );
        serial.worker( serialInput, serialInput, mask, indices, 0, 0 );
    }

    reductionHarness parallel;
    prepareWorkerReduction( parallel, 4, 5, 5, { 2 } );
    parallel.m_rightReason = true;
    parallel.m_rightReasonRadius = 100;
    mx::improc::eigenCube<float> parallelInput = input;
    {
        OpenMPThreadGuard threads( 3 );
        parallel.worker( parallelInput, parallelInput, mask, indices, 0, 0 );
    }

    requireApprox( serial.m_psfsub[0].cube(), inputMatrix );
    requireApprox( parallel.m_psfsub[0].cube(), inputMatrix );
    requireApprox( serial.m_psfsub[0].cube(), parallel.m_psfsub[0].cube() );
    REQUIRE( std::isfinite( serial.t_psf ) );
    REQUIRE( std::isfinite( parallel.t_psf ) );
    REQUIRE( serial.t_psf >= 0 );
    REQUIRE( parallel.t_psf >= 0 );
}

/// Verify KLIPreduction::worker safely rethrows a per-target empty-library failure after its OpenMP region.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP worker rejects an empty selected library", "[KLIPreduction][worker][exclude][error][OpenMP]" )
{
    reductionHarness reduction;
    prepareWorkerReduction( reduction, 2, 4, 4, { 1 } );
    reduction.m_excludeMethod = mx::improc::HCI::exclude::imno;

    mx::improc::eigenCube<float> images( 2, 1, 4 );
    images.image( 0 ) << 3, 0.2;
    images.image( 1 ) << 0.1, 2;
    images.image( 2 ) << 1, -0.3;
    images.image( 3 ) << -0.5, 0.8;
    reductionT::imageT mask;
    std::vector<size_t> indices{ 0, 1 };

    OpenMPThreadGuard threads( 3 );
    REQUIRE_THROWS_WITH( reduction.worker( images, images, mask, indices, 99, 0 ),
                         Catch::Matchers::Contains( "has no admissible references" ) );
}

/// Verify KLIPreduction::worker rejects empty or inconsistent geometry and propagates a master-solver failure.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP worker adaptive validation failures", "[KLIPreduction][worker][basis][validation][error]" )
{
    reductionT::imageT mask;
    std::vector<size_t> indices{ 0, 1 };

    SECTION( "empty image collection" )
    {
        reductionHarness reduction;
        prepareWorkerReduction( reduction, 2, 1, 2, { 1 } );
        mx::improc::eigenCube<float> empty;
        mx::improc::eigenCube<float> target( 2, 1, 1 );
        target.cube().setOnes();

        REQUIRE_THROWS_WITH( reduction.worker( empty, target, mask, indices, 0, 0 ),
                             Catch::Matchers::Contains( "KLIP requires reference and target images" ) );
    }

    SECTION( "inconsistent image geometry" )
    {
        reductionHarness reduction;
        prepareWorkerReduction( reduction, 2, 1, 2, { 1 } );
        mx::improc::eigenCube<float> references( 2, 1, 2 );
        references.cube().setOnes();
        mx::improc::eigenCube<float> target( 3, 1, 1 );
        target.cube().setOnes();

        REQUIRE_THROWS_WITH( reduction.worker( references, target, mask, indices, 0, 0 ),
                             Catch::Matchers::Contains( "invalid KLIP worker image dimensions" ) );
    }

    SECTION( "master eigensolver allocation failure" )
    {
        reductionHarness reduction;
        prepareWorkerReduction( reduction, 2, 4, 4, { 1 } );
        mx::improc::eigenCube<float> images( 2, 1, 4 );
        images.image( 0 ) << 3, 0.2;
        images.image( 1 ) << 0.1, 2;
        images.image( 2 ) << 1, -0.3;
        images.image( 3 ) << -0.5, 0.8;

        EigenAllocatorGuard allocatorFailure;
        REQUIRE_THROWS_WITH( reduction.worker( images, images, mask, indices, 0, 0 ),
                             Catch::Matchers::Contains( "KLIP eigensolver failed with status -1" ) );
    }
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

    reduction.m_rightReason = true;
    reduction.m_rightReasonRadius = 1;
    reduction.m_diagnosticDirectory = directory.file( "right-reason" ).string();
    reduction.m_psfsub[0].cube().setZero();
    images.image( 0 ) << 1, -1, 0;
    images.image( 1 ) << 1, 1, -2;
    reduction.worker( images, images, mask, indices, 0, 0 );

    REQUIRE( std::filesystem::exists( directory.file( "right-reason/rrMask.fits" ) ) );
    REQUIRE( std::filesystem::exists( directory.file( "right-reason/projMat.fits" ) ) );
    REQUIRE( std::filesystem::exists( directory.file( "right-reason/projMatrr.fits" ) ) );
    REQUIRE( reduction.m_psfsub[0].cube().isFinite().all() );
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

    reduction.m_RDIimageMJD = { 1, 2 };
    reduction.m_imageMJD = { std::numeric_limits<double>::quiet_NaN() };
    REQUIRE_THROWS( reduction.worker( references, targets, mask, indices, 0, 0 ) );

    reduction.m_includeMethod = static_cast<mx::improc::HCI::include>( 99 );
    reduction.m_includeRefNum = 0;
    REQUIRE_THROWS( reduction.worker( references, targets, mask, indices, 0, 0 ) );

    reduction.m_includeMethod = mx::improc::HCI::include::angle;
    reduction.m_includeRefNum = 1;
    reduction.m_RDIderotF.m_angles.clear();
    reduction.m_derotF.m_angles.clear();
    REQUIRE_THROWS( reduction.worker( references, targets, mask, indices, 0, 0 ) );

    reduction.m_RDIderotF.m_angles = { 0, std::numeric_limits<double>::quiet_NaN() };
    reduction.m_derotF.m_angles = { 0 };
    REQUIRE_THROWS( reduction.worker( references, targets, mask, indices, 0, 0 ) );

    reduction.m_RDIderotF.m_angles = { 0, 1 };
    reduction.m_derotF.m_angles = { std::numeric_limits<double>::quiet_NaN() };
    REQUIRE_THROWS( reduction.worker( references, targets, mask, indices, 0, 0 ) );

    reduction.m_includeMethod = mx::improc::HCI::include::all;
    reduction.m_includeRefNum = 0;
    reduction.m_meanSubMethod = static_cast<mx::improc::HCI::meanSub>( 99 );
    REQUIRE_THROWS( reduction.worker( references, targets, mask, indices, 0, 0 ) );
}

/// Verify KLIPreduction::regions converts each exclusion unit and restores configuration after worker failure.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP region exclusion conversion", "[KLIPreduction][regions][exclude][validation]" )
{
    const std::vector<mx::improc::HCI::exclude> methods{ mx::improc::HCI::exclude::pixel,
                                                         mx::improc::HCI::exclude::angle,
                                                         mx::improc::HCI::exclude::imno };

    for( const mx::improc::HCI::exclude method : methods )
    {
        reductionHarness reduction;
        prepareRegionReduction( reduction );
        reduction.m_Nmodes = { 1 };
        reduction.m_meanSubMethod = static_cast<mx::improc::HCI::meanSub>( 99 );
        reduction.m_excludeMethod = method;
        reduction.m_excludeMethodMax = method;
        reduction.m_minDPx = 0.1F;
        reduction.m_maxDPx = 2.0F;
        reduction.m_derotF.m_angles = { 0, 0.5, 1.0 };

        REQUIRE_THROWS( reduction.regions( 1, 2, 0, 360 ) );
        REQUIRE( reduction.m_excludeMethod == method );
        REQUIRE( reduction.m_excludeMethodMax == method );
    }
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

    prepareRegionReduction( reduction );
    reduction.m_padSize = -1;
    REQUIRE_THROWS( reduction.regions( 0, 2, 0, 360 ) );

    prepareRegionReduction( reduction );
    reduction.m_padSize = 4;
    reduction.m_refIms.resize( 2, 3, 1 );
    REQUIRE_THROWS( reduction.regions( 0, 2, 0, 360 ) );

    reductionHarness unloaded;
    unloaded.m_Nmodes = { 1 };
    unloaded.m_padSize = 0;
    REQUIRE_THROWS( unloaded.regions( 0, 2, 0, 360 ) );

    reductionHarness missingRDI;
    prepareRegionReduction( missingRDI );
    missingRDI.m_RDIfilesRead = false;
    missingRDI.m_RDIfileList = { "/definitely/not/a/reference.fits" };
    REQUIRE_THROWS( missingRDI.regions( 0, 2, 0, 360 ) );
}

/// Verify KLIPreduction::regions derives an image size and honors preprocess-only early completion.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP region early completion", "[KLIPreduction][regions][preprocess]" )
{
    reductionHarness inferredSize;
    prepareRegionReduction( inferredSize );
    inferredSize.m_Nmodes = { 1 };
    inferredSize.m_imSize = 0;
    inferredSize.m_padSize = 0;
    inferredSize.m_preProcess_only = true;
    inferredSize.m_skipPreProcess = false;

    REQUIRE( inferredSize.regions( 0, 2, 0, 360 ) == 0 );
    REQUIRE( inferredSize.m_imSize == 4 );
    REQUIRE( inferredSize.m_psfsub.empty() );
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

    REQUIRE( reduction.m_psfsubValidity.empty() );
    REQUIRE( reduction.finalProcess() == 0 );
    REQUIRE( reduction.m_psfsubValidity.empty() );

    const reductionT::imageT negative = reductionT::imageT::Constant( 2, 2, -2 );
    const reductionT::imageT positive = reductionT::imageT::Constant( 2, 2, 2 );
    REQUIRE( reduction.m_psfsub[0].image( 0 ).isApprox( negative ) );
    REQUIRE( reduction.m_psfsub[0].image( 1 ).isZero() );
    REQUIRE( reduction.m_psfsub[0].image( 2 ).isApprox( positive ) );
    REQUIRE( reduction.m_finim.rows() == 2 );
    REQUIRE( reduction.m_finim.cols() == 2 );
    REQUIRE( reduction.m_finim.planes() == 1 );
    REQUIRE( reduction.m_finim.image( 0 ).isZero() );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::KLIPreduction<float, mx::improc::ADIDerotator<float, mx::verbose::vv>, double, mx::verbose::vv>
        doxygenReduction;
    doxygenReduction.finalProcess();
#endif
    // clang-format on
}

/// Verify KLIPreduction::appendReductionHeader and finalProcess preserve complete, ordered KLIP FITS provenance.
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

    reductionHarness::fitsHeaderT reductionHeader;
    reduction.appendReductionHeader( reductionHeader );
    std::vector<std::string> reductionKeywords;
    for( auto card = reductionHeader.begin(); card != reductionHeader.end(); ++card )
    {
        if( !card->keyword().empty() )
        {
            reductionKeywords.push_back( card->keyword() );
        }
    }
    REQUIRE( reductionKeywords == std::vector<std::string>{ "MEAN SUB METHOD",
                                                            "PIXTS NORM METHOD",
                                                            "NMODES",
                                                            "RIGHT REASON",
                                                            "RIGHT REASON RADIUS",
                                                            "REGMINR",
                                                            "REGMAXR",
                                                            "REGMINQ",
                                                            "REGMAXQ",
                                                            "EXMTHDMN",
                                                            "MINDPX",
                                                            "EXMTHDMX",
                                                            "MAXDPX",
                                                            "INMTHDMX",
                                                            "INCLREFN" } );

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

    size_t adiPosition = header.size();
    size_t klipPosition = header.size();
    size_t position = 0;
    for( auto card = header.begin(); card != header.end(); ++card, ++position )
    {
        if( card->keyword() == "POSTMEDS" )
            adiPosition = position;
        else if( card->keyword() == "MEAN SUB METHOD" )
            klipPosition = position;
    }
    REQUIRE( adiPosition < klipPosition );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::KLIPreduction<float, mx::improc::ADIDerotator<float, mx::verbose::vv>, double, mx::verbose::vv>
        doxygenReduction;
    doxygenReduction.appendReductionHeader( reductionHeader );
#endif
    // clang-format on
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

/// Verify KLIPreduction::processPSFSub rejects malformed saved mode metadata.
/** \ingroup KLIPreduction_unit_tests */
TEST_CASE( "KLIP saved mode-list parse failure", "[KLIPreduction][processPSFSub][validation]" )
{
    TestDirectory directory;
    reductionT::imageT image( 2, 2 );
    image.setOnes();
    mx::fits::fitsHeader<mx::verbose::vv> header;
    header.append<int>( "REDUCTION", 0, "reduction index" );
    header.append<int>( "IMAGE", 0, "image index" );
    header.append<std::string>( "NMODES", "not-an-integer", "KLIP mode counts" );

    mx::fits::fitsFile<float, mx::verbose::vv> writer;
    REQUIRE( writer.write( directory.file( "malformed_000_00000.fits" ).string(), image, header ) ==
             mx::error_t::noerror );

    reductionHarness reduction;
    REQUIRE_THROWS( reduction.processPSFSub( directory.path().string(), "malformed" ) );
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

    mask.resize( 1, 0 );
    REQUIRE_THROWS( reduction.meanSubtract( references, targets, mask, norms ) );
    REQUIRE( references.image( 0 ).isApprox( original.image( 0 ) ) );

    mask.resize( 2, 2 );
    mask.setOnes();
    REQUIRE_THROWS( reduction.meanSubtract( references, targets, mask, norms ) );
    REQUIRE( references.image( 0 ).isApprox( original.image( 0 ) ) );
}

} // namespace KLIPreduction_test
} // namespace unitTest
