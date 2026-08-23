/** \file HCIobservation_output_test.cpp
 * \brief Tests HCIobservation FITS headers and output products.
 * \author Jared R. Males
 */

#include "../catch2/catch.hpp"

#include "HCIobservation_test_fixture.hpp"

#include <fstream>

namespace unitTest
{
namespace HCIobservation_output_test
{

/// \cond HCIobservation_test_harness
namespace
{
HCIobservationTestHarness::imageT indexedImage( int rows, int cols, float offset = 0 )
{
    HCIobservationTestHarness::imageT image( rows, cols );
    for( int row = 0; row < rows; ++row )
    {
        for( int col = 0; col < cols; ++col )
        {
            image( row, col ) = offset + 10 * row + col;
        }
    }
    return image;
}

void readFitsImage( const std::filesystem::path &path,
                    HCIobservationTestHarness::imageT &image,
                    HCIobservationTestHarness::fitsHeaderT &header )
{
    HCIobservationTestHarness::fitsFileT file;
    const mx::error_t result = file.read( image, header, path.string() );
    if( result != mx::error_t::noerror )
    {
        throw mx::exception<mx::verbose::vv>( result, "reading output FITS image" );
    }
}

class currentPathGuard
{
  public:
    /// Change to a temporary current directory and remember the original location.
    explicit currentPathGuard( const std::filesystem::path &path /**< [in] temporary current directory */ )
        : m_previous( std::filesystem::current_path() )
    {
        std::filesystem::current_path( path );
    }

    /// Restore the original current directory.
    ~currentPathGuard()
    {
        std::filesystem::current_path( m_previous );
    }

  private:
    /// Current directory active before construction.
    std::filesystem::path m_previous;
};

} // namespace
/// \endcond

/// Verify HCIobservation::stdFitsHeader appends the configured reduction metadata without discarding input cards.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation standard FITS header", "[HCIobservation][output][stdFitsHeader]" )
{
    HCIobservationTestHarness observation;
    observation.m_deleteFront = 2;
    observation.m_deleteBack = 3;
    observation.m_qualityFile = "/tmp/quality-values.dat";
    observation.m_qualityThreshold = 0.75F;
    observation.m_Nims = 4;
    observation.m_imSize = 17;
    observation.m_coaddMethod = mx::improc::HCI::coadd::mean;
    observation.m_coaddMaxImno = 5;
    observation.m_coaddMaxTime = 2.5F;
    observation.m_coaddMaxAngle = 0.75F;
    observation.m_maskFile = "mask.fits";
    observation.m_preProcess_beforeCoadd = true;
    observation.m_preProcess_mask = false;
    observation.m_preProcess_subradprof = true;
    observation.m_preProcess_azUSM_azHalfWidth = 4;
    observation.m_preProcess_azUSM_maxAz = 35;
    observation.m_preProcess_azUSM_radHalfWidth = 2;
    observation.m_preProcess_medianUSM_fwhm = 3;
    observation.m_preProcess_gaussUSM_fwhm = 1.5F;
    observation.m_preProcess_meanSubMethod = mx::improc::HCI::meanSub::medianImage;
    observation.m_preProcess_pixelTSNormMethod = mx::improc::HCI::pixelTSNorm::rms;

    HCIobservationTestHarness::fitsHeaderT header;
    header.append<int>( "INPUT", 42, "original card" );
    observation.stdFitsHeader( header );

    REQUIRE( header["INPUT"].Int() == 42 );
    REQUIRE( header["FDELFRNT"].Int() == 2 );
    REQUIRE( header["FDELBACK"].Int() == 3 );
    REQUIRE( header["QFILE"].String() == "quality-values.dat" );
    REQUIRE( header["QTHRESH"].value<float>() == Approx( 0.75F ) );
    REQUIRE( header["NUMIMS"].Int() == 4 );
    REQUIRE( header["IMSIZE"].Int() == 17 );
    REQUIRE( header["COADMTHD"].String().starts_with( "mean" ) );
    REQUIRE( header["COADIMNO"].Int() == 5 );
    REQUIRE( header["COADTIME"].value<float>() == Approx( 2.5F ) );
    REQUIRE( header["COADANGL"].value<float>() == Approx( 0.75F ) );
    REQUIRE( header["MASKFILE"].String() == "mask.fits" );
    REQUIRE( header["PREPROC BEFORE"].Int() == 1 );
    REQUIRE( header["PREPROC MASK"].Int() == 0 );
    REQUIRE( header["PREPROC SUBRADPROF"].Int() == 1 );
    REQUIRE( header["PREPROC AZUSM AZWIDTH"].value<float>() == Approx( 4 ) );
    REQUIRE( header["PREPROC AZUSM MAXAZ"].value<float>() == Approx( 35 ) );
    REQUIRE( header["PREPROC AZUSM RADWIDTH"].value<float>() == Approx( 2 ) );
    REQUIRE( header["PREPROC MEDIANUSM FWHM"].value<float>() == Approx( 3 ) );
    REQUIRE( header["PREPROC GAUSSUSM FWHM"].value<float>() == Approx( 1.5F ) );
    REQUIRE( header["PREPROC MEANSUB METHOD"].String() == "medianImage" );
    REQUIRE( header["PREPROC PIXELTSNORM METHOD"].String() == "rms" );

    HCIobservationTestHarness noCoadd;
    HCIobservationTestHarness::fitsHeaderT noCoaddHeader;
    noCoadd.stdFitsHeader( noCoaddHeader );
    REQUIRE( noCoaddHeader["COADMTHD"].String() == "none" );
    REQUIRE( noCoaddHeader["COADIMNO"].Int() == 0 );
    REQUIRE( noCoaddHeader["COADTIME"].value<float>() == 0 );
    REQUIRE( noCoaddHeader["COADANGL"].value<float>() == 0 );
}

/// Verify HCIobservation::outputPreProcessed writes rectangular target images, copied headers, and nested paths.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation target preprocessed output", "[HCIobservation][output][preprocessed][target]" )
{
    TestDirectory directory;
    HCIobservationTestHarness observation;
    observation.m_preProcess_outputPrefix = directory.file( "nested/target_" ).string();
    observation.m_Nims = 2;
    observation.m_Nrows = 2;
    observation.m_Ncols = 3;
    observation.m_imSize = 2;
    observation.m_tgtIms.resize( 2, 3, 2 );
    observation.m_tgtIms.image( 0 ) = indexedImage( 2, 3 );
    observation.m_tgtIms.image( 1 ) = indexedImage( 2, 3, 100 );
    observation.m_heads.resize( 2 );
    observation.m_heads[0].append<int>( "ORIGINAL", 10, "original target card" );
    observation.m_heads[1].append<int>( "ORIGINAL", 20, "original target card" );

    observation.outputPreProcessed();

    for( int plane = 0; plane < 2; ++plane )
    {
        const auto path = directory.file( std::format( "nested/target_{:06}.fits", plane ) );
        REQUIRE( std::filesystem::exists( path ) );

        HCIobservationTestHarness::imageT image;
        HCIobservationTestHarness::fitsHeaderT header;
        readFitsImage( path, image, header );
        REQUIRE( image.rows() == 2 );
        REQUIRE( image.cols() == 3 );
        REQUIRE( image.isApprox( observation.m_tgtIms.image( plane ) ) );
        REQUIRE( header["ORIGINAL"].Int() == 10 * ( plane + 1 ) );
        REQUIRE( header["NUMIMS"].Int() == 2 );
    }

    observation.m_heads.pop_back();
    REQUIRE_THROWS( observation.outputPreProcessed() );
}

/// Verify HCIobservation::outputRDIPreProcessed writes distinct reference filenames, data, and copied headers.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation RDI preprocessed output", "[HCIobservation][output][preprocessed][RDI]" )
{
    TestDirectory directory;
    HCIobservationTestHarness observation;
    observation.m_preProcess_outputPrefix = directory.file( "nested/preprocessed_" ).string();
    observation.m_Nims = 1;
    observation.m_refIms.resize( 3, 2, 2 );
    observation.m_refIms.image( 0 ) = indexedImage( 3, 2, 50 );
    observation.m_refIms.image( 1 ) = indexedImage( 3, 2, 150 );
    observation.m_RDIheads.resize( 2 );
    observation.m_RDIheads[0].append<int>( "REFERENCE", 1, "original RDI card" );
    observation.m_RDIheads[1].append<int>( "REFERENCE", 2, "original RDI card" );

    observation.outputRDIPreProcessed();

    for( int plane = 0; plane < 2; ++plane )
    {
        const auto path = directory.file( std::format( "nested/preprocessed_RDI_{:06}.fits", plane ) );
        REQUIRE( std::filesystem::exists( path ) );

        HCIobservationTestHarness::imageT image;
        HCIobservationTestHarness::fitsHeaderT header;
        readFitsImage( path, image, header );
        REQUIRE( image.rows() == 3 );
        REQUIRE( image.cols() == 2 );
        REQUIRE( image.isApprox( observation.m_refIms.image( plane ) ) );
        REQUIRE( header["REFERENCE"].Int() == plane + 1 );
    }

    observation.m_RDIheads.clear();
    REQUIRE_THROWS( observation.outputRDIPreProcessed() );
}

/// Verify HCIobservation preprocessed-output methods reject empty inputs and propagate directory failures.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation preprocessed output validation", "[HCIobservation][output][preprocessed][validation]" )
{
    TestDirectory directory;
    HCIobservationTestHarness noPrefix;
    REQUIRE_NOTHROW( noPrefix.outputPreProcessed() );
    REQUIRE_NOTHROW( noPrefix.outputRDIPreProcessed() );

    HCIobservationTestHarness emptyTarget;
    emptyTarget.m_preProcess_outputPrefix = directory.file( "empty-target_" ).string();
    REQUIRE_THROWS( emptyTarget.outputPreProcessed() );

    HCIobservationTestHarness emptyRDI;
    emptyRDI.m_preProcess_outputPrefix = directory.file( "empty-reference_" ).string();
    REQUIRE_THROWS( emptyRDI.outputRDIPreProcessed() );

    const auto blocker = directory.file( "not-a-directory" );
    writeTextFile( blocker, "regular file" );
    HCIobservationTestHarness blocked;
    blocked.m_preProcess_outputPrefix = ( blocker / "target_" ).string();
    blocked.m_tgtIms.resize( 1, 1, 1 );
    blocked.m_tgtIms.image( 0 ).setZero();
    blocked.m_heads.resize( 1 );
    REQUIRE_THROWS( blocked.outputPreProcessed() );

    HCIobservationTestHarness targetWriteFailure;
    targetWriteFailure.m_preProcess_outputPrefix = directory.file( "target-write_" ).string();
    targetWriteFailure.m_tgtIms.resize( 1, 1, 1 );
    targetWriteFailure.m_tgtIms.image( 0 ).setZero();
    targetWriteFailure.m_heads.resize( 1 );
    std::filesystem::create_directory( directory.file( "target-write_000000.fits" ) );
    writeTextFile( directory.file( "target-write_000000.fits/marker" ), "keep directory non-empty" );
    REQUIRE_THROWS( targetWriteFailure.outputPreProcessed() );

    HCIobservationTestHarness blockedRDI;
    blockedRDI.m_preProcess_outputPrefix = ( blocker / "reference_" ).string();
    blockedRDI.m_refIms.resize( 1, 1, 1 );
    blockedRDI.m_refIms.image( 0 ).setZero();
    blockedRDI.m_RDIheads.resize( 1 );
    REQUIRE_THROWS( blockedRDI.outputRDIPreProcessed() );

    HCIobservationTestHarness rdiWriteFailure;
    rdiWriteFailure.m_preProcess_outputPrefix = directory.file( "reference-write_" ).string();
    rdiWriteFailure.m_refIms.resize( 1, 1, 1 );
    rdiWriteFailure.m_refIms.image( 0 ).setZero();
    rdiWriteFailure.m_RDIheads.resize( 1 );
    std::filesystem::create_directory( directory.file( "reference-write_RDI_000000.fits" ) );
    writeTextFile( directory.file( "reference-write_RDI_000000.fits/marker" ), "keep directory non-empty" );
    REQUIRE_THROWS( rdiWriteFailure.outputRDIPreProcessed() );
}

/// Verify HCIobservation::writeFinim writes sequential and exact cube products with complete headers.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation final-image output", "[HCIobservation][output][writeFinim]" )
{
    TestDirectory directory;
    HCIobservationTestHarness observation;
    observation.m_outputDir = directory.file( "nested/final" ).string();
    observation.m_finimName = "combined_";
    observation.m_Nims = 3;
    observation.m_imSize = 2;
    observation.m_combineMethod = mx::improc::HCI::combine::mean;
    observation.m_weightFile = "weights.dat";
    observation.m_minGoodFract = 0.5F;
    observation.m_finim.resize( 2, 3, 2 );
    observation.m_finim.image( 0 ) = indexedImage( 2, 3, 20 );
    observation.m_finim.image( 1 ) = indexedImage( 2, 3, 120 );

    HCIobservationTestHarness::fitsHeaderT additional;
    additional.append<int>( "EXTRA", 77, "additional output card" );
    const auto occupiedAuxiliaryPath = directory.file( "nested/final/combined_0045_outputs" );
    std::filesystem::create_directories( occupiedAuxiliaryPath );
    writeTextFile( occupiedAuxiliaryPath / "optimizer-summary.yaml", "reserved by an auxiliary-only run" );
    writeTextFile( directory.file( "nested/final/combined_0046.fits" ), "later retained product" );
    const auto firstPath = directory.file( "nested/final/combined_0047.fits" );
    REQUIRE( observation.finalImageOutputPath() == firstPath.string() );
    REQUIRE( observation.finalImageOutputPath() == firstPath.string() );
    HCIobservationTestHarness::fitsHeaderT constructedHeader;
    observation.finalImageHeader( constructedHeader, &additional );
    REQUIRE( constructedHeader["COMBINATION METHOD"].String().starts_with( "mean" ) );
    REQUIRE( constructedHeader["EXTRA"].Int() == 77 );
    observation.writeFinim( &additional );

    REQUIRE( std::filesystem::exists( firstPath ) );
    HCIobservationTestHarness::cubeT cube;
    HCIobservationTestHarness::fitsHeaderT header;
    HCIobservationTestHarness::fitsFileT file;
    REQUIRE( file.read( cube, header, firstPath.string() ) == mx::error_t::noerror );
    REQUIRE( cube.rows() == 2 );
    REQUIRE( cube.cols() == 3 );
    REQUIRE( cube.planes() == 2 );
    REQUIRE( cube.image( 0 ).isApprox( observation.m_finim.image( 0 ) ) );
    REQUIRE( cube.image( 1 ).isApprox( observation.m_finim.image( 1 ) ) );
    REQUIRE( header["COMBINATION METHOD"].String().starts_with( "mean" ) );
    REQUIRE( header["WEIGHT FILE"].String() == "weights.dat" );
    REQUIRE( header["MIN GOOD FRACTION"].value<float>() == Approx( 0.5F ) );
    REQUIRE( header["EXTRA"].Int() == 77 );

    observation.writeFinim();
    REQUIRE( std::filesystem::exists( directory.file( "nested/final/combined_0048.fits" ) ) );

    const auto resolvedPath = directory.file( "nested/final/resolved-product.fits" );
    observation.writeFinimAtPath( resolvedPath.string(), &additional );
    REQUIRE( std::filesystem::exists( resolvedPath ) );

    observation.m_exactFinimName = true;
    observation.m_finimName = "exact-product.bin";
    observation.m_finim.image( 0 ).setConstant( 31 );
    observation.writeFinim();
    observation.m_finim.image( 0 ).setConstant( 47 );
    observation.writeFinim();

    const auto exactPath = directory.file( "nested/final/exact-product.bin" );
    REQUIRE( std::filesystem::exists( exactPath ) );
    REQUIRE( file.read( cube, exactPath.string() ) == mx::error_t::noerror );
    REQUIRE( cube.image( 0 ).minCoeff() == 47 );
    REQUIRE( cube.image( 0 ).maxCoeff() == 47 );

    observation.m_combineMethod = mx::improc::HCI::combine::sigmaMean;
    observation.m_sigmaThreshold = 2.25F;
    observation.m_finimName = "sigma-product.fits";
    observation.writeFinim();
    HCIobservationTestHarness::fitsHeaderT sigmaHeader;
    REQUIRE( file.read( cube, sigmaHeader, directory.file( "nested/final/sigma-product.fits" ).string() ) ==
             mx::error_t::noerror );
    REQUIRE( sigmaHeader["SIGMA THRESHOLD"].Double() == Approx( 2.25 ) );

    HCIobservationTestHarness empty;
    empty.m_outputDir = directory.file( "unused" ).string();
    REQUIRE_THROWS( empty.writeFinim() );

    const auto blocker = directory.file( "final-blocker" );
    writeTextFile( blocker, "regular file" );
    HCIobservationTestHarness blocked;
    blocked.m_outputDir = ( blocker / "output" ).string();
    blocked.m_finim.resize( 1, 1, 1 );
    blocked.m_finim.image( 0 ).setZero();
    REQUIRE_THROWS( blocked.writeFinim() );

    HCIobservationTestHarness writeFailure;
    writeFailure.m_outputDir = directory.path().string();
    writeFailure.m_finimName = "existing-directory";
    writeFailure.m_exactFinimName = true;
    writeFailure.m_finim.resize( 1, 1, 1 );
    writeFailure.m_finim.image( 0 ).setZero();
    std::filesystem::create_directory( directory.file( "existing-directory" ) );
    writeTextFile( directory.file( "existing-directory/marker" ), "keep directory non-empty" );
    REQUIRE_THROWS( writeFailure.writeFinim() );
}

/// Verify HCIobservation::outputPSFSub writes every reduction/image, copied headers, indices, and weight sidecar.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation PSF-subtracted output", "[HCIobservation][output][outputPSFSub]" )
{
    TestDirectory directory;
    HCIobservationTestHarness observation;
    observation.m_outputDir = directory.file( "nested/psf-output" ).string();
    observation.m_PSFSubPrefix = "products/psf";
    observation.m_Nims = 2;
    observation.m_imSize = 2;
    observation.m_comboWeights = { 0.25F, 0.75F };
    observation.m_heads.resize( 2 );
    observation.m_heads[0].append<int>( "SOURCE", 10, "source image" );
    observation.m_heads[1].append<int>( "SOURCE", 20, "source image" );
    observation.m_psfsub.resize( 2 );
    for( size_t reduction = 0; reduction < observation.m_psfsub.size(); ++reduction )
    {
        observation.m_psfsub[reduction].resize( 2, 3, 2 );
        for( int plane = 0; plane < 2; ++plane )
        {
            observation.m_psfsub[reduction].image( plane ) =
                indexedImage( 2, 3, 100 * static_cast<float>( reduction ) + 10 * plane );
        }
    }

    HCIobservationTestHarness::fitsHeaderT additional;
    additional.append<int>( "EXTRA", 88, "additional output card" );
    observation.outputPSFSub( &additional );

    HCIobservationTestHarness::fitsFileT file;
    for( size_t reduction = 0; reduction < observation.m_psfsub.size(); ++reduction )
    {
        for( int plane = 0; plane < 2; ++plane )
        {
            const auto path =
                directory.file( std::format( "nested/psf-output/products/psf_{:03}_{:05}.fits", reduction, plane ) );
            REQUIRE( std::filesystem::exists( path ) );

            HCIobservationTestHarness::imageT image;
            HCIobservationTestHarness::fitsHeaderT header;
            REQUIRE( file.read( image, header, path.string() ) == mx::error_t::noerror );
            REQUIRE( image.isApprox( observation.m_psfsub[reduction].image( plane ) ) );
            REQUIRE( header["SOURCE"].Int() == 10 * ( plane + 1 ) );
            REQUIRE( header["EXTRA"].Int() == 88 );
            REQUIRE( header["REDUCTION"].Int() == static_cast<int>( reduction ) );
            REQUIRE( header["IMAGE"].Int() == plane );
        }
    }

    const auto weightPath = directory.file( "nested/psf-output/products/psfweights.dat" );
    REQUIRE( std::filesystem::exists( weightPath ) );
    std::ifstream weights( weightPath );
    std::string firstLine;
    std::string secondLine;
    REQUIRE( static_cast<bool>( std::getline( weights, firstLine ) ) );
    REQUIRE( static_cast<bool>( std::getline( weights, secondLine ) ) );
    REQUIRE( firstLine.find( "psf_000_00000.fits 0.25" ) != std::string::npos );
    REQUIRE( secondLine.find( "psf_000_00001.fits 0.75" ) != std::string::npos );

    HCIobservationTestHarness loaded;
    loaded.readPSFSub( directory.file( "nested/psf-output/products" ).string(), "psf" );
    REQUIRE( loaded.m_filesRead );
    REQUIRE( loaded.m_psfsub.size() == observation.m_psfsub.size() );
    REQUIRE( loaded.m_Nrows == 2 );
    REQUIRE( loaded.m_Ncols == 3 );
    REQUIRE( loaded.m_Nims == 2 );
    REQUIRE( loaded.m_heads.size() == 2 );
    REQUIRE( loaded.m_hookEvents == std::vector<std::string>{ "target-read" } );
    for( size_t reduction = 0; reduction < loaded.m_psfsub.size(); ++reduction )
    {
        REQUIRE( loaded.m_psfsub[reduction].cube().isApprox( observation.m_psfsub[reduction].cube() ) );
    }

    observation.m_comboWeights.pop_back();
    REQUIRE_THROWS( observation.outputPSFSub() );

    observation.m_comboWeights.clear();
    observation.m_heads.clear();
    REQUIRE_THROWS( observation.outputPSFSub() );

    observation.m_heads.resize( 2 );
    observation.m_psfsub[1].resize( 1, 1, 1 );
    REQUIRE_THROWS( observation.outputPSFSub() );

    HCIobservationTestHarness empty;
    REQUIRE_THROWS( empty.outputPSFSub() );

    HCIobservationTestHarness direct;
    direct.m_PSFSubPrefix = directory.file( "direct-psf" ).string();
    direct.m_psfsub.resize( 1 );
    direct.m_psfsub[0].resize( 1, 2, 1 );
    direct.m_psfsub[0].image( 0 ).setConstant( 9 );
    direct.m_heads.resize( 1 );
    direct.outputPSFSub();
    REQUIRE( std::filesystem::exists( directory.file( "direct-psf_000_00000.fits" ) ) );

    {
        currentPathGuard currentDirectory( directory.path() );
        HCIobservationTestHarness relative;
        relative.m_PSFSubPrefix = "relative-psf";
        relative.m_psfsub.resize( 1 );
        relative.m_psfsub[0].resize( 1, 1, 1 );
        relative.m_psfsub[0].image( 0 ).setZero();
        relative.m_heads.resize( 1 );
        relative.outputPSFSub();
    }
    REQUIRE( std::filesystem::exists( directory.file( "relative-psf_000_00000.fits" ) ) );

    HCIobservationTestHarness emptyCube;
    emptyCube.m_psfsub.resize( 1 );
    REQUIRE_THROWS( emptyCube.outputPSFSub() );

    const auto blocker = directory.file( "psf-blocker" );
    writeTextFile( blocker, "regular file" );
    HCIobservationTestHarness blocked;
    blocked.m_outputDir = ( blocker / "output" ).string();
    blocked.m_PSFSubPrefix = "psf";
    blocked.m_psfsub.resize( 1 );
    blocked.m_psfsub[0].resize( 1, 1, 1 );
    blocked.m_psfsub[0].image( 0 ).setZero();
    blocked.m_heads.resize( 1 );
    REQUIRE_THROWS( blocked.outputPSFSub() );

    HCIobservationTestHarness weightOpenFailure;
    weightOpenFailure.m_outputDir = directory.path().string();
    weightOpenFailure.m_PSFSubPrefix = "blocked-weights/psf";
    weightOpenFailure.m_psfsub.resize( 1 );
    weightOpenFailure.m_psfsub[0].resize( 1, 1, 1 );
    weightOpenFailure.m_psfsub[0].image( 0 ).setZero();
    weightOpenFailure.m_heads.resize( 1 );
    weightOpenFailure.m_comboWeights = { 1 };
    std::filesystem::create_directories( directory.file( "blocked-weights/psfweights.dat" ) );
    REQUIRE_THROWS( weightOpenFailure.outputPSFSub() );

    HCIobservationTestHarness weightWriteFailure;
    weightWriteFailure.m_PSFSubPrefix = directory.file( "full" ).string();
    weightWriteFailure.m_psfsub.resize( 1 );
    weightWriteFailure.m_psfsub[0].resize( 1, 1, 1 );
    weightWriteFailure.m_psfsub[0].image( 0 ).setZero();
    weightWriteFailure.m_heads.resize( 1 );
    weightWriteFailure.m_comboWeights = { 1 };
    std::filesystem::create_symlink( "/dev/full", directory.file( "fullweights.dat" ) );
    REQUIRE_THROWS( weightWriteFailure.outputPSFSub() );

    HCIobservationTestHarness imageWriteFailure;
    imageWriteFailure.m_outputDir = directory.path().string();
    imageWriteFailure.m_PSFSubPrefix = "blocked-image/psf";
    imageWriteFailure.m_psfsub.resize( 1 );
    imageWriteFailure.m_psfsub[0].resize( 1, 1, 1 );
    imageWriteFailure.m_psfsub[0].image( 0 ).setZero();
    imageWriteFailure.m_heads.resize( 1 );
    std::filesystem::create_directories( directory.file( "blocked-image/psf_000_00000.fits" ) );
    writeTextFile( directory.file( "blocked-image/psf_000_00000.fits/marker" ), "keep directory non-empty" );
    REQUIRE_THROWS( imageWriteFailure.outputPSFSub() );
}

/// Verify HCIobservation::readPSFSub rejects incomplete grids, malformed names, and mismatched header indices.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation saved PSF-subtracted validation", "[HCIobservation][readPSFSub][validation]" )
{
    TestDirectory directory;
    const auto makeSavedImage = [&directory]( const std::string &name, int reduction, int imageIndex )
    {
        HCIobservationTestHarness::imageT image = indexedImage( 2, 3, 10 * reduction + imageIndex );
        HCIobservationTestHarness::fitsHeaderT header;
        header.append<int>( "REDUCTION", reduction, "reduction index" );
        header.append<int>( "IMAGE", imageIndex, "image index" );
        writeFitsImage( directory.file( name ), image, &header );
    };

    makeSavedImage( "saved_000_00000.fits", 0, 0 );
    makeSavedImage( "saved_000_00001.fits", 0, 1 );
    makeSavedImage( "saved_001_00000.fits", 1, 0 );

    HCIobservationTestHarness incomplete;
    REQUIRE_THROWS( incomplete.readPSFSub( directory.path().string(), "saved" ) );
    REQUIRE( incomplete.m_psfsub.empty() );

    makeSavedImage( "saved_001_00001.fits", 9, 1 );
    HCIobservationTestHarness mismatchedHeader;
    REQUIRE_THROWS( mismatchedHeader.readPSFSub( directory.path().string(), "saved" ) );
    REQUIRE( mismatchedHeader.m_psfsub.empty() );

    writeTextFile( directory.file( "saved_bad_00002.fits" ), "not a FITS file" );
    HCIobservationTestHarness malformedName;
    REQUIRE_THROWS( malformedName.readPSFSub( directory.path().string(), "saved" ) );

    HCIobservationTestHarness missing;
    REQUIRE_THROWS( missing.readPSFSub( directory.file( "missing" ).string(), "saved" ) );
    REQUIRE_THROWS( missing.readPSFSub( directory.path().string(), "" ) );

    const auto emptyDirectory = directory.file( "empty" );
    std::filesystem::create_directory( emptyDirectory );
    REQUIRE_THROWS( missing.readPSFSub( emptyDirectory.string(), "saved" ) );

    const auto unreadableDirectory = directory.file( "unreadable" );
    std::filesystem::create_directory( unreadableDirectory );
    std::filesystem::permissions( unreadableDirectory, std::filesystem::perms::none );
    REQUIRE_THROWS( missing.readPSFSub( unreadableDirectory.string(), "saved" ) );
    std::filesystem::permissions( unreadableDirectory, std::filesystem::perms::owner_all );

    const auto danglingDirectory = directory.file( "dangling" );
    std::filesystem::create_directory( danglingDirectory );
    std::filesystem::create_symlink( danglingDirectory / "missing-target", danglingDirectory / "saved_0_0.fits" );
    REQUIRE_THROWS( missing.readPSFSub( danglingDirectory.string(), "saved" ) );

    const auto malformedDirectory = directory.file( "malformed" );
    std::filesystem::create_directory( malformedDirectory );
    writeTextFile( malformedDirectory / "saved_bad.fits", "not a FITS file" );
    REQUIRE_THROWS( missing.readPSFSub( malformedDirectory.string(), "saved" ) );

    const auto duplicateDirectory = directory.file( "duplicate" );
    std::filesystem::create_directory( duplicateDirectory );
    HCIobservationTestHarness::imageT duplicateImage = indexedImage( 1, 1 );
    HCIobservationTestHarness::fitsHeaderT duplicateHeader;
    duplicateHeader.append<int>( "REDUCTION", 0, "reduction index" );
    duplicateHeader.append<int>( "IMAGE", 0, "image index" );
    writeFitsImage( duplicateDirectory / "saved_0_0.fits", duplicateImage, &duplicateHeader );
    writeFitsImage( duplicateDirectory / "saved_00_00.fits", duplicateImage, &duplicateHeader );
    REQUIRE_THROWS( missing.readPSFSub( duplicateDirectory.string(), "saved" ) );

    const auto invalidFitsDirectory = directory.file( "invalid-fits" );
    std::filesystem::create_directory( invalidFitsDirectory );
    writeTextFile( invalidFitsDirectory / "saved_0_0.fits", "not a FITS file" );
    REQUIRE_THROWS( missing.readPSFSub( invalidFitsDirectory.string(), "saved" ) );

    const auto inconsistentDirectory = directory.file( "inconsistent" );
    std::filesystem::create_directory( inconsistentDirectory );
    HCIobservationTestHarness::imageT firstImage = indexedImage( 1, 1 );
    HCIobservationTestHarness::imageT secondImage = indexedImage( 2, 1 );
    HCIobservationTestHarness::fitsHeaderT firstHeader;
    firstHeader.append<int>( "REDUCTION", 0, "reduction index" );
    firstHeader.append<int>( "IMAGE", 0, "image index" );
    HCIobservationTestHarness::fitsHeaderT secondHeader;
    secondHeader.append<int>( "REDUCTION", 0, "reduction index" );
    secondHeader.append<int>( "IMAGE", 1, "image index" );
    writeFitsImage( inconsistentDirectory / "saved_0_0.fits", firstImage, &firstHeader );
    writeFitsImage( inconsistentDirectory / "saved_0_1.fits", secondImage, &secondHeader );
    REQUIRE_THROWS( missing.readPSFSub( inconsistentDirectory.string(), "saved" ) );

    // clang-format off
#ifdef __DOXY_ONLY__
    mx::improc::HCIobservation<float, mx::verbose::vv>::readPSFSub( "", "" );
#endif
    // clang-format on
}

/// Verify the HCIobservation read-to-write workflow preserves independently calculated data and essential metadata.
/** \ingroup HCIobservation_unit_tests */
TEST_CASE( "HCIobservation golden read-to-write workflow", "[HCIobservation][output][workflow]" )
{
    TestDirectory directory;
    HCIobservationTestHarness observation;
    observation.m_dateKeyword.clear();
    observation.m_coaddMethod = mx::improc::HCI::coadd::mean;
    observation.m_coaddMaxImno = 2;
    observation.m_preProcess_meanSubMethod = mx::improc::HCI::meanSub::meanImage;
    observation.m_combineMethod = mx::improc::HCI::combine::median;

    const std::vector<float> values{ 1, 1, 3, 3, 8, 8 };
    for( size_t index = 0; index < values.size(); ++index )
    {
        HCIobservationTestHarness::imageT image( 2, 3 );
        image.setConstant( values[index] );
        const auto path = directory.file( std::format( "input_{:02}.fits", index ) );
        writeFitsImage( path, image );
        observation.m_fileList.push_back( path.string() );
    }

    observation.readFiles();
    REQUIRE( observation.m_tgtIms.planes() == 3 );
    REQUIRE( observation.m_tgtIms.image( 0 ).minCoeff() == Approx( -3 ) );
    REQUIRE( observation.m_tgtIms.image( 1 ).minCoeff() == Approx( -1 ) );
    REQUIRE( observation.m_tgtIms.image( 2 ).minCoeff() == Approx( 4 ) );

    observation.m_psfsub = { observation.m_tgtIms };
    observation.combineFinim();
    REQUIRE( observation.m_finim.planes() == 1 );
    REQUIRE( observation.m_finim.image( 0 ).minCoeff() == Approx( -1 ) );
    REQUIRE( observation.m_finim.image( 0 ).maxCoeff() == Approx( -1 ) );

    observation.m_outputDir = directory.file( "nested/workflow" ).string();
    observation.m_finimName = "golden.fits";
    observation.m_exactFinimName = true;
    observation.writeFinim();

    HCIobservationTestHarness::cubeT output;
    HCIobservationTestHarness::fitsHeaderT header;
    HCIobservationTestHarness::fitsFileT file;
    const auto outputPath = directory.file( "nested/workflow/golden.fits" );
    REQUIRE( file.read( output, header, outputPath.string() ) == mx::error_t::noerror );
    REQUIRE( output.rows() == 2 );
    REQUIRE( output.cols() == 2 );
    REQUIRE( output.planes() == 1 );
    REQUIRE( output.image( 0 ).minCoeff() == Approx( -1 ) );
    REQUIRE( output.image( 0 ).maxCoeff() == Approx( -1 ) );
    REQUIRE( header["NUMIMS"].Int() == 3 );
    REQUIRE( header["COADMTHD"].String().starts_with( "mean" ) );
    REQUIRE( header["COADIMNO"].Int() == 2 );
    REQUIRE( header["PREPROC MEANSUB METHOD"].String().starts_with( "meanImage" ) );
    REQUIRE( header["COMBINATION METHOD"].String().starts_with( "median" ) );
}

} // namespace HCIobservation_output_test
} // namespace unitTest
