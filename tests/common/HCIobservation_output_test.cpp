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
    observation.writeFinim( &additional );

    const auto firstPath = directory.file( "nested/final/combined_0000.fits" );
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
    REQUIRE( std::filesystem::exists( directory.file( "nested/final/combined_0001.fits" ) ) );

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
