/** \file hciGaussianReference.cpp
 * \brief Exports the production hciAnalyze Gaussian filter before annular SNR normalization.
 */

#define HCIREDUCE_HCIANALYZE_NO_MAIN
#include "src/apps/hciAnalyze.cpp"

namespace
{

/// Expose the production Gaussian operation for a saved-image comparison.
class GaussianReference : public hciAnalyze
{
  public:
    using hciAnalyze::filterCube;
};

} // namespace

/// Read a science cube and export the same low-pass operation selected by filter.lpfGaussFW.
int main( int argc, /**< [in] argument count */
          char **argv /**< [in] input FITS, new output FITS, and positive FWHM in pixels */ )
{
    try
    {
        if( argc != 4 )
        {
            throw std::invalid_argument( "usage: hciGaussianReference INPUT.fits OUTPUT.fits FWHM" );
        }
        size_t parsed{ 0 };
        const float fwhm = std::stof( argv[3], &parsed );
        if( parsed != std::string( argv[3] ).size() || !std::isfinite( fwhm ) || fwhm <= 0 ||
            std::filesystem::exists( argv[2] ) )
        {
            throw std::invalid_argument( "require positive finite FWHM and a new output path" );
        }
        hciAnalyze::cubeT cube, invalid;
        hciAnalyze::fitsHeaderT header;
        mx::fits::fitsFile<float, mx::verbose::vv> file;
        if( file.read( cube, header, argv[1] ) != mx::error_t::noerror )
        {
            throw std::runtime_error( "cannot read input FITS" );
        }
        mx::improc::zeroNaNCube( cube, &invalid );
        GaussianReference application;
        application.filterCube( cube, invalid, 0, fwhm );
        for( int plane = 0; plane < cube.planes(); ++plane )
        {
            for( int column = 0; column < cube.cols(); ++column )
            {
                for( int row = 0; row < cube.rows(); ++row )
                {
                    if( invalid.image( plane )( row, column ) != 0 )
                    {
                        cube.image( plane )( row, column ) = std::numeric_limits<float>::quiet_NaN();
                    }
                }
            }
        }
        if( header.append( "LPFGFW", fwhm, "hciAnalyze Gaussian low-pass FWHM [pix]" ) != mx::error_t::noerror ||
            header.append( "GAUSREF", 1, "Gaussian filter before annular SNR normalization" ) != mx::error_t::noerror ||
            file.write( argv[2], cube, header ) != mx::error_t::noerror )
        {
            throw std::runtime_error( "cannot write Gaussian reference FITS" );
        }
        return EXIT_SUCCESS;
    }
    catch( const std::exception &error )
    {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
