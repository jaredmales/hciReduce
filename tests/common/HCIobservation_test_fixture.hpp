/** \file HCIobservation_test_fixture.hpp
 * \brief Shared test harness and temporary-input helpers for HCIobservation tests.
 * \author Jared R. Males
 */

#ifndef tests_common_HCIobservation_test_fixture_hpp
#define tests_common_HCIobservation_test_fixture_hpp

#include <filesystem>
#include <string>
#include <vector>

#include "src/common/HCIobservation.hpp"

/// \cond HCIobservation_test_harness
namespace unitTest
{

/// Protected-access harness for focused HCIobservation tests.
struct HCIobservationTestHarness : public mx::improc::HCIobservation<float, mx::verbose::vv>
{
    using baseT = mx::improc::HCIobservation<float, mx::verbose::vv>;
    using cubeT = mx::improc::eigenCube<float>;

    using baseT::load_fileList;
    using baseT::m_coaddKeywords;
    using baseT::m_coaddMaxImno;
    using baseT::m_coaddMaxTime;
    using baseT::m_coaddMethod;
    using baseT::m_combineMethod;
    using baseT::m_comboWeights;
    using baseT::m_dateIsISO8601;
    using baseT::m_dateKeyword;
    using baseT::m_dateUnit;
    using baseT::m_deleteBack;
    using baseT::m_deleteFront;
    using baseT::m_directory;
    using baseT::m_doOutputPSFSub;
    using baseT::m_exactFinimName;
    using baseT::m_extension;
    using baseT::m_fileList;
    using baseT::m_fileListFile;
    using baseT::m_filesDeleted;
    using baseT::m_filesRead;
    using baseT::m_finimName;
    using baseT::m_heads;
    using baseT::m_imageMJD;
    using baseT::m_imSize;
    using baseT::m_keywords;
    using baseT::m_mask;
    using baseT::m_maskCube;
    using baseT::m_maskFile;
    using baseT::m_minGoodFract;
    using baseT::m_Ncols;
    using baseT::m_Nims;
    using baseT::m_Npix;
    using baseT::m_Nrows;
    using baseT::m_outputDir;
    using baseT::m_pixelTSSigma;
    using baseT::m_prefix;
    using baseT::m_preProcess_azUSM_azHalfWidth;
    using baseT::m_preProcess_azUSM_maxAz;
    using baseT::m_preProcess_azUSM_radHalfWidth;
    using baseT::m_preProcess_beforeCoadd;
    using baseT::m_preProcess_gaussUSM_fwhm;
    using baseT::m_preProcess_mask;
    using baseT::m_preProcess_meanSubMethod;
    using baseT::m_preProcess_medianUSM_fwhm;
    using baseT::m_preProcess_only;
    using baseT::m_preProcess_outputPrefix;
    using baseT::m_preProcess_pixelTSNormMethod;
    using baseT::m_preProcess_subradprof;
    using baseT::m_PSFSubPrefix;
    using baseT::m_qualityFile;
    using baseT::m_qualityThreshold;
    using baseT::m_RDIdateIsISO8601;
    using baseT::m_RDIdateKeyword;
    using baseT::m_RDIdateUnit;
    using baseT::m_RDIdeleteBack;
    using baseT::m_RDIdeleteFront;
    using baseT::m_RDIdirectory;
    using baseT::m_RDIextension;
    using baseT::m_RDIfileList;
    using baseT::m_RDIfileListFile;
    using baseT::m_RDIfilesDeleted;
    using baseT::m_RDIfilesRead;
    using baseT::m_RDIheads;
    using baseT::m_RDIimageMJD;
    using baseT::m_RDIkeywords;
    using baseT::m_RDImask;
    using baseT::m_RDImaskCube;
    using baseT::m_RDImaskFile;
    using baseT::m_RDImaskUseInput;
    using baseT::m_RDIprefix;
    using baseT::m_RDIqualityFile;
    using baseT::m_RDIqualityThreshold;
    using baseT::m_refIms;
    using baseT::m_sigmaThreshold;
    using baseT::m_skipPreProcess;
    using baseT::m_tgtIms;
    using baseT::m_thresholdOnly;
    using baseT::m_weightFile;

    /// Ordered record of target/RDI post-read and post-coadd hook calls.
    std::vector<std::string> m_hookEvents;

    /// Record the target post-read hook.
    void postReadFiles() override;

    /// Record the target post-coadd hook.
    void postCoadd() override;

    /// Record the RDI post-read hook.
    void postRDIReadFiles() override;

    /// Record the RDI post-coadd hook.
    void postRDICoadd() override;
};

/// A unique temporary directory removed when the fixture leaves scope.
class TestDirectory
{
  public:
    /// Create a unique directory below the operating-system temporary directory.
    TestDirectory();

    /// Remove the exact temporary directory created by this fixture.
    ~TestDirectory();

    TestDirectory( const TestDirectory & ) = delete;
    TestDirectory &operator=( const TestDirectory & ) = delete;
    TestDirectory( TestDirectory && ) = delete;
    TestDirectory &operator=( TestDirectory && ) = delete;

    /// Return the temporary directory path.
    const std::filesystem::path &path() const;

    /// Return a path below the temporary directory.
    std::filesystem::path file( const std::string &name /**< [in] relative file name */ ) const;

  private:
    /// Exact path owned by this fixture.
    std::filesystem::path m_path;
};

/// Pin OpenMP's default team size for deterministic numerical tests.
class OpenMPThreadGuard
{
  public:
    /// Save the current setting and request the supplied number of threads.
    explicit OpenMPThreadGuard( int threads /**< [in] positive OpenMP thread count */ = 1 );

    /// Restore the saved OpenMP thread setting.
    ~OpenMPThreadGuard();

    OpenMPThreadGuard( const OpenMPThreadGuard & ) = delete;
    OpenMPThreadGuard &operator=( const OpenMPThreadGuard & ) = delete;

  private:
    /// OpenMP thread setting active before construction.
    int m_previousThreads;
};

/// Write text exactly as supplied to a test file.
void writeTextFile( const std::filesystem::path &path, /**< [in] output file path */
                    const std::string &contents /**< [in] complete file contents */ );

/// Write a small floating-point FITS image and optional header.
void writeFitsImage( const std::filesystem::path &path,              /**< [in] output FITS path */
                     const HCIobservationTestHarness::imageT &image, /**< [in] image pixels */
                     HCIobservationTestHarness::fitsHeaderT *header /**< [in] optional FITS header */ = nullptr );

/// Write a small floating-point FITS cube and optional header.
void writeFitsCube( const std::filesystem::path &path,            /**< [in] output FITS path */
                    const HCIobservationTestHarness::cubeT &cube, /**< [in] cube pixels */
                    HCIobservationTestHarness::fitsHeaderT *header /**< [in] optional FITS header */ = nullptr );

} // namespace unitTest
/// \endcond

#endif // tests_common_HCIobservation_test_fixture_hpp
