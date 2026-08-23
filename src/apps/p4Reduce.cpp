/** \file p4Reduce.cpp
 * \brief Defines the p4Reduce command-line application.
 * \author Jared R. Males
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include <mx/app/application.hpp>
#include <mx/ioutils/fits/fitsFile.hpp>

#include "../common/ADIDerotator.hpp"
#include "../common/P4NegativeOptimizer.hpp"
#include "../common/P4Reduction.hpp"

/// Command-line application for Pixel Prediction Post-Processing.
/** The initial application supports one configured target-only P4 reduction. `basic` and `normal` are exact aliases;
 * saved-product postprocessing remains unavailable until the P4 validity-checkpoint format is defined.
 */
class p4Reduce : public mx::app::application
{
  protected:
    std::string m_mode{ "basic" };  ///< Execution mode: basic or normal.

    bool m_showTiming{ false };     ///< Whether to print the completed reduction timing report.

    mx::improc::P4Reductionf m_obs; ///< Configured P4 observation and reduction state.

    /** \name P4 Negative-Companion Optimization - Data
     * @{
     */

    /// One completed delete-one-time-block optimizer refit.
    struct OptimizerJackknifeSample
    {
        std::size_t block{ 0 };               ///< Zero-based omitted block index.

        std::size_t omittedBegin{ 0 };        ///< First omitted physical frame in the half-open interval.

        std::size_t omittedEnd{ 0 };          ///< Physical frame one past the omitted interval.

        std::size_t retainedFrames{ 0 };      ///< Number of physical frames retained by the refit.

        bool converged{ false };              ///< Whether every optimizer convergence condition was satisfied.

        bool denseAgreement{ false };         ///< Whether the contrast result agreed with its dense validation basin.

        mx::improc::P4LocalTrial trial;       ///< Best finite polar trial from the resample.

        double rowDelta{ 0 };                 ///< Fitted Cartesian row displacement from the configured trial.

        double columnDelta{ 0 };              ///< Fitted Cartesian column displacement from the configured trial.

        double merit{ 0 };                    ///< Best aperture mean-square residual in the resample.

        std::size_t evaluationCount{ 0 };     ///< Number of distinct local evaluations used by the refit.

        double evaluationElapsedSeconds{ 0 }; ///< Sum of local-evaluation elapsed seconds for the refit.

        std::string status;                   ///< Stable optimizer completion status for the refit.
    };

    /// Aggregate state and statistics for one requested jackknife uncertainty calculation.
    struct OptimizerJackknifeResult
    {
        std::size_t requestedBlocks{ 0 }; ///< Configured contiguous delete-one-block count; zero disables resampling.

        bool complete{ false };           ///< Whether every requested resample converged and produced covariance.

        std::size_t evaluationCount{ 0 }; ///< Total distinct local evaluations across uncertainty refits.

        double evaluationElapsedSeconds{ 0 }; ///< Sum of elapsed local-evaluation time across uncertainty refits.

        mx::improc::ReductionTiming timing;   ///< Aggregate algorithm timing across uncertainty refits.

        std::vector<OptimizerJackknifeSample> samples; ///< One diagnostic record per requested omitted block.

        std::optional<mx::improc::P4JackknifeStatistics> statistics; ///< Complete ordinary jackknife statistics.

        std::string status{ "disabled" };                            ///< Stable uncertainty completion status.
    };

    bool m_optimizeEnabled{ false };               ///< Whether to run local negative-companion optimization.

    double m_optimizeModeFraction{ 0.15 };         ///< Exact configured P4 mode fraction used by the merit function.

    double m_optimizeApertureRadius{ 0 };          ///< Merit radius; zero selects twice `p4.psfRadius`.

    bool m_optimizeFitPosition{ false };           ///< Whether to fit Cartesian astrometry with contrast.

    double m_optimizeContrastLower{ -0.05 };       ///< Inclusive lower signed contrast bound.

    double m_optimizeContrastUpper{ 0 };           ///< Inclusive upper signed contrast bound.

    std::size_t m_optimizeMaxEvaluations{ 64 };    ///< Maximum number of distinct local P4 calls.

    std::size_t m_optimizeValidationSamples{ 21 }; ///< Uniform dense validation-grid size.

    double m_optimizeParameterTolerance{ 1e-5 };   ///< Absolute contrast convergence tolerance.

    double m_optimizePositionTolerance{ 1e-3 };    ///< Absolute Cartesian convergence tolerance in pixels.

    double m_optimizeMeritTolerance{ 1e-6 };       ///< Combined absolute/relative merit tolerance.

    double m_optimizePositionBound{ 1 };           ///< Symmetric Cartesian position bound in pixels.

    std::size_t m_optimizeUncertaintyBlocks{ 8 };  ///< Contiguous delete-one-block jackknife count; zero disables it.

    std::string m_optimizeOutputPrefix{ "p4Negative_" }; ///< Prefix inside the final image's auxiliary directory.

    /// @}

    /// Return the exact configured output-plane index selected by the optimizer.
    std::size_t optimizerModeIndex() const
    {
        for( std::size_t index = 0; index < m_obs.m_modeFractions.size(); ++index )
        {
            const double fraction = static_cast<double>( m_obs.m_modeFractions[index] );
            const double tolerance = 8 * static_cast<double>( std::numeric_limits<float>::epsilon() ) *
                                     std::max( { 1.0, std::abs( fraction ), std::abs( m_optimizeModeFraction ) } );
            if( std::abs( fraction - m_optimizeModeFraction ) <= tolerance )
            {
                return index;
            }
        }
        throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                              "p4Optimize.modeFraction is not a configured p4.modeFractions entry" );
    }

    /// Return the configured or PSF-radius-derived optimizer aperture radius.
    double optimizerApertureRadius() const
    {
        return m_optimizeApertureRadius > 0 ? m_optimizeApertureRadius : 2 * static_cast<double>( m_obs.m_psfRadius );
    }

    /// Create the optimizer output directory and return its complete product prefix.
    std::filesystem::path optimizerPrefixPath() const
    {
        namespace fs = std::filesystem;
        const fs::path outputDirectory( m_obs.auxiliaryOutputDirectory() );
        const fs::path prefixPath = outputDirectory / m_optimizeOutputPrefix;
        std::error_code directoryError;
        fs::create_directories( prefixPath.parent_path(), directoryError );
        if( directoryError )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::fileoerr,
                                                  "could not create P4 optimizer output directory: " +
                                                      directoryError.message() );
        }
        return prefixPath;
    }

    /// Remove a configuration fragment left by an earlier optimizer run with the same prefix.
    static void removeOptimizerConfiguration(
        const std::filesystem::path &prefixPath /**< [in] complete optimizer product prefix */ )
    {
        const std::filesystem::path configurationPath( prefixPath.string() + "best.conf" );
        std::error_code removalError;
        std::filesystem::remove( configurationPath, removalError );
        if( removalError )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::fileoerr,
                                                  "could not remove stale P4 optimizer configuration fragment " +
                                                      configurationPath.string() + ": " + removalError.message() );
        }
    }

    /// Accumulate every reduction timing field from one optimizer refit.
    static void
    accumulateOptimizerTiming( mx::improc::ReductionTiming &total, /**< [in,out] aggregate uncertainty timing */
                               const mx::improc::ReductionTiming &value /**< [in] timing from one completed refit */ )
    {
        total.geometryElapsedSeconds += value.geometryElapsedSeconds;
        total.regressionElapsedSeconds += value.regressionElapsedSeconds;
        total.samplingWorkerSeconds += value.samplingWorkerSeconds;
        total.sameImageSamplingWorkerSeconds += value.sameImageSamplingWorkerSeconds;
        total.temporalSamplingWorkerSeconds += value.temporalSamplingWorkerSeconds;
        total.gramWorkerSeconds += value.gramWorkerSeconds;
        total.eigensolveWorkerSeconds += value.eigensolveWorkerSeconds;
        total.modeWorkerSeconds += value.modeWorkerSeconds;
        total.projectionWorkerSeconds += value.projectionWorkerSeconds;
        total.psfWorkerSeconds += value.psfWorkerSeconds;
        total.psfReconstructionElapsedSeconds += value.psfReconstructionElapsedSeconds;
    }

    /// Return the ordered physical frames retained after deleting one contiguous block.
    std::vector<int>
    optimizerJackknifeFrames( std::size_t block,      /**< [in] zero-based block to omit */
                              std::size_t blockCount, /**< [in] total equal-count partition block count */
                              std::size_t &begin,     /**< [out] first omitted physical frame */
                              std::size_t &end /**< [out] physical frame one past the omitted interval */ ) const
    {
        const std::size_t frameCount = m_obs.targetFrameCount();
        if( blockCount < 2 || block >= blockCount || blockCount > frameCount )
        {
            throw mx::exception<mx::verbose::vv>(
                mx::error_t::invalidconfig,
                "P4 jackknife blocks must be between two and the loaded frame count" );
        }
        begin = block * frameCount / blockCount;
        end = ( block + 1 ) * frameCount / blockCount;
        std::vector<int> includedFrames;
        includedFrames.reserve( frameCount - ( end - begin ) );
        for( std::size_t frame = 0; frame < frameCount; ++frame )
        {
            if( frame < begin || frame >= end )
            {
                includedFrames.push_back( static_cast<int>( frame ) );
            }
        }
        return includedFrames;
    }

    /// Run the configured contrast-only delete-one-time-block jackknife refits.
    OptimizerJackknifeResult runContrastJackknife(
        const mx::improc::P4ContrastOptimizerConfig &configuration, /**< [in] point-fit optimizer controls */
        const mx::improc::P4LocalTrial &initialTrial,               /**< [in] common configured initial trial */
        std::size_t modeIndex,                                      /**< [in] selected output-plane index */
        double apertureRadius /**< [in] effective merit-aperture radius */ )
    {
        OptimizerJackknifeResult jackknife;
        jackknife.requestedBlocks = m_optimizeUncertaintyBlocks;
        if( jackknife.requestedBlocks == 0 )
        {
            return jackknife;
        }

        std::vector<mx::improc::P4JackknifeEstimate> estimates;
        estimates.reserve( jackknife.requestedBlocks );
        jackknife.samples.reserve( jackknife.requestedBlocks );
        for( std::size_t block = 0; block < jackknife.requestedBlocks; ++block )
        {
            std::size_t begin{ 0 };
            std::size_t end{ 0 };
            const std::vector<int> includedFrames =
                optimizerJackknifeFrames( block, jackknife.requestedBlocks, begin, end );
            std::cerr << "P4 jackknife block " << block + 1 << " / " << jackknife.requestedBlocks
                      << ": omitting frames [" << begin << ',' << end << ")\n";
            const auto evaluate = [this, &initialTrial, &includedFrames]( double contrast )
            {
                mx::improc::P4LocalTrial trial = initialTrial;
                trial.contrast = contrast;
                return m_obs.evaluateLocal( trial, includedFrames );
            };
            const mx::improc::P4ContrastOptimizationResult refit =
                mx::improc::optimizeP4Contrast( configuration, evaluate, modeIndex, apertureRadius );
            mx::improc::P4LocalTrial trial = initialTrial;
            trial.contrast = refit.bestContrast;
            const bool converged = refit.converged && refit.denseAgreement;
            jackknife.samples.push_back( { block,
                                           begin,
                                           end,
                                           includedFrames.size(),
                                           converged,
                                           refit.denseAgreement,
                                           trial,
                                           0,
                                           0,
                                           refit.bestMerit,
                                           refit.evaluationCount,
                                           refit.evaluationElapsedSeconds,
                                           refit.status } );
            jackknife.evaluationCount += refit.evaluationCount;
            jackknife.evaluationElapsedSeconds += refit.evaluationElapsedSeconds;
            accumulateOptimizerTiming( jackknife.timing, refit.timing );
            if( converged )
            {
                estimates.push_back( { 0, 0, refit.bestContrast } );
            }
        }
        jackknife.complete = estimates.size() == jackknife.requestedBlocks;
        if( jackknife.complete )
        {
            jackknife.statistics = mx::improc::p4JackknifeStatistics( estimates, false );
            jackknife.status = "complete";
        }
        else
        {
            jackknife.status = "one or more block refits did not converge";
        }
        return jackknife;
    }

    /// Run the configured joint-position delete-one-time-block jackknife refits.
    OptimizerJackknifeResult runPositionJackknife(
        const mx::improc::P4ContrastOptimizerConfig &configuration, /**< [in] point-fit optimizer controls */
        const mx::improc::P4LocalTrial &referenceTrial,             /**< [in] configured reporting reference */
        const mx::improc::P4LocalTrial &searchCenterTrial,          /**< [in] converged full-data search center */
        std::size_t modeIndex,                                      /**< [in] selected output-plane index */
        double apertureRadius /**< [in] effective merit-aperture radius */ )
    {
        OptimizerJackknifeResult jackknife;
        jackknife.requestedBlocks = m_optimizeUncertaintyBlocks;
        if( jackknife.requestedBlocks == 0 )
        {
            return jackknife;
        }

        std::vector<mx::improc::P4JackknifeEstimate> estimates;
        estimates.reserve( jackknife.requestedBlocks );
        jackknife.samples.reserve( jackknife.requestedBlocks );
        const mx::improc::P4CartesianOffset referenceOffset = mx::improc::p4TrialCartesianOffset( referenceTrial );
        for( std::size_t block = 0; block < jackknife.requestedBlocks; ++block )
        {
            std::size_t begin{ 0 };
            std::size_t end{ 0 };
            const std::vector<int> includedFrames =
                optimizerJackknifeFrames( block, jackknife.requestedBlocks, begin, end );
            std::cerr << "P4 jackknife block " << block + 1 << " / " << jackknife.requestedBlocks
                      << ": omitting frames [" << begin << ',' << end << ")\n";
            const auto evaluate = [this, &includedFrames]( const mx::improc::P4LocalTrial &trial )
            { return m_obs.evaluateLocal( trial, includedFrames ); };
            const mx::improc::P4PositionContrastOptimizationResult refit =
                mx::improc::optimizeP4PositionContrast( configuration,
                                                        searchCenterTrial,
                                                        m_optimizePositionBound,
                                                        evaluate,
                                                        modeIndex,
                                                        apertureRadius );
            const bool converged = refit.converged && refit.denseAgreement;
            const mx::improc::P4CartesianOffset fittedOffset = mx::improc::p4TrialCartesianOffset( refit.bestTrial );
            const double rowDelta = fittedOffset.row - referenceOffset.row;
            const double columnDelta = fittedOffset.column - referenceOffset.column;
            jackknife.samples.push_back( { block,
                                           begin,
                                           end,
                                           includedFrames.size(),
                                           converged,
                                           refit.denseAgreement,
                                           refit.bestTrial,
                                           rowDelta,
                                           columnDelta,
                                           refit.bestMerit,
                                           refit.evaluationCount,
                                           refit.evaluationElapsedSeconds,
                                           refit.status } );
            jackknife.evaluationCount += refit.evaluationCount;
            jackknife.evaluationElapsedSeconds += refit.evaluationElapsedSeconds;
            accumulateOptimizerTiming( jackknife.timing, refit.timing );
            if( converged )
            {
                estimates.push_back( { rowDelta, columnDelta, refit.bestTrial.contrast } );
            }
        }
        jackknife.complete = estimates.size() == jackknife.requestedBlocks;
        if( jackknife.complete )
        {
            jackknife.statistics = mx::improc::p4JackknifeStatistics( estimates, true );
            jackknife.status = "complete";
        }
        else
        {
            jackknife.status = "one or more block refits did not converge";
        }
        return jackknife;
    }

    /// Append jackknife controls and completed uncertainty statistics to one FITS header.
    void appendJackknifeHeader(
        mx::improc::P4Reductionf::fitsHeaderT &header, /**< [in,out] header receiving uncertainty cards */
        const OptimizerJackknifeResult &jackknife /**< [in] completed, disabled, or failed uncertainty outcome */ )
        const
    {
        header.append<unsigned long long>( "P4 OPT JK BLOCKS",
                                           static_cast<unsigned long long>( jackknife.requestedBlocks ),
                                           "delete-one-time-block jackknife count" );
        header.append<int>( "P4 OPT JK COMPLETE", jackknife.complete ? 1 : 0, "jackknife uncertainty completion" );
        header.append<unsigned long long>( "P4 OPT JK EVALUATIONS",
                                           static_cast<unsigned long long>( jackknife.evaluationCount ),
                                           "local evaluations across jackknife refits" );
        header.append<std::string>( "P4 OPT JK STATUS", jackknife.status, "jackknife completion status" );
        if( !jackknife.statistics )
        {
            return;
        }
        const auto &statistics = *jackknife.statistics;
        header.append<double>( "P4 OPT JK ROW ERROR", statistics.rowStandardError, "jackknife row error [pixels]" );
        header.append<double>( "P4 OPT JK COLUMN ERROR",
                               statistics.columnStandardError,
                               "jackknife column error [pixels]" );
        header.append<double>( "P4 OPT JK CONTRAST ERROR",
                               statistics.contrastStandardError,
                               "jackknife signed-contrast error" );
        constexpr std::array<const char *, 9> covarianceNames{ "P4 OPT JK COV RR",
                                                               "P4 OPT JK COV RC",
                                                               "P4 OPT JK COV RX",
                                                               "P4 OPT JK COV CR",
                                                               "P4 OPT JK COV CC",
                                                               "P4 OPT JK COV CX",
                                                               "P4 OPT JK COV XR",
                                                               "P4 OPT JK COV XC",
                                                               "P4 OPT JK COV XX" };
        for( std::size_t index = 0; index < covarianceNames.size(); ++index )
        {
            header.append<double>( covarianceNames[index], statistics.covariance[index], "jackknife covariance" );
        }
    }

    /// Write one row per delete-one-time-block refit when uncertainty resampling is enabled.
    void writeJackknifeTable(
        const OptimizerJackknifeResult &jackknife, /**< [in] uncertainty controls and refit diagnostics */
        const std::filesystem::path &prefixPath /**< [in] complete optimizer product prefix */ ) const
    {
        if( jackknife.requestedBlocks == 0 )
        {
            return;
        }
        const std::filesystem::path path( prefixPath.string() + "jackknife.csv" );
        std::ofstream stream( path );
        if( !stream )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::fileoerr,
                                                  "could not open P4 jackknife table " + path.string() );
        }
        stream << "block,omitted_begin,omitted_end,retained_frames,converged,dense_agreement,row_delta,"
                  "column_delta,separation,position_angle,contrast,merit,evaluations,elapsed_seconds,status\n"
               << std::setprecision( 17 );
        for( const OptimizerJackknifeSample &sample : jackknife.samples )
        {
            stream << sample.block << ',' << sample.omittedBegin << ',' << sample.omittedEnd << ','
                   << sample.retainedFrames << ',' << ( sample.converged ? 1 : 0 ) << ','
                   << ( sample.denseAgreement ? 1 : 0 ) << ',' << sample.rowDelta << ',' << sample.columnDelta << ','
                   << sample.trial.separation << ',' << sample.trial.positionAngle << ',' << sample.trial.contrast
                   << ',' << sample.merit << ',' << sample.evaluationCount << ',' << sample.evaluationElapsedSeconds
                   << ",\"" << sample.status << "\"\n";
        }
        stream.close();
        if( !stream )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::fileoerr,
                                                  "could not write P4 jackknife table " + path.string() );
        }
    }

    /// Append optimizer configuration and outcome provenance to one FITS header.
    /** \tparam resultT contrast-only or joint optimizer result type */
    template <typename resultT>
    void appendOptimizerHeader(
        mx::improc::P4Reductionf::fitsHeaderT &header, /**< [in,out] header receiving optimizer cards */
        const resultT &result,                         /**< [in] completed optimizer outcome */
        const mx::improc::P4LocalTrial &bestTrial,     /**< [in] fitted polar trial */
        bool fitPosition,                              /**< [in] whether astrometry was optimized */
        bool positionConverged,                        /**< [in] position-search convergence */
        bool contrastConverged,                        /**< [in] final contrast-search convergence */
        double rowDelta,                               /**< [in] fitted Cartesian row displacement */
        double columnDelta,                            /**< [in] fitted Cartesian column displacement */
        double apertureRadius,                         /**< [in] effective merit-aperture radius */
        double apertureRow,                            /**< [in] fixed full-image aperture-center row */
        double apertureColumn,                         /**< [in] fixed full-image aperture-center column */
        const OptimizerJackknifeResult &jackknife /**< [in] uncertainty controls and outcome */ ) const
    {
        header.append<int>( "P4 OPTIMIZER", 1, "negative-companion optimizer enabled" );
        header.append<int>( "P4 OPT FIT POSITION", fitPosition ? 1 : 0, "joint Cartesian astrometry enabled" );
        header.append<double>( "P4 OPT MODE FRACTION", m_optimizeModeFraction, "optimized P4 mode fraction" );
        header.append<double>( "P4 OPT APERTURE RADIUS", apertureRadius, "uniform L2 aperture radius [pixels]" );
        header.append<std::string>( "P4 OPT APERTURE FRAME", "FIXED_INITIAL_SKY", "merit-aperture convention" );
        header.append<double>( "P4 OPT APERTURE CENTER ROW", apertureRow, "fixed aperture-center row [pixels]" );
        header.append<double>( "P4 OPT APERTURE CENTER COLUMN",
                               apertureColumn,
                               "fixed aperture-center column [pixels]" );
        header.append<double>( "P4 OPT CONTRAST LOWER", m_optimizeContrastLower, "lower contrast bound" );
        header.append<double>( "P4 OPT CONTRAST UPPER", m_optimizeContrastUpper, "upper contrast bound" );
        header.append<double>( "P4 OPT CONTRAST TOLERANCE",
                               m_optimizeParameterTolerance,
                               "contrast convergence tolerance" );
        header.append<double>( "P4 OPT POSITION BOUND", m_optimizePositionBound, "Cartesian position bound [pixels]" );
        header.append<double>( "P4 OPT POSITION TOLERANCE",
                               m_optimizePositionTolerance,
                               "Cartesian convergence tolerance [pixels]" );
        header.append<double>( "P4 OPT MERIT TOLERANCE", m_optimizeMeritTolerance, "merit convergence tolerance" );
        header.append<double>( "P4 OPT BEST SEPARATION", bestTrial.separation, "fitted separation [pixels]" );
        header.append<double>( "P4 OPT BEST PA", bestTrial.positionAngle, "fitted position angle [degrees]" );
        header.append<double>( "P4 OPT BEST ROW DELTA", rowDelta, "fitted Cartesian row displacement [pixels]" );
        header.append<double>( "P4 OPT BEST COLUMN DELTA",
                               columnDelta,
                               "fitted Cartesian column displacement [pixels]" );
        header.append<double>( "P4 OPT BEST CONTRAST", bestTrial.contrast, "best evaluated signed contrast" );
        header.append<double>( "P4 OPT BEST MERIT", result.bestMerit, "best aperture mean-square residual" );
        header.append<int>( "P4 OPT POSITION CONVERGED", positionConverged ? 1 : 0, "Cartesian minimizer convergence" );
        header.append<int>( "P4 OPT CONTRAST CONVERGED", contrastConverged ? 1 : 0, "final contrast convergence" );
        header.append<int>( "P4 OPT CONVERGED", result.converged ? 1 : 0, "complete optimizer convergence" );
        header.append<int>( "P4 OPT DENSE AGREEMENT", result.denseAgreement ? 1 : 0, "dense-grid basin agreement" );
        header.append<unsigned long long>( "P4 OPT EVALUATIONS",
                                           static_cast<unsigned long long>( result.evaluationCount ),
                                           "distinct local evaluations" );
        header.append<std::string>( "P4 OPT STATUS", result.status, "optimizer completion status" );
        appendJackknifeHeader( header, jackknife );
    }

    /// Write optimizer summary, best local FITS cubes, and a converged configuration fragment.
    /** \tparam resultT contrast-only or joint optimizer result type */
    template <typename resultT>
    void
    writeOptimizerFinalOutputs( const resultT &result,                        /**< [in] completed optimizer outcome */
                                const mx::improc::P4LocalTrial &initialTrial, /**< [in] configured initial trial */
                                const mx::improc::P4LocalTrial &bestTrial,    /**< [in] fitted polar trial */
                                bool fitPosition,       /**< [in] whether astrometry was optimized */
                                bool positionConverged, /**< [in] position-search convergence */
                                bool contrastConverged, /**< [in] final contrast-search convergence */
                                double rowDelta,        /**< [in] fitted Cartesian row displacement */
                                double columnDelta,     /**< [in] fitted Cartesian column displacement */
                                double apertureRadius,  /**< [in] effective merit-aperture radius */
                                double apertureRow,     /**< [in] fixed full-image aperture-center row */
                                double apertureColumn,  /**< [in] fixed full-image aperture-center column */
                                const OptimizerJackknifeResult &jackknife,
                                /**< [in] uncertainty controls and outcome */
                                const std::filesystem::path &prefixPath /**< [in] complete optimizer product prefix */ )
    {
        namespace fs = std::filesystem;

        removeOptimizerConfiguration( prefixPath );
        writeJackknifeTable( jackknife, prefixPath );

        const fs::path summaryPath( prefixPath.string() + "summary.yaml" );
        std::ofstream summaryStream( summaryPath );
        if( !summaryStream )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::fileoerr,
                                                  "could not open P4 optimizer summary " + summaryPath.string() );
        }
        summaryStream << std::setprecision( 17 ) << "p4NegativeOptimization:\n"
                      << "  status: \"" << result.status << "\"\n"
                      << "  converged: " << ( result.converged ? "true" : "false" ) << "\n"
                      << "  fitPosition: " << ( fitPosition ? "true" : "false" ) << "\n"
                      << "  positionConverged: " << ( positionConverged ? "true" : "false" ) << "\n"
                      << "  contrastConverged: " << ( contrastConverged ? "true" : "false" ) << "\n"
                      << "  denseAgreement: " << ( result.denseAgreement ? "true" : "false" ) << "\n"
                      << "  initial:\n"
                      << "    separation: " << initialTrial.separation << "\n"
                      << "    positionAngle: " << initialTrial.positionAngle << "\n"
                      << "    contrast: " << initialTrial.contrast << "\n"
                      << "  fitted:\n"
                      << "    separation: " << bestTrial.separation << "\n"
                      << "    positionAngle: " << bestTrial.positionAngle << "\n"
                      << "    contrast: " << bestTrial.contrast << "\n"
                      << "    rowDelta: " << rowDelta << "\n"
                      << "    columnDelta: " << columnDelta << "\n"
                      << "  modeFraction: " << m_optimizeModeFraction << "\n"
                      << "  apertureRadius: " << apertureRadius << "\n"
                      << "  apertureFrame: \"fixed-initial-sky\"\n"
                      << "  apertureCenterRow: " << apertureRow << "\n"
                      << "  apertureCenterColumn: " << apertureColumn << "\n"
                      << "  positionBound: " << m_optimizePositionBound << "\n"
                      << "  positionTolerance: " << m_optimizePositionTolerance << "\n"
                      << "  contrastTolerance: " << m_optimizeParameterTolerance << "\n"
                      << "  meritTolerance: " << m_optimizeMeritTolerance << "\n"
                      << "  contrastBounds: [" << m_optimizeContrastLower << ", " << m_optimizeContrastUpper << "]\n"
                      << "  evaluationCount: " << result.evaluationCount << "\n"
                      << "  bestMerit: " << result.bestMerit << "\n"
                      << "  evaluationElapsedSeconds: " << result.evaluationElapsedSeconds << "\n"
                      << "  timing:\n"
                      << "    geometryElapsedSeconds: " << result.timing.geometryElapsedSeconds << "\n"
                      << "    regressionElapsedSeconds: " << result.timing.regressionElapsedSeconds << "\n"
                      << "    samplingWorkerSeconds: " << result.timing.samplingWorkerSeconds << "\n"
                      << "    gramWorkerSeconds: " << result.timing.gramWorkerSeconds << "\n"
                      << "    eigensolveWorkerSeconds: " << result.timing.eigensolveWorkerSeconds << "\n"
                      << "    projectionWorkerSeconds: " << result.timing.projectionWorkerSeconds << "\n";
        summaryStream << "  jackknife:\n"
                      << "    requestedBlocks: " << jackknife.requestedBlocks << "\n"
                      << "    complete: " << ( jackknife.complete ? "true" : "false" ) << "\n"
                      << "    status: \"" << jackknife.status << "\"\n"
                      << "    evaluationCount: " << jackknife.evaluationCount << "\n"
                      << "    evaluationElapsedSeconds: " << jackknife.evaluationElapsedSeconds << "\n";
        if( jackknife.statistics )
        {
            const auto &statistics = *jackknife.statistics;
            summaryStream << "    mean:\n"
                          << "      rowDelta: " << statistics.mean.rowDelta << "\n"
                          << "      columnDelta: " << statistics.mean.columnDelta << "\n"
                          << "      contrast: " << statistics.mean.contrast << "\n"
                          << "    standardError:\n"
                          << "      rowDelta: " << statistics.rowStandardError << "\n"
                          << "      columnDelta: " << statistics.columnStandardError << "\n"
                          << "      contrast: " << statistics.contrastStandardError << "\n"
                          << "    covarianceParameterOrder: [rowDelta, columnDelta, contrast]\n"
                          << "    covariance:\n"
                          << "      - [" << statistics.covariance[0] << ", " << statistics.covariance[1] << ", "
                          << statistics.covariance[2] << "]\n"
                          << "      - [" << statistics.covariance[3] << ", " << statistics.covariance[4] << ", "
                          << statistics.covariance[5] << "]\n"
                          << "      - [" << statistics.covariance[6] << ", " << statistics.covariance[7] << ", "
                          << statistics.covariance[8] << "]\n";
        }
        summaryStream.close();
        if( !summaryStream )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::fileoerr,
                                                  "could not write P4 optimizer summary " + summaryPath.string() );
        }

        m_obs.m_fakeSep = { static_cast<float>( bestTrial.separation ) };
        m_obs.m_fakePA = { static_cast<float>( bestTrial.positionAngle ) };
        m_obs.m_fakeContrast = { static_cast<float>( bestTrial.contrast ) };
        m_obs.m_finim = result.bestEvaluation.residual;
        m_obs.m_localFinalValidity = result.bestEvaluation.validity;

        mx::improc::P4Reductionf::fitsHeaderT optimizerHeader;
        m_obs.appendReductionHeader( optimizerHeader );
        appendOptimizerHeader( optimizerHeader,
                               result,
                               bestTrial,
                               fitPosition,
                               positionConverged,
                               contrastConverged,
                               rowDelta,
                               columnDelta,
                               apertureRadius,
                               apertureRow,
                               apertureColumn,
                               jackknife );
        optimizerHeader.append<std::string>( "P4 PRODUCT ROLE", "LOCAL_RESIDUAL", "P4 product role" );
        const fs::path residualPath( prefixPath.string() + "best.fits" );
        m_obs.writeFinimAtPath( residualPath.string(), &optimizerHeader );

        mx::improc::P4Reductionf::fitsHeaderT validityAdditional;
        m_obs.appendReductionHeader( validityAdditional );
        appendOptimizerHeader( validityAdditional,
                               result,
                               bestTrial,
                               fitPosition,
                               positionConverged,
                               contrastConverged,
                               rowDelta,
                               columnDelta,
                               apertureRadius,
                               apertureRow,
                               apertureColumn,
                               jackknife );
        validityAdditional.append<std::string>( "P4 PRODUCT ROLE", "LOCAL_VALIDITY", "P4 product role" );
        mx::improc::P4Reductionf::fitsHeaderT validityHeader;
        m_obs.finalImageHeader( validityHeader, &validityAdditional );
        const fs::path validityPath( prefixPath.string() + "best_validity.fits" );
        mx::fits::fitsFile<float, mx::verbose::vv> writer;
        const mx::error_t validityResult =
            writer.write( validityPath.string(), m_obs.m_localFinalValidity, validityHeader );
        if( validityResult != mx::error_t::noerror )
        {
            throw mx::exception<mx::verbose::vv>( validityResult,
                                                  "could not write P4 optimizer validity " + validityPath.string() );
        }
        std::cerr << "P4 optimizer validity written to: " << validityPath.string() << '\n';

        if( result.converged && result.denseAgreement && ( jackknife.requestedBlocks == 0 || jackknife.complete ) )
        {
            const fs::path configurationPath( prefixPath.string() + "best.conf" );
            std::ofstream configurationStream( configurationPath );
            if( !configurationStream )
            {
                throw mx::exception<mx::verbose::vv>( mx::error_t::fileoerr,
                                                      "could not open P4 optimizer configuration fragment " +
                                                          configurationPath.string() );
            }
            configurationStream << std::setprecision( 17 ) << "[fake]\n"
                                << "sep=" << bestTrial.separation << '\n'
                                << "PA=" << bestTrial.positionAngle << '\n'
                                << "contrast=" << bestTrial.contrast << '\n';
            configurationStream.close();
            if( !configurationStream )
            {
                throw mx::exception<mx::verbose::vv>( mx::error_t::fileoerr,
                                                      "could not write P4 optimizer configuration fragment " +
                                                          configurationPath.string() );
            }
        }
    }

    /// Write contrast-only merit samples and all common optimizer products.
    void writeOptimizerOutputs(
        const mx::improc::P4ContrastOptimizationResult &result, /**< [in] completed contrast-only outcome */
        const mx::improc::P4LocalTrial &initialTrial,           /**< [in] configured initial trial */
        double apertureRadius,                                  /**< [in] effective merit-aperture radius */
        const OptimizerJackknifeResult &jackknife /**< [in] uncertainty controls and outcome */ )
    {
        const std::filesystem::path prefixPath = optimizerPrefixPath();
        const std::filesystem::path meritPath( prefixPath.string() + "merit.csv" );
        std::ofstream meritStream( meritPath );
        if( !meritStream )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::fileoerr,
                                                  "could not open P4 optimizer merit table " + meritPath.string() );
        }
        meritStream << "evaluation,stage,contrast,merit,elapsed_seconds\n" << std::setprecision( 17 );
        for( std::size_t index = 0; index < result.samples.size(); ++index )
        {
            const auto &sample = result.samples[index];
            meritStream << index << ',' << ( sample.denseValidation ? "dense" : "brent" ) << ',' << sample.contrast
                        << ',' << sample.merit << ',' << sample.elapsedSeconds << '\n';
        }
        meritStream.close();
        if( !meritStream )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::fileoerr,
                                                  "could not write P4 optimizer merit table " + meritPath.string() );
        }

        mx::improc::P4LocalTrial bestTrial = initialTrial;
        bestTrial.contrast = result.bestContrast;
        writeOptimizerFinalOutputs( result,
                                    initialTrial,
                                    bestTrial,
                                    false,
                                    false,
                                    result.converged,
                                    0,
                                    0,
                                    apertureRadius,
                                    result.bestEvaluation.sourceRow,
                                    result.bestEvaluation.sourceColumn,
                                    jackknife,
                                    prefixPath );
    }

    /// Write joint position/contrast merit samples and all common optimizer products.
    void writeOptimizerOutputs(
        const mx::improc::P4PositionContrastOptimizationResult &result, /**< [in] completed joint outcome */
        const mx::improc::P4LocalTrial &initialTrial,                   /**< [in] configured initial trial */
        double apertureRadius,                                          /**< [in] effective merit-aperture radius */
        const OptimizerJackknifeResult &jackknife /**< [in] uncertainty controls and outcome */ )
    {
        const std::filesystem::path prefixPath = optimizerPrefixPath();
        const std::filesystem::path meritPath( prefixPath.string() + "merit.csv" );
        std::ofstream meritStream( meritPath );
        if( !meritStream )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::fileoerr,
                                                  "could not open P4 optimizer merit table " + meritPath.string() );
        }
        meritStream << "evaluation,stage,row_delta,column_delta,separation,position_angle,contrast,merit,"
                       "elapsed_seconds\n"
                    << std::setprecision( 17 );
        for( std::size_t index = 0; index < result.samples.size(); ++index )
        {
            const auto &sample = result.samples[index];
            meritStream << index << ',' << sample.stage << ',' << sample.rowDelta << ',' << sample.columnDelta << ','
                        << sample.trial.separation << ',' << sample.trial.positionAngle << ',' << sample.trial.contrast
                        << ',' << sample.merit << ',' << sample.elapsedSeconds << '\n';
        }
        meritStream.close();
        if( !meritStream )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::fileoerr,
                                                  "could not write P4 optimizer merit table " + meritPath.string() );
        }

        writeOptimizerFinalOutputs( result,
                                    initialTrial,
                                    result.bestTrial,
                                    true,
                                    result.positionConverged,
                                    result.contrastConverged,
                                    result.bestRowDelta,
                                    result.bestColumnDelta,
                                    apertureRadius,
                                    result.apertureRow,
                                    result.apertureColumn,
                                    jackknife,
                                    prefixPath );
    }

    /// Print the accumulated optimizer timing hierarchy.
    /** \tparam resultT contrast-only or joint optimizer result type */
    template <typename resultT>
    void printOptimizerTiming( const resultT &result /**< [in] completed optimizer outcome */ ) const
    {
        std::cout << "P4 optimizer times:\n"
                  << "  Local evaluations: " << result.evaluationElapsedSeconds << " elapsed sec\n"
                  << "    Geometry: " << result.timing.geometryElapsedSeconds << " elapsed sec\n"
                  << "    P4 algorithm: " << result.timing.regressionElapsedSeconds << " elapsed sec\n"
                  << "      Sampling: " << result.timing.samplingWorkerSeconds << " worker sec\n"
                  << "      Gram construction: " << result.timing.gramWorkerSeconds << " worker sec\n"
                  << "      EigenDecomposition: " << result.timing.eigensolveWorkerSeconds << " worker sec\n"
                  << "      Projection/residual: " << result.timing.projectionWorkerSeconds << " worker sec\n";
    }

    /// Print uncertainty-refit totals below the point-estimate timing hierarchy.
    void printJackknifeTiming(
        const OptimizerJackknifeResult &jackknife /**< [in] completed or partial uncertainty outcome */ ) const
    {
        if( jackknife.requestedBlocks == 0 )
        {
            return;
        }
        std::cout << "  Jackknife local evaluations: " << jackknife.evaluationElapsedSeconds << " elapsed sec\n"
                  << "    Geometry: " << jackknife.timing.geometryElapsedSeconds << " elapsed sec\n"
                  << "    P4 algorithm: " << jackknife.timing.regressionElapsedSeconds << " elapsed sec\n"
                  << "      Sampling: " << jackknife.timing.samplingWorkerSeconds << " worker sec\n"
                  << "      Gram construction: " << jackknife.timing.gramWorkerSeconds << " worker sec\n"
                  << "      EigenDecomposition: " << jackknife.timing.eigensolveWorkerSeconds << " worker sec\n"
                  << "      Projection/residual: " << jackknife.timing.projectionWorkerSeconds << " worker sec\n";
    }

    /// Print the configured jackknife completion state and available standard errors.
    void printJackknifeResult( const OptimizerJackknifeResult
                                   &jackknife /**< [in] completed, disabled, or failed uncertainty outcome */ ) const
    {
        if( jackknife.requestedBlocks == 0 )
        {
            return;
        }
        std::cout << "P4 jackknife uncertainty: " << jackknife.status << ", blocks=" << jackknife.requestedBlocks
                  << ", evaluations=" << jackknife.evaluationCount;
        if( jackknife.statistics )
        {
            const auto &statistics = *jackknife.statistics;
            if( m_optimizeFitPosition )
            {
                std::cout << ", rowError=" << statistics.rowStandardError
                          << ", columnError=" << statistics.columnStandardError;
            }
            std::cout << ", contrastError=" << statistics.contrastStandardError;
        }
        std::cout << '\n';
    }

    /// Execute the opt-in contrast-only or joint optimizer and persist only its final products.
    int executeOptimizer()
    {
        const std::size_t modeIndex = optimizerModeIndex();
        const double apertureRadius = optimizerApertureRadius();
        const mx::improc::P4LocalTrial initialTrial{ static_cast<double>( m_obs.m_fakeSep[0] ),
                                                     static_cast<double>( m_obs.m_fakePA[0] ),
                                                     static_cast<double>( m_obs.m_fakeContrast[0] ) };
        mx::improc::P4ContrastOptimizerConfig optimizerConfiguration;
        optimizerConfiguration.contrastLower = m_optimizeContrastLower;
        optimizerConfiguration.contrastUpper = m_optimizeContrastUpper;
        optimizerConfiguration.maxEvaluations = m_optimizeMaxEvaluations;
        optimizerConfiguration.validationSamples = m_optimizeValidationSamples;
        optimizerConfiguration.parameterTolerance = m_optimizeParameterTolerance;
        optimizerConfiguration.positionTolerance = m_optimizePositionTolerance;
        optimizerConfiguration.meritTolerance = m_optimizeMeritTolerance;

        if( m_optimizeFitPosition )
        {
            const auto evaluate = [this]( const mx::improc::P4LocalTrial &trial )
            { return m_obs.evaluateLocal( trial ); };
            const mx::improc::P4PositionContrastOptimizationResult result =
                mx::improc::optimizeP4PositionContrast( optimizerConfiguration,
                                                        initialTrial,
                                                        m_optimizePositionBound,
                                                        evaluate,
                                                        modeIndex,
                                                        apertureRadius );
            OptimizerJackknifeResult jackknife;
            jackknife.requestedBlocks = m_optimizeUncertaintyBlocks;
            if( result.converged && result.denseAgreement )
            {
                jackknife = runPositionJackknife( optimizerConfiguration,
                                                  initialTrial,
                                                  result.bestTrial,
                                                  modeIndex,
                                                  apertureRadius );
            }
            else if( jackknife.requestedBlocks != 0 )
            {
                jackknife.status = "skipped because the point fit did not converge";
            }
            writeOptimizerOutputs( result, initialTrial, apertureRadius, jackknife );
            std::cout << std::setprecision( 10 ) << "P4 negative-companion optimization: " << result.status
                      << ", separation=" << result.bestTrial.separation << ", PA=" << result.bestTrial.positionAngle
                      << ", contrast=" << result.bestTrial.contrast << ", merit=" << result.bestMerit
                      << ", evaluations=" << result.evaluationCount << '\n';
            printJackknifeResult( jackknife );
            if( m_showTiming )
            {
                printOptimizerTiming( result );
                printJackknifeTiming( jackknife );
            }
            return result.converged && result.denseAgreement && ( jackknife.requestedBlocks == 0 || jackknife.complete )
                       ? 0
                       : EXIT_FAILURE;
        }

        const auto evaluate = [this, &initialTrial]( double contrast )
        {
            mx::improc::P4LocalTrial trial = initialTrial;
            trial.contrast = contrast;
            return m_obs.evaluateLocal( trial );
        };
        const mx::improc::P4ContrastOptimizationResult result =
            mx::improc::optimizeP4Contrast( optimizerConfiguration, evaluate, modeIndex, apertureRadius );
        OptimizerJackknifeResult jackknife;
        jackknife.requestedBlocks = m_optimizeUncertaintyBlocks;
        if( result.converged && result.denseAgreement )
        {
            jackknife = runContrastJackknife( optimizerConfiguration, initialTrial, modeIndex, apertureRadius );
        }
        else if( jackknife.requestedBlocks != 0 )
        {
            jackknife.status = "skipped because the point fit did not converge";
        }
        writeOptimizerOutputs( result, initialTrial, apertureRadius, jackknife );
        std::cout << std::setprecision( 10 ) << "P4 negative-companion optimization: " << result.status
                  << ", contrast=" << result.bestContrast << ", merit=" << result.bestMerit
                  << ", evaluations=" << result.evaluationCount << '\n';
        printJackknifeResult( jackknife );
        if( m_showTiming )
        {
            printOptimizerTiming( result );
            printJackknifeTiming( jackknife );
        }
        return result.converged && result.denseAgreement && ( jackknife.requestedBlocks == 0 || jackknife.complete )
                   ? 0
                   : EXIT_FAILURE;
    }

  public:
    /// Construct the application with its standard configuration paths.
    p4Reduce()
    {
        m_configPathGlobal_env = "P4REDUCE_GLOBAL_CONFIG";
        m_configPathLocal = "p4Reduce.conf";
        m_requireConfigPathLocal = false;
        config.m_sources = true;
    }

    /// Register application dispatch and P4 reduction configuration targets.
    void setupConfig() override
    {
        config.add( "mode",
                    "",
                    "mode",
                    mx::app::argType::Required,
                    "",
                    "mode",
                    false,
                    "string",
                    "The mode of operation: basic or normal (exact aliases; default basic)" );

        config.add( "showTiming",
                    "",
                    "showTiming",
                    mx::app::argType::True,
                    "",
                    "showTiming",
                    false,
                    "bool",
                    "print the completed reduction timing report to standard output" );

        config.add( "p4Optimize.enabled",
                    "",
                    "p4Optimize.enabled",
                    mx::app::argType::True,
                    "p4Optimize",
                    "enabled",
                    false,
                    "bool",
                    "run bounded negative-companion optimization through the P4 local path" );
        config.add( "p4Optimize.modeFraction",
                    "",
                    "p4Optimize.modeFraction",
                    mx::app::argType::Required,
                    "p4Optimize",
                    "modeFraction",
                    false,
                    "double",
                    "exact configured P4 mode fraction used by the optimizer; default 0.15" );
        config.add( "p4Optimize.apertureRadius",
                    "",
                    "p4Optimize.apertureRadius",
                    mx::app::argType::Required,
                    "p4Optimize",
                    "apertureRadius",
                    false,
                    "double",
                    "positive L2 aperture radius; zero selects twice p4.psfRadius" );
        config.add( "p4Optimize.fitPosition",
                    "",
                    "p4Optimize.fitPosition",
                    mx::app::argType::True,
                    "p4Optimize",
                    "fitPosition",
                    false,
                    "bool",
                    "fit bounded Cartesian position as well as contrast" );
        config.add( "p4Optimize.contrastLower",
                    "",
                    "p4Optimize.contrastLower",
                    mx::app::argType::Required,
                    "p4Optimize",
                    "contrastLower",
                    false,
                    "double",
                    "inclusive lower signed contrast bound; default -0.05" );
        config.add( "p4Optimize.contrastUpper",
                    "",
                    "p4Optimize.contrastUpper",
                    mx::app::argType::Required,
                    "p4Optimize",
                    "contrastUpper",
                    false,
                    "double",
                    "inclusive nonpositive upper signed contrast bound; default 0" );
        config.add( "p4Optimize.maxEvaluations",
                    "",
                    "p4Optimize.maxEvaluations",
                    mx::app::argType::Required,
                    "p4Optimize",
                    "maxEvaluations",
                    false,
                    "size_t",
                    "maximum distinct local P4 calls; default 64" );
        config.add( "p4Optimize.validationSamples",
                    "",
                    "p4Optimize.validationSamples",
                    mx::app::argType::Required,
                    "p4Optimize",
                    "validationSamples",
                    false,
                    "size_t",
                    "number of uniform dense-grid validation samples; default 21" );
        config.add( "p4Optimize.parameterTolerance",
                    "",
                    "p4Optimize.parameterTolerance",
                    mx::app::argType::Required,
                    "p4Optimize",
                    "parameterTolerance",
                    false,
                    "double",
                    "absolute contrast convergence tolerance; default 1e-5" );
        config.add( "p4Optimize.positionTolerance",
                    "",
                    "p4Optimize.positionTolerance",
                    mx::app::argType::Required,
                    "p4Optimize",
                    "positionTolerance",
                    false,
                    "double",
                    "absolute Cartesian-pixel convergence tolerance; default 1e-3" );
        config.add( "p4Optimize.meritTolerance",
                    "",
                    "p4Optimize.meritTolerance",
                    mx::app::argType::Required,
                    "p4Optimize",
                    "meritTolerance",
                    false,
                    "double",
                    "relative dense-grid merit agreement tolerance; default 1e-6" );
        config.add( "p4Optimize.positionBound",
                    "",
                    "p4Optimize.positionBound",
                    mx::app::argType::Required,
                    "p4Optimize",
                    "positionBound",
                    false,
                    "double",
                    "symmetric Cartesian position bound around the initial trial in pixels; default 1" );
        config.add( "p4Optimize.uncertaintyBlocks",
                    "",
                    "p4Optimize.uncertaintyBlocks",
                    mx::app::argType::Required,
                    "p4Optimize",
                    "uncertaintyBlocks",
                    false,
                    "size_t",
                    "contiguous delete-one-block jackknife count; zero disables uncertainty; default 8" );
        config.add( "p4Optimize.outputPrefix",
                    "",
                    "p4Optimize.outputPrefix",
                    mx::app::argType::Required,
                    "p4Optimize",
                    "outputPrefix",
                    false,
                    "string",
                    "optimizer product prefix inside the final image's _outputs directory; default p4Negative_" );

        m_obs.setupConfig( config );
    }

    /// Load application and P4 reduction configuration values and reject unrecognized inputs.
    void loadConfig() override
    {
        m_obs.loadConfig( config );
        config( m_mode, "mode" );
        config( m_showTiming, "showTiming" );
        config( m_optimizeEnabled, "p4Optimize.enabled" );
        config( m_optimizeModeFraction, "p4Optimize.modeFraction" );
        config( m_optimizeApertureRadius, "p4Optimize.apertureRadius" );
        config( m_optimizeFitPosition, "p4Optimize.fitPosition" );
        config( m_optimizeContrastLower, "p4Optimize.contrastLower" );
        config( m_optimizeContrastUpper, "p4Optimize.contrastUpper" );
        config( m_optimizeMaxEvaluations, "p4Optimize.maxEvaluations" );
        config( m_optimizeValidationSamples, "p4Optimize.validationSamples" );
        config( m_optimizeParameterTolerance, "p4Optimize.parameterTolerance" );
        config( m_optimizePositionTolerance, "p4Optimize.positionTolerance" );
        config( m_optimizeMeritTolerance, "p4Optimize.meritTolerance" );
        config( m_optimizePositionBound, "p4Optimize.positionBound" );
        config( m_optimizeUncertaintyBlocks, "p4Optimize.uncertaintyBlocks" );
        config( m_optimizeOutputPrefix, "p4Optimize.outputPrefix" );

        bool unusedPrinted{ false };
        for( const auto &[name, target] : config.m_targets )
        {
            if( target.used )
            {
                continue;
            }

            if( !unusedPrinted )
            {
                std::cerr << "****************************************************\n";
                std::cerr << "WARNING: unused config options (this is a programmer error):\n";
                unusedPrinted = true;
            }
            std::cerr << "   " << name << '\n';
        }

        bool invalidInput{ false };
        if( !config.m_unusedConfigs.empty() )
        {
            std::cerr << "****************************************************\n";
            std::cerr << "ERROR: unrecognized config options:\n";
            for( const auto &[name, target] : config.m_unusedConfigs )
            {
                std::cerr << "   ";
                if( target.keyword.empty() )
                {
                    std::cerr << name;
                }
                else
                {
                    if( !target.section.empty() )
                    {
                        std::cerr << target.section << '.';
                    }
                    std::cerr << target.keyword;
                }
                if( config.m_sources && !target.sources.empty() )
                {
                    std::cerr << " [" << target.sources.front() << ']';
                }
                std::cerr << '\n';
            }
            std::cerr << "****************************************************\n";
            invalidInput = true;
        }

        if( !config.nonOptions.empty() )
        {
            std::cerr << "****************************************************\n";
            std::cerr << "ERROR: unrecognized command line arguments\n";
            invalidInput = true;
        }

        if( invalidInput )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                  "unrecognized p4Reduce configuration input" );
        }
    }

    /// Validate application-level dispatch without duplicating P4 scientific validation.
    void checkConfig() override
    {
        if( m_mode != "basic" && m_mode != "normal" )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                  "invalid p4Reduce mode: " + m_mode +
                                                      "; supported modes are basic and normal" );
        }
        if( !m_optimizeEnabled )
        {
            return;
        }
        if( m_obs.m_localStampSize <= 0 )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                  "p4Optimize.enabled requires positive p4.localStampSize" );
        }
        if( m_obs.m_regressionFrame != mx::improc::P4RegressionFrame::detector )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                  "p4Optimize.enabled supports detector-frame P4 only" );
        }
        if( m_obs.m_fakeSep.size() != 1 || m_obs.m_fakePA.size() != 1 || m_obs.m_fakeContrast.size() != 1 )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                  "p4Optimize.enabled requires one configured fake trial" );
        }
        if( !mx::math::isFinite( m_optimizeContrastLower ) || !mx::math::isFinite( m_optimizeContrastUpper ) ||
            m_optimizeContrastLower >= m_optimizeContrastUpper || m_optimizeContrastUpper > 0 )
        {
            throw mx::exception<mx::verbose::vv>(
                mx::error_t::invalidconfig,
                "p4Optimize contrast bounds must be finite, ordered, and nonpositive" );
        }
        const double initialContrast = static_cast<double>( m_obs.m_fakeContrast[0] );
        if( initialContrast < m_optimizeContrastLower || initialContrast > m_optimizeContrastUpper )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                  "configured fake.contrast is outside p4Optimize contrast bounds" );
        }
        if( m_optimizeValidationSamples < 3 || m_optimizeMaxEvaluations < m_optimizeValidationSamples + 3 )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                  "p4Optimize.maxEvaluations cannot contain the dense scan" );
        }
        if( m_optimizeFitPosition && ( m_optimizePositionBound <= 0 ||
                                       m_optimizeMaxEvaluations < mx::improc::p4PositionContrastMinimumEvaluations(
                                                                      m_optimizeValidationSamples ) ) )
        {
            throw mx::exception<mx::verbose::vv>(
                mx::error_t::invalidconfig,
                "joint p4Optimize requires positive positionBound and enough evaluations for both contrast scans" );
        }
        if( !mx::math::isFinite( m_optimizeParameterTolerance ) || m_optimizeParameterTolerance <= 0 ||
            ( m_optimizeFitPosition &&
              ( !mx::math::isFinite( m_optimizePositionTolerance ) || m_optimizePositionTolerance <= 0 ) ) ||
            !mx::math::isFinite( m_optimizeMeritTolerance ) || m_optimizeMeritTolerance < 0 ||
            !mx::math::isFinite( optimizerApertureRadius() ) || optimizerApertureRadius() <= 0 ||
            !mx::math::isFinite( m_optimizePositionBound ) || m_optimizePositionBound < 0 ||
            ( m_optimizeUncertaintyBlocks != 0 && m_optimizeUncertaintyBlocks < 2 ) || m_optimizeOutputPrefix.empty() ||
            std::filesystem::path( m_optimizeOutputPrefix ).is_absolute() )
        {
            throw mx::exception<mx::verbose::vv>( mx::error_t::invalidconfig,
                                                  "one or more p4Optimize controls are invalid" );
        }
        const double positionSupport = m_optimizeFitPosition ? m_optimizePositionBound : 0;
        if( optimizerApertureRadius() + positionSupport > static_cast<double>( m_obs.m_localStampSize / 2 ) )
        {
            throw mx::exception<mx::verbose::vv>(
                mx::error_t::invalidconfig,
                "p4.localStampSize half-width must be at least p4Optimize.apertureRadius plus positionBound" );
        }
        static_cast<void>( optimizerModeIndex() );
    }

    /// Execute one fully configured P4 reduction.
    /** \returns zero on success, including a configured preprocessing-only stop.
     * \throws mx::exception when configuration, input, reduction, or output fails.
     */
    int execute() override
    {
        if( m_optimizeEnabled )
        {
            return executeOptimizer();
        }
        const int result = m_obs.reduce();
        if( result == 0 && m_showTiming )
        {
            m_obs.dump_times();
        }
        return result;
    }
};

/// Run p4Reduce with command-line exception reporting.
/** The application framework's help status is normalized to `EXIT_SUCCESS`; any reduction status is otherwise
 * propagated. A standard exception is unwound, printed with a help hint, and returned as `EXIT_FAILURE`.
 *
 * \returns `EXIT_SUCCESS` on help or successful execution, otherwise `EXIT_FAILURE`.
 */
int p4ReduceMain( int argc, /**< [in] standard command-line argument count */
                  char **argv /**< [in] standard command-line argument vector */ )
{
    const std::string argv0 = argc > 0 && argv != nullptr && argv[0] != nullptr ? argv[0] : "p4Reduce";
    p4Reduce reduction;

    try
    {
        const int result = reduction.main( argc, argv );
        return result == 1 ? EXIT_SUCCESS : result;
    }
    catch( const std::exception &exception )
    {
        std::vector<std::string> messages;
        mx::unwind_exceptions( messages, exception );
        mx::print_exceptions( messages, std::format( "{}: exception(s) encountered during execution", argv0 ) );
        std::cerr << std::format( "\nTo get help try: {} -h\n\n", argv0 );
        return EXIT_FAILURE;
    }
}

#ifndef HCIREDUCE_P4REDUCE_NO_MAIN
/// p4Reduce process entry point.
int main( int argc, /**< [in] standard command-line argument count */
          char **argv /**< [in] standard command-line argument vector */ )
{
    return p4ReduceMain( argc, argv );
}
#endif
