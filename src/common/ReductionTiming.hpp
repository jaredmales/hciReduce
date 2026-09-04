/** \file ReductionTiming.hpp
 * \brief Defines reusable reduction timing records.
 * \author Jared R. Males
 */

#ifndef ReductionTiming_hpp
#define ReductionTiming_hpp

namespace mx
{
namespace improc
{

/// Instance-owned elapsed and aggregate-worker timing values for one reduction algorithm.
/** Worker values sum independent thread-local measurements and can therefore exceed elapsed wall-clock time when work
 * runs in parallel.
 *
 * \ingroup programming_library
 */
struct ReductionTiming
{
    double geometryElapsedSeconds{ 0 };         ///< Elapsed serial geometry/setup time in seconds.

    double regressionElapsedSeconds{ 0 };       ///< Elapsed algorithm-worker phase time in seconds.

    double samplingWorkerSeconds{ 0 };          ///< Aggregate worker time used to assemble local input matrices.

    double sameImageSamplingWorkerSeconds{ 0 }; ///< Aggregate worker time for same-image target and OR sampling.

    double temporalSamplingWorkerSeconds{ 0 };  ///< Aggregate worker time for temporal-summary sampling.

    double gramWorkerSeconds{ 0 };              ///< Aggregate worker time used to form normal-equation Gram matrices.

    double eigensolveWorkerSeconds{ 0 };        ///< Aggregate worker time spent in eigensolvers and rank selection.

    double baseFactorWorkerSeconds{ 0 };        ///< Aggregate worker time used to form reusable base singular systems.

    double deletionWorkerSeconds{ 0 };          ///< Aggregate worker time used by factor-deletion backends.

    double explicitFallbackWorkerSeconds{ 0 };  ///< Aggregate worker time used by explicit numerical fallbacks.

    double modeWorkerSeconds{ 0 };              ///< Aggregate worker time used to construct modal bases.

    double projectionWorkerSeconds{ 0 };        ///< Aggregate worker time used to apply modes and construct residuals.

    double psfWorkerSeconds{ 0 }; ///< Aggregate worker time used to construct compact frozen-model PSF stamps.

    double psfReconstructionElapsedSeconds{ 0 }; ///< Elapsed final spatially variable PSF reconstruction time.

    /// Reset every timing value to zero before a new reduction.
    void reset()
    {
        geometryElapsedSeconds = 0;
        regressionElapsedSeconds = 0;
        samplingWorkerSeconds = 0;
        sameImageSamplingWorkerSeconds = 0;
        temporalSamplingWorkerSeconds = 0;
        gramWorkerSeconds = 0;
        eigensolveWorkerSeconds = 0;
        baseFactorWorkerSeconds = 0;
        deletionWorkerSeconds = 0;
        explicitFallbackWorkerSeconds = 0;
        modeWorkerSeconds = 0;
        projectionWorkerSeconds = 0;
        psfWorkerSeconds = 0;
        psfReconstructionElapsedSeconds = 0;
    }
};

} // namespace improc
} // namespace mx

#endif // ReductionTiming_hpp
