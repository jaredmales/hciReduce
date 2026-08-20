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
    double geometryElapsedSeconds{ 0 };   ///< Elapsed serial geometry/setup time in seconds.

    double regressionElapsedSeconds{ 0 }; ///< Elapsed algorithm-worker phase time in seconds.

    double samplingWorkerSeconds{ 0 };    ///< Aggregate worker time used to assemble local input matrices.

    double gramWorkerSeconds{ 0 };        ///< Aggregate worker time used to form normal-equation Gram matrices.

    double eigensolveWorkerSeconds{ 0 };  ///< Aggregate worker time spent in eigensolvers and rank selection.

    double modeWorkerSeconds{ 0 };        ///< Aggregate worker time used to construct modal bases.

    double projectionWorkerSeconds{ 0 };  ///< Aggregate worker time used to apply modes and construct residuals.

    /// Reset every timing value to zero before a new reduction.
    void reset()
    {
        geometryElapsedSeconds = 0;
        regressionElapsedSeconds = 0;
        samplingWorkerSeconds = 0;
        gramWorkerSeconds = 0;
        eigensolveWorkerSeconds = 0;
        modeWorkerSeconds = 0;
        projectionWorkerSeconds = 0;
    }
};

} // namespace improc
} // namespace mx

#endif // ReductionTiming_hpp
