/** \file p4DecompositionBenchmark.cpp
 * \brief Benchmarks resident FP64 symmetric eigendecomposition on CPU and GPU.
 * \author Jared R. Males
 */

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <Eigen/Dense>

#include <cuda_runtime_api.h>
#include <cusolverDn.h>
#include <omp.h>

#include <mx/math/eigenLapack.hpp>

extern "C"
{
    /// Set the process-wide OpenBLAS worker count.
    void openblas_set_num_threads( int numThreads /**< [in] requested OpenBLAS workers */ );

    /// Return the process-wide OpenBLAS worker count.
    int openblas_get_num_threads();
}

namespace
{

/// Column-major FP64 matrix matching LAPACK and cuSOLVER storage.
using matrixT = Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>;

/// Runtime controls for one benchmark invocation.
struct BenchmarkOptions
{
    int dimension{ 620 };       ///< Order of every symmetric matrix.
    int predictorCount{ 2480 }; ///< Columns used to synthesize the Gram matrix.
    int matrixCount{ 64 };      ///< Number of independent decompositions per repetition.
    int streamCount{ 4 };       ///< CUDA streams used in the streamed-batch measurement.
    int cpuThreadCount{ 0 };    ///< OpenMP workers; zero selects the available maximum.
    int repetitions{ 3 };       ///< Timed repetitions of each benchmark variant.
    int device{ 0 };            ///< CUDA logical device.
    std::uint64_t seed{ 42 };   ///< Deterministic random seed.
    bool cpuOnly{ false };      ///< Run CPU measurements without requiring a CUDA device.
};

/// Summary statistics for repeated wall-clock measurements.
struct TimingSummary
{
    double minimumSeconds{ 0 }; ///< Minimum repetition time.
    double meanSeconds{ 0 };    ///< Arithmetic-mean repetition time.
    double stddevSeconds{ 0 };  ///< Population standard deviation.
};

/// CPU-owned state reused for one decomposition at a time.
struct CpuSlot
{
    matrixT matrix;                       ///< Destructive LAPACK input.
    matrixT eigenvectors;                 ///< LAPACK eigenvector output.
    matrixT eigenvalues;                  ///< LAPACK eigenvalue output.
    mx::math::syevrMem<double> workspace; ///< Reusable mxlib/LAPACK workspace.
};

/// Throw a contextual exception for a CUDA runtime failure.
void checkCuda( cudaError_t status,          /**< [in] CUDA status to validate */
                std::string_view operation ) /**< [in] operation producing the status */
{
    if( status != cudaSuccess )
    {
        throw std::runtime_error( std::string( operation ) + ": " + cudaGetErrorString( status ) );
    }
}

/// Throw a contextual exception for a cuSOLVER failure.
void checkCusolver( cusolverStatus_t status,     /**< [in] cuSOLVER status to validate */
                    std::string_view operation ) /**< [in] operation producing the status */
{
    if( status != CUSOLVER_STATUS_SUCCESS )
    {
        throw std::runtime_error( std::string( operation ) + ": cuSOLVER status " +
                                  std::to_string( static_cast<int>( status ) ) );
    }
}

/// Owning device allocation used by the benchmark.
template <typename valueT>
class DeviceBuffer
{
  public:
    /// Construct an empty device allocation.
    DeviceBuffer() = default;

    /// Device allocations have unique ownership.
    DeviceBuffer( const DeviceBuffer & ) = delete;

    /// Device allocations have unique ownership.
    DeviceBuffer &operator=( const DeviceBuffer & ) = delete;

    /// Release the device allocation.
    ~DeviceBuffer()
    {
        if( m_data )
        {
            cudaFree( m_data );
        }
    }

    /// Allocate storage for a number of values.
    void allocate( std::size_t count /**< [in] number of values to allocate */ )
    {
        if( count == 0 )
        {
            return;
        }
        checkCuda( cudaMalloc( reinterpret_cast<void **>( &m_data ), count * sizeof( valueT ) ), "cudaMalloc" );
        m_count = count;
    }

    /// Access the beginning of the allocation.
    valueT *data()
    {
        return m_data;
    }

    /// Access the beginning of the allocation.
    const valueT *data() const
    {
        return m_data;
    }

    /// Return the allocation capacity in values.
    std::size_t size() const
    {
        return m_count;
    }

  private:
    valueT *m_data{ nullptr }; ///< Device pointer owned by this object.
    std::size_t m_count{ 0 };  ///< Allocation capacity in values.
};

/// Owning pinned-host byte allocation for cuSOLVER workspace.
class PinnedBuffer
{
  public:
    /// Construct an empty pinned allocation.
    PinnedBuffer() = default;

    /// Pinned allocations have unique ownership.
    PinnedBuffer( const PinnedBuffer & ) = delete;

    /// Pinned allocations have unique ownership.
    PinnedBuffer &operator=( const PinnedBuffer & ) = delete;

    /// Release the pinned allocation.
    ~PinnedBuffer()
    {
        if( m_data )
        {
            cudaFreeHost( m_data );
        }
    }

    /// Allocate the requested number of pinned bytes.
    void allocate( std::size_t bytes /**< [in] requested capacity in bytes */ )
    {
        if( bytes == 0 )
        {
            return;
        }
        checkCuda( cudaMallocHost( &m_data, bytes ), "cudaMallocHost" );
        m_bytes = bytes;
    }

    /// Access the allocation.
    void *data()
    {
        return m_data;
    }

    /// Return the allocation capacity in bytes.
    std::size_t size() const
    {
        return m_bytes;
    }

  private:
    void *m_data{ nullptr };  ///< Pinned host pointer owned by this object.
    std::size_t m_bytes{ 0 }; ///< Allocation capacity in bytes.
};

/// One independently streamed batched cuSOLVER call.
class CudaBatch
{
  public:
    /// Construct an uninitialized stream batch.
    CudaBatch() = default;

    /// CUDA batches have unique ownership.
    CudaBatch( const CudaBatch & ) = delete;

    /// CUDA batches have unique ownership.
    CudaBatch &operator=( const CudaBatch & ) = delete;

    /// Release CUDA handles and the stream after all queued work completes.
    ~CudaBatch()
    {
        if( m_stream )
        {
            cudaStreamSynchronize( m_stream );
        }
        if( m_params )
        {
            cusolverDnDestroyParams( m_params );
        }
        if( m_solver )
        {
            cusolverDnDestroy( m_solver );
        }
        if( m_stream )
        {
            cudaStreamDestroy( m_stream );
        }
    }

    /// Create a stream, solver state, and workspaces for one contiguous matrix batch.
    void initialize( int dimension,            /**< [in] symmetric matrix order */
                     std::size_t matrixOffset, /**< [in] first matrix in the device arrays */
                     int matrixCount,          /**< [in] number of matrices in this batch */
                     double *matrices,         /**< [in,out] complete device matrix array */
                     double *eigenvalues )     /**< [out] complete device eigenvalue array */
    {
        m_dimension = dimension;
        m_matrixOffset = matrixOffset;
        m_matrixCount = matrixCount;

        const std::uint64_t packedValues = static_cast<std::uint64_t>( dimension ) * dimension * matrixCount;
        if( packedValues > static_cast<std::uint64_t>( std::numeric_limits<std::int32_t>::max() ) )
        {
            throw std::invalid_argument( "cuSOLVER batch exceeds n*lda*batchSize <= INT32_MAX" );
        }

        checkCuda( cudaStreamCreateWithFlags( &m_stream, cudaStreamNonBlocking ), "cudaStreamCreateWithFlags" );
        checkCusolver( cusolverDnCreate( &m_solver ), "cusolverDnCreate" );
        checkCusolver( cusolverDnSetStream( m_solver, m_stream ), "cusolverDnSetStream" );
        checkCusolver( cusolverDnCreateParams( &m_params ), "cusolverDnCreateParams" );

        const std::size_t matrixStride = static_cast<std::size_t>( dimension ) * dimension;
        checkCusolver( cusolverDnXsyevBatched_bufferSize( m_solver,
                                                          m_params,
                                                          CUSOLVER_EIG_MODE_VECTOR,
                                                          CUBLAS_FILL_MODE_LOWER,
                                                          dimension,
                                                          CUDA_R_64F,
                                                          matrices + matrixOffset * matrixStride,
                                                          dimension,
                                                          CUDA_R_64F,
                                                          eigenvalues + matrixOffset * dimension,
                                                          CUDA_R_64F,
                                                          &m_deviceWorkspaceBytes,
                                                          &m_hostWorkspaceBytes,
                                                          matrixCount ),
                       "cusolverDnXsyevBatched_bufferSize" );

        m_deviceWorkspace.allocate( m_deviceWorkspaceBytes );
        m_hostWorkspace.allocate( m_hostWorkspaceBytes );
    }

    /// Queue the eigendecomposition for this stream batch.
    void execute( double *matrices,    /**< [in,out] complete device matrix array */
                  double *eigenvalues, /**< [out] complete device eigenvalue array */
                  int *info )          /**< [out] complete device status array */
    {
        const std::size_t matrixStride = static_cast<std::size_t>( m_dimension ) * m_dimension;
        checkCusolver( cusolverDnXsyevBatched( m_solver,
                                               m_params,
                                               CUSOLVER_EIG_MODE_VECTOR,
                                               CUBLAS_FILL_MODE_LOWER,
                                               m_dimension,
                                               CUDA_R_64F,
                                               matrices + m_matrixOffset * matrixStride,
                                               m_dimension,
                                               CUDA_R_64F,
                                               eigenvalues + m_matrixOffset * m_dimension,
                                               CUDA_R_64F,
                                               m_deviceWorkspace.data(),
                                               m_deviceWorkspaceBytes,
                                               m_hostWorkspace.data(),
                                               m_hostWorkspaceBytes,
                                               info + m_matrixOffset,
                                               m_matrixCount ),
                       "cusolverDnXsyevBatched" );
    }

    /// Wait until all work in this batch's stream completes.
    void synchronize()
    {
        checkCuda( cudaStreamSynchronize( m_stream ), "cudaStreamSynchronize" );
    }

    /// Return the device workspace capacity in bytes.
    std::size_t deviceWorkspaceBytes() const
    {
        return m_deviceWorkspaceBytes;
    }

    /// Return the host workspace capacity in bytes.
    std::size_t hostWorkspaceBytes() const
    {
        return m_hostWorkspaceBytes;
    }

  private:
    int m_dimension{ 0 };                      ///< Symmetric matrix order.
    std::size_t m_matrixOffset{ 0 };           ///< First matrix owned by this batch.
    int m_matrixCount{ 0 };                    ///< Matrices owned by this batch.
    cudaStream_t m_stream{ nullptr };          ///< Nonblocking CUDA stream.
    cusolverDnHandle_t m_solver{ nullptr };    ///< cuSOLVER handle bound to the stream.
    cusolverDnParams_t m_params{ nullptr };    ///< cuSOLVER generic parameters.
    DeviceBuffer<std::byte> m_deviceWorkspace; ///< Device solver workspace.
    PinnedBuffer m_hostWorkspace;              ///< Pinned host solver workspace.
    std::size_t m_deviceWorkspaceBytes{ 0 };   ///< Queried device workspace size.
    std::size_t m_hostWorkspaceBytes{ 0 };     ///< Queried host workspace size.
};

/// Print command-line help.
void printHelp( const char *program /**< [in] executable name */ )
{
    std::cout << "Usage: " << program << " [options]\n"
              << "  --dimension N       matrix order (default 620)\n"
              << "  --predictors K      synthetic Gram predictor count (default 2480)\n"
              << "  --matrices B        decompositions per repetition (default 64)\n"
              << "  --streams S         streams in streamed GPU run (default 4)\n"
              << "  --cpu-threads N     OpenMP CPU workers (default available maximum)\n"
              << "  --repetitions N     timed repetitions (default 3)\n"
              << "  --device N          CUDA logical device (default 0)\n"
              << "  --seed N            deterministic random seed (default 42)\n"
              << "  --cpu-only          run CPU measurements without a CUDA device\n"
              << "  --help              show this message\n";
}

/// Parse a positive integer command-line value.
int parsePositiveInt( const char *value,        /**< [in] value text */
                      std::string_view option ) /**< [in] option name for diagnostics */
{
    const std::string text( value );
    std::size_t parsedCharacters{ 0 };
    const long long parsed = std::stoll( text, &parsedCharacters );
    if( parsedCharacters != text.size() || parsed <= 0 || parsed > std::numeric_limits<int>::max() )
    {
        throw std::invalid_argument( std::string( option ) + " requires a positive integer" );
    }
    return static_cast<int>( parsed );
}

/// Parse a nonnegative integer command-line value.
int parseNonnegativeInt( const char *value,        /**< [in] value text */
                         std::string_view option ) /**< [in] option name for diagnostics */
{
    const std::string text( value );
    std::size_t parsedCharacters{ 0 };
    const long long parsed = std::stoll( text, &parsedCharacters );
    if( parsedCharacters != text.size() || parsed < 0 || parsed > std::numeric_limits<int>::max() )
    {
        throw std::invalid_argument( std::string( option ) + " requires a nonnegative integer" );
    }
    return static_cast<int>( parsed );
}

/// Parse all benchmark command-line options.
BenchmarkOptions parseOptions( int argc,             /**< [in] argument count */
                               char **argv,          /**< [in] argument vector */
                               bool &helpRequested ) /**< [out] whether help was requested */
{
    BenchmarkOptions options;
    helpRequested = false;

    for( int argument = 1; argument < argc; ++argument )
    {
        const std::string_view option( argv[argument] );
        if( option == "--help" )
        {
            helpRequested = true;
            continue;
        }
        if( option == "--cpu-only" )
        {
            options.cpuOnly = true;
            continue;
        }
        if( argument + 1 >= argc )
        {
            throw std::invalid_argument( std::string( option ) + " requires a value" );
        }
        const char *value = argv[++argument];
        if( option == "--dimension" )
        {
            options.dimension = parsePositiveInt( value, option );
        }
        else if( option == "--predictors" )
        {
            options.predictorCount = parsePositiveInt( value, option );
        }
        else if( option == "--matrices" )
        {
            options.matrixCount = parsePositiveInt( value, option );
        }
        else if( option == "--streams" )
        {
            options.streamCount = parsePositiveInt( value, option );
        }
        else if( option == "--cpu-threads" )
        {
            options.cpuThreadCount = parsePositiveInt( value, option );
        }
        else if( option == "--repetitions" )
        {
            options.repetitions = parsePositiveInt( value, option );
        }
        else if( option == "--device" )
        {
            options.device = parseNonnegativeInt( value, option );
        }
        else if( option == "--seed" )
        {
            const std::string text( value );
            std::size_t parsedCharacters{ 0 };
            options.seed = std::stoull( text, &parsedCharacters );
            if( parsedCharacters != text.size() )
            {
                throw std::invalid_argument( "--seed requires an unsigned integer" );
            }
        }
        else
        {
            throw std::invalid_argument( "unknown option: " + std::string( option ) );
        }
    }

    if( options.predictorCount < options.dimension )
    {
        throw std::invalid_argument( "--predictors must be at least --dimension" );
    }
    options.streamCount = std::min( options.streamCount, options.matrixCount );
    return options;
}

/// Calculate timing summary statistics.
TimingSummary summarize( const std::vector<double> &seconds /**< [in] repetition times */ )
{
    TimingSummary summary;
    summary.minimumSeconds = *std::min_element( seconds.begin(), seconds.end() );
    summary.meanSeconds = std::accumulate( seconds.begin(), seconds.end(), 0.0 ) / seconds.size();
    double variance{ 0 };
    for( const double value : seconds )
    {
        const double difference = value - summary.meanSeconds;
        variance += difference * difference;
    }
    summary.stddevSeconds = std::sqrt( variance / seconds.size() );
    return summary;
}

/// Construct a deterministic batch of Gram-like symmetric matrices.
std::vector<double> makeMatrices( const BenchmarkOptions &options /**< [in] matrix dimensions and seed */ )
{
    std::mt19937_64 generator( options.seed );
    std::normal_distribution<double> normal( 0.0, 1.0 );
    matrixT predictors( options.dimension, options.predictorCount );
    for( Eigen::Index column = 0; column < predictors.cols(); ++column )
    {
        for( Eigen::Index row = 0; row < predictors.rows(); ++row )
        {
            predictors( row, column ) = normal( generator );
        }
    }

    matrixT gram =
        ( predictors.matrix() * predictors.matrix().transpose() / static_cast<double>( options.predictorCount ) )
            .array();
    gram.matrix().diagonal().array() += 1.0e-8;

    const std::size_t matrixValues = static_cast<std::size_t>( options.dimension ) * options.dimension;
    std::vector<double> matrices( matrixValues * options.matrixCount );
    for( int matrix = 0; matrix < options.matrixCount; ++matrix )
    {
        std::copy( gram.data(), gram.data() + matrixValues, matrices.data() + matrixValues * matrix );
        for( int diagonal = 0; diagonal < options.dimension; ++diagonal )
        {
            matrices[matrixValues * matrix + static_cast<std::size_t>( diagonal ) * ( options.dimension + 1 )] +=
                static_cast<double>( matrix ) * 1.0e-12;
        }
    }
    return matrices;
}

/// Perform one full CPU symmetric eigendecomposition through the production mxlib API.
int solveCpu( CpuSlot &slot /**< [in,out] matrix, outputs, and reusable workspace */ )
{
    return mx::math::eigenSYEVR( slot.eigenvectors,
                                 slot.eigenvalues,
                                 slot.matrix,
                                 0,
                                 static_cast<int>( slot.matrix.rows() ),
                                 'L',
                                 &slot.workspace );
}

/// Time serial or OpenMP-parallel resident CPU eigendecompositions.
TimingSummary benchmarkCpu( const std::vector<double> &matrices, /**< [in] resident column-major input batch */
                            const BenchmarkOptions &options,     /**< [in] dimensions and repetitions */
                            int threadCount )                    /**< [in] simultaneous CPU decompositions */
{
    const int workers = std::max( 1, std::min( threadCount, options.matrixCount ) );
    const std::size_t matrixValues = static_cast<std::size_t>( options.dimension ) * options.dimension;
    std::vector<std::unique_ptr<CpuSlot>> slots;
    slots.reserve( workers );
    for( int worker = 0; worker < workers; ++worker )
    {
        auto slot = std::make_unique<CpuSlot>();
        slot->matrix.resize( options.dimension, options.dimension );
        slot->eigenvectors.resize( options.dimension, options.dimension );
        slot->eigenvalues.resize( options.dimension, 1 );
        std::copy( matrices.data(), matrices.data() + matrixValues, slot->matrix.data() );
        if( solveCpu( *slot ) != 0 )
        {
            throw std::runtime_error( "CPU eigensolver warmup failed" );
        }
        slots.push_back( std::move( slot ) );
    }

    std::vector<double> repetitionTimes;
    repetitionTimes.reserve( options.repetitions );
    for( int repetition = 0; repetition < options.repetitions; ++repetition )
    {
        double repetitionSeconds{ 0 };
        for( int first = 0; first < options.matrixCount; first += workers )
        {
            const int active = std::min( workers, options.matrixCount - first );
            for( int worker = 0; worker < active; ++worker )
            {
                const double *source = matrices.data() + matrixValues * static_cast<std::size_t>( first + worker );
                std::copy( source, source + matrixValues, slots[worker]->matrix.data() );
            }

            std::atomic<int> solverStatus{ 0 };
            const auto begin = std::chrono::steady_clock::now();
#pragma omp parallel for num_threads( active ) schedule( static, 1 )
            for( int worker = 0; worker < active; ++worker )
            {
                const int status = solveCpu( *slots[worker] );
                if( status != 0 )
                {
                    solverStatus.store( status, std::memory_order_relaxed );
                }
            }
            repetitionSeconds += std::chrono::duration<double>( std::chrono::steady_clock::now() - begin ).count();
            if( solverStatus.load( std::memory_order_relaxed ) != 0 )
            {
                throw std::runtime_error( "CPU eigensolver failed with status " +
                                          std::to_string( solverStatus.load( std::memory_order_relaxed ) ) );
            }
        }
        repetitionTimes.push_back( repetitionSeconds );
    }
    return summarize( repetitionTimes );
}

/// Create equally sized contiguous stream batches, differing by at most one matrix.
std::vector<std::unique_ptr<CudaBatch>> makeCudaBatches( int streamCount,                 /**< [in] number of streams */
                                                         const BenchmarkOptions &options, /**< [in] dimensions */
                                                         double *matrices,     /**< [in,out] device matrices */
                                                         double *eigenvalues ) /**< [out] device eigenvalues */
{
    std::vector<std::unique_ptr<CudaBatch>> batches;
    batches.reserve( streamCount );
    std::size_t offset{ 0 };
    for( int stream = 0; stream < streamCount; ++stream )
    {
        const int count = options.matrixCount / streamCount + ( stream < options.matrixCount % streamCount ? 1 : 0 );
        auto batch = std::make_unique<CudaBatch>();
        batch->initialize( options.dimension, offset, count, matrices, eigenvalues );
        offset += static_cast<std::size_t>( count );
        batches.push_back( std::move( batch ) );
    }
    return batches;
}

/// Queue and synchronize all CUDA stream batches once.
void executeCudaBatches( std::vector<std::unique_ptr<CudaBatch>> &batches, /**< [in,out] stream batches */
                         double *matrices,                                 /**< [in,out] device matrices */
                         double *eigenvalues,                              /**< [out] device eigenvalues */
                         int *info )                                       /**< [out] device statuses */
{
    for( auto &batch : batches )
    {
        batch->execute( matrices, eigenvalues, info );
    }
    for( auto &batch : batches )
    {
        batch->synchronize();
    }
}

/// Time one resident GPU batch layout without including any host/device or device/device copies.
TimingSummary benchmarkGpu( const BenchmarkOptions &options,      /**< [in] dimensions and repetitions */
                            int streamCount,                      /**< [in] concurrent stream count */
                            const DeviceBuffer<double> &pristine, /**< [in] resident pristine device matrices */
                            DeviceBuffer<double> &working,        /**< [in,out] destructive solver matrices */
                            DeviceBuffer<double> &eigenvalues,    /**< [out] device eigenvalues */
                            DeviceBuffer<int> &info,              /**< [out] device solver statuses */
                            std::size_t &deviceWorkspaceBytes,    /**< [out] total workspace across streams */
                            std::size_t &hostWorkspaceBytes )     /**< [out] total workspace across streams */
{
    const std::size_t matrixBytes = pristine.size() * sizeof( double );
    auto batches = makeCudaBatches( streamCount, options, working.data(), eigenvalues.data() );
    deviceWorkspaceBytes = 0;
    hostWorkspaceBytes = 0;
    for( const auto &batch : batches )
    {
        deviceWorkspaceBytes += batch->deviceWorkspaceBytes();
        hostWorkspaceBytes += batch->hostWorkspaceBytes();
    }

    checkCuda( cudaMemcpy( working.data(), pristine.data(), matrixBytes, cudaMemcpyDeviceToDevice ),
               "warmup matrix reset" );
    executeCudaBatches( batches, working.data(), eigenvalues.data(), info.data() );

    std::vector<double> repetitionTimes;
    repetitionTimes.reserve( options.repetitions );
    for( int repetition = 0; repetition < options.repetitions; ++repetition )
    {
        checkCuda( cudaMemcpy( working.data(), pristine.data(), matrixBytes, cudaMemcpyDeviceToDevice ),
                   "timed matrix reset" );
        const auto begin = std::chrono::steady_clock::now();
        executeCudaBatches( batches, working.data(), eigenvalues.data(), info.data() );
        repetitionTimes.push_back( std::chrono::duration<double>( std::chrono::steady_clock::now() - begin ).count() );
    }

    std::vector<int> hostInfo( options.matrixCount );
    checkCuda( cudaMemcpy( hostInfo.data(), info.data(), hostInfo.size() * sizeof( int ), cudaMemcpyDeviceToHost ),
               "copying cuSOLVER status" );
    for( int matrix = 0; matrix < options.matrixCount; ++matrix )
    {
        if( hostInfo[matrix] != 0 )
        {
            throw std::runtime_error( "GPU eigensolver matrix " + std::to_string( matrix ) +
                                      " returned info=" + std::to_string( hostInfo[matrix] ) );
        }
    }
    return summarize( repetitionTimes );
}

/// Solve the first input matrix on the CPU outside the timed measurements.
std::unique_ptr<CpuSlot> makeCpuReference( const std::vector<double> &matrices, /**< [in] resident input batch */
                                           int dimension )                      /**< [in] matrix order */
{
    auto reference = std::make_unique<CpuSlot>();
    const std::size_t matrixValues = static_cast<std::size_t>( dimension ) * dimension;
    reference->matrix.resize( dimension, dimension );
    std::copy( matrices.data(), matrices.data() + matrixValues, reference->matrix.data() );
    if( solveCpu( *reference ) != 0 )
    {
        throw std::runtime_error( "CPU reference eigensolver failed" );
    }
    return reference;
}

/// Report solver accuracy for the first GPU matrix without timing its device-to-host copy.
void reportValidation( const std::vector<double> &matrices,         /**< [in] original input matrices */
                       const CpuSlot &cpuReference,                 /**< [in] CPU eigenpairs */
                       const DeviceBuffer<double> &gpuEigenvectors, /**< [in] overwritten GPU matrices */
                       const DeviceBuffer<double> &gpuEigenvalues,  /**< [in] GPU eigenvalue arrays */
                       int dimension )                              /**< [in] matrix order */
{
    const std::size_t matrixValues = static_cast<std::size_t>( dimension ) * dimension;
    std::vector<double> eigenvectors( matrixValues );
    std::vector<double> eigenvalues( dimension );
    checkCuda( cudaMemcpy( eigenvectors.data(),
                           gpuEigenvectors.data(),
                           matrixValues * sizeof( double ),
                           cudaMemcpyDeviceToHost ),
               "copying validation eigenvectors" );
    checkCuda( cudaMemcpy( eigenvalues.data(),
                           gpuEigenvalues.data(),
                           static_cast<std::size_t>( dimension ) * sizeof( double ),
                           cudaMemcpyDeviceToHost ),
               "copying validation eigenvalues" );

    const Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> input( matrices.data(),
                                                                                         dimension,
                                                                                         dimension );
    const Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> vectors( eigenvectors.data(),
                                                                                           dimension,
                                                                                           dimension );
    const Eigen::Map<const Eigen::VectorXd> values( eigenvalues.data(), dimension );

    double maximumRelativeEigenvalueError{ 0 };
    for( int index = 0; index < dimension; ++index )
    {
        const double scale = std::max( 1.0, std::abs( cpuReference.eigenvalues( index, 0 ) ) );
        maximumRelativeEigenvalueError =
            std::max( maximumRelativeEigenvalueError,
                      std::abs( values( index ) - cpuReference.eigenvalues( index, 0 ) ) / scale );
    }

    const Eigen::MatrixXd residual = input.matrix() * vectors - vectors * values.asDiagonal();
    const Eigen::MatrixXd orthogonality =
        vectors.transpose() * vectors - Eigen::MatrixXd::Identity( dimension, dimension );
    const double normalizedResidual = residual.norm() / ( input.norm() * vectors.norm() );
    const double normalizedOrthogonality = orthogonality.norm() / dimension;

    std::cout << "  Validation (untimed transfers):\n"
              << "    Maximum relative eigenvalue difference: " << maximumRelativeEigenvalueError << "\n"
              << "    GPU normalized residual: " << normalizedResidual << "\n"
              << "    GPU normalized orthogonality error: " << normalizedOrthogonality << "\n";
}

/// Print one timing result in the benchmark report.
void printTiming( std::string_view name,       /**< [in] displayed result name */
                  const TimingSummary &timing, /**< [in] timing statistics */
                  int matrixCount )            /**< [in] decompositions per repetition */
{
    std::cout << "    " << name << ":\n"
              << "      Minimum: " << timing.minimumSeconds << " sec\n"
              << "      Mean: " << timing.meanSeconds << " sec\n"
              << "      Standard deviation: " << timing.stddevSeconds << " sec\n"
              << "      Throughput at mean: " << matrixCount / timing.meanSeconds << " matrices/sec\n";
}

/// Run only the resident CPU measurements for local validation or CPU-only hosts.
int runCpuOnlyBenchmark( const BenchmarkOptions &options /**< [in] validated benchmark controls */ )
{
    const int cpuThreads = options.cpuThreadCount > 0 ? options.cpuThreadCount : omp_get_max_threads();
    openblas_set_num_threads( 1 );
    std::cout << std::setprecision( 10 ) << "Generating resident benchmark matrices...\n";
    const std::vector<double> matrices = makeMatrices( options );
    const TimingSummary cpuSerial = benchmarkCpu( matrices, options, 1 );
    const TimingSummary cpuParallel = benchmarkCpu( matrices, options, cpuThreads );

    std::cout << "P4 resident eigendecomposition benchmark:\n"
              << "  Matrix dimension: " << options.dimension << "\n"
              << "  Synthetic predictor count: " << options.predictorCount << "\n"
              << "  Matrices per repetition: " << options.matrixCount << "\n"
              << "  Repetitions: " << options.repetitions << "\n"
              << "  CPU:\n"
              << "    OpenBLAS threads per solve: " << openblas_get_num_threads() << "\n"
              << "    Parallel decomposition workers: " << std::min( cpuThreads, options.matrixCount ) << "\n";
    printTiming( "Serial", cpuSerial, options.matrixCount );
    printTiming( "OpenMP batch", cpuParallel, options.matrixCount );
    std::cout << "  GPU: skipped by --cpu-only\n";
    return 0;
}

/// Run the complete resident decomposition benchmark.
int runBenchmark( const BenchmarkOptions &options /**< [in] validated benchmark controls */ )
{
    int deviceCount{ 0 };
    checkCuda( cudaGetDeviceCount( &deviceCount ), "cudaGetDeviceCount" );
    if( options.device >= deviceCount )
    {
        throw std::invalid_argument( "requested CUDA device " + std::to_string( options.device ) + " but only " +
                                     std::to_string( deviceCount ) + " device(s) are available" );
    }
    checkCuda( cudaSetDevice( options.device ), "cudaSetDevice" );

    cudaDeviceProp deviceProperties{};
    checkCuda( cudaGetDeviceProperties( &deviceProperties, options.device ), "cudaGetDeviceProperties" );

    const int cpuThreads = options.cpuThreadCount > 0 ? options.cpuThreadCount : omp_get_max_threads();
    openblas_set_num_threads( 1 );

    std::cout << std::setprecision( 10 ) << "Generating resident benchmark matrices...\n";
    const std::vector<double> matrices = makeMatrices( options );
    const std::size_t matrixValues = static_cast<std::size_t>( options.dimension ) * options.dimension;
    const std::size_t matrixBytes = matrices.size() * sizeof( double );

    DeviceBuffer<double> pristine;
    DeviceBuffer<double> working;
    DeviceBuffer<double> eigenvalues;
    DeviceBuffer<int> info;
    pristine.allocate( matrices.size() );
    working.allocate( matrices.size() );
    eigenvalues.allocate( static_cast<std::size_t>( options.matrixCount ) * options.dimension );
    info.allocate( options.matrixCount );
    checkCuda( cudaMemcpy( pristine.data(), matrices.data(), matrixBytes, cudaMemcpyHostToDevice ),
               "preloading benchmark matrices" );

    const TimingSummary cpuSerial = benchmarkCpu( matrices, options, 1 );
    const TimingSummary cpuParallel = benchmarkCpu( matrices, options, cpuThreads );

    std::size_t singleDeviceWorkspace{ 0 };
    std::size_t singleHostWorkspace{ 0 };
    const TimingSummary gpuSingle =
        benchmarkGpu( options, 1, pristine, working, eigenvalues, info, singleDeviceWorkspace, singleHostWorkspace );

    std::size_t streamedDeviceWorkspace{ 0 };
    std::size_t streamedHostWorkspace{ 0 };
    const TimingSummary gpuStreamed = benchmarkGpu( options,
                                                    options.streamCount,
                                                    pristine,
                                                    working,
                                                    eigenvalues,
                                                    info,
                                                    streamedDeviceWorkspace,
                                                    streamedHostWorkspace );

    std::cout << "P4 resident eigendecomposition benchmark:\n"
              << "  Matrix dimension: " << options.dimension << "\n"
              << "  Synthetic predictor count: " << options.predictorCount << "\n"
              << "  Matrices per repetition: " << options.matrixCount << "\n"
              << "  Repetitions: " << options.repetitions << "\n"
              << "  Matrix batch size: " << static_cast<double>( matrixBytes ) / ( 1024.0 * 1024.0 ) << " MiB\n"
              << "  One matrix size: " << static_cast<double>( matrixValues * sizeof( double ) ) / ( 1024.0 * 1024.0 )
              << " MiB\n"
              << "  CPU:\n"
              << "    OpenBLAS threads per solve: " << openblas_get_num_threads() << "\n"
              << "    Parallel decomposition workers: " << std::min( cpuThreads, options.matrixCount ) << "\n";
    printTiming( "Serial", cpuSerial, options.matrixCount );
    printTiming( "OpenMP batch", cpuParallel, options.matrixCount );
    std::cout << "  GPU:\n"
              << "    Device: " << options.device << " (" << deviceProperties.name << ")\n"
              << "    Timed host/device transfers: 0\n"
              << "    Single-batch device workspace: "
              << static_cast<double>( singleDeviceWorkspace ) / ( 1024.0 * 1024.0 ) << " MiB\n"
              << "    Single-batch host workspace: " << static_cast<double>( singleHostWorkspace ) / ( 1024.0 * 1024.0 )
              << " MiB\n";
    printTiming( "One stream, one batch", gpuSingle, options.matrixCount );
    std::cout << "    Streamed batch count: " << options.streamCount << "\n"
              << "    Streamed device workspace total: "
              << static_cast<double>( streamedDeviceWorkspace ) / ( 1024.0 * 1024.0 ) << " MiB\n"
              << "    Streamed host workspace total: "
              << static_cast<double>( streamedHostWorkspace ) / ( 1024.0 * 1024.0 ) << " MiB\n";
    printTiming( "Concurrent stream batches", gpuStreamed, options.matrixCount );

    const std::unique_ptr<CpuSlot> cpuReference = makeCpuReference( matrices, options.dimension );
    reportValidation( matrices, *cpuReference, working, eigenvalues, options.dimension );
    return 0;
}

} // namespace

/// Program entry point for the P4 eigendecomposition benchmark.
int main( int argc,     /**< [in] command-line argument count */
          char **argv ) /**< [in] command-line argument vector */
{
    try
    {
        bool helpRequested{ false };
        const BenchmarkOptions options = parseOptions( argc, argv, helpRequested );
        if( helpRequested )
        {
            printHelp( argv[0] );
            return 0;
        }
        if( options.cpuOnly )
        {
            return runCpuOnlyBenchmark( options );
        }
        return runBenchmark( options );
    }
    catch( const std::exception &exception )
    {
        std::cerr << "p4DecompositionBenchmark: " << exception.what() << "\n";
        return 1;
    }
}
