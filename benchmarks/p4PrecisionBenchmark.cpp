/** \file p4PrecisionBenchmark.cpp
 * \brief Benchmarks CPU precision policies for the uncentered in-sample P4 PCA kernel.
 * \author Jared R. Males
 */

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/SVD>

#include <mx/math/eigenLapack.hpp>

namespace
{

/// Clock used for every benchmark interval.
using benchmarkClockT = std::chrono::steady_clock;

/// Dynamic column-major array used by the benchmark pipeline.
template <typename scalarT>
using matrixT = Eigen::Array<scalarT, Eigen::Dynamic, Eigen::Dynamic>;

/// Dynamic column vector used by the benchmark pipeline.
template <typename scalarT>
using vectorT = Eigen::Array<scalarT, Eigen::Dynamic, 1>;

/// Frozen acceptance-contract schema emitted by this benchmark.
constexpr std::string_view acceptanceSchema{ "hcireduce-pca-precision-acceptance/v1" };

/// Versioned long-form CSV schema emitted by this benchmark.
constexpr std::string_view csvSchema{ "hcireduce-pca-precision-metric/v2" };

/// SHA-256 of the only acceptance contract this executable is authorized to bind.
constexpr std::string_view acceptanceSha256{ "255f70885da886504be1ab9ae43ef5e587cee1ef6454a946116e40da10d1135c" };

/// Incremental SHA-256 calculator used for contract, executable, and canonical source-value identities.
class Sha256
{
  public:
    /// Construct a calculator with the SHA-256 initial state.
    Sha256();

    /// Add bytes to the digest.
    void update( const void *data,   /**< [in] bytes to hash */
                 std::size_t size ); /**< [in] number of bytes */

    /// Finish the digest and return lowercase hexadecimal text.
    std::string finish();

  private:
    /// Apply the SHA-256 compression function to the buffered block.
    void transform();

    /// Current SHA-256 chaining state.
    std::array<std::uint32_t, 8> m_state;

    /// Pending input block.
    std::array<std::uint8_t, 64> m_buffer{};

    /// Number of pending bytes in `m_buffer`.
    std::size_t m_bufferSize{ 0 };

    /// Total unpadded input length in bits.
    std::uint64_t m_bitCount{ 0 };
};

/// Rotate a 32-bit word right by a compile-time-independent count.
std::uint32_t rotateRight( std::uint32_t value, /**< [in] word to rotate */
                           unsigned count )     /**< [in] rotation count */
{
    return ( value >> count ) | ( value << ( 32U - count ) );
}

Sha256::Sha256()
    : m_state{ 0x6a09e667U, 0xbb67ae85U, 0x3c6ef372U, 0xa54ff53aU, 0x510e527fU, 0x9b05688cU, 0x1f83d9abU, 0x5be0cd19U }
{
}

void Sha256::update( const void *data, std::size_t size )
{
    const auto *bytes = static_cast<const std::uint8_t *>( data );
    if( size > ( std::numeric_limits<std::uint64_t>::max() - m_bitCount ) / 8U )
    {
        throw std::overflow_error( "SHA-256 input is too large" );
    }
    m_bitCount += static_cast<std::uint64_t>( size ) * 8U;
    for( std::size_t index = 0; index < size; ++index )
    {
        m_buffer[m_bufferSize++] = bytes[index];
        if( m_bufferSize == m_buffer.size() )
        {
            transform();
            m_bufferSize = 0;
        }
    }
}

void Sha256::transform()
{
    static constexpr std::array<std::uint32_t, 64> constants{
        0x428a2f98U, 0x71374491U, 0xb5c0fbcfU, 0xe9b5dba5U, 0x3956c25bU, 0x59f111f1U, 0x923f82a4U, 0xab1c5ed5U,
        0xd807aa98U, 0x12835b01U, 0x243185beU, 0x550c7dc3U, 0x72be5d74U, 0x80deb1feU, 0x9bdc06a7U, 0xc19bf174U,
        0xe49b69c1U, 0xefbe4786U, 0x0fc19dc6U, 0x240ca1ccU, 0x2de92c6fU, 0x4a7484aaU, 0x5cb0a9dcU, 0x76f988daU,
        0x983e5152U, 0xa831c66dU, 0xb00327c8U, 0xbf597fc7U, 0xc6e00bf3U, 0xd5a79147U, 0x06ca6351U, 0x14292967U,
        0x27b70a85U, 0x2e1b2138U, 0x4d2c6dfcU, 0x53380d13U, 0x650a7354U, 0x766a0abbU, 0x81c2c92eU, 0x92722c85U,
        0xa2bfe8a1U, 0xa81a664bU, 0xc24b8b70U, 0xc76c51a3U, 0xd192e819U, 0xd6990624U, 0xf40e3585U, 0x106aa070U,
        0x19a4c116U, 0x1e376c08U, 0x2748774cU, 0x34b0bcb5U, 0x391c0cb3U, 0x4ed8aa4aU, 0x5b9cca4fU, 0x682e6ff3U,
        0x748f82eeU, 0x78a5636fU, 0x84c87814U, 0x8cc70208U, 0x90befffaU, 0xa4506cebU, 0xbef9a3f7U, 0xc67178f2U,
    };
    std::array<std::uint32_t, 64> words{};
    for( std::size_t index = 0; index < 16; ++index )
    {
        const std::size_t offset = 4 * index;
        words[index] = ( static_cast<std::uint32_t>( m_buffer[offset] ) << 24U ) |
                       ( static_cast<std::uint32_t>( m_buffer[offset + 1] ) << 16U ) |
                       ( static_cast<std::uint32_t>( m_buffer[offset + 2] ) << 8U ) |
                       static_cast<std::uint32_t>( m_buffer[offset + 3] );
    }
    for( std::size_t index = 16; index < words.size(); ++index )
    {
        const std::uint32_t previous15 = words[index - 15];
        const std::uint32_t previous2 = words[index - 2];
        const std::uint32_t sigma0 =
            rotateRight( previous15, 7 ) ^ rotateRight( previous15, 18 ) ^ ( previous15 >> 3U );
        const std::uint32_t sigma1 = rotateRight( previous2, 17 ) ^ rotateRight( previous2, 19 ) ^ ( previous2 >> 10U );
        words[index] = words[index - 16] + sigma0 + words[index - 7] + sigma1;
    }

    std::uint32_t a = m_state[0];
    std::uint32_t b = m_state[1];
    std::uint32_t c = m_state[2];
    std::uint32_t d = m_state[3];
    std::uint32_t e = m_state[4];
    std::uint32_t f = m_state[5];
    std::uint32_t g = m_state[6];
    std::uint32_t h = m_state[7];
    for( std::size_t index = 0; index < words.size(); ++index )
    {
        const std::uint32_t sum1 = rotateRight( e, 6 ) ^ rotateRight( e, 11 ) ^ rotateRight( e, 25 );
        const std::uint32_t choice = ( e & f ) ^ ( ~e & g );
        const std::uint32_t temporary1 = h + sum1 + choice + constants[index] + words[index];
        const std::uint32_t sum0 = rotateRight( a, 2 ) ^ rotateRight( a, 13 ) ^ rotateRight( a, 22 );
        const std::uint32_t majority = ( a & b ) ^ ( a & c ) ^ ( b & c );
        const std::uint32_t temporary2 = sum0 + majority;
        h = g;
        g = f;
        f = e;
        e = d + temporary1;
        d = c;
        c = b;
        b = a;
        a = temporary1 + temporary2;
    }
    m_state[0] += a;
    m_state[1] += b;
    m_state[2] += c;
    m_state[3] += d;
    m_state[4] += e;
    m_state[5] += f;
    m_state[6] += g;
    m_state[7] += h;
}

std::string Sha256::finish()
{
    const std::uint64_t originalBitCount = m_bitCount;
    const std::uint8_t marker{ 0x80U };
    update( &marker, 1 );
    const std::uint8_t zero{ 0 };
    while( m_bufferSize != 56 )
    {
        update( &zero, 1 );
    }
    std::array<std::uint8_t, 8> length{};
    for( std::size_t index = 0; index < length.size(); ++index )
    {
        length[length.size() - 1 - index] = static_cast<std::uint8_t>( originalBitCount >> ( 8U * index ) );
    }
    update( length.data(), length.size() );

    std::ostringstream output;
    output << std::hex << std::setfill( '0' );
    for( const std::uint32_t word : m_state )
    {
        output << std::setw( 8 ) << word;
    }
    return output.str();
}

/// Calculate the SHA-256 digest of a complete file.
std::string sha256File( const std::string &path /**< [in] file to read */ )
{
    std::ifstream input( path, std::ios::binary );
    if( !input )
    {
        throw std::runtime_error( "cannot open file for SHA-256: " + path );
    }
    Sha256 digest;
    std::array<char, 32768> buffer{};
    while( input )
    {
        input.read( buffer.data(), static_cast<std::streamsize>( buffer.size() ) );
        const std::streamsize count = input.gcount();
        if( count > 0 )
        {
            digest.update( buffer.data(), static_cast<std::size_t>( count ) );
        }
    }
    if( !input.eof() )
    {
        throw std::runtime_error( "failed while reading file for SHA-256: " + path );
    }
    return digest.finish();
}

/// Add a 64-bit unsigned integer to a digest using canonical big-endian bytes.
void hashUnsigned( Sha256 &digest,       /**< [in,out] digest to update */
                   std::uint64_t value ) /**< [in] integer to encode */
{
    std::array<std::uint8_t, 8> bytes{};
    for( std::size_t index = 0; index < bytes.size(); ++index )
    {
        bytes[bytes.size() - 1 - index] = static_cast<std::uint8_t>( value >> ( 8U * index ) );
    }
    digest.update( bytes.data(), bytes.size() );
}

/// Add canonical IEEE-754 binary32 values in Eigen storage order to a digest.
void hashFloatValues( Sha256 &digest,          /**< [in,out] digest to update */
                      const float *values,     /**< [in] contiguous values */
                      std::size_t valueCount ) /**< [in] number of values */
{
    static_assert( sizeof( float ) == sizeof( std::uint32_t ) );
    for( std::size_t index = 0; index < valueCount; ++index )
    {
        std::uint32_t bits{ 0 };
        std::memcpy( &bits, values + index, sizeof( bits ) );
        std::array<std::uint8_t, 4> bytes{
            static_cast<std::uint8_t>( bits >> 24U ),
            static_cast<std::uint8_t>( bits >> 16U ),
            static_cast<std::uint8_t>( bits >> 8U ),
            static_cast<std::uint8_t>( bits ),
        };
        digest.update( bytes.data(), bytes.size() );
    }
}

/// Deterministic input families supported by the benchmark.
enum class InputCase
{
    wellConditioned, ///< Existing pseudorandom, well-conditioned smoke input.
    rankBoundary,    ///< Diagonal spectrum with exact zeros and nonzero-threshold boundary values.
    rankZero         ///< Dense exact dependencies evaluated with the production zero threshold.
};

/// Runtime controls and declared engineering smoke-test limits.
struct BenchmarkOptions
{
    int samples{ 96 };                                 ///< Number of temporal samples, T.
    int predictors{ 384 };                             ///< Number of local predictors, K.
    int modes{ 16 };                                   ///< Number of leading modes applied to the target.
    int warmups{ 1 };                                  ///< Untimed pipeline executions per candidate.
    int repetitions{ 6 };                              ///< Timed executions per candidate; six balances rotated order.
    std::uint64_t seed{ 42 };                          ///< Seed for the deterministic FP32 source values.
    double gramRelativeTolerance{ 5e-5 };              ///< Accepted FP32 Gram relative Frobenius error.
    double eigenvalueRelativeTolerance{ 5e-4 };        ///< Accepted retained-eigenvalue relative L2 error.
    double projectorRelativeTolerance{ 2e-2 };         ///< Accepted retained-projector relative Frobenius error.
    double residualRelativeTolerance{ 2e-3 };          ///< Accepted final-residual relative L2 error.
    double solverResidualTolerance{ 2e-4 };            ///< Accepted normalized eigensolver residual.
    double orthogonalityTolerance{ 2e-3 };             ///< Accepted retained-factor orthogonality defect.
    InputCase inputCase{ InputCase::wellConditioned }; ///< Deterministic source family to benchmark.
    double rankRelativeThreshold{ 1e-3 };              ///< Relative eigenvalue threshold for rank diagnostics.
    double rankAmbiguityMultiplier{ 8 };               ///< Multiplier on the precision-scaled ambiguity band.
    std::string acceptanceContractPath;                ///< Runtime path whose content must match the frozen contract.
    std::string expectedAcceptanceSha256{ acceptanceSha256 }; ///< Optional caller pin, constrained to the frozen hash.
    std::string caseId;                                       ///< Stable caller-declared case identifier.
    int innerBlasLapackThreads{ -1 };                         ///< Declared inner BLAS/LAPACK thread count.
    std::string csvOutputPath;                                ///< Optional native CSV artifact path.
    std::string summaryOutputPath;                            ///< Optional native human-summary artifact path.
};

/// Verified identities shared by every metric row in one invocation.
struct RunIdentity
{
    std::string acceptanceContractPath;   ///< Supplied contract path.
    std::string acceptanceContractSha256; ///< Digest verified against the compiled-in frozen digest.
    std::string caseSetId;                ///< Acceptance required-case-set identifier.
    std::string caseId;                   ///< Stable caller-declared case identifier.
    std::string sourceValueSha256;        ///< Digest of canonical source-value serialization.
    std::string executableSha256;         ///< Digest of the running executable.
};

/// One normalized value with its unprotected numerator and protected denominator.
struct NormalizedMetric
{
    double numerator{ 0 };   ///< Norm or scalar in the formula numerator.
    double denominator{ 0 }; ///< Protected scalar denominator used by the formula.
    double value{ 0 };       ///< Quotient of numerator and denominator.
};

/// Independent FP64 thin-SVD reference for ordinary validation gates.
struct NumericalOracle
{
    Eigen::MatrixXd gram;            ///< Direct FP64 Gram from original FP32 source values.
    Eigen::VectorXd eigenvalues;     ///< Retained squared singular values in descending order.
    Eigen::MatrixXd gramModeColumns; ///< Retained left or right singular vectors for the selected Gram.
    Eigen::VectorXd residual;        ///< Retained-left-subspace target residual.
    double targetNorm{ 0 };          ///< Original target L2 norm.
    double residualNorm{ 0 };        ///< Oracle residual L2 norm.
};

/// Long-form acceptance metric description emitted as one CSV record.
struct MetricRecord
{
    std::string record{ "acceptance_metric" }; ///< Record kind within the versioned CSV schema.
    std::string gateId;                        ///< Frozen acceptance gate identifier.
    std::string oracleId;                      ///< Oracle identifier from the acceptance contract.
    std::string metricName;                    ///< Unambiguous metric name.
    std::string numeratorName;                 ///< Formula numerator name, or empty when not meaningful.
    double numeratorValue{ std::numeric_limits<double>::quiet_NaN() }; ///< Formula numerator, or NaN for null.
    std::string numeratorNotApplicableReason; ///< Required explanation when the numerator is null.
    std::string denominatorName;              ///< Formula denominator name, or empty when not meaningful.
    double denominatorValue{ std::numeric_limits<double>::quiet_NaN() }; ///< Protected denominator, or NaN for null.
    std::string denominatorNotApplicableReason;               ///< Required explanation when the denominator is null.
    double value{ std::numeric_limits<double>::quiet_NaN() }; ///< Contract metric value, or NaN for null.
    std::string units;                                        ///< Units spelled as required by the contract.
    std::string aggregationScope;                             ///< Scope over which this raw value was calculated.
    std::string limitOperator;                                ///< Frozen comparison operator, or `null`.
    double limitValue{ std::numeric_limits<double>::quiet_NaN() }; ///< Frozen numeric limit, or NaN for null.
    bool applicable{ true };            ///< Whether this invocation can evaluate the gate metric.
    std::string notApplicableReason;    ///< Required reason when `applicable` is false.
    int passed{ -1 };                   ///< One for pass, zero for fail, or negative for null.
    std::string gateStatus{ "frozen" }; ///< Gate status from the acceptance contract.
    std::string stage{ "validation" };  ///< Validation or timed diagnostic stage.
    double diagnosticValue{ std::numeric_limits<double>::quiet_NaN() }; ///< Non-gating diagnostic value.
    std::string diagnosticUnits;                                        ///< Units of the non-gating diagnostic value.
};

/// One pipeline execution's exclusive stage timings.
struct StageTimings
{
    double conversionSeconds{ 0 };         ///< FP32-source conversion or copy time.
    double gramCrossProductSeconds{ 0 };   ///< Gram and predictor-target cross-product time.
    double eigensolveSeconds{ 0 };         ///< Solver-input conversion, SYEVR, and output conversion time.
    double projectionResidualSeconds{ 0 }; ///< Mode projection and residual update time.
    double totalSeconds{ 0 };              ///< Complete timed pipeline wall time.
};

/// Shared FP32 values presented to every precision candidate.
struct SourceData
{
    matrixT<float> predictors;        ///< T-by-K predictor matrix generated once per invocation.
    vectorT<float> target;            ///< T-sample target generated once per invocation.
    int designedStructuralRank{ -1 }; ///< Known exact source rank, or negative when unspecified.
};

/// Type-erased numerical outputs produced by a separate untimed validation execution.
struct NumericalResult
{
    std::string candidate;        ///< Stable precision-policy identifier.
    matrixT<double> gram;         ///< Selected normal-equation matrix.
    matrixT<double> eigenvectors; ///< Complete ascending Gram eigenvectors.
    vectorT<double> eigenvalues;  ///< Complete ascending Gram eigenvalues.
    vectorT<double> residual;     ///< Residual after the requested leading modes.
};

/// Threshold-boundary values for one ascending normalized eigenvalue sequence.
struct RankBoundary
{
    double ratioBelow{ std::numeric_limits<double>::quiet_NaN() }; ///< Largest ratio at or below the threshold.
    double ratioAbove{ std::numeric_limits<double>::quiet_NaN() }; ///< Smallest ratio strictly above the threshold.
    double eigengap{ std::numeric_limits<double>::quiet_NaN() };   ///< Gap from ratioBelow to ratioAbove.
};

/// Independent FP64 singular-value rank reference for a rank-sensitive input.
struct RankOracle
{
    Eigen::VectorXd normalizedEigenvalues;  ///< Ascending squared singular values normalized by the largest.
    RankBoundary boundary;                  ///< Oracle ratios immediately surrounding the threshold.
    int thresholdRank{ -1 };                ///< Oracle threshold rank, structurally classified at zero tolerance.
    int observedStructuralRank{ -1 };       ///< FP64 SVD rank at its roundoff-scale zero threshold.
    int designedStructuralRank{ -1 };       ///< Exact structural rank established by source construction.
    int exactNullity{ -1 };                 ///< Dimension minus the independently observed structural rank.
    double smallestPositiveRatio{ 0 };      ///< Smallest nonzero oracle lambda/lambda-max ratio.
    double rawLargestNullRatio{ 0 };        ///< Largest raw FP64 SVD ratio classified as exact-null roundoff.
    double gramFrobeniusRatio{ 0 };         ///< Gram Frobenius norm divided by lambda-max.
    double absoluteGramFrobeniusRatio{ 0 }; ///< Absolute-product Gram Frobenius bound divided by lambda-max.
    double minimumThresholdDistance{ 0 };   ///< Smallest realized absolute distance from the rank threshold.
    double lambdaMaximum{ 0 };              ///< Largest FP64 direct-SVD squared singular value.
};

/// Ambiguity-aware numerical-rank diagnostics for one precision candidate.
struct RankMetrics
{
    bool evaluated{ false };                    ///< Whether rank-sensitive diagnostics were requested.
    int oracleRank{ -1 };                       ///< Independent FP64 rank at the requested threshold.
    int candidateRank{ -1 };                    ///< Candidate rank using the production strict comparison.
    int oracleStructuralRank{ -1 };             ///< Independent FP64 SVD structural-rank check.
    int designedStructuralRank{ -1 };           ///< Exact structural rank established by source construction.
    int oracleExactNullity{ -1 };               ///< Independently observed count of exact-zero singular values.
    double oracleSmallestPositiveRatio{ 0 };    ///< Smallest nonzero oracle lambda/lambda-max ratio.
    double oracleRawLargestNullRatio{ 0 };      ///< Largest raw SVD ratio assigned to the exact nullspace.
    double oracleGramFrobeniusRatio{ 0 };       ///< Oracle Gram Frobenius norm divided by lambda-max.
    double oracleAbsoluteGramRatio{ 0 };        ///< Absolute-product Gram bound divided by lambda-max.
    double oracleMinimumThresholdDistance{ 0 }; ///< Smallest realized absolute distance from the rank threshold.
    int ambiguousEigenvalues{ 0 };              ///< Oracle eigenvalues inside this candidate's ambiguity band.
    double ambiguityHalfWidth{ 0 };             ///< Normalized precision-scaled half-width around the threshold.
    RankBoundary oracleBoundary;                ///< Independent FP64 ratios surrounding the threshold.
    RankBoundary candidateBoundary;             ///< Candidate ratios surrounding the threshold.
    bool rankDiffersFromOracle{ false };        ///< Whether the raw candidate and oracle rank counts differ.
    bool agreementOutsideAmbiguity{ true };     ///< Whether every unambiguous retain/reject decision agrees.
    double calculationGamma{ 0 };               ///< Gamma-n term for candidate Gram accumulation.
    double eigensolverGamma{ 0 };               ///< Gamma-n term for candidate eigensolve dimension.
    double calculationBandComponent{ 0 };       ///< Scaled candidate calculation contribution to the band.
    double eigensolverBandComponent{ 0 };       ///< Scaled eigensolver contribution to the band.
    double oracleLambdaMaximum{ 0 };            ///< FP64 direct-SVD lambda maximum.
    std::string ambiguousIndicesDescending;     ///< Semicolon-delimited one-based descending ambiguous indices.
    std::string oracleSpectrumDescending;       ///< Semicolon-delimited normalized FP64 oracle spectrum.
    std::string candidateSpectrumDescending;    ///< Semicolon-delimited normalized candidate spectrum.
    std::string oracleDecisionsDescending;      ///< Semicolon-delimited strict oracle retain decisions.
    std::string candidateDecisionsDescending;   ///< Semicolon-delimited strict candidate retain decisions.
};

/// Numerical comparisons against P4-D64 and intrinsic solver checks.
struct ValidationMetrics
{
    NormalizedMetric gramRelativeError;          ///< Relative Frobenius Gram error.
    NormalizedMetric eigenvalueRelativeError;    ///< Relative L2 retained-eigenvalue error.
    NormalizedMetric projectorRelativeError;     ///< Relative Frobenius retained-projector error.
    NormalizedMetric residualOverOracleResidual; ///< Residual error divided by the FP64 oracle residual norm.
    NormalizedMetric residualOverTarget;         ///< Residual error divided by the original target norm.
    double residualMaximumError{ 0 };            ///< Maximum absolute final-residual error.
    double gramMaximumError{ 0 };                ///< Maximum absolute Gram-entry error.
    NormalizedMetric solverResidual;             ///< Normalized retained eigensolver residual.
    double orthogonalityDefect{ 0 };             ///< Frobenius defect of Q-transpose-Q.
    double maximumPrincipalAngle{ 0 };           ///< Largest retained-subspace principal angle in radians.
    RankMetrics rank;                            ///< Optional ambiguity-aware numerical-rank diagnostics.
    bool passed{ false };                        ///< Whether every declared engineering limit passed.
};

/// Reusable storage for one calculation/eigensolver scalar combination.
template <typename calculationT, typename eigensolverT>
struct CandidateWorkspace
{
    matrixT<calculationT> predictors;                 ///< Candidate calculation-type predictor matrix.
    vectorT<calculationT> target;                     ///< Candidate calculation-type target.
    matrixT<calculationT> gram;                       ///< Candidate calculation-type normal equations.
    vectorT<calculationT> crossProduct;               ///< Predictor-target product for the K-less-than-T branch.
    matrixT<eigensolverT> solverInput;                ///< Destructive eigensolver-type covariance input.
    matrixT<eigensolverT> solverEigenvectors;         ///< Eigensolver-type full eigenvectors.
    matrixT<eigensolverT> solverEigenvalues;          ///< Eigensolver-type full eigenvalues.
    mx::math::syevrMem<eigensolverT> solverWorkspace; ///< Reusable production mxlib SYEVR storage.
    matrixT<calculationT> eigenvectors;               ///< Calculation-type full eigenvectors.
    matrixT<calculationT> eigenvalues;                ///< Calculation-type full eigenvalues.
    vectorT<calculationT> residual;                   ///< Candidate residual scratch and output.
    vectorT<calculationT> projectedMode;              ///< Candidate single-mode projection scratch.
};

/// Repeated timing and worst observed validation values for one candidate.
struct CandidateSummary
{
    std::array<std::vector<double>, 5> stageSeconds; ///< Raw samples for four stages plus total.
    ValidationMetrics maximumMetrics;                ///< Component-wise maximum validation values.
    bool passed{ true };                             ///< Whether every repetition passed.
};

/// Return elapsed seconds since a steady-clock time point.
double secondsSince( const benchmarkClockT::time_point &start /**< [in] beginning of the interval */ )
{
    return std::chrono::duration<double>( benchmarkClockT::now() - start ).count();
}

/// Return the stable machine-readable name of an input family.
std::string_view inputCaseName( InputCase inputCase /**< [in] input family selector */ )
{
    if( inputCase == InputCase::wellConditioned )
    {
        return "well-conditioned";
    }
    return inputCase == InputCase::rankBoundary ? "rank-boundary" : "rank-zero";
}

/// Parse a supported deterministic input-family name.
InputCase parseInputCase( std::string_view value, /**< [in] complete input-family name */
                          std::string_view name ) /**< [in] option name for diagnostics */
{
    if( value == "well-conditioned" )
    {
        return InputCase::wellConditioned;
    }
    if( value == "rank-boundary" )
    {
        return InputCase::rankBoundary;
    }
    if( value == "rank-zero" )
    {
        return InputCase::rankZero;
    }
    throw std::invalid_argument( std::string( name ) + " must be well-conditioned, rank-boundary, or rank-zero" );
}

/// Parse a nonnegative or positive command-line integer into `int`.
int parseInteger( std::string_view value, /**< [in] complete decimal text */
                  std::string_view name,  /**< [in] option name for diagnostics */
                  bool allowZero )        /**< [in] whether zero is accepted */
{
    std::size_t consumed{ 0 };
    long long parsed{ 0 };
    try
    {
        parsed = std::stoll( std::string( value ), &consumed );
    }
    catch( const std::exception & )
    {
        throw std::invalid_argument( std::string( name ) + " requires an integer" );
    }

    if( consumed != value.size() || parsed > std::numeric_limits<int>::max() || parsed < 0 ||
        ( !allowZero && parsed == 0 ) )
    {
        throw std::invalid_argument( std::string( name ) + " is outside its accepted integer range" );
    }
    return static_cast<int>( parsed );
}

/// Parse a complete unsigned decimal seed.
std::uint64_t parseSeed( std::string_view value, /**< [in] complete decimal text */
                         std::string_view name ) /**< [in] option name for diagnostics */
{
    if( value.empty() || value.front() == '-' )
    {
        throw std::invalid_argument( std::string( name ) + " requires an unsigned integer" );
    }
    std::size_t consumed{ 0 };
    unsigned long long parsed{ 0 };
    try
    {
        parsed = std::stoull( std::string( value ), &consumed );
    }
    catch( const std::exception & )
    {
        throw std::invalid_argument( std::string( name ) + " requires an unsigned integer" );
    }

    if( consumed != value.size() )
    {
        throw std::invalid_argument( std::string( name ) + " requires an unsigned integer" );
    }
    return static_cast<std::uint64_t>( parsed );
}

/// Parse a complete finite nonnegative floating-point tolerance.
double parseTolerance( std::string_view value, /**< [in] complete decimal text */
                       std::string_view name ) /**< [in] option name for diagnostics */
{
    std::size_t consumed{ 0 };
    double parsed{ 0 };
    try
    {
        parsed = std::stod( std::string( value ), &consumed );
    }
    catch( const std::exception & )
    {
        throw std::invalid_argument( std::string( name ) + " requires a floating-point value" );
    }

    if( consumed != value.size() || !std::isfinite( parsed ) || parsed < 0 )
    {
        throw std::invalid_argument( std::string( name ) + " requires a finite nonnegative value" );
    }
    return parsed;
}

/// Return whether text is a lowercase or uppercase 64-digit hexadecimal SHA-256.
bool isSha256( std::string_view value /**< [in] text to validate */ )
{
    return value.size() == 64 && std::all_of( value.begin(),
                                              value.end(),
                                              []( char digit )
                                              {
                                                  return ( digit >= '0' && digit <= '9' ) ||
                                                         ( digit >= 'a' && digit <= 'f' ) ||
                                                         ( digit >= 'A' && digit <= 'F' );
                                              } );
}

/// Return whether a case identifier is stable and safe for manifests and CSV consumers.
bool isCaseId( std::string_view value /**< [in] identifier to validate */ )
{
    return !value.empty() && std::all_of( value.begin(),
                                          value.end(),
                                          []( char character )
                                          {
                                              return ( character >= 'a' && character <= 'z' ) ||
                                                     ( character >= 'A' && character <= 'Z' ) ||
                                                     ( character >= '0' && character <= '9' ) || character == '.' ||
                                                     character == '_' || character == ':' || character == '-';
                                          } );
}

/// Return whether two supplied paths resolve to the same filesystem location.
bool samePath( const std::string &first,   /**< [in] first nonempty path */
               const std::string &second ) /**< [in] second nonempty path */
{
    std::error_code firstError;
    std::error_code secondError;
    const std::filesystem::path firstPath = std::filesystem::weakly_canonical( first, firstError );
    const std::filesystem::path secondPath = std::filesystem::weakly_canonical( second, secondError );
    return !firstError && !secondError && firstPath == secondPath;
}

/// Print command-line usage and the declared default numerical smoke gates.
void printUsage( std::ostream &output, /**< [out] destination stream */
                 std::string_view executable /**< [in] invoked program name */ )
{
    output << "Usage: " << executable << " [options]\n"
           << "  --samples T                    Temporal sample count (default 96)\n"
           << "  --predictors K                 Predictor count (default 384)\n"
           << "  --modes N                      Leading retained modes (default 16)\n"
           << "  --warmups N                    Untimed executions per candidate (default 1)\n"
           << "  --repetitions N                Timed executions; multiples of 3 balance order (default 6)\n"
           << "  --seed N                       Deterministic FP32 source seed (default 42)\n"
           << "  --gram-rel-tol X               Frozen Gram engineering smoke gate (5e-5)\n"
           << "  --eigenvalue-rel-tol X         Frozen eigenvalue engineering smoke gate (5e-4)\n"
           << "  --projector-rel-tol X          Frozen projector engineering smoke gate (2e-2)\n"
           << "  --residual-rel-tol X           Frozen residual engineering smoke gate (2e-3)\n"
           << "  --solver-residual-tol X        Frozen eigensolver engineering smoke gate (2e-4)\n"
           << "  --orthogonality-tol X          Frozen orthogonality engineering smoke gate (2e-3)\n"
           << "  --input-case NAME              well-conditioned (default), rank-boundary, or rank-zero\n"
           << "  --rank-rel-tol X               Relative rank threshold (default 1e-3; rank-zero requires 0)\n"
           << "  --rank-band-multiplier X       Frozen precision-band multiplier for rank gates (8)\n"
           << "  --acceptance-contract PATH     Required frozen acceptance-contract file\n"
           << "  --acceptance-contract-sha256 H Optional caller pin; must equal the compiled frozen hash\n"
           << "  --case-id ID                   Required stable [A-Za-z0-9._:-]+ run-case identifier\n"
           << "  --inner-blas-lapack-threads N  Required declared inner thread count; v1 requires 1\n"
           << "  --csv-output PATH              Write CSV directly to PATH instead of stdout\n"
           << "  --summary-output PATH          Write the human summary to PATH instead of stderr\n"
           << "  --help                         Show this message\n\n"
           << "Contract-bound long-form CSV v2 is written to stdout by default; the human summary is written to "
              "stderr by default.\n"
           << "The benchmark hashes the supplied contract, running executable, and canonical FP32 source values.\n"
           << "Defaults are capability smoke gates, not final scientific acceptance tolerances.\n"
           << "Scope is the uncentered, in-sample Gram/eigensolve/projection microbenchmark; rank selection, "
              "centering, held-out fits, coefficients, probes, and product materialization are intentionally "
              "excluded from timing. Rank-sensitive analysis is an untimed, bounded diagnostic; its "
              "precision-scaled ambiguity band is not a science tolerance. Mode application remains folded into "
              "projection/residual timing.\n";
}

/// Parse benchmark controls and reject unknown or incomplete options.
BenchmarkOptions parseOptions( int argc,    /**< [in] command-line argument count */
                               char **argv, /**< [in] command-line argument vector */
                               bool &showHelp /**< [out] whether help was requested */ )
{
    BenchmarkOptions options;
    showHelp = false;
    for( int argument = 1; argument < argc; ++argument )
    {
        const std::string_view name( argv[argument] );
        if( name == "--help" )
        {
            showHelp = true;
            continue;
        }
        if( argument + 1 >= argc )
        {
            throw std::invalid_argument( std::string( name ) + " requires a value" );
        }
        const std::string_view value( argv[++argument] );
        if( name == "--samples" )
        {
            options.samples = parseInteger( value, name, false );
        }
        else if( name == "--predictors" )
        {
            options.predictors = parseInteger( value, name, false );
        }
        else if( name == "--modes" )
        {
            options.modes = parseInteger( value, name, false );
        }
        else if( name == "--warmups" )
        {
            options.warmups = parseInteger( value, name, true );
        }
        else if( name == "--repetitions" )
        {
            options.repetitions = parseInteger( value, name, false );
        }
        else if( name == "--seed" )
        {
            options.seed = parseSeed( value, name );
        }
        else if( name == "--gram-rel-tol" )
        {
            options.gramRelativeTolerance = parseTolerance( value, name );
        }
        else if( name == "--eigenvalue-rel-tol" )
        {
            options.eigenvalueRelativeTolerance = parseTolerance( value, name );
        }
        else if( name == "--projector-rel-tol" )
        {
            options.projectorRelativeTolerance = parseTolerance( value, name );
        }
        else if( name == "--residual-rel-tol" )
        {
            options.residualRelativeTolerance = parseTolerance( value, name );
        }
        else if( name == "--solver-residual-tol" )
        {
            options.solverResidualTolerance = parseTolerance( value, name );
        }
        else if( name == "--orthogonality-tol" )
        {
            options.orthogonalityTolerance = parseTolerance( value, name );
        }
        else if( name == "--input-case" )
        {
            options.inputCase = parseInputCase( value, name );
        }
        else if( name == "--rank-rel-tol" )
        {
            options.rankRelativeThreshold = parseTolerance( value, name );
        }
        else if( name == "--rank-band-multiplier" )
        {
            options.rankAmbiguityMultiplier = parseTolerance( value, name );
        }
        else if( name == "--acceptance-contract" )
        {
            options.acceptanceContractPath = value;
        }
        else if( name == "--acceptance-contract-sha256" )
        {
            options.expectedAcceptanceSha256 = value;
        }
        else if( name == "--case-id" )
        {
            options.caseId = value;
        }
        else if( name == "--inner-blas-lapack-threads" )
        {
            options.innerBlasLapackThreads = parseInteger( value, name, false );
        }
        else if( name == "--csv-output" )
        {
            options.csvOutputPath = value;
        }
        else if( name == "--summary-output" )
        {
            options.summaryOutputPath = value;
        }
        else
        {
            throw std::invalid_argument( "unknown option: " + std::string( name ) );
        }
    }

    if( showHelp )
    {
        return options;
    }
    if( options.acceptanceContractPath.empty() )
    {
        throw std::invalid_argument( "--acceptance-contract is required for contract-bound CSV v2 output" );
    }
    if( !isSha256( options.expectedAcceptanceSha256 ) || options.expectedAcceptanceSha256 != acceptanceSha256 )
    {
        throw std::invalid_argument( "--acceptance-contract-sha256 does not match the compiled frozen contract" );
    }
    if( !isCaseId( options.caseId ) )
    {
        throw std::invalid_argument( "--case-id is required and must contain only [A-Za-z0-9._:-]" );
    }
    if( options.innerBlasLapackThreads != 1 )
    {
        throw std::invalid_argument( "--inner-blas-lapack-threads 1 is required by the v1 acceptance contract" );
    }
    if( options.gramRelativeTolerance != 5e-5 || options.eigenvalueRelativeTolerance != 5e-4 ||
        options.projectorRelativeTolerance != 2e-2 || options.residualRelativeTolerance != 2e-3 ||
        options.solverResidualTolerance != 2e-4 || options.orthogonalityTolerance != 2e-3 )
    {
        throw std::invalid_argument( "contract-bound mode requires the frozen v1 engineering smoke limits" );
    }
    if( options.inputCase != InputCase::wellConditioned && options.rankAmbiguityMultiplier != 8.0 )
    {
        throw std::invalid_argument( "rank-sensitive contract-bound mode requires --rank-band-multiplier 8" );
    }
    if( !options.csvOutputPath.empty() && !options.summaryOutputPath.empty() &&
        samePath( options.csvOutputPath, options.summaryOutputPath ) )
    {
        throw std::invalid_argument( "--csv-output and --summary-output must identify different files" );
    }
    if( ( !options.csvOutputPath.empty() && samePath( options.csvOutputPath, options.acceptanceContractPath ) ) ||
        ( !options.summaryOutputPath.empty() &&
          samePath( options.summaryOutputPath, options.acceptanceContractPath ) ) )
    {
        throw std::invalid_argument( "output paths must not overwrite the acceptance contract" );
    }

    const int maximumModes = std::min( options.samples, options.predictors );
    if( options.modes > maximumModes )
    {
        throw std::invalid_argument( "--modes must not exceed min(T,K)" );
    }
    if( options.inputCase != InputCase::wellConditioned )
    {
        constexpr int reservedBoundaryValues{ 8 };
        constexpr long long maximumSourceElements{ 1048576 };
        constexpr int maximumGramDimension{ 256 };
        if( maximumModes < reservedBoundaryValues + 4 )
        {
            throw std::invalid_argument( "rank-sensitive inputs require min(T,K) >= 12" );
        }
        if( options.modes > maximumModes - reservedBoundaryValues )
        {
            throw std::invalid_argument( "rank-sensitive --modes must not enter the eight-value diagnostic tail" );
        }
        if( maximumModes > maximumGramDimension ||
            static_cast<long long>( options.samples ) * options.predictors > maximumSourceElements )
        {
            throw std::invalid_argument( "rank-sensitive inputs are bounded to min(T,K) <= 256 and T*K <= 1048576" );
        }
        if( options.inputCase == InputCase::rankBoundary &&
            ( options.rankRelativeThreshold <= 0 || options.rankRelativeThreshold >= 0.25 ) )
        {
            throw std::invalid_argument( "rank-boundary --rank-rel-tol must be in (0,0.25)" );
        }
        const double oracleZeroScale =
            std::max( options.samples, options.predictors ) * std::numeric_limits<double>::epsilon();
        const double minimumBoundaryThreshold = 4e4 * oracleZeroScale * oracleZeroScale;
        if( options.inputCase == InputCase::rankBoundary && options.rankRelativeThreshold <= minimumBoundaryThreshold )
        {
            throw std::invalid_argument(
                "rank-boundary --rank-rel-tol is too small to separate its near-positive mode from the bounded "
                "FP64 structural-rank threshold" );
        }
        if( options.inputCase == InputCase::rankZero && options.rankRelativeThreshold != 0 )
        {
            throw std::invalid_argument( "rank-zero requires --rank-rel-tol 0" );
        }
        if( options.rankAmbiguityMultiplier <= 0 )
        {
            throw std::invalid_argument( "rank-sensitive --rank-band-multiplier must be positive" );
        }
    }
    return options;
}

/// Generate one exactly representable pseudorandom FP32 value in [-1,1).
float nextUnitFloat( std::mt19937 &generator /**< [in,out] deterministic integer generator */ )
{
    const std::uint32_t mantissa = generator() >> 8;
    return static_cast<float>( mantissa ) * ( 1.0f / 8388608.0f ) - 1.0f;
}

/// Return Higham's gamma-n accumulation factor using a conservative machine-epsilon unit.
template <typename scalarT>
double gammaN( int operations /**< [in] number of rounded operations in the accumulation model */ )
{
    const double scaledEpsilon = operations * static_cast<double>( std::numeric_limits<scalarT>::epsilon() );
    if( scaledEpsilon >= 1 )
    {
        throw std::invalid_argument( "rank ambiguity model requires n*epsilon < 1" );
    }
    return scaledEpsilon / ( 1.0 - scaledEpsilon );
}

/// Return a candidate-specific diagnostic ambiguity half-width in lambda/lambda-max units.
/** For `makeRankBoundarySourceData`, the diagonal embedding leaves at most one nonzero product in each Gram entry,
 * so the calculation term is gamma-1. `makeRankZeroSourceData` uses dense dot products and therefore the full inner
 * dimension. The Gram-formation term is scaled by the FP64 oracle's Frobenius norm of the absolute-product Gram,
 * while the eigensolver gamma-N term is scaled by the ordinary Gram Frobenius norm. Both are normalized by
 * lambda-max, and the factor of two covers normalization by the computed largest eigenvalue. The configurable
 * multiplier absorbs solver constants. This remains a diagnostic exclusion band for these bounded constructions,
 * not a general forward-error theorem or science tolerance. Placement-only calls use unit scales for near-boundary
 * values and the diagonal sqrt(N) worst case for the far value; validation supplies measured oracle scales.
 */
template <typename calculationT, typename eigensolverT>
double rankAmbiguityHalfWidth(
    const BenchmarkOptions &options,        /**< [in] dimensions and band multiplier */
    double calculationMagnitudeScale = 1.0, /**< [in] absolute-product Gram norm divided by lambda-max */
    double solverMagnitudeScale = 1.0 )     /**< [in] Gram norm divided by lambda-max */
{
    const int gramDimension = std::min( options.samples, options.predictors );
    const int innerDimension = options.samples <= options.predictors ? options.predictors : options.samples;
    const int calculationOperations = options.inputCase == InputCase::rankZero ? innerDimension : 1;
    const double halfWidth = 2.0 * options.rankAmbiguityMultiplier *
                             ( gammaN<calculationT>( calculationOperations ) * calculationMagnitudeScale +
                               gammaN<eigensolverT>( gramDimension ) * solverMagnitudeScale );
    if( !std::isfinite( halfWidth ) )
    {
        throw std::invalid_argument( "rank ambiguity half-width is not finite" );
    }
    return halfWidth;
}

/// Construct one deterministic P4-like FP32 predictor and correlated target pair.
SourceData makeSourceData( const BenchmarkOptions &options /**< [in] dimensions and seed */ )
{
    SourceData source;
    source.predictors.resize( options.samples, options.predictors );
    source.target.resize( options.samples );
    std::seed_seq seedSequence{ static_cast<std::uint32_t>( options.seed ),
                                static_cast<std::uint32_t>( options.seed >> 32 ) };
    std::mt19937 generator( seedSequence );

    for( Eigen::Index column = 0; column < source.predictors.cols(); ++column )
    {
        const float position = static_cast<float>( column + 1 ) / static_cast<float>( source.predictors.cols() );
        const float scale = 0.75f + 0.5f * position;
        for( Eigen::Index row = 0; row < source.predictors.rows(); ++row )
        {
            source.predictors( row, column ) = scale * nextUnitFloat( generator );
        }
    }

    const Eigen::Index last = source.predictors.cols() - 1;
    const Eigen::Index upperMiddle = 3 * source.predictors.cols() / 4;
    const Eigen::Index middle = source.predictors.cols() / 2;
    for( Eigen::Index row = 0; row < source.target.rows(); ++row )
    {
        source.target( row ) = 0.55f * source.predictors( row, last ) - 0.30f * source.predictors( row, upperMiddle ) +
                               0.20f * source.predictors( row, middle ) + 0.08f * nextUnitFloat( generator );
    }
    return source;
}

/// Construct a bounded FP32 matrix with exact zeros and eigenvalues straddling the requested rank threshold.
/** The diagonal embedding gives the FP32 source a known exact structural rank. Its realized spectrum is still
 * reconstructed independently from the FP32 values with an FP64 SVD before rank decisions are checked. Near-tail
 * values use the unit-scale nominal band, while the far-above value uses the worst possible sqrt(N) Frobenius scale
 * for a normalized diagonal Gram. Casting the requested square roots to FP32 does not generally realize exact
 * equality with the threshold, so the oracle also reports the minimum realized threshold distance.
 */
SourceData
makeRankBoundarySourceData( const BenchmarkOptions &options /**< [in] dimensions, seed, and rank controls */ )
{
    constexpr int boundaryValues{ 8 };
    const int dimension = std::min( options.samples, options.predictors );
    const int robustValues = dimension - boundaryValues;
    const int innerDimension = options.samples <= options.predictors ? options.predictors : options.samples;
    const double nominalFloatBand = rankAmbiguityHalfWidth<float, float>( options );
    const double outerFloatBand = rankAmbiguityHalfWidth<float, float>( options,
                                                                        std::sqrt( static_cast<double>( dimension ) ),
                                                                        std::sqrt( static_cast<double>( dimension ) ) );
    const double threshold = options.rankRelativeThreshold;
    const double robustFloor = std::min( 0.75, std::max( 0.05, 8.0 * threshold ) );
    if( threshold + 4.0 * outerFloatBand >= 1 )
    {
        throw std::invalid_argument( "rank-boundary ambiguity multiplier leaves no normalized far-above value" );
    }

    std::vector<double> normalizedEigenvalues( static_cast<std::size_t>( dimension ), 0 );
    for( int value = 0; value < robustValues; ++value )
    {
        const double fraction = robustValues == 1 ? 0 : static_cast<double>( value ) / ( robustValues - 1 );
        normalizedEigenvalues[static_cast<std::size_t>( value )] = std::exp( fraction * std::log( robustFloor ) );
    }

    const std::array<double, boundaryValues> boundarySpectrum{
        threshold + 4.0 * outerFloatBand,
        threshold + 0.5 * nominalFloatBand,
        threshold,
        std::max( 0.5 * threshold, threshold - 0.5 * nominalFloatBand ),
        std::max( 0.25 * threshold, threshold - 4.0 * outerFloatBand ),
        1e-4 * threshold,
        0,
        0,
    };
    std::copy( boundarySpectrum.begin(), boundarySpectrum.end(), normalizedEigenvalues.begin() + robustValues );

    SourceData source;
    source.predictors = matrixT<float>::Zero( options.samples, options.predictors );
    for( int diagonal = 0; diagonal < dimension; ++diagonal )
    {
        source.predictors( diagonal, diagonal ) =
            static_cast<float>( std::sqrt( normalizedEigenvalues[static_cast<std::size_t>( diagonal )] ) );
    }
    source.designedStructuralRank = dimension - 2;

    source.target.resize( options.samples );
    std::seed_seq seedSequence{ static_cast<std::uint32_t>( options.seed ),
                                static_cast<std::uint32_t>( options.seed >> 32 ),
                                static_cast<std::uint32_t>( innerDimension ) };
    std::mt19937 generator( seedSequence );
    for( int row = 0; row < options.samples; ++row )
    {
        source.target( row ) = 0.05f * nextUnitFloat( generator );
    }
    for( int diagonal = 0; diagonal < robustValues; ++diagonal )
    {
        const float coefficient = diagonal % 2 == 0 ? 0.30f : -0.20f;
        source.target( diagonal ) += coefficient * source.predictors( diagonal, diagonal );
    }
    return source;
}

/// Construct a dense FP32 source with exact dependencies and one near-dependent independent direction.
/** Four rows or columns are exact copies of independent source vectors, according to the selected Gram orientation.
 * One remaining independent row or column is scaled down by 1e-4. The FP64 SVD oracle verifies the intended exact
 * structural rank before any candidate's strict lambda-greater-than-zero rank is interpreted.
 */
SourceData makeRankZeroSourceData( const BenchmarkOptions &options /**< [in] dimensions and deterministic seed */ )
{
    constexpr int exactDependencies{ 4 };
    SourceData source = makeSourceData( options );
    const bool useTemporalGram = options.samples <= options.predictors;
    const int dimension = std::min( options.samples, options.predictors );
    const int independentValues = dimension - exactDependencies;
    if( useTemporalGram )
    {
        source.predictors.row( independentValues - 1 ) *= 1e-4f;
        for( int dependency = 0; dependency < exactDependencies; ++dependency )
        {
            source.predictors.row( independentValues + dependency ) = source.predictors.row( dependency );
        }
    }
    else
    {
        source.predictors.col( independentValues - 1 ) *= 1e-4f;
        for( int dependency = 0; dependency < exactDependencies; ++dependency )
        {
            source.predictors.col( independentValues + dependency ) = source.predictors.col( dependency );
        }
    }
    source.designedStructuralRank = independentValues;
    return source;
}

/// Construct the selected deterministic source family without entering a timed region.
SourceData makeSelectedSourceData( const BenchmarkOptions &options /**< [in] source-family controls */ )
{
    if( options.inputCase == InputCase::rankBoundary )
    {
        return makeRankBoundarySourceData( options );
    }
    if( options.inputCase == InputCase::rankZero )
    {
        return makeRankZeroSourceData( options );
    }
    return makeSourceData( options );
}

/// Hash the exact source arrays using a canonical versioned binary32 serialization.
std::string sourceValueSha256( const SourceData &source /**< [in] source arrays to identify */ )
{
    Sha256 digest;
    constexpr std::string_view encoding{ "hcireduce-p4-source-binary32-big-endian-column-major/v1" };
    digest.update( encoding.data(), encoding.size() );
    hashUnsigned( digest, static_cast<std::uint64_t>( source.predictors.rows() ) );
    hashUnsigned( digest, static_cast<std::uint64_t>( source.predictors.cols() ) );
    hashFloatValues( digest, source.predictors.data(), static_cast<std::size_t>( source.predictors.size() ) );
    hashUnsigned( digest, static_cast<std::uint64_t>( source.target.rows() ) );
    hashFloatValues( digest, source.target.data(), static_cast<std::size_t>( source.target.size() ) );
    return digest.finish();
}

/// Verify immutable files and construct the identities shared by every emitted row.
RunIdentity makeRunIdentity( const BenchmarkOptions &options, /**< [in] declared contract and case controls */
                             const SourceData &source,        /**< [in] generated source values */
                             std::string_view executable )    /**< [in] argv-zero fallback executable path */
{
    RunIdentity identity;
    identity.acceptanceContractPath = options.acceptanceContractPath;
    identity.acceptanceContractSha256 = sha256File( options.acceptanceContractPath );
    if( identity.acceptanceContractSha256 != acceptanceSha256 )
    {
        throw std::runtime_error( "acceptance contract content does not match compiled SHA-256 " +
                                  std::string( acceptanceSha256 ) );
    }
    identity.caseSetId = options.inputCase == InputCase::wellConditioned
                             ? "p4_well_separated"
                             : ( options.inputCase == InputCase::rankBoundary ? "rank_nonzero_boundary"
                                                                              : "rank_zero_exact_dependencies" );
    identity.caseId = options.caseId;
    identity.sourceValueSha256 = sourceValueSha256( source );
    try
    {
        identity.executableSha256 = sha256File( "/proc/self/exe" );
    }
    catch( const std::exception & )
    {
        identity.executableSha256 = sha256File( std::string( executable ) );
    }
    return identity;
}

/// Construct the direct FP64 thin-SVD oracle required by the engineering smoke gates.
NumericalOracle makeNumericalOracle( const SourceData &source, /**< [in] exact original FP32 source values */
                                     const BenchmarkOptions &options /**< [in] Gram orientation and retained modes */ )
{
    const Eigen::MatrixXd predictors = source.predictors.matrix().cast<double>();
    const Eigen::VectorXd target = source.target.matrix().cast<double>();
    Eigen::JacobiSVD<Eigen::MatrixXd> decomposition( predictors, Eigen::ComputeThinU | Eigen::ComputeThinV );
    if( decomposition.info() != Eigen::Success )
    {
        throw std::runtime_error( "independent FP64 JacobiSVD numerical oracle failed" );
    }

    NumericalOracle oracle;
    if( options.samples <= options.predictors )
    {
        oracle.gram = predictors * predictors.transpose();
        oracle.gramModeColumns = decomposition.matrixU().leftCols( options.modes );
    }
    else
    {
        oracle.gram = predictors.transpose() * predictors;
        oracle.gramModeColumns = decomposition.matrixV().leftCols( options.modes );
    }
    oracle.eigenvalues = decomposition.singularValues().head( options.modes ).array().square();
    const Eigen::MatrixXd leftModes = decomposition.matrixU().leftCols( options.modes );
    oracle.residual = target - leftModes * ( leftModes.transpose() * target );
    oracle.targetNorm = target.norm();
    oracle.residualNorm = oracle.residual.norm();
    return oracle;
}

/// Execute the uncentered P4 Gram/eigensolve/projection microkernel for one precision policy.
/** The only mxlib numerical call is the production `eigenSYEVR` path. For the mixed policy, conversion to and from
 * FP64 is included in the eigensolve stage so every other timed operation remains FP32. This intentionally stops
 * short of the complete direct-path contract: rank selection, output-plane materialization, coefficients, probes,
 * and input validation are not measured here.
 */
template <typename calculationT, typename eigensolverT>
StageTimings
executeCandidate( std::string_view candidate,      /**< [in] stable precision-policy identifier */
                  const SourceData &source,        /**< [in] common FP32 source values */
                  const BenchmarkOptions &options, /**< [in] dimensions and retained-mode count */
                  CandidateWorkspace<calculationT, eigensolverT> &state /**< [in,out] reusable candidate storage */ )
{
    const bool useTemporalGram = options.samples <= options.predictors;
    const int gramDimension = std::min( options.samples, options.predictors );
    const auto totalStart = benchmarkClockT::now();

    const auto conversionStart = benchmarkClockT::now();
    state.predictors = source.predictors.template cast<calculationT>();
    state.target = source.target.template cast<calculationT>();
    const double conversionSeconds = secondsSince( conversionStart );

    const auto gramStart = benchmarkClockT::now();
    if( useTemporalGram )
    {
        state.gram = ( state.predictors.matrix() * state.predictors.matrix().transpose() ).array();
        state.crossProduct.resize( 0 );
    }
    else
    {
        state.gram = ( state.predictors.matrix().transpose() * state.predictors.matrix() ).array();
        state.crossProduct = ( state.predictors.matrix().transpose() * state.target.matrix() ).array();
    }
    const double gramSeconds = secondsSince( gramStart );

    const auto eigensolveStart = benchmarkClockT::now();
    MXLAPACK_INT solverStatus{ 0 };
    if constexpr( std::is_same_v<calculationT, eigensolverT> )
    {
        solverStatus = mx::math::eigenSYEVR( state.eigenvectors,
                                             state.eigenvalues,
                                             state.gram,
                                             0,
                                             gramDimension,
                                             'L',
                                             &state.solverWorkspace );
    }
    else
    {
        state.solverInput = state.gram.template cast<eigensolverT>();
        solverStatus = mx::math::eigenSYEVR( state.solverEigenvectors,
                                             state.solverEigenvalues,
                                             state.solverInput,
                                             0,
                                             gramDimension,
                                             'L',
                                             &state.solverWorkspace );
    }
    if( solverStatus != 0 )
    {
        throw std::runtime_error( std::string( candidate ) + " eigenSYEVR failed with status " +
                                  std::to_string( solverStatus ) );
    }
    if constexpr( !std::is_same_v<calculationT, eigensolverT> )
    {
        state.eigenvectors = state.solverEigenvectors.template cast<calculationT>();
        state.eigenvalues = state.solverEigenvalues.block( 0, 0, gramDimension, 1 ).template cast<calculationT>();
    }
    const double eigensolveSeconds = secondsSince( eigensolveStart );

    const auto projectionStart = benchmarkClockT::now();
    state.residual = state.target;
    for( int retainedMode = 1; retainedMode <= options.modes; ++retainedMode )
    {
        const int eigenIndex = gramDimension - retainedMode;
        const calculationT eigenvalue = state.eigenvalues( eigenIndex );
        if( !std::isfinite( static_cast<double>( eigenvalue ) ) || eigenvalue <= calculationT{ 0 } )
        {
            throw std::runtime_error( std::string( candidate ) + " produced a nonpositive retained eigenvalue" );
        }

        const vectorT<calculationT> eigenvector = state.eigenvectors.col( eigenIndex );
        if( useTemporalGram )
        {
            const calculationT coefficient = eigenvector.matrix().dot( state.target.matrix() );
            state.projectedMode = eigenvector * coefficient;
        }
        else
        {
            const calculationT coefficient = eigenvector.matrix().dot( state.crossProduct.matrix() ) / eigenvalue;
            state.projectedMode = ( state.predictors.matrix() * eigenvector.matrix() ).array();
            state.projectedMode *= coefficient;
        }
        state.residual -= state.projectedMode;
    }
    const double projectionSeconds = secondsSince( projectionStart );
    const double totalSeconds = secondsSince( totalStart );

    return { conversionSeconds, gramSeconds, eigensolveSeconds, projectionSeconds, totalSeconds };
}

/// Capture an execution's outputs and independently rebuild its pre-SYEVR Gram outside timed regions.
template <typename calculationT, typename eigensolverT>
NumericalResult captureNumericalResult(
    std::string_view candidate,      /**< [in] stable precision-policy identifier */
    const SourceData &source,        /**< [in] common FP32 source values */
    const BenchmarkOptions &options, /**< [in] dimensions and Gram orientation */
    const CandidateWorkspace<calculationT, eigensolverT> &state /**< [in] completed validation execution */ )
{
    const matrixT<calculationT> predictors = source.predictors.template cast<calculationT>();
    matrixT<calculationT> validationGram;
    if( options.samples <= options.predictors )
    {
        validationGram = ( predictors.matrix() * predictors.matrix().transpose() ).array();
    }
    else
    {
        validationGram = ( predictors.matrix().transpose() * predictors.matrix() ).array();
    }

    NumericalResult result;
    result.candidate = candidate;
    result.gram = validationGram.template cast<double>();
    result.eigenvectors = state.eigenvectors.template cast<double>();
    result.eigenvalues = state.eigenvalues.block( 0, 0, state.eigenvalues.rows(), 1 ).template cast<double>();
    result.residual = state.residual.template cast<double>();
    return result;
}

/// Return a retained eigenvector block ordered from the largest eigenvalue downward.
Eigen::MatrixXd retainedEigenvectors( const NumericalResult &result, /**< [in] complete ascending eigensystem */
                                      int modes )                    /**< [in] number of leading vectors to return */
{
    const Eigen::Index dimension = result.eigenvectors.cols();
    Eigen::MatrixXd retained( result.eigenvectors.rows(), modes );
    for( int mode = 0; mode < modes; ++mode )
    {
        retained.col( mode ) = result.eigenvectors.matrix().col( dimension - 1 - mode );
    }
    return retained;
}

/// Return retained eigenvalues ordered from largest to smallest.
Eigen::VectorXd retainedEigenvalues( const NumericalResult &result, /**< [in] complete ascending eigensystem */
                                     int modes )                    /**< [in] number of leading values to return */
{
    const Eigen::Index dimension = result.eigenvalues.rows();
    Eigen::VectorXd retained( modes );
    for( int mode = 0; mode < modes; ++mode )
    {
        retained( mode ) = result.eigenvalues( dimension - 1 - mode );
    }
    return retained;
}

/// Capture a norm ratio together with the exact numerator and protected denominator used.
template <typename derivedT>
NormalizedMetric normalizedMetric( const Eigen::MatrixBase<derivedT> &difference, /**< [in] numerator values */
                                   double referenceNorm ) /**< [in] already-evaluated denominator norm */
{
    NormalizedMetric metric;
    metric.numerator = difference.norm();
    metric.denominator = std::max( referenceNorm, std::numeric_limits<double>::min() );
    metric.value = metric.numerator / metric.denominator;
    return metric;
}

/// Calculate the largest principal angle between two equal-width subspaces.
double maximumPrincipalAngle( const Eigen::MatrixXd &candidateModes, /**< [in] candidate orthonormal columns */
                              const Eigen::MatrixXd &oracleModes )   /**< [in] oracle orthonormal columns */
{
    Eigen::JacobiSVD<Eigen::MatrixXd> decomposition( candidateModes.transpose() * oracleModes );
    if( decomposition.info() != Eigen::Success || decomposition.singularValues().rows() == 0 )
    {
        throw std::runtime_error( "principal-angle SVD failed" );
    }
    const double minimumCosine = std::clamp( decomposition.singularValues().minCoeff(), 0.0, 1.0 );
    return std::acos( minimumCosine );
}

/// Find the normalized eigenvalue ratios immediately below and above a strict relative threshold.
RankBoundary characterizeRankBoundary( const Eigen::VectorXd &normalizedEigenvalues, /**< [in] ascending ratios */
                                       double threshold /**< [in] strict relative retain threshold */ )
{
    RankBoundary boundary;
    for( Eigen::Index index = 0; index < normalizedEigenvalues.rows(); ++index )
    {
        const double ratio = normalizedEigenvalues( index );
        if( ratio <= threshold )
        {
            boundary.ratioBelow = ratio;
        }
        else if( !std::isfinite( boundary.ratioAbove ) )
        {
            boundary.ratioAbove = ratio;
        }
    }
    if( std::isfinite( boundary.ratioBelow ) && std::isfinite( boundary.ratioAbove ) )
    {
        boundary.eigengap = boundary.ratioAbove - boundary.ratioBelow;
    }
    return boundary;
}

/// Build an independent FP64 singular-value oracle for a bounded rank-sensitive input.
RankOracle makeRankOracle( const SourceData &source, /**< [in] common FP32 source values */
                           const BenchmarkOptions &options /**< [in] rank threshold and dimensions */ )
{
    const Eigen::MatrixXd predictors = source.predictors.cast<double>().matrix();
    Eigen::JacobiSVD<Eigen::MatrixXd> decomposition( predictors );
    if( decomposition.info() != Eigen::Success )
    {
        throw std::runtime_error( "independent FP64 JacobiSVD rank oracle failed" );
    }

    const Eigen::VectorXd singularValues = decomposition.singularValues();
    if( singularValues.rows() == 0 || !singularValues.allFinite() || singularValues( 0 ) <= 0 )
    {
        throw std::runtime_error( "independent FP64 JacobiSVD rank oracle returned an invalid spectrum" );
    }

    RankOracle oracle;
    const double largestSquared = singularValues( 0 ) * singularValues( 0 );
    oracle.lambdaMaximum = largestSquared;
    const double zeroThreshold =
        std::max( options.samples, options.predictors ) * std::numeric_limits<double>::epsilon() * singularValues( 0 );
    oracle.observedStructuralRank = static_cast<int>( ( singularValues.array() > zeroThreshold ).count() );
    oracle.designedStructuralRank = source.designedStructuralRank;
    oracle.exactNullity = static_cast<int>( singularValues.rows() ) - oracle.observedStructuralRank;
    if( oracle.observedStructuralRank <= 0 )
    {
        throw std::runtime_error( "independent FP64 oracle found no structurally positive singular value" );
    }

    oracle.normalizedEigenvalues.resize( singularValues.rows() );
    for( Eigen::Index ascending = 0; ascending < singularValues.rows(); ++ascending )
    {
        const double singularValue = singularValues( singularValues.rows() - 1 - ascending );
        const double rawRatio = singularValue * singularValue / largestSquared;
        if( ascending < oracle.exactNullity )
        {
            oracle.rawLargestNullRatio = std::max( oracle.rawLargestNullRatio, rawRatio );
            oracle.normalizedEigenvalues( ascending ) = 0;
        }
        else
        {
            oracle.normalizedEigenvalues( ascending ) = rawRatio;
        }
    }
    oracle.smallestPositiveRatio = oracle.normalizedEigenvalues( oracle.exactNullity );
    oracle.gramFrobeniusRatio = oracle.normalizedEigenvalues.norm();
    const Eigen::MatrixXd absolutePredictors = predictors.cwiseAbs();
    Eigen::MatrixXd absoluteGram;
    if( options.samples <= options.predictors )
    {
        absoluteGram = absolutePredictors * absolutePredictors.transpose();
    }
    else
    {
        absoluteGram = absolutePredictors.transpose() * absolutePredictors;
    }
    oracle.absoluteGramFrobeniusRatio = absoluteGram.norm() / largestSquared;
    oracle.minimumThresholdDistance =
        ( oracle.normalizedEigenvalues.array() - options.rankRelativeThreshold ).abs().minCoeff();
    oracle.thresholdRank =
        options.inputCase == InputCase::rankZero
            ? oracle.observedStructuralRank
            : static_cast<int>( ( oracle.normalizedEigenvalues.array() > options.rankRelativeThreshold ).count() );
    oracle.boundary = characterizeRankBoundary( oracle.normalizedEigenvalues, options.rankRelativeThreshold );
    if( oracle.observedStructuralRank != oracle.designedStructuralRank )
    {
        throw std::runtime_error( "independent FP64 oracle disagrees with the designed exact structural rank" );
    }
    return oracle;
}

/// Serialize an ascending normalized spectrum in one-based descending order.
std::string serializeSpectrumDescending( const Eigen::VectorXd &values /**< [in] ascending normalized spectrum */ )
{
    std::ostringstream output;
    output << std::setprecision( 17 );
    for( Eigen::Index descending = 0; descending < values.rows(); ++descending )
    {
        if( descending != 0 )
        {
            output << ';';
        }
        output << values( values.rows() - 1 - descending );
    }
    return output.str();
}

/// Serialize strict threshold decisions for an ascending normalized spectrum in descending order.
std::string serializeDecisionsDescending( const Eigen::VectorXd &values, /**< [in] ascending normalized spectrum */
                                          double threshold )             /**< [in] strict threshold */
{
    std::ostringstream output;
    for( Eigen::Index descending = 0; descending < values.rows(); ++descending )
    {
        if( descending != 0 )
        {
            output << ';';
        }
        output << ( values( values.rows() - 1 - descending ) > threshold ? '1' : '0' );
    }
    return output.str();
}

/// Compare one candidate's strict rank decisions with the FP64 oracle outside its precision-scaled ambiguity band.
template <typename calculationT, typename eigensolverT>
RankMetrics validateRank( const NumericalResult &candidate, /**< [in] candidate ascending eigensystem */
                          const RankOracle &oracle,         /**< [in] independent FP64 SVD spectrum */
                          const BenchmarkOptions &options /**< [in] threshold and ambiguity controls */ )
{
    if( candidate.eigenvalues.rows() != oracle.normalizedEigenvalues.rows() )
    {
        throw std::runtime_error( candidate.candidate + " rank spectrum has the wrong dimension" );
    }
    if( !candidate.eigenvalues.matrix().allFinite() )
    {
        throw std::runtime_error( candidate.candidate + " rank spectrum contains a nonfinite eigenvalue" );
    }
    for( Eigen::Index index = 1; index < candidate.eigenvalues.rows(); ++index )
    {
        if( candidate.eigenvalues( index ) < candidate.eigenvalues( index - 1 ) )
        {
            throw std::runtime_error( candidate.candidate + " rank spectrum is not ascending" );
        }
    }
    const double largestEigenvalue = candidate.eigenvalues( candidate.eigenvalues.rows() - 1 );
    if( !std::isfinite( largestEigenvalue ) || largestEigenvalue <= 0 )
    {
        throw std::runtime_error( candidate.candidate + " rank spectrum has no positive maximum" );
    }

    const Eigen::VectorXd candidateRatios = candidate.eigenvalues.matrix() / largestEigenvalue;
    const long double candidateThreshold =
        static_cast<long double>( options.rankRelativeThreshold ) * static_cast<long double>( largestEigenvalue );
    RankMetrics metrics;
    metrics.evaluated = true;
    metrics.oracleRank = oracle.thresholdRank;
    metrics.candidateRank = 0;
    for( Eigen::Index index = 0; index < candidate.eigenvalues.rows(); ++index )
    {
        if( static_cast<long double>( candidate.eigenvalues( index ) ) > candidateThreshold )
        {
            ++metrics.candidateRank;
        }
    }
    metrics.oracleStructuralRank = oracle.observedStructuralRank;
    metrics.designedStructuralRank = oracle.designedStructuralRank;
    metrics.oracleExactNullity = oracle.exactNullity;
    metrics.oracleSmallestPositiveRatio = oracle.smallestPositiveRatio;
    metrics.oracleRawLargestNullRatio = oracle.rawLargestNullRatio;
    metrics.oracleGramFrobeniusRatio = oracle.gramFrobeniusRatio;
    metrics.oracleAbsoluteGramRatio = oracle.absoluteGramFrobeniusRatio;
    metrics.oracleMinimumThresholdDistance = oracle.minimumThresholdDistance;
    const int gramDimension = std::min( options.samples, options.predictors );
    const int innerDimension = options.samples <= options.predictors ? options.predictors : options.samples;
    const int calculationOperations = options.inputCase == InputCase::rankZero ? innerDimension : 1;
    metrics.calculationGamma = gammaN<calculationT>( calculationOperations );
    metrics.eigensolverGamma = gammaN<eigensolverT>( gramDimension );
    metrics.calculationBandComponent = metrics.calculationGamma * oracle.absoluteGramFrobeniusRatio;
    metrics.eigensolverBandComponent = metrics.eigensolverGamma * oracle.gramFrobeniusRatio;
    metrics.oracleLambdaMaximum = oracle.lambdaMaximum;
    metrics.ambiguityHalfWidth = rankAmbiguityHalfWidth<calculationT, eigensolverT>( options,
                                                                                     oracle.absoluteGramFrobeniusRatio,
                                                                                     oracle.gramFrobeniusRatio );
    metrics.oracleBoundary = oracle.boundary;
    metrics.candidateBoundary = characterizeRankBoundary( candidateRatios, options.rankRelativeThreshold );
    metrics.rankDiffersFromOracle = metrics.candidateRank != metrics.oracleRank;
    metrics.oracleSpectrumDescending = serializeSpectrumDescending( oracle.normalizedEigenvalues );
    metrics.candidateSpectrumDescending = serializeSpectrumDescending( candidateRatios );
    metrics.oracleDecisionsDescending =
        serializeDecisionsDescending( oracle.normalizedEigenvalues, options.rankRelativeThreshold );
    metrics.candidateDecisionsDescending =
        serializeDecisionsDescending( candidateRatios, options.rankRelativeThreshold );

    for( Eigen::Index index = 0; index < oracle.normalizedEigenvalues.rows(); ++index )
    {
        const double oracleRatio = oracle.normalizedEigenvalues( index );
        if( std::abs( oracleRatio - options.rankRelativeThreshold ) <= metrics.ambiguityHalfWidth )
        {
            ++metrics.ambiguousEigenvalues;
            const Eigen::Index descendingIndex = oracle.normalizedEigenvalues.rows() - index;
            if( !metrics.ambiguousIndicesDescending.empty() )
            {
                metrics.ambiguousIndicesDescending += ';';
            }
            metrics.ambiguousIndicesDescending += std::to_string( descendingIndex );
            continue;
        }
        const bool oracleRetained = oracleRatio > options.rankRelativeThreshold;
        const bool candidateRetained = static_cast<long double>( candidate.eigenvalues( index ) ) > candidateThreshold;
        metrics.agreementOutsideAmbiguity = metrics.agreementOutsideAmbiguity && oracleRetained == candidateRetained;
    }
    return metrics;
}

/// Calculate invariant, eigensolver, and final-residual validation metrics.
ValidationMetrics validateResult( const NumericalResult &candidate, /**< [in] result to evaluate */
                                  const NumericalOracle &oracle,    /**< [in] direct FP64 thin-SVD oracle */
                                  const BenchmarkOptions &options,  /**< [in] retained modes and tolerances */
                                  const RankMetrics &rankMetrics /**< [in] optional ambiguity-aware rank result */ )
{
    ValidationMetrics metrics;
    const Eigen::MatrixXd candidateGram = candidate.gram.matrix();
    const Eigen::MatrixXd gramDifference = candidateGram - oracle.gram;
    metrics.gramRelativeError = normalizedMetric( gramDifference, oracle.gram.norm() );
    metrics.gramMaximumError = gramDifference.cwiseAbs().maxCoeff();

    const Eigen::MatrixXd candidateVectors = retainedEigenvectors( candidate, options.modes );
    const Eigen::VectorXd candidateValues = retainedEigenvalues( candidate, options.modes );
    metrics.eigenvalueRelativeError =
        normalizedMetric( candidateValues - oracle.eigenvalues, oracle.eigenvalues.norm() );

    const Eigen::MatrixXd candidateProjector = candidateVectors * candidateVectors.transpose();
    const Eigen::MatrixXd oracleProjector = oracle.gramModeColumns * oracle.gramModeColumns.transpose();
    metrics.projectorRelativeError = normalizedMetric( candidateProjector - oracleProjector, oracleProjector.norm() );
    metrics.maximumPrincipalAngle = maximumPrincipalAngle( candidateVectors, oracle.gramModeColumns );

    const Eigen::VectorXd residualDifference = candidate.residual.matrix() - oracle.residual;
    metrics.residualOverOracleResidual = normalizedMetric( residualDifference, oracle.residualNorm );
    metrics.residualOverTarget = normalizedMetric( residualDifference, oracle.targetNorm );
    metrics.residualMaximumError = residualDifference.cwiseAbs().maxCoeff();

    const Eigen::MatrixXd solverDifference =
        candidateGram * candidateVectors - candidateVectors * candidateValues.asDiagonal();
    const double solverScale = candidateGram.norm() * candidateVectors.norm();
    metrics.solverResidual = normalizedMetric( solverDifference, solverScale );
    const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity( options.modes, options.modes );
    metrics.orthogonalityDefect = ( candidateVectors.transpose() * candidateVectors - identity ).norm();
    metrics.rank = rankMetrics;

    metrics.passed = metrics.gramRelativeError.value <= options.gramRelativeTolerance &&
                     metrics.eigenvalueRelativeError.value <= options.eigenvalueRelativeTolerance &&
                     metrics.projectorRelativeError.value <= options.projectorRelativeTolerance &&
                     metrics.residualOverOracleResidual.value <= options.residualRelativeTolerance &&
                     metrics.residualOverTarget.value <= options.residualRelativeTolerance &&
                     metrics.solverResidual.value <= options.solverResidualTolerance &&
                     metrics.orthogonalityDefect <= options.orthogonalityTolerance &&
                     ( !metrics.rank.evaluated || metrics.rank.agreementOutsideAmbiguity );
    return metrics;
}

/// Format a finite floating-point CSV field or the literal `null`.
std::string csvNumber( double value /**< [in] value, with NaN representing null */ )
{
    if( std::isnan( value ) )
    {
        return "null";
    }
    std::ostringstream output;
    output << std::setprecision( 17 ) << value;
    return output.str();
}

/// Quote one CSV field when RFC-4180 delimiters are present.
std::string csvField( std::string_view value /**< [in] unescaped field text */ )
{
    if( value.find_first_of( ",\"\r\n" ) == std::string_view::npos )
    {
        return std::string( value );
    }
    std::string escaped{ "\"" };
    for( const char character : value )
    {
        if( character == '\"' )
        {
            escaped += '\"';
        }
        escaped += character;
    }
    escaped += '\"';
    return escaped;
}

/// Emit a complete CSV row with consistent escaping.
void emitFields( std::ostream &output,                    /**< [out] machine-readable destination */
                 const std::vector<std::string> &fields ) /**< [in] ordered unescaped field values */
{
    for( std::size_t field = 0; field < fields.size(); ++field )
    {
        if( field != 0 )
        {
            output << ',';
        }
        output << csvField( fields[field] );
    }
    output << '\n';
}

/// Print the versioned long-form CSV schema used for contract-bound metric rows.
void emitCsvHeader( std::ostream &output /**< [out] machine-readable destination */ )
{
    output << "record_schema,record,acceptance_schema,acceptance_contract_path,acceptance_contract_sha256,"
              "case_set_id,case_id,gate_id,candidate_id,oracle_id,source_value_sha256,executable_sha256,algorithm,"
              "repetition,execution_position,outer_worker_count,inner_blas_lapack_threads,seed,warmups,"
              "total_repetitions,input_case,scenario,samples,predictors,targets,pixels,references,gram_orientation,"
              "gram_dimension,inner_product_dimension,requested_modes,realized_modes,calculation_precision,"
              "eigensolver_precision,metric_name,numerator_name,numerator_value,numerator_not_applicable_reason,"
              "denominator_name,denominator_value,denominator_not_applicable_reason,value,units,aggregation_scope,"
              "limit_operator,limit_value,applicable,not_applicable_reason,passed,gate_status,stage,diagnostic_value,"
              "diagnostic_units,benchmark_scope,rank_threshold,strict_comparison,fp32_epsilon,fp64_epsilon,"
              "rank_band_multiplier,calculation_gamma,eigensolver_gamma,calculation_band_component,"
              "eigensolver_band_component,oracle_lambda_max,oracle_threshold_rank,candidate_threshold_rank,"
              "designed_structural_rank,oracle_structural_rank,exact_nullity,oracle_smallest_positive_ratio,"
              "oracle_raw_largest_null_ratio,oracle_gram_frobenius_ratio,oracle_absolute_gram_frobenius_ratio,"
              "rank_ambiguity_half_width,ambiguous_indices_descending,ambiguous_count,"
              "oracle_ratio_below_threshold,oracle_ratio_above_threshold,candidate_ratio_below_threshold,"
              "candidate_ratio_above_threshold,oracle_boundary_eigengap,candidate_boundary_eigengap,"
              "oracle_normalized_spectrum_descending,candidate_normalized_spectrum_descending,"
              "oracle_decisions_descending,candidate_decisions_descending,requested_mode_cluster_membership,"
              "outside_band_agreement,affected_modes,product_sensitivity_record,ambiguity_resolution,"
              "promotion_ready,gram_triangle,source_value_encoding\n";
}

/// Reject internally inconsistent null, applicability, denominator, and pass-state representations.
void validateMetricRecord( const MetricRecord &metric /**< [in] record about to be emitted */ )
{
    const bool numeratorNull = std::isnan( metric.numeratorValue );
    const bool denominatorNull = std::isnan( metric.denominatorValue );
    if( metric.numeratorName.empty() != numeratorNull ||
        ( numeratorNull && metric.numeratorNotApplicableReason.empty() ) )
    {
        throw std::logic_error( "metric numerator representation is inconsistent for " + metric.metricName );
    }
    if( !numeratorNull && !std::isfinite( metric.numeratorValue ) )
    {
        throw std::logic_error( "metric numerator is nonfinite for " + metric.metricName );
    }
    if( metric.denominatorName.empty() != denominatorNull ||
        ( denominatorNull && metric.denominatorNotApplicableReason.empty() ) )
    {
        throw std::logic_error( "metric denominator representation is inconsistent for " + metric.metricName );
    }
    if( !denominatorNull && ( !std::isfinite( metric.denominatorValue ) || metric.denominatorValue <= 0 ) )
    {
        throw std::logic_error( "metric denominator is not finite and positive for " + metric.metricName );
    }
    if( !metric.applicable )
    {
        if( !std::isnan( metric.value ) || metric.passed >= 0 || metric.notApplicableReason.empty() )
        {
            throw std::logic_error( "unavailable metric representation is inconsistent for " + metric.metricName );
        }
        return;
    }
    if( !std::isfinite( metric.value ) )
    {
        throw std::logic_error( "applicable metric is nonfinite for " + metric.metricName );
    }
    if( metric.gateStatus == "frozen-report-only" ? metric.passed >= 0 : metric.passed < 0 )
    {
        throw std::logic_error( "metric pass-state representation is inconsistent for " + metric.metricName );
    }
}

/// Emit one acceptance-bound long-form metric record.
void emitMetricRecord( std::ostream &output,            /**< [out] machine-readable destination */
                       const RunIdentity &identity,     /**< [in] verified run identities */
                       std::string_view candidate,      /**< [in] precision-policy identifier */
                       const BenchmarkOptions &options, /**< [in] dimensions and declared controls */
                       const RankMetrics &rank,         /**< [in] optional rank evidence */
                       const MetricRecord &metric,      /**< [in] gate-specific metric evidence */
                       int repetition,                  /**< [in] zero-based execution repetition */
                       int executionPosition )          /**< [in] balanced-order position, or negative */
{
    validateMetricRecord( metric );
    const bool temporalGram = options.samples <= options.predictors;
    const std::string_view calculationPrecision = candidate == "P4-D64" ? "fp64" : "fp32";
    const std::string_view eigensolverPrecision = candidate == "P4-F32" ? "fp32" : "fp64";
    std::string_view benchmarkScope{ "well_conditioned_uncentered_in_sample_microbenchmark" };
    if( options.inputCase == InputCase::rankBoundary )
    {
        benchmarkScope = "rank_boundary_uncentered_in_sample_microbenchmark";
    }
    else if( options.inputCase == InputCase::rankZero )
    {
        benchmarkScope = "rank_zero_uncentered_in_sample_microbenchmark";
    }
    const auto integerOrNull = []( int value ) { return value < 0 ? std::string( "null" ) : std::to_string( value ); };
    const auto booleanOrNull = []( int value )
    { return value < 0 ? std::string( "null" ) : std::string( value == 0 ? "false" : "true" ); };
    const std::string nullReason = metric.applicable ? "null" : metric.notApplicableReason;
    const std::string numeratorName = metric.numeratorName.empty() ? "null" : metric.numeratorName;
    const std::string denominatorName = metric.denominatorName.empty() ? "null" : metric.denominatorName;
    const std::string rankValue = rank.evaluated ? csvNumber( options.rankRelativeThreshold ) : "null";
    const std::string ambiguityResolution =
        !rank.evaluated ? "not_applicable" : ( rank.ambiguousEigenvalues == 0 ? "not_required" : "unresolved" );
    const std::string affectedModes =
        !rank.evaluated ? "not_applicable" : "none_requested_modes_outside_diagnostic_tail";
    const std::string productSensitivity = !rank.evaluated ? "not_applicable" : "not_measured_by_kernel_microbenchmark";

    emitFields( output,
                { std::string( csvSchema ),
                  metric.record,
                  std::string( acceptanceSchema ),
                  identity.acceptanceContractPath,
                  identity.acceptanceContractSha256,
                  identity.caseSetId,
                  identity.caseId,
                  metric.gateId,
                  std::string( candidate ),
                  metric.oracleId,
                  identity.sourceValueSha256,
                  identity.executableSha256,
                  "p4",
                  std::to_string( repetition ),
                  std::to_string( executionPosition ),
                  "1",
                  std::to_string( options.innerBlasLapackThreads ),
                  std::to_string( options.seed ),
                  std::to_string( options.warmups ),
                  std::to_string( options.repetitions ),
                  std::string( inputCaseName( options.inputCase ) ),
                  "direct_uncentered_in_sample",
                  std::to_string( options.samples ),
                  std::to_string( options.predictors ),
                  "null",
                  "null",
                  "null",
                  temporalGram ? "temporal_T_le_K" : "predictor_T_gt_K",
                  std::to_string( std::min( options.samples, options.predictors ) ),
                  std::to_string( temporalGram ? options.predictors : options.samples ),
                  std::to_string( options.modes ),
                  std::to_string( options.modes ),
                  std::string( calculationPrecision ),
                  std::string( eigensolverPrecision ),
                  metric.metricName,
                  numeratorName,
                  csvNumber( metric.numeratorValue ),
                  metric.numeratorNotApplicableReason.empty() ? "null" : metric.numeratorNotApplicableReason,
                  denominatorName,
                  csvNumber( metric.denominatorValue ),
                  metric.denominatorNotApplicableReason.empty() ? "null" : metric.denominatorNotApplicableReason,
                  csvNumber( metric.value ),
                  metric.units,
                  metric.aggregationScope,
                  metric.limitOperator.empty() ? "null" : metric.limitOperator,
                  csvNumber( metric.limitValue ),
                  metric.applicable ? "true" : "false",
                  nullReason,
                  booleanOrNull( metric.passed ),
                  metric.gateStatus,
                  metric.stage,
                  csvNumber( metric.diagnosticValue ),
                  metric.diagnosticUnits.empty() ? "null" : metric.diagnosticUnits,
                  std::string( benchmarkScope ),
                  rankValue,
                  rank.evaluated ? "lambda_i>tau*lambda_1" : "null",
                  rank.evaluated ? csvNumber( std::numeric_limits<float>::epsilon() ) : "null",
                  rank.evaluated ? csvNumber( std::numeric_limits<double>::epsilon() ) : "null",
                  rank.evaluated ? csvNumber( options.rankAmbiguityMultiplier ) : "null",
                  rank.evaluated ? csvNumber( rank.calculationGamma ) : "null",
                  rank.evaluated ? csvNumber( rank.eigensolverGamma ) : "null",
                  rank.evaluated ? csvNumber( rank.calculationBandComponent ) : "null",
                  rank.evaluated ? csvNumber( rank.eigensolverBandComponent ) : "null",
                  rank.evaluated ? csvNumber( rank.oracleLambdaMaximum ) : "null",
                  rank.evaluated ? integerOrNull( rank.oracleRank ) : "null",
                  rank.evaluated ? integerOrNull( rank.candidateRank ) : "null",
                  rank.evaluated ? integerOrNull( rank.designedStructuralRank ) : "null",
                  rank.evaluated ? integerOrNull( rank.oracleStructuralRank ) : "null",
                  rank.evaluated ? integerOrNull( rank.oracleExactNullity ) : "null",
                  rank.evaluated ? csvNumber( rank.oracleSmallestPositiveRatio ) : "null",
                  rank.evaluated ? csvNumber( rank.oracleRawLargestNullRatio ) : "null",
                  rank.evaluated ? csvNumber( rank.oracleGramFrobeniusRatio ) : "null",
                  rank.evaluated ? csvNumber( rank.oracleAbsoluteGramRatio ) : "null",
                  rank.evaluated ? csvNumber( rank.ambiguityHalfWidth ) : "null",
                  rank.evaluated && !rank.ambiguousIndicesDescending.empty() ? rank.ambiguousIndicesDescending
                                                                             : ( rank.evaluated ? "none" : "null" ),
                  rank.evaluated ? std::to_string( rank.ambiguousEigenvalues ) : "null",
                  rank.evaluated ? csvNumber( rank.oracleBoundary.ratioBelow ) : "null",
                  rank.evaluated ? csvNumber( rank.oracleBoundary.ratioAbove ) : "null",
                  rank.evaluated ? csvNumber( rank.candidateBoundary.ratioBelow ) : "null",
                  rank.evaluated ? csvNumber( rank.candidateBoundary.ratioAbove ) : "null",
                  rank.evaluated ? csvNumber( rank.oracleBoundary.eigengap ) : "null",
                  rank.evaluated ? csvNumber( rank.candidateBoundary.eigengap ) : "null",
                  rank.evaluated ? rank.oracleSpectrumDescending : "null",
                  rank.evaluated ? rank.candidateSpectrumDescending : "null",
                  rank.evaluated ? rank.oracleDecisionsDescending : "null",
                  rank.evaluated ? rank.candidateDecisionsDescending : "null",
                  rank.evaluated ? "none_no_requested_mode_cluster_case" : "not_applicable",
                  rank.evaluated ? ( rank.agreementOutsideAmbiguity ? "true" : "false" ) : "null",
                  affectedModes,
                  productSensitivity,
                  ambiguityResolution,
                  "false",
                  "full_symmetric",
                  "ieee754_binary32_big_endian_column_major_v1" } );
}

/// Construct one applicable normalized gate record.
MetricRecord normalizedGate( std::string gateId,             /**< [in] frozen gate identifier */
                             std::string oracleId,           /**< [in] contract oracle identifier */
                             std::string metricName,         /**< [in] unambiguous metric name */
                             std::string numeratorName,      /**< [in] formula numerator identity */
                             std::string denominatorName,    /**< [in] protected denominator identity */
                             const NormalizedMetric &metric, /**< [in] calculated metric evidence */
                             double limit )                  /**< [in] frozen upper limit */
{
    MetricRecord record;
    record.gateId = std::move( gateId );
    record.oracleId = std::move( oracleId );
    record.metricName = std::move( metricName );
    record.numeratorName = std::move( numeratorName );
    record.numeratorValue = metric.numerator;
    record.denominatorName = std::move( denominatorName );
    record.denominatorValue = metric.denominator;
    record.value = metric.value;
    record.units = "dimensionless";
    record.aggregationScope = "single_pre_timing_validation_execution";
    record.limitOperator = "less_than_or_equal";
    record.limitValue = limit;
    record.passed = metric.value <= limit ? 1 : 0;
    return record;
}

/// Emit all applicable ordinary or rank-specific validation records for one candidate.
void emitValidationRecords( std::ostream &output,             /**< [out] machine-readable destination */
                            const RunIdentity &identity,      /**< [in] verified run identities */
                            std::string_view candidate,       /**< [in] precision-policy identifier */
                            const ValidationMetrics &metrics, /**< [in] validation evidence */
                            const BenchmarkOptions &options ) /**< [in] dimensions and limits */
{
    if( metrics.rank.evaluated )
    {
        MetricRecord rankRecord;
        rankRecord.gateId = "exact.rank_unambiguous";
        rankRecord.oracleId = "analytic_spectrum";
        rankRecord.metricName = "outside_ambiguity_retain_reject_agreement";
        rankRecord.numeratorNotApplicableReason = "exact_boolean_metric_has_no_numerator";
        rankRecord.denominatorNotApplicableReason = "exact_boolean_metric_has_no_denominator";
        rankRecord.value = metrics.rank.agreementOutsideAmbiguity ? 1.0 : 0.0;
        rankRecord.units = "none";
        rankRecord.aggregationScope = "every_oracle_eigenvalue_outside_candidate_specific_band";
        rankRecord.limitOperator = "exact_true";
        rankRecord.passed = metrics.rank.agreementOutsideAmbiguity ? 1 : 0;
        emitMetricRecord( output, identity, candidate, options, metrics.rank, rankRecord, 0, -1 );
        return;
    }

    const auto emit = [&]( const MetricRecord &record )
    { emitMetricRecord( output, identity, candidate, options, metrics.rank, record, 0, -1 ); };
    emit( normalizedGate( "smoke.gram_relative_frobenius",
                          "fp64-direct-Gram-from-original-fp32-values",
                          "gram_relative_frobenius",
                          "norm_frobenius_candidate_minus_fp64_gram",
                          "protected_norm_fp64_gram_frobenius",
                          metrics.gramRelativeError,
                          options.gramRelativeTolerance ) );
    emit( normalizedGate( "smoke.retained_eigenvalue_relative_l2",
                          "fp64_direct_thin_svd",
                          "retained_eigenvalue_relative_l2",
                          "norm_l2_candidate_minus_oracle_eigenvalues",
                          "protected_norm_oracle_retained_eigenvalues_l2",
                          metrics.eigenvalueRelativeError,
                          options.eigenvalueRelativeTolerance ) );
    emit( normalizedGate( "smoke.retained_projector_relative_frobenius",
                          "fp64_direct_thin_svd",
                          "retained_projector_relative_frobenius",
                          "norm_frobenius_candidate_minus_oracle_projector",
                          "protected_norm_oracle_projector_frobenius",
                          metrics.projectorRelativeError,
                          options.projectorRelativeTolerance ) );
    emit( normalizedGate( "smoke.residual_over_oracle_residual",
                          "fp64_direct_thin_svd",
                          "residual_error_over_oracle_residual",
                          "norm_l2_candidate_minus_oracle_residual",
                          "protected_norm_oracle_residual_l2",
                          metrics.residualOverOracleResidual,
                          options.residualRelativeTolerance ) );
    emit( normalizedGate( "smoke.residual_over_target",
                          "fp64_direct_thin_svd",
                          "residual_error_over_target",
                          "norm_l2_candidate_minus_oracle_residual",
                          "protected_norm_original_target_l2",
                          metrics.residualOverTarget,
                          options.residualRelativeTolerance ) );
    emit( normalizedGate( "smoke.normalized_solver_residual",
                          "intrinsic-candidate-eigensystem",
                          "normalized_solver_residual",
                          "norm_frobenius_gram_q_minus_q_lambda",
                          "protected_norm_gram_times_norm_q",
                          metrics.solverResidual,
                          options.solverResidualTolerance ) );

    MetricRecord orthogonality;
    orthogonality.gateId = "smoke.retained_factor_orthogonality";
    orthogonality.oracleId = "mathematical-identity";
    orthogonality.metricName = "retained_factor_orthogonality";
    orthogonality.numeratorName = "norm_frobenius_q_transpose_q_minus_identity";
    orthogonality.numeratorValue = metrics.orthogonalityDefect;
    orthogonality.denominatorNotApplicableReason = "formula_has_no_denominator";
    orthogonality.value = metrics.orthogonalityDefect;
    orthogonality.units = "dimensionless";
    orthogonality.aggregationScope = "retained_factor";
    orthogonality.limitOperator = "less_than_or_equal";
    orthogonality.limitValue = options.orthogonalityTolerance;
    orthogonality.passed = metrics.orthogonalityDefect <= options.orthogonalityTolerance ? 1 : 0;
    emit( orthogonality );

    MetricRecord gramMaximum;
    gramMaximum.gateId = "smoke.gram_maximum_absolute";
    gramMaximum.oracleId = "fp64-direct-Gram-from-original-fp32-values";
    gramMaximum.metricName = "gram_maximum_absolute";
    gramMaximum.numeratorName = "maximum_absolute_gram_entry_error";
    gramMaximum.numeratorValue = metrics.gramMaximumError;
    gramMaximum.denominatorNotApplicableReason = "maximum_absolute_metric_has_no_denominator";
    gramMaximum.value = metrics.gramMaximumError;
    gramMaximum.units = "squared_input_units";
    gramMaximum.aggregationScope = "all_gram_entries";
    gramMaximum.passed = -1;
    gramMaximum.gateStatus = "frozen-report-only";
    emit( gramMaximum );

    MetricRecord residualMaximum;
    residualMaximum.gateId = "smoke.residual_maximum_absolute";
    residualMaximum.oracleId = "fp64_direct_thin_svd";
    residualMaximum.metricName = "residual_maximum_absolute";
    residualMaximum.numeratorName = "maximum_absolute_residual_error";
    residualMaximum.numeratorValue = metrics.residualMaximumError;
    residualMaximum.denominatorNotApplicableReason = "maximum_absolute_metric_has_no_denominator";
    residualMaximum.value = metrics.residualMaximumError;
    residualMaximum.units = "input_intensity_units";
    residualMaximum.aggregationScope = "all_target_samples";
    residualMaximum.passed = -1;
    residualMaximum.gateStatus = "frozen-report-only";
    emit( residualMaximum );

    MetricRecord principalAngle;
    principalAngle.gateId = "smoke.maximum_principal_angle";
    principalAngle.oracleId = "fp64_direct_thin_svd";
    principalAngle.metricName = "maximum_principal_angle";
    principalAngle.numeratorNotApplicableReason = "angle_metric_has_no_numerator";
    principalAngle.denominatorNotApplicableReason = "angle_metric_has_no_denominator";
    principalAngle.value = metrics.maximumPrincipalAngle;
    principalAngle.units = "radians";
    principalAngle.aggregationScope = "retained_subspace";
    principalAngle.passed = -1;
    principalAngle.gateStatus = "frozen-report-only";
    emit( principalAngle );

    MetricRecord productionMapping;
    productionMapping.gateId = "smoke.production_mapping";
    productionMapping.oracleId = "production_api_same_policy";
    productionMapping.metricName = "production_mapping";
    productionMapping.numeratorNotApplicableReason = "metric_unavailable";
    productionMapping.denominatorNotApplicableReason = "metric_unavailable";
    productionMapping.units = "dimensionless";
    productionMapping.aggregationScope = "single_pre_timing_validation_execution";
    productionMapping.applicable = false;
    productionMapping.notApplicableReason = "benchmark_local_p4_staging_is_not_the_production_precision_adapter";
    productionMapping.passed = -1;
    emit( productionMapping );
}

/// Emit five non-gating stage diagnostics without claiming the owner-TBD end-to-end performance gate.
void emitTimingRecords( std::ostream &output,             /**< [out] machine-readable destination */
                        const RunIdentity &identity,      /**< [in] verified run identities */
                        std::string_view candidate,       /**< [in] precision-policy identifier */
                        const StageTimings &timing,       /**< [in] exclusive and total execution timings */
                        const ValidationMetrics &metrics, /**< [in] rank evidence shared by this candidate */
                        const BenchmarkOptions &options,  /**< [in] dimensions and controls */
                        int repetition,                   /**< [in] zero-based timed repetition */
                        int executionPosition )           /**< [in] balanced-order position */
{
    const std::array<std::string_view, 5> stageNames{ "conversion",
                                                      "gram_cross_product",
                                                      "eigensolve",
                                                      "projection_residual",
                                                      "total" };
    const std::array<double, 5> stageSeconds{ timing.conversionSeconds,
                                              timing.gramCrossProductSeconds,
                                              timing.eigensolveSeconds,
                                              timing.projectionResidualSeconds,
                                              timing.totalSeconds };
    for( std::size_t stage = 0; stage < stageNames.size(); ++stage )
    {
        MetricRecord record;
        record.record = "diagnostic_timing";
        record.gateId = "performance.end_to_end_wall_time";
        record.oracleId = "algorithm-authoritative-production";
        record.metricName = std::string( "microbenchmark_stage_seconds:" ) + std::string( stageNames[stage] );
        record.numeratorNotApplicableReason = "owner_tbd_end_to_end_metric_unavailable";
        record.denominatorNotApplicableReason = "owner_tbd_end_to_end_metric_unavailable";
        record.units = "seconds";
        record.aggregationScope = "single_timed_microbenchmark_execution";
        record.applicable = false;
        record.notApplicableReason =
            "kernel_stage_timing_is_not_owner_designated_end_to_end_performance_evidence_and_limit_is_owner_tbd";
        record.passed = -1;
        record.gateStatus = "owner_tbd";
        record.stage = stageNames[stage];
        record.diagnosticValue = stageSeconds[stage];
        record.diagnosticUnits = "seconds";
        emitMetricRecord( output, identity, candidate, options, metrics.rank, record, repetition, executionPosition );
    }
}

/// Update one candidate's timing samples and component-wise worst validation values.
void accumulateSummary( CandidateSummary &summary,  /**< [in,out] aggregate candidate result */
                        const StageTimings &timing, /**< [in] current execution timings */
                        const ValidationMetrics &metrics /**< [in] current validation result */ )
{
    const std::array<double, 5> seconds{ timing.conversionSeconds,
                                         timing.gramCrossProductSeconds,
                                         timing.eigensolveSeconds,
                                         timing.projectionResidualSeconds,
                                         timing.totalSeconds };
    for( std::size_t stage = 0; stage < seconds.size(); ++stage )
    {
        summary.stageSeconds[stage].push_back( seconds[stage] );
    }
    summary.maximumMetrics = metrics;
    summary.passed = summary.passed && metrics.passed;
}

/// Return the median of one nonempty collection of raw timing samples.
double median( std::vector<double> values /**< [in] samples copied for sorting */ )
{
    std::sort( values.begin(), values.end() );
    const std::size_t middle = values.size() / 2;
    if( values.size() % 2 == 0 )
    {
        return 0.5 * ( values[middle - 1] + values[middle] );
    }
    return values[middle];
}

/// Return an environment value or a stable unset marker.
std::string_view environmentValue( const char *name /**< [in] process environment variable name */ )
{
    const char *value = std::getenv( name );
    return value ? std::string_view( value ) : std::string_view( "<unset>" );
}

/// Print concise median timings, validation maxima, gates, and thread controls.
void printHumanSummary( std::ostream &output,            /**< [out] human-readable destination */
                        const BenchmarkOptions &options, /**< [in] dimensions and declared gates */
                        const std::array<CandidateSummary, 3> &summaries /**< [in] aggregate candidate results */ )
{
    const std::array<std::string_view, 3> names{ "P4-D64", "P4-M32D64", "P4-F32" };
    output << "P4 CPU precision microbenchmark (" << inputCaseName( options.inputCase )
           << ", uncentered, in-sample): T=" << options.samples << ", K=" << options.predictors
           << ", orientation=" << ( options.samples <= options.predictors ? "temporal" : "predictor" )
           << ", retained modes=" << options.modes << ", warmups=" << options.warmups
           << ", repetitions=" << options.repetitions << '\n';
    output << "Thread environment: OPENBLAS_NUM_THREADS=" << environmentValue( "OPENBLAS_NUM_THREADS" )
           << ", OMP_NUM_THREADS=" << environmentValue( "OMP_NUM_THREADS" ) << '\n';
    output << "Candidate order rotates once per repetition; "
           << ( options.repetitions % 3 == 0 ? "this run contains complete balanced triplets."
                                             : "this run ends with an incomplete, unbalanced triplet." )
           << '\n';
    output << "Median milliseconds:\n"
           << std::left << std::setw( 13 ) << "candidate" << std::right << std::setw( 12 ) << "convert"
           << std::setw( 14 ) << "gram/cross" << std::setw( 12 ) << "eigensolve" << std::setw( 14 ) << "projection"
           << std::setw( 12 ) << "total" << '\n';
    output << std::fixed << std::setprecision( 3 );
    for( std::size_t candidate = 0; candidate < summaries.size(); ++candidate )
    {
        output << std::left << std::setw( 13 ) << names[candidate] << std::right;
        for( const auto &stage : summaries[candidate].stageSeconds )
        {
            output << std::setw( 12 ) << 1000.0 * median( stage );
        }
        output << '\n';
    }

    output << std::scientific << std::setprecision( 3 );
    output << "Single pre-timing validation (gram, eigenvalue, projector, residual/oracle, residual/target, solver "
              "residual, orthogonality):\n";
    for( std::size_t candidate = 0; candidate < summaries.size(); ++candidate )
    {
        const ValidationMetrics &metrics = summaries[candidate].maximumMetrics;
        output << "  " << names[candidate] << ": " << metrics.gramRelativeError.value << ", "
               << metrics.eigenvalueRelativeError.value << ", " << metrics.projectorRelativeError.value << ", "
               << metrics.residualOverOracleResidual.value << ", " << metrics.residualOverTarget.value << ", "
               << metrics.solverResidual.value << ", " << metrics.orthogonalityDefect << " ["
               << ( summaries[candidate].passed ? "PASS" : "FAIL" ) << "]\n";
    }
    output << "Engineering smoke gates (not final science tolerances): " << options.gramRelativeTolerance << ", "
           << options.eigenvalueRelativeTolerance << ", " << options.projectorRelativeTolerance << ", "
           << options.residualRelativeTolerance << ", " << options.solverResidualTolerance << ", "
           << options.orthogonalityTolerance << '\n';
    if( options.inputCase != InputCase::wellConditioned )
    {
        output << "Rank diagnostics use strict lambda/lambda-max > " << options.rankRelativeThreshold
               << "; only decisions outside each construction-specific precision band are gated.\n";
        if( options.inputCase == InputCase::rankZero )
        {
            output << "At zero tolerance the oracle classifies its exact nullspace with the FP64 SVD structural "
                      "threshold; candidate ranks use the production strict eigenvalue > 0 comparison.\n";
        }
        output << "Independent FP64 SVD oracle rank=" << summaries[0].maximumMetrics.rank.oracleRank
               << ", structural rank=" << summaries[0].maximumMetrics.rank.oracleStructuralRank << " (designed "
               << summaries[0].maximumMetrics.rank.designedStructuralRank
               << "), exact nullity=" << summaries[0].maximumMetrics.rank.oracleExactNullity
               << ", smallest positive ratio=" << summaries[0].maximumMetrics.rank.oracleSmallestPositiveRatio
               << ", raw largest null ratio=" << summaries[0].maximumMetrics.rank.oracleRawLargestNullRatio
               << ", Gram norm scales=" << summaries[0].maximumMetrics.rank.oracleGramFrobeniusRatio << '/'
               << summaries[0].maximumMetrics.rank.oracleAbsoluteGramRatio
               << ", nearest threshold distance=" << summaries[0].maximumMetrics.rank.oracleMinimumThresholdDistance
               << ", boundary ratios=[" << summaries[0].maximumMetrics.rank.oracleBoundary.ratioBelow << ", "
               << summaries[0].maximumMetrics.rank.oracleBoundary.ratioAbove
               << "], eigengap=" << summaries[0].maximumMetrics.rank.oracleBoundary.eigengap << '\n';
        output << "Candidate rank diagnostics (rank, band half-width, ambiguous values, candidate boundary gap, "
                  "unambiguous agreement):\n";
        for( std::size_t candidate = 0; candidate < summaries.size(); ++candidate )
        {
            const RankMetrics &rank = summaries[candidate].maximumMetrics.rank;
            output << "  " << names[candidate] << ": " << rank.candidateRank << ", " << rank.ambiguityHalfWidth << ", "
                   << rank.ambiguousEigenvalues << ", " << rank.candidateBoundary.eigengap << " ["
                   << ( rank.agreementOutsideAmbiguity ? "PASS_UNAMBIGUOUS" : "FAIL_UNAMBIGUOUS" );
            if( rank.rankDiffersFromOracle )
            {
                output << "; RAW_RANK_DIFFERS";
            }
            output << "]\n";
        }
        output << "The ambiguity model applies only to this bounded "
               << ( options.inputCase == InputCase::rankBoundary ? "diagonal spectral" : "dense exact-dependency" )
               << " construction; it is not a general error bound or final science tolerance.\n";
    }
}

} // namespace

/// Run the CPU precision benchmark and return nonzero when execution or validation fails.
int main( int argc, /**< [in] command-line argument count */
          char **argv /**< [in] command-line argument vector */ )
{
    try
    {
        bool showHelp{ false };
        const BenchmarkOptions options = parseOptions( argc, argv, showHelp );
        if( showHelp )
        {
            printUsage( std::cout, argv[0] );
            return 0;
        }

        std::ofstream csvFile;
        std::ofstream summaryFile;
        std::ostream *csvOutput = &std::cout;
        std::ostream *summaryOutput = &std::cerr;
        if( !options.csvOutputPath.empty() )
        {
            csvFile.open( options.csvOutputPath );
            if( !csvFile )
            {
                throw std::runtime_error( "cannot open --csv-output path: " + options.csvOutputPath );
            }
            csvOutput = &csvFile;
        }
        if( !options.summaryOutputPath.empty() )
        {
            summaryFile.open( options.summaryOutputPath );
            if( !summaryFile )
            {
                throw std::runtime_error( "cannot open --summary-output path: " + options.summaryOutputPath );
            }
            summaryOutput = &summaryFile;
        }

        const SourceData source = makeSelectedSourceData( options );
        const RunIdentity identity = makeRunIdentity( options, source, argv[0] );
        const NumericalOracle numericalOracle = makeNumericalOracle( source, options );
        std::array<NumericalResult, 3> numericalResults;
        {
            CandidateWorkspace<double, double> validationD64Workspace;
            CandidateWorkspace<float, double> validationMixedWorkspace;
            CandidateWorkspace<float, float> validationF32Workspace;
            executeCandidate( "P4-D64", source, options, validationD64Workspace );
            numericalResults[0] = captureNumericalResult( "P4-D64", source, options, validationD64Workspace );
            executeCandidate( "P4-M32D64", source, options, validationMixedWorkspace );
            numericalResults[1] = captureNumericalResult( "P4-M32D64", source, options, validationMixedWorkspace );
            executeCandidate( "P4-F32", source, options, validationF32Workspace );
            numericalResults[2] = captureNumericalResult( "P4-F32", source, options, validationF32Workspace );
        }

        std::array<RankMetrics, 3> rankMetrics;
        if( options.inputCase != InputCase::wellConditioned )
        {
            const RankOracle oracle = makeRankOracle( source, options );
            rankMetrics[0] = validateRank<double, double>( numericalResults[0], oracle, options );
            rankMetrics[1] = validateRank<float, double>( numericalResults[1], oracle, options );
            rankMetrics[2] = validateRank<float, float>( numericalResults[2], oracle, options );
        }

        std::array<ValidationMetrics, 3> validationMetrics;
        for( std::size_t candidate = 0; candidate < validationMetrics.size(); ++candidate )
        {
            validationMetrics[candidate] =
                validateResult( numericalResults[candidate], numericalOracle, options, rankMetrics[candidate] );
        }

        CandidateWorkspace<double, double> d64Workspace;
        CandidateWorkspace<float, double> mixedWorkspace;
        CandidateWorkspace<float, float> f32Workspace;
        for( int warmup = 0; warmup < options.warmups; ++warmup )
        {
            executeCandidate( "P4-D64", source, options, d64Workspace );
            executeCandidate( "P4-M32D64", source, options, mixedWorkspace );
            executeCandidate( "P4-F32", source, options, f32Workspace );
        }

        emitCsvHeader( *csvOutput );
        std::array<CandidateSummary, 3> summaries;
        const std::array<std::string_view, 3> candidateNames{ "P4-D64", "P4-M32D64", "P4-F32" };
        for( std::size_t candidate = 0; candidate < validationMetrics.size(); ++candidate )
        {
            emitValidationRecords( *csvOutput,
                                   identity,
                                   candidateNames[candidate],
                                   validationMetrics[candidate],
                                   options );
        }
        for( int repetition = 0; repetition < options.repetitions; ++repetition )
        {
            std::array<StageTimings, 3> timings;
            std::array<int, 3> executionPositions;
            for( int offset = 0; offset < 3; ++offset )
            {
                const int candidate = ( repetition + offset ) % 3;
                executionPositions[static_cast<std::size_t>( candidate )] = offset;
                if( candidate == 0 )
                {
                    timings[0] = executeCandidate( "P4-D64", source, options, d64Workspace );
                }
                else if( candidate == 1 )
                {
                    timings[1] = executeCandidate( "P4-M32D64", source, options, mixedWorkspace );
                }
                else
                {
                    timings[2] = executeCandidate( "P4-F32", source, options, f32Workspace );
                }
            }

            for( std::size_t candidate = 0; candidate < timings.size(); ++candidate )
            {
                emitTimingRecords( *csvOutput,
                                   identity,
                                   candidateNames[candidate],
                                   timings[candidate],
                                   validationMetrics[candidate],
                                   options,
                                   repetition,
                                   executionPositions[candidate] );
                accumulateSummary( summaries[candidate], timings[candidate], validationMetrics[candidate] );
            }
        }

        printHumanSummary( *summaryOutput, options, summaries );
        csvOutput->flush();
        summaryOutput->flush();
        if( !*csvOutput )
        {
            throw std::runtime_error( "failed while writing CSV output" );
        }
        if( !*summaryOutput )
        {
            throw std::runtime_error( "failed while writing summary output" );
        }
        const bool passed = std::all_of( summaries.begin(),
                                         summaries.end(),
                                         []( const CandidateSummary &summary ) { return summary.passed; } );
        if( !passed )
        {
            *summaryOutput << "Numerical validation failed; see CSV metrics and declared engineering smoke gates.\n";
            summaryOutput->flush();
            if( !*summaryOutput )
            {
                throw std::runtime_error( "failed while writing summary output" );
            }
            return 2;
        }
        return 0;
    }
    catch( const std::exception &error )
    {
        std::cerr << "p4PrecisionBenchmark: " << error.what() << '\n';
        return 1;
    }
}
