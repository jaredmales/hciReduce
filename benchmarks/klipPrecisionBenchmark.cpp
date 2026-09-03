/** \file klipPrecisionBenchmark.cpp
 * \brief Benchmarks CPU precision policies for direct KLIP mode construction and subtraction.
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
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/SVD>

#include "../src/common/ADIDerotator.hpp"
#include "../src/common/KLIPreduction.hpp"

namespace
{

/// Clock used for every benchmark interval.
using benchmarkClockT = std::chrono::steady_clock;

/// Dynamic column-major array used by the benchmark pipeline.
template <typename scalarT>
using matrixT = Eigen::Array<scalarT, Eigen::Dynamic, Eigen::Dynamic>;

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

/// Supported KLIP execution topology.
enum class BenchmarkScenario : std::uint8_t
{
    sharedReference, ///< Construct one reference-Gram basis and apply it to every target.
    targetPixel      ///< Construct a new pixel-Gram basis independently for every target.
};

/// Runtime controls and declared engineering smoke-test limits.
struct BenchmarkOptions
{
    BenchmarkScenario scenario{ BenchmarkScenario::sharedReference }; ///< KLIP basis-reuse topology.
    int targets{ 32 };                                                ///< Number of target vectors, T.
    int pixels{ 1024 };                                               ///< Pixels in the vectorized region, P.
    int references{ 64 };                                             ///< References in the selected library, R.
    int modes{ 16 };                                                  ///< Number of leading KL modes subtracted.
    int warmups{ 1 };                                                 ///< Untimed pipeline executions per candidate.
    int repetitions{ 6 };                 ///< Timed executions per candidate; even counts balance alternating order.
    std::uint64_t seed{ 42 };             ///< Seed for deterministic FP32 source values.
    double gramRelativeTolerance{ 5e-5 }; ///< Accepted FP32 Gram relative Frobenius error.
    double eigenvalueRelativeTolerance{ 5e-4 }; ///< Accepted retained-eigenvalue relative L2 error.
    double projectorRelativeTolerance{ 2e-2 };  ///< Accepted retained-projector relative Frobenius error.
    double residualRelativeTolerance{ 2e-3 };   ///< Accepted final-residual relative Frobenius error.
    double productionRelativeTolerance{ 5e-5 }; ///< Accepted difference from calcKLModesAdaptive output.
    double orthogonalityTolerance{ 2e-3 };      ///< Accepted retained-mode orthogonality defect.
    double solverResidualTolerance{ 2e-4 };     ///< Accepted normalized eigensolver residual.
    std::string acceptanceContractPath;         ///< Runtime path whose content must match the frozen contract.
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
    double gramSeconds{ 0 };             ///< Selected reference- or pixel-Gram construction time.
    double eigensolveSeconds{ 0 };       ///< Gram conversion, SYEVR, and eigenpair conversion time.
    double modeConstructionSeconds{ 0 }; ///< Reference projection or pixel-mode materialization time.
    double subtractionSeconds{ 0 };      ///< Coefficient, PSF, and residual construction time.
    double totalSeconds{ 0 };            ///< Complete pipeline wall time, including small loop overheads.
};

/// Shared FP32 values presented to both precision candidates.
struct SourceData
{
    matrixT<float> references; ///< P-by-R vectorized reference library.
    matrixT<float> targets;    ///< P-by-T target vectors.
};

/// Reusable storage for one eigensolver scalar policy.
template <typename eigensolverT>
struct CandidateWorkspace
{
    matrixT<float> gram;                              ///< FP32 Gram matrix constructed by the KLIP calculation path.
    mx::math::syevrMem<eigensolverT> solverWorkspace; ///< Reusable production mxlib SYEVR storage.
    matrixT<float> eigenvectors;                      ///< Selected eigenvectors returned to FP32 calculation storage.
    matrixT<float> eigenvalues;                       ///< Selected eigenvalues returned to FP32 calculation storage.
    matrixT<float> modes;                             ///< KL modes stored by row in ascending eigenvalue order.
    matrixT<float> coefficients;                      ///< Mode coefficients for one target.
    matrixT<float> psf;                               ///< Reconstructed PSF for one target.
    matrixT<float> residual;                          ///< Reusable one-target residual scratch.
    matrixT<float> residuals;                         ///< P-by-T final residual vectors.
};

/// Numerical results retained from one separate untimed execution.
struct NumericalResult
{
    std::string candidate;        ///< Stable precision-policy identifier.
    matrixT<double> gram;         ///< FP32 candidate Gram promoted only for validation.
    matrixT<double> eigenvectors; ///< Selected candidate Gram factors in ascending eigenvalue order.
    matrixT<double> eigenvalues;  ///< Selected ascending eigenvalues.
    matrixT<double> modes;        ///< Candidate KL modes stored by ascending eigenvalue.
    matrixT<double> residuals;    ///< Candidate residuals for every target.
};

/// Independent FP64 thin-SVD reference calculated once per invocation.
struct OracleResult
{
    Eigen::MatrixXd gram;        ///< Complete direct FP64 Gram matrix.
    Eigen::VectorXd eigenvalues; ///< Leading squared singular values in descending order.
    Eigen::MatrixXd modeColumns; ///< Leading left singular vectors stored by descending column.
    Eigen::MatrixXd residuals;   ///< FP64 projected residuals for every target.
    double targetNorm{ 0 };      ///< Frobenius norm used to scale residual errors.
    double residualNorm{ 0 };    ///< Frobenius norm of the FP64 oracle residuals.
};

/// Numerical comparisons against the FP64 SVD oracle and production implementation.
struct ValidationMetrics
{
    double gramRelativeError{ 0 };       ///< Relative Frobenius Gram error against direct FP64 arithmetic.
    double eigenvalueRelativeError{ 0 }; ///< Relative L2 retained-eigenvalue error against squared singular values.
    double projectorRelativeError{ 0 };  ///< Relative Frobenius retained-projector error against FP64 thin SVD.
    double residualRelativeError{ 0 };   ///< FP64-oracle residual error normalized by target Frobenius norm.
    double residualMaximumError{ 0 };    ///< Maximum absolute final-residual error against the FP64 oracle.
    double productionProjectorRelativeError{ 0 }; ///< Projector difference from calcKLModesAdaptive.
    double productionResidualRelativeError{ 0 };  ///< Residual difference from calcKLModesAdaptive modes.
    double orthogonalityDefect{ 0 };              ///< Frobenius defect of the retained spatial modes.
    NormalizedMetric gramEvidence;                ///< Gram numerator, denominator, and quotient.
    NormalizedMetric eigenvalueEvidence;          ///< Retained-eigenvalue numerator, denominator, and quotient.
    NormalizedMetric projectorEvidence;           ///< Retained-projector numerator, denominator, and quotient.
    NormalizedMetric residualOverOracleResidual;  ///< Residual error divided by oracle residual norm.
    NormalizedMetric residualOverTarget;          ///< Residual error divided by original target norm.
    NormalizedMetric productionProjectorEvidence; ///< Staged-to-production projector metric evidence.
    NormalizedMetric productionResidualEvidence;  ///< Staged-to-production residual metric evidence.
    NormalizedMetric solverResidualEvidence;      ///< Candidate eigensystem backward-residual evidence.
    double gramMaximumError{ 0 };                 ///< Maximum absolute authoritative Gram-entry error.
    double maximumPrincipalAngle{ 0 };            ///< Largest retained-subspace principal angle in radians.
    bool passed{ false };                         ///< Whether every declared engineering smoke gate passed.
};

/// Repeated timing and validation values for one candidate.
struct CandidateSummary
{
    std::array<std::vector<double>, 5> stageSeconds; ///< Raw samples for four exclusive stages plus total.
    ValidationMetrics metrics;                       ///< Validation values from the single untimed execution.
    bool passed{ true };                             ///< Whether the candidate passed its validation execution.
};

/// Return elapsed seconds since a steady-clock time point.
double secondsSince( const benchmarkClockT::time_point &start /**< [in] beginning of the interval */ )
{
    return std::chrono::duration<double>( benchmarkClockT::now() - start ).count();
}

/// Return the stable command-line and CSV spelling of a benchmark scenario.
std::string_view scenarioName( BenchmarkScenario scenario /**< [in] scenario to identify */ )
{
    if( scenario == BenchmarkScenario::sharedReference )
    {
        return "shared-reference";
    }
    return "target-pixel";
}

/// Parse an exact benchmark-scenario spelling.
BenchmarkScenario parseScenario( std::string_view value /**< [in] command-line scenario value */ )
{
    if( value == "shared-reference" )
    {
        return BenchmarkScenario::sharedReference;
    }
    if( value == "target-pixel" )
    {
        return BenchmarkScenario::targetPixel;
    }
    throw std::invalid_argument( "--scenario must be shared-reference or target-pixel" );
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
void printUsage( std::ostream &output,         /**< [out] destination stream */
                 std::string_view executable ) /**< [in] invoked program name */
{
    output << "Usage: " << executable << " [options]\n"
           << "  --scenario NAME                shared-reference or target-pixel (default shared-reference)\n"
           << "  --targets T                    Target-vector count (default 32)\n"
           << "  --pixels P                     Vectorized region pixels (default 1024)\n"
           << "  --references R                 Selected reference count (default 64)\n"
           << "  --modes N                      Leading subtracted modes (default 16)\n"
           << "  --warmups N                    Positive untimed executions per candidate (default 1)\n"
           << "  --repetitions N                Timed executions; even counts balance order (default 6)\n"
           << "  --seed N                       Deterministic FP32 source seed (default 42)\n"
           << "  --gram-rel-tol X               Frozen Gram engineering smoke gate (5e-5)\n"
           << "  --eigenvalue-rel-tol X         Frozen eigenvalue engineering smoke gate (5e-4)\n"
           << "  --projector-rel-tol X          Frozen projector engineering smoke gate (2e-2)\n"
           << "  --residual-rel-tol X           Frozen residual engineering smoke gate (2e-3)\n"
           << "  --production-rel-tol X         Frozen production-agreement smoke gate (5e-5)\n"
           << "  --orthogonality-tol X          Frozen orthogonality engineering smoke gate (2e-3)\n"
           << "  --solver-residual-tol X        Frozen eigensolver engineering smoke gate (2e-4)\n"
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
           << "Defaults are engineering capability smoke gates, not final scientific acceptance tolerances.\n"
           << "shared-reference requires R <= P and constructs one basis for all T targets. target-pixel requires "
              "P < R and repeats a direct adaptive refit for each target. The latter uses one selected in-memory "
              "library so reference selection/collapse and OpenMP scheduling remain outside this numerical "
              "microbenchmark. Candidate matrices are reused after warmup, so transient per-call allocations are "
              "also excluded. One M-mode residual is materialized per target; production multi-mode-list output "
              "materialization is excluded. This is arithmetic staging, not complete KLIP worker timing.\n";
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
        if( name == "--scenario" )
        {
            options.scenario = parseScenario( value );
        }
        else if( name == "--targets" )
        {
            options.targets = parseInteger( value, name, false );
        }
        else if( name == "--pixels" )
        {
            options.pixels = parseInteger( value, name, false );
        }
        else if( name == "--references" )
        {
            options.references = parseInteger( value, name, false );
        }
        else if( name == "--modes" )
        {
            options.modes = parseInteger( value, name, false );
        }
        else if( name == "--warmups" )
        {
            options.warmups = parseInteger( value, name, false );
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
        else if( name == "--production-rel-tol" )
        {
            options.productionRelativeTolerance = parseTolerance( value, name );
        }
        else if( name == "--orthogonality-tol" )
        {
            options.orthogonalityTolerance = parseTolerance( value, name );
        }
        else if( name == "--solver-residual-tol" )
        {
            options.solverResidualTolerance = parseTolerance( value, name );
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
        options.productionRelativeTolerance != 5e-5 || options.solverResidualTolerance != 2e-4 ||
        options.orthogonalityTolerance != 2e-3 )
    {
        throw std::invalid_argument( "contract-bound mode requires the frozen v1 engineering smoke limits" );
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

    if( options.modes > std::min( options.pixels, options.references ) )
    {
        throw std::invalid_argument( "--modes must not exceed min(P,R)" );
    }
    if( options.scenario == BenchmarkScenario::sharedReference && options.references > options.pixels )
    {
        throw std::invalid_argument( "shared-reference requires R <= P" );
    }
    if( options.scenario == BenchmarkScenario::targetPixel && options.pixels >= options.references )
    {
        throw std::invalid_argument( "target-pixel requires P < R" );
    }
    return options;
}

/// Generate one exactly representable pseudorandom FP32 value in [-1,1).
float nextUnitFloat( std::mt19937 &generator /**< [in,out] deterministic integer generator */ )
{
    const std::uint32_t mantissa = generator() >> 8;
    return static_cast<float>( mantissa ) * ( 1.0f / 8388608.0f ) - 1.0f;
}

/// Construct one deterministic, well-conditioned FP32 reference library and correlated targets.
SourceData makeSourceData( const BenchmarkOptions &options /**< [in] dimensions and seed */ )
{
    SourceData source;
    source.references.resize( options.pixels, options.references );
    source.targets.resize( options.pixels, options.targets );
    std::seed_seq seedSequence{ static_cast<std::uint32_t>( options.seed ),
                                static_cast<std::uint32_t>( options.seed >> 32 ) };
    std::mt19937 generator( seedSequence );

    for( Eigen::Index reference = 0; reference < source.references.cols(); ++reference )
    {
        const float position = static_cast<float>( reference + 1 ) / static_cast<float>( source.references.cols() );
        const float scale = 0.75f + 0.5f * position;
        for( Eigen::Index pixel = 0; pixel < source.references.rows(); ++pixel )
        {
            source.references( pixel, reference ) = scale * nextUnitFloat( generator );
        }
    }

    for( Eigen::Index target = 0; target < source.targets.cols(); ++target )
    {
        const Eigen::Index first = target % source.references.cols();
        const Eigen::Index second = ( 5 * target + source.references.cols() / 3 ) % source.references.cols();
        const Eigen::Index third = ( 11 * target + 2 * source.references.cols() / 3 ) % source.references.cols();
        for( Eigen::Index pixel = 0; pixel < source.targets.rows(); ++pixel )
        {
            source.targets( pixel, target ) =
                0.55f * source.references( pixel, first ) - 0.30f * source.references( pixel, second ) +
                0.20f * source.references( pixel, third ) + 0.08f * nextUnitFloat( generator );
        }
    }
    return source;
}

/// Hash the exact source arrays using a canonical versioned binary32 serialization.
std::string sourceValueSha256( const SourceData &source /**< [in] source arrays to identify */ )
{
    Sha256 digest;
    constexpr std::string_view encoding{ "hcireduce-klip-source-binary32-big-endian-column-major/v1" };
    digest.update( encoding.data(), encoding.size() );
    hashUnsigned( digest, static_cast<std::uint64_t>( source.references.rows() ) );
    hashUnsigned( digest, static_cast<std::uint64_t>( source.references.cols() ) );
    hashFloatValues( digest, source.references.data(), static_cast<std::size_t>( source.references.size() ) );
    hashUnsigned( digest, static_cast<std::uint64_t>( source.targets.rows() ) );
    hashUnsigned( digest, static_cast<std::uint64_t>( source.targets.cols() ) );
    hashFloatValues( digest, source.targets.data(), static_cast<std::size_t>( source.targets.size() ) );
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
    identity.caseSetId = "klip_well_separated";
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

/// Construct the selected FP32 Gram matrix used by the adaptive KLIP branch.
void constructGram( matrixT<float> &gram,             /**< [out] selected symmetric Gram matrix */
                    const matrixT<float> &references, /**< [in] P-by-R selected reference library */
                    BenchmarkScenario scenario )      /**< [in] branch whose shape has already been validated */
{
    if( scenario == BenchmarkScenario::sharedReference )
    {
        mx::math::eigenSYRK( gram, references );
    }
    else
    {
        gram = ( references.matrix() * references.matrix().transpose() ).array();
    }
}

/// Solve for the largest Gram eigenpairs with the production SYEVR range and conversion policy.
/** The FP32 covariance copy, including its FP64 promotion for KLIP-M32D64, is part of this eigensolve stage. The
 * reference-Gram branch also repeats the same-type `mem->cvd`, `mem->evecsd`, and `mem->evalsd` assignments made when
 * `calcKLModes` aliases those arrays through `calcEigenVecs`. Its normalization and validity checks occur in
 * eigensolver precision before conversion to FP32. The pixel-Gram branch performs one covariance conversion and
 * converts raw eigenpairs before its mode-stage validity checks.
 */
template <typename eigensolverT>
void solveGram( std::string_view candidate, /**< [in] stable precision-policy identifier */
                int modes,                  /**< [in] number of largest eigenpairs requested */
                BenchmarkScenario scenario, /**< [in] adaptive branch controlling normalization precision */
                CandidateWorkspace<eigensolverT> &state /**< [in,out] Gram, eigenpairs, and solver storage */ )
{
    const int gramDimension = state.gram.rows();
    state.solverWorkspace.cvd = state.gram.template cast<eigensolverT>();
    if( scenario == BenchmarkScenario::sharedReference )
    {
        state.solverWorkspace.cvd = state.solverWorkspace.cvd.template cast<eigensolverT>();
    }
    const MXLAPACK_INT status = mx::math::eigenSYEVR( state.solverWorkspace.evecsd,
                                                      state.solverWorkspace.evalsd,
                                                      state.solverWorkspace.cvd,
                                                      gramDimension - modes,
                                                      gramDimension,
                                                      'L',
                                                      &state.solverWorkspace );
    if( status != 0 )
    {
        throw std::runtime_error( std::string( candidate ) + " eigenSYEVR failed with status " +
                                  std::to_string( status ) );
    }

    if( scenario == BenchmarkScenario::sharedReference )
    {
        state.solverWorkspace.evecsd = state.solverWorkspace.evecsd.template cast<eigensolverT>();
        state.solverWorkspace.evalsd = state.solverWorkspace.evalsd.template cast<eigensolverT>();
        for( int mode = 0; mode < modes; ++mode )
        {
            bool validMode = mx::math::isFinite( state.solverWorkspace.evalsd( mode ) ) &&
                             state.solverWorkspace.evalsd( mode ) > eigensolverT{ 0 };
            if( validMode )
            {
                state.solverWorkspace.evecsd.col( mode ) /= std::sqrt( state.solverWorkspace.evalsd( mode ) );
            }
            for( int reference = 0; validMode && reference < state.solverWorkspace.evecsd.rows(); ++reference )
            {
                validMode = mx::math::isFinite( state.solverWorkspace.evecsd( reference, mode ) );
            }
            if( !validMode )
            {
                state.solverWorkspace.evecsd.col( mode ).setZero();
            }
        }
    }

    state.eigenvectors = state.solverWorkspace.evecsd.template cast<float>();
    state.eigenvalues = state.solverWorkspace.evalsd.template cast<float>();
}

/// Materialize spatial KL modes exactly as the selected production adaptive branch does.
void constructModes( matrixT<float> &modes,              /**< [out] M-by-P spatial KL modes */
                     const matrixT<float> &eigenvectors, /**< [in] selected Gram eigenvectors */
                     const matrixT<float> &eigenvalues,  /**< [in] selected Gram eigenvalues */
                     const matrixT<float> &references,   /**< [in] P-by-R selected reference library */
                     BenchmarkScenario scenario )        /**< [in] selected adaptive Gram branch */
{
    const int modeCount = eigenvectors.cols();
    if( scenario == BenchmarkScenario::sharedReference )
    {
        modes.resize( modeCount, references.rows() );
        mx::math::gemm<float>( CblasColMajor,
                               CblasTrans,
                               CblasTrans,
                               modeCount,
                               references.rows(),
                               references.cols(),
                               1.0F,
                               eigenvectors.data(),
                               eigenvectors.rows(),
                               references.data(),
                               references.rows(),
                               0.0F,
                               modes.data(),
                               modes.rows() );
        return;
    }

    modes = eigenvectors.matrix().transpose().array();
    for( int mode = 0; mode < modeCount; ++mode )
    {
        bool validMode = mx::math::isFinite( eigenvalues( mode ) ) && eigenvalues( mode ) > 0;
        for( int pixel = 0; validMode && pixel < modes.cols(); ++pixel )
        {
            validMode = mx::math::isFinite( modes( mode, pixel ) );
        }
        if( !validMode )
        {
            modes.row( mode ).setZero();
        }
    }
}

/// Apply all available leading modes to one target using the production KLIP accumulation order.
template <typename targetT>
void subtractModes( matrixT<float> &residual,     /**< [out] P-by-1 target residual */
                    matrixT<float> &coefficients, /**< [in,out] reusable 1-by-M coefficients */
                    matrixT<float> &psf,          /**< [in,out] reusable 1-by-P PSF estimate */
                    const matrixT<float> &modes,  /**< [in] M-by-P modes in ascending eigenvalue order */
                    const targetT &target )       /**< [in] contiguous P-by-1 target expression */
{
    coefficients.resize( 1, modes.rows() );
    for( int mode = 0; mode < coefficients.size(); ++mode )
    {
        coefficients( mode ) = modes.row( mode ).matrix().dot( target.col( 0 ).matrix() );
    }

    psf = coefficients( coefficients.size() - 1 ) * modes.row( modes.rows() - 1 );
    for( int mode = coefficients.size() - 2; mode >= 0; --mode )
    {
        psf += coefficients( mode ) * modes.row( mode );
    }
    residual = target - psf.transpose();
}

/// Execute direct adaptive KLIP mode construction and subtraction for one precision policy.
/** This benchmark-local staging mirrors `calcKLModesAdaptive`: FP32 Gram construction, `calcEigenVecs`-equivalent
 * promotion and largest-range SYEVR, branch-specific spatial-mode construction, and the worker's coefficient/PSF
 * accumulation order. A separate untimed validation calls the production function directly to guard that mapping.
 */
template <typename eigensolverT>
StageTimings executeCandidate( std::string_view candidate,      /**< [in] stable precision-policy identifier */
                               const SourceData &source,        /**< [in] common FP32 source values */
                               const BenchmarkOptions &options, /**< [in] topology, dimensions, and retained modes */
                               CandidateWorkspace<eigensolverT> &state /**< [in,out] reusable candidate storage */ )
{
    StageTimings timing;
    state.residuals.resize( options.pixels, options.targets );
    const auto totalStart = benchmarkClockT::now();

    const int fitCount = options.scenario == BenchmarkScenario::sharedReference ? 1 : options.targets;
    for( int fit = 0; fit < fitCount; ++fit )
    {
        const auto gramStart = benchmarkClockT::now();
        constructGram( state.gram, source.references, options.scenario );
        timing.gramSeconds += secondsSince( gramStart );

        const auto eigensolveStart = benchmarkClockT::now();
        solveGram( candidate, options.modes, options.scenario, state );
        timing.eigensolveSeconds += secondsSince( eigensolveStart );

        const auto modeStart = benchmarkClockT::now();
        constructModes( state.modes, state.eigenvectors, state.eigenvalues, source.references, options.scenario );
        timing.modeConstructionSeconds += secondsSince( modeStart );

        const int firstTarget = options.scenario == BenchmarkScenario::sharedReference ? 0 : fit;
        const int pastLastTarget = options.scenario == BenchmarkScenario::sharedReference ? options.targets : fit + 1;
        const auto subtractionStart = benchmarkClockT::now();
        for( int target = firstTarget; target < pastLastTarget; ++target )
        {
            subtractModes( state.residual, state.coefficients, state.psf, state.modes, source.targets.col( target ) );
            state.residuals.col( target ) = state.residual;
        }
        timing.subtractionSeconds += secondsSince( subtractionStart );
    }

    timing.totalSeconds = secondsSince( totalStart );
    return timing;
}

/// Capture one candidate's completed outputs in FP64 validation storage.
template <typename eigensolverT>
NumericalResult captureNumericalResult(
    std::string_view candidate,                        /**< [in] stable precision-policy identifier */
    const CandidateWorkspace<eigensolverT> &workspace, /**< [in] completed untimed candidate execution */
    BenchmarkScenario scenario )                       /**< [in] Gram branch controlling valid matrix storage */
{
    NumericalResult result;
    result.candidate = candidate;
    result.gram.setZero( workspace.gram.rows(), workspace.gram.cols() );
    if( scenario == BenchmarkScenario::sharedReference )
    {
        for( Eigen::Index column = 0; column < workspace.gram.cols(); ++column )
        {
            for( Eigen::Index row = column; row < workspace.gram.rows(); ++row )
            {
                result.gram( row, column ) = workspace.gram( row, column );
            }
        }
    }
    else
    {
        result.gram = workspace.gram.template cast<double>();
    }
    result.eigenvectors = workspace.eigenvectors.template cast<double>();
    if( scenario == BenchmarkScenario::sharedReference )
    {
        for( Eigen::Index mode = 0; mode < result.eigenvectors.cols(); ++mode )
        {
            result.eigenvectors.col( mode ) *= std::sqrt( static_cast<double>( workspace.eigenvalues( mode ) ) );
        }
    }
    result.eigenvalues = workspace.eigenvalues.topRows( workspace.eigenvectors.cols() ).template cast<double>();
    result.modes = workspace.modes.template cast<double>();
    result.residuals = workspace.residuals.template cast<double>();
    return result;
}

/// Return retained spatial modes as descending columns for projector calculations.
Eigen::MatrixXd retainedModeColumns( const matrixT<double> &modes /**< [in] ascending M-by-P spatial modes */ )
{
    Eigen::MatrixXd columns( modes.cols(), modes.rows() );
    for( Eigen::Index mode = 0; mode < modes.rows(); ++mode )
    {
        columns.col( mode ) = modes.matrix().row( modes.rows() - 1 - mode ).transpose();
    }
    return columns;
}

/// Return retained eigenvalues in descending order.
Eigen::VectorXd retainedEigenvalues( const matrixT<double> &eigenvalues /**< [in] ascending selected eigenvalues */ )
{
    Eigen::VectorXd retained( eigenvalues.rows() );
    for( Eigen::Index mode = 0; mode < eigenvalues.rows(); ++mode )
    {
        retained( mode ) = eigenvalues( eigenvalues.rows() - 1 - mode );
    }
    return retained;
}

/// Apply a matrix of descending orthonormal mode columns to every target in FP64.
Eigen::MatrixXd projectedResiduals( const Eigen::MatrixXd &modeColumns, /**< [in] P-by-M orthonormal modes */
                                    const Eigen::MatrixXd &targets )    /**< [in] P-by-T target vectors */
{
    return targets - modeColumns * ( modeColumns.transpose() * targets );
}

/// Apply FP32 modes to every target using the production KLIP worker accumulation order.
matrixT<float> productionOrderedResiduals( const matrixT<float> &modes, /**< [in] M-by-P production modes */
                                           const matrixT<float> &targets /**< [in] P-by-T FP32 target vectors */ )
{
    matrixT<float> residuals( targets.rows(), targets.cols() );
    matrixT<float> residual;
    matrixT<float> coefficients;
    matrixT<float> psf;
    for( Eigen::Index target = 0; target < targets.cols(); ++target )
    {
        subtractModes( residual, coefficients, psf, modes, targets.col( target ) );
        residuals.col( target ) = residual;
    }
    return residuals;
}

/// Construct the independent FP64 thin-SVD oracle once for both precision candidates.
OracleResult makeOracle( const SourceData &source, /**< [in] exact common FP32 inputs */
                         const BenchmarkOptions &options /**< [in] Gram branch and retained modes */ )
{
    const Eigen::MatrixXd references = source.references.matrix().cast<double>();
    const Eigen::MatrixXd targets = source.targets.matrix().cast<double>();
    OracleResult oracle;
    if( options.scenario == BenchmarkScenario::sharedReference )
    {
        oracle.gram = references.transpose() * references;
    }
    else
    {
        oracle.gram = references * references.transpose();
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> decomposition( references, Eigen::ComputeThinU );
    oracle.modeColumns = decomposition.matrixU().leftCols( options.modes );
    oracle.eigenvalues = decomposition.singularValues().head( options.modes ).array().square();
    oracle.residuals = projectedResiduals( oracle.modeColumns, targets );
    oracle.targetNorm = targets.norm();
    oracle.residualNorm = oracle.residuals.norm();
    return oracle;
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

/// Return compact retained-projector distance evidence without materializing either projector.
NormalizedMetric projectorMetric( const Eigen::MatrixXd &candidateModes,  /**< [in] candidate mode columns */
                                  const Eigen::MatrixXd &referenceModes ) /**< [in] reference mode columns */
{
    const Eigen::MatrixXd candidateGram = candidateModes.transpose() * candidateModes;
    const Eigen::MatrixXd referenceGram = referenceModes.transpose() * referenceModes;
    const Eigen::MatrixXd crossGram = candidateModes.transpose() * referenceModes;
    const long double distanceSquared = static_cast<long double>( candidateGram.squaredNorm() ) +
                                        static_cast<long double>( referenceGram.squaredNorm() ) -
                                        2.0L * static_cast<long double>( crossGram.squaredNorm() );
    NormalizedMetric metric;
    metric.numerator = std::sqrt( static_cast<double>( std::max( 0.0L, distanceSquared ) ) );
    metric.denominator = std::max( referenceGram.norm(), std::numeric_limits<double>::min() );
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

/// Return lower-triangle Gram evidence and maximum absolute error under the canonical KLIP rule.
NormalizedMetric
lowerTriangleMetric( const Eigen::MatrixXd &candidate, /**< [in] candidate Gram with valid lower data */
                     const Eigen::MatrixXd &reference, /**< [in] complete reference Gram matrix */
                     double &maximumError )            /**< [out] maximum absolute lower error */
{
    double differenceSquared{ 0 };
    double referenceSquared{ 0 };
    maximumError = 0;
    for( Eigen::Index column = 0; column < reference.cols(); ++column )
    {
        for( Eigen::Index row = column; row < reference.rows(); ++row )
        {
            const double difference = candidate( row, column ) - reference( row, column );
            differenceSquared += difference * difference;
            referenceSquared += reference( row, column ) * reference( row, column );
            maximumError = std::max( maximumError, std::abs( difference ) );
        }
    }
    NormalizedMetric metric;
    metric.numerator = std::sqrt( differenceSquared );
    metric.denominator = std::max( std::sqrt( referenceSquared ), std::numeric_limits<double>::min() );
    metric.value = metric.numerator / metric.denominator;
    return metric;
}

/// Calculate oracle and production-agreement metrics for one candidate.
template <typename eigensolverT>
ValidationMetrics validateResult( const NumericalResult &candidate, /**< [in] staged result to evaluate */
                                  const SourceData &source,         /**< [in] exact common FP32 inputs */
                                  const BenchmarkOptions &options,  /**< [in] dimensions and engineering gates */
                                  const OracleResult &oracle )      /**< [in] independent FP64 thin-SVD reference */
{
    const Eigen::MatrixXd candidateModes = retainedModeColumns( candidate.modes );
    const Eigen::VectorXd candidateEigenvalues = retainedEigenvalues( candidate.eigenvalues );
    const Eigen::MatrixXd candidateResiduals = candidate.residuals.matrix();

    matrixT<float> productionCovariance;
    if( options.scenario == BenchmarkScenario::sharedReference )
    {
        mx::math::eigenSYRK( productionCovariance, source.references );
    }
    matrixT<float> productionModes;
    mx::math::syevrMem<eigensolverT> productionWorkspace;
    const MXLAPACK_INT productionStatus = mx::improc::calcKLModesAdaptive<eigensolverT>( productionModes,
                                                                                         productionCovariance,
                                                                                         source.references,
                                                                                         options.modes,
                                                                                         &productionWorkspace );
    if( productionStatus != 0 )
    {
        throw std::runtime_error( candidate.candidate + " production calcKLModesAdaptive failed with status " +
                                  std::to_string( productionStatus ) );
    }
    const Eigen::MatrixXd productionModeColumns = retainedModeColumns( productionModes.template cast<double>() );
    const Eigen::MatrixXd productionResiduals =
        productionOrderedResiduals( productionModes, source.targets ).matrix().cast<double>();

    ValidationMetrics metrics;
    if( options.scenario == BenchmarkScenario::sharedReference )
    {
        metrics.gramEvidence = lowerTriangleMetric( candidate.gram.matrix(), oracle.gram, metrics.gramMaximumError );
    }
    else
    {
        const Eigen::MatrixXd gramDifference = candidate.gram.matrix() - oracle.gram;
        metrics.gramEvidence = normalizedMetric( gramDifference, oracle.gram.norm() );
        metrics.gramMaximumError = gramDifference.cwiseAbs().maxCoeff();
    }
    metrics.gramRelativeError = metrics.gramEvidence.value;
    metrics.eigenvalueEvidence =
        normalizedMetric( candidateEigenvalues - oracle.eigenvalues, oracle.eigenvalues.norm() );
    metrics.eigenvalueRelativeError = metrics.eigenvalueEvidence.value;
    metrics.projectorEvidence = projectorMetric( candidateModes, oracle.modeColumns );
    metrics.projectorRelativeError = metrics.projectorEvidence.value;
    metrics.maximumPrincipalAngle = maximumPrincipalAngle( candidateModes, oracle.modeColumns );
    const Eigen::MatrixXd residualDifference = candidateResiduals - oracle.residuals;
    metrics.residualOverOracleResidual = normalizedMetric( residualDifference, oracle.residualNorm );
    metrics.residualOverTarget = normalizedMetric( residualDifference, oracle.targetNorm );
    metrics.residualRelativeError = metrics.residualOverTarget.value;
    metrics.residualMaximumError = residualDifference.cwiseAbs().maxCoeff();
    metrics.productionProjectorEvidence = projectorMetric( candidateModes, productionModeColumns );
    metrics.productionProjectorRelativeError = metrics.productionProjectorEvidence.value;
    metrics.productionResidualEvidence =
        normalizedMetric( candidateResiduals - productionResiduals, oracle.targetNorm );
    metrics.productionResidualRelativeError = metrics.productionResidualEvidence.value;
    metrics.orthogonalityDefect =
        ( candidateModes.transpose() * candidateModes - Eigen::MatrixXd::Identity( options.modes, options.modes ) )
            .norm();

    Eigen::MatrixXd candidateGram = candidate.gram.matrix();
    if( options.scenario == BenchmarkScenario::sharedReference )
    {
        candidateGram = candidateGram.selfadjointView<Eigen::Lower>();
    }
    const Eigen::MatrixXd solverDifference =
        candidateGram * candidate.eigenvectors.matrix() -
        candidate.eigenvectors.matrix() * candidate.eigenvalues.matrix().asDiagonal();
    metrics.solverResidualEvidence =
        normalizedMetric( solverDifference, candidateGram.norm() * candidate.eigenvectors.matrix().norm() );

    metrics.passed = metrics.gramRelativeError <= options.gramRelativeTolerance &&
                     metrics.eigenvalueRelativeError <= options.eigenvalueRelativeTolerance &&
                     metrics.projectorRelativeError <= options.projectorRelativeTolerance &&
                     metrics.residualRelativeError <= options.residualRelativeTolerance &&
                     metrics.productionProjectorRelativeError <= options.productionRelativeTolerance &&
                     metrics.productionResidualRelativeError <= options.productionRelativeTolerance &&
                     metrics.solverResidualEvidence.value <= options.solverResidualTolerance &&
                     metrics.orthogonalityDefect <= options.orthogonalityTolerance;
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
                       const MetricRecord &metric,      /**< [in] gate-specific metric evidence */
                       int repetition,                  /**< [in] zero-based execution repetition */
                       int executionPosition )          /**< [in] balanced-order position, or negative */
{
    validateMetricRecord( metric );
    const bool referenceGram = options.scenario == BenchmarkScenario::sharedReference;
    const std::string_view eigensolverPrecision = candidate == "KLIP-F32" ? "fp32" : "fp64";
    const auto booleanOrNull = []( int value )
    { return value < 0 ? std::string( "null" ) : std::string( value == 0 ? "false" : "true" ); };
    const std::string nullReason = metric.applicable ? "null" : metric.notApplicableReason;
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
                  "klip",
                  std::to_string( repetition ),
                  std::to_string( executionPosition ),
                  "1",
                  std::to_string( options.innerBlasLapackThreads ),
                  std::to_string( options.seed ),
                  std::to_string( options.warmups ),
                  std::to_string( options.repetitions ),
                  "well-conditioned",
                  std::string( scenarioName( options.scenario ) ),
                  "null",
                  "null",
                  std::to_string( options.targets ),
                  std::to_string( options.pixels ),
                  std::to_string( options.references ),
                  referenceGram ? "reference_R_le_P" : "pixel_P_lt_R",
                  std::to_string( referenceGram ? options.references : options.pixels ),
                  std::to_string( referenceGram ? options.pixels : options.references ),
                  std::to_string( options.modes ),
                  std::to_string( options.modes ),
                  "fp32",
                  std::string( eigensolverPrecision ),
                  metric.metricName,
                  metric.numeratorName.empty() ? "null" : metric.numeratorName,
                  csvNumber( metric.numeratorValue ),
                  metric.numeratorNotApplicableReason.empty() ? "null" : metric.numeratorNotApplicableReason,
                  metric.denominatorName.empty() ? "null" : metric.denominatorName,
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
                  "uncentered_direct_in_memory_reused_buffers_no_transient_allocations_no_selection_no_openmp",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "null",
                  "not_applicable",
                  "null",
                  "not_applicable",
                  "not_applicable",
                  "not_applicable",
                  "false",
                  referenceGram ? "authoritative_lower_triangle_including_diagonal" : "full_symmetric",
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

/// Emit every available engineering validation record for one candidate.
void emitValidationRecords( std::ostream &output,             /**< [out] machine-readable destination */
                            const RunIdentity &identity,      /**< [in] verified run identities */
                            std::string_view candidate,       /**< [in] precision-policy identifier */
                            const ValidationMetrics &metrics, /**< [in] validation evidence */
                            const BenchmarkOptions &options ) /**< [in] dimensions and limits */
{
    const auto emit = [&]( const MetricRecord &record )
    { emitMetricRecord( output, identity, candidate, options, record, 0, -1 ); };
    emit( normalizedGate( "smoke.gram_relative_frobenius",
                          "fp64-direct-Gram-from-original-fp32-values",
                          "gram_relative_frobenius",
                          "norm_frobenius_candidate_minus_fp64_gram",
                          "protected_norm_fp64_gram_frobenius",
                          metrics.gramEvidence,
                          options.gramRelativeTolerance ) );
    emit( normalizedGate( "smoke.retained_eigenvalue_relative_l2",
                          "fp64_direct_thin_svd",
                          "retained_eigenvalue_relative_l2",
                          "norm_l2_candidate_minus_oracle_eigenvalues",
                          "protected_norm_oracle_retained_eigenvalues_l2",
                          metrics.eigenvalueEvidence,
                          options.eigenvalueRelativeTolerance ) );
    emit( normalizedGate( "smoke.retained_projector_relative_frobenius",
                          "fp64_direct_thin_svd",
                          "retained_projector_relative_frobenius",
                          "norm_frobenius_candidate_minus_oracle_projector",
                          "protected_norm_oracle_projector_frobenius",
                          metrics.projectorEvidence,
                          options.projectorRelativeTolerance ) );
    MetricRecord residualOracle = normalizedGate( "smoke.residual_over_oracle_residual",
                                                  "fp64_direct_thin_svd",
                                                  "residual_error_over_oracle_residual",
                                                  "norm_frobenius_candidate_minus_oracle_residual",
                                                  "protected_norm_oracle_residual_frobenius",
                                                  metrics.residualOverOracleResidual,
                                                  options.residualRelativeTolerance );
    if( metrics.residualOverOracleResidual.denominator == std::numeric_limits<double>::min() )
    {
        residualOracle.numeratorName.clear();
        residualOracle.numeratorValue = std::numeric_limits<double>::quiet_NaN();
        residualOracle.numeratorNotApplicableReason = "metric_inapplicable";
        residualOracle.denominatorName.clear();
        residualOracle.denominatorValue = std::numeric_limits<double>::quiet_NaN();
        residualOracle.denominatorNotApplicableReason = "oracle_residual_norm_is_zero";
        residualOracle.value = std::numeric_limits<double>::quiet_NaN();
        residualOracle.applicable = false;
        residualOracle.notApplicableReason = "oracle_residual_norm_is_zero";
        residualOracle.passed = -1;
    }
    emit( residualOracle );
    emit( normalizedGate( "smoke.residual_over_target",
                          "fp64_direct_thin_svd",
                          "residual_error_over_target",
                          "norm_frobenius_candidate_minus_oracle_residual",
                          "protected_norm_original_target_frobenius",
                          metrics.residualOverTarget,
                          options.residualRelativeTolerance ) );
    emit( normalizedGate( "smoke.normalized_solver_residual",
                          "intrinsic-candidate-eigensystem",
                          "normalized_solver_residual",
                          "norm_frobenius_gram_q_minus_q_lambda",
                          "protected_norm_gram_times_norm_q",
                          metrics.solverResidualEvidence,
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
    orthogonality.aggregationScope = "retained_spatial_modes";
    orthogonality.limitOperator = "less_than_or_equal";
    orthogonality.limitValue = options.orthogonalityTolerance;
    orthogonality.passed = metrics.orthogonalityDefect <= options.orthogonalityTolerance ? 1 : 0;
    emit( orthogonality );

    emit( normalizedGate( "smoke.production_mapping",
                          "production_api_same_policy",
                          "production_mapping_projector",
                          "norm_frobenius_staged_minus_production_projector",
                          "protected_norm_production_projector_frobenius",
                          metrics.productionProjectorEvidence,
                          options.productionRelativeTolerance ) );
    emit( normalizedGate( "smoke.production_mapping",
                          "production_api_same_policy",
                          "production_mapping_residual",
                          "norm_frobenius_staged_minus_production_residual",
                          "protected_norm_original_target_frobenius",
                          metrics.productionResidualEvidence,
                          options.productionRelativeTolerance ) );

    MetricRecord gramMaximum;
    gramMaximum.gateId = "smoke.gram_maximum_absolute";
    gramMaximum.oracleId = "fp64-direct-Gram-from-original-fp32-values";
    gramMaximum.metricName = "gram_maximum_absolute";
    gramMaximum.numeratorName = "maximum_absolute_gram_entry_error";
    gramMaximum.numeratorValue = metrics.gramMaximumError;
    gramMaximum.denominatorNotApplicableReason = "maximum_absolute_metric_has_no_denominator";
    gramMaximum.value = metrics.gramMaximumError;
    gramMaximum.units = "squared_input_units";
    gramMaximum.aggregationScope = options.scenario == BenchmarkScenario::sharedReference
                                       ? "authoritative_lower_triangle_including_diagonal"
                                       : "all_gram_entries";
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
    residualMaximum.aggregationScope = "all_targets_and_pixels";
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
    principalAngle.aggregationScope = "retained_spatial_subspace";
    principalAngle.passed = -1;
    principalAngle.gateStatus = "frozen-report-only";
    emit( principalAngle );
}

/// Emit five non-gating stage diagnostics without claiming the owner-TBD end-to-end performance gate.
void emitTimingRecords( std::ostream &output,            /**< [out] machine-readable destination */
                        const RunIdentity &identity,     /**< [in] verified run identities */
                        std::string_view candidate,      /**< [in] precision-policy identifier */
                        const StageTimings &timing,      /**< [in] four exclusive and one total timing */
                        const BenchmarkOptions &options, /**< [in] dimensions and controls */
                        int repetition,                  /**< [in] zero-based timed repetition */
                        int executionPosition )          /**< [in] balanced-order position */
{
    const std::array<std::string_view, 5> stageNames{ "gram",
                                                      "eigensolve",
                                                      "mode_construction",
                                                      "subtraction",
                                                      "total" };
    const std::array<double, 5> stageSeconds{ timing.gramSeconds,
                                              timing.eigensolveSeconds,
                                              timing.modeConstructionSeconds,
                                              timing.subtractionSeconds,
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
        emitMetricRecord( output, identity, candidate, options, record, repetition, executionPosition );
    }
}

/// Update one candidate's timing samples.
void accumulateSummary( CandidateSummary &summary, /**< [in,out] aggregate candidate result */
                        const StageTimings &timing /**< [in] current execution timings */ )
{
    const std::array<double, 5> seconds{ timing.gramSeconds,
                                         timing.eigensolveSeconds,
                                         timing.modeConstructionSeconds,
                                         timing.subtractionSeconds,
                                         timing.totalSeconds };
    for( std::size_t stage = 0; stage < seconds.size(); ++stage )
    {
        summary.stageSeconds[stage].push_back( seconds[stage] );
    }
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

/// Print concise median timings, validation metrics, gates, and thread controls.
void printHumanSummary( std::ostream &output,            /**< [out] human-readable destination */
                        const BenchmarkOptions &options, /**< [in] dimensions and declared gates */
                        const std::array<CandidateSummary, 2> &summaries /**< [in] aggregate candidate results */ )
{
    const std::array<std::string_view, 2> names{ "KLIP-M32D64", "KLIP-F32" };
    output << "KLIP CPU precision microbenchmark (well-conditioned, uncentered, direct): scenario="
           << scenarioName( options.scenario ) << ", T=" << options.targets << ", P=" << options.pixels
           << ", R=" << options.references << ", modes=" << options.modes << ", warmups=" << options.warmups
           << ", repetitions=" << options.repetitions << '\n';
    output << "Thread environment: OPENBLAS_NUM_THREADS=" << environmentValue( "OPENBLAS_NUM_THREADS" )
           << ", OMP_NUM_THREADS=" << environmentValue( "OMP_NUM_THREADS" ) << '\n';
    output << "Candidate order alternates once per repetition; "
           << ( options.repetitions % 2 == 0 ? "this run is balanced." : "this run ends unbalanced." ) << '\n';
    output << "Scope: arithmetic stages use warmed reusable matrices and one M-mode residual per target; transient "
              "per-call allocation, multi-mode-list output materialization, reference selection/collapse, and OpenMP "
              "scheduling are excluded.\n";
    output << "Median milliseconds:\n"
           << std::left << std::setw( 15 ) << "candidate" << std::right << std::setw( 12 ) << "gram" << std::setw( 12 )
           << "eigensolve" << std::setw( 14 ) << "modes" << std::setw( 14 ) << "subtraction" << std::setw( 12 )
           << "total" << '\n';
    output << std::fixed << std::setprecision( 3 );
    for( std::size_t candidate = 0; candidate < summaries.size(); ++candidate )
    {
        output << std::left << std::setw( 15 ) << names[candidate] << std::right;
        for( const auto &stage : summaries[candidate].stageSeconds )
        {
            output << std::setw( 12 ) << 1000.0 * median( stage );
        }
        output << '\n';
    }

    output << std::scientific << std::setprecision( 3 );
    output << "Single untimed validation (gram, eigenvalue, SVD projector, residual/oracle, residual/target, "
              "production projector, production residual, solver residual, orthogonality):\n";
    for( std::size_t candidate = 0; candidate < summaries.size(); ++candidate )
    {
        const ValidationMetrics &metrics = summaries[candidate].metrics;
        output << "  " << names[candidate] << ": " << metrics.gramRelativeError << ", "
               << metrics.eigenvalueRelativeError << ", " << metrics.projectorRelativeError << ", "
               << metrics.residualOverOracleResidual.value << ", " << metrics.residualOverTarget.value << ", "
               << metrics.productionProjectorRelativeError << ", " << metrics.productionResidualRelativeError << ", "
               << metrics.solverResidualEvidence.value << ", " << metrics.orthogonalityDefect << " ["
               << ( summaries[candidate].passed ? "PASS" : "FAIL" ) << "]\n";
    }
    output << "Engineering smoke gates (not final science tolerances): " << options.gramRelativeTolerance << ", "
           << options.eigenvalueRelativeTolerance << ", " << options.projectorRelativeTolerance << ", "
           << options.residualRelativeTolerance << ", " << options.productionRelativeTolerance << ", "
           << options.solverResidualTolerance << ", " << options.orthogonalityTolerance << '\n';
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

        const SourceData source = makeSourceData( options );
        const RunIdentity identity = makeRunIdentity( options, source, argv[0] );
        const OracleResult oracle = makeOracle( source, options );
        std::array<ValidationMetrics, 2> validationMetrics;
        {
            CandidateWorkspace<double> mixedValidationWorkspace;
            CandidateWorkspace<float> fp32ValidationWorkspace;
            executeCandidate( "KLIP-M32D64", source, options, mixedValidationWorkspace );
            executeCandidate( "KLIP-F32", source, options, fp32ValidationWorkspace );
            validationMetrics[0] = validateResult<double>(
                captureNumericalResult( "KLIP-M32D64", mixedValidationWorkspace, options.scenario ),
                source,
                options,
                oracle );
            validationMetrics[1] =
                validateResult<float>( captureNumericalResult( "KLIP-F32", fp32ValidationWorkspace, options.scenario ),
                                       source,
                                       options,
                                       oracle );
        }

        CandidateWorkspace<double> mixedWorkspace;
        CandidateWorkspace<float> fp32Workspace;
        for( int warmup = 0; warmup < options.warmups; ++warmup )
        {
            executeCandidate( "KLIP-M32D64", source, options, mixedWorkspace );
            executeCandidate( "KLIP-F32", source, options, fp32Workspace );
        }

        emitCsvHeader( *csvOutput );
        std::array<CandidateSummary, 2> summaries;
        for( std::size_t candidate = 0; candidate < summaries.size(); ++candidate )
        {
            summaries[candidate].metrics = validationMetrics[candidate];
            summaries[candidate].passed = validationMetrics[candidate].passed;
        }
        const std::array<std::string_view, 2> names{ "KLIP-M32D64", "KLIP-F32" };
        for( std::size_t candidate = 0; candidate < validationMetrics.size(); ++candidate )
        {
            emitValidationRecords( *csvOutput, identity, names[candidate], validationMetrics[candidate], options );
        }
        for( int repetition = 0; repetition < options.repetitions; ++repetition )
        {
            std::array<StageTimings, 2> timings;
            std::array<int, 2> executionPositions;
            for( int offset = 0; offset < 2; ++offset )
            {
                const int candidate = ( repetition + offset ) % 2;
                executionPositions[static_cast<std::size_t>( candidate )] = offset;
                if( candidate == 0 )
                {
                    timings[0] = executeCandidate( names[0], source, options, mixedWorkspace );
                }
                else
                {
                    timings[1] = executeCandidate( names[1], source, options, fp32Workspace );
                }
            }

            for( std::size_t candidate = 0; candidate < timings.size(); ++candidate )
            {
                emitTimingRecords( *csvOutput,
                                   identity,
                                   names[candidate],
                                   timings[candidate],
                                   options,
                                   repetition,
                                   executionPositions[candidate] );
                accumulateSummary( summaries[candidate], timings[candidate] );
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
            *summaryOutput << "Numerical validation failed; see CSV metrics and engineering smoke gates.\n";
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
        std::cerr << "klipPrecisionBenchmark: " << error.what() << '\n';
        return 1;
    }
}
