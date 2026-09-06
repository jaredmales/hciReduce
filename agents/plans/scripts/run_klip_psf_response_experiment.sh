#!/usr/bin/env bash

set -euo pipefail

export LC_ALL=C

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "${script_dir}/../../.." && pwd)
roc_working_dir="${repo_root}/working/roc"
base_config=${BASE_CONFIG:-"${script_dir}/klipReduce_afLepNaco_psf_response.conf"}
psf_file=${PSF_FILE:-/home/jrmales/Source/mxWork/NACO/AFLep/2011-10-21/out/psf_reg_median.fits}
klipreduce_bin=${KLIPREDUCE_BIN:-klipReduce}
psf_stamp_size=${PSF_STAMP_SIZE:-11}
experiment_dir=${EXPERIMENT_DIR:-"${roc_working_dir}/klip_psf_response_$(date -u +%Y%m%dT%H%M%SZ)"}
dry_run=false
analyze_only=false

all_cases=(
    science_only
    reference_dr1_a16
    radial_dr2_a4
    radial_dr2_a8
    radial_dr2_a16
    radial_dr4_a16
    radial_dr8_a16
)

usage()
{
    cat <<EOF
Usage: $(basename "$0") [options] [case ...]

Run the AF Lep/NACO KLIP sparse frozen-basis response experiment. With no case
arguments, all cases are run sequentially. Completed cases are skipped.

Options:
  --list                 List the available cases and exit.
  --dry-run              Validate inputs and print commands without reducing.
  --analyze-only DIR     Analyze completed cases in DIR without reducing.
  -h, --help             Show this help.

Environment overrides:
  KLIPREDUCE_BIN         klipReduce executable (default: klipReduce from PATH)
  BASE_CONFIG            response-compatible configuration (default: ${base_config})
  PSF_FILE               centered post-preprocessing PSF (default: ${psf_file})
  PSF_STAMP_SIZE         response stamp width (default: ${psf_stamp_size})
  EXPERIMENT_DIR         fixed output directory, useful when resuming
  OMP_NUM_THREADS        OpenMP worker limit passed through to klipReduce

Examples:
  OMP_NUM_THREADS=32 nohup $(basename "$0") > klip_psf_driver.log 2>&1 &
  EXPERIMENT_DIR=/data/klip-psf $(basename "$0") science_only reference_dr1_a16 radial_dr2_a16
  $(basename "$0") --analyze-only /data/klip-psf
EOF
}

list_cases()
{
    cat <<'EOF'
science_only       normal KLIP reduction without analytic response measurement
reference_dr1_a16  radii 6.5:1:59.5 pixels, 16 angles per radius
radial_dr2_a4      radii 7:2:59 pixels, 4 angles per radius
radial_dr2_a8      radii 7:2:59 pixels, 8 angles per radius
radial_dr2_a16     radii 7:2:59 pixels, 16 angles per radius
radial_dr4_a16     radii 8:4:56 plus 59 pixels, 16 angles per radius
radial_dr8_a16     radii 10:8:58 pixels, 16 angles per radius
EOF
}

radii_sequence()
{
    local first_radius=$1
    local step=$2
    local final_radius=$3
    seq -s, "${first_radius}" "${step}" "${final_radius}"
}

case_parameters()
{
    local case_name=$1
    case_response=true
    case_radii=
    case_angles=0
    case "${case_name}" in
        science_only)
            case_response=false
            ;;
        reference_dr1_a16)
            case_radii=$(radii_sequence 6.5 1 59.5)
            case_angles=16
            ;;
        radial_dr2_a4)
            case_radii=$(radii_sequence 7 2 59)
            case_angles=4
            ;;
        radial_dr2_a8)
            case_radii=$(radii_sequence 7 2 59)
            case_angles=8
            ;;
        radial_dr2_a16)
            case_radii=$(radii_sequence 7 2 59)
            case_angles=16
            ;;
        radial_dr4_a16)
            case_radii="$(radii_sequence 8 4 56),59"
            case_angles=16
            ;;
        radial_dr8_a16)
            case_radii=$(radii_sequence 10 8 58)
            case_angles=16
            ;;
        *)
            printf 'Unknown experiment case: %s\n' "${case_name}" >&2
            list_cases >&2
            exit 2
            ;;
    esac
}

shell_join()
{
    local argument
    for argument in "$@"; do
        printf '%q ' "${argument}"
    done
    printf '\n'
}

case_complete()
{
    local case_directory=$1
    [[ -f "${case_directory}/complete" && -f "${case_directory}/finim.fits" ]]
}

selected_cases=()
while (($#)); do
    case "$1" in
        --list)
            list_cases
            exit 0
            ;;
        --dry-run)
            dry_run=true
            shift
            ;;
        --analyze-only)
            if (($# < 2)); then
                printf '%s\n' '--analyze-only requires an experiment directory.' >&2
                exit 2
            fi
            analyze_only=true
            experiment_dir=$2
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --*)
            printf 'Unknown option: %s\n' "$1" >&2
            usage >&2
            exit 2
            ;;
        *)
            selected_cases+=("$1")
            shift
            ;;
    esac
done

if [[ "${analyze_only}" == true ]]; then
    exec python3 "${script_dir}/compare_klip_psf_response.py" "${experiment_dir}"
fi

if ((${#selected_cases[@]} == 0)); then
    selected_cases=("${all_cases[@]}")
fi

[[ -r "${base_config}" ]] || { printf 'Base configuration is not readable: %s\n' "${base_config}" >&2; exit 1; }
[[ -r "${psf_file}" ]] || { printf 'PSF template is not readable: %s\n' "${psf_file}" >&2; exit 1; }
[[ -x /usr/bin/time ]] || { printf '%s\n' '/usr/bin/time is required.' >&2; exit 1; }
command -v python3 >/dev/null 2>&1 || { printf '%s\n' 'python3 is required for result comparison.' >&2; exit 1; }
python3 -c 'import astropy.io.fits, numpy' >/dev/null 2>&1 || {
    printf '%s\n' 'The Python astropy and numpy packages are required for result comparison.' >&2
    exit 1
}
command -v "${klipreduce_bin}" >/dev/null 2>&1 || {
    printf 'klipReduce executable was not found: %s\n' "${klipreduce_bin}" >&2
    printf '%s\n' 'Set KLIPREDUCE_BIN to the newly built executable.' >&2
    exit 1
}

help_text=$("${klipreduce_bin}" --help 2>&1)
if [[ "${help_text}" != *"--klip.psfSampleRadii"* || "${help_text}" != *"--klip.outputPSFModels"* ]]; then
    printf 'klipReduce does not expose the sparse response options: %s\n' "${klipreduce_bin}" >&2
    printf '%s\n' 'Build this hciReduce checkout and set KLIPREDUCE_BIN to that executable.' >&2
    exit 1
fi
if grep -Eq '^[[:space:]]*(psfFile|psfStampSize|psfSampleRadii|psfSamplesPerRadius|outputPSFModels)[[:space:]]*=' \
    "${base_config}"; then
    printf 'The base configuration must not enable KLIP PSF measurement: %s\n' "${base_config}" >&2
    exit 1
fi

for case_name in "${selected_cases[@]}"; do
    case_parameters "${case_name}"
done

mkdir -p "${experiment_dir}"
base_snapshot="${experiment_dir}/klipReduce_afLepNaco_psf_response.base.conf"
if [[ -e "${base_snapshot}" ]]; then
    if ! cmp -s "${base_config}" "${base_snapshot}"; then
        printf 'The base configuration differs from the snapshot in %s\n' "${experiment_dir}" >&2
        printf '%s\n' 'Choose a new EXPERIMENT_DIR to preserve experiment provenance.' >&2
        exit 1
    fi
else
    cp -- "${base_config}" "${base_snapshot}"
fi

binary_path=$(command -v "${klipreduce_bin}")
binary_sha256=$(sha256sum "${binary_path}" | awk '{ print $1 }')
psf_sha256=$(sha256sum "${psf_file}" | awk '{ print $1 }')
provenance_file="${experiment_dir}/provenance.txt"
if [[ -e "${provenance_file}" ]]; then
    recorded_binary_sha256=$(awk -F= '$1 == "klipreduce_sha256" { print $2 }' "${provenance_file}")
    recorded_psf_sha256=$(awk -F= '$1 == "psf_file_sha256" { print $2 }' "${provenance_file}")
    recorded_stamp_size=$(awk -F= '$1 == "psf_stamp_size" { print $2 }' "${provenance_file}")
    recorded_workers=$(awk -F= '$1 == "omp_num_threads" { print $2 }' "${provenance_file}")
    if [[ "${recorded_binary_sha256}" != "${binary_sha256}" || "${recorded_psf_sha256}" != "${psf_sha256}" ||
          "${recorded_stamp_size}" != "${psf_stamp_size}" ||
          "${recorded_workers}" != "${OMP_NUM_THREADS:-unlimited}" ]]; then
        printf 'The binary, PSF, stamp size, or worker count differs from the provenance in %s\n' \
            "${experiment_dir}" >&2
        printf '%s\n' 'Choose a new EXPERIMENT_DIR rather than mixing experiment inputs.' >&2
        exit 1
    fi
else
    {
        printf 'created_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        printf 'hostname=%s\n' "$(hostname)"
        printf 'klipreduce_path=%s\n' "${binary_path}"
        printf 'klipreduce_sha256=%s\n' "${binary_sha256}"
        printf 'base_config=%s\n' "${base_config}"
        printf 'base_config_sha256=%s\n' "$(sha256sum "${base_config}" | awk '{ print $1 }')"
        printf 'psf_file=%s\n' "${psf_file}"
        printf 'psf_file_sha256=%s\n' "${psf_sha256}"
        printf 'psf_stamp_size=%s\n' "${psf_stamp_size}"
        printf 'omp_num_threads=%s\n' "${OMP_NUM_THREADS:-unlimited}"
        printf 'hciReduce_commit=%s\n' "$(git -C "${repo_root}" rev-parse HEAD 2>/dev/null || printf unknown)"
    } > "${provenance_file}"
fi

summary_file="${experiment_dir}/runs.tsv"
if [[ ! -e "${summary_file}" ]]; then
    printf 'case\tstatus\twall_seconds\tstarted_utc\tfinished_utc\n' > "${summary_file}"
fi

printf 'Experiment directory: %s\n' "${experiment_dir}"
printf 'Base configuration: %s\n' "${base_config}"
printf 'klipReduce: %s\n' "${binary_path}"
printf 'PSF template: %s\n' "${psf_file}"
printf 'OMP_NUM_THREADS: %s\n' "${OMP_NUM_THREADS:-unlimited}"

for case_name in "${selected_cases[@]}"; do
    case_parameters "${case_name}"
    case_dir="${experiment_dir}/${case_name}"

    if case_complete "${case_dir}"; then
        printf '\n[%s] completion marker and final image exist; skipping.\n' "${case_name}"
        continue
    fi
    if [[ -e "${case_dir}/run.log" || -e "${case_dir}/finim.fits" ]]; then
        printf '\n[%s] has incomplete prior output in %s; refusing to overwrite it.\n' "${case_name}" "${case_dir}" >&2
        printf '%s\n' 'Move that case directory aside or choose a new EXPERIMENT_DIR.' >&2
        exit 1
    fi

    mkdir -p "${case_dir}"
    command_line=(
        "${klipreduce_bin}"
        --config "${base_config}"
        --output.directory "${case_dir}"
        --output.fileName finim.fits
        --output.exactFName=true
        --showTiming=true
    )
    if [[ "${case_response}" == true ]]; then
        command_line+=(
            --klip.psfFile "${psf_file}"
            --klip.psfStampSize "${psf_stamp_size}"
            --klip.psfSampleRadii "${case_radii}"
            --klip.psfSamplesPerRadius "${case_angles}"
            --klip.outputPSFModels=true
            --klip.psfOutputPrefix klipPSF_
        )
    fi

    shell_join "${command_line[@]}" > "${case_dir}/command.txt"
    cp -- "${base_config}" "${case_dir}/klipReduce_afLepNaco_psf_response.base.conf"
    {
        printf 'case=%s\n' "${case_name}"
        printf 'response_enabled=%s\n' "${case_response}"
        printf 'sample_radii=%s\n' "${case_radii}"
        printf 'samples_per_radius=%s\n' "${case_angles}"
        printf 'psf_file=%s\n' "${psf_file}"
        printf 'psf_stamp_size=%s\n' "${psf_stamp_size}"
        printf 'omp_num_threads=%s\n' "${OMP_NUM_THREADS:-unlimited}"
    } > "${case_dir}/case.env"

    printf '\n[%s]\n' "${case_name}"
    shell_join "${command_line[@]}"
    if [[ "${dry_run}" == true ]]; then
        continue
    fi

    started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
    set +e
    /usr/bin/time -f 'wall_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\nmaximum_rss_kib=%M' \
        -o "${case_dir}/resource_usage.txt" \
        "${command_line[@]}" 2>&1 | tee "${case_dir}/run.log"
    command_status=${PIPESTATUS[0]}
    set -e
    finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
    wall_seconds=$(awk -F= '$1 == "wall_seconds" { print $2 }' "${case_dir}/resource_usage.txt")
    printf '%s\t%s\t%s\t%s\t%s\n' \
        "${case_name}" "${command_status}" "${wall_seconds:-unknown}" "${started_utc}" "${finished_utc}" \
        >> "${summary_file}"

    if ((command_status != 0)); then
        printf '[%s] klipReduce failed with status %d.\n' "${case_name}" "${command_status}" >&2
        exit "${command_status}"
    fi
    if [[ ! -f "${case_dir}/finim.fits" ]]; then
        printf '[%s] klipReduce did not publish its final image.\n' "${case_name}" >&2
        exit 1
    fi
    if [[ "${case_response}" == true ]]; then
        python3 - "${case_dir}" <<'PY'
import sys
from pathlib import Path

from astropy.io import fits

case_directory = Path(sys.argv[1])
product_directory = case_directory / "finim_outputs"
responses = sorted(product_directory.glob("klipPSF_mode*_radial_response.fits"))
validities = sorted(product_directory.glob("klipPSF_mode*_radial_validity.fits"))
if not responses or len(responses) != len(validities):
    raise SystemExit(f"incomplete KLIP response products in {product_directory}")
for response in responses:
    header = fits.getheader(response)
    if int(header.get("KLIP PSF PRODUCT SCHEMA", 0)) != 1:
        raise SystemExit(f"unexpected KLIP response schema in {response}")
    if str(header.get("KLIP PSF ACCUMULATION", "")).strip() != "WORKER_SUM":
        raise SystemExit(f"expected bounded worker-sum accumulation in {response}")
PY
    fi
    printf 'completed_utc=%s\n' "${finished_utc}" > "${case_dir}/complete"
done

if [[ "${dry_run}" == true ]]; then
    printf '\nDry run complete; no reductions were started.\n'
    exit 0
fi

has_sparse_case=false
for case_name in "${all_cases[@]:2}"; do
    if case_complete "${experiment_dir}/${case_name}"; then
        has_sparse_case=true
        break
    fi
done
if case_complete "${experiment_dir}/reference_dr1_a16" && [[ "${has_sparse_case}" == true ]]; then
    python3 "${script_dir}/compare_klip_psf_response.py" "${experiment_dir}"
else
    printf '\nA completed response reference and sparse case are not both present yet; skipping comparison.\n'
fi
