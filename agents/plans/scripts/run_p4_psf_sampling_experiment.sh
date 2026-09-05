#!/usr/bin/env bash

set -euo pipefail

export LC_ALL=C

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "${script_dir}/../../.." && pwd)
roc_working_dir="${repo_root}/working/roc"
base_config=${BASE_CONFIG:-"${roc_working_dir}/p4Reduce_afLepNaco.conf"}
psf_file=${PSF_FILE:-/home/jrmales/Source/mxWork/NACO/AFLep/2011-10-21/out/psf_reg_median.fits}
p4reduce_bin=${P4REDUCE_BIN:-p4Reduce}
psf_stamp_size=${PSF_STAMP_SIZE:-11}
psf_filter=${PSF_FILTER:-true}
experiment_dir=${EXPERIMENT_DIR:-"${roc_working_dir}/p4_psf_sampling_$(date -u +%Y%m%dT%H%M%SZ)"}
dry_run=false
analyze_only=false

all_cases=(
    dense
    radial_dr1_a16
    radial_dr2_a4
    radial_dr2_a8
    radial_dr2_a16
    radial_dr2_a32
    radial_dr4_a16
    radial_dr8_a16
)

usage()
{
    cat <<EOF
Usage: $(basename "$0") [options] [case ...]

Run the AF Lep/NACO P4 frozen-response sampling experiment. With no case
arguments, all cases are run sequentially. Completed cases are skipped.

Options:
  --list                 List the available cases and exit.
  --dry-run              Print commands without running p4Reduce.
  --analyze-only DIR     Analyze completed cases in DIR without reducing.
  -h, --help             Show this help.

Environment overrides:
  P4REDUCE_BIN           p4Reduce executable (default: p4Reduce from PATH)
  BASE_CONFIG            standard configuration (default: ${base_config})
  PSF_FILE               centered post-preprocessing PSF (default: ${psf_file})
  PSF_STAMP_SIZE         response stamp width (default: ${psf_stamp_size})
  PSF_FILTER             write filtered science products (default: ${psf_filter})
  EXPERIMENT_DIR         fixed output directory, useful when resuming
  OMP_NUM_THREADS        OpenMP worker limit passed through to p4Reduce

Examples:
  nohup $(basename "$0") > psf_sampling_driver.log 2>&1 &
  EXPERIMENT_DIR=/data/psf-test $(basename "$0") dense radial_dr2_a16
  $(basename "$0") --analyze-only /data/psf-test
EOF
}

list_cases()
{
    cat <<'EOF'
dense             exact response at every owned P4 search pixel
radial_dr1_a16    radii 0.5:1:59.5 pixels, 16 requested angles per radius
radial_dr2_a4     radii 1:2:59 pixels, 4 requested angles per radius
radial_dr2_a8     radii 1:2:59 pixels, 8 requested angles per radius
radial_dr2_a16    radii 1:2:59 pixels, 16 requested angles per radius
radial_dr2_a32    radii 1:2:59 pixels, 32 requested angles per radius
radial_dr4_a16    radii 2:4:58 pixels, 16 requested angles per radius
radial_dr8_a16    radii 4:8:52 plus 58 pixels, 16 requested angles per radius
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
    case "${case_name}" in
        dense)
            case_radii=
            case_angles=0
            ;;
        radial_dr1_a16)
            case_radii=$(radii_sequence 0.5 1 59.5)
            case_angles=16
            ;;
        radial_dr2_a4)
            case_radii=$(radii_sequence 1 2 59)
            case_angles=4
            ;;
        radial_dr2_a8)
            case_radii=$(radii_sequence 1 2 59)
            case_angles=8
            ;;
        radial_dr2_a16)
            case_radii=$(radii_sequence 1 2 59)
            case_angles=16
            ;;
        radial_dr2_a32)
            case_radii=$(radii_sequence 1 2 59)
            case_angles=32
            ;;
        radial_dr4_a16)
            case_radii=$(radii_sequence 2 4 58)
            case_angles=16
            ;;
        radial_dr8_a16)
            case_radii="$(radii_sequence 4 8 52),58"
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

manifest_complete()
{
    local manifest_path=$1
    [[ -f "${manifest_path}" ]] || return 1
    python3 - "${manifest_path}" <<'PY'
import sys
from astropy.io import fits

try:
    complete = int(fits.getheader(sys.argv[1]).get("P4 PSF COMPLETE", 0)) == 1
except Exception:
    complete = False
raise SystemExit(0 if complete else 1)
PY
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
    exec python3 "${script_dir}/compare_p4_psf_sampling.py" "${experiment_dir}"
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
command -v "${p4reduce_bin}" >/dev/null 2>&1 || {
    printf 'p4Reduce executable was not found: %s\n' "${p4reduce_bin}" >&2
    printf '%s\n' 'Set P4REDUCE_BIN to the newly built executable.' >&2
    exit 1
}

help_text=$("${p4reduce_bin}" --help 2>&1)
if [[ "${help_text}" != *"--p4.psfSampleRadii"* ]]; then
    printf 'p4Reduce does not expose --p4.psfSampleRadii: %s\n' "${p4reduce_bin}" >&2
    printf '%s\n' 'Build this hciReduce checkout and set P4REDUCE_BIN to that executable.' >&2
    exit 1
fi
if grep -Eq '^[[:space:]]*(psfSampleRadii|psfSamplesPerRadius)[[:space:]]*=' "${base_config}"; then
    printf 'The dense reference requires no active sparse-sampling keys in %s\n' "${base_config}" >&2
    printf '%s\n' 'Comment out those keys or point BASE_CONFIG at the unsampled standard configuration.' >&2
    exit 1
fi

for case_name in "${selected_cases[@]}"; do
    case_parameters "${case_name}"
done

mkdir -p "${experiment_dir}"
base_snapshot="${experiment_dir}/p4Reduce_afLepNaco.base.conf"
if [[ -e "${base_snapshot}" ]]; then
    if ! cmp -s "${base_config}" "${base_snapshot}"; then
        printf 'The base configuration differs from the snapshot in %s\n' "${experiment_dir}" >&2
        printf '%s\n' 'Choose a new EXPERIMENT_DIR to preserve experiment provenance.' >&2
        exit 1
    fi
else
    cp -- "${base_config}" "${base_snapshot}"
fi

binary_path=$(command -v "${p4reduce_bin}")
binary_sha256=$(sha256sum "${binary_path}" | awk '{ print $1 }')
psf_sha256=$(sha256sum "${psf_file}" | awk '{ print $1 }')
provenance_file="${experiment_dir}/provenance.txt"
if [[ -e "${provenance_file}" ]]; then
    recorded_binary_sha256=$(awk -F= '$1 == "p4reduce_sha256" { print $2 }' "${provenance_file}")
    recorded_psf_sha256=$(awk -F= '$1 == "psf_file_sha256" { print $2 }' "${provenance_file}")
    if [[ "${recorded_binary_sha256}" != "${binary_sha256}" || "${recorded_psf_sha256}" != "${psf_sha256}" ]]; then
        printf 'The binary or PSF differs from the provenance recorded in %s\n' "${experiment_dir}" >&2
        printf '%s\n' 'Choose a new EXPERIMENT_DIR rather than mixing experiment inputs.' >&2
        exit 1
    fi
else
    {
        printf 'created_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        printf 'hostname=%s\n' "$(hostname)"
        printf 'p4reduce_path=%s\n' "${binary_path}"
        printf 'p4reduce_sha256=%s\n' "${binary_sha256}"
        printf 'base_config=%s\n' "${base_config}"
        printf 'base_config_sha256=%s\n' "$(sha256sum "${base_config}" | awk '{ print $1 }')"
        printf 'psf_file=%s\n' "${psf_file}"
        printf 'psf_file_sha256=%s\n' "${psf_sha256}"
        printf 'hciReduce_commit=%s\n' "$(git -C "${repo_root}" rev-parse HEAD 2>/dev/null || printf unknown)"
    } > "${provenance_file}"
fi

summary_file="${experiment_dir}/runs.tsv"
if [[ ! -e "${summary_file}" ]]; then
    printf 'case\tstatus\twall_seconds\tstarted_utc\tfinished_utc\n' > "${summary_file}"
fi

printf 'Experiment directory: %s\n' "${experiment_dir}"
printf 'Base configuration: %s\n' "${base_config}"
printf 'p4Reduce: %s\n' "${binary_path}"
printf 'PSF template: %s\n' "${psf_file}"
printf 'OMP_NUM_THREADS: %s\n' "${OMP_NUM_THREADS:-unlimited}"

for case_name in "${selected_cases[@]}"; do
    case_parameters "${case_name}"
    case_dir="${experiment_dir}/${case_name}"
    manifest="${case_dir}/finim_outputs/p4PSF_manifest.fits"

    if manifest_complete "${manifest}"; then
        printf '\n[%s] complete manifest exists; skipping.\n' "${case_name}"
        continue
    fi
    if [[ -e "${case_dir}/run.log" || -e "${case_dir}/finim.fits" ]]; then
        printf '\n[%s] has incomplete prior output in %s; refusing to overwrite it.\n' "${case_name}" "${case_dir}" >&2
        printf '%s\n' 'Move that case directory aside or choose a new EXPERIMENT_DIR.' >&2
        exit 1
    fi

    mkdir -p "${case_dir}"
    command_line=(
        "${p4reduce_bin}"
        --config "${base_config}"
        --p4.psfFile "${psf_file}"
        --p4.psfStampSize "${psf_stamp_size}"
        --p4.outputPSFModels=true
        --p4.psfFilter="${psf_filter}"
        --p4.psfOutputPrefix p4PSF_
        --output.directory "${case_dir}"
        --output.fileName finim.fits
        --output.exactFName=true
        --showTiming=true
    )
    if [[ -n "${case_radii}" ]]; then
        command_line+=(
            --p4.psfSampleRadii "${case_radii}"
            --p4.psfSamplesPerRadius "${case_angles}"
        )
    fi

    shell_join "${command_line[@]}" > "${case_dir}/command.txt"
    cp -- "${base_config}" "${case_dir}/p4Reduce_afLepNaco.base.conf"
    {
        printf 'case=%s\n' "${case_name}"
        printf 'sample_radii=%s\n' "${case_radii}"
        printf 'samples_per_radius=%s\n' "${case_angles}"
        printf 'psf_file=%s\n' "${psf_file}"
        printf 'psf_stamp_size=%s\n' "${psf_stamp_size}"
        printf 'psf_filter=%s\n' "${psf_filter}"
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
        printf '[%s] p4Reduce failed with status %d.\n' "${case_name}" "${command_status}" >&2
        exit "${command_status}"
    fi
    if [[ ! -e "${manifest}" ]]; then
        printf '[%s] p4Reduce exited successfully but did not publish %s\n' "${case_name}" "${manifest}" >&2
        exit 1
    fi
done

if [[ "${dry_run}" == true ]]; then
    printf '\nDry run complete; no reductions were started.\n'
    exit 0
fi

has_sparse_case=false
for case_name in "${all_cases[@]:1}"; do
    if manifest_complete "${experiment_dir}/${case_name}/finim_outputs/p4PSF_manifest.fits"; then
        has_sparse_case=true
        break
    fi
done
if manifest_complete "${experiment_dir}/dense/finim_outputs/p4PSF_manifest.fits" && [[ "${has_sparse_case}" == true ]]; then
    python3 "${script_dir}/compare_p4_psf_sampling.py" "${experiment_dir}"
else
    printf '\nA completed dense reference and sparse case are not both present yet; skipping comparison.\n'
fi
