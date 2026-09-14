#!/usr/bin/env bash

set -euo pipefail

export LC_ALL=C
export PYTHONDONTWRITEBYTECODE=1

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "${script_dir}/../../.." && pwd)
roc_working_dir="${repo_root}/working/roc"
source_experiment=${SOURCE_EXPERIMENT:-"${roc_working_dir}/klip_cpp_response_20260913T222913Z"}
science_file=${SCIENCE_FILE:-"${source_experiment}/science_only/finim.fits"}
default_sparse_manifest="${source_experiment}/radial_ld_refit4_filter/finim_outputs/klipPSF_manifest.fits"
sparse_manifest=${SPARSE_RESPONSE_MANIFEST:-"${default_sparse_manifest}"}
experiment_dir=${EXPERIMENT_DIR:-"${roc_working_dir}/klip_exact_response_$(date -u +%Y%m%dT%H%M%SZ)"}
base_config=${BASE_CONFIG:-"${script_dir}/klipReduce_afLepNaco_psf_response.conf"}
psf_file=${PSF_FILE:-/home/jrmales/Source/mxWork/NACO/AFLep/2011-10-21/out/psf_reg_median.fits}
klipreduce_bin=${KLIPREDUCE_BIN:-klipReduce}
hcianalyze_bin=${HCIANALYZE_BIN:-hciAnalyze}
mode_counts=${MODE_COUNTS:-125,150,175,200,225,250,300,350}
planet_sep=${PLANET_SEP:-11.73795339222688}
planet_pa=${PLANET_PA:-262.1667995998323}
planet_contrast=${PLANET_CONTRAST:-0.004763925929356391}
response_fraction=${RESPONSE_FRACTION:-0.25}
psf_stamp_size=${PSF_STAMP_SIZE:-11}
psf_filter_min_good_fract=${PSF_FILTER_MIN_GOOD_FRACT:-1}
gaussian_fwhms=${GAUSSIAN_FWHMS:-3.6}
lambda_d=${LAMBDA_D:-3.6}
planet_radius=${PLANET_RADIUS:-5}
snr_min_radius=${SNR_MIN_RADIUS:-6}
snr_max_radius=${SNR_MAX_RADIUS:-60}
snr_aperture_radius=${SNR_APERTURE_RADIUS:-2}
dry_run=false

usage()
{
    cat <<EOF
Usage: $(basename "$0") [--dry-run]

Measure an exact local KLIP response at the nearest integer pixel to the known
candidate. The response is the paired end-to-end difference

  [KLIP(data + epsilon PSF) - KLIP(data - epsilon PSF)] / (2 epsilon).

Only the two fixed-location KLIP reductions are run; no negative-planet
position or contrast optimization is performed. The result is packaged as an
EXACT_AZIMUTHAL manifest and compared in hciAnalyze with the unfiltered,
sparse-response, and Gaussian-smoothed cases.

Environment overrides:
  SOURCE_EXPERIMENT          completed native sparse-response experiment
  SCIENCE_FILE               standalone unfiltered KLIP cube
  SPARSE_RESPONSE_MANIFEST   compatible sparse KLIP response manifest
  EXPERIMENT_DIR             output directory (default: ${experiment_dir})
  BASE_CONFIG                maintained AF Lep/NACO KLIP configuration
  PSF_FILE                   centered perturbation PSF
  KLIPREDUCE_BIN             klipReduce executable (default: ${klipreduce_bin})
  HCIANALYZE_BIN             hciAnalyze executable (default: ${hcianalyze_bin})
  MODE_COUNTS                exact KL mode vector (default: ${mode_counts})
  PLANET_SEP, PLANET_PA      known candidate coordinates
  PLANET_CONTRAST            known candidate contrast
  RESPONSE_FRACTION          epsilon/PLANET_CONTRAST (default: ${response_fraction})
  PSF_STAMP_SIZE             odd exact-response stamp width (default: ${psf_stamp_size})
  PSF_FILTER_MIN_GOOD_FRACT  minimum hciAnalyze filter support
  GAUSSIAN_FWHMS             comma-separated comparison widths
  LAMBDA_D                   pixels per lambda/D
  PLANET_RADIUS              SNR noise-exclusion radius
  SNR_MIN_RADIUS             inner radial-noise bound
  SNR_MAX_RADIUS             outer radial-noise bound
  SNR_APERTURE_RADIUS        signal aperture radius
  OMP_NUM_THREADS            OpenMP worker limit passed through to klipReduce

Example:
  OMP_NUM_THREADS=48 nohup $(basename "$0") > klip_exact_response_driver.log 2>&1 &
EOF
}

while (($#)); do
    case "$1" in
        --dry-run)
            dry_run=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            printf 'Unknown option: %s\n' "$1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

shell_join()
{
    local argument
    for argument in "$@"; do
        printf '%q ' "${argument}"
    done
    printf '\n'
}

absolute_file()
{
    local input_file=$1
    local input_directory
    [[ -f "${input_file}" ]] || {
        printf 'Required file is not readable: %s\n' "${input_file}" >&2
        exit 1
    }
    input_directory=$(cd -- "$(dirname -- "${input_file}")" && pwd)
    printf '%s/%s\n' "${input_directory}" "$(basename -- "${input_file}")"
}

run_timed()
{
    local stage_directory=$1
    shift
    local command_line=("$@")
    printf '\n[%s]\n' "${stage_directory#"${experiment_dir}/"}"
    shell_join "${command_line[@]}"
    if [[ "${dry_run}" == true ]]; then
        return
    fi
    mkdir -p "${stage_directory}"
    shell_join "${command_line[@]}" > "${stage_directory}/command.txt"
    set +e
    /usr/bin/time -f 'wall_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\nmaximum_rss_kib=%M' \
        -o "${stage_directory}/resource_usage.txt" \
        "${command_line[@]}" 2>&1 | tee "${stage_directory}/run.log"
    local command_status=${PIPESTATUS[0]}
    set -e
    if ((command_status != 0)); then
        printf '[%s] failed with status %d.\n' "${stage_directory#"${experiment_dir}/"}" \
            "${command_status}" >&2
        exit "${command_status}"
    fi
    [[ -s "${stage_directory}/finim.fits" ]] || {
        printf '[%s] did not publish finim.fits.\n' "${stage_directory#"${experiment_dir}/"}" >&2
        exit 1
    }
    : > "${stage_directory}/complete"
}

command -v "${klipreduce_bin}" >/dev/null 2>&1 || {
    printf 'klipReduce executable was not found: %s\n' "${klipreduce_bin}" >&2
    exit 1
}
command -v "${hcianalyze_bin}" >/dev/null 2>&1 || {
    printf 'hciAnalyze executable was not found: %s\n' "${hcianalyze_bin}" >&2
    exit 1
}
science_file=$(absolute_file "${science_file}")
sparse_manifest=$(absolute_file "${sparse_manifest}")
base_config=$(absolute_file "${base_config}")
psf_file=$(absolute_file "${psf_file}")

if [[ -e "${experiment_dir}" ]]; then
    printf 'Experiment output already exists; refusing to overwrite: %s\n' "${experiment_dir}" >&2
    exit 1
fi

mapfile -t exact_values < <(
    python3 - "${science_file}" "${planet_sep}" "${planet_pa}" "${planet_contrast}" "${response_fraction}" <<'PY'
import math
import sys

from astropy.io import fits

shape = fits.getdata(sys.argv[1], memmap=True).shape
if len(shape) != 3:
    raise SystemExit("SCIENCE_FILE must contain a three-dimensional KLIP cube")
separation, pa, contrast, fraction = map(float, sys.argv[2:])
if not all(math.isfinite(value) for value in (separation, pa, contrast, fraction)):
    raise SystemExit("planet geometry, contrast, and response fraction must be finite")
if separation < 0 or contrast <= 0 or fraction <= 0:
    raise SystemExit("planet separation, contrast, and response fraction are invalid")
columns, rows = shape[1:]
center_row = 0.5 * (rows - 1)
center_column = 0.5 * (columns - 1)
pa_radians = math.radians(pa)
requested_row = center_row - separation * math.sin(pa_radians)
requested_column = center_column + separation * math.cos(pa_radians)
target_row = math.floor(requested_row + 0.5)
target_column = math.floor(requested_column + 0.5)
delta_row = target_row - center_row
delta_column = target_column - center_column
exact_separation = math.hypot(delta_row, delta_column)
exact_pa = math.degrees(math.atan2(-delta_row, delta_column)) % 360
epsilon = contrast * fraction
for value in (target_row, target_column, exact_separation, exact_pa, epsilon):
    print(f"{value:.17g}")
PY
)
if ((${#exact_values[@]} != 5)); then
    printf 'Could not resolve the exact detector-pixel response coordinate.\n' >&2
    exit 1
fi
exact_row=${exact_values[0]}
exact_column=${exact_values[1]}
exact_sep=${exact_values[2]}
exact_pa=${exact_values[3]}
response_contrast=${exact_values[4]}

printf 'Nearest exact candidate pixel: row=%s column=%s separation=%s PA=%s\n' \
    "${exact_row}" "${exact_column}" "${exact_sep}" "${exact_pa}"
printf 'Paired perturbation half-amplitude: %s\n' "${response_contrast}"

if [[ "${dry_run}" == false ]]; then
    mkdir -p "${experiment_dir}"
    printf '%s\n' \
        "science_file=${science_file}" \
        "sparse_manifest=${sparse_manifest}" \
        "mode_counts=${mode_counts}" \
        "planet_sep=${planet_sep}" \
        "planet_pa=${planet_pa}" \
        "planet_contrast=${planet_contrast}" \
        "response_fraction=${response_fraction}" \
        "response_contrast=${response_contrast}" \
        "exact_row=${exact_row}" \
        "exact_column=${exact_column}" \
        "exact_sep=${exact_sep}" \
        "exact_pa=${exact_pa}" \
        > "${experiment_dir}/settings.env"
fi

for sign in plus minus; do
    signed_contrast=${response_contrast}
    if [[ "${sign}" == minus ]]; then
        signed_contrast="-${response_contrast}"
    fi
    command_line=(
        "${klipreduce_bin}"
        --config "${base_config}"
        --klip.Nmodes "${mode_counts}"
        --psfResponse.file ""
        --psfResponse.outputModels=false
        --psfResponse.filter=false
        --planet.sep "${planet_sep}"
        --planet.PA "${planet_pa}"
        --planet.contrast "${planet_contrast}"
        --fake.method single
        --fake.fileName "${psf_file}"
        --fake.sep "${exact_sep}"
        --fake.PA "${exact_pa}"
        --fake.contrast "${signed_contrast}"
        --fake.subtractPlanet=false
        --output.directory "${experiment_dir}/runs/${sign}"
        --output.fileName finim.fits
        --output.exactFName=true
        --showTiming=true
    )
    run_timed "${experiment_dir}/runs/${sign}" "${command_line[@]}"
done

product_directory="${experiment_dir}/exact_response"
build_command=(
    python3 "${script_dir}/build_klip_exact_response.py"
    "${science_file}"
    "${experiment_dir}/runs/plus/finim.fits"
    "${experiment_dir}/runs/minus/finim.fits"
    "${sparse_manifest}"
    "${product_directory}"
    --separation "${exact_sep}"
    --pa "${exact_pa}"
    --perturbation "${response_contrast}"
    --stamp-size "${psf_stamp_size}"
    --minimum-support "${psf_filter_min_good_fract}"
)
printf '\n[build exact response]\n'
shell_join "${build_command[@]}"
if [[ "${dry_run}" == true ]]; then
    printf '\n[hciAnalyze comparison]\n'
    printf 'EXACT_RESPONSE_MANIFEST=%q %q\n' \
        "${product_directory}/klipExact_manifest.fits" \
        "${script_dir}/run_klip_hciAnalyze_filter_comparison.sh"
    printf '\nDry run complete; no KLIP reductions were started.\n'
    exit 0
fi
shell_join "${build_command[@]}" > "${experiment_dir}/build_command.txt"
"${build_command[@]}" 2>&1 | tee "${experiment_dir}/build.log"

comparison_directory="${experiment_dir}/hciAnalyze"
EXPERIMENT_DIR="${comparison_directory}" \
SCIENCE_FILE="${science_file}" \
RESPONSE_MANIFEST="${sparse_manifest}" \
EXACT_RESPONSE_MANIFEST="${product_directory}/klipExact_manifest.fits" \
HCIANALYZE_BIN="${hcianalyze_bin}" \
GAUSSIAN_FWHMS="${gaussian_fwhms}" \
LAMBDA_D="${lambda_d}" \
PLANET_SEP="${planet_sep}" \
PLANET_PA="${planet_pa}" \
PLANET_CONTRAST="${planet_contrast}" \
PLANET_RADIUS="${planet_radius}" \
SNR_MIN_RADIUS="${snr_min_radius}" \
SNR_MAX_RADIUS="${snr_max_radius}" \
SNR_APERTURE_RADIUS="${snr_aperture_radius}" \
    "${script_dir}/run_klip_hciAnalyze_filter_comparison.sh"

: > "${experiment_dir}/complete"
printf '\nKLIP exact-response validation complete: %s\n' \
    "${comparison_directory}/hciAnalyze_filter_summary.md"
