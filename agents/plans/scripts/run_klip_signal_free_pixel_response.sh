#!/usr/bin/env bash

set -euo pipefail

export LC_ALL=C
export PYTHONDONTWRITEBYTECODE=1

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "${script_dir}/../../.." && pwd)
roc_working_dir="${repo_root}/working/roc"
source_experiment=${SOURCE_EXPERIMENT:-"${roc_working_dir}/klip_cpp_response_20260913T222913Z"}
science_file=${SCIENCE_FILE:-"${source_experiment}/science_only/finim.fits"}
sparse_manifest=${SPARSE_RESPONSE_MANIFEST:-"${source_experiment}/radial_ld_refit4_filter/finim_outputs/klipPSF_manifest.fits"}
experiment_dir=${EXPERIMENT_DIR:-"${roc_working_dir}/klip_signal_free_pixel_response_$(date -u +%Y%m%dT%H%M%SZ)"}
base_config=${BASE_CONFIG:-"${script_dir}/klipReduce_afLepNaco_psf_response.conf"}
psf_file=${PSF_FILE:-/home/jrmales/Source/mxWork/NACO/AFLep/2011-10-21/out/psf_reg_median.fits}
klipreduce_bin=${KLIPREDUCE_BIN:-klipReduce}
hcianalyze_bin=${HCIANALYZE_BIN:-hciAnalyze}
mode_counts=${MODE_COUNTS:-125,150,175,200,225,250,300,350}
optimizer_mode=${OPTIMIZER_MODE:-200}
planet_sep=${PLANET_SEP:-11.73795339222688}
planet_pa=${PLANET_PA:-262.1667995998323}
planet_contrast=${PLANET_CONTRAST:-0.004763925929356391}
position_bound=${POSITION_BOUND:-1}
contrast_lower=${CONTRAST_LOWER:-0}
contrast_upper=${CONTRAST_UPPER:-0.01}
optimizer_aperture_radius=${OPTIMIZER_APERTURE_RADIUS:-5}
optimizer_max_evaluations=${OPTIMIZER_MAX_EVALUATIONS:-192}
optimizer_parameter_tolerance=${OPTIMIZER_PARAMETER_TOLERANCE:-0.0005}
optimizer_merit_tolerance=${OPTIMIZER_MERIT_TOLERANCE:-0.00001}
response_fraction=${RESPONSE_FRACTION:-1}
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

Run the full KLIP signal-free response oracle:

  1. Fit separation, PA, and positive contrast by repeatedly injecting a
     negative fake and minimizing the fixed-aperture final-image L2 residual.
  2. Put that positive fit in [planet] and apply it once with
     fake.subtractPlanet=true, producing the signal-free baseline.
  3. At every eligible integer search pixel, measure

       [KLIP(D-Pbest+epsilon Pxy) - KLIP(D-Pbest-epsilon Pxy)] / (2 epsilon).

The native response run writes a schema-2 PIXEL_EXACT manifest. The final stage
uses hciAnalyze to apply that response to the original planet-bearing science
cube and compares it with the sparse response, no filter, and Gaussian smooth.

Environment overrides:
  EXPERIMENT_DIR                 output directory (default: ${experiment_dir})
  BASE_CONFIG                    maintained AF Lep/NACO KLIP config
  PSF_FILE                       centered injection PSF
  KLIPREDUCE_BIN                 klipReduce executable
  HCIANALYZE_BIN                 hciAnalyze executable
  SOURCE_EXPERIMENT              completed native sparse-response experiment
  SCIENCE_FILE                   original planet-bearing final-image cube
  SPARSE_RESPONSE_MANIFEST       sparse response used in the comparison
  MODE_COUNTS                    dense-response KL modes (default: ${mode_counts})
  OPTIMIZER_MODE                 one KL mode used by the negative fit (default: ${optimizer_mode})
  PLANET_SEP, PLANET_PA          initial negative-fit position
  PLANET_CONTRAST                initial positive source contrast
  POSITION_BOUND                 Cartesian fit half-width (default: ${position_bound})
  CONTRAST_LOWER, CONTRAST_UPPER positive contrast bounds
  OPTIMIZER_APERTURE_RADIUS      fixed merit-aperture radius
  OPTIMIZER_MAX_EVALUATIONS      Powell evaluation limit
  OPTIMIZER_PARAMETER_TOLERANCE  normalized Powell coordinate tolerance
  OPTIMIZER_MERIT_TOLERANCE      relative Powell merit tolerance
  RESPONSE_FRACTION              epsilon/fitted contrast (default: ${response_fraction})
  PSF_STAMP_SIZE                 odd per-pixel response width (default: ${psf_stamp_size})
  PSF_FILTER_MIN_GOOD_FRACT      later hciAnalyze filter-support threshold
  GAUSSIAN_FWHMS                 comma-separated comparison widths
  LAMBDA_D, PLANET_RADIUS        hciAnalyze S/N controls
  SNR_MIN_RADIUS, SNR_MAX_RADIUS hciAnalyze noise-annulus bounds
  SNR_APERTURE_RADIUS            hciAnalyze signal aperture
  OMP_NUM_THREADS                OpenMP worker limit passed through to klipReduce

The 6--60 pixel annulus contains about 11,200 integer pixels. At the measured
2.934 seconds per paired-trial reduction, its roughly 22,400 trials should take
about 18 hours. An interrupted native pixel-response stage cannot yet resume
within that stage; completed optimizer evaluations are restartable.

Example:
  OMP_NUM_THREADS=48 nohup $(basename "$0") > klip_signal_free_pixel_response_driver.log 2>&1 &
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
    mkdir -p "${stage_directory}"
    shell_join "${command_line[@]}" > "${stage_directory}/command.txt"
    printf '\n[%s]\n' "${stage_directory#"${experiment_dir}/"}"
    shell_join "${command_line[@]}"
    if [[ "${dry_run}" == true ]]; then
        return
    fi
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
}

pixel_response_complete()
{
    local case_directory=$1
    local manifest="${case_directory}/finim_outputs/klipPSF_manifest.fits"
    [[ -s "${manifest}" ]] || return 1
    python3 - "${manifest}" "${case_directory}" "${mode_counts}" "${fitted_sep}" "${fitted_pa}" \
        "${fitted_contrast}" "${response_contrast}" "${psf_file}" <<'PY'
import math
import pathlib
import sys

import numpy as np
from astropy.io import fits

manifest_path, case_text, modes_text, sep_text, pa_text, contrast_text, response_text, psf_text = sys.argv[1:]
case = pathlib.Path(case_text)

def vector(header, keyword):
    text = str(header.get(keyword, "")).strip()
    return [] if not text else [float(token) for token in text.split(",")]

try:
    header = fits.getheader(manifest_path)
    count = int(str(header["KLIP PSF MEASUREMENT COUNT"]).strip())
    trials = int(str(header["KLIP PSF REFIT TRIAL COUNT"]).strip())
    modes = [int(token) for token in str(header["NMODES"]).split(",")]
    coordinates = np.asarray(fits.getdata(case / "finim_outputs/klipPSF_coordinates.fits"))
    if coordinates.ndim == 2 and coordinates.shape[0] == 4:
        coordinates = coordinates.T
    required = [
        case / "finim.fits",
        *[case / f"finim_outputs/klipPSF_mode{index:03d}_pixel_response.fits" for index in range(len(modes))],
        *[case / f"finim_outputs/klipPSF_mode{index:03d}_pixel_validity.fits" for index in range(len(modes))],
    ]
    sep = vector(header, "PLANETSEP")
    pa = vector(header, "PLANETPA")
    contrast = vector(header, "PLANETCONT")
    complete = (
        int(header.get("KLIP PSF COMPLETE", 0)) == 1
        and int(header.get("KLIP PSF PRODUCT SCHEMA", 0)) == 2
        and str(header.get("KLIP PSF PRODUCT", "")).strip() == "MANIFEST"
        and str(header.get("KLIP PSF SPATIAL MODEL", "")).strip() == "PIXEL_EXACT"
        and str(header.get("KLIP PSF RESPONSE METHOD", "")).strip() == "refitDifference"
        and str(header.get("KLIP PSF ACCUMULATION", "")).strip() == "PAIRED_FINAL_DIFFERENCE"
        and int(header.get("KLIP PSF SAMPLE EVERY PIXEL", 0)) == 1
        and math.isclose(float(header.get("KLIP PSF SAMPLE AVOID RADIUS", -1)), 0, rel_tol=0, abs_tol=1e-12)
        and math.isclose(float(header.get("KLIP PSF REFIT CONTRAST", -1)), float(response_text), rel_tol=1e-6, abs_tol=1e-9)
        and count > 0
        and trials == 2 * count
        and coordinates.shape == (count, 4)
        and modes == [int(token) for token in modes_text.split(",")]
        and len(sep) == len(pa) == len(contrast) == 1
        and math.isclose(sep[0], float(sep_text), rel_tol=0, abs_tol=5e-4)
        and math.isclose(pa[0] % 360, float(pa_text) % 360, rel_tol=0, abs_tol=5e-4)
        and math.isclose(contrast[0], float(contrast_text), rel_tol=0, abs_tol=5e-9)
        and pathlib.Path(str(header["FAKEFILE"]).strip()).resolve() == pathlib.Path(psf_text).resolve()
        and all(path.is_file() for path in required)
    )
except Exception:
    complete = False
raise SystemExit(0 if complete else 1)
PY
}

[[ -x /usr/bin/time ]] || { printf '%s\n' '/usr/bin/time is required.' >&2; exit 1; }
command -v python3 >/dev/null 2>&1 || { printf '%s\n' 'python3 is required.' >&2; exit 1; }
python3 -c 'import astropy.io.fits, numpy, scipy.optimize' >/dev/null 2>&1 || {
    printf '%s\n' 'The Python astropy, numpy, and scipy packages are required.' >&2
    exit 1
}
command -v "${klipreduce_bin}" >/dev/null 2>&1 || {
    printf 'klipReduce executable was not found: %s\n' "${klipreduce_bin}" >&2
    exit 1
}
command -v "${hcianalyze_bin}" >/dev/null 2>&1 || {
    printf 'hciAnalyze executable was not found: %s\n' "${hcianalyze_bin}" >&2
    exit 1
}
base_config=$(absolute_file "${base_config}")
psf_file=$(absolute_file "${psf_file}")
science_file=$(absolute_file "${science_file}")
sparse_manifest=$(absolute_file "${sparse_manifest}")
help_text=$("${klipreduce_bin}" --help 2>&1)
for required_option in --psfResponse.sampleEvery --psfResponse.refitContra --fake.subtractPlanet; do
    if [[ "${help_text}" != *"${required_option}"* ]]; then
        printf 'klipReduce does not expose required option %s; build this checkout first.\n' \
            "${required_option}" >&2
        exit 1
    fi
done

if [[ "${dry_run}" == false ]]; then
    mkdir -p "${experiment_dir}"
fi
optimizer_directory="${experiment_dir}/negative_optimizer"
optimizer_command=(
    python3 "${script_dir}/optimize_klip_negative_planet.py"
    "${base_config}"
    "${psf_file}"
    "${optimizer_directory}"
    --klipreduce-bin "${klipreduce_bin}"
    --mode-count "${optimizer_mode}"
    --initial-separation "${planet_sep}"
    --initial-pa "${planet_pa}"
    --initial-contrast "${planet_contrast}"
    --position-bound "${position_bound}"
    --contrast-lower "${contrast_lower}"
    --contrast-upper "${contrast_upper}"
    --aperture-radius "${optimizer_aperture_radius}"
    --maximum-evaluations "${optimizer_max_evaluations}"
    --parameter-tolerance "${optimizer_parameter_tolerance}"
    --merit-tolerance "${optimizer_merit_tolerance}"
)
printf '\n[negative optimizer]\n'
shell_join "${optimizer_command[@]}"
if [[ "${dry_run}" == true ]]; then
    printf '\n[pixel response]\n'
    printf 'Uses fitted [planet] values, fake.subtractPlanet=true, and refitDifference at every pixel.\n'
    printf '\nDry run complete; fitted values are needed to print the exact response command.\n'
    exit 0
fi
"${optimizer_command[@]}" 2>&1 | tee "${optimizer_directory}/driver.log"

readarray -t fitted_planet < <(python3 - "${optimizer_directory}/summary.json" <<'PY'
import json
import math
import sys

summary = json.load(open(sys.argv[1], encoding="utf-8"))
best = summary.get("best", {})
values = [best.get("separation"), best.get("position_angle"), best.get("contrast")]
if summary.get("complete") is not True or summary.get("status") != "converged":
    raise SystemExit("negative-planet fit is not complete and converged")
if not all(isinstance(value, (int, float)) and math.isfinite(value) for value in values):
    raise SystemExit("negative-planet fit has invalid best parameters")
if values[0] < 0 or values[2] <= 0:
    raise SystemExit("negative-planet fit has invalid separation or contrast")
for value in values:
    print(f"{value:.17g}")
PY
)
if ((${#fitted_planet[@]} != 3)); then
    printf 'Could not read the converged negative-planet fit.\n' >&2
    exit 1
fi
fitted_sep=${fitted_planet[0]}
fitted_pa=${fitted_planet[1]}
fitted_contrast=${fitted_planet[2]}
response_contrast=$(python3 - "${fitted_contrast}" "${response_fraction}" <<'PY'
import math
import sys

contrast, fraction = map(float, sys.argv[1:])
value = contrast * fraction
if not math.isfinite(value) or value <= 0:
    raise SystemExit("RESPONSE_FRACTION must produce a positive finite perturbation")
print(f"{value:.17g}")
PY
)
printf 'Best positive planet: separation=%s PA=%s contrast=%s\n' \
    "${fitted_sep}" "${fitted_pa}" "${fitted_contrast}"
printf 'Paired finite-difference half-amplitude: %s\n' "${response_contrast}"

settings_path="${experiment_dir}/settings.env"
settings_lines=(
    "base_config=${base_config}"
    "psf_file=${psf_file}"
    "science_file=${science_file}"
    "sparse_manifest=${sparse_manifest}"
    "mode_counts=${mode_counts}"
    "optimizer_mode=${optimizer_mode}"
    "fitted_sep=${fitted_sep}"
    "fitted_pa=${fitted_pa}"
    "fitted_contrast=${fitted_contrast}"
    "response_fraction=${response_fraction}"
    "response_contrast=${response_contrast}"
    "psf_stamp_size=${psf_stamp_size}"
)
if [[ -e "${settings_path}" ]]; then
    mapfile -t saved_settings < "${settings_path}"
    [[ "${saved_settings[*]}" == "${settings_lines[*]}" ]] || {
        printf 'Existing experiment settings differ: %s\n' "${settings_path}" >&2
        exit 1
    }
else
    printf '%s\n' "${settings_lines[@]}" > "${settings_path}"
fi

provenance_path="${experiment_dir}/provenance.txt"
if [[ ! -e "${provenance_path}" ]]; then
    {
        printf 'created_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        printf 'hostname=%s\n' "$(hostname)"
        printf 'hciReduce_commit=%s\n' "$(git -C "${repo_root}" rev-parse HEAD)"
        printf 'klipreduce_path=%s\n' "$(command -v "${klipreduce_bin}")"
        printf 'hcianalyze_path=%s\n' "$(command -v "${hcianalyze_bin}")"
        sha256sum "$(command -v "${klipreduce_bin}")" "$(command -v "${hcianalyze_bin}")" \
            "${base_config}" "${psf_file}" "${science_file}" "${sparse_manifest}"
    } > "${provenance_path}"
fi

response_case="${experiment_dir}/signal_free_pixel_response"
response_manifest="${response_case}/finim_outputs/klipPSF_manifest.fits"
if pixel_response_complete "${response_case}"; then
    printf '\n[signal_free_pixel_response] completed response exists; skipping.\n'
elif [[ -e "${response_case}/run.log" || -e "${response_manifest}" ]]; then
    printf 'Incomplete pixel-response stage exists and cannot be resumed: %s\n' "${response_case}" >&2
    exit 1
else
    response_command=(
        "${klipreduce_bin}"
        --config "${base_config}"
        --klip.Nmodes "${mode_counts}"
        --planet.sep "${fitted_sep}"
        --planet.PA "${fitted_pa}"
        --planet.contrast "${fitted_contrast}"
        --fake.method single
        --fake.fileName "${psf_file}"
        --fake.subtractPlanet=true
        --psfResponse.file "${psf_file}"
        --psfResponse.stampSize "${psf_stamp_size}"
        --psfResponse.sampleEveryPixel=true
        --psfResponse.method refitDifference
        --psfResponse.sampleAvoidRadius 0
        --psfResponse.refitContrast "${response_contrast}"
        --psfResponse.outputModels=true
        --psfResponse.filter=false
        --psfResponse.filterMinGoodFract "${psf_filter_min_good_fract}"
        --psfResponse.outputPrefix klipPSF_
        --output.directory "${response_case}"
        --output.fileName finim.fits
        --output.exactFName=true
        --showTiming=true
    )
    run_timed "${response_case}" "${response_command[@]}"
    pixel_response_complete "${response_case}" || {
        printf 'Native KLIP pixel response did not publish a complete schema-2 product.\n' >&2
        exit 1
    }
    : > "${response_case}/complete"
fi

comparison_directory="${experiment_dir}/hciAnalyze"
if [[ -s "${comparison_directory}/hciAnalyze_filter_summary.md" ]]; then
    printf '\n[hciAnalyze] completed comparison exists; skipping.\n'
else
    [[ ! -e "${comparison_directory}" ]] || {
        printf 'Incomplete hciAnalyze comparison exists; refusing to overwrite: %s\n' \
            "${comparison_directory}" >&2
        exit 1
    }
    EXPERIMENT_DIR="${comparison_directory}" \
    SCIENCE_FILE="${science_file}" \
    RESPONSE_MANIFEST="${sparse_manifest}" \
    EXACT_RESPONSE_MANIFEST="${response_manifest}" \
    HCIANALYZE_BIN="${hcianalyze_bin}" \
    GAUSSIAN_FWHMS="${gaussian_fwhms}" \
    LAMBDA_D="${lambda_d}" \
    PLANET_SEP="${fitted_sep}" \
    PLANET_PA="${fitted_pa}" \
    PLANET_CONTRAST="${fitted_contrast}" \
    PLANET_RADIUS="${planet_radius}" \
    SNR_MIN_RADIUS="${snr_min_radius}" \
    SNR_MAX_RADIUS="${snr_max_radius}" \
    SNR_APERTURE_RADIUS="${snr_aperture_radius}" \
        "${script_dir}/run_klip_hciAnalyze_filter_comparison.sh"
fi

: > "${experiment_dir}/complete"
printf '\nKLIP signal-free pixel-response oracle complete: %s\n' \
    "${comparison_directory}/hciAnalyze_filter_summary.md"
