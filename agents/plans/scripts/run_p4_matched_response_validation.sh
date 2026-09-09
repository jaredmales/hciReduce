#!/usr/bin/env bash

set -euo pipefail

export LC_ALL=C

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "${script_dir}/../../.." && pwd)
roc_working_dir="${repo_root}/working/roc"
base_config=${BASE_CONFIG:-"${roc_working_dir}/p4Reduce_afLepNaco.conf"}
psf_file=${PSF_FILE:-/home/jrmales/Source/mxWork/NACO/AFLep/2011-10-21/out/psf_reg_median.fits}
p4reduce_bin=${P4REDUCE_BIN:-p4Reduce}
experiment_dir=${EXPERIMENT_DIR:-"${roc_working_dir}/p4_matched_response_$(date -u +%Y%m%dT%H%M%SZ)"}

planet_sep=${PLANET_SEP:-11.782}
planet_pa=${PLANET_PA:-262.051}
planet_contrast=${PLANET_CONTRAST:-0.004878}
mode_fraction=${MODE_FRACTION:-0.15}
position_bound=${POSITION_BOUND:-1}
psf_stamp_size=${PSF_STAMP_SIZE:-11}
psf_sample_avoid_radius=${PSF_SAMPLE_AVOID_RADIUS:-5}
noise_exclusion_radius=${NOISE_EXCLUSION_RADIUS:-5}
noise_min_radius=${NOISE_MIN_RADIUS:-6}
noise_max_radius=${NOISE_MAX_RADIUS:-60}
lambda_d=${LAMBDA_D:-3.6}

optimizer_aperture_radius=${OPTIMIZER_APERTURE_RADIUS:-5}
optimizer_local_stamp_size=${OPTIMIZER_LOCAL_STAMP_SIZE:-13}
optimizer_contrast_lower=${OPTIMIZER_CONTRAST_LOWER:--0.05}
optimizer_contrast_upper=${OPTIMIZER_CONTRAST_UPPER:-0}
optimizer_max_evaluations=${OPTIMIZER_MAX_EVALUATIONS:-192}
optimizer_validation_samples=${OPTIMIZER_VALIDATION_SAMPLES:-21}
optimizer_uncertainty_blocks_explicit=${OPTIMIZER_UNCERTAINTY_BLOCKS+x}
optimizer_uncertainty_blocks=${OPTIMIZER_UNCERTAINTY_BLOCKS:-0}

sparse_response_reference=${SPARSE_RESPONSE_DIR:-}
optimizer_reference=${OPTIMIZER_PRODUCTS_DIR:-}
signal_free_reference=${SIGNAL_FREE_RESPONSE_DIR:-}
dry_run=false

usage()
{
    cat <<EOF
Usage: $(basename "$0") [--dry-run]

Run the AF Lep/NACO validation sequence for a sparse response-backed matched
filter. Completed stages are skipped, so EXPERIMENT_DIR can be resumed.

Stages:
  1. Run a target-aware region2/a4 response on the original science data.
  2. Fit position and contrast from that fixed sparse matched response.
  3. Run the finite-amplitude local P4 negative-companion optimizer.
  4. Subtract the fitted planet once and calculate a dense signal-free response.
  5. Apply that oracle response to the same original science image and compare fits.

Environment overrides:
  EXPERIMENT_DIR                 output directory (default: ${experiment_dir})
  BASE_CONFIG                    standard ROC P4 configuration (default: ${base_config})
  P4REDUCE_BIN                   p4Reduce executable (default: ${p4reduce_bin})
  PSF_FILE                       centered P4-input PSF template (default: ${psf_file})
  PLANET_SEP                     initial separation (default: ${planet_sep})
  PLANET_PA                      initial PA east of north (default: ${planet_pa})
  PLANET_CONTRAST                initial positive contrast (default: ${planet_contrast})
  MODE_FRACTION                  exact P4 mode used by every fit (default: ${mode_fraction})
  POSITION_BOUND                 Cartesian fit half-width in pixels (default: ${position_bound})
  PSF_SAMPLE_AVOID_RADIUS        sparse detector-sample exclusion (default: ${psf_sample_avoid_radius})
  NOISE_EXCLUSION_RADIUS         matched-filter noise exclusion (default: ${noise_exclusion_radius})
  NOISE_MIN_RADIUS               inner noise-profile radius (default: ${noise_min_radius})
  NOISE_MAX_RADIUS               outer noise-profile radius (default: ${noise_max_radius})
  LAMBDA_D                       pixels per lambda/D for small-sample correction (default: ${lambda_d})
  OPTIMIZER_UNCERTAINTY_BLOCKS   exact-fit jackknife blocks; zero disables (default: ${optimizer_uncertainty_blocks})
  SPARSE_RESPONSE_DIR            reuse a completed sparse case directory
  OPTIMIZER_PRODUCTS_DIR         reuse a directory containing a converged p4Negative_* point fit
  SIGNAL_FREE_RESPONSE_DIR       reuse a completed, planet-subtracted dense case directory
  OMP_NUM_THREADS                OpenMP worker limit passed through to p4Reduce

Example:
  OMP_NUM_THREADS=48 nohup $(basename "$0") > p4_matched_response_driver.log 2>&1 &

Reuse the already completed exact optimizer if desired:
  OPTIMIZER_PRODUCTS_DIR=/path/to/finim_outputs \
    OMP_NUM_THREADS=48 nohup $(basename "$0") > p4_matched_response_driver.log 2>&1 &
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

absolute_directory()
{
    local directory=$1
    [[ -d "${directory}" ]] || { printf 'Directory is not readable: %s\n' "${directory}" >&2; exit 1; }
    (cd -- "${directory}" && pwd)
}

manifest_complete()
{
    local case_directory=$1
    local manifest="${case_directory}/finim_outputs/p4PSF_manifest.fits"
    [[ -f "${manifest}" ]] || return 1
    python3 - "${manifest}" <<'PY'
import sys
from astropy.io import fits

try:
    complete = int(fits.getheader(sys.argv[1]).get("P4 PSF COMPLETE", 0)) == 1
except Exception:
    complete = False
raise SystemExit(0 if complete else 1)
PY
}

dense_complete()
{
    local case_directory=$1
    local manifest="${case_directory}/finim_outputs/p4PSF_manifest.fits"
    manifest_complete "${case_directory}" || return 1
    python3 - "${manifest}" <<'PY'
import sys
from astropy.io import fits

try:
    header = fits.getheader(sys.argv[1])
    complete = (
        str(header.get("P4 PSF SPATIAL MODEL", "")).strip() == "PER_PIXEL"
        and int(header.get("P4 PSF RADII PER REGION", 0)) == 0
        and int(header.get("P4 PSF SAMPLES PER RADIUS", 0)) == 0
        and not str(header.get("P4 PSF SAMPLE RADII", "")).strip()
        and int(header.get("P4 PSF MODEL OUTPUT", 0)) == 1
    )
except Exception:
    complete = False
raise SystemExit(0 if complete else 1)
PY
}

sparse_complete()
{
    local case_directory=$1
    local manifest="${case_directory}/finim_outputs/p4PSF_manifest.fits"
    manifest_complete "${case_directory}" || return 1
    python3 - "${manifest}" <<'PY'
import sys
from astropy.io import fits

try:
    header = fits.getheader(sys.argv[1])
    complete = (
        int(header.get("P4 PSF PRODUCT SCHEMA", 0)) == 6
        and str(header.get("P4 PSF SPATIAL MODEL", "")).strip() == "REGION_TARGET_RADIAL_LINEAR"
        and str(header.get("P4 PSF COMPOSITION", "")).strip() == "TARGET_PIXEL"
        and int(header.get("P4 PSF RADII PER REGION", 0)) == 2
        and int(header.get("P4 PSF SAMPLES PER RADIUS", 0)) == 4
        and float(header.get("P4 PSF SAMPLE AVOID RADIUS", 0)) > 0
        and int(header.get("P4 PSF MODEL OUTPUT", 0)) == 1
    )
except Exception:
    complete = False
raise SystemExit(0 if complete else 1)
PY
}

optimizer_summary_path()
{
    local products=$1
    [[ -d "${products}" ]] || return 1
    python3 - "${products}" <<'PY'
import math
import pathlib
import re
import sys

products = pathlib.Path(sys.argv[1])
candidates = sorted(products.glob("p4Negative*_summary.yaml"))
if len(candidates) != 1:
    raise SystemExit(1)
summary = candidates[0]
text = summary.read_text(encoding="utf-8")
status = re.search(r'^  status: "([^"]+)"$', text, flags=re.MULTILINE)
converged = re.search(r"^  converged: (true|false)$", text, flags=re.MULTILINE)
dense = re.search(r"^  denseAgreement: (true|false)$", text, flags=re.MULTILINE)
fitted = re.search(r"^  fitted:\n(?P<body>(?:    [^\n]*\n)+)", text, flags=re.MULTILINE)
if (
    status is None
    or status.group(1) != "converged"
    or converged is None
    or converged.group(1) != "true"
    or dense is None
    or dense.group(1) != "true"
    or fitted is None
):
    raise SystemExit(1)
values = []
for key in ("separation", "positionAngle", "contrast"):
    match = re.search(rf"^    {key}: ([+\-0-9.eE]+)$", fitted.group("body"), flags=re.MULTILINE)
    if match is None:
        raise SystemExit(1)
    values.append(float(match.group(1)))
if not all(math.isfinite(value) for value in values) or values[0] < 0 or values[2] >= 0:
    raise SystemExit(1)
print(summary)
PY
}

optimizer_point_complete()
{
    local products=$1
    optimizer_summary_path "${products}" >/dev/null
}

fit_complete()
{
    local fit_directory=$1
    [[ -s "${fit_directory}/summary.json" && -s "${fit_directory}/surface.csv" ]]
}

run_timed()
{
    local stage_directory=$1
    shift
    local command_line=("$@")
    mkdir -p "${stage_directory}"
    shell_join "${command_line[@]}" > "${stage_directory}/command.txt"
    printf '\n[%s]\n' "$(basename "${stage_directory}")"
    shell_join "${command_line[@]}"
    if [[ "${dry_run}" == true ]]; then
        return 0
    fi
    set +e
    /usr/bin/time -f 'wall_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\nmaximum_rss_kib=%M' \
        -o "${stage_directory}/resource_usage.txt" \
        "${command_line[@]}" 2>&1 | tee "${stage_directory}/run.log"
    local command_status=${PIPESTATUS[0]}
    set -e
    if ((command_status != 0)); then
        printf '[%s] failed with status %d.\n' "$(basename "${stage_directory}")" "${command_status}" >&2
    fi
    return "${command_status}"
}

[[ -r "${base_config}" ]] || { printf 'Base configuration is not readable: %s\n' "${base_config}" >&2; exit 1; }
[[ -r "${psf_file}" ]] || { printf 'PSF template is not readable: %s\n' "${psf_file}" >&2; exit 1; }
[[ -x /usr/bin/time ]] || { printf '%s\n' '/usr/bin/time is required.' >&2; exit 1; }
command -v python3 >/dev/null 2>&1 || { printf '%s\n' 'python3 is required.' >&2; exit 1; }
python3 -c 'import astropy.io.fits, numpy' >/dev/null 2>&1 || {
    printf '%s\n' 'The Python astropy and numpy packages are required.' >&2
    exit 1
}
command -v "${p4reduce_bin}" >/dev/null 2>&1 || {
    printf 'p4Reduce executable was not found: %s\n' "${p4reduce_bin}" >&2
    exit 1
}

help_text=$("${p4reduce_bin}" --help 2>&1)
for required_option in --p4.modeFractions --p4.psfRadiiPerRegion --p4.psfSamplingMode \
    --p4Optimize.enabled --fake.subtractPlanet; do
    if [[ "${help_text}" != *"${required_option}"* ]]; then
        printf 'p4Reduce does not expose required option %s; build this checkout first.\n' "${required_option}" >&2
        exit 1
    fi
done
if grep -Eq '^[[:space:]]*(psfSampleRadii|psfRadiiPerRegion|psfSamplesPerRadius|psfSampleArcStep)[[:space:]]*=' "${base_config}"; then
    printf 'Base configuration contains an active sparse PSF grid: %s\n' "${base_config}" >&2
    printf '%s\n' 'Use the unsampled standard ROC configuration so the dense oracle is unambiguous.' >&2
    exit 1
fi

negative_initial_contrast=$(python3 - "${planet_contrast}" <<'PY'
import math
import sys

value = float(sys.argv[1])
if not math.isfinite(value) or value <= 0:
    raise SystemExit("PLANET_CONTRAST must be finite and positive")
print(f"{-value:.17g}")
PY
)

mkdir -p "${experiment_dir}"
experiment_dir=$(absolute_directory "${experiment_dir}")
if [[ -z "${optimizer_uncertainty_blocks_explicit}" && -f "${experiment_dir}/settings.txt" ]]; then
    saved_uncertainty_blocks=$(awk -F= '$1 == "optimizer_uncertainty_blocks" { print $2 }' \
        "${experiment_dir}/settings.txt")
    if [[ "${saved_uncertainty_blocks}" =~ ^[0-9]+$ ]]; then
        optimizer_uncertainty_blocks=${saved_uncertainty_blocks}
    fi
fi
base_snapshot="${experiment_dir}/p4Reduce_afLepNaco.base.conf"
if [[ -e "${base_snapshot}" ]]; then
    cmp -s "${base_config}" "${base_snapshot}" || {
        printf 'Base configuration differs from the experiment snapshot: %s\n' "${base_snapshot}" >&2
        exit 1
    }
else
    cp -- "${base_config}" "${base_snapshot}"
fi

binary_path=$(command -v "${p4reduce_bin}")
binary_sha256=$(sha256sum "${binary_path}" | awk '{ print $1 }')
base_config_sha256=$(sha256sum "${base_config}" | awk '{ print $1 }')
psf_file_sha256=$(sha256sum "${psf_file}" | awk '{ print $1 }')
settings_content=$(printf '%s\n' \
    "p4reduce_path=${binary_path}" \
    "p4reduce_sha256=${binary_sha256}" \
    "base_config=${base_config}" \
    "base_config_sha256=${base_config_sha256}" \
    "psf_file=${psf_file}" \
    "psf_file_sha256=${psf_file_sha256}" \
    "planet_sep=${planet_sep}" \
    "planet_pa=${planet_pa}" \
    "planet_contrast=${planet_contrast}" \
    "mode_fraction=${mode_fraction}" \
    "position_bound=${position_bound}" \
    "psf_stamp_size=${psf_stamp_size}" \
    "psf_sample_avoid_radius=${psf_sample_avoid_radius}" \
    "noise_exclusion_radius=${noise_exclusion_radius}" \
    "noise_min_radius=${noise_min_radius}" \
    "noise_max_radius=${noise_max_radius}" \
    "lambda_d=${lambda_d}" \
    "optimizer_aperture_radius=${optimizer_aperture_radius}" \
    "optimizer_local_stamp_size=${optimizer_local_stamp_size}" \
    "optimizer_contrast_lower=${optimizer_contrast_lower}" \
    "optimizer_contrast_upper=${optimizer_contrast_upper}" \
    "optimizer_max_evaluations=${optimizer_max_evaluations}" \
    "optimizer_validation_samples=${optimizer_validation_samples}" \
    "optimizer_uncertainty_blocks=${optimizer_uncertainty_blocks}" \
    "sparse_response_reference=${sparse_response_reference}" \
    "optimizer_reference=${optimizer_reference}" \
    "signal_free_reference=${signal_free_reference}" \
    "omp_num_threads=${OMP_NUM_THREADS:-unlimited}")
settings_file="${experiment_dir}/settings.txt"
if [[ -e "${settings_file}" ]]; then
    [[ "$(<"${settings_file}")" == "${settings_content}" ]] || {
        printf 'Experiment settings differ from the saved settings: %s\n' "${settings_file}" >&2
        printf '%s\n' 'Choose a new EXPERIMENT_DIR or restore the original environment overrides.' >&2
        exit 1
    }
else
    printf '%s\n' "${settings_content}" > "${settings_file}"
fi

provenance_file="${experiment_dir}/provenance.txt"
if [[ ! -e "${provenance_file}" ]]; then
    {
        printf 'created_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        printf 'hostname=%s\n' "$(hostname)"
        printf 'hciReduce_commit=%s\n' "$(git -C "${repo_root}" rev-parse HEAD 2>/dev/null || printf unknown)"
        printf '%s\n' "${settings_content}"
    } > "${provenance_file}"
fi

printf 'Experiment directory: %s\n' "${experiment_dir}"
printf 'P4 executable: %s\n' "${binary_path}"
printf 'Mode fraction: %s\n' "${mode_fraction}"
printf 'OMP_NUM_THREADS: %s\n' "${OMP_NUM_THREADS:-unlimited}"

if [[ -n "${sparse_response_reference}" ]]; then
    sparse_case=$(absolute_directory "${sparse_response_reference}")
    sparse_complete "${sparse_case}" || {
        printf 'SPARSE_RESPONSE_DIR is not a completed schema-6 region2/a4 case: %s\n' "${sparse_case}" >&2
        exit 1
    }
else
    sparse_case="${experiment_dir}/sparse_response"
    if sparse_complete "${sparse_case}"; then
        printf '\n[sparse_response] completed response exists; skipping.\n'
    elif [[ -e "${sparse_case}/run.log" || -e "${sparse_case}/finim.fits" ]]; then
        printf 'Incomplete sparse-response stage exists; refusing to overwrite: %s\n' "${sparse_case}" >&2
        exit 1
    else
        sparse_command=(
            "${p4reduce_bin}"
            --config "${base_config}"
            --input.imSize 256
            --p4.modeFractions "${mode_fraction}"
            --planet.sep "${planet_sep}"
            --planet.PA "${planet_pa}"
            --planet.contrast "${planet_contrast}"
            --p4.psfFile "${psf_file}"
            --p4.psfStampSize "${psf_stamp_size}"
            --p4.outputPSFModels=true
            --p4.psfFilter=true
            --p4.psfOutputPrefix p4PSF_
            --p4.psfSamplingMode detectorLocal
            --p4.psfRadiiPerRegion 2
            --p4.psfSamplesPerRadius 4
            --p4.psfSampleAvoidRadius "${psf_sample_avoid_radius}"
            --p4Optimize.enabled=false
            --output.directory "${sparse_case}"
            --output.fileName finim.fits
            --output.exactFName=true
            --showTiming=true
        )
        run_timed "${sparse_case}" "${sparse_command[@]}"
        if [[ "${dry_run}" == false ]] && ! sparse_complete "${sparse_case}"; then
            printf 'Sparse response did not publish the required schema-6 manifest.\n' >&2
            exit 1
        fi
    fi
fi

sparse_fit="${experiment_dir}/sparse_fit"
if [[ "${dry_run}" == false ]]; then
    if fit_complete "${sparse_fit}"; then
        printf '\n[sparse_fit] completed fit exists; skipping.\n'
    else
        [[ ! -e "${sparse_fit}" ]] || {
            printf 'Incomplete sparse-fit directory exists; refusing to overwrite: %s\n' "${sparse_fit}" >&2
            exit 1
        }
        sparse_fit_command=(
            python3 "${script_dir}/fit_p4_matched_response.py"
            "${sparse_case}/finim.fits"
            "${sparse_case}/finim_outputs/p4PSF_manifest.fits"
            "${sparse_fit}"
            --mode-fraction "${mode_fraction}"
            --initial-separation "${planet_sep}"
            --initial-pa "${planet_pa}"
            --position-bound "${position_bound}"
            --noise-exclusion-radius "${noise_exclusion_radius}"
            --noise-min-radius "${noise_min_radius}"
            --noise-max-radius "${noise_max_radius}"
            --lambda-d "${lambda_d}"
        )
        run_timed "${sparse_fit}" "${sparse_fit_command[@]}"
    fi
fi

if [[ -n "${optimizer_reference}" ]]; then
    optimizer_products=$(absolute_directory "${optimizer_reference}")
    if ! optimizer_point_complete "${optimizer_products}" && \
        optimizer_point_complete "${optimizer_products}/finim_outputs"; then
        optimizer_products="${optimizer_products}/finim_outputs"
    fi
    optimizer_point_complete "${optimizer_products}" || {
        printf 'OPTIMIZER_PRODUCTS_DIR does not contain a converged p4Negative point fit: %s\n' \
            "${optimizer_products}" >&2
        exit 1
    }
else
    optimizer_case="${experiment_dir}/exact_optimizer"
    optimizer_products="${optimizer_case}/finim_outputs"
    if optimizer_point_complete "${optimizer_products}"; then
        printf '\n[exact_optimizer] converged point-fit products exist; skipping.\n'
    elif [[ -e "${optimizer_case}/run.log" || -e "${optimizer_products}" ]]; then
        printf 'Incomplete exact-optimizer stage exists; refusing to overwrite: %s\n' "${optimizer_case}" >&2
        exit 1
    else
        optimizer_command=(
            "${p4reduce_bin}"
            --config "${base_config}"
            --input.imSize 256
            --p4.modeFractions "${mode_fraction}"
            --planet.sep "${planet_sep}"
            --planet.PA "${planet_pa}"
            --planet.contrast "${planet_contrast}"
            --fake.method single
            --fake.fileName "${psf_file}"
            --fake.sep "${planet_sep}"
            --fake.PA "${planet_pa}"
            --fake.contrast "${negative_initial_contrast}"
            --fake.subtractPlanet=false
            --p4.localStampSize "${optimizer_local_stamp_size}"
            --p4.psfFile ""
            --p4.outputPSFModels=false
            --p4.psfFilter=false
            --p4Optimize.enabled=true
            --p4Optimize.modeFraction "${mode_fraction}"
            --p4Optimize.apertureRadius "${optimizer_aperture_radius}"
            --p4Optimize.fitPosition=true
            --p4Optimize.contrastLower "${optimizer_contrast_lower}"
            --p4Optimize.contrastUpper "${optimizer_contrast_upper}"
            --p4Optimize.maxEvaluations "${optimizer_max_evaluations}"
            --p4Optimize.validationSamples "${optimizer_validation_samples}"
            --p4Optimize.parameterTolerance 1e-5
            --p4Optimize.positionTolerance 0.001
            --p4Optimize.meritTolerance 1e-6
            --p4Optimize.positionBound "${position_bound}"
            --p4Optimize.uncertaintyBlocks "${optimizer_uncertainty_blocks}"
            --p4Optimize.outputPrefix p4Negative_
            --output.directory "${optimizer_case}"
            --output.fileName finim.fits
            --output.exactFName=true
            --showTiming=true
        )
        optimizer_run_status=0
        run_timed "${optimizer_case}" "${optimizer_command[@]}" || optimizer_run_status=$?
        if [[ "${dry_run}" == false ]]; then
            if ! optimizer_point_complete "${optimizer_products}"; then
                printf 'Exact optimizer did not publish a converged point fit.\n' >&2
                if ((optimizer_run_status != 0)); then
                    exit "${optimizer_run_status}"
                fi
                exit 1
            fi
            if ((optimizer_run_status != 0)); then
                printf '%s\n' \
                    'Exact point fit converged; continuing despite incomplete optional uncertainty products.' >&2
            fi
        fi
    fi
fi

if [[ "${dry_run}" == true ]] && ! optimizer_point_complete "${optimizer_products}"; then
    printf '\nExecution-dependent fit and signal-free commands are omitted in dry-run mode.\n'
    exit 0
fi

optimizer_summary=$(optimizer_summary_path "${optimizer_products}")
readarray -t fitted_planet < <(python3 - "${optimizer_summary}" <<'PY'
import math
import re
import sys

text = open(sys.argv[1], encoding="utf-8").read()
fitted = re.search(r"^  fitted:\n(?P<body>(?:    [^\n]*\n)+)", text, flags=re.MULTILINE)
if fitted is None:
    raise SystemExit("optimizer summary has no fitted point")
values = []
for key in ("separation", "positionAngle", "contrast"):
    match = re.search(rf"^    {key}: ([+\-0-9.eE]+)$", fitted.group("body"), flags=re.MULTILINE)
    if match is None:
        raise SystemExit(f"optimizer summary has no fitted {key}")
    values.append(float(match.group(1)))
values[2] = -values[2]
if not all(math.isfinite(value) for value in values) or values[0] < 0 or values[2] <= 0:
    raise SystemExit("fitted optimizer values are not a valid positive planet")
for value in values:
    print(f"{value:.17g}")
PY
)
fitted_sep=${fitted_planet[0]}
fitted_pa=${fitted_planet[1]}
fitted_contrast=${fitted_planet[2]}
printf '\nExact fitted planet: separation=%s PA=%s contrast=%s\n' \
    "${fitted_sep}" "${fitted_pa}" "${fitted_contrast}"

if [[ -n "${signal_free_reference}" ]]; then
    signal_free_case=$(absolute_directory "${signal_free_reference}")
    dense_complete "${signal_free_case}" || {
        printf 'SIGNAL_FREE_RESPONSE_DIR is not a completed dense P4 response case: %s\n' \
            "${signal_free_case}" >&2
        exit 1
    }
else
    signal_free_case="${experiment_dir}/signal_free_oracle"
    if dense_complete "${signal_free_case}"; then
        printf '\n[signal_free_oracle] completed response exists; skipping.\n'
    elif [[ -e "${signal_free_case}/run.log" || -e "${signal_free_case}/finim.fits" ]]; then
        printf 'Incomplete signal-free stage exists; refusing to overwrite: %s\n' "${signal_free_case}" >&2
        exit 1
    else
        signal_free_command=(
            "${p4reduce_bin}"
            --config "${base_config}"
            --input.imSize 256
            --p4.modeFractions "${mode_fraction}"
            --planet.sep "${fitted_sep}"
            --planet.PA "${fitted_pa}"
            --planet.contrast "${fitted_contrast}"
            --fake.method single
            --fake.fileName "${psf_file}"
            --fake.subtractPlanet=true
            --p4.localStampSize 0
            --p4.psfFile "${psf_file}"
            --p4.psfStampSize "${psf_stamp_size}"
            --p4.outputPSFModels=true
            --p4.psfFilter=true
            --p4.psfOutputPrefix p4PSF_
            --p4Optimize.enabled=false
            --output.directory "${signal_free_case}"
            --output.fileName finim.fits
            --output.exactFName=true
            --showTiming=true
        )
        run_timed "${signal_free_case}" "${signal_free_command[@]}"
        if [[ "${dry_run}" == false ]] && ! dense_complete "${signal_free_case}"; then
            printf 'Signal-free oracle did not publish a completed dense response manifest.\n' >&2
            exit 1
        fi
    fi
fi

oracle_fit="${experiment_dir}/signal_free_fit"
if fit_complete "${oracle_fit}"; then
    printf '\n[signal_free_fit] completed fit exists; skipping.\n'
else
    [[ ! -e "${oracle_fit}/run.log" && ! -e "${oracle_fit}/summary.json" ]] || {
        printf 'Incomplete signal-free-fit directory exists; refusing to overwrite: %s\n' "${oracle_fit}" >&2
        exit 1
    }
    oracle_fit_command=(
        python3 "${script_dir}/fit_p4_matched_response.py"
        "${sparse_case}/finim.fits"
        "${signal_free_case}/finim_outputs/p4PSF_manifest.fits"
        "${oracle_fit}"
        --mode-fraction "${mode_fraction}"
        --initial-separation "${planet_sep}"
        --initial-pa "${planet_pa}"
        --position-bound "${position_bound}"
        --noise-exclusion-radius "${noise_exclusion_radius}"
        --noise-min-radius "${noise_min_radius}"
        --noise-max-radius "${noise_max_radius}"
        --lambda-d "${lambda_d}"
    )
    run_timed "${oracle_fit}" "${oracle_fit_command[@]}"
fi

if [[ "${dry_run}" == true ]]; then
    printf '\nComparison requires completed fit products and is omitted in dry-run mode.\n'
    exit 0
fi

python3 "${script_dir}/compare_p4_matched_response.py" \
    "${sparse_fit}/summary.json" \
    "${oracle_fit}/summary.json" \
    "${optimizer_products}" \
    "${experiment_dir}"

printf '\nValidation complete. Summary: %s\n' "${experiment_dir}/fit_comparison.md"
