#!/usr/bin/env bash

set -euo pipefail

export LC_ALL=C
export PYTHONDONTWRITEBYTECODE=1

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "${script_dir}/../../.." && pwd)
roc_working_dir="${repo_root}/working/roc"
reference_experiment=${REFERENCE_EXPERIMENT:-"${roc_working_dir}/p4_matched_response_20260907T225744Z"}
base_config=${BASE_CONFIG:-"${roc_working_dir}/p4Reduce_afLepNaco.conf"}
psf_file=${PSF_FILE:-/home/jrmales/Source/mxWork/NACO/AFLep/2011-10-21/out/psf_reg_median.fits}
p4reduce_bin=${P4REDUCE_BIN:-p4Reduce}
experiment_dir=${EXPERIMENT_DIR:-"${roc_working_dir}/p4_mean_combine_$(date -u +%Y%m%dT%H%M%SZ)"}
mode_fraction=${MODE_FRACTION:-0.15}
dry_run=false

usage()
{
    cat <<EOF
Usage: $(basename "$0") [--dry-run]

Test whether final-image sigma clipping explains the P4 analytic-response flux
discrepancy. The driver runs only two selected-mode reductions:

  1. Original data combined with an arithmetic mean.
  2. Exact fitted planet subtracted, then combined with an arithmetic mean.

It projects their finite difference onto the existing dense signal-free response
from REFERENCE_EXPERIMENT and compares it with the corresponding sigmaMean pair.
Neither the exact optimizer nor the PSF response field is recalculated.

Environment overrides:
  REFERENCE_EXPERIMENT   completed matched-response experiment
                        (default: ${reference_experiment})
  EXPERIMENT_DIR         output directory (default: ${experiment_dir})
  BASE_CONFIG            standard ROC P4 configuration (default: ${base_config})
  P4REDUCE_BIN           p4Reduce executable (default: ${p4reduce_bin})
  PSF_FILE               centered input PSF used for planet subtraction
                        (default: ${psf_file})
  MODE_FRACTION          sole P4 mode to calculate (default: ${mode_fraction})
  OMP_NUM_THREADS        OpenMP worker limit passed through to p4Reduce

Example:
  OMP_NUM_THREADS=48 nohup $(basename "$0") > p4_mean_combine_driver.log 2>&1 &
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

optimizer_summary_path()
{
    local reference=$1
    python3 - "${reference}" <<'PY'
import math
import pathlib
import re
import sys

products = pathlib.Path(sys.argv[1]) / "exact_optimizer" / "finim_outputs"
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

mean_reduction_complete()
{
    local final_image=$1
    local expects_subtraction=$2
    [[ -s "${final_image}" ]] || return 1
    python3 - "${final_image}" "${mode_fraction}" "${expects_subtraction}" <<'PY'
import sys

import numpy as np
from astropy.io import fits

path, requested_text, expects_subtraction_text = sys.argv[1:]
requested = float(requested_text)
expects_subtraction = expects_subtraction_text == "true"
try:
    header = fits.getheader(path)
    data = np.asarray(fits.getdata(path))
    modes = np.asarray([float(value) for value in str(header["P4 MODE FRACTIONS"]).split(",")])
    matches = np.flatnonzero(
        np.abs(modes - requested)
        <= 8
        * np.finfo(np.float32).eps
        * np.maximum.reduce((np.ones_like(modes), np.abs(modes), np.full_like(modes, abs(requested))))
    )
    complete = (
        str(header.get("COMBINATION METHOD", "")).strip() == "mean"
        and matches.size == 1
        and modes.size == 1
        and data.size > 0
        and np.any(np.isfinite(data))
        and bool(str(header.get("FAKEFILE", "")).strip()) == expects_subtraction
    )
except Exception:
    complete = False
raise SystemExit(0 if complete else 1)
PY
}

reference_response_complete()
{
    local manifest=$1
    [[ -s "${manifest}" ]] || return 1
    python3 - "${manifest}" "${mode_fraction}" <<'PY'
import sys

import numpy as np
from astropy.io import fits

path, requested_text = sys.argv[1:]
requested = float(requested_text)
try:
    header = fits.getheader(path)
    modes = np.asarray([float(value) for value in str(header["P4 MODE FRACTIONS"]).split(",")])
    matches = np.flatnonzero(
        np.abs(modes - requested)
        <= 8
        * np.finfo(np.float32).eps
        * np.maximum.reduce((np.ones_like(modes), np.abs(modes), np.full_like(modes, abs(requested))))
    )
    complete = (
        int(header.get("P4 PSF COMPLETE", 0)) == 1
        and str(header.get("P4 PSF SPATIAL MODEL", "")).strip() == "PER_PIXEL"
        and str(header.get("P4 PSF COMBINATION", "")).strip() == "mean"
        and int(header.get("P4 PSF MODE COUNT", 0)) == 1
        and modes.size == 1
        and matches.size == 1
    )
except Exception:
    complete = False
raise SystemExit(0 if complete else 1)
PY
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
for required_option in --combine.method --fake.subtractPlanet --p4.modeFractions --p4Optimize.enabled; do
    if [[ "${help_text}" != *"${required_option}"* ]]; then
        printf 'p4Reduce does not expose required option %s; build this checkout first.\n' "${required_option}" >&2
        exit 1
    fi
done

reference_experiment=$(absolute_directory "${reference_experiment}")
optimizer_summary=$(optimizer_summary_path "${reference_experiment}") || {
    printf 'Reference experiment lacks a converged exact point fit: %s\n' "${reference_experiment}" >&2
    exit 1
}
response_manifest="${reference_experiment}/signal_free_oracle/finim_outputs/p4PSF_manifest.fits"
reference_response_complete "${response_manifest}" || {
    printf 'Reference response is not a completed one-mode, dense, mean-combined product: %s\n' \
        "${response_manifest}" >&2
    exit 1
}
reference_original="${reference_experiment}/sparse_response/finim.fits"
reference_signal_free="${reference_experiment}/signal_free_oracle/finim.fits"
response_directory="${reference_experiment}/signal_free_oracle/finim_outputs"
response_coordinates="${response_directory}/p4PSF_coordinates.fits"
response_model="${response_directory}/p4PSF_model_0000.fits"
response_validity="${response_directory}/p4PSF_validity_0000.fits"
for required_product in "${reference_original}" "${reference_signal_free}" "${response_coordinates}" \
    "${response_model}" "${response_validity}"; do
    [[ -r "${required_product}" ]] || {
        printf 'Reference experiment product is missing: %s\n' "${required_product}" >&2
        exit 1
    }
done

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

mkdir -p "${experiment_dir}"
experiment_dir=$(absolute_directory "${experiment_dir}")
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
settings_content=$(printf '%s\n' \
    "p4reduce_path=${binary_path}" \
    "p4reduce_sha256=$(sha256sum "${binary_path}" | awk '{ print $1 }')" \
    "base_config=${base_config}" \
    "base_config_sha256=$(sha256sum "${base_config}" | awk '{ print $1 }')" \
    "psf_file=${psf_file}" \
    "psf_file_sha256=$(sha256sum "${psf_file}" | awk '{ print $1 }')" \
    "reference_experiment=${reference_experiment}" \
    "reference_original_sha256=$(sha256sum "${reference_original}" | awk '{ print $1 }')" \
    "reference_signal_free_sha256=$(sha256sum "${reference_signal_free}" | awk '{ print $1 }')" \
    "response_manifest_sha256=$(sha256sum "${response_manifest}" | awk '{ print $1 }')" \
    "response_coordinates_sha256=$(sha256sum "${response_coordinates}" | awk '{ print $1 }')" \
    "response_model_sha256=$(sha256sum "${response_model}" | awk '{ print $1 }')" \
    "response_validity_sha256=$(sha256sum "${response_validity}" | awk '{ print $1 }')" \
    "optimizer_summary_sha256=$(sha256sum "${optimizer_summary}" | awk '{ print $1 }')" \
    "fitted_sep=${fitted_sep}" \
    "fitted_pa=${fitted_pa}" \
    "fitted_contrast=${fitted_contrast}" \
    "mode_fraction=${mode_fraction}" \
    "combination=mean" \
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
printf 'Reference experiment: %s\n' "${reference_experiment}"
printf 'Exact fitted planet: separation=%s PA=%s contrast=%s\n' \
    "${fitted_sep}" "${fitted_pa}" "${fitted_contrast}"
printf 'Mode fraction: %s\n' "${mode_fraction}"
printf 'OMP_NUM_THREADS: %s\n' "${OMP_NUM_THREADS:-unlimited}"

original_case="${experiment_dir}/mean_original"
original_image="${original_case}/finim.fits"
if mean_reduction_complete "${original_image}" false; then
    printf '\n[mean_original] completed reduction exists; skipping.\n'
elif [[ -e "${original_case}/run.log" || -e "${original_image}" ]]; then
    printf 'Incomplete mean-original stage exists; refusing to overwrite: %s\n' "${original_case}" >&2
    exit 1
else
    original_command=(
        "${p4reduce_bin}"
        --config "${base_config}"
        --input.imSize 256
        --p4.modeFractions "${mode_fraction}"
        --planet.sep "${fitted_sep}"
        --planet.PA "${fitted_pa}"
        --planet.contrast "${fitted_contrast}"
        --fake.fileName ""
        --fake.subtractPlanet=false
        --p4.localStampSize 0
        --p4.psfFile ""
        --p4.outputPSFModels=false
        --p4.psfFilter=false
        --p4Optimize.enabled=false
        --combine.method mean
        --output.directory "${original_case}"
        --output.fileName finim.fits
        --output.exactFName=true
        --showTiming=true
    )
    run_timed "${original_case}" "${original_command[@]}"
    if [[ "${dry_run}" == false ]] && ! mean_reduction_complete "${original_image}" false; then
        printf 'Mean-original reduction did not publish the expected final image.\n' >&2
        exit 1
    fi
fi

signal_free_case="${experiment_dir}/mean_signal_free"
signal_free_image="${signal_free_case}/finim.fits"
if mean_reduction_complete "${signal_free_image}" true; then
    printf '\n[mean_signal_free] completed reduction exists; skipping.\n'
elif [[ -e "${signal_free_case}/run.log" || -e "${signal_free_image}" ]]; then
    printf 'Incomplete mean-signal-free stage exists; refusing to overwrite: %s\n' "${signal_free_case}" >&2
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
        --p4.psfFile ""
        --p4.outputPSFModels=false
        --p4.psfFilter=false
        --p4Optimize.enabled=false
        --combine.method mean
        --output.directory "${signal_free_case}"
        --output.fileName finim.fits
        --output.exactFName=true
        --showTiming=true
    )
    run_timed "${signal_free_case}" "${signal_free_command[@]}"
    if [[ "${dry_run}" == false ]] && ! mean_reduction_complete "${signal_free_image}" true; then
        printf 'Mean-signal-free reduction did not publish the expected final image.\n' >&2
        exit 1
    fi
fi

if [[ "${dry_run}" == true ]]; then
    printf '\nComparison requires completed mean-combined final images and is omitted in dry-run mode.\n'
    exit 0
fi

python3 "${script_dir}/compare_p4_mean_combine.py" \
    "${reference_experiment}" \
    "${original_image}" \
    "${signal_free_image}" \
    "${experiment_dir}" \
    --mode-fraction "${mode_fraction}"

printf '\nMean-combination diagnostic complete. Summary: %s\n' \
    "${experiment_dir}/mean_combine_comparison.md"
