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
experiment_dir=${EXPERIMENT_DIR:-"${roc_working_dir}/p4_refit_difference_$(date -u +%Y%m%dT%H%M%SZ)"}
mode_fraction=${MODE_FRACTION:-0.15}
position_bound=${POSITION_BOUND:-1}
psf_stamp_size=${PSF_STAMP_SIZE:-11}
psf_sample_avoid_radius=${PSF_SAMPLE_AVOID_RADIUS:-5}
noise_exclusion_radius=${NOISE_EXCLUSION_RADIUS:-5}
noise_min_radius=${NOISE_MIN_RADIUS:-6}
noise_max_radius=${NOISE_MAX_RADIUS:-60}
lambda_d=${LAMBDA_D:-3.6}
refit_contrast=${REFIT_CONTRAST:-}
dry_run=false

usage()
{
    cat <<EOF
Usage: $(basename "$0") [--dry-run]

Measure a sparse P4 response with paired +epsilon/-epsilon coefficient refits,
fit the original final image with that response, and compare the result with
the prior frozen-response fit and exact negative-planet optimizer.

Completed stages are skipped, so EXPERIMENT_DIR can be resumed.

Environment overrides:
  REFERENCE_EXPERIMENT       completed matched-response experiment
                            (default: ${reference_experiment})
  EXPERIMENT_DIR             output directory (default: ${experiment_dir})
  BASE_CONFIG                standard ROC P4 configuration (default: ${base_config})
  P4REDUCE_BIN               p4Reduce executable (default: ${p4reduce_bin})
  PSF_FILE                   centered P4-input PSF (default: ${psf_file})
  REFIT_CONTRAST             central-difference half-amplitude; defaults to
                            the exact optimizer's positive planet contrast
  MODE_FRACTION              sole P4 mode to calculate (default: ${mode_fraction})
  PSF_SAMPLE_AVOID_RADIUS    detector radius excluded around known planets
                            (default: ${psf_sample_avoid_radius})
  OMP_NUM_THREADS            OpenMP worker limit passed through to p4Reduce

Example:
  OMP_NUM_THREADS=48 nohup $(basename "$0") > p4_refit_difference_driver.log 2>&1 &
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
    return "${command_status}"
}

response_complete()
{
    local manifest=$1
    local expected_contrast=$2
    [[ -s "${manifest}" ]] || return 1
    python3 - "${manifest}" "${mode_fraction}" "${expected_contrast}" <<'PY'
import math
import sys

import numpy as np
from astropy.io import fits

path, requested_text, contrast_text = sys.argv[1:]
requested = float(requested_text)
contrast = float(contrast_text)
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
        and int(header.get("P4 PSF PRODUCT SCHEMA", 0)) == 7
        and str(header.get("P4 PSF RESPONSE", "")).strip() == "REFIT_CENTRAL_DIFFERENCE"
        and str(header.get("P4 PSF COEFFICIENT SCOPE", "")).strip() == "PAIRED_REFIT"
        and str(header.get("P4 PSF SAMPLING MODE", "")).strip() == "refitDifference"
        and int(header.get("P4 PSF REFIT FIT COUNT", 0)) > 0
        and math.isclose(float(header.get("P4 PSF REFIT CONTRAST", math.nan)), contrast, rel_tol=2e-6)
        and modes.size == 1
        and matches.size == 1
    )
except Exception:
    complete = False
raise SystemExit(0 if complete else 1)
PY
}

fit_complete()
{
    local summary=$1
    local manifest=$2
    [[ -s "${summary}" ]] || return 1
    python3 - "${summary}" "${manifest}" "${mode_fraction}" <<'PY'
import json
import pathlib
import sys

summary_path, manifest_path, requested_text = sys.argv[1:]
try:
    summary = json.loads(pathlib.Path(summary_path).read_text(encoding="utf-8"))
    complete = (
        summary.get("fit", {}).get("status") == "converged"
        and pathlib.Path(summary["manifest"]).resolve() == pathlib.Path(manifest_path).resolve()
        and abs(float(summary["mode_fraction"]) - float(requested_text)) < 1e-6
    )
except Exception:
    complete = False
raise SystemExit(0 if complete else 1)
PY
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
for required_option in --p4.psfFile --p4.psfSamplingMode --p4.psfRefitContrast; do
    if [[ "${help_text}" != *"${required_option}"* ]]; then
        printf 'p4Reduce does not expose required option %s; build this checkout first.\n' "${required_option}" >&2
        exit 1
    fi
done

reference_experiment=$(absolute_directory "${reference_experiment}")
reference_science="${reference_experiment}/sparse_response/finim.fits"
reference_fit="${reference_experiment}/sparse_fit/summary.json"
optimizer_products="${reference_experiment}/exact_optimizer/finim_outputs"
optimizer_summary="${optimizer_products}/p4Negative_summary.yaml"
for required_product in "${reference_science}" "${reference_fit}" "${optimizer_summary}"; do
    [[ -r "${required_product}" ]] || {
        printf 'Reference experiment product is missing: %s\n' "${required_product}" >&2
        exit 1
    }
done

readarray -t exact_planet < <(python3 - "${optimizer_summary}" <<'PY'
import math
import re
import sys

text = open(sys.argv[1], encoding="utf-8").read()
fitted = re.search(r"^  fitted:\n(?P<body>(?:    [^\n]*\n)+)", text, flags=re.MULTILINE)
if fitted is None or not re.search(r"^  converged: true$", text, flags=re.MULTILINE):
    raise SystemExit("exact optimizer did not converge")
values = []
for key in ("separation", "positionAngle", "contrast"):
    match = re.search(rf"^    {key}: ([+\-0-9.eE]+)$", fitted.group("body"), flags=re.MULTILINE)
    if match is None:
        raise SystemExit(f"exact optimizer summary lacks {key}")
    values.append(float(match.group(1)))
values[2] = -values[2]
if values[0] < 0 or values[2] <= 0 or not all(map(math.isfinite, values)):
    raise SystemExit("exact optimizer point is invalid")
for value in values:
    print(f"{value:.17g}")
PY
)
exact_sep=${exact_planet[0]}
exact_pa=${exact_planet[1]}
exact_contrast=${exact_planet[2]}
if [[ -z "${refit_contrast}" ]]; then
    refit_contrast=${exact_contrast}
fi
python3 - "${refit_contrast}" <<'PY'
import math
import sys

value = float(sys.argv[1])
if not math.isfinite(value) or value <= 0:
    raise SystemExit("REFIT_CONTRAST must be finite and positive")
PY

readarray -t sampling_planet < <(python3 - "${reference_science}" <<'PY'
import math
import sys

from astropy.io import fits

header = fits.getheader(sys.argv[1])
values = []
for key in ("PLANETSEP", "PLANETPA", "PLANETCONT"):
    entries = [float(value) for value in str(header[key]).split(",")]
    if len(entries) != 1:
        raise SystemExit(f"reference response does not contain one {key} value")
    values.append(entries[0])
if values[0] < 0 or values[2] < 0 or not all(map(math.isfinite, values)):
    raise SystemExit("reference planet metadata is invalid")
for value in values:
    print(f"{value:.17g}")
PY
)
sampling_sep=${sampling_planet[0]}
sampling_pa=${sampling_planet[1]}
sampling_contrast=${sampling_planet[2]}

mkdir -p "${experiment_dir}"
experiment_dir=$(absolute_directory "${experiment_dir}")
base_snapshot="${experiment_dir}/p4Reduce_afLepNaco.base.conf"
if [[ -e "${base_snapshot}" ]]; then
    cmp -s "${base_config}" "${base_snapshot}" || {
        printf 'Base configuration differs from experiment snapshot: %s\n' "${base_snapshot}" >&2
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
    "exact_sep=${exact_sep}" \
    "exact_pa=${exact_pa}" \
    "exact_contrast=${exact_contrast}" \
    "refit_contrast=${refit_contrast}" \
    "sampling_sep=${sampling_sep}" \
    "sampling_pa=${sampling_pa}" \
    "sampling_contrast=${sampling_contrast}" \
    "mode_fraction=${mode_fraction}" \
    "psf_stamp_size=${psf_stamp_size}" \
    "psf_sample_avoid_radius=${psf_sample_avoid_radius}" \
    "omp_num_threads=${OMP_NUM_THREADS:-unlimited}")
settings_file="${experiment_dir}/settings.txt"
if [[ -e "${settings_file}" ]]; then
    [[ "$(<"${settings_file}")" == "${settings_content}" ]] || {
        printf 'Experiment settings differ from saved settings: %s\n' "${settings_file}" >&2
        exit 1
    }
else
    printf '%s\n' "${settings_content}" > "${settings_file}"
fi
if [[ ! -e "${experiment_dir}/provenance.txt" ]]; then
    {
        printf 'created_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        printf 'hostname=%s\n' "$(hostname)"
        printf 'hciReduce_commit=%s\n' "$(git -C "${repo_root}" rev-parse HEAD 2>/dev/null || printf unknown)"
        printf '%s\n' "${settings_content}"
    } > "${experiment_dir}/provenance.txt"
fi

printf 'Experiment directory: %s\n' "${experiment_dir}"
printf 'Exact planet: sep=%s PA=%s contrast=%s\n' "${exact_sep}" "${exact_pa}" "${exact_contrast}"
printf 'Refit-difference half-amplitude: %s\n' "${refit_contrast}"

response_case="${experiment_dir}/refit_response"
response_manifest="${response_case}/finim_outputs/p4PSF_manifest.fits"
if response_complete "${response_manifest}" "${refit_contrast}"; then
    printf '\n[refit_response] completed response exists; skipping.\n'
elif [[ -e "${response_case}/run.log" || -e "${response_manifest}" ]]; then
    printf 'Incomplete refit response exists; refusing to overwrite: %s\n' "${response_case}" >&2
    exit 1
else
    response_command=(
        "${p4reduce_bin}"
        --config "${base_config}"
        --input.imSize 256
        --p4.modeFractions "${mode_fraction}"
        --planet.sep "${sampling_sep}"
        --planet.PA "${sampling_pa}"
        --planet.contrast "${sampling_contrast}"
        --fake.fileName ""
        --fake.subtractPlanet=false
        --p4.psfFile "${psf_file}"
        --p4.psfStampSize "${psf_stamp_size}"
        --p4.outputPSFModels=true
        --p4.psfFilter=false
        --p4.psfOutputPrefix p4PSF_
        --p4.psfSamplingMode refitDifference
        --p4.psfRefitContrast "${refit_contrast}"
        --p4.psfRadiiPerRegion 2
        --p4.psfSamplesPerRadius 4
        --p4.psfSampleAvoidRadius "${psf_sample_avoid_radius}"
        --p4Optimize.enabled=false
        --output.directory "${response_case}"
        --output.fileName finim.fits
        --output.exactFName=true
        --showTiming=true
    )
    run_timed "${response_case}" "${response_command[@]}"
    if [[ "${dry_run}" == false ]] && ! response_complete "${response_manifest}" "${refit_contrast}"; then
        printf 'Refit-difference stage did not publish the expected response manifest.\n' >&2
        exit 1
    fi
fi

fit_case="${experiment_dir}/refit_fit"
fit_summary="${fit_case}/summary.json"
if [[ "${dry_run}" == true ]]; then
    printf '\nMatched-response fitting requires completed products; omitted in dry-run mode.\n'
    exit 0
fi
if fit_complete "${fit_summary}" "${response_manifest}"; then
    printf '\n[refit_fit] completed fit exists; skipping.\n'
elif [[ -e "${fit_summary}" || -e "${fit_case}/amplitude.fits" ]]; then
    printf 'Incomplete refit-response fit exists; refusing to overwrite: %s\n' "${fit_case}" >&2
    exit 1
else
    python3 "${script_dir}/fit_p4_matched_response.py" \
        "${response_case}/finim.fits" \
        "${response_manifest}" \
        "${fit_case}" \
        --mode-fraction "${mode_fraction}" \
        --initial-separation "${exact_sep}" \
        --initial-pa "${exact_pa}" \
        --position-bound "${position_bound}" \
        --noise-exclusion-radius "${noise_exclusion_radius}" \
        --noise-min-radius "${noise_min_radius}" \
        --noise-max-radius "${noise_max_radius}" \
        --lambda-d "${lambda_d}"
fi

python3 "${script_dir}/compare_p4_refit_difference.py" \
    "${reference_experiment}" \
    "${response_manifest}" \
    "${fit_summary}" \
    "${experiment_dir}"

printf '\nP4 refit-difference experiment complete. Summary: %s\n' \
    "${experiment_dir}/refit_difference_comparison.md"
