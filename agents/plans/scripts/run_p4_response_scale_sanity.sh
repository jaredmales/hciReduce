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
experiment_dir=${EXPERIMENT_DIR:-"${roc_working_dir}/p4_response_scale_sanity_$(date -u +%Y%m%dT%H%M%SZ)"}
mode_fraction=${MODE_FRACTION:-0.15}
position_bound=${POSITION_BOUND:-1}
psf_stamp_size=${PSF_STAMP_SIZE:-11}
psf_sample_avoid_radius=${PSF_SAMPLE_AVOID_RADIUS:-5}
noise_exclusion_radius=${NOISE_EXCLUSION_RADIUS:-5}
noise_min_radius=${NOISE_MIN_RADIUS:-6}
noise_max_radius=${NOISE_MAX_RADIUS:-60}
lambda_d=${LAMBDA_D:-3.6}
over_subtraction_factor=${OVER_SUBTRACTION_FACTOR:-4}
dry_run=false

usage()
{
    cat <<EOF
Usage: $(basename "$0") [--dry-run]

Run two independent sanity checks for the P4 matched-response flux scale:

  1. Refit P4 after subtracting the original matched-response contrast and,
     separately, four times the exact negative-fit contrast. The first case
     should retain a positive planet; the second should be strongly negative.
  2. Multiply the response-input PSF by the exact fitted contrast before
     calculating a new sparse response, fit in those scaled-template units,
     and convert the result back to physical contrast.

Completed stages are skipped, so EXPERIMENT_DIR can be resumed.

Environment overrides:
  REFERENCE_EXPERIMENT       completed matched-response experiment
                            (default: ${reference_experiment})
  EXPERIMENT_DIR             output directory (default: ${experiment_dir})
  BASE_CONFIG                standard ROC P4 configuration (default: ${base_config})
  P4REDUCE_BIN               p4Reduce executable (default: ${p4reduce_bin})
  PSF_FILE                   unscaled centered P4-input PSF (default: ${psf_file})
  MODE_FRACTION              sole P4 mode to calculate (default: ${mode_fraction})
  OVER_SUBTRACTION_FACTOR    exact-contrast multiplier, greater than one
                            (default: ${over_subtraction_factor})
  PSF_SAMPLE_AVOID_RADIUS    sparse detector-sample exclusion
                            (default: ${psf_sample_avoid_radius})
  OMP_NUM_THREADS            OpenMP worker limit passed through to p4Reduce

Example:
  OMP_NUM_THREADS=48 nohup $(basename "$0") > p4_response_scale_sanity_driver.log 2>&1 &
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
    if ((command_status != 0)); then
        printf '[%s] failed with status %d.\n' "$(basename "${stage_directory}")" "${command_status}" >&2
    fi
    return "${command_status}"
}

final_reduction_complete()
{
    local final_image=$1
    local expected_contrast=$2
    [[ -s "${final_image}" ]] || return 1
    python3 - "${final_image}" "${mode_fraction}" "${expected_contrast}" <<'PY'
import sys

import numpy as np
from astropy.io import fits

path, requested_text, expected_contrast_text = sys.argv[1:]
requested = float(requested_text)
expected_contrast = float(expected_contrast_text)
try:
    header = fits.getheader(path)
    data = np.asarray(fits.getdata(path))
    modes = np.asarray([float(value) for value in str(header["P4 MODE FRACTIONS"]).split(",")])
    contrasts = np.asarray([float(value) for value in str(header["PLANETCONT"]).split(",")])
    matches = np.flatnonzero(
        np.abs(modes - requested)
        <= 8
        * np.finfo(np.float32).eps
        * np.maximum.reduce((np.ones_like(modes), np.abs(modes), np.full_like(modes, abs(requested))))
    )
    complete = (
        modes.size == 1
        and matches.size == 1
        and contrasts.size == 1
        and np.isclose(contrasts[0], expected_contrast, rtol=2e-5, atol=0)
        and bool(str(header.get("FAKEFILE", "")).strip())
        and data.size > 0
        and np.any(np.isfinite(data))
    )
except Exception:
    complete = False
raise SystemExit(0 if complete else 1)
PY
}

scaled_response_complete()
{
    local case_directory=$1
    local scaled_template=$2
    local manifest="${case_directory}/finim_outputs/p4PSF_manifest.fits"
    [[ -s "${manifest}" ]] || return 1
    python3 - "${manifest}" "${mode_fraction}" "${scaled_template}" <<'PY'
import pathlib
import sys
import numpy as np
from astropy.io import fits

path, requested_text, template_text = sys.argv[1:]
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
        and int(header.get("P4 PSF PRODUCT SCHEMA", 0)) == 6
        and str(header.get("P4 PSF SPATIAL MODEL", "")).strip() == "REGION_TARGET_RADIAL_LINEAR"
        and str(header.get("P4 PSF NORMALIZATION", "")).strip() == "STORED"
        and pathlib.Path(str(header.get("P4 PSF TEMPLATE", ""))).resolve() == pathlib.Path(template_text).resolve()
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
for required_option in --fake.subtractPlanet --p4.psfFile --p4.psfSamplingMode --p4Optimize.enabled; do
    if [[ "${help_text}" != *"${required_option}"* ]]; then
        printf 'p4Reduce does not expose required option %s; build this checkout first.\n' "${required_option}" >&2
        exit 1
    fi
done

reference_experiment=$(absolute_directory "${reference_experiment}")
reference_original="${reference_experiment}/sparse_response/finim.fits"
reference_manifest="${reference_experiment}/sparse_response/finim_outputs/p4PSF_manifest.fits"
reference_fit_summary="${reference_experiment}/sparse_fit/summary.json"
optimizer_summary="${reference_experiment}/exact_optimizer/finim_outputs/p4Negative_summary.yaml"
for required_product in "${reference_original}" "${reference_manifest}" "${reference_fit_summary}" \
    "${optimizer_summary}" "${reference_experiment}/signal_free_oracle/finim.fits"; do
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

readarray -t response_planet < <(python3 - "${reference_fit_summary}" <<'PY'
import json
import math
import pathlib
import sys

fit = json.loads(pathlib.Path(sys.argv[1]).read_text(encoding="utf-8"))["fit"]
values = [float(fit["separation"]), float(fit["position_angle"]), float(fit["contrast"])]
if fit.get("status") != "converged" or values[0] < 0 or values[2] <= 0 or not all(map(math.isfinite, values)):
    raise SystemExit("reference matched-response fit is invalid")
for value in values:
    print(f"{value:.17g}")
PY
)
response_sep=${response_planet[0]}
response_pa=${response_planet[1]}
response_contrast=${response_planet[2]}
readarray -t sampling_planet < <(python3 - "${reference_original}" <<'PY'
import math
import sys

from astropy.io import fits

header = fits.getheader(sys.argv[1])
values = []
for key in ("PLANETSEP", "PLANETPA", "PLANETCONT"):
    entries = [float(value) for value in str(header[key]).split(",")]
    if len(entries) != 1:
        raise SystemExit(f"reference sparse response does not have one {key} value")
    values.append(entries[0])
if values[0] < 0 or values[2] < 0 or not all(map(math.isfinite, values)):
    raise SystemExit("reference sparse-response planet metadata is invalid")
for value in values:
    print(f"{value:.17g}")
PY
)
sampling_sep=${sampling_planet[0]}
sampling_pa=${sampling_planet[1]}
sampling_contrast=${sampling_planet[2]}
over_contrast=$(python3 - "${exact_contrast}" "${over_subtraction_factor}" <<'PY'
import sys
contrast, factor = map(float, sys.argv[1:])
if contrast <= 0 or factor <= 1:
    raise SystemExit("invalid over-subtraction request")
print(f"{contrast * factor:.17g}")
PY
)

mkdir -p "${experiment_dir}"
experiment_dir=$(absolute_directory "${experiment_dir}")
scaled_psf="${experiment_dir}/psf_scaled_to_exact_contrast.fits"
python3 - "${psf_file}" "${scaled_psf}" "${exact_contrast}" <<'PY'
import pathlib
import sys
import warnings

import numpy as np
from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning

source_text, output_text, scale_text = sys.argv[1:]
warnings.simplefilter("ignore", VerifyWarning)
source = pathlib.Path(source_text).resolve()
output = pathlib.Path(output_text)
scale = float(scale_text)
header = fits.getheader(source)
data = np.asarray(fits.getdata(source), dtype=np.float32)
expected = data * np.float32(scale)
if output.exists():
    existing = np.asarray(fits.getdata(output), dtype=np.float32)
    existing_header = fits.getheader(output)
    valid = (
        np.array_equal(existing, expected)
        and abs(float(existing_header.get("HCIREDUCE PSF SCALE", float("nan"))) - scale) <= 1e-12 * scale
    )
    if not valid:
        raise SystemExit(f"existing scaled PSF does not match the requested product: {output}")
else:
    header["HIERARCH HCIREDUCE PSF SCALE"] = (scale, "multiplicative sanity-check scale")
    fits.PrimaryHDU(data=expected, header=header).writeto(output)
PY

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
    "scaled_psf=${scaled_psf}" \
    "scaled_psf_sha256=$(sha256sum "${scaled_psf}" | awk '{ print $1 }')" \
    "reference_experiment=${reference_experiment}" \
    "reference_manifest_sha256=$(sha256sum "${reference_manifest}" | awk '{ print $1 }')" \
    "reference_fit_summary_sha256=$(sha256sum "${reference_fit_summary}" | awk '{ print $1 }')" \
    "optimizer_summary_sha256=$(sha256sum "${optimizer_summary}" | awk '{ print $1 }')" \
    "exact_sep=${exact_sep}" \
    "exact_pa=${exact_pa}" \
    "exact_contrast=${exact_contrast}" \
    "response_sep=${response_sep}" \
    "response_pa=${response_pa}" \
    "response_contrast=${response_contrast}" \
    "sampling_sep=${sampling_sep}" \
    "sampling_pa=${sampling_pa}" \
    "sampling_contrast=${sampling_contrast}" \
    "over_subtraction_factor=${over_subtraction_factor}" \
    "over_contrast=${over_contrast}" \
    "mode_fraction=${mode_fraction}" \
    "omp_num_threads=${OMP_NUM_THREADS:-unlimited}")
settings_file="${experiment_dir}/settings.txt"
if [[ -e "${settings_file}" ]]; then
    [[ "$(<"${settings_file}")" == "${settings_content}" ]] || {
        printf 'Experiment settings differ from saved settings: %s\n' "${settings_file}" >&2
        printf '%s\n' 'Choose a new EXPERIMENT_DIR or restore the original overrides.' >&2
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
printf 'Response-fit subtraction: sep=%s PA=%s contrast=%s\n' \
    "${response_sep}" "${response_pa}" "${response_contrast}"
printf 'Exact planet: sep=%s PA=%s contrast=%s\n' "${exact_sep}" "${exact_pa}" "${exact_contrast}"
printf 'Over-subtraction contrast (%sx exact): %s\n' "${over_subtraction_factor}" "${over_contrast}"
printf 'Scaled response PSF: %s\n' "${scaled_psf}"

response_subtraction_case="${experiment_dir}/response_fit_subtraction"
response_subtraction_image="${response_subtraction_case}/finim.fits"
if final_reduction_complete "${response_subtraction_image}" "${response_contrast}"; then
    printf '\n[response_fit_subtraction] completed reduction exists; skipping.\n'
elif [[ -e "${response_subtraction_case}/run.log" || -e "${response_subtraction_image}" ]]; then
    printf 'Incomplete response-fit subtraction exists; refusing to overwrite: %s\n' \
        "${response_subtraction_case}" >&2
    exit 1
else
    response_subtraction_command=(
        "${p4reduce_bin}"
        --config "${base_config}"
        --input.imSize 256
        --p4.modeFractions "${mode_fraction}"
        --planet.sep "${response_sep}"
        --planet.PA "${response_pa}"
        --planet.contrast "${response_contrast}"
        --fake.method single
        --fake.fileName "${psf_file}"
        --fake.subtractPlanet=true
        --p4.localStampSize 0
        --p4.psfFile ""
        --p4.outputPSFModels=false
        --p4.psfFilter=false
        --p4Optimize.enabled=false
        --output.directory "${response_subtraction_case}"
        --output.fileName finim.fits
        --output.exactFName=true
        --showTiming=true
    )
    run_timed "${response_subtraction_case}" "${response_subtraction_command[@]}"
    if [[ "${dry_run}" == false ]] && \
        ! final_reduction_complete "${response_subtraction_image}" "${response_contrast}"; then
        printf 'Response-fit subtraction did not publish the expected final image.\n' >&2
        exit 1
    fi
fi

over_subtraction_case="${experiment_dir}/four_times_exact_subtraction"
over_subtraction_image="${over_subtraction_case}/finim.fits"
if final_reduction_complete "${over_subtraction_image}" "${over_contrast}"; then
    printf '\n[four_times_exact_subtraction] completed reduction exists; skipping.\n'
elif [[ -e "${over_subtraction_case}/run.log" || -e "${over_subtraction_image}" ]]; then
    printf 'Incomplete over-subtraction exists; refusing to overwrite: %s\n' "${over_subtraction_case}" >&2
    exit 1
else
    over_subtraction_command=(
        "${p4reduce_bin}"
        --config "${base_config}"
        --input.imSize 256
        --p4.modeFractions "${mode_fraction}"
        --planet.sep "${exact_sep}"
        --planet.PA "${exact_pa}"
        --planet.contrast "${over_contrast}"
        --fake.method single
        --fake.fileName "${psf_file}"
        --fake.subtractPlanet=true
        --p4.localStampSize 0
        --p4.psfFile ""
        --p4.outputPSFModels=false
        --p4.psfFilter=false
        --p4Optimize.enabled=false
        --output.directory "${over_subtraction_case}"
        --output.fileName finim.fits
        --output.exactFName=true
        --showTiming=true
    )
    run_timed "${over_subtraction_case}" "${over_subtraction_command[@]}"
    if [[ "${dry_run}" == false ]] && ! final_reduction_complete "${over_subtraction_image}" "${over_contrast}"; then
        printf 'Over-subtraction did not publish the expected final image.\n' >&2
        exit 1
    fi
fi

scaled_response_case="${experiment_dir}/scaled_response"
scaled_manifest="${scaled_response_case}/finim_outputs/p4PSF_manifest.fits"
if scaled_response_complete "${scaled_response_case}" "${scaled_psf}"; then
    printf '\n[scaled_response] completed response exists; skipping.\n'
elif [[ -e "${scaled_response_case}/run.log" || -e "${scaled_manifest}" ]]; then
    printf 'Incomplete scaled response exists; refusing to overwrite: %s\n' "${scaled_response_case}" >&2
    exit 1
else
    scaled_response_command=(
        "${p4reduce_bin}"
        --config "${base_config}"
        --input.imSize 256
        --p4.modeFractions "${mode_fraction}"
        --planet.sep "${sampling_sep}"
        --planet.PA "${sampling_pa}"
        --planet.contrast "${sampling_contrast}"
        --fake.fileName ""
        --fake.subtractPlanet=false
        --p4.psfFile "${scaled_psf}"
        --p4.psfStampSize "${psf_stamp_size}"
        --p4.outputPSFModels=true
        --p4.psfFilter=true
        --p4.psfOutputPrefix p4PSF_
        --p4.psfSamplingMode detectorLocal
        --p4.psfRadiiPerRegion 2
        --p4.psfSamplesPerRadius 4
        --p4.psfSampleAvoidRadius "${psf_sample_avoid_radius}"
        --p4Optimize.enabled=false
        --output.directory "${scaled_response_case}"
        --output.fileName finim.fits
        --output.exactFName=true
        --showTiming=true
    )
    run_timed "${scaled_response_case}" "${scaled_response_command[@]}"
    if [[ "${dry_run}" == false ]] && ! scaled_response_complete "${scaled_response_case}" "${scaled_psf}"; then
        printf 'Scaled response did not publish the expected response manifest.\n' >&2
        exit 1
    fi
fi

scaled_fit="${experiment_dir}/scaled_fit"
scaled_fit_summary="${scaled_fit}/summary.json"
if [[ "${dry_run}" == true ]]; then
    printf '\nScaled response fit and comparison require completed reductions; omitted in dry-run mode.\n'
    exit 0
fi
if fit_complete "${scaled_fit_summary}" "${scaled_manifest}"; then
    printf '\n[scaled_fit] completed fit exists; skipping.\n'
elif [[ -e "${scaled_fit_summary}" || -e "${scaled_fit}/amplitude.fits" ]]; then
    printf 'Incomplete scaled fit exists; refusing to overwrite: %s\n' "${scaled_fit}" >&2
    exit 1
else
    python3 "${script_dir}/fit_p4_matched_response.py" \
        "${reference_original}" \
        "${scaled_manifest}" \
        "${scaled_fit}" \
        --mode-fraction "${mode_fraction}" \
        --initial-separation "${exact_sep}" \
        --initial-pa "${exact_pa}" \
        --position-bound "${position_bound}" \
        --noise-exclusion-radius "${noise_exclusion_radius}" \
        --noise-min-radius "${noise_min_radius}" \
        --noise-max-radius "${noise_max_radius}" \
        --lambda-d "${lambda_d}"
fi

python3 "${script_dir}/compare_p4_response_scale_sanity.py" \
    "${reference_experiment}" \
    "${response_subtraction_image}" \
    "${over_subtraction_image}" \
    "${scaled_manifest}" \
    "${scaled_fit_summary}" \
    "${experiment_dir}" \
    --mode-fraction "${mode_fraction}" \
    --template-scale "${exact_contrast}" \
    --over-subtraction-factor "${over_subtraction_factor}"

printf '\nP4 response-scale sanity experiment complete. Summary: %s\n' \
    "${experiment_dir}/response_scale_sanity.md"
