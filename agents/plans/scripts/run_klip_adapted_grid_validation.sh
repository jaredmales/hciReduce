#!/usr/bin/env bash

set -euo pipefail

export LC_ALL=C
export PYTHONDONTWRITEBYTECODE=1

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "${script_dir}/../../.." && pwd)
roc_working_dir="${repo_root}/working/roc"
reference_experiment=${REFERENCE_EXPERIMENT:-"${roc_working_dir}/klip_matched_response_20260913T012040Z"}
finite_response_experiment=${FINITE_RESPONSE_EXPERIMENT:-"${roc_working_dir}/klip_finite_response_20260913T150323Z"}
candidate_experiment=${CANDIDATE_EXPERIMENT:-"${roc_working_dir}/klip_central_response_20260913T154659Z"}
reference_case_name=${REFERENCE_CASE:-radial_ld_fixed16_filter}
base_config=${BASE_CONFIG:-"${script_dir}/klipReduce_afLepNaco_psf_response.conf"}
psf_file=${PSF_FILE:-/home/jrmales/Source/mxWork/NACO/AFLep/2011-10-21/out/psf_reg_median.fits}
klipreduce_bin=${KLIPREDUCE_BIN:-klipReduce}
experiment_dir=${EXPERIMENT_DIR:-"${roc_working_dir}/klip_adapted_grid_$(date -u +%Y%m%dT%H%M%SZ)"}
mode_counts=${MODE_COUNTS:-}
planet_sep=${PLANET_SEP:-11.73795339222688}
planet_pa=${PLANET_PA:-262.1667995998323}
planet_contrast=${PLANET_CONTRAST:-0.004763925929356391}
candidate_avoid_radius=${CANDIDATE_AVOID_RADIUS:-5}
samples_per_radius=${SAMPLES_PER_RADIUS:-4}
amplitude_fraction=${AMPLITUDE_FRACTION:-1.0}
candidate_amplitude_fraction=${CANDIDATE_AMPLITUDE_FRACTION:-1.0}
minimum_support_fraction=${MINIMUM_SUPPORT_FRACTION:-1}
position_bound=${POSITION_BOUND:-1}
noise_exclusion_radius=${NOISE_EXCLUSION_RADIUS:-5}
noise_min_radius=${NOISE_MIN_RADIUS:-6}
noise_max_radius=${NOISE_MAX_RADIUS:-60}
lambda_d=${LAMBDA_D:-3.6}
dry_run=false

usage()
{
    cat <<EOF
Usage: $(basename "$0") [--dry-run]

Build a complete candidate-avoiding paired-refit KLIP response grid on the
accepted 3.6-pixel radial nodes. The driver measures positive and negative
perturbations, rotates and averages the clear responses at each radius,
publishes normalized-filter products, and fits every configured KL mode.

Environment overrides:
  REFERENCE_EXPERIMENT       completed frozen-response experiment
  FINITE_RESPONSE_EXPERIMENT completed one-sided response experiment used by the measurement analyzer
  CANDIDATE_EXPERIMENT       completed central-response experiment containing the candidate oracle
  REFERENCE_CASE             frozen-response case (default: ${reference_case_name})
  EXPERIMENT_DIR             output directory (default: ${experiment_dir})
  BASE_CONFIG                maintained AF Lep/NACO KLIP configuration
  KLIPREDUCE_BIN             klipReduce executable
  PSF_FILE                   centered input PSF used for perturbations
  MODE_COUNTS                comma-separated KL modes; default is reference NMODES
  PLANET_SEP, PLANET_PA      known candidate coordinates
  PLANET_CONTRAST            fiducial perturbation contrast
  CANDIDATE_AVOID_RADIUS     excluded distance around the candidate (default: ${candidate_avoid_radius})
  SAMPLES_PER_RADIUS         requested uniform angular samples (default: ${samples_per_radius})
  AMPLITUDE_FRACTION         epsilon/PLANET_CONTRAST (default: ${amplitude_fraction})
  CANDIDATE_AMPLITUDE_FRACTION
                             existing candidate-oracle fraction (default: ${candidate_amplitude_fraction})
  MINIMUM_SUPPORT_FRACTION   matched-filter usable-stamp fraction (default: ${minimum_support_fraction})
  POSITION_BOUND             fit half-width in pixels (default: ${position_bound})
  NOISE_EXCLUSION_RADIUS     fitted-source noise exclusion (default: ${noise_exclusion_radius})
  NOISE_MIN_RADIUS           inner noise radius (default: ${noise_min_radius})
  NOISE_MAX_RADIUS           outer noise radius (default: ${noise_max_radius})
  LAMBDA_D                   pixels per lambda/D (default: ${lambda_d})
  OMP_NUM_THREADS            OpenMP worker limit passed through to klipReduce

With the defaults, three of 60 requested locations are excluded around AF Lep b,
leaving 57 samples and 114 KLIP reductions.

Example:
  OMP_NUM_THREADS=48 nohup $(basename "$0") > klip_adapted_grid_driver.log 2>&1 &
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

[[ -d "${reference_experiment}/${reference_case_name}" ]] || {
    printf 'Reference response case is not readable: %s\n' \
        "${reference_experiment}/${reference_case_name}" >&2
    exit 1
}
for required_product in \
    "${finite_response_experiment}/empirical_response.fits" \
    "${candidate_experiment}/klip_central_response.json" \
    "${candidate_experiment}/central_response.fits"; do
    [[ -r "${required_product}" ]] || {
        printf 'Required response product is missing: %s\n' "${required_product}" >&2
        exit 1
    }
done
python3 -c 'import astropy.io.fits, numpy' >/dev/null 2>&1 || {
    printf '%s\n' 'The Python astropy and numpy packages are required.' >&2
    exit 1
}

reference_case="${reference_experiment}/${reference_case_name}"
reference_response="${reference_case}/finim_outputs/klipPSF_mode000_radial_response.fits"
[[ -r "${reference_response}" ]] || {
    printf 'Reference radial response is missing: %s\n' "${reference_response}" >&2
    exit 1
}
if [[ -z "${mode_counts}" ]]; then
    mode_counts=$(python3 - "${reference_case}/finim.fits" <<'PY'
import sys
from astropy.io import fits
print(str(fits.getheader(sys.argv[1])["NMODES"]))
PY
)
fi

sample_specs=$(python3 - "${reference_response}" "${planet_sep}" "${planet_pa}" \
    "${candidate_avoid_radius}" "${samples_per_radius}" <<'PY'
import math
import re
import sys

from astropy.io import fits

path, planet_sep_text, planet_pa_text, avoid_text, count_text = sys.argv[1:]
planet_sep = float(planet_sep_text)
planet_pa = float(planet_pa_text)
avoid = float(avoid_text)
count = int(count_text)
if not all(math.isfinite(value) for value in (planet_sep, planet_pa, avoid)):
    raise SystemExit("planet coordinates and avoidance radius must be finite")
if planet_sep < 0 or avoid < 0 or count <= 0:
    raise SystemExit("planet separation, avoidance radius, or angular sample count is invalid")
radii = [float(token) for token in str(fits.getheader(path)["KLIP PSF SAMPLE RADII"]).split(",")]
specs = []
counts = []
for radius_index, radius in enumerate(radii):
    retained = 0
    for angle_index in range(count):
        pa = 360.0 * angle_index / count
        delta = math.radians(pa - planet_pa)
        distance = math.sqrt(
            radius * radius + planet_sep * planet_sep
            - 2 * radius * planet_sep * math.cos(delta)
        )
        if distance < avoid:
            continue
        label = f"r{radius_index:02d}_a{angle_index:02d}"
        specs.append(f"{label}:{radius:.17g}:{pa:.17g}:clear")
        retained += 1
    if retained == 0:
        raise SystemExit(f"candidate avoidance removed every sample at radius {radius}")
    counts.append(retained)
if not specs or any(re.fullmatch(r"[A-Za-z0-9_]+:[^;]+", spec) is None for spec in specs):
    raise SystemExit("failed to construct valid adapted-grid sample specifications")
print(";".join(specs))
PY
)

central_command=("${script_dir}/run_klip_central_response_validation.sh")
if [[ "${dry_run}" == true ]]; then
    central_command+=(--dry-run)
fi
EXPERIMENT_DIR="${experiment_dir}" \
REFERENCE_EXPERIMENT="${reference_experiment}" \
FINITE_RESPONSE_EXPERIMENT="${finite_response_experiment}" \
REFERENCE_CASE="${reference_case_name}" \
BASE_CONFIG="${base_config}" \
PSF_FILE="${psf_file}" \
KLIPREDUCE_BIN="${klipreduce_bin}" \
MODE_COUNTS="${mode_counts}" \
PLANET_SEP="${planet_sep}" \
PLANET_PA="${planet_pa}" \
PLANET_CONTRAST="${planet_contrast}" \
CANDIDATE_AVOID_RADIUS="${candidate_avoid_radius}" \
AMPLITUDE_FRACTIONS="${amplitude_fraction}" \
SAMPLE_SPECS="${sample_specs}" \
    "${central_command[@]}"

if [[ "${dry_run}" == true ]]; then
    printf '\nDry run complete; adapted-grid construction requires the paired reductions.\n'
    exit 0
fi

adapted_case="${experiment_dir}/adapted_grid"
adapted_summary="${adapted_case}/adapted_grid_summary.json"
adapted_manifest="${adapted_case}/finim_outputs/klipPSF_manifest.fits"
if python3 - "${adapted_summary}" "${adapted_manifest}" <<'PY'
import json
import pathlib
import sys

from astropy.io import fits

summary_path, manifest_path = map(pathlib.Path, sys.argv[1:])
try:
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    header = fits.getheader(manifest_path)
    complete = (
        int(summary.get("schema", 0)) == 1
        and int(header.get("KLIP PSF COMPLETE", 0)) == 1
        and str(header.get("KLIP PSF RESPONSE METHOD", "")).strip() == "PAIRED_CENTRAL_REFIT"
    )
except Exception:
    complete = False
raise SystemExit(0 if complete else 1)
PY
then
    printf '\n[adapted_grid] completed products exist; skipping.\n'
elif [[ -e "${adapted_case}" ]]; then
    printf 'Incomplete adapted-grid products exist; refusing to overwrite: %s\n' "${adapted_case}" >&2
    exit 1
else
    build_command=(
        python3 "${script_dir}/build_klip_adapted_grid.py"
        "${reference_case}"
        "${experiment_dir}"
        "${adapted_case}"
        --candidate-experiment "${candidate_experiment}"
        --amplitude-fraction "${amplitude_fraction}"
        --candidate-amplitude-fraction "${candidate_amplitude_fraction}"
        --minimum-radius "${noise_min_radius}"
        --maximum-radius "${noise_max_radius}"
        --minimum-support-fraction "${minimum_support_fraction}"
    )
    shell_join "${build_command[@]}" > "${experiment_dir}/adapted_grid_command.txt"
    "${build_command[@]}" 2>&1 | tee "${experiment_dir}/adapted_grid.log"
fi

IFS=, read -r -a selected_modes <<< "${mode_counts}"
for mode_count in "${selected_modes[@]}"; do
    fit_directory="${experiment_dir}/adapted_fit_mode${mode_count}"
    if [[ -s "${fit_directory}/summary.json" && -s "${fit_directory}/surface.csv" ]]; then
        printf '\nCompleted adapted-grid fit exists; skipping: %s\n' "${fit_directory}/summary.json"
    elif [[ -s "${fit_directory}/failed.txt" ]]; then
        printf '\nRecorded failed adapted-grid fit exists; skipping: %s\n' "${fit_directory}/failed.txt"
    elif [[ -e "${fit_directory}" ]]; then
        printf 'Incomplete adapted-grid fit exists; refusing to overwrite: %s\n' "${fit_directory}" >&2
        exit 1
    else
        fit_command=(
            python3 "${script_dir}/fit_klip_matched_response.py"
            "${adapted_case}"
            "${fit_directory}"
            --mode-count "${mode_count}"
            --initial-separation "${planet_sep}"
            --initial-pa "${planet_pa}"
            --position-bound "${position_bound}"
            --noise-exclusion-radius "${noise_exclusion_radius}"
            --noise-min-radius "${noise_min_radius}"
            --noise-max-radius "${noise_max_radius}"
            --lambda-d "${lambda_d}"
        )
        mkdir -p "${fit_directory}"
        shell_join "${fit_command[@]}" > "${fit_directory}/command.txt"
        set +e
        "${fit_command[@]}" 2>&1 | tee "${fit_directory}/run.log"
        fit_status=${PIPESTATUS[0]}
        set -e
        if ((fit_status != 0)); then
            printf 'exit_status=%d\n' "${fit_status}" > "${fit_directory}/failed.txt"
            printf 'KLIP adapted-grid fit for mode %s failed with status %d; continuing.\n' \
                "${mode_count}" "${fit_status}" >&2
        fi
    fi
done

python3 - "${experiment_dir}" "${planet_contrast}" <<'PY'
import csv
import json
import pathlib
import sys

experiment = pathlib.Path(sys.argv[1])
fiducial = float(sys.argv[2])
rows = []
failures = []
for path in sorted(experiment.glob("adapted_fit_mode*/summary.json")):
    summary = json.loads(path.read_text(encoding="utf-8"))
    fit = summary["fit"]
    rows.append(
        {
            "mode_count": int(summary["mode_count"]),
            "status": str(fit["status"]),
            "separation": float(fit["separation"]),
            "position_angle": float(fit["position_angle"]),
            "contrast": float(fit["contrast"]),
            "contrast_fraction_of_fiducial": float(fit["contrast"]) / fiducial,
            "contrast_standard_error": float(fit["contrast_standard_error"]),
            "snr": float(fit["snr"]),
        }
    )
for path in sorted(experiment.glob("adapted_fit_mode*/failed.txt")):
    failures.append(path.parent.name.removeprefix("adapted_fit_mode"))
if not rows and not failures:
    raise SystemExit("no completed or failed adapted-grid fits were found")
rows.sort(key=lambda row: row["mode_count"])
csv_path = experiment / "klip_adapted_grid_fit_summary.csv"
with csv_path.open("w", encoding="utf-8", newline="") as stream:
    fieldnames = [
        "mode_count",
        "status",
        "separation",
        "position_angle",
        "contrast",
        "contrast_fraction_of_fiducial",
        "contrast_standard_error",
        "snr",
    ]
    writer = csv.DictWriter(stream, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
markdown_path = experiment / "klip_adapted_grid_fit_summary.md"
with markdown_path.open("w", encoding="utf-8") as stream:
    stream.write("# KLIP adapted-grid matched-response summary\n\n")
    stream.write("| KL modes | Status | Separation | PA | Contrast | Fiducial fraction | Contrast SE | S/N |\n")
    stream.write("|---:|---|---:|---:|---:|---:|---:|---:|\n")
    for row in rows:
        stream.write(
            f"| {row['mode_count']} | {row['status']} | {row['separation']:.8g} | "
            f"{row['position_angle']:.8g} | {row['contrast']:.8g} | "
            f"{row['contrast_fraction_of_fiducial']:.8g} | "
            f"{row['contrast_standard_error']:.8g} | {row['snr']:.8g} |\n"
        )
    if failures:
        stream.write("\nFailed mode fits (see each run log): " + ", ".join(failures) + "\n")
print(f"Wrote {markdown_path}")
PY

printf '\nKLIP adapted-grid validation complete: %s\n' \
    "${experiment_dir}/klip_adapted_grid_fit_summary.md"
