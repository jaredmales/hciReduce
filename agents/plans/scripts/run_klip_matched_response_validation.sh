#!/usr/bin/env bash

set -euo pipefail

export LC_ALL=C

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "${script_dir}/../../.." && pwd)
roc_working_dir="${repo_root}/working/roc"
experiment_dir=${EXPERIMENT_DIR:-"${roc_working_dir}/klip_matched_response_$(date -u +%Y%m%dT%H%M%SZ)"}
base_config=${BASE_CONFIG:-"${script_dir}/klipReduce_afLepNaco_psf_response.conf"}
psf_file=${PSF_FILE:-/home/jrmales/Source/mxWork/NACO/AFLep/2011-10-21/out/psf_reg_median.fits}
klipreduce_bin=${KLIPREDUCE_BIN:-klipReduce}
psf_stamp_size=${PSF_STAMP_SIZE:-11}
psf_filter_min_good_fract=${PSF_FILTER_MIN_GOOD_FRACT:-1}
mode_counts=${MODE_COUNTS:-125,150,175,200,225,250,300,350}
response_case=${RESPONSE_CASE:-radial_ld_fixed16_filter}
planet_sep=${PLANET_SEP:-11.73795339222688}
planet_pa=${PLANET_PA:-262.1667995998323}
planet_contrast=${PLANET_CONTRAST:-0.004763925929356391}
candidate_avoid_radius=${CANDIDATE_AVOID_RADIUS:-5}
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

Run the first AF Lep/NACO KLIP matched-response validation. A science-only
control verifies that enabling the accepted 3.6-pixel radial, fixed-16-angle
response grid leaves the ordinary final cube unchanged. The response run
publishes production normalized-filter cubes and fits position/contrast from
those products without rerunning KLIP.

Environment overrides:
  EXPERIMENT_DIR             output directory (default: ${experiment_dir})
  BASE_CONFIG                response-compatible KLIP config (default: ${base_config})
  KLIPREDUCE_BIN             klipReduce executable (default: ${klipreduce_bin})
  PSF_FILE                   centered input PSF template (default: ${psf_file})
  PSF_STAMP_SIZE             response stamp width (default: ${psf_stamp_size})
  PSF_FILTER_MIN_GOOD_FRACT  minimum usable filter-stamp fraction (default: ${psf_filter_min_good_fract})
  MODE_COUNTS                comma-separated exact KL mode counts to fit (default: ${mode_counts})
  RESPONSE_CASE              response case run by the maintained driver (default: ${response_case})
  PLANET_SEP                 initial separation in pixels (default: ${planet_sep})
  PLANET_PA                  initial PA east of north (default: ${planet_pa})
  PLANET_CONTRAST            paired-response half-amplitude when applicable (default: ${planet_contrast})
  CANDIDATE_AVOID_RADIUS     paired-grid avoidance radius (default: ${candidate_avoid_radius})
  POSITION_BOUND             Cartesian fit half-width in pixels (default: ${position_bound})
  NOISE_EXCLUSION_RADIUS     matched-filter noise exclusion (default: ${noise_exclusion_radius})
  NOISE_MIN_RADIUS           inner noise-profile radius (default: ${noise_min_radius})
  NOISE_MAX_RADIUS           outer noise-profile radius (default: ${noise_max_radius})
  LAMBDA_D                   pixels per lambda/D (default: ${lambda_d})
  OMP_NUM_THREADS            OpenMP worker limit passed through to klipReduce

Example:
  OMP_NUM_THREADS=48 nohup $(basename "$0") > klip_matched_response_driver.log 2>&1 &
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

runner_command=("${script_dir}/run_klip_psf_response_experiment.sh")
if [[ "${dry_run}" == true ]]; then
    runner_command+=(--dry-run)
fi
runner_command+=(science_only "${response_case}")

EXPERIMENT_DIR="${experiment_dir}" \
BASE_CONFIG="${base_config}" \
KLIPREDUCE_BIN="${klipreduce_bin}" \
PSF_FILE="${psf_file}" \
PSF_STAMP_SIZE="${psf_stamp_size}" \
PSF_FILTER_MIN_GOOD_FRACT="${psf_filter_min_good_fract}" \
PLANET_SEP="${planet_sep}" \
PLANET_PA="${planet_pa}" \
PLANET_CONTRAST="${planet_contrast}" \
CANDIDATE_AVOID_RADIUS="${candidate_avoid_radius}" \
    "${runner_command[@]}"

if [[ "${dry_run}" == true ]]; then
    printf '\nDry run complete; matched-response fitting requires completed products.\n'
    exit 0
fi

case_directory="${experiment_dir}/${response_case}"
python3 - "${experiment_dir}/science_only/finim.fits" "${case_directory}/finim.fits" <<'PY'
import sys

import numpy as np
from astropy.io import fits

control = np.asarray(fits.getdata(sys.argv[1]))
response = np.asarray(fits.getdata(sys.argv[2]))
if control.shape != response.shape or not np.array_equal(control, response, equal_nan=True):
    raise SystemExit("enabling KLIP response measurement changed the ordinary final science cube")
print("KLIP science-only and response-enabled final cubes are elementwise identical")
PY

IFS=, read -r -a selected_modes <<< "${mode_counts}"
if ((${#selected_modes[@]} == 0)); then
    printf '%s\n' 'MODE_COUNTS must contain at least one exact KL mode count.' >&2
    exit 2
fi
for mode_count in "${selected_modes[@]}"; do
    fit_directory="${experiment_dir}/matched_fit_mode${mode_count}"
    summary_path="${fit_directory}/summary.json"
    if [[ -s "${summary_path}" && -s "${fit_directory}/surface.csv" ]]; then
        printf '\nCompleted KLIP matched-response fit exists; skipping: %s\n' "${summary_path}"
    elif [[ -s "${fit_directory}/failed.txt" ]]; then
        printf '\nRecorded failed KLIP matched-response fit exists; skipping: %s\n' "${fit_directory}/failed.txt"
    elif [[ -e "${fit_directory}" ]]; then
        printf 'Incomplete matched-response fit exists; refusing to overwrite: %s\n' "${fit_directory}" >&2
        exit 1
    else
        fit_command=(
            python3 "${script_dir}/fit_klip_matched_response.py"
            "${case_directory}"
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
        printf '%q ' "${fit_command[@]}" > "${fit_directory}/command.txt"
        printf '\n' >> "${fit_directory}/command.txt"
        set +e
        "${fit_command[@]}" 2>&1 | tee "${fit_directory}/run.log"
        fit_status=${PIPESTATUS[0]}
        set -e
        if ((fit_status != 0)); then
            printf 'exit_status=%d\n' "${fit_status}" > "${fit_directory}/failed.txt"
            printf 'KLIP matched-response fit for mode %s failed with status %d; continuing.\n' \
                "${mode_count}" "${fit_status}" >&2
        fi
    fi
done

python3 - "${experiment_dir}" <<'PY'
import csv
import json
import pathlib
import sys

experiment = pathlib.Path(sys.argv[1])
rows = []
failures = []
for path in sorted(experiment.glob("matched_fit_mode*/summary.json")):
    summary = json.loads(path.read_text(encoding="utf-8"))
    fit = summary["fit"]
    rows.append(
        {
            "mode_count": int(summary["mode_count"]),
            "status": str(fit["status"]),
            "separation": float(fit["separation"]),
            "position_angle": float(fit["position_angle"]),
            "contrast": float(fit["contrast"]),
            "contrast_standard_error": float(fit["contrast_standard_error"]),
            "snr": float(fit["snr"]),
        }
    )
for path in sorted(experiment.glob("matched_fit_mode*/failed.txt")):
    failures.append(path.parent.name.removeprefix("matched_fit_mode"))
if not rows and not failures:
    raise SystemExit("no completed or failed KLIP matched-response fits were found")
rows.sort(key=lambda row: row["mode_count"])
csv_path = experiment / "klip_matched_response_summary.csv"
with csv_path.open("w", encoding="utf-8", newline="") as stream:
    fieldnames = [
        "mode_count",
        "status",
        "separation",
        "position_angle",
        "contrast",
        "contrast_standard_error",
        "snr",
    ]
    writer = csv.DictWriter(stream, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
markdown_path = experiment / "klip_matched_response_summary.md"
with markdown_path.open("w", encoding="utf-8") as stream:
    stream.write("# KLIP matched-response summary\n\n")
    stream.write("| KL modes | Status | Separation | PA | Contrast | Contrast SE | S/N |\n")
    stream.write("|---:|---|---:|---:|---:|---:|---:|\n")
    for row in rows:
        stream.write(
            f"| {row['mode_count']} | {row['status']} | {row['separation']:.8g} | "
            f"{row['position_angle']:.8g} | {row['contrast']:.8g} | "
            f"{row['contrast_standard_error']:.8g} | {row['snr']:.8g} |\n"
        )
    if failures:
        stream.write("\nFailed mode fits (see each run log): " + ", ".join(failures) + "\n")
print(f"Wrote {markdown_path}")
PY

printf '\nKLIP matched-response validation complete: %s\n' \
    "${experiment_dir}/klip_matched_response_summary.md"
