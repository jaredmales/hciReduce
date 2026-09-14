#!/usr/bin/env bash

set -euo pipefail

export LC_ALL=C

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "${script_dir}/../../.." && pwd)
roc_working_dir="${repo_root}/working/roc"
source_experiment=${SOURCE_EXPERIMENT:-"${roc_working_dir}/klip_cpp_response_20260913T222913Z"}
science_file=${SCIENCE_FILE:-"${source_experiment}/science_only/finim.fits"}
default_response_manifest="${source_experiment}/radial_ld_refit4_filter/finim_outputs/klipPSF_manifest.fits"
response_manifest=${RESPONSE_MANIFEST:-"${default_response_manifest}"}
exact_response_manifest=${EXACT_RESPONSE_MANIFEST:-}
experiment_dir=${EXPERIMENT_DIR:-"${roc_working_dir}/klip_hciAnalyze_filters_$(date -u +%Y%m%dT%H%M%SZ)"}
hcianalyze_bin=${HCIANALYZE_BIN:-hciAnalyze}
lambda_d=${LAMBDA_D:-3.6}
gaussian_fwhms=${GAUSSIAN_FWHMS:-3.6}
planet_sep=${PLANET_SEP:-11.73795339222688}
planet_pa=${PLANET_PA:-262.1667995998323}
planet_contrast=${PLANET_CONTRAST:-0.004763925929356391}
planet_radius=${PLANET_RADIUS:-5}
snr_min_radius=${SNR_MIN_RADIUS:-6}
snr_max_radius=${SNR_MAX_RADIUS:-60}
snr_aperture_radius=${SNR_APERTURE_RADIUS:-2}
dry_run=false

usage()
{
    cat <<EOF
Usage: $(basename "$0") [--dry-run]

Apply hciAnalyze to one KLIP science cube with no spatial filter, Gaussian
low-pass smoothing, and the persisted sparse KLIP response. Each case receives
an identical copy of the input cube so its finim_snr.fits output is retained.

Environment overrides:
  SOURCE_EXPERIMENT   completed KLIP response experiment
                      (default: ${source_experiment})
  SCIENCE_FILE        standalone unfiltered KLIP cube
                      (default: ${science_file})
  RESPONSE_MANIFEST   complete sparse-response manifest; its sibling response
                      and validity cubes must remain beside it
                      (default: ${response_manifest})
  EXACT_RESPONSE_MANIFEST
                      optional exact-azimuthal or per-pixel response manifest;
                      when set, add an exact_response comparison case
  EXPERIMENT_DIR      comparison output directory (default: ${experiment_dir})
  HCIANALYZE_BIN      hciAnalyze executable (default: ${hcianalyze_bin})
  LAMBDA_D            pixels per lambda/D (default: ${lambda_d})
  GAUSSIAN_FWHMS      comma-separated low-pass Gaussian FWHMs in pixels
                      (default: ${gaussian_fwhms})
  PLANET_SEP          candidate separation in pixels (default: ${planet_sep})
  PLANET_PA           candidate PA east of north (default: ${planet_pa})
  PLANET_CONTRAST     reported candidate contrast (default: ${planet_contrast})
  PLANET_RADIUS       candidate noise-exclusion radius (default: ${planet_radius})
  SNR_MIN_RADIUS      inner radial-noise bound (default: ${snr_min_radius})
  SNR_MAX_RADIUS      outer radial-noise bound (default: ${snr_max_radius})
  SNR_APERTURE_RADIUS signal aperture radius (default: ${snr_aperture_radius})

Examples:
  $(basename "$0")
  GAUSSIAN_FWHMS=2.4,3.0,3.6,4.2 $(basename "$0")
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

command -v "${hcianalyze_bin}" >/dev/null 2>&1 || {
    printf 'hciAnalyze executable was not found: %s\n' "${hcianalyze_bin}" >&2
    exit 1
}
science_file=$(absolute_file "${science_file}")
response_manifest=$(absolute_file "${response_manifest}")
if [[ -n "${exact_response_manifest}" ]]; then
    exact_response_manifest=$(absolute_file "${exact_response_manifest}")
fi

if [[ -e "${experiment_dir}" ]]; then
    printf 'Experiment output already exists; refusing to overwrite: %s\n' "${experiment_dir}" >&2
    exit 1
fi

common_options=(
    --lambdaD "${lambda_d}"
    --planet.sep "${planet_sep}"
    --planet.PA "${planet_pa}"
    --planet.contrast "${planet_contrast}"
    --planet.R "${planet_radius}"
    --snr.minRad "${snr_min_radius}"
    --snr.maxRad "${snr_max_radius}"
    --snr.apertureR "${snr_aperture_radius}"
    --filter.hpfGaussFW 0
    --diagnostics=true
)

run_case()
{
    local case_name=$1
    local filter_description=$2
    shift 2
    local filter_options=("$@")
    local case_directory="${experiment_dir}/${case_name}"
    local case_file="${case_directory}/finim.fits"
    local command_line=(
        "${hcianalyze_bin}"
        --file "${case_file}"
        "${common_options[@]}"
        "${filter_options[@]}"
    )

    printf '\n[%s] %s\n' "${case_name}" "${filter_description}"
    shell_join "${command_line[@]}"
    if [[ "${dry_run}" == true ]]; then
        return
    fi

    mkdir -p "${case_directory}"
    cp -- "${science_file}" "${case_file}"
    shell_join "${command_line[@]}" > "${case_directory}/command.txt"
    printf 'filter=%s\n' "${filter_description}" > "${case_directory}/case.env"
    "${command_line[@]}" > "${case_directory}/results.txt" 2> "${case_directory}/diagnostics.txt"
    [[ -s "${case_directory}/finim_snr.fits" ]] || {
        printf '[%s] hciAnalyze did not publish finim_snr.fits.\n' "${case_name}" >&2
        exit 1
    }
    cat "${case_directory}/results.txt"
    : > "${case_directory}/complete"
}

run_case unfiltered "none" --filter.lpfGaussFW 0
run_case sparse_response "sparse KLIP response" \
    --filter.lpfGaussFW 0 \
    --filter.psfResponse "${response_manifest}"
if [[ -n "${exact_response_manifest}" ]]; then
    run_case exact_response "exact KLIP response" \
        --filter.lpfGaussFW 0 \
        --filter.psfResponse "${exact_response_manifest}"
fi

IFS=, read -r -a gaussian_widths <<< "${gaussian_fwhms}"
if ((${#gaussian_widths[@]} == 0)); then
    printf 'GAUSSIAN_FWHMS must contain at least one positive FWHM.\n' >&2
    exit 2
fi
for gaussian_fwhm in "${gaussian_widths[@]}"; do
    gaussian_tag=${gaussian_fwhm//./p}
    gaussian_tag=${gaussian_tag//-/_minus_}
    run_case "gaussian_fwhm_${gaussian_tag}" "Gaussian FWHM ${gaussian_fwhm} pixels" \
        --filter.lpfGaussFW "${gaussian_fwhm}"
done

if [[ "${dry_run}" == true ]]; then
    printf '\nDry run complete; no hciAnalyze calls were started.\n'
    exit 0
fi

python3 - "${experiment_dir}" <<'PY'
import csv
import pathlib
import sys

experiment = pathlib.Path(sys.argv[1])
rows = []
for result_path in sorted(experiment.glob("*/results.txt")):
    lines = [line.split() for line in result_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not lines or lines[0] != ["plane", "mode", "signal", "separation", "position_angle", "contrast", "snr"]:
        raise SystemExit(f"unexpected hciAnalyze result format: {result_path}")
    for values in lines[1:]:
        if len(values) != len(lines[0]):
            raise SystemExit(f"malformed hciAnalyze result row: {result_path}")
        row = dict(zip(lines[0], values))
        row["case"] = result_path.parent.name
        row["snr"] = float(row["snr"])
        rows.append(row)

if not rows:
    raise SystemExit("no hciAnalyze results were found")

baseline = {
    (row["plane"], row["mode"], row["signal"]): row["snr"]
    for row in rows
    if row["case"] == "unfiltered"
}
for row in rows:
    key = (row["plane"], row["mode"], row["signal"])
    row["ratio_to_unfiltered"] = row["snr"] / baseline[key]

rows.sort(key=lambda row: (int(row["plane"]), row["case"], int(row["signal"])))
csv_path = experiment / "hciAnalyze_filter_summary.csv"
with csv_path.open("w", encoding="utf-8", newline="") as stream:
    fieldnames = ["case", "plane", "mode", "signal", "snr", "ratio_to_unfiltered"]
    writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)

markdown_path = experiment / "hciAnalyze_filter_summary.md"
with markdown_path.open("w", encoding="utf-8") as stream:
    stream.write("# KLIP hciAnalyze filter comparison\n\n")
    stream.write("| Case | Plane | KL modes | Signal | S/N | Ratio to unfiltered |\n")
    stream.write("|---|---:|---:|---:|---:|---:|\n")
    for row in rows:
        stream.write(
            f"| {row['case']} | {row['plane']} | {row['mode']} | {row['signal']} | "
            f"{row['snr']:.8g} | {row['ratio_to_unfiltered']:.8g} |\n"
        )

print(f"Wrote {markdown_path}")
PY

printf '\nKLIP hciAnalyze filter comparison complete: %s\n' \
    "${experiment_dir}/hciAnalyze_filter_summary.md"
