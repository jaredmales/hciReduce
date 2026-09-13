#!/usr/bin/env bash

set -euo pipefail

export LC_ALL=C
export PYTHONDONTWRITEBYTECODE=1

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "${script_dir}/../../.." && pwd)
roc_working_dir="${repo_root}/working/roc"
reference_experiment=${REFERENCE_EXPERIMENT:-"${roc_working_dir}/klip_matched_response_20260913T012040Z"}
base_config=${BASE_CONFIG:-"${script_dir}/klipReduce_afLepNaco_psf_response.conf"}
psf_file=${PSF_FILE:-/home/jrmales/Source/mxWork/NACO/AFLep/2011-10-21/out/psf_reg_median.fits}
klipreduce_bin=${KLIPREDUCE_BIN:-klipReduce}
experiment_dir=${EXPERIMENT_DIR:-"${roc_working_dir}/klip_finite_response_$(date -u +%Y%m%dT%H%M%SZ)"}
mode_counts=${MODE_COUNTS:-}
planet_sep=${PLANET_SEP:-11.73795339222688}
planet_pa=${PLANET_PA:-262.1667995998323}
planet_contrast=${PLANET_CONTRAST:-0.004763925929356391}
reference_case_name=${REFERENCE_CASE:-radial_ld_fixed16_filter}
dry_run=false

usage()
{
    cat <<EOF
Usage: $(basename "$0") [--dry-run]

Test the KLIP basis-derivative hypothesis with one complete negative-companion
reduction. The driver reuses the frozen sparse response and original science
cube from REFERENCE_EXPERIMENT, subtracts the accepted P4 exact-negative planet
before refitting KLIP, and compares the end-to-end removed signal with the
frozen-basis response for every retained KL mode.

Environment overrides:
  REFERENCE_EXPERIMENT  completed KLIP matched-response experiment
                        (default: ${reference_experiment})
  REFERENCE_CASE        response case beneath the reference experiment
                        (default: ${reference_case_name})
  EXPERIMENT_DIR        output directory (default: ${experiment_dir})
  BASE_CONFIG           maintained AF Lep/NACO KLIP configuration
                        (default: ${base_config})
  KLIPREDUCE_BIN        klipReduce executable (default: ${klipreduce_bin})
  PSF_FILE              centered input PSF used for subtraction
                        (default: ${psf_file})
  MODE_COUNTS           comma-separated KL modes; default is reference NMODES
  PLANET_SEP            subtracted separation in pixels (default: ${planet_sep})
  PLANET_PA             subtracted PA east of north (default: ${planet_pa})
  PLANET_CONTRAST       positive subtracted contrast (default: ${planet_contrast})
  OMP_NUM_THREADS       OpenMP worker limit passed through to klipReduce

Example:
  OMP_NUM_THREADS=48 nohup $(basename "$0") > klip_finite_response_driver.log 2>&1 &
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

reference_complete()
{
    local case_directory=$1
    python3 - "${case_directory}" <<'PY'
import pathlib
import sys

import numpy as np
from astropy.io import fits

case = pathlib.Path(sys.argv[1])
manifest = case / "finim_outputs" / "klipPSF_manifest.fits"
required = (case / "finim.fits", case / "finim_filtered.fits", manifest)
try:
    header = fits.getheader(manifest)
    modes = [int(token) for token in str(header["NMODES"]).split(",")]
    complete = (
        all(path.is_file() for path in required)
        and int(header.get("KLIP PSF COMPLETE", 0)) == 1
        and int(header.get("KLIP PSF PRODUCT SCHEMA", 0)) == 1
        and int(header.get("KLIP PSF FILTER", 0)) == 1
        and bool(modes)
        and len(set(modes)) == len(modes)
        and all((case / "finim_outputs" / f"klipPSF_mode{index:03d}_radial_response.fits").is_file()
                for index in range(len(modes)))
        and np.asarray(fits.getdata(case / "finim.fits")).size > 0
    )
except Exception:
    complete = False
raise SystemExit(0 if complete else 1)
PY
}

negative_complete()
{
    local final_image=$1
    python3 - "${final_image}" "${mode_counts}" "${planet_sep}" "${planet_pa}" "${planet_contrast}" \
        "${psf_file}" <<'PY'
import math
import pathlib
import sys

import numpy as np
from astropy.io import fits

path, modes_text, sep_text, pa_text, contrast_text, psf_text = sys.argv[1:]

def vector(header, keyword):
    return [float(token) for token in str(header[keyword]).split(",")]

try:
    header = fits.getheader(path)
    data = np.asarray(fits.getdata(path))
    modes = [int(token) for token in str(header["NMODES"]).split(",")]
    requested_modes = [int(token) for token in modes_text.split(",")]
    separation = vector(header, "PLANETSEP")
    pa = vector(header, "PLANETPA")
    contrast = vector(header, "PLANETCONT")
    complete = (
        modes == requested_modes
        and data.ndim in (2, 3)
        and data.size > 0
        and np.any(np.isfinite(data))
        and len(separation) == len(pa) == len(contrast) == 1
        and math.isclose(separation[0], float(sep_text), rel_tol=0, abs_tol=5e-4)
        and math.isclose(pa[0], float(pa_text), rel_tol=0, abs_tol=5e-4)
        and math.isclose(contrast[0], float(contrast_text), rel_tol=0, abs_tol=5e-9)
        and pathlib.Path(str(header["FAKEFILE"]).strip()).resolve() == pathlib.Path(psf_text).resolve()
        and not str(header.get("FAKESEP", "")).strip()
    )
except Exception:
    complete = False
raise SystemExit(0 if complete else 1)
PY
}

comparison_complete()
{
    local directory=$1
    python3 - "${directory}" "${mode_counts}" "${planet_sep}" "${planet_pa}" "${planet_contrast}" <<'PY'
import json
import math
import pathlib
import sys

directory, modes_text, sep_text, pa_text, contrast_text = sys.argv[1:]
directory = pathlib.Path(directory)
required = (
    directory / "klip_finite_response_comparison.json",
    directory / "klip_finite_response_comparison.csv",
    directory / "klip_finite_response_comparison.md",
    directory / "empirical_response.fits",
    directory / "frozen_response.fits",
    directory / "empirical_minus_scaled_frozen.fits",
)
try:
    summary = json.loads(required[0].read_text(encoding="utf-8"))
    injection = summary["negative_injection"]
    modes = [int(row["mode_count"]) for row in summary["modes"]]
    complete = (
        all(path.is_file() for path in required)
        and modes == [int(token) for token in modes_text.split(",")]
        and math.isclose(float(injection["separation"]), float(sep_text), rel_tol=0, abs_tol=1e-12)
        and math.isclose(float(injection["position_angle"]), float(pa_text) % 360, rel_tol=0, abs_tol=1e-12)
        and math.isclose(float(injection["contrast"]), float(contrast_text), rel_tol=0, abs_tol=1e-15)
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
command -v "${klipreduce_bin}" >/dev/null 2>&1 || {
    printf 'klipReduce executable was not found: %s\n' "${klipreduce_bin}" >&2
    exit 1
}
help_text=$("${klipreduce_bin}" --help 2>&1)
for required_option in --klip.Nmodes --planet.sep --planet.PA --planet.contrast --fake.subtractPlanet; do
    if [[ "${help_text}" != *"${required_option}"* ]]; then
        printf 'klipReduce does not expose required option %s; build this checkout first.\n' \
            "${required_option}" >&2
        exit 1
    fi
done

reference_experiment=$(absolute_directory "${reference_experiment}")
reference_case="${reference_experiment}/${reference_case_name}"
reference_complete "${reference_case}" || {
    printf 'Reference KLIP response case is incomplete: %s\n' "${reference_case}" >&2
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
python3 - "${mode_counts}" "${planet_sep}" "${planet_pa}" "${planet_contrast}" <<'PY'
import math
import sys
try:
    modes = [int(token) for token in sys.argv[1].split(",")]
    values = [float(token) for token in sys.argv[2:]]
except ValueError:
    raise SystemExit("mode counts and planet parameters must be numeric")
if not modes or len(set(modes)) != len(modes) or any(mode <= 0 for mode in modes):
    raise SystemExit("MODE_COUNTS must contain unique positive integers")
if not all(math.isfinite(value) for value in values) or values[0] < 0 or values[2] <= 0:
    raise SystemExit("planet separation, PA, and positive contrast must be finite and valid")
PY

mkdir -p "${experiment_dir}"
settings_path="${experiment_dir}/settings.env"
settings_lines=(
    "reference_experiment=${reference_experiment}"
    "reference_case=${reference_case_name}"
    "base_config=$(cd -- "$(dirname -- "${base_config}")" && pwd)/$(basename -- "${base_config}")"
    "psf_file=$(cd -- "$(dirname -- "${psf_file}")" && pwd)/$(basename -- "${psf_file}")"
    "mode_counts=${mode_counts}"
    "planet_separation=${planet_sep}"
    "planet_position_angle=${planet_pa}"
    "planet_contrast=${planet_contrast}"
    "subtract_planet=true"
    "omp_num_threads=${OMP_NUM_THREADS:-}"
)
if [[ -e "${settings_path}" ]]; then
    mapfile -t saved_settings < "${settings_path}"
    [[ "${saved_settings[*]}" == "${settings_lines[*]}" ]] || {
        printf 'Existing experiment settings do not match this invocation: %s\n' "${settings_path}" >&2
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
        sha256sum "$(command -v "${klipreduce_bin}")" "${base_config}" "${psf_file}" \
            "${reference_case}/finim.fits" "${reference_case}/finim_outputs/klipPSF_manifest.fits"
    } > "${provenance_path}"
fi

negative_case="${experiment_dir}/negative_planet"
negative_final="${negative_case}/finim.fits"
if negative_complete "${negative_final}"; then
    printf '\n[negative_planet] completed reduction exists; skipping.\n'
elif [[ -e "${negative_case}/run.log" || -e "${negative_final}" ]]; then
    printf 'Incomplete negative-planet stage exists; refusing to overwrite: %s\n' "${negative_case}" >&2
    exit 1
else
    negative_command=(
        "${klipreduce_bin}"
        --config "${base_config}"
        --klip.Nmodes "${mode_counts}"
        --klip.psfFile ""
        --klip.outputPSFModels=false
        --klip.psfFilter=false
        --planet.sep "${planet_sep}"
        --planet.PA "${planet_pa}"
        --planet.contrast "${planet_contrast}"
        --fake.method single
        --fake.fileName "${psf_file}"
        --fake.subtractPlanet=true
        --output.directory "${negative_case}"
        --output.fileName finim.fits
        --output.exactFName=true
        --showTiming=true
    )
    run_timed "${negative_case}" "${negative_command[@]}"
    if [[ "${dry_run}" == false ]]; then
        negative_complete "${negative_final}" || {
            printf 'Negative-planet reduction did not publish the expected final cube.\n' >&2
            exit 1
        }
        : > "${negative_case}/complete"
    fi
fi

if [[ "${dry_run}" == true ]]; then
    printf '\nDry run complete; finite-response comparison requires the negative reduction.\n'
    exit 0
fi

if comparison_complete "${experiment_dir}"; then
    printf '\n[comparison] completed analysis exists; skipping.\n'
else
    for product in klip_finite_response_comparison.json klip_finite_response_comparison.csv \
        klip_finite_response_comparison.md empirical_response.fits frozen_response.fits \
        empirical_minus_scaled_frozen.fits; do
        [[ ! -e "${experiment_dir}/${product}" ]] || {
            printf 'Incomplete comparison products exist; refusing to overwrite: %s\n' "${experiment_dir}" >&2
            exit 1
        }
    done
    comparison_command=(
        python3 "${script_dir}/compare_klip_finite_response.py"
        "${reference_experiment}"
        "${negative_final}"
        "${experiment_dir}"
        --reference-case "${reference_case_name}"
        --separation "${planet_sep}"
        --pa "${planet_pa}"
        --contrast "${planet_contrast}"
    )
    shell_join "${comparison_command[@]}" > "${experiment_dir}/comparison_command.txt"
    "${comparison_command[@]}" 2>&1 | tee "${experiment_dir}/comparison.log"
    comparison_complete "${experiment_dir}" || {
        printf 'KLIP finite-response comparison did not publish a complete result.\n' >&2
        exit 1
    }
fi

printf '\nKLIP finite-response validation complete: %s\n' \
    "${experiment_dir}/klip_finite_response_comparison.md"
