#!/usr/bin/env bash

set -euo pipefail

export LC_ALL=C
export PYTHONDONTWRITEBYTECODE=1

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "${script_dir}/../../.." && pwd)
roc_working_dir="${repo_root}/working/roc"
reference_experiment=${REFERENCE_EXPERIMENT:-"${roc_working_dir}/klip_matched_response_20260913T012040Z"}
finite_response_experiment=${FINITE_RESPONSE_EXPERIMENT:-"${roc_working_dir}/klip_finite_response_20260913T150323Z"}
base_config=${BASE_CONFIG:-"${script_dir}/klipReduce_afLepNaco_psf_response.conf"}
psf_file=${PSF_FILE:-/home/jrmales/Source/mxWork/NACO/AFLep/2011-10-21/out/psf_reg_median.fits}
klipreduce_bin=${KLIPREDUCE_BIN:-klipReduce}
experiment_dir=${EXPERIMENT_DIR:-"${roc_working_dir}/klip_central_response_$(date -u +%Y%m%dT%H%M%SZ)"}
reference_case_name=${REFERENCE_CASE:-radial_ld_fixed16_filter}
mode_counts=${MODE_COUNTS:-}
planet_sep=${PLANET_SEP:-11.73795339222688}
planet_pa=${PLANET_PA:-262.1667995998323}
planet_contrast=${PLANET_CONTRAST:-0.004763925929356391}
candidate_avoid_radius=${CANDIDATE_AVOID_RADIUS:-5}
amplitude_fractions=${AMPLITUDE_FRACTIONS:-0.25,0.5,1.0}
default_sample_specs="candidate:11.73795339222688:262.1667995998323:candidate;r07p8_pa000:7.8:0:clear;r07p8_pa090:7.8:90:clear;r07p8_pa180:7.8:180:clear;r11p4_pa000:11.4:0:clear;r11p4_pa090:11.4:90:clear;r11p4_pa180:11.4:180:clear;r33p0_pa000:33:0:clear;r33p0_pa090:33:90:clear;r33p0_pa180:33:180:clear;r33p0_pa270:33:270:clear;r58p2_pa000:58.2:0:clear;r58p2_pa090:58.2:90:clear;r58p2_pa180:58.2:180:clear;r58p2_pa270:58.2:270:clear"
sample_specs=${SAMPLE_SPECS:-"${default_sample_specs}"}
dry_run=false
analyze_only=false

usage()
{
    cat <<EOF
Usage: $(basename "$0") [--dry-run] [--analyze-only]

Measure paired central-difference KLIP responses at a bounded set of sky
locations and perturbation amplitudes. Each response is

    [KLIP(data + epsilon * PSF) - KLIP(data - epsilon * PSF)] / (2 epsilon).

The default clear samples cover the inner and outer region boundaries, the
planet's radial neighborhood, and a mid-radius control while remaining at
least CANDIDATE_AVOID_RADIUS pixels from the known planet. The candidate itself
is included as a separately labelled diagnostic.

Environment overrides:
  REFERENCE_EXPERIMENT       completed frozen-response experiment
                             (default: ${reference_experiment})
  FINITE_RESPONSE_EXPERIMENT completed one-sided finite-response experiment
                             (default: ${finite_response_experiment})
  REFERENCE_CASE             response case beneath REFERENCE_EXPERIMENT
                             (default: ${reference_case_name})
  EXPERIMENT_DIR             output directory (default: ${experiment_dir})
  BASE_CONFIG                maintained AF Lep/NACO KLIP configuration
  KLIPREDUCE_BIN             klipReduce executable
  PSF_FILE                   centered input PSF used for perturbations
  MODE_COUNTS                comma-separated KL modes; default is reference NMODES
  PLANET_SEP, PLANET_PA      known candidate coordinates
  PLANET_CONTRAST            fiducial contrast used to scale epsilon
  CANDIDATE_AVOID_RADIUS     minimum clear-sample distance in pixels
                             (default: ${candidate_avoid_radius})
  AMPLITUDE_FRACTIONS        positive epsilon/PLANET_CONTRAST values
                             (default: ${amplitude_fractions})
  SAMPLE_SPECS               semicolon-separated label:sep:PA:role entries;
                             role is candidate or clear
  OMP_NUM_THREADS            OpenMP worker limit passed through to klipReduce

The defaults launch 90 short KLIP reductions: 15 positions, three amplitudes,
and two signs. On the 2026-09-13 ROC timing this should take about 4.5 minutes.

Example:
  OMP_NUM_THREADS=48 nohup $(basename "$0") > klip_central_response_driver.log 2>&1 &
EOF
}

while (($#)); do
    case "$1" in
        --dry-run)
            dry_run=true
            shift
            ;;
        --analyze-only)
            analyze_only=true
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
    printf '\n[%s]\n' "${stage_directory#"${experiment_dir}/"}"
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
        printf '[%s] failed with status %d.\n' "${stage_directory#"${experiment_dir}/"}" \
            "${command_status}" >&2
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
        and int(header.get("KLIP PSF FILTER", 0)) == 1
        and bool(modes)
        and all((case / "finim_outputs" / f"klipPSF_mode{index:03d}_radial_response.fits").is_file()
                for index in range(len(modes)))
        and np.asarray(fits.getdata(case / "finim.fits")).size > 0
    )
except Exception:
    complete = False
raise SystemExit(0 if complete else 1)
PY
}

perturbation_complete()
{
    local final_image=$1
    local sample_sep=$2
    local sample_pa=$3
    local signed_contrast=$4
    python3 - "${final_image}" "${mode_counts}" "${sample_sep}" "${sample_pa}" "${signed_contrast}" \
        "${planet_sep}" "${planet_pa}" "${planet_contrast}" "${psf_file}" <<'PY'
import math
import pathlib
import sys

import numpy as np
from astropy.io import fits

path, modes_text, sep_text, pa_text, contrast_text, planet_sep_text, planet_pa_text, planet_contrast_text, psf_text = sys.argv[1:]

def vector(header, keyword):
    return [float(token) for token in str(header[keyword]).split(",")]

try:
    header = fits.getheader(path)
    data = np.asarray(fits.getdata(path))
    modes = [int(token) for token in str(header["NMODES"]).split(",")]
    fake_sep = vector(header, "FAKESEP")
    fake_pa = vector(header, "FAKEPA")
    fake_contrast = vector(header, "FAKECONT")
    known_sep = vector(header, "PLANETSEP")
    known_pa = vector(header, "PLANETPA")
    known_contrast = vector(header, "PLANETCONT")
    complete = (
        modes == [int(token) for token in modes_text.split(",")]
        and data.ndim in (2, 3)
        and data.size > 0
        and np.any(np.isfinite(data))
        and len(fake_sep) == len(fake_pa) == len(fake_contrast) == 1
        and len(known_sep) == len(known_pa) == len(known_contrast) == 1
        and math.isclose(fake_sep[0], float(sep_text), rel_tol=0, abs_tol=5e-4)
        and math.isclose(fake_pa[0], float(pa_text) % 360, rel_tol=0, abs_tol=5e-4)
        and math.isclose(fake_contrast[0], float(contrast_text), rel_tol=0, abs_tol=5e-9)
        and math.isclose(known_sep[0], float(planet_sep_text), rel_tol=0, abs_tol=5e-4)
        and math.isclose(known_pa[0], float(planet_pa_text) % 360, rel_tol=0, abs_tol=5e-4)
        and math.isclose(known_contrast[0], float(planet_contrast_text), rel_tol=0, abs_tol=5e-9)
        and pathlib.Path(str(header["FAKEFILE"]).strip()).resolve() == pathlib.Path(psf_text).resolve()
    )
except Exception:
    complete = False
raise SystemExit(0 if complete else 1)
PY
}

comparison_complete()
{
    local directory=$1
    python3 - "${directory}" <<'PY'
import json
import pathlib
import sys

directory = pathlib.Path(sys.argv[1])
required = (
    directory / "klip_central_response.json",
    directory / "klip_central_response.csv",
    directory / "klip_central_response.md",
    directory / "central_response.fits",
    directory / "frozen_response.fits",
    directory / "central_minus_scaled_frozen.fits",
)
try:
    summary = json.loads(required[0].read_text(encoding="utf-8"))
    complete = all(path.is_file() for path in required) and int(summary.get("schema", 0)) == 1
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
for required_option in --klip.Nmodes --fake.sep --fake.PA --fake.contrast --planet.sep; do
    if [[ "${help_text}" != *"${required_option}"* ]]; then
        printf 'klipReduce does not expose required option %s; build this checkout first.\n' \
            "${required_option}" >&2
        exit 1
    fi
done

reference_experiment=$(absolute_directory "${reference_experiment}")
finite_response_experiment=$(absolute_directory "${finite_response_experiment}")
reference_case="${reference_experiment}/${reference_case_name}"
reference_complete "${reference_case}" || {
    printf 'Reference KLIP response case is incomplete: %s\n' "${reference_case}" >&2
    exit 1
}
[[ -r "${finite_response_experiment}/empirical_response.fits" ]] || {
    printf 'Finite-response empirical cube is missing: %s\n' "${finite_response_experiment}" >&2
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

mkdir -p "${experiment_dir}"
experiment_dir=$(cd -- "${experiment_dir}" && pwd)
settings_path="${experiment_dir}/settings.env"
settings_lines=(
    "reference_experiment=${reference_experiment}"
    "finite_response_experiment=${finite_response_experiment}"
    "reference_case=${reference_case_name}"
    "base_config=$(cd -- "$(dirname -- "${base_config}")" && pwd)/$(basename -- "${base_config}")"
    "psf_file=$(cd -- "$(dirname -- "${psf_file}")" && pwd)/$(basename -- "${psf_file}")"
    "mode_counts=${mode_counts}"
    "planet_separation=${planet_sep}"
    "planet_position_angle=${planet_pa}"
    "planet_contrast=${planet_contrast}"
    "candidate_avoid_radius=${candidate_avoid_radius}"
    "amplitude_fractions=${amplitude_fractions}"
    "sample_specs=${sample_specs}"
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

sample_manifest="${experiment_dir}/samples.json"
python3 - "${sample_manifest}" "${sample_specs}" "${amplitude_fractions}" "${planet_sep}" "${planet_pa}" \
    "${planet_contrast}" "${candidate_avoid_radius}" <<'PY'
import json
import math
import pathlib
import re
import sys

path = pathlib.Path(sys.argv[1])
specs_text, fractions_text = sys.argv[2:4]
planet_sep, planet_pa, planet_contrast, avoid_radius = map(float, sys.argv[4:])
if not all(math.isfinite(value) for value in (planet_sep, planet_pa, planet_contrast, avoid_radius)):
    raise SystemExit("planet and avoidance parameters must be finite")
if planet_sep < 0 or planet_contrast <= 0 or avoid_radius < 0:
    raise SystemExit("planet separation/contrast and avoidance radius are invalid")

fractions = [float(token) for token in fractions_text.split(",")]
if not fractions or not all(math.isfinite(value) and value > 0 for value in fractions):
    raise SystemExit("AMPLITUDE_FRACTIONS must contain positive finite values")
if fractions != sorted(set(fractions)):
    raise SystemExit("AMPLITUDE_FRACTIONS must be unique and strictly increasing")

samples = []
for index, spec in enumerate(specs_text.split(";")):
    fields = spec.split(":")
    if len(fields) != 4:
        raise SystemExit(f"invalid SAMPLE_SPECS entry: {spec}")
    label, separation_text, pa_text, role = fields
    if re.fullmatch(r"[A-Za-z0-9_]+", label) is None or role not in ("candidate", "clear"):
        raise SystemExit(f"invalid sample label or role: {spec}")
    separation = float(separation_text)
    pa = float(pa_text) % 360
    if not math.isfinite(separation) or separation < 0 or not math.isfinite(pa):
        raise SystemExit(f"invalid sample coordinates: {spec}")
    angle = math.radians(pa - planet_pa)
    candidate_distance = math.sqrt(
        separation * separation + planet_sep * planet_sep
        - 2 * separation * planet_sep * math.cos(angle)
    )
    if role == "clear" and candidate_distance < avoid_radius:
        raise SystemExit(
            f"clear sample {label} is only {candidate_distance:.6g} pixels from the candidate"
        )
    if role == "candidate" and (
        not math.isclose(separation, planet_sep, rel_tol=0, abs_tol=1e-12)
        or not math.isclose(pa, planet_pa % 360, rel_tol=0, abs_tol=1e-12)
    ):
        raise SystemExit("the candidate sample must match PLANET_SEP and PLANET_PA")
    samples.append(
        {
            "index": index,
            "label": label,
            "separation": separation,
            "position_angle": pa,
            "role": role,
            "candidate_distance": candidate_distance,
        }
    )
if not samples or len({sample["label"] for sample in samples}) != len(samples):
    raise SystemExit("SAMPLE_SPECS must contain unique labelled samples")
if sum(sample["role"] == "candidate" for sample in samples) > 1:
    raise SystemExit("SAMPLE_SPECS may contain at most one candidate sample")

manifest = {
    "schema": 1,
    "planet": {
        "separation": planet_sep,
        "position_angle": planet_pa % 360,
        "contrast": planet_contrast,
    },
    "candidate_avoid_radius": avoid_radius,
    "amplitude_fractions": fractions,
    "samples": samples,
}
if path.exists():
    if json.loads(path.read_text(encoding="utf-8")) != manifest:
        raise SystemExit(f"existing sample manifest does not match this invocation: {path}")
else:
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY

provenance_path="${experiment_dir}/provenance.txt"
if [[ ! -e "${provenance_path}" ]]; then
    {
        printf 'created_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        printf 'hostname=%s\n' "$(hostname)"
        printf 'hciReduce_commit=%s\n' "$(git -C "${repo_root}" rev-parse HEAD)"
        printf 'klipreduce_path=%s\n' "$(command -v "${klipreduce_bin}")"
        sha256sum "$(command -v "${klipreduce_bin}")" "${base_config}" "${psf_file}" \
            "${reference_case}/finim.fits" \
            "${reference_case}/finim_outputs/klipPSF_manifest.fits" \
            "${finite_response_experiment}/empirical_response.fits"
    } > "${provenance_path}"
fi

if [[ "${analyze_only}" == false ]]; then
    while IFS=$'\t' read -r sample_index sample_label sample_sep sample_pa sample_role fraction fraction_tag epsilon; do
        for sign in plus minus; do
            signed_contrast=${epsilon}
            if [[ "${sign}" == minus ]]; then
                signed_contrast="-${epsilon}"
            fi
            stage_directory="${experiment_dir}/runs/${sample_label}/fraction_${fraction_tag}/${sign}"
            final_image="${stage_directory}/finim.fits"
            if [[ -f "${final_image}" ]] && \
                perturbation_complete "${final_image}" "${sample_sep}" "${sample_pa}" "${signed_contrast}"; then
                printf '\n[%s] completed reduction exists; skipping.\n' \
                    "${stage_directory#"${experiment_dir}/"}"
                continue
            fi
            if [[ -e "${stage_directory}/run.log" || -e "${final_image}" ]]; then
                printf 'Incomplete perturbation stage exists; refusing to overwrite: %s\n' \
                    "${stage_directory}" >&2
                exit 1
            fi
            command_line=(
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
                --fake.sep "${sample_sep}"
                --fake.PA "${sample_pa}"
                --fake.contrast "${signed_contrast}"
                --fake.subtractPlanet=false
                --output.directory "${stage_directory}"
                --output.fileName finim.fits
                --output.exactFName=true
                --showTiming=true
            )
            run_timed "${stage_directory}" "${command_line[@]}"
            if [[ "${dry_run}" == false ]]; then
                perturbation_complete "${final_image}" "${sample_sep}" "${sample_pa}" \
                    "${signed_contrast}" || {
                    printf 'Perturbation reduction did not publish the expected final cube.\n' >&2
                    exit 1
                }
                : > "${stage_directory}/complete"
            fi
        done
    done < <(python3 - "${sample_manifest}" <<'PY'
import json
import pathlib
import sys

manifest = json.loads(pathlib.Path(sys.argv[1]).read_text(encoding="utf-8"))
contrast = float(manifest["planet"]["contrast"])
for sample in manifest["samples"]:
    for fraction in manifest["amplitude_fractions"]:
        tag = format(float(fraction), ".12g").replace(".", "p")
        print(
            sample["index"], sample["label"], sample["separation"], sample["position_angle"],
            sample["role"], fraction, tag, contrast * float(fraction), sep="\t"
        )
PY
)
fi

if [[ "${dry_run}" == true ]]; then
    printf '\nDry run complete; central-response comparison requires the perturbation reductions.\n'
    exit 0
fi

if comparison_complete "${experiment_dir}"; then
    printf '\n[comparison] completed analysis exists; skipping.\n'
else
    for product in klip_central_response.json klip_central_response.csv klip_central_response.md \
        central_response.fits frozen_response.fits central_minus_scaled_frozen.fits; do
        [[ ! -e "${experiment_dir}/${product}" ]] || {
            printf 'Incomplete comparison products exist; refusing to overwrite: %s\n' "${experiment_dir}" >&2
            exit 1
        }
    done
    comparison_command=(
        python3 "${script_dir}/compare_klip_central_response.py"
        "${reference_experiment}"
        "${experiment_dir}"
        --reference-case "${reference_case_name}"
        --finite-response-experiment "${finite_response_experiment}"
    )
    shell_join "${comparison_command[@]}" > "${experiment_dir}/comparison_command.txt"
    "${comparison_command[@]}" 2>&1 | tee "${experiment_dir}/comparison.log"
    comparison_complete "${experiment_dir}" || {
        printf 'KLIP central-response comparison did not publish a complete result.\n' >&2
        exit 1
    }
fi

printf '\nKLIP central-response validation complete: %s\n' \
    "${experiment_dir}/klip_central_response.md"
