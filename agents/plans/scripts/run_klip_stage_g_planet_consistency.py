#!/usr/bin/env python3
"""Run exact-contrast KLIP injections that match the Stage-F planet aperture."""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter
from scipy.stats import t as student_t

SCRIPT_PATH = Path(__file__).resolve()
SCRIPT_DIRECTORY = SCRIPT_PATH.parent
FROZEN_DEPENDENCIES = (SCRIPT_PATH.parents[2] / "software"
                       if len(SCRIPT_PATH.parents) > 2 else Path("/nonexistent"))
if (FROZEN_DEPENDENCIES / "run_klip_stage_c_development.py").is_file():
    sys.path.insert(0, str(FROZEN_DEPENDENCIES))
else:
    sys.path.insert(0, str(SCRIPT_DIRECTORY))

import prepare_klip_stage_c_development as preparation  # noqa: E402
import run_klip_covariance_stage_a as stage  # noqa: E402
import run_klip_stage_c_development as development  # noqa: E402
import run_klip_stage_e_validation as validation  # noqa: E402
import run_klip_stage_f_planet as stage_f  # noqa: E402


STAGE_DIRECTORY = "stage_g_planet_consistency"
SITE_COUNT = 12
NOMINAL_RADIUS = 12.0
PRIMARY_MODE = 200
ANALYZER_REPORTING_APERTURE_RADIUS = 60.0
PSF_BLUR_FWHM = {"nominal": 0.0, "blur0p9": 0.9, "blur1p8": 1.8}
ARM_DEFINITIONS = (
    {"name": "nominal_integer", "psf": "nominal", "phase": "integer"},
    {"name": "nominal_analysis_phase", "psf": "nominal", "phase": "analysis"},
    {"name": "nominal_optimized_phase", "psf": "nominal", "phase": "optimized"},
    {"name": "blur0p9_optimized_phase", "psf": "blur0p9", "phase": "optimized"},
    {"name": "blur1p8_optimized_phase", "psf": "blur1p8", "phase": "optimized"},
)
PRIMARY_ARM = "nominal_optimized_phase"
PAIR_DEFINITIONS = (
    ("gaussian_fwhm3p6_minus_exact_identity", "gaussian_fwhm3p6", "exact_identity"),
    ("gaussian_fwhm3p6_minus_sparse_identity", "gaussian_fwhm3p6", "sparse_identity"),
    ("gaussian_fwhm3p6_minus_gaussian_fwhm2p4", "gaussian_fwhm3p6", "gaussian_fwhm2p4"),
    ("raw_rectangular_m0p3_minus_exact_identity", "raw_rectangular_m0p3", "exact_identity"),
    ("radial_hann_m0p1_trunc0p75_minus_exact_identity",
     "radial_hann_m0p1_trunc0p75", "exact_identity"),
)


def read(path: Path) -> object:
    """Read one JSON document."""
    return json.loads(path.read_text(encoding="utf-8"))


def stage_g_root(root: Path) -> Path:
    """Return the isolated Stage-G product directory."""
    return root / STAGE_DIRECTORY


def source_from_polar(shape: tuple[int, int], separation: float,
                      position_angle: float) -> tuple[float, float]:
    """Return a floating row and column from hciReduce separation and PA."""
    metadata = {"known_planet": {"separation": separation,
                                  "position_angle": position_angle}}
    return tuple(float(value) for value in development.raw.planet_position(shape, metadata))


def polar_from_source(shape: tuple[int, int], row: float,
                      column: float) -> tuple[float, float]:
    """Return hciReduce separation and PA for a floating row and column."""
    center_row = 0.5 * (shape[1] - 1)
    center_column = 0.5 * (shape[0] - 1)
    delta_row = row - center_row
    delta_column = column - center_column
    separation = math.hypot(delta_row, delta_column)
    position_angle = math.degrees(math.atan2(-delta_row, delta_column)) % 360.0
    return separation, position_angle


def fractional_phase(source: tuple[float, float]) -> tuple[float, float]:
    """Return the source displacement from its nearest native pixel."""
    return source[0] - round(source[0]), source[1] - round(source[1])


def sample_summary(values: list[float]) -> dict[str, float | int]:
    """Return finite sample statistics with the sample standard deviation."""
    array = np.asarray(values, dtype=np.float64)
    stage.require(len(array) >= 2 and np.all(np.isfinite(array)),
                  "sample summary requires at least two finite values")
    return {
        "count": int(len(array)),
        "mean": float(np.mean(array)),
        "sample_standard_deviation": float(np.std(array, ddof=1)),
        "minimum": float(np.min(array)),
        "maximum": float(np.max(array)),
    }


def normal_predictive(summary: dict[str, float | int], observed: float) \
        -> dict[str, float]:
    """Return a diagnostic normal predictive comparison for one observation."""
    count = int(summary["count"])
    mean = float(summary["mean"])
    deviation = float(summary["sample_standard_deviation"])
    stage.require(count >= 2 and deviation > 0 and np.isfinite(observed),
                  "normal predictive comparison is undefined")
    standardized = (observed - mean) / deviation
    predictive_t = standardized / math.sqrt(1.0 + 1.0 / count)
    return {
        "observed": observed,
        "standardized_deviation": standardized,
        "predictive_t": predictive_t,
        "degrees_of_freedom": count - 1,
        "two_sided_p": float(2 * student_t.sf(abs(predictive_t), df=count - 1)),
    }


def make_psf_controls(source: Path, destination: Path) -> list[dict[str, object]]:
    """Create flux-preserving nominal and Gaussian-broadened injection PSFs."""
    data, header = fits.getdata(source, header=True)
    data = np.asarray(data, dtype=np.float64)
    stage.require(data.ndim == 2 and np.all(np.isfinite(data)) and np.sum(data) > 0,
                  "nominal injection PSF is invalid")
    destination.mkdir()
    original_sum = float(np.sum(data))
    records = []
    for name, fwhm in PSF_BLUR_FWHM.items():
        path = destination / f"{name}.fits"
        if fwhm == 0:
            shutil.copy2(source, path)
        else:
            sigma = fwhm / math.sqrt(8 * math.log(2))
            control = gaussian_filter(data, sigma=sigma, mode="constant", cval=0.0,
                                      truncate=4.0)
            stage.require(np.all(np.isfinite(control)) and np.sum(control) > 0,
                          f"invalid broadened PSF {name}")
            control *= original_sum / float(np.sum(control))
            control_header = header.copy()
            control_header["PSFCTRL"] = name
            control_header["PSFFWHM"] = fwhm
            control_header["PSFSUM0"] = original_sum
            control_header["PSFSUM1"] = float(np.sum(control))
            fits.writeto(path, control.astype(np.float32), control_header)
        stored = np.asarray(fits.getdata(path), dtype=np.float64)
        stage.require(np.isclose(float(np.sum(stored)), original_sum, rtol=2e-6, atol=0),
                      f"PSF control {name} did not preserve total flux")
        records.append({
            "name": name,
            "additional_gaussian_fwhm_pixels": fwhm,
            "sum": float(np.sum(stored)),
            "peak": float(np.max(stored)),
            "fingerprint": stage.fingerprint(path),
        })
    return records


def select_sites(protocol: dict[str, object]) -> list[dict[str, object]]:
    """Select twelve score-blind radius-12 sites unused for positive injections."""
    candidates = [dict(site) for site in protocol["sites"]
                  if float(site["nominal_radius"]) == NOMINAL_RADIUS and
                  site["role"] in {"calibration", "heldout_null"}]
    stage.require(len(candidates) >= SITE_COUNT,
                  "insufficient non-positive radius-12 sites")
    chosen = preparation.maximin(candidates, SITE_COUNT)
    stage.require(len(chosen) == SITE_COUNT and
                  len({(int(site["row"]), int(site["column"])) for site in chosen}) == SITE_COUNT,
                  "Stage-G site selection is incomplete")
    return chosen


def task_command(output: Path, parent: dict[str, object], task: dict[str, object],
                 psf_paths: dict[str, Path]) -> list[str]:
    """Construct one exact-contrast signal-free injection command."""
    source = task["source"]
    known = parent["known_planet"]
    paths = parent["paths"]
    return [
        str(paths["klipreduce"]), "--config", str(output / "reduction.conf"),
        "--input.directory=", "--input.fileList", str(output / "inputs.txt"),
        "--klip.Nmodes", ",".join(map(str, stage.MODES)),
        "--psfResponse.file=", "--psfResponse.outputModels=false",
        "--psfResponse.filter=false",
        "--planet.sep", str(known["separation"]),
        "--planet.PA", str(known["position_angle"]),
        "--planet.contrast", str(known["contrast"]),
        "--fake.method", "single",
        "--fake.fileName", str(psf_paths[str(task["psf_control"])]),
        "--fake.sep", str(source["separation"]),
        "--fake.PA", str(source["position_angle"]),
        "--fake.contrast", str(task["contrast"]),
        "--fake.subtractPlanet=true",
        "--output.directory", str(output / "reductions" / str(task["name"])),
        "--output.fileName", "finim.fits", "--output.exactFName=true",
        "--showTiming=true",
    ]


def verify_stage_f_boundary(root: Path) \
        -> tuple[dict[str, object], Path, dict[str, object]]:
    """Verify Stage E, the frozen Stage-F manifest, and completed planet products."""
    protocol, _, analyzer = stage_f.verify_stage_e(root)
    stage_f.verify_complete(root)
    output = stage_f.stage_f_root(root)
    manifest = read(output / "manifest.json")
    stage.verify(manifest["input_records"] + manifest["software_records"] +
                 manifest["model_receipts"])
    stage.require(stage_f.model_receipts(root) == manifest["model_receipts"],
                  "Stage-F calibration model products changed")
    receipt = read(output / "complete.json")
    stage.require(receipt["known_planet_opened"] and
                  not receipt["stage_e_policy_changed"] and
                  not receipt["method_selection_performed"],
                  "Stage-F completion boundary changed")
    return protocol, analyzer, manifest


def preflight_geometries(parent: dict[str, object],
                         tasks: list[dict[str, object]]) -> list[dict[str, object]]:
    """Require exact response and baseline-fit support for every unique aperture."""
    baseline = np.asarray(fits.getdata(parent["paths"]["baseline"], memmap=True),
                          dtype=np.float64)
    parent_response = Path(str(parent["parent_response"]))
    paths = development.response47.product_paths(parent_response)
    coordinates = np.asarray(fits.getdata(paths["coordinates"], memmap=True),
                             dtype=np.float64).T
    lookup = {(int(row), int(column)): index
              for index, (row, column, _, _) in enumerate(coordinates)}
    mode_index = stage.MODES.index(PRIMARY_MODE)
    responses = np.asarray(fits.getdata(paths["responses"][mode_index], memmap=True),
                           dtype=np.float64)
    validities = np.asarray(fits.getdata(paths["validities"][mode_index], memmap=True),
                            dtype=np.float64) > 0.5
    optimized = source_from_polar((128, 128), float(parent["known_planet"]["separation"]),
                                  float(parent["known_planet"]["position_angle"]))
    unique: dict[tuple[float, float], dict[str, object]] = {}
    for task in tasks:
        source = task["source"]
        unique[(float(source["row"]), float(source["column"]))] = task
    result = []
    for index, task in enumerate(unique.values(), start=1):
        settings = dict(task["settings"])
        source, aperture, searches, bins = stage_f.analysis_geometry((128, 128), settings)
        expected = task["source"]
        stage.require(np.allclose(source, [expected["row"], expected["column"]],
                                  rtol=0, atol=2e-12),
                      "Stage-G source geometry does not round trip")
        stage.require(all(query in lookup for query in searches),
                      "Stage-G aperture leaves exact-response coordinates")
        minimum_support = development.SUPPORT * development.SUPPORT
        width = int(parent["training"]["candidate_specific_half_width_by_radius"][
            str(stage_f.policy_radius(float(expected["separation"])))])
        for query in searches:
            source_index = lookup[query]
            template = development.raw.crop(responses[source_index].T, development.SUPPORT)
            validity = development.raw.crop(validities[source_index].T, development.SUPPORT)
            fitted = stage_f.fit_planet_query(
                baseline[mode_index], query, searches, template, validity, optimized, width)
            minimum_support = min(minimum_support,
                                  int(np.count_nonzero(fitted["support"])))
        result.append({
            "source_row_column": list(source),
            "separation": float(expected["separation"]),
            "position_angle": float(expected["position_angle"]),
            "aperture_pixels": int(np.count_nonzero(aperture)),
            "queries": len(searches),
            "radial_bins": bins,
            "minimum_exact_support_pixels": minimum_support,
            "training_half_width": width,
        })
        print(f"preflight geometry {index}/{len(unique)}: "
              f"{task['base_site']} {task['phase']}", flush=True)
    return result


def prepare(root: Path, config: Path) -> None:
    """Freeze the exact-contrast aperture and PSF-control campaign."""
    root = root.resolve()
    config = config.resolve()
    parent, analyzer, _ = verify_stage_f_boundary(root)
    stage.require(sorted(os.sched_getaffinity(0)) == parent["resources"]["cpu_affinity"],
                  "Stage-G prepare must use the frozen response CPU affinity")
    output = stage_g_root(root)
    stage.require(not output.exists(), f"Stage-G output already exists: {output}")
    settings = stage_f.parse_analysis_config(config)
    stage.require(settings == stage_f.EXPECTED_SETTINGS,
                  "Stage-G requires the reviewed working/analyze.conf settings")
    output.mkdir()
    (output / "software").mkdir()
    shutil.copy2(config, output / "analyze.conf")
    shutil.copy2(root / "reduction.conf", output / "reduction.conf")
    shutil.copy2(root / "inputs.txt", output / "inputs.txt")

    frozen_runner = output / "software" / SCRIPT_PATH.name
    shutil.copy2(SCRIPT_PATH, frozen_runner)
    stage_f_source = stage_f.stage_f_root(root) / "software" / "run_klip_stage_f_planet.py"
    stage.require(stage_f_source.is_file(), "frozen Stage-F runner is missing")
    shutil.copy2(stage_f_source, output / "software" / stage_f_source.name)

    nominal_psf = Path(str(parent["paths"]["psf"]))
    psf_records = make_psf_controls(nominal_psf, output / "psfs")
    psf_paths = {str(record["name"]): Path(str(record["fingerprint"]["path"]))
                 for record in psf_records}

    analysis_source = source_from_polar((128, 128), settings["separation"],
                                        settings["position_angle"])
    optimized_source = source_from_polar(
        (128, 128), float(parent["known_planet"]["separation"]),
        float(parent["known_planet"]["position_angle"]))
    phases = {
        "integer": (0.0, 0.0),
        "analysis": fractional_phase(analysis_source),
        "optimized": fractional_phase(optimized_source),
    }
    sites = select_sites(parent)
    tasks = []
    for site_index, site in enumerate(sites):
        for arm in ARM_DEFINITIONS:
            phase = phases[str(arm["phase"])]
            row = float(site["row"]) + phase[0]
            column = float(site["column"]) + phase[1]
            separation, position_angle = polar_from_source((128, 128), row, column)
            source = source_from_polar((128, 128), separation, position_angle)
            stage.require(np.allclose(source, [row, column], rtol=0, atol=2e-12),
                          "Stage-G source polar conversion is not reversible")
            arm_settings = dict(settings)
            arm_settings.update(separation=separation, position_angle=position_angle)
            tasks.append({
                "name": f"site{site_index:02d}__{arm['name']}",
                "site_index": site_index,
                "base_site": site["name"],
                "base_row_column": [int(site["row"]), int(site["column"])],
                "arm": arm["name"],
                "psf_control": arm["psf"],
                "phase": arm["phase"],
                "phase_offset_row_column": list(phase),
                "source": {"row": source[0], "column": source[1],
                           "separation": separation, "position_angle": position_angle},
                "contrast": stage.PLANET_CONTRAST,
                "settings": arm_settings,
            })
    stage.require(len(tasks) == SITE_COUNT * len(ARM_DEFINITIONS),
                  "Stage-G task count changed")
    for task in tasks:
        task["command"] = task_command(output, parent, task, psf_paths)
    geometry = preflight_geometries(parent, tasks)

    protocol = {
        "schema": 1,
        "stage": "KLIP Stage G planet-consistency injections",
        "purpose": "direct optimized-contrast radius-12 injections using the Stage-F aperture, subpixel phases, and predeclared PSF broadening controls",
        "parent_root": str(root),
        "stage_f_root": str(stage_f.stage_f_root(root)),
        "primary_mode": PRIMARY_MODE,
        "nominal_radius": NOMINAL_RADIUS,
        "planet_contrast": stage.PLANET_CONTRAST,
        "site_count": SITE_COUNT,
        "site_selection": "deterministic angular maximin from radius-12 calibration and held-out-null centers; no positive score used",
        "sites": sites,
        "phase_offsets_row_column": {name: list(value) for name, value in phases.items()},
        "arms": list(ARM_DEFINITIONS),
        "primary_arm": PRIMARY_ARM,
        "psf_controls": psf_records,
        "tasks": tasks,
        "methods": list(stage_f.PLANET_METHODS),
        "primary_endpoint": "paired mode-200 aperture-maximum and nearest-pixel SNR differences: covariance minus exact identity",
        "secondary_endpoints": [
            "Gaussian 3.6 minus exact identity and Gaussian 2.4",
            "integer versus configured and optimized subpixel phase",
            "0.25 and 0.5 lambda/D additional Gaussian PSF broadening at optimized phase",
            "finite positive-response fidelity to the nominal exact response",
        ],
        "analysis_contract": {
            "aperture_radius_pixels": settings["aperture_radius"],
            "hcianalyze_reporting_aperture_radius_pixels":
                ANALYZER_REPORTING_APERTURE_RADIUS,
            "hcianalyze_reporting_aperture_affects_annular_normalization": False,
            "source_exclusion_radius_pixels": settings["source_radius"],
            "known_planet_and_trial_excluded_from_annular_profile": True,
            "weights_fit_from_signal_free_baseline": True,
            "entire_injection_aperture_excluded_from_training": True,
            "no_parameter_selected_from_planet_snr": True,
        },
        "preflight_geometries": geometry,
        "resources": parent["resources"],
        "paths": {
            "klipreduce": parent["paths"]["klipreduce"],
            "hcianalyze": str(analyzer),
            "baseline": parent["paths"]["baseline"],
            "parent_response": parent["parent_response"],
            "sparse_manifest": parent["paths"]["sparse_manifest"],
        },
    }
    stage.write_json(output / "protocol.json", protocol)
    model_receipts = stage_f.model_receipts(root)
    input_records = [
        stage.fingerprint(output / "analyze.conf"),
        stage.fingerprint(output / "reduction.conf"),
        stage.fingerprint(output / "inputs.txt"),
        stage.fingerprint(root / "protocol.json"),
        stage.fingerprint(root / "validation_complete.json"),
        stage.fingerprint(stage_f.stage_f_root(root) / "complete.json"),
        stage.fingerprint(stage_f.stage_f_root(root) / "manifest.json"),
        stage.fingerprint(stage_f.stage_f_root(root) / "planet" / "results.json"),
        stage.fingerprint(output / "protocol.json"),
    ] + [record["fingerprint"] for record in psf_records]
    software_records = [
        stage.fingerprint(frozen_runner),
        stage.fingerprint(output / "software" / stage_f_source.name),
        stage.fingerprint(Path(str(parent["paths"]["klipreduce"]))),
        stage.fingerprint(analyzer),
    ]
    manifest = {
        "schema": 1,
        "purpose": "immutable Stage-G exact-contrast planet-consistency campaign",
        "root": str(root),
        "input_records": input_records,
        "software_records": software_records,
        "model_receipts": model_receipts,
        "stage_f_manifest_sha256": stage.fingerprint(
            stage_f.stage_f_root(root) / "manifest.json")["sha256"],
        "stage_f_policy_changed": False,
        "method_selection_performed": False,
        "task_count": len(tasks),
    }
    stage.write_json(output / "manifest.json", manifest)
    stage.write_json(output / "state.json", {
        "status": "prepared", "completed_reductions": 0,
        "completed_analyses": 0, "task_count": len(tasks),
        "stage_f_policy_changed": False, "method_selection_performed": False,
    })
    print(output / "manifest.json", flush=True)
    print(f"taskset -c 12-27 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "
          f"MKL_NUM_THREADS=1 python3 {frozen_runner} run {root} --workers 4",
          flush=True)


def replace_fingerprint(records: list[dict[str, object]], path: Path) -> None:
    """Replace the unique manifest fingerprint for one repaired file."""
    matches = [index for index, record in enumerate(records)
               if Path(str(record["path"])) == path]
    stage.require(len(matches) == 1, f"manifest has no unique record for {path}")
    records[matches[0]] = stage.fingerprint(path)


def repair_analysis_aperture(root: Path) -> None:
    """Upgrade a failed initial Stage-G analysis while retaining its reductions."""
    root = root.resolve()
    output = stage_g_root(root)
    protocol_path = output / "protocol.json"
    manifest_path = output / "manifest.json"
    frozen_runner = output / "software" / SCRIPT_PATH.name
    stage.require(output.is_dir() and protocol_path.is_file() and manifest_path.is_file(),
                  "prepared Stage-G output is missing")
    stage.require(not (output / "complete.json").exists(),
                  "completed Stage-G output must not be repaired")
    manifest = read(manifest_path)
    protocol = read(protocol_path)
    stage.verify(manifest["input_records"] + manifest["software_records"] +
                 manifest["model_receipts"])
    stage.require(len(protocol["tasks"]) == SITE_COUNT * len(ARM_DEFINITIONS),
                  "Stage-G repair task count changed")
    receipts = [output / "reductions" / str(task["name"]) / "complete.json"
                for task in protocol["tasks"]]
    stage.require(all(path.is_file() for path in receipts),
                  "Stage-G repair requires all reductions to be complete")
    for path in receipts:
        verify_task_receipt(path)

    repair_directory = output / "repairs" / "analysis_aperture"
    stage.require(not repair_directory.exists(),
                  f"Stage-G analysis-aperture repair already exists: {repair_directory}")
    repair_directory.mkdir(parents=True)
    for path in (protocol_path, manifest_path, frozen_runner):
        shutil.copy2(path, repair_directory / (path.name + ".before"))
    before = {path.name: stage.fingerprint(repair_directory / (path.name + ".before"))
              for path in (protocol_path, manifest_path, frozen_runner)}

    contract = protocol["analysis_contract"]
    stage.require(float(contract["aperture_radius_pixels"]) ==
                  stage_f.EXPECTED_SETTINGS["aperture_radius"],
                  "Stage-G scientific aperture changed before repair")
    contract["hcianalyze_reporting_aperture_radius_pixels"] = (
        ANALYZER_REPORTING_APERTURE_RADIUS)
    contract["hcianalyze_reporting_aperture_affects_annular_normalization"] = False
    contract["analysis_repair"] = (
        "The reporting aperture is 60 pixels so the exclusion-only known-planet "
        "entry has valid pixels; source measurements retain the reviewed 3-pixel aperture."
    )
    stage.write_json(protocol_path, protocol)
    shutil.copy2(SCRIPT_PATH, frozen_runner)
    replace_fingerprint(manifest["input_records"], protocol_path)
    replace_fingerprint(manifest["software_records"], frozen_runner)
    manifest["repairs"] = [{
        "name": "analysis_aperture",
        "reason": (
            "hciAnalyze also scores exclusion entries and rejected the masked "
            "known-planet aperture"
        ),
        "reductions_reused": len(receipts),
        "scientific_aperture_radius_pixels": contract["aperture_radius_pixels"],
        "hcianalyze_reporting_aperture_radius_pixels":
            ANALYZER_REPORTING_APERTURE_RADIUS,
        "before": before,
    }]
    stage.write_json(manifest_path, manifest)
    stage.verify(manifest["input_records"] + manifest["software_records"] +
                 manifest["model_receipts"])
    stage.write_json(output / "state.json", {
        "status": "analysis_repair_prepared", "completed_reductions": len(receipts),
        "completed_analyses": 0, "task_count": len(protocol["tasks"]),
        "stage_f_policy_changed": False, "method_selection_performed": False,
    })
    print(f"repaired Stage-G analysis aperture; retained {len(receipts)} reductions",
          flush=True)
    print(f"taskset -c 12-27 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "
          f"MKL_NUM_THREADS=1 python3 {frozen_runner} run {root} --workers 4",
          flush=True)


def enable(root: Path) -> tuple[dict[str, object], dict[str, object], Path]:
    """Verify the Stage-G manifest and return its protocol and analyzer."""
    root = root.resolve()
    output = stage_g_root(root)
    protocol = read(output / "protocol.json")
    manifest = read(output / "manifest.json")
    contract = protocol["analysis_contract"]
    stage.require(protocol["primary_mode"] == PRIMARY_MODE and
                  protocol["primary_arm"] == PRIMARY_ARM and
                  tuple(protocol["methods"]) == stage_f.PLANET_METHODS and
                  len(protocol["tasks"]) == SITE_COUNT * len(ARM_DEFINITIONS) and
                  float(contract["hcianalyze_reporting_aperture_radius_pixels"]) ==
                  ANALYZER_REPORTING_APERTURE_RADIUS and
                  not contract["hcianalyze_reporting_aperture_affects_annular_normalization"] and
                  not manifest["stage_f_policy_changed"] and
                  not manifest["method_selection_performed"],
                  "Stage-G frozen contract changed")
    stage.verify(manifest["input_records"] + manifest["software_records"] +
                 manifest["model_receipts"])
    stage.require(stage_f.model_receipts(root) == manifest["model_receipts"],
                  "Stage-G calibration model products changed")
    stage.require(stage.fingerprint(SCRIPT_PATH) in manifest["software_records"],
                  "Stage-G runner is outside the frozen software receipt")
    parent, analyzer, _ = verify_stage_f_boundary(root)
    stage.require(analyzer == Path(str(protocol["paths"]["hcianalyze"])),
                  "Stage-G analyzer no longer matches the Stage-F boundary")
    return protocol, parent, analyzer


def validate_reduction(path: Path) -> dict[str, object]:
    """Validate one Stage-G injected KLIP cube."""
    data, header = fits.getdata(path, header=True, memmap=True)
    stage.require(np.asarray(data).shape == (len(stage.MODES), 128, 128) and
                  stage.read_modes(header) == stage.MODES,
                  "Stage-G reduction schema changed")
    return stage.fingerprint(path)


def verify_task_receipt(path: Path) -> None:
    """Verify one complete reduction or analysis receipt."""
    receipt = read(path)
    stage.require(receipt["status"] == "complete", f"incomplete receipt: {path}")
    stage.verify(receipt["products"])


def reduce_all(root: Path, protocol: dict[str, object], parent: dict[str, object]) -> None:
    """Run or verify every frozen Stage-G injection reduction."""
    output = stage_g_root(root)
    reductions = output / "reductions"
    reductions.mkdir(exist_ok=True)
    tasks = list(protocol["tasks"])
    for index, task in enumerate(tasks, start=1):
        directory = reductions / str(task["name"])
        complete = directory / "complete.json"
        if complete.exists():
            verify_task_receipt(complete)
        else:
            if directory.exists():
                development.archive_incomplete(output, directory, "reductions")
            directory.mkdir()
            elapsed = stage.run_command([str(value) for value in task["command"]], directory,
                                        directory / "run.log",
                                        development.environment(parent, True))
            product = validate_reduction(directory / "finim.fits")
            stage.write_json(complete, {
                "status": "complete", "task": task["name"],
                "elapsed_seconds": elapsed,
                "products": [product, stage.fingerprint(directory / "run.log"),
                             stage.fingerprint(directory / "run.command.json")],
            })
        stage.write_json(output / "state.json", {
            "status": "reducing", "completed_reductions": index,
            "completed_analyses": len(list((output / "analysis").glob("*/complete.json")))
                if (output / "analysis").exists() else 0,
            "task_count": len(tasks), "stage_f_policy_changed": False,
            "method_selection_performed": False,
        })
        print(f"Stage-G reduction {index}/{len(tasks)}: {task['name']}", flush=True)


def production_snr(protocol: dict[str, object], parent: dict[str, object], analyzer: Path,
                   task: dict[str, object], maps: np.ndarray, header: fits.Header,
                   output: Path) -> tuple[np.ndarray, dict[str, object], list[str]]:
    """Run production annular normalization with known-planet and trial exclusions."""
    flat = maps.reshape((-1, *maps.shape[-2:])).astype(np.float32)
    working_header = header.copy()
    working_header["HCI FILTER LABELS"] = ",".join(
        f"m{mode}_{method}" for mode in stage.MODES for method in stage_f.PLANET_METHODS)
    fits.writeto(output / "amplitudes.fits", flat, working_header)
    known = parent["known_planet"]
    source = task["source"]
    settings = task["settings"]
    command = [
        str(analyzer), "--config", str(stage_g_root(Path(str(protocol["parent_root"]))) /
                                         "analyze.conf"),
        "--file=amplitudes.fits", f"--lambdaD={settings['lambda_d']}",
        f"--planet.sep={known['separation']},{source['separation']}",
        f"--planet.PA={known['position_angle']},{source['position_angle']}",
        f"--planet.R={known['exclusion_radius']},{settings['source_radius']}",
        f"--snr.apertureR={ANALYZER_REPORTING_APERTURE_RADIUS}",
        f"--snr.minRad={settings['minimum_radius']}",
        f"--snr.maxRad={settings['maximum_radius']}",
        "--filter.psfResponse=", "--filter.lpfGaussFW=0", "--filter.hpfGaussFW=0",
        "--noise.model=identity", "--noise.only=false", "--noise.outputDiagnostics=false",
    ]
    with (output / "analysis.log").open("w", encoding="utf-8") as log:
        subprocess.run(command, cwd=output, env=development.environment(parent, False),
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    snr, snr_header = fits.getdata(output / "amplitudes_snr.fits", header=True)
    snr = np.asarray(snr, dtype=np.float64).reshape(maps.shape)
    stage.require(int(snr_header["SNRAPER"]) ==
                  int(ANALYZER_REPORTING_APERTURE_RADIUS) and
                  int(snr_header["SNRMINR"]) == int(settings["minimum_radius"]) and
                  int(snr_header["SNRMAXR"]) == int(settings["maximum_radius"]) and
                  int(snr_header["SNRMEAN"]) == 1 and int(snr_header["SNRSMALL"]) == 1,
                  "Stage-G hciAnalyze settings changed")

    known_position = source_from_polar(maps.shape[-2:], float(known["separation"]),
                                       float(known["position_angle"]))
    trial_position = (float(source["row"]), float(source["column"]))
    exclusions = [
        (known_position[0], known_position[1], float(known["exclusion_radius"])),
        (trial_position[0], trial_position[1], float(settings["source_radius"])),
    ]
    radius_map = development.image_radius(maps.shape[-2:])
    maximum = 0.0
    records = {}
    for mode_index, mode in enumerate(stage.MODES):
        mode_record = {}
        for method_index, method in enumerate(stage_f.PLANET_METHODS):
            expected = development.annular_oracle(maps[mode_index, method_index], exclusions)
            expected[(radius_map < float(settings["minimum_radius"])) |
                     (radius_map > float(settings["maximum_radius"]))] = np.nan
            valid = np.isfinite(maps[mode_index, method_index]) & np.isfinite(expected)
            stage.require(np.any(valid) and np.all(np.isfinite(snr[mode_index, method_index][valid])) and
                          np.allclose(snr[mode_index, method_index][valid], expected[valid],
                                      rtol=2e-6, atol=2e-6),
                          f"Stage-G annular oracle mismatch for {mode} {method}")
            error = float(np.max(np.abs(snr[mode_index, method_index][valid] - expected[valid])))
            maximum = max(maximum, error)
            mode_record[method] = {"maximum_error": error,
                                   "checked_pixels": int(np.count_nonzero(valid))}
        records[str(mode)] = mode_record
    oracle = {"rtol": 2e-6, "atol": 2e-6,
              "maximum_production_oracle_error": maximum,
              "known_planet_exclusion_row_column_radius": [*known_position,
                                                              known["exclusion_radius"]],
              "trial_exclusion_row_column_radius": [*trial_position,
                                                        settings["source_radius"]],
              "methods_by_mode": records}
    return snr, oracle, command


def response_fidelity(root: Path, parent: dict[str, object], task: dict[str, object],
                      positive: np.ndarray, diagnostics: dict[str, object]) \
        -> dict[str, object]:
    """Compare the finite primary-mode injected response with the nominal exact response."""
    mode_index = stage.MODES.index(PRIMARY_MODE)
    baseline = np.asarray(fits.getdata(parent["paths"]["baseline"], memmap=True),
                          dtype=np.float64)
    response_paths = development.response47.product_paths(Path(str(parent["parent_response"])))
    coordinates = np.asarray(fits.getdata(response_paths["coordinates"], memmap=True),
                             dtype=np.float64).T
    lookup = {(int(row), int(column)): index
              for index, (row, column, _, _) in enumerate(coordinates)}
    responses = np.asarray(fits.getdata(response_paths["responses"][mode_index], memmap=True),
                           dtype=np.float64)
    validities = np.asarray(fits.getdata(response_paths["validities"][mode_index], memmap=True),
                            dtype=np.float64) > 0.5
    searches = [tuple(map(int, query)) for query in diagnostics["aperture_queries"]]
    source = (float(task["source"]["row"]), float(task["source"]["column"]))
    query = min(searches, key=lambda value: math.hypot(value[0] - source[0],
                                                        value[1] - source[1]))
    source_index = lookup[query]
    template = development.raw.crop(responses[source_index].T, development.SUPPORT)
    validity = development.raw.crop(validities[source_index].T, development.SUPPORT)
    optimized = source_from_polar((128, 128), float(parent["known_planet"]["separation"]),
                                  float(parent["known_planet"]["position_angle"]))
    nominal = stage_f.policy_radius(math.hypot(query[0] - 63.5, query[1] - 63.5))
    width = int(parent["training"]["candidate_specific_half_width_by_radius"][str(nominal)])
    fitted = stage_f.fit_planet_query(baseline[mode_index], query, searches, template,
                                      validity, optimized, width)
    delta = ((development.stamp(positive[mode_index], query) -
              development.stamp(baseline[mode_index], query)) / float(task["contrast"]))
    metrics = validation.response_fidelity(delta, template, fitted)
    return {"mode": PRIMARY_MODE, "query_row_column": list(query),
            "source_offset_from_query_pixels": [source[0] - query[0], source[1] - query[1]],
            "metrics": metrics}


def analyze_task(root_value: str, task_name: str) -> str:
    """Analyze one Stage-G injection with the complete Stage-F aperture contract."""
    root = Path(root_value)
    protocol, parent, analyzer = enable(root)
    task = next(record for record in protocol["tasks"] if record["name"] == task_name)
    output = stage_g_root(root) / "analysis" / task_name
    complete = output / "complete.json"
    if complete.exists():
        verify_task_receipt(complete)
        return str(output / "results.json")
    if output.exists():
        development.archive_incomplete(stage_g_root(root), output, "analysis")
    output.mkdir(parents=True)
    started = time.monotonic()
    reduction = stage_g_root(root) / "reductions" / task_name / "finim.fits"
    validate_reduction(reduction)
    positive, header = fits.getdata(reduction, header=True, memmap=True)
    positive = np.asarray(positive, dtype=np.float64)
    settings = dict(task["settings"])
    maps, responses, policy_map, diagnostics, aperture, bins = stage_f.build_amplitude_maps(
        root, parent, positive, settings)
    map_header = header.copy()
    map_header["HCI RADIAL BINS"] = ",".join(map(str, bins))
    fits.writeto(output / "responses.fits",
                 responses.reshape((-1, *responses.shape[-2:])).astype(np.float32), map_header)
    fits.writeto(output / "policy_radius.fits", policy_map.astype(np.float32), map_header)
    stage.write_json(output / "diagnostics.json", diagnostics)
    snr, oracle, command = production_snr(protocol, parent, analyzer, task, maps,
                                           map_header, output)
    stage.write_json(output / "annular_verification.json", oracle)
    stage.write_json(output / "command.json", command)
    result = stage_f.summarize_planet(maps, responses, snr, diagnostics, aperture, settings)
    result.update({
        "purpose": "Stage-G exact-contrast injection analyzed with the Stage-F aperture",
        "known_planet_opened": True,
        "injection_task": {key: value for key, value in task.items() if key != "command"},
        "response_fidelity": response_fidelity(root, parent, task, positive, diagnostics),
        "annular_oracle_maximum_error": oracle["maximum_production_oracle_error"],
        "elapsed_seconds": time.monotonic() - started,
    })
    stage.write_json(output / "results.json", result)
    products = [stage.fingerprint(output / name) for name in (
        "amplitudes.fits", "amplitudes_snr.fits", "responses.fits", "policy_radius.fits",
        "diagnostics.json", "annular_verification.json", "analysis.log", "command.json",
        "results.json")]
    stage.write_json(complete, {"status": "complete", "task": task_name,
                                "products": products})
    verify_task_receipt(complete)
    return str(output / "results.json")


def analyze_all(root: Path, protocol: dict[str, object], workers: int) -> None:
    """Run or verify every Stage-G aperture analysis."""
    output = stage_g_root(root)
    analysis = output / "analysis"
    analysis.mkdir(exist_ok=True)
    tasks = list(protocol["tasks"])
    pending = []
    for task in tasks:
        reduction_receipt = output / "reductions" / str(task["name"]) / "complete.json"
        verify_task_receipt(reduction_receipt)
        complete = analysis / str(task["name"]) / "complete.json"
        if complete.exists():
            verify_task_receipt(complete)
        else:
            pending.append(str(task["name"]))
    completed = len(tasks) - len(pending)
    if pending:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(analyze_task, str(root), name): name for name in pending}
            for future in as_completed(futures):
                name = futures[future]
                future.result()
                completed += 1
                stage.write_json(output / "state.json", {
                    "status": "analyzing", "completed_reductions": len(tasks),
                    "completed_analyses": completed, "task_count": len(tasks),
                    "stage_f_policy_changed": False, "method_selection_performed": False,
                })
                print(f"Stage-G analysis {completed}/{len(tasks)}: {name}", flush=True)


def paired_summary(records: list[dict[str, object]], first: str, second: str,
                   field: str) -> dict[str, object]:
    """Summarize one paired method difference over common injection sites."""
    values = [float(record["methods_by_mode"][str(PRIMARY_MODE)][first][field]) -
              float(record["methods_by_mode"][str(PRIMARY_MODE)][second][field])
              for record in records]
    result = sample_summary(values)
    result["first_method"] = first
    result["second_method"] = second
    result["field"] = field
    result["values_by_site"] = {
        str(record["injection_task"]["base_site"]): value
        for record, value in zip(records, values)
    }
    return result


def summarize(root: Path, protocol: dict[str, object]) -> None:
    """Aggregate Stage-G arms and compare them with the Stage-F planet."""
    output = stage_g_root(root)
    records = []
    for task in protocol["tasks"]:
        complete = output / "analysis" / str(task["name"]) / "complete.json"
        verify_task_receipt(complete)
        records.append(read(output / "analysis" / str(task["name"]) / "results.json"))
    by_arm: dict[str, list[dict[str, object]]] = defaultdict(list)
    for record in records:
        by_arm[str(record["injection_task"]["arm"])].append(record)
    stage.require(set(by_arm) == {str(arm["name"]) for arm in ARM_DEFINITIONS} and
                  all(len(values) == SITE_COUNT for values in by_arm.values()),
                  "Stage-G summary lost an arm or site")
    for values in by_arm.values():
        values.sort(key=lambda record: int(record["injection_task"]["site_index"]))

    planet = read(stage_f.stage_f_root(root) / "planet" / "results.json")
    planet_methods = planet["methods_by_mode"][str(PRIMARY_MODE)]
    arms = {}
    for arm_name, arm_records in sorted(by_arm.items()):
        methods = {}
        for method in stage_f.PLANET_METHODS:
            nearest = sample_summary([
                float(record["methods_by_mode"][str(PRIMARY_MODE)][method]["nearest_pixel_snr"])
                for record in arm_records])
            aperture = sample_summary([
                float(record["methods_by_mode"][str(PRIMARY_MODE)][method]["aperture_maximum_snr"])
                for record in arm_records])
            nearest["planet_comparison"] = normal_predictive(
                nearest, float(planet_methods[method]["nearest_pixel_snr"]))
            aperture["planet_comparison"] = normal_predictive(
                aperture, float(planet_methods[method]["aperture_maximum_snr"]))
            methods[method] = {"nearest_pixel_snr": nearest,
                               "aperture_maximum_snr": aperture}
        pairs = {}
        for name, first, second in PAIR_DEFINITIONS:
            pair = {}
            for field in ("nearest_pixel_snr", "aperture_maximum_snr"):
                summary = paired_summary(arm_records, first, second, field)
                observed = (float(planet_methods[first][field]) -
                            float(planet_methods[second][field]))
                summary["planet_comparison"] = normal_predictive(summary, observed)
                pair[field] = summary
            pairs[name] = pair
        fidelity = {}
        for metric_name in ("unweighted", "raw_rectangular_m0p3", "radial_hann_m0p1"):
            fidelity[metric_name] = {}
            for value_name in ("cosine", "projection_scale", "best_scaled_relative_residual"):
                fidelity[metric_name][value_name] = sample_summary([
                    float(record["response_fidelity"]["metrics"][metric_name][value_name])
                    for record in arm_records])
        arms[arm_name] = {"methods": methods, "paired_differences": pairs,
                          "response_fidelity": fidelity}

    by_site_arm = {(str(record["injection_task"]["base_site"]),
                    str(record["injection_task"]["arm"])): record for record in records}
    arm_effects = {}
    reference_arm = "nominal_optimized_phase"
    for arm in (str(value["name"]) for value in ARM_DEFINITIONS):
        if arm == reference_arm:
            continue
        method_effects = {}
        for method in stage_f.PLANET_METHODS:
            fields = {}
            for field in ("nearest_pixel_snr", "aperture_maximum_snr"):
                values = []
                for site in (str(site["name"]) for site in protocol["sites"]):
                    first = by_site_arm[(site, arm)]
                    second = by_site_arm[(site, reference_arm)]
                    values.append(float(first["methods_by_mode"][str(PRIMARY_MODE)][method][field]) -
                                  float(second["methods_by_mode"][str(PRIMARY_MODE)][method][field]))
                fields[field] = sample_summary(values)
            method_effects[method] = fields
        arm_effects[f"{arm}_minus_{reference_arm}"] = method_effects

    result = {
        "schema": 1,
        "purpose": "direct optimized-contrast planet/injection consistency and PSF-model controls",
        "primary_mode": PRIMARY_MODE,
        "site_count_per_arm": SITE_COUNT,
        "primary_arm": PRIMARY_ARM,
        "planet_source": planet["source_row_column"],
        "planet_methods": planet_methods,
        "arms": arms,
        "paired_arm_effects": arm_effects,
        "limits": [
            "All injections share one observing sequence and correlated residual field.",
            "Gaussian broadening is a controlled PSF perturbation, not a unique physical model.",
            "Normal predictive probabilities are diagnostic; paired values and ranks are primary.",
        ],
    }
    stage.write_json(output / "results.json", result)

    csv_path = output / "results.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["arm", "method", "nearest_mean", "nearest_sample_std",
                         "aperture_mean", "aperture_sample_std", "planet_nearest",
                         "planet_aperture"])
        for arm_name, arm in arms.items():
            for method, values in arm["methods"].items():
                nearest = values["nearest_pixel_snr"]
                aperture = values["aperture_maximum_snr"]
                writer.writerow([arm_name, method, nearest["mean"],
                                 nearest["sample_standard_deviation"], aperture["mean"],
                                 aperture["sample_standard_deviation"],
                                 planet_methods[method]["nearest_pixel_snr"],
                                 planet_methods[method]["aperture_maximum_snr"]])

    lines = [
        "# KLIP Stage-G exact-contrast injection result", "",
        "Twelve score-blind radius-12 sites are measured in each arm. Values below are",
        "mode-200 aperture-maximum SNR means ± sample standard deviations; the planet",
        "column is the unchanged Stage-F value.", "",
        "| Arm | Method | Injections | Planet | Planet standardized deviation |",
        "| :--- | :--- | ---: | ---: | ---: |",
    ]
    for arm_name, arm in arms.items():
        for method in stage_f.PLANET_METHODS:
            summary = arm["methods"][method]["aperture_maximum_snr"]
            comparison = summary["planet_comparison"]
            lines.append(f"| {arm_name} | {stage_f.METHOD_LABELS[method]} | "
                         f"{summary['mean']:.4f} ± {summary['sample_standard_deviation']:.4f} | "
                         f"{comparison['observed']:.4f} | "
                         f"{comparison['standardized_deviation']:+.2f} |")
    lines.extend(["", "## Paired covariance endpoint", "",
                  "| Arm | Difference | Injections, mean ± sample std | Planet | Standardized difference |",
                  "| :--- | :--- | ---: | ---: | ---: |"])
    for arm_name, arm in arms.items():
        for name in ("raw_rectangular_m0p3_minus_exact_identity",
                     "radial_hann_m0p1_trunc0p75_minus_exact_identity"):
            summary = arm["paired_differences"][name]["aperture_maximum_snr"]
            comparison = summary["planet_comparison"]
            lines.append(f"| {arm_name} | {name} | {summary['mean']:+.4f} ± "
                         f"{summary['sample_standard_deviation']:.4f} | "
                         f"{comparison['observed']:+.4f} | "
                         f"{comparison['standardized_deviation']:+.2f} |")
    lines.extend(["", "## Paired arm effects", "",
                  "Each cell is the mean aperture-maximum SNR change relative to the nominal",
                  "optimized-phase arm at the same 12 base sites.", "",
                  "| Arm minus nominal optimized phase | Gaussian 3.6 | Exact identity | Raw covariance | Radial covariance |",
                  "| :--- | ---: | ---: | ---: | ---: |"])
    effect_methods = ("gaussian_fwhm3p6", "exact_identity",
                      "raw_rectangular_m0p3", "radial_hann_m0p1_trunc0p75")
    for comparison, effects in arm_effects.items():
        values = [float(effects[method]["aperture_maximum_snr"]["mean"])
                  for method in effect_methods]
        lines.append("| " + comparison + " | " +
                     " | ".join(f"{value:+.4f}" for value in values) + " |")
    lines.extend(["", "## Nominal-template response fidelity", "",
                  "| Arm | Unweighted cosine | Unweighted residual | Raw-covariance residual | Radial-covariance residual |",
                  "| :--- | ---: | ---: | ---: | ---: |"])
    for arm_name, arm in arms.items():
        fidelity = arm["response_fidelity"]
        cosine = fidelity["unweighted"]["cosine"]
        unweighted = fidelity["unweighted"]["best_scaled_relative_residual"]
        raw = fidelity["raw_rectangular_m0p3"]["best_scaled_relative_residual"]
        radial = fidelity["radial_hann_m0p1"]["best_scaled_relative_residual"]
        lines.append(f"| {arm_name} | {cosine['mean']:.6f} ± "
                     f"{cosine['sample_standard_deviation']:.6f} | "
                     f"{unweighted['mean']:.4f} ± {unweighted['sample_standard_deviation']:.4f} | "
                     f"{raw['mean']:.4f} ± {raw['sample_standard_deviation']:.4f} | "
                     f"{radial['mean']:.4f} ± {radial['sample_standard_deviation']:.4f} |")
    lines.extend(["", "The nominal optimized-phase arm is the primary endpoint. Other phase",
                  "and PSF arms are predeclared controls and do not select a filter or tune a",
                  "PSF from the observed planet SNR.", ""])
    (output / "results.md").write_text("\n".join(lines), encoding="utf-8")

    task_receipts = [stage.fingerprint(output / "analysis" / str(task["name"]) /
                                       "complete.json") for task in protocol["tasks"]]
    reduction_receipts = [stage.fingerprint(output / "reductions" / str(task["name"]) /
                                            "complete.json") for task in protocol["tasks"]]
    products = [stage.fingerprint(output / name)
                for name in ("results.json", "results.csv", "results.md")]
    stage.write_json(output / "complete.json", {
        "status": "complete", "task_count": len(protocol["tasks"]),
        "stage_f_policy_changed": False, "method_selection_performed": False,
        "products": products, "reduction_receipts": reduction_receipts,
        "analysis_receipts": task_receipts,
    })
    stage.write_json(output / "state.json", {
        "status": "complete", "completed_reductions": len(protocol["tasks"]),
        "completed_analyses": len(protocol["tasks"]), "task_count": len(protocol["tasks"]),
        "stage_f_policy_changed": False, "method_selection_performed": False,
    })


def verify_complete(root: Path) -> None:
    """Verify the complete Stage-G receipt and all nested task products."""
    output = stage_g_root(root)
    receipt = read(output / "complete.json")
    stage.require(receipt["status"] == "complete" and
                  not receipt["stage_f_policy_changed"] and
                  not receipt["method_selection_performed"],
                  "Stage-G completion receipt changed")
    stage.verify(receipt["products"] + receipt["reduction_receipts"] +
                 receipt["analysis_receipts"])
    for record in receipt["reduction_receipts"] + receipt["analysis_receipts"]:
        verify_task_receipt(Path(str(record["path"])))


def run(root: Path, workers: int) -> None:
    """Run the resumable reductions, aperture analyses, and frozen summary."""
    root = root.resolve()
    protocol, parent, _ = enable(root)
    output = stage_g_root(root)
    if (output / "complete.json").exists():
        verify_complete(root)
        print(output / "results.md", flush=True)
        return
    stage.require(workers >= 1, "workers must be positive")
    reduce_all(root, protocol, parent)
    analyze_all(root, protocol, workers)
    summarize(root, protocol)
    verify_complete(root)
    print(output / "results.md", flush=True)


def check(config: Path) -> None:
    """Check geometry conversions, phase definitions, PSF controls, and method pairs."""
    settings = stage_f.parse_analysis_config(config.resolve())
    stage.require(settings == stage_f.EXPECTED_SETTINGS,
                  "Stage-G check requires the reviewed analysis settings")
    analysis = source_from_polar((128, 128), settings["separation"],
                                 settings["position_angle"])
    optimized = source_from_polar((128, 128), stage.PLANET_SEPARATION, stage.PLANET_PA)
    analysis_phase = fractional_phase(analysis)
    optimized_phase = fractional_phase(optimized)
    stage.require(np.allclose(analysis_phase, [0.16879332335545882, -0.12935152033411157],
                              rtol=0, atol=1e-12) and
                  np.allclose(optimized_phase, [-0.2770699293779586, 0.4859637006126931],
                              rtol=0, atol=1e-12),
                  "Stage-G source phases changed")
    for source in (analysis, optimized, (52.25, 61.75)):
        separation, angle = polar_from_source((128, 128), *source)
        stage.require(np.allclose(source_from_polar((128, 128), separation, angle), source,
                                  rtol=0, atol=2e-12),
                      "Stage-G source conversion failed")
    stage.require(len(ARM_DEFINITIONS) == len({arm["name"] for arm in ARM_DEFINITIONS}) and
                  PRIMARY_ARM in {arm["name"] for arm in ARM_DEFINITIONS} and
                  set(PSF_BLUR_FWHM) == {arm["psf"] for arm in ARM_DEFINITIONS} and
                  all(first in stage_f.PLANET_METHODS and second in stage_f.PLANET_METHODS
                      for _, first, second in PAIR_DEFINITIONS),
                  "Stage-G arm or method contract changed")
    synthetic = np.zeros((31, 31), dtype=np.float64)
    synthetic[15, 15] = 1.0
    for fwhm in (0.9, 1.8):
        blurred = gaussian_filter(synthetic, sigma=fwhm / math.sqrt(8 * math.log(2)),
                                  mode="constant", cval=0.0, truncate=4.0)
        stage.require(np.isclose(np.sum(blurred), 1.0, rtol=2e-7, atol=2e-7) and
                      0 < np.max(blurred) <= 1,
                      "Stage-G Gaussian PSF control failed")
    print("KLIP Stage-G planet-consistency checks passed", flush=True)


def parser() -> argparse.ArgumentParser:
    """Build the Stage-G command-line parser."""
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="action", required=True)
    check_parser = subparsers.add_parser("check")
    check_parser.add_argument("--config", type=Path, default=Path("working/analyze.conf"))
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("root", type=Path)
    prepare_parser.add_argument("--config", type=Path, default=Path("working/analyze.conf"))
    repair_parser = subparsers.add_parser("repair-analysis-aperture")
    repair_parser.add_argument("root", type=Path)
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("root", type=Path)
    run_parser.add_argument("--workers", type=int, default=4)
    return result


def main() -> None:
    """Dispatch the requested Stage-G action."""
    arguments = parser().parse_args()
    if arguments.action == "check":
        check(arguments.config)
    elif arguments.action == "prepare":
        prepare(arguments.root, arguments.config)
    elif arguments.action == "repair-analysis-aperture":
        repair_analysis_aperture(arguments.root)
    else:
        run(arguments.root, arguments.workers)


if __name__ == "__main__":
    main()
