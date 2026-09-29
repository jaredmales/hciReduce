#!/usr/bin/env python3
"""Reanalyze the KLIP planet and matched injections with cubic-shifted responses."""
from __future__ import annotations

import argparse
import csv
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
from scipy.stats import t as student_t

SCRIPT_PATH = Path(__file__).resolve()
SCRIPT_DIRECTORY = SCRIPT_PATH.parent
FROZEN_DEPENDENCIES = (SCRIPT_PATH.parents[2] / "software"
                       if len(SCRIPT_PATH.parents) > 2 else Path("/nonexistent"))
if (FROZEN_DEPENDENCIES / "run_klip_stage_h_subpixel_response.py").is_file():
    sys.path.insert(0, str(FROZEN_DEPENDENCIES))
else:
    sys.path.insert(0, str(SCRIPT_DIRECTORY))

import run_klip_stage_h_subpixel_response as stage_h  # noqa: E402


stage_g = stage_h.stage_g
stage_f = stage_h.stage_f
development = stage_h.development
stage = stage_h.stage

STAGE_DIRECTORY = "stage_i_shifted_planet"
PRIMARY_MODE = stage_h.PRIMARY_MODE
MODE_INDEX = stage.MODES.index(PRIMARY_MODE)
NOMINAL_RADIUS = stage_g.NOMINAL_RADIUS
SITE_COUNT = stage_g.SITE_COUNT
ANALYZER_REPORTING_APERTURE_RADIUS = stage_g.ANALYZER_REPORTING_APERTURE_RADIUS
SHIFTED_METHODS = (
    "shifted_exact_identity",
    "shifted_raw_rectangular_m0p3",
    "shifted_radial_hann_m0p1_trunc0p75",
)
METHODS = (
    "gaussian_fwhm2p4",
    "gaussian_fwhm3p6",
    "integer_exact_identity",
    "shifted_exact_identity",
    "integer_raw_rectangular_m0p3",
    "shifted_raw_rectangular_m0p3",
    "integer_radial_hann_m0p1_trunc0p75",
    "shifted_radial_hann_m0p1_trunc0p75",
)
CURRENT_METHODS = {
    "gaussian_fwhm2p4": "gaussian_fwhm2p4",
    "gaussian_fwhm3p6": "gaussian_fwhm3p6",
    "integer_exact_identity": "exact_identity",
    "integer_raw_rectangular_m0p3": "raw_rectangular_m0p3",
    "integer_radial_hann_m0p1_trunc0p75": "radial_hann_m0p1_trunc0p75",
}
SHIFTED_TO_CURRENT = {
    "shifted_exact_identity": "exact_identity",
    "shifted_raw_rectangular_m0p3": "raw_rectangular_m0p3",
    "shifted_radial_hann_m0p1_trunc0p75": "radial_hann_m0p1_trunc0p75",
}
PAIRS = (
    ("shifted_identity_minus_integer_identity",
     "shifted_exact_identity", "integer_exact_identity"),
    ("shifted_raw_minus_integer_raw",
     "shifted_raw_rectangular_m0p3", "integer_raw_rectangular_m0p3"),
    ("shifted_radial_minus_integer_radial",
     "shifted_radial_hann_m0p1_trunc0p75",
     "integer_radial_hann_m0p1_trunc0p75"),
    ("shifted_raw_minus_shifted_identity",
     "shifted_raw_rectangular_m0p3", "shifted_exact_identity"),
    ("shifted_radial_minus_shifted_identity",
     "shifted_radial_hann_m0p1_trunc0p75", "shifted_exact_identity"),
    ("gaussian_fwhm3p6_minus_shifted_identity",
     "gaussian_fwhm3p6", "shifted_exact_identity"),
    ("gaussian_fwhm3p6_minus_gaussian_fwhm2p4",
     "gaussian_fwhm3p6", "gaussian_fwhm2p4"),
)


def read(path: Path) -> object:
    """Read one JSON document."""
    return json.loads(path.read_text(encoding="utf-8"))


def output_root(root: Path) -> Path:
    """Return the isolated Stage-I output directory."""
    return root / STAGE_DIRECTORY


def response_paths(parent: dict[str, object]) -> dict[str, object]:
    """Return the promoted exact-response products used by Stage I."""
    return development.response47.product_paths(Path(str(parent["parent_response"])))


def source_settings(separation: float, position_angle: float) -> dict[str, float]:
    """Return the Stage-F settings at one explicitly selected source position."""
    settings = dict(stage_f.EXPECTED_SETTINGS)
    settings["separation"] = float(separation)
    settings["position_angle"] = float(position_angle)
    return settings


def verify_parents(root: Path) \
        -> tuple[dict[str, object], dict[str, object], dict[str, object], Path]:
    """Verify Stages F through H and return their frozen protocols."""
    stage_h.verify_complete(root)
    stage_g.verify_complete(root)
    stage_f.verify_complete(root)
    parent = read(root / "protocol.json")
    protocol_g = read(stage_g.stage_g_root(root) / "protocol.json")
    protocol_h = read(stage_h.output_root(root) / "protocol.json")
    analyzer = Path(str(protocol_g["paths"]["hcianalyze"]))
    stage.require(analyzer.is_file() and
                  int(protocol_g["primary_mode"]) == PRIMARY_MODE and
                  str(protocol_g["primary_arm"]) == stage_g.PRIMARY_ARM and
                  int(protocol_h["primary_mode"]) == PRIMARY_MODE and
                  np.allclose(protocol_h["phase_offset_row_column"],
                              [-0.2770699293779586, 0.4859637006126931],
                              rtol=0, atol=1e-12),
                  "Stage-I parent boundary changed")
    return parent, protocol_g, protocol_h, analyzer


def primary_injection_tasks(protocol: dict[str, object]) -> list[dict[str, object]]:
    """Return the ordered Stage-G optimized-phase injection tasks."""
    tasks = [dict(task) for task in protocol["tasks"]
             if str(task["arm"]) == stage_g.PRIMARY_ARM]
    tasks.sort(key=lambda task: int(task["site_index"]))
    stage.require(len(tasks) == SITE_COUNT and
                  len({str(task["base_site"]) for task in tasks}) == SITE_COUNT,
                  "Stage-I lost a Stage-G primary site")
    return tasks


def prepare(root: Path) -> None:
    """Freeze the shifted-template planet and injection reanalysis."""
    root = root.resolve()
    parent, protocol_g, protocol_h, analyzer = verify_parents(root)
    stage.require(sorted(os.sched_getaffinity(0)) == parent["resources"]["cpu_affinity"],
                  "Stage-I prepare must use the frozen response CPU affinity")
    output = output_root(root)
    stage.require(not output.exists(), f"Stage-I output already exists: {output}")
    output.mkdir()
    (output / "software").mkdir()
    shutil.copy2(stage_g.stage_g_root(root) / "analyze.conf", output / "analyze.conf")

    frozen_runner = output / "software" / SCRIPT_PATH.name
    shutil.copy2(SCRIPT_PATH, frozen_runner)
    dependency_sources = (
        stage_h.output_root(root) / "software" / stage_h.SCRIPT_PATH.name,
        stage_g.stage_g_root(root) / "software" / stage_g.SCRIPT_PATH.name,
        stage_f.stage_f_root(root) / "software" / stage_f.SCRIPT_PATH.name,
    )
    stage.require(all(path.is_file() for path in dependency_sources),
                  "Stage-I frozen dependency set is incomplete")
    for source in dependency_sources:
        shutil.copy2(source, output / "software" / source.name)

    science, _, _ = stage_f.original_science_path(parent)
    known = parent["known_planet"]
    optimized_settings = source_settings(float(known["separation"]),
                                         float(known["position_angle"]))
    optimized_source = stage_g.source_from_polar(
        (128, 128), optimized_settings["separation"],
        optimized_settings["position_angle"])
    optimized_phase = stage_g.fractional_phase(optimized_source)
    stage.require(np.allclose(optimized_phase, protocol_h["phase_offset_row_column"],
                              rtol=0, atol=2e-12),
                  "Stage-I optimized planet phase changed")

    tasks = [{
        "name": "planet",
        "kind": "planet",
        "base_site": "optimized_planet",
        "science": str(science),
        "parent_amplitudes": str(stage_f.stage_f_root(root) / "planet" /
                                 "amplitudes.fits"),
        "parent_policy": str(stage_f.stage_f_root(root) / "planet" /
                             "policy_radius.fits"),
        "parent_results": str(stage_f.stage_f_root(root) / "planet" /
                              "results.json"),
        "source": {
            "row": optimized_source[0], "column": optimized_source[1],
            "separation": optimized_settings["separation"],
            "position_angle": optimized_settings["position_angle"],
        },
        "settings": optimized_settings,
    }]
    for task in primary_injection_tasks(protocol_g):
        name = str(task["name"])
        tasks.append({
            "name": name,
            "kind": "injection",
            "base_site": str(task["base_site"]),
            "science": str(stage_g.stage_g_root(root) / "reductions" / name /
                           "finim.fits"),
            "parent_amplitudes": str(stage_g.stage_g_root(root) / "analysis" /
                                     name / "amplitudes.fits"),
            "parent_policy": str(stage_g.stage_g_root(root) / "analysis" /
                                 name / "policy_radius.fits"),
            "parent_results": str(stage_g.stage_g_root(root) / "analysis" /
                                  name / "results.json"),
            "source": dict(task["source"]),
            "settings": dict(task["settings"]),
        })
    stage.require(len(tasks) == SITE_COUNT + 1 and
                  all(np.allclose(stage_g.fractional_phase(
                      (float(task["source"]["row"]), float(task["source"]["column"]))),
                      optimized_phase, rtol=0, atol=2e-12)
                      for task in tasks),
                  "Stage-I task phases are not identical")

    unit_directory = (root / "calibration" / "units" /
                      development.unit_name(NOMINAL_RADIUS, PRIMARY_MODE))
    unit_results = read(unit_directory / "results.json")
    stage.require(int(unit_results["generic_positions"]) == 245 and
                  list(unit_results["required_bins"]) == [10, 11, 12, 13],
                  "Stage-I radius-12 calibration coverage changed")
    paths = response_paths(parent)
    records = [
        stage.fingerprint(stage_h.output_root(root) / "complete.json"),
        stage.fingerprint(stage_h.output_root(root) / "results.json"),
        stage.fingerprint(stage_g.stage_g_root(root) / "complete.json"),
        stage.fingerprint(stage_g.stage_g_root(root) / "protocol.json"),
        stage.fingerprint(stage_g.stage_g_root(root) / "results.json"),
        stage.fingerprint(stage_f.stage_f_root(root) / "complete.json"),
        stage.fingerprint(stage_f.stage_f_root(root) / "planet" / "results.json"),
        stage.fingerprint(Path(str(parent["paths"]["baseline"]))),
        stage.fingerprint(unit_directory / "complete.json"),
        stage.fingerprint(unit_directory / "models.npz"),
        stage.fingerprint(unit_directory / "results.json"),
        stage.fingerprint(paths["manifest"]),
        stage.fingerprint(paths["coordinates"]),
        stage.fingerprint(paths["responses"][MODE_INDEX]),
        stage.fingerprint(paths["validities"][MODE_INDEX]),
        stage.fingerprint(output / "analyze.conf"),
    ]
    for task in tasks:
        for key in ("science", "parent_amplitudes", "parent_policy",
                    "parent_results"):
            records.append(stage.fingerprint(Path(str(task[key]))))

    protocol = {
        "schema": 1,
        "stage": "KLIP Stage I shifted-template planet closure",
        "purpose": (
            "apply the promoted cubic phase shift to the optimized planet and "
            "twelve matched injections with production annular SNR"
        ),
        "parent_root": str(root),
        "primary_mode": PRIMARY_MODE,
        "nominal_radius": NOMINAL_RADIUS,
        "phase_offset_row_column": list(optimized_phase),
        "methods": list(METHODS),
        "pairs": [list(pair) for pair in PAIRS],
        "tasks": tasks,
        "primary_endpoints": [
            "optimized-nearest-pixel shifted covariance minus shifted identity SNR",
            "optimized-nearest-pixel Gaussian 3.6 minus shifted identity SNR",
            "planet comparison with the twelve injection paired-difference distributions",
        ],
        "secondary_endpoints": [
            "shifted-minus-integer SNR for each response filter",
            "nearest-pixel amplitude for integer and shifted response filters",
            "Gaussian 3.6 minus Gaussian 2.4 fixed-pixel SNR",
        ],
        "analysis_contract": {
            "source_position": "optimized negative-companion position",
            "source_phase_identical_for_planet_and_injections": True,
            "candidate_statistic": "nearest native pixel to optimized source",
            "generic_noise_weights": "radius-12 mode-200 signal-free baseline",
            "generic_calibration_positions": 245,
            "generic_calibration_bins": [10, 11, 12, 13],
            "candidate_weights": "full source-aperture exclusion on signal-free baseline",
            "annular_normalization": "production hciAnalyze with independent oracle",
            "known_planet_and_trial_excluded_from_noise": True,
            "aperture_maximum_retested": False,
            "new_klip_reductions": 0,
        },
        "resources": parent["resources"],
        "paths": {
            "baseline": parent["paths"]["baseline"],
            "parent_response": parent["parent_response"],
            "analyzer": str(analyzer),
            "unit_directory": str(unit_directory),
        },
    }
    stage.write_json(output / "protocol.json", protocol)
    records.append(stage.fingerprint(output / "protocol.json"))
    software_records = [
        stage.fingerprint(frozen_runner),
        *[stage.fingerprint(output / "software" / source.name)
          for source in dependency_sources],
        stage.fingerprint(analyzer),
    ]
    stage.write_json(output / "manifest.json", {
        "schema": 1,
        "purpose": "immutable Stage-I shifted-template planet closure",
        "root": str(root),
        "input_records": records,
        "software_records": software_records,
        "task_count": len(tasks),
        "new_klip_reductions": 0,
        "method_selection_performed": False,
        "planet_score_used_for_selection": False,
    })
    stage.write_json(output / "state.json", {
        "status": "prepared", "calibration_complete": False,
        "completed_analyses": 0, "task_count": len(tasks),
    })
    print(output / "manifest.json", flush=True)
    print(f"taskset -c 12-27 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "
          f"MKL_NUM_THREADS=1 python3 {frozen_runner} run {root}", flush=True)


def enable(root: Path) -> tuple[dict[str, object], dict[str, object], Path]:
    """Verify the frozen Stage-I boundary."""
    root = root.resolve()
    output = output_root(root)
    protocol = read(output / "protocol.json")
    manifest = read(output / "manifest.json")
    stage.require(int(protocol["primary_mode"]) == PRIMARY_MODE and
                  float(protocol["nominal_radius"]) == NOMINAL_RADIUS and
                  tuple(protocol["methods"]) == METHODS and
                  tuple(tuple(pair) for pair in protocol["pairs"]) == PAIRS and
                  len(protocol["tasks"]) == SITE_COUNT + 1 and
                  not manifest["method_selection_performed"] and
                  not manifest["planet_score_used_for_selection"] and
                  int(manifest["new_klip_reductions"]) == 0,
                  "Stage-I frozen contract changed")
    stage.verify(manifest["input_records"] + manifest["software_records"])
    stage.require(stage.fingerprint(SCRIPT_PATH) in manifest["software_records"],
                  "Stage-I runner is outside the frozen software receipt")
    stage_h.verify_complete(root)
    parent = read(root / "protocol.json")
    analyzer = Path(str(protocol["paths"]["analyzer"]))
    return protocol, parent, analyzer


def precision_weight(covariance: np.ndarray, template: np.ndarray,
                     support: np.ndarray, scale: np.ndarray,
                     kind: str, cutoff: float) -> np.ndarray:
    """Return the normalized precision weight for one dense covariance."""
    covariance = np.asarray(covariance, dtype=np.float64)
    template = np.asarray(template, dtype=np.float64)
    support = np.asarray(support, dtype=bool)
    scale = np.asarray(scale, dtype=np.float64)
    selected = np.flatnonzero(support.ravel())
    submatrix = covariance[np.ix_(selected, selected)]
    eigenvalues, eigenvectors = np.linalg.eigh(submatrix)
    target = float(np.mean(np.diag(submatrix)))
    stage.require(np.all(eigenvalues > 0) and target > 0,
                  "Stage-I precision covariance is invalid")
    transformed = (template / scale).ravel()[selected]
    coefficients = eigenvectors.T @ transformed
    threshold = cutoff * target
    if kind == "full":
        precision = 1 / eigenvalues
    elif kind == "truncate":
        precision = np.where(eigenvalues >= threshold, 1 / eigenvalues, 0)
    else:
        raise RuntimeError(f"unsupported Stage-I precision kind: {kind}")
    inverse_template = eigenvectors @ (precision * coefficients)
    energy = float(transformed @ inverse_template)
    stage.require(np.isfinite(energy) and energy > 0,
                  "Stage-I precision template has no energy")
    transformed_weight = np.zeros(template.size, dtype=np.float64)
    transformed_weight[selected] = inverse_template / energy
    weight = transformed_weight / scale.ravel()
    stage.require(np.isclose(weight @ template.ravel(), 1,
                             rtol=3e-10, atol=3e-12),
                  "Stage-I precision weight lost unit response")
    return weight


def shifted_weights(template: np.ndarray, shifted: np.ndarray,
                    support: np.ndarray, fitted: dict[str, object]) -> np.ndarray:
    """Build shifted identity, raw, and radial weights and replay current weights."""
    identity = development.normalized_weight(shifted, support)
    raw_covariance = fitted["raw_detail"]["covariance"]
    radial_covariance = fitted["radial_detail"]["covariance"]
    raw = precision_weight(raw_covariance, shifted, support,
                           np.ones(shifted.shape), "full", 0)
    radial = precision_weight(radial_covariance, shifted, support,
                              fitted["scale"], "truncate", 0.75)
    replay_raw = precision_weight(raw_covariance, template, support,
                                  np.ones(template.shape), "full", 0)
    replay_radial = precision_weight(radial_covariance, template, support,
                                     fitted["scale"], "truncate", 0.75)
    stage.require(np.allclose(replay_raw,
                              fitted["weights"]["raw_rectangular_m0p3"],
                              rtol=3e-10, atol=3e-12) and
                  np.allclose(replay_radial,
                              fitted["weights"]["radial_hann_m0p1_trunc0p75"],
                              rtol=3e-10, atol=3e-12),
                  "Stage-I precision replay differs from the parent fit")
    return np.stack((identity, raw, radial))


def verify_receipt(path: Path) -> None:
    """Verify one completed Stage-I receipt."""
    receipt = read(path)
    stage.require(receipt["status"] == "complete", f"incomplete receipt: {path}")
    stage.verify(receipt["products"])


def calibrate(root: Path, protocol: dict[str, object],
              parent: dict[str, object]) -> dict[str, np.ndarray]:
    """Build phase-shifted generic weights for the radius-12 noise annuli."""
    output = output_root(root) / "calibration"
    complete = output / "complete.json"
    if complete.exists():
        verify_receipt(complete)
        with np.load(output / "weights.npz", allow_pickle=False) as source:
            return {name: np.array(source[name], copy=True)
                    for name in source.files}
    if output.exists():
        development.archive_incomplete(output_root(root), output, "calibration")
    output.mkdir()
    unit = development.load_unit(root, NOMINAL_RADIUS, PRIMARY_MODE)
    baseline_cube = fits.getdata(protocol["paths"]["baseline"], memmap=True)
    baseline = np.asarray(baseline_cube[MODE_INDEX], dtype=np.float64)
    paths = response_paths(parent)
    responses = fits.getdata(paths["responses"][MODE_INDEX], memmap=True)
    validities = fits.getdata(paths["validities"][MODE_INDEX], memmap=True)
    planet = stage_g.source_from_polar(
        baseline.shape, float(parent["known_planet"]["separation"]),
        float(parent["known_planet"]["position_angle"]))
    phase = tuple(float(value) for value in protocol["phase_offset_row_column"])
    width = int(parent["training"]["candidate_specific_half_width_by_radius"][
        str(NOMINAL_RADIUS)])
    positions = np.asarray(unit["positions"], dtype=np.int16)
    source_indices = np.asarray(unit["source_indices"], dtype=np.int32)
    stage.require(len(positions) == len(source_indices) == 245,
                  "Stage-I generic calibration position count changed")
    stage.require(len({tuple(map(int, value)) for value in positions}) ==
                  len(positions),
                  "Stage-I generic calibration positions are not unique")
    weights = np.full((len(positions), len(SHIFTED_METHODS),
                       stage_f.SUPPORT**2), np.nan, dtype=np.float64)
    replay_errors = []
    valid_covariance = 0
    for local, (query_value, source_index) in enumerate(
            zip(positions, source_indices), start=1):
        query = tuple(map(int, query_value))
        full_template = np.asarray(responses[int(source_index)],
                                   dtype=np.float64).T
        full_validity = (np.asarray(validities[int(source_index)],
                                    dtype=np.float64).T > 0.5)
        template = development.raw.crop(full_template, stage_f.SUPPORT)
        support = development.raw.crop(full_validity, stage_f.SUPPORT)
        shifted, shifted_support = stage_h.shifted_template(
            full_template, full_validity, phase)
        stage.require(np.array_equal(support, shifted_support) and np.all(support),
                      "Stage-I generic shifted support changed")
        generic_support = support & development.candidate_mask(query, planet)
        stage.require(generic_support[stage_f.HALF, stage_f.HALF],
                      "Stage-I generic source anchor is excluded")
        weights[local - 1, 0] = development.normalized_weight(
            shifted, generic_support)
        stored_identity = unit["reference_weights"][
            local - 1,
            development.REFERENCE_WEIGHT_METHODS.index("exact_identity")]
        current_identity = development.normalized_weight(
            template, generic_support)
        replay_errors.append(float(np.max(np.abs(current_identity -
                                                   stored_identity))))
        stored_raw = unit["covariance_weights"][
            local - 1,
            development.COVARIANCE_METHODS.index("raw_rectangular_m0p3")]
        stored_radial = unit["covariance_weights"][
            local - 1,
            development.COVARIANCE_METHODS.index(
                "radial_hann_m0p1_trunc0p75")]
        if np.all(np.isfinite(stored_raw)) and np.all(np.isfinite(stored_radial)):
            searches = [(query[0] + delta_row, query[1] + delta_column)
                        for delta_row, delta_column
                        in development.footprint.SEARCH_OFFSETS]
            fitted = development.fit_query(
                baseline, query, searches, template, support, planet, width)
            stage.require(np.array_equal(fitted["support"], generic_support),
                          "Stage-I fitted generic support changed")
            calculated = shifted_weights(
                template, shifted, generic_support, fitted)
            weights[local - 1] = calculated
            replay_errors.extend([
                float(np.max(np.abs(
                    fitted["weights"]["raw_rectangular_m0p3"] - stored_raw))),
                float(np.max(np.abs(
                    fitted["weights"]["radial_hann_m0p1_trunc0p75"] -
                    stored_radial))),
            ])
            valid_covariance += 1
        if local % 25 == 0 or local == len(positions):
            print(f"Stage-I shifted calibration {local}/{len(positions)}",
                  flush=True)
    stage.require(valid_covariance > 0 and
                  float(np.max(replay_errors)) <= 3e-10,
                  "Stage-I generic calibration did not replay frozen weights")
    np.savez_compressed(output / "weights.npz", positions=positions,
                        source_indices=source_indices, shifted_weights=weights)
    stage.write_json(output / "results.json", {
        "mode": PRIMARY_MODE,
        "nominal_radius": NOMINAL_RADIUS,
        "phase_offset_row_column": list(phase),
        "positions": len(positions),
        "valid_covariance_positions": valid_covariance,
        "maximum_parent_weight_replay_error": float(np.max(replay_errors)),
    })
    products = [stage.fingerprint(output / "weights.npz"),
                stage.fingerprint(output / "results.json")]
    stage.write_json(complete, {"status": "complete", "products": products})
    verify_receipt(complete)
    with np.load(output / "weights.npz", allow_pickle=False) as source:
        return {name: np.array(source[name], copy=True)
                for name in source.files}


def load_parent_maps(task: dict[str, object]) \
        -> tuple[np.ndarray, np.ndarray, fits.Header]:
    """Load one parent mode-200 amplitude map and its policy map."""
    flat, header = fits.getdata(task["parent_amplitudes"],
                                header=True, memmap=True)
    flat = np.asarray(flat, dtype=np.float64)
    expected = len(stage.MODES) * len(stage_f.PLANET_METHODS)
    stage.require(flat.shape == (expected, 128, 128),
                  "Stage-I parent amplitude schema changed")
    maps = flat.reshape((len(stage.MODES), len(stage_f.PLANET_METHODS),
                         128, 128))[MODE_INDEX]
    policy = np.asarray(fits.getdata(task["parent_policy"], memmap=True),
                        dtype=np.float64)
    stage.require(policy.shape == (len(stage.MODES), 128, 128),
                  "Stage-I parent policy schema changed")
    return maps, policy[MODE_INDEX], header


def candidate_weights(baseline: np.ndarray, parent: dict[str, object],
                      protocol: dict[str, object], task: dict[str, object],
                      responses: np.ndarray, validities: np.ndarray,
                      coordinate_lookup: dict[tuple[int, int], int]) \
        -> tuple[tuple[int, int], np.ndarray, np.ndarray]:
    """Return current and shifted source-safe weights at the optimized nearest pixel."""
    source, _, searches, _ = stage_f.analysis_geometry(
        baseline.shape, dict(task["settings"]))
    query = min(searches, key=lambda value: math.hypot(
        value[0] - source[0], value[1] - source[1]))
    stage.require(query in coordinate_lookup,
                  "Stage-I candidate lacks an exact response")
    displacement = (source[0] - query[0], source[1] - query[1])
    stage.require(np.allclose(displacement,
                              protocol["phase_offset_row_column"],
                              rtol=0, atol=2e-12),
                  "Stage-I candidate phase changed")
    source_index = coordinate_lookup[query]
    full_template = np.asarray(responses[source_index], dtype=np.float64).T
    full_validity = (np.asarray(validities[source_index],
                                dtype=np.float64).T > 0.5)
    template = development.raw.crop(full_template, stage_f.SUPPORT)
    support = development.raw.crop(full_validity, stage_f.SUPPORT)
    shifted, shifted_support = stage_h.shifted_template(
        full_template, full_validity, displacement)
    stage.require(np.array_equal(support, shifted_support) and np.all(support),
                  "Stage-I candidate shifted support changed")
    optimized = stage_g.source_from_polar(
        baseline.shape, float(parent["known_planet"]["separation"]),
        float(parent["known_planet"]["position_angle"]))
    radius_map = development.image_radius(baseline.shape)
    nominal = stage_f.policy_radius(float(radius_map[query[1], query[0]]))
    stage.require(nominal == NOMINAL_RADIUS,
                  "Stage-I candidate left the radius-12 policy")
    width = int(parent["training"]["candidate_specific_half_width_by_radius"][
        str(nominal)])
    fitted = stage_f.fit_planet_query(
        baseline, query, searches, template, support, optimized, width)
    current = np.stack((
        development.normalized_weight(template, support),
        fitted["weights"]["raw_rectangular_m0p3"],
        fitted["weights"]["radial_hann_m0p1_trunc0p75"],
    ))
    shifted_values = shifted_weights(template, shifted, support, fitted)
    return query, current, shifted_values


def production_snr(protocol: dict[str, object], task: dict[str, object],
                   analyzer: Path, maps: np.ndarray, header: fits.Header,
                   directory: Path) -> tuple[np.ndarray, dict[str, object]]:
    """Run production annular SNR and verify it with the independent oracle."""
    working_header = header.copy()
    working_header["HCI FILTER LABELS"] = ",".join(METHODS)
    fits.writeto(directory / "amplitudes.fits",
                 maps.astype(np.float32), working_header)
    known = read(Path(str(protocol["parent_root"])) / "protocol.json")[
        "known_planet"]
    signals = [known]
    if str(task["kind"]) == "injection":
        signals.append({
            "separation": task["source"]["separation"],
            "position_angle": task["source"]["position_angle"],
            "exclusion_radius": task["settings"]["source_radius"],
        })
    command = [
        str(analyzer), "--config",
        str(output_root(Path(str(protocol["parent_root"]))) / "analyze.conf"),
        "--file=amplitudes.fits",
        f"--lambdaD={task['settings']['lambda_d']}",
        "--planet.sep=" + ",".join(str(value["separation"])
                                   for value in signals),
        "--planet.PA=" + ",".join(str(value["position_angle"])
                                  for value in signals),
        "--planet.R=" + ",".join(str(value["exclusion_radius"])
                                 for value in signals),
        f"--snr.apertureR={ANALYZER_REPORTING_APERTURE_RADIUS}",
        f"--snr.minRad={task['settings']['minimum_radius']}",
        f"--snr.maxRad={task['settings']['maximum_radius']}",
        "--filter.psfResponse=", "--filter.lpfGaussFW=0",
        "--filter.hpfGaussFW=0", "--noise.model=identity",
        "--noise.only=false", "--noise.outputDiagnostics=false",
    ]
    with (directory / "analysis.log").open("w", encoding="utf-8") as log:
        subprocess.run(command, cwd=directory,
                       env=development.environment(
                           read(Path(str(protocol["parent_root"])) /
                                "protocol.json"), False),
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    snr, snr_header = fits.getdata(directory / "amplitudes_snr.fits",
                                   header=True)
    snr = np.asarray(snr, dtype=np.float64)
    stage.require(snr.shape == maps.shape and
                  int(snr_header["SNRAPER"]) ==
                  int(ANALYZER_REPORTING_APERTURE_RADIUS) and
                  int(snr_header["SNRMINR"]) ==
                  int(task["settings"]["minimum_radius"]) and
                  int(snr_header["SNRMAXR"]) ==
                  int(task["settings"]["maximum_radius"]) and
                  int(snr_header["SNRMEAN"]) == 1 and
                  int(snr_header["SNRSMALL"]) == 1,
                  "Stage-I hciAnalyze contract changed")
    exclusions = [
        (*stage_g.source_from_polar(
            maps.shape[-2:], float(value["separation"]),
            float(value["position_angle"])),
         float(value["exclusion_radius"]))
        for value in signals
    ]
    radius_map = development.image_radius(maps.shape[-2:])
    maximum = 0.0
    for method_index, method in enumerate(METHODS):
        expected = development.annular_oracle(
            maps[method_index], exclusions)
        expected[(radius_map < float(task["settings"]["minimum_radius"])) |
                 (radius_map > float(task["settings"]["maximum_radius"]))] = np.nan
        valid = np.isfinite(maps[method_index]) & np.isfinite(expected)
        stage.require(np.any(valid) and
                      np.allclose(snr[method_index][valid], expected[valid],
                                  rtol=2e-6, atol=2e-6),
                      f"Stage-I annular oracle mismatch for {method}")
        maximum = max(maximum, float(np.max(np.abs(
            snr[method_index][valid] - expected[valid]))))
    stage.write_json(directory / "command.json", command)
    oracle = {
        "maximum_production_oracle_error": maximum,
        "exclusions_row_column_radius": [list(value)
                                          for value in exclusions],
    }
    stage.write_json(directory / "annular_verification.json", oracle)
    return snr, oracle


def analyze_task(root: Path, protocol: dict[str, object],
                 parent: dict[str, object], analyzer: Path,
                 task: dict[str, object], calibration: dict[str, np.ndarray],
                 coordinate_lookup: dict[tuple[int, int], int],
                 responses: np.ndarray, validities: np.ndarray) -> None:
    """Measure one planet or injection at the optimized nearest pixel."""
    directory = output_root(root) / "analysis" / str(task["name"])
    complete = directory / "complete.json"
    if complete.exists():
        verify_receipt(complete)
        return
    if directory.exists():
        development.archive_incomplete(output_root(root), directory,
                                       "analysis")
    directory.mkdir(parents=True)
    cube, science_header = fits.getdata(task["science"],
                                        header=True, memmap=True)
    stage.require(np.asarray(cube).shape ==
                  (len(stage.MODES), 128, 128) and
                  stage.read_modes(science_header) == stage.MODES,
                  "Stage-I science cube schema changed")
    science = np.asarray(cube[MODE_INDEX], dtype=np.float64)
    baseline_cube = fits.getdata(protocol["paths"]["baseline"], memmap=True)
    stage.require(np.asarray(baseline_cube).shape ==
                  (len(stage.MODES), 128, 128),
                  "Stage-I baseline cube schema changed")
    baseline = np.asarray(baseline_cube[MODE_INDEX], dtype=np.float64)
    parent_maps, policy, _ = load_parent_maps(task)
    maps = np.stack([
        parent_maps[stage_f.PLANET_METHODS.index("gaussian_fwhm2p4")],
        parent_maps[stage_f.PLANET_METHODS.index("gaussian_fwhm3p6")],
        parent_maps[stage_f.PLANET_METHODS.index("exact_identity")],
        parent_maps[stage_f.PLANET_METHODS.index("exact_identity")],
        parent_maps[stage_f.PLANET_METHODS.index("raw_rectangular_m0p3")],
        parent_maps[stage_f.PLANET_METHODS.index("raw_rectangular_m0p3")],
        parent_maps[stage_f.PLANET_METHODS.index(
            "radial_hann_m0p1_trunc0p75")],
        parent_maps[stage_f.PLANET_METHODS.index(
            "radial_hann_m0p1_trunc0p75")],
    ])
    shifted_indices = [METHODS.index(name) for name in SHIFTED_METHODS]
    positions = np.asarray(calibration["positions"], dtype=np.int16)
    weights = np.asarray(calibration["shifted_weights"],
                         dtype=np.float64)
    stage.require(positions.shape == (245, 2) and
                  weights.shape == (245, len(SHIFTED_METHODS),
                                    stage_f.SUPPORT**2),
                  "Stage-I shifted calibration schema changed")
    replaced_mask = np.zeros(science.shape, dtype=bool)
    replaced = 0
    for local, query_value in enumerate(positions):
        row, column = map(int, query_value)
        if not np.isclose(policy[column, row], NOMINAL_RADIUS,
                          rtol=0, atol=1e-12):
            continue
        data = development.stamp(science, (row, column)).ravel()
        for shifted_local, method_index in enumerate(shifted_indices):
            weight = weights[local, shifted_local]
            maps[method_index, column, row] = (
                weight @ data if np.all(np.isfinite(weight)) else np.nan)
        replaced_mask[column, row] = True
        replaced += 1
    stage.require(replaced > 100,
                  "Stage-I shifted generic map lost radius-12 coverage")

    query, current_candidate, shifted_candidate = candidate_weights(
        baseline, parent, protocol, task, responses, validities,
        coordinate_lookup)
    data = development.stamp(science, query).ravel()
    current_indices = [
        METHODS.index("integer_exact_identity"),
        METHODS.index("integer_raw_rectangular_m0p3"),
        METHODS.index("integer_radial_hann_m0p1_trunc0p75"),
    ]
    for local, method_index in enumerate(current_indices):
        maps[method_index, query[1], query[0]] = (
            current_candidate[local] @ data)
    for local, method_index in enumerate(shifted_indices):
        maps[method_index, query[1], query[0]] = (
            shifted_candidate[local] @ data)

    source = (float(task["source"]["row"]), float(task["source"]["column"]))
    query_radius = math.hypot(query[0] - 63.5, query[1] - 63.5)
    bracket = (math.floor(query_radius - 0.5),
               math.floor(query_radius - 0.5) + 1)
    radius_map = development.image_radius(science.shape)
    native_column, native_row = np.indices(science.shape)
    noise = np.ones(science.shape, dtype=bool)
    known = parent["known_planet"]
    exclusions = [(*stage_g.source_from_polar(
        science.shape, float(known["separation"]),
        float(known["position_angle"])), float(known["exclusion_radius"]))]
    if str(task["kind"]) == "injection":
        exclusions.append((*source, float(task["settings"]["source_radius"])))
    for row, column, radius in exclusions:
        noise &= np.hypot(native_row - row, native_column - column) > radius + 0.5
    relevant = noise & (((radius_map > bracket[0]) &
                         (radius_map <= bracket[0] + 1)) |
                        ((radius_map > bracket[1]) &
                         (radius_map <= bracket[1] + 1)))
    relevant &= np.isfinite(
        maps[METHODS.index("integer_exact_identity")])
    stage.require(np.any(relevant) and
                  np.allclose(policy[relevant], NOMINAL_RADIUS,
                              rtol=0, atol=1e-12) and
                  all(np.array_equal(np.isfinite(maps[shifted]),
                                     np.isfinite(maps[current]))
                      for shifted, current in zip(
                          shifted_indices, current_indices)),
                  "Stage-I shifted noise map does not match radius-12 support")

    snr, oracle = production_snr(
        protocol, task, analyzer, maps, science_header, directory)
    methods = {
        method: {
            "amplitude": float(maps[index, query[1], query[0]]),
            "snr": float(snr[index, query[1], query[0]]),
        }
        for index, method in enumerate(METHODS)
    }
    replay = {}
    if str(task["kind"]) == "injection":
        parent_result = read(Path(str(task["parent_results"])))
        parent_mode = parent_result["methods_by_mode"][str(PRIMARY_MODE)]
        for current, parent_name in CURRENT_METHODS.items():
            expected_query = tuple(
                parent_mode[parent_name]["nearest_pixel_row_column"])
            expected = float(parent_mode[parent_name]["nearest_pixel_snr"])
            observed = methods[current]["snr"]
            stage.require(tuple(query) == expected_query and
                          np.isclose(observed, expected,
                                     rtol=2e-6, atol=2e-6),
                          f"Stage-I current replay failed for {current}")
            replay[current] = {
                "parent": expected, "observed": observed,
                "absolute_error": abs(observed - expected),
            }
    result = {
        "schema": 1,
        "task": task["name"],
        "kind": task["kind"],
        "base_site": task["base_site"],
        "source_row_column": list(source),
        "query_row_column": list(query),
        "query_radius": query_radius,
        "phase_offset_row_column": [source[0] - query[0],
                                    source[1] - query[1]],
        "generic_positions_replaced": replaced,
        "noise_bracket_lower_bins": list(bracket),
        "methods": methods,
        "current_parent_replay": replay,
        "annular_oracle_maximum_error":
            oracle["maximum_production_oracle_error"],
    }
    stage.write_json(directory / "results.json", result)
    products = [stage.fingerprint(directory / name) for name in (
        "amplitudes.fits", "amplitudes_snr.fits", "analysis.log",
        "command.json", "annular_verification.json", "results.json")]
    stage.write_json(complete, {"status": "complete", "products": products})
    verify_receipt(complete)


def analyze_all(root: Path, protocol: dict[str, object],
                parent: dict[str, object], analyzer: Path,
                calibration: dict[str, np.ndarray]) -> None:
    """Analyze or verify the planet and all twelve injections."""
    output = output_root(root)
    (output / "analysis").mkdir(exist_ok=True)
    paths = response_paths(parent)
    coordinates = np.asarray(
        fits.getdata(paths["coordinates"], memmap=True),
        dtype=np.float64).T
    lookup = {(int(row), int(column)): index
              for index, (row, column, _, _) in enumerate(coordinates)}
    responses = fits.getdata(paths["responses"][MODE_INDEX], memmap=True)
    validities = fits.getdata(paths["validities"][MODE_INDEX], memmap=True)
    for index, task in enumerate(protocol["tasks"], start=1):
        analyze_task(root, protocol, parent, analyzer, task,
                     calibration, lookup, responses, validities)
        stage.write_json(output / "state.json", {
            "status": "analyzing", "calibration_complete": True,
            "completed_analyses": index,
            "task_count": len(protocol["tasks"]),
        })
        print(f"Stage-I analysis {index}/{len(protocol['tasks'])}: "
              f"{task['name']}", flush=True)


def sample_summary(values: list[float]) -> dict[str, float | int]:
    """Summarize the twelve finite injection values."""
    array = np.asarray(values, dtype=np.float64)
    stage.require(array.size == SITE_COUNT and np.all(np.isfinite(array)),
                  "Stage-I injection summary lost a site")
    return {
        "count": int(array.size),
        "mean": float(np.mean(array)),
        "sample_standard_deviation": float(np.std(array, ddof=1)),
        "minimum": float(np.min(array)),
        "maximum": float(np.max(array)),
        "values": [float(value) for value in array],
    }


def planet_comparison(summary: dict[str, object],
                      observed: float) -> dict[str, float | int]:
    """Compare the planet with an exchangeable injection sample."""
    values = np.asarray(summary["values"], dtype=np.float64)
    deviation = float(summary["sample_standard_deviation"])
    standardized = ((observed - float(summary["mean"])) / deviation
                    if deviation > 0 else math.nan)
    predictive = standardized / math.sqrt(1 + 1 / len(values))
    return {
        "planet": observed,
        "standardized_deviation": standardized,
        "predictive_t": predictive,
        "degrees_of_freedom": len(values) - 1,
        "two_sided_p": float(2 * student_t.sf(
            abs(predictive), df=len(values) - 1)),
        "injections_below_or_equal": int(np.count_nonzero(values <= observed)),
        "injections_above_or_equal": int(np.count_nonzero(values >= observed)),
        "exchangeable_lower_tail_rank":
            float((1 + np.count_nonzero(values <= observed)) /
                  (len(values) + 1)),
        "exchangeable_upper_tail_rank":
            float((1 + np.count_nonzero(values >= observed)) /
                  (len(values) + 1)),
    }


def summarize(root: Path, protocol: dict[str, object]) -> None:
    """Aggregate fixed-pixel planet and injection comparisons."""
    output = output_root(root)
    records = {
        str(task["name"]): read(output / "analysis" /
                                str(task["name"]) / "results.json")
        for task in protocol["tasks"]
    }
    planet = records["planet"]
    injections = [records[str(task["name"])]
                  for task in protocol["tasks"]
                  if str(task["kind"]) == "injection"]
    methods = {}
    for method in METHODS:
        snr = sample_summary([
            float(record["methods"][method]["snr"])
            for record in injections])
        amplitude = sample_summary([
            float(record["methods"][method]["amplitude"])
            for record in injections])
        snr["planet_comparison"] = planet_comparison(
            snr, float(planet["methods"][method]["snr"]))
        amplitude["planet_comparison"] = planet_comparison(
            amplitude, float(planet["methods"][method]["amplitude"]))
        methods[method] = {"snr": snr, "amplitude": amplitude}
    pairs = {}
    for name, first, second in PAIRS:
        values = [
            float(record["methods"][first]["snr"]) -
            float(record["methods"][second]["snr"])
            for record in injections
        ]
        summary = sample_summary(values)
        observed = (float(planet["methods"][first]["snr"]) -
                    float(planet["methods"][second]["snr"]))
        summary["first_method"] = first
        summary["second_method"] = second
        summary["planet_comparison"] = planet_comparison(summary, observed)
        pairs[name] = summary
    result = {
        "schema": 1,
        "purpose": "optimized-nearest-pixel cubic-shifted KLIP response closure",
        "primary_mode": PRIMARY_MODE,
        "phase_offset_row_column": protocol["phase_offset_row_column"],
        "planet": planet,
        "injection_count": len(injections),
        "methods": methods,
        "paired_differences": pairs,
        "maximum_annular_oracle_error": max(
            float(record["annular_oracle_maximum_error"])
            for record in records.values()),
        "maximum_parent_replay_error": max(
            (float(value["absolute_error"])
             for record in injections
             for value in record["current_parent_replay"].values()),
            default=0),
        "limits": [
            "The planet position is the independently optimized negative-companion position.",
            "The statistic is the nearest native pixel at the common optimized phase.",
            "Aperture-maximum closure remains the completed Stage-G endpoint and is not retested.",
        ],
    }
    stage.write_json(output / "results.json", result)

    with (output / "results.csv").open(
            "w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "comparison", "injection_mean", "injection_sample_std",
            "injection_minimum", "injection_maximum", "planet",
            "planet_standardized_deviation", "lower_tail_rank",
            "upper_tail_rank",
        ])
        for name, value in pairs.items():
            comparison = value["planet_comparison"]
            writer.writerow([
                name, value["mean"], value["sample_standard_deviation"],
                value["minimum"], value["maximum"],
                comparison["planet"], comparison["standardized_deviation"],
                comparison["exchangeable_lower_tail_rank"],
                comparison["exchangeable_upper_tail_rank"],
            ])

    lines = [
        "# KLIP Stage-I shifted-template planet closure", "",
        "All values below are mode-200 SNR differences at the nearest native",
        "pixel to the common optimized-phase source position. Injection values",
        "are means +/- sample standard deviations across twelve sites.", "",
        "| Comparison | Injections | Range | Planet | Planet deviation |",
        "| :--- | ---: | :--- | ---: | ---: |",
    ]
    labels = {
        "shifted_identity_minus_integer_identity":
            "Shifted identity - integer identity",
        "shifted_raw_minus_integer_raw":
            "Shifted raw covariance - integer raw covariance",
        "shifted_radial_minus_integer_radial":
            "Shifted radial covariance - integer radial covariance",
        "shifted_raw_minus_shifted_identity":
            "Shifted raw covariance - shifted identity",
        "shifted_radial_minus_shifted_identity":
            "Shifted radial covariance - shifted identity",
        "gaussian_fwhm3p6_minus_shifted_identity":
            "Gaussian 3.6 - shifted identity",
        "gaussian_fwhm3p6_minus_gaussian_fwhm2p4":
            "Gaussian 3.6 - Gaussian 2.4",
    }
    for name, _, _ in PAIRS:
        value = pairs[name]
        comparison = value["planet_comparison"]
        lines.append(
            f"| {labels[name]} | {value['mean']:+.4f} +/- "
            f"{value['sample_standard_deviation']:.4f} | "
            f"[{value['minimum']:+.4f}, {value['maximum']:+.4f}] | "
            f"{comparison['planet']:+.4f} | "
            f"{comparison['standardized_deviation']:+.2f} SD |"
        )
    lines.extend([
        "", "## Absolute fixed-pixel SNR", "",
        "| Method | Injection SNR | Planet SNR | Planet deviation |",
        "| :--- | ---: | ---: | ---: |",
    ])
    for method in METHODS:
        value = methods[method]["snr"]
        comparison = value["planet_comparison"]
        lines.append(
            f"| {method} | {value['mean']:.4f} +/- "
            f"{value['sample_standard_deviation']:.4f} | "
            f"{comparison['planet']:.4f} | "
            f"{comparison['standardized_deviation']:+.2f} SD |"
        )
    lines.extend([
        "", "The report is descriptive. The predeclared paired differences and",
        "exchangeable ranks determine whether response registration resolves the",
        "fixed-pixel planet discrepancy.", "",
    ])
    (output / "results.md").write_text("\n".join(lines), encoding="utf-8")


def finish(root: Path, protocol: dict[str, object]) -> None:
    """Write and verify the final Stage-I receipt."""
    output = output_root(root)
    analysis_receipts = [
        stage.fingerprint(output / "analysis" / str(task["name"]) /
                          "complete.json")
        for task in protocol["tasks"]
    ]
    products = [stage.fingerprint(output / name)
                for name in ("results.json", "results.csv", "results.md")]
    calibration = stage.fingerprint(output / "calibration" / "complete.json")
    stage.write_json(output / "complete.json", {
        "status": "complete",
        "calibration_receipt": calibration,
        "analysis_receipts": analysis_receipts,
        "products": products,
        "method_selection_performed": False,
    })
    stage.write_json(output / "state.json", {
        "status": "complete", "calibration_complete": True,
        "completed_analyses": len(protocol["tasks"]),
        "task_count": len(protocol["tasks"]),
    })
    verify_complete(root)


def verify_complete(root: Path) -> None:
    """Verify the final Stage-I receipt and all nested products."""
    output = output_root(root)
    receipt = read(output / "complete.json")
    stage.require(receipt["status"] == "complete" and
                  not receipt["method_selection_performed"],
                  "Stage-I completion receipt changed")
    stage.verify([receipt["calibration_receipt"], *receipt["analysis_receipts"],
                  *receipt["products"]])
    verify_receipt(Path(str(receipt["calibration_receipt"]["path"])))
    for record in receipt["analysis_receipts"]:
        verify_receipt(Path(str(record["path"])))


def run(root: Path) -> None:
    """Run the resumable shifted calibration and fixed-pixel closure."""
    root = root.resolve()
    protocol, parent, analyzer = enable(root)
    output = output_root(root)
    if (output / "complete.json").exists():
        verify_complete(root)
        print(output / "results.md", flush=True)
        return
    calibration = calibrate(root, protocol, parent)
    analyze_all(root, protocol, parent, analyzer, calibration)
    summarize(root, protocol)
    finish(root, protocol)
    print(output / "results.md", flush=True)


def check() -> None:
    """Check precision reconstruction, pairs, and optimized source geometry."""
    template = np.arange(1, stage_f.SUPPORT**2 + 1,
                         dtype=np.float64).reshape(
                             stage_f.SUPPORT, stage_f.SUPPORT)
    support = np.ones(template.shape, dtype=bool)
    covariance = np.eye(template.size, dtype=np.float64)
    identity = development.normalized_weight(template, support)
    precision = precision_weight(
        covariance, template, support, np.ones(template.shape), "full", 0)
    partial_support = support.copy()
    partial_support[0, :3] = False
    partial_identity = development.normalized_weight(template, partial_support)
    partial_precision = precision_weight(
        covariance, template, partial_support, np.ones(template.shape),
        "full", 0)
    stage.require(np.allclose(identity, precision,
                              rtol=2e-14, atol=2e-14) and
                  np.allclose(partial_identity, partial_precision,
                              rtol=2e-14, atol=2e-14) and
                  len(METHODS) == len(set(METHODS)) and
                  all(first in METHODS and second in METHODS
                      for _, first, second in PAIRS),
                  "Stage-I precision or method contract failed")
    source = stage_g.source_from_polar(
        (128, 128), stage.PLANET_SEPARATION, stage.PLANET_PA)
    stage.require(np.allclose(stage_g.fractional_phase(source),
                              [-0.2770699293779586,
                               0.4859637006126931],
                              rtol=0, atol=1e-12),
                  "Stage-I optimized source phase changed")
    print("KLIP Stage-I shifted-planet checks passed", flush=True)


def parser() -> argparse.ArgumentParser:
    """Build the Stage-I command-line parser."""
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="action", required=True)
    subparsers.add_parser("check")
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("root", type=Path)
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("root", type=Path)
    return result


def main() -> None:
    """Dispatch the requested Stage-I action."""
    arguments = parser().parse_args()
    if arguments.action == "check":
        check()
    elif arguments.action == "prepare":
        prepare(arguments.root)
    else:
        run(arguments.root)


if __name__ == "__main__":
    main()
