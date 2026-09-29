#!/usr/bin/env python3
"""Measure shifted KLIP filters across the complete frozen mode-fraction grid."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
from astropy.io import fits

SCRIPT_PATH = Path(__file__).resolve()
SCRIPT_DIRECTORY = SCRIPT_PATH.parent
FROZEN_DEPENDENCIES = (SCRIPT_PATH.parents[2] / "software"
                       if len(SCRIPT_PATH.parents) > 2 else Path("/nonexistent"))
if (FROZEN_DEPENDENCIES / "run_klip_stage_i_shifted_planet.py").is_file():
    sys.path.insert(0, str(FROZEN_DEPENDENCIES))
else:
    sys.path.insert(0, str(SCRIPT_DIRECTORY))

import run_klip_stage_i_shifted_planet as stage_i  # noqa: E402


stage_h = stage_i.stage_h
stage_g = stage_i.stage_g
stage_f = stage_i.stage_f
development = stage_i.development
stage = stage_i.stage

STAGE_DIRECTORY = "stage_j_mode_dependence"
MODES = tuple(stage.MODES)
MODE_COUNT = len(MODES)
NOMINAL_RADIUS = stage_i.NOMINAL_RADIUS
SITE_COUNT = stage_i.SITE_COUNT
METHODS = stage_i.METHODS
CURRENT_METHODS = stage_i.CURRENT_METHODS
SHIFTED_METHODS = stage_i.SHIFTED_METHODS
PAIRS = stage_i.PAIRS
ANALYZER_REPORTING_APERTURE_RADIUS = stage_i.ANALYZER_REPORTING_APERTURE_RADIUS
PRIMARY_PAIR_NAMES = (
    "shifted_raw_minus_shifted_identity",
    "shifted_radial_minus_shifted_identity",
    "gaussian_fwhm3p6_minus_shifted_identity",
)
DISPLAY_METHODS = (
    "gaussian_fwhm2p4",
    "gaussian_fwhm3p6",
    "shifted_exact_identity",
    "shifted_raw_rectangular_m0p3",
    "shifted_radial_hann_m0p1_trunc0p75",
)


def read(path: Path) -> object:
    """Read one JSON document."""
    return json.loads(path.read_text(encoding="utf-8"))


def output_root(root: Path) -> Path:
    """Return the isolated Stage-J output directory."""
    return root / STAGE_DIRECTORY


def response_paths(parent: dict[str, object]) -> dict[str, object]:
    """Return the promoted exact-response products used by Stage J."""
    return development.response47.product_paths(Path(str(parent["parent_response"])))


def verify_parents(root: Path) \
        -> tuple[dict[str, object], dict[str, object], dict[str, object], Path]:
    """Verify Stage I and return its protocol, root protocol, and analyzer."""
    stage_i.verify_complete(root)
    protocol_i = read(stage_i.output_root(root) / "protocol.json")
    parent = read(root / "protocol.json")
    protocol_g = read(stage_g.stage_g_root(root) / "protocol.json")
    analyzer = Path(str(protocol_i["paths"]["analyzer"]))
    stage.require(analyzer.is_file() and
                  tuple(protocol_i["methods"]) == METHODS and
                  len(protocol_i["tasks"]) == SITE_COUNT + 1 and
                  int(protocol_i["primary_mode"]) == 200 and
                  tuple(stage.read_modes(fits.getheader(
                      protocol_i["tasks"][0]["science"]))) == MODES and
                  str(protocol_g["primary_arm"]) == stage_g.PRIMARY_ARM,
                  "Stage-J parent boundary changed")
    return protocol_i, parent, protocol_g, analyzer


def prepare(root: Path) -> None:
    """Freeze the all-mode shifted-template planet and injection reanalysis."""
    root = root.resolve()
    protocol_i, parent, _, analyzer = verify_parents(root)
    stage.require(sorted(os.sched_getaffinity(0)) == parent["resources"]["cpu_affinity"],
                  "Stage-J prepare must use the frozen response CPU affinity")
    output = output_root(root)
    stage.require(not output.exists(), f"Stage-J output already exists: {output}")
    output.mkdir()
    (output / "software").mkdir()
    shutil.copy2(stage_i.output_root(root) / "analyze.conf", output / "analyze.conf")

    frozen_runner = output / "software" / SCRIPT_PATH.name
    shutil.copy2(SCRIPT_PATH, frozen_runner)
    dependency_names = (
        stage_i.SCRIPT_PATH.name,
        stage_h.SCRIPT_PATH.name,
        stage_g.SCRIPT_PATH.name,
        stage_f.SCRIPT_PATH.name,
    )
    dependency_sources = tuple(
        stage_i.output_root(root) / "software" / name
        for name in dependency_names
    )
    stage.require(all(path.is_file() for path in dependency_sources),
                  "Stage-J frozen dependency set is incomplete")
    for source in dependency_sources:
        shutil.copy2(source, output / "software" / source.name)

    tasks = [dict(task) for task in protocol_i["tasks"]]
    stage.require(tasks[0]["kind"] == "planet" and
                  sum(task["kind"] == "injection" for task in tasks) == SITE_COUNT,
                  "Stage-J task set changed")
    paths = response_paths(parent)
    records = [
        stage.fingerprint(stage_i.output_root(root) / "complete.json"),
        stage.fingerprint(stage_i.output_root(root) / "protocol.json"),
        stage.fingerprint(stage_i.output_root(root) / "results.json"),
        stage.fingerprint(Path(str(parent["paths"]["baseline"]))),
        stage.fingerprint(paths["manifest"]),
        stage.fingerprint(paths["coordinates"]),
        stage.fingerprint(output / "analyze.conf"),
    ]
    unit_positions = {}
    for mode_index, mode in enumerate(MODES):
        unit_directory = (root / "calibration" / "units" /
                          development.unit_name(NOMINAL_RADIUS, mode))
        unit_results = read(unit_directory / "results.json")
        stage.require(int(unit_results["generic_positions"]) == 245 and
                      list(unit_results["required_bins"]) == [10, 11, 12, 13],
                      f"Stage-J radius-12 calibration coverage changed for {mode}")
        unit_positions[str(mode)] = int(unit_results["generic_positions"])
        records.extend([
            stage.fingerprint(unit_directory / "complete.json"),
            stage.fingerprint(unit_directory / "models.npz"),
            stage.fingerprint(unit_directory / "results.json"),
            stage.fingerprint(paths["responses"][mode_index]),
            stage.fingerprint(paths["validities"][mode_index]),
        ])
    for task in tasks:
        for key in ("science", "parent_amplitudes", "parent_policy",
                    "parent_results"):
            records.append(stage.fingerprint(Path(str(task[key]))))

    protocol = {
        "schema": 1,
        "stage": "KLIP Stage J shifted-template mode dependence",
        "purpose": (
            "measure site-specific and planet KLIP SNR curves over the complete "
            "frozen mode-fraction grid without planet-based mode selection"
        ),
        "parent_root": str(root),
        "modes": list(MODES),
        "mode_fraction_by_code": {str(mode): mode / 1000 for mode in MODES},
        "nominal_radius": NOMINAL_RADIUS,
        "phase_offset_row_column": protocol_i["phase_offset_row_column"],
        "methods": list(METHODS),
        "pairs": [list(pair) for pair in PAIRS],
        "primary_pair_names": list(PRIMARY_PAIR_NAMES),
        "tasks": tasks,
        "primary_endpoints": [
            "modewise shifted covariance minus shifted identity SNR",
            "independently mode-scanned method maxima and paired maxima differences",
            "injection-selected modes with leave-one-site-out injection evaluation",
        ],
        "secondary_endpoints": [
            "absolute modewise SNR",
            "within-site SNR standard deviation and range across modes",
            "best-mode histograms across injection sites",
            "shifted-minus-integer response effects by mode",
        ],
        "selection_contract": {
            "planet_selects_mode": False,
            "modewise_results_retained": True,
            "mode_scan_applied_identically_to_planet_and_each_injection": True,
            "injection_selected_planet_mode_uses_all_twelve_injections": True,
            "injection_selected_injection_modes_use_leave_one_site_out": True,
            "method_families_select_modes_independently": True,
            "primary_mode_replaced": False,
        },
        "analysis_contract": {
            "source_position": "optimized negative-companion position",
            "candidate_statistic": "nearest native pixel at common optimized phase",
            "generic_noise_weights": "mode-specific radius-12 signal-free baseline",
            "generic_calibration_positions_by_mode": unit_positions,
            "candidate_weights": "mode-specific complete source-aperture exclusion",
            "annular_normalization": "production hciAnalyze with independent oracle",
            "known_planet_and_trial_excluded_from_noise": True,
            "new_klip_reductions": 0,
        },
        "resources": parent["resources"],
        "paths": {
            "baseline": parent["paths"]["baseline"],
            "parent_response": parent["parent_response"],
            "analyzer": str(analyzer),
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
        "purpose": "immutable Stage-J all-mode shifted-template comparison",
        "root": str(root),
        "input_records": records,
        "software_records": software_records,
        "task_count": len(tasks),
        "mode_count": MODE_COUNT,
        "new_klip_reductions": 0,
        "method_selection_performed": False,
        "planet_score_used_for_mode_selection": False,
    })
    stage.write_json(output / "state.json", {
        "status": "prepared", "calibration_modes_complete": 0,
        "completed_analyses": 0, "task_count": len(tasks),
    })
    print(output / "manifest.json", flush=True)
    print(f"taskset -c 12-27 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "
          f"MKL_NUM_THREADS=1 python3 {frozen_runner} run {root}", flush=True)


def enable(root: Path) \
        -> tuple[dict[str, object], dict[str, object], Path]:
    """Verify the frozen Stage-J boundary."""
    root = root.resolve()
    output = output_root(root)
    protocol = read(output / "protocol.json")
    manifest = read(output / "manifest.json")
    stage.require(tuple(protocol["modes"]) == MODES and
                  tuple(protocol["methods"]) == METHODS and
                  tuple(tuple(pair) for pair in protocol["pairs"]) == PAIRS and
                  len(protocol["tasks"]) == SITE_COUNT + 1 and
                  int(manifest["mode_count"]) == MODE_COUNT and
                  int(manifest["new_klip_reductions"]) == 0 and
                  not manifest["method_selection_performed"] and
                  not manifest["planet_score_used_for_mode_selection"],
                  "Stage-J frozen contract changed")
    stage.verify(manifest["input_records"] + manifest["software_records"])
    stage.require(stage.fingerprint(SCRIPT_PATH) in manifest["software_records"],
                  "Stage-J runner is outside the frozen software receipt")
    stage_i.verify_complete(root)
    parent = read(root / "protocol.json")
    analyzer = Path(str(protocol["paths"]["analyzer"]))
    return protocol, parent, analyzer


def verify_receipt(path: Path) -> None:
    """Verify one completed Stage-J receipt."""
    receipt = read(path)
    stage.require(receipt["status"] == "complete", f"incomplete receipt: {path}")
    stage.verify(receipt["products"])


def calibrate(root: Path, protocol: dict[str, object],
              parent: dict[str, object]) -> dict[str, np.ndarray]:
    """Build phase-shifted radius-12 weights for every frozen KLIP mode."""
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

    baseline_cube = fits.getdata(protocol["paths"]["baseline"], memmap=True)
    stage.require(np.asarray(baseline_cube).shape ==
                  (MODE_COUNT, 128, 128),
                  "Stage-J baseline cube schema changed")
    paths = response_paths(parent)
    planet = stage_g.source_from_polar(
        (128, 128), float(parent["known_planet"]["separation"]),
        float(parent["known_planet"]["position_angle"]))
    phase = tuple(float(value) for value in protocol["phase_offset_row_column"])
    width = int(parent["training"]["candidate_specific_half_width_by_radius"][
        str(NOMINAL_RADIUS)])
    positions_by_mode = np.empty((MODE_COUNT, 245, 2), dtype=np.int16)
    sources_by_mode = np.empty((MODE_COUNT, 245), dtype=np.int32)
    weights = np.full((MODE_COUNT, 245, len(SHIFTED_METHODS),
                       stage_f.SUPPORT**2), np.nan, dtype=np.float64)
    replay_by_mode = {}
    valid_by_mode = {}

    for mode_index, mode in enumerate(MODES):
        unit = development.load_unit(root, NOMINAL_RADIUS, mode)
        positions = np.asarray(unit["positions"], dtype=np.int16)
        source_indices = np.asarray(unit["source_indices"], dtype=np.int32)
        stage.require(positions.shape == (245, 2) and
                      source_indices.shape == (245,) and
                      len({tuple(map(int, value)) for value in positions}) == 245,
                      f"Stage-J calibration schema changed for {mode}")
        positions_by_mode[mode_index] = positions
        sources_by_mode[mode_index] = source_indices
        baseline = np.asarray(baseline_cube[mode_index], dtype=np.float64)
        responses = fits.getdata(paths["responses"][mode_index], memmap=True)
        validities = fits.getdata(paths["validities"][mode_index], memmap=True)
        replay_errors = []
        valid_covariance = 0
        for local, (query_value, source_index) in enumerate(
                zip(positions, source_indices), start=1):
            query = tuple(map(int, query_value))
            full_template = np.asarray(
                responses[int(source_index)], dtype=np.float64).T
            full_validity = (np.asarray(
                validities[int(source_index)], dtype=np.float64).T > 0.5)
            template = development.raw.crop(full_template, stage_f.SUPPORT)
            support = development.raw.crop(full_validity, stage_f.SUPPORT)
            shifted, shifted_support = stage_h.shifted_template(
                full_template, full_validity, phase)
            stage.require(np.array_equal(support, shifted_support) and
                          np.all(support),
                          f"Stage-J shifted support changed for {mode}")
            generic_support = support & development.candidate_mask(query, planet)
            stage.require(generic_support[stage_f.HALF, stage_f.HALF],
                          f"Stage-J generic source anchor excluded for {mode}")
            weights[mode_index, local - 1, 0] = (
                development.normalized_weight(shifted, generic_support))
            stored_identity = unit["reference_weights"][
                local - 1,
                development.REFERENCE_WEIGHT_METHODS.index("exact_identity")]
            current_identity = development.normalized_weight(
                template, generic_support)
            replay_errors.append(float(np.max(np.abs(
                current_identity - stored_identity))))
            stored_raw = unit["covariance_weights"][
                local - 1,
                development.COVARIANCE_METHODS.index(
                    "raw_rectangular_m0p3")]
            stored_radial = unit["covariance_weights"][
                local - 1,
                development.COVARIANCE_METHODS.index(
                    "radial_hann_m0p1_trunc0p75")]
            if np.all(np.isfinite(stored_raw)) and np.all(np.isfinite(stored_radial)):
                searches = [
                    (query[0] + delta_row, query[1] + delta_column)
                    for delta_row, delta_column
                    in development.footprint.SEARCH_OFFSETS
                ]
                fitted = development.fit_query(
                    baseline, query, searches, template, support, planet, width)
                stage.require(np.array_equal(fitted["support"], generic_support),
                              f"Stage-J fitted support changed for {mode}")
                calculated = stage_i.shifted_weights(
                    template, shifted, generic_support, fitted)
                weights[mode_index, local - 1] = calculated
                replay_errors.extend([
                    float(np.max(np.abs(
                        fitted["weights"]["raw_rectangular_m0p3"] - stored_raw))),
                    float(np.max(np.abs(
                        fitted["weights"]["radial_hann_m0p1_trunc0p75"] -
                        stored_radial))),
                ])
                valid_covariance += 1
            if local % 50 == 0 or local == len(positions):
                print(f"Stage-J shifted calibration mode {mode}: "
                      f"{local}/{len(positions)}", flush=True)
        maximum = float(np.max(replay_errors))
        stage.require(valid_covariance > 0 and maximum <= 3e-10,
                      f"Stage-J parent replay failed for {mode}")
        replay_by_mode[str(mode)] = maximum
        valid_by_mode[str(mode)] = valid_covariance
        stage.write_json(output_root(root) / "state.json", {
            "status": "calibrating",
            "calibration_modes_complete": mode_index + 1,
            "completed_analyses": 0,
            "task_count": len(protocol["tasks"]),
        })

    np.savez_compressed(output / "weights.npz",
                        positions=positions_by_mode,
                        source_indices=sources_by_mode,
                        shifted_weights=weights)
    stage.write_json(output / "results.json", {
        "modes": list(MODES),
        "nominal_radius": NOMINAL_RADIUS,
        "phase_offset_row_column": list(phase),
        "positions_per_mode": 245,
        "valid_covariance_positions_by_mode": valid_by_mode,
        "maximum_parent_weight_replay_error_by_mode": replay_by_mode,
        "maximum_parent_weight_replay_error": max(replay_by_mode.values()),
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
    """Load all parent amplitude maps and radius-policy maps."""
    flat, header = fits.getdata(task["parent_amplitudes"],
                                header=True, memmap=True)
    flat = np.asarray(flat, dtype=np.float64)
    expected = MODE_COUNT * len(stage_f.PLANET_METHODS)
    stage.require(flat.shape == (expected, 128, 128),
                  "Stage-J parent amplitude schema changed")
    maps = flat.reshape((MODE_COUNT, len(stage_f.PLANET_METHODS),
                         128, 128))
    policy = np.asarray(fits.getdata(task["parent_policy"], memmap=True),
                        dtype=np.float64)
    stage.require(policy.shape == (MODE_COUNT, 128, 128),
                  "Stage-J parent policy schema changed")
    return maps, policy, header


def candidate_weights(baseline: np.ndarray, parent: dict[str, object],
                      protocol: dict[str, object], task: dict[str, object],
                      responses: np.ndarray, validities: np.ndarray,
                      coordinate_lookup: dict[tuple[int, int], int]) \
        -> tuple[tuple[int, int], np.ndarray, np.ndarray]:
    """Return current and shifted source-safe weights for one mode."""
    source, _, searches, _ = stage_f.analysis_geometry(
        baseline.shape, dict(task["settings"]))
    query = min(searches, key=lambda value: math.hypot(
        value[0] - source[0], value[1] - source[1]))
    stage.require(query in coordinate_lookup,
                  "Stage-J candidate lacks an exact response")
    displacement = (source[0] - query[0], source[1] - query[1])
    stage.require(np.allclose(displacement,
                              protocol["phase_offset_row_column"],
                              rtol=0, atol=2e-12),
                  "Stage-J candidate phase changed")
    source_index = coordinate_lookup[query]
    full_template = np.asarray(responses[source_index], dtype=np.float64).T
    full_validity = (np.asarray(validities[source_index],
                                dtype=np.float64).T > 0.5)
    template = development.raw.crop(full_template, stage_f.SUPPORT)
    support = development.raw.crop(full_validity, stage_f.SUPPORT)
    shifted, shifted_support = stage_h.shifted_template(
        full_template, full_validity, displacement)
    stage.require(np.array_equal(support, shifted_support) and np.all(support),
                  "Stage-J candidate shifted support changed")
    optimized = stage_g.source_from_polar(
        baseline.shape, float(parent["known_planet"]["separation"]),
        float(parent["known_planet"]["position_angle"]))
    radius_map = development.image_radius(baseline.shape)
    nominal = stage_f.policy_radius(float(radius_map[query[1], query[0]]))
    stage.require(nominal == NOMINAL_RADIUS,
                  "Stage-J candidate left the radius-12 policy")
    width = int(parent["training"]["candidate_specific_half_width_by_radius"][
        str(nominal)])
    fitted = stage_f.fit_planet_query(
        baseline, query, searches, template, support, optimized, width)
    current = np.stack((
        development.normalized_weight(template, support),
        fitted["weights"]["raw_rectangular_m0p3"],
        fitted["weights"]["radial_hann_m0p1_trunc0p75"],
    ))
    shifted_values = stage_i.shifted_weights(
        template, shifted, support, fitted)
    return query, current, shifted_values


def production_snr(protocol: dict[str, object], parent: dict[str, object],
                   task: dict[str, object], analyzer: Path,
                   maps: np.ndarray, header: fits.Header,
                   directory: Path) -> tuple[np.ndarray, dict[str, object]]:
    """Run all-mode production annular SNR and verify the independent oracle."""
    flat = maps.reshape((-1, *maps.shape[-2:])).astype(np.float32)
    working_header = header.copy()
    working_header["HCI FILTER LABELS"] = ",".join(
        f"m{mode}_{method}" for mode in MODES for method in METHODS)
    fits.writeto(directory / "amplitudes.fits", flat, working_header)
    known = parent["known_planet"]
    signals = [known]
    if str(task["kind"]) == "injection":
        signals.append({
            "separation": task["source"]["separation"],
            "position_angle": task["source"]["position_angle"],
            "exclusion_radius": task["settings"]["source_radius"],
        })
    command = [
        str(analyzer), "--config", str(output_root(
            Path(str(protocol["parent_root"]))) / "analyze.conf"),
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
                       env=development.environment(parent, False),
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    snr, snr_header = fits.getdata(directory / "amplitudes_snr.fits",
                                   header=True)
    snr = np.asarray(snr, dtype=np.float64).reshape(maps.shape)
    stage.require(snr.shape == maps.shape and
                  int(snr_header["SNRAPER"]) ==
                  int(ANALYZER_REPORTING_APERTURE_RADIUS) and
                  int(snr_header["SNRMINR"]) ==
                  int(task["settings"]["minimum_radius"]) and
                  int(snr_header["SNRMAXR"]) ==
                  int(task["settings"]["maximum_radius"]) and
                  int(snr_header["SNRMEAN"]) == 1 and
                  int(snr_header["SNRSMALL"]) == 1,
                  "Stage-J hciAnalyze contract changed")
    exclusions = [
        (*stage_g.source_from_polar(
            maps.shape[-2:], float(value["separation"]),
            float(value["position_angle"])),
         float(value["exclusion_radius"]))
        for value in signals
    ]
    radius_map = development.image_radius(maps.shape[-2:])
    maximum = 0.0
    checked = 0
    for mode_index, mode in enumerate(MODES):
        for method_index, method in enumerate(METHODS):
            expected = development.annular_oracle(
                maps[mode_index, method_index], exclusions)
            expected[(radius_map < float(task["settings"]["minimum_radius"])) |
                     (radius_map > float(task["settings"]["maximum_radius"]))] = np.nan
            valid = (np.isfinite(maps[mode_index, method_index]) &
                     np.isfinite(expected))
            stage.require(np.any(valid) and
                          np.allclose(snr[mode_index, method_index][valid],
                                      expected[valid], rtol=2e-6, atol=2e-6),
                          f"Stage-J annular oracle mismatch for {mode} {method}")
            maximum = max(maximum, float(np.max(np.abs(
                snr[mode_index, method_index][valid] - expected[valid]))))
            checked += int(np.count_nonzero(valid))
    stage.write_json(directory / "command.json", command)
    oracle = {
        "maximum_production_oracle_error": maximum,
        "checked_pixels": checked,
        "exclusions_row_column_radius": [list(value) for value in exclusions],
    }
    stage.write_json(directory / "annular_verification.json", oracle)
    return snr, oracle


def analyze_task(root: Path, protocol: dict[str, object],
                 parent: dict[str, object], analyzer: Path,
                 task: dict[str, object], calibration: dict[str, np.ndarray],
                 coordinate_lookup: dict[tuple[int, int], int],
                 responses: list[np.ndarray],
                 validities: list[np.ndarray]) -> None:
    """Measure one planet or injection over every KLIP mode."""
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
    stage.require(np.asarray(cube).shape == (MODE_COUNT, 128, 128) and
                  tuple(stage.read_modes(science_header)) == MODES,
                  "Stage-J science cube schema changed")
    science = np.asarray(cube, dtype=np.float64)
    baseline_cube = fits.getdata(protocol["paths"]["baseline"], memmap=True)
    stage.require(np.asarray(baseline_cube).shape ==
                  (MODE_COUNT, 128, 128),
                  "Stage-J baseline cube schema changed")
    parent_maps, policy, _ = load_parent_maps(task)
    maps = np.stack([
        parent_maps[:, stage_f.PLANET_METHODS.index("gaussian_fwhm2p4")],
        parent_maps[:, stage_f.PLANET_METHODS.index("gaussian_fwhm3p6")],
        parent_maps[:, stage_f.PLANET_METHODS.index("exact_identity")],
        parent_maps[:, stage_f.PLANET_METHODS.index("exact_identity")],
        parent_maps[:, stage_f.PLANET_METHODS.index("raw_rectangular_m0p3")],
        parent_maps[:, stage_f.PLANET_METHODS.index("raw_rectangular_m0p3")],
        parent_maps[:, stage_f.PLANET_METHODS.index(
            "radial_hann_m0p1_trunc0p75")],
        parent_maps[:, stage_f.PLANET_METHODS.index(
            "radial_hann_m0p1_trunc0p75")],
    ], axis=1)
    shifted_indices = [METHODS.index(name) for name in SHIFTED_METHODS]
    current_indices = [
        METHODS.index("integer_exact_identity"),
        METHODS.index("integer_raw_rectangular_m0p3"),
        METHODS.index("integer_radial_hann_m0p1_trunc0p75"),
    ]
    source = (float(task["source"]["row"]), float(task["source"]["column"]))
    radius_map = development.image_radius(science.shape[-2:])
    native_column, native_row = np.indices(science.shape[-2:])
    known = parent["known_planet"]
    exclusions = [(*stage_g.source_from_polar(
        science.shape[-2:], float(known["separation"]),
        float(known["position_angle"])), float(known["exclusion_radius"]))]
    if str(task["kind"]) == "injection":
        exclusions.append((*source, float(task["settings"]["source_radius"])))
    noise = np.ones(science.shape[-2:], dtype=bool)
    for row, column, radius in exclusions:
        noise &= np.hypot(native_row - row, native_column - column) > radius + 0.5

    query_by_mode = {}
    replaced_by_mode = {}
    for mode_index, mode in enumerate(MODES):
        positions = np.asarray(calibration["positions"][mode_index],
                               dtype=np.int16)
        weights = np.asarray(calibration["shifted_weights"][mode_index],
                             dtype=np.float64)
        stage.require(positions.shape == (245, 2) and
                      weights.shape == (245, len(SHIFTED_METHODS),
                                        stage_f.SUPPORT**2),
                      f"Stage-J shifted calibration schema changed for {mode}")
        replaced_mask = np.zeros(science.shape[-2:], dtype=bool)
        replaced = 0
        for local, query_value in enumerate(positions):
            row, column = map(int, query_value)
            if not np.isclose(policy[mode_index, column, row], NOMINAL_RADIUS,
                              rtol=0, atol=1e-12):
                continue
            data = development.stamp(
                science[mode_index], (row, column)).ravel()
            for shifted_local, method_index in enumerate(shifted_indices):
                weight = weights[local, shifted_local]
                maps[mode_index, method_index, column, row] = (
                    weight @ data if np.all(np.isfinite(weight)) else np.nan)
            replaced_mask[column, row] = True
            replaced += 1

        query, current_candidate, shifted_candidate = candidate_weights(
            np.asarray(baseline_cube[mode_index], dtype=np.float64),
            parent, protocol, task, responses[mode_index],
            validities[mode_index], coordinate_lookup)
        query_by_mode[str(mode)] = list(query)
        data = development.stamp(science[mode_index], query).ravel()
        for local, method_index in enumerate(current_indices):
            maps[mode_index, method_index, query[1], query[0]] = (
                current_candidate[local] @ data)
        for local, method_index in enumerate(shifted_indices):
            maps[mode_index, method_index, query[1], query[0]] = (
                shifted_candidate[local] @ data)

        query_radius = math.hypot(query[0] - 63.5, query[1] - 63.5)
        bracket = (math.floor(query_radius - 0.5),
                   math.floor(query_radius - 0.5) + 1)
        relevant = noise & (((radius_map > bracket[0]) &
                             (radius_map <= bracket[0] + 1)) |
                            ((radius_map > bracket[1]) &
                             (radius_map <= bracket[1] + 1)))
        relevant &= np.isfinite(
            maps[mode_index, METHODS.index("integer_exact_identity")])
        stage.require(np.any(relevant) and
                      np.all(replaced_mask[relevant]) and
                      np.allclose(policy[mode_index][relevant], NOMINAL_RADIUS,
                                  rtol=0, atol=1e-12) and
                      all(np.array_equal(
                          np.isfinite(maps[mode_index, shifted]),
                          np.isfinite(maps[mode_index, current]))
                          for shifted, current in zip(
                              shifted_indices, current_indices)),
                      f"Stage-J shifted noise support changed for {mode}")
        replaced_by_mode[str(mode)] = replaced

    unique_queries = {tuple(value) for value in query_by_mode.values()}
    stage.require(len(unique_queries) == 1,
                  "Stage-J source query changed across modes")
    query = next(iter(unique_queries))
    snr, oracle = production_snr(
        protocol, parent, task, analyzer, maps, science_header, directory)
    methods_by_mode = {}
    for mode_index, mode in enumerate(MODES):
        methods_by_mode[str(mode)] = {
            method: {
                "amplitude": float(maps[mode_index, method_index,
                                        query[1], query[0]]),
                "snr": float(snr[mode_index, method_index,
                                 query[1], query[0]]),
            }
            for method_index, method in enumerate(METHODS)
        }

    replay = {}
    if str(task["kind"]) == "injection":
        parent_result = read(Path(str(task["parent_results"])))
        for mode in MODES:
            parent_mode = parent_result["methods_by_mode"][str(mode)]
            mode_replay = {}
            for current, parent_name in CURRENT_METHODS.items():
                expected_query = tuple(
                    parent_mode[parent_name]["nearest_pixel_row_column"])
                expected = float(parent_mode[parent_name]["nearest_pixel_snr"])
                observed = methods_by_mode[str(mode)][current]["snr"]
                stage.require(tuple(query) == expected_query and
                              np.isclose(observed, expected,
                                         rtol=2e-6, atol=2e-6),
                              f"Stage-J replay failed for {mode} {current}")
                mode_replay[current] = {
                    "parent": expected,
                    "observed": observed,
                    "absolute_error": abs(observed - expected),
                }
            replay[str(mode)] = mode_replay
    result = {
        "schema": 1,
        "task": task["name"],
        "kind": task["kind"],
        "base_site": task["base_site"],
        "source_row_column": list(source),
        "query_row_column": list(query),
        "phase_offset_row_column": [source[0] - query[0],
                                    source[1] - query[1]],
        "generic_positions_replaced_by_mode": replaced_by_mode,
        "methods_by_mode": methods_by_mode,
        "current_parent_replay_by_mode": replay,
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
    responses = [fits.getdata(path, memmap=True)
                 for path in paths["responses"]]
    validities = [fits.getdata(path, memmap=True)
                  for path in paths["validities"]]
    for index, task in enumerate(protocol["tasks"], start=1):
        analyze_task(root, protocol, parent, analyzer, task,
                     calibration, lookup, responses, validities)
        stage.write_json(output / "state.json", {
            "status": "analyzing",
            "calibration_modes_complete": MODE_COUNT,
            "completed_analyses": index,
            "task_count": len(protocol["tasks"]),
        })
        print(f"Stage-J analysis {index}/{len(protocol['tasks'])}: "
              f"{task['name']}", flush=True)


def sample_summary(values: list[float]) -> dict[str, object]:
    """Summarize the twelve finite injection values."""
    return stage_i.sample_summary(values)


def planet_comparison(summary: dict[str, object], observed: float) \
        -> dict[str, object]:
    """Compare a planet statistic with its exchangeable injection sample."""
    return stage_i.planet_comparison(summary, observed)


def method_curve(record: dict[str, object], method: str) -> np.ndarray:
    """Return one task's ordered all-mode SNR curve."""
    values = np.asarray([
        record["methods_by_mode"][str(mode)][method]["snr"]
        for mode in MODES
    ], dtype=np.float64)
    stage.require(values.shape == (MODE_COUNT,) and np.all(np.isfinite(values)),
                  f"Stage-J lost a finite curve for {method}")
    return values


def scanned_method(record: dict[str, object], method: str) \
        -> tuple[float, int]:
    """Return the maximum SNR and lowest maximizing mode for one task."""
    values = method_curve(record, method)
    index = int(np.argmax(values))
    return float(values[index]), MODES[index]


def leave_one_out_scores(curves: np.ndarray) \
        -> tuple[np.ndarray, np.ndarray]:
    """Select a mode on the other sites and score each held-out site."""
    curves = np.asarray(curves, dtype=np.float64)
    stage.require(curves.shape == (SITE_COUNT, MODE_COUNT) and
                  np.all(np.isfinite(curves)),
                  "Stage-J leave-one-out curve matrix is invalid")
    scores = np.empty(SITE_COUNT, dtype=np.float64)
    modes = np.empty(SITE_COUNT, dtype=np.int32)
    for site in range(SITE_COUNT):
        training = np.delete(curves, site, axis=0)
        selected = int(np.argmax(np.mean(training, axis=0)))
        scores[site] = curves[site, selected]
        modes[site] = MODES[selected]
    return scores, modes


def summarize(root: Path, protocol: dict[str, object]) -> None:
    """Aggregate modewise, scanned, and injection-selected comparisons."""
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
    stage.require(len(injections) == SITE_COUNT,
                  "Stage-J summary lost an injection")

    modewise = {}
    for mode in MODES:
        methods = {}
        for method in METHODS:
            summary = sample_summary([
                float(record["methods_by_mode"][str(mode)][method]["snr"])
                for record in injections])
            summary["planet_comparison"] = planet_comparison(
                summary,
                float(planet["methods_by_mode"][str(mode)][method]["snr"]))
            methods[method] = summary
        pairs = {}
        for name, first, second in PAIRS:
            values = [
                float(record["methods_by_mode"][str(mode)][first]["snr"]) -
                float(record["methods_by_mode"][str(mode)][second]["snr"])
                for record in injections
            ]
            summary = sample_summary(values)
            observed = (
                float(planet["methods_by_mode"][str(mode)][first]["snr"]) -
                float(planet["methods_by_mode"][str(mode)][second]["snr"])
            )
            summary["first_method"] = first
            summary["second_method"] = second
            summary["planet_comparison"] = planet_comparison(summary, observed)
            pairs[name] = summary
        modewise[str(mode)] = {"methods": methods,
                               "paired_differences": pairs}

    scanned_methods = {}
    scanned_values = {}
    for method in METHODS:
        injection_scans = [scanned_method(record, method)
                           for record in injections]
        planet_scan = scanned_method(planet, method)
        values = [value for value, _ in injection_scans]
        summary = sample_summary(values)
        summary["planet_comparison"] = planet_comparison(
            summary, planet_scan[0])
        summary["planet_selected_mode"] = planet_scan[1]
        summary["injection_selected_modes"] = [mode for _, mode in injection_scans]
        summary["injection_mode_histogram"] = {
            str(mode): count
            for mode, count in sorted(Counter(
                mode for _, mode in injection_scans).items())
        }
        scanned_methods[method] = summary
        scanned_values[method] = {
            "injections": np.asarray(values, dtype=np.float64),
            "planet": planet_scan[0],
        }
    scanned_pairs = {}
    for name, first, second in PAIRS:
        values = (scanned_values[first]["injections"] -
                  scanned_values[second]["injections"])
        summary = sample_summary(values.tolist())
        observed = (float(scanned_values[first]["planet"]) -
                    float(scanned_values[second]["planet"]))
        summary["first_method"] = first
        summary["second_method"] = second
        summary["planet_comparison"] = planet_comparison(summary, observed)
        scanned_pairs[name] = summary

    selected_methods = {}
    selected_values = {}
    for method in METHODS:
        curves = np.stack([method_curve(record, method)
                           for record in injections])
        scores, modes = leave_one_out_scores(curves)
        training_means = np.mean(curves, axis=0)
        planet_mode_index = int(np.argmax(training_means))
        planet_mode = MODES[planet_mode_index]
        planet_score = float(method_curve(planet, method)[planet_mode_index])
        summary = sample_summary(scores.tolist())
        summary["planet_comparison"] = planet_comparison(summary, planet_score)
        summary["planet_mode_selected_from_all_injections"] = planet_mode
        summary["all_injection_mean_by_mode"] = {
            str(mode): float(value)
            for mode, value in zip(MODES, training_means)
        }
        summary["leave_one_out_selected_modes"] = [int(value) for value in modes]
        summary["leave_one_out_mode_histogram"] = {
            str(mode): count
            for mode, count in sorted(Counter(map(int, modes)).items())
        }
        selected_methods[method] = summary
        selected_values[method] = {
            "injections": scores,
            "planet": planet_score,
        }
    selected_pairs = {}
    for name, first, second in PAIRS:
        values = (selected_values[first]["injections"] -
                  selected_values[second]["injections"])
        summary = sample_summary(values.tolist())
        observed = (float(selected_values[first]["planet"]) -
                    float(selected_values[second]["planet"]))
        summary["first_method"] = first
        summary["second_method"] = second
        summary["planet_comparison"] = planet_comparison(summary, observed)
        selected_pairs[name] = summary

    variability = {}
    for method in METHODS:
        injection_curves = [method_curve(record, method)
                            for record in injections]
        standard_deviations = [float(np.std(curve, ddof=1))
                               for curve in injection_curves]
        ranges = [float(np.max(curve) - np.min(curve))
                  for curve in injection_curves]
        planet_curve = method_curve(planet, method)
        deviation_summary = sample_summary(standard_deviations)
        range_summary = sample_summary(ranges)
        deviation_summary["planet_comparison"] = planet_comparison(
            deviation_summary, float(np.std(planet_curve, ddof=1)))
        range_summary["planet_comparison"] = planet_comparison(
            range_summary, float(np.max(planet_curve) - np.min(planet_curve)))
        variability[method] = {
            "within_site_sample_standard_deviation": deviation_summary,
            "within_site_range": range_summary,
        }

    replay_errors = [
        float(value["absolute_error"])
        for record in injections
        for mode in record["current_parent_replay_by_mode"].values()
        for value in mode.values()
    ]
    result = {
        "schema": 1,
        "purpose": "all-mode shifted-template KLIP planet consistency",
        "modes": list(MODES),
        "mode_fraction_by_code": protocol["mode_fraction_by_code"],
        "injection_count": len(injections),
        "planet": planet,
        "modewise": modewise,
        "independent_mode_scan": {
            "methods": scanned_methods,
            "paired_differences": scanned_pairs,
        },
        "injection_selected_mode": {
            "methods": selected_methods,
            "paired_differences": selected_pairs,
        },
        "mode_variability": variability,
        "maximum_annular_oracle_error": max(
            float(record["annular_oracle_maximum_error"])
            for record in records.values()),
        "maximum_parent_replay_error": max(replay_errors, default=0),
        "limits": [
            "All modewise curves are retained; no mode is chosen from the planet for inference.",
            "Independent scans apply the same maximum operation to every planet and injection curve.",
            "Injection-selected planet modes use all twelve injections; "
            "injection scores use leave-one-site-out selection.",
            "The twelve injections share an observing sequence and correlated residual field.",
        ],
    }
    stage.write_json(output / "results.json", result)

    with (output / "mode_curves.csv").open(
            "w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["task", "kind", "base_site", "mode_code",
                         "mode_fraction", "method", "snr", "amplitude"])
        for task in protocol["tasks"]:
            record = records[str(task["name"])]
            for mode in MODES:
                for method in METHODS:
                    value = record["methods_by_mode"][str(mode)][method]
                    writer.writerow([
                        task["name"], task["kind"], task["base_site"], mode,
                        mode / 1000, method, value["snr"], value["amplitude"],
                    ])

    with (output / "mode_summary.csv").open(
            "w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["mode_code", "mode_fraction", "method",
                         "injection_mean", "injection_sample_std", "planet",
                         "planet_standardized_deviation", "lower_tail_rank",
                         "upper_tail_rank"])
        for mode in MODES:
            for method in METHODS:
                value = modewise[str(mode)]["methods"][method]
                comparison = value["planet_comparison"]
                writer.writerow([
                    mode, mode / 1000, method, value["mean"],
                    value["sample_standard_deviation"], comparison["planet"],
                    comparison["standardized_deviation"],
                    comparison["exchangeable_lower_tail_rank"],
                    comparison["exchangeable_upper_tail_rank"],
                ])

    lines = [
        "# KLIP Stage-J mode-dependence result", "",
        "Values are optimized-nearest-pixel SNR at the common source phase.",
        "Injection entries are means +/- sample standard deviations over twelve sites.",
        "", "## Modewise absolute SNR", "",
        "| Mode fraction | Method | Injections | Planet | Planet deviation |",
        "| :--- | :--- | ---: | ---: | ---: |",
    ]
    for mode in MODES:
        for method in DISPLAY_METHODS:
            value = modewise[str(mode)]["methods"][method]
            comparison = value["planet_comparison"]
            lines.append(
                f"| {mode / 1000:.3f} | {method} | "
                f"{value['mean']:.4f} +/- "
                f"{value['sample_standard_deviation']:.4f} | "
                f"{comparison['planet']:.4f} | "
                f"{comparison['standardized_deviation']:+.2f} SD |"
            )
    lines.extend([
        "", "## Modewise primary paired differences", "",
        "| Mode fraction | Comparison | Injections | Planet | Planet deviation |",
        "| :--- | :--- | ---: | ---: | ---: |",
    ])
    for mode in MODES:
        for name in PRIMARY_PAIR_NAMES:
            value = modewise[str(mode)]["paired_differences"][name]
            comparison = value["planet_comparison"]
            lines.append(
                f"| {mode / 1000:.3f} | {name} | "
                f"{value['mean']:+.4f} +/- "
                f"{value['sample_standard_deviation']:.4f} | "
                f"{comparison['planet']:+.4f} | "
                f"{comparison['standardized_deviation']:+.2f} SD |"
            )
    lines.extend([
        "", "## Independent maximum over modes", "",
        "Each task and method selects its own maximum before planet/injection comparison.",
        "", "| Method | Injections | Planet | Planet mode | Planet deviation |",
        "| :--- | ---: | ---: | ---: | ---: |",
    ])
    for method in DISPLAY_METHODS:
        value = scanned_methods[method]
        comparison = value["planet_comparison"]
        lines.append(
            f"| {method} | {value['mean']:.4f} +/- "
            f"{value['sample_standard_deviation']:.4f} | "
            f"{comparison['planet']:.4f} | "
            f"{value['planet_selected_mode'] / 1000:.3f} | "
            f"{comparison['standardized_deviation']:+.2f} SD |"
        )
    lines.extend([
        "", "## Injection-selected modes", "",
        "The planet mode maximizes the twelve-injection mean. Each injection is",
        "scored at the mode selected by the other eleven sites.", "",
        "| Method | Injections, leave-one-out | Planet | Planet mode | Planet deviation |",
        "| :--- | ---: | ---: | ---: | ---: |",
    ])
    for method in DISPLAY_METHODS:
        value = selected_methods[method]
        comparison = value["planet_comparison"]
        lines.append(
            f"| {method} | {value['mean']:.4f} +/- "
            f"{value['sample_standard_deviation']:.4f} | "
            f"{comparison['planet']:.4f} | "
            f"{value['planet_mode_selected_from_all_injections'] / 1000:.3f} | "
            f"{comparison['standardized_deviation']:+.2f} SD |"
        )
    lines.extend([
        "", "The complete per-site curves, paired scans, mode histograms, and",
        "within-site mode variability are retained in results.json and mode_curves.csv.", "",
    ])
    (output / "results.md").write_text("\n".join(lines), encoding="utf-8")


def finish(root: Path, protocol: dict[str, object]) -> None:
    """Write and verify the final Stage-J receipt."""
    output = output_root(root)
    analysis_receipts = [
        stage.fingerprint(output / "analysis" / str(task["name"]) /
                          "complete.json")
        for task in protocol["tasks"]
    ]
    products = [stage.fingerprint(output / name) for name in (
        "results.json", "mode_curves.csv", "mode_summary.csv", "results.md")]
    calibration = stage.fingerprint(output / "calibration" / "complete.json")
    stage.write_json(output / "complete.json", {
        "status": "complete",
        "calibration_receipt": calibration,
        "analysis_receipts": analysis_receipts,
        "products": products,
        "method_selection_performed": False,
        "planet_score_used_for_mode_selection": False,
    })
    stage.write_json(output / "state.json", {
        "status": "complete", "calibration_modes_complete": MODE_COUNT,
        "completed_analyses": len(protocol["tasks"]),
        "task_count": len(protocol["tasks"]),
    })
    verify_complete(root)


def verify_complete(root: Path) -> None:
    """Verify the final Stage-J receipt and all nested products."""
    output = output_root(root)
    receipt = read(output / "complete.json")
    stage.require(receipt["status"] == "complete" and
                  not receipt["method_selection_performed"] and
                  not receipt["planet_score_used_for_mode_selection"],
                  "Stage-J completion receipt changed")
    stage.verify([receipt["calibration_receipt"],
                  *receipt["analysis_receipts"], *receipt["products"]])
    verify_receipt(Path(str(receipt["calibration_receipt"]["path"])))
    for record in receipt["analysis_receipts"]:
        verify_receipt(Path(str(record["path"])))


def run(root: Path) -> None:
    """Run the resumable all-mode shifted-template comparison."""
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
    """Check the mode grid, scan, and leave-one-out selection contracts."""
    stage.require(MODES == (125, 150, 175, 200, 225, 250, 300, 350) and
                  len(METHODS) == len(set(METHODS)) and
                  all(first in METHODS and second in METHODS
                      for _, first, second in PAIRS),
                  "Stage-J mode or method contract changed")
    curves = np.tile(np.arange(MODE_COUNT, dtype=np.float64),
                     (SITE_COUNT, 1))
    curves[:, 3] += 20
    scores, modes = leave_one_out_scores(curves)
    stage.require(np.all(modes == MODES[3]) and
                  np.all(scores == curves[:, 3]),
                  "Stage-J leave-one-out selection failed")
    synthetic = {
        "methods_by_mode": {
            str(mode): {method: {"snr": float(index), "amplitude": 0.0}
                        for method in METHODS}
            for index, mode in enumerate(MODES)
        }
    }
    value, mode = scanned_method(synthetic, METHODS[0])
    stage.require(value == MODE_COUNT - 1 and mode == MODES[-1],
                  "Stage-J independent scan failed")
    print("KLIP Stage-J mode-dependence checks passed", flush=True)


def parser() -> argparse.ArgumentParser:
    """Build the Stage-J command-line parser."""
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="action", required=True)
    subparsers.add_parser("check")
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("root", type=Path)
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("root", type=Path)
    return result


def main() -> None:
    """Dispatch the requested Stage-J action."""
    arguments = parser().parse_args()
    if arguments.action == "check":
        check()
    elif arguments.action == "prepare":
        prepare(arguments.root)
    else:
        run(arguments.root)


if __name__ == "__main__":
    main()
