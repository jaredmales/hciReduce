#!/usr/bin/env python3
"""Measure the known KLIP planet after immutable Stage-E validation."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter

SCRIPT_PATH = Path(__file__).resolve()
SCRIPT_DIRECTORY = SCRIPT_PATH.parent
FROZEN_DEPENDENCIES = (SCRIPT_PATH.parents[2] / "software"
                       if len(SCRIPT_PATH.parents) > 2 else Path("/nonexistent"))
if (FROZEN_DEPENDENCIES / "run_klip_stage_e_validation.py").is_file():
    sys.path.insert(0, str(FROZEN_DEPENDENCIES))
else:
    sys.path.insert(0, str(SCRIPT_DIRECTORY))

import freeze_klip_stage_d_policy as policy_module  # noqa: E402
import prepare_klip_stage_c_development as preparation  # noqa: E402
import run_klip_covariance_stage_a as stage  # noqa: E402
import run_klip_stage_c_development as development  # noqa: E402
import run_klip_stage_e_validation as validation  # noqa: E402


PLANET_METHODS = (
    "native",
    "gaussian_fwhm2p4",
    "gaussian_fwhm3p6",
    "exact_identity",
    "sparse_identity",
    "exact_identity_lpf1p8",
    "exact_identity_lpf2p7",
    "exact_isotropic_fitted_mean",
    "raw_rectangular_m0p3",
    "radial_hann_m0p1_trunc0p75",
)
METHOD_LABELS = {
    "native": "Native",
    "gaussian_fwhm2p4": "Gaussian 2.4",
    "gaussian_fwhm3p6": "Gaussian 3.6",
    "exact_identity": "Exact identity",
    "sparse_identity": "Sparse identity",
    "exact_identity_lpf1p8": "Exact response LPF 1.8",
    "exact_identity_lpf2p7": "Exact response LPF 2.7",
    "exact_isotropic_fitted_mean": "Exact identity, fitted mean",
    "raw_rectangular_m0p3": "Raw rectangular PSD",
    "radial_hann_m0p1_trunc0p75": "Radial Hann, truncation 0.75",
}
SUPPORT = development.SUPPORT
HALF = development.HALF
EXPECTED_SETTINGS = {
    "lambda_d": 3.6,
    "separation": 11.782,
    "position_angle": 262.051,
    "source_radius": 7.0,
    "minimum_radius": 6.0,
    "maximum_radius": 60.0,
    "aperture_radius": 3.0,
}


def read(path: Path) -> object:
    """Read one JSON document."""
    return json.loads(path.read_text(encoding="utf-8"))


def stage_f_root(root: Path) -> Path:
    """Return the isolated Stage-F product directory."""
    return root / "stage_f_planet"


def parse_analysis_config(path: Path) -> dict[str, float]:
    """Read the small hciAnalyze guidance file with its sectionless root key."""
    values: dict[str, str] = {}
    section = ""
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1].strip().lower()
            continue
        stage.require("=" in line, f"invalid analysis configuration line: {raw_line}")
        key, value = (part.strip() for part in line.split("=", 1))
        values[f"{section}.{key.lower()}" if section else key.lower()] = value
    required = ("lambdad", "planet.sep", "planet.pa", "planet.r",
                "snr.minrad", "snr.aperturer")
    stage.require(all(key in values for key in required),
                  "working/analyze.conf guidance is incomplete")
    result = {
        "lambda_d": float(values["lambdad"]),
        "separation": float(values["planet.sep"]),
        "position_angle": float(values["planet.pa"]),
        "source_radius": float(values["planet.r"]),
        "minimum_radius": float(values["snr.minrad"]),
        "maximum_radius": float(values.get("snr.maxrad", "60")),
        "aperture_radius": float(values["snr.aperturer"]),
    }
    stage.require(all(np.isclose(result[key], expected, rtol=0, atol=1e-12)
                      for key, expected in EXPECTED_SETTINGS.items()),
                  "working/analyze.conf guidance changed from the reviewed Stage-F contract")
    return result


def original_science_path(protocol: dict[str, object]) -> tuple[Path, Path, Path]:
    """Follow the promoted response provenance to the original science cube."""
    response_protocol_path = Path(str(protocol["parent_response"])) / "protocol.json"
    response_protocol = read(response_protocol_path)
    stage_a_root = Path(str(response_protocol["parent_stage_a"]))
    stage_a_protocol_path = stage_a_root / "protocol.json"
    stage_a_protocol = read(stage_a_protocol_path)
    science = Path(str(stage_a_protocol["paths"]["original_science"]))
    return science, response_protocol_path, stage_a_protocol_path


def verify_stage_e(root: Path) -> tuple[dict[str, object], dict[str, object], Path]:
    """Recursively verify completed Stage E before exposing the planet."""
    protocol, _ = preparation.load_experiment(root)
    policy_completion = read(root / "policy_complete.json")
    stage.require(policy_completion["status"] == "complete",
                  "Stage-D policy is incomplete")
    stage.verify(policy_completion["products"])
    policy_manifest = read(root / "policy_manifest.json")
    stage.verify(policy_manifest["input_records"] + policy_manifest["software_records"] +
                 policy_manifest["policy_records"])
    validation.verify_repair_records(policy_manifest.get("repair_records", []))
    local_validation_hash = stage.fingerprint(Path(validation.__file__))["sha256"]
    stage.require(local_validation_hash in {record["sha256"]
                                            for record in policy_manifest["software_records"]},
                  "Stage-F imported a validation implementation outside the frozen policy")
    policy = read(root / "policy" / "policy.json")
    stage.require(policy["frozen_before_heldout_or_validation"] and
                  policy["primary_mode"] == policy_module.PRIMARY_MODE and
                  tuple(policy["methods"]["ordered"]) == validation.METHODS and
                  tuple(policy["methods"]["covariance_candidates"]) ==
                  validation.COVARIANCE_CANDIDATES and not policy["joint_mode_maximization"],
                  "frozen Stage-E policy changed")
    analyzer = Path(str(read(root / "development_manifest.json")["hcianalyze"]["path"]))
    stage.require(analyzer.is_file(), "frozen hciAnalyze executable is missing")
    validation.verify_models(root)
    validation.verify_heldout(root)
    receipt = read(root / "validation_complete.json")
    stage.require(receipt["status"] == "complete" and receipt["policy_unchanged"] and
                  not receipt["known_planet_opened"],
                  "Stage-E completion receipt is not the unopened frozen endpoint")
    stage.verify(receipt["products"])
    for task in protocol["validation_tasks_unopened"]:
        name = str(task["name"])
        for parent in ("validation_reductions", "validation_analysis"):
            nested = read(root / parent / name / "complete.json")
            stage.verify(nested["products"])
    state = read(root / "state.json")
    stage.require(state["status"] == "validation_complete" and
                  not state["known_planet_opened"],
                  "Stage-E state changed before Stage F")
    return protocol, policy, analyzer


def model_receipts(root: Path) -> list[dict[str, object]]:
    """Fingerprint the frozen calibration units used by Stage F."""
    records = []
    for radius in preparation.RADII:
        for mode in stage.MODES:
            path = (root / "calibration" / "units" /
                    development.unit_name(radius, mode) / "complete.json")
            stage.verify(read(path)["products"])
            records.append(stage.fingerprint(path))
    return records


def prepare(root: Path, config: Path) -> None:
    """Freeze Stage-F inputs and software after validating the Stage-E boundary."""
    root = root.resolve()
    config = config.resolve()
    protocol, policy, analyzer = verify_stage_e(root)
    settings = parse_analysis_config(config)
    output = stage_f_root(root)
    stage.require(not output.exists(), f"Stage-F output already exists: {output}")
    science, response_protocol, stage_a_protocol = original_science_path(protocol)
    for path in (science, response_protocol, stage_a_protocol, analyzer):
        stage.require(path.is_file(), f"missing Stage-F input: {path}")
    science_record = stage.fingerprint(science)
    stage_a_manifest_path = stage_a_protocol.parent / "manifest.json"
    stage_a_manifest = read(stage_a_manifest_path)
    expected = [record for record in stage_a_manifest["fixed_records"]
                if Path(str(record["path"])) == science]
    stage.require(len(expected) == 1 and expected[0] == science_record,
                  "original science cube differs from the Stage-A frozen inventory")
    cube, header = fits.getdata(science, header=True, memmap=True)
    stage.require(tuple(cube.shape) == (len(stage.MODES), 128, 128) and
                  stage.read_modes(header) == stage.MODES,
                  "original science cube mode ordering or shape changed")

    output.mkdir()
    (output / "software").mkdir()
    frozen_runner = output / "software" / SCRIPT_PATH.name
    shutil.copy2(SCRIPT_PATH, frozen_runner)
    frozen_config = output / "analyze.conf"
    shutil.copy2(config, frozen_config)
    inputs = [
        science_record,
        stage.fingerprint(frozen_config),
        stage.fingerprint(root / "protocol.json"),
        stage.fingerprint(root / "policy_complete.json"),
        stage.fingerprint(root / "policy_manifest.json"),
        stage.fingerprint(root / "validation_complete.json"),
        stage.fingerprint(root / "validation_results.json"),
        stage.fingerprint(root / "heldout_results.json"),
        stage.fingerprint(response_protocol),
        stage.fingerprint(stage_a_protocol),
        stage.fingerprint(stage_a_manifest_path),
    ]
    software = [stage.fingerprint(frozen_runner), stage.fingerprint(analyzer)]
    manifest = {
        "schema": 1,
        "purpose": "descriptive known-planet endpoint after immutable KLIP Stage-E validation",
        "stage_e_policy_changed": False,
        "method_selection_performed": False,
        "known_planet_opened": False,
        "root": str(root),
        "settings": settings,
        "methods": list(PLANET_METHODS),
        "validated_methods": list(policy["methods"]["ordered"]),
        "additional_predeclared_controls": [method for method in PLANET_METHODS
                                               if method not in policy["methods"]["ordered"]],
        "radial_policy": "nearest frozen Stage-C nominal radius; ties use the smaller radius",
        "training_exclusion": "optimized-planet disk plus the union of all 3-pixel analysis-aperture 11x11 footprints",
        "candidate_support": "exact or sparse response validity without the training planet mask",
        "source_config": stage.fingerprint(config),
        "input_records": inputs,
        "software_records": software,
        "model_receipts": model_receipts(root),
    }
    stage.write_json(output / "manifest.json", manifest)
    stage.write_json(output / "state.json", {
        "status": "prepared", "known_planet_opened": False,
        "stage_e_policy_changed": False, "method_selection_performed": False,
    })
    print(output / "manifest.json", flush=True)
    print(f"python3 {frozen_runner} run {root}", flush=True)


def enable(root: Path) -> tuple[dict[str, object], Path, dict[str, float], Path]:
    """Verify the frozen Stage-F manifest and return its required inputs."""
    root = root.resolve()
    protocol, _, analyzer = verify_stage_e(root)
    output = stage_f_root(root)
    manifest = read(output / "manifest.json")
    stage.require(not manifest["known_planet_opened"] and
                  not manifest["stage_e_policy_changed"] and
                  not manifest["method_selection_performed"] and
                  tuple(manifest["methods"]) == PLANET_METHODS,
                  "Stage-F frozen contract changed")
    stage.verify(manifest["input_records"] + manifest["software_records"] +
                 manifest["model_receipts"])
    stage.require(model_receipts(root) == manifest["model_receipts"],
                  "Stage-F calibration model products changed")
    stage.require(stage.fingerprint(SCRIPT_PATH) in manifest["software_records"],
                  "Stage-F runner is outside the frozen software receipt")
    settings = parse_analysis_config(output / "analyze.conf")
    stage.require(settings == manifest["settings"], "Stage-F analysis settings changed")
    science, _, _ = original_science_path(protocol)
    return protocol, analyzer, settings, science


def policy_radius(value: float) -> float:
    """Map an arbitrary candidate radius to the nearest frozen policy radius."""
    return min((float(radius) for radius in preparation.RADII),
               key=lambda radius: (abs(radius - value), radius))


def planet_reference_weights(template: np.ndarray, validity: np.ndarray,
                             sparse_template: np.ndarray, sparse_validity: np.ndarray) \
        -> tuple[dict[str, np.ndarray], dict[str, float], dict[str, int]]:
    """Build predeclared identity controls without masking the source itself."""
    exact_support = np.asarray(validity, dtype=bool)
    sparse_support = np.asarray(sparse_validity, dtype=bool)
    stage.require(exact_support[HALF, HALF] and sparse_support[HALF, HALF],
                  "planet response anchor is invalid")
    templates = {
        "exact_identity": np.asarray(template, dtype=np.float64),
        "sparse_identity": np.asarray(sparse_template, dtype=np.float64),
    }
    for name, fwhm in (("exact_identity_lpf1p8", 1.8),
                       ("exact_identity_lpf2p7", 2.7)):
        sigma = fwhm / math.sqrt(8 * math.log(2))
        templates[name] = gaussian_filter(template, sigma=sigma, mode="constant",
                                          cval=0.0, truncate=4.0)
    supports = {
        "exact_identity": exact_support,
        "sparse_identity": sparse_support,
        "exact_identity_lpf1p8": exact_support,
        "exact_identity_lpf2p7": exact_support,
    }
    weights = {name: development.normalized_weight(value, supports[name])
               for name, value in templates.items()}
    responses = {name: float(weight @ template.ravel())
                 for name, weight in weights.items()}
    counts = {name: int(np.count_nonzero(supports[name])) for name in weights}
    return weights, responses, counts


def fit_planet_query(image: np.ndarray, query: tuple[int, int],
                     searches: list[tuple[int, int]], template: np.ndarray,
                     validity: np.ndarray, optimized_planet: tuple[float, float],
                     width: int) -> dict[str, object]:
    """Fit frozen PSD models while keeping the source pixels in candidate support."""
    finite = np.isfinite(image)
    rings = development.raw.training_stencils(finite, query, searches,
                                               optimized_planet, SUPPORT)
    raw_samples, raw_counts = development.combined_samples(image, rings, width)
    rectangular = development.raw.regularize(
        development.raw.fit_periodogram_base(raw_samples, SUPPORT, "rectangular"), 0.3)

    excluded = development.radial.exclusion_mask(image.shape, searches, SUPPORT,
                                                  optimized_planet)
    counts = development.radial.profile_counts(image, excluded)
    stage.require(np.all(counts >= development.radial.MINIMUM_PROFILE_PIXELS),
                  "Stage-F radial profile lost frozen support")
    scale_map, profile = development.radial.variance_profile(image, excluded)
    standardized = image / scale_map
    radial_samples, radial_counts = development.combined_samples(standardized, rings, width)
    hann = development.raw.regularize(
        development.raw.fit_periodogram_base(radial_samples, SUPPORT, "hann"), 0.1)

    support = np.asarray(validity, dtype=bool)
    data = development.stamp(image, query)
    scale = development.stamp(scale_map, query)
    stage.require(support[HALF, HALF] and np.all(np.isfinite(data[support])) and
                  np.all(np.isfinite(scale[support])) and np.all(scale[support] > 0),
                  "Stage-F candidate support is invalid")
    raw_weights, _, raw_detail = development.precision_grid(
        rectangular, template, support, np.ones(template.shape),
        (("raw_rectangular_m0p3", "full", 0.0),))
    radial_weights, _, radial_detail = development.precision_grid(
        hann, template, support, scale, development.RADIAL_POLICIES)
    for record in raw_detail["policies"].values():
        record["samples"] = len(raw_samples)
        record["split_samples"] = raw_counts
    for record in radial_detail["policies"].values():
        record["samples"] = len(radial_samples)
        record["split_samples"] = radial_counts
    return {
        "weights": raw_weights | radial_weights,
        "isotropic_mean": np.asarray(rectangular["mean"], dtype=np.float64),
        "support": support,
        "scale": scale,
        "raw_detail": raw_detail,
        "radial_detail": radial_detail,
        "profile_minimum_pixels": int(np.min(counts)),
        "profile": profile,
    }


def analysis_geometry(shape: tuple[int, int], settings: dict[str, float]) \
        -> tuple[tuple[float, float], np.ndarray, list[tuple[int, int]], list[int]]:
    """Build the reviewed source aperture and its required one-pixel radial bins."""
    metadata = {"known_planet": {"separation": settings["separation"],
                                  "position_angle": settings["position_angle"]}}
    source = development.raw.planet_position(shape, metadata)
    column_grid, row_grid = np.indices(shape)
    distance = np.hypot(row_grid - source[0], column_grid - source[1])
    aperture = distance <= settings["aperture_radius"] + 0.5
    searches = [(int(row), int(column)) for column, row in np.argwhere(aperture)]
    searches.sort()
    radius_map = development.image_radius(shape)
    bins: set[int] = set()
    for row, column in searches:
        lower = math.floor(float(radius_map[column, row]) - 0.5)
        bins.update((lower, lower + 1))
    stage.require(searches and min(bins) >= 5 and max(bins) <= 25,
                  "Stage-F aperture left frozen radial coverage")
    return source, aperture, searches, sorted(bins)


def stitched_maps(root: Path, image: np.ndarray, mode: int) -> tuple[np.ndarray, np.ndarray]:
    """Stitch frozen generic maps using the nearest predeclared radial policy."""
    indices = [development.METHODS.index(method) for method in PLANET_METHODS]
    output = np.full((len(PLANET_METHODS), *image.shape), np.nan, dtype=np.float64)
    chosen_radius = np.full(image.shape, np.nan, dtype=np.float64)
    best = np.full(image.shape, np.inf, dtype=np.float64)
    radius_map = development.image_radius(image.shape)
    for nominal in preparation.RADII:
        unit = development.load_unit(root, nominal, mode)
        directory = root / "calibration" / "units" / development.unit_name(nominal, mode)
        all_maps = development.positive_maps(image, nominal, mode, unit, directory)
        values = all_maps[indices]
        common = np.all(np.isfinite(values), axis=0)
        distance = np.abs(radius_map - float(nominal))
        take = common & ((distance < best) |
                         (np.isclose(distance, best) & (float(nominal) < chosen_radius)))
        output[:, take] = values[:, take]
        chosen_radius[take] = float(nominal)
        best[take] = distance[take]
    return output, chosen_radius


def build_amplitude_maps(root: Path, protocol: dict[str, object], science: np.ndarray,
                         settings: dict[str, float]) \
        -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, object], np.ndarray, list[int]]:
    """Apply frozen maps and source-safe aperture weights to the planet cube."""
    source, aperture, searches, bins = analysis_geometry(science.shape[-2:], settings)
    modes = len(stage.MODES)
    maps = np.full((modes, len(PLANET_METHODS), *science.shape[-2:]), np.nan,
                   dtype=np.float64)
    responses = np.full_like(maps, np.nan)
    policy_map = np.full((modes, *science.shape[-2:]), np.nan, dtype=np.float64)
    diagnostics: dict[str, object] = {"source_row_column": list(source),
                                      "aperture_queries": [], "modes": {}}

    baseline = np.asarray(fits.getdata(protocol["paths"]["baseline"], memmap=True),
                          dtype=np.float64)
    parent = Path(str(protocol["parent_response"]))
    paths = development.response47.product_paths(parent)
    coordinates = np.asarray(fits.getdata(paths["coordinates"], memmap=True),
                             dtype=np.float64).T
    coordinate_lookup = {(int(row), int(column)): index
                         for index, (row, column, _, _) in enumerate(coordinates)}
    optimized_planet = development.raw.planet_position(
        science.shape[-2:], {"known_planet": protocol["known_planet"]})
    radius_map = development.image_radius(science.shape[-2:])
    method_indices = {name: index for index, name in enumerate(PLANET_METHODS)}
    gaussian_definitions = (("gaussian_fwhm2p4", 2.4), ("gaussian_fwhm3p6", 3.6))

    for mode_index, mode in enumerate(stage.MODES):
        image = np.asarray(science[mode_index], dtype=np.float64)
        maps[mode_index], policy_map[mode_index] = stitched_maps(root, image, mode)
        exact_responses = np.asarray(fits.getdata(paths["responses"][mode_index], memmap=True),
                                     dtype=np.float64)
        exact_validities = np.asarray(
            fits.getdata(paths["validities"][mode_index], memmap=True),
            dtype=np.float64) > 0.5
        sparse_cube, sparse_validity, sparse_radii = development.sparse_products(
            Path(str(protocol["paths"]["sparse_manifest"])), mode_index)
        gaussian_images = {name: preparation.gaussian_map(image, fwhm)
                           for name, fwhm in gaussian_definitions}
        mode_diagnostics = {}
        for row, column in searches:
            query = (row, column)
            stage.require(query in coordinate_lookup,
                          f"Stage-F aperture query lacks exact response: {query}")
            source_index = coordinate_lookup[query]
            template = development.raw.crop(
                np.asarray(exact_responses[source_index], dtype=np.float64).T, SUPPORT)
            validity = development.raw.crop(
                np.asarray(exact_validities[source_index], dtype=bool).T, SUPPORT)
            center_row = 0.5 * (image.shape[1] - 1)
            center_column = 0.5 * (image.shape[0] - 1)
            actual_radius = math.hypot(row - center_row, column - center_column)
            angle = math.atan2(row - center_row, column - center_column)
            sparse_template, sparse_support = stage.evaluate_sparse_response(
                sparse_cube, sparse_validity, sparse_radii, actual_radius, angle)
            reference_weights, reference_responses, support_counts = planet_reference_weights(
                template, validity, sparse_template, sparse_support)
            nominal = policy_radius(float(radius_map[column, row]))
            width = int(protocol["training"]["candidate_specific_half_width_by_radius"][
                str(nominal)])
            fitted = fit_planet_query(baseline[mode_index], query, searches, template,
                                      validity, optimized_planet, width)
            data = development.stamp(image, query).ravel()

            maps[mode_index, method_indices["native"], column, row] = image[column, row]
            responses[mode_index, method_indices["native"], column, row] = template[HALF, HALF]
            for name, fwhm in gaussian_definitions:
                maps[mode_index, method_indices[name], column, row] = gaussian_images[name][column, row]
                sigma = fwhm / math.sqrt(8 * math.log(2))
                filtered = gaussian_filter(template, sigma=sigma, mode="constant",
                                           cval=0.0, truncate=4.0)
                responses[mode_index, method_indices[name], column, row] = filtered[HALF, HALF]
            for name, weight in reference_weights.items():
                maps[mode_index, method_indices[name], column, row] = weight @ data
                responses[mode_index, method_indices[name], column, row] = reference_responses[name]
            identity = reference_weights["exact_identity"]
            mean_method = "exact_isotropic_fitted_mean"
            maps[mode_index, method_indices[mean_method], column, row] = (
                identity @ (data - fitted["isotropic_mean"]))
            responses[mode_index, method_indices[mean_method], column, row] = (
                identity @ template.ravel())
            for name in ("raw_rectangular_m0p3", "radial_hann_m0p1_trunc0p75"):
                weight = fitted["weights"][name]
                maps[mode_index, method_indices[name], column, row] = weight @ data
                responses[mode_index, method_indices[name], column, row] = weight @ template.ravel()
            policy_map[mode_index, column, row] = nominal
            mode_diagnostics[f"{row},{column}"] = {
                "query_row_column": [row, column],
                "radius": float(radius_map[column, row]),
                "policy_radius": nominal,
                "training_half_width": width,
                "exact_support_pixels": int(np.count_nonzero(validity)),
                "reference_support_pixels": support_counts,
                "radial_profile_minimum_pixels": fitted["profile_minimum_pixels"],
                "raw_samples": fitted["raw_detail"]["policies"][
                    "raw_rectangular_m0p3"]["samples"],
                "radial_samples": fitted["radial_detail"]["policies"][
                    "radial_hann_m0p1_trunc0p75"]["samples"],
                "raw_retained_modes": fitted["raw_detail"]["policies"][
                    "raw_rectangular_m0p3"]["retained_modes"],
                "radial_retained_modes": fitted["radial_detail"]["policies"][
                    "radial_hann_m0p1_trunc0p75"]["retained_modes"],
            }
        diagnostics["modes"][str(mode)] = mode_diagnostics

    stage.require(np.all(np.isfinite(maps[:, :, aperture])) and
                  np.all(np.isfinite(responses[:, :, aperture])),
                  "Stage-F aperture lacks common method support")
    diagnostics["aperture_queries"] = [list(query) for query in searches]
    diagnostics["aperture_pixels"] = len(searches)
    diagnostics["required_radial_bins"] = bins
    diagnostics["optimized_training_exclusion_row_column"] = list(optimized_planet)
    return maps, responses, policy_map, diagnostics, aperture, bins


def run_hcianalyze(root: Path, protocol: dict[str, object], analyzer: Path,
                   maps: np.ndarray, header: fits.Header, output: Path) \
        -> tuple[np.ndarray, fits.Header, list[str]]:
    """Run hciAnalyze with the copied user guidance and no additional filter."""
    flat = maps.reshape((-1, *maps.shape[-2:])).astype(np.float32)
    working_header = header.copy()
    working_header["HCI FILTER LABELS"] = ",".join(
        f"m{mode}_{method}" for mode in stage.MODES for method in PLANET_METHODS)
    fits.writeto(output / "amplitudes.fits", flat, working_header)
    command = [str(analyzer), "--config", str(stage_f_root(root) / "analyze.conf"),
               "--file=amplitudes.fits", "--filter.psfResponse=",
               "--filter.lpfGaussFW=0", "--filter.hpfGaussFW=0",
               "--noise.model=identity", "--noise.only=false",
               "--noise.outputDiagnostics=false"]
    with (output / "analysis.log").open("w", encoding="utf-8") as log:
        subprocess.run(command, cwd=output, env=development.environment(protocol, False),
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    snr, snr_header = fits.getdata(output / "amplitudes_snr.fits", header=True)
    return np.asarray(snr, dtype=np.float64).reshape(maps.shape), snr_header, command


def verify_annular_snr(maps: np.ndarray, snr: np.ndarray,
                       settings: dict[str, float], source: tuple[float, float]) \
        -> dict[str, object]:
    """Verify production SNR against the independent one-pixel annular oracle."""
    exclusions = [(source[0], source[1], settings["source_radius"])]
    radius_map = development.image_radius(maps.shape[-2:])
    methods = {}
    maximum = 0.0
    for mode_index, mode in enumerate(stage.MODES):
        mode_result = {}
        for method_index, method in enumerate(PLANET_METHODS):
            expected = development.annular_oracle(maps[mode_index, method_index], exclusions)
            expected[(radius_map < settings["minimum_radius"]) |
                     (radius_map > settings["maximum_radius"])] = np.nan
            valid = np.isfinite(maps[mode_index, method_index]) & np.isfinite(expected)
            stage.require(np.any(valid) and np.all(np.isfinite(snr[mode_index, method_index][valid])) and
                          np.allclose(snr[mode_index, method_index][valid], expected[valid],
                                      rtol=2e-6, atol=2e-6),
                          f"Stage-F annular oracle mismatch for {mode} {method}")
            error = float(np.max(np.abs(snr[mode_index, method_index][valid] - expected[valid])))
            maximum = max(maximum, error)
            mode_result[method] = {"maximum_error": error,
                                   "checked_pixels": int(np.count_nonzero(valid))}
        methods[str(mode)] = mode_result
    return {"rtol": 2e-6, "atol": 2e-6,
            "maximum_production_oracle_error": maximum, "methods_by_mode": methods}


def method_support(method: str, detail: dict[str, object]) -> int | None:
    """Return the response support count relevant to one aperture method."""
    if method == "native":
        return 1
    if method.startswith("gaussian_"):
        return None
    if method == "sparse_identity":
        return int(detail["reference_support_pixels"]["sparse_identity"])
    return int(detail["exact_support_pixels"])


def summarize_planet(maps: np.ndarray, responses: np.ndarray, snr: np.ndarray,
                     diagnostics: dict[str, object], aperture: np.ndarray,
                     settings: dict[str, float]) -> dict[str, object]:
    """Measure nearest-pixel and aperture-maximum planet statistics."""
    source = tuple(float(value) for value in diagnostics["source_row_column"])
    searches = [tuple(map(int, query)) for query in diagnostics["aperture_queries"]]
    nearest = min(searches, key=lambda query: math.hypot(query[0] - source[0],
                                                         query[1] - source[1]))
    methods_by_mode = {}
    for mode_index, mode in enumerate(stage.MODES):
        mode_result = {}
        for method_index, method in enumerate(PLANET_METHODS):
            values = np.where(aperture, snr[mode_index, method_index], np.nan)
            peak_flat = int(np.nanargmax(values))
            peak_column, peak_row = np.unravel_index(peak_flat, values.shape)
            peak = (int(peak_row), int(peak_column))
            nearest_amplitude = float(maps[mode_index, method_index, nearest[1], nearest[0]])
            peak_amplitude = float(maps[mode_index, method_index, peak[1], peak[0]])
            nearest_response = float(responses[mode_index, method_index, nearest[1], nearest[0]])
            peak_response = float(responses[mode_index, method_index, peak[1], peak[0]])
            stage.require(nearest_response != 0 and peak_response != 0,
                          f"zero Stage-F response for {mode} {method}")
            detail = diagnostics["modes"][str(mode)][f"{peak[0]},{peak[1]}"]
            mode_result[method] = {
                "label": METHOD_LABELS[method],
                "nearest_pixel_row_column": list(nearest),
                "nearest_pixel_snr": float(snr[mode_index, method_index, nearest[1], nearest[0]]),
                "nearest_pixel_amplitude": nearest_amplitude,
                "nearest_pixel_response_to_exact": nearest_response,
                "nearest_pixel_contrast_estimate": nearest_amplitude / nearest_response,
                "aperture_peak_row_column": list(peak),
                "aperture_maximum_snr": float(snr[mode_index, method_index, peak[1], peak[0]]),
                "aperture_peak_amplitude": peak_amplitude,
                "aperture_peak_response_to_exact": peak_response,
                "aperture_peak_contrast_estimate": peak_amplitude / peak_response,
                "peak_offset_from_nominal_pixels": math.hypot(peak[0] - source[0],
                                                                peak[1] - source[1]),
                "valid_aperture_pixels": int(np.count_nonzero(np.isfinite(values))),
                "peak_support_pixels": method_support(method, detail),
                "peak_policy_radius": detail["policy_radius"],
                "peak_training_half_width": detail["training_half_width"],
            }
        methods_by_mode[str(mode)] = mode_result
    return {
        "purpose": "descriptive known-planet endpoint; no method selection",
        "stage_e_policy_changed": False,
        "method_selection_performed": False,
        "known_planet_opened": True,
        "settings": settings,
        "source_row_column": list(source),
        "aperture_pixels": int(np.count_nonzero(aperture)),
        "methods_by_mode": methods_by_mode,
    }


def validation_primary(root: Path) -> tuple[dict[str, dict[str, object]], dict[str, int]]:
    """Return frozen faint-source and held-out primary-mode results."""
    result = read(root / "validation_results.json")
    heldout = read(root / "heldout_results.json")
    rows = [row for row in result["summary_rows"]
            if int(row["mode"]) == policy_module.PRIMARY_MODE and
            float(row["target_source_snr"]) == 3.0]
    summary = {}
    for method in validation.METHODS:
        selected = [row for row in rows if row["method"] == method]
        stage.require(len(selected) == len(preparation.RADII),
                      f"Stage-F closure lost validation method {method}")
        summary[method] = {
            "recoveries": sum(int(row["recoveries"]) for row in selected),
            "sites": sum(int(row["sites"]) for row in selected),
            "mean_maximum_snr": float(np.mean([row["mean_search_snr"] for row in selected])),
        }
    nulls = {method: int(value["frozen_exceedances"])
             for method, value in heldout["sensitivity"][str(policy_module.PRIMARY_MODE)].items()}
    return summary, nulls


def write_tables(root: Path, output: Path, result: dict[str, object],
                 oracle: dict[str, object]) -> None:
    """Write compact machine-readable and Markdown closure tables."""
    validation_summary, nulls = validation_primary(root)
    rows = []
    for mode in stage.MODES:
        for method in PLANET_METHODS:
            value = result["methods_by_mode"][str(mode)][method]
            injection = validation_summary.get(method) if mode == policy_module.PRIMARY_MODE else None
            rows.append({
                "mode": mode,
                "method": method,
                "validation_snr3_recoveries": injection["recoveries"] if injection else "",
                "validation_snr3_sites": injection["sites"] if injection else "",
                "validation_snr3_mean_maximum_snr": injection["mean_maximum_snr"] if injection else "",
                "heldout_exceedances": nulls.get(method, "") if mode == policy_module.PRIMARY_MODE else "",
                "planet_nearest_pixel_snr": value["nearest_pixel_snr"],
                "planet_aperture_maximum_snr": value["aperture_maximum_snr"],
                "planet_peak_row": value["aperture_peak_row_column"][0],
                "planet_peak_column": value["aperture_peak_row_column"][1],
                "planet_peak_offset_pixels": value["peak_offset_from_nominal_pixels"],
                "planet_peak_amplitude": value["aperture_peak_amplitude"],
                "planet_peak_response_to_exact": value["aperture_peak_response_to_exact"],
                "planet_peak_contrast_estimate": value["aperture_peak_contrast_estimate"],
                "planet_peak_support_pixels": ("" if value["peak_support_pixels"] is None
                                                   else value["peak_support_pixels"]),
            })
    with (output / "results.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    primary = result["methods_by_mode"][str(policy_module.PRIMARY_MODE)]
    lines = [
        "# KLIP Stage-F known-planet closure", "",
        "This is a descriptive single-planet endpoint performed after the Stage-E policy and validation result were immutable. It does not select or rescue a method.", "",
        "## Mode-200 closure", "",
        "| Method | Validation SNR-3 recovery | Held-out exceedances | Planet nearest-pixel SNR | Planet aperture-max SNR | Peak row, column | Peak offset (pixels) | Peak contrast estimate |",
        "| :--- | :---: | :---: | ---: | ---: | :---: | ---: | ---: |",
    ]
    for method in PLANET_METHODS:
        injection = validation_summary.get(method)
        recovery = (f"{injection['recoveries']} / {injection['sites']}"
                    if injection is not None else "not in Stage E")
        held = str(nulls[method]) if method in nulls else "not in Stage E"
        value = primary[method]
        row, column = value["aperture_peak_row_column"]
        lines.append(
            f"| {METHOD_LABELS[method]} | {recovery} | {held} | "
            f"{value['nearest_pixel_snr']:.4f} | {value['aperture_maximum_snr']:.4f} | "
            f"{row}, {column} | {value['peak_offset_from_nominal_pixels']:.3f} | "
            f"{value['aperture_peak_contrast_estimate']:.6g} |")
    lines.extend(["", "## Planet aperture-maximum SNR by KL mode", "",
                  "| Method | " + " | ".join(map(str, stage.MODES)) + " |",
                  "| :--- | " + " | ".join("---:" for _ in stage.MODES) + " |"])
    for method in PLANET_METHODS:
        values = [result["methods_by_mode"][str(mode)][method]["aperture_maximum_snr"]
                  for mode in stage.MODES]
        lines.append(f"| {METHOD_LABELS[method]} | " +
                     " | ".join(f"{value:.4f}" for value in values) + " |")
    lines.extend([
        "", "## Verification", "",
        f"- Analysis aperture: {result['aperture_pixels']} native pixels within 3.5 pixels of the `working/analyze.conf` center.",
        f"- Maximum production-versus-independent annular-SNR error: `{oracle['maximum_production_oracle_error']:.7g}`.",
        "- Generic annular maps use the nearest frozen Stage-C radial policy. Aperture weights are fit only on the signal-free baseline and exclude the complete analysis aperture from training.",
        "- Exact and sparse response identity, response-smoothed identity, fitted-mean identity, and both covariance filters retain the planet pixels in candidate support.",
        "", "The injection and held-out columns remain the method-selection evidence. The planet columns are a consistency check on one real source.", "",
    ])
    (output / "results.md").write_text("\n".join(lines), encoding="utf-8")
    stage.write_json(output / "results.json", result)


def archive_partial(output: Path) -> None:
    """Archive an incomplete planet attempt before a clean replay."""
    planet = output / "planet"
    if not planet.exists():
        return
    interrupted = output / "interrupted"
    interrupted.mkdir(exist_ok=True)
    index = 1
    while (interrupted / f"attempt_{index:04d}").exists():
        index += 1
    planet.rename(interrupted / f"attempt_{index:04d}")


def verify_complete(root: Path) -> None:
    """Verify a completed Stage-F receipt and every direct product."""
    output = stage_f_root(root)
    receipt = read(output / "complete.json")
    stage.require(receipt["status"] == "complete" and receipt["known_planet_opened"] and
                  not receipt["stage_e_policy_changed"] and
                  not receipt["method_selection_performed"],
                  "Stage-F completion receipt changed")
    stage.verify(receipt["products"])


def run(root: Path) -> None:
    """Apply every frozen finalist and permanent reference to the known planet."""
    root = root.resolve()
    protocol, analyzer, settings, science_path = enable(root)
    output = stage_f_root(root)
    if (output / "complete.json").exists():
        verify_complete(root)
        print(output / "planet" / "results.md", flush=True)
        return
    archive_partial(output)
    planet_output = output / "planet"
    planet_output.mkdir()
    stage.write_json(output / "state.json", {
        "status": "planet_analyzing", "known_planet_opened": True,
        "stage_e_policy_changed": False, "method_selection_performed": False,
    })
    started = time.monotonic()
    science, header = fits.getdata(science_path, header=True, memmap=True)
    science = np.asarray(science, dtype=np.float64)
    maps, responses, policy_map, diagnostics, aperture, bins = build_amplitude_maps(
        root, protocol, science, settings)
    map_header = header.copy()
    map_header["HCI RADIAL BINS"] = ",".join(map(str, bins))
    fits.writeto(planet_output / "responses.fits",
                 responses.reshape((-1, *responses.shape[-2:])).astype(np.float32), map_header)
    fits.writeto(planet_output / "policy_radius.fits", policy_map.astype(np.float32), map_header)
    stage.write_json(planet_output / "diagnostics.json", diagnostics)
    snr, snr_header, command = run_hcianalyze(root, protocol, analyzer, maps,
                                               map_header, planet_output)
    stage.require(int(snr_header["SNRAPER"]) == int(settings["aperture_radius"]) and
                  int(snr_header["SNRMINR"]) == int(settings["minimum_radius"]) and
                  int(snr_header["SNRMAXR"]) == int(settings["maximum_radius"]) and
                  int(snr_header["SNRMEAN"]) == 1 and int(snr_header["SNRSMALL"]) == 1,
                  "Stage-F hciAnalyze settings changed")
    source = tuple(float(value) for value in diagnostics["source_row_column"])
    oracle = verify_annular_snr(maps, snr, settings, source)
    stage.write_json(planet_output / "annular_verification.json", oracle)
    stage.write_json(planet_output / "command.json", command)
    result = summarize_planet(maps, responses, snr, diagnostics, aperture, settings)
    result["elapsed_seconds"] = time.monotonic() - started
    result["annular_oracle_maximum_error"] = oracle["maximum_production_oracle_error"]
    write_tables(root, planet_output, result, oracle)
    products = [stage.fingerprint(planet_output / name) for name in (
        "amplitudes.fits", "amplitudes_snr.fits", "responses.fits",
        "policy_radius.fits", "diagnostics.json", "annular_verification.json",
        "analysis.log", "command.json", "results.json", "results.csv", "results.md")]
    stage.write_json(output / "complete.json", {
        "status": "complete", "known_planet_opened": True,
        "stage_e_policy_changed": False, "method_selection_performed": False,
        "products": products,
    })
    stage.write_json(output / "state.json", {
        "status": "closure_complete", "known_planet_opened": True,
        "stage_e_policy_changed": False, "method_selection_performed": False,
    })
    verify_complete(root)
    print(planet_output / "results.md", flush=True)


def check(config: Path) -> None:
    """Check method mappings, analysis guidance, and radius-policy boundaries."""
    settings = parse_analysis_config(config.resolve())
    stage.require(len(PLANET_METHODS) == len(set(PLANET_METHODS)) and
                  set(PLANET_METHODS).issubset(development.METHODS) and
                  tuple(validation.METHODS) == policy_module.SELECTED_METHODS,
                  "Stage-F method mapping changed")
    stage.require(policy_radius(8.75) == 7.5 and policy_radius(8.76) == 10.0 and
                  policy_radius(11.0) == 10.0 and policy_radius(11.01) == 12.0 and
                  policy_radius(14.0) == 12.0 and policy_radius(14.01) == 16.0,
                  "Stage-F nearest-radius tie policy changed")
    template = np.arange(1, SUPPORT * SUPPORT + 1, dtype=np.float64).reshape(SUPPORT, SUPPORT)
    sparse = template[::-1].copy()
    validity = np.ones(template.shape, dtype=bool)
    weights, responses, counts = planet_reference_weights(template, validity, sparse, validity)
    stage.require(all(np.isfinite(weight).all() for weight in weights.values()) and
                  all(np.isfinite(value) and value != 0 for value in responses.values()) and
                  all(value == SUPPORT * SUPPORT for value in counts.values()) and
                  settings == EXPECTED_SETTINGS,
                  "Stage-F reference-weight control changed")
    source, _, searches, bins = analysis_geometry((128, 128), settings)
    query = min(searches, key=lambda value: math.hypot(value[0] - source[0],
                                                       value[1] - source[1]))
    optimized = development.raw.planet_position(
        (128, 128), {"known_planet": {"separation": stage.PLANET_SEPARATION,
                                       "position_angle": stage.PLANET_PA}})
    generator = np.random.default_rng(260926)
    image = generator.normal(size=(128, 128))
    coordinate = np.arange(-HALF, HALF + 1)
    column, row = np.meshgrid(coordinate, coordinate)
    model = np.exp(-(np.square(row) + np.square(column)) / (2 * 1.5 ** 2))
    fitted = fit_planet_query(image, query, searches, model, validity, optimized, 20)
    stage.require(len(searches) == 39 and bins == list(range(8, 16)) and
                  np.count_nonzero(fitted["support"]) == SUPPORT * SUPPORT,
                  "Stage-F source-safe full-aperture fit changed")
    print("KLIP Stage-F planet checks passed", flush=True)


def parser() -> argparse.ArgumentParser:
    """Build the Stage-F command-line parser."""
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="action", required=True)
    check_parser = subparsers.add_parser("check")
    check_parser.add_argument("--config", type=Path, default=Path("working/analyze.conf"))
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("root", type=Path)
    prepare_parser.add_argument("--config", type=Path, default=Path("working/analyze.conf"))
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("root", type=Path)
    return result


def main() -> None:
    """Dispatch the requested Stage-F action."""
    args = parser().parse_args()
    if args.action == "check":
        check(args.config)
    elif args.action == "prepare":
        prepare(args.root, args.config)
    else:
        run(args.root)


if __name__ == "__main__":
    main()
