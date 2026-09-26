#!/usr/bin/env python3
"""Run held-out nulls and fresh KLIP validation under the frozen Stage-D policy."""
from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parent))
import freeze_klip_stage_d_policy as policy_module  # noqa: E402
import prepare_klip_stage_c_development as preparation  # noqa: E402
import run_klip_covariance_stage_a as stage  # noqa: E402
import run_klip_stage_c_development as development  # noqa: E402


METHODS = policy_module.SELECTED_METHODS
COVARIANCE_CANDIDATES = policy_module.COVARIANCE_CANDIDATES
WEIGHT_METHODS = (
    "exact_identity",
    "sparse_identity",
    "raw_rectangular_m0p3",
    "radial_hann_m0p1_trunc0p75",
)
SUPPORT = development.SUPPORT
MODEL_ROLES = ("validation", "heldout_null")


def read(path: Path) -> object:
    """Read one JSON document."""
    return json.loads(path.read_text(encoding="utf-8"))


def radius_key(radius: float) -> str:
    """Return the Stage-C JSON key for one nominal radius."""
    return str(float(radius))


def same_content(first: dict[str, object], second: dict[str, object]) -> bool:
    """Return whether two fingerprints describe identical bytes."""
    return (int(first["bytes"]) == int(second["bytes"]) and
            str(first["sha256"]) == str(second["sha256"]))


def record_at(path: Path, content: dict[str, object]) -> dict[str, object]:
    """Describe known fingerprint content at a different resolved path."""
    return {"path": str(path.resolve()), "bytes": int(content["bytes"]),
            "sha256": str(content["sha256"])}


def verify_repair_records(records: list[dict[str, object]]) -> None:
    """Verify repair receipts and their preserved prior-state artifacts."""
    stage.verify(records)
    for record in records:
        repair = read(Path(str(record["path"])))
        stage.verify([repair["previous_runner"], repair["previous_policy_manifest"],
                      repair["previous_policy_completion"]])
        archive = repair.get("repair_detail", {}).get("archive_manifest")
        if archive is not None:
            stage.verify([archive])


def archive_validation_analysis(root: Path) -> dict[str, object]:
    """Archive every partial Stage-E analysis after verifying retained upstream products."""
    protocol = read(root / "protocol.json")
    verify_models(root)
    verify_heldout(root)
    for task in protocol["validation_tasks_unopened"]:
        receipt = read(root / "validation_reductions" / str(task["name"]) / "complete.json")
        stage.verify(receipt["products"])
    stage.require(not (root / "validation_complete.json").exists() and
                  not any((root / name).exists() for name in
                          ("validation_results.json", "validation_results.csv",
                           "validation_results.md")),
                  "final validation products exist; refusing analysis-only repair")
    analysis = root / "validation_analysis"
    stage.require(analysis.is_dir(), "partial validation analysis is missing")
    receipts = sorted(analysis.glob("*/complete.json"))
    for receipt_path in receipts:
        stage.verify(read(receipt_path)["products"])
    directories = [path for path in analysis.iterdir() if path.is_dir()]
    archive_root = root / "interrupted" / "stage_e_analysis_repairs"
    archive_root.mkdir(parents=True, exist_ok=True)
    attempt = 1
    while (archive_root / f"attempt_{attempt:04d}").exists():
        attempt += 1
    destination = archive_root / f"attempt_{attempt:04d}"
    shutil.move(analysis, destination)
    records = [stage.fingerprint(path) for path in sorted(destination.rglob("*"))
               if path.is_file()]
    archive_manifest = destination / "archive_manifest.json"
    stage.write_json(archive_manifest, {
        "schema": 1,
        "reason": "replace the original Stage-C masked-fidelity helper with its recorded repair",
        "completed_analyses": len(receipts),
        "incomplete_analyses": len(directories) - len(receipts),
        "files": records,
    })
    return {"completed_analyses": len(receipts),
            "incomplete_analyses": len(directories) - len(receipts),
            "archive_manifest": stage.fingerprint(archive_manifest)}


def repair_frozen_runner(root: Path, frozen: Path) -> None:
    """Apply a provenance-preserving Stage-E runner repair when its guard permits."""
    source = Path(__file__).resolve()
    stage.require(frozen.is_file(), "frozen Stage-E validation runner is missing")
    manifest_path = root / "policy_manifest.json"
    completion_path = root / "policy_complete.json"
    manifest = read(manifest_path)
    expected = [record for record in manifest["software_records"]
                if Path(str(record["path"])).resolve() == frozen.resolve()]
    stage.require(len(expected) == 1, "policy manifest lacks one frozen validation runner")
    source_record = stage.fingerprint(source)
    observed = stage.fingerprint(frozen)
    if same_content(observed, expected[0]) and same_content(source_record, observed):
        return

    stage.verify(manifest["input_records"] + manifest["policy_records"])
    verify_repair_records(manifest.get("repair_records", []))
    for record in manifest["software_records"]:
        if record != expected[0]:
            stage.verify([record])
    stage.require(same_content(observed, expected[0]),
                  "frozen validation runner changed outside a guarded repair")
    state = read(root / "state.json")
    completion = read(completion_path)
    stage.require(completion["status"] == "complete",
                  "policy completion receipt is invalid")

    if state["status"] == "policy_frozen":
        stage.require(not state["heldout_nulls_opened"] and
                      not state["validation_products_opened"] and
                      not completion["heldout_nulls_opened"] and
                      not completion["validation_products_opened"] and
                      not (root / "validation_models").exists() and
                      not (root / "heldout_analysis").exists() and
                      not (root / "validation_reductions").exists() and
                      not (root / "validation_analysis").exists(),
                      "entry-point repair requires an unopened policy")
        repair_name = "stage_e_entry_unpack_20260926"
        reason = "correct Stage-E main() to unpack the two-value preparation loader"
        repair_detail = {"retained_products": "no Stage-E products existed"}
    elif state["status"] == "validation_analyzing":
        stage.require(state["heldout_nulls_opened"] and
                      state["validation_products_opened"],
                      "analysis repair requires opened Stage-E products")
        repair_name = "stage_e_masked_fidelity_20260926"
        reason = ("use the recorded Stage-C masked-support fidelity repair; detection "
                  "maps, held-out scores, models, and reductions are unchanged")
        repair_detail = archive_validation_analysis(root)
        stage.write_json(root / "state.json", {
            "status": "validation_reductions_complete",
            "development_reductions": len(read(root / "protocol.json")["development_tasks"]),
            "validation_reductions": len(read(root / "protocol.json")[
                "validation_tasks_unopened"]),
            "heldout_nulls_opened": True, "validation_products_opened": True,
            "positive_analysis_started": True})
    else:
        raise RuntimeError("Stage-E runner differs outside a supported repair state")

    repair_directory = root / "policy_repairs" / repair_name
    stage.require(not repair_directory.exists(), f"repair directory already exists: {repair_name}")
    repair_directory.mkdir(parents=True)
    previous_runner = repair_directory / "previous_run_klip_stage_e_validation.py"
    previous_manifest = repair_directory / "previous_policy_manifest.json"
    previous_completion = repair_directory / "previous_policy_complete.json"
    shutil.copy2(frozen, previous_runner)
    shutil.copy2(manifest_path, previous_manifest)
    shutil.copy2(completion_path, previous_completion)
    replacement = record_at(frozen, source_record)
    repair_path = repair_directory / "repair.json"
    stage.write_json(repair_path, {
        "schema": 1, "reason": reason, "scientific_policy_changed": False,
        "heldout_nulls_opened": bool(state["heldout_nulls_opened"]),
        "validation_products_opened": bool(state["validation_products_opened"]),
        "previous_runner": stage.fingerprint(previous_runner),
        "replacement_runner": replacement,
        "previous_policy_manifest": stage.fingerprint(previous_manifest),
        "previous_policy_completion": stage.fingerprint(previous_completion),
        "repair_detail": repair_detail,
    })

    temporary = frozen.with_suffix(frozen.suffix + ".replacement")
    shutil.copy2(source, temporary)
    temporary.replace(frozen)
    stage.require(stage.fingerprint(frozen) == replacement,
                  "replacement validation runner fingerprint changed")
    manifest["software_records"] = [replacement if record == expected[0] else record
                                    for record in manifest["software_records"]]
    repair_record = stage.fingerprint(repair_path)
    manifest["repair_records"] = [*manifest.get("repair_records", []), repair_record]
    stage.write_json(manifest_path, manifest)
    products = [*manifest["policy_records"], *manifest["repair_records"],
                stage.fingerprint(manifest_path)]
    stage.write_json(completion_path, {
        "status": "complete",
        "heldout_nulls_opened": bool(completion["heldout_nulls_opened"]),
        "validation_products_opened": bool(completion["validation_products_opened"]),
        "products": products})
    stage.verify(products + manifest["input_records"] + manifest["software_records"])


def ensure_frozen_runner(root: Path) -> None:
    """Execute the policy-frozen runner, applying the guarded pre-exposure repair once."""
    source = Path(__file__).resolve()
    frozen = root / "software" / source.name
    repair_frozen_runner(root, frozen)
    if source != frozen.resolve():
        os.execv(sys.executable, [sys.executable, str(frozen), *sys.argv[1:]])


def enable(root: Path) -> tuple[dict[str, object], dict[str, object], Path]:
    """Verify the immutable policy and return its protocol, policy, and analyzer."""
    protocol, _ = preparation.load_experiment(root)
    completion = read(root / "policy_complete.json")
    stage.require(completion["status"] == "complete", "Stage-D policy is incomplete")
    stage.verify(completion["products"])
    manifest = read(root / "policy_manifest.json")
    stage.verify(manifest["input_records"] + manifest["software_records"] +
                 manifest["policy_records"])
    verify_repair_records(manifest.get("repair_records", []))
    own = stage.fingerprint(Path(__file__))
    stage.require(own in manifest["software_records"],
                  "validation runner is outside the frozen policy")
    policy = read(root / "policy" / "policy.json")
    stage.require(policy["frozen_before_heldout_or_validation"] and
                  policy["primary_mode"] == policy_module.PRIMARY_MODE and
                  tuple(policy["methods"]["ordered"]) == METHODS and
                  tuple(policy["methods"]["covariance_candidates"]) == COVARIANCE_CANDIDATES and
                  not policy["joint_mode_maximization"],
                  "frozen validation policy changed")
    analyzer = Path(str(read(root / "development_manifest.json")["hcianalyze"]["path"]))
    stage.require(analyzer.is_file(), "frozen hciAnalyze executable is missing")
    return protocol, policy, analyzer


def model_unit_name(radius: float, mode: int) -> str:
    """Return one validation-model unit directory name."""
    return f"r{development.radius_tag(radius)}_mode{mode}"


def model_sites(protocol: dict[str, object], radius: float) -> list[dict[str, object]]:
    """Return validation and held-out sites in a fixed order."""
    sites = [site for site in protocol["sites"]
             if float(site["nominal_radius"]) == float(radius) and site["role"] in MODEL_ROLES]
    sites.sort(key=lambda site: (MODEL_ROLES.index(str(site["role"])), int(site["role_index"])))
    expected = sum(int(protocol["sites_per_radius_by_role"][role]) for role in MODEL_ROLES)
    stage.require(len(sites) == expected, f"validation model lost radius {radius:g} sites")
    return sites


def calculate_model_unit(root_value: str, radius_value: float, mode_index: int) -> str:
    """Freeze baseline-derived weights for all unopened sites at one radius and mode."""
    root = Path(root_value)
    protocol = read(root / "protocol.json")
    radius_value = float(radius_value)
    mode = stage.MODES[mode_index]
    directory = root / "validation_models" / "units" / model_unit_name(radius_value, mode)
    directory.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()

    baseline_cube = fits.getdata(protocol["paths"]["baseline"], memmap=True)
    image = np.asarray(baseline_cube[mode_index], dtype=np.float64)
    parent = Path(str(protocol["parent_response"]))
    paths = development.response47.product_paths(parent)
    coordinates = np.asarray(fits.getdata(paths["coordinates"], memmap=True), dtype=np.float64).T
    coordinate_lookup = {(int(row), int(column)): index
                         for index, (row, column, _, _) in enumerate(coordinates)}
    responses = np.asarray(fits.getdata(paths["responses"][mode_index], memmap=True), dtype=np.float64)
    validities = np.asarray(fits.getdata(paths["validities"][mode_index], memmap=True),
                            dtype=np.float64) > 0.5
    planet = development.raw.planet_position(image.shape,
                                              {"known_planet": protocol["known_planet"]})
    gaussian_maps = {name: preparation.gaussian_map(image, fwhm)
                     for name, fwhm in development.GAUSSIAN_METHODS if name in METHODS}
    unit = development.load_unit(root, radius_value, mode)
    lookup = {tuple(map(int, query)): index for index, query in enumerate(unit["positions"])}
    sites = model_sites(protocol, radius_value)
    positions = np.empty((len(sites), 5, 2), dtype=np.int16)
    weights = np.empty((len(sites), 5, len(WEIGHT_METHODS), SUPPORT * SUPPORT),
                       dtype=np.float64)
    responses_out = np.empty((len(sites), 5, len(METHODS)), dtype=np.float64)
    baseline_out = np.empty((len(sites), 5, len(METHODS)), dtype=np.float64)
    sigmas = np.empty((len(sites), 5, len(COVARIANCE_CANDIDATES)), dtype=np.float64)
    diagnostics = {}
    width = int(protocol["training"]["candidate_specific_half_width_by_radius"][
        radius_key(radius_value)])

    for site_index, site in enumerate(sites):
        searches = [(int(site["row"]) + delta_row, int(site["column"]) + delta_column)
                    for delta_row, delta_column in development.footprint.SEARCH_OFFSETS]
        details = []
        for search_index, query in enumerate(searches):
            stage.require(query in lookup and query in coordinate_lookup,
                          f"site query is outside the calibrated annulus: {site['name']}")
            local = lookup[query]
            source = coordinate_lookup[query]
            template = development.raw.crop(np.asarray(responses[source], dtype=np.float64).T,
                                             SUPPORT)
            validity = development.raw.crop(np.asarray(validities[source], dtype=bool).T, SUPPORT)
            fitted = development.fit_query(image, query, searches, template, validity, planet, width)
            data = development.stamp(image, query).ravel()
            exact = unit["reference_weights"][
                local, development.REFERENCE_WEIGHT_METHODS.index("exact_identity")]
            sparse = unit["reference_weights"][
                local, development.REFERENCE_WEIGHT_METHODS.index("sparse_identity")]
            raw_weight = fitted["weights"]["raw_rectangular_m0p3"]
            radial_weight = fitted["weights"]["radial_hann_m0p1_trunc0p75"]
            local_weights = (exact, sparse, raw_weight, radial_weight)
            positions[site_index, search_index] = query
            weights[site_index, search_index] = np.stack(local_weights)
            for method_index, method in enumerate(METHODS):
                if method == "native":
                    amplitude = float(image[query[1], query[0]])
                    response = float(template[development.HALF, development.HALF])
                elif method in gaussian_maps:
                    amplitude = float(gaussian_maps[method][query[1], query[0]])
                    response = float(unit["source_responses"][local,
                                                                  development.METHODS.index(method)])
                else:
                    weight = local_weights[WEIGHT_METHODS.index(method)]
                    amplitude = float(weight @ data)
                    response = float(weight @ template.ravel())
                stage.require(np.isfinite(amplitude) and np.isfinite(response) and response != 0,
                              f"invalid frozen model for {site['name']} {mode} {method}")
                baseline_out[site_index, search_index, method_index] = amplitude
                responses_out[site_index, search_index, method_index] = response
            sigmas[site_index, search_index, 0] = fitted["sigmas"]["raw_rectangular_m0p3"]
            sigmas[site_index, search_index, 1] = fitted["sigmas"][
                "radial_hann_m0p1_trunc0p75"]
            details.append({
                "query_row_column": list(query),
                "profile_minimum_pixels": fitted["profile_minimum_pixels"],
                "raw": fitted["raw_detail"]["policies"]["raw_rectangular_m0p3"],
                "radial": fitted["radial_detail"]["policies"][
                    "radial_hann_m0p1_trunc0p75"],
            })
        diagnostics[str(site["name"])] = details

    np.savez_compressed(directory / "models.npz",
                        site_names=np.asarray([str(site["name"]) for site in sites]),
                        positions=positions, weights=weights, source_responses=responses_out,
                        baseline_amplitudes=baseline_out, conditional_sigmas=sigmas)
    stage.write_json(directory / "diagnostics.json", {
        "radius": radius_value,
        "mode": mode,
        "baseline_only": True,
        "positive_products_read": False,
        "training_half_width": width,
        "sites": diagnostics,
        "elapsed_seconds": time.monotonic() - started,
    })
    products = [stage.fingerprint(directory / name)
                for name in ("models.npz", "diagnostics.json")]
    stage.write_json(directory / "complete.json", {"status": "complete", "products": products})
    return str(directory / "complete.json")


def load_model(root: Path, radius: float, mode: int) -> dict[str, np.ndarray]:
    """Load one completed validation-model unit."""
    directory = root / "validation_models" / "units" / model_unit_name(radius, mode)
    stage.verify(read(directory / "complete.json")["products"])
    with np.load(directory / "models.npz", allow_pickle=False) as source:
        return {name: np.array(source[name], copy=True) for name in source.files}


def verify_models(root: Path) -> None:
    """Verify the complete score-blind validation-model receipt."""
    receipt = read(root / "validation_models_complete.json")
    stage.require(receipt["status"] == "complete" and receipt["baseline_only"] and
                  not receipt["heldout_scores_read"] and not receipt["validation_products_opened"],
                  "validation model receipt is invalid")
    stage.verify(receipt["products"])
    paths = sorted((root / "validation_models" / "units").glob("*/complete.json"))
    stage.require(len(paths) == len(preparation.RADII) * len(stage.MODES),
                  "validation model receipt count changed")
    for path in paths:
        stage.verify(read(path)["products"])


def models(root: Path, workers: int) -> None:
    """Build baseline-only candidate weights before reading any held-out or positive score."""
    protocol, _, _ = enable(root)
    if (root / "validation_models_complete.json").exists():
        verify_models(root)
        print("Stage-E validation models already complete.", flush=True)
        return
    state = read(root / "state.json")
    stage.require(state["status"] == "policy_frozen" and
                  not state["heldout_nulls_opened"] and
                  not state["validation_products_opened"],
                  "validation models require a fresh frozen policy")
    model_root = root / "validation_models" / "units"
    model_root.mkdir(parents=True, exist_ok=True)
    pending, receipts = [], []
    for radius in preparation.RADII:
        for mode_index, mode in enumerate(stage.MODES):
            directory = model_root / model_unit_name(radius, mode)
            receipt = directory / "complete.json"
            if receipt.exists():
                stage.verify(read(receipt)["products"])
                receipts.append(receipt)
            else:
                if directory.exists():
                    development.archive_incomplete(root, directory, "validation_models")
                pending.append((radius, mode_index))
    if pending:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(calculate_model_unit, str(root), radius, mode_index):
                       (radius, stage.MODES[mode_index]) for radius, mode_index in pending}
            for completed, future in enumerate(as_completed(futures), start=1):
                receipts.append(Path(future.result()))
                radius, mode = futures[future]
                print(f"validation model {completed}/{len(pending)}: radius {radius:g}, mode {mode}",
                      flush=True)
    stage.require(len(receipts) == len(preparation.RADII) * len(stage.MODES),
                  "validation model unit count changed")
    products = [stage.fingerprint(path) for path in sorted(receipts)]
    stage.write_json(root / "validation_models_complete.json", {
        "status": "complete",
        "baseline_only": True,
        "heldout_scores_read": False,
        "validation_products_opened": False,
        "products": products,
    })
    stage.write_json(root / "state.json", {
        "status": "validation_models_complete",
        "development_reductions": len(protocol["development_tasks"]),
        "heldout_nulls_opened": False,
        "validation_products_opened": False,
        "positive_analysis_started": True,
    })


def selected_baseline_maps(root: Path, radius: float, mode: int) -> np.ndarray:
    """Load only the policy-selected generic baseline amplitude maps."""
    directory = root / "calibration" / "units" / development.unit_name(radius, mode)
    values = np.asarray(fits.getdata(directory / "baseline_amplitudes.fits"), dtype=np.float64)
    return values[[development.METHODS.index(method) for method in METHODS]]


def selected_positive_maps(image: np.ndarray, root: Path, radius: float, mode: int) -> np.ndarray:
    """Apply the Stage-C generic weights and retain only policy-selected maps."""
    unit = development.load_unit(root, radius, mode)
    directory = root / "calibration" / "units" / development.unit_name(radius, mode)
    values = development.positive_maps(image, radius, mode, unit, directory)
    return values[[development.METHODS.index(method) for method in METHODS]]


def production_snr(protocol: dict[str, object], analyzer: Path, maps: np.ndarray,
                   header: fits.Header, site: dict[str, object], directory: Path) \
        -> tuple[np.ndarray, float, list[dict[str, object]]]:
    """Run production annular normalization and verify the independent oracle."""
    directory.mkdir(parents=True, exist_ok=True)
    planet = protocol["known_planet"]
    command = [str(analyzer), "--file=amplitudes.fits", f"--lambdaD={stage.LAMBDA_D}",
               f"--planet.sep={planet['separation']},{site['separation']}",
               f"--planet.PA={planet['position_angle']},{site['position_angle']}",
               f"--planet.R={planet['exclusion_radius']},{planet['exclusion_radius']}",
               "--snr.apertureR=60", "--snr.minRad=0", "--snr.maxRad=60",
               "--filter.psfResponse=", "--filter.lpfGaussFW=0", "--filter.hpfGaussFW=0",
               "--noise.model=identity", "--noise.only=false", "--noise.outputDiagnostics=false"]
    flat = maps.reshape((-1, *maps.shape[-2:])).astype(np.float32)
    working_header = header.copy()
    working_header["HCI FILTER LABELS"] = ",".join(
        f"m{mode}_{method}" for mode in stage.MODES for method in METHODS)
    with tempfile.TemporaryDirectory(dir=directory) as temporary_value:
        temporary = Path(temporary_value)
        fits.writeto(temporary / "amplitudes.fits", flat, working_header)
        with (directory / "analysis.log").open("w", encoding="utf-8") as log:
            subprocess.run(command, cwd=temporary, env=development.environment(protocol, False),
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        snr, snr_header = fits.getdata(temporary / "amplitudes_snr.fits", header=True)
    snr = np.asarray(snr, dtype=np.float64).reshape(maps.shape)
    stage.require(int(snr_header["SNRMEAN"]) == 1 and int(snr_header["SNRSMALL"]) == 1,
                  "hciAnalyze annular SNR contract changed")
    planet_position = development.raw.planet_position(maps.shape[-2:],
                                                       {"known_planet": planet})
    exclusions = [(planet_position[0], planet_position[1], float(planet["exclusion_radius"])),
                  (float(site["row"]), float(site["column"]),
                   float(planet["exclusion_radius"]))]
    checked = [(int(site["row"]) + delta_row, int(site["column"]) + delta_column)
               for delta_row, delta_column in development.footprint.SEARCH_OFFSETS]
    radii = development.image_radius(maps.shape[-2:])
    maximum_error = 0.0
    fallbacks = []
    for mode_index, mode in enumerate(stage.MODES):
        for method_index, method in enumerate(METHODS):
            expected = development.annular_oracle(maps[mode_index, method_index], exclusions)
            supported = development.annular_oracle(maps[mode_index, method_index], exclusions,
                                                    supported_edge=True)
            for row, column in checked:
                observed = float(snr[mode_index, method_index, column, row])
                oracle = float(expected[column, row])
                if np.isfinite(oracle):
                    stage.require(np.isfinite(observed) and
                                  np.isclose(observed, oracle, rtol=2e-6, atol=2e-6),
                                  f"annular oracle mismatch for {method}")
                    maximum_error = max(maximum_error, abs(observed - oracle))
                    continue
                radius = float(radii[column, row])
                replacement = float(supported[column, row])
                stage.require(radius < stage.SNR_MIN_RADIUS + 0.5 and observed == 0 and
                              np.isfinite(replacement),
                              f"unsupported non-boundary annular SNR for {method}")
                snr[mode_index, method_index, column, row] = replacement
                fallbacks.append({
                    "mode": mode, "method": method, "row": row, "column": column,
                    "radius": radius, "production_value": observed, "replacement": replacement,
                    "policy": "nearest supported one-pixel annular mean and sample deviation; exact small-sample correction at candidate radius",
                })
    stage.write_json(directory / "command.json", command)
    stage.write_json(directory / "annular_verification.json", {
        "maximum_production_oracle_error": maximum_error,
        "boundary_fallbacks": fallbacks,
    })
    return snr, maximum_error, fallbacks


def site_model(model: dict[str, np.ndarray], site_name: str) -> tuple[int, np.ndarray]:
    """Return the unique index and search positions for one frozen site model."""
    names = [str(value) for value in model["site_names"]]
    stage.require(names.count(site_name) == 1, f"site model is missing {site_name}")
    index = names.index(site_name)
    return index, np.asarray(model["positions"][index], dtype=np.int64)


def heldout_site(root_value: str, site_name: str) -> str:
    """Measure one held-out null under unchanged policy thresholds."""
    root = Path(root_value)
    protocol = read(root / "protocol.json")
    policy = read(root / "policy" / "policy.json")
    analyzer = Path(str(read(root / "development_manifest.json")["hcianalyze"]["path"]))
    site = next(site for site in protocol["sites"] if site["name"] == site_name)
    stage.require(site["role"] == "heldout_null", "held-out analysis received another role")
    radius = float(site["nominal_radius"])
    directory = root / "heldout_analysis" / site_name
    complete = directory / "complete.json"
    if complete.exists():
        stage.verify(read(complete)["products"])
        return str(directory / "scores.json")
    if directory.exists():
        development.archive_incomplete(root, directory, "heldout_analysis")
    directory.mkdir(parents=True)
    _, header = fits.getdata(protocol["paths"]["baseline"], header=True, memmap=True)
    maps = []
    raw_amplitudes = []
    for mode in stage.MODES:
        selected = selected_baseline_maps(root, radius, mode)
        model = load_model(root, radius, mode)
        index, positions = site_model(model, site_name)
        amplitudes = np.asarray(model["baseline_amplitudes"][index], dtype=np.float64)
        for search_index, (row, column) in enumerate(positions):
            selected[:, int(column), int(row)] = amplitudes[search_index]
        maps.append(selected)
        raw_amplitudes.append(amplitudes)
    maps_array = np.stack(maps)
    snr, maximum_error, fallbacks = production_snr(protocol, analyzer, maps_array, header,
                                                    site, directory)
    searches = [(int(site["row"]) + delta_row, int(site["column"]) + delta_column)
                for delta_row, delta_column in development.footprint.SEARCH_OFFSETS]
    modes = {}
    for mode_index, mode in enumerate(stage.MODES):
        methods = {}
        for method_index, method in enumerate(METHODS):
            values = [float(snr[mode_index, method_index, column, row])
                      for row, column in searches]
            amplitudes = [float(raw_amplitudes[mode_index][search_index, method_index])
                          for search_index in range(5)]
            threshold = float(policy["thresholds"][radius_key(radius)][str(mode)][method])
            stage.require(np.all(np.isfinite(values)) and np.all(np.isfinite(amplitudes)),
                          f"held-out site lacks {method} support")
            methods[method] = {
                "snr_pixels": values,
                "search_score": max(values),
                "amplitude_pixels": amplitudes,
                "threshold": threshold,
                "exceeds_threshold": max(values) > threshold,
            }
        modes[str(mode)] = methods
    result = {"site": site, "methods_by_mode": modes,
              "annular_oracle_maximum_error": maximum_error,
              "annular_boundary_fallbacks": fallbacks}
    stage.write_json(directory / "scores.json", result)
    products = [stage.fingerprint(directory / name) for name in
                ("scores.json", "analysis.log", "command.json", "annular_verification.json")]
    stage.write_json(complete, {"status": "complete", "products": products})
    return str(directory / "scores.json")


def heldout_sensitivity(root: Path, protocol: dict[str, object],
                        records: list[dict[str, object]]) -> dict[str, object]:
    """Evaluate held-out exceedances over the fixed paired calibration resamples."""
    with np.load(root / "policy" / "resamples.npz", allow_pickle=False) as source:
        indices = np.asarray(source["bootstrap_indices"], dtype=np.int64)
    calibration = read(root / "calibration" / "thresholds.json")
    lookup = {str(record["site"]["name"]): record for record in records}
    result = {}
    for mode in stage.MODES:
        methods = {}
        for method in METHODS:
            counts = np.zeros(policy_module.BOOTSTRAP_REPLICATES, dtype=np.int16)
            frozen_count = 0
            by_radius = {}
            for radius_index, radius in enumerate(preparation.RADII):
                maxima = np.asarray(calibration[radius_key(radius)][str(mode)][method]["site_maxima"],
                                    dtype=np.float64)
                resampled = np.max(maxima[indices[radius_index]], axis=1)
                held_sites = policy_module.selected_sites(protocol, radius, "heldout_null")
                values = np.asarray([
                    lookup[str(site["name"])]["methods_by_mode"][str(mode)][method]["search_score"]
                    for site in held_sites], dtype=np.float64)
                count = np.sum(values[:, None] > resampled[None, :], axis=0)
                counts += count.astype(np.int16)
                threshold = float(calibration[radius_key(radius)][str(mode)][method]["threshold"])
                frozen_radius_count = int(np.sum(values > threshold))
                frozen_count += frozen_radius_count
                by_radius[radius_key(radius)] = {
                    "frozen_exceedances": frozen_radius_count,
                    "bootstrap_exceedances": policy_module.quantiles(count),
                }
            methods[method] = {"frozen_exceedances": frozen_count,
                               "bootstrap_exceedances": policy_module.quantiles(counts),
                               "by_radius": by_radius}
        result[str(mode)] = methods
    return result


def verify_heldout(root: Path) -> dict[str, object]:
    """Verify the complete held-out-null receipt."""
    receipt = read(root / "heldout_complete.json")
    stage.require(receipt["status"] == "complete" and receipt["policy_unchanged"],
                  "held-out receipt is invalid")
    stage.verify(receipt["products"])
    paths = sorted((root / "heldout_analysis").glob("*/complete.json"))
    stage.require(len(paths) == len(preparation.RADII) *
                  preparation.ROLE_COUNTS["heldout_null"],
                  "held-out receipt count changed")
    for path in paths:
        stage.verify(read(path)["products"])
    return read(root / "heldout_results.json")


def heldout(root: Path, workers: int) -> None:
    """Expose and analyze the preassigned held-out nulls without retuning."""
    protocol, _, _ = enable(root)
    verify_models(root)
    if (root / "heldout_complete.json").exists():
        verify_heldout(root)
        print(root / "heldout_results.json", flush=True)
        return
    state = read(root / "state.json")
    stage.require(state["status"] in {"validation_models_complete", "heldout_running"} and
                  not state["validation_products_opened"],
                  "held-out nulls require completed score-blind models and unopened validation")
    analysis = root / "heldout_analysis"
    analysis.mkdir(exist_ok=True)
    sites = [site for site in protocol["sites"] if site["role"] == "heldout_null"]
    paths = []
    stage.write_json(root / "state.json", {
        "status": "heldout_running", "development_reductions": len(protocol["development_tasks"]),
        "heldout_nulls_opened": True, "validation_products_opened": False,
        "positive_analysis_started": True})
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(heldout_site, str(root), str(site["name"])): str(site["name"])
                   for site in sites}
        for completed, future in enumerate(as_completed(futures), start=1):
            paths.append(Path(future.result()))
            print(f"held-out null {completed}/{len(sites)}: {futures[future]}", flush=True)
    records = [read(path) for path in sorted(paths)]
    sensitivity = heldout_sensitivity(root, protocol, records)
    result = {
        "purpose": "KLIP Stage-E held-out null result under frozen Stage-D policy",
        "policy_changed": False,
        "sites": len(records),
        "records": records,
        "sensitivity": sensitivity,
        "primary_mode": policy_module.PRIMARY_MODE,
    }
    stage.write_json(root / "heldout_results.json", result)
    products = [stage.fingerprint(root / "heldout_results.json")]
    products.extend(stage.fingerprint(root / "heldout_analysis" / str(site["name"]) /
                                      "complete.json") for site in sites)
    stage.write_json(root / "heldout_complete.json", {
        "status": "complete", "policy_unchanged": True, "products": products})
    stage.write_json(root / "state.json", {
        "status": "heldout_complete", "development_reductions": len(protocol["development_tasks"]),
        "heldout_nulls_opened": True, "validation_products_opened": False,
        "positive_analysis_started": True})
    print(root / "heldout_results.json", flush=True)


def validate_reduction(path: Path) -> dict[str, object]:
    """Validate one positive KLIP cube and return its fingerprint."""
    data, header = fits.getdata(path, header=True, memmap=True)
    stage.require(np.asarray(data).shape == (len(stage.MODES), 128, 128) and
                  stage.read_modes(header) == stage.MODES,
                  "validation reduction schema changed")
    return stage.fingerprint(path)


def reduce(root: Path) -> None:
    """Run or verify all 108 frozen validation reductions."""
    protocol, _, _ = enable(root)
    verify_models(root)
    verify_heldout(root)
    if (root / "validation_complete.json").exists():
        receipt = read(root / "validation_complete.json")
        stage.verify(receipt["products"])
        print("Stage-E validation already complete.", flush=True)
        return
    commands = read(root / "policy" / "validation_commands.json")
    stage.require(len(commands) == int(protocol["expected_future_validation_reductions"]),
                  "validation command count changed")
    reductions = root / "validation_reductions"
    reductions.mkdir(exist_ok=True)
    stage.write_json(root / "state.json", {
        "status": "validation_reducing", "development_reductions": len(protocol["development_tasks"]),
        "validation_reductions": len(list(reductions.glob("*/complete.json"))),
        "heldout_nulls_opened": True, "validation_products_opened": True,
        "positive_analysis_started": True})
    for index, record in enumerate(commands, start=1):
        name = str(record["task"])
        directory = reductions / name
        complete = directory / "complete.json"
        if complete.exists():
            stage.verify(read(complete)["products"])
        else:
            if directory.exists():
                development.archive_incomplete(root, directory, "validation_reductions")
            directory.mkdir()
            elapsed = stage.run_command([str(value) for value in record["command"]], directory,
                                        directory / "run.log", development.environment(protocol, True))
            product = validate_reduction(directory / "finim.fits")
            stage.write_json(complete, {
                "status": "complete", "elapsed_seconds": elapsed,
                "products": [product, stage.fingerprint(directory / "run.log"),
                             stage.fingerprint(directory / "run.command.json")],
            })
        print(f"validation reduction {index}/{len(commands)}: {name}", flush=True)
    stage.write_json(root / "state.json", {
        "status": "validation_reductions_complete",
        "development_reductions": len(protocol["development_tasks"]),
        "validation_reductions": len(commands),
        "heldout_nulls_opened": True, "validation_products_opened": True,
        "positive_analysis_started": True})


def response_fidelity(delta: np.ndarray, template: np.ndarray,
                      fitted: dict[str, object]) -> dict[str, object]:
    """Compare a finite response using the repaired full-space support convention."""
    mask = np.asarray(fitted["support"], dtype=bool).ravel()
    raw_delta = np.asarray(delta, dtype=np.float64).ravel()[mask]
    raw_template = np.asarray(template, dtype=np.float64).ravel()[mask]

    def metric(data: np.ndarray, model_template: np.ndarray,
               covariance: np.ndarray) -> dict[str, float]:
        covariance = np.asarray(covariance, dtype=np.float64)
        stage.require(covariance.shape == (mask.size, mask.size),
                      "fidelity covariance and response support differ")
        submatrix = covariance[np.ix_(mask, mask)]
        inverse_data = np.linalg.solve(submatrix, data)
        inverse_template = np.linalg.solve(submatrix, model_template)
        cross = float(model_template @ inverse_data)
        template_energy = float(model_template @ inverse_template)
        data_energy = float(data @ inverse_data)
        projection = cross / template_energy
        cosine = cross / math.sqrt(template_energy * data_energy)
        residual = data - projection * model_template
        residual_energy = float(residual @ np.linalg.solve(submatrix, residual))
        return {"cosine": cosine, "projection_scale": projection,
                "best_scaled_relative_residual": math.sqrt(residual_energy / data_energy)}

    result = {"unweighted": metric(raw_delta, raw_template,
                                     np.eye(mask.size, dtype=np.float64))}
    result["raw_rectangular_m0p3"] = metric(
        raw_delta, raw_template, fitted["raw_detail"]["covariance"])
    scale = np.asarray(fitted["scale"], dtype=np.float64).ravel()[mask]
    result["radial_hann_m0p1"] = metric(
        raw_delta / scale, raw_template / scale, fitted["radial_detail"]["covariance"])
    return result


def analyze_task(root_value: str, task_name: str) -> str:
    """Analyze one unopened validation positive using baseline-frozen selected filters."""
    root = Path(root_value)
    protocol = read(root / "protocol.json")
    analyzer = Path(str(read(root / "development_manifest.json")["hcianalyze"]["path"]))
    task = next(task for task in protocol["validation_tasks_unopened"] if task["name"] == task_name)
    site = next(site for site in protocol["sites"] if site["name"] == task["site"])
    radius = float(site["nominal_radius"])
    directory = root / "validation_analysis" / task_name
    complete = directory / "complete.json"
    if complete.exists():
        stage.verify(read(complete)["products"])
        return str(directory / "measurements.json")
    if directory.exists():
        development.archive_incomplete(root, directory, "validation_analysis")
    directory.mkdir(parents=True)
    started = time.monotonic()
    source = root / "validation_reductions" / task_name / "finim.fits"
    validate_reduction(source)
    positive_cube, header = fits.getdata(source, header=True, memmap=True)
    baseline_cube = fits.getdata(protocol["paths"]["baseline"], memmap=True)
    parent = Path(str(protocol["parent_response"]))
    paths = development.response47.product_paths(parent)
    coordinates = np.asarray(fits.getdata(paths["coordinates"], memmap=True), dtype=np.float64).T
    coordinate_lookup = {(int(row), int(column)): index
                         for index, (row, column, _, _) in enumerate(coordinates)}
    planet = development.raw.planet_position((128, 128),
                                              {"known_planet": protocol["known_planet"]})
    searches = [(int(site["row"]) + delta_row, int(site["column"]) + delta_column)
                for delta_row, delta_column in development.footprint.SEARCH_OFFSETS]
    width = int(site["training_half_width"])
    maps, source_responses, baseline_amplitudes = [], [], []
    fidelity_records, diagnostics = {}, {}
    for mode_index, mode in enumerate(stage.MODES):
        positive = np.asarray(positive_cube[mode_index], dtype=np.float64)
        baseline = np.asarray(baseline_cube[mode_index], dtype=np.float64)
        selected = selected_positive_maps(positive, root, radius, mode)
        model = load_model(root, radius, mode)
        site_index, positions = site_model(model, str(site["name"]))
        stage.require(np.array_equal(positions, np.asarray(searches)),
                      "validation model search order changed")
        local_weights = np.asarray(model["weights"][site_index], dtype=np.float64)
        local_responses = np.asarray(model["source_responses"][site_index], dtype=np.float64)
        local_baseline = np.asarray(model["baseline_amplitudes"][site_index], dtype=np.float64)
        for search_index, query in enumerate(searches):
            data = development.stamp(positive, query).ravel()
            for weighted_index, method in enumerate(WEIGHT_METHODS):
                selected[METHODS.index(method), query[1], query[0]] = (
                    local_weights[search_index, weighted_index] @ data)
        maps.append(selected)
        source_responses.append(local_responses)
        baseline_amplitudes.append(local_baseline)

        query = searches[0]
        source_index = coordinate_lookup[query]
        response_cube = np.asarray(fits.getdata(paths["responses"][mode_index], memmap=True),
                                   dtype=np.float64)
        validity_cube = np.asarray(fits.getdata(paths["validities"][mode_index], memmap=True),
                                   dtype=np.float64) > 0.5
        template = development.raw.crop(response_cube[source_index].T, SUPPORT)
        validity = development.raw.crop(validity_cube[source_index].T, SUPPORT)
        fitted = development.fit_query(baseline, query, searches, template, validity, planet, width)
        for weighted_index, method in enumerate(WEIGHT_METHODS[2:], start=2):
            stage.require(np.allclose(fitted["weights"][method], local_weights[0, weighted_index],
                                     rtol=3e-10, atol=3e-12),
                          f"recomputed baseline weight changed for {method}")
        delta = ((development.stamp(positive, query) - development.stamp(baseline, query)) /
                 float(task["contrast"]))
        fidelity_records[str(mode)] = response_fidelity(delta, template, fitted)
        diagnostics[str(mode)] = {
            "profile_minimum_pixels": fitted["profile_minimum_pixels"],
            "raw": fitted["raw_detail"]["policies"]["raw_rectangular_m0p3"],
            "radial": fitted["radial_detail"]["policies"][
                "radial_hann_m0p1_trunc0p75"],
        }
    maps_array = np.stack(maps)
    snr, maximum_error, fallbacks = production_snr(protocol, analyzer, maps_array, header,
                                                    site, directory)
    response_array = np.stack(source_responses)
    baseline_array = np.stack(baseline_amplitudes)
    modes = {}
    for mode_index, mode in enumerate(stage.MODES):
        method_records = {}
        for method_index, method in enumerate(METHODS):
            values = [float(snr[mode_index, method_index, column, row])
                      for row, column in searches]
            amplitudes = [float(maps_array[mode_index, method_index, column, row])
                          for row, column in searches]
            stage.require(np.all(np.isfinite(values)) and np.all(np.isfinite(amplitudes)),
                          f"validation task lacks {method} support")
            best = int(np.argmax(values))
            response = float(response_array[mode_index, 0, method_index])
            center_estimate = amplitudes[0] / response
            baseline_center = float(baseline_array[mode_index, 0, method_index]) / response
            method_records[method] = {
                "snr_pixels": values, "search_score": values[best],
                "center_snr": values[0], "best_search_index": best,
                "localization_error_pixels": math.hypot(
                    *development.footprint.SEARCH_OFFSETS[best]),
                "amplitude_pixels": amplitudes,
                "source_response_center": response,
                "center_contrast_estimate": center_estimate,
                "center_relative_contrast_error": (
                    center_estimate - float(task["contrast"])) / float(task["contrast"]),
                "positive_minus_baseline_throughput": (
                    center_estimate - baseline_center) / float(task["contrast"]),
            }
        modes[str(mode)] = method_records
    result = {
        "task": task, "site": site, "source": stage.fingerprint(source),
        "methods_by_mode": modes,
        "annular_oracle_maximum_error": maximum_error,
        "annular_boundary_fallbacks": fallbacks,
        "response_fidelity": fidelity_records,
        "fit_diagnostics": diagnostics,
        "elapsed_seconds": time.monotonic() - started,
    }
    stage.write_json(directory / "measurements.json", result)
    products = [stage.fingerprint(directory / name) for name in
                ("measurements.json", "analysis.log", "command.json", "annular_verification.json")]
    stage.write_json(complete, {"status": "complete", "products": products})
    return str(directory / "measurements.json")


def validation_bootstrap(measurements: list[dict[str, object]], method: str,
                         target: float) -> dict[str, float]:
    """Bootstrap the paired mean-SNR change within each radius."""
    rows = []
    for radius in preparation.RADII:
        chosen = [record for record in measurements
                  if float(record["site"]["nominal_radius"]) == float(radius) and
                  float(record["task"]["target_source_snr"]) == float(target)]
        chosen.sort(key=lambda record: int(record["site"]["role_index"]))
        stage.require(len(chosen) == preparation.ROLE_COUNTS["validation"],
                      f"validation bootstrap lost radius {radius:g}")
        rows.append([float(record["methods_by_mode"][str(policy_module.PRIMARY_MODE)][method][
                                   "search_score"]) -
                     float(record["methods_by_mode"][str(policy_module.PRIMARY_MODE)][
                         policy_module.PRIMARY_COMPARATOR]["search_score"])
                     for record in chosen])
    differences = np.asarray(rows, dtype=np.float64)
    generator = np.random.default_rng(policy_module.VALIDATION_BOOTSTRAP_SEED)
    indices = generator.integers(0, differences.shape[1],
                                 size=(policy_module.BOOTSTRAP_REPLICATES,
                                       differences.shape[0], differences.shape[1]))
    sampled = np.take_along_axis(differences[None, :, :], indices, axis=2)
    means = np.mean(sampled, axis=(1, 2))
    result = policy_module.quantiles(means)
    result["observed_mean"] = float(np.mean(differences))
    result["paired_sites"] = int(differences.size)
    return result


def summarize(root: Path, protocol: dict[str, object], policy: dict[str, object],
              measurements: list[dict[str, object]]) -> None:
    """Write all-mode validation tables and evaluate the preregistered gate."""
    rows = []
    for mode in stage.MODES:
        for radius in preparation.RADII:
            for target in preparation.TARGET_SNRS:
                chosen = [record for record in measurements
                          if float(record["site"]["nominal_radius"]) == float(radius) and
                          float(record["task"]["target_source_snr"]) == float(target)]
                stage.require(len(chosen) == preparation.ROLE_COUNTS["validation"],
                              "validation summary lost paired sites")
                for method in METHODS:
                    values = [record["methods_by_mode"][str(mode)][method] for record in chosen]
                    threshold = float(policy["thresholds"][radius_key(radius)][str(mode)][method])
                    rows.append({
                        "mode": mode, "radius": radius, "target_source_snr": target,
                        "method": method, "sites": len(values), "threshold": threshold,
                        "recoveries": sum(float(value["search_score"]) > threshold
                                          for value in values),
                        "mean_search_snr": float(np.mean([value["search_score"] for value in values])),
                        "median_search_snr": float(np.median([value["search_score"] for value in values])),
                        "mean_center_snr": float(np.mean([value["center_snr"] for value in values])),
                        "mean_localization_error_pixels": float(np.mean([
                            value["localization_error_pixels"] for value in values])),
                        "median_relative_contrast_error": float(np.median([
                            value["center_relative_contrast_error"] for value in values])),
                        "relative_contrast_error_scatter": float(np.std([
                            value["center_relative_contrast_error"] for value in values], ddof=1)),
                        "mean_positive_minus_baseline_throughput": float(np.mean([
                            value["positive_minus_baseline_throughput"] for value in values])),
                    })
    heldout_result = verify_heldout(root)
    heldout_primary = heldout_result["sensitivity"][str(policy_module.PRIMARY_MODE)]
    primary = [row for row in rows if row["mode"] == policy_module.PRIMARY_MODE]
    gate = {}
    comparator = policy_module.PRIMARY_COMPARATOR
    for candidate in COVARIANCE_CANDIDATES:
        bootstrap = {str(float(target)): validation_bootstrap(measurements, candidate, target)
                     for target in preparation.TARGET_SNRS}
        recovery = {}
        throughput_radius = {}
        candidate_total = comparator_total = 0
        all_levels_at_least = True
        support_complete = True
        pooled_throughput = {}
        for target in preparation.TARGET_SNRS:
            candidate_rows = [row for row in primary
                              if row["target_source_snr"] == target and row["method"] == candidate]
            comparator_rows = [row for row in primary
                               if row["target_source_snr"] == target and row["method"] == comparator]
            stage.require(len(candidate_rows) == len(comparator_rows) == len(preparation.RADII),
                          "primary validation rows are incomplete")
            candidate_recovery = sum(int(row["recoveries"]) for row in candidate_rows)
            comparator_recovery = sum(int(row["recoveries"]) for row in comparator_rows)
            recovery[str(float(target))] = {"candidate": candidate_recovery,
                                            "comparator": comparator_recovery}
            candidate_total += candidate_recovery
            comparator_total += comparator_recovery
            all_levels_at_least &= candidate_recovery >= comparator_recovery
            support_complete &= all(int(row["sites"]) == 6 for row in candidate_rows + comparator_rows)
            radius_values = {radius_key(float(row["radius"])):
                             float(row["mean_positive_minus_baseline_throughput"])
                             for row in candidate_rows}
            throughput_radius[str(float(target))] = radius_values
            pooled_throughput[str(float(target))] = float(np.mean(list(radius_values.values())))
        heldout_candidate = int(heldout_primary[candidate]["frozen_exceedances"])
        heldout_comparator = int(heldout_primary[comparator]["frozen_exceedances"])
        heldout_pass = heldout_candidate <= heldout_comparator
        recovery_pass = all_levels_at_least and candidate_total > comparator_total
        snr_success = any(bootstrap[str(float(target))]["p02p5"] > 0
                          for target in (3.0, 5.0))
        no_snr_loss = all(bootstrap[str(float(target))]["p97p5"] >= 0
                          for target in preparation.TARGET_SNRS)
        pooled_range = policy_module.THROUGHPUT_POOLED_RANGE
        radius_range = policy_module.THROUGHPUT_RADIUS_RANGE
        throughput_pass = (all(pooled_range[0] <= value <= pooled_range[1]
                               for value in pooled_throughput.values()) and
                           all(radius_range[0] <= value <= radius_range[1]
                               for target_values in throughput_radius.values()
                               for value in target_values.values()))
        conditions = {"heldout": heldout_pass, "recovery": recovery_pass,
                      "paired_snr": snr_success and no_snr_loss,
                      "throughput": throughput_pass, "common_support": support_complete}
        gate[candidate] = {
            "accepted": all(conditions.values()), "conditions": conditions,
            "heldout_exceedances": {"candidate": heldout_candidate,
                                     "comparator": heldout_comparator},
            "recoveries": recovery, "total_recoveries": {"candidate": candidate_total,
                                                            "comparator": comparator_total},
            "paired_snr_bootstrap": bootstrap,
            "pooled_mean_throughput": pooled_throughput,
            "radius_level_mean_throughput": throughput_radius,
        }

    result = {
        "purpose": "KLIP Stage-E fresh validation under immutable Stage-D policy",
        "policy_changed": False, "selection_performed": False,
        "known_planet_opened": False, "measurements": len(measurements),
        "summary_rows": rows, "acceptance": gate,
        "maximum_annular_oracle_error": max(float(record["annular_oracle_maximum_error"])
                                              for record in measurements),
    }
    stage.write_json(root / "validation_results.json", result)
    with (root / "validation_results.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    lines = ["# KLIP Stage-E fresh validation", "",
             "All values use the immutable Stage-D methods, thresholds, sites, contrasts, and baseline-derived "
             "weights. Mode 200 is primary; no maximization over methods or modes was performed.", "",
             "## Mode-200 pooled comparison", "",
             "| Target | Method | Recovered / 36 | Mean maximum SNR | Mean throughput |",
             "| ---: | :--- | ---: | ---: | ---: |"]
    for target in preparation.TARGET_SNRS:
        for method in METHODS:
            chosen = [row for row in primary
                      if row["target_source_snr"] == target and row["method"] == method]
            lines.append(f"| {target:g} | {method} | {sum(row['recoveries'] for row in chosen)} | "
                         f"{np.mean([row['mean_search_snr'] for row in chosen]):.4f} | "
                         f"{np.mean([row['mean_positive_minus_baseline_throughput'] for row in chosen]):.4f} |")
    lines.extend(["", "## Frozen covariance acceptance gate", "",
                  "| Candidate | Accepted | Held-out (candidate / G3.6) | Recovery totals (candidate / G3.6) | SNR gate | Throughput | Support |",
                  "| :--- | :---: | :---: | :---: | :---: | :---: | :---: |"])
    for candidate in COVARIANCE_CANDIDATES:
        item = gate[candidate]
        lines.append(f"| {candidate} | {'yes' if item['accepted'] else 'no'} | "
                     f"{item['heldout_exceedances']['candidate']} / {item['heldout_exceedances']['comparator']} | "
                     f"{item['total_recoveries']['candidate']} / {item['total_recoveries']['comparator']} | "
                     f"{'pass' if item['conditions']['paired_snr'] else 'fail'} | "
                     f"{'pass' if item['conditions']['throughput'] else 'fail'} | "
                     f"{'pass' if item['conditions']['common_support'] else 'fail'} |")
    lines.extend(["", "The other seven KL modes are retained in the CSV and JSON as separate generalization "
                  "checks. They are not combined into a maximum statistic."])
    (root / "validation_results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def analyze(root: Path, workers: int) -> None:
    """Analyze all fresh validation positives and write the immutable result."""
    protocol, policy, _ = enable(root)
    verify_models(root)
    verify_heldout(root)
    for task in protocol["validation_tasks_unopened"]:
        stage.verify(read(root / "validation_reductions" / str(task["name"]) /
                          "complete.json")["products"])
    if (root / "validation_complete.json").exists():
        receipt = read(root / "validation_complete.json")
        stage.verify(receipt["products"])
        for path in sorted((root / "validation_analysis").glob("*/complete.json")):
            stage.verify(read(path)["products"])
        print(root / "validation_results.md", flush=True)
        return
    analysis_root = root / "validation_analysis"
    analysis_root.mkdir(exist_ok=True)
    tasks = [str(task["name"]) for task in protocol["validation_tasks_unopened"]]
    paths = []
    stage.write_json(root / "state.json", {
        "status": "validation_analyzing", "development_reductions": len(protocol["development_tasks"]),
        "validation_reductions": len(tasks), "heldout_nulls_opened": True,
        "validation_products_opened": True, "positive_analysis_started": True})
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(analyze_task, str(root), name): name for name in tasks}
        for completed, future in enumerate(as_completed(futures), start=1):
            paths.append(Path(future.result()))
            print(f"validation analysis {completed}/{len(tasks)}: {futures[future]}", flush=True)
    measurements = [read(path) for path in sorted(paths)]
    summarize(root, protocol, policy, measurements)
    products = [stage.fingerprint(root / name) for name in
                ("validation_results.json", "validation_results.csv", "validation_results.md")]
    products.extend(stage.fingerprint(root / "validation_analysis" / str(task["name"]) /
                                      "complete.json")
                    for task in protocol["validation_tasks_unopened"])
    stage.write_json(root / "validation_complete.json", {
        "status": "complete", "policy_unchanged": True,
        "known_planet_opened": False, "products": products})
    stage.write_json(root / "state.json", {
        "status": "validation_complete", "development_reductions": len(protocol["development_tasks"]),
        "validation_reductions": len(tasks), "heldout_nulls_opened": True,
        "validation_products_opened": True, "positive_analysis_started": True,
        "known_planet_opened": False})
    print(root / "validation_results.md", flush=True)


def run(root: Path, workers: int) -> None:
    """Run score-blind models, held-out nulls, validation reductions, and analysis in order."""
    models(root, workers)
    heldout(root, workers)
    reduce(root)
    analyze(root, workers)


def check() -> None:
    """Check selected mappings, bootstrap pairing, and acceptance comparisons."""
    stage.require(tuple(METHODS) == policy_module.SELECTED_METHODS and
                  set(WEIGHT_METHODS).issubset(METHODS) and
                  set(COVARIANCE_CANDIDATES).issubset(WEIGHT_METHODS) and
                  len(METHODS) == 7, "Stage-E selected-method mapping changed")
    differences = np.arange(36, dtype=np.float64).reshape(6, 6)
    generator = np.random.default_rng(policy_module.VALIDATION_BOOTSTRAP_SEED)
    indices = generator.integers(0, 6, size=(policy_module.BOOTSTRAP_REPLICATES, 6, 6))
    sampled = np.take_along_axis(differences[None, :, :], indices, axis=2)
    stage.require(sampled.shape == (policy_module.BOOTSTRAP_REPLICATES, 6, 6) and
                  np.all(np.isfinite(np.mean(sampled, axis=(1, 2)))),
                  "radius-stratified bootstrap changed")
    stage.require(not (2 > 2) and (3 > 2), "strict threshold comparison changed")
    template = np.arange(SUPPORT * SUPPORT, dtype=np.float64).reshape(SUPPORT, SUPPORT) + 1
    mask = np.ones(template.shape, dtype=bool)
    mask[0, 0] = False
    fitted = {"support": mask, "scale": np.ones(template.shape),
              "raw_detail": {"covariance": np.eye(mask.size)},
              "radial_detail": {"covariance": np.eye(mask.size)}}
    fidelity = response_fidelity(template, template, fitted)
    stage.require(all(np.isclose(value["cosine"], 1) and
                      np.isclose(value["projection_scale"], 1) and
                      value["best_scaled_relative_residual"] <= 1e-14
                      for value in fidelity.values()),
                  "masked-support response fidelity changed")
    print("KLIP Stage-E validation checks passed", flush=True)


def parser() -> argparse.ArgumentParser:
    """Build the Stage-E command-line parser."""
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="action", required=True)
    subparsers.add_parser("check")
    for action in ("models", "heldout", "reduce", "analyze", "run"):
        current = subparsers.add_parser(action)
        current.add_argument("root", type=Path)
        if action in {"models", "heldout", "analyze", "run"}:
            current.add_argument("--workers", type=int, default=4)
    return result


def main() -> None:
    """Dispatch the Stage-E action."""
    arguments = parser().parse_args()
    if arguments.action == "check":
        check()
        return
    root = arguments.root.resolve()
    ensure_frozen_runner(root)
    preparation.load_experiment(root)
    workers = int(getattr(arguments, "workers", 1))
    stage.require(workers >= 1, "worker count must be positive")
    if arguments.action == "models":
        models(root, workers)
    elif arguments.action == "heldout":
        heldout(root, workers)
    elif arguments.action == "reduce":
        reduce(root)
    elif arguments.action == "analyze":
        analyze(root, workers)
    else:
        run(root, workers)


if __name__ == "__main__":
    main()
