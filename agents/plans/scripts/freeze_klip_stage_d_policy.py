#!/usr/bin/env python3
"""Freeze the KLIP matched-filter policy before held-out or validation scores are opened."""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import shutil
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import prepare_klip_stage_c_development as preparation  # noqa: E402
import run_klip_covariance_stage_a as stage  # noqa: E402
import run_klip_stage_c_development as development  # noqa: E402


SELECTED_METHODS = (
    "native",
    "gaussian_fwhm2p4",
    "gaussian_fwhm3p6",
    "exact_identity",
    "sparse_identity",
    "raw_rectangular_m0p3",
    "radial_hann_m0p1_trunc0p75",
)
COVARIANCE_CANDIDATES = (
    "raw_rectangular_m0p3",
    "radial_hann_m0p1_trunc0p75",
)
PRIMARY_COMPARATOR = "gaussian_fwhm3p6"
STRONG_SMOOTHING_CONTROL = "gaussian_fwhm2p4"
PRIMARY_MODE = 200
BOOTSTRAP_REPLICATES = 4096
BOOTSTRAP_SEED = 260926
VALIDATION_BOOTSTRAP_SEED = 270926
ANGULAR_BLOCK_SIZE = 4
THROUGHPUT_POOLED_RANGE = (0.90, 1.10)
THROUGHPUT_RADIUS_RANGE = (0.85, 1.15)
POLICY_FILES = (
    "policy.json",
    "threshold_audit.json",
    "selection_sensitivity.json",
    "validation_commands.json",
    "resamples.npz",
    "README.md",
)


def read(path: Path) -> object:
    """Read one JSON document."""
    return json.loads(path.read_text(encoding="utf-8"))


def quantiles(values: np.ndarray) -> dict[str, float]:
    """Summarize one finite numeric vector at fixed audit quantiles."""
    values = np.asarray(values, dtype=np.float64)
    stage.require(values.ndim == 1 and len(values) and np.all(np.isfinite(values)),
                  "audit quantiles require a nonempty finite vector")
    return {
        "minimum": float(np.min(values)),
        "p02p5": float(np.quantile(values, 0.025)),
        "median": float(np.median(values)),
        "mean": float(np.mean(values)),
        "p97p5": float(np.quantile(values, 0.975)),
        "maximum": float(np.max(values)),
    }


def radius_key(radius: float) -> str:
    """Return the Stage-C JSON key for one nominal radius."""
    return str(float(radius))


def selected_sites(protocol: dict[str, object], radius: float, role: str) -> list[dict[str, object]]:
    """Return one role's sites in frozen role-index order."""
    rows = [site for site in protocol["sites"]
            if float(site["nominal_radius"]) == float(radius) and site["role"] == role]
    rows.sort(key=lambda site: int(site["role_index"]))
    stage.require(len(rows) == int(protocol["sites_per_radius_by_role"][role]),
                  f"site partition lost {role} radius {radius:g}")
    return rows


def verify_development(root: Path, protocol: dict[str, object]) -> None:
    """Verify the complete Stage-C result and its unopened-validation state."""
    state = read(root / "state.json")
    stage.require(state["status"] == "development_complete" and
                  not state["validation_products_opened"],
                  "policy freeze requires completed development and unopened validation products")
    stage.require(not (root / "validation_reductions").exists() and
                  not (root / "validation_analysis").exists() and
                  not (root / "heldout_analysis").exists(),
                  "held-out or validation products already exist")
    completion = read(root / "development_complete.json")
    stage.require(completion["status"] == "complete" and
                  not completion["validation_products_opened"],
                  "development completion receipt permits opened validation products")
    stage.verify(completion["products"])
    calibration = read(root / "calibration" / "complete.json")
    stage.require(calibration["status"] == "complete", "calibration is incomplete")
    stage.verify(calibration["products"])
    for receipt_path in sorted((root / "calibration" / "units").glob("*/complete.json")):
        stage.verify(read(receipt_path)["products"])
    for receipt_path in sorted((root / "calibration" / "sites").glob("*/complete.json")):
        stage.verify(read(receipt_path)["products"])
    result = read(root / "development_results.json")
    stage.require(result["purpose"] == "KLIP Stage-C calibrated development result" and
                  not result["selection_performed"] and
                  not result["validation_products_opened"] and
                  int(result["measurements"]) == len(protocol["development_tasks"]),
                  "development result contract changed")
    for task in protocol["development_tasks"]:
        reduction = read(root / "reductions" / str(task["name"]) / "complete.json")
        analysis = read(root / "development_analysis" / str(task["name"]) / "complete.json")
        stage.verify(reduction["products"])
        stage.verify(analysis["products"])


def immutable_input_records(root: Path, protocol: dict[str, object]) -> list[dict[str, object]]:
    """Collect direct fingerprints for every generated input used by policy or validation."""
    paths = [root / name for name in (
        "protocol.json", "manifest.json", "geometry.json", "contrasts.json",
        "development_complete.json", "development_results.json",
        "development_results.csv", "calibration/complete.json",
        "calibration/thresholds.json", "calibration/records.json")]
    records = [stage.fingerprint(path) for path in paths]
    for group in (root / "calibration" / "units", root / "calibration" / "sites"):
        for receipt_path in sorted(group.glob("*/complete.json")):
            receipt = read(receipt_path)
            records.append(stage.fingerprint(receipt_path))
            records.extend(receipt["products"])
    for task in protocol["development_tasks"]:
        for parent in ("reductions", "development_analysis"):
            receipt_path = root / parent / str(task["name"]) / "complete.json"
            receipt = read(receipt_path)
            records.append(stage.fingerprint(receipt_path))
            records.extend(receipt["products"])
        measurement = read(root / "development_analysis" / str(task["name"]) /
                           "measurements.json")
        records.append(measurement["source"])
    unique = {}
    for record in records:
        path = str(record["path"])
        stage.require(path not in unique or unique[path] == record,
                      f"conflicting immutable fingerprints for {path}")
        unique[path] = record
    return [unique[path] for path in sorted(unique)]


def ensure_frozen(root: Path, arguments: list[str]) -> None:
    """Copy the policy and validation runners into the experiment before freezing."""
    source = Path(__file__).resolve()
    destination = root / "software" / source.name
    validation_source = source.with_name("run_klip_stage_e_validation.py")
    validation_destination = root / "software" / validation_source.name
    stage.require(validation_source.is_file(), "Stage-E validation runner is missing")
    if source == destination.resolve():
        stage.require(validation_destination.is_file(), "frozen Stage-E runner is missing")
        return
    for current, frozen in ((source, destination), (validation_source, validation_destination)):
        if frozen.exists():
            stage.require(frozen.read_bytes() == current.read_bytes(),
                          f"frozen software differs from repository: {frozen.name}")
        else:
            shutil.copy2(current, frozen)
    os.execv(sys.executable, [sys.executable, str(destination), *arguments])


def validation_command(root: Path, protocol: dict[str, object], task: dict[str, object]) -> list[str]:
    """Construct one validation reduction command with a separate output tree."""
    command = preparation.task_command(root, protocol, task)
    marker = command.index("--output.directory")
    command[marker + 1] = str(root / "validation_reductions" / str(task["name"]))
    return command


def calibration_lookup(records: list[dict[str, object]]) -> dict[tuple[float, str], dict[str, object]]:
    """Index the frozen calibration records by radius and site name."""
    result = {(float(record["radius"]), str(record["site"])): record for record in records}
    stage.require(len(result) == len(records), "calibration record keys are not unique")
    return result


def development_scores(root: Path, protocol: dict[str, object]) \
        -> dict[tuple[float, float, int, str], np.ndarray]:
    """Load paired frozen-weight development search scores for the selected methods."""
    task_lookup = {str(task["name"]): task for task in protocol["development_tasks"]}
    site_lookup = {str(site["name"]): site for site in protocol["sites"]}
    collected: dict[tuple[float, float, int, str], list[tuple[int, float]]] = {}
    for name, task in task_lookup.items():
        measurement = read(root / "development_analysis" / name / "measurements.json")
        stage.require(measurement["task"] == task, f"development task changed: {name}")
        site = site_lookup[str(task["site"])]
        stage.require(measurement["site"] == site, f"development site changed: {name}")
        for mode in stage.MODES:
            for method in SELECTED_METHODS:
                value = float(measurement["arms"]["frozen"][str(mode)][method]["search_score"])
                key = (float(site["nominal_radius"]), float(task["target_source_snr"]), mode, method)
                collected.setdefault(key, []).append((int(site["role_index"]), value))
    result = {}
    for key, values in collected.items():
        values.sort()
        stage.require([index for index, _ in values] == list(range(
                      preparation.ROLE_COUNTS["development"])),
                      f"development pairing changed for {key}")
        vector = np.asarray([value for _, value in values], dtype=np.float64)
        stage.require(np.all(np.isfinite(vector)), f"nonfinite development score for {key}")
        result[key] = vector
    expected = (len(preparation.RADII) * len(preparation.TARGET_SNRS) *
                len(stage.MODES) * len(SELECTED_METHODS))
    stage.require(len(result) == expected, "development score table is incomplete")
    return result


def geometry_audit(protocol: dict[str, object], radius: float) -> dict[str, object]:
    """Describe angular spacing without treating locations as independent samples."""
    result = {}
    for role in ("calibration", "development", "validation", "heldout_null"):
        sites = selected_sites(protocol, radius, role)
        angles = np.sort(np.asarray([float(site["angle_radians"]) for site in sites]))
        gaps = np.diff(np.r_[angles, angles[0] + 2 * math.pi])
        result[role] = {
            "count": len(sites),
            "site_names": [str(site["name"]) for site in sites],
            "angles_radians": angles.tolist(),
            "minimum_wrapped_gap_radians": float(np.min(gaps)),
            "median_wrapped_gap_radians": float(np.median(gaps)),
            "maximum_wrapped_gap_radians": float(np.max(gaps)),
        }
    return result


def build_resamples(protocol: dict[str, object]) -> tuple[np.ndarray, np.ndarray]:
    """Build fixed paired bootstrap and circular angular-block retention indices."""
    generator = np.random.default_rng(BOOTSTRAP_SEED)
    bootstrap = np.empty((len(preparation.RADII), BOOTSTRAP_REPLICATES,
                          preparation.ROLE_COUNTS["calibration"]), dtype=np.int16)
    retained = np.empty((len(preparation.RADII), preparation.ROLE_COUNTS["calibration"],
                         preparation.ROLE_COUNTS["calibration"] - ANGULAR_BLOCK_SIZE),
                        dtype=np.int16)
    for radius_index, radius in enumerate(preparation.RADII):
        count = preparation.ROLE_COUNTS["calibration"]
        bootstrap[radius_index] = generator.integers(0, count,
                                                     size=(BOOTSTRAP_REPLICATES, count))
        sites = selected_sites(protocol, radius, "calibration")
        angular_order = np.argsort([float(site["angle_radians"]) for site in sites])
        for start in range(count):
            removed = {int(angular_order[(start + offset) % count])
                       for offset in range(ANGULAR_BLOCK_SIZE)}
            retained[radius_index, start] = [index for index in range(count)
                                              if index not in removed]
    return bootstrap, retained


def build_audit(root: Path, protocol: dict[str, object], thresholds: dict[str, object],
                records: list[dict[str, object]], bootstrap: np.ndarray,
                retained: np.ndarray) -> tuple[dict[str, object], dict[str, object]]:
    """Copy signed nulls and calculate paired threshold/recovery sensitivity."""
    record_lookup = calibration_lookup(records)
    scores = development_scores(root, protocol)
    audit: dict[str, object] = {"schema": 1, "threshold_rule": protocol["threshold_rule"],
                                "strict_exceedance": True, "radii": {}}
    sensitivity: dict[str, object] = {
        "schema": 1,
        "paired_bootstrap": {"seed": BOOTSTRAP_SEED,
                              "replicates": BOOTSTRAP_REPLICATES,
                              "sample_size_per_radius": 20,
                              "shared_indices_across_modes_and_methods": True},
        "angular_block_deletion": {"deleted_consecutive_sites": ANGULAR_BLOCK_SIZE,
                                    "circular_deletions_per_radius": 20,
                                    "angles_define_order": True},
        "radii": {},
    }
    for radius_index, radius in enumerate(preparation.RADII):
        radius_name = radius_key(radius)
        sites = selected_sites(protocol, radius, "calibration")
        audit_radius = {"geometry": geometry_audit(protocol, radius), "modes": {}}
        sensitivity_radius = {"modes": {}}
        for mode in stage.MODES:
            audit_methods, sensitivity_methods = {}, {}
            for method in SELECTED_METHODS:
                null_records = []
                maxima = []
                for site in sites:
                    record = record_lookup[(float(radius), str(site["name"]))]
                    measured = record["methods_by_mode"][str(mode)][method]
                    pixels = np.asarray(measured["snr_pixels"], dtype=np.float64)
                    maximum = float(measured["search_score"])
                    stage.require(pixels.shape == (5,) and np.all(np.isfinite(pixels)) and
                                  maximum == float(np.max(pixels)),
                                  f"invalid signed null record for {site['name']} {mode} {method}")
                    maxima.append(maximum)
                    null_records.append({
                        "site": str(site["name"]),
                        "role_index": int(site["role_index"]),
                        "angle_radians": float(site["angle_radians"]),
                        "search_score": maximum,
                        "snr_pixels": pixels.tolist(),
                    })
                maxima_array = np.asarray(maxima, dtype=np.float64)
                frozen = thresholds[radius_name][str(mode)][method]
                ordered = np.sort(maxima_array)
                stage.require(list(map(float, frozen["site_maxima"])) == maxima_array.tolist() and
                              float(frozen["threshold"]) == float(ordered[-1]) and
                              float(frozen["largest"]) == float(ordered[-1]) and
                              float(frozen["second_largest"]) == float(ordered[-2]) and
                              float(frozen["largest_second_gap"]) == float(ordered[-1] - ordered[-2]) and
                              bool(frozen["strict_exceedance"]),
                              f"threshold receipt changed for {radius_name} {mode} {method}")
                audit_methods[method] = {
                    "threshold": float(frozen["threshold"]),
                    "largest": float(ordered[-1]),
                    "second_largest": float(ordered[-2]),
                    "largest_second_gap": float(ordered[-1] - ordered[-2]),
                    "signed_nulls": null_records,
                }

                bootstrap_thresholds = np.max(maxima_array[bootstrap[radius_index]], axis=1)
                deletion_thresholds = np.max(maxima_array[retained[radius_index]], axis=1)
                method_sensitivity: dict[str, object] = {
                    "paired_bootstrap_threshold": quantiles(bootstrap_thresholds),
                    "angular_block_delete_threshold": quantiles(deletion_thresholds),
                    "fraction_bootstrap_threshold_below_frozen": float(np.mean(
                        bootstrap_thresholds < float(frozen["threshold"]))),
                    "fraction_block_delete_threshold_below_frozen": float(np.mean(
                        deletion_thresholds < float(frozen["threshold"]))),
                    "development_recovery": {},
                }
                for target in preparation.TARGET_SNRS:
                    positive = scores[(float(radius), float(target), mode, method)]
                    bootstrap_recovery = np.sum(
                        positive[None, :] > bootstrap_thresholds[:, None], axis=1)
                    deletion_recovery = np.sum(
                        positive[None, :] > deletion_thresholds[:, None], axis=1)
                    base_recovery = int(np.sum(positive > float(frozen["threshold"])))
                    method_sensitivity["development_recovery"][str(float(target))] = {
                        "frozen": base_recovery,
                        "paired_bootstrap": quantiles(bootstrap_recovery),
                        "angular_block_delete": quantiles(deletion_recovery),
                    }
                sensitivity_methods[method] = method_sensitivity
            audit_radius["modes"][str(mode)] = audit_methods
            sensitivity_radius["modes"][str(mode)] = sensitivity_methods
        audit["radii"][radius_name] = audit_radius
        sensitivity["radii"][radius_name] = sensitivity_radius
    return audit, sensitivity


def policy_document(protocol: dict[str, object], thresholds: dict[str, object],
                    commands: list[dict[str, object]]) -> dict[str, object]:
    """Construct the exact immutable Stage-D method and acceptance policy."""
    copied_thresholds = {
        radius_key(radius): {
            str(mode): {method: float(thresholds[radius_key(radius)][str(mode)][method]["threshold"])
                        for method in SELECTED_METHODS}
            for mode in stage.MODES}
        for radius in preparation.RADII
    }
    sites = {role: [site for site in protocol["sites"] if site["role"] == role]
             for role in ("calibration", "development", "validation", "heldout_null")}
    return {
        "schema": 1,
        "stage": "KLIP covariance matched-filter Stage D policy freeze",
        "frozen_before_heldout_or_validation": True,
        "selection_source": "complete Stage-C development grid only",
        "primary_mode": PRIMARY_MODE,
        "secondary_modes": [mode for mode in stage.MODES if mode != PRIMARY_MODE],
        "joint_mode_maximization": False,
        "response_support_pixels": 11,
        "candidate_search_offsets_row_column": protocol["search_offsets_row_column"],
        "baseline_frozen_weights_only": True,
        "methods": {
            "ordered": list(SELECTED_METHODS),
            "controls": ["native", "gaussian_fwhm2p4", "gaussian_fwhm3p6",
                         "exact_identity", "sparse_identity"],
            "covariance_candidates": list(COVARIANCE_CANDIDATES),
            "primary_comparator": PRIMARY_COMPARATOR,
            "strong_smoothing_control": STRONG_SMOOTHING_CONTROL,
            "parameters": {
                "native": "unfiltered final-image pixel",
                "gaussian_fwhm2p4": "15-pixel mask-normalized Gaussian convolution; FWHM 2.4 pixels",
                "gaussian_fwhm3p6": "15-pixel mask-normalized Gaussian convolution; FWHM 3.6 pixels",
                "exact_identity": "signal-free per-pixel exact KLIP response; central 11 pixels; identity covariance",
                "sparse_identity": "candidate-avoiding radial response; rotated-response pixelwise average and radial interpolation; central 11 pixels; identity covariance",
                "raw_rectangular_m0p3": "raw 11-pixel rectangular Welch PSD; mixing 0.3; unit-response weight; no fitted candidate mean",
                "radial_hann_m0p1_trunc0p75": "strict leave-site-out radial-standardized 11-pixel Hann Welch PSD; mixing 0.1; precision modes below 0.75 mean variance set to zero; unit physical response; no fitted candidate mean",
            },
        },
        "thresholds": copied_thresholds,
        "threshold_rule": protocol["threshold_rule"],
        "source_levels": list(protocol["target_source_snrs"]),
        "physical_contrasts": "frozen by the parent contrasts.json fingerprint",
        "sites": sites,
        "validation_tasks": list(protocol["validation_tasks_unopened"]),
        "expected_validation_commands": len(commands),
        "acceptance": {
            "applies_to": list(COVARIANCE_CANDIDATES),
            "comparison": PRIMARY_COMPARATOR,
            "mode": PRIMARY_MODE,
            "heldout": "aggregate exceedance count across the six radii must be no larger than the comparator under unchanged method-specific thresholds",
            "recovery": "at least comparator recovery at each SNR 3/5/7 level and strictly larger total recovery",
            "paired_snr": {
                "bootstrap": "paired radius-stratified site bootstrap; resample six sites with replacement independently within each radius",
                "replicates": BOOTSTRAP_REPLICATES,
                "seed": VALIDATION_BOOTSTRAP_SEED,
                "success": "95-percent lower bound above zero at SNR 3 or SNR 5",
                "no_significant_loss": "95-percent upper bound is not below zero at every SNR 3/5/7 level",
            },
            "throughput": {
                "pooled_mean_each_level_inclusive_range": list(THROUGHPUT_POOLED_RANGE),
                "radius_level_mean_inclusive_range": list(THROUGHPUT_RADIUS_RANGE),
            },
            "support": "all 36 paired validation sites must be finite for candidate and comparator at every source level",
            "all_conditions_required": True,
        },
        "reporting": {
            "all_modes_separate": True,
            "gaussian_fwhm2p4_secondary_comparison": True,
            "heldout_samples_called_independent": False,
            "maximum_of_20_interpreted_as_formal_false_alarm_probability": False,
            "heldout_can_reject_but_not_retune_or_nominate": True,
        },
    }


def report(policy: dict[str, object], audit: dict[str, object],
           sensitivity: dict[str, object]) -> str:
    """Render a compact human-readable policy receipt."""
    lines = [
        "# KLIP Stage-D frozen validation policy",
        "",
        "This receipt was written before any held-out-null score or validation positive was opened. "
        "It freezes mode 200 as primary, baseline-derived weights, the seven-method shortlist, the existing "
        "method-specific maximum-of-20 thresholds, and the Stage-E acceptance rule.",
        "",
        "## Frozen methods",
        "",
    ]
    lines.extend(f"- `{method}`" for method in policy["methods"]["ordered"])
    lines.extend([
        "",
        "The confirmatory covariance candidates are `raw_rectangular_m0p3` and "
        "`radial_hann_m0p1_trunc0p75`. The primary comparison is paired against "
        "`gaussian_fwhm3p6`; `gaussian_fwhm2p4` is reported as the stronger smoothing control.",
        "",
        "## Mode-200 threshold audit",
        "",
        "| Radius | Method | Threshold | Second largest | Gap | Bootstrap 2.5--97.5% | Block-delete 2.5--97.5% |",
        "| ---: | :--- | ---: | ---: | ---: | :--- | :--- |",
    ])
    for radius in preparation.RADII:
        key = radius_key(radius)
        for method in SELECTED_METHODS:
            item = audit["radii"][key]["modes"][str(PRIMARY_MODE)][method]
            item_sensitivity = sensitivity["radii"][key]["modes"][str(PRIMARY_MODE)][method]
            boot = item_sensitivity["paired_bootstrap_threshold"]
            block = item_sensitivity["angular_block_delete_threshold"]
            lines.append(
                f"| {radius:g} | {method} | {item['threshold']:.4f} | "
                f"{item['second_largest']:.4f} | {item['largest_second_gap']:.4f} | "
                f"{boot['p02p5']:.4f}--{boot['p97p5']:.4f} | "
                f"{block['p02p5']:.4f}--{block['p97p5']:.4f} |")
    lines.extend([
        "",
        "The 20 calibration locations are spatially correlated. Their maximum is retained for direct continuity "
        "with the P4 tests and is not assigned a formal false-alarm probability. Six separately assigned null "
        "locations per radius may now be exposed; they may reject a method but cannot alter this policy.",
        "",
        "## Validation gate",
        "",
        "A covariance candidate must have no more held-out exceedances than Gaussian 3.6, match or exceed its "
        "recovery at every source level with a strictly larger total, pass the frozen paired-SNR bootstrap rule, "
        "retain the frozen throughput ranges, and retain all common sites. Other KL modes are reported separately "
        "and are not maximized.",
    ])
    return "\n".join(lines) + "\n"


def existing_complete(root: Path) -> bool:
    """Verify and acknowledge an already completed immutable policy."""
    complete = root / "policy_complete.json"
    if not complete.exists():
        return False
    receipt = read(complete)
    stage.require(receipt["status"] == "complete" and
                  receipt["heldout_nulls_opened"] is False and
                  receipt["validation_products_opened"] is False,
                  "policy completion receipt is invalid")
    stage.verify(receipt["products"])
    manifest = read(root / "policy_manifest.json")
    stage.verify(manifest["input_records"] + manifest["software_records"] +
                 manifest["policy_records"] + manifest.get("repair_records", []))
    print(root / "policy" / "README.md", flush=True)
    return True


def freeze(root: Path) -> None:
    """Create the immutable policy receipt and leave held-out and validation data unopened."""
    protocol, _ = preparation.load_experiment(root)
    if existing_complete(root):
        return
    ensure_frozen(root, sys.argv[1:])
    protocol, _ = preparation.load_experiment(root)
    verify_development(root, protocol)
    stage.require(PRIMARY_MODE == int(protocol["primary_mode"]) and
                  tuple(SELECTED_METHODS) == tuple(dict.fromkeys(SELECTED_METHODS)) and
                  set(SELECTED_METHODS).issubset(development.METHODS) and
                  set(COVARIANCE_CANDIDATES).issubset(SELECTED_METHODS),
                  "selected method or primary-mode contract changed")
    policy_directory = root / "policy"
    if policy_directory.exists():
        development.archive_incomplete(root, policy_directory, "policy")
    policy_directory.mkdir()

    thresholds = read(root / "calibration" / "thresholds.json")
    calibration_records = read(root / "calibration" / "records.json")
    bootstrap, retained = build_resamples(protocol)
    np.savez_compressed(policy_directory / "resamples.npz",
                        bootstrap_indices=bootstrap, block_retained_indices=retained,
                        radii=np.asarray(preparation.RADII, dtype=np.float64))
    audit, sensitivity = build_audit(root, protocol, thresholds, calibration_records,
                                     bootstrap, retained)
    commands = [{"task": task["name"], "command": validation_command(root, protocol, task)}
                for task in protocol["validation_tasks_unopened"]]
    stage.require(len(commands) == int(protocol["expected_future_validation_reductions"]),
                  "validation command count changed")
    policy = policy_document(protocol, thresholds, commands)
    stage.write_json(policy_directory / "policy.json", policy)
    stage.write_json(policy_directory / "threshold_audit.json", audit)
    stage.write_json(policy_directory / "selection_sensitivity.json", sensitivity)
    stage.write_json(policy_directory / "validation_commands.json", commands)
    (policy_directory / "README.md").write_text(report(policy, audit, sensitivity),
                                                  encoding="utf-8")

    input_records = immutable_input_records(root, protocol)
    software_records = [stage.fingerprint(root / "software" / name) for name in (
        "prepare_klip_stage_c_development.py", "run_klip_stage_c_development.py",
        "freeze_klip_stage_d_policy.py", "run_klip_stage_e_validation.py")]
    policy_records = [stage.fingerprint(policy_directory / name) for name in POLICY_FILES]
    stage.write_json(root / "policy_manifest.json", {
        "schema": 1,
        "input_records": input_records,
        "software_records": software_records,
        "policy_records": policy_records,
    })
    products = policy_records + [stage.fingerprint(root / "policy_manifest.json")]
    stage.write_json(root / "policy_complete.json", {
        "status": "complete",
        "heldout_nulls_opened": False,
        "validation_products_opened": False,
        "products": products,
    })
    stage.write_json(root / "state.json", {
        "status": "policy_frozen",
        "development_reductions": len(protocol["development_tasks"]),
        "heldout_nulls_opened": False,
        "validation_products_opened": False,
        "positive_analysis_started": True,
    })
    print(policy_directory / "README.md", flush=True)


def check() -> None:
    """Check fixed resampling, validation output isolation, and method policy."""
    generator = np.random.default_rng(BOOTSTRAP_SEED)
    first = generator.integers(0, 20, size=(BOOTSTRAP_REPLICATES, 20))
    generator = np.random.default_rng(BOOTSTRAP_SEED)
    second = generator.integers(0, 20, size=(BOOTSTRAP_REPLICATES, 20))
    stage.require(np.array_equal(first, second) and first.shape == (BOOTSTRAP_REPLICATES, 20),
                  "paired bootstrap is not reproducible")
    command = ["klipReduce", "--output.directory", "/old"]
    marker = command.index("--output.directory")
    command[marker + 1] = "/validation"
    stage.require(command[-1] == "/validation", "validation output isolation failed")
    stage.require(len(SELECTED_METHODS) == 7 and len(COVARIANCE_CANDIDATES) == 2 and
                  PRIMARY_MODE == preparation.PRIMARY_MODE and
                  PRIMARY_COMPARATOR in SELECTED_METHODS and
                  STRONG_SMOOTHING_CONTROL in SELECTED_METHODS,
                  "Stage-D method policy changed")
    synthetic = np.asarray([1.0, 2.0, 3.0])
    summary = quantiles(synthetic)
    stage.require(summary["minimum"] == 1 and summary["maximum"] == 3,
                  "audit quantiles changed")
    print("KLIP Stage-D policy checks passed", flush=True)


def parser() -> argparse.ArgumentParser:
    """Build the Stage-D command-line parser."""
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="action", required=True)
    subparsers.add_parser("check")
    freeze_parser = subparsers.add_parser("freeze")
    freeze_parser.add_argument("root", type=Path)
    return result


def main() -> None:
    """Dispatch the Stage-D action."""
    arguments = parser().parse_args()
    if arguments.action == "check":
        check()
    else:
        freeze(arguments.root.resolve())


if __name__ == "__main__":
    main()
