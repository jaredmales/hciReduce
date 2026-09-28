#!/usr/bin/env python3
"""Test shifted and fully regenerated KLIP responses at the planet's subpixel phase."""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time

import numpy as np
from astropy.io import fits
from scipy.ndimage import shift as cubic_shift

SCRIPT_PATH = Path(__file__).resolve()
SCRIPT_DIRECTORY = SCRIPT_PATH.parent
FROZEN_DEPENDENCIES = (SCRIPT_PATH.parents[2] / "software"
                       if len(SCRIPT_PATH.parents) > 2 else Path("/nonexistent"))
if (FROZEN_DEPENDENCIES / "run_klip_stage_g_planet_consistency.py").is_file():
    sys.path.insert(0, str(FROZEN_DEPENDENCIES))
else:
    sys.path.insert(0, str(SCRIPT_DIRECTORY))

import run_klip_stage_g_planet_consistency as stage_g  # noqa: E402


stage = stage_g.stage
development = stage_g.development
validation = stage_g.validation
stage_f = stage_g.stage_f

STAGE_DIRECTORY = "stage_h_subpixel_response"
PRIMARY_MODE = stage_g.PRIMARY_MODE
PRIMARY_ARM = stage_g.PRIMARY_ARM
SITE_COUNT = stage_g.SITE_COUNT
SHIFT_ORDER = 3
SHIFT_MODE = "nearest"
INTERPOLATION_GUARD = stage_f.SUPPORT + 4
MODEL_NAMES = ("integer", "shifted", "paired_exact")
METRIC_NAMES = ("unweighted", "raw_rectangular_m0p3", "radial_hann_m0p1")
VALUE_NAMES = ("cosine", "projection_scale", "best_scaled_relative_residual")


def read(path: Path) -> object:
    """Read one JSON document."""
    return json.loads(path.read_text(encoding="utf-8"))


def output_root(root: Path) -> Path:
    """Return the isolated Stage-H product directory."""
    return root / STAGE_DIRECTORY


def primary_tasks(protocol: dict[str, object]) -> list[dict[str, object]]:
    """Return the ordered Stage-G nominal optimized-phase tasks."""
    tasks = [dict(task) for task in protocol["tasks"]
             if str(task["arm"]) == PRIMARY_ARM]
    tasks.sort(key=lambda task: int(task["site_index"]))
    stage.require(len(tasks) == SITE_COUNT and
                  len({str(task["base_site"]) for task in tasks}) == SITE_COUNT,
                  "Stage-H parent tasks lost a primary-arm site")
    return tasks


def replace_option(command: list[str], option: str, value: str) -> None:
    """Replace the unique value following one two-token command option."""
    matches = [index for index, item in enumerate(command) if item == option]
    stage.require(len(matches) == 1 and matches[0] + 1 < len(command),
                  f"command lacks a unique {option}")
    command[matches[0] + 1] = value


def negative_command(output: Path, task: dict[str, object]) -> list[str]:
    """Turn one frozen Stage-G positive command into its paired negative trial."""
    command = [str(value) for value in task["command"]]
    contrast = float(task["contrast"])
    stage.require(np.isfinite(contrast) and contrast > 0,
                  "Stage-G task contrast is not positive")
    replace_option(command, "--fake.contrast", repr(-contrast))
    replace_option(command, "--output.directory",
                   str(output / "negative" / str(task["name"])))
    return command


def stage_g_boundary(root: Path) -> tuple[dict[str, object], dict[str, object]]:
    """Verify completed Stage G and return its protocol plus Stage-C protocol."""
    stage_g.verify_complete(root)
    parent_output = stage_g.stage_g_root(root)
    protocol = read(parent_output / "protocol.json")
    manifest = read(parent_output / "manifest.json")
    stage.verify(manifest["input_records"] + manifest["software_records"] +
                 manifest["model_receipts"] + manifest.get("repair_records", []))
    stage.require(str(Path(str(protocol["parent_root"])).resolve()) == str(root.resolve()) and
                  int(protocol["primary_mode"]) == PRIMARY_MODE and
                  str(protocol["primary_arm"]) == PRIMARY_ARM,
                  "Stage-G parent boundary changed")
    parent = read(root / "protocol.json")
    stage.require(Path(str(parent["paths"]["baseline"])).is_file(),
                  "Stage-C signal-free baseline is missing")
    return protocol, parent


def response_product_records(parent: dict[str, object]) -> list[dict[str, object]]:
    """Fingerprint the exact-response coordinate, response, and validity products."""
    paths = development.response47.product_paths(Path(str(parent["parent_response"])))
    _, header = fits.getdata(paths["manifest"], header=True)
    stage.require(str(header["KLIP PSF SPATIAL MODEL"]).strip() == "PIXEL_EXACT" and
                  str(header["KLIP PSF RESPONSE METHOD"]).strip() == "refitDifference" and
                  int(header["KLIP PSF STAMP SIZE"]) == 47,
                  "Stage-H parent response contract changed")
    products = [paths["manifest"], paths["coordinates"],
                *paths["responses"], *paths["validities"]]
    stage.require(all(Path(path).is_file() for path in products),
                  "Stage-H exact-response product set is incomplete")
    return [stage.fingerprint(Path(path)) for path in products]


def prepare(root: Path) -> None:
    """Freeze the twelve paired negative reductions and response comparison."""
    root = root.resolve()
    protocol_g, parent = stage_g_boundary(root)
    stage.require(sorted(os.sched_getaffinity(0)) == parent["resources"]["cpu_affinity"],
                  "Stage-H prepare must use the frozen response CPU affinity")
    output = output_root(root)
    stage.require(not output.exists(), f"Stage-H output already exists: {output}")
    output.mkdir()
    (output / "software").mkdir()

    frozen_runner = output / "software" / SCRIPT_PATH.name
    shutil.copy2(SCRIPT_PATH, frozen_runner)
    parent_runner = stage_g.stage_g_root(root) / "software" / stage_g.SCRIPT_PATH.name
    parent_stage_f = stage_f.stage_f_root(root) / "software" / stage_f.SCRIPT_PATH.name
    stage.require(parent_runner.is_file() and parent_stage_f.is_file(),
                  "Stage-H parent frozen software is incomplete")
    for source in (parent_runner, parent_stage_f):
        shutil.copy2(source, output / "software" / source.name)

    tasks = []
    input_records = [
        stage.fingerprint(stage_g.stage_g_root(root) / "complete.json"),
        stage.fingerprint(stage_g.stage_g_root(root) / "protocol.json"),
        stage.fingerprint(stage_g.stage_g_root(root) / "manifest.json"),
        stage.fingerprint(stage_g.stage_g_root(root) / "results.json"),
        stage.fingerprint(Path(str(parent["paths"]["baseline"]))),
    ] + response_product_records(parent)
    for task in primary_tasks(protocol_g):
        name = str(task["name"])
        positive = stage_g.stage_g_root(root) / "reductions" / name / "finim.fits"
        receipt = stage_g.stage_g_root(root) / "reductions" / name / "complete.json"
        stage_g.verify_task_receipt(receipt)
        stage_g.validate_reduction(positive)
        record = {
            "name": name,
            "site_index": int(task["site_index"]),
            "base_site": str(task["base_site"]),
            "source": dict(task["source"]),
            "settings": dict(task["settings"]),
            "contrast": float(task["contrast"]),
            "phase_offset_row_column": list(task["phase_offset_row_column"]),
            "positive": str(positive),
        }
        record["negative_command"] = negative_command(output, task)
        tasks.append(record)
        input_records.extend([stage.fingerprint(positive), stage.fingerprint(receipt)])

    expected_phase = tuple(float(value)
                           for value in protocol_g["phase_offsets_row_column"]["optimized"])
    stage.require(np.allclose(expected_phase, [-0.2770699293779586, 0.4859637006126931],
                              rtol=0, atol=1e-12) and
                  all(np.allclose(task["phase_offset_row_column"], expected_phase,
                                  rtol=0, atol=2e-12) for task in tasks),
                  "Stage-H optimized subpixel phase changed")

    protocol = {
        "schema": 1,
        "stage": "KLIP Stage H subpixel-response closure",
        "purpose": (
            "compare the current integer response and a cubic-shifted integer response "
            "with a full paired KLIP response regenerated at the planet-like subpixel phase"
        ),
        "parent_root": str(root),
        "parent_stage_g": str(stage_g.stage_g_root(root)),
        "primary_mode": PRIMARY_MODE,
        "primary_arm": PRIMARY_ARM,
        "site_count": SITE_COUNT,
        "modes": list(stage.MODES),
        "phase_offset_row_column": list(expected_phase),
        "shift": {
            "implementation": "scipy.ndimage.shift on the 47-pixel integer response",
            "order": SHIFT_ORDER,
            "mode": SHIFT_MODE,
            "prefilter": True,
            "crop_after_shift_pixels": stage_f.SUPPORT,
            "required_central_validity_guard_pixels": INTERPOLATION_GUARD,
        },
        "response_targets": {
            "central": "(positive-negative)/(2*contrast)",
            "positive": "(positive-baseline)/contrast",
            "negative": "(baseline-negative)/contrast",
        },
        "models": list(MODEL_NAMES),
        "metrics": list(METRIC_NAMES),
        "values": list(VALUE_NAMES),
        "primary_endpoints": [
            "mode-200 shifted-minus-integer cosine against the paired central response",
            "mode-200 integer-minus-shifted best-scaled residual against the paired central response",
            "mode-200 reduction in absolute projection-scale error against the paired central response",
        ],
        "secondary_endpoints": [
            "the same comparisons against each one-sided response and in every retained KL mode",
            "paired-exact template fidelity to the positive response as the finite-amplitude ceiling",
        ],
        "interpretation": (
            "cosine is the retained optimal matched-filter SNR in the stated covariance metric; "
            "projection scale is noiseless contrast gain; the paired-exact model is an oracle "
            "response target and does not use annular detection noise"
        ),
        "tasks": tasks,
        "resources": parent["resources"],
        "paths": {
            "klipreduce": parent["paths"]["klipreduce"],
            "baseline": parent["paths"]["baseline"],
            "parent_response": parent["parent_response"],
        },
    }
    stage.write_json(output / "protocol.json", protocol)
    input_records.append(stage.fingerprint(output / "protocol.json"))
    software_records = [
        stage.fingerprint(frozen_runner),
        stage.fingerprint(output / "software" / parent_runner.name),
        stage.fingerprint(output / "software" / parent_stage_f.name),
        stage.fingerprint(Path(str(parent["paths"]["klipreduce"]))),
    ]
    manifest = {
        "schema": 1,
        "purpose": "immutable KLIP Stage-H subpixel-response campaign",
        "root": str(root),
        "input_records": input_records,
        "software_records": software_records,
        "task_count": len(tasks),
        "positive_reductions_reused": len(tasks),
        "negative_reductions_required": len(tasks),
        "method_selection_performed": False,
        "planet_score_used_for_selection": False,
    }
    stage.write_json(output / "manifest.json", manifest)
    stage.write_json(output / "state.json", {
        "status": "prepared", "completed_negative_reductions": 0,
        "completed_analyses": 0, "task_count": len(tasks),
    })
    print(output / "manifest.json", flush=True)
    print(f"taskset -c 12-27 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "
          f"MKL_NUM_THREADS=1 python3 {frozen_runner} run {root}", flush=True)


def enable(root: Path) -> tuple[dict[str, object], dict[str, object]]:
    """Verify the frozen Stage-H boundary and return its protocol and parent."""
    root = root.resolve()
    output = output_root(root)
    protocol = read(output / "protocol.json")
    manifest = read(output / "manifest.json")
    stage.require(int(protocol["primary_mode"]) == PRIMARY_MODE and
                  str(protocol["primary_arm"]) == PRIMARY_ARM and
                  int(protocol["site_count"]) == SITE_COUNT and
                  list(protocol["modes"]) == stage.MODES and
                  tuple(protocol["models"]) == MODEL_NAMES and
                  len(protocol["tasks"]) == SITE_COUNT and
                  not manifest["method_selection_performed"] and
                  not manifest["planet_score_used_for_selection"],
                  "Stage-H frozen contract changed")
    stage.verify(manifest["input_records"] + manifest["software_records"])
    stage.require(stage.fingerprint(SCRIPT_PATH) in manifest["software_records"],
                  "Stage-H runner is outside the frozen software receipt")
    stage_g.verify_complete(root)
    parent = read(root / "protocol.json")
    return protocol, parent


def verify_task_receipt(path: Path) -> None:
    """Verify one completed Stage-H reduction or analysis receipt."""
    receipt = read(path)
    stage.require(receipt["status"] == "complete", f"incomplete receipt: {path}")
    stage.verify(receipt["products"])


def reduce_all(root: Path, protocol: dict[str, object], parent: dict[str, object]) -> None:
    """Run or verify the twelve paired negative KLIP reductions."""
    output = output_root(root)
    (output / "negative").mkdir(exist_ok=True)
    for index, task in enumerate(protocol["tasks"], start=1):
        directory = output / "negative" / str(task["name"])
        complete = directory / "complete.json"
        if complete.exists():
            verify_task_receipt(complete)
        else:
            if directory.exists():
                development.archive_incomplete(output, directory, "negative")
            directory.mkdir()
            started = time.monotonic()
            stage.run_command([str(value) for value in task["negative_command"]], directory,
                              directory / "run.log", development.environment(parent, True))
            product = directory / "finim.fits"
            stage_g.validate_reduction(product)
            stage.write_json(complete, {
                "status": "complete", "task": task["name"],
                "elapsed_seconds": time.monotonic() - started,
                "products": [stage.fingerprint(product),
                             stage.fingerprint(directory / "run.log"),
                             stage.fingerprint(directory / "run.command.json")],
            })
            verify_task_receipt(complete)
        stage.write_json(output / "state.json", {
            "status": "reducing", "completed_negative_reductions": index,
            "completed_analyses": len(list((output / "analysis").glob("*/complete.json")))
            if (output / "analysis").exists() else 0,
            "task_count": len(protocol["tasks"]),
        })
        print(f"Stage-H negative reduction {index}/{len(protocol['tasks'])}: "
              f"{task['name']}", flush=True)


def shifted_template(response: np.ndarray, validity: np.ndarray,
                     displacement: tuple[float, float]) -> tuple[np.ndarray, np.ndarray]:
    """Cubic-shift a full response, then return the central reviewed support."""
    response = np.asarray(response, dtype=np.float64)
    validity = np.asarray(validity, dtype=bool)
    stage.require(response.shape == validity.shape and response.ndim == 2 and
                  response.shape[0] == response.shape[1] and
                  response.shape[0] >= INTERPOLATION_GUARD,
                  "shifted response has an invalid full-stamp schema")
    guard = development.raw.crop(validity, INTERPOLATION_GUARD)
    stage.require(np.all(guard),
                  "shifted response lacks the central cubic-interpolation guard")
    shifted = cubic_shift(response, displacement, order=SHIFT_ORDER,
                          mode=SHIFT_MODE, prefilter=True)
    template = development.raw.crop(shifted, stage_f.SUPPORT)
    support = development.raw.crop(validity, stage_f.SUPPORT)
    stage.require(np.all(np.isfinite(template[support])) and np.all(support),
                  "shifted response lost reviewed 11-pixel support")
    return template, support


def response_lookup(parent: dict[str, object]) \
        -> tuple[dict[tuple[int, int], int], list[np.ndarray], list[np.ndarray]]:
    """Memory-map the integer exact-response coordinate lookup and mode cubes."""
    paths = development.response47.product_paths(Path(str(parent["parent_response"])))
    coordinates = np.asarray(fits.getdata(paths["coordinates"], memmap=True),
                             dtype=np.float64).T
    lookup = {(int(row), int(column)): index
              for index, (row, column, _, _) in enumerate(coordinates)}
    responses = [fits.getdata(path, memmap=True) for path in paths["responses"]]
    validities = [fits.getdata(path, memmap=True) for path in paths["validities"]]
    return lookup, responses, validities


def analyze_task(root: Path, protocol: dict[str, object], parent: dict[str, object],
                 task: dict[str, object], lookup: dict[tuple[int, int], int],
                 responses: list[np.ndarray], validities: list[np.ndarray]) -> None:
    """Compare integer, shifted, and paired exact responses for one site."""
    output = output_root(root) / "analysis" / str(task["name"])
    complete = output / "complete.json"
    if complete.exists():
        verify_task_receipt(complete)
        return
    if output.exists():
        development.archive_incomplete(output_root(root), output, "analysis")
    output.mkdir(parents=True)
    positive, positive_header = fits.getdata(task["positive"], header=True, memmap=True)
    negative_path = output_root(root) / "negative" / str(task["name"]) / "finim.fits"
    negative, negative_header = fits.getdata(negative_path, header=True, memmap=True)
    baseline, baseline_header = fits.getdata(protocol["paths"]["baseline"],
                                             header=True, memmap=True)
    positive = np.asarray(positive, dtype=np.float64)
    negative = np.asarray(negative, dtype=np.float64)
    baseline = np.asarray(baseline, dtype=np.float64)
    stage.require(positive.shape == negative.shape == baseline.shape ==
                  (len(stage.MODES), 128, 128) and
                  stage.read_modes(positive_header) ==
                  stage.read_modes(negative_header) ==
                  stage.read_modes(baseline_header) == stage.MODES,
                  "Stage-H response cubes changed schema")

    settings = dict(task["settings"])
    source, _, searches, _ = stage_f.analysis_geometry((128, 128), settings)
    expected_source = (float(task["source"]["row"]), float(task["source"]["column"]))
    stage.require(np.allclose(source, expected_source, rtol=0, atol=2e-12),
                  "Stage-H source geometry changed")
    query = min(searches, key=lambda value: math.hypot(value[0] - source[0],
                                                        value[1] - source[1]))
    stage.require(query in lookup, "Stage-H nearest query lacks an exact response")
    displacement = (source[0] - query[0], source[1] - query[1])
    stage.require(np.allclose(displacement, protocol["phase_offset_row_column"],
                              rtol=0, atol=2e-12),
                  "Stage-H nearest-query phase changed")
    source_index = lookup[query]
    optimized = stage_g.source_from_polar(
        (128, 128), float(parent["known_planet"]["separation"]),
        float(parent["known_planet"]["position_angle"]))
    radius_map = development.image_radius((128, 128))
    nominal = stage_f.policy_radius(float(radius_map[query[1], query[0]]))
    width = int(parent["training"]["candidate_specific_half_width_by_radius"][str(nominal)])
    contrast = float(task["contrast"])
    mode_results = {}
    primary_arrays = None
    for mode_index, mode in enumerate(stage.MODES):
        full_integer = np.asarray(responses[mode_index][source_index],
                                  dtype=np.float64).T
        full_validity = np.asarray(
            validities[mode_index][source_index], dtype=np.float64).T > 0.5
        integer = development.raw.crop(full_integer, stage_f.SUPPORT)
        integer_support = development.raw.crop(full_validity, stage_f.SUPPORT)
        shifted, shifted_support = shifted_template(full_integer, full_validity, displacement)
        central = ((development.stamp(positive[mode_index], query) -
                    development.stamp(negative[mode_index], query)) / (2 * contrast))
        one_sided_positive = ((development.stamp(positive[mode_index], query) -
                               development.stamp(baseline[mode_index], query)) / contrast)
        one_sided_negative = ((development.stamp(baseline[mode_index], query) -
                               development.stamp(negative[mode_index], query)) / contrast)
        support = (integer_support & shifted_support & np.isfinite(central) &
                   np.isfinite(one_sided_positive) & np.isfinite(one_sided_negative))
        stage.require(np.all(support) and np.count_nonzero(support) == stage_f.SUPPORT**2,
                      f"Stage-H common support changed for mode {mode}")
        fitted = stage_f.fit_planet_query(baseline[mode_index], query, searches,
                                          integer, support, optimized, width)
        models = {"integer": integer, "shifted": shifted,
                  "paired_exact": central}
        targets = {"central": central, "positive": one_sided_positive,
                   "negative": one_sided_negative}
        metrics = {
            target_name: {
                model_name: validation.response_fidelity(target, model, fitted)
                for model_name, model in models.items()
            }
            for target_name, target in targets.items()
        }
        self_metrics = metrics["central"]["paired_exact"]
        stage.require(all(np.isclose(value["cosine"], 1, rtol=0, atol=2e-12) and
                          np.isclose(value["projection_scale"], 1, rtol=0, atol=2e-12) and
                          value["best_scaled_relative_residual"] <= 3e-12
                          for value in self_metrics.values()),
                      f"Stage-H paired response self-check failed for mode {mode}")
        mode_results[str(mode)] = {
            "support_pixels": int(np.count_nonzero(support)),
            "policy_radius": nominal,
            "training_half_width": width,
            "metrics": metrics,
        }
        if mode == PRIMARY_MODE:
            primary_arrays = np.stack((central, one_sided_positive,
                                       one_sided_negative, integer, shifted))

    stage.require(primary_arrays is not None, "Stage-H primary response was not retained")
    header = positive_header.copy()
    header["HIERARCH STAGE H MODE"] = PRIMARY_MODE
    header["HIERARCH STAGE H LABELS"] = "central,positive,negative,integer,shifted"
    fits.writeto(output / "mode200_responses.fits",
                 primary_arrays.astype(np.float32), header)
    result = {
        "schema": 1,
        "task": task["name"],
        "site_index": task["site_index"],
        "base_site": task["base_site"],
        "source_row_column": list(source),
        "query_row_column": list(query),
        "phase_offset_row_column": list(displacement),
        "contrast": contrast,
        "modes": mode_results,
    }
    stage.write_json(output / "results.json", result)
    products = [stage.fingerprint(output / "results.json"),
                stage.fingerprint(output / "mode200_responses.fits")]
    stage.write_json(complete, {"status": "complete", "task": task["name"],
                                "products": products})
    verify_task_receipt(complete)


def analyze_all(root: Path, protocol: dict[str, object], parent: dict[str, object]) -> None:
    """Analyze or verify every paired subpixel response."""
    output = output_root(root)
    (output / "analysis").mkdir(exist_ok=True)
    lookup, responses, validities = response_lookup(parent)
    for index, task in enumerate(protocol["tasks"], start=1):
        verify_task_receipt(output / "negative" / str(task["name"]) / "complete.json")
        analyze_task(root, protocol, parent, task, lookup, responses, validities)
        stage.write_json(output / "state.json", {
            "status": "analyzing",
            "completed_negative_reductions": len(protocol["tasks"]),
            "completed_analyses": index, "task_count": len(protocol["tasks"]),
        })
        print(f"Stage-H analysis {index}/{len(protocol['tasks'])}: {task['name']}",
              flush=True)


def sample_summary(values: list[float]) -> dict[str, float | int]:
    """Summarize a finite sample with its sample standard deviation."""
    array = np.asarray(values, dtype=np.float64)
    stage.require(array.size == SITE_COUNT and np.all(np.isfinite(array)),
                  "Stage-H summary lost a finite site")
    return {"count": int(array.size), "mean": float(np.mean(array)),
            "sample_standard_deviation": float(np.std(array, ddof=1)),
            "minimum": float(np.min(array)), "maximum": float(np.max(array)),
            "values": [float(value) for value in array]}


def summarize(root: Path, protocol: dict[str, object]) -> None:
    """Aggregate the paired response comparisons and write review tables."""
    output = output_root(root)
    records = []
    for task in protocol["tasks"]:
        complete = output / "analysis" / str(task["name"]) / "complete.json"
        verify_task_receipt(complete)
        records.append(read(output / "analysis" / str(task["name"]) / "results.json"))
    by_mode = {}
    for mode in stage.MODES:
        targets = {}
        for target in ("central", "positive", "negative"):
            models = {}
            for model in MODEL_NAMES:
                metrics = {}
                for metric in METRIC_NAMES:
                    values = {}
                    for value in VALUE_NAMES:
                        values[value] = sample_summary([
                            float(record["modes"][str(mode)]["metrics"][target]
                                         [model][metric][value])
                            for record in records
                        ])
                    metrics[metric] = values
                models[model] = metrics
            targets[target] = models
        improvements = {}
        for target in ("central", "positive", "negative"):
            improvements[target] = {}
            for metric in METRIC_NAMES:
                integer_values = [record["modes"][str(mode)]["metrics"][target]
                                  ["integer"][metric] for record in records]
                shifted_values = [record["modes"][str(mode)]["metrics"][target]
                                  ["shifted"][metric] for record in records]
                improvements[target][metric] = {
                    "shifted_minus_integer_cosine": sample_summary([
                        float(second["cosine"]) - float(first["cosine"])
                        for first, second in zip(integer_values, shifted_values)
                    ]),
                    "integer_minus_shifted_residual": sample_summary([
                        float(first["best_scaled_relative_residual"]) -
                        float(second["best_scaled_relative_residual"])
                        for first, second in zip(integer_values, shifted_values)
                    ]),
                    "integer_minus_shifted_absolute_projection_error": sample_summary([
                        abs(float(first["projection_scale"]) - 1) -
                        abs(float(second["projection_scale"]) - 1)
                        for first, second in zip(integer_values, shifted_values)
                    ]),
                }
        by_mode[str(mode)] = {"targets": targets, "shift_improvements": improvements}

    result = {
        "schema": 1,
        "purpose": "KLIP integer, shifted, and paired-exact subpixel response comparison",
        "primary_mode": PRIMARY_MODE,
        "site_count": SITE_COUNT,
        "phase_offset_row_column": protocol["phase_offset_row_column"],
        "by_mode": by_mode,
        "limits": [
            "The paired-exact response is an oracle generated from the same observing sequence.",
            ("Fidelity cosine is matched-filter SNR retention under the stated "
             "covariance; no annular detection SNR is re-estimated."),
            ("A successful cubic shift supports a cheap phase model; a remaining "
             "residual motivates a measured phase grid."),
        ],
    }
    stage.write_json(output / "results.json", result)

    with (output / "results.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["target", "model", "metric", "cosine_mean", "cosine_std",
                         "projection_mean", "projection_std", "residual_mean",
                         "residual_std"])
        primary = by_mode[str(PRIMARY_MODE)]["targets"]
        for target in ("central", "positive", "negative"):
            for model in MODEL_NAMES:
                for metric in METRIC_NAMES:
                    values = primary[target][model][metric]
                    writer.writerow([
                        target, model, metric,
                        values["cosine"]["mean"],
                        values["cosine"]["sample_standard_deviation"],
                        values["projection_scale"]["mean"],
                        values["projection_scale"]["sample_standard_deviation"],
                        values["best_scaled_relative_residual"]["mean"],
                        values["best_scaled_relative_residual"]["sample_standard_deviation"],
                    ])

    lines = [
        "# KLIP Stage-H subpixel-response result", "",
        f"Twelve planet-like phase responses at offset "
        f"`({protocol['phase_offset_row_column'][0]:+.4f}, "
        f"{protocol['phase_offset_row_column'][1]:+.4f})` pixels are compared at mode "
        f"{PRIMARY_MODE}. Values are means +/- sample standard deviations.", "",
        "The paired central response is `(positive-negative)/(2*contrast)`. The positive",
        "target is the actual one-sided Stage-G injection response. Cosine is the retained",
        "optimal matched-filter SNR fraction under the stated covariance metric.", "",
        "## Paired central-response target", "",
        ("| Covariance metric | Model | Cosine | Projection scale | "
         "Best-scaled residual |"),
        "| :--- | :--- | ---: | ---: | ---: |",
    ]
    primary = by_mode[str(PRIMARY_MODE)]["targets"]
    labels = {"unweighted": "Identity", "raw_rectangular_m0p3": "Raw rectangular",
              "radial_hann_m0p1": "Radial Hann"}
    model_labels = {"integer": "Current integer", "shifted": "Cubic shifted",
                    "paired_exact": "Paired exact"}
    for metric in METRIC_NAMES:
        for model in ("integer", "shifted", "paired_exact"):
            values = primary["central"][model][metric]
            lines.append(
                f"| {labels[metric]} | {model_labels[model]} | "
                f"{values['cosine']['mean']:.6f} +/- "
                f"{values['cosine']['sample_standard_deviation']:.6f} | "
                f"{values['projection_scale']['mean']:.5f} +/- "
                f"{values['projection_scale']['sample_standard_deviation']:.5f} | "
                f"{values['best_scaled_relative_residual']['mean']:.4f} +/- "
                f"{values['best_scaled_relative_residual']['sample_standard_deviation']:.4f} |"
            )
    lines.extend(["", "## Actual positive-response target", "",
                  ("| Covariance metric | Model | Cosine | Projection scale | "
                   "Best-scaled residual |"),
                  "| :--- | :--- | ---: | ---: | ---: |"])
    for metric in METRIC_NAMES:
        for model in ("integer", "shifted", "paired_exact"):
            values = primary["positive"][model][metric]
            lines.append(
                f"| {labels[metric]} | {model_labels[model]} | "
                f"{values['cosine']['mean']:.6f} +/- "
                f"{values['cosine']['sample_standard_deviation']:.6f} | "
                f"{values['projection_scale']['mean']:.5f} +/- "
                f"{values['projection_scale']['sample_standard_deviation']:.5f} | "
                f"{values['best_scaled_relative_residual']['mean']:.4f} +/- "
                f"{values['best_scaled_relative_residual']['sample_standard_deviation']:.4f} |"
            )
    lines.extend(["", "## Interpretation guide", "",
                  ("- If the shifted model approaches the paired-exact model, the present "
                   "mismatch is predominantly registration and a cheap phase-shift model "
                   "is sufficient."),
                  ("- If a substantial residual remains, the KLIP response shape itself "
                   "depends on phase and the next implementation should measure a "
                   "fractional-phase response grid."),
                  ("- The paired-exact versus positive rows bound finite-amplitude "
                   "asymmetry that neither integer nor shifted registration can remove."),
                  ""])
    (output / "results.md").write_text("\n".join(lines), encoding="utf-8")


def verify_complete(root: Path) -> None:
    """Verify the final Stage-H receipt and every nested task receipt."""
    output = output_root(root)
    receipt = read(output / "complete.json")
    stage.require(receipt["status"] == "complete" and
                  not receipt["method_selection_performed"],
                  "Stage-H completion receipt changed")
    stage.verify(receipt["products"] + receipt["negative_receipts"] +
                 receipt["analysis_receipts"])
    for record in receipt["negative_receipts"] + receipt["analysis_receipts"]:
        verify_task_receipt(Path(str(record["path"])))


def finish(root: Path, protocol: dict[str, object]) -> None:
    """Write and verify the final immutable Stage-H receipt."""
    output = output_root(root)
    negative_receipts = [stage.fingerprint(
        output / "negative" / str(task["name"]) / "complete.json")
        for task in protocol["tasks"]]
    analysis_receipts = [stage.fingerprint(
        output / "analysis" / str(task["name"]) / "complete.json")
        for task in protocol["tasks"]]
    products = [stage.fingerprint(output / name)
                for name in ("results.json", "results.csv", "results.md")]
    stage.write_json(output / "complete.json", {
        "status": "complete", "negative_receipts": negative_receipts,
        "analysis_receipts": analysis_receipts, "products": products,
        "method_selection_performed": False,
    })
    stage.write_json(output / "state.json", {
        "status": "complete", "completed_negative_reductions": SITE_COUNT,
        "completed_analyses": SITE_COUNT, "task_count": SITE_COUNT,
    })
    verify_complete(root)


def run(root: Path) -> None:
    """Run the resumable negative reductions and response comparison."""
    root = root.resolve()
    protocol, parent = enable(root)
    output = output_root(root)
    if (output / "complete.json").exists():
        verify_complete(root)
        print(output / "results.md", flush=True)
        return
    reduce_all(root, protocol, parent)
    analyze_all(root, protocol, parent)
    summarize(root, protocol)
    finish(root, protocol)
    print(output / "results.md", flush=True)


def check() -> None:
    """Check command mutation, subpixel-shift direction, and fidelity conventions."""
    command = ["klipReduce", "--fake.contrast", "0.2", "--output.directory", "/old"]
    task = {"command": command, "contrast": 0.2, "name": "trial"}
    changed = negative_command(Path("/new"), task)
    stage.require(changed[2] == "-0.2" and changed[4].endswith("negative/trial") and
                  command[2] == "0.2" and command[4] == "/old",
                  "Stage-H negative command mutation failed")
    full = np.zeros((47, 47), dtype=np.float64)
    full[23, 23] = 1
    validity = np.ones(full.shape, dtype=bool)
    displacement = (-0.25, 0.5)
    shifted, support = shifted_template(full, validity, displacement)
    center = stage_f.SUPPORT // 2
    stage.require(np.all(support) and
                  np.sum(shifted[center - 1, :]) > np.sum(shifted[center + 1, :]) and
                  np.sum(shifted[:, center + 1]) > np.sum(shifted[:, center - 1]),
                  "Stage-H cubic shift direction changed")
    fitted = {"support": np.ones((stage_f.SUPPORT, stage_f.SUPPORT), dtype=bool),
              "scale": np.ones((stage_f.SUPPORT, stage_f.SUPPORT)),
              "raw_detail": {"covariance": np.eye(stage_f.SUPPORT**2)},
              "radial_detail": {"covariance": np.eye(stage_f.SUPPORT**2)}}
    template = np.arange(1, stage_f.SUPPORT**2 + 1, dtype=np.float64).reshape(
        stage_f.SUPPORT, stage_f.SUPPORT)
    metrics = validation.response_fidelity(template, template, fitted)
    stage.require(all(np.isclose(value["cosine"], 1) and
                      np.isclose(value["projection_scale"], 1) and
                      value["best_scaled_relative_residual"] <= 1e-14
                      for value in metrics.values()),
                  "Stage-H response fidelity convention changed")
    print("KLIP Stage-H subpixel-response checks passed", flush=True)


def parser() -> argparse.ArgumentParser:
    """Build the Stage-H command-line parser."""
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="action", required=True)
    subparsers.add_parser("check")
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("root", type=Path)
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("root", type=Path)
    return result


def main() -> None:
    """Dispatch the requested Stage-H action."""
    arguments = parser().parse_args()
    if arguments.action == "check":
        check()
    elif arguments.action == "prepare":
        prepare(arguments.root)
    else:
        run(arguments.root)


if __name__ == "__main__":
    main()
